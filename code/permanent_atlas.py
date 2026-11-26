"""Enumerate all maximal permanent catalogues of square and triangle targets."""

from argparse import ArgumentParser
from collections import defaultdict
from functools import cache
from itertools import combinations
import json
from pathlib import Path

DIRECTIONS = ((1, 0), (-1, 0), (0, 1), (0, -1))


COUNTS = {
    'triangles': (1, 1, 1, 3, 4, 12, 24, 66, 160, 448, 1186),
    'squares': (1, 1, 2, 5, 12, 35, 108, 369, 1285, 4655, 17073),
}


@cache
def square_shape(cells):
    images = []
    for _ in range(4):
        x0, y0 = min(cells)
        images.append(tuple(sorted((x - x0, y - y0) for x, y in cells)))
        cells = tuple((-y, x) for x, y in cells)
    return min(images)


@cache
def triangle_shape(cells):
    images = []
    for _ in range(6):
        x0 = 3 * (min(x for x, _ in cells) // 3)
        y0 = 3 * (min(y for _, y in cells) // 3)
        images.append(tuple(sorted((x - x0, y - y0) for x, y in cells)))
        cells = tuple((-y, x + y) for x, y in cells)
    return min(images)


def neighbors(kind, cell):
    x, y = cell
    if kind == 'squares':
        return tuple((x + dx, y + dy) for dx, dy in DIRECTIONS)
    step = 1 if x % 3 == 1 else -1
    return ((x + step, y + step), (x - 2 * step, y + step), (x + step, y - 2 * step))


@cache
def proper(kind, cells):
    return square_shape(cells) if kind == 'squares' else triangle_shape(cells)


@cache
def mirror(kind, cells):
    # Swapping the two coordinate axes reflects either lattice.
    return proper(kind, tuple(sorted((y, x) for x, y in cells)))


def free(kind, cells):
    return min(proper(kind, cells), mirror(kind, cells))


def free_shapes(kind, limit):
    unit = ((0, 0),) if kind == 'squares' else ((1, 1),)
    shapes = {unit}
    for stock in range(1, limit + 1):
        assert len(shapes) == COUNTS[kind][stock - 1], (kind, stock, len(shapes))
        yield stock, sorted(shapes)
        if stock < limit:
            next_shapes = set()
            for shape in shapes:
                border = {point for cell in shape for point in neighbors(kind, cell)} - set(shape)
                next_shapes.update(free(kind, tuple(sorted((*shape, point)))) for point in border)
            shapes = next_shapes


def connected_parts(kind, hand):
    adjacency = []
    for cell in hand:
        adjacent = set(neighbors(kind, cell))
        adjacency.append(sum(1 << i for i, point in enumerate(hand) if point in adjacent))
    masks = {1 << i: bits for i, bits in enumerate(adjacency)}
    for mask in range(1, 1 << len(hand)):
        if mask not in masks:
            continue
        frontier = masks[mask] & ~mask
        while frontier:
            bit = frontier & -frontier
            frontier ^= bit
            masks[mask | bit] = (masks[mask] | adjacency[bit.bit_length() - 1]) & ~(mask | bit)
    return {mask: proper(kind, tuple(point for i, point in enumerate(hand) if mask & (1 << i)))
            for mask in masks}


def geometry(kind, target):
    placements = [connected_parts(kind, hand) for hand in sorted({target, mirror(kind, target)})]
    shapes = sorted({body for group in placements for body in group.values()}, key=lambda s: (len(s), s))
    ids = {body: i for i, body in enumerate(shapes)}
    joins = defaultdict(set)
    bags = set()
    full = (1 << len(target)) - 1
    for group in placements:
        by_mask = {mask: ids[body] for mask, body in group.items()}
        by_cell = [[] for _ in target]
        for mask, result in by_mask.items():
            count = mask.bit_count()
            if 2 <= count <= len(target) // 2:
                for position in range(len(target)):
                    if mask & (1 << position):
                        by_cell[position].append((mask, result))
            left = (mask - 1) & mask
            while left:
                right = mask ^ left
                if left < right and left in by_mask and right in by_mask:
                    pair = tuple(sorted((by_mask[left], by_mask[right])))
                    joins[pair].add(result)
                left = (left - 1) & mask

        @cache
        def pack(remaining):
            if not remaining:
                return {()}
            bit = remaining & -remaining
            result = pack(remaining ^ bit).copy()
            for mask, part in by_cell[bit.bit_length() - 1]:
                if mask & remaining == mask:
                    result.update(tuple(sorted((*bag, part))) for bag in pack(remaining ^ mask))
            return result

        bags.update(pack(full))
        pack.cache_clear()
    return shapes, joins, bags


def small_catalogues(sizes, stock, joins, packable):
    parents = defaultdict(list)
    for pair, results in joins.items():
        for result in results:
            if sizes[result] <= stock // 2:
                parents[result].append(pair)
    small = [i for i, size in enumerate(sizes) if 1 < size <= stock // 2]
    selected = {0}

    def search(start, bags):
        # From four cells onward every rule has the unique dimer.
        if not small or small[0] in selected:
            yield frozenset(selected), bags
        for position in range(start, len(small)):
            part = small[position]
            if not any(a in selected and b in selected for a, b in parents[part]):
                continue
            additions = []
            for bag, mass in bags:
                extended = bag
                for total in range(mass + sizes[part], stock + 1, sizes[part]):
                    extended += (part,)
                    if extended not in packable:
                        break
                    additions.append((extended, total))
                else:
                    continue
                break
            else:
                selected.add(part)
                yield from search(position + 1, bags + additions)
                selected.remove(part)

    yield from search(0, [((), 0)])


def maximal_catalogues(shapes, joins, packable, stock):
    if stock == 1:
        return [(0,)]
    sizes = list(map(len, shapes))
    large = [i for i, size in enumerate(sizes) if stock // 2 < size < stock]
    terminals = {i for i, size in enumerate(sizes) if size == stock}
    parents = defaultdict(list)
    for pair, outcomes in joins.items():
        for body in outcomes:
            if sizes[body] > stock // 2:
                parents[body].append(pair)
    results = []
    for small, bags in small_catalogues(sizes, stock, joins, packable):
        # Residual inventories that cannot progress using only small joins.
        waiting = defaultdict(list)
        for count in {stock, *(stock - sizes[i] for i in large)}:
            for bag, mass in bags:
                if mass > count:
                    continue
                inventory = (0,) * (count - mass) + bag
                pairs = set(combinations(inventory, 2))
                if not any(joins.get(pair, set()) & small for pair in pairs):
                    waiting[count].append((set(inventory), pairs))

        viable = terminals.copy()
        for body in reversed(large):
            if all(any(joins.get((part, body), set()) & viable for part in residual)
                   for residual, _ in waiting[stock - sizes[body]]):
                viable.add(body)
        if not all(any(joins.get(pair, set()) & viable for pair in pairs)
                   for _, pairs in waiting[stock]):
            continue

        # Keep the large types actually constructible from the selected small ones.
        reachable = set(small)
        for body in sorted(viable):
            if any(a in reachable and b in reachable for a, b in parents[body]):
                reachable.add(body)
        assert reachable & terminals
        results.append(frozenset(reachable))

    maximal = []
    for candidate in sorted(set(results), key=len, reverse=True):
        if not any(candidate <= other for other in maximal):
            maximal.append(candidate)
    return sorted(tuple(sorted(catalogue)) for catalogue in maximal)


def record(kind, stock, index, target):
    shapes, joins, packable = geometry(kind, target)
    catalogues = maximal_catalogues(shapes, joins, packable, stock)
    active = sorted({part for catalogue in catalogues for part in catalogue})
    indices = {i: j for j, i in enumerate(active)}
    mirror_indices = {i: shapes.index(mirror(kind, shapes[i])) for i in active}
    mirror_groups = {min(catalogue, tuple(sorted(mirror_indices[i] for i in catalogue)))
                     for catalogue in catalogues}
    return {
        'id': f'{kind[0].upper()}{stock:02}-{index:04}',
        'shape': target,
        'parts': [shapes[i] for i in active],
        'factorizations': [[indices[i] for i in catalogue] for catalogue in catalogues],
        'mirror_groups': len(mirror_groups),
    }


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--kind', choices=COUNTS, action='append')
    parser.add_argument('--limit', type=int, choices=range(1, 12), default=10)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    for kind in args.kind or COUNTS:
        destination = args.output / kind
        destination.mkdir(parents=True, exist_ok=True)
        for stock, targets in free_shapes(kind, args.limit):
            records = [record(kind, stock, index, target) for index, target in enumerate(targets, 1)]
            (destination / f'{stock:02}.json').write_text(json.dumps(records, separators=(',', ':')) + '\n')
            print(kind, stock, len(records), sum(bool(row['factorizations']) for row in records),
                  sum(len(row['factorizations']) for row in records), flush=True)


if __name__ == '__main__':
    main()
