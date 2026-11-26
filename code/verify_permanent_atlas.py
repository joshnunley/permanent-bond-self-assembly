"""Audit the published stock tables with direct inventory and packing tests."""

from argparse import ArgumentParser
from itertools import combinations
import json
from pathlib import Path

from permanent_atlas import COUNTS, free_shapes, geometry, mirror, neighbors


def inventories(stock, sizes, parts):
    if stock == 0:
        yield ()
    for i, part in enumerate(parts):
        if sizes[part] <= stock:
            for rest in inventories(stock - sizes[part], sizes, parts[i:]):
                yield (part, *rest)


def exact(sizes, joins, catalogue):
    built = {0}
    for part in sorted(catalogue - built):
        if not any(part in joins.get((a, b), ()) for a in built for b in built if a <= b):
            return False
        built.add(part)
    return all(len(state) == 1 or any(joins.get(pair, set()) & catalogue
                                     for pair in combinations(state, 2))
               for state in inventories(max(sizes), sizes, sorted(catalogue)))


def stock_theorem(kind, target, sizes, joins):
    stock = len(target)
    if stock <= 7:
        def match(cells, singles):
            if not cells:
                return singles == 0
            cell = min(cells)
            rest = cells - {cell}
            return (singles and match(rest, singles - 1)) or any(
                match(rest - {other}, singles) for other in rest.intersection(neighbors(kind, cell))
            )
        return bool(match(set(target), stock % 2))

    def join(a, b):
        return joins.get(tuple(sorted((a, b))), set())

    terminals = {part for part, size in enumerate(sizes) if size == stock}
    unit, domino = 0, 1
    fours = join(domino, domino)

    def packs_three(a, b, c):
        return any(join(part, last) & terminals for left, right, last in ((a, b, c), (a, c, b), (b, c, a))
                   for part in join(left, right))

    if stock == 8:
        return any(join(four, four) & terminals for four in fours)
    if stock == 10:
        return any(packs_three(four, four, domino) for four in fours)
    return any(packs_three(four, four, unit) for four in fours) or any(
        packs_three(trimer, trimer, trimer) and
        any(packs_three(five, domino, domino) and packs_three(five, trimer, unit)
            for five in join(trimer, domino))
        for trimer in join(unit, domino)
    )


def direct_maximal(sizes, joins):
    maximal = set()
    for count in range(len(sizes) - 1, 0, -1):
        for rest in combinations(range(1, len(sizes)), count):
            candidate = frozenset({0, *rest})
            if sizes[max(candidate)] == max(sizes) and not any(candidate <= other for other in maximal):
                if exact(sizes, joins, candidate):
                    maximal.add(candidate)
    if len(sizes) == 1:
        maximal.add(frozenset({0}))
    return maximal


def verify(kind, data):
    checked = 0
    for stock, targets in free_shapes(kind, 10):
        rows = json.loads((data / kind / f'{stock:02}.json').read_text())
        assert len(rows) == len(targets)
        brute = 0
        for index, (row, target) in enumerate(zip(rows, targets), 1):
            assert row['id'] == f'{kind[0].upper()}{stock:02}-{index:04}'
            assert tuple(map(tuple, row['shape'])) == target
            shapes, joins, _ = geometry(kind, target)
            sizes = list(map(len, shapes))
            shape_ids = {shape: i for i, shape in enumerate(shapes)}
            parts = [shape_ids[tuple(map(tuple, part))] for part in row['parts']]
            catalogues = {frozenset(parts[i] for i in group) for group in row['factorizations']}
            assert len(catalogues) == len(row['factorizations'])
            assert set(parts) == set().union(*catalogues)
            assert bool(catalogues) == stock_theorem(kind, target, sizes, joins), row['id']
            for group in catalogues:
                assert exact(sizes, joins, group), row['id']
                assert not any(group < other for other in catalogues)
            mirrored = {i: shape_ids[mirror(kind, shape)] for i, shape in enumerate(shapes)}
            pairs = {min(tuple(sorted(group)), tuple(sorted(mirrored[i] for i in group))) for group in catalogues}
            assert len(pairs) == row['mirror_groups']
            assert all(frozenset(mirrored[i] for i in group) in catalogues for group in catalogues)
            if len(shapes) <= 12:
                assert catalogues == direct_maximal(sizes, joins), row['id']
                brute += 1
        checked += len(rows)
        print(f'{kind} {stock}: {len(rows)} statuses and all listed catalogues checked, {brute} full subset searches',
              flush=True)
    return checked


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('data', type=Path)
    parser.add_argument('--kind', choices=COUNTS, action='append')
    args = parser.parse_args()
    print(f'{sum(verify(kind, args.data) for kind in args.kind or COUNTS)} targets checked')
