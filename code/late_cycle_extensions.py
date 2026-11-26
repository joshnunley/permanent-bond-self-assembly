"""Certify extension of late cycle candidates over points and one-angle cosets."""

from collections import Counter
from fractions import Fraction
from functools import cache
from itertools import combinations
from math import lcm
import sys

import numpy as np

from exact_core import Bodies, CycField, family_search, lemma15_candidates
from formal_multi import FormalBodies, Ring, formal_lemma15
from late_cycle_next import extend_cosets
from late_cycles import admissible, dot, first_bodies, heading_chain, univariate_roots
from special import Polys


@cache
def compare_unit(field, vector):
    """Use the certified interval branch for every non-unit length test."""
    if vector == field.zero:
        return -1
    if vector in field.unit:
        return 0
    return field._cmp1_hp(vector)


class ExactBodies(Bodies):
    def valid_new(self, cells, x):
        return all(compare_unit(self.F, self.F.sub(x, c)) >= 0 for c in cells)


class WholeTorusRing(Ring):
    def vals(self, arr):
        return np.zeros(arr.shape[: -self.r - 1])

    def rads(self, arr):
        # The older candidate generator uses these only to skip impossible
        # length matches. Over the whole torus, log every non-formal match.
        return np.full(arr.shape[: -self.r - 1], np.inf)


class FormalRelaxation(FormalBodies):
    """Permit variable overlaps, and log equalities over the whole torus."""

    def touching_pairs(self, arr, vals, rads):
        pairs = []
        for i, j in combinations(range(len(arr)), 2):
            g = arr[j] - arr[i]
            root = self.ring.unit_of(g)
            if root is not None:
                k, p = root
                pairs.extend([(i, j, k, p), (j, i, tuple(-x for x in k), -p)])
            elif len(self.ring.support(g)) >= 2:
                self.log("contact", g)
        return pairs

    def valid_new_arr(self, arr, x):
        for a in arr:
            support = self.ring.monomials(x - a)
            if not support or (len(support) == 1 and compare_unit(self.ring.K, support[0][1]) < 0):
                return False
        return True

    def lookup(self, img, t):
        ids = []
        for row in img:
            i = self.pos[t].get(row.tobytes())
            if i is None:
                # Any successful specialization must place this first missing point.
                for a in self.arr[t]:
                    g = row - a
                    if len(self.ring.support(g)) >= 2:
                        self.log("zero", g)
                return None
            ids.append(i)
        return ids


def chain_points(word, offset, basis):
    """A point is a sum of roots indexed by exponent vector and phase."""
    points = [Counter()]
    count = len(offset)
    for h in heading_chain(word):
        a = tuple(dot(h[:count], v) for v in basis)
        p = (dot(h[:count], offset) + Fraction(h[-1], 6)) % 1
        points.append(points[-1] + Counter({(a, p): 1}))
    return points


def starting_bodies(candidate):
    words, neighbor, (offset, basis), heading = candidate
    starts = [[Counter()]]
    starts.extend(chain_points(w, offset, basis) for w in words[2:])
    old = starts[-1]
    x = Counter(old[neighbor])
    x[heading] -= 1
    starts.append(old + [x])
    return starts


def coefficient_order(starts):
    return lcm(2, *(p.denominator for pts in starts for g in pts for (_, p), c in g.items() if c))


def integer_bounds(field):
    # The searches work with at most fifteen cells. Sums, unit multiplications
    # and the products in Polys.norm are bounded by this conservative margin.
    height = max(abs(c) for row in field.pw for c in row)
    assert 10**8 * field.phi**3 * height**5 < 2**63, (field.Q, height)


def register(registry, starts, is_formal, n0):
    if is_formal:
        ring = registry.ring
        order = ring.M

        def encode(poly):
            return ring.elem({(*a, int(p * order) % order): c for (a, p), c in poly.items() if c})
    else:
        field = registry.F
        order = field.Q

        def encode(poly):
            return field.from_poly({int(p * order) % order: c for (_, p), c in poly.items() if c})

    ids = [None]
    for size, pts in enumerate(starts, 1):
        arr = [encode(g) for g in pts]
        if is_formal:
            valid = registry.valid_set(arr)
        else:
            valid = all(registry.valid_new(arr[:i], p) for i, p in enumerate(arr))
        if not valid:
            return None
        i = registry.add(arr)
        if 3 <= size < n0 and not registry.is_chain(i):
            return None
        if not is_formal and size == n0 and registry.is_chain(i):
            return None
        ids.append(i)
    return ids


def certify_fourteen_cell_ring(registry, ids):
    field = registry.F
    assert field.Q % 14 == 0 and len(ids) == 15
    turn = field.Q // 14
    cells = [field.zero]
    for step in range(13):
        cells.append(field.add(cells[-1], field.pw[step * turn]))
    for size, body in enumerate(ids[1:], 1):
        arc = registry.add(cells[:size])
        assert registry.shape(arc) == registry.shape(body), size


def run(registry, starts, n0, limit, is_formal):
    ids = register(registry, starts, is_formal, n0)
    if ids is None:
        return 0
    ext = formal_lemma15(registry).ext if is_formal else lemma15_candidates(registry).ext

    def candidates(family, size):
        new_size = size + 1
        if new_size < len(ids):
            return [ids[new_size]]
        a = min(4, new_size - 4)
        assert a >= 3
        return {
            child
            for y in family[size]
            for p in {a, new_size - a}
            for child in ext(y, family[p][0])
        }

    reach, final = family_search(
        registry, limit, candidates, start=[[None, (ids[1],)]]
    )
    assert not final, (n0, limit, starts, final)
    reached = max(reach)
    if reached == 14:
        assert n0 == 14 and not is_formal
        certify_fourteen_cell_ring(registry, ids)
    return reached


def evaluate(starts, point):
    result = []
    for pts in starts:
        body = []
        for g in pts:
            value = Counter()
            for (a, p), c in g.items():
                value[(), (dot(a, point) + p) % 1] += c
            body.append(value)
        result.append(body)
    return result


@cache
def exact_point(starts, n0, limit):
    starts = [[dict(g) for g in pts] for pts in starts]
    field = CycField(coefficient_order(starts))
    integer_bounds(field)
    return run(ExactBodies(field), starts, n0, limit, False)


def point_check(starts, value, n0, limit):
    evaluated = evaluate(starts, value)
    key = tuple(tuple(tuple(sorted((k, c) for k, c in g.items() if c)) for g in pts) for pts in evaluated)
    return exact_point(key, n0, limit)


def one_angle(registry, starts, n0, limit, offset, basis, counts):
    reached = run(registry, starts, n0, limit, True)
    counts["one angle", reached] += 1
    polynomials = Polys(registry.ring)
    exceptional = set()
    for (kind, _), data in registry.logs.items():
        g, _ = polynomials.relation(kind, data)
        if not g:
            # For the ratio log the denominator is a difference of two
            # units. Its nonzero roots, if any, all lie on the unit circle.
            # Identically equal moduli force its quotient with the numerator
            # to be a monomial times a constant of modulus one. That constant
            # would have to be a root of unity in the coefficient field at
            # any admissible match, which would have made the ratio formal.
            assert kind == "ratio", kind
            continue
        exceptional.update(univariate_roots(registry.ring.M, tuple(sorted(g.items()))))
    for t in sorted(exceptional):
        turns = tuple((a + b * t) % 1 for a, b in zip(offset, basis[0]))
        if admissible(len(offset), turns, ()):
            reached = point_check(starts, (t,), n0, limit)
            counts["special point", reached] += 1


def substitute(poly, offset, basis):
    result = Counter()
    for (a, p), c in poly.items():
        result[tuple(dot(a, v) for v in basis), (dot(a, offset) + p) % 1] += c
    return result


def apply_extension(registry, body, starts, extension):
    neighbor, point, columns, heading, turns, turn_columns = extension
    ring = registry.ring
    canonical = []
    for g in registry.arr[body]:
        poly = {
            (a, Fraction(j, ring.M)): c
            for a, coeff in ring.monomials(g)
            for j, c in enumerate(coeff)
            if c
        }
        canonical.append(substitute(poly, point, columns))
    refined = [[substitute(g, point, columns) for g in pts] for pts in starts]
    x = Counter(canonical[neighbor])
    x[heading] -= 1
    return refined + [canonical + [x]], turns, turn_columns


def check_case(starts, n0, limit, words, offset, basis, counts):
    dimension = len(basis)
    if dimension == 0:
        reached = point_check(starts, (), n0, limit)
        counts["fixed point", reached] += 1
        return
    assert dimension in (1, 2)
    ring = WholeTorusRing(coefficient_order(starts), 32 if dimension == 1 else 16, dimension)
    integer_bounds(ring.K)
    registry = FormalRelaxation(ring)
    if dimension == 1:
        one_angle(registry, starts, n0, limit, offset, basis, counts)
        return

    ids = register(registry, starts, True, n0)
    if ids is None:
        counts["two angles", 0] += 1
        return
    assert n0 in (6, 7) and len(starts) <= 7, (n0, len(starts))
    body = ids[-1]
    edges = {(a, b): (k, Fraction(p % ring.M, ring.M)) for a, b, k, p in registry.pairs[body]}
    extensions = extend_cosets(edges, len(starts), n0, words, offset, basis)
    counts["two angles", len(starts)] += 1
    for extension in sorted(extensions):
        refined, turns, columns = apply_extension(registry, body, starts, extension)
        check_case(refined, n0, limit, words, turns, columns, counts)


def check(n0):
    assert 6 <= n0 <= 14
    limit = 15
    first = first_bodies(n0)
    sizes = Counter(len(candidate[2][1]) for candidate in first)
    counts = Counter()
    for candidate in sorted(first):
        offset, basis = candidate[2]
        check_case(starting_bodies(candidate), n0, limit, candidate[0], offset, basis, counts)
    print(f"First non-chain size {n0}. Candidate counts by dimension: {sorted(sizes.items())}")
    for (kind, reached), count in sorted(counts.items()):
        print(f"  {kind:13s} max size {reached:2d}: {count}")
    print(f"No prefix reaches size {limit}.")
    if n0 == 14:
        print("Every prefix reaching size 14 consists of arcs of a regular 14-cell ring.")


if __name__ == "__main__":
    for n in map(int, sys.argv[1:]):
        check(n)
