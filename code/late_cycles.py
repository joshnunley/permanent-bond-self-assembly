"""Exact direction and closure constraints on the first non-chain body."""

from collections import Counter, defaultdict
from fractions import Fraction
from functools import cache
from itertools import combinations
from math import gcd, lcm

from sympy import factorint, isprime

from one_turn import chain_prefixes, cw, n_sizes
from special import field, orders_up_to_degree, vanishes_at

ZERO = (0, 0, 0, 0)
HALF = (0, 0, 0, 3)


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def turn(letter):
    sign = 1 if letter > 0 else -1
    return tuple(sign * (abs(letter) == i) for i in (1, 2, 3)) + (0,)


def heading_chain(word):
    """Headings as three turn coefficients and a constant in sixths of a turn."""
    headings = [ZERO]
    for letter in word:
        headings.append(add(headings[-1], turn(letter)))
    return headings


@cache
def graph_splits(m):
    """Contact graphs of a path of length m-1 plus vertex m."""
    vertices = set(range(m + 1))
    result = []
    for degree in range(1, min(m, 6) + 1):
        for neighbors in combinations(range(m), degree):
            if degree == 1 and neighbors[0] in (0, m - 1):
                continue
            adj = [set() for _ in vertices]
            for i in range(m - 1):
                adj[i].add(i + 1)
                adj[i + 1].add(i)
            for i in neighbors:
                adj[i].add(m)
                adj[m].add(i)

            splits = [[] for _ in range((m + 1) // 2 + 1)]
            for left in range(m):
                for right in range(left + 1, m + 1):
                    other = vertices.difference(range(left, right))
                    ends = sorted(v for v in other if len(adj[v] & other) == 1)
                    if len(other) > 1 and (
                        len(ends) != 2 or any(len(adj[v] & other) > 2 for v in other)
                    ):
                        continue
                    path = [ends[0] if ends else m]
                    while next_vertices := (adj[path[-1]] & other) - set(path):
                        assert len(next_vertices) == 1
                        path.append(next_vertices.pop())
                    if len(path) == len(other):
                        splits[min(len(path), right - left)].append((path, left, right))
            if all(splits[1:]):
                result.append((neighbors, splits))
    return result


def split_equations(words, neighbors, splits):
    """Disjunctions of equations h[j] - h[i] = value for the new headings."""
    old_word = words[-1]
    old = heading_chain(old_word)
    m = len(old) + 1
    directions = {}
    for i, a in enumerate(old):
        directions[i, i + 1] = (0, a)
        directions[i + 1, i] = (0, add(a, HALF))
    for j, i in enumerate(neighbors, 1):
        directions[m, i] = (j, ZERO)
        directions[i, m] = (j, HALF)

    groups = []
    for j, (a, b) in enumerate(zip(neighbors, neighbors[1:]), 1):
        if b == a + 1:  # An equilateral triangle on the old edge.
            groups.append({
                ((0, j, add(old[a], (0, 0, 0, 3 + sign))),
                 (0, j + 1, add(old[a], (0, 0, 0, -sign))))
                for sign in (-1, 1)
            })
        if b == a + 2:  # The other common neighbor of the two old endpoints.
            groups.append({((0, j, add(old[a + 1], HALF)), (0, j + 1, old[a]))})

    for witnesses in splits[2:]:
        alternatives = set()
        for path, left, right in witnesses:
            interval = old_word[left : right - 2] if right - left >= 3 else ()
            if cw.shape(interval) != words[right - left]:
                continue
            w = words[len(path)]
            path_angles = [directions[a, b] for a, b in zip(path, path[1:])]
            for turns in {w, cw.neg(w), cw.rev(w), cw.neg(cw.rev(w))}:
                eqs = []
                for (i, a), (j, b), letter in zip(path_angles, path_angles[1:], turns):
                    value = add(sub(a, b), turn(letter))
                    if i != j or value != ZERO:
                        eqs.append((i, j, value))
                alternatives.add(tuple(eqs))
        if not alternatives:
            return None
        if () not in alternatives:
            groups.append(alternatives)
    return sorted(groups, key=len)


def normalize(a):
    return (*a[:3], a[3] % 6)


def angle_systems(groups, degree):
    """Solve for the new headings up to root headings and old-angle equations."""
    systems = set()

    def visit(k, roots, weights, relations):
        if k == len(groups):
            systems.add((tuple(roots), tuple(weights), frozenset(relations)))
            return
        for equations in groups[k]:
            rr, ww, cc = roots[:], weights[:], set(relations)
            for i, j, value in equations:
                p, q = rr[i], rr[j]
                value = normalize(add(value, sub(ww[i], ww[j])))
                if p == q:
                    coeffs, phase = value[:3], value[3]
                    if not any(coeffs):
                        if phase:
                            break
                    elif phase == 0 and max(map(abs, coeffs)) == 1 and sum(map(abs, coeffs)) <= 2:
                        # A turn would be zero, or two turn sizes would coincide.
                        break
                    else:
                        cc.add(min(value, normalize(sub(ZERO, value))))
                else:
                    if p > q:
                        p, q = q, p
                        value = sub(ZERO, value)
                    for j, root in enumerate(rr):
                        if root == q:
                            rr[j] = p
                            ww[j] = normalize(add(ww[j], value))
            else:
                visit(k + 1, rr, ww, cc)

    visit(0, list(range(degree + 1)), [ZERO] * (degree + 1), set())
    return systems


def bezout(a, b):
    if not b:
        return (1, 0, a) if a > 0 else (-1, 0, -a)
    u, v, g = bezout(b, a % b)
    return v, u - (a // b) * v, g


def whole_torus(dimension):
    return (Fraction(0),) * dimension, tuple(
        tuple(int(i == j) for i in range(dimension)) for j in range(dimension)
    )


@cache
def impose(offset, basis, row, rhs):
    """Intersect a torus coset with row*theta = rhs modulo one."""
    coeffs = [dot(row, v) for v in basis]
    rhs = (rhs - dot(row, offset)) % 1
    if not any(coeffs):
        return () if rhs else ((offset, basis),)
    vectors = list(basis)
    for k in range(1, len(basis)):
        a, b = coeffs[0], coeffs[k]
        if b:
            s, t, g = bezout(a, b)
            u, v = vectors[0], vectors[k]
            vectors[0] = tuple(s * x + t * y for x, y in zip(u, v))
            vectors[k] = tuple((-b // g) * x + (a // g) * y for x, y in zip(u, v))
            coeffs[0] = g
    g = coeffs[0]
    if g < 0:
        vectors[0] = tuple(-x for x in vectors[0])
        g = -g
    return tuple(
        (tuple((x + ((rhs + k) / g) * y) % 1 for x, y in zip(offset, vectors[0])),
         tuple(vectors[1:]))
        for k in range(g)
    )


@cache
def admissible(count, offset, basis):
    """Reject cosets forcing an invalid turn or coincident turn sizes."""
    for i in range(count):
        a = offset[i]
        A = tuple(v[i] for v in basis)
        if not any(A) and not (0 < a < Fraction(1, 3) or Fraction(2, 3) < a < 1):
            return False
        for j in range(i):
            b = offset[j]
            B = tuple(v[j] for v in basis)
            if (a == b and A == B) or ((a + b) % 1 == 0 and A == tuple(-x for x in B)):
                return False
    return True


@cache
def finite_field_root(order):
    """A certified primitive root of the given order in a prime field."""
    p = order + 1
    while p < 2**64 and not isprime(p):
        p += order
    assert p < 2**64
    divisors = factorint(order)
    for a in range(2, p):
        x = pow(a, (p - 1) // order, p)
        if all(pow(x, order // d, p) != 1 for d in divisors):
            assert pow(x, order, p) == 1
            return x, p
    raise AssertionError(order)


@cache
def univariate_roots(order, polynomial):
    g = dict(polynomial)
    span = max(g)[0] - min(g)[0]
    if not span:
        return ()
    terms = [(k[0], j, c) for k, cs in polynomial for j, c in enumerate(cs) if c]
    roots = []
    for q in orders_up_to_degree(order, span):
        total_order = lcm(order, q)
        root, p = finite_field_root(total_order)
        powers = [pow(root, i, p) for i in range(total_order)]
        for k in range(q):
            if gcd(k, q) != 1:
                continue
            a, b = k * (total_order // q), total_order // order
            if sum(c * powers[(i * a + j * b) % total_order] for i, j, c in terms) % p:
                continue
            t = Fraction(k, q)
            if vanishes_at(order, g, (t,)):
                roots.append(t)
    return tuple(roots)


@cache
def vanishing_coefficient(phases):
    if not phases:
        return True
    order = lcm(*(f.denominator for f in phases))
    if order == 1:
        return False
    exps = [int(f * order) % order for f in phases]
    x, p = finite_field_root(order)
    if sum(pow(x, k, p) for k in exps) % p:
        return False
    co = [0] * order
    for k in exps:
        co[k] += 1
    return vanishes_at(order, {(): co}, ())


@cache
def minimal_patterns(size):
    """All rotations anchored at one term, from Table 1 of Poonen–Rubinstein."""
    assert 2 <= size <= 10
    patterns = set()

    def record(roots):
        assert len(set(roots)) == size
        for origin in roots:
            patterns.add(tuple(sorted((r - origin) % 1 for r in roots if r != origin)))

    for prime in (2, 3, 5, 7):
        count = size - prime
        if not 0 <= count < prime or (count and prime < 5):
            continue
        for replaced in combinations(range(prime), count):
            roots = [
                (Fraction(j, prime) + shift) % 1
                for j in range(prime)
                for shift in ((Fraction(1, 6), Fraction(-1, 6)) if j in replaced else (0,))
            ]
            record(roots)
    if size == 10:
        # R_7:R_5, with the shared term rotated to 1.
        record([Fraction(j, 7) for j in range(1, 7)] + [
            (Fraction(1, 2) + Fraction(j, 5)) % 1 for j in range(1, 5)
        ])
    return tuple(sorted(patterns))


@cache
def coefficients(offset, basis, terms):
    groups = defaultdict(list)
    for row, phase in terms:
        exp = tuple(dot(row, v) for v in basis)
        groups[exp].append((dot(row, offset) + phase) % 1)
    return {
        exp for exp, phases in groups.items()
        if not vanishing_coefficient(tuple(sorted(phases)))
    }


@cache
def solve_sum(turn_count, offset, basis, terms):
    """All torsion zeros of the unit sum, as a finite union of cosets."""
    active = coefficients(offset, basis, terms)
    if not active:
        return ((offset, basis),)
    if len(active) == 1:
        return ()
    terms = tuple(
        (a, p) for a, p in terms if tuple(dot(a, v) for v in basis) in active
    )
    result = set()
    pivot, origin = terms[0]

    def match(off, bas, rest, phases):
        if not phases:
            result.update(solve_sum(turn_count, off, bas, rest))
            return
        remaining = coefficients(off, bas, terms)
        if not remaining:
            result.add((off, bas))
            return
        if len(remaining) == 1:
            return
        for j, (a, p) in enumerate(rest):
            if j and (a, p) == rest[j - 1]:
                continue
            for oo, bb in impose(off, bas, sub(a, pivot), (phases[0] + origin - p) % 1):
                if admissible(turn_count, oo, bb):
                    match(oo, bb, rest[:j] + rest[j + 1 :], phases[1:])

    for size in range(2, len(terms) + 1):
        if len(terms) - size == 1:
            continue
        for pattern in minimal_patterns(size):
            match(offset, basis, terms[1:], pattern)
    return tuple(result)


def unit_length(terms):
    """The univariate Laurent polynomial |sum terms|^2 - 1."""
    groups = defaultdict(Counter)
    groups[0][0] = -1
    for (k, p), a in terms.items():
        for (j, q), b in terms.items():
            groups[k - j][(p - q) % 6] += a * b
    K = field(6)
    result = {}
    for k, coeffs in groups.items():
        v = K.from_poly(coeffs)
        if any(v):
            result[k,] = v
    return tuple(sorted(result.items()))


@cache
def one_turn_placements(steps, neighbors, i, heading):
    """Specialize the anchored cell using only its necessary unit lengths."""
    h, _, _, phase = heading
    # At turn zero the old centres are collinear with unit spacing.
    # Two required contacts at least three indices apart cannot both hold.
    assert neighbors[-1] - neighbors[0] >= 3
    polynomials = set()
    for j in neighbors:
        terms = Counter({(h, phase): 1})
        for step in steps[min(i, j) : max(i, j)]:
            terms[step[0], 0] += 1 if j > i else -1
        p = unit_length(terms)
        if p:
            polynomials.add(p)
    assert polynomials
    g = min(polynomials, key=lambda p: p[-1][0][0] - p[0][0][0])
    result = []
    for t in univariate_roots(6, g):
        if not admissible(1, (t,), ()):
            continue
        if all(vanishes_at(6, dict(p), (t,)) for p in polynomials):
            result.append((i, ((t,), ()), ((), (h * t + Fraction(phase, 6)) % 1)))
    return tuple(result)


def one_turn_anchors(words, steps, witnesses):
    """Read a new heading from the part containing x in the split 4+(N-4)."""
    old_word = words[-1]
    x = len(steps) + 1
    anchors = set()

    def old(a, b):
        return steps[a] if a < b else add(steps[b], HALF)

    for path, left, right in witnesses:
        if cw.shape(old_word[left : right - 2]) != words[right - left]:
            continue
        w = words[len(path)]
        variants = {w, cw.neg(w), cw.rev(w), cw.neg(cw.rev(w))}
        for j, (a, b, c) in enumerate(zip(path, path[1:], path[2:])):
            if a == x:
                anchors.update((b, normalize(sub(old(b, c), turn(v[j])))) for v in variants)
                break
            if c == x:
                anchors.update((b, normalize(add(add(old(a, b), turn(v[j])), HALF))) for v in variants)
                break
        else:
            raise AssertionError(path)
    return anchors


def first_bodies(n0):
    """Candidate first bodies as a chain plus one specified unit displacement."""
    m = n0 - 1
    prefixes = chain_prefixes(m)[m]
    words = {
        (None, (), (), *(cw.shape(F[i][0]) for i in range(3, m + 1)))
        for F in prefixes
        if n_sizes(F[m][0]) >= 1
    }
    result = set()
    for seq in sorted(words):
        r = n_sizes(seq[-1])
        steps = tuple(heading_chain(seq[-1]))
        for neighbors, splits in graph_splits(m):
            if r == 1 and n0 >= 8:
                for i, heading in one_turn_anchors(seq, steps, splits[4]):
                    for item in one_turn_placements(steps, neighbors, i, heading):
                        result.add((seq, *item))
                continue
            groups = split_equations(seq, neighbors, splits)
            if groups is None:
                continue
            for roots, weights, constraints in angle_systems(groups, len(neighbors)):
                free = sorted(set(roots) - {0})
                states = [whole_torus(r + len(free))]
                for eq in sorted(constraints):
                    row = eq[:r] + (0,) * len(free)
                    states = [
                        (off, bas)
                        for a, b in states
                        for off, bas in impose(a, b, row, Fraction(-eq[-1], 6))
                        if admissible(r, off, bas)
                    ]
                if not states:
                    continue

                def heading(j):
                    row = weights[j][:r] + tuple(int(roots[j] == u) for u in free)
                    return row, Fraction(weights[j][-1], 6)

                cycles = set()
                for j, (a, b) in enumerate(zip(neighbors, neighbors[1:]), 1):
                    old = [(e[:r] + (0,) * len(free), Fraction(e[-1], 6)) for e in steps[a:b]]
                    row, p = heading(j + 1)
                    cycles.add(tuple(sorted([*old, heading(j), (row, (p + Fraction(1, 2)) % 1)])))
                for cycle in sorted(cycles, key=len):
                    states = {item for a, b in states for item in solve_sum(r, a, b, cycle)}
                    if not states:
                        break
                row, p = heading(1)
                for off, bas in states:
                    case = off[:r], tuple(v[:r] for v in bas)
                    new_heading = tuple(dot(row, v) for v in bas), (dot(row, off) + p) % 1
                    result.add((seq, neighbors[0], case, new_heading))
    return result
