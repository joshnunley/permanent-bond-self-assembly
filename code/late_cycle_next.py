"""Reduce the two-angle cases at sizes six and seven by one-cell extensions."""

from collections import deque
from fractions import Fraction
from itertools import combinations

from late_cycles import admissible, dot, impose, solve_sum, sub, whole_torus
from one_turn import cw


def induced_path(adj, vertices):
    if len(vertices) == 1:
        return tuple(vertices)
    degree = {v: len(adj[v] & vertices) for v in vertices}
    ends = sorted(v for v in vertices if degree[v] == 1)
    if len(ends) != 2 or any(d == 0 or d > 2 for d in degree.values()):
        return None
    path = [ends[0]]
    while choices := (adj[path[-1]] & vertices) - set(path):
        assert len(choices) == 1
        path.append(choices.pop())
    return tuple(path) if len(path) == len(vertices) else None


def extension_graphs(edges, size, n0):
    """Keep the splits for which both parts predate the first non-chain body."""
    vertices = set(range(size + 1))
    for degree in range(1, min(size, 6) + 1):
        for neighbors in combinations(range(size), degree):
            adj = [set() for _ in vertices]
            for a, b in edges:
                adj[a].add(b)
            for i in neighbors:
                adj[i].add(size)
                adj[size].add(i)
            groups = []
            for k in range(max(2, size + 2 - n0), (size + 1) // 2 + 1):
                paths = []
                for p in combinations(range(size + 1), k):
                    if 2 * k == size + 1 and size in p:
                        continue
                    A = set(p)
                    u = induced_path(adj, A)
                    v = induced_path(adj, vertices - A)
                    if u and v:
                        paths.append((u, v))
                if not paths:
                    break
                groups.append(paths)
            else:
                yield neighbors, sorted(groups, key=len)


def old_turns(word, offset, basis):
    def encode(letter):
        if not letter:
            return (0,) * len(basis), Fraction(0)
        k = abs(letter) - 1
        sign = 1 if letter > 0 else -1
        return tuple(sign * b[k] for b in basis), (sign * offset[k]) % 1

    return {
        tuple(encode(t) for t in w)
        for w in (word, cw.rev(word), cw.neg(word), cw.rev(cw.neg(word)))
    }


def path_equations(path, turns, directions):
    headings = [directions[a, b] for a, b in zip(path, path[1:])]
    observed = [(sub(b, a), (q - p) % 1) for (a, p), (b, q) in zip(headings, headings[1:])]
    return {
        tuple(sorted((sub(a, b), (q - p) % 1) for (a, p), (b, q) in zip(observed, w) if (a, p) != (b, q)))
        for w in turns
    }


def shortest_headings(edges, size, a, b):
    queue = deque([[a]])
    seen = {a}
    while queue:
        path = queue.popleft()
        if path[-1] == b:
            return [edges[i, j] for i, j in zip(path, path[1:])]
        for w in range(size):
            if (path[-1], w) in edges and w not in seen:
                seen.add(w)
                queue.append(path + [w])
    raise AssertionError((a, b))


def compose(offset, basis, inner_offset, inner_basis):
    """Compose two affine torus parametrizations, with basis vectors as columns."""
    rows = tuple(tuple(v[i] for v in basis) for i in range(len(offset)))
    return (
        tuple((a + dot(row, inner_offset)) % 1 for a, row in zip(offset, rows)),
        tuple(tuple(dot(row, b) for row in rows) for b in inner_basis),
    )


def extend_cosets(edges, size, n0, words, offset, basis):
    """Parameterize positions of a new cell, using all splits into two old chains."""
    k = len(basis)
    assert k == 2 and size + 1 <= 8
    old = {i: old_turns(w, offset, basis) for i, w in enumerate(words[1:], 1)}
    result = set()

    for neighbors, groups in extension_graphs(edges, size, n0):
        degree = len(neighbors)
        directions = {key: (a + (0,) * degree, p) for key, (a, p) in edges.items()}
        expected = {
            i: {tuple((a + (0,) * degree, p) for a, p in w) for w in turns}
            for i, turns in old.items()
        }
        for j, i in enumerate(neighbors, k):
            row = tuple(int(t == j) for t in range(k + degree))
            directions[size, i] = row, Fraction(0)
            directions[i, size] = row, Fraction(1, 2)

        choices = []
        for split_opts in groups:
            alternatives = {
                tuple(sorted(left + right))
                for a, b in split_opts
                for left in path_equations(a, expected[len(a)], directions)
                for right in path_equations(b, expected[len(b)], directions)
            }
            if () not in alternatives:
                choices.append(alternatives)

        states = {whole_torus(k + degree)}
        for alternatives in sorted(choices, key=len):
            next_states = set()
            for equations in alternatives:
                active = states
                for row, phase in equations:
                    active = {
                        item for a, b in active for item in impose(a, b, row, phase)
                    }
                    if not active:
                        break
                next_states.update(active)
            states = next_states
            if not states:
                break

        for a, b in zip(neighbors, neighbors[1:]):
            old_path = shortest_headings(directions, size, a, b)
            row, p = directions[size, b]
            law = tuple(sorted([*old_path, directions[size, a], (row, (p + Fraction(1, 2)) % 1)]))
            states = {item for off, bas in states for item in solve_sum(0, off, bas, law)}
            if not states:
                break

        for off, bas in states:
            point = off[:k]
            columns = tuple(v[:k] for v in bas)
            turns, turn_columns = compose(offset, basis, point, columns)
            if admissible(len(offset), turns, turn_columns):
                heading = tuple(v[k] for v in bas), off[k]
                result.add((neighbors[0], point, columns, heading, turns, turn_columns))
    return result
