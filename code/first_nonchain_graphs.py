"""Contact graphs that the first non-chain body can have (proofs/late-cycles.md). Graphs only.

Let S(N0) be the first body that is not a chain. Then S(N0) = Y + x with Y = c_0 .. c_{M-1} a chain (M = N0 - 1),
and for every a + b = N0 the body splits into a copy of S(a) and a copy of S(b), which are chains. So the contact
graph G = (path c_0 .. c_{M-1}) + (vertex x joined to c_i, i in I) has this property:
    (*) for every a = 1 .. N0 // 2 the vertex set splits into two sets of sizes a and N0 - a, each inducing a path,
        with at least one edge between them.
This script lists, for each N0, the sets I (up to reversing the path) for which G is not a path and (*) holds.
A cell touches at most 6 others, so |I| <= 6; for |I| >= 2 consecutive members of I that are adjacent on the path
form a triangle with x.

It prints the survivors with |I| <= 2 (I = {0, M-1} is a ring; I = {0, j} is a cycle of length j + 2 with one tail),
counts the others, and checks the three statements of Lemma 12:
  (a) |I| = 1 (the body is a tree) only for N0 = 5 and N0 = 7;
  (b) for N0 >= 8: |I| >= 2, and x lies on a cycle of length at least N0 // 2 + 1;
  (c) for N0 >= 10: x touches an end cell of the chain.

usage: python3 first_nonchain_graphs.py [N0_max]   (default 14, which covers every case by Theorem 11)
"""
import sys
from itertools import combinations


def induced_path(adj, S):
    """Is the subgraph induced on the vertex set S a path (connected, max degree 2, |S| - 1 edges)?"""
    S = set(S)
    if not S:
        return False
    if len(S) == 1:
        return True
    deg = {v: len(adj[v] & S) for v in S}
    if max(deg.values()) > 2 or sum(deg.values()) != 2 * (len(S) - 1):
        return False
    start = next(iter(S))
    seen = {start}
    stack = [start]
    while stack:
        v = stack.pop()
        for w in adj[v] & S:
            if w not in seen:
                seen.add(w)
                stack.append(w)
    return len(seen) == len(S)


def split_sizes(M, I):
    """Sizes a for which path + x splits into two induced paths of sizes a and M + 1 - a (touching).
    x is in part A. By connectivity the other part B is an interval [p, q] of the path."""
    x = M
    adj = {v: set() for v in range(M + 1)}
    for v in range(M - 1):
        adj[v].add(v + 1); adj[v + 1].add(v)
    for i in I:
        adj[i].add(x); adj[x].add(i)
    sizes = set()
    V = set(range(M + 1))
    for p in range(M):
        for q in range(p, M):
            B = set(range(p, q + 1))
            A = V - B
            if induced_path(adj, A) and induced_path(adj, B) and any(adj[v] & B for v in A):
                sizes.add(len(A)); sizes.add(len(B))
    return sizes


def longest_cycle_through_x(M, I):
    I = sorted(I)
    return max(I) - min(I) + 2 if len(I) >= 2 else 0


def main(N0max):
    ok_tree, ok_cycle, ok_end = True, True, True
    for N0 in range(4, N0max + 1):
        M = N0 - 1
        out = []
        for d in range(1, 7):
            for I in combinations(range(M), d):
                if I > tuple(M - 1 - i for i in reversed(I)):
                    continue                                   # up to reversing the path
                if d == 1 and I[0] in (0, M - 1):
                    continue                                   # a path again
                need = set(range(1, N0 // 2 + 1))
                if need <= split_sizes(M, I):
                    out.append(I)
        by_d = {d: [I for I in out if len(I) == d] for d in range(1, 7)}
        for I in out:
            if len(I) == 1 and N0 not in (5, 7):
                ok_tree = False
            if N0 >= 8 and (len(I) < 2 or longest_cycle_through_x(M, I) < N0 // 2 + 1):
                ok_cycle = False
            if N0 >= 10 and 0 not in I and (M - 1) not in I:
                ok_end = False
        print(f"N0 = {N0:2d}: survivors by |I|: " + ", ".join(f"{d}: {len(by_d[d])}" for d in range(1, 7)))
        print(f"         |I| = 1: {by_d[1]};  |I| = 2: {by_d[2]}")
    print(f"For 4 <= N0 <= {N0max}:")
    print(f"  (a) |I| = 1 (a tree) only for N0 = 5 and 7: {'holds' if ok_tree else 'FAILS'}")
    print(f"  (b) N0 >= 8: |I| >= 2, and x lies on a cycle of length >= N0 // 2 + 1: {'holds' if ok_cycle else 'FAILS'}")
    print(f"  (c) N0 >= 10: x touches an end cell of the chain: {'holds' if ok_end else 'FAILS'}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 14)
