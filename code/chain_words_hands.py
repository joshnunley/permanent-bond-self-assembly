"""Chain prefixes with hands allowed: how long can the bodies of a merge-closed family stay chains without the
small ones being constant or alternating?  (proofs/all-side-counts.md)

Setting. F(1), ..., F(K) satisfy the definition of a merge-closed family for every split a + b <= K, and every body
in them is a chain (its contact graph is a path). Nothing is assumed about F(N) for N > K. Then nothing rules out
both hands of a body at N > K/2 (the "no two hands" argument needs F(2N) to be a chain), so, unlike chain_words.py,
this search lets F(N) be one class or a mirror pair, exactly as merge_closed.py does for bodies.

Words. A chain of N cells is its word of N - 2 turns. Two chains are properly congruent iff their words are equal
or one is the other reversed and negated; the mirror image negates the word. Property 3 for a + b = N and a chain T:
T splits into a copy of A and a copy of B only by one cut, after a cells or after b cells.

Letters. Only equalities between turns matter. A letter is 0 (straight on) or +-k, k = 1, 2, ... a turn size in
order of first appearance. No bound on the number of sizes is assumed: at each step the new letter is 0, +-an old
size, or one new size (its sign is immaterial: flipping the sign of a size that occurs nowhere else is a renaming).

Search. Depth-first over all sequences F(1..N). F(N + 1) is {T}, {mirror T} or {T, mirror T} for a candidate T that
is a member of F(N) plus one cell at either end (property 3 with b = 1 gives at least one such member).

Output. For each N: how many prefixes reach N, how many distinct shape sequences they are, and how many of those
are not "simple" up to level m, for m = 6, 8, 11, 14. A prefix is simple up to m when S(m) has a constant or
alternating word (then so has every S(i), i <= m, because S(i) is a sub-chain of S(m)). For N = 10..13 it prints
the shape sequences whose S(N) is neither constant nor alternating (the exceptions of Lemma 10).

usage: python3 chain_words_hands.py [K_max]   (default 30)
"""
import sys
import time


def neg(w):
    return tuple(-x for x in w)


def rev(w):
    return tuple(reversed(w))


def canon(w):
    """Class under rotation and translation."""
    return min(w, rev(neg(w)))


def shape(w):
    """Class under all isometries."""
    return min(w, rev(w), neg(w), rev(neg(w)))


def sub(w, i, j):
    """Word of the sub-chain of cells i..j (1-based, inclusive) of the chain with word w."""
    return w[i - 1:j - 2] if j - i >= 2 else ()


def simple(w):
    """Constant (including all zeros) or alternating word."""
    if len(w) <= 1:
        return True
    if all(x == w[0] for x in w):
        return True
    return w[0] != 0 and all(w[i] == -w[i - 1] for i in range(1, len(w)))


def ok(opt, F, M):
    """Property 3 for every split a + b = M, with F(M) = opt."""
    for a in range(1, M // 2 + 1):
        b = M - a
        for A in F[a]:
            for B in F[b]:
                good = False
                for t in opt:
                    if (canon(sub(t, 1, a)) == A and canon(sub(t, a + 1, M)) == B) or \
                       (canon(sub(t, 1, b)) == B and canon(sub(t, b + 1, M)) == A):
                        good = True
                        break
                if not good:
                    return False
    return True


def search(Kmax, one_hand=False):
    reach = [0] * (Kmax + 1)
    seqs = [set() for _ in range(Kmax + 1)]          # shape sequences reaching N

    def dfs(F, N):
        reach[N] += 1
        seqs[N].add(tuple(shape(F[i][0]) for i in range(3, N + 1)))
        if N == Kmax:
            return
        M = N + 1
        cands = set()
        for w in F[N]:
            sizes = sorted({abs(x) for x in w if x})
            new = (sizes[-1] + 1) if sizes else 1
            letters = [0, new] + [s * k for k in sizes for s in (1, -1)]
            for v in (w, rev(neg(w))):
                for s in letters:
                    t = canon(v + (s,))
                    cands.add(min(t, canon(neg(t))))
        for t in sorted(cands):
            tm = canon(neg(t))
            options = [(t,)] if t == tm else ([(t,), (tm,)] if one_hand else [(t,), (tm,), (t, tm)])
            for opt in options:
                if ok(opt, F, M):
                    dfs(F + [opt], M)

    dfs([None, ((),), ((),)], 2)
    return reach, seqs


def report(Kmax, one_hand):
    t0 = time.time()
    reach, seqs = search(Kmax, one_hand)
    kind = "one class at each N (as chain_words.py, but with no bound on the number of turn sizes)" if one_hand \
        else "one class or a mirror pair at each N"
    print(f"chain prefixes, {kind}; K_max = {Kmax}  ({time.time() - t0:.0f} s)")
    print("  N: prefixes reaching N / distinct shape sequences / of these, not simple up to m = 6, 8, 11, 14")
    for N in range(3, Kmax + 1):
        ns = {m: sum(1 for s in seqs[N] if m <= N and not simple(s[m - 3])) for m in (6, 8, 11, 14)}
        print(f"  {N:3d}: {reach[N]:7d} / {len(seqs[N]):6d} / " + "  ".join(
            f"m={m}: {ns[m] if m <= N else '-'}" for m in (6, 8, 11, 14)))
    for m in (6, 8, 11, 14):
        K = next((N for N in range(m, Kmax + 1) if all(simple(s[m - 3]) for s in seqs[N])), None)
        print(f"  least K <= {Kmax} such that every prefix reaching K is simple up to {m}: {K}")
    if not one_hand:
        fmt = lambda w: " ".join(f"{x:+d}" if x else "0" for x in w) if w else "(no turn)"
        for N in range(10, min(Kmax, 13) + 1):
            exc = sorted(s for s in seqs[N] if not simple(s[-1]))
            print(f"  shape sequences reaching N = {N} with S({N}) neither constant nor alternating: {len(exc)}")
            for s in exc:
                print(f"    S({N}) = {fmt(s[-1])}")
                print("      earlier: " + "; ".join(f"S({i + 3}) = {fmt(w)}" for i, w in enumerate(s[:-1]) if i + 3 >= 5))
    last = sorted({s[-1] for s in seqs[Kmax]})
    print(f"  shapes S({Kmax}) of the prefixes reaching {Kmax}: {len(last)}")
    for w in last[:12]:
        print("    " + " ".join(f"{x:+d}" if x else "0" for x in w))
    sys.stdout.flush()


if __name__ == "__main__":
    Kmax = int(sys.argv[1]) if len(sys.argv) > 1 else 30
    report(min(Kmax, 24), one_hand=True)
    report(Kmax, one_hand=False)
