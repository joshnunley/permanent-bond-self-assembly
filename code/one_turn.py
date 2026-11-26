"""The first non-chain body when the chain before it has one turn size: reduction to finitely many turns, and the
exact searches at those turns.  (proofs/late-cycles.md)

Setting. F is merge-closed, S(N0) is the first non-chain body, N0 >= 8, and all nonzero turns of the chain
S(N0 - 1) have the same size tau (0 < tau < pi). Put z = exp(i tau). (For N0 >= 11 this always holds, by the word
model of ../nonchain/chain_words_hands.py: see the table printed by step "words".)

Step "formal" (all n at once). The chains S(1), ..., S(N0 - 1) are words in the letters 0, +1, -1, and their centres
are Laurent polynomials in z with integer coefficients. By Lemma 15 of nonchain_attempt.md, S(N0) = Y + x with
x = h(c): either
  (i) x = y1 + rho * (c - a)(z) or y1 + rho * (c - a)(1/z), where a, b are adjacent cells of S_p (p = 3 or N0 - 3),
      y1, y2 are adjacent cells of Y, and rho = (y2 - y1)/(b - a), resp. (y2 - y1)/conj(b - a), is +-z^k; or
  (ii) x is the middle of a copy of S(3) whose ends are the two end cells of Y; then
      |c_last - c_first|^2 = |1 + z^w|^2 (w = +-1 or 0 the letter of S(3)), a nonzero Laurent polynomial
      equation E(z) = 0.
In case (i), x - c_j =: g_j(z) is a Laurent polynomial for every cell c_j of Y, and x touches c_j exactly when
G_j(z) := g_j(z) g_j(1/z) - 1 = 0. G_j is identically 0 exactly when g_j = +-z^a ("formal contact").
The script lists every G_j that is not identically 0 and every E, over all single-turn-size chain prefixes of the
word model with 7 <= N0 - 1 <= 13, and finds all roots of unity z that are roots of one of them: the set Z*.
For z not in Z*, x has formal contacts only; two formal contacts are with cells c_i, c_{i+2}; so the contact graph
of S(N0) is a tree or has a single 4-ring, which contradicts Lemma 12 of nonchain_attempt.md for N0 >= 8.

Step "exact" (one z in Z* at a time). For z = exp(2 pi i k / d) the script runs the complete search in exact
arithmetic in Z[zeta_Q], Q = lcm(d, 2), under the coin condition: chains with turns 0, +-tau up to N = 5, then
all Lemma 15 candidates (rotations by roots of unity of the field), chains only up to N = 7. A chain prefix that
reaches N = 14 is not continued (Theorem 11: N0 <= 14). It reports the prefixes whose first non-chain body is
S(N0), 8 <= N0 <= 14, and the largest N they reach. The claim is that none reaches N_max.

usage: python3 -B one_turn.py words
       python3 -B one_turn.py remaining              (the chain prefixes with two or three turn sizes)
       python3 -B one_turn.py formal                 (writes out/one_turn_special.txt)
       python3 -B one_turn.py exact d k [N_max]      (one turn; default N_max 20)
       python3 -B one_turn.py exact_all [N_max] [part/parts]   (every (d, k) of out/one_turn_special.txt)
"""
import math
import sys
import time
from collections import defaultdict

import numpy as np

import chain_words_hands as cw

from exact_core import CycField, Bodies, family_search, lemma15_candidates          # noqa: E402


# ---------------------------------------------------------------------------------------------------------------
# the word model: all chain prefixes F(1..K) (copied from chain_words_hands.search, but returning the prefixes)
# ---------------------------------------------------------------------------------------------------------------

def chain_prefixes(Kmax):
    out = defaultdict(list)

    def dfs(F, N):
        out[N].append(list(F))
        if N == Kmax:
            return
        M = N + 1
        cands = set()
        for w in F[N]:
            sizes = sorted({abs(x) for x in w if x})
            new = (sizes[-1] + 1) if sizes else 1
            letters = [0, new] + [s * k for k in sizes for s in (1, -1)]
            for v in (w, cw.rev(cw.neg(w))):
                for s in letters:
                    t = cw.canon(v + (s,))
                    cands.add(min(t, cw.canon(cw.neg(t))))
        for t in sorted(cands):
            tm = cw.canon(cw.neg(t))
            options = [(t,)] if t == tm else [(t,), (tm,), (t, tm)]
            for opt in options:
                if cw.ok(opt, F, M):
                    dfs(F + [opt], M)

    dfs([None, ((),), ((),)], 2)
    return out


def n_sizes(w):
    return len({abs(x) for x in w if x})


def words_report(Kmax=13):
    pre = chain_prefixes(Kmax)
    print("chain prefixes F(1..K) of the word model (one class or a mirror pair at each N), by the number of turn "
          "sizes of S(K)")
    print("  K: prefixes; with 0, 1, 2, 3 turn sizes")
    for K in range(3, Kmax + 1):
        cnt = defaultdict(int)
        for F in pre[K]:
            cnt[n_sizes(F[K][0])] += 1
        print(f"  {K:2d}: {len(pre[K]):5d}; " + ", ".join(str(cnt[i]) for i in range(4)))
    print("So for K >= 10 every prefix has at most one turn size.")


def remaining_report(Kmax=9):
    """The chain prefixes that Theorem 4.1 does not exclude: shape sequences S(3..K) with two or three turn sizes
    in S(K), 5 <= K <= 9 (K = N0 - 1)."""
    pre = chain_prefixes(Kmax)
    fmt = lambda w: " ".join(f"{x:+d}" if x else "0" for x in w)
    print("chain prefixes with two or three turn sizes: what is left of (G)")
    print("(shape sequences S(3) | S(4) | ... | S(K), K = N0 - 1; letters are turn sizes in order of appearance)")
    for K in range(5, Kmax + 1):
        seqs = {2: set(), 3: set()}
        for F in pre[K]:
            ns = n_sizes(F[K][0])
            if ns >= 2:
                seqs[ns].add(tuple(cw.shape(F[i][0]) for i in range(3, K + 1)))
        print(f"K = {K} (N0 = {K + 1}): {len(seqs[2])} shape sequences with two sizes, {len(seqs[3])} with three")
        for ns in (3, 2):
            for sq in sorted(seqs[ns]):
                print(f"  [{ns}] " + " | ".join(fmt(w) for w in sq))


# ---------------------------------------------------------------------------------------------------------------
# Laurent polynomials in z as dicts {exponent: coefficient}
# ---------------------------------------------------------------------------------------------------------------

def padd(p, q, s=1):
    r = dict(p)
    for e, c in q.items():
        v = r.get(e, 0) + s * c
        if v:
            r[e] = v
        else:
            r.pop(e, None)
    return r


def pmul(p, q):
    r = defaultdict(int)
    for e1, c1 in p.items():
        for e2, c2 in q.items():
            r[e1 + e2] += c1 * c2
    return {e: c for e, c in r.items() if c}


def pshift(p, k, s=1):
    return {e + k: s * c for e, c in p.items()}


def pinv(p):
    """p(1/z): the complex conjugate on |z| = 1."""
    return {-e: c for e, c in p.items()}


def chain_poly(word):
    """Centres c_0..c_{N-1} and direction exponents e_0..e_{N-2} of the chain with the given word (letters 0, +-1):
    c_{i+1} = c_i + z^{e_i}, e_0 = 0, e_i = e_{i-1} + word[i-1]."""
    cs = [{}, {0: 1}]
    es = [0]
    for s in word:
        es.append(es[-1] + s)
        cs.append(padd(cs[-1], {es[-1]: 1}))
    return cs, es


def norm_key(g):
    """Key of g up to g -> +-z^a g(z^{+-1}); G = g(z) g(1/z) - 1 has the same roots of unity for all of these."""
    def one(h):
        m = min(h)
        items = tuple(sorted((e - m, c) for e, c in h.items()))
        neg = tuple((e, -c) for e, c in items)
        return min(items, neg)
    return min(one(g), one(pinv(g)))


TOT = {}


def totients(n):
    phi = list(range(n + 1))
    for p in range(2, n + 1):
        if phi[p] == p:
            for k in range(p, n + 1, p):
                phi[k] -= phi[k] // p
    return phi


def cyclo_coeffs(d, cache={}):
    if d not in cache:
        from sympy import Poly, cyclotomic_poly
        from sympy.abc import x
        cache[d] = [int(c) for c in Poly(cyclotomic_poly(d, x), x).all_coeffs()]      # high -> low
    return cache[d]


def poly_rem_zero(coeffs_high_low, d):
    """Is the integer polynomial divisible by the d-th cyclotomic polynomial? Exact integer long division."""
    div = cyclo_coeffs(d)
    a = list(coeffs_high_low)
    n = len(div) - 1
    for i in range(len(a) - n):
        c = a[i]
        if c:
            for j in range(1, n + 1):
                a[i + j] -= c * div[j]
            a[i] = 0
    return not any(a)


class RootFinder:
    """All d such that the d-th cyclotomic polynomial divides a given integer polynomial. A float evaluation at
    exp(2 pi i / d) screens the candidates d (phi(d) <= degree); a value of modulus > 1e-6 is certainly not 0,
    because the rounding error is below 1e-10 for the polynomials met here (at most 400 terms of size <= 100); the
    remaining d are decided by exact division."""

    def __init__(self, maxdeg):
        self.maxdeg = maxdeg
        lim = 2 * maxdeg * maxdeg + 10
        phi = totients(lim)
        self.ds = np.array([d for d in range(1, lim + 1) if phi[d] <= maxdeg])
        self.phis = np.array([phi[d] for d in self.ds])
        self.theta = 2 * np.pi / self.ds

    def orders(self, P):
        """P: dict exponent -> coeff (Laurent). Returns the list of d with Phi_d | P."""
        if not P:
            raise ValueError("zero polynomial")
        m = min(P)
        deg = max(P) - m
        assert deg <= self.maxdeg, deg
        es = np.array([e - m for e in P])
        cs = np.array([P[e] for e in P], dtype=float)
        vals = np.exp(1j * np.outer(self.theta, es)) @ cs
        out = []
        coeffs = None
        for idx in np.nonzero((np.abs(vals) < 1e-6) & (self.phis <= deg))[0]:
            d = int(self.ds[idx])
            if coeffs is None:
                coeffs = [0] * (deg + 1)
                for e, c in P.items():
                    coeffs[deg - (e - m)] = c
            if poly_rem_zero(coeffs, d):
                out.append(d)
        return out


def formal_step(Kmin=7, Kmax=13, out_path=None):
    t0 = time.time()
    pre = chain_prefixes(Kmax)
    RF = RootFinder(120)
    cache = {}
    special = defaultdict(lambda: defaultdict(int))           # d -> N0 -> number of (candidate, cell) pairs
    stats = defaultdict(lambda: defaultdict(int))

    def roots_of_contact(g):
        key = norm_key(g)
        if key not in cache:
            G = padd(pmul(g, pinv(g)), {0: 1}, -1)
            cache[key] = RF.orders(G) if G else None
        return cache[key]

    for K in range(Kmin, Kmax + 1):
        N0 = K + 1
        combos = set()
        for F in pre[K]:
            if n_sizes(F[K][0]) > 1:
                continue
            for Y in F[K]:
                for p in {3, K - 2}:
                    combos.add((Y, F[p][0]))
        stats[N0]["(Y, S_p) pairs"] = len(combos)
        ring_done = set()
        for (Y, S) in sorted(combos):
            Yc, Ye = chain_poly(Y)
            Sc, Se = chain_poly(S)
            p = len(Sc)
            M = len(Yc)
            # case (ii): x between the two ends of Y, as the middle of a copy of S(3) (which has a nonzero turn)
            if p == 3 and (Y, S) not in ring_done:
                ring_done.add((Y, S))
                D = padd(Yc[-1], Yc[0], -1)
                ends = padd(Sc[2], Sc[0], -1)                      # 1 + z^w, w the letter of S(3)
                E = padd(pmul(D, pinv(D)), pmul(ends, pinv(ends)), -1)
                assert E, "ring equation is identically zero"
                for d in RF.orders(E):
                    special[d][N0] += 1
                stats[N0]["ring equations"] += 1
            # case (i)
            for flip in (False, True):
                Sx = [pinv(c) for c in Sc] if flip else Sc
                for ci in range(p):
                    pair = None
                    for a in range(p - 1):
                        if a != ci and a + 1 != ci:
                            pair = (a, a + 1)
                            break
                    if pair is None:
                        continue                       # p = 3 and c the middle cell: case (ii)
                    a, b = pair
                    ea = -Se[a] if flip else Se[a]     # b - a = z^ea
                    rel = padd(Sx[ci], Sx[a], -1)
                    for j in range(M - 1):
                        for (i1, sgn) in ((j, 1), (j + 1, -1)):      # (y1, y2) = (c_j, c_j+1) or (c_j+1, c_j)
                            # y2 - y1 = sgn z^{Ye[j]};  rho = sgn z^{Ye[j] - ea}
                            x = padd(Yc[i1], pshift(rel, Ye[j] - ea, sgn))
                            stats[N0]["candidates"] += 1
                            formal = []
                            coincide = False
                            nonformal = []
                            for t in range(M):
                                g = padd(x, Yc[t], -1)
                                if not g:
                                    coincide = True
                                    break
                                if len(g) == 1 and abs(next(iter(g.values()))) == 1:
                                    formal.append(t)
                                else:
                                    nonformal.append(g)
                            if coincide:
                                stats[N0]["candidates equal to a cell of Y"] += 1
                                continue
                            assert len(formal) <= 2 and (len(formal) < 2 or formal[1] - formal[0] == 2), formal
                            stats[N0][f"candidates with {len(formal)} formal contacts"] += 1
                            for g in nonformal:
                                ds = roots_of_contact(g)
                                stats[N0]["contact equations"] += 1
                                for d in ds:
                                    special[d][N0] += 1
    lines = []
    lines.append(f"one turn size, first non-chain body S(N0), {Kmin + 1} <= N0 <= {Kmax + 1}: the special turns "
                 f"({time.time() - t0:.0f} s)")
    for N0 in sorted(stats):
        lines.append(f"  N0 = {N0}: " + "; ".join(f"{k}: {v}" for k, v in stats[N0].items()))
    lines.append(f"  distinct contact polynomials (up to g -> +-z^a g(z^+-1)): {len(cache)}")
    ds = sorted(special)
    lines.append(f"  orders d of the roots of unity z in Z*: {len(ds)} values, largest {ds[-1] if ds else None}")
    npairs = 0
    for d in ds:
        ks = [k for k in range(1, d // 2 + 1) if math.gcd(k, d) == 1 and 2 * k != d]
        npairs += len(ks)
        lines.append(f"  d = {d}: arises for N0 in {sorted(special[d])}; turns 2 pi k / d, k in {ks}")
    lines.append(f"  turns to search (z and 1/z are mirror images; z = 1 and z = -1 are not turns): {npairs}")
    text = "\n".join(lines)
    print(text)
    if out_path:
        with open(out_path, "w") as f:
            f.write(text + "\n")
    return ds


# ---------------------------------------------------------------------------------------------------------------
# exact search at one turn
# ---------------------------------------------------------------------------------------------------------------

def chain_order(Rg, i):
    """Cells of a chain body in path order (indices)."""
    nb = Rg.contact_graph(i)
    N = len(Rg.bodies[i])
    if N == 1:
        return [0]
    ends = [t for t in range(N) if len(nb[t]) == 1]
    order = [ends[0]]
    while len(order) < N:
        order.append(next(t for t in nb[order[-1]] if t not in order[-2:]))
    return order


def exact_turn(d, k, Nmax=20, verbose=False):
    """Exact search at z = exp(2 pi i k / d). Returns a dict with the outcome."""
    Q = d if d % 2 == 0 else 2 * d
    zexp = k * (Q // d)
    F = CycField(Q)
    Rg = Bodies(F)
    l15 = lemma15_candidates(Rg)
    chain_cache = {}

    def all_chains(Fam):
        key = tuple(Fam[1:])
        if key not in chain_cache:
            chain_cache[key] = all(Rg.is_chain(i) for opt in Fam[1:] for i in opt)
        return chain_cache[key]

    def chain_ext(i):
        cells = Rg.bodies[i]
        out = set()
        if len(cells) == 1:
            return {Rg.add([F.zero, F.one])}
        order = chain_order(Rg, i)
        for (e, e1) in ((order[0], order[1]), (order[-1], order[-2])):
            kk = F.unit[F.sub(cells[e], cells[e1])]
            for s in (0, 1, -1):
                xn = F.add(cells[e], F.pw[(kk + s * zexp) % Q])
                if xn in Rg.pos[i] or not Rg.valid_new(cells, xn):
                    continue
                t = Rg.add(cells + (xn,))
                if Rg.is_chain(t):
                    out.add(t)
        return out

    def cands(Fam, N):
        M = N + 1
        chains = all_chains(Fam)
        if chains and N >= 14:
            return []
        if M <= 5:
            out = set()
            for A in Fam[N]:
                out |= chain_ext(A)
            return out
        out = l15(Fam, N)
        if M <= 7:
            out = [t for t in out if Rg.is_chain(t)]
        return out

    res = {"d": d, "k": k, "Q": Q, "by_N0": defaultdict(lambda: [0, 0]), "reach": None}

    def on_prefix(Fam, N):
        if N >= 8 and not Rg.is_chain(Fam[N][0]):
            N0 = next(M for M in range(3, N + 1) if not Rg.is_chain(Fam[M][0]))
            if N == N0:
                res["by_N0"][N0][0] += 1
            res["by_N0"][N0][1] = max(res["by_N0"][N0][1], N)

    # start: the single cell; S(2); S(3) is found by chain_ext (turn 0 or +-tau)
    reach, finals = family_search(Rg, Nmax, cands, on_prefix=on_prefix)
    res["reach"] = dict(reach)
    res["alive"] = sum(1 for Fam in finals if not all_chains(Fam))
    res["bodies"] = len(Rg.bodies)
    res["hp"] = F.hp_calls
    return res


def fmt(res):
    s = f"d = {res['d']}, k = {res['k']} (turn {360 * res['k'] / res['d']:.3f} deg, Q = {res['Q']}): "
    if not res["by_N0"]:
        s += "no prefix has a first non-chain body S(N0) with N0 >= 8"
    else:
        s += "; ".join(f"N0 = {N0}: {c} prefixes, last alive at N = {m}" for N0, (c, m) in sorted(res["by_N0"].items()))
    s += f"; non-chain prefixes alive at N_max: {res['alive']}"
    return s


def read_special(path):
    pairs = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("d = "):
                d = int(line.split(":")[0].split("=")[1])
                ks = eval(line.split("k in ")[1])
                pairs += [(d, k) for k in ks]
    return pairs


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "words":
        words_report()
    elif cmd == "remaining":
        remaining_report()
    elif cmd == "formal":
        formal_step(out_path=os.path.join(HERE, "out", "one_turn_special.txt"))
    elif cmd == "exact":
        d, k = int(sys.argv[2]), int(sys.argv[3])
        Nmax = int(sys.argv[4]) if len(sys.argv) > 4 else 20
        t0 = time.time()
        r = exact_turn(d, k, Nmax)
        print(fmt(r), f"({time.time() - t0:.1f} s; reach {r['reach']})")
    elif cmd == "exact_all":
        Nmax = int(sys.argv[2]) if len(sys.argv) > 2 else 20
        part, parts = (int(v) for v in sys.argv[3].split("/")) if len(sys.argv) > 3 else (0, 1)
        pairs = read_special(os.path.join(HERE, "out", "one_turn_special.txt"))
        pairs = [pq for i, pq in enumerate(pairs) if i % parts == part]
        t0 = time.time()
        bad = 0
        for (d, k) in pairs:
            r = exact_turn(d, k, Nmax)
            bad += r["alive"]
            print(fmt(r))
            sys.stdout.flush()
        print(f"part {part}/{parts}: {len(pairs)} turns searched to N_max = {Nmax}; non-chain prefixes alive at "
              f"N_max: {bad}  ({time.time() - t0:.0f} s)")
