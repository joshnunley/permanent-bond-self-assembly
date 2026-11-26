"""Special points of a box: the torsion points at which a relation logged by the formal search holds.
(proofs/late-cycles.md)

Each logged relation gives a Laurent polynomial G that vanishes wherever the relation holds:
    zero g:              G = g
    contact g:           G = g * conj(g) - 1
    ratio (d, u1, u2):   G = d * conj(d) - (u1 - u2) * conj(u1 - u2)
(conj replaces t by 1/t and conjugates the coefficients; it is complex conjugation on the torus.)

One free angle (r = 1). The special points are among the roots of unity t0 in the interval with G(t0) = 0. Their
order d satisfies phi(lcm(M, d)) <= phi(M) * span(G). They are screened in floating point and tested exactly in
Z[zeta_Q], Q = lcm(M, d, 2).

Two free angles (r = 2). A special point t either
  (1) lies on a torsion coset {t^lam = eta} (lam a primitive integer vector, eta a root of unity) on which G
      vanishes identically; these cosets are found exactly: lam is parallel to a difference of two exponents of G,
      and eta is a common root of the class polynomials (the sums of the terms of G whose exponents agree modulo
      lam, as polynomials in eta); or
  (2) satisfies t^delta, t^delta' in mu_m for two independent differences delta, delta' of exponents of the short
      relation, where m = lcm(M, product of the primes <= K) and K is the number of roots of unity in the short
      relation (Lemma 5.4, from Mann's theorem). These points are enumerated in the box, screened in floating
      point (margin 1e-7) and tested exactly: in Z[zeta_Q] when phi(Q) <= 600, and otherwise by exact division by
      the Q-th cyclotomic polynomial (vanishes_at).
"""
import math
from collections import defaultdict
from fractions import Fraction
from itertools import combinations

import numpy as np

from exact_core import CycField

_TOT = {}


def _totients(n):
    if n not in _TOT:
        phi = list(range(n + 1))
        for p in range(2, n + 1):
            if phi[p] == p:
                for k in range(p, n + 1, p):
                    phi[k] -= phi[k] // p
        _TOT[n] = phi
    return _TOT[n]


def totient(n):
    res = n
    p = 2
    m = n
    while p * p <= m:
        if m % p == 0:
            while m % p == 0:
                m //= p
            res -= res // p
        p += 1
    if m > 1:
        res -= res // m
    return res


_ORDERS = {}


def orders_up_to_degree(M, deg):
    """All d >= 1 with phi(lcm(M, d)) <= phi(M) * deg: the possible orders of a root of unity that is a root of a
    nonzero polynomial of degree <= deg over Q(zeta_M). (phi(L) >= sqrt(L / 2) bounds the search.)"""
    key = (M, deg)
    if key not in _ORDERS:
        bound = totient(M) * deg
        lim = 2 * bound * bound + 2
        phi = _totients(lim)
        out = []
        for d in range(1, lim + 1):
            L = d * M // math.gcd(d, M)
            # if L > lim then phi(L) >= sqrt(L / 2) > bound, so d is not wanted
            if L <= lim and phi[L] <= bound:
                out.append(d)
        _ORDERS[key] = out
    return _ORDERS[key]


_FIELDS = {}


def field(Q):
    if Q not in _FIELDS:
        _FIELDS[Q] = CycField(Q)
    return _FIELDS[Q]


def primorial(K):
    out = 1
    for p in range(2, K + 1):
        if all(p % q for q in range(2, int(p ** 0.5) + 1)):
            out *= p
    return out


# ---------------------------------------------------------------------------------------------------------------
# polynomials as dicts {exponent tuple: coefficient tuple in Z[zeta_M]}
# ---------------------------------------------------------------------------------------------------------------

class Polys:
    def __init__(self, ring):
        self.ring = ring
        self.K = ring.K
        self.r = ring.r

    def of(self, g):
        return dict(self.ring.monomials(g))

    def conj(self, P):
        return {tuple(-v for v in a): self.K.conj(c) for a, c in P.items()}

    def mul(self, P, Q):
        out = {}
        for a, c in P.items():
            for b, e in Q.items():
                key = tuple(x + y for x, y in zip(a, b))
                prod = self.K.mul(c, e)
                if key in out:
                    out[key] = self.K.add(out[key], prod)
                else:
                    out[key] = prod
        return {a: c for a, c in out.items() if any(c)}

    def sub(self, P, Q):
        out = dict(P)
        for a, c in Q.items():
            out[a] = self.K.sub(out[a], c) if a in out else tuple(-v for v in c)
        return {a: c for a, c in out.items() if any(c)}

    def norm(self, P):
        return self.mul(P, self.conj(P))

    def relation(self, typ, payload):
        """(G, short): G the polynomial that vanishes where the relation holds; short = (exponent differences of
        the short relation, number of roots of unity in it)."""
        one = {(0,) * self.r: self.K.one}
        if typ == "zero":
            g = self.of(payload)
            G = g
            exps = list(g)
            K = sum(sum(1 for v in c if v) for c in g.values())
            delta = {tuple(x - y for x, y in zip(a, b)) for a in exps for b in exps if a != b}
        elif typ == "contact":
            g = self.of(payload)
            G = self.sub(self.norm(g), one)
            exps = list(g)
            K = sum(sum(1 for v in c if v) for c in g.values()) + 1
            delta = {tuple(x - y for x, y in zip(a, b)) for a in exps for b in exps if a != b}
        else:
            d, u1, u2 = (self.of(v) for v in payload)
            e = self.sub(u1, u2)
            G = self.sub(self.norm(d), self.norm(e))
            exps = list(d)
            (b1,), (b2,) = list(u1), list(u2)
            eps = tuple(x - y for x, y in zip(b1, b2))
            K = sum(sum(1 for v in c if v) for c in d.values()) + 2
            dd = {tuple(x - y for x, y in zip(a, b)) for a in exps for b in exps}
            delta = set(dd) | {tuple(x + y for x, y in zip(eps, v)) for v in dd}
            delta.discard((0,) * self.r)
        return G, (delta, K)

    def float_terms(self, P):
        """(exponent array (n, r), complex coefficient array (n,))."""
        ex = np.array(list(P.keys()), dtype=float).reshape(len(P), self.r)
        co = np.array([np.dot(c, self.ring.emb) for c in P.values()])
        return ex, co


def specialize(M, P, F, logs_t):
    """Value in F = Q(zeta_Q) of the polynomial P at t_i = zeta_Q^(logs_t[i]). Q is a multiple of M."""
    tM = F.Q // M
    poly = defaultdict(int)
    for a, c in P.items():
        base = sum(x * y for x, y in zip(a, logs_t))
        for j, v in enumerate(c):
            if v:
                poly[(j * tM + base) % F.Q] += v
    return F.from_poly(poly)


_PHI = {}


def cyclotomic_coeffs(Q):
    """Coefficients of the Q-th cyclotomic polynomial, low degree first (Python integers)."""
    if Q not in _PHI:
        from sympy import Poly, cyclotomic_poly
        from sympy.abc import x
        _PHI[Q] = [int(c) for c in Poly(cyclotomic_poly(Q, x), x).all_coeffs()][::-1]
    return _PHI[Q]


def vanishes_at(M, P, fracs):
    """Exact test: is the polynomial P (dict {exponent tuple: coefficient in Z[zeta_M]}) zero at the torsion point
    t_i = exp(2 pi i fracs[i])? With Q = lcm(M, denominators), P(t) = A(zeta_Q) for an integer polynomial A of
    degree < Q, and A(zeta_Q) = 0 exactly when the Q-th cyclotomic polynomial divides A. The division is done with
    Python integers, so the test is exact for every Q (it is used when Q is too large for CycField)."""
    Q = M
    for f in fracs:
        Q = Q * f.denominator // math.gcd(Q, f.denominator)
    tM = Q // M
    lt = [int(f * Q) % Q for f in fracs]
    A = [0] * Q
    for a, c in P.items():
        base = sum(x * y for x, y in zip(a, lt))
        for j, v in enumerate(c):
            if v:
                A[(base + j * tM) % Q] += v
    phi = cyclotomic_coeffs(Q)
    n = len(phi) - 1                        # phi is monic of degree n
    for i in range(Q - 1, n - 1, -1):
        c = A[i]
        if c:
            for j in range(n + 1):
                A[i - n + j] -= c * phi[j]
    return not any(A[:n])


# ---------------------------------------------------------------------------------------------------------------
# r = 1
# ---------------------------------------------------------------------------------------------------------------

def find_special_1(ring, logs, lo, hi):
    """All (d, k), gcd(k, d) = 1, 0 <= k < d, with lo <= 2 pi k / d <= hi modulo 2 pi (up to 1e-9), at which the
    polynomial G of some logged relation vanishes."""
    M = ring.M
    PS = Polys(ring)
    special = set()
    for (typ, _), payload in logs.items():
        G, _ = PS.relation(typ, payload)
        assert G, "a logged relation holds identically"
        exps = [a[0] for a in G]
        span = max(max(exps) - min(exps), 1)
        ex, co = PS.float_terms(G)
        cands = []
        for d in orders_up_to_degree(M, span):
            k0 = math.ceil((lo - 1e-9) * d / (2 * math.pi))
            k1 = math.floor((hi + 1e-9) * d / (2 * math.pi))
            for k in range(k0, k1 + 1):
                if math.gcd(k, d) == 1:
                    cands.append((d, k % d))
        if not cands:
            continue
        th = np.array([2 * math.pi * k / d for d, k in cands])
        vals = np.abs(np.exp(1j * np.outer(th, ex[:, 0])) @ co)
        for idx in np.flatnonzero(vals < 1e-7):
            d, k = cands[idx]
            if (d, k) in special:
                continue
            Q = d * M // math.gcd(d, M)
            if Q % 2:
                Q *= 2
            F = field(Q)
            if not any(specialize(M, G, F, (k * (Q // d),))):
                special.add((d, k))
    return special


# ---------------------------------------------------------------------------------------------------------------
# r = 2
# ---------------------------------------------------------------------------------------------------------------

def _ext_gcd(a, b):
    if b == 0:
        return (1, 0) if a > 0 else (-1, 0)
    x, y = _ext_gcd(b, a % b)
    return (y, x - (a // b) * y)


def coset_basis(a, b):
    """(c, d) with a d - b c = 1, for a primitive vector (a, b)."""
    x, y = _ext_gcd(a, b)           # a x + b y = +-1 * gcd
    g = a * x + b * y
    assert abs(g) == 1
    x, y = x * g, y * g
    c, d = -y, x
    assert a * d - b * c == 1
    return c, d


def primitive(v):
    g = math.gcd(abs(v[0]), abs(v[1]))
    a, b = v[0] // g, v[1] // g
    if a < 0 or (a == 0 and b < 0):
        a, b = -a, -b
    return (a, b)


def hnf2(v1, v2):
    """Hermite normal form (a, b, d) of the lattice spanned by two independent integer vectors: basis (a, b), (0, d)
    with a, d > 0 and 0 <= b < d."""
    (a, b), (c, d) = v1, v2
    while c != 0:
        q = a // c
        a, b, c, d = c, d, a - q * c, b - q * d
    if a < 0:
        a, b = -a, -b
    d = abs(d)
    return (a, b % d, d)


def find_special_2(ring, logs, lo, hi, max_phi=600, stats=None):
    """Special cosets and special points of the box [lo, hi] (radians) for two free angles.
    Returns (cosets, points): cosets = {(a, b, d, k)}: t1^a t2^b = exp(2 pi i k / d), (a, b) primitive;
    points = {(f1, f2)}: theta_i = 2 pi f_i, f_i Fractions in [0, 1)."""
    M = ring.M
    PS = Polys(ring)
    lo = np.array(lo)
    hi = np.array(hi)
    cosets = set()
    points = set()
    corners = np.array([[lo[0], lo[1]], [lo[0], hi[1]], [hi[0], lo[1]], [hi[0], hi[1]]])
    for (typ, _), payload in logs.items():
        G, (delta, K) = PS.relation(typ, payload)
        assert G, "a logged relation holds identically"
        # ---- (1) codimension 1: cosets on which G vanishes identically ----
        mine = set()
        Gex = list(G)
        dirs = {primitive(tuple(x - y for x, y in zip(p, q))) for p, q in combinations(Gex, 2)}
        for (a, b) in dirs:
            c, d_ = coset_basis(a, b)
            classes = defaultdict(list)
            for al in Gex:
                classes[a * al[1] - b * al[0]].append((d_ * al[0] - c * al[1], G[al]))
            if any(len(v) < 2 for v in classes.values()):
                continue
            span = min(max(k for k, _ in v) - min(k for k, _ in v) for v in classes.values())
            rng = corners @ np.array([a, b], dtype=float)
            l0, h0 = rng.min(), rng.max()
            cands = []
            for d in orders_up_to_degree(M, span):
                k0 = math.ceil((l0 - 1e-9) * d / (2 * math.pi))
                k1 = math.floor((h0 + 1e-9) * d / (2 * math.pi))
                for k in range(k0, k1 + 1):
                    if math.gcd(k, d) == 1:
                        cands.append((d, k % d))
            if not cands:
                continue
            th = np.array([2 * math.pi * k / d for d, k in cands])
            ok = np.ones(len(cands), dtype=bool)
            for v in classes.values():
                ex = np.array([k for k, _ in v], dtype=float)
                co = np.array([np.dot(cf, ring.emb) for _, cf in v])
                ok &= np.abs(np.exp(1j * np.outer(th, ex)) @ co) < 1e-7
            for idx in np.flatnonzero(ok):
                d, k = cands[idx]
                Q = d * M // math.gcd(d, M)
                if Q % 2:
                    Q *= 2
                F = field(Q)
                good = True
                for v in classes.values():
                    P = {(kk,): cf for kk, cf in v}
                    if any(specialize(M, P, F, (k * (Q // d),))):
                        good = False
                        break
                if good:
                    mine.add((a, b, d, k))
        cosets |= mine
        # ---- (2) isolated points ----
        m = M * primorial(K) // math.gcd(M, primorial(K))
        if stats is not None:
            stats["K"][K] += 1
        exG, coG = PS.float_terms(G)
        seen_lat = set()
        dl = sorted(delta)
        for v1, v2 in combinations(dl, 2):
            det = v1[0] * v2[1] - v1[1] * v2[0]
            if det == 0:
                continue
            key = hnf2(v1, v2)                       # visit each lattice once
            if key in seen_lat:
                continue
            seen_lat.add(key)
            r1 = corners @ np.array(v1, dtype=float)
            r2 = corners @ np.array(v2, dtype=float)
            p1 = np.arange(math.ceil((r1.min() - 1e-9) * m / (2 * math.pi)),
                           math.floor((r1.max() + 1e-9) * m / (2 * math.pi)) + 1)
            p2 = np.arange(math.ceil((r2.min() - 1e-9) * m / (2 * math.pi)),
                           math.floor((r2.max() + 1e-9) * m / (2 * math.pi)) + 1)
            if len(p1) == 0 or len(p2) == 0:
                continue
            P1, P2 = np.meshgrid(p1, p2, indexing="ij")
            P1 = P1.ravel()
            P2 = P2.ravel()
            # theta / 2 pi = adj (p1, p2) / (m det)
            n1 = v2[1] * P1 - v1[1] * P2
            n2 = -v2[0] * P1 + v1[0] * P2
            den = m * det
            t1 = 2 * math.pi * n1 / den
            t2 = 2 * math.pi * n2 / den
            inbox = (t1 >= lo[0] - 1e-9) & (t1 <= hi[0] + 1e-9) & (t2 >= lo[1] - 1e-9) & (t2 <= hi[1] + 1e-9)
            if not inbox.any():
                continue
            n1, n2, t1, t2 = n1[inbox], n2[inbox], t1[inbox], t2[inbox]
            vals = np.abs(np.exp(1j * (np.outer(t1, exG[:, 0]) + np.outer(t2, exG[:, 1]))) @ coG)
            for idx in np.flatnonzero(vals < 1e-7):
                f1 = Fraction(int(n1[idx]), den) % 1
                f2 = Fraction(int(n2[idx]), den) % 1
                if (f1, f2) in points:
                    continue
                # on one of the cosets found for this relation?
                if any((a * f1 + b * f2 - Fraction(k, d)).denominator == 1 for (a, b, d, k) in mine):
                    continue
                Q = M
                for f in (f1, f2):
                    Q = Q * f.denominator // math.gcd(Q, f.denominator)
                if Q % 2:
                    Q *= 2
                if totient(Q) <= max_phi:
                    F = field(Q)
                    lt = (int(f1 * Q), int(f2 * Q))
                    if any(specialize(M, G, F, lt)):
                        continue
                elif not vanishes_at(M, G, (f1, f2)):
                    continue
                points.add((f1, f2))
    return cosets, points
