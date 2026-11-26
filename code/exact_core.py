"""Exact arithmetic for bodies whose centres lie in Z[zeta], zeta = exp(2 pi i / Q), and the search for merge-closed
families in that arithmetic.  (proofs/fixed-side-count.md)

Why this is exact.
- A centre is an integer vector: its coordinates in the basis 1, zeta, ..., zeta^(phi(Q) - 1) of Z[zeta]. Sums,
  products with powers of zeta and complex conjugates are computed with integers, reduced by the cyclotomic polynomial.
- Q is even, so the roots of unity of the field Q(zeta) are exactly the Q powers of zeta. An element d of Z[zeta] has
  |d| = 1 exactly when it is a root of unity (Kronecker), that is, when it is one of
  these Q vectors. So "two cells touch" and "two cells coincide" are decided by comparing integer vectors.
- "Centres less than 1 apart" (the coin form of overlap) is decided with integer bounds on the real and imaginary
  parts of d. The enclosing trigonometric intervals are computed once per field. Any unresolved comparison is
  repeated at higher interval precision. Floating-point approximations never discard a candidate here.

Bodies. A body is a sorted tuple of centres. Two bodies of two or more cells are the same up to rotation and
translation exactly when their canonical forms agree; the canonical form is the least sorted tuple over all ways of
moving an ordered touching pair of cells to (0, 1). (Any rotation-and-translation carrying one connected body onto
another carries touching pairs to touching pairs, so this is a complete invariant; the rotation is a ratio of two
roots of unity.) The mirror image is the complex conjugate.

Overlap is the coin condition only: centres at least 1 apart. Regular n-gon cells with inradius 1/2 contain these
coins, so every n-gon body passes; a search under the coin condition keeps a superset of the n-gon families.
"""
from collections import defaultdict
from functools import cached_property

import numpy as np
from sympy import Poly, cyclotomic_poly, totient
from sympy.abc import x as _x


class _LazyR:
    """R[k] = matrix of multiplication by zeta^k (built on first use)."""

    def __init__(self, F):
        self.F = F
        self.cache = {}

    def __getitem__(self, k):
        F = self.F
        k %= F.Q
        m = self.cache.get(k)
        if m is None:
            m = np.array([F.pw[(j + k) % F.Q] for j in range(F.phi)], dtype=np.int64)
            self.cache[k] = m
        return m


class CycField:
    def __init__(self, Q):
        assert Q % 2 == 0 and Q >= 2
        self.Q = Q
        self.phi = int(totient(Q))
        co = [int(c) for c in Poly(cyclotomic_poly(Q, _x), _x).all_coeffs()][::-1]      # low -> high, monic
        phi = self.phi
        pw = []
        v = [0] * phi
        v[0] = 1
        for _ in range(Q):
            pw.append(tuple(v))
            carry = v[-1]
            v = [0] + v[:-1]
            if carry:
                for i in range(phi):
                    v[i] -= carry * co[i]
        assert tuple(v) == pw[0]
        self.pw = pw                                                   # pw[k] = zeta^k
        self.R = _LazyR(self)                                          # R[k]: row vector a @ R[k] = a * zeta^k
        self.C = np.array([pw[(-j) % Q] for j in range(phi)], dtype=np.int64)          # complex conjugation
        self.unit = {pw[k]: k for k in range(Q)}                       # root of unity -> exponent
        self.zero = (0,) * phi
        self.one = pw[0]
        self.hp_calls = 0

    # ---- elements (tuples of ints) ----
    def add(self, a, b):
        return tuple(p + q for p, q in zip(a, b))

    def sub(self, a, b):
        return tuple(p - q for p, q in zip(a, b))

    def mulz(self, a, k):
        """a * zeta^k."""
        return tuple((np.array(a, dtype=np.int64) @ self.R[k % self.Q]).tolist())

    def conj(self, a):
        return tuple((np.array(a, dtype=np.int64) @ self.C).tolist())

    def conj(self, a):
        return tuple((np.array(a, dtype=np.int64) @ self.C).tolist())

    def mul(self, a, b):
        """Product of two elements."""
        out = np.zeros(self.phi, dtype=np.int64)
        bb = np.array(b, dtype=np.int64)
        for j, c in enumerate(a):
            if c:
                out += c * (bb @ self.R[j])
        return tuple(out.tolist())

    def from_poly(self, d):
        """sum of c * zeta^e for a dict {e: c}."""
        v = [0] * self.phi
        for e, c in d.items():
            p = self.pw[e % self.Q]
            for i in range(self.phi):
                v[i] += c * p[i]
        return tuple(v)

    def cmp1(self, d):
        """Compare |d| with 1 for d in Z[zeta]: returns 0 (root of unity, |d| = 1), -1 (|d| < 1) or +1 (|d| > 1)."""
        if d in self.unit:
            return 0
        scale, powers = self._embedding_bounds
        lo_x = hi_x = lo_y = hi_y = 0
        for c, (a, b, u, v) in zip(d, powers):
            if not c:
                continue
            c = int(c)
            if c < 0:
                a, b, u, v = b, a, v, u
            lo_x += c * a
            hi_x += c * b
            lo_y += c * u
            hi_y += c * v
        x_sq = lo_x * lo_x, hi_x * hi_x
        y_sq = lo_y * lo_y, hi_y * hi_y
        if max(x_sq) + max(y_sq) < scale * scale:
            return -1
        low_x = 0 if lo_x <= 0 <= hi_x else min(x_sq)
        low_y = 0 if lo_y <= 0 <= hi_y else min(y_sq)
        if low_x + low_y > scale * scale:
            return 1
        return self._cmp1_hp(d)

    @cached_property
    def _embedding_bounds(self):
        """Outward integer enclosures for scale times each basis coordinate."""
        from mpmath import iv

        scale = 10**32
        previous = iv.dps
        try:
            iv.dps = 60
            powers = []
            for j in range(self.phi):
                angle = 2 * iv.pi * j / self.Q
                x, y = scale * iv.cos(angle), scale * iv.sin(angle)
                # int truncates toward zero, so add a full unit in each outward direction.
                powers.append((int(x.a) - 1, int(x.b) + 1, int(y.a) - 1, int(y.b) + 1))
            return scale, tuple(powers)
        finally:
            iv.dps = previous

    def _cmp1_hp(self, d):
        import mpmath
        self.hp_calls += 1
        for dps in (60, 200, 1000):
            mpmath.iv.dps = dps
            xr = mpmath.iv.mpf(0)
            yi = mpmath.iv.mpf(0)
            for j, c in enumerate(d):
                if c:
                    ang = 2 * mpmath.iv.pi * j / self.Q
                    xr += c * mpmath.iv.cos(ang)
                    yi += c * mpmath.iv.sin(ang)
            val = xr * xr + yi * yi - 1
            if val.b < 0:
                return -1
            if val.a > 0:
                return 1
        raise ArithmeticError("could not separate |d| from 1")


class Bodies:
    """Registry of bodies up to rotation and translation, in canonical position, with cached mirrors and splits."""

    def __init__(self, F):
        self.F = F
        self.ids = {}
        self.bodies = []       # id -> canonical sorted tuple of centres
        self.arr = []          # id -> numpy array
        self.pos = []          # id -> {centre: index}
        self.pairs = []        # id -> list of (i, j, k): ordered touching pairs, centre_j - centre_i = zeta^k
        self.mir = {}
        self.split_cache = {}
        self.place_cache = {}
        self.canon_cache = {}

    # ---- canonical form ----
    def touching_pairs(self, cells):
        F = self.F
        out = []
        N = len(cells)
        for i in range(N):
            ci = cells[i]
            for j in range(N):
                if i != j:
                    k = F.unit.get(tuple(p - q for p, q in zip(cells[j], ci)))
                    if k is not None:
                        out.append((i, j, k))
        return out

    def canon(self, cells):
        """Canonical sorted tuple of a set of centres (assumed connected by contacts)."""
        cells = tuple(sorted(cells))
        hit = self.canon_cache.get(cells)
        if hit is not None:
            return hit
        F = self.F
        N = len(cells)
        if N == 1:
            return (F.zero,)
        arr = np.array(cells, dtype=np.int64)
        best = None
        for (i, j, k) in self.touching_pairs(cells):
            img = (arr - arr[i]) @ F.R[(-k) % F.Q]
            key = tuple(sorted(map(tuple, img.tolist())))
            if best is None or key < best:
                best = key
        if best is None:
            raise ValueError("no touching pair")
        if len(self.canon_cache) < 300_000:
            self.canon_cache[cells] = best
        return best

    def add(self, cells):
        key = self.canon(cells)
        i = self.ids.get(key)
        if i is None:
            i = len(self.bodies)
            self.ids[key] = i
            self.bodies.append(key)
            self.arr.append(np.array(key, dtype=np.int64))
            self.pos.append({c: t for t, c in enumerate(key)})
            self.pairs.append(self.touching_pairs(key))
        return i

    def mirror(self, i):
        if i not in self.mir:
            F = self.F
            img = (self.arr[i] @ F.C).tolist()
            j = self.add([tuple(r) for r in img])
            self.mir[i] = j
            self.mir[j] = i
        return self.mir[i]

    def shape(self, i):
        return min(i, self.mirror(i))

    # ---- geometry ----
    def valid_new(self, cells, xnew):
        """Coin condition for a new centre against the centres of a body: no centre coincides or is closer than 1."""
        F = self.F
        for c in cells:
            d = tuple(p - q for p, q in zip(xnew, c))
            if d not in F.unit and F.cmp1(d) < 0:
                return False
        return True

    def contact_graph(self, i):
        nb = defaultdict(set)
        for (a, b, _) in self.pairs[i]:
            nb[a].add(b)
        return nb

    def is_chain(self, i):
        N = len(self.bodies[i])
        if N <= 2:
            return True
        nb = self.contact_graph(i)
        return sorted(len(nb[t]) for t in range(N)) == [1, 1] + [2] * (N - 2)

    def degrees(self, i):
        nb = self.contact_graph(i)
        return sorted(len(nb[t]) for t in range(len(self.bodies[i])))

    # ---- placements and splits ----
    def placements(self, a, t):
        """All sets of cells of body t that are rotated-and-translated copies of body a (as frozensets of indices)."""
        key = (a, t)
        hit = self.place_cache.get(key)
        if hit is not None:
            return hit
        F = self.F
        A = self.arr[a]
        T = self.bodies[t]
        pos = self.pos[t]
        out = set()
        if len(A) == 1:
            out = {frozenset([i]) for i in range(len(T))}
        elif len(A) <= len(T):
            # the canonical form of A has cells at 0 and 1
            for (i, j, k) in self.pairs[t]:
                img = (A @ F.R[k] + self.arr[t][i]).tolist()
                ids = []
                for r in img:
                    q = pos.get(tuple(r))
                    if q is None:
                        break
                    ids.append(q)
                else:
                    out.add(frozenset(ids))
        self.place_cache[key] = out
        return out

    def splits(self, t, a, b):
        """Is body t the disjoint union of a copy of a and a copy of b (rotations and translations)? The two parts
        touch because t is connected."""
        key = (t, a, b)
        hit = self.split_cache.get(key)
        if hit is not None:
            return hit
        T = self.bodies[t]
        res = False
        if len(self.bodies[a]) + len(self.bodies[b]) == len(T):
            bkey = self.bodies[b]
            for P in self.placements(a, t):
                rest = [T[i] for i in range(len(T)) if i not in P]
                if len(rest) == 1:
                    res = len(bkey) == 1
                else:
                    try:
                        res = self.canon(rest) == bkey
                    except ValueError:
                        res = False
                if res:
                    break
        self.split_cache[key] = res
        return res


def family_search(Rg, Nmax, candidates, start=None, on_prefix=None, allow=None, max_prefixes=None):
    """Depth-first search over all prefixes F(1..N), N <= Nmax, that satisfy the definition of a merge-closed family
    for every split a + b <= N. F[N] is a tuple of body ids: one body, or a body
    and its mirror image.

    candidates(F, N) -> iterable of body ids: the bodies T of N + 1 cells that may be members of F(N + 1) (any
        superset of the true candidates up to mirror images is allowed).
    start: a list of prefixes to start from (each a list F with F[0] = None); default: the single cell.
    allow(M, body id) -> bool: optional restriction on the members of F(M).
    on_prefix(F, N): called for every prefix reached.
    Returns (reach, finals): reach[N] = number of prefixes reaching N, finals = the prefixes reaching Nmax."""
    reach = defaultdict(int)
    finals = []
    count = [0]

    def dfs(F, N):
        reach[N] += 1
        count[0] += 1
        if max_prefixes is not None and count[0] > max_prefixes:
            raise RuntimeError("too many prefixes")
        if on_prefix is not None:
            on_prefix(F, N)
        if N == Nmax:
            finals.append(list(F))
            return
        M = N + 1
        cands = set()
        for T in candidates(F, N):
            cands.add(min(T, Rg.mirror(T)))
        for T in sorted(cands):
            if allow is not None and not allow(M, T):
                continue
            Tm = Rg.mirror(T)
            options = [(T,)] if T == Tm else [(T,), (Tm,), (T, Tm)]
            for opt in options:
                ok = True
                for a in range(M // 2, 0, -1):
                    for A in F[a]:
                        for B in F[M - a]:
                            if not any(Rg.splits(t, A, B) for t in opt):
                                ok = False
                                break
                        if not ok:
                            break
                    if not ok:
                        break
                if ok:
                    dfs(F + [opt], M)

    if start is None:
        start = [[None, (Rg.add([Rg.F.zero]),)]]
    for F0 in start:
        dfs(list(F0), len(F0) - 1)
    return reach, finals


# ---------------------------------------------------------------------------------------------------------------
# candidate generators
# ---------------------------------------------------------------------------------------------------------------

def ngon_candidates(Rg, n):
    """Candidates for regular n-gon cells: Q = 2n, zeta = exp(i pi / n). A cell whose contact with a neighbour is in
    direction zeta^m has edge normals zeta^(m + 2e). The new cell is placed across any free edge of any cell, and
    must satisfy the coin condition; if it lands at distance exactly 1 from another cell, the direction between them
    must be an edge normal of both (otherwise the two n-gons overlap)."""
    F = Rg.F
    assert F.Q == 2 * n
    cache = {}

    def parities(i):
        """Parity of the exponents of the edge normals of each cell (cell 0 of the canonical form has a neighbour
        at +1 = zeta^0, or is alone). None if the body is not consistent."""
        N = len(Rg.bodies[i])
        par = [None] * N
        nb = defaultdict(list)
        for (a, b, k) in Rg.pairs[i]:
            nb[a].append((b, k))
        z = Rg.pos[i][F.zero]
        par[z] = 0
        stack = [z]
        while stack:
            a = stack.pop()
            for (b, k) in nb[a]:
                if (k - par[a]) % 2:
                    return None
                pb = (k + n) % 2
                if par[b] is None:
                    par[b] = pb
                    stack.append(b)
                elif par[b] != pb:
                    return None
        return par

    def ext(i):
        if i in cache:
            return cache[i]
        cells = Rg.bodies[i]
        pos = Rg.pos[i]
        par = parities(i)
        out = set()
        if par is not None:
            for a, c in enumerate(cells):
                for e in range(n):
                    m = (par[a] + 2 * e) % F.Q
                    xn = F.add(c, F.pw[m])
                    if xn in pos:
                        continue
                    px = (m + n) % 2
                    good = True
                    for b, c2 in enumerate(cells):
                        d = F.sub(c2, xn)
                        k = F.unit.get(d)
                        if k is not None:
                            if (k - px) % 2 or (k + n - par[b]) % 2:
                                good = False
                                break
                        elif F.cmp1(d) < 0:
                            good = False
                            break
                    if good:
                        out.add(Rg.add(cells + (xn,)))
        cache[i] = sorted(out)
        return cache[i]

    def candidates(Fam, N):
        out = set()
        for A in Fam[N]:
            out.update(ext(A))
        return out

    return candidates


def lemma15_candidates(Rg):
    """Candidates by Lemma 15 of nonchain_attempt.md, for N + 1 >= 6, in exact arithmetic.

    S(N + 1) = Y + x with Y a member of F(N), and for the split (3, N - 2) the part P that contains x is h(S_p),
    p in {3, N - 2}, h an isometry, S_p a member of F(p). All cells of P except x are cells of Y. So x = h(c) for a
    cell c of S_p, where h carries the other cells of S_p into Y. h(w) = y1 + rho (w - a) or y1 + rho conj(w - a)
    with rho a root of unity of the field (see the note, Lemma 2.3). We fix two other cells a, b of S_p (a touching
    pair if there is one) and run over all ordered pairs (y1, y2) of cells of Y with y2 - y1 = rho (b - a), resp.
    rho conj(b - a)."""
    F = Rg.F
    Q = F.Q
    cache = {}

    def ext(Yid, Sid):
        key = (Yid, Sid)
        if key in cache:
            return cache[key]
        Y = Rg.bodies[Yid]
        Ypos = Rg.pos[Yid]
        Yarr = Rg.arr[Yid]
        S = Rg.bodies[Sid]
        p = len(S)
        diffs = defaultdict(list)                      # y2 - y1 -> [index of y1]
        for i in range(len(Y)):
            for j in range(len(Y)):
                if i != j:
                    diffs[tuple((Yarr[j] - Yarr[i]).tolist())].append(i)
        touching = defaultdict(set)
        for (a, b, _) in Rg.pairs[Sid]:
            touching[a].add(b)
        out = set()
        for flip in (False, True):
            Sarr = Rg.arr[Sid] @ F.C if flip else Rg.arr[Sid]
            for ci in range(p):
                others = [t for t in range(p) if t != ci]
                pair = None
                for a in others:
                    for b in others:
                        if b in touching[a]:
                            pair = (a, b)
                            break
                    if pair:
                        break
                if pair is None:
                    pair = (others[0], others[1])
                a, b = pair
                e = Sarr[b] - Sarr[a]
                rel = Sarr - Sarr[a]
                for k in range(Q):
                    d = tuple((e @ F.R[k]).tolist())
                    for i1 in diffs.get(d, ()):
                        img = (rel @ F.R[k] + Yarr[i1]).tolist()
                        xn = tuple(img[ci])
                        if xn in Ypos:
                            continue
                        if any(tuple(img[t]) not in Ypos for t in others):
                            continue
                        if Rg.valid_new(Y, xn):
                            out.add(Rg.add(Y + (xn,)))
        cache[key] = sorted(out)
        return cache[key]

    def candidates(Fam, N):
        M = N + 1
        assert M >= 6
        out = set()
        for Yid in Fam[N]:
            for p in {3, M - 3}:
                out.update(ext(Yid, Fam[p][0]))
        return out

    candidates.ext = ext
    return candidates
