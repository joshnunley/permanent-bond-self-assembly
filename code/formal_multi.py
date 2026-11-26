"""Bodies with r free angles (r = 1 or 2): centres are Laurent polynomials in t_1..t_r, t_i = exp(i theta_i), with
coefficients in Z[zeta_M]; theta ranges over a box B.  (proofs/late-cycles.md)

The formal search. Equalities are decided formally (identities of Laurent polynomials):
  - two cells coincide when their polynomials are equal;
  - two cells touch when their difference is a unit zeta^j t^k of the ring ("formal contact");
  - rotations are multiplications by such units.
The only inequality, "two centres are less than 1 apart", is decided on the whole box by the value at the midpoint
and a Lipschitz bound (|d/dtheta_i of sum c_e t^e| <= sum |e_i| |c_e|): certified "< 1 on all of B" kills a
candidate; anything else keeps it. So for each theta in B the formal search keeps a superset of the prefixes of the
true search at theta, PROVIDED no relevant equality holds at theta that does not hold formally (Lemma 5.2 of the
note). Every place where such a non-formal equality could matter is logged, as a relation that holds at theta:
  - "contact" g:  g is a difference of two cells of one body, not a unit, and |g| = 1 is possible on B;
                  relation: g = omega for some root of unity omega;
  - "zero" g:     g = p - q, p an image under a formal rotation-and-translation and q a cell it is looked up
                  against, p != q formally, p = q possible on B; relation: g = 0;
  - "ratio" (d, u1, u2): in the Lemma 15 generator, d = y2 - y1 for two cells of Y, u1 - u2 the difference of the
                  two ends of S(3) seen from its middle cell (u1, u2 units), d not a unit multiple of u1 - u2,
                  |d| = |u1 - u2| possible on B; relation: d = rho (u1 - u2) for some root of unity rho.
The special points of B are the torsion points (all theta_i rational multiples of pi) at which a logged relation
holds. They are found by special.py; each is covered by a search of lower dimension.

Elements are numpy integer arrays of shape (W,)*r + (phi,): entry [e_1 + E, .., e_r + E] is the coefficient of
t^e as a vector in the basis 1, zeta_M, .., zeta_M^(phi-1). The window E must be large enough; this is asserted.
"""
import cmath
import math
from collections import defaultdict

import numpy as np

from exact_core import CycField

TOL = 1e-9


class Ring:
    def __init__(self, M, E, r=1):
        self.K = CycField(M)
        self.M = M
        self.E = E
        self.r = r
        self.W = 2 * E + 1
        self.phi = self.K.phi
        self.shape = (self.W,) * r + (self.phi,)
        self.axes = tuple(range(-r - 1, -1))
        self.exps = np.arange(-E, E + 1)
        self.emb = np.array([cmath.exp(2j * math.pi * j / M) for j in range(self.phi)])
        self.absgrid = []
        for i in range(r):
            sh = [1] * r
            sh[i] = self.W
            self.absgrid.append(np.abs(self.exps).astype(float).reshape(sh))
        self.set_box([0.0] * r, [0.0] * r)

    def set_box(self, lo, hi):
        self.lo = np.array(lo, dtype=float)
        self.hi = np.array(hi, dtype=float)
        self.mid = (self.lo + self.hi) / 2
        self.h = (self.hi - self.lo) / 2
        self.ph = [np.exp(1j * self.exps * self.mid[i]) for i in range(self.r)]

    def zero(self):
        return np.zeros(self.shape, dtype=np.int64)

    def elem(self, poly):
        """sum of c * zeta_M^j * t^k for poly = {(k_1, .., k_r, j): c}."""
        a = self.zero()
        for key, c in poly.items():
            k, j = key[:-1], key[-1]
            assert len(k) == self.r and all(abs(v) <= self.E for v in k), "exponent window too small"
            a[tuple(v + self.E for v in k)] += c * np.array(self.K.pw[j % self.M], dtype=np.int64)
        return a

    def support(self, g):
        """Exponent vectors with nonzero coefficient, as an array (n, r) in lexicographic order."""
        return np.argwhere(g.any(axis=-1)) - self.E

    def monomials(self, g):
        return [(tuple(int(v) for v in a), tuple(g[tuple(a + self.E)].tolist())) for a in self.support(g)]

    def unit_of(self, g):
        """(k, j) if g = zeta^j t^k (k a tuple), else None."""
        nz = np.argwhere(g.any(axis=-1))
        if len(nz) != 1:
            return None
        j = self.K.unit.get(tuple(g[tuple(nz[0])].tolist()))
        if j is None:
            return None
        return (tuple(int(v) - self.E for v in nz[0]), j)

    def mul_unit(self, A, k, j):
        """A * zeta^j t^k for an array A of shape (...,) + self.shape."""
        if any(k):
            for i in range(self.r):
                if k[i]:
                    other = tuple(a for a in range(A.ndim) if a != A.ndim - 1 - self.r + i)
                    sup = np.flatnonzero(A.any(axis=other))
                    if len(sup):
                        assert 0 <= sup[0] + k[i] and sup[-1] + k[i] < self.W, "exponent window too small"
            A = np.roll(A, k, axis=self.axes)
        if j % self.M:
            A = A @ self.K.R[j]
        return A

    def conj(self, A):
        return np.flip(A, axis=self.axes) @ self.K.C

    def vals(self, A):
        V = A @ self.emb
        for i in reversed(range(self.r)):
            V = V @ self.ph[i]
        return V

    def rads(self, A):
        """Radius of the disc around the midpoint value that contains the values on the whole box."""
        absV = np.abs(A @ self.emb)
        ax = tuple(range(-self.r, 0))
        out = 0.0
        for i in range(self.r):
            out = out + (absV * self.absgrid[i]).sum(axis=ax) * self.h[i]
        return out

    def unit_ratio(self, d, e):
        """(k, j) with d = zeta^j t^k e, or None."""
        sd = np.argwhere(d.any(axis=-1))
        se = np.argwhere(e.any(axis=-1))
        if len(sd) != len(se) or len(sd) == 0:
            return None
        k = sd[0] - se[0]
        if not np.array_equal(sd, se + k):
            return None
        cd = d[tuple(sd.T)]
        ce = e[tuple(se.T)]
        for j in range(self.M):
            if np.array_equal(cd, ce @ self.K.R[j]):
                return (tuple(int(v) for v in k), j)
        return None

    def shift_norm(self, g):
        """g times the monomial that moves its lexicographically lowest exponent to 0."""
        nz = np.argwhere(g.any(axis=-1))
        shift = self.E - nz[0]
        assert ((nz + shift >= 0) & (nz + shift < self.W)).all(), "exponent window too small"
        return np.roll(g, tuple(int(v) for v in shift), axis=tuple(range(self.r)))


class FormalBodies:
    """The registry of exact_core.Bodies, for formal bodies on the box set in the ring."""

    def __init__(self, ring):
        self.ring = ring
        self.ids = {}
        self.arr = []
        self.pos = []
        self.pairs = []        # (i, j, k, jj): cell_j - cell_i = zeta^jj t^k
        self.vals = []
        self.rads = []
        self.mir = {}
        self.place_cache = {}
        self.split_cache = {}
        self.logs = {}

    def log(self, typ, g, u1=None, u2=None):
        r = self.ring
        if typ == "ratio":
            self.logs.setdefault((typ, g.tobytes() + u1.tobytes() + u2.tobytes()), (g, u1, u2))
        else:
            g = r.shift_norm(g)
            self.logs.setdefault((typ, g.tobytes()), g)

    def touching_pairs(self, arr, vals, rads):
        r = self.ring
        m = np.abs(vals[None, :] - vals[:, None])
        thr = rads[None, :] + rads[:, None] + TOL
        out = []
        for i, j in np.argwhere(np.abs(m - 1) <= thr):
            if i == j:
                continue
            g = arr[j] - arr[i]
            u = r.unit_of(g)
            if u is not None:
                out.append((int(i), int(j), u[0], u[1]))
            elif abs(m[i, j] - 1) <= r.rads(g) + TOL:
                self.log("contact", g)
        return out

    def canon(self, arr):
        r = self.ring
        N = len(arr)
        if N == 1:
            z = np.zeros_like(arr)
            return (z[0].tobytes(),), z
        vals = r.vals(arr)
        rads = r.rads(arr)
        best = None
        bestimg = None
        for (i, j, k, jj) in self.touching_pairs(arr, vals, rads):
            img = r.mul_unit(arr - arr[i], tuple(-v for v in k), -jj)
            keys = [row.tobytes() for row in img]
            order = sorted(range(N), key=keys.__getitem__)
            key = tuple(keys[t] for t in order)
            if best is None or key < best:
                best = key
                bestimg = img[order]
        if best is None:
            raise ValueError("no touching pair")
        return best, bestimg

    def add(self, cells):
        arr = np.array(cells, dtype=np.int64) if not isinstance(cells, np.ndarray) else cells
        key, img = self.canon(arr)
        i = self.ids.get(key)
        if i is None:
            r = self.ring
            i = len(self.arr)
            self.ids[key] = i
            self.arr.append(img)
            self.pos.append({k: t for t, k in enumerate(key)})
            v = r.vals(img)
            ra = r.rads(img)
            self.vals.append(v)
            self.rads.append(ra)
            self.pairs.append(self.touching_pairs(img, v, ra))
        return i

    def mirror(self, i):
        if i not in self.mir:
            j = self.add(self.ring.conj(self.arr[i]))
            self.mir[i] = j
            self.mir[j] = i
        return self.mir[i]

    def shape(self, i):
        return min(i, self.mirror(i))

    def contact_graph(self, i):
        nb = defaultdict(set)
        for (a, b, _, _) in self.pairs[i]:
            nb[a].add(b)
        return nb

    def is_chain(self, i):
        N = len(self.arr[i])
        if N <= 2:
            return True
        nb = self.contact_graph(i)
        return sorted(len(nb[t]) for t in range(N)) == [1, 1] + [2] * (N - 2)

    def degrees(self, i):
        nb = self.contact_graph(i)
        return sorted(len(nb[t]) for t in range(len(self.arr[i])))

    def valid_new_arr(self, arr, x):
        """False only if x is certainly (for every theta in the box) closer than 1 to some cell of arr."""
        r = self.ring
        diffs = x[None] - arr
        v = np.abs(r.vals(diffs))
        L = r.rads(diffs) + TOL
        if np.any(v + L < 1):
            return False
        for idx in np.flatnonzero(v - L <= 1):
            g = diffs[idx]
            if r.unit_of(g) is None:
                self.log("contact", g)
        return True

    def valid_set(self, cells):
        arr = np.array(cells, dtype=np.int64)
        if len({c.tobytes() for c in arr}) < len(arr):
            return False
        return all(self.valid_new_arr(arr[:i], arr[i]) for i in range(1, len(arr)))

    def lookup(self, img, t):
        """Indices in body t of the cells img (formal equality), or None. Logs near misses."""
        pos = self.pos[t]
        ids = [pos.get(row.tobytes()) for row in img]
        if None not in ids:
            return ids
        r = self.ring
        miss = [a for a, q in enumerate(ids) if q is None]
        sub = img[miss]
        D = np.abs(r.vals(sub)[:, None] - self.vals[t][None, :])
        thr = r.rads(sub)[:, None] + self.rads[t][None, :] + TOL
        for a, b in np.argwhere(D <= thr):
            g = sub[a] - self.arr[t][b]
            if D[a, b] <= r.rads(g) + TOL:
                self.log("zero", g)
        return None

    def placements(self, a, t):
        key = (a, t)
        hit = self.place_cache.get(key)
        if hit is not None:
            return hit
        r = self.ring
        A = self.arr[a]
        T = self.arr[t]
        out = set()
        if len(A) == 1:
            out = {frozenset([i]) for i in range(len(T))}
        elif len(A) <= len(T):
            for (i, j, k, jj) in self.pairs[t]:
                ids = self.lookup(r.mul_unit(A, k, jj) + T[i], t)
                if ids is not None:
                    out.add(frozenset(ids))
        self.place_cache[key] = out
        return out

    def rest_is(self, t, rest, b):
        B = self.arr[b]
        if len(B) != len(rest):
            return False
        if len(B) == 1:
            return True
        r = self.ring
        T = self.arr[t]
        for (i, j, k, jj) in self.pairs[t]:
            if i in rest and j in rest:
                ids = self.lookup(r.mul_unit(B, k, jj) + T[i], t)
                if ids is not None and set(ids) == rest:
                    return True
        return False

    def splits(self, t, a, b):
        key = (t, a, b)
        hit = self.split_cache.get(key)
        if hit is not None:
            return hit
        res = False
        nT = len(self.arr[t])
        if len(self.arr[a]) + len(self.arr[b]) == nT:
            full = set(range(nT))
            for P in self.placements(a, t):
                if self.rest_is(t, full - P, b):
                    res = True
                    break
        self.split_cache[key] = res
        return res


def formal_lemma15(Rg):
    """The Lemma 15 candidates (see exact_core.lemma15_candidates), formally. Rotations are units of the ring."""
    r = Rg.ring
    cache = {}

    def ext(Yid, Sid):
        key = (Yid, Sid)
        if key in cache:
            return cache[key]
        Y = Rg.arr[Yid]
        Ypos = Rg.pos[Yid]
        S = Rg.arr[Sid]
        p = len(S)
        NY = len(Y)
        out = set()
        tpairs = {}
        for (a, b, k, j) in Rg.pairs[Sid]:
            tpairs[(a, b)] = (k, j)

        def process(img, ci, others):
            x = img[ci]
            if x.tobytes() in Ypos:
                return
            if Rg.lookup(img[others], Yid) is None:
                return
            if Rg.valid_new_arr(Y, x):
                out.add(Rg.add(np.concatenate([Y, x[None]])))

        for flip in (False, True):
            Sarr = r.conj(S) if flip else S
            for ci in range(p):
                others = [t for t in range(p) if t != ci]
                pair = None
                for a in others:
                    for b in others:
                        if (a, b) in tpairs:
                            pair = (a, b)
                            break
                    if pair:
                        break
                if pair is not None:
                    a, b = pair
                    ke, je = tpairs[pair]
                    if flip:
                        ke, je = tuple(-v for v in ke), -je
                    rel = Sarr - Sarr[a]
                    for (i1, i2, kd, jd) in Rg.pairs[Yid]:
                        process(r.mul_unit(rel, tuple(x - y for x, y in zip(kd, ke)), jd - je) + Y[i1], ci, others)
                else:
                    # S is a chain of 3 cells and ci is its middle cell
                    a, b = others[0], others[1]
                    u1 = Sarr[b] - Sarr[ci]
                    u2 = Sarr[a] - Sarr[ci]
                    e = u1 - u2
                    rel = Sarr - Sarr[a]
                    ve = abs(r.vals(e))
                    re_ = r.rads(e)
                    yv = Rg.vals[Yid]
                    yr = Rg.rads[Yid]
                    for i1 in range(NY):
                        for i2 in range(NY):
                            if i1 == i2:
                                continue
                            if abs(abs(yv[i2] - yv[i1]) - ve) > yr[i1] + yr[i2] + re_ + TOL:
                                continue
                            d = Y[i2] - Y[i1]
                            u = r.unit_ratio(d, e)
                            if u is not None:
                                process(r.mul_unit(rel, u[0], u[1]) + Y[i1], ci, others)
                            elif abs(abs(r.vals(d)) - ve) <= r.rads(d) + re_ + TOL:
                                Rg.log("ratio", d, u1, u2)
        cache[key] = sorted(out)
        return cache[key]

    def candidates(Fam, N):
        M = N + 1
        out = set()
        for Yid in Fam[N]:
            for p in {3, M - 3}:
                if p >= 3:
                    out.update(ext(Yid, Fam[p][0]))
        return out

    candidates.ext = ext
    return candidates
