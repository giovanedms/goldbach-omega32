"""Rigorous (ball arithmetic, python-flint arb) data for T3: odd N, Omega(ab) <= 32, N = sqrt(8) z^16.

Every number that enters a certificate is either an exact integer/dyadic, or an arb ball, or a float64
that is a *proven* upper (hi) / lower (lo) bound of an arb ball (converted by hi()/lo() below, which
check the direction with an exact arb comparison).

Per-prime region (p < P0 = 1e6): exact per-prime ball sums.
Block region (P0 <= p < XMAX): blocks [A,B) with exact (count n, sum S1) from primes_blocks.py; for a convex
f on [A,B] (f = ell or e_delta, convexity proven in RIGOROSO.md sec. 2):
    n f(S1/n)  <=  sum_{p in [A,B)} f(p)  <=  n f(A) + (f(B)-f(A)) (S1 - nA)/(B-A)      (Jensen / chord)

Notation (P = primes not in Q):
  ell_p = log 2 (p=2), -log(1-2/p) (p odd)   ->  L(x) = sum_{p<x, p in P} ell_p,  V(z) >= exp(-L(z))
  e_d(p) = log(1 + 8 p^-d)                    ->  E_d(x) = sum_{p<x, p in P} e_d(p)
  log K(z) = max(0, sup_{w<z} [L(z) - L(w) - 2 log log z + 2 log log w])   (L(w) := sum_{p<w})
"""
import math, os
import numpy as np
from flint import arb, ctx

ctx.prec = 128
HERE = os.path.dirname(os.path.abspath(__file__))
P0 = 10**6
N0 = 117465180365547648498934439
DELTAS_NUM = list(range(26, 128))                 # delta = j/128 (dyadic, exact)
DELTAS = [j / 128 for j in DELTAS_NUM]
CQ = {(3,): arb(19) / 6, (3, 5): arb(11)}         # |r_d| <= C_Q 2^omega(d)   (paper Lemma 4.3 / RIGOROSO.md sec. 3)


def hi(x):
    """float64 f with f >= every point of ball x (checked)."""
    f = math.nextafter(float(x.upper()), math.inf)
    while not (arb(f) >= x):
        f = math.nextafter(f, math.inf)
    return f


def lo(x):
    f = math.nextafter(float(x.lower()), -math.inf)
    while not (arb(f) <= x):
        f = math.nextafter(f, -math.inf)
    return f


def ell(p):
    p = arb(p)
    return -(1 - 2 / p).log()


def loglog(x):
    return arb(x).log().log()


def logF(s, k):
    """log F(s,k), F = 1 - (s+3)/(2 e^k) (2 e k/(s-3))^((s-3)/2); s,k arb. Returns None if not certified s>=2k+3, F>0."""
    if not (s >= 2 * k + 3):
        return None
    t = (s + 3) / (2 * k.exp()) * (((s - 3) / 2) * (2 * arb.const_e() * k / (s - 3)).log()).exp()
    F = 1 - t
    if not (F > 0):
        return None
    return F.log()


def logF_np(s, k):
    with np.errstate(all="ignore"):
        F = 1 - (s + 3) / (2 * np.exp(k)) * (2 * math.e * k / (s - 3)) ** ((s - 3) / 2)
        return np.where((s >= 2 * k + 3) & (F > 0), np.log(F), -np.inf)


def z0_lower():
    """u0 = float <= z0* = (N0^2/8)^(1/32), checked exactly: 8 u0^32 <= N0^2."""
    from fractions import Fraction
    x = ((arb(N0) ** 2 / 8).log() / 32).exp()
    u0 = lo(x)
    assert 8 * Fraction(u0) ** 32 <= Fraction(N0) ** 2
    return u0


class RData:
    def __init__(self, Q, xmax, verbose=True):
        self.Q = tuple(sorted(Q)); self.xmax = xmax
        from primes_blocks import primes_below
        pa = primes_below(P0)
        p = pa[~np.isin(pa, np.array(self.Q))]
        self.p = p; n = len(p)
        # ---- per-prime region ------------------------------------------------------------------
        Lc = [arb(0)]
        for q in p:
            Lc.append(Lc[-1] + (arb(2).log() if q == 2 else ell(int(q))))
        self.Lc_end = Lc[-1]
        self.Lhi_lt = np.array([hi(x) for x in Lc])           # index i: sum over first i primes of P
        self.Llo_lt = np.array([lo(x) for x in Lc])
        ll = [loglog(int(q)) for q in p]
        Hhi = np.array([hi(2 * ll[i] - arb(self.Llo_lt[i])) for i in range(n)])
        self.Hpm_hi = np.maximum.accumulate(Hhi)              # bound for sup over w <= p_i
        self.Kq_hi = np.array([hi(arb(self.Lhi_lt[i + 1]) - 2 * ll[i] + arb(self.Hpm_hi[i])) for i in range(n)])
        if verbose: print(f"# per-prime: {n} primes in P below {P0}", flush=True)
        self.Ehi_lt = np.zeros((len(DELTAS), n + 1))
        self.E_end = []
        for jd, dn in enumerate(DELTAS_NUM):
            d = arb(dn) / 128
            acc = arb(0); row = [hi(acc)]
            for q in p:
                acc = acc + (1 + 8 * arb(int(q)) ** (-d)).log()
                row.append(hi(acc))
            self.Ehi_lt[jd] = row; self.E_end.append(acc)
        if verbose: print("# per-prime E sums done", flush=True)
        # ---- block region -------------------------------------------------------------------------
        bl = np.load(os.path.join(HERE, "blocks.npz"))
        b, bn, bs = bl["bounds"], bl["n"], bl["s1"]
        k = int(np.searchsorted(b, xmax)); assert b[k] == xmax and b[0] == P0
        b, bn, bs = b[:k + 1], bn[:k], bs[:k]
        self.b = b
        nb = len(bn)
        fl = [ell(int(x)) for x in b]
        LBh, LBl = [self.Lc_end], [self.Lc_end]
        HB = []
        for j in range(nb):
            A, B, m, S = int(b[j]), int(b[j + 1]), int(bn[j]), int(bs[j])
            if m:
                up = m * fl[j] + (fl[j + 1] - fl[j]) * (arb(S - m * A) / (B - A))
                dn_ = m * ell(arb(S) / m)
                HB.append(hi(2 * loglog(B) - arb(lo(LBl[-1]))))
            else:
                up = dn_ = arb(0); HB.append(-math.inf)
            LBh.append(LBh[-1] + up); LBl.append(LBl[-1] + dn_)
        self.LBhi = np.array([hi(x) for x in LBh]); self.LBlo = np.array([lo(x) for x in LBl])
        self.HBpm = np.maximum(np.maximum.accumulate(np.array(HB)), self.Hpm_hi[-1])   # sup_w for w < b_{j+1}
        self.EBhi = np.zeros((len(DELTAS), nb + 1))
        for jd, dn in enumerate(DELTAS_NUM):
            d = arb(dn) / 128
            fe = [(1 + 8 * arb(int(x)) ** (-d)).log() for x in b]
            acc = self.E_end[jd]; row = [hi(acc)]
            for j in range(nb):
                A, B, m, S = int(b[j]), int(b[j + 1]), int(bn[j]), int(bs[j])
                if m:
                    acc = acc + m * fe[j] + (fe[j + 1] - fe[j]) * (arb(S - m * A) / (B - A))
                row.append(hi(acc))
            self.EBhi[jd] = row
        if verbose: print(f"# blocks: {nb} up to {xmax:.3g}", flush=True)

    # ---- quantities at a point Y that is a block boundary (for the analytic range) ----------------
    def at_boundary(self, Y):
        j = int(np.searchsorted(self.b, Y)); assert self.b[j] == Y
        L_hi = self.LBhi[j]
        logKex_hi = max(0.0, hi(arb(L_hi) - 2 * loglog(Y) + arb(self.HBpm[j - 1])))
        return L_hi, logKex_hi

    # ---- cells ------------------------------------------------------------------------------------
    def cells(self, u0, eps1=1e-4):
        """returns dict of arrays: u, v (floats; exact dyadics), Lhi(v), Ehi(v,delta) [cells x deltas], logK_hi."""
        G = int(math.ceil(math.log(P0 / u0) / math.log1p(eps1)))
        z = u0 * (1 + eps1) ** np.arange(G + 1); z[-1] = float(P0)
        z = np.unique(z)
        u1, v1 = z[:-1], z[1:]
        a = np.searchsorted(self.p, u1, side="right")      # #P-primes <= u
        bb = np.searchsorted(self.p, v1, side="left")      # #P-primes <  v
        L1 = self.Lhi_lt[bb]; E1 = self.Ehi_lt[:, bb].T
        K1 = np.zeros(len(u1))
        for i in range(len(u1)):
            c = hi(arb(self.Lhi_lt[a[i]]) - 2 * loglog(float(u1[i])) + arb(self.Hpm_hi[a[i] - 1]))
            if bb[i] > a[i]:
                c = max(c, float(self.Kq_hi[a[i]:bb[i]].max()))
            K1[i] = max(0.0, c)
        b = self.b.astype(float)
        u2, v2 = b[:-1], b[1:]
        L2 = self.LBhi[1:]; E2 = self.EBhi[:, 1:].T
        K2 = np.array([max(0.0, hi(arb(self.LBhi[j + 1]) - 2 * loglog(int(self.b[j])) + arb(self.HBpm[j])))
                       for j in range(len(u2))])
        return dict(u=np.concatenate((u1, u2)), v=np.concatenate((v1, v2)), L=np.concatenate((L1, L2)),
                    E=np.concatenate((E1, E2)), logK=np.concatenate((K1, K2)))


def logX_lo_arb(u, Q):
    """log of a lower bound for X valid for all odd N > sqrt(8) u^16 (RIGOROSO.md sec. 1, 3)."""
    N = arb(8).sqrt() * arb(u) ** 16
    if Q == (3,):
        return ((N - 1) / 3).log()
    return ((N - 1) / 5 - 3).log()


def certify(u, v, L_hi, E_hi, logK_hi, s, jd, Q):
    """Rigorous lower bound (float) of log(main) - log(error) for every z in (u,v], or -inf."""
    k = 2 + arb(logK_hi)
    lf = logF(arb(s), k)
    if lf is None:
        return -math.inf
    d = arb(DELTAS_NUM[jd]) / 128
    m = logX_lo_arb(u, Q) - arb(L_hi) + lf - (2 * CQ[Q]).log() - arb(s) * d * arb(v).log() - arb(E_hi)
    return lo(m)


def optimise(u, v, L, E, logK, Q, ns=160, chunk=4000):
    """float64 search of (s, delta) per cell; returns s (float), jd (int), approx margin."""
    n = len(u)
    xf, xs = (1 / 3, 0.0) if Q == (3,) else (1 / 5, 3.0)
    logX = np.log((math.sqrt(8) * u ** 16 - 1) * xf - xs)
    lC = math.log(2 * float(CQ[Q].mid()))
    lv = np.log(v)
    best = np.full(n, -np.inf); bs = np.zeros(n); bj = np.zeros(n, int)
    frac = np.linspace(1e-7, 10, ns)
    for c0 in range(0, n, chunk):
        sl = slice(c0, min(n, c0 + chunk)); k = 2 + logK[sl]
        S = (2 * k + 3)[:, None] + frac[None, :]
        LF = logF_np(S, k[:, None])
        base = logX[sl] - L[sl] - lC
        for jd, d in enumerate(DELTAS):
            vals = LF - S * (d * lv[sl])[:, None]
            j = np.argmax(vals, axis=1); r = vals[np.arange(len(j)), j] + base - E[sl, jd]
            up = r > best[sl]; idx = np.arange(sl.start, sl.stop)[up]
            best[idx] = r[up]; bs[idx] = S[np.arange(len(j)), j][up]; bj[idx] = jd
    # golden refinement of s for the chosen delta
    k = 2 + logK; d = np.array(DELTAS)[bj]; h = 10 / (ns - 1)
    a_ = np.maximum(bs - h, 2 * k + 3 + 1e-7); b_ = bs + h
    g = lambda s: logF_np(s, k) - s * d * lv
    gr = (math.sqrt(5) - 1) / 2
    for _ in range(40):
        c1 = b_ - gr * (b_ - a_); c2 = a_ + gr * (b_ - a_)
        f1, f2 = g(c1), g(c2)
        m = f1 > f2
        b_ = np.where(m, c2, b_); a_ = np.where(m, a_, c1)
    s_ref = (a_ + b_) / 2
    better = g(s_ref) > g(bs)
    bs = np.where(better, s_ref, bs)
    approx = logX - L - lC + g(bs) - E[np.arange(n), bj]
    return bs, bj, approx


def logK_an_hi(Y, logKex_hi):
    """log K valid for every z >= Y (RIGOROSO.md sec. 4, Lemma A1):
       K_an = max(1, K_exact(Y)) * exp(2/log^2 Y + 1/(Y-2) + 1/(Y-1))   (Rosser-Schoenfeld Thm 5, both sides)."""
    Y = arb(Y)
    return hi(arb(max(0.0, logKex_hi)) + 2 / Y.log() ** 2 + 1 / (Y - 2) + 1 / (Y - 1))
