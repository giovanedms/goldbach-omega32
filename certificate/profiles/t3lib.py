"""
T3 / Campbell (arXiv 2608.09488) -- exploratory (NON-rigorous, double precision) model of the
Friedlander-Iwaniec Lambda^- Lambda^2 lower-bound sieve inequality used in Sec. 4-5 of the paper,
generalised to: target M (Omega(ab) <= M), parity of N, pre-sieve set Q, per-interval K,
continuous s, delta grid.

Notation follows the paper:
  Theorem 2.1: S(A,z) >= X V(z) F(s,k) - 2 R4(A,D),  k = kappa + log K, s = log D / log z, s >= 2k+3
  F(s,k) = 1 - (s+3)/(2 e^k) * (2 e k/(s-3))^((s-3)/2)                       (2.2)
  2 R4(A,z^s) <= 2 C_Q z^{s delta} prod_{p<z, p notin Q} (1 + 8/p^delta)      (4.4), C_{3} = 19/6
  X V(z) >= X_lo * (1/2) prod_{p<z, p odd, p notin Q} (1 - 2/p)
Everything here is float64; it is an exploration tool, not a proof.
"""
import math
import numpy as np

E = math.e


def primes_upto(n):
    """All primes < n (numpy, odd-only sieve)."""
    if n <= 2:
        return np.zeros(0, dtype=np.int64)
    m = (n - 1) // 2  # index i <-> 2i+1, i < m  (numbers < n)
    s = np.ones(m, dtype=bool)
    s[0] = False
    r = int(math.isqrt(n)) + 1
    for i in range(1, (r - 1) // 2 + 1):
        if s[i]:
            p = 2 * i + 1
            s[(p * p - 1) // 2::p] = False
    odd = 2 * np.nonzero(s)[0].astype(np.int64) + 1
    odd = odd[odd < n]
    return np.concatenate(([2], odd)).astype(np.int64)


def F(s, k):
    """Friedlander-Iwaniec sieve factor (2.2). Vectorised in s."""
    s = np.asarray(s, dtype=float)
    with np.errstate(all="ignore"):
        val = 1.0 - (s + 3.0) / (2.0 * math.exp(k)) * ((2.0 * E * k) / (s - 3.0)) ** ((s - 3.0) / 2.0)
    val = np.where(s >= 2 * k + 3, val, -np.inf)
    return val


def C_Q(Q):
    """Constant C with |r_d| <= C 2^{omega(d)} for d>1 (d | P(z)).
    Q={3}: paper, Lemma 4.3: 19/6.  General Q (our derivation, see SPEC sec. 2.4):
    r_d = sum_{e | prod Q} mu(e) (rho_{de} - g(d) rho_e);  |rho_l| <= 2^{omega(l)}, rho_1 = 0, g(d) <= 1/2
    => |r_d| <= 3^{|Q|} 2^{omega(d)} + (3^{|Q|}-1)/2 <= (3^{|Q|} + (3^{|Q|}-1)/4) 2^{omega(d)}."""
    Q = tuple(sorted(Q))
    if Q == (3,):
        return 19.0 / 6.0
    t = 3 ** len(Q)
    return t + (t - 1) / 4.0


class PrimeData:
    """Cumulative prime data up to LIM for a fixed pre-sieve set Q."""

    def __init__(self, lim, Q=(3,), deltas=None, sample_fn=None):
        """sample_fn(p_all) -> array of z-values at which the delta-sums are needed (memory saver);
        if None, full cumulative arrays are kept (needs ~8 bytes x #primes x #deltas)."""
        self.lim = lim
        self.Q = tuple(sorted(Q))
        p = primes_upto(lim)
        self.p_all = p
        keep = ~np.isin(p, np.array(self.Q, dtype=np.int64))
        self.p = p[keep]
        pf = self.p.astype(float)
        self.pf = pf  # float copy: searchsorted with float keys on an int array would cast the whole array per call
        # ell_p = log (1-g(p))^{-1} worst case: p=2 -> log 2; odd -> -log(1-2/p)
        ell = np.where(self.p == 2, math.log(2.0), -np.log1p(-2.0 / pf))
        self.ell = ell
        self.cumL = np.concatenate(([0.0], np.cumsum(ell)))  # cumL[i] = sum_{j<i} ell_j  (primes before index i)
        # log W after including primes up to index i-1:  W = (1/2) prod_{odd p notin Q}(1-2/p) = exp(-cumL)
        # (because ell_2 = log 2 and ell_p = -log(1-2/p))  => logW(v) = -L(v)
        self.deltas = np.array(deltas if deltas is not None else np.round(np.arange(0.20, 0.991, 0.01), 4))
        self.cumE = {}
        self.sidx = None
        if sample_fn is not None:
            zs = np.asarray(sample_fn(p), dtype=float)
            self.sidx = np.unique(np.searchsorted(self.p, zs, side="left"))
        for d in self.deltas:
            e = np.log1p(8.0 * pf ** (-d))
            full = np.concatenate(([0.0], np.cumsum(e)))
            self.cumE[float(d)] = full if self.sidx is None else full[self.sidx].copy()
            del full, e
        # dimension-condition prefix max:  H_i = max_{w prime, index<=i} (2 log log w - L(w)),  L(w)=cumL[index(w)]
        with np.errstate(divide="ignore"):
            llw = 2.0 * np.log(np.log(pf))
        Hw = llw - self.cumL[:-1]
        Hw[0] = 2.0 * math.log(math.log(2.0)) - 0.0  # w = 2
        # w=q for q in Q is dominated (same product as next prime, smaller log w) -- ignore
        self.Hpref = np.maximum.accumulate(Hw)
        # logK at z -> p_i^+ : L(p_i^+) - 2 log log p_i + max(0, Hpref_i)
        self.logK_at = np.maximum(0.0, self.cumL[1:] - llw + self.Hpref)

    def idx_lt(self, x):
        """number of (kept) primes < x."""
        return int(np.searchsorted(self.pf, x, side="left"))

    def logL(self, x):
        return self.cumL[self.idx_lt(x)]

    def logK_interval(self, u, v):
        """log of the smallest K valid for every z in (u, v] (worst-case g)."""
        i0 = self.idx_lt(u + 1e-12)  # primes <= u  -> z = u^+ uses primes <= u
        # candidate z->u^+ : primes p <= u included
        j_u = int(np.searchsorted(self.pf, u, side="right")) - 1
        vals = []
        if j_u >= 0:
            vals.append(self.logK_at[j_u] - 2 * math.log(math.log(u)) + 2 * math.log(math.log(self.p[j_u])))
        j_v = self.idx_lt(v)  # primes < v
        if j_v > j_u + 1:
            vals.append(float(self.logK_at[j_u + 1:j_v].max()))
        return max(vals) if vals else 0.0

    def logE(self, v, d):
        i = self.idx_lt(v)
        if self.sidx is None:
            return self.cumE[float(d)][i]
        k = int(np.searchsorted(self.sidx, i))
        assert k < len(self.sidx) and self.sidx[k] == i, "v not a sample point"
        return self.cumE[float(d)][k]


def best_log_ratio(pd, u, v, logX, logK=None, s_grid=None, kappa=2.0):
    """max over (s, delta) of log(main) - log(error) on the interval z in (u,v].
    logX = log of lower bound for X (at the smallest N of the interval).
    Returns (log_ratio, s, delta, K)."""
    if logK is None:
        logK = pd.logK_interval(u, v)
    k = kappa + logK
    if s_grid is None:
        s_grid = np.linspace(2 * k + 3, 2 * k + 3 + 14, 7001)
    Fs = F(s_grid, k)
    ok = Fs > 0
    s_ok = s_grid[ok]
    logF = np.log(Fs[ok])
    logW = -pd.logL(v)
    lv = math.log(v)
    C = C_Q(pd.Q)
    best = (-math.inf, None, None)
    for d in pd.deltas:
        le = math.log(2 * C) + pd.logE(v, d)
        vals = logF - s_ok * d * lv
        j = int(np.argmax(vals))
        lr = logX + logW + vals[j] - le
        if lr > best[0]:
            best = (lr, float(s_ok[j]), float(d))
    return best[0], best[1], best[2], math.exp(logK)
