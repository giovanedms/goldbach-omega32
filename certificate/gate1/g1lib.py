"""Gate 1 -- independent re-implementation (does NOT import t3lib). float64, exploratory."""
import math, numpy as np

def sieve(n):
    b = np.ones(n, dtype=np.bool_); b[:2] = False
    for i in range(2, int(n**0.5) + 1):
        if b[i]: b[i*i::i] = False
    return np.flatnonzero(b).astype(np.int64)

def logF(s, k):
    s = np.asarray(s, float)
    with np.errstate(all="ignore"):
        F = 1 - (s + 3) / (2 * np.exp(k)) * (2 * math.e * k / (s - 3)) ** ((s - 3) / 2)
        return np.where((s >= 2 * k + 3) & (F > 0), np.log(F), -np.inf)

class Model:
    """Q = excluded pre-sieve primes. Arrays over kept primes P (p not in Q)."""
    def __init__(self, lim, Q):
        pa = sieve(lim)
        self.p = pa[~np.isin(pa, Q)]
        pf = self.p.astype(float); self.pf = pf
        ell = np.where(self.p == 2, math.log(2), -np.log1p(-2 / pf))   # log (1-g)^-1 worst case
        self.L = np.concatenate(([0.0], np.cumsum(ell)))               # L[j]= sum over idx<j
        llp = np.log(np.log(pf))
        H = np.maximum.accumulate(2 * llp - self.L[:-1])               # max_w<=p_i (2loglog w - L(w))
        self.H = H
        self.logK_q = np.maximum(0, self.L[1:] - 2 * llp + H)         # z -> p_i^+
    def n_le(self, x): return int(np.searchsorted(self.pf, x, side="right"))
    def n_lt(self, x): return int(np.searchsorted(self.pf, x, side="left"))
    def logK_uplus(self, u):
        i = self.n_le(u)
        return max(0.0, self.L[i] - 2 * math.log(math.log(u)) + self.H[i - 1])
    def logK_brute(self, z):
        """direct sup over w (prime w and w=2) for z (not prime): check."""
        P = self.pf[self.pf < z]; best = 0.0
        for j in range(len(P)):
            val = (self.L[len(P)] - self.L[j]) + 2*math.log(math.log(P[j])) - 2*math.log(math.log(z))
            best = max(best, val)
        return best
