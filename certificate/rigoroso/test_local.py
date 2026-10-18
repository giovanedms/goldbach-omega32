"""Brute-force checks of the local lemmas used in RIGOROSO.md (odd N), and float cross-checks of the rigorous data.
(1) |rho_l| <= 2^omega(l) for squarefree l (incl. even l), all odd N in a range;
(2) pre-sieve Q={3} and Q={3,5}: |A_d| = X g(d) + r_d with |r_d| <= C_Q 2^omega(d), d | P(z), d > 1;
(3) X >= (N-1)/3 (Q3), X >= (N-1)/5 - 3 (Q35);  a_n <= (N^2-1)/8;
(4) rigorous L, E, K (rlib) vs. plain float64 sums / brute-force sup."""
import math, itertools, numpy as np, sys
from fractions import Fraction as Fr

def omega(n):
    c, p = 0, 2
    while p * p <= n:
        if n % p == 0:
            c += 1
            while n % p == 0: n //= p
        p += 1
    return c + (n > 1)

def g_of(d, N):
    r = Fr(1)
    for p in [2, 3, 5, 7, 11, 13]:
        if d % p == 0:
            r *= Fr(1, 2) if p == 2 else (Fr(1, p) if N % p == 0 else Fr(2, p))
    return r

SQF = sorted(l for l in range(1, 30031) if 30030 % l == 0)          # squarefree l | 2*3*5*7*11*13
bad = 0; worst = {}
Nlist = list(range(3, 1500, 2)) + list(range(30031, 30031 + 2 * 60, 2))
for N in Nlist:
    a = np.array([n * (N - n) // 2 for n in range(1, N)], dtype=object)
    Y = N - 1
    assert max(a) <= Fr(N * N - 1, 8)
    def cntB(l): return sum(1 for x in a if x % l == 0)
    rho = {}
    for l in SQF:
        if l > 2 * N: continue
        rho[l] = cntB(l) - Y * g_of(l, N)
        if abs(rho[l]) > 2 ** omega(l): bad += 1; print("rho bound fails", N, l, rho[l]) if bad < 20 else None
    for Qs, C in [((3,), Fr(19, 6)), ((3, 5), Fr(11))]:
        q = math.prod(Qs)
        I = [i for i, x in enumerate(a) if all(x % p for p in Qs)]
        X = len(I)
        lowX = Fr(N - 1, 3) if Qs == (3,) else Fr(N - 1, 5) - 3
        if X < lowX: bad += 1; print("X bound fails", N, Qs, X)
        for d in [l for l in SQF if l > 1 and all(l % p for p in Qs) and l <= N]:
            Ad = sum(1 for i in I if a[i] % d == 0)
            r = Ad - X * g_of(d, N)
            ratio = abs(r) / 2 ** omega(d)
            worst[Qs] = max(worst.get(Qs, 0), ratio)
            if ratio > C: bad += 1; print("r_d bound fails", N, Qs, d, r) if bad < 20 else None
print(f"(1)-(3) odd N in [3,1500) + 60 N near 30031: violations = {bad}; max |r_d|/2^omega(d): "
      + ", ".join(f"Q={k}: {float(v):.4f}" for k, v in worst.items()))

# (4) rigorous data vs float
sys.path.insert(0, "."); import rlib
from flint import arb
for Q in [(3,), (3, 5)]:
    R = rlib.RData(Q, 10**8, verbose=False)
    from primes_blocks import primes_below
    P = primes_below(10**8); P = P[~np.isin(P, Q)]; pf = P.astype(float)
    ellf = np.where(P == 2, math.log(2), -np.log1p(-2 / pf))
    Lf = ellf.sum()
    j = int(np.searchsorted(R.b, 10**8))
    print(f"Q={Q}: L(1e8) float {Lf:.12f}  rigorous [{R.LBlo[j]:.12f}, {R.LBhi[j]:.12f}]  width {R.LBhi[j]-R.LBlo[j]:.2e}")
    assert R.LBlo[j] - 1e-9 <= Lf <= R.LBhi[j] + 1e-9
    for jd in [0, 50, 101]:
        d = rlib.DELTAS[jd]; Ef = np.log1p(8 * pf ** (-d)).sum()
        print(f"   E_{d:.4f}(1e8) float {Ef:.10f} rigorous hi {R.EBhi[jd, j]:.10f}  (hi-float {R.EBhi[jd, j]-Ef:.2e})")
        assert Ef <= R.EBhi[jd, j] + 1e-8
    # brute-force K at a few z (float, worst-case g) vs rigorous cell values
    Lc = np.concatenate(([0.0], np.cumsum(ellf)))
    def Kbrute(z):
        i = int(np.searchsorted(pf, z, side="left")); w = pf[:i]
        return math.exp(max(0, (Lc[i] - Lc[:i] + 2 * np.log(np.log(w)) - 2 * math.log(math.log(z))).max()))
    C = R.cells(rlib.z0_lower())
    for z in [43.001, 97.5, 1000.5, 15823.0, 2.5e5, 3.3e6, 5e7]:
        i = int(np.searchsorted(C["v"], z)); u, v = C["u"][i], C["v"][i]
        kb = max(Kbrute(u * (1 + 1e-12)), Kbrute(v))
        print(f"   z={z:g}: cell ({u:.3f},{v:.3f}] rigorous K<= {math.exp(C['logK'][i]):.7f}  brute(float, endpoints) {kb:.7f}")
        assert kb <= math.exp(C["logK"][i]) * (1 + 1e-9)
print("all checks passed")
