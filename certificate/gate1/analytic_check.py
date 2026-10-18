"""Independent Case-3 (analytic range) check for odd N, Omega<=32, N = sqrt8 z^16 (and paper r=16 variant A)."""
import sys, math, numpy as np
sys.path.insert(0, sys.path[0]); from g1lib import *
from scipy.optimize import minimize_scalar
def worst_log_ratio(logz, logmain_fn, k, C, B1Q, phi):
    lz = logz
    best = -1e9
    for s in np.linspace(2*k+3+1e-9, 17, 1500):
        lf = logF(s, k)
        if not np.isfinite(lf): continue
        # min over alpha of 8 phi(a) - a s
        r = minimize_scalar(lambda a: 8*phi(a) - a*s, bounds=(1e-4, 3), method="bounded")
        lerr = math.log(2*C) + s*lz + 8*math.log(lz) + 8*B1Q + r.fun
        best = max(best, logmain_fn(lz) + lf - lerr)
    return best
# Variant A: paper constants, r=16 (z=N^{1/16}), s=14.66, alpha=.3246
def A(log10z):
    lz = log10z*math.log(10); logN = 16*lz
    main = logN - math.log(3) + math.log(1.241) - 2*math.log(lz) + math.log(0.0563)
    err = math.log(2*0.563) + 14.66*lz + 8*math.log(lz)
    return main - err
for e in [11.5, 11.9, 12.0, 12.1, 12.25]:
    print(f"A: log10 z={e}: log ratio={A(e):+.3f}  (N=10^{16*e:.1f})")
# root of A
lo, hi = 10, 13
for _ in range(60):
    mid = (lo+hi)/2; (lo, hi) = (lo, mid) if A(mid) > 0 else (mid, hi)
print(f"A closes at log10 z = {hi:.3f}, N = 10^{16*hi:.1f}")
# Variant Q35 (C-style): K split at 1e8, cV from exact product to 1e8, phi = e^a - 1
m = Model(10**8 + 1, [3, 5])
iz = m.n_lt(1e8)
logK8 = max(0.0, m.L[iz] - 2*math.log(math.log(1e8)) + m.H[iz-1])
K = math.exp(logK8 + 2/math.log(1e8)**2 + 2.1/(1e8-2)); k = 2 + math.log(K)
cV = math.exp(-m.L[iz]) * math.log(1e8)**2 * math.exp(-2/math.log(1e8)**2 - 2.1/(1e8-2))
B = 0.2614972128
print(f"Q35: K={K:.5f} cV={cV:.5f}")
for e in [8.0, 8.5, 9, 10, 12, 20, 40]:
    lz0 = e*math.log(10)
    B1Q = B + 1/(2*lz0**2) - 1/3 - 1/5
    mainf = lambda lz: 16*lz + 0.5*math.log(8) + math.log(1/5) + math.log(cV) - 2*math.log(lz)
    print(f"Q35 C-style: log10 z={e}: log ratio={worst_log_ratio(lz0, mainf, k, 11.0, B1Q, lambda a: math.exp(a)-1):+.3f}")
