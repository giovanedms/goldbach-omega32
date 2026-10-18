"""Independent full cover of z in (z0*, 1e8] for N odd, Omega<=32 (N = sqrt8 z^16), by elementary
geometric cells (z_g, z_{g+1}], ratio 1+eps. Each cell must have margin>0 on its own.
usage: cover_cells.py Q(3 | 3,5) eps"""
import sys, math, time, numpy as np
sys.path.insert(0, sys.path[0]); from g1lib import *
Q = [int(x) for x in sys.argv[1].split(",")]; eps = float(sys.argv[2])
LIM = 10**8; N0 = 117465180365547648498934439
t = time.time(); m = Model(LIM + 1, Q); print("# model", round(time.time()-t,1), "s, primes kept", len(m.p), flush=True)
z0 = (N0**2 / 8) ** (1/32)
n = int(math.ceil(math.log(LIM / z0) / math.log1p(eps)))
z = z0 * (1 + eps) ** np.arange(n + 1); z[-1] = LIM          # last endpoint exactly 1e8
u, v = z[:-1], z[1:]
if Q == [3]:
    C, xf, xsub = 19/6, 1/3, 0.0
else:
    C, xf, xsub = 11.0, 1/5, 2/3 + 2 + 4                       # |rho3|+|rho5|+|rho15|
N = np.maximum(math.sqrt(8) * u**16, float(N0))
logX = np.log((N - 1) * xf - xsub)
a = np.searchsorted(m.pf, u, side="right")                     # primes <= u
b = np.searchsorted(m.pf, v, side="left")                      # primes < v
logKu = np.maximum(0, m.L[a] - 2*np.log(np.log(u)) + m.H[a-1])
logK = logKu.copy()
ne = b > a
cm = np.array([m.logK_q[x:y].max() for x, y in zip(a[ne], b[ne])]) if ne.sum() < 400000 else None
logK[ne] = np.maximum(logK[ne], cm)
if len(sys.argv) > 3 and sys.argv[3] == 'fixed':
    logK[:] = math.log(float(sys.argv[4])); print('# K fixed', sys.argv[4])
k = 2 + logK
logW = -m.L[b]
lv = np.log(v)
deltas = np.round(np.arange(0.20, 0.9951, 0.005), 4)
best = np.full(n, -np.inf); bs = np.zeros(n); bd = np.zeros(n)
cum = {}
for d in deltas:
    e = np.concatenate(([0.0], np.cumsum(np.log1p(8 * m.pf ** (-d)))))
    cum[d] = e[b]
NS = 600; frac = np.linspace(0, 12, NS)
CH = 8000
for c0 in range(0, n, CH):
    sl = slice(c0, min(n, c0 + CH))
    S = (2*k[sl] + 3)[:, None] + frac[None, :]
    LF = logF(S, k[sl][:, None])
    base = logX[sl] + logW[sl] - math.log(2*C)
    for d in deltas:
        vals = LF - S * (d * lv[sl])[:, None]
        j = np.argmax(vals, axis=1); r = vals[np.arange(len(j)), j] + base - cum[d][sl]
        upd = r > best[sl]
        idx = np.arange(sl.start, sl.stop)[upd]
        best[idx] = r[upd]; bs[idx] = S[np.arange(len(j)), j][upd]; bd[idx] = d
i = int(np.argmin(best))
print(f"# Q={Q} eps={eps} z0*={z0:.6f} cells={n} last v={v[-1]:.1f}")
print(f"# min cell margin = {best[i]:.4f} on ({u[i]:.3f},{v[i]:.3f}] s={bs[i]:.3f} d={bd[i]} K={math.exp(logK[i]):.5f}")
print(f"# cells with margin<=0: {(best<=0).sum()}")
for lo, hi in [(z0,1e2),(1e2,1e3),(1e3,1e4),(1e4,1e5),(1e5,1e6),(1e6,1e7),(1e7,1e8+1)]:
    msk = (u >= lo) & (u < hi); j = np.flatnonzero(msk)[np.argmin(best[msk])]
    print(f"  z in [{lo:.3g},{hi:.3g}): min margin {best[j]:.4f} at {u[j]:.1f} s={bs[j]:.3f} d={bd[j]} K={math.exp(logK[j]):.5f}")
