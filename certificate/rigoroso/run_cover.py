"""Rigorous cover of z in (u0, XMAX] (u0 <= z0*), odd N, Omega<=32, N = sqrt8 z^16.
usage: run_cover.py Q(3|3,5) XMAX Kmode
  Kmode = interval            : K computed (rigorous upper bound) separately on every cell
          uniform:<K>         : one K for every cell (e.g. uniform:1.146 = paper Lemma 4.1(1))
          uniform:sup         : one K = rigorous sup of the per-cell K over the whole finite range
Every cell (u,v] gets a certified lower bound (arb) of log(main/error) valid for ALL z in (u,v].
Writes out_cover_<Q>_<Kmode>.txt and data_<Q>.json (boundary data for analytic.py)."""
import sys, math, json, time
import numpy as np
from flint import arb
import rlib
from rlib import RData, certify, optimise, z0_lower, hi, DELTAS

Q = tuple(int(x) for x in sys.argv[1].split(",")); XMAX = int(float(sys.argv[2])); mode = sys.argv[3]
tag = f"Q{''.join(map(str, Q))}_{mode.replace(':', '')}"
out = open(f"out_cover_{tag}.txt", "w")
def say(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.write(s + "\n"); out.flush()
t0 = time.time()
R = RData(Q, XMAX)
u0 = z0_lower()
C = R.cells(u0)
u, v, L, E, logK = C["u"], C["v"], C["L"], C["E"], C["logK"]
say(f"# Q={Q} XMAX={XMAX:.3g} mode={mode} u0={u0!r} (<= z0*, checked exactly) cells={len(u)} "
    f"first u={u[0]!r} last v={v[-1]!r}  build {time.time()-t0:.0f}s")
assert np.all(u[1:] == v[:-1]) and u[0] == u0 and v[-1] == XMAX      # contiguous cover
logK_cell_max = float(logK.max()); jm = int(np.argmax(logK))
say(f"# per-cell K (rigorous upper bounds): max K = {math.exp(logK_cell_max):.6f} on ({u[jm]:.4f},{v[jm]:.4f}]")
if mode == "envelope":
    # K(z) replaced by sup_{z'' >= z} K(z''): per-cell suffix max, and the analytic tail bound K_an(XMAX) for z >= XMAX
    Lh, lK = R.at_boundary(XMAX); tail = rlib.logK_an_hi(XMAX, lK)
    logK = np.maximum(np.maximum.accumulate(logK[::-1])[::-1], tail)
    say(f"# envelope K: tail K_an(XMAX) = {math.exp(tail):.6f}")
if mode.startswith("uniform"):
    val = mode.split(":")[1]
    lk = logK_cell_max if val == "sup" else hi((arb(val)).log())
    logK = np.full(len(u), lk)
    say(f"# uniform K = exp({lk!r}) = {math.exp(lk):.6f}")
s, jd, approx = optimise(u, v, L, E, logK, Q)
t1 = time.time()
cert = np.array([certify(float(u[i]), float(v[i]), float(L[i]), float(E[i, jd[i]]), float(logK[i]),
                         float(s[i]), int(jd[i]), Q) for i in range(len(u))])
say(f"# certification {time.time()-t1:.0f}s; max |float approx - certified| = {np.nanmax(np.abs(approx - cert)[np.isfinite(cert)]):.2e}")
nf = int((cert <= 0).sum())
i = int(np.argmin(cert))
say(f"RESULT cells={len(u)} failing(cert<=0)={nf}")
say(f"RESULT min certified margin = {cert[i]:.6f} on ({u[i]:.4f},{v[i]:.4f}] s={s[i]:.5f} delta={DELTAS[jd[i]]:.6f} K={math.exp(logK[i]):.6f}")
for lo_, hi_ in [(u0, 1e2), (1e2, 1e3), (1e3, 1e4), (1e4, 1e5), (1e5, 1e6), (1e6, 1e7), (1e7, 1e8), (1e8, 1e9), (1e9, 1e10)]:
    msk = (v > lo_) & (u < hi_)
    if msk.any():
        j = np.flatnonzero(msk)[np.argmin(cert[msk])]
        say(f"  z in ({lo_:.3g},{hi_:.3g}]: min cert {cert[j]:+.5f} on ({u[j]:.3f},{v[j]:.3f}] s={s[j]:.4f} "
            f"d={DELTAS[jd[j]]:.4f} K={math.exp(logK[j]):.6f}  fails={int((cert[msk]<=0).sum())}")
if nf:
    bad = np.flatnonzero(cert <= 0)
    say(f"# failing z-range: [{u[bad[0]]:.4f}, {v[bad[-1]]:.4f}]")
np.savez(f"cells_{tag}.npz", u=u, v=v, cert=cert, s=s, jd=jd, logK=logK)
if mode == "interval":
    data = {"Q": Q, "logK_cell_max": logK_cell_max, "boundaries": {}}
    for Y in [10**8, 10**9, 10**10]:
        if Y <= XMAX:
            Lh, lK = R.at_boundary(Y)
            data["boundaries"][str(Y)] = {"L_hi": Lh, "logKexact_hi": lK}
            say(f"# boundary Y={Y:.0e}: L_hi={Lh!r} (W(Y)>={math.exp(-Lh):.8g})  K_exact(Y)<={math.exp(lK):.6f}")
    json.dump(data, open(f"data_Q{''.join(map(str, Q))}.json", "w"), indent=1)
say(f"# total {time.time()-t0:.0f}s")
