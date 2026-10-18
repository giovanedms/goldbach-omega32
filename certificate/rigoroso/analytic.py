"""Rigorous analytic range z >= Y, odd N, Omega<=32, N = sqrt8 z^16 (RIGOROSO.md sec. 4).
For t = log z >= t_Y = log Y, with fixed s, alpha (eta = alpha/t, delta = 1 - eta):
  M(t) >= 16t + log sqrt8 - log qX + log(1 - c1/N) + log cV - 10 log t + log F(s,k) - log(2 C_Q)
          - s t + alpha s - 8 (B - sigma_Q) - 4/t^2 - 8 (e^alpha - 1)
  cV = exp(-L(Y)) log^2 Y exp(-2/log^2 Y - 1/(Y-2) - 1/(Y-1)),   k = 2 + log K_an(Y)
dM/dt >= (16 - s) - 10/t + 8/t^3 > 0 whenever (16 - s) t_Y > 10, so M(t) >= M(t_Y) for all t >= t_Y.
usage: analytic.py Q Y [K]   (K omitted: K_an(Y) from data_<Q>.json; K given: max(K, K_an(Y)) is used)
"""
import sys, json, math
import numpy as np
from flint import arb
import rlib
from rlib import hi, lo, logF, logF_np, CQ

Q = tuple(int(x) for x in sys.argv[1].split(",")); Y = int(float(sys.argv[2]))
Kuser = sys.argv[3] if len(sys.argv) > 3 else None
d = json.load(open(f"data_Q{''.join(map(str, Q))}.json"))["boundaries"][str(Y)]
L_hi, lKex = d["L_hi"], d["logKexact_hi"]
lKan = rlib.logK_an_hi(Y, lKex)
lK = lKan if Kuser is None else max(lKan, hi(arb(Kuser).log()))
Ya = arb(Y); tY = Ya.log()
B_up = arb("0.2614972128477")                    # Mertens constant 0.26149721284764278... (upper bound)
sig = sum(arb(1) / q for q in Q)
qX, c1 = (3, 1) if Q == (3,) else (5, 16)
lcV = -arb(L_hi) + 2 * tY.log() - 2 / tY ** 2 - 1 / (Ya - 2) - 1 / (Ya - 1)

def M_arb(s, al):
    k = 2 + arb(lK); s = arb(s); al = arb(al)
    lf = logF(s, k)
    if lf is None or not ((16 - s) * tY > 10) or not (al > 0) or not (al < tY):
        return None
    N = arb(8).sqrt() * (16 * tY).exp()
    return (N.log() - arb(qX).log() + (1 - c1 / N).log() + lcV - 10 * tY.log() + lf - (2 * CQ[Q]).log()
            - s * tY + al * s - 8 * (B_up - sig) - 4 / tY ** 2 - 8 * (al.exp() - 1))

# float search: alpha*(s) = log(s/8); s on a grid
t = float(tY.mid()); k = 2 + lK
ss = np.linspace(2 * k + 3 + 1e-6, 16 - 10 / t - 1e-6, 20001)
al = np.log(ss / 8)
val = (16 * t + 0.5 * math.log(8) - math.log(qX) - L_hi + 2 * math.log(t) - 10 * math.log(t) + logF_np(ss, k)
       - math.log(2 * float(CQ[Q].mid())) - ss * t + al * ss - 8 * (0.2614972128 - float(sig.mid())) - 8 * (np.exp(al) - 1))
j = int(np.argmax(val)); s_opt, a_opt = float(ss[j]), float(al[j])
m = M_arb(s_opt, a_opt)
print(f"Q={Q} Y={Y:.0e} K_exact(Y)<={math.exp(lKex):.6f} K_an(Y)<={math.exp(lKan):.6f} K_used={math.exp(lK):.6f} "
      f"cV>={math.exp(lo(lcV)):.6f}")
print(f"  s={s_opt:.5f} alpha={a_opt:.5f} (16-s)log Y={(16-s_opt)*t:.3f} (>10 needed)  float M(Y)={val[j]:+.5f}")
print(f"  CERTIFIED M(t) >= {lo(m):+.6f} for all z >= Y" if m is not None and lo(m) > 0 else
      f"  NOT certified: M lower = {lo(m) if m is not None else None}")
