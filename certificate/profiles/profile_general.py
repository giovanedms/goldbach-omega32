"""Pointwise margin profile (tiny intervals) for N = c z^r, pre-sieve set Q.
usage: profile_general.py Q LIM c r label"""
import math, sys, numpy as np
from t3lib import *
Q = tuple(int(x) for x in sys.argv[1].split(",")); LIM=int(float(sys.argv[2]))
c=float(eval(sys.argv[3])); r=float(sys.argv[4]); label=sys.argv[5]
N0 = 117465180365547648498934439
pd = PrimeData(LIM, Q=Q, deltas=np.round(np.arange(0.20,0.991,0.02),4))
qf=1.0
for q in Q: qf*=(1-2/q)
z0=(N0/c)**(1/r)
def logX(u):
    N=max(N0+1,c*u**r); return math.log((N-1)*qf-(3**len(Q)-1))
zs=np.unique(np.concatenate(([z0*1.0000001],np.geomspace(max(z0,41),LIM*0.99,60))))
res=[]
for u in zs:
    lr,s,d,K=best_log_ratio(pd,u,u*1.000001,logX(u)); res.append((lr,u,s,d,K))
mn=min(res)
print(f"{label}: Q={Q} c={c:.3f} r={r} z0={z0:.3f}  min logratio={mn[0]:.3f} at z={mn[1]:.4g} (s={mn[2]:.2f},d={mn[3]:.2f},K={mn[4]:.4f}); at z0: {res[0][0]:.2f}; at 1e8: {res[-1][0]:.2f}")
neg=[x for x in res if x[0]<=0]
if neg: print(f"   negative on z in [{min(x[1] for x in neg):.4g}, {max(x[1] for x in neg):.4g}] (grid)")
