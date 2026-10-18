"""Copy of the profile_general.py logic used for Figure F4. Saves every grid point.
Run (ulimit -v 4000000): python F4_gen_profiles.py ; output profiles_dump.csv = F4_profiles_float.csv."""
import math, sys, numpy as np
from F4_t3lib import *
N0 = 117465180365547648498934439
LIM=10**8
out=[]
for label,Q,c,r in [("M32_Q3",(3,),2**1.5,16.0),("M32_Q35",(3,5),2**1.5,16.0),("M31_Q3",(3,),2**1.5,15.5),("M31_Q35",(3,5),2**1.5,15.5)]:
    pd = PrimeData(LIM, Q=Q, deltas=np.round(np.arange(0.20,0.991,0.02),4))
    qf=1.0
    for q in Q: qf*=(1-2/q)
    z0=(N0/c)**(1/r)
    def logX(u):
        N=max(N0+1,c*u**r); return math.log((N-1)*qf-(3**len(Q)-1))
    zs=np.unique(np.concatenate(([z0*1.0000001],np.geomspace(max(z0,41),LIM*0.99,240))))
    for u in zs:
        lr,s,d,K=best_log_ratio(pd,u,u*1.000001,logX(u)); out.append((label,u,lr))
    del pd; import gc; gc.collect()
    print(label,min(x[2] for x in out if x[0]==label),flush=True)
with open("profiles_dump.csv","w") as f:
    f.write("case,z,margin\n")
    for l,u,lr in out: f.write(f"{l},{float(u)!r},{float(lr)!r}\n")
