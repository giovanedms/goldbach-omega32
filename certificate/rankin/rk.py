import math
def primes(n):
    s=bytearray([1])*(n+1); s[0]=s[1]=0
    for i in range(2,int(n**.5)+1):
        if s[i]: s[i*i::i]=bytearray(len(s[i*i::i]))
    return [i for i in range(n+1) if s[i]]
for z in [43.5, 1000.0, 16000.0, 1e5]:
    P=[p for p in primes(int(z)) if p<z and p!=3]
    g={p:(0.5 if p==2 else 2/p) for p in P}
    # local K at level z: sup over w (prime endpoints and w=2)
    K=0; prod=1.0
    for p in reversed(P):
        prod*=1/(1-g[p]); K=max(K,prod*(math.log(p)/math.log(z))**2)
    k=2+math.log(K)
    worst=-1e9
    for c in [0.1*i for i in range(1,40)]:
        eps=c/math.log(z)
        lhs=sum(g[p]*(p**eps-1) for p in P)
        worst=max(worst,lhs-k*(math.exp(c)-1))
    print(f"z={z:g} K_local={K:.5f} max_c[LHS - k(e^c-1)]={worst:.4f} (deve ser <=0)")
