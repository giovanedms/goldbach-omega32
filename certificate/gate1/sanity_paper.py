import sys, math, numpy as np
sys.path.insert(0, sys.path[0]); from g1lib import *
m = Model(10**6, [3])
# (1) paper Lemma 4.1(1) case 3: sup over z0<=z<286 of K  ~ 1.1458 at z->43+, w=5
z0 = 117465180365547648498934439**(1/16.5)
vals = [m.logK_uplus(z0)] + [m.logK_q[i] for i in range(len(m.p)) if z0 <= m.p[i] < 286]
print("K sup z0..286 =", math.exp(max(vals)))
for z in [40.5, 43.5, 100.3, 1000.5, 16150.5]:
    print("z", z, "K prefix", math.exp(m.logK_uplus(z)), "brute", math.exp(m.logK_brute(z)))
# (2) paper Case 2 S_lower on (97,151] delta .6, s 14.8, K 1.146 -> 1.0217e29
def S_paper(u, v, s, d):
    k = 2 + math.log(1.146)
    W = math.exp(-m.L[m.n_lt(v)])
    E = np.sum(np.log1p(8 * m.pf[m.pf < v] ** (-d)))
    return (u**16.5 - 1)/3 * W * math.exp(logF(s, k)) - 19/3 * v**(s*d) * math.exp(E)
print("S(97,151;14.8,.6) =", S_paper(97, 151, 14.8, 0.6), "(author 1.021694840443e+29)")
print("S(z0,97;14.8,.2) =", S_paper(z0, 97, 14.8, 0.2), "(author 1.876202819369e+22)")
