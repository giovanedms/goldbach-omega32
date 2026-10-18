"""Cross-check: forca bruta independente dos lemas locais do pre-crivo Q={3,5} (N impar).
Confere: (i) identidade r_d = sum_{e|15} mu(e)(rho_{de} - g(d) rho_e) exatamente;
(ii) |rho_l| <= 2^omega(l); (iii) |r_d| <= 11*2^omega(d); (iv) X >= (N-1)/5 - 3;
(v) max de |rho_3|+|rho_5|+|rho_15| (periodico em N mod 60) -> constante C_Q mais fina."""
from fractions import Fraction as Fr
from itertools import combinations
import math

PR = [2, 7, 11, 13, 17]


def g(d, N):
    r = Fr(1)
    for p in [2, 3, 5, 7, 11, 13, 17, 19]:
        if d % p == 0:
            r *= Fr(1, 2) if p == 2 else (Fr(1, p) if N % p == 0 else Fr(2, p))
    return r


def om(d):
    return sum(1 for p in [2, 3, 5, 7, 11, 13, 17, 19] if d % p == 0)


ds = [1]
for k in range(1, len(PR) + 1):
    for c in combinations(PR, k):
        ds.append(math.prod(c))
viol = 0; worst = Fr(0); worstrho = Fr(0); maxsum = Fr(0); minXgap = None
Ns = list(range(3, 1201, 2)) + list(range(4001, 4401, 2)) + list(range(30031, 30031 + 120, 2))
for N in Ns:
    a = [n * (N - n) // 2 for n in range(1, N)]
    Y = N - 1
    rho = {}
    def B(l): return sum(1 for x in a if x % l == 0)
    for e in [1, 3, 5, 15]:
        for d in ds:
            l = d * e
            if l not in rho:
                rho[l] = B(l) - Y * g(l, N)
                worstrho = max(worstrho, abs(rho[l]) / 2 ** om(l))
                if abs(rho[l]) > 2 ** om(l): viol += 1
    I = [x for x in a if x % 3 and x % 5]
    X = len(I)
    if X < Fr(N - 1, 5) - 3: viol += 1
    gap = X - (Fr(N - 1, 5) - 3); minXgap = gap if minXgap is None else min(minXgap, gap)
    maxsum = max(maxsum, abs(rho[3]) + abs(rho[5]) + abs(rho[15]))
    mu = {1: 1, 3: -1, 5: -1, 15: 1}
    for d in ds[1:]:
        Ad = sum(1 for x in I if x % d == 0)
        rd = Ad - X * g(d, N)
        formula = sum(mu[e] * (rho[d * e] - g(d, N) * rho[e]) for e in mu)
        assert rd == formula, (N, d)
        worst = max(worst, abs(rd) / 2 ** om(d))
        if abs(rd) > 11 * 2 ** om(d): viol += 1
print(f"N testados={len(Ns)}  violacoes={viol}  identidade r_d exata OK")
print(f"max|rho_l|/2^omega = {float(worstrho):.4f}; max|r_d|/2^omega(d) = {float(worst):.4f}")
print(f"max(|rho3|+|rho5|+|rho15|) = {maxsum} = {float(maxsum):.4f}  -> C_Q fino = 9 + max/4 = {float(9 + maxsum / 4):.4f}")
print(f"min X - ((N-1)/5-3) = {minXgap}")
