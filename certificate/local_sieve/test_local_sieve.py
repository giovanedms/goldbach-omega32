"""Checks of the Lambda^- Lambda^2 lower bound under the LOCAL dimension condition.

(A1)-(A3) are exact rational checks (with the weights psi(d) replaced by rational approximations);
(A4), (A4a) are floating-point checks with tolerance 1e-9 relative to V(z)^{-1} (resp. |mid|); (A5), (A6) use
absolute tolerances 1e-9 / 1e-12; (A6) checks
Proposition 4.4 and Theorem 4.1 in the range s >= 2k+3 (kappa = 2, smallest K allowed at z).
Checks every link of Appendix A for the actual sequence
a_n = n(N-n)/2 (N odd) or n(N-n) (N even), pre-sieved by Q, and for random densities:

  (A1) S(A,z) >= sum_d lambda_d |A_d|                     (combinatorial inequality, exact)
  (A2) sum_d lambda_d g(d) == V(D,z) from identity (2.3)  (algebraic identity, exact)
  (A3) |theta_e| <= 1  and  |lambda_d| <= 2 tau_4(d)         (weights)
  (A4) J^2 V(D,z) >= V(z)^{-1} - (log(Y_0 z^2)/log z) I(Y_0/z,z) (main-term bound; floats with margin)
  (A5) S(A,z) >= X V(z) (1 - c I V) - 2 R_4               (the assembled inequality)

Usage: python test_local_sieve.py [n_random]
"""
from __future__ import annotations

import itertools
import math
import random
import sys
from fractions import Fraction as Fr


def primes_upto(n: int) -> list[int]:
    s = bytearray([1]) * (n + 1)
    s[0:2] = b"\x00\x00"
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]


def squarefree_divisors(ps: list[int]):
    """All (d, tuple_of_primes) for squarefree d composed of ps."""
    out = []
    for r in range(len(ps) + 1):
        for c in itertools.combinations(ps, r):
            out.append((math.prod(c), c))
    return out


def weights(ps, g, Y, z):
    """Selberg-type y_d = J^{-1} w(d), w = min(1, log(Y/d)/log z); returns rho, y, J (floats->Fractions)."""
    h = {p: g[p] / (1 - g[p]) for p in ps}
    divs = squarefree_divisors(ps)
    hd = {d: math.prod((h[p] for p in c), start=Fr(1)) for d, c in divs}

    def w(d):  # exact enough: rational approximation of a real weight, consistency is what matters
        if d >= Y:
            return Fr(0)
        return Fr(min(1.0, math.log(Y / d) / math.log(z))).limit_denominator(10 ** 12)

    wd = {d: w(d) for d, _ in divs}
    J = sum(hd[d] * wd[d] for d, _ in divs)
    y = {d: wd[d] / J for d, _ in divs}
    gd = {d: math.prod((g[p] for p in c), start=Fr(1)) for d, c in divs}
    # Moebius inversion: g(e) rho_e = mu(e) sum_{e | d} h(d) y_d
    rho = {}
    for e, ce in divs:
        se = set(ce)
        tot = sum(hd[d] * y[d] for d, c in divs if se.issubset(c))
        rho[e] = (-1) ** len(ce) * tot / gd[e]
    return divs, h, hd, gd, y, rho, J


def lambdas(divs, rho, ps):
    """lambda_d = lambda+_d - sum_{p|d} (lambda+_d + lambda+_{d/p}),  lambda+_d = sum_{[d1,d2]=d} rho rho."""
    fac = {d: set(c) for d, c in divs}
    lp = {d: Fr(0) for d, _ in divs}
    nz = [d for d, _ in divs if rho[d] != 0]
    for d1 in nz:
        for d2 in nz:
            m = d1 * d2 // math.gcd(d1, d2)
            lp[m] += rho[d1] * rho[d2]
    lam = {}
    for d, c in divs:
        v = lp[d] - sum(lp[d] + lp[d // p] for p in c)
        lam[d] = v
    return lam


def check_identity_and_weights(ps, g, Y, z, tag):
    divs, h, hd, gd, y, rho, J = weights(ps, g, Y, z)
    assert rho[1] == 1, f"{tag}: rho_1 != 1"
    lam = lambdas(divs, rho, ps)
    # (A2) sum lambda_d g(d) vs identity (2.3)
    lhs = sum(lam[d] * gd[d] for d, _ in divs)
    rhs = sum(hd[d] * y[d] ** 2 for d, _ in divs)
    for d, c in divs:
        for p in ps:
            if p not in c:
                rhs -= g[p] * hd[d] * (y[d] - y[d * p]) ** 2
    assert lhs == rhs, f"{tag}: identity (2.3) fails: {float(lhs)} vs {float(rhs)}"
    # (A3)
    for d, c in divs:
        assert abs(rho[d]) <= 1, f"{tag}: |rho_{d}| = {float(rho[d])} > 1"
        assert abs(lam[d]) <= 2 * 4 ** len(c), f"{tag}: |lambda_{d}| too big"
    # (A4a) the sharper intermediate step: J^2 V(D,z) >= sum_{d<=Y/z} h - sum_{Y/z<d<Yz} h(d) log d/log z
    mid = sum(float(hd[d]) for d, _ in divs if d <= Y / z) - sum(
        float(hd[d]) * math.log(d) / math.log(z) for d, _ in divs if Y / z < d < Y * z)
    assert float(J * J * lhs) >= mid - 1e-9 * max(1.0, abs(mid)), f"{tag}: (A4a) fails"
    # (A4) J^2 V(D,z) >= V^{-1} - c I(Y/z)
    Vinv = math.prod((1 + h[p] for p in ps), start=Fr(1))
    I = sum(hd[d] for d, _ in divs if d >= Y / z)
    c = math.log(Y * z * z) / math.log(z)
    a4 = float(J * J * lhs) - (float(Vinv) - c * float(I))
    assert a4 >= -1e-9 * float(Vinv), f"{tag}: (A4) fails by {a4}"
    # V(D,z) >= V(z)(1 - c I V) needs J <= 1/V(z) AND a non-negative bracket (Appendix A, two cases)
    V = 1 / Vinv
    assert J <= Vinv, f"{tag}: J > 1/V(z)"
    bracket = 1 - c * float(I) * float(V)
    if bracket >= 0:
        assert float(lhs) >= float(V) * bracket - 1e-12, f"{tag}: V(D,z) >= V(1-cIV) fails"
    return divs, gd, lam, J, lhs, float(V), float(I), c


def local_K(ps, g, z):
    """Smallest K >= 1 with prod_{w<=p<z}(1-g(p))^{-1} <= K (log z/log w)^2 for all 2 <= w < z (kappa = 2).
    The sup over w in [p_i, p_{i+1}) is attained at w = p_i (product constant, (log w)^2 increasing)."""
    best = 1.0
    cands = sorted(set([2] + list(ps)))
    for w in cands:
        if w >= z:
            continue
        prod = 1.0
        for p in ps:
            if w <= p < z:
                prod /= (1 - float(g[p]))
        best = max(best, prod * (math.log(w) / math.log(z)) ** 2)
    return best


def F_sk(s, k):
    return 1 - (s + 3) / (2 * math.exp(k)) * (2 * math.e * k / (s - 3)) ** ((s - 3) / 2)


def goldbach_case(N: int, Q: tuple[int, ...], z: float, s: float):
    """Exact test on a_n for one N with sifting primes P = primes < z not in Q."""
    odd = N % 2 == 1
    a = {n: (n * (N - n) // 2 if odd else n * (N - n)) for n in range(1, N)}
    I_Q = [n for n in range(1, N) if all(a[n] % q for q in Q)]
    X = len(I_Q)
    ps = [p for p in primes_upto(int(z) + 1) if p < z and p not in Q]
    D = z ** s
    Y = math.sqrt(D / z)
    # local density g (Campbell 3.): g(2)=1/2, g(p)=2/p (p not | N), 1/p (p | N)
    g = {}
    for p in ps:
        if p == 2:
            g[p] = Fr(1, 2)
        else:
            g[p] = Fr(1, p) if N % p == 0 else Fr(2, p)
    tag = f"N={N} Q={Q} z={z} s={s}"
    divs, gd, lam, J, VDz, V, I, c = check_identity_and_weights(ps, g, Y, z, tag)
    # exact |A_d| and r_d
    Ad = {d: sum(1 for n in I_Q if a[n] % d == 0) for d, _ in divs}
    r = {d: Ad[d] - X * gd[d] for d, _ in divs}
    assert r[1] == 0
    P = math.prod(ps)
    S = sum(1 for n in I_Q if math.gcd(a[n], P) == 1)
    # (A1) S >= sum lambda_d |A_d|
    lower = sum(lam[d] * Ad[d] for d, _ in divs)
    assert S >= lower, f"{tag}: (A1) fails S={S} < {float(lower)}"
    # support: lambda_d != 0 only for d < D
    for d, _ in divs:
        if lam[d] != 0:
            assert d < D, f"{tag}: lambda supported at d={d} >= D"
    # (A5) assembled inequality with the crude bracket (1 - c I V)
    R4 = sum(4 ** len(cc) * abs(r[d]) for d, cc in divs if 1 < d < D)
    bracket = 1 - c * I * V
    a5 = X * V * bracket - 2 * float(R4)
    if bracket < 0:                       # trivial case: S >= 0 > a5
        assert a5 < 0
    assert S >= a5 - 1e-9, f"{tag}: (A5) fails S={S} < {a5}"
    # (A6) in the range of Theorem 4.1 (kappa = 2, local K at this z, s >= 2k+3):
    #      Proposition 4.4 bound c*I*V <= 1 - F(s,k), and the theorem S >= X V F - 2 R_4
    K = local_K(ps, g, z); k = 2 + math.log(K)
    checkedF = False
    if s >= 2 * k + 3:
        F = F_sk(s, k)
        assert c * I * V <= (1 - F) + 1e-12, f"{tag}: Prop 4.4 bound fails"
        assert S >= X * V * F - 2 * float(R4) - 1e-9, f"{tag}: Theorem 4.1 fails"
        checkedF = 2 if F > 0 else 1
    return S, float(lower), a5, checkedF


def main() -> int:
    nrand = int(sys.argv[1]) if len(sys.argv) > 1 else 300
    random.seed(20261018)
    # random densities on random prime sets
    for i in range(nrand):
        ps = sorted(random.sample(primes_upto(60), random.randint(1, 6)))
        g = {p: Fr(random.randint(1, 97), 100) for p in ps}
        z = max(ps) + random.random() * 5 + 0.01
        s = random.uniform(1.0, 6.0)
        Y = math.sqrt(z ** s / z)
        check_identity_and_weights(ps, g, Y, z, f"rand#{i}")
    print(f"random densities: {nrand} cases, (A2)(A3)(A4) OK")
    # actual Goldbach sequences
    n = nF = npos = 0
    for N in list(range(41, 400, 2)) + list(range(42, 400, 2)) + [1001, 2023, 4096, 9999]:
        for Q in ((3,), (3, 5)):
            for z in (6.5, 11.5, 14.0, 18.0):
                for s in (2.0, 3.0, 4.5, 7.5, 9.0, 11.0, 15.0, 17.0, 20.0, 24.0):
                    _, _, _, cf = goldbach_case(N, Q, z, s)
                    n += 1; nF += cf > 0; npos += cf == 2
    print(f"Goldbach sequences: {n} cases, (A1)-(A5) OK; "
          f"(A6) Proposition 4.4 and Theorem 4.1 checked in {nF} cases with s >= 2k+3, {npos} of them with F > 0")
    return 0


if __name__ == "__main__":
    sys.exit(main())
