"""Cross-check: reverificacao independente de celulas (nao importa rlib/primes_blocks).
Para cada celula (u,v] escolhida, recalcula do zero com crivo proprio + somas por primo (sem Jensen/corda):
  margem = log X_lo(u) - L(v) + log F(s, 2+log K_cel) - log(2 C_Q) - s*delta*log v - E_delta(v)
K_cel (cenario I) = sup_{z in (u,v]} sup_{w<z} prod_{w<=p<z,p in P}(1-2/p)^-1 (log w/log z)^2, via candidatos
z -> u+ e z -> q+ (q primo de P em (u,v)), w primo de P com w <= candidato.
Somas em float64 com fsum por bloco (erro << 1e-9); resultado comparado com cert do autor com folga 1e-7.
usage: python3 crosscheck_cells.py <seed>"""
import sys, math
import numpy as np
import mpmath as mp

mp.mp.dps = 40
import os; BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "rigoroso") + os.sep  # (copia: aponta para ../rigoroso)


def sieve(n):
    s = np.ones(n, dtype=bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.flatnonzero(s)


def logF(s, k):
    s = mp.mpf(s); k = mp.mpf(k)
    assert s >= 2 * k + 3
    F = 1 - (s + 3) / (2 * mp.e ** k) * (2 * mp.e * k / (s - 3)) ** ((s - 3) / 2)
    assert F > 0
    return mp.log(F)


def check(tag, Q, C, idx_list, PR):
    d = np.load(BASE + f"cells_{tag}.npz")
    P = PR[~np.isin(PR, np.array(Q))]
    pf = P.astype(float)
    ell = np.where(P == 2, math.log(2), -np.log1p(-2.0 / np.maximum(pf, 3)))
    Lpre = np.concatenate(([0.0], np.cumsum(ell)))            # Lpre[i] = sum of first i primes of P
    H = 2 * np.log(np.log(pf)) - Lpre[:-1]                     # w = P[i]: 2 loglog w - L(w), L(w)=sum p<w
    Hmax = np.maximum.accumulate(H)
    for i in idx_list:
        u, v, cert, s, jd, lKa = (float(d["u"][i]), float(d["v"][i]), float(d["cert"][i]), float(d["s"][i]),
                                  int(d["jd"][i]), float(d["logK"][i]))
        delta = (26 + jd) / 128
        a = int(np.searchsorted(pf, u, side="right"))          # primes <= u
        b = int(np.searchsorted(pf, v, side="left"))           # primes < v
        # K_cel
        cands = [(a, math.log(u))] + [(j + 1, math.log(pf[j])) for j in range(a, b)]   # (#primes included, log z)
        lK = 0.0
        for (n_inc, lz) in cands:
            if n_inc == 0: continue
            val = Lpre[n_inc] - 2 * math.log(lz) + Hmax[n_inc - 1]
            lK = max(lK, val)
        lK_c = lK + 1e-9                                       # folga
        Lv = math.fsum(ell[:b]) + 1e-9
        Ev = math.fsum(np.log1p(8 * pf[:b] ** (-delta))) + 1e-9 * (1 + b * 1e-6)
        N = mp.sqrt(8) * mp.mpf(u) ** 16
        logX = mp.log((N - 1) / 3) if Q == (3,) else mp.log((N - 1) / 5 - 3)
        mI = logX - Lv + logF(s, 2 + lK_c) - mp.log(2 * C) - mp.mpf(s) * delta * mp.log(v) - Ev
        mA = logX - Lv + logF(s, 2 + lKa) - mp.log(2 * C) - mp.mpf(s) * delta * mp.log(v) - Ev
        print(f"{tag} cell#{i} ({u:.6f},{v:.6f}] s={s:.5f} d={delta:.5f} #cand={len(cands)}")
        print(f"   K_cel check={math.exp(lK):.9f}  K main={math.exp(lKa):.9f}  (main>=check: {lKa >= lK})")
        print(f"   margin check (K check) = {float(mI):+.9f} | check w/ K main = {float(mA):+.9f} | cert main = {cert:+.9f}"
              f" | main<=check(K main)+1e-7: {cert <= float(mA) + 1e-7}")


if __name__ == "__main__":
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 20261017
    PR = sieve(10 ** 8 + 1)
    rng = np.random.default_rng(seed)
    d = np.load(BASE + "cells_Q35_interval.npz"); n = len(d["u"])
    lowmask = np.flatnonzero(d["v"] <= 1e6); highmask = np.flatnonzero(d["u"] >= 1e6)
    picks = [int(rng.choice(lowmask)), int(rng.choice(lowmask)), int(rng.choice(highmask)), int(np.argmin(d["cert"]))]
    print("seed", seed, "picks Q35_interval", picks)
    check("Q35_interval", (3, 5), 11, picks, PR)
    d3 = np.load(BASE + "cells_Q3_interval.npz")
    m3 = np.flatnonzero(d3["v"] <= 1e8)
    p3 = [int(np.argmin(d3["cert"])), int(rng.choice(m3))]
    print("picks Q3_interval", p3)
    check("Q3_interval", (3,), mp.mpf(19) / 6, p3, PR)
    # envelope: margem do autor com K envelope armazenado (K so verificado >= K_cel)
    de = np.load(BASE + "cells_Q3_envelope.npz")
    pe = [int(np.argmin(de["cert"]))]
    print("picks Q3_envelope", pe)
    check("Q3_envelope", (3,), mp.mpf(19) / 6, pe, PR)
