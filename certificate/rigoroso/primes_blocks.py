"""Exact prime data for the rigorous T3 verification.

  * primes_below(n): all primes < n (numpy odd-only Eratosthenes), exact.
  * block boundaries b_0 = P0 < b_1 < ... (integers, geometric spacing, anchors 1e6..1e10 included exactly)
  * for each block [b_j, b_{j+1}): n_j = #primes, S1_j = sum of primes  (exact int64; segmented sieve)
Output: blocks.npz.  Cross-checks: pi(10^k) against known values; blocks also checked against the
primesieve CLI (independent C++ implementation) in check_blocks.py.

usage: primes_blocks.py [XMAX=1e10]
"""
import sys, math, time
import numpy as np

P0 = 10**6
ANCHORS = [10**6, 10**7, 10**8, 10**9, 10**10]
EPS = {10**6: 2e-4, 10**7: 5e-4, 10**8: 1e-3, 10**9: 1e-3}   # relative spacing on [anchor, next anchor)


def primes_below(n):
    """All primes p < n, int64."""
    if n <= 2:
        return np.zeros(0, np.int64)
    m = (n + 1) // 2                      # index i <-> 2i+1 ; need 2i+1 < n
    s = np.ones(m, dtype=bool); s[0] = False
    for i in range(1, (math.isqrt(n) - 1) // 2 + 1):
        if s[i]:
            p = 2 * i + 1
            s[p * p // 2::p] = False
    odd = 2 * np.flatnonzero(s).astype(np.int64) + 1
    odd = odd[odd < n]
    return np.concatenate((np.array([2], np.int64), odd))


def make_bounds(xmax):
    b = []
    for a, c in zip(ANCHORS[:-1], ANCHORS[1:]):
        if a >= xmax:
            break
        e = EPS[a]
        g = int(math.ceil(math.log(c / a) / math.log1p(e)))
        x = np.round(a * (1 + e) ** np.arange(g)).astype(np.int64)
        b.extend(x.tolist())
    b.append(min(xmax, ANCHORS[-1]))
    b = np.unique(np.array(b, np.int64))
    for a in ANCHORS:
        assert a > xmax or a in b, a
    return b


def block_stats(bounds, seg=10**8):
    lo_all, hi_all = int(bounds[0]), int(bounds[-1])
    base = primes_below(math.isqrt(hi_all) + 2)[1:]          # odd base primes
    nb = len(bounds) - 1
    cnt = np.zeros(nb, np.int64); s1 = np.zeros(nb, np.int64)
    lo = lo_all
    while lo < hi_all:
        hi = min(lo + seg, hi_all)
        # odd numbers in [lo, hi): o = lo0 + 2i
        lo0 = lo if lo % 2 == 1 else lo + 1
        m = (hi - lo0 + 1) // 2
        s = np.ones(m, dtype=bool)
        for p in base:
            p = int(p)
            if p * p >= hi:
                break
            st = max(p * p, ((lo0 + p - 1) // p) * p)
            if st % 2 == 0:
                st += p
            s[(st - lo0) // 2::p] = False
        pr = lo0 + 2 * np.flatnonzero(s).astype(np.int64)
        pr = pr[pr < hi]
        if lo <= 2 < hi:
            raise ValueError("block region must start above 2")
        cs = np.concatenate(([0], np.cumsum(pr)))
        pos = np.searchsorted(pr, bounds, side="left")            # #primes(seg) < b
        cnt += np.diff(pos)
        s1 += np.diff(cs[pos])
        lo = hi
    return cnt, s1


if __name__ == "__main__":
    xmax = int(float(sys.argv[1])) if len(sys.argv) > 1 else 10**10
    t = time.time()
    small = primes_below(P0)
    b = make_bounds(xmax)
    n, s1 = block_stats(b)
    print(f"blocks={len(b)-1} range [{b[0]},{b[-1]}) time {time.time()-t:.1f}s")
    known = {10**6: 78498, 10**7: 664579, 10**8: 5761455, 10**9: 50847534, 10**10: 455052511}
    cum = len(small) + np.concatenate(([0], np.cumsum(n)))
    for x, v in known.items():
        if x <= b[-1]:
            j = int(np.searchsorted(b, x)); got = int(cum[j]) if x > P0 else len(small)
            print(f"pi({x:.0e}) = {got}  known {v}  {'OK' if got == v else 'MISMATCH'}")
            assert got == v
    np.savez("blocks.npz", bounds=b, n=n, s1=s1)
    print("saved blocks.npz")
