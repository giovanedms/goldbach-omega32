"""Independent cross-check of blocks.npz against the primesieve C++ CLI (~/src/primesieve/build/primesieve).
Counts for 60 random blocks + all blocks of 4 decade-spanning unions; exact sums for 12 random blocks (primesieve -p)."""
import subprocess, numpy as np, os
PS = os.path.expanduser("~/src/primesieve/build/primesieve")
d = np.load("blocks.npz"); b, n, s1 = d["bounds"], d["n"], d["s1"]
rng = np.random.default_rng(12345)
def ps_count(lo, hi):          # primes in [lo, hi)  -> primesieve counts [start, stop] inclusive
    out = subprocess.run([PS, str(lo), str(hi - 1), "-c"], capture_output=True, text=True).stdout
    return int([l for l in out.splitlines() if l.startswith("Primes")][0].split()[-1])
bad = 0
for j in rng.choice(len(n), 60, replace=False):
    c = ps_count(int(b[j]), int(b[j + 1])); bad += c != n[j]
for (x, y) in [(10**6, 10**7), (10**7, 10**8), (10**8, 10**9), (10**9, 10**10)]:
    i, k = np.searchsorted(b, x), np.searchsorted(b, y)
    c = ps_count(x, y); bad += c != n[i:k].sum()
    print(f"[{x:.0e},{y:.0e}) primesieve {c}  blocks {n[i:k].sum()}")
for j in rng.choice(len(n), 12, replace=False):
    out = subprocess.run([PS, str(int(b[j])), str(int(b[j + 1]) - 1), "-p"], capture_output=True, text=True).stdout
    s = sum(int(t) for t in out.split()); bad += s != int(s1[j])
print("mismatches:", bad)
assert bad == 0
