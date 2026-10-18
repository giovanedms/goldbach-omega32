"""Re-run the certificate in a scratch directory and compare every output with the shipped one.

usage: python certificate/verify_all.py [--quick] [--keep]
  --quick  skip the slowest step (local_sieve/test_local_sieve.py, about 1.5 min)
  --keep   keep the scratch directory (its path is printed)

The shipped files are never touched: every script runs with the scratch directory as the
current directory. Comparisons are exact, except that
  * run times ("# total 40s", "build 43s", "certification 2s") are masked;
  * Python warning lines ("test_local.py:60: RuntimeWarning ...", plus the source line
    that follows) are dropped: they go to stderr and their position in the log varies;
  * out_analytic.txt was written before analytic.py printed "(A-s)" instead of "(16-s)"
    (A = 16 for odd N); the label is normalised, every number is compared as is.
check_blocks.py runs only if the primesieve CLI is found (PATH or PRIMESIEVE=...).
Exit status 0 iff SHA256SUMS verifies and every comparison passes.  Total time about 5 min.
"""
import hashlib, os, re, shutil, subprocess, sys, tempfile, time
import numpy as np

sys.stdout.reconfigure(line_buffering=True)

CERT = os.path.dirname(os.path.abspath(__file__))
RIG = os.path.join(CERT, "rigoroso")
PY = sys.executable
TIME = re.compile(r"\b\d+(\.\d+)?s\b")
WARN = re.compile(r"^\S*\.py:\d+: \w*Warning:.*\n(?:[ \t]+.*\n)?", re.M)
fails = []


def check(name, ok):
    print(f"  {'OK  ' if ok else 'FAIL'} {name}")
    if not ok:
        fails.append(name)


def norm(text, label=False):
    text = TIME.sub("<t>", WARN.sub("", text))
    return text.replace("(A-s)", "(16-s)") if label else text


def run(args, out=None, env=None):
    t = time.time()
    r = subprocess.run([PY, *args], cwd=WORK, capture_output=True, text=True,
                       env={**os.environ, **(env or {})})
    if out:
        with open(os.path.join(WORK, out), "w") as f:
            f.write(r.stdout + r.stderr)
    print(f"  ran  {' '.join(os.path.relpath(a, CERT) if a.startswith(CERT) else a for a in args)}"
          f"{' [' + ' '.join(f'{k}={v}' for k, v in env.items()) + ']' if env else ''}  ({time.time() - t:.0f}s)")
    check(f"exit status 0: {os.path.basename(args[0])} {' '.join(args[1:])}", r.returncode == 0)
    return r


def same_text(new, old, label=False):
    a = open(os.path.join(WORK, new)).read()
    b = open(old).read()
    check(f"{os.path.relpath(old, CERT)}", norm(a, label) == norm(b, label))


def same_npz(new, old):
    a, b = np.load(os.path.join(WORK, new)), np.load(old)
    check(f"{os.path.relpath(old, CERT)}", set(a) == set(b) and all(np.array_equal(a[k], b[k]) for k in a))


def same_bytes(new, old):
    h = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
    check(f"{os.path.relpath(old, CERT)}", h(os.path.join(WORK, new)) == h(old))


if __name__ == "__main__":
    quick, keep = "--quick" in sys.argv, "--keep" in sys.argv
    WORK = tempfile.mkdtemp(prefix="goldbach_verify_")
    print(f"scratch directory: {WORK}")

    print("[1] hashes")
    r = subprocess.run(["sha256sum", "-c", "--quiet", "SHA256SUMS"], cwd=CERT, capture_output=True, text=True)
    check("SHA256SUMS", r.returncode == 0)

    print("[2] blocks.npz against primesieve")
    if os.environ.get("PRIMESIEVE") or shutil.which("primesieve"):
        run([os.path.join(RIG, "check_blocks.py")], "out_check_blocks.txt")
        same_text("out_check_blocks.txt", os.path.join(RIG, "out_check_blocks.txt"))
    else:
        print("  skip primesieve CLI not found")

    print("[3] finite cover")
    covers = [("3,5", "1e8", "interval", {}), ("3", "1e10", "interval", {}), ("3,5", "1e8", "envelope", {}),
              ("3", "1e8", "uniform:1.146", {}), ("3,5", "1e8", "uniform:sup", {}),
              ("3,5", "1e8", "interval", {"GB_PARITY": "even"})]
    for q, x, mode, env in covers:
        run([os.path.join(RIG, "run_cover.py"), q, x, mode], env=env)
        tag = f"Q{q.replace(',', '')}_{mode.replace(':', '')}{'_even' if env else ''}"
        same_text(f"out_cover_{tag}.txt", os.path.join(RIG, f"out_cover_{tag}.txt"))
        if not mode.startswith("uniform"):
            same_npz(f"cells_{tag}.npz", os.path.join(RIG, f"cells_{tag}.npz"))
    for j in ["data_Q35.json", "data_Q3.json", "data_Q35_even.json"]:
        same_bytes(j, os.path.join(RIG, j))

    print("[4] analytic range")
    odd = [("3,5", "1e8"), ("3,5", "1e8", "1.08997"), ("3", "1e8"), ("3", "1e9"), ("3", "1e10"),
           ("3", "1e10", "1.146"), ("3", "1e9", "1.097"), ("3", "1e10", "1.097")]
    text = "".join(run([os.path.join(RIG, "analytic.py"), *a]).stdout for a in odd)
    open(os.path.join(WORK, "out_analytic.txt"), "w").write(text)
    same_text("out_analytic.txt", os.path.join(RIG, "out_analytic.txt"), label=True)
    run([os.path.join(RIG, "analytic.py"), "3,5", "1e8"], "out_analytic_even.txt", env={"GB_PARITY": "even"})
    same_text("out_analytic_even.txt", os.path.join(RIG, "out_analytic_even.txt"))

    print("[5] local lemmas, float vs rigorous")
    run([os.path.join(RIG, "test_local.py")], "out_test_local.txt")
    same_text("out_test_local.txt", os.path.join(RIG, "out_test_local.txt"))

    print("[6] separate re-verification (AI-written)")
    for s in ["local", "cells", "rs"]:
        run([os.path.join(CERT, "crosscheck", f"crosscheck_{s}.py")], f"out_crosscheck_{s}.txt")
        same_text(f"out_crosscheck_{s}.txt", os.path.join(CERT, "crosscheck", f"out_crosscheck_{s}.txt"))

    if not quick:
        print("[7] algebra of Appendix A")
        run([os.path.join(CERT, "local_sieve", "test_local_sieve.py"), "300"], "out_test_local_sieve.txt")
        same_text("out_test_local_sieve.txt", os.path.join(CERT, "local_sieve", "out_test_local_sieve.txt"))

    if not keep:
        shutil.rmtree(WORK)
    print(f"\n{'ALL CHECKS PASSED' if not fails else f'{len(fails)} FAILED: ' + '; '.join(fails)}")
    sys.exit(1 if fails else 0)
