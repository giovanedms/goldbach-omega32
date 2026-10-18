# goldbach-omega32

Code, outputs and paper for **"Goldbach-type sums a+b=N with Ω(ab) ≤ 32"**
(Giovane de Morais, Instituto Tecnológico de Aeronáutica, ITA).

## Result

**Theorem.** Every integer N ≥ 2 can be written as N = a + b with positive
integers a, b and **Ω(ab) ≤ 32**, where Ω counts prime factors with multiplicity.

This improves the bound Ω(ab) ≤ 33 of P. J. Campbell, *A sharper explicit result for
the sum of two almost primes*, [arXiv:2608.09488](https://arxiv.org/abs/2608.09488),
whose verification scripts are at
<https://github.com/PeterJCampbell1/almost-prime-sum-verifications>.
The proof stays inside Campbell's framework; the new ingredients are a per-interval
dimension constant, per-interval sieve parameters and pre-sieving by 3 and 5. The
intermediate range (117 415 intervals for odd N, 118 322 for even N) is certified in
ball arithmetic by the code here. For completeness, Appendix A of the paper gives a
self-contained proof of the explicit Friedlander–Iwaniec lower bound with the dimension
condition imposed only at the sifting level; `certificate/local_sieve/test_local_sieve.py`
checks the algebra of that proof numerically (partly in exact rational arithmetic, partly
in floating point).

The paper is a **preprint, submitted to the Journal of Number Theory** (Computational Section) on 18 October 2026; it has not been peer reviewed.
Archived at Zenodo: [10.5281/zenodo.23456153](https://doi.org/10.5281/zenodo.23456153) (all versions); release v1.0.0: [10.5281/zenodo.23456154](https://doi.org/10.5281/zenodo.23456154).

## Contents

| Path | What it is |
|---|---|
| `paper/` | `main.tex`, `main.pdf`, `refs.bib`, `figs/` (PDF figures and `make_figs.py`) |
| `certificate/rigoroso/` | The certificate: exact prime-block data (`primes_blocks.py`, `blocks.npz`), ball-arithmetic library `rlib.py`, finite cover `run_cover.py`, analytic range `analytic.py`, local-lemma tests `test_local.py`, outputs `out_*.txt`, `data_*.json`, `cells_*.npz` |
| `certificate/crosscheck/` | Separate re-verification (AI-written, internal consistency) (`crosscheck_local.py`, `crosscheck_cells.py`, `crosscheck_rs.py`) and outputs |
| `certificate/gate1/` | Separate float64 reimplementation (AI-written, internal consistency), outputs |
| `certificate/local_sieve/` | `test_local_sieve.py`, numerical check of the algebra of Appendix A (Theorem 4.1), and its output |
| `certificate/rankin/` | `rk.py`, numeric check of the Rankin step (paper Proposition 4.4) |
| `certificate/campbell_repro/` | Outputs and timings of Campbell's three scripts (regression for Ω ≤ 33). His code is not included, see `THIRD_PARTY_NOTICES.md` |
| `certificate/SHA256SUMS` | SHA-256 of every file in `certificate/` |

## Requirements

- Python ≥ 3.12 (developed with 3.14.4)
- `python-flint` (0.9.0), `numpy` (2.5.4), `mpmath` (1.3.0, `crosscheck/`), `scipy` (1.18.1, `gate1/`), `matplotlib` (3.11.2, figures only); exact pins in `requirements.txt`
- [`primesieve`](https://github.com/kimwalisch/primesieve) CLI, optional, only for `check_blocks.py` (tested with 12.17); it is looked up on `PATH`, or set `PRIMESIEVE=/path/to/primesieve`

```bash
python -m venv .venv && . .venv/bin/activate
pip install -r requirements.txt
```

No GPU needed. Times below are on one CPU core of the author's workstation.

## Verifying everything in one command

```bash
python certificate/verify_all.py          # about 5 min; --quick skips Appendix A test (~1.5 min)
```

It checks `SHA256SUMS`, re-runs every step below in a temporary directory (the shipped
files are never touched) and compares each output with the shipped one: the cell files
(`cells_*.npz`) and `data_*.json` must be identical, the text outputs identical up to run
times and Python warning lines. It exits with status 0 only if every check passes. On a
fresh clone with a fresh virtual environment (18 October 2026) all checks passed; with
unpinned packages (`mpmath` 1.4.1) the outputs were also identical.

## Reproducing step by step

All scripts locate `blocks.npz` relative to their own path. Step 3 overwrites
`out_*.txt`, `data_*.json` and `cells_*.npz` in the **current directory**; to compare
without touching the shipped files, run from a scratch directory
(`python /path/to/certificate/rigoroso/run_cover.py ...`).

```bash
cd certificate/rigoroso
# 1) (optional) regenerate blocks.npz: segmented sieve up to 1e10 (~38 s, ~257 MB RAM)
python primes_blocks.py 1e10
# 2) check blocks.npz against primesieve (< 1 s; needs the primesieve CLI, see Requirements)
python check_blocks.py
# 3) finite cover
python run_cover.py 3,5 1e8 interval      # main chain, ~20-40 s
python run_cover.py 3 1e10 interval       # alternative chain, ~43 s
python run_cover.py 3,5 1e8 envelope      # scenario E, ~60 s
python run_cover.py 3 1e8 uniform:1.146   # scenario U (expected to fail), ~40-55 s
python run_cover.py 3,5 1e8 uniform:sup   # scenario U, sup of the per-cell K (expected to fail), ~40 s
GB_PARITY=even python run_cover.py 3,5 1e8 interval   # even N, ~70 s (writes *_even outputs)
# 4) analytic range (reads data_Q*.json from the current directory, <1 s each).
#    The two runs used by the proof:
python analytic.py 3,5 1e8                       # odd N
GB_PARITY=even python analytic.py 3,5 1e8        # even N  -> out_analytic_even.txt
#    out_analytic.txt is the concatenated output of these eight runs (Table 2; an optional
#    third argument K replaces the dimension constant by max(K, K_an(Y))):
for a in "3,5 1e8" "3,5 1e8 1.08997" "3 1e8" "3 1e9" "3 1e10" "3 1e10 1.146" "3 1e9 1.097" "3 1e10 1.097"; do
  python analytic.py $a; done
# 5) brute force of the local lemmas and float vs. rigorous comparison (~40 s)
python test_local.py
# 6) separate (AI-written) re-verification (~10 s)
python ../crosscheck/crosscheck_local.py && python ../crosscheck/crosscheck_cells.py && python ../crosscheck/crosscheck_rs.py
# 7) numerical check of the algebra of Appendix A (~2 min)
python ../local_sieve/test_local_sieve.py 300
```

The `uniform` scenarios are meant to fail (they show that a single global constant
does not suffice); their `.npz` cell files are not shipped, only the `.txt` outputs.

`out_analytic.txt` was written before `analytic.py` printed the label `(A-s)` (with
A = 16 for odd N, 16.5 for even N) instead of `(16-s)`; every number in it is reproduced
exactly by the current code. It is kept as is because its SHA-256 is quoted in the paper.

Auxiliary checks (not part of the proof; float64, outputs in the same directories):

```bash
cd certificate
python gate1/cover_cells.py 3,5 1e-4 && python gate1/cover_cells.py 3 1e-4   # ~1 min each
python gate1/analytic_check.py && python gate1/sanity_paper.py
python rankin/rk.py                                                          # Proposition 4.4, numerically
python profiles/profile_general.py 3 1e8 "8**0.5" 16 "odd M=32"              # one line of out_profiles.txt
```

Their shipped outputs hold stdout only; NumPy `RuntimeWarning` lines (division by zero
at p = 2, handled by `np.where`) go to stderr.

Figures: `cd paper/figs && python make_figs.py` (reads `certificate/rigoroso/`).
Paper: `cd paper && pdflatex main && bibtex main && pdflatex main && pdflatex main`
(class `elsarticle`).

## Verifying hashes

```bash
cd certificate
sha256sum -c SHA256SUMS
```

Two shipped outputs (`out_primes_blocks.txt`, `out_test_local.txt`) had the
author's local absolute paths replaced by neutral ones before publication. Their
hashes in `SHA256SUMS` refer to the published versions. The numerical content is unchanged.

## Licence and citation

Code: MIT (`LICENSE`). Paper text and figures: all rights reserved by the author
unless stated otherwise. See `CITATION.cff` and `THIRD_PARTY_NOTICES.md`.

## Result tables

`results/` contains the numerical tables of the paper as CSV, extracted from the certificate
outputs by `python results/extract_tables.py` (run from the repository root):

| File | Paper | Content |
|---|---|---|
| `table1_decades.csv` | Table 1 | minimal certified margin per decade, Q = {3,5}, per-cell K |
| `table2_analytic.csv` | Table 2 | analytic range: constants, s, alpha and certified M_lo(log Z) |
| `table3_margins.csv` | Table 3 | minimal certified margins of the finite range, scenarios I / E / U |
| `table4_omega31.csv` | Table 4 | exploratory floating-point profiles for Omega(ab) <= 31 (not certified) |

`certificate/profiles/` holds the script and output of the Omega(ab) <= 31 profiles
(`python profile_general.py Q LIM c r label`, float64, a grid of 60 points).
