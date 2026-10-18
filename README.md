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
ball arithmetic by the code here. Appendix A of the paper proves that the explicit
Friedlander–Iwaniec lower bound holds with the dimension condition imposed only at the
sifting level, which is the reading Campbell's argument requires;
`certificate/local_sieve/test_local_sieve.py` tests every step of that proof in exact
rational arithmetic.

The paper is a **preprint, to be submitted**. Status and DOI will be added here.

## Contents

| Path | What it is |
|---|---|
| `paper/` | `main.tex`, `main.pdf`, `refs.bib`, `figs/` (PDF figures and `make_figs.py`) |
| `certificate/rigoroso/` | The certificate: exact prime-block data (`primes_blocks.py`, `blocks.npz`), ball-arithmetic library `rlib.py`, finite cover `run_cover.py`, analytic range `analytic.py`, local-lemma tests `test_local.py`, outputs `out_*.txt`, `data_*.json`, `cells_*.npz` |
| `certificate/crosscheck/` | Independent re-verification (`crosscheck_local.py`, `crosscheck_cells.py`, `crosscheck_rs.py`) and outputs |
| `certificate/gate1/` | Independent float64 reimplementation, outputs |
| `certificate/local_sieve/` | `test_local_sieve.py`, exact rational test of every step of Appendix A (Theorem 4.1) |
| `certificate/rankin/` | `rk.py`, numeric check of the Rankin step (paper §4.4) |
| `certificate/campbell_repro/` | Outputs and timings of Campbell's three scripts (regression for Ω ≤ 33). His code is not included, see `THIRD_PARTY_NOTICES.md` |
| `certificate/SHA256SUMS` | SHA-256 of every file in `certificate/` |

## Requirements

- Python ≥ 3.12 (developed with 3.14.4)
- `python-flint` (0.9.0), `numpy` (2.5.4), `mpmath` (1.3.0, `crosscheck/`), `scipy` (1.18.1, `gate1/`), `matplotlib` (figures only)
- [`primesieve`](https://github.com/kimwalisch/primesieve) CLI, optional, only for `check_blocks.py` (tested with 12.17)

```bash
pip install python-flint numpy mpmath scipy matplotlib
```

No GPU needed. Times below are on one CPU core of the author's workstation.

## Reproducing

All scripts locate `blocks.npz` relative to their own path. Step 3 overwrites
`out_*.txt`, `data_*.json` and `cells_*.npz` in the **current directory**; to compare
without touching the shipped files, run from a scratch directory
(`python /path/to/certificate/rigoroso/run_cover.py ...`).

```bash
cd certificate/rigoroso
# 1) (optional) regenerate blocks.npz: segmented sieve up to 1e10 (~38 s, ~257 MB RAM)
python primes_blocks.py 1e10
# 2) check blocks.npz against primesieve (< 1 s)
python check_blocks.py
# 3) finite cover
python run_cover.py 3,5 1e8 interval      # main chain, ~20-40 s
python run_cover.py 3 1e10 interval       # alternative chain, ~43 s
python run_cover.py 3,5 1e8 envelope      # scenario E, ~60 s
python run_cover.py 3 1e8 uniform:1.146   # scenario U (expected to fail), ~40-55 s
GB_PARITY=even python run_cover.py 3,5 1e8 interval   # even N, ~70 s (writes *_even outputs)
# 4) analytic range (reads data_Q*.json from the current directory, <1 s each)
python analytic.py 3,5 1e8
python analytic.py 3 1e9
GB_PARITY=even python analytic.py 3,5 1e8
# 5) brute force of the local lemmas and float vs. rigorous comparison (~40 s)
python test_local.py
# 6) independent re-verification (~10 s)
python ../crosscheck/crosscheck_local.py && python ../crosscheck/crosscheck_cells.py && python ../crosscheck/crosscheck_rs.py
# 7) exact rational test of the sieve inequality of Appendix A (~1 min)
python ../local_sieve/test_local_sieve.py 300
```

The `uniform` scenarios are meant to fail (they show that a single global constant
does not suffice); their `.npz` cell files are not shipped, only the `.txt` outputs.
`certificate/gate1/cover_cells.py` (float64 reimplementation) takes about 3 min.

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
