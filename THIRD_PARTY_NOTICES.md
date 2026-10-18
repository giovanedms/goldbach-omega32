# Third-party notices

No source file in this repository was copied or derived from the code of
P. J. Campbell, "Verification scripts for explicit almost-prime computations"
(https://github.com/PeterJCampbell1/almost-prime-sum-verifications, release
`arxiv-v1`, commit `2755a58b6eb13c0b43fb55826c7bf26e403d80ec`, MIT licence).

The directory `certificate/campbell_repro/` contains only the **text outputs and
timings** produced by running the three scripts of that repository
(`K_presieved`, `V_presieved_lower_bound`, `case2_presieved`), as a regression
check of the Omega(ab) <= 33 framework. The scripts themselves are not included.
If any of Campbell's scripts are added to this repository in the future, his
copyright notice and the full MIT licence text from his repository must be
kept alongside them.

The paper (`paper/`) builds on P. J. Campbell, "A sharper explicit result for
the sum of two almost primes", arXiv:2608.09488, and cites it in `refs.bib`.

Python dependencies (python-flint, numpy, mpmath, scipy, matplotlib) and the
optional `primesieve` tool are used as installed packages under their own licences.
