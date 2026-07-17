# Jensen-Window PF Newman Polymath-15 Lambda 0.1965 Provenance Audit

Date: 2026-07-17

Status: published-bound provenance and candidate-certificate audit. The
candidate is reproduced in part but remains unrefereed. This is not a
proof of `Lambda <= 0` or RH.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_lambda01965_provenance_audit.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_lambda01965_provenance_audit.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_lambda01965_provenance_audit.py
```

## Established Record

Rodgers-Tao prove `Lambda >= 0`. Polymath15 proved `Lambda <= 0.22`,
and Platt-Trudgian Theorem 1 rigorously verifies RH through the exact
height `3000175332800`. Their Corollary 2 gives

```text
0 <= Lambda <= 1/5 = 0.2
```

This is the peer-reviewed published interval used by the corpus.

## Candidate Record

The Zenodo deposit `10.5281/zenodo.20724170`, dated
`2026-06-17`, claims

```text
t0=177/1000
y0^2=39/1000
t0+y0^2/2=393/2000=0.1965
X=6000000185827
T_PT-X/2=350479773/2
```

The PDF calls itself `DRAFT v0.3 - not for distribution`; the bundle
README says `PRE-FREEZE DRAFT`, and the deposit is a preprint. The paper
also discloses that the result and certificate campaign are AI-derived.

## Integrity And Rerun

```text
archive MD5=da7ab6a1bd1e43bb26d080704927cba1
manifest entries=505
manifest bytes hashed=213197118
manifest bad=0
manifest missing=0
portable invocations=22
portable pass=22
portable fail=0
portable measured seconds=271.18
```

The local pass covered both assembly lines, both binding lines, all
three normalization seams and their second lines, both site-glue lines,
both Euler-3 tail lines, the three block-uniform finite-box modes, the
y-reduction theorem, and the criterion/error-term spines.

## Nonpromotion Boundary

The following compiled producer packages were hash-checked but not
live-rerun on this Windows machine:

```text
record/record_package_197
certified1965/grid_full
certified1965/grid_tbox
```

Nor has this audit supplied an independent line-by-line review linking
every implementation formula to Polymath15. The correct corpus label is
therefore `reproduced unrefereed candidate`, not established record.

## Consequence

If accepted, the candidate changes only

```text
[0,1/5] -> [0,393/2000]
absolute shrink=7/2000
relative shrink=7/400
```

It does not alter the qualitative target: rule out every positive
boundary collision. No `Lambda <= 0` or RH conclusion is supplied.
