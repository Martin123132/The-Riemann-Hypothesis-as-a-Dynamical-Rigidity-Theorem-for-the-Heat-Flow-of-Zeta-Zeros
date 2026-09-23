# Uniform approximation and a successful conditional-moment correction

Copyright (c) 2026 TWO HANDS NETWORK LTD. All rights reserved.
Private local research. 17 September 2026.

## Result

One approximation obstruction is removed analytically. A second obstruction
was diagnosed on the existing N=2980 arithmetic data and avoided by a fixed,
complete height partition. This is not a proof of RH or an unbounded sign theorem.

| Same N2980 retained polynomial and original gate | Certified score enclosure |
| --- | --- |
| Whole-window fixed-pair moment bound | (-0.541799, -0.0309742) |
| Same moment method in 128 fixed complete height blocks | (1.2072164, 1.5787755) |
| Direct gated transverse minimum | (1.5293473, 1.6178036) |

The direct column is a tighter enclosure of the existing retained-polynomial
quantity, not a newly improved physical margin. All 65536 height children
remain included. No source terms, gate-uncertain cells or negative pieces
were discarded. The finite test uses the full 2976-term source on all four
bands; it must not be confused with the N64 tapered hybrid or an asymptotic
multiscale-ladder execution.

## What changed mathematically

For every N>=64, the inherited rational-polynomial comparisons admit explicit
Bernstein refinement and unnormalized log-sum lower comparisons ell_j with

    0 <= b_j-ell_j < 5/d^2,  d=ceil(log(N+1)).

The complete moment functional changes by at most 5 B m/d^2, with B bands
and m the original gate mass. Under the inherited m=O(sqrt d) and B=O(log d)
inputs, this is O(log d/d^(3/2)), below the d^-1 collar scale. The construction
is expensive: about N d^4 controls per band. It is not an efficient all-height
algorithm. UNIFORM_DERIVATION.md proves the bound and states its hypotheses.

The N2980 test then shows that approximation accuracy is not the only issue:
the unsmoothed whole-window moment functional itself is negative. Its
dispersion penalty is greater than 1.64355, while the actual outer-to-minimum
competition is less than 0.04083. Negative gap means alone do not prevent
this loss.

Conditioning the moments on a fixed partition preserves more height
dependence. Cauchy-Schwarz proves that it can only lower the exact cap.
Here the cap drops from more than 1.66626 to less than 0.338528, restoring
a positive score. The 128 blocks were fixed before the refined results were
read, not selected for their signs. All blocks use the original global
normalization. Smoothing stability still costs 5 B m/d^2, with no extra
factor of 128. See CONDITIONAL_DERIVATION.md and PARTITION_PLAN.md.

## Verification performed

- Both complete refinement builds ran at 160 and 224 bits. Their 65536-child
  mathematical ledgers agree byte-for-byte; result payloads agree exactly.
- Independent Fraction reaggregation checks all 65536 children, all three
  first/second moments and all 128 conditional blocks.
- A separately coded power-basis transformation checks 81920 control
  intervals on 256 fixed children. Independent 75-digit mpmath.iv checks
  their weight and gate enclosures. This is a fixed subset of control/phase
  checks, not an independently recomputed source sum at every child.
- The earlier 8192-cell hard/smooth test and its inconclusive enclosures are
  retained. Its two Arb builds agree, with 256 independent smooth-band checks.
- All 18 new named test families passed, in addition to the uniform theorem's
  150 rational-panel controls, 1000 stability controls and exact constants.
- Source jets, third-derivative bounds and their original provenance remain
  inherited. No Dirichlet or auxiliary-source sum was reevaluated this turn.

The delivery's CLEAN_REPLAY_RECEIPT.json records fresh extraction and complete
224-bit refinement reconstruction with independent acceptance. It is external
to the immutable package so its creation cannot alter the package seal.

## What remains

The next theorem must bound the ORIGINAL shared-height arithmetic in the
conditional first and second moments, on a defensible unbounded family.
It must control partition complexity and exceed every attached source,
collar, ladder and remaining payment. Merely making the partition arbitrarily
fine is not that theorem. One failed finite global-moment score does not
exclude success on an unbounded subsequence.

The recommended next step is a source-derived, phase-local moment estimate
for a deterministic partition, retaining the exact middle transverse region
and original gate. No automatic height scan, Forge study or publication is
queued. No RH proof or new zeta-zero count is claimed.

## Preservation

Historical packages are read-only. The session checkpoint is updated only
after delivery replay, with all 1333 inherited pins and previous checkpoint
objects checked unchanged. GitHub, the desktop input folder and unrelated
jobs are untouched. No provider calls or subagents were used.
