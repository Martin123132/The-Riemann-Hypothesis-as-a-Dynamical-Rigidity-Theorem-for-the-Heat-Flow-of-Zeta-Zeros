# Hardy Exact-Binary128 Point-Ball Defect Gate

Date: 2026-08-06
Status: rigorous saved-point source-expression certificate; not a proof and no uniform analytic error theorem

## Exact Inputs

The telemetry records each proof-relevant `real(kind=16)` value as a 32-hex-digit payload. The gate decodes each payload as an exact IEEE-754 binary128 dyadic rational and constructs the corresponding exact Arb point. Decimal telemetry is not used for arithmetic.

The phase constant remains source-derived: `p=4*ATAN(1)`, `tpp=2*p`, `tpm=-tpp`. Its saved binary payload is decoded directly; no separate value of pi is inserted.

## Certificate

At `96` and `160` decimal digits, all `64` low/high precision defect balls overlap, and all `64` high-precision balls exclude zero.

The smallest rigorous defect-magnitude lower endpoint is `0.0323006316806926616947992985969913945123004819`. The largest upper endpoint for the gap between exact-real Arb evaluation and the logged binary128 computation is `4.82442692200603594697156922365979315831003267e-31`. The largest gap-to-defect ratio upper endpoint is `6.65954162642018133423951821226155835735435017e-30`.

This establishes that the observed source-convention local discrepancy is genuine at the saved points and is not a binary128 accumulation artifact. It strengthens the case for enclosing the whole finite defect directly, where the W1-W5 cancellation is retained.

## Next Target

Promote coefficients, multiplier, and analytic qq from exact point payloads to outward-rounded functions on a real cell or fixed complex disk, while proving the transformed-level adapter and every selector margin.

## Proof Boundary

This gate rigorously encloses a source-convention finite expression at 64 saved exact binary128 input points and proves those point defect balls exclude zero. It does not prove that the source-convention transformed levels equal the paper's exact equation-120 objects, enclose the analytic W1-W5 corrections, control any real interval or complex disk, exclude selector or hierarchy transitions, bound the outer Hardy representation error, establish physical-height complexity, or imply a carrier value, Lambda<=0, PF-infinity, RH, or a prize-level conclusion.
