# Hardy Chain Telemetry Scout

Date: 2026-08-06
Status: finite low-height diagnostic; not a proof and no external error theorem

## Evaluator Boundary

The accepted resumable evaluator remains byte-pinned. The separately compiled telemetry derivative matches its full-precision checkpoint journal after removing only the run provenance id, and enabling telemetry leaves that derivative's journal and displayed values exactly unchanged.

The phase constant is not inserted as an unexplained `pi`. The accepted source computes `p = 4*atan(1)`, then `tpm = -2*p`; each chain records that actual `tpm`, and every independent reconstruction consumes the recorded value.

## Independent Reconstruction

The scout contains `64` chains and `64` recurrence rows. Every finite transformed level is re-summed independently at 50, 80, 120 decimal digits from its logged length and coefficients.

Initial lengths range from `104` to `104` terms before the inclusive endpoint adjustment; kernel lengths range from `1` to `2`. The maximum 80-to-120-digit ladder drift is `5.24035173954193484932302495016e-78`.

The conjugation-aware affine accumulation identity reconstructs the observed root error with maximum residual `2.24102370007601542184427111411e-33`.

## Route Comparison

- **Direct local defect:** observed norm summary `{"max": "0.156894051603196055299913871966", "median": "0.0802833246094686657057218036195", "min": "0.032300631680692661694799298597", "p90": "0.120779983562291635523604033461"}`. The W1-to-defect ratio summary is `{"max": "60.4384859580801730776741810526", "median": "7.71260735581735522807835758342", "min": "0.822322920146988894126499770064", "p90": "29.4858785286205303478489534799"}`. This is the leading certificate candidate because it preserves cancellations instead of bounding five source components separately.
- **Exact first-parent shell:** direct term range `[105, 105]` and exact elimination on `64` chains. Because every saved chain has `MIT=2`, this shell is the original 105-term sum; it validates the bypass identity but does not yet demonstrate a cheaper scalable shell.
- **Termwise W1:** observed `t5` summary `{"max": "3.95734872730712205026976064125", "median": "0.714338542587184886965713392359", "min": "0.064690256617508805640891124486", "p90": "2.18576831454898235027042134944"}`. It remains necessary as an analytic fallback, but these samples supply no uniform W1 bound.

Provisional order: direct finite-defect enclosure, first-parent exact shell, then termwise W1 where direct cost forces it.

## Proof Boundary

This scout proves evaluator non-interference for one accepted low-height fixture, independently reconstructs its finite transformed sums at an increasing precision ladder, verifies the conjugation-aware finite defect accumulation numerically, and measures exact-shell costs. It does not prove a transformed-level interval adapter, uniform W1 or special-function bounds, selector/transition exclusion, outer Hardy representation error, physical-height complexity, evaluator C6 or disk control, a carrier value, Lambda<=0, PF-infinity, RH, or a prize-level conclusion.
