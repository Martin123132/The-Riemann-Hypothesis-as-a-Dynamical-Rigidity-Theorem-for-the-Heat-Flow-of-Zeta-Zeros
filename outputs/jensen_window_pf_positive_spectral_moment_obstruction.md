# Jensen-Window PF Positive Spectral Moment Obstruction

Date: 2026-07-06

Status: ordinary positive spectral moment obstruction. This is not a proof of
the all-order column recurrence, Schur positivity, Jensen-window PF-infinity,
Jensen hyperbolicity, RH, or `Lambda <= 0`.

## Purpose

The positive-readout target keeps a spectral-transform subroute alive. This
note separates an impossible interpretation from the still-live one.

For:

```text
H(t) = 1 + g_1*t + g_2*t^2 + ...
E(t) = 1/H(-t)
mu_m = [t^m]E(t)
```

we have:

```text
mu_0 = 1
mu_1 = g_1
mu_2 = g_1^2 - g_2
Delta_2 = det [[mu_0, mu_1], [mu_1, mu_2]] = -g_2
```

So if `g_2>0`, the sequence `mu_m` cannot be ordinary power moments of a
positive measure.

Machine-readable obstruction:

```text
work/rh_compute/results/jensen_window_pf_positive_spectral_moment_obstruction.json
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_positive_spectral_moment_obstruction.py
```

Current result:

```text
validated Jensen-window PF positive spectral moment obstruction: 3 symbolic rows, 735 finite Delta2 obstruction rows, 0 issues
```

## Symbolic Rows

```text
psm_01_degree2_hankel_obstruction:
  Delta_2 = -g_2, so any window with g_2>0 fails the 2x2 positive moment
  Hankel condition.

psm_02_signed_hankel_signature_conflict:
  The signed J-fraction finite grid records (-1)^(r(r-1)/2) Delta_r>0. At
  r=2 this is Delta_2<0, incompatible with ordinary positive moments.

psm_03_transform_escape_hatch:
  A spectral transform can remain live only if it is not the ordinary moment
  representation mu_m=int x^m dnu(x) of a positive measure.
```

## Finite Evidence

The finite check reads:

```text
outputs/jensen_window_pf_reciprocal_signed_j_fraction_scout.md
work/rh_compute/results/jensen_window_pf_reciprocal_signed_j_fraction_scout.json
work/rh_compute/results/jensen_window_pf_reciprocal_signed_j_fraction_lamgrid_n0_n20_d2_d8_dps520.jsonl
```

It validates:

```text
735 finite Delta_2 obstruction rows
expected_delta_sign = -1
all checked Delta_2 rows negative
```

The source scout validates:

```text
validated Jensen-window PF reciprocal signed J-fraction scout: 2 symbolic rows, 3675 signed Hankel rows, 2940 signed-lambda rows, 0 issues
```

## Consequence

Reject:

```text
mu_m = int x^m dnu(x), where nu is a positive measure
```

as a proof route for the raw reciprocal coefficients.

Still live:

```text
nonordinary positive transform
positive kernel identity
non-power basis coefficient extraction
Xi/Phi-specific positive scalar functional
```

provided the final Taylor coefficients are exactly `mu_m` and the proof does
not assume endpoint PF, Jensen hyperbolicity, Laguerre-Polya membership, RH,
or `Lambda <= 0`.

## Boundary

Passing this checker means the ordinary positive moment version of the
positive spectral transform has been rejected. It does not reject nonordinary
positive transforms, positive kernels, Xi/Phi-specific readouts, or any future
exact positive scalar functional whose coefficient extraction is not ordinary
power moments of the raw `mu_m` sequence.
