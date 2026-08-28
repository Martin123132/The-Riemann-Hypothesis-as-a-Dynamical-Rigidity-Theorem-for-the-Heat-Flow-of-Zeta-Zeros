# Mordell parameter boxes and direct phase-scale audit

Date: 2026-08-24

Status: interval-certificate and method-barrier note; central Mordell boxes,
recurrence-wall splitting, and two microscopic complete-current boxes are
certified, but scalable recursion and physical quadrature remain open

## Parameter boxes

Let `Z=[z_-,z_+]` lie in `[-1/2,1/2]` and let
`T=[tau_-,tau_+]` have `tau_->0`.  Substituting the Arb balls `Z,T` directly
in Kuznetsov's rotated integral encloses `h(Z,T)` and `h_z(Z,T)` on `[0,Y]`.
The Section 11.464 tails are evaluated at `tau_-`, where they are largest:

```text
T_h <=coth(pi Y/sqrt(2))*erfc(sqrt(pi tau_-)Y)/sqrt(tau_-),
T_hz<=2*coth(pi Y/sqrt(2))*exp(-pi tau_- Y^2)/tau_-. (PB1)
```

Target boxes crossing a half integer are split there.  On each piece one
fixed exact unit recurrence transports the central enclosure; the branch
balls are then joined by interval hull.  Thus a recurrence wall is covered,
not sampled around or silently assigned to one side.                    (PB2)

Production at 320 bits,
`Y=10`, and tolerance
`1e-55` certifies
3 target boxes, including one crossing `z=1/2`.
All 27 rigorous pointwise probes are
contained in their parent boxes.  The largest target-box radii are

```text
h:   0.01045107631944119930267333984375,
h_z: 0.0419043576694093644618988037109375.       (PB3)
```

The physical map also certifies two joined endpoint-current boxes while
checking that `r=nearest(x(A+2s)/2)` and `m=floor(Kx)` remain fixed.  The
corner box crosses the Mordell recurrence wall and is split automatically.

## Complete-current microboxes

For fixed `r,m`, write

```text
T(w,sigma)=sum_(k=0)^m (r+k)exp(2pi i(wk+sigma k^2)).
```

Around an exact rational centre `(w_0,sigma_0)`, the elementary chord bound
`|exp(iu)-exp(iv)|<=|u-v|` gives

```text
|T-T_0|<=2pi(delta_w S_1+delta_sigma S_2),
S_1=sum_(k=0)^m(r+k)k,
S_2=sum_(k=0)^m(r+k)k^2.                             (PB4)
```

The centre `T_0` is enclosed by exact rational-period compression.  The
prefactor is treated as `q exp(i pi phi)`, where

```text
q=2/x^(3/2),  w=(A+2s)/2-r/x,
sigma=-1/(2x),  phi=1/4+r(A+2s)-r^2/x.              (PB5)
```

Combining (PB4) with interval endpoint currents preserves the exact centre
cancellation and bounds only its variation.  At half-width `1e-28`, both
full-roster boxes contain their independently enclosed exact centre source
currents.  Their largest absolute complete radius is
`0.0067461658982210792601108551025390625`.

## Scaling decision

The ten-row width audit applies exactly the same proved first-difference
bound at half-widths `1e-12,1e-16,1e-20,1e-24,1e-28`.  Its largest and
smallest main-variation bounds are respectively
`9.267943404748744140625e+12` and
`0.004266662147212761137249348308841945254244`.  The `1e-24` rows
already have main-variation bounds above `42.66662147212760913816964603029191493988`.

This is a barrier for this direct first-difference box strategy, not an
impossibility theorem for all interval methods.  It shows that tiling the
physical domain with direct source-length boxes is not a credible next step.
The next implementation must recurse the truncated theta current, reduce its
active length, and propagate boxes through the transformed main and joined
endpoint currents before attempting quadrature.

Machine-audited companion:

```text
outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_quadratic_theta_mordell_parameter_box_phase_scale_gate.md
work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_quadratic_theta_mordell_parameter_box_phase_scale_gate.json
work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_quadratic_theta_mordell_parameter_box_phase_scale_gate.py
work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_quadratic_theta_mordell_parameter_box_phase_scale_gate.py
```

Primary-source boundary: the rotated integral and exact unit recurrence are
from Alexey Kuznetsov, *Computing the truncated theta function via Mordell
integral*, arXiv:1306.4081v2.  The box substitution, wall partition, and
first-difference phase audit are derived here.  Practical Gauss--Laguerre
quadrature is not used.

Pi provenance: every `pi` in this gate comes from the inherited Fourier phase
or rotated Mordell representation.  No fitted circle or polygon constant is
introduced.

This gate proves three central/transported Mordell parameter boxes, two
physical joined endpoint boxes, and two microscopic complete-current boxes.
It does not prove a scalable parameter-box current evaluator, the small-`tau`
branch, recursive interval control, physical quadrature, the non-A bound,
joined `R_after_A`, `R_Dir`, `Q_K-T`, an all-height theorem, `Lambda<=0`,
PF-infinity, RH, or a prize-level conclusion.
