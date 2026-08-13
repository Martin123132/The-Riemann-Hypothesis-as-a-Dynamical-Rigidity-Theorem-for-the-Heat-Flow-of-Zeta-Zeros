# Inner/outer corrected-defect projection

Date: 2026-08-13

Status: rigorous branch-separated corrected-defect enclosure on the ordinary
top corridor; this is not yet a bound for the paired companion

Split the carrier-suppressed corrected statistic at event zero:

```text
W_corr=W_inner+W_outer,
W_inner=sum_(m=39696)^39894 a_m(G_ex,m-G_opp,m),
W_outer=sum_(m=39895)^40094 a_m(G_ex,m-G_opp,m),       (IO1)
a_39894=0.
```

The beta-minus-four parts were evaluated as separate 198- and 200-mode
Fourier polynomials on `65536` midpoint panels at five height nodes.  No
modewise integral or post-summation triangle was used.  Separate derivative
majorants give

```text
inner y error  <[2.454609350472676196041354767812238484839399626274357672898069852439335e-6 +/- 3.12e-75],
inner t error  <[2.121363330522976683876470243108722683582376598644426669961748914864063e-6 +/- 2.56e-75],
outer y error  <[2.498098263254372095657565935104497394855112389617512986482162129603088e-6 +/- 2.84e-75],
outer t error  <[2.142780599949190184685105944798293601737018747265965574099400180709426e-6 +/- 2.82e-75].                 (IO2)
```

The exact Kummer-ODE profile correction is split with the two certified
whole-kernel norms, giving bounds

```text
inner exact finite-t correction <[4.846799008581719365689356468453006519279549878379020540582005561970243e-12 +/- 3.08e-83],
outer exact finite-t correction <[7.149933274481865640424026552806462455707074888775419215574443400097031e-12 +/- 3.18e-82]. (IO3)
```

After restoring the common carrier and paired projection, the rigorous
component intervals are

```text
[-0.0001203897766118226554484106714916196648123880167197652086037612573236531 +/- 2.04e-75] < R_inner
 < [-6.843790506728551139731184939265476041102790808132825802653609266545410e-5 +/- 3.40e-75],
[-0.0001106889260878981874183396287309648315629755044380029225474474189586277 +/- 3.22e-74] < R_outer
 < [-5.799699853739199647304659546700487497225383934159084706115962090689216e-5 +/- 3.38e-75].      (IO4)
```

Their sum is the previously certified `R_corr`.  The outer component cancels
exactly in Section 11.406, so the operative fold increment is

```text
Delta Q_fold=R_inner+R_pair,
R_pair=2Re[e^(-i*pi/8)J_pair].                        (IO5)
```

In fact (IO4) certifies `R_inner<0` uniformly.  Consequently the explicit
sufficient target for a negative fold increment is

```text
R_pair<[6.843790506728551139731184939265476041102790808132825802653609266545409e-5 +/- 3.45e-75].                                      (IO6)
```

No sign for the fold increment is claimed until (IO6) is proved on the same
corridor.

Pi provenance: the projection factor and corridor width come from
`beta^3=pi*C^2/8`, `hY=2pi`, inverse Airy Fourier normalization, and exact
odd-square half-domain reflection.  No fitted constant is introduced.

Machine-audited companion:

```text
outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_inner_outer_corrected_defect_projection_gate.md
work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_inner_outer_corrected_defect_projection_gate.json
work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_inner_outer_corrected_defect_projection_gate.py
work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_inner_outer_corrected_defect_projection_gate.py
```

No paired-companion bound, signed fold increment, fixed-state residual,
all-corridor telescope, complete `Q_K-T` or `T_upper`, height-uniform theorem,
`Lambda<=0`, PF-infinity, RH, or prize-level conclusion is proved.
