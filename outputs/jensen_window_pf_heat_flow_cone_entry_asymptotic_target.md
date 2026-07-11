# Jensen-Window PF Heat-Flow Cone-Entry Asymptotic Target

Date: 2026-07-06

Status: open forward-flow theorem target. This is not a proof of propagated
monotone contractions, Jensen-window PF-infinity, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_heat_flow_cone_entry_asymptotic_target`.

Proof boundary: full ratio-cone entry is now proved at `lambda=-100`. This
artifact states the remaining infinite or collared finite forward-flow
legitimacy theorem. It does not prove `jwpf_06` or establish `Lambda <= 0`.

Machine-readable target:

```text
work/rh_compute/results/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.json
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_heat_flow_cone_entry_asymptotic_target.py
```

Current result:

```text
validated Jensen-window PF heat-flow cone-entry asymptotic target: 8 rows, 0 issues, 1 ready-to-apply rows, 0 live routes
```

## 2026-07-10 Upper-Wall Handoff

The later interval theorem certificate
`outputs/jensen_window_pf_kernel_mellin_upper_wall_certificate.md` proves
`x_k(lambda)<=1` for every real `lambda` and every `k>=1`. Together with the
separate Cauchy-Schwarz boundary-threshold lemma, both pointwise cone walls
are now exact. This target remains open only because full entry also requires
the adjacent-`k` inequality `x_(k+1)>=x_k`; the fixed-`k=22` real-`T`
certificate does not provide that all-`k` bridge.

## 2026-07-10 T=1156 Counterexample Handoff

The interval certificate
`outputs/jensen_window_pf_negative_lambda_t1156_monotone_wall_counterexample_certificate.md`
proves that the actual zeta sequence at `lambda=-1156` has
`x_121<x_120`. Thus the fixed-`k=22`, `T>=1156` certificate cannot be
promoted into full cone entry at `T=1156`, and a theorem asserting the
adjacent wall for all real heat parameters is false.

The parameter-specific route succeeds at `lambda=-100`: a repaired finite
prefix and zeta-specific saddle/tail theorem splice into a full all-`k` entry
time. Applying cone invariance forward still requires a rigorous infinite or
collared finite-flow argument.

The exact kernel summand-shift lemma
`outputs/jensen_window_pf_negative_lambda_kernel_summand_shift_lemma.md`
now sharpens the `lambda=-100`, `k>=300` tail. It rewrites every kernel
summand as a shifted copy of `phi_1` and proves that the complete shifted
`n=2..20` contribution on `v<=3/2` is below `2.122e-29` of the dominant
moment; `n>=21` starts beyond that split. The later dominance and paired
remainder theorems discharge the remaining tail burden.

## 2026-07-10 First-Summand Reduction

The later all-k certificate
`outputs/jensen_window_pf_negative_lambda_first_summand_dominance_certificate.md`
discharges that perturbation bound in stronger form:

```text
M_k=M_k^(1)*(1+delta_k), 0<=delta_k<=2/k^6, k>=300,
|L_k-L_k^(1)|<=16/(k-1)^6,              k>=301.
```

The finite Arb certificate
`outputs/jensen_window_pf_negative_lambda_m100_k320_collar_extension_certificate.md`
also promotes the repaired source through adjacent index `k=318`. The paired
compact and ray certificates prove

```text
L_k^(1)>=1/(4*k^2), k>=319,
```

for the first Newman summand at `lambda=-100`. Exact saddle geometry and nine
60-digit samples through `k=20000` are recorded in
`outputs/jensen_window_pf_negative_lambda_first_summand_saddle_wall_target.md`,
and the exact bridge
`outputs/jensen_window_pf_negative_lambda_first_summand_cumulant_bridge.md`
reduces this wall to the sufficient continuous estimate

```text
kappa_3,t(2*log(U))>=-37/(50*t^2), t>=318.
```

The Gamma contribution and rational transfer are exact. The interval theorem
`outputs/jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate.md`
certifies the expansion through fifth order. The paired theorem
`outputs/jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate.md`
then proves the remaining floor on the full compact mode range
`0.9264<=u<=5`. The analytic ray theorem
`outputs/jensen_window_pf_negative_lambda_first_summand_paired_remainder_ray_certificate.md`
proves

```text
seventh-order normalized saddle remainder>=-3/250, u>=5.
```

Thus the global cumulant hypothesis and lambda=-100 adjacent wall are closed;
the remaining issue in this target is downstream cone/all-order propagation.

## Entry Target

The ratio cone from the conditional invariance lemma is:

```text
x_k=(A_{k+1}/A_k)/(A_k/A_{k-1})
(2*k-1)/(2*k+1) <= x_k <= 1
x_{k+1} >= x_k
```

The exact conditional lemma is:

```text
outputs/jensen_window_pf_heat_flow_ratio_cone_invariance_lemma.md
validated Jensen-window PF heat-flow ratio cone invariance lemma: 6 exact rows, 315 lower rows, 315 upper rows, 310 monotone rows, 0 issues
```

The actual zeta coefficients now enter the full cone at `lambda=-100`:

```text
outputs/jensen_window_pf_negative_lambda_m100_full_cone_entry_certificate.md
```

To propagate that entry forward, one still needs:

```text
1. analytic legitimacy for the infinite or collared finite heat-flow argument;
2. explicit control preventing a first exit from escaping to k=infinity;
3. no endpoint PF, Jensen hyperbolicity, Laguerre-Polya membership, RH, or
   Lambda <= 0 as input.
```

## Formal Asymptotic Route

For `lambda=-T`, `T -> +infty`, write:

```text
Phi(u)=c0+c2*u^2+c4*u^4+c6*u^6+...
s=k+1/2
a=c2/c0
b=c4/c0
c=c6/c0
```

The fixed-k formal expansion has shape:

```text
A_k(-T)=C*4^(-k)*T^(-k-1/2)
        *(1+a*s/T+b*s*(s+1)/T^2+c*s*(s+1)*(s+2)/T^3+...)
```

So:

```text
log x_k
  = (2*b-a^2)/T^2
    + 2*(a^3*k-3*a*b*k-a*b+3*c*k+3*c)/T^3
    + O(T^-4)

log(x_{k+1}/x_k)
  = 2*(a^3-3*a*b+3*c)/T^3
    + O(T^-4)
```

The leading Gaussian term gives `x_k -> 1`, so it lands on the upper wall. A
strict cone-entry proof must use subleading Phi Taylor data and error control.

The live asymptotic route needed the local signs:

```text
2*b-a^2 < 0
2*(a^3-3*a*b+3*c) > 0
```

These local signs are now certified with explicit Phi-tail bounds:

```text
outputs/jensen_window_pf_phi_taylor_cone_entry_sign_scout.md
validated Jensen-window PF Phi Taylor cone-entry sign scout: 4 coefficient balls, 2 certified signs, 0 ready-to-apply rows, 0 issues
```

These local Taylor signs are now diagnostic history: the repaired-prefix plus
analytic-tail route proves entry directly at `lambda=-100`.

## Target Rows

```text
hfcet_01_exact_cone_entry_statement:
  ready theorem; full infinite ratio-cone entry at lambda=-100

hfcet_02_fixed_k_asymptotic_formal_reduction:
  formal fixed-k asymptotic route; insufficient for the full cone

hfcet_03_gaussian_leading_degeneracy:
  leading asymptotic gives x_k -> 1 and is too weak alone

hfcet_04_phi_taylor_sign_route:
  diagnostic fixed-k route superseded by the repaired-prefix/analytic-tail proof

hfcet_05_uniform_k_or_tail_collar_gap:
  hard gap: justify infinite-dimensional or collared finite forward flow
  from the now-proved lambda=-100 entry

hfcet_06_finite_grid_support:
  finite evidence only from existing Arb cache and negative-lambda prefix
  certificate

hfcet_07_forbidden_endpoint_shortcuts:
  endpoint PF/LP/RH/Newman-direction inputs remain circular

hfcet_08_conditional_application:
  cone entry plus conditional invariance would prove the covered
  monotone-contraction theorem
```

## Evidence And Boundaries

Finite sanity support:

```text
outputs/jensen_window_pf_monotone_contraction_stress.md
validated Jensen-window PF monotone contraction stress: 2875 rows, 2875 positive rows, 0 issues

outputs/jensen_window_pf_negative_lambda_cone_entry_prefix_scout.md
validated Jensen-window PF negative-lambda cone-entry prefix scout: 69 coefficient rows, 63 lower-wall rows, 63 upper-wall rows, 60 monotone-gap rows, 0 issues

outputs/jensen_window_pf_negative_lambda_finite_collar_contract.md
validated Jensen-window PF negative-lambda finite-collar contract: active depth K=19, 57 active lower rows, 57 active upper rows, 57 active monotone rows, 6 collar lower/upper rows, 3 collar monotone rows, 0 issues
```

The negative-lambda prefix scout certifies the actual ACB-enclosed ratios for
`lambda=-25,-50,-100`, lower/upper walls through `k<=21`, and monotone gaps
through `k<=20`. It remains finite evidence only; it does not supply the all-`k`
tail theorem or finite-collar flow theorem.

The finite-collar contract extracts the usable active depth from that prefix:
`K=19`, with first collar `x_20` and second collar `x_21`. It records that
pushing to `K=20` needs further coefficient enclosures or a genuine tail theorem.

The analytic handoff is now stated as:

```text
outputs/jensen_window_pf_negative_lambda_defect_tail_theorem_target.md
validated Jensen-window PF negative-lambda defect-tail theorem target: 8 rows, 0 issues, 0 ready-to-apply rows, 2 live routes
```

It requires `0 <= d_k <= 2/(2*k+1)` for `k>=22` and `d_(k+1)<=d_k`
for `k>=21`, where `d_k=1-x_k`.

Boundary and invariance algebra:

```text
outputs/jensen_window_pf_heat_flow_boundary_threshold_lemma.md
validated Jensen-window PF heat-flow boundary threshold lemma: 5 exact rows, 315 strong-threshold rows, 315 heat-threshold rows, 0 issues

outputs/jensen_window_pf_heat_flow_ratio_cone_invariance_lemma.md
validated Jensen-window PF heat-flow ratio cone invariance lemma: 6 exact rows, 315 lower rows, 315 upper rows, 310 monotone rows, 0 issues
```

The target remains open because none of these artifacts proves the zeta
coefficient sequence enters the full infinite cone at any suitable lambda.
