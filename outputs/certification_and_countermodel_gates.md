# Certification And Countermodel Gates

Date: 2026-07-10

Status: proof-safety gate ledger. This is not a proof of RH or `Lambda <= 0`; it defines certification levels, manifest checks, and countermodel boundaries.

Purpose:

Define hard gates for the RH dynamical-rigidity proof programme so numerical evidence, conditional stability, and local zero-flow intuition do not get mistaken for a proof of `Lambda <= 0`.

## Gate A: Numerical Evidence Levels

### Level 0: Notebook Evidence

Examples:

```text
floating roots appear real
plots look stable
determinants print with expected signs
```

Allowed use:

```text
exploration only
```

Not allowed:

```text
manuscript theorem support
```

### Level 1: Reproducible High-Precision Evidence

Examples from current work:

```text
repro_hankel_15c: 585/585 signed determinants positive
repro_jensen_fragile_d16_d20_n50_n60: 55/55 Jensen cases passed
```

Allowed use:

```text
finite evidence, stress-test result, route selection
```

Required metadata:

```text
dps
n_sum
cutoff
lambda grid
degree/shift ranges
JSONL or CSV row logs
failure rows preserved
```

Not allowed:

```text
claiming exact finite certification
claiming all-degree or all-shift behavior
claiming RH evidence beyond diagnostic status
```

### Level 2: Independent Stability Evidence

Examples from current work:

```text
high-precision rerun of fragile Hankel block:
lambda=0, m=9..12, shift=0..8, dps=160, 36/36 passed
```

Allowed use:

```text
robust finite evidence
```

Still not allowed:

```text
exact certification of transcendental moment determinants
```

### Level 3: Ball/Interval Propagation Evidence

Current scaffold:

```text
work/rh_compute/scripts/arb_hankel_sign_probe.py
work/rh_compute/scripts/arb_toeplitz_pf_probe.py
```

These probe whether determinant signs survive a chosen coefficient error radius using `python-flint` / Arb balls. The Toeplitz/PF probe also has a structural-zero front gate for identically zero upper-triangular Toeplitz minors.

Allowed use:

```text
conditional certificate:
if cached coefficient balls enclose exact A_k(lambda), then determinant sign is certified
```

Not allowed unless moment enclosures are proved:

```text
unconditional determinant certificate
```

### Level 4: Full Finite Certificate

Required:

```text
rigorous enclosure of each A_k(lambda)
validated truncation/tail bound for Phi series
validated integration/quadrature error bound
ball determinant sign separated from zero
for Jensen: root count stable under coefficient balls, not just rationalized decimal coefficients
```

Only Level 4 can support wording like:

```text
certified finite computation
```

Update added 2026-07-03:

The truncation/cutoff tail part has a first analytic bound scaffold:

```text
work/rh_compute/scripts/coefficient_tail_bounds.py
```

For the `repro_hankel_15c` settings:

```text
n_sum = 100
cutoff = 8
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}
k <= 32
```

the worst omitted `A_k` tail is about `1.91e-13597`.

For the fragile Jensen cache settings:

```text
n_sum = 120
cutoff = 9
lambda = 0
k <= 80
```

the worst omitted `A_k` tail is about `1.53e-19597`.

These are far below the `1e-80` Arb propagation radius. The remaining Level-4 blocker is validated quadrature for the retained finite sum on `[0, cutoff]`.

Further update added 2026-07-03:

Selected retained-integral coefficient enclosures now exist:

```text
work/rh_compute/scripts/acb_coefficient_enclosures.py
```

using `python-flint`'s rigorous `acb.integral`.

Current certified coefficient ranges:

```text
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}, k = 0..32
n_sum = 100, cutoff = 8
```

Each certified `A_k` ball is contained in the cache ball `A_k(cache) +/- 1e-80`. This promotes the corresponding structural Arb Toeplitz/PF probes to finite certificates in the ranges listed in `outputs/finite_toeplitz_certificate_status.md`.

Signed-Hankel update added 2026-07-04:

The signed-Hankel grid:

```text
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}
m = 0..19
s = 0..24
```

is now checked directly against the rigorous `A_ball` enclosure rows:

```text
work/rh_compute/scripts/arb_hankel_enclosure_sign_probe.py
python work/rh_compute/scripts/check_hankel_certificate_manifest.py
```

The checker validates:

```text
2,500 signed-Hankel finite certificates
```

This is a finite determinant certificate only; it does not cross Gate D without an all-order bridge theorem.

Expanded certified Toeplitz/PF ranges now include:

```text
lambda grid, N = 10, orders <= 4
lambda grid, N = 12, orders <= 4
lambda grid, N = 14, orders <= 4
lambda grid, N = 16, orders <= 4
lambda grid, N = 18, orders <= 4
lambda grid, N = 10, orders <= 5
lambda grid, N = 12, orders <= 5
lambda grid, N = 14, orders <= 5
lambda grid, N = 16, orders <= 5
lambda grid, N = 18, orders <= 5
lambda grid, N = 20, orders <= 4
lambda grid, N = 20, orders <= 5
lambda grid, N = 22, orders <= 4
lambda grid, N = 22, orders <= 5
lambda grid, N = 24, orders <= 4
lambda grid, N = 26, orders <= 4
lambda grid, N = 28, orders <= 4
lambda grid, N = 30, orders <= 4
lambda grid, N = 24, orders <= 5
```

GPU/Torch prefilter note added 2026-07-04:

```text
work/rh_compute/scripts/torch_toeplitz_prefilter.py
```

This is only a non-rigorous floating-point screen. It can prioritize Arb work but cannot promote any claim to "certified finite computation." The proof boundary is recorded in:

```text
outputs/gpu_prefilter_and_certificate_boundary.md
```

Edrei-log diagnostic added 2026-07-04:

```text
work/rh_compute/scripts/arb_edrei_log_sign_probe.py
work/rh_compute/scripts/check_edrei_log_sign_manifest.py
work/rh_compute/scripts/arb_edrei_power_hankel_probe.py
work/rh_compute/scripts/check_edrei_power_hankel_manifest.py
work/rh_compute/scripts/edrei_power_hankel_midpoint_scout.py
work/rh_compute/scripts/check_edrei_power_hankel_frontier_manifest.py
work/rh_compute/scripts/check_edrei_power_hankel_boundary_manifest.py
work/rh_compute/scripts/edrei_moment_quadrature_scout.py
work/rh_compute/scripts/check_edrei_quadrature_scout_manifest.py
outputs/edrei_log_sign_diagnostic.md
outputs/edrei_moment_quadrature_scout.md
```

Using rigorous `c_ball` coefficient enclosures, the manifest validates:

```text
320 finite Edrei-log sign diagnostics
4,205 finite Edrei power-Hankel diagnostics
```

These are finite necessary-condition diagnostics for the entire coefficient PF route. They do not prove PF-infinity or Laguerre-Polya membership. The promoted power-Hankel set now includes the uniform five-lambda frontier through the finite staircase `2m+s <= 57`, `s >= 1`: from `m = 0, s = 1..57` down to `m = 28, s = 1`. The `n = 51` layer needed a tighter `lambda = 1e-6`, `k = 0..51`, `dps = 220`, `abs_tol = 1e-140` coefficient repair; the `n = 53` and later layers needed the same tighter coefficient treatment through `k = 57` for all five lambdas, followed by `dps = 2400` Edrei-log reruns. The boundary checker validates both the retired inconclusive rows and the repaired positive rows. The midpoint-only scout at `m = 20, s = 1` was mixed before tightening and is explicitly non-rigorous; only Arb-separated determinants are promoted. The midpoint frontier checker validates that warning artifact separately so future notes do not confuse center-only signs with certified signs.

The Edrei moment-recurrence scout is a finite Arb diagnostic. It validates
one positive recurrence range, orders `2..12`, and one broader frontier range,
orders `2..20`, where the persisted Edrei-log enclosures become inconclusive
from order `13` onward. This is a theorem-search frontier note, not a proof
of the all-order Edrei representation.

Promoted finite certificates and non-promotion guards now have manifest gates:

```text
python work/rh_compute/scripts/check_core_proof_programme_gates.py
python work/rh_compute/scripts/check_output_reference_integrity.py
python work/rh_compute/scripts/check_output_status_manifest.py
python work/rh_compute/scripts/check_proof_claim_ledger.py
python work/rh_compute/scripts/check_signed_hankel_jensen_dependency_graph.py
python work/rh_compute/scripts/check_result_language_boundaries.py
python work/rh_compute/scripts/check_toeplitz_certificate_manifest.py
python work/rh_compute/scripts/check_toeplitz_jacobi_trudi_map.py
python work/rh_compute/scripts/check_hankel_sign_consistency_reduction_audit.py
python work/rh_compute/scripts/check_arb_hankel_sign_consistency_reduction_manifest.py
python work/rh_compute/scripts/check_arb_shifted_hankel_sign_consistency_manifest.py
python work/rh_compute/scripts/check_jensen_hankel_bridge_algebra.py
python work/rh_compute/scripts/check_signed_hankel_jensen_bridge_target.py
python work/rh_compute/scripts/check_jensen_window_pf_bridge_obligations.py
python work/rh_compute/scripts/check_jensen_window_pf_theorem_machinery_fit_matrix.py
python work/rh_compute/scripts/check_jensen_window_pf_structural_ansatz_matrix.py
python work/rh_compute/scripts/check_jensen_window_pf_cauchy_binet_low_degree_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_log_concavity_frontier_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_ratio_condition_scout.py
python work/rh_compute/scripts/check_edrei_log_sign_manifest.py
python work/rh_compute/scripts/check_edrei_power_hankel_manifest.py
python work/rh_compute/scripts/check_edrei_power_hankel_frontier_manifest.py
python work/rh_compute/scripts/check_edrei_power_hankel_boundary_manifest.py
python work/rh_compute/scripts/check_edrei_quadrature_scout_manifest.py
```

The checker validates the advertised positive summaries, empty problem-row files,
available stderr logs, and the beta negative control before a finite range is
treated as part of the proof programme.

The Grussler-Damm reshaped-Hankel sign-consistency route now has both an
exact-rationalized cache point audit and an Arb/enclosure-backed finite
certificate for the five-lambda grid, `k = 2..7`, `N = 20`. The Arb manifest
validates `689,795/689,795` finite minors as positive and separated from zero.
This is a finite certificate only; it does not prove the all-order
sign-consistency bridge.

The shifted reshaped-Hankel Arb manifest adds finite shift coverage:
five-lambda grid, `n = 0..20`, `k = 2..5`, `N = 18`, validating
`1,322,685/1,322,685` finite shifted minors. This is still a finite rectangle, not
the all-shift theorem required by the bridge target.

The order-6 shifted Arb manifest adds the next-order finite slice:
five-lambda grid, `n = 0..20`, `k = 6`, `N = 16`, validating
`840,840/840,840` finite shifted minors using coefficient enclosures through
`A_40`. This is a finite next-order stress test, not an all-order theorem.

The order-7 shifted Arb manifest adds another bounded next-order slice:
five-lambda grid, `n = 0..20`, `k = 7`, `N = 15`, validating
`675,675/675,675` finite shifted minors using coefficient enclosures through
`A_40`. This is finite evidence only.

The order-8 shifted Arb manifest adds the next bounded slice:
five-lambda grid, `n = 0..20`, `k = 8`, `N = 14`, validating
`315,315/315,315` finite shifted minors using coefficient enclosures through
`A_40`. This is finite evidence only.

The consolidated shifted staircase checker is:

```text
python work/rh_compute/scripts/check_arb_shifted_hankel_staircase_manifest.py
```

It validates all promoted shifted slices together:

```text
3,154,515/3,154,515 finite shifted minors
k = 2..5 at N = 18
k = 6 at N = 16
k = 7 at N = 15
k = 8 at N = 14
```

This is a consolidated finite-certificate manifest, not an all-order theorem.

The Jensen/Hankel bridge algebra gate records the exact degree-2 identity
between the `m = 1` signed-Hankel determinant and the degree-2 Jensen
discriminant, and an exact positive rational degree-3 countermodel showing
that finite low-order reshaped-Hankel signs cannot be promoted into Jensen
hyperbolicity.

The signed-Hankel/Jensen bridge target checker requires the active theorem
target to remain all-order and all-shift, with the finite Arb certificate
listed only as evidence and the degree-3 obstruction listed as a kill gate.

The signed-Hankel/Jensen dependency graph documents the allowed relation
between finite evidence, countermodel gates, open theorem targets, and the
not-proved Newman-direction conclusion:

```text
outputs/signed_hankel_jensen_dependency_graph.md
work/rh_compute/results/signed_hankel_jensen_dependency_graph.json
python work/rh_compute/scripts/check_signed_hankel_jensen_dependency_graph.py
```

It is a dependency hygiene gate only, not a proof bridge.

The Jensen-window PF bridge obligation ledger decomposes the central open
bridge target into exact, finite, open, conditional, rejected, and
route-separated obligations:

```text
outputs/jensen_window_pf_bridge_obligations.md
work/rh_compute/results/jensen_window_pf_bridge_obligations.json
python work/rh_compute/scripts/check_jensen_window_pf_bridge_obligations.py
```

It validates `10` obligations with `3` still open. This is theorem-search
hygiene only, not proof of the bridge.

The Jensen-window PF theorem machinery fit matrix audits total-positivity,
PF-sequence, zero-preserver, sign-regularity, downstream Laguerre-Polya, and
rejected-shortcut theorem families against the central `jwpf_06` bridge:

```text
outputs/jensen_window_pf_theorem_machinery_fit_matrix.md
work/rh_compute/results/jensen_window_pf_theorem_machinery_fit_matrix.json
python work/rh_compute/scripts/check_jensen_window_pf_theorem_machinery_fit_matrix.py
```

It validates `7` rows with `0` ready-to-apply rows. This is source-anchored
theorem-search hygiene only.

The Jensen-window PF structural ansatz matrix turns the theorem-machinery
audit into exact hard tests for candidate proof mechanisms:

```text
outputs/jensen_window_pf_structural_ansatz_matrix.md
work/rh_compute/results/jensen_window_pf_structural_ansatz_matrix.json
python work/rh_compute/scripts/check_jensen_window_pf_structural_ansatz_matrix.py
```

It validates `6` candidate, blocked, and rejected ansatz rows with `0`
ready-to-apply rows. This is structural proof-search hygiene only.

The Jensen-window PF Cauchy-Binet low-degree scout checks the live kernel
ansatz against the selected exact degree-2/3/4 formulas:

```text
outputs/jensen_window_pf_cauchy_binet_low_degree_scout.md
work/rh_compute/results/jensen_window_pf_cauchy_binet_low_degree_scout.json
python work/rh_compute/scripts/check_jensen_window_pf_cauchy_binet_low_degree_scout.py
```

It validates `15` formula rows with nonnegative Bernstein coefficients under
adjacent log-concavity and records `0` kernel identities found.

The Jensen-window PF log-concavity frontier scout extends that diagnostic to
larger contiguous minors:

```text
outputs/jensen_window_pf_log_concavity_frontier_scout.md
work/rh_compute/results/jensen_window_pf_log_concavity_frontier_scout.json
python work/rh_compute/scripts/check_jensen_window_pf_log_concavity_frontier_scout.py
```

It validates `14` contiguous rows and locates first Bernstein-certificate
failures at degree 3 size 6 and degree 4 size 5, plus first exact countermodel
negatives at degree 3 size 8 and degree 4 size 6.

The Jensen-window PF ratio-condition scout rejects nearby strengthened
conditions as standalone bridges:

```text
outputs/jensen_window_pf_ratio_condition_scout.md
work/rh_compute/results/jensen_window_pf_ratio_condition_scout.json
python work/rh_compute/scripts/check_jensen_window_pf_ratio_condition_scout.py
```

It validates `7` candidate rows with `4` rejected by exact countermodel and
`1` rejected by constructed positive extension.

The kernel Mellin upper-wall theorem and its adjacent-wall countermodel form
a matched promotion gate:

```text
outputs/jensen_window_pf_kernel_mellin_upper_wall_certificate.md
outputs/jensen_window_pf_log_concave_mellin_monotone_wall_countermodel.md
python work/rh_compute/scripts/check_jensen_window_pf_kernel_mellin_upper_wall_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_log_concave_mellin_monotone_wall_countermodel.py
```

The first proves `x_k(lambda)<=1` for every real `lambda` and every `k>=1`
from a global interval-backed log-concavity theorem and Berwald-Borell. The
second gives an exact log-concave Mellin witness with `x_2<x_1`, preventing
that valid upper-wall theorem from being promoted to the still-open
adjacent-`k` wall.

The actual zeta kernel now also has a parameter-specific promotion gate:

```text
outputs/jensen_window_pf_negative_lambda_t1156_monotone_wall_counterexample_certificate.md
work/rh_compute/results/acb_enclosures_lambda_m1156_k119_k122_dps250.jsonl
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_t1156_monotone_wall_counterexample_certificate.py
```

Four rigorous coefficient enclosures at `lambda=-1156` certify
`x_121-x_120<-1.68e-8`. This leaves the fixed-`k=22` theorem intact but
blocks its promotion to an all-`k` cone theorem at `T=1156`. Cone entry is
therefore redirected to a moderate fixed negative parameter with a finite
collar plus a separate eventual-`k` theorem.

The first exact reduction on that revised route is:

```text
outputs/jensen_window_pf_negative_lambda_kernel_summand_shift_lemma.md
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_kernel_summand_shift_lemma.py
```

It proves every Newman-kernel summand is a shifted `n=1` copy and certifies
the complete compact `n=2..20` perturbation below `2.122e-29` for
`lambda=-100`, `k>=300`.

The compact reduction, forward-flow handoff, and first higher-minor audit are
now covered by twenty-two linked gates:

```text
outputs/jensen_window_pf_negative_lambda_first_summand_dominance_certificate.md
outputs/jensen_window_pf_negative_lambda_m100_k320_collar_extension_certificate.md
outputs/jensen_window_pf_negative_lambda_first_summand_saddle_wall_target.md
outputs/jensen_window_pf_negative_lambda_first_summand_cumulant_bridge.md
outputs/jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate.md
outputs/jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate.md
outputs/jensen_window_pf_negative_lambda_first_summand_paired_remainder_ray_certificate.md
outputs/jensen_window_pf_negative_lambda_m100_full_cone_entry_certificate.md
outputs/jensen_window_pf_heat_flow_infinite_cone_invariance_certificate.md
outputs/jensen_window_pf_defect_complete_monotonicity_scout.md
outputs/jensen_window_pf_heat_flow_jensen_hierarchy_lemma.md
outputs/jensen_window_pf_rank_two_boundary_family_lemma.md
outputs/jensen_window_pf_multiplier_counting_measure_target.md
outputs/jensen_window_pf_mellin_multiplier_power_sum_obstruction.md
outputs/jensen_window_pf_cubic_reciprocal_defect_invariance_lemma.md
outputs/jensen_window_pf_cubic_m100_tail_entry_certificate.md
outputs/jensen_window_pf_cubic_forward_uniform_tail_certificate.md
outputs/jensen_window_pf_quartic_boundary_flow_obstruction.md
outputs/jensen_window_pf_quartic_double_root_threshold_lemma.md
outputs/jensen_window_pf_quartic_quintic_polar_contact_lemma.md
outputs/jensen_window_pf_cofinal_degree_polar_closure_lemma.md
outputs/jensen_window_pf_cofinal_scaling_limit_equivalence_gate.md
outputs/jensen_window_pf_polar_heat_collision_cascade_lemma.md
outputs/jensen_window_pf_scaled_double_zero_boundary_layer_lemma.md
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_first_summand_dominance_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_m100_k320_collar_extension_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_first_summand_saddle_wall_target.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_first_summand_cumulant_bridge.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_first_summand_paired_remainder_ray_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_m100_full_cone_entry_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_heat_flow_infinite_cone_invariance_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_defect_complete_monotonicity_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_heat_flow_jensen_hierarchy_lemma.py
python work/rh_compute/scripts/check_jensen_window_pf_rank_two_boundary_family_lemma.py
python work/rh_compute/scripts/check_jensen_window_pf_multiplier_counting_measure_target.py
python work/rh_compute/scripts/check_jensen_window_pf_mellin_multiplier_power_sum_obstruction.py
python work/rh_compute/scripts/check_jensen_window_pf_cubic_reciprocal_defect_invariance_lemma.py
python work/rh_compute/scripts/check_jensen_window_pf_cubic_m100_tail_entry_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_cubic_forward_uniform_tail_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_quartic_boundary_flow_obstruction.py
python work/rh_compute/scripts/check_jensen_window_pf_quartic_double_root_threshold_lemma.py
python work/rh_compute/scripts/check_jensen_window_pf_quartic_quintic_polar_contact_lemma.py
python work/rh_compute/scripts/check_jensen_window_pf_cofinal_degree_polar_closure_lemma.py
python work/rh_compute/scripts/check_jensen_window_pf_cofinal_scaling_limit_equivalence_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_polar_heat_collision_cascade_lemma.py
python work/rh_compute/scripts/check_jensen_window_pf_scaled_double_zero_boundary_layer_lemma.py
```

The first discharges the shifted far tail and proves
`M_k=M_k^(1)(1+delta_k)`, `0<=delta_k<=2/k^6`, for every integer `k>=300`.
The second extends the certified `lambda=-100` adjacent wall through `k=318`.
The third records exact saddle geometry and, using the paired remainder
theorems, proves `L_k^(1)>=1/(4*k^2)` for every integer `k>=319`.
The fourth proves the exact cumulant/Gamma reduction to the sufficient estimate
`kappa_3,t(2*log(U))>=-37/(50*t^2)` for every real `t>=318`.
The fifth certifies the leading saddle term below `13/20`, its cubic
Edgeworth correction below `1/100`, and its fifth-order correction below
`1/1000`. The sixth proves the seventh-order normalized remainder floor
`-79/1000` throughout `0.9264<=u<=5`, reaching `t>1.5241916613e10`.
The seventh proves `H_t>=-3/250` on `u>=5`, closing the global cumulant
hypothesis and lambda=-100 adjacent wall.
The eighth merges the repaired prefix and composes it with the exact pointwise
walls and analytic adjacent tail, proving full ratio-cone entry at lambda=-100.
The ninth proves an infinite-sequence first-crossing maximum principle for the
heat ODE and propagates that full cone to every finite lambda>=-100, including
lambda=0. The all-order PF/Jensen bridge remains open.
The tenth certifies both defect and negative-log-contraction alternating
hierarchies through order 8 on the cached grid, retains 838 high-order
inconclusives, and supplies the exact degree-3 discriminant countermodel. The
eleventh records the exact coupled
degree-shift Jensen heat hierarchy and the one-atom Hausdorff full-cone
boundary witness with normalized frontier derivative `329728/2109375>0`,
showing that complete defect monotonicity plus the local heat ODE still does
not give cubic invariance.
The twelfth proves an all-degree rank-two boundary factorization and supplies
integer pointwise-product multiplier closure. The exact `-937/3456` mixture
and `<-27/125` fractional-power cubic guards reject arbitrary positive measure
weights in the resulting canonical-product ansatz.
The thirteenth records the open sufficient unit-atomic factorization target,
with exact convergence, Laplace-kernel, and failure-gate requirements.
The fourteenth Arb-encloses the natural Mellin interpolation's log moments and
finds three negative power-sum Hankel determinants. It closes that continuous-
index strengthening but explicitly leaves the integer-only product target open.
The fifteenth factors the cubic discriminant in reciprocal defect coordinates,
proves the shift-coupled inward boundary condition, and certifies 318
lambda=-100 prefix margins plus 310 nonnegative-grid margins.
The sixteenth certifies negative third log-cumulant on 4,074 compact Arb blocks
and an explicit analytic ray, derives a reciprocal-defect floor, and closes
`q_(k+1)-q_k<1` for every `k>=319`. Together with the prefix this proves all
shifted cubics hyperbolic at lambda=-100.
The seventeenth derives the exact one-sided q-increment flow, proves the
weighted source cap `sqrt(k)*g_k'<=7*r_k`, and uses an explicit coercive
supersolution to establish compact-forward `O_L(k^(-1/2))` tail decay. This
closes the cubic first-crossing handoff and proves every shifted cubic
hyperbolic at lambda=0. Degree four and the all-degree bridge remain open.
The eighteenth derives the normalized quartic discriminant and supplies an
exact rational hyperbolic boundary point with four positive ratio margins and
three strict cubic margins but `Q'/r_1<0`. It blocks promotion of the cubic
cone unchanged to degree four without claiming any failure of the actual zeta
trajectory.
The nineteenth converts that obstruction into the exact branch-aware inward
law `(3*a^2-4*a+p)*(u-U)<=0`, proves shift-independence of `U`, and identifies
the necessary triple-root equality and tangent factor. A globally closed
quartic contraction cone is still open.
The twentieth proves the adjacent-degree polar identity, shows that a quartic
double root under a hyperbolic quintic extension forces a quintic triple root,
and identifies `u=U` as exactly that contact value. The identity extends to all
degrees, but its infinite top-down hierarchy remains open.
The twenty-first iterates polar descent and proves that an unbounded sequence
of hyperbolic terminal degrees at each shift suffices for every finite degree.
The finite Sturm data now reach every degree from 3 through 12 on the bounded
grid, totaling 1050 certified rows. This is still not a cofinal sequence, so
the terminal theorem remains open.
The twenty-second proves that the fixed-shift cofinal degree condition is
equivalent, through the Jensen scaling limit, to Laguerre-Polya membership.
It blocks using cofinal degrees as a weaker assumed endpoint or promoting the
bounded degree-3 through degree-12 ladder into that limit.
The twenty-third proves the general polar heat-collision cascade: a fixed-degree
multiple root under a hyperbolic upper tower forces linearly increasing
multiplicity and hence an exponential-polynomial entire function. A genuine
Xi-type LP boundary must therefore escape to unbounded Jensen degree.
The twenty-fourth resolves the first two orders of that escaping layer. It
proves the universal quadratic, the D^(-3/2) root gap, the rho/(8D) collision
shift, and the D^(-2) regularized external-field correction. The required
degree- and zero-height-uniform remainder estimate remains open.
The twenty-fifth identifies that external correction exactly with the
regularized Newman root field, derives the positive pair stiffness and the
second-order squared-gap law, and supplies opposite-sign finite LP
countermodels. The remaining sign and tail estimate is therefore Xi-specific.
The twenty-sixth rewrites each unshifted Jensen polynomial as a positive
integral of generalized Laguerre kernels. Fixed scales and every half-integer
Gamma law are hyperbolic in all degrees, but exact positive-mixture and
log-concave-density countermodels show that integration alone is not a
preserver. A common-interlacing or total-positive Xi connection remains open.
The twenty-seventh computes the classical continuum root field as
`-pi/8+O(1/x)`, matching the positive-time quantile drift, and imports the
published fixed-time high-zero simplicity theorem to confine any multiple zero
below `exp(C/t)`. Exact bounded-location even-lattice deformations force either
unbounded field sign, so only a lambda-uniform reciprocal-gap estimate or an
analytic compact-region exclusion can close this reduction.
The twenty-eighth integrates the uniform zero-count discrepancy against the
field kernel and proves that every mesoscopic and far contribution is
`O(1/log c)` outside `H=log(4c)^2`. The remaining field is the inverse-square
weighted odd local count. An exact even backward-heat polynomial has field
`-pi/8` and drift `-pi/4` while still undergoing positive pair birth, proving
that field balance alone is not a collision obstruction.

The core runner is documented in:

```text
outputs/core_proof_programme_gates.md
```

Current core runner status:

```text
validated 193/193 core proof-programme gates in the integrated non-slow run.
The current runner has 200 gates: 193 non-slow and 7 slow.
```

Current manifest status:

```text
95 promoted positive Toeplitz/PF certificate summaries
1 beta negative control
```

Countermodel library added 2026-07-04:

```text
outputs/countermodel_library.md
work/rh_compute/scripts/countermodel_gate_examples.py
```

The executable gate validates:

```text
local heat-equation square-root birth with g' = 4/g
exact rational finite moment-recurrence prefix trap:
  recurrence order <= 12
  moments 0..23 preserved
  positive edited moment 24 breaks the next Hankel/moment gate
exact rational Stieltjes multiplier trap:
  positive measure 10 delta_0 + delta_1 + delta_2
  leading Hankel determinants 12, 51, 40, 0
  c_k = mu_k/(2k)! gives c_1^2 - c_0 c_2 = -1/4
exact finite shifted-principal signed-Hankel grid trap:
  base a_k = 1/k!
  m <= 4, shifts <= 8 preserved
  positive a_17 breaks the next shifted m = 1 signed-Hankel/Jensen gate
exact finite Schur/Toeplitz prefix trap:
  base h_k = 1/k!, generating function exp(z)
  h_0..h_6 and 2,940 finite Toeplitz/Schur tests preserved
  positive h_7 breaks the Jacobi-Trudi skew-Schur determinant s_(6,6)
finite Jensen-window rectangle trap:
  current Jensen-window PF/Sturm coefficient inputs A_0..A_25 preserved
  positive A_26 breaks degree-2, shift-24 Jensen hyperbolicity
positive one-term extensions of the current coefficient prefixes that break:
  order-2 Toeplitz/PF
  m = 1 signed-Hankel
  degree-2 Jensen hyperbolicity
```

These examples are not claims about the actual zeta coefficients. They are proof-safety tests: any proposed bridge lemma that uses only local repulsion, only a finite coefficient prefix, only a finite Schur/Toeplitz shape prefix, only a finite shifted-principal Hankel grid, only a finite moment/recurrence prefix, or only Stieltjes/Hankel positivity where Toeplitz PF is required cannot prove `Lambda <= 0`.

## Gate B: Circularity Checks

No argument may assume:

```text
all zeros are real at lambda = 0
```

when the conclusion is:

```text
Lambda <= 0
```

Reason:

RH is equivalent to `Lambda <= 0` in the de Bruijn-Newman framework.

Allowed conditional statement:

```text
If all zeros are real at lambda_0 and the comparison hypotheses hold,
then real-rootedness persists for larger lambda.
```

Forbidden promotion:

```text
This proves RH when lambda_0 = 0.
```

## Gate C: Local-Repulsion Countermodel

Any proof step relying only on local repulsion must be tested against:

```text
F(lambda,t) = (t-a)^2 - 2(lambda-lambda_*).
```

This satisfies:

```text
F_lambda = -F_tt
```

and has:

```text
lambda < lambda_*: complex conjugate pair
lambda = lambda_*: double real zero
lambda > lambda_*: two real zeros
```

On the real-zero side:

```text
g' = 4/g.
```

Therefore:

```text
local repulsion is compatible with positive square-root birth.
```

Forbidden claim:

```text
The local law g' ~ 4/g excludes positive Newman boundary.
```

Allowed claim:

```text
The local law describes the already-real side after birth.
```

The executable check is:

```text
python work/rh_compute/scripts/countermodel_gate_examples.py
```

## Gate D: Signed-Hankel Route Claims

Current evidence supports:

```text
finite signed Hankel sign pattern for tested m,s,lambda
finite Jensen hyperbolicity for tested d,n,lambda
```

Current evidence does not support:

```text
signed Hankel signs imply all Jensen polynomials are hyperbolic
signed Hankel signs imply H_lambda is Laguerre-Polya
signed Hankel signs imply Lambda <= 0
finite shifted-principal Hankel grids imply all-order sign consistency
```

To cross this gate, we need a theorem of the form:

```text
signed/indefinite Hankel sign-regularity
  + Xi/Phi-specific analytic hypotheses
  => all-degree/all-shift Jensen hyperbolicity
```

or:

```text
signed/indefinite Hankel sign-regularity
  + Xi/Phi-specific analytic hypotheses
  => H_lambda in Laguerre-Polya
```

or:

```text
signed/indefinite Hankel sign-regularity
  + Xi/Phi-specific analytic hypotheses
  => no positive Newman boundary.
```

## Gate E: Theorem-Search Fit/Misfit

Naive route rejected:

```text
A_k is PF-infinity via Toeplitz total positivity.
```

Reason:

The bundle reports robust negative Toeplitz minors for `A_k`.

Correction added 2026-07-03:

This rejects only the naive PF-infinity claim for the exponential/Jensen coefficients `A_k`. The ordinary Taylor coefficients of

```text
F_lambda(s) = sum A_k(lambda) s^k / k!
```

are

```text
c_k(lambda) = A_k(lambda) / k!.
```

The ASW/Edrei PF-sequence route should be tested on `c_k`, not on `A_k`. A later finite exact-rational audit found no negative Toeplitz minors for `c_k` in the tested ranges, so `c_k` PF-infinity is a live route but still requires all-order proof or certified theorem reduction.

Live route:

```text
A_k may encode a signed/indefinite Hankel structure, while c_k = A_k/k!
may be the ordinary coefficient sequence for a classical PF-infinity route.
Both require a missing theorem or certified all-order positivity.
```

Required before manuscript promotion:

```text
identify the transformation or state the missing theorem explicitly
```

## Gate F: Result Language

Use:

```text
reproduced
stress-tested
stable under higher precision
conditional on coefficient enclosure
finite evidence
route hypothesis
missing bridge theorem
```

Do not use:

```text
proved RH
certified, unless Level 4 is met
obstructed complex zeros globally, unless positive Newman boundary is ruled out
rank-2 proof, unless full infinite-dimensional theorem exists
```
