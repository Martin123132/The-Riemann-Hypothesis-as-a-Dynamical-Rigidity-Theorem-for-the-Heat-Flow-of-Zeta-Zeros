# Proof Claim Ledger

Date: 2026-07-11

Status: claim-classification ledger. This is not a proof of RH or `Lambda <= 0`; it records which claims are exact, finite, diagnostic, open, rejected, or hygiene gates.

## Purpose

The proof programme now has many artifacts. This ledger keeps the main claims
separated by type so finite evidence and diagnostics cannot be silently
promoted into an all-order theorem.

Machine-readable ledger:

```text
work/rh_compute/results/proof_claim_ledger.json
```

Validator:

```text
python work/rh_compute/scripts/check_proof_claim_ledger.py
```

Current result:

```text
validated proof-claim ledger: 182 claims, 0 issues, 13 open theorem targets
```

## Categories

```text
exact_lemma:
  exact noncircular identities already available for proof development

finite_certificate:
  promoted finite interval or manifest-backed certificates

interval_certificate:
  promoted interval-backed theorems on a complete real interval or half-line

diagnostic:
  finite necessary-condition or theorem-search diagnostics

algebraic_reindexing:
  exact algebraic translations such as Toeplitz/Jacobi-Trudi

countermodel_gate:
  executable proof-safety obstructions

theorem_target:
  open bridge theorem needed before any proof promotion

forbidden_promotion:
  invalid proof step rejected by countermodel or logic

hygiene_gate:
  reproducibility, language, status, and reference integrity checks
```

## Current Open Targets

```text
target_direct_coefficient_pf:
  prove all Toeplitz minors of c_k(0)=mu_{2k}(0)/(2k)! are nonnegative

target_signed_hankel_jensen_bridge:
  prove that all-order signed-Hankel/sign-consistency structure for A_k,
  with Xi/Phi-specific analytic hypotheses, implies all-degree/all-shift
  Jensen hyperbolicity, Laguerre-Polya membership, or Lambda <= 0
  specification: outputs/signed_hankel_jensen_bridge_target.md

target_jensen_window_pf_bridge:
  prove that every binomially weighted Jensen window
  B^{d,n,0}_j=binom(d,j) A_{n+j}(0) is finite PF-infinity, or prove this
  from all-order shifted reshaped-Hankel sign consistency
  specification: outputs/jensen_window_pf_bridge_target.md

target_schur_positive_specialization:
  construct a positive Schur/Edrei-Thoma specialization h_k -> d_k(0)
  specification: outputs/positive_schur_specialization_target.md

target_positive_determinant_integral:
  derive a positive determinant integral formula for every structurally
  nonzero Toeplitz minor

target_edrei_log_power_representation:
  prove an all-order positive Edrei log-power representation for H_0
```

Each is explicitly recorded as `open_target`. The checker rejects theorem
targets that lack a current blocker or required upgrade.

## New Finite Diagnostic

```text
hankel_sign_consistency_reduction_point_audit:
  exact-rationalized cache point audit for the Grussler-Damm reshaped-Hankel
  finite condition; validates five lambdas, k=2..5, N=18; not an interval
  certificate or all-order sign-consistency proof

hankel_sign_consistency_reduction_finite_certificate:
  Arb/enclosure-backed finite certificate for the same reshaped-Hankel
  frontier; validates 689,795 finite minors for k=2..7, N=20; not an all-order
  sign-consistency theorem

shifted_hankel_sign_consistency_finite_certificate:
  Arb/enclosure-backed finite certificate for shifted reshaped-Hankel blocks;
  staircase manifest validates 3,154,515 finite minors for shifts n=0..20:
  k=2..5 at N=18, k=6 at N=16, k=7 at N=15, and k=8 at N=14; not an
  all-shift or all-order sign-consistency theorem

jensen_hankel_bridge_algebra_gate:
  exact degree-2 signed-Hankel/Jensen identity plus a positive rational
  degree-3 countermodel; blocks promotion from finite low-order reshaped signs
  to Jensen hyperbolicity

jensen_window_pf_obligation_algebra_gate:
  exact low-degree Jensen-window PF obligation algebra; degree 2 matches the
  signed-Hankel threshold, while degree 3 and degree 4 introduce additional
  banded Toeplitz obligations and finite low-order countermodel failures

arb_jensen_window_pf_obligation_diagnostic:
  Arb/enclosure-backed finite diagnostic for selected degree-3/4
  Jensen-window contiguous Toeplitz determinants; validates 1,470/1,470
  finite determinants for the five-lambda grid and shifts n=0..20

arb_jensen_window_sturm_hyperbolicity_diagnostic:
  Arb/Sturm finite diagnostic for selected degree-3/4 Jensen-window
  positive-root counts; validates 210/210 finite root-count rows for the
  five-lambda grid and shifts n=0..20

arb_jensen_window_sturm_d5_hyperbolicity_diagnostic:
  Arb/Sturm finite diagnostic for selected degree-5 Jensen-window
  positive-root counts; validates 105/105 finite root-count rows for the
  five-lambda grid and shifts n=0..20

arb_jensen_window_sturm_d6_d12_hyperbolicity_diagnostic:
  Arb/Sturm finite diagnostic for degree-6 through degree-12 Jensen-window
  positive-root counts; validates 735/735 finite rows on the five-lambda grid
  and shifts n=0..20, bringing the degree-3 through degree-12 total to
  1050/1050 with no failed or inconclusive rows

jensen_window_sturm_pf_finite_consequence:
  finite window-by-window PF consequence of the Arb/Sturm root-count
  manifests; validates 1050/1050 checked Jensen windows across degrees
  d=3..12, five lambdas, and shifts n=0..20

jensen_window_pf_bridge_obligation_ledger:
  theorem-obligation hygiene gate for target_jensen_window_pf_bridge;
  validates 10 exact, finite, open, conditional, rejected, and route-separated
  obligations, with 3 open obligations and finite rows blocked from closing
  the target

jensen_window_pf_theorem_machinery_fit_matrix:
  theorem-search hygiene gate for jwpf_06; validates 7 total-positivity,
  PF, zero-preserver, sign-regularity, downstream, and rejected-shortcut
  theorem-family rows, with 0 ready-to-apply rows and explicit fatal gaps

jensen_window_pf_sign_regular_transfer_gap_matrix:
  finite theorem-search diagnostic combining exact degree-2 signed-Hankel
  contact, degree-3/4 countermodel gates, and all-order/binomial/shift
  uniformity requirements into 9 transfer rows, with 0 ready-to-apply rows

jensen_window_pf_factorial_multiplier_split_audit:
  finite theorem-search diagnostic separating the valid conditional
  Pólya-Schur factorial multiplier step from the false raw moment-window input
  route; validates 315 raw degree-2 anti-hyperbolic rows and 315 normalized
  degree-2 positive rows, with 0 ready-to-apply rows

jensen_window_pf_structural_ansatz_matrix:
  structural proof-search hygiene gate for jwpf_06; validates 6 candidate,
  blocked, and rejected ansatz rows against exact degree-2/3/4 hard tests and
  the finite countermodel kill gate, with 0 ready-to-apply rows

jensen_window_pf_cauchy_binet_low_degree_scout:
  symbolic theorem-search diagnostic for the live Cauchy-Binet ansatz;
  validates 15 degree-2/3/4 hard formulas by adjacent-log-concavity ratio
  parametrization with nonnegative Bernstein coefficients, finds 0 kernel
  identities, and keeps the larger-minor countermodel warning active

jensen_window_pf_log_concavity_frontier_scout:
  symbolic frontier diagnostic extending the low-degree scout to larger
  contiguous Jensen-window Toeplitz minors; validates 14 rows and locates first
  Bernstein-certificate failures at degree 3 size 6 and degree 4 size 5, and
  first exact log-concave countermodel negatives at degree 3 size 8 and degree
  4 size 6

jensen_window_pf_ratio_condition_scout:
  ratio-condition theorem-search diagnostic; validates 7 candidate rows,
  rejecting adjacent log-concavity, decreasing ratio contractions,
  second-order log-concavity, and selected low-degree Bernstein positivity by
  exact countermodel, and contraction log-concavity by constructed positive
  extension

jensen_window_pf_contraction_log_concavity_scout:
  ratio-condition rejection diagnostic; validates a positive rational
  extension satisfying x2^2 >= x1*x3 while the degree 3 size 8 and degree 4
  size 6 contiguous Jensen-window minors remain negative

jensen_window_pf_schur_shape_contract:
  exact shape-contract diagnostic; validates the Jensen-window
  Jacobi-Trudi/Schur shape map on a bounded N=8, order<=5 grid, recording
  15,709 finite-band nonzero bounded shapes and 2 hard frontier column-shape
  obligations with mixed-sign h-monomial expansions

jensen_window_pf_column_recurrence_contract:
  column-shape recurrence diagnostic; validates the elementary-symmetric
  recurrence C_m=h_0^m e_m and confirms the degree 3 size 8 and degree 4 size
  6 hard frontier recurrence rows match the exact negative rational
  countermodel values

jensen_window_pf_column_recurrence_finite_coverage:
  finite coverage diagnostic; validates that checked zeta-window grids support
  the column recurrence target with 1,470 direct positive Arb determinant rows,
  210 hard recurrence rows, 315 checked Sturm/PF windows, and 12,600 positive
  recurrence-only stress rows, without promoting this to an all-order
  recurrence theorem

arb_jensen_window_column_recurrence_stress:
  finite Arb stress diagnostic; validates 12,600 positive Jensen-window column
  recurrence rows for degrees d=3..8, sizes m=1..20, five lambda values, and
  shifts n=0..20

jensen_window_pf_reciprocal_coefficient_extended_stress:
  extended reciprocal coefficient diagnostic; validates 72,600 positive
  normalized [t^m]1/H(-t) rows for degrees d=2..12, sizes m=1..40, five
  lambda values, and shifts n=0..32

jensen_window_pf_reciprocal_positivity_route_matrix:
  theorem-search hygiene gate for the column recurrence target
  [t^m]1/H(-t)>=0; validates 9 reciprocal-positivity route rows, with 0
  ready-to-apply rows, 3 live representation candidates, a rejected standard
  positive Stieltjes/Jacobi fraction route, rejected generic ratio shortcuts,
  and finite stress rows kept finite

jensen_window_pf_reciprocal_fraction_scout:
  reciprocal continued-fraction sign diagnostic; validates 3 symbolic rows
  and 735 finite Arb zeta-window sign rows, showing that the standard positive
  S-fraction/J-fraction route for E(t)=1/H(-t) has the wrong first nontrivial
  sign while leaving signed or modified fractions open

jensen_window_pf_reciprocal_signed_j_fraction_scout:
  reciprocal signed J-fraction diagnostic; validates 2 symbolic rows, 3,675
  finite signed reciprocal-Hankel determinant rows, and 2,940 finite
  signed-lambda rows, supporting the signed modified fraction target without
  proving the all-order theorem

jensen_window_pf_reciprocal_signed_jacobi_beta_scout:
  reciprocal signed Jacobi beta diagnostic; validates 3 symbolic rows and
  3,675 finite beta rows, with 2,940 positive rows, 630 negative beta_1 rows,
  and 105 terminal degree-2 zero-containing rows, sharpening the signed
  modified fraction target without proving the all-order theorem

jensen_window_pf_reciprocal_motzkin_path_obstruction_scout:
  reciprocal Motzkin path obstruction diagnostic; validates 3 symbolic rows,
  735 finite mu_2 cancellation rows, and 630 finite beta_1 diagonal obstruction
  rows, rejecting the raw ordinary J-fraction Motzkin path model as manifest
  positivity while leaving modified signed models open

jensen_window_pf_reciprocal_motzkin_parity_lift_obstruction_scout:
  reciprocal Motzkin parity-lift obstruction diagnostic; validates 3 symbolic
  rows and 5,145 finite same-length mixed-sign witness rows, rejecting global
  length-parity signs and diagonal sign conjugation as repairs of the raw
  ordinary J-fraction Motzkin path model while leaving state-space modified
  models open

jensen_window_pf_signed_j_fraction_theorem_target:
  theorem-target hygiene gate for the signed modified J-fraction route;
  validates 7 fit/misfit rows and states the missing theorem needed to convert
  all-order signed reciprocal-Hankel/Jacobi signature into coefficientwise
  nonnegativity of E(t)=1/H(-t), with 0 ready-to-apply rows

jensen_window_pf_modified_signed_model_target:
  modified signed-model hygiene gate; validates 9 model rows, rejecting the
  raw ordinary Motzkin/J-fraction model, diagonal sign conjugation, and global
  length-parity sign repairs while leaving 4 modified candidates live only
  conditionally, with 0 ready-to-apply rows

jensen_window_pf_state_space_sign_lift_obstruction_scout:
  derived finite obstruction diagnostic; validates 3 symbolic rows and 735
  mu_2 rows showing that an absolute-value sign-state cover of raw Motzkin
  paths overshoots the actual reciprocal coefficient by 2*kappa_1>0, while
  leaving genuinely modified state-space doubled models open

jensen_window_pf_oscillatory_resolvent_fit_matrix:
  oscillatory/resolvent theorem-search hygiene gate; validates 8 fit/misfit
  rows with 0 ready-to-apply rows, rejects ordinary entrywise Jacobi powers,
  diagonal similarity, absolute-value majorants, classical oscillatory
  spectral conclusions, indefinite moment language, and finite signed
  patterns as standalone coefficient-positivity proofs, and leaves only
  positive spectral-transform and Xi/Phi positive-kernel routes live
  conditionally

jensen_window_pf_positive_readout_theorem_target:
  positive-readout theorem-target hygiene gate; validates 8 candidate rows
  with 0 ready-to-apply rows and 2 live foundational routes, isolating the
  exact positive scalar readout obligation to either a noncircular positive
  spectral transform or an Xi/Phi-specific positive resolvent kernel while
  rejecting wrappers, endpoint factorization, finite quadrature, raw signed
  readouts, and absolute-value majorants as standalone proofs

jensen_window_pf_positive_spectral_moment_obstruction:
  positive spectral moment obstruction diagnostic; validates 3 symbolic rows
  and 735 finite Delta_2 obstruction rows, showing that the raw reciprocal
  coefficients mu_m cannot be ordinary power moments of a positive measure,
  while leaving nonordinary positive transforms and Xi/Phi positive kernels
  open

jensen_window_pf_nonordinary_positive_transform_ansatz_matrix:
  nonordinary positive-transform ansatz diagnostic; validates 8 ansatz rows
  with 0 ready-to-apply rows and 3 live ansatz rows, narrowing the surviving
  positive-readout search to non-power positive functionals, Xi/Phi positive
  kernels, or genuinely modified exact state-space transfer models

jensen_window_pf_nonpower_functional_low_degree_scout:
  nonpower functional low-degree diagnostic; validates 7 scout rows with
  0 ready-to-apply rows and 1 live contract row, recomputing reciprocal
  coefficient formulas through mu_6 and signed composition counts through
  m=8 to show that npt_04 needs a genuine positive cone, basis, and
  functional absorbing signed low-degree cancellations

jensen_window_pf_nonpower_functional_cone_candidate_matrix:
  nonpower functional cone-candidate diagnostic; validates 8 cone rows with
  0 ready-to-apply rows and 2 live cone rows, rejecting raw g-coordinate,
  standalone ratio/log-concavity, tautological, and endpoint PF/LP cones while
  leaving Xi/Phi kernel and Cauchy-Binet/determinant-integral cone routes live

jensen_window_pf_cauchy_binet_cone_frontier_matrix:
  Cauchy-Binet cone frontier diagnostic; validates 8 frontier rows with
  0 ready-to-apply rows and 2 live frontier rows, showing that the route
  cannot rest on selected low-degree Bernstein certificates or adjacent
  log-concavity and must instead construct a positive determinant integral
  for the hard column frontier or an all-shape Cauchy-Binet/Andreief kernel

jensen_window_pf_monotone_contraction_frontier_scout:
  monotone-contraction frontier diagnostic; validates 2 exact hard-row
  certificates with 88 positive Bernstein coefficients and 210 finite zeta
  diagnostic rows, showing that the first hard column-frontier polynomials are
  positive under x1 <= x2 <= x3 while the rational log-concavity countermodel
  violates that sharper condition

jensen_window_pf_monotone_contraction_column_extension_scout:
  bounded exact column-extension diagnostic; validates 25 degree-3/4/5 column
  rows with 3,329 positive Bernstein coefficients under monotone contractions,
  including 3 rows beyond the original first hard frontier and a degree-5 band
  through m=8

jensen_window_pf_monotone_contraction_sparse_degree6_scout:
  bounded exact sparse degree-6 diagnostic; validates 10 degree-6 column rows
  through m=10 with 63,347 strictly positive Bernstein coefficients, 0 zero
  coefficients, and 0 negative coefficients under monotone contractions

jensen_window_pf_monotone_contraction_sparse_degree7_frontier_scout:
  bounded exact sparse degree-7 frontier diagnostic; validates 9 degree-7
  column rows through m=9 with 670,891 strictly positive Bernstein
  coefficients, and records the first one-shot global Bernstein certificate
  obstruction at m=10 with 126 negative Bernstein coefficients and minimum
  -4928

jensen_window_pf_monotone_contraction_sparse_degree7_subdivision_scout:
  bounded exact degree-7 m=10 subdivision diagnostic; repairs the one-shot
  global Bernstein obstruction by splitting s0 into [0,1/2], [1/2,3/4], and
  [3/4,1], certifying 785,400 strictly positive slab Bernstein coefficients

jensen_window_pf_monotone_contraction_all_m_counterexample:
  exact countermodel gate; gives the shift-0 infinite static-cone witness
  x=(19/20,19/20,1,1,1,1,1,...) satisfying every pointwise wall and monotone
  contraction, while its normalized degree-7 m=11 column recurrence is
  negative; the propagated static cone alone is not an all-m theorem

jensen_window_pf_monotone_contraction_theorem_target:
  validated on the needed heat regime; the infinite maximum principle proves
  x_(k+1)>=x_k for every k and finite lambda>=-100. The former all-real-lambda
  statement remains false at lambda=-1156, and contractions alone remain
  insufficient for all-m column recurrence positivity

jensen_window_pf_heat_flow_cone_entry_asymptotic_target:
  validated infinite cone-propagation closure; full entry at lambda=-100 and
  the adjacent-defect maximum principle propagate both pointwise walls and
  x_(k+1)>=x_k to every finite lambda>=-100, including lambda=0

jensen_window_pf_phi_taylor_cone_entry_sign_scout:
  finite Taylor sign certificate; validates 4 Phi coefficient balls and the
  two local sign combinations 2*b-a^2<0 and
  2*(a^3-3*a*b+3*c)>0 needed by the fixed-k cone-entry asymptotic route,
  while leaving the uniform-in-k or collared finite cone-entry theorem open

jensen_window_pf_negative_lambda_cone_entry_prefix_scout:
  finite negative-lambda prefix certificate; validates the canonical 69-row
  ACB prefix through lower/upper ratio-cone walls k<=21 and monotone gaps
  k<=20, the extended k30 93-row prefix through lower/upper walls k<=29
  and monotone gaps k<=28, the extended k50 153-row prefix through
  lower/upper walls k<=49 and monotone gaps k<=48, the extended k60
  183-row prefix through lower/upper walls k<=59 and monotone gaps k<=58,
  the extended k80 243-row prefix through lower/upper walls k<=79 and
  monotone gaps k<=78, the extended k100 303-row prefix through
  lower/upper walls k<=99 and monotone gaps k<=98, the extended k150
  453-row prefix through lower/upper walls k<=149 and monotone gaps k<=148,
  and the extended k200 603-row prefix through lower/upper walls k<=199
  and monotone gaps k<=198 for lambda=-25,-50,-100,
  while leaving the all-k tail theorem or finite-collar flow theorem open

jensen_window_pf_negative_lambda_finite_collar_contract:
  finite-collar accounting diagnostic; extracts from the negative-lambda
  prefix and exact ratio-cone collar rule that the canonical usable active
  finite depth is K=19 with collars x_20,x_21, the extended k30 usable
  active finite depth is K=27 with collars x_28,x_29, the extended k50
  usable active finite depth is K=47 with collars x_48,x_49, the
  extended k60 usable active finite depth is K=57 with collars x_58,x_59,
  the extended k80 usable active finite depth is K=77 with collars
  x_78,x_79, the extended k100 usable active finite depth is K=97 with
  collars x_98,x_99, and the extended k150 usable active finite depth is
  K=147 with collars x_148,x_149, and the extended k200 usable active
  finite depth is K=197 with collars x_198,x_199

jensen_window_pf_negative_lambda_tail_barrier_scout:
  finite theorem-search diagnostic; rewrites the missing negative-lambda tail
  as defect inequalities, certifies 63 one-third-width buffer rows and 60
  defect-monotone rows on the canonical finite prefix, certifies 87
  one-third-width buffer rows and 84 defect-monotone rows on the extended k30
  prefix, certifies 147 cone-buffer rows and 144 defect-monotone rows on the
  extended k50 prefix, certifies 177 cone-buffer rows and 174 defect-monotone
  rows on the extended k60 prefix, certifies 237 cone-buffer rows and 234
  defect-monotone rows on the extended k80 prefix, certifies 297 cone-buffer
  rows and 294 defect-monotone rows on the extended k100 prefix, certifies
  447 cone-buffer rows and 444 defect-monotone rows on the extended k150
  prefix, certifies 597 cone-buffer rows and 594 defect-monotone rows on the
  extended k200 prefix, records the k200 one-third-width buffer frontier at
  179/597 rows, and rejects the
  scaled-defect nonincreasing shortcut without promoting the result to an
  all-k tail theorem

jensen_window_pf_negative_lambda_scaled_defect_frontier_scout:
  finite theorem-search diagnostic; validates 597 exact-cone rows on the k200
  prefix, records that the fixed half-width buffer holds on only 521/597 rows
  with first failures at lambda=-50 and k=191 and lambda=-25 and k=133,
  records that the one-third-width buffer holds on only 179/597 rows with
  first failure at lambda=-25 and k=31, and keeps scaled-defect nonincrease
  rejected without promoting any buffer into an all-k theorem

jensen_window_pf_negative_lambda_defect_recurrence_scout:
  finite theorem-search diagnostic; validates the buffered sufficient condition
  on 63 finite rows, validates defect monotonicity on 60 rows, and rejects the
  direct width-preserving recurrence on all 60 checked adjacent rows without
  promoting the result to an all-k tail theorem

jensen_window_pf_negative_lambda_log_curvature_bridge:
  finite theorem-search diagnostic and exact algebraic reduction; translates
  the buffered defect route into the sufficient conditions
  0<=B_k<=2/(3*(2*k+1)) and B_(k+1)<=B_k, validates 63 simple log-buffer rows
  and 60 curvature-monotone rows on the finite prefix, and identifies the
  monotone part with the open Delta^3 log A theorem target

jensen_window_pf_negative_lambda_bounded_log_curvature_target:
  historical finite diagnostic; records the formerly proposed fixed wall
  B_k=-Delta^2 log A_k<=2/(3*(2*k+1)) and its raw moment-curvature equivalent
  on the old k<=22 prefix, but is no longer a live theorem target after the
  repaired k300 obstruction

jensen_window_pf_negative_lambda_bounded_log_curvature_k300_obstruction:
  finite obstruction gate; rewrites the old wall as C_k=(2*k+1)*B_k<=2/3,
  validates that only 179/897 repaired k300 rows satisfy it while 718/897
  fail it with 0 inconclusive rows, and redirects the live route toward the
  zeta-specific raw corridor or linear curvature barrier

jensen_window_pf_negative_lambda_gaussian_curvature_matrix:
  finite theorem-search diagnostic and exact baseline comparison; identifies
  the bounded log-curvature target as a controlled positive deficit from the
  Gaussian raw-moment baseline, validates 63 positive-deficit and
  bounded-deficit rows, and rejects positive Gaussian scale-mixture arguments
  as an upper-wall proof template

jensen_window_pf_negative_lambda_signed_gaussian_perturbation_matrix:
  finite theorem-search diagnostic; packages the surviving fixed-k signed
  Gaussian perturbation route, using the certified Phi Taylor signs to record
  positive leading Gaussian deficit and positive monotone correction while
  preserving the uniform-remainder gap

jensen_window_pf_negative_lambda_uniform_remainder_target:
  open theorem target; states the two-scale uniform remainder theorem needed
  to promote the signed Gaussian perturbation route, records that the
  leading-only cutoff floors 0, 3, and 16 for T=25,50,100 all fall below the
  k>=22 tail start, and leaves both local/mesoscopic k/T remainders and the
  far-tail moving-saddle theorem open

jensen_window_pf_negative_lambda_taylor_moment_budget:
  finite theorem-search diagnostic; derives the exact Gaussian-moment
  normalization for the local q/T route, validates 7 k=22 tail-start samples,
  rejects the low-order Taylor truncation as a finite proof model for
  lambda=-25,-50,-100, and isolates the positivity and log-curvature Taylor
  remainder obligations

jensen_window_pf_negative_lambda_high_order_taylor_scout:
  finite theorem-search diagnostic; automatically generates local Phi Taylor
  polynomials, certifies 8 even coefficient rows through c14, and validates
  35 k=22 truncation rows showing invalid normalizers, upper-wall violations,
  and overbound rows across truncation degrees 6,8,10,12,14

jensen_window_pf_negative_lambda_defect_tail_theorem_target:
  open theorem target; with finite evidence through x_149, states the missing
  all-k defect-tail theorem: prove 0<=d_k<=2/(2*k+1) for k>=150 and
  d_(k+1)<=d_k for k>=149 without endpoint PF/RH/Newman assumptions

jensen_window_pf_negative_lambda_half_width_tail_target:
  finite diagnostic/rejected-route record; finite-rejects the fixed
  half-width scaled-defect route at k150: s_k<=1/2 holds on only 430/447
  checked rows, first fails at lambda=-25 and k=133, and must be replaced by
  an exact-cone or adaptive scaled-defect target plus the separate monotone
  defect bridge

jensen_window_pf_negative_lambda_adaptive_scaled_defect_target:
  open theorem target; replaces the rejected fixed half-width route with the
  exact cone 0<=s_k<=1 or an explicit adaptive envelope from k>=200, plus the
  separate monotone defect bridge from k>=199. The k200 frontier validates
  597/597 exact-cone rows while half-width fails on 76 rows and one-third
  fails on 418 rows, but this remains finite evidence, not an all-k theorem

jensen_window_pf_negative_lambda_adaptive_envelope_matrix:
  finite theorem-search diagnostic; validates the k200 monotone-envelope
  pattern supporting the adaptive target: 594/594 adjacent k-increase rows,
  398/398 cross-lambda order rows, maximum checked scaled defect
  0.5376643171065356005 at lambda=-25 and k=199, and 76 finite half-width
  failures. It does not prove k-uniform or continuous-lambda monotonicity

jensen_window_pf_negative_lambda_adaptive_envelope_obligations:
  exact algebraic obligation diagnostic; validates 9 rows separating the
  adaptive target into exact lower-threshold input, the open upper wall
  x_k<=1, the open monotone bridge x_(k+1)>=x_k, the exact scaled
  k-monotone identity 2+(2*k+1)*x_k-(2*k+3)*x_(k+1)>=0, and the rejected
  fixed half-width/one-third buffers

jensen_window_pf_negative_lambda_raw_moment_bridge_matrix:
  exact finite theorem-search diagnostic; translates the adaptive envelope
  route into raw moment ratios R_k=M_(k+1)*M_(k-1)/M_k^2. It validates the
  exact identities x_k=((2*k-1)/(2*k+1))*R_k, 0<=s_k<=1 iff
  1<=R_k<=(2*k+1)/(2*k-1), and the adaptive R_(k+1) corridor, then checks
  597/597 raw-cone rows and 594/594 corridor rows on the k200 prefix while
  preserving the 76 half-width and 418 one-third fixed-buffer failures

jensen_window_pf_negative_lambda_raw_ratio_decrement_corridor_scout:
  exact finite theorem-search diagnostic; rewrites the raw adaptive corridor
  as the decrement recurrence
  2*(R_k-1)/(2*k+1)<=R_k-R_(k+1)<=4*R_k/(2*k+1)^2. It validates 594/594
  decrement-corridor rows, 591/591 theta-k monotonicity rows, and 396/396
  theta lambda-order rows on the k200 prefix, while two exact raw-cone
  monotone counterexamples block using raw decrease as a shortcut

jensen_window_pf_negative_lambda_k300_precision_repair_audit:
  finite precision-repair diagnostic; extends the raw-ratio decrement route
  to a repaired k300 stress. The broad dps160/cutoff6 run produces
  lambda=-100 high-k precision alarms, while local dps220/cutoff7 repairs
  over k220..250 and k245..320 restore 897/897 raw wall rows, 894/894
  decrement-corridor rows, 891/891 theta-k monotonicity rows, and 596/596
  theta lambda-order rows. This remains finite evidence only

jensen_window_pf_negative_lambda_raw_log_decrement_bridge:
  exact finite theorem-search diagnostic; rewrites the raw decrement corridor
  in log-ratio coordinates with p_k=log(R_k) and
  delta_k=p_(k+1)-p_k. It proves the exact equivalence to the two-sided
  delta_k bounds, validates 897/897 raw-log wall rows and 894/894
  log-corridor rows on repaired k300 data, and keeps two exact raw-cone
  counterexamples blocking raw-log decrease as a shortcut

jensen_window_pf_negative_lambda_coefficient_curvature_corridor_bridge:
  exact finite theorem-search diagnostic; rewrites the raw/log decrement
  route in coefficient-curvature variables
  B_k=-log(((2*k-1)/(2*k+1))*R_k). It proves the exact equivalence to
  log((2*k+3)/(2+(2*k+1)*exp(-B_k)))<=B_(k+1)<=B_k, validates 897/897
  B-wall rows and 894/894 curvature-corridor rows on repaired k300 data,
  and blocks monotone curvature or raw walls as standalone shortcuts

jensen_window_pf_negative_lambda_linear_curvature_barrier_scout:
  exact finite theorem-search diagnostic; proves the exact sufficient lemma
  L_k(B)<=((2*k+1)/(2*k+3))*B for B>=0, reducing the nonlinear lower
  curvature barrier to the stronger linear target
  B_(k+1)>=((2*k+1)/(2*k+3))*B_k. It validates 897/897 B-wall rows and
  894/894 linear-barrier rows on repaired k300 data while keeping the
  analogous defect-width recurrence rejected

jensen_window_pf_negative_lambda_scaled_curvature_monotonicity_target:
  open theorem target; names the replacement all-k theorem
  C_(k+1)>=C_k, exactly equivalent to
  B_(k+1)>=((2*k+1)/(2*k+3))*B_k. It validates 894/894 repaired k300
  finite stress rows, records that the retired fixed wall C_k<=2/3 fails on
  718/897 rows, and must still be paired with the monotone-contraction upper
  side B_(k+1)<=B_k before the raw-corridor route can close

jensen_window_pf_negative_lambda_scaled_curvature_log_ceiling_bridge:
  exact finite theorem-search diagnostic; rewrites C_(k+1)>=C_k as the
  affine log-ratio ceiling
  delta_k<=h_(k+1)-((2*k+1)/(2*k+3))*h_k-2*p_k/(2*k+3). It validates
  894/894 repaired k300 scaled-ceiling rows, 894/894 scaled-log-corridor
  rows, and records that the affine ceiling is sharper than the nonlinear
  raw-log upper wall, with only a tight high-k slack left in the finite data

jensen_window_pf_negative_lambda_relative_gaussian_curvature_bridge:
  exact finite theorem-search diagnostic; rewrites B_k as the negative second
  difference of the relative-Gaussian log moment sequence f_k and rewrites
  C_(k+1)>=C_k as the weighted four-point inequality
  (2*k+1)*f_(k-1)-(6*k+5)*f_k+(6*k+7)*f_(k+1)-(2*k+3)*f_(k+2)>=0. It
  validates repaired k300 B-positive, B-decrease, C-increase, and
  C-lambda-order rows while leaving the weighted all-k theorem open

jensen_window_pf_negative_lambda_relative_gaussian_taylor_stencil_scout:
  finite theorem-search diagnostic; certifies positive fixed-k Taylor leading
  signs for B_k, B_k-B_(k+1), and C_(k+1)-C_k, but validates that finite
  Taylor truncations are unstable proof objects: only 4/35 sampled truncation
  rows are all-positive, so a uniform Taylor-tail remainder theorem is still
  required

jensen_window_pf_negative_lambda_relative_gaussian_stencil_remainder_obligations:
  exact finite theorem-search diagnostic; decomposes the Taylor-tail log
  error into exact B, companion, and weighted-gap epsilon stencils, records
  4/35 positive baseline rows and 31 blocked baseline rows, and identifies
  the weakest finite half-margin budget 1.166490564421582442E-8. This
  sharpens, but does not prove, the uniform remainder theorem

jensen_window_pf_negative_lambda_relative_gaussian_pointwise_tail_budget:
  exact finite theorem-search diagnostic; converts the exact epsilon-stencil
  margins into pointwise log-tail and multiplicative relative-tail tolerances,
  with weakest half-safety log-tail envelope 1.458113205526978052E-9 and
  relative-tail ratio bound 1.458113204463930993E-9. These are required
  tolerances, not proved analytic tail estimates

jensen_window_pf_negative_lambda_relative_gaussian_next_increment_stencil_stress:
  finite theorem-search diagnostic; shows that the known next Taylor increment
  exceeds the crude pointwise half-safety budget on both tested positive
  baselines while the structured B, companion, and weighted-gap stencil signs
  are preserved. This rejects the current pointwise triangle shortcut and
  sharpens the live direct signed stencil-tail route

jensen_window_pf_negative_lambda_relative_gaussian_degree16_stencil_continuation:
  finite theorem-search diagnostic; computes the degree-16 Taylor coefficient
  and tests all four current positive baselines one Taylor step further. Three
  preserve structured stencil signs, while the degree-14 T=1000 baseline fails
  through the companion stencil and the degree-14 T=2000 baseline survives,
  sharpening the route toward an explicit large-T or q/T collar

jensen_window_pf_negative_lambda_relative_gaussian_degree16_collar_scan:
  finite theorem-search diagnostic; scans every integer T from 900 to 2200 at
  k=22 for the degree-16 continuation. The M=7 baseline first has positive
  stencils at T=918, the M=8 continuation first preserves signs at T=1156, the
  half-safety condition first holds at T=1483, and the pointwise budget fails
  on all 1283 baseline-positive rows

jensen_window_pf_negative_lambda_relative_gaussian_degree16_real_t_collar_scout:
  finite theorem-search diagnostic; certifies, for the rationalized fixed-k
  degree-16 finite surrogate only, that the four normalizers and all three
  structured stencil signs persist on the real half-line T>=1156. The remaining
  upgrade is interval coefficient control plus signed infinite-tail stencil
  bounds on the same collar

jensen_window_pf_negative_lambda_relative_gaussian_degree16_arb_real_t_collar_certificate:
  finite theorem-search diagnostic; upgrades the fixed-k degree-16 real-T
  surrogate collar to Arb coefficient-ratio balls. Bernstein coefficients
  certify F_21 through F_24, the stripped B product, the stripped companion
  product, and the stripped weighted-gap derivative numerator on
  0<=u<=1/1156, while leaving the infinite residual Taylor-tail theorem open

jensen_window_pf_negative_lambda_relative_gaussian_degree40_arb_collar_ladder_stress:
  finite theorem-search diagnostic; repeats the Arb real-T collar test for
  every even finite surrogate degree 16 through 40 at fixed k=22. All tested
  levels certify positive F_21 through F_24 normalizers and all three
  structured stencils on T>=1156 with zero Bernstein failures, while leaving
  the residual Taylor-tail bound beyond degree 40 and all-k upgrade open

jensen_window_pf_negative_lambda_relative_gaussian_degree40_residual_tail_budget:
  exact finite theorem-search diagnostic; converts the degree-40 Arb collar
  margins into sufficient fixed-k residual targets on 0<=u<=1/1156:
  |R_i(u)|<=5.382819486765314521E-01*u^3 and
  |R_i'(u)|<=9.315354075509573936E-03*u for i=21..24. The finite
  degree-42..80 tail profile consumes less than 0.1% of those budgets, but no
  analytic majorant for the infinite residual tail is proved

jensen_window_pf_negative_lambda_relative_gaussian_formal_tail_obstruction_scout:
  finite theorem-search obstruction; extends the formal residual-term profile
  through degree 240 and shows that after decreasing to a least-term region
  near j=103, the formal value and derivative terms grow again from j=103 to
  j=104 for F_21 through F_24. This rejects monotone/geometric infinite
  formal-tail summation and fixed-radius Cauchy termwise summation as proof
  templates, while leaving the actual asymptotic remainder theorem open

jensen_window_pf_negative_lambda_relative_gaussian_asymptotic_remainder_target:
  exact theorem-search diagnostic; converts the degree-40 residual budget and
  formal-tail obstruction into explicit sufficient analytic-remainder
  constants. A common 1000x first-omitted-term theorem would fit inside the
  half-safety value and derivative budgets for F_21 through F_24 on
  0<=u<=1/1156, with common multiplier limit about 1419.939. This is target
  calibration only; the first-omitted-term and least-term remainder theorems
  remain open

jensen_window_pf_negative_lambda_relative_gaussian_actual_endpoint_remainder_scout:
  floating endpoint theorem-search diagnostic; evaluates the actual
  relative-Gaussian multiplier at T=1156 by generalized Gauss-Laguerre
  quadrature. At the selected order N=320, the value and derivative residuals
  for F_21 through F_24 are below one first omitted formal term and below 0.1%
  of the degree-40 half-safety budgets. This is endpoint evidence only; it is
  not interval-certified and does not prove a uniform collar remainder theorem

jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout:
  floating cancellation-reduced theorem-search diagnostic; subtracts the
  degree-40 polynomial inside the Gamma expectation and samples
  T=1156,1500,2000,5000,10000 for F_21 through F_24. Across four quadrature
  orders, all sampled value and derivative residuals are below one first
  omitted formal term, with worst value ratio about 0.9707100590 and worst
  derivative ratio about 0.9693567775. This is finite floating grid evidence
  only, not interval-certified and not a uniform collar theorem

jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target:
  open numerical-certification target; converts the cancellation-reduced
  grid slack into explicit intervalization obligations. A future certificate
  with total ratio error below 1.0e-2 would keep the finite grid below one
  first omitted formal term, but Laguerre node/weight intervals, Phi n-tail
  bounds, quadrature-remainder bounds, coefficient propagation, and a
  grid-to-collar bridge remain open

jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_bound_scout:
  analytic padded-range theorem-search diagnostic; bounds the n>30 tails in
  Phi, Phi', and Phi(0) on 0<=x<=1 far below the 2.0e-3 per-source
  intervalization cap after normalizing by the diagnostic Phi(0)>=0.44 proxy.
  This only narrows the nlrgit_03 tail source; interval Laguerre-node range
  and Phi(0) lower-bound certificates, plus quadrature, coefficient, rounding,
  and grid-to-collar obligations, remain open

jensen_window_pf_negative_lambda_relative_gaussian_node_c0_range_certificate:
  exact/Arb theorem-search diagnostic; certifies the two finite-grid side
  conditions used by the padded Phi-tail scout. Gershgorin/AM-GM bounds every
  recorded generalized Laguerre node by 809<T_min=1156, hence x<=1, and the
  n=1 term in Phi(0) is Arb-certified above 0.44. This does not certify
  individual Laguerre weights, quadrature remainders, coefficient propagation,
  rounding, or the grid-to-collar bridge

jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_grid_certificate:
  finite-grid Phi-tail source certificate; composes the padded n>30 Phi,
  Phi', and Phi(0) tail majorants with the certified side conditions x<=1 and
  Phi(0)>=0.44. The three omitted-n tail sources sit far below the 2.0e-3
  per-source intervalization cap for the recorded finite grid. This does not
  certify finite n<=30 node evaluations, Laguerre weights, quadrature error,
  coefficient propagation, rounding, or the grid-to-collar bridge

jensen_window_pf_negative_lambda_relative_gaussian_quadrature_ladder_scout:
  high-order floating theorem-search diagnostic; recomputes the worst
  cancellation-reduced grid row T=10000, F_21 at quadrature orders 96 through
  320. All sampled ratios remain below one first omitted formal term, with
  order-spread below 1e-14, calibrating a future rigorous quadrature-radius
  target of 1e-6. This is not a quadrature-remainder theorem or interval
  weight certificate

jensen_window_pf_negative_lambda_relative_gaussian_quadrature_remainder_route_matrix:
  exact/formula route diagnostic; records the classical N=320 generalized
  Gauss-Laguerre remainder prefactor for alpha=41/2, about 2.79e-159,
  translates the 1e-6 quadrature ratio-radius target into explicit value and
  derivative 640th-derivative supremum caps, and verifies that adding that
  cap to the worst-row finite-plus-tail budget still keeps both ratios below
  one first omitted term. This is not a derivative bound, not interval
  adaptive integration, and not a quadrature-remainder theorem

jensen_window_pf_negative_lambda_relative_gaussian_worst_row_far_tail_split_certificate:
  Arb worst-row far-tail split diagnostic; certifies the finite n<=30
  cancellation-reduced continuum tail y>=200 for T=10000, F_21 using
  monotone Phi majorants and upper incomplete-Gamma moment bounds. The value
  and derivative tail ratios are about 8.66e-26 and 4.50e-27, far below the
  1e-6 quadrature ratio-radius target. Compact integration on 0<=y<=200,
  aggregation, all-row coverage, and grid-to-collar coverage remain open

jensen_window_pf_negative_lambda_relative_gaussian_worst_row_compact_interval_integration_scout:
  Arb worst-row compact-interval diagnostic; imports the y>=200 far-tail
  split and tests six raw Arb panel hulls on the remaining compact interval
  0<=y<=200 for T=10000, F_21. The raw value and derivative width-to-cap
  ratios are about 1.38e39 and 4.28e37, so plain interval-Riemann hulls are
  rejected. The live route is now a local Taylor/Chebyshev panel model with
  exact Gamma-weighted moments; compact integration remains open

jensen_window_pf_negative_lambda_relative_gaussian_worst_row_chebyshev_panel_moment_scout:
  high-precision floating Chebyshev panel-moment diagnostic; integrates
  Chebyshev interpolants against incomplete-Gamma panel moments on the same
  six compact panels 0<=y<=200. Consecutive degree deltas from 16->20 onward
  are below the unscaled quadrature caps in both value and derivative
  channels, with reference degree-32 estimates about -6.58e-34 and -1.38e-32.
  This calibrates the local-model route but does not provide Arb coefficient
  balls or interpolation-remainder bounds

jensen_window_pf_negative_lambda_relative_gaussian_worst_row_arb_chebyshev_interpolant_moment_scout:
  Arb Chebyshev interpolant-moment diagnostic; promotes the floating
  Chebyshev ladder to Arb-enclosed interpolant coefficients and
  incomplete-Gamma panel moments for degrees 16, 20, 24, and 32. All
  consecutive interpolant Cauchy deltas are below the unscaled quadrature
  caps, but this still certifies only interpolant arithmetic; the
  interpolation remainder between the interpolants and the true compact core
  remains open

jensen_window_pf_negative_lambda_relative_gaussian_worst_row_interpolation_remainder_route_matrix:
  interpolation-remainder route diagnostic; quantifies the remaining compact
  panel theorem after Arb interpolant arithmetic. The heaviest Gamma panel is
  20<=y<=50 with mass upper about 0.6021615933248953, and the artifact records
  20 Bernstein degree/rho budgets plus 16 minimal-degree rows for the analytic
  domain and sup-norm certificates needed to bound true-function interpolation
  remainders. This is route sizing only; no analytic-domain, sup-norm,
  endpoint, Taylor-model, or compact-integral certificate is proved

jensen_window_pf_negative_lambda_relative_gaussian_endpoint_parity_repair_matrix:
  endpoint parity-repair diagnostic; makes the y=0 branch obligation explicit
  for the compact interpolation route. Arb series expansion of the finite
  n<=30 Phi truncation records 8 low-order odd Taylor rows through degree 15,
  with first odd coefficient absolute upper about 1.50e-1300. This does not
  prove endpoint analyticity: the admissible repairs are exact infinite-kernel
  evenness plus certified tail charge, or an x-variable endpoint-panel
  certificate

jensen_window_pf_negative_lambda_relative_gaussian_endpoint_x_panel_route_matrix:
  endpoint x-panel route diagnostic; quantifies the first-panel repair route
  after the parity matrix. The change y=T*x^2 maps 0<=y<=1 to 0<=x<=0.01,
  the transformed Gamma density has x^42 power for alpha=20.5, and the
  first-panel mass upper is about 1.62e-21. The artifact records 18
  Bernstein degree/rho budgets for a future x-domain and sup-norm certificate,
  but does not prove the exact x moments, x-panel interpolation remainder, or
  compact interval certificate

jensen_window_pf_negative_lambda_relative_gaussian_endpoint_x_moment_taylor_certificate:
  finite Arb endpoint certificate for the T=10000, F_21 worst row. It expands
  the finite n<=30 cancellation-reduced value and derivative cores through
  degree 64, integrates all 65 Taylor rows by exact transformed Gamma moments,
  and bounds the true-function tail by a rigorous |z|<=0.2 Cauchy majorant.
  The resulting first-panel integrals are certified negative, with both
  absolute-to-global-cap ratios below 4.10e-47. This closes only 0<=y<=1 for
  the finite core; the separate n>30 tail, five y>=1 panels, all-row coverage,
  and uniform collar theorem remain open

jensen_window_pf_negative_lambda_relative_gaussian_worst_row_compact_x_moment_taylor_certificate:
  finite Arb compact certificate for the T=10000, F_21 worst row. The complete
  range 0<=y<=200 maps into x<=sqrt(0.02), inside the certified |z|<=0.38
  disk. A degree-128 Taylor model is integrated by 129 exact transformed
  Gamma moments, with true-function remainder radii consuming less than
  4.0e-10 and 1.3e-9 of the value and derivative caps. This replaces the
  six-panel interpolation obligation for the finite n<=30 core. The separate
  n>30 Phi tail, y>200 far tail, all-row coverage, and uniform collar theorem
  remain open

jensen_window_pf_negative_lambda_relative_gaussian_worst_row_full_expectation_certificate:
  complete true-expectation certificate for the T=10000, F_21 row. It
  composes the degree-128 compact finite core, the finite-core y>=200 tail,
  and a new Arb global n>=31 Phi/Phi'/Phi(0) normalization correction. The
  full value and derivative balls are negative, with first-omitted ratio
  uppers about 0.9707101 and 0.9693568 and margins above 0.029. No worst-row
  integral source remains open and generalized Gauss-Laguerre quadrature is
  not used in the final row certificate. The other finite-grid rows,
  all-source aggregation, and finite-grid-to-collar theorem remain open

jensen_window_pf_negative_lambda_relative_gaussian_all_row_direct_expectation_certificate:
  complete recorded-grid direct expectation certificate. A uniform
  degree-384 exact-moment/Cauchy core on x<=8/25, rowwise incomplete-Gamma
  real-tail bounds, and the global n>=31 normalization correction certify all
  20 value and all 20 derivative expectation balls negative and below one
  first omitted term. The worst ratios remain at T=10000, F_21, below
  0.970711 and 0.969357. This retires the recorded-grid integration source,
  but signed finite-degree stencil aggregation, interval coverage in T, the
  real-T collar, sign-regularity bridge, RH, and Lambda <= 0 remain open

jensen_window_pf_negative_lambda_relative_gaussian_recorded_grid_stencil_composition_certificate:
  fixed-k=22 recorded-grid stencil certificate. The 20 expectation balls fit
  the exact rational budgets |R_i|<=(1/2)u^3 and |R_i'|<=(9/1000)u, using
  less than 5.7e-4 of either budget. A fresh Arb perturbation ledger retains
  positive normalizer, B-product, companion-product, and weighted-gap
  derivative margins, certifying all five recorded T systems. The source
  floating thresholds are not proof inputs. Real-T interval coverage,
  remaining k windows, cone entry, sign regularity, RH, and Lambda <= 0
  remain open

jensen_window_pf_negative_lambda_relative_gaussian_finite_collar_segment_stencil_certificate:
  real-interval certificate on 1156<=T<=10000. A two-regime compact Cauchy
  and incomplete-Gamma majorant proves the exact rational residual budgets
  uniformly for F_21 through F_24, then composes them with the degree-40 Arb
  perturbation ledger. The fixed-k=22 stencil system is certified throughout
  the bounded segment. The ray T>10000, remaining k windows, cone entry, sign
  regularity, RH, and Lambda <= 0 remain open

jensen_window_pf_negative_lambda_relative_gaussian_full_kernel_evenness_cauchy_lemma:
  exact analytic lemma. The Jacobi theta transformation and differential
  operator covariance prove the full infinite Riemann kernel even and
  analytic on |z|<pi/8. Subtracting normalized even Taylor terms through
  degree 40 leaves an exact order-42 zero and x^42-factored Cauchy bounds.
  This is not inferred from finite-truncation near-evenness and is not by
  itself a Gamma-expectation, cone-entry, RH, or Lambda <= 0 theorem

jensen_window_pf_negative_lambda_relative_gaussian_full_real_T_fixed_k_stencil_certificate:
  full real-interval certificate for fixed k=22. The order-42 factor, a full
  infinite-kernel disk majorant, full real-tail majorants, geometric n-tail
  bounds, and upper-Gamma hazard monotonicity certify the residual budgets on
  T>=10000. Together with the bounded segment, all four normalizers and the
  three target stencils are certified for every real T>=1156. All-k coverage,
  cone entry, sign regularity, RH, and Lambda <= 0 remain open

jensen_window_pf_negative_lambda_relative_gaussian_first_omitted_denominator_certificate:
  Arb denominator-side theorem-search diagnostic; certifies the first omitted
  coefficient r_21 is negative and bounded away from zero, then proves every
  value and derivative first-omitted denominator on the recorded finite grid is
  positive. The worst lows occur at T=10000, F_21, translating a 1e-6
  ratio-radius target into scaled absolute-radius caps about 6.78e-28 and
  1.42e-30. This does not certify the residual numerator, quadrature error,
  rounding, or the grid-to-collar bridge

jensen_window_pf_negative_lambda_relative_gaussian_coefficient_core_certificate:
  Arb coefficient-core propagation diagnostic; rebuilds coefficient-ratio
  balls r_0 through r_21 and propagates their radii through exact Gamma
  moments on the recorded finite grid. The worst value and derivative
  coefficient-radius ratios to the first omitted denominators are about
  8.61e-81 and 5.52e-83 at T=10000, F_21, far below the 1e-6 ratio-radius
  target. This retires the coefficient-ball source for finite-grid
  intervalization only; Phi node evaluation, Laguerre node/weight intervals,
  quadrature error, rounding aggregation, and grid-to-collar coverage remain
  open

jensen_window_pf_negative_lambda_relative_gaussian_worst_row_laguerre_root_bracket_certificate:
  Arb worst-row Laguerre node diagnostic; certifies 320 disjoint
  sign-changing brackets for the roots of L_320^(41/2) at the T=10000, F_21
  quadrature row. The widest bracket has width about 4.15e-17. The same row
  records 30 zero floating SciPy weights, so the remaining proof work is a
  non-floating Christoffel-weight interval certificate before this can become
  a quadrature-row interval enclosure

jensen_window_pf_negative_lambda_relative_gaussian_worst_row_christoffel_weight_midpoint_scout:
  Arb worst-row Christoffel-weight midpoint diagnostic; evaluates the
  generalized Laguerre weight formula at midpoints of the certified N=320
  root brackets. It gives 320 positive non-floating midpoint weights, repairs
  30 SciPy double-weight underflows, and matches Gamma(43/2) to relative error
  about 1.8e-18. Direct interval denominator evaluation still contains zero on
  all 320 rows, so the midpoint artifact itself remains non-promoted; the
  follow-up interval certificate closes the weight source for this one row

jensen_window_pf_negative_lambda_relative_gaussian_worst_row_christoffel_weight_interval_certificate:
  Arb worst-row Christoffel-weight interval diagnostic; evaluates
  L_321^(41/2) on each certified L_320^(41/2) root bracket by a centered
  Taylor enclosure. It separates zero in all 320 denominator intervals,
  certifies 320 positive Christoffel-weight intervals, repairs the 30 SciPy
  double-weight underflows as interval weights, and has a weight-sum interval
  containing Gamma(43/2). Phi/Phi' node evaluation, quadrature remainder,
  rounding aggregation, all-row coverage, and grid-to-collar coverage remain
  open

jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_part_weighted_sum_interval_certificate:
  Arb worst-row finite-part weighted-sum diagnostic; refines the T=10000,
  F_21, N=320 Laguerre node brackets to 120 bisection steps, evaluates the
  finite n<=30 Phi/Phi' terms on those intervals, sums them with certified
  Christoffel weights, and subtracts the polynomial part through exact Gamma
  moments. The value and derivative residual ratios are certified below one
  first omitted term. The n>30 tail composition, quadrature remainder,
  rounding aggregation, all-row coverage, and grid-to-collar coverage remain
  open

jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_plus_tail_budget_certificate:
  Arb worst-row finite-plus-tail budget diagnostic; composes the certified
  finite n<=30 weighted-sum ratio uppers with the full 2.0e-3 Phi-tail
  source cap for T=10000, F_21, N=320. The composed value and derivative
  ratio uppers are 0.9853957992836557769419015895036210773888 and
  0.9714055674762067320093698741711260875260, both below one first omitted
  term. Quadrature remainder, rounding aggregation, all-row coverage, and
  grid-to-collar coverage remain open

jensen_window_pf_negative_lambda_raw_moment_obstruction_matrix:
  exact countermodel gate; validates three positive two-atom Stieltjes moment
  witnesses showing that generic raw moment positivity/log-convexity cannot
  prove the upper raw wall, scaled-upper corridor side, or monotone-bridge
  lower corridor side. It blocks generic moment-sequence promotion while
  leaving the zeta-specific all-k theorem open

jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target:
  open theorem target now localized at lambda=-100. Berwald-Borell proves the
  upper raw wall, the repaired finite collar reaches k=318, and all n>=2
  moments are controlled by 2/k^6. Tail corridor occupancy reduces to the
  dominant first-summand saddle wall

jensen_window_pf_heat_flow_monotone_closure_scout:
  finite heat-flow closure diagnostic; validates 4 exact lambda-flow algebra
  rows, 315 finite Arb threshold rows, and 305 finite Arb flow-bracket rows,
  using the exact boundary-threshold lemma to isolate the remaining global
  flow-invariance and adjacent-log-concavity gaps without promoting them to an
  analytic theorem

jensen_window_pf_heat_flow_ratio_cone_invariance_lemma:
  exact conditional ratio-cone invariance lemma; proves inward-pointing
  boundary algebra for (2*k-1)/(2*k+1)<=x_k<=1 and x_{k+1}>=x_k under the
  heat-flow ratio ODE, while leaving the zeta cone-entry and infinite-flow
  legitimacy theorems open

jensen_window_pf_heat_flow_boundary_threshold_lemma:
  exact boundary-threshold lemma; proves from Phi positivity and
  Cauchy-Schwarz raw-moment log-convexity that
  x_k >= (2*k-1)/(2*k+1), hence the heat-flow boundary threshold
  x_k >= (2*k-1)/(2*k+5), while leaving adjacent log-concavity and global
  monotone contractions open

jensen_window_pf_kernel_mellin_upper_wall_certificate:
  interval-backed theorem certificate; proves y->Phi(sqrt(y)) strictly
  log-concave by a 200-subinterval Arb/Cauchy compact proof plus an analytic
  full-kernel ray. Berwald-Borell then proves x_k(lambda)<=1 for every real
  lambda and every k>=1. The adjacent-k wall x_(k+1)>=x_k remains open

jensen_window_pf_log_concave_mellin_monotone_wall_countermodel:
  exact interval countermodel gate; the log-concave density
  exp(-5y)*1_[0,1](y) has Gamma-normalized Mellin contractions x_1 and x_2
  inside the upper wall but x_2<x_1 by more than 0.027. Generic
  log-concavity therefore cannot prove the remaining adjacent-k wall

jensen_window_pf_negative_lambda_t1156_monotone_wall_counterexample_certificate:
  interval zeta-kernel counterexample gate; rigorous ACB enclosures of
  A_119..A_122 at lambda=-1156 certify x_121-x_120<-1.68e-8 while both
  contractions remain inside the Mellin upper wall. This blocks all-k cone
  entry at T=1156 and any universal promotion of the fixed-k T>=1156 theorem,
  while leaving moderate-parameter finite-collar plus eventual-k routes open

jensen_window_pf_negative_lambda_kernel_summand_shift_lemma:
  exact all-n kernel-shift lemma with a compact interval theorem; proves
  phi_n(u)=n^(-1/2)phi_1(u+(log n)/2) and, at lambda=-100 and k>=300,
  bounds the complete shifted n=2..20 contribution on v<=3/2 below
  2.122e-29 of the n=1 moment. The later all-k dominance theorem discharges
  the shifted far tail; the dominant n=1 adjacent wall remains open

jensen_window_pf_negative_lambda_first_summand_dominance_certificate:
  all-k analytic kernel-tail theorem; exact ratio monotonicity, strict
  first-integrand concavity, and 15 Arb-positive propagation gates prove
  M_k=M_k^(1)(1+delta_k), 0<=delta_k<=2/k^6 for every k>=300. The adjacent
  log-wall perturbation is at most 16/(k-1)^6 for k>=301

jensen_window_pf_negative_lambda_m100_k320_collar_extension_certificate:
  finite Arb collar extension; promotes the repaired A_245..A_320 source into
  76 positive coefficients, 74 ratio-cone rows, and 73 adjacent-wall rows.
  Nineteen newly exposed rows extend the lambda=-100 collar through k=318,
  with minimum new log gap above 3.709e-6

jensen_window_pf_negative_lambda_first_summand_cumulant_bridge:
  exact conditional reduction; writes the cubic moment difference as an
  average of kappa_3,t(2 log U), extracts the exact positive Gamma wall, and
  proves that kappa_3,t(2 log U)>=-37/(50 t^2) for t>=318 is sufficient for
  L_k^(1)>=1/(4 k^2) for every k>=319. The uniform cumulant bound remains open

jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate:
  global interval theorem for the leading cumulant asymptotic; symbolic
  differentiation, 40,740 Arb compact intervals, and an exact analytic ray
  prove caps 13/20, 1/100, and 1/1000 through fifth order for all t>=318.
  The remaining sufficient target is a seventh-order floor -79/1000

jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate:
  large compact all-real interval theorem; exact paired moments, 40,736 Arb
  eighth-derivative envelope intervals, rigorous midpoint errors, and explicit
  two-sided tails prove the normalized remainder floor -79/1000 throughout
  0.9264<=u<=5, reaching t>1.5241916613e10. The same 4,074 Arb blocks certify
  strict negative third log-cumulant, with maximum upper endpoint below
  -9.13e-6. The separate analytic ray certificate completes the half-line

jensen_window_pf_negative_lambda_first_summand_paired_remainder_ray_certificate:
  analytic asymptotic-ray theorem; an adaptive sqrt(8 log q) window,
  first-order Gaussian moment comparison, and explicit tails prove
  H_t>=-3/250 on u>=5. Combined with the compact theorem this closes the
  global cumulant hypothesis for every t>=318

jensen_window_pf_negative_lambda_first_summand_saddle_wall_target:
  validated dominant-summand wall closure; exact strict-concavity geometry gives
  one saddle in an all-k bracket, and nine 60-digit samples through k=20000
  satisfy L_k^(1)>1/(4*k^2) with increasing k^2 L_k^(1). The exact cumulant
  bridge and paired compact/ray theorems prove L_k^(1)>=1/(4*k^2) for every
  k>=319. The dominance bound and finite collar then prove the full
  lambda=-100 adjacent wall

jensen_window_pf_negative_lambda_m100_full_cone_entry_certificate:
  interval/analytic theorem composition; a precedence-merged repaired source
  certifies 321 positive coefficients, 319 pointwise cone rows, and 318
  adjacent prefix rows. Exact all-k pointwise walls and the analytic k>=319
  adjacent tail prove full infinite ratio-cone entry at lambda=-100; the
  separate maximum-principle certificate closes forward-flow legitimacy

jensen_window_pf_heat_flow_infinite_cone_invariance_certificate:
  exact infinite-index maximum principle; in adjacent-defect coordinates the
  heat ODE is cooperative with nonnegative source and bounded potential.
  Uniform h_k->0 yields a finite active set for every negative spatial
  minimum, and the Dini-minimum argument propagates the full cone from
  lambda=-100 to every finite lambda>=-100, including lambda=0

jensen_window_pf_defect_complete_monotonicity_scout:
  finite Arb diagnostic with exact guard; certifies 3284 defect and 3288
  negative-log-contraction alternating differences, both complete through
  order 8, while retaining 838 high-order interval inconclusives. The exact
  x=(1/2,1,1,...) Hausdorff-defect witness has cubic Jensen discriminant
  -27/16, so these patterns are not an all-shape bridge

jensen_window_pf_heat_flow_jensen_hierarchy_lemma:
  exact shift-degree heat hierarchy; identifies partial_lambda J_(d,n) with a
  shifted-window operator and with a first/second z-derivative operator on
  J_(d+1,n). The one-atom Hausdorff defects d_k=(4/9)(9/25)^(k-1) give the
  exact full-cone cubic boundary x=5/9, y=21/25, z=589/625 with
  (partial_lambda F)/r_0=329728/2109375>0; even complete defect monotonicity
  plus the local heat ODE does not provide an invariant higher-minor cone

jensen_window_pf_rank_two_boundary_family_lemma:
  exact all-degree benchmark; the sequence A_k=a^(k-1)(c+kb)/u factors every
  shifted Jensen window into one simple and one repeated negative root. In
  alpha coordinates, finite pointwise products are multiplier sequences and
  give a discrete canonical-product formula for -log x_k. The -937/3456
  mixture and <-27/125 fractional-power guards reject arbitrary positive
  measure weights; a genuine counting-measure representation remains open

jensen_window_pf_cubic_reciprocal_defect_invariance_lemma:
  exact degree-3 coordinate and conditional flow theorem with finite Arb entry;
  q_k=(1-x_k)^(-1/2) turns cubic hyperbolicity into
  q_(k+1)-q_k<=1. At saturation, the heat field points inward exactly when the
  next increment satisfies the same bound. All 318 lambda=-100 prefix margins
  and 310 nonnegative-grid margins are certified. The later cubic tail-entry
  certificate closes the all-k lambda=-100 tail; forward-uniformity and every
  higher-degree minor cone remain open

jensen_window_pf_cubic_m100_tail_entry_certificate:
  all-k degree-3 entry theorem at lambda=-100; 4,074 compact negative-skewness
  blocks plus an analytic ray prove kappa_3<0 for every t>=318. Two-sided
  adjacent log walls imply d_m>=1/(5m+1), q_(k+1)-q_k<1 for every k>=319,
  and q_(k+1)-q_k->0. Together with the 318 prefix margins this proves every
  shifted cubic Jensen polynomial hyperbolic at lambda=-100. The subsequent
  forward-uniform certificate propagates degree 3 through lambda=0; all higher
  degrees remain open

jensen_window_pf_cubic_forward_uniform_tail_certificate:
  exact all-shift degree-3 forward theorem; the q-increment vector field is
  one-sided cooperative and obeys sqrt(k)*g_k'<=7*r_k at the weighted barrier.
  The entry estimate sqrt(k)*g_k<12 and an explicit coercive supersolution give
  sup_[-100,L] g_k<=(12+7*R_L*(L+100))/sqrt(k) for every finite L. This closes
  the infinite first-crossing handoff and proves all shifted cubics hyperbolic
  at lambda=0; degree 4, PF-infinity, and RH remain open

jensen_window_pf_quartic_boundary_flow_obstruction:
  exact countermodel gate; derives Disc(J_4)=256*x^6*y^2*Q and the heat
  derivative of Q. A rational hyperbolic double-root boundary point satisfies
  all local ratio walls and three strict neighboring cubic inequalities, but
  has Q'/r_1=-13108711376416987159336748097/
  20606742971316325673502124987495<0. Therefore the propagated cubic cone is
  not by itself a quartic invariant; this does not show failure of the zeta
  trajectory and instead requires a new coupled degree-4 condition

jensen_window_pf_quartic_double_root_threshold_lemma:
  exact local quartic boundary theorem; for
  P=(1+a*w)^2(1+b*w)(1+c*w), p=bc, the heat value at the double root has sign
  u-U with U=(-a^2+2*a+p)(3*a^2-5*a+5*p)/(6*p^2). The complete inward
  condition is (3*a^2-4*a+p)*(u-U)<=0. On the triple-root stratum first-order
  viability requires u=U and the first variation retains (1+a*w)^2. A closed
  contraction-coordinate quartic invariant and higher-time tangency remain open

jensen_window_pf_quartic_quintic_polar_contact_lemma:
  exact adjacent-degree theorem; normalized windows obey
  P_d=P_(d+1)-w*P_(d+1)'/(d+1). For a positive-root hyperbolic quintic, a
  double root of the quartic polar derivative forces a quintic triple root.
  In the quartic boundary coordinates P_5(-1/a) is a nonzero prefactor times
  -(u-U), so the heat threshold is exactly quintic triple contact. This points
  to a general adjacent-degree hierarchy but does not close its infinite
  top-down dependence

jensen_window_pf_cofinal_degree_polar_closure_lemma:
  exact all-degree reduction; repeated zero-pole polar descent proves that one
  hyperbolic terminal degree D closes every lower degree at the same shift.
  Therefore an unbounded sequence D_j of hyperbolic terminal degrees is enough
  for all finite degrees, and cofinal terminal sequences at every shift suffice
  for the complete Jensen family. Current Sturm evidence reaches only degrees
  3 through 12 on a finite grid; no cofinal zeta terminal theorem is proved

jensen_window_pf_cofinal_scaling_limit_equivalence_gate:
  exact noncircularity theorem; P_(D,n)(z/D) converges locally uniformly to the
  fixed-shift entire function F_n. Hence an unbounded hyperbolic degree
  subsequence forces F_n into the Laguerre-Polya class, and Jensen's theorem
  gives the converse. At shift zero the cofinal terminal theorem is already an
  endpoint-equivalent statement, while the bounded 1050-row Sturm ladder does
  not enter the limit

jensen_window_pf_polar_heat_collision_cascade_lemma:
  exact collision-cascade theorem; at a multiplicity-m root the adjacent
  degree controls the full low heat jet. If every higher degree is hyperbolic,
  polar lifting raises the common-root multiplicity once per degree. Such a
  cascade forces the entire coefficient function to be exponential times a
  bounded-degree polynomial. Therefore a non-exponential-polynomial LP
  boundary is strict in every fixed Jensen degree, and failure can approach it
  only through degrees tending to infinity

jensen_window_pf_scaled_double_zero_boundary_layer_lemma:
  exact degree-scaled boundary theorem; the first two Jensen corrections are
  -z^2*F''/(2D) and (z^3*F'''/3+z^4*F''''/8)/D^2. Near a nondegenerate double
  zero rho<0, the joint layer is eta^2+8*rho*tau-rho^2. The finite pair has an
  unscaled D^(-3/2) gap, and its D^(-2) collision correction contains
  F'''/F''=3U'/U, the regularized global root external field. The nested
  finite-degree eventual thresholds exhaust Lambda, but the required degree-
  and zero-height-uniform remainder theorem remains open

jensen_window_pf_newman_root_external_field_lemma:
  exact root-field translation; at a squared double zero the D^(-2) Jensen
  correction is governed by E_x=sum_(j!=*)m_j/(x_j^2-x^2), while its signed
  Newman stiffness is the strictly positive sum K_x=sum_y 1/(x-y)^2. The
  squared pair gap satisfies q'=8-4qS and has correction -16*K_x*(t-t_*)^2.
  Two finite LP products realize opposite field signs, so generic LP or
  coefficient positivity cannot replace the open Xi-specific balance theorem

jensen_window_pf_newman_classical_field_balance_gate:
  exact reference-field and compactness gate; the Riemann-von Mangoldt
  continuum field tends to -pi/8, matching the -pi/4 high-positive-time zero
  drift. Published asymptotics make every sufficiently high fixed-time zero
  simple, so a positive boundary collision lies below exp(C/Lambda). Exact
  even-lattice perturbations show that bounded classical-location error still
  allows either unbounded field sign; a lambda-uniform reciprocal-gap estimate
  or compact-region exclusion remains open

jensen_window_pf_newman_local_odd_count_reduction_lemma:
  exact Stieltjes localization theorem; the published uniform zero count
  controls every field contribution outside radius H=log(4c)^2, leaving
  B(c)=-pi/8+S_H+O(1/log c), where S_H is the inverse-square weighted odd
  local counting discrepancy. An exact even backward-heat polynomial has the
  same -pi/8 field and -pi/4 drift while still producing positive square-root
  birth. Thus both a lambda-uniform Xi odd-count theorem and a separate global
  collision obstruction remain necessary

jensen_window_pf_newman_boundary_energy_direction_gate:
  exact collision-energy and directionality gate; every double-zero birth has
  renormalized pair energy asymptotic to 1/(8*(t-t_*)), so finite integrated
  energy down to the boundary would exclude collision. The exact classical-
  field model has positive stiffness and positive cubic gap jet while its
  energy decreases after birth from a nonintegrable trace. The published
  Rodgers-Tao theorem assumes Lambda<0 and starts at Lambda/4, away from the
  boundary, so an Xi-specific boundary-trace estimate remains open

jensen_window_pf_newman_positive_boundary_attainment_lemma:
  exact compactness composition; under Lambda>0 the absolute strip bound and
  uniform positive-time high-zero theorem trap nonreal zeros below the boundary
  in one compact rectangle, forcing a finite real multiple zero of H_Lambda.
  Thus Lambda<=0 is equivalent to simplicity for every t>0. The universal
  multiplicity-m Hermite split has ordered cluster energy
  m(m-1)/(8*(t-Lambda))+O((t-Lambda)^(-1/2)), so every multiplicity has a
  nonintegrable endpoint trace; Xi simplicity or endpoint control remains open

jensen_window_pf_newman_strict_laguerre_correlation_target:
  open RH-equivalent theorem target; the first Laguerre expression is exactly
  the Fourier transform of K_(1,t)(v)=integral phi_t(s+v)phi_t(s-v)s^2 ds.
  Positive-boundary attainment and Wiener imply Lambda<=0 iff these kernels'
  translates are dense in L1 for every 0<t<=1/2. An exact smooth strictly
  log-concave positive-definite convolution has double Fourier zeros, so the
  target needs stronger Xi-specific correlation or total-positive structure

jensen_window_pf_newman_correlation_hierarchy_gaussian_mixture_gate:
  exact hierarchy and route-elimination gate; all Xi correlations evolve by
  partial_t K_n=2v^2 K_n+2K_(n+1), and a multiple boundary root has a rigid
  hierarchy contact, with F_2=3F_1'' at a double root. Complete monotonicity in
  v^2 would force strict Fourier positivity, but the exact double-exponential
  Phi tail makes K_(1,t) decay faster than every Gaussian and therefore rules
  out Gaussian mixtures and direct PF-infinity membership. A separate
  55-digit scout detects the corresponding negative log-convexity minor

jensen_window_pf_newman_theta_summand_spectral_square_gate:
  exact theta/Mellin and finite-block obstruction; Phi=(R''-R)/8 with
  R'(0)=-1/2 gives H_t=1/16+D_t[C_t]/8 for an explicit second-order D_t and
  reduces L_t to one curvature expression in A_t=D_t[C_t]. The bilateral
  shifted-profile transform reconstructs xi itself. Every finite theta
  truncation retains a positive odd endpoint jet A_N and consequently has
  first Laguerre tail -2A_N^2/x^6+O(x^-8)<0. The infinite modular cancellation
  is essential; finite self/cross spectral squares cannot prove the target

jensen_window_pf_newman_gasper_fake_xi_remainder_gate:
  exact Gasper comparison and route-elimination gate; the scaled fake-Xi block
  is a genuine real-zero benchmark. Arb midpoint certificates with explicit
  tails prove that the best scalar absolute-remainder budget at x=25 exceeds
  one at t=0 and t=1/2; 80-digit values are 2.5852... and 2.5692.... Exactly,
  the kernel ratio Phi/Psi is greater than one at the origin and tends to one
  at infinity, so it cannot be a nonconstant positive cosh-mixture. The
  one-block triangle estimate and direct Cardon convolution transfer are
  closed; sign-aware and multi-block routes remain

jensen_window_pf_newman_gasper_residual_two_block_gate:
  exact and interval-certified two-block obstruction; Phi-Psi_9 is strictly
  positive, but the transform of the residual fails the first Laguerre test.
  More generally, tail positivity confines the 9/5 coefficient beta to
  0<=beta<=pi-3/(2pi). Acb Cauchy derivative bounds and Arb convexity cover
  that full interval with L[R_beta](66)<0 or L[R_beta](50)<0. Larger beta has
  a negative residual tail, and the tail-matched multiplier also fails the
  standard universal-factor imaginary-zero test. Signed coupled blocks remain

jensen_window_pf_newman_classical_three_block_residual_gate:
  exact and interval-certified classical three-block obstruction; the Polya
  P2 and de Bruijn real-zero kernels both leave strictly positive residuals,
  but Arb certifies negative first Laguerre values at x=86 and x=52. Tail and
  origin inequalities confine every nonnegative 9/5/1 block with a globally
  nonnegative residual to a rational triangle. Three Acb spectral certificates
  cover all 64,908 closed parameter boxes in that triangle. Gasper's published
  single-shift square does not sign the mixed products of a multi-shift sum;
  signed higher blocks or a new coupled square remain open

jensen_window_pf_newman_signed_universal_factor_residual_gate:
  exact and interval-certified signed universal-factor obstruction; the full
  signed 9/5/1 Polya multiplier condition is equivalent to one quartic having
  all roots in [0,1]. Exact endpoint and critical-point constraints plus the
  Xi origin bound place every globally positive residual candidate in one
  rational rectangle. Acb residual tests at x=86 and x=122 and exact quartic
  and critical discriminants classify 4,094 adaptive leaves grown from 3,416
  base boxes, with maximum depth six and no unresolved box. Independent LP
  residuals are closed for the whole standard signed three-shift multiplier
  cone; higher shifts and genuinely coupled mixed-term squares remain open

jensen_window_pf_laguerre_scale_mixture_gate:
  exact kernel theorem and preservation guard; each unshifted Jensen window is
  a positive scale integral of L_D^(-1/2), and each fixed scale is hyperbolic.
  Exact two-atom and exponential-density countermodels show that positive
  mixing and log-concavity do not preserve this property. Half-integer Gamma
  mixing is nevertheless hyperbolic in every degree by Euler-Jacobi
  factorization, isolating a possible Xi-specific connection theorem

jensen_window_pf_multiplier_counting_measure_target:
  open sufficient bridge target; asks for an exact convergent unit-atomic
  product of elementary multiplier sequences for the normalized zeta
  coefficients, equivalently a discrete representation of -log x_k. Finite
  log-complete-monotonicity is necessary evidence only, exact guards reject
  fractional weights and convex mixtures, and the natural Mellin interpolation
  is blocked by a negative power-sum Hankel determinant. Integer-only product
  equality remains open without an interpolation-uniqueness theorem

jensen_window_pf_mellin_multiplier_power_sum_obstruction:
  interval-certified continuous-interpolation guard; eleven Arb-enclosed Phi
  log moments produce nine candidate power sums and six shifted Hankel tests.
  Three determinants are strictly negative, including the shift-2 size-4 value
  near -2.588644974358276e-19. This rejects the natural Gamma-normalized Mellin
  product, but not a product identity asserted only at integer indices

jensen_window_pf_monotone_contraction_stress:
  finite monotone-contraction stress diagnostic; validates 2875 finite
  Arb-classified zeta-window rows across degrees d=3..12 and the k<=64 cache,
  all satisfying adjacent log-concavity and increasing ratio contractions
  with no zero-containing required positivity gaps

signed_hankel_jensen_dependency_graph:
  dependency hygiene gate connecting finite evidence, countermodel gates, open
  theorem targets, and the `lambda_le_0_goal` node; validates that finite and
  diagnostic nodes only support open targets and have no direct proving edge
  to the not-proved conclusion

target_jensen_window_pf_bridge:
  theorem target that reformulates all-degree/all-shift Jensen hyperbolicity
  as finite PF-infinity of every binomially weighted window
  B^{d,n,0}_j=binom(d,j) A_{n+j}(0); open and not proved

countermodel_gates:
  validates 11 proof-safety examples, including local heat birth, finite
  coefficient-prefix promotion, finite Jensen-window rectangle promotion,
  finite Schur-prefix promotion, finite signed-Hankel grid promotion, finite
  moment-recurrence promotion, and Stieltjes/Hankel-to-Toeplitz promotion traps
```

## Boundary

Passing the ledger checker means the corpus has not confused finite
certificates, diagnostics, or algebraic translations with the missing bridge
theorem. It does not prove any open target.
