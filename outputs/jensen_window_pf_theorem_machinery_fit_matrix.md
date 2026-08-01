# Jensen-Window PF Theorem Machinery Fit Matrix

Date: 2026-07-24

Status: theorem-search fit matrix. This is not a proof of Jensen-window
PF-infinity, Jensen hyperbolicity, Laguerre-Polya membership, RH, or
`Lambda <= 0`; it audits which total-positivity and zero-preserver theorem
families could or could not close `jwpf_06_sign_regular_to_jensen_pf_conversion`.

## Purpose

The obligation ledger identifies the central open bridge:

```text
jwpf_06_sign_regular_to_jensen_pf_conversion
```

This matrix asks a sharper question: which known theorem machinery can turn
proved sign-regular information about the actual `A_k(0)` into total
nonnegativity of every binomially weighted Jensen-window Toeplitz matrix?

Machine-readable matrix:

```text
work/rh_compute/results/jensen_window_pf_theorem_machinery_fit_matrix.json
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_theorem_machinery_fit_matrix.py
```

Current result:

```text
validated Jensen-window PF theorem machinery fit matrix: 11 rows, 0 issues, 0 ready-to-apply rows
```

## Required Bridge Features

A theorem capable of closing `jwpf_06` must:

```text
output every binomially weighted Jensen-window Toeplitz matrix for all d,n
handle binomial weights binom(d,j)
handle all shifts n uniformly
start from proved noncircular hypotheses about the actual A_k(0)
avoid assuming Jensen hyperbolicity, Laguerre-Polya membership, RH, or Lambda <= 0
```

## Source-Anchored Rows

```text
tm_01_asw_edrei_pf_sequence_characterization:
  endpoint_equivalence_only
  ASW/Edrei plus Polya-Schur and derivative closure identify all-order
  coefficient PF-infinity with every Jensen window being PF-infinity. They do
  not convert signed-Hankel data into that common endpoint theorem.

tm_08_sokal_log_derivative_stieltjes_criterion:
  endpoint_equivalence_only
  for the entire normalized H, the Stieltjes moment property of
  a_r=p_(r+1)=(-1)^r[z^r]H'/H is exactly equivalent to LP+, coefficient
  PF-infinity, and the Jensen endpoint. It does not transfer the failed
  original A_k signed-Hankel antecedent into that endpoint.

tm_09_krein_stieltjes_pick_phi_kernel:
  possible_if_new_kernel_representation
  the endpoint is exactly P_Phi(z)=-Im(F_0'(z)conj(F_0(z)))>0 throughout the
  upper half-plane, equivalently positivity of a polarized Phi double integral.
  A noncircular Phi-derived positive resolvent would also close the endpoint.

tm_10_edrei_radial_heat_backward_boundary:
  route_mismatch
  the exact Edrei-moment heat hierarchy and finite-rank boundary-flux identity
  identify a canonical high-shift collision sensor. An exact
  orthogonal-polynomial Schur quotient proves that repeated atoms force
  exponential logarithmic determinant growth, whereas simple atoms give only
  polynomial growth. The remaining Phi-specific target is the subexponential
  fixed-size high-shift flux estimate; generic cone invariance is false, and
  independent algebraic-accuracy raw-moment estimates are exponentially
  under-resolved at collision.

tm_11_suzuki_arithmetic_hankel_canonical_system:
  possible_if_new_kernel_representation
  Suzuki constructs the arithmetic Hankel operator unconditionally on a local
  interval and gives an exact global criterion. Along the continuous
  self-adjoint truncation path, the first live hypothesis is equivalently
  pointwise ||K[t]||<1 for every finite t on a cofinal omega sequence. The
  published form is det(I+/-K[t])!=0 for all t. The causal-multiplier audit
  proves that this cofinal determinant family already implies RH by
  shifted-zero accumulation, so the published terminal premise is redundant
  at sequence level. At fixed omega, this is equivalent to reduced innerness
  and exact horizontal zero cancellation. Suzuki's scalar criterion therefore
  leaves an explicit L2 Jordan-totient summatory residual, or the stronger
  eventual-sign condition, as the direct arithmetic target. Suzuki's full
  hierarchy permits eventual sign at any logarithmic smoothing order on one
  cofinal omega sequence. The generalized Muntz identity adds an unsmoothed
  finite target: uniform boundedness of natural Mobius partial sums in their
  exact weighted all-height Mellin-Plancherel norm. An invertible tail-Hardy
  operator identifies this norm with Burnol's weighted natural
  fractional-part approximants at epsilon=2omega. Reciprocal cells make that
  norm equivalent to a discrete energy Q_(omega,N) and force the
  reciprocal-zeta rate N^(-1/2-omega). The exact finite-tail split further
  isolates limiting energy and a summable dyadic Mobius square function. A
  stationary Gram form also isolates every possible growth in one signed
  Mobius off-diagonal. Its Ornstein-Uhlenbeck leading term is exactly a
  reciprocal-Mobius tail energy with a cofinal RH-equivalent uniform bound.
  The sign-indefinite remainder blocks comparison with the full Gram. All
  antecedent estimates remain open. Weighted Hardy/Copson identities make
  the leading energy a dyadic Mertens mean square, equivalently one signed
  origin-anchored cumulative Mobius-correlation target; averaged-Chowla
  does not control that anchored slice. The weighted-prefix transfer shows
  that the first q block misses a constant anchor mode, while pre-collapse
  Vaughan gives the exact signed handoff B_1(X)-TI_X+TII_X. Across all
  cutoffs, adjacent reciprocal tails reconstruct the anchor and identify
  its exact cost as `sum_N P_(alpha/2,N)/N`; a critical scalar model proves
  `sup_N P_(omega,N)` is one logarithm short. Generic Hardy and
  logarithmic-spacing estimates still cost `N^(1-alpha)`. Exact Neumann
  cosine diagonalization proves all dyadic tail modes above
  `ceil(sqrt(K))` summable and confines the RH-equivalent burden to the
  mean and `O(sqrt(K))` smooth modes. Davenport's arbitrary logarithmic
  savings remain power-short. The local-path transference identifies the
  mean as `sqrt(K)r_K-tau_(K,0)` and the nonconstant modes as an anchored
  Mertens-path DCT. Its exact non-Mobius model has zero dyadic endpoint tails
  and zero means but a divergent first mode, so generic coefficient or
  large-sieve machinery cannot close the gate. See
  `outputs/jensen_window_pf_mertens_local_path_cosine_transference.md`.
  Passing to the quotient by constants makes the Abel map uniformly
  invertible and identifies the nonconstant gate with Brownian-bridge
  kernel `min(i,j)-ij/K`; its exact pre-collapse handoff is
  `C_K=-TI_K+TII_K`. The diagonal is summable, but existing averaged or
  separately absolute machinery remains one power short. See
  `outputs/jensen_window_pf_mertens_centered_bridge_vaughan_handoff.md`.
  The invertible scale filter `u_K-(1/2)u_(2K)` then cancels the affine
  mean's infinite plateau and combines its compact tent square with the
  Brownian bridge as one local kernel. Its two-interval handoff is
  `O_(loc,K)=-TI_(loc,K)+TII_(loc,K)`; this is a sharper theorem-search
  coordinate, not a Type I/II estimate. The reciprocal-weighted version
  also has an equivalent positive compact tail-lattice square because the
  affine rank-one term fills the bridge's constant direction exactly. A
  uniformly invertible finite Abel map then gives an equivalent unweighted
  Mobius suffix-square target with one bounded-coefficient future anchor.
  Its sufficient calibration is H_K=O(K^(2+epsilon)), epsilon<alpha,
  versus the trivial O(K^3). See
  `outputs/jensen_window_pf_mertens_affine_tent_bridge_handoff.md`.
  The ordinary suffix energy can furthermore be written as the norm of
  a bounded K-dimensional Mobius feature sum with explicit Gram diagonal
  at most (5/4)K^2. Componentwise finite Vaughan on the two dyadic
  intervals gives B_K=-V_(I,K)+V_(II,K) before squaring. This makes a
  vector-valued Type I/II square-function or variational estimate a
  precise theorem-search target while preserving the signed difference.
  See
  `outputs/jensen_window_pf_mertens_ordinary_vector_vaughan_handoff.md`.
  Diagonalizing the min-kernel gives the mixed-boundary half-odd sine
  modes theta_(K,r)=(2r-1)pi/(2K+1). The target is uniformly equivalent
  to sum_K K^(-alpha)sum_r|X_(K,r)|^2/(2r-1)^2, with each X_(K,r) a
  bounded ordinary-Mobius test and
  X_(K,r)=-T_(I,K,r)+T_(II,K,r) before squaring. This is the sharp
  additive-twist/large-sieve interface; generic Parseval and Davenport
  remain one power short. See
  `outputs/jensen_window_pf_mertens_mixed_boundary_sine_vaughan_handoff.md`.
  Splitting each mode into its current block and endpoint anchor proves
  that the modes above R_(alpha,K)=ceil(K^(1-alpha/2)) have normalized
  tail at most 2K^(-1)+8K^(-alpha/2), hence are dyadically summable.
  For each fixed positive alpha the additive-twist theorem search may
  therefore be restricted exactly to O(K^(1-alpha/2)) low modes. The
   reduction is not uniform at alpha=0 and proves no low-mode gain. See
   `outputs/jensen_window_pf_mertens_spectral_anchor_high_mode_reduction.md`.
   Summing those retained modes gives an exact entrywise-positive
   ordinary-Mobius kernel with
   `P_(K,infinity)(i,j)=pi^2(2K+1)^(-2)min(i,j)`, a uniformly bounded
   diagonal below `7pi^2/24`, and monotone endpoint Abel weights.
   Finite Vaughan preserves the joint cross term `-2<u_I,G u_II>`.
   Thus the theorem-search target is now precisely a weighted signed
   low-kernel off-diagonal or joint bilinear estimate; positivity and
   bounded trace do not prove it. See
   `outputs/jensen_window_pf_mertens_truncated_half_odd_kernel_handoff.md`.
   The off-diagonal now has an exact shift decomposition. Every fixed
   shift has `O(1/K)` height and
   `Var_(n in I_(K,h))W_(K,h)(n)<6rho_K`; every fixed base point has a
   dual `O(1/K)` variation bound. Abel summation exposes maximal
   shifted correlations and local Mertens intervals, but either
   axiswise absolute route remains power-short even under conditional
   square-root input. The theorem search must therefore preserve
    joint cancellation across shifts, base points, Vaughan types, and
    possibly scales. See
    `outputs/jensen_window_pf_mertens_shift_kernel_variation_handoff.md`.
    Zero extension and two-dimensional Abel summation now keep both
    arithmetic coordinates together. The resulting mixed kernel variation is
    below `3*pi^2*(1+log(2R))/K`, and weighted Cauchy-Schwarz isolates the
    sufficient curvature energy
    `E_(alpha,K)=O_epsilon(K^(1+epsilon))`. This estimate remains open; its
    stronger maximal-prefix version already contains RH-scale block-Mertens
    input, while an explicit block-diagonal array blocks promotion from
    separate row/column square-root prefixes. See
    `outputs/jensen_window_pf_mertens_planar_abel_handoff.md`.
    Lifting each triangular-band pair to an edge indicator writes this
    curvature energy as one exact Gram quadratic form. Its complete diagonal is
    below `(27/2)*pi^2*K*(1+log(2R))`, so the sufficient energy target is
    equivalent to the signed edge-pair bound
    `|C_(alpha,K)|=O_epsilon(K^(1+epsilon))`. One-vertex collisions and
    four-distinct-endpoint pairs remain jointly signed, and supported finite
    Vaughan symmetrization retains their cancellation. Lewko-Lewko's
    variational BDH theorem is linear and character-averaged,
    Tao-Teravainen's quantitative Mobius uniformity has iterated-log decay,
    and averaged Chowla is power-short at the natural curvature scale; none
    provides this one-vector quadratic nested-edge estimate. See
    `outputs/jensen_window_pf_mertens_planar_edge_gram_vaughan_handoff.md`.
    The complete endpoint-incidence audit now gives exact formulas
    `C=C_3+C_4` and a perfect cancellation witness
    `C_3=210, C_4=-210, C=0`. Fair random cuts reproduce the strong
    energy gate exactly, so they do not lower its arithmetic cost.
    Before absolute mixed curvature, the current-current interior has
    the genuinely weaker signed projection
    `O_(CC,int,K)=(4*pi^2/N^3)sum_(r<=R)sigma_(K,r)Y_(K,r)`.
    A global modewise flux now closes the transition, current-future,
    and future-future cells exactly:
    `O=O_(CC,int)+sum_(r<=R)Z_(K,r)/q_r^2`. The joined `Z_(K,r)`
    projection has boundary part
    `Z_(K,r)` is the mandatory joined difference
    `phi_r(K-1)U_r-e_r V_r`; bounding its two terms separately discards
    a real shoulder cancellation. Orthogonal completion also gives the
    lossless low-mode energy
    `E_perp+c(beta+gamma/c)^2`. A further suffix-Abel reindexing gives
    `R_B=(E_B-D_B)/2+C_delta`, with `E_B>=0`, `D_B<5pi^2/12`, and
    `|R_B-E_B/2|=O(1+log R)`. Thus the complete boundary/future
    remainder is subpower exactly when `E_B` is subpower. Uniform odd-sine
    control also gives `E_B<=pi(1/2+2pi)K^(-2)H_K` for the earlier
    affine-tent energy. Square-wave bulk localization leaves only an
    `O(K^(alpha/2))` terminal collar, and bounded adjacent suffix increments
    recover it. Hence cofinal all-epsilon `E_B` control is RH-equivalent.
    A bounded terminal witness has `K^2 E_B/H_K->0`, so this is not a
    uniform reverse norm comparison. The separated `Y_(K,r)` plus `E_B`
    route is therefore not a weaker shortcut. The remaining candidate join
    also closes exactly:
    `J_B=O_(CC,int)+E_B/2=L/2+(D_B-D_L)/2-C_delta`.
    Both diagonals are bounded and `C_delta=O(log R)`, so the join returns
    the lossless criterion up to a summable defect. The `Y_(K,r)`, `E_B`,
    and lossless estimates remain open.
    Menon's 2026 almost-all short-interval results
    (`arXiv:2607.15574`) remain logarithmic and have the wrong averaging
    quantifiers for this anchored quadratic target. See
    `outputs/jensen_window_pf_mertens_planar_incidence_antidiagonal_handoff.md`,
    `outputs/jensen_window_pf_mertens_planar_boundary_flux_anchor_handoff.md`,
    `outputs/jensen_window_pf_mertens_planar_boundary_suffix_energy_handoff.md`,
    `outputs/jensen_window_pf_mertens_planar_boundary_suffix_localization_gate.md`,
    and
    `outputs/jensen_window_pf_mertens_planar_joined_energy_equivalence_gate.md`.
    The downstream float64 diagnostic records 28 actual-Mobius rows through
    `K=1024`; on that finite grid `K*V_K<6.571`, `E_K/K<0.831`, and direct
    versus double-Abel discrepancies are below `6.7e-16`. This is
    non-falsification only. See
    `outputs/jensen_window_pf_mertens_planar_curvature_energy_scout.md`.
    A symmetric
  four-zero gate further proves that fixed-shift rational innerness and
  finite causal H2 energy do not generically control reciprocal square
  energy on the same shifted line.

tm_02_schoenberg_pf_functions_variation_diminishing:
  possible_if_new_kernel_representation
  promising only if we construct a PF kernel or variation-diminishing operator
  that produces the actual Jensen-window Toeplitz minors.

tm_03_karlin_basic_composition_cauchy_binet:
  possible_if_new_kernel_representation
  promising only if each Jensen-window Toeplitz minor becomes a Cauchy-Binet
  sum/integral of nonnegative determinant products.

tm_04_polya_schur_multiplier_preservers:
  possible_if_preserver_hypotheses_are_proved
  relevant to binomial weights, but it still needs an input family of windows
  already proved to satisfy the appropriate real-root/PF hypothesis.

tm_05_gantmacher_krein_sign_regular_matrices:
  possible_if_preserver_hypotheses_are_proved
  closest in spirit to signed-Hankel evidence, but no transfer theorem from the
  signed/indefinite Hankel pattern to Jensen-window Toeplitz TN is identified.

tm_06_laguerre_polya_jensen_limit:
  conditional_downstream_only
  useful only after all-degree/all-shift Jensen hyperbolicity is proved.

tm_07_finite_grid_or_rh_assuming_shortcuts:
  rejected_circular
  finite grids, local repulsion, or RH-assuming arguments are proof-safety traps.
```

The exact endpoint equivalence used in `tm_01` is machine-audited in:

```text
outputs/jensen_window_pf_coefficient_pf_equivalence_gate.md
work/rh_compute/results/jensen_window_pf_coefficient_pf_equivalence_gate.json
python work/rh_compute/scripts/check_jensen_window_pf_coefficient_pf_equivalence_gate.py
```

It corrects the earlier route separation while leaving the `jwpf_06`
signed-Hankel/Xi-to-endpoint implication open.

The exact logarithmic-derivative/Stieltjes endpoint equivalence used in
`tm_08` is machine-audited in:

```text
outputs/jensen_window_pf_edrei_stieltjes_equivalence_gate.md
work/rh_compute/results/jensen_window_pf_edrei_stieltjes_equivalence_gate.json
python work/rh_compute/scripts/check_jensen_window_pf_edrei_stieltjes_equivalence_gate.py
```

Its live antecedent is an all-order proof that `a_r=p_(r+1)` is Stieltjes,
equivalently both Hankel columns `s=1,2`, a nonnegative S-fraction, or a direct
Stieltjes representation of `H_0'/H_0`. The finite staircase does not verify
that hypothesis.

The exact Pick-kernel and radial-heat refinements used in `tm_09` and `tm_10`
are machine-audited in:

```text
outputs/jensen_window_pf_phi_pick_kernel_target.md
work/rh_compute/results/jensen_window_pf_phi_pick_kernel_target.json
python work/rh_compute/scripts/check_jensen_window_pf_phi_pick_kernel_target.py
outputs/jensen_window_pf_edrei_heat_flow_boundary_gate.md
work/rh_compute/results/jensen_window_pf_edrei_heat_flow_boundary_gate.json
python work/rh_compute/scripts/check_jensen_window_pf_edrei_heat_flow_boundary_gate.py
outputs/jensen_window_pf_edrei_hankel_boundary_flux_gate.md
work/rh_compute/results/jensen_window_pf_edrei_hankel_boundary_flux_gate.json
python work/rh_compute/scripts/check_jensen_window_pf_edrei_hankel_boundary_flux_gate.py
outputs/jensen_window_pf_edrei_raw_moment_collision_resolution_gate.md
work/rh_compute/results/jensen_window_pf_edrei_raw_moment_collision_resolution_gate.json
python work/rh_compute/scripts/check_jensen_window_pf_edrei_raw_moment_collision_resolution_gate.py
```

They sharpen the surviving Phi sign/operator target, reject generic backward
cone propagation, derive the exact orthogonal-polynomial quotient behind the
collision rate, prove that its raw-moment subtraction condition number grows
like `(beta/q)^s`, and isolate the exact equivalent rigidity criterion

```text
limsup_(s to infinity)
  max(1,partial_lambda log D_(r,s)(lambda))^(1/s) <= 1
```

for every fixed `r` at a candidate positive Newman boundary. Neither proves
this Phi-specific subexponential estimate. With Lemma 11.9 and the collision
sensor, the all-`r` criterion is exactly equivalent to `Lambda<=0`; it is a
rigidity reformulation, not a weaker theorem or a proved Xi/Phi antecedent.
The collision-resolution gate further shows that a symmetric gap `epsilon`
is not resolved until
`s=2*log(beta/epsilon)/log(beta/q)+O(1)`. A viable moment route must preserve
the determinant cancellation structurally or provide gap-adapted
exponentially small remainders. Adjoining any fixed larger simple prefix
preserves this crossover at arbitrary determinant rank.

The exact xi directional coordinate and Suzuki arithmetic-Hankel refinement
used in `tm_11` are machine-audited in:

```text
outputs/jensen_window_pf_xi_pick_suzuki_hankel_bridge.md
work/rh_compute/results/jensen_window_pf_xi_pick_suzuki_hankel_bridge.json
python work/rh_compute/scripts/check_jensen_window_pf_xi_pick_suzuki_hankel_bridge.py
outputs/jensen_window_pf_suzuki_spectral_frontier.md
work/rh_compute/results/jensen_window_pf_suzuki_spectral_frontier.json
python work/rh_compute/scripts/check_jensen_window_pf_suzuki_spectral_frontier.py
outputs/jensen_window_pf_suzuki_determinant_only_reduction.md
work/rh_compute/results/jensen_window_pf_suzuki_determinant_only_reduction.json
python work/rh_compute/scripts/check_jensen_window_pf_suzuki_determinant_only_reduction.py
outputs/jensen_window_pf_suzuki_fixed_omega_phase_diagram.md
work/rh_compute/results/jensen_window_pf_suzuki_fixed_omega_phase_diagram.json
python work/rh_compute/scripts/check_jensen_window_pf_suzuki_fixed_omega_phase_diagram.py
outputs/jensen_window_pf_suzuki_jordan_totient_sign_scout.md
work/rh_compute/results/jensen_window_pf_suzuki_jordan_totient_sign_scout.json
python work/rh_compute/scripts/check_jensen_window_pf_suzuki_jordan_totient_sign_scout.py
outputs/jensen_window_pf_suzuki_cofinal_monotonicity_hierarchy.md
work/rh_compute/results/jensen_window_pf_suzuki_cofinal_monotonicity_hierarchy.json
python work/rh_compute/scripts/check_jensen_window_pf_suzuki_cofinal_monotonicity_hierarchy.py
outputs/jensen_window_pf_suzuki_jordan_error_kernel_reduction.md
work/rh_compute/results/jensen_window_pf_suzuki_jordan_error_kernel_reduction.json
python work/rh_compute/scripts/check_jensen_window_pf_suzuki_jordan_error_kernel_reduction.py
outputs/jensen_window_pf_suzuki_cofinal_l2_hierarchy.md
work/rh_compute/results/jensen_window_pf_suzuki_cofinal_l2_hierarchy.json
python work/rh_compute/scripts/check_jensen_window_pf_suzuki_cofinal_l2_hierarchy.py
outputs/jensen_window_pf_jordan_muntz_causal_energy_bridge.md
work/rh_compute/results/jensen_window_pf_jordan_muntz_causal_energy_bridge.json
python work/rh_compute/scripts/check_jensen_window_pf_jordan_muntz_causal_energy_bridge.py
outputs/jensen_window_pf_jordan_muntz_burnol_hardy_intertwiner.md
work/rh_compute/results/jensen_window_pf_jordan_muntz_burnol_hardy_intertwiner.json
python work/rh_compute/scripts/check_jensen_window_pf_jordan_muntz_burnol_hardy_intertwiner.py
outputs/jensen_window_pf_burnol_cell_energy_tail_obstruction.md
work/rh_compute/results/jensen_window_pf_burnol_cell_energy_tail_obstruction.json
python work/rh_compute/scripts/check_jensen_window_pf_burnol_cell_energy_tail_obstruction.py
outputs/jensen_window_pf_burnol_tail_discrepancy_dyadic_reduction.md
work/rh_compute/results/jensen_window_pf_burnol_tail_discrepancy_dyadic_reduction.json
python work/rh_compute/scripts/check_jensen_window_pf_burnol_tail_discrepancy_dyadic_reduction.py
outputs/jensen_window_pf_weighted_fractional_autocorrelation_gram_bridge.md
work/rh_compute/results/jensen_window_pf_weighted_fractional_autocorrelation_gram_bridge.json
python work/rh_compute/scripts/check_jensen_window_pf_weighted_fractional_autocorrelation_gram_bridge.py
outputs/jensen_window_pf_weighted_autocorrelation_ou_tail_energy_reduction.md
work/rh_compute/results/jensen_window_pf_weighted_autocorrelation_ou_tail_energy_reduction.json
python work/rh_compute/scripts/check_jensen_window_pf_weighted_autocorrelation_ou_tail_energy_reduction.py
outputs/jensen_window_pf_ou_mertens_mean_square_reduction.md
work/rh_compute/results/jensen_window_pf_ou_mertens_mean_square_reduction.json
python work/rh_compute/scripts/check_jensen_window_pf_ou_mertens_mean_square_reduction.py
outputs/jensen_window_pf_mertens_weighted_prefix_affine_defect_reduction.md
work/rh_compute/results/jensen_window_pf_mertens_weighted_prefix_affine_defect_reduction.json
python work/rh_compute/scripts/check_jensen_window_pf_mertens_weighted_prefix_affine_defect_reduction.py
outputs/jensen_window_pf_mertens_anchor_logarithmic_tail_energy_reduction.md
work/rh_compute/results/jensen_window_pf_mertens_anchor_logarithmic_tail_energy_reduction.json
python work/rh_compute/scripts/check_jensen_window_pf_mertens_anchor_logarithmic_tail_energy_reduction.py
outputs/jensen_window_pf_mertens_dyadic_cosine_mode_bottleneck.md
work/rh_compute/results/jensen_window_pf_mertens_dyadic_cosine_mode_bottleneck.json
python work/rh_compute/scripts/check_jensen_window_pf_mertens_dyadic_cosine_mode_bottleneck.py
outputs/jensen_window_pf_fixed_shift_inner_reciprocal_boundary_separation_gate.md
work/rh_compute/results/jensen_window_pf_fixed_shift_inner_reciprocal_boundary_separation_gate.json
python work/rh_compute/scripts/check_jensen_window_pf_fixed_shift_inner_reciprocal_boundary_separation_gate.py
```

The Fredholm series makes the signed-Hankel contact literal, but a rank-one
totally nonnegative kernel can still have `det(I-K)=0`. Likewise,
real-boundary unimodularity does not permit a high-contour inverse Fourier
integral to be shifted through possible poles. The exact spectral-frontier
audit further proves that Hilbert-Schmidt domination diverges, Suzuki's
`M_v^2 exp(4vt)<1` certificate has a finite time ceiling, and a contraction
gap uniform in `t` is impossible in the desired HB case. The route therefore
needs either a cancellation-sensitive signed quadratic-form barrier
preventing a finite `+1` or `-1` crossing, or a direct arithmetic proof of
Suzuki's `L2` summatory residual. Eventual one-sign is a stronger sufficient
target. The 6000-point positive sign scout supplies regression anchors, not
tail control. The cofinal hierarchy makes `k=2` a legitimate first target:
prove eventual sign for a logarithmic antiderivative on one explicit
sequence `omega_j->0`.

The Abel reduction shows why ordinary asymptotics do not suffice: the full
positive residue main term cancels against the Suzuki kernel. What remains
is a signed Mobius-error convolution, and absolute estimates retain an
`x^(1/2-omega)` power loss at every finite smoothing order.

Central Taylor subtraction supplies an exact all-order energy formulation:
`H_(omega,k)(exp(t))-P_(omega,k)(t)` is in `L2(0,infinity)` if and only if
the fixed-shift quotient is inner. Thus one such estimate along any
`omega_j->0` is RH-equivalent. A boundary-unimodular quotient with a
right-half-plane pole gives an exponentially growing residual, so trace
modulus alone does not prove the energy estimate. In fact, the regularized
boundary trace already has finite `L2` norm; the pole appears as an
anti-causal inverse-Fourier component. Any Plancherel route must therefore
prove positive-time support, not merely finite boundary energy. The analogous
Riesz-smoothed classical totient theorem assumes RH and reciprocal
zeta-derivative bounds and does not transfer the missing antecedent.

The unsmoothed generalized Muntz identity now gives the cleanest finite
version. If `E_(omega,N)` is the natural Mobius partial sum, then

```text
||E_(omega,N)||_H^2
 =1/(2*pi) integral_R W_omega(t)
  |zeta(1/2+omega+i*t)M_N(1/2+omega+i*t)|^2 dt.
```

Pointwise convergence and Fatou show that
`sup_N ||E_(omega,N)||_H<infinity` supplies the required error energy.
Burnol and Balazard-Saias machinery gives the converse convergence under RH,
so the cofinal uniform bound is an internally audited RH-equivalent theorem
candidate. The finite Brownian max kernel is totally nonnegative but gives
only a lower bound of zero. Scalar `l2` coefficients do not make the dilation
orbit Bessel, and the audited general Nyman sufficiency does not transfer
because `R_omega` lacks its compact-support hypothesis.

The exact tail-Hardy identity sharpens that literature contact:

```text
E_(omega,N)(1/t)
 =-(I-omega*H_star)[t^(-omega)f_(2omega,N)(t)],

(1-2omega)||t^(-omega)f_(2omega,N)||_2
 <=||E_(omega,N)||_H
 <=(1+2omega)||t^(-omega)f_(2omega,N)||_2.
```

Thus the Jordan and Burnol finite norm problems are equivalent. Burnol's
published cofinal theorem and Balazard-Saias conditional convergence
source-back the RH-equivalent composition, but neither supplies the
unconditional uniform bound.

The reciprocal-cell reduction identifies the exact discrete form

```text
Q_(omega,N)
 =sum_(k>=1)k^(2omega-2)
  |sum_(d<=N)mu(d)d^(-2omega){k/d}|^2.
```

Uniform boundedness of `Q_(omega,N)` is equivalent to the Burnol norm.
The stable prefix `k<=N` forces

```text
|A_(omega,N)-1/zeta(1+2omega)|=O(N^(-1/2-omega)),
```

so Abel summation makes `Re(s)>1/2+omega` zero-free. Any successful
estimate must therefore recover this scalar rate as well as supply the
weighted short-multiplicative-interval cancellation. Neither is proved.

The omitted-divisor identity sharpens the last phrase. Uniform
`Q_(omega,N)` is exactly equivalent to

```text
Q_(omega,infinity)<infinity,
r_(omega,N)=O(N^(-1/2-omega)),
sup_N sum_(j>=0)(2^j*N)^(2omega-2)
 V_(omega,N)(2^j*N)<infinity.
```

The quotient partition `q=floor(k/d)` expresses `V_(omega,N)` through
weighted Mobius increments on `(k/(q+1),k/q]`. A critical block estimate
without a summable gain across `j` is not enough.

The weighted fractional-part autocorrelation gives a second exact view:

```text
Gamma_(alpha,N)
 =sum_(d,e<=N)mu(d)mu(e)(de)^(-(1+alpha)/2)
  C_alpha(log(d/e)).
```

The kernel is positive definite with zeta spectral density, and its diagonal
is uniformly bounded. Positive definiteness still supplies no upper bound
for the signed Mobius off-diagonal; kernel positivity and a bounded Gram diagonal
are therefore nonpromotion guards, not the missing theorem.

## Primary Source Anchors

The matrix records source anchors for the theorem families, including:

```text
https://www.pnas.org/doi/10.1073/pnas.37.5.303
https://link.springer.com/article/10.1007/BF02786970
https://www.cambridge.org/core/journals/canadian-journal-of-mathematics/article/proof-of-a-conjecture-of-schoenberg-on-the-generating-function-of-a-totally-positive-sequence/7CC0A2ADBCB2ADF112F737EDE229766C
https://link.springer.com/article/10.1007/BF02790092
https://annals.math.princeton.edu/wp-content/uploads/annals-v170-n1-p14-p.pdf
https://doi.org/10.1016/j.jmaa.2022.126432
https://web.math.ku.dk/~berg/manus/castellon.pdf
https://doi.org/10.1016/j.jfa.2021.109116
```

These are used as theorem-search anchors, not as claims that their hypotheses
are already satisfied by the zeta heat-flow coefficients.

## Current Best Route

The best current routes are still structural:

```text
construct a positive kernel, planar network, production matrix, determinant
integral, or sign-regular-to-Toeplitz transfer theorem that produces every
Jensen-window Toeplitz minor for B^{d,n,0}_j = binom(d,j) A_{n+j}(0)

or prove the polarized global sign P_Phi(z)>0 for every upper-half-plane z,
construct a Phi-derived positive resolvent for H_0'/H_0, or equivalently prove
both Edrei-log Hankel columns s=1,2 at every order

or close Suzuki's arithmetic-Hankel criterion by preventing every finite
first +/-1 crossing on a cofinal omega sequence through a noncircular signed
quadratic-form or determinant factorization; the determinant-only
causal-multiplier reduction removes the terminal premise at cofinal-sequence
level; equivalently prove Suzuki's explicit L2 Jordan-totient summatory
residual, or the stronger eventual-sign condition, on a sequence omega->0
at any selected logarithmic smoothing orders k_j>=1
through cancellation or order in the exact Jordan-error kernel, not
positive main-term domination; in the sharp finite formulation, prove
sup_N ||E_(omega,N)||_H<infinity by a uniform weighted all-height
natural-mollifier estimate, equivalently prove
sup_N ||t^(-omega)f_(2omega,N)||_2<infinity in Burnol coordinates,
equivalently sup_N Q_(omega,N)<infinity in reciprocal cells
```

The next proof-search action is therefore not another finite grid by itself.
For the Suzuki route, attack the mean and modes
`1<=r<ceil(sqrt(K))` in the exact dyadic reciprocal-tail cosine
decomposition. Their explicit one-hump and smooth sine weights require a
power gain beyond Davenport, not another arbitrary-log estimate. In
parallel, decompose the exact pre-collapse Vaughan forms and seek

```text
-TI_X+TII_X=O_epsilon(X^(2+epsilon))
```

while preserving signed cancellation across both shift and terminal length.
Applying Vaughan after collapse is endogenous and recycles the Mertens
energy. Averaged-Chowla and almost-all interval estimates do not control the
anchored slice. Continue blocking the sign-indefinite remainder by
logarithmic separation and require a summable gain over the critical
post-prefix scale. Limiting energy and the forced `N^(-1/2-omega)`
reciprocal-zeta rate remain separate required inputs.
For the direct PF route it remains a symbolic theorem search for an identity
that survives the degree-3 countermodel and handles the binomial weights.

## Structural Ansatz Workbench

The next layer of this audit is:

```text
outputs/jensen_window_pf_structural_ansatz_matrix.md
work/rh_compute/results/jensen_window_pf_structural_ansatz_matrix.json
python work/rh_compute/scripts/check_jensen_window_pf_structural_ansatz_matrix.py
```

It validates:

```text
validated Jensen-window PF structural ansatz matrix: 6 ansatz rows, 0 issues, 0 ready-to-apply rows
```

The workbench tests positive Cauchy-Binet, planar-network/production-matrix,
determinant-integral, binomial-preserver, direct-Hankel, and finite-grid
ansatz rows against exact degree-2/3/4 low-degree formulas and the finite
countermodel kill gate.

The structural rows now also have a bounded Schur/Jacobi-Trudi shape contract:

```text
outputs/jensen_window_pf_schur_shape_contract.md
work/rh_compute/results/jensen_window_pf_schur_shape_contract.json
python work/rh_compute/scripts/check_jensen_window_pf_schur_shape_contract.py
```

It records the finite-band shape obligations that any positive Schur, network,
production-matrix, Cauchy-Binet, or determinant-integral theorem would need to
cover.

The Cauchy-Binet row also has a low-degree symbolic scout:

```text
outputs/jensen_window_pf_cauchy_binet_low_degree_scout.md
work/rh_compute/results/jensen_window_pf_cauchy_binet_low_degree_scout.json
python work/rh_compute/scripts/check_jensen_window_pf_cauchy_binet_low_degree_scout.py
```

It records `15` formula rows with nonnegative Bernstein coefficients under
adjacent log-concavity, but `0` kernel identities found.

## Boundary

Passing this checker means the theorem-search map does not overstate any known
machinery. It records no `ready_to_apply` row, and it keeps endpoint
equivalences, possible structural routes, downstream limiting arguments, and
rejected shortcuts separated.
