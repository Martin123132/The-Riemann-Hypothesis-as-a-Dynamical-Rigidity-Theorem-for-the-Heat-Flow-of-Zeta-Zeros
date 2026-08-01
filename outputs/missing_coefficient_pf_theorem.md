# Missing Coefficient-PF Bridge Theorem

Date: 2026-07-03

Status: theorem target. This is the exact kind of statement that would turn the current finite Toeplitz/PF certificates into a possible proof route. It is not currently proved.

## Object

Let

```text
c_k(lambda) = mu_{2k}(lambda) / (2k)!
mu_{2k}(lambda) = integral_R u^(2k) exp(lambda u^2) Phi(u) du
F_lambda(s) = sum_{k >= 0} c_k(lambda) s^k.
```

The current certified finite evidence concerns upper-triangular Toeplitz minors of:

```text
T_lambda[i,j] = c_{j-i}(lambda), if j >= i,
T_lambda[i,j] = 0, otherwise.
```

## Theorem Target A: Direct Coefficient PF-Infinity

For `lambda = 0`, prove:

```text
all minors of T_0 are nonnegative.
```

Equivalently:

```text
{c_k(0)}_{k >= 0} is a PF-infinity sequence.
```

If true, ASW/Edrei gives:

```text
F_0 is in the restricted Laguerre-Polya class.
```

Reference anchors:

```text
Aissen-Schoenberg-Whitney I: https://doi.org/10.1007/BF02786970
Edrei II: https://doi.org/10.1007/BF02786971
```

This would be a direct coefficient-side route to the RH conclusion in the transformed variable.

## Theorem Target B: Uniform Heat-Interval Coefficient PF-Infinity

For all:

```text
0 <= lambda <= lambda_1
```

or for all:

```text
lambda in [0, epsilon]
```

prove:

```text
{c_k(lambda)} is PF-infinity.
```

This would be stronger than needed at `lambda = 0`, but it may be more natural if the heat deformation gives a monotonic or variation-diminishing structure.

Hard gate:

```text
No proof may assume real-rootedness of F_0 or H_0.
```

## Theorem Target C: Determinantal Integral Formula

For every finite row set:

```text
0 <= r_1 < ... < r_m
```

and column set:

```text
0 <= q_1 < ... < q_m
```

prove either structural zero or:

```text
det(c_{q_b-r_a}(lambda))_{a,b=1}^m
  = integral_{\Omega} W_{R,Q,lambda}(x) d\nu(x),
```

where:

```text
W_{R,Q,lambda}(x) >= 0
```

for the required lambda range.

This would be the cleanest bridge because it upgrades every Toeplitz minor in one stroke.

Current blocker:

```text
No positive integrand W_{R,Q,lambda} is known.
```

## Theorem Target D: Schur-Positivity Specialization

Using the Jacobi-Trudi correspondence, normalize:

```text
d_k(lambda) = c_k(lambda) / c_0(lambda)
```

and prove that the specialization:

```text
h_k -> d_k(lambda)
```

is Schur-positive:

```text
s_{\lambda/\mu}(d_0,d_1,d_2,...) >= 0
```

for all skew shapes `lambda / mu` corresponding to upper-triangular Toeplitz minors.

This is equivalent in spirit to coefficient PF-infinity, but it may expose the right algebraic structure: a positive specialization, a measure on partitions, a determinantal point process, or a production-matrix model.

The exact Toeplitz-minor-to-skew-shape map is recorded in:

```text
work/rh_compute/scripts/toeplitz_jacobi_trudi_map.py
work/rh_compute/scripts/check_toeplitz_jacobi_trudi_map.py
outputs/toeplitz_jacobi_trudi_bridge_note.md
```

The validator checks the `N = 10`, orders `<= 5` grid:

```text
124129 minors
39094 nonzero maps
85035 structural zeros
```

This verifies the algebraic translation used by Target D; it does not prove
that the resulting specialization is Schur-positive.

Current blocker:

```text
No known positive specialization has been identified for the zeta heat-kernel coefficients.
```

Noncircular sufficient targets include:

```text
1. identify d_k = c_k/c_0 as h_k(alpha) for a nonnegative specialization alpha;
2. construct a planar-network model whose path matrix has entries d_{q-r};
3. construct a production-matrix or continued-fraction model that implies
   Schur positivity before any appeal to ASW/Edrei.
```

## Theorem Target E: Multiplier/Closure Factorization

Find a factorization:

```text
c_k(lambda) = a_k(lambda) b_k
```

where:

```text
{a_k(lambda)} is PF-infinity,
{b_k} is a PF-preserving multiplier sequence,
```

or where the corresponding generating function factorizes as a product of restricted Laguerre-Polya functions.

The obvious factor:

```text
b_k = 1 / (2k)!
```

has a restricted Laguerre-Polya generating function:

```text
sum b_k s^k = cosh(sqrt(s)).
```

But this does not help unless the remaining moment sequence has the right PF property.

Current blocker:

```text
Stieltjes/Hankel moment positivity of mu_{2k} does not imply Toeplitz PF-infinity.
```

Executable countermodel:

```text
python work/rh_compute/scripts/countermodel_gate_examples.py
```

The exact Stieltjes multiplier trap uses the positive measure
`10 delta_0 + delta_1 + delta_2`, whose leading Hankel determinants are
`12, 51, 40, 0`. After normalization `c_k = mu_k/(2k)!`, however:

```text
c_0 = 12
c_1 = 3/2
c_2 = 5/24
c_1^2 - c_0 c_2 = -1/4
```

Thus Target E cannot be closed by ordinary moment positivity plus the
factorial denominator. It needs a genuine PF-preserver theorem with verified
hypotheses or an explicit restricted Laguerre-Polya factorization.

## Theorem Target F: Edrei Log-Power-Sum Representation

Normalize:

```text
d_k(lambda) = c_k(lambda) / c_0(lambda)
H_lambda(z) = sum_{k >= 0} d_k(lambda) z^k
log H_lambda(z) = sum_{n >= 1} ell_n(lambda) z^n
q_n(lambda) = n ell_n(lambda)
```

Set

```text
p_n(lambda) = (-1)^(n-1) q_n(lambda)
a_r(lambda) = p_(r+1)(lambda)
              = (-1)^r [z^r] H_lambda'(z)/H_lambda(z).
```

Sokal's logarithmic-derivative criterion applies because `H_lambda` is entire
and `H_lambda(0)=1`. It gives the exact all-order equivalence

```text
H_lambda is Laguerre-Polya type I
  <=> (a_r(lambda))_(r>=0) is a Stieltjes moment sequence.
```

Equivalently, there is a unique measure of the Edrei zero form for which

```text
H_lambda'(z)/H_lambda(z)
  = integral_[0,infinity) dnu_lambda(x)/(1+xz),
nu_lambda
  = gamma*delta_0 + sum_j m_j*beta_j*delta_(beta_j),
gamma>=0, beta_j>0, m_j in Z_(>=1), sum_j m_j*beta_j<infinity.
```

Entireness supplies the discreteness and integer-residue part once the
all-order Stieltjes property is proved. Thus no separate product-reconstruction
conjecture is needed. Sign alternation `p_n>=0` alone remains only necessary.

A clean strict sufficient target is

```text
det(p_(i+j+1)(0))_(i,j=0)^m > 0,
det(p_(i+j+2)(0))_(i,j=0)^m > 0
```

for every `m>=0`. By the Stieltjes moment criterion these are the two required
Hankel columns. The exact indexing and primary-source gate are recorded in
`outputs/jensen_window_pf_edrei_stieltjes_equivalence_gate.md`.

Finite status:

```text
work/rh_compute/scripts/arb_edrei_log_sign_probe.py
work/rh_compute/scripts/check_edrei_log_sign_manifest.py
work/rh_compute/scripts/arb_edrei_power_hankel_probe.py
work/rh_compute/scripts/check_edrei_power_hankel_manifest.py
outputs/edrei_log_sign_diagnostic.md
```

Using rigorous `c_ball` enclosures through `k = 64`, the manifest validates:

```text
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}
n = 1..64
320/320 finite log-sign diagnostics pass
```

The same coefficient enclosures also test finite pieces of the exact
all-order Stieltjes criterion. For:

```text
p_n(lambda) = (-1)^(n-1) q_n(lambda)
```

the finite shifted Hankel/power-sum diagnostic validates:

```text
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}
m = 0,    s = 1..57
m = 1,    s = 1..55
m = 2,    s = 1..53
m = 3,    s = 1..51
m = 4,    s = 1..49
m = 5,    s = 1..47
m = 6,    s = 1..45
m = 7,    s = 1..43
m = 8,    s = 1..41
m = 9,    s = 1..39
m = 10,   s = 1..37
m = 11,   s = 1..35
m = 12,   s = 1..33
m = 13,   s = 1..31
m = 14,   s = 1..29
m = 15,   s = 1..27
m = 16,   s = 1..25
m = 17,   s = 1..23
m = 18,   s = 1..21
m = 19,   s = 1..19
m = 20,   s = 1..17
m = 21,   s = 1..15
m = 22,   s = 1..13
m = 23,   s = 1..11
m = 24,   s = 1..9
m = 25,   s = 1..7
m = 26,   s = 1..5
m = 27,   s = 1..3
m = 28,   s = 1
4,205/4,205 finite Edrei power-Hankel diagnostics pass
```

A larger scout to `s = 48` found no negative determinants but did encounter
high-shift interval-width inconclusives beginning at `m = 2, s = 47` and
moving inward for larger `m`. The top edge is now promoted through
`m = 24, s = 1`, and the lower high-shift wedge is promoted through
`2m+s <= 49` for `m <= 8`, after tighter coefficient enclosures through
`k = 49` and Edrei-log rows through `n = 49`. Higher-shift frontiers
are not promoted.

The next uniform promoted layer uses coefficient enclosures through `k = 57`
and Edrei-log rows through `n = 57`, now completing the five-lambda staircase
`2m+s <= 57`. The `n = 51` layer needed a tighter `lambda = 1e-6` repair; the
`n = 53` and later layers needed tighter `k <= 57`, `dps = 220`,
`abs_tol = 1e-140` coefficient inputs and `dps = 2400` Edrei-log reruns for all
five lambdas.
The checker
`work/rh_compute/scripts/check_edrei_power_hankel_boundary_manifest.py`
validates both the retired blockers and the repaired positive rows.

Current blocker:

```text
finite log-sign alternation, finite power-Hankel positivity, and finite recurrence positivity do not prove that a_r=p_(r+1) is a Stieltjes moment sequence at all orders; no all-m proof of the two required Hankel columns or nonnegative S-fraction is known.
```

The direct ratio handoff is now exact rather than schematic. Krein's
Stieltjes/Pick characterization and Sokal's criterion give

```text
H_0 in LP+
  <=> P_Phi(z):=-Im(F_0'(z)*conj(F_0(z)))>0
      for every Im(z)>0.
```

Polarization writes `P_Phi` as a concrete double integral against
`Phi(u)Phi(v)`. An equivalent sufficient route is a noncircular identity

```text
H_0'(z)/H_0(z)=gamma+<v,(I+z*T)^(-1)v>
```

with `T` positive self-adjoint and derived from Phi. Neither the global sign
nor such an operator is known. The exact target and two-scale mixture guard
are in `outputs/jensen_window_pf_phi_pick_kernel_target.md`.

The same endpoint now has an explicit published arithmetic-operator form.
With

```text
E_(omega,nu)(r)=xi(1/2+omega-ir)^nu
```

and Suzuki's truncated Hankel operator `K_(omega,nu)[t]`, a concrete
equivalent target is:

```text
choose omega_n decreasing to 0 and nu_n*omega_n>1;
prove det(I+/-K_(omega_n,nu_n)[t])!=0 for every t>=0;
```

Suzuki's published criterion also lists a terminal `J`-limit. The internally
audited causal-multiplier reduction shows that this terminal premise is
redundant for a strictly decreasing cofinal determinant family: all finite
forms extend to one bounded causal multiplier, and any off-line zero would
force distinct shifted zeros accumulating at itself. This does not prove the
displayed determinant condition. Suzuki proves only a local interval
unconditionally. The Fredholm expansion writes every coefficient as an
integral of a continuum Hankel determinant, but total nonnegativity alone is
insufficient: a rank-one TN kernel can have eigenvalue one and `det(I-K)=0`.
Horizontal xi-modulus growth and real-boundary unimodularity are also
theorem-mismatched shortcuts. The exact route and guards are in
`outputs/jensen_window_pf_xi_pick_suzuki_hankel_bridge.md` and
`outputs/jensen_window_pf_suzuki_determinant_only_reduction.md`.

The truncations form a Hilbert-Schmidt-continuous self-adjoint path starting
at zero, with nondecreasing operator norm. Therefore all-`t` determinant
nonvanishing is exactly pointwise `||K[t]||<1` for every finite `t`, or
equivalently exclusion of a finite first `+1` or `-1` crossing. This does not
ask for one gap uniform in `t`; such a gap is impossible in the desired
Hermite-Biehler case. The exact Hilbert-Schmidt formula also diverges for
every nonzero kernel, and Suzuki's sup-contour local certificate has a finite
time ceiling. The corrected spectral target and route guards are in
`outputs/jensen_window_pf_suzuki_spectral_frontier.md`.

The Edrei moments also satisfy an exact radial heat hierarchy, but the LP+
model `(1+beta*z)^2` exits every shifted `2x2` Stieltjes wall immediately
under backward heat. Thus generic backward cone invariance is unavailable;
an Xi/Phi-specific all-order no-escape theorem would be required. See
`outputs/jensen_window_pf_edrei_heat_flow_boundary_gate.md`.

Arb recurrence reconstruction scout:

```text
work/rh_compute/scripts/edrei_moment_quadrature_scout.py
outputs/edrei_moment_quadrature_scout.md
```

Using the shifted sequence `a_n = p_{n+1}`, finite Arb recurrence data is
positive for orders `2..12` across the five-lambda grid. A broader order
`2..20` scout becomes interval-inconclusive from order `13` onward, with no
negative rows. This is a finite constructive diagnostic for Target F, not an
all-order representation proof.

Countermodel gate:

```text
python work/rh_compute/scripts/countermodel_gate_examples.py
```

The exact rational moment-recurrence trap preserves the finite recurrence
prefix through order `12` using factorial moments and then chooses a positive
next even moment that makes the next Hankel determinant negative. This blocks
any Target F proof attempt that tries to promote finite recurrence positivity
to the all-order Stieltjes property without a new infinite theorem.

## Numerical Evidence Boundary

The finite certificates are currently evidence for the theorem target, not a proof of it. They can be used to:

```text
1. falsify bad theorem targets;
2. identify likely structural zero patterns;
3. calibrate conjectures;
4. guide which all-order theorem is plausible.
```

They cannot be used to infer the theorem by extrapolation.

## Required Referee-Grade Upgrade

Before this route can be promoted in a manuscript, we need one of:

```text
1. a proof of Target A;
2. a proof of Target B;
3. a proof of Target C;
4. a proof of Target D;
5. a valid closure/factorization proof under Target E;
6. an all-order Stieltjes-moment proof for a_r=p_(r+1), equivalently a
   nonnegative Stieltjes continued fraction or the two all-order Hankel columns
   under Target F;
7. Suzuki's cofinal all-t Fredholm nonvanishing, equivalently pointwise
   strict contractivity of every selected finite truncation;
8. an external theorem whose hypotheses are verified exactly for c_k(lambda).
```

Until then, the honest status is:

```text
certified finite evidence for a coefficient PF-infinity conjecture.
```

The route triage and kill gates are recorded in:

```text
outputs/coefficient_pf_decision_tree.md
```

The focused positive Schur-specialization contract is recorded in:

```text
outputs/positive_schur_specialization_target.md
python work/rh_compute/scripts/check_positive_schur_specialization_target.py
```

The theorem-family fit/misfit audit is recorded in:

```text
outputs/sign_regularity_theorem_fit_matrix.md
```

Its current conclusion is:

```text
ASW/Edrei, Polya-Schur/Jensen, and Sokal's logarithmic-derivative criterion
give equivalent endpoint languages; ordinary Hankel moment positivity of the
raw moments does not imply Toeplitz PF-infinity, while Stieltjes positivity of
the nonlinear Edrei-log sequence a_r=p_(r+1) is endpoint-equivalent.
```

The three noncircular bridge targets now worth prioritizing are:

```text
1. construct a positive Schur/Edrei-Thoma specialization h_k -> d_k(0);
2. prove that a_r=p_(r+1) is a Stieltjes moment sequence, for example through
   the two all-order Hankel columns or a nonnegative S-fraction;
3. prevent every finite +/-1 crossing in Suzuki's arithmetic-Hankel
   truncation path by a signed quadratic-form or determinant factorization
   without contour, Hermite-Biehler, or RH circularity;
4. derive a positive determinant integral formula for every Toeplitz minor.
```
