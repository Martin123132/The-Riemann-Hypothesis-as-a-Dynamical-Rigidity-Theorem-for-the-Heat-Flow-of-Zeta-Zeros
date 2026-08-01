# Countermodel Library

Date: 2026-07-09

Status: proof-safety artifact. This is not evidence against RH. It records small models, finite-prefix traps, and finite-grid traps that block invalid bridge lemmas in the RH dynamical-rigidity programme.

Executable gate:

```text
python work/rh_compute/scripts/countermodel_gate_examples.py
```

Current result:

```text
validated 11 countermodel gate examples
```

## Purpose

The programme is trying to prove the missing Newman direction:

```text
Lambda <= 0.
```

The dangerous failure mode is to turn one of these into a proof:

```text
local zero repulsion
finite signed-Hankel evidence
finite Jensen hyperbolicity
finite Jensen-window PF/Sturm rectangles
finite Toeplitz/PF certificates
finite Schur/Toeplitz shape-prefix evidence
finite Edrei moment-recurrence evidence
Stieltjes/Hankel moment positivity without Toeplitz total positivity
```

Each item is useful, but none is logically enough by itself. This file gives concrete gates that proposed proof steps must pass before they can enter a manuscript as more than evidence or a conditional claim.

## Gate 1: Local Heat Birth

The exact model:

```text
F_tau(t) = t^2 - 2 tau
```

satisfies the same Newman heat sign:

```text
partial_tau F = -2
partial_tt F = 2
partial_tau F = - partial_tt F
```

Its zeros are:

```text
tau < 0:  t = +/- i sqrt(-2 tau)
tau = 0:  double real zero at t = 0
tau > 0:  t = +/- sqrt(2 tau)
```

On the real-zero side, the gap is:

```text
g(tau) = 2 sqrt(2 tau)
```

and:

```text
g'(tau) = 4/g(tau).
```

Blocked proof step:

```text
The local law g' = 4/g excludes positive Newman birth.
```

Correct use:

```text
The local law describes the already-real side after a square-root birth.
```

Therefore any proof of `Lambda <= 0` must use global Xi structure, a sign-regularity theorem, or another nonlocal invariant. Local repulsion alone cannot do the job.

## Gate 2: Finite Prefix Is Not All-Order PF

Let:

```text
c_k(lambda) = A_k(lambda) / k!
```

be the ordinary coefficient sequence used in the coefficient-PF route. Suppose a finite prefix through `c_K` has passed every certified Toeplitz/PF test we have run.

That still cannot prove PF-infinity. Preserve the whole prefix through `K`, and define one positive next coefficient by:

```text
c_{K+1} = 2 c_K^2 / c_{K-1}.
```

Then the order-2 Toeplitz minor:

```text
det [[c_K,   c_{K+1}],
     [c_{K-1}, c_K]]
= c_K^2 - c_{K-1} c_{K+1}
= -c_K^2
< 0.
```

This does not say the actual zeta coefficient sequence fails. It says a proof step using only a finite prefix is invalid unless it adds a structural all-order theorem.

Blocked proof step:

```text
The certified finite Toeplitz/PF ledger makes c_k PF-infinity.
```

Correct use:

```text
The finite Toeplitz/PF ledger is falsification pressure and theorem-search evidence.
```

## Gate 2A: Finite Schur Prefix Is Not A Positive Specialization

The positive Schur-specialization target studies:

```text
h_k -> d_k(0)
```

and tries to prove all skew-Schur evaluations nonnegative. A finite set of
Schur or Toeplitz checks still cannot prove such a specialization exists.

The exact gate starts with:

```text
h_k = 1/k!
```

whose generating function is `exp(z)`, a clean restricted PF-infinity model.
It preserves:

```text
h_0, h_1, ..., h_6
```

and exactly checks a finite Toeplitz/Schur grid:

```text
N = 7
orders <= 4
2,940 finite tests
1,204 structurally nonzero positive minors
1,736 structural zeros
```

Then it chooses the next positive complete-homogeneous coordinate:

```text
h_7 = 1/2160
```

At the first untested Jacobi-Trudi shape:

```text
lambda = (6,6)
mu = (0,0)
```

the determinant becomes:

```text
s_(6,6) = det [[h_6, h_7],
               [h_5, h_6]]
        = -1/518400 < 0.
```

Blocked proof step:

```text
A finite Schur/Toeplitz shape ledger proves h_k -> d_k is a positive
specialization.
```

Correct use:

```text
Finite Schur/Toeplitz checks test proposed formulas. A proof needs an
all-order positive specialization, planar network, production matrix,
continued fraction, positive determinant integral, or equivalent theorem.
```

## Gate 3: Finite Signed-Hankel Is Not All-Order Signed Regularity

For:

```text
D_{m,s} = det(A_{i+j+s})_{i,j=0}^m
sigma(m) = (-1)^(m(m+1)/2)
```

the observed signed-Hankel condition is:

```text
sigma(m) D_{m,s} > 0.
```

For `m = 1`, this says:

```text
-(A_s A_{s+2} - A_{s+1}^2) > 0.
```

Given any positive prefix through `A_K`, preserve it and define:

```text
A_{K+1} = 2 A_K^2 / A_{K-1}.
```

At shift `s = K-1`:

```text
D_{1,K-1} = A_{K-1} A_{K+1} - A_K^2
          = A_K^2
          > 0
```

so:

```text
sigma(1) D_{1,K-1} = -A_K^2 < 0.
```

Blocked proof step:

```text
The finite signed-Hankel certificate grid proves all-order signed regularity.
```

Correct use:

```text
The finite signed-Hankel grid supports a conjectural all-order signed-regularity target.
```

## Gate 4: Finite Jensen Hyperbolicity Is Not Jensen Criterion

For degree 2 and shift `K-1`, the Jensen polynomial has the form:

```text
P_{2,K-1}(x)
= A_{K-1} + 2 A_K x + A_{K+1} x^2.
```

With the same positive one-term extension:

```text
A_{K+1} = 2 A_K^2 / A_{K-1}
```

the discriminant becomes:

```text
4(A_K^2 - A_{K-1} A_{K+1})
= -4 A_K^2
< 0.
```

So the next degree-2 Jensen polynomial is not hyperbolic, despite every earlier finite check being preserved.

Blocked proof step:

```text
Many finite Jensen/Sturm passes imply all Jensen polynomials are hyperbolic.
```

Correct use:

```text
Finite Jensen/Sturm passes guide theorem search and catch numerical failures.
```

## Gate 4A: Finite Consecutive Signed-Hankel Grid Is Not All Shifts

The executable gate also includes an exact finite-grid trap independent of
the zeta coefficients. Start with:

```text
a_k = 1/k!
```

For the shifted-principal signed-Hankel determinants:

```text
sigma(m) det(a_{i+j+s})_{i,j=0}^m
sigma(m) = (-1)^(m(m+1)/2)
```

the exact rational check validates the whole grid:

```text
m = 0..4
s = 0..8
45/45 signed determinants > 0
coefficients used: a_0..a_16
minimum signed grid value about 4.048066611965355946e-53
```

Now preserve all those coefficients and choose the next positive coefficient:

```text
a_17 = 1/167382319104000.
```

At the next untested shift:

```text
s = 15
```

the `m = 1` signed-Hankel value becomes negative:

```text
-2.284340357080483672e-27
```

and the degree-2 Jensen discriminant at the same shift is negative:

```text
-9.137361428321934688e-27.
```

Blocked proof step:

```text
The finite shifted-principal signed-Hankel grid proves all shifts, all-order
sign-regularity, or Jensen hyperbolicity.
```

Correct use:

```text
The finite grid is a certified finite diagnostic. A promotion needs a known
Hankel sign-consistency reduction plus an all-order proof, or a new bridge
theorem matching the zeta coefficient sequence.
```

## Gate 4B: Current Jensen-Window Rectangle Is Not All Shifts

The current Arb Jensen-window diagnostics are deliberately finite:

```text
PF obligation checks:
  degrees d = 3,4
  shifts n = 0..20
  coefficients used through A_24

Sturm/root-count checks:
  degrees d = 3,4,5
  shifts n = 0..20
  coefficients used through A_25
```

The executable gate preserves every coefficient those finite Jensen-window
checks can see:

```text
A_0, A_1, ..., A_25
```

for all five lambda rows, then chooses a positive next coefficient:

```text
A_26 = 2 A_25^2 / A_24.
```

At the next untested degree-2 Jensen window, shift `24`, the discriminant is:

```text
4(A_25^2 - A_24 A_26) = -4 A_25^2 < 0.
```

So all existing finite Jensen-window PF/Sturm inputs are preserved, but the
next degree-2 Jensen polynomial is not hyperbolic.

Blocked proof step:

```text
The finite Jensen-window PF/Sturm rectangle proves all-shift Jensen
hyperbolicity.
```

Correct use:

```text
The finite Jensen-window rectangle is a strong stress test and normalization
check. A proof still needs an all-degree/all-shift theorem for the actual
coefficient sequence.
```

## Gate 5: Finite Moment Recurrence Is Not An Edrei Representation

The Edrei reconstruction scout checks finite recurrence data for the shifted
moment sequence:

```text
a_n = p_{n+1}.
```

Positive recurrence data through any fixed order still does not prove an
all-order positive measure or Edrei zero-parameter representation.

The executable gate uses exact rational arithmetic. Start with the factorial
moments:

```text
m_n = n!
```

These are genuine Stieltjes moments of `exp(-x) dx` on `[0, infinity)`, so the
finite recurrence/Hankel prefix is honestly positive. For the current default
order `12`, preserve:

```text
m_0, m_1, ..., m_23
```

Then choose the next even moment `m_24` as a positive number below the exact
Schur-complement threshold for the next Hankel matrix. The gate reports:

```text
all preserved leading Hankel determinants > 0
adversarial m_24 > 0
next Hankel determinant < 0
```

Blocked proof step:

```text
The finite Arb recurrence scout through order 12 proves the all-order Edrei
moment representation.
```

Correct use:

```text
The recurrence scout is a constructive finite diagnostic and precision
frontier. It becomes a proof only after an all-order moment theorem,
positive parameter construction, or analytic recurrence formula is supplied.
```

## Gate 6: Stieltjes Moment Positivity Is Not Coefficient PF

The coefficient-PF route uses:

```text
c_k = mu_k / (2k)!.
```

A tempting but invalid bridge is:

```text
mu_k is a Stieltjes moment sequence
and 1/(2k)! has a restricted Laguerre-Polya generating function
therefore c_k is PF-infinity.
```

The executable gate blocks this with an exact positive measure:

```text
10 delta_0 + delta_1 + delta_2 on [0, infinity).
```

Its moments begin:

```text
mu_0..mu_6 = 12, 3, 5, 9, 17, 33, 65
```

and the leading Hankel determinants are:

```text
size 1..4 = 12, 51, 40, 0
```

so the moment sequence is Stieltjes/Hankel-nonnegative. But after the
coefficient-route normalization:

```text
c_0 = 12
c_1 = 3/2
c_2 = 5/24
```

the first order-2 Toeplitz/PF minor is:

```text
c_1^2 - c_0 c_2 = -1/4 < 0.
```

Blocked proof step:

```text
Stieltjes/Hankel positivity of the moments mu_k, together with the
factor 1/(2k)!, proves coefficient PF-infinity.
```

Correct use:

```text
Moment positivity can motivate the route, but a valid proof needs a
Toeplitz-total-positivity theorem, a positive determinant integral formula,
or an explicit restricted Laguerre-Polya factorization.
```

## Gate 7: Arithmetic-Hankel Shortcuts Do Not Close Suzuki's Criterion

The exact Xi/Pick/Suzuki audit adds three independent proof-safety witnesses:

```text
F_*(z)=z^2+6z+25:
  horizontal growth of M_*(w)=F_*(w^2) holds on Re(w)>1,
  but the required hyperbolic Pick direction is negative at an exact point;

K(x+y)=1 on L2(0,1):
  the Hankel kernel is totally nonnegative and rank one,
  but det(I-K)=0 because its nonzero eigenvalue is one;

Theta_*(r)=(r+i)/(r-i):
  |Theta_*(u)|=1 for every real u,
  but the upper-half-plane pole contributes a residue under contour motion.
```

These block promotion from horizontal xi-modulus monotonicity, Hankel-kernel
total nonnegativity, or real-boundary unimodularity into Suzuki's all-`t`
Fredholm gate. The later determinant-only audit shows that the terminal
premise is redundant for a cofinal all-time family, but these witnesses still
leave that surviving spectral gate open. The executable exact audit is:

```text
outputs/jensen_window_pf_xi_pick_suzuki_hankel_bridge.md
python work/rh_compute/scripts/check_jensen_window_pf_xi_pick_suzuki_hankel_bridge.py
```

## Gate 8: Suzuki Norm Certificates Must Respect Spectral Flow

The exact truncation-path audit adds four route guards:

```text
isolated-time determinant guard:
  A=2P has ||A||=2 while det(I-A)=-1 and det(I+A)=3;
  determinant nonvanishing implies contraction only when the whole
  continuous path from t=0 is included;

Hilbert-Schmidt guard:
  ||K[t]||_HS^2=integral_0^(2t)(2t-s)|K(s)|^2 ds,
  which diverges for every nonzero supported continuous kernel;

high-contour guard:
  Suzuki's sufficient test M_v^2 exp(4vt)<1 has a finite t ceiling for
  every fixed zeta pair (omega,nu), because M_v has polynomial scale;

uniform-gap guard:
  under the desired Hermite-Biehler hypothesis, ||K[t]||<1 for every
  finite t but ||K[t]|| tends to one.
```

These reject an isolated determinant check, a global Hilbert-Schmidt
estimate, indefinite repetition of the local sup-contour argument, and one
epsilon-sized gap uniform in `t`. The surviving target is to prevent a
finite first `+1` or `-1` crossing by a cancellation-sensitive signed
quadratic form. The executable audit is:

```text
outputs/jensen_window_pf_suzuki_spectral_frontier.md
python work/rh_compute/scripts/check_jensen_window_pf_suzuki_spectral_frontier.py
```

## Gate 9: One Fixed Suzuki Pair Cannot Replace Cofinality

All-time determinant nonvanishing for one fixed `(omega,nu)` extends the
supported kernel to a bounded causal multiplier and removes upper-half-plane
poles of the xi quotient. At a denominator zero `rho`, however, this yields
only the cancellation rule

```text
xi(rho-2*omega)=0.
```

One such shifted zero is not a contradiction. A strictly decreasing
`omega_n->0` sequence is essential: it produces distinct zeros
`rho-2*omega_n` accumulating at `rho`, which an entire nonzero xi function
cannot have. Thus neither one fixed family nor a finite omega list may be
promoted to the determinant-only RH equivalence. The exact reduction and
guard are:

```text
outputs/jensen_window_pf_suzuki_determinant_only_reduction.md
python work/rh_compute/scripts/check_jensen_window_pf_suzuki_determinant_only_reduction.py
```

## Gate 10: Fixed-Shift Cancellation And Finite Sign Do Not Prove RH

The fixed-`omega` phase diagram has an exact cancellation witness:

```text
F(z)=(z^2+1)*(z^2+9),
F(z-i)/F(z+i)=(z-4i)/(z+4i).
```

The quotient is inner in the upper half-plane even though `F` has nonreal
zeros. A single good shift can therefore conceal off-axis zeros by exact
equal-height horizontal cancellation. Multiplicity is essential:

```text
F_def(z)=(z^2+1)*(z^2+9)^2,
F_def(z-i)/F_def(z+i)
 =(z-4i)^2*(z+2i)/((z-2i)*(z+4i)^2),
```

which leaves an uncancelled pole at `z=2i`. Pole-freeness and boundary
modulus are still insufficient without growth control: `exp(-i*z)` has
modulus one on the real line and no poles, but is unbounded in the upper
half-plane.

Suzuki's scalar eventual-sign target has a separate finite-prefix guard.
All 6000 sampled values of `sqrt(x)*h_omega^<1>(x)` are positive for five
shifts down to `omega=1/32` and `x<=5000`, but the primitive weight itself is
signed:

```text
g_(1/2)^<1>(exp(-4))=-39.803349588<0,
g_(1/2)^<1>(1/2)=1.174060048>0.
```

Hence coefficient positivity is not a termwise proof, and no finite grid
may be promoted to eventual sign or the required `L2` tail. The exact and
finite audits are:

```text
outputs/jensen_window_pf_suzuki_fixed_omega_phase_diagram.md
python work/rh_compute/scripts/check_jensen_window_pf_suzuki_fixed_omega_phase_diagram.py
outputs/jensen_window_pf_suzuki_jordan_totient_sign_scout.md
python work/rh_compute/scripts/check_jensen_window_pf_suzuki_jordan_totient_sign_scout.py
outputs/jensen_window_pf_suzuki_cofinal_monotonicity_hierarchy.md
python work/rh_compute/scripts/check_jensen_window_pf_suzuki_cofinal_monotonicity_hierarchy.py
outputs/jensen_window_pf_suzuki_jordan_error_kernel_reduction.md
python work/rh_compute/scripts/check_jensen_window_pf_suzuki_jordan_error_kernel_reduction.py
```

There is also an exact main-term guard. Abel summation shows that every
Suzuki smoothing kernel annihilates the positive residue main term of the
cumulative Jordan-totient mass. The remaining quantity is a signed
Mobius-error convolution. Elementary absolute estimates give only
`O(x^(1/2-omega)*(1+log x))`, with the same power at every finite smoothing
order. A proof may therefore use neither positive main-term domination nor
naive absolute values to establish eventual sign.

## Gate 11: Boundary-Unimodular Energy Is Not Half-Plane Analyticity

For `a>0`, consider the right-half-plane quotients

```text
Q_good(z)=(a-z)/(a+z),
Q_bad(z)=(a+z)/(a-z).
```

Both have modulus one on the imaginary axis. The first is inner; the second
has a pole at `z=a`. If `T_(k-1)` is the degree-`k-1` Taylor polynomial at
zero, exact division gives

```text
[Q_good(z)-T_(k-1)(z)]/z^k
 =2*(-1)^k/[a^(k-1)*(a+z)],

[Q_bad(z)-T_(k-1)(z)]/z^k
 =2/[a^(k-1)*(a-z)].
```

The first inverse Laplace transform decays like `exp(-a*t)` and is in
`L2(0,infinity)`. The second actual positive-time inverse grows like
`exp(a*t)` and is not. Yet the second quotient's regularized boundary trace
still has finite `L2` norm: its inverse Fourier transform is the
anti-causal function `2*a^(-(k-1))*exp(a*t)*1_(t<0)`. Thus the missing
Hardy condition is causality, not a finite boundary Plancherel norm. An
`L2` Suzuki residual detects right-half-plane pole-freeness; boundary
modulus one by itself does not provide it. A finite residual-energy cutoff
also cannot distinguish a long delayed tail from genuine integrability.

The checker verifies the Taylor-remainder identities exactly through six
smoothing orders, checks both boundary moduli, and integrates the decaying
model energy:

```text
outputs/jensen_window_pf_suzuki_cofinal_l2_hierarchy.md
work/rh_compute/results/jensen_window_pf_suzuki_cofinal_l2_hierarchy.json
python work/rh_compute/scripts/check_jensen_window_pf_suzuki_cofinal_l2_hierarchy.py
```

The model is a proof-safety guard, not a surrogate for xi. It rejects
trace-only and finite-energy promotions while leaving the cofinal
Jordan-error `L2` estimate open.

## Gate 12: Natural Mobius Form Does Not Supply Its Uniform Norm

The generalized Muntz identity writes the unsmoothed Jordan error as

```text
E_omega(x)
 =sum_(d>=1)mu(d)d^(-omega)R_omega(x/d).
```

This exact form has three distinct nonpromotion guards.

First, the scalar coefficients of the unitary dilation expansion are in
`l2`, but this alone gives no Bessel bound. The elementary Hilbert-space
model `v_d=v` with coefficients `1/d` has unit vectors and square-summable
coefficients while its partial sums diverge.

Second, the finite energy kernel

```text
K_X(u,v)=1/max(1,u,v)-1/X
```

is positive semidefinite and its ordered finite matrices are totally
nonnegative. That proves only `J_omega(X)>=0`, not a uniform upper bound.
The discrete, cross, and continuum pieces each have the same potentially
divergent scale, so their signed cancellation remains essential.

Third, the boundary transforms `1/(a+z)` and `1/(a-z)` have the same squared
boundary norm `1/(2a)`. The former is causal and decaying; the latter is
anti-causal, while its positive-time inverse grows. Consequently, even the
exact finite all-height Plancherel formula cannot be promoted to the Hardy
estimate without causality.

The exact Burnol-Hardy intertwiner adds a fourth guard. It gives

```text
(1-2omega)||t^(-omega)f_(2omega,N)||_2
 <=||E_(omega,N)||_H
 <=(1+2omega)||t^(-omega)f_(2omega,N)||_2,
```

so it transfers the norm problem into published Nyman coordinates but does
not bound the arithmetic fractional-part sum. Bounded invertibility of the
operator is not boundedness of its input family.

The reciprocal-cell reduction adds a fifth guard. Uniform boundedness is
equivalent to

```text
sup_N Q_(omega,N)<infinity,

Q_(omega,N)
 =sum_(k>=1)k^(2omega-2)
  |sum_(d<=N)mu(d)d^(-2omega){k/d}|^2.
```

The stable divisor prefix then forces

```text
|sum_(d<=N)mu(d)d^(-(1+2omega))-1/zeta(1+2omega)|
 =O(N^(-1/2-omega)).
```

Thus any proposed generic norm proof must also imply this explicit
RH-strength reciprocal-zeta tail rate and the remaining weighted
short-multiplicative-interval cancellation. An operator norm, positive
kernel, or coefficient `l2` estimate that cannot recover this rate has not
closed the arithmetic gate.

The finite-tail split adds a sixth guard. It gives the exact criterion

```text
Q_(omega,infinity)<infinity,
r_(omega,N)=O(N^(-1/2-omega)),
sup_N sum_(j>=0)(2^j*N)^(2omega-2)
 V_(omega,N)(2^j*N)<infinity.
```

The critical block estimate

```text
V_(omega,N)(K)=O(K^(2-2omega))
```

is not summable across dyadic `K=2^j*N` and therefore cannot be promoted
to the needed post-prefix energy. A proof needs a summable gain or a
different cancellation mechanism. Likewise, Báez-Duarte's unweighted
natural-approximation lower bounds and the unweighted fractional-part
autocorrelation kernel cannot be transferred by deleting the present
weights.

The stationary-Gram reduction adds a seventh guard. The weighted
autocorrelation kernel is positive definite and has nonnegative spectral
density, while

```text
Diag_(alpha,N)
 <=C_alpha(0)zeta(1+alpha)
```

is uniformly bounded. Neither fact bounds the signed Mobius off-diagonal.
Indeed, the pointwise positive-definite estimate
`|C_alpha(u)|<=C_alpha(0)` gives only

```text
Gamma_(alpha,N)=O(N^(1-alpha)),
```

which diverges. Kernel positivity, spectral positivity, and a bounded
diagonal therefore cannot be promoted to the uniform Burnol norm.

The Ornstein-Uhlenbeck reduction adds an eighth guard. The leading kernel
has the exact Gram identity

```text
G_(alpha,N)
 =sum_(k=1)^N [k^(1-alpha)-(k-1)^(1-alpha)]
  |sum_(d=k)^N mu(d)/d|^2.
```

Uniform boundedness of this energy for every member of one cofinal sequence
`alpha_j->0` is equivalent to RH, but this does not compare it with the
full Gram. The natural remainder density

```text
[|zeta((1-alpha)/2+it)|^2
 -(1-alpha)L_alpha]
/[((1-alpha)/2)^2+t^2]
```

is negative at `alpha=1/2`, `t=0`, so the remainder kernel is not positive
semidefinite. A positivity argument cannot discard it or order the full
Gram against the leading energy. Generic Montgomery-Vaughan spacing on
`log n` costs order `N^(1-alpha)`, and the multiplicative Hilbert matrix has
a product kernel rather than the required ratio kernel. A valid promotion
therefore needs a direct Mobius-specific remainder estimate.

The OU/Mertens reduction adds a ninth guard. If

```text
r_k=sum_(n>=k)mu(n)/n,
B_alpha=sum_(k>=1)k^(-alpha)|r_k|^2,
M_alpha=sum_(k>=1)M(k)^2/k^(2+alpha),
```

then explicit weighted Hardy and Copson bounds give

```text
B_alpha<infinity iff M_alpha<infinity.
```

The OU weights are comparable with `k^(-alpha)`, so this is not a new
generic Hilbert-space shortcut: it is the classical weighted Mertens
mean-square problem in exact coordinates. Equivalently,

```text
RH iff
sum_(K<=k<2K) M(k)^2=O_epsilon(K^(2+epsilon))
for every epsilon>0.
```

Moreover, with `C_h(Y)=sum_(m<=Y)mu(m)mu(m+h)`, exact pair reindexing
writes the off-diagonal part of `sum_(k<=X)M(k)^2` as

```text
2 sum_(h<X) sum_(Y<=X-h) C_h(Y).
```

Thus the remaining correlation target is signed, origin-anchored, and
cumulative in both the shift and terminal point. Generic Hardy applied
directly to the Mobius prefix pays
`sum mu(n)^2 n^(-alpha)`, which diverges. Averaged Chowla estimates average
all shift tuples, while almost-all short-interval estimates may discard the
exceptional interval anchored at the origin; neither theorem supplies this
anchored cumulative estimate. Even a terminal absolute estimate
`sum_h |C_h(Y)|=o(Y^2)` loses an extra factor when summed over `Y`. A valid
argument must preserve signed cancellation across the joint `(h,Y)` region
rather than take absolute values row by row.

The weighted-prefix/affine-defect reduction adds a tenth guard. For

```text
A_alpha(k)=sum_(n<=k)mu(n)n^(-alpha),
```

the continuous energy

```text
integral_1^infinity A_alpha(x)^2*x^(alpha-2)dx
```

is an exact Volterra-Hardy coordinate for the weighted Mertens energy and
has the same reciprocal-zeta boundary spectrum. Yet on the first
post-prefix block the q-discrepancy sees only

```text
A_alpha(k)-A_alpha(N)-k*r_N.
```

A generic coefficient sequence with nonzero constant prefix from `N`
onward has `r_N=0` and zero post-prefix discrepancy but positive absolute
weighted-prefix energy. Thus the dyadic discrepancy and reciprocal-tail
rate cannot be promoted at one cutoff to the Mertens energy.

The all-cutoff tail reduction adds the sharper logarithmic guard. Compatible
adjacent tails reconstruct the prefix exactly, and

```text
sum_N N^(alpha-2)|A_alpha(N)|^2<infinity
 iff
sum_N P_(alpha/2,N)/N<infinity.
```

The critical scalar sequence

```text
gamma=(1+alpha)/2,
r_N=(N+1)^(-gamma),
a_N=N[r_(N-1)-r_N]
```

satisfies `0<a_N<=gamma*N^(-alpha)` and has uniformly bounded stable-prefix
energy `P_N`, but `sum_N P_N/N` diverges. Equivalently, its weighted-prefix
energy has a positive harmonic asymptotic. Thus

```text
sup_N P_(omega,N)<infinity
```

cannot be promoted generically to the required all-cutoff anchor energy.
This model is not Mobius; it leaves a Mobius-specific logarithmic
mean-square gain open.

The fixed-shift spectral route has a separate promotion trap. The real even
quartic

```text
F(z)=((z-omega)^2+T^2)((z+omega)^2+T^2)
```

has a symmetric off-axis zero quartet, but

```text
Q(z)=F(z-omega)/F(z+omega)
```

cancels its boundary-zero factor and becomes rational inner in the right
half-plane. In fact `||Q-1||_(H2)^2=8omega`. Nevertheless,

```text
integral_R dt/[(1+t^2)|F(omega+it)|^2]=infinity
```

because the shifted line still passes through the simple zeros. Thus
finite limiting causal energy cannot be promoted generically to the
reciprocal boundary energy required by the logarithmic stable-prefix
criterion. The model blocks only this same-shift comparison; it does not
block the cofinal full-Burnol implication.

Vaughan decomposition has a second promotion trap. After collapsing the
anchored correlation, its test is

```text
f_X(n)=(X-n+1)M(n-1),
```

and generic Cauchy-Schwarz feeds `sum M(k)^2` back into the estimate.
Applying Vaughan before collapse avoids that circularity and gives the
exact bounded-test aggregate

```text
B_1(X)-TI_X+TII_X.
```

The new theorem target is a signed `O_epsilon(X^(2+epsilon))` estimate for
`-TI_X+TII_X`, with cancellation retained jointly over shift and terminal
length. Rowwise absolute values are still forbidden.

The executable audit is:

```text
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
outputs/jensen_window_pf_fixed_shift_inner_reciprocal_boundary_separation_gate.md
work/rh_compute/results/jensen_window_pf_fixed_shift_inner_reciprocal_boundary_separation_gate.json
python work/rh_compute/scripts/check_jensen_window_pf_fixed_shift_inner_reciprocal_boundary_separation_gate.py
```

It leaves open the equivalent uniform bounds
`sup_N ||E_(omega,N)||_H<infinity` and
`sup_N Q_(omega,N)<infinity` on a cofinal shift sequence.

## Current Cache-Based Gate Run

The executable gate reads:

```text
work/rh_compute/results/repro_hankel_15c_coefficients.jsonl
```

and applies the finite-prefix traps to the five lambda rows:

```text
0
1e-6
1e-4
1e-2
1e-1
```

With the current cache, the prefix is preserved through:

```text
K = 32
```

For each lambda, a positive one-term extension at `K+1` breaks:

```text
order-2 Toeplitz/PF
m = 1 signed-Hankel
degree-2 Jensen hyperbolicity
```

The same run also validates the exact rational moment-recurrence trap:

```text
recurrence order <= 12
moments 0..23 preserved
positive edited moment 24
next Hankel/moment gate breaks
```

It validates the exact finite shifted-principal signed-Hankel grid trap:

```text
base a_k = 1/k!
m <= 4, shifts <= 8 preserved
positive a_17 breaks the next shifted m = 1 signed-Hankel/Jensen gate
```

It validates the exact finite Schur-prefix trap:

```text
base h_k = 1/k!
h_0..h_6 preserved
2,940 finite Toeplitz/Schur tests preserved
positive h_7 breaks s_(6,6)
```

It validates the finite Jensen-window rectangle trap:

```text
current Jensen-window PF/Sturm coefficient inputs A_0..A_25 preserved
positive A_26 breaks degree-2, shift-24 Jensen hyperbolicity
```

It also validates the Stieltjes multiplier trap:

```text
positive measure = 10 delta_0 + delta_1 + delta_2
leading Hankel determinants = 12, 51, 40, 0
c_1^2 - c_0 c_2 = -1/4
```

Again, these are proof-safety models, not claims about the actual zeta coefficients. Their role is logical: they show that finite-prefix evidence, finite shifted-principal grids, finite recurrence evidence, and generic Stieltjes/Hankel positivity cannot be promoted into an all-order theorem by wording alone.

## Manuscript Rule

Additional signed-model gate:

```text
The raw ordinary Motzkin/J-fraction path model, diagonal sign conjugation,
global length-parity sign repairs, and absolute-value sign-state covers cannot
be used as manifest positivity proofs for E(t)=1/H(-t). A surviving signed
route must be a genuinely modified state-space, production/Riordan/network,
oscillatory resolvent, or Xi/Phi positive-cancellation construction, and it
must still prove all-order coefficientwise nonnegativity.
```

Before accepting any proposed bridge lemma, ask:

```text
Does the lemma fail on the local heat-birth model?
Does the lemma rely only on a finite coefficient prefix?
Does the lemma rely only on a finite Schur/Toeplitz shape prefix?
Does the lemma rely only on a finite shifted-principal signed-Hankel grid?
Does the lemma rely only on a finite moment or recurrence prefix?
Does the lemma prove only Stieltjes/Hankel positivity when Toeplitz PF is required?
Does the lemma silently assume PF-infinity, Jensen hyperbolicity, or Laguerre-Polya membership?
Does the lemma assume RH at lambda = 0?
```

If yes, it cannot be used to prove `Lambda <= 0`.

The only acceptable promotions are:

```text
finite certificate
conditional theorem with explicit hypotheses
all-order theorem with noncircular structural proof
```

Executable result-language scan:

```text
python work/rh_compute/scripts/check_output_reference_integrity.py
python work/rh_compute/scripts/check_output_status_manifest.py
python work/rh_compute/scripts/check_proof_claim_ledger.py
python work/rh_compute/scripts/check_signed_hankel_jensen_dependency_graph.py
python work/rh_compute/scripts/check_jensen_window_pf_bridge_obligations.py
python work/rh_compute/scripts/check_jensen_window_pf_theorem_machinery_fit_matrix.py
python work/rh_compute/scripts/check_jensen_window_pf_sign_regular_transfer_gap_matrix.py
python work/rh_compute/scripts/check_jensen_window_pf_factorial_multiplier_split_audit.py
python work/rh_compute/scripts/check_jensen_window_pf_structural_ansatz_matrix.py
python work/rh_compute/scripts/check_jensen_window_pf_schur_shape_contract.py
python work/rh_compute/scripts/check_jensen_window_pf_column_recurrence_contract.py
python work/rh_compute/scripts/check_jensen_window_pf_column_recurrence_finite_coverage.py
python work/rh_compute/scripts/check_arb_jensen_window_column_recurrence_stress.py
python work/rh_compute/scripts/check_jensen_window_pf_reciprocal_positivity_route_matrix.py
python work/rh_compute/scripts/check_jensen_window_pf_reciprocal_fraction_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_reciprocal_signed_j_fraction_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_signed_j_fraction_theorem_target.py
python work/rh_compute/scripts/check_jensen_window_pf_modified_signed_model_target.py
python work/rh_compute/scripts/check_jensen_window_pf_oscillatory_resolvent_fit_matrix.py
python work/rh_compute/scripts/check_jensen_window_pf_positive_readout_theorem_target.py
python work/rh_compute/scripts/check_jensen_window_pf_positive_spectral_moment_obstruction.py
python work/rh_compute/scripts/check_jensen_window_pf_nonordinary_positive_transform_ansatz_matrix.py
python work/rh_compute/scripts/check_jensen_window_pf_nonpower_functional_low_degree_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_nonpower_functional_cone_candidate_matrix.py
python work/rh_compute/scripts/check_jensen_window_pf_cauchy_binet_cone_frontier_matrix.py
python work/rh_compute/scripts/check_jensen_window_pf_monotone_contraction_frontier_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_monotone_contraction_column_extension_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_monotone_contraction_sparse_degree6_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_monotone_contraction_sparse_degree7_frontier_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_monotone_contraction_sparse_degree7_subdivision_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_monotone_contraction_all_m_counterexample.py
python work/rh_compute/scripts/check_jensen_window_pf_monotone_contraction_theorem_target.py
python work/rh_compute/scripts/check_jensen_window_pf_heat_flow_monotone_closure_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_heat_flow_boundary_threshold_lemma.py
python work/rh_compute/scripts/check_jensen_window_pf_kernel_mellin_upper_wall_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_log_concave_mellin_monotone_wall_countermodel.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_t1156_monotone_wall_counterexample_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_kernel_summand_shift_lemma.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_first_summand_dominance_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_m100_k320_collar_extension_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_first_summand_saddle_wall_target.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_first_summand_cumulant_bridge.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_heat_flow_ratio_cone_invariance_lemma.py
python work/rh_compute/scripts/check_jensen_window_pf_heat_flow_cone_entry_asymptotic_target.py
python work/rh_compute/scripts/check_jensen_window_pf_phi_taylor_cone_entry_sign_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_cone_entry_prefix_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_finite_collar_contract.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_tail_barrier_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_defect_recurrence_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_log_curvature_bridge.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_bounded_log_curvature_target.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_gaussian_curvature_matrix.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_signed_gaussian_perturbation_matrix.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_uniform_remainder_target.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_taylor_moment_budget.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_high_order_taylor_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_defect_tail_theorem_target.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_half_width_tail_target.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_adaptive_scaled_defect_target.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_adaptive_envelope_matrix.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_adaptive_envelope_obligations.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_raw_moment_bridge_matrix.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_raw_ratio_decrement_corridor_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_k300_precision_repair_audit.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_raw_log_decrement_bridge.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_coefficient_curvature_corridor_bridge.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_linear_curvature_barrier_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_scaled_curvature_monotonicity_target.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_scaled_curvature_log_ceiling_bridge.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_curvature_bridge.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_taylor_stencil_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_stencil_remainder_obligations.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_pointwise_tail_budget.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_next_increment_stencil_stress.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_degree16_stencil_continuation.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_degree16_collar_scan.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_degree16_real_t_collar_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_degree16_arb_real_t_collar_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_degree40_arb_collar_ladder_stress.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_degree40_residual_tail_budget.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_formal_tail_obstruction_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_asymptotic_remainder_target.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_actual_endpoint_remainder_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_bound_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_node_c0_range_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_grid_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_quadrature_ladder_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_quadrature_remainder_route_matrix.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_worst_row_far_tail_split_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_first_omitted_denominator_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_coefficient_core_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_worst_row_laguerre_root_bracket_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_worst_row_christoffel_weight_midpoint_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_worst_row_christoffel_weight_interval_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_part_weighted_sum_interval_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_plus_tail_budget_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_raw_moment_obstruction_matrix.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target.py
python work/rh_compute/scripts/check_jensen_window_pf_monotone_contraction_stress.py
python work/rh_compute/scripts/check_jensen_window_pf_state_space_sign_lift_obstruction_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_cauchy_binet_low_degree_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_log_concavity_frontier_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_ratio_condition_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_contraction_log_concavity_scout.py
python work/rh_compute/scripts/check_result_language_boundaries.py
```

Current result:

```text
validated output references: scanned 425 markdown files, 4383 path references, 0 missing required paths, 3 planned missing deliverables
validated output artifact statuses: scanned 425 markdown files, 0 status issues
validated proof-claim ledger: 359 claims, 0 issues, 9 open theorem targets
validated signed-Hankel/Jensen dependency graph with 0 issues
validated Jensen-window PF coefficient-PF equivalence gate: 13 rows, 0 issues, 3 exact coefficient identities, 4 classical/closure steps, 1 seven-way equivalence, 3 guards, 1 open structural handoff
validated Jensen-window PF Edrei-Stieltjes equivalence gate: 14 rows, 0 issues, 12 exact indexing checks, 7 exact Hankel checks, 1 unified endpoint, 3 finite/nonpromotion guards, 1 open Xi/Phi handoff
validated Jensen-window PF Edrei heat-flow boundary gate: 11 rows, 0 issues, 4 exact flow identities, 9 rank-one orientation checks, 2 exact heat witnesses, 1 rejected generic backward-invariance shortcut, 1 open Xi/Phi handoff
validated Jensen-window PF Phi Pick-kernel target: 12 rows, 0 issues, 7 exact identities, 3 independent polarization checks, 1 exact mixture guard, 2 open structural routes
validated Jensen-window PF Xi Pick/Suzuki Hankel bridge: 16 rows, 0 issues, 7 exact coordinate/guard identities, 4 published operator steps, 3 exact countermodels, 1 open global gate
validated Jensen-window PF Suzuki spectral frontier: 16 rows, 0 issues, 10 exact path/reduction identities, 5 route guards, 1 open global obligation
validated Suzuki determinant-only reduction: 16 rows, 0 issues, 12 exact bridge steps, 2 theorem/corollary candidates, 1 route guard, 1 open arithmetic gate
validated Suzuki fixed-omega phase diagram: 21 rows, 0 issues, 3 exact countermodel guards, 2 published scalar targets, 1 open arithmetic gate
validated Suzuki Jordan-totient sign scout: 6000 samples, 5 omega values, 0 negative rows, 0 issues
validated Suzuki cofinal monotonicity hierarchy: 21 rows, 0 issues, 3 signed/nonpromotion guards, 1 cofinal equivalence candidate, 1 open arithmetic gate
validated Suzuki Jordan-error kernel reduction: 16 rows, 0 issues, 8 exact identities/bounds, 2 route guards, 1 open cancellation gate
validated Suzuki cofinal L2 hierarchy: 22 rows, 0 issues, 3 countermodel guards, 1 literature-fit guard, 1 cofinal equivalence candidate, 1 open arithmetic energy gate
validated Jordan-Muntz causal-energy bridge: 32 rows, 0 issues, 5 theorem candidates, 5 proof guards, 1 open all-height mollifier gate
validated Jordan-Muntz/Burnol Hardy intertwiner: 18 rows, 0 issues, 11 exact identities, 1 source-backed cofinal theorem candidate, 2 nonpromotion guards, 1 open natural-mollifier gate
validated Burnol cell-energy/tail obstruction: 18 rows, 0 issues, 10 exact identities, 4 norm-reduction steps, 1 reciprocal-zeta tail obstruction, 1 open short-interval gate
validated Burnol tail-discrepancy/dyadic reduction: 20 rows, 0 issues, 12 exact identities, 4 equivalence steps, 2 literature guards, 1 open dyadic square-function gate
validated weighted fractional autocorrelation Gram bridge: 20 rows, 0 issues, 12 exact identities, 3 proof guards, 1 open signed off-diagonal gate
validated weighted-autocorrelation OU tail-energy reduction: 22 rows, 0 issues, 13 exact identities, 4 proof guards, 1 open comparison gate
validated OU/Mertens mean-square reduction: 25 rows, 0 issues, 16 exact reductions, 4 proof guards, 1 open anchored mean-square gate
validated Mertens weighted-prefix/affine-defect reduction: 34 rows, 0 issues, 24 exact reductions, 4 proof guards, 2 open handoff gates
validated Jensen-window PF bridge obligations: 16 obligations, 0 issues, 3 open obligations
validated Jensen-window PF theorem machinery fit matrix: 11 rows, 0 issues, 0 ready-to-apply rows
validated Jensen-window PF sign-regular transfer gap matrix: 9 transfer rows, 2 countermodel gates, 3 open requirements, 3 rejected shortcuts, 0 ready-to-apply rows, 0 issues
validated Jensen-window PF factorial multiplier split audit: 5 exact rows, 315 raw degree-2 anti-hyperbolic rows, 315 normalized degree-2 positive rows, 0 ready-to-apply rows, 0 issues
validated Jensen-window PF structural ansatz matrix: 6 ansatz rows, 0 issues, 0 ready-to-apply rows
validated Jensen-window PF Schur shape contract: 4 grid rows, 0 issues, 2 frontier rows
validated Jensen-window PF column recurrence contract: 4 degree rows, 0 issues, 2 hard frontier rows
validated Jensen-window PF column recurrence finite coverage: 1470 direct positive rows, 210 hard recurrence rows, 315 Sturm/PF windows, 0 issues
validated 12600 Arb Jensen-window column recurrence stress rows with 0 issues
validated Jensen-window PF reciprocal positivity route matrix: 9 rows, 0 issues, 0 ready-to-apply rows
validated Jensen-window PF reciprocal fraction scout: 3 symbolic rows, 735 finite rows, 0 issues
validated Jensen-window PF reciprocal signed J-fraction scout: 2 symbolic rows, 3675 signed Hankel rows, 2940 signed-lambda rows, 0 issues
validated Jensen-window PF signed J-fraction theorem target: 7 fit rows, 0 issues, 0 ready-to-apply rows
validated Jensen-window PF modified signed-model target: 9 model rows, 0 issues, 0 ready-to-apply rows, 4 live modified candidates
validated Jensen-window PF oscillatory resolvent fit matrix: 8 fit rows, 0 issues, 0 ready-to-apply rows
validated Jensen-window PF positive readout theorem target: 8 candidate rows, 0 issues, 0 ready-to-apply rows, 2 live foundational routes
validated Jensen-window PF positive spectral moment obstruction: 3 symbolic rows, 735 finite Delta2 obstruction rows, 0 issues
validated Jensen-window PF nonordinary positive transform ansatz matrix: 8 ansatz rows, 0 issues, 0 ready-to-apply rows, 3 live ansatz rows
validated Jensen-window PF nonpower functional low-degree scout: 7 scout rows, 0 issues, 0 ready-to-apply rows, 1 live contract rows
validated Jensen-window PF nonpower functional cone candidate matrix: 8 cone rows, 0 issues, 0 ready-to-apply rows, 2 live cone rows
validated Jensen-window PF Cauchy-Binet cone frontier matrix: 8 frontier rows, 0 issues, 0 ready-to-apply rows, 2 live frontier rows
validated Jensen-window PF monotone contraction frontier scout: 2 exact rows, 88 Bernstein coefficients, 210 finite zeta rows, 0 issues
validated Jensen-window PF monotone-contraction column extension scout: 25 column rows, 3329 Bernstein coefficients, 3 beyond-frontier rows, 0 negative Bernstein rows, 0 issues
validated Jensen-window PF monotone-contraction sparse degree-6 scout: 10 degree-6 rows, 63347 Bernstein coefficients, m<=10, 0 negative Bernstein rows, 0 zero Bernstein rows, 0 issues
validated Jensen-window PF monotone-contraction sparse degree-7 frontier scout: 9 positive rows, 1 certificate-obstruction row, 932691 Bernstein coefficients, first obstruction m=10, 126 negative Bernstein coefficients, 0 zero Bernstein coefficients, 0 issues
validated Jensen-window PF monotone-contraction sparse degree-7 subdivision scout: 3 dyadic slabs, 785400 slab Bernstein coefficients, 0 negative slab coefficients, 0 zero slab coefficients, repaired m=10 obstruction, 0 issues
validated Jensen-window PF monotone-contraction all-m counterexample: degree 7, m=11, exact monotone witness, negative normalized value, 0 issues
validated Jensen-window PF monotone contraction theorem target: 9 rows, 0 issues, 0 ready-to-apply rows, 2 live routes
validated Jensen-window PF heat-flow monotone closure scout: 4 exact rows, 315 threshold rows, 305 flow-bracket rows, 0 issues
validated Jensen-window PF heat-flow boundary threshold lemma: 5 exact rows, 315 strong-threshold rows, 315 heat-threshold rows, 0 issues
validated Jensen-window PF kernel Mellin upper-wall certificate: 8 rows, 0 issues, 200 positive compact intervals, 1 positive analytic ray, 1 remaining open cone clause, 0 ready-to-apply rows
validated Jensen-window PF log-concave Mellin monotone-wall countermodel: 6 rows, 0 issues, 2 upper-wall contractions, 1 monotone-wall violation
validated Jensen-window PF T=1156 monotone-wall counterexample certificate: 7 rows, 0 issues, 4 coefficient enclosures, 1 zeta monotone-wall violation
validated Jensen-window PF negative-lambda kernel summand-shift lemma: 8 rows, 0 issues, 6 exact rows, 1 compact interval row, 1 open far-tail row, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda first-summand dominance certificate: 10 rows, 0 issues, 4 exact rows, 5 interval rows, 15 positive analytic gates, 1 open dominant-wall row, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda -100 k320 collar extension certificate: 6 rows, 0 issues, 76 positive coefficients, 74 cone rows, 73 adjacent-wall rows, 19 new extension rows, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda first-summand saddle-wall target: 9 rows, 0 issues, 9 positive samples, 9 quarter-k2 samples, 9 bracketed saddles, 1 open requirement, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda first-summand cumulant bridge: 8 rows, 0 issues, 4 exact identities, 1 conditional bridge, 9 positive samples, 1 open requirement, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda first-summand leading-saddle certificate: 8 rows, 0 issues, 40740 positive leading intervals, 40740 positive cubic-correction intervals, 40740 positive fifth-correction intervals, 3 positive analytic ray gates, 9 positive seventh-remainder samples, 1 open remainder, 0 ready-to-apply rows
validated Jensen-window PF heat-flow ratio cone invariance lemma: 6 exact rows, 315 lower rows, 315 upper rows, 310 monotone rows, 0 issues
validated Jensen-window PF heat-flow cone-entry asymptotic target: 8 rows, 0 issues, 0 ready-to-apply rows, 1 live routes
validated Jensen-window PF Phi Taylor cone-entry sign scout: 4 coefficient balls, 2 certified signs, 0 ready-to-apply rows, 0 issues
validated Jensen-window PF negative-lambda cone-entry prefix scout: 69 coefficient rows, 63 lower-wall rows, 63 upper-wall rows, 60 monotone-gap rows, 0 issues
validated Jensen-window PF negative-lambda cone-entry prefix scout: 183 coefficient rows, 177 lower-wall rows, 177 upper-wall rows, 174 monotone-gap rows, 0 issues
validated Jensen-window PF negative-lambda cone-entry prefix scout: 243 coefficient rows, 237 lower-wall rows, 237 upper-wall rows, 234 monotone-gap rows, 0 issues
validated Jensen-window PF negative-lambda cone-entry prefix scout: 303 coefficient rows, 297 lower-wall rows, 297 upper-wall rows, 294 monotone-gap rows, 0 issues
validated Jensen-window PF negative-lambda cone-entry prefix scout: 453 coefficient rows, 447 lower-wall rows, 447 upper-wall rows, 444 monotone-gap rows, 0 issues
validated Jensen-window PF negative-lambda cone-entry prefix scout: 603 coefficient rows, 597 lower-wall rows, 597 upper-wall rows, 594 monotone-gap rows, 0 issues
validated Jensen-window PF negative-lambda finite-collar contract: active depth K=19, 57 active lower rows, 57 active upper rows, 57 active monotone rows, 6 collar lower/upper rows, 3 collar monotone rows, 0 issues
validated Jensen-window PF negative-lambda finite-collar contract: active depth K=57, 171 active lower rows, 171 active upper rows, 171 active monotone rows, 6 collar lower/upper rows, 3 collar monotone rows, 0 issues
validated Jensen-window PF negative-lambda finite-collar contract: active depth K=77, 231 active lower rows, 231 active upper rows, 231 active monotone rows, 6 collar lower/upper rows, 3 collar monotone rows, 0 issues
validated Jensen-window PF negative-lambda finite-collar contract: active depth K=97, 291 active lower rows, 291 active upper rows, 291 active monotone rows, 6 collar lower/upper rows, 3 collar monotone rows, 0 issues
validated Jensen-window PF negative-lambda finite-collar contract: active depth K=147, 441 active lower rows, 441 active upper rows, 441 active monotone rows, 6 collar lower/upper rows, 3 collar monotone rows, 0 issues
validated Jensen-window PF negative-lambda finite-collar contract: active depth K=197, 591 active lower rows, 591 active upper rows, 591 active monotone rows, 6 collar lower/upper rows, 3 collar monotone rows, 0 issues
validated Jensen-window PF negative-lambda tail-barrier scout: 63 cone-buffer rows, 60 defect-monotone rows, 63 one-third-buffer rows, 60 scaled-defect increase rows, 1 rejected candidate, 0 issues
validated Jensen-window PF negative-lambda tail-barrier scout: 177 cone-buffer rows, 174 defect-monotone rows, 139 one-third-buffer rows, 174 scaled-defect increase rows, 1 rejected candidate, 0 issues
validated Jensen-window PF negative-lambda scaled-defect frontier scout: 177 scaled rows, 177 cone rows, 177 half-width rows, 139 one-third rows, 38 one-third failures, 174 scaled-increase rows, 0 issues
validated Jensen-window PF negative-lambda tail-barrier scout: 237 cone-buffer rows, 234 defect-monotone rows, 159 one-third-buffer rows, 234 scaled-defect increase rows, 1 rejected candidate, 0 issues
validated Jensen-window PF negative-lambda scaled-defect frontier scout: 237 scaled rows, 237 cone rows, 237 half-width rows, 159 one-third rows, 78 one-third failures, 234 scaled-increase rows, 0 issues
validated Jensen-window PF negative-lambda tail-barrier scout: 297 cone-buffer rows, 294 defect-monotone rows, 179 one-third-buffer rows, 294 scaled-defect increase rows, 1 rejected candidate, 0 issues
validated Jensen-window PF negative-lambda scaled-defect frontier scout: 297 scaled rows, 297 cone rows, 297 half-width rows, 179 one-third rows, 118 one-third failures, 294 scaled-increase rows, 0 issues
validated Jensen-window PF negative-lambda tail-barrier scout: 447 cone-buffer rows, 444 defect-monotone rows, 179 one-third-buffer rows, 444 scaled-defect increase rows, 1 rejected candidate, 0 issues
validated Jensen-window PF negative-lambda scaled-defect frontier scout: 447 scaled rows, 447 cone rows, 430 half-width rows, 179 one-third rows, 268 one-third failures, 444 scaled-increase rows, 0 issues
validated Jensen-window PF negative-lambda tail-barrier scout: 597 cone-buffer rows, 594 defect-monotone rows, 179 one-third-buffer rows, 594 scaled-defect increase rows, 1 rejected candidate, 0 issues
validated Jensen-window PF negative-lambda scaled-defect frontier scout: 597 scaled rows, 597 cone rows, 521 half-width rows, 179 one-third rows, 418 one-third failures, 594 scaled-increase rows, 0 issues
validated Jensen-window PF negative-lambda defect-recurrence scout: 63 buffered rows, 60 defect-monotone rows, 60 width-recurrence rejections, 1 live sufficient routes, 0 issues
validated Jensen-window PF negative-lambda log-curvature bridge: 63 simple log-buffer rows, 63 exact defect-buffer rows, 60 curvature-monotone rows, 5 bridge rows, 0 issues
validated Jensen-window PF negative-lambda bounded log-curvature target: 8 rows, 0 issues, 0 ready-to-apply rows, 2 live routes, 63 raw-threshold rows
validated Jensen-window PF negative-lambda bounded log-curvature k300 obstruction: 7 rows, 0 issues, 718 two-thirds failures, 894 scaled-curvature increase rows, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda Gaussian curvature matrix: 7 matrix rows, 63 positive-deficit rows, 63 bounded-deficit rows, 63 raw-threshold rows, 0 issues
validated Jensen-window PF negative-lambda signed Gaussian perturbation matrix: 8 matrix rows, 2 certified Taylor signs, 1 fixed-k activation estimates, 0 ready-to-apply rows, 0 issues
validated Jensen-window PF negative-lambda uniform remainder target: 8 rows, 0 issues, 0 ready-to-apply rows, 2 open requirements, 3 leading-scale rows
validated Jensen-window PF negative-lambda Taylor moment budget: 9 budget rows, 7 tail-start samples, 4 invalid truncation rows, 2 bounded truncation rows, 0 ready-to-apply rows, 0 issues
validated Jensen-window PF negative-lambda high-order Taylor scout: 8 coefficient rows, 35 truncation rows, 9 invalid normalizers, 2 upper-wall violations, 3 overbound rows, 0 ready-to-apply rows, 0 issues
validated Jensen-window PF negative-lambda defect-tail theorem target: 8 rows, 0 issues, 0 ready-to-apply rows, 2 live routes
validated Jensen-window PF negative-lambda half-width tail target: 9 rows, 0 issues, 0 ready-to-apply rows, 0 live routes, 430 half-width rows, 17 half-width failures
validated Jensen-window PF negative-lambda adaptive scaled-defect target: 8 rows, 0 issues, 2 live routes, 597 exact-cone rows, 76 half-width failures
validated Jensen-window PF negative-lambda adaptive envelope matrix: 7 matrix rows, 0 issues, 594 k-increase rows, 398 lambda-order rows, 76 half-width failures
validated Jensen-window PF negative-lambda adaptive envelope obligations: 9 obligation rows, 0 issues, 3 exact rows, 3 open requirements, 1 rejected routes
validated Jensen-window PF negative-lambda raw-moment bridge matrix: 8 matrix rows, 0 issues, 597 raw-cone rows, 594 corridor rows, 76 half-width failures
validated Jensen-window PF negative-lambda raw-ratio decrement-corridor scout: 9 rows, 0 issues, 594 decrement-corridor rows, 591 theta-k-monotone rows, 2 exact counterexamples, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda k300 precision-repair audit: 7 rows, 0 issues, 894 repaired decrement-corridor rows, 891 repaired theta-k-monotone rows, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda raw-log decrement bridge: 8 rows, 0 issues, 894 log-corridor rows, 894 log-decrease rows, 2 exact counterexamples, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda coefficient-curvature corridor bridge: 9 rows, 0 issues, 894 curvature-corridor rows, 894 monotone-curvature rows, 2 exact counterexamples, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda linear curvature-barrier scout: 8 rows, 0 issues, 894 linear-barrier rows, 894 monotone-curvature rows, 2 exact counterexamples, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda scaled-curvature monotonicity target: 10 rows, 0 issues, 2 live routes, 894 scaled-curvature increase rows, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda scaled-curvature log-ceiling bridge: 8 rows, 0 issues, 894 scaled-ceiling rows, 894 scaled-log-corridor rows, 894 ceiling-dominance rows, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda relative-Gaussian curvature bridge: 8 rows, 0 issues, 897 B-positive rows, 894 B-decrease rows, 894 C-increase rows, 598 C-lambda-order rows, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda relative-Gaussian Taylor stencil scout: 8 rows, 0 issues, 3 certified leading-sign rows, 35 truncation rows, 4 all-positive stencil rows, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda relative-Gaussian stencil remainder obligations: 9 rows, 0 issues, 4 positive baseline rows, 31 blocked baseline rows, 4 exact stencil rows, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda relative-Gaussian pointwise tail budget: 9 rows, 0 issues, 4 positive baseline rows, 4 budget rows, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda relative-Gaussian next-increment stencil stress: 8 rows, 0 issues, 2 tested next-increment rows, 2 pointwise budget failures, 2 stencil-sign-preserving rows, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda relative-Gaussian degree-16 stencil continuation: 7 rows, 0 issues, 4 tested continuation rows, 3 stencil-sign-preserving rows, 1 stencil-sign-failure rows, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda relative-Gaussian degree-16 collar scan: 7 rows, 0 issues, 1301 scan rows, 1045 continuation-positive rows, 718 half-safety rows, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda relative-Gaussian degree-16 real-T collar scout: 8 rows, 0 issues, 4 positive normalizer rows, 3 certified surrogate stencil rows, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda relative-Gaussian degree-16 Arb real-T collar certificate: 8 rows, 0 issues, 4 positive normalizer rows, 3 certified stencil rows, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda relative-Gaussian degree-40 Arb collar ladder stress: 8 rows, 0 issues, 13 degree levels, max degree 40, 0 failed Bernstein rows, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda relative-Gaussian degree-40 residual tail budget: 8 rows, 0 issues, 5 budget inequalities, 4 finite tail profile rows, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda relative-Gaussian formal-tail obstruction scout: 8 rows, 0 issues, 4 profile rows, 4 formal-tail turnaround rows, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda relative-Gaussian asymptotic remainder target: 6 rows, 0 issues, 4 first-omitted rows, 4 optimized-window rows, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda relative-Gaussian actual endpoint remainder scout: 6 rows, 0 issues, 4 endpoint rows, 5 quadrature orders, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda relative-Gaussian cancellation-reduced remainder grid scout: 6 rows, 0 issues, 20 grid rows, 5 T values, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda relative-Gaussian intervalization target: 6 rows, 0 issues, 8 obligations, 5 open requirements, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda relative-Gaussian Phi tail bound scout: 6 rows, 0 issues, 3 tail bounds below 1e-1000, 2 conditional requirements, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda relative-Gaussian node-c0 range certificate: 5 rows, 0 issues, 16 Laguerre bound rows, 2 certified side conditions, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda relative-Gaussian Phi-tail grid certificate: 6 rows, 0 issues, 3 certified tail sources, 2 certified side conditions, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda relative-Gaussian quadrature ladder scout: 5 rows, 0 issues, 7 ladder rows, 320 reference order, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda relative-Gaussian quadrature-remainder route matrix: 7 rows, 0 issues, derivative order 640, 2 derivative-sup caps, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda relative-Gaussian worst-row far-tail split certificate: 7 rows, 0 issues, split y=200, 2 tail ratios below quadrature cap, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda relative-Gaussian worst-row compact-interval integration scout: 7 rows, 0 issues, 6 panels, plain interval Riemann rejected, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda relative-Gaussian worst-row Chebyshev panel-moment scout: 7 rows, 0 issues, 5 degrees, 4 Cauchy pairs, 3 cap-safe pairs, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda relative-Gaussian worst-row Arb Chebyshev interpolant-moment scout: 7 rows, 0 issues, 4 degrees, 3 Cauchy pairs, 3 cap-safe pairs, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda relative-Gaussian first-omitted denominator certificate: 6 rows, 0 issues, 20 denominator rows, 2 ratio-cap rows, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda relative-Gaussian coefficient-core propagation certificate: 6 rows, 0 issues, 22 coefficient rows, 20 propagation rows, 2 intervalization rows, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda relative-Gaussian worst-row Laguerre root-bracket certificate: 6 rows, 0 issues, 320 root brackets, 30 zero floating weights, 2 intervalization rows, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda relative-Gaussian worst-row Christoffel-weight midpoint scout: 6 rows, 0 issues, 320 midpoint weights, 30 repaired floating underflows, 320 direct interval obstructions, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda relative-Gaussian worst-row Christoffel-weight interval certificate: 6 rows, 0 issues, 320 interval weights, 0 Taylor denominator obstructions, 30 repaired floating underflows, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda relative-Gaussian worst-row finite-part weighted-sum interval certificate: 6 rows, 0 issues, 320 refined nodes, 320 interval weights, 2 below-one ratios, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda relative-Gaussian worst-row finite-plus-tail budget certificate: 6 rows, 0 issues, 2 composed ratios, 3 tail sources, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda raw-moment obstruction matrix: 7 matrix rows, 0 issues, 3 exact counterexamples, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda zeta-specific raw-corridor target: 9 rows, 0 issues, 2 live routes, 2 rejected shortcuts, 0 ready-to-apply rows
validated Jensen-window PF monotone contraction stress: 2875 rows, 2875 positive rows, 0 issues
validated Jensen-window PF state-space sign-lift obstruction scout: 3 symbolic rows, 735 mu2 sign-lift obstruction rows, 0 issues
validated Jensen-window PF Cauchy-Binet low-degree scout: 15 formula rows, 0 issues, 0 kernel identities found
validated Jensen-window PF log-concavity frontier scout: 14 contiguous rows, 0 issues
validated Jensen-window PF ratio-condition scout: 7 candidate rows, 0 issues, 4 rejected by countermodel, 1 rejected by construction
validated Jensen-window PF contraction-log-concavity scout: 1 rejected by construction, 0 issues, 2 negative frontier rows
validated result-language boundaries: scanned 423 markdown files, 0 overclaims
```
