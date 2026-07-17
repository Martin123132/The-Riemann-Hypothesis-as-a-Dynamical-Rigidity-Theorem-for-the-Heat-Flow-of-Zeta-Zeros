# Formal Core

Date: 2026-07-13

Status: exact-lemma ledger. This is not a proof of RH or `Lambda <= 0`; it records noncircular identities and separates them from conjectural bridges.

Project:

```text
RH Dynamical Rigidity / de Bruijn-Newman Zero-Flow Programme
```

Purpose:

Collect the exact, non-circular mathematical statements that can safely be used in later proof development. Numerical observations, conjectural bridges, and RH-equivalent assumptions are deliberately separated from exact lemmas.

## 0. Conventions

We use the standard de Bruijn-Newman family

```text
H_lambda(z) = integral_0^infty exp(lambda u^2) Phi(u) cos(z u) du.
```

It satisfies the heat-type equation

```text
partial_lambda H_lambda(z) = - partial_zz H_lambda(z).
```

The Newman parameter is `lambda`.

The de Bruijn-Newman constant `Lambda` is defined so that `H_lambda` has only real zeros in the real-zero regime. With the usual convention:

```text
RH is equivalent to Lambda <= 0.
```

Rodgers-Tao proves:

```text
Lambda >= 0.
```

Therefore this programme must prove:

```text
Lambda <= 0
```

without assuming RH at `lambda = 0`.

## 1. Zero Tracking

### Lemma 1.1: Simple Zero Tracking

Let `F(lambda,z)` be analytic in both variables and real-valued on the real axis for real `lambda`. Suppose

```text
F(lambda_0, x_0) = 0,
partial_z F(lambda_0, x_0) != 0.
```

Then there is a unique differentiable real branch `x(lambda)` near `lambda_0` such that

```text
x(lambda_0) = x_0,
F(lambda, x(lambda)) = 0.
```

The branch satisfies

```text
x'(lambda) = - F_lambda(lambda,x(lambda)) / F_z(lambda,x(lambda)).
```

Proof:

Immediate from the implicit function theorem and differentiation of

```text
F(lambda, x(lambda)) = 0.
```

### Corollary 1.2: Newman Zero Velocity

For `F = H_lambda`, since

```text
F_lambda = -F_zz,
```

any simple real zero branch satisfies

```text
x_i'(lambda) = F_zz(lambda,x_i) / F_z(lambda,x_i).
```

This is exact wherever the zero remains simple.

## 2. Local Product Expansion And Singular Velocity

### Lemma 2.1: Local Hadamard/Weierstrass Velocity Form

Suppose that near a finite cluster of simple real zeros

```text
x_1(lambda), ..., x_r(lambda)
```

the function can be written as

```text
F(lambda,z) = U(lambda,z) product_{k=1}^r (z - x_k(lambda)),
```

where `U` is analytic and nonzero in the neighbourhood.

Then at a zero `x_i`,

```text
F_zz/F_z (x_i)
= 2 sum_{1 <= j <= r, j != i} 1/(x_i - x_j)
  + 2 partial_z log U(x_i).
```

For the Newman flow,

```text
x_i'
= 2 sum_{j != i} 1/(x_i - x_j)
  + E_i,
```

where

```text
E_i = 2 partial_z log U(x_i)
```

is the local external field contributed by all structure outside the finite cluster.

Proof sketch:

Write `P(z) = product_k (z - x_k)`. At a simple zero `x_i`,

```text
P_zz/P_z (x_i) = 2 sum_{j != i} 1/(x_i - x_j).
```

The logarithmic derivative of the nonvanishing analytic factor contributes `2 U_z/U`.

## 3. Adjacent Gap Equation In The Finite Pairwise Model

Consider the reduced finite model

```text
x_i' = 2 sum_{j != i} 1/(x_i - x_j),
```

with ordered real points

```text
x_1 < x_2 < ... < x_n.
```

Let

```text
g_i = x_{i+1} - x_i.
```

Then

```text
g_i'
= 4/g_i
  - 2 g_i sum_{j != i,i+1}
      1 / ((x_{i+1}-x_j)(x_i-x_j)).
```

For every outside point `x_j`, the product

```text
(x_{i+1}-x_j)(x_i-x_j)
```

is positive. Therefore the pairwise tail in this reduced model is compressive:

```text
tail <= 0.
```

Consequence:

A bridge lemma asserting nonnegative pairwise tail is false for the reduced finite model. Any useful gap comparison theorem must use additional Xi-specific structure, a controlled-negativity estimate, or a different global principle.

## 4. Finite Cluster Non-Collapse In The Repulsive Direction

### Lemma 4.1: Internal Variance Growth

For the internal reduced flow on `r` points

```text
x_i' = 2 sum_{j != i} 1/(x_i - x_j),
```

define the cluster mean

```text
m = (1/r) sum_i x_i
```

and internal variance

```text
V = sum_i (x_i - m)^2.
```

Then the internal pairwise contribution satisfies

```text
V'_internal = 2r(r-1).
```

Proof:

Since `sum_i (x_i-m)=0`,

```text
V' = 2 sum_i (x_i-m) x_i'.
```

Pairing terms `i,j` gives one unit contribution per ordered interaction pair, hence

```text
V'_internal
= 4 sum_{i<j} 1
= 2r(r-1).
```

### Lemma 4.2: Non-Collapse With Locally Lipschitz External Field

Suppose a finite cluster evolves by

```text
x_i'
= 2 sum_{j != i} 1/(x_i - x_j)
  + E(x_i,lambda),
```

and suppose `E` is locally Lipschitz in `x` across the cluster, with Lipschitz constant `L`.

Then

```text
V' >= 2r(r-1) - 2LV.
```

Hence a positive-variance cluster cannot collapse to `V=0` in finite forward `lambda` time in the repulsive direction.

Interpretation:

This is a genuine local theorem about already-real finite clusters. It does not exclude birth of a real pair from a complex conjugate pair at a positive Newman boundary.

## 5. Square-Root Boundary Normal Form

### Lemma 5.1: Local Double-Zero Normal Form

Let `F(lambda,z)` be analytic and real on the real axis. Suppose

```text
F(lambda_*, t_*) = 0,
F_z(lambda_*, t_*) = 0,
F_zz(lambda_*, t_*) != 0,
```

and the parameter is transverse:

```text
F_lambda(lambda_*, t_*) != 0.
```

Then, locally,

```text
F(lambda,z)
= U(lambda,z) [ (z - t_* - r(lambda))^2 - s(lambda) ],
```

where `U` is analytic and nonzero, and `s(lambda_*)=0` with `s'(lambda_*) != 0`.

The local zeros are

```text
z_+(lambda), z_-(lambda)
= t_* + r(lambda) +/- sqrt(s(lambda)).
```

Thus the transition between a conjugate pair and two real zeros is square-root type.

## 6. Minimal Obstruction To The Local-Repulsion Proof Strategy

Consider

```text
F(lambda,t) = (t-a)^2 - 2(lambda-lambda_*).
```

Then

```text
F_lambda = -F_tt.
```

The zeros are:

```text
lambda < lambda_*:  t = a +/- i sqrt(2(lambda_*-lambda))
lambda = lambda_*:  t = a, double real zero
lambda > lambda_*:  t = a +/- sqrt(2(lambda-lambda_*))
```

On the real-zero side,

```text
g(lambda) = 2 sqrt(2(lambda-lambda_*))
```

and

```text
g' = 4/g.
```

Therefore the local repulsive law is compatible with positive square-root birth. It cannot, by itself, prove `Lambda <= 0`.

## 7. Exact Implication Needed For RH

A non-circular proof through the Newman route must establish:

```text
There is no lambda_* > 0 and real t_* such that
H_{lambda_*}(t_*) = 0
partial_t H_{lambda_*}(t_*) = 0.
```

Equivalently, it must rule out a positive Newman boundary without assuming all zeros are real at `lambda = 0`.

## 8. Conditional Stability Statement

The local cluster and gap machinery may support conditional statements of the form:

```text
If all zeros are already real at some lambda_0,
and if the relevant external-field/comparison hypotheses hold,
then real-rootedness persists for larger lambda.
```

This is useful, but it is not a proof of RH if `lambda_0 = 0`, because that assumption is RH itself.

## 9. Numerical Evidence Status

The existing numerical material supports:

- reduced finite-window rank-2 structure;
- Hadamard alignment of the leading mode;
- density-balancing secondary mode;
- signed Hankel/Jensen patterns in tested finite ranges.

These are evidence and guideposts. They are not currently exact theorems about the full de Bruijn-Newman family.

## 10. Next Formal Target

The next non-circular theorem search should focus on one of:

```text
A. no positive Newman boundary;
B. Laguerre-Polya membership for H_lambda in the required lambda range;
C. signed total positivity / sign-regularity sufficient for real-rootedness;
D. Xi-specific comparison principle controlling the global tail.
```

Recommended first attack:

```text
C. signed total positivity / sign-regularity.
```

Reason:

The local zero-flow route has a known obstruction: positive square-root birth is locally compatible with the same repulsive law. The signed Hankel/Jensen evidence points toward a different global mechanism.

## 11. Unbounded-Degree Boundary Refinements

### Lemma 11.1: Canonical Root Field At A Double Zero

At a real-zero Newman time, put `s=-z^2` and write

```text
F_lambda(s)=F_lambda(0)*product_j (1+s/x_j(lambda)^2).
```

If `x>0` is a double zero, `rho=-x^2`, and

```text
F_*(s)=(s-rho)^2 U(s),
```

then the regularized squared-zero field and stiffness are

```text
E_x=U'(rho)/U(rho)=sum_(j!=*) m_j/(x_j^2-x^2),
K_s=-E_x'(rho)=sum_(j!=*) m_j/(x_j^2-x^2)^2>0,
F_*'''(rho)/F_*''(rho)=3E_x.
```

In signed Newman variables, if `H_*(z)=(z-x)^2V(z)`, then

```text
B_x=V'(x)/V(x)=1/x-2xE_x,
K_x=sum_(signed y outside the pair) 1/(x-y)^2
   =1/(2x^2)+2E_x+4x^2K_s>0.
```

### Lemma 11.2: Pair Gap And Stiffness

For an already-real pair with center `c`, gap `g`, and `q=g^2`, the Newman
zero ODE gives

```text
c' -> 2B_x,
g'=4/g-2g*sum_y 1/((x_+-y)(x_--y)),
q'=8-4q*S(c,q).
```

At a finite double-zero birth time `t_*`,

```text
q(t)=8(t-t_*)-16K_x(t-t_*)^2+O((t-t_*)^3).
```

The stiffness correction is coercive, but its field component has no generic
sign: `(1+s)^2(1+s/2)` gives `E_x=1`, while `(1+s)^2(1+2s)` gives `E_x=-2`.
Both are positive-coefficient Laguerre-Polya polynomials. Thus an argument for
Xi must prove a paired-cancellation and far-tail estimate specific to its zero
measure.

### Lemma 11.3: Laguerre Scale Kernel

For the established moment normalization

```text
A_j(lambda)=j!/(2j)!*integral_0^infinity W_lambda(u)u^(2j)du,
```

the unshifted degree-`D` Jensen polynomial is

```text
J_D(w)=integral_0^infinity W_lambda(u)K_D(w,u^2)du,
K_D(w,y)={}_1F_1(-D;1/2;-wy/4)
        =D!/(1/2)_D*L_D^(-1/2)(-wy/4).
```

For every `y>0`, `K_D(w,y)` has `D` simple negative roots.

### Gate 11.4: Integration Is Not A Generic Preserver

The positive measure

```text
(9/10)delta_(1/4)+(1/10)delta_(5/2)
```

produces at `D=2`

```text
J_2(w)=1+19w/40+109w^2/1920,
disc(J_2)=-7/4800.
```

Even the log-concave exponential density produces at `D=3`

```text
J_3(w)=1+3w/2+w^2/2+w^3/20,
disc(J_3)=-1/200.
```

Hence neither positive scale mixing nor log-concavity can be used as the
missing all-degree theorem.

### Lemma 11.5: Half-Integer Gamma Benchmark

If `Y` has the Gamma law with shape `alpha` and scale `theta`, then

```text
P_(D,alpha,theta)(w)={}_2F_1(-D,alpha;1/2;-theta*w/4).
```

For `alpha=m+1/2`, `m>=0`, Euler's transformation and the Jacobi identity
show that all `D` roots are real and negative. There are `min(D,m)` simple
roots in `(-4/theta,0)` and an endpoint root `-4/theta` of multiplicity
`max(D-m,0)`.

This is an exact all-degree benchmark, not an Xi representation theorem. A
closing kernel route still needs an Xi-specific common-interlacing,
total-positive connection, or direct variation-diminishing theorem.

### Lemma 11.6: Classical Field And Fixed-Time Compactness

For the Riemann-von Mangoldt reference density

```text
rho_0(y)=log(y/(4*pi))/(4*pi), y>=a>4*pi,
```

the signed continuum field is

```text
B_rho(x)=2x*PV integral_a^infinity rho_0(y)/(x^2-y^2)dy
        =-pi/8+O_a(1/x).
```

The positive-time counting density adds `t/(16y)`, changing the field only by

```text
t/(8x)*log(sqrt(x^2-a^2)/a)=O_a(t*log(x/a)/x).
```

Implicit differentiation of the corresponding quantile law gives

```text
dx_n/dt=-pi*x_n*log(x_n/(4*pi))
        /(4*x_n*log(x_n/(4*pi))+pi*t)
        -> -pi/4,
```

consistent with `dx_n/dt=2B_n` and the field limit `-pi/8`.

For each fixed `t>0`, published positive-time asymptotics make all zeros above
a scale `exp(C/t)` real, simple, and uniquely localized. Thus, conditional on
`Lambda>0`, every multiple boundary zero is confined to
`|x|<exp(C/Lambda)`. This does not prove `Lambda<=0`: the compact region grows
without bound as `Lambda->0+`.

Macroscopic location control is insufficient. In a signed integer lattice with
double roots at `+/-n`, move the even neighboring pair `+/-(n+1)` to
`+/-(n+epsilon)`; the field at `+n` tends to `-infinity`. Moving
`+/-(n-1)` to `+/-(n-epsilon)` makes it tend to `+infinity`. Every root moves
by less than one. A closing perturbative argument therefore needs a
lambda-uniform reciprocal-gap weighted discrepancy estimate, not merely a
Riemann-von Mangoldt counting bound.

### Lemma 11.7: Local Odd-Count Reduction

At an even positive double zero `c`, remove the positive double atom and let
`mu` count the other positive roots. Then

```text
B(c)=1/c+PV integral K_c(y)dmu(y),
K_c(y)=2c/(c^2-y^2)=1/(c-y)+1/(c+y).
```

Let the reference density be the derivative of the positive-time counting
law, and let `F` be its cumulative discrepancy from `mu`. If

```text
|F(y)|<=E on [a,2c],
|F(y)|<=C*log(2+y) for y>=2c,
0<H<=c/2,
```

then Stieltjes integration by parts gives the explicit outer estimate

```text
|field discrepancy outside (c-H,c+H)|
  <=5E/H+2C*(log(4c)+1)/c.
```

The published positive-time count gives `E=O(log c)` uniformly for
`0<t<=1/2`. Taking `H=log(4c)^2` and using the continuum field yields

```text
B(c)=-pi/8+S_H(c)+O(1/log(4c)+log(4c)^3/c),
S_H(c)=sum_(0<|y-c|<H) 1/(c-y).
```

If

```text
D_c(u)=# outsider roots in [c-u,c)-# outsider roots in (c,c+u],
```

then exactly

```text
S_H=D_c(H)/H+integral_0^H D_c(u)/u^2 du.
```

Thus the unresolved field statistic is local and odd under reflection about
the collision. A global `O(log c)` counting error cannot control its
inverse-square weight near `u=0`.

Field control is still not collision exclusion. Let

```text
a^2=1+16/(8+pi),
P(z)=(z^2-1)^2(z^2-a^2),
F_tau=exp(-tau*d_z^2)P.
```

This is an exact even solution of `F_tau=-F_zz`. At the double zero `z=1`,

```text
B(1)=-pi/8,
z_+/-=1+/-sqrt(2*tau)-pi*tau/4+O(tau^(3/2)).
```

It has the classical field and drift while still producing a real pair for
positive `tau`. Any closing Newman route therefore needs both Xi-specific
local odd-count control and an additional global mechanism forbidding birth.

### Lemma 11.8: Boundary Energy Direction Gate

For a real pair born at `t_*`, let `g=x_+-x_-`, `q=g^2`, and

```text
K=sum_(outsider y) 1/(c-y)^2.
```

The zero ODE gives

```text
q'=8-4q*S,
q(t)=8*(t-t_*)-16*K*(t-t_*)^2+O((t-t_*)^3).
```

In the exact even model of Lemma 11.7,

```text
K=(pi^2+24*pi+160)/64>0,
q(tau)=8*tau-16*K*tau^2+C3*tau^3+O(tau^(7/2)),
C3=(pi^4+48*pi^3+800*pi^2+5632*pi+14848)/32>0.
```

Thus the classical field, positive stiffness, negative quadratic gap
curvature, and positive cubic gap jet are all compatible with positive birth.

For a fixed reference gap `Delta>0`, the Rodgers-Tao renormalized interaction
uses

```text
V(r)=r^(-2)+2r-3,
E_pair=Delta^(-2)*V(g/Delta)
      =1/g^2+2g/Delta^3-3/Delta^2.
```

Since `V(1)=V'(1)=0` and `V''(r)=6/r^4>0`, this interaction is nonnegative.
At a collision,

```text
E_pair(t)=1/(8*(t-t_*))+O(1),
dE_pair/dt=-1/(8*(t-t_*)^2)+O((t-t_*)^(-1/2)),
integral_(t_*)^(t_*+epsilon) E_pair(t)dt=+infinity.
```

Hence nonnegative forward energy decrease does not exclude birth: the pair
relaxes from an infinite boundary trace. On the other hand, any Xi-specific
renormalized energy containing this interaction with a weight bounded below
would exclude collision if it were locally time-integrable all the way down
to `t_*`.

The published Rodgers-Tao integrated-energy theorem is not such a boundary
estimate. It is proved under the contradiction assumption `Lambda<0` and
integrates from `Lambda/4` to `0`, a distance `3|Lambda|/4` after the boundary.
It controls relaxation in an already-real interval and uses zeta information
at `t=0`; it is unavailable in the opposite hypothetical regime `Lambda>0`.
The open energy handoff is therefore a noncircular Xi boundary-trace theorem
on `(Lambda,Lambda+epsilon)`, not an extension of the existing interior
estimate by terminology alone.

### Lemma 11.9: Positive Boundary Attainment

Two published positive-time estimates remove the high-index escape loophole
under the contradiction hypothesis `Lambda>0`:

```text
all zeros of H_t satisfy |Im z|<=sqrt(1-2t)<=1 for 0<t<=1/2;
if Re z>=exp(C/t) and H_t(z)=0, then Im z=0,
```

where `C` is absolute. Choose `t_n<Lambda` with `t_n->Lambda` and
`t_n>=Lambda/2`. Since `t_n<Lambda`, each `H_(t_n)` has a nonreal zero `z_n`.
The two estimates give the fixed compact bound

```text
|Re z_n|<X_Lambda=exp(2C/Lambda),
|Im z_n|<=1.
```

After passage to a subsequence, `z_n->z_*`. The Fourier integral gives local
uniform convergence `H_(t_n)->H_Lambda`, so `H_Lambda(z_*)=0`; this zero is
real because all zeros at `t=Lambda` are real. Both `z_n` and its distinct
complex conjugate converge to `z_*`. Rouche zero counting therefore gives

```text
mult_(z_*) H_Lambda>=2.
```

Consequently,

```text
Lambda>0 => H_Lambda has a finite real multiple zero c,
             |c|<=exp(2C/Lambda).
```

Together with simplicity for every `t>Lambda` and the published bound
`Lambda<=1/5`, this yields the sharper exact equivalence

```text
Lambda<=0
iff H_t has only simple zeros for every 0<t<=1/5.
```

If the boundary zero has multiplicity `m`, the published local heat theorem
splits it as

```text
x_a(Lambda+tau)=c+sqrt(2tau)*lambda_a+O(tau),
```

where the `lambda_a` are the simple roots of the probabilists' Hermite
polynomial `He_m`. Its ODE gives

```text
sum_(a!=b) 1/(lambda_a-lambda_b)^2=m(m-1)/4.
```

Hence the ordered cluster energy satisfies

```text
sum_(a!=b) 1/(x_a-x_b)^2
  =m(m-1)/(8tau)+O(tau^(-1/2)),
```

and every multiplicity `m>=2` has logarithmically nonintegrable boundary
energy. Thus the energy criterion of Lemma 11.8 covers every hypothetical
positive Newman boundary, not only an attained double zero. The remaining
theorem is Xi-specific positive-time simplicity or finite-truncation endpoint
energy integrability; neither published interior relaxation nor compactness
proves that theorem.

### Lemma 11.10: Strict Laguerre Correlation Criterion

Platt-Trudgian Corollary 2 gives the published bound `Lambda<=1/5`.
Consequently, if `Lambda>0`, then the attained boundary lies in `(0,1/5]`.
It therefore suffices to work on this smaller interval. For `0<t<=1/5`, set

```text
phi_t(u)=exp(t*u^2)*Phi(u),
H_t(x)=integral_0^infinity phi_t(u)cos(xu)du,
L_t(x)=H_t'(x)^2-H_t(x)H_t''(x).
```

Lemma 11.9 and the simple real-zero factorization give the exact equivalence

```text
Lambda<=0
iff L_t(x)>0 for every real x and every 0<t<=1/5.
```

Indeed, the forward direction has

```text
L_t/H_t^2=sum_j 1/(x-r_j)^2
```

away from the simple roots, while at a root `L_t(r_j)=H_t'(r_j)^2>0`.
Conversely, `Lambda>0` forces a multiple real zero `c` of `H_Lambda`, where
`L_Lambda(c)=0`.

Define the first correlation kernel

```text
K_(1,t)(v)=integral_R phi_t(s+v)phi_t(s-v)s^2 ds.
```

The midpoint/difference change of variables gives the normalization-exact
identity

```text
L_t(x)=integral_R K_(1,t)(v)cos(2xv)dv,
Fourier[K_(1,t)](xi)=L_t(xi/2).
```

Also `L_t(0)>0`, since `H_t(0)>0`, `H_t'(0)=0`, and `H_t''(0)<0`.
Wiener's `L^1` Tauberian theorem therefore yields another exact equivalence:

```text
Lambda<=0
iff for every 0<t<=1/5, the translations of K_(1,t)
    are dense in L^1(R).
```

Ordinary positive definiteness is insufficient. At a hypothetical boundary,
`H_Lambda` is Laguerre-Polya and `K_(1,Lambda)` is positive definite, but its
Fourier transform still vanishes at `2c`.

Strict log-concavity does not repair this generically. Let

```text
T_a(y)=(a-|y|)_+,
G_sigma(x)=exp(-x^2/(2sigma^2)),
K_(a,sigma)=T_a*G_sigma.
```

If `p_x(y)` is proportional to `T_a(y)G_sigma(x-y)`, then

```text
(log K)''=Var_x(Y)/sigma^4-1/sigma^2
          <=(a^2-sigma^2)/sigma^4<0       (sigma>a).
```

Thus `K` is smooth, positive, even, and strictly log-concave. It is also
positive definite, because

```text
Fourier[K](xi)=4*sqrt(2*pi)*sigma*exp(-sigma^2*xi^2/2)
               *sin(a*xi/2)^2/xi^2>=0.
```

Nevertheless this transform has a double zero at every
`xi=2*pi*n/a`, `n!=0`, so the translations are not dense. The live theorem
must prove Xi-specific zero-freeness of the first correlation transform using
structure beyond generic log-concavity and positive definiteness.

### Lemma 11.11: Correlation Hierarchy and Gaussian-Mixture Gate

For every integer `n>=0`, define

```text
K_(n,t)(v)=integral_R phi_t(s+v)phi_t(s-v)s^(2n) ds,
F_(n,t)(xi)=Fourier[K_(n,t)](xi).
```

Since

```text
(s+v)^2+(s-v)^2=2(s^2+v^2),
```

differentiation under the super-exponentially convergent integral and Fourier
transformation give the exact hierarchy

```text
partial_t K_(n,t)=2v^2 K_(n,t)+2K_(n+1,t),
partial_t F_(n,t)=-2 partial_xi^2 F_(n,t)+2F_(n+1,t).
```

If the generalized Laguerre expressions are normalized by

```text
|H_t(x+iy)|^2=sum_(n>=0)L_(n,t)(x)y^(2n),
```

then the correlation formula is

```text
L_(n,t)(x)=2^(2n-1)/(2n)! F_(n,t)(2x).
```

Consequently, a real root `c` of multiplicity `m` has the hierarchy signature

```text
F_(n,t)(2c)=0                                      for n<m,
F_(m,t)(2c)=(2m)!/2^(2m-1)*(H_t^(m)(c)/m!)^2>0.
```

For a double root this specializes to

```text
F_(1,t)(2c)=partial_xi F_(1,t)(2c)=0,
partial_xi^2 F_(1,t)(2c)=H_t''(c)^2/4,
F_(2,t)(2c)=3H_t''(c)^2/4=3 partial_xi^2 F_(1,t)(2c),
partial_t F_(1,t)(2c)=H_t''(c)^2.
```

Thus the exact hierarchy is compatible with a nonnegative quadratic boundary
touch; it does not itself give a favorable maximum-principle contradiction.

There is a strong sufficient criterion. Set

```text
g_t(r)=K_(1,t)(sqrt(r)).
```

If `g_t` were nonzero and completely monotone, Bernstein-Schoenberg would give

```text
g_t(r)=integral_[0,infinity) exp(-a r)dmu_t(a),
Fourier[K_(1,t)](xi)
 =sqrt(pi)*integral_(a>0)a^(-1/2)exp(-xi^2/(4a))dmu_t(a)>0.
```

This would prove the strict criterion in Lemma 11.10. The hypothesis is,
however, exactly impossible for the Xi correlation. For `u>=0`, the defining
Phi series gives

```text
0<Phi(u)<=C_Phi exp(9u-pi exp(4u)),
C_Phi=2pi^2 sum_(n>=1)n^4 exp(-pi(n^2-1))<infinity.
```

For `v>=0`, use

```text
|s+v|+|s-v|=2max(|s|,v),
exp(4|s+v|)+exp(4|s-v|)>=exp(4v)+exp(4|s|).
```

Then, uniformly for `0<=t<=1/2`,

```text
0<K_(1,t)(v)
 <=C_Phi^2 M exp(2t v^2+18v-pi exp(4v)),
M=integral_R s^2 exp(s^2+18|s|-pi exp(4|s|))ds<infinity.
```

Hence

```text
lim_(r->infinity) -log(g_t(r))/r=+infinity.
```

A nonzero completely monotone function has a positive Laplace measure and
therefore satisfies `g(r)>=mu([0,A])exp(-Ar)` for some finite `A`. The two
statements contradict each other, so `g_t` is not completely monotone for any
`0<=t<=1/2`.

The local numerical behavior is deliberately deceptive. Composite
Gauss-Legendre quadrature and an independent 55-digit mpmath run find
alternating signs `(-1)^n g_t^(n)(0)>0` through `n=8`, but already

```text
g_0(0)g_0''(0)-g_0'(0)^2
 =-0.0000078630150226764958198468357098...,
g_(1/2)(0)g_(1/2)''(0)-g_(1/2)'(0)^2
 =-0.0000079006312495514071411797904742....
```

Thus the first log-convexity minor required of a Laplace mixture is negative.
These values are diagnostics; the super-Gaussian tail is the exact global
obstruction.

Two related generic routes are also closed. A differentiable even kernel has
`K'(0)=0`, so it cannot be both nonconstant decreasing and convex on the
positive half-line; the classical Polya convexity criterion does not apply.
Also, Schoenberg's classification represents a nondegenerate PF-infinity
function through Gaussian and one-sided exponential convolution factors. An
even such function has a Gaussian lower tail if the Gaussian factor is
present, and an exponential lower tail otherwise. The two-sided
super-Gaussian bound therefore rules out direct PF-infinity membership of
`K_(1,t)`.

The remaining correlation route must be tail-compatible: either construct an
Xi/theta-summand spectral-square decomposition with no common real spectral
zero, or prove a relative hierarchy coercivity estimate that excludes the
universal contact `F_1=F_1'=0`, `F_2=3F_1''>0`. This is an open handoff, not a
proof of strict Fourier positivity, Wiener density, RH, or `Lambda<=0`.

### Lemma 11.11A: Positive-Time Strong Log-Concavity

Let `ell(u)=log Phi(u)`. Csordas and Varga proved that

```text
q(r)=ell(sqrt(r))
```

is strictly concave for `r>0`. For `u>0`,

```text
q'(u^2)=ell'(u)/(2u),
q''(u^2)=(u ell''(u)-ell'(u))/(4u^3)<0.
```

Since `Phi` is even, `lim_(u->0+) ell'(u)/u=ell''(0)`. Strict decrease of
`q'` and the displayed second identity therefore give, successively,

```text
ell'(u)/u<ell''(0),
ell''(u)<ell'(u)/u<ell''(0).
```

Evenness extends the curvature ceiling to every real `u`. With

```text
kappa=-ell''(0)=-Phi''(0)/Phi(0),
```

the exact conclusion is

```text
(log Phi)''(u)<=-kappa                 (u in R).
```

A 256-bit Arb sum through theta index eight, with explicit geometric bounds
for every omitted term, certifies

```text
Phi(0)   =0.4466969004671234440869846670547091...,
Phi''(0)=-33.46100154940651373358024870154595...,
kappa    =74.90761971801328467366520780869287...,
kappa/2-1/5>37.2538098590066423368326039043.
```

Consequently, for `phi_t(u)=exp(tu^2)Phi(u)`,

```text
(log phi_t)''(u)<=-(kappa-2t)<0       (0<=t<=1/5).
```

Thus every kernel in the exact positive-time target window is positive, even,
strictly decreasing on the positive half-line, super-exponentially decaying,
and uniformly strongly log-concave.

This curvature propagates through the entire correlation hierarchy. On
`s>0`, the logarithmic Hessian of

```text
phi_t(s+v)phi_t(s-v)s^(2n)
```

has quadratic form

```text
a(p+q)^2+b(p-q)^2-2n p^2/s^2,
a=(log phi_t)''(s+v), b=(log phi_t)''(s-v).
```

It is bounded above by `-2(kappa-2t)(p^2+q^2)`. Strong Prekopa
marginalization then yields

```text
(log K_(n,t))''(v)<=-2(kappa-2t)      (n>=0, 0<=t<=1/5).
```

This is a substantial shape theorem but not a Fourier zero theorem. Indeed,

```text
Fourier[K_(0,t)](xi)=2H_t(xi/2)^2.
```

At `t=0`, Hardy's theorem supplies real zeros of `H_0`, so this transform has
double real zeros although `K_(0,0)` has the full strong-log-concavity,
positive-definiteness, admissibility, and Xi-tail package. The artifact and
independent checker are

```text
outputs/jensen_window_pf_newman_positive_time_strong_logconcavity_gate.md
python work/rh_compute/scripts/check_jensen_window_pf_newman_positive_time_strong_logconcavity_gate.py
```

### Lemma 11.11B: Weighted Strong-Log-Concavity Countermodel

Even passing from `K_0` to the first `s^2`-weighted correlation does not make
generic strong log-concavity sufficient. For `delta,epsilon>=0`, define

```text
phi_(delta,epsilon)(u)
 =exp(-u^2-delta*u^4-epsilon(cosh(4u)-1))(1+u^4/10).
```

Put `z=u^2/sqrt(10)`. Direct differentiation and
`z/(1+z^2)<=1/2` give

```text
(log phi_(delta,epsilon))''
 =-2-12delta*u^2-16epsilon*cosh(4u)
   +4/sqrt(10)*z(3-z^2)/(1+z^2)^2
 <=-(2-6/sqrt(10))-12delta*u^2-16epsilon*cosh(4u)<0.
```

For `delta>0` or `epsilon>0`, this is a smooth positive even strictly
decreasing admissible kernel. When `epsilon>0`, it has theta-type
double-exponential decay.

At `delta=1/10`, `epsilon=1/1000`, the stronger root-variable shape property
also holds. Direct differentiation gives

```text
d^2/dr^2[-r-r^2/10+log(1+r^2/10)]
 =-r^2(r^2+30)/(5(r^2+10)^2)<0,
d^2/dr^2 cosh(4sqrt(r))
 =[4sqrt(r)cosh(4sqrt(r))-sinh(4sqrt(r))]/r^(3/2)>0.
```

For the last sign, put `y=sqrt(r)`; the numerator vanishes at zero and its
derivative is `16y sinh(4y)>0`. Subtracting the theta term therefore makes the
root profile still more concave. The explicit countermodel matches ordinary
uniform strong log-concavity, strict concavity of `log(phi(sqrt(r)))`, and the
Xi theta tail class. At the endpoint `delta=epsilon=0`, Gaussian
differentiation evaluates its full Fourier transform exactly:

```text
F_0(x)=sqrt(pi)exp(-x^2/4)(x^4-12x^2+172)/160.
```

Its first Laguerre expression is

```text
L_1[F_0](x)
 =pi exp(-x^2/2)
  *(x^8-16x^6+440x^4-7680x^2+37840)/51200,
L_1[F_0](3)=-743pi exp(-9/2)/51200<0.
```

A 256-bit Acb integral on `[0,6]`, with explicit upper incomplete-gamma
enclosures for all omitted tails, directly certifies at `delta=1/10`,
`epsilon=1/1000`, and `x=21/5` that

```text
F       =[ 0.010079077726118413514324903055736354411454698489725216
           +/- 3.24e-55],
F'      =[-0.02883692696821660760997548984142684809891905250938752
           +/- 3.43e-54],
F''     =[ 0.1036847681860763492844176464243140394504271561765499
           +/- 3.77e-53],
L_1[F]  =[-0.000213478480591774964507646766853475076598648832979809
           +/- 5.09e-55]<0.
```

The exact first-correlation identity now gives

```text
L_1[F_(delta,epsilon)](x)=4 Fourier[K_(1,delta,epsilon)](2x),
Fourier[K_(1,1/10,1/1000)](42/5)<0.
```

Thus an admissible input with Xi-style root-variable concavity, uniform strong
log-concavity, and theta-tail decay can have a first weighted correlation that
is not even positive definite. This closes generic
Xi-style-shape-and-tail-plus-weighting promotion. It does not decide the Xi
correlation, whose surviving route must use arithmetic theta structure, a
modularly grouped weighted square, or relative hierarchy coercivity. The exact
artifact is

```text
outputs/jensen_window_pf_newman_weighted_strong_logconcavity_countermodel_gate.md
python work/rh_compute/scripts/check_jensen_window_pf_newman_weighted_strong_logconcavity_countermodel_gate.py
```

### Lemma 11.11C: Strict-Laguerre Monotonicity Counter-Gate

The following stronger condition is sufficient but false for Xi. Put

```text
L_t(x)=H_t'(x)^2-H_t(x)H_t''(x),
M_t(x)=-L_t'(x)=H_t(x)H_t'''(x)-H_t'(x)H_t''(x).
```

Because every Fourier jet tends to zero at positive infinity, strict
`M_t(x)>0` for all `x>0` gives

```text
L_t(x)=integral_x^infinity M_t(y)dy>0.
```

Uniformly for `0<t<=1/5`, this would imply `Lambda<=0` by Lemma 11.10. A
multiple boundary root `c>0` would instead have `H_t(c)=H_t'(c)=0` and hence
`M_t(c)=0`, so the strict margin directly excludes the universal contact.

Differentiating the first-correlation transform gives the exact sine form

```text
M_t(x)=4 integral_0^infinity v*K_(1,t)(v)*sin(2xv)dv.
```

Define the positive decreasing tail

```text
Q_t(s)=integral_s^infinity v*K_(1,t)(v)dv.
```

Fubini and `sin(2xv)=2x integral_0^v cos(2xs)ds` give

```text
M_t(x)=8x integral_0^infinity Q_t(s)cos(2xs)ds.
```

This does not fall under the classical positive-decreasing-convex criterion.
Indeed `Q_t'=-sK_(1,t)` and `Q_t''=-K_(1,t)-sK_(1,t)'`, so
`Q_t''(0)=-K_(1,t)(0)<0`. More sharply, if `ell=log K_(1,t)`, then
`h(s)=1+s ell'(s)` decreases strictly from `1` to negative infinity because
`h'=ell'+s ell''<0` and strong log-concavity gives `ell'(s)<=-cs`. Hence
`Q_t''=-K_(1,t)h` changes sign exactly once, from negative to positive. The
tail is concave then convex, and its transform sign remains Xi-specific.

The elementary fixed-interior complex-shift criterion is also unavailable.
At every regular shift `p>0`, evenness and real analyticity imply that
`Q_t(ip)` is real. Hence `-Im Q_t(s+ip)` starts at zero and cannot be both
nonnegative and nonincreasing unless it vanishes identically.

For the endpoint-subtracted theta primitive `A_t=D_t[C_t]=8H_t-1/2`, the same
target is

```text
64M_t(x)=(A_t(x)+1/2)A_t'''(x)-A_t'(x)A_t''(x)>0.
```

A reproducible scout finds `L_t'(x)<0` on 12000 dense rows with
`x=j/50<=40`, six times through `t=1/2`, and on 20 selected 75-digit rows
through `x=200`; the `t=0` rows are independently checked from Xi derivatives.
Those grids miss the close-pair stress region. At the explicit rational point
`x=1401016343/100000`, a 256-bit Arb Xi jet through third order certifies

```text
M_0(x)/A_0(x)^2=[-3.24186746975685391605462931182063307428944779983857697073123186e-6 +/- 5.42e-69],
M_0(x)/L_0(x)=[-0.0001463088540165360168258875503621433728403782920777681331619380172 +/- 5.15e-68].
```

Thus `M_0(x)<0`. The theta tail makes
`|u|^3 exp(u^2/5)Phi(u)` integrable, so dominated convergence makes the
heat-flow jet through third order continuous in `t` on `[0,1/5]`. The strict
negative sign therefore persists for all sufficiently small positive `t`.
This rejects the condition
`M_t(x)>0` for every `x>0` and every `0<t<=1/5`.

The mechanism is exact. If
`F(m+y)=(y^2-a^2)G(m+y)`, then

```text
M[F](m)=-4a^2 G(m)G'(m)+a^4(G(m)G'''(m)-G'(m)G''(m)).
```

Unlike first-Laguerre positivity, the stronger sign depends at leading order
on the slope of the regular factor near a close pair. The global sine-transform
and coupled boundary-contour branch is therefore retired. This does not reject
`L_t>0`; the surviving Newman endgames are direct strict Laguerre positivity
or corrected `C1` double-zero transversality. The artifact is

```text
outputs/jensen_window_pf_newman_strict_laguerre_monotonicity_scout.md
python work/rh_compute/scripts/check_jensen_window_pf_newman_strict_laguerre_monotonicity_scout.py
```

### Lemma 11.12: Theta-Summand Spectral-Square Gate

For `u>=0`, write

```text
phi_n(u)=(2pi^2 n^4 exp(9u)-3pi n^2 exp(5u))
         *exp(-pi n^2 exp(4u)),
Phi(u)=sum_(n>=1)phi_n(u).
```

Lemma 8.16 gives the arithmetic shift

```text
phi_n(u)=n^(-1/2)phi_1(u+(log n)/2).
```

There is also an exact differential primitive. If

```text
h(u)=exp(u-pi exp(4u)),
R(u)=sum_(n>=1)exp(u-pi n^2 exp(4u)),
```

then

```text
phi_1=(h''-h)/8,
Phi=(R''-R)/8.
```

Let `theta(a)=sum_(n in Z)exp(-pi n^2 a)`. Differentiating
`theta(a)=a^(-1/2)theta(1/a)` at `a=1` gives

```text
theta'(1)=-theta(1)/4,
R'(0)=-1/2.
```

Thus, for

```text
C(x)=integral_0^infinity R(u)cos(xu)du,
```

two integrations by parts give the exact undeformed transform identity

```text
H_0(x)=1/16-(1+x^2)C(x)/8.
```

This is a cancellation formula rather than a positive square: the boundary
constant cancels the complete algebraic high-frequency expansion of the
second term.

The endpoint-subtracted formula extends exactly to every target time. Define

```text
C_t(x)=integral_0^infinity exp(tu^2)R(u)cos(xu)du,
D_t=-4t^2 partial_x^2+4tx partial_x+(2t-1-x^2).
```

For `w(u)=exp(tu^2)cos(xu)`, two integrations by parts use `w'(0)=0` and
`R'(0)=-1/2` to give

```text
H_t(x)=1/16+D_t[C_t](x)/8
      =1/16+(-4t^2 C_t''+4tx C_t'+(2t-1-x^2)C_t)/8.
```

Writing `A_t=D_t[C_t]`, the strict Laguerre expression becomes

```text
L_t(x)=(A_t'(x)^2-(A_t(x)+1/2)A_t''(x))/64.
```

A multiple real zero is exactly a contact
`A_t(c)=-1/2`, `A_t'(c)=0`. The infinite endpoint cancellation is therefore
built into one explicit all-time curvature target.

The corresponding bilateral profile audit is also exact. With
`s=(1+iz)/2`,

```text
Qhat(z)=integral_R phi_1(u)exp(izu)du
       =-(1+z^2)/32*pi^(-(1+iz)/4)*Gamma((1+iz)/4),
sum_(n>=1)n^(-1/2)exp(-iz log(n)/2)=zeta(s),
Qhat(z)zeta(s)=xi(s)/4.
```

The sum converges absolutely first in `Im z<-1`, and the last identity then
continues analytically. Hence the direct shifted-profile spectral assembly
reconstructs the completed zeta function itself; it is not a new zero-free
factorization.

Termwise Laguerre positivity fails more strongly. Define

```text
H_(n,t)(x)=integral_0^infinity exp(tu^2)phi_n(u)cos(xu)du,
L[f]=f'^2-ff'',
B[f,g]=2f'g'-fg''-gf''.
```

Then

```text
L[sum_n H_n]=sum_n L[H_n]+sum_(m<n)B[H_m,H_n].
```

Neither the self nor cross terms have a fixed sign. More decisively, every
finite theta block fails globally. For

```text
f_(N,t)(u)=exp(tu^2)sum_(n=1)^N phi_n(u),
```

direct differentiation gives

```text
phi_n'(0)=pi n^2 exp(-pi n^2)
          *(-8pi^2 n^4+30pi n^2-15).
```

The roots of the final quadratic in `q=pi n^2` are
`(15+-sqrt(105))/8`; hence `phi_1'(0)>0` and `phi_n'(0)<0` for `n>=2`.
Since the infinite kernel is even, `sum_n phi_n'(0)=Phi'(0)=0`. Therefore

```text
A_N=f_(N,t)'(0)=sum_(n=1)^N phi_n'(0)
   =-sum_(n>N)phi_n'(0)>0,
```

independently of `t`. Repeated integration by parts yields, for every fixed
finite `N` and real `t`,

```text
H_(N,t)(x)=-A_N/x^2+O_(N,t)(x^-4),
H_(N,t)'(x)=2A_N/x^3+O(x^-5),
H_(N,t)''(x)=-6A_N/x^4+O(x^-6),
L[H_(N,t)](x)=-2A_N^2/x^6+O_(N,t)(x^-8)<0
```

for all sufficiently large `x`. Thus no finite positive grouping of theta
summands can establish the global first-Laguerre criterion. Its failure moves
to higher frequency as `N` grows because `A_N->0`, but disappears only after
the exact infinite modular cancellation.

Independent 60-digit quadrature locates the obstruction at moderate
frequencies:

```text
B[H_(1,0),H_(2,0)](16.75)
 =-5.455807042295599276470773299858...e-9,
L[H_(1,0)](39.5)
 =-1.435395530583386003271352965275...e-11,
L[H_0](39.5)
 =+8.112479065720435571302196294182...e-12.
```

At `t=1/2` the same signs persist. These values are diagnostics; the finite-N
asymptotic is exact.

The live theorem must therefore operate after infinite endpoint cancellation.
The exact deformed formula reduces it to

```text
A_t'(x)^2-(A_t(x)+1/2)A_t''(x)>0
```

for every real `x` and `0<t<=1/2`, proved from the positive theta primitive
`R` and its modular identity. Finite summand squares, pairwise-positive cross
terms, and direct Mellin reconstruction are closed routes. No noncircular
curvature estimate of this strength is currently proved. This is not a proof
of RH or `Lambda<=0`.

### Lemma 11.13: Gasper Fake-Xi Remainder and Convolution Gate

Define the scaled Polya fake-Xi kernel and transform by

```text
Psi(u)=4pi^2 exp(-2pi cosh(4u))cosh(9u),
P_t(x)=integral_0^infinity exp(tu^2)Psi(u)cos(xu)du.
```

In the variable used here, Gasper's model has the exact normalization

```text
P_0(x)=Xi_star(x/2)/8
      =pi^2/2*(K_(9/4+ix/4)(2pi)+K_(9/4-ix/4)(2pi)).
```

Gasper's integral-of-squares theorem proves that `P_0` has only real zeros;
the classical heat universal-factor theorem preserves this for `t>=0`.
This supplies a genuine comparison model, but not a result about the true
transform `H_t`.

Let

```text
L[f]=f'^2-ff'',
B[f,g]=2f'g'-fg''-gf'',
E_alpha=H-alpha P.
```

Then the exact remainder decomposition is

```text
L[H]=alpha^2 L[P]+B[alpha P,E_alpha]+L[E_alpha],
B[alpha P,E_alpha]=alpha B[P,H]-2alpha^2 L[P],
L[E_alpha]=L[H]-alpha B[P,H]+alpha^2 L[P].
```

The elementary absolute-margin route would require

```text
R_alpha=(|B[alpha P,E_alpha]|+|L[E_alpha]|)/(alpha^2 L[P])<1.
```

At a fixed point suppose `L[P]>0`, `L[H]>0`, `B[P,H]>0`, and set

```text
rho=B[P,H]^2/(L[P]L[H]).
```

If `0<rho<2`, direct completion of the square gives

```text
inf_(alpha>0) R_alpha=3-rho,
alpha_star=L[H]/B[P,H],
R_alpha-(3-rho)
 =(alpha B[P,H]-L[H])^2/(alpha^2 L[P]L[H])
```

in the minimizing sign branch; the other branch is no smaller. Arb midpoint
quadrature on `[0,1]`, with interval second-derivative remainders and explicit
`n>=9` and `u>=1` tails, certifies the sign chamber and a minimum strictly
above one at `x=25` for both endpoint times. Independent 80-digit quadrature
gives

```text
t=0:   rho=0.4147964651353511451...,
       inf R_alpha=2.5852035348646488549...,
t=1/2: rho=0.4307610900175452454...,
       inf R_alpha=2.5692389099824547546....
```

The displayed digits are independent cross-checks of the Arb balls. The
certificate rejects every positive scalar normalization of this one-block
triangle estimate; it does not reject a sign-aware estimate because the mixed
and residual terms can cancel.

The direct positive-convolution transfer fails exactly. For `u>=0`, put

```text
M(u)=Phi(u)/Psi(u).
```

The origin values satisfy

```text
Phi(0)>pi(2pi-3)exp(-pi)>4pi^2 exp(-2pi)=Psi(0),
```

so `M(0)>1`. On the other hand,

```text
Psi(u)=2pi^2 exp(9u-pi exp(4u)-pi exp(-4u))(1+exp(-18u)),
phi_1(u)=2pi^2 exp(9u-pi exp(4u))
         *(1-3exp(-4u)/(2pi)),
```

and every `n>=2` theta term is superexponentially smaller than the first.
Therefore

```text
lim_(u->infinity) M(u)=1.
```

If a finite positive measure represented this ratio as

```text
M(u)=integral_R cosh(us)dmu(s),
```

then a finite limit at infinity would force `mu(R\{0})=0` by monotone
convergence, making `M` constant. This contradicts `M(0)>1` and `M(u)->1`.
Thus the direct Cardon/Polya positive-cosh multiplier from fake Xi to Xi does
not exist.

The surviving Gasper route must preserve signs or use several
special-function blocks with a proved common real-zero/interlacing structure.
The next admissible target is a multi-block or sign-aware mixed-Laguerre
estimate performed after the full modular endpoint cancellation. A scalar
absolute budget and a one-factor positive-cosh transfer are closed routes.
This is not strict Laguerre positivity, RH, or `Lambda<=0`.

### Lemma 11.14: Positive Gasper Residual Two-Block Obstruction

The first fake-Xi subtraction has a positive kernel residual. Put

```text
Psi_9(u)=4pi^2 exp(-2pi cosh(4u))cosh(9u),
q=exp(-4u),
a=3/(2pi).
```

For `u>=0`,

```text
phi_1(u)/Psi_9(u)=exp(pi q)(1-aq)/(1+q^(9/2)).
```

Since `0<q<=1`, `pi>3`, and `a<1/2`,

```text
exp(pi q)(1-aq)
 >(1+pi q)(1-aq)
 =1+(pi-a)q-(3/2)q^2
 >1+q
 >=1+q^(9/2).
```

All higher theta summands are also positive on the positive half-line.
Consequently

```text
Phi(u)>phi_1(u)>Psi_9(u)
```

for every finite `u>=0`. Thus positivity of the residual kernel is not the
missing theorem.

Consider the full nonnegative 9/5 two-block family

```text
Psi_beta(u)=4pi^2 exp(-2pi cosh(4u))
            *(cosh(9u)+beta cosh(5u)),
P_beta=P_(9/4)+beta P_(5/4),
R_beta=H_0-P_beta,
beta>=0.
```

Its residual kernel has the exact leading tail

```text
Phi(u)-Psi_beta(u)
 ~2pi^2*(pi-3/(2pi)-beta)*exp(5u-pi exp(4u)).
```

Therefore a residual kernel that remains nonnegative in the tail must satisfy

```text
0<=beta<=beta_tail,
beta_tail=pi-3/(2pi).
```

Write `E_0=H_0-P_(9/4)` and `Q=P_(5/4)`. The first Laguerre expression is the
exact quadratic

```text
L[R_beta]=L[E_0]-beta B[E_0,Q]+beta^2 L[Q].
```

Acb evaluation of the exact xi and Bessel formulas, combined with central
differences and Cauchy derivative bounds, certifies this quadratic uniformly
in `beta`. The derivative enclosure uses step `h=10^-6`, Cauchy radius `0.04`,
and an Acb square of radius `0.1`; if `M` bounds the function on the square,
the added first- and second-derivative errors are respectively
`h^2 M/r^3` and `2h^2 M/r^4`.

At `x=66`, Arb gives

```text
c0=[-1.483e-21 +/- 5.70e-25],
c1=[-6.213e-21 +/- 4.82e-25],
c2=[+2.8475e-21 +/- 3.10e-26].
```

The quadratic is strictly convex and negative at both endpoints of
`0<=beta<=23/10`. At `x=50`, the corresponding coefficients are

```text
c0=[+1.411e-15 +/- 1.31e-19],
c1=[-1.4587e-15 +/- 6.49e-20],
c2=[+3.4559e-16 +/- 4.59e-21].
```

This quadratic is strictly convex and negative at both endpoints of
`23/10<=beta<=beta_tail`. Convexity therefore proves

```text
L[R_beta](66)<0  for 0<=beta<=23/10,
L[R_beta](50)<0  for 23/10<=beta<=pi-3/(2pi).
```

For larger `beta`, the residual kernel is eventually negative by the tail
formula. Hence no nonnegative `beta` can produce both a globally nonnegative
residual kernel and a residual transform satisfying even the first Laguerre
inequality on the whole real axis. In particular, Xi cannot be proved by
splitting off a positive 9/5 Gasper block and treating the positive remainder
as a second independent Laguerre-Polya block.

The tail-matched base also fails the standard universal-multiplier route. If

```text
U_beta(z)=cosh(9z)+beta cosh(5z)
         =cosh(z) Q_beta(cosh(z)^2),
```

then

```text
Q_beta(y)=256y^4-576y^3+(432+16beta)y^2
          -(120+20beta)y+9+5beta,
disc Q_beta=-2^24*(20beta^5-64beta^4+184beta^3
                   -432beta^2+972beta-729).
```

The final quintic is strictly increasing and positive for `beta>=11/10`;
therefore the discriminant is negative and `U_beta` has off-imaginary zeros.
Since `beta_tail>5/2`, the tail-matched multiplier is outside the standard
Laguerre-Polya universal-factor class.

The live Gasper handoff is now necessarily signed and coupled: at least one
block must be signed, or the mixed terms must enter one common square identity
after the full modular cancellation. Positive one-block domination and a
positive 9/5 block plus independently LP residual are closed architectures.
This is not strict Laguerre positivity, RH, or `Lambda<=0`.

### Lemma 11.15: Classical Three-Block Residual Obstruction

In the corpus normalization, set

```text
P_a(x)=pi^2/2*(K_(a+ix/4)(2pi)+K_(a-ix/4)(2pi)),
a=3/(2pi),
b=pi-a.
```

The two classical comparison transforms are

```text
B_P2=P_(9/4)-a P_(5/4),
B_dB=P_(9/4)+b P_(5/4)+P_(1/4).
```

Polya and de Bruijn proved, respectively, that these transforms have only
real zeros. Those established source facts are not re-proved here. What is
new is an exhaustive audit of their residuals and of the full positive
three-block family.

Put `q=exp(-4u)` and `r=sqrt(q)`. After division by
`2pi^2 exp(9u-pi exp(4u))` and multiplication by `exp(pi q)`, the general
block is

```text
S_(beta,gamma)(q)
 =1+q^(9/2)+beta q(1+q^(5/2))+gamma q^2(1+q^(1/2)).
```

For the signed P2 parameters `beta=-a`, `gamma=0`, the elementary bound
`exp(pi q)>1+pi q` gives

```text
exp(pi q)(1-aq)-S_P2(q)
 >q*(pi-3q/2-q^(7/2)+a q^(5/2))
 >=q*(pi-5/2)>0.
```

Every higher theta summand is positive on `u>=0`, hence

```text
Phi(u)-Psi_P2(u)>0
```

for every finite `u>=0`.

For the de Bruijn parameters `beta=b`, `gamma=1`, use the fourth Taylor
lower bound for `exp(pi q)`. The resulting difference is `q^2 p(r)`, where

```text
p(r)=pi^2/2-5/2-r+(-3pi/4+pi^3/6)r^2
     +(3-2pi^2)r^3/(2pi)+(-pi^2/4+pi^4/24)r^4
     -r^5-pi^3 r^6/16.
```

All seven degree-six Bernstein coefficients of `p` on `[0,1]` are
rigorously positive; their Arb lower margins range from `0.2356...` to
`2.5075...`. Therefore

```text
Phi(u)-Psi_dB(u)>0
```

for every finite `u>=0` as well.

Now consider all `beta,gamma>=0`. The exact tail expansion begins with

```text
exp(pi q)(1-aq)=1+bq+((pi^2-3)/2)q^2+O(q^3).
```

Thus eventual residual nonnegativity requires `beta<=b`. Using `pi<22/7`
gives `b<821/308<27/10`. At the origin, global residual nonnegativity also
requires

```text
4pi^2 exp(-2pi)(1+beta+gamma)<=Phi(0).
```

Arb summation of the first eight exact theta terms, followed by the geometric
tail majorant

```text
2pi^2*9^4*exp(-81pi)/(1-(10/9)^4 exp(-19pi)),
```

certifies

```text
Phi(0)/(4pi^2 exp(-2pi))<61/10.
```

Consequently every positive three-block candidate with a globally
nonnegative residual lies in the rational triangle

```text
0<=beta<27/10, gamma>=0, beta+gamma<51/10.
```

With `E=H_0-P_(9/4)`, `P5=P_(5/4)`, and `P1=P_(1/4)`, the exact bivariate
Laguerre expression is

```text
L[R_(beta,gamma)]
 =L[E]-beta B[E,P5]-gamma B[E,P1]+beta^2 L[P5]
  +beta gamma B[P5,P1]+gamma^2 L[P1].
```

Acb Cauchy jet enclosures at `x=48,52,86`, using central-difference step
`10^-6`, Cauchy radius `0.04`, and Acb box radius `0.1`, certify this
quadratic on a closed `1/80` parameter mesh. The boxes are

```text
I_i=[i/80,(i+1)/80], J_j=[j/80,(j+1)/80],
0<=i<216, 0<=j<408-i.
```

All `64,908` boxes are strictly negative at at least one of the three
spectral points. With outward-rounded stored coefficient balls, `12,929`
boxes are assigned to `x=48`, `5,833` to `x=52`, and `46,146` to `x=86`.
In particular,

```text
L[H_0-B_P2](86)=[-2.34e-27 +/- 5.35e-30],
L[H_0-B_dB](52)=[-4.8e-17 +/- 1.15e-19].
```

Hence both classical real-zero bases leave positive residual kernels whose
transforms fail a necessary Laguerre-Polya inequality. More generally, no
nonnegative-coefficient 9/5/1 block with a globally nonnegative residual can
treat that residual as an independent Laguerre-Polya component.

Gasper's published square proves the single-shift theorem for
`F_(a,c)=K_(i(z-ic))(a)+K_(i(z+ic))(a)`: at a zero it obtains equality of two
Bessel moduli and converts their difference into one positive square
integral. A multi-shift linear combination has mixed products between
distinct shifts, and that single-shift identity does not sign them. This is a
scope guard, not a proof that no new matrix-valued or coupled square exists.

The remaining Gasper handoff is therefore a signed higher-block or genuinely
coupled mixed-term identity after full modular cancellation. This is not
strict Laguerre positivity, RH, or `Lambda<=0`.

### Lemma 11.16: Signed Universal-Factor Residual Obstruction

The standard Polya universal-factor hypothesis can be audited on the full
signed 9/5/1 coefficient plane. Define

```text
U_(beta,gamma)(z)=cosh(9z)+beta cosh(5z)+gamma cosh(z).
```

The Chebyshev identities give

```text
U_(beta,gamma)(z)=cosh(z) Q_(beta,gamma)(cosh(z)^2),

Q_(beta,gamma)(y)
 =256y^4-576y^3+(432+16beta)y^2-(120+20beta)y
  +9+5beta+gamma.
```

The zeros of `U_(beta,gamma)` are all imaginary if and only if all four
zeros of `Q_(beta,gamma)` are real and lie in `[0,1]`. Indeed,
`cosh(z)^2=r` has only imaginary solutions exactly for `0<=r<=1`.

This root condition and global residual nonnegativity imply a small rational
outer rectangle. First,

```text
Q'(1)=40+12beta>=0,
```

so `beta>=-10/3>-17/5`. Rolle's theorem also forces `Q'` to have three real
zeros. Its discriminant is

```text
disc_y Q'=-2^22 f(beta),
f(beta)=32beta^3-297beta^2+1053beta-1215.
```

The derivative `f'` has discriminant `-51516`, hence `f` is strictly
increasing, while `f(11/5)=607/125>0`. Therefore `beta<11/5`.

Put `s=beta+gamma`. If the roots of `Q` are `r_j` in `[0,1]`, then

```text
Q(1)=1+s=256 product_j(1-r_j)>=0,
```

so `s>=-1`. The origin residual inequality from Lemma 11.15 gives
`s<51/10`. Every admissible signed candidate is consequently contained in

```text
-17/5<=beta<11/5, -1<=s=beta+gamma<51/10.
```

The second algebraic exclusion uses

```text
disc_y Q=2^24 Delta(beta,gamma),

Delta=-20beta^5+16beta^4 gamma+64beta^4+124beta^3 gamma
      -184beta^3-128beta^2 gamma^2-384beta^2 gamma+432beta^2
      -48beta gamma^2+1080beta gamma-972beta
      +256gamma^3-459gamma^2-486gamma+729.
```

Thus `f(beta)>0` excludes three real critical points, and `Delta<0`
excludes four real quartic roots.

On the rational rectangle, use the exact residual quadratic from Lemma 11.15
and Acb Cauchy jets at `x=86` and `x=122`. Start with `56*61=3,416` boxes of
side `1/10` in `(beta,s)`. A box is accepted, in order, when

```text
L[R_(beta,gamma)](86)<0,
L[R_(beta,gamma)](122)<0,
f(beta)>0,
Delta(beta,gamma)<0.
```

Otherwise bisect it in both coordinates. Outward-rounded Arb arithmetic at
170 decimal digits terminates after `4,094` leaves with maximum depth six:

```text
x=86 Laguerre leaves:       2329,
x=122 Laguerre leaves:      1281,
critical-discriminant leaves: 240,
quartic-discriminant leaves:  244,
unresolved leaves:             0.
```

The closest accepted bounds remain strict: the Laguerre upper bounds are at
most about `-7.87e-32` and `-1.56e-38`, the critical-factor lower bound is
about `+0.419`, and the quartic-discriminant upper bound is about `-0.198`.

It follows that every point in the necessary rectangle either violates the
first residual Laguerre inequality or cannot satisfy the imaginary-zero
multiplier hypothesis. Therefore no signed 9/5/1 base certified by the
standard Polya universal-factor theorem can be paired with a globally
nonnegative residual treated as an independent Laguerre-Polya component.

This closes the full signed three-shift universal-factor decomposition, not
all signed comparison arguments. Higher shifts and a new matrix-valued square
that controls mixed terms together remain open. This is not strict Laguerre
positivity, RH, or `Lambda<=0`.

### Lemma 11.17: Lambda=-100 Scaled-Curvature Growth

Let `M_t^(1)` denote the continuous first-summand moment family at
`lambda=-100`, and put

```text
F(t)=log M_t^(1),
H(t)=log Gamma(t+1/2)-F(t),
B(t)=H(t+1)-2H(t)+H(t-1),
C(t)=(2t+1)B(t).
```

The centered second difference has the exact tent representation

```text
B(t)=integral_[-1,1] (1-|s|) H''(t+s) ds.
```

Writing

```text
Q(r)=2H''(r)+(2r+1)H'''(r),
```

differentiation gives

```text
C'(t)=integral_[-1,1] (1-|s|)
       *(Q(t+s)-2sH'''(t+s)) ds.
```

It is consequently sufficient to prove the buffered pointwise inequality

```text
Q(r)-2|H'''(r)| > 64/(r-3)^5, r>=318.             (11.17.1)
```

For the compact log-variable mode range `0.9264<=u<=5`, the standardized
first-summand density is integrated by paired interval composite Simpson
quadrature. The fourth-derivative Simpson remainder is evaluated directly
from the potential jet through order eight. The low range uses the established
`|y|<=6`, `V^(8)/a^4<=1/50000` envelope. On `2<=u<=5`, `3,000` Arb intervals
prove the sharper envelope

```text
sup_|y|<=8 V^(8)(x_t+y/sqrt(a_t))/a_t^4 <= 1/10^10.
```

All `16,074` compact mode blocks prove (11.17.1). The weakest outward-rounded
full-kernel transfer margin occurs on `4.9998<=u<=5` and is

```text
2.7939059038581400959435161115718783e-13.
```

On the ray `u>=5`, put `q=pi exp(4u)`, `t=V'`, `a=V''`, and `b=V'''`.
After the substitutions `u=5+v` and `q=10^9+Q`, clearing positive
denominators in

```text
u^2 t ((2t+1)b/a^3-2/a)-1/5
```

gives a polynomial with `66` strictly positive coefficients. Hence its
leading value is at least `1/5`. The already certified paired-ray moment
bounds

```text
|m1+alpha/2|<=13/q,
|m2-1|<=36/q,
|kappa_3(Y)+alpha|<=120/q
```

together with `t<=3uq`, `a>=(39/10)u^2q`, the elementary Hurwitz-zeta
bounds for the Gamma derivatives, and the `H'''` buffer consume at most

```text
2106009471479000117/29250000000000000000
```

of that normalized margin. The remaining ray margin is

```text
3743990528520999883/29250000000000000000 > 1/10,
```

which proves (11.17.1) on `u>=5` and also dominates `64/(r-3)^5`.

For the complete kernel write

```text
M_k=M_k^(1)(1+delta_k), 0<=delta_k<=2/k^6.
```

If `e_k=log(1+delta_k)`, direct expansion of the two adjacent scaled
curvatures gives

```text
|(C_(k+1)-C_k)-(C1_(k+1)-C1_k)|
 <=2[(2k+3)/(k+2)^6 +(6k+7)/(k+1)^6
      +(6k+5)/k^6 +(2k+1)/(k-1)^6]
 <=64/(k-1)^5,
```

for `k>=319`. Thus (11.17.1) proves strict complete-kernel growth on that
tail. Independently, the precedence-merged repaired Arb source proves

```text
C_(k+1)-C_k>0, 1<=k<=318,
```

with minimum lower margin
`7.8681970512268396360017168088540193e-4` at `k=318`. Therefore the actual
zeta heat-flow coefficients satisfy

```text
C_(k+1)>=C_k for every integer k>=1 at lambda=-100.       (11.17.2)
```

This is an interval and analytic theorem at one heat parameter. It supplies
neither PF-infinity, any higher-degree Jensen cone, RH, nor `Lambda<=0`.

### Lemma 11.18: Lambda=-100 Raw-Ratio Corridor

For the complete moments define

```text
R_k=M_(k+1)M_(k-1)/M_k^2,
x_k=((2k-1)/(2k+1))R_k,
B_k=-log x_k,
C_k=(2k+1)B_k.
```

The full ratio-cone entry theorem gives, at `lambda=-100`,

```text
B_k>=0, B_(k+1)<=B_k, k>=1.                              (11.18.1)
```

Lemma 11.17 gives

```text
B_(k+1)>=((2k+1)/(2k+3))B_k.                             (11.18.2)
```

For `B>=0`, set `alpha=(2k+1)/(2k+3)`. The function

```text
alpha B-log((2k+3)/(2+(2k+1)exp(-B)))
```

vanishes at zero and has nonnegative derivative, since the derivative of the
logarithmic term is `y/(2+y)` with `0<=y<=(2k+1)`. Hence

```text
log((2k+3)/(2+(2k+1)exp(-B))) <= alpha B.                 (11.18.3)
```

Combining (11.18.1)-(11.18.3) yields the full coefficient-curvature corridor.
Returning to `R_k` gives, for every `k>=1`,

```text
((2k-1)(2k+3)/(2k+1)^2)R_k
 <=R_(k+1)
 <=(2+(2k-1)R_k)/(2k+1),
```

at `lambda=-100`. This closes the zeta-specific raw-corridor target at that
parameter. It supplies neither an all-order Toeplitz/Jensen bridge,
PF-infinity, RH, nor `Lambda<=0`.

### Lemma 11.19: Lambda=-100 Adaptive-Defect Closure

At `lambda=-100`, define

```text
x_k=((2k-1)/(2k+1))R_k,
d_k=1-x_k,
s_k=((2k+1)/2)d_k.
```

The full ratio-cone entry theorem gives

```text
(2k-1)/(2k+1)<=x_k<=1,
x_(k+1)>=x_k,
```

for every `k>=1`. These inequalities are exactly

```text
0<=d_k<=2/(2k+1),
d_(k+1)<=d_k,
0<=s_k<=1.                                           (11.19.1)
```

Moreover,

```text
s_(k+1)-s_k
 =(2+(2k-1)R_k-(2k+1)R_(k+1))/2.
```

The upper raw-corridor wall in Lemma 11.18 therefore gives

```text
s_(k+1)>=s_k for every k>=1.                         (11.19.2)
```

Thus the defect-tail statement and the exact adaptive scaled-defect statement
needed by the one-entry-parameter route hold from `k=1` at `lambda=-100`.
The older targets asked for the stronger simultaneous theorem at
`lambda=-25,-50,-100`; that simultaneous statement is not proved and is not
needed once full entry is established at `lambda=-100`. This lemma supplies
neither a higher-degree minor cone, PF-infinity, RH, nor `Lambda<=0`.

### Lemma 11.20: Multiplier Hausdorff Uniqueness

Suppose a positive multiset `{alpha_j}` satisfies

```text
sum_j alpha_j^(-2)<infinity
```

and the normalized contraction curvatures have the elementary multiplier
form

```text
y_k=-log x_k=sum_j -log(1-1/(k+alpha_j)^2), k>=1.       (11.20.1)
```

For `z>1`, Frullani's identity gives

```text
-log(1-z^(-2))
 =integral_0^infinity exp(-z*t)q(t)dt,
q(t)=2(cosh(t)-1)/t>0.                                 (11.20.2)
```

Indeed, the numerator after multiplication by `exp(-z*t)` is

```text
exp(-(z-1)t)+exp(-(z+1)t)-2exp(-zt),
```

and the two Frullani integrals sum to
`log(z/(z-1))+log(z/(z+1))`.

Put

```text
S(t)=sum_j exp(-alpha_j*t).
```

The summability hypothesis makes the counting measure locally finite and
`S(t)` finite for every `t>0`. Tonelli's theorem and (11.20.2) give

```text
y_k=integral_0^infinity exp(-k*t)q(t)S(t)dt.
```

After `r=exp(-t)`, this becomes the Hausdorff moment representation

```text
y_k=integral_[0,1] r^(k-1)dnu(r),
dnu(r)=q(-log r)S(-log r)dr.                            (11.20.3)
```

The measure is finite because its total mass is `y_1`. Since polynomials are
dense in `C([0,1])`, the moments `y_1,y_2,...` determine `nu` uniquely. In
particular,

```text
(-1)^m Delta^m y_k
 =integral_[0,1] r^(k-1)(1-r)^m dnu(r)>=0.              (11.20.4)
```

Conversely, the integer multiplier target holds exactly when this unique
Hausdorff measure is absolutely continuous and its density has the form

```text
dnu/dr=q(-log r)sum_j r^alpha_j                         (11.20.5)
```

for a unit-integer-multiplicity counting measure satisfying the stated
summability condition. If (11.20.5) exists, division by `q(t)>0` and
uniqueness of Laplace transforms determine the multiset `{alpha_j}` uniquely.

This does not identify the product interpolation with the natural Mellin
interpolation. The exact function

```text
f(s)=sin(2*pi*s)
```

vanishes at every integer and satisfies

```text
f(s+1)-2f(s)+f(s-1)=0.
```

Thus integer agreement, even together with agreement of centered continuous
log curvature, permits a nontrivial periodic interpolation factor. A growth
or canonical-interpolation theorem is still required before the continuous
Mellin power-sum obstruction can be applied to the integer-only product.
Lemma 11.20 proves uniqueness and an exact recovery characterization, not
existence of the unit counting measure, PF-infinity, RH, or `Lambda<=0`.

### Lemma 11.21: Conditional Leading-Atom Bounds

Assume the unit counting representation in Lemma 11.20 exists, and define

```text
a_m=(-1)^m Delta^m y_1.
```

For one atom put

```text
g_m(alpha)
 =integral_0^1 r^(alpha-1)(1-r)^(m+2)/(-log r) dr.
```

Then Tonelli's theorem gives

```text
a_m=sum_j g_m(alpha_j).                                (11.21.1)
```

Differentiation under the integral is exact:

```text
partial_alpha g_m(alpha)
 =-integral_0^1 r^(alpha-1)(1-r)^(m+2)dr
 =-Beta(alpha,m+3)<0.                                  (11.21.2)
```

Let `beta_m` be the unique solution of `g_m(beta_m)=a_m`. If
`alpha_min=min_j alpha_j`, then (11.21.1)-(11.21.2) imply

```text
alpha_min>=beta_m.                                     (11.21.3)
```

Using the 220-digit lambda-zero coefficient enclosures through `A_57`, with
all derived Arb operations at 250 digits, order six gives

```text
g_6(4.863538496)-a_6>0,
g_6(4.863538497)-a_6<0.
```

The respective outward-rounded margins are about `3.6340e-14` and
`-1.4440e-13`, with radii below `6e-130`. Hence

```text
4.863538496<beta_6<4.863538497,
alpha_min>4.863538496.                                 (11.21.4)
```

At the lower endpoint, all other orders `0<=m<=55`, `m!=6`, satisfy
`g_m(4.863538496)-a_m<0`; thus order six is the strongest bound in the
available finite triangle.

More generally, if `N(A)` counts target atoms at or below `A`, positivity and
monotonicity give

```text
N(A)g_m(A)<=a_m.                                       (11.21.5)
```

At `A=11/2`, order three yields the certified ratio

```text
a_3/g_3(11/2)<1.632<2,
```

so

```text
N(11/2)<=1.                                            (11.21.6)
```

Equations (11.21.4) and (11.21.6) are conditional constraints on any
multiplier counting measure. They supply neither an atom below `11/2` nor a
construction of the measure, PF-infinity, RH, or `Lambda<=0`.

### Lemma 11.22: Unit-Atomic Multiplier Obstruction

Retain the atom kernel and first difference moments from Lemma 11.21. At
lambda zero, Arb proves

```text
g_6(4.863538496)-a_6>0.                                (11.22.1)
```

If the proposed unit-atomic product existed, any atom with
`alpha<=4.863538496` would contribute at least the left kernel in (11.22.1),
already exceeding the complete moment `a_6`. Hence every target atom would
have to satisfy

```text
alpha_j>4.863538496.                                   (11.22.2)
```

Now put

```text
R_m(alpha)=g_(m+1)(alpha)/g_m(alpha).
```

Under the probability measure proportional to the integrand defining
`g_m(alpha)`, one has

```text
R_m(alpha)=E_alpha[1-r],
R_m'(alpha)=Cov_alpha(1-r,log r)<0.                    (11.22.3)
```

The strict sign follows from the double-integral covariance identity, because
`1-r` is strictly decreasing and `log r` is strictly increasing on `(0,1)`.
Equations (11.21.1), (11.22.2), and (11.22.3) would therefore force the
weighted-average bound

```text
a_7/a_6<R_6(4.863538496).                              (11.22.4)
```

The same 250-digit Arb source instead proves

```text
a_7/a_6-R_6(4.863538496)
 =[0.00076825009751290047 +/- 8.00e-126] >0.           (11.22.5)
```

Thus (11.22.4) and (11.22.5) contradict each other. The normalized
lambda-zero zeta coefficient sequence has no convergent elementary
multiplier product with positive atoms of unit integer multiplicity and
`sum_j alpha_j^(-2)<infinity`.

The unit-multiplicity hypothesis is essential in the cutoff step. Lemma 11.22
rejects this sufficient multiplier subclass, not every multiplier sequence,
the other all-order Jensen/PF bridges, PF-infinity, RH, or `Lambda<=0`.

### Lemma 11.23: Theta/Bessel Higher-Shift Regularization Gate

For `n>=1`, set

```text
phi_n(u)
 =(2*pi^2*n^4*exp(9u)-3*pi*n^2*exp(5u))
  *exp(-pi*n^2*exp(4u)),
g_n(u)=(phi_n(u)+phi_n(-u))/2.
```

The theta functional equation makes the complete kernel even, so the two
absolutely convergent pointwise theta series give

```text
Phi(u)=sum_(n>=1)g_n(u).                              (11.23.1)
```

For each fixed `n`, symmetrize before expanding the two remaining
exponentials. Pairing the `exp(9u)` coefficient at index `m` with the
`exp(5u)` coefficient at index `m-1` gives

```text
g_n(u)=exp(-2*pi*n^2*cosh(4u))
       *sum_(m>=0)c_(n,m)cosh((9-4m)u),              (11.23.2)

c_(n,m)=pi^m*n^(2m)*(2*pi^2*n^4-3m)/m!.             (11.23.3)
```

The separate positive exponential series furnish an integrable absolute
majorant on `[0,infinity)` for every fixed `n`. Thus (11.23.2), and its
termwise cosine transform, are ordinary absolutely justified identities at
fixed `n`. With

```text
P_(n,m)(x)
 =1/8*(K_(9/4-m+ix/4)(2*pi*n^2)
       +K_(9/4-m-ix/4)(2*pi*n^2)),
```

the standard Bessel integral gives

```text
I_n(x):=integral_0^infinity g_n(u)cos(xu)du
       =sum_(m>=0)c_(n,m)P_(n,m)(x).                 (11.23.4)
```

The coefficient sign is exact:

```text
sign(c_(n,m))=sign(2*pi^2*n^4-3m).                  (11.23.5)
```

In particular, `3<pi<22/7` puts `2*pi^2/3` strictly between `6` and
`7`; hence `c_(1,m)>0` for `0<=m<=6` and `c_(1,m)<0` for every
`m>=7`. The higher-shift expansion is intrinsically signed.

The fixed-`n` transforms cannot be summed over `n` in the ordinary sense.
Indeed, the translate identity

```text
phi_n(u)=n^(-1/2)phi_1(u+(log n)/2)
```

and the bilateral profile transform

```text
Qhat(x)=-(1+x^2)/32*pi^(-(1+ix)/4)*Gamma((1+ix)/4)
```

give at zero frequency

```text
I_n(0)=Qhat(0)/(2*sqrt(n)),
Qhat(0)=-Gamma(1/4)/(32*pi^(1/4))<0.                 (11.23.6)
```

Consequently,

```text
sum_(n>=1)I_n(0)=-infinity,                          (11.23.7)
```

whereas `H_0(0)=integral_0^infinity Phi(u)du` is finite. Thus (11.23.1)
is not interchangeable with the half-line Fourier integral. Analytic
continuation of the formal arithmetic sum only reconstructs

```text
Qhat(x)zeta((1+ix)/2)=xi((1+ix)/2)/4,
```

and supplies no positive block summation theorem.

For any finite or legitimately regularized family, the remaining mixed-term
problem has the exact matrix form

```text
L[sum_j a_jF_j]=sum_(i,j)a_i*a_j*M_(i,j),
M_(i,j)=F_i'*F_j'-(F_i*F_j''+F_j*F_i'')/2.          (11.23.8)
```

Therefore a viable higher-shift proof must perform the theta modular
cancellation before transformation and prove positivity of the coupled
quadratic form, or establish a noncircular sign-preserving renormalized
summation theorem. Lemma 11.23 closes only naive ordinary termwise
higher-shift summation; it proves neither strict Laguerre positivity, Wiener
density, RH, nor `Lambda<=0`.

### Lemma 11.24: Theta Cell Renormalization

Extend the theta index continuously. For `a>0`, define

```text
phi_a(u)
 =(2*pi^2*a^4*exp(9u)-3*pi*a^2*exp(5u))
  *exp(-pi*a^2*exp(4u)),
g_a(u)=(phi_a(u)+phi_a(-u))/2.
```

The translate law remains exact:

```text
phi_a(u)=a^(-1/2)phi_1(u+(log a)/2).                 (11.24.1)
```

Put `q=pi*exp(4u)`. The Gaussian moments

```text
integral_0^infinity a^2*exp(-q*a^2)da
 =sqrt(pi)/(4*q^(3/2)),
integral_0^infinity a^4*exp(-q*a^2)da
 =3*sqrt(pi)/(8*q^(5/2))
```

cancel exactly in the defining linear combination. Therefore, for every real
`u`,

```text
integral_0^infinity phi_a(u)da
 =integral_0^infinity g_a(u)da=0.                   (11.24.2)
```

Define the continuum-cell residual

```text
r_n(u)=g_n(u)-integral_(n-1)^n g_a(u)da.             (11.24.3)
```

The cell intervals partition `(0,infinity)`, so Lemma 11.23 and (11.24.2)
give the exact pointwise identity

```text
Phi(u)=sum_(n>=1)r_n(u).                             (11.24.4)
```

Let `v=u+(log a)/2`. Differentiating (11.24.1) gives

```text
partial_a phi_a(u)
 =a^(-3/2)*(phi_1'(v)-phi_1(v))/2.                  (11.24.5)
```

Translation in `u` and the rapid decay of `phi_1` imply, for every integer
`k>=0` and `a>=1`,

```text
integral_R |u|^k*|partial_a g_a(u)|du
 <=C_k*a^(-3/2)*(1+|log a|^k).                     (11.24.6)
```

For `n>=2`, write `r_n` as the integral of `g_n-g_a` across its cell and use
(11.24.6). Near `a=0`, (11.24.1) instead gives the integrable majorant
`C_k*a^(-1/2)*(1+|log a|^k)`. Hence

```text
sum_(n>=1) integral_R |u|^k*|r_n(u)|du<infinity
for every integer k>=0.                              (11.24.7)
```

Set

```text
J_n(x)=integral_0^infinity r_n(u)cos(xu)du.
```

Equation (11.24.7) permits termwise transformation and differentiation to
every fixed order, absolutely and uniformly on the real axis:

```text
H_0^(k)(x)=sum_(n>=1)J_n^(k)(x).                     (11.24.8)
```

This renormalization has an explicit arithmetic transform. For
`s=(1+ix)/2`, put

```text
e_n(s)
 =n^(-s)-(n^(1-s)-(n-1)^(1-s))/(1-s).              (11.24.9)
```

The profile transform in Lemma 11.23 gives

```text
J_n(x)
 =1/4*(Qhat(x)e_n(s)+Qhat(-x)e_n(conjugate(s))).    (11.24.10)
```

The cell terms telescope. Euler's convergent sum-integral limit gives, for
`0<Re(s)<1`,

```text
sum_(n>=1)e_n(s)
 =lim_(N->infinity)
   (sum_(n=1)^N n^(-s)-N^(1-s)/(1-s))
 =zeta(s).                                           (11.24.11)
```

Thus (11.24.8)-(11.24.11) reconstruct `H_0` through an ordinary convergent
series, without using any information about zeta zeros.

At zero frequency,

```text
J_n(0)
 =Qhat(0)/2*(n^(-1/2)-2*(sqrt(n)-sqrt(n-1))).       (11.24.12)
```

The cell integral of `a^(-1/2)` is strictly larger than its right-endpoint
value, while `Qhat(0)<0`. Consequently,

```text
J_n(0)>0 for every n>=1.                             (11.24.13)
```

Normal convergence through the second jet also gives the exact coupled
Laguerre matrix

```text
L[H_0](x)=sum_(m,n>=1)M_(m,n)(x),                    (11.24.14)

M_(m,n)
 =J_m'*J_n'-(J_m*J_n''+J_n*J_m'')/2,
```

with an absolutely and locally uniformly convergent double series. No sign of
this matrix sum is asserted.

Finally, the endpoint decomposition is not compatible block by block with
positive Newman time. Expanding the negative-side profile gives

```text
r_n(u)=-(pi/2)*(3n-1)exp(-5u)+O_n(exp(-9u))          (11.24.15)
```

because

```text
n^2-integral_(n-1)^n a^2 da=n-1/3.
```

For every `t>0`, `exp(tu^2)r_n(u)` is therefore nonintegrable. A direct
positive-time proof must cancel the full modular tail inside each admissible
group before applying the Gaussian weight, or work with the already
endpoint-subtracted theta primitive. Lemma 11.24 supplies a rigorous endpoint
renormalization and matrix identity, not the positive-time strict Laguerre
theorem, Wiener density, RH, or `Lambda<=0`.

### Lemma 11.25: Reciprocal-Gamma Mixture Sign Gate

Write the even moment sequence as

```text
m_k(lambda)=mu_(2k)(lambda)
           =integral_0^infinity t^k rho_lambda(dt),
A_k(lambda)=gamma_k*m_k(lambda),
gamma_k=k!/(2k)!.
```

The duplication formula gives

```text
gamma_k=sqrt(pi)/(4^k*Gamma(k+1/2)).                 (11.25.1)
```

Karlin's reciprocal-gamma sign-regularity theorem states that, if a minor has
size `nu+1` and `alpha>-nu`, then every minor of

```text
[1/Gamma(alpha+i+j)]_(i,j>=0)
```

with increasing row and column sets has strict sign

```text
(-1)^(nu*(nu+1)/2).                                  (11.25.2)
```

This is Lemma 9.2 of S. Karlin, *Total positivity, absorption probabilities
and applications*, Trans. Amer. Math. Soc. 111 (1964), 33-107. Applying
(11.25.2) with `alpha=n+1/2` and using positive row and column scalings proves
that, for every `t>0` and `n>=0`,

```text
K_(n,t)(i,j)=gamma_(n+i+j)*t^(n+i+j)                (11.25.3)
```

is strictly sign-regular in every order, with signature
`epsilon_k=(-1)^(k(k-1)/2)`.

For the leading contiguous `k` by `k` minor one may evaluate the determinant
explicitly:

```text
det[gamma_(i+j)]_(i,j=0..k-1)
 =(-1)^(k(k-1)/2)*pi^(k/2)*4^(-k(k-1))
  *prod_(r=1)^(k-1)r!
  /prod_(i=0)^(k-1)Gamma(k-1/2+i).                  (11.25.4)
```

Thus the reciprocal-gamma factor itself has exactly the signature seen in the
finite shifted-Hankel staircase.

Positive scale mixing is the missing step. Multilinearity in the determinant
rows gives, for `0<=j_1<...<j_k`,

```text
R_(k,n)(j_1,...,j_k)
 =integral_((0,infinity)^k)
   det[gamma_(n+i+j_l)*t_i^(n+i+j_l)]
       _(i=0..k-1,l=1..k)
   prod_(i=0)^(k-1)rho_lambda(dt_i).                (11.25.5)
```

The scale variables in (11.25.5) are independent. Karlin's theorem signs the
common-scale diagonal `t_0=...=t_(k-1)`, not the full integrand. Already for
`k=2`, `n=0`, and columns `(0,1)`, symmetrization in the two row scales gives

```text
D_sym(t_0,t_1)
 =(t_0^2-6*t_0*t_1+t_1^2)/24.                      (11.25.6)
```

It satisfies

```text
D_sym(1,1)=-1/6,
D_sym(1,10)=41/24.
```

Hence the natural determinant integrand has no fixed sign. The failure is not
merely pointwise. For the positive measure

```text
rho=10*delta_(1/100)+delta_1+delta_2,
```

one obtains

```text
(m_0,m_1,m_2)=(12,31/10,5001/1000),
(A_0,A_1,A_2)=(12,31/20,1667/4000),
A_0*A_2-A_1^2=5197/2000>0.                           (11.25.7)
```

The required order-two signed-Hankel sign is negative. Thus positive mixing
of fixed-scale strictly sign-regular reciprocal-gamma kernels is not closed,
even when every scale is strictly positive.

For the actual Xi measure, the completed order-two cone has a useful exact
interpretation. Define the `n`-tilted probability law

```text
P_n(dt)=t^n*rho_lambda(dt)/m_n(lambda)
```

and its squared coefficient of variation

```text
CV_n^2
 =Var_(P_n)(T)/E_(P_n)(T)^2
 =m_n*m_(n+2)/m_(n+1)^2-1.
```

Then

```text
A_(n+1)^2>=A_n*A_(n+2)
iff m_(n+1)^2/(m_n*m_(n+2))>=(2n+1)/(2n+3)
iff CV_n^2<=2/(2n+1).                                (11.25.8)
```

The full ratio-cone entry at `lambda=-100` and infinite heat-flow invariance
therefore prove (11.25.8) for every shift at `lambda=0`. This closes the
complete order-two signed-Hankel layer by an Xi-specific concentration
theorem.

Higher orders remain open. They require a hierarchy controlling the
sign-changing compound integrand in (11.25.5), a coupling that keeps the row
scales sufficiently close to the common-scale diagonal, or a different
positive compound-kernel factorization. Lemma 11.25 proves the fixed-scale
all-order theorem, the exact mixture obstruction, and the Xi order-two
concentration result; it does not prove higher-order mixture sign regularity,
the signed-Hankel/Jensen bridge, PF-infinity, RH, or `Lambda<=0`.

### Lemma 11.26: Reciprocal-Defect Compound Order Three

Let

```text
rho_k=A_(k+1)/A_k,
x_k=rho_k/rho_(k-1),
d_k=1-x_k,
q_k=d_k^(-1/2).
```

For the contiguous order-three signed-Hankel minor

```text
D_(3,n)=det[A_(n+i+j)]_(i,j=0..2),
```

direct elimination of the coefficient ratios gives

```text
D_(3,n)=A_n^3*rho_n^6*x_(n+1)^3*F_n,                (11.26.1)

F_n
 =x_(n+1)*x_(n+2)^2*x_(n+3)
  -x_(n+1)*x_(n+2)^2
  -x_(n+2)^2*x_(n+3)+2*x_(n+2)-1.
```

All factors outside `F_n` are positive. Substitution of `x_k=1-d_k`
collapses the frontier to

```text
F_n
 =d_(n+1)*d_(n+3)*x_(n+2)^2-d_(n+2)^2.             (11.26.2)
```

The required order-three signature is negative. Hence

```text
D_(3,n)<0
iff d_(n+2)^2>x_(n+2)^2*d_(n+1)*d_(n+3).           (11.26.3)
```

In reciprocal-defect coordinates define

```text
C_n=q_(n+1)*q_(n+3)-q_(n+2)^2+1.                   (11.26.4)
```

The exact positive-factor identity is

```text
d_(n+2)^2-x_(n+2)^2*d_(n+1)*d_(n+3)
 =C_n*(q_(n+1)*q_(n+3)+q_(n+2)^2-1)
  /(q_(n+1)^2*q_(n+2)^4*q_(n+3)^2).                (11.26.5)
```

Inside the full ratio cone the second factor in (11.26.5) is positive, so

```text
D_(3,n)<0 iff C_n>0.                                (11.26.6)
```

This is unit-buffered log-concavity of the reciprocal defects. If

```text
a_n=q_(n+2)-q_(n+1),
b_n=q_(n+3)-q_(n+2),
```

then

```text
C_n=1-a_n*b_n+q_(n+2)*(b_n-a_n).                   (11.26.7)
```

Consequently `0<=a_n<=b_n<=1` is sufficient for `C_n>=0`, with strict
positivity unless `a_n=b_n=1`. On the reciprocal-square boundary family
`q_k=alpha+beta*k`, `0<=beta<1`, one has the strict constant margin

```text
C_n=1-beta^2>0.                                     (11.26.8)
```

This sufficient increment theorem is not implied by the previously proved
cones. Consider

```text
q_1=10, q_2=109/10, q_3=58/5.
```

Then

```text
(d_1,d_2,d_3)=(1/100,100/11881,25/3364),
(x_1,x_2,x_3)=(99/100,11781/11881,3339/3364),
(s_1,s_2,s_3)=(3/200,250/11881,175/6728).
```

The contractions lie strictly inside their full ratio-cone walls and
increase; the defects decrease; the scaled defects increase and remain below
one; and the reciprocal-defect increments are `9/10` and `7/10`, both
strictly below one. The two neighboring cubic Jensen frontier values are

```text
-828039/1411581610000,
-1610484375/1597415764323856,
```

so both cubic windows are strictly hyperbolic. Nevertheless

```text
C_0=-181/100<0.
```

The positive coefficient prefix generated by these contractions satisfies

```text
det[A_(i+j)]_(i,j=0..2)>0
```

with exact value

```text
4106267526339/1899424214416000000.
```

Thus the full ratio cone, decreasing defect, increasing scaled defect, and
strict cubic Jensen cone do not imply the contiguous order-three
signed-Hankel sign. The first higher compound concentration target is the
all-shift theorem `C_n>0`, initially at `lambda=-100`, together with a
forward-invariance argument to `lambda=0`. Lemma 11.26 supplies the exact
coordinate and a strict promotion countermodel; it does not prove this new
curvature theorem, arbitrary order-three column sets, the all-order
signed-Hankel antecedent, the Jensen bridge, RH, or `Lambda<=0`.

### Lemma 11.27: Contiguous Order-Three Entry at Lambda=-100

At `lambda=-100`, retain the coordinates of Lemma 11.26 and put

```text
s_k=((2*k+1)/2)*d_k.
```

The repaired coefficient enclosures through `A_321`, evaluated with
1024-bit Arb arithmetic, prove

```text
C_n>0,                                                     (11.27.1)
d_(n+2)^2-x_(n+2)^2*d_(n+1)*d_(n+3)>0,
0<=n<=317,
```

and the splice anchor

```text
s_319>251/500.                                             (11.27.2)
```

The all-index adaptive-defect theorem at this heat parameter gives

```text
d_(j+1)<=d_j,
s_(j+1)>=s_j,
j>=1.                                                      (11.27.3)
```

It remains to control the tail without extrapolating the finite data. For a
center `k>=320`, set

```text
m=2*k+1>=641,
u_k=q_k-q_(k-1),
v_k=q_(k+1)-q_k.
```

By (11.27.3), `u_k,v_k>=0`. The scaled-defect part of (11.27.3) and
(11.27.2) imply, for every `j>=319`,

```text
d_j>251/(250*(2*j+1)),
q_j^2<250*(2*j+1)/251.                                    (11.27.4)
```

Moreover,

```text
(2*k-1)*d_(k-1)<=(2*k+1)*d_k,
```

and therefore

```text
u_k<=q_k*r_m,
r_m=1-sqrt(1-2/m).                                        (11.27.5)
```

For every `m>=641`,

```text
r_m<1/m+1/m^2.                                            (11.27.6)
```

Indeed the right side in the equivalent square-root inequality is positive,
and after squaring the difference is

```text
(1-2/m)-(1-1/m-1/m^2)^2
 =(m^2-2*m-1)/m^4>0.
```

Combining (11.27.4)--(11.27.6), and applying the same estimate at `m+2`,
gives

```text
q_k*u_k<(250/251)*(1+1/m),                                (11.27.7)
u_k*v_k<(250/251)*(1/m)*(1+1/m)^2.                        (11.27.8)
```

For (11.27.8), both `u_k^2` and `v_k^2` are bounded by the displayed
right side; the elementary function `(1/m)*(1+1/m)^2` decreases for
positive `m`.

The exact increment form from Lemma 11.26 now yields

```text
C_(k-2)
 =1-u_k*v_k+q_k*(v_k-u_k)
 >=1-u_k*v_k-q_k*u_k
 >1-(250/251)*(1+2/m+2/m^2+1/m^3).                       (11.27.9)
```

The bracket in (11.27.9) decreases with `m`. At `m=641`, exact rational
arithmetic gives

```text
C_(k-2)>57613471/66107054971>0.                           (11.27.10)
```

Thus (11.27.10) covers `n=k-2>=318`, while (11.27.1) covers
`0<=n<=317`. Consequently

```text
C_n(-100)>0,
D_(3,n)(-100)<0,
for every integer n>=0.                                  (11.27.11)
```

This is an all-shift entry theorem for the contiguous order-three compound
sign at one heat parameter. It does not prove forward invariance to
`lambda=0`, noncontiguous order-three minors, any higher compound order, the
all-order signed-Hankel/Jensen bridge, PF-infinity, RH, or `Lambda<=0`.

### Lemma 11.28: Contiguous Order-Three Forward Invariance

Along the actual positive heat-flow coefficient trajectory, the reciprocal
defects satisfy the exact one-sided equation

```text
q_j'
 =r_j*((2*j-1)*q_j
       -(2*j+3)*(q_j^3-q_j)/q_(j+1)^2).             (11.28.1)
```

For `k=n+2`, abbreviate

```text
a=q_(k-1), b=q_k, c=q_(k+1), d=q_(k+2),
C_n=a*c-b^2+1,
C_(n+1)=b*d-c^2+1.
```

Using `r_(k-1)=r_k/(1-b^(-2))` and
`r_(k+1)=r_k*(1-c^(-2))`, direct differentiation of `C_n` in (11.28.1)
gives the cooperative identity

```text
C_n'/r_k=alpha_k*C_(n+1)+beta_k*C_n,                (11.28.2)

alpha_k
 =(2*k+5)*a*(b*d+c^2-1)/(c*d^2),                   (11.28.3)

beta_k=-N_k/(c^2*(b^2-1)),                          (11.28.4)

N_k
 =(2*k+1)*a^2*c^2+(2*k+1)*a*b^2*c-(2*k+1)*a*c
  +(4*k+6)*b^4+(-4*k+2)*b^2*c^2-(4*k+6)*b^2.
```

The strict moment upper wall gives `q_j>1`. Hence every denominator in
(11.28.2)--(11.28.4) is nonzero and

```text
alpha_k>0.                                            (11.28.5)
```

In particular, at `C_n=0`,

```text
C_n'/r_k=alpha_k*C_(n+1).                            (11.28.6)
```

Thus the finite-dimensional boundary algebra points inward whenever the next
compound margin is nonnegative. It remains to justify the infinite system.

Fix a finite `L>=-100`. The completed ratio cone gives

```text
q_j^2>=(2*j+1)/2,
r_j<=r_0.
```

The cubic forward-uniform theorem supplies a constant `B_L` such that, on
`-100<=lambda<=L`,

```text
0<=g_j:=q_(j+1)-q_j<=B_L/sqrt(j)                    (11.28.7)
```

on the spatial tail. The finitely many omitted indices are harmless by
continuity. Summing (11.28.7) and using compactness of one fixed `q_K` gives

```text
sqrt(j)<=q_j<=Q_L*sqrt(j)                            (11.28.8)
```

after increasing `Q_L` if necessary.

The coefficient cancellation in (11.28.2) is essential. Put

```text
k=h^(-2),
b=y/h,
a=(y-u*h^2)/h,
c=(y+v*h^2)/h,
d=(y+(v+w)*h^2)/h,                                  (11.28.9)
```

where `u,v,w` are the neighboring increments multiplied by `sqrt(k)`.
By (11.28.7)--(11.28.8), `(y,u,v,w)` stays in a compact box with `y>=1`.
Substitution of (11.28.9) into (11.28.3)--(11.28.4) and exact cancellation
give

```text
lim_(h->0) h^2*alpha_k=4,                            (11.28.10)

lim_(h->0)(alpha_k+beta_k)
 =2*(u+2*v-3*w)/y.                                  (11.28.11)
```

After cancellation, the only denominator factors are, up to sign,

```text
(y^2-h^2)*(y+v*h^2)^2*(y+(v+w)*h^2)^2.
```

They stay uniformly away from zero on the compact tail box. Equations
(11.28.10)--(11.28.11), together with continuity at the finitely many initial
indices, prove

```text
alpha_k=O_L(k),
alpha_k+beta_k=O_L(1).                               (11.28.12)
```

Also, from

```text
C_n=1-g_(k-1)*g_k+q_k*(g_k-g_(k-1)),
```

equations (11.28.7)--(11.28.8) imply `C_n=O_L(1)` uniformly in `n`.
Define

```text
z_n=C_n/(n+1).                                       (11.28.13)
```

Then `z_n->0` uniformly on the compact heat interval. Dividing (11.28.2) by
`n+1=k-1` gives

```text
z_n'
 =r_k*alpha_k*(k/(k-1))*z_(n+1)+r_k*beta_k*z_n.      (11.28.14)
```

By (11.28.12) and `r_k<=r_0`, the effective diagonal coefficient

```text
r_k*(alpha_k+beta_k+alpha_k/(k-1))                  (11.28.15)
```

has a finite upper bound `M_L`. Set
`y_n(lambda)=exp(-M_L*(lambda+100))*z_n(lambda)`. If the spatial infimum of
the `y_n` is negative, uniform tail decay makes it occur at a finite index.
At an active minimum, `y_(n+1)>=y_n`; (11.28.5), (11.28.14), and
(11.28.15) therefore give the nonnegative upper-right derivative required by
the same connected-component minimum argument used for the ratio cone. A
negative component cannot form.

Lemma 11.27 supplies `C_n(-100)>0` for every `n`. Hence

```text
C_n(lambda)>=0
```

for every finite `lambda>=-100`. Strictness at each fixed index follows from
variation of constants in (11.28.2): the initial term is strictly positive,
`alpha_k>=0`, and `C_(n+1)>=0`. Consequently

```text
C_n(lambda)>0 for every n>=0 and finite lambda>=-100,
D_(3,n)(0)<0 for every n>=0.                         (11.28.16)
```

Lemma 11.28 closes the complete shifted contiguous order-three layer at
`lambda=0`. It does not prove noncontiguous order-three minors, any compound
order four or higher, the all-order signed-Hankel/Jensen bridge, PF-infinity,
RH, or `Lambda<=0`.

### Lemma 11.29: Noncontiguous Order-Three Secant Transfer

Fix `n>=0` and a heat parameter on the trajectory of Lemma 11.28. For a
column offset `j`, positive column scaling gives

```text
(A_(n+j),A_(n+j+1),A_(n+j+2))^T
 =A_(n+j)*(1,u_j,v_j)^T,                             (11.29.1)

u_j=r_(n+j),
v_j=r_(n+j)*r_(n+j+1).
```

The strict upper ratio wall gives

```text
u_(j+1)<u_j.                                         (11.29.2)
```

View `P_j=(u_j,v_j)` as a planar point and define its successive edge slope

```text
sigma_j=(v_(j+1)-v_j)/(u_(j+1)-u_j).                 (11.29.3)
```

Direct determinant algebra gives

```text
det[(1,u_(j+l),v_(j+l))^T]_(l=0,1,2)
 =(u_(j+1)-u_j)*(u_(j+2)-u_(j+1))
  *(sigma_(j+1)-sigma_j).                            (11.29.4)
```

The first two factors in (11.29.4) are both negative. Lemma 11.28 proves the
strict negative sign of every contiguous order-three determinant. Therefore

```text
sigma_(j+1)<sigma_j                                  (11.29.5)
```

for every `j>=0`.

For arbitrary integers `a<b`, let

```text
S_(a,b)=(v_b-v_a)/(u_b-u_a).
```

Telescoping numerator and denominator shows that `S_(a,b)` is the weighted
average of `sigma_a,...,sigma_(b-1)` with strictly positive weights
`u_l-u_(l+1)`. If `a<b<c`, (11.29.5) implies that every edge slope in the
first interval is larger than every edge slope in the second, and hence

```text
S_(a,b)>S_(b,c).                                     (11.29.6)
```

Another application of (11.29.4), now to the three nonconsecutive points,
gives

```text
det[(1,u_l,v_l)^T]_(l=a,b,c)
 =(u_b-u_a)*(u_c-u_b)*(S_(b,c)-S_(a,b))<0.           (11.29.7)
```

Restoring the positive factors from (11.29.1) proves

```text
R_(3,n)(j_1,j_2,j_3)<0                               (11.29.8)
```

for every `n>=0` and `0<=j_1<j_2<j_3` at `lambda=0` (indeed at every finite
`lambda>=-100`). For two columns, the same normalization gives

```text
R_(2,n)(j_1,j_2)
 =A_(n+j_1)*A_(n+j_2)
  *(r_(n+j_2)-r_(n+j_1))<0.                          (11.29.9)
```

In particular, `R_(2,n)(j_1,j_2)<0` for every `j_1<j_2`.

Thus the complete arbitrary-column reshaped-Hankel layers of orders two and
three have the required signatures. The planar secant argument has no
automatic order-four analogue. Lemma 11.29 does not prove compound order four
or higher, the all-order sign-regular-to-Jensen transfer, PF-infinity, RH, or
`Lambda<=0`.

### Lemma 11.30: Contiguous Order-Four Condensation Frontier

Write

```text
H_(m,n)=det[A_(n+i+j)]_(i,j=0..m-1).
```

Desnanot-Jacobi condensation applied to the contiguous `4` by `4` Hankel
block gives

```text
H_(4,n)*H_(2,n+2)
 =H_(3,n)*H_(3,n+2)-H_(3,n+1)^2.                    (11.30.1)
```

Lemmas 11.25 and 11.29 give `H_(2,j)<0` and `H_(3,j)<0`. Put

```text
T_j=-H_(3,j)>0.
```

Since the order-four signature is positive, (11.30.1) gives the exact
equivalence

```text
H_(4,n)>0
iff T_(n+1)^2>T_n*T_(n+2).                          (11.30.2)
```

Thus contiguous order four is strict log-concavity of the magnitudes of the
completed contiguous order-three determinants.

In the defect coordinates of Lemma 11.26, define the positive order-three gap

```text
G_n
 =d_(n+2)^2-x_(n+2)^2*d_(n+1)*d_(n+3)>0.           (11.30.3)
```

The determinant factorization gives

```text
T_n=A_n^3*r_n^6*x_(n+1)^3*G_n.                      (11.30.4)
```

Taking two consecutive ratios in (11.30.4) yields

```text
T_(n+1)^2/(T_n*T_(n+2))
 =G_(n+1)^2/(x_(n+3)^3*G_n*G_(n+2)).                (11.30.5)
```

Consequently the exact new contiguous order-four target is

```text
G_(n+1)^2>x_(n+3)^3*G_n*G_(n+2)                    (11.30.6)
```

for every `n>=0`.

The repaired `lambda=-100` coefficient enclosures through `A_322` prove
(11.30.6) for all `0<=n<=316`. All `317` Arb margins are strictly positive;
the smallest outward-rounded lower bound is greater than

```text
0.004717895386893703020409137714
```

at `n=316`. This is a finite prefix theorem, not the missing all-index tail.

The lower compound layers do not imply (11.30.6). Consider

```text
s_1=3/10, s_2=2/5, s_3=41/100, s_4=49/100, s_5=11/20.  (11.30.7)
```

Then

```text
(d_1,...,d_5)=(1/5,4/25,41/350,49/450,1/10),
(x_1,...,x_5)=(4/5,21/25,309/350,401/450,9/10),
(G_0,G_1,G_2)=(1417/156250,10943/76562500,603553/236250000).
```

The contractions satisfy every strict pointwise ratio wall and increase; the
defects decrease; and the scaled defects increase. The four exact cubic
frontiers are

```text
-319/15625,
-89697/10937500,
-15257191/2756250000,
-9791/2250000,
```

so every available cubic window is strict. The positive coefficient prefix
obtained from `A_0=r_0=1` has every available order-two, contiguous
order-three, and arbitrary-column order-three determinant with the required
sign. However,

```text
G_1^2-x_3^3*G_0*G_2
 =-933340447356927/58618164062500000000<0,          (11.30.8)

H_(4,0)
 =-6608596712914764288/582076609134674072265625<0.  (11.30.9)
```

The required order-four sign is positive. Thus even the completed ratio,
cubic, and order-three layers cannot be promoted to order four without the
new buffered gap-log-concavity theorem (11.30.6).

Lemma 11.30 supplies the exact contiguous order-four coordinate, a 317-row
interval prefix, and a strict lower-order promotion countermodel. It does not
prove the all-index `lambda=-100` tail, forward order-four invariance,
arbitrary-column order four, the all-order bridge, PF-infinity, RH, or
`Lambda<=0`.

The tail obligation admits a sharper quantitative form. Put `k=n+3` and

```text
P_n=log(G_n*G_(n+2)/G_(n+1)^2).                     (11.30.10)
```

Equation (11.30.6) is equivalent to

```text
P_n<-3*log(x_k).                                     (11.30.11)
```

At `lambda=-100`, the increasing scaled-defect theorem and the certified
anchor `s_319>251/500` give, for `k>=319`,

```text
-3*log(x_k)>3*d_k>=753/(250*(2*k+1)).                (11.30.12)
```

For the reduced tail `k>=320`, moreover,

```text
4/k^2<753/(250*(2*k+1)), k>=320.                    (11.30.13)
```

After clearing denominators, the numerator in (11.30.13) is
`753*k^2-1000*(2*k+1)`; it equals `76466200>0` at `k=320` and is increasing
thereafter. Hence the single analytic estimate

```text
P_n<=4/(n+3)^2 for every n>=317                     (11.30.14)
```

is sufficient to close the complete `lambda=-100` order-four tail and is
proved in Lemma 11.31 below. The repaired Arb
prefix proves the stronger finite bound

```text
P_n<2/(n+3)^2, 0<=n<=316,                           (11.30.15)
```

on all 317 available rows. Equation (11.30.14), not extrapolation of
(11.30.15), is the analytic tail theorem proved below.

### Lemma 11.31: Order-Four First-Summand Curvature Reduction

Let `B_1(t)` be the continuous first-summand curvature at `lambda=-100`,
and set

```text
x_1(t)=exp(-B_1(t)),
d_1(t)=1-x_1(t),
ell_1(t)=log d_1(t),
J_1(t)=2B_1(t)-ell_1(t-1)+2ell_1(t)-ell_1(t+1).     (11.31.1)
```

The continuous order-three gap centered at `t` is

```text
g_1(t)=d_1(t)^2-x_1(t)^2*d_1(t-1)*d_1(t+1).
```

Equation (11.31.1) gives the exact stable factorization

```text
g_1(t)=d_1(t)^2*(1-exp(-J_1(t))).                    (11.31.2)
```

Thus the cancellation between two terms of order `t^-2` is represented by
the positive coordinates `d_1(t)` and `J_1(t)`, each of order `t^-1`.

Put

```text
phi(z)=1/(exp(z)-1),
chi(z)=exp(z)/(exp(z)-1)^2.
```

For `ell=log(1-exp(-B))`, direct differentiation gives

```text
ell'=phi(B)*B',
ell''=phi(B)*B''-chi(B)*(B')^2.                     (11.31.3)
```

For `r=0,1,2`, differentiation of (11.31.1) gives

```text
J_1^(r)(t)
 =2B_1^(r)(t)-ell_1^(r)(t-1)+2ell_1^(r)(t)-ell_1^(r)(t+1).
```

Consequently, wherever `J_1(t)>0`,

```text
K_1(t):=(log g_1(t))''
 =2ell_1''(t)+phi(J_1(t))*J_1''(t)
  -chi(J_1(t))*(J_1'(t))^2.                         (11.31.4)
```

There is a stable same-point localization of (11.31.4). Define

```text
j_0=2B_1-ell_1'',
j_1=2B_1'-ell_1''',
j_2=2B_1''-ell_1'''',

E_r=(1/12)*sup_(|s|<=1)|ell_1^(r+4)(t+s)|,
r=0,1,2.                                             (11.31.4a)
```

The tent representation of the centered second difference and Taylor's
theorem give

```text
|J_1^(r)(t)-j_r(t)|<=E_r(t), r=0,1,2.               (11.31.4b)
```

Both `phi` and `chi` decrease on `(0,infinity)`. Therefore, if `j_0>E_0`,

```text
K_1(t)<=U(t),

U(t)=2ell_1''
 +phi(j_0-E_0)*max(j_2+E_2,0)
 -chi(j_0+E_0)*max(|j_1|-E_1,0)^2.                  (11.31.4c)
```

Consequently the pointwise inequalities

```text
j_0(t)>E_0(t), U(t)<=7/(2t^2)                       (11.31.4d)
```

are sufficient for the remaining curvature ceiling. Unlike direct enclosure
of `B_1(t-1),B_1(t),B_1(t+1)`, (11.31.4a)-(11.31.4d) preserve the common
real parameter in all leading terms; only explicit fourth- through
sixth-derivative envelopes retain shifted dependence.

The required positivity of `J_1` has a quantitative lower bound at the
integer points used below. If

```text
M_j=M_j^(1)*(1+delta_j), 0<=delta_j<=2/j^6,
```

and `e_j=log(1+delta_j)`, then

```text
|B_j-B_j^(1)|
 <=a_j:=2*((j-1)^(-6)+2*j^(-6)+(j+1)^(-6)).         (11.31.5)
```

The scaled-defect anchor and the pointwise cone give, for integer `m>=319`,

```text
251/(250*(2m+1))<=d_m<=2/(2m+1).
```

Using `d<=-log(1-d)<=d/(1-d)` and (11.31.5) gives

```text
1/(2m+1)<=B_1(m)<=3/(2m-1), m>=319.                (11.31.6)
```

The inequalities needed to absorb `a_m` are much weaker than their exact
values; for example `a_m<=1/(250*(2m+1))` and `a_m<=1/(2m-1)` follow at
`m=319` and strengthen thereafter.

Two additional interval-Simpson blocks on

```text
0.925<=u<=0.926, 0.926<=u<=0.9264
```

prove `H'''<0` and the positive shifted scaled-curvature buffer below the
old compact endpoint. Their minimum outward-rounded buffer is greater than
`5.59e-4`, and `V'(x(0.925))<316`. Combined with the compact and ray
certificates used in Lemma 11.17, this proves `B_1'(t)<0` and
`((2t+1)B_1(t))'>0` already for `t>=317`. The cumulant bridge proves, for
real `t>=319`,

```text
B_1(t)-B_1(t+1)>=1/(4t^2).                          (11.31.7)
```

Using the floor and ceiling integers in (11.31.6), monotonicity gives

```text
1/(2t+3)<=B_1(t)<=3/(2t-3), t>=319.                (11.31.8)
```

For `t>=319`, scaled-curvature growth between `t-1` and `t` gives

```text
B_1(t-1)-B_1(t)<=2B_1(t)/(2t-1).
```

Since `phi` is decreasing and `phi(B)<=1/B`,

```text
ell_1(t-1)-ell_1(t)<=2/(2t-1).                      (11.31.9)
```

Also `exp(z)<=1/(1-z)` for `0<=z<1`. Equations (11.31.7)-(11.31.8)
therefore give

```text
ell_1(t)-ell_1(t+1)
 >=phi(B_1(t))/(4t^2)
 >=(t-3)/(6t^2).                                    (11.31.10)
```

Combining (11.31.1) and (11.31.8)-(11.31.10) yields

```text
J_1(t)
 >=2/(2t+3)-2/(2t-1)+(t-3)/(6t^2)
 >=1/(7t), t>=319.                                  (11.31.11)
```

After `t=319+n`, the numerator of the last comparison is

```text
4n^3+3412n^2+955637n+87486770,
```

which is positive coefficientwise for `n>=0`. As an independent cross-check,
the repaired Arb coefficients and (11.31.5) certify

```text
J_1(319)-1/(7*319)>0.002695906801650624679.          (11.31.12)
```

Thus

```text
J_1(t)>=1/(7t) for every real t>=319.               (11.31.12a)
```

The remaining global first-summand curvature estimate is

```text
K_1(t)<=7/(2t^2) for every real t>=319.             (11.31.13)
```

A deterministic paired interval-Simpson certificate proves a substantial
real-parameter part of (11.31.13). On `0.9255<=u<=2.00002`, `107,452`
adjacent mode tiles of width `10^-5` enclose `H^(2),...,H^(8)` with
outward-rounded Arb arithmetic. The localized inequality
(11.31.4a)-(11.31.4d) assembles these into `1,073` positive central blocks
covering

```text
319<=t<=V'(2)=37850.3222102113848162951108925613....  (11.31.13a)
```

The weakest curvature margin is
`1.14426650120038260530285332880335545e-10`, and the largest certified
upper bound for `t^2 U(t)` is
`3.337233557567905543142604940455149<7/2`. Thus (11.31.13a) is an interval
theorem. The finite exact-corridor and analytic-ray arguments below complete
the remaining mode tail.

The ray has an exact higher-cumulant contract. Put `q=pi*exp(4u)`,
`epsilon=q^(-1/2)`, and write the standardized potential as

```text
R_epsilon(y)=sum_(r=3)^8 L_r*epsilon^(r-2)*y^r/r!.
```

Exact tilted-Gaussian partition algebra through `epsilon^6` gives the leading
signature

```text
kappa_3~-epsilon, kappa_4~2epsilon^2, kappa_5~-6epsilon^3,
kappa_6~24epsilon^4, kappa_7~-120epsilon^5, kappa_8~720epsilon^6.
```

The complete formal polynomials are recorded in the Gaussian-cumulant ray
target. A sufficient explicit corridor proposed there is

```text
2/5<=q*(kappa_2-1)<=4/5,
1<=(-1)^r*kappa_r*q^(r/2-1)/(r-2)!<=c_r, r=3,...,8,
(c_3,...,c_8)=(6/5,27/20,3/2,17/10,2,5/2).
```

Assuming these inequalities throughout each full `t+-2` collar, Arb boxes
clear (11.31.4d) at seven representative ray modes from `u=2` to `u=20`.
This is a conditional finite compatibility test, not a ray theorem.

The epsilon-six formal model itself is now closed globally. A deterministic
`1,800,000`-block Arb cover proves all seven formal corridors on `2<=u<=20`,
with weakest margin
`0.0184369217593030043644830584526062`. On `u>=20`, coefficient-positive
leading-model gates give a `1/(10u)` corridor buffer, fourteen exact
potential-jet sign gates prove

```text
|V^(r)-qP_r(u)|<=1000u^8, r=2,...,8,
```

and explicit polynomial norms give

```text
|R_r^[6]-F_r|<=22000000u^6/q<1/(20u).
```

Thus the exact epsilon-six formal cumulant polynomial satisfies the candidate
corridors for every `u>=2`.

The first omitted formal layer has now also been closed globally. Exact
tilted-Gaussian algebra through `epsilon^8`, audited term for term against all
42 stored epsilon-six coefficients, gives

```text
scaled(kappa_r^[8]-kappa_r^[6])=q^-3 C_r, r=2,3,4,
scaled(kappa_r^[8]-kappa_r^[6])=q^-2 C_r, r=5,6,
scaled(kappa_r^[8]-kappa_r^[6])=q^-1 C_r, r=7,8.
```

A centered sixth-order Arb Taylor cover of 1,800 blocks proves explicit signed
bounds for all seven `C_r` on `2<=u<=20`. On `u>=20`, fourteen leading buffer
gates, four new order-nine/ten potential-jet gates, and polynomial norm
transfer prove the same bounds analytically. In particular,

```text
-1/20<C_2<0, -3/20<C_3<0, -1/2<C_4<0,
0<C_5<2, 0<C_6<16/5, 0<C_7<5, 0<C_8<61/10.
```

The epsilon-eight correction is therefore smaller than `1/1000` throughout
`2<=u<=20` and smaller than `1/(100u)` on `u>=20`. The formal recurrence has
also been extended through `epsilon^10` and audited against all 56 stored
epsilon-eight coefficients. Its second-next scaled hierarchy is `q^-4` in
orders two through four, `q^-3` in orders five and six, and `q^-2` in orders
seven and eight. A 3,600-block centered Arb Taylor cover and a coefficient-
positive asymptotic transfer prove the signed second-next bounds globally.

The complete epsilon-ten correction is below `10^-7` on `2<=u<=20` and below
`1/(1000u)` on `u>=20`. Exact composition therefore reduces the density
theorem to the sharper sufficient budgets

```text
scaled |kappa_r-kappa_r^[10]|<9/1000,   2<=u<=20,
scaled |kappa_r-kappa_r^[10]|<1/(100u), u>=20.
```

The reserved final corridor margins are `79999/10000000` on the finite segment
and `29/(1000u)` on the ray.

Those exact-density budgets are now proved. The cancellation-preserving
complex-disk factorization is

```text
A_u(z)=exp(-z^2/2)/sqrt(2*pi) * integral_R exp(z*y-W_u(y))dy,
K_u(z)=z^2/2+log A_u(z)-log A_u(0).
```

Exact partition/logarithm recurrences through epsilon ten pass 70 cumulant
coefficient audits. Formal partition nonvanishing, relative logarithm bounds,
and Cauchy's estimate reduce all seven cumulant errors to

```text
|A_u-P_u^[10]|<1/(100000*q^3),       2<=u<=20,
|A_u-P_u^[10]|<1/(20000*u*q^3),      u>=20.          (11.31.13b)
```

Use the adaptive cutoff `Y=1+sqrt(32 log q)`. A degree-thirty formal-density
majorant and exact Gaussian hazards close both formal tails. Two positive-
coefficient curvature inequalities give

```text
(39/10)*v^2*q_v<=V''(x_v),   V''(x_u)<=5*u^2*q,
W_u''(y)>59319/100000>1/2,   |y|<=Y,
```

and global potential monotonicity then closes both exact-density tails.

The remaining central subtraction preserves rather than discards the Gaussian
cancellations. Centered Arb Taylor covers of 5,400 partition blocks prove

```text
||Z_11||_1<2, ||Z_12||_1<2,
||Z_13||_1<21/10, ||Z_14||_1<12/5,
```

while 5,430 shifted-jet blocks give `L_17<4000` on the finite collar. A
222-term Bell-polynomial Taylor remainder retains the full `q^(-15/2)` factor,
and an exact seventeenth-order potential remainder supplies the final density
comparison. The resulting central budget ratios are

```text
finite ratio <0.363964,
ray ratio    <3.29*10^-63.
```

Together with all four tails, (11.31.13b) is proved. The complex-disk Cauchy
contract therefore proves the simultaneous exact-minus-epsilon-ten cumulant
budgets through order eight, and the reserve arithmetic proves all seven
alternating-factorial exact cumulant corridors for every `u>=2`.

The exact corridor-to-curvature composition is now also closed. First, a
nonuniform rational cover retaining the common `q`, `t`, curvature,
Hurwitz-zeta, and corridor dependence proves (11.31.4d) on `2<=u<=20`.
It has 20,700 mode blocks and 41,400 shifted-collar gates; its largest
certified upper bound is

```text
t^2*U(t)<3.491488695857048<7/2.                     (11.31.13c)
```

The remaining ray admits a short dimensionless proof. On the enlarged collar
`u>=19`, coefficient-positive potential geometry and monotonicity of `q/u`
give

```text
q>=10^33,
u*q<=t<=(201/100)*u*q,
V''(x_u)>=(361/100)*u^2*q,
1/t<=10^-30.                                        (11.31.13d)
```

Moreover `t(20)-t(19)>2`, so every `t+s`, `|s|<=1`, belonging to a central
mode `u>=20` remains in that enlarged collar. For `2<=r<=8`, put

```text
x_r=(-1)^r*t^(r-1)*H^(r)/(r-2)!.
```

Convex midpoint quadrature and the integral lower bound for the Hurwitz zeta
function give

```text
(t/(t+1/2))^(r-1)
 <=(r-1)*t^(r-1)*zeta(r,t+1/2)<=1.
```

Combining this with (11.31.13d) and the exact cumulant corridors proves

```text
0<x_r<=1 (2<=r<=8), x_2>=97/100, x_3>=24/25.       (11.31.13e)
```

Write

```text
ell=log(1-exp(-B))=log B+R(B),
R(B)=log((1-exp(-B))/B).
```

The convergent product expansion

```text
R(B)=-B/2+sum_(m>=1)(-1)^(m+1)*zeta(2m)
      *B^(2m)/(m*(2*pi)^(2m))
```

and the elementary bounds `zeta(2m)<2`, `2*pi>6` imply
`|R^(j)(B)|<1` for `1<=j<=6`. The exact partial-Bell identity

```text
B_(n,k)(1!,2!,...)=(n!/k!)*binom(n-1,k-1)
```

then gives

```text
t^2*ell''<=23/20,
|t^r*ell^(r)|<30000, 2<=r<=6,
t^(r+1)*E_r<1/1000, r=0,1,2.                       (11.31.13f)
```

Consequently the scaled localized quantities satisfy

```text
t*(j_0-E_0)>193/100,       t*(j_0+E_0)<201/100,
t^3*max(j_2+E_2,0)<401/100,
t^2*max(|j_1|-E_1,0)>191/100.                       (11.31.13g)
```

For `0<z<=1/1000`, use

```text
z/(exp(z)-1)<=1,
z^2*exp(z)/(exp(z)-1)^2>=exp(-z)>=1-z>=999/1000.
```

The negative square term in (11.31.4c) is essential and yields the strict
rational comparison

```text
t^2*U(t)
 <2*(23/20)+401/193-(999/1000)*(191/201)^2
 =3011223637/866377000
 <7/2,                                                (11.31.13h)
```

with reserve `21095863/866377000`. Thus (11.31.4d) holds on `u>=20`.
Together with (11.31.13c) and the compact theorem (11.31.13a), this proves
the formerly open global estimate

```text
K_1(t)<=7/(2t^2) for every real t>=319.              (11.31.13)
```

The exact consequence of the global estimate is as follows. For
`k=n+3>=320`, the tent identity gives

```text
P_n^(1)
 =integral_[-1,1](1-|s|)*K_1(k+s) ds
 <=(7/2)*[-log(1-1/k^2)]
 <=7/(2*(k^2-1))
 <=18/(5k^2).                                       (11.31.14)
```

The final rational margin is

```text
18/(5k^2)-7/(2*(k^2-1))
 =(k^2-36)/(10k^2*(k^2-1))>=0.
```

The full-kernel transfer is already unconditional. Along the segment from
`B_1` to `B`, (11.31.5) and the coefficient cone give `B_theta,j>1/(4j)`.
Hence

```text
|ell_j-ell_j^(1)|<=32j/(j-1)^6,
|J_j-J_j^(1)|<=176j/(j-2)^6<=1/(14j).              (11.31.15)
```

Together with (11.31.11)-(11.31.12), the derivative bounds for
`log(1-exp(-J))` imply

```text
|log g_j-log g_j^(1)|<=2528j^2/(j-2)^6, j>=319.
```

Therefore

```text
|P_n-P_n^(1)|
 <=10112*(k+1)^2/(k-3)^6
 <=2/(5k^2), k>=320.                                (11.31.16)
```

After `k=320+n`, the last inequality is equivalent to positivity of

```text
n^6+1902n^5+1482055n^4+604691300n^3
 +135889991935n^2+15877422036942n+748002501678169.
```

All coefficients are positive. Consequently (11.31.13) combines with
(11.31.14)-(11.31.16) to give

```text
P_n<=4/(n+3)^2, n>=317.                              (11.31.17)
```

This is exactly the tail required in Lemma 11.30.

### Lemma 11.32: All-Shift Contiguous Order-Four Entry at Lambda=-100

For `k=n+3>=320`, the scaled-defect anchor gives

```text
-3*log(x_k)>3*d_k>=753/(250*(2*k+1)).                (11.32.1)
```

The penalty ceiling (11.31.17) lies strictly below this buffer:

```text
4/k^2<753/(250*(2*k+1)), k>=320.                    (11.32.2)
```

After `k=320+m`, the cleared numerator of (11.32.2) is

```text
753*m^2+479920*m+76466200,
```

whose coefficients are positive. Therefore

```text
P_n<-3*log(x_(n+3)), n>=317.                        (11.32.3)
```

By the exact sign equivalence (11.30.11), (11.32.3) proves
`H_(4,n)(-100)>0` throughout the analytic tail. The repaired 1024-bit Arb
prefix in Lemma 11.30 proves the same sign for all `0<=n<=316`. Hence

```text
H_(4,n)(-100)>0 for every integer n>=0.             (11.32.4)
```

This is a complete all-shift contiguous order-four entry theorem at one heat
parameter. It does not prove forward order-four invariance through
`lambda=0`, arbitrary-column order four, compound order five or higher, the
all-order signed-Hankel/Jensen bridge, PF-infinity, RH, or `Lambda<=0`. The
next live target is the exact order-four compound-flow law and a weighted
infinite-index maximum principle starting from (11.32.4).

### Lemma 11.33: Cooperative Order-Four Flow Reduction

The normalized coefficients obey

```text
A_j'=(4j+2)A_(j+1).                                 (11.33.1)
```

Let `delta` denote the unweighted shift derivation
`delta(A_j)=A_(j+1)`. Direct multilinearity of a contiguous Hankel
determinant gives the affine-weight identity

```text
H_(m,n)'=(4(n+2m-2)+2)*delta(H_(m,n)).              (11.33.2)
```

For orders three and four, the relevant Plucker identity is

```text
H_(3,n+1)*delta(H_(4,n))
 =H_(3,n)*H_(4,n+1)
  +delta(H_(3,n+1))*H_(4,n).                       (11.33.3)
```

Put `T_n=-H_(3,n)>0` and `Q_n=H_(4,n)`. Equations
(11.33.1)--(11.33.3) give the exact one-sided system

```text
Q_n'=a_n*Q_(n+1)+b_n*Q_n,                          (11.33.4)

a_n=(4n+26)*T_n/T_(n+1)>0,
b_n=((4n+26)/(4n+22))*(log T_(n+1))'.              (11.33.5)
```

In particular, at `Q_n=0`, the derivative has the sign of `Q_(n+1)`.
Thus the local order-four boundary algebra is cooperative and does not
require an order-five sign.

There is a bounded stable normalization. With the positive order-three gaps
from Lemma 11.30, define

```text
F_n=G_(n+1)^2-x_(n+3)^3*G_n*G_(n+2),

S_n=A_n^4*rho_n^12*x_(n+1)^8*x_(n+2)^4/d_(n+3).
```

Exact coefficient-ratio elimination gives

```text
H_(4,n)=S_n*F_n, S_n>0.                            (11.33.6)
```

Inside the completed ratio and order-three cones, `0<G_j<=1`; hence
`|F_n|<=1`. Positive rescaling of (11.33.4) yields

```text
F_n'=alpha_n*F_(n+1)+beta_n*F_n,                   (11.33.7)

alpha_n
 =(4n+26)*rho_(n+2)*x_(n+3)^4
  *(d_(n+3)/d_(n+4))*(G_n/G_(n+1))>0,

beta_n
 =((4n+26)/(4n+22))*(log T_(n+1))'-(log S_n)'.     (11.33.8)
```

Set `z_n=F_n/(n+1)`. Then `z_n->0` uniformly on every compact heat interval,
and

```text
z_n'
 =alpha_n*((n+2)/(n+1))*z_(n+1)+beta_n*z_n.        (11.33.9)
```

Consequently the standard exponential weighted-minimum argument would
propagate (11.32.4) through every finite heat time after proving only

```text
sup_(-100<=lambda<=L,n>=0)
 [beta_n+alpha_n*(n+2)/(n+1)]<infinity              (11.33.10)
```

for each finite `L>=-100`. Lemma 11.33 proves the local cooperative flow,
positive boundary orientation, bounded tail coordinate, and exact infinite-
system reduction. By itself it does not prove (11.33.10), forward order-four
invariance, arbitrary-column order four, higher compounds, PF-infinity, RH,
or `Lambda<=0`. At this stage the direct weighted-infimum route has the sole
blocker (11.33.10); Theorem 11.37 later bypasses it with a uniform eventual
tail and finite confinement.

### Lemma 11.34: Lambda-Zero Eventual Order-Four Positivity

Let `gamma(k)` be the standard Xi coefficients

```text
xi(1/2+z)=sum_(k>=0) gamma(k)*z^(2k)/k!.
```

The identities `H_0(x)=xi((1+i*x)/2)/8` and
`A_k(0)=2k!/(2k)! integral_0^infinity Phi(u)u^(2k)du` give

```text
A_k(0)=gamma(k)/4^(k+1),                            (11.34.1)

H_(4,n)[A(0)]=4^(-4n-16) H_(4,n)[gamma].           (11.34.2)
```

Use Theorem 2.1 of Griffin, Ono, Rolen, Thorner, Tripp, and Wagner,
`Jensen Polynomials for the Riemann Xi Function`. For integers `1<=j<M`, it
gives the convergent ratio expansion

```text
log(gamma(M-j)/gamma(M))
 =-sum_(m>=1) G_m(M)*Delta(M)^(2m-2)*j^m,           (11.34.3)

Delta(M)~1/sqrt(2M),
G_m(M)->2^(m-1)/(m(m-1)),
|G_m(M)|<<_C(2C)^m.                                 (11.34.4)
```

In particular, `G_2(M)->1`. Put `M=n+6`, `h=Delta(M)^2`, and in the
`(i,j)` entry of the order-four determinant set `r=6-i-j`. The `m=1`
term in (11.34.3) factors positively from rows and columns. The remaining
matrix has entries

```text
K_(i,j)(h)
 =exp(-sum_(m>=2)G_m(M)h^(m-1)(6-i-j)^m).           (11.34.5)
```

Exact truncated-series determinant algebra over all 24 permutations gives

```text
[h^0,...,h^6] det K
 =[0,0,0,0,0,0,768*G_2(M)^6].                      (11.34.6)
```

Thus the terms containing `G_3,...,G_7` cancel from the first nonzero
coefficient. The uniform bound in (11.34.4) makes the contribution from
`m>=8` equal to `O(h^7)` entrywise for the fixed shifts `0<=r<=6`.
Consequently

```text
H_(4,n)[gamma]
 =gamma(M)^4*exp(-12G_1(M))
  *(768G_2(M)^6*Delta(M)^12+o(Delta(M)^12)).        (11.34.7)
```

Every factor outside the parentheses is positive and `G_2(M)->1`. Hence
there exists a finite `N_H4` such that

```text
H_(4,n)(0)>0 for every n>=N_H4.                     (11.34.8)
```

Independently, a direct 24,576-bit Arb expansion of `xi(1/2+z)` gives
outward-rounded coefficient balls through `A_506(0)` and proves

```text
H_(4,n)(0)>0 for every 0<=n<=500.                   (11.34.9)
```

All 501 stable margins are also strictly positive; the smallest occurs at
`n=500`, with

```text
F_500=5.48126835651812768442696180364636...e-20.
```

Lemma 11.34 is an exact eventual theorem plus a rigorous finite prefix. By
itself it does not make `N_H4` explicit and therefore does not prove the
all-shift lambda-zero order-four sign. A direct splice would require an
explicit remainder constant in (11.34.7) and certification of any indices
`501<=n<N_H4`. The uniform compact-heat route below closes the all-shift
contiguous order-four statement without constructing that direct splice.

### Lemma 11.35: Uniform Higher-Theta Suppression On Compact Heat

Write `T=-lambda`, so `0<=T<=100`, and let `A_k^(1)(-T)` denote the first
theta-summand contribution to `A_k(-T)`. Define

```text
delta_k(T)=A_k(-T)/A_k^(1)(-T)-1.                    (11.35.1)
```

After normalizing the first summand, this is the expectation of the
pointwise higher-summand ratio `epsilon` under

```text
dmu_(k,T)(u)=Z_(k,T)^(-1)u^(2k)exp(-T*u^2)Phi_1(u)du,
delta_k(T)=E_(mu_(k,T))[epsilon(U)].                 (11.35.2)
```

The function `epsilon(u)` is decreasing, whereas `u^2` is increasing. Exact
differentiation and the two-copy covariance identity give

```text
delta_k'(T)=-Cov_(mu_(k,T))(epsilon(U),U^2)>=0.      (11.35.3)
```

Hence `0<=delta_k(T)<=delta_k(100)`. The already proved `T=100` low- and
high-region estimates are

```text
B_low(k)=-8k/(5(log k+4/5))+(5/2)log k+1
         +pi(exp(2/5)-1)sqrt(k),
epsilon(a(k))<=17exp(-3pi sqrt(k)),  a(k)=log(k)/8.
```

For every fixed `p>0`,

```text
lim_(k->infinity)(B_low(k)+p log k)log k/k=-8/5,
lim_(k->infinity)(log 17-3pi sqrt(k)+p log k)/sqrt(k)=-3pi.
```

Both limits are negative. Therefore, for every `p>0`, there is a `K_p`
such that

```text
sup_(0<=T<=100) delta_k(T)<=2k^(-p),  k>=K_p.        (11.35.4)
```

For `e_k(T)=log(1+delta_k(T))`, every fixed forward or backward difference
stencil satisfies

```text
|Delta^m e_k(T)|<=2^m max_(0<=j<=m)e_(k+j)(T).      (11.35.5)
```

Thus all fixed local logarithmic differences of the higher-theta correction
are uniformly superpolynomial in `k` on `0<=T<=100`. In particular, they are
negligible to every fixed order needed by the seven-node order-four ratio
expansion.

### Lemma 11.36: Uniform First-Summand Heat-Tilt Asymptotics

Theorem 5.2 of Cormac O'Sullivan, `Zeros of Jensen polynomials and
asymptotics for the Riemann xi function`, gives arbitrary-order saddle
asymptotics for

```text
I_alpha(f;n)=integral_1^infinity (log t)^n exp(-alpha*t)f(t)dt.  (11.36.1)
```

Its Section 5 explicitly includes logarithmic Gaussian multipliers. For the
compact family

```text
f_T(t)=exp(-T(log t)^2/16),  0<=T<=100,              (11.36.2)
```

one has the exact suitability quotient

```text
f_T(t(1+x))/f_T(t)
 =exp(-T((v/8)log(1+x)+log(1+x)^2/16)),  t=exp(v).  (11.36.3)
```

On the complex disk used in the published proof, its exponent is bounded
uniformly in `T`. Cauchy estimates therefore make every fixed-order
suitability coefficient and remainder uniform over the compact family. The
exact `t=exp(4u)` substitution writes the first theta summand as

```text
M_k^(1)(T)=C_k(2pi^2 I_pi(t^(5/4)f_T;2k)
                    -3pi I_pi(t^(1/4)f_T;2k)).      (11.36.4)
```

Applying the all-order expansion to both integrals and combining them gives,
for every fixed truncation order `R`,

```text
log R_T^(1)(k)
 =-T*w^2/16+sum_(r=1)^(R-1)Q_(r,T)(w)/k^r
   +O_R(w^(3R)/k^R),                                (11.36.5)
w=W(2k/pi),
```

uniformly for `0<=T<=100`. Here `R_T^(1)(k)` is the first-summand heat
ratio, so the heat-independent leading factor has cancelled.

The Lambert recurrence and the integral formula for finite differences are

```text
w'(k)=w/(k(1+w)),
Delta^m F(k)=integral_[0,1]^m F^(m)(k+s_1+...+s_m)ds. (11.36.6)
```

Exact differentiation through order seven yields
`d^m(w^2)/dk^m=O(w/k^m)`. Taking `R>m` in (11.36.5), and differentiating
the polynomial correction and remainder bounds, proves

```text
Delta_k^m log R_T^(1)(k)=O(log(k)/k^m),
2<=m<=7, uniformly for 0<=T<=100.                   (11.36.7)
```

The same statement holds for backward stencils. This is the complete
first-summand heat correction required by the compact-heat order-four tail.

### Theorem 11.37: Contiguous Order-Four Forward Invariance Through Zero

Combine the lambda-zero graded Xi ratio expansion of Lemma 11.34 with
Lemmas 11.35 and 11.36. On the seven determinant nodes, exact Newton
interpolation writes a normalized logarithmic ratio as

```text
q_M(j)=sum_(m=1)^6 D_m(M) binom(j,m),  0<=j<=6.     (11.37.1)
```

The coefficient of `j^r` in (11.37.1) uses only `D_m` with `m>=r`, and its
diagonal term is `D_r/r!`. From (11.36.7) and Lemma 11.35, the heat
correction at degree `r` is `O(log M/M^r)`. Since `h(M)~1/(2M)`,

```text
c_r(T,M)/h(M)^(r-1)=O(log M/M)=o(1),  2<=r<=6,      (11.37.2)
```

uniformly in `T`. The degree-seven ratio contract therefore has the same
leading `G_2` limit, uniformly:

```text
G_2(T,M)->1,  0<=T<=100.                             (11.37.3)
```

The exact determinant calculation (11.34.6) is universal in the graded
ratio coefficients. It now gives

```text
det K_T=768G_2(T,M)^6h(M)^6+o(h(M)^6)               (11.37.4)
```

uniformly on `0<=T<=100`. All factors removed in the row and column
normalization are positive. Compact uniformity and (11.37.3) imply the
existence of one finite, possibly ineffective, integer `N` such that

```text
H_(4,n)(lambda)>0 for n>=N and -100<=lambda<=0.      (11.37.5)
```

Let `Q_n` be the stable order-four margin from Lemma 11.33; it has the same
sign as `H_(4,n)` and obeys

```text
Q_n'=alpha_n Q_(n+1)+beta_n Q_n,  alpha_n>0.         (11.37.6)
```

For each fixed `n`, variation of constants gives

```text
Q_n(lambda)=E_n(lambda)(Q_n(-100)
 +integral_(-100)^lambda E_n(s)^(-1)alpha_n(s)Q_(n+1)(s)ds),   (11.37.7)
```

where `E_n(lambda)>0`. Equation (11.37.5) supplies the positive boundary
`Q_N` throughout the heat interval, while the all-shift entry theorem gives
`Q_n(-100)>0` for every `n`. Backward induction in (11.37.7), from `N-1`
to zero, proves

```text
H_(4,n)(lambda)>0
for every integer n>=0 and every -100<=lambda<=0.    (11.37.8)
```

In particular, every contiguous order-four Hankel minor is positive at
`lambda=0`. The theorem does not cover noncontiguous order-four minors,
compound order at least five, PF-infinity, the all-degree Jensen bridge, RH,
or `Lambda<=0`.

Machine-audited companions:

```text
outputs/jensen_window_pf_uniform_superpolynomial_first_summand_dominance.md
outputs/jensen_window_pf_uniform_first_summand_heat_tilt_asymptotic_theorem.md
outputs/jensen_window_pf_compound_order4_uniform_heat_forward_invariance_certificate.md
```

### Theorem 11.38: Arbitrary-Column Transfer Through Order Four

Put

```text
epsilon_k=(-1)^binom(k,2).                           (11.38.1)
```

The completed lambda-zero layers give

```text
epsilon_k H_(k,s)(0)>0
for every s>=0 and k=1,2,3,4.                       (11.38.2)
```

Fix integers `n,N>=0` and reverse the first `N+1` Hankel columns in the
four-row block:

```text
B_(i,q)^(n,N)=A_(n+i+N-q)(0),
0<=i<=3, 0<=q<=N.                                   (11.38.3)
```

Take a solid `k` by `k` minor with row start `r` and column start `c`, where
`1<=k<=4`. Reversing its columns gives the exact identity

```text
det B[r:r+k,c:c+k]
 =epsilon_k H_(k,n+r+N-c-k+1)(0).                  (11.38.4)
```

The mapped shift is nonnegative because `c<=N-k+1`. Therefore (11.38.2)
makes every solid minor of `B^(n,N)` positive. In particular, every initial
minor, meaning a solid minor touching the first row or first column, is
positive.

The rectangular initial-minor criterion of Gasca and Pena, `Total positivity
and Neville elimination`, Linear Algebra and its Applications 165 (1992),
25-44, states that a real rectangular matrix is strictly totally positive if
and only if all its initial minors are positive. Applying it to (11.38.3)
shows that every minor of `B^(n,N)` is positive.

Now let `0<=j_1<...<j_k<=N`. The columns
`N-j_k<...<N-j_1` occur in increasing order in `B`, and direct reversal gives

```text
det B[0:k | N-j_k,...,N-j_1]
 =epsilon_k R_(k,n)(j_1,...,j_k).                  (11.38.5)
```

The left side is positive. Taking `N=j_k` proves, in particular,

```text
R_(4,n)(j_1,j_2,j_3,j_4)>0
for every n>=0 and 0<=j_1<j_2<j_3<j_4.             (11.38.6)
```

The argument does not depend on the number four. For every fixed `m`, an
`m`-row version of (11.38.3) proves the structural implication

```text
[epsilon_k H_(k,s)>0 for all s and 1<=k<=m]
 =>[epsilon_k R_(k,n)(j_1,...,j_k)>0
    for all n, all increasing columns, and 1<=k<=m]. (11.38.7)
```

All lower layers in (11.38.2) are essential. The positive rational sequence

```text
(10,9,29,18,21,25,3,16)
```

has `H_(4,0)=288076>0` and `H_(4,1)=264875>0`, but
`R_(4,0)(0,1,3,4)=-231169<0`; it also has the wrong lower sign
`H_(2,0)=209>0`. Thus contiguous order four alone cannot be promoted by this
criterion.

Theorem 11.38 completes the consecutive-row, arbitrary-column signed-Hankel
structure through order four. It also removes arbitrary-column transfer as a
separate burden at every future fixed order. The first new layer is contiguous
order five; PF-infinity, the all-degree Jensen bridge, RH, and `Lambda<=0`
remain open.

Machine-audited companion:

```text
outputs/jensen_window_pf_order4_noncontiguous_total_positivity_transfer.md
```

### Theorem 11.39: Uniform Order-Five Tail And Endpoint Reduction

The published Xi ratio theorem and O'Sullivan's all-order suitable-multiplier
expansion are not restricted to the seven differences used at order four.
For the compact family `f_T(t)=exp(-T(log t)^2/16)`, the same Cauchy audit
through degree eleven and exact differentiation of `w=W(2k/pi)` give

```text
Delta_k^r log(A_k^(1)(-T)/A_k^(1)(0))
 =O(log(k)/k^r), 2<=r<=11,
uniformly for 0<=T<=100.                             (11.39.1)
```

Lemma 11.35 removes the higher theta summands to every fixed order. Exact
nine-node Newton interpolation then transfers (11.39.1) to the graded ratio
coefficients on the order-five determinant stencil.

Put `M=n+8` and `h=Delta(M)^2`. For entries

```text
K_(i,j)(h)=exp(-sum_(r>=2)G_r h^(r-1)(8-i-j)^r),
0<=i,j<=4,                                           (11.39.2)
```

exact truncated-series algebra over all 120 permutations gives

```text
[h^0,...,h^10]det K
 =[0,0,0,0,0,0,0,0,0,0,294912G_2^10].              (11.39.3)
```

Here `294912=2^10(1!2!3!4!)`; all terms containing `G_3,...,G_11`
cancel from the first nonzero coefficient. Since the heat correction preserves
`G_2(T,M)->1` uniformly, (11.39.3) proves the existence of one finite,
possibly ineffective, integer `N_5` such that

```text
H_(5,n)(lambda)>0
for n>=N_5 and -100<=lambda<=0.                      (11.39.4)
```

The exact Desnanot-Jacobi identity at the new layer is

```text
H_(5,n)H_(3,n+2)
 =H_(4,n)H_(4,n+2)-H_(4,n+1)^2.                    (11.39.5)
```

The completed signs `H_(3,n)<0` and `H_(4,n)>0` therefore give the endpoint
equivalence

```text
H_(5,n)(-100)>0
 iff H_(4,n+1)(-100)^2
      >H_(4,n)(-100)H_(4,n+2)(-100).                (11.39.6)
```

Let `delta(A_j)=A_(j+1)`. Affine-weight differentiation and the adjacent
Plucker identity give

```text
H_(5,n)'=(4n+34)delta(H_(5,n)),

H_(4,n+1)delta(H_(5,n))
 =H_(4,n)H_(5,n+1)+delta(H_(4,n+1))H_(5,n).         (11.39.7)
```

Thus, for `Q_n=H_(5,n)`,

```text
Q_n'=a_n Q_(n+1)+b_n Q_n,
a_n=(4n+34)H_(4,n)/H_(4,n+1)>0,
b_n=((4n+34)/(4n+30))(log H_(4,n+1))'.              (11.39.8)
```

The order-five flow is cooperative and requires no order-six sign. If the
endpoint condition in (11.39.6) is proved for every `n`, the uniform tail
(11.39.4) provides a positive boundary at `N_5`; variation of constants and
finite backward induction then give

```text
[H_(5,n)(-100)>0 for every n]
 =>[H_(5,n)(lambda)>0
    for every n and every -100<=lambda<=0].          (11.39.9)
```

The endpoint condition is genuinely new. The positive rational sequence

```text
(1,1,1/2,1/6,1/24,1/120,1/720,1/5040,1/42000)
```

has every available strict signed contiguous minor through order four, but

```text
H_(5,0)=-1/3657830400000<0,
H_(4,1)^2-H_(4,0)H_(4,2)
 =-1/3792438558720000000<0.                          (11.39.10)
```

Therefore Theorem 11.39 is an unconditional uniform-tail theorem and an exact
conditional propagation theorem, not an all-shift order-five theorem. Its
sole open input is the all-shift endpoint inequality (11.39.6) at
`lambda=-100`. PF-infinity, the all-degree Jensen bridge, RH, and `Lambda<=0`
remain open.

Machine-audited companion:

```text
outputs/jensen_window_pf_compound_order5_uniform_tail_flow_reduction.md
```

### Corollary 11.40: Rigorous Order-Five Entry Prefix

Retain the contraction and stable order-four coordinates

```text
x_k=A_(k-1)A_(k+1)/A_k^2,  d_k=1-x_k,
G_n=d_(n+2)^2-x_(n+2)^2d_(n+1)d_(n+3),
F_n=G_(n+1)^2-x_(n+3)^3G_nG_(n+2).                 (11.40.1)
```

Define the stable order-five margin

```text
J_n=d_(n+3)d_(n+5)F_(n+1)^2
    -x_(n+4)^4d_(n+4)^2F_nF_(n+2).                 (11.40.2)
```

Exact coefficient-ratio elimination gives

```text
H_(5,n)=W_n J_n,                                    (11.40.3)

W_n=A_n^5 rho_n^20 x_(n+1)^15 x_(n+2)^10 x_(n+3)^5
    /(d_(n+3)d_(n+4)^2d_(n+5)G_(n+2)).             (11.40.4)
```

Every factor in `W_n` is positive in the completed ratio, order-three, and
order-four cones. Thus `J_n` has exactly the sign of `H_(5,n)`.

Seven cached coefficient sources, merged in recorded precedence order and
hashed in the certificate, give 1024-bit outward-rounded Arb balls for

```text
A_k(-100)>0, 0<=k<=324.                              (11.40.5)
```

Direct interval evaluation of (11.40.1)--(11.40.2) proves

```text
J_n(-100)>0,
relative_margin_n(-100)>0,
H_(5,n)(-100)>0,
for every 0<=n<=316.                                (11.40.6)
```

The weakest rows occur at `n=316`, with

```text
J_316=[1.195506752987560996592413174e-45 +/- 8.91e-73],
relative_316
 =[0.00626902754512557813924681177 +/- 7.83e-30].   (11.40.7)
```

Consequently the endpoint target in Theorem 11.39 is reduced to

```text
J_n(-100)>0 for every n>=317.                        (11.40.8)
```

Corollary 11.40 is a rigorous finite prefix, not an analytic tail or an
all-shift order-five theorem. The eventual theorem (11.39.4) has a
non-effective threshold and does not by itself splice (11.40.6) to all shifts.

Machine-audited companion:

```text
outputs/jensen_window_pf_compound_order5_m100_prefix_certificate.md
```

### Lemma 11.41: Scalar Curvature Target For The Order-Five Tail

Set `k=n+4` and, inside the completed positive lower cone, define

```text
C_n=log(F_nF_(n+2)/F_(n+1)^2)
    +log(d_k^2/(d_(k-1)d_(k+1))).                   (11.41.1)
```

Dividing (11.40.2) by its positive second term and taking logarithms gives the
exact equivalence

```text
J_n>0 iff C_n<-4log x_k.                             (11.41.2)
```

In finite-difference notation,

```text
C_n=Delta^2 log F_n-Delta^2 log d_(n+3).             (11.41.3)
```

The lambda `-100` scaled-defect anchor already proves, for `k>=320`,

```text
-4log x_k>=4d_k>=502/(125(2k+1)).                   (11.41.4)
```

For every `k>=321`, exact rational arithmetic gives

```text
100/k^2<502/(125(2k+1)).                             (11.41.5)
```

Indeed, after putting `k=321+m`, the cleared difference is

```text
502m^2+297284m+43689082,                             (11.41.6)
```

whose coefficients are positive. Consequently the single estimate

```text
C_n<=100/(n+4)^2 for every n>=317                    (11.41.7)
```

would imply `J_n(-100)>0` on the complete missing tail. Together with
Corollary 11.40 it would prove all-shift entry at `lambda=-100`; Theorem 11.39
would then propagate contiguous order five through `lambda=0`, and Theorem
11.38 would supply arbitrary columns.

At the finite splice, rigorous evaluation gives

```text
320^2 C_316=3.5869277550969014082...<100.            (11.41.8)
```

Thus the proposed constant has a factor above 27 of numerical reserve at the
boundary. Lemma 11.41 does not prove (11.41.7); it reduces the entire remaining
order-five endpoint problem to that zeta-specific stable log-curvature ceiling.

Machine-audited companion:

```text
outputs/jensen_window_pf_compound_order5_m100_tail_curvature_reduction.md
```

### Lemma 11.42: First-Summand Bridge For The Order-Five Curvature

Retain the continuous first-summand coordinates `B`, `d=1-exp(-B)`,
`ell=log d`, and

```text
g=d^2-x^2 d(t-1)d(t+1),  h=log g.                 (11.42.1)
```

The second stable gap and its logarithmic coordinate are

```text
R(t)=3B(t)-h(t-1)+2h(t)-h(t+1),
f(t)=g(t)^2-x(t)^3g(t-1)g(t+1)
    =g(t)^2(1-exp(-R(t))),
q(t)=log(f(t)/d(t))
    =2h(t)-ell(t)+log(1-exp(-R(t))).               (11.42.2)
```

Thus, for `k=n+4`, the first-summand version of (11.41.1) satisfies

```text
C_n^(1)=Delta^2 q(k)
       =integral_[-1,1](1-|s|)q''(k+s) ds.         (11.42.3)
```

Writing

```text
phi(z)=1/(exp z-1),
chi(z)=exp z/(exp z-1)^2,                           (11.42.4)
```

exact differentiation, with the negative square retained, gives

```text
q''=2h''-ell''+phi(R)R''-chi(R)(R')^2.             (11.42.5)
```

The completed first-summand gap theorem and the full-kernel bound
`0<=delta_j<=2/j^6` imply the stable floors

```text
min(J_j,J_j^(1))>=1/(8j),  j>=319,
min(R_j,R_j^(1))>=7/(5j),  j>=320.                 (11.42.6)
```

Consequently every logarithmic perturbation is taken away from zero. The
explicit Lipschitz chain in the companion certificate gives, without a
conditional hypothesis,

```text
|C_n-C_n^(1)|<=37/k^2,  k=n+4>=321.                (11.42.7)
```

After clearing denominators and setting `k=321+m`, the reserve in (11.42.7)
is a degree-52 polynomial with all 53 coefficients positive. In particular,
the single continuous estimate

```text
q''(t)<=60/t^2 for every real t>=320               (11.42.8)
```

would imply

```text
C_n^(1)<63/k^2,
C_n<100/k^2,  k>=321.                              (11.42.9)
```

Lemma 11.42 is an exact first/full transfer. It does not assert (11.42.8).

Machine-audited companion:

```text
outputs/jensen_window_pf_compound_order5_first_summand_curvature_bridge.md
```

### Theorem 11.43: Global First-Summand Nested-Curvature Bound

The continuous estimate (11.42.8) holds. Its proof is the union of three
rigorous interval theorems.

First, one common `t+-3` collar and outward-rounded Taylor jets through
`H^(8)` evaluate the nested stable logarithms without subtracting raw
determinants. The `107452`-tile order-four cache, with SHA-256

```text
d721a22738543dd2f62181a31732b26666d13eb3f8c4f1c8c46ead3e84ada4cf,
```

assembles into 36 central blocks and proves

```text
q''(t)<=60/t^2 for 320<=t<=V'(2).                  (11.43.1)
```

The largest certified scaled upper on this compact range is

```text
t^2q''(t)<=22.0266810828795017766530766838....     (11.43.2)
```

Second, 100 additional interval-Simpson tiles close the mode-two collar, and
1850 exact-cumulant-corridor blocks prove

```text
q''(t)<=60/t^2 for every mode 2<=u<=20.            (11.43.3)
```

The largest scaled upper there is

```text
11.9132395331211577148958456918....                (11.43.4)
```

Third, on the enlarged `u>=19` collar put

```text
x_r=(-1)^r t^(r-1)H^(r)/(r-2)!.
```

The proved normalized boxes are

```text
0<x_r<=1, 2<=r<=8,
x_2>=97/100, x_3>=24/25, 1/t<=10^(-30).            (11.43.5)
```

The separation `t(20)-t(19)>6.84*10^36` keeps every `t+-3` point in that
collar. For `b_r=t^(r+1)B^(r)/r!`, exact rational tent bounds give

```text
b_0 in [969/1000,1001/1000],
b_1 in [-1001/1000,-959/1000],
(-1)^r b_r in [0,1001/1000], 2<=r<=6.              (11.43.6)
```

After analytically removing the irrelevant `log(1/t)` constants, the stable
logarithms are evaluated through

```text
R_0(w)=log((1-exp(-w))/w).                          (11.43.7)
```

The exact partial-Bell identity bounds the normalized coefficient correction
of order `r` by `2^(r-1)C*10^(-30)` for `r>=1`, with `C=2` at the defect
layer and `C=100` at the nested layers. A single outward-rounded evaluation
of the resulting dimensionless box proves

```text
t^2q''(t)<=9.15835270172636955976486206055...<10
for every mode u>=20.                              (11.43.8)
```

The ranges (11.43.1), (11.43.3), and (11.43.8) cover every real `t>=320`.
Therefore

```text
q''(t)<=60/t^2 for every real t>=320.               (11.43.9)
```

Here `q` denotes the first-summand coordinate of Lemma 11.42. The theorem is
a rigorous interval composition, not an inference from the finite scout
values.

Machine-audited companions:

```text
outputs/jensen_window_pf_compound_order5_nested_curvature_compact_certificate.md
outputs/jensen_window_pf_compound_order5_nested_curvature_finite_ray_certificate.md
outputs/jensen_window_pf_compound_order5_nested_curvature_asymptotic_ray_certificate.md
```

### Theorem 11.44: All-Shift Order-Five Entry At Lambda Minus One Hundred

For `k=n+4>=321`, Theorem 11.43 and the tent identity give

```text
C_n^(1)
 <=60[-log(1-1/k^2)]
 <=60/(k^2-1)
 <63/k^2.                                          (11.44.1)
```

Combining (11.44.1) with the unconditional transfer (11.42.7) yields

```text
C_n<100/k^2 for every n>=317.                      (11.44.2)
```

Equations (11.41.2)--(11.41.5) now imply

```text
J_n(-100)>0,
H_(5,n)(-100)>0,
for every n>=317.                                  (11.44.3)
```

The rigorous prefix (11.40.6) covers the complementary shifts. Hence

```text
H_(5,n)(-100)>0 for every integer n>=0.            (11.44.4)
```

This discharges the former endpoint target
`target_compound_order5_m100_entry`. It is an all-shift contiguous theorem at
one heat parameter, not yet the forward or arbitrary-column theorem.

Machine-audited companion:

```text
outputs/jensen_window_pf_compound_order5_m100_entry_certificate.md
```

### Theorem 11.45: Fixed Order Five On The Full Heat Interval

Theorem 11.44 supplies positive endpoint data for every shift. The uniform
eventual tail (11.39.4) and the cooperative system (11.39.8) confine a first
loss of positivity to finitely many coordinates. For the integrating factor
`E_n(lambda)>0`, variation of constants reads

```text
Q_n(lambda)=E_n(lambda)
 [Q_n(-100)+integral_(-100)^lambda
  E_n(s)^(-1)a_n(s)Q_(n+1)(s) ds].                 (11.45.1)
```

Backward induction from the positive uniform tail makes the integrand in
(11.45.1) nonnegative. Therefore

```text
H_(5,n)(lambda)>0
for every n>=0 and every -100<=lambda<=0.           (11.45.2)
```

Apply the fixed-order Gasca-Pena initial-minor transfer of Theorem 11.38
pointwise in `lambda`. The completed contiguous layers through order five
then give

```text
R_(5,n)(j_1,j_2,j_3,j_4,j_5;lambda)>0             (11.45.3)
```

for every `n>=0`, every `0<=j_1<j_2<j_3<j_4<j_5`, and every
`-100<=lambda<=0`.

Thus the consecutive-row arbitrary-column signed-Hankel structure is closed
through fixed order five on the complete heat interval. No claim is made for
order six or higher, PF-infinity, the all-degree Jensen bridge, RH, or
`Lambda<=0`. The next fixed-order problem is a cancellation-preserving
order-six entry coordinate and flow reduction.

Machine-audited companion:

```text
outputs/jensen_window_pf_compound_order5_uniform_heat_forward_invariance_certificate.md
```

### Theorem 11.46: Uniform Signed Order-Six Tail And Endpoint Reduction

Put

```text
epsilon_m=(-1)^binom(m,2),
Q_(m,n)=epsilon_m H_(m,n).                          (11.46.1)
```

Thus `Q_(6,n)=-H_(6,n)`. The all-order suitable-multiplier theorem, the
compact heat-tilt audit through order sixteen, and Lemma 11.35 give the
graded order-six ratio contract uniformly for `-100<=lambda<=0`.

For the normalized six-by-six determinant, exact algebra over all 720
permutations and 684 weighted monomials gives

```text
[h^0,...,h^15]det K
 =[0,...,0,-1132462080 G_2^15].                    (11.46.2)
```

Every contribution containing `G_3,...,G_16` cancels from the first nonzero
coefficient. The orientation in (11.46.1) changes its sign, so there is one
finite, possibly ineffective, integer `N_6` such that

```text
Q_(6,n)(lambda)=-H_(6,n)(lambda)>0
for n>=N_6 and -100<=lambda<=0.                    (11.46.3)
```

Desnanot-Jacobi gives the signed recursion

```text
Q_(m,n)Q_(m-2,n+2)
 =Q_(m-1,n+1)^2-Q_(m-1,n)Q_(m-1,n+2).             (11.46.4)
```

At order six this is

```text
Q_(6,n)H_(4,n+2)
 =H_(5,n+1)^2-H_(5,n)H_(5,n+2).                   (11.46.5)
```

Since orders four and five are complete, signed order-six endpoint entry is
equivalent to strict log concavity of the positive `H_5` sequence.

The affine heat derivative and adjacent Plucker identity similarly give

```text
Q_(6,n)'=a_n Q_(6,n+1)+b_n Q_(6,n),
a_n=(4n+42)H_(5,n)/H_(5,n+1)>0,
b_n=((4n+42)/(4n+38))(log H_(5,n+1))'.             (11.46.6)
```

Thus the order-six flow is cooperative over the completed order-five cone.
The uniform tail (11.46.3) and variation of constants imply

```text
[Q_(6,n)(-100)>0 for every n]
 =>[Q_(6,n)(lambda)>0 for every n
    and every -100<=lambda<=0].                    (11.46.7)
```

This endpoint condition is genuinely new. The exact rational countermodel in
the companion certificate has every available strict signed contiguous minor
through order five but `Q_(6,0)<0`. Therefore the lower cone alone cannot
promote order six.

Machine-audited companion:

```text
outputs/jensen_window_pf_compound_order6_uniform_tail_flow_reduction.md
```

### Corollary 11.47: Rigorous Signed Order-Six Prefix

Inside the completed order-five cone define

```text
K_n=H_(5,n+1)^2/[H_(5,n)H_(5,n+2)]-1.             (11.47.1)
```

Equation (11.46.5) gives the exact positive-scale identity

```text
Q_(6,n)=H_(5,n)H_(5,n+2)K_n/H_(4,n+2).            (11.47.2)
```

The merged coefficient sources provide 1024-bit outward-rounded balls for
`A_0(-100),...,A_326(-100)`. Direct stable evaluation proves

```text
K_n(-100)>0,
Q_(6,n)(-100)>0,
for every 0<=n<=316.                               (11.47.3)
```

The weakest relative lower bound occurs at `n=316` and exceeds

```text
0.0078096078288102208346.                          (11.47.4)
```

Consequently the endpoint target in Theorem 11.46 is reduced to the analytic
tail `n>=317`.

Machine-audited companion:

```text
outputs/jensen_window_pf_compound_order6_m100_prefix_certificate.md
```

### Lemma 11.48: Canonical Stable Curvature For The Order-Six Tail

Retain the continuous stable coordinates of Lemma 11.42:

```text
x=exp(-B),  d=1-x,
g=d^2-x^2d(t-1)d(t+1),  h=log g,
f=g^2-x^3g(t-1)g(t+1),  q=log(f/d).                (11.48.1)
```

Add the third gap and its stable logarithm,

```text
S(t)=4B(t)-q(t-1)+2q(t)-q(t+1),
p(t)=2q(t)-h(t)+log(1-exp(-S(t))).                 (11.48.2)
```

Exact elimination of all coefficient-ratio prefactors gives the canonical
identity

```text
H_(5,n)=A_(n+4)^5 exp(p(n+4)).                     (11.48.3)
```

Set `k=n+5` and

```text
P_k=p(k-1)-2p(k)+p(k+1).                           (11.48.4)
```

Then

```text
D_n=log[H_(5,n)H_(5,n+2)/H_(5,n+1)^2]
   =5log x_k+P_k,                                  (11.48.5)

Q_(6,n)>0 iff D_n<0.                               (11.48.6)
```

The completed endpoint defect theorem gives, for `k>=320`,

```text
-5log x_k>=5d_k>=251/[50(2k+1)].                  (11.48.7)
```

For every `k>=322`, exact rational arithmetic proves

```text
320/k^2<251/[50(2k+1)].                            (11.48.8)
```

After setting `k=322+m`, the cleared numerator is

```text
251m^2+129644m+15704684,                           (11.48.9)
```

whose coefficients are positive. Hence

```text
P_k<=320/k^2 for every k>=322                      (11.48.10)
```

is sufficient for `Q_(6,n)(-100)>0` on every `n>=317`.

Machine-audited companion:

```text
outputs/jensen_window_pf_compound_order6_m100_tail_curvature_reduction.md
```

### Lemma 11.49: First-Summand Bridge For The Order-Six Curvature

The endpoint first-summand theorem sharpens the complete-to-first moment
defect to

```text
0<=delta_k=(M_k-M_k^(1))/M_k^(1)<2/k^7,
k>=316.                                            (11.49.1)
```

This is a half-line theorem: exact derivative monotonicity reduces both
endpoint comparisons to `k=316`, where outward-rounded Arb logarithmic gates
are strictly negative.

The third stable gap remains away from zero in both kernels:

```text
min(S_j,S_j^(1))>=3/(2j),  j>=321.                 (11.49.2)
```

For orientation, the exact perturbation chain is

```text
a_j=2[(j-1)^(-7)+2j^(-7)+(j+1)^(-7)],
L_j=4j a_j,
U_j=2a_j+L_(j-1)+2L_j+L_(j+1),
V_j=2L_j+8jU_j,
W_j=3a_j+V_(j-1)+2V_j+V_(j+1),
E_j=2V_j+(5j/7)W_j+L_j,
Z_j=4a_j+E_(j-1)+2E_j+E_(j+1),
Y_j=2E_j+V_j+(2j/3)Z_j.                            (11.49.3)
```

Triangle and mean-value inequalities inside the proved floors give
`|p_j-p_j^(1)|<=Y_j`. Clearing denominators and setting `k=322+m` produces a
degree-75 polynomial with all 76 coefficients positive, proving

```text
|P_k-P_k^(1)|<100/k^2,  k>=322.                   (11.49.4)
```

At the splice the exact scaled transfer is

```text
322^2 transfer_error
 =97.4494017121864439638...<100.                   (11.49.5)
```

The tent identity gives

```text
P_k^(1)=integral_[-1,1](1-|s|)p_1''(k+s) ds.       (11.49.6)
```

Therefore

```text
p_1''(t)<=200/t^2 for every t>=321
 =>P_k^(1)<201/k^2
 =>P_k<301/k^2<320/k^2,  k>=322.                  (11.49.7)
```

Machine-audited companions:

```text
outputs/jensen_window_pf_negative_lambda_first_summand_power7_dominance_extension.md
outputs/jensen_window_pf_compound_order6_first_summand_curvature_bridge.md
```

### Theorem 11.50: Global Third-Nested First-Summand Curvature

The continuous hypothesis in (11.49.7) holds. Four centered stable layers
require a common `t+-4` cover and derivatives `H^(2),...,H^(10)`.

First, the hashed compact H cache, an aligned 107452-row H9/H10 extension,
and one complete right-collar tile assemble into 38 adaptive blocks. Every
block proves positive `J`, `R`, and `S` coordinates and

```text
p_1''(t)<=200/t^2 for 321<=t<=V'(2).               (11.50.1)
```

The largest compact scaled upper is

```text
t^2p_1''(t)<197.8700.                              (11.50.2)
```

Second, the exact signed cumulant corridors through order eight are extended
by the coarse but rigorous bounds

```text
|kappa_r|q^(r/2-1)/(r-2)!<50000,
r=9,10, u>=2.                                      (11.50.3)
```

The epsilon-ten partition recurrence has 45 terms at order nine and 64 at
order ten. Rational termwise majorization and the existing unit-disk Cauchy
residual prove (11.50.3).

A rigorous mode-two quadrature collar and 17999 rational mode blocks then use
dimensionless stable arithmetic to prove

```text
p_1''(t)<=200/t^2 for every saddle mode 2<=u<=20.  (11.50.4)
```

The largest finite-ray scaled upper is below `21.337`.

Third, on the enlarged `u>=19` collar, the normalized signed H boxes through
order eight and (11.50.3) imply

```text
b_0 in [969/1000,1001/1000],
b_1 in [-1001/1000,-959/1000],
(-1)^r b_r in [0,1001/1000], 2<=r<=6,
|b_7|<2500, |b_8|<2500.                            (11.50.5)
```

Put `z=1/t`. At every stable layer evaluate

```text
log(1-exp(-zv))=log z+log v+R_0(zv),
R_0(w)=log((1-exp(-w))/w).                         (11.50.6)
```

The convergent defect series and exact partial-Bell identity bound all
analytic correction coefficients uniformly on `0<=z<=10^(-30)`. One
outward-rounded dimensionless interval evaluation of the complete four-layer
recursion gives

```text
t^2p_1''(t)
 <22.7683610696345567704...<100<200
for every mode u>=20.                              (11.50.7)
```

The ranges (11.50.1), (11.50.4), and (11.50.7) cover every real `t>=321`.
Consequently

```text
p_1''(t)<=200/t^2 for every real t>=321.            (11.50.8)
```

Machine-audited companions:

```text
outputs/jensen_window_pf_compound_order6_high_cumulant_coarse_corridor.md
outputs/jensen_window_pf_compound_order6_nested_curvature_compact_certificate.md
outputs/jensen_window_pf_compound_order6_nested_curvature_finite_ray_certificate.md
outputs/jensen_window_pf_compound_order6_nested_curvature_asymptotic_ray_certificate.md
```

### Theorem 11.51: All-Shift Signed Order-Six Entry At Lambda Minus One Hundred

Theorem 11.50 and (11.49.6) give

```text
P_k^(1)<201/k^2,  k>=322.                          (11.51.1)
```

Combining (11.51.1) with (11.49.4) yields

```text
P_k<301/k^2<320/k^2,  k>=322.                     (11.51.2)
```

Equations (11.48.5)--(11.48.8) therefore prove

```text
Q_(6,n)(-100)=-H_(6,n)(-100)>0
for every n>=317.                                  (11.51.3)
```

The rigorous prefix (11.47.3) covers the complementary shifts. Hence

```text
Q_(6,n)(-100)=-H_(6,n)(-100)>0
for every integer n>=0.                            (11.51.4)
```

This discharges `target_compound_order6_m100_entry`. It is an all-shift
contiguous theorem at one heat parameter, not by itself the forward or
arbitrary-column theorem.

Machine-audited companion:

```text
outputs/jensen_window_pf_compound_order6_m100_entry_certificate.md
```

### Theorem 11.52: Fixed Order Six On The Full Heat Interval

Theorem 11.51 supplies positive signed endpoint data for every shift. The
uniform eventual tail (11.46.3) and cooperative system (11.46.6) confine a
first loss to finitely many coordinates. Variation of constants and backward
induction from the positive tail prove

```text
Q_(6,n)(lambda)=-H_(6,n)(lambda)>0
for every n>=0 and every -100<=lambda<=0.           (11.52.1)
```

Together with Theorem 11.45, all signed contiguous layers through order six
are now complete on the heat interval. Apply the fixed-order Gasca-Pena
initial-minor transfer of Theorem 11.38 pointwise in `lambda`. It gives

```text
epsilon_6 R_(6,n)(j_1,j_2,j_3,j_4,j_5,j_6;lambda)>0
                                                               (11.52.2)
```

for every `n>=0`, every increasing column sextuple, and every
`-100<=lambda<=0`.

Thus the consecutive-row arbitrary-column signed-Hankel structure is closed
through fixed order six on the complete heat interval. No claim is made for
order seven or higher, PF-infinity, the all-degree Jensen bridge, RH, or
`Lambda<=0`. The next fixed-order problem is a cancellation-preserving
order-seven entry coordinate and its corresponding analytic tail.

Machine-audited companion:

```text
outputs/jensen_window_pf_compound_order6_uniform_heat_forward_invariance_certificate.md
```

### Theorem 11.53: Graded-Kernel Vandermonde At Every Fixed Order

Fix `m>=1`, put

```text
D=binom(m,2),  L=2(m-1),  epsilon_m=(-1)^D,
K_(i,j)=exp(-sum_(r>=2)G_r h^(r-1)(L-i-j)^r).     (11.53.1)
```

Writing `z_i=L-i` and `y_j=j`, row and column gauges with constant term one
reduce the determinant to the mixed kernel

```text
M_h(z,y)=exp(sum_(a,b>=1)c_(a,b)h^(a+b-1)z^a y^b),
c_(1,1)=2G_2.                                     (11.53.2)
```

If

```text
M_h(z,y)=sum_(p,q>=0)c_(p,q)(h)z^p y^q,           (11.53.3)
```

then a product of `ell` mixed factors with bidegree `(p,q)` has exact
`h`-degree `p+q-ell`. Since `ell<=min(p,q)`,

```text
ord_h c_(p,q)>=max(p,q),
[h^k]c_(k,k)=(2G_2)^k/k!.                         (11.53.4)
```

In the Cauchy-Binet expansion, increasing exponent tuples `P,Q` satisfy

```text
sum_i max(p_i,q_(pi(i)))>=D.                      (11.53.5)
```

Equality forces `P=Q=(0,1,...,m-1)`. For this block,

```text
sum_i max(i,pi(i))=D+(1/2)sum_i|i-pi(i)|,         (11.53.6)
```

so only the identity coefficient permutation contributes. The two remaining
alternants are Vandermonde determinants, and therefore

```text
[h^D]det K
 =epsilon_m 2^D prod_(j=1)^(m-1)j! G_2^D.        (11.53.7)
```

Consequently the signed minor `Q_(m,n)=epsilon_m H_(m,n)` has a positive
universal first term at every fixed order.

The heat-tilt quantifier is also available at arbitrary fixed difference
order. O'Sullivan's suitable-multiplier theorem may be truncated at every
fixed order, while Cauchy estimates make the compact Gaussian-log family
uniform for `0<=T<=100`. If `w=W(2k/pi)`, define

```text
F_0(w)=w^2,
F_(s+1)(w)=-sF_s(w)+wF_s'(w)/(1+w).               (11.53.8)
```

Then

```text
d^s(w^2)/dk^s=k^(-s)F_s(w),  F_s(w)=O_s(w).       (11.53.9)
```

For fixed `s`, choose a saddle truncation `R>s`. The finite difference of the
uniform remainder is bounded directly by

```text
2^s O_R(w^(3R)/k^R)=o(w/k^s),                    (11.53.10)
```

without differentiating an unspecified remainder. Hence

```text
Delta^s log R_T^(1)(k)=O_s(log(k)/k^s)            (11.53.11)
```

uniformly for `0<=T<=100`. Composing this with the all-order Xi ratio input
and superpolynomial higher-theta suppression proves the correctly ordered
quantifiers

```text
for every fixed m there exists N_m such that
Q_(m,n)(lambda)>0 for n>=N_m and -100<=lambda<=0. (11.53.12)
```

The threshold may depend on `m` and is not effective here. Equation
(11.53.12) is neither a common-tail theorem uniform in `m` nor an all-shift
or PF-infinity theorem.

Machine-audited companion:

```text
outputs/jensen_window_pf_graded_kernel_vandermonde_all_order_lemma.md
```

### Lemma 11.54: Order-Seven Coordinate And Cooperative Flow

At `m=7`, `D=21` and `epsilon_7=-1`. Theorem 11.53 gives

```text
det K=-52183852646400 G_2^21 h^21+o(h^21),
Q_7=-H_7=positive_scale*
       (52183852646400 G_2^21 h^21+o(h^21)).       (11.54.1)
```

Thus there is a compact-uniform eventual positive `Q_7` tail. Signed
Desnanot-Jacobi gives the exact endpoint coordinate

```text
Q_(7,n)H_(5,n+2)
 =Q_(6,n+1)^2-Q_(6,n)Q_(6,n+2).                  (11.54.2)
```

Since the order-five and signed order-six layers are complete and positive,
order-seven entry is exactly strict log-concavity of the endpoint `Q_6`
sequence.

The affine-Hankel and adjacent Plucker identities give

```text
Q_(7,n)'=(4n+50)delta(Q_(7,n)),                   (11.54.3)
Q_n'=a_n Q_(n+1)+b_n Q_n,                        (11.54.4)
a_n=(4n+50)Q_(6,n)/Q_(6,n+1)>0,                  (11.54.5)
b_n=((4n+50)/(4n+46))(log Q_(6,n+1))'.           (11.54.6)
```

Theorem 11.53 confines a first loss to finitely many coordinates. Variation
of constants and backward induction therefore prove the conditional theorem

```text
[Q_(7,n)(-100)>0 for every n>=0]
 =>[Q_(7,n)(lambda)>0 for every n>=0
    and every -100<=lambda<=0].                   (11.54.7)
```

The fixed-order initial-minor transfer would then supply every
arbitrary-column signed order-seven minor.

This endpoint input cannot be inferred from the completed lower cone. The
exact sequence

```text
1/0!,1/1!,...,1/11!,1/480100000                  (11.54.8)
```

has every available signed contiguous layer through order six positive, but

```text
Q6_1^2-Q6_0 Q6_2
 =-34201/59436732624881410374275156859658650098073600000000000000000,
Q7_0=-34201/88845982207870628726518579200000000000. (11.54.9)
```

Hence a new zeta-specific endpoint theorem is indispensable.

Machine-audited companion:

```text
outputs/jensen_window_pf_compound_order7_uniform_tail_flow_reduction.md
```

### Theorem 11.55: Rigorous Order-Seven Endpoint Prefix

Inside the completed lower cone, define

```text
L_n=Q_(6,n+1)^2/[Q_(6,n)Q_(6,n+2)]-1.            (11.55.1)
```

Equation (11.54.2) becomes

```text
Q_(7,n)=Q_(6,n)Q_(6,n+2)L_n/H_(5,n+2),           (11.55.2)
```

so `Q_(7,n)>0` is exactly `L_n>0`. The computation evaluates the completed
stable `H_4`, `H_5`, and `Q_6` factorizations rather than subtracting a raw
seven-by-seven determinant.

Outward-rounded 2048-bit Arb coefficient balls cover

```text
A_k(-100)>0 for 0<=k<=326.                        (11.55.3)
```

A dedicated retained-integral repair covers `A_179,...,A_190` at 220 decimal
digits with the established analytic summand-tail and cutoff-tail bounds.
Every formerly inconclusive row is thereby replaced by an enclosure theorem.
Stable interval evaluation proves

```text
L_n(-100)>0 and Q_(7,n)(-100)>0
for every 0<=n<=314.                              (11.55.4)
```

The weakest row is

```text
L_314 lower
 =0.009383261154960606031722710481968...>9/1000. (11.55.5)
```

This is a rigorous finite prefix. It leaves precisely the endpoint tail
`n>=315`; the non-effective threshold in (11.53.12) cannot by itself splice
the two ranges.

Machine-audited companion:

```text
outputs/jensen_window_pf_compound_order7_m100_prefix_certificate.md
```

### Lemma 11.56: Fourth-Nested Curvature Reduction For Order Seven

At `lambda=-100`, retain the completed stable hierarchy

```text
x=exp(-B), d=1-x,
g=d^2-x^2 d(t-1)d(t+1), h=log g,
f=g^2-x^3 g(t-1)g(t+1), q=log(f/d),
S=4B-q(t-1)+2q-q(t+1),
p=2q-h+log(1-exp(-S)).                             (11.56.1)
```

Signed condensation adds one stable gap and logarithm:

```text
T(t)=5B(t)-p(t-1)+2p(t)-p(t+1),                  (11.56.2)
r(t)=2p(t)-q(t)+log(1-exp(-T(t))),                (11.56.3)
Q_(6,n)=A_(n+5)^6 exp(r(n+5)).                    (11.56.4)
```

For `k=n+6`, set

```text
R_k=r(k-1)-2r(k)+r(k+1),
E_n=log[Q_(6,n)Q_(6,n+2)/Q_(6,n+1)^2]
   =6log x_k+R_k.                                 (11.56.5)
```

Then

```text
L_n=exp(-E_n)-1,
Q_(7,n)>0 iff E_n<0.                              (11.56.6)
```

The completed endpoint defect theorem gives

```text
-6log x_k>=6d_k>=753/[125(2k+1)],  k>=320.        (11.56.7)
```

For every `k>=321`, exact rational arithmetic proves

```text
900/k^2<753/[125(2k+1)].                          (11.56.8)
```

After `k=321+m`, the cleared numerator is

```text
753m^2+258426m+5252373,                           (11.56.9)
```

which is positive coefficientwise. Therefore the single estimate

```text
R_k<=900/k^2 for every k>=321                     (11.56.10)
```

would prove `Q_(7,n)(-100)>0` for every `n>=315` and, together with Theorem
11.55, all-shift endpoint entry. Theorems 11.57--11.58 below close the
complete-to-first-summand transfer, and Theorems 11.59--11.60 close the
continuous fourth-nested curvature theorem through saddle mode `u=2`. The
remaining analytic task is the finite and asymptotic mode ray `u>=2`.

Equation (11.56.10) is still open. Consequently Theorem 11.54 cannot yet be
promoted to all-shift order seven, PF-infinity, the all-degree Jensen bridge,
RH, or `Lambda<=0`.

Machine-audited companion:

```text
outputs/jensen_window_pf_compound_order7_m100_tail_curvature_reduction.md
```

### Theorem 11.57: Rebalanced Inverse-Eighth-Power Dominance

The completed first-summand theorem used the split `a(k)=log(k)/8`. Its
high-region estimate had far more reserve than the low-region probability.
Rebalance the same exact summand-ratio argument by setting

```text
a(k)=log(k)/10,
b(k)=a(k)+1/10,
c(k)=b(k)+1/100,
q(a(k))=pi k^(2/5).                                (11.57.1)
```

The inherited theorem `epsilon(0)<0.0022` and the same tilted-integral
comparison give the two candidate bounds

```text
epsilon(a(k))<=17exp(-3pi k^(2/5)),                (11.57.2)
epsilon(0)P_k(u<a(k))<k^(-8).                     (11.57.3)
```

At `k=300`, outward-rounded Arb arithmetic proves

```text
log 17-3pi k^(2/5)+8log k<-43.8190,               (11.57.4)
B_10(k)+log(0.0022)+8log k<-22.4542,              (11.57.5)
17-high_prefactor>0.17979.                         (11.57.6)
```

All fourteen endpoint and derivative gates are strict. The six ratio-log
derivatives have the forms

```text
3/5+1/(L-alpha)-2/(L+beta)>0,
1+1/(L-alpha)-2/(L+beta)>0,                       (11.57.7)
```

for the three relevant pairs `(alpha,beta)`, so every comparison strengthens
on `k>=300`. Consequently

```text
epsilon(a(k))<k^(-8),
epsilon(0)P_k(u<a(k))<k^(-8),                     (11.57.8)
0<=delta_k=(M_k-M_k^(1))/M_k^(1)<2/k^8           (11.57.9)
```

for every integer `k>=300`. Since `log(1+z)<=z`, the adjacent wall satisfies

```text
|B_j-B_j^(1)|<=a_j
 =2[(j-1)^(-8)+2j^(-8)+(j+1)^(-8)],  j>=301.     (11.57.10)
```

This is a lambda-minus-one-hundred moment-tail theorem. It does not itself
prove any fourth-nested curvature inequality.

Machine-audited companion:

```text
outputs/jensen_window_pf_negative_lambda_first_summand_power8_rebalanced_dominance_extension.md
```

### Lemma 11.58: First-Summand Bridge For The Order-Seven Curvature

Write the fourth stable layer as

```text
T(t)=5B(t)-p(t-1)+2p(t)-p(t+1),
r(t)=2p(t)-q(t)+log(1-exp(-T(t))).                (11.58.1)
```

The completed order-six first/full bounds imply, for `j>=322`,

```text
T_j^(1)>=5/(2j+1)-201/j^2,
T_j>=251/[50(2j+1)]-301/j^2.                     (11.58.2)
```

After `j=322+m`, subtracting `3/(2j)` leaves respectively

```text
4m^2+1769m+154480,
101m^2+34869m+740684,                             (11.58.3)
```

up to positive denominators. For `j=320,321`, the rigorous order-six prefix
gives `T_j=log(1+K_(j-5))`; using `log(1+x)>=x/(1+x)` and (11.57.10) leaves
first/full margins above `0.00307`. Hence

```text
min(T_j,T_j^(1))>=3/(2j),  j>=320.                (11.58.4)
```

The complete perturbation chain is

```text
a_j=2[(j-1)^(-8)+2j^(-8)+(j+1)^(-8)],
L_j=4j a_j,
U_j=2a_j+L_(j-1)+2L_j+L_(j+1),
V_j=2L_j+8jU_j,
W_j=3a_j+V_(j-1)+2V_j+V_(j+1),
E_j=2V_j+(5j/7)W_j+L_j,
Z_j=4a_j+E_(j-1)+2E_j+E_(j+1),
Y_j=2E_j+V_j+(2j/3)Z_j,
O_j=5a_j+Y_(j-1)+2Y_j+Y_(j+1),
N_j=2Y_j+E_j+(2j/3)O_j.                          (11.58.5)
```

Triangle and mean-value inequalities inside (11.58.4) give
`|T_j-T_j^(1)|<=O_j` and `|r_j-r_j^(1)|<=N_j`. Exact rational
common-denominator arithmetic then proves

```text
|R_k-R_k^(1)|
 <=N_(k-1)+2N_k+N_(k+1)<262/k^2,  k>=321.       (11.58.6)
```

After `k=321+m`, the shifted numerator has degree 102 and all 103
coefficients are positive. At the splice,

```text
321^2 transfer_error
 =261.334434117711045553...<262.                  (11.58.7)
```

The tent identity now reduces the full endpoint theorem to one continuous
first-summand estimate:

```text
r_1''(t)<=600/t^2 for every real t>=320           (11.58.8)
 =>R_k^(1)<601/k^2
 =>R_k<863/k^2<900/k^2,  k>=321.                 (11.58.9)
```

Thus (11.58.8) closes every endpoint shift `n>=315`, splices to Theorem
11.55, and activates the conditional flow theorem (11.54.7), once its four
continuous pieces are composed. Sections 11.59--11.63 prove those pieces;
Sections 11.64--11.65 record the endpoint and heat-flow compositions. The
all-order PF-infinity and Jensen bridges remain separate.

Machine-audited companions:

```text
outputs/jensen_window_pf_negative_lambda_first_summand_power8_rebalanced_dominance_extension.md
outputs/jensen_window_pf_compound_order7_first_summand_curvature_bridge.md
```

## 11.59 Order-seven shifted-jet compact bridge

The compact part of (11.58.8) is now a rigorous interval theorem. Fix a
rational anchor `a` and write `B=Delta^2 H` for the first stable gap. For any
gap coordinate `C` and logarithmic layer `f`, use the normalized Taylor
coefficients

```text
c_n=a^(n+1) C^(n)/n!,
f_n=a^n f^(n)/n!.                                  (11.59.1)
```

On the common collar `[a-5,b+5]`, the unit-tent identity
`Delta^2 f(t)=integral_(-1)^1(1-|v|)f''(t+v)dv` first gives

```text
b_n in a^(n+1) Hull(H^(n+2))/n!.
```

If `f` is the current stable logarithmic layer and `m=2,3,4,5` is the
leading factor at the next gap, exact coefficient scaling gives

```text
c_n=m*b_n-((n+1)(n+2)/a)*f_(n+2).                 (11.59.2)
```

Outward-rounded stable-log series propagate (11.59.2) through `J,R,S,T`.
Each centered difference consumes one unit of collar, so `H^(2),...,H^(21)`
on `[a-5,b+5]` enclose every normalized coefficient through order eleven on
the central block `[a,b]`. Strict positivity of the zeroth `B,J,R,S,T`
intervals certifies the domain of every logarithm.

At the exact anchor, degree-thirty exact-potential panel Taylor quadrature
encloses `H^(0),...,H^(10)` at all eleven shifts `a+j`, `-5<=j<=5`. Exact
centered-difference series algebra gives the point jet
`r_n(a)=r^(n)(a)/n!` through order ten. If

```text
rho_11=a^11*r^(11)/11!,
```

then for `0<=h<=b-a`, Taylor's theorem and the common collar give

```text
r''(a+h)
 <=2r_2(a)+sum_(n=3)^10 n(n-1)|r_n(a)|h^(n-2)
   +(110|rho_11|/a^2)(h/a)^9.                     (11.59.3)
```

The rational partition starts at `a_0=320` and uses

```text
a_(i+1)=min(1000,(1/10)floor(10*(161/160)*a_i)),  (11.59.4)
```

so `(a_(i+1)-a_i)/a_i<=1/160`. A hashed `1/100` t-cache supplies the first
73 collars and an overlapping hashed `1/10` cache supplies the remaining 113.
All 186 blocks pass. The weakest normalized final-gap floor and the largest
scaled curvature upper are

```text
min T_dimensionless >2.0695966777590944,
max_i a_(i+1)^2*sup_(block i) r''
 <50.910519421521962<600,                          (11.59.5)
```

with the largest upper on `500<=t<=503.1`. Therefore

```text
r_1''(t)<=600/t^2 for every real 320<=t<=1000.     (11.59.6)
```

This discharges the first compact part of (11.58.8). The next theorem extends
the result from `t=1000` through saddle mode `u=2`, and Sections 11.62--11.65
close the remaining fixed-order-seven pieces. This theorem alone supplies no
PF-infinity, all-degree Jensen bridge, RH, or `Lambda<=0` conclusion.

Machine-audited companion:

```text
outputs/jensen_window_pf_compound_order7_shifted_jet_t320_t1000_certificate.md
```

## 11.60 Order-seven nested-curvature compact certificate

The exact-mode compact method continues (11.59.6) without a Taylor remainder
in `t`. On each rational mode tile, rigorous interval quadrature encloses

```text
H^(2),...,H^(12)
```

on one common `t+-5` collar. The aligned source grid has mode width `10^(-5)`
and reaches beyond `u=2`; one additional right-collar tile reaches
`u=2.00003`. The H2--H8, H9--H10, and H11--H12 files contain exactly 107452
aligned rows apiece.

For a local series `F(z)=sum_n F_n z^n`, put

```text
B_n=H^(n+2)/n!,
d=1-exp(-B), ell=log(d).                           (11.60.1)
```

The four cancellation-preserving stable layers are propagated by

```text
J_n=2B_n-(n+1)(n+2)ell_(n+2),
h=2ell+log(1-exp(-J)),
R_n=3B_n-(n+1)(n+2)h_(n+2),
q=2h-ell+log(1-exp(-R)),
S_n=4B_n-(n+1)(n+2)q_(n+2),
p=2q-h+log(1-exp(-S)),
T_n=5B_n-(n+1)(n+2)p_(n+2),
r=2p-q+log(1-exp(-T)).                            (11.60.2)
```

Every operation in (11.60.1)--(11.60.2) is outward rounded on the same
derivative hull. Strict lower bounds for `B,d,J,R,S,T` certify the logarithm
domains. Adaptive bisection gives 82 contiguous rational mode blocks from the
tile containing `t=1000` through `u=2`; every central block has more than five
units of cached `t` collar on both sides. The extremal diagnostics are

```text
max t^2*sup r_1''(t)
 <358.732317558974206057490510776<600,
min T>0.0001024760187815426029870334694.          (11.60.3)
```

The first block starts below `t=1000`, and the shifted-jet theorem (11.59.6)
contains the splice point. The final block ends at mode two, where

```text
V'(2)=37850.32221021138481629511089256....
```

Consequently

```text
r_1''(t)<=600/t^2 for every real 1000<=t<=V'(2), (11.60.4)
r_1''(t)<=600/t^2 for every real 320<=t<=V'(2).  (11.60.5)
```

The finite and asymptotic rays are proved separately in Sections 11.62 and
11.63. This compact theorem alone does not prove endpoint entry, PF-infinity,
the all-degree Jensen bridge, RH, or `Lambda<=0`.

Machine-audited companion:

```text
outputs/jensen_window_pf_compound_order7_nested_curvature_compact_certificate.md
```

## 11.61 Order-seven high-cumulant corridors

The finite and asymptotic saddle arguments require `H` derivatives through
order twelve. For `r=11,12`, apply the exact epsilon-ten tilted-Gaussian
partition recurrence to the already certified potential box
`L_3,...,L_12`. No potential derivative of order thirteen or fourteen enters.
Termwise rational substitution in all 72 formal terms gives

```text
|scaled kappa_r^[10]|<14000,   2<=u<=20,
|scaled kappa_r^[10]|<700000,  u>=20.              (11.61.1)
```

The exact Cauchy factor is `r(r-1)<=12*11=132`. In both ranges the certified
unit-disk residual, multiplied by 132, is below one. Consequently

```text
|kappa_r| q^(r/2-1)/(r-2)!<14001,
    r=11,12, 2<=u<=20,
|kappa_r| q^(r/2-1)/(r-2)!<700001,
    r=11,12, u>=20.                                (11.61.2)
```

This is an exact coarse corridor theorem, not a sampled asymptotic fit.

Machine-audited companion:

```text
outputs/jensen_window_pf_compound_order7_high_cumulant_coarse_corridor.md
```

## 11.62 Order-seven finite saddle ray

Combine the exact normalized cumulant corridors through order eight, the
order-nine and order-ten cap `50000`, and (11.61.2). They enclose
`H^(2),...,H^(12)` on the strict `t+-5` collar required by (11.60.2). A short
aligned quadrature extension covers the first mode block
`2<=u<=2001/1000`; the rest of the interval is partitioned into rational
blocks of width `1/1000`.

On every block, outward-rounded dimensionless series arithmetic proves
strict positivity of `J,R,S,T` and encloses the fourth-layer curvature. All
17,999 ray blocks and the initial collar pass. The extremal diagnostics are

```text
max t^2*sup r_1''(t)
 <73.5426280196174861068325277832<600,
min T_dimensionless
 >3.854837399921988165453005686.                   (11.62.1)
```

Therefore

```text
r_1''(t)<=600/t^2 for every saddle mode 2<=u<=20.  (11.62.2)
```

The rational cover and the exact initial collar overlap, so (11.62.2) is a
continuous interval theorem and not point sampling.

Machine-audited companion:

```text
outputs/jensen_window_pf_compound_order7_nested_curvature_finite_ray_certificate.md
```

## 11.63 Order-seven asymptotic saddle ray

Put `z=1/t`. The exact low-cumulant ray theorem and the high-cumulant
corridors imply normalized `H` boxes through order twelve; in particular,

```text
|b_7|,|b_8|<2500,
|b_9|,|b_10|<40000.                                (11.63.1)
```

At each stable layer evaluate the small argument without cancellation by

```text
log(1-exp(-zv))=log(z)+log(v)+R(zv),
R(w)=log((1-exp(-w))/w).                           (11.63.2)
```

The convergent defect series and the exact partial-Bell identity bound the
coefficient error by `C*10^(-30)` in degree zero and by
`2^(r-1)C*10^(-30)` in degree `r>=1`. One outward-rounded evaluation over the
whole interval `0<=z<=10^(-30)` gives

```text
J_0>1.9379999996,
R_0>2.9069999995,
S_0>3.8759999992,
T_0>4.8449999991,
t^2 r_1''(t)<55.540195343064144254<100<600.        (11.63.3)
```

Thus

```text
r_1''(t)<=600/t^2 for every saddle mode u>=20.     (11.63.4)
```

Equations (11.59.6), (11.60.4), (11.62.2), and (11.63.4), including their
shared endpoints, cover every real `t>=320`. Hence the formerly conditional
first-summand target is now a theorem:

```text
r_1''(t)<=600/t^2 for every real t>=320.           (11.63.5)
```

Machine-audited companion:

```text
outputs/jensen_window_pf_compound_order7_nested_curvature_asymptotic_ray_certificate.md
```

## 11.64 All-shift order-seven entry at lambda=-100

Insert (11.63.5) into the exact tent identity of Lemma 11.58. For every
integer `k>=321`,

```text
R_k^(1)<600[-log(1-1/k^2)]<601/k^2.                (11.64.1)
```

The complete-to-first transfer (11.58.6) then gives

```text
R_k<R_k^(1)+|R_k-R_k^(1)|
   <601/k^2+262/k^2
   =863/k^2<900/k^2.                               (11.64.2)
```

The exact tail reduction of Theorem 11.56 converts (11.64.2), with `k=n+6`,
into

```text
Q_(7,n)(-100)>0 for every n>=315.                  (11.64.3)
```

The 1024-bit prefix theorem (Theorem 11.55) supplies the complementary range
`0<=n<=314`. Since `epsilon_7=-1`, prefix and tail prove

```text
Q_(7,n)(-100)=-H_(7,n)(-100)>0
for every integer n>=0.                            (11.64.4)
```

This is the all-shift endpoint theorem at one heat parameter. Its proof uses
the complete Newman kernel in (11.64.2); it is not a first-summand-only
conclusion.

Machine-audited companion:

```text
outputs/jensen_window_pf_compound_order7_m100_entry_certificate.md
```

## 11.65 Fixed order-seven heat-flow closure

The all-fixed-order Vandermonde theorem gives a compact-uniform eventual tail:

```text
there exists N_7 such that Q_(7,n)(lambda)>0
for n>=N_7 and -100<=lambda<=0.                    (11.65.1)
```

Inside the completed order-six cone, signed Desnanot-Jacobi and the affine
Hankel heat derivative give the cooperative system

```text
Q_n'=a_n Q_(n+1)+b_n Q_n,
a_n=(4n+50)Q_(6,n)/Q_(6,n+1)>0,                   (11.65.2)
b_n=((4n+50)/(4n+46))(log Q_(6,n+1))'.
```

The uniform tail (11.65.1) confines a hypothetical first loss to finitely
many shifts. On that finite set, the integrating-factor identity is

```text
Q_n(lambda)=E_n(lambda)[Q_n(-100)
 +integral_(-100)^lambda E_n(s)^(-1)a_n(s)Q_(n+1)(s) ds].
                                                               (11.65.3)
```

Starting from the positive tail and inducting backward, the integral in
(11.65.3) is nonnegative, while (11.64.4) makes its initial term strictly
positive. Therefore no first loss exists and

```text
Q_(7,n)(lambda)=-H_(7,n)(lambda)>0
for every n>=0 and -100<=lambda<=0.                (11.65.4)
```

Together with the completed lower layers, (11.65.4) supplies every signed
initial minor through order seven. Applying the fixed-order Gasca-Pena
initial-minor criterion pointwise in `lambda` proves every consecutive-row,
arbitrary-column reshaped-Hankel minor of order at most seven has its required
sign on the same heat interval.

This closes fixed order seven only. The countermodel in the order-seven flow
reduction proves that lower-layer positivity alone does not generate the next
order. Fixed order eight therefore requires the new endpoint input constructed
in Sections 11.67--11.70; higher orders still require their own input or a
genuinely uniform-in-order theorem. PF-infinity, the all-degree Jensen bridge,
RH, and `Lambda<=0` remain unproved.

Machine-audited companions:

```text
outputs/jensen_window_pf_compound_order7_uniform_tail_flow_reduction.md
outputs/jensen_window_pf_compound_order7_uniform_heat_forward_invariance_certificate.md
```

## 11.66 Fixed order-eight propagation reduction

The all-fixed-order Vandermonde theorem specializes at `m=8` with

```text
D=binom(8,2)=28,
epsilon_8=1,
[h^28]Q_8=33664847019245568000*G_2^28>0.           (11.66.1)
```

Consequently there is a finite, possibly order-dependent and non-effective,
threshold `N_8` such that

```text
Q_(8,n)(lambda)>0 for n>=N_8 and -100<=lambda<=0. (11.66.2)
```

Signed Desnanot-Jacobi over the completed order-six and order-seven cone gives

```text
Q_(8,n)Q_(6,n+2)
 =Q_(7,n+1)^2-Q_(7,n)Q_(7,n+2).                  (11.66.3)
```

Since `Q_(6,n+2)>0`, endpoint order-eight positivity is exactly strict
log-concavity of the endpoint `Q_7` sequence. The affine heat derivative and
adjacent Plucker identity specialize to

```text
Q_(8,n)'=(4n+58)delta(Q_(8,n)),
Q_n'=a_n Q_(n+1)+b_n Q_n,
a_n=(4n+58)Q_(7,n)/Q_(7,n+1)>0,
b_n=((4n+58)/(4n+54))(log Q_(7,n+1))'.            (11.66.4)
```

Thus (11.66.2), finite confinement, and variation of constants prove the
conditional implication

```text
[Q_(8,n)(-100)>0 for every n]
 =>[Q_(8,n)(lambda)>0 for every n and -100<=lambda<=0],
                                                               (11.66.5)
```

after which the fixed-order initial-minor theorem transfers the result to
arbitrary columns.

This reduction genuinely needs new endpoint input. Let

```text
a_j=1/j! for 0<=j<=13,
a_14=1/87120000000.                                (11.66.6)
```

Exact rational determinant evaluation proves every available signed
contiguous minor through order seven is strictly positive, whereas

```text
Q_(7,1)^2-Q_(7,0)Q_(7,2)
 =-463/10246317406459624992735028091299392073078206541205704645017600000000000000000000000,
Q_(8,0)
 =-463/69210459277322870541707704182767616000000000000000<0.  (11.66.7)
```

Hence completed lower layers do not imply order eight. The sole new
fixed-order-eight input is the zeta-specific theorem

```text
Q_(8,n)(-100)>0 for every n>=0,
```

equivalently strict log-concavity of `Q_(7,n)(-100)`. Sections 11.67--11.70
discharge this zeta-specific input. All higher orders, PF-infinity, the
all-degree Jensen bridge, RH, and `Lambda<=0` remain open.

Machine-audited companion:

```text
outputs/jensen_window_pf_compound_order8_uniform_tail_flow_reduction.md
```

## 11.67 Strengthened finite order-eight endpoint prefix

Inside the completed positive order-seven cone, define the stable relative
log-concavity margin

```text
M_n=Q_(7,n+1)^2/(Q_(7,n)Q_(7,n+2))-1.            (11.67.1)
```

Equation (11.66.3) becomes

```text
Q_(8,n)
 =Q_(7,n)Q_(7,n+2)M_n/Q_(6,n+2),                 (11.67.2)
```

so `Q_(8,n)>0` is exactly `M_n>0`. Stable 2048-bit Arb evaluation from
1,257 outward-rounded endpoint coefficient balls, including twelve repaired
coefficients and 930 retained-integral extension rows, proves

```text
M_n(-100)>0 and Q_(8,n)(-100)>0
for every 0<=n<=1242.                             (11.67.3)
```

The weakest certified row is the final shift:

```text
M_1242
 =[0.00355609273479377717793692780206927653148465783678139703539852
   +/- 1.36E-63],
inf M_1242
 =3.55609273479377717793692780206927653148465783678139703539851E-3
 >1/300.                                          (11.67.4)
```

Thus the remaining endpoint problem starts at

```text
Q_(8,n)(-100)>0 for every n>=1243.                (11.67.5)
```

Machine-audited companion:

```text
outputs/jensen_window_pf_compound_order8_m100_prefix_certificate.md
```

## 11.68 Exact order-eight endpoint-tail bridge

Write `k=n+7`. The fifth nested stable logarithm and its centered second
difference are

```text
U(t)=6B(t)-r(t-1)+2r(t)-r(t+1),
s(t)=2r(t)-p(t)+log(1-exp(-U(t))),
W_k=s(k-1)-2s(k)+s(k+1).                          (11.68.1)
```

Exact stable factorization gives

```text
E_n=7log(x_k)+W_k,
M_n=exp(-E_n)-1,
Q_(8,n)(-100)>0 iff E_n<0.                        (11.68.2)
```

The completed coefficient-defect buffer proves

```text
-7log(x_k)>=1757/(250(2k+1)),
4300/k^2<1757/(250(2k+1)), k>=1250.              (11.68.3)
```

Hence `W_k<=4300/k^2` suffices for the complete tail (11.67.5). The
rebalanced retained-summand theorem supplies

```text
0<=delta_k=(M_k-M_k^(1))/M_k^(1)<2/k^9,
k>=300,                                           (11.68.4)
```

and exact five-layer error propagation, positive gap floors, and a
coefficient-positive degree-133 rational inequality give

```text
|W_k-W_k^(1)|<190/k^2, k>=1250.                  (11.68.5)
```

Finally, the tent identity

```text
W_k^(1)=integral_[-1,1](1-|v|)s_1''(k+v)dv       (11.68.6)
```

proves the conditional implication

```text
s_1''(t)<=4000/t^2 for t>=999
 => W_k^(1)<4001/k^2
 => W_k<4191/k^2<4300/k^2,
k>=1250.                                          (11.68.7)
```

Machine-audited companions:

```text
outputs/jensen_window_pf_compound_order8_m100_tail_curvature_reduction.md
outputs/jensen_window_pf_negative_lambda_first_summand_power9_rebalanced_dominance_extension.md
outputs/jensen_window_pf_compound_order8_first_summand_curvature_bridge.md
```

## 11.69 Continuous fifth-nested curvature theorem

The epsilon-ten formal logarithm has vanishing normalized cumulants in
orders thirteen and fourteen. Combining this with the exact-minus-formal
unit-disk residual and the Cauchy factor `14*13=182` proves

```text
|kappa_j|q^(j/2-1)/(j-2)!<1,
j=13,14 and every saddle mode u>=2.               (11.69.1)
```

These caps complete the `H^(2),...,H^(14)` input required by the fifth
stable layer. Four outward-rounded covers then prove

```text
s_1''(t)<=2000/t^2, 699<=t<=999,                 (11.69.2)
s_1''(t)<=4000/t^2, 999<=t<=V'(2),               (11.69.3)
s_1''(t)<=4000/t^2, 2<=u<=20,                    (11.69.4)
t^2 s_1''(t)<134.49<200<4000, u>=20.             (11.69.5)
```

Here (11.69.2) uses 185 common-collar shifted-jet blocks, (11.69.3) uses 96
adaptive compact blocks, and (11.69.4) uses 17,999 contiguous rational mode
blocks. For (11.69.5), exact normalized `H2-H14` boxes and one interval over
`0<=1/t<=10^-30` keep every stable coordinate positive. The mode-to-`t`
cover is exact, so (11.69.3)--(11.69.5) yield

```text
s_1''(t)<=4000/t^2 for every real t>=999.         (11.69.6)
```

Machine-audited companions:

```text
outputs/jensen_window_pf_compound_order8_high_cumulant_coarse_corridor.md
outputs/jensen_window_pf_compound_order8_shifted_jet_t699_t999_certificate.md
outputs/jensen_window_pf_compound_order8_nested_curvature_compact_certificate.md
outputs/jensen_window_pf_compound_order8_nested_curvature_finite_ray_certificate.md
outputs/jensen_window_pf_compound_order8_nested_curvature_asymptotic_ray_certificate.md
```

## 11.70 All-shift order-eight endpoint entry

Substituting (11.69.6) into (11.68.7) proves `W_k<4300/k^2` for every
`k>=1250`. Equations (11.68.2)--(11.68.3) therefore give

```text
Q_(8,n)(-100)=H_(8,n)(-100)>0 for every n>=1243. (11.70.1)
```

Splicing (11.70.1) with the rigorous finite prefix (11.67.3) proves the
zeta-specific endpoint theorem

```text
Q_(8,n)(-100)=H_(8,n)(-100)>0
for every integer n>=0.                           (11.70.2)
```

Machine-audited companion:

```text
outputs/jensen_window_pf_compound_order8_m100_entry_certificate.md
```

## 11.71 Fixed order-eight heat propagation

The compact-uniform eventual theorem (11.66.2), the completed lower cone,
and the exact cooperative recursion (11.66.4) confine any hypothetical first
loss to finitely many shifts. On that finite set, variation of constants is

```text
Q_n(lambda)=E_n(lambda)[Q_n(-100)
 +integral_(-100)^lambda E_n(v)^(-1)a_n(v)Q_(n+1)(v)dv].
                                                               (11.71.1)
```

Backward induction from the positive uniform tail makes the integral
nonnegative, while (11.70.2) supplies the strictly positive initial term.
Consequently

```text
Q_(8,n)(lambda)=H_(8,n)(lambda)>0
for every n>=0 and -100<=lambda<=0.               (11.71.2)
```

Together with the lower layers, (11.71.2) supplies every signed initial
minor through order eight. The fixed-order Gasca-Pena criterion then proves
every consecutive-row, arbitrary-column reshaped-Hankel minor through order
eight has the required sign on the same heat interval.

This closes fixed order eight only. Order nine and higher still require new
zeta-specific endpoint input or a genuinely uniform-in-order mechanism.
PF-infinity, the all-degree Jensen bridge, RH, and `Lambda<=0` remain open.

Machine-audited companion:

```text
outputs/jensen_window_pf_compound_order8_uniform_heat_forward_invariance_certificate.md
```

## 11.72 Signed order-nine coordinate, tail, and finite endpoint prefix

At `m=9`, `D=36` and `epsilon_9=1`. The all-fixed-order graded-kernel
theorem specializes to

```text
det K=347485857744891213250560000*G_2^36*h^36+o(h^36),
there exists N_9 such that Q_(9,n)(lambda)>0
for n>=N_9 and -100<=lambda<=0.                  (11.72.1)
```

Signed Desnanot-Jacobi gives the exact endpoint coordinate

```text
Q_(9,n)Q_(7,n+2)
 =Q_(8,n+1)^2-Q_(8,n)Q_(8,n+2).                 (11.72.2)
```

Since orders seven and eight are complete, `Q_(9,n)>0` is equivalent to
strict log concavity of the positive `Q_8` sequence. The general affine
Hankel and adjacent Plucker identities specialize to

```text
Q_(9,n)'=(4n+66)delta(Q_(9,n)),                  (11.72.3)
Q_n'=a_n Q_(n+1)+b_n Q_n,
a_n=(4n+66)Q_(8,n)/Q_(8,n+1)>0,
b_n=((4n+66)/(4n+62))(log Q_(8,n+1))'.           (11.72.4)
```

The lower cone alone cannot supply endpoint entry. The exact rational
sequence `(1,1,1/2!,...,1/15!,1/20926400000000)` has every available signed
contiguous minor through order eight strictly positive but `Q_(9,0)<0`.

For the actual heat-deformed zeta coefficients define the stable margin

```text
M_n=Q_(8,n+1)^2/(Q_(8,n)Q_(8,n+2))-1.
```

Outward-rounded 2048-bit coefficient balls through `A_1256`, including 38
retained-integral repairs, prove

```text
Q_(9,n)(-100)>0 for every 0<=n<=1240.             (11.72.5)
```

The weakest relative margin is the final cached row and remains strictly
above `1/250`.

Machine-audited companions:

```text
outputs/jensen_window_pf_compound_order9_uniform_tail_flow_reduction.md
outputs/jensen_window_pf_compound_order9_m100_prefix_certificate.md
```

## 11.73 Continuous order-nine endpoint tail

Put `k=n+8`. At the sixth stable layer, exact cancellation gives

```text
E_n=log(Q_(8,n)Q_(8,n+2)/Q_(8,n+1)^2)
   =8log(x_k)+Y_k,
Y_k=w(k-1)-2w(k)+w(k+1),                         (11.73.1)
Q_(9,n)(-100)>0 iff E_n<0.                       (11.73.2)
```

The coefficient-defect buffer proves

```text
-8log(x_k)>=1004/(125(2k+1)), k>=320,            (11.73.3)
4900/k^2<1004/(125(2k+1)), k>=1249.              (11.73.4)
```

The rebalanced complete-to-first-summand transfer through all six stable
logarithms gives

```text
|Y_k-Y_k^(1)|<550/k^2, k>=1251.                  (11.73.5)
```

The continuous first-summand theorem is certified on four exactly joined
ranges:

```text
w_1''(t)<=4200/t^2, 1250<=t<=5700,
w_1''(t)<=4200/t^2, 5700<=t<=V'(2),
w_1''(t)<=4200/t^2, 2<=u<=20,
t^2 w_1''(t)<500<4200, u>=20.                    (11.73.6)
```

The first range uses 279 root segments and 874 accepted localized blocks.
The remaining ranges use 108 compact blocks, 17,999 rational finite-ray
blocks, and one normalized asymptotic interval. Their largest scaled upper is

```text
4199.18548221032914747349727361016012002522663395440335471778<4200.
```

Thus

```text
w_1''(t)<=4200/t^2 for every real t>=1250.        (11.73.7)
```

Tent integration and (11.73.5) now give

```text
Y_k^(1)<4201/k^2,
Y_k<4201/k^2+550/k^2=4751/k^2<4900/k^2,
k>=1251.                                         (11.73.8)
```

Equations (11.73.1)--(11.73.4) prove `Q_(9,n)(-100)>0` for every
`n>=1243`. Separate retained-integral balls for `A_1257,A_1258` prove the two
splice signs `n=1241,1242`. Together with (11.72.5), this yields

```text
Q_(9,n)(-100)=H_(9,n)(-100)>0
for every integer n>=0.                           (11.73.9)
```

Machine-audited companions:

```text
outputs/jensen_window_pf_compound_order9_first_summand_curvature_bridge.md
outputs/jensen_window_pf_compound_order9_localized_lower_bridge_certificate.md
outputs/jensen_window_pf_compound_order9_nested_curvature_compact_certificate.md
outputs/jensen_window_pf_compound_order9_nested_curvature_finite_ray_certificate.md
outputs/jensen_window_pf_compound_order9_nested_curvature_asymptotic_ray_certificate.md
outputs/jensen_window_pf_compound_order9_first_summand_curvature_certificate.md
outputs/jensen_window_pf_compound_order9_m100_finite_splice_certificate.md
outputs/jensen_window_pf_compound_order9_m100_entry_certificate.md
```

## 11.74 Fixed order-nine heat propagation

The uniform eventual tail (11.72.1), completed order-eight cone, endpoint
theorem (11.73.9), and cooperative recursion (11.72.4) give, by variation of
constants and finite backward induction in the shift,

```text
Q_(9,n)(lambda)=H_(9,n)(lambda)>0
for every n>=0 and -100<=lambda<=0.               (11.74.1)
```

The fixed-order initial-minor criterion then proves every consecutive-row,
arbitrary-column signed reshaped-Hankel minor through order nine has the
required sign on the same heat interval. This closes fixed order nine only;
it does not prove the all-order endpoint hierarchy, the Jensen-window PF
bridge, RH, or `Lambda<=0`.

Machine-audited companion:

```text
outputs/jensen_window_pf_compound_order9_uniform_heat_forward_invariance_certificate.md
```

## 11.75 All-order endpoint-to-heat reduction

The repeated fixed-order propagation has an exact arbitrary-order form. Put

```text
H_(m,n)=det[A_(n+i+j)]_(0<=i,j<m), H_(0,n)=1,
epsilon_m=(-1)^binom(m,2), Q_(m,n)=epsilon_m H_(m,n). (11.75.1)
```

For each fixed `m`, Theorem 11.53 supplies a possibly order-dependent
threshold:

```text
for every fixed m exists N_m such that
Q_(m,n)(lambda)>0 for n>=N_m and -100<=lambda<=0. (11.75.2)
```

No common threshold in `m` is needed. To see this, first derive the general
flow. Let

```text
C_j=(A_(n+i+j))_(0<=i<m).
```

Under the sequence-shift derivation, determinant multilinearity gives

```text
delta(H_(m,n))=det[C_0,...,C_(m-2),C_m].          (11.75.3)
```

Every earlier shifted column duplicates its right neighbor. Splitting the
coefficient heat law

```text
A_r'=(4r+2)A_(r+1)
```

into column and row parts leaves only the final shifted column and final
shifted row:

```text
column contribution=[4(n+m-1)+2]delta(H_(m,n)),
row contribution=4(m-1)delta(H_(m,n)).            (11.75.4)
```

Therefore, at every finite order,

```text
Q_(m,n)'=c_(m,n)delta(Q_(m,n)),
c_(m,n)=4n+8m-6.                                  (11.75.5)
```

Apply the three-term flag-Plucker relation to the `m` by `(m+1)` matrix with
columns `C_0,...,C_m` and the corresponding maximal minors of its first
`m-1` rows. The Hankel index map gives

```text
Q_(m-1,n+1)delta(Q_(m,n))
 =Q_(m-1,n)Q_(m,n+1)
  +delta(Q_(m-1,n+1))Q_(m,n).                    (11.75.6)
```

Since `c_(m-1,n+1)=c_(m,n)-4`, equations (11.75.5)--(11.75.6) yield

```text
Q_(m,n)'=a_(m,n)Q_(m,n+1)+b_(m,n)Q_(m,n),
a_(m,n)=c_(m,n)Q_(m-1,n)/Q_(m-1,n+1)>0,
b_(m,n)=c_(m,n)/(c_(m,n)-4)
         *(log Q_(m-1,n+1))'.                    (11.75.7)
```

Fix `m` and assume the preceding layer is positive throughout the heat
interval and that `Q_(m,n)(-100)>0` for every shift. Choose this order's own
`N_m` from (11.75.2). For `n=N_m-1,...,0`, variation of constants gives

```text
Q_(m,n)(lambda)=E_(m,n)(lambda)
 [Q_(m,n)(-100)+integral_(-100)^lambda
  E_(m,n)(s)^(-1)a_(m,n)(s)Q_(m,n+1)(s) ds],
E_(m,n)>0.                                        (11.75.8)
```

Finite backward induction in `n` makes the integral nonnegative and proves
the complete order-`m` layer. Ordinary induction in `m`, starting from the
completed base (11.74.1), therefore proves the exact equivalence

```text
[Q_(m,n)(-100)>0 for every m>=10,n>=0]
iff
[Q_(m,n)(lambda)>0 for every m>=10,n>=0
 and -100<=lambda<=0].                            (11.75.9)
```

The reverse direction is restriction to `lambda=-100`. There is no exchange
of the quantifiers `forall m` and `exists N_m` in the forward direction: one
finite order is completed before the next order is considered.

Signed Desnanot-Jacobi also gives, for every `m>=2`,

```text
Q_(m,n)Q_(m-2,n+2)
 =Q_(m-1,n+1)^2-Q_(m-1,n)Q_(m-1,n+2),            (11.75.10)
```

because

```text
binom(m,2)+binom(m-2,2)-2binom(m-1,2)=1.
```

At this stage of the reduction, the candidate static endpoint antecedent was
the infinite hierarchy

```text
Q_(m,n)(-100)>0 for every m>=10 and n>=0,         (11.75.11)
```

equivalently continuation of strict log concavity through every signed
layer. Conditional on (11.75.11), the fixed-order initial-minor theorem
applies at every finite `m` and every heat parameter, giving all-order
consecutive-row arbitrary-column signed-Hankel positivity. Section 11.78
shows, however, that (11.75.11) is false for the actual endpoint sequence.

Even if the antecedent had held, it would not by itself prove PF-infinity of
every binomially weighted Jensen window. The sign-regular-to-Jensen-PF
conversion remains a separate open theorem; RH and `Lambda<=0` remain open.

Machine-audited companion:

```text
outputs/jensen_window_pf_all_order_endpoint_heat_reduction.md
```

## 11.76 Normalized deep-Schur endpoint coordinate

The static hierarchy (11.75.11) has a sharper exact coordinate. Since the
rigorous endpoint enclosures give `A_0(-100)>0`, put

```text
h_k=A_k(-100)/A_0(-100), k>=0,
h_k=0, k<0.                                           (11.76.1)
```

Then `h_0=1`, as required for an ordinary complete-homogeneous Schur
specialization. Reverse the columns of `H_(m,n)` and transpose. With
`N=n+m-1`, the resulting entries are

```text
A_(N-i+j)(-100), 1<=i,j<=m.
```

The column reversal has sign `epsilon_m=(-1)^binom(m,2)`, so Jacobi-Trudi
gives the exact rectangular identity

```text
Q_(m,n)(-100)
 =det[A_(N-i+j)(-100)]_(1<=i,j<=m)
 =A_0(-100)^m s_((N^m))(h), N=n+m-1.              (11.76.2)
```

In particular, the candidate endpoint hierarchy is minimally the
rectangular statement

```text
s_((N^m))(h)>0 for every m>=10 and N>=m-1.         (11.76.3)
```

The same coordinate covers arbitrary columns. For
`0<=j_0<...<j_(m-1)`, define

```text
R_(m,n)(j_0,...,j_(m-1))
 =det[A_(n+i+j_l)(-100)]_(0<=i,l<m).               (11.76.4)
```

Set, for `0<=q<m`,

```text
lambda_(q+1)=n+j_(m-1-q)+q.                        (11.76.5)
```

Then

```text
lambda_r-lambda_(r+1)
 =j_(m-r)-j_(m-r-1)-1>=0,
lambda_m=n+j_0+m-1>=m-1,

epsilon_m R_(m,n)(j_0,...,j_(m-1))
 =A_0(-100)^m s_lambda(h).                         (11.76.6)
```

Conversely, every partition

```text
lambda_1>=...>=lambda_m>=m-1
```

has the canonical inverse

```text
n=lambda_m-(m-1),
j_l=lambda_(m-l)-lambda_m+l, 0<=l<m.              (11.76.7)
```

Here `j_0=0`, and
`j_l-j_(l-1)=lambda_(m-l)-lambda_(m-l+1)+1>=1`.
Thus consecutive-row arbitrary-column signed minors are in bijection with
the deep partition cone

```text
D_m={lambda_1>=...>=lambda_m>=m-1}.                (11.76.8)
```

The depth boundary is exact rather than conventional:

```text
min_(1<=i,j<=m)(lambda_i-i+j)=lambda_m-m+1>=0.     (11.76.9)
```

Hence `D_m` is precisely the Jacobi-Trudi region that never uses an
artificial negative-index zero. Combining (11.76.2)--(11.76.8) with the
completed base through order nine and the fixed-order initial-minor theorem
proves

```text
[Q_(m,n)(-100)>0 for every m>=10,n>=0]
iff
[s_lambda(h)>0 for every m>=10,lambda in D_m].     (11.76.10)
```

The reverse implication restricts to rectangles. For the forward
implication, the endpoint hierarchy plus the known lower base supplies every
contiguous layer; the finite initial-minor theorem then supplies every
arbitrary column set at each fixed order. Thus, within this conditional
route, only the rectangles in (11.76.3) are independent analytic
obligations. Section 11.78 rigorously refutes that candidate statement at its
first previously uncompleted order.

This coordinate is not full Polya-frequency positivity. Exact rational
interval propagation of the rigorous endpoint balls for `A_0,...,A_3`
proves

```text
s_(1,1,1)(h)=h_1^3-2h_1h_2+h_3
 <-4.8484*10^(-11)<0.                              (11.76.11)
```

Therefore the actual normalized endpoint sequence is not `PF_3` and not
`PF-infinity`. The failed partition `(1,1,1)` is outside `D_3`, because its
smallest part is `1<2`. Full Schur positivity consequently imposes an
additional false inequality and cannot be substituted for the deep endpoint
target. Equation (11.76.11) alone does not decide the infinite deep cone;
the direct order-ten test in Section 11.78 does.

The audited primary literature supplies the already-used finite
initial-minor transfer and supports the Hankel/Toeplitz reversal signature.
Edrei's classification addresses the strictly larger full-PF inequality
family, while published eventual total positivity by matrix powers uses a
different coordinate. No direct closing theorem occurs in the bounded
audited primary-source set; this is not an exhaustive literature claim.

Theorem 11.75 would propagate (11.76.3), equivalently (11.76.10), over
`-100<=lambda<=0`; the later counterexample means this conditional cannot be
used for the actual endpoint sequence. The binomially weighted Jensen-window
problem remains separate. Neither RH nor `Lambda<=0` follows from this
coordinate or its failure.

Machine-audited companion:

```text
outputs/jensen_window_pf_endpoint_deep_schur_coordinate.md
```

## 11.77 Rectangular Toda coordinate and boundary obstructions

The candidate rectangles in (11.76.3) satisfy an exact discrete Toda
relation. Write

```text
tau_(m,N)=s_((N^m))(h)
          =Q_(m,N-m+1)(-100)/A_0(-100)^m.          (11.77.1)
```

Applying Desnanot-Jacobi to the rectangular Jacobi-Trudi matrix gives

```text
tau_(m+1,N) tau_(m-1,N)
 =tau_(m,N)^2-tau_(m,N-1) tau_(m,N+1).            (11.77.2)
```

Therefore, whenever `tau_(m-1,N)>0`,

```text
tau_(m+1,N)>0
iff tau_(m,N)^2>tau_(m,N-1) tau_(m,N+1).          (11.77.3)
```

Thus the next rectangular layer is exactly strict width-log-concavity of the
current row. Equation (11.77.2) is not a positivity induction: its right-hand
side is subtractive, so the known positive factors do not determine its sign.
The proposed analytic input was (11.77.3) for every row `m>=9` and
admissible width, equivalently (11.76.3) for `m>=10`. Section 11.78 shows
that this strict width-log-concavity already fails for `m=9` at four widths.

A natural shifted-tail shortcut fails at the Jacobi-Trudi boundary. Define

```text
b_k^(s)=h_(s+k)/h_s, k>=0,
b_k^(s)=0, k<0.                                   (11.77.4)
```

Then

```text
s_mu(b^(s))=h_s^(-r) s_(mu+(s^r))(h)              (11.77.5)
```

is valid for a length-`r` shape only when `mu_r>=r-1`, because only then do
all tail Jacobi-Trudi indices avoid the artificial negative-index zero. The
smallest exact failure is

```text
s_(0,0)(b^(1))-h_1^(-2)s_(1,1)(h)
 =h_0 h_2/h_1^2>0.                                (11.77.6)
```

Consequently, deep endpoint positivity cannot be promoted to ordinary
moving-tail PF by (11.77.5); the attempted translation already presupposes
the same depth condition.

There is also no generic implication from even full strict Schur positivity
to Jensen hyperbolicity. Consider the specialization

```text
H(z)=exp(z/100)/((1-z)(1-2z)).                     (11.77.7)
```

Its Plancherel component gives, for every partition `lambda`,

```text
s_lambda[X+Pl_epsilon]
 >=s_lambda[Pl_epsilon]
 =f^lambda epsilon^|lambda|/|lambda|!>0.           (11.77.8)
```

Nevertheless its shift-zero cubic Jensen polynomial has exact discriminant
`-222484532394597/2000000000000<0`; after multiplication by a positive
constant, the polynomial is

```text
6000000+54180000*x+126540900*x^2+90420901*x^3.     (11.77.9)
```

It has one real zero and a nonreal conjugate pair. Hence unweighted Schur/PF
positivity, even in strict all-shape form, is insufficient for the
binomially weighted Jensen conclusion. The remaining bridge must use extra
structure specific to the actual Xi/Phi coefficients.

Equations (11.77.2), (11.77.6), and (11.77.9) isolate three distinct issues.
The first candidate route is rejected in Section 11.78; the boundary and
generic-Schur shortcuts are rejected here. The surviving task is a weaker
Xi/Phi-specific Jensen-window theorem. Neither RH nor `Lambda<=0` is
established here.

Machine-audited companion:

```text
outputs/jensen_window_pf_deep_schur_toda_boundary_gate.md
```

## 11.78 First-open-order endpoint counterexample and route correction

The first endpoint layer not covered by the completed base through order nine
can be decided directly. At `lambda=-100`, the signed Desnanot-Jacobi identity
specialized to order ten is

```text
Q_(10,n) Q_(8,n+2)
 =Q_(9,n+1)^2-Q_(9,n)Q_(9,n+2).                  (11.78.1)
```

Because the order-eight and order-nine layers are strictly positive for every
shift, define

```text
L_n=Q_(9,n+1)^2/[Q_(9,n)Q_(9,n+2)]-1.            (11.78.2)
```

Then `sign Q_(10,n)=sign L_n`. Independent `4096`-bit Arb evaluations of the
direct `10` by `10` Hankel determinant, the Jacobi-Trudi determinant, and the
Toda residual certify

```text
L_0<0, L_1<0, L_2<0, L_3<0, L_4>0.               (11.78.3)
```

The enclosing midpoints, with radii below `6*10^(-66)`, are respectively

```text
-0.032923664115365083326
-0.034086141684355165823
-0.023431440617663769611
-0.005349707526967545098
 0.013273398491067584882.                         (11.78.4)
```

Since `epsilon_10=(-1)^45=-1` and

```text
Q_(10,n)(-100)=A_0(-100)^10 s_(((n+9)^10))(h),    (11.78.5)
```

equations (11.78.3)--(11.78.5) rigorously give

```text
s_((N^10))(h)<0 for N=9,10,11,12,
s_((13^10))(h)>0.                                 (11.78.6)
```

All four negative shapes are inside the required deep cone because
`N>=10-1`. A cancellation-preserving Arb condensation scan further certifies
`L_n>0` for every `4<=n<=1240`, with no inconclusive rows, but that finite
positive stretch cannot repair an all-shift assertion. Consequently

```text
s_((N^m))(h)>0 for every m>=10,N>=m-1             (11.78.7)
```

is false for the actual endpoint sequence. Equivalently, the antecedent in
(11.75.9) and the deep-cone statement in (11.76.10) fail. The conditional
endpoint-to-heat equivalence, the rectangular Toda identity, and all results
through order nine remain valid.

This is a counterexample to the proposed all-order signed-Hankel/deep-Schur
route, not to RH, Jensen hyperbolicity, or `Lambda<=0`. The programme must
therefore seek a weaker Xi/Phi-specific condition that controls the
binomially weighted Jensen windows without assuming all-shift signed-Hankel
positivity.

Machine-audited companion:

```text
outputs/jensen_window_pf_endpoint_order10_counterexample.md
```

## 11.79 Delayed order-ten recovery at lambda zero

The failure in Section 11.78 is an endpoint failure at four shifts, not a
claim that order ten must remain negative under the heat flow.  The positive
endpoint ray beginning at `n=4` can be completed without discarding those
four signs.

For the first Newman summand, independent interval covers of the localized,
compact, finite-saddle, and asymptotic-saddle ranges compose to give

```text
z_1''(t)<=4200/t^2 for every real t>=1251.         (11.79.1)
```

The exact seven-layer complete-to-first-summand transfer gives

```text
|Z_k-Z_k^(1)|<10/k^2 for every integer k>=1252.   (11.79.2)
```

Equations (11.79.1)--(11.79.2), the seventh stable-coordinate endpoint
factorization, and the two retained-integral splice signs at `n=1241,1242`
join the direct positive block to the analytic tail.  They prove the complete
endpoint sign chart

```text
Q_(10,n)(-100)<0 for n=0,1,2,3,
Q_(10,n)(-100)>0 for every integer n>=4.           (11.79.3)
```

The completed order-nine heat layer supplies the positive lower coefficient
in the cooperative equation.  Together with the fixed-order eventual tail,
variation of constants permits a finite descending induction beginning at
any delayed shift `n0`.  Specializing to `m=10,n0=4`, (11.79.3) therefore
implies

```text
Q_(10,n)(lambda)>0 for every n>=4
and every -100<=lambda<=0.                         (11.79.4)
```

No sign for `n<4` is used in this propagation.  Independently, outward-rounded
520-digit Arb determinant rebuilds at lambda zero prove

```text
Q_(10,n)(0)>0 for n=0,1,2,3.                       (11.79.5)
```

The disjoint statements (11.79.4)--(11.79.5) now close the fixed layer:

```text
Q_(10,n)(0)>0 for every integer n>=0.              (11.79.6)
```

Combining (11.79.6) with the completed lower layers and the fixed-order
initial-minor transfer yields

```text
epsilon_k R_(k,n)(j_1,...,j_k)(0)>0
for 1<=k<=10, n>=0, and 0<=j_1<...<j_k.            (11.79.7)
```

Thus the lambda-zero signed-Hankel and arbitrary-column frontier advances
from order nine to order ten even though the all-shift endpoint hierarchy is
false.  This is a fixed-order theorem only: it proves neither order eleven
nor PF-infinity, and the Xi/Phi-specific Jensen bridge remains open.  In
particular, neither RH nor `Lambda<=0` is established here.

Machine-audited companions:

```text
outputs/jensen_window_pf_compound_order10_first_summand_curvature_certificate.md
outputs/jensen_window_pf_compound_order10_m100_delayed_entry_certificate.md
outputs/jensen_window_pf_delayed_cooperative_heat_tail_lemma.md
outputs/jensen_window_pf_compound_order10_lambda0_prefix_certificate.md
outputs/jensen_window_pf_compound_order10_lambda0_completion_certificate.md
```

## 11.80 Order-eleven reduction to one continuum theorem

The same delayed-ray strategy advances one further order without assuming a
false all-shift endpoint cone.  Direct stable-coordinate and `11` by `11`
determinant certificates give

```text
Q_(11,n)(-100)>0 for every integer 0<=n<=1242.     (11.80.1)
```

For `k=n+10`, the next exact factorization is governed by

```text
X(t)=9B(t)-z(t-1)+2z(t)-z(t+1),
y(t)=2z(t)-w(t)+log(1-exp(-X(t))),
E_n=10 log(x_k)+Y_k.                              (11.80.2)
```

The inherited order-ten theorem supplies a positive inverse-linear floor for
`X`.  An exact eighteen-row power envelope then proves

```text
|Y_k-Y_k^(1)|<37/k^2 for every integer k>=1253.   (11.80.3)
```

Consequently the single continuous target

```text
y_1''(t)<=6000/t^2 for every real t>=1252         (11.80.4)
```

would imply `Q_(11,n)(-100)>0` for every `n>=1243`; this meets (11.80.1)
without an index gap.  The shifted cooperative theorem would then propagate
the `n>=4` ray through the heat interval, while independent Arb determinants
already prove `Q_(11,n)(0)>0` for `n=0,1,2,3`.

Thus every algebraic, finite, and heat-flow input for fixed order eleven is
ready except (11.80.4).  Eight rigorous point jets support the target and a
localized H24 interval core has been validated, but neither is a continuous
half-line proof.  Order eleven, PF-infinity, the Jensen bridge, RH, and
`Lambda<=0` therefore remain unproved.

Machine-audited companions:

```text
outputs/jensen_window_pf_compound_order11_m100_prefix_certificate.md
outputs/jensen_window_pf_compound_order11_lambda0_prefix_certificate.md
outputs/jensen_window_pf_compound_order11_first_summand_point_scout.md
outputs/jensen_window_pf_compound_order11_curvature_bridge_target.md
```
