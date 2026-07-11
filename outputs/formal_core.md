# Formal Core

Date: 2026-07-11

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

Together with simplicity for every `t>Lambda`, this yields the exact
equivalence

```text
Lambda<=0
iff H_t has only simple zeros for every t>0.
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

For `0<t<=1/2`, set

```text
phi_t(u)=exp(t*u^2)*Phi(u),
H_t(x)=integral_0^infinity phi_t(u)cos(xu)du,
L_t(x)=H_t'(x)^2-H_t(x)H_t''(x).
```

Lemma 11.9 and the simple real-zero factorization give the exact equivalence

```text
Lambda<=0
iff L_t(x)>0 for every real x and every 0<t<=1/2.
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
iff for every 0<t<=1/2, the translations of K_(1,t)
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
