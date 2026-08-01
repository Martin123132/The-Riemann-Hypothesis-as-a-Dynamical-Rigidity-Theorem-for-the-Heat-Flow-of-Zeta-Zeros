# Score-Abel Radial Dimension-Lift Gate

Date: 2026-07-25

Status: exact radial dimension ladder for every shifted Jensen
window, with a sharp two-Gaussian nonpromotion gate. This is
not a proof of PF-infinity, RH, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_newman_score_abel_radial_dimension_lift_gate.json
python work/rh_compute/scripts/jensen_window_pf_newman_score_abel_radial_dimension_lift_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_score_abel_radial_dimension_lift_gate.py
```

Current result:

```text
validated score-Abel radial dimension-lift gate: 15 rows, 0 issues, 1 imported Abel coordinate, 1 radial probability ladder, 1 planar marginal, 1 Bessel derivative identity, 1 own-dimension positive-definiteness theorem, 1 dimension walk, 1 Jensen identification, 1 Schoenberg scope guard, 3 exact Gaussian countermodel rows, 2 imported Xi low-degree closures, 1 imported degree-361 closure, 1 nonpromotion gate
```

## Radial Probability Ladder

Start from the exact Abel representation already obtained:

```text
Let mathcal_F_t(z)=2*H_t(sqrt(-z))=sum_(k>=0) A_k(t)*z^k/k!. For the established Abel density r_t, M_k=integral_0^infinity v^k*r_t(v)dv=k!*A_k.
```

For every shift `n>=0`, define

```text
For n>=0 and y in R^(2n+2), define g_(n,t)(y)=r_t(|y|^2/4)/((4*pi)^(n+1)*A_n(t)). Then g_(n,t)>=0 and integral_R^(2n+2) g_(n,t)(y)dy=1.
Under g_(n,t), V=|Y|^2/4 has probability law dmu_(n,t)(v)=v^n*r_t(v)dv/M_n(t). Thus each shift n uses the n-th size bias of the same Abel measure.
```

Indeed, the sphere area in dimension `2n+2` and the change
`v=|y|^2/4` leave the factor `1/n!`; the remaining radial
moment is `M_n=n!*A_n`. Thus the density has mass one.

At the bottom of the ladder there is a direct geometric
interpretation:

```text
For n=0, g_(0,t) is the isotropic planar Abel/Radon lift of f_t/A_0: integral_R g_(0,t)(u,y)dy=f_t(u)/A_0. A direct proof swaps the Abel integrals and uses integral_u^s r*dr/(sqrt(r^2-u^2)*sqrt(s^2-r^2))=pi/2.
```

So the Abel transform was not merely a moment device: it is the
isotropic planar density whose one-coordinate marginal is the
normalized Newman kernel.

## Bessel-Fourier Identity

Termwise differentiation gives

```text
mathcal_F_t^(n)(z)=1/n!*integral_0^infinity v^n*0F1(;n+1;v*z)*r_t(v)dv, and 0F1(;n+1;-v*x^2)=n!*(x*sqrt(v))^(-n)*J_n(2*x*sqrt(v)), with the continuous value at x=0.
```

The right side is precisely the spherical characteristic
kernel in dimension `2n+2`. Therefore

```text
With the Fourier convention exp(i*xi dot y), hat(g_(n,t))(xi)=mathcal_F_t^(n)(-|xi|^2)/A_n(t). Consequently every normalized derivative is radial positive definite in its own dimension 2n+2.
If phi_(n,t)(x)=mathcal_F_t^(n)(-x^2)/A_n(t), then phi_(n+1,t)(x)=-A_n(t)*phi_(n,t)'(x)/(2*x*A_(n+1)(t)), x>0, with continuous extension at 0.
```

This makes every derivative a characteristic function, but in
its own increasing dimension and under its own size-biased law.

## Jensen Windows

```text
P_(D,n)(w)=sum_(k=0)^D binom(D,k)A_(n+k)w^k is exactly the degree-D Jensen polynomial of mathcal_F_t^(n), equivalently of the dimension-(2n+2) radial characteristic profile after normalization.
```

The radial ladder therefore captures every shift, not only the
unshifted window. It also tells us exactly why positive
definiteness alone is not the missing theorem.

## Schoenberg Scope Guard

```text
The all-dimensions Schoenberg theorem concerns one fixed radial profile in every dimension. Here both phi_n and its size-biased radial law change with n. It therefore gives no complete-monotonicity or Laguerre-Polya promotion; indeed phi_(0,0) has real Xi zeros and cannot be a nonzero Gaussian mixture valid in every dimension.
```

The fixed-profile all-dimension hypothesis is absent. Treating
the sequence `phi_n` as one Schoenberg profile would silently
replace the actual size-biased family by a stronger false
premise.

## Two-Gaussian Newman Countermodel

The failure can be made exact while retaining the same heat-flow
deformation:

```text
For c>=1, 0<=t<=1/5, let a=c^2-t, b=2*c^2-t and f^G_(t,c)(u)=exp(-a*u^2)+exp(-b*u^2)=exp(t*u^2)*(exp(-c^2*u^2)+exp(-2*c^2*u^2)). Then partial_t f^G=u^2*f^G.
Writing q=exp(-c^2*u^2), (log f^G)''=-2*a-2*c^2*q/(1+q)+4*c^4*u^2*q/(1+q)^2 <= -((2-4/e)*c^2-2/5)<0. Hence the family is uniformly strongly log-concave, and the certified curvature can be made arbitrarily large by increasing c.
r^G_(t,c)(v)=4*sqrt(pi*a)*exp(-4*a*v)+4*sqrt(pi*b)*exp(-4*b*v). It is positive on v>0 and therefore generates the same normalized radial lifts and dimension walk in every shift.
```

Thus the countermodel has every radial lift above, full Abel
support, and an arbitrarily large certified strong-concavity
constant. Nevertheless,

```text
mathcal_G_(t,c)(z)=sqrt(pi/a)*exp(z/(4*a))+sqrt(pi/b)*exp(z/(4*b)).
The complete zero set of mathcal_G_(t,c) is z_k=2*a*b*log(a/b)/c^2+4*pi*i*a*b*(2*k+1)/c^2, k in Z. Every zero is nonreal.
Put p=sqrt(pi/a), q=sqrt(pi/b), alpha=1/(4*a), beta=1/(4*b), so A_k=p*alpha^k+q*beta^k. Then disc P_(2,n)=4*(A_(n+1)^2-A_n*A_(n+2))=-4*p*q*(alpha*beta)^n*(alpha-beta)^2<0 for every n>=0.
```

So it fails the degree-two Jensen test for every shift. The
obstruction is stronger than a single accidental polynomial:
the entire generating function has a complete nonreal zero
lattice.

## Quadratic And Cubic Reconciliation

The Gaussian countermodel does not conflict with the existing
Xi low-degree theorems:

```text
The Gaussian model does not satisfy the Xi-specific input used below: y maps to f^G_(t,c)(sqrt(y))=exp(-a*y)+exp(-b*y), whose logarithmic second derivative is (a-b)^2*exp(-(a+b)*y)/(exp(-a*y)+exp(-b*y))^2>0. It is strictly log-convex, not log-concave, in y.
The established Xi theorem says y maps to Phi(sqrt(y)) is strictly log-concave on [0,infinity). Multiplication by exp(t*y) preserves this for every real t. If E_k(t)=integral_0^infinity y^(k-1/2)*exp(t*y)*Phi(sqrt(y))dy, then A_k(t)=sqrt(pi)*4^(-k)*E_k(t)/Gamma(k+1/2). Berwald-Borell therefore gives A_(n+1)^2>=A_n*A_(n+2) for every n>=0 and real t. Since M_n=n!*A_n, equivalently M_(n+1)^2>=(n+1)/(n+2)*M_n*M_(n+2), so every shifted degree-two Jensen polynomial is hyperbolic.
The validated reciprocal-defect entry and forward-uniform heat theorem separately prove that every shifted degree-three Jensen polynomial is hyperbolic at every finite t>=-100, in particular throughout 0<=t<=1/5. This conclusion uses the Xi coefficient heat trajectory and is not a consequence of the radial ladder or of quadratic log-concavity.
```

Thus the Abel concentration inequality in degree two is not
open, and the separate reciprocal-defect heat theorem also
closes degree three throughout the target heat interval.

## Degree-361 Reconciliation

```text
A separate Xi real-zero-band and sector theorem now closes every shifted Jensen layer through degree 361 uniformly on 0<=t<=1/5. Thus the outer-root quartic threshold u<=U(a,p) and the adjacent quintic obstruction are no longer open for the actual Xi heat flow on this interval. The radial ladder remains an exact coordinate, not the source of that theorem. The first unproved direct degree is 362, but one more bounded degree is not the terminal need: the live target is an unbounded cofinal degree sequence or an all-degree zero-preserving theorem. A separate strong-log-concave local Mellin countermodel satisfies all neighboring ratio and cubic signs while failing its quartic, and an exact length-ten positive rational segment satisfies the strengthened scaled-defect and reciprocal-increment corridors together with all 120, 126, and 56 supported signed-Hankel conditions of orders two, three, and four across multiple shifts while still having u>U and a nonhyperbolic adjacent quintic. A downstream exact Bernstein certificate proves that every signed order-three/order-four continuation of this fixed prefix through x_12 fails the necessary increasing-x_13 compatibility margin, so this witness cannot become an infinite countermodel. A distinct exact rational tail from the same outer contact nevertheless clears every stated scalar gate and all 364, 715, and 792 finite signed minors of orders two through four through x_13, where its compatibility margin is positive. Its own x_14 compatibility is negative, but only for that selected tail. A later length-14 survivor passes 4,043 supported signed minors through order eight while its adjacent quintic is nonhyperbolic. These remain valid generic nonpromotion guards, while the Xi sector theorem excludes them through degree 361 by theta-kernel geometry. The surviving condition must control unbounded degree.
```

The radial theorem is therefore useful as a coordinate and a
theorem-search filter. The Xi sector theorem now closes the
earlier quartic and quintic frontier, and indeed every shifted
degree through 361, without promoting the generic radial facts.
The direct finite-degree frontier begins at degree 362. The
mathematically relevant next step is an unbounded cofinal
sequence or a genuine all-degree theorem.

References: https://dlmf.nist.gov/10.39 and
https://doi.org/10.2307/1968466; for modern dimension-walk
context see https://arxiv.org/abs/2408.11612.
