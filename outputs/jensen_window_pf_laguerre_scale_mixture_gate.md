# Jensen-Window PF Laguerre Scale-Mixture Gate

Date: 2026-07-11

Status: exact Laguerre scale-mixture reduction with sharp preservation
guards. This is not a proof of PF-infinity, RH, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_laguerre_scale_mixture_gate.json
python work/rh_compute/scripts/jensen_window_pf_laguerre_scale_mixture_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_laguerre_scale_mixture_gate.py
```

Current result:

```text
validated Jensen-window PF Laguerre scale-mixture gate: 10 rows, 0 issues, 3 exact kernel identities, 1 individual-kernel hyperbolicity theorem, 1 positive-mixture countermodel, 1 Gamma reduction, 1 half-integer all-degree theorem, 1 log-concavity countermodel, 1 open Xi-specific handoff
```

## Exact Kernel

With the established normalization

```text
A_j(lambda)=j!/(2j)!*integral_0^infinity W_lambda(u)u^(2j)du, W_lambda(u)=2*exp(lambda*u^2)*Phi(u)>0
J_D(w)=integral_0^infinity W_lambda(u)K_D(w,u^2)du, K_D(w,y)=sum_(j=0)^D (D)_j*(w*y)^j/(2j)!
K_D(w,y)={}_1F_1(-D;1/2;-w*y/4)
K_D(w,y)=D!/(1/2)_D*L_D^(-1/2)(-w*y/4)
```

For y>0, K_D(w,y) has D simple roots w=-4*ell_(D,r)/y<0, where ell_(D,r) are the positive roots of L_D^(-1/2).

Thus every scale entering the integral is individually hyperbolic.
That observation is exact, but the next promotion is false.

## Preservation Guards

The positive two-atom measure

```text
(9/10)*delta_(1/4)+(1/10)*delta_(5/2) in y=u^2
J_2(w)=109*w**2/1920 + 19*w/40 + 1
disc_w J_2=-7/4800<0
```

is already a counterexample to positive-mixture preservation. Even
log-concavity is insufficient: the exponential density gives

```text
J_3(w)=w**3/20 + w**2/2 + 3*w/2 + 1
disc_w J_3=-1/200<0.
```

So the route cannot close from positivity, moment existence, or
log-concavity of the squared-scale density alone.

## Half-Integer Gamma Theorem

For `Y~Gamma(alpha,theta)`,

```text
P_(D,alpha,theta)(w)={}_2F_1(-D,alpha;1/2;-theta*w/4)
```

When `alpha=m+1/2`, Euler's transformation and the Jacobi identity give

```text
P=(1-X)^(D-m)*{}_2F_1(-m,D+1/2;1/2;X), the cofactor being proportional to P_m^(-1/2,D-m)(1-2X)
P=D!/(1/2)_D*P_D^(-1/2,m-D)(1-2X)
All D roots are real and negative: min(D,m) are simple in (-4/theta,0), and w=-4/theta has multiplicity max(D-m,0).
```

This is a genuine all-degree hyperbolic benchmark family. It matches
the half-integer power structure of an even analytic weight after the
pushforward `y=u^2`, but positive mixing of these blocks is not enough.

## Live Handoff

Exploit the actual pushforward density W_lambda(sqrt(y))/(2*sqrt(y)) through a half-integer Gamma/Laguerre expansion with a weighted total-positive connection or direct variation-diminishing rule. Because the density has full support y>0 and every component root scales as 1/y, the complete fixed-scale family has no common interlacer in degree D>=2; positivity of the expansion coefficients alone is also insufficient. The fractional Abel kernel has both signs of 2x2 minors, so any variation-diminishing argument must use a weighted Xi-specific composition. The Abel measure now also gives every shift n an exact radial probability in dimension 2n+2 and a classical dimension walk. A two-Gaussian Newman flow has that entire radial ladder and arbitrarily strong log-concavity while every shifted quadratic Jensen discriminant is negative. The model is log-convex in the squared variable and therefore lies outside the Xi root-variable log-concavity theorem. That Xi theorem already closes every shifted quadratic, and the reciprocal-defect heat theorem closes every shifted cubic on the target interval. Hence radial positive definiteness cannot replace the missing Xi-specific degree-four or all-degree connection.

The useful question is now whether the actual Phi pushforward carries
a weighted total-positive Gamma/Laguerre connection that the two
countermodels cannot possess. Global common interlacing of all scales
is impossible by full support. The remaining Xi-specific theorem
target is not a generic measure argument, and it cannot come from
total positivity of the bare Abel kernel.

The exact full-support guard and score/Beta Abel coordinate are in
`outputs/jensen_window_pf_newman_one_sided_phase_moment_bridge_gate.md`.
Their all-shift radial dimension ladder and exact two-Gaussian
nonpromotion theorem are in
`outputs/jensen_window_pf_newman_score_abel_radial_dimension_lift_gate.md`.

References: https://dlmf.nist.gov/13.6, https://dlmf.nist.gov/15.8,
and https://dlmf.nist.gov/18.16.
