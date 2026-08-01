# Newman Theta-Curvature Probability/Operator Gate

Date: 2026-07-24

Status: exact theta-curvature probability/operator reduction with
a generic double-contact guard. This is not a proof of
`Lambda <= 0` or RH.

```text
work/rh_compute/results/jensen_window_pf_newman_theta_curvature_probability_operator_gate.json
python work/rh_compute/scripts/jensen_window_pf_newman_theta_curvature_probability_operator_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_theta_curvature_probability_operator_gate.py
```

Current result:

```text
validated Newman theta-curvature probability/operator gate: 20 rows, 0 issues, 3 theta-primitive identities, 2 monotone-convex inequalities, 1 probability law, 1 fixed-weight theta mixture, 1 uniform C2 dominant-block budget, 1 dominant-mass endpoint-jet guard, 2 transform identities, 1 oscillator factorization, 1 Doob diffusion identity, 1 Sturm nodal-loss guard, 1 asymptotic drift-curvature guard, 1 characteristic contact reduction, 1 componentwise contact decomposition, 1 explicit C1 tail disjunction, 1 endpoint Laplace comparison, 1 generic double-contact guard, 1 open Xi transversality handoff
```

## Curvature Probability

For the positive theta primitive, set

```text
R(u)=sum_(n>=1) exp(u-pi*n^2*exp(4u)), u>=0
R''(u)-R(u)=8*Phi(u)
R'(0)=-1/2
S_t(u)=exp(t*u^2)*R(u), 0<=t<=1/5
```

Writing `y=4*pi*n^2*exp(4u)`, every summand obeys

```text
s_(n,t)'/s_(n,t)=1+2*t*u-y
s_(n,t)''/s_(n,t)=(y-1-2*t*u)^2+2*t-4*y
1+2*t*u<=1+(2/5)u<=y/12, hence y-1-2*t*u>=11*y/12
s_(n,t)''/s_(n,t)>=(121/144)*y^2-4*y=y*((121/144)*y-4)>0
```

Thus `S_t` is positive, strictly decreasing, and strictly convex
uniformly for `0<=t<=1/5`. Its curvature is a probability law:

```text
dmu_t(v)=2*S_t''(v)dv on [0,infinity)
mu_t([0,infinity))=2*(S_t'(infinity)-S_t'(0))=1
S_t(u)=(1/2)*integral_[0,infinity)(v-u)_+ dmu_t(v)
```

The probability law has a fixed arithmetic component split:

```text
dmu_(n,t)(v)=2*s_(n,t)''(v)dv
w_n=mu_(n,t)([0,infinity))=2*(4*pi*n^2-1)*exp(-pi*n^2)
The weights w_n are independent of t and sum_n w_n=1.
dnu_(n,t)=dmu_(n,t)/w_n, A_t(x)=sum_(n>=1)w_n*a_(n,t)(x)
sum_(n>=2)w_n<1/2800, hence w_1>2799/2800
w_1=2*(4*pi-1)*exp(-pi)=0.99965638867482904249949954238374
```

Integration by parts gives a uniform dominant-block budget:

```text
integral v*dmu_(n,t)=2*exp(-pi*n^2)
integral v^2*dmu_(n,t)=4*integral_0^infinity s_(n,t)(v)dv
integral v^3*dmu_(n,t)=12*integral_0^infinity v*s_(n,t)(v)dv
s_(n,t)(v)<=exp(-pi*n^2)*exp(-(4*pi*n^2-1)*v)
|A_t-w_1*a_(1,t)|<1/2800
|A_t'-w_1*a_(1,t)'|<1/137200
|A_t''-w_1*a_(1,t)''|<1/3361400
|A_t'''-w_1*a_(1,t)'''|<1/54900000
```

The mass hierarchy cannot be promoted to a global one-block
approximation. Modular evenness gives the exact endpoint identities

```text
R^((2k+1))(0)=R'(0)=-1/2 for every k>=0
s_(1,0)^(5)(0)=exp(-pi)*(-1024*pi^5+11520*pi^4-33920*pi^3+26400*pi^2-3124*pi+1)>300
sum_(n>=2)s_(n,0)^(5)(0)=-1/2-s_(1,0)^(5)(0)<-601/2
```

Thus a tail carrying less than `1/2800` of the probability mass
cancels a fifth endpoint jet of magnitude above `300`.

The polynomial coefficients in `J_t` still amplify fixed absolute
errors at high frequency, so this is a rigorous compact-frequency
handoff rather than a global positivity proof.

## Positive Primitive Transform

Let

```text
C_t(x)=integral_0^infinity S_t(u)*cos(x*u)du
A_t(x)=integral_[0,infinity)cos(x*v)dmu_t(v)
```

Twice integrating the triangular mixture gives

```text
C_t(x)=(1/(2*x^2))*integral(1-cos(x*v))dmu_t(v)=(1-A_t(x))/(2*x^2), x!=0
C_t(x)>0 for every real x; for x!=0, 0<x^2*C_t(x)<1
```

This is a genuine strict Fourier-positivity theorem, but for the
theta primitive rather than for `K_(1,t)`.

## Operator And Contact

The endpoint-subtracted operator has the exact factorization

```text
D_t=-4*t^2*d_x^2+4*t*x*d_x+(2*t-1-x^2)
Q_t=2*t*d_x-x
D_t=-(Q_t^2+1)
H_t(x)=1/16+D_t[C_t](x)/8=(1-2*(Q_t^2+1)C_t(x))/16
```

Because both numerator and denominator solve the same backward
heat equation, the positive-normalizer ratio also has an exact
evolution law:

```text
partial_t H_t=-partial_x^2 H_t and partial_t C_t=-partial_x^2 C_t
G_t=8*H_t/C_t
partial_t G_t=-partial_x^2 G_t-2*(partial_x log C_t)*partial_x G_t=-C_t^(-2)*partial_x(C_t^2*partial_x G_t)
For tau=T-t, partial_tau G=C^(-2)*partial_x(C^2*partial_x G).
A_t f=-C_t^(-2)*partial_x(C_t^2*partial_x f) is nonnegative in L2(C_t^2 dx): <f,A_t f>_(C_t^2)=integral C_t^2*(f')^2 dx.
```

This weighted Sturm form does not itself exclude nodal loss. The
exact calibration model

```text
C_t(x)=1, H_t(x)=(x^2+a-2t)/8, G_t(x)=x^2+a-2t
At t_*=a/2, G has a double zero at 0; for t>t_* it has two simple real zeros, while for t<t_* it has none.
Under tau=T-t, the positive heat semigroup loses two real nodes at the multiple zero. This is allowed by the Sturm zero-number theorem.
```

shows that a nonnegative generator and the zero-number theorem
permit two nodes to disappear at a multiple zero. An Xi-specific
conserved nodal flux or quantitative transversality is still
required.

The modular endpoint jets also determine the Doob weight at high
frequency:

```text
S_t'(0)=-1/2 and S_t'''(0)=R'''(0)+6*t*R'(0)=-1/2-3*t
C_t(x)=1/(2*x^2)-(1/2+3*t)/x^4+O_t(x^-6)
(log C_t)''=2/x^2-6*(1+6*t)/x^4+O(x^-6)>0 for all sufficiently large x, uniformly on 0<=t<=1/5.
2*(log C_t)'=-4/x+4*(1+6*t)/x^3+O(x^-5)
```

Thus `C_t` is eventually log-convex, uniformly on the target
window. A global contracting-drift or positive Bakry-Emery
curvature argument is therefore unavailable.

Since `C_t>0`, it supplies a global positive normalizer. More
concretely, with `J_t=16*x^4*H_t`,

```text
J_t=4*t^2*x^2*A_t''-4*t*x*(x^2+4*t)*A_t'+(x^4+(6*t+1)*x^2+24*t^2)*A_t-((6*t+1)*x^2+24*t^2)
Because H_t(0)>0, a boundary multiple zero c is nonzero; H_t(c)=H_t'(c)=0 iff J_t(c)=J_t'(c)=0.
```

The same identity can be read as the explicit oscillatory
expectation

```text
J_t=E_mu[((x^4+(6t+1)x^2+24t^2)-4t^2x^2V^2)*cos(xV)+4tx(x^2+4t)V*sin(xV)]-((6t+1)x^2+24t^2)
```

The fixed theta mixture makes the contact equation additive:

```text
A_(n,t)(x)=integral cos(xv)dmu_(n,t)(v), C_(n,t)=(w_n-A_(n,t))/(2*x^2), H_(n,t)=w_n/16+D_t[C_(n,t)]/8
J_t=sum_(n>=1)J_(n,t), J_(n,t)=16*x^4*H_(n,t)
P_t=x^4+(6*t+1)*x^2+24*t^2, R_t=(6*t+1)*x^2+24*t^2
J_(n,t)=4*t^2*x^2*A_(n,t)''-4*t*x*(x^2+4*t)*A_(n,t)'+P_t*A_(n,t)-R_t*w_n
J_(n,t)'=4*t^2*x^2*A_(n,t)'''-4*t*x*(x^2+2*t)*A_(n,t)''+(P_t-4*t*(3*x^2+4*t))*A_(n,t)'+(4*x^3+2*(6*t+1)*x)*A_(n,t)-2*(6*t+1)*x*w_n
C_(1,0)(x)=Re[(1/4)*pi^(-(1+i*x)/4)*Gamma((1+i*x)/4,pi)]
```

Writing `J_(1,t)` for the unnormalized first component, the
remaining theta tail obeys the explicit error bars

```text
delta_0=1/2800, delta_1=1/137200, delta_2=1/3361400, delta_3=1/54900000
B_0(t,x)=4*t^2*x^2*delta_2+4*t*x*(x^2+4*t)*delta_1+(x^4+2*((6*t+1)*x^2+24*t^2))*delta_0
B_1(t,x)=4*t^2*x^2*delta_3+4*t*x*(x^2+2*t)*delta_2+|P_t-4*t*(3*x^2+4*t)|*delta_1+(4*x^3+4*(6*t+1)*x)*delta_0
|J_t-J_(1,t)|<B_0(t,x) and |J_t'-J_(1,t)'|<B_1(t,x), x>0
```

Consequently the following is a rigorous pointwise exclusion test:

```text
At any 0<t<=1/5 and x>0, if |J_(1,t)(x)|>B_0(t,x) or |J_(1,t)'(x)|>B_1(t,x), then (J_t(x),J_t'(x))!=(0,0).
If the pointwise disjunction holds for every 0<t<=1/5 and x>0, then the boundary-contact criterion gives Lambda<=0.
```

The implication is proved, but the raw moment bars grow too
quickly to serve as an all-frequency estimate. The open theorem
requires a compact certificate joined to a modularly coupled
high-frequency argument.

At the undeformed endpoint this collapses to a Laplace reference:

```text
16*x^2*H_0(x)=(1+x^2)*A_0(x)-1
1/(1+x^2) is the cosine characteristic function of the unit exponential law on [0,infinity).
H_0(c)=H_0'(c)=0 iff A_0(c)=1/(1+c^2) and A_0'(c)=-2*c/(1+c^2)^2.
```

Thus a multiple endpoint zero is exactly tangential contact between
the theta-curvature characteristic function and that of the unit
exponential law.

## Nonpromotion Guard

The positive normalization is not itself the missing theorem.
Let

```text
f=(1-|.|)_+ * exp(-(.^2)/8), a smooth positive even Gaussian-tail kernel
g(u)=-4*integral_0^infinity (exp(-|u-v|)+exp(-(u+v)))*f(v)dv
g''-g=8*f, g'(0)=0, g<0
R_epsilon(u)=exp(-u)/2+epsilon*g(u). For sufficiently small epsilon>0, R_epsilon is positive and decreasing, R_epsilon''=R_epsilon+8*epsilon*f>0, and R_epsilon'(0)=-1/2.
```

The kernel has a double Fourier zero at `2*pi`. For small positive
`epsilon`, its Neumann lift has all the curvature-probability and
primitive-transform positivity properties above, while the
associated transform still has that multiple zero. Hence generic
probability positivity, convexity, and `C_t>0` cannot be promoted
to the Newman conclusion.

## Live Handoff

Prove uniformly for 0<t<=1/5 that the Xi theta-curvature characteristic expression J_t has no common real zero with J_t'. Equivalently, use the explicit positive probability density dmu_t=2*d_u^2[exp(tu^2)R(u)]du to obtain a quantitative C1 separation for the displayed oscillatory expectation. The generic Neumann-lift guard forbids promotion from probability positivity, convexity, or C_t>0 alone. Any successful estimate must preserve the theta arithmetic and remain uniform as t tends to zero. The separate compact-transversality certificate now proves the displayed first-block disjunction on 1/4<=x<=38 and combines it with an exact origin collar. The live obligation is therefore x>38, where the raw B_0/B_1 bars must be abandoned in favor of the corrected Riemann-Siegel phase-aware partition.
