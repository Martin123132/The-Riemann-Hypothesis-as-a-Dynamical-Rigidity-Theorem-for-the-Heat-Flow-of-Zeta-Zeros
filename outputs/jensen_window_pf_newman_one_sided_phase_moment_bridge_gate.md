# Jensen-Window PF Newman One-Sided Phase/Moment Bridge Gate

Date: 2026-07-25

Status: exact zero-free one-sided phase lift and score-probability
bridge to the signed-Hankel coefficients and their positive Abel
measure. The Xi joint-avoidance theorem remains open;
this is not a proof of `Lambda <= 0` or RH.

```text
work/rh_compute/results/jensen_window_pf_newman_one_sided_phase_moment_bridge_gate.json
python work/rh_compute/scripts/jensen_window_pf_newman_one_sided_phase_moment_bridge_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_one_sided_phase_moment_bridge_gate.py
```

## Score Probability

Import the proved positive-time kernel shape:

```text
For 0<=t<=1/5, f_t(u)=exp(t*u^2)*Phi(u) is positive, even, C-infinity, strictly decreasing for u>0, and super-exponentially decaying; (log f_t)''<=-(kappa-2t), kappa=-Phi''(0)/Phi(0)>74.9076
dnu_t(u)=-f_t'(u)du/f_t(0), u>0, is a nondegenerate probability law with a continuous positive density
```

The probability is not an auxiliary model. It is obtained directly
from the logarithmic slope of the positive half-kernel.

## Zero-Free Complex Lift

Integration by parts gives

```text
F_t(x)=integral_0^infinity f_t(u)exp(i*x*u)du=H_t(x)+i*Y_t(x); chi_t(x)=integral_0^infinity exp(i*x*u)dnu_t(u); F_t(x)=i*f_t(0)*(1-chi_t(x))/x for x>0
Writing chi_t=C_t+i*B_t gives H_t=f_t(0)*B_t/x and Y_t=f_t(0)*(1-C_t)/x>0 for x>0
```

Because the score law has a continuous positive density, `C_t(x)<1`
for every `x>0`. Hence `Y_t(x)>0`, so this exact complex lift never
vanishes on the positive real axis.

## Exact Contact System

```text
For Theta_t=arg(F_t) in (0,pi), at H_t(x)=0 one has Theta_t=pi/2, Theta_t'=-H_t'/Y_t=-B_t'/(1-C_t); B_t=E_nu[sin(xU)] and B_t'=E_nu[U*cos(xU)]; therefore H_t=H_t'=0 iff B_t=B_t'=0
T_L[H_t]=H_t^2+(H_t'/L)^2=f_t(0)^2/x^2*(B_t^2+(B_t'-B_t/x)^2/L^2)
```

This removes the exceptional `E=0` branch from the phase reduction.
The surviving obligation is the joint small-ball problem for the
sine characteristic component and its derivative.

## Signed-Hankel Coordinate

The same probability law encodes the coefficient workstream:

```text
If mu_(2k)(t)=integral_R u^(2k)f_t(u)du and c_k(t)=mu_(2k)(t)/(2k)!, then E_nu[U^(2k+1)]=(2k+1)*mu_(2k)/(2*f_t(0)), c_k=2*f_t(0)*E_nu[U^(2k+1)]/(2k+1)!, A_k=k!*c_k, and 2*H_t(sqrt(-z))=sum_(k>=0)c_k*z^k
For q=s*(1-s) and P_(D,n)(w)=sum_(k=0)^D binom(D,k)A_(n+k)w^k, P_(D,n)(w)=2*f_t(0)*D!/(D+n)!*integral_0^1 E_nu[U^(2n+1)*q^n*L_D^(n)(-q*U^2*w)]ds. For every U>0 and 0<s<1, the fixed kernel has D simple negative w-roots; the remaining issue is preservation under this specific positive score/Beta scale mixture
```

Thus the phase-critical and signed-Hankel programmes are two
coordinates on the same score law. The identity itself does not
promote the positive Laguerre scale mixture or the finite signed
minors to PF-infinity.

## Score-Beta Abel Measure

The two-variable score/Beta mixture can be pushed to one
positive measure on the Laguerre scale:

```text
Let S be uniform on (0,1), independent of U~nu_t, Q=S*(1-S), and define the finite measure rho_t(E)=2*f_t(0)*E_nu[U*1_(Q*U^2 in E)]. Then rho_t has density r_t(v)=4*integral_(2*sqrt(v))^infinity (-f_t'(u))/sqrt(u^2-4*v)du, v>0, and M_k(t)=integral_0^infinity v^k*r_t(v)dv=k!*A_k(t)
2*H_t(sqrt(-z))=integral_0^infinity I_0(2*sqrt(v*z))*r_t(v)dv and P_(D,n)(w)=D!/(D+n)!*integral_0^infinity v^n*L_D^(n)(-v*w)*r_t(v)dv
Writing R_t(v)=integral_v^infinity r_t(s)ds, partial_t r_t(v)=4*v*r_t(v)-2*R_t(v), equivalently M_k'=2*(2*k+1)*M_(k+1)/(k+1). Moreover disc P_(2,n)=4*(A_(n+1)^2-A_n*A_(n+2))>=0 iff M_(n+1)^2/(M_n*M_(n+2))>=(n+1)/(n+2); positivity of rho_t alone gives only the opposite-side Cauchy-Schwarz upper bound M_(n+1)^2/(M_n*M_(n+2))<=1
```

This makes the first Jensen obstruction a sharp concentration
threshold for the tilted Abel measure. Ordinary moment positivity
supplies only Cauchy-Schwarz log-convexity, so it cannot be used
as the missing all-degree preservation theorem.

## Interlacing Guard

```text
Because -f_t'(u)>0 for u>0, r_t(v)>0 for every v>0. If 0<ell_1<...<ell_D are the roots of L_D^(n), then the roots of L_D^(n)(-v*w) are -ell_D/v<...<-ell_1/v. For D>=2 the family {L_D^(n)(-v*w):v>0} has no common interlacer: as v tends to 0 every root tends to -infinity, while as v tends to infinity every root tends to 0 from below. Hence a direct global common-interlacing mixture theorem cannot close the Xi target; a weighted Xi-specific total-positive or variation-diminishing connection would still be sufficient
```

The support geometry therefore retires global common interlacing
of the component family. It does not rule out a theorem using the
specific Abel weights and a stronger total-positive connection.

The bare Abel operator cannot supply that connection by itself:

```text
For K(x,y)=(y-x)_+^(-1/2), the 2x2 minor at x=(0,1), y=(2,3) is 1/2-1/sqrt(3)<0, whereas the minor at x=(0,2), y=(1,3) is 1>0. Thus the bare Abel kernel is neither TP_2 nor sign-regular of order 2; any variation-diminishing closure must use the Xi weight or a larger composed kernel
```

## Closed Score Flow

The probability coordinate remains closed under Newman time:

```text
Writing m_t(u)=-f_t'(u)/f_t(0) and bar_nu_t(u)=integral_u^infinity m_t(v)dv=f_t(u)/f_t(0), partial_t m_t=u^2*m_t-2*u*bar_nu_t. Equivalently, partial_t chi_t=-chi_t''+2*chi_t'/x+2*(1-chi_t)/x^2 and partial_t B_t=-B_t''+2*B_t'/x-2*B_t/x^2
```

This is an exact radial backward-heat equation for the sine
observable. It does not by itself prevent collision; the missing
input remains Xi-specific sign regularity or arithmetic phase
separation.

## Normalizer-Compatible Jet

The score coordinate is quantitatively compatible with the exact
normalization used by the corrected Riemann-Siegel branch:

```text
Let mathcal_A_t=|M_t((1-i*x)/2)|, Z_t=H_t/mathcal_A_t, g_t=x*mathcal_A_t/f_t(0), and a_t=(log mathcal_A_t)'. Then B_t=g_t*Z_t and B_t'-B_t/x=g_t*(Z_t'+a_t*Z_t). On L>=50, 0<tL<=25, |a_t|/L<1/2, so with m_-=(9-sqrt(17))/8 and m_+=(9+sqrt(17))/8, m_-*g_t^2*T_L[Z_t] <= B_t^2+(B_t'-B_t/x)^2/L^2 <= m_+*g_t^2*T_L[Z_t]. Equivalently, m_-*mathcal_A_t^2*T_L[Z_t] <= T_L[H_t] <= m_+*mathcal_A_t^2*T_L[Z_t]
```

Thus the asymptotically flat raw phase is a conditioning effect,
not a new obstruction after the lift amplitude is restored. The
existing corrected lower-bound target remains unchanged.

## Conditioning Guard

Repeated integration by parts gives

```text
Uniformly for 0<=t<=1/5, Y_t(x)=f_t(0)/x-f_t''(0)/x^3+O(x^-5)=f_t(0)/x*(1+(kappa-2t)/x^2+O(x^-4)), while H_t and H_t' are rapidly decreasing; hence pi/2-Theta_t and Theta_t' are O_A(x^-A) for every A
```

The exact one-sided phase is globally nonvanishing but asymptotically flat. No fixed polynomial absolute lower bound for |Theta_t'| can be transferred from the corrected Riemann-Siegel phase; quantitative margins must retain their own lift amplitude.

The corrected Riemann-Siegel phase and this exact phase therefore
cannot share an absolute velocity threshold without an explicit
amplitude conversion.

## Shape Countermodel

```text
For f=K_(1,2)=(1-|.|)_+*exp(-(.^2)/8), f is positive, even, strictly decreasing and strongly log-concave, while Fourier[f](xi)=8*sqrt(2*pi)*exp(-2*xi^2)*sin(xi/2)^2/xi^2 has double zeros at xi=2*pi*n; its score probability therefore has B=B'=0 at those points
lim_(xi->2*pi) Fourier[f](xi)/(xi-2*pi)^2=sqrt(2*pi)*exp(-8*pi^2)/(2*pi^2)>0
```

So probability positivity, monotonicity, and strong log-concavity
still do not prove transversality.

## Numerical Crosscheck

| label | t | x | Re F | Im F | B | phase velocity |
|---|---:|---:|---:|---:|---:|---:|
| origin_scale | 0.0 | 1.0000000000000000000000000000000000 | 0.0617821234887595375475711022950 | 0.00536226124433314528087727480169 | 0.138308825120909150722710536595 | 0.0865011865313348111105123548922 |
| first_xi_crossing | 0.0 | 28.269450283469387580914503967124941 | 3.65266341309597454617387947171e-50 | 0.0181612061379581274434711380032 | -5.61758462520972539318165601487e-49 | 0.00475849139201126423575890656899 |
| positive_time_edge | 0.20000000000000000000 | 38.000000000000000000000000000000000 | -0.0000154478620675650912327473739218 | 0.0124455239128418939102969267037 | -0.00131413215080205734711010188599 | -0.000600591980544819600376217770342 |

The quadrature checks the exact lift and first-jet identities at
two endpoint rows, including the first Xi crossing, and one
positive-time row. These are diagnostics only.

## Live Target

```text
Prove for the Xi score laws nu_t that (B_t(x),B_t'(x))!=(0,0) for every x>38 and 0<t<=1/5, using theta arithmetic or an all-order PF/sign-regular theorem; probability, strong log-concavity, and the zero-free complex lift alone are insufficient
```

A successful next theorem must use the theta arithmetic of this
specific score law or close the all-order PF/sign-regular route.
