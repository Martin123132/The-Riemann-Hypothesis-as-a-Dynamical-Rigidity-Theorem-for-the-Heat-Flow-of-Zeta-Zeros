# Newman First-Order Centered Real-Residual Reduction

Date: 2026-07-26

Status: exact saddle-centered scalar reduction with certified
nuisance budgets. The Xi arithmetic theorem remains open; this
is not a proof of RH.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_first_order_centered_real_residual_reduction.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_real_residual_reduction.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_real_residual_reduction.py
```

## Centered Crossing Scalar

```text
E_[1]=X+iY, R_a=P_a+iQ_a, lambda_a=u_a+i*v_a, E_[1],x=lambda_a*E_[1]+R_a
U=u_a*X-v_a*Y+P_a=u_a*X+S_a
S_a=P_a-v_a*Y=Re(R_a)-v_a*Im(E_[1])=U-u_a*X
W_[1]=v_a*(X^2+Y^2)+Q_a*X-P_a*Y=X*(v_a*X+Q_a)-S_a*Y
At X=0, U=S_a and W_[1]=-S_a*Y.
```

Unlike `W_[1]*Y`, `S_a` still classifies a simple crossing when
`E_[1]=0`; no complex-main zero is deleted.

## Critical Frame

```text
s=(1-i*x)/2, L=log(x/(4*pi)), a^2=x/(4*pi)+t/16, alpha(s)=1/(2s)+1/(s-1)+(1/2)Log(s/(2*pi))
Re(alpha)=L/2+(1/4)log(1+x^(-2))-1/(1+x^2)
Im(alpha)=-(1/2)atan(x)+3x/(1+x^2), x>0
C_x=Re(alpha'(s))=(7x^2-5)/(x^2+1)^2; D_x=Im(alpha'(s))=x*(x^2+5)/(x^2+1)^2
phi'/phi=i*beta_t', beta_t'=-(1/2)Re(alpha+(t/2)alpha*alpha')
s_*'=t*D_x/4+i*(-1/2-t*C_x/4)
u_a=-(t/4)D_x*log(a); delta_a=log(a)-Re(alpha); v_a=delta_a/2+(t/4)*(C_x*delta_a+Im(alpha)*D_x)
0<D_x<=2/x and log(a)<L imply |u_a|<=tL/(2x)<=25/(2x)<exp(-L)
Put delta_0=1/(1+x^2)-(1/4)log(1+x^(-2)), z=pi*t/(4x), delta_t=(1/2)log(1+z). Then v_a=delta_0/2+(1/4)(log(1+z)-z)+(t/4)*C_x*(delta_0+delta_t)+(t/4)*(Im(alpha)*D_x+pi/(4x)).
|v_a|<3/x^2<exp(-2L)
```

The order-`t/x` pieces in `v_a` cancel after centering at the
moving Riemann-Siegel saddle. This is why `v_a` is `O(x^-2)`.

## Absolute Budgets

```text
sum_(n=1)^N |e_n|<=50*exp(L/4)
|alpha_n|<2L, |alpha'|<=7/x, t/4+t^2|alpha_n|^2/8<313, hence |d_n|<1
|K(T_0)|/A_t<20*exp(-L/4), exp(t*pi^2/64)<2, and |F+F'''/(12*pi^2*a)|<6 imply |g_0|<exp(L/4)
|E_[1]|<101*exp(L/4)
|d_(n,x)|<4223/x^2 and therefore |sum e_n*d_(n,x)|<1500*exp(-7L/4)<1e-7*exp(-5L/4)
```

## Core Arithmetic Scalar

```text
M_(1,a)=sum_(n=1)^N log(n/a)f_n, D_(1,x)=sum_(n=1)^N e_n*d_(n,x), G_a=g_(0,x)-lambda_a*g_0
R_a=-s_*'*M_(1,a)+D_(1,x)+G_a
A_a=-(1/2)Im(M_(1,a))-(t/4)*(D_x*Re(M_(1,a))+C_x*Im(M_(1,a)))+Re(G_a)
P_a=A_a+Re(D_(1,x))
For |X|<=50000*exp(-5L/4), |U-A_a|<1e-6*exp(-5L/4).
```

Thus the non-negligible term is one centered logarithmic moment
plus the retained endpoint defect.

## Endpoint Defect

```text
H_a=F(p)+F'''(p)/(12*pi^2*a), g_0=-kappa_N*H_a, K(T_0)=M_0(iT_0)U*exp(pi*i/8)
a_x=1/(8*pi*a), p_x=-1/(4*pi*a), H_(a,x)=p_x*(F'(p)+F''''(p)/(12*pi^2*a))-a_x*F'''(p)/(12*pi^2*a^2)
mu_a=kappa_(N,x)/kappa_N-lambda_a=K_x/K-(M_t(s))_x/M_t(s)+s_*'*log(a)
K_x/K=-pi/8+1/(4T_0)+1/(2(T_0+i)); (M_t(s))_x/M_t(s)=(-i/2)*(alpha+(t/2)alpha*alpha')
G_a=-kappa_N*(H_(a,x)+mu_a*H_a)
```

## Contact And Winding Target

```text
A full contact forces |X|<50000*exp(-5L/4) and |U|<100000L*exp(-5L/4).
On q=2tL^2>=1, prove that |X|<=50000*exp(-5L/4) implies |A_a|>(100000L+1)*exp(-5L/4).
The scalar theorem and |U-A_a|<1e-6*exp(-5L/4) give |U|>100000L*exp(-5L/4), excluding contact.
At X=0, U=S_a. Thus every simple upward crossing is exactly X=0,S_a>0, including E_[1]=0.
Under the sufficient scalar theorem, sign(U)=sign(A_a) at every X=0 crossing; no exceptional complex-main class is needed.
kappa_j=N_(X=0,A_a>0;t_j)-N_(X=0,A_a>0;t_(j+1))+I_i(V_j)+I_i(finite shoulders and joins), and it is enough to prove this integer is <1.
```

## Route Audit

```text
At X=0, W_[1]=-S_a*Y. The Wronskian sign loses the E_[1]=0 crossing, while S_a and the certified core A_a do not.
The Wronskian magnitude disjunction is a stronger sufficient contact theorem, not a reduction of the scalar slope difficulty.
The RH-level input is now a lower bound and signed crossing budget for A_a, an explicit centered logarithmic moment plus the exact C_0+C_1/a endpoint defect.
```

## Remaining Cells

```text
Prove the scalar band theorem and its oriented A_a-positive crossing budget on q>=1 in the live 0<tL<c_*+o(1) layer.
Treat q<1 by a multiplicity-compatible Hermite or degree chart.
Close bounded-L, L_epsilon, vertical-connector, core-to-main, and chart-join phase cells separately.
```

This artifact proves the saddle-centered coordinate identities, critical-frame formulas, u_a and v_a bounds, corrected-main and d_(n,x) budgets, core-scalar approximation, explicit endpoint defect, all-crossing sign unification, and successor-count substitution. It does not prove the q>=1 scalar lower bound or signed count, the q<1 multiplicity-compatible theorem, finite phase cells, one-sided successor winding, contact exclusion, Lambda<=0, RH, PF-infinity, or a Clay-prize conclusion.
