# Newman Polymath-15 Critical First-Order Signed-Contact Reduction

Date: 2026-07-26

Status: exact first-order contact reduction with certified error
constants. This is not a proof of contact exclusion, `Lambda <= 0`,
or RH.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_critical_first_order_signed_contact_reduction.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_critical_first_order_signed_contact_reduction.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_critical_first_order_signed_contact_reduction.py
```

## Corrected Complex Main

On each radius-`1/L` proof disk, select one prescribed-`N`
analytic lift. The adjacent-lift discrepancy is already part of
the certified residual.

```text
s=(1-i*x)/2, s_*=s+t*alpha(s)/2, phi=M_t(s)/|M_t(s)|, e_n=phi*exp(t*log(n)^2/4-s_*log(n))
alpha_n=alpha(s)-log(n), d_n=1/(6s)+alpha'(s)*(t/4+t^2*alpha_n^2/8), f_n=e_n*(1+d_n)
d_(n,x)=(-i/2)*[-1/(6s^2)+alpha''(s)*(t/4+t^2*alpha_n^2/8)+(t^2/4)*alpha_n*alpha'(s)^2]
T_0=x/2+pi*t/8, a=sqrt(T_0/(2*pi)), p=1-2*(a-N), F=C_0, kappa_N=(-1)^N*exp(t*pi^2/64)*M_0(iT_0)*U*exp(pi*i/8)/|M_t(s)|, g_0=-kappa_N*(F(p)+F'''(p)/(12*pi^2*a))
E_[1]=g_0+sum_(n=1)^N f_n, J_[1]=2*Re(E_[1])
```

The previous peel can now be absorbed without changing the zero
system:

```text
J_[1]=J_[0]+DeltaJ, r_[0]=DeltaJ+r_[1], so J_[0]+r_[0]=J_[1]+r_[1] and the value/derivative contact systems are identical
```

## Rate-Free Contact Identity

Use real component values and derivatives, so isolated component
zeros do not create singular logarithmic rates.

```text
For components z_j in {g_0,f_1,...,f_N}, put c_j=Re(z_j), d_j=Re(z_(j,x)); P={c_j>0}, N_-={c_j<0}, Z={c_j=0}; M_+=sum_P c_j, M_-=sum_N_(-c_j), h_+=sum_P d_j/M_+, h_-=sum_N_(-d_j)/M_-, D_0=sum_Z d_j, with an empty mean set to 0
X_[1]=Re(E_[1])=M_+-M_-, U_[1]=Re(E_[1],x)=M*(h_+-h_-)+(X_[1]/2)*(h_++h_-)+D_0, M=(M_++M_-)/2
At H=H_x=0, X_[1]=-r_[1]/2 and U_[1]=-r_[1],x/2; hence M*(h_+-h_-)+D_0=-r_[1],x/2+(r_[1]/4)*(h_++h_-)
```

## Certified Contact Box

The global first-order theorem supplies

```text
|X_[1]|<50000*exp(-5L/4), |U_[1]|<100000*L*exp(-5L/4)
```

Therefore the exact open frequency-layer theorem is

```text
|X_[1]|<=50000*exp(-5L/4) implies |U_[1]|>100000*L*exp(-5L/4)
```

A stronger partition form is

```text
Whenever |X_[1]|<=50000*exp(-5L/4), prove |M*(h_+-h_-)+D_0|>25000*exp(-5L/4)*(4L+|h_++h_-|)
```

Neither inequality is proved here.

## Saddle-Centered Observable

Centering at the moving Riemann-Siegel saddle gives

```text
lambda_a=phi'/phi-s_*'*log(a)=u_a+i*v_a; R_a=-s_*'*sum_(n=1)^N log(n/a)f_n+sum_(n=1)^N e_n*d_(n,x)+(g_(0,x)-lambda_a*g_0); E_[1],x=lambda_a*E_[1]+R_a
At contact, Re(R_a)-v_a*Im(E_[1])=-(r_[1],x-u_a*r_[1])/2, so |Re(R_a)-v_a*Im(E_[1])|<(100000L+50000|u_a|)*exp(-5L/4)
```

This is the proof-facing arithmetic object: a centered logarithmic
moment, the explicit `d_n` derivative, and one endpoint defect.

## Domain Split

```text
L=log(x/(4*pi))>=50, 0<tL<=25, q=2*t*L^2; the live frequency layer has q>=1 and 0<tL<=4911678521/1933561194+o(1)
For every fixed epsilon>0, tL>=4911678521/1933561194+epsilon is already closed for sufficiently large L by the oscillatory-zeta theorem
Prove the signed band only on q>=1 in the live asymptotic frequency layer; use a multiplicity-compatible parabolic/Hermite or boundary-degree theorem on q<1; terminate the bounded-L shoulder
```

The restriction to `q>=1` is logical, not cosmetic:

```text
H_tau(y)=y^2-2tau has no positive-time contact but its root slopes are 2*sqrt(2tau)->0; a fixed positive slope floor uniform to tau=0 would impose endpoint simplicity
```

A slope floor uniform to `t=0` would demand endpoint simplicity,
which RH does not require. The ultra-small parabolic layer therefore
needs a multiplicity-compatible Hermite or boundary-degree argument.

## Proof Boundary

This artifact proves the corrected component representation, old/new repartition equivalence, rate-free near-crossing identity, explicit contact box, and saddle-centered observable. It does not prove the strict Xi crossing-band inequality, the q<1 multiplicity-compatible theorem, the bounded-L shoulder, contact exclusion, Lambda<=0, RH, PF-infinity, or a Clay-prize result.
