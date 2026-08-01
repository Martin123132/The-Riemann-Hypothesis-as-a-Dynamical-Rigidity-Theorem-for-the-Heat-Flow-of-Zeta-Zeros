# Jensen-Window PF Newman Polymath-15 Critical Dirichlet First Correction Gate

Date: 2026-07-26

Status: exact first signed finite-sum correction and guarded
second-order remainder programme. This is not a proof of
`Lambda <= 0` or RH.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_critical_dirichlet_first_correction_gate.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_critical_dirichlet_first_correction_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_critical_dirichlet_first_correction_gate.py
```

## Exact Centered Ratio

[Polymath 15](https://arxiv.org/abs/1904.12438), equation (40), gives

```text
alpha_n=alpha(s)-log(n), dmu(v)=pi^(-1/2)exp(-v^2)dv
r_(t,n)(s)=exp(-t*alpha_n^2/4)*integral exp(-sqrt(t)*v*alpha_n)r_(0,n)(s+w(v))dmu(v), w(v)=sqrt(t)*v+t*alpha_n/2
r_(0,n)(s)=M_0(s)*n^(-s)*GammaStar(s/2), GammaStar(z)=Gamma(z)/(exp(-z)z^z*sqrt(2*pi/z))
```

Substitution and exact cancellation yield

```text
m_(t,n)(s)=M_0(s)*n^(-s)*exp(t*alpha_n^2/4)
-t*alpha_n^2/4-sqrt(t)*v*alpha_n+alpha_n*w=t*alpha_n^2/4
r_(t,n)(s)/m_(t,n)(s)=integral GammaStar((s+w)/2)*exp(log(M_0(s+w))-log(M_0(s))-alpha(s)w)dmu(v)
```

This ratio contains every term discarded by Proposition 6.1.

## First Signed Term

The centered Gaussian moment is

```text
integral w(v)^2 dmu(v)=t/2+t^2*alpha_n^2/4
```

The scaled-gamma coefficient `1/(6s)` and the quadratic
`log M_0` term therefore give

```text
d_(t,n)(s)=1/(6s)+alpha'(s)*(t/4+t^2*alpha_n^2/8)
m_[1](t,n;s)=m_(t,n)(s)*(1+d_(t,n)(s))
```

This is signed and componentwise. It is the finite-sum analogue
of retaining the first omitted Riemann-Siegel endpoint coefficient.

## Remainder Anatomy

The [DLMF log-gamma expansion](https://dlmf.nist.gov/5.11) gives

```text
log GammaStar(z)=1/(12z)+R_G(z), |R_G(z)|<=sec(arg(z)/2)^4/(360*|z|^3)
For z=(s+w)/2 with Im(s+w)>0, |R_G(z)|<=16*|s+w|/(45*Im(s+w)^4)
alpha'(s)=-1/(2s^2)-1/(s-1)^2+1/(2s); alpha''(s)=1/s^3+2/(s-1)^3-1/(2s^2)
R_M=log(M_0(s+w))-log(M_0(s))-alpha(s)w-alpha'(s)w^2/2, |R_M|<=|w|^3*sup_(0<=u<=1)|alpha''(s+u*w)|/6
rho(v)=R_M+R_G+1/(6(s+w))-1/(6s)
```

A global cubic Taylor bound cannot simply be exponentiated under
the whole Gaussian. The correct proof split is

```text
Split |v|<=V and |v|>V with V=T^(1/3), T=Im(s): use the second-order expansions on the central part and the published Proposition 6.1 quadratic Gaussian envelope on the tail
```

The central range uses the signed second-order expansion; the tail
reuses the published quadratic Gaussian envelope.

## Scale Forecast

The theorem to certify per component is

```text
|r_(t,n)/m_(t,n)-(1+d_(t,n))|<=C_D/T^2 uniformly for n<=N on the critical radius-1/L collar
sum of the two absolute finite-main coefficient masses <=50*exp(L/4), T>=2*pi*exp(L-1)
C_D/T^2 times 50*exp(L/4)=O(C_D*exp(-7L/4))
```

Thus the peeled Dirichlet sums should be smaller than the new
endpoint remainder:

```text
The peeled Dirichlet residual is O(exp(-7L/4)); the RS-C_1 endpoint and cutoff residuals are O(exp(-5L/4)) and therefore set the next global C^1 scale
```

The exponent comparison is conditional until an explicit `C_D`
is proved.

## Signed Contact Handoff

```text
DeltaJ=(DeltaA+DeltaB-DeltaC)/A_t, J_[1]=J_[0]+DeltaJ, r_[0]=DeltaJ+r_[1]
Prove explicit constants eta_0,eta_1 with |r_[1]|<=eta_0*exp(-5L/4), |r_[1],x|<=eta_1*L*exp(-5L/4) on L>=50, 0<tL<=25, including adjacent cutoffs, then use r_0=DeltaJ in the signed contact-normal inequality
```

This is the next quantitative job: certify the central Gaussian
moments, tail, first derivative, and cutoff collars with explicit
constants. Only the residual then pays an absolute-value budget.
