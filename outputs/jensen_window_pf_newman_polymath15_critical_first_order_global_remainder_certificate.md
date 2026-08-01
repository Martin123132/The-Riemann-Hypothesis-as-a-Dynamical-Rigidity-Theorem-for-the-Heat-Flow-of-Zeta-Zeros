# Newman Polymath-15 Critical First-Order Global Remainder Certificate

Date: 2026-07-26

Status: global first-order endpoint and cutoff remainder on the
`L>=50`, `0<tL<=25` critical ray. This is not a proof of contact
exclusion, `Lambda<=0`, RH, or a Clay-prize conclusion.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_critical_first_order_global_remainder_certificate.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_critical_first_order_global_remainder_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_critical_first_order_global_remainder_certificate.py
```

## Exact First-Order Cancellation

[Polymath 15](https://arxiv.org/pdf/1904.12438), Propositions 6.2 and 6.3,
give the endpoint heat integral and its coefficient/remainder
bounds. The extracted coefficient is

```text
C_1(p,sigma)=F'''(p)/(12*pi^2)+i*(2*sigma-1)*F'(p)/(4*pi), F=C_0
For sigma=(1-y)/2 and p'(T)=-1/(2*pi*a), (C_1(p,sigma)-C_1(p,1/2))/a=(i*y/2)*F'(p)*p'(T)=-i*y*F'(p)/(4*pi*a)
```

Thus the off-axis `F'` term is not discarded. It is exactly the
first complex displacement of `C0(p(T(z)))`. This cancellation is
what permits a second-order holomorphic collar.

## Fourth-Derivative Strip

The old removable-disc argument is enlarged by a radius-`1/10`
Cauchy disk:

```text
On |Re p|<=101/100, |Im p|<=1/100, Cauchy radius 1/10 stays in |Re p|<=111/100, |Im p|<=11/100; |F|<5 gives |F'|<50, |F''|<1000, |F'''|<30000, |F''''|<1200000
direct expanded C0 bound = [4.545006240074531751161629560727004915340271656380749415359425070964983 +/- 4.82e-70] < 5
p displacement = [4.420669830983571318809063394042258698086734524135207993367751675303341e-14 +/- 2.77e-84] < 1/100
```

## Heat Prefactor

Writing the residual multiplier as `P(u)`, the `C0` and retained
`C1` pieces satisfy

```text
|C_0|*integral|P-1|dmu<=2/T
|C_1(p,u)|<=300+10|u| and a^-1*integral|P-1||C_1|dmu<=1000/T
scaled C0 constant = [3.333364444716298570224192436663702224903351057681065122711056240564068 +/- 3.23e-71] < 4
tilted C1 moment = [418.3887741438575101841332060857839436867705243105808336548409006595615 +/- 3.19e-68] < 500
```

This includes the gamma/Stirling and `log M0` defects from the
proof of Proposition 6.3.

## Positive Shift

```text
For u>=0 choose K=2; the C_2/a^2 and RS_2 terms contribute less than 9/T after Gaussian integration
positive a^-2 coefficient = [0.6537465017405402999969178080450667596372026498570120002237305091122994 +/- 2.97e-71] < 1
```

## Negative Shift

A fixed `K=2` is invalid when `u<0`. Instead use

```text
For u<0 choose K_u=max(floor(-u)+3,floor(T_0/pi)); the k>=2 low series is split into even/odd Gamma progressions, and the published extreme tail delta_3<=2*10^-30/a^14; the total negative contribution is less than 1/T
negative low series = [0.01226709802781813630862695378997646815059986111588753861772092344133046 +/- 3.75e-72] < 1/50
combined T-scaled tail = [4.189457731676133389787674308186841860067643027501996545237719259736133 +/- 4.89e-70] < 10
```

The extreme range is the same Fubini-Tonelli tail bounded by
`2*10^-30/a^14` in the published proof.

## Holomorphic Lift

```text
K(T)=-(sqrt(pi)/8)*exp(-pi*T/4)*(T^(3/2)+i*T^(1/2)); (T^(3/2)+i*T^(1/2))/a=sqrt(2*pi)*(T+i)
E_[1],N(z)=(-1)^N exp(t*pi^2/64)*(G_[1],N(z)+G_[1],N^#(z)), where G_[1],N=K(T(z))*(F(p_N(T(z)))+F'''(p_N(T(z)))/(12*pi^2*a(T(z))))
For A(T)=T^(3/2)+i*T^(1/2), |theta|<=1/L and T>=10^6: |A'|/|A(T)|<=2/T, |A''|/|A(T)|<=1/T^2, |p'|<=(2*pi*(T-|theta|))^-1/2, |p''|<=(2*sqrt(2*pi)*(T-|theta|)^(3/2))^-1; therefore T|S_0''|/|A(T)|<2100, T|S_1'|/|A(T)|<10133, and the two-component paper/lift mismatch is <1000/T
two-component lift constant = [406.5509466460494330466387380103983765396614506777735504436333361948041 +/- 1.39e-68] < 1000
```

Combining the heat expansion and lift comparison gives

```text
The two heat endpoints, after the C_1 lift is retained, have bracket error <=4000/T on the doubled transition collar
The published endpoint prefactor, the collar loss, and sup|B_t(z)|/A_t(x)<2 give |K(T)|/A_t(x)<2*exp(7/4)*exp(-L/4)<20*exp(-L/4); also T>pi*exp(L)
|R_endpoint|/A_t(x)<30000*exp(-5L/4)
```

## Fixed Cutoffs

The finite Dirichlet second-order theorem is lower order:

```text
|R_D|/A_t(x)<8000000*exp(-7L/4), |partial_x R_D|/A_t(x)<8000000*L*exp(-7L/4)
On a fixed-N disk, sup |R_[1]|/A_t(x)<40000*exp(-5L/4)
L=50 absorption = [0.0001111035509197121647572941099686948552831980830416404446683823996870509 +/- 4.52e-74] < 1
```

## Adjacent Cutoffs

The inherited `C0` mismatch, entering corrected Dirichlet term,
and parity-matched endpoint `C1` term give

```text
The inherited C_0/main-block adjacent mismatch is <1000*exp(-5L/4)
|d_(t,n)|<100/T for the entering Dirichlet term; critical C_1 parity gives exact endpoint value matching and its adjacent derivative difference is O(T^-1); together with the old mismatch, |Delta lift_[1]|/A_t(x)<20000*exp(-5L/4)
From |alpha'(s)|<=1/(T-3), t|alpha_n|<=27, d=1/(6s)+alpha'(s)(t/4+t^2*alpha_n^2/8) gives T|d|<100. The entering two-sharp block is <200*exp(-L/4), costing <7000*exp(-5L/4). C_1 parity makes the endpoint difference zero at a=m; the F'''' derivative bound costs <6000*exp(-5L/4) on the doubled collar
saved adjacent sum = [10335.90604091243822572862489824958285031082414467031985879389928188035 +/- 2.70e-66] < 20000
```

## Global First Jet

One fixed analytic lift is now available on every Cauchy disk.
Cauchy's estimate and the real normalizer derivative yield

```text
For every critical disk, including cutoff crossings, |r_[1](x)|<100000*exp(-5L/4) and |partial_x r_[1](x)|<200000*L*exp(-5L/4)
A radius-1/L Cauchy estimate gives |R_[1]'|/A_t<100000*L*exp(-5L/4); with |partial_x log A_t|<L/2 this gives |partial_x r_[1]|<150000*L*exp(-5L/4)<200000*L*exp(-5L/4)
```

## Remaining Theorem

The approximation side no longer sets the critical-ray scale.
The next obligation is

```text
Insert eta_0=100000 and eta_1=200000 in the exact signed contact-normal inequality; the strict Xi arithmetic reversal is the next open theorem target
```

That signed arithmetic inequality, the bounded-`L` shoulder, and
the final exhaustion remain open.
