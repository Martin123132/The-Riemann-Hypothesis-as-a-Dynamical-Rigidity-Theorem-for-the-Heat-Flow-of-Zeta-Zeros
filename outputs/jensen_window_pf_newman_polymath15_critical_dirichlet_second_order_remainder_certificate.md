# Newman Critical Dirichlet Second-Order Remainder Certificate

Date: 2026-07-26

Status: explicit fixed-cutoff finite-Dirichlet theorem; not a proof
of contact exclusion, `Lambda <= 0`, or RH.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_critical_dirichlet_second_order_remainder_certificate.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_critical_dirichlet_second_order_remainder_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_critical_dirichlet_second_order_remainder_certificate.py
```

## Uniform Shift

```text
L>=50, 0<t<=1/2, tL<=25, T>=2*pi*exp(L-1), |sigma|<=1, n<=N, and the radius-1/L disk lies in one prescribed-N cell
On the critical collar, |alpha(s)|<=L/2+(1+pi)/2+2/T and log(n)<=L/2+1; hence t|alpha_n|<=27
w=sqrt(t)v+c, c=t*alpha_n/2, |c|<=27/2
```

The actual lower bound for `T` at `L=50` is
`11984211960000554356237.4704588507437795648004701854030141285`. The proof below only
uses `T>=10^6`.

## Central Range

```text
V=T^(1/3); central: |v|<=V; tail: |v|>V
For T>=10^6, |v|<=V implies |w|<=V and Im(s+u*w)>=T/2 for 0<=u<=1
|alpha'(s)|<=1/T and sup_(0<=u<=1)|alpha''(s+u*w)|<=3/T^2
g=1/(6s)+alpha'(s)w^2/2; rho=R_M+R_G+1/(6(s+w))-1/(6s); |rho|<=(|w|^3/2+|w|/3+12/T)/T^2
|g+rho|<0.006, |exp(g+rho)-1-g|<=|rho|+|g+rho|^2
E|w|^3<10980, E|w|<15, E|w|^4<307330; E|rho|<5496/T^2, E|g|^2<153666/T^2, E|rho|^2<T^-2
integral_(|v|<=V)|exp(g+rho)-1-g|dmu <=312830/T^2
```

The pointwise exponent is below
`0.00500116666666666666666666666666666666666666666666666666666667` at the
artificial audit endpoint and decreases thereafter.

## Gaussian Tail

The proof of [Polymath 15 Proposition 6.1](https://arxiv.org/abs/1904.12438)
supplies the quadratic envelope

```text
The pointwise envelope in the proof of Polymath Proposition 6.1 gives |I(v)|dmu<=pi^(-1/2)exp(92/T)*exp(-(1-0.26/T)v^2)dv
For T>=10^6: integral_tail |I|dmu<=T^-2, mu(tail)<=T^-3, E[v^2 1_tail]<=T^-1, and integral_tail |I-1-g|dmu<=2/T^2
```

The three logarithmic endpoint margins are respectively
`-9877.53637167713074180570825027127433957372711751249529615204`,
`-9963.73100345501996914278385040072670203420602358839728816922`, and
`-9980.88617207548768921843883624106690977111568795660956093761`; their
checked derivatives are negative. Their proof role is supplied
by the exact rational guards saved alongside this audit:

```text
At T=10^6, T^(1/3)=100 and T^(2/3)=10000. The exact rational bounds log(T)<14, log(2)<1, pi>157/50, and e<68/25 make all three tail logarithms negative, their derivative upper bounds negative, the central exponent <3/500, and the fixed-cell constant <8000000
```

## Per-Term Theorem

```text
|r_(t,n)(s)/m_(t,n)(s)-(1+d_(t,n)(s))|<=400000/T^2
```

The central and tail constants sum to less than `400000`; no
asymptotic `O` constant is hidden in this statement.

## Fixed-Cell Transfer

```text
Using the two-block coefficient mass <=50*exp(L/4), sup_collar |B_t|/A_t<2, and T>=2*pi*exp(L-1): |R_D|/A_t<8000000*exp(-7L/4)
Cauchy on the radius-1/L disk gives |partial_x R_D|/A_t<8000000*L*exp(-7L/4)
8000000*exp(-7L/4)<0.000112*exp(-5L/4) for L>=50
```

Thus the two finite Dirichlet blocks no longer determine the
global error exponent after their first signed corrections are
retained.

## Remaining Boundary

```text
The fixed-N Dirichlet residual is closed. A global corrected remainder still needs the heat-integrated endpoint a^-2 bound and one adjacent-cutoff comparison for the new signed lift.
```

The next quantitative theorem is the heat-integrated endpoint
`a^-2` bound, followed by its adjacent-cutoff lift. Only after
those steps can the signed contact-normal inequality be tested.
