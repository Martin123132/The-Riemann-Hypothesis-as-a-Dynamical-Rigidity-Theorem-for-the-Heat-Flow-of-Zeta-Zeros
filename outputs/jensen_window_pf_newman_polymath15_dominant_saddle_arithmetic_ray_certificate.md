# Jensen-Window PF Newman Polymath-15 Dominant-Saddle Arithmetic Ray Certificate

Date: 2026-07-16

Status: rigorous positive normalized curvature for the complete finite
Riemann-Siegel main sum on the dominant-saddle ray. This is not yet a
proof for `H_t`, `Lambda <= 0`, or RH.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_dominant_saddle_arithmetic_ray_certificate.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_dominant_saddle_arithmetic_ray_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_dominant_saddle_arithmetic_ray_certificate.py
```

The input finite sum and cutoff come from D. H. J. Polymath,
[Effective approximation of heat flow evolution of the Riemann xi function](https://arxiv.org/abs/1904.12438),
Theorem 1.3.

## Ray

```text
0<t<=1/2, L=log(x/(4*pi)), t*L>=24; hence L>=48 and x=4*pi*exp(L)>10^21
```

The implication `L>=48` places the entire argument beyond
`x=[8817490397450238834356.185569219276955454895729125369982449075043686890 +/- 4.51e-49]`, rigorously larger than `10^21`.
The large threshold is a theorem-building choice, not a numerical
claim that smaller points fail.

## Coefficients

On this ray, elementary estimates for the published normalizer give

```text
Re(alpha)=L/2+(1/4)log(1+x^-2)-1/(1+x^2)
|alpha|<=L/2+1, |alpha'|<=2/x, |alpha''|<=3/x^2
g'=alpha+(t/2)alpha*alpha', L/2-x^-2-(L/4+1/2)/x<=Re(g')<=L/2+(4x^2)^-1+(L/4+1/2)/x; |g''|<=2/x+(3L/8+7/4)/x^2<=4/x
0.24L<=b=|beta'|<=0.26L, a=|s_*'|<=0.51, c=|beta''|<=10^-21, d=|s_*''|<=10^-42, |V|<=10^-21
```

Here `g=log(M_t)`. The first line is exact. The modulus bounds use
`|s|=|s-1|=sqrt(1+x^2)/2`; substituting them into
`g'=alpha+(t/2)alpha*alpha'` and
`g''=alpha'+(t/2)(alpha'^2+alpha*alpha'')` gives the displayed
normalizer bounds. Since `L>=48` and `x>10^21`, they imply the
stated decimal rational envelopes for `b,a,c,d,V`.

The published lower bound for `Re(s_*)`, together with the saddle
cutoff, then proves

```text
p_n>=1/2+(t/4)log((x/(4*pi))/n)-t/(2x^2)>=7/2-[1/(512*exp(L))+1/(4x^2)]>349/100
|exp((t/4)log(n)^2)*n^(-s_*)|<=n^(-3.49), 2<=n<=N
```

The exponent loss from replacing the exact power by `3.49` is
enclosed by `[2.783523599103388879466733144014579482500946841747077570921460357700027e-24 +/- 4.14e-94]`, below `1/100`.

## Arb Moments

For `p=349/100`, outward-rounded Arb series evaluation gives

```text
S0=sum_(n>=2)n^-p       = [0.1278709998316298986307925385245536726455110461461403941676865286585381 +/- 1.49e-72] < 16/125
S1=sum_(n>=2)log(n)n^-p = [0.1143098929764890452589236997639857299254191266796851536249812640499371 +/- 9.21e-72] < 23/200
S2=sum_(n>=2)log(n)^2n^-p = [0.1198423145124326358478067088669221942024306397216094225232479589965527 +/- 4.88e-71] < 121/1000
```

Exact logarithmic differentiation therefore yields

```text
Q0=2*S0; Q1=2*(b*S0+a*S1); Q2=2*((b^2+c)*S0+(2*b*a+d)*S1+a^2*S2)
Q0 <= 0.256
Q1 <= 0.06901*L
Q2 <= 0.01861*L^2
```

These bounds use the infinite zeta tails, so they cover every finite
cutoff without computation at the enormous endpoint scale.

## Curvature

The leading phase and the exact perturbation identity give

```text
C[P0]>=4b^2-2c-4|V|>0.23L^2, P0=2cos(beta)
Err<=4bQ1+Q1^2+2Q2+Q0(2c+2b^2)+Q0Q2+|V|(4Q0+Q0^2)
```

The final budget is exact rational arithmetic:

```text
phase floor / L^2 > 23/100
tail error / L^2  < 77/500
strict margin / L^2 = 19/250
```

Hence, with `N` held fixed when differentiating its finite analytic
formula,

```text
C_t[P_(N,t)]>0.076L^2 for every point on the ray, with N held fixed when differentiating the finite analytic main sum
```

## Remaining Gap

This closes the arithmetic main sum, not the exact heat-flow function.
The next job is to put the published remainder on a fixed-`N` complex
collar, transfer it through two derivatives, and show its normalized
curvature error is below `0.076L^2`. Adjacent-`N` overlap is still
needed at the discrete cutoff transitions. The complementary region
`0<tL<24` remains the coupled near-zero-time boundary layer.
