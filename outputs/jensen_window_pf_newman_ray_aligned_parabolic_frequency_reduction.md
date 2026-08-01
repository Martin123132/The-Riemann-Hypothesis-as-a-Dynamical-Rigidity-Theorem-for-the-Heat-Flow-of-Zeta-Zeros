# Newman Ray-Aligned Parabolic-Frequency Reduction

Date: 2026-07-26

Status: exact exhaustion redesign, exact positive scale, and conditional
one-antecedent cofinal reduction. The Xi descendant estimate remains open;
this is not a proof of `Lambda<=0` or RH.

## Why Change The Exhaustion

The linear-radius schedule used previously,

```text
t_j=1/(5j), R_j=j+38,
```

forces each new right strip into `tL->0`, even though the existing
high-frequency theorem already proves strict first-Laguerre positivity on
`L>=50`, `tL>=25`. The diagonal-exhaustion theorem permits any time floors
tending to zero and any radii tending to infinity. They need not be coupled
linearly.

Choose instead

```text
c=25, a=100,
t_j=c/(a+j),
L_j=log(R_j/(4pi))=a+j+1,
R_j=4pi exp(L_j).                                  (1)
```

The base time is `t_0=1/4>1/5`, so the published `Lambda<=1/5` bound makes
the base rectangle `[1/4,1/2]x[0,R_0]` contact-free.

## Every New Strip Is Already Closed

Let

```text
P_j=[t_j,1/2]x[0,R_j],
C_j^out=[t_(j+1),t_j]x[38,R_j],
S_j=[t_(j+1),1/2]x[R_j,R_(j+1)].
```

For every point of `S_j`,

```text
tL(x)>=t_(j+1)L_j
      =[25/(a+j+1)](a+j+1)=25.                    (2)
```

Also `L(x)>=L_j>=101>50`. The checked dominant-saddle global-ray theorem
therefore gives

```text
H_x^2-H H_xx>0                                    (3)
```

throughout every new strip. Hence a new strip cannot contain a first-jet
contact. This removes the formerly open all-stage right-strip cone or
weighted slope-gap theorem from the induction antecedents.

The schedule is cofinal:

```text
t_j->0, R_j->infinity.                             (4)
```

The independently certified compact theorem handles `0<=x<=38` for all
relevant times. What remains is exactly the outer old descendant collar
`C_j^out`.

## Parabolic-Frequency Scale

For `L=log(x/(4pi))>0`, define

```text
s_pf(t,x)=(L^2+(2t)^(-1))^(-1/2)
         =sqrt(2t)/sqrt(1+2tL^2),                  (5)
V_pf=(H,s_pf H_x).
```

Put `q=2tL^2` and compare (5) with the parabolic scale
`s_par=sqrt(2t)` and frequency scale `s_freq=1/L`:

```text
alpha=s_pf/s_par=1/sqrt(1+q),
beta=s_pf/s_freq=sqrt(q/(1+q)),
alpha^2+beta^2=1.                                 (6)
```

Consequently

```text
min(s_par,s_freq)/sqrt(2)<=s_pf<=min(s_par,s_freq). (7)
```

The transition `q=1` has `alpha=beta=1/sqrt(2)`. There is no badly
conditioned chart interface.

The exact scale derivatives are

```text
partial_t log(s_pf)=1/[2t(1+q)],
partial_t log(alpha)=-L^2/(1+q),
partial_t log(beta)=1/[2t(1+q)],
partial_x log(s_pf)=-2tL/[x(1+q)].                 (8)
```

At entry, (1)-(2) give

```text
q_entry=2t_(j+1)L_j^2=50L_j>=5050,                (9)
```

so the scale is almost exactly `1/L`. For any fixed descendant coordinate
`x`, `q=2tL(x)^2` tends to zero with `t`, and the same scale becomes almost
exactly `sqrt(2t)`, the multiplicity-compatible Hermite coordinate.

## Relative Two-Chart Handoff

Let `V_par=(H,s_par H_x)` and `V_freq=(H,s_freq H_x)`. Since
`V_pf=diag(1,alpha)V_par=diag(1,beta)V_freq`, the exact rescaling inequality
gives

```text
q<=1:
  kappa_pf<=sqrt(2) kappa_par+L^2/(1+q);

q>=1:
  kappa_pf<=sqrt(2) kappa_freq+1/[2t(1+q)].        (10)
```

Here `kappa_par` or `kappa_freq` denotes a proved relative bound for the
corresponding Xi jet. Equation (10) does not assert either bound. It shows
that once such bounds are established in their natural regimes, the smooth
interface costs at most `sqrt(2)`.

The schedule also has

```text
log(t_j/t_(j+1))
 =log((a+j+1)/(a+j))->0,                           (11)
```

so a Hermite-type `K/t` estimate has vanishing one-step logarithmic cost.

## Conditional Cofinal Theorem

Suppose that on every outer old collar `C_j^out` there is an integrable
`kappa_j(t)` such that

```text
||partial_t V_pf(t,x)||<=kappa_j(t)||V_pf(t,x)||.  (12)
```

Gronwall transports every old nonzero jet through `C_j^out`. Equation (3)
certifies `S_j`, and the positive scale preserves contact sets and degree.
Induction then certifies every `P_j`; cofinality (4) gives positive-time
simplicity and therefore `Lambda<=0`.

This composition is exact, but (12) for Xi is not proved.

The rectangular condition has an equivalent direct wedge formulation. Put

```text
tau(x)=min(1/4,25/L(x)),             x>=38.         (13)
```

At `t=tau(x)`, the first jet is nonzero: when `L(x)<=100`, this follows
from `1/4>Lambda`; when `L(x)>=100`, it follows from `tL=25` and the
dominant-saddle theorem. Therefore it is sufficient to prove, for every
fixed `x>=38`,

```text
||partial_t V_pf(t,x)||<=kappa_x(t)||V_pf(t,x)||,
0<t<=tau(x),    integral_delta^tau(x) kappa_x<infinity
for every delta>0.                                  (14)
```

Gronwall then excludes a contact at every positive `t`. The exponential
rectangles (1) are a finite-stage exhaustion of this one wedge target.

There is an essential nonpromotion guard. For any `C^1` vector path `V` on
`[delta,tau]`, a nonzero endpoint together with an integrable relative bound
implies nonvanishing by Gronwall. Conversely, if `V` is already nonvanishing,
then

```text
kappa(t)=||partial_t V(t)||/||V(t)||
```

is continuous and supplies such a bound. Thus bare existence of `kappa_x`
is equivalent to the desired positive-time noncontact statement. A useful
Xi theorem must construct an a priori majorant from independently bounded
arithmetic or analytic quantities, without dividing by the unknown jet norm.
The scale `s_pf` conditions the two asymptotic charts; it does not itself
provide that majorant.

## Explicit Schedule Diagnostics

| stage | t_j | L_j | t_(j+1) | min tL | entry q | beta^2 | unit-K log cost |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | 1/4 | 101 | 25/101 | 25 | 5050 | 0.99980201940209856 | 0.009950330853168092 |
| 1 | 25/101 | 102 | 25/102 | 25 | 5100 | 0.99980396000784155 | 0.0098522964430116395 |
| 10 | 5/22 | 111 | 25/111 | 25 | 5550 | 0.99981985227886871 | 0.0090498355199178562 |
| 100 | 1/8 | 201 | 25/201 | 25 | 10050 | 0.99990050741219783 | 0.0049875415110389679 |
| 1000 | 1/44 | 1101 | 25/1101 | 25 | 55050 | 0.99998183502570348 | 0.00090867793621811225 |
| 10000 | 1/404 | 10101 | 25/10101 | 25 | 505050 | 0.99999802000194038 | 9.900499983335809e-05 |

## Sharper Asymptotic Variant

For any fixed `epsilon>0`, the same construction may replace `25` by
`c_*+epsilon`, where

```text
c_*=4911678521/1933561194.
```

The oscillatory-zeta theorem then closes every new strip once its
existential `L_epsilon` threshold is passed. The `c=25`, `L>=50`
construction above is preferred as the fully explicit reduction.

## Proof Boundary

The exhaustion identities, new-strip composition, blended scale,
two-chart conditioning, derivative formulas, and conditional induction are
rigorous. The result does not prove the relative Xi bound (12), contact-free
old collars at every stage, `Lambda<=0`, RH, PF-infinity, or the Clay prize.

Machine-audited files:

```text
work/rh_compute/results/jensen_window_pf_newman_ray_aligned_parabolic_frequency_reduction.json
work/rh_compute/scripts/jensen_window_pf_newman_ray_aligned_parabolic_frequency_reduction.py
work/rh_compute/scripts/check_jensen_window_pf_newman_ray_aligned_parabolic_frequency_reduction.py
```
