# Jensen-Window PF Newman Strict-Laguerre Monotonicity Scout

Date: 2026-07-17

Status: exact sufficient condition and finite scouts, retired by a rigorous Xi Lehmer-point counterexample. This is not a proof of RH or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_newman_strict_laguerre_monotonicity_scout.json
python work/rh_compute/scripts/jensen_window_pf_newman_strict_laguerre_monotonicity_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_strict_laguerre_monotonicity_scout.py
```

Current result:

```text
validated Jensen-window PF Newman strict-Laguerre monotonicity scout: 9 rows, 0 issues, 5 exact target identities, 2 exact classical-route obstructions, 6 dense time rows, 20 high-frequency rows, 1 theta-tail nonpromotion guard, 1 Arb Xi monotonicity rejection
```

## Exact Sufficient Candidate

Define

```text
L_t(x)=H_t'(x)^2-H_t(x)*H_t''(x)
M_t(x)=-L_t'(x)=H_t(x)*H_t'''(x)-H_t'(x)*H_t''(x)
```

Then

```text
If M_t(x)>0 for every x>0 and L_t(x)->0 as x->infinity, then L_t(x)=integral_x^infinity M_t(y)dy>0 for every x>=0.
M_t(x)>0 for every x>0 and every 0<t<=1/5 implies Lambda<=0.
```

This condition is stronger than strict Laguerre positivity. The implication is
exact, but the condition itself is false for Xi, as certified below.

## Correlation And Theta Forms

```text
M_t(x)=2*integral_R v*K_(1,t)(v)*sin(2*x*v)dv=4*integral_0^infinity v*K_(1,t)(v)*sin(2*x*v)dv
M_t(x)=((A_t(x)+1/2)*A_t'''(x)-A_t'(x)*A_t''(x))/64
```

The candidate analytic problem was positivity of a sine transform of
`v*K_(1,t)(v)`, or equivalently the displayed third-order theta-primitive
Wronskian.

Writing the sine kernel as an integrated cosine gives a second exact form:

```text
Q_t(s)=integral_s^infinity v*K_(1,t)(v)dv
M_t(x)=8*x*integral_0^infinity Q_t(s)*cos(2*x*s)ds
Q_t'(s)=-s*K_(1,t)(s), Q_t''(s)=-K_(1,t)(s)-s*K_(1,t)'(s)
```

For ell=log K_(1,t), h(s)=1+s*ell'(s) decreases strictly from 1 to -infinity because h'=ell'+s*ell''<0 and strong log-concavity gives ell'(s)<=-c*s. Hence Q_t''=-K_(1,t)*h changes sign exactly once, from negative to positive.
In particular, `Q_t''(0)=-K_(1,t)(0)<0, so the classical positive-decreasing-convex cosine-transform criterion cannot establish M_t>0.`
Thus the classical convex-kernel positivity theorem does not close this target.
Moreover, At every regular fixed shift p>0, evenness and real analyticity give Q_t(i*p) real. Hence -Im Q_t(s+i*p) starts at zero and cannot be both nonnegative and nonincreasing unless it is identically zero, so Tuck's elementary interior-shift criterion cannot apply nontrivially.
A viable complex-contour proof needs a boundary correction or a frequency-adaptive modular shift; the theta series identifies abs(Im z)=pi/8 as its natural termwise-control boundary.
The relevant primary source is
https://doi.org/10.1017/S0004972700047511.

## Dense Moderate Scout

| t | x rows | largest L_t' | at x | smallest L_t' | at x |
|---:|---:|---:|---:|---:|---:|
| 0 | 2000 | -4.11895579540490916e-12 | 40 | -3.92227242966496190e-06 | 6.9800000000000004 |
| 0.05 | 2000 | -4.10209881280951840e-12 | 40 | -3.93257078956562750e-06 | 6.9800000000000004 |
| 0.1 | 2000 | -4.08531227661249307e-12 | 40 | -3.94290383904955908e-06 | 6.9800000000000004 |
| 0.15 | 2000 | -4.06859688927927284e-12 | 40 | -3.95327171929068172e-06 | 6.9800000000000004 |
| 0.2 | 2000 | -4.05195335919097325e-12 | 40 | -3.96367457213198842e-06 | 6.9800000000000004 |
| 0.5 | 2000 | -3.95364154848317632e-12 | 40 | -4.02687440479802134e-06 | 6.96 |

The grid is `x=j/50`, `1<=j<=2000`. The maximum relative coarse/fine delta is
`1.74849836598498601e-10`. This is finite numerical
evidence, not continuum control.

## High-Frequency Scout

| t | x | L_t'(x) |
|---:|---:|---:|
| 0 | 40 | -0.00000000000411895579448632007353931142173811409819378510379368486431220 |
| 0 | 80 | -2.22290026634735121755276383842595196021310280009484941075692e-24 |
| 0 | 120 | -1.00634691469262228721287653498564011344956118045731557868619e-37 |
| 0 | 160 | -3.08599394685849317548546681814828971355035472342185972563336e-50 |
| 0 | 200 | -2.80325572279256566791574436050465122218566479313945104835900e-63 |
| 0.1 | 40 | -0.00000000000408531227568601710680446354610037081003233688207207262763809 |
| 0.1 | 80 | -2.25448648844558094273777661902921104780226057385782605875543e-24 |
| 0.1 | 120 | -1.05166087076946782841663807007500747979753539813664345885561e-37 |
| 0.1 | 160 | -3.33865621804258982945894909739915565636680745808989797800181e-50 |
| 0.1 | 200 | -2.95693562058735799011796257092637269516453798352928137855205e-63 |
| 0.2 | 40 | -0.00000000000405195335825762182855874424741942694559485048139509307942653 |
| 0.2 | 80 | -2.28693492558438958447517508980957774054812303590026709903936e-24 |
| 0.2 | 120 | -1.10273211172287721198649066493683243815671194037688303661425e-37 |
| 0.2 | 160 | -3.60653975990912179067078815377070050198252857511712950816567e-50 |
| 0.2 | 200 | -3.12046794559790362545710302003715582745416360579597774896448e-63 |
| 0.5 | 40 | -0.00000000000395364154752945574189771550160078020293911202566157252451126 |
| 0.5 | 80 | -2.38939093337290294731578567667251763100712112554904703997572e-24 |
| 0.5 | 120 | -1.29451695714042784331138514857963727849703957707716825872298e-37 |
| 0.5 | 160 | -4.50719488796955858889112471749374956993360634916541629078307e-50 |
| 0.5 | 200 | -3.67574637578504169896540983775428968999536768792266071728697e-63 |

The 75-digit maximum relative coarse/fine delta is
`5.85979491689557237184563e-46`; the `t=0` rows are independently
checked against derivatives of `xi((1+i*x)/2)/8`, with maximum relative delta
`1.31007236842613583495048e-46`.

These selected rows miss the close-pair stress region and are retained as a
warning against promoting a regular frequency grid.

## Exact Xi Counterexample

At the rational Lehmer stress point `x=1401016343/100000`, Arb certifies

```text
[-3.24186746975685391605462931182063307428944779983857697073123186e-6 +/- 5.42e-69]
[-0.0001463088540165360168258875503621433728403782920777681331619380172 +/- 5.15e-68]
M_0(x)=H_0(x)H_0'''(x)-H_0'(x)H_0''(x)=-L_0'(x)<0
```

The theta tail makes |u|^3*exp(u^2/5)*Phi(u) integrable. Dominated convergence therefore makes the first three x-derivatives of H_t continuous in t on [0,1/5]. The strict value M_0(x)<0 persists for all sufficiently small positive t, rejecting the universal assertion M_t(x)>0 for every x>0 and every 0<t<=1/5.

This rejects only the stronger strict-monotonicity route. It does not reject L_t(x)>0, double-zero transversality, Lambda<=0, or RH.

## Nonpromotion Guard

For the explicit theta-tail countermodel,

```text
phi_(1/10,1/1000)(u)=exp(-u^2-u^4/10-(cosh(4u)-1)/1000)*(1+u^4/10)
L(21/5)=[-0.000213478480591774964507646766853475076598648832979809 +/- 5.09E-55]<0.
```

L(0)>0, the certified L(21/5)<0, and L(x)->0 imply that L'(y)>0 for some y>21/5. Thus uniform strong log-concavity, strict root-variable concavity, and theta-type decay do not imply the monotonicity target.

The generic countermodel had already shown that shape and tail hypotheses do
not supply monotonicity. The exact Xi certificate now closes the stronger
route itself.

## Decision

Do not pursue a global positive sine-transform theorem for `v*K_(1,t)` or a
coupled contour sign for `-L_t'`. The surviving Newman routes are direct
strict Laguerre positivity and, more economically near Lehmer pairs, exclusion
of simultaneous zeros of `(H_t,H_t')` through the corrected `C1`
transversality target.
