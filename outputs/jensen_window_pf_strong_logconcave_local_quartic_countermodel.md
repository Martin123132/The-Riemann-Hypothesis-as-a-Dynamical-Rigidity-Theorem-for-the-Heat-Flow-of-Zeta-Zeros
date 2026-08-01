# Jensen-Window PF Strong-Log-Concave Local Quartic Countermodel

Date: 2026-07-25

Status: exact interval local quartic countermodel gate. This is not a proof
or disproof of degree-four Xi hyperbolicity, PF-infinity, RH, or
`Lambda <= 0`; the witness is not the Xi kernel.

```text
work/rh_compute/results/jensen_window_pf_strong_logconcave_local_quartic_countermodel.json
python work/rh_compute/scripts/jensen_window_pf_strong_logconcave_local_quartic_countermodel.py
python work/rh_compute/scripts/check_jensen_window_pf_strong_logconcave_local_quartic_countermodel.py
```

Current result:

```text
validated Jensen-window PF strong-log-concave local quartic countermodel: 10 rows, 0 issues, 6 positive Mellin moments, 4 local ratio-wall contractions, 3 monotone gaps, 3 strict cubic margins, 1 negative quartic frontier, 1 full-support approximation theorem, 1 Xi-specific handoff
```

## Strong-Log-Concave Mellin Witness

On the convex support `[1/20,1]`, take

```text
f(y)=exp(-3*y-y^2/10)*1_[1/20,1](y)
A_k=integral_(1/20)^1 y^(k-1/2)f(y)dy/Gamma(k+1/2).
```

The potential `3*y+y^2/10` has second derivative `1/5`, so this is
a `1/5`-strongly log-concave squared-variable density. Arb encloses
`A_0,...,A_5` strictly above zero at 512-bit precision.

## Low-Degree Data

For `x_k=A_(k-1)A_(k+1)/A_k^2`, Arb gives

```text
x_1=[0.546260487327330744593299720444906878822025179257789770409088 +/- 4.00E-61]
x_2=[0.818498181607290144339963533721499951426211436768873476867139 +/- 2.03E-61]
x_3=[0.851638198486083412271372927980468278727325989879406655860049 +/- 3.94E-61]
x_4=[0.863362882556287015914332781935164391033742149130755272172689 +/- 4.36E-61]
```

All four satisfy

```text
(2*k-1)/(2*k+1)<x_k<1, 1<=k<=4,
x_1<x_2<x_3<x_4.
```

The three consecutive cubic frontier balls are

```text
F_shift_0=[-0.0237343900999162294985373243981561163330896523419611882190875 +/- 7.45E-63]<0
F_shift_1=[-0.0159417188925092933088746142017136146072988159580009747333247 +/- 3.58E-62]<0
F_shift_2=[-0.0110064304334018563944386285045908398412396185477237133973711 +/- 1.92E-62]<0
```

Thus the local quadratic walls, Stieltjes lower walls, monotone
contractions, and shifted cubic tests are all strict.

## Quartic Failure

The normalized quartic obeys

```text
J_4(w)=1+4*w+6*x_1*w^2+4*x_1^2*x_2*w^3+x_1^3*x_2^2*x_3*w^4
Disc(J_4)=256*x_1^6*x_2^2*Q(x_1,x_2,x_3)
Q=[-0.00587537323825915999878047191451202719150512396456697171798624 +/- 2.83E-63]<0.
```

Hence the shift-zero quartic has negative discriminant and a nonreal
conjugate pair. This is
a moment-realized strengthening of the abstract quartic boundary-flow
obstruction: the currently reconciled local low-degree facts still do
not provide the missing quartic invariant.

## Full-Support Guard

For finite `L>0`, define on `y>=0`

```text
f_L(y)=exp(-3*y-y^2/10-L*((1/20-y)_+ +(y-1)_+)), y>=0, L<infinity
```

Its potential is `1/5`-strongly convex and the density is positive on
the full half-line. As `L` tends to infinity, its first six normalized
moments converge to those above by dominated convergence. Every sign
used here is strict and depends continuously on those moments, so all
of the displayed countermodel inequalities persist for sufficiently
large finite `L`. Full support and strong log-concavity therefore do not
repair the local implication.

## Proof Boundary

The witness is not a Newman heat trajectory and does not satisfy or
challenge the global Xi all-shift theorem as a dynamical statement.
It blocks only a local generic promotion. The surviving target is a
heat-compatible branch-aware quartic invariant, an Xi-specific weighted
theta/Laguerre composition, or a noncircular all-degree theorem.

```text
outputs/jensen_window_pf_kernel_mellin_upper_wall_certificate.md
outputs/jensen_window_pf_cubic_forward_uniform_tail_certificate.md
outputs/jensen_window_pf_quartic_double_root_threshold_lemma.md
outputs/jensen_window_pf_quartic_boundary_flow_obstruction.md
```

Summary:

A 1/5-strongly log-concave Gamma-normalized Mellin density satisfies all local ratio walls and three consecutive strict cubic tests while its shift-zero quartic discriminant is negative. Full-support strongly log-concave approximants preserve the strict signs, so the quartic bridge needs genuinely global or Xi-specific structure beyond the reconciled low-degree facts.
