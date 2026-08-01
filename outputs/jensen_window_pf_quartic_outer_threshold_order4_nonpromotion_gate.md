# Quartic Outer-Threshold Order-Four Nonpromotion Gate

Date: 2026-07-25

Status: exact finite local scope gate. This is not a proof of
quartic invariance, PF-infinity, `Lambda<=0`, or RH, and it is not
an Xi counterexample or an infinite sign-regular countermodel.

```text
work/rh_compute/results/jensen_window_pf_quartic_outer_threshold_order4_nonpromotion_gate.json
python work/rh_compute/scripts/jensen_window_pf_quartic_outer_threshold_order4_nonpromotion_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_quartic_outer_threshold_order4_nonpromotion_gate.py
```

Current result:

```text
validated quartic outer-threshold order-four nonpromotion gate: 11 rows, 0 issues, 1 exact outer quartic contact, 8 strict ratio coordinates, 7 strict scaled-defect steps, 7 reciprocal-defect increment bounds, 7 strict cubic frontiers, 120 signed order-two minors, 126 signed order-three minors, 56 signed order-four minors, 1 failed outer threshold, 1 negative quintic discriminant, 1 forbidden promotion, 1 closed downstream handoff
```

## Exact Contact

Take

```text
(1+3719*w/5000)^2*(1+4943*w/5000)*(1+7619*w/5000)
2*a+b+c=4.
```

Its contraction coordinates, followed by two exact extensions, are

```text
x_2=24154639/25000000
x_3=567331181410000/583446585220321
x_4=909681023616163/930852655574884
x_5=2453/2500
x_6=123/125
x_7=19719/20000
x_8=987309/1000000
x_9=98823997/100000000
```

Regard these as `x_2,...,x_9`. Set `A_1=A_2=1`,
`r_j=r_(j-1)x_j`, and `A_(j+1)=A_j r_j`. This gives a positive
rational segment `A_1,...,A_10`. Every contraction is strictly
below one, strictly above its pointwise wall
`(2j-1)/(2j+1)`, and the eight contractions increase strictly.
In addition, the scaled defects
`s_j=(2j+1)(1-x_j)/2` increase at all seven steps, and exact
one-squaring certificates prove
`(1-x_(j+1))^(-1/2)-(1-x_j)^(-1/2)<1` at all seven steps.

All seven supported adjacent cubic frontiers are strict:

```text
F_1=-27542245385772/364654115762700625<0
F_2=-5309435528978030980920000/135775700809608417375013504441<0
F_3=-107160123180070672535747959690079/5415541664942585129778322584100000000<0
F_4=-1270039/97656250000<0
F_5=-60359031/6250000000000<0
F_6=-3063049172838759/400000000000000000000<0
F_7=-64094836187510322510671/10000000000000000000000000000<0
```

## Signed-Hankel Audit

For every shift and strictly increasing column set whose entries
remain in `A_1,...,A_10`, exact rational enumeration gives

```text
order 2: 120/120 required signed margins positive
order 3: 126/126 required signed margins positive
order 4:  56/56  required signed margins positive
```

The closest exact margins are

```text
order 2: 13166581882944804486567406270206319940986855650073384441926610688804320054395501922559959737635920383146323215605942671964340489148392192866833144671276689492841391468685276480984938647233643/4443772403973784816914603573505985700810904260387334510596701875329017639160156250000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
order 3: 401842677616759001417844284705536349262609673435336911656193470737864254413284387236276052794891810718116479047412703806484845227708746363462717226905208817560048207068828998661229273650648326849806293708252121503/7414474041680064294814205000954856187581221331083473911284887265082943486049771308898925781250000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
order 4: 12574434445944550439275379655690051061196184174134414753333789708491269843317673678924926107957838214818897153882034836078104890312654672824881100214789922453434934294792092805044099063684374408577034531211864875149/10419431193301813813878870292956905925436505008974599765569024611977511085569858551025390625000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
```

Thus the witness includes every local arbitrary-column order-three
test and the available strict contiguous order-four test, not merely
the boundary determinant used to select the outer branch.

## Failed Outer Threshold

At the quartic contact,

```text
C=(a-b)(a-c)=5967/31250>0
U(a,p)=1391468765486968/1418329604954089
u-U=490607234960317/3545824012385222500>0.
```

The repeated root is on the live outer branch, but the required
`u<=U(a,p)` condition fails strictly. The heat vector is outward.

The adjacent quintic exposes the missing structure exactly:

```text
P_5(-1/a)=-2*p^2*(u-U)/(a^2*B)
             =-490607234960317/1054952331831895000<0
Disc(P_5)=-1051319854565846275082407616343548911197300752504299706357671105291155412668332062800253267699328157259919/1032932599056736970859207646198272705078125000000000000000000000000000000000000000000000000000000000000000000000000<0.
```

## Consequence

The completed local signed-Hankel information through order four
cannot by itself prove the remaining quartic inward threshold.
This does not challenge the Xi theorems: the witness is finite and
is neither the Xi sequence nor a Newman trajectory. It instead
forces the next proof step to use genuinely nonlocal/far-column
information, a degree-five polar closure, theta arithmetic, or
another Xi-specific heat-compatible invariant.

The immediate continuation problem is exact rather than merely
numerical. Appending `x_10` preserves the new terminal order-three
and order-four signs whenever

```text
122559637567354810219/123942620623363810219 < x_10
x_10 < 2855598054244597686471023128949393139897289/2887821089634174590758480845749303731397289
corridor width=18955815844931742539348666080562500000/36649337448547309731315880413404413655162994699>0.
```

The rational point `x_10=4944208739/5000000000` lies
strictly inside and also preserves the next pointwise wall,
scaled-defect increase, reciprocal-increment bound, and cubic
frontier. The exact downstream obstruction now settles this
particular all-length question: after arbitrary strict corridor
choices through `x_12`, the best possible compatibility margin
for an increasing `x_13` and the next order-four sign is still
negative. Thus this fixed prefix cannot be an infinite
countermodel. A uniform theorem over all outer-contact prefixes,
or a genuinely global Xi, degree-five, or theta-specific
constraint, remains open.

```text
outputs/jensen_window_pf_quartic_outer_branch_length13_obstruction.md
```
