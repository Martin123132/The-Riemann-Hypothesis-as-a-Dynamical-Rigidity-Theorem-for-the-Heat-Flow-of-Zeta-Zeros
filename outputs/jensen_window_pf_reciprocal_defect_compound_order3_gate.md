# Jensen-Window PF Reciprocal-Defect Compound Order-Three Gate

Date: 2026-07-12

Status: complete all-column order-three theorem with strict abstract
countermodel gate. This is not a proof of RH or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_reciprocal_defect_compound_order3_gate`.

```text
work/rh_compute/results/jensen_window_pf_reciprocal_defect_compound_order3_gate.json
python work/rh_compute/scripts/jensen_window_pf_reciprocal_defect_compound_order3_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_reciprocal_defect_compound_order3_gate.py
```

Current result:

```text
validated Jensen-window PF reciprocal-defect compound order-three gate: 11 rows, 0 issues, 2 exact coordinate identities, 2 exact sign equivalences, 1 sufficient increment theorem, 1 exact boundary benchmark, 1 strict cone countermodel, 1 all-shift lambda=-100 entry theorem, 1 full forward propagation theorem, 1 arbitrary-column order-three theorem, 1 open order-four handoff
```

## Exact Coordinate

For the contiguous minor

```text
D_(3,n)=det[A_(n+i+j)]_(i,j=0..2)
D_(3,n)=A_n^3*rho_n^6*x_(n+1)^3*F_n
F_n=x_(n+1)*x_(n+2)^2*x_(n+3)-x_(n+1)*x_(n+2)^2-x_(n+2)^2*x_(n+3)+2*x_(n+2)-1
```

all factors outside `F_n` are positive. In defect variables this
collapses to

```text
F_n=d_(n+1)*d_(n+3)*x_(n+2)^2-d_(n+2)^2
D_(3,n)<0 iff d_(n+2)^2>x_(n+2)^2*d_(n+1)*d_(n+3)
```

and in reciprocal defects it becomes

```text
D_(3,n)<0 iff C_n:=q_(n+1)*q_(n+3)-q_(n+2)^2+1>0
d_(n+2)^2-x_(n+2)^2*d_(n+1)*d_(n+3)=C_n*(q_(n+1)*q_(n+3)+q_(n+2)^2-1)/(q_(n+1)^2*q_(n+2)^4*q_(n+3)^2)
```

This is the first explicit higher compound concentration coordinate.

## Increment Budget

Set consecutive reciprocal-defect increments `a_n,b_n`. Then

```text
a_n=q_(n+2)-q_(n+1), b_n=q_(n+3)-q_(n+2)
C_n=1-a_n*b_n+q_(n+2)*(b_n-a_n)
If b_n<a_n, prove q_(n+2)*(a_n-b_n)+a_n*b_n<1.
```

One sufficient theorem is immediate:

```text
0<=a_n<=b_n<=1
C_n>=1-a_n*b_n>=0
For q_k=alpha+beta*k with 0<=beta<1, C_n=1-beta^2>0.
```

Actual Xi increments need not be nondecreasing, so the exact curvature
budget remains the live statement.

## Strict Countermodel

Take

```text
q_1,q_2,q_3=['10', '109/10', '58/5']
d_1,d_2,d_3=['1/100', '100/11881', '25/3364']
x_1,x_2,x_3=['99/100', '11781/11881', '3339/3364']
s_1,s_2,s_3=['3/200', '250/11881', '175/6728']
q increments=['9/10', '7/10']
```

The resulting positive coefficient prefix is

```text
A_0,...,A_4=['1', '1', '99/100', '115465581/118810000', '449662131526605921/474856053604000000']
two cubic frontier values=['-828039/1411581610000', '-1610484375/1597415764323856']
C_0=-181/100
det[A_(i+j)]_(i,j=0..2)=4106267526339/1899424214416000000>0
```

It lies strictly inside the full ratio cone, has decreasing defects,
increasing scaled defects, reciprocal-defect increments below one, and
strictly hyperbolic neighboring cubic Jensen windows. Nevertheless the
order-three signed-Hankel minor has the wrong sign. This blocks promotion
from every currently proved scalar/cubic cone to the new compound margin.

## Lambda=-100 Entry

A repaired Arb prefix through n=317 and an exact scaled-defect tail prove C_n(-100)>0 for every n>=0.

```text
C_n(-100)>0 and D_(3,n)(-100)<0 for every n>=0,
uniform analytic tail lower bound=57613471/66107054971.
```

The finite prefix is spliced to an exact scaled-defect tail theorem;
this is no longer a finite-only observation.

## Forward Invariance

The exact cooperative system C_n'/r_(n+2)=alpha_n*C_(n+1)+beta_n*C_n, the cubic forward-uniform tail, and a weighted maximum principle propagate C_n>0 through lambda=0.

```text
C_n(lambda)>0 for every n>=0 and finite lambda>=-100
D_(3,n)(0)<0 for every n>=0.
```

## Arbitrary Columns

Positive column scaling embeds the three-row Hankel columns as planar points. Contiguous negativity makes their edge slopes strictly decrease, and secant averaging transfers the sign to every strictly increasing column triple.

```text
R_(3,n)(j_1,j_2,j_3)<0 for every n>=0 and 0<=j_1<j_2<j_3 at lambda=0
```

## Live Handoff

Find coordinates, entry theorems, and forward-invariance laws for contiguous order four; then transfer to arbitrary columns and continue to higher compound orders.

All-column reshaped-Hankel order three is closed. Compound order four
and every higher order remain open.
