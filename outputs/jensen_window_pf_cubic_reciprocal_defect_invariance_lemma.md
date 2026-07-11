# Jensen-Window PF Cubic Reciprocal-Defect Invariance Lemma

Date: 2026-07-10

Status: exact cubic reciprocal-defect invariance reduction with finite entry
certificates. This is not a proof of all-degree Jensen hyperbolicity, PF-infinity,
RH, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_cubic_reciprocal_defect_invariance_lemma.json
python work/rh_compute/scripts/jensen_window_pf_cubic_reciprocal_defect_invariance_lemma.py
python work/rh_compute/scripts/check_jensen_window_pf_cubic_reciprocal_defect_invariance_lemma.py
```

Current result:

```text
validated Jensen-window PF cubic reciprocal-defect invariance lemma: 10 rows, 0 issues, 6 exact coordinate rows, 1 conditional maximum principle, 318 lambda=-100 prefix margins, 310 nonnegative-grid margins, 0 failed or inconclusive, 1 open tail handoff
```

## Exact Cubic Coordinate

Set `d_k=1-x_k` and `q_k=d_k^(-1/2)`. For the normalized cubic
frontier `F=x^2*y^2-6*x*y+4*x+4*y-3`, exact factorization gives

```text
q_k^4*q_(k+1)^4*F
 =(q_k-q_(k+1)-1)*(q_k-q_(k+1)+1)
  *(q_k+q_(k+1)-1)*(q_k+q_(k+1)+1).
```

Inside the full ratio cone, `q_(k+1)>=q_k>=1`. Since the cubic
discriminant is `-27*x^2*F`, degree-3 hyperbolicity is exactly

```text
0 <= q_(k+1)-q_k <= 1.
```

## Boundary Flow

At the nontrivial boundary, write `t=sqrt(1-x)`. Then

```text
x=1-t^2,
y=(1+2*t)/(1+t)^2,
d_j=t^2/(1+(j-1)*t)^2.
```

This is exactly the reciprocal-square defect law from the rank-two boundary
family. Substitution into the heat ratio ODE gives

```text
partial_lambda F/r_(n+1)=C_n(t)*(z-Z(t)),
C_n(t)=8*t^3*(2*n+7)*(1-t)*(1+2*t)^2/(1+t)^3 > 0,
Z(t)=(1+t)*(1+3*t)/(1+2*t)^2,
1-Z(t)=t^2/(1+2*t)^2.
```

Therefore, when `q_(k+1)-q_k=1`, the cubic frontier points inward
exactly when `q_(k+2)-q_(k+1)<=1`. The threshold is independent of the
shift. If the increments tend uniformly to zero in the spatial tail, a
finite-active-set first-crossing argument propagates the whole cubic cone.

## Finite Arb Entry

At `lambda=-100`, all `318` repaired-prefix margins
through `k=318` are strictly positive. The minimum is

```text
k=1: [0.801865140205956003946528256176261697433175013244485162830701160752896247320145010016140934856495888 +/- 8.66e-100]
```

Across the five cached nonnegative heat times, all `310`
available margins are strictly positive. The global finite minimum is

```text
lambda=0.1, k=1:
[0.673764497340926782808244698920510186865940542132052945121961008099664296554818336635 +/- 5.89e-85]
```

## Remaining Theorem

The cubic propagation route now needs one specific analytic input:

```text
q_(k+1)-q_k -> 0 uniformly on compact lambda intervals,
and q_(k+1)-q_k<1 for every k>=319 at lambda=-100.
```

That would promote the existing finite entry and exact boundary algebra to all
degree-3 shifts. Higher degrees still need additional minor coordinates, so this
is a controlled cubic advance rather than the missing all-order bridge.
