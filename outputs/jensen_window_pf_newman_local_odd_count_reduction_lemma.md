# Jensen-Window PF Newman Local Odd-Count Reduction Lemma

Date: 2026-07-11

Status: exact local odd-count reduction with counting and birth
guards. This is not a proof of RH or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_newman_local_odd_count_reduction_lemma.json
python work/rh_compute/scripts/jensen_window_pf_newman_local_odd_count_reduction_lemma.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_local_odd_count_reduction_lemma.py
```

Current result:

```text
validated Jensen-window PF Newman local odd-count reduction lemma: 10 rows, 0 issues, 3 exact Stieltjes identities, 1 explicit outer bound, 1 log-squared localization theorem, 1 odd-count formula, 3 finite reciprocal-gap checks, 1 published uniform counting input, 1 exact classical-field birth countermodel, 1 open Xi collision-exclusion handoff
```

## Stieltjes Reduction

For a positive double zero `c`,

```text
At an even double zero c>0, remove the positive double atom and let mu count the other positive roots. Then B(c)=1/c+PV integral K_c(y)dmu(y), K_c(y)=2c/(c^2-y^2)=1/(c-y)+1/(c+y).
K_c'(y)=4c*y/(c^2-y^2)^2>0
```

Let `F` be the difference between the positive outsider-root count
and the positive-time reference count. Integration by parts gives

```text
If |F(y)|<=E on [a,2c] and |F(y)|<=C*log(2+y) for y>=2c, then for 0<H<=c/2 the field discrepancy outside (c-H,c+H) is at most 5*E/H+2*C*(log(4c)+1)/c.
```

The published theorem supplies `E=O(log(4c))` uniformly for
`0<t<=1/2`. Choosing `H=log(4c)^2` and using the classical
continuum field therefore gives

```text
Uniformly for 0<t<=1/2 at any large positive double zero c, B(c)=-pi/8+S_H(c)+O(1/log(4c)+log(4c)^3/c), H=log(4c)^2.
```

All mesoscopic and far-field uncertainty is now `O(1/log c)`.
The unresolved quantity is genuinely local.

Primary counting source: https://arxiv.org/abs/1904.12438.

## Odd Local Statistic

Define the left-minus-right outsider count at radius `u` by `D_c(u)`.
Then

```text
Let D_c(u)=#roots in [c-u,c)-#roots in (c,c+u]. If H avoids roots, S_H=D_c(H)/H+integral_0^H D_c(u)/u^2 du.
For paired left/right distances ell_m,r_m, 1/ell_m-1/r_m=(r_m-ell_m)/(ell_m*r_m); unmatched roots retain their signed reciprocal terms.
W_H=|D_c(H)|/H+integral_0^H |D_c(u)|/u^2 du controls |S_H|.
```

This is the strict local target. A global `O(log c)` counting error
does not control the inverse-square weight at small `u`.

## Birth Guard

Set

```text
a^2=1 + 16/(pi + 8)
P(z)=(z^2-1)^2*(z^2-a^2)
F_tau=exp(-tau*d_z^2)P
B(1)=-pi/8, center velocity=-pi/4
z_+/-=1+/-sqrt(2*tau)-pi*tau/4+O(tau^(3/2)), tau>0
```

This is an exact even solution of `F_tau=-F_zz`. It has the
classical field `-pi/8` and drift `-pi/4`, yet a real pair is born
for `tau>0`. Therefore even a successful odd-count theorem would
control the field but would not, by itself, exclude a positive Newman
boundary.

## Live Handoff

Prove an Xi-specific lambda-uniform bound on the signed odd-count integral or paired reciprocal-gap discrepancy, and then add a genuinely global mechanism that turns that control into collision exclusion; field balance alone is compatible with positive birth.

The next proof mechanism must use an additional Xi-specific global
constraint beyond the first regularized root field.
