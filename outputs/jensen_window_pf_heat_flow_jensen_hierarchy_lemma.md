# Jensen-Window PF Heat-Flow Jensen Hierarchy Lemma

Date: 2026-07-10

Status: exact heat-flow hierarchy lemma and open higher-minor handoff. This is not a
proof of Jensen hyperbolicity, PF-infinity, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_heat_flow_jensen_hierarchy_lemma`.

```text
work/rh_compute/results/jensen_window_pf_heat_flow_jensen_hierarchy_lemma.json
python work/rh_compute/scripts/jensen_window_pf_heat_flow_jensen_hierarchy_lemma.py
python work/rh_compute/scripts/check_jensen_window_pf_heat_flow_jensen_hierarchy_lemma.py
```

Current result:

```text
validated Jensen-window PF heat-flow Jensen hierarchy lemma: 9 rows, 0 issues, 5 exact hierarchy identities, 2 cubic countermodels, 1 open higher-minor handoff, 0 ready-to-apply rows
```

## Exact Hierarchy

For

```text
J_(d,n)(z)=sum_(j=0)^d binom(d,j)*A_(n+j)*z^j
A_k'=2*(2*k+1)*A_(k+1),
```

coefficient comparison gives

```text
J_(d+1,n)=J_(d,n)+z*J_(d,n+1),
partial_z J_(d,n)=d*J_(d-1,n+1),
partial_lambda J_(d,n)=(4*n+2)*J_(d,n+1)+4*z*partial_z J_(d,n+1),
partial_lambda J_(d,n)=((4*n+2)*partial_z+4*z*partial_z^2)*J_(d+1,n)/(d+1).
```

The heat evolution is therefore coupled across shifts and degrees. These identities do
not supply a closed scalar PDE for one fixed Jensen window.

## Cubic Gate

The normalized degree-3 discriminant is

```text
Disc(1+3*z+3*x*z^2+x^2*y*z^3)
  = -27*x^2*(x^2*y^2-6*x*y+4*x+4*y-3).
```

At `x=1/2,y=1`, the infinite contraction sequence lies in the full static ratio cone
and has completely monotone defects, but the discriminant is `-27/16`. Thus a
static higher-minor proof needs an additional invariant.

There is also an exact dynamic boundary guard. At shift `n=0`, set

```text
d_k=(4/9)*(9/25)^(k-1),
x=5/9, y=21/25, z=589/625,
F=x^2*y^2-6*x*y+4*x+4*y-3=0,
(partial_lambda F)/r_0=329728/2109375>0.
```

Since the discriminant is a positive factor times `-F`, its heat derivative is
negative at this boundary point. The defects form the one-atom Hausdorff moment
sequence `(4/9)*delta_(9/25)` and satisfy every full-cone wall. Thus even complete
defect monotonicity plus the local heat ODE does not supply forward cubic invariance. A viable
higher-minor flow proof needs an additional shift-coupled invariant.

## Handoff

The lambda derivative of a Jensen window is exactly coupled to shift n+1 and, equivalently, to degree d+1. This gives a concrete hierarchy for higher-minor theorem search. The -27/16 static cubic guard and the exact one-atom Hausdorff boundary-exit guard 329728/2109375 prove that neither the propagated static cone nor complete defect monotonicity plus the local heat ODE supplies the missing invariant structure.

```text
outputs/jensen_window_pf_heat_flow_infinite_cone_invariance_certificate.md
outputs/jensen_window_pf_defect_complete_monotonicity_scout.md
outputs/jensen_window_pf_bridge_target.md
```
