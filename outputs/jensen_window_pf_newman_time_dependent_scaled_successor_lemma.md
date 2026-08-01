# Newman Time-Dependent Scaled-Jet Successor Lemma

Date: 2026-07-25

Status: exact generalized successor and relative-transport lemma with an open two-regime Xi antecedent; not a proof of an all-stage successor, `Lambda<=0`, or RH.

## Positive Scaled Jet

Let `H_t=-H_xx` and choose any continuous scale

```text
s(t,x)>0,                    V_s=(H,sH_x).          (1)
```

Then `V_s=0` exactly when `(H,H_x)=0`. On any contact-free boundary,

```text
A_r=diag(1,(1-r)+r s),       0<=r<=1,               (2)
det(A_r)=(1-r)+r s>0
```

is an explicit homotopy from the ordinary first jet to (1). Hence positive
scaling preserves the contact set, local index, boundary winding, and degree.

## Exact Transport Derivative

Unlike the earlier scale `1/ell(x)`, the new scale may depend on time.
Direct differentiation gives

```text
partial_t V_s=(-H_xx,s_t H_x-sH_xxx).              (3)
```

Thus the original additive phase-cell theorem remains valid verbatim after
its derivative bound is replaced by the norm of (3). If

```text
V_s(t_j,I) subset C,          d=dist(0,C)>0,
M>=sup_collar ||partial_t V_s||,
```

then

```text
||V_s(t,x)||>=d-|t-t_j|M.                           (4)
```

## Relative Transport

The time-dependent coordinate also exposes a stronger interface. Suppose
on one descendant collar

```text
||partial_t V_s(t,x)||<=kappa(t)||V_s(t,x)||,       (5)
```

where `kappa` is integrable over that positive-time collar. Gronwall in
reverse time gives

```text
||V_s(t,x)||
 >=exp(-integral_t^t_j kappa(r)dr)||V_s(t_j,x)||.   (6)
```

Consequently an old cell with clearance `d_j` has descendant clearance

```text
d_j exp(-integral_t^t_j kappa)>0.                  (7)
```

No additive ratio below one and no floor uniform as `t` tends to zero is
needed. A finite relative integral on each individual collar is enough.

## Hermite Specialization

For the multiplicity-compatible scale

```text
s=sqrt(2t),                   s_t=1/s,              (8)
```

the exact Hermite benchmark supplies

```text
kappa_m(t)=C_m/(2t).                                 (9)
```

On the cofinal times `t_j=1/(5j)`, (6) becomes

```text
||V_s(t_(j+1),x)||
 >=(j/(j+1))^(C_m/2)||V_s(t_j,x)||,                (10)
```

and the one-step logarithmic cost
`(C_m/2)log(1+1/j)` tends to zero. Every factor in (10) is positive,
although the endpoint limit may have zero clearance.

## Scale Interface

Entry-strip arithmetic may be cleaner in a frequency scale such as `1/L`,
whereas old-cell multiplicity transport is cleaner in `sqrt(2t)`. Any two
positive scales have the interpolation

```text
s_r=(1-r)s_entry+r s_descendant>0.                 (11)
```

Equation (2) applied to (11) proves that switching coordinates cannot create
a contact or alter degree. A rigorous implementation still has to bound the
transport derivative through a chosen overlap; the topology itself is exact.

## Conditional Successor

For

```text
P_j=[t_j,1/5]x[0,R_j],
C_j=[t_(j+1),t_j]x[0,R_j],
S_j=[t_(j+1),1/5]x[R_j,R_(j+1)],
```

a multiplicity-compatible all-stage theorem would follow from:

```text
1. one certified finite base P_J;
2. a strict first-jet half-plane cone on each new strip S_j in a positive
   entry scale;
3. an integrable relative bound (5) on every old-cell descendant collar,
   possibly after a positive scale switch;
4. explicit finite interfaces and source/remainder bounds.
```

The ordinary successor proof then applies to the scaled boundary paths,
and (2) transfers the resulting zero degree back to the ordinary first jet.

## Proof Boundary

The contact equivalence, positive-scaling degree homotopy, derivative (3),
additive and relative transport inequalities, cell-clearance factor,
Hermite specialization, cofinal factor, and positive scale switch are
exact. The Xi entry cone, relative descendant estimate, overlap bounds,
finite transition, and all-`j` composition are open. This lemma proves no
Q209 theorem, `Lambda<=0`, RH, PF-infinity, or Clay-prize conclusion.

Machine-audited files:

```text
work/rh_compute/results/jensen_window_pf_newman_time_dependent_scaled_successor_lemma.json
work/rh_compute/scripts/jensen_window_pf_newman_time_dependent_scaled_successor_lemma.py
work/rh_compute/scripts/check_jensen_window_pf_newman_time_dependent_scaled_successor_lemma.py
```
