# Newman Adiabatic Phase-Cell Successor Lemma

Date: 2026-07-25

Status: exact conditional successor theorem with finite Q208
calibration. The Xi all-j antecedent remains open; this is
not a proof of `Lambda<=0` or RH.

## Shell Decomposition

```text
P_j=[t_j,1/5]x[0,R_j]
C_j=[t_(j+1),t_j]x[0,R_j]
S_j=[t_(j+1),1/5]x[R_j,R_(j+1)]
P_(j+1)=P_j union C_j union S_j
delta_j=1/(5*j*(j + 1))
```

A successor is not merely a new corner: lowering the bottom
time adds a collar across the complete old x-range. The
one-unit right strip is the only spatially local addition.

## Heat-Jet Transport

For a positive scale `ell(x)` independent of time,

```text
V_ell=(H,H_x/ell(x)), ell(x)>0
partial_t V_ell=(-H_xx,-H_xxx/ell)
||partial_t V_ell||_2^2=h_xx**2 + h_xxx**2/ell**2
```

On a phase cell `I_(j,k)`, prove

```text
V_ell(t_j,I_(j,k)) subset C_(j,k)
d_(j,k)=dist(0,C_(j,k))>0
M_(j,k)>=sup sqrt(H_xx^2+(H_xxx/ell)^2) on the collar
delta_j*M_(j,k)<d_(j,k)
```

Then reverse triangle and convex-cell thickening give

```text
||V_ell(t,x)||_2>=d_(j,k)-delta_j*M_(j,k)>0.
```

Thus the actual heat flow is an origin-free homotopy between
the old and new bottom paths.

## Right Strip

It is enough to find a fixed vector `q_j` and `eta_j>0`
with

```text
q_j dot V_ell(t,x)>=eta_j on S_j.
```

The strip then lies in an open half-plane and contains no
contact. Q208 calibrates this with `q=(0,1)` and `F_x>0`
on all 20 stored right-strip cells.

## Successor Theorem

```text
P_j contact-free
and every bottom cell satisfies delta_j*M_(j,k)<d_(j,k)
and the new right strip satisfies one strict half-plane gate
imply P_(j+1) contact-free with zero first-jet degree.
```

A finite base plus these hypotheses for every later stage
would certify the cofinal exhaustion and hence `Lambda<=0`.

## Live Obligation

Construct a rigorous collar certifier for sqrt(H_xx^2+(H_xxx/ell)^2), first on the Q207-to-Q208 collar as a calibration, and compare delta_j*M_(j,k) with exact phase-cell distances. Then seek an analytic all-j bound before using Q209 as anything beyond a falsification case.

## Proof Boundary

The shell decomposition, heat-jet transport identity, cell-thickening criterion, half-plane criterion, and successor induction are exact. No Q207-to-Q208 transport budget, no all-j bottom-collar estimate, no all-j right-strip cone, no cofinal Xi theorem, no Lambda<=0, and no RH proof is supplied.
