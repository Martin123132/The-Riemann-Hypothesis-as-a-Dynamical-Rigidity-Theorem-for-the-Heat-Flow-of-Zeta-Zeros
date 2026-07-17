# Jensen-Window PF Compound Order-Four Forward-Flow Reduction

Date: 2026-07-13

Status: exact cooperative order-four flow reduction with one open
spatial-tail coefficient bound. This is not a proof of forward
order-four invariance, PF-infinity, RH, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_compound_order4_forward_flow_reduction.json
python work/rh_compute/scripts/jensen_window_pf_compound_order4_forward_flow_reduction.py
python work/rh_compute/scripts/check_jensen_window_pf_compound_order4_forward_flow_reduction.py
```

## Affine Hankel Flow

The coefficient heat equation is

```text
A_j'=(4*j+2)*A_(j+1)
```

If `delta(A_j)=A_(j+1)`, exact determinant algebra gives

```text
H_(m,n)'=(4*(n+2*m-2)+2)*delta(H_(m,n))
H_(3,n+1)*delta(H_(4,n))=H_(3,n)*H_(4,n+1)+delta(H_(3,n+1))*H_(4,n)
```

Put `T_n=-H_(3,n)>0` and `Q_n=H_(4,n)`. The completed order-three
theorem therefore turns the order-four flow into

```text
Q_n'=a_n*Q_(n+1)+b_n*Q_n,
a_n=(4*n+26)*H_(3,n)/H_(3,n+1) >0,
b_n=((4*n+26)/(4*n+22))*H_(3,n+1)'/H_(3,n+1).
```

At `Q_n=0`, the vector field points inward whenever `Q_(n+1)>=0`.
Thus the finite-dimensional boundary algebra is cooperative; order five
does not enter this local boundary identity.

## Stable Margin

The positive ratio factorization is

```text
F_n=G_(n+1)^2-x_(n+3)^3*G_n*G_(n+2)
S_n=A_n^4*rho_n^12*x_(n+1)^8*x_(n+2)^4/d_(n+3)
H_(4,n)=S_n*F_n
```

Inside the completed ratio and order-three cones, `0<G_j<=1`, so
`|F_n|<=1`. Positive rescaling preserves the one-sided system:

```text
F_n'=alpha_n*F_(n+1)+beta_n*F_n
alpha_n=(4*n+26)*rho_(n+2)*x_(n+3)^4*(d_(n+3)/d_(n+4))*(G_n/G_(n+1)) >0,
beta_n=((4*n+26)/(4*n+22))*(log T_(n+1))'-(log S_n)'.
```

## Infinite-Index Target

Set `z_n=F_n/(n+1)`. Boundedness gives `z_n->0` uniformly, so a
negative spatial infimum is attained. The exact weighted flow is

```text
z_n'=alpha_n*((n+2)/(n+1))*z_(n+1)+beta_n*z_n
```

and the standard exponential minimum argument closes as soon as, for
every finite `L>=-100`,

```text
for every finite L>=-100, sup_(-100<=lambda<=L,n>=0) [beta_n+alpha_n*(n+2)/(n+1)]<infinity
```

This is now the sole forward-propagation blocker. It is a scalar spatial
tail theorem, not a missing local sign identity. A proof needs quantitative
control of neighboring reciprocal-defect increments and the positive
order-three gaps on compact heat intervals.

```text
outputs/jensen_window_pf_compound_order4_m100_entry_certificate.md
outputs/jensen_window_pf_compound_order3_forward_invariance_certificate.md
outputs/signed_hankel_jensen_dependency_graph.md
outputs/formal_core.md
```
