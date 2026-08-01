# Newman First-Order Centered Absolute-Phase Anchor Reduction

Date: 2026-07-26

Correction notice (2026-07-31): any endpoint-specific declaration that
`H_a` is real is superseded by the complex endpoint source-normalization
gate and Formal Core Section 11.146.  The branch-free anchor and determinant
collapse remain exact.

Status: exact branch-free phase normalization and a decisive
determinant-collapse guard. The direct anchored projection
theorem remains open; this is not a proof of RH.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_first_order_centered_absolute_phase_anchor_reduction.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_absolute_phase_anchor_reduction.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_absolute_phase_anchor_reduction.py
```

## First-Coefficient Anchor

The n=1 term removes every logarithmic phase:

```text
f_1=phi*(1+d_1), |f_1|=|1+d_1|, |d_1|<2189/x<1/2
eta=f_1/|f_1|=phi*(1+d_1)/|1+d_1|
```

At `L=50`, the certified `|d_1|` upper bound is
`3.35978744122717451e-20`, so the anchor is uniformly
nonzero without choosing a branch of `arg`.

Relative coefficients are therefore exact:

```text
q_n=f_n/f_1=exp[t*log(n)^2/4-s_*log(n)]*(1+d_n)/(1+d_1), q_1=1
```

## Endpoint In The Same Frame

The established endpoint phase cancellation gives

```text
beta=(sqrt(pi)/8)*exp(t*pi^2/64-pi*T_0/4)*T_0^(1/2)>0; kappa_N=(-1)^(N+1)*beta*(T_0+i)/|M_t(s)|; g_0=(-1)^N*beta*(T_0+i)*H_a/|M_t(s)|; G_a=(-1)^N*beta*(T_0+i)*(H_(a,x)+mu_a*H_a)/|M_t(s)|
```

After division by the first coefficient, the absolute
normalizer phase cancels algebraically:

```text
r_0=g_0/f_1=(-1)^N*beta*(T_0+i)*H_a/[M_t(s)*(1+d_1)]; r_A=G_a/f_1=(-1)^N*beta*(T_0+i)*(H_(a,x)+mu_a*H_a)/[M_t(s)*(1+d_1)]
```

Define the two aggregate shapes by

```text
P_0=sum_(n=1)^N q_n, P_1=sum_(n=1)^N ell_n*q_n, ell_n=log(n/a); Z_0=P_0+r_0=E_[1]/f_1; Z_A=-s_*'*P_1+r_A=C_a/f_1, C_a=-s_*'*M_(1,a)+G_a, A_a=Re(C_a)
```

Then the retained real value and centered scalar are exactly

```text
(X,A_a)^T=S(Z_0,Z_A)*(Re(f_1),Im(f_1))^T, S=[[Re(Z_0),-Im(Z_0)],[Re(Z_A),-Im(Z_A)]]; X/|f_1|=Re(eta*Z_0), A_a/|f_1|=Re(eta*Z_A)
```

## Exact Contact Split

The contact system now has a branch-free geometric form:

```text
If Z_0!=0, then X=A_a=0 iff Re(eta*Z_0)=0 and Im(conj(Z_0)*Z_A)=0. If Z_0=0, then X=A_a=0 iff Re(eta*Z_A)=0.
```

The second clause is essential. It retains the complex-main-
zero class that any division by `Z_0` would delete.

For the ordinary branch, the natural shape determinant is

```text
Delta_anchor=det(S)=-Im(conj(Z_0)*Z_A); |f_1|^2*Delta_anchor=-Im(C_a*conj(E_[1]))
||(X,A_a)||_2>=|f_1|*|Delta_anchor|/sqrt(|Z_0|^2+|Z_A|^2) when |Z_0|^2+|Z_A|^2>0
```

## Determinant Collapse

The determinant is not a new independent arithmetic route.
Substituting the exact centered derivative decomposition gives

```text
E_[1],x=lambda_a*E_[1]+C_a+D_(1,x), lambda_a=u_a+i*v_a; hence |f_1|^2*Delta_anchor=v_a*|E_[1]|^2+Im(D_(1,x)*conj(E_[1]))-W_[1]
```

At a real-part crossing this reduces to

```text
At X=0, |f_1|^2*Delta_anchor=A_a*Y. If E_[1]=iY!=0, then A_a>0 iff Delta_anchor*Y>0.
```

but at a complex-main zero,

```text
At E_[1]=0 one has Z_0=Delta_anchor=W_[1]=0, while A_a=|f_1|*Re(eta*Z_A) can be nonzero. The determinant cannot classify this crossing.
```

Thus a lower bound for `Delta_anchor` is a stronger sufficient
Wronskian theorem with the same blind spot. It cannot replace
the direct centered scalar.

## Chart-Covariance Guard

```text
The anchor f_1 is independent of the prescribed cutoff. |Re(eta*Delta Z_0)|<2200*exp(-5L/4) and |Re(eta*Delta Z_A)|<10000*exp(-7L/4). The adjacent theorem controls these real projections, not the complex shape increments; Q_N can have a large imaginary part. Therefore Delta_anchor has no certified adjacent-chart absolute bound.
```

This matters because the adjacent cancellation is real, not
complex-absolute. Treating `Delta_anchor` as chart invariant
would discard the cancellation already proved in the endpoint
recurrence.

## Surviving Theorem

The exact q>=1 target is now

```text
On q=2tL^2>=1 in the canonical N=floor(a) chart, prove |Re(eta*Z_0)|<=50000*exp(-5L/4)/|f_1| implies |Re(eta*Z_A)|>(100000L+1)*exp(-5L/4)/|f_1|, then bound N_(Re(eta*Z_0)=0,Re(eta*Z_A)>0) below one successor turn. This formulation retains Z_0=0.
```

A logically complete proof may split the domain as follows:

```text
A possible proof split is: on |Z_0|>=tau use anchored phase-difference information together with the exact determinant collapse; on |Z_0|<tau prove the direct real projection Re(eta*Z_A) bound. The second branch is essential and cannot be inferred from Delta_anchor.
```

The independent q<1 obligation is

```text
Construct the separate q=2tL^2<1 multiplicity-compatible local chart and finite phase-cell theorem; the present first-order anchor does not provide it.
```

This reduction proves no anchored projection lower bound,
signed crossing budget, q<1 chart, finite phase closure,
contact exclusion, `Lambda<=0`, PF-infinity, RH, or
Clay-prize conclusion.
