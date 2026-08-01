# Complex Endpoint Source-Normalization Corrigendum

Date: 2026-07-31

Status: exact correction gate and nonpromotion audit. This is
not a proof of the Xi contact gap, `Lambda<=0`, or RH.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_first_order_centered_complex_endpoint_source_normalization_gate.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_complex_endpoint_source_normalization_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_complex_endpoint_source_normalization_gate.py
```

Current result:

```text
validated complex endpoint source-normalization gate: 16 rows, 1 nonreal source witness, 4 corrected Cartesian projections, 2 contact minors, 2 rank guards, 7 historical artifacts quarantined, 0 signed lower bounds, 0 Abel gaps, 0 winding bounds
```

## Source Correction

Published source: D.H.J. Polymath, Effective approximation of heat flow evolution of the Riemann xi function, and a new upper bound for the de Bruijn-Newman constant, arXiv:1904.12438, equation (53), C_0(p).

https://arxiv.org/abs/1904.12438

```text
C_0(p)={exp(pi*i*(p^2/2+3/8))-i*sqrt(2)*cos(pi*p/2)}/{2*cos(pi*p)}
H_a(p)=C_0(p)+C_0'''(p)/(12*pi^2*a)
```

The endpoint function C_0(p)=F(p) is not the Vaughan low coefficient C_0(n) from the preceding decomposition.

C_0 and H_a are generally complex on the real p-axis; the source gives |C_0(p)|<=1/2, not C_0(p) in R.

At the exact midpoint,

```text
C_0(-p)=C_0(p), hence C_0'''(0)=0
C_0(0)=sqrt(2-sqrt(2))/4+i{sqrt(2+sqrt(2))/4-sqrt(2)/2}
Im C_0(0)<0 because sqrt(2+sqrt(2))<2*sqrt(2)
For every a>0, H_a(0)=C_0(0) is nonreal.
```

p=0 occurs at every half-integer saddle a=N+1/2, including arbitrarily large cutoffs.

## Correct Endpoint

```text
kappa in R, H_a=H_R+iH_I, J_a=J_R+iJ_I, e=kappa(T_0+i)H_a, g=kappa(T_0+i)J_a
Re e=kappa(T_0H_R-H_I), Im e=kappa(T_0H_I+H_R)
Re g=kappa(T_0J_R-J_I), Im g=kappa(T_0J_I+J_R)
g*conj(e)=kappa^2(T_0^2+1)J_a*conj(H_a)
E_N=kappa(T_0+i)[J_a-s_*'log(a)H_a]
```

The centered endpoint atom is

```text
alpha=c*u_N, c_0=Re e, d_0=Re g-alpha*Re e
c_0=kappa(T_0H_R-H_I)
d_0=kappa{T_0J_R-J_I-c*u_N(T_0H_R-H_I)}
h_0=(T_0J_R-J_I)/(T_0H_R-H_I)-c*u_N when T_0H_R-H_I!=0
The division-free endpoint fibre is T_0H_R-H_I=0; H_a=0 is only a subfibre.
```

## Odd Fibre

```text
For the physical endpoint g_0=A_N+iB_N and S_odd=X_odd+iY_odd, Re(S_odd*conj(g_0))=A_NX_odd+B_NY_odd.
|S_odd+g_0|^2=(X_odd+A_N)^2+(Y_odd+B_N)^2
For D=sqrt(A_N^2+B_N^2)>0, P=(A_NX_odd+B_NY_odd)/D and Q=(-B_NX_odd+A_NY_odd)/D give |S_odd+g_0|^2=(P+D)^2+Q^2.
For D=0 retain |S_odd|^2 without division.
```

The rotation follows the actual complex endpoint vector. A
fixed `(T_0,1)` endpoint direction is not source-derived.

## Contact Minors

```text
For z_n=X_n+iY_n, d_n=c(u_n-u_N)X_n-b*u_nY_n.
K_(0n)=c_0d_n-d_0X_n=kappa{(T_0H_R-H_I)u_n(cX_n-bY_n)-(T_0J_R-J_I)X_n}.
K_(nm)=X_nd_m-d_nX_m=c*log(n/m)X_nX_m-b{u_mX_nY_m-u_nY_nX_m}.
```

Neither determinant has a sign from the currently proved amplitude, distance, or source-normalization data.

## Vaughan Rank

```text
For r=(R_0,R_I,R_II) and sigma=(1,-1,-1), Re mathfrak B_N=sigma*r.
sigma*sigma^T=[[1,-1,-1],[-1,1,1],[-1,1,1]] has rank 1 and eigenvalues 3,0,0.
After pullback to the endpoint/carrier atoms the square is d*d^T, the original centered-slope rank-one kernel.
Adding the value row gives c*c^T+d*d^T of rank at most 2; on c*x=0 only the rank-at-most-one slope square remains.
```

Exact Vaughan recombination changes coordinates but cannot create a contact margin. A new Xi restriction on the actual coefficient curve is still required.

## Supersession

Only endpoint rows that used H_a in R are superseded. The branch-free anchor, complex factorization, abstract component transport, Abel/Mangoldt identities, adjacent cutoff recurrence, and endpoint-composed Vaughan identity remain exact.

- `absolute_phase_anchor_reduction: real-H declaration only`
- `direct_projection_regime_reduction: endpoint projections`
- `carrier_kernel_abel_prefix_reduction: endpoint Hermitian product`
- `endpoint_first_pivot_odd_small_ball_guard: real-H sentence only`
- `endpoint_odd_fibre_correlation_feasibility_guard: fixed endpoint direction`
- `contact_signed_transport_reduction: endpoint atom and exceptional fibre`
- `interior_projective_current_gate: endpoint effective-slope text`

Use the corrected Cartesian endpoint atom in every later contact or Type-I/II argument. Do not infer a fixed (T_0,1) endpoint direction, and do not divide merely under H_a!=0.

## Pi Provenance

Every pi in this correction is inherited from the published completed-zeta/Riemann-Siegel endpoint normalization. The Cartesian repair introduces no new pi.

## Boundary

This artifact corrects the recurrent endpoint's Cartesian projections, Hermitian product, centered atom, exceptional fibre, odd-fibre rotation, and contact minors. It proves that the exact Vaughan split has no new rank capable of supplying the missing contact margin. It does not prove an Xi phase law, signed endpoint/carrier correlation estimate, Abel-scalar gap, successor winding cap, contact exclusion, Q209, the cofinal descendant theorem, Lambda<=0, PF-infinity, RH, or a prize-level conclusion.
