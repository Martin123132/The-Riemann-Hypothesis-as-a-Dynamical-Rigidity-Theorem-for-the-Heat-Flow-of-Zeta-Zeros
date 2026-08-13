# Morse-Fresnel Observation-Image Compression Gate

Date: 2026-08-02

Status: exact eight-observation compression and degree-five polynomial-amplitude
closure; `0 coefficient-aware interior bounds`, `0 signed flow bounds`, and
this is not a proof of RH.

## Observation Image

For a=U_0^Tp and b=U_1^Tp, g_p=2M_xi p+ell_xi=U_0(2Jb)+U_1(2Ja+t_T).

For x=U_0^Tr and z=U_1^Tr, g_p^Tr+r^TM_xi r=2a^TJz+2b^TJx+t_T^Tz+2x^TJz.

The remainder affects the flow through x and z only: four real observation errors and their four frequency derivatives. Any remainder in ker(U_0^T) intersect ker(U_1^T) is exactly invisible.

Writing the four base errors as `x=(x_V,x_N,x_A,x_Q)` and their frequency
derivatives as `z=(z_V,z_N,z_A,z_Q)`, the exact error is

```text
The exact flow error is (N_p+A_(T,x))z_V+(X_T+V_p)z_N-(Q_p+X_(T,x))z_A-(A_T+A_p)z_Q+N_(p,xi)x_V+V_(p,xi)x_N-Q_(p,xi)x_A-A_(p,xi)x_Q+x_Vz_N+x_Nz_V-x_Az_Q-x_Qz_A.
```

## Polynomial-Amplitude Closure

For P(lambda)=sum_(j=0)^5 p_jlambda^j, set A_P(u)=P(log u)exp(t(log u)^2/4-sigma log u). The Morse maps A to b, c, c', and the outside C_rA operator linearly.

Define R[P]=kappa sum_T e(phi_r(alpha/r)) integral c'_(P,r)e(-y^2/2)+kappa^2 sum_(r notin T) integral e(phi_r)C_rA_P. Then R[P]=sum_j p_j rho_j^(osc).

The four base observation polynomials are P_V=C, P_N=N, P_A=RC, and P_Q=Q, with degrees at most 2,4,3,3.

The frequency-derivative row of P is dot(P)(lambda)=i(lambda-log a)P(lambda), exactly matching delta_a(r).

The base and lifted observation amplitudes have degrees at most 2,4,3,3 and 3,5,4,4. No degree above five is introduced.

Let E[P]=nu c_xi R[P], with the fixed source normalizer and common phase. For X in {V,N,A,Q}, set x_X=Re E[P_X] and z_X=Re E[i(lambda-log a)P_X].

## One Linear Functional

```text
P_lin=N_(p,xi)P_V+V_(p,xi)P_N-Q_(p,xi)P_A-A_(p,xi)P_Q+i(lambda-log a){(N_p+A_(T,x))P_V+(X_T+V_p)P_N-(Q_p+X_(T,x))P_A-(A_T+A_p)P_Q}.
```

The entire linear remainder interaction is exactly g_p^Tr=Re E[P_lin], where deg(P_lin)<=5.

The quadratic remainder is

```text
The quadratic remainder is exactly r^TM_xi r=x_Vz_N+x_Nz_V-x_Az_Q-x_Qz_A. It contains four signed base/derivative pairings and no hidden moment directions.
```

Large componentwise moment errors can be harmless if all eight visible projections vanish. Conversely, smallness of only the six raw moments is not the theorem dictated by the flow.

## Tail Propagation

If P(lambda)=sum_(j=0)^5p_jlambda^j, the known outside-tail bounds imply |T[P]|<2h^2 sum_(j=0)^5 K_j|p_j|, with K=(6,14,55,336,2738,27936).

This is a coefficient-dependent propagation template. The physical coefficients of P_lin and the eight base/lifted rows are not yet bounded here, and the c' interior is untouched.

## Nonpromotion Guard

An exact algebraic specialization has P_lin of degree five with leading coefficient 6i and P_lin(log a)=-21. Thus the compression does not force a universal factor lambda-log a or any automatic endpoint zero. This witness is not asserted to be a physical Xi state.

## Route Decision

Estimate the single grouped interior functional associated with P_lin before taking absolute values. Treat the four quadratic observation pairings at their normalized scale. Do not allocate independent budgets to rho_0,...,rho_5.

## Next Target

Expand the physical coefficients of P_lin from C,N,RC,Q and the retained carrier observations. Derive the exact c'_(P_lin,r) Morse amplitude, audit its discrete mode variation and saddle zeros, and prove or falsify a grouped bound after the common phase, reciprocal pairing, endpoints, and terminal row are retained.

## Proof Boundary

This proves exact observation-image compression, polynomial-amplitude closure through degree five, collapse of the full linear remainder to one complex Morse functional, four exact quadratic pairings, and one coefficient-dependent outside-tail template. It does not evaluate P_lin on the physical chart, bound its c-prime interior, control the four quadratic pairings, prove an h^2 flow error, evaluate the retained carrier, prove a signed flow or Phi_B bound, exclude contact, establish a retained aggregate or Xi theorem, prove Q209, Lambda<=0, PF-infinity, RH, or a prize-level conclusion.

## Reproduce

```text
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_observation_image_compression_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_observation_image_compression_gate.py
```
