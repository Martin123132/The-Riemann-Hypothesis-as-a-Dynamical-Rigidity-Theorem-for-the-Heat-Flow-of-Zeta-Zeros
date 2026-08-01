# Endpoint-Composed Six-Moment Saddle Flow

Date: 2026-08-01

Status: exact saddle-flow reduction with `0 signed flow bounds` and `0 Phi_B
bounds`. This is not a proof of RH.

## Auxiliary Current

```text
For T_0-epsilon<=xi<=T_0 put z_n(xi)=omega_a A_n exp(-i*xi*u_n). Then z_(n,xi)=-i*u_n*z_n and conj(z_n)_(xi)=i*u_n*conj(z_n).

Define Psi(xi)=Re sum_(n<=B)L_nz_n(xi)+(1/2)Re sum_(n,m<=B)[H^+_(n,m)z_n(xi)conj(z_m(xi))+T^+_(n,m)z_n(xi)z_m(xi)]. At xi=Omega, Psi(Omega)=h^2*Phi_B exactly.
```

Only xi=Omega is the physical carrier state. Intermediate xi and xi=T_0 form an exact auxiliary homotopy and are not asserted to satisfy the terminal source recurrence separately.

## Exact Flow

```text
Psi'(xi)=Re sum_n[-i*u_n*L_n]z_n+(1/2)Re sum_(n,m){[-i(u_n-u_m)H^+_(n,m)]z_nconj(z_m)+[-i(u_n+u_m)T^+_(n,m)]z_nz_m}.

h^2*Phi_B-Psi(T_0)=-epsilon*integral_0^1 Psi'(T_0-theta*epsilon)dtheta.
```

If H^+_(m,n)=conj(H^+_(n,m)), then -i(u_n-u_m)H^+_(n,m) is Hermitian. If T^+_(m,n)=T^+_(n,m), then -i(u_n+u_m)T^+_(n,m) is transpose-symmetric. The complete saddle transfer remains one endpoint-linear, one phase-difference, and one phase-sum family. No new correlation family or four-variable expansion is created.

## Leading Pair Kernel

```text
For M_0=f+ig, M_1=p_1+iq_1, M_2=p_2+iq_2, and q_3=Im(M_3), put K=q_1p_2-fq_3+2u_x(q_1g-fp_1)=4P_bulk^(0)'(xi).

With Delta_nm=u_n-u_m and Sigma_nm=u_n+u_m, K=-(1/4)sum_(n,m)[Delta_nm^2 Sigma_nm Im(z_nz_m)+Delta_nm Sigma_nm^2 Im(z_nconj(z_m))]-u_x sum_(n,m)Sigma_nm Re(z_nz_m).

The Hermitian and cubic transpose pieces vanish on n=m. The exact diagonal is K_diag=-2u_x sum_n u_n Re(z_n^2), which has no fixed sign.
```

This pair identity is exactly the xi derivative of the ideal kernels H_0^+=-Sigma_nm^2/8 and T_0^+=-Delta_nm^2/8-i*u_x/2. The identity preserves both physical phase families. Its diagonal coefficient is small, but no diagonal or aggregate sign is inferred without the actual phases.

## Six-Moment Closure

```text
In ell_n=log(n), C,D have degree at most 2, R,R_x degree at most 1, Q degree at most 3, N and L degree at most 4. Multiplication by u_n, u_n-u_m, or u_n+u_m raises each one-variable degree to at most 5.

Therefore Psi'(xi) closes exactly on the six correction-free moments G_0(xi),...,G_5(xi), with no seventh moment and no new correction denominator.
```

The exact degree audit is

```json
{
  "L_n": 4,
  "u_n_L_n": 5,
  "H_plus_ell_n": 4,
  "H_plus_ell_m": 4,
  "flow_H_ell_n": 5,
  "flow_H_ell_m": 5,
  "flow_T_ell_n": 5,
  "flow_T_ell_m": 5
}
```

## Order-Five Mangoldt Form

```text
Put c_xi=exp[i(Omega-xi)log(a)], s_xi=sigma-i*xi, and q_n^(0)(xi)=exp[t log(n)^2/4-s_xi log(n)]. Then z_n(xi)=(eta/S_a)c_xi q_n^(0)(xi). The factor c_xi is common to every carrier and stays outside each divisor sum.

G_5(xi)=(eta/S_a)c_xi sum_(dm<=B)Lambda(d)log(dm)^4q_(dm)^(0)(xi)=(eta*c_xi/(2S_a))sum_(dm<=B)[Lambda(d)+Lambda(m)]log(dm)^4q_(dm)^(0)(xi).

The exact factor q_(dm)^(0)=q_d^(0)q_m^(0)exp[(t/2)log(d)log(m)] is unchanged; xi only changes the common logarithmic phase.

For D=floor(sqrt(B)), the order-five sum has the same exact D-by-D square and two hyperbolic wings as orders one through four.

For every prime sqrt(B)<p<=B, the order-five symmetric {1,p} minor is [[0,log(p)^5/2],[log(p)^5/2,0]] and has determinant -log(p)^10/4<0.
```

The formal audit checks `205` one-sided
and `205` symmetric rows.

## Exact Target

```text
Prove directly Psi(T_0)-epsilon integral_0^1Psi'(T_0-theta*epsilon)dtheta<=h^2/400, preferably h^2/800, for the actual endpoint and six-moment Xi vector.
```

Apply a uniform signed Type-I/II, Vaughan, or reciprocal-Poisson estimate to the endpoint-linear, Hermitian flow, and transpose flow families together. The new arithmetic input is only the order-five logarithmic moment; no new correlation geometry is needed.

## Pi Provenance

No new `pi` is introduced. The existing `u_x=h^2/(8pi)` and saddle phase
inherit the completed-zeta, Riemann--Siegel, and Poisson normalizations
recorded in Formal Core Section 11.118.

## Boundary

This proves the exact endpoint-composed frequency-flow and transfer identities, preservation of Hermitian and transpose symmetry, the leading pair-flow formula, degree-five closure, and the order-five Mangoldt/balanced-hyperbola extension. It does not prove a signed flow estimate, saddle-proxy upper bound, Phi_B upper bound, contact exclusion, retained aggregate or Xi-level current theorem, Q209, Lambda<=0, PF-infinity, RH, or a prize-level result.

## Reproduce

```text
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_saddle_flow_reduction.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_saddle_flow_reduction.py
```
