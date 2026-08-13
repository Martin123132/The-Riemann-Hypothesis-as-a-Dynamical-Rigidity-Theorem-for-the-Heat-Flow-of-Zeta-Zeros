# Reciprocal Terminal Full-Support Homotopy Gate

Date: 2026-08-02

Status: exact all-carrier homotopy, terminal reciprocal strip, and outer-gap reduction proved; edge-composed signed bound open; not a proof of RH.

This is not a proof of RH. It changes the auxiliary transform so the terminal carriers and the bulk are one moving block, while the genuine endpoint edge remains explicit.

## Current Reassembly

omega=(V,mathcal N,A,Q) and Q(omega)=omega^T J omega=Vmathcal N-AQ.

omega_T=omega_E+omega_C(Omega), where omega_C is the carrier tail B<n<=N and omega_E is the genuine endpoint edge.

Psi_B=Q(omega_B)+2omega_T^T J omega_B and Psi_B(Omega)=h^2Phi_B.

tilde(Psi)(xi)=Q(omega_E+omega_B(xi)+omega_C(xi)).

tilde(Psi)(Omega)=Q(omega_T+omega_B(Omega))=h^2P_ret.

tilde(Psi)'=2(omega_E+omega_B+omega_C)^T J(dot(omega_B)+dot(omega_C)).

The exact comparison with the previous frozen-terminal flow is

```text
tilde(Psi)'-Psi_B'=2[omega_C(xi)-omega_C(Omega)]^TJdot(omega_B)+2[omega_E+omega_B+omega_C(xi)]^TJdot(omega_C).
```

## Full-Support Observation Image

G_j^[N]=G_j^[B]+G_j^(C), 0<=j<=5, with G_j^(C) summed over B<n<=N under the same auxiliary rotation.

For Y=y_B+y_C, tilde(Psi)'=Y^TM_xiY+ell_E^TY, where ell_E=U_1(2Jomega_E).

The terminal carriers now live in Y; only the genuine endpoint edge remains in the external linear row.

For Y=p+r, the linear remainder is one degree-at-most-five P_lin^[N]=F_alpha+i(lambda-log a)F_beta with beta=(mathcal N_p+mathcal N_E,V_p+V_E,-Q_p-Q_E,-A_p-A_E).

The nonlinear remainder remains epsilon_Vepsilon_(N,xi)+epsilon_Nepsilon_(V,xi)-epsilon_Aepsilon_(Q,xi)-epsilon_Qepsilon_(A,xi).

## Terminal Interval

H^[N][A]=H^[B][A]+sum_(q=B+1)^N A(q)e(alpha_Plog q).

M_N[A]-M_B[A]=(1/2)A(B)e_B+sum_(q=B+1)^(N-1)A(q)e_q+(1/2)A(N)e_N.

I_(A,N)(r)=I_(A,B)(r)+integral_B^N A(u)e(alpha_Plog u-ru)du.

The B endpoint half is transferred into the terminal interval; it must not be retained simultaneously as a hard endpoint of both transforms.

## Reciprocal Terminal Strip

The cells 1<=q<=N tile (alpha_P/(N+1/2),2alpha_P].

The cells B+1<=q<=N tile (alpha_P/(N+1/2),alpha_P/(B+1/2)].

Their diagonal inversions return q=B+1,...,N-1 with full weight and q=N with half weight; the upper hard endpoint supplies the other half at N.

The full-support scalar transform is still Fourier inversion. Its value is the alignment of the terminal strip and endpoint edge, not scalar smallness.

## Outer Complement

For r<=alpha_P/(N+1/2) and 1<=u<=N, alpha_P/u-r>=alpha_P/[2N(N+1/2)]>1/4 on the physical segment.

For r>2alpha_P and 1<=u<=N, r-alpha_P/u>alpha_P.

For n=N-m in the terminal block, u_n<2hK and |z_n(xi)-z_n(Omega)|<7h^3K|z_n(Omega)|/(4L^2)<h^3K|z_n(Omega)|.

These gaps prove nonstationarity only. They do not bound the endpoint-composed outer complement at h^2 scale.

## Handoff

The next theorem must estimate the full-support observation kernel, not a scalar reciprocal cell. Compose the true endpoint edge with the two outer complement pieces, retain the terminal-strip carriers, the endpoint-linear, Hermitian, and transpose channels, and the four quadratic residual pairings, and only then take a modulus. The direct target is

```text
tilde(Psi)(T_0)-epsilon integral_0^1
  tilde(Psi)'(T_0-theta epsilon)dtheta < 0.
```

## Pi Provenance

`alpha_P=xi/(2pi)` and `e(x)=exp(2pi i x)` come from the fixed completed-zeta and Fourier conventions. No circle, polygon, or fitted geometric constant is added.

## Proof Boundary

This gate proves exact current polarization, the all-carrier physical anchor, six-moment closure on `1<=n<=N`, edge-only external linearization, terminal interval half-weight transfer, terminal reciprocal-strip tiling, two outer phase-gradient gaps, and a terminal phase-motion bound. It proves no edge-complement estimate, signed full-current theorem, quadratic residual bound, `Phi_B` upper bound, contact exclusion, retained aggregate or Xi theorem, `Q209`, cofinal descendant theorem, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion.
