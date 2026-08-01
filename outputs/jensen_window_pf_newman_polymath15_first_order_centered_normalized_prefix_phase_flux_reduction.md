# Newman Normalized-Prefix Phase-Flux Reduction

Date: 2026-07-27

Status: exact normalized-prefix and first-jet flux formulas,
with a certified Xi spiralling law and a decisive recrossing
guard. The successor flux theorem remains open. This is not a proof of RH.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_first_order_centered_normalized_prefix_phase_flux_reduction.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_normalized_prefix_phase_flux_reduction.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_normalized_prefix_phase_flux_reduction.py
```

## Pi Provenance

```text
The pi in a^2=T_0/(2*pi)=x/(4*pi)+t/16 is inherited from xi(s)=s*(s-1)*pi^(-s/2)*Gamma(s/2)*zeta(s)/2 and the Riemann-Siegel saddle. The cutoff width 4*pi*(2N+1) is its algebraic cell difference. A phase cell uses the same constant because exp(i*(theta+2*pi))=exp(i*theta). No circle or prefix polygon is selected.
```

The same universal constant is obtained from every Euclidean
circle, but no circle is selected in this calculation. The
prefix polygon is only a complex partial-sum path.

## Normalized Polygon

```text
Because z_n=eta*q_n, e=eta*r_0, and g=eta*r_A, G_k=conj(eta)*F_k=sum_(n=1)^k q_n with G_0=0 and G_1=1. Thus Z_0=r_0+G_N=conj(eta)*W_0. The prefix polygon is the piecewise-linear path through G_0,G_1,...,G_N; it has no role in defining pi.
Let V_tilde_N=r_A-s_*'*u_N*r_0 and B_N=V_tilde_N+s_*'*sum_(k=1)^(N-1)h_k*G_k. Then Z_A=s_*'*u_N*Z_0+B_N=conj(eta)*W_A. For Pi_eta(w)=Re(eta*w), mathsf_X=Pi_eta(Z_0), mathsf_A=Pi_eta(Z_A), and mathsf_A=Pi_eta(B_N)-b*u_N*Im(eta*Z_0)+c*u_N*mathsf_X. Hence at mathsf_X=0 the normalized formula is exactly mathcal_C_N=mathsf_A, including Z_0=0.
```

The normalization removes a nonzero unit phase; it does not
discard the absolute projector or the endpoint.

## Carrier Currents

```text
On a fixed-N cell put epsilon_n=d_(n,x)/(1+d_n). Then gamma_n=q_(n,x)/q_n=-s_*'*log(n)+epsilon_n-epsilon_1, gamma_1=0, rho_n=Re(gamma_n), vartheta_n=Im(gamma_n), and G_(k,x)=sum_(n=1)^k gamma_n*q_n.
Since |epsilon_n|<8446/x^2, rho_n=-c*log(n)+Re(epsilon_n-epsilon_1). On q>=1, c=t*D_x/4>1/(8L^2*x). Therefore for n>=2, rho_n<-log(2)/(16L^2*x)<0; at L=50 the discarded coefficient-current error is less than 1.5e-14 of this certified inward floor.
For 1<=n<m<=N, vartheta_m-vartheta_n=(-b)*log(m/n)+Im(epsilon_m-epsilon_n)>=log(m/n)/2-16892/x^2>0. In a complete fixed-t cutoff cell, h_n=log((n+1)/n)>1/N and Delta x=4*pi*(2N+1), so every adjacent relative phase advances by more than 4*pi. The worst L=50 integrated error is below 8.3e-20 of the surplus over 4*pi.
```

At the L=50 boundary the inward error-to-floor ratio is
`1.49617241664434464e-14`, and the
adjacent integrated error-to-4pi-surplus ratio is
`8.25272247692305813e-20`.

The prefix currents themselves are

```text
Writing q_n=R_n*exp(i*theta_n), delta_nm=theta_n-theta_m, gamma_n=rho_n+i*vartheta_n, Re(G_(k,x)*conj(G_k))=sum_n rho_n*R_n^2+sum_(n<m)R_n*R_m*[(rho_n+rho_m)*cos(delta_nm)-(vartheta_n-vartheta_m)*sin(delta_nm)]; Im(G_(k,x)*conj(G_k))=sum_n vartheta_n*R_n^2+sum_(n<m)R_n*R_m*[(rho_n-rho_m)*sin(delta_nm)+(vartheta_n+vartheta_m)*cos(delta_nm)].
The diagonal angular terms are positive, but both pair kernels contain unrestricted sine and cosine factors. Strict inward motion and ordered angular currents therefore do not assign a sign to a prefix radial or angular current.
```

## Recrossing Guard

```text
{
  "coefficients": "q_2=exp(-epsilon*theta)*exp(i*theta), q_8=(1/2)*exp(-epsilon*theta)*exp(3*i*theta), epsilon>0",
  "downward_count": 3,
  "downward_zeros": "pi/3,2pi/3,3pi/2",
  "endpoint_relation": "Gamma_epsilon(2pi)=exp(-2pi*epsilon)*Gamma_epsilon(0) for Gamma_epsilon=X_epsilon+i*partial_theta X_epsilon",
  "factorization": "cos(theta)+(1/2)cos(3theta)=cos(theta)*(2cos(theta)^2-1/2)",
  "first_jet_winding": -3,
  "interpretation": "Both carriers move strictly inward, their angular currents are 1 and 3, the second amplitude is half the first, and the relative phase makes two turns. Nevertheless one phase cell has six simple crossings, three upward crossings, zero net signed crossing count, and first-jet winding -3. The simple zeros persist when the faster angular current is increased slightly, so a strict-more-than-two-turn version also exists.",
  "projection": "X_epsilon(theta)=exp(-epsilon*theta)*[cos(theta)+(1/2)cos(3theta)]",
  "signed_count": 0,
  "upward_count": 3,
  "upward_zeros": "pi/2,4pi/3,5pi/3",
  "zero_count": 6,
  "zeros": "pi/3,pi/2,2pi/3,4pi/3,3pi/2,5pi/3"
}
```

This exact model has the same leading logarithmic frequency
ratio `log(8)/log(2)=3` and amplitude ratio
`8^(-1/2)/2^(-1/2)=1/2`. It is a route guard, not an Xi
counterexample: the actual endpoint and full coefficient chain
are deliberately absent.

The first arithmetic completion test gives

```text
{
  "chebyshev_polynomial": "For c=cos(theta), 2sqrt(2)*C_3=4sqrt(2)c^3+4c^2-sqrt(2)c-2=:P(c).",
  "completed_block": "C_3(theta)=cos(theta)+2^(-1/2)cos(2theta)+(1/2)cos(3theta)",
  "completed_upward_count": 1,
  "completed_zero_count": 2,
  "discriminant": "-1696",
  "endpoint_signs": "P(-1)=2-3sqrt(2)<0<P(1)=2+3sqrt(2)",
  "interpretation": "Restoring the arithmetically required intermediate power removes the sparse model's extra two upward crossings. This does not prove the full Xi chain, but it shows that exact multiplicative completion can supply rigidity that generic amplitude and phase ordering miss.",
  "next_hypothesis": "Test complete prime-power chains, then odd-part times dyadic chain decompositions, with the heat-quadratic and d_n perturbations retained. A theorem must also control sums of different chains and the endpoint; blockwise one-turn behavior alone is not promotable.",
  "real_root_count": 1,
  "required_intermediate": "The actual t=0 dyadic chain also contains n=4 with q_4/q_2=2^(-1/2) and phase 2theta.",
  "sparse_block": "cos(theta)+(1/2)cos(3theta), from the leading n=2 and n=8 coefficients after scaling by q_2",
  "unit_interval_root_count": 1
}
```

Thus the generic shortcut is false, but the actual
multiplicative chain repairs its smallest sparse obstruction.
The full prime-power, cross-chain, heat, and endpoint theorem
is still open.

```text
For a regular real scalar X on [A,B], N_up-N_down=[sgn(X(B))-sgn(X(A))]/2. This endpoint flux controls only the net signed count. In the exact recrossing model it is zero while N_up=3.
```

So a signed endpoint flux cannot replace the one-sided count.

## First-Jet Flux

```text
For ell>0 define Gamma_ell=mathsf_X+i*mathsf_A/ell. Where Gamma_ell!=0, partial_x arg(Gamma_ell)={mathsf_X*partial_x(mathsf_A/ell)-(mathsf_A/ell)*partial_x mathsf_X}/{mathsf_X^2+(mathsf_A/ell)^2}. At mathsf_X=0, mathsf_A>0, and partial_x mathsf_X>0, the positive-imaginary-ray intersection has sign -1. Thus the first-jet argument flux, with the prescribed connector and half-open conventions, is the correct object for the upward count.
Let eta_x=i*omega_eta*eta, Z_(0,x)=r_(0,x)+sum_n gamma_n*q_n, V_tilde_(N,x)=r_(A,x)-s_*''*u_N*r_0-s_*'*u_(N,x)*r_0-s_*'*u_N*r_(0,x), u_(N,x)=1/(8*pi*a^2)=1/(4*T_0), and B_(N,x)=V_tilde_(N,x)+s_*''*sum_k h_k*G_k+s_*'*sum_k h_k*sum_(n<=k)gamma_n*q_n. Then Z_(A,x)=(s_*''*u_N+s_*'*u_(N,x))*Z_0+s_*'*u_N*Z_(0,x)+B_(N,x), mathsf_X_x=Pi_eta(Z_(0,x))-omega_eta*Im(eta*Z_0), and mathsf_A_x=Pi_eta(Z_(A,x))-omega_eta*Im(eta*Z_A). This supplies an exact O(N) integrand for the first-jet phase flux.
```

This is the surviving exact O(N) phase-flux representation.

## Live Theorem

```text
Retain the pointwise all-fiber lower bound |mathcal_C_N|>A_L+epsilon_term on |mathsf_X|<=delta_L. For the integer, do not count turns of an individual carrier. Insert the normalized O(N) flux into the complete successor boundary composition and prove 0<=kappa_j=(2*pi)^(-1)*Delta_arg(Gamma_j)<1 after top, bottom, vertical, chart-join, and finite-shoulder terms are combined. The lower inequality is topological; the strict upper inequality is the open arithmetic theorem.
The q=2tL^2<1 layer still requires a separate multiplicity-compatible parabolic/Hermite first-jet flux chart. No uniform positive slope at t=0 is introduced.
```

No Xi prefix lower bound, successor flux upper bound, q<1
chart, finite shoulder closure, contact exclusion,
`Lambda<=0`, PF-infinity, RH, or Clay-prize conclusion is
claimed.
