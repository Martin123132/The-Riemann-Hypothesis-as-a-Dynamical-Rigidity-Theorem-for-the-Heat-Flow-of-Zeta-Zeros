# Newman Joined Dyadic Odd-Prefix First-Jet Reduction

Date: 2026-07-27

Status: exact joined reduction and hypothetical homotopy-base
theorem; not a proof of `Lambda<=0`, RH, PF-infinity, or a
Clay-prize result.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_first_order_centered_joined_dyadic_odd_prefix_first_jet_reduction.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_joined_dyadic_odd_prefix_first_jet_reduction.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_joined_dyadic_odd_prefix_first_jet_reduction.py
```

## Pi Provenance

The pi in a^2=x/(4*pi)+t/16 is inherited from xi(s)=s*(s-1)*pi^(-s/2)*Gamma(s/2)*zeta(s)/2 and the Riemann-Siegel saddle. The cutoff-cell difference is 4*pi*(2N+1), while 2*pi is the period of exp(i*theta). The dyadic base 2, odd prefixes, phase lift, and phase-synchronized hypothetical do not define pi.

## Deterministic Heat Shift

```text
For A_t(n;s)=exp[t*log(n)^2/4-s*log(n)], h=log(2), A_t(2^k*m;s)=A_t(2^k;s)*A_t(m;s-t*k*h/2).
```

This factorization is stronger than merely writing the heat
quadratic as a Gaussian expectation: each dyadic layer is an
exactly shifted odd prefix before any averaging.

## Joined Odd Prefixes

```text
Put A_t(n;s)=exp[t*log(n)^2/4-s*log(n)], c_n=(1+d_n)/(1+d_1), epsilon_n=d_(n,x)/(1+d_n), delta_n=epsilon_n-epsilon_1, and q_n=A_t(n;s_*)c_n. Then q_(n,x)/q_n=-s_*'*log(n)+delta_n.
Let h=log(2), K=floor(log_2 N), M_k=floor(N/2^k), and for j>=0 define O_(k,j)=sum_(m<=M_k,m odd)(log m)^j*A_t(m;s_*-t*k*h/2)c_(2^k*m). Then, for H_r=sum_(n<=N)(log n)^r q_n, H_r=sum_(k=0)^K A_t(2^k;s_*)*sum_(j=0)^r binom(r,j)(k*h)^(r-j)O_(k,j). This is an exact deterministic heat-shifted odd-prefix identity, not a conditional or averaged approximation.
Define E_(k,j) by replacing c_(2^k*m) in O_(k,j) with delta_(2^k*m)c_(2^k*m), and put D_r=sum_k A_t(2^k;s_*)*sum_(j=0)^r binom(r,j)(k*h)^(r-j)E_(k,j). Then D_r=sum_(n<=N)(log n)^r delta_n q_n exactly.
```

Unique dyadic valuation is doing the bookkeeping. At a cutoff
jump, exactly one odd-prefix layer gains the entering integer.

## Phase Lift

```text
Write sigma_*=Re(s_*), omega_*=-Im(s_*), a_(k,m)=exp[t*log(2^k*m)^2/4-sigma_*log(2^k*m)]>0, and define P_N(z,xi,mu)=sum_(k=0)^K z^k*sum_(m<=M_k,m odd)a_(k,m)exp(i*xi*log m)*(1+mu*d_(2^k*m))/(1+mu*d_1). For z_*=exp(i*omega_*h), P_N(z_*,omega_*,1)=H_0. With D_h=h*z*partial_z-i*partial_xi, D_h^r P_N(z_*,omega_*,1)=H_r.
For D_h=h*z*partial_z-i*partial_xi, D_h[z^k*exp(i*xi*log(m))]=log(2^k*m)*z^k*exp(i*xi*log(m)); D_h^r supplies the r-th total logarithmic moment. Also partial_mu P_N is obtained termwise from partial_mu[(1+mu*d_n)/(1+mu*d_1)]=(d_n-d_1)/(1+mu*d_1)^2. The physical fixed-chart derivative is H_(0,x)=-s_*'*D_h P_N+D_0 at the actual point.
```

The variables `z` and `xi` separate the dyadic phase from the
odd Mellin phase. The physical Xi prefix sits on the linked
curve `z=exp(i*omega_*log(2))`, `xi=omega_*`.

## Five-Current First Jet

```text
The fixed-chart logarithmic moments obey H_(r,x)=-s_*'*H_(r+1)+D_r. Hence the endpoint-complete first jet closes on H_0,H_1,H_2,D_0,D_1 together with r_0,r_A and their x derivatives; no individual carrier phase or Gaussian fibre has to be discarded.
Let A=log(a), u_x=A_x=1/(8*pi*a^2), P_1=H_1-AH_0, Z_0=r_0+H_0, and Z_A=r_A-s_*'P_1. Then Z_(0,x)=r_(0,x)-s_*'H_1+D_0, P_(1,x)=-s_*'(H_2-AH_1)+D_1-AD_0-u_xH_0, and Z_(A,x)=r_(A,x)-s_*''P_1-s_*'P_(1,x). These formulas include the endpoint before any absolute value or winding estimate is taken.
For eta_x=i*omega_eta*eta and Pi_eta(w)=Re(eta*w), mathsf_X=Pi_eta(Z_0), mathsf_A=Pi_eta(Z_A), mathsf_X_x=Pi_eta(Z_(0,x))-omega_eta*Im(eta*Z_0), and mathsf_A_x=Pi_eta(Z_(A,x))-omega_eta*Im(eta*Z_A). For ell>0, Gamma_ell=mathsf_X+i*mathsf_A/ell and partial_x arg Gamma_ell=[mathsf_X*partial_x(mathsf_A/ell)-(mathsf_A/ell)*mathsf_X_x]/[mathsf_X^2+(mathsf_A/ell)^2].
```

Thus the exact first jet needs only `H_0,H_1,H_2,D_0,D_1`
plus endpoint value and derivative data. This is a reduction
of organization, not a sign theorem.

## Cutoff And Endpoint

```text
At N->N+1 write n=N+1=2^v*m with m odd. Exactly one odd-prefix layer, v, gains one term, so Delta H_r=(log n)^r q_n and Delta D_r=(log n)^r delta_n q_n. With j_0=kappa_N*J_a/f_1 and j_A=kappa_N*(J_(a,x)+mu_a*J_a)/f_1, Delta Z_0=q_n+j_0=Q_N/f_1 and Delta Z_A=-s_*'*log(n/a)q_n+j_A. Differentiating gives Delta Z_(0,x)=gamma_nq_n+(j_0)_x and Delta Z_(A,x)=(j_A)_x-s_*''log(n/a)q_n-s_*'[-u_xq_n+log(n/a)gamma_nq_n].
The imported exact adjacent recurrence is therefore the cutoff boundary law of the joined representation, not a separate approximation. Its real projections retain |Delta X|<1100exp(-5L/4) and |Delta A_a|<5000exp(-7L/4); no complex jump bound is invented.
```

The terminal recurrence is part of the joined representation.
It is not appended after a bulk winding calculation.

## Phase-Synchronized Hypothetical

The deliberately deformed case `xi=mu=0` can be solved exactly:

```text
Set xi=0 and mu=0 while keeping z free. Then P_fr(z)=sum_(k=0)^K B_kz^k with B_k=sum_(m<=M_k,m odd)a_(k,m)>0. For every m retained in layer k+1, the certified correction-free dyadic ratio is below 2^(-49/100), and the odd index set shrinks. Thus 0<B_(k+1)<2^(-49/100)B_k. Enestrom-Kakeya puts every zero outside |z|=2^(49/100)>1. Therefore P_fr is unit-disk zero-free with winding zero, and zP_fr has winding one. This is a proved theorem for the stipulated phase-synchronized hypothetical, not for the actual Xi odd phases.
If rho_1,...,rho_K are the roots of P_fr, then |rho_j|>R_0=2^(49/100). Hence on |z|=1, |P_fr(z)|=B_K*product_j|z-rho_j|>B_K*(R_0-1)^K. This explicit margin is positive but may be too small for the full odd-phase homotopy.
```

This is the useful thought experiment. Synchronizing all odd
phases produces a zero-free joined polynomial. The physical
problem is therefore an explicit transport question from that
proved base to the actual Mellin phases, rather than an
undefined appeal to multiplicativity.

## Homotopy Target

```text
The exact deformation from the proved base to the actual finite prefix is (xi,mu):(0,0)->(omega_*,1) in P_N, with z=z_*. A viable new theorem must exclude P_N(z_*,xi,mu)=0 along a specified path or control every boundary crossing using the D_h first and second moments. The endpoint r_0 and centered jet r_A-s_*'P_1 must then be inserted before the successor argument is counted. Crude L1 phase perturbation is not assumed to be smaller than the phase-frozen margin.
```

A successful theorem may require exponential-polynomial
Hermite-Biehler/de Branges machinery, a signed two-parameter
Jacobian, or a linked-curve argument-principle estimate.

## Route Guards

### Averaging does not preserve winding

```json
{
  "equal_weight_mean": "(A+B)/2",
  "fibre_windings": "wind(A)=wind(B)=1 on |z|=1",
  "first_fibre": "A(z)=z*(1+r*z+r^2*z^2)",
  "interpretation": "Averaging boundary-nonzero winding-one fibres can produce winding three. A per-Gaussian-fibre theorem cannot be averaged without additional Gaussian/Xi structure.",
  "mean_inner_discriminant": "-19/2",
  "mean_inner_root_modulus": "sqrt(2/5)<1",
  "mean_winding": 3,
  "second_fibre": "B(z)=-(4/5)*z*(1+r*z)"
}
```

### The second logarithmic moment is indispensable

```json
{
  "interpretation": "For unrestricted joined coefficients, value and first logarithmic moment do not determine the second moment entering the first-jet derivative. Xi-specific structure must control H_2 rather than delete it.",
  "moments": "H_0=0, H_1=0, H_2=2h^2",
  "support": "logarithmic nodes 0,h,2h",
  "weights": "1,-2,1"
}
```

These are generic exact guards, not Xi counterexamples.

## Live Theorem

```text
On q=2tL^2>=1, prove an Xi-specific transport theorem from the phase-frozen base to (z_*,omega_*,1), strong enough to yield |mathsf_X|<=delta_L => |mathcal_C_N|>A_L+epsilon_term and the complete composed successor bound 0<=kappa_j<1. Candidate mechanisms include a Hermite-Biehler/de Branges exponential-polynomial criterion, a signed two-parameter Jacobian, or a Rouche/argument-principle estimate that uses the linked dyadic and odd phases rather than independent block windings.
```

The separate small-`q` obligation is

```text
The q=2tL^2<1 layer remains a separate multiplicity-compatible parabolic/Hermite first-jet chart. The joined dyadic identities remain valid there, but no uniform endpoint simplicity or positive slope is inferred.
```

## Boundary

The heat-shift factorization, joined odd-prefix moments, five-current derivative closure, phase lift, cutoff jump, endpoint first jet, phase-frozen zero-free theorem and margin, and two generic route guards are exact or imported from certified prerequisites. The actual odd-phase/correction homotopy, joined Xi lower bound, strict successor flux upper bound, q<1 chart, finite shoulders and connectors, contact exclusion, Lambda<=0, PF-infinity, RH, and a Clay-prize conclusion remain open.
