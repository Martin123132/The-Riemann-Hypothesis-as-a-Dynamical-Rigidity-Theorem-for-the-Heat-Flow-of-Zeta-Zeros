# Complete Prime-Power-Chain Occupation Gate

Date: 2026-07-28

Status: exact chain reduction and route obstruction. This is not a proof of RH.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_first_order_centered_complete_prime_power_chain_occupation_gate.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_complete_prime_power_chain_occupation_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_complete_prime_power_chain_occupation_gate.py
```

## Coefficient Contract

The active coefficients are q_n=exp[t*log(n)^2/4-s_*log(n)]*(1+d_n)/(1+d_1), q_1=1. They carry coefficient +1; there is no Mobius factor in this Xi prefix.

## Exact Chain Algebra

```text
Fix p prime, p does not divide m, ell=log(p), n_k=m*p^k, and R_m=floor(log_p(N/m)). With B_m=phi*exp[t*log(m)^2/4-s_*log(m)]/|1+d_1|, w_(m,p)=exp[ell*(t*log(m)/2-s_*)], and a_(m,k)=exp[t*ell^2*k^2/4]*(1+d_(m*p^k)), one has z_(m*p^k)=eta*q_(m*p^k)=B_m*a_(m,k)*w_(m,p)^k exactly.
Q_(m,p)(w)=sum_(k=0)^R_m a_(m,k)w^k, Z_(m,p)=B_m*Q_(m,p)(w_(m,p)), and K_(m,p)=B_m*w_(m,p)*Q_(m,p)'(w_(m,p))=sum_k k*z_(m*p^k).
Put kappa_0=log(N/m), u_0=log(a/m), C_(m,p)=Re Z_(m,p), and s_*'=c_s+i*b. Then D_(m,p)=Re[(c_s*kappa_0+i*b*u_0)Z_(m,p)-ell*s_*'*K_(m,p)].
On every fixed-N chart and fixed chain, D_(m,p)=partial_x C_(m,p)+c_s*log(N)*C_(m,p)-R_(m,p), where R_(m,p)=sum_k r_(m*p^k). Thus a chain contact still reduces to transversality, not positivity.
```

## Complete-Chain Threshold Guard

Let A_0,A_1,A_2>0 with A_2<A_0, put theta_0=arccos(A_2/A_0), varphi=(pi-theta_0)/2, and theta_k=theta_0+k*varphi. Then theta_2=pi, c_0=A_2, c_2=-A_2, and c_1<0.

For u_0>2ell, u_k=u_0-k*ell, kappa_k=kappa_0-k*ell, beta=-b>=1/2, and 0<=c_s<1/4, h_k=c_s*kappa_k+beta*u_k*tan(theta_k) satisfies h_1<h_2<h_0. Indeed h_1-h_2<c_s*ell-beta*u_1<ell/4-ell/2<0, while h_0-h_2=2c_s*ell+beta*u_0*tan(theta_0)>0.

For P_(m,p)(s)=sum_(k=0)^2 c_k*1_(s<h_k), P=c_1<0 below h_1, P=c_0+c_2=0 on (h_1,h_2), and P=c_0=A_2>0 on (h_2,h_0). Thus neither a fixed sign nor a strict positive lower modulus follows from complete-chain amplitudes and affine phases.

```text
A=(1,2^(-1/2),1/2), theta=(pi/3,2pi/3,pi), u=(5ell/2,3ell/2,ell/2), c_s=0, b=-1
h_1=-3sqrt(3)ell/2<h_2=0<h_0=5sqrt(3)ell/2
-sqrt(2)/4 below h_1; 0 on (h_1,h_2); 1/2 on (h_2,h_0)
ell*(10sqrt(3)+3sqrt(6))/8>0
```

Adding pi to every theta_k leaves every h_k unchanged but sends every c_k, P_(m,p), C_(m,p), and D_(m,p) to its negative. Any composable sign theorem must control the p-free external phase B_m; internal chain contraction cannot do so.

For fixed correction phases delta_k, the two affine parameters alpha,varphi can still set alpha+delta_0=theta_0 and alpha+2varphi+delta_2=pi. The middle phase is (theta_0+pi)/2+delta_1-(delta_0+delta_2)/2. Small d_n corrections therefore do not create a coefficient-only sign; this does not assert that the actual Xi trajectory realizes the guard phases.

## Pi Provenance

The pi in the Xi saddle and a^2=x/(4*pi)+t/16 comes from the completed-zeta Gamma normalization. The pi in theta_2=pi is separately the half-period of exp(i*theta); no polygon or prime chain defines pi.

## P-Free Rejoining

```text
For every fixed prime p, each 1<=n<=N has the unique form n=m*p^k with p not dividing m. Therefore mu=sum_(p not dividing m)mu_(m,p), P=sum_m P_(m,p), C_bulk=sum_m C_(m,p), and D_bulk=sum_m D_(m,p) exactly.
Let M=floor(N/p). The number of p-free complete chains with R_m=0 is S_p(N)=N-2M+floor(M/p). For p=2, S_2(N)>=floor(N/4). Thus a fixed-prime chain split leaves order-N singleton chains with no internal completion at all.
For the correction-free Dirichlet prefix S_M(s)=sum_(n<=M)n^(-s), O_M^(p)(s)=S_M(s)-p^(-s)S_floor(M/p)(s). If R=floor(log_p N), then sum_(k=0)^R p^(-ks)O_floor(N/p^k)^(p)(s)=S_N(s) by exact telescoping.
```

The subtraction in O_M^(p)=S_M-p^(-s)S_floor(M/p) is a one-prime p-free inclusion identity. It does not turn the original carrier coefficients into Mobius coefficients and supplies no positivity.

The Gaussian identity writes the heat factor as an expectation of tilted Dirichlet prefixes. The p-free decomposition telescopes for every tilt before taking the expectation. This exactly reconstructs the original prefix; expectation cannot be passed through the phase-dependent indicators 1_(s<h_n).

## Route Decision

Complete prime-power chains remain a useful exact grouping, but chainwise threshold positivity is false. The p-free external phases, order-N singleton chains, sharp cutoff, and recurrent endpoint must be controlled together. The next target is the endpoint-complete joined-prefix first-jet boundary flux, beginning with the first unresolved compact boundary beyond Q_207.

## Boundary

Proved: exact physical chain factorization, Euler moment, chain transport, three-level threshold guard, external-phase reversal, p-free partition, singleton count, and telescoping rejoin. Open: a joined Xi boundary margin, winding cap, q<1 chart, finite-height cofinal closure, contact exclusion, Lambda<=0, PF-infinity, RH, and any Clay-prize conclusion.
