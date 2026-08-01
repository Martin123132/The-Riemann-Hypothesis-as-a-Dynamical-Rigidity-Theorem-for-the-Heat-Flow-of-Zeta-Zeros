# Signed Occupation Transport Reduction

Date: 2026-07-28

Status: exact transport, normalization, residual budget, and
transversality equivalence; not a proof of the Xi threshold-mass
bound, contact exclusion, `Lambda<=0`, PF-infinity, RH, or a
Clay-prize result.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_first_order_centered_signed_occupation_transport_reduction.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_signed_occupation_transport_reduction.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_signed_occupation_transport_reduction.py
```

## Pi Provenance

The pi in a^2=x/(4*pi)+t/16 and u_x=1/(8*pi*a^2) remains the completed-zeta and Riemann-Siegel saddle constant. The occupation variable s, Dirac masses, and integrating factor N^(sigma_*) introduce no circle, fitted pi, or geometric normalization.

## Distributional Transport

Fix one N=floor(a) chart on L>=50, 0<tL<=25, and q=2*t*L^2>=1. The distributional slope law is local to a pole-free interior projective chart. The summed division-free carrier identity applies to every 1<=n<=N, including the terminal carrier.

```text
On a pole-free fixed-N chart let I be a fixed set of finite interior slopes and define mu(s;x)=sum_(i in I)c_i*delta(s-h_i), F(s;x)=sum_(i in I)c_i*h_(i,x)*delta(s-h_i), and S(s;x)=sum_(i in I)c_(i,x)*delta(s-h_i).
For every compactly supported smooth test phi, partial_x sum_i c_i*phi(h_i)-sum_i c_i*h_(i,x)*phi'(h_i)=sum_i c_(i,x)*phi(h_i). Equivalently, partial_x mu+partial_s F=S in distributions.
Put P(s;x)=sum_(i in I)c_i*1_(s<h_i). Then partial_x P=F+sum_(i in I)c_(i,x)*1_(s<h_i) as a distribution in s.
```

## Xi Source Collapse

```text
Put epsilon_i=d_(i,x)/(1+d_i)=xi_i+i*chi_i. For z_i=eta*q_i=X_i+iY_i, varrho_i=-c_s*log(i)+xi_i-xi_1 and nu_i=b*u_i+e_i with e_i=v_a+chi_i. Hence X_(i,x)=varrho_i*X_i-nu_i*Y_i.
Define r_i=(xi_i-xi_1)X_i-e_iY_i. Since k_i=log(N/i), the exact division-free identity is X_(i,x)=d_i-c_s*log(N)*X_i+r_i, or d_i=X_(i,x)+c_s*log(N)*X_i-r_i.
Where X_i!=0 and h_i=d_i/X_i, X_(i,x)=[h_i-c_s*log(N)]X_i+r_i. Therefore S=[s-c_s*log(N)]mu+R, where R(s;x)=sum_i r_i*delta(s-h_i).
```

## Common Normalization

```text
Let Omega_N=N^(sigma_*), so Omega_(N,x)=c_s*log(N)*Omega_N. With bar_mu=Omega_N*mu, bar_F=Omega_N*F, and bar_R=Omega_N*R, the exact normalized law is partial_x bar_mu+partial_s bar_F=s*bar_mu+bar_R.
Let M_>(s;x)=sum_i d_i*1_(s<h_i) and R_>(s;x)=sum_i r_i*1_(s<h_i). Then partial_x P=F+M_>-c_s*log(N)*P+R_>, so partial_x(Omega_N*P)=Omega_N*(F+M_>+R_>).
```

The common N^(sigma_*) factor removes the shared frame drift exactly, but the universal reaction s*bar_mu and the signed flux bar_F remain. This normalization is not a positivity transformation.

## Residual Budget

The source bounds |epsilon_i|<8446/x^2 and |e_i|<8449/x^2 give |r_i|<18888*|z_i|/x^2, because 16892^2+8449^2<18888^2.

The corrected adjacent amplitude ratio obeys |z_(n+1)|/|z_n|<exp[-(2/5)log((n+1)/n)]. Since |z_1|=1, |z_n|<n^(-2/5), and sum_(n=1)^N|z_n|<(5/3)N^(3/5)<(10/3)exp(3L/10).

```text
Using x^2>144*exp(2L), |R_bulk|:=|sum_(n=1)^N r_n|<438*exp(-17L/10)<10^(-7)*exp(-5L/4) for L>=50. At L=50 the last inequality follows after squaring from (19/7)^45>(438*10^7)^2 and e>19/7; its ratio then decreases with L.
```

## First Moment

```text
Summing the division-free carrier identity over 1<=n<=N, including the terminal carrier, gives D_bulk=(C_bulk)_x+c_s*log(N)*C_bulk-R_bulk, and Omega_N*D_bulk=partial_x(Omega_N*C_bulk)-Omega_N*R_bulk.
Put E_0=d_0-c_(0,x)-c_s*log(N)c_0. Since mathsf_X=c_0+C_bulk and mathcal_C_N=d_0+D_bulk, exactly mathcal_C_N=partial_x mathsf_X+c_s*log(N)mathsf_X+E_0-R_bulk=Omega_N^(-1)partial_x(Omega_N*mathsf_X)+E_0-R_bulk.
```

R_bulk is uniformly tiny, but E_0 is not assigned a separate sign or magnitude. Only the endpoint-complete combination E_0-R_bulk is compared with the anchored first-jet identity; the terminal carrier remains included division-free in the bulk sum.

The anchored identity and mathcal_C_N=mathsf_A-c_s*u_N*mathsf_X give mathcal_C_N=partial_x mathsf_X+(xi_1-u_a-c_s*u_N)mathsf_X+v_a*mathsf_Y-Re(D_(1,x))/|f_1|. Since log(a)=log(N)+u_N, comparison yields E_0-R_bulk=(xi_1-u_a-c_s*log(a))mathsf_X+v_a*mathsf_Y-Re(D_(1,x))/|f_1|.

On |mathsf_X|<=delta_L the existing endpoint-complete bounds make |partial_x mathsf_X-mathcal_C_N|/A_L less than 2.03e-14. Thus the occupation first moment recovers the certified crossing-orientation handoff, conditional on the still-open pointwise gap.

## Route Guards

The strict theorem h_(i,x)<0 does not make F a one-sign measure because c_i is signed. With c_1=1,c_2=-1, the negative rates (-1,-2) give total flux +1, while (-2,-1) give -1. These are exact local transport jets, not asserted Xi states. In the correction-free normalized local law take (c_1,c_2)=(1,-1), h_(1,x)=h_(2,x)=-1, and zero residual. The slope pairs (1,0), (0,1), and (1,1) all have contact mass c_1+c_2=0 but first moments +1,-1,0. Individual clockwise motion therefore does not sign the normalized first moment. These are local transport guards, not Xi carriers.

If c_i->0 while d_i->d_*!=0, then h_i=d_i/c_i escapes to projective infinity. The compactly supported measure c_i*delta_(h_i) disappears locally while its first moment c_i*h_i=d_i has the finite limit d_*. That boundary mass is D_perp and must be retained division-free or joined in the reciprocal chart.

The existing exact backward-heat family has mathcal_C=partial_x mathsf_X and an arbitrarily strong contact-band gap, yet has any prescribed number of upward crossings. Therefore even a successful pointwise occupation first-moment bound would not by itself prove the one-turn successor theorem.

## Route Decision

Promote the normalized occupation law as an exact organizational lemma, but do not treat it as new coercivity. Its correction density is negligible, while its common frame drift cancels exactly; nevertheless its flux remains signed, and its first moment is precisely the already-open real transversality coordinate. A future gain must prove an Xi-specific signed threshold-mass estimate, most plausibly using complete multiplicative chains and the recurrent endpoint, or return to the endpoint-complete phase-flux boundary contract. Generic transport monotonicity cannot supply either theorem.

The q<1 layer remains a separate multiplicity-compatible parabolic/Hermite or boundary-degree theorem. The affine occupation chart is not continued through t=0 multiplicities by assumption.

## Boundary

This proves the weak and cumulative signed occupation transport laws, exact Xi source collapse, common N^(sigma_*) normalization, quantitative summed residual bound, division-free bulk first moment, endpoint-complete equivalence with the known real transversality identity, signed-flux guards, and projective-infinity boundary term. It does not prove a signed Xi threshold-mass estimate, the pointwise Abel-scalar gap, an edge sign, a one-turn successor bound, q<1 closure, finite-height effectivity, contact exclusion, Lambda<=0, PF-infinity, RH, or a Clay-prize conclusion.
