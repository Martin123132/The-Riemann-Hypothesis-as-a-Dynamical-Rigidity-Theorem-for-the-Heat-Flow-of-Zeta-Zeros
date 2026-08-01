# Newman Abel-Scalar Shear-Flux Reduction

Date: 2026-07-28

Status: exact linked-point topological reduction and route guard;
not a proof of the Xi pointwise gap, one-turn phase budget,
`Lambda<=0`, PF-infinity, RH, or a Clay-prize result.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_first_order_centered_abel_scalar_shear_flux_reduction.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_abel_scalar_shear_flux_reduction.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_abel_scalar_shear_flux_reduction.py
```

## Pi Provenance

The pi in a^2=x/(4pi)+t/16 and the Xi endpoint remains inherited from the completed-zeta normalization and the Riemann-Siegel saddle. The 2*pi in a winding normalization is the period of exp(i*theta). The pi used in the generic sine guard is only the standard sine period and introduces no new constant into the Xi formulas.

## Exact Shear

```text
On every fixed-N q>=1 chart, put alpha=c*u_N. The exact all-fiber identity is mathsf_A=mathcal_C_N+alpha*mathsf_X.
If |mathsf_X|<=delta_L implies |mathcal_C_N|>A_L+epsilon_term, then (mathsf_X,mathcal_C_N) is nonzero everywhere on that q>=1 arc: inside the band the second coordinate is nonzero, and outside it the first coordinate is nonzero.
For 0<=s<=1 and ell>0 define Gamma_(s,ell)=mathsf_X+i*(mathcal_C_N+s*alpha*mathsf_X)/ell. A zero would force mathsf_X=mathcal_C_N=0. Thus the pointwise gap gives a nonvanishing homotopy from Psi_ell=mathsf_X+i*mathcal_C_N/ell to Gamma_ell=mathsf_X+i*mathsf_A/ell.
[[1,0],[s*alpha/ell,1]] with determinant 1.
Changing ell through positive values is the diagonal orientation-preserving map diag(1,1/ell). Once (mathsf_X,mathcal_C_N) is nonzero, it preserves winding and positive-imaginary-ray intersection signs.
```

## Reduced Proxy

```text
The q>=1 signed integer may therefore be computed from the division-free Abel proxy Psi_ell=mathsf_X+i*mathcal_C_N/ell. This reduction includes W_0=0 and removes the terminal shear c*u_N*mathsf_X from the topological coordinate.
Inside a fixed-N chart, alpha_x=Re(s_*'')*u_N+c/(4*T_0) and (mathcal_C_N)_x=(mathsf_A)_x-alpha_x*mathsf_X-alpha*(mathsf_X)_x. The existing endpoint-complete formulas for (mathsf_X)_x and (mathsf_A)_x retain H_0,H_1,H_2,D_0,D_1 and all endpoint derivatives, so the reduced flux remains an O(N) prefix computation.
partial_s arg(Psi_ell)={ell*mathsf_X*(mathcal_C_N)_s-mathsf_X*mathcal_C_N*ell_s-ell*mathcal_C_N*(mathsf_X)_s}/{ell^2*mathsf_X^2+mathcal_C_N^2}.
At mathsf_X=0 and mathcal_C_N!=0, partial_s arg(Psi_ell)=-ell*(mathsf_X)_s/mathcal_C_N.
```

## Crossing Orientation

On |mathsf_X|<=delta_L, |partial_x mathsf_X-mathcal_C_N| is at most epsilon_term+2*exp(-L)*delta_L+[101*exp(-7L/4)+1e-7*exp(-5L/4)]/|f_1|.

After division by A_L, the four terms are below exp(-3L/2)/(4L), exp(-L)/L, 101*exp(-L/2)/(100000L+1), and 1e-7/(100000L+1). They decrease for L>=50 and sum to less than 2.03e-14 at L=50.

Under the pointwise gap, |mathcal_C_N|>A_L+epsilon_term>A_L exceeds the certified orientation error. Hence sign(partial_x mathsf_X)=sign(mathcal_C_N) throughout the contact band. At mathsf_X=0, mathcal_C_N>0 is exactly an upward crossing and contributes -1 on an increasing-x edge.

## Zero Fibre

At W_0=0 one has mathsf_X=mathsf_Y=0 and Psi_ell=i*mathcal_C_N/ell. The pointwise gap keeps this value nonzero; no ratio, Wronskian, or Schur coordinate is used.

The shear is exact inside each canonical fixed-N chart. The already proved adjacent real-projection homotopies must still be used at cutoff changes; no complex chart-invariance of mathcal_C_N is inferred.

## Backward-Heat Guard

```text
For integer m>=1, X_m(x,t)=exp(m^2*t)*sin(m*x+pi/4) satisfies partial_t X_m=-partial_x^2 X_m. Put mathcal_C_m=partial_x X_m and alpha=0.
At t=0, on |X_m|<=delta<1, |mathcal_C_m|>=m*sqrt(1-delta^2). The pointwise crossing gap can be made arbitrarily large.
On the half-open interval 0<=x<2*pi there are exactly m upward crossings, and wind[X_m+i*mathcal_C_m]=-m. The argument-flux numerator is -m^2*exp(2m^2*t)<0.
```

The backward-heat Fourier family proves that an arbitrarily strong pointwise crossing gap, exact identity mathcal_C=partial_x mathsf_X, and heat evolution do not bound the one-sided crossing count. It is a generic route guard, not an Xi counterexample.

## Minimal Linked-Point Route

```text
Pointwise Xi target: prove |mathsf_X|<=delta_L implies |mathcal_C_N|>A_L+epsilon_term in each prescribed canonical q>=1 chart, retaining W_0=0 and the endpoint.
Independent signed Xi target: compose Psi_ell through the q>=1 horizontal arcs, adjacent cutoff homotopies, vertical connector, finite shoulders, and chart joins, and prove Delta_(partial D_j)arg(Psi_j)<2*pi. Since the transferred integer kappa_j is already nonnegative, this gives 0<=kappa_j<1 and hence kappa_j=0.
The q=2tL^2<1 multiplicity-compatible parabolic/Hermite first-jet chart and its finite connectors remain separate.
```

The smallest surviving linked-point route consists of two independent Xi theorems: the Abel-scalar contact-band gap and the endpoint-complete Abel-proxy phase budget below one turn. Full Schur disk stability and full-circle three-cylinder zero-freeness are stronger than required.

## Boundary

The all-fiber shear identity, determinant-one homotopy, positive-scale invariance, reduced argument flux, crossing orientation handoff, zero-fibre inclusion, and generic backward-heat many-crossing guard are exact or certified. The pointwise Xi Abel-scalar gap and the Xi phase budget below one turn are independent open theorems. No q<1 closure, complete boundary composition, contact exclusion, Lambda<=0, PF-infinity, RH proof, or Clay-prize conclusion is asserted.
