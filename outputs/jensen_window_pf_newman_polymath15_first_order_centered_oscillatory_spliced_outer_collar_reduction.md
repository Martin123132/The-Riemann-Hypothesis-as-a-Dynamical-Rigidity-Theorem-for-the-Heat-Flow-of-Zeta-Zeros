# Newman Oscillatory-Spliced Outer-Collar Reduction

Date: 2026-07-28

Status: exact conditional domain and degree reduction. This is not a
proof of the Xi Abel gap, `Lambda<=0`, PF-infinity, RH, or a Clay-prize
result.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_first_order_centered_oscillatory_spliced_outer_collar_reduction.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_oscillatory_spliced_outer_collar_reduction.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_oscillatory_spliced_outer_collar_reduction.py
```

## Pi Provenance

The pi in L=log(x/(4*pi)) and in the successor radius R_j=4*pi*exp(L_j) is inherited unchanged from the completed-zeta normalization and the Riemann-Siegel saddle. The splice introduces no circle, polygon, fitted pi, or new normalization.

## Current Frontier

```text
The pinned July-2026 ANTEDB audit, including the four Tao-Trudgian-Yang and two Cushing post-2023 pairs, leaves the exact pointwise-beta contact unchanged at c_*=4911678521/1933561194. Its active contact is alpha_*=62831/155153, beta_*=220633/620612.
```

## Epsilon Collar

```text
Fix 0<epsilon<25-c_*. The oscillatory-zeta theorem gives a finite L_epsilon such that the exact Xi first Laguerre quantity is positive whenever 0<t<=1/2, L>=L_epsilon, and c=tL>=c_*+epsilon. Put B_epsilon=max(50,L_epsilon).
For L_j=101+j, t_b=25/L_j, t_t=25/(L_j-1), and D_j=[0,4*pi*exp(L_j)]x[t_b,t_t], choose J_epsilon=max(0,ceil(B_epsilon)-101). For j>=J_epsilon define L_(epsilon,*)(t)=max(B_epsilon,(2t)^(-1/2)). Then Omega_(j,epsilon)^out is exactly the closed part of D_j with L>=B_epsilon and q=2tL^2>=1. Its complement Omega_(j,epsilon)^in contains q<1 and the fixed L<=B_epsilon shoulder.
The raised interface is continuous and piecewise analytic, with at most one corner at t=1/(2B_epsilon^2). Since L_j>=B_epsilon and sqrt(L_j/50)<L_j, it lies inside both successor horizontal edges. On its q=1 branch, c=1/(2L)<=1/100<c_*. Its L=B_epsilon branch may cross the oscillatory seam, where exact noncontact is already known.
```

For the concrete bookkeeping choice `epsilon=1/100`, the new closed
Abel cap is

```text
c<= 246550706647/96678059700 = 2.550223984760009...
```

instead of `c<=25`. The theorem remains valid for every smaller fixed
positive epsilon, with its own finite threshold.

## Three-Regime Cover

```text
Write c_epsilon=c_*+epsilon and split the raised outer collar into A_(j,epsilon)={c<=c_epsilon}, Z_(j,epsilon)={c_epsilon<=c<=25}, and G_(j,epsilon)={c>=25}. The first set lies in the checked L>=50, q>=1, 0<c<=25 first-order domain. The exact oscillatory-zeta theorem gives noncontact on Z, and the dominant-saddle theorem gives noncontact on G. The closed sets meet at c=c_epsilon and c=25 and cover the collar.
The only open seam tests are therefore q=1, c=c_epsilon, L=B_epsilon, W_0=0, the recurrent endpoint, and adjacent cutoff overlaps. The q=1 seam is strictly inside the low-c Abel region. At c=c_epsilon and c=25 the globally continuous exact Xi jet, rather than an unjoined chart coordinate, is used. The known exact noncontact theorems include both equality seams.
```

## Conditional Excision

```text
Assume only on A_(j,epsilon) that every prescribed canonical chart satisfies |mathsf_X|<=delta_L implies |mathcal_C_N|>A_L+epsilon_term, including W_0=0, endpoints, equality boundaries, and cutoff overlaps. The checked first-order transfer makes the exact Xi jet J_H=(H_t,partial_xH_t) nonzero on A. The oscillatory and dominant theorems make it nonzero on Z and G. Hence J_H is nonzero on all of Omega_(j,epsilon)^out.
The raised inner and outer interface copies have opposite orientations, so partial D_j=partial Omega_(j,epsilon)^in+partial Omega_(j,epsilon)^out. Conditional outer noncontact gives deg(J_H,Omega_(j,epsilon)^out,0)=0 and therefore kappa_j=wind(J_H(partial Omega_(j,epsilon)^in),0) for every j>=J_epsilon.
```

## Exact Consequence

```text
For each fixed admissible epsilon, the sufficient outer Abel theorem is needed only on the closed wedge L>=B_epsilon, q>=1, 0<c<=4911678521/1933561194+epsilon. Requiring it through c=25 is a valid but unnecessary strengthening.
Because epsilon is arbitrary, every fixed scaled ray c>c_* is eventually closed. Thus the asymptotically live outer arithmetic envelope is c<=c_*+o(1). This statement does not manufacture a uniform explicit epsilon(L), and it does not include the q<1 parabolic layer.
```

## Route Guards

```text
The current ANTEDB pointwise-beta machinery cannot lower this seam: the exact TY1/TY2 contact survives the six post-2023 pairs and twelve audited beta-to-beta iterations. A lower seam needs an improved beta bound near (alpha_*,beta_*) or cancellation beyond pointwise beta majorants.

L_epsilon is existential in the imported theorem. The logical reduction is valid with that finite constant, but a finite-height computer certificate needs an effective upper bound for B_epsilon. Until then, L<B_epsilon and the finitely many early successor slabs remain explicit inner or finite-shoulder obligations; they may not be silently declared checked.
```

## Route Decision

Use the oscillatory-spliced collar as the primary outer route. Attack the Abel scalar only in the low-c wedge and develop the multiplicity-compatible q<1/bounded-L inner degree theorem separately. Retain the c<=25 collar and the H_j<3*pi/2 boundary ledger as stronger fallbacks.

## Open Theorems

Outer:

```text
For one fixed 0<epsilon<25-c_* and its B_epsilon=max(50,L_epsilon), prove uniformly on L>=B_epsilon, q=2tL^2>=1, 0<tL<=c_*+epsilon that |mathsf_X|<=delta_L implies |mathcal_C_N|>A_L+epsilon_term. Include q=1, W_0=0, the recurrent endpoint, equality boundaries, and every adjacent-cutoff overlap.
```

Inner:

```text
Construct one multiplicity-compatible exact or rigorously dominated successor proxy on q<1 together with the fixed L<=B_epsilon shoulder, core, early slabs, and all joins. Prove zero degree, or boundary nonvanishing and winding strictly below one. Effectivize B_epsilon before claiming a finite exhaustive certificate.
```

## Boundary

This artifact proves the epsilon-quantified raised-collar geometry, three-regime cover, exact seam allocation, conditional degree excision, and the reduction of the outer Abel burden from c<=25 to c<=c_*+epsilon. It does not prove the Xi Abel-scalar gap, an effective L_epsilon, the q<1/bounded-L inner theorem, complete boundary nonvanishing, contact exclusion, Lambda<=0, PF-infinity, RH, or a Clay-prize conclusion.
