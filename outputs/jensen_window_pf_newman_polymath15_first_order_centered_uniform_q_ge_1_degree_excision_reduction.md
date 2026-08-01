# Newman Uniform-q>=1 Degree-Excision Reduction

Date: 2026-07-28

Status: exact conditional outer-collar degree reduction. This is not a
proof of the uniform Xi Abel gap, the inner successor theorem,
`Lambda<=0`, PF-infinity, RH, or a Clay-prize result.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_first_order_centered_uniform_q_ge_1_degree_excision_reduction.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_uniform_q_ge_1_degree_excision_reduction.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_uniform_q_ge_1_degree_excision_reduction.py
```

## Pi Provenance

The pi in L=log(x/(4*pi)) and X_*(t)=4*pi*exp(L_*(t)) is the same universal constant already fixed by the completed-zeta normalization and Riemann-Siegel saddle. The 2*pi in winding is the period of exp(i*theta). The excision introduces no circle, polygon, fitted pi, or new normalization.

## Successor Collar

```text
For L_j=101+j put t_b=25/L_j, t_t=25/(L_j-1), R_j=4*pi*exp(L_j), and D_j={(x,t):0<=x<=R_j, t_b<=t<=t_t}, with the standard (x,t) orientation.
For L(x)=log(x/(4*pi)), define L_*(t)=max(50,(2t)^(-1/2)) and X_*(t)=4*pi*exp(L_*(t)). Then Omega_j^out={(x,t) in D_j:x>=X_*(t)} is exactly the closed part of D_j with L>=50 and q=2tL^2>=1.
Let Omega_j^in be the closure of D_j\Omega_j^out. It contains the q<1 layer, every L<50 bounded/core shoulder, and their finite joins. The common interface is I_j={(X_*(t),t):t_b<=t<=t_t}.
At t_b and t_t the q=1 logarithmic thresholds are sqrt(L_j/50) and sqrt((L_j-1)/50), respectively. Both are strictly below L_j for L_j>=101, so the outer collar is nonempty on both horizontal edges.
The interface is continuous and piecewise analytic. Its only possible corner is where (2t)^(-1/2)=50, namely t=1/5000. Thus Omega_j^in and Omega_j^out are compact Lipschitz domains even when one successor slab contains that corner.
```

## Two-Regime Cover

```text
Split Omega_j^out into A_j={tL<=25} and G_j={tL>=25}. On A_j the still-open uniform Abel-scalar gap lies in the certified L>=50, q>=1, 0<tL<=25 first-order remainder domain. On G_j the already proved dominant-saddle theorem gives L_t(x)>0 for the exact Xi flow. The two closed pieces meet at tL=25 and cover the entire outer collar.
The whole moving interface lies strictly on the tL<25 side: if L_*=50 then tL_*=50t<=25/2, while on q=1 one has tL_*=sqrt(t/2)<1. Thus no dominant/Abel seam is placed on the degree-excision interface.
```

The first-order remainder is used only on `tL<=25`. The already proved
dominant-saddle theorem supplies exact noncontact on `tL>=25`; the two
pieces meet at equality.

## Conditional Noncontact

```text
Assume the closed-domain Xi estimate |mathsf_X|<=delta_L implies |mathcal_C_N|>A_L+epsilon_term in every prescribed canonical chart throughout A_j, including q=1, tL=25, W_0=0, endpoints, and cutoff overlaps. The checked shear/scale and first-order remainder transfers then keep the exact jet J_H=(H_t,partial_x H_t) nonzero on A_j. The dominant theorem keeps J_H nonzero on G_j. Hence J_H is nonzero on all of Omega_j^out.
mathcal_C_N is used only inside a canonical Abel chart to certify nonvanishing of the exact jet. It is not declared a global complex coordinate. Degree and interface cancellation are performed with the globally continuous exact Xi jet J_H; adjacent-cutoff and positive-normalizer homotopies retain their existing endpoint and overlap obligations.
```

## Oriented Excision

```text
Orient the outer boundary as bottom outer left-to-right, right edge bottom-to-top, top outer right-to-left, and I_j top-to-bottom. Orient the inner boundary as bottom inner left-to-right, I_j bottom-to-top, top inner right-to-left, and the symmetry axis top-to-bottom. The two copies of I_j cancel, so partial D_j=partial Omega_j^in+partial Omega_j^out as oriented one-chains.
If the conditional outer noncontact theorem holds and the remaining inner boundary contract makes all displayed windings defined, then deg(J_H,Omega_j^out,0)=wind(J_H(partial Omega_j^out),0)=0. Therefore kappa_j=wind(J_H(partial D_j),0)=wind(J_H(partial Omega_j^in),0). All q>=1 and dominant-ray horizontal turns cancel through the zero-degree outer collar; no separate signed outer phase bound is needed.
```

## Route Guards

```text
For integer m>=1, X_m=exp(m^2t)sin(mx+pi/4) and C_m=partial_xX_m give a globally nonzero exact backward-heat jet. On [0,2*pi] each forward horizontal path winds -m, but the reversed upper path winds +m and the two vertical paths cancel. Thus the full zero-free strip has degree zero for every m. Arbitrarily many turns on one horizontal edge do not obstruct two-dimensional excision.

F(t,x)=x^2-2t satisfies F_t=-F_xx. On [-1,1]x[-1/4,1/4], J_F=(F,F_x) is boundary-nonzero but has the interior zero (x,t)=(0,0), whose Jacobian determinant in (x,t) coordinates is +4. Hence its boundary winding is +1. Boundary nonvanishing alone cannot justify outer excision; nonvanishing on the whole closed collar is essential.
```

The first guard shows why arbitrarily many turns on one horizontal edge do
not obstruct a zero-free two-dimensional collar. The second shows why
boundary-only nonvanishing is too weak for the same conclusion.

## Route Fork

```text
Boundary-only route: retain the proved connector cap |C_j|<pi/2 and prove the joined complementary inequality H_j<3*pi/2. Uniform-collar route: prove the Abel gap on all of A_j, combine it with dominant noncontact on G_j, excise Omega_j^out, and prove only the inner-region successor degree theorem. The second route asks for a stronger pointwise domain statement but removes the independent many-turn outer phase theorem.
Promote the uniform-collar route as primary because the Abel target is already pointwise and its formulas are uniform in (x,t) on the source domain. Retain the connector quarter-turn theorem and H_j<3*pi/2 ledger as a valid fallback if the gap can be proved only on boundary arcs.
```

## Open Theorems

Outer:

```text
Prove uniformly on the closed set L>=50, q=2tL^2>=1, 0<tL<=25 that |mathsf_X|<=delta_L implies |mathcal_C_N|>A_L+epsilon_term, with equality boundaries, W_0=0, the recurrent endpoint, and every cutoff overlap included. This is the only new outer arithmetic theorem.
```

Inner:

```text
On Omega_j^in, construct one multiplicity-compatible exact or rigorously dominated proxy covering q<1, L<50, the core, finite shoulders, and all joins. Prove boundary nonvanishing and an oriented winding below one (or direct zero degree) for every j. Positive local contact degree then forces the localized successor integer to vanish.
```

## Boundary

The successor geometry, moving q=1/L=50 interface, two-regime cover, oriented chain cancellation, conditional degree excision, paired many-turn cancellation guard, boundary-only contact guard, and route comparison are exact or inherited from checked sources. This does not prove the uniform Xi Abel-scalar gap, its exact-H noncontact consequence, the inner q<1/bounded-L successor theorem, complete boundary nonvanishing, contact exclusion, Lambda<=0, PF-infinity, RH, or a Clay-prize conclusion.
