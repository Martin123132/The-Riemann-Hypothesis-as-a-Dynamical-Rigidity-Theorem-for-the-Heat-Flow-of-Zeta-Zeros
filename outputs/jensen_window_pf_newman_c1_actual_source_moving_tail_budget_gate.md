# Newman C1 Actual-Source Moving-Tail Budget Gate

Date: 2026-08-04

Status: the complete observation-level moving-tail package is below rho A_T/1000000; the nonterminal quadratic package and signed main remain open; not a proof of RH.

## Ownership and definition

The all-carrier transfer differs from the frozen-terminal transfer only through the integral of D_move; this observation-level term is outside the later relative-lift and quadratic packages and is inserted exactly once.

R_move=-varepsilon integral_0^1D_move(T_0-theta varepsilon)dtheta and E_move=2R_move/|nu|^2, where varepsilon=T_0-Omega.

## Source-free envelopes

For q_n=exp[t(log n)^2/4-sigma log n], q=1 gives sigma>1/2 and t(log n)^2/4<1/8, hence q_n<(8/7)n^(-1/2).

For every 1<=n<=N, |P_V(log n)|,|P_mathcalN(log n)|,|P_A(log n)|,|P_Q(log n)|<P_*:=R_*^2/2, R_*=(L+3)/2.

S_0=sum_(n<=N)q_n<(16/7)h^(-1/2).

S_1=sum_(n<=N)u_nq_n<(16/7)R_*h^(-1/2).

For T_1=sum_(B<n<=N)u_nq_n, K^18<=a and K>=2 give n>=a/2, q_n<2h^(1/2), and T_1<4h^(3/2)K^2.

Uniformly on the auxiliary segment, ||omega_B+omega_C||_infinity<=|nu|P_*S_0, ||dot(omega_B)||_infinity<=|nu|P_*S_1, ||dot(omega_C)||_infinity<=|nu|P_*T_1, and ||Delta omega_C||_infinity<=varepsilon|nu|P_*T_1.

## Common-unit budget

|D_move|<4varepsilon|nu|^2P_*^2T_1S_1+4|nu|^2P_*^2S_0T_1+(12/5)|nu|P_*T_1.

|E_move|<8varepsilon P_*^2(varepsilon T_1S_1+S_0T_1)+(24/5)varepsilon P_*(B_0T_0)T_1.

Using varepsilon<7h^2/(8L^2), |E_move|<16h^3R_*^4K^2/L^2+14h^5R_*^5K^2/L^4+(42/5)(B_0T_0)h^(7/2)R_*^2K^2/L^2.

Using K^2<=h^(-1/9), B_0T_0/rho<61/100, rho>(329/3300)Lexp(-L/4), A_T>14300, and exp(L/4)<h^(-1/2), the three ratios are below [52800/(329*14300)]h^(43/18)R_*^4/L^3, [46200/(329*14300)]h^(79/18)R_*^5/L^5, and [2562/(500*14300)]h^(61/18)R_*^2/L^2.

After h<exp(-L/2), the logarithmic derivatives are -43/36+4/(L+3)-3/L, -79/36+5/(L+3)-5/L, and -61/36+2/(L+3)-2/L; all are negative for L>=50.

At L=50, exp(-25)<1/72000000000 follows from the 61-term positive Taylor floor for exp(25). Replacing powers 43/18,79/18,61/18 by 2,4,3 gives an exact rational sum below 1/1000000.

|E_move|<rho A_T/1000000.

|E_R+E_term+E_aa+E_Delta+E_aff^C+E_move|<(140379929/429000000)rho A_T.

If |E_quad|<=delta_1rho A_T, adverse-witness negativity requires E_FAN<-[280040071/429000000-delta_1]rho A_T.

## Gate rows

| id | role | state | claim |
|---|---|---|---|
| mtb_01_ownership | package ownership | proved | The moving-tail defect is a single observation-level transfer package. |
| mtb_02_definition | exact definition | proved | The package has a fixed sign and scalar normalization. |
| mtb_03_weight | source-free weight | proved | Every correction-free carrier weight has a square-root majorant. |
| mtb_04_observation | observation box | proved | All four carrier observation polynomials share one source-free envelope. |
| mtb_05_mass | all-carrier mass | proved | The complete carrier mass is O(h^(-1/2)). |
| mtb_06_moment | all-carrier moment | proved | The complete first distance moment is O(R_*h^(-1/2)). |
| mtb_07_terminal | terminal moment | proved | The moving terminal derivative has an h^(3/2)K^2 envelope. |
| mtb_08_displacement | terminal displacement | proved | The terminal observation displacement gains one more varepsilon. |
| mtb_09_bilinear | J bilinear | proved | The four-coordinate current has a dimension-fixed sup-norm bound. |
| mtb_10_pointwise | pointwise defect | proved | Both moving-tail terms have complete pointwise source envelopes. |
| mtb_11_integrated | integrated defect | proved | Physical integration and source normalization give three explicit positive terms. |
| mtb_12_main | carrier ratio | proved | The leading all-carrier/terminal interaction is anchor-negligible. |
| mtb_13_delta | displacement ratio | proved | The terminal-displacement interaction is anchor-negligible. |
| mtb_14_edge | edge ratio | proved | The genuine-edge/moving-terminal interaction is source-normalized and anchor-negligible. |
| mtb_15_close | moving budget | proved | The complete moving-tail package costs less than one millionth of the anchor. |
| mtb_16_target | live theorem | open | Only the nonterminal quadratic package remains in the absolute-secondary frontier. |

## Next action

Retain the exact four nonterminal quadratic pairings and their Hermitian/transpose phase structure. Derive a common-unit E_quad bound without componentwise row-norm promotion; then prove or falsify E_FAN<-[280040071/429000000-delta_1]rho A_T on every adverse actual-source witness arc.

## Pi provenance

No new geometric pi is introduced. The only pi-dependent rate used in the observation box is u_(N,x)=h^2/(8pi), inherited from the completed-zeta/Riemann-Siegel saddle coordinate; pi>3 supplies a rational upper bound.

## Proof boundary

This gate proves a fixed-q1 source-free observation envelope and the common-unit budget |E_move|<rho A_T/1000000 for the complete observation-level moving-tail defect. It proves no nonterminal quadratic budget, signed finite-near-affine inequality, complete-current sign, all-q transport, contact exclusion, Q209 shell, cofinal successor, Lambda<=0, PF-infinity, RH, or prize-level conclusion.

built Newman C1 actual-source moving-tail budget gate: 16 rows, 0 issues, |E_move|<rho*A_T/1000000, 1 remaining absolute secondary package, 1 live signed near-affine target
