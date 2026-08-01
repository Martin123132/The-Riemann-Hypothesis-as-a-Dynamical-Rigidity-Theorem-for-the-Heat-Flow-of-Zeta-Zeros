# First-Order Rectangular Boundary-Degree Reduction

Date: 2026-07-29

Status: exact boundary and degree reduction with open Xi Abel and
winding targets. This is not a proof of contact exclusion,
Lambda<=0, or RH.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_first_order_rectangular_boundary_degree_reduction.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_rectangular_boundary_degree_reduction.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_rectangular_boundary_degree_reduction.py
```

Current result:

```text
validated Newman first-order rectangular boundary-degree reduction: 22 rows, 2 error-box coordinates, 2 rectangular homotopies, 1 positive-index degree transfer, 1 boundary-only Abel target, 2 route branches, 2 nonpromotion guards, 0 contact exclusions
```

## Error-Matched Boundary Geometry

The first-order remainder is componentwise:

```text
|r_[1]|<B_0(L), |partial_x r_[1]/L|<B_1(L)
```

Define

```text
M_L(J)=max(|J_[1]|/B_0(L),|partial_x J_[1]/L|/B_1(L))
M_L(J)=G_L(X,U)
```

Then `M_L(J)>1` is exactly the box-optimal signed-band condition.
On a boundary it gives the componentwise Rouche homotopy

```text
V_s=V_J+s*V_r, M_L(V_s)>=M_L(V_J)-s*M_L(V_r)>0
```

so the first-order main and exact normalized Xi jet have the same
boundary winding. The older Euclidean target remains sufficient but
is strictly stronger.

## Error Whitening

```text
D_L=diag(B_0(L)^(-1),B_1(L)^(-1))
det(D_L)=1/(B_0(L)*B_1(L))>0
D_L maps the certified error box into (-1,1)^2
```

The positive diagonal homotopy preserves winding. This packages the
unequal C1 value and derivative errors without replacing their
rectangle by a larger ellipse.

## Abel Boundary Target

The division-free target is unchanged but is required only on the
selected boundary arcs:

```text
|mathsf_X|<=delta_L => |mathcal_C_N|>A_L+epsilon_term
U=u_N*S+sum_(k=1)^(N-1)log((k+1)/k)*F_k
```

The target implies M_L(J)>1 on the selected boundary, including W_0=0.

Adjacent charts already lie in a strict sub-box:

```text
|Delta J_[1]|/B_0(L)<1/5
|partial_x Delta J_[1]/L|/B_1(L)<3/20
M_L(J+s*Delta J)>1-s/5>=4/5
```

## Degree Transfer

The positive normalizer shear preserves degree, and every real heat
contact has local index `floor(m/2)>0` in the standard `(x,t)`
orientation. Consequently boundary winding zero excludes all
contacts. On one successor the weaker integer trap is enough:

```text
0<=kappa_j=wind(V_J)<1 => kappa_j=0
```

Boundary nonvanishing alone is not enough: `F=x^2-2t` has a
boundary-nonzero rectangle with one positive-index contact.

## Hard Bottom Interval

On the ray schedule,

```text
t_j=25/(100+j), L_j=101+j, R_j=4*pi*exp(L_j)
q=50*L^2/(100+j)
t_j*L=25*L/(100+j)
```

the remaining low-c, q>=1 bottom interval is exactly

```text
max(B_epsilon,sqrt((100+j)/50))<=L<=min(L_j,(c_*+epsilon)*(100+j)/25)
```

Oscillatory-zeta closes the raised middle and the dominant theorem
closes t_jL>=25. q<1 and L<B_epsilon remain separate.

## Route Fork

Route A:

Prove the Abel gap on the full two-dimensional raised outer wedge; excise it at zero degree.

Route B:

Prove the Abel gap only on joined successor boundaries and H_j<3*pi/2; use C_j<pi/2 and the integer trap.

Recommended next bounded target:

Attempt the one-dimensional hard-bottom Abel theorem first, preserving the full joined phase ledger.

Open handoff:

For every successor j and fixed epsilon, prove the normalized Abel gap on each nonempty hard bottom interval max(B_epsilon,sqrt((100+j)/50))<=L<=min(L_j,(c_*+epsilon)(100+j)/25), retaining W_0=0, the recurrent endpoint, q=1, equality boundaries, and adjacent-cutoff homotopies. Then prove the completely joined horizontal phase contribution H_j<3*pi/2. If the same Abel estimate extends to the full raised outer wedge, use the uniform-collar excision route instead.

Pi provenance is unchanged: `pi` comes from the completed-zeta,
Riemann-Siegel, Fourier, and phase-period normalizations already
recorded in the source chain.

## Boundary

This artifact proves the rectangular first-order boundary homotopy, error-whitening and winding invariance, adjacent-chart sub-box join, positive-index degree transfer, exact hard-bottom interval, and the two-route comparison. It does not prove the Xi Abel gap, horizontal phase bound, q<1 or bounded-L inner theorem, complete boundary nonvanishing, contact exclusion, Lambda<=0, PF-infinity, RH, or a Clay-prize conclusion.
