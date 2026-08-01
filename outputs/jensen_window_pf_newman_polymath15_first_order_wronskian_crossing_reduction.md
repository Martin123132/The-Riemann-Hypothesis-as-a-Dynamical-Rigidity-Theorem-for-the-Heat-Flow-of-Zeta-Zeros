# Newman First-Order Wronskian Crossing Reduction

Date: 2026-07-26

Status: exact Wronskian and oriented-crossing reduction. The
Xi arithmetic theorem remains open; this is not a proof of RH.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_first_order_wronskian_crossing_reduction.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_wronskian_crossing_reduction.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_wronskian_crossing_reduction.py
```

## First-Order Wronskian

```text
E_[1]=X+iY, E_[1],x=U+iV, W_[1]=Im(E_[1],x*conj(E_[1]))
W_[1]=V*X-U*Y
J_[1]=2X, J_[1],x=2U
At a full contact, X=-r_[1]/2 and U=-r_[1],x/2.
2*|W_[1]|<exp(-5L/4)*(100000*|E_[1],x|+200000*L*|E_[1]|)
For every point in the live layer, prove either |X|>50000*exp(-5L/4) or 2*|W_[1]|>exp(-5L/4)*(100000*|E_[1],x|+200000*L*|E_[1]|).
```

## Boundary Crossings

```text
At X=0, W_[1]=-U*Y.
At X=0 and E_[1]!=0, X^2+(U/L)^2>12500000000*exp(-5L/2) iff |W_[1]|>50000*sqrt(5)*L*|E_[1]|*exp(-5L/4).
Where E_[1]!=0, theta_[1],x=W_[1]/|E_[1]|^2.
At X=0, E_[1]=iY!=0, an upward zero U>0 is equivalent to W_[1]*Y<0.
If E_[1]=0 then X=Y=W_[1]=0; a proxy crossing with U>0 is not classified by W_[1]*Y and must be counted separately.
N_up=N_(X=0,E_[1]!=0,W_[1]*Y<0)+N_(E_[1]=0,U>0), with half-open endpoint conventions.
```

A Wronskian magnitude bound supplies separation, but the winding
also needs signed crossings and the explicit `E_[1]=0` class.

## Successor Composition

```text
Insert the decomposed N_up on both horizontal edges into kappa_j=N_up(t_j;R_j)-N_up(t_(j+1);R_j)+I_i(V_j), then add the finite shoulder/join intersection cells.
It is enough to prove the resulting oriented count is <1; positive contact degree then forces kappa_j=0.
The absolute Wronskian margin proves boundary nonvanishing. The signed W_[1]*Y counts and exceptional complex-main zeros are additional data needed for the one-sided winding bound.
```

## Route Guards

```text
The aggregate phase derivative is not known to have one sign; stored corrected crossings exhibit both signs below the asymptotic L>=50 layer.
The instantaneous component-rate Hermitian form has inertia (1,1,m-2) when rates differ, so it is not positive semidefinite.
Positive amplitudes and ordered same-sign component speeds admit an exact double-crossing countermodel.
Finite sign diagnostics shape the route but neither prove nor disprove the L>=50 heat-coupled theorem.
```

## Remaining Xi Input

```text
On q=2tL^2>=1 in the live 0<tL<c_*+o(1) layer, prove the first-order Wronskian disjunction and a signed crossing/phase budget including E_[1]=0 exceptional crossings.
Treat q<1 by a multiplicity-compatible Hermite or degree chart; do not demand a fixed phase-speed floor at t=0.
Close bounded L, L_epsilon, vertical connector, and chart-join phase cells separately.
```

## Live Handoff

```text
Use the Wronskian magnitude to certify first-jet separation at real-part crossings, and use the signs of W_[1]*Im(E_[1]) plus the explicit E_[1]=0 count to upper-bound the successor winding. Do not replace this signed arithmetic theorem by phase monotonicity or component-rate ordering.
```

This artifact proves the first-order contact Wronskian bound, crossing-margin equivalence, upward-sign identity, exceptional-zero split, and successor-count composition. It does not prove the q>=1 Wronskian theorem, q<1 multiplicity-compatible theorem, finite shoulders, one-sided successor winding, contact exclusion, Lambda<=0, RH, PF-infinity, or a Clay-prize conclusion.
