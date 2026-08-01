# Newman Polymath-15 First-Order Cofinal Boundary Reduction

Date: 2026-07-26

Status: exact boundary reduction with a certified first-order
jet budget. This is not a proof of the cofinal boundary theorem,
`Lambda <= 0`, or RH.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_first_order_cofinal_boundary_reduction.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_cofinal_boundary_reduction.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_cofinal_boundary_reduction.py
```

## First-Order Boundary Budget

```text
V_r=(r_[1],r_[1],x/L), ||V_r||_2^2<50000000000*exp(-5L/2)
If V_J=(J_[1],J_[1],x/L) satisfies J_[1]^2+(J_[1],x/L)^2>50000000000*exp(-5L/2) pointwise on a closed boundary, then V_J+s*V_r is nonzero there for 0<=s<=1 and wind(V_Z)=wind(V_J)
For J_[1]=2X_[1] and J_[1],x=2U_[1], it is enough that X_[1]^2+(U_[1]/L)^2>12500000000*exp(-5L/2)
```

The comparison with the old corrected budget is

```text
Relative to the old 32000000*exp(-3L/2) squared budget, the new/old ratio is (3125/2)*exp(-L), decreasing in L; at L>=50 it is <5*10^-19
```

The exact rational `L=50` upper audit is `156250000000000000000000000000000000000000000000000000/369988485035126972924700782451696644186473100389722973815184405301748249`
or `4.223104402429321692e-19`.

## Ray-Aligned Boundary

```text
t_j=25/(100+j), L_j=101+j, R_j=4*pi*exp(L_j), P_j=[t_j,1/2]x[0,R_j]
On x=R_j and t_j<=t<=1/2, tL_j>=25*(101+j)/(100+j)>25; the dominant-saddle theorem makes this boundary contact-free
For S_j=[t_(j+1),1/2]x[R_j,R_(j+1)], tL>=t_(j+1)L_j=25, so every new strip is already closed
```

The compact core and top edge are inherited from the certified
base and published positive-time bound. The only new arithmetic
work is on the bottom boundary and finite joins.

## Bottom Split

```text
On the bottom t=t_j, q=2t_jL^2=50L^2/(100+j); q<1 iff L<sqrt((100+j)/50), and q>=1 beyond that point
On a fixed bottom, c=t_jL; every fixed c>4911678521/1933561194 is asymptotically closed after the oscillatory-zeta handoff, subject to its epsilon-dependent finite shoulder
The remaining proof work is one-dimensional but has three pieces: the q>=1 critical bottom arc, a multiplicity-compatible q<1 bottom arc, and finite bounded-L/epsilon shoulders, together with the winding computation
```

This reduces the hard two-dimensional wedge to selected
one-dimensional boundary arcs, without claiming those arcs are
already controlled.

## Continuity And Winding

```text
On a doubled cutoff collar, the two adjacent first-order lifts satisfy |Delta J_[1]|<20000*exp(-5L/4)
Cauchy on radius 1/L and |partial_x log A_t|<L/2 give |partial_x Delta J_[1]|<30000*L*exp(-5L/4)
||V_(J,N+1)-V_(J,N)||_2^2<1300000000*exp(-5L/2)=(13/500)*50000000000*exp(-5L/2)
If the standard boundary margin ||V_(J,N)||_2^2>50000000000*exp(-5L/2) holds, then V_(J,N)+s*(V_(J,N+1)-V_(J,N)) is nonzero for 0<=s<=1; the local cutoff lifts therefore join continuously without a separate arithmetic theorem
Boundary nonvanishing transfers winding but does not determine it; because the transferred contact degree is a nonnegative integer, it is enough to prove wind(V_J)<1 on every successor boundary
The ray-aligned P_0 has zero winding because t_0=1/4>1/5>=Lambda; Q208 remains a compact calibration. For every successor, prove the first-order boundary domination and continuous one-sided phase bound only on its boundary; the same-sign positive contact index in the standard (x,t) orientation then excludes all interior contacts
```

The cutoff joins now cost no independent theorem once the standard
margin holds. Positivity of the contact degree means the resulting
integer needs only a strict upper bound below one.

## Proof Boundary

This artifact proves the first-order vector remainder budget, boundary-Rouche transfer condition, exponential improvement over the old budget, ray-schedule right-edge closure, and exact bottom q split. It also proves that the same boundary margin joins adjacent cutoff lifts by a nonvanishing overlap homotopy. It does not prove the q>=1 boundary margin, q<1 multiplicity-compatible margin, one-sided proxy phase bound, bounded-L shoulders, the all-stage cofinal boundary theorem, contact exclusion, Lambda<=0, RH, PF-infinity, or a Clay-prize result.
