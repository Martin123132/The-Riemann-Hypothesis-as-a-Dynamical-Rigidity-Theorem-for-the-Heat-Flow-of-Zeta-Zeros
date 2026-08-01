# Ray-Aligned Energy Cell-Localization Gate

Date: 2026-07-29

Status: exact localization and endpoint-ownership audit. This
is not a proof of contact exclusion, Lambda<=0, or RH.

```text
work/rh_compute/results/jensen_window_pf_newman_ray_aligned_energy_cell_localization_gate.json
python work/rh_compute/scripts/jensen_window_pf_newman_ray_aligned_energy_cell_localization_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_ray_aligned_energy_cell_localization_gate.py
```

Current result:

```text
validated ray-aligned energy cell-localization gate: 20 rows, 2 exact local balances, 2 exact contact floors, 2 conditional no-contact criteria, 2 outer anchors, 0 cofinal internal anchor registries, 1 finite successor instance, 2 nonpromotion guards, 0 pointwise contact exclusions
```

## Cell Identity

For

```text
38=x_0<x_1<...<x_N=R_j
I_k=[x_(k-1),x_k], 1<=k<=N
h_k=x_k-x_(k-1)>0
```

the fixed-slice balance gives

```text
D_k(t)=d/dt integral_(x_(k-1))^(x_k) E_pf dx+2[J_pf(t,x_k)-J_pf(t,x_(k-1))]
D_k(t)=2 integral_(I_k) B_pf(t,x) dx>=0
```

A contact in `I_k` therefore forces

```text
V_pf(t,c)=0 in I_k => D_k(t)>=2[e_(k-1)(t)+e_k(t)]^2/h_k
V_pf(t,c)=0 in I_k => D_k(t)>=T_k(t), T_k=2[m_(k-1)(t)+m_k(t)]^2/h_k
m_(k-1)>0, m_k=0 => T_k=2m_(k-1)^2/h_k; the right-sided analogue is identical
```

The `1/h_k` gain is real. It is useful only when at least one
endpoint of that cell has an independently certified positive
current-time margin.

## Telescoping

Summing the local identities gives

```text
sum_(k=1)^N D_k(t)=d/dt integral_38^R_j E_pf dx+2[J_pf(t,R_j)-J_pf(t,38)]=2 integral_38^R_j B_pf dx=:D_total(t)
```

Every internal current J_pf(t,x_k), 1<=k<N, appears once with sign + and once with sign -.

The exact summed criterion is

```text
If T_k(t)>0 for every k and D_total(t)<min_(1<=k<=N) T_k(t), then V_pf(t,x) has no zero on [38,R_j].
```

A contact in I_r gives D_r>=T_r, while D_total=sum_k D_k>=D_r because every D_k>=0.

The tempting replacement of the minimum by a sum is false:

```text
T_1=1, T_2=100, D_1=1, D_2=0 gives D_total=1<T_1+T_2=101 while D_1>=T_1.
```

This is an exact implication counterexample for the budget logic, not a claimed realization by the Xi field.

## Endpoint Ownership

The positive-margin nodes form a vertex-cover problem on the
partition path:

```text
{x_0=38,x_N=R_j}
The two owned outer nodes cover both cells only when N<=2. For N>=3 at least one interior edge has neither endpoint owned; bounded-width localization therefore needs a cofinal interior anchor skeleton.
```

Source audit:

- Compact: At x=38 and 0<t<=1/5, the compact certificate gives |J_t'(38)|>20577/350 and the inherited margin sqrt(E_pf(t,38))>=m_(L,j)>0.
- Dominant: The dominant theorem proves strict normalized Laguerre curvature on tL>=25. The endpoint gate exposes the quantitative margin at x=R_j. Uniformly over a descendant collar, the dominant region shrinks to x=R_j at the bottom time, so it does not populate a bounded-cell skeleton through the old interior.
- Oscillatory: For each epsilon>0 the oscillatory theorem gives an existential L_epsilon above c_*+epsilon, c_*=4911678521/1933561194, but no practical L_epsilon or effective cofinal endpoint-margin registry.
- Critical C1: For L>=50 and 0<tL<25, the corrected global C1 remainder satisfies r^2+(r'/L)^2<32000000 exp(-3L/2), but the finite-main lower bound T_L[J]>32000000 exp(-3L/2) remains open.
- Finite phase cells: The Q208 top and bottom certificates own fixed paths t=1/5 and t=1/1040 on 0<=x<=246. The Q207-to-Q208 transport owns one finite successor on 38<=x<=245.
- Conditional successor: The all-stage phase-cell lemma would transfer margins if delta_j M_(j,k)<d_(j,k) on every cell, but that Xi antecedent is precisely open.

The finite inventory is exact:

```text
Q208 bottom: t=1/1040, 492 cells, 0 unresolved
Q208 top: t=1/5, 492 cells, 0 unresolved
Q207-Q208 transport: 1 successor, 525 panels
```

These are finite certificates, not a cofinal current-time
anchor registry.

## Exact Tradeoff

Keeping the thresholds cellwise preserves their 1/h_k strength but also keeps the internal fluxes in D_k.

Summing cancels every internal flux but compresses all contact information to the weakest threshold min_k T_k.

Telescoping is exact bookkeeping, not a substitute for current-time interior endpoint ownership.

## Route Decision

Retain the local and globally telescoped criteria as exact conditional lemmas, but park the energy-localization route. The audited sources do not supply a cofinal current-time interior anchor vertex cover, and attempting to transport one reintroduces the open all-j phase-cell estimate. Return to the critical contact-normal C1 arithmetic target.

Open handoff:

Prove the corrected finite-main inequality T_L[J]>32000000 exp(-3L/2) for L>=50 and 0<tL<25, with cutoff transitions retained. First derive a cutoff-stable phase/amplitude decomposition of T_L[J] and identify the smallest signed cancellation statement that would force the lower bound; do not divide by the unknown first-jet norm.

Pi provenance is unchanged: every occurrence comes from the
established coordinate `L=log(x/(4*pi))`, where `pi` is the
ordinary circle constant.

## Boundary

This artifact proves the cell-local balance, local contact floors, one-sided threshold, exact internal-flux telescoping, global minimum-threshold criterion, vertex-cover ownership condition, and source ownership audit. It does not prove a cofinal interior anchor registry, any strict Xi local or global bulk upper budget, the corrected finite-main C1 lower bound, Q209, the cofinal descendant theorem, Lambda<=0, PF-infinity, RH, or a Clay-prize conclusion.
