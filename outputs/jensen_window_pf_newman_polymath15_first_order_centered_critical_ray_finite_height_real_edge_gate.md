# Newman Polymath-15 Critical-Ray Finite-Height Real-Edge Gate

Date: 2026-07-31

Status: effective sign for the retained first-order real-edge model on the full critical range. This is not a proof of an Xi-level edge sign, Lambda <= 0, or RH.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_first_order_centered_critical_ray_finite_height_real_edge_gate.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_critical_ray_finite_height_real_edge_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_critical_ray_finite_height_real_edge_gate.py
```

## Domain Extension

For L>=50 and 0<tL<=25,

```text
0<t<=1/2,          h=1/a<exp(-25)<1/72000000000.
```

The exact stable formulas from the q=1 gate depend on t but do not use q=1. This gate redoes every t-dependent majorant on the enlarged box rather than extrapolating the q=1 constants.

## Enlarged Source Bounds

The two visibly changed inputs are

```text
rho/h^2 <= 1/32,
t/4+t^2|w|^2/8 < 1/2.
```

After the same constant and removable-quotient cancellations, the full normalized chain gives

```text
|alpha'|<h^2,            |alpha''|<h^4,
|d|<h^2/4,               |d_x|<h^4/32,
|W|<4h,                  |W_x|<h^2,
|Z+Q|<8h,                |Z_x-i h theta Q/2|<6h^2,
|mu|<h^2,                |mu_x|<h^4.
```

The checker stores and re-evaluates all 27 normalized rational inequalities. The 4096-box C_0 derivative cover and both factored Cauchy charts are recomputed as well.

## Four-Jet Transfer

The complete normalized edge still satisfies

```text
c=A+h e_0,                       |e_0|<100,
d_edge=hB+h^2 e_1,               |e_1|<250,
c_x=hB+h^2 e_2,                  |e_2|<250,
d_(edge,x)=h^2M+h^3 e_3,         |e_3|<900.
```

Therefore the same division-free perturbation budget holds:

```text
|a^2 J_edge/S_a^2-K_edge(p)|<5500h<1/10000000
a^2 J_edge/S_a^2<-3749/10000000<0
```

uniformly for L>=50, 0<tL<=25, and -1<=p<=1.

## q>=1 Consequence

Since q=2tL^2, every q>=1 point in the critical range is included. Thus the pointwise prescribed-chart q>=1 extension requested after the q=1 gate is closed for the retained first-order model.

## Boundary and Handoff

A pointwise sign on each prescribed-N chart does not itself sign the transition when the cutoff changes. The next gate is the adjacent-cutoff signed splice. The omitted higher-order Xi/source remainder and its insertion into the cumulative contact scalar remain separate. No Xi-level edge sign, Abel-scalar gap, winding cap, contact exclusion, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.
