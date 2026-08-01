# Adjacent-Cutoff Real-Edge Projective Splice Gate

Date: 2026-07-31

Status: effective signed adjacent-cutoff splice for the retained first-order
real-edge model. This is not a proof of an Xi-level edge theorem, an Abel
gap, `Lambda <= 0`, or RH. It is not an Xi-level result.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_first_order_centered_adjacent_cutoff_real_edge_projective_splice_gate.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_adjacent_cutoff_real_edge_projective_splice_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_adjacent_cutoff_real_edge_projective_splice_gate.py
```

## Cutoff Charts

L>=50, 0<tL<=25, h=1/a<1/72000000000; at an integer cutoff a=m compare p_minus=-1 in the N=m-1 chart with p_plus=1 in the N=m chart.

The adjacent recurrence gives `kappa_m=-kappa_(m-1)`. Therefore

```text
S_minus=kappa_(m-1)T_0, S_plus=kappa_mT_0=-S_minus; v_minus=(C_minus/S_minus,D_minus/S_minus), v_plus=(C_plus/S_plus,D_plus/S_plus). Real nonzero rescaling does not change the projective points, and c_minus,c_plus<0 select their common short lift.
```

Projective points are unchanged by either real nonzero normalizer. The
proved bound `c_minus,c_plus<-1/4` aligns the two representatives in the
same half-plane and selects the short projective connector.

## Exact Endpoint Wedge

At the old endpoint `p_minus=-1` and new endpoint `p_plus=1`,

```text
A_minus=-sqrt(2+sqrt(2))/4
A_plus =-sqrt(2+sqrt(2))/4
B_minus=-3sqrt(2-sqrt(2))/16
B_plus =-sqrt(2-sqrt(2))/16

A_minus*B_plus-B_minus*A_plus=-sqrt(2)/32.
```

These radicals come from exact evaluation of the completed-zeta and
Riemann-Siegel `C_0` phases at `p=+-1`. No fitted circle or polygon defines
`pi` or the radicals.

## Finite-Height Budget

The critical-ray theorem supplies, independently at both endpoints,

```text
c_sigma=A_sigma+h e0_sigma, |e0_sigma|<100; d_sigma=h B_sigma+h^2 e1_sigma, |e1_sigma|<250, sigma in {minus,plus}.
```

Direct expansion gives

```text
Delta_cut=c_minus*d_plus-d_minus*c_plus=-sqrt(2)h/32+h^2[A(e1_plus-e1_minus)+e0_minus B_plus-B_minus e0_plus]+h^3[e0_minus e1_plus-e1_minus e0_plus].
```

Using `|A|<1/2`, `|B_minus|,|B_plus|<1`, and the saved defect caps,
the two remainder coefficients are at most `500` and `50000`. Since

```text
h<1/72000000000,
500h+50000h^2<1/600,
sqrt(2)/32>1/24,
```

the strict cutoff sign is

```text
Delta_cut/h<-1/25,
Delta_cut<-h/25<0.
```

## Signed Join

For the affine connector

```text
v_s=(1-s)v_minus+s v_plus; c_s*d_(s)-d_s*c_(s)=Delta_cut,
c_minus<-1/4, c_plus<-1/4, c_s<-1/4.
```

Hence

```text
For v_s=(c_s,d_s), partial_s arg(c_s+i d_s)=Delta_cut/(c_s^2+d_s^2)<0.
```

The apparent missing `rho_xx` or second adjacent x-jet is not needed for
this connector: the within-cell x-current uses `d_(edge,x)`, whereas the
cutoff homotopy current differentiates `(c_s,d_s)` with respect to `s` and
uses only the two certified row values.

## Consequence and Boundary

Together with the full critical-ray within-cell sign, this closes the
pointwise cells and every adjacent cutoff by strict clockwise projective
arcs for the retained first-order edge model, including every `q>=1` point
in the stated critical range.

The omitted higher-order Xi/source remainder has not been transferred to
this projective connector, and the edge block has not been inserted into
the complete cumulative contact/minor scalar. No Xi-level edge sign,
cumulative-minor estimate, Abel-scalar gap, successor winding cap, contact
exclusion, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion follows.
