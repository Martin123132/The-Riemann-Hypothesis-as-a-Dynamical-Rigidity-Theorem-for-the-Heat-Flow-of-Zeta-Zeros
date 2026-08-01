# Newman Dominant-Ray Connector Phase Cap

Date: 2026-07-28

Status: exact right-connector phase cap and reduced horizontal phase
ledger; not a proof of the remaining Xi boundary inequalities,
`Lambda<=0`, PF-infinity, RH, or a Clay-prize result.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_first_order_centered_dominant_ray_connector_phase_cap.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_dominant_ray_connector_phase_cap.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_dominant_ray_connector_phase_cap.py
```

## Pi Provenance

The pi in x=4*pi*exp(L), the completed-zeta normalizer, and the Riemann-Siegel saddle is the same universal Euclidean constant already present in xi. The 2*pi phase turn is the period of exp(i*theta). The classical rational bounds 3<pi<22/7 are used only to turn the source-derived angular bounds into strict rational inequalities; no circle or polygon is selected by the connector argument.

## Ray Geometry

```text
For L=L_j=101+j and x=R_j=4*pi*exp(L), the positively oriented right edge of D_j traverses I_L=[25/L,25/(L-1)]. Its length is Delta t=25/[L*(L-1)], and t*L>=25 throughout.
Put Z_t=H_t/A_t and V(t)=(Z_t(x),partial_x Z_t(x)/L). The positive normalizer and its derivative act on the full first jet by a positive determinant shear, so V may be used on the connector without changing the closed-boundary degree.
```

## One-Saddle Ellipse

```text
At fixed x, M_t(s)=exp[t*alpha(s)^2/4]*M_0(s), so beta_t=arg M_t(s) obeys partial_t beta=Im(alpha^2)/4=Re(alpha)*Im(alpha)/2. For x=4*pi*exp(L), Re(alpha)<L/2 and |Im(alpha)|<pi/4; hence |Delta beta|<25*pi/[16*(L-1)]<=pi/64.
With b=-partial_x beta and lambda=b/L, the certified bounds give 6/25<=lambda<=13/50. The one-saddle normalized jet is V_0(t)=(2*cos(beta),2*lambda*sin(beta)), so |V_0|>=12/25.
```

## Uniform Cone

```text
The finite arithmetic tail and exact-H collar give |V_1-V_(0,1)|<=32/125+1/8000=2049/8000 and |V_2-V_(0,2)|<=6901/100000+7/160000=55243/800000. Their Euclidean norm is strictly below 4/15, so |V-V_0|/|V_0|<5/9 uniformly, including cutoff transitions.
Because sin(3/5)>3/5-(3/5)^3/6=141/250>5/9, the relative angle delta(t)=arg V(t)-arg V_0(t) has one continuous branch with |delta(t)|<3/5 on the whole connector. Equivalently V/V_0 remains in a disk contained in one open half-plane, so no hidden relative turn is possible.
```

## Connector Cap

```text
For theta_0=arg(cos(beta)+i*lambda*sin(beta)), |d theta_0|<=(1/lambda)|d beta|+[1/(2*lambda)]|d lambda|. At fixed x, lambda is affine in t; its certified range has width 1/50. Therefore |Delta theta_0|<25*pi/384+1/24<331/1344, using pi<22/7.
A continuous argument lift of V on I_L has total angular range below 6/5+331/1344=9719/6720<3/2<pi/2. In particular |Delta_(right D_j)arg V|<pi/2 for every j>=0. The exact dominant-ray connector therefore costs strictly less than one quarter turn.
```

The bound is uniform across every cutoff transition. It controls the
continuous argument lift of the complete connector, not just its endpoints.

## Phase Ledger

```text
On the symmetry axis H_t(0)>0 and partial_x H_t(0)=0. Its normalized first-jet path lies on the positive real ray and contributes zero phase.
Insert the exact normalized connector into the boundary proxy through the already required nonvanishing homotopies, assigning the two endpoint tracks once to the complementary path. Write C_j for the connector phase and H_j for every other oriented piece: top and bottom q>=1 Abel arcs, q<1 arcs, finite/core shoulders, adjacent-cutoff and chart joins, and the two connector endpoint tracks. Then 2*pi*kappa_j=H_j+C_j, |C_j|<pi/2, and kappa_j is a nonnegative integer.
The strict signed inequality H_j<3*pi/2 is sufficient for H_j+C_j<2*pi, hence 0<=kappa_j<1 and kappa_j=0. The vertical connector is no longer an open arithmetic phase theorem; the remaining phase budget is the completely joined horizontal composite, including its endpoint tracks.
```

The vertical arithmetic connector has therefore been discharged. The two
connector endpoint tracks remain in the horizontal composite so that no
open-path phase is silently lost.

## Seam Guard

The cap is proved for the exact normalized Xi first jet on the right edge, not for an unjoined Abel chart. Replacing the horizontal exact jet by the Abel proxy requires the contact-band gap and explicit nonvanishing endpoint tracks. Those tracks are part of H_j and may not be discarded or charged twice.

## Endpoint-Only Guard

```text
Endpoint agreement alone cannot cap an open-path phase: G_m(u)=exp(2*pi*i*m*u), 0<=u<=1, has G_m(0)=G_m(1)=1 but Delta arg G_m=2*pi*m. The connector proof instead uses the uniform pointwise cone |arg(V/V_0)|<3/5, which forbids this hidden winding.
```

This is a generic path countermodel, not an Xi counterexample.

## Remaining Target

```text
Prove the Abel-scalar pointwise gap on every q>=1 horizontal chart and prove H_j<3*pi/2 after the top, bottom, q<1, finite shoulder, cutoff/chart, and connector-seam pieces are joined. No absolute-value sum over local carrier fluxes is licensed.
```

## Boundary

The ray geometry, normalizer time-phase identity, one-saddle ellipse, uniform exact-H cone, strict connector phase cap, zero axis phase, reduced horizontal ledger, and endpoint-only winding guard are exact or inherited from certified source bounds. This does not prove the Xi Abel-scalar gap, the 3*pi/2 horizontal-composite inequality, q<1 closure, finite shoulders or endpoint tracks, complete boundary composition, contact exclusion, Lambda<=0, PF-infinity, RH, or a Clay-prize conclusion.
