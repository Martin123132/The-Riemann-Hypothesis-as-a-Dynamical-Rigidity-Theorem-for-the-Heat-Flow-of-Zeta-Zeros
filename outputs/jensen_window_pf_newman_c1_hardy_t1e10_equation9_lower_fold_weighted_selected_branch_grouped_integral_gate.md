# Grouped selected-branch leading-defect integral

Date: 2026-08-13

Status: rigorous direct grouped midpoint/height enclosure; valid but too
coarse for the local normalized headroom

The selected leading-defect correction is kept as one finite integral,

```text
W_sel(mu)=Integral_0^Y e^(iy^2/(4beta))
 sum_(m!=q) a_m(mu)e^(-id_my)C_sel,m(lambda,y)dy.    (SG1)
```

At each height the two mode halves are evaluated as whole 198- and 200-term
Fourier polynomials.  No modewise integral or absolute-value sum is taken.
Composite midpoint quadrature on `65536` panels is enclosed using a
uniform second-derivative bound; five height nodes are transported over their
adjacent subintervals by an analytic `lambda`-derivative majorant.

The resulting ordinary-top-corridor certificate is

```text
max node modulus <[3.165275519133648783211842985660279667739141524876319869628766999713402e-5 +/- 1.80e-75],
y-quadrature error <[4.971215543876200470358556210857949815761673655338850843099502637910145e-6 +/- 5.44e-75],
height transport error <[4.264143930472166868561576187907016285319395345910392244061149095573490e-6 +/- 5.14e-75],
sup_mu |W_sel(mu)|<[4.088811466568485517103856225536776277847248425001244178344832173061765e-5 +/- 1.35e-74].       (SG2)
```

The direct grouped values remain near `3.2e-5`, but the present analytic
derivative envelope is deliberately conservative.  The bound is therefore
a certified interface, not yet a closure of the previous local headroom.
Its next sharpening should exploit the exact Airy equation inside the
grouped kernel rather than a global `sum|a_m|` derivative majorant.

Pi provenance: the height radius is `pi/(16t*)`, the fold scale obeys
`beta^3=pi C^2/8`, and the mode phases use `hY=2pi`.  No unrelated circle
normalization enters.

Machine-audited companion:

```text
outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_weighted_selected_branch_grouped_integral_gate.md
work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_weighted_selected_branch_grouped_integral_gate.json
work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_weighted_selected_branch_grouped_integral_gate.py
work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_weighted_selected_branch_grouped_integral_gate.py
```

No exact finite-integral amplitude remainder, local-headroom closure,
completed source/initial-data splice, complete `Q_K-T` or `T_upper`,
all-corridor theorem, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion
is proved.
