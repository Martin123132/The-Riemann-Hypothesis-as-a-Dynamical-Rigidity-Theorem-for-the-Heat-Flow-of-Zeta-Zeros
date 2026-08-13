# Fixed-B selector-boundary nonzero completed-strip height cell

Date: 2026-08-12

Status: rigorous event-free nonzero height cell; this is not a proof of the
ordinary-Morse join, `T_upper`, `Lambda<=0`, or RH

For fixed `A=159577`, the fold scale is independent of height:

```text
beta^3=pi*A^2/8=t*,        lambda=(t*-t)/beta.         (HC1)
```

The exact Fourier period `hY=2pi` is therefore also fixed.  The all-mode
midpoint plus endpoint half-current continues to collapse the exact and
canonical completed strips to `y=0` for every real `t` in

```text
t*-pi/16 <= t <= t*+pi/16,
|lambda| <= [9.11373501846528081509557975247001911858917670407189115712808383608025308676795678667665259e-5 +/- 3.95e-94].          (HC2)
```

Differentiating the already completed exact-minus-canonical fold integral at
the center gives

```text
delta'(0)=[-1.04128031117923623624828467845055206223170838695394799573070766403050938992304846923848771e-7 +/- 1.70e-28].              (HC3)
```

Its formal first term is `-(13/60)2pi Ai(0)/beta^2`; the rigorous-to-formal
ratio is `[1.00000000000000238134256845544948335352878969830187376108265873142499913228675723075866699 +/- 1.63e-21]`.

A 960-panel interval majorant retains the exact-minus-
canonical difference before taking moduli and proves

```text
sup_|lambda|<=L |delta''(lambda)|
 <= [1.19408956840079163281886326451694192660319643158126144480545378979677138818265906313031738e-6 +/- 2.42e-94].             (HC4)
```

Taylor's theorem applied to the completed object now gives

```text
|delta(lambda)|
 <= [9.49612500096965344123853751189785591899931291950937710064324691677659606545530133694636546e-12 +/- 2.94e-29],
equation-(9) normalized completed-strip error
 <= [8.80375772537209782260368813063480746412771804002560215479153369790086744800590867653722710e-13 +/- 2.72e-30] < 8.9e-13. (HC5)
```

The cell is geometrically clean.  The nearest old-selector internal event is
`pi/8` below `t*`, leaving margin
`[0.196349540849362077403915211454968930262323087460944113810934037019238525392888062414252177 +/- 4.63e-91]`.  The first adjacent-selector event is
mode `40094` at `3103pi/8` above `t*`.
The fixed source endpoint remains `B=5122423` throughout.

This is a genuine nonzero selector-boundary comparison, not a sampled endpoint
check.  It still does not supply the Airy-to-ordinary-Morse remainder match at
the external event buffers.

Machine-audited companion:

```text
outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_fixed_B_selector_boundary_nonzero_height_completed_strip_cell_gate.md
work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_fixed_B_selector_boundary_nonzero_height_completed_strip_cell_gate.json
work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_t1e10_equation9_fixed_B_selector_boundary_nonzero_height_completed_strip_cell_gate.py
work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_fixed_B_selector_boundary_nonzero_height_completed_strip_cell_gate.py
```

No complete `T_upper`, `Lambda<=0`, RH, or prize-level conclusion is proved.
