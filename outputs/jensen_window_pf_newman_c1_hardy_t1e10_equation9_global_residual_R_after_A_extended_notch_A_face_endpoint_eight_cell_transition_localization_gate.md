# Eight-cell localization of the A-face endpoint transition

Date: 2026-08-23

Status: exact certificate localizing the special-function endpoint layer to
eight midpoint cells; the eight-cell interval integral remains open

Split the 42-cell strip exactly as

```text
[39894.5,39936.5]
 =[39894.5,39902.5] union [39902.5,39936.5].          (EL1)
```

The first interval contains precisely the midpoint modes `39895,...,39902`;
the second contains precisely `39903,...,39936`.  Therefore

```text
D_42(x)=D_(tr,8)(x)+D_(tail,34)(x),
E_42(x)=-D_42(x).                                    (EL2)
```

On the complete endpoint layer `x_16<=x<=1/2`, the smallest normal-tail
coordinate in the 34-cell block occurs at `(y,x)=(39902.5,1/2)` and is

```text
-q_A=2*(39902.5-A/4)=33/2=16.5.                     (EL3)
```

Hence the six-current theorem applies uniformly to all 34 tail cells.  Its
pointwise error increases with `x` and decreases with `y`, giving

```text
sup_tail |G_A-G_6|
 <=[8.1760625933336825095042782997168188560880707726369730416112722211163328484654990e-11 +/- 9.60e-92],

|D_(tail,34)-D_(tail,34),6|
 <=[5.5597225634669041064629092438074368221398881253931416682956651103591063369565393e-9 +/- 1.35e-89].     (EL4)
```

The factor in the second line is `68=34+34`: 34 point masses plus strip
length 34.  The reduced tail defect is again the explicit rational-log
formula of Section 11.454, now restricted to the 34-mode roster.

After the exact equation-(9) physical normalization, replacing the endpoint
tail by its six-current form costs at most

```text
[5.4079206670534129761782425685682766729826526391215214490413340655583402084037140e-15 +/- 3.04e-95]. (EL5)
```

At the layer boundaries the rational-log tail itself has certified moduli

```text
|D_(tail,34),6(x_16)|
 =[1.8536422985122183777039986720042979954039468151006936879745742793146201538092856 +/- 1.81e-80],

|D_(tail,34),6(1/2)|
 =[7.4576156232238033141740191713568373170657409005217230142649409762505771640026088 +/- 3.23e-80]. (EL6)
```

Only

```text
D_(tr,8)(x)
 =sum_(m=39895)^39902 G_A(m,x)
  -integral_(39894.5)^39902.5 G_A(y,x)dy             (EL7)
```

still requires exact Fresnel enclosure on the endpoint layer.  This is a
compact rectangle in the normal variables: eight unit cells and
`1/4<=39894.5-A*x/2<=7.9993`.  No remote mode or 34-cell tail must enter that
special-function calculation.

Floating route telemetry, not used as proof, gives the complete endpoint
layer candidate `-0.000079920799443454167892760164806031637091850371775049` with
coarse/fine relative-integral difference
`2.6562e-12`.

Machine-audited companion:

```text
outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_endpoint_eight_cell_transition_localization_gate.md
work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_endpoint_eight_cell_transition_localization_gate.json
work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_endpoint_eight_cell_transition_localization_gate.py
work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_endpoint_eight_cell_transition_localization_gate.py
```

Pi provenance: every `pi` in (EL1)--(EL7) is inherited from the canonical
Fresnel primitive, Gaussian Abel factor, equation-(9) normalization, and
exact Kummer/Fourier phase.  The integer eight-cell split is forced by the
next half-integer edge with endpoint `|q_A|>16`, not fitted to telemetry.

Proof boundary: exact eight-cell/34-cell endpoint decomposition, uniform
six-current tail replacement below `5.56e-9`, and physical replacement cost
below `6e-15` only.  No interval value for the eight-cell transition,
complete endpoint-layer value, full signed A-face bound, quantitative
`R_after_A`, `R_Dir`, or `Q_K-T` estimate, all-height theorem, `Lambda<=0`,
PF-infinity, RH, or prize-level conclusion is proved.
