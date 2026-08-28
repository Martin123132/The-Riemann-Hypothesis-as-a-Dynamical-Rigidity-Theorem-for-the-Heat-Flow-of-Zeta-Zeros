# Tiny target-correction K on the first height subcell

Date: 2026-08-28

Status: rigorous exact target-correction Hardy-operator interval; complete
`K_T` assembly remains open

The A-free packet contains the finite correction

```text
P_corr=(C_G-H)D_T,
D_T=sum_(m=622)^39894 m^(-1/2-it),
C_G=(1+exp(-2*pi*t))^(-1/2).                       (TC1)
```

For `h=H'/H` and `L_t=d_t+i theta'-h`, the unsimplified product rule is

```text
L_t[P_corr]=(C_G'-H')D_T+(C_G-H)D_T'
             +(i theta'-h)(C_G-H)D_T.              (TC2)
```

The two apparent `H'D_T` terms cancel exactly, leaving

```text
L_t[P_corr]
 =(C_G-H)(D_T'+i theta'D_T)+(C_G'-h C_G)D_T

 =sum_(m=622)^39894
  [C_G'-h C_G+i(C_G-H)(theta'-log m)]m^(-1/2-it). (TC3)
```

The production evaluator sums the final line of (TC3) mode by mode on
`I_1=[10^10-10^-4,10^10+10^-4]`, so it never forms the separately widened
`D_T'` and `i theta'D_T` pieces.  Arb obtains

```text
K_corr: real [5.500468498557898249657668840437704115237943948644485969454565405637985e-22 +/- 2.21e-21]
        imag [-1.000516952322587066568675222693703882081594161118321505652697656047088e-21 +/- 2.20e-21],

Hardy_t[K_corr]/H=[-3.950973934851031966947260177334543158889079705886615556664764881134033e-22 +/- 6.20e-21].       (TC4)
```

The compact weighted sum overlaps both the simplified two-term assembly and
the raw three-term product rule.  The independent checker raises precision,
reverses all 39,273 modes, switches to the Gamma-duplication formula for
`H'/H`, and repeats all three forms.  A changed six-mode analytic model also
compares (TC3) with a direct numerical derivative of `(C-H)D`.

Pi provenance: `pi` comes only from the exact thermal Gamma coefficient,
Riemann--Siegel phase, and inherited Fourier/Gamma normalization.  No fitted
geometric constant supplies `pi`.

Proof boundary: (TC1)--(TC4) certify only the exact tiny target-correction
operator on `I_1`.  Final assembly with the transition--upper-arc,
lower-plus-ordinary, and positive-real-tail packets remains open.  No wider
`Q_K-T` sign interval, full event-cell theorem, wall handoff, all-height
theorem, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion is proved.
