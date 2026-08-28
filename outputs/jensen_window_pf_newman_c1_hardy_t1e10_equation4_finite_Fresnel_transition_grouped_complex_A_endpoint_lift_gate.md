# Grouped complex A endpoint lift on the transition block

Date: 2026-08-27

Status: rigorous grouped complex endpoint lift; transition Fresnel/Gamma block
not yet evaluated.

The local affine and local exact-minus-affine caches retain the raw complex
endpoint contribution for every mode `39853..39936`.  The exact-exterior gate
retains their grouped complex rational contour value, while the endpoint-tail
gate supplies a rigorous complex modulus allowance for replacing the exact
tail.  Therefore, before physical projection,

```text
R_A^raw=R_A,loc^aff+R_A,loc^nonlinear+R_A,ext^exact,
X_A,endpoint=c_t exp(-i*pi/8)R_A^raw.
```

The certified complex preprojection is

```text
X_A,endpoint = [-0.06463340946584340200520588922629059434888262632560741854752705790299494 +/- 1.51e-12]
             + i [-0.008301313013576003740619392437225776118147834166152429962928301166425154 +/- 1.51e-12].
```

It reconstructs the previously certified physical endpoint through

```text
2 Re X_A,endpoint = [-0.1292668189316868040104117784525811886977652526512148370950541158059899 +/- 3.02e-12].
```

Using the exact `H(t)` and Riemann-Siegel phase gives the natural grouped lift

```text
mathcal A_A^nat=H(t)exp(-i theta(t))X_A,endpoint
               =sum_(m=39853)^39936 m^(-s)C_A,m,

Hardy_t[mathcal A_A^nat]=H(t)A_endpoint.
```

Thus the entire transition contribution can be bounded as one block:

```text
sum_(m=39853)^39936 m^(-s)(beta_m-C_G)
 -mathcal A_A^nat,
```

which is exactly the modewise sum of `m^(-s)(beta_m-C_G-C_A,m)`.  Separate
full complex exterior integration for 84 modes is not logically required.

Pi provenance: `pi` comes from the exact endpoint phase, bi-Morse
normalization, paper rotation, and Riemann-Siegel phase.  No fitted constant
is introduced.

Proof boundary: rigorous grouped complex A endpoint lift at `t=10^10` only.
This gate does not evaluate or bound the joined Fresnel/Gamma transition
packet, prove the `P_W` Mordell join, enclose `J_Z` or `D_K`, or prove a non-A,
all-height, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion.
