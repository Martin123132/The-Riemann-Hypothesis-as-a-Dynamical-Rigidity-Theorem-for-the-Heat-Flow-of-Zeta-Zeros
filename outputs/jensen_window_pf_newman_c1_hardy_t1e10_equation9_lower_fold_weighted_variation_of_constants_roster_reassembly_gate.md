# Weighted variation-of-constants roster reassembly

Date: 2026-08-13

Status: exact branchwise finite-Fubini reassembly and conditional packet
coefficients proved; completed source and initial-data cancellation open

Choose event zero as reference and write

```text
r_j=r_0+j Delta_r,
F_sigma(s)=exp(i chi(d(s)))[R_0,sigma+R_Y,sigma+Q_sigma]. (VR1)
```

For either selected Hankel branch, Section 11.388 gives

```text
H_sigma(r_j)=-partial_s K(r_j,r_0)H_sigma(r_0)
 +K(r_j,r_0)H_sigma,r(r_0)
 +Integral_(r_0)^(r_j)K(r_j,s)F_sigma(s)ds.          (VR2)
```

Keep the certified leading-defect weights `a_j` attached.  Finite Fubini
reassembly on the positive branch is exactly

```text
sum_(j=1)^200 a_j H_+(r_j)
 =-D_+ H_+(r_0)+A_+ H_+,r(r_0)
 +Integral_(r_0)^(r_200)
    [sum_(j:r_j>=s) a_jK(r_j,s)]F_+(s)ds,            (VR3)
```

whereas orientation on the negative branch gives

```text
sum_(j=-198)^-1 a_j H_-(r_j)
 =-D_- H_-(r_0)+A_- H_-,r(r_0)
 -Integral_(r_-198)^(r_0)
    [sum_(j:r_j<=s) a_jK(r_j,s)]F_-(s)ds.            (VR4)
```

Here `D_sigma=sum a_j partial_sK(r_j,r_0)` and
`A_sigma=sum a_jK(r_j,r_0)`.  At fixed `s`, (VR3) uses a suffix of
`1..200`; (VR4) uses a prefix of `-198..-1`.  The exact source-cell incidence
counts are `20100` and
`19701`, respectively.  An independent
rational finite-array checker verifies both orientations.

The inverse gauge remains common:

```text
sum a_jG_sigma(d_j)=exp(i beta^3)sum a_jH_sigma(r_j). (VR5)
```

Thus no modewise phase is introduced after (VR3)--(VR4).  Sections
11.390--11.391 supply the branchwise coefficients

```text
|A_-| and negative source packet
 <[3.773903448169599551191850315432134375919658704338404930639827597634795645143518324339296821328192359e-7 +/- 1.40e-107],
|A_+| and positive source packet
 <[3.803311209208802542424917587689813444070251655385279973081910541362153303059285741228013472911439103e-7 +/- 3.15e-107],

|D_-|<[0.0008130582038663737342069787758298923635702310623270155310617545506295678791616509437420557386391222821 +/- 4.34e-104],
|D_+|<[0.0008193938777114114838080827418321118813518608278440911741065844719169044886212161734318185531535825621 +/- 3.74e-104].       (VR6)
```

Equations (VR3)--(VR6) expose the remaining wall precisely.  The `K`
source coefficients are at the `10^-7` scale, but the derivative packets
multiply the branch initial data at the `8.2e-4` scale.  Neither branch
initial value may be discarded.  The event-zero term, the beta-minus-four
identity

```text
H4+sum G4+C4=g4_lambda(0),                            (VR7)
```

and the exact-minus-beta-four source must be rewritten in the same
`H_sigma(r_0),H_sigma,r(r_0),F_sigma` coordinates before any final norm.

This is an exact reassembly for the leading stationary-defect weights only.
It does not identify the finite-integral amplitude remainder with `a_j`, nor
does it prove the cancellation demanded after (VR7).

Pi provenance: (VR5) uses `beta^3=pi C^2/8` and integer mode parity; all
kernel normalizations are inherited from the Airy Wronskian.

Machine-audited companion:

```text
outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_weighted_variation_of_constants_roster_reassembly_gate.md
work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_weighted_variation_of_constants_roster_reassembly_gate.json
work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_weighted_variation_of_constants_roster_reassembly_gate.py
work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_weighted_variation_of_constants_roster_reassembly_gate.py
```

No completed initial-data/source cancellation, exact finite-integral
amplitude remainder, grouped 398-mode splice, complete `Q_K-T` or
`T_upper`, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion is
proved.
