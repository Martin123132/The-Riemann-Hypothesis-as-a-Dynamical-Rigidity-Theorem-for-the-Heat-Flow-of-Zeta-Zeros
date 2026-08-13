# Airy Green kernel for the selected-branch detuning equation

Date: 2026-08-13

Status: exact gauge conjugation, variation-of-constants kernel, and uniform
pointwise kernel envelope proved; no grouped Green-operator bound

For the detuning equation (11.387.2), put

```text
chi(d)=2beta^2 d+beta d^2,
r=beta^2+lambda+2beta d,
G(d)=exp[-i chi(d)]H(r).                              (GK1)
```

Direct differentiation gives the exact conjugation

```text
L_beta G=exp[-i chi(d)](H_rr+rH).                     (GK2)
```

Thus if `L_beta G=S`, the transformed equation is simply

```text
H_rr+rH=exp[i chi(d(r))]S(d(r)).                     (GK3)
```

Its causal Green kernel is the Airy determinant

```text
K(r,s)=pi[Ai(-r)Bi(-s)-Bi(-r)Ai(-s)],                (GK4)
K(s,s)=0,       partial_r K(s,s)=1.                  (GK5)
```

Consequently, for any reference point `r_0`,

```text
H(r)=-partial_s K(r,r_0)H(r_0)+K(r,r_0)H_r(r_0)
 +Integral_(r_0)^r K(r,s)e^[i chi(d(s))]S(d(s))ds.   (GK6)
```

For `d_m=beta(4m/C-1)`, the Green coordinate is affine:

```text
r_m=beta^2(8m/C-1)+lambda,
Delta_r=8beta^2/C,
r_(q+j)=r_q+j Delta_r,
r_(q-j)=r_q-j Delta_r.                               (GK7)
```

The event-ordered pairs therefore form an exactly symmetric Airy lattice,
not merely an approximately regular sample.  On the full top corridor,

```text
[4595460.49636281631441140967198588125771989228195278633114051663533798837583542906822287403 +/- 4.49e-83] <= r <= [4688073.02641084492514309318692142001202428332278035181173676342487014138405818579998314640 +/- 4.73e-83],
Delta_r=[232.694798886661458660881172825585318541936808144908828028964517489628567492347238423600029 +/- 2.12e-87].              (GK8)
```

DLMF 9.8.20 and its signed first-neglected-term remainder imply, at order
zero and for positive `r`,

```text
M(-r)^2<=1/(pi sqrt(r)).                              (GK9)
```

Using `M(-r)^2=Ai(-r)^2+Bi(-r)^2` in the determinant gives

```text
|K(r,s)|<=pi M(-r)M(-s)<=(rs)^(-1/4)
          <[0.000466482634806340821896175549566115059023874252868881411276861381617387754052114340490062689 +/- 2.61e-93]<4.67e-4. (GK10)
```

The pointwise bound alone is not summed over the `r` interval: doing so by
absolute values would throw away the equally spaced event-pair phases.  The
next estimate must insert the compressed sources (11.387.7)--(11.387.8) and
the grouped forcing projection into (GK6), then use the Airy phase and finite
Abel differences on the symmetric lattice.

Pi provenance: the Airy Wronskian contributes the `pi` in (GK4); the lattice
scale comes from `beta^3=pi C^2/8` and the exact Kummer/Fourier detuning.
The modulus theorem is sourced from `https://dlmf.nist.gov/9.8`.

Machine-audited companion:

```text
outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_selected_beta4_detuning_airy_green_kernel_gate.md
work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_selected_beta4_detuning_airy_green_kernel_gate.json
work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_selected_beta4_detuning_airy_green_kernel_gate.py
work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_selected_beta4_detuning_airy_green_kernel_gate.py
```

No grouped Green-operator estimate, endpoint-completion cancellation,
finite-integral 398-mode splice, all-corridor continuation, complete
`Q_K-T` or `T_upper`, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion
is proved.
