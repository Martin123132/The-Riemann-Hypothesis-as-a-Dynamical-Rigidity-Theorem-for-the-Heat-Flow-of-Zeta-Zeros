# Contact-Centered Mangoldt/Abel Equivalence Gate

Date: 2026-07-30

Status: exact contact centering, Abel duality, and balanced hyperbola
decomposition. The required signed lower bound remains open. This is not
a proof of contact exclusion, `Lambda<=0`, PF-infinity, RH, or a
Clay-prize result.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_first_order_centered_mangoldt_contact_centering_gate.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_mangoldt_contact_centering_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_mangoldt_contact_centering_gate.py
```

## Pi Provenance

The `pi` in `a^2=x/(4*pi)+t/16` is still the completed-zeta and
Riemann-Siegel normalization. Contact centering, Abel summation, the
Mangoldt identity, and the balanced hyperbola introduce no new `pi`.

## One Physical Mangoldt Moment

```text
Put Z_Lambda=eta[H_1+varrho_1 H_2+varrho_2 H_3]. Then Z_Lambda=sum_(n<=N)log(n)z_n=sum_(dm<=N)Lambda(d)z_(dm)=(1/2)sum_(dm<=N)[Lambda(d)+Lambda(m)]z_(dm).
```

The new coefficients are called `varrho_1,varrho_2`; the established
symbol `rho_1=partial_x log|f_1|` keeps its earlier meaning.

```text
With S=sum_(n<=N)z_n and L_a=log(a), U=L_a S-Z_Lambda=L_a(W_0-e)-Z_Lambda.

Define E_N=g-s_*'L_a e and mathfrak B_N=E_N-s_*'Z_Lambda+s_*'L_a W_0. Then mathfrak B_N=g+s_*'U exactly.
```

Thus the four correction-free moments do not remain four independent
objects at contact. They combine into one physical logarithmic derivative
plus the already present value `W_0`.

## Contact Fibres

```text
The endpoint-complete scalar is mathcal_C_N=Re(mathfrak B_N)-c u_N mathsf_X. Equivalently, for W_0=mathsf_X+i mathsf_Y and u_N=L_a-log(N), mathcal_C_N=Re[E_N-s_*'Z_Lambda]+c log(N)mathsf_X-b L_a mathsf_Y.

At W_0=0, mathfrak B_N=E_N-s_*'Z_Lambda and mathcal_C_N=Re[E_N-s_*'Z_Lambda]. No division by W_0, H_a, or a carrier is used.

At mathsf_X=0, mathcal_C_N=Re[E_N-s_*'Z_Lambda]-b L_a mathsf_Y. The unrestricted imaginary fibre remains explicit.
```

The zero fibre is genuinely simpler, but the real-crossing fibre still
contains `mathsf_Y`; no sector or sign has been assumed.

## Abel Duality

```text
For F_k=sum_(n<=k)z_n and h_k=log((k+1)/k), Z_Lambda=log(N)S-sum_(k=1)^(N-1)h_kF_k. Hence at W_0=0, E_N-s_*'Z_Lambda=V_N+s_*'sum_k h_kF_k, where V_N=g-s_*'u_Ne.
```

The independent rational audit uses
`5` carriers and has zero mismatches.
This identity proves that the prefix and Mangoldt coordinates are two
representations of the same jet, not two sources of evidence.

## Endpoint

```text
With kappa=(-1)^N B_0 real, e=kappa(T_0+i)H_a and g=kappa(T_0+i)J_a, E_N=kappa(T_0+i)[J_a-s_*'L_a H_a]. This identity remains division-free when H_a=0.
```

The formula does not divide by `H_a`, so the endpoint-zero fibre remains
inside the theorem statement.

## Balanced Hyperbola

```text
Let D=floor(sqrt(N)). Since no pair dm<=N has d>D and m>D, Z_Lambda equals (1/2)sum_(d,m<=D)[Lambda(d)+Lambda(m)]z_(dm)+sum_(d<=D<m<=N/d)[Lambda(d)+Lambda(m)]z_(dm).
```

The checker reconstructs this split at
`8` cutoffs with exact rational-complex
carriers and zero coefficient mismatches. The small square has no hidden
cutoff because `floor(sqrt(N))^2<=N`; every remaining pair lies in the
written wing. This is an exact organization, not a cancellation bound.

## Live Handoff

```text
On |mathsf_X|<=delta_L, the existing shear bound |c u_N mathsf_X|<epsilon_term shows that |Re(mathfrak B_N)|>A_L+2epsilon_term is sufficient for the Abel gap. This estimate is not proved here.

Treat Z_Lambda, the Abel-prefix sum, and the centered jet as three exact coordinates of one object. Apply a Type-I/II or Vaughan decomposition to Z_Lambda only after composing E_N and s_*'L_aW_0, and retain the balanced-wing boundary. Do not count the coordinates as independent evidence.
```

The next candidate must estimate the combined endpoint, Mangoldt moment,
and `W_0` term before absolute values. It must retain the wing boundary,
`W_0=0`, `mathsf_X=0`, `H_a=0`, `q=1`, cutoff equality, and adjacent
charts.

## Boundary

This gate proves the contact-centering identity, one physical Mangoldt
moment, exact Abel-prefix duality, recurrent endpoint factorization,
zero-fibre and real-crossing formulas, notation separation, and balanced
hyperbola decomposition. It proves no signed lower bound, Abel-scalar gap,
horizontal successor winding cap, contact exclusion, Q209, cofinal
descendant theorem, `Lambda<=0`, PF-infinity, RH, or prize-level result.
