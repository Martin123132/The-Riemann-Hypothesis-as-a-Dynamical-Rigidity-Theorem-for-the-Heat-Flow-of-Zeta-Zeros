# Gaussian B-face window and summable Fresnel remainder

Date: 2026-08-13

Status: exact asymptotic-remainder lemma plus saved-height interval
certificate; not a proof of the complete B-face estimate

Use the exact lower B trace saddle and scale

```text
xi=sqrt(H_B)(x-x_B^-),
H_B=Phi_B''(x_B^-),
ell_B=2sqrt(H_B)/B.                                   (GW1)
```

At `t=10^10`, interval arithmetic gives

```text
H_B=[84898495600365202.65624792695616691950048595680296737501208996420099071288000101841395161978689913904 +/- 3.77e-81],
ell_B=[113.7639659575591168838447929362502704780399209513182261552805317195249336337848032059547189216373739 +/- 2.55e-96].        (GW2)
```

The crossing center of mode `m` is
`xi_m=ell_B(m-m_B^-)`.  The four nearest relevant centers are

```text
xi_620=[-177.0171417448578063781396492298604493793266025390878566682446753263823761012805979467823443629558730 +/- 1.57e-93],
xi_621=[-63.25317578729868949429485629361017890128668158776963051296414360685744246749579474082762544131849914 +/- 1.57e-93],
xi_622=[50.51079017026042738954993664264009157675323936354859564231638811266749116628900846512709348031887473 +/- 1.57e-93],
xi_623=[164.2747561278195442733947295788903620547931603148668217975969198321924248000738116710818124019562486 +/- 1.57e-93].                                    (GW3)
```

Therefore the closed window `|xi|<=70` contains exactly the B crossings
`621` and `622`.  Throughout it, every other positive mode has

```text
|q_B|>75,                                             (GW4)
```

with the nearest left and right margins certified respectively by
`[85.43999157462733963055533524235702944614969244547663418841046936771545431281029092258699643828991365 +/- 1.26e-93]` and
`[75.19229444540905128800747316578856386911775020640150220728567722317996809083189999128509256389458062 +/- 1.25e-93]`.

For `E(r)=exp(i*pi*r^2/2)`, repeated integration by parts gives

```text
T(r)=iE/(pi*r)+E/(pi^2*r^3)-3iE/(pi^3*r^5)+R_T,
|R_T|<=30/(pi^4*r^7).                                 (GW5)
```

Substitution into the sign-adapted current from Section 11.362, followed by
the exact completed phase, parity, and external `1/x` factor, gives

```text
B/[i*pi(Bx-2m)]
-m/[2*pi^2(c-m)^3]
+3i*m*x/[4*pi^3(c-m)^5] + R_ext,m,  c=Bx/2,          (GW6)

|R_ext,m|<=15*m*x^2/[4*pi^4*|c-m|^7].                (GW7)
```

Unlike the leading rational terms, (GW7) is absolutely summable.  Uniformly
on `|xi|<=70`, summing every positive nonlocal mode gives

```text
sum_(1<=m<=39894, m notin {621,622})|R_ext,m|
 <= [7.470760628709106569060729801440952825966599100021607237679652159711519742590894136914089943366345949e-6 +/- 8.40e-100],

sum_(m>=1, m notin {621,622})|R_ext,m|
 <= [7.470760628709106569060729806369077706582819368595599541019052558707095371060331968954272470302183300e-6 +/- 8.40e-100] < 7.48e-6.   (GW8)
```

The infinite tail in (GW8) is bounded by monotone integral comparison; it is
not numerically truncated.

This closes the nonlocal third-order Fresnel remainder on the natural B
window.  It deliberately does not sum the three pole-bearing terms in
(GW6) separately.  Those must be pole-subtracted with modes 621 and 622 and
combined with the exact step current, negative/outer completion, and the
grouped endpoint difference.

Pi provenance: every `pi` is inherited from the equation-(9) Fresnel phase,
the exact face trace, or integration by parts.  No fitted constant is used.

Proof boundary: exact window roster, nonlocal `q_B` separation, three-term
Fresnel expansion, and its absolutely summed remainder only.  The three
rational lattice sums, exact two-mode current, outside-window tail, and
A-fold splice remain open.  This result does not establish the
complete paired residual, `T_upper`, `Lambda<=0`, PF-infinity, RH, or a
prize-level conclusion.
