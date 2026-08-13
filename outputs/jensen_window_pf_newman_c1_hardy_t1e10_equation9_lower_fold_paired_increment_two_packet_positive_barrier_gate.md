# Paired selector-increment two-packet positive barrier

Date: 2026-08-13

Status: rigorous positive lower bound for the selector fold increment; this
rejects the signed-increment route and is not a fixed-state residual theorem

For the exact common-profile coefficient

```text
G_ex,(39895+n)=integral_0^1 g_ex(s)e^(-2pi*i*n*s)ds,
```

the paired fold increment has exactly two 200-frequency packets:

```text
D_fold/kappa_C=sum_(n=-200)^-1 G_ex,(39895+n)
              +sum_(n=-79789)^-79590 G_ex,(39895+n),  (PB1)
kappa_C=e^(i beta^3)2sqrt(2).
```

The first packet is evaluated directly from the beta-minus-four profile by
a `65536`-panel Arb midpoint rule.  A `2048`-panel
rectangle atlas bounds the profile and its first two `s` derivatives, giving

```text
[5.245425776827566236943266463129423405943725528079303759875527567736141 +/- 4.02e-70] < Re low beta4
 < [5.642923301672220093734239897989780405654863843125257915178877459062819 +/- 3.37e-70],
low midpoint error <[0.1887212926607737079446085741196316248555691575229770776516749456633393 +/- 2.15e-71].    (PB2)
```

The remote packet is integrated by parts twice.  Integer endpoint phases are
one, so its two boundary currents are explicit and the remaining integral is
bounded by `sup|g4''| sum(2pi n)^-2`:

```text
[-0.02513041023683791850132711520964267513559868834668563379032829443852351 +/- 1.24e-72] < Re high beta4
 < [0.02513288995314875341584634548573750084004698860416141639834642723189638 +/- 3.64e-72],
high second-current remainder <[0.02513164045435434051349233569368362358097488467903680634433736083520995 +/- 2.71e-72]. (PB3)
```

Finally the exact Kummer-ODE theorem gives a uniform profile error.  The two
length-200 Dirichlet packets have the same `L1` majorant, so the complete
exact-minus-beta-minus-four packet error is

```text
canonical packet error <[4.217228401895724188898160646781440387291937587674597294825046991137283e-6 +/- 4.50e-76].      (PB4)
```

Combining (PB1)--(PB4), restoring the carrier, and applying the exact paired
projection proves uniformly on the ordinary top corridor

```text
20<[29.53042617185773784003499392716790466468402423517815945484063879176039 +/- 3.49e-69]<Delta Q_fold
 <[32.06339161020397093705532841772994983869768915680920807149458612209941 +/- 3.17e-69]<35.                  (PB5)
```

Thus the selector increment is one large removed source contribution, not a
small signed correction.  Its positive scale is more than
`[431492.2577309249985754096063820386162038636065218406361821937004900736 +/- 3.50e-65]` times the certified inner
corrected margin.  Telescoping these increments without simultaneous target
and ownership updates cannot close the fixed-state `Q_K-T` residual.

Pi provenance: every frequency and projection factor comes from the exact
Kummer/Fourier character, `beta^3=pi*C^2/8`, and odd-square half-domain
reflection.  No fitted constant is used.

Machine-audited companion:

```text
outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_paired_increment_two_packet_positive_barrier_gate.md
work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_paired_increment_two_packet_positive_barrier_gate.json
work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_paired_increment_two_packet_positive_barrier_gate.py
work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_paired_increment_two_packet_positive_barrier_gate.py
```

No fixed-state fold residual, moving-target telescope, complete `Q_K-T` or
`T_upper`, height-uniform theorem, `Lambda<=0`, PF-infinity, RH, or
prize-level conclusion is proved.
