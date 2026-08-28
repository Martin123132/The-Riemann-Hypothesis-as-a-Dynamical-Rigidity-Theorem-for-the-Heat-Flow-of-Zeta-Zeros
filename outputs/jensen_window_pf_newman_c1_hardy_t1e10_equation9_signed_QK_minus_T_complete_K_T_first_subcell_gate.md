# Complete K_T on the first height subcell

Date: 2026-08-28

Status: rigorous complete `K_T` and strictly positive joined height derivative
on the first nonzero-radius subcell; not a proof of any wider or global theorem

On `I_1=[10^10-10^-4,10^10+10^-4]`, exact split-contour ownership gives

```text
K_T=K_(L+O)+(K_tr+K_U)+K_tail+K_corr.              (KT1)
```

The first two parentheses are cancellation-preserving joins.  The positive-real
tail is inserted as a centered complex rectangle obtained from its rigorous
modulus bound, while the exact finite correction is retained with its phase.
Arb obtains the complete complex enclosure

```text
K_T: real [0.001075832723884759782215980917140674724944806468367439016634794394864449 +/- 2.13e-5]
     imag [-0.0007637304928410179862501573958117249009178211545193283492887941611183215 +/- 2.14e-5].         (KT2)
```

The exact A-free derivative identity is

```text
(Q_K-T)'=Hardy_t[K_T]/H.                            (KT3)
```

There are two production enclosures of its right side.  Direct projection of
(KT2) gives

```text
Hardy_t[K_T]/H=[0.0007213559574665850959718227386474609375000000000000000000000000000000000 +/- 6.33e-5].  (KT4)
```

Exact linearity also permits the four already stabilized projected packets to
be added without reintroducing the large shared phase dependency:

```text
lower plus ordinary : [0.0007211635411294992081820964813232421875000000000000000000000000000000000 +/- 6.22e-5]
transition plus arc : [1.928065307995447362487961839860872714780271053314209000000000000000000e-7 +/- 1.08e-6]
positive tail bound : [8.894313554516310380487661293621294175740210658602407276e-1873216240 +/- 3.43e-1873216295]
tiny correction     : [-3.950973934851031966947260177334543158889079705886615556664764881134033e-22 +/- 6.21e-21]

Hardy_t[K_T]/H=[0.0007213563476602987525232478840221250780767520093718771050110920294113384 +/- 6.32e-5] > 0. (KT5)
```

As a deliberately coarser guard, the transition--arc midpoint is discarded
and only its certified absolute envelope `1.3e-6` is retained.  Even then,

```text
(Q_K-T)' in [0.0007211635411294992077869990878381389908052739822665456841110920294113384 +/- 6.35e-5] > 0. (KT6)
```

The independent checker changes precision and summation order, rebuilds the
phase and positive `H` transport, verifies direct/projected linearity on an
altered four-packet fixture, and repeats the envelope-only positivity test.

Pi provenance: every `pi` is inherited from the exact Mellin--Fresnel kernel,
Fourier spacing, Gamma normalization, Riemann--Siegel phase, or `p=t/(2*pi)`.
No circle, polygon, visual symmetry, or fitted geometric constant supplies
`pi`.

Proof boundary: (KT1)--(KT6) prove the complete nonzero-radius `K_T` enclosure
and `(Q_K-T)'>0` only on `I_1`.  The previously certified `Q_K-T<0` interval is
therefore increasing there, but no larger sign interval follows without a new
adaptive cover.  No full event-cell theorem, event-wall handoff, all-height
transport, equation-(4) global infinite-series identity, `Lambda<=0`,
PF-infinity, RH, or prize-level conclusion is proved.
