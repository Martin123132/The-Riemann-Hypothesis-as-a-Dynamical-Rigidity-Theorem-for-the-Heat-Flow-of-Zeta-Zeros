# Complete lower-fold detuning-lattice quadrature

Date: 2026-08-11

Status: all 84 finite Airy-Fresnel detuning values and their grouped sum
validated; not a proof of the height-uniform fold theorem

For every mode `m=39853..39936`, the cache contains an independent complex
ball for

```text
G_64(d_m)=integral_0^64 Ai(-lambda-y)
 exp(i*y^2/(4beta)-i*d_m*y)dy.                         (DL1)
```

All 84 rows are deterministic, append-only, and individually resumable.
Their exact grouped canonical value is

```text
2pi sum_m G_64(d_m)
 =[-43.33293318556235591049231631821292024307103273683162131356787883669535 +/- 2.26e-22]
  +i*[1.428190029203347133044642455140469773370106139214021263248523619533087 +/- 2.26e-22].              (DL2)
```

Both components overlap the independently computed Airy-first value

```text
[-43.33293318556235591049231631678778100000000000000000000000000000000000 +/- 1.06e-23]
 +i*[1.428190029203347133044642455139145500000000000000000000000000000000000 +/- 1.06e-23]. (DL3)
```

The midpoint differences are

```text
real: [-1.425139243071032736831621313567878836695354949965609320997597695526888e-27 +/- 3.43e-97],
imag: [1.324273370106139214021263248523619533087310347588099861778655215605112e-30 +/- 2.02e-100].              (DL4)
```

This cross-check uses two different organizations of the same canonical
fold: 84 separate detuning transforms in (DL2), and one Airy integral with
the complete finite Dirichlet kernel in (DL3).  It supplies rigorous lattice
state data for the detuning ODE from Section 11.324.

The discrete adjacent variation is recorded as

```text
sum_m |G_64(d_(m+1))-G_64(d_m)|
 =[3.545089444128317770033231384557728650032219030649061887675619327570521 +/- 9.02e-23].              (DL5)
```

It is diagnostic state information, not a termwise final-error bound.  The
next step is to propagate the ODE over detuning and height cells while
preserving the grouped finite sum, then match the outer-Morse limit.

Proof boundary: complete rigorous canonical detuning lattice at one saved
height and `Y=64` only.  No finite-t outer-chart join, complete `T_upper`,
height-uniform source theorem, `Lambda<=0`, RH, or prize-level conclusion is
proved.
