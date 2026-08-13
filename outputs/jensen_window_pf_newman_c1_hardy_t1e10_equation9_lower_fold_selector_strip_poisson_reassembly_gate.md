# Exact selector-strip finite-Poisson reassembly

Date: 2026-08-12

Status: exact fixed-upper-endpoint selector reassembly validated; not a proof
of the 399-mode approximation bound or complete source selector splice

Write the old odd source roster as

```text
alpha=C+2u,       u=0,...,L,
f_C(u)=(C+2u)exp(i*pi*x*(C+2u)^2/4).                  (SS1)
```

When the lower selector advances from `C` to `C+2`, the new roster is the
shifted old roster `u=1,...,L`.  Therefore exactly one source lattice point,
`alpha=C`, is removed.  The counts here are

```text
C=159577,       C+2=159579,
old terms=2481423,
new terms=2481422.             (SS2)
```

For integer Fourier mode `m`, shifting `u=v+1` leaves the character unchanged.
The mode difference is the one-cell strip

```text
S_m(x)=integral_0^1 f_C(u)exp(-2pi*i*m*u)du
 =(-1)^m/2 integral_C^(C+2)
   alpha exp(i*pi[x*alpha^2/4-m*alpha])dalpha.         (SS3)
```

The endpoint half-current changes by

```text
H_C(x)=[f_C(0)-f_C(1)]/2.                              (SS4)
```

Finite Poisson on the one-cell interval gives

```text
lim_(M->infinity) sum_(m=-M)^M S_m(x)
 =[f_C(0)+f_C(1)]/2,                                   (SS5)

H_C+lim_sym sum_m S_m=f_C(0)
 =C exp(i*pi*x*C^2/4).                                 (SS6)
```

After multiplication by the Kummer weight and x integration, (SS6) is exactly
the single source term removed when the lower odd selector advances.  Thus
the complete fixed-`B` finite-Poisson representation reassembles exactly.

This also shows why the 399 positive turning modes are insufficient by
themselves.  The selector identity requires the complete symmetric mode sum,
including zero and negative modes, together with the endpoint half-current.
Separating the common odd-endpoint current produces a false divergent object.

The adjacent selector begins with the 399-mode positive turning roster
`39696..40094`,
but its quantitative Airy/Morse estimate must be embedded in (SS6), not used
as a replacement for it.

Proof boundary: exact selector reassembly for a fixed upper endpoint `B` and
the complete symmetric Poisson sum only.  No 399-mode quantitative remainder,
simultaneous upper-endpoint transition, complete source splice, `T_upper`,
`Lambda<=0`, RH, or prize-level conclusion is proved.
