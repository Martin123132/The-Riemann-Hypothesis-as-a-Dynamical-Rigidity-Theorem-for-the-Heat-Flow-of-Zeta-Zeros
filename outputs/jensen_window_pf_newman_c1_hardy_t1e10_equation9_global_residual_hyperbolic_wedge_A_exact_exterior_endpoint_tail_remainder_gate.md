# Four-term exact-endpoint tail remainder on the A exterior

Date: 2026-08-14

Status: rigorous replacement-error bound for the full exact exterior;
not an exterior value or complete A theorem

For the exact phase-stripped endpoint tail of Section 11.437, write

```text
B_m=b_0+b_1+b_2+b_3+R_4,
b_0=h_A/(iF_m'),   b_(j+1)=-b_j'/(iF_m').            (ET1)
```

Repeated exact integration by parts gives

```text
|R_4(x)| <= integral_0^x |b_3'(u)|du.                (ET2)
```

The derivative is the explicit rational function

```text
b_3'(x)=240 i x^(7/2) P_9(x,m)
 /[pi^4(Ax-2m)^8(Ax+2m)^8],                         (ET3)
```

where the real and imaginary coefficients of `P_9` are bounded term by term.
The interval calculation preserves the `x^(7/2)` vanishing on
`[0,1e-12]`, then uses 4096 deterministic real
slabs for every mode through its exact cutoff.

After integrating the cumulative remainder against the complete exterior
weight `x^(-7/4)(1-x)^(-1/4)`, summing all 84 modes, restoring `1/(2pi)`,
and applying the equation-(9) physical normalization,

```text
max_m sup_x |R_4,m(x)| <=
 [5.5879916143724273713729929564479817431033120576626607257410681549e-11 +/- 4.01e-76],

complete canonical exterior replacement error <=
 [6.4469271173775255601697624215408380018545346108520865288770891074e-13 +/- 7.15e-69],

complete physical exterior replacement error <=
 [2.2823525500438195706919927350738271666944727746985178417424047566e-15 +/- 2.53e-71]. (ET4)
```

An independent checker repeats the entire calculation at 70 digits with
8192 slabs per mode and obtains a tighter nested upper bound.

Pi provenance: every `pi` in (ET1)--(ET4) comes from the exact endpoint
phase, endpoint driver, and equation-(9) normalization.  No fitted constant
is introduced.

Proof boundary: uniform four-term endpoint-tail replacement error over the
full exact A exterior for modes 39853..39936 only.  No value of the resulting
rational common-phase integral, outer-IBP remainder bound, compact nonlinear
amplitude, complete A endpoint theorem, `R_Dir` estimate, RH, or prize-level
conclusion is proved.
