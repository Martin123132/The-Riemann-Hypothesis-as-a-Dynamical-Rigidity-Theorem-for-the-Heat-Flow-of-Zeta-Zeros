# Fixed-regulator embedding of the translated A-face channel

Date: 2026-08-24

Status: exact saved-height channel embedding and complement definition;
the non-A complement is not bounded

Let

```text
T={622,...,39894},       U={622,...,39936},
Q=U minus T={39895,...,39936}.                    (EM1)
```

The post-A finite-regulator remainder `mathfrak R_A,(M,epsilon)` is the
common object of Section 11.443.  Sections 11.449--11.453 prove that, for
every `epsilon>0`, translating the upper notch edge from `39894.5` to
`39936.5` produces exactly

```text
C_U,epsilon-C_T,epsilon
 =-sum_(m in Q) exp(-pi*epsilon*m^2)exp(-2*pi*i*m*s),

h_U,epsilon-h_T,epsilon
 =-q_epsilon(y)1_(39894.5<y<39936.5).               (EM2)
```

After the exact A-face lift and removal of its common unit-modulus phase,
the paired edge-translation channel is therefore

```text
E_(42,epsilon)(x)
 =[-sum_(m=39895)^39936 f_(epsilon,x)(m)]
   -[-integral_(39894.5)^39936.5 f_(epsilon,x)(y)dy]. (EM3)
```

Thus the finite block has coefficient `-1` and the translated strip has
coefficient `+1`.  Both pieces in (EM3) are one channel.  Neither may be
estimated or charged as a second copy of the 42-mode row.

For every common cutoff `M>=B=5122421`, let
`mathfrak E_(A,42),epsilon` denote the exact lifted channel (EM3), and define

```text
mathfrak R_nonA,(M,epsilon)
 :=mathfrak R_A,(M,epsilon)-mathfrak E_(A,42),epsilon.

mathfrak R_A,(M,epsilon)
 =mathfrak E_(A,42),epsilon
  +mathfrak R_nonA,(M,epsilon).                     (EM4)
```

The first term is independent of `M`: its finite support ends at `39936`
and its continuous strip is compact.  Equation (EM4) is an analytic-channel
decomposition of the existing post-A remainder.  It does **not** delete a
coefficient row.  In particular, modes `39895..39936` retain the unchanged
six-class summand

```text
P_-m+B_m+B_-m.                                      (EM5)
```

The A notch remains removed, and no `P_m`, `A_m`, or `A_-m` atom is
reintroduced.

The Abel-zero physical split is legitimate.  On the finite strip,
`exp(-pi*epsilon*y^2)` tends uniformly to one.  On `0<x<=x_16`, the exact
six-current formula has a factor `x` in every retained term and its rigorous
remainder has a factor `x^6`; hence `G_A(y,x)=O(x)` uniformly on the strip.
The physical density is consequently

```text
x^(-5/4)(1-x)^(-1/4) O(x)=O(x^(-1/4)),              (EM6)
```

which is integrable at zero.  The interval `x_16<=x<=1/2` is compact.
Dominated convergence, followed by linearity of the already fixed common
limit order, gives

```text
I_(A,42)=lim_(epsilon down 0)
          mathcal P_t[mathfrak E_(A,42),epsilon],

R_nonA=lim_(epsilon down 0)lim_(M to infinity)
        mathcal P_t[mathfrak R_nonA,(M,epsilon)],

R_after_A=I_(A,42)+R_nonA.                          (EM7)
```

Section 11.458 proves `|I_(A,42)|<0.00364`.  Therefore the inherited target
`R_after_A<0.0368147039947` follows from the sufficient, deliberately
strong condition

```text
|R_nonA| < 0.0331747039947,

because R_after_A
 <=|I_(A,42)|+|R_nonA|<0.0368147039947.             (EM8)
```

This gate derives (EM8); it does not prove it.

Machine-audited companion:

```text
outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_fixed_regulator_embedding_gate.md
work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_fixed_regulator_embedding_gate.json
work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_fixed_regulator_embedding_gate.py
work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_fixed_regulator_embedding_gate.py
```

Pi provenance: no new occurrence of `pi` is introduced.  Equation (EM2)
inherits the Gaussian Abel/Fourier normalization, (EM3) inherits the
canonical Fresnel/Kummer phase, and (EM6)--(EM8) introduce no `pi`.

Proof boundary: the translated 42-mode A-face channel is embedded exactly in
the saved-height common-regulator `R_after_A`, its physical Abel-zero split is
justified, and the sufficient non-A target is derived.  The non-A complement,
joined `R_after_A`, `R_Dir`, `Q_K-T`, all-height theorem, `Lambda<=0`,
PF-infinity, RH, and every prize-level conclusion remain unproved.
