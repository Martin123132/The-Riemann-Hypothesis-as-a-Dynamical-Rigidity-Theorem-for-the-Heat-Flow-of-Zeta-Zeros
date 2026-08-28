# Post-A extended-notch normal form

Date: 2026-08-23

Status: exact post-A projector compression and transition-free upper spectral
edge certified; joined phase-adapted bound open

Let

```text
T={622,...,39894},
W_A={39853,...,39936},
U={622,...,39936}=T union W_A.                     (EN1)
```

The two halves of `W_A` contain 42 modes each: `39853..39894` is already in
the original target, while `39895..39936` lies outside it.  Therefore the
indicator identity

```text
chi_T+chi_A(1-chi_T)=chi_U                            (EN2)
```

holds coefficientwise.  Subtracting the certified finite A transition from
the common projector kernel before any norm gives, for every inherited common
cutoff and regulator,

```text
Delta_(M,epsilon)-mathcal A_(M,epsilon)
 =H_x+sum_(m=-M)^M w_m I_m
  -sum_(m=622)^39936 w_m P_m
  -sum_(m=39853)^39936 w_m(A_m+A_-m).                (EN3)
```

Thus the positive full-line projector is no longer a six-row mask: it is one
contiguous extended notch.  Equivalently,

```text
Delta-mathcal A
 =H_x+integral_0^L f_x(u)[D_(M,epsilon)-G_(U,epsilon)]du
  +sum_(m in U)w_m(A_m+B_m)
  -sum_(m in W_A)w_m(A_m+A_-m).                      (EN4)
```

The endpoint rows in (EN4) reduce exactly to

```text
622..39852:    A_m+B_m,
39853..39936:  B_m-A_-m.                             (EN5)
```

The one-cell fold applies directly to the source integral in (EN4).  The
two-jet reassembly gives the equivalent common-kernel form

```text
H_x+E_x
 +2*pi*i sum_(m=1)^M m*w_m[hat Phi_x(m)-hat Phi_x(-m)]
 -sum_(m=622)^39936 w_m P_m
 -sum_(m=39853)^39936 w_m(A_m+A_-m).                 (EN6)
```

The B window and analytic B exterior remain separate subtractions:

```text
mathfrak R_A=(Delta-mathcal A)-chi_W Btr-O.           (EN7)
```

No B atom or negative Gamma bulk has been removed by (EN1)--(EN7).

The extended target has half-integer edges

```text
c=621.5,             d_A=39936.5.                    (EN8)
```

Hence its exact Gaussian interpolant

```text
h_A,epsilon(y)=e^(-pi*epsilon*y^2)
               [1_(y<c)+1_(y>d_A)]                  (EN9)
```

has integer samples equal to the complement of `U`.  Poisson summation, the
two closed boundary currents, and the pole-free centred-cell recombination
from Sections 11.444--11.445 apply with `d` replaced by `d_A`: both edges are
half-integers and therefore have the same integer phase `(-1)^k`.

This edge translation has a concrete geometric effect.  The A half-boundary
is `A/4=39894.25`.  The old and new upper-edge gaps are

```text
39894.5-A/4=0.25,
39936.5-A/4=42.25.                                   (EN10)
```

For the exact A normal coordinate at `x=1/2`,

```text
q_A(m,1/2)=(A-4m)/2.                                 (EN11)
```

The deleted integer window runs from `q_A=82.5` to `q_A=-83.5`.  The nearest
surviving integers have `q_A=84.5` and `q_A=-85.5`, while the new upper
half-integer edge has `q_A=-84.5` instead of the old `-0.5`.  The certified A
transition has therefore moved the upper spectral edge out of the A
half-boundary collar rather than merely subtracting a scalar after the fact.

Machine-audited companion:

```text
outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_normal_form_gate.md
work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_normal_form_gate.json
work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_normal_form_gate.py
work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_normal_form_gate.py
```

Pi provenance: every `pi` in (EN1)--(EN11) is inherited from the Gaussian
Abel regulator, integer Fourier character, and Kummer phase.  No fitted or
geometric occurrence is introduced.

Proof boundary: exact finite-regulator post-A coefficient compression,
extended-notch fold/modular inheritance, and rational upper-edge geometry
only.  No quantitative bound for the remaining joined kernel, no complete
`R_after_A`, `R_Dir`, or `Q_K-T` estimate, no all-height theorem,
`Lambda<=0`, PF-infinity, RH, or prize-level conclusion is proved.
