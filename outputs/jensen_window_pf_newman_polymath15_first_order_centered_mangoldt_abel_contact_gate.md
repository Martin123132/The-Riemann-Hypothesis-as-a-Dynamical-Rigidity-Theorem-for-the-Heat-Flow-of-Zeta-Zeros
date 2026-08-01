# Corrected Mangoldt-Abel Contact Gate

Date: 2026-07-30

Status: exact first-order moment and symmetric divisor-pair reduction
with a nonpromotion guard. This is not a proof of the Xi Abel gap,
contact exclusion, `Lambda<=0`, PF-infinity, RH, or a Clay-prize result.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_first_order_centered_mangoldt_abel_contact_gate.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_mangoldt_abel_contact_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_mangoldt_abel_contact_gate.py
```

## Pi Provenance

The `pi` in `a^2=x/(4*pi)+t/16` remains the completed-zeta and
Riemann-Siegel normalization. No new `pi` is introduced by the
Mangoldt convolution, a circle, a polygon, or the legacy images.

## Exact First-Order Collapse

```text
Put B=alpha'(s)t^2/8, a_0=1+1/(6s)+alpha'(s)t/4+B alpha(s)^2, a_1=-2B alpha(s), and a_2=B. Then, exactly, 1+d_n=a_0+a_1 log(n)+a_2 log(n)^2.

Since a_0=1+d_1 and |d_1|<1/2, put varrho_1=a_1/a_0 and varrho_2=a_2/a_0. With q_n^(0)=exp[t log(n)^2/4-s_*log(n)] and eta=f_1/|f_1|, the anchored carrier is z_n=eta q_n^(0)[1+varrho_1 log(n)+varrho_2 log(n)^2].
```

The retained correction is therefore not an uncontrolled perturbation:
it changes the correction-free Dirichlet chain by exactly two logarithmic
moments.

```text
For H_j=sum_(n<=N)log(n)^j q_n^(0), j=0,1,2,3, S=eta[H_0+varrho_1 H_1+varrho_2 H_2] and U=eta[log(a)H_0+(varrho_1 log(a)-1)H_1+(varrho_2 log(a)-varrho_1)H_2-varrho_2 H_3].
```

This is an exact collapse to four finite logarithmic moments.
The recurrent endpoint is still present through `e` and `g`; it has not
been absorbed into a bulk estimate.

## Symmetric Mangoldt Form

Using `sum_(d|n)Lambda(d)=log(n)` and swapping `d,m` gives

```text
For j=1,2,3, the finite moment is exactly H_j=sum_(dm<=N)Lambda(d)log(dm)^(j-1)q_(dm)^(0)=(1/2)sum_(dm<=N)[Lambda(d)+Lambda(m)]log(dm)^(j-1)q_(dm)^(0).

For dm<=N, q_(dm)^(0)=q_d^(0)q_m^(0)exp[(t/2)log(d)log(m)].
```

Thus the actual symmetry is transpose symmetry on the multiplicative
hyperbola `dm<=N`. The divisor audit checked
`128` exact coefficient rows
through `n=128` with zero mismatches.

The untruncated heat coupling has the exact square expansion

```text
For t>=0 the untruncated heat coupling G_t(d,m)=exp[(t/2)log(d)log(m)] is positive semidefinite, because sum_(d,m)y_d y_m G_t(d,m)=sum_(r>=0)(t/2)^r/r! [sum_d y_d log(d)^r]^2.
```

but the retained moment is a complex bilinear form:

```text
Combining the heat factorization with the preceding identity writes H_j=sum_(dm<=N)K_(j,t,N)(d,m)q_d^(0)q_m^(0), where K is real and transpose-symmetric. This is a complex bilinear form, not a Hermitian form.
```

## Contact Normal Form

```text
With e=eta r_0=g_0/|f_1|, g=eta r_A=G_a/|f_1|, W_0=e+S, mathsf_X=Re(W_0), s_*'=c+ib, and u_N=log(a/N), the exact endpoint-complete scalar is mathcal_C_N=Re[g+s_*'U]-c u_N mathsf_X, with S and U given by the four-moment formulas.
```

On the hard bottom intervals the already certified contact-band shear
therefore gives the sufficient handoff

```text
On |mathsf_X|<=delta_L, the already certified terminal shear obeys |c u_N mathsf_X|<epsilon_term. Therefore a signed endpoint-plus-four-moment estimate |Re[g+s_*'U]|>A_L+2epsilon_term is sufficient for the existing Abel gap, but is not proved here.
```

This is a theorem target, not a proved inequality.

## Symmetry Guard

The symmetric kernel is not positive semidefinite:

```text
If sqrt(N)<p<=N is prime, then for every j=1,2,3 the {1,p} principal minor of K_(j,t,N) is [[0,log(p)^j/2],[log(p)^j/2,0]], with determinant -log(p)^(2j)/4<0. Bertrand's postulate supplies such a p for every N>=2, with N=2,3,4 checked directly.
```

The machine witness records 8 cutoffs and
24 strict negative minors. The sign is
independent of the heat parameter because `log(1)=0`.
These are negative prime-edge minors of the actual symmetric kernel.

There is a second guard: the sum uses `q_d^(0)q_m^(0)`, not
`q_d^(0)conjugate(q_m^(0))`. Even a positive real kernel would not make
the real part of this complex bilinear form positive.

## Relation To The Older Symmetry

```text
The exact first-jet symmetry is multiplicative and hyperbolically truncated: dm<=N with additive weight Lambda(d)+Lambda(m). The legacy radial square uses d^2+m^2 and product weight Lambda(d)Lambda(m). The latter is a different, second-order object and is not the Xi contact scalar.
```

The images were useful in pointing at symmetry, but the exact symmetry
that reaches the Xi contact problem is multiplicative and first order.
The radial `Lambda(d)Lambda(m)` square is a separate second-order object.

## Route Decision

```text
Retain the symmetric Mangoldt form as a Type-I/II or bilinear cancellation coordinate. Do not promote transpose symmetry, entrywise positivity, or the untruncated heat Gram to a contact lower bound. A successful estimate must keep the complex phase, hyperbolic cutoff, recurrent endpoint, first-order coefficients, W_0=0, and adjacent charts joined.
```

The promising use of this normal form is a cancellation-preserving
Vaughan/Type-I/II decomposition of the three symmetric moments, composed
with the endpoint before absolute values. Another positivity search on
the symmetric matrix is ruled out by the exact prime-edge minor.

## Boundary

This gate proves the exact quadratic correction polynomial, four-moment Abel collapse, symmetric Mangoldt divisor-pair identity, heat coupling, endpoint-complete contact normal form, and indefinite prime-edge guard. It proves no signed endpoint-plus-moment lower bound, Abel-scalar gap, horizontal successor winding cap, q<1 or bounded-L closure, contact exclusion, Lambda<=0, PF-infinity, RH, or prize-level conclusion.
