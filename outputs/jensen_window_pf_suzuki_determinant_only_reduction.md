# Jensen-Window PF Suzuki Determinant-Only Reduction

Date: 2026-07-23

Status: internally audited exact reduction and theorem candidate.
It is not a proof of RH or `Lambda <= 0`.
It sharpens the logical form of Suzuki's criterion but does not establish the
still-open all-time determinant premise. Independent expert review remains
required before publication.

```text
work/rh_compute/results/jensen_window_pf_suzuki_determinant_only_reduction.json
python work/rh_compute/scripts/jensen_window_pf_suzuki_determinant_only_reduction.py
python work/rh_compute/scripts/check_jensen_window_pf_suzuki_determinant_only_reduction.py
```

## Statement

For zeta, Suzuki's Theorem 2.4 remains equivalent to RH after
condition (4), the terminal canonical-kernel limit, is deleted:

```text
RH iff there exist strictly decreasing omega_n->0 and positive integers nu_n with nu_n*omega_n>1 such that det(I+/-K_(omega_n,nu_n)[t])!=0 for every n and every t>=0
```

The reverse implication is not a proof that these determinants
are nonzero. It is a reduction showing that this one all-time
operator condition would already be enough.

## 1. From Every Truncation To One Global Form

The spectral-frontier lemma gives

```text
det(I+/-K[t])!=0 for every t>=0 implies ||K[t]||<1 for every finite t
```

On the dense subspace of compactly supported L2 functions define

```text
D=L2_c(R)={compactly supported L2 functions}; B(f,g)=integral_R integral_R K(x+y)f(y)conj(g(x)) dy dx
```

For any fixed `f,g`, choose `t` above both supports. Then

```text
For f,g in D and t above both supports: B(f,g)=<K[t]f,g>, hence |B(f,g)|<=||f||_2||g||_2
```

The bound is independent of the chosen support cutoff. Therefore

```text
B extends uniquely to a bounded Hermitian form on L2(R), represented by a self-adjoint H with ||H||<=1
```

No limit of the pointwise strict constants is claimed; their
common non-strict upper bound one is exactly what the extension
requires.

## 2. Reflection Produces A Causal Multiplier

Let `J` be reflection. A change of variables gives

```text
Jf(x)=f(-x), G=H*J; for f in D, (Gf)(x)=integral_R K(x-y)f(y)dy=(K*f)(x)
```

The equality holds first against compact tests. Since both sides
represent the same distribution and `G` is bounded on L2, it
identifies Suzuki's locally defined kernel with the global
convolution distribution. Its support now supplies

```text
G commutes with every translation and P_a G=P_a G P_a for P_a=1_(-infinity,a); thus G is bounded, translation-invariant, and causal
```

The standard L2 Paley-Wiener multiplier theorem therefore gives

```text
There is M in H-infinity({Re(s)>0}) with ||M||_infinity=||G||<=1 and L(Gf)(s)=M(s)L(f)(s) for causal f
```

This is the step that must not be replaced by an unjustified
real-axis contour shift: boundedness is obtained from all finite
forms before the multiplier theorem is invoked.

## 3. Identify The Arithmetic Transfer Function

Suzuki proves absolute convergence of the kernel transform in a
high half-plane. There, ordinary Fubini calculation gives

```text
For Re(s)>c: M(s)=integral_0^infinity K(x)e^(-s*x)dx=Theta_(omega,nu)(i*s)
```

A concrete division guard is

```text
For q=1_[0,1], L(q)(s)=(1-e^(-s))/s!=0 when Re(s)>0, so the high-strip multiplier identity can be divided by L(q)(s)
```

Rotating the transfer function to the upper half-plane yields

```text
Theta_tilde(z)=M(-i*z) belongs to H-infinity(C+) with ||Theta_tilde||_infinity<=1 and agrees with [xi(1/2-omega-i*z)/xi(1/2+omega-i*z)]^nu on Im(z)>c
```

The two meromorphic functions agree on a nonempty open strip, so

```text
The meromorphic identity theorem makes every apparent pole of the xi quotient in C+ removable
```

This conclusion derives inner-type analyticity from the all-time
finite-section bound; it does not assume HB, RH, or a contour
shift across unknown poles.

## 4. Pole Cancellation Becomes A Zero Shift

Let `rho=beta+i*gamma` be a zero of xi with
`beta>1/2+omega`. The corresponding denominator zero occurs at

```text
z=-gamma+i*(beta-1/2-omega) in C+.
```

Pole removal then gives

```text
If xi(rho)=0 and Re(rho)>1/2+omega, pole removal at z=-Im(rho)+i*(Re(rho)-1/2-omega) forces xi(rho-2*omega)=0
```

The exact coordinate check uses `rho=3/4+7i`,
`omega=1/8`, and `z=-7+i/8`. The denominator
argument is `3/4+7i` and the numerator
argument is `1/2+7i`, a real shift of
`1/4=2*omega`.

For one omega this is only cancellation, not a contradiction.
For a strictly decreasing cofinal sequence, however,

```text
If omega_n decreases strictly to 0, an off-line zero rho would force distinct zeros rho-2*omega_n -> rho for all large n, contradicting isolation of zeros of nonzero entire xi
```

Thus xi has no zeros with real part greater than one half. Its
functional equation excludes the reflected left half, proving
the reverse implication in the displayed equivalence. Conversely,
RH gives Suzuki's HB hypothesis and hence all-time strict
contractivity and determinant nonvanishing.

## Terminal Condition

The resulting logical consequence is

```text
Within Suzuki's cofinal criterion, the terminal condition J_(omega_n,nu_n)(t;z,z)->0 follows after the determinant conditions imply RH; it is not an additional logical premise
```

This is sequence-level redundancy. It is not the stronger and
generally unjustified claim that one fixed all-time determinant
family directly forces its terminal limit. The guard is

```text
For one fixed omega, holomorphic continuation permits a denominator zero to be canceled by a numerator zero at rho-2*omega; cofinal distinct shifts are essential
```

## Surviving Proof Obligation

The sharpened Suzuki route now has one arithmetic gate:

```text
The reduction does not establish the antecedent: one still must prove ||K_(omega_n,nu_n)[t]||<1 for every finite t on a cofinal omega_n sequence without assuming RH or HB
```

Finite grids, the unconditional local interval, Hilbert-Schmidt
domination, and an RH/HB-derived isometry do not close this gate.

## Sources

- Masatoshi Suzuki, `Hamiltonians arising from L-functions in the Selberg class`, Theorems 2.3-2.4 and Propositions 4.1-4.4 in arXiv v3: https://arxiv.org/abs/1606.05726
- Chris Guiver, Hartmut Logemann, and Mark R. Opmeer, `Operator-valued multiplier theorems for causal translation-invariant operators with applications to control theoretic input-output stability`, classical L2 equivalence in equations (1.1) and (1.3): https://link.springer.com/article/10.1007/s00498-024-00387-4
- `outputs/jensen_window_pf_suzuki_spectral_frontier.md`

## Proof Boundary

This artifact proves, modulo the cited standard causal-multiplier theorem, that a cofinal sequence of Suzuki all-time determinant conditions implies RH and hence makes Suzuki's terminal J-limit premise redundant in that cofinal criterion. It does not prove that the determinant conditions hold for zeta, does not rule out common-zero cancellation for one fixed omega, and does not prove RH or Lambda<=0. The determinant-only equivalence is a new corpus theorem candidate and should receive independent expert review before being represented as a published result.
