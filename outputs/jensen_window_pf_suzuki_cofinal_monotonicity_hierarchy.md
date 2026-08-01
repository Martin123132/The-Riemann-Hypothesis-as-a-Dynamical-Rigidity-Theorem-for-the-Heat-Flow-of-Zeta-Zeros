# Jensen-Window PF Suzuki Cofinal Monotonicity Hierarchy

Date: 2026-07-23

Status: internally audited theorem-candidate reduction with one open
arithmetic gate. This is not a proof of eventual sign, RH, or
`Lambda <= 0`. Independent expert review is required.

```text
work/rh_compute/results/jensen_window_pf_suzuki_cofinal_monotonicity_hierarchy.json
python work/rh_compute/scripts/jensen_window_pf_suzuki_cofinal_monotonicity_hierarchy.py
python work/rh_compute/scripts/check_jensen_window_pf_suzuki_cofinal_monotonicity_hierarchy.py
```

## Published Hierarchy

For `0<omega<1/2`, Suzuki defines

```text
c_omega(n)
 =n^omega*sum_(d|n)mu(d)d^(-2omega)
 =n^omega*product_(p|n)(1-p^(-2omega))>0.
```

Starting from his explicit beta-integral weight `g_(omega,1)`, put

```text
g_(omega,k+1)(x)
 =integral_x^1 sqrt(y/x)*g_(omega,k)(y)dy/y,

H_(omega,k)(x)
 =x^(-1/2)*sum_(n<=x)c_omega(n)*g_(omega,k)(n/x).
```

The smoothing is a logarithmic antiderivative:

```text
H_(omega,k+1)(x)=integral_1^x H_(omega,k)(y)dy/y.
```

Suzuki's Mellin identity is

```text
integral_1^infinity H_(omega,k)(x)x^(1/2-s)dx/x
 =[xi(s-omega)/xi(s+omega)]/(s-1/2)^k,              (CMH.1)
```

initially for `Re(s)>1+omega`. At level `k=1`,

```text
H_(omega,1)(x)=sqrt(x)*h_omega^<1>(x)
```

in the canonical-system notation used by the fixed-omega phase diagram.

## Cofinal Equivalence

If one `H_(omega,k)` has a single sign for every sufficiently large `x`,
Landau's theorem applied to (CMH.1), together with the absence of real xi
zeros, removes every uncancelled pole of the reduced quotient in
`Re(s)>1/2`. Suzuki's high-strip estimate and Phragmen-Lindelof then give
meromorphic innerness. Therefore

```text
eventual one-sign of H_(omega,k) => D(omega).        (CMH.2)
```

The fixed-omega phase theorem makes cofinality decisive:

```text
If omega_j decreases to zero and H_(omega_j,k_j) is eventually
one-signed for every j, for any integers k_j>=1, then RH.       (CMH.3)
```

Conversely, Suzuki proves under RH that, for each fixed
`0<omega<1/2`,

```text
H_(omega,1)(x)=1+O_omega(x^(-B_omega))
```

for some `B_omega>0`, while every fixed `k>=2` has a positive leading
term `(log x)^(k-1)/(k-1)!`. Thus the sharpened equivalence is

```text
RH
iff there exist omega_j->0 and integers k_j>=1 such that every
    H_(omega_j,k_j) has one eventual sign.            (CMH.4)
```

Suzuki's published statement asks for all shifts in an interval. The
weakening to one cofinal sequence is the new corpus consequence and depends
on the internally audited fixed-omega phase theorem.

## Signed-Weight Guard

Higher smoothing does not make coefficient positivity a termwise proof.
For `0<omega<1/2`, the explicit first weight satisfies

```text
g_(omega,1)(x)~-a_omega*x^(omega-1),  a_omega>0,
```

as `x->0+`. The recursion gives, for every `k>=1`,

```text
g_(omega,k)(x)
 ~-a_omega*(1/2-omega)^(-(k-1))*x^(omega-1)
```

near zero. Near one, `g_(omega,k)` is positive and vanishes like a
positive constant times `(1-x)^(omega+k-1)`. Every smoothing weight is
therefore genuinely signed. The positive arithmetic coefficients cannot
close (CMH.4) term by term.

## Surviving Target

The broadened arithmetic obligation is

```text
Choose an explicit omega_j->0 and any convenient smoothing orders k_j>=1.
Prove eventual one-sign behavior of H_(omega_j,k_j) for every j by direct
arithmetic inequalities, without assuming a zero-free half-plane.
```

This is weaker than proving the unsmoothed `k=1` target at every shift.
The `k=2` level is the first natural focus because it is one logarithmic
antiderivative of the sampled function while retaining the same Landau
consequence.

Finite positivity remains nonpromotable, and one fixed shift remains
insufficient because exact horizontal zero cancellation can hide off-line
zeros.

## Sources

- Masatoshi Suzuki, `On monotonicity of certain weighted summatory functions associated with L-functions`: https://arxiv.org/abs/1204.1823
- Masatoshi Suzuki, `A canonical system of differential equations arising from the Riemann zeta-function`: https://arxiv.org/abs/1204.1827
- `outputs/jensen_window_pf_suzuki_fixed_omega_phase_diagram.md`

## Proof Boundary

The published hierarchy and Landau implication are exact source-backed
inputs. The cofinal weakening is an internally audited theorem candidate
obtained by composing them with the fixed-omega phase theorem. No eventual
sign is proved at any subcritical cofinal family, and this artifact does not
prove RH or Lambda<=0.
