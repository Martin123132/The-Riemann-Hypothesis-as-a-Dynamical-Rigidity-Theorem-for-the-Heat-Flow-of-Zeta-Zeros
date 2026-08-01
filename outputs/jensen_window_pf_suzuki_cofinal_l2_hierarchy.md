# Jensen-Window PF Suzuki Cofinal L2 Hierarchy

Date: 2026-07-23

Status: internally audited theorem-candidate reduction with one open
arithmetic energy gate. This is not a proof of an `L2` residual, RH, or
`Lambda <= 0`. Independent expert review is required.

```text
work/rh_compute/results/jensen_window_pf_suzuki_cofinal_l2_hierarchy.json
python work/rh_compute/scripts/jensen_window_pf_suzuki_cofinal_l2_hierarchy.py
python work/rh_compute/scripts/check_jensen_window_pf_suzuki_cofinal_l2_hierarchy.py
```

## Central Subtraction

For `0<omega<1/2`, put `z=s-1/2` and

```text
Q_omega(z)=xi(1/2+z-omega)/xi(1/2+z+omega).
```

The functional equation gives `Q_omega(0)=1`. Suzuki's smoothing hierarchy
has the Laplace form

```text
integral_0^infinity H_(omega,k)(exp(t))*exp(-z*t)dt
 =Q_omega(z)/z^k,                                    (SL2.1)
```

initially for `Re(z)>1/2+omega`.

Write the Taylor expansion at the central point as

```text
Q_omega(z)=sum_(j>=0)q_(omega,j)z^j,
T_(omega,k-1)(z)=sum_(j=0)^(k-1)q_(omega,j)z^j,

P_(omega,k)(t)
 =sum_(j=0)^(k-1)
  q_(omega,j)*t^(k-1-j)/(k-1-j)!.
```

Since the Laplace transform of `t^m/m!` is `z^(-m-1)`, the residual

```text
r_(omega,k)(t)
 =H_(omega,k)(exp(t))-P_(omega,k)(t)
```

satisfies the exact identity

```text
Laplace(r_(omega,k))(z)
 =[Q_omega(z)-T_(omega,k-1)(z)]/z^k.                (SL2.2)
```

At `k=1`, the subtraction is simply `P_(omega,1)=1`. At higher levels
the leading term is `t^(k-1)/(k-1)!`.

The first broadened target is completely explicit. The functional equation
gives

```text
q_(omega,1)
 =-2*xi'(1/2+omega)/xi(1/2+omega),

P_(omega,2)(t)=t+q_(omega,1),

r_(omega,2)(t)
 =-q_(omega,1)
  +integral_0^t (H_(omega,1)(exp(u))-1)du.            (SL2.2a)
```

Thus `k=2` asks for square integrability of a normalized cumulative
level-one discrepancy, not merely for ordinary smoothing.

## Fixed-Shift Equivalence

Suppose first that `r_(omega,k)` belongs to `L2(0,infinity)`. The
half-plane Paley-Wiener theorem makes its Laplace transform an `H2`
function in `Re(z)>0`. Equation (SL2.2), initially valid farther right,
then gives

```text
Q_omega(z)
 =T_(omega,k-1)(z)+z^k*Laplace(r_(omega,k))(z).
```

Thus the reduced quotient has no pole in the right half-plane. Boundary
unimodularity, Suzuki's high-strip estimate, and Phragmen-Lindelof promote
it to an inner function.

Conversely, suppose `Q_omega` is inner. Then

```text
F_(omega,k)(z)
 =[Q_omega(z)-T_(omega,k-1)(z)]/z^k
```

is bounded near `z=0` by Taylor cancellation. On the imaginary boundary
it is `O(1/|t|)` at infinity because `Q_omega` is bounded and the Taylor
polynomial has degree `k-1`. Hence `F_(omega,k)` belongs to `H2`.
Paley-Wiener and transform uniqueness recover the actual residual in
`L2(0,infinity)`. Therefore, for every fixed `omega` and every `k>=1`,

```text
r_(omega,k) in L2(0,infinity)
 iff Q_omega is inner
 iff D(omega).                                       (SL2.3)
```

This extends the published level-one `L2` criterion to the complete
logarithmic smoothing hierarchy.

## Cofinal Criterion

The fixed-omega phase theorem says that a sequence of shifts in `D`
accumulating at zero forces RH. Under RH every shifted quotient is inner.
Combining this with (SL2.3) gives the sharpened equivalence

```text
RH
iff there exist omega_j->0 and integers k_j>=1 such that
    H_(omega_j,k_j)(exp(t))-P_(omega_j,k_j)(t)
    belongs to L2(0,infinity) for every j.            (SL2.4)
```

The smoothing order may vary with the shift. This creates an energy target
parallel to, and weaker in shape than, eventual one-sign behavior.

## Arithmetic Form

The exact Jordan-error reduction turns the residual into

```text
x^(-1/2)*integral_0^1
 E_omega(x*u)W_(omega,k)(u)du/u
 -P_(omega,k)(log x).                                (SL2.5)
```

Thus the surviving arithmetic task is to prove, for one explicit
`omega_j->0` and convenient `k_j`, that (SL2.5) has finite squared norm
against `dx/x`, without assuming a zero-free half-plane. At level one this
is

```text
integral_1^infinity |H_(omega,1)(x)-1|^2 dx/x
 <infinity.                                          (SL2.6)
```

At level two it is the concrete condition

```text
integral_1^infinity |
 x^(-1/2)*integral_0^1 E_omega(x*u)W_(omega,2)(u)du/u
 -log(x)+2*xi'(1/2+omega)/xi(1/2+omega)
 |^2 dx/x < infinity.                                (SL2.7)
```

## Countermodel Guard

For `a>0`, the right-half-plane inner function

```text
Q_good(z)=(a-z)/(a+z)
```

has Taylor coefficients `q_0=1` and
`q_j=2*(-1)^j/a^j` for `j>=1`. Exact subtraction gives

```text
[Q_good(z)-T_(k-1)(z)]/z^k
 =2*(-1)^k/[a^(k-1)*(a+z)],
```

whose inverse Laplace transform is a decaying exponential in `L2`.

In contrast,

```text
Q_bad(z)=(a+z)/(a-z)
```

also has modulus one on the imaginary axis, but has a pole at `z=a`.
Its actual positive-time residual contains `exp(a*t)` and is not in `L2`.

There is a sharper warning. The regularized boundary trace

```text
B_(omega,k)(y)
 =[Q_omega(i*y)-T_(omega,k-1)(i*y)]/(i*y)^k
```

is already in `L2(R)` without RH: Taylor cancellation controls `y=0` and
boundary unimodularity gives an `O(1/|y|)` tail. For `Q_bad`, its inverse
Fourier transform is the anti-causal decaying function

```text
2*a^(-(k-1))*exp(a*t)*1_(t<0).
```

For `Q_good`, the inverse Fourier transform is supported on `t>=0`.
Paley-Wiener therefore makes positive-time support, not finite boundary
energy, the decisive condition. A Plancherel calculation of the boundary
norm alone is automatic and cannot prove (SL2.3).

One fixed shift is also insufficient because exact horizontal zero
cancellation can hide off-line zeros. A finite numerical energy integral
cannot certify the infinite tail.

## Riesz-Smoothing Comparison

Das, Lang, Wan, and Xu study positive logarithmic Riesz means of the
classical Euler-totient sum. Their square-root bound after two smoothing
operations assumes RH and a reciprocal-zeta-derivative moment bound. That
result is useful contour technology, but it neither identifies Suzuki's
signed archimedean kernel nor proves (SL2.5). In Suzuki's transform the
Jordan residue main term is annihilated exactly and the central Taylor
polynomial supplies the target profile.

## Sources

- Masatoshi Suzuki, `On monotonicity of certain weighted summatory functions associated with L-functions`: https://arxiv.org/abs/1204.1823
- Masatoshi Suzuki, `A canonical system of differential equations arising from the Riemann zeta-function`: https://arxiv.org/abs/1204.1827
- Sanjana Das, Hannah Lang, Hamilton Wan, and Nancy Xu, `The Distribution of Error Terms of Smoothed Summatory Totient Functions`: https://arxiv.org/abs/2207.07722
- `outputs/jensen_window_pf_suzuki_fixed_omega_phase_diagram.md`
- `outputs/jensen_window_pf_suzuki_jordan_error_kernel_reduction.md`

## Proof Boundary

Equations (SL2.1)-(SL2.3) use published Suzuki identities, standard
Paley-Wiener theory, and the internally audited fixed-shift innerness
handoff. Equation (SL2.4) is a new corpus theorem candidate requiring
independent review. No residual in (SL2.5) or (SL2.6) is proved
unconditionally, so this artifact does not prove RH or `Lambda <= 0`.
