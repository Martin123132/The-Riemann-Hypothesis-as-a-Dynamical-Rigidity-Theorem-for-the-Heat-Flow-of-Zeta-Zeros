# Jensen-Window PF Suzuki Jordan-Error Kernel Reduction

Date: 2026-07-23

Status: exact arithmetic reduction with one open cancellation gate. This is
not a proof of an eventual sign, RH, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_suzuki_jordan_error_kernel_reduction.json
python work/rh_compute/scripts/jensen_window_pf_suzuki_jordan_error_kernel_reduction.py
python work/rh_compute/scripts/check_jensen_window_pf_suzuki_jordan_error_kernel_reduction.py
```

## Cumulative Jordan Mass

For `0<omega<1/2`,

```text
c_omega(n)
 =sum_(d*m=n)mu(d)d^(-omega)m^omega,

C_omega(x)=sum_(n<=x)c_omega(n),

A_omega=1/((1+omega)*zeta(1+2omega)),

E_omega(x)=C_omega(x)-A_omega*x^(1+omega).           (JEK.1)
```

The constant `A_omega` is the residue main-term coefficient coming from
`zeta(s-omega)/zeta(s+omega)` at `s=1+omega`.

For Suzuki's level-`k` weight, set

```text
W_(omega,k)(u)=-u*g_(omega,k)'(u).
```

The logarithmic smoothing recursion gives

```text
W_(omega,k)
 =g_(omega,k-1)+(1/2)g_(omega,k).                    (JEK.2)
```

Stieltjes summation, using `g_(omega,k)(1)=0`, yields

```text
H_(omega,k)(x)
 =x^(-1/2)*integral_0^x
  C_omega(t)W_(omega,k)(t/x)dt/t.                   (JEK.3)
```

## Exact Main-Term Annihilation

The archimedean Mellin transform of `g_(omega,k)` contains the factor
`s-omega-1`. At `s=1+omega`,

```text
integral_0^1 g_(omega,k)(u)u^omega du=0.
```

Integration by parts, with vanishing endpoint products, then gives

```text
integral_0^1 W_(omega,k)(u)u^omega du=0.             (JEK.4)
```

Therefore the entire positive residue main term in (JEK.1) contributes
exactly zero to (JEK.3):

```text
H_(omega,k)(x)
 =x^(-1/2)*integral_0^1
  E_omega(x*u)W_(omega,k)(u)du/u.                   (JEK.5)
```

This is the sharp arithmetic handoff. Eventual positivity cannot come from
ordinary Jordan-totient main-term domination; it is wholly a signed
summatory-error cancellation problem.

## Exact Mobius Decomposition

Put

```text
P_omega(y)=sum_(m<=y)m^omega,
R_omega(y)=P_omega(y)-y^(1+omega)/(1+omega).
```

Using the coefficient convolution and the absolutely convergent identity
`sum mu(d)d^(-(1+2omega))=1/zeta(1+2omega)` gives

```text
E_omega(x)
 =sum_(d<=x)mu(d)d^(-omega)R_omega(x/d)
  -x^(1+omega)/(1+omega)
   *sum_(d>x)mu(d)d^(-(1+2omega)).                   (JEK.6)
```

Both terms require Mobius cancellation and neither has a fixed sign.

If `R_omega(y)=-y^(1+omega)/(1+omega)` is retained for `0<y<1`, the
tail is absorbed and (JEK.6) is equivalently the full generalized Muntz
identity

```text
E_omega(x)
 =sum_(d>=1)mu(d)d^(-omega)R_omega(x/d).
```

The finite natural-dilation energy of this full form is developed in
`outputs/jensen_window_pf_jordan_muntz_causal_energy_bridge.md`.

## Absolute-Bound Barrier

The elementary power-sum estimate

```text
R_omega(y)=O_omega(y^omega)
```

and absolute summation in (JEK.6) give

```text
E_omega(x)=O_omega(x^(1-omega)).                     (JEK.7)
```

Every finite smoothing level has

```text
W_(omega,k)(u)=O_(omega,k)(u^(omega-1))
```

near zero. Substitution into (JEK.5), split at `u=1/x`, yields only

```text
H_(omega,k)(x)
 =O_(omega,k)(x^(1/2-omega)*(1+log x)).              (JEK.8)
```

This is far above the RH-conditional scale `1` for `k=1` or
`(log x)^(k-1)` for higher `k`. Increasing `k` does not improve the power in
this naive absolute estimate because the near-zero exponent of every
smoothed weight is the same.

## Surviving Target

For one explicit sequence `omega_j->0` and convenient smoothing levels
`k_j`, prove directly that

```text
integral_0^1
 E_(omega_j)(x*u)W_(omega_j,k_j)(u)du/u
```

has one eventual sign after the `x^(-1/2)` normalization. A successful
argument must use cancellation or order structure beyond (JEK.7); positive
residue-main-term domination and absolute-value estimates are ruled out.

## Sources

- Masatoshi Suzuki, `On monotonicity of certain weighted summatory functions associated with L-functions`: https://arxiv.org/abs/1204.1823
- `outputs/jensen_window_pf_suzuki_cofinal_monotonicity_hierarchy.md`

## Proof Boundary

Equations (JEK.1)-(JEK.8) are exact identities or elementary bounds. They
identify the signed Mobius-error convolution that must be controlled and
reject two inadequate proof routes. They do not prove the required
cancellation, eventual sign, RH, or Lambda<=0.
