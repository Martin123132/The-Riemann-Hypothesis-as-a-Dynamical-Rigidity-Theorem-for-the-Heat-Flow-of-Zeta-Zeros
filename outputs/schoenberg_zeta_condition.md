# Schoenberg/Zeta Condition and Coefficient Normalization Audit

Date: 2026-07-03

Status: theorem-map correction and finite diagnostic evidence. This is not a proof of RH and not a rigorous interval certificate.

## 1. Exact Theorem Shape To Keep Separate

Schoenberg's total-positivity route is not a statement about the direct heat kernel moments by themselves.

The relevant chain in Groechenig's exposition is:

```text
Polya frequency function Lambda
  <=> 1/Psi is a two-sided Laplace transform of Lambda
  <=> Psi is in the Laguerre-Polya class, with Psi(0) > 0.
```

Read backwards, this says:

```text
Psi is Laguerre-Polya
  <=> Fourier/Laplace transform of 1/Psi is a Polya frequency function.
```

For zeta, Groechenig's RH-equivalent condition is:

```text
RH
  <=> inverse Fourier transform of 1/xi(1/2 + tau) is a Polya frequency function.
```

There is also a one-sided version using the even factorization:

```text
Xi(s) = Xi_1(-s^2).

RH
  <=> 1/Xi_1 is a Laplace transform of a one-sided Polya frequency function.
```

This is an inverse-function/kernel condition. It is not automatically the same as total positivity of the original Newman heat kernel `Phi_lambda`.

## 2. Coefficient Normalization Issue

Our moment expansion has the form:

```text
mu_2k(lambda) = integral_R u^(2k) exp(lambda*u^2) Phi(u) du
A_k(lambda) = mu_2k(lambda) * k! / (2k)!
```

For the even entire function:

```text
H_lambda(z) = integral_R exp(i*z*u) Phi_lambda(u) du
```

the even expansion can be written as:

```text
H_lambda(z) = F_lambda(-z^2)
F_lambda(s) = sum_{k >= 0} mu_2k(lambda) s^k / (2k)!
            = sum_{k >= 0} A_k(lambda) s^k / k!
```

So `A_k(lambda)` are exponential/Jensen coefficients, usually denoted `beta_k` in Jensen-polynomial statements.

The ordinary Taylor coefficients of `F_lambda` are instead:

```text
c_k(lambda) = A_k(lambda) / k!
             = mu_2k(lambda) / (2k)!.
```

This matters because the ASW/Edrei PF-sequence theorem for the restricted Laguerre-Polya class is stated for ordinary Taylor coefficients `c_k`, not for exponential/Jensen coefficients `A_k`.

## 3. Consequence For The Old Toeplitz Obstruction

The bundle's earlier conclusion was:

```text
negative Toeplitz minors damage the simple PF-infinity coefficient route.
```

That is true for the tested sequence `A_k`.

But it does not by itself rule out the ASW/Edrei coefficient route, because that route should test:

```text
c_k = A_k / k!
```

not `A_k`.

This is a significant correction to the theorem map.

## 4. Finite Exact-Rational Audit

New script:

```text
work/rh_compute/scripts/toeplitz_pf_audit.py
```

Source cache:

```text
work/rh_compute/results/repro_hankel_15c_coefficients.jsonl
dps = 90
n_sum = 100
cutoff = 8
```

The script builds the upper-triangular Toeplitz matrix:

```text
T[i,j] = c[j-i] if j >= i, otherwise 0.
```

It rationalizes the cached coefficients to 70 digits and computes finite minors exactly over QQ.

### Lambda = 0, N = 10, Order <= 4

Run:

```text
toeplitz_pf_beta_taylor_N10_o4
```

Results:

```text
sequence beta = A_k:
  tests: 60,625
  positive: 15,101
  zero: 40,935
  negative: 4,589

sequence taylor = A_k/k!:
  tests: 60,625
  positive: 19,690
  zero: 40,935
  negative: 0
```

The old reported obstruction is exactly reproduced for `A_k`:

```text
rows = 0,1,2
cols = 1,2,3
det(beta minor) = -1.71119821042992898147325e-11
```

For the ordinary Taylor sequence, the same minor is positive:

```text
det(taylor minor) = 7.59800565413986900638776e-11
```

Smallest nonzero minor seen in this finite `taylor` pass:

```text
rows = 0,1,2,3
cols = 6,7,8,9
det approximately 9.68569470272923019694506e-75
positive
```

This is a warning: the finite pass is sign-clean, but high-order minors can be extremely small.

### Lambda = 0, Wider Finite Checks

```text
taylor, N = 10, order <= 5:
  tests: 124,129
  negative: 0

taylor, N = 12, order <= 4:
  tests: 297,925
  negative: 0
```

### Lambda Grid, N = 10, Order <= 4

For `taylor = A_k/k!`, no negative minors were found for:

```text
lambda = 0
lambda = 1e-6
lambda = 1e-4
lambda = 1e-2
lambda = 1e-1
```

Each lambda-grid run tested 60,625 finite minors through order 4.

### Arb Propagation Probe

A follow-up Arb probe was added:

```text
work/rh_compute/scripts/arb_toeplitz_pf_probe.py
```

For `taylor = A_k/k!`, `lambda = 0`, `N = 10`, orders `<= 4`, and assumed absolute coefficient radii on `A_k` of `1e-80`, `1e-90`, and `1e-100`, the probe found:

```text
negative intervals: 0
positive intervals: 19,690
exact Arb zeros: 37,635
inconclusive zero-containing intervals: 3,300
```

The exact-rational point audit has 40,935 zero minors in the same range, so the 3,300 inconclusive rows appear to be interval dependency/overestimation around point-zero minors. This requires a structural-zero classifier or a dependency-preserving determinant method before it can become a certificate.

Update:

The structural-zero classifier has now been added. It expands the Toeplitz determinant combinatorially in formal symbols `a_d` and detects minors whose determinant polynomial is identically zero. With this front gate, the same Arb run becomes:

```text
positive intervals: 19,690
structural zeros: 40,935
negative intervals: 0
inconclusive zero-containing intervals: 0
```

The same structural Arb method is clean for:

```text
taylor, lambda = 0, N = 10, orders <= 5
taylor, lambda = 0, N = 12, orders <= 4
taylor, lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}, N = 10, orders <= 4
```

The remaining certificate gap is no longer zero-minor dependency in these ranges; it is proving that the input coefficient balls enclose the exact analytic coefficients.

## 5. Updated Route

The naive statement:

```text
A_k is PF-infinity.
```

is false for the tested data.

The corrected statement to investigate is:

```text
c_k = A_k/k! is PF-infinity.
```

If proved for the exact coefficient sequence of `F_0 = Xi_1`, this would place `F_0` in the restricted Laguerre-Polya class and would therefore imply RH through the classical coefficient theorem.

That means this route is potentially aligned with a known theorem ecosystem, but the missing work is enormous:

```text
finite rationalized minors
  -> all exact Toeplitz minors of c_k are nonnegative
  -> F_0 is restricted Laguerre-Polya
  -> RH.
```

At present, only the first diagnostic arrow has evidence.

## 6. Gates

Do not claim:

```text
the Toeplitz route is dead
```

without specifying that the negative minors were for `A_k`, not for `A_k/k!`.

Do not claim:

```text
c_k is PF-infinity
```

until there is an all-order proof or a rigorous theorem reducing the infinite check to certified finite/structural data.

Do not claim:

```text
Schoenberg proves our kernel is totally positive
```

because the direct Schoenberg zeta condition concerns the inverse transform of `1/Xi` or `1/Xi_1`, not the original heat kernel moments.

## 7. Next Tasks

1. Add interval-enclosed coefficient generation for `c_k = mu_2k/(2k)!`.
2. Replace assumed Arb coefficient radii with interval-enclosed moment computation.
3. Re-run Toeplitz/PF audits with proved coefficient balls and the structural-zero front gate.
4. Search ASW/Edrei/Karlin variants for sufficient conditions that prove PF-infinity from structured moment or recurrence data.
5. Compare the `c_k` Toeplitz minors against Turan and higher-order log-concavity inequalities.
6. Decide whether the heat-flow deformation preserves or destroys PF-infinity for `c_k(lambda)` as lambda varies.

## Sources

- Groechenig, "Schoenberg's Theory of Totally Positive Functions and the Riemann Zeta Function": https://arxiv.org/abs/2007.12889
- Aissen, Schoenberg, Whitney, "On the generating functions of totally positive sequences I" (classical ASW theorem, cited by Groechenig).
- Edrei's theorem on totally positive sequences and generating functions, as referenced in the total-positivity literature.
