# Jensen-Window PF Factorial Multiplier Split Audit

Date: 2026-07-06

Status: finite theorem-search diagnostic. This is not a proof of
Jensen-window PF-infinity, Jensen hyperbolicity, Laguerre-Polya
membership, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_factorial_multiplier_split_audit`.

Proof boundary: this artifact separates a valid conditional
Pólya-Schur multiplier step from the missing input theorem. It does not
prove the sign-regular-to-Jensen transfer theorem.

Machine-readable artifact:

```text
work/rh_compute/results/jensen_window_pf_factorial_multiplier_split_audit.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_factorial_multiplier_split_audit.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_factorial_multiplier_split_audit.py
```

Current result:

```text
validated Jensen-window PF factorial multiplier split audit: 5 exact rows, 315 raw degree-2 anti-hyperbolic rows, 315 normalized degree-2 positive rows, 0 ready-to-apply rows, 0 issues
```

## Exact Split

For the raw moments and normalized Jensen coefficients:

```text
A_k = gamma_k * mu_{2k}
gamma_k = k!/(2*k)!
sum gamma_k*z^k/k! = cosh(sqrt(z))
```

The shifted sequences `gamma_{n+j}` are conditional multiplier
sequences because they come from derivatives of `cosh(sqrt(z))`.
Thus they preserve real-rootedness if the raw shifted moment-window
polynomials are already real-rooted.

But raw moment-window hyperbolicity is the wrong input theorem:

```text
Delta_2(M_{2,n}) = 4*(mu_{2*n+2}^2 - mu_{2*n}*mu_{2*n+4}) <= 0
```

Cauchy-Schwarz gives raw moment log-convexity, which is the opposite
degree-2 sign from real-rooted positive-coefficient Jensen windows.

The factorial normalization changes the degree-2 threshold to:

```text
A_{n+1}^2 >= A_n*A_{n+2}
iff mu_{2*n+2}^2/(mu_{2*n}*mu_{2*n+4}) >= (2*n+1)/(2*n+3)
```

## Finite Arb Sanity Check

Using the existing coefficient enclosures:

```text
lambdas: 0.0, 0.000001, 0.0001, 0.01, 0.1
max coefficient index: 64
raw degree-2 rows: 315
raw degree-2 anti-hyperbolic rows: 315
normalized degree-2 rows: 315
normalized degree-2 positive rows: 315
normalized threshold positive rows: 315
```

Sharpest finite margins:

```text
max raw discriminant: -3.575201093265539228E-80 at lambda=0.0, n=62
min normalized discriminant: 2.277193819536730908E-328 at lambda=0.0, n=62
min normalized threshold margin: 7.510922618285697771E-3 at lambda=0.0, n=62
```

## Route Consequence

The factorial multiplier route is not dead, but it is conditional. It
must be paired with a genuinely normalized Xi/Phi-specific theorem,
such as a ratio-cone argument, determinant integral, positive kernel,
or Schur-positive specialization. Raw moment positivity is not enough.

The subsequent reciprocal-gamma mixture audit proves the complete
fixed-scale sign-regularity theorem and isolates the remaining failure:
positive scale mixtures are not sign-regularity preserving, even at
order two. For Xi, the completed order-two ratio cone is exactly the
tilted relative-variance bound `CV_n^2<=2/(2n+1)`; higher compound
concentration remains open.

Integration:

```text
outputs/coefficient_pf_bridge_obstruction.md
outputs/jensen_window_pf_heat_flow_boundary_threshold_lemma.md
outputs/jensen_window_pf_sign_regular_transfer_gap_matrix.md
outputs/jensen_window_pf_theorem_machinery_fit_matrix.md
outputs/jensen_window_pf_reciprocal_gamma_mixture_sign_gate.md
```

Summary:

The factorial multiplier gamma_k=k!/(2*k)! is not the missing bridge by itself. It gives a clean conditional preservation step for shifted windows, but the natural raw moment-window input theorem is anti-aligned at degree 2. The useful route is therefore a normalized Xi/Phi-specific ratio-cone, determinant integral, or Schur-positive theorem, not raw moment positivity.
