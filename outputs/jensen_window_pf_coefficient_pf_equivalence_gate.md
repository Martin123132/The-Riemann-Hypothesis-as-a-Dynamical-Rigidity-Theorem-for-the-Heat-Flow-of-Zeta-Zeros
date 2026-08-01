# Jensen-Window PF Coefficient-PF Equivalence Gate

Date: 2026-07-23

Status: exact theorem-equivalence and proof-route correction. This is not
a proof of coefficient PF-infinity, the signed-Hankel-to-PF bridge,
Laguerre-Polya membership, RH, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_coefficient_pf_equivalence_gate.json
python work/rh_compute/scripts/jensen_window_pf_coefficient_pf_equivalence_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_coefficient_pf_equivalence_gate.py
```

Current result:

```text
validated Jensen-window PF coefficient-PF equivalence gate: 13 rows, 0 issues, 3 exact coefficient identities, 4 classical/closure steps, 1 seven-way equivalence, 3 guards, 1 open structural handoff
```

## Exact Normalization

For the corpus normalization,

A_k>=0 for every k, A_0>0, and F is entire; all hold for the corpus zeta coefficient sequence.

```text
c_k=A_k/k! and F(z)=sum_(k>=0)c_k*z^k=sum_(k>=0)A_k*z^k/k!
F_n(z)=F^(n)(z)=sum_(j>=0)A_(n+j)*z^j/j!
```

Define the diagonal tail operator by

```text
T_n[z^j]=A_(n+j)*z^j
T_n[(1+z)^d]=P_(d,n)(z)=sum_(j=0)^d C(d,j)*A_(n+j)*z^j
```

Thus the binomially weighted Jensen polynomial is exactly the
Polya-Schur test polynomial for the shifted sequence.

## Exact Equivalence

The Polya-Schur theorem, ASW/Edrei characterization, closure under
differentiation, and the finite ASW theorem give

```text
c is PF-infinity <=> F is Laguerre-Polya type I <=> A is a multiplier sequence <=> every P_(d,0) is hyperbolic <=> every shifted A_(n+j) is a multiplier sequence <=> every P_(d,n) is hyperbolic <=> every B^(d,n) is finite PF-infinity
```

The all-order coefficient-PF route and the all-degree/all-shift
Jensen-window route are therefore not two independent endpoint
theorems. They are exact formulations of the same endpoint theorem.

## Programme Correction

What remains separate is the present finite evidence. Finite Toeplitz
certificates for `c_k` do not prove the infinite PF statement, and finite
Jensen-window certificates do not prove every degree and shift.

Likewise, the equivalence supplies no transfer from the signed-Hankel
certificates to coefficient PF-infinity. That structural implication is
still the live bridge, now with one unambiguous endpoint target rather
than two nominally separate ones.

## Positive-Mixture Guard

The positive Phi moment integral can be viewed as a positive mixture of
dilations. Coefficient positivity or positive mixing does not make that
diagonal operator a hyperbolicity preserver. Polya-Schur requires the
operator's multiplier-sequence/type-I hypothesis, which is an endpoint
statement and cannot be smuggled in through positivity language.

## Machine Audit

The generator and independent checker verify

```text
normalization identities: 12
derivative-tail identities: 42
diagonal-window identities: 80
```

These finite exact checks audit the indexing and factorial normalization.
The arbitrary-order equivalence comes from the cited classical theorems.

## Primary Sources

```text
https://annals.math.princeton.edu/wp-content/uploads/annals-v170-n1-p14-p.pdf
https://doi.org/10.1007/BF02786970
https://doi.org/10.1007/BF02786971
https://doi.org/10.1073/pnas.37.5.303
```

## Boundary

Passing this checker proves only that the target formulations and guards
are stated consistently. It does not prove the all-order antecedent for
the zeta coefficients and therefore does not prove RH or `Lambda <= 0`.
