# Jensen-Window PF Newman Polymath-15 Critical Transversality Target

Date: 2026-07-17

Status: exact corrected `C1` reduction of the positive-boundary
collision problem. This is not a proof of `Lambda <= 0` or RH.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_critical_transversality_target.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_critical_transversality_target.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_critical_transversality_target.py
```

## Double-Zero Reduction

On the real axis write

```text
H_t=A_t*Z_t with A_t>0; H_t'=A_t*((log A_t)'*Z_t+Z_t')
H_t=H_t'=0 if and only if Z_t=Z_t'=0
```

Thus the Newman endgame does not require a global lower bound for the
full Laguerre curvature. It is enough to keep the normalized value and
first derivative from vanishing simultaneously.

## Corrected First Jet

Use the refined Polymath-15 split

```text
Z_t=J_(N,t)+r_ref, J_(N,t)=P_(N,t)-Q_(N,t), P_(N,t)=2Re(exp(i*beta_t)D_(N,t))
exp(i*beta)D=X+iY, exp(i*beta)D'=U+iV1, B=beta'
J=2X-Q, J'=2*(U-BY)-Q'
```

The scale-adapted distance from a double zero is exactly

```text
T_L[J]=(2X-Q)^2+((2*(U-BY)-Q')/L)^2
```

This is positive by construction and avoids the indefinite second-jet
terms in the corrected curvature formula.

## C1 Exclusion Lemma

Suppose the refined remainder satisfies

```text
|r_ref|<=epsilon_0, |r_ref'|<=L*epsilon_1
```

At an exact double zero `Z=Z'=0`, so `J=-r_ref` and `J'=-r_ref'`.
Consequently

```text
At a double zero, T_L[J]=r_ref^2+(r_ref'/L)^2 <=epsilon_0^2+epsilon_1^2
```

Therefore the concrete live target is

```text
T_L[J]>epsilon_0^2+epsilon_1^2 uniformly on the critical region
```

A radius-`1/L` collar now needs only one Cauchy derivative:

```text
A radius-1/L collar transfers the refined scalar remainder through one derivative with one factor L; no second-derivative remainder is needed
```

## Global Composition

```text
Combine the oscillatory-zeta theorem for every fixed tL>=c_*+epsilon at sufficiently large L, the critical transversality target on the residual high-frequency layer 0<tL<=c_*+o(1), and compact no-double-zero certificates for every bounded-L remainder, where c_*=4911678521/1933561194
If all three regions are closed, positive-boundary attainment gives Lambda<=0, hence RH
```

The Lehmer stress gate explains why this first-jet target is preferable
to a fixed `kappa*L^2` curvature floor: close pairs make the exact
curvature margin genuinely small. The lower bound above remains open
and is the new proof-facing obligation, not an established sign theorem.
