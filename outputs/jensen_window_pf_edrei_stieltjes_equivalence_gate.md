# Jensen-Window PF Edrei-Stieltjes Equivalence Gate

Date: 2026-07-23

Status: exact endpoint-equivalence and theorem-route correction. This is not
a proof of the all-order Stieltjes property, coefficient PF-infinity, Jensen
hyperbolicity for zeta, RH, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_edrei_stieltjes_equivalence_gate.json
python work/rh_compute/scripts/jensen_window_pf_edrei_stieltjes_equivalence_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_edrei_stieltjes_equivalence_gate.py
```

Current result:

```text
validated Jensen-window PF Edrei-Stieltjes equivalence gate: 14 rows, 0 issues, 12 exact indexing checks, 7 exact Hankel checks, 1 unified endpoint, 3 finite/nonpromotion guards, 1 open Xi/Phi handoff
```

## Exact Indexing

```text
H_lambda(z)=F_lambda(z)/F_lambda(0)=sum_(k>=0)d_k(lambda)z^k, H_lambda(0)=1, and H_lambda is entire
log H(z)=sum_(n>=1)ell_n*z^n, q_n=n*ell_n, p_n=(-1)^(n-1)*q_n
a_r:=(-1)^r*[z^r](H'(z)/H(z))=p_(r+1)
```

The last line is the important off-by-one check: the sequence tested by the
existing moment-recurrence scout, `a_r=p_(r+1)`, is exactly Sokal's signed
Taylor-coefficient sequence for `H'/H`.

## Classical Equivalence

Sokal's Proposition 6 applies because `H` is entire and `H(0)=1`:

```text
For entire H with H(0)!=0: H in LP+ <=> (a_r)_(r>=0) is a Stieltjes moment sequence
H'(z)/H(z)=integral_[0,infinity)dnu(x)/(1+x*z)=gamma+sum_j m_j*beta_j/(1+beta_j*z)
```

Entireness is decisive. It forces the representing measure to be the unique
discrete zero measure

```text
a_r=integral_[0,infinity)x^r*dnu(x), nu=gamma*delta_0+sum_j m_j*beta_j*delta_(beta_j), gamma>=0, beta_j>0, m_j in Z_(>=1), sum_j m_j*beta_j<infinity
```

so no extra integrality or product-reconstruction conjecture is needed after
the all-order Stieltjes property has actually been proved.

Combining this with the coefficient-PF/Jensen gate gives

```text
c is PF-infinity <=> H is LP+ <=> a_r=p_(r+1) is Stieltjes <=> every shifted Jensen window is finite PF-infinity
```

This is an equivalence of endpoint targets, not evidence that the endpoint
holds.

## Exact Hankel Target

The Stieltjes moment criterion is

```text
a is Stieltjes <=> (a_(i+j))_(i,j>=0) and (a_(i+j+1))_(i,j>=0) are positive semidefinite on every finite block
```

A clean strict sufficient target in the corpus notation is therefore

```text
D_(m,1)=det(p_(i+j+1))_(i,j=0)^m>0 and D_(m,2)=det(p_(i+j+2))_(i,j=0)^m>0 for every m>=0 => a is Stieltjes
```

Only the two columns `s=1,2` are structurally required. The wider finite
staircase remains useful stress evidence, but it is not an all-order proof.
The exact continued-fraction alternative is

```text
a is Stieltjes <=> sum_(r>=0)a_r*t^r has a formal Stieltjes continued fraction with nonnegative coefficients, including terminating cases
```

## Existing Finite Evidence

```text
320 finite Edrei-log sign rows pass
4205 finite power-Hankel rows pass
55 finite recurrence rows through order 12 are positive
orders 13..20 are interval-inconclusive, not negative
```

These are finite pieces of an endpoint-equivalent Stieltjes target. They remain
diagnostics and cannot be promoted by extrapolation.

## Xi/Phi Handoff

The moment integral gives

```text
F_lambda(z)=integral_R exp(lambda*u^2)*Phi(u)*cosh(u*sqrt(z))*du
H_lambda'(z)/H_lambda(z)=[integral_R exp(lambda*u^2)*Phi(u)*u*sinh(u*sqrt(z))/(2*sqrt(z))*du]/[integral_R exp(lambda*u^2)*Phi(u)*cosh(u*sqrt(z))*du]
```

A noncircular proof that this ratio is a Stieltjes function at `lambda=0`
would prove the common endpoint. Positivity of `Phi`, or the fact that every
fixed-scale `cosh(u*sqrt(z))` lies in `LP+`, is insufficient: positive
mixtures need not preserve the required zero or logarithmic-derivative class.

## Object-Separation Guard

```text
Hankel matrices of a_r=p_(r+1) are nonlinear Edrei-log objects and are not the original signed-Hankel matrices of A_k
```

Thus the known order-ten failure in the original signed-Hankel hierarchy does
not refute this nonlinear Edrei-log Hankel target. It only forbids conflating
the two objects.

## Machine Audit

The generator and independent checker use a rational finite Edrei product to
verify the indexing and both Hankel columns:

```text
indexing checks: 12
Hankel determinant checks: 7
s=1 signs by sizes 1..4: ['286/315', '424/59535', '160/992436543', '0']
s=2 signs by sizes 1..3: ['358/3969', '200/15752961', '0']
```

The terminal zeros are the expected finite-support rank termination, not
failures of the Stieltjes property.

## Primary Source

Sokal, Proposition 6 and Proposition 7:

```text
https://doi.org/10.1016/j.jmaa.2022.126432
https://discovery.ucl.ac.uk/10150527/1/Sokal_1-s2.0-S0022247X22004462-main.pdf
```

## Boundary

This artifact proves an exact equivalence of all-order target formulations, identifies a strict two-column Hankel sufficient target, and records an exact Xi/Phi logarithmic-derivative handoff. It does not prove the all-order Stieltjes property, coefficient PF-infinity, Jensen hyperbolicity for zeta, RH, or Lambda <= 0.
