# Jensen-Window PF Burnol Cell-Energy/Tail Obstruction

Date: 2026-07-23

Status: exact cell-energy and discrete-discrepancy reduction with one
reciprocal-zeta tail obstruction and one open short-interval gate. This is
not a proof of RH or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_burnol_cell_energy_tail_obstruction.json
python work/rh_compute/scripts/jensen_window_pf_burnol_cell_energy_tail_obstruction.py
python work/rh_compute/scripts/check_jensen_window_pf_burnol_cell_energy_tail_obstruction.py
```

## Arithmetic Cells

Fix `0<omega<1/2` and define

```text
A_(omega,N)
 =sum_(d<=N)mu(d)d^(-(1+2omega)),

b_(omega,N)(n)
 =sum_(d|n,d<=N)mu(d)d^(-2omega),

B_(omega,N)(k)
 =sum_(n<=k)b_(omega,N)(n)
 =sum_(d<=N)mu(d)d^(-2omega)floor(k/d).             (BCET.1)
```

For Burnol's finite fractional-part function,

```text
f_(2omega,N)(t)
 =sum_(d<=N)mu(d)d^(-2omega){1/(d*t)},
```

the reciprocal coordinate `x=1/t` gives

```text
f_(2omega,N)(1/x)
 =A_(omega,N)x-B_(omega,N)(floor(x)).               (BCET.2)
```

Consequently,

```text
||t^(-omega)f_(2omega,N)||_2^2
 =sum_(k>=0) integral_k^(k+1)
  x^(2omega-2)|A_(omega,N)x-B_(omega,N)(k)|^2 dx.   (BCET.3)
```

The zero cell is exactly

```text
A_(omega,N)^2/(1+2omega).
```

For `k>=1`, write `Delta_a(k)=(k+1)^a-k^a`. Direct integration gives

```text
J_(omega,N)(k)
 =A_N^2*Delta_(1+2omega)(k)/(1+2omega)
  -2*A_N*B_N(k)*Delta_(2omega)(k)/(2omega)
  +B_N(k)^2*Delta_(2omega-1)(k)/(2omega-1).         (BCET.4)
```

Every `J_(omega,N)(k)` is a nonnegative cell integral even though the three
expanded terms need not be nonnegative separately.

## Discrete Energy

Define the summatory discrepancy

```text
S_(omega,N)(k)
 =A_(omega,N)k-B_(omega,N)(k)
 =sum_(d<=N)mu(d)d^(-2omega){k/d}.                  (BCET.5)
```

On `x=k+y`, the cell error is simply

```text
A_N*x-B_N(k)=S_N(k)+A_N*y,  0<=y<1.
```

Because `2omega-2` lies in `(-2,-1)`, the cell weight is comparable to
`k^(2omega-2)`, uniformly in `k>=1`. Also

```text
sup_N |A_(omega,N)|<=zeta(1+2omega).
```

The elementary inequalities

```text
|S+A*y|^2<=2S^2+2A^2,

S^2<=2|S+A*y|^2+2A^2*y^2
```

therefore prove the exact uniform-boundedness equivalence

```text
sup_N ||t^(-omega)f_(2omega,N)||_2<infinity

iff

sup_N Q_(omega,N)<infinity,                         (BCET.6)

Q_(omega,N)
 =sum_(k>=1)k^(2omega-2)S_(omega,N)(k)^2.
```

This is the discrete short-multiplicative-interval form of the natural
mollifier gate.

## Forced Reciprocal-Zeta Rate

Absolute convergence gives

```text
A_(omega,infinity)=1/zeta(1+2omega).
```

For every `k<=N`, all divisors of every `n<=k` already occur, so

```text
B_(omega,N)(k)=B_(omega,infinity)(k)
```

and hence

```text
S_(omega,N)(k)-S_(omega,infinity)(k)
 =k*(A_(omega,N)-A_(omega,infinity)).               (BCET.7)
```

If `Q_(omega,N)<=C` uniformly, Fatou first gives
`Q_(omega,infinity)<=C`. Taking the weighted `l2` norm of (BCET.7) on
`1<=k<=N` and using the triangle inequality gives

```text
|A_(omega,N)-1/zeta(1+2omega)|
 *sqrt(sum_(k<=N)k^(2omega))
 <=2*sqrt(C).
```

Since

```text
sum_(k<=N)k^(2omega)
 >=N^(1+2omega)/(1+2omega),
```

the uniform norm forces

```text
|A_(omega,N)-1/zeta(1+2omega)|
 <=2*sqrt(C*(1+2omega))*N^(-1/2-omega).             (BCET.8)
```

Thus the norm target already contains a square-root reciprocal-zeta tail
estimate before the remaining short-interval cancellation is considered.

## Zero-Free Consequence

Put `s_0=1+2omega` and

```text
r_N=A_(omega,infinity)-A_(omega,N).
```

The bound `r_N=O(N^(-1/2-omega))`, followed by Abel summation, makes

```text
sum_(n>=1)mu(n)n^(-s)
```

converge analytically for

```text
Re(s)>s_0-(1/2+omega)=1/2+omega.                    (BCET.9)
```

It agrees with `1/zeta(s)` where `Re(s)>1`, so analytic continuation makes
the open half-plane in (BCET.9) zero-free. No boundary-line assertion is
made. On a cofinal `omega_j->0` sequence, these open half-planes already
force RH.

This gives a direct arithmetic proof of the hard direction of the uniform
partial-norm criterion. It also shows why a generic Hilbert-space or
total-positivity estimate cannot close the problem: it would have to imply
the explicit reciprocal-zeta rate (BCET.8).

## Open Gate

The remaining target can be stated entirely discretely:

```text
sup_N sum_(k>=1) k^(2omega-2)
 |sum_(d<=N)mu(d)d^(-2omega){k/d}|^2
 <infinity                                           (BCET.10)
```

for every shift in one explicit cofinal sequence. A proof must supply both
the scalar tail rate (BCET.8) and cancellation in the weighted
short-multiplicative-interval discrepancies. Neither component is proved
here.

## Proof Boundary

Equations (BCET.1)-(BCET.9) are finite reindexings, elementary integrals,
weighted `l2` estimates, or Abel-summation consequences.
Equation (BCET.10) is open. The artifact exposes an RH-strength necessary condition;
it does not establish that condition, RH, PF-infinity, or `Lambda <= 0`.
