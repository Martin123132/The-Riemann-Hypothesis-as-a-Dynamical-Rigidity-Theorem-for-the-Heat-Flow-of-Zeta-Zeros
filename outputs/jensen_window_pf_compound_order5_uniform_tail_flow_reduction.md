# Jensen-Window PF Compound Order-Five Uniform Tail And Flow Reduction

Date: 2026-07-13

Status: exact uniform eventual order-five tail and conditional forward
reduction with one open `lambda=-100` entry. This is not a proof of
all-shift order five, PF-infinity, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_compound_order5_uniform_tail_flow_reduction`.

```text
work/rh_compute/results/jensen_window_pf_compound_order5_uniform_tail_flow_reduction.json
python work/rh_compute/scripts/jensen_window_pf_compound_order5_uniform_tail_flow_reduction.py
python work/rh_compute/scripts/check_jensen_window_pf_compound_order5_uniform_tail_flow_reduction.py
```

## Compact-Heat Asymptotic

O'Sullivan's all-order suitable-multiplier theorem applies to the same
compact family through every fixed order needed here. The exact
suitability audit and Lambert recurrence now run through order eleven:

```text
Delta_k^m log(A_k^(1)(-T)/A_k^(1)(0))=O(log(k)/k^m), m=2,...,11, uniformly for 0<=T<=100
all fixed local differences of log(A_k(-T)/A_k^(1)(-T)) are O_p,m(k^-p) uniformly
```

Nine-node Newton interpolation transfers those estimates to the graded
ratio coefficients. For `M=n+8` and `h=Delta(M)^2`, the exact
120-permutation determinant expansion is

```text
[h^0,...,h^10] det K=[0,0,0,0,0,0,0,0,0,0,294912*G_2^10]
294912*G2^10*h^10
```

The coefficient is the positive Vandermonde constant
`2^10*(1!*2!*3!*4!)=294912`; `G_3,...,G_11` cancel from the first
nonzero term. Since `G_2->1` uniformly,

```text
there exists N_5 such that H_(5,n)(lambda)>0 for every n>=N_5 and -100<=lambda<=0
```

## Exact Entry Coordinate

Desnanot-Jacobi gives

```text
H_(5,n)*H_(3,n+2)=H_(4,n)*H_(4,n+2)-H_(4,n+1)^2
```

The completed signs `H_(3,n)<0` and `H_(4,n)>0` make the order-five
entry target exactly

```text
H_(4,n+1)(-100)^2>H_(4,n)(-100)*H_(4,n+2)(-100), every n>=0.
```

Thus the new analytic problem is strict log-concavity of the actual
positive order-four determinant sequence at `lambda=-100`.

## Cooperative Flow

Exact affine-Hankel differentiation and the adjacent Plucker identity
give

```text
Q_n'=a_n*Q_(n+1)+b_n*Q_n, a_n=(4*n+34)*H_(4,n)/H_(4,n+1)>0, b_n=((4*n+34)/(4*n+30))*(log H_(4,n+1))'
```

The off-diagonal coefficient is positive because order four is already
complete. The flow needs no order-six sign. If the open entry target is
proved, the uniform tail confines the remaining system to finitely many
indices, and variation of constants gives

```text
[H_(5,n)(-100)>0 for all n] => [H_(5,n)(lambda)>0 for all n and -100<=lambda<=0]
```

## Lower-Layer Kill Gate

Order five is not a formal consequence of the lower signed layers. The
positive rational sequence

```text
['1', '1', '1/2', '1/6', '1/24', '1/120', '1/720', '1/5040', '1/42000']
```

has every available signed contiguous minor through order four strict,
but

```text
H_(5,0)=-1/3657830400000<0,
H_(4,1)^2-H_(4,0)*H_(4,2)=-1/3792438558720000000<0.
```

## Open Input

```text
prove H_(5,n)(-100)>0 for every integer n>=0
equivalently prove strict log-concavity of H_(4,n)(-100) for every n>=0
```

This is now the sole order-five propagation input. It is not supplied by
finite evidence, the lower layers, or the eventual tail.

```text
outputs/jensen_window_pf_compound_order4_uniform_heat_forward_invariance_certificate.md
outputs/jensen_window_pf_order4_noncontiguous_total_positivity_transfer.md
outputs/signed_hankel_jensen_bridge_target.md
outputs/formal_core.md
```
