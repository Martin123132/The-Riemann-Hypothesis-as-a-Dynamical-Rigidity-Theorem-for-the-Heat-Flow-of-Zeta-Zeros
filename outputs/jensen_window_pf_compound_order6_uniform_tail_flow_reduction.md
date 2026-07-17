# Jensen-Window PF Compound Order-Six Uniform Tail And Flow Reduction

Date: 2026-07-13

Status: exact uniform eventual signed order-six tail and conditional
forward reduction with one open lambda=-100 entry. This is not a proof
of all-shift order six, PF-infinity, RH, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_compound_order6_uniform_tail_flow_reduction.json
python work/rh_compute/scripts/jensen_window_pf_compound_order6_uniform_tail_flow_reduction.py
python work/rh_compute/scripts/check_jensen_window_pf_compound_order6_uniform_tail_flow_reduction.py
```

## Universal Signed Tail

The published all-order Xi ratio expansion, a compact heat-tilt audit
through order sixteen, and superpolynomial higher-theta suppression
reduce the normalized six-by-six determinant to

```text
K_(i,j)(h)=exp(-sum_(m>=2)G_m*h^(m-1)*(10-i-j)^m)
[h^0,...,h^15] det K=[0,...,0,-1132462080*G_2^15].
```

The exact audit checks `720` permutations
against `684` weighted monomials,
or `492480` exact permutation-monomial
terms. Every contribution containing `G_3,...,G_16` cancels from the
first nonzero coefficient. Since `Q_6=-H_6`, this proves

```text
Q_(6,n)(lambda)=positive_scale(lambda,n)*(1132462080*G_2(lambda,n)^15*h(lambda,n)^15+o(h(lambda,n)^15)) uniformly
there exists N_6 such that Q_(6,n)(lambda)=-H_(6,n)(lambda)>0 for every n>=N_6 and -100<=lambda<=0
```

The threshold is finite but non-effective.

## Fixed-Order Recursion

Put `epsilon_m=(-1)^binom(m,2)` and `Q_(m,n)=epsilon_m*H_(m,n)`.
Desnanot-Jacobi and the orientation identity
`epsilon_m*epsilon_(m-2)=-1` give

```text
Q_(m,n)*Q_(m-2,n+2)=Q_(m-1,n+1)^2-Q_(m-1,n)*Q_(m-1,n+2)
Q_(m,n)>0 iff Q_(m-1,n+1)^2>Q_(m-1,n)*Q_(m-1,n+2), provided Q_(m-2,n+2)>0
```

Thus every new signed contiguous layer is strict log-concavity of the
preceding signed layer. The affine heat derivative and adjacent Plucker
identity similarly give the general cooperative recursion

```text
Q_(m,n)'=c_(m,n)*delta(Q_(m,n)), c_(m,n)=4*n+8*m-6
Q_(m-1,n+1)*delta(Q_(m,n))=Q_(m-1,n)*Q_(m,n+1)+delta(Q_(m-1,n+1))*Q_(m,n)
Q_(m,n)'=a_(m,n)*Q_(m,n+1)+b_(m,n)*Q_(m,n), a_(m,n)=c_(m,n)*Q_(m-1,n)/Q_(m-1,n+1)>0, b_(m,n)=c_(m,n)/(c_(m,n)-4)*(log Q_(m-1,n+1))'
```

The positive off-diagonal coefficient uses only the completed order
`m-1` sign; no order `m+1` sign is required.

## Order-Six Specialization

At order six the endpoint coordinate and flow are

```text
H_(6,n)*H_(4,n+2)=H_(5,n)*H_(5,n+2)-H_(5,n+1)^2
Q_(6,n)=-H_(6,n)>0 iff H_(5,n+1)^2>H_(5,n)*H_(5,n+2)
Q_n'=a_n*Q_(n+1)+b_n*Q_n, a_n=(4*n+42)*H_(5,n)/H_(5,n+1)>0, b_n=((4*n+42)/(4*n+38))*(log H_(5,n+1))'
```

The completed `H_5>0` layer makes `a_n>0`. The uniform signed tail
then confines a hypothetical first loss to finitely many shifts, so
variation of constants proves the conditional theorem

```text
[Q_(6,n)(-100)>0 for all n] => [Q_(6,n)(lambda)>0 for all n and -100<=lambda<=0]
completed signed contiguous layers through order six imply every arbitrary-column signed layer through order six
```

## Countermodel Gate

The rational sequence

```text
(1,1,1/2,1/6,1/24,1/120,1/720,1/5040,1/40320,1/362880,1/3590000)
```

has every available strict signed contiguous minor through order five,
but

```text
H_(6,0)=3247/6004150801779916800000000>0,
Q_(6,0)=-3247/6004150801779916800000000<0,
H_(5,1)^2-H_(5,0)H_(5,2)=-3247/8784866131573981512204288000000000000<0.
```

Order six is therefore genuinely new and cannot be promoted from the
completed lower cone.

## Remaining Endpoint Target

The remaining fixed-order problem is all-shift signed order-six entry:

```text
Q_(6,n)(-100)>0 for every integer n>=0, equivalently H_(5,n+1)(-100)^2>H_(5,n)(-100)*H_(5,n+2)(-100)
```

A separate rigorous prefix certificate handles the currently available
finite coefficient range. The analytic all-shift tail remains open.

```text
outputs/jensen_window_pf_compound_order5_uniform_heat_forward_invariance_certificate.md
outputs/jensen_window_pf_uniform_first_summand_heat_tilt_asymptotic_theorem.md
outputs/jensen_window_pf_uniform_superpolynomial_first_summand_dominance.md
outputs/signed_hankel_jensen_bridge_target.md
outputs/formal_core.md
```
