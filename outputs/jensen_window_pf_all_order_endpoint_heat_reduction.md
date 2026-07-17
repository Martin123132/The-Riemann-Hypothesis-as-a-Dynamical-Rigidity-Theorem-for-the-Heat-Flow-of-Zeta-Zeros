# Jensen-Window PF All-Order Endpoint-To-Heat Reduction

Date: 2026-07-16

Status: exact all-order endpoint-to-heat equivalence, with the signed
base completed through order nine and the candidate continuation
rejected by a rigorous order-ten endpoint counterexample. This is not
a proof of the separate Jensen-window PF bridge, RH, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_all_order_endpoint_heat_reduction.json
python work/rh_compute/scripts/jensen_window_pf_all_order_endpoint_heat_reduction.py
python work/rh_compute/scripts/check_jensen_window_pf_all_order_endpoint_heat_reduction.py
```

## Coordinates And Eventual Tail

Put

```text
H_(m,n)=det[A_(n+i+j)]_(0<=i,j<m), H_(0,n)=1, epsilon_m=(-1)^binom(m,2), Q_(m,n)=epsilon_m*H_(m,n)
```

The graded-kernel Vandermonde theorem has the correctly ordered
quantifiers

```text
for every fixed m exists N_m such that for every n>=N_m and -100<=lambda<=0, Q_(m,n)(lambda)>0
```

The threshold may depend on `m`. That is sufficient below because the
argument fixes one finite `m` before choosing and using `N_m`.

## Arbitrary-Order Heat Algebra

Let `C_j=(A_(n+i+j))_(0<=i<m)`. Under the sequence-shift derivation,
multilinearity gives

```text
delta(H_(m,n))=det[C_0,...,C_(m-2),C_m], C_j=(A_(n+i+j))_(0<=i<m)
```

Every earlier shifted column duplicates its right neighbor. Splitting
`4(n+i+j)+2` into a column part and a row part leaves only the final
shifted column and final shifted row:

```text
column part=[4*(n+m-1)+2]*delta(H); row part=4*(m-1)*delta(H)
Q_(m,n)'=c_(m,n)*delta(Q_(m,n)), c_(m,n)=4*n+8*m-6
```

For the `m` by `(m+1)` matrix with columns `C_0,...,C_m`, apply the
three-term flag-Plucker relation to its maximal minors and the maximal
minors of its first `m-1` rows. With the Hankel index map this is

```text
Q_(m-1,n+1)*delta(Q_(m,n))=Q_(m-1,n)*Q_(m,n+1)+delta(Q_(m-1,n+1))*Q_(m,n)
```

Since `c_(m-1,n+1)=c_(m,n)-4`, division by the positive preceding
layer yields

```text
Q_(m,n)'=a_(m,n)*Q_(m,n+1)+b_(m,n)*Q_(m,n), a_(m,n)=c_(m,n)*Q_(m-1,n)/Q_(m-1,n+1)>0, b_(m,n)=c_(m,n)/(c_(m,n)-4)*(log Q_(m-1,n+1))'
```

Thus the order-`m` system is one-sided cooperative whenever order
`m-1` has already been completed.

## Tail-To-Endpoint Propagation

Fix `m`. Choose its own `N_m` from the eventual-tail theorem. For
`n=N_m-1,...,0`, variation of constants reads

```text
Q_(m,n)(lambda)=E_(m,n)(lambda)*(Q_(m,n)(-100)+integral_(-100)^lambda E_(m,n)(s)^(-1)*a_(m,n)(s)*Q_(m,n+1)(s)ds), E_(m,n)>0
```

The endpoint term is strictly positive. Backward induction in `n`
makes the integrand nonnegative because `a_(m,n)>0` and the next shift
is positive throughout the interval. Hence

```text
[Q_(k,n)>0 on the heat interval for every k<m,n and Q_(m,n)(-100)>0 for every n] => [Q_(m,n)(lambda)>0 for every n and -100<=lambda<=0]
```

There is no exchange of `forall m` and `exists N_m` in this proof.
Ordinary induction completes one order before moving to the next.

## Exact All-Order Reduction

The fixed-order programme now supplies

```text
Q_(m,n)(lambda)>0 for 1<=m<=9, every n>=0, and -100<=lambda<=0
```

Iterating the single-layer lemma from `m=10` proves the exact
equivalence

```text
[Q_(m,n)(-100)>0 for every m>=10,n>=0] iff [Q_(m,n)(lambda)>0 for every m>=10,n>=0,-100<=lambda<=0]
```

The reverse implication is simply restriction to `lambda=-100`. Once
the contiguous hierarchy is complete, the fixed-order initial-minor
criterion applies separately at each finite `m` and each `lambda`,
giving

```text
all-order contiguous interval positivity => all-order consecutive-row arbitrary-column signed-Hankel positivity
```

## Rejected Endpoint Antecedent

Desnanot-Jacobi and the orientation identity
`binom(m,2)+binom(m-2,2)-2binom(m-1,2)=1` give

```text
Q_(m,n)*Q_(m-2,n+2)=Q_(m-1,n+1)^2-Q_(m-1,n)*Q_(m-1,n+2)
Q_(m,n)>0 iff Q_(m-1,n+1)^2>Q_(m-1,n)*Q_(m-1,n+2), provided Q_(m-2,n+2)>0
```

The conditional all-order route would require the static statement

```text
Q_(m,n)(-100)>0 for every integer m>=10 and n>=0
```

equivalently the continuation from order nine of the strict
log-concavity hierarchy in the positive signed coordinates. The
order-ten counterexample now proves that this antecedent is false:

```text
Q_(10,n)(-100)<0 for n=0,1,2,3
```

Thus the endpoint/interval equivalence remains exact, but it cannot
serve as a proof route for the actual sequence.

## Separate Bridge

Abstract all-order shifted signed-Hankel positivity would not
automatically establish Jensen hyperbolicity:

```text
all-order shifted signed-Hankel positivity does not by itself prove PF-infinity of every binomially weighted Jensen window
```

The binomially weighted Jensen-window problem remains open, but its
replacement Xi/Phi-specific theorem cannot assume the rejected
all-shift signed-Hankel hierarchy.

## Machine Audit

The artifact records `255` exact orientation-parity checks
and independent symbolic determinant specializations at `4`
orders. Those finite checks audit the implementation; arbitrary-order
validity comes from the determinant cancellation and the two general
determinant identities proved above.
