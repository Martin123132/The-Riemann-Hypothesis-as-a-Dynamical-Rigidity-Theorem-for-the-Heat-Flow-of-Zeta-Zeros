# Delayed Cooperative Heat-Tail Lemma

Date: 2026-07-16

Status: exact shifted heat-flow propagation lemma. This is not a proof of RH or `Lambda <= 0`.

## Lemma

Fix an order `m` and a starting shift `n0`. Assume the lower layer is
positive throughout the heat interval, the fixed-order eventual tail
holds, and `Q_(m,n)(-100)>0` for every `n>=n0`.

Choose `N=max(n0,N_m)`. The eventual-tail theorem supplies positivity
for every `n>=N`. At `n=N-1`, variation of constants has a positive
initial term and a positive forcing term from `Q_(m,N)`. Repeat at
`n=N-2,...,n0`. This is a finite descending induction.

```text
Q_(m,n)(lambda)=E_(m,n)(lambda)*(Q_(m,n)(-100)+integral_(-100)^lambda E_(m,n)(s)^(-1)*a_(m,n)(s)*Q_(m,n+1)(s)ds), E_(m,n)>0
[Q_(m-1,n)(lambda)>0 for every n>=n0 on -100<=lambda<=0, the fixed-order m eventual tail holds, and Q_(m,n)(-100)>0 for every n>=n0] => [Q_(m,n)(lambda)>0 for every n>=n0 and -100<=lambda<=0]
```

The proof uses no `Q_(m,n)` with `n<n0`.

## Order Ten

The completed order-nine layer gives `a_(10,n)>0`, so setting `n0=4` gives

```text
[Q_(10,n)(-100)>0 for every integer n>=4] implies Q_(10,n)(lambda)>0 for every n>=4 and -100<=lambda<=0
```

This remains conditional on proving the delayed endpoint premise. It
does not assign signs to shifts `0,1,2,3` at negative lambda and is not
PF-infinity, RH, or `Lambda<=0`.

## Order Eleven

Applying the same theorem at `m=11,n0=4` gives

```text
[Q_(10,n)(lambda)>0 for every n>=4 and -100<=lambda<=0 and Q_(11,n)(-100)>0 for every n>=4] implies Q_(11,n)(lambda)>0 for every n>=4 and -100<=lambda<=0
```

This second specialization is conditional on the order-ten heat ray
and the order-eleven endpoint ray; it does not assert either premise.
