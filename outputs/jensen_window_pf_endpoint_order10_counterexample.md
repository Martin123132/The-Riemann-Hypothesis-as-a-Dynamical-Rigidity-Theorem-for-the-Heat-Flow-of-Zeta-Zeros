# Endpoint Order-Ten Counterexample

Date: 2026-07-16

Status: rigorous first-open-order counterexample. This rejects the
proposed all-order deep rectangular endpoint hierarchy; it is not a
counterexample to RH or Jensen hyperbolicity. This is not a proof of
RH or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_endpoint_order10_counterexample.json
python work/rh_compute/scripts/jensen_window_pf_endpoint_order10_counterexample.py
python work/rh_compute/scripts/check_jensen_window_pf_endpoint_order10_counterexample.py
```

## Exact Coordinate

```text
epsilon_10=(-1)^45=-1, Q_(10,n)=-H_(10,n)
Q_(10,n)*Q_(8,n+2)=Q_(9,n+1)^2-Q_(9,n)*Q_(9,n+2)
L_n=Q_(9,n+1)^2/(Q_(9,n)*Q_(9,n+2))-1
inside the positive order-eight/order-nine cone, Q_(10,n)>0 iff L_n>0
Q_(10,n)(-100)=A_0(-100)^10*s_(((n+9)^10))(h)
```

Orders eight and nine are already positive at every endpoint shift.
Therefore the sign of `L_n` is exactly the sign of `Q_(10,n)`.

## Rigorous Failure

Direct `4096`-bit Arb determinants give:

```text
n=0, N=9, ((9^10)): L_n=[-0.0329236641153650833261588690420219786527056477281666156649796573329067 +/- 4.60E-71]; sign=negative
n=1, N=10, ((10^10)): L_n=[-0.03408614168435516582333818798571984202088637651796258403447069775257 +/- 1.42E-69]; sign=negative
n=2, N=11, ((11^10)): L_n=[-0.0234314406176637696109214723258548906569834163387079762522048107971 +/- 5.84E-68]; sign=negative
n=3, N=12, ((12^10)): L_n=[-0.005349707526967545098317802014181744287132566047650903189376754622 +/- 2.67E-67]; sign=negative
n=4, N=13, ((13^10)): L_n=[0.01327339849106758488246975707860360074010336617198100333757491474 +/- 5.32E-66]; sign=positive
```

Thus the four required deep rectangles
`(9^10),(10^10),(11^10),(12^10)` are strictly negative. The
next rectangle `(13^10)` is strictly positive. These are deep
shapes: their smallest part satisfies `N>=10-1`.

For every direct row, the checker independently verifies the raw
Hankel determinant, the Jacobi-Trudi determinant, the Toda identity,
and overlap with the stable condensation coordinate.

## Finite Sign Map

The cancellation-preserving chain through `Q_9` classifies
`1241` endpoint rows:

```text
negative: n=[0, 1, 2, 3]
positive: 4<=n<=1240
inconclusive: 0
```

The finite positive stretch and the fixed-order eventual-tail theorem
do not repair the all-shift statement: one negative required rectangle
is enough to refute it.

## Consequence

```text
s_((N^m))(h)>0 for every m>=10,N>=m-1 is false
```

The endpoint-to-heat theorem and Toda identity remain valid exact
conditional statements, but their all-order positivity antecedent is
false for the actual endpoint sequence. The programme must no longer
present completion of the deep cone as the remaining theorem.

What survives is the separate problem:

```text
seek a weaker Xi/Phi-specific route to Jensen hyperbolicity that does not assume all-shift signed-Hankel positivity
```

This is a route correction, not negative evidence for RH itself.
