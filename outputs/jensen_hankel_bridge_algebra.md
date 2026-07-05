# Jensen/Hankel Bridge Algebra Gate

Date: 2026-07-05

Status: exact algebra and countermodel gate. This is not a proof of RH or
`Lambda <= 0`; it records a low-degree identity and a finite obstruction to
promoting finite reshaped-Hankel signs into a Jensen theorem.

## Purpose

The signed-Hankel route needs a theorem of the form:

```text
all-order reshaped-Hankel sign-consistency for A_k(lambda)
  plus Xi/Phi-specific hypotheses
  => all-degree/all-shift Jensen hyperbolicity
```

The current Arb certificate gives finite evidence for the left-hand side. This
gate records what is exact at degree 2 and where finite low-order evidence
already stops being a theorem.

Executable artifacts:

```text
python work/rh_compute/scripts/jensen_hankel_bridge_algebra.py
python work/rh_compute/scripts/check_jensen_hankel_bridge_algebra.py
work/rh_compute/results/jensen_hankel_bridge_algebra.json
```

Current result:

```text
validated Jensen/Hankel bridge algebra gate: degree2 identity and degree3 finite countermodel with 0 issues
```

## Degree 2 Exact Identity

For:

```text
P_{2,n}(x) = A_n + 2 A_{n+1} x + A_{n+2} x^2
```

the discriminant is:

```text
Delta_2 = 4(A_{n+1}^2 - A_n A_{n+2})
```

and:

```text
Delta_2
= -4 det [[A_n, A_{n+1}],
          [A_{n+1}, A_{n+2}]]
```

So the `m = 1` signed-Hankel condition is exactly the degree-2 Jensen
hyperbolicity condition for positive `A_n`.

## Degree 3 Countermodel Gate

For:

```text
P_3(x) = A_0 + 3 A_1 x + 3 A_2 x^2 + A_3 x^3
```

the cubic discriminant is:

```text
81 A_1^2 A_2^2 - 108 A_0 A_2^3 - 108 A_1^3 A_3
- 27 A_0^2 A_3^2 + 162 A_0 A_1 A_2 A_3
```

The exact rational sequence:

```text
A_0..A_4 = 1, 33/40, 429/640, 4719/12800, 4719/35840
```

has positive coefficients and passes the finite reshaped-Hankel sign checks
for `k = 2,3`, `N = 3`:

```text
k=2 columns (0,1): det = -33/3200
k=2 columns (0,2): det = -4719/25600
k=2 columns (1,2): det = -297297/2048000
k=3 columns (0,1,2): det = -281710143/9175040000
```

The expected sign for both `k=2` and `k=3` is negative, so these finite
reshaped signs pass. But the degree-3 Jensen discriminant is:

```text
-2476526481/3276800000
```

which is negative. Thus the degree-3 Jensen polynomial has a nonreal pair and
is not hyperbolic.

## Boundary

Blocked proof step:

```text
finite reshaped-Hankel sign-consistency for low k,N implies Jensen hyperbolicity
```

Correct use:

```text
finite reshaped-Hankel certificates are theorem-search evidence;
the proof still needs an all-order sign-consistency-to-Jensen bridge
with Xi/Phi-specific hypotheses.
```
