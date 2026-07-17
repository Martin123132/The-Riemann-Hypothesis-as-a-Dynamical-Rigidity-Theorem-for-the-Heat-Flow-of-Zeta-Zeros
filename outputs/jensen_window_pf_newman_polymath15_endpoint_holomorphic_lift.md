# Jensen-Window PF Newman Polymath-15 Endpoint Holomorphic Lift

Date: 2026-07-17

Status: exact holomorphic-lift reduction for the refined endpoint
correction. This is not a proof of `Lambda <= 0` or RH.
It is not a `C1` remainder certificate.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_endpoint_holomorphic_lift.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_endpoint_holomorphic_lift.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_endpoint_holomorphic_lift.py
```

The endpoint term is imported from the A+B-C approximation in
[Polymath 15](https://arxiv.org/abs/1904.12438). Its displayed off-axis real-part formula
is not itself holomorphic, so it cannot be fed directly into Cauchy's
estimate. The following lift repairs that point.

## Exact Phase Cancellation

Direct substitution of `M_0` and `U` gives

```text
M_0(i*T)*U(T)*exp(pi*i/8)=-(sqrt(pi)/8)*exp(-pi*T/4)*(T^(3/2)+i*T^(1/2))
F_N(T)=exp(-pi*T/4)S_N(T), S_N(T)=-(sqrt(pi)/8)*(T^(3/2)+i*T^(1/2))*C_0(p_N(T))
```

The apparent high-frequency phase is gone exactly. Independent
80-digit checks at `T=100,1000,10000` have maximum relative
difference `1.041882638615332398161565e-77`.

## Holomorphic Lift

For fixed `N`, set

```text
T(z)=z/2+pi*t/8, a(z)=sqrt(T(z)/(2*pi)), p_N(z)=1-2*(a(z)-N)
F_N^#(z)=conj(F_N(conj(z)))
C_hat_(N,t)(z)=(-1)^N*exp(t*pi^2/64)*(F_N(z)+F_N^#(z))
```

Then

```text
C_hat_(N,t)(x)=C_(N,t)(x) for every real x in the fixed-N cell
```

so the lift changes nothing in the corrected real-axis main.

## Defect Equation

The exponential derivative cancels:

```text
(d/dz+pi/8)F_N(z)=(1/2)*exp(-pi*T(z)/4)*S_N'(T(z))
```

Comparing the lift with the paper's off-axis extension yields

```text
C_paper(x+iy)=exp(-pi*i*y/8)*C_(N,t)(x)
D=C_hat-C_paper satisfies D(x,0)=0 and partial_y D+i*pi*D/8=i*(partial_z+pi/8)C_hat
|D(x,y)|<=|y|*exp(pi*|y|/8)*sup_|v|<=|y| |(partial_z+pi/8)C_hat(x+iv)|
```

Thus the discrepancy is controlled by the slow derivative, not by
the full endpoint block. Since `p_N'(T)=O(T^-1/2)`, bounded `C_0`
and `C_0'` give

```text
If C_0 and C_0' are bounded on the required p-strip, then |S_N'(T)|=O(T), one saddle factor below |S_N(T)|=O(T^(3/2))
After division by A_t(x), the radius-1/L lift defect is O(exp(-3L/4)/L) uniformly for bounded c=tL
```

## Remaining Local Bound

A 2,001-point high-precision real grid reports
`max |C_0'|=0.6908594309245895837673515` at
`p=-1.0`. This is only a diagnostic.
The next obligation is

```text
Certify explicit complex-strip bounds for C_0 and C_0' and combine them with the published e_A+e_B+e_C constants
```

Closing that finite-dimensional strip estimate would make the refined
endpoint correction compatible with a rigorous `C1` Cauchy transfer.
It would still leave the RH-level corrected first-jet lower bound open.
