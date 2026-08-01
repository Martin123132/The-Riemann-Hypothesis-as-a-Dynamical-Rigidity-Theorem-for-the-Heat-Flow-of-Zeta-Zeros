# Newman Single-Carrier Adiabatic Benchmark

Date: 2026-07-25

Status: exact conditional thought-experiment theorem and rejection guard; not a proof of `Lambda<=0` or RH.

## Hypothetical Carrier

Let

```text
E=A*exp(i*theta), A>0,
r=E_x/E=u+i*v,
J=2*Re(E),
V_ell=(J,J_x/ell), ell>0.
```

Then

```text
V_ell=2*A [[1,0],[u/ell,-v/ell]]
              (cos(theta),sin(theta))^T.          (1)
```

The matrix in (1) has determinant `-v/ell` and squared
Frobenius norm `1+(u^2+v^2)/ell^2`. Therefore

```text
||V_ell||_2
 >=2*A*|v|/sqrt(ell^2+u^2+v^2).                 (2)
```

This is the clean hypothetical geometry: a nonzero phase
speed supplies a phase-independent first-jet floor.

## Heat Transport

If the carrier itself obeys `E_t=-E_xx`, then

```text
E_xx/E=r^2+r_x,
E_xxx/E=r^3+3*r*r_x+r_xx.                       (3)
```

Assume fixed constants `c>0`, `C1,C2,C3>=0` satisfy

```text
|v|>=c*ell,
|r|<=C1*ell,
|r_x|<=C2*ell^2,
|r_xx|<=C3*ell^3.                               (4)
```

Put

```text
B2=C1^2+C2,
B3=C1^3+3*C1*C2+C3,
K=sqrt(1+C1^2)*sqrt(B2^2+B3^2)/c.
```

Equations (2)-(4) give the exact conditional bound

```text
||partial_t V_ell||_2/||V_ell||_2<=K*ell^2.     (5)
```

Gronwall then transports a nonzero old jet multiplicatively:

```text
||V_ell(t,x)||_2
 >=exp(-K*ell(x)^2*|t-t_j|)*||V_ell(t_j,x)||_2. (6)
```

For the cofinal shells,

```text
delta_j=1/(5*j*(j+1)),
R_j=j+38,
ell_j=max(1,log(R_j/(4*pi))),
delta_j*ell_j^2 -> 0.                           (7)
```

Thus the ideal single-carrier thought experiment really does
make the successor increasingly adiabatic.

## Why It Is Not Xi Yet

The exact four-carrier model

```text
E_toy=exp(-3ix)-exp(-4ix)+i*exp(-ix)
      -(i/2)*exp(-2ix)
```

has positive amplitudes and distinct ordered phase speeds, but

```text
E_toy(0)=i/2,
E_toy'(0)=i,
Re(E_toy)(0)=Re(E_toy'(0))=0,
Re(E_toy)''(0)=7,
Im(E_toy'(0)/E_toy(0))=0.                       (8)
```

So summing individually good carriers can destroy the phase-speed
floor exactly at a real double zero. Generic positive amplitudes,
ordered speeds, and omitted-tail domination cannot establish (4).

## Xi Handoff

For the corrected Riemann-Siegel main `E=X+iY` and
`E'=U+iV`, the live theorem remains the crossing-restricted
arithmetic estimate

```text
X^2+(U/L)^2
 >8000000*exp(-3*L/2)                            (9)
```

through the residual `tL->0` layer, or an equivalent
Xi-specific restriction that excludes the negative Wronskian
direction on the crossing set. Equation (7) explains why such a
theorem would make the adiabatic successor work; it does not prove
that theorem.

## Proof Boundary

The carrier algebra, singular-value floor, logarithmic-jet
condition number, cofinal `delta_j*ell_j^2` limit, Gronwall
transport, and four-carrier countermodel are exact. The uniform
carrier hypotheses are not proved for Xi. The crossing small-ball
estimate, all-`j` successor theorem, `Lambda<=0`, RH, PF-infinity,
and the Clay prize remain open.

Machine-audited files:

```text
work/rh_compute/results/jensen_window_pf_newman_single_carrier_adiabatic_benchmark.json
work/rh_compute/scripts/jensen_window_pf_newman_single_carrier_adiabatic_benchmark.py
work/rh_compute/scripts/check_jensen_window_pf_newman_single_carrier_adiabatic_benchmark.py
```
