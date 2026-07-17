# Jensen-Window PF Newman Polymath-15 Critical Lehmer Margin Gate

Date: 2026-07-17

Status: rigorous small-margin counter-gate and corrected finite
diagnostics. This is not a proof of `Lambda <= 0` or RH.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_critical_lehmer_margin_gate.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_critical_lehmer_margin_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_critical_lehmer_margin_gate.py
```

## Exact Stress Point

At the explicit rational point `x=14010.16343`, Arb evaluates the
exact `H_0=xi((1+i*x)/2)/8` through third order and divides its
Laguerre expression by the exact normalizer amplitude. It certifies

```text
[0.022157698462940963514405593163006144166470992886153587019597917503836 +/- 3.42e-70]
[0.00045007184792928814038031537868985932326046414288328445442403532767930 +/- 5.92e-72]
0<C[H_0/A_0]/L^2<1/1000<3/40
```

Thus the `3/40*L^2` margin proved on `tL>=25` cannot simply be
continued into the critical layer.

## Strict Monotonicity Counter-Gate

The same Arb jet certifies

```text
[-3.24186746975685391605462931182063307428944779983857697073123186e-6 +/- 5.42e-69]
[-0.0001463088540165360168258875503621433728403782920777681331619380172 +/- 5.15e-68]
M_0(x)=H_0(x)H_0'''(x)-H_0'(x)H_0''(x)=-L_0'(x)<0
```

The theta tail makes `|u|^3*exp(u^2/5)*Phi(u)` integrable, so
dominated convergence makes the heat-flow jet through third order
continuous in `t` on `[0,1/5]`. The strict negative value at `t=0` therefore
persists for all sufficiently small positive `t`. This rigorously
retires the proposed global condition `M_t(x)>0` for every `x>0` and
every `0<t<=1/5`. It does not challenge `L_t(x)>0` itself.

## Why The Margin Shrinks

For a close pair with local factorization
`F(m+y)=(y^2-a^2)G(m+y)`, exact differentiation gives

```text
If F(m+y)=(y^2-a^2)G(m+y), then L[F](m)=2*a^2*G(m)^2+a^4*(G'(m)^2-G(m)G''(m))
If F(m+y)=(y^2-a^2)G(m+y), then M[F](m)=-4*a^2*G(m)*G'(m)+a^4*(G(m)*G'''(m)-G'(m)*G''(m))
```

The leading curvature is quadratic in the half-gap, while the stronger
monotonicity sign also depends on the slope of the regular factor and has
no fixed close-pair sign. Under the Newman
heat sign, a double zero is nevertheless transversal:

```text
For partial_t H=-partial_x^2 H, partial_t L=H_xx^2-2H_x H_xxx+H H_xxxx; at H=H_x=0 this equals H_xx^2
H_tau(y)=y^2-2*tau has L[H_tau](y)=2*y^2+4*tau
```

## Corrected Main

Independent 60- and 80-digit runs of the corrected Polymath-15 main
give:

| c=tL | N | C[J]/L^2 | J | J'/L | J^2+(J'/L)^2 |
|---:|---:|---:|---:|---:|---:|
| 0 | 33 | 0.000444842135 | -0.00392103038 | -0.000213787658 | 1.54201844e-5 |
| 0.01 | 33 | 0.0013460776 | -0.0118525497 | -0.000658322054 | 0.000140916323 |
| 0.1 | 33 | 0.00941328768 | -0.0821325233 | -0.00445932838 | 0.006765637 |
| 0.5 | 33 | 0.0441316203 | -0.371500915 | -0.0172154693 | 0.138309302 |

The numerical rows are diagnostics. Their role is to show the shape
of the viable target: a fixed curvature floor is too blunt near close
pairs, while the two-component first-jet norm directly measures
distance from a double zero and requires only a `C1` error transfer.

No finite or point computation here is promoted to a uniform sign,
absence of positive-time collisions, `Lambda <= 0`, or RH.
