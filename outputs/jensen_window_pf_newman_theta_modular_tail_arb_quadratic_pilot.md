# Newman Theta Modular-Tail Arb Quadratic Pilot

Date: 2026-07-24

Status: exact weighted reduction plus one rigorous finite compact
Arb integral. This is not an exact infinite-tail bound and not a
proof of `Lambda<=0`, RH, or a Clay-prize result.

```text
work/rh_compute/results/jensen_window_pf_newman_theta_modular_tail_arb_quadratic_pilot.json
python work/rh_compute/scripts/jensen_window_pf_newman_theta_modular_tail_arb_quadratic_pilot.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_theta_modular_tail_arb_quadratic_pilot.py
```

## Exact Weighted Reduction

For a real function `f` and nonnegative weight `W`,

```text
int_I W|f| <= (int_I W)^(1/2)(int_I W*f^2)^(1/2)
```

This is exactly Cauchy-Schwarz. Applied term by term to the
positive heat-envelope weights
`W_b=exp(Tu^2)binom(k,b)P_b(T,u)`, it replaces each nonsmooth
absolute-value integral by a weight mass and an analytic
quadratic integral.

## Rigorous Finite Compact Pilot

At 96-bit Arb precision, the certified quantity is

```text
integral_0^(11/5) exp(u^2/5)*[D^9 sum_(n=7)^16 b_n(u)]^2 du
in [1.4163309784009929435546867526867157e-20 +/- 3.65e-30]
```

The enclosure has 31 relative
accuracy bits. The integrand is entire and the Arb callback
therefore has no branch-cut qualification.

## Non-Promotion Guard

This is only the finite arithmetic block `n=7..16` on
`u in [0,11/5]`. It omits `n>=17`, `u>11/5`, the remaining
heat-polynomial weights and derivative orders, and all retained
first-jet lower bounds.

## Next Certificate

Build the finite `(b,k,N,t)` quadratic matrix, attach analytic
bounds for the omitted arithmetic and outer-u tails, and compare
the resulting `d0,d1` balls with rigorous retained `J,J'` lower
balls.

No strict Laguerre conclusion, `Lambda<=0`, RH, or Clay-prize
conclusion is claimed.
