# Centered Mangoldt Endpoint-Composed Vaughan Gate

Date: 2026-07-31

Status: exact Vaughan handoff and cancellation-guard artifact; not a proof
of a signed centered-jet lower bound, Abel gap, `Lambda<=0`, PF-infinity,
RH, or a Clay-prize result.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_first_order_centered_mangoldt_endpoint_composed_vaughan_gate.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_mangoldt_endpoint_composed_vaughan_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_mangoldt_endpoint_composed_vaughan_gate.py
```

## Why The Extra Convolution Matters

The physical observable is `sum_(n<=N)log(n)z_n`, not
`sum_(n<=N)Lambda(n)z_n`.  Therefore Vaughan is applied to the Mangoldt
divisor in `log=Lambda*1`, and the final convolution by `1` must remain.

```text
Let lambda_0=Lambda*1_(n<=U), lambda_1=Lambda-lambda_0, mu_0=mu*1_(n<=V), and mu_1=mu-mu_0. Since Lambda=mu*log, log=Lambda*1, and mu*1=delta, Lambda=lambda_0+mu_0*log-lambda_0*mu_0*1+mu_1*lambda_1*1. Therefore log=(lambda_0*1)+(mu_0*log*1-lambda_0*mu_0*1*1)+(mu_1*lambda_1*1*1).
```

Put

```text
Define C_0=lambda_0*1, C_I=mu_0*log*1-lambda_0*mu_0*1*1, and C_II=mu_1*lambda_1*1*1. Then C_0+C_I+C_II=log coefficientwise and Z_Lambda=Z_0+Z_I+Z_II, where Z_j=sum_(n<=N)C_j(n)z_n.
```

This is an exact coefficient identity for every finite `N,U,V`.

## Endpoint Composition

```text
With R_(0;U,V)=E_N+s_*'log(a)W_0-s_*'Z_0, mathfrak B_N=R_(0;U,V)-s_*'Z_I-s_*'Z_II exactly. The endpoint and W_0 term are composed before any real projection or absolute value.

At W_0=0, R_(0;U,V)=E_N-s_*'Z_0. At mathsf_X=0, mathcal C_N=Re[R_(0;U,V)-s_*'Z_I-s_*'Z_II]. No division by W_0, H_a, or a Vaughan component is used.
```

The decomposition has not changed the live scalar or deleted either
exceptional fibre.

## First Stronger Candidate

```text
Put R_0=Re R_(0;U,V), R_I=Re(s_*'Z_I), and R_II=Re(s_*'Z_II). The reverse triangle inequality gives |Re mathfrak B_N|>=|R_0|-|R_I|-|R_II|. Hence |R_0|-|R_I|-|R_II|>A_L+2epsilon_term is sufficient for the Abel gap, but no such uniform inequality is proved.
```

This criterion is rigorously sufficient, but it is deliberately stronger
than the target.  Squaring does not permit the pieces to be separated:

```text
Squaring the exact real scalar retains (R_0-R_I-R_II)^2=R_0^2+R_I^2+R_II^2-2R_0R_I-2R_0R_II+2R_IR_II. The three cross terms are part of the theorem. Separate Type-I/II norms or post-decomposition absolute values do not control them.
```

## Canonical Cubic Split

```text
The canonical finite choice U=V=K_N=floor(N^(1/3)) is an exact cubic-balanced coordinate. It creates no estimate and no sign by itself.

If N+1 is not a perfect cube, K_(N+1)=K_N and each Vaughan component gains only C_j(N+1)z_(N+1). Their sum is log(N+1)z_(N+1).

If N+1=r^3, K changes from r-1 to r. For m<=N put tau_j(m)=C_j^(r)(m)-C_j^(r-1)(m). Then sum_j tau_j(m)=0 coefficientwise, so all historical parameter transfers cancel. The entering coefficients still sum to log(r^3).
```

The exact audit checks `214` transitions:
`209` ordinary and
`5` perfect-cube transitions.  It
also contains `1` simultaneous
square/cube overlap, at `n=64`, and has
`0` coefficient mismatches.

## Exact Cancellation Guards

```text
The exact 1,5,25 chain has W_0=0, zero endpoint, Z_I=log(5), Z_II=-log(5), and Z_Lambda=0. Thus signed Type-I/II cancellation can be complete while separate absolute projections pay log(5)/8.
```

The countermodel uses `s_*'=1/16-i/2`.  The two real projections are
`+log(5)/16` and `-log(5)/16`; the exact signed sum is zero, while separate
absolute values cost `log(5)/8`.

```text
The exact 1,11 prime-edge model has W_0=0 and Z_I=log(11). Choosing the allowed formal endpoint slope g=s_*'log(11) gives mathfrak B_N=0. A proof must use the actual Xi coupling of the endpoint to the carriers.
```

These are route countermodels.  They prove no physical Xi failure.

## Surviving Theorem

```text
The new decomposition identifies a precise missing input: an Xi-specific pointwise signed correlation theorem for the endpoint core and the joint Type-I/II pair on the contact band. Mean-square Vaughan estimates, componentwise upper bounds, and generic amplitude ordering do not imply it.

Retain the cubic-balanced Vaughan split as a proof-facing coordinate and falsification harness. The next candidate must state a joint signed endpoint/Type-I/Type-II inequality for the actual logarithmic carriers and survive W_0=0, mathsf_X=0, H_a=0, q=1, prime edges, cube transfers, the n=64 square/cube overlap, and adjacent real projection.
```

The earlier Mobius Vaughan gates concern a different bridge problem and
do not supply this pointwise signed theorem.

## Pi Provenance

No new `pi` appears in the Vaughan identity.  The only `pi` retained by
the endpoint-complete Xi coordinate remains inherited from the
completed-zeta and Riemann-Siegel normalization.

## Boundary

This gate proves the finite Vaughan convolution identity for the physical
logarithmic moment, endpoint composition, exceptional-fibre formulas,
reverse-triangle sufficient criterion, cubic-root cutoff transport, and
two exact cancellation guards.  It proves no pointwise Xi-specific signed
correlation estimate, signed centered-jet lower bound, Abel-scalar gap,
horizontal successor winding cap, contact exclusion, Q209, cofinal
descendant theorem, `Lambda<=0`, PF-infinity, RH, or prize-level result.
