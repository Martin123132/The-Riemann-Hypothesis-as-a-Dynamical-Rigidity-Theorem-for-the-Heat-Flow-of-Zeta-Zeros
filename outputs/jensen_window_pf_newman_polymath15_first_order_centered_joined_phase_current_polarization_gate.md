# Joined Phase-Current Polarization Gate

Date: 2026-07-30

Status: exact joined-current identity and nonpromotion gate. This is not
a proof of a joined Xi Abel-gap theorem, `Lambda<=0`, PF-infinity, or RH.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_first_order_centered_joined_phase_current_polarization_gate.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_joined_phase_current_polarization_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_joined_phase_current_polarization_gate.py
```

## Exact Polarization

For a finite endpoint-complete decomposition `W=sum_a w_a`, define

```text
K(W)=Im((D_j W) conjugate(W)).
```

Direct expansion gives

```text
K(W)
 =sum_a Im((D_j w_a)conjugate(w_a))
  +sum_(a<b) Im((D_j w_a)conjugate(w_b)
                +(D_j w_b)conjugate(w_a)).
```

The checker reconstructs this identity over exact rational complex
arithmetic. No quotient or nonvanishing assumption is used.

## Contact Identity

Write `W=mathsf_X+i mathsf_Y` and retain the certified orientation
defect

```text
E_j=D_j mathsf_X-x mathcal_C_N.
```

Then

```text
K=mathsf_X D_j mathsf_Y
  -x mathsf_Y mathcal_C_N
  -mathsf_Y E_j.
```

At `mathsf_X=0`,

```text
K=-mathsf_Y(x mathcal_C_N+E_j).
```

If `mathsf_Y=0`, then `W=0` and `K=0`; phase current cannot replace
the division-free Abel scalar on that fibre.

## Prime-Power Insertion

For one physical chain `W_m=eta B_m Q_m`, let

```text
tau_m=-x log(p) b>0,
e_m=x(beta_m'-b log(m)),
eta_L=i Omega_eta eta.
```

The exact fixed-ray formula gives

```text
K_m
 =|B_m|^2[
    tau_m J_ray,m
    +(Omega_eta+e_m-tau_m)|Q_m|^2].
```

Thus the proved positive currents for `p=2,M>=8`,
`p=3,M>=4`, and `p>=5,M>=2` enter the joined current only through
the displayed diagonal term. External rates, the recurrent endpoint,
short chains, singletons, and every cross current remain.

## Exact Guards

For

```text
W_a(theta)=exp(i theta)+a exp(3i theta),
```

the joined current at `theta=pi/2` is exactly

```text
K=(1-a)(1-3a).
```

Both self currents are positive: `1` and `3a^2`. Nevertheless:

```text
a=1/3:
  Re(W)=0, Re(W')=0,
  cos(theta)+(1/3)cos(3theta)=(4/3)cos(theta)^3,
  K=0/1;

a=1/2:
  K=-1/4<0.
```

The cubic row is a real-projection contact produced by exact cross-current
cancellation. The negative row shows that individual starlikeness is not
closed under addition. These are generic route guards, not Xi
counterexamples or phase-attainment claims.

## Route Decision

The long-chain theorems remain valid and useful diagonal input, but they
do not by themselves lower-bound the endpoint-complete Abel scalar or
the joined phase current. A viable continuation must either:

1. prove an Xi-specific bound for the endpoint, external-rate, short,
   singleton, and pair cross-current aggregate; or
2. prove the division-free Abel-scalar gap directly, retaining `W_0=0`.

## Pi Provenance

The `pi/2` in the guard is the ordinary quarter-period of the complex
exponential. Any `pi` in the Xi saddle scale remains inherited from
completed-zeta normalization; no polygon or prime block defines it.

## Boundary

This gate proves the exact joined-current polarization, contact identity, prime-power diagonal insertion, and generic cubic-contact and negative-current guards. It proves no Xi-specific cross-current bound, endpoint or singleton absorption, Abel-scalar gap, horizontal successor winding cap, contact exclusion, Q209, cofinal descendant theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion.
