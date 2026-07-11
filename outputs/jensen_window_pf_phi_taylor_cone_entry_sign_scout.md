# Jensen-Window PF Phi Taylor Cone-Entry Sign Scout

Date: 2026-07-06

Status: finite Taylor sign certificate. This is not a proof of zeta cone entry,
monotone contractions, Jensen-window PF-infinity, RH, or the Newman-direction
goal.

Artifact kind: `jensen_window_pf_phi_taylor_cone_entry_sign_scout`.

Proof boundary: this artifact certifies only the local Taylor-combination
signs used by the fixed-k large-negative-lambda asymptotic route. It does not
control the asymptotic remainder uniformly in `k`, does not prove the full
infinite or collared finite ratio cone, does not prove `jwpf_06`, and does not
settle `Lambda <= 0`.

Machine-readable certificate:

```text
work/rh_compute/results/jensen_window_pf_phi_taylor_cone_entry_sign_scout.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_phi_taylor_cone_entry_sign_scout.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_phi_taylor_cone_entry_sign_scout.py
```

Current result:

```text
validated Jensen-window PF Phi Taylor cone-entry sign scout: 4 coefficient balls, 2 certified signs, 0 ready-to-apply rows, 0 issues
```

## Role

The cone-entry asymptotic target isolated the fixed-k sign requirements:

```text
Phi(u)=c0+c2*u^2+c4*u^4+c6*u^6+...
a=c2/c0
b=c4/c0
c=c6/c0
2*b-a^2 < 0
2*(a^3-3*a*b+3*c) > 0
```

This certificate supplies those local signs with an explicit finite sum plus
tail bound. It does not supply the missing uniform-in-`k` or collared finite
cone coverage.

Source target:

```text
outputs/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.md
```

## Coefficient Enclosures

```text
c0: [0.4466969004671234440869846670547091132204 +/- 2.44e-41] (positive)
c2: [-16.73050077470325686679012435077297959908 +/- 2.45e-39] (negative)
c4: [270.7209444851414291692625428649721885115 +/- 4.48e-38] (positive)
c6: [-2254.921151385262962904875830379485325230 +/- 3.67e-37] (negative)
```

Parameters:

```text
precision_bits = 256
tail_cutoff_n = 12
tail_safety_factor = 4
```

Tail method:

```text
For each polynomial degree m, bound sum_{n>N} n^(2*m)*exp(-pi*n^2) by first/(1-r), where first=(N+1)^(2*m)*exp(-pi*(N+1)^2) and r=((N+2)/(N+1))^(2*m)*exp(-pi*(2*N+3)).
```

## Certified Sign Combinations

```text
ptces_01_upper_wall_sign:
  2*b-a^2=(2*c4*c0-c2^2)/c0^2
  enclosure: [-190.6865836820046823899449624711272165733 +/- 5.91e-39]
  certified sign: negative
  boundary: Local Taylor sign only; it does not give a uniform-in-k cone-entry theorem or a controlled asymptotic remainder.

ptces_02_monotone_wall_sign:
  2*(a^3-3*a*b+3*c)
  enclosure: [825.9976249272062741326759643225349914950 +/- 2.19e-38]
  certified sign: positive
  boundary: Local Taylor sign only; it does not give a uniform-in-k cone-entry theorem or a controlled asymptotic remainder.
```

## Consequence

The local Phi Taylor signs are no longer the immediate obstruction in the
fixed-k cone-entry route. The remaining hard work is to convert the fixed-k
asymptotic indication into a theorem with controlled remainders and either
uniform-in-`k` coverage or a rigorous finite-prefix collar plus tail argument.
