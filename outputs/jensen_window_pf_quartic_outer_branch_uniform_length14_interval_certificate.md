# Quartic Outer-Branch Uniform Length-14 Interval Certificate

Date: 2026-07-25

Status: rigorous one-contact uniform length-fourteen obstruction;
not a proof of the Xi quartic threshold, `Lambda <= 0`, or RH.

```text
work/rh_compute/results/jensen_window_pf_quartic_outer_branch_uniform_length14_interval_certificate.json
work/rh_compute/results/jensen_window_pf_quartic_outer_branch_uniform_length14_interval_events_endpoint.jsonl
python work/rh_compute/scripts/jensen_window_pf_quartic_outer_branch_uniform_length14_interval_certificate.py --resume
python work/rh_compute/scripts/check_jensen_window_pf_quartic_outer_branch_uniform_length14_interval_certificate.py
```

## Exact Reduction

For the fixed outward contact through `x_5`, the first
scaled-defect inequality is exactly

```text
13*d_6>11*d_5  iff  0<y_6<Y_6,
Y_6=9865549779980282617079853707566678613195505108205132800/26486057908576767522100028185953787483471843406523808729.
```

All later signed order-three/order-four corridors have
`0<y_k<1`. Cancellation of the gap recurrence gives positive
monomials in the `y_k` and inverse powers of earlier `x_k`, so
the rational recurrence extends continuously to the closed box

```text
[0,Y_6] x [0,1]^7.
```

## Rigorous Cover

The append-only event log is a complete binary partition of
that box. On each leaf, 256-bit Arb automatic differentiation
fixes every coordinate having a sign-definite derivative at its
maximizing face and applies a mean-value enclosure to the
remaining coordinates.

```text
events: 74947
splits: 37473
certified leaves: 37474
maximum depth: 23
largest certified upper endpoint: -1.5503268049844423e-10
pending boxes: 0
```

Every upper endpoint is strictly negative. Therefore
`Delta_14<0` throughout the closed proving box.

The exact optimizer-indicated corner independently gives

```text
Delta_14(Y_6,1,...,1)=[-1.251561117811508282044237455837485634774282332162030473736406679569461969442e-6 +/- 3.68e-82]<0.
```

## Scope

This is uniform over every tail from this one fixed outer
contact satisfying the first scaled step. It upgrades two
fixed-tail examples to a complete same-contact obstruction at
length 14. It is not uniform over the complete `(a,p,u)` outer
contact family and does not prove the Xi quartic threshold,
PF-infinity, `Lambda <= 0`, or RH.
