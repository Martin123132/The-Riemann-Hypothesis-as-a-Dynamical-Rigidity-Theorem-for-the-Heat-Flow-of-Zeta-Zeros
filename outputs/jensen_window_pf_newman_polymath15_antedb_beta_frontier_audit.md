# Jensen-Window PF Newman Polymath-15 ANTEDB Beta-Frontier Audit

Date: 2026-07-17

Status: current-source frontier audit and route nonpromotion gate. The
required cancellation improvement remains open. This is not a proof of
`Lambda <= 0` or RH.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_antedb_beta_frontier_audit.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_antedb_beta_frontier_audit.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_antedb_beta_frontier_audit.py
```

## Source Snapshot

The audit uses ANTEDB commit `99668603896af86e6cda90ed6755cf3116aab0ac` dated 2026-07-11 and
the 2025 Tao-Trudgian-Yang paper. It includes their four new exponent
pairs and the two Cushing 2025 pairs recorded by the current blueprint.
This is an as-of-source audit, not a claim about future literature.

## Exact Coordinate

```text
For tau=N^(2+o(1)) and block length M=N^r=tau^(alpha+o(1)), alpha=r/2 and Phi(r)=2*beta(alpha); hence c_req(alpha)=(4*beta(alpha)-2*alpha)/(alpha*(1-alpha))
```

Thus the current pointwise phase input must be optimized directly as a
beta profile; restricting attention to named exponent pairs could miss a
stronger piecewise beta estimate.

## Current Hull

The independently recomputed rational pair hull, including the six
post-2023 pairs and Heath-Brown sequence through `m=100`, gives

```text
line count=108
r_*=125662/155153
active=TY1,TY2
c_*=4911678521/1933561194=2.540223984760009...
omitted-tail trivial ceiling=9802/4851=2.020614306328592...
```

The new pairs improve intermediate radii, but the global maximum remains
the old `TY1/TY2` contact.

## Direct Beta Audit

The documented one-pass ANTEDB pipeline was reproduced with Python
`3.12.13`, pycddlib `2.1.8.post1`, and
`PYTHONHASHSEED=0`. Two seeded runs returned
68 final pieces. At the
critical point, its left, point, and right pieces meet exactly:

```text
alpha_*=62831/155153
beta_*=220633/620612
18/199+(521/796)*alpha
4742/38463+(88225/153852)*alpha
569/2800+(1053/2800)*alpha
c_req(alpha_*)=4911678521/1933561194
```

Intermediate partition counts and raw row hashes depend on set/tie ordering, so they are not promoted. Two PYTHONHASHSEED=0 runs returned 68 final pieces and the same exact rational contact.

A separate beta-only iteration audit then applied the proved ANTEDB
van der Corput transform twelve times:

```text
Twelve successive applications of the proved ANTEDB beta-to-beta van der Corput lemma improved other intervals but left the exact critical contact and global threshold unchanged on every pass
piece counts=[62, 69, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81]
contact beta on every pass=220633/620612
global threshold on every pass=4911678521/1933561194
```

Finite twelve-pass theorem-search audit only; it does not prove stabilization under infinitely many passes or closure under the full beta/exponent-pair transform system
The broader exploratory beta/exponent-pair fixed-point run was terminated
and contributes no mathematical claim here.

## Consequence

```text
The post-2023 exponent pairs and the pinned current ANTEDB one-pass direct-beta envelope improve intermediate radii but do not lower the TY1/TY2 contact; the audited pointwise-bound frontier remains c_*
Lowering c_* still requires a beta bound strictly below 220633/620612 at alpha=62831/155153, compatible improvement across the exposed neighborhood, or cancellation beyond pointwise beta majorants; this audit supplies none of those inputs
```

At `c=2`, the unchanged local phase deficit is
`3133668399/48144906818`. This remains a route benchmark,
not a proof of cancellation, inner Wronskian separation, or RH.
