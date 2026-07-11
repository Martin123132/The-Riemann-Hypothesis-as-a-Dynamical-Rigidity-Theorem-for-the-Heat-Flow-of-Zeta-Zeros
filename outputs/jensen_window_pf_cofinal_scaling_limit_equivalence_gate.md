# Jensen-Window PF Cofinal Scaling-Limit Equivalence Gate

Date: 2026-07-10

Status: exact scaling-limit and proof-route guard. This is not a proof
of Laguerre-Polya membership, PF-infinity, RH, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_cofinal_scaling_limit_equivalence_gate.json
python work/rh_compute/scripts/jensen_window_pf_cofinal_scaling_limit_equivalence_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_cofinal_scaling_limit_equivalence_gate.py
```

Current result:

```text
validated Jensen-window PF cofinal scaling-limit equivalence gate: 10 rows, 0 issues, 2 exact scaling identities, 2 analytic limit steps, 1 cofinal-to-LP theorem, 1 LP-to-all-degrees theorem, 1 fixed-shift equivalence, 2 non-promotion guards, 1 open independent-route handoff
```

## Scaling Limit

For a fixed shift `n`,

```text
P_(D,n)(z/D)=sum_(j=0)^D ((D)_j/D^j)*A_(n+j)*z^j/j!
(D)_j/D^j=product_(m=0)^(j-1)(1-m/D)->1 for fixed j
```

Since `F_n` is entire, the convergence is locally uniform:

```text
P_(D,n)(z/D)->F_n(z) locally uniformly because F_n is entire
```

If an unbounded subsequence of the scaled polynomials is hyperbolic,
the Laguerre-Polya closure theorem places the nonzero limit `F_n` in
the Laguerre-Polya class. Conversely, Jensen's theorem gives every
degree once `F_n` is Laguerre-Polya. Hence

```text
cofinal degree hyperbolicity at fixed n <=> F_n is Laguerre-Polya
```

## Proof-Safety Consequence

At shift zero, proving the cofinal terminal sequence is already an
endpoint theorem, not a weaker corollary of familiar fixed-degree
asymptotics. The 1,050 finite Sturm rows through degree 12 do not enter
this limit because their degrees are bounded.

The cofinal formulation can still be useful if heat-flow or kernel
structure proves it independently. It cannot be assumed from
Laguerre-Polya membership or from the desired RH endpoint without
becoming circular.
