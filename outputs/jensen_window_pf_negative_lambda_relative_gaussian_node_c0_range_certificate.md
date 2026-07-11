# Jensen-Window PF Negative-Lambda Relative-Gaussian Node-C0 Range Certificate

Date: 2026-07-07

Status: node-range and Phi0 lower certificate. This is not a proof
of a finite-grid interval certificate, a uniform residual estimate,
scaled-curvature monotonicity, cone entry, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_node_c0_range_certificate`.

Proof boundary: this artifact certifies only the two side conditions
needed by the padded Phi-tail scout: `x<=1` for the recorded
Laguerre-node grid and `Phi(0)>=0.44`. It does not certify
individual nodes, weights, quadrature remainders, coefficient
propagation, rounding, or a grid-to-collar bridge.

Machine-readable artifact:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_node_c0_range_certificate.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_node_c0_range_certificate.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_node_c0_range_certificate.py
```

Current result:

```text
validated Jensen-window PF negative-lambda relative-Gaussian node-c0 range certificate: 5 rows, 0 issues, 16 Laguerre bound rows, 2 certified side conditions, 0 ready-to-apply rows
```

## Node Range

```text
min T: 1156
max quadrature order: 192
max index: 24
worst Laguerre node upper bound: 809
worst x^2 upper bound: 809/1156
worst x^2 upper bound decimal: 6.998269896193771711E-01
worst x^2 slack to one: 347/1156
node range x<=1 certified: True
```

The exact row bound uses the Laguerre Jacobi matrix, Gershgorin,
and AM-GM. For the worst recorded case, `N=192` and
`alpha=47/2`, every node is bounded by `809`, while `T>=1156`.

## Phi0 Lower Bound

```text
certified c0 lower: 0.44
n=1 Phi(0) term lower: 0.44572697131799474430335802786658539152643783769655419085980516650442890855659827
n=1 margin over c0 lower: 0.0057269713179947443033580278665853915264378376965541908598051534749906705348328164
c0 lower certified by n=1 term: True
```

Because

```text
2*pi^2*n^4 - 3*pi*n^2 = pi*n^2*(2*pi*n^2-3)>0 for n>=1, using pi>3
```

the full `Phi(0)` is at least the certified `n=1` contribution.

Integration:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout.json
outputs/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_bound_scout.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_bound_scout.json
outputs/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.json
outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md
outputs/signed_hankel_jensen_dependency_graph.md
```

Summary:

The finite relative-Gaussian grid satisfies the two concrete side conditions used by the padded Phi-tail scout: Gershgorin/AM-GM gives every recorded Laguerre node <=809<T_min=1156, hence x<=1, and Arb certifies the n=1 term in Phi(0) already exceeds 0.44. This supports the tail scout but leaves weights, quadrature error, coefficient propagation, rounding, and the grid-to-collar bridge open.
