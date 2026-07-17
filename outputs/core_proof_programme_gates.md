# Core Proof-Programme Gate Runner

Date: 2026-07-17

Status: reproducibility ledger. This is not a proof of PF-infinity,
Laguerre-Polya membership, RH, or `Lambda <= 0`.

## Purpose

The corpus now has many individual finite-evidence manifests and
countermodel guards. The core runner gives a single command that checks the
main promoted finite certificates, dependency hygiene graph, and the
non-promotion gates used in the referee-facing notes.

Script:

```text
work/rh_compute/scripts/check_core_proof_programme_gates.py
```

Command:

```text
python work/rh_compute/scripts/check_core_proof_programme_gates.py
python work/rh_compute/scripts/check_core_proof_programme_gates.py --skip-slow
```

Current result:

The runner has `349` gates: `333` non-slow and `16` slow. The latest full
single-process umbrella execution before the order-nine integration completed
`301/301` with no failures in `1348.4` seconds. The fifteen order-nine and
all-order gates have also passed focused validation, including the two slow
order-nine cache/transfer checks. They complete the order-nine curvature and
endpoint chain, signed contiguous/arbitrary-column closure through fixed order
nine, and the all-order endpoint-to-heat reduction. The deep-Schur coordinate,
Toda/boundary gate, and rigorous order-ten counterexample now reject the
former all-order endpoint antecedent and raised the ledger to 295 claims. The
current non-slow umbrella run completed `304/304` with no failures in `1104.5`
seconds, including the newly registered counterexample gate. The subsequently
registered Polymath-15 and Edrei gates, the expanded 304-claim ledger, and the
new Xi strict-monotonicity rejection have passed focused validation. Sixteen
order-ten cache, curvature, delayed-heat, and lambda-zero completion gates are
now registered. The complete `328/328` non-slow umbrella run passed with no
failures in `969.3` seconds. Five subsequently registered order-eleven gates
have passed focused validation, and the claim ledger now contains 324 claims.

```text
OK core gate: countermodel proof-safety gates
OK core gate: result-language boundary scan
OK core gate: output artifact status manifest
OK core gate: output reference integrity
OK core gate: proof-claim ledger
OK core gate: signed-Hankel/Jensen dependency graph
OK core gate: signed-Hankel finite certificates
OK core gate: Toeplitz/PF finite certificates
OK core gate: Toeplitz/Jacobi-Trudi reindexing
OK core gate: Hankel sign-consistency reduction point audits
OK core gate: Hankel sign-consistency reduction finite certificates
OK core gate: shifted Hankel staircase finite certificates
OK core gate: Jensen/Hankel bridge algebra gate
OK core gate: Jensen-window PF obligation algebra gate
OK core gate: Arb Jensen-window PF obligation finite diagnostics
OK core gate: Arb Jensen-window Sturm hyperbolicity finite diagnostics
OK core gate: Arb Jensen-window Sturm degree-5 hyperbolicity finite diagnostics
OK core gate: Arb Jensen-window Sturm degree-6 through degree-12 hyperbolicity finite diagnostics
OK core gate: finite Sturm-to-PF Jensen-window consequences
OK core gate: signed-Hankel/Jensen bridge target specification
OK core gate: Jensen-window PF bridge target
OK core gate: Jensen-window PF bridge obligation ledger
OK core gate: Jensen-window PF theorem machinery fit matrix
OK core gate: Jensen-window PF sign-regular transfer gap matrix
OK core gate: Jensen-window PF factorial multiplier split audit
OK core gate: Jensen-window PF reciprocal-gamma mixture sign gate
OK core gate: Jensen-window PF reciprocal-defect compound order-three gate
OK core gate: Jensen-window PF negative-lambda -100 compound order-three entry certificate
OK core gate: Jensen-window PF compound order-three forward-invariance certificate
OK core gate: Jensen-window PF order-three noncontiguous secant-transfer lemma
OK core gate: Jensen-window PF compound order-four condensation gate
OK core gate: Jensen-window PF compound order-four first-summand curvature bridge
OK core gate: Jensen-window PF compound order-four localized-curvature compact certificate
OK core gate: Jensen-window PF compound order-four Gaussian cumulant ray target
OK core gate: Jensen-window PF compound order-four formal cumulant corridor certificate
OK core gate: Jensen-window PF compound order-four formal cumulant asymptotic ray certificate
OK core gate: Jensen-window PF compound order-four formal cumulant next-parity certificate
OK core gate: Jensen-window PF compound order-four next-parity finite certificate
OK core gate: Jensen-window PF compound order-four next-parity asymptotic ray certificate
OK core gate: Jensen-window PF compound order-four exact cumulant remainder budget
OK core gate: Jensen-window PF compound order-four formal cumulant second-next parity certificate
OK core gate: Jensen-window PF compound order-four second-next finite certificate
OK core gate: Jensen-window PF compound order-four second-next asymptotic ray certificate
OK core gate: Jensen-window PF compound order-four exact cumulant complex-disk contract
OK core gate: Jensen-window PF compound order-four exact cumulant formal tails
OK core gate: Jensen-window PF compound order-four exact cumulant exact tails
OK core gate: Jensen-window PF compound order-four exact cumulant corridor theorem
OK core gate: Jensen-window PF compound order-four finite partition extension
OK core gate: Jensen-window PF compound order-four exact cumulant central residual
OK core gate: Jensen-window PF compound order-four finite exact-corridor curvature theorem
OK core gate: Jensen-window PF compound order-four exact-corridor curvature ray theorem
OK core gate: Jensen-window PF compound order-four lambda=-100 entry theorem
OK core gate: Jensen-window PF compound order-four forward-flow reduction
OK core gate: Arb Xi lambda-zero order-four prefix certificate
OK core gate: Jensen-window PF compound order-four lambda-zero eventual positivity
OK core gate: Jensen-window PF lambda-zero first-summand dominance transfer
OK core gate: Jensen-window PF compound order-four eventual-tail forward reduction
OK core gate: Jensen-window PF compound order-four uniform-heat tail reduction
OK core gate: Jensen-window PF uniform superpolynomial first-summand dominance
OK core gate: Jensen-window PF uniform first-summand heat-tilt asymptotic theorem
OK core gate: Jensen-window PF compound order-four uniform-heat forward invariance
OK core gate: Jensen-window PF order-four noncontiguous total-positivity transfer
OK core gate: Jensen-window PF compound order-five uniform-tail and flow reduction
OK core gate: Jensen-window PF compound order-five lambda=-100 prefix
OK core gate: Jensen-window PF compound order-five lambda=-100 tail curvature reduction
OK core gate: Jensen-window PF compound order-five first-summand curvature bridge
OK core gate: Jensen-window PF compound order-five nested-curvature compact certificate
OK core gate: Jensen-window PF compound order-five nested-curvature finite ray
OK core gate: Jensen-window PF compound order-five nested-curvature asymptotic ray
OK core gate: Jensen-window PF compound order-five lambda=-100 entry
OK core gate: Jensen-window PF compound order-five uniform-heat forward invariance
OK core gate: Jensen-window PF compound order-six uniform-tail and flow reduction
OK core gate: Jensen-window PF compound order-six lambda=-100 prefix
OK core gate: Jensen-window PF compound order-six lambda=-100 tail curvature reduction
OK core gate: Jensen-window PF inverse-seventh-power first-summand dominance
OK core gate: Jensen-window PF compound order-six first/full curvature bridge
OK core gate: Jensen-window PF compound order-six high-cumulant corridor
OK core gate: Jensen-window PF compound order-six nested-curvature compact certificate
OK core gate: Jensen-window PF compound order-six nested-curvature finite ray
OK core gate: Jensen-window PF compound order-six nested-curvature asymptotic ray
OK core gate: Jensen-window PF compound order-six lambda=-100 entry
OK core gate: Jensen-window PF compound order-six uniform-heat forward invariance
OK core gate: Jensen-window PF graded-kernel all-order Vandermonde lemma
OK core gate: Jensen-window PF compound order-seven uniform-tail and flow reduction
OK core gate: Jensen-window PF compound order-seven lambda=-100 prefix
OK core gate: Jensen-window PF compound order-seven lambda=-100 tail curvature reduction
OK core gate: Jensen-window PF inverse-eighth-power rebalanced first-summand dominance
OK core gate: Jensen-window PF compound order-seven first/full curvature bridge
OK core gate: Jensen-window PF compound order-seven shifted-jet compact bridge
OK core gate: Jensen-window PF compound order-seven nested-curvature compact certificate
OK core gate: Jensen-window PF compound order-seven high-cumulant coarse corridor
OK core gate: Jensen-window PF compound order-seven nested-curvature finite-ray certificate
OK core gate: Jensen-window PF compound order-seven nested-curvature asymptotic-ray certificate
OK core gate: Jensen-window PF compound order-seven lambda=-100 entry
OK core gate: Jensen-window PF compound order-seven uniform-heat forward invariance
OK core gate: Jensen-window PF compound order-eight uniform-tail and flow reduction
OK core gate: Jensen-window PF compound order-eight lambda=-100 prefix
OK core gate: Jensen-window PF compound order-eight lambda=-100 tail-curvature reduction
OK core gate: Jensen-window PF inverse-ninth-power first-summand dominance
OK core gate: Jensen-window PF compound order-eight first/full curvature bridge
OK core gate: Jensen-window PF compound order-eight high-cumulant coarse corridor
OK core gate: Jensen-window PF compound order-eight shifted-jet certificate
OK core gate: Jensen-window PF compound order-eight nested-curvature compact certificate
OK core gate: Jensen-window PF compound order-eight nested-curvature finite-ray certificate
OK core gate: Jensen-window PF compound order-eight nested-curvature asymptotic-ray certificate
OK core gate: Jensen-window PF compound order-eight lambda=-100 entry
OK core gate: Jensen-window PF compound order-eight uniform-heat forward invariance
OK core gate: Jensen-window PF compound order-nine uniform-tail and flow reduction
OK core gate: Jensen-window PF compound order-nine lambda=-100 prefix
OK core gate: Jensen-window PF compound order-nine lambda=-100 tail-curvature reduction
OK core gate: Jensen-window PF compound order-nine first/full curvature bridge
OK core gate: Jensen-window PF compound order-nine high-cumulant coarse corridor
OK core gate: Jensen-window PF compound order-nine exact-point H0-H8 cache
OK core gate: Jensen-window PF compound order-nine localized lower bridge
OK core gate: Jensen-window PF compound order-nine nested-curvature compact certificate
OK core gate: Jensen-window PF compound order-nine nested-curvature finite-ray certificate
OK core gate: Jensen-window PF compound order-nine nested-curvature asymptotic-ray certificate
OK core gate: Jensen-window PF compound order-nine global first-summand curvature
OK core gate: Jensen-window PF compound order-nine finite endpoint splice
OK core gate: Jensen-window PF compound order-nine lambda=-100 entry
OK core gate: Jensen-window PF compound order-nine uniform-heat forward invariance
OK core gate: Jensen-window PF all-order endpoint-to-heat reduction
OK core gate: Jensen-window PF endpoint deep-Schur coordinate
OK core gate: Jensen-window PF deep-Schur Toda/boundary gate
OK core gate: Jensen-window PF endpoint order-ten counterexample
OK core gate: Jensen-window PF order-moment transport fit gate
OK core gate: Jensen-window PF structural ansatz matrix
OK core gate: Jensen-window PF Schur shape contract
OK core gate: Jensen-window PF column recurrence contract
OK core gate: Jensen-window PF column recurrence finite coverage
OK core gate: Arb Jensen-window column recurrence stress
OK core gate: Jensen-window PF reciprocal coefficient extended stress
OK core gate: Jensen-window PF reciprocal positivity route matrix
OK core gate: Jensen-window PF reciprocal fraction scout
OK core gate: Jensen-window PF reciprocal signed J-fraction scout
OK core gate: Jensen-window PF reciprocal signed Jacobi beta scout
OK core gate: Jensen-window PF reciprocal Motzkin path obstruction scout
OK core gate: Jensen-window PF reciprocal Motzkin parity-lift obstruction scout
OK core gate: Jensen-window PF signed J-fraction theorem target
OK core gate: Jensen-window PF modified signed-model target
OK core gate: Jensen-window PF oscillatory resolvent fit matrix
OK core gate: Jensen-window PF positive readout theorem target
OK core gate: Jensen-window PF positive spectral moment obstruction
OK core gate: Jensen-window PF nonordinary positive transform ansatz matrix
OK core gate: Jensen-window PF nonpower functional low-degree scout
OK core gate: Jensen-window PF nonpower functional cone candidate matrix
OK core gate: Jensen-window PF Cauchy-Binet cone frontier matrix
OK core gate: Jensen-window PF monotone contraction frontier scout
OK core gate: Jensen-window PF monotone-contraction column extension scout
OK core gate: Jensen-window PF monotone-contraction sparse degree-6 scout
OK core gate: Jensen-window PF monotone-contraction sparse degree-7 frontier scout
OK core gate: Jensen-window PF monotone-contraction sparse degree-7 subdivision scout
OK core gate: Jensen-window PF monotone-contraction all-m counterexample
OK core gate: Jensen-window PF monotone contraction theorem target
OK core gate: Jensen-window PF heat-flow monotone closure scout
OK core gate: Jensen-window PF heat-flow boundary threshold lemma
OK core gate: Jensen-window PF kernel Mellin upper-wall certificate
OK core gate: Jensen-window PF log-concave Mellin monotone-wall countermodel
OK core gate: Jensen-window PF T=1156 monotone-wall counterexample certificate
OK core gate: Jensen-window PF negative-lambda kernel summand-shift lemma
OK core gate: Jensen-window PF negative-lambda first-summand dominance certificate
OK core gate: Jensen-window PF negative-lambda -100 k320 collar extension certificate
OK core gate: Jensen-window PF negative-lambda -100 full cone-entry certificate
OK core gate: Jensen-window PF infinite heat-flow cone-invariance certificate
OK core gate: Jensen-window PF defect complete-monotonicity scout
OK core gate: Jensen-window PF multiplier complete-monotonicity frontier scout
OK core gate: Jensen-window PF multiplier Hausdorff-uniqueness bridge
OK core gate: Jensen-window PF multiplier leading-atom bound certificate
OK core gate: Jensen-window PF multiplier unit-atomic obstruction certificate
OK core gate: Jensen-window PF heat-flow Jensen hierarchy lemma
OK core gate: Jensen-window PF cubic reciprocal-defect invariance lemma
OK core gate: Jensen-window PF cubic lambda=-100 tail-entry certificate
OK core gate: Jensen-window PF cubic forward-uniform tail certificate
OK core gate: Jensen-window PF quartic boundary-flow obstruction
OK core gate: Jensen-window PF quartic double-root threshold lemma
OK core gate: Jensen-window PF quartic-quintic polar-contact lemma
OK core gate: Jensen-window PF cofinal-degree polar-closure lemma
OK core gate: Jensen-window PF cofinal scaling-limit equivalence gate
OK core gate: Jensen-window PF polar heat-collision cascade lemma
OK core gate: Jensen-window PF scaled double-zero boundary-layer lemma
OK core gate: Jensen-window PF Newman root external-field lemma
OK core gate: Jensen-window PF Newman classical-field balance gate
OK core gate: Jensen-window PF Newman local odd-count reduction lemma
OK core gate: Jensen-window PF Newman boundary-energy direction gate
OK core gate: Jensen-window PF Newman positive-boundary attainment lemma
OK core gate: Jensen-window PF Newman strict-Laguerre correlation target
OK core gate: Jensen-window PF Newman Polymath-15 oscillatory zeta handoff theorem
OK core gate: Jensen-window PF Newman Polymath-15 cancellation/zero-free wall gate
OK core gate: Jensen-window PF Newman Polymath-15 Gaussian/Legendre duality gate
OK core gate: Jensen-window PF Newman Polymath-15 ANTEDB beta-frontier audit
OK core gate: Jensen-window PF Newman Polymath-15 critical scaled coercivity target
OK core gate: Jensen-window PF Newman correlation hierarchy Gaussian-mixture gate
OK core gate: Jensen-window PF Newman theta-summand spectral-square gate
OK core gate: Jensen-window PF Newman theta/Bessel higher-shift regularization gate
OK core gate: Jensen-window PF Newman theta cell-renormalization gate
OK core gate: Jensen-window PF Newman Gasper fake-Xi remainder gate
OK core gate: Jensen-window PF Newman Gasper residual two-block gate
OK core gate: Jensen-window PF Newman classical three-block residual gate
OK core gate: Jensen-window PF Newman signed universal-factor residual gate
OK core gate: Jensen-window PF Laguerre scale-mixture gate
OK core gate: Jensen-window PF rank-two boundary-family lemma
OK core gate: Jensen-window PF multiplier counting-measure target
OK core gate: Jensen-window PF Mellin multiplier power-sum obstruction
OK core gate: Jensen-window PF negative-lambda first-summand saddle-wall closure
OK core gate: Jensen-window PF negative-lambda first-summand cumulant bridge
OK core gate: Jensen-window PF negative-lambda first-summand leading-saddle certificate
OK core gate: Jensen-window PF negative-lambda first-summand paired-remainder certificate
OK core gate: Jensen-window PF negative-lambda first-summand paired-remainder ray certificate
OK core gate: Jensen-window PF heat-flow ratio cone invariance lemma
OK core gate: Jensen-window PF heat-flow cone-entry asymptotic target
OK core gate: Jensen-window PF Phi Taylor cone-entry sign scout
OK core gate: Jensen-window PF negative-lambda cone-entry prefix scout
OK core gate: Jensen-window PF negative-lambda cone-entry prefix k30 scout
OK core gate: Jensen-window PF negative-lambda cone-entry prefix k50 scout
OK core gate: Jensen-window PF negative-lambda cone-entry prefix k60 scout
OK core gate: Jensen-window PF negative-lambda cone-entry prefix k80 scout
OK core gate: Jensen-window PF negative-lambda cone-entry prefix k100 scout
OK core gate: Jensen-window PF negative-lambda cone-entry prefix k150 scout
OK core gate: Jensen-window PF negative-lambda cone-entry prefix k200 scout
OK core gate: Jensen-window PF negative-lambda finite-collar contract
OK core gate: Jensen-window PF negative-lambda finite-collar k30 contract
OK core gate: Jensen-window PF negative-lambda finite-collar k50 contract
OK core gate: Jensen-window PF negative-lambda finite-collar k60 contract
OK core gate: Jensen-window PF negative-lambda finite-collar k80 contract
OK core gate: Jensen-window PF negative-lambda finite-collar k100 contract
OK core gate: Jensen-window PF negative-lambda finite-collar k150 contract
OK core gate: Jensen-window PF negative-lambda finite-collar k200 contract
OK core gate: Jensen-window PF negative-lambda tail-barrier scout
OK core gate: Jensen-window PF negative-lambda tail-barrier k30 scout
OK core gate: Jensen-window PF negative-lambda tail-barrier k50 scout
OK core gate: Jensen-window PF negative-lambda tail-barrier k60 scout
OK core gate: Jensen-window PF negative-lambda tail-barrier k80 scout
OK core gate: Jensen-window PF negative-lambda tail-barrier k100 scout
OK core gate: Jensen-window PF negative-lambda tail-barrier k150 scout
OK core gate: Jensen-window PF negative-lambda tail-barrier k200 scout
OK core gate: Jensen-window PF negative-lambda scaled-defect frontier k50 scout
OK core gate: Jensen-window PF negative-lambda scaled-defect frontier k60 scout
OK core gate: Jensen-window PF negative-lambda scaled-defect frontier k80 scout
OK core gate: Jensen-window PF negative-lambda scaled-defect frontier k100 scout
OK core gate: Jensen-window PF negative-lambda scaled-defect frontier k150 scout
OK core gate: Jensen-window PF negative-lambda scaled-defect frontier k200 scout
OK core gate: Jensen-window PF negative-lambda defect-recurrence scout
OK core gate: Jensen-window PF negative-lambda log-curvature bridge
OK core gate: Jensen-window PF negative-lambda bounded log-curvature target
OK core gate: Jensen-window PF negative-lambda bounded log-curvature k300 obstruction
OK core gate: Jensen-window PF negative-lambda Gaussian curvature matrix
OK core gate: Jensen-window PF negative-lambda signed Gaussian perturbation matrix
OK core gate: Jensen-window PF negative-lambda uniform remainder target
OK core gate: Jensen-window PF negative-lambda Taylor moment budget
OK core gate: Jensen-window PF negative-lambda high-order Taylor scout
OK core gate: Jensen-window PF negative-lambda defect-tail theorem target
OK core gate: Jensen-window PF negative-lambda half-width tail target
OK core gate: Jensen-window PF negative-lambda adaptive scaled-defect target
OK core gate: Jensen-window PF negative-lambda adaptive envelope matrix
OK core gate: Jensen-window PF negative-lambda adaptive envelope obligations
OK core gate: Jensen-window PF negative-lambda raw-moment bridge matrix
OK core gate: Jensen-window PF negative-lambda raw-ratio decrement-corridor scout
OK core gate: Jensen-window PF negative-lambda k300 precision-repair audit
OK core gate: Jensen-window PF negative-lambda raw-log decrement bridge
OK core gate: Jensen-window PF negative-lambda coefficient-curvature corridor bridge
OK core gate: Jensen-window PF negative-lambda linear curvature-barrier scout
OK core gate: Jensen-window PF negative-lambda scaled-curvature monotonicity target
OK core gate: Jensen-window PF negative-lambda scaled-curvature continuous bridge
OK core gate: Jensen-window PF negative-lambda scaled-curvature log-ceiling bridge
OK core gate: Jensen-window PF negative-lambda relative-Gaussian curvature bridge
OK core gate: Jensen-window PF negative-lambda relative-Gaussian Taylor stencil scout
OK core gate: Jensen-window PF negative-lambda relative-Gaussian stencil remainder obligations
OK core gate: Jensen-window PF negative-lambda relative-Gaussian pointwise tail budget
OK core gate: Jensen-window PF negative-lambda relative-Gaussian next-increment stencil stress
OK core gate: Jensen-window PF negative-lambda relative-Gaussian degree-16 stencil continuation
OK core gate: Jensen-window PF negative-lambda relative-Gaussian degree-16 collar scan
OK core gate: Jensen-window PF negative-lambda relative-Gaussian degree-16 real-T collar scout
OK core gate: Jensen-window PF negative-lambda relative-Gaussian degree-16 Arb real-T collar certificate
OK core gate: Jensen-window PF negative-lambda relative-Gaussian degree-40 Arb collar ladder stress
OK core gate: Jensen-window PF negative-lambda relative-Gaussian degree-40 residual tail budget
OK core gate: Jensen-window PF negative-lambda relative-Gaussian formal-tail obstruction scout
OK core gate: Jensen-window PF negative-lambda relative-Gaussian asymptotic remainder target
OK core gate: Jensen-window PF negative-lambda relative-Gaussian actual endpoint remainder scout
OK core gate: Jensen-window PF negative-lambda relative-Gaussian cancellation-reduced remainder grid scout
OK core gate: Jensen-window PF negative-lambda relative-Gaussian intervalization target
OK core gate: Jensen-window PF negative-lambda relative-Gaussian Phi tail bound scout
OK core gate: Jensen-window PF negative-lambda relative-Gaussian node-c0 range certificate
OK core gate: Jensen-window PF negative-lambda relative-Gaussian Phi-tail grid certificate
OK core gate: Jensen-window PF negative-lambda relative-Gaussian quadrature ladder scout
OK core gate: Jensen-window PF negative-lambda relative-Gaussian quadrature-remainder route matrix
OK core gate: Jensen-window PF negative-lambda relative-Gaussian worst-row far-tail split certificate
OK core gate: Jensen-window PF negative-lambda relative-Gaussian worst-row compact-interval integration scout
OK core gate: Jensen-window PF negative-lambda relative-Gaussian worst-row Chebyshev panel-moment scout
OK core gate: Jensen-window PF negative-lambda relative-Gaussian worst-row Arb Chebyshev interpolant-moment scout
OK core gate: Jensen-window PF negative-lambda relative-Gaussian worst-row interpolation-remainder route matrix
OK core gate: Jensen-window PF negative-lambda relative-Gaussian endpoint parity-repair matrix
OK core gate: Jensen-window PF negative-lambda relative-Gaussian endpoint x-panel route matrix
OK core gate: Jensen-window PF negative-lambda relative-Gaussian endpoint x-moment Taylor certificate
OK core gate: Jensen-window PF negative-lambda relative-Gaussian worst-row compact x-moment Taylor certificate
OK core gate: Jensen-window PF negative-lambda relative-Gaussian worst-row full-expectation certificate
OK core gate: Jensen-window PF negative-lambda relative-Gaussian all-row direct expectation certificate
OK core gate: Jensen-window PF negative-lambda relative-Gaussian recorded-grid stencil composition certificate
OK core gate: Jensen-window PF negative-lambda relative-Gaussian finite-collar-segment stencil certificate
OK core gate: Jensen-window PF negative-lambda relative-Gaussian full-kernel evenness/Cauchy lemma
OK core gate: Jensen-window PF negative-lambda relative-Gaussian full-real-T fixed-k stencil certificate
OK core gate: Jensen-window PF negative-lambda relative-Gaussian first-omitted denominator certificate
OK core gate: Jensen-window PF negative-lambda relative-Gaussian coefficient-core propagation certificate
OK core gate: Jensen-window PF negative-lambda relative-Gaussian worst-row Laguerre root-bracket certificate
OK core gate: Jensen-window PF negative-lambda relative-Gaussian worst-row Christoffel-weight midpoint scout
OK core gate: Jensen-window PF negative-lambda relative-Gaussian worst-row Christoffel-weight interval certificate
OK core gate: Jensen-window PF negative-lambda relative-Gaussian worst-row finite-part weighted-sum interval certificate
OK core gate: Jensen-window PF negative-lambda relative-Gaussian worst-row finite-plus-tail budget certificate
OK core gate: Jensen-window PF negative-lambda raw-moment obstruction matrix
OK core gate: Jensen-window PF negative-lambda zeta-specific raw-corridor target
OK core gate: Jensen-window PF negative-lambda -100 raw-corridor certificate
OK core gate: Jensen-window PF negative-lambda -100 adaptive-defect certificate
OK core gate: Jensen-window PF monotone contraction stress
OK core gate: Jensen-window PF state-space sign-lift obstruction scout
OK core gate: Jensen-window PF Cauchy-Binet low-degree scout
OK core gate: Jensen-window PF log-concavity frontier scout
OK core gate: Jensen-window PF ratio-condition scout
OK core gate: Jensen-window PF contraction-log-concavity scout
OK core gate: sign-regularity theorem fit matrix
OK core gate: positive Schur-specialization target note
OK core gate: Edrei-log sign diagnostics
OK core gate: Edrei power-Hankel diagnostics
OK core gate: Edrei midpoint frontier non-promotion guard
OK core gate: Edrei power-Hankel boundary repair manifest
OK core gate: Edrei moment-recurrence scout manifest
FULL EXECUTION: validated 301/301 registered core proof-programme gates
```

## Included Gates

```text
countermodel proof-safety gates:
  validates 11 exact/logical traps, including the finite Jensen-window
  rectangle extension gate

result-language boundary scan:
  scans output markdown for obvious unqualified global-proof overclaims

output artifact status manifest:
  requires every output markdown file to declare Date, Status, artifact kind,
  and proof-boundary language near the top

output reference integrity:
  checks markdown references to scripts, output notes, and concrete result
  files; current run has 0 missing required paths and 3 planned deliverables

proof-claim ledger:
  validates 273 classified claims, including 7 open theorem targets that remain
  explicitly unproved

signed-Hankel/Jensen dependency graph:
  validates outputs/signed_hankel_jensen_dependency_graph.md and
  work/rh_compute/results/signed_hankel_jensen_dependency_graph.json,
  checking that finite signed-Hankel/Jensen evidence feeds into open theorem
  targets and not directly into the Lambda <= 0 conclusion

signed-Hankel finite certificates:
  validates 2,500 finite signed-Hankel certificates

Toeplitz/PF finite certificates:
  validates 95 promoted positive summaries and 1 negative control

Toeplitz/Jacobi-Trudi reindexing:
  validates N=10, orders<=5 exact map over 124,129 minors
  records 39,094 nonzero maps and 85,035 structural zeros

Hankel sign-consistency reduction point audits:
  validates 20 exact-rationalized finite reshaped-Hankel point audits
  for the five-lambda grid, orders k=2..5, and N=18

Hankel sign-consistency reduction finite certificates:
  validates 689,795 Arb/enclosure-backed reshaped-Hankel finite certificates
  for the five-lambda grid, orders k=2..7, and N=20

shifted Hankel staircase finite certificates:
  validates 3,154,515 Arb/enclosure-backed shifted reshaped-Hankel finite
  certificates for the five-lambda grid and shifts n=0..20:
  k=2..5 at N=18, k=6 at N=16, k=7 at N=15, and k=8 at N=14

Jensen/Hankel bridge algebra gate:
  validates the exact degree-2 signed-Hankel/Jensen identity and an exact
  degree-3 finite countermodel blocking low-order finite-sign promotion

Jensen-window PF obligation algebra gate:
  validates exact low-degree Jensen-window PF obligations: degree 2 matches the
  signed-Hankel threshold, while degree 3 and degree 4 introduce additional
  banded Toeplitz obligations and finite low-order countermodel failures

Arb Jensen-window PF obligation finite diagnostics:
  validates 1,470 selected Arb/enclosure-backed degree-3/4 Jensen-window
  contiguous Toeplitz determinants for the five-lambda grid and shifts n=0..20

Arb Jensen-window Sturm hyperbolicity finite diagnostics:
  validates 210 Arb/Sturm degree-3/4 Jensen-window positive-root counts for
  Q_{d,n,lambda}(y)=P_{d,n,lambda}(-y) on the five-lambda grid and shifts
  n=0..20

Arb Jensen-window Sturm degree-5 hyperbolicity finite diagnostics:
  validates 105 Arb/Sturm degree-5 Jensen-window positive-root counts for
  Q_{5,n,lambda}(y)=P_{5,n,lambda}(-y) on the five-lambda grid and shifts
  n=0..20

Arb Jensen-window Sturm degree-6 through degree-12 hyperbolicity finite diagnostics:
  validates 735 Arb/Sturm positive-root counts across every degree d=6..12,
  the five-lambda grid, and shifts n=0..20, bringing the d=3..12 finite total
  to 1050 certified rows

finite Sturm-to-PF Jensen-window consequences:
  validates the finite Polya-frequency consequence for 1050 checked
  Jensen windows across degrees d=3..12, five lambdas, and shifts n=0..20

signed-Hankel/Jensen bridge target specification:
  validates that the active bridge target is stated as an all-order,
  all-shift open theorem with explicit low-degree obstructions

Jensen-window PF bridge target:
  validates the total-positivity reformulation that every binomially weighted
  Jensen window B^{d,n,0}_j=binom(d,j) A_{n+j}(0) must be finite PF-infinity,
  while preserving that this is an open theorem target

Jensen-window PF bridge obligation ledger:
  validates a 10-row obligation decomposition for the Jensen-window PF bridge:
  exact equivalence/contact rows, finite evidence rows, 3 open obligations,
  a conditional limiting row, and rejection/route-separation guards

Jensen-window PF theorem machinery fit matrix:
  validates a 7-row source-anchored theorem-family audit for jwpf_06,
  including ASW/Edrei, Schoenberg, Karlin/Cauchy-Binet,
  Polya-Schur/preserver, sign-regular matrix, downstream Laguerre-Polya, and
  rejected shortcut rows, with 0 ready-to-apply rows

Jensen-window PF sign-regular transfer gap matrix:
  validates a 9-row transfer-gap diagnostic combining exact degree-2
  signed-Hankel contact, degree-3/4 countermodel gates, 3 open all-order
  requirements, and 3 rejected shortcuts, with 0 ready-to-apply rows

Jensen-window PF factorial multiplier split audit:
  validates a 5-row route-separation audit for gamma_k=k!/(2*k)!,
  recording 315 raw degree-2 anti-hyperbolic rows and 315 normalized
  degree-2 positive rows, with 0 ready-to-apply rows

Jensen-window PF reciprocal-gamma mixture sign gate:
  validates Karlin's all-order fixed-scale reciprocal-gamma signature, an
  exact independent-row scale integral, two exact positive-mixture failure
  gates, and the equivalence between the completed Xi order-two cone and the
  tilted concentration bound CV_n^2<=2/(2n+1)

Jensen-window PF reciprocal-defect compound order-three gate:
  validates the exact contiguous 3x3 reciprocal-defect margin C_n, its
  increment curvature budget and boundary benchmark, plus a strict rational
  countermodel satisfying every existing ratio/adaptive/cubic cone while
  reversing the required signed-Hankel sign; composes the separate entry and
  forward-invariance theorems for the actual Xi trajectory

Jensen-window PF negative-lambda -100 compound order-three entry certificate:
  validates 318 Arb prefix compound margins, the exact scaled-defect tail
  bound C_n>57613471/66107054971 for n>=318, and the resulting all-shift
  contiguous 3x3 entry theorem at lambda=-100

Jensen-window PF compound order-three forward-invariance certificate:
  validates the cooperative flow C_n'/r_(n+2)=alpha_n*C_(n+1)+beta_n*C_n,
  compact-tail coefficient cancellation, and a weighted infinite maximum
  principle proving D_(3,n)(0)<0 for every shift; the following secant theorem
  closes noncontiguous order three

Jensen-window PF order-three noncontiguous secant-transfer lemma:
  validates the planar orientation identity and weighted-secant transfer,
  proving every arbitrary-column reshaped-Hankel sign in orders two and three
  at lambda=0; downstream total positivity now closes arbitrary columns at
  order four as well

Jensen-window PF compound order-four condensation gate:
  validates the Desnanot-Jacobi/log-concavity coordinate, 317 positive
  lambda=-100 prefix margins, and a strict rational lower-order countermodel
  with H_(4,0)<0; downstream the all-index G-gap tail, flow, and
  arbitrary-column transfer are closed

Jensen-window PF compound order-four first-summand curvature bridge:
  validates the stable gap factorization, the real J_1(t)>=1/(7t) floor,
  the same-point derivative localization, and the exact 2/(5k^2) full-kernel
  perturbation budget; downstream exact-corridor theorems now discharge the
  former analytic u>=2 condition

Jensen-window PF compound order-four localized-curvature compact certificate:
  validates a 107,452-tile outward-rounded Arb cache and reconstructs 1,073
  positive real-parameter blocks proving K_1(t)<=7/(2t^2) for
  319<=t<=V'(2); downstream finite and ray theorems close the rest

Jensen-window PF compound order-four Gaussian cumulant ray target:
  validates the exact epsilon-six cumulant polynomials through order eight,
  their alternating factorial leading signature, explicit candidate
  corridors, two proved global formal-corridor components, and seven positive
  conditional full-collar boxes; exact-minus-formal errors remain open

Jensen-window PF compound order-four formal cumulant corridor certificate:
  validates 1,800,000 adjacent Arb blocks covering 2<=u<=20 and all seven
  epsilon-six formal cumulant corridors, with weakest margin above 0.01843

Jensen-window PF compound order-four formal cumulant asymptotic ray certificate:
  validates seven coefficient-positive buffered leading-corridor gates,
  fourteen exact potential-jet remainder sign gates, and the transfer
  |R_r^[6]-F_r|<1/(20u), closing the formal model on u>=20

Jensen-window PF compound order-four finite exact-corridor curvature theorem:
  validates 20,700 correlated mode blocks and 41,400 shifted-collar gates,
  proving j_0>E_0 and t^2 U(t)<7/2 throughout 2<=u<=20

Jensen-window PF compound order-four exact-corridor curvature ray theorem:
  validates the coefficient-positive geometry, seven normalized H boxes,
  six logarithmic-defect derivative bounds, and the strict rational ceiling
  3011223637/866377000<7/2; finite plus ray close corridor-to-U globally

Jensen-window PF compound order-four lambda=-100 entry theorem:
  validates the 317-row repaired prefix, analytic penalty tail, strict
  scaled-defect buffer, and H_(4,n)(-100)>0 for every n>=0; the uniform-tail
  theorem below now supplies the formerly open forward propagation

Jensen-window PF compound order-four forward-flow reduction:
  validates three exact affine-Hankel/Plucker identities, the positive
  one-sided order-four determinant and stable-gap systems, uniform tail
  attainment, and reduction of forward propagation to one effective-diagonal
  coefficient ceiling; the finite-confinement reduction below provides an
  alternative route that does not require that ceiling

Arb Xi lambda-zero order-four prefix certificate:
  validates 507 outward-rounded coefficient balls from a direct 24,576-bit
  Xi series and proves every raw H_4 determinant and stable margin strictly
  positive on 0<=n<=500, with zero inconclusive rows

Jensen-window PF compound order-four lambda-zero eventual positivity:
  validates the exact Xi normalization, seven determinant-series coefficients,
  universal leading term 768*G_2^6*h^6, rigorous H_4 positivity for all
  0<=n<=500, and eventual H_4 positivity; its direct effective finite splice
  remains open but is no longer needed for the all-shift conclusion

Jensen-window PF lambda-zero first-summand dominance transfer:
  validates covariance monotonicity across 0<=T<=100, transfers the certified
  first-summand error from T=100 through lambda zero, and preserves the exact
  order-four perturbation budget; its direct curvature-tail handoff remains a
  separate, superseded route

Jensen-window PF compound order-four eventual-tail forward reduction:
  validates variation of constants and finite backward induction, reducing
  all-shift propagation to a uniform eventual positive tail without any
  spatial coefficient ceiling

Jensen-window PF compound order-four uniform-heat tail reduction:
  validates the universal 768*G_2^6*h^6 determinant term under a compact-
  uniform degree-seven ratio contract and isolates the first-summand heat tilt
  as the last input to that contract

Jensen-window PF uniform superpolynomial first-summand dominance:
  validates uniform superpolynomial suppression of all higher theta summands
  and every fixed local logarithmic difference on 0<=T<=100

Jensen-window PF uniform first-summand heat-tilt asymptotic theorem:
  validates the compact-uniform suitability audit for O'Sullivan's published
  all-order saddle theorem, the exact first-summand integral reduction, and
  the required difference bounds through order seven

Jensen-window PF compound order-four uniform-heat forward invariance:
  composes the uniform ratio contract, universal determinant cancellation,
  all-shift lambda=-100 entry, and finite cooperative propagation to prove
  every contiguous H_(4,n)(lambda)>0 on -100<=lambda<=0

Jensen-window PF order-four noncontiguous total-positivity transfer:
  validates four exact reversal orders, 240 solid-block maps, the rectangular
  Gasca-Pena initial-minor criterion, and 1,020 rational benchmark minors;
  proves every arbitrary-column order-four sign and the fixed-order transfer
  from contiguous layers to arbitrary columns

Jensen-window PF compound order-five uniform-tail and flow reduction:
  the separately run slow gate validates twelve compact-family suitability
  coefficients, eleven Lambert orders, nine-node Newton transfer, the exact
  120-permutation term 294912*G_2^10*h^10, cooperative order-five flow, and
  finite confinement to lambda=-100 entry

Jensen-window PF compound order-five lambda=-100 prefix:
  validates the exact H_(5,n)=W_n*J_n stable factorization, 325 hashed
  outward-rounded coefficient balls, 317 positive J_n and relative margins,
  and H_(5,n)(-100)>0 for every 0<=n<=316

Jensen-window PF compound order-five lambda=-100 tail curvature reduction:
  validates J_n>0 iff C_n<-4log(x_(n+4)), the defect anchor, and the
  coefficient-positive comparison showing C_n<=100/(n+4)^2 is sufficient on
  n>=317; the downstream nested-curvature chain now supplies that ceiling

Jensen-window PF compound order-five first-summand curvature bridge:
  validates two positive nested stable floors, the exact q'' identity, and a
  degree-52 coefficient-positive full-kernel perturbation reserve; splits the
  scalar budget as 37/k^2 plus a 63/k^2 first-summand contribution

Jensen-window PF compound order-five nested-curvature compact certificate:
  validates 36 outward-rounded blocks built from the hashed 107452-tile
  H-derivative cache; proves q_1''(t)<=60/t^2 throughout 320<=t<=V'(2)

Jensen-window PF compound order-five nested-curvature finite ray:
  validates 100 rigorous mode-two extension tiles and 1850 exact-cumulant
  corridor blocks; proves the same curvature ceiling for every 2<=u<=20

Jensen-window PF compound order-five nested-curvature asymptotic ray:
  validates the normalized-H collar, analytic stable-log majorant, and one
  dimensionless interval box; proves t^2q_1''(t)<9.159 for every u>=20

Jensen-window PF compound order-five lambda=-100 entry:
  composes the three continuous ranges, the 37+63 full-kernel transfer, and
  the rigorous prefix to prove H_(5,n)(-100)>0 for every integer n>=0

Jensen-window PF compound order-five uniform-heat forward invariance:
  composes endpoint entry, the uniform eventual tail, cooperative variation
  of constants, and the fixed-order initial-minor theorem; proves contiguous
  and arbitrary-column signed-Hankel order five on the full heat segment

Jensen-window PF compound order-six uniform-tail and flow reduction:
  validates 17 compact-family suitability coefficients, 16 Lambert orders,
  ten Newton coefficients, 720 determinant permutations, and 684 weighted
  monomials; proves the uniform signed tail and the exact cooperative
  order-six recursion conditional only on lambda=-100 entry

Jensen-window PF compound order-six lambda=-100 prefix:
  validates 327 coefficient records, 317 strictly positive relative H_5
  margins and Q_6 signs, and 39 repaired coefficients; proves the endpoint
  prefix through n=316 with weakest certified margin above 7/1000

Jensen-window PF compound order-six lambda=-100 tail curvature reduction:
  validates the two exact stable factorizations that reduce the analytic tail
  from n>=317 to one global first-summand scalar-curvature ceiling

Jensen-window PF inverse-seventh-power first-summand dominance:
  validates six exact rows and eleven positive gates; supplies the complete
  inverse-seventh-power first-summand dominance transfer from k>=316

Jensen-window PF compound order-six first/full curvature bridge:
  validates the exact first-summand and full-kernel curvature transfer,
  including the cancellation-preserving perturbation budget used by the
  endpoint tail theorem

Jensen-window PF compound order-six high-cumulant corridor:
  validates 109 exact terms for formal orders nine and ten and proves the two
  high-cumulant corridors used on the long ray

Jensen-window PF compound order-six nested-curvature compact certificate:
  validates 38 outward-rounded blocks and proves p_1''(t)<=200/t^2 throughout
  the compact interval 321<=t<=V'(2)

Jensen-window PF compound order-six nested-curvature finite ray:
  validates 17,999 rational mode blocks and proves the same curvature ceiling
  throughout the complete ray segment 2<=u<=20

Jensen-window PF compound order-six nested-curvature asymptotic ray:
  validates the normalized-H collar and one dimensionless interval theorem;
  proves t^2*p_1''(t)<22.769 for every u>=20

Jensen-window PF compound order-six lambda=-100 entry:
  composes the compact, finite-ray, and asymptotic curvature ranges with the
  full-kernel transfer, analytic n>=317 tail, and rigorous prefix to prove
  Q_(6,n)(-100)=-H_(6,n)(-100)>0 for every n>=0

Jensen-window PF compound order-six uniform-heat forward invariance:
  composes endpoint entry, the uniform signed tail, cooperative recursion,
  the completed order-five lower cone, and the fixed-order initial-minor
  theorem; proves contiguous and arbitrary-column signed-Hankel order six on
  the full heat segment and hands the fixed-order frontier to order seven

Jensen-window PF graded-kernel all-order Vandermonde lemma:
  validates 12 order rows, 46,233 permutation stress cases, and 169
  coefficient-valuation cells; proves the universal positive signed leading
  term and one compact-uniform eventual tail at every fixed order, with an
  order-dependent non-effective threshold

Jensen-window PF compound order-seven uniform-tail and flow reduction:
  specializes the universal tail to the positive signed coefficient
  52183852646400*G_2^21, validates the exact condensation coordinate and
  cooperative recursion, and records both the conditional forward theorem
  and exact lower-cone countermodel

Jensen-window PF compound order-seven lambda=-100 prefix:
  validates 327 coefficients, 317 positive Q_6 values, 315 strict relative
  margins and Q_7 signs, and twelve repaired coefficients; proves the rigorous
  endpoint prefix through n=314 with weakest margin above 9/1000

Jensen-window PF compound order-seven lambda=-100 tail curvature reduction:
  validates the fourth stable factorization, exact centered-curvature
  coordinate, and positive rational comparison; reduces every missing n>=315
  endpoint sign to the ceiling R_k<=900/k^2 for k>=321, later discharged by
  the continuous and full-kernel composition

Jensen-window PF inverse-eighth-power rebalanced first-summand dominance:
  validates fourteen strict endpoint and half-line derivative gates for the
  log(k)/10 split; proves the complete-to-first moment defect is below 2/k^8
  for every k>=300

Jensen-window PF compound order-seven first/full curvature bridge:
  validates the global fourth-gap floor, the full four-layer perturbation
  chain, and a degree-102 shifted numerator with 103 positive coefficients;
  proves |R_k-R_k^(1)|<262/k^2 and reduces the endpoint tail to the single
  continuous target r_1''(t)<=600/t^2 on t>=320

Jensen-window PF compound order-seven shifted-jet compact bridge:
  validates 186 contiguous rational blocks, eleven exact shifted point jets
  per anchor, and dimensionless H2-H21 common-collar remainders; proves
  r_1''(t)<=600/t^2 continuously on 320<=t<=1000 with worst scaled upper
  below 50.911

Jensen-window PF compound order-seven nested-curvature compact certificate:
  validates 107452 aligned H2-H12 mode tiles, a strict t+-5 collar, and 82
  contiguous adaptive rational blocks; proves r_1''(t)<=600/t^2 continuously
  on 1000<=t<=V'(2) with worst scaled upper below 358.733

Jensen-window PF compound order-seven high-cumulant coarse corridor:
  validates 72 exact formal terms and the unit-disk residual transfer; proves
  normalized kappa_11 and kappa_12 caps below 14001 on 2<=u<=20 and below
  700001 on u>=20

Jensen-window PF compound order-seven nested-curvature finite-ray certificate:
  validates the exact mode-two collar and 17,999 rational mode blocks; proves
  r_1''(t)<=600/t^2 for every 2<=u<=20 with worst scaled upper below 73.543

Jensen-window PF compound order-seven nested-curvature asymptotic-ray certificate:
  validates normalized H2-H12 boxes, the stable-log defect majorant, and one
  interval over 0<=1/t<=10^(-30); proves t^2r_1''(t)<55.541 for all u>=20

Jensen-window PF compound order-seven lambda=-100 entry:
  composes all four continuous ranges, the exact tent and full-kernel
  transfers, the analytic n>=315 tail, and the rigorous prefix; proves
  Q_(7,n)(-100)=-H_(7,n)(-100)>0 for every n>=0

Jensen-window PF compound order-seven uniform-heat forward invariance:
  composes endpoint entry, the uniform eventual tail, cooperative recursion,
  the completed order-six cone, and the fixed-order initial-minor theorem;
  proves contiguous and arbitrary-column signed-Hankel order seven on the
  full heat segment

Jensen-window PF compound order-eight uniform-tail and flow reduction:
  specializes the universal signed term to order eight, validates the exact
  Q_7 log-concavity coordinate and cooperative recursion, and records both
  the conditional forward theorem and an exact lower-cone countermodel

Jensen-window PF compound order-eight lambda=-100 prefix:
  validates 1,257 endpoint coefficients, 1,245 positive Q_7 values, and all
  1,243 relative Q_7 margins and Q_8 signs; proves the rigorous endpoint
  prefix through n=1242 with the weakest margin above 1/300

Jensen-window PF compound order-eight lambda=-100 tail-curvature reduction:
  validates the stable endpoint factorization, the exact logarithmic buffer,
  and the sufficient ceiling W_k<=4300/k^2 for every k>=1250

Jensen-window PF inverse-ninth-power first-summand dominance:
  validates all 14 analytic gates and proves the retained-summand defect
  0<=delta_k<2/k^9 for every integer k>=300

Jensen-window PF compound order-eight first/full curvature bridge:
  validates the fifth stable-gap floor, 134 positive transfer coefficients,
  and |W_k-W_k^(1)|<190/k^2 for every k>=1250

Jensen-window PF compound order-eight high-cumulant coarse corridor:
  validates the vanishing epsilon-ten formal terms and exact normalized
  cumulant caps in orders 13 and 14 on every saddle mode u>=2

Jensen-window PF compound order-eight shifted-jet certificate:
  validates 185 contiguous common-collar blocks and proves
  s_1''(t)<=2000/t^2 for every real 699<=t<=999

Jensen-window PF compound order-eight nested-curvature compact certificate:
  validates 96 adaptive H2-H14 blocks and proves s_1''(t)<=4000/t^2 from
  t=999 through the mode-two saddle boundary

Jensen-window PF compound order-eight nested-curvature finite-ray certificate:
  validates the exact mode-two collar and 17,999 rational mode blocks; proves
  s_1''(t)<=4000/t^2 for every saddle mode 2<=u<=20

Jensen-window PF compound order-eight nested-curvature asymptotic-ray certificate:
  validates normalized H2-H14 boxes and one interval over 0<=1/t<=10^(-30);
  proves t^2s_1''(t)<134.49 for every u>=20

Jensen-window PF compound order-eight lambda=-100 entry:
  composes all continuous ranges, tent integration, the full-kernel transfer,
  analytic n>=1243 tail, and rigorous prefix; proves
  Q_(8,n)(-100)=H_(8,n)(-100)>0 for every n>=0

Jensen-window PF compound order-eight uniform-heat forward invariance:
  composes endpoint entry, the uniform eventual tail, cooperative recursion,
  the completed order-seven cone, and the fixed-order initial-minor theorem;
  proves contiguous and arbitrary-column signed-Hankel order eight on the
  full heat segment

Jensen-window PF compound order-nine uniform-tail and flow reduction:
  specializes the universal signed term to order nine, validates the Q_8
  log-concavity coordinate and cooperative recursion, and retains an exact
  lower-cone countermodel

Jensen-window PF compound order-nine lambda=-100 prefix:
  validates 1,257 endpoint coefficients, 38 repaired rows, and every Q_9 sign
  through n=1240, with the weakest relative margin above 1/250

Jensen-window PF compound order-nine lambda=-100 tail-curvature reduction:
  validates the sixth stable factorization and reduces the endpoint tail to
  Y_k<=4900/k^2 for k>=1249

Jensen-window PF compound order-nine first/full curvature bridge:
  validates the sixth-gap floor, 169 positive transfer coefficients, and
  |Y_k-Y_k^(1)|<550/k^2 for k>=1251

Jensen-window PF compound order-nine high-cumulant coarse corridor:
  validates the exact fifteenth- and sixteenth-cumulant corridors on every
  saddle mode u>=2

Jensen-window PF compound order-nine exact-point H0-H8 cache:
  validates all 8,929 hash-bound jet rows used by the localized lower bridge

Jensen-window PF compound order-nine localized lower bridge:
  validates 279 root segments and 874 accepted blocks; proves
  w_1''(t)<=4200/t^2 on 1250<=t<=5700

Jensen-window PF compound order-nine nested-curvature compact certificate:
  validates 108 H2-H16 blocks and proves the same ceiling from t=5700 through
  the mode-two saddle boundary

Jensen-window PF compound order-nine nested-curvature finite-ray certificate:
  validates 17,999 rational mode blocks and proves the ceiling for 2<=u<=20

Jensen-window PF compound order-nine nested-curvature asymptotic-ray certificate:
  validates one normalized asymptotic interval and proves
  t^2*w_1''(t)<324.906<500 for every u>=20

Jensen-window PF compound order-nine global first-summand curvature:
  composes the four exact ranges and proves w_1''(t)<=4200/t^2 for every real
  t>=1250

Jensen-window PF compound order-nine finite endpoint splice:
  validates retained-integral balls for A_1257,A_1258 and proves the two
  missing signs n=1241,1242

Jensen-window PF compound order-nine lambda=-100 entry:
  composes the global curvature theorem, tent and full-kernel transfers,
  analytic tail, splice, and prefix; proves Q_(9,n)(-100)>0 for every n>=0

Jensen-window PF compound order-nine uniform-heat forward invariance:
  composes endpoint entry, the eventual tail, cooperative recursion, and the
  fixed-order initial-minor theorem; proves contiguous and arbitrary-column
  order nine on the full heat segment

Jensen-window PF all-order endpoint-to-heat reduction:
  validates the arbitrary-order affine and flag-Plucker identities, the
  order-dependent-tail quantifiers, and the exact equivalence between the
  heat hierarchy and the endpoint antecedent Q_(m,n)(-100)>0 for every
  m>=10,n>=0; the equivalence remains exact although the antecedent is false

Jensen-window PF endpoint deep-Schur coordinate:
  validates the normalized rectangular and arbitrary-column Jacobi-Trudi
  maps, the inverse map onto deep partitions, a rigorous endpoint PF_3
  counterexample, and the bounded primary-literature fit audit; records that
  the candidate hierarchy s_((N^m))(h)>0 for every m>=10,N>=m-1 is rejected

Jensen-window PF deep-Schur Toda/boundary gate:
  validates the exact rectangular Toda identity, 251 shifted-tail boundary
  checks, and 138 strict-Schur checks; blocks the zero-reset moving-tail PF
  shortcut and generic strict-Schur-to-Jensen implication, and diagnoses the
  rejected order-ten rectangle layer while leaving a weaker Xi/Phi bridge open

Jensen-window PF endpoint order-ten counterexample:
  independently checks direct Hankel, Jacobi-Trudi, Toda, and stable
  condensation coordinates; certifies four negative required deep rectangles
  at N=9,10,11,12 and rejects the proposed all-order endpoint hierarchy

Jensen-window PF order-moment transport fit gate:
  validates the exact Gamma-average reparametrization and the positive origin
  derivative that excludes complete monotonicity; blocks promotion of an
  ordinary positive-Hankel transport theorem into reciprocal sign regularity

Jensen-window PF structural ansatz matrix:
  validates a 6-row structural proof-search workbench for jwpf_06, checking
  positive Cauchy-Binet, planar-network/production-matrix, determinant-integral,
  preserver, direct-Hankel, and finite-grid ansatz rows against the exact
  degree-2/3/4 hard tests and finite countermodel kill gate

Jensen-window PF Schur shape contract:
  validates a bounded finite-support Jacobi-Trudi shape contract for
  h_j=binom(d,j)A_{n+j}; records 15,709 finite-band nonzero bounded shapes
  and 2 hard frontier column-shape obligations

Jensen-window PF column recurrence contract:
  validates the elementary-symmetric recurrence C_m=h_0^m e_m for the
  column-shape obligations and confirms the 2 hard frontier recurrence rows
  match exact negative rational countermodel values

Jensen-window PF column recurrence finite coverage:
  validates finite zeta-window support for the recurrence target: 1,470 direct
  positive Arb determinant rows, 210 hard recurrence rows, and 315 checked
  Sturm/PF windows, without promoting this to an all-order theorem

Arb Jensen-window column recurrence stress:
  validates 12,600 positive Arb recurrence rows for degrees d=3..8, sizes
  m=1..20, five lambda values, and shifts n=0..20

Jensen-window PF reciprocal coefficient extended stress:
  validates 72,600 positive normalized reciprocal coefficient rows for
  degrees d=2..12, sizes m=1..40, five lambda values, and shifts n=0..32

Jensen-window PF reciprocal positivity route matrix:
  validates a 9-row theorem-search matrix for the column recurrence target
  [t^m]1/H(-t)>=0, separating endpoint language, live renewal/signed-fraction/
  production-matrix candidates, standard positive fraction mismatch, Kaluza
  mismatch, rejected generic ratio shortcuts, finite recurrence stress, and
  the still-missing all-Schur lift

Jensen-window PF reciprocal fraction scout:
  validates 3 symbolic continued-fraction sign rows and 735 finite Arb
  zeta-window sign rows, rejecting the standard positive Stieltjes/Jacobi
  S-fraction or J-fraction route for E(t)=1/H(-t) while leaving signed or
  modified fraction routes open

Jensen-window PF reciprocal signed J-fraction scout:
  validates 2 symbolic signature rows, 3,675 finite signed reciprocal-Hankel
  determinant rows, and 2,940 finite signed-lambda rows, supporting the signed
  modified J-fraction target while preserving that no all-order theorem is
  proved

Jensen-window PF reciprocal signed Jacobi beta scout:
  validates 3 symbolic rows and 3,675 finite beta rows, with 2,940 positive
  rows, 630 negative beta_1 rows, and 105 terminal degree-2 zero-containing
  rows, sharpening the signed modified J-fraction target while preserving that
  no all-order theorem is proved

Jensen-window PF reciprocal Motzkin path obstruction scout:
  validates 3 symbolic rows, 735 finite mu_2 cancellation rows, and 630 finite
  beta_1 diagonal obstruction rows, rejecting the raw ordinary J-fraction
  Motzkin path model as manifest positivity while leaving modified signed
  models open

Jensen-window PF reciprocal Motzkin parity-lift obstruction scout:
  validates 3 symbolic rows and 5,145 finite same-length mixed-sign witness
  rows, rejecting global length-parity signs and diagonal sign conjugation as
  repairs of the raw ordinary J-fraction Motzkin path model while leaving
  state-space modified models open

Jensen-window PF signed J-fraction theorem target:
  validates a 7-row target specification for the missing theorem that would
  convert all-order signed reciprocal-Hankel/Jacobi signature into
  coefficientwise nonnegativity of E(t)=1/H(-t), with 0 ready-to-apply rows

Jensen-window PF modified signed-model target:
  validates a 9-row audit separating dead signed-model repairs from live
  conditional modified candidates: raw ordinary Motzkin/J-fraction positivity,
  diagonal sign conjugation, and global length-parity repairs are rejected,
  while state-space, modified production/Riordan/network, oscillatory
  resolvent, and Xi/Phi positive-cancellation routes remain conditional only

Jensen-window PF oscillatory resolvent fit matrix:
  validates an 8-row theorem-search matrix for the oscillatory/resolvent
  subroute, rejecting entrywise Jacobi powers, diagonal similarity,
  absolute-value majorants, classical oscillatory spectral conclusions,
  indefinite moment language, and finite signed patterns as standalone
  coefficient-positivity proofs, while leaving positive spectral-transform and
  Xi/Phi positive-kernel variants live conditionally

Jensen-window PF positive readout theorem target:
  validates an 8-row target for the exact positive scalar readout still needed
  after the oscillatory/resolvent fit matrix, leaving only a noncircular
  positive spectral transform or Xi/Phi-specific positive resolvent kernel
  live as foundational routes, now sharpened by the nonordinary transform
  ansatz matrix, and rejecting abstract wrappers, endpoint factorization,
  finite quadrature, raw signed readouts, and absolute-value majorants as
  standalone proofs

Jensen-window PF positive spectral moment obstruction:
  validates a 3-row symbolic and 735-row finite obstruction showing that the
  raw reciprocal coefficients mu_m cannot be ordinary power moments of a
  positive measure because Delta_2=-g_2<0, while leaving nonordinary positive
  transforms and Xi/Phi positive kernels open

Jensen-window PF nonordinary positive transform ansatz matrix:
  validates 8 ansatz rows with 0 ready-to-apply rows and 3 live ansatz rows,
  narrowing the surviving positive-readout search to non-power positive
  functionals, Xi/Phi positive kernels, or genuinely modified exact
  state-space transfer models

Jensen-window PF nonpower functional low-degree scout:
  validates 7 scout rows with 0 ready-to-apply rows and 1 live contract row,
  recomputing reciprocal coefficient formulas through mu_6 and signed
  composition counts through m=8 to show that the non-power functional route
  needs a genuine positive cone, basis, and functional absorbing signed
  low-degree cancellations

Jensen-window PF nonpower functional cone candidate matrix:
  validates 8 cone rows with 0 ready-to-apply rows and 2 live cone rows,
  rejecting raw g-coordinate, standalone ratio/log-concavity, tautological,
  and endpoint PF/LP cones while leaving Xi/Phi kernel and
  Cauchy-Binet/determinant-integral cone routes live

Jensen-window PF Cauchy-Binet cone frontier matrix:
  validates 8 frontier rows with 0 ready-to-apply rows and 2 live frontier
  rows, showing that the determinant-integral route cannot rest on selected
  low-degree Bernstein certificates or adjacent log-concavity and must instead
  construct a hard column-frontier determinant integral or all-shape
  Cauchy-Binet/Andreief kernel

Jensen-window PF monotone contraction frontier scout:
  validates 2 exact hard-row certificates with 88 positive Bernstein
  coefficients and 210 finite zeta diagnostic rows, showing that the first
  hard column-frontier polynomials are positive under x1 <= x2 <= x3 while
  the rational log-concavity countermodel violates that sharper condition

Jensen-window PF monotone-contraction column extension scout:
  validates 25 bounded degree-3/4/5 column rows with 3,329 positive Bernstein
  coefficients under monotone contractions, including 3 rows beyond the
  original first hard frontier and a degree-5 band through m=8

Jensen-window PF monotone-contraction sparse degree-6 scout:
  validates 10 bounded degree-6 column rows through m=10 with 63,347 strictly
  positive Bernstein coefficients under monotone contractions, using a sparse
  exact common-denominator Bernstein transform

Jensen-window PF monotone-contraction sparse degree-7 frontier scout:
  validates 9 bounded degree-7 column rows through m=9 with 670,891 strictly
  positive Bernstein coefficients, and records a degree-7 m=10 one-shot
  global Bernstein certificate obstruction with 126 negative coefficients

Jensen-window PF monotone-contraction sparse degree-7 subdivision scout:
  repairs the degree-7 m=10 one-shot global Bernstein obstruction by splitting
  s0 into 3 dyadic slabs and validating 785,400 strictly positive slab
  Bernstein coefficients

Jensen-window PF monotone-contraction all-m counterexample:
  validates an exact shift-0 infinite witness satisfying all pointwise
  lower/upper walls and monotone contractions, but with a negative degree-7
  m=11 normalized column recurrence; the propagated static cone alone is not
  an all-m column recurrence theorem

Jensen-window PF monotone contraction theorem target:
  is superseded by the validated full-cone entry and infinite heat-flow
  invariance certificates, which prove Delta^3 log A_{k-1}(lambda)>=0 for
  every finite lambda>=-100; the all-order PF/Jensen bridge remains open

Jensen-window PF heat-flow monotone closure scout:
  validates 4 exact lambda-flow algebra rows, 315 finite Arb threshold rows,
  and 305 finite Arb flow-bracket rows, isolating q >= (2*k-1)/(2*k+5) as a
  boundary condition without proving a closed differential inequality

Jensen-window PF heat-flow boundary threshold lemma:
  validates a 5-row exact lemma proving x_k >= (2*k-1)/(2*k+1), hence the
  heat-flow threshold x_k >= (2*k-1)/(2*k+5), from Phi positivity and
  Cauchy-Schwarz raw-moment log-convexity

Jensen-window PF kernel Mellin upper-wall certificate:
  validates global strict log-concavity of y->Phi(sqrt(y)) and applies
  Berwald-Borell to prove x_k(lambda)<=1 for every real lambda and k>=1

Jensen-window PF log-concave Mellin monotone-wall countermodel:
  validates an abstract log-concave density with x_2<x_1, blocking generic
  promotion from the Mellin upper wall to adjacent-k monotonicity

Jensen-window PF T=1156 monotone-wall counterexample certificate:
  validates four ACB actual-kernel coefficient enclosures and certifies
  x_121-x_120<-1.68e-8 at lambda=-1156, blocking all-k cone entry there and
  universal promotion of the fixed-k T>=1156 theorem

Jensen-window PF negative-lambda kernel summand-shift lemma:
  validates the exact identity phi_n(u)=n^(-1/2)phi_1(u+(log n)/2) and an
  Arb compact-tail bound below 2.122e-29 for n=2..20 at lambda=-100,
  k>=300, v<=3/2, isolating the dominant n=1 saddle

Jensen-window PF negative-lambda first-summand dominance certificate:
  validates the exact summand-ratio monotonicity and interval tail estimates
  proving M_k=M_k^(1)(1+delta_k), 0<=delta_k<=2/k^6, at lambda=-100 for
  every integer k>=300; consequently |L_k-L_k^(1)|<=16/(k-1)^6 for k>=301

Jensen-window PF negative-lambda -100 k320 collar extension certificate:
  validates 76 positive coefficients, 74 ratio-cone rows, and 73 adjacent-wall
  rows from the repaired ACB source, extending the certified finite collar
  through x_319>x_318 (adjacent index k=318)

Jensen-window PF negative-lambda -100 full cone-entry certificate:
  validates a precedence-merged repaired source with 321 positive coefficients,
  319 pointwise cone rows, and 318 adjacent prefix rows; exact all-k pointwise
  walls and the analytic adjacent tail prove full cone entry at lambda=-100

Jensen-window PF infinite heat-flow cone-invariance certificate:
  validates the exact defect evolution, a uniform zeroth-order coefficient
  bound, spatial-tail attainment, and the Dini derivative of the finite active
  minimum set; propagates the full ratio cone to every finite lambda>=-100,
  including lambda=0, while leaving the all-order PF/Jensen bridge open

Jensen-window PF defect complete-monotonicity scout:
  validates 3284 defect and 3288 negative-log-contraction Arb-positive
  alternating differences on five cached heat times, with both channels
  complete through order 8 and 838 high-order inconclusives; an exact -27/16
  cubic witness blocks promotion to the all-shape bridge

Jensen-window PF multiplier complete-monotonicity frontier scout:
  evaluates the tighter dps220 A_0..A_57 sources at 250-digit Arb precision
  and certifies all 7980 alternating log-contraction differences through the
  complete finite order-55 triangle, with no inconclusive or negative interval

Jensen-window PF multiplier Hausdorff-uniqueness bridge:
  proves the exact Hausdorff moment representation and uniqueness theorem,
  characterizes the unit-atomic density needed by the multiplier target, and
  retains a periodic-interpolation guard against importing the continuous
  Mellin obstruction without an additional uniqueness theorem

Jensen-window PF multiplier leading-atom bound certificate:
  conditionally brackets the strongest finite atom-kernel root at order 6,
  proves alpha_min>4.863538496 and N(11/2)<=1 for any target multiset, and
  leaves existence of an atom below 11/2 explicitly open

Jensen-window PF multiplier unit-atomic obstruction certificate:
  combines the order-6 atom cutoff with exact ratio monotonicity and an Arb
  consecutive-moment gap above 7.68e-4, ruling out the unit-atomic elementary
  multiplier product while leaving general multiplier sequences open

Jensen-window PF heat-flow Jensen hierarchy lemma:
  validates five exact degree-shift heat identities, the static -27/16 cubic
  guard, and a one-atom Hausdorff full-cone boundary witness with normalized
  frontier derivative 329728/2109375>0; even complete defect monotonicity plus
  the local heat ODE needs a stronger coupled higher-minor invariant

Jensen-window PF cubic reciprocal-defect invariance lemma:
  factors the cubic discriminant in q_k=(1-x_k)^(-1/2), making degree-3
  hyperbolicity exactly q_(k+1)-q_k<=1; the saturated heat boundary is inward
  iff the next shift satisfies the same condition. Arb certifies 318
  lambda=-100 prefix margins and 310 nonnegative-grid margins; the later
  tail-entry certificate closes the all-k lambda=-100 tail

Jensen-window PF cubic lambda=-100 tail-entry certificate:
  validates strict negative third log-cumulant on 4,074 compact Arb blocks and
  an explicit analytic ray, then proves d_m>=1/(5*m+1) and
  q_(k+1)-q_k<1 for every k>=319. Combined with 318 prefix margins this proves
  every shifted cubic Jensen polynomial hyperbolic at lambda=-100; the next
  certificate closes the uniform forward-heat spatial tail

Jensen-window PF cubic forward-uniform tail certificate:
  validates the exact one-sided q-increment flow and weighted source cap
  sqrt(k)*g_k'<=7*r_k. The entry estimate sqrt(k)*g_k<12 and an explicit
  coercive supersolution prove sup_[-100,L] g_k=O_L(k^(-1/2)), legitimizing
  the infinite first-crossing principle and proving every shifted cubic Jensen
  polynomial hyperbolic at lambda=0; degree four and all higher degrees remain
  open

Jensen-window PF quartic boundary-flow obstruction:
  validates the normalized quartic discriminant and exact heat derivative,
  together with a rational hyperbolic double-root boundary point satisfying
  four positive ratio margins and three strict cubic margins but having
  Q'/r_1<0. This blocks unchanged promotion of the cubic cone to degree four
  without claiming failure of the actual zeta trajectory

Jensen-window PF quartic double-root threshold lemma:
  validates symmetric double-root coordinates and the exact branch-aware
  condition (3*a^2-4*a+p)*(u-U)<=0, with U independent of the shift index.
  On the triple-root stratum first-order viability requires u=U and the first
  variation retains (1+a*w)^2; global quartic closure remains open

Jensen-window PF quartic-quintic polar-contact lemma:
  validates P_d=P_(d+1)-w*P_(d+1)'/(d+1), strict simplicity of nonroot polar
  zeros, and the multiplicity drop at shared roots. Consequently a quartic
  double root under a hyperbolic quintic extension forces quintic triple
  contact exactly at u=U; the noncircular infinite-degree closure remains open

Jensen-window PF cofinal-degree polar-closure lemma:
  validates one-step interlacing, finite polar descent, and the exact reduction
  from all finite degrees to an unbounded sequence of hyperbolic terminal
  degrees at each shift. The evidence audit records 1050 finite Sturm rows at
  degrees 3 through 12 while keeping 2,875 degree-3 through degree-12
  contraction rows separate from hyperbolicity; no cofinal sequence is proved

Jensen-window PF cofinal scaling-limit equivalence gate:
  validates P_(D,n)(z/D)->F_n(z) locally uniformly and the equivalence between
  cofinal fixed-shift hyperbolicity and Laguerre-Polya membership. This blocks
  treating a shift-zero cofinal theorem as a weaker endpoint surrogate or
  promoting the bounded degree-3 through degree-12 Sturm ladder into the limit

Jensen-window PF polar heat-collision cascade lemma:
  validates the general multiple-root heat jet and polar multiplicity lift.
  An all-hyperbolic upper tower would force a fixed-root multiplicity growing
  once per degree, hence an exponential-polynomial entire function. A genuine
  non-exponential-polynomial LP boundary is therefore strict in every fixed
  degree, and any failure must escape to unbounded degree

Jensen-window PF scaled double-zero boundary-layer lemma:
  validates the first two 1/D Jensen corrections and the universal local model
  eta^2+8*rho*tau-rho^2. It proves the D^(-3/2) finite root gap, the rho/(8D)
  collision shift, and the D^(-2) correction involving the regularized global
  root field U'(rho)/U(rho). The degree- and zero-height-uniform remainder
  estimate needed for Lambda<=0 remains open

Jensen-window PF Newman root external-field lemma:
  validates the canonical-product formulas for the regularized squared-zero
  field E_x and positive signed stiffness K_x, together with the exact pair
  center/gap flow and q=8*dt-16*K_x*dt^2+O(dt^3). Opposite-sign finite LP
  countermodels prove that generic coefficient positivity cannot control the
  D^(-2) collision term; the Xi-specific height-uniform balance remains open

Jensen-window PF Newman classical-field balance gate:
  validates the exact weighted perturbation identity and the -pi/8 continuum
  field of the Riemann-von Mangoldt density, matching the -pi/4 positive-time
  quantile drift. Published high-zero asymptotics confine fixed-time multiple
  zeros below exp(C/t), while exact bounded-location countermodels force either
  field sign. A lambda-uniform reciprocal-gap or compact-exclusion theorem is
  still required

Jensen-window PF Newman local odd-count reduction lemma:
  validates the Stieltjes outer-field estimate and composes it with the
  published uniform zero count. Outside H=log(4c)^2 the error is O(1/log c),
  leaving only the inverse-square weighted odd local count. An exact even
  heat-flow polynomial has the classical field and drift but positive pair
  birth, so local field balance and a separate Xi collision obstruction are
  both still required

Jensen-window PF Newman boundary-energy direction gate:
  validates the higher pair-gap jet and the universal 1/(8*tau) collision
  interaction. Nonnegative forward energy decay is compatible with birth from
  an infinite trace, while finite endpoint-integrated energy would exclude it.
  The published Rodgers-Tao theorem assumes Lambda<0 and stays away from the
  boundary, so it cannot be reversed into the missing upper direction

Jensen-window PF Newman positive-boundary attainment lemma:
  composes the absolute strip and uniform high-zero theorems to force any
  hypothetical positive boundary to have a finite real multiple zero. It
  validates the equivalent positive-time simplicity target and the exact
  multiplicity-m Hermite cluster coefficient m(m-1)/(8*tau), leaving Xi
  simplicity or finite-truncation endpoint-energy integrability open

Jensen-window PF Newman strict-Laguerre correlation target:
  validates the exact first-correlation representation of
  H_t'^2-H_t*H_t'' and the RH-equivalent Wiener translate-density criterion
  for every 0<t<=1/2. An exact smooth strictly log-concave positive-definite
  convolution has double Fourier zeros, so generic kernel shape cannot replace
  the open Xi-specific zero-free correlation theorem

Jensen-window PF Newman Polymath-15 oscillatory zeta handoff theorem:
  validates an eleven-exponent-pair envelope with exact threshold
  c_*=4911678521/1933561194 and the zeta-jet/remainder composition proving
  exact first-Laguerre positivity for every fixed c>=c_*+epsilon at
  sufficiently large L. It supplies no practical L_epsilon and does not cover
  the scaled strip at or below c_*

Jensen-window PF Newman Polymath-15 cancellation/zero-free wall gate:
  validates the exact weighted-block cancellation frontier, the conditional
  c=2 zeta 1-line handoff, and the fixed c<2 nonpromotion wall. It records the
  precise local exponent deficit without asserting a stronger exponential-sum
  theorem or the open inner Wronskian separation

Jensen-window PF Newman Polymath-15 Gaussian/Legendre duality gate:
  validates the exact finite Gaussian-shift identity and Legendre equivalence,
  including q_*=4800718975/7734244776 equality and the exact c=2 deficit. It
  blocks same-input semigroup repackaging without asserting stronger
  cancellation or inner Wronskian separation

Jensen-window PF Newman Polymath-15 ANTEDB beta-frontier audit:
  validates the exact beta-coordinate objective against pinned ANTEDB commit
  99668603896af86e6cda90ed6755cf3116aab0ac. Six post-2023 pairs and the
  documented one-pass direct-beta envelope retain the exact c_* contact, as do
  twelve finite beta-only van der Corput passes. The audit does not promote raw
  optimizer partitions, infinite transform closure, stronger cancellation, or
  inner Wronskian separation

Jensen-window PF Newman Polymath-15 critical scaled coercivity target:
  validates the corrected finite Riemann-Siegel curvature identities and
  refined C2 handoff on 0<=c<=c_*+o(1), while keeping the phase-sensitive
  coercivity inequality explicitly open and subordinate to the global
  strict-correlation target

Jensen-window PF Newman correlation hierarchy Gaussian-mixture gate:
  validates the exact K_n/F_n heat hierarchy and its universal multiple-root
  contact. Complete monotonicity in v^2 would close strict Fourier positivity,
  but the exact double-exponential Phi tail rejects Gaussian mixtures and
  direct PF-infinity membership; independent 55-digit quadrature also detects
  the first negative log-convexity minor. A tail-compatible spectral-square or
  hierarchy-coercivity theorem remains open

Jensen-window PF Newman theta-summand spectral-square gate:
  validates Phi=(R''-R)/8, R'(0)=-1/2, the all-time identity
  H_t=1/16+D_t[C_t]/8, and the bilateral Gamma/Mellin reconstruction of xi.
  Every finite theta truncation retains A_N>0 and therefore has eventual
  Laguerre tail -2A_N^2/x^6<0; 60-digit witnesses also reject pairwise
  self/cross positivity. One infinite-theta curvature estimate remains

Jensen-window PF Newman Gasper fake-Xi remainder gate:
  validates the exact fake-Xi normalization and Laguerre remainder algebra.
  A pointwise invariant optimizes over every positive scalar, and Arb midpoint
  certificates with explicit tails put its absolute budget above one at x=25
  for t=0 and t=1/2. Exactly, Phi/Psi is nonconstant with finite tail limit
  and cannot be a positive cosh mixture. One-block scalar domination and
  direct Cardon convolution are closed, while sign-aware and multi-block
  routes remain open

Jensen-window PF Newman Gasper residual two-block gate:
  validates exact positivity of Phi-Psi_9, the complete 9/5 residual-tail
  parameter bound, and a beta-quadratic Laguerre identity. Acb Cauchy jets and
  Arb convexity cover the whole tail-admissible interval with a negative first
  Laguerre value at x=50 or x=66. Larger beta has a negative residual tail;
  only signed, coupled multi-block identities remain open

Jensen-window PF Newman classical three-block residual gate:
  validates the full nonnegative 9/5/1 residual obstruction, including the
  two classical real-zero benchmarks and complete Acb/Arb coverage of 64,908
  rational parameter boxes

Jensen-window PF Newman signed universal-factor residual gate:
  validates the signed 9/5/1 standard universal-factor obstruction through an
  exact quartic reduction and a complete 4,094-leaf adaptive interval cover

Jensen-window PF Newman theta/Bessel higher-shift regularization gate:
  validates the exact fixed-index higher-shift expansion and sign threshold,
  then proves at zero frequency that its arithmetic Bessel-transform sum
  diverges; modular grouping or sign-preserving renormalization is required

Jensen-window PF Newman theta cell-renormalization gate:
  validates a normally convergent endpoint cell subtraction, Euler-zeta
  transform assembly, positive zero-frequency blocks, and an absolutely
  convergent coupled Laguerre matrix; its nonzero exp(-5u) block tails prove
  that the construction cannot be deformed termwise to positive Newman time

Jensen-window PF Laguerre scale-mixture gate:
  validates the exact Kummer/Laguerre kernel in every unshifted Jensen moment
  integral. Each fixed scale is hyperbolic, but exact positive-mixture and
  log-concave-density countermodels block integration as a generic preserver.
  Half-integer Gamma mixing is hyperbolic in every degree by Jacobi
  factorization; an Xi-specific zero-preserving connection theorem is open

Jensen-window PF rank-two boundary-family lemma:
  validates an exact all-degree/all-shift model whose Jensen windows factor
  into one simple and one repeated negative root, plus multiplier closure for
  finite pointwise products; the -937/3456 mixture and <-27/125 fractional-
  power guards leave only a discrete counting-measure factorization route

Jensen-window PF multiplier counting-measure target:
  validates the exact sufficient unit-atomic multiplier-product theorem shape,
  its contraction-log and Laplace-kernel reductions, the finite necessary
  zeta sign evidence, two exact continuous-weight failure gates, and the
  separate Mellin-interpolation obstruction; no zeta counting measure is
  constructed

Jensen-window PF Mellin multiplier power-sum obstruction:
  validates eleven Arb-enclosed Phi log moments, nine continuous-product
  power-sum candidates, and six shifted Hankel determinants; three are
  strictly negative, ruling out the natural continuous Mellin product while
  leaving integer-only coefficient equality open

Jensen-window PF negative-lambda first-summand saddle-wall closure:
  validates exact strict saddle concavity, the uniform saddle bracket, and nine
  high-precision cross-checks through k=20000; the paired theorems prove
  L_k^(1)>=1/(4*k^2) for every integer k>=319 and close the adjacent wall

Jensen-window PF negative-lambda first-summand cumulant bridge:
  validates four exact identities and the conditional reduction from
  kappa_3,t(2*log(U))>=-37/(50*t^2), t>=318, to the quantitative dominant
  wall; the paired remainder theorem discharges the uniform special-kernel
  cumulant estimate, while nine high-precision samples remain cross-checks

Jensen-window PF negative-lambda first-summand leading-saddle certificate:
  validates seven symbolic potential derivatives, three sets of 40,740 Arb
  intervals, and three exact ray gates proving caps 13/20, 1/100, and 1/1000
  through fifth order; the open cumulant step is reduced to a seventh-order
  normalized remainder floor -79/1000

Jensen-window PF negative-lambda first-summand paired-remainder certificate:
  validates exact paired standardized moments, 40,736 Arb eighth-derivative
  envelope intervals, rigorous midpoint errors, and explicit two-sided tails;
  proves the floor -79/1000 and strict negative third log-cumulant on all 4,074
  compact blocks covering 0.9264<=u<=5

Jensen-window PF negative-lambda first-summand paired-remainder ray certificate:
  validates the adaptive sqrt(8 log q) window, Gaussian moment comparison,
  and analytic tails proving H_t>=-3/250 on u>=5; composes with the compact
  theorem to close the global cumulant hypothesis for every t>=318

Jensen-window PF heat-flow ratio cone invariance lemma:
  validates a 6-row exact conditional invariance lemma for the ratio cone
  (2*k-1)/(2*k+1)<=x_k<=1 and x_{k+1}>=x_k; its former zeta cone-entry and
  infinite-flow hypotheses are now discharged by the paired certificates

Jensen-window PF heat-flow cone-entry asymptotic target:
  records full ratio-cone entry at lambda=-100 and infinite-dimensional
  forward propagation to lambda=0 as validated; the remaining hard gap is
  downstream all-shape PF/Jensen propagation

Jensen-window PF Phi Taylor cone-entry sign scout:
  validates 4 Phi Taylor coefficient balls and the two local sign combinations
  2*b-a^2<0 and 2*(a^3-3*a*b+3*c)>0 needed by the fixed-k cone-entry route,
  while leaving the uniform-in-k or collared finite theorem open

Jensen-window PF negative-lambda cone-entry prefix scout:
  validates 69 ACB coefficient rows and finite ratio-cone prefix inequalities
  for lambda=-25,-50,-100: lower and upper walls through k<=21 and monotone
  gaps through k<=20, while leaving the all-k tail or finite-collar theorem open

Jensen-window PF negative-lambda cone-entry prefix k30 scout:
  validates 93 ACB coefficient rows and finite ratio-cone prefix inequalities
  for lambda=-25,-50,-100: lower and upper walls through k<=29 and monotone
  gaps through k<=28, while leaving the all-k tail or finite-collar theorem open

Jensen-window PF negative-lambda cone-entry prefix k50 scout:
  validates 153 ACB coefficient rows and finite ratio-cone prefix inequalities
  for lambda=-25,-50,-100: lower and upper walls through k<=49 and monotone
  gaps through k<=48, while leaving the all-k tail or finite-collar theorem open

Jensen-window PF negative-lambda cone-entry prefix k60 scout:
  validates 183 ACB coefficient rows and finite ratio-cone prefix inequalities
  for lambda=-25,-50,-100: lower and upper walls through k<=59 and monotone
  gaps through k<=58, while leaving the all-k tail or finite-collar theorem open

Jensen-window PF negative-lambda cone-entry prefix k80 scout:
  validates 243 ACB coefficient rows and finite ratio-cone prefix inequalities
  for lambda=-25,-50,-100: lower and upper walls through k<=79 and monotone
  gaps through k<=78, while leaving the all-k tail or finite-collar theorem open

Jensen-window PF negative-lambda cone-entry prefix k100 scout:
  validates 303 ACB coefficient rows and finite ratio-cone prefix inequalities
  for lambda=-25,-50,-100: lower and upper walls through k<=99 and monotone
  gaps through k<=98, while leaving the all-k tail or finite-collar theorem open

Jensen-window PF negative-lambda cone-entry prefix k150 scout:
  validates 453 ACB coefficient rows and finite ratio-cone prefix inequalities
  for lambda=-25,-50,-100: lower and upper walls through k<=149 and monotone
  gaps through k<=148, while leaving the all-k tail or finite-collar theorem open

Jensen-window PF negative-lambda cone-entry prefix k200 scout:
  validates 603 ACB coefficient rows and finite ratio-cone prefix inequalities
  for lambda=-25,-50,-100: lower and upper walls through k<=199 and monotone
  gaps through k<=198, while leaving the all-k tail or finite-collar theorem open

Jensen-window PF negative-lambda finite-collar contract:
  validates the finite-collar accounting for the current negative-lambda
  prefix: active depth K=19, first collar x_20, and second collar x_21

Jensen-window PF negative-lambda finite-collar k30 contract:
  validates the finite-collar accounting for the extended k30 prefix:
  active depth K=27, first collar x_28, and second collar x_29

Jensen-window PF negative-lambda finite-collar k50 contract:
  validates the finite-collar accounting for the extended k50 prefix:
  active depth K=47, first collar x_48, and second collar x_49

Jensen-window PF negative-lambda finite-collar k60 contract:
  validates the finite-collar accounting for the extended k60 prefix:
  active depth K=57, first collar x_58, and second collar x_59

Jensen-window PF negative-lambda finite-collar k80 contract:
  validates the finite-collar accounting for the extended k80 prefix:
  active depth K=77, first collar x_78, and second collar x_79

Jensen-window PF negative-lambda finite-collar k100 contract:
  validates the finite-collar accounting for the extended k100 prefix:
  active depth K=97, first collar x_98, and second collar x_99

Jensen-window PF negative-lambda finite-collar k150 contract:
  validates the finite-collar accounting for the extended k150 prefix:
  active depth K=147, first collar x_148, and second collar x_149

Jensen-window PF negative-lambda finite-collar k200 contract:
  validates the finite-collar accounting for the extended k200 prefix:
  active depth K=197, first collar x_198, and second collar x_199

Jensen-window PF negative-lambda tail-barrier scout:
  validates the defect-form tail-barrier diagnostics on the current
  negative-lambda prefix: 63 cone-buffer rows, 60 defect-monotone rows, a
  certified one-third-width buffer, and one rejected scaled-defect shortcut

Jensen-window PF negative-lambda tail-barrier k30 scout:
  validates the defect-form tail-barrier diagnostics on the extended k30
  prefix: 87 cone-buffer rows, 84 defect-monotone rows, a certified
  one-third-width buffer, and one rejected scaled-defect shortcut

Jensen-window PF negative-lambda tail-barrier k50 scout:
  validates the defect-form tail-barrier diagnostics on the extended k50
  prefix: 147 cone-buffer rows and 144 defect-monotone rows; the stronger
  one-third-width buffer holds on 128/147 rows, marking a finite frontier,
  and the scaled-defect shortcut remains rejected

Jensen-window PF negative-lambda tail-barrier k60 scout:
  validates the defect-form tail-barrier diagnostics on the extended k60
  prefix: 177 cone-buffer rows and 174 defect-monotone rows; the stronger
  one-third-width buffer holds on 139/177 rows, sharpening the finite frontier,
  and the scaled-defect shortcut remains rejected

Jensen-window PF negative-lambda tail-barrier k80 scout:
  validates the defect-form tail-barrier diagnostics on the extended k80
  prefix: 237 cone-buffer rows and 234 defect-monotone rows; the stronger
  one-third-width buffer holds on 159/237 rows, sharpening the finite frontier,
  and the scaled-defect shortcut remains rejected

Jensen-window PF negative-lambda tail-barrier k100 scout:
  validates the defect-form tail-barrier diagnostics on the extended k100
  prefix: 297 cone-buffer rows and 294 defect-monotone rows; the stronger
  one-third-width buffer holds on 179/297 rows, sharpening the finite frontier,
  and the scaled-defect shortcut remains rejected

Jensen-window PF negative-lambda tail-barrier k150 scout:
  validates the defect-form tail-barrier diagnostics on the extended k150
  prefix: 447 cone-buffer rows and 444 defect-monotone rows; the stronger
  one-third-width buffer holds on 179/447 rows, sharpening the finite frontier,
  and the scaled-defect shortcut remains rejected

Jensen-window PF negative-lambda tail-barrier k200 scout:
  validates the defect-form tail-barrier diagnostics on the extended k200
  prefix: 597 cone-buffer rows and 594 defect-monotone rows; the stronger
  one-third-width buffer holds on 179/597 rows, sharpening the finite frontier,
  and the scaled-defect shortcut remains rejected

Jensen-window PF negative-lambda scaled-defect frontier k50 scout:
  validates the scaled-defect frontier on the extended k50 prefix: exact cone
  and half-width buffers hold on 147/147 rows, the one-third buffer holds on
  128/147 rows with 19 failures, and scaled-defect nonincrease remains rejected

Jensen-window PF negative-lambda scaled-defect frontier k60 scout:
  validates the scaled-defect frontier on the extended k60 prefix: exact cone
  and half-width buffers hold on 177/177 rows, the one-third buffer holds on
  139/177 rows with 38 failures, and scaled-defect nonincrease remains rejected

Jensen-window PF negative-lambda scaled-defect frontier k80 scout:
  validates the scaled-defect frontier on the extended k80 prefix: exact cone
  and half-width buffers hold on 237/237 rows, the one-third buffer holds on
  159/237 rows with 78 failures, and scaled-defect nonincrease remains rejected

Jensen-window PF negative-lambda scaled-defect frontier k100 scout:
  validates the scaled-defect frontier on the extended k100 prefix: exact cone
  and half-width buffers hold on 297/297 rows, the one-third buffer holds on
  179/297 rows with 118 failures, and scaled-defect nonincrease remains rejected

Jensen-window PF negative-lambda scaled-defect frontier k150 scout:
  validates the scaled-defect frontier on the extended k150 prefix: exact cone
  holds on 447/447 rows, fixed half-width holds on 430/447 rows with first
  failure at lambda=-25 and k=133, the one-third buffer holds on 179/447 rows
  with 268 failures, and scaled-defect nonincrease remains rejected

Jensen-window PF negative-lambda scaled-defect frontier k200 scout:
  validates the scaled-defect frontier on the extended k200 prefix: exact cone
  holds on 597/597 rows, fixed half-width holds on 521/597 rows with first
  failures at lambda=-50 and k=191 and lambda=-25 and k=133, the one-third
  buffer holds on 179/597 rows with 418 failures, and scaled-defect
  nonincrease remains rejected

Jensen-window PF negative-lambda defect-recurrence scout:
  validates the recurrence-route diagnostic on the current negative-lambda
  prefix: 63 buffered rows, 60 defect-monotone rows, and 60 rejections of the
  direct width-preserving recurrence, while leaving the all-k theorem open

Jensen-window PF negative-lambda log-curvature bridge:
  validates the log-curvature reduction for the buffered route: 63 simple
  log-buffer rows, 63 exact defect-buffer rows, 60 curvature-monotone rows,
  and 5 bridge rows, while leaving the bounded-curvature theorem open

Jensen-window PF negative-lambda bounded log-curvature target:
  validates the historical k<=22 compatibility target for the fixed wall
  B_k=-Delta^2 log A_k<=2/(3*(2*k+1)), including 63 finite raw-threshold rows,
  now classified as a diagnostic rather than a live theorem target

Jensen-window PF negative-lambda bounded log-curvature k300 obstruction:
  validates the repaired k300 obstruction to the fixed wall
  C_k=(2*k+1)*B_k<=2/3, with 179/897 rows passing, 718/897 rows failing,
  0 inconclusive rows, and 894/894 scaled-curvature increase rows supporting
  the replacement route only as finite evidence

Jensen-window PF negative-lambda Gaussian curvature matrix:
  validates 7 Gaussian-baseline matrix rows, including 63 positive-deficit
  rows, 63 bounded-deficit rows, and a rejected positive Gaussian
  scale-mixture template for the upper-wall proof

Jensen-window PF negative-lambda signed Gaussian perturbation matrix:
  validates 8 fixed-k signed-perturbation rows, 2 certified Taylor signs, and
  1 finite-depth activation scale while preserving the uniform-remainder gap

Jensen-window PF negative-lambda uniform remainder target:
  validates 8 open-target rows isolating the two-scale theorem needed to
  promote the signed perturbation route, including 2 open analytic
  requirements and 3 leading-only scale diagnostics below the k>=22 tail start

Jensen-window PF negative-lambda Taylor moment budget:
  validates 9 Taylor-moment budget rows and 7 k=22 tail-start samples,
  rejecting low-order Taylor truncation as a finite proof model while
  isolating the required positivity and log-curvature Taylor-tail estimates

Jensen-window PF negative-lambda high-order Taylor scout:
  validates 8 certified even Taylor coefficient rows through c14 and 35
  high-order truncation rows, rejecting higher finite Taylor truncation as a
  replacement for an analytic infinite-tail remainder theorem

Jensen-window PF negative-lambda defect-tail theorem target:
  validates 8 open-target rows stating the all-k defect-tail handoff:
  0<=d_k<=2/(2*k+1) for k>=150 and d_(k+1)<=d_k for k>=149, with 2 live
  noncircular proof routes and 0 ready-to-apply rows

Jensen-window PF negative-lambda half-width tail target:
  validates 9 finite diagnostic rows rejecting the fixed half-width scaled-defect
  route: 430/447 finite half-width rows, 17 half-width failures, 0 live
  half-width proof routes, and a separate monotone defect bridge requirement

Jensen-window PF negative-lambda adaptive scaled-defect target:
  validates 8 open-target rows replacing the rejected fixed half-width route:
  597 exact-cone rows, 76 half-width failures, 2 live noncircular proof routes,
  and a separate monotone defect bridge requirement

Jensen-window PF negative-lambda adaptive envelope matrix:
  validates 7 finite theorem-search rows isolating the monotone-envelope
  route: 594/594 k-increase rows, 398/398 lambda-order rows, 76 half-width
  failures, and no ready-to-apply proof rows

Jensen-window PF negative-lambda adaptive envelope obligations:
  validates 9 exact algebraic obligation rows: 3 exact rows, 3 open
  requirements, 1 rejected route, and no ready-to-apply proof rows

Jensen-window PF negative-lambda raw-moment bridge matrix:
  validates 8 exact/finite theorem-search rows translating the adaptive route
  into raw moment ratios: 597/597 raw-cone rows, 594/594 adaptive corridor
  rows, 76 half-width failures, and no all-k proof claim

Jensen-window PF negative-lambda raw-ratio decrement-corridor scout:
  validates 9 exact/finite theorem-search rows rewriting the raw adaptive
  corridor as a decrement recurrence: 594/594 decrement-corridor rows,
  591/591 theta-k monotonicity rows, 396/396 theta lambda-order rows, and
  2 exact counterexamples blocking raw monotone-decrease shortcuts

Jensen-window PF negative-lambda k300 precision-repair audit:
  validates 7 precision-repair rows extending the raw-ratio decrement route to
  repaired k300 stress: 894/894 decrement-corridor rows, 891/891 theta-k
  monotonicity rows, 596/596 theta lambda-order rows, and no ready-to-apply
  proof rows

Jensen-window PF negative-lambda raw-log decrement bridge:
  validates 8 exact/finite theorem-search rows rewriting the decrement
  corridor as log-ratio inequalities: 894/894 log-corridor rows, 894/894
  log-decrease rows, 2 exact counterexamples, and no ready-to-apply proof rows

Jensen-window PF negative-lambda coefficient-curvature corridor bridge:
  validates 9 exact/finite theorem-search rows rewriting the decrement
  corridor as coefficient-curvature inequalities: 894/894 curvature-corridor
  rows, 894/894 monotone-curvature rows, 2 exact counterexamples, and no
  ready-to-apply proof rows

Jensen-window PF negative-lambda linear curvature-barrier scout:
  validates 8 exact/finite theorem-search rows reducing the nonlinear lower
  curvature barrier to a stronger linear target: 894/894 linear-barrier rows,
  894/894 monotone-curvature rows, 2 exact counterexamples, and no
  ready-to-apply proof rows

Jensen-window PF negative-lambda scaled-curvature monotonicity target:
  validates 10 historical target rows that isolated the replacement theorem
  C_(k+1)>=C_k, with 894/894 repaired k300 finite stress rows; the following
  continuous-bridge gate now discharges that target at lambda=-100

Jensen-window PF negative-lambda scaled-curvature continuous bridge:
  validates 10 exact/interval/analytic rows, 16,074 continuous compact mode
  blocks, 318 repaired prefix gaps, and one exact u>=5 ray; proves
  C_(k+1)>=C_k for every k>=1 at lambda=-100 with 0 open requirements

Jensen-window PF negative-lambda scaled-curvature log-ceiling bridge:
  validates 8 exact/finite rows rewriting C_(k+1)>=C_k as an affine
  delta_k ceiling, with 894/894 scaled-ceiling rows, 894/894
  scaled-log-corridor rows, and 894/894 ceiling-dominance rows

Jensen-window PF negative-lambda relative-Gaussian curvature bridge:
  validates 8 exact/finite rows rewriting C_(k+1)>=C_k as a weighted
  relative-Gaussian four-point inequality, with 897/897 B-positive rows,
  894/894 B-decrease rows, 894/894 C-increase rows, and 598/598
  C-lambda-order rows

Jensen-window PF negative-lambda relative-Gaussian Taylor stencil scout:
  validates 8 theorem-search rows with 3 certified fixed-k leading-sign rows,
  35 finite truncation rows, and only 4 all-positive stencil rows, keeping the
  uniform Taylor-tail remainder theorem open

Jensen-window PF negative-lambda relative-Gaussian stencil remainder obligations:
  validates 9 exact/finite rows decomposing the Taylor-tail error into B,
  companion, and weighted-gap epsilon stencils, with 4 positive baseline rows,
  31 blocked baseline rows, 4 exact stencil rows, and 0 ready-to-apply rows

Jensen-window PF negative-lambda relative-Gaussian pointwise tail budget:
  validates 9 exact/finite rows converting stencil half-margin budgets into
  pointwise log-tail and multiplicative relative-tail tolerances, with 4
  positive baseline rows, 4 budget rows, and 0 ready-to-apply rows

Jensen-window PF negative-lambda relative-Gaussian next-increment stencil stress:
  validates 8 finite theorem-search rows showing 2 tested next-increment rows,
  2 pointwise budget failures, 2 stencil-sign-preserving rows, and 0
  ready-to-apply rows

Jensen-window PF negative-lambda relative-Gaussian degree-16 stencil continuation:
  validates 7 finite theorem-search rows computing the degree-16 continuation,
  with 4 tested continuation rows, 3 stencil-sign-preserving rows, 1
  stencil-sign-failure row, and 0 ready-to-apply rows

Jensen-window PF negative-lambda relative-Gaussian degree-16 collar scan:
  validates 7 finite theorem-search rows scanning integer T=900..2200 at k=22,
  with 1301 scan rows, 1045 continuation-positive rows, 718 half-safety rows,
  pointwise-budget failure on all 1283 baseline-positive rows, and 0
  ready-to-apply rows

Jensen-window PF negative-lambda relative-Gaussian degree-16 real-T collar scout:
  validates 8 finite theorem-search rows certifying a rationalized fixed-k
  degree-16 surrogate collar on the real half-line T>=1156, with 4 positive
  normalizer rows, 3 certified surrogate stencil rows, and 0 ready-to-apply
  rows

Jensen-window PF negative-lambda relative-Gaussian degree-16 Arb real-T collar certificate:
  validates 8 finite theorem-search rows upgrading the fixed-k degree-16
  real-T surrogate collar to Arb coefficient-ratio balls, with 4 positive
  normalizer rows, 3 Bernstein-certified stencil rows, and 0 ready-to-apply
  rows

Jensen-window PF negative-lambda relative-Gaussian degree-40 Arb collar ladder stress:
  validates 8 finite theorem-search rows stress-testing the same T>=1156
  collar for every even finite surrogate degree 16 through 40, with 13 degree
  levels, zero Bernstein failures, and 0 ready-to-apply rows

Jensen-window PF negative-lambda relative-Gaussian degree-40 residual tail budget:
  validates 8 exact/finite theorem-search rows converting the degree-40 collar
  margins into sufficient fixed-k residual targets, with 5 budget
  inequalities, 4 finite degree-42..80 tail-profile rows, and 0
  ready-to-apply rows

Jensen-window PF negative-lambda relative-Gaussian formal-tail obstruction scout:
  validates 8 finite obstruction rows extending the formal residual profile
  through degree 240, rejecting monotone/geometric infinite formal-tail
  summation after the j=103 turnaround, and leaving the actual asymptotic
  remainder theorem open

Jensen-window PF negative-lambda relative-Gaussian asymptotic remainder target:
  validates 6 exact theorem-search rows calibrating sufficient analytic
  remainder constants: a 1000x first-omitted-term theorem would fit the
  fixed-k degree-40 half-safety budget, while an optimized finite window
  target keeps the actual remainder theorem open

Jensen-window PF negative-lambda relative-Gaussian actual endpoint remainder scout:
  validates 6 floating endpoint theorem-search rows comparing the actual
  T=1156 relative-Gaussian quadrature residuals against first omitted terms
  and degree-40 budgets, with 4 endpoint rows, 5 quadrature orders, and 0
  ready-to-apply rows

Jensen-window PF negative-lambda relative-Gaussian cancellation-reduced remainder grid scout:
  validates 6 floating cancellation-reduced theorem-search rows subtracting
  the degree-40 polynomial inside the Gamma expectation across 20 grid rows,
  5 T values, 4 indices, and 4 quadrature orders, with 0 ready-to-apply rows

Jensen-window PF negative-lambda relative-Gaussian intervalization target:
  validates 6 open numerical-certification rows translating the
  cancellation-reduced finite-grid slack into 8 obligations, 5 open
  requirements, and 0 ready-to-apply rows

Jensen-window PF negative-lambda relative-Gaussian Phi tail bound scout:
  validates 6 analytic padded-range tail rows bounding the n>30 Phi, Phi',
  and Phi(0) tails below 1e-1000 after diagnostic normalization, with 2
  conditional requirements and 0 ready-to-apply rows

Jensen-window PF negative-lambda relative-Gaussian node-c0 range certificate:
  validates 5 exact/Arb node-c0 rows, with 16 Laguerre Gershgorin bound rows,
  2 certified side conditions for the Phi-tail scout, and 0 ready-to-apply
  rows

Jensen-window PF negative-lambda relative-Gaussian quadrature ladder scout:
  validates 5 high-order floating quadrature-ladder rows on the worst
  cancellation-reduced grid row, with 7 ladder rows through reference order
  320 and 0 ready-to-apply rows

Jensen-window PF negative-lambda relative-Gaussian quadrature-remainder route matrix:
  validates 7 route-matrix rows converting the N=320 Laguerre remainder
  prefactor into derivative-supremum and interval-integration obligations,
  with derivative order 640 and 0 ready-to-apply rows

Jensen-window PF negative-lambda relative-Gaussian worst-row far-tail split certificate:
  validates 7 Arb far-tail split rows for T=10000, F_21, certifying y>=200
  below the quadrature cap while leaving compact integration on 0<=y<=200 open

Jensen-window PF negative-lambda relative-Gaussian worst-row compact-interval integration scout:
  validates 7 compact-interval diagnostic rows for T=10000, F_21, with 6 raw
  Arb panels on 0<=y<=200; plain interval Riemann is rejected and 0 rows are
  ready-to-apply

Jensen-window PF negative-lambda relative-Gaussian worst-row Chebyshev panel-moment scout:
  validates 7 high-precision floating Chebyshev panel-moment rows, with 5
  degrees, 4 Cauchy pairs, 3 cap-safe pairs, and 0 ready-to-apply rows

Jensen-window PF negative-lambda relative-Gaussian worst-row Arb Chebyshev interpolant-moment scout:
  validates 7 Arb Chebyshev interpolant-moment rows, with 4 degrees, 3 Cauchy
  pairs, 3 cap-safe pairs, and 0 ready-to-apply rows

Jensen-window PF negative-lambda relative-Gaussian worst-row interpolation-remainder route matrix:
  validates 8 route-matrix rows, with 6 panel masses, 20 Bernstein budgets,
  16 minimal-degree rows, and 0 ready-to-apply rows

Jensen-window PF negative-lambda relative-Gaussian endpoint parity-repair matrix:
  validates 7 endpoint parity-repair rows, with 8 odd Taylor rows through
  order 15 and 0 ready-to-apply rows

Jensen-window PF negative-lambda relative-Gaussian endpoint x-panel route matrix:
  validates 7 endpoint x-panel route rows, with x interval 0<=x<=0.01,
  18 Bernstein budgets, and 0 ready-to-apply rows

Jensen-window PF negative-lambda relative-Gaussian endpoint x-moment Taylor certificate:
  validates 7 first-panel certificate rows, with 65 exact transformed-moment
  rows, one certified finite-core endpoint panel, 5 open compact panels, and
  0 ready-to-apply rows

Jensen-window PF negative-lambda relative-Gaussian worst-row compact x-moment Taylor certificate:
  validates 7 compact certificate rows, with 129 exact transformed-moment
  rows, one certified finite-core interval 0<=y<=200, 0 open compact panels,
  and 0 ready-to-apply rows

Jensen-window PF negative-lambda relative-Gaussian worst-row full-expectation certificate:
  validates 8 complete one-row composition rows, with 3 rigorous sources,
  2 complete true first-omitted ratios below one, 0 open worst-row integral
  sources, no final quadrature dependence, and 0 ready-to-apply rows

Jensen-window PF negative-lambda relative-Gaussian first-omitted denominator certificate:
  validates 6 Arb denominator-side rows, with 20 positive first-omitted
  denominator rows, 2 ratio-cap conversion rows, and 0 ready-to-apply rows

Jensen-window PF negative-lambda relative-Gaussian coefficient-core propagation certificate:
  validates 6 coefficient-core propagation rows, with 22 Arb coefficient
  rows, 20 exact-Gamma propagation rows, 2 intervalization handoff rows, and
  0 ready-to-apply rows

Jensen-window PF negative-lambda relative-Gaussian worst-row Laguerre root-bracket certificate:
  validates 6 worst-row node-bracket rows, with 320 Arb sign-changing root
  brackets for L_320^(41/2), 30 zero floating SciPy weights recorded as an
  underflow boundary, 2 intervalization handoff rows, and 0 ready-to-apply rows

Jensen-window PF negative-lambda relative-Gaussian worst-row Christoffel-weight midpoint scout:
  validates 6 worst-row weight-midpoint rows, with 320 positive Arb midpoint
  weights, 30 repaired SciPy floating underflows, 320 direct interval
  denominator obstructions, and 0 ready-to-apply rows

Jensen-window PF negative-lambda relative-Gaussian worst-row Christoffel-weight interval certificate:
  validates 6 worst-row weight-interval rows, with 320 positive Christoffel
  weight intervals, 0 Taylor denominator obstructions, 30 repaired SciPy
  floating underflows, a weight-sum interval containing Gamma(43/2), and 0
  ready-to-apply rows

Jensen-window PF negative-lambda relative-Gaussian worst-row finite-part weighted-sum interval certificate:
  validates 6 worst-row finite-part weighted-sum rows, with 320 refined
  nodes, 320 interval weights, 2 below-one first-omitted ratios, and 0
  ready-to-apply rows

Jensen-window PF negative-lambda relative-Gaussian worst-row finite-plus-tail budget certificate:
  validates 6 worst-row finite-plus-tail budget rows, with 2 composed
  first-omitted ratios below one, 3 certified tail sources, and 0
  ready-to-apply rows

Jensen-window PF negative-lambda raw-moment obstruction matrix:
  validates 7 exact countermodel-gate rows with 3 positive two-atom Stieltjes
  witnesses blocking generic raw moment/log-convexity promotion into the upper
  raw wall or adaptive corridor

Jensen-window PF negative-lambda zeta-specific raw-corridor target:
  validates 9 historical target rows naming the actual zeta raw-wall/corridor
  theorem, with 2 live routes and 2 rejected shortcuts; the following exact
  composition gate now discharges that target at lambda=-100

Jensen-window PF negative-lambda -100 raw-corridor certificate:
  validates 6 exact theorem-composition rows combining full ratio-cone entry,
  scaled-curvature growth, and the linear-to-nonlinear calculus lemma into the
  complete raw-ratio decrement corridor for every k>=1 at lambda=-100

Jensen-window PF negative-lambda -100 adaptive-defect certificate:
  validates 8 exact theorem-composition rows; proves the all-k defect cone,
  decreasing defect, exact scaled-defect cone, and increasing scaled defect at
  lambda=-100, while leaving the stronger simultaneous three-lambda statement
  unproved

Jensen-window PF monotone contraction stress:
  validates 2875 finite Arb-classified zeta-window rows across degrees d=3..12
  and the k<=64 cache, all satisfying adjacent log-concavity and increasing
  ratio contractions with no zero-containing required positivity gaps

Jensen-window PF state-space sign-lift obstruction scout:
  validates 3 symbolic rows and 735 derived mu_2 rows showing that an
  absolute-value sign-state cover of raw Motzkin paths overshoots the actual
  reciprocal coefficient by 2*kappa_1>0; this rejects only the cheap cover and
  leaves genuinely modified state-space doubled models open

Jensen-window PF Cauchy-Binet low-degree scout:
  validates a 15-row symbolic ratio/Bernstein scout for the live Cauchy-Binet
  ansatz; all selected degree-2/3/4 hard formulas are certified by adjacent
  log-concavity, but 0 kernel identities are found

Jensen-window PF log-concavity frontier scout:
  validates a 14-row larger-contiguous-minor frontier scout; adjacent
  log-concavity first loses the simple Bernstein certificate at degree 3 size
  6 and degree 4 size 5, while the exact rational countermodel first becomes
  negative at degree 3 size 8 and degree 4 size 6

Jensen-window PF ratio-condition scout:
  validates a 7-row ratio-condition audit, rejecting adjacent log-concavity,
  decreasing ratio contractions, second-order log-concavity, and selected
  low-degree Bernstein positivity by exact countermodel, and contraction
  log-concavity by constructed positive extension

Jensen-window PF contraction-log-concavity scout:
  validates a positive rational extension satisfying x2^2 >= x1*x3 while the
  degree 3 size 8 and degree 4 size 6 contiguous Jensen-window minors remain
  negative

sign-regularity theorem fit matrix:
  validates that ASW/Edrei, Schur-positive, Schoenberg, Jensen, Hankel, and
  signed-Hankel route notes preserve fit/misfit boundaries and kill gates

positive Schur-specialization target note:
  validates that the coefficient-PF Schur target states the exact
  Jacobi-Trudi equivalence, noncircular sufficient constructions, and
  rejection gates for finite or circular proof attempts

Edrei-log sign diagnostics:
  validates 320 finite Edrei-log sign diagnostics

Edrei power-Hankel diagnostics:
  validates 4,205 finite Edrei power-Hankel diagnostics

Edrei midpoint frontier non-promotion guard:
  validates 5 non-rigorous midpoint frontier scouts

Edrei power-Hankel boundary repair manifest:
  validates 2 retired inconclusive blocker rows and 3 repaired positives

Edrei moment-recurrence scout manifest:
  validates 1 positive Arb recurrence scout and 1 inconclusive frontier scout
```

## Boundary

Passing this runner means the advertised finite ledgers and proof-safety
guards are internally reproducible, and that the exact Toeplitz-to-skew-Schur
reindexing used by the coefficient-PF route still passes its algebraic check.
It also checks for obvious forbidden result-language drift in the markdown
notes and enforces explicit artifact classification metadata. It does not
convert finite certificates, moment diagnostics, midpoint scouts,
counterexamples, reindexing identities, dependency graphs, claim
classification, clean wording, clean metadata, or intact file references into
an all-order bridge theorem.

The remaining theorem burden is narrower but still open:

```text
contiguous and arbitrary-column signed-Hankel layers through fixed order nine
are closed on the heat segment from -100 to zero. The arbitrary-order heat
algebra and induction give an exact conditional endpoint-to-heat equivalence.
After normalizing h_k=A_k(-100)/A_0(-100), however, the proposed static
rectangle antecedent fails: four required order-ten rectangles are strictly
negative. The endpoint sequence also fails a full PF_3 inequality outside
the deep cone, the zero-reset moving-tail translation fails on shallow
shapes, and strict full unweighted Schur positivity has an exact cubic Jensen
counterexample. Construct a noncircular, weaker Xi/Phi-specific kernel,
determinant, or direct Jensen mechanism that tolerates those failures and
still proves Laguerre-Polya membership or all-degree/all-shift Jensen
hyperbolicity. On
the direct Newman route, retain theta modular cancellation before the
spectral transform and establish positivity of the coupled mixed-term matrix
for every 0<t<=1/2. The universal all-real-lambda monotone-wall formulation
is false at lambda=-1156.
```

## Useful Options

For quick checks that skip gates marked slow:

```text
python work/rh_compute/scripts/check_core_proof_programme_gates.py --skip-slow
```

The default per-gate timeout is `600` seconds. This accommodates observed
runtime variance in the 320-node finite-part interval certificate, which can
exceed the former 300-second default without any change in its result.

Current quick result:

```text
latest completed non-slow umbrella coverage: 328/328 core proof-programme gates
the current runner has 349 gates: 333 non-slow and 16 slow
```

The current 324-claim ledger expectation passes. The latest non-slow umbrella
run completed `328/328` gates in `969.3` seconds, including the later
Polymath-15 and Edrei gates, the Xi strict-monotonicity rejection, and the
sixteen order-ten completion gates. Five added order-eleven gates pass focused
validation; a fresh `333/333` umbrella timing is pending. The latest full umbrella run before
order-nine integration completed all former `301/301` gates in one execution,
and the fifteen order-nine/all-order gates have passed focused validation,
including both newly marked slow checks. Earlier
reduction and finite-prefix artifacts retain their historical local
boundaries; their open handoffs are supplied by separately validated
downstream entry and uniform-heat certificates.

For machine-readable output:

```text
python work/rh_compute/scripts/check_core_proof_programme_gates.py --json
```
