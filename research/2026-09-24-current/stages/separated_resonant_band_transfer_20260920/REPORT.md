# Separated resonant band: current-to-response replacement

20 September 2026. TWO HANDS NETWORK LTD. Private local RH research.

## Outcome

The surviving third-mode current has a new directly proved simplification.
On the deterministic product band

    |n1*n2*n3 - (N+r_c)^2| <= (N+r_c) floor(d^(1/4)),
    d=ceil(log(N+1)),

retain only the already cube-root-separated triples. Their complete
five-carrier assigned current can be replaced by their signed response:

    J_band = T_band + O_I(d^(3/4)) = T_band + o_I(log N).

For the narrower half-width N+r_c the payment is O_I(sqrt(d)). These
bands are nested alternatives, not two additive deletions. Both contain
the entire moving |P-(N+r)^2|<= (N+r)/2 region on the original window.

The proof retains the exact Euler rotation, original common inverse,
physical-height Jacobian and both sharp endpoints. It uses a ternary-
divisor short-interval upper bound solely for an absolute error payment.
The signed coefficient sequence is not assumed multiplicative.

This removes the selected inverse-derivative expression from the remaining
normal form. It does NOT remove the response or prove its sign. Outside
products and the surviving first-mode channel are unchanged. The original
0.0011250999 signed logarithmic target and RH remain open.

## Why the sign is still a separate problem

A fixed actual-profile control used N=101^12 and the exact product P=N^2.
There are seven distinct power-factorizations using powers101^2 through
101^11, with exponents summing24. All are cube-root-separated, coprime
to30 and inside the tested d^2-to-M core. At x=1/10,1,4, their coefficient
products include both signs, even though P, the carrier phase and the
companion multiplier are shared. Ordered multiplicities were retained.

This is not a full-source sign calculation or a claim about their net
asymptotic mass. It rules out assigning a favourable sign to every term
merely because the bare resonant cosine is positive.

## Checks

- 19 primary exact check families, including all five carrier identities.
- 600 independently implemented exact polynomial integration cases,
  with nonzero endpoint and nonzero density-derivative rejection guards.
- 63 overlapping primary Arb fields at112 and176bits.
- Three additional profile integrals from its defining measure, not the
  inherited positive-series implementation, at each precision.
- Twelve finite full ternary-divisor band masses at each precision.
- 15 further overlapping interval fields from those checks.

The sparse11-term source-shaped checks use actual height phases, actual
profile coefficients and the exact2/3/5 Euler product. They are not an
evaluation of the full retained RH source. The two interval precisions
share FLINT/Arb; the analytic theorem is not machine-formalized. The
Shiu theorem is a cited analytic input, not verified by finite sampling.

Two failed test executions are retained: an empty sparse control band and
an exclusive-output filename collision. Neither is counted as a pass;
no mathematical counterexample was repaired. See ITERATIONS.md.

## Next target

Use J_reduced=J_sep-J_band+T_band. Seek the joint signed bound for this
response, the outside-product third-mode current and the surviving first
mode. Keep shared product coefficients collected before assigning signs.
No endpoint-uniform claim or favourable numerical constant follows from
the present O_I bound.

The closure seal and clean replay are recorded externally beside the
package. Historical experiments and research packages remain unchanged.
No provider execution, Forge work, Git push or external publication.
