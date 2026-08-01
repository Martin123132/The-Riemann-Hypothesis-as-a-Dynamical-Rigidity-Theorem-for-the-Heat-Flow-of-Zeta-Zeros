#!/usr/bin/env python3
"""Build the planar incidence and signed anti-diagonal handoff."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_mertens_planar_incidence_antidiagonal_handoff.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_mertens_planar_incidence_antidiagonal_handoff.md"
)


def row(
    row_id: str,
    role: str,
    status: str,
    statement: str,
    proof_boundary: str,
) -> dict[str, str]:
    return {
        "id": row_id,
        "role": role,
        "status": status,
        "statement": statement,
        "proof_boundary": proof_boundary,
    }


def build_rows() -> list[dict[str, str]]:
    return [
        row(
            "peia_01_parent_energy",
            "exact_equivalence",
            "available_exact",
            "Corollary 11.22Z.6 gives E=D+C on the edge lift, "
            "with D=O_epsilon(K^(1+epsilon)) and one open signed "
            "off-diagonal edge-pair form C.",
            "This inherits, but does not prove, the edge-pair estimate.",
        ),
        row(
            "peia_02_tail_kernel",
            "exact_definition",
            "available_exact",
            "T_K(X,G)=sum_(x>=X,H>=G)nu_(x,H), so "
            "G_K((p,q),(r,s))=T_K(max(p,r),max(q-p,s-r)).",
            "The curvature measure is nu=|Delta Wtilde| only on the "
            "strong sufficient energy route.",
        ),
        row(
            "peia_03_threshold_graph",
            "exact_definition",
            "available_exact",
            "The threshold graph Gamma_(x,H) has vertices "
            "K+1,...,4K-1 and edges p<q with p<=x and q-p<=H.",
            "The orientation is retained only to encode base and shift.",
        ),
        row(
            "peia_04_adjacency_quadratic",
            "exact_identity",
            "available_exact",
            "If A_t is the symmetric adjacency matrix of Gamma_t, then "
            "S_t=(1/2)<mu,A_t mu>.",
            "A_t is indefinite and this identity supplies no sign bound.",
        ),
        row(
            "peia_05_threshold_diagonal",
            "exact_definition",
            "available_exact",
            "Q_t=sum_({u,v} in Gamma_t)mu(u)^2mu(v)^2; "
            "D=sum_t nu_t Q_t.",
            "Q_t is the equal-edge contribution, not the full energy.",
        ),
        row(
            "peia_06_neighbor_field",
            "exact_definition",
            "available_exact",
            "R_t(v)=sum_(u adjacent to v)mu(u) and "
            "Q_t(v)=sum_(u adjacent to v)mu(u)^2.",
            "The neighbor intervals depend on both x and H.",
        ),
        row(
            "peia_07_full_collision",
            "exact_identity",
            "available_exact",
            "The ordered one-vertex collision at threshold t is "
            "C_3(t)=sum_v mu(v)^2(R_t(v)^2-Q_t(v)).",
            "This includes same-base, same-tip, and chain collisions.",
        ),
        row(
            "peia_08_sorted_triple",
            "exact_identity",
            "available_exact",
            "After summing nu, C_3 is the explicit sorted-triple sum "
            "in (PEIA.4), with the three tail weights "
            "T(a,c-a), T(b,max(b-a,c-b)), and T(b,c-a).",
            "The repeated endpoint introduces mu(v)^2; no Mobius "
            "cancellation is inferred.",
        ),
        row(
            "peia_09_sorted_quadruple",
            "exact_identity",
            "available_exact",
            "The four-distinct-endpoint contribution C_4 is the sorted "
            "quadruple sum in (PEIA.5), one term for each perfect "
            "matching of a<b<c<d.",
            "This is a weighted four-point Mobius form.",
        ),
        row(
            "peia_10_endpoint_expansion",
            "exact_identity",
            "available_exact",
            "The complete signed edge off-diagonal is C=C_3+C_4.",
            "Neither stratum is asserted to be small separately.",
        ),
        row(
            "peia_11_cancellation_witness",
            "countermodel_gate",
            "guard_validated",
            "For K=4, x=11, H=10, and signs +1 on 5..12 and -1 on "
            "13..15, there are 49 edges, S=7, Q=49, C=0, "
            "C_3=210, and C_4=-210.",
            "This bounded synthetic sequence proves exact cancellation "
            "between strata; it is not a statement about Mobius.",
        ),
        row(
            "peia_12_stratum_separation_guard",
            "proof_guard",
            "guard_active",
            "Separate bounds for C_3 and C_4 are sufficient but not "
            "equivalent to the signed C bound.",
            "The cancellation witness forbids presenting separated "
            "stratum estimates as a lossless reformulation.",
        ),
        row(
            "peia_13_random_cut",
            "exact_definition",
            "available_exact",
            "For independent fair vertex colors chi, J_t(chi) is the "
            "sum of c_e over threshold edges whose endpoints have "
            "different colors.",
            "The coloring is auxiliary and does not change the Mobius "
            "labels.",
        ),
        row(
            "peia_14_cut_probabilities",
            "exact_identity",
            "available_exact",
            "A fixed edge is cut with probability 1/2, while any two "
            "distinct edges are both cut with probability 1/4, whether "
            "they collide or are disjoint.",
            "The second probability is the reason both signed strata "
            "remain joined.",
        ),
        row(
            "peia_15_cut_second_moment",
            "exact_identity",
            "available_exact",
            "At every threshold, 4 E_chi[J_t(chi)^2]=S_t^2+Q_t.",
            "This is an identity, not an upper bound for S_t.",
        ),
        row(
            "peia_16_integrated_cut_identity",
            "exact_identity",
            "available_exact",
            "For J_K(chi)=sum_t nu_t J_t(chi)^2, "
            "E=4 E_chi[J_K(chi)]-D.",
            "The already-closed D is the only correction.",
        ),
        row(
            "peia_17_cut_bilinear_form",
            "exact_identity",
            "available_exact",
            "With f_chi=mu 1_(chi=0) and g_chi=mu 1_(chi=1), "
            "J_t=B_t(f_chi,g_chi)+B_t(g_chi,f_chi)"
            "=2<f_chi,A_t g_chi>.",
            "The random support masks destroy multiplicative structure "
            "unless a theorem explicitly tolerates them.",
        ),
        row(
            "peia_18_cut_target_equivalence",
            "exact_equivalence",
            "available_exact",
            "The energy target is equivalent to an "
            "O_epsilon(K^(1+epsilon)) bound for the average cut "
            "bilinear energy.",
            "The equivalence uses the closed diagonal and provides no "
            "automatic gain.",
        ),
        row(
            "peia_19_uniform_coloring_guard",
            "proof_guard",
            "guard_active",
            "A bound uniform over all colorings is sufficient but is "
            "strictly stronger than the required average.",
            "Do not replace the exact average by a worst-case coloring "
            "without recording the loss.",
        ),
        row(
            "peia_20_biased_sign_family",
            "exact_definition",
            "available_exact",
            "Let eta_v be independent signs with E eta_v=m and put "
            "z=m^2; P_K(z) is the expected nu-weighted square after "
            "multiplying every edge by eta_u eta_v.",
            "Only 0<=z<=1 has a direct biased-sign interpretation.",
        ),
        row(
            "peia_21_chaos_polynomial",
            "exact_identity",
            "available_exact",
            "P_K(z)=D+z C_3+z^2 C_4.",
            "The original energy is the boundary value P_K(1).",
        ),
        row(
            "peia_22_randomization_no_gain",
            "proof_guard",
            "guard_active",
            "Mean-zero randomization gives P_K(0)=D, but no positivity "
            "principle controls P_K(1) from interior biased values.",
            "Randomization isolates the chaos levels and does not prove "
            "the deterministic boundary estimate.",
        ),
        row(
            "peia_23_partial_fourier_data",
            "exact_definition",
            "available_exact",
            "F_x(theta)=sum_(K<p<=x)mu(p)e(p theta), "
            "F_*(theta)=sum_(K<q<4K)mu(q)e(q theta), and "
            "K_H(theta)=sum_(h=1)^H e(h theta).",
            "Here e(theta)=exp(2 pi i theta).",
        ),
        row(
            "peia_24_partial_prefix_fourier",
            "exact_identity",
            "available_exact",
            "S_K(x,H)=integral_0^1 F_x(theta) "
            "conj(F_*(theta)) K_H(theta)dtheta.",
            "This keeps the varying base prefix and is bilinear rather "
            "than a full-prefix modulus square.",
        ),
        row(
            "peia_25_full_prefix_specialization",
            "exact_identity",
            "available_exact",
            "At x=4K-2, the real form reduces to the inherited "
            "|F_K|^2 odd-correlation identity.",
            "The varying-x problem is not removed by this specialization.",
        ),
        row(
            "peia_26_current_interior_index",
            "exact_definition",
            "available_exact",
            "For the current-current interior, i>=1, H>=1, "
            "i+H<=K-2, and s=2i+H ranges from 3 to 2K-5.",
            "Transition and future cells lie outside this index set.",
        ),
        row(
            "peia_27_antidiagonal_curvature",
            "exact_identity",
            "available_exact",
            "Delta Wtilde_K(K+i,H)=kappa_(K,R)(2i+H), where kappa_s "
            "is the signed double integral of D_R in (PEIA.13).",
            "The energy route replaces kappa_s by |kappa_s|; the direct "
            "route retains its sign.",
        ),
        row(
            "peia_28_smoothed_odd_spectrum",
            "exact_identity",
            "available_exact",
            "kappa_s=(4 pi^2/N^3)sum_(r<=R)sigma_r "
            "cos(q_r(s+3/2)pi/N), with q_r=2r-1 and the two sinc "
            "factors sigma_r in (PEIA.14).",
            "This is an exact finite odd-frequency expansion.",
        ),
        row(
            "peia_29_curvature_parseval",
            "exact_identity",
            "available_exact",
            "sum_(s=0)^(2N-1)kappa_s^2="
            "(16 pi^4/N^5)sum_(r<=R)sigma_r^2"
            "<=16 pi^4 R/N^5.",
            "The identity is over a complete 2N-point period.",
        ),
        row(
            "peia_30_line_aggregate",
            "exact_definition",
            "available_exact",
            "A_K(s)=sum_(2i+H=s,i+H<=K-2)S_K(K+i,H).",
            "This sums prefixes along one current-interior "
            "anti-diagonal.",
        ),
        row(
            "peia_31_edge_multiplicity",
            "exact_identity",
            "available_exact",
            "A_K(s)=sum_(a,g)mu(K+a)mu(K+a+g)N_s(a,g), with "
            "N_s(a,g)=[floor((s-g)/2)-max(a,s-K+2)+1]_+.",
            "The edge range is a>=1, g>=1, a+g<=K-2.",
        ),
        row(
            "peia_32_signed_interior_factorization",
            "exact_identity",
            "available_exact",
            "The direct current-current interior contribution is "
            "O_CC,int=sum_s kappa_s A_K(s).",
            "This is one component of O_(alpha,K), not the complete "
            "transition/future decomposition.",
        ),
        row(
            "peia_33_projection_coefficients",
            "exact_definition",
            "available_exact",
            "Y_(K,r)=sum_(s=3)^(2K-5)A_K(s)"
            "cos(q_r(s+3/2)pi/N).",
            "Only the first R odd frequencies enter the actual curvature.",
        ),
        row(
            "peia_34_spectral_projection",
            "exact_identity",
            "available_exact",
            "O_CC,int=(4 pi^2/N^3)sum_(r<=R)sigma_r Y_(K,r).",
            "No absolute value over r or s has been taken.",
        ),
        row(
            "peia_35_projection_cauchy",
            "exact_bound",
            "available_exact",
            "|O_CC,int|^2<=(16 pi^4 R/N^6)"
            "sum_(r<=R)|Y_(K,r)|^2.",
            "This Cauchy step is over the retained odd frequencies only.",
        ),
        row(
            "peia_36_projection_sufficient_gate",
            "conditional_target",
            "open",
            "For every epsilon>0, "
            "sum_(r<=R)|Y_(K,r)|^2="
            "O_epsilon(N^6 R^(-1)K^epsilon) would make "
            "O_CC,int subpower on dyadic K.",
            "This target covers only the signed current-current interior "
            "and is not proved.",
        ),
        row(
            "peia_37_bessel_fallback",
            "exact_bound",
            "available_exact",
            "Odd-frequency Bessel gives "
            "sum_(r<=R)|Y_(K,r)|^2<=N sum_s|A_K(s)|^2.",
            "A separate arithmetic line-energy estimate is still needed.",
        ),
        row(
            "peia_38_bessel_power_guard",
            "proof_guard",
            "guard_active",
            "Geometry-only bounds for A_K(s) do not reach the projection "
            "target; the all-one sequence retains polynomial excess.",
            "Bessel orthogonality is an interface, not the missing Mobius "
            "cancellation.",
        ),
        row(
            "peia_39_energy_vs_signed_route",
            "proof_guard",
            "guard_active",
            "The cut bilinearization is exactly equivalent to the strong "
            "energy route, whereas the anti-diagonal projection acts "
            "before absolute mixed curvature and is genuinely weaker.",
            "The two theorem targets must not be conflated.",
        ),
        row(
            "peia_40_remaining_blocks",
            "proof_guard",
            "guard_active",
            "The transition, current-future, and future-future cells from "
            "Corollary 11.22Z.5 remain separate open components on the "
            "direct signed route.",
            "Controlling O_CC,int alone does not prove the complete "
            "off-diagonal criterion.",
        ),
        row(
            "peia_40a_recent_short_interval_guard",
            "literature_guard",
            "guard_active",
            "Menon's 2026 almost-all short-interval Fourier and averaged "
            "Chowla bounds give logarithmic-scale decay after averaging "
            "over base points or shifts.",
            "They do not supply the anchored quadratic odd-frequency "
            "projection bound for Y_(K,r), nor the transition/future "
            "estimates required here.",
        ),
        row(
            "peia_41_finite_validation",
            "finite_validation",
            "validated_finite",
            "Independent finite checks reproduce the collision and "
            "matching formulas, exact cancellation witness, cut moments, "
            "biased chaos polynomial, partial Fourier identity, "
            "anti-diagonal multiplicities, spectrum, and Parseval law.",
            "Finite identities do not prove either asymptotic target.",
        ),
        row(
            "peia_42_proof_boundary",
            "proof_guard",
            "guard_active",
            "The incidence expansion, random-cut equivalence, and signed "
            "current-interior anti-diagonal projection are exact.",
            "No projection estimate, remaining-block estimate, full "
            "energy bound, RH, PF-infinity, Lambda<=0, or Clay-prize "
            "result is proved.",
        ),
    ]


NOTE_BODY = """# Jensen-Window PF Mertens Planar Incidence/Anti-Diagonal Handoff

Date: 2026-07-24

Status: exact incidence, cut-bilinear, and signed current-interior
anti-diagonal reductions with one open projection gate and open
transition/future blocks. This is not a proof of RH, PF-infinity,
or `Lambda <= 0`.

```text
work/rh_compute/results/
jensen_window_pf_mertens_planar_incidence_antidiagonal_handoff.json
python work/rh_compute/scripts/
check_jensen_window_pf_mertens_planar_incidence_antidiagonal_handoff.py
```

## Complete Endpoint Incidence Expansion

Retain the edge set, labels, threshold features, and positive curvature
measure from Corollary 11.22Z.6. Write

```text
T_K(X,G):=sum_(x>=X,H>=G)nu_(x,H).             (PEIA.1)
```

Thus the edge Gram kernel is

```text
G_K((p,q),(r,s))
 =T_K(max(p,r),max(q-p,s-r)).                  (PEIA.2)
```

At a fixed threshold let `Gamma_t` be the simple graph whose edges are
the active pairs, and put

```text
R_t(v):=sum_(u adjacent to v)mu(u),
Q_t(v):=sum_(u adjacent to v)mu(u)^2.
```

The full ordered collision contribution, including same-base,
same-tip, and chain collisions, is

```text
C_3(t)
 =sum_v mu(v)^2[R_t(v)^2-Q_t(v)].             (PEIA.3)
```

After summing over the curvature measure, sort the three endpoints as
`a<b<c`. The complete collision stratum is

```text
C_3=2 sum_(a<b<c) [
 mu(a)^2mu(b)mu(c) T_K(a,c-a)
 +mu(a)mu(b)^2mu(c) T_K(b,max(b-a,c-b))
 +mu(a)mu(b)mu(c)^2 T_K(b,c-a)
].                                            (PEIA.4)
```

For four distinct sorted endpoints `a<b<c<d`, the three perfect
matchings give

```text
C_4=2 sum_(a<b<c<d) mu(a)mu(b)mu(c)mu(d) [
 T_K(c,max(b-a,d-c))
 +T_K(b,max(c-a,d-b))
 +T_K(b,d-a)
].                                            (PEIA.5)
```

Hence the complete signed edge-pair form is exactly `C=C_3+C_4`.
This endpoint expansion is lossless, but estimating the two lines
separately is not.

An exact bounded-sign witness makes that warning concrete. Take
`K=4`, threshold `(x,H)=(11,10)`, and coefficients

```text
c_v=+1 for 5<=v<=12,
c_v=-1 for 13<=v<=15.
```

Then

```text
edge count=49, S=7, Q=49,
C=S^2-Q=0, C_3=210, C_4=-210.                 (PEIA.6)
```

Thus collision and disjoint strata can cancel perfectly. This is a
synthetic proof-safety example, not a Mobius estimate.

## Exact Random-Cut Bilinearization

Color the vertices independently and fairly by `chi_v in {0,1}` and
let

```text
J_t(chi)
 :=sum_({u,v} in Gamma_t)
   mu(u)mu(v)1_(chi_u!=chi_v).                (PEIA.7)
```

One edge is cut with probability `1/2`; any two distinct edges are
both cut with probability `1/4`, whether they share a vertex or not.
Therefore

```text
4 E_chi[J_t(chi)^2]=S_t^2+Q_t.               (PEIA.8)
```

After multiplying by `nu_t` and summing,

```text
E_(alpha,K)
 =4 E_chi[sum_t nu_t J_t(chi)^2]
  -D_(alpha,K).                               (PEIA.9)
```

If

```text
f_chi=mu 1_(chi=0),  g_chi=mu 1_(chi=1),
```

then the cut sum is the separated-color bilinear form

```text
J_t(chi)
 =B_t(f_chi,g_chi)+B_t(g_chi,f_chi)
 =2<f_chi,A_t g_chi>.                         (PEIA.10)
```

Because the diagonal is closed, the planar energy target is
equivalent to the average cut-bilinear square-function bound at
`K^(1+epsilon)`. This is a cleaner bilinear interface, but it is not
a gain: a bound uniform over all colorings is stronger, and the
random masks do not preserve a divisor identity automatically.

A second exact audit uses independent signs with mean `m`. Put
`z=m^2` and let `P_K(z)` be the expected curvature-weighted square
after every edge label is multiplied by `eta_u eta_v`. Then

```text
P_K(z)=D+z C_3+z^2 C_4.                      (PEIA.11)
```

Mean-zero randomization gives only `P_K(0)=D`. Positivity for
`0<=z<1` does not control the deterministic boundary value
`P_K(1)=E`; randomization separates the chaos levels but does not
prove the edge gate.

## Partial-Prefix Fourier Identity

Use `e(theta)=exp(2*pi*i*theta)` and define

```text
F_x(theta)=sum_(K<p<=x)mu(p)e(p theta),
F_*(theta)=sum_(K<q<4K)mu(q)e(q theta),
K_H(theta)=sum_(h=1)^H e(h theta).
```

Orthogonality gives the exact varying-prefix formula

```text
S_K(x,H)
 =integral_0^1 F_x(theta)conj(F_*(theta))
   K_H(theta)dtheta.                          (PEIA.12)
```

Unlike the full-prefix identity, this retains one partial Mobius
polynomial and one fixed full polynomial. It identifies a bilinear
Carleson-type interface without asserting that an existing theorem
has the needed arithmetic hypotheses.

## Signed Current-Interior Anti-Diagonal

Now return to the signed direct pairing before replacing mixed
curvature by its absolute value. In the current-current interior,

```text
i>=1, H>=1, i+H<=K-2, s=2i+H,

kappa_(K,R)(s)
 :=(2/N) integral_0^(pi/N) integral_0^(2pi/N)
   D_R(s*pi/N+u+v)dvdu.                      (PEIA.13)
```

Corollary 11.22Z.5 gives
`Delta Wtilde_K(K+i,H)=kappa_(K,R)(2i+H)`. Thus the signed curvature
depends on the single anti-diagonal coordinate `s`.

For `q_r=2r-1`, put

```text
sigma_(K,r)
 :=sinc(q_r*pi/(2N))sinc(q_r*pi/N).
```

Termwise integration of the finite Dirichlet kernel gives

```text
kappa_(K,R)(s)
 =(4*pi^2/N^3) sum_(r=1)^R sigma_(K,r)
   cos(q_r(s+3/2)pi/N).                      (PEIA.14)
```

Odd-frequency orthogonality over a complete `2N`-point period yields

```text
sum_(s=0)^(2N-1)kappa_(K,R)(s)^2
 =(16*pi^4/N^5)sum_(r=1)^R sigma_(K,r)^2
 <=16*pi^4*R/N^5.                            (PEIA.15)
```

Define the anti-diagonal prefix aggregate

```text
A_K(s)
 :=sum_(i,H>=1,i+H<=K-2,2i+H=s)S_K(K+i,H).
                                                        (PEIA.16)
```

For an edge `(K+a,K+a+g)`, its multiplicity in that line is

```text
N_s(a,g)
 :=[floor((s-g)/2)-max(a,s-K+2)+1]_+,

A_K(s)
 =sum_(a,g>=1,a+g<=K-2)
   mu(K+a)mu(K+a+g)N_s(a,g).                 (PEIA.17)
```

Consequently the complete signed current-current interior part of the
double-Abel pairing is

```text
O_(CC,int,K)
 =sum_(s=3)^(2K-5)kappa_(K,R)(s)A_K(s).      (PEIA.18)
```

Put

```text
Y_(K,r)
 :=sum_(s=3)^(2K-5)A_K(s)
   cos(q_r(s+3/2)pi/N).                      (PEIA.19)
```

Equations (PEIA.14) and (PEIA.18) give the exact low odd-frequency
projection

```text
O_(CC,int,K)
 =(4*pi^2/N^3)sum_(r=1)^R sigma_(K,r)Y_(K,r).
                                                        (PEIA.20)
```

No absolute value over cells or frequencies has entered. Cauchy only
at this final `R`-mode interface gives

```text
|O_(CC,int,K)|^2
 <=(16*pi^4*R/N^6)
   sum_(r=1)^R|Y_(K,r)|^2.                  (PEIA.21)
```

Therefore the concrete sufficient current-interior target

```text
sum_(r=1)^R|Y_(K,r)|^2
 =O_epsilon(N^6*R^(-1)*K^epsilon)           (PEIA.22)
```

for every `epsilon>0` would make this component subpower. Since
`R=ceil(K^(1-alpha/2))`, its scale is
`K^(5+alpha/2+epsilon)`.

The complete odd cosine family gives only

```text
sum_(r=1)^R|Y_(K,r)|^2
 <=N sum_s|A_K(s)|^2.                       (PEIA.23)
```

This Bessel interface does not supply the required Mobius line-energy
gain by geometry alone. The direct projection target is nevertheless
strictly weaker than the absolute planar energy gate: it retains the
Dirichlet signs and asks only for the actual `R` odd-frequency
projection.

## Recent Short-Interval Literature Guard

Menon, *Improved bounds for multiplicative functions in almost all
short intervals* (arXiv:2607.15574, 2026), improves Fourier-uniformity
and averaged Chowla estimates for Mobius/Liouville in almost all short
intervals. The gains relevant here remain logarithmic and are obtained
after averaging over interval bases or shifts. The paper also identifies
roughly inverse-logarithmic decay as the natural limit of the
Matomaki-Radziwill method for its averaged Chowla application.

Those results are valuable context, but they do not imply (PEIA.22):
`Y_(K,r)` is an anchored quadratic edge aggregate with nested prefix
multiplicities, and the present theorem target asks for a power-scale
square-function estimate. They also do not treat the transition,
current-future, or future-future pieces. No almost-all short-interval
statement may therefore be promoted to the missing projection theorem
without an explicit transfer argument.

## Open Gates And Proof Boundary

There are now two deliberately separate theorem-search interfaces:

1. The strong sufficient energy route is exactly equivalent to an
   average random-cut bilinear square function, with the closed
   diagonal correction.
2. The weaker signed route reduces the current-current interior to
   (PEIA.20)-(PEIA.22), before absolute mixed curvature.

The second item does not yet include the transition,
current-future, or future-future cells. Those remain the explicit
boundary and anchor obligations from Corollary 11.22Z.5 and the
earlier truncated-kernel handoffs.

The endpoint incidence formulas, perfect-cancellation witness,
random-cut identity, biased-sign polynomial, partial Fourier formula,
anti-diagonal multiplicity, smoothed odd spectrum, Parseval law, and
projection inequality are exact and independently checked. No
projection estimate, remaining-block estimate, full curvature-energy
bound, full Burnol bound, RH, PF-infinity, `Lambda<=0`, or Clay-prize
result is proved.
"""


def build_payload() -> dict:
    rows = build_rows()
    return {
        "kind": (
            "jensen_window_pf_mertens_planar_"
            "incidence_antidiagonal_handoff"
        ),
        "date": "2026-07-24",
        "status": (
            "exact endpoint incidence, random-cut bilinear, and signed "
            "current-interior anti-diagonal reductions with one open "
            "projection gate"
        ),
        "source": {
            "parent": (
                "outputs/"
                "jensen_window_pf_mertens_planar_edge_gram_vaughan_"
                "handoff.md"
            ),
            "planar_parent": (
                "outputs/jensen_window_pf_mertens_planar_abel_handoff.md"
            ),
            "formal_core": "outputs/formal_core.md",
        },
        "literature": {
            "menon_2026": "https://arxiv.org/abs/2607.15574",
        },
        "parameters": {
            "vertex_interval": "K<n<4K",
            "current_interior": "i>=1, H>=1, i+H<=K-2",
            "antidiagonal": "s=2i+H",
            "mode_cutoff": "R=ceil(K^(1-alpha/2))",
            "normalization": "N=2K+1",
        },
        "rows": rows,
        "summary": {
            "row_count": 43,
            "exact_reduction_count": 32,
            "proof_guard_count": 7,
            "literature_guard_count": 1,
            "countermodel_count": 1,
            "finite_validation_count": 1,
            "open_projection_gate_count": 1,
            "endpoint_expansion_proved": True,
            "perfect_stratum_cancellation_witness_proved": True,
            "cut_bilinear_equivalence_proved": True,
            "partial_prefix_fourier_identity_proved": True,
            "antidiagonal_projection_proved": True,
            "projection_estimate_proved": False,
            "remaining_blocks_proved": False,
            "curvature_energy_proved": False,
            "rh_proved": False,
            "pf_infinity_proved": False,
            "lambda_le_zero_proved": False,
        },
        "proof_boundary": (
            "Exact incidence and current-interior reductions only. "
            "No missing arithmetic projection or boundary/future "
            "estimate is proved."
        ),
    }


def render_note(payload: dict) -> str:
    lines = [NOTE_BODY.strip(), "", "## Claim Ledger", ""]
    lines.append("| ID | Role | Status | Statement | Boundary |")
    lines.append("|---|---|---|---|---|")
    for item in payload["rows"]:
        cells = [
            f"`{item['id']}`",
            f"`{item['role']}`",
            f"`{item['status']}`",
            item["statement"],
            item["proof_boundary"],
        ]
        cells = [cell.replace("|", "\\|") for cell in cells]
        lines.append("| " + " | ".join(cells) + " |")
    summary = payload["summary"]
    lines.extend(
        [
            "",
            "Summary:",
            "",
            f"- rows: {summary['row_count']}",
            f"- exact reductions: {summary['exact_reduction_count']}",
            f"- proof guards: {summary['proof_guard_count']}",
            f"- literature guards: {summary['literature_guard_count']}",
            f"- countermodels: {summary['countermodel_count']}",
            f"- open projection gates: "
            f"{summary['open_projection_gate_count']}",
            "",
        ]
    )
    return "\n".join(lines)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_RESULT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    payload = build_payload()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2) + "\n",
        encoding="utf-8",
    )
    args.note.write_text(render_note(payload), encoding="utf-8")
    summary = payload["summary"]
    print(
        "built planar incidence/anti-diagonal handoff: "
        f"{summary['row_count']} rows, "
        f"{summary['exact_reduction_count']} exact reductions, "
        f"{summary['proof_guard_count']} proof guards, "
        f"{summary['countermodel_count']} countermodel, "
        f"{summary['open_projection_gate_count']} open projection gate"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
