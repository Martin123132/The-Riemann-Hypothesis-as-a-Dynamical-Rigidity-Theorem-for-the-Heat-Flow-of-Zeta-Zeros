#!/usr/bin/env python3
"""Build the Anthropic Weil/rank-trace density-bridge audit."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
import os
from pathlib import Path
import re

import mpmath as mp
from pypdf import PdfReader
import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = "jensen_window_pf_newman_anthropic_weil_density_bridge_gate"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"

PAPER = REPO_ROOT / "work/rh_compute/external/anthropic_zeta23_20260810/claude_more_than_two_thirds_zeta_zeros.pdf"
INFORMAL_NOTE = REPO_ROOT / "work/rh_compute/external/anthropic_zeta23_20260810/anthropic_informal_note_67_percent.pdf"
LOCAL_DEPENDENCIES = {
    "fixed_root_tangency": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_fixed_root_shift_tangency_rigidity_gate.json",
    "polar_collision_cascade": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_polar_heat_collision_cascade_lemma.json",
    "cofinal_boundary_target": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_cofinal_boundary_degree_transfer_target.json",
    "Airy_core": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_endpoint_logit_airy_core_reduction_gate.json",
}
EXPECTED_HASHES = {
    PAPER: "6792988e6cd0e17690621ce898abd5d534f98407741bc7cb14bbe7d07c77d72f",
    INFORMAL_NOTE: "45e0330ad37965e5531fa1f4f11e5bebcae147a5237a3e5b3d029efa7ddf759d",
}


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    certificate: str
    proof_boundary: str


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(text, encoding="utf-8")
    os.replace(temporary, path)


def compact_pdf_text(path: Path) -> str:
    reader = PdfReader(path)
    text = "\n".join(page.extract_text() or "" for page in reader.pages)
    return re.sub(r"[^a-z0-9]+", "", text.lower())


def audit_sources() -> dict:
    for path, expected in EXPECTED_HASHES.items():
        require(path.is_file(), f"missing external source: {path}")
        require(file_hash(path) == expected, f"external source hash drift: {path}")

    paper_text = compact_pdf_text(PAPER)
    note_text = compact_pdf_text(INFORMAL_NOTE)
    for token in (
        "morethantwothirdsofthezerosoftheriemannzetafunctionlieonthecriticalline",
        "ranktraceinequality",
        "rhitselfisoutofreachofthemechanism",
        "06725",
    ):
        require(token in paper_text, f"paper token missing: {token}")
    for token in ("67ofthezeroesareontheline", "guinandweil"):
        require(token in note_text, f"informal-note token missing: {token}")

    return {
        "release_date": "2026-08-10",
        "announcement": "https://www.anthropic.com/research/riemann-zeta",
        "paper": {
            "path": str(PAPER.relative_to(REPO_ROOT)).replace("\\", "/"),
            "sha256": EXPECTED_HASHES[PAPER],
            "url": "https://www-cdn.anthropic.com/564f962e60643842f5fcb4a17c9dbc8f608f1c37.pdf",
        },
        "informal_note": {
            "path": str(INFORMAL_NOTE.relative_to(REPO_ROOT)).replace("\\", "/"),
            "sha256": EXPECTED_HASHES[INFORMAL_NOTE],
            "url": "https://www-cdn.anthropic.com/23455459f8832d06bb175cc0f88d019aed962ef8.pdf",
        },
        "Lean_companion": "https://github.com/anthropics/zeta-23-lean",
        "acceptance_guard": "Newly released external result. This gate does not substitute for broad independent scrutiny or reprove every analytic estimate.",
    }


def constant_certificate() -> dict:
    theta = 1 / sp.sqrt(2)
    c1 = sp.sqrt(2) * sp.tan(theta) / (1 + theta * sp.tan(theta))
    critical = sp.simplify(2 - 1 / c1)
    alternate = sp.Rational(3, 2) - sp.cot(theta) / sp.sqrt(2)
    require(sp.simplify(sp.trigsimp(critical - alternate)) == 0, "constant formulas disagree")

    mp.mp.dps = 100
    theta_mp = 1 / mp.sqrt(2)
    c1_mp = mp.sqrt(2) * mp.tan(theta_mp) / (1 + theta_mp * mp.tan(theta_mp))
    critical_mp = 2 - 1 / c1_mp
    exceptional_mp = 1 - critical_mp
    pair_mp = exceptional_mp / 2
    require(mp.mpf("0.6725007") < critical_mp < mp.mpf("0.6725008"), "critical proportion interval failed")
    require(mp.mpf("0.3274992") < exceptional_mp < mp.mpf("0.3274993"), "exceptional cap interval failed")
    require(mp.mpf("0.1637496") < pair_mp < mp.mpf("0.1637497"), "off-line pair cap interval failed")

    return {
        "theta": mp.nstr(theta_mp, 80),
        "Montgomery_Taylor_c1": mp.nstr(c1_mp, 80),
        "simple_on_line_lower_constant": mp.nstr(critical_mp, 80),
        "exceptional_fraction_upper_cap": mp.nstr(exceptional_mp, 80),
        "off_line_pair_fraction_threshold": mp.nstr(pair_mp, 80),
        "exact_formula": "c_*=3/2-cot(1/sqrt(2))/sqrt(2)=2-1/c_1*",
        "exceptional_definition": "E(T)=N(T,2T)-N_0^s(T,2T)",
        "asymptotic_consequence": "limsup E(T)/N(T,2T)<=delta_*=1-c_*",
    }


def rank_trace_certificate() -> dict:
    c, r, b = sp.symbols("c r b", positive=True)
    tr_p = c * r / 2
    tr_q = c * b
    frob = c**2 * r / 4 + c**2 * b
    rhs = c * tr_p - c**2 * r / 4 + 2 * c * tr_q - c**2 * b
    require(sp.expand(frob - rhs) == 0, "rank-trace equality model failed")

    trace_a, count_n, frob_a = sp.symbols("trace_A N frob_A", real=True)
    assembled = 4 * trace_a - 2 * count_n - frob_a
    require(
        sp.expand(assembled.subs({trace_a: count_n, frob_a: sp.Rational(4, 3) * count_n}) - sp.Rational(2, 3) * count_n) == 0,
        "two-thirds assembly failed",
    )
    return {
        "lemma": "For P>=0, rank(P)<=r and n_+(Q)<=b, ||P+Q||_F^2 >= c tr(P)-c^2 r/4+2c tr(Q)-c^2 b.",
        "c_equals_2": "r>=2tr(P)+4tr(Q)-4b-||P+Q||_F^2.",
        "equality_model": "P=(c/2)Pi_1 and Q=c Pi_2 for orthogonal projections of ranks r and b.",
        "two_thirds_assembly": "tr(A)~N, ||A||_F^2~4N/3, and tr(P)+2b<=N give rank(P)>=(2/3-o(1))N.",
        "scope": "The finite-dimensional algebra is checked here. The paper's explicit-formula and prime-side asymptotics are source-audited, not independently reproved in full.",
    }


def density_bridge_certificate(constants: dict) -> dict:
    t = sp.symbols("t", positive=True)
    dominant_zero_count = t * sp.log(t / (2 * sp.pi)) / (2 * sp.pi)
    require(sp.limit(t / dominant_zero_count, t, sp.oo) == 0, "linear family should have zero relative density")
    require(sp.limit(sp.log(t) / dominant_zero_count, t, sp.oo) == 0, "logarithmic family should have zero relative density")
    require(sp.limit(1 / dominant_zero_count, t, sp.oo) == 0, "finite family should have zero relative density")

    return {
        "required_contradiction": (
            "If one off-line zero exists, prove that for some epsilon>0 and an unbounded sequence T_j, "
            "N_off(T_j,2T_j)/N(T_j,2T_j)>=delta_*+epsilon."
        ),
        "equivalent_pair_target": (
            "If C_off counts symmetric off-line pairs with multiplicity, force "
            "C_off(T_j,2T_j)/N(T_j,2T_j)>delta_*/2."
        ),
        "delta_star": constants["exceptional_fraction_upper_cap"],
        "pair_threshold": constants["off_line_pair_fraction_threshold"],
        "scale_guard": "N(T,2T)~T log(T)/(2 pi). O(1), O(log T), O(T^alpha) for alpha<1, and even O(T) descendants all have zero relative density.",
        "fixed_root_assessment": "The fixed-root shift chain propagates a common Jensen root and a finite Hankel defect; it does not manufacture distinct zeta zeros in height.",
        "polar_cascade_assessment": "The polar cascade increases Jensen degree and multiplicity at one scaled root; degree growth is not zero-count growth.",
        "cofinal_contact_assessment": "The cofinal boundary programme seeks direct absence of contacts. Its local index is sign-definite but presently has no lower count proportional to T log T.",
        "Airy_assessment": "The current Airy branch controls a grouped source remainder. It supplies neither a dyadic off-line-zero count nor an evolution law for a Weil compression.",
        "bridge_closed": False,
    }


def dependency_audit() -> dict[str, dict[str, str]]:
    audited: dict[str, dict[str, str]] = {}
    for name, path in LOCAL_DEPENDENCIES.items():
        require(path.is_file(), f"missing local dependency: {path}")
        audited[name] = {
            "path": str(path.relative_to(REPO_ROOT)).replace("\\", "/"),
            "sha256": file_hash(path),
        }
    return audited


def build_rows(constants: dict, rank_trace: dict, bridge: dict) -> list[GateRow]:
    rows = [
        GateRow("awd_01_source", "external source", "source_audited", "The released paper and informal note are preserved by SHA-256.", "Two PDF snapshots and source URLs.", "External acceptance remains provisional on the release date."),
        GateRow("awd_02_statement", "theorem statement", "source_audited", "The paper claims an unconditional 0.672500... lower proportion of simple on-line zeros.", constants["exact_formula"], "The full analytic proof is not independently reconstructed here."),
        GateRow("awd_03_constant", "constant", "proved", "The two displayed formulas for the Montgomery-Taylor constant agree.", constants["simple_on_line_lower_constant"], "Constant algebra only."),
        GateRow("awd_04_exceptional", "logical consequence", "proved", "The theorem implies limsup exceptional fraction at most delta_*.", constants["exceptional_fraction_upper_cap"], "Conditional on the external theorem being correct."),
        GateRow("awd_05_rank_trace", "linear algebra", "proved", "The rank-trace equality model and two-thirds assembly are exact.", rank_trace["c_equals_2"], rank_trace["scope"]),
        GateRow("awd_06_density_target", "candidate bridge", "open", "A contradiction needs off-line density above delta_* on an unbounded subsequence.", bridge["required_contradiction"], "Infinitely many or positive but smaller density is insufficient."),
        GateRow("awd_07_pair_target", "candidate bridge", "open", "A pair-counting version needs density above delta_*/2.", bridge["equivalent_pair_target"], "Pair multiplicities and window normalization must be retained."),
        GateRow("awd_08_scale", "non-promotion guard", "proved", "Sub-T-log-T descendant counts have zero relative density.", bridge["scale_guard"], "Degree and shift counts cannot be relabelled as zeta-zero counts."),
        GateRow("awd_09_fixed_root", "local corpus comparison", "guard_validated", "The fixed-root and polar cascades do not close the density target.", bridge["fixed_root_assessment"], bridge["polar_cascade_assessment"]),
        GateRow("awd_10_contact", "local corpus comparison", "guard_validated", "The cofinal contact route remains direct contact exclusion, not density amplification.", bridge["cofinal_contact_assessment"], "Closing that route would already prove Lambda<=0 without the proportion theorem."),
        GateRow("awd_11_Airy", "route separation", "guard_validated", "The Airy amplitude/tail obligation is unchanged.", bridge["Airy_assessment"], "No Weil-form evolution theorem is inferred."),
        GateRow("awd_12_decision", "route decision", "parked", "Record the theorem as a strong external constraint but do not pivot the active proof route.", "Continue grouped Airy amplitude and contour tails; reopen only with a genuine T log T amplification or heat-dependent inertia estimate.", "No RH, Lambda<=0, or prize-level conclusion follows."),
    ]
    require(len(rows) == 12, "row count drift")
    return rows


def render_note(artifact: dict) -> str:
    c = artifact["constant_certificate"]
    r = artifact["rank_trace_certificate"]
    b = artifact["density_bridge_certificate"]
    return f"""# Anthropic Weil Rank-Trace Density-Bridge Audit

Date: 2026-08-10

Status: external theorem statement and sources audited; central constant and finite-dimensional rank-trace assembly independently checked; density bridge open; not a proof of RH.

## Exact External Claim

The newly released paper claims unconditionally that

```text
liminf N_0^s(T,2T)/N(T,2T) >= c_*
c_*={c['exact_formula']}
c_*={c['simple_on_line_lower_constant']}
```

Here `N_0^s` counts simple zeros on the critical line and `N` counts all
nontrivial zeros with multiplicity.  The result is a newly released external
claim.  The source paper, informal note, and public Lean companion make it a
serious input, but this local gate does not replace broad independent scrutiny
or independently reprove every analytic estimate.

## Rank-Trace Core

{r['lemma']}

At `c=2`,

```text
{r['c_equals_2']}
```

The equality model is `{r['equality_model']}`  Combining the paper's
zero-side block structure with its first and second moments gives

```text
{r['two_thirds_assembly']}
```

This gate checks that algebra and the optimized constant.  It treats the
explicit-formula localization and prime-side moment evaluation as audited
external inputs rather than claiming a second complete proof.

## Exceptional Cap

Put

```text
E(T)=N(T,2T)-N_0^s(T,2T).
```

The claimed theorem gives

```text
limsup E(T)/N(T,2T) <= delta_*=1-c_*
delta_*={c['exceptional_fraction_upper_cap']}
delta_*/2={c['off_line_pair_fraction_threshold']}
```

Every off-line zero is exceptional, but exceptional also includes a repeated
on-line zero.  Therefore an off-line family can contradict the theorem only
by exceeding the same upper cap.

## Required Density Amplification

The exact candidate bridge is:

```text
{b['required_contradiction']}
```

If symmetric pairs are counted instead, the target is:

```text
{b['equivalent_pair_target']}
```

This is a very large requirement.  Since

```text
N(T,2T) ~ T log(T)/(2 pi),
```

{b['scale_guard']}  In particular, one collision per height unit is still
zero density.  An infinite descendant chain is nowhere near sufficient.

## Comparison With Our Corpus

- Fixed-root shift tangency: {b['fixed_root_assessment']}
- Polar collision cascade: {b['polar_cascade_assessment']}
- Cofinal contact route: {b['cofinal_contact_assessment']}
- Active Airy route: {b['Airy_assessment']}

Thus none of the current exact cascades closes the new density target.  Calling
Jensen-degree growth, shift growth, or local contact index a density of zeta
zeros would be a category error.

## Heat-Flow Inertia Option

The other possible interface would be a heat-dependent family of Weil
compressions whose negative inertia is quantitatively observable from the
source side.  At `t=0`, the classical explicit formula supplies the prime-side
representation used in the paper.  Away from `t=0`, our present heat-flow
source formulas do not supply the corresponding trace and Hilbert-Schmidt
moment theorem.  A naive congruence transport would preserve inertia and is
incompatible with the collision mechanism, so no such transport is promoted.

A useful future theorem would need either:

1. a differential inequality for negative spectral mass of a heat-dependent
   Weil compression with source-side errors under control; or
2. a uniform observability theorem forcing any off-line pair to contribute a
   non-negligible negative eigenvalue across enough independent windows.

Neither theorem is currently available, and sparse zeros arbitrarily close to
the line are precisely what first and second aggregate moments fail to see.

## Route Decision

Do not pivot the main proof programme.  Preserve this result as a strong new
external constraint and a precisely stated optional bridge.  Continue the
grouped Airy amplitude, contour-tail, interior-join, and `T_upper` assembly.
Reopen the density branch only if a mechanism produces order `T log T`
distinct off-line zeros in dyadic windows, or if a genuine heat-dependent
negative-inertia estimate appears.

## Pi Provenance

The proportion constant `c_*` contains no `pi`; it comes from the optimized
Montgomery-Taylor cosine window with parameter `theta=1/sqrt(2)`.  The `pi` in
the density scale comes from the Riemann-von Mangoldt main term
`N(T,2T)~T log(T)/(2 pi)`.  In the paper's Gabor compression, `2 pi` also comes
from the Fourier sampling spacing `h=2 pi/L`.  These are separate appearances
and are not inserted as free normalization constants.

## Proof Boundary

This audit proves the constant identities, rank-trace equality model,
exceptional-cap logic, and zero-density scale guard.  It does not independently
reprove every analytic estimate in the external paper, establish community
acceptance, prove a defect-density amplification theorem, transport Weil
inertia through the heat flow, complete the grouped Airy remainder, prove
`Lambda<=0`, prove RH, or establish a prize-level conclusion.
"""


def main() -> int:
    source_audit = audit_sources()
    constants = constant_certificate()
    rank_trace = rank_trace_certificate()
    bridge = density_bridge_certificate(constants)
    dependencies = dependency_audit()
    rows = build_rows(constants, rank_trace, bridge)
    artifact = {
        "kind": KIND,
        "date": "2026-08-10",
        "status": "external theorem and central rank-trace mechanism audited; density bridge remains open and active Airy route unchanged",
        "passed": True,
        "decision": {
            "external_statement_source_audited": True,
            "full_external_analytic_proof_independently_reproved": False,
            "community_acceptance_treated_as_settled": False,
            "rank_trace_core_checked": True,
            "exceptional_cap_derived": True,
            "density_amplification_proved": False,
            "heat_dependent_Weil_inertia_transport_proved": False,
            "direct_RH_bridge": False,
            "main_Airy_route_changed": False,
            "side_branch": "parked_after_exact_nonbridge_reduction",
        },
        "source_audit": source_audit,
        "constant_certificate": constants,
        "rank_trace_certificate": rank_trace,
        "density_bridge_certificate": bridge,
        "dependencies": dependencies,
        "rows": [asdict(row) for row in rows],
        "proof_boundary": "The source snapshots, theorem statement, optimized constant, central rank-trace algebra, exceptional-cap consequence, and density-scale non-promotion guard are checked. No full independent analytic reproof, community acceptance claim, density amplification, heat-dependent Weil-inertia theorem, grouped Airy remainder, Lambda<=0, RH, or prize-level conclusion is proved.",
    }
    atomic_write(RESULT_PATH, json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    atomic_write(NOTE_PATH, render_note(artifact))
    print(
        "built Anthropic Weil density-bridge audit: "
        f"c*={constants['simple_on_line_lower_constant'][:18]}, "
        f"delta*={constants['exceptional_fraction_upper_cap'][:18]}, "
        f"{len(rows)} rows, bridge_closed={bridge['bridge_closed']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
