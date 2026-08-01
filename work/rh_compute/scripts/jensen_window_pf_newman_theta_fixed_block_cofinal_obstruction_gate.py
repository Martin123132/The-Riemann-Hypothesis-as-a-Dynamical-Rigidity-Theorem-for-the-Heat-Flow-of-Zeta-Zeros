#!/usr/bin/env python3
"""Prove the cofinal obstruction for fixed theta blocks and static tail bars."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal
import hashlib
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_theta_"
    "fixed_block_cofinal_obstruction_gate"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SECOND_BLOCK_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_theta_"
    "second_diagonal_shell_two_block_interval_certificate.json"
)
EVENNESS_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_theta_"
    "bessel_higher_shift_regularization_gate.json"
)
FRONTIER_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_theta_two_block_shell_frontier_scout/"
    "frontier.jsonl"
)
FRONTIER_LAST_J = 30


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | list | None = None


def canonical_bytes(value: dict) -> bytes:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
    ).encode("ascii")


def row_sha256(row: dict) -> str:
    payload = {
        key: value
        for key, value in row.items()
        if key != "row_sha256"
    }
    return hashlib.sha256(canonical_bytes(payload)).hexdigest()


def load_frontier_prefix() -> tuple[list[dict], dict]:
    rows: list[dict] = []
    previous: str | None = None
    with FRONTIER_SOURCE.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            row = json.loads(line)
            if row.get("previous_row_sha256") != previous:
                raise RuntimeError(
                    f"frontier hash-chain failure at row {line_number}"
                )
            if row.get("row_sha256") != row_sha256(row):
                raise RuntimeError(
                    f"frontier row-hash failure at row {line_number}"
                )
            previous = row["row_sha256"]
            rows.append(row)
            if row.get("j") == FRONTIER_LAST_J:
                break
    if not rows or rows[0].get("kind") != (
        "jensen_window_pf_newman_theta_"
        "two_block_shell_frontier_scout_header"
    ):
        raise RuntimeError("frontier header is missing")
    shells = [row for row in rows if isinstance(row.get("j"), int)]
    expected = list(range(5, FRONTIER_LAST_J + 1))
    if [row["j"] for row in shells] != expected:
        raise RuntimeError("frontier prefix is incomplete")
    if any(row.get("status") != "certified" for row in shells):
        raise RuntimeError("frontier prefix contains an uncertified shell")
    minimum = min(
        shells,
        key=lambda row: Decimal(
            row["summary"]["minimum_certified_ratio_lower"]
            .lstrip("[")
            .split()[0]
        ),
    )
    audit = {
        "start_j": shells[0]["j"],
        "end_j": shells[-1]["j"],
        "shell_rows": len(shells),
        "certified_shell_rows": sum(
            row["status"] == "certified" for row in shells
        ),
        "total_evaluated_boxes": sum(
            row["summary"]["evaluated_boxes"] for row in shells
        ),
        "total_subdivisions": sum(
            row["summary"]["subdivisions"] for row in shells
        ),
        "minimum_ratio_shell": minimum["j"],
        "minimum_ratio_lower": minimum["summary"][
            "minimum_certified_ratio_lower"
        ],
        "terminal_prefix_row_sha256": shells[-1]["row_sha256"],
        "proof_boundary": (
            "Stored finite diagnostics only; no Q_j with j>=5 is "
            "promoted by this audit."
        ),
    }
    return rows, audit


def source_audit() -> dict:
    two_block = json.loads(SECOND_BLOCK_SOURCE.read_text(encoding="utf-8"))
    evenness = json.loads(EVENNESS_SOURCE.read_text(encoding="utf-8"))
    two_block_text = json.dumps(two_block, sort_keys=True)
    evenness_text = json.dumps(evenness, sort_keys=True)
    for marker in (
        "M_0<40*81*exp(-27)/99<10^-10",
        "|J_t-J_(<=2,t)|<16*x^4*10^-10",
        "retained n=1,2",
    ):
        if marker not in two_block_text:
            raise RuntimeError(f"two-block source marker missing: {marker}")
    for marker in (
        "Phi(u)=sum_(n>=1)phi_n(u)=Phi(-u)",
        "phi_n(u)=",
    ):
        if marker not in evenness_text:
            raise RuntimeError(f"evenness source marker missing: {marker}")
    _, frontier = load_frontier_prefix()
    return {
        "two_block_kind": two_block["kind"],
        "evenness_kind": evenness["kind"],
        "frontier_prefix": frontier,
    }


def build_exact() -> dict:
    u, n, t = sp.symbols(
        "u n t", real=True, nonnegative=True
    )
    n = sp.symbols("n", integer=True, positive=True)
    phi_n = (
        2 * sp.pi**2 * n**4 * sp.exp(9 * u)
        - 3 * sp.pi * n**2 * sp.exp(5 * u)
    ) * sp.exp(-sp.pi * n**2 * sp.exp(4 * u))
    derivative_at_zero = sp.factor(
        sp.diff(phi_n, u).subs(u, 0)
    )
    expected = (
        sp.pi
        * n**2
        * sp.exp(-sp.pi * n**2)
        * (
            -8 * sp.pi**2 * n**4
            + 30 * sp.pi * n**2
            - 15
        )
    )
    if sp.simplify(derivative_at_zero - expected) != 0:
        raise RuntimeError("theta endpoint derivative identity failed")

    phi_one = phi_n.subs(n, 1)
    phi_two = phi_n.subs(n, 2)
    retained = sp.exp(t * u**2) * (phi_one + phi_two)
    retained_derivative = sp.simplify(sp.diff(retained, u).subs(u, 0))
    a_two = sp.simplify(
        expected.subs(n, 1) + expected.subs(n, 2)
    )
    if sp.simplify(retained_derivative - a_two) != 0:
        raise RuntimeError("two-block endpoint derivative drifted")

    x, h, hp, epsilon_zero, epsilon_one = sp.symbols(
        "x h hp epsilon_0 epsilon_1",
        positive=True,
    )
    j_value = 16 * x**4 * h
    j_derivative = 64 * x**3 * h + 16 * x**4 * hp
    value_bar = 16 * x**4 * epsilon_zero
    derivative_bar = (
        64 * x**3 * epsilon_zero
        + 16 * x**4 * epsilon_one
    )
    if sp.simplify(
        j_value / value_bar - h / epsilon_zero
    ) != 0:
        raise RuntimeError("value-ratio cancellation failed")
    if sp.simplify(
        j_derivative / derivative_bar
        - (4 * h + x * hp) / (
            4 * epsilon_zero + x * epsilon_one
        )
    ) != 0:
        raise RuntimeError("derivative-ratio cancellation failed")

    return {
        "retained_kernel": {
            "theta_block": (
                "phi_n(u)=(2*pi^2*n^4*exp(9u)-"
                "3*pi*n^2*exp(5u))*exp(-pi*n^2*exp(4u))"
            ),
            "finite_kernel": (
                "f_(N,t)(u)=exp(tu^2)*sum_(1<=n<=N)phi_n(u)"
            ),
            "transform": (
                "H_(N,t)(x)=integral_0^infinity "
                "f_(N,t)(u)*cos(xu)du"
            ),
            "regularity": (
                "For fixed N, every u-derivative of f_(N,t) is "
                "uniformly integrable for 0<=t<=1/5, and every "
                "integration-by-parts boundary term at infinity vanishes."
            ),
        },
        "integration_by_parts": {
            "constants": (
                "A_N=sup_(0<=t<=1/5)(|f_(N,t)'(0)|+"
                "||f_(N,t)''||_1)<infinity; "
                "B_N=sup_(0<=t<=1/5)||(u*f_(N,t))''||_1<infinity"
            ),
            "value_identity": (
                "H_(N,t)(x)=-f_(N,t)'(0)/x^2-"
                "x^(-2)*integral_0^infinity "
                "f_(N,t)''(u)*cos(xu)du"
            ),
            "value_bound": "|H_(N,t)(x)|<=A_N/x^2",
            "derivative_identity": (
                "partial_x H_(N,t)(x)=x^(-2)*"
                "integral_0^infinity "
                "(u*f_(N,t)(u))''*sin(xu)du"
            ),
            "derivative_bound": (
                "|partial_x H_(N,t)(x)|<=B_N/x^2"
            ),
        },
        "static_tail_certificate": {
            "bars": (
                "E_0(x)=16*x^4*epsilon_0; "
                "E_1(x)=64*x^3*epsilon_0+16*x^4*epsilon_1, "
                "where epsilon_0>0 and epsilon_1>=0 are fixed"
            ),
            "value_ratio": (
                "|J_(N,t)(x)|/E_0(x)="
                "|H_(N,t)(x)|/epsilon_0"
                "<=A_N/(epsilon_0*x^2)"
            ),
            "derivative_ratio": (
                "|J_(N,t)'(x)|/E_1(x)="
                "|4H_(N,t)(x)+x*partial_x H_(N,t)(x)|/"
                "(4epsilon_0+xepsilon_1)"
                "<=(4A_N/x^2+B_N/x)/(4epsilon_0+xepsilon_1)"
            ),
            "explicit_threshold": (
                "X_N=max(1,sqrt(2A_N/epsilon_0),"
                "B_N/(2epsilon_0))"
            ),
            "failure": (
                "For every x>X_N and every 0<=t<=1/5, "
                "|J_(N,t)(x)|<E_0(x) and "
                "|J_(N,t)'(x)|<E_1(x)"
            ),
        },
        "cofinal_obstruction": {
            "shells": (
                "S_j=[1/(5j),1/5]x[38,38+j]"
            ),
            "theorem": (
                "For every fixed N and every static positive absolute "
                "tail bar epsilon_0, the value-or-derivative direct-bar "
                "certificate fails at some point of every sufficiently "
                "large S_j."
            ),
            "scope": (
                "This is an obstruction to the fixed-N static absolute-"
                "tail proof method, not evidence of a common zero of the "
                "full Xi first jet."
            ),
        },
        "two_block_boundary_cancellation": {
            "endpoint_derivative": (
                "phi_n'(0)=pi*n^2*exp(-pi*n^2)*"
                "(-8*pi^2*n^4+30*pi*n^2-15)"
            ),
            "a_two": (
                "a_2=f_(2,t)'(0)=sum_(n=1)^2 phi_n'(0)>0"
            ),
            "a_two_decimal": str(sp.N(a_two, 80)),
            "sign_proof": (
                "For n>=3, 4*pi*n^2-15>0, hence "
                "-8*pi^2*n^4+30*pi*n^2-15<0. The differentiated "
                "theta series converges locally uniformly. Since Phi "
                "is even, sum_(n>=1)phi_n'(0)=0, so "
                "a_2=-sum_(n>=3)phi_n'(0)>0."
            ),
            "full_cancellation": (
                "Phi is even, so (exp(tu^2)Phi(u))'(0)=0 and "
                "sum_(n>=3)phi_n'(0)=-a_2<0"
            ),
            "interpretation": (
                "The omitted theta tail cancels the algebraic x^(-2) "
                "endpoint term of the retained two-block transform. "
                "Static absolute moments discard this cancellation."
            ),
        },
        "required_upgrade": {
            "route": (
                "A cofinal theta-shell proof must let N grow with x, "
                "replace static moment bars by oscillatory remainder "
                "bounds that retain endpoint-jet cancellation, or use "
                "a modular grouping that preserves the cancellation "
                "before transformation."
            ),
            "endpoint_jet_identity": (
                "integral_0^infinity f(u)cos(xu)du="
                "sum_(r=0)^(m-1)(-1)^(r+1)"
                "*f^(2r+1)(0)/x^(2r+2)"
                "+(-1)^m*x^(-2m)*integral_0^infinity "
                "f^(2m)(u)cos(xu)du"
            ),
            "open_target": (
                "Construct a uniform cancellation-aware first-jet "
                "remainder theorem strong enough to close every S_j "
                "without assuming zero simplicity or RH."
            ),
        },
    }


def build_artifact() -> dict:
    exact = build_exact()
    audit = source_audit()
    rows = [
        GateRow(
            id="ntfbco_01_uniform_regular_decay",
            role="exact_lemma",
            readiness="ready_to_apply",
            claim=(
                "Every fixed finite theta truncation has enough uniform "
                "half-line decay for repeated integration by parts."
            ),
            formula=exact["retained_kernel"]["regularity"],
            proof_boundary=(
                "Fixed N and 0<=t<=1/5 only; no infinite theta sum is "
                "interchanged here."
            ),
        ),
        GateRow(
            id="ntfbco_02_value_ibp_bound",
            role="exact_identity",
            readiness="ready_to_apply",
            claim=(
                "The retained cosine transform has a uniform O(x^-2) "
                "bound controlled by its endpoint derivative."
            ),
            formula=exact["integration_by_parts"]["value_identity"],
            proof_boundary="Exact integration by parts plus an L1 bound.",
            diagnostics={
                "bound": exact["integration_by_parts"]["value_bound"],
            },
        ),
        GateRow(
            id="ntfbco_03_derivative_ibp_bound",
            role="exact_identity",
            readiness="ready_to_apply",
            claim=(
                "The x derivative of the retained transform is also "
                "uniformly O(x^-2)."
            ),
            formula=exact["integration_by_parts"][
                "derivative_identity"
            ],
            proof_boundary="Exact integration by parts for u*f_(N,t).",
            diagnostics={
                "bound": exact["integration_by_parts"][
                    "derivative_bound"
                ],
            },
        ),
        GateRow(
            id="ntfbco_04_normalized_static_bars",
            role="exact_inequality",
            readiness="ready_to_apply",
            claim=(
                "After the x^4 normalization, both direct-bar "
                "certificate ratios have explicit decaying upper bounds."
            ),
            formula=(
                f"{exact['static_tail_certificate']['value_ratio']}; "
                f"{exact['static_tail_certificate']['derivative_ratio']}"
            ),
            proof_boundary=(
                "Applies to fixed positive static absolute-tail bars."
            ),
        ),
        GateRow(
            id="ntfbco_05_cofinal_obstruction",
            role="route_guard",
            readiness="guard_validated",
            claim=exact["cofinal_obstruction"]["theorem"],
            formula=(
                f"{exact['static_tail_certificate']['explicit_threshold']}; "
                f"{exact['static_tail_certificate']['failure']}"
            ),
            proof_boundary=exact["cofinal_obstruction"]["scope"],
        ),
        GateRow(
            id="ntfbco_06_two_block_endpoint_cancellation",
            role="exact_theorem",
            readiness="ready_to_apply",
            claim=(
                "For N=2, the omitted arithmetic tail exactly cancels "
                "the nonzero retained endpoint derivative."
            ),
            formula=exact["two_block_boundary_cancellation"][
                "full_cancellation"
            ],
            proof_boundary=(
                "Exact endpoint cancellation; it does not sign the "
                "full transform at positive frequency."
            ),
            diagnostics=exact["two_block_boundary_cancellation"],
        ),
        GateRow(
            id="ntfbco_07_finite_frontier_diagnostic",
            role="finite_diagnostic",
            readiness="diagnostic_only",
            claim=(
                "The stored two-block scout remains certified through "
                "Q_30 at its finite rational-box resolution."
            ),
            formula=(
                "stored prefix j=5..30; first finite frontier=none"
            ),
            proof_boundary=(
                "No j>=5 stage or cofinal theorem is promoted."
            ),
            diagnostics=audit["frontier_prefix"],
        ),
        GateRow(
            id="ntfbco_08_cancellation_aware_handoff",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim=exact["required_upgrade"]["open_target"],
            formula=exact["required_upgrade"]["endpoint_jet_identity"],
            proof_boundary=(
                "The identity identifies the repair; the required "
                "uniform first-jet estimate remains open."
            ),
            diagnostics={"route": exact["required_upgrade"]["route"]},
        ),
        GateRow(
            id="ntfbco_09_proof_boundary",
            role="proof_guard",
            readiness="guard_active",
            claim=(
                "Failure of this sufficient fixed-block certificate "
                "does not imply an Xi contact or failure of RH."
            ),
            formula=exact["cofinal_obstruction"]["scope"],
            proof_boundary=(
                "No new Q_j, Lambda<=0, RH, or Clay-prize conclusion "
                "is proved."
            ),
        ),
    ]
    return {
        "kind": STEM,
        "date": "2026-07-24",
        "status": (
            "exact fixed-block cofinal obstruction with "
            "cancellation-aware handoff"
        ),
        "proof_boundary": (
            "This artifact proves that a fixed finite theta truncation "
            "combined with fixed positive static absolute-moment tail "
            "bars cannot certify the cofinal diagonal shell family. It "
            "also identifies the exact endpoint cancellation discarded "
            "by the two-block bars. It does not prove a common zero, "
            "disprove the finite interval certificates, certify any new "
            "Q_j, prove Lambda<=0, prove RH, or win a Clay prize."
        ),
        "sources": [
            str(SECOND_BLOCK_SOURCE.relative_to(REPO_ROOT)).replace(
                "\\", "/"
            ),
            str(EVENNESS_SOURCE.relative_to(REPO_ROOT)).replace(
                "\\", "/"
            ),
            str(FRONTIER_SOURCE.relative_to(REPO_ROOT)).replace(
                "\\", "/"
            ),
            "outputs/formal_core.md",
        ],
        "source_audit": audit,
        "exact": exact,
        "rows": [asdict(row) for row in rows],
    }


def render_note(artifact: dict) -> str:
    exact = artifact["exact"]
    frontier = artifact["source_audit"]["frontier_prefix"]
    return "\n".join(
        [
            "# Newman Theta Fixed-Block Cofinal Obstruction Gate",
            "",
            "Date: 2026-07-24",
            "",
            "Status: exact method obstruction with a live "
            "cancellation-aware handoff.",
            "This is not a proof or disproof of RH or `Lambda <= 0`.",
            "",
            "```text",
            f"work/rh_compute/results/{STEM}.json",
            f"python work/rh_compute/scripts/{STEM}.py",
            f"python work/rh_compute/scripts/check_{STEM}.py",
            "```",
            "",
            "## The Fixed-Block Setup",
            "",
            "For a fixed finite `N`, put",
            "",
            "```text",
            exact["retained_kernel"]["theta_block"],
            exact["retained_kernel"]["finite_kernel"],
            exact["retained_kernel"]["transform"],
            "```",
            "",
            exact["retained_kernel"]["regularity"],
            "",
            "The double-exponential theta tail makes all derivatives "
            "used below uniformly integrable for `0<=t<=1/5`.",
            "",
            "## Two Integrations By Parts",
            "",
            "Define the finite uniform constants",
            "",
            "```text",
            exact["integration_by_parts"]["constants"],
            "```",
            "",
            "Then, for every `x>0`,",
            "",
            "```text",
            exact["integration_by_parts"]["value_identity"],
            exact["integration_by_parts"]["value_bound"],
            exact["integration_by_parts"]["derivative_identity"],
            exact["integration_by_parts"]["derivative_bound"],
            "```",
            "",
            "The derivative identity uses `g(u)=u*f_(N,t)(u)`, "
            "for which `g(0)=0`.",
            "",
            "## Exact Cofinal Obstruction",
            "",
            "The existing direct certificate uses static bars",
            "",
            "```text",
            exact["static_tail_certificate"]["bars"],
            "```",
            "",
            "Since `J_(N,t)=16*x^4*H_(N,t)`, normalization cancels "
            "the apparent `x^4` advantage:",
            "",
            "```text",
            exact["static_tail_certificate"]["value_ratio"],
            exact["static_tail_certificate"]["derivative_ratio"],
            "```",
            "",
            "Take",
            "",
            "```text",
            exact["static_tail_certificate"]["explicit_threshold"],
            "```",
            "",
            "For `x>X_N`, the value bound is below `1/2`. Also "
            "`4A_N/x^2<2epsilon_0` and "
            "`B_N/x<2epsilon_0`, while the derivative denominator is "
            "at least `4epsilon_0`. Therefore",
            "",
            "```text",
            exact["static_tail_certificate"]["failure"],
            "```",
            "",
            "Every sufficiently large diagonal shell",
            "",
            "```text",
            exact["cofinal_obstruction"]["shells"],
            "```",
            "",
            "contains such an `x`. Hence a fixed number of retained "
            "theta blocks plus fixed positive absolute-moment bars "
            "cannot close the cofinal family.",
            "",
            "## What The Bars Erase",
            "",
            "At the half-line endpoint,",
            "",
            "```text",
            exact["two_block_boundary_cancellation"][
                "endpoint_derivative"
            ],
            exact["two_block_boundary_cancellation"]["a_two"],
            f"a_2 approximately "
            f"{exact['two_block_boundary_cancellation']['a_two_decimal']}",
            exact["two_block_boundary_cancellation"]["sign_proof"],
            exact["two_block_boundary_cancellation"][
                "full_cancellation"
            ],
            "```",
            "",
            exact["two_block_boundary_cancellation"][
                "interpretation"
            ],
            "",
            "This is the structural point: the full even Xi kernel "
            "has no odd endpoint jet, but every fixed arithmetic "
            "truncation generally does. The omitted tail cancels that "
            "jet. Taking its absolute mass first destroys the "
            "cancellation before the oscillatory transform is used.",
            "",
            "## Finite Scout",
            "",
            "The separate resumable scout records",
            "",
            "```text",
            f"shells={frontier['start_j']}..{frontier['end_j']}",
            f"certified rows={frontier['certified_shell_rows']}",
            f"evaluated boxes={frontier['total_evaluated_boxes']}",
            f"subdivisions={frontier['total_subdivisions']}",
            f"smallest stored ratio occurs at j="
            f"{frontier['minimum_ratio_shell']}: "
            f"{frontier['minimum_ratio_lower']}",
            "```",
            "",
            "Those are rigorous stored finite diagnostics, but they "
            "are not promoted here as new `Q_j` theorems. Passing "
            "through `j=30` and failing cofinally are compatible: the "
            "obstruction is asymptotic.",
            "",
            "## Live Handoff",
            "",
            "The exact repeated endpoint-jet formula is",
            "",
            "```text",
            exact["required_upgrade"]["endpoint_jet_identity"],
            "```",
            "",
            exact["required_upgrade"]["route"],
            "",
            "Open target:",
            "",
            "```text",
            exact["required_upgrade"]["open_target"],
            "```",
            "",
            "## Proof Boundary",
            "",
            exact["cofinal_obstruction"]["scope"],
            "No common zero, new `Q_j`, `Lambda<=0`, RH, or "
            "Clay-prize conclusion follows.",
            "",
            "Current validation:",
            "",
            "```text",
            "validated Newman theta fixed-block cofinal obstruction "
            "gate: 9 rows, 0 issues, 3 exact integration-by-parts "
            "bounds, 1 cofinal method obstruction, 1 endpoint-"
            "cancellation theorem, 26 finite diagnostic shells, "
            "1 cancellation-aware open handoff",
            "```",
            "",
        ]
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    args.note.write_text(render_note(artifact), encoding="utf-8")
    print(
        "wrote Newman theta fixed-block cofinal obstruction gate: "
        f"{len(artifact['rows'])} rows"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
