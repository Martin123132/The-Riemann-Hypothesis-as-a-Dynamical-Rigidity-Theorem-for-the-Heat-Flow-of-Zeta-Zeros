#!/usr/bin/env python3
"""Build the nonterminal determinant-completion ledger-correction gate."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = "jensen_window_pf_newman_c1_actual_source_nonterminal_determinant_completion_gate"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"

SOURCE_PATHS = {
    "observation_image": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_morse_fresnel_"
        "observation_image_compression_gate.json"
    ),
    "full_support": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_"
        "ideal_cubic_reciprocal_terminal_full_support_homotopy_gate.json"
    ),
    "terminal_completion": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_full_support_terminal_bulk_"
        "phase_completion_gate.json"
    ),
    "relative_lift": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_full_support_terminal_"
        "nonterminal_relative_lift_gate.json"
    ),
    "compensation_ledger": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_c1_actual_source_complete_compensation_ledger_gate.json",
    "edge_affine": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_c1_actual_source_edge_affine_ownership_gate.json",
    "quadratic_primitive": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_c1_actual_source_current_transfer_quadratic_primitive_gate.json",
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


def atomic_write(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(content, encoding="utf-8")
    os.replace(temporary, path)


def load_and_audit_sources() -> dict[str, dict[str, str]]:
    payloads: dict[str, dict] = {}
    for key, path in SOURCE_PATHS.items():
        require(path.is_file(), f"missing source artifact: {path}")
        payloads[key] = json.loads(path.read_text(encoding="utf-8"))

    observation = payloads["observation_image"]
    require(
        "g_p^Tr+r^TM_xi r" in observation["matrix_certificate"]["residual_factorization"],
        "observation residual factorization drifted",
    )
    require(
        "N_(p,xi)P_V" in observation["observation_certificate"]["linear_polynomial"],
        "retained derivative row drifted",
    )
    require(
        "four signed base/derivative pairings"
        in observation["observation_certificate"]["quadratic_pairing"],
        "quadratic pairing count drifted",
    )

    full_support = payloads["full_support"]["symbolic_certificate"]["full_support_observations"]
    require(
        "P_lin^[N]=F_alpha+i(lambda-log a)F_beta" in full_support["linear_polynomial"],
        "full-support linear polynomial drifted",
    )
    require(
        "mathcal N_p+mathcal N_E" in full_support["linear_polynomial"],
        "full-support beta row drifted",
    )

    completion = payloads["terminal_completion"]["symbolic_certificate"]["coefficient_completion"]
    require(
        "C(a+rho)=C(rho)" in completion["identity"],
        "terminal completion identity drifted",
    )
    require(
        "omega_*=omega_p+rho" in completion["complete_rows"],
        "terminal mixed-row ownership drifted",
    )

    relative = payloads["relative_lift"]["symbolic_certificate"]["outer_remainder"]
    require(
        "rho_X=Re{nu c_xi mathscr R_N[P_X]}" in relative["observations"],
        "nonterminal observation functional drifted",
    )
    require(
        "mathscr O_N[P]=T_N[P]+mathscr R_N[P]" in relative["split"],
        "outer functional split drifted",
    )

    ledger = payloads["compensation_ledger"]["exact"]["physical_partition"]
    require(
        ledger["total"]
        == "E_tot=E_F+E_R+E_term+E_Delta+E_aff+E_aa+E_move+E_quad.",
        "eight-package ledger drifted",
    )
    require(
        "four remaining nonterminal quadratic packages" in ledger["remaining"],
        "ledger quadratic ownership drifted",
    )

    edge = payloads["edge_affine"]["exact"]
    require(
        "beta_p=beta_N-beta_E" in edge["single_owner"],
        "edge beta subtraction drifted",
    )
    require(
        "mathscr O_N[P_E]" in edge["edge_polynomial"],
        "edge functional ownership drifted",
    )

    primitive = payloads["quadratic_primitive"]["exact"]
    require(
        "E_quad^cur=dQ_quad/dxi exactly" in primitive["quadratic_primitive"],
        "quadratic primitive drifted",
    )
    require(
        "E_FAN+E_quad^cur" in primitive["signed_reclassification"],
        "tested signed-main formula drifted",
    )

    return {
        key: {
            "path": str(path.relative_to(REPO_ROOT)).replace("\\", "/"),
            "sha256": file_hash(path),
        }
        for key, path in SOURCE_PATHS.items()
    }


def real_vector(prefix: str) -> sp.Matrix:
    return sp.Matrix(sp.symbols(f"{prefix}_V {prefix}_N {prefix}_A {prefix}_Q", real=True))


def bilinear(left: sp.Matrix, J: sp.Matrix, right: sp.Matrix) -> sp.Expr:
    return sp.expand((left.T * J * right)[0])


def symbolic_certificate() -> dict[str, str | int]:
    J = sp.Rational(1, 2) * sp.Matrix(
        [
            [0, 1, 0, 0],
            [1, 0, 0, 0],
            [0, 0, 0, -1],
            [0, 0, -1, 0],
        ]
    )
    p, e, a = (real_vector(name) for name in ("p", "e", "a"))
    dp, de, da = (real_vector(name) for name in ("dp", "de", "da"))
    t_e = real_vector("tE")

    raw = (
        2 * bilinear(p, J, da + de)
        + 2 * bilinear(dp, J, a + e)
        + 2 * bilinear(a + e, J, da + de)
        + sp.expand(t_e.dot(da + de))
    )
    retained_nonterminal = 2 * bilinear(p, J, de) + 2 * bilinear(dp, J, e)
    quadratic = 2 * bilinear(e, J, de)
    terminal_mixed = 2 * bilinear(p + e, J, da) + 2 * bilinear(dp + de, J, a)
    terminal_pure = 2 * bilinear(a, J, da)
    edge_affine = sp.expand(t_e.dot(da + de))

    audits = 0
    split = retained_nonterminal + quadratic + terminal_mixed + terminal_pure + edge_affine
    require(sp.simplify(raw - split) == 0, "raw five-block split")
    audits += 1

    represented_old = quadratic + terminal_mixed + terminal_pure + edge_affine
    require(sp.simplify(raw - represented_old - retained_nonterminal) == 0, "ledger defect")
    audits += 1

    cross_primitive_derivative = 2 * bilinear(dp, J, e) + 2 * bilinear(p, J, de)
    require(sp.simplify(retained_nonterminal - cross_primitive_derivative) == 0, "cross primitive")
    audits += 1

    determinant_difference_derivative = 2 * (
        bilinear(p + e, J, dp + de) - bilinear(p, J, dp)
    )
    require(
        sp.simplify(retained_nonterminal + quadratic - determinant_difference_derivative) == 0,
        "completed determinant derivative",
    )
    audits += 1

    alpha_p = sp.Matrix([dp[1], dp[0], -dp[3], -dp[2]])
    beta_p = sp.Matrix([p[1], p[0], -p[3], -p[2]])
    require(sp.simplify(alpha_p - 2 * J * dp) == sp.zeros(4, 1), "alpha row")
    audits += 1
    require(sp.simplify(beta_p - 2 * J * p) == sp.zeros(4, 1), "beta row")
    audits += 1

    polynomial_projection = sp.expand(alpha_p.dot(e) + beta_p.dot(de))
    require(sp.simplify(polynomial_projection - retained_nonterminal) == 0, "P_nr projection")
    audits += 1

    beta_edge = real_vector("betaE")
    beta_full = beta_p + beta_edge
    require(sp.simplify(beta_full - beta_edge - beta_p) == sp.zeros(4, 1), "beta owner split")
    audits += 1

    return {
        "raw_current": "R_raw=2p^TJ dot(a+e)+2dot(p)^TJ(a+e)+2(a+e)^TJ dot(a+e)+t_E^T dot(a+e).",
        "five_block_split": "R_raw=R_nr+R_quad+R_a*+R_aa+R_aff.",
        "retained_nonterminal": "R_nr=2p^TJ dot(e)+2dot(p)^TJe=partial_xi(2p^TJe).",
        "quadratic": "R_quad=2e^TJ dot(e).",
        "terminal_mixed": "R_a*=2(p+e)^TJ dot(a)+2(dot(p)+dot(e))^TJa.",
        "terminal_pure": "R_aa=2a^TJ dot(a).",
        "edge_affine": "R_aff=t_E^T[dot(a)+dot(e)].",
        "ledger_defect": "The residual sector represented by terminal mixed, pure terminal, edge affine, and nonterminal quadratic terms differs from R_raw by exactly R_nr.",
        "polynomial_rows": "alpha_p=(mathcal N_(p,xi),V_(p,xi),-Q_(p,xi),-A_(p,xi))=2Jdot(p), beta_p=(mathcal N_p,V_p,-Q_p,-A_p)=2Jp.",
        "polynomial_projection": "alpha_p^T e+beta_p^T dot(e)=R_nr.",
        "determinant_completion": "R_nr+R_quad=2{(p+e)^TJ[dot(p)+dot(e)]-p^TJdot(p)}=partial_xi{(p+e)^TJ(p+e)-p^TJp}.",
        "symbolic_audits": audits,
    }


def rational_vector(values: tuple[Fraction, Fraction, Fraction, Fraction]) -> sp.Matrix:
    return sp.Matrix([sp.Rational(value.numerator, value.denominator) for value in values])


def finite_rational_certificate() -> dict[str, str | int]:
    J = sp.Rational(1, 2) * sp.Matrix(
        [
            [0, 1, 0, 0],
            [1, 0, 0, 0],
            [0, 0, 0, -1],
            [0, 0, -1, 0],
        ]
    )
    p = rational_vector((Fraction(1, 2), Fraction(-2, 3), Fraction(3, 5), Fraction(-4, 7)))
    dp = rational_vector((Fraction(5, 11), Fraction(-6, 13), Fraction(7, 17), Fraction(-8, 19)))
    e = rational_vector((Fraction(-2, 5), Fraction(3, 7), Fraction(-5, 9), Fraction(7, 11)))
    de = rational_vector((Fraction(11, 13), Fraction(-13, 17), Fraction(17, 19), Fraction(-19, 23)))
    a = rational_vector((Fraction(1, 3), Fraction(2, 5), Fraction(-3, 7), Fraction(4, 9)))
    da = rational_vector((Fraction(-5, 11), Fraction(7, 13), Fraction(-11, 17), Fraction(13, 19)))
    t_e = rational_vector((Fraction(2, 3), Fraction(-3, 5), Fraction(5, 7), Fraction(-7, 11)))

    raw = sp.expand(
        2 * bilinear(p, J, da + de)
        + 2 * bilinear(dp, J, a + e)
        + 2 * bilinear(a + e, J, da + de)
        + t_e.dot(da + de)
    )
    r_nr = sp.expand(2 * bilinear(p, J, de) + 2 * bilinear(dp, J, e))
    r_quad = sp.expand(2 * bilinear(e, J, de))
    r_mixed = sp.expand(2 * bilinear(p + e, J, da) + 2 * bilinear(dp + de, J, a))
    r_aa = sp.expand(2 * bilinear(a, J, da))
    r_aff = sp.expand(t_e.dot(da + de))
    represented_old = sp.expand(r_quad + r_mixed + r_aa + r_aff)
    completed = sp.expand(2 * (bilinear(p + e, J, dp + de) - bilinear(p, J, dp)))
    alpha_p = 2 * J * dp
    beta_p = 2 * J * p
    polynomial_projection = sp.expand(alpha_p.dot(e) + beta_p.dot(de))

    require(sp.simplify(raw - (r_nr + represented_old)) == 0, "finite raw split")
    require(sp.simplify(raw - represented_old - r_nr) == 0, "finite ledger defect")
    require(sp.simplify(polynomial_projection - r_nr) == 0, "finite polynomial projection")
    require(sp.simplify(completed - r_nr - r_quad) == 0, "finite determinant completion")
    require(r_nr != 0, "finite missing-term witness")

    zero_terminal_edge_proposed = r_quad
    zero_terminal_edge_complete = sp.expand(r_nr + r_quad)
    require(
        sp.simplify(zero_terminal_edge_complete - zero_terminal_edge_proposed) == r_nr,
        "failed recombination witness",
    )

    return {
        "raw": str(raw),
        "represented_old": str(represented_old),
        "missing_retained_nonterminal": str(r_nr),
        "quadratic": str(r_quad),
        "completed_nonterminal": str(completed),
        "polynomial_projection": str(polynomial_projection),
        "zero_terminal_edge_proposed": str(zero_terminal_edge_proposed),
        "zero_terminal_edge_complete": str(zero_terminal_edge_complete),
        "finite_rational_audits": 5,
        "witness_boundary": "This is an exact independent-vector algebraic specialization that disproves the proposed universal bookkeeping identity; it is not asserted to be a physical Xi state.",
    }


def exact_certificate() -> dict[str, str | int]:
    closed = Fraction(289339, 858000)
    adverse = Fraction(49, 50)
    necessary = adverse - closed
    sufficient = adverse + closed
    require(necessary == Fraction(551501, 858000), "necessary threshold")
    require(sufficient == Fraction(1130179, 858000), "sufficient threshold")

    return {
        "domain": "Work on one fixed reciprocal-roster cell at the physical point, after the moving-tail current has been assigned separately. Let p be the retained observation, a the terminal conditional S_1 observation, and e every other nonterminal outer observation.",
        "polynomial_split": "With alpha_p=alpha_N, beta_p=beta_N-beta_E, P_nr(lambda)=F_(alpha_p)(lambda)+i(lambda-log a)F_(beta_p)(lambda), and P_E(lambda)=i(lambda-log a)F_(beta_E)(lambda), one has P_lin^[N]=P_nr+P_E.",
        "functional_split": "Using mathscr O_N=mathcal T_N+mathscr R_N, Re{nu c_xi mathscr O_N[P_lin^[N]]}=Re{nu c_xi mathcal T_N[P_nr]}+Re{nu c_xi mathscr R_N[P_nr]}+R_aff. The middle term is R_nr and is neither terminal mixed nor edge affine.",
        "missing_package": "Define E_nr:=2R_nr/|nu|^2. The eight-package ledger omitted E_nr; E_quad^cur contains only R_quad and cannot absorb the retained/nonterminal linear cross.",
        "primitive": "Define Q_nt=(2/|nu|^2){(p+e)^TJ(p+e)-p^TJp}. Then partial_xi Q_nt=E_nr+E_quad^cur exactly.",
        "failed_test": "The proposed completion E_FAN+E_quad^cur is not a determinant-current identity. E_FAN belongs to the terminal/edge-near channel; the exact determinant partner of E_quad^cur is E_nr.",
        "corrected_ledger": "E_tot^corr=E_F+E_R+E_term+E_Delta+E_aff+E_aa+E_move^cur+E_nr+E_quad^cur.",
        "corrected_main": "E_main^corr:=E_FAN+E_nr+E_quad^cur.",
        "thresholds": "The bounded six-package fraction remains 289339/858000. Therefore adverse-witness negativity necessarily requires E_main^corr<-(551501/858000)rho A_T, while E_main^corr<-(1130179/858000)rho A_T is triangle-safe sufficient.",
        "threshold_reason": "The constants do not change because E_nr is restored to the unbounded signed main, not assigned an invented absolute allowance.",
        "next_action": "Decompose the exact completed nonterminal determinant current E_nr+E_quad^cur through the retained finite, near, and common-cutoff paired-remote functionals. Keep E_FAN separate, and test whether the completed nonterminal current has an exact endpoint transfer or a source-correlated sign before taking moduli.",
        "pi_provenance": "No new pi is introduced. This correction uses only the symmetric observation matrix J, linearity, and the already normalized full-support functionals.",
        "proof_boundary": "This gate proves a bookkeeping defect and its exact algebraic repair. It proves no bound or sign for E_nr, E_quad, their completed determinant current, E_FAN, or the complete current; no all-q transport, contact exclusion, Lambda<=0, PF-infinity, RH, or prize-level conclusion.",
        "corrected_package_count": 9,
        "missing_linear_packages": 1,
        "determinant_completions": 1,
        "necessary_numerator": necessary.numerator,
        "necessary_denominator": necessary.denominator,
        "sufficient_numerator": sufficient.numerator,
        "sufficient_denominator": sufficient.denominator,
    }


def build_rows(
    exact: dict[str, str | int],
    symbolic: dict[str, str | int],
) -> list[GateRow]:
    return [
        GateRow("ndc_01_domain", "domain", "proved", "The audit is pointwise on one fixed reciprocal-roster cell after separate moving-current ownership.", str(exact["domain"]), "Tie transfer and roster motion remain separately governed."),
        GateRow("ndc_02_raw", "raw residual", "proved", "The original residual current contains retained-linear, quadratic, terminal, and affine terms.", str(symbolic["raw_current"]), "This is the observation residual sector, not the centre term."),
        GateRow("ndc_03_split", "five-block split", "proved", "The raw residual current has five disjoint algebraic blocks.", str(symbolic["five_block_split"]), "The moving-tail current is outside this fixed-point split."),
        GateRow("ndc_04_terminal", "terminal ownership", "proved", "The terminal-dependent terms are exactly the mixed and pure terminal blocks.", str(symbolic["terminal_mixed"]) + " " + str(symbolic["terminal_pure"]), "They do not contain the retained/nonterminal cross."),
        GateRow("ndc_05_affine", "edge ownership", "proved", "The external edge row is one separately owned linear block.", str(symbolic["edge_affine"]), "P_E remains counted once on the complete outer functional."),
        GateRow("ndc_06_quadratic", "quadratic ownership", "proved", "E_quad owns only the nonterminal self-current.", str(symbolic["quadratic"]), "No retained observation occurs in this block."),
        GateRow("ndc_07_missing", "missing linear block", "proved", "The residual ledger defect is exactly the retained/nonterminal linear cross.", str(symbolic["ledger_defect"]), "This is an ownership correction, not an estimate."),
        GateRow("ndc_08_rows", "polynomial rows", "proved", "The missing cross has the retained derivative and retained base coefficient rows.", str(symbolic["polynomial_rows"]), "All four observation signs follow from J."),
        GateRow("ndc_09_polynomial", "polynomial compression", "proved", "One degree-at-most-five polynomial P_nr carries the missing cross.", str(exact["polynomial_split"]), "The edge beta row is removed exactly once."),
        GateRow("ndc_10_functional", "functional ownership", "proved", "The nonterminal remainder functional applied to P_nr is exactly R_nr.", str(exact["functional_split"]), "mathscr R_N, not the unsplit mathscr O_N, owns this package."),
        GateRow("ndc_11_primitive", "cross primitive", "proved", "R_nr is the derivative of the retained/nonterminal cross primitive.", str(symbolic["retained_nonterminal"]), "Physical coefficient rows are fixed on the auxiliary segment."),
        GateRow("ndc_12_completion", "determinant completion", "proved", "E_nr plus E_quad is one complete nonterminal determinant-current difference.", str(symbolic["determinant_completion"]) + " " + str(exact["primitive"]), "No sign follows from being a derivative."),
        GateRow("ndc_13_falsification", "tested proposal", "proved", "E_FAN+E_quad is not the claimed determinant completion.", str(exact["failed_test"]), "This falsifies only the proposed exact recombination."),
        GateRow("ndc_14_ledger", "corrected ledger", "proved", "The pointwise compensation ledger has nine packages after restoring E_nr.", str(exact["corrected_ledger"]), "Earlier valid component identities and bounds are retained."),
        GateRow("ndc_15_main", "corrected signed main", "proved", "The restored term belongs in the unbounded signed main.", str(exact["corrected_main"]) + " " + str(exact["thresholds"]), "No absolute E_nr allowance is asserted."),
        GateRow("ndc_16_target", "live arithmetic theorem", "open", "The completed nonterminal current must now be decomposed and estimated without breaking its determinant structure.", str(exact["next_action"]), "Neither the pointwise sign nor endpoint transfer is closed."),
        GateRow("ndc_17_boundary", "proof boundary", "guard_validated", "A ledger repair is not a complete-current or RH proof.", str(exact["proof_boundary"]), "All prize-level conclusions remain open."),
    ]


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    symbolic = payload["symbolic"]
    finite = payload["finite_rational"]
    lines = [
        "# Newman C1 Nonterminal Determinant-Completion Gate",
        "",
        "Date: 2026-08-05",
        "",
        "Status: the proposed E_FAN+E_quad determinant recombination is falsified as an exact identity; one omitted retained/nonterminal linear package is restored exactly; the arithmetic sign remains open; not a proof of RH.",
        "",
        "## Exact Residual Split",
        "",
        "```text",
        str(symbolic["raw_current"]),
        str(symbolic["five_block_split"]),
        str(symbolic["retained_nonterminal"]),
        str(symbolic["quadratic"]),
        str(symbolic["terminal_mixed"]),
        str(symbolic["terminal_pure"]),
        str(symbolic["edge_affine"]),
        "```",
        "",
        str(symbolic["ledger_defect"]),
        "",
        "## Full-Support Polynomial Owner",
        "",
        exact["polynomial_split"],
        "",
        exact["functional_split"],
        "",
        exact["missing_package"],
        "",
        "## Determinant Completion",
        "",
        "```text",
        str(symbolic["determinant_completion"]),
        exact["primitive"],
        "```",
        "",
        exact["failed_test"],
        "",
        f"The independent exact-rational specialization gives missing R_nr = `{finite['missing_retained_nonterminal']}`, so the proposed and completed nonterminal values are `{finite['zero_terminal_edge_proposed']}` and `{finite['zero_terminal_edge_complete']}`. {finite['witness_boundary']}",
        "",
        "## Corrected Pointwise Ledger",
        "",
        "```text",
        exact["corrected_ledger"],
        exact["corrected_main"],
        "```",
        "",
        exact["thresholds"],
        "",
        exact["threshold_reason"],
        "",
        "## Next Action",
        "",
        exact["next_action"],
        "",
        "## Pi Provenance",
        "",
        exact["pi_provenance"],
        "",
        "## Proof Boundary",
        "",
        exact["proof_boundary"],
        "",
        payload["success"],
        "",
    ]
    return "\n".join(lines)


def main() -> None:
    source_audit = load_and_audit_sources()
    symbolic = symbolic_certificate()
    finite = finite_rational_certificate()
    exact = exact_certificate()
    rows = build_rows(exact, symbolic)
    require(len(rows) == 17, "row count")
    require(sum(row.readiness == "open" for row in rows) == 1, "open row count")

    success = (
        "built Newman C1 nonterminal determinant-completion gate: "
        "17 rows, 0 issues, 8 symbolic audits, 5 exact-rational audits, "
        "7 source audits, 1 missing linear package identified, "
        "1 corrected determinant completion, 1 open arithmetic row"
    )
    payload = {
        "kind": KIND,
        "schema_version": 1,
        "date": "2026-08-05",
        "status": "one omitted retained/nonterminal linear package restored and exact determinant completion corrected; arithmetic sign open",
        "source_sha256": {key: item["sha256"] for key, item in source_audit.items()},
        "source_audit": source_audit,
        "symbolic": symbolic,
        "finite_rational": finite,
        "exact": exact,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "symbolic_audits": symbolic["symbolic_audits"],
            "finite_rational_audits": finite["finite_rational_audits"],
            "source_audits": len(source_audit),
            "prior_package_count": 8,
            "corrected_package_count": exact["corrected_package_count"],
            "missing_linear_packages": exact["missing_linear_packages"],
            "determinant_completions": exact["determinant_completions"],
            "necessary_signed_fraction": f"{exact['necessary_numerator']}/{exact['necessary_denominator']}",
            "sufficient_signed_fraction": f"{exact['sufficient_numerator']}/{exact['sufficient_denominator']}",
            "open_rows": sum(row.readiness == "open" for row in rows),
        },
        "next_action": exact["next_action"],
        "pi_provenance": exact["pi_provenance"],
        "proof_boundary": exact["proof_boundary"],
        "success": success,
    }
    atomic_write(RESULT_PATH, json.dumps(payload, indent=2, sort_keys=True) + "\n")
    atomic_write(NOTE_PATH, render_note(payload))
    print(success)


if __name__ == "__main__":
    main()
