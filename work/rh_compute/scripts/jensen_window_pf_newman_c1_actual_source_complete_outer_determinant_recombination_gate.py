#!/usr/bin/env python3
"""Build the complete outer-determinant recombination gate."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = "jensen_window_pf_newman_c1_actual_source_complete_outer_determinant_recombination_gate"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"

SOURCE_PATHS = {
    "nonterminal_completion": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_c1_actual_source_nonterminal_determinant_completion_gate.json",
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
    "edge_affine": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_c1_actual_source_edge_affine_ownership_gate.json",
    "compensation_ledger": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_c1_actual_source_complete_compensation_ledger_gate.json",
    "current_transfer": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_c1_actual_source_current_transfer_quadratic_primitive_gate.json",
    "raw_center": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_c1_actual_source_raw_center_sign_obstruction_gate.json",
    "full_support": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_"
        "ideal_cubic_reciprocal_terminal_full_support_homotopy_gate.json"
    ),
    "joined_recombination": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_full_support_joined_lower_"
        "interior_abel_recombination_gate.json"
    ),
    "symmetric_outer": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_full_support_ideal_cubic_"
        "symmetric_outer_pairing_gate.json"
    ),
    "finite_cell": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_full_support_finite_band_"
        "dirichlet_cell_reduction_gate.json"
    ),
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

    nonterminal = payloads["nonterminal_completion"]["exact"]
    require(
        nonterminal["corrected_ledger"].startswith("E_tot^corr=E_F+E_R"),
        "corrected pointwise ledger drifted",
    )
    require(
        "partial_xi Q_nt=E_nr+E_quad^cur" in nonterminal["primitive"],
        "nonterminal primitive drifted",
    )

    terminal = payloads["terminal_completion"]["symbolic_certificate"]["coefficient_completion"]
    require("C(a+rho)=C(rho)" in terminal["identity"], "terminal split drifted")
    require("Delta_a=2omega_*^TJdot a" in terminal["complete_rows"], "terminal current drifted")

    relative = payloads["relative_lift"]["symbolic_certificate"]
    require(
        "mathscr O_N[P]=T_N[P]+mathscr R_N[P]" in relative["outer_remainder"]["split"],
        "outer split drifted",
    )
    require(
        "actual mixed real current" in relative["mixed_current"]["projection"],
        "mixed-current projection drifted",
    )

    edge = payloads["edge_affine"]["exact"]
    require("t_E=2Jomega_E" in edge["row_split"], "edge row drifted")
    require("mathscr O_N[P_E]" in edge["edge_polynomial"], "edge functional drifted")

    ledger = payloads["compensation_ledger"]["exact"]
    require("E_W=E_F+E_R+E_term" in ledger["ideal_partition"]["projection"], "ideal partition drifted")
    require("linear actual-source projection" in ledger["physical_partition"]["correction"], "physical correction drifted")

    current = payloads["current_transfer"]["exact"]
    require(current["current_budget"] == "|E_move^cur|<rho A_T/100.", "moving-current budget drifted")

    center = payloads["raw_center"]["exact"]["sign_obstruction"]
    require(">49/100" in center["positive_point"], "positive center witness drifted")
    require("<-49/100" in center["negative_point"], "negative center witness drifted")
    require("roster interiors" in center["tie_guard"], "center tie guard drifted")

    full_support = payloads["full_support"]["symbolic_certificate"]
    require("G_j^[N]=G_j^[B]+G_j^(C)" in full_support["full_support_observations"]["moments"], "full-support carrier drifted")
    require("M_N[A]-M_B[A]" in full_support["terminal_extension"]["starred_sum"], "starred carrier drifted")

    joined = payloads["joined_recombination"]["symbolic_certificate"]["reverse_recombination"]
    require("R_N[P]=-D_N[P]-T_N[P]" in joined["objects"], "R/D/T identity drifted")
    require("R_N[P]=O_N[P]-T_N[P]" in joined["outer"], "R/O/T identity drifted")

    outer = payloads["symmetric_outer"]["symbolic_certificate"]
    require("mathscr O_N=E_band-tau_band=-mathscr D_N" in outer["recombination"]["exterior"], "outer/defect identity drifted")
    require("Define mathscr N_N[P]" in outer["cutoff"]["near"], "near split drifted")
    require("paired series mathscr C_N[P]" in outer["cutoff"]["remote"], "remote split drifted")

    finite = payloads["finite_cell"]["symbolic_certificate"]
    require("mathscr D_N[P]=-mathscr M_N[P]" in finite["cell_reduction"]["defect"], "finite defect drifted")
    require("half-open jump laws" in finite["tie_transport"]["moving_boundary_guard"], "tie law drifted")

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


def increment_current(
    base: sp.Matrix,
    increment: sp.Matrix,
    base_dot: sp.Matrix,
    increment_dot: sp.Matrix,
    J: sp.Matrix,
) -> sp.Expr:
    return sp.expand(
        2 * bilinear(base, J, increment_dot)
        + 2 * bilinear(base_dot, J, increment)
        + 2 * bilinear(increment, J, increment_dot)
    )


def symbolic_certificate() -> dict[str, str | int]:
    J = sp.Rational(1, 2) * sp.Matrix(
        [
            [0, 1, 0, 0],
            [1, 0, 0, 0],
            [0, 0, 0, -1],
            [0, 0, -1, 0],
        ]
    )
    p, a, o, near, remote, edge = (
        real_vector(name) for name in ("p", "a", "o", "near", "remote", "edge")
    )
    dp, da, do, dnear, dremote = (
        real_vector(name) for name in ("dp", "da", "do", "dnear", "dremote")
    )
    e = o - a
    de = do - da

    r_nr = 2 * bilinear(p, J, de) + 2 * bilinear(dp, J, e)
    r_quad = 2 * bilinear(e, J, de)
    r_mixed = 2 * bilinear(p + e, J, da) + 2 * bilinear(dp + de, J, a)
    r_aa = 2 * bilinear(a, J, da)
    residual_outer = increment_current(p, o, dp, do, J)

    audits = 0
    require(sp.simplify(a + e - o) == sp.zeros(4, 1), "terminal plus nonterminal")
    audits += 1
    require(sp.simplify(r_nr + r_quad + r_mixed + r_aa - residual_outer) == 0, "outer determinant recombination")
    audits += 1

    t_edge = 2 * J * edge
    affine = sp.expand(t_edge.dot(do))
    require(sp.simplify(affine - 2 * bilinear(edge, J, do)) == 0, "edge affine")
    audits += 1

    base = edge + p
    geometric = sp.expand(residual_outer + affine)
    shifted = increment_current(base, o, dp, do, J)
    require(sp.simplify(geometric - shifted) == 0, "edge-shifted determinant")
    audits += 1

    primitive_derivative = 2 * (
        bilinear(base + o, J, dp + do) - bilinear(base, J, dp)
    )
    require(sp.simplify(geometric - primitive_derivative) == 0, "geometric primitive")
    audits += 1

    near_remote_substitution = {o[index]: near[index] + remote[index] for index in range(4)}
    near_remote_substitution.update(
        {do[index]: dnear[index] + dremote[index] for index in range(4)}
    )
    geometric_nr = sp.expand(geometric.subs(near_remote_substitution, simultaneous=True))
    near_current = increment_current(base, near, dp, dnear, J)
    remote_current = increment_current(base + near, remote, dp + dnear, dremote, J)
    require(sp.simplify(geometric_nr - near_current - remote_current) == 0, "near/remote determinant split")
    audits += 1

    defect, ddefect = real_vector("defect"), real_vector("ddefect")
    finite_substitution = {o[index]: -defect[index] for index in range(4)}
    finite_substitution.update({do[index]: -ddefect[index] for index in range(4)})
    geometric_finite = sp.expand(geometric.subs(finite_substitution, simultaneous=True))
    finite_current = increment_current(base, -defect, dp, -ddefect, J)
    require(sp.simplify(geometric_finite - finite_current) == 0, "finite-defect current")
    audits += 1

    band, dband = real_vector("band"), real_vector("dband")
    defect_to_band = {defect[index]: band[index] - p[index] for index in range(4)}
    defect_to_band.update(
        {ddefect[index]: dband[index] - dp[index] for index in range(4)}
    )
    finite_band = sp.expand(finite_current.subs(defect_to_band, simultaneous=True))
    finite_expected = 2 * (
        bilinear(edge + 2 * p - band, J, 2 * dp - dband)
        - bilinear(edge + p, J, dp)
    )
    require(sp.simplify(finite_band - finite_expected) == 0, "band/carrier primitive")
    audits += 1

    reverse_near = increment_current(base, remote, dp, dremote, J)
    reverse_near += increment_current(base + remote, near, dp + dremote, dnear, J)
    require(sp.simplify(near_current + remote_current - reverse_near) == 0, "allocation order total")
    audits += 1

    return {
        "observation_owners": "p_X=Re{nu c_xi mathscr M_N[P_X]}, a_X=Re{nu c_xi mathcal T_N[P_X]}, e_X=Re{nu c_xi mathscr R_N[P_X]}, and o_X=a_X+e_X=Re{nu c_xi mathscr O_N[P_X]}.",
        "terminal_cancellation": "Because mathscr O_N=mathcal T_N+mathscr R_N, one has a+e=o exactly before any modulus.",
        "outer_current": "R_outer=2p^TJdot(o)+2dot(p)^TJo+2o^TJdot(o)=partial_xi{(p+o)^TJ(p+o)-p^TJp}.",
        "package_recombination": "R_nr+R_quad+R_(a*)+R_aa=R_outer.",
        "edge_shift": "With b=omega_E+p and t_E=2Jomega_E, R_geom=R_outer+R_aff=partial_xi{(b+o)^TJ(b+o)-b^TJb}.",
        "near_remote": "For o=n+c with n the finite near observation and c the common-cutoff paired-remote observation, R_geom=partial_xi{Q(b+n)-Q(b)}+partial_xi{Q(b+n+c)-Q(b+n)}.",
        "canonical_order": "The near-first order assigns the near/remote cross to the remote increment, matching the certified rule that the finite mean-one near kernel is composed first and the paired remote tail is attached afterward. Reversing the order changes only cross ownership, not the total.",
        "finite_defect": "Since mathscr O_N=-mathscr D_N, write d_X=Re{nu c_xi mathscr D_N[P_X]}; then R_geom=partial_xi{Q(b-d)-Q(b)}.",
        "band_form": "Because mathscr D_N=mathscr F_N-mathscr M_N, b-d=omega_E+2p-f, where f_X=Re{nu c_xi mathscr F_N[P_X]}. Thus the completed geometric current has a purely finite reciprocal-band representation.",
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
    a = rational_vector((Fraction(1, 3), Fraction(2, 5), Fraction(-3, 7), Fraction(4, 9)))
    da = rational_vector((Fraction(-5, 11), Fraction(7, 13), Fraction(-11, 17), Fraction(13, 19)))
    near = rational_vector((Fraction(-2, 5), Fraction(3, 7), Fraction(-5, 9), Fraction(7, 11)))
    dnear = rational_vector((Fraction(11, 13), Fraction(-13, 17), Fraction(17, 19), Fraction(-19, 23)))
    remote = rational_vector((Fraction(3, 11), Fraction(-5, 13), Fraction(7, 17), Fraction(-11, 19)))
    dremote = rational_vector((Fraction(-13, 23), Fraction(17, 29), Fraction(-19, 31), Fraction(23, 37)))
    edge = rational_vector((Fraction(2, 3), Fraction(-3, 5), Fraction(5, 7), Fraction(-7, 11)))
    o, do = near + remote, dnear + dremote
    e, de = o - a, do - da
    base = edge + p

    r_nr = 2 * bilinear(p, J, de) + 2 * bilinear(dp, J, e)
    r_quad = 2 * bilinear(e, J, de)
    r_mixed = 2 * bilinear(p + e, J, da) + 2 * bilinear(dp + de, J, a)
    r_aa = 2 * bilinear(a, J, da)
    outer = increment_current(p, o, dp, do, J)
    affine = 2 * bilinear(edge, J, do)
    geometric = sp.expand(outer + affine)
    shifted = increment_current(base, o, dp, do, J)
    near_current = increment_current(base, near, dp, dnear, J)
    remote_current = increment_current(base + near, remote, dp + dnear, dremote, J)
    defect, ddefect = -o, -do
    finite_current = increment_current(base, -defect, dp, -ddefect, J)
    band, dband = p + defect, dp + ddefect
    band_current = 2 * (
        bilinear(edge + 2 * p - band, J, 2 * dp - dband)
        - bilinear(edge + p, J, dp)
    )

    require(sp.simplify(r_nr + r_quad + r_mixed + r_aa - outer) == 0, "finite package recombination")
    require(sp.simplify(geometric - shifted) == 0, "finite edge shift")
    require(sp.simplify(geometric - near_current - remote_current) == 0, "finite near remote")
    require(sp.simplify(geometric - finite_current) == 0, "finite defect")
    require(sp.simplify(geometric - band_current) == 0, "finite band")
    require(geometric != 0, "nonzero finite witness")

    return {
        "outer": str(outer),
        "affine": str(affine),
        "geometric": str(geometric),
        "near_current": str(near_current),
        "remote_current": str(remote_current),
        "finite_defect_current": str(finite_current),
        "finite_band_current": str(band_current),
        "finite_rational_audits": 6,
        "witness_boundary": "This exact rational-vector audit checks the identities independently of the builder's symbolic names; it is not asserted to be a physical Xi state.",
    }


def exact_certificate() -> dict[str, str | int]:
    calibrated_center = Fraction(49, 50)
    moving = Fraction(1, 100)
    necessary = calibrated_center - moving
    sufficient = calibrated_center + moving
    require(necessary == Fraction(97, 100), "necessary geometric threshold")
    require(sufficient == Fraction(99, 100), "sufficient geometric threshold")

    return {
        "domain": "Work pointwise on one fixed reciprocal-roster cell with physical coefficient rows fixed and the genuine edge omega_E fixed. Roster ties are handled by complete mode transfer rather than derivatives of floor functions.",
        "terminal_package": "The normalized terminal mixed current is E_(a*)=E_F+E_R+E_term+E_Delta, and E_aa is the pure terminal current. Therefore E_(a*)+E_aa+E_nr+E_quad^cur is the normalized derivative of the complete outer determinant difference based on o=a+e.",
        "geometric_package": "Define Q_geom=(2/|nu|^2){(omega_E+p+o)^TJ(omega_E+p+o)-(omega_E+p)^TJ(omega_E+p)} and E_geom^cur=partial_xi Q_geom.",
        "package_collapse": "E_geom^cur=E_F+E_R+E_term+E_Delta+E_aa+E_nr+E_quad^cur+E_aff. Hence E_tot^corr=E_geom^cur+E_move^cur.",
        "near_primitive": "Q_near=(2/|nu|^2){Q(omega_E+p+n)-Q(omega_E+p)}.",
        "remote_primitive": "Q_remote=(2/|nu|^2){Q(omega_E+p+n+c)-Q(omega_E+p+n)}. It contains the near/remote cross and remote self-determinant, so the existing linear remote budget cannot be reused for it without a new proof.",
        "near_remote_completion": "Q_geom=Q_near+Q_remote exactly for mathscr O_N=mathscr N_N+mathscr C_N, with mathscr C_N defined only as a common-cutoff paired-remote observation.",
        "finite_cell": "With d_X=Re{nu c_xi mathscr D_N[P_X]} and mathscr O_N=-mathscr D_N, Q_geom=(2/|nu|^2){Q(omega_E+p-d)-Q(omega_E+p)}. Since mathscr D_N=-mathscr M_N+H_N{f_P(1)-f_P(N)}+V_N[P], this contains no conditional infinite series.",
        "tie_rule": "At either reciprocal boundary, transfer the complete mode and its lift between mathscr F_N and mathscr O_N. The full observation and Q_geom are tie invariant; differentiating m_N or n_N as floor functions is forbidden.",
        "current_criterion": "Sigma_full=-2rho A_TB+E_geom^cur+E_move^cur. On a roster-interior continuous arc carrying the certified values B>49/100 and B<-49/100, the intermediate-value theorem supplies a calibrated point B=-49/100. At that point, |E_move^cur|<rho A_T/100 makes E_geom^cur<-(97/100)rho A_T necessary for negativity, while E_geom^cur<-(99/100)rho A_T is sufficient independently of the moving-current sign.",
        "threshold_interpretation": "These 97/100 and 99/100 constants replace the prior signed-main bookkeeping only for the newly completed geometric package. The previously proved component budgets remain valid fallbacks but are not spent in this joined determinant formulation.",
        "next_action": "Choose between two equivalent arithmetic charts. In the finite chart, insert the Dirichlet-cell formula for mathscr D_N and derive the reciprocal-stationary decomposition of the determinant Q(omega_E+p-d)-Q(omega_E+p). In the near/remote chart, compose the mean-one finite near kernel with the retained carrier before bounding the paired-remote cross and self-determinant. Test both charts on the physical coefficient rows before selecting Type-I/II, Vaughan, adjacent-recurrence, or reciprocal-pair machinery.",
        "pi_provenance": "No new pi is introduced. The determinant recombination uses only J and linearity. Pi in the finite and paired-remote representations is inherited from e(x)=exp(2pi i x), kappa=1/(2pi i), and the finite Dirichlet kernel.",
        "proof_boundary": "This gate proves exact package recombination, edge completion, canonical near/paired-remote and finite-cell determinant forms, tie ownership, and geometric-current threshold arithmetic. It proves no sign or bound for E_geom, Q_near, Q_remote, the finite-cell determinant, or the complete current; no all-q transport, contact exclusion, Lambda<=0, PF-infinity, RH, or prize-level conclusion.",
        "collapsed_packages": 8,
        "remaining_pointwise_packages": 2,
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
        GateRow("odr_01_domain", "domain", "proved", "The recombination is fixed-cell and tie-complete.", str(exact["domain"]), "No moving floor derivative is used."),
        GateRow("odr_02_owners", "source owners", "proved", "The retained, terminal, nonterminal, and complete outer observations have exact functional owners.", str(symbolic["observation_owners"]), "All source phases remain fixed."),
        GateRow("odr_03_terminal", "terminal cancellation", "proved", "The terminal conditional and nonterminal remainder reconstruct the complete outer functional.", str(symbolic["terminal_cancellation"]), "This is functional addition, not a numerical cancellation estimate."),
        GateRow("odr_04_outer", "outer determinant", "proved", "All terminal and nonterminal residual-current blocks form one outer determinant derivative.", str(symbolic["outer_current"]), "The external edge is not yet included in this row."),
        GateRow("odr_05_packages", "package recombination", "proved", "The restored nonterminal blocks and terminal blocks exhaust the outer determinant current.", str(symbolic["package_recombination"]) + " " + str(exact["terminal_package"]), "Previously valid package identities are retained."),
        GateRow("odr_06_edge", "edge completion", "proved", "The genuine-edge affine row shifts the determinant base by omega_E.", str(symbolic["edge_shift"]), "omega_E is fixed on the auxiliary segment."),
        GateRow("odr_07_geom", "geometric primitive", "proved", "Eight pointwise packages are one exact geometric determinant current.", str(exact["geometric_package"]), "A derivative identity supplies no sign."),
        GateRow("odr_08_collapse", "ledger collapse", "proved", "The corrected nine-package ledger reduces to geometric current plus moving current.", str(exact["package_collapse"]), "E_move remains separately owned."),
        GateRow("odr_09_outer_split", "canonical outer split", "proved", "The complete outer functional is finite near plus common-cutoff paired remote.", str(exact["near_remote_completion"]), "One-sided remote infinities remain inadmissible."),
        GateRow("odr_10_near", "near determinant", "proved", "The retained carrier and mean-one near kernel form the first determinant increment.", str(exact["near_primitive"]), "No near sign is proved."),
        GateRow("odr_11_remote", "remote determinant", "proved", "The paired remote tail forms the second determinant increment after the near join.", str(exact["remote_primitive"]), "The old linear remote bound does not close this nonlinear increment."),
        GateRow("odr_12_order", "cross ownership", "proved", "Near-first ordering assigns the cross term canonically without changing the total.", str(symbolic["canonical_order"]), "This is a bookkeeping convention tied to the certified outer split."),
        GateRow("odr_13_defect", "finite-cell form", "proved", "The same geometric determinant is a returned-carrier minus finite-cell-defect increment.", str(symbolic["finite_defect"]) + " " + str(exact["finite_cell"]), "The finite defect is not known small."),
        GateRow("odr_14_band", "finite band", "proved", "The completed determinant has a representation containing no conditional infinite series.", str(symbolic["band_form"]), "Reciprocal-stationary cancellation remains to be proved."),
        GateRow("odr_15_ties", "tie transport", "guard_validated", "Complete mode transfer preserves the geometric primitive across roster boundaries.", str(exact["tie_rule"]), "Cellwise derivatives must be patched by jump laws."),
        GateRow("odr_16_threshold", "signed criterion", "proved", "At the calibrated adverse level, only the moving-current allowance is detached from the geometric target.", str(exact["current_criterion"]) + " " + str(exact["threshold_interpretation"]), "Neither geometric threshold inequality is proved."),
        GateRow("odr_17_target", "live arithmetic theorem", "open", "A physical signed estimate for one of the two equivalent determinant charts remains open.", str(exact["next_action"]), "No chart has yet supplied the 99/100 sufficient sign."),
        GateRow("odr_18_boundary", "proof boundary", "guard_validated", "Exact determinant completion is not a complete-current or RH proof.", str(exact["proof_boundary"]), "All prize-level conclusions remain open."),
    ]


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    symbolic = payload["symbolic"]
    finite = payload["finite_rational"]
    lines = [
        "# Newman C1 Complete Outer-Determinant Recombination Gate",
        "",
        "Date: 2026-08-05",
        "",
        "Status: eight residual and affine packages recombined into one geometric determinant current; moving current remains separate; signed arithmetic open; not a proof of RH.",
        "",
        "## Observation Owners",
        "",
        symbolic["observation_owners"],
        "",
        symbolic["terminal_cancellation"],
        "",
        "## Complete Determinant",
        "",
        "```text",
        symbolic["outer_current"],
        symbolic["package_recombination"],
        symbolic["edge_shift"],
        exact["geometric_package"],
        exact["package_collapse"],
        "```",
        "",
        f"The independent exact-rational audit gives outer `{finite['outer']}`, affine `{finite['affine']}`, and geometric total `{finite['geometric']}`. {finite['witness_boundary']}",
        "",
        "## Near And Paired Remote",
        "",
        "```text",
        exact["near_primitive"],
        exact["remote_primitive"],
        exact["near_remote_completion"],
        "```",
        "",
        symbolic["canonical_order"],
        "",
        "## Finite-Cell Form",
        "",
        symbolic["finite_defect"],
        "",
        symbolic["band_form"],
        "",
        exact["finite_cell"],
        "",
        exact["tie_rule"],
        "",
        "## Corrected Signed Criterion",
        "",
        exact["current_criterion"],
        "",
        exact["threshold_interpretation"],
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
    require(len(rows) == 18, "row count")
    require(sum(row.readiness == "open" for row in rows) == 1, "open row count")

    success = (
        "built Newman C1 complete outer-determinant recombination gate: "
        "18 rows, 0 issues, 9 symbolic audits, 6 exact-rational audits, "
        "11 source audits, 8 packages collapsed, 2 pointwise packages remain, "
        "1 finite-cell determinant form, 1 open arithmetic row"
    )
    payload = {
        "kind": KIND,
        "schema_version": 1,
        "date": "2026-08-05",
        "status": "complete outer and edge determinant recombination proved; moving current and signed arithmetic remain open",
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
            "collapsed_packages": exact["collapsed_packages"],
            "remaining_pointwise_packages": exact["remaining_pointwise_packages"],
            "near_remote_determinant_splits": 1,
            "finite_cell_determinant_forms": 1,
            "necessary_geometric_fraction": f"{exact['necessary_numerator']}/{exact['necessary_denominator']}",
            "sufficient_geometric_fraction": f"{exact['sufficient_numerator']}/{exact['sufficient_denominator']}",
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
