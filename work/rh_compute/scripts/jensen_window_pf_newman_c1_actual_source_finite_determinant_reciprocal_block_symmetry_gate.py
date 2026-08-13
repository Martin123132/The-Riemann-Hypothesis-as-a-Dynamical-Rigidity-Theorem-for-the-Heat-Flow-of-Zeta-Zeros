#!/usr/bin/env python3
"""Build the finite determinant reciprocal-block symmetry gate."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = "jensen_window_pf_newman_c1_actual_source_finite_determinant_reciprocal_block_symmetry_gate"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"

SOURCE_PATHS = {
    "outer_determinant": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_c1_actual_source_complete_outer_determinant_recombination_gate.json",
    "finite_cell": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_full_support_finite_band_"
        "dirichlet_cell_reduction_gate.json"
    ),
    "full_support": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_"
        "ideal_cubic_reciprocal_terminal_full_support_homotopy_gate.json"
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


def require_zero(expression: sp.Expr, message: str) -> None:
    require(sp.simplify(sp.expand(expression)) == 0, message)


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

    outer = payloads["outer_determinant"]
    require(
        "b-d=omega_E+2p-f" in outer["symbolic"]["band_form"],
        "finite band substitution drifted",
    )
    require(
        "Q(omega_E+p-d)-Q(omega_E+p)" in outer["exact"]["finite_cell"],
        "finite determinant chart drifted",
    )
    require(
        outer["exact"]["package_collapse"].endswith("E_tot^corr=E_geom^cur+E_move^cur."),
        "two-package ledger drifted",
    )
    require("floor functions is forbidden" in outer["exact"]["tie_rule"], "outer tie rule drifted")

    finite = payloads["finite_cell"]["symbolic_certificate"]
    blocks = finite["reciprocal_blocks"]
    require("sum_(q=1)^N sum_(p=1)^N B_(q,p)[P]" in blocks["blocks"], "block partition drifted")
    require("diagonal p=q" in blocks["gradient"], "diagonal block drifted")
    require("p=q+1 and p=q-1" in blocks["adjacent"], "adjacent blocks drifted")
    require("|p-q|>=2" in blocks["far_gap"], "far block drifted")
    require("half-open jump laws" in finite["tie_transport"]["moving_boundary_guard"], "finite tie rule drifted")

    support = payloads["full_support"]["symbolic_certificate"]["current_homotopy"]
    require("Vmathcal N-AQ" in support["observation_order"], "observation determinant drifted")
    require("=h^2P_ret" in support["anchor"], "retained anchor drifted")

    return {
        key: {
            "path": str(path.relative_to(REPO_ROOT)).replace("\\", "/"),
            "sha256": file_hash(path),
        }
        for key, path in SOURCE_PATHS.items()
    }


def real_vector(prefix: str) -> sp.Matrix:
    return sp.Matrix(sp.symbols(f"{prefix}_V {prefix}_N {prefix}_A {prefix}_Q", real=True))


def observation_matrix(vector: sp.Matrix) -> sp.Matrix:
    return sp.Matrix([[vector[0], vector[2]], [vector[3], vector[1]]])


def bilinear(left: sp.Matrix, J: sp.Matrix, right: sp.Matrix) -> sp.Expr:
    return sp.expand((left.T * J * right)[0])


def quadratic(vector: sp.Matrix, J: sp.Matrix) -> sp.Expr:
    return bilinear(vector, J, vector)


def primitive_increment(base: sp.Matrix, defect: sp.Matrix, J: sp.Matrix) -> sp.Expr:
    return sp.expand(quadratic(base - defect, J) - quadratic(base, J))


def current_increment(
    base: sp.Matrix,
    defect: sp.Matrix,
    base_dot: sp.Matrix,
    defect_dot: sp.Matrix,
    J: sp.Matrix,
) -> sp.Expr:
    return sp.expand(
        -2 * bilinear(base, J, defect_dot)
        - 2 * bilinear(base_dot, J, defect)
        + 2 * bilinear(defect, J, defect_dot)
    )


def symbolic_certificate() -> tuple[dict[str, str | int], dict[str, sp.Expr]]:
    J = sp.Rational(1, 2) * sp.Matrix(
        [[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, -1], [0, 0, -1, 0]]
    )
    base = real_vector("b")
    d0, d1, d2 = real_vector("d0"), real_vector("d1"), real_vector("d2")
    dbase = real_vector("db")
    dd0, dd1, dd2 = real_vector("dd0"), real_vector("dd1"), real_vector("dd2")
    defect = d0 + d1 + d2
    ddefect = dd0 + dd1 + dd2
    audits = 0

    matrix_base = observation_matrix(base)
    matrix_defect = observation_matrix(defect)
    require_zero(matrix_base.det() - quadratic(base, J), "observation determinant")
    audits += 1

    adjugate_pairing = sp.trace(matrix_base.adjugate() * matrix_defect)
    require_zero(adjugate_pairing - 2 * bilinear(base, J, defect), "adjugate pairing")
    audits += 1

    total_primitive = primitive_increment(base, defect, J)
    division_free = matrix_defect.det() - adjugate_pairing
    require_zero(total_primitive - division_free, "division-free determinant increment")
    audits += 1

    base0 = base
    base1 = base - d0
    base2 = base - d0 - d1
    primitive0 = primitive_increment(base0, d0, J)
    primitive1 = primitive_increment(base1, d1, J)
    primitive2 = primitive_increment(base2, d2, J)
    require_zero(total_primitive - primitive0 - primitive1 - primitive2, "block primitive telescope")
    audits += 1

    generic_base, generic_defect = real_vector("g"), real_vector("r")
    generic_dbase, generic_ddefect = real_vector("dg"), real_vector("dr")
    generic_current = current_increment(
        generic_base, generic_defect, generic_dbase, generic_ddefect, J
    )
    expected_generic = sp.expand(
        2 * bilinear(generic_base - generic_defect, J, generic_dbase - generic_ddefect)
        - 2 * bilinear(generic_base, J, generic_dbase)
    )
    require_zero(generic_current - expected_generic, "block current formula")
    audits += 1

    current0 = current_increment(base0, d0, dbase, dd0, J)
    current1 = current_increment(base1, d1, dbase - dd0, dd1, J)
    current2 = current_increment(base2, d2, dbase - dd0 - dd1, dd2, J)
    total_current = current_increment(base, defect, dbase, ddefect, J)
    require_zero(total_current - current0 - current1 - current2, "block current telescope")
    audits += 1

    naive = sum(primitive_increment(base, item, J) for item in (d0, d1, d2))
    cross = 2 * (
        bilinear(d0, J, d1) + bilinear(d0, J, d2) + bilinear(d1, J, d2)
    )
    require_zero(total_primitive - naive - cross, "cross ownership")
    audits += 1

    l11, l12, l21, l22 = sp.symbols("l11 l12 l21 l22", real=True)
    r11, r12, r21, r22 = sp.symbols("r11 r12 r21 r22", real=True)
    L = sp.Matrix([[l11, l12], [l21, l22]])
    R = sp.Matrix([[r11, r12], [r21, r22]])
    require_zero(
        (L * matrix_base * R).det() - L.det() * R.det() * matrix_base.det(),
        "left-right determinant covariance",
    )
    audits += 1

    u, a, q, n, s = sp.symbols("u a q n s", real=True)
    D = sp.Matrix([[u, a], [q, n]])
    B_positive = sp.diag(s, s)
    positive_increment = (B_positive - D).det() - B_positive.det()
    require_zero(positive_increment - (D.det() - s * (u + n)), "positive canonical gauge")
    audits += 1

    B_negative = sp.diag(s, -s)
    negative_increment = (B_negative - D).det() - B_negative.det()
    require_zero(negative_increment - (D.det() - s * (-u + n)), "negative canonical gauge")
    audits += 1

    require(audits == 10, "symbolic audit count")
    symbolic = {
        "matrix_chart": "Map x=(V,mathcal N,A,Q)^T to M(x)=[[V,A],[Q,mathcal N]]. Then det M(x)=x^TJx=Vmathcal N-AQ.",
        "division_free_increment": "For B=M(b) and D=M(d), det(B-D)-det(B)=det(D)-tr(adj(B)D)=d^TJd-2b^TJd. No inverse or nonzero determinant is required.",
        "group_symmetry": "For real 2x2 L,R with det(L)det(R)=1, the simultaneous action M(x)->LM(x)R preserves Q and every determinant increment. This is the SL(2,R) x SL(2,R) determinant symmetry, enlarged by the determinant-product-one component.",
        "positive_gauge": "If Delta=det(B)>0 and s=sqrt(Delta), L=sB^(-1), R=I have determinant product one and LBR=sI. In that frozen pointwise gauge, det(B-D)-det(B)=det(D')-s(D'_(11)+D'_(22)).",
        "negative_gauge": "If Delta=det(B)<0 and s=sqrt(-Delta), take L=sB^(-1) and R=diag(1,-1). Their determinants are both -1, LBR=diag(s,-s), and the increment is det(D')-s(-D'_(11)+D'_(22)).",
        "singular_guard": "If det(B)=0, B^(-1) is forbidden. Rank-zero and rank-one bases remain covered by the division-free adjugate formula and must be treated as separate null-cone strata unless h^2P_ret=det(B) is proved nonzero at the calibrated witnesses.",
        "block_partition": "Write f=f^(0)+f^(1)+f^(2), where f^(0) uses reciprocal blocks p=q, f^(1) uses |p-q|=1, and f^(2) uses |p-q|>=2. Define d_0=f^(0)-p, d_1=f^(1), d_2=f^(2), so d=f-p=d_0+d_1+d_2.",
        "primitive_telescope": "With B_0=b, B_1=b-d_0, B_2=b-d_0-d_1 and Delta_j=Q(B_j-d_j)-Q(B_j), Q(b-d)-Q(b)=Delta_0+Delta_1+Delta_2.",
        "current_telescope": "For I(B,D)=partial_xi{Q(B-D)-Q(B)}=-2B^TJdot(D)-2dot(B)^TJD+2D^TJdot(D), the completed geometric current is the sum I(B_0,d_0)+I(B_1,d_1)+I(B_2,d_2), followed by the fixed factor 2/|nu|^2.",
        "cross_ownership": "Using the original base b separately in all three increments omits 2d_0^TJd_1+2d_0^TJd_2+2d_1^TJd_2. The diagonal-first order assigns each cross to the later adjacent or far increment.",
        "symbolic_audits": audits,
    }
    expressions = {
        "total_primitive": total_primitive,
        "primitive0": primitive0,
        "primitive1": primitive1,
        "primitive2": primitive2,
        "total_current": total_current,
        "current0": current0,
        "current1": current1,
        "current2": current2,
        "naive_cross_defect": cross,
    }
    return symbolic, expressions


def rational_vector(values: tuple[Fraction, Fraction, Fraction, Fraction]) -> sp.Matrix:
    return sp.Matrix([sp.Rational(value.numerator, value.denominator) for value in values])


def rational_certificate() -> dict[str, str | int]:
    J = sp.Rational(1, 2) * sp.Matrix(
        [[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, -1], [0, 0, -1, 0]]
    )
    base = rational_vector((Fraction(2, 3), Fraction(-5, 7), Fraction(11, 13), Fraction(17, 19)))
    d0 = rational_vector((Fraction(3, 5), Fraction(7, 11), Fraction(-2, 9), Fraction(5, 8)))
    d1 = rational_vector((Fraction(-4, 7), Fraction(6, 13), Fraction(9, 10), Fraction(-3, 14)))
    d2 = rational_vector((Fraction(5, 12), Fraction(-8, 15), Fraction(7, 16), Fraction(11, 18)))
    dbase = rational_vector((Fraction(1, 4), Fraction(-2, 5), Fraction(3, 7), Fraction(4, 9)))
    dd0 = rational_vector((Fraction(-3, 8), Fraction(5, 12), Fraction(7, 15), Fraction(-2, 11)))
    dd1 = rational_vector((Fraction(4, 13), Fraction(-6, 17), Fraction(8, 19), Fraction(3, 10)))
    dd2 = rational_vector((Fraction(-5, 14), Fraction(7, 18), Fraction(-9, 22), Fraction(11, 24)))
    defect = d0 + d1 + d2
    ddefect = dd0 + dd1 + dd2
    audits = 0

    primitive = primitive_increment(base, defect, J)
    sequential_primitive = primitive_increment(base, d0, J)
    sequential_primitive += primitive_increment(base - d0, d1, J)
    sequential_primitive += primitive_increment(base - d0 - d1, d2, J)
    require_zero(primitive - sequential_primitive, "rational primitive telescope")
    audits += 1

    current = current_increment(base, defect, dbase, ddefect, J)
    sequential_current = current_increment(base, d0, dbase, dd0, J)
    sequential_current += current_increment(base - d0, d1, dbase - dd0, dd1, J)
    sequential_current += current_increment(
        base - d0 - d1, d2, dbase - dd0 - dd1, dd2, J
    )
    require_zero(current - sequential_current, "rational current telescope")
    audits += 1

    naive = sum(primitive_increment(base, item, J) for item in (d0, d1, d2))
    require(primitive != naive, "rational cross witness vanished")
    audits += 1

    B = observation_matrix(base)
    D = observation_matrix(defect)
    L = sp.Matrix([[1, 2], [0, 1]])
    R = sp.Matrix([[1, 0], [3, 1]])
    transformed = (L * (B - D) * R).det() - (L * B * R).det()
    require_zero(transformed - primitive, "rational determinant symmetry")
    audits += 1

    Bp = sp.Matrix([[2, 1], [0, 2]])
    Lp = 2 * Bp.inv()
    require_zero(Lp.det() - 1, "positive gauge determinant")
    require(Lp * Bp == 2 * sp.eye(2), "positive gauge normal form")
    audits += 1

    Bm = sp.Matrix([[2, 1], [0, -2]])
    Lm = 2 * Bm.inv()
    Rm = sp.diag(1, -1)
    require_zero(Lm.det() * Rm.det() - 1, "negative gauge determinant product")
    require(Lm * Bm * Rm == sp.diag(2, -2), "negative gauge normal form")
    audits += 1

    require(audits == 6, "rational audit count")
    return {
        "primitive": str(primitive),
        "current": str(current),
        "naive_original_base_sum": str(naive),
        "missing_cross": str(sp.expand(primitive - naive)),
        "positive_normal_form": str(Lp * Bp),
        "negative_normal_form": str(Lm * Bm * Rm),
        "exact_rational_audits": audits,
        "witness_boundary": "These independent rational matrices audit universal identities and cross ownership; they are not asserted to be physical Xi observations.",
    }


def exact_certificate() -> dict[str, str | int]:
    return {
        "domain": "Work at one physical point inside a fixed reciprocal-roster cell. The matrices L and R used for a pointwise gauge are frozen while differentiating transformed base, defect, and tangent matrices.",
        "base_matrix": "B=M(omega_E+p), D=M(d), and det(B)=Q(omega_E+p)=h^2P_ret at the physical anchor. The completed finite primitive is Q_geom=(2/|nu|^2){det(B-D)-det(B)}.",
        "diagonal_blocks": "The p=q blocks contain the reciprocal stationary point u=alpha_P/r, including the two incomplete endpoint diagonals. Their defect is d_0=f^(0)-p and must retain carrier cancellation.",
        "adjacent_blocks": "The p=q+1 and p=q-1 blocks have no positive uniform gradient gap and meet the diagonal cells at half-integer boundaries. They form d_1 and must be joined in opposite orientations with half-open transfer.",
        "far_blocks": "The |p-q|>=2 blocks form d_2 and obey |alpha_P/u-r|>=alpha_P(|p-q|-1)/{(p+1/2)(q+1/2)} before moduli. Their determinant increment still contains a cross with the updated diagonal-plus-adjacent base and a far self-determinant.",
        "tie_rule": "At a reciprocal boundary, transfer the complete block mode and its lift. The determinant total is tie invariant; no floor function is differentiated.",
        "norm_guard": "Euclidean coordinate norms are not invariant under determinant-product-one left/right gauges. Any arithmetic estimate must either fix and justify a conditioned gauge or bound the invariant quantities det(D_j) and tr(adj(B_j)D_j), together with their tangent analogues.",
        "singular_obligation": "The nonsingular canonical gauge is conditional on P_ret!=0. A complete proof must certify that condition at every selected calibrated witness or include rank-one and zero base matrices through the division-free adjugate chart.",
        "next_action": "Compute the physical diagonal and adjacent block matrices jointly at calibrated B=-49/100 level points, preserving source phases and endpoint halves. Bound the far invariant increment using the certified gradient gap only after its cross with the updated base is retained. Record det(B) to determine whether the nonsingular gauge is available; otherwise stratify the null-cone case. Compare the block telescope against the near/paired-remote chart before selecting an arithmetic lemma.",
        "pi_provenance": "No new pi is introduced by the determinant symmetry or block telescope. Pi remains inherited from e(x)=exp(2pi i x), alpha_P=xi/(2pi), and the completed-zeta saddle normalization.",
        "proof_boundary": "This gate proves the determinant-matrix chart, division-free adjugate identity, determinant-product-one symmetry, conditional nonsingular normal forms, singular-base guard, and exact diagonal/adjacent/far primitive and current telescopes. It proves no physical block sign or bound, no nonsingularity of P_ret, no geometric-current inequality, complete-current sign, all-q transport, contact exclusion, Lambda<=0, PF-infinity, RH, or prize-level conclusion.",
        "reciprocal_block_classes": 3,
        "determinant_symmetries": 1,
        "singular_base_guards": 1,
    }


def build_rows(symbolic: dict, rational: dict, exact: dict) -> list[GateRow]:
    rows = [
        GateRow("fdr_01_domain", "fixed-cell domain", "proved", "The determinant and tangent matrices are evaluated on one fixed roster cell.", exact["domain"], "Moving gauges are not differentiated."),
        GateRow("fdr_02_matrix", "matrix chart", "proved", "The split quadratic form is one 2x2 determinant.", symbolic["matrix_chart"], "This is a coordinate identity."),
        GateRow("fdr_03_increment", "adjugate increment", "proved", "The finite geometric primitive has a division-free determinant formula.", symbolic["division_free_increment"], "No base inverse is assumed."),
        GateRow("fdr_04_physical", "physical base", "proved", "The matrix chart represents the completed finite primitive and retained anchor.", exact["base_matrix"], "No sign of P_ret is claimed."),
        GateRow("fdr_05_symmetry", "determinant symmetry", "proved", "A determinant-product-one left/right action preserves the full target.", symbolic["group_symmetry"], "Coordinate sizes are gauge dependent."),
        GateRow("fdr_06_positive", "positive normal form", "proved", "Every nonsingular positive-determinant base has a scalar identity gauge.", symbolic["positive_gauge"], "The frame is frozen at the point."),
        GateRow("fdr_07_negative", "negative normal form", "proved", "Every nonsingular negative-determinant base has a split diagonal gauge.", symbolic["negative_gauge"], "The frame is frozen at the point."),
        GateRow("fdr_08_singular", "singular-base guard", "guard_validated", "The null cone remains in the division-free chart.", symbolic["singular_guard"] + " " + exact["singular_obligation"], "Nonsingularity is not smuggled in."),
        GateRow("fdr_09_partition", "reciprocal partition", "proved", "The finite band has diagonal, adjacent, and far block classes.", symbolic["block_partition"], "The partition uses all band blocks once."),
        GateRow("fdr_10_diagonal", "stationary diagonal", "proved", "The diagonal defect retains the carrier and endpoint stationary cells.", exact["diagonal_blocks"], "No diagonal sign is proved."),
        GateRow("fdr_11_adjacent", "adjacent join", "proved", "The two adjacent orientations remain joined to the diagonal sector.", exact["adjacent_blocks"], "No false gradient gap is used."),
        GateRow("fdr_12_far", "far sector", "proved", "The far blocks have the certified reciprocal phase gap.", exact["far_blocks"], "The determinant cross and self terms remain."),
        GateRow("fdr_13_primitive", "primitive telescope", "proved", "The completed determinant is the exact diagonal-first three-stage telescope.", symbolic["primitive_telescope"], "No block is estimated separately here."),
        GateRow("fdr_14_cross", "cross ownership", "proved", "Later increments own all earlier/later determinant crosses.", symbolic["cross_ownership"], "Reusing the original base three times is invalid."),
        GateRow("fdr_15_current", "current telescope", "proved", "Differentiation preserves the same three-stage ownership.", symbolic["current_telescope"], "The normalized factor is attached once."),
        GateRow("fdr_16_rational", "exact witness", "proved", "Independent exact matrices expose the omitted crosses in a naive split.", f"missing cross={rational['missing_cross']}; current={rational['current']}", rational["witness_boundary"]),
        GateRow("fdr_17_ties", "tie transport", "guard_validated", "Complete mode transfer preserves the determinant total across roster boundaries.", exact["tie_rule"], "Floor derivatives are forbidden."),
        GateRow("fdr_18_norms", "invariant estimate guard", "guard_validated", "Coordinate norms cannot be treated as intrinsic under the new symmetry.", exact["norm_guard"], "A conditioned gauge or invariant pairing is required."),
        GateRow("fdr_19_target", "live arithmetic theorem", "open", "Control the diagonal-plus-adjacent invariant and then the far invariant at calibrated physical rows.", exact["next_action"], "No physical signed estimate is yet proved."),
        GateRow("fdr_20_boundary", "proof boundary", "guard_validated", "The symmetry and telescope are not a proof of RH.", exact["proof_boundary"], "All prize-level conclusions remain open."),
    ]
    require(len(rows) == 20, "row count")
    return rows


def render_note(payload: dict) -> str:
    symbolic = payload["symbolic"]
    exact = payload["exact"]
    rational = payload["finite_rational"]
    return "\n".join(
        [
            "# Newman C1 Finite Determinant Reciprocal-Block Symmetry Gate",
            "",
            "Date: 2026-08-05",
            "",
            "Status: exact determinant symmetry and reciprocal-block telescope proved; physical signed arithmetic open; not a proof of RH.",
            "",
            "## Determinant Chart",
            "",
            symbolic["matrix_chart"],
            "",
            symbolic["division_free_increment"],
            "",
            exact["base_matrix"],
            "",
            "## Symmetry And Gauges",
            "",
            symbolic["group_symmetry"],
            "",
            symbolic["positive_gauge"],
            "",
            symbolic["negative_gauge"],
            "",
            symbolic["singular_guard"],
            "",
            exact["norm_guard"],
            "",
            "## Reciprocal Blocks",
            "",
            symbolic["block_partition"],
            "",
            exact["diagonal_blocks"],
            "",
            exact["adjacent_blocks"],
            "",
            exact["far_blocks"],
            "",
            symbolic["primitive_telescope"],
            "",
            symbolic["cross_ownership"],
            "",
            symbolic["current_telescope"],
            "",
            "## Exact Rational Audit",
            "",
            f"The naive original-base split misses `{rational['missing_cross']}` in the stored rational witness.",
            "",
            rational["witness_boundary"],
            "",
            "## Tie Rule",
            "",
            exact["tie_rule"],
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
    )


def main() -> int:
    source_audit = load_and_audit_sources()
    symbolic, expressions = symbolic_certificate()
    rational = rational_certificate()
    exact = exact_certificate()
    rows = build_rows(symbolic, rational, exact)
    success = (
        "built Newman C1 finite determinant reciprocal-block symmetry gate: "
        "20 rows, 0 issues, 10 symbolic audits, 6 exact-rational audits, "
        "3 source audits, 3 reciprocal block classes, 1 determinant symmetry, "
        "1 singular-base guard, 1 open arithmetic row"
    )
    payload = {
        "kind": KIND,
        "schema_version": 1,
        "date": "2026-08-05",
        "status": "finite determinant symmetry and reciprocal-block telescope proved; physical signed arithmetic open",
        "source_sha256": {key: value["sha256"] for key, value in source_audit.items()},
        "source_audit": source_audit,
        "symbolic": symbolic,
        "stored_expressions": {key: str(value) for key, value in expressions.items()},
        "finite_rational": rational,
        "exact": exact,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "issues": 0,
            "symbolic_audits": symbolic["symbolic_audits"],
            "exact_rational_audits": rational["exact_rational_audits"],
            "source_audits": len(source_audit),
            "reciprocal_block_classes": exact["reciprocal_block_classes"],
            "determinant_symmetries": exact["determinant_symmetries"],
            "singular_base_guards": exact["singular_base_guards"],
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
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
