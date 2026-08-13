#!/usr/bin/env python3
"""Build the fixed-root shift-tangency rigidity and Hankel bridge gate."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = "jensen_window_pf_fixed_root_shift_tangency_rigidity_gate"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"

SOURCE_PATHS = {
    "heat_hierarchy": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_heat_flow_jensen_hierarchy_lemma.json",
    "polar_contact": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_quartic_quintic_polar_contact_lemma.json",
    "coefficient_pf": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_coefficient_pf_equivalence_gate.json",
    "degree_cascade": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_polar_heat_collision_cascade_lemma.json",
    "order4_m100": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_compound_order4_m100_entry_certificate.json",
    "order4_forward": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_compound_order4_uniform_heat_forward_invariance_certificate.json",
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


def require_zero(expression: sp.Expr, label: str) -> None:
    value = sp.simplify(sp.factor(sp.together(expression)))
    if value != 0 and value.equals(0) is not True:
        raise RuntimeError(f"{label}: {value}")


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(text, encoding="utf-8")
    os.replace(temporary, path)


def load_sources() -> dict[str, dict]:
    payloads: dict[str, dict] = {}
    for key, path in SOURCE_PATHS.items():
        require(path.exists(), f"missing source artifact: {path}")
        payloads[key] = json.loads(path.read_text(encoding="utf-8"))
    require(len(payloads["heat_hierarchy"].get("rows", [])) == 9, "heat hierarchy drift")
    require(len(payloads["polar_contact"].get("rows", [])) == 10, "polar contact drift")
    require(len(payloads["coefficient_pf"].get("rows", [])) == 13, "coefficient-PF drift")
    require(len(payloads["degree_cascade"].get("rows", [])) == 10, "degree cascade drift")
    require(
        payloads["order4_m100"].get("exact", {}).get("all_shift_entry")
        == "H_(4,n)(-100)>0, every n>=0",
        "order-four entry drift",
    )
    require(len(payloads["order4_forward"].get("rows", [])) == 9, "order-four forward drift")
    require(
        payloads["order4_forward"].get("exact", {}).get("all_interval_theorem")
        == "H_(4,n)(lambda)>0 for every integer n>=0 and every lambda in [-100,0]",
        "order-four interval theorem drift",
    )
    return payloads


def source_audit() -> dict[str, dict[str, str]]:
    return {
        key: {
            "path": str(path.relative_to(REPO_ROOT)).replace("\\", "/"),
            "sha256": file_hash(path),
        }
        for key, path in SOURCE_PATHS.items()
    }


def jensen_polynomial(values: list[sp.Expr], degree: int, shift: int, z: sp.Symbol) -> sp.Expr:
    return sp.expand(
        sum(
            sp.binomial(degree, j) * values[shift + j] * z**j
            for j in range(degree + 1)
        )
    )


def derivative_and_transfer_certificate() -> dict[str, str | int]:
    z = sp.symbols("z")
    values = [sp.Rational((k + 2) * (3 * k + 1), 7 * (k + 1)) for k in range(20)]
    derivative_checks = 0
    for degree in range(1, 8):
        for shift in range(3):
            polynomial = jensen_polynomial(values, degree, shift, z)
            for order in range(degree + 1):
                direct = sp.diff(polynomial, z, order)
                expected = sp.factorial(degree) / sp.factorial(degree - order) * sum(
                    sp.binomial(degree - order, j)
                    * z**j
                    * values[shift + order + j]
                    for j in range(degree - order + 1)
                )
                require_zero(direct - expected, "Jensen derivative identity")
                derivative_checks += 1

    adjacent_checks = 0
    for degree in range(1, 7):
        for shift in range(4):
            upper = jensen_polynomial(values, degree + 1, shift, z)
            lower = jensen_polynomial(values, degree, shift, z)
            shifted = jensen_polynomial(values, degree, shift + 1, z)
            require_zero(upper - lower - z * shifted, "adjacent shift identity")
            adjacent_checks += 1

    transfer_checks = 0
    for root in (sp.Rational(-2), sp.Rational(-3, 2)):
        for multiplicity in range(1, 6):
            lower_unit = -root + (multiplicity + 2) * z + 2 * z**2
            upper_unit = 1 + (multiplicity + 1) * z + z**2
            lower = (z - root) ** multiplicity * lower_unit
            upper = (z - root) ** (multiplicity + 1) * upper_unit
            quotient, remainder = sp.div(sp.expand(upper - lower), z)
            require_zero(remainder, "transfer divisibility")
            transfer_checks += 1
            for order in range(multiplicity):
                require_zero(
                    sp.diff(quotient, z, order).subs(z, root),
                    "transferred root jet",
                )
                transfer_checks += 1

    return {
        "derivative": "For J_A^(d,n)(z)=sum_(j=0)^d binom(d,j)A_(n+j)z^j, (J_A^(d,n))^(s)(r)=d!/(d-s)! sum_(j=0)^(d-s)binom(d-s,j)r^jA_(n+s+j).",
        "adjacent": "Pascal's identity gives J_A^(d+1,n)=J_A^(d,n)+zJ_A^(d,n+1).",
        "transfer": "Let r!=0. If r has multiplicity at least m in J_A^(d,n) and at least m+1 in J_A^(d+1,n), then J_A^(d,n+1)=(J_A^(d+1,n)-J_A^(d,n))/z has multiplicity at least m at the same r.",
        "polar_dependency": "The existing polar-contact lemma supplies the m-to-(m+1) lift only when the adjacent degree-(d+1) window is negative-root hyperbolic. The same-root shift transfer is therefore conditional on that higher-degree hypothesis; the recurrence theorem itself merely assumes the common root.",
        "zero_guard": "The condition r!=0 is essential. Division by z can lower multiplicity at r=0, so the zero-root channel remains separate.",
        "derivative_checks": derivative_checks,
        "adjacent_checks": adjacent_checks,
        "transfer_checks": transfer_checks,
    }


def fixed_root_certificate() -> dict[str, str | int]:
    x = sp.symbols("x", integer=True, nonnegative=True)
    classification_checks = 0
    recurrence_checks = 0
    for degree in range(2, 8):
        for multiplicity in range(1, degree + 1):
            p = degree - multiplicity + 1
            q = sp.Rational(degree + 2, degree + 3)
            root = -1 / q
            polynomial = sum((h + 1) * x**h for h in range(p))

            def sequence(index: int) -> sp.Expr:
                return sp.expand(q**index * polynomial.subs(x, index))

            for shift in range(3):
                jensen = sum(
                    sp.binomial(degree, j) * sequence(shift + j) * sp.Symbol("z") ** j
                    for j in range(degree + 1)
                )
                z_symbol = sp.Symbol("z")
                for order in range(multiplicity):
                    require_zero(
                        sp.diff(jensen, z_symbol, order).subs(z_symbol, root),
                        "classification converse",
                    )
                    classification_checks += 1
            for shift in range(4):
                recurrence = sum(
                    sp.binomial(p, j)
                    * root**j
                    * sequence(shift + multiplicity - 1 + j)
                    for j in range(p + 1)
                )
                require_zero(recurrence, "highest root recurrence")
                recurrence_checks += 1

    return {
        "recurrence": "If one r!=0 is a root of multiplicity at least m in every shifted degree-d Jensen window, then with M=m-1 and p=d-m+1, (1+rE)^p A_(n+M)=0 for every n>=0.",
        "tail": "With q=-1/r, every recurrence solution has A_(M+k)=q^kP(k), deg(P)<=p-1=d-m. The lower derivative conditions determine the finite prefix consistently.",
        "differential": "For F(z)=sum_(n>=0)A_nz^n/n!, the all-shift root equations are (1+rD)^(d-s)F^(s)=0 for 0<=s<m. In particular F(z)=exp(qz)Q(z), q=-1/r, with deg(Q)<=d-m.",
        "converse": "Conversely, F=exp(qz)Q with deg(Q)<=d-m gives A_n=q^nP(n) with the same degree bound, and every J_A^(d,n) has r=-1/q as a root of multiplicity at least m.",
        "normalization": "The rational OGF sum A_nz^n has a finite pole for a nonzero geometric-polynomial sequence, but the entire zeta generating function is the EGF F=sum A_nz^n/n!. The correct zeta obstruction is the exponential-polynomial classification, not OGF entireness.",
        "indexing": "For d=4,m=2, the recurrence derived without relabelling is A_(n+1)+3rA_(n+2)+3r^2A_(n+3)+r^3A_(n+4)=0. Writing A_n through A_(n+3) requires an explicit shift of the recurrence index.",
        "classification_checks": classification_checks,
        "recurrence_checks": recurrence_checks,
    }


def persistent_certificate() -> dict[str, str | int]:
    n, z = sp.symbols("n z", integer=True, nonnegative=True)
    q = sp.Rational(1, 10)
    root = sp.Integer(-10)

    def value(index: sp.Expr) -> sp.Expr:
        return (33 * index**2 + 3 * index + 4) / (
            sp.Integer(4) * sp.Integer(10) ** index
        )

    checks = 0
    recurrence = sp.simplify(
        value(n)
        + 3 * root * value(n + 1)
        + 3 * root**2 * value(n + 2)
        + root**3 * value(n + 3)
    )
    require_zero(recurrence, "persistent recurrence")
    checks += 1

    for shift in range(11):
        polynomial = sp.expand(
            sum(sp.binomial(4, j) * value(shift + j) * z**j for j in range(5))
        )
        require_zero(polynomial.subs(z, root), "persistent common root")
        require_zero(sp.diff(polynomial, z).subs(z, root), "persistent double root")
        require(sp.diff(polynomial, z, 2).subs(z, root) != 0, "persistent exact multiplicity")
        checks += 3

    ogf = (1 + sp.Rational(7, 10) * z + sp.Rational(17, 200) * z**2) / (1 - z / 10) ** 3
    egf = (1 + sp.Rational(9, 10) * z + sp.Rational(33, 400) * z**2) * sp.exp(z / 10)
    for index in range(11):
        require_zero(sp.expand(ogf.series(z, 0, 12).removeO()).coeff(z, index) - value(index), "persistent OGF")
        require_zero(sp.diff(egf, z, index).subs(z, 0) - value(index), "persistent EGF")
        checks += 2

    numerator = 1 + sp.Rational(7, 10) * z + sp.Rational(17, 200) * z**2
    require(numerator.subs(z, 10) == sp.Rational(33, 2), "persistent pole numerator")
    checks += 1

    base_polynomial = (33 * n**2 + 3 * n + 4) / 4
    newton = [
        base_polynomial.subs(n, 0),
        (base_polynomial.subs(n, 1) - base_polynomial.subs(n, 0)),
        (
            base_polynomial.subs(n, 2)
            - 2 * base_polynomial.subs(n, 1)
            + base_polynomial.subs(n, 0)
        ),
    ]
    require(newton == [1, 9, sp.Rational(33, 2)], "persistent Newton data")
    checks += 3

    for shift in range(6):
        matrix = sp.Matrix(4, 4, lambda i, j: value(shift + i + j))
        require_zero(matrix.det(), "persistent Hankel rank")
        checks += 1

    q_symbol, a, b, c = sp.symbols("q a b c", nonzero=True)
    q_prime = 4 * q_symbol**2
    a_prime = 10 * q_symbol * a
    b_prime = 8 * q_symbol * a + 6 * q_symbol * b
    c_prime = 2 * q_symbol * (a + b + c)
    polynomial = a * n**2 + b * n + c
    parameter_flow = (
        n * (q_prime / q_symbol) * polynomial
        + a_prime * n**2
        + b_prime * n
        + c_prime
    )
    hierarchy = 2 * (2 * n + 1) * q_symbol * (
        a * (n + 1) ** 2 + b * (n + 1) + c
    )
    require_zero(parameter_flow - hierarchy, "geometric-quadratic hierarchy closure")
    checks += 1

    return {
        "sequence": "A_n=(33n^2+3n+4)/(4*10^n) has the exact common quartic double root r=-10 and satisfies the corrected cubic recurrence at every shift.",
        "ogf": "sum A_nz^n=(1+7z/10+17z^2/200)/(1-z/10)^3, whose numerator equals 33/2 at z=10.",
        "egf": "sum A_nz^n/n!=(1+9z/10+33z^2/400)exp(z/10), which is entire and is exactly the classified exponential-polynomial form.",
        "hankel": "Every contiguous order-four Hankel determinant of this rank-three sequence vanishes.",
        "hierarchy": "The geometric-quadratic class is invariant under A_n'=2(2n+1)A_(n+1), with q'=4q^2, a'=10qa, b'=8qa+6qb, and c'=2q(a+b+c).",
        "scope": "This exact witness validates the recurrence and the normalization correction. The separately reported canonical continuation and its A_10 value are not promoted here because its native source data were not supplied.",
        "checks": checks,
    }


def hankel_certificate() -> dict[str, str | int]:
    hankel_checks = 0
    r = sp.symbols("r", nonzero=True)
    for order in range(1, 6):
        symbols = sp.symbols(f"a0:{2 * order + 1}")
        matrix = sp.Matrix(order + 1, order + 1, lambda i, j: symbols[i + j])
        coefficients = sp.Matrix([sp.binomial(order, j) * r**j for j in range(order + 1)])
        product = matrix * coefficients
        for row in range(order + 1):
            expected = sum(
                sp.binomial(order, j) * r**j * symbols[row + j]
                for j in range(order + 1)
            )
            require_zero(product[row] - expected, "Hankel recurrence vector")
            hankel_checks += 1

        q = sp.Rational(order + 2, order + 3)
        degree = order - 1
        polynomial = lambda k: sum((h + 2) * k**h for h in range(degree + 1))
        sequence = lambda k: q**k * polynomial(k)
        sample = sp.Matrix(order + 1, order + 1, lambda i, j: sequence(2 + i + j))
        require_zero(sample.det(), "geometric-polynomial Hankel rank")
        hankel_checks += 1

    return {
        "matrix": "Let c_p(r)=(binom(p,0),binom(p,1)r,...,r^p)^T and H_(p+1,M+n)=[A_(M+n+i+j)]_(i,j=0)^p. The vector of p+1 consecutive recurrence defects is H_(p+1,M+n)c_p(r).",
        "rank": "If p+1 consecutive fixed-root recurrences are exact, H_(p+1,M+n)c_p(r)=0. Since the first component of c_p(r) is 1, det H_(p+1,M+n)=0.",
        "first_boundary": "Let p=d-m+1. If J_A^(d,n) has one nonzero multiplicity-m root and the p consecutive extensions J_A^(d+1,n),...,J_A^(d+1,n+p-1) are negative-root hyperbolic, repeated polar lifting and adjacent transfer give the same multiplicity-m root in shifts n through n+p. Therefore det H_(p+1,m-1+n)=0. At a first loss of global Jensen hyperbolicity, every fixed higher-degree window remains hyperbolic by closure, so a nonzero determinant excludes that degree-d collision.",
        "quartic": "For d=4,m=2, p=3 and M=1. Four consecutive fixed-root double tangencies force H_(4,n+1)=det[A_(n+1+i+j)]_(i,j=0)^3=0.",
        "m100": "The certified anchor theorem H_(4,k)(-100)>0 for every k>=0 excludes four consecutive quartic double tangencies at one fixed nonzero root at lambda=-100.",
        "interval": "The certified forward-invariance theorem strengthens the anchor to H_(4,k)(lambda)>0 for every k>=0 and every lambda in [-100,0]. Hence the same fixed-root four-chain is impossible throughout that heat interval.",
        "conditional": "At any lambda in [-100,0], if a quartic double root at shift n has negative-root hyperbolic quintic extensions at shifts n,n+1,n+2, polar lifting and adjacent transfer create the four forbidden quartic tangencies. The initial tangency and all three higher-degree hypotheses therefore cannot coexist anywhere on the interval.",
        "checks": hankel_checks,
        "m100_applications": 1,
        "heat_interval_applications": 1,
    }


def finite_difference(values: list[sp.Expr], order: int, index: int) -> sp.Expr:
    return sp.simplify(
        sum(
            (-1) ** (order - j) * sp.binomial(order, j) * values[index + j]
            for j in range(order + 1)
        )
    )


def near_rigidity_certificate() -> dict[str, str | int]:
    newton_checks = 0
    bound_checks = 0
    values = [sp.Rational(3 * k**4 - 2 * k**3 + 5 * k + 7, k + 2) for k in range(14)]
    for order in range(1, 6):
        for k in range(order, 14):
            polynomial = sum(
                sp.binomial(k, s) * finite_difference(values, s, 0)
                for s in range(order)
            )
            remainder = sum(
                sp.binomial(k - 1 - j, order - 1)
                * finite_difference(values, order, j)
                for j in range(k - order + 1)
            )
            require_zero(values[k] - polynomial - remainder, "Newton remainder")
            newton_checks += 1
            maximum = max(
                abs(finite_difference(values, order, j))
                for j in range(k - order + 1)
            )
            require(abs(remainder) <= sp.binomial(k, order) * maximum, "Newton bound")
            bound_checks += 1

    r0, r1 = sp.symbols("r0 r1")
    drift_checks = 0
    for order in range(1, 7):
        for power in range(1, order + 1):
            quotient = sum(r0 ** (power - 1 - ell) * r1**ell for ell in range(power))
            require_zero(r0**power - r1**power - (r0 - r1) * quotient, "root drift factor")
            drift_checks += 1

    return {
        "normalized_defect": "For epsilon_n=(1+rE)^pA_(M+n), q=-1/r, and B_n=q^(-n)A_(M+n), one has Delta^pB_n=(-1)^p q^(-n)epsilon_n exactly.",
        "newton": "Let P_(p-1)(k)=sum_(s=0)^(p-1)binom(k,s)Delta^sB_0. Then B_k-P_(p-1)(k)=sum_(j=0)^(k-p)binom(k-1-j,p-1)Delta^pB_j.",
        "bound": "Consequently |B_k-P_(p-1)(k)|<=binom(k,p)max_(0<=j<=k-p)|q|^(-j)|epsilon_j|. For quartic double roots the amplification factor is binom(k,3).",
        "hankel_margin": "For the p+1 defect vector epsilon=Hc, ||epsilon||_2>=sigma_min(H)||c||_2>=|det H| ||c||_2/||H||_2^p. Hence max_i|epsilon_i| is at least this quantity divided by sqrt(p+1).",
        "drift": "If row i is exact at r_i but tested at r_*, its residual is sum_(j=1)^pbinom(p,j)(r_*^j-r_i^j)A_(M+n+i+j). On |r_i|,|r_*|<=R, use |r_*^j-r_i^j|<=jR^(j-1)|r_*-r_i| and compare with the Hankel-margin lower bound.",
        "open": "A zeta-specific near-chain theorem still needs certified singular-value or determinant/norm margins on the relevant heat interval, a compact annulus excluding r=0 and infinity, and control of the adjacent-root drift. Finite closeness to an exponential-polynomial does not contradict entireness by itself.",
        "newton_checks": newton_checks,
        "bound_checks": bound_checks,
        "drift_checks": drift_checks,
    }


def handoff_certificate() -> dict[str, str]:
    return {
        "decision": "Promote the fixed-shift theorem in EGF normalization and use finite Hankel margins before invoking an infinite-chain pole argument.",
        "two_axis": "The existing polar cascade controls a fixed shift while degree grows. This gate controls a fixed degree while shift grows. For a non-exponential-polynomial zeta source, an exact obstruction cannot remain in either fixed-root axis indefinitely.",
        "next": "Build interval lower bounds for sigma_min(H_(4,n)) or a determinant-plus-norm substitute on the heat interval where tangency propagation is needed. Combine those margins with the root-drift residual and the adjacent quintic transfer hypothesis.",
        "separation": "This is a Jensen/PF branch. It does not alter the distinct Section 11.193 carrier-plus-near-kernel obligation or license mixing their estimates.",
        "reserve": "No approximate-chain exclusion on the zeta heat interval, all-degree or all-shift hyperbolicity theorem, PF-infinity, Lambda<=0, RH, or prize-level conclusion is proved.",
    }


def build_rows(certificate: dict) -> list[GateRow]:
    transfer = certificate["transfer"]
    fixed = certificate["fixed_root"]
    persistent = certificate["persistent"]
    hankel = certificate["hankel"]
    near = certificate["near_rigidity"]
    handoff = certificate["handoff"]
    rows = [
        GateRow("frt_01_sources", "source chain", "proved", "The six parent artifacts are current and hash-pinned.", "Parent row counts, the order-four anchor, and the full heat-interval theorem are checked.", "No parent result is strengthened."),
        GateRow("frt_02_derivative", "Jensen derivative", "proved", "Every derivative has the exact shifted binomial form.", transfer["derivative"], "Finite polynomial identity."),
        GateRow("frt_03_adjacent", "adjacent shift", "proved", "Adjacent degrees expose the next equal-degree shift.", transfer["adjacent"], "No hyperbolicity follows from this identity."),
        GateRow("frt_04_transfer", "same-root transfer", "proved", "A lifted multiplicity transfers the same nonzero root to the next shift.", transfer["transfer"], transfer["zero_guard"]),
        GateRow("frt_05_polar_dependency", "higher-degree dependency", "guard_validated", "The multiplicity lift requires the parent hyperbolicity hypothesis.", transfer["polar_dependency"], "A repeated tangency chain is not automatic without this input."),
        GateRow("frt_06_zero_root", "zero-root guard", "guard_validated", "The transfer theorem excludes r=0.", transfer["zero_guard"], "The zero-root channel remains open."),
        GateRow("frt_07_recurrence", "common-root recurrence", "proved", "A fixed common root forces a repeated-root recurrence.", fixed["recurrence"], "The theorem assumes every indicated shift."),
        GateRow("frt_08_indexing", "quartic indexing", "proved", "The quartic recurrence starts at A_(n+1) before relabelling.", fixed["indexing"], "No silent shift is allowed."),
        GateRow("frt_09_tail", "recurrence solution", "proved", "The normalized coefficient tail is polynomial.", fixed["tail"], "This is exact fixed-root rigidity."),
        GateRow("frt_10_differential", "EGF equation", "proved", "The all-shift equations become constant-coefficient ODEs.", fixed["differential"], "The EGF normalization is essential."),
        GateRow("frt_11_classification", "exponential-polynomial", "proved", "The entire Jensen generating function is exponential-polynomial.", fixed["differential"], "This class is entire, so entireness alone is not a contradiction."),
        GateRow("frt_12_converse", "classification converse", "proved", "Every classified source has the common-root chain.", fixed["converse"], "Multiplicity may exceed the stated lower bound in degenerate cases."),
        GateRow("frt_13_ogf", "ordinary series", "proved", "The ordinary series has a finite pole for a nonzero geometric-polynomial sequence.", persistent["ogf"], fixed["normalization"]),
        GateRow("frt_14_normalization", "zeta normalization", "guard_validated", "The OGF pole is not the zeta entireness contradiction.", fixed["normalization"], "Use the EGF exponential-polynomial classification."),
        GateRow("frt_15_persistent", "persistent recurrence", "proved", "The positive persistent source realizes the quartic recurrence exactly.", persistent["sequence"], persistent["scope"]),
        GateRow("frt_16_persistent_root", "persistent common root", "proved", "Every tested quartic has exact double root -10.", persistent["sequence"], "The symbolic recurrence proves persistence beyond the tests."),
        GateRow("frt_17_persistent_ogf", "persistent OGF", "proved", "The reported pole numerator is exact.", persistent["ogf"], "This is the ordinary series only."),
        GateRow("frt_18_persistent_egf", "persistent EGF", "proved", "The same source has an entire exponential-polynomial EGF.", persistent["egf"], "This is the normalization correction witness."),
        GateRow("frt_19_hierarchy", "hierarchy closure", "proved", "The geometric-quadratic family closes under the coefficient flow.", persistent["hierarchy"], "Closure does not imply membership of the zeta source."),
        GateRow("frt_20_source_scope", "source boundary", "guard_validated", "Only reproducible supplied data are promoted.", persistent["scope"], "The canonical A_10 claim awaits its native source."),
        GateRow("frt_21_hankel_vector", "Hankel bridge", "proved", "Consecutive recurrence defects form one Hankel matrix-vector product.", hankel["matrix"], "No determinant sign is inferred generically."),
        GateRow("frt_22_rank", "finite rank defect", "proved", "p+1 exact recurrences force a finite Hankel determinant to vanish.", hankel["rank"], "This is stronger than requiring an infinite chain."),
        GateRow("frt_23_quartic_hankel", "quartic specialization", "proved", "Four fixed-root quartic double tangencies force H_4=0.", hankel["quartic"], "The root must be common and nonzero."),
        GateRow("frt_24_heat_interval", "heat-interval application", "proved", "All-shift H_4 positivity excludes such four-chains throughout [-100,0].", hankel["interval"], "The root must remain common across the four rows."),
        GateRow("frt_25_conditional_transfer", "finite first-boundary reduction", "proved", "p hyperbolic upper extensions force one finite Hankel rank defect.", hankel["first_boundary"], "This does not supply all-order Hankel nonvanishing or the higher-degree hypotheses away from a first boundary."),
        GateRow("frt_26_normalized_defect", "near recurrence", "proved", "Recurrence defect is a normalized p-th finite difference.", near["normalized_defect"], "Fixed reference root only."),
        GateRow("frt_27_newton", "Newton reconstruction", "proved", "The near-chain remainder has an exact discrete Green kernel.", near["newton"], "No asymptotic notation is used."),
        GateRow("frt_28_near_bound", "near-rigidity bound", "proved", "Small fixed-root defects force a finite block near a geometric-polynomial tail.", near["bound"], "The binomial amplification must be retained."),
        GateRow("frt_29_hankel_margin", "quantitative obstruction", "proved", "A nonzero singular-value margin forces at least one tangency defect.", near["hankel_margin"], "Numerical interval margins are not supplied here."),
        GateRow("frt_30_root_drift", "drifting root", "open", "Control root drift by comparing row coefficients to one reference root.", near["drift"], near["open"]),
        GateRow("frt_31_zeta_near", "zeta near-chain", "open", "Supply heat-uniform Hankel margins and a compact root annulus.", handoff["next"], handoff["reserve"]),
        GateRow("frt_32_boundary", "proof boundary", "guard_validated", "This gate is not a proof of RH.", handoff["reserve"], handoff["separation"]),
    ]
    require(len(rows) == 32, "row count drifted")
    return rows


def render_note(artifact: dict) -> str:
    cert = artifact["symbolic_certificate"]
    transfer = cert["transfer"]
    fixed = cert["fixed_root"]
    persistent = cert["persistent"]
    hankel = cert["hankel"]
    near = cert["near_rigidity"]
    handoff = cert["handoff"]
    return f"""# Fixed-Root Shift-Tangency Rigidity Gate

Date: 2026-08-03

Status: exact fixed-root recurrence, EGF classification, finite Hankel obstruction, and fixed-root near-rigidity theorem; drifting-root zeta estimate open; not a proof of RH.

## Adjacent Transfer

{transfer['derivative']}

{transfer['adjacent']}

{transfer['transfer']}

{transfer['polar_dependency']}

{transfer['zero_guard']}

## Fixed-Root Classification

{fixed['recurrence']}

{fixed['indexing']}

{fixed['tail']}

{fixed['differential']}

{fixed['converse']}

{fixed['normalization']}

## Persistent Rational Witness

{persistent['sequence']}

{persistent['ogf']}

{persistent['egf']}

{persistent['hankel']}

{persistent['hierarchy']}

{persistent['scope']}

## Finite Hankel Bridge

{hankel['matrix']}

{hankel['rank']}

{hankel['first_boundary']}

{hankel['quartic']}

{hankel['m100']}

{hankel['interval']}

{hankel['conditional']}

## Quantitative Near-Rigidity

{near['normalized_defect']}

{near['newton']}

{near['bound']}

{near['hankel_margin']}

{near['drift']}

{near['open']}

## Route Decision

{handoff['decision']}

{handoff['two_axis']}

{handoff['next']}

{handoff['separation']}

## Pi Provenance

This gate introduces no pi. Every identity is binomial, finite-difference, rational-function, polynomial, differential-operator, or Hankel algebra.

## Proof Boundary

{handoff['reserve']} This gate is not a proof of RH.
"""


def main() -> int:
    load_sources()
    certificate = {
        "transfer": derivative_and_transfer_certificate(),
        "fixed_root": fixed_root_certificate(),
        "persistent": persistent_certificate(),
        "hankel": hankel_certificate(),
        "near_rigidity": near_rigidity_certificate(),
        "handoff": handoff_certificate(),
    }
    rows = build_rows(certificate)
    counts = {
        "rows": len(rows),
        "derivative_checks": certificate["transfer"]["derivative_checks"],
        "adjacent_shift_checks": certificate["transfer"]["adjacent_checks"],
        "same_root_transfer_checks": certificate["transfer"]["transfer_checks"],
        "classification_checks": certificate["fixed_root"]["classification_checks"],
        "recurrence_checks": certificate["fixed_root"]["recurrence_checks"],
        "persistent_witness_checks": certificate["persistent"]["checks"],
        "hankel_checks": certificate["hankel"]["checks"],
        "m100_applications": certificate["hankel"]["m100_applications"],
        "heat_interval_applications": certificate["hankel"]["heat_interval_applications"],
        "newton_identity_checks": certificate["near_rigidity"]["newton_checks"],
        "newton_bound_checks": certificate["near_rigidity"]["bound_checks"],
        "root_drift_factor_checks": certificate["near_rigidity"]["drift_checks"],
        "zeta_near_chain_bounds": 0,
        "rh_conclusions": 0,
    }
    artifact = {
        "kind": KIND,
        "date": "2026-08-03",
        "status": "fixed-root shift recurrence, EGF classification, finite Hankel heat-interval obstruction, and fixed-root near-rigidity proved; drifting-root zeta estimate open",
        "counts": counts,
        "symbolic_certificate": certificate,
        "source_audit": source_audit(),
        "rows": [asdict(row) for row in rows],
        "proof_boundary": "This proves the exact adjacent-shift identity and conditional same-root transfer, fixed-root recurrence, EGF exponential-polynomial classification, persistent witness normalization, finite Hankel rank obstruction, fixed-root four-chain exclusion throughout lambda in [-100,0], and fixed-root near-rigidity formula. It proves no heat-uniform drifting-root exclusion, all-degree or all-shift hyperbolicity theorem, PF-infinity, Lambda<=0, RH, or prize-level conclusion.",
    }
    atomic_write(RESULT_PATH, json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    atomic_write(NOTE_PATH, render_note(artifact))
    print(
        "built fixed-root shift-tangency rigidity gate: "
        f"{counts['rows']} rows, {counts['derivative_checks']} derivative checks, "
        f"{counts['adjacent_shift_checks']} adjacent shifts, "
        f"{counts['same_root_transfer_checks']} transfer checks, "
        f"{counts['classification_checks']} classification checks, "
        f"{counts['hankel_checks']} Hankel checks, "
        f"{counts['newton_identity_checks']} Newton identities, "
        f"{counts['m100_applications']} lambda=-100 anchor, "
        f"{counts['heat_interval_applications']} heat-interval application, "
        f"{counts['zeta_near_chain_bounds']} zeta near-chain bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
