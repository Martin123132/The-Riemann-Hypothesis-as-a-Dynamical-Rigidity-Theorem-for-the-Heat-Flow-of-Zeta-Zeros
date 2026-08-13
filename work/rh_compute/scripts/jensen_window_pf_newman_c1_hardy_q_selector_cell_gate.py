#!/usr/bin/env python3
"""Build selector-stable local q-input cells around the saved cubic roster."""

from __future__ import annotations

from decimal import Decimal, ROUND_CEILING, localcontext
from fractions import Fraction
import hashlib
import json
import math
import os
from pathlib import Path
import sys
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
FLINT_ROOT = Path(
    os.environ.get(
        "RH_PYTHON_FLINT_ROOT",
        r"C:\Users\ollet\Documents\Codex\third_party\python_flint_0_8_0",
    )
)
sys.path.insert(0, str(FLINT_ROOT))

try:
    import flint
    from flint import arb, ctx
except ImportError as exc:  # pragma: no cover
    raise RuntimeError(f"python-flint runtime unavailable at {FLINT_ROOT}") from exc


SOURCE = (
    REPO_ROOT
    / "work/rh_compute/external/hardy_fastcode/resumable/zeta14cubicmult_resumable.f90"
)
TELEMETRY = (
    REPO_ROOT
    / "work/rh_compute/external/hardy_fastcode/fixtures/chain_telemetry/t1e10/enabled/chain_telemetry.jsonl"
)
CHILD_ADAPTER = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_cubic_child_legendre_adapter_gate.json"
)
RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_q_selector_cell_gate.json"
)
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_q_selector_cell_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = REPO_ROOT / "work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_q_selector_cell_gate.py"
SOURCE_SHA256 = "0fb64f090c21b185f27edc1c9194254cf471e74bfd6af3b88cb7f0cb40ec0c3d"
EXPECTED_TPM_HEX = "C001921FB54442D18469898CC51701B8"
TOL = Fraction(1, 10**12)
ERF_BOUNDARIES = tuple(
    Fraction(value)
    for value in (
        0,
        Fraction(1, 10**14),
        Fraction(4, 5),
        Fraction(9, 8),
        Fraction(11, 8),
        Fraction(13, 8),
        Fraction(15, 8),
        Fraction(17, 8),
        Fraction(19, 8),
        Fraction(21, 8),
    )
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def fraction_record(value: Fraction) -> dict[str, int]:
    return {"numerator": value.numerator, "denominator": value.denominator}


def decimal_text(value: Fraction, digits: int = 70) -> str:
    with localcontext() as context:
        context.prec = digits
        return format(Decimal(value.numerator) / Decimal(value.denominator), ".50E")


def upper_decimal(value: Fraction, digits: int = 90) -> str:
    with localcontext() as context:
        context.prec = digits
        context.rounding = ROUND_CEILING
        return format(Decimal(value.numerator) / Decimal(value.denominator), ".80E")


def binary128_fraction(payload: str) -> Fraction:
    require(isinstance(payload, str) and len(payload) == 32, "invalid binary128 payload")
    bits = int(payload, 16)
    sign = -1 if bits >> 127 else 1
    exponent = (bits >> 112) & 0x7FFF
    fraction = bits & ((1 << 112) - 1)
    require(exponent != 0x7FFF, "non-finite binary128 payload")
    if exponent == 0:
        if fraction == 0:
            return Fraction(0)
        mantissa = fraction
        binary_exponent = 1 - 16383 - 112
    else:
        mantissa = (1 << 112) + fraction
        binary_exponent = exponent - 16383 - 112
    if binary_exponent >= 0:
        return Fraction(sign * mantissa * 2**binary_exponent, 1)
    return Fraction(sign * mantissa, 2 ** (-binary_exponent))


def fraction_ball(value: Fraction) -> arb:
    return arb(value.numerator) / arb(value.denominator)


def interval_ball(center: Fraction, radius: Fraction) -> arb:
    require(radius > 0, "nonpositive interval radius")
    return fraction_ball(center) + arb(f"[0 +/- {upper_decimal(radius)}]")


def bound_fraction(value: arb) -> Fraction:
    mantissa, exponent = value.man_exp()
    mantissa = int(mantissa)
    exponent = int(exponent)
    return Fraction(mantissa * 2**exponent, 1) if exponent >= 0 else Fraction(mantissa, 2 ** (-exponent))


def bounds(value: arb) -> tuple[Fraction, Fraction]:
    return bound_fraction(value.lower()), bound_fraction(value.upper())


def nearest_integer(value: Fraction) -> int:
    if value >= 0:
        return math.floor(value + Fraction(1, 2))
    return math.ceil(value - Fraction(1, 2))


def source_locations() -> dict[str, int]:
    lines = SOURCE.read_text(encoding="utf-8").splitlines()
    tokens = {
        "q_sqrt_xr": "sx=sqrt(xr(nit))",
        "q_endpoint": "t3=0.5*(1.0+endpoint)",
        "q_frac_lower": "con1=1.0-fracL(nit)",
        "q_psi_lower": "call psi(con2,ps2)",
        "q_psi_phase": "call psi(phicoeff(1,nit)+1.0,ps2)",
        "q_erf_first": "call erf(z,6,cr1)",
        "q_linear_branch": "if (phicoeff(1,nit).gt.0.0) then",
        "q_jbot": "jbot=ceiling(phicoeff(1,nit))",
        "q_ip_cap": "if (ip1.gt.(0.5*L(nit))) then",
        "q_lower_denominator": "con1=i*1.0-phicoeff(1,nit)",
        "q_upper_denominator": "c1=L(nit)+fracL(nit)",
        "q_newton_branch": "if (m1(nit).gt.5) then",
        "psi_zero_guard": "IF (XX.LT.1e-12) THEN",
        "psi_negative_integer_guard": "IF (X.LT.0.0.AND.ABS(X-NINT(X)).LT.1e-12) THEN",
        "psi_small_large_split": "IF (XX.LT.4.0) THEN",
        "erf_first_taylor_split": "IF (ZA.GE.0.8.AND.ZA.LT.1.125) THEN",
        "erf_asymptotic_split": "ELSE IF (ZA.GE.2.625) THEN",
    }
    locations: dict[str, int] = {}
    for name, token in tokens.items():
        matches = [index for index, line in enumerate(lines, 1) if token in line]
        require(matches, f"source q-selector token missing: {token}")
        locations[name] = matches[0]
    return locations


def load_chains() -> dict[int, dict[str, Any]]:
    chains: dict[int, dict[str, Any]] = {}
    for line in TELEMETRY.read_text(encoding="utf-8").splitlines():
        record = json.loads(line)
        chain_id = int(record.get("chain", 0))
        if record["type"] == "chain":
            chains[chain_id] = {"header": record, "levels": {}}
        elif record["type"] == "level":
            chains[chain_id]["levels"][int(record["level"])] = record
    require(sorted(chains) == list(range(1, 65)), "q-cell chain sequence drift")
    return chains


def psi_boundaries() -> tuple[Fraction, ...]:
    values = {Fraction(0)}
    for half_index in range(-12, 13):
        values.add(Fraction(half_index, 2))
    for special in (Fraction(1, 2), Fraction(1), Fraction(2), Fraction(3)):
        for sign in (-1, 1):
            values.add(sign * (special - TOL))
            values.add(sign * (special + TOL))
    for integer in range(-6, 1):
        values.add(Fraction(integer) - TOL)
        values.add(Fraction(integer) + TOL)
    return tuple(sorted(values))


PSI_BOUNDARIES = psi_boundaries()


def boundary_margin(value: arb, boundaries: tuple[Fraction, ...]) -> Fraction | None:
    lower, upper = bounds(value)
    margins: list[Fraction] = []
    for boundary in boundaries:
        if lower <= boundary <= upper:
            return None
        margins.append(lower - boundary if boundary < lower else boundary - upper)
    return min(margins)


def interval_contains(value: arb, lower_target: Fraction, upper_target: Fraction) -> bool:
    lower, upper = bounds(value)
    return lower_target < lower <= upper < upper_target


def selector_cell(
    chain_id: int,
    header: dict[str, Any],
    level_one: dict[str, Any],
    level_two: dict[str, Any],
) -> dict[str, Any]:
    a1, a2, a3 = [binary128_fraction(value) for value in level_one["coefficients_hex"]]
    y = binary128_fraction(level_one["xr_hex"])
    frac = binary128_fraction(level_one["frac_length_hex"])
    require(y == 2 * a2 and y > 0, f"chain {chain_id} parent xr drift")
    require(0 < frac < 1, f"chain {chain_id} fracL drift")
    require(a1 != 0 and abs(a1) < Fraction(1, 2), f"chain {chain_id} a1 branch drift")
    child_length = int(level_two["length"])
    require(child_length in (1, 2), f"chain {chain_id} child length drift")

    center_y = y
    s1 = Fraction(1, 2) / center_y
    s2 = -a3 / center_y**3
    x = -a1
    raw_r1 = 2 * s1 * x + 3 * s2 * x**2
    raw_r2 = s1 + 3 * s2 * x
    n1 = nearest_integer(raw_r1)
    ix = nearest_integer(2 * raw_r2)
    raw_linear_sign = 1 if raw_r1 - n1 > 0 else -1
    child_xr_sign = 1 if 2 * raw_r2 - ix > 0 else -1
    jbot = math.ceil(a1)

    base_radii = {
        "a1": min(abs(a1), Fraction(1, 2) - abs(a1)) / 16,
        "y": y / 16,
        "a3": abs(a3) / 4,
        "fracL": min(frac, 1 - frac) / 16,
    }
    require(all(radius > 0 for radius in base_radii.values()), "degenerate base q cell")

    ctx.dps = 180
    ctx.threads = 1
    p_hat = -binary128_fraction(header["tpm_hex"]) / 2
    sp = fraction_ball(p_hat).sqrt() + arb(f"[0 +/- {upper_decimal(Fraction(1, 2**110))}]")

    accepted: dict[str, Any] | None = None
    for halvings in range(0, 61):
        factor = Fraction(1, 2**halvings)
        radii = {name: radius * factor for name, radius in base_radii.items()}
        a1_ball = interval_ball(a1, radii["a1"])
        y_ball = interval_ball(y, radii["y"])
        a3_ball = interval_ball(a3, radii["a3"])
        frac_ball = interval_ball(frac, radii["fracL"])

        if not interval_contains(a1_ball, Fraction(-1, 2), Fraction(1, 2)):
            continue
        if not interval_contains(y_ball, Fraction(0), Fraction(1)):
            continue
        if not interval_contains(frac_ball, Fraction(0), Fraction(1)):
            continue
        if boundary_margin(a1_ball, (Fraction(0),)) is None:
            continue

        source_x = -a1_ball
        source_s1 = 1 / (2 * y_ball)
        source_s2 = -a3_ball / y_ball**3
        r1_ball = 2 * source_s1 * source_x + 3 * source_s2 * source_x**2
        r2_ball = source_s1 + 3 * source_s2 * source_x
        if not interval_contains(r1_ball, Fraction(n1, 1) - Fraction(1, 2), Fraction(n1, 1) + Fraction(1, 2)):
            continue
        linear_residual = r1_ball - n1
        if raw_linear_sign > 0:
            if bounds(linear_residual)[0] <= 0:
                continue
        elif bounds(linear_residual)[1] >= 0:
            continue
        if not interval_contains(2 * r2_ball, Fraction(ix, 1) - Fraction(1, 2), Fraction(ix, 1) + Fraction(1, 2)):
            continue
        child_xr = 2 * r2_ball - ix
        if child_xr_sign > 0:
            if bounds(child_xr)[0] <= 0:
                continue
        elif bounds(child_xr)[1] >= 0:
            continue

        psi_arguments = {
            "one_minus_fracL": 1 - frac_ball,
            "length_plus_fracL_plus_one": child_length + frac_ball + 1,
            "a1_plus_one": a1_ball + 1,
            "length_plus_one_minus_a1": child_length + 1 - a1_ball,
            "jbot_minus_a1": jbot - a1_ball,
            "jbot_minus_length_fracL": jbot - (child_length + frac_ball),
        }
        psi_margins = {name: boundary_margin(value, PSI_BOUNDARIES) for name, value in psi_arguments.items()}
        if any(margin is None for margin in psi_margins.values()):
            continue

        q_con1 = a1_ball if a1 > 0 else a1_ball + 1
        z_first = sp * (1 - frac_ball) / y_ball.sqrt()
        z_second = sp * q_con1 / y_ball.sqrt()
        erf_margins = {
            "frac_endpoint_argument": boundary_margin(z_first, ERF_BOUNDARIES),
            "linear_endpoint_argument": boundary_margin(z_second, ERF_BOUNDARIES),
        }
        if any(margin is None for margin in erf_margins.values()):
            continue

        denominator_lowers: list[Fraction] = []
        discriminant_lowers: list[Fraction] = []
        branch_parameter_uppers: list[Fraction] = []
        h1 = 3 * a3_ball / y_ball**2
        for index in range(child_length + 1):
            u = index - a1_ball
            denominator = 1 + h1 * u - (arb(3) / 2) * h1**2 * u**2
            d_lower, d_upper = bounds(denominator)
            if d_lower <= 0 <= d_upper:
                break
            denominator_lowers.append(min(abs(d_lower), abs(d_upper)))
            branch_parameter = 12 * a3_ball * u / y_ball**2
            discriminant = 1 + branch_parameter
            disc_lower, _ = bounds(discriminant)
            if disc_lower <= 0:
                break
            discriminant_lowers.append(disc_lower)
            z_lower, z_upper = bounds(branch_parameter)
            branch_parameter_uppers.append(max(abs(z_lower), abs(z_upper)))
        if len(denominator_lowers) != child_length + 1:
            continue

        accepted = {
            "halvings": halvings,
            "scale": fraction_record(factor),
            "radii": {
                name: {"exact": fraction_record(radius), "decimal": decimal_text(radius)}
                for name, radius in radii.items()
            },
            "child_selectors": {
                "n1": n1,
                "ix": ix,
                "linear_residual_sign": raw_linear_sign,
                "child_xr_sign": child_xr_sign,
                "source_conjugates": bool(level_two["conjugate"]),
                "minimum_linear_selector_boundary_gap": decimal_text(
                    boundary_margin(r1_ball, (Fraction(n1, 1) - Fraction(1, 2), Fraction(n1, 1) + Fraction(1, 2)))
                ),
                "minimum_quadratic_selector_boundary_gap": decimal_text(
                    boundary_margin(2 * r2_ball, (Fraction(ix, 1) - Fraction(1, 2), Fraction(ix, 1) + Fraction(1, 2)))
                ),
            },
            "q_branches": {
                "a1_sign": 1 if a1 > 0 else -1,
                "jbot": jbot,
                "ip1": 1,
                "newton_used": False,
                "psi_call_count": len(psi_arguments),
                "erf_call_count": 2,
                "minimum_psi_boundary_gap": decimal_text(min(margin for margin in psi_margins.values() if margin is not None)),
                "minimum_erf_boundary_gap": decimal_text(min(margin for margin in erf_margins.values() if margin is not None)),
            },
            "analytic_margins": {
                "minimum_denominator_abs_lower": decimal_text(min(denominator_lowers)),
                "minimum_stationary_discriminant_lower": decimal_text(min(discriminant_lowers)),
                "maximum_branch_parameter_abs_upper": decimal_text(max(branch_parameter_uppers)),
                "minimum_frac_endpoint_margin": decimal_text(min(frac - radii["fracL"], 1 - frac - radii["fracL"])),
                "minimum_parent_xr_lower": decimal_text(y - radii["y"]),
            },
        }
        break

    require(accepted is not None, f"chain {chain_id} no selector-stable cell after 60 halvings")
    return {
        "chain": chain_id,
        "block": int(header["block"]),
        "sum_index": int(header["sum_index"]),
        "branch": int(header["branch"]),
        "child_length": child_length,
        "center": {
            "a1_hex": level_one["coefficients_hex"][0],
            "xr_hex": level_one["xr_hex"],
            "a3_hex": level_one["coefficients_hex"][2],
            "fracL_hex": level_one["frac_length_hex"],
        },
        **accepted,
    }


def minimum_decimal(rows: list[dict[str, Any]], path: tuple[str, ...]) -> str:
    values = []
    for row in rows:
        value: Any = row
        for key in path:
            value = value[key]
        values.append(Decimal(value))
    return format(min(values), ".45E")


def build_artifact() -> dict[str, Any]:
    for path in (SOURCE, TELEMETRY, CHILD_ADAPTER, CHECKER):
        require(path.is_file(), f"missing q-selector dependency: {path}")
    require(file_hash(SOURCE) == SOURCE_SHA256, "accepted source hash drift")
    require(flint.__version__ == "0.8.0", "python-flint version drift")
    child = json.loads(CHILD_ADAPTER.read_text(encoding="utf-8"))
    require(child["saved_fixture"]["chain_count"] == 64, "child adapter dependency drift")
    chains = load_chains()
    rows = [
        selector_cell(
            chain_id,
            chains[chain_id]["header"],
            chains[chain_id]["levels"][1],
            chains[chain_id]["levels"][2],
        )
        for chain_id in sorted(chains)
    ]

    require(all(row["sum_index"] == (row["chain"] + 1) // 2 for row in rows), "sum-index roster drift")
    return {
        "kind": "jensen_window_pf_newman_c1_hardy_q_selector_cell_gate",
        "status": "rigorous_local_q_selector_cells_with_physical_coverage_and_analytic_values_open",
        "source": {
            "path": relative(SOURCE),
            "sha256": file_hash(SOURCE),
            "locations": source_locations(),
        },
        "cell_contract": {
            "coordinates": ["a1", "xr=2*a2", "a3", "fracL"],
            "discrete_labels": "L(0)=104, L(1) in {1,2}, m1(1)=m1(2)=3, ip=3",
            "construction": (
                "Start from one-sixteenth of the elementary sign/integer margins (one-quarter for a3), "
                "then halve all four radii until Arb proves every selector and analytic denominator stable."
            ),
            "preserved": [
                "child n1 and ix nearest-integer selectors",
                "child linear-residual sign and xr orientation/conjugation",
                "q a1-sign, jbot, ip1, and no-Newton branches",
                "all six PSI source branch paths",
                "both ERF source branch paths",
                "fracL endpoint separation and parent xr positivity",
                "cubic stationary discriminant positivity",
                "centered kernel denominator nonvanishing",
            ],
            "not_preserved_or_proved": (
                "No physical rae/t trajectory is yet enclosed in these derived-input boxes, and no source "
                "PSI/ERF approximation error or analytic q value is bounded."
            ),
        },
        "aggregate": {
            "cell_count": len(rows),
            "branch_counts": {str(branch): sum(row["branch"] == branch for row in rows) for branch in (1, 2)},
            "sum_index_range": [min(row["sum_index"] for row in rows), max(row["sum_index"] for row in rows)],
            "maximum_halvings": max(row["halvings"] for row in rows),
            "halving_histogram": {
                str(value): sum(row["halvings"] == value for row in rows)
                for value in sorted({row["halvings"] for row in rows})
            },
            "minimum_radii": {
                name: minimum_decimal(rows, ("radii", name, "decimal"))
                for name in ("a1", "y", "a3", "fracL")
            },
            "minimum_linear_selector_boundary_gap": minimum_decimal(rows, ("child_selectors", "minimum_linear_selector_boundary_gap")),
            "minimum_quadratic_selector_boundary_gap": minimum_decimal(rows, ("child_selectors", "minimum_quadratic_selector_boundary_gap")),
            "minimum_psi_boundary_gap": minimum_decimal(rows, ("q_branches", "minimum_psi_boundary_gap")),
            "minimum_erf_boundary_gap": minimum_decimal(rows, ("q_branches", "minimum_erf_boundary_gap")),
            "minimum_denominator_abs_lower": minimum_decimal(rows, ("analytic_margins", "minimum_denominator_abs_lower")),
            "minimum_stationary_discriminant_lower": minimum_decimal(rows, ("analytic_margins", "minimum_stationary_discriminant_lower")),
            "minimum_frac_endpoint_margin": minimum_decimal(rows, ("analytic_margins", "minimum_frac_endpoint_margin")),
            "minimum_parent_xr_lower": minimum_decimal(rows, ("analytic_margins", "minimum_parent_xr_lower")),
        },
        "rows": rows,
        "next_handoff": {
            "first": (
                "Derive interval coefficient maps from the physical rae parameter into these four q inputs "
                "and prove that each target roster segment lies in one certified box or an explicit finite partition."
            ),
            "second": (
                "Replace source PSI and six-term ERF evaluations by rigorous special-function balls on each "
                "occupied branch, then enclose t1 through t5 and qq without discarding their correlation."
            ),
            "falsification_rule": (
                "Reject a box immediately if a physical path crosses any recorded child, PSI, ERF, sign, "
                "integer, denominator, or stationary-discriminant boundary."
            ),
        },
        "sources": {
            "telemetry": {"path": relative(TELEMETRY), "sha256": file_hash(TELEMETRY)},
            "child_adapter": {"path": relative(CHILD_ADAPTER), "sha256": file_hash(CHILD_ADAPTER)},
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
            "python_flint": {"version": flint.__version__, "flint_version": flint.__FLINT_VERSION__},
        },
        "proof_boundary": (
            "This gate proves sixty-four nonzero-radius local boxes in derived q-input coordinates on which the "
            "recorded integer, parity, orientation, PSI, ERF, loop, denominator, and cubic stationary selectors "
            "remain fixed. It does not prove that a physical coefficient trajectory enters or remains in any box, "
            "bound source special-function approximation errors, enclose analytic W1-W5 or qq, control the outer "
            "Hardy representation, evaluate a physical carrier, prove a determinant sign, Lambda<=0, PF-infinity, "
            "RH, or a prize-level conclusion."
        ),
    }


def render_note(artifact: dict[str, Any]) -> str:
    aggregate = artifact["aggregate"]
    return "\n".join(
        [
            "# Hardy q Selector-Cell Gate",
            "",
            "Date: 2026-08-06",
            "Status: rigorous local selector cells; not a proof and physical coverage/analytic q open",
            "",
            "## Cell Contract",
            "",
            "Each cell is a nonzero-radius box in `(a1, xr=2*a2, a3, fracL)` with fixed discrete labels. Arb interval evaluation preserves the transformed-child integer/parity selectors, q sign and loop branches, all six `PSI` branch paths, both `ERF` branch paths, denominator nonvanishing, and the cubic stationary branch.",
            "",
            f"All `{aggregate['cell_count']}` saved chains admit such a cell. The largest shrink count is `{aggregate['maximum_halvings']}`; the halving histogram is `{json.dumps(aggregate['halving_histogram'], sort_keys=True)}`.",
            "",
            "## Tightest Certified Margins",
            "",
            f"- Child linear selector: `{aggregate['minimum_linear_selector_boundary_gap']}`",
            f"- Child quadratic selector: `{aggregate['minimum_quadratic_selector_boundary_gap']}`",
            f"- `PSI` branch boundary: `{aggregate['minimum_psi_boundary_gap']}`",
            f"- `ERF` branch boundary: `{aggregate['minimum_erf_boundary_gap']}`",
            f"- Kernel denominator: `{aggregate['minimum_denominator_abs_lower']}`",
            f"- Stationary discriminant: `{aggregate['minimum_stationary_discriminant_lower']}`",
            "",
            "## Critical Scope Boundary",
            "",
            artifact["cell_contract"]["not_preserved_or_proved"],
            "",
            "The next proof obligation is to map actual `rae` intervals into these boxes. Only then is it meaningful to replace the source `PSI` and six-term `ERF` approximations with rigorous special-function balls and enclose the correlated sum `q=t3+t5+conjg(t1+t2+t4)`.",
            "",
            "## Proof Boundary",
            "",
            artifact["proof_boundary"],
            "",
        ]
    )


def main() -> int:
    artifact = build_artifact()
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    NOTE.parent.mkdir(parents=True, exist_ok=True)
    result_tmp = RESULT.with_suffix(".json.tmp")
    note_tmp = NOTE.with_suffix(".md.tmp")
    result_tmp.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    note_tmp.write_text(render_note(artifact), encoding="utf-8")
    result_tmp.replace(RESULT)
    note_tmp.replace(NOTE)
    print(
        "built Hardy q selector-cell gate: "
        f"{artifact['aggregate']['cell_count']} cells, maximum "
        f"{artifact['aggregate']['maximum_halvings']} halvings, "
        "6 PSI and 2 ERF paths per cell"
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
