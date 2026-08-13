#!/usr/bin/env python3
"""Build the exact cubic level-one-to-level-two Legendre adapter gate."""

from __future__ import annotations

from decimal import Decimal, localcontext
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
    from flint import acb, arb, ctx
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
INITIAL_ADAPTER = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_initial_level_adapter_gate.json"
)
RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_cubic_child_legendre_adapter_gate.json"
)
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_cubic_child_legendre_adapter_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = (
    REPO_ROOT
    / "work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_cubic_child_legendre_adapter_gate.py"
)
SOURCE_SHA256 = "0fb64f090c21b185f27edc1c9194254cf471e74bfd6af3b88cb7f0cb40ec0c3d"
EXPECTED_TPM_HEX = "C001921FB54442D18469898CC51701B8"
PRECISION_LADDER = (96, 160)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def fraction_record(value: Fraction) -> dict[str, int]:
    return {"numerator": value.numerator, "denominator": value.denominator}


def decimal_text(value: Fraction, digits: int = 55) -> str:
    with localcontext() as context:
        context.prec = digits
        return format(Decimal(value.numerator) / Decimal(value.denominator), ".45E")


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
        return Fraction(sign * mantissa * (2**binary_exponent), 1)
    return Fraction(sign * mantissa, 2 ** (-binary_exponent))


def binary128_ball(payload: str) -> arb:
    value = binary128_fraction(payload)
    exponent = value.denominator.bit_length() - 1
    require(value.denominator == 2**exponent, "binary128 value is not dyadic")
    return arb((value.numerator, -exponent))


def binary128_complex(payload: list[str]) -> acb:
    require(isinstance(payload, list) and len(payload) == 2, "invalid complex payload")
    return acb(binary128_ball(payload[0]), binary128_ball(payload[1]))


def fraction_ball(value: Fraction) -> arb:
    return arb(value.numerator) / arb(value.denominator)


def nearest_integer(value: Fraction) -> int:
    if value >= 0:
        return math.floor(value + Fraction(1, 2))
    return math.ceil(value - Fraction(1, 2))


def dyadic(value: arb) -> list[int]:
    mantissa, exponent = value.man_exp()
    return [int(mantissa), int(exponent)]


def arb_record(value: arb, digits: int = 48) -> dict[str, Any]:
    lower = value.lower()
    upper = value.upper()
    return {
        "display": value.str(digits),
        "lower_dyadic": dyadic(lower),
        "upper_dyadic": dyadic(upper),
        "lower_decimal": lower.str(digits, radius=False),
        "upper_decimal": upper.str(digits, radius=False),
    }


def acb_record(value: acb, digits: int = 44) -> dict[str, Any]:
    return {
        "display": value.str(digits),
        "real": arb_record(value.real, digits),
        "imag": arb_record(value.imag, digits),
        "abs": arb_record(abs(value), digits),
        "contains_zero": value.contains(0),
        "relative_accuracy_bits": int(value.rel_accuracy_bits()),
    }


def source_locations() -> dict[str, int]:
    lines = SOURCE.read_text(encoding="utf-8").splitlines()
    tokens = {
        "scheme_quadratic_scale": "y=2*phicoeff(2,j)",
        "legendre_quadratic": "see(1)=1/(2*y)",
        "legendre_cubic": "see(2)=-phicoeff(3,j)/(y**3)",
        "first_omitted_quartic": "small=9*(phicoeff(3,j)**2)/(2*(y**5))",
        "translation_center": "x=-phicoeff(1,j)",
        "binomial_translation": "phicoeff(i,j+1)=sum",
        "raised_order_zero_extension": "phicoeff(m1(j),j-1)=0.0",
        "raise_test": "smalle=(x0**(m1(j)+1))*small",
        "dual_quadratic_integer": "ix=nint(2*phicoeff(2,j))",
        "dual_quadratic_remainder": "phicoeff(2,j)=0.5*xr(j)",
        "dual_linear_mod_one": "phicoeff(1,j)=phicoeff(1,j)-anint(phicoeff(1,j),dp1)",
        "dual_odd_parity": "if (mod(abs(ix),2).eq.1) then",
        "dual_orientation": "if (xr(j).gt.0.0) then",
        "kernel_h1": "h1=3*phicoeff(3,MIT-1)/(xr(MIT-1)**2)",
        "kernel_denominator": "denn=1-phicoeff(1,MIT-1)*(h1-phicoeff(1,MIT-1)*(6*h2-1.5*h3))",
        "kernel_multiplier": "c1=exp(pp)/(EPI4*sqrt(x))",
    }
    locations: dict[str, int] = {}
    for name, token in tokens.items():
        matches = [index for index, line in enumerate(lines, 1) if token in line]
        require(matches, f"source child-adapter token missing: {token}")
        locations[name] = matches[0]
    return locations


def load_chains() -> dict[int, dict[str, Any]]:
    chains: dict[int, dict[str, Any]] = {}
    for line in TELEMETRY.read_text(encoding="utf-8").splitlines():
        record = json.loads(line)
        chain_id = int(record.get("chain", 0))
        if record["type"] == "chain":
            chains[chain_id] = {"header": record, "levels": {}, "steps": {}}
        elif record["type"] == "level":
            chains[chain_id]["levels"][int(record["level"])] = record
        elif record["type"] == "recurrence":
            chains[chain_id]["steps"][int(record["nit"])] = record
    require(sorted(chains) == list(range(1, 65)), "child-adapter chain sequence drift")
    for chain_id, data in chains.items():
        require(data["header"]["mit"] == 2, f"chain {chain_id} MIT drift")
        require(sorted(data["levels"]) == [1, 2], f"chain {chain_id} level gap")
        require(sorted(data["steps"]) == [1], f"chain {chain_id} recurrence gap")
    return chains


def cubic_legendre_raw(parent: list[Fraction]) -> tuple[list[Fraction], Fraction]:
    require(len(parent) == 3, "cubic parent required")
    a1, a2, a3 = parent
    y = 2 * a2
    require(y != 0, "degenerate quadratic coefficient")
    s1 = Fraction(1, 2) / y
    s2 = -a3 / (y**3)
    x = -a1
    raw = [
        s1 * x**2 + s2 * x**3,
        2 * s1 * x + 3 * s2 * x**2,
        s1 + 3 * s2 * x,
        s2,
    ]
    h4 = 9 * a3**2 / (2 * y**5)
    return raw, h4


def normalize_raw(raw: list[Fraction]) -> dict[str, Any]:
    require(len(raw) == 4, "cubic raw coefficient count drift")
    n1 = nearest_integer(raw[1])
    linear_residual = raw[1] - n1
    ix = nearest_integer(2 * raw[2])
    xr = 2 * raw[2] - ix
    sigma = 1 if linear_residual >= 0 else -1
    linear = linear_residual
    if ix % 2:
        linear -= Fraction(sigma, 2)
    pre_orientation = [linear, xr / 2, raw[3]]
    orientation = 1 if xr > 0 else -1
    stored = [orientation * value for value in pre_orientation]
    return {
        "n1": n1,
        "ix": ix,
        "xr": xr,
        "sigma": sigma,
        "conjugate": xr > 0,
        "pre_orientation": pre_orientation,
        "stored": stored,
        "linear_selector_margin": Fraction(1, 2) - abs(linear_residual),
        "quadratic_selector_margin": Fraction(1, 2) - abs(xr),
        "linear_sign_margin": abs(linear_residual),
        "orientation_margin": abs(xr),
    }


def reconstruct_source_raw(
    expected_raw: list[Fraction], level_two: dict[str, Any]
) -> tuple[list[Fraction], dict[str, Any]]:
    stored = [binary128_fraction(value) for value in level_two["coefficients_hex"]]
    require(len(stored) == 3, "saved child is not cubic")
    orientation = 1 if level_two["conjugate"] else -1
    pre = [orientation * value for value in stored]
    xr_logged = orientation * binary128_fraction(level_two["xr_hex"])
    require(2 * pre[1] == xr_logged, "logged child xr/coefficient mismatch")

    exact = normalize_raw(expected_raw)
    ix = int(exact["ix"])
    n1 = int(exact["n1"])
    sigma = int(exact["sigma"])
    source_raw = [
        binary128_fraction(level_two["phi0_hex"]),
        pre[0] + n1 + (Fraction(sigma, 2) if ix % 2 else 0),
        pre[1] + Fraction(ix, 2),
        pre[2],
    ]
    require(nearest_integer(2 * source_raw[2]) == ix, "source quadratic selector changed")
    require(nearest_integer(source_raw[1]) == n1, "source linear selector changed")
    require(2 * source_raw[2] - ix == xr_logged, "source raw xr reconstruction failed")
    require((xr_logged > 0) == bool(level_two["conjugate"]), "source orientation mismatch")
    return source_raw, exact


def phase(coefficients: list[Fraction], index: int) -> Fraction:
    total = Fraction(0)
    for coefficient in reversed(coefficients):
        total = (total + coefficient) * index
    return total


def denominator_forms(parent: list[Fraction], index: int) -> tuple[Fraction, Fraction]:
    a1, a2, a3 = parent[:3]
    y = 2 * a2
    h1 = 3 * a3 / y**2
    h2 = parent[3] / y**3 if len(parent) > 3 else Fraction(0)
    h3 = h1**2
    denn = 1 - a1 * (h1 - a1 * (6 * h2 - Fraction(3, 2) * h3))
    rminorcor = h1 + 3 * a1 * (h3 - 4 * h2)
    rminorcor2 = Fraction(3, 2) * h3 - 6 * h2
    source = denn + index * (rminorcor - index * rminorcor2)
    u = Fraction(index) - a1
    centered = 1 + h1 * u + (6 * h2 - Fraction(3, 2) * h1**2) * u**2
    return source, centered


def weighted_kernel(
    level_one: dict[str, Any],
    level_two: dict[str, Any],
    step: dict[str, Any],
    tpm_hex: str,
    dps: int,
) -> tuple[acb, acb]:
    ctx.dps = dps
    ctx.threads = 1
    parent = [binary128_fraction(value) for value in level_one["coefficients_hex"]]
    coefficients = [binary128_ball(value) for value in level_two["coefficients_hex"]]
    tpm = binary128_ball(tpm_hex)
    total = acb(0)
    for index in range(int(level_two["length"]) + 1):
        angle = arb(0)
        for coefficient in reversed(coefficients):
            angle = (angle + coefficient) * index
        denominator, centered = denominator_forms(parent, index)
        require(denominator == centered and denominator != 0, "kernel denominator identity failed")
        total += acb(arb(0), tpm * angle).exp() / fraction_ball(denominator)
    if level_two["subtract_one"]:
        total -= 1
    if level_two["conjugate"]:
        total = total.conjugate()
    logged = binary128_complex(step["state_before_hex"])
    return total, total - logged


def exact_legendre_tail(
    parent: list[Fraction], index: int, dps: int
) -> tuple[arb, arb, arb]:
    ctx.dps = dps
    ctx.threads = 1
    a1, a2, a3 = [fraction_ball(value) for value in parent]
    y = 2 * a2
    u = index - a1
    z = 12 * a3 * u / y**2
    root = (1 + z).sqrt()
    stationary_x = 2 * u / (y * (1 + root))
    exact = u * stationary_x - (y / 2) * stationary_x**2 - a3 * stationary_x**3
    truncated = u**2 / (2 * y) - a3 * u**3 / y**3
    return exact - truncated, z, 1 + z


def polynomial_power_coefficient(series: list[Fraction], power: int, order: int) -> Fraction:
    polynomial = [Fraction(0)] * (order + 1)
    polynomial[0] = Fraction(1)
    for _ in range(power):
        product = [Fraction(0)] * (order + 1)
        for left, left_value in enumerate(polynomial):
            if left_value == 0:
                continue
            for right in range(1, min(len(series), order - left + 1)):
                product[left + right] += left_value * series[right]
        polynomial = product
    return polynomial[order]


def formal_legendre_coefficients(
    parent: dict[int, Fraction], maximum_degree: int
) -> dict[int, Fraction]:
    require(parent.get(2, Fraction(0)) != 0, "formal parent has zero quadratic term")
    derivative = {
        degree - 1: degree * coefficient
        for degree, coefficient in parent.items()
        if degree >= 2 and coefficient != 0
    }
    linear = derivative[1]
    inverse = [Fraction(0)] * maximum_degree
    for order in range(1, maximum_degree):
        remainder = Fraction(0)
        for power, coefficient in derivative.items():
            if power >= 2:
                remainder += coefficient * polynomial_power_coefficient(inverse, power, order)
        target = Fraction(1) if order == 1 else Fraction(0)
        inverse[order] = (target - remainder) / linear
    return {degree: inverse[degree - 1] / degree for degree in range(2, maximum_degree + 1)}


def synthetic_normalization_audit() -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    checks = 0
    for a1 in (Fraction(-2, 5), Fraction(-1, 7), Fraction(0), Fraction(1, 7), Fraction(2, 5)):
        for a2 in (
            Fraction(-1, 5),
            Fraction(-3, 20),
            Fraction(1, 8),
            Fraction(3, 20),
            Fraction(1, 5),
        ):
            for a3 in (Fraction(-1, 500), Fraction(0), Fraction(1, 500)):
                raw, _ = cubic_legendre_raw([a1, a2, a3])
                normalized = normalize_raw(raw)
                if normalized["quadratic_selector_margin"] == 0:
                    continue
                for index in range(17):
                    raw_value = raw[0] + phase(raw[1:], index)
                    reduced_value = raw[0] + phase(normalized["pre_orientation"], index)
                    require((raw_value - reduced_value).denominator == 1, "synthetic parity failure")
                    checks += 1
                rows.append(
                    {
                        "a1": fraction_record(a1),
                        "a2": fraction_record(a2),
                        "a3": fraction_record(a3),
                        "ix": normalized["ix"],
                        "n1": normalized["n1"],
                        "sigma": normalized["sigma"],
                        "xr_sign": 1 if normalized["xr"] > 0 else -1,
                        "conjugate": normalized["conjugate"],
                    }
                )
    require(rows, "no synthetic normalization fixtures")
    return {
        "fixture_count": len(rows),
        "integer_phase_checks": checks,
        "ix_parities": sorted({abs(int(row["ix"])) % 2 for row in rows}),
        "xr_signs": sorted({int(row["xr_sign"]) for row in rows}),
        "conjugation_values": sorted({bool(row["conjugate"]) for row in rows}),
        "sample_rows": rows[:12],
    }


def formal_series_audit() -> dict[str, Any]:
    cases = 0
    coefficient_checks = 0
    zero_extension_checks = 0
    for degree in range(3, 8):
        for seed in range(1, 5):
            parent = {1: Fraction((-1) ** seed, 7 + seed), 2: Fraction(2 + seed, 11 + seed)}
            for order in range(3, degree + 1):
                parent[order] = Fraction((-1) ** (order + seed), 40 + 3 * order + seed)
            coefficients = formal_legendre_coefficients(parent, min(10, degree + 3))
            extended = dict(parent)
            extended[degree + 1] = Fraction(0)
            extended_coefficients = formal_legendre_coefficients(extended, min(10, degree + 3))
            require(coefficients == extended_coefficients, "zero-extension prefix failure")
            zero_extension_checks += len(coefficients)

            y = 2 * parent[2]
            require(coefficients[2] == Fraction(1, 2) / y, "formal h2 failure")
            require(coefficients[3] == -parent[3] / y**3, "formal h3 failure")
            expected_h4 = 9 * parent[3] ** 2 / (2 * y**5)
            if degree >= 4:
                expected_h4 -= parent[4] / y**4
            require(coefficients[4] == expected_h4, "formal h4 failure")
            coefficient_checks += 3
            if degree >= 4:
                expected_h5 = -27 * parent[3] ** 3 / y**7 + 12 * parent[3] * parent[4] / y**6
                if degree >= 5:
                    expected_h5 -= parent[5] / y**5
                require(coefficients[5] == expected_h5, "formal h5 failure")
                coefficient_checks += 1
            cases += 1
    return {
        "case_count": cases,
        "coefficient_checks": coefficient_checks,
        "zero_extension_checks": zero_extension_checks,
        "degrees": list(range(3, 8)),
        "conclusion": (
            "The formal inverse/Legendre coefficients are prefix-stable when a raised-order "
            "parent is formed by appending an exact zero coefficient. This proves the mathematical "
            "zero-extension operation; source formulas above quartic still require their own audit."
        ),
    }


def choose_arb(values: list[arb], maximum: bool) -> arb:
    require(values, "empty Arb choice")
    return (max if maximum else min)(values, key=lambda value: float(value.mid()))


def build_artifact() -> dict[str, Any]:
    for path in (SOURCE, TELEMETRY, INITIAL_ADAPTER, CHECKER):
        require(path.is_file(), f"missing child-adapter dependency: {path}")
    require(file_hash(SOURCE) == SOURCE_SHA256, "accepted source hash drift")
    require(flint.__version__ == "0.8.0", "python-flint version drift")
    initial = json.loads(INITIAL_ADAPTER.read_text(encoding="utf-8"))
    require(initial["saved_fixture"]["coefficient_map_exact_count"] == 64, "initial adapter drift")

    chains = load_chains()
    rows: list[dict[str, Any]] = []
    all_gaps: list[Fraction] = []
    truncation_indicators: list[Fraction] = []
    selector_margins: dict[str, list[Fraction]] = {
        "linear": [],
        "quadratic": [],
        "linear_sign": [],
        "orientation": [],
    }
    denominator_values: list[Fraction] = []
    kernel_gaps: list[arb] = []
    legendre_tails: list[arb] = []
    discriminants: list[arb] = []
    z_values: list[arb] = []
    parity_checks = 0
    denominator_checks = 0
    ix_histogram: dict[str, int] = {}

    for chain_id in sorted(chains):
        data = chains[chain_id]
        header = data["header"]
        level_one = data["levels"][1]
        level_two = data["levels"][2]
        step = data["steps"][1]
        require(header["tpm_hex"] == EXPECTED_TPM_HEX, f"chain {chain_id} tpm drift")
        require(int(level_one["degree"]) == int(level_two["degree"]) == 3, "saved degree drift")

        parent = [binary128_fraction(value) for value in level_one["coefficients_hex"]]
        require(2 * parent[1] == binary128_fraction(level_one["xr_hex"]), "parent y/xr drift")
        expected_raw, h4 = cubic_legendre_raw(parent)
        source_raw, normalized = reconstruct_source_raw(expected_raw, level_two)
        gaps = [abs(left - right) for left, right in zip(source_raw, expected_raw)]
        all_gaps.extend(gaps)
        ix_key = str(normalized["ix"])
        ix_histogram[ix_key] = ix_histogram.get(ix_key, 0) + 1

        child_length = int(level_two["length"])
        indicator = abs(h4 * (Fraction(child_length) - parent[0]) ** 4)
        truncation_indicators.append(indicator)
        require(indicator < Fraction(1, 1000), f"chain {chain_id} source raise criterion")

        for name, key in (
            ("linear", "linear_selector_margin"),
            ("quadratic", "quadratic_selector_margin"),
            ("linear_sign", "linear_sign_margin"),
            ("orientation", "orientation_margin"),
        ):
            selector_margins[name].append(normalized[key])
        require(all(values[-1] > 0 for values in selector_margins.values()), "saved selector tie")

        for index in range(child_length + 1):
            raw_value = expected_raw[0] + phase(expected_raw[1:], index)
            reduced_value = expected_raw[0] + phase(normalized["pre_orientation"], index)
            require((raw_value - reduced_value).denominator == 1, "saved child parity failure")
            source_denominator, centered_denominator = denominator_forms(parent, index)
            require(source_denominator == centered_denominator, "saved denominator centering failure")
            require(source_denominator != 0, "saved denominator zero")
            denominator_values.append(abs(source_denominator))
            parity_checks += 1
            denominator_checks += 1

            tail, z_value, discriminant = exact_legendre_tail(parent, index, PRECISION_LADDER[-1])
            legendre_tails.append(abs(tail))
            z_values.append(abs(z_value))
            discriminants.append(discriminant)

        low_kernel, low_gap = weighted_kernel(
            level_one, level_two, step, header["tpm_hex"], PRECISION_LADDER[0]
        )
        high_kernel, high_gap = weighted_kernel(
            level_one, level_two, step, header["tpm_hex"], PRECISION_LADDER[-1]
        )
        require(low_kernel.overlaps(high_kernel), f"chain {chain_id} kernel precision drift")
        kernel_gaps.append(abs(high_gap))

        rows.append(
            {
                "chain": chain_id,
                "branch": int(header["branch"]),
                "parent_length": int(level_one["length"]),
                "child_length": child_length,
                "parent_degree": int(level_one["degree"]),
                "child_degree": int(level_two["degree"]),
                "ix": int(normalized["ix"]),
                "n1": int(normalized["n1"]),
                "sigma": int(normalized["sigma"]),
                "xr_sign": 1 if normalized["xr"] > 0 else -1,
                "source_conjugates_child": bool(level_two["conjugate"]),
                "source_subtracts_one": bool(level_two["subtract_one"]),
                "maximum_exact_coefficient_gap": fraction_record(max(gaps)),
                "maximum_exact_coefficient_gap_decimal": decimal_text(max(gaps)),
                "first_omitted_h4": fraction_record(h4),
                "source_raise_indicator": fraction_record(indicator),
                "source_raise_indicator_decimal": decimal_text(indicator),
                "kernel_source_roundoff_gap": acb_record(high_gap),
            }
        )

    maximum_gap = max(all_gaps)
    maximum_indicator = max(truncation_indicators)
    maximum_kernel_gap = choose_arb([value.upper() for value in kernel_gaps], maximum=True)
    maximum_tail = choose_arb([value.upper() for value in legendre_tails], maximum=True)
    maximum_z = choose_arb([value.upper() for value in z_values], maximum=True)
    minimum_discriminant = choose_arb([value.lower() for value in discriminants], maximum=False)

    synthetic = synthetic_normalization_audit()
    formal = formal_series_audit()
    denominator_synthetic_checks = 0
    for a1 in (Fraction(-2, 5), Fraction(1, 7), Fraction(3, 8)):
        for a2 in (Fraction(1, 8), Fraction(3, 20)):
            for a3 in (Fraction(-1, 500), Fraction(1, 500)):
                for a4 in (Fraction(0), Fraction(1, 700)):
                    parent = [a1, a2, a3, a4]
                    for index in range(-3, 4):
                        source_form, centered_form = denominator_forms(parent, index)
                        require(source_form == centered_form, "synthetic denominator identity failure")
                        denominator_synthetic_checks += 1

    return {
        "kind": "jensen_window_pf_newman_c1_hardy_cubic_child_legendre_adapter_gate",
        "status": "exact_cubic_child_adapter_with_finite_bit_and_analytic_q_tail_open",
        "source": {
            "path": relative(SOURCE),
            "sha256": file_hash(SOURCE),
            "locations": source_locations(),
        },
        "exact_cubic_lemma": {
            "parent_phase": "F(x)=a1*x+a2*x^2+a3*x^3 and y=2*a2",
            "stationary_coordinate": "u=F'(x)-a1=y*x+3*a3*x^2",
            "legendre_dual": "H(u)=u*X(u)-a2*X(u)^2-a3*X(u)^3 with u=y*X+3*a3*X^2",
            "source_truncation": "H3(u)=u^2/(2*y)-a3*u^3/y^3",
            "first_omitted_term": "[u^4]H(u)=9*a3^2/(2*y^5)",
            "translation": "phifactors is the binomial expansion H3(k-a1)=r0+r1*k+r2*k^2+r3*k^3",
            "raw_coefficients": {
                "r0": "s1*x^2+s2*x^3",
                "r1": "2*s1*x+3*s2*x^2",
                "r2": "s1+3*s2*x",
                "r3": "s2",
                "definitions": "x=-a1, s1=1/(2*y), s2=-a3/y^3",
            },
            "integer_normalization": (
                "After r1 is reduced modulo one and 2*r2=ix+xr is split, the raw and "
                "normalized phases differ by the integer-valued expression n1*k+ix*k^2/2 "
                "and, for odd ix, sigma*k/2."
            ),
            "parity_reason": "ix*k^2+sigma*k is even for odd ix because sigma is +1 or -1",
            "orientation": (
                "If xr>0 the stored negative-phase child is conjugated; if xr<=0 its nonconstant "
                "coefficients are negated. Both routes recover the same positive normalized phase."
            ),
            "constant_phase": (
                "r0 is not negated. The recurrence multiplier exp(i*tpp*r0)/(EPI4*sqrt(abs(xr_parent))) "
                "restores the translated Legendre constant."
            ),
        },
        "kernel_denominator_lemma": {
            "definitions": "h1=3*a3/y^2, h2=a4/y^3 (zero for the saved cubic rows), u=k-a1",
            "centered_form": "D(k)=1+h1*u+(6*h2-3*h1^2/2)*u^2",
            "source_form": "D(k)=denn+k*(rminorcor-k*rminorcor2)",
            "conclusion": "The centered and source forms are identical polynomials, including raised-order a4.",
            "saved_checks": denominator_checks,
            "synthetic_raised_order_checks": denominator_synthetic_checks,
            "minimum_saved_absolute_denominator": fraction_record(min(denominator_values)),
            "minimum_saved_absolute_denominator_decimal": decimal_text(min(denominator_values)),
        },
        "saved_fixture": {
            "chain_count": len(rows),
            "branch_counts": {
                str(branch): sum(row["branch"] == branch for row in rows) for branch in (1, 2)
            },
            "mit_histogram": {"2": len(rows)},
            "degree_pair_histogram": {"3_to_3": len(rows)},
            "parent_length_range": [min(row["parent_length"] for row in rows), max(row["parent_length"] for row in rows)],
            "child_length_range": [min(row["child_length"] for row in rows), max(row["child_length"] for row in rows)],
            "ix_histogram": ix_histogram,
            "conjugated_child_count": sum(row["source_conjugates_child"] for row in rows),
            "subtract_one_count": sum(row["source_subtracts_one"] for row in rows),
            "integer_phase_checks": parity_checks,
            "maximum_exact_raw_coefficient_gap": fraction_record(maximum_gap),
            "maximum_exact_raw_coefficient_gap_decimal": decimal_text(maximum_gap),
            "minimum_selector_margins": {
                key: {
                    "exact": fraction_record(min(values)),
                    "decimal": decimal_text(min(values)),
                }
                for key, values in selector_margins.items()
            },
            "source_raise_threshold": fraction_record(Fraction(1, 1000)),
            "maximum_source_raise_indicator": fraction_record(maximum_indicator),
            "maximum_source_raise_indicator_decimal": decimal_text(maximum_indicator),
            "rows": rows,
        },
        "rigorous_point_reconstruction": {
            "precision_ladder_dps": list(PRECISION_LADDER),
            "weighted_kernel_formula": (
                "Apply the logged negative phase and exact centered D(k), subtract one when flagged, "
                "then conjugate when flagged."
            ),
            "maximum_weighted_kernel_source_roundoff_gap_upper": arb_record(maximum_kernel_gap),
            "maximum_exact_legendre_tail_at_saved_child_indices": arb_record(maximum_tail),
            "maximum_dimensionless_branch_parameter_abs": arb_record(maximum_z),
            "minimum_stationary_discriminant": arb_record(minimum_discriminant),
            "scope": (
                "These are rigorous balls at the finitely many saved child indices. They are not a "
                "uniform real-cell or complex-disk enclosure and do not enclose analytic qq."
            ),
        },
        "synthetic_normalization": synthetic,
        "formal_raised_order": formal,
        "route_conclusion": {
            "closed": (
                "The exact cubic formal Legendre transform through degree three, its translation, "
                "integer/parity normalization, orientation, constant phase, and quadratic denominator centering."
            ),
            "saved_source_case": (
                "All 64 MIT=2 chains select the same exact mathematical adapter with positive selector margins; "
                "their finite source coefficient and weighted-kernel gaps are explicitly certified."
            ),
            "new_observation": (
                "The source raise-order test uses transformed length L(1)=1 or 2, not original length L(0)=104. "
                "Its h4*(L(1)-a1)^4 indicator is below 0.001 on every saved chain."
            ),
            "open": (
                "The 0.001 first-omitted-term test is a heuristic, not a tail theorem. Source formulas above "
                "quartic, selector-stable interval transport, analytic W1-W5/qq, and the outer Hardy remainder remain open."
            ),
            "next_target": (
                "Turn the exact cubic stationary branch into an outward-rounded interval cell enclosure, then "
                "bound the full Legendre tail and analytic qq on that same selector-stable domain."
            ),
        },
        "sources": {
            "telemetry": {"path": relative(TELEMETRY), "sha256": file_hash(TELEMETRY)},
            "initial_adapter": {"path": relative(INITIAL_ADAPTER), "sha256": file_hash(INITIAL_ADAPTER)},
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
            "python_flint": {"version": flint.__version__, "flint_version": flint.__FLINT_VERSION__},
        },
        "proof_boundary": (
            "This gate proves the exact mathematical cubic level-one-to-level-two Legendre/parity/orientation "
            "adapter, the exact centered denominator identity, and finite exact-input point reconstructions for "
            "sixty-four saved chains. It does not turn the source's first-omitted-term criterion into a rigorous "
            "tail bound, audit all source formulas above quartic, enclose analytic W1-W5 or qq on a real cell or "
            "disk, control selector transitions or the outer Hardy representation, evaluate a physical carrier, "
            "prove a determinant sign, Lambda<=0, PF-infinity, RH, or a prize-level conclusion."
        ),
    }


def render_note(artifact: dict[str, Any]) -> str:
    saved = artifact["saved_fixture"]
    denominator = artifact["kernel_denominator_lemma"]
    rigorous = artifact["rigorous_point_reconstruction"]
    return "\n".join(
        [
            "# Hardy Cubic Child Legendre Adapter Gate",
            "",
            "Date: 2026-08-06",
            "Status: exact cubic transformed-child adapter; not a proof and analytic q/tail enclosure open",
            "",
            "## Exact Transform",
            "",
            "For `F(x)=a1*x+a2*x^2+a3*x^3`, put `y=2*a2` and solve `u=y*X+3*a3*X^2` on the branch `X(0)=0`. Its Legendre dual satisfies",
            "",
            "`H(u)=u^2/(2*y)-a3*u^3/y^3+[9*a3^2/(2*y^5)]u^4+...`.",
            "",
            "The source's `scheme` constructs the displayed cubic truncation, and `phifactors` is exactly its binomial translation `H3(k-a1)`. The subsequent integer, parity, and orientation operations preserve the mathematical integer-index phase. The untranslated constant is carried by `exp(i*tpp*phi0)` in the recurrence multiplier.",
            "",
            "## Denominator",
            "",
            f"The source kernel denominator is exactly `{denominator['centered_form']}`. The gate performs `{denominator['saved_checks']}` saved checks and `{denominator['synthetic_raised_order_checks']}` synthetic raised-order checks. The minimum saved absolute denominator is `{denominator['minimum_saved_absolute_denominator_decimal']}`.",
            "",
            "## Saved Chains",
            "",
            f"All `{saved['chain_count']}` chains have `MIT=2` and degree pair `3 -> 3`. The original parent has `{saved['parent_length_range']}` terms under the source's inclusive convention, while the transformed child length is only `{saved['child_length_range']}`. The raise-order test therefore uses `L(1)`, not `L(0)`.",
            "",
            f"The largest exact raw-coefficient source gap is `{saved['maximum_exact_raw_coefficient_gap_decimal']}`. The largest first-omitted-term indicator is `{saved['maximum_source_raise_indicator_decimal']}`, below the source threshold `0.001`.",
            "",
            "## Rigorous Point Audit",
            "",
            f"The maximum exact Legendre-tail ball over all saved child indices is `{rigorous['maximum_exact_legendre_tail_at_saved_child_indices']['display']}`. The maximum weighted-kernel/source gap is `{rigorous['maximum_weighted_kernel_source_roundoff_gap_upper']['display']}`. These are finite exact-input point statements, not interval theorems.",
            "",
            "## What Remains",
            "",
            artifact["route_conclusion"]["open"],
            "",
            artifact["route_conclusion"]["next_target"],
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
        "built Hardy cubic child Legendre adapter gate: 64 saved chains, "
        f"{artifact['saved_fixture']['integer_phase_checks']} parity checks, "
        f"{artifact['kernel_denominator_lemma']['saved_checks']} denominator checks"
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
