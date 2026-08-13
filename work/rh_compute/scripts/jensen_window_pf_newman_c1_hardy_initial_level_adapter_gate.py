#!/usr/bin/env python3
"""Build the exact first-level MGS coefficient/parity adapter gate."""

from __future__ import annotations

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
POINT_BALL = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_point_ball_defect_gate.json"
)
RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_initial_level_adapter_gate.json"
)
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_initial_level_adapter_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = REPO_ROOT / "work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_initial_level_adapter_gate.py"
SOURCE_SHA256 = "0fb64f090c21b185f27edc1c9194254cf471e74bfd6af3b88cb7f0cb40ec0c3d"
EXPECTED_TPM_HEX = "C001921FB54442D18469898CC51701B8"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def fraction_record(value: Fraction) -> dict[str, int]:
    return {"numerator": value.numerator, "denominator": value.denominator}


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


def binary128_arb(payload: str) -> arb:
    value = binary128_fraction(payload)
    denominator_exponent = value.denominator.bit_length() - 1
    require(value.denominator == 2**denominator_exponent, "binary128 denominator is not dyadic")
    return arb((value.numerator, -denominator_exponent))


def nearest_integer(value: Fraction) -> int:
    if value >= 0:
        return math.floor(value + Fraction(1, 2))
    return math.ceil(value - Fraction(1, 2))


def preliminary_coefficients(initial: list[Fraction], ix: int, xr: Fraction) -> tuple[list[Fraction], int]:
    sigma = 1 if initial[0] >= 0 else -1
    result = list(initial)
    result[1] = xr / 2
    if ix % 2:
        result[0] = initial[0] - Fraction(sigma, 2)
    return result, sigma


def phase_integer(ix: int, sigma: int, index: int) -> Fraction:
    if ix % 2:
        return Fraction(ix * index * index + sigma * index, 2)
    return Fraction(ix * index * index, 2)


def source_locations() -> dict[str, int]:
    lines = SOURCE.read_text(encoding="utf-8").splitlines()
    tokens = {
        "initial_length": "L(0)=N",
        "twice_quadratic": "x0=2*initialcoeff(2)",
        "nearest_integer": "ix=nint(x0)",
        "quadratic_remainder": "phicoeff(2,1)=xr(1)/2.0",
        "odd_parity": "if (mod(ix,2).eq.1) then",
        "positive_remainder_conjugation": "if (xr(1).gt.0.0) then",
        "negative_remainder_negation": "phicoeff(1,1)=-phicoeff(1,1)",
    }
    locations: dict[str, int] = {}
    for name, token in tokens.items():
        matches = [index for index, line in enumerate(lines, 1) if token in line]
        require(matches, f"source adapter token missing: {token}")
        locations[name] = matches[0]
    return locations


def load_chains() -> dict[int, dict[str, Any]]:
    chains: dict[int, dict[str, Any]] = {}
    for line in TELEMETRY.read_text(encoding="utf-8").splitlines():
        record = json.loads(line)
        if record["type"] == "chain":
            chains[int(record["chain"])] = {"header": record, "level_one": None}
        elif record["type"] == "level" and int(record["level"]) == 1:
            chains[int(record["chain"])]["level_one"] = record
    require(sorted(chains) == list(range(1, 65)), "adapter chain sequence drift")
    require(all(data["level_one"] is not None for data in chains.values()), "level-one record missing")
    return chains


def audit_saved_chain(chain_id: int, data: dict[str, Any]) -> dict[str, Any]:
    header = data["header"]
    level = data["level_one"]
    initial = [binary128_fraction(value) for value in header["initial_coefficients_hex"]]
    logged = [binary128_fraction(value) for value in level["coefficients_hex"]]
    require(len(initial) == len(logged) == int(header["base_degree"]), "saved degree drift")

    conjugates = bool(level["conjugate"])
    preliminary = logged if conjugates else [-value for value in logged]
    xr_abs = binary128_fraction(level["xr_hex"])
    xr = xr_abs if conjugates else -xr_abs
    ix_fraction = 2 * initial[1] - xr
    require(ix_fraction.denominator == 1, f"chain {chain_id} nonintegral ix")
    ix = int(ix_fraction)
    require(nearest_integer(2 * initial[1]) == ix, f"chain {chain_id} NINT mismatch")

    expected, sigma = preliminary_coefficients(initial, ix, xr)
    require(preliminary == expected, f"chain {chain_id} level-one coefficient map mismatch")
    require((xr > 0) == conjugates, f"chain {chain_id} orientation mismatch")

    phase_values = [phase_integer(ix, sigma, index) for index in range(int(header["initial_length"]) + 1)]
    require(all(value.denominator == 1 for value in phase_values), f"chain {chain_id} parity failure")
    max_phase = max((abs(value) for value in phase_values), default=Fraction(0))
    return {
        "chain": chain_id,
        "branch": int(header["branch"]),
        "length": int(header["initial_length"]),
        "ix": ix,
        "sigma": sigma,
        "xr_sign": 1 if xr > 0 else (-1 if xr < 0 else 0),
        "source_conjugates": conjugates,
        "initial_coefficients_hex": header["initial_coefficients_hex"],
        "level_one_coefficients_hex": level["coefficients_hex"],
        "integer_phase_count": len(phase_values),
        "maximum_integer_phase_shift": fraction_record(max_phase),
        "mathematical_adapter_exact": True,
        "binary_tpp_periodicity_bound_zero": max_phase == 0,
    }


def synthetic_parity_fixtures() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for ix in range(-4, 5):
        for sigma in (-1, 1):
            for xr in (Fraction(-1, 4), Fraction(1, 4)):
                a1 = Fraction(sigma, 4)
                initial = [a1, Fraction(ix, 2) + xr / 2, Fraction(1, 7)]
                require(nearest_integer(2 * initial[1]) == ix, "synthetic NINT mismatch")
                transformed, observed_sigma = preliminary_coefficients(initial, ix, xr)
                require(observed_sigma == sigma, "synthetic sign mismatch")
                phase_values = [phase_integer(ix, sigma, index) for index in range(65)]
                require(all(value.denominator == 1 for value in phase_values), "synthetic parity failure")
                for index, phase in enumerate(phase_values):
                    original_value = sum(initial[j] * (index ** (j + 1)) for j in range(3))
                    transformed_value = sum(
                        transformed[j] * (index ** (j + 1)) for j in range(3)
                    )
                    require(original_value - transformed_value == phase, "synthetic phase identity failure")
                rows.append(
                    {
                        "ix": ix,
                        "sigma": sigma,
                        "xr": fraction_record(xr),
                        "checked_indices": 65,
                        "maximum_integer_phase_shift": fraction_record(
                            max(abs(value) for value in phase_values)
                        ),
                    }
                )
    return rows


def arb_record(value: arb) -> dict[str, Any]:
    lower = value.lower()
    upper = value.upper()
    lower_m, lower_e = lower.man_exp()
    upper_m, upper_e = upper.man_exp()
    return {
        "display": value.str(55),
        "lower_dyadic": [int(lower_m), int(lower_e)],
        "upper_dyadic": [int(upper_m), int(upper_e)],
        "lower_decimal": lower.str(55, radius=False),
        "upper_decimal": upper.str(55, radius=False),
    }


def build_artifact() -> dict[str, Any]:
    for path in (SOURCE, TELEMETRY, POINT_BALL, CHECKER):
        require(path.is_file(), f"missing initial-adapter artifact: {path}")
    require(file_hash(SOURCE) == SOURCE_SHA256, "accepted source hash drift")
    require(flint.__version__ == "0.8.0", "python-flint version drift")
    point_ball = json.loads(POINT_BALL.read_text(encoding="utf-8"))
    require(point_ball["aggregate"]["defect_zero_exclusion_count"] == 64, "point-ball dependency drift")

    chains = load_chains()
    saved_rows = [audit_saved_chain(chain_id, chains[chain_id]) for chain_id in sorted(chains)]
    synthetic_rows = synthetic_parity_fixtures()
    ix_histogram: dict[str, int] = {}
    for row in saved_rows:
        key = str(row["ix"])
        ix_histogram[key] = ix_histogram.get(key, 0) + 1

    tpm_payloads = {data["header"]["tpm_hex"] for data in chains.values()}
    require(tpm_payloads == {EXPECTED_TPM_HEX}, "saved tpm payload drift")
    ctx.dps = 160
    ctx.threads = 1
    tpp_hat = -binary128_arb(EXPECTED_TPM_HEX)
    phase_normalization_delta = tpp_hat - 2 * arb.pi()
    delta_abs = abs(phase_normalization_delta)

    return {
        "kind": "jensen_window_pf_newman_c1_hardy_initial_level_adapter_gate",
        "status": "exact_mathematical_adapter_with_open_binary_phase_budget",
        "source": {
            "path": relative(SOURCE),
            "sha256": file_hash(SOURCE),
            "locations": source_locations(),
        },
        "exact_lemma": {
            "original_phase": "A(k)=sum_(j=1)^m a_j k^j",
            "reduced_phase": "B(k)=sum_(j=1)^m b_j k^j",
            "quadratic_split": "2*a_2=ix+xr with ix=nint(2*a_2) and |xr|<=1/2",
            "linear_shift": "b_1=a_1 for even ix; b_1=a_1-sigma/2 for odd ix, sigma=sign_nonnegative(a_1)",
            "phase_integer": "A(k)-B(k)=ix*k^2/2 for even ix and (ix*k^2+sigma*k)/2 for odd ix",
            "parity_reason": "k^2+k and k^2-k are even for every integer k",
            "orientation": "If xr>0, conjugate the negative-phase sum; if xr<=0, negate B before the negative-phase sum.",
            "conclusion": "With mathematical 2*pi, the adapted level-one finite sum equals the original finite generalized Gauss sum exactly.",
        },
        "saved_fixture": {
            "chain_count": len(saved_rows),
            "branch_counts": {
                str(branch): sum(row["branch"] == branch for row in saved_rows)
                for branch in (1, 2)
            },
            "ix_histogram": ix_histogram,
            "coefficient_map_exact_count": len(saved_rows),
            "integer_phase_checks": sum(row["integer_phase_count"] for row in saved_rows),
            "zero_integer_phase_shift_count": sum(
                row["binary_tpp_periodicity_bound_zero"] for row in saved_rows
            ),
            "rows": saved_rows,
        },
        "synthetic_parity_coverage": {
            "fixture_count": len(synthetic_rows),
            "ix_values": list(range(-4, 5)),
            "a1_signs": [-1, 1],
            "xr_signs": [-1, 1],
            "checked_indices_per_fixture": 65,
            "rows": synthetic_rows,
        },
        "binary_phase_normalization": {
            "tpm_hex": EXPECTED_TPM_HEX,
            "definition": "tpp_hat=-decode_binary128(tpm_hex)",
            "delta": "Delta=tpp_hat-2*pi",
            "delta_ball": arb_record(phase_normalization_delta),
            "delta_abs_ball": arb_record(delta_abs),
            "generic_bound": "|G_hat(A)-G_hat(B)| <= |Delta| * sum_(k=0)^N |A(k)-B(k)|",
            "saved_bound": "zero because ix=0 on all 64 saved chains, so A(k)-B(k)=0 exactly",
            "physical_requirement": "Record or bound ix and N on every physical branch before using the adapter with binary tpp.",
        },
        "route_conclusion": {
            "closed": "The mathematical first-level coefficient/parity adapter and both source orientation branches.",
            "saved_source_case": "All 64 rows have ix=0, so the binary phase-normalization drift vanishes identically at level one.",
            "open": "Nonzero-ix physical branches require the explicit Delta-weighted phase budget; deeper transformed levels and analytic qq remain open.",
            "next_target": "Prove the level-one-to-level-two transformed-child adapter and build an interval analytic-qq enclosure on one selector-stable cell.",
        },
        "sources": {
            "telemetry": {"path": relative(TELEMETRY), "sha256": file_hash(TELEMETRY)},
            "point_ball": {"path": relative(POINT_BALL), "sha256": file_hash(POINT_BALL)},
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "proof_boundary": (
            "This gate proves the exact mathematical first-level parity/orientation adapter "
            "and validates the source coefficient map on sixty-four saved binary128 chains. "
            "It exposes, but does not uniformly close, the binary phase-normalization budget "
            "when ix is nonzero. It does not prove the deeper transformed-child adapter, "
            "analytic W1-W5 or special-function bounds, a real-cell or disk enclosure, "
            "transition exclusion, outer Hardy representation control, physical-height "
            "complexity, a carrier value, Lambda<=0, PF-infinity, RH, or a prize-level conclusion."
        ),
    }


def render_note(artifact: dict[str, Any]) -> str:
    saved = artifact["saved_fixture"]
    binary = artifact["binary_phase_normalization"]
    return "\n".join(
        [
            "# Hardy Initial-Level Parity Adapter Gate",
            "",
            "Date: 2026-08-06",
            "Status: exact mathematical level-one adapter; not a proof and physical nonzero-ix binary phase budget open",
            "",
            "## Exact Lemma",
            "",
            "Write `2*a2 = ix + xr`, where `ix=nint(2*a2)`. For even `ix`, the reduced phase differs by `ix*k^2/2`. For odd `ix`, the source shifts the linear coefficient by `-sigma/2`, and the difference is `(ix*k^2+sigma*k)/2`. This is an integer because `k^2+k` and `k^2-k` are even.",
            "",
            "The source's conjugation branch for positive `xr` and coefficient-negation branch for nonpositive `xr` both convert its negative-phase internal sum back to the same positive reduced phase. With mathematical `2*pi`, the adapted level-one sum is therefore exactly the original generalized Gauss sum.",
            "",
            "## Saved Chains",
            "",
            f"All `{saved['chain_count']}` saved chains reproduce the exact source coefficient map. The checker performs `{saved['integer_phase_checks']}` integer-phase checks. Their `ix` histogram is `{json.dumps(saved['ix_histogram'], sort_keys=True)}`.",
            "",
            "Because every saved row has `ix=0`, its reduced polynomial equals the original polynomial exactly. The source's binary approximation to `2*pi` therefore creates no level-one periodicity drift on this fixture.",
            "",
            "## New Physical Budget",
            "",
            f"The source binary constant satisfies `Delta=tpp_hat-2*pi` with `{binary['delta_ball']['display']}`. For nonzero `ix`, a safe generic bound is `{binary['generic_bound']}`. Physical chains must record or bound both `ix` and `N`; mathematical periodicity alone is not a bit-level source bound.",
            "",
            "## Next Target",
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
        "built Hardy initial-level adapter gate: 64 exact source maps, "
        "6720 saved parity checks, 36 synthetic sign/parity fixtures and 1 open binary phase budget"
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
