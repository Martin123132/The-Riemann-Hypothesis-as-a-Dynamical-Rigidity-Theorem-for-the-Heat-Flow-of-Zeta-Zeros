#!/usr/bin/env python3
"""Independently check the sparse joined source-owned-carrier pilot."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
import math
import os
from pathlib import Path
import sys
from typing import Any


for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))
from mpmath import mp


STEM = (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_"
    "non_A_source_owned_carrier_sparse_pilot_gate"
)
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
ROWS = REPO_ROOT / f"work/rh_compute/results/{STEM}_rows.jsonl"
BUILDER = Path(__file__).with_name(Path(__file__).name.removeprefix("check_"))

A = 159_577
B = 5_122_421
L = 2_481_422
HEIGHT = 10_000_000_000
FIRST = 622
ORDINARY_LAST = 39_852
A_WINDOW_FIRST = 39_853
A_WINDOW_LAST = 39_936
SOURCE_MASS = (L + 1) * (A + B) // 2
CHECK_DPS = 90
CHECK_LABELS = (
    "first_B_entry",
    "closest_pair_B_20223",
    "quarter",
    "midpoint",
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def mp_fraction(value: Fraction) -> Any:
    return mp.mpf(value.numerator) / value.denominator


def mp_cis_pi_fraction(value: Fraction) -> Any:
    reduced = value % 2
    return mp.exp(mp.j * mp.pi * reduced.numerator / reduced.denominator)


def float_cis_pi_fraction(value: Fraction) -> complex:
    reduced = value % 2
    angle = math.pi * reduced.numerator / reduced.denominator
    return complex(math.cos(angle), math.sin(angle))


def parse_complex(record: dict[str, str]) -> Any:
    return mp.mpc(record["real"], record["imag"])


def direct_lower_fresnel(q: Any) -> Any:
    return mp.fresnelc(q) + mp.j * mp.fresnels(q) + (1 + mp.j) / 2


def independent_source(x: Fraction, block_size: int = 193) -> complex:
    """Replay the endpoint-complete source with a different float block layout."""

    second_ratio = float_cis_pi_fraction(2 * x)
    block_s0_real: list[float] = []
    block_s0_imag: list[float] = []
    block_s1_real: list[float] = []
    block_s1_imag: list[float] = []
    length = L + 1
    for start in range(0, length, block_size):
        stop = min(length, start + block_size)
        phase = float_cis_pi_fraction(x * (start * start + A * start))
        ratio = float_cis_pi_fraction(x * (A + 2 * start + 1))
        local_s0_real: list[float] = []
        local_s0_imag: list[float] = []
        local_s1_real: list[float] = []
        local_s1_imag: list[float] = []
        for index in range(start, stop):
            local_s0_real.append(phase.real)
            local_s0_imag.append(phase.imag)
            local_s1_real.append(index * phase.real)
            local_s1_imag.append(index * phase.imag)
            phase *= ratio
            ratio *= second_ratio
        block_s0_real.append(math.fsum(local_s0_real))
        block_s0_imag.append(math.fsum(local_s0_imag))
        block_s1_real.append(math.fsum(local_s1_real))
        block_s1_imag.append(math.fsum(local_s1_imag))
    s0 = complex(math.fsum(block_s0_real), math.fsum(block_s0_imag))
    s1 = complex(math.fsum(block_s1_real), math.fsum(block_s1_imag))
    base = float_cis_pi_fraction(x * A * A / 4)
    return base * (A * s0 + 2 * s1)


def independent_carrier(x: Fraction) -> Any:
    """Use direct rational phases, descending modes, and mpmath Fresnel functions."""

    x_mp = mp_fraction(x)
    scale = mp.sqrt(2 / x_mp) / x_mp
    total = 1 + mp.j
    endpoint_a = mp_cis_pi_fraction(x * A * A / 4)

    p_terms: list[Any] = []
    for mode in range(A_WINDOW_LAST, FIRST - 1, -1):
        phase = mp_cis_pi_fraction(Fraction(mode, 1) - Fraction(mode * mode, 1) / x)
        p_terms.append(phase * mode * scale * total)
    full_p = mp.fsum(p_terms)

    a_pair_terms: list[Any] = []
    for mode in range(A_WINDOW_LAST, A_WINDOW_FIRST - 1, -1):
        phase = mp_cis_pi_fraction(Fraction(mode, 1) - Fraction(mode * mode, 1) / x)
        coefficient = phase * mode * scale
        q_minus = mp.sqrt(x_mp / 2) * (A - 2 * mode / x_mp)
        q_plus = mp.sqrt(x_mp / 2) * (A + 2 * mode / x_mp)
        a_pair_terms.append(
            -2 * endpoint_a / (mp.j * mp.pi * x_mp)
            + coefficient * (direct_lower_fresnel(q_plus) - direct_lower_fresnel(q_minus))
        )
    return -full_p - mp.fsum(a_pair_terms)


def physical_density(x: Fraction, value: Any) -> Any:
    x_mp = mp_fraction(x)
    normalization = 2 * (mp.pi / (32 * HEIGHT)) ** mp.mpf("0.25")
    weight = mp.exp(mp.j * HEIGHT / 2 * mp.log((1 - x_mp) / x_mp))
    weight /= (x_mp * (1 - x_mp)) ** mp.mpf("0.25")
    return normalization * mp.re(mp.exp(-mp.j * mp.pi / 8) * weight * value)


def main() -> int:
    mp.dps = CHECK_DPS
    artifact = load_json(RESULT)
    require(artifact.get("passed") is True, "builder artifact did not pass")
    require(artifact["kind"] == STEM, "artifact kind drift")
    require(artifact["decision"]["common_pointwise_mathematical_level_respected"] is True, "level audit failed")
    require(artifact["decision"]["physical_x_quadrature_completed"] is False, "pilot overclaimed quadrature")
    require(artifact["decision"]["non_A_bound_proved"] is False, "pilot overclaimed non-A bound")
    require(artifact["decision"]["rh_implication"] is False, "pilot overclaimed RH")
    require(artifact["decision"]["naive_value_space_panel_route_selected"] is False, "raw panel route guard lost")
    require(
        artifact["decision"]["physical_transform_ownership_ledger_selected_next"] is True,
        "physical-transform ledger route not selected",
    )

    require(file_hash(BUILDER) == artifact["sources"]["builder"]["sha256"], "builder hash drift")
    require(file_hash(Path(__file__)) == artifact["sources"]["checker"]["sha256"], "checker hash drift")
    require(file_hash(ROWS) == artifact["sources"]["rows"]["sha256"], "row-cache hash drift")
    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(file_hash(path) == dependency["sha256"], f"dependency hash drift: {path}")

    rows = artifact["rows"]
    require(len(rows) == 18, "pilot row count drift")
    require(all(row["passed"] is True for row in rows), "stored row failure")
    require(
        [Fraction(row["x"]) for row in rows] == sorted(Fraction(row["x"]) for row in rows),
        "pilot rows are not ordered",
    )
    require(len({row["config_signature"] for row in rows}) == 1, "mixed row-cache signatures")

    jsonl_rows = [json.loads(line) for line in ROWS.read_text(encoding="utf-8").splitlines() if line.strip()]
    signature = rows[0]["config_signature"]
    matching = {row["label"]: row for row in jsonl_rows if row.get("config_signature") == signature}
    require(set(matching) >= {row["label"] for row in rows}, "cache is missing a promoted row")

    # Independent coefficient arithmetic in the order P_m, A_m, A_-m, B_m.
    abel = (-1, -1, 0, -1)
    ordinary_endpoint = (0, 1, 0, 1)
    high_endpoint = (0, 0, -1, 1)
    require(tuple(a + b for a, b in zip(abel, ordinary_endpoint)) == (-1, 0, 0, 0), "ordinary vector drift")
    require(tuple(a + b for a, b in zip(abel, high_endpoint)) == (-1, -1, -1, 0), "A-window vector drift")

    by_label = {row["label"]: row for row in rows}
    for label in CHECK_LABELS:
        row = by_label[label]
        x = Fraction(row["x"])
        stored_source = parse_complex(row["source"]["value"])
        checked_source_complex = independent_source(x)
        checked_source = mp.mpc(repr(checked_source_complex.real), repr(checked_source_complex.imag))
        require(
            abs(checked_source - stored_source) / SOURCE_MASS <= mp.mpf("2e-10"),
            f"independent source recurrence mismatch at {label}",
        )

        checked_carrier = independent_carrier(x)
        stored_carrier = parse_complex(row["selected_representation"]["joined_carrier"])
        carrier_scale = max(mp.mpf(1), abs(checked_carrier), abs(stored_carrier))
        require(
            abs(checked_carrier - stored_carrier) / carrier_scale <= mp.mpf("1e-38"),
            f"independent joined-carrier mismatch at {label}",
        )

        checked_assembly = checked_source + checked_carrier
        stored_density = mp.mpf(row["selected_representation"]["normalized_physical_density"])
        checked_density = physical_density(x, checked_assembly)
        density_tolerance = mp.mpf("2e-10") * SOURCE_MASS
        density_tolerance *= 2 * (mp.pi / (32 * HEIGHT)) ** mp.mpf("0.25")
        density_tolerance /= (mp_fraction(x) * (1 - mp_fraction(x))) ** mp.mpf("0.25")
        require(abs(checked_density - stored_density) <= density_tolerance, f"physical density mismatch at {label}")

    altered = [row["altered_endpoint_Fresnel_cross_check"] for row in rows if row["altered_endpoint_Fresnel_cross_check"]]
    require(len(altered) == 8, "altered cross-check count drift")
    require(all(item["passed"] is True for item in altered), "altered representation failure")
    require(max(mp.mpf(item["relative_discrepancy"]) for item in altered) <= mp.mpf("1e-42"), "altered discrepancy drift")

    periodic = [row["periodic_source_cross_check"] for row in rows if row["periodic_source_cross_check"]]
    require(len(periodic) == 9, "periodic source check count drift")
    require(all(item["passed"] is True for item in periodic), "periodic source check failure")
    require(
        artifact["mathematical_level_audit"]["forbidden_mixing_avoided"] is True,
        "pointwise/integrated level guard lost",
    )
    print("independently checked sparse joined source-owned-carrier pilot", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
