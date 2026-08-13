#!/usr/bin/env python3
"""Independently check the shifted-Hardy diagnostic bridge gate."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import re

import mpmath as mp
import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = "jensen_window_pf_newman_c1_shifted_hardy_diagnostic_bridge_gate"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_complex(text: str) -> complex:
    match = re.fullmatch(
        r"\(([^ ]+)\s*([+-])\s*([^j]+)j\)", text.strip()
    )
    require(match is not None, f"cannot parse complex value: {text}")
    sign = 1.0 if match.group(2) == "+" else -1.0
    return complex(float(match.group(1)), sign * float(match.group(3)))


def parse_real_part(text: str) -> float:
    match = re.match(r"\(([^ ]+)", text.strip())
    require(match is not None, f"cannot parse real part: {text}")
    return float(match.group(1))


def check_sources(payload: dict) -> int:
    audit = payload["source_audit"]
    require(set(audit) == set(payload["source_sha256"]), "source key drift")
    require(len(audit) == 16, "source artifact count")
    for key, item in audit.items():
        path = REPO_ROOT / item["path"]
        require(path.is_file(), f"missing source artifact {key}")
        digest = file_hash(path)
        require(digest == item["sha256"], f"source hash drift {key}")
        require(digest == payload["source_sha256"][key], f"source map drift {key}")
    return len(audit)


def parse_single(path: Path) -> float:
    values = re.findall(
        r"Grand total of Hardy function Z\(t\)=\s*([-+0-9.Ee]+)",
        path.read_text(encoding="utf-8"),
    )
    require(len(values) == 1, "single raw-output count")
    return float(values[0])


def parse_multi(path: Path) -> np.ndarray:
    values = re.findall(
        r"Grand total of Hardy function Z\(t[^=]*=\s*([-+0-9.Ee]+)",
        path.read_text(encoding="utf-8"),
    )
    require(len(values) == 15, f"multi raw-output count: {path}")
    return np.asarray([float(value) for value in values])


def leading_correction(t_value: mp.mpf) -> float:
    a_value = mp.sqrt(t_value / (2 * mp.pi))
    n_value = int(mp.floor(a_value))
    fraction = a_value - n_value
    numerator = mp.cos(
        2 * mp.pi * (fraction * (fraction - 1) - mp.mpf(1) / 16)
    )
    denominator = mp.cos(2 * mp.pi * fraction)
    return float(
        (-1) ** (n_value - 1)
        * a_value ** (-mp.mpf("0.5"))
        * numerator
        / denominator
    )


def independent_central_main(t_value: mp.mpf) -> tuple[int, float]:
    n_value = int(mp.floor(mp.sqrt(t_value / (2 * mp.pi))))
    total = 0j
    for integer in range(1, n_value + 1):
        total += complex(mp.e ** (1j * t_value * mp.log(integer))) / math.sqrt(
            integer
        )
    phase = complex(mp.e ** (-1j * mp.siegeltheta(t_value)))
    return n_value, 2.0 * float(np.real(phase * total))


def check_fixture_values(payload: dict) -> int:
    audited = payload["source_audit"]
    single_path = REPO_ROOT / audited["fixture:single_t1e10_et05"]["path"]
    multi_path = REPO_ROOT / audited["fixture:multi_t1e10_et005"]["path"]
    h002_path = REPO_ROOT / audited["fixture:multi_t1e10_et005_h002"]["path"]
    h004_path = REPO_ROOT / audited["fixture:multi_t1e10_et005_h004"]["path"]
    t12_path = REPO_ROOT / audited["fixture:multi_t1e12_et005"]["path"]

    single = parse_single(single_path)
    exact = mp.siegelz(mp.mpf("1e10"))
    stored = payload["calibration"]["single_t1e10"]
    require(abs(single - stored["fast_value"]) < 5.0e-13, "stored single value")
    require(
        abs(float(abs(exact - single)) - stored["absolute_error"]) < 1.0e-15,
        "single independent error",
    )

    n_value, central_main = independent_central_main(mp.mpf("1e10"))
    require(n_value == 39894, "central Riemann-Siegel cutoff")
    grand = parse_multi(multi_path)
    fast_main = grand[7] - leading_correction(mp.mpf("1e10"))
    central_error = fast_main - central_main
    stored_error = payload["calibration"]["t1e10_et005"]["errors"][7]
    require(abs(central_error - stored_error) < 5.0e-10, "central direct-main error")

    require(len(parse_multi(h002_path)) == 15, "0.02 control values")
    require(len(parse_multi(h004_path)) == 15, "0.04 control values")
    require(len(parse_multi(t12_path)) == 15, "T=10^12 values")
    calibration = payload["calibration"]
    require(calibration["t1e10_et005"]["maximum_absolute"] < 0.01, "accepted grid")
    require(calibration["t1e12_et005"]["maximum_absolute"] < 0.005, "scaled grid")
    require(calibration["t1e10_h002"]["maximum_absolute"] > 1.0, "0.02 rejection")
    require(calibration["t1e10_h004"]["maximum_absolute"] > 3.0, "0.04 rejection")
    return 9


def polynomial(values: list[dict[str, str]]) -> np.ndarray:
    return np.asarray(
        [complex(float(value["re"]), float(value["im"])) for value in values]
    )


def independent_kernel_check(payload: dict) -> int:
    audited = payload["source_audit"]
    roots_payload = json.loads(
        (REPO_ROOT / audited["source:calibrated_roots"]["path"]).read_text(
            encoding="utf-8"
        )
    )
    adapter_payload = json.loads(
        (REPO_ROOT / audited["source:physical_adapter"]["path"]).read_text(
            encoding="utf-8"
        )
    )
    root = roots_payload["precision_ladder"][-1]["root"]
    physical = adapter_payload["physical_rows"][0]
    stored_root = payload["physical_kernels"][0]
    require(stored_root["terminal_peel"] == 0, "unexpected terminal peel")
    require(stored_root["N"] == int(physical["N"]), "stored physical N")

    t_value = mp.mpf(physical["chart"]["T_0"])
    n_value = int(physical["N"])
    heat_time = float(physical["chart"]["t"])
    sigma = parse_real_part(root["s_star"])
    omega = parse_complex(root["eta"])
    if n_value % 2:
        omega = -omega
    line = complex(mp.e ** (-1j * mp.siegeltheta(t_value)))
    require(abs((omega / line).imag) < 1.0e-14, "stationary endpoint phase lock")

    shift_values = np.asarray(stored_root["shifts"])
    require(np.array_equal(shift_values, np.arange(-7, 8) * 0.01), "shift grid")
    phases = np.asarray(
        [
            complex(
                mp.e
                ** (
                    -1j
                    * mp.siegeltheta(t_value + mp.mpf(str(float(delta))))
                )
            )
            for delta in shift_values
        ]
    )
    rng = np.random.default_rng(20260805)
    upper = math.log(n_value)
    points = np.unique(
        np.concatenate(([0.0, upper], upper * rng.random(3001) ** 2))
    )
    matrix = 2.0 * phases[None, :] * np.exp(
        1j * points[:, None] * shift_values[None, :]
    )
    weight = np.exp(
        heat_time * points**2 / 4.0 - (sigma - 0.5) * points
    )

    stored_by_key = {
        (row["family"], row["name"]): row for row in stored_root["rows"]
    }
    checks = (("base_rows", "P_V"), ("tangent_rows", "P_N"))
    for family, name in checks:
        stored = stored_by_key[(family, name)]
        coefficients = np.asarray(stored["real_hardy_coefficients"])
        require(len(coefficients) == 15, f"coefficient count {family}:{name}")
        source_polynomial = polynomial(physical["polynomials"][family][name])
        truth = (
            omega
            * weight
            * np.polynomial.polynomial.polyval(points, source_polynomial)
        )
        prediction = matrix @ coefficients
        scale = max(1.0, float(np.max(np.abs(truth))))
        error = float(np.max(np.abs(prediction - truth)) / scale)
        require(error < 1.0e-8, f"independent kernel error {family}:{name}")

    all_rows = [
        row for root_row in payload["physical_kernels"] for row in root_row["rows"]
    ]
    require(len(all_rows) == 24, "all physical kernel rows")
    require(all(len(row["real_hardy_coefficients"]) == 15 for row in all_rows), "all coefficient counts")
    maximum_l1 = max(row["coefficient_l1"] for row in all_rows)
    require(
        abs(maximum_l1 - payload["physical_kernel_summary"]["maximum_coefficient_l1"])
        < 1.0e-6,
        "coefficient summary",
    )
    require(maximum_l1 > 1.0e6, "conditioning guard")
    return 6


def check_rows(payload: dict) -> None:
    suffixes = (
        "pin",
        "port",
        "single",
        "multi",
        "scale",
        "phase",
        "kernel",
        "endpoint",
        "physical",
        "condition",
        "h002",
        "h004",
        "runtime",
        "boundary",
        "handoff",
    )
    expected = [
        f"shd_{index:02d}_{suffix}"
        for index, suffix in enumerate(suffixes, start=1)
    ]
    require([row["id"] for row in payload["rows"]] == expected, "row ids")
    require(sum(row["readiness"] == "open" for row in payload["rows"]) == 1, "open row")
    require(sum(row["readiness"] == "falsified" for row in payload["rows"]) == 2, "falsified controls")


def main() -> None:
    mp.mp.dps = 80
    require(RESULT_PATH.is_file(), "missing result")
    require(NOTE_PATH.is_file(), "missing note")
    payload = json.loads(RESULT_PATH.read_text(encoding="utf-8"))
    require(payload["kind"] == KIND, "kind")
    require(payload["schema_version"] == 1, "schema")
    require(payload["upstream"]["commit"] == "2e16dac3206b707052c3ac4cacdf3d1a2325e636", "upstream commit")
    require(payload["upstream"]["license"] == "GPL-3.0", "upstream license")
    source_audits = check_sources(payload)
    fixture_audits = check_fixture_values(payload)
    kernel_audits = independent_kernel_check(payload)
    check_rows(payload)

    summary = payload["summary"]
    require(summary["rows"] == 15, "summary rows")
    require(summary["accepted_shift_steps"] == [0.01], "accepted shifts")
    require(summary["rejected_shift_steps"] == [0.02, 0.04], "rejected shifts")
    require(summary["physical_kernel_rows"] == 24, "physical rows")
    require(summary["terminal_peel"] == 0, "terminal peel")
    require(summary["physical_values_evaluated"] == 0, "physical values")
    require(summary["interval_enclosures"] == 0, "intervals")
    require(summary["open_rows"] == 1, "open summary")

    note = NOTE_PATH.read_text(encoding="utf-8")
    require(payload["success"] in note, "note success")
    require("0.02" in note and "0.04" in note, "rejected controls in note")
    require("No terminal band peel" in note, "endpoint correction note")
    require("not a proof" in note.lower(), "proof boundary note")

    print(
        "validated Newman C1 shifted-Hardy diagnostic bridge gate: "
        f"15 rows, 0 issues, {source_audits} source/hash audits, "
        f"{fixture_audits} independent fixture audits, {kernel_audits} kernel audits, "
        "24 physical kernel fits, shift_step 0.01 retained, 0.02 and 0.04 rejected, "
        "0 physical values and 1 open resumable-evaluator row"
    )


if __name__ == "__main__":
    main()
