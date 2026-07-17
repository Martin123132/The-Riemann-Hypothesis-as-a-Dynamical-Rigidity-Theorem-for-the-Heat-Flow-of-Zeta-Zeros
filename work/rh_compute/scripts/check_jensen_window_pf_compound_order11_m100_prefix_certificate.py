#!/usr/bin/env python3
"""Validate the lambda=-100 signed order-eleven prefix certificate."""

from __future__ import annotations

from decimal import Decimal
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
for candidate in (SCRIPT_DIR, VENDOR):
    if candidate.exists() and str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

import flint  # noqa: E402

import jensen_window_pf_compound_order11_m100_prefix_certificate as generator  # noqa: E402
import jensen_window_pf_compound_order9_m100_prefix_certificate as prefix  # noqa: E402
import jensen_window_pf_endpoint_order10_counterexample as endpoint  # noqa: E402


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def validate_repair_overlap() -> list[str]:
    issues = []
    old = prefix.load_source(generator.SOURCE_PATHS[0])
    repair = prefix.load_source(generator.PRECISION_REPAIR)
    if set(repair) != set(range(134, 153)):
        issues.append("precision repair index range changed")
        return issues
    for index, repaired in repair.items():
        if index not in old:
            issues.append(f"old source misses repair index {index}")
        elif not bool(repaired.overlaps(old[index])):
            issues.append(f"precision repair is disjoint at k={index}")
    extension = prefix.load_source(generator.ENDPOINT_EXTENSION)
    if set(extension) != {1261, 1262}:
        issues.append("endpoint extension index range changed")
    return issues


def validate_direct_rows(
    artifact: dict,
    values: dict[int, flint.arb],
) -> list[str]:
    issues = []
    rows = artifact.get("direct_audit", [])
    if [row.get("n") for row in rows] != list(range(5)):
        return ["direct audit rows are not n=0..4"]
    for row in rows:
        n = int(row["n"])
        q9 = endpoint.direct_signed_hankel(values, 9, n + 2)
        q10_left = endpoint.direct_signed_hankel(values, 10, n)
        q10_center = endpoint.direct_signed_hankel(values, 10, n + 1)
        q10_right = endpoint.direct_signed_hankel(values, 10, n + 2)
        q11 = endpoint.direct_signed_hankel(values, 11, n)
        schur = endpoint.direct_deep_schur(values, 11, n + 10)
        numerator = q10_center**2 - q10_left * q10_right
        checks = {
            "q9_positive": endpoint.sign_class(q9) == "positive",
            "numerator_positive": endpoint.sign_class(numerator) == "positive",
            "q11_positive": endpoint.sign_class(q11) == "positive",
            "deep_schur_positive": endpoint.sign_class(schur) == "positive",
            "toda_residual_contains_zero": bool((q11 * q9 - numerator).contains(0)),
            "schur_residual_contains_zero": bool(
                (schur - q11 / values[0] ** 11).contains(0)
            ),
        }
        for name, passed in checks.items():
            if not passed:
                issues.append(f"direct audit {name} failed at n={n}")
        if not q11.overlaps(flint.arb(row["Q11_ball"])):
            issues.append(f"stored direct Q11 ball is disjoint at n={n}")
        if not numerator.overlaps(flint.arb(row["condensation_numerator_ball"])):
            issues.append(f"stored direct numerator is disjoint at n={n}")
    return issues


def main() -> int:
    path = generator.DEFAULT_OUT
    artifact = load_json(path)
    issues = []
    if artifact.get("kind") != (
        "jensen_window_pf_compound_order11_m100_prefix_certificate"
    ):
        issues.append("bad artifact kind")
    if artifact.get("status") != (
        "rigorous lambda=-100 signed order-eleven prefix through n=1242"
    ):
        issues.append("bad artifact status")
    expected_exact = {
        "orientation": "epsilon_11=(-1)^55=-1, Q_(11,n)=-H_(11,n)",
        "condensation": (
            "Q_(11,n)*Q_(9,n+2)=Q_(10,n+1)^2-"
            "Q_(10,n)*Q_(10,n+2)"
        ),
        "relative_coordinate": (
            "M_n=Q_(10,n+1)^2/(Q_(10,n)*Q_(10,n+2))-1"
        ),
        "exceptional_prefix": (
            "Q_(10,n+1)^2-Q_(10,n)*Q_(10,n+2)>0 for 0<=n<=3"
        ),
        "positive_cone_block": "M_n>0 for every 4<=n<=1242",
        "finite_prefix": "Q_(11,n)(-100)>0 for every 0<=n<=1242",
        "analytic_handoff": "prove Q_(11,n)(-100)>0 for every n>=1243",
    }
    if artifact.get("exact") != expected_exact:
        issues.append("exact identity block changed")

    sources = artifact.get("sources", {})
    for name in ("order10_endpoint", "precision_repair", "endpoint_extension"):
        source = sources.get(name, {})
        source_path = REPO_ROOT / source.get("path", "missing")
        if not source_path.exists():
            issues.append(f"missing source {name}")
        elif source.get("sha256") != generator.sha256(source_path):
            issues.append(f"source hash changed for {name}")
        summary_path_text = source.get("summary_path")
        if summary_path_text:
            summary_path = REPO_ROOT / summary_path_text
            if not summary_path.exists():
                issues.append(f"missing summary {name}")
            elif source.get("summary_sha256") != generator.sha256(summary_path):
                issues.append(f"summary hash changed for {name}")

    issues.extend(validate_repair_overlap())
    flint.ctx.prec = generator.PRECISION_BITS
    values, diagnostics = generator.load_coefficients()
    if len(values) != 1263 or set(values) != set(range(1263)):
        issues.append("coefficient coverage is not A_0..A_1262")
    if artifact.get("source_diagnostics") != diagnostics:
        issues.append("source diagnostics changed")

    rebuilt = generator.stable_order11_prefix(values)
    rebuilt.pop("internal")
    finite = artifact.get("finite", {})
    if finite != rebuilt:
        issues.append("full stable prefix rebuild differs from artifact")
    finite_rows = finite.get("rows", [])
    if [row.get("n") for row in finite_rows] != list(range(1243)):
        issues.append("finite rows are not contiguous n=0..1242")
    else:
        for row in finite_rows:
            n = int(row["n"])
            try:
                numerator_lower = Decimal(row["condensation_numerator_lower"])
                q11_lower = Decimal(row["Q11_lower"])
            except Exception as exc:
                issues.append(f"unparseable lower bound at n={n}: {exc}")
                continue
            if numerator_lower <= 0 or q11_lower <= 0:
                issues.append(f"nonpositive stored Q11 bound at n={n}")
            if row.get("Q11_sign") != "positive_by_signed_condensation":
                issues.append(f"bad Q11 sign label at n={n}")
            if n < 4:
                if row.get("relative_Q10_margin_lower") is not None:
                    issues.append(f"exceptional row has relative margin at n={n}")
            else:
                try:
                    relative_lower = Decimal(row["relative_Q10_margin_lower"])
                except Exception as exc:
                    issues.append(f"unparseable relative margin at n={n}: {exc}")
                else:
                    if relative_lower <= 0:
                        issues.append(f"nonpositive relative margin at n={n}")

    q10_map = finite.get("Q10_sign_map", {})
    if q10_map != {
        "negative_indices": [0, 1, 2, 3],
        "positive_range": [4, 1244],
        "inconclusive_indices": [],
    }:
        issues.append("bad Q10 sign map")
    if finite.get("minimum_relative_n") != 1242:
        issues.append("minimum relative margin is not at n=1242")
    elif Decimal(finite.get("minimum_relative_lower", "0")) <= Decimal("0.005"):
        issues.append("minimum relative margin lower bound is not above 0.005")

    issues.extend(validate_direct_rows(artifact, values))
    summary = artifact.get("summary", {})
    expected_summary = {
        "rows": 7,
        "ready_rows": 6,
        "open_rows": 1,
        "coefficients": 1263,
        "precision_repair_rows": 19,
        "endpoint_extension_rows": 2,
        "Q11_numerators": 1243,
        "positive_Q11_rows": 1243,
        "negative_Q11_rows": 0,
        "inconclusive_Q11_rows": 0,
        "exceptional_Q11_rows": 4,
        "positive_cone_Q11_rows": 1239,
        "direct_Q11_checks": 5,
        "finite_prefix_theorems": 1,
        "open_analytic_tail_targets": 1,
    }
    if summary != expected_summary:
        issues.append("summary counts changed")

    if issues:
        for issue in issues:
            print(f"issue: {issue}")
        return 1
    print(
        "validated order-eleven lambda=-100 prefix: "
        "1263 coefficients, 1243 positive Q11 rows, "
        "5 direct audits, 0 issues"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
