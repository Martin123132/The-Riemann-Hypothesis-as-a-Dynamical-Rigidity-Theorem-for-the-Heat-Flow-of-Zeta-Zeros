#!/usr/bin/env python3
"""Certify a long lambda-zero order-four prefix from an Arb Xi series."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
import math
import os
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402


MAX_N = 500
MAX_K = MAX_N + 6
PRECISION_BITS = 24_576
SERIES_ORDER = 2 * MAX_K + 4
SERIALIZED_DIGITS = 180
CHECK_PRECISION_BITS = 2048
DEFAULT_CACHE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "arb_xi_lambda0_order4_prefix_coefficients_n0_n506_bits24576.jsonl"
)
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "arb_xi_lambda0_order4_prefix_certificate.json"
)
DEFAULT_NOTE = REPO_ROOT / "outputs/arb_xi_lambda0_order4_prefix_certificate.md"


@dataclass(frozen=True)
class PrefixRow:
    n: int
    H4_ball: str
    F_ball: str
    F_over_G1_squared_ball: str


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def compact_ball(value: flint.arb, digits: int = SERIALIZED_DIGITS) -> str:
    text = value.str(digits)
    rounded = flint.arb(text)
    if not rounded.contains(value):
        raise RuntimeError("serialized Arb ball does not contain its source")
    return text


def xi_series_coefficients() -> list[flint.arb]:
    flint.ctx.prec = PRECISION_BITS
    flint.ctx.cap = SERIES_ORDER
    flint.ctx.threads = min(8, os.cpu_count() or 1)
    s = flint.acb_series([flint.acb("0.5"), 1], prec=SERIES_ORDER)
    xi = (
        (s * (s - 1) / 2)
        * ((-s / 2) * flint.acb.pi().log()).exp()
        * (s / 2).gamma()
        * s.zeta()
    )
    coefficients = []
    for k in range(MAX_K + 1):
        even = xi[2 * k]
        if not even.imag.contains(0):
            raise RuntimeError(f"Xi coefficient has nonzero imaginary part at k={k}")
        coefficient = even.real * math.factorial(k) / (flint.arb(4) ** (k + 1))
        if not coefficient > 0:
            raise RuntimeError(f"A_{k}(0) is not sign-separated: {coefficient}")
        coefficients.append(coefficient)
    return coefficients


def write_cache(path: Path, coefficients: list[flint.arb]) -> None:
    metadata = {
        "kind": "arb_xi_lambda0_order4_prefix_metadata",
        "method": "direct Acb power series for xi(1/2+z)",
        "source_formula": (
            "xi(s)=s*(s-1)*pi^(-s/2)*Gamma(s/2)*zeta(s)/2; "
            "A_k(0)=k!*[z^(2k)]xi(1/2+z)/4^(k+1)"
        ),
        "precision_bits": PRECISION_BITS,
        "series_order": SERIES_ORDER,
        "max_n": MAX_N,
        "max_k": MAX_K,
        "serialized_digits": SERIALIZED_DIGITS,
        "flint_version": flint.__version__,
    }
    lines = [json.dumps(metadata, sort_keys=True)]
    for k, coefficient in enumerate(coefficients):
        lines.append(
            json.dumps(
                {
                    "kind": "arb_xi_lambda0_coefficient",
                    "k": k,
                    "A_ball": compact_ball(coefficient),
                },
                sort_keys=True,
            )
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def load_cache(path: Path) -> tuple[dict, dict[int, flint.arb]]:
    metadata: dict | None = None
    coefficients: dict[int, flint.arb] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            row = json.loads(line)
            if row.get("kind") == "arb_xi_lambda0_order4_prefix_metadata":
                metadata = row
            elif row.get("kind") == "arb_xi_lambda0_coefficient":
                coefficients[int(row["k"])] = flint.arb(row["A_ball"])
    if metadata is None:
        raise RuntimeError("Xi prefix cache metadata is missing")
    missing = [k for k in range(MAX_K + 1) if k not in coefficients]
    if missing:
        raise RuntimeError(f"Xi prefix cache is missing coefficients: {missing[:10]}")
    return metadata, coefficients


def stable_margin(coefficients: dict[int, flint.arb], n: int) -> tuple[flint.arb, flint.arb]:
    ratios = [coefficients[k + 1] / coefficients[k] for k in range(n, n + 6)]
    contractions = [ratios[k + 1] / ratios[k] for k in range(5)]
    defects = [1 - value for value in contractions]
    gaps = [
        defects[center] ** 2
        - contractions[center] ** 2
        * defects[center - 1]
        * defects[center + 1]
        for center in range(1, 4)
    ]
    margin = gaps[1] ** 2 - contractions[2] ** 3 * gaps[0] * gaps[2]
    return margin, margin / gaps[1] ** 2


def build_artifact_from_cache(cache_path: Path = DEFAULT_CACHE) -> dict:
    flint.ctx.prec = CHECK_PRECISION_BITS
    metadata, coefficients = load_cache(cache_path)
    expected_metadata = {
        "precision_bits": PRECISION_BITS,
        "series_order": SERIES_ORDER,
        "max_n": MAX_N,
        "max_k": MAX_K,
        "serialized_digits": SERIALIZED_DIGITS,
    }
    for key, expected in expected_metadata.items():
        if metadata.get(key) != expected:
            raise RuntimeError(f"Xi prefix metadata mismatch for {key}")
    rows: list[PrefixRow] = []
    for n in range(MAX_N + 1):
        determinant = flint.arb_mat(
            [[coefficients[n + i + j] for j in range(4)] for i in range(4)]
        ).det()
        margin, normalized = stable_margin(coefficients, n)
        if not determinant > 0:
            raise RuntimeError(f"H4 sign not separated at n={n}: {determinant}")
        if not margin > 0:
            raise RuntimeError(f"stable H4 margin not separated at n={n}: {margin}")
        rows.append(
            PrefixRow(
                n=n,
                H4_ball=compact_ball(determinant, 100),
                F_ball=compact_ball(margin, 100),
                F_over_G1_squared_ball=compact_ball(normalized, 100),
            )
        )
    minimum_margin = min(rows, key=lambda row: float(flint.arb(row.F_ball).mid()))
    minimum_normalized = min(
        rows,
        key=lambda row: float(flint.arb(row.F_over_G1_squared_ball).mid()),
    )
    relative_radius_rows = []
    for k in (0, 64, 128, 256, 384, MAX_K):
        value = coefficients[k]
        relative_radius_rows.append(
            {
                "k": k,
                "A_ball": compact_ball(value, 100),
                "relative_radius_upper": compact_ball(value.rad() / value.abs_lower(), 30),
            }
        )
    return {
        "kind": "arb_xi_lambda0_order4_prefix_certificate",
        "date": "2026-07-13",
        "status": "rigorous direct-Xi-series lambda-zero contiguous order-four prefix through n=500",
        "proof_boundary": (
            "This is a finite Arb certificate for 0<=n<=500 only. It is not an "
            "all-shift order-four theorem, PF-infinity, RH, or Lambda<=0."
        ),
        "method": metadata,
        "cache": {
            "path": str(cache_path.relative_to(REPO_ROOT)).replace("\\", "/"),
            "sha256": sha256(cache_path),
            "coefficient_rows": len(coefficients),
        },
        "finite": {
            "lambda": "0",
            "n_range": [0, MAX_N],
            "coefficient_range": [0, MAX_K],
            "all_coefficients_positive": True,
            "all_H4_positive": True,
            "all_stable_margins_positive": True,
            "minimum_margin_n": minimum_margin.n,
            "minimum_margin_ball": minimum_margin.F_ball,
            "minimum_normalized_n": minimum_normalized.n,
            "minimum_normalized_ball": minimum_normalized.F_over_G1_squared_ball,
            "representative_coefficient_precision": relative_radius_rows,
            "rows": [asdict(row) for row in rows],
        },
        "summary": {
            "coefficient_rows": MAX_K + 1,
            "prefix_rows": MAX_N + 1,
            "positive_coefficients": MAX_K + 1,
            "positive_H4_rows": MAX_N + 1,
            "positive_stable_margin_rows": MAX_N + 1,
            "inconclusive_rows": 0,
        },
        "generator": "work/rh_compute/scripts/arb_xi_lambda0_order4_prefix_certificate.py",
        "checker": "work/rh_compute/scripts/check_arb_xi_lambda0_order4_prefix_certificate.py",
    }


def write_note(path: Path, artifact: dict) -> None:
    finite = artifact["finite"]
    method = artifact["method"]
    lines = [
        "# Arb Xi Lambda-Zero Order-Four Prefix Certificate",
        "",
        "Date: 2026-07-13",
        "",
        "Status: rigorous direct-Xi-series finite certificate through `n=500`.",
        "This is not an all-shift order-four theorem, PF-infinity, RH, or",
        "`Lambda <= 0`.",
        "",
        "```text",
        "work/rh_compute/results/arb_xi_lambda0_order4_prefix_certificate.json",
        "work/rh_compute/results/arb_xi_lambda0_order4_prefix_coefficients_n0_n506_bits24576.jsonl",
        "python work/rh_compute/scripts/arb_xi_lambda0_order4_prefix_certificate.py",
        "python work/rh_compute/scripts/check_arb_xi_lambda0_order4_prefix_certificate.py",
        "```",
        "",
        "## Direct Series",
        "",
        "Python-flint/Arb expands the exact identity",
        "",
        "```text",
        "xi(s)=s*(s-1)*pi^(-s/2)*Gamma(s/2)*zeta(s)/2",
        "A_k(0)=k!*[z^(2k)]xi(1/2+z)/4^(k+1).",
        "```",
        "",
        f"The expansion uses `{method['precision_bits']}` bits and series order",
        f"`{method['series_order']}`. Each serialized decimal Arb ball is rounded",
        "outward and checked to contain the original high-precision ball.",
        "",
        "## Certified Prefix",
        "",
        "Recomputing the raw `4x4` Hankel determinant and the stable gap",
        "factorization from the cached coefficient balls proves",
        "",
        "```text",
        "A_k(0)>0 for every 0<=k<=506",
        "H_(4,n)(0)>0 for every 0<=n<=500",
        "F_n(0)>0 for every 0<=n<=500",
        "```",
        "",
        f"The smallest recorded stable margin occurs at `n={finite['minimum_margin_n']}`:",
        "",
        "```text",
        f"F_500={finite['minimum_margin_ball']}",
        f"F_500/G_501^2={finite['minimum_normalized_ball']}",
        "```",
        "",
        "This extends the previous rigorous endpoint prefix from `n=58` to",
        "`n=500`. It remains finite evidence; it composes with the separate",
        "eventual Xi-asymptotic theorem only after that theorem's threshold is",
        "made explicit.",
        "",
        "```text",
        "outputs/jensen_window_pf_compound_order4_lambda0_eventual_positivity_certificate.md",
        "outputs/formal_core.md",
        "```",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    coefficients = xi_series_coefficients()
    write_cache(args.cache, coefficients)
    artifact = build_artifact_from_cache(args.cache)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    write_note(args.note, artifact)
    print(
        "built Arb Xi lambda-zero order-four prefix: "
        f"{artifact['summary']['coefficient_rows']} coefficients, "
        f"{artifact['summary']['positive_H4_rows']} positive H4 rows, "
        "0 inconclusive"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
