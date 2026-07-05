#!/usr/bin/env python3
"""Run the core proof-programme reproducibility gates.

This is a compact smoke/ledger runner for the RH dynamical-rigidity corpus.
It does not prove PF-infinity, Laguerre-Polya membership, RH, or Lambda <= 0.
It checks that the promoted finite-evidence manifests and executable
countermodel guards still match the advertised proof-programme status.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
import subprocess
import sys
from time import perf_counter


REPO_ROOT = Path(__file__).resolve().parents[3]


@dataclass(frozen=True)
class GateSpec:
    name: str
    command: tuple[str, ...]
    expected: tuple[str, ...]
    category: str
    slow: bool = False


GATES: tuple[GateSpec, ...] = (
    GateSpec(
        name="countermodel proof-safety gates",
        command=("work/rh_compute/scripts/countermodel_gate_examples.py",),
        expected=("validated 9 countermodel gate examples",),
        category="non-promotion guards",
    ),
    GateSpec(
        name="result-language boundary scan",
        command=("work/rh_compute/scripts/check_result_language_boundaries.py",),
        expected=("validated result-language boundaries: scanned", "0 overclaims"),
        category="non-promotion guards",
    ),
    GateSpec(
        name="output artifact status manifest",
        command=("work/rh_compute/scripts/check_output_status_manifest.py",),
        expected=("validated output artifact statuses: scanned", "0 status issues"),
        category="non-promotion guards",
    ),
    GateSpec(
        name="output reference integrity",
        command=("work/rh_compute/scripts/check_output_reference_integrity.py",),
        expected=("validated output references: scanned", "0 missing required paths"),
        category="non-promotion guards",
    ),
    GateSpec(
        name="proof-claim ledger",
        command=("work/rh_compute/scripts/check_proof_claim_ledger.py",),
        expected=("validated proof-claim ledger: 19 claims, 0 issues, 5 open theorem targets",),
        category="non-promotion guards",
    ),
    GateSpec(
        name="signed-Hankel finite certificates",
        command=("work/rh_compute/scripts/check_hankel_certificate_manifest.py",),
        expected=("validated 2500 signed-Hankel finite certificates",),
        category="promoted finite evidence",
    ),
    GateSpec(
        name="Toeplitz/PF finite certificates",
        command=("work/rh_compute/scripts/check_toeplitz_certificate_manifest.py",),
        expected=("validated 95 promoted positive certificate summaries and 1 negative control",),
        category="promoted finite evidence",
    ),
    GateSpec(
        name="Toeplitz/Jacobi-Trudi reindexing",
        command=("work/rh_compute/scripts/check_toeplitz_jacobi_trudi_map.py",),
        expected=("validated Toeplitz/Jacobi-Trudi reindexing: N=10, orders<=5, 124129 minors, 39094 nonzero maps, 85035 structural zeros",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Hankel sign-consistency reduction point audits",
        command=("work/rh_compute/scripts/check_hankel_sign_consistency_reduction_audit.py",),
        expected=("validated 20 reshaped Hankel sign-consistency point audits with 0 issues",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Hankel sign-consistency reduction finite certificates",
        command=("work/rh_compute/scripts/check_arb_hankel_sign_consistency_reduction_manifest.py",),
        expected=("validated 689795 Arb reshaped-Hankel sign-consistency finite certificates with 0 issues",),
        category="promoted finite evidence",
    ),
    GateSpec(
        name="shifted Hankel sign-consistency finite certificates",
        command=("work/rh_compute/scripts/check_arb_shifted_hankel_sign_consistency_manifest.py",),
        expected=("validated 818805 shifted Arb reshaped-Hankel sign-consistency finite certificates with 0 issues",),
        category="promoted finite evidence",
    ),
    GateSpec(
        name="Jensen/Hankel bridge algebra gate",
        command=("work/rh_compute/scripts/check_jensen_hankel_bridge_algebra.py",),
        expected=("validated Jensen/Hankel bridge algebra gate: degree2 identity and degree3 finite countermodel with 0 issues",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="signed-Hankel/Jensen bridge target specification",
        command=("work/rh_compute/scripts/check_signed_hankel_jensen_bridge_target.py",),
        expected=("validated signed-Hankel/Jensen bridge target specification with 0 issues",),
        category="open theorem target hygiene",
    ),
    GateSpec(
        name="Edrei-log sign diagnostics",
        command=("work/rh_compute/scripts/check_edrei_log_sign_manifest.py",),
        expected=("validated 320 finite Edrei-log sign diagnostics",),
        category="promoted finite evidence",
    ),
    GateSpec(
        name="Edrei power-Hankel diagnostics",
        command=("work/rh_compute/scripts/check_edrei_power_hankel_manifest.py",),
        expected=("validated 4205 finite Edrei power-Hankel diagnostics",),
        category="promoted finite evidence",
        slow=True,
    ),
    GateSpec(
        name="Edrei midpoint frontier non-promotion guard",
        command=("work/rh_compute/scripts/check_edrei_power_hankel_frontier_manifest.py",),
        expected=("validated 5 non-rigorous Edrei midpoint frontier scouts",),
        category="non-promotion guards",
    ),
    GateSpec(
        name="Edrei power-Hankel boundary repair manifest",
        command=("work/rh_compute/scripts/check_edrei_power_hankel_boundary_manifest.py",),
        expected=("validated 2 retired inconclusive blocker rows and 3 repaired positive boundary rows",),
        category="promoted finite evidence",
    ),
    GateSpec(
        name="Edrei moment-recurrence scout manifest",
        command=("work/rh_compute/scripts/check_edrei_quadrature_scout_manifest.py",),
        expected=("validated 1 positive Arb recurrence scout and 1 inconclusive frontier scout",),
        category="finite theorem-search diagnostics",
    ),
)


def command_for(spec: GateSpec) -> list[str]:
    return [sys.executable, *spec.command]


def tail(text: str, lines: int = 12) -> str:
    parts = text.splitlines()
    return "\n".join(parts[-lines:])


def run_gate(spec: GateSpec, timeout: int) -> dict:
    start = perf_counter()
    completed = subprocess.run(
        command_for(spec),
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
        timeout=timeout,
        check=False,
    )
    elapsed = perf_counter() - start
    combined = completed.stdout + "\n" + completed.stderr
    missing = [needle for needle in spec.expected if needle not in combined]
    ok = completed.returncode == 0 and not missing
    return {
        "name": spec.name,
        "category": spec.category,
        "command": " ".join(command_for(spec)),
        "returncode": completed.returncode,
        "elapsed_seconds": round(elapsed, 3),
        "ok": ok,
        "missing_expected": missing,
        "stdout_tail": tail(completed.stdout),
        "stderr_tail": tail(completed.stderr),
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--timeout", type=int, default=300, help="Per-gate timeout in seconds.")
    parser.add_argument("--skip-slow", action="store_true", help="Skip gates marked slow.")
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON summary.")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    specs = [spec for spec in GATES if not (args.skip_slow and spec.slow)]
    results = [run_gate(spec, args.timeout) for spec in specs]
    ok = all(result["ok"] for result in results)

    if args.json:
        print(json.dumps({"ok": ok, "gates": results}, indent=2, sort_keys=True))
    else:
        for result in results:
            status = "OK" if result["ok"] else "FAIL"
            print(
                f"{status} core gate: {result['name']} "
                f"({result['category']}, {result['elapsed_seconds']}s)"
            )
            if not result["ok"]:
                if result["missing_expected"]:
                    print(f"  missing expected: {result['missing_expected']}")
                if result["stdout_tail"]:
                    print("  stdout tail:")
                    print(result["stdout_tail"])
                if result["stderr_tail"]:
                    print("  stderr tail:")
                    print(result["stderr_tail"])
        print(f"validated {sum(1 for result in results if result['ok'])}/{len(results)} core proof-programme gates")

    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
