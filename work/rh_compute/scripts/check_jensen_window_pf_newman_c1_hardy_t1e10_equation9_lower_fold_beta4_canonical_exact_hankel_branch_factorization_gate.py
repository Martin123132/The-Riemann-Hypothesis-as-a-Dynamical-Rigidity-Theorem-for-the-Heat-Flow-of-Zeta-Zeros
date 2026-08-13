#!/usr/bin/env python3
"""Independently check the corrected Airy/Hankel branch factorization."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))

import mpmath as mp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_beta4_canonical_exact_hankel_branch_factorization_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    require(RESULT.is_file() and NOTE.is_file(), "missing result or note")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact["kind"] == STEM and artifact["passed"] is True, "gate identity drift")
    for record in artifact["dependencies"].values():
        path = REPO_ROOT / record["path"]
        require(path.is_file() and file_hash(path) == record["sha256"], f"dependency hash drift: {path}")
    for record in artifact["sources"].values():
        path = REPO_ROOT / record["path"]
        require(path.is_file() and file_hash(path) == record["sha256"], f"source hash drift: {path}")

    mp.mp.dps = 78
    for x in (mp.mpf("0.00007"), mp.mpf("0.031"), mp.mpf("0.9"), mp.mpf("8.5"), mp.mpf("93")):
        xi = mp.mpf(2) * x**mp.mpf("1.5") / 3
        a_bessel = mp.sqrt(x) / 3 * (mp.besselj(-mp.mpf(1) / 3, xi) + mp.besselj(mp.mpf(1) / 3, xi))
        ax_bessel = x / 3 * (mp.besselj(-mp.mpf(2) / 3, xi) - mp.besselj(mp.mpf(2) / 3, xi))
        require(abs(a_bessel - mp.airyai(-x)) < mp.mpf("1e-68"), f"Bessel Ai identity drift at X={x}")
        require(abs(ax_bessel + mp.airyai(-x, 1)) < mp.mpf("1e-68"), f"Bessel A_X identity drift at X={x}")

        h1_a = mp.sqrt(x) / (2 * mp.sqrt(3)) * mp.exp(1j * mp.pi / 6) * mp.hankel1(mp.mpf(1) / 3, xi)
        h2_a = mp.sqrt(x) / (2 * mp.sqrt(3)) * mp.exp(-1j * mp.pi / 6) * mp.hankel2(mp.mpf(1) / 3, xi)
        h1_x = -x / (2 * mp.sqrt(3)) * mp.exp(-1j * mp.pi / 6) * mp.hankel1(mp.mpf(2) / 3, xi)
        h2_x = -x / (2 * mp.sqrt(3)) * mp.exp(1j * mp.pi / 6) * mp.hankel2(mp.mpf(2) / 3, xi)
        require(abs(h1_a + h2_a - a_bessel) < mp.mpf("1e-68"), "Hankel Ai split drift")
        require(abs(h1_x + h2_x - ax_bessel) < mp.mpf("1e-68"), "Hankel derivative split drift")
        require(abs(h2_a - mp.conj(h1_a)) < mp.mpf("1e-68"), "Ai branch conjugation drift")
        require(abs(h2_x - mp.conj(h1_x)) < mp.mpf("1e-68"), "A_X branch conjugation drift")

    for nu, xi in ((mp.mpf(1) / 3, mp.mpf("0.4")), (mp.mpf(2) / 3, mp.mpf("7.25"))):
        h1 = lambda value: mp.hankel1(nu, value)
        h2 = lambda value: mp.hankel2(nu, value)
        wronskian = h1(xi) * mp.diff(h2, xi) - mp.diff(h1, xi) * h2(xi)
        require(abs(wronskian + 4j / (mp.pi * xi)) < mp.mpf("1e-68"), "Hankel Wronskian drift")

    decision = artifact["decision"]
    require(decision["exact_Hankel_branch_factorization_proved"] is True, "exact split flag lost")
    require(decision["raw_Hankel_functions_singular_at_X_zero"] is True, "raw-Hankel zero guard lost")
    require(decision["prefactored_branches_have_finite_positive_limits"] is True, "finite branch-limit guard lost")
    require(decision["Airy_or_explicit_limit_required_at_X_zero"] is True, "fold-collar interpretation guard lost")
    require(decision["ordinary_logistic_branch_identification_proved"] is False, "ordinary bridge overclaim")
    require(decision["corridor_remainder_bound_proved"] is False, "corridor overclaim")
    require(decision["rh_implication"] is False, "RH overclaim")
    note = NOTE.read_text(encoding="utf-8")
    for token in ("exact Bessel connection", "not fitted phases", "finite one-sided branch limits", "no RH"):
        require(token in note, f"note boundary token missing: {token}")
    print("independently checked exact Airy/Hankel split: 5 positive-X samples, 2 Wronskians, 0 corridor claims", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
