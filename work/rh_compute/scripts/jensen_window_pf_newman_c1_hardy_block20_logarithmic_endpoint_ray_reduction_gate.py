#!/usr/bin/env python3
"""Record the exact logarithmic endpoint-ray reduction and branch diagnostics."""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_ROOT))

import mpmath as mp

import jensen_window_pf_newman_c1_hardy_block20_native_q_required_compensation_gate as native_q
import jensen_window_pf_newman_c1_hardy_q_selector_cell_gate as cells


PAPER = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/Lewis_Brereton_2607.15310.pdf"
SECTOR = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_cubic_sector_contour_repair_gate.json"
)
EQ69 = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_exact_eq69_saddle_family_gate.json"
)
RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_logarithmic_endpoint_ray_reduction_gate.json"
)
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_block20_logarithmic_endpoint_ray_reduction_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = (
    REPO_ROOT
    / "work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_block20_logarithmic_endpoint_ray_reduction_gate.py"
)
WITNESS_CHAINS = (1, 2, 3, 36)
DIAGNOSTIC_DPS = 50


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def decimal(value: Fraction) -> str:
    return cells.decimal_text(value)


def harmonic(number: int) -> Fraction:
    return sum((Fraction(1, index) for index in range(1, number + 1)), Fraction(0))


def diagnostic(data: dict[str, Any], chain: int) -> dict[str, Any]:
    mp.mp.dps = DIAGNOSTIC_DPS
    parent = data["levels"][1]
    fractions = [cells.binary128_fraction(value) for value in parent["coefficients_hex"]]
    phi1, phi2, phi3 = [mp.mpf(value.numerator) / value.denominator for value in fractions]
    length = int(parent["length"])
    xi = phi1 + 2 * phi2 * length + 3 * phi3 * length**2
    transformed_length = int(mp.floor(xi))
    first_tail_index = transformed_length + 1
    first_saddle_index = 0 if phi1 < 0 else 1
    pi = mp.pi
    imaginary = 1j

    def phase(z: Any) -> Any:
        return phi1 * z + phi2 * z**2 + phi3 * z**3

    def derivative(z: Any) -> Any:
        return phi1 + 2 * phi2 * z + 3 * phi3 * z**2

    def kernel(u: Any, start: int) -> Any:
        exponential = mp.exp(2 * pi * imaginary * u)
        value = -mp.log1p(-exponential)
        for index in range(1, start):
            value -= exponential**index / index
        return value

    def ray(endpoint: int, angle: Any, phase_sign: int, start: int) -> Any:
        direction = mp.exp(imaginary * angle)

        def integrand(radius: Any) -> Any:
            u = radius * direction
            z = endpoint + u
            return (
                derivative(z)
                * mp.exp(phase_sign * 2 * pi * imaginary * phase(z))
                * kernel(u, start)
                * direction
            )

        return mp.quad(integrand, [0, mp.mpf("0.25"), 1, 4, 16, 64, 160, 320])

    b_angle = 5 * pi / 6 if phi3 < 0 else pi / 2
    c_angle = pi / 2 if phi3 < 0 else pi / 6
    b_source = ray(0, b_angle, -1, first_tail_index) - ray(
        length, b_angle, -1, first_tail_index
    )
    c_positive = -(
        ray(0, c_angle, 1, 1) - ray(length, c_angle, 1, 1)
    )
    c_source = mp.conj(c_positive)
    endpoint = mp.exp(-2 * pi * imaginary * phase(length))
    endpoint_half = (1 + endpoint) / 2
    a_endpoint = (
        imaginary
        * (endpoint - 1)
        * mp.mpf(harmonic(transformed_length).numerator)
        / (2 * pi * harmonic(transformed_length).denominator)
    )
    breakpoints = list(range(0, length + 1, 13))
    require(breakpoints[-1] == length, "log-ray diagnostic real breakpoints drift")
    zero_mode = (
        mp.quad(lambda y: mp.exp(-2 * pi * imaginary * phase(y)), breakpoints)
        if phi1 > 0
        else 0
    )
    formula = endpoint_half + a_endpoint + b_source + c_source + zero_mode
    parent_sum = mp.fsum(
        mp.exp(-2 * pi * imaginary * phase(index)) for index in range(length + 1)
    )
    saddle_sum = mp.fsum(
        mp.quad(
            lambda y, index=index: mp.exp(
                2 * pi * imaginary * (index * y - phase(y))
            ),
            breakpoints,
        )
        for index in range(first_saddle_index, transformed_length + 1)
    )
    target = parent_sum - saddle_sum
    gap = formula - target
    return {
        "chain": chain,
        "phi1_sign": "positive" if phi1 > 0 else "negative",
        "phi3_sign": "positive" if phi3 > 0 else "negative",
        "transformed_length": transformed_length,
        "first_saddle_index": first_saddle_index,
        "first_tail_index": first_tail_index,
        "b_angle": "5*pi/6" if phi3 < 0 else "pi/2",
        "c_angle": "pi/2" if phi3 < 0 else "pi/6",
        "formula_real": mp.nstr(mp.re(formula), 40),
        "formula_imag": mp.nstr(mp.im(formula), 40),
        "target_real": mp.nstr(mp.re(target), 40),
        "target_imag": mp.nstr(mp.im(target), 40),
        "absolute_gap": mp.nstr(abs(gap), 20),
        "rigorous": False,
    }


def build_note(artifact: dict[str, Any]) -> str:
    aggregate = artifact["aggregate"]
    witness_lines = "\n".join(
        f"chain {row['chain']:>3}: Phi1 {row['phi1_sign']}, Phi3 {row['phi3_sign']}, |gap|={row['absolute_gap']}"
        for row in artifact["diagnostic_witnesses"]
    )
    return f"""# Hardy block-20 logarithmic endpoint-ray reduction gate

Date: 2026-08-07

Status: exact analytic reduction with non-rigorous branch diagnostics; not a proof of W2--W4 error control or RH

## Summed ray kernel

For `Im(u)>0` and integer `m>=1`, define

```text
K_m(u)=sum_(n=m)^infinity exp(2*pi*i*n*u)/n
      =-Log(1-exp(2*pi*i*u))-sum_(n=1)^(m-1) exp(2*pi*i*n*u)/n.  (1)
```

The principal logarithm is analytic here because
`|exp(2*pi*i*u)|<1`.  Near the ray origin, `K_m(u)=-Log(-2*pi*i*u)+O(1)`,
so its logarithmic singularity is integrable.  This is the exact summed
version of the absolute-interchange statement in the sector-contour gate.

Put `E_-(z)=exp(-2*pi*i*F(z))`, `E_+(z)=exp(2*pi*i*F(z))`, and

```text
J^-_m(x;theta)=integral_0^infinity F'(x+r*e^(i*theta))
               E_-(x+r*e^(i*theta)) K_m(r*e^(i*theta)) e^(i*theta) dr,

J^+_m(x;theta)=the same expression with E_+ in place of E_-.
```

With `L=floor(xi)`, `m=L+1`, the exact source-oriented nonsaddle pieces are

```text
B^-=J^-_m(0;theta_b)-J^-_m(N;theta_b),
C^-=conj(J^+_1(N;theta_c)-J^+_1(0;theta_c)),           (2)
```

where `theta_b` and `theta_c` are the sign-aware angles from the preceding
gate.  Integer `N` is essential: it makes `K_m(N+u)=K_m(u)`.

## Exact recurrence remainder identity

Let

```text
P^-=sum_(k=0)^N E_-(k),
I69^-=sum_(n=j0)^L integral_0^N exp(2*pi*i*(n*y-F(y)))dy,
j0=ceil(Phi1) in {{0,1}},
H_L=sum_(n=1)^L 1/n.
```

Then direct integration by parts in `(67a)`, followed by (2), gives

```text
P^-=I69^-+Q_exact^-,

Q_exact^-=(1+E_-(N))/2 + i*(E_-(N)-1)*H_L/(2*pi)
          +B^-+C^-+1_(Phi1>0)*integral_0^N E_-(y)dy.  (3)
```

Thus the infinite nonsaddle families are reduced to four endpoint-ray
integrals, one optional finite zero-mode integral, and explicit endpoint
terms.  No truncation in `n` appears in (3).

## Roster and diagnostics

The exact selector audit covers {aggregate['call_count']} calls, with
`L=1 or 2`, `m=2 or 3`, and all four `(sign(Phi1),sign(Phi3))` combinations.
The following high-precision mpmath integrations test orientation only:

```text
{witness_lines}
```

These four values are non-rigorous diagnostics and are not used to certify
the exact theorem.  Their purpose is to catch conjugation, endpoint, and
zero-mode mistakes before implementing Arb quadrature.

## Pi provenance and boundary

Pi in (1)--(3) is inherited from the paper's Fourier exponential and fixes
both the logarithmic kernel and phase orientation.  It is not inserted as a
fitted normalization.

The reduction is exact under the already certified contour and interchange
hypotheses.  It does not yet rigorously enclose the logarithmic endpoint
integrals, compare them with W2--W4 on all 374 calls, establish a
height-uniform recurrence, control the outer Hardy remainder, or prove
`Lambda<=0`, PF-infinity, RH, or a prize-level conclusion.
"""


def main() -> int:
    for path in (PAPER, SECTOR, EQ69, CHECKER):
        require(path.is_file(), f"missing log-ray dependency: {path}")
    sector = json.loads(SECTOR.read_text(encoding="utf-8"))
    require(
        sector["status"] == "exact_sign_aware_cubic_nonsaddle_contour_existence_and_rotation_lemma",
        "log-ray sector gate not admitted",
    )
    recursive = native_q.load_recursive_chains()
    sector_rows = {int(row["chain"]): row for row in sector["rows"]}
    require(set(recursive) == set(sector_rows), "log-ray sector roster drift")

    counts: Counter[str] = Counter()
    rows: list[dict[str, Any]] = []
    for chain in sorted(recursive):
        data = recursive[chain]
        parent = data["levels"][1]
        child = data["levels"][2]
        phi1, phi2, phi3 = (
            cells.binary128_fraction(value) for value in parent["coefficients_hex"]
        )
        length = int(parent["length"])
        require(length == 104, "log-ray parent length drift")
        xi = phi1 + 2 * phi2 * length + 3 * phi3 * length**2
        transformed_length = xi.numerator // xi.denominator
        require(transformed_length == int(child["length"]) in (1, 2), "log-ray transformed length drift")
        first_tail_index = transformed_length + 1
        first_saddle_index = 0 if phi1 < 0 else 1
        require(first_saddle_index in (0, 1), "log-ray first saddle index drift")
        sign_key = f"phi1_{'positive' if phi1 > 0 else 'negative'}_phi3_{'positive' if phi3 > 0 else 'negative'}"
        counts[sign_key] += 1
        counts[f"L_{transformed_length}"] += 1
        sector_row = sector_rows[chain]
        rows.append(
            {
                "chain": chain,
                "sum_index": int(data["header"]["sum_index"]),
                "branch": int(data["header"]["branch"]),
                "phi1_sign": "positive" if phi1 > 0 else "negative",
                "phi3_sign": "positive" if phi3 > 0 else "negative",
                "transformed_length": transformed_length,
                "first_saddle_index": first_saddle_index,
                "first_tail_index": first_tail_index,
                "harmonic_exact": {
                    "numerator": harmonic(transformed_length).numerator,
                    "denominator": harmonic(transformed_length).denominator,
                    "decimal": decimal(harmonic(transformed_length)),
                },
                "b_angle": sector_row["b_contract"]["angle"],
                "c_angle": sector_row["c_contract"]["angle"],
            }
        )

    for key in (
        "phi1_positive_phi3_positive",
        "phi1_positive_phi3_negative",
        "phi1_negative_phi3_positive",
        "phi1_negative_phi3_negative",
    ):
        require(counts[key] > 0, f"log-ray missing sign branch: {key}")
    diagnostics = [diagnostic(recursive[chain], chain) for chain in WITNESS_CHAINS]
    diagnostic_signs = {
        (row["phi1_sign"], row["phi3_sign"]) for row in diagnostics
    }
    require(len(diagnostic_signs) == 4, "log-ray diagnostic sign coverage drift")
    for row in diagnostics:
        require(mp.mpf(row["absolute_gap"]) < mp.mpf("1e-12"), f"chain {row['chain']} log-ray diagnostic gap")

    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_block20_logarithmic_endpoint_ray_reduction_gate",
        "status": "exact_logarithmic_endpoint_ray_reduction_with_nonrigorous_orientation_diagnostics",
        "scope": {
            "block": 20,
            "recursive_call_count": 374,
            "diagnostic_witness_count": len(diagnostics),
            "diagnostic_decimal_digits": DIAGNOSTIC_DPS,
        },
        "kernel_identity": (
            "K_m(u)=sum_{n=m}^infinity exp(2*pi*i*n*u)/n="
            "-Log(1-exp(2*pi*i*u))-sum_{n=1}^{m-1}exp(2*pi*i*n*u)/n, Im(u)>0"
        ),
        "source_oriented_identity": {
            "b": "B^-=J^-_m(0;theta_b)-J^-_m(N;theta_b)",
            "c": "C^-=conj(J^+_1(N;theta_c)-J^+_1(0;theta_c))",
            "a_endpoint": "i*(E_-(N)-1)*H_L/(2*pi)",
            "zero_mode": "1_(Phi1>0)*integral_0^N E_-(y)dy",
            "complete": "P^-=I69^-+(1+E_-(N))/2+a_endpoint+B^-+C^-+zero_mode",
        },
        "analytic_contract": {
            "branch": "principal Log; 1-exp(2*pi*i*u) remains in the open right half-plane for Im(u)>0",
            "origin": "K_m(u)=-Log(-2*pi*i*u)+O(1), an integrable logarithmic singularity",
            "infinity": "the sign-aware cubic sector bound dominates the polynomial amplitude",
            "integer_endpoint": "N=104 implies exp(2*pi*i*n*(N+u))=exp(2*pi*i*n*u)",
        },
        "aggregate": {
            "call_count": len(rows),
            "length_one_count": counts["L_1"],
            "length_two_count": counts["L_2"],
            "sign_combination_counts": {
                key: counts[key]
                for key in sorted(counts)
                if key.startswith("phi1_")
            },
        },
        "rows": rows,
        "diagnostic_witnesses": diagnostics,
        "paper": {
            "path": relative(PAPER),
            "sha256": file_hash(PAPER),
            "pages": [23, 27, 28, 29, 30],
            "equations": [67, 68, 69, 82, 83, 88, 90, 92, 95, 96],
        },
        "sources": {
            "sector_gate": {"path": relative(SECTOR), "sha256": file_hash(SECTOR)},
            "equation_69_gate": {"path": relative(EQ69), "sha256": file_hash(EQ69)},
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "next_obligation": (
            "Implement singularity-subtracted Arb quadrature for the four logarithmic endpoint rays and optional zero "
            "mode, then compare the resulting exact nonsaddle balls with W2--W4 on all 374 calls."
        ),
        "proof_boundary": (
            "The logarithmic reduction is exact under the admitted contour and sum-interchange lemma; the four mpmath "
            "checks are non-rigorous orientation diagnostics. This does not rigorously enclose the endpoint integrals, "
            "prove W2--W4 error control or a height-uniform recurrence, control the outer Hardy representation, or prove "
            "Lambda<=0, PF-infinity, RH, or a prize-level conclusion."
        ),
    }
    RESULT.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    NOTE.write_text(build_note(artifact), encoding="utf-8")
    print("built logarithmic endpoint-ray reduction gate: 374 calls, 4 diagnostic branches")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
