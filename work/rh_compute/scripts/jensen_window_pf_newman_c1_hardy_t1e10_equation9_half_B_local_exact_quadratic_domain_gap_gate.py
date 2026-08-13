#!/usr/bin/env python3
"""Certify the exact Gaussian/Fresnel domain gap for B modes 621 and 622."""

from __future__ import annotations

import hashlib
import json
import os
import sys
import time
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import arb, ctx
import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_local_exact_quadratic_domain_gap_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "pair_Morse": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_signed_pair_logistic_morse_transform_gate.json",
    "discrete_gap": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_face_discrete_joint_gradient_gap_gate.json",
}

PRECISION = 100
T = 10_000_000_000
B = 5_122_421


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing dependency: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def set_low_priority() -> str:
    try:
        import psutil

        process = psutil.Process()
        if os.name == "nt":
            process.nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
            return "below_normal"
        process.nice(10)
        return "nice_10"
    except Exception as exc:  # pragma: no cover
        return f"unavailable:{type(exc).__name__}"


def symbolic_certificate() -> dict[str, str]:
    y, u = sp.symbols("y u", real=True)
    phase = -y**2 + sp.pi * u**2 / 2
    gradient = [sp.diff(phase, y), sp.diff(phase, u)]
    require(gradient == [-2 * y, sp.pi * u], "quadratic gradient failed")
    return {
        "outer_coordinate": "y=sqrt(t)*s/2, s=sgn(v-1)sqrt(2[v-1-log v])",
        "outer_inverse": "v=-W_-1(-exp(-1-s^2/2)) for s>=0 and v=-W_0(-exp(-1-s^2/2)) for s<=0",
        "normal_coordinate": "u=q_B(m,x)=sqrt(x/2)(B-2m/x)",
        "exact_phase": "Psi_m-Psi_m(x_m)+normal completion=-y^2+pi*u^2/2",
        "phase_gradient": "(-2y,pi*u)",
        "mode_621_domain": "u<=q_621(y); q_621(0)<0, so the Gaussian origin is excluded",
        "mode_622_domain": "u>=q_622(y); q_622(0)>0, so the Gaussian origin is excluded",
        "minimax_rule": "On each nearest boundary arc, 2|y| and pi|q_m(y)| vary oppositely; the minimum max component is their unique equality point.",
    }


def interval_certificate() -> dict[str, Any]:
    ctx.dps = PRECISION
    ctx.threads = 1
    t = arb(T)
    endpoint = arb(B)
    pi = arb.pi()

    def coordinates(mode: int, x: arb) -> tuple[arb, arb]:
        m = arb(mode)
        r = t / (2 * pi * m**2)
        v = (1 - x) / (r * x)
        defect = 2 * (v - 1 - v.log())
        require(defect > 0, "outer Morse defect is not positive")
        s = defect.sqrt() if v > 1 else -defect.sqrt()
        y = t.sqrt() * s / 2
        q = (endpoint * x - 2 * m) / (2 * x).sqrt()
        return y, q

    def bracket(mode: int, low_text: str, high_text: str, side: str) -> dict[str, Any]:
        low, high = arb(low_text), arb(high_text)
        y_low, q_low = coordinates(mode, low)
        y_high, q_high = coordinates(mode, high)
        if side == "left":
            f_low = -2 * y_low + pi * q_low
            f_high = -2 * y_high + pi * q_high
            require(f_low < 0 and f_high > 0, "left exact balance bracket failed")
            margin = min(-2 * y_low, -2 * y_high)
            q_sign = q_high
            require(y_low < 0 and y_high < 0 and q_sign < 0, "left-domain signs failed")
        else:
            f_low = 2 * y_low - pi * q_low
            f_high = 2 * y_high - pi * q_high
            require(f_low > 0 and f_high < 0, "right exact balance bracket failed")
            margin = min(2 * y_low, 2 * y_high)
            q_sign = q_low
            require(y_low > 0 and y_high > 0 and q_sign > 0, "right-domain signs failed")
        y_mid, q_mid = coordinates(mode, (low + high) / 2)
        euclidean_gradient = (4 * y_mid**2 + pi**2 * q_mid**2).sqrt()
        return {
            "mode": mode,
            "x_bracket": [low_text, high_text],
            "balance_sign_balls": [f_low.str(PRECISION, more=True), f_high.str(PRECISION, more=True)],
            "y_mid_ball": y_mid.str(PRECISION, more=True),
            "q_mid_ball": q_mid.str(PRECISION, more=True),
            "max_gradient_component_lower_ball": margin.str(PRECISION, more=True),
            "euclidean_gradient_mid_ball": euclidean_gradient.str(PRECISION, more=True),
        }

    mode_621 = bracket(621, "0.00024238523658", "0.00024238523660", "left")
    mode_622 = bracket(622, "0.00024291643709", "0.00024291643710", "right")
    require(arb(mode_621["max_gradient_component_lower_ball"]) > arb("57.17"), "mode-621 exact gap below 57.17")
    require(arb(mode_622["max_gradient_component_lower_ball"]) > arb("45.63"), "mode-622 exact gap below 45.63")

    def saddle_q(mode: int) -> arb:
        m = arb(mode)
        x = 2 * pi * m**2 / (t + 2 * pi * m**2)
        return (endpoint * x - 2 * m) / (2 * x).sqrt()

    q621 = saddle_q(621)
    q622 = saddle_q(622)
    require(q621 < -arb("50.44") and q622 > arb("40.28"), "saddle normal margins drift")

    return {
        "height": T,
        "endpoint": B,
        "mode_621": mode_621,
        "mode_622": mode_622,
        "q_621_at_outer_saddle_ball": q621.str(PRECISION, more=True),
        "q_622_at_outer_saddle_ball": q622.str(PRECISION, more=True),
        "local_exact_max_gradient_component_lower_bound": "45.63",
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["interval_certificate"]
    a, b = c["mode_621"], c["mode_622"]
    return f"""# Exact quadratic-domain gap for the local B modes

Date: 2026-08-13

Status: exact-coordinate interval certificate; not a proof of the complete
B-face estimate

Use the global outer Morse coordinate and the exact endpoint Fresnel
coordinate

```text
y=sqrt(t)s/2,
u=q_B(m,x)=sqrt(x/2)(B-2m/x).                         (QD1)
```

The complete two-dimensional phase defect is exactly

```text
-y^2+pi*u^2/2,   grad=(-2y,pi*u).                    (QD2)
```

There is no phase Taylor remainder.  For mode 621 the residual domain is on
the side `u<=q_621(y)` and `q_621(0)<0`; for mode 622 it is on the side
`u>=q_622(y)` and `q_622(0)>0`.  Thus both exact domains exclude the Gaussian
saddle `(0,0)`.

At the outer saddles,

```text
q_621(0)={c['q_621_at_outer_saddle_ball']},
q_622(0)={c['q_622_at_outer_saddle_ball']}.            (QD3)
```

On the nearest boundary arc the two gradient components vary oppositely, so
the minimax point is the unique balance `2|y|=pi|q_m(y)|`.  Arb brackets the
mode-621 balance by

```text
x in [{a['x_bracket'][0]}, {a['x_bracket'][1]}],
y={a['y_mid_ball']}, q={a['q_mid_ball']},
max(2|y|,pi|q|)>{a['max_gradient_component_lower_ball']}>57.17. (QD4)
```

For mode 622,

```text
x in [{b['x_bracket'][0]}, {b['x_bracket'][1]}],
y={b['y_mid_ball']}, q={b['q_mid_ball']},
max(2|y|,pi|q|)>{b['max_gradient_component_lower_ball']}>45.63. (QD5)
```

These are the exact nonstationary margins for the only two B-local integer
modes.  They are stronger than the common face-coordinate bound because
they use the true quadratic phase scales.  The tangent half-plane model can
be nearly degenerate only after replacing the curved boundary by a remote
line; neither exact local domain contains a joint saddle.

Pi provenance: `pi` comes from the equation-(9) outer and endpoint phases.
No fitted constant is used.

Proof boundary: exact phase coordinates and local-domain gradient margins
only.  No transformed-amplitude derivative bound, integration-by-parts
constant, complete B estimate, A-fold splice, complete paired residual,
`T_upper`, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion is
established.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")
    require(
        dependencies["pair_Morse"]["decision"]["positive_and_negative_modes_share_outer_Morse_phase"] is True,
        "pair-Morse dependency drift",
    )
    require(
        dependencies["discrete_gap"]["decision"]["all_positive_integer_face_minimax_margin_above_22p40"] is True,
        "discrete-gap dependency drift",
    )

    artifact = {
        "kind": STEM,
        "status": "exact_quadratic_B_local_domains_exclude_joint_saddle_with_gradient_margin_above_45p63",
        "passed": True,
        "symbolic_certificate": symbolic_certificate(),
        "interval_certificate": interval_certificate(),
        "decision": {
            "outer_and_normal_phase_is_exactly_quadratic": True,
            "mode_621_residual_domain_excludes_Gaussian_saddle": True,
            "mode_622_residual_domain_excludes_Gaussian_saddle": True,
            "local_exact_max_gradient_component_above_45p63": True,
            "tangent_half_plane_degeneracy_is_exact_domain_degeneracy": False,
            "local_nonstationary_integral_bound_completed": False,
            "complete_B_face_estimate_proved": False,
            "complete_T_upper_proved": False,
            "rh_implication": False,
        },
        "dependencies": {
            name: {"path": relative(path), "sha256": file_hash(path)}
            for name, path in DEPENDENCIES.items()
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "runtime": {
            "elapsed_seconds": round(time.time() - started, 3),
            "workers": 1,
            "flint_threads": 1,
            "process_priority": priority,
        },
        "next_obligation": "Express the two exact local residuals as Gaussian integrals over their curved excluded-saddle domains. Split where 2|y|=pi|q|, apply outer integration by parts on the crossing side and the three-term Fresnel tail on the saddle side, and certify the transformed amplitude and boundary terms.",
        "proof_boundary": "Exact local B-domain phase and gradient gaps only. No amplitude constants, local integral bound, complete B estimate, A-fold splice, complete paired residual, complete T_upper, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified exact quadratic-domain gaps for B modes 621 and 622", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
