#!/usr/bin/env python3
"""Build the Fourier/score-moment regular-field representation and guards."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
import os
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = "jensen_window_pf_newman_fourier_moment_regular_field_gate"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"

SOURCE_PATHS = {
    "regular_field_nonclosure": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_regular_field_first_loss_nonclosure_gate.json",
    "score_bridge": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_one_sided_phase_moment_bridge_gate.json",
    "strong_logconcavity": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_positive_time_strong_logconcavity_gate.json",
    "weighted_shape_guard": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_weighted_strong_logconcavity_countermodel_gate.json",
    "theta_spectral_guard": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_theta_summand_spectral_square_gate.json",
}

MIN_MULTIPLICITY = 2
MAX_AUDIT_MULTIPLICITY = 16
AUDIT_PHASES = (sp.Rational(5, 4) * sp.pi, sp.Rational(3, 2) * sp.pi, sp.Rational(7, 4) * sp.pi)


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    certificate: str
    proof_boundary: str


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def expression_hash(expression: sp.Expr) -> str:
    return hashlib.sha256(sp.srepr(sp.expand(expression)).encode("utf-8")).hexdigest()


def atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(text, encoding="utf-8")
    os.replace(temporary, path)


def load_sources() -> dict[str, dict]:
    payloads: dict[str, dict] = {}
    for key, path in SOURCE_PATHS.items():
        require(path.exists(), f"missing source artifact: {path}")
        payloads[key] = json.loads(path.read_text(encoding="utf-8"))

    require(
        "impose no bound or sign on ell"
        in payloads["regular_field_nonclosure"].get("symbolic_certificate", {}).get("nonclosure", ""),
        "regular-field nonclosure contract drift",
    )
    require(
        "dnu_t(u)=-f_t'(u)du/f_t(0)"
        in payloads["score_bridge"].get("exact", {}).get("score_probability", ""),
        "score probability contract drift",
    )
    require(
        "H_t=f_t(0)*B_t/x"
        in payloads["score_bridge"].get("exact", {}).get("cartesian_lift", ""),
        "score lift contract drift",
    )
    require(
        "(log Phi)''(u)<=-kappa"
        in payloads["strong_logconcavity"].get("exact", {}).get("uniform_phi_curvature", ""),
        "strong-log-concavity contract drift",
    )
    require(
        "theta-type double-exponential decay"
        in payloads["weighted_shape_guard"].get("exact", {}).get("nonpromotion_decision", ""),
        "weighted shape guard drift",
    )
    require(
        "infinite theta/modular endpoint cancellation"
        in payloads["theta_spectral_guard"].get("exact", {}).get("nonpromotion_decision", ""),
        "theta spectral guard drift",
    )
    return payloads


def source_audit() -> dict[str, dict[str, str]]:
    return {
        key: {
            "path": str(path.relative_to(REPO_ROOT)).replace("\\", "/"),
            "sha256": file_hash(path),
        }
        for key, path in SOURCE_PATHS.items()
    }


def sinc(z: sp.Expr) -> sp.Expr:
    return sp.sin(z) / z


def heat_monomial(m: int, z: sp.Symbol, direction: int) -> sp.Expr:
    return sp.expand(
        sum(
            direction**order * sp.diff(z**m, z, 2 * order) / sp.factorial(order)
            for order in range(m // 2 + 1)
        )
    )


def build_moment_audit() -> list[dict]:
    x, u = sp.symbols("x u", real=True)
    rows: list[dict] = []
    for m in range(MIN_MULTIPLICITY, MAX_AUDIT_MULTIPLICITY + 1):
        direct = sp.expand_trig(sp.diff(sp.cos(x * u), x, m))
        score = sp.expand_trig(sp.diff(sp.sin(x * u), x, m))
        if m % 2 == 0:
            r = m // 2
            direct_target = (-1) ** r * u**m * sp.cos(x * u)
            score_target = (-1) ** r * u**m * sp.sin(x * u)
            direct_ratio = "B_c=-S_(m+1)/[(m+1)C_m]"
            score_ratio = "B_c=Cnu_(m+1)/[(m+1)Snu_m]-1/c"
        else:
            r = (m - 1) // 2
            direct_target = (-1) ** (r + 1) * u**m * sp.sin(x * u)
            score_target = (-1) ** r * u**m * sp.cos(x * u)
            direct_ratio = "B_c=C_(m+1)/[(m+1)S_m]"
            score_ratio = "B_c=-Snu_(m+1)/[(m+1)Cnu_m]-1/c"
        require(sp.simplify(direct - direct_target) == 0, f"direct moment parity drift m={m}")
        require(sp.simplify(score - score_target) == 0, f"score moment parity drift m={m}")
        rows.append(
            {
                "multiplicity": m,
                "direct_derivative_sha256": expression_hash(direct),
                "score_derivative_sha256": expression_hash(score),
                "direct_field_ratio": direct_ratio,
                "score_field_ratio": score_ratio,
            }
        )
    return rows


def build_score_jet_audit() -> list[dict]:
    x, c, a, b = sp.symbols("x c a b", nonzero=True)
    h = x - c
    rows: list[dict] = []
    for m in range(MIN_MULTIPLICITY, MAX_AUDIT_MULTIPLICITY + 1):
        score = sp.expand(h**m * (1 + a * h + b * h**2))
        transform = sp.cancel(score / x)
        score_ratio = sp.factor(
            sp.diff(score, x, m + 1).subs(x, c)
            / ((m + 1) * sp.diff(score, x, m).subs(x, c))
        )
        transform_ratio = sp.factor(
            sp.diff(transform, x, m + 1).subs(x, c)
            / ((m + 1) * sp.diff(transform, x, m).subs(x, c))
        )
        ell = sp.factor((2 * c * transform_ratio - m) / 4)
        require(score_ratio == a, f"score regular ratio drift m={m}")
        require(sp.simplify(transform_ratio - (a - 1 / c)) == 0, f"score/transform conversion drift m={m}")
        require(sp.simplify(ell - (c * a / 2 - sp.Rational(m + 2, 4))) == 0, f"score ell drift m={m}")
        rows.append(
            {
                "multiplicity": m,
                "score_ratio": str(score_ratio),
                "transform_ratio": str(transform_ratio),
                "ell": str(ell),
            }
        )
    return rows


def build_sinc_audit() -> list[dict]:
    x = sp.symbols("x")
    c = sp.Integer(1)
    rows: list[dict] = []
    for m in range(2, 7):
        for phase in AUDIT_PHASES:
            transform = sinc(sp.pi * x) ** m * sinc(phase * x)
            derivative_m = sp.simplify(sp.diff(transform, x, m).subs(x, c))
            derivative_next = sp.simplify(sp.diff(transform, x, m + 1).subs(x, c))
            field = sp.simplify(derivative_next / ((m + 1) * derivative_m))
            target = sp.simplify(phase / sp.tan(phase) - m - 1)
            ell = sp.simplify((2 * field - m) / 4)
            target_ell = sp.simplify((2 * phase / sp.tan(phase) - 3 * m - 2) / 4)
            require(derivative_m != 0, f"sinc contact denominator drift m={m}, phase={phase}")
            require(sp.simplify(field - target) == 0, f"sinc field drift m={m}, phase={phase}")
            require(sp.simplify(ell - target_ell) == 0, f"sinc ell drift m={m}, phase={phase}")
            rows.append(
                {
                    "multiplicity": m,
                    "phase": str(phase),
                    "field": str(field),
                    "ell": str(ell),
                    "transform_sha256": expression_hash(transform),
                    "contact_derivative_sha256": expression_hash(derivative_m),
                }
            )
    return rows


def build_gaussian_audit() -> list[dict]:
    x = sp.symbols("x")
    c = sp.Integer(1)
    phase = sp.Rational(5, 4) * sp.pi
    rows: list[dict] = []
    for m in range(2, 7):
        radius = m * sp.pi + phase
        variance = (m + 2) ** 2 * sp.pi**2 + 1
        base = sinc(sp.pi * x) ** m * sinc(phase * x)
        transform = sp.exp(-variance * x**2 / 2) * base
        derivative_m = sp.simplify(sp.diff(transform, x, m).subs(x, c))
        derivative_next = sp.simplify(sp.diff(transform, x, m + 1).subs(x, c))
        field = sp.simplify(derivative_next / ((m + 1) * derivative_m))
        target = sp.simplify(phase / sp.tan(phase) - m - 1 - variance)
        curvature_margin = sp.simplify(variance - radius**2)
        require(sp.simplify(field - target) == 0, f"Gaussian-smoothed field drift m={m}")
        require(curvature_margin.is_positive is True, f"strong curvature margin drift m={m}")
        rows.append(
            {
                "multiplicity": m,
                "phase": str(phase),
                "support_radius": str(radius),
                "gaussian_variance": str(variance),
                "curvature_margin": str(curvature_margin),
                "field": str(field),
                "transform_sha256": expression_hash(transform),
            }
        )
    return rows


def build_hermite_audit() -> list[dict]:
    z = sp.symbols("z")
    rows: list[dict] = []
    for m in range(MIN_MULTIPLICITY, MAX_AUDIT_MULTIPLICITY + 1):
        forward = heat_monomial(m, z, -1)
        backward = heat_monomial(m, z, 1)
        forward_roots = int(sp.Poly(forward, z).count_roots(-sp.oo, sp.oo))
        backward_roots = int(sp.Poly(backward, z).count_roots(-sp.oo, sp.oo))
        require(forward_roots == m, f"forward Hermite drift m={m}")
        require(backward_roots == m % 2, f"backward Hermite drift m={m}")
        rows.append(
            {
                "multiplicity": m,
                "forward_real_roots": forward_roots,
                "backward_real_roots": backward_roots,
                "forward_sha256": expression_hash(forward),
                "backward_sha256": expression_hash(backward),
            }
        )
    return rows


def build_certificate() -> dict[str, str]:
    return {
        "xi_kernel": "f_t(u)=exp(tu^2)Phi(u), H_t(x)=integral_0^infinity f_t(u)cos(xu)du.",
        "complex_moments": (
            "M_k(c,t)=integral_0^infinity u^k f_t(u)exp(icu)du and "
            "H_t^(k)(c)=Re[i^k M_k(c,t)]."
        ),
        "multiple_contact": (
            "At a multiplicity-m zero c, Re[i^jM_j]=0 for 0<=j<m, "
            "Re[i^mM_m]!=0, and B_c=Re[i^(m+1)M_(m+1)]/"
            "((m+1)Re[i^mM_m])."
        ),
        "direct_parity": (
            "For even m, B_c=-S_(m+1)/[(m+1)C_m]; for odd m, "
            "B_c=C_(m+1)/[(m+1)S_m], where C_k=int u^k f cos(cu) and "
            "S_k=int u^k f sin(cu)."
        ),
        "ell_from_field": "ell=(2cB_c-m)/4.",
        "score_probability": (
            "dnu_t(u)=-f_t'(u)du/f_t(0), S_t(x)=E_nu[sin(xU)], and "
            "H_t(x)=f_t(0)S_t(x)/x for x>0."
        ),
        "score_contact": (
            "A multiplicity-m zero of H_t at c>0 is exactly a multiplicity-m "
            "zero of S_t; B_c=S_t^(m+1)(c)/[(m+1)S_t^(m)(c)]-1/c."
        ),
        "score_ell": (
            "ell=c*S_t^(m+1)(c)/[2(m+1)S_t^(m)(c)]-(m+2)/4."
        ),
        "score_parity": (
            "For even m, B_c=Cnu_(m+1)/[(m+1)Snu_m]-1/c; for odd m, "
            "B_c=-Snu_(m+1)/[(m+1)Cnu_m]-1/c."
        ),
        "compact_model": (
            "H_(m,y,c)(x)=sinc(pi*x/c)^m*sinc(y*x/c), pi<y<2pi. It is the "
            "characteristic function of the convolution of m Uniform[-pi/c,pi/c] "
            "laws and one Uniform[-y/c,y/c] law."
        ),
        "compact_field": (
            "At x=c, B_c=(y*cot(y)-m-1)/c and "
            "ell=(2y*cot(y)-3m-2)/4."
        ),
        "phase_bijection": (
            "y*cot(y) decreases continuously from +infinity to -infinity on "
            "(pi,2pi), so the compact positive-kernel model realizes every real ell."
        ),
        "compact_shape": (
            "The frequency density is even, compactly supported, and log-concave "
            "because it is a convolution of interval indicators."
        ),
        "compact_heat": (
            "For every tau, H_tau(x)=int p(u)exp(tau*u^2)exp(ixu)du. For tau>=0, "
            "exp(-tau D_x^2) preserves the Laguerre-Polya class; for tau<0 the "
            "multiplicity-m local limit exp(D_z^2)z^m has a nonreal pair. Thus the "
            "compact model has exact threshold zero."
        ),
        "gaussian_model": (
            "Multiplying H_(m,y,c) by exp(-q*x^2/2) convolves its compact kernel "
            "with N(0,q) and changes B_c by -q*c."
        ),
        "strong_curvature": (
            "If the compact summand Y lies in [-R,R], then (log p_q)''="
            "Var(Y|X=u)/q^2-1/q<=-(q-R^2)/q^2. Choosing "
            "q=((m+2)^2*pi^2+1)/c^2>R^2 makes the kernel smooth, positive, even, "
            "strictly decreasing, and uniformly strongly log-concave."
        ),
        "gaussian_arbitrary_field": (
            "The smoothed model has ell=[2y*cot(y)-3m-2-2q*c^2]/4; the same phase "
            "bijection realizes every real ell and gives a continuous positive score law."
        ),
        "nonclosure": (
            "Positive Fourier representation, compact log-concavity, or smooth positive "
            "uniform strong log-concavity plus local first-loss dynamics do not bound ell."
        ),
        "theta_guard": (
            "The strong model has Gaussian rather than Xi double-exponential tails. The "
            "existing theta-tail weighted countermodel separately forbids shape-and-tail "
            "promotion, but does not by itself realize this exact contact."
        ),
        "xi_handoff": (
            "A surviving field theorem must use the simultaneous Xi theta/modular source "
            "and multiple-contact moment equations, or bypass ell through global first-jet winding."
        ),
        "uniform_handoff": (
            "The cofinal Jensen route still independently requires a remainder uniform in degree."
        ),
    }


def build_rows(c: dict[str, str]) -> list[GateRow]:
    return [
        GateRow("fmr_01_sources", "source chain", "proved", "Five parent contracts are current and hash-pinned.", "Regular-field, score, strong-shape, theta-tail, and modular guards are audited.", "No parent result is strengthened."),
        GateRow("fmr_02_complex", "complex moments", "proved", "Every spatial derivative is one real projection of an adjacent Fourier moment.", c["complex_moments"], "Differentiation is justified by the Xi super-exponential tail."),
        GateRow("fmr_03_contact", "multiple contact", "proved", "Multiplicity and the regular field are exact adjacent-moment conditions.", c["multiple_contact"], "The denominator is nonzero by exact multiplicity."),
        GateRow("fmr_04_parity", "parity split", "proved", "The cosine/sine ratio and its sign depend on the parity of m.", c["direct_parity"], "No absolute-value replacement is allowed."),
        GateRow("fmr_05_ell", "s-coordinate field", "proved", "The Jensen field is recovered exactly from the signed-variable moment ratio.", c["ell_from_field"], "Uses rho=-c^2 and the established field conversion."),
        GateRow("fmr_06_score", "score probability", "proved", "The same contact is encoded by the positive Xi score law.", c["score_probability"], "This imports the proved monotone-kernel integration by parts."),
        GateRow("fmr_07_score_contact", "score jet", "proved", "The full multiplicity transfers through the nonzero factor 1/x.", c["score_contact"], "Only c>0 is used."),
        GateRow("fmr_08_score_ell", "score field ratio", "proved", "The exact Xi-field target is one adjacent score-moment ratio.", c["score_ell"] + " " + c["score_parity"], "The ratio remains signed and oscillatory."),
        GateRow("fmr_09_compact", "positive Fourier model", "proved", "A sinc product supplies a positive compact frequency kernel with an m-fold contact.", c["compact_model"], "This is a countermodel family, not Xi."),
        GateRow("fmr_10_compact_shape", "kernel geometry", "proved", "The compact kernel is even and log-concave.", c["compact_shape"], "It is not smooth at every spline knot."),
        GateRow("fmr_11_compact_field", "arbitrary field", "proved", "The positive-kernel family realizes every real ell.", c["compact_field"] + " " + c["phase_bijection"], "The tuning phase is uniquely selected in (pi,2pi)."),
        GateRow("fmr_12_compact_heat", "exact first loss", "proved", "The compact model has a full all-time heat flow and threshold zero.", c["compact_heat"], "This uses Laguerre-Polya preservation and local Hermite splitting."),
        GateRow("fmr_13_gaussian", "smooth score model", "proved", "Gaussian smoothing keeps the contact and leaves the field freely tunable.", c["gaussian_model"] + " " + c["gaussian_arbitrary_field"], "Its positive-time heat interval ends at the Gaussian integrability threshold."),
        GateRow("fmr_14_curvature", "uniform strong log-concavity", "proved", "The smoothed frequency kernel has an explicit negative curvature ceiling.", c["strong_curvature"], "The tail is Gaussian, not theta-type."),
        GateRow("fmr_15_score_nonclosure", "generic score no-go", "proved", "A continuous positive score probability and strong kernel curvature do not bound ell.", c["nonclosure"], "This does not combine the exact contact with the Xi theta tail."),
        GateRow("fmr_16_theta_guard", "shape/tail guard", "guard_validated", "Generic Xi-like shape and tail information is separately insufficient for correlation positivity.", c["theta_guard"], "The two countermodels are not silently merged into one model."),
        GateRow("fmr_17_pi", "normalization provenance", "proved", "Pi enters only as the first nonzero zero of the uniform characteristic sine.", "sinc(pi*x/c) vanishes at x=c because sin(pi)=0; the Xi Fourier convention itself adds no hidden pi here.", "This pi is from the explicit uniform model, not inserted into the Xi moment identity."),
        GateRow("fmr_18_xi", "Xi arithmetic handoff", "open", "Any closing field restriction must use simultaneous Xi modular arithmetic and contact.", c["xi_handoff"], "No actual Xi field bound is proved."),
        GateRow("fmr_19_uniform", "uniform Jensen handoff", "open", "The local moment formula does not control the cofinal Jensen error.", c["uniform_handoff"], "No degree-uniform remainder is proved."),
        GateRow("fmr_20_boundary", "proof boundary", "guard_validated", "This gate is not a proof of RH.", "It proves the exact moment coordinates and generic positive-kernel/score nonclosure.", "No Xi field bound, collision exclusion, Lambda<=0, RH, or prize-level conclusion is proved."),
    ]


def render_note(artifact: dict) -> str:
    counts = artifact["counts"]
    return f"""# Newman Fourier-Moment Regular-Field Gate

Date: 2026-08-03

Status: exact Xi Fourier/score-moment representation and generic positive-kernel nonclosure; the simultaneous Xi theta-contact theorem and degree-uniform Jensen remainder remain open; not a proof of RH.

## Xi Moment Coordinate

Write

```text
f_t(u)=exp(tu^2)Phi(u),
H_t(x)=integral_0^infinity f_t(u)cos(xu)du,
M_k(c,t)=integral_0^infinity u^k f_t(u)exp(icu)du.
```

The theta tail permits differentiation under the integral, giving

```text
H_t^(k)(c)=Re[i^k M_k(c,t)].
```

At a zero `c>0` of exact multiplicity `m`,

```text
Re[i^j M_j]=0                  (0<=j<m),
Re[i^m M_m]!=0,
B_c=Re[i^(m+1)M_(m+1)]/((m+1)Re[i^m M_m]),
ell=(2cB_c-m)/4.
```

With `C_k=int u^k f_t(u)cos(cu)du` and `S_k=int u^k f_t(u)sin(cu)du`, this becomes

```text
m even: B_c=-S_(m+1)/[(m+1)C_m],
m odd:  B_c= C_(m+1)/[(m+1)S_m].
```

The parity signs are part of the theorem and cannot be replaced by absolute values.

## Positive Score Coordinate

The proved Xi score probability is

```text
dnu_t(u)=-f_t'(u)du/f_t(0),
S_t(x)=E_nu[sin(xU)],
H_t(x)=f_t(0)S_t(x)/x.
```

Thus `S_t` has the same multiplicity at `c`, and

```text
B_c=S_t^(m+1)(c)/[(m+1)S_t^(m)(c)]-1/c,
ell=c*S_t^(m+1)(c)/[2(m+1)S_t^(m)(c)]-(m+2)/4.
```

For even `m`, the score ratio is `Cnu_(m+1)/Snu_m`; for odd `m`, it is `-Snu_(m+1)/Cnu_m`. This is the exact Xi-only quantity a field theorem would have to constrain.

## Positive-Kernel Countermodel

Let `sinc(z)=sin(z)/z`. For `pi<y<2pi`, define

```text
H_(m,y,c)(x)=sinc(pi*x/c)^m sinc(y*x/c).
```

This is the characteristic function of a sum of `m` independent uniform variables on `[-pi/c,pi/c]` and one on `[-y/c,y/c]`. Its frequency kernel is therefore even, nonnegative, compactly supported, and log-concave. At `x=c`,

```text
B_c=(y*cot(y)-m-1)/c,
ell=(2y*cot(y)-3m-2)/4.
```

The map `y*cot(y)` decreases continuously from positive infinity to negative infinity on `(pi,2pi)`. Hence this positive-Fourier-kernel family realizes every real `ell`.

Each sinc factor is Laguerre-Polya. Compact support makes

```text
H_tau(x)=integral p(u)exp(tau*u^2)exp(ixu)du
```

defined for every real `tau`. The factorized operator `exp(-tau D_x^2)` preserves the Laguerre-Polya class for `tau>=0`, while the negative-time local polynomial `exp(D_z^2)z^m` has a nonreal pair. The exact heat threshold is therefore zero.

## Strongly Log-Concave Score Guard

Multiply the transform by `exp(-q*x^2/2)`. In frequency this convolves the compact law with a Gaussian of variance `q`. If its compact summand lies in `[-R,R]`, then

```text
(log p_q)''(u)=Var(Y|X=u)/q^2-1/q
              <=-(q-R^2)/q^2.
```

Choosing

```text
q=((m+2)^2*pi^2+1)/c^2
```

makes `q>R^2` uniformly for every `pi<y<2pi`. The kernel is now smooth, positive, even, strictly decreasing, and uniformly strongly log-concave; its score law has a continuous positive density. Its field is

```text
ell=[2y*cot(y)-3m-2-2q*c^2]/4,
```

which still realizes every real value. Thus Fourier positivity, score-probability positivity, and uniform strong log-concavity, even together with a local first-loss heat split, do not constrain the field.

This smoothed model has Gaussian rather than Xi double-exponential tails. Separately, the promoted weighted countermodel proves that Xi-like strong/root-variable log-concavity and theta-type decay do not imply the needed weighted-correlation sign. These statements are kept separate: no single countermodel here is claimed to combine the exact contact with every Xi tail property.

## What Remains

The surviving target is narrower than before. It must use the simultaneous Xi theta/modular source and the displayed multiplicity-m moment equations, or bypass the field through global first-jet winding. Finite theta blocks and termwise spectral squares are already closed by their endpoint defects. The cofinal Jensen route also still requires a remainder uniform in degree.

Pi provenance is explicit: `pi` appears in the countermodel because `sin(pi)=0`, the first nonzero zero of the uniform characteristic function. No `pi` is inserted into the Xi moment identities themselves.

## Audit

The builder and independent checker cover {counts['multiplicities']} multiplicities, {counts['sinc_models']} exact sinc contacts, {counts['gaussian_models']} Gaussian-smoothed contacts, and {counts['hermite_checks']} forward/backward Hermite splits.

```text
python work/rh_compute/scripts/jensen_window_pf_newman_fourier_moment_regular_field_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_fourier_moment_regular_field_gate.py
```

Current result:

```text
validated Fourier-moment regular-field gate: 20 rows, 5 sources, 15 multiplicities, 15 direct moment checks, 15 score moment checks, 15 score-jet checks, 15 sinc models, 5 Gaussian models, 5 strong-curvature checks, 15 Hermite checks, 2 open handoffs, 0 Xi field bounds, 0 degree-uniform bounds
```
"""


def main() -> None:
    load_sources()
    certificate = build_certificate()
    moment_audit = build_moment_audit()
    score_jet_audit = build_score_jet_audit()
    sinc_audit = build_sinc_audit()
    gaussian_audit = build_gaussian_audit()
    hermite_audit = build_hermite_audit()
    rows = build_rows(certificate)
    artifact = {
        "kind": KIND,
        "date": "2026-08-03",
        "status": "exact Fourier and score moment field representation with generic positive-kernel nonclosure",
        "proof_boundary": (
            "This gate proves the parity-correct Xi Fourier/score-moment representation "
            "of the regular field and exact positive-kernel and strongly-log-concave "
            "families realizing every real field. It proves no bound for the actual Xi "
            "moment ratio, no simultaneous theta-contact theorem, no degree-uniform "
            "Jensen remainder, no Xi collision exclusion, Lambda<=0, RH, or prize-level conclusion."
        ),
        "source_audit": source_audit(),
        "symbolic_certificate": certificate,
        "moment_audit": moment_audit,
        "score_jet_audit": score_jet_audit,
        "sinc_audit": sinc_audit,
        "gaussian_audit": gaussian_audit,
        "hermite_audit": hermite_audit,
        "rows": [asdict(row) for row in rows],
        "counts": {
            "rows": len(rows),
            "sources": len(SOURCE_PATHS),
            "multiplicities": len(moment_audit),
            "direct_moment_checks": len(moment_audit),
            "score_moment_checks": len(moment_audit),
            "score_jet_checks": len(score_jet_audit),
            "sinc_models": len(sinc_audit),
            "gaussian_models": len(gaussian_audit),
            "strong_curvature_checks": len(gaussian_audit),
            "hermite_checks": len(hermite_audit),
            "open_handoffs": 2,
            "actual_xi_field_bounds": 0,
            "degree_uniform_bounds": 0,
            "rh_conclusions": 0,
        },
    }
    atomic_write(RESULT_PATH, json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    atomic_write(NOTE_PATH, render_note(artifact))
    counts = artifact["counts"]
    print(
        "built Fourier-moment regular-field gate: "
        f"{counts['rows']} rows, {counts['sources']} sources, "
        f"{counts['multiplicities']} multiplicities, "
        f"{counts['direct_moment_checks']} direct moment checks, "
        f"{counts['score_moment_checks']} score moment checks, "
        f"{counts['score_jet_checks']} score-jet checks, "
        f"{counts['sinc_models']} sinc models, "
        f"{counts['gaussian_models']} Gaussian models, "
        f"{counts['strong_curvature_checks']} strong-curvature checks, "
        f"{counts['hermite_checks']} Hermite checks, "
        f"{counts['open_handoffs']} open handoffs, "
        f"{counts['actual_xi_field_bounds']} Xi field bounds, "
        f"{counts['degree_uniform_bounds']} degree-uniform bounds"
    )


if __name__ == "__main__":
    main()
