#!/usr/bin/env python3
"""Independently check the signed compact Morse--Fresnel overlap integral."""

from __future__ import annotations

import json
import os
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
for candidate in (VENDOR, SCRIPT_ROOT):
    sys.path.insert(0, str(candidate))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import acb, arb, ctx


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_hyperbolic_morse_fresnel_compact_integral_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
BUILDER = REPO_ROOT / f"work/rh_compute/scripts/{STEM}.py"
OVERLAP_GATE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_hyperbolic_morse_fresnel_overlap_core_gate.json"

PRECISION = 100
TOLERANCE = "3e-26"
C = 159_577
MODE = 39_894
Y_CUTOFF = 64
CORE_U = 200
PANELS_PER_HALF = 48
SERIES_DEGREE = 28
SERIES_RADIUS = "0.01"

ROOTS = {
    "lower_face": (
        "-0.002825759532961365076401160417159114247825864097320687324089683704073666283224381896263434854721777473",
        "0.002831092859949456681693112196296051440370287003489608130887631666653078012463105775307871867951718149",
    ),
    "event": (
        "-0.002825759532933649445131464450396983743915033772561613810770360267229064685468657752846686488110707596",
        "0.002831092859921636330917475537313707362085208216652466697578982972472672141177079747032833042459072342",
    ),
    "upper_face": (
        "-0.002825759532905933813862584261780894515755158231241335126296618978417123373928547289289614360234877368",
        "0.002831092859893815980142658768797893714991216035497558963491765267230760723646577968637109603399278758",
    ),
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


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


def parse_complex(record: dict[str, str]) -> acb:
    return acb(arb(record["real_ball"]), arb(record["imag_ball"]))


class IndependentFace:
    def __init__(self, label: str, offset: int) -> None:
        self.label = label
        self.i = acb(0, 1)
        self.pi = arb.pi()
        self.c = arb(C)
        self.m = arb(MODE)
        self.tstar = self.pi * self.c**2 / 8
        self.beta = self.tstar ** (arb(1) / 3)
        self.tau = self.pi * self.m * (self.c - 2 * self.m)
        self.t = self.tau + arb(offset) * self.pi / 16
        self.ymax = arb(Y_CUTOFF) / (2 * self.beta).sqrt()
        self.r = self.t / (2 * self.pi * self.m**2)
        self.den0 = 1 + self.r
        self.endpoint0 = (self.tau - self.t) / (self.pi * self.m)
        self.slope = self.t / (self.pi * self.m)
        self.q_at_zero = self.endpoint0 * (self.pi / (2 * self.den0)).sqrt()
        self.linear = 2 * (self.t - self.tau) / (self.c * self.pi.sqrt())
        self.f_const = (self.i * self.pi / 4).exp() * (self.pi / 2).sqrt()
        self.f_rotate = (-self.i * self.pi / 4).exp() / arb(2).sqrt()
        self.delta_lo = arb(ROOTS[label][0])
        self.delta_hi = arb(ROOTS[label][1])
        z0 = self.beta * self.r.log() / 2
        eta = self.tstar / self.t
        exact_carrier = self.t * (z0 / self.beta - eta * (z0 / self.beta).tanh())
        d = -self.beta / self.c
        canonical_carrier = (self.tstar - self.t) * d / self.beta - d**3 / 3
        self.carrier_phase = exact_carrier - canonical_carrier
        self.carrier = (self.i * self.carrier_phase).exp()

    def k(self, delta: acb) -> acb:
        value = acb(0)
        for n in range(SERIES_DEGREE, -1, -1):
            value = value * delta + arb((-1) ** n) / (n + 2)
        radius = arb(SERIES_RADIUS)
        error = radius ** (SERIES_DEGREE + 1) / ((SERIES_DEGREE + 3) * (1 - radius))
        return value + acb(arb(0, error), arb(0, error))

    def primitive(self, q: acb) -> acb:
        return self.f_const * (self.f_rotate * q).erf()

    def h0(self, a: acb, b: acb) -> acb:
        root = (2 * a).sqrt()
        lower = b / root
        upper = root * self.ymax + lower
        return (-self.i * b**2 / (4 * a)).exp() * (self.primitive(upper) - self.primitive(lower)) / root

    def h1(self, a: acb, b: acb, h0: acb) -> acb:
        edge = (self.i * (a * self.ymax**2 + b * self.ymax)).exp()
        return -(b * h0 + self.i * (edge - 1)) / (2 * a)

    def integrand(self, delta: acb, analytic: bool) -> acb:
        k = self.k(delta)
        den = self.den0 + self.r * delta
        x = 1 / den
        u = delta * (self.t * k).sqrt()
        q = (self.endpoint0 - self.slope * delta) * (self.pi / (2 * den)).sqrt()
        rho = (2 * x).sqrt()

        stable_numerator = (
            -self.pi * self.den0 * self.endpoint0 * self.slope
            - self.pi * (self.endpoint0 * self.endpoint0) * self.r / 2
            + self.t * self.den0 * delta * (self.r - k * den)
        )
        exact_boundary = delta * stable_numerator / (2 * self.den0 * den)
        velocity = (1 + delta) * (2 * k).sqrt()
        amplitude0 = self.r ** (arb(1) / 4) * (1 + delta) ** (-arb(1) / 4) * velocity * rho
        eps = (self.r * x / self.t).sqrt()
        exact_h0 = self.h0(x, rho * q)
        exact_h1 = self.h1(x, rho * q, exact_h0)
        exact = (
            (self.i * exact_boundary).exp()
            * amplitude0
            * ((1 + q * eps) * exact_h0 + rho * eps * exact_h1)
        )

        model_boundary = self.linear * u + u**2 / (2 * self.c) + u**3 / (3 * self.c * self.pi.sqrt())
        model = (self.i * model_boundary).exp() * self.h0(acb(arb(1) / 2), -u)
        du_ddelta = (self.t / 2).sqrt() / velocity
        return (self.carrier * exact - model) * du_ddelta

    def u_real(self, delta: arb) -> arb:
        kval = self.k(acb(delta))
        require(kval.imag.contains(0), "real-axis k enclosure missed zero imaginary part")
        return delta * (self.t * kval.real).sqrt()


def bounds(face: IndependentFace, index: int) -> tuple[arb, arb]:
    if index < PANELS_PER_HALF:
        return (
            face.delta_lo * (PANELS_PER_HALF - index) / PANELS_PER_HALF,
            face.delta_lo * (PANELS_PER_HALF - index - 1) / PANELS_PER_HALF,
        )
    j = index - PANELS_PER_HALF
    return face.delta_hi * j / PANELS_PER_HALF, face.delta_hi * (j + 1) / PANELS_PER_HALF


def integrate(face: IndependentFace) -> acb:
    total = acb(0)
    for index in range(2 * PANELS_PER_HALF):
        left, right = bounds(face, index)
        total += acb.integral(
            face.integrand,
            left,
            right,
            abs_tol=arb(TOLERANCE),
            rel_tol=arb(TOLERANCE),
            eval_limit=450_000,
            depth_limit=48,
        )
    u_error = abs(face.u_real(face.delta_lo) + CORE_U).upper() + abs(face.u_real(face.delta_hi) - CORE_U).upper()
    require(u_error < arb("1e-52"), f"{face.label} independent endpoint error too large")
    return acb(arb(total.real, 3 * u_error), arb(total.imag, 3 * u_error))


def main() -> None:
    ctx.dps = PRECISION
    priority = set_low_priority()
    require(RESULT.is_file(), "missing compact-integral result")
    require(BUILDER.is_file(), "missing compact-integral builder")
    require(OVERLAP_GATE.is_file(), "missing overlap-core dependency")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("status") == "signed_compact_hyperbolic_Morse_Fresnel_overlap_integral_certified", "unexpected result status")
    require(artifact["claims"]["signed_compact_difference_integrated_before_absolute_values"] is True, "signed claim missing")
    require(artifact["claims"]["outer_u_tail_proved"] is False, "outer-tail boundary was promoted")
    require(artifact["claims"]["RH_proved"] is False, "RH boundary was promoted")

    stored = {row["label"]: row for row in artifact["certified_faces"]}
    faces = [IndependentFace("lower_face", -1), IndependentFace("event", 0), IndependentFace("upper_face", 1)]
    maxima = arb(0)
    for face in faces:
        value = integrate(face)
        recorded = parse_complex(stored[face.label]["normalized_compact_difference_ball"])
        require(value.real.overlaps(recorded.real), f"{face.label} real result does not overlap")
        require(value.imag.overlaps(recorded.imag), f"{face.label} imaginary result does not overlap")
        require(abs(value) < arb("0.00049"), f"{face.label} independent normalized target failed")
        require(arb(2).sqrt() / face.pi * abs(value) < arb("0.000221"), f"{face.label} independent physical target failed")
        maxima = max(maxima, abs(value).upper())

    note = (REPO_ROOT / artifact["artifacts"]["note"]["path"]).read_text(encoding="utf-8")
    for token in (
        "cancellation-preserving compact integral certified",
        "|Delta I_core| < 0.00049",
        "does not bound `|u|>200`",
    ):
        require(token in note, f"missing note token: {token}")
    print(
        f"validated signed compact overlap integral independently: faces=3, max |difference|={maxima}; priority={priority}",
        flush=True,
    )


if __name__ == "__main__":
    main()
