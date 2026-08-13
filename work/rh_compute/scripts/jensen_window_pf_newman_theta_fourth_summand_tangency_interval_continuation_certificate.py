#!/usr/bin/env python3
"""Interval-certify the fourth-theta-summand contact continuation branch.

The certificate uses high-order Taylor models for the oscillatory transforms
and a parameter-uniform Krawczyk operator on each lambda slab.  Pointwise ACB
quadrature preserves the cancellation at x about 136; positive moment
integrals bound every Taylor remainder.
"""

from __future__ import annotations

import argparse
import ctypes
from fractions import Fraction
import hashlib
import json
import math
import os
from pathlib import Path
import sys
from time import perf_counter

for variable in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ.setdefault(variable, "1")

import psutil


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
if str(VENDOR) not in sys.path:
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402
from flint import acb, arb  # noqa: E402


STEM = (
    "jensen_window_pf_newman_theta_fourth_summand_"
    "tangency_interval_continuation_certificate"
)
PARENT_STEM = (
    "jensen_window_pf_newman_theta_fourth_summand_"
    "tangency_homotopy_scout"
)
PARENT_RESULT = REPO_ROOT / "work/rh_compute/results" / f"{PARENT_STEM}.json"
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"

SCHEMA_VERSION = 1
DATE = "2026-08-03"
PRECISION_BITS = 192
ABS_TOL = "1e-34"
INTEGRATION_CUTOFF = Fraction(2)
TAIL_RADIUS = "1e-1500"
TAYLOR_TIME_ORDER = 10
TAYLOR_X_ORDER = 28
MAX_FIELD_DERIVATIVE = 3
TIME_SAFETY_RADIUS = Fraction(1, 250)  # 0.004
X_SAFETY_RADIUS = Fraction(3, 1000)  # 0.003


class ResourcePark(RuntimeError):
    """Raised only after a completed chart when the daytime CPU gate fires."""


class CpuMonitor:
    def __init__(self) -> None:
        self.samples: list[float] = []
        self.consecutive_high = 0

    def sample(self) -> None:
        value = float(psutil.cpu_percent(interval=0.5))
        self.samples.append(value)
        self.consecutive_high = self.consecutive_high + 1 if value > 75.0 else 0
        if self.consecutive_high >= 2:
            raise ResourcePark("two consecutive daytime CPU samples exceeded 75%")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def request_below_normal_priority() -> bool:
    if os.name != "nt":
        try:
            psutil.Process().nice(5)
            return True
        except (psutil.Error, PermissionError, OSError):
            return False
    try:
        kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
        kernel32.GetCurrentProcess.restype = ctypes.c_void_p
        kernel32.SetPriorityClass.argtypes = [ctypes.c_void_p, ctypes.c_uint32]
        kernel32.SetPriorityClass.restype = ctypes.c_int
        return bool(
            kernel32.SetPriorityClass(kernel32.GetCurrentProcess(), 0x00004000)
        )
    except (AttributeError, OSError):
        return False


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json_atomic(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def write_text_atomic(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(text)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def fraction_from_text(value: str) -> Fraction:
    return Fraction(value)


def fraction_text(value: Fraction, digits: int = 80) -> str:
    sign = "-" if value < 0 else ""
    value = abs(value)
    integer, remainder = divmod(value.numerator, value.denominator)
    if remainder == 0:
        return f"{sign}{integer}"
    places: list[str] = []
    for _ in range(digits):
        remainder *= 10
        digit, remainder = divmod(remainder, value.denominator)
        places.append(str(digit))
        if remainder == 0:
            break
    return f"{sign}{integer}." + "".join(places)


def arb_fraction(value: Fraction) -> arb:
    return arb(value.numerator) / value.denominator


def symmetric_ball(radius: Fraction) -> arb:
    return arb(0, arb_fraction(radius))


def arb_text(value: arb) -> str:
    return value.str(70, radius=True, more=True)


def ball_endpoints(value: arb) -> dict[str, str]:
    return {
        "ball": arb_text(value),
        "lower": arb_text(value.lower()),
        "upper": arb_text(value.upper()),
    }


def matrix_inverse_2x2(matrix: list[list[arb]]) -> list[list[arb]]:
    a, b = matrix[0]
    c, d = matrix[1]
    determinant = a * d - b * c
    require(not determinant.contains(0), "point preconditioner matrix is singular")
    return [[d / determinant, -b / determinant], [-c / determinant, a / determinant]]


def matmul_2x2(left: list[list[arb]], right: list[list[arb]]) -> list[list[arb]]:
    return [
        [sum((left[i][k] * right[k][j] for k in range(2)), arb(0)) for j in range(2)]
        for i in range(2)
    ]


def matvec_2(left: list[list[arb]], right: list[arb]) -> list[arb]:
    return [
        sum((left[i][k] * right[k] for k in range(2)), arb(0))
        for i in range(2)
    ]


def arb_power(value: arb, exponent: int) -> arb:
    result = arb(1)
    for _ in range(exponent):
        result *= value
    return result


def load_parent_rows() -> tuple[dict, list[dict]]:
    parent = json.loads(PARENT_RESULT.read_text(encoding="utf-8"))
    rows = parent["continuation"]["rows"]
    require(len(rows) >= 12, "parent continuation does not cover lambda=0 through 0.22")
    for index, row in enumerate(rows[:12]):
        require(
            fraction_from_text(row["lambda"]) == Fraction(index, 50),
            f"unexpected parent lambda at row {index}",
        )
    return parent, rows[:12]


def exact_tail_audit(max_moment: int) -> dict:
    require(max_moment <= 60, "tail audit currently covers moments only through 60")
    exp_two_lower = sum(Fraction(2**k, math.factorial(k)) for k in range(6))
    exp_eight_lower = sum(Fraction(8**k, math.factorial(k)) for k in range(16))
    require(exp_two_lower > 7, "e^2>7 rational audit failed")
    require(exp_eight_lower > 2000, "e^8>2000 rational audit failed")
    require(math.factorial(60) < 10**82, "60!<10^82 integer audit failed")
    require(2560 * 10**82 * 6561**2 < 10**94, "tail prefactor audit failed")
    # alpha*y0 > (3-1/56)*2000 > 5964.  The prefactor is below
    # 10^94 < exp(282), while ln(10)<3 follows from exp(3)>10.
    # Hence every component tail is below exp(-5682)<10^-1800.
    require(Fraction(167, 56) * 2000 > 5964, "tail exponent audit failed")
    require(5964 - 282 > 3 * 1800, "decimal tail conversion audit failed")
    return {
        "range": "u>=2, t<=1/2, 0<=moment<=60, 1<=n<=4",
        "component_bound": "integral tail < 10^-1800",
        "radius_used_per_transform": TAIL_RADIUS,
        "derivation": (
            "u^m<=m!e^u, u^2<=e^(4u)/(4e^2), e^2>7, "
            "3<pi<22/7, 2000<e^8<3^8=6561, and y=e^(4u)"
        ),
        "weighted_family_note": (
            "The n=1,2,3 plus lambda_c*n=4 sum has total tail below "
            "4*10^-1800, hence below the attached 10^-1500 radius."
        ),
    }


class TransformFamily:
    def __init__(
        self,
        lambda_center: Fraction,
        time_center: Fraction,
        time_high: Fraction,
        x_center: Fraction,
        weighted: bool,
    ) -> None:
        self.lambda_center = lambda_center
        self.time_center = time_center
        self.time_high = time_high
        self.x_center = x_center
        self.weighted = weighted
        self.pi = acb.pi()
        self.zero = acb(0)
        self.cutoff = acb(2)
        self.abs_tol = arb(ABS_TOL)
        self.tail_real = arb(f"[0 +/- {TAIL_RADIUS}]")
        self.tail_complex = acb(self.tail_real, self.tail_real)
        self.oscillatory_cache: dict[int, acb] = {}
        self.moment_cache: dict[tuple[Fraction, int], arb] = {}

    def phi_n(self, u: acb, index: int) -> acb:
        n = acb(index)
        a = self.pi * n * n
        e4 = (4 * u).exp()
        return a * (5 * u).exp() * (2 * a * e4 - 3) * (-a * e4).exp()

    def phi(self, u: acb) -> acb:
        if not self.weighted:
            return self.phi_n(u, 4)
        return (
            self.phi_n(u, 1)
            + self.phi_n(u, 2)
            + self.phi_n(u, 3)
            + acb(arb_fraction(self.lambda_center)) * self.phi_n(u, 4)
        )

    def oscillatory_integral(self, moment: int) -> acb:
        cached = self.oscillatory_cache.get(moment)
        if cached is not None:
            return cached
        time_ball = acb(arb_fraction(self.time_center))
        x_ball = acb(arb_fraction(self.x_center))

        def integrand(u: acb, analytic: bool) -> acb:
            del analytic
            return (
                u**moment
                * (time_ball * u * u).exp()
                * self.phi(u)
                * (acb(0, 1) * x_ball * u).exp()
            )

        value = acb.integral(
            integrand,
            self.zero,
            self.cutoff,
            abs_tol=self.abs_tol,
            rel_tol=self.abs_tol,
            deg_limit=64,
            eval_limit=200000,
            depth_limit=64,
        ) + self.tail_complex
        self.oscillatory_cache[moment] = value
        return value

    def moment_upper(self, moment: int, time: Fraction) -> arb:
        key = (time, moment)
        cached = self.moment_cache.get(key)
        if cached is not None:
            return cached
        time_ball = acb(arb_fraction(time))

        def integrand(u: acb, analytic: bool) -> acb:
            del analytic
            return u**moment * (time_ball * u * u).exp() * self.phi(u)

        retained = acb.integral(
            integrand,
            self.zero,
            self.cutoff,
            abs_tol=self.abs_tol,
            rel_tol=self.abs_tol,
            deg_limit=64,
            eval_limit=200000,
            depth_limit=64,
        ).real
        value = retained + self.tail_real
        require(value.upper() > 0, "positive moment integral has no positive upper bound")
        self.moment_cache[key] = value
        return value.upper()

    def point_jet(self, derivative_order: int) -> arb:
        return (acb(0, 1) ** derivative_order * self.oscillatory_integral(derivative_order)).real

    def transform_box(
        self,
        derivative_order: int,
        time_radius: Fraction,
        x_radius: Fraction,
    ) -> tuple[arb, dict[str, str]]:
        tau = symmetric_ball(time_radius)
        xi = symmetric_ball(x_radius)
        polynomial = arb(0)
        for time_order in range(TAYLOR_TIME_ORDER + 1):
            for x_order in range(TAYLOR_X_ORDER + 1):
                moment = derivative_order + 2 * time_order + x_order
                coefficient = (
                    acb(0, 1) ** (derivative_order + x_order)
                    * self.oscillatory_integral(moment)
                ).real
                polynomial += (
                    coefficient
                    * arb_power(tau, time_order)
                    * arb_power(xi, x_order)
                    / (math.factorial(time_order) * math.factorial(x_order))
                )

        time_moment = derivative_order + 2 * (TAYLOR_TIME_ORDER + 1)
        time_remainder = (
            arb_power(arb_fraction(time_radius), TAYLOR_TIME_ORDER + 1)
            / math.factorial(TAYLOR_TIME_ORDER + 1)
            * self.moment_upper(time_moment, self.time_high)
        )

        x_remainder = arb(0)
        for time_order in range(TAYLOR_TIME_ORDER + 1):
            moment = derivative_order + 2 * time_order + TAYLOR_X_ORDER + 1
            x_remainder += (
                arb_power(arb_fraction(time_radius), time_order)
                / math.factorial(time_order)
                * arb_power(arb_fraction(x_radius), TAYLOR_X_ORDER + 1)
                / math.factorial(TAYLOR_X_ORDER + 1)
                * self.moment_upper(moment, self.time_center)
            )

        remainder = (time_remainder + x_remainder).upper()
        result = polynomial + arb(0, remainder)
        return result, {
            "time_remainder_upper": arb_text(time_remainder.upper()),
            "x_remainder_upper": arb_text(x_remainder.upper()),
            "total_remainder_upper": arb_text(remainder),
        }


class ContinuationChart:
    def __init__(self, index: int, left: dict, right: dict) -> None:
        self.index = index
        self.lambda_low = fraction_from_text(left["lambda"])
        self.lambda_high = fraction_from_text(right["lambda"])
        self.lambda_center = (self.lambda_low + self.lambda_high) / 2
        left_time = fraction_from_text(left["time"])
        right_time = fraction_from_text(right["time"])
        left_x = fraction_from_text(left["x"])
        right_x = fraction_from_text(right["x"])
        self.time_center = (left_time + right_time) / 2
        self.x_center = (left_x + right_x) / 2
        self.time_radius = abs(right_time - left_time) / 2 + TIME_SAFETY_RADIUS
        self.x_radius = abs(right_x - left_x) / 2 + X_SAFETY_RADIUS
        self.time_high = self.time_center + self.time_radius
        require(
            Fraction(0) <= self.lambda_center <= Fraction(1),
            "lambda center left the homotopy range",
        )
        require(
            self.time_high <= Fraction(1, 2),
            "chart exceeds the audited tail time range",
        )

        self.weighted = TransformFamily(
            self.lambda_center,
            self.time_center,
            self.time_high,
            self.x_center,
            weighted=True,
        )
        self.fourth = TransformFamily(
            self.lambda_center,
            self.time_center,
            self.time_high,
            self.x_center,
            weighted=False,
        )
        self.weighted_boxes: list[arb] = []
        self.fourth_boxes: list[arb] = []
        self.remainders: list[dict] = []
        for derivative in range(MAX_FIELD_DERIVATIVE + 1):
            weighted_box, weighted_remainder = self.weighted.transform_box(
                derivative, self.time_radius, self.x_radius
            )
            fourth_box, fourth_remainder = self.fourth.transform_box(
                derivative, self.time_radius, self.x_radius
            )
            self.weighted_boxes.append(weighted_box)
            self.fourth_boxes.append(fourth_box)
            self.remainders.append(
                {"weighted": weighted_remainder, "fourth": fourth_remainder}
            )

        self.lambda_radius = (self.lambda_high - self.lambda_low) / 2

    def field_boxes(self, delta_lambda: arb | None = None) -> list[arb]:
        if delta_lambda is None:
            delta_lambda = symmetric_ball(self.lambda_radius)
        return [
            self.weighted_boxes[q] + delta_lambda * self.fourth_boxes[q]
            for q in range(MAX_FIELD_DERIVATIVE + 1)
        ]

    def point_field(self, delta_lambda: arb | None = None) -> list[arb]:
        if delta_lambda is None:
            delta_lambda = symmetric_ball(self.lambda_radius)
        return [
            self.weighted.point_jet(q) + delta_lambda * self.fourth.point_jet(q)
            for q in range(MAX_FIELD_DERIVATIVE + 1)
        ]

    def krawczyk(self, delta_lambda: arb | None = None) -> dict:
        field_box = self.field_boxes(delta_lambda)
        point_field = self.point_field(delta_lambda)
        jacobian_box = [
            [-field_box[2], field_box[1]],
            [-field_box[3], field_box[2]],
        ]
        point_jets = [self.weighted.point_jet(q).mid() for q in range(4)]
        point_jacobian = [
            [-point_jets[2], point_jets[1]],
            [-point_jets[3], point_jets[2]],
        ]
        preconditioner = matrix_inverse_2x2(point_jacobian)
        product = matmul_2x2(preconditioner, jacobian_box)
        defect = [
            [arb(int(i == j)) - product[i][j] for j in range(2)]
            for i in range(2)
        ]
        correction = [-value for value in matvec_2(preconditioner, point_field[:2])]
        delta_box = [symmetric_ball(self.time_radius), symmetric_ball(self.x_radius)]
        defect_term = matvec_2(defect, delta_box)
        image = [correction[i] + defect_term[i] for i in range(2)]
        radii = [arb_fraction(self.time_radius), arb_fraction(self.x_radius)]
        margins = [radii[i] - image[i].abs_upper() for i in range(2)]
        passed = all(margin.lower() > 0 for margin in margins)
        return {
            "passed": passed,
            "image_offsets": [ball_endpoints(value) for value in image],
            "strict_interior_margins": [ball_endpoints(value) for value in margins],
            "preconditioner": [[arb_text(value) for value in row] for row in preconditioner],
            "jacobian_boxes": [[ball_endpoints(value) for value in row] for row in jacobian_box],
            "field_at_center_parameter_box": [ball_endpoints(value) for value in point_field],
        }

    def build(self) -> dict:
        krawczyk = self.krawczyk()
        field_box = self.field_boxes()
        fourth_positive = self.fourth_boxes[0].lower() > 0
        curvature_negative = field_box[2].upper() < 0
        require(krawczyk["passed"], f"chart {self.index} Krawczyk inclusion failed")
        require(fourth_positive, f"chart {self.index} does not prove H4>0")
        require(curvature_negative, f"chart {self.index} does not prove F_xx<0")
        return {
            "index": self.index,
            "lambda": {
                "low": fraction_text(self.lambda_low),
                "center": fraction_text(self.lambda_center),
                "high": fraction_text(self.lambda_high),
                "radius": fraction_text(self.lambda_radius),
            },
            "root_box": {
                "time_center": fraction_text(self.time_center),
                "time_radius": fraction_text(self.time_radius),
                "time_ball": arb_text(
                    arb_fraction(self.time_center) + symmetric_ball(self.time_radius)
                ),
                "x_center": fraction_text(self.x_center),
                "x_radius": fraction_text(self.x_radius),
                "x_ball": arb_text(
                    arb_fraction(self.x_center) + symmetric_ball(self.x_radius)
                ),
            },
            "krawczyk": krawczyk,
            "signs": {
                "H4": ball_endpoints(self.fourth_boxes[0]),
                "F_xx": ball_endpoints(field_box[2]),
                "H4_strictly_positive": fourth_positive,
                "F_xx_strictly_negative": curvature_negative,
                "dt_dlambda_strictly_negative": fourth_positive and curvature_negative,
                "identity": "dt/dlambda=H_4/F_xx at a contact",
            },
            "taylor_remainders": self.remainders,
        }


def pilot_chart(index: int) -> dict:
    _, rows = load_parent_rows()
    require(0 <= index < len(rows) - 1, "pilot chart index is out of range")
    started = perf_counter()
    chart = ContinuationChart(index, rows[index], rows[index + 1])
    result = chart.build()
    result["elapsed_seconds"] = perf_counter() - started
    return result


def connector_record(left: ContinuationChart, right: ContinuationChart) -> dict:
    require(left.lambda_high == right.lambda_low, "connector lambda endpoints differ")
    endpoint_delta = arb_fraction(left.lambda_radius)
    endpoint_krawczyk = left.krawczyk(endpoint_delta)
    require(endpoint_krawczyk["passed"], f"endpoint Krawczyk failed after chart {left.index}")
    root_time = arb_fraction(left.time_center) + arb(
        endpoint_krawczyk["image_offsets"][0]["ball"]
    )
    root_x = arb_fraction(left.x_center) + arb(
        endpoint_krawczyk["image_offsets"][1]["ball"]
    )
    next_time = arb_fraction(right.time_center) + symmetric_ball(right.time_radius)
    next_x = arb_fraction(right.x_center) + symmetric_ball(right.x_radius)
    margins = {
        "time_lower": root_time.lower() - next_time.lower(),
        "time_upper": next_time.upper() - root_time.upper(),
        "x_lower": root_x.lower() - next_x.lower(),
        "x_upper": next_x.upper() - root_x.upper(),
    }
    passed = all(value.lower() > 0 for value in margins.values())
    require(passed, f"connector {left.index}->{right.index} containment failed")
    return {
        "left_chart": left.index,
        "right_chart": right.index,
        "lambda": fraction_text(left.lambda_high),
        "passed": passed,
        "left_endpoint_root_enclosure": {
            "time": ball_endpoints(root_time),
            "x": ball_endpoints(root_x),
        },
        "strict_containment_margins_in_right_chart": {
            key: ball_endpoints(value) for key, value in margins.items()
        },
        "logic": (
            "The left endpoint root lies strictly inside the next chart box; "
            "the next chart's Krawczyk uniqueness therefore identifies both roots."
        ),
    }


def endpoint_root(chart: ContinuationChart, right: bool) -> dict:
    delta = arb_fraction(chart.lambda_radius)
    if not right:
        delta = -delta
    krawczyk = chart.krawczyk(delta)
    require(krawczyk["passed"], f"chart {chart.index} endpoint Krawczyk failed")
    time_ball = arb_fraction(chart.time_center) + arb(
        krawczyk["image_offsets"][0]["ball"]
    )
    x_ball = arb_fraction(chart.x_center) + arb(
        krawczyk["image_offsets"][1]["ball"]
    )
    return {
        "lambda": fraction_text(chart.lambda_high if right else chart.lambda_low),
        "time": ball_endpoints(time_ball),
        "x": ball_endpoints(x_ball),
        "krawczyk": krawczyk,
    }


def certify_crossing(parent: dict, containing_chart: ContinuationChart) -> dict:
    source = parent["boundary_crossing"]
    lambda_center = fraction_from_text(source["lambda"])
    x_center = fraction_from_text(source["x"])
    lambda_radius = Fraction(1, 10**12)
    x_radius = Fraction(1, 10**10)
    time_center = Fraction(0)
    weighted = TransformFamily(
        lambda_center,
        time_center,
        time_center,
        x_center,
        weighted=True,
    )
    fourth = TransformFamily(
        lambda_center,
        time_center,
        time_center,
        x_center,
        weighted=False,
    )
    weighted_boxes: list[arb] = []
    fourth_boxes: list[arb] = []
    for derivative in range(3):
        weighted_box, _ = weighted.transform_box(derivative, Fraction(0), x_radius)
        fourth_box, _ = fourth.transform_box(derivative, Fraction(0), x_radius)
        weighted_boxes.append(weighted_box)
        fourth_boxes.append(fourth_box)

    delta_lambda = symmetric_ball(lambda_radius)
    field_boxes = [
        weighted_boxes[q] + delta_lambda * fourth_boxes[q] for q in range(3)
    ]
    # Here lambda is one of the two Newton unknowns, so G(z_0) is evaluated
    # at the point lambda_center.  Its box variation belongs only in DG(Z).
    point_field = [weighted.point_jet(q) for q in range(2)]
    point_weighted = [weighted.point_jet(q).mid() for q in range(3)]
    point_fourth = [fourth.point_jet(q).mid() for q in range(2)]
    jacobian_box = [
        [fourth_boxes[0], field_boxes[1]],
        [fourth_boxes[1], field_boxes[2]],
    ]
    point_jacobian = [
        [point_fourth[0], point_weighted[1]],
        [point_fourth[1], point_weighted[2]],
    ]
    preconditioner = matrix_inverse_2x2(point_jacobian)
    product = matmul_2x2(preconditioner, jacobian_box)
    defect = [
        [arb(int(i == j)) - product[i][j] for j in range(2)]
        for i in range(2)
    ]
    correction = [-value for value in matvec_2(preconditioner, point_field)]
    delta_box = [symmetric_ball(lambda_radius), symmetric_ball(x_radius)]
    image = [
        correction[i] + matvec_2(defect, delta_box)[i]
        for i in range(2)
    ]
    radii = [arb_fraction(lambda_radius), arb_fraction(x_radius)]
    margins = [radii[i] - image[i].abs_upper() for i in range(2)]
    passed = all(value.lower() > 0 for value in margins)
    require(passed, "crossing Krawczyk inclusion failed")

    lambda_ball = arb_fraction(lambda_center) + image[0]
    x_ball = arb_fraction(x_center) + image[1]
    chart_lambda = arb_fraction(containing_chart.lambda_center) + symmetric_ball(
        containing_chart.lambda_radius
    )
    chart_time = arb_fraction(containing_chart.time_center) + symmetric_ball(
        containing_chart.time_radius
    )
    chart_x = arb_fraction(containing_chart.x_center) + symmetric_ball(
        containing_chart.x_radius
    )
    branch_match = (
        chart_lambda.contains(lambda_ball)
        and chart_time.contains(arb(0))
        and chart_x.contains(x_ball)
    )
    require(branch_match, "crossing enclosure is not inside the continuation chart")
    require(fourth_boxes[0].lower() > 0, "crossing does not prove H4>0")
    require(field_boxes[2].upper() < 0, "crossing does not prove F_xx<0")
    return {
        "passed": passed,
        "lambda_center": fraction_text(lambda_center),
        "lambda_radius": fraction_text(lambda_radius),
        "lambda_enclosure": ball_endpoints(lambda_ball),
        "time": "0",
        "x_center": fraction_text(x_center),
        "x_radius": fraction_text(x_radius),
        "x_enclosure": ball_endpoints(x_ball),
        "krawczyk_image_offsets": [ball_endpoints(value) for value in image],
        "strict_interior_margins": [ball_endpoints(value) for value in margins],
        "jacobian_boxes": [[ball_endpoints(value) for value in row] for row in jacobian_box],
        "H4_strictly_positive": True,
        "F_xx_strictly_negative": True,
        "dt_dlambda_strictly_negative": True,
        "matched_to_continuation_chart": containing_chart.index,
        "branch_match_passed": branch_match,
    }


def build_payload(
    baseline: list[float], priority_lowered: bool, progress: bool = False
) -> tuple[dict, bool]:
    started = perf_counter()
    parent, rows = load_parent_rows()
    monitor = CpuMonitor()
    charts: list[ContinuationChart] = []
    chart_records: list[dict] = []
    parked = False
    for index in range(len(rows) - 1):
        chart = ContinuationChart(index, rows[index], rows[index + 1])
        record = chart.build()
        charts.append(chart)
        chart_records.append(record)
        if progress:
            print(
                f"certified chart {index}: lambda={record['lambda']['low']}.."
                f"{record['lambda']['high']}",
                flush=True,
            )
        try:
            monitor.sample()
        except ResourcePark:
            parked = True
            break

    complete = len(charts) == len(rows) - 1 and not parked
    connectors: list[dict] = []
    crossing: dict | None = None
    endpoints: dict | None = None
    if complete:
        connectors = [
            connector_record(charts[index], charts[index + 1])
            for index in range(len(charts) - 1)
        ]
        initial = endpoint_root(charts[0], right=False)
        final = endpoint_root(charts[-1], right=True)
        require(arb(initial["time"]["ball"]).lower() > 0, "initial branch time is not positive")
        require(arb(final["time"]["ball"]).upper() < 0, "final branch time is not negative")
        endpoints = {
            "initial": initial,
            "final": final,
            "initial_time_strictly_positive": True,
            "final_time_strictly_negative": True,
        }
        crossing = certify_crossing(parent, charts[-1])

    max_moment = MAX_FIELD_DERIVATIVE + 2 * TAYLOR_TIME_ORDER + TAYLOR_X_ORDER + 1
    elapsed = perf_counter() - started
    payload = {
        "kind": STEM,
        "schema_version": SCHEMA_VERSION,
        "date": DATE,
        "status": (
            "rigorous interval continuation certificate complete"
            if complete
            else "parked after resource threshold"
        ),
        "source_sha256": sha256_path(Path(__file__).resolve()),
        "parent": {
            "result": str(PARENT_RESULT.relative_to(REPO_ROOT)).replace("\\", "/"),
            "result_sha256": sha256_path(PARENT_RESULT),
        },
        "resource_policy": {
            "mode": "daytime",
            "active_compute_workers": 1,
            "thread_caps": 1,
            "below_normal_priority_applied": priority_lowered,
            "baseline_cpu_percent": baseline,
            "baseline_mean_percent": sum(baseline) / len(baseline),
            "runtime_cpu_percent": monitor.samples,
            "resource_parked": parked,
            "elapsed_seconds": elapsed,
        },
        "interval_configuration": {
            "precision_bits": PRECISION_BITS,
            "absolute_integration_tolerance": ABS_TOL,
            "integration_interval": [0, 2],
            "tail_radius": TAIL_RADIUS,
            "taylor_time_order": TAYLOR_TIME_ORDER,
            "taylor_x_order": TAYLOR_X_ORDER,
            "maximum_moment": max_moment,
            "time_safety_radius": fraction_text(TIME_SAFETY_RADIUS),
            "x_safety_radius": fraction_text(X_SAFETY_RADIUS),
        },
        "exact_setup": {
            "homotopy": "F_lambda=H_1+H_2+H_3+lambda*H_4",
            "contact": "G=(F_lambda,partial_x F_lambda)=(0,0)",
            "heat_identities": "F_t=-F_xx and (F_x)_t=-F_xxx",
            "contact_jacobian": "[[-F_xx,F_x],[-F_xxx,F_xx]]",
            "implicit_identity_at_contact": "dt/dlambda=H_4/F_xx",
            "parametric_krawczyk_logic": (
                "For every fixed lambda in a slab, the uniform Krawczyk image "
                "lies strictly inside the same (t,x) box, proving one unique contact there."
            ),
            "pi_provenance": (
                "pi is fixed by the Jacobi theta normalization sum_n exp(-pi*n^2*y), "
                "with y=exp(4u); changing pi changes the zeta/Xi kernel."
            ),
        },
        "tail_audit": exact_tail_audit(max_moment),
        "continuation": {
            "lambda_coverage": ["0", "0.22"],
            "expected_charts": len(rows) - 1,
            "certified_charts": len(chart_records),
            "unresolved_charts": (len(rows) - 1) - len(chart_records),
            "charts": chart_records,
            "connectors": connectors,
            "all_krawczyk_inclusions_passed": complete,
            "all_H4_positive": complete,
            "all_F_xx_negative": complete,
            "all_dt_dlambda_negative": complete,
            "branch_chain_connected": complete and len(connectors) == len(charts) - 1,
            "endpoints": endpoints,
        },
        "boundary_crossing": crossing,
        "theorem_decision": {
            "proved_if_complete": (
                "A unique regular contact branch in the certified local chain exists for "
                "0<=lambda<=0.22, has F_xx<0 and dt/dlambda<0 throughout, and crosses "
                "t=0 exactly once in the attached crossing enclosure."
            ),
            "proof_boundary": (
                "This certifies one connected contact branch for the four-term finite theta "
                "homotopy. It does not count all contacts outside these boxes, control lambda>0.22, "
                "add n>=5, exclude a complete-Xi contact, prove Lambda<=0, RH, or the Clay theorem."
            ),
        },
    }
    return payload, complete


def render_note(payload: dict) -> str:
    continuation = payload["continuation"]
    crossing = payload["boundary_crossing"]
    lines = [
        "# Fourth-Summand Tangency Interval Continuation Certificate",
        "",
        f"Date: {DATE}",
        "",
        "Status: rigorous interval continuation certificate for one finite theta homotopy.",
        "It is not a proof of RH and does not represent the complete theta sum.",
        "",
        "```text",
        f"work/rh_compute/results/{STEM}.json",
        f"python work/rh_compute/scripts/{STEM}.py --progress",
        f"python work/rh_compute/scripts/check_{STEM}.py",
        "```",
        "",
        "## Certified Statement",
        "",
        payload["theorem_decision"]["proved_if_complete"],
        "",
        f"The certificate uses `{continuation['certified_charts']}` overlapping parameter",
        f"charts and `{len(continuation['connectors'])}` strict endpoint connectors.",
        "Every chart proves `H_4>0`, `F_xx<0`, and therefore",
        "`dt/dlambda=H_4/F_xx<0` on its unique contact.",
        "",
        "## Crossing",
        "",
        f"The unique `t=0` crossing lies in lambda ball `{crossing['lambda_enclosure']['ball']}`",
        f"and x ball `{crossing['x_enclosure']['ball']}`.",
        "",
        "## Method",
        "",
        "Direct interval substitution erases the high-frequency cancellation. The certificate",
        "instead evaluates point oscillatory moments with ACB, encloses each parameter box by a",
        f"time-order `{TAYLOR_TIME_ORDER}` and frequency-order `{TAYLOR_X_ORDER}` Taylor model,",
        "adds positive-moment remainders and the explicit cutoff tail, and applies a uniform",
        "two-dimensional Krawczyk inclusion on every lambda slab.",
        "",
        "## Pi Provenance",
        "",
        payload["exact_setup"]["pi_provenance"],
        "",
        "## Proof Boundary",
        "",
        payload["theorem_decision"]["proof_boundary"],
        "",
    ]
    return "\n".join(lines)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--pilot-chart", type=int)
    parser.add_argument("--baseline-seconds", type=int, default=5)
    parser.add_argument("--baseline-max-percent", type=float, default=60.0)
    parser.add_argument("--progress", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    priority_lowered = request_below_normal_priority()
    baseline = [float(psutil.cpu_percent(interval=1.0)) for _ in range(args.baseline_seconds)]
    baseline_mean = sum(baseline) / len(baseline)
    print(
        "interval continuation baseline CPU: "
        f"samples={','.join(f'{value:.1f}' for value in baseline)}, "
        f"mean={baseline_mean:.1f}%",
        flush=True,
    )
    if baseline_mean > args.baseline_max_percent:
        print("interval continuation deferred: daytime baseline is already busy")
        return 2
    flint.ctx.prec = PRECISION_BITS
    max_moment = MAX_FIELD_DERIVATIVE + 2 * TAYLOR_TIME_ORDER + TAYLOR_X_ORDER + 1
    tail_audit = exact_tail_audit(max_moment)
    if args.pilot_chart is not None:
        result = pilot_chart(args.pilot_chart)
        print(json.dumps({"pilot": result, "tail_audit": tail_audit}, indent=2))
        return 0
    payload, complete = build_payload(baseline, priority_lowered, args.progress)
    write_json_atomic(args.out, payload)
    if complete:
        write_text_atomic(args.note, render_note(payload))
    print(
        "interval continuation: "
        f"charts={payload['continuation']['certified_charts']}/"
        f"{payload['continuation']['expected_charts']}, "
        f"connectors={len(payload['continuation']['connectors'])}, "
        f"crossing={payload['boundary_crossing'] is not None}, "
        f"parked={payload['resource_policy']['resource_parked']}, "
        f"elapsed={payload['resource_policy']['elapsed_seconds']:.3f}s"
    )
    return 0 if complete else 3


if __name__ == "__main__":
    sys.exit(main())
