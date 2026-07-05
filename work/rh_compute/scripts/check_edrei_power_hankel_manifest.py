#!/usr/bin/env python3
"""Validate finite Edrei-log power-Hankel diagnostic summaries.

This is a necessary-condition diagnostic for the positive power-sum
interpretation of the Edrei-log coefficients.  It is not a proof of
PF-infinity, Laguerre-Polya membership, RH, or Lambda <= 0.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path


@dataclass(frozen=True)
class PowerHankelSpec:
    lam: str
    safe_lam: str
    m_min: int
    m_max: int
    shift_min: int
    shift_max: int
    rows: int
    needed_max_n: int
    dps: int = 340
    suffix: str = ""
    coefficient_frontier: bool = False
    source_log_max_n: int = 0
    source_log_dps: int = 1000
    source_log_suffix: str = "_frontier_tol1e-120"
    boundary_coefficient_frontier: bool = False

    @property
    def stem(self) -> str:
        return (
            f"arb_edrei_power_hankel_lam{self.safe_lam}"
            f"_m{self.m_min}_m{self.m_max}"
            f"_s{self.shift_min}_s{self.shift_max}_dps{self.dps}{self.suffix}"
        )


LAMBDAS: tuple[tuple[str, str], ...] = (
    ("0", "0"),
    ("1e-6", "1em6"),
    ("1e-4", "1em4"),
    ("1e-2", "1em2"),
    ("1e-1", "1em1"),
)

COEFF_LAMBDA_TEXT = {
    "0": "0.0",
    "1e-6": "0.000001",
    "1e-4": "0.0001",
    "1e-2": "0.01",
    "1e-1": "0.1",
}

SPECS: tuple[PowerHankelSpec, ...] = tuple(
    spec
    for lam, safe_lam in LAMBDAS
    for spec in (
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=0,
            m_max=8,
            shift_min=1,
            shift_max=29,
            rows=261,
            needed_max_n=45,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=9,
            m_max=9,
            shift_min=1,
            shift_max=26,
            rows=26,
            needed_max_n=44,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=10,
            m_max=10,
            shift_min=1,
            shift_max=23,
            rows=23,
            needed_max_n=43,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=11,
            m_max=11,
            shift_min=1,
            shift_max=20,
            rows=20,
            needed_max_n=42,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=12,
            m_max=12,
            shift_min=1,
            shift_max=17,
            rows=17,
            needed_max_n=41,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=13,
            m_max=13,
            shift_min=1,
            shift_max=15,
            rows=15,
            needed_max_n=41,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=14,
            m_max=14,
            shift_min=1,
            shift_max=12,
            rows=12,
            needed_max_n=40,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=15,
            m_max=15,
            shift_min=1,
            shift_max=10,
            rows=10,
            needed_max_n=40,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=16,
            m_max=16,
            shift_min=1,
            shift_max=8,
            rows=8,
            needed_max_n=40,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=17,
            m_max=17,
            shift_min=1,
            shift_max=5,
            rows=5,
            needed_max_n=39,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=18,
            m_max=18,
            shift_min=1,
            shift_max=2,
            rows=2,
            needed_max_n=38,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=19,
            m_max=19,
            shift_min=1,
            shift_max=1,
            rows=1,
            needed_max_n=39,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=20,
            m_max=20,
            shift_min=1,
            shift_max=1,
            rows=1,
            needed_max_n=41,
            dps=1000,
            suffix="_frontier_tol1e-120",
            coefficient_frontier=True,
            source_log_max_n=41,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=11,
            m_max=11,
            shift_min=21,
            shift_max=21,
            rows=1,
            needed_max_n=43,
            dps=1000,
            suffix="_frontier_tol1e-120",
            coefficient_frontier=True,
            source_log_max_n=43,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=12,
            m_max=12,
            shift_min=18,
            shift_max=19,
            rows=2,
            needed_max_n=43,
            dps=1000,
            suffix="_frontier_tol1e-120",
            coefficient_frontier=True,
            source_log_max_n=43,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=13,
            m_max=13,
            shift_min=16,
            shift_max=17,
            rows=2,
            needed_max_n=43,
            dps=1000,
            suffix="_frontier_tol1e-120",
            coefficient_frontier=True,
            source_log_max_n=43,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=14,
            m_max=14,
            shift_min=13,
            shift_max=15,
            rows=3,
            needed_max_n=43,
            dps=1000,
            suffix="_frontier_tol1e-120",
            coefficient_frontier=True,
            source_log_max_n=43,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=15,
            m_max=15,
            shift_min=11,
            shift_max=13,
            rows=3,
            needed_max_n=43,
            dps=1000,
            suffix="_frontier_tol1e-120",
            coefficient_frontier=True,
            source_log_max_n=43,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=16,
            m_max=16,
            shift_min=9,
            shift_max=11,
            rows=3,
            needed_max_n=43,
            dps=1000,
            suffix="_frontier_tol1e-120",
            coefficient_frontier=True,
            source_log_max_n=43,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=17,
            m_max=17,
            shift_min=6,
            shift_max=9,
            rows=4,
            needed_max_n=43,
            dps=1000,
            suffix="_frontier_tol1e-120",
            coefficient_frontier=True,
            source_log_max_n=43,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=18,
            m_max=18,
            shift_min=3,
            shift_max=7,
            rows=5,
            needed_max_n=43,
            dps=1000,
            suffix="_frontier_tol1e-120",
            coefficient_frontier=True,
            source_log_max_n=43,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=19,
            m_max=19,
            shift_min=2,
            shift_max=5,
            rows=4,
            needed_max_n=43,
            dps=1000,
            suffix="_frontier_tol1e-120",
            coefficient_frontier=True,
            source_log_max_n=43,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=20,
            m_max=20,
            shift_min=2,
            shift_max=2,
            rows=1,
            needed_max_n=42,
            dps=1000,
            suffix="_frontier_tol1e-120",
            coefficient_frontier=True,
            source_log_max_n=43,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=20,
            m_max=20,
            shift_min=3,
            shift_max=3,
            rows=1,
            needed_max_n=43,
            dps=1000,
            suffix="_frontier_tol1e-120",
            coefficient_frontier=True,
            source_log_max_n=43,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=21,
            m_max=21,
            shift_min=1,
            shift_max=1,
            rows=1,
            needed_max_n=43,
            dps=1000,
            suffix="_frontier_tol1e-120",
            coefficient_frontier=True,
            source_log_max_n=43,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=9,
            m_max=9,
            shift_min=27,
            shift_max=27,
            rows=1,
            needed_max_n=45,
            dps=1000,
            suffix="_frontier_tol1e-120",
            coefficient_frontier=True,
            source_log_max_n=45,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=10,
            m_max=10,
            shift_min=24,
            shift_max=25,
            rows=2,
            needed_max_n=45,
            dps=1000,
            suffix="_frontier_tol1e-120",
            coefficient_frontier=True,
            source_log_max_n=45,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=11,
            m_max=11,
            shift_min=22,
            shift_max=23,
            rows=2,
            needed_max_n=45,
            dps=1000,
            suffix="_frontier_tol1e-120",
            coefficient_frontier=True,
            source_log_max_n=45,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=12,
            m_max=12,
            shift_min=20,
            shift_max=21,
            rows=2,
            needed_max_n=45,
            dps=1000,
            suffix="_frontier_tol1e-120",
            coefficient_frontier=True,
            source_log_max_n=45,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=13,
            m_max=13,
            shift_min=18,
            shift_max=19,
            rows=2,
            needed_max_n=45,
            dps=1000,
            suffix="_frontier_tol1e-120",
            coefficient_frontier=True,
            source_log_max_n=45,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=14,
            m_max=14,
            shift_min=16,
            shift_max=17,
            rows=2,
            needed_max_n=45,
            dps=1000,
            suffix="_frontier_tol1e-120",
            coefficient_frontier=True,
            source_log_max_n=45,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=15,
            m_max=15,
            shift_min=14,
            shift_max=15,
            rows=2,
            needed_max_n=45,
            dps=1000,
            suffix="_frontier_tol1e-120",
            coefficient_frontier=True,
            source_log_max_n=45,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=16,
            m_max=16,
            shift_min=12,
            shift_max=13,
            rows=2,
            needed_max_n=45,
            dps=1000,
            suffix="_frontier_tol1e-120",
            coefficient_frontier=True,
            source_log_max_n=45,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=17,
            m_max=17,
            shift_min=10,
            shift_max=11,
            rows=2,
            needed_max_n=45,
            dps=1000,
            suffix="_frontier_tol1e-120",
            coefficient_frontier=True,
            source_log_max_n=45,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=18,
            m_max=18,
            shift_min=8,
            shift_max=9,
            rows=2,
            needed_max_n=45,
            dps=1000,
            suffix="_frontier_tol1e-120",
            coefficient_frontier=True,
            source_log_max_n=45,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=19,
            m_max=19,
            shift_min=6,
            shift_max=7,
            rows=2,
            needed_max_n=45,
            dps=1000,
            suffix="_frontier_tol1e-120",
            coefficient_frontier=True,
            source_log_max_n=45,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=20,
            m_max=20,
            shift_min=4,
            shift_max=5,
            rows=2,
            needed_max_n=45,
            dps=1000,
            suffix="_frontier_tol1e-120",
            coefficient_frontier=True,
            source_log_max_n=45,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=21,
            m_max=21,
            shift_min=2,
            shift_max=3,
            rows=2,
            needed_max_n=45,
            dps=1000,
            suffix="_frontier_tol1e-120",
            coefficient_frontier=True,
            source_log_max_n=45,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=22,
            m_max=22,
            shift_min=1,
            shift_max=1,
            rows=1,
            needed_max_n=45,
            dps=1000,
            suffix="_frontier_tol1e-120",
            coefficient_frontier=True,
            source_log_max_n=45,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=9,
            m_max=9,
            shift_min=28,
            shift_max=29,
            rows=2,
            needed_max_n=47,
            dps=1000,
            suffix="_frontier_tol1e-120",
            coefficient_frontier=True,
            source_log_max_n=47,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=10,
            m_max=10,
            shift_min=26,
            shift_max=27,
            rows=2,
            needed_max_n=47,
            dps=1000,
            suffix="_frontier_tol1e-120",
            coefficient_frontier=True,
            source_log_max_n=47,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=11,
            m_max=11,
            shift_min=24,
            shift_max=25,
            rows=2,
            needed_max_n=47,
            dps=1000,
            suffix="_frontier_tol1e-120",
            coefficient_frontier=True,
            source_log_max_n=47,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=12,
            m_max=12,
            shift_min=22,
            shift_max=23,
            rows=2,
            needed_max_n=47,
            dps=1000,
            suffix="_frontier_tol1e-120",
            coefficient_frontier=True,
            source_log_max_n=47,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=13,
            m_max=13,
            shift_min=20,
            shift_max=21,
            rows=2,
            needed_max_n=47,
            dps=1000,
            suffix="_frontier_tol1e-120",
            coefficient_frontier=True,
            source_log_max_n=47,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=14,
            m_max=14,
            shift_min=18,
            shift_max=19,
            rows=2,
            needed_max_n=47,
            dps=1000,
            suffix="_frontier_tol1e-120",
            coefficient_frontier=True,
            source_log_max_n=47,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=15,
            m_max=15,
            shift_min=16,
            shift_max=17,
            rows=2,
            needed_max_n=47,
            dps=1000,
            suffix="_frontier_tol1e-120",
            coefficient_frontier=True,
            source_log_max_n=47,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=16,
            m_max=16,
            shift_min=14,
            shift_max=15,
            rows=2,
            needed_max_n=47,
            dps=1000,
            suffix="_frontier_tol1e-120",
            coefficient_frontier=True,
            source_log_max_n=47,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=17,
            m_max=17,
            shift_min=12,
            shift_max=13,
            rows=2,
            needed_max_n=47,
            dps=1000,
            suffix="_frontier_tol1e-120",
            coefficient_frontier=True,
            source_log_max_n=47,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=18,
            m_max=18,
            shift_min=10,
            shift_max=11,
            rows=2,
            needed_max_n=47,
            dps=1000,
            suffix="_frontier_tol1e-120",
            coefficient_frontier=True,
            source_log_max_n=47,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=19,
            m_max=19,
            shift_min=8,
            shift_max=9,
            rows=2,
            needed_max_n=47,
            dps=1000,
            suffix="_frontier_tol1e-120",
            coefficient_frontier=True,
            source_log_max_n=47,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=20,
            m_max=20,
            shift_min=6,
            shift_max=7,
            rows=2,
            needed_max_n=47,
            dps=1000,
            suffix="_frontier_tol1e-120",
            coefficient_frontier=True,
            source_log_max_n=47,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=21,
            m_max=21,
            shift_min=4,
            shift_max=5,
            rows=2,
            needed_max_n=47,
            dps=1000,
            suffix="_frontier_tol1e-120",
            coefficient_frontier=True,
            source_log_max_n=47,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=22,
            m_max=22,
            shift_min=2,
            shift_max=3,
            rows=2,
            needed_max_n=47,
            dps=1000,
            suffix="_frontier_tol1e-120",
            coefficient_frontier=True,
            source_log_max_n=47,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=23,
            m_max=23,
            shift_min=1,
            shift_max=1,
            rows=1,
            needed_max_n=47,
            dps=1000,
            suffix="_frontier_tol1e-120",
            coefficient_frontier=True,
            source_log_max_n=47,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=10,
            m_max=10,
            shift_min=28,
            shift_max=29,
            rows=2,
            needed_max_n=49,
            dps=1000,
            suffix="_frontier_tol1e-120",
            coefficient_frontier=True,
            source_log_max_n=49,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=11,
            m_max=11,
            shift_min=26,
            shift_max=27,
            rows=2,
            needed_max_n=49,
            dps=1000,
            suffix="_frontier_tol1e-120",
            coefficient_frontier=True,
            source_log_max_n=49,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=12,
            m_max=12,
            shift_min=24,
            shift_max=25,
            rows=2,
            needed_max_n=49,
            dps=1000,
            suffix="_frontier_tol1e-120",
            coefficient_frontier=True,
            source_log_max_n=49,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=13,
            m_max=13,
            shift_min=22,
            shift_max=23,
            rows=2,
            needed_max_n=49,
            dps=1000,
            suffix="_frontier_tol1e-120",
            coefficient_frontier=True,
            source_log_max_n=49,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=14,
            m_max=14,
            shift_min=20,
            shift_max=21,
            rows=2,
            needed_max_n=49,
            dps=1000,
            suffix="_frontier_tol1e-120",
            coefficient_frontier=True,
            source_log_max_n=49,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=15,
            m_max=15,
            shift_min=18,
            shift_max=19,
            rows=2,
            needed_max_n=49,
            dps=1000,
            suffix="_frontier_tol1e-120",
            coefficient_frontier=True,
            source_log_max_n=49,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=16,
            m_max=16,
            shift_min=16,
            shift_max=17,
            rows=2,
            needed_max_n=49,
            dps=1000,
            suffix="_frontier_tol1e-120",
            coefficient_frontier=True,
            source_log_max_n=49,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=17,
            m_max=17,
            shift_min=14,
            shift_max=15,
            rows=2,
            needed_max_n=49,
            dps=1000,
            suffix="_frontier_tol1e-120",
            coefficient_frontier=True,
            source_log_max_n=49,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=18,
            m_max=18,
            shift_min=12,
            shift_max=13,
            rows=2,
            needed_max_n=49,
            dps=1000,
            suffix="_frontier_tol1e-120",
            coefficient_frontier=True,
            source_log_max_n=49,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=19,
            m_max=19,
            shift_min=10,
            shift_max=11,
            rows=2,
            needed_max_n=49,
            dps=1000,
            suffix="_frontier_tol1e-120",
            coefficient_frontier=True,
            source_log_max_n=49,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=20,
            m_max=20,
            shift_min=8,
            shift_max=9,
            rows=2,
            needed_max_n=49,
            dps=1000,
            suffix="_frontier_tol1e-120",
            coefficient_frontier=True,
            source_log_max_n=49,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=21,
            m_max=21,
            shift_min=6,
            shift_max=7,
            rows=2,
            needed_max_n=49,
            dps=1000,
            suffix="_frontier_tol1e-120",
            coefficient_frontier=True,
            source_log_max_n=49,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=22,
            m_max=22,
            shift_min=4,
            shift_max=5,
            rows=2,
            needed_max_n=49,
            dps=1000,
            suffix="_frontier_tol1e-120",
            coefficient_frontier=True,
            source_log_max_n=49,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=23,
            m_max=23,
            shift_min=2,
            shift_max=3,
            rows=2,
            needed_max_n=49,
            dps=1000,
            suffix="_frontier_tol1e-120",
            coefficient_frontier=True,
            source_log_max_n=49,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=24,
            m_max=24,
            shift_min=1,
            shift_max=1,
            rows=1,
            needed_max_n=49,
            dps=1000,
            suffix="_frontier_tol1e-120",
            coefficient_frontier=True,
            source_log_max_n=49,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=0,
            m_max=8,
            shift_min=30,
            shift_max=33,
            rows=36,
            needed_max_n=49,
            dps=1000,
            suffix="_frontier_tol1e-120",
            coefficient_frontier=True,
            source_log_max_n=49,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=0,
            m_max=7,
            shift_min=34,
            shift_max=35,
            rows=16,
            needed_max_n=49,
            dps=1000,
            suffix="_frontier_tol1e-120",
            coefficient_frontier=True,
            source_log_max_n=49,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=0,
            m_max=6,
            shift_min=36,
            shift_max=37,
            rows=14,
            needed_max_n=49,
            dps=1000,
            suffix="_frontier_tol1e-120",
            coefficient_frontier=True,
            source_log_max_n=49,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=0,
            m_max=5,
            shift_min=38,
            shift_max=39,
            rows=12,
            needed_max_n=49,
            dps=1000,
            suffix="_frontier_tol1e-120",
            coefficient_frontier=True,
            source_log_max_n=49,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=0,
            m_max=4,
            shift_min=40,
            shift_max=41,
            rows=10,
            needed_max_n=49,
            dps=1000,
            suffix="_frontier_tol1e-120",
            coefficient_frontier=True,
            source_log_max_n=49,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=0,
            m_max=3,
            shift_min=42,
            shift_max=43,
            rows=8,
            needed_max_n=49,
            dps=1000,
            suffix="_frontier_tol1e-120",
            coefficient_frontier=True,
            source_log_max_n=49,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=0,
            m_max=2,
            shift_min=44,
            shift_max=45,
            rows=6,
            needed_max_n=49,
            dps=1000,
            suffix="_frontier_tol1e-120",
            coefficient_frontier=True,
            source_log_max_n=49,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=0,
            m_max=1,
            shift_min=46,
            shift_max=47,
            rows=4,
            needed_max_n=49,
            dps=1000,
            suffix="_frontier_tol1e-120",
            coefficient_frontier=True,
            source_log_max_n=49,
        ),
        PowerHankelSpec(
            lam=lam,
            safe_lam=safe_lam,
            m_min=0,
            m_max=0,
            shift_min=48,
            shift_max=49,
            rows=2,
            needed_max_n=49,
            dps=1000,
            suffix="_frontier_tol1e-120",
            coefficient_frontier=True,
            source_log_max_n=49,
        ),
    )
)

FRONTIER_51_BLOCKS: tuple[tuple[int, int, int, int, int], ...] = (
    (0, 0, 50, 51, 2),
    (1, 1, 48, 49, 2),
    (2, 2, 46, 47, 2),
    (3, 3, 44, 45, 2),
    (4, 4, 42, 43, 2),
    (5, 5, 40, 41, 2),
    (6, 6, 38, 39, 2),
    (7, 7, 36, 37, 2),
    (8, 8, 34, 35, 2),
    (9, 9, 30, 33, 4),
    (10, 10, 30, 31, 2),
    (11, 11, 28, 29, 2),
    (12, 12, 26, 27, 2),
    (13, 13, 24, 25, 2),
    (14, 14, 22, 23, 2),
    (15, 15, 20, 21, 2),
    (16, 16, 18, 19, 2),
    (17, 17, 16, 17, 2),
    (18, 18, 14, 15, 2),
    (19, 19, 12, 13, 2),
    (20, 20, 10, 11, 2),
    (21, 21, 8, 9, 2),
    (22, 22, 6, 7, 2),
    (23, 23, 4, 5, 2),
)

FRONTIER_51_REPAIR_BLOCKS: tuple[tuple[int, int, int, int, int], ...] = (
    (24, 24, 2, 3, 2),
    (25, 25, 1, 1, 1),
)

FRONTIER_53_BLOCKS: tuple[tuple[int, int, int, int, int], ...] = (
    (0, 0, 52, 53, 2),
    (1, 1, 50, 51, 2),
    (2, 2, 48, 49, 2),
    (3, 3, 46, 47, 2),
    (4, 4, 44, 45, 2),
    (5, 5, 42, 43, 2),
    (6, 6, 40, 41, 2),
    (7, 7, 38, 39, 2),
    (8, 8, 36, 37, 2),
    (9, 9, 34, 35, 2),
    (10, 10, 32, 33, 2),
    (11, 11, 30, 31, 2),
    (12, 12, 28, 29, 2),
    (13, 13, 26, 27, 2),
    (14, 14, 24, 25, 2),
    (15, 15, 22, 23, 2),
    (16, 16, 20, 21, 2),
    (17, 17, 18, 19, 2),
    (18, 18, 16, 17, 2),
    (19, 19, 14, 15, 2),
    (20, 20, 12, 13, 2),
    (21, 21, 10, 11, 2),
    (22, 22, 8, 9, 2),
    (23, 23, 6, 7, 2),
    (24, 24, 4, 5, 2),
    (25, 25, 2, 3, 2),
    (26, 26, 1, 1, 1),
)

FRONTIER_55_BLOCKS: tuple[tuple[int, int, int, int, int], ...] = (
    (0, 0, 54, 55, 2),
    (1, 1, 52, 53, 2),
    (2, 2, 50, 51, 2),
    (3, 3, 48, 49, 2),
    (4, 4, 46, 47, 2),
    (5, 5, 44, 45, 2),
    (6, 6, 42, 43, 2),
    (7, 7, 40, 41, 2),
    (8, 8, 38, 39, 2),
    (9, 9, 36, 37, 2),
    (10, 10, 34, 35, 2),
    (11, 11, 32, 33, 2),
    (12, 12, 30, 31, 2),
    (13, 13, 28, 29, 2),
    (14, 14, 26, 27, 2),
    (15, 15, 24, 25, 2),
    (16, 16, 22, 23, 2),
    (17, 17, 20, 21, 2),
    (18, 18, 18, 19, 2),
    (19, 19, 16, 17, 2),
    (20, 20, 14, 15, 2),
    (21, 21, 12, 13, 2),
    (22, 22, 10, 11, 2),
    (23, 23, 8, 9, 2),
    (24, 24, 6, 7, 2),
    (25, 25, 4, 5, 2),
    (26, 26, 2, 3, 2),
    (27, 27, 1, 1, 1),
)

FRONTIER_57_BLOCKS: tuple[tuple[int, int, int, int, int], ...] = (
    (0, 0, 56, 57, 2),
    (1, 1, 54, 55, 2),
    (2, 2, 52, 53, 2),
    (3, 3, 50, 51, 2),
    (4, 4, 48, 49, 2),
    (5, 5, 46, 47, 2),
    (6, 6, 44, 45, 2),
    (7, 7, 42, 43, 2),
    (8, 8, 40, 41, 2),
    (9, 9, 38, 39, 2),
    (10, 10, 36, 37, 2),
    (11, 11, 34, 35, 2),
    (12, 12, 32, 33, 2),
    (13, 13, 30, 31, 2),
    (14, 14, 28, 29, 2),
    (15, 15, 26, 27, 2),
    (16, 16, 24, 25, 2),
    (17, 17, 22, 23, 2),
    (18, 18, 20, 21, 2),
    (19, 19, 18, 19, 2),
    (20, 20, 16, 17, 2),
    (21, 21, 14, 15, 2),
    (22, 22, 12, 13, 2),
    (23, 23, 10, 11, 2),
    (24, 24, 8, 9, 2),
    (25, 25, 6, 7, 2),
    (26, 26, 4, 5, 2),
    (27, 27, 2, 3, 2),
    (28, 28, 1, 1, 1),
)

SPECS = SPECS + tuple(
    PowerHankelSpec(
        lam=lam,
        safe_lam=safe_lam,
        m_min=m_min,
        m_max=m_max,
        shift_min=shift_min,
        shift_max=shift_max,
        rows=rows,
        needed_max_n=51,
        dps=1000,
        suffix="_frontier_tol1e-120",
        coefficient_frontier=True,
        source_log_max_n=51,
    )
    for lam, safe_lam in LAMBDAS
    for m_min, m_max, shift_min, shift_max, rows in FRONTIER_51_BLOCKS
)

SPECS = SPECS + tuple(
    PowerHankelSpec(
        lam=lam,
        safe_lam=safe_lam,
        m_min=m_min,
        m_max=m_max,
        shift_min=shift_min,
        shift_max=shift_max,
        rows=rows,
        needed_max_n=51,
        dps=1000,
        suffix="_frontier_tol1e-120",
        coefficient_frontier=True,
        source_log_max_n=51,
    )
    for lam, safe_lam in LAMBDAS
    if lam != "1e-6"
    for m_min, m_max, shift_min, shift_max, rows in FRONTIER_51_REPAIR_BLOCKS
)

SPECS = SPECS + tuple(
    PowerHankelSpec(
        lam="1e-6",
        safe_lam="1em6",
        m_min=m_min,
        m_max=m_max,
        shift_min=shift_min,
        shift_max=shift_max,
        rows=rows,
        needed_max_n=51,
        dps=2400,
        suffix="_boundary_tol1e-140",
        coefficient_frontier=True,
        source_log_max_n=51,
        source_log_dps=2400,
        source_log_suffix="_boundary_tol1e-140",
        boundary_coefficient_frontier=True,
    )
    for m_min, m_max, shift_min, shift_max, rows in FRONTIER_51_REPAIR_BLOCKS
)

SPECS = SPECS + tuple(
    PowerHankelSpec(
        lam=lam,
        safe_lam=safe_lam,
        m_min=m_min,
        m_max=m_max,
        shift_min=shift_min,
        shift_max=shift_max,
        rows=rows,
        needed_max_n=53,
        dps=2400,
        suffix="_boundary_tol1e-140",
        coefficient_frontier=True,
        source_log_max_n=53,
        source_log_dps=2400,
        source_log_suffix="_boundary_tol1e-140",
        boundary_coefficient_frontier=True,
    )
    for lam, safe_lam in LAMBDAS
    for m_min, m_max, shift_min, shift_max, rows in FRONTIER_53_BLOCKS
)

SPECS = SPECS + tuple(
    PowerHankelSpec(
        lam=lam,
        safe_lam=safe_lam,
        m_min=m_min,
        m_max=m_max,
        shift_min=shift_min,
        shift_max=shift_max,
        rows=rows,
        needed_max_n=55,
        dps=2400,
        suffix="_boundary_tol1e-140",
        coefficient_frontier=True,
        source_log_max_n=55,
        source_log_dps=2400,
        source_log_suffix="_boundary_tol1e-140",
        boundary_coefficient_frontier=True,
    )
    for lam, safe_lam in LAMBDAS
    for m_min, m_max, shift_min, shift_max, rows in FRONTIER_55_BLOCKS
)

SPECS = SPECS + tuple(
    PowerHankelSpec(
        lam=lam,
        safe_lam=safe_lam,
        m_min=m_min,
        m_max=m_max,
        shift_min=shift_min,
        shift_max=shift_max,
        rows=rows,
        needed_max_n=57,
        dps=2400,
        suffix="_boundary_tol1e-140",
        coefficient_frontier=True,
        source_log_max_n=57,
        source_log_dps=2400,
        source_log_suffix="_boundary_tol1e-140",
        boundary_coefficient_frontier=True,
    )
    for lam, safe_lam in LAMBDAS
    for m_min, m_max, shift_min, shift_max, rows in FRONTIER_57_BLOCKS
)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_spec(results_dir: Path, spec: PowerHankelSpec) -> None:
    if spec.coefficient_frontier:
        coeff_summary_path = (
            results_dir
            / f"acb_enclosures_edrei_frontier_lam{spec.safe_lam}_k0_k41_dps180_tol1e-120_summary.json"
        )
        coeff_summary = load_json(coeff_summary_path)
        coeff_checks = {
            "kind": "acb_coefficient_enclosure_summary",
            "k_min": 0,
            "k_max": 41,
            "dps": 180,
            "abs_tol": "1e-120",
            "rows": 42,
        }
        for key, expected in coeff_checks.items():
            actual = coeff_summary.get(key)
            if actual != expected:
                raise AssertionError(
                    f"{coeff_summary_path.name}: {key} {actual!r} != {expected!r}"
                )
        if coeff_summary.get("lambdas") != [COEFF_LAMBDA_TEXT[spec.lam]]:
            raise AssertionError(
                f"{coeff_summary_path.name}: unexpected lambdas {coeff_summary.get('lambdas')!r}"
            )

        log_max_n = spec.source_log_max_n or spec.needed_max_n
        if log_max_n > 41:
            coeff_tail_path = (
                results_dir
                / "acb_enclosures_edrei_frontier_lamgrid_k42_k43_dps180_tol1e-120_summary.json"
            )
            coeff_tail_summary = load_json(coeff_tail_path)
            coeff_tail_checks = {
                "kind": "acb_coefficient_enclosure_summary",
                "k_min": 42,
                "k_max": 43,
                "dps": 180,
                "abs_tol": "1e-120",
                "rows": 10,
            }
            for key, expected in coeff_tail_checks.items():
                actual = coeff_tail_summary.get(key)
                if actual != expected:
                    raise AssertionError(
                        f"{coeff_tail_path.name}: {key} {actual!r} != {expected!r}"
                    )
        if log_max_n > 43:
            coeff_tail2_path = (
                results_dir
                / "acb_enclosures_edrei_frontier_lamgrid_k44_k45_dps180_tol1e-120_summary.json"
            )
            coeff_tail2_summary = load_json(coeff_tail2_path)
            coeff_tail2_checks = {
                "kind": "acb_coefficient_enclosure_summary",
                "k_min": 44,
                "k_max": 45,
                "dps": 180,
                "abs_tol": "1e-120",
                "rows": 10,
            }
            for key, expected in coeff_tail2_checks.items():
                actual = coeff_tail2_summary.get(key)
                if actual != expected:
                    raise AssertionError(
                        f"{coeff_tail2_path.name}: {key} {actual!r} != {expected!r}"
                    )
        if log_max_n > 45:
            coeff_tail3_path = (
                results_dir
                / "acb_enclosures_edrei_frontier_lamgrid_k46_k47_dps180_tol1e-120_summary.json"
            )
            coeff_tail3_summary = load_json(coeff_tail3_path)
            coeff_tail3_checks = {
                "kind": "acb_coefficient_enclosure_summary",
                "k_min": 46,
                "k_max": 47,
                "dps": 180,
                "abs_tol": "1e-120",
                "rows": 10,
            }
            for key, expected in coeff_tail3_checks.items():
                actual = coeff_tail3_summary.get(key)
                if actual != expected:
                    raise AssertionError(
                        f"{coeff_tail3_path.name}: {key} {actual!r} != {expected!r}"
                    )
        if log_max_n > 47:
            coeff_tail4_path = (
                results_dir
                / "acb_enclosures_edrei_frontier_lamgrid_k48_k49_dps180_tol1e-120_summary.json"
            )
            coeff_tail4_summary = load_json(coeff_tail4_path)
            coeff_tail4_checks = {
                "kind": "acb_coefficient_enclosure_summary",
                "k_min": 48,
                "k_max": 49,
                "dps": 180,
                "abs_tol": "1e-120",
                "rows": 10,
            }
            for key, expected in coeff_tail4_checks.items():
                actual = coeff_tail4_summary.get(key)
                if actual != expected:
                    raise AssertionError(
                        f"{coeff_tail4_path.name}: {key} {actual!r} != {expected!r}"
                    )
        if log_max_n > 49:
            coeff_tail5_path = (
                results_dir
                / "acb_enclosures_edrei_frontier_lamgrid_k50_k51_dps180_tol1e-120_summary.json"
            )
            coeff_tail5_summary = load_json(coeff_tail5_path)
            coeff_tail5_checks = {
                "kind": "acb_coefficient_enclosure_summary",
                "k_min": 50,
                "k_max": 51,
                "dps": 180,
                "abs_tol": "1e-120",
                "rows": 10,
            }
            for key, expected in coeff_tail5_checks.items():
                actual = coeff_tail5_summary.get(key)
                if actual != expected:
                    raise AssertionError(
                        f"{coeff_tail5_path.name}: {key} {actual!r} != {expected!r}"
                    )
        if spec.boundary_coefficient_frontier:
            boundary_k_max = log_max_n
            boundary_sources = [(0, boundary_k_max)]
            if boundary_k_max == 55:
                boundary_sources = [(0, 53), (54, 55)]
            if boundary_k_max == 57:
                boundary_sources = [(0, 53), (54, 55), (56, 57)]
            for k_min, k_max in boundary_sources:
                coeff_boundary_path = (
                    results_dir
                    / (
                        f"acb_enclosures_edrei_boundary_lam{spec.safe_lam}_"
                        f"k{k_min}_k{k_max}_dps220_tol1e-140_summary.json"
                    )
                )
                coeff_boundary_summary = load_json(coeff_boundary_path)
                coeff_boundary_checks = {
                    "kind": "acb_coefficient_enclosure_summary",
                    "k_min": k_min,
                    "k_max": k_max,
                    "dps": 220,
                    "abs_tol": "1e-140",
                    "rows": k_max - k_min + 1,
                }
                for key, expected in coeff_boundary_checks.items():
                    actual = coeff_boundary_summary.get(key)
                    if actual != expected:
                        raise AssertionError(
                            f"{coeff_boundary_path.name}: {key} {actual!r} != {expected!r}"
                        )
                if coeff_boundary_summary.get("lambdas") != [COEFF_LAMBDA_TEXT[spec.lam]]:
                    raise AssertionError(
                        f"{coeff_boundary_path.name}: unexpected lambdas "
                        f"{coeff_boundary_summary.get('lambdas')!r}"
                    )

        log_summary_path = (
            results_dir
            / (
                f"arb_edrei_log_sign_lam{spec.safe_lam}_n1_n{log_max_n}_"
                f"dps{spec.source_log_dps}{spec.source_log_suffix}_summary.json"
            )
        )
        log_summary = load_json(log_summary_path)
        log_checks = {
            "kind": "arb_edrei_log_sign_probe_summary",
            "lam": spec.lam,
            "max_n": log_max_n,
            "dps": spec.source_log_dps,
            "rows": log_max_n,
            "failed_or_inconclusive": 0,
            "all_ok": True,
        }
        for key, expected in log_checks.items():
            actual = log_summary.get(key)
            if actual != expected:
                raise AssertionError(
                    f"{log_summary_path.name}: {key} {actual!r} != {expected!r}"
                )

    summary_path = results_dir / f"{spec.stem}_summary.json"
    rows_path = results_dir / f"{spec.stem}.jsonl"
    summary = load_json(summary_path)
    checks = {
        "kind": "arb_edrei_power_hankel_probe_summary",
        "lam": spec.lam,
        "rows": spec.rows,
        "ok": spec.rows,
        "failed_or_inconclusive": 0,
        "m_min": spec.m_min,
        "m_max": spec.m_max,
        "shift_min": spec.shift_min,
        "shift_max": spec.shift_max,
        "needed_max_n": spec.needed_max_n,
        "dps": spec.dps,
        "all_ok": True,
    }
    for key, expected in checks.items():
        actual = summary.get(key)
        if actual != expected:
            raise AssertionError(f"{summary_path.name}: {key} {actual!r} != {expected!r}")
    if not bool(summary.get("rigorous")):
        raise AssertionError(f"{summary_path.name}: rigorous flag is false")
    counts = summary.get("counts", {})
    if counts.get("positive") != spec.rows or len(counts) != 1:
        raise AssertionError(f"{summary_path.name}: unexpected counts {counts!r}")

    if not rows_path.exists():
        raise FileNotFoundError(f"missing row log {rows_path.name}")
    row_count = 0
    with rows_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            row = json.loads(line)
            row_count += 1
            if row.get("kind") != "arb_edrei_power_hankel_probe":
                raise AssertionError(f"{rows_path.name}: bad row kind")
            if row.get("lam") != spec.lam:
                raise AssertionError(f"{rows_path.name}: lambda mismatch")
            if row.get("classification") != "positive":
                raise AssertionError(f"{rows_path.name}: non-positive row {row!r}")
            if row.get("contains_zero"):
                raise AssertionError(f"{rows_path.name}: row contains zero {row!r}")
            if not row.get("ok"):
                raise AssertionError(f"{rows_path.name}: non-ok row {row!r}")
    if row_count != spec.rows:
        raise AssertionError(f"{rows_path.name}: row count {row_count} != {spec.rows}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", default="work/rh_compute/results", type=Path)
    args = parser.parse_args()

    for spec in SPECS:
        validate_spec(args.results_dir, spec)
        print(
            "OK Edrei power-Hankel diagnostic: "
            f"lambda={spec.lam}, m={spec.m_min}..{spec.m_max}, "
            f"shifts={spec.shift_min}..{spec.shift_max}"
        )
    print(f"validated {sum(spec.rows for spec in SPECS)} finite Edrei power-Hankel diagnostics")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
