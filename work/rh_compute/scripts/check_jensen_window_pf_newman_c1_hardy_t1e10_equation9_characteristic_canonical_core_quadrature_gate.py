#!/usr/bin/env python3
"""Independently check the characteristic canonical-core quadrature gate."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))

from flint import arb, acb, ctx


RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_characteristic_canonical_core_quadrature_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_characteristic_canonical_core_quadrature_gate.md"
FORMULA_VERSION = "characteristic_core_degree17_tanh_v1"
T = 10_000_000_000
C = 159577
MODE_LO = 39853
MODE_HI = 39936
Y_BOUND = 64
PANELS = 32

COEFFICIENTS = {
    1: Fraction(1), 3: Fraction(-1, 3), 5: Fraction(2, 15),
    7: Fraction(-17, 315), 9: Fraction(62, 2835),
    11: Fraction(-1382, 155925), 13: Fraction(21844, 6081075),
    15: Fraction(-929569, 638512875), 17: Fraction(6404582, 10854718875),
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_complex(record: dict[str, str]) -> acb:
    return acb(arb(record["real_ball"]), arb(record["imag_ball"]))


def overlap(left: acb, right: acb) -> bool:
    return (left.real - right.real).contains(arb(0)) and (left.imag - right.imag).contains(arb(0))


class IndependentCore:
    def __init__(self) -> None:
        self.i = acb(0, 1); self.pi = arb.pi(); self.t = arb(T); self.c = arb(C); self.y = arb(Y_BOUND)
        self.eta = self.pi*self.c**2/(8*self.t); self.beta = (self.t*self.eta)**(arb(1)/3)
        self.lam = (self.eta-1)*self.t/self.beta; self.sigma = 4*self.beta/(self.pi*self.c)
        self.h = 4*self.beta/self.c; self.cy = self.sigma/self.c
        self.ds = [self.h*(arb(m)-self.c/4) for m in range(MODE_LO,MODE_HI+1)]
        self.k = (self.i*self.pi/4).exp()*(self.pi/2).sqrt(); self.r = (-self.i*self.pi/4).exp()/arb(2).sqrt()

    def f0(self,q): return self.k*(self.r*q).erf()

    def h0(self,a,p,analytic):
        root=(2*a).sqrt(); q0=-p*root/(2*a); q1=(self.y-p/(2*a))*root
        return (-self.i*p*p/(4*a)).exp()/root*(self.f0(q1)-self.f0(q0))

    def h1(self,a,p,h0):
        return p*h0/(2*a)-self.i*((self.i*(a*self.y**2-p*self.y)).exp()-1)/(2*a)

    def poly(self,z):
        return sum((arb(q.numerator)/q.denominator*z**k for k,q in COEFFICIENTS.items()),acb(0))

    def remainder_poly(self,z):
        return sum((-arb(q.numerator)/q.denominator*z**k for k,q in COEFFICIENTS.items() if k>=3),acb(0))

    def canonical(self,z,analytic):
        a=1/(4*self.beta); inner=sum((self.h0(a,z+d,analytic) for d in self.ds),acb(0))
        return (self.i*(z**3/3-self.lam*z)).exp()*inner

    def finite(self,z,analytic):
        w=z/self.beta; tau=self.poly(w); x=(1-tau)/2; a=x/(2*self.beta); common=self.beta*tau; inner=acb(0)
        for d in self.ds:
            p=common+d; h0=self.h0(a,p,analytic); inner += h0+self.cy*self.h1(a,p,h0)
        phase=-self.lam*z+self.beta**3*self.remainder_poly(w)
        return w.cosh()**(-arb(3)/2)*(self.i*phase).exp()*inner

    def airy(self,y,_):
        kernel=sum(((-self.i*d*y).exp() for d in self.ds),acb(0))
        return (-(self.lam+y)).airy_ai()*(self.i*y*y/(4*self.beta)).exp()*kernel


def main() -> None:
    require(RESULT.is_file() and NOTE.is_file(), "missing result or note")
    artifact=json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "gate not passed")
    for record in artifact["dependencies"].values():
        path=REPO_ROOT/record["path"]; require(path.is_file() and file_hash(path)==record["sha256"],f"dependency hash drift: {path}")
    for record in artifact["sources"].values():
        path=REPO_ROOT/record["path"]; require(path.is_file() and file_hash(path)==record["sha256"],f"source hash drift: {path}")

    cache_path=REPO_ROOT/artifact["dependencies"]["cache"]["path"]
    rows=[json.loads(line) for line in cache_path.read_text(encoding="utf-8").splitlines() if line.strip()]
    rows=[row for row in rows if row.get("formula_version")==FORMULA_VERSION]
    require(len(rows)==64,"cache row count drift")
    keys={(row["integral_kind"],int(row["panel_index"])) for row in rows}
    require(len(keys)==64,"duplicate or missing cache key")
    by_key={(row["integral_kind"],int(row["panel_index"])):row for row in rows}
    canonical=sum((parse_complex(by_key[("canonical",i)]["value"]) for i in range(PANELS)),acb(0))
    finite=sum((parse_complex(by_key[("finite_t_degree17",i)]["value"]) for i in range(PANELS)),acb(0))
    values=artifact["certified_values"]
    require(overlap(canonical,parse_complex(values["canonical_rectangle_ball"])),"canonical sum drift")
    require(overlap(finite,parse_complex(values["finite_t_degree17_rectangle_ball"])),"finite-t sum drift")

    ctx.dps=40
    core=IndependentCore()
    sqrt_lam=core.lam.sqrt(); plus=sqrt_lam-core.ds[-1]; minus=sqrt_lam+core.ds[0]
    require(plus>arb("0.0058") and minus>arb("0.032"),"independent stationary margins failed")

    full=acb(0)
    for left in range(0,Y_BOUND,4):
        full += acb.integral(core.airy,arb(left),arb(left+4),abs_tol=arb("1e-20"),rel_tol=arb("1e-20"),eval_limit=200000,depth_limit=40)
    full*=2*core.pi
    require(overlap(full,parse_complex(values["full_z_line_compact_y_ball"])),"Airy-first value drift")

    for kind,index,function in (("canonical",0,core.canonical),("finite_t_degree17",31,core.finite)):
        left=arb(-4)+arb(8*index)/PANELS; right=arb(-4)+arb(8*(index+1))/PANELS
        fresh=acb.integral(function,left,right,abs_tol=arb("1e-11"),rel_tol=arb("1e-11"),eval_limit=500000,depth_limit=50)
        cached=parse_complex(by_key[(kind,index)]["value"])
        require(overlap(fresh,cached),f"independent panel drift: {kind} {index}")

    correction=arb(values["finite_t_relative_correction_ball"])
    ztail=arb(values["canonical_abs_z_gt_4_relative_to_rectangle_ball"])
    require(correction<arb("6.6e-6"),"relative correction threshold failed")
    require(ztail<arb("0.064"),"z-tail threshold failed")
    require(artifact["decision"]["canonical_z_tail_negligible"] is False,"tail guard drift")
    require(artifact["decision"]["outer_normal_saddle_join_proved"] is False,"proof boundary drift")
    note=NOTE.read_text(encoding="utf-8")
    for token in ("I_M(infinity,Y)","<6.6e-6","<0.064","not a proof"):
        require(token in note,f"note token missing: {token}")
    print("validated characteristic canonical core independently: 64 cache rows, correction<6.6e-6, z-tail<0.064")


if __name__ == "__main__":
    main()
