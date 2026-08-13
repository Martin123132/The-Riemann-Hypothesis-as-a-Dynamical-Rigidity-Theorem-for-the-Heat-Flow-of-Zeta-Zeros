# Ideal-Cubic Double-Morse Hyperbolic-Transport Gate

Date: 2026-08-03

Status: exact transport profiles, hyperbolic phase, carrier ridge, transverse kernel, and ideal derivative identities proved; signed ridge-completed ideal estimate open; not a proof of RH.

## Transport Profiles

Define h(w)=w-1-log w, mathfrak R(w)=2h(w)/(w-1), and mathfrak B(w)=2h(w)/(w-1)^2, with removable values mathfrak R(1)=0 and mathfrak B(1)=1. Also z(w)J(w)=w mathfrak R(w).

The exact logarithmic Jacobian identity is zJ partial_w log J=1-mathfrak B(w). It converts the weighted flux transport into two elementary profile functions.

At w=1+d, mathfrak R=d-2d^2/3+d^3/2-2d^4/5+..., while mathfrak B=1-2d/3+d^2/2-2d^3/5+d^4/3+.... Thus the transport has a removable first-order zero at the double saddle.

mathfrak B(w)>0 for w>0, while mathfrak R(w) has the sign of w-1. A coefficient-blind ridge transport therefore changes orientation across the reciprocal saddle.

## Exact Interior Transport

The boundary-flux transport vector is mathcal L=eta partial_eta-y partial_y=z(rho)J(rho)partial_rho-z(v)J(v)partial_v.

For G_P=exp(S(x))P(x)J(rho)J(v)/rho, x=log(pv/rho), the exact identity is mathcal L G_P=exp(S)J(rho)J(v)/rho times {[mathfrak B(v)-mathfrak B(rho)-mathfrak R(rho)]P-[mathfrak R(rho)+mathfrak R(v)](P'+gP)}.

At rho=v=1 both coefficients vanish. To first order they are -(delta_rho+2delta_v)/3 and -(delta_rho+delta_v), so the interior transport gains one local displacement but not two.

On the carrier ridge v=rho, mathcal L G_P=-exp(S(log p))J(rho)^2 mathfrak R(rho){P+2(P'+gP)}(log p)/rho. The phase is constant there, so this trace must be recombined with the returned carrier rather than bounded as an oscillatory cell error.

## Hyperbolic Coordinates

Put s=(eta+y)/sqrt(2) and d=(eta-y)/sqrt(2). Then (eta^2-y^2)/2=sd and deta dy has unit absolute Jacobian.

The transport becomes mathcal L=d partial_s+s partial_d. It differentiates transversely on the carrier ridge d=0, where the phase e(sd) is constant.

For fixed s!=0 and finite transverse limits d_-<d_+, integral_(d_-)^(d_+)e(sd)dd=kappa{e(sd_+)-e(sd_-)}/s, with removable value d_+-d_- at s=0. This is the exact finite transverse Fourier kernel.

Subtract the ridge amplitude before estimating the transverse residual. For |s| away from zero, integrate the residual in d with the exact 1/s kernel; keep a central s-neighborhood and reciprocal side traces joined to the returned carrier. A two-dimensional modulus would erase this structure.

## Ideal Cubics

The exact ideal derivatives are (P_H^0)'=-i(x-u_N)(3x+u_N)/4 and (P_T^0)'=-i(x+u_N)(3x-u_N)/4+u_(N,x). They enter the transport only through P'+gP.

At x=-u_N, P_H^0=0 but (P_H^0)'=-iu_N^2. Also P_T^0=-2u_Nu_(N,x) and (P_T^0)'+g_NP_T^0=u_(N,x)(1-2g_Nu_N). Thus the Hermitian value null does not delete the terminal-neighborhood interior transport.

Neither P_H^0+2{(P_H^0)'+gP_H^0} nor its transpose analogue vanishes identically. The nonoscillatory ridge transport is therefore a live joined-carrier term, not a polynomial null.

## Handoff

Use the exact hyperbolic transverse kernel after subtracting and recombining the carrier ridge. Do not estimate eta G_eta-y G_y as a two-dimensional absolute remainder and do not infer that the Hermitian terminal value null removes its derivative transport.

Write the complete ideal Hermitian and transpose joins in (s,d) coordinates, isolate the ridge amplitude G_P(s,0), and prove the exact cancellation or retained coefficient supplied by both returned carriers and the transpose terminal block. Then split the residual into a central |s| chart and an outer transverse 1/s kernel, keeping physical rails, reciprocal sides, and far compression joined.

## Pi Provenance

The finite transverse kernel uses only e(x)=exp(2pi i x) and kappa=1/(2pi i). The profile functions come from w-1-log w and the exact Morse Jacobian; no geometric or fitted constant is inserted.

## Proof Boundary

No signed ridge-completed ideal join, transverse residual bound, far aggregate, h^2 reserve, completed current, or RH-level theorem is proved. This gate proves no signed completed-current bound, `Phi_B` theorem, contact exclusion, retained aggregate or Xi theorem, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion.
