# Jensen-Window PF Newman Signed Universal-Factor Residual Gate

Date: 2026-07-11

Status: exact and interval-certified signed universal-factor residual
obstruction. This is not a proof or disproof of RH or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_newman_signed_universal_factor_residual_gate`.

```text
work/rh_compute/results/jensen_window_pf_newman_signed_universal_factor_residual_gate.json
python work/rh_compute/scripts/jensen_window_pf_newman_signed_universal_factor_residual_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_signed_universal_factor_residual_gate.py
```

Current result:

```text
validated Jensen-window PF Newman signed universal-factor residual gate: 8 rows, 0 issues, 1 exact multiplier reduction, 1 rational parameter rectangle, 2 discriminant exclusions, 2 Acb spectral certificates, 4094 adaptive leaves, 3416 base boxes, maximum depth 6, 1 exhaustive signed universal-factor obstruction, 1 coupled handoff
```

## Multiplier Cone

```text
U_(beta,gamma)(z)=cosh(9z)+beta*cosh(5z)+gamma*cosh(z)
U_(beta,gamma)(z)=cosh(z)*Q_(beta,gamma)(cosh(z)^2)
Q_(beta,gamma)(y)=5*beta + gamma + 256*y**4 - 576*y**3 + y**2*(16*beta + 432) - y*(20*beta + 120) + 9
U_(beta,gamma) has only imaginary zeros iff every zero of Q_(beta,gamma) is real and lies in [0,1].
```

Thus the standard Polya universal-factor hypothesis is a concrete
quartic root-location condition, including signed coefficients.

## Necessary Rectangle

```text
If Q has roots in [0,1], Q'(1)=40+12*beta>=0, so beta>=-10/3>-17/5.
Rolle gives three real roots of Q', hence f(beta)<=0 for f=32beta^3-297beta^2+1053beta-1215. Since disc(f')=-51516, f is strictly increasing, and f(11/5)=607/125>0, beta<11/5.
Q(1)=1+beta+gamma=256*product_j(1-r_j)>=0, so s=beta+gamma>=-1.
Global residual nonnegativity at u=0 and the Arb origin bound give s<51/10.
-17/5<=beta<11/5, -1<=s=beta+gamma<51/10
```

The upper sum bound is the independently certified Xi-kernel origin
inequality from the classical three-block gate.

## Adaptive Certificate

The residual Laguerre expression is

```text
L[R_(beta,gamma)]=L[E]-beta*B[E,P5]-gamma*B[E,P1]+beta^2*L[P5]+beta*gamma*B[P5,P1]+gamma^2*L[P1]
```

Central differences use h=1e-6. Acb boxes of radius 0.1 bound each component on Cauchy circles of radius 0.04; the added f' and f'' errors are h^2*M/r^3 and 2*h^2*M/r^4.
The adaptive classification is:

```text
initial boxes=3416
adaptive leaves=4094
maximum depth=6
L(x=86) leaves=2329
L(x=122) leaves=1281
critical-discriminant leaves=240
quartic-discriminant leaves=244
unresolved boxes=0
leaf sha256=fbc1f39873be290c8b2b976396d5cf778d8604c008bfa053fa397f9e4c57b69e
```

f(beta)>0 makes disc(Q')<0, while Delta(beta,gamma)<0 makes disc(Q)<0; either condition excludes four roots in [0,1].
Every leaf is therefore either spectrally incompatible with an
independently Laguerre-Polya residual or algebraically outside the
standard imaginary-zero multiplier cone.

## Boundary

This exhausts signed 9/5/1 bases certified by the standard Polya universal-factor hypothesis when their residual kernel is globally nonnegative and is treated as an independent Laguerre-Polya block. It does not reject higher shifts or a coupled signed identity.
The next route must add higher shifts or derive one genuinely coupled
matrix square that retains and controls the mixed signs.

References: https://arxiv.org/abs/0801.2996,
https://arxiv.org/abs/1502.06844, and
https://doi.org/10.1007/BF02565336.
