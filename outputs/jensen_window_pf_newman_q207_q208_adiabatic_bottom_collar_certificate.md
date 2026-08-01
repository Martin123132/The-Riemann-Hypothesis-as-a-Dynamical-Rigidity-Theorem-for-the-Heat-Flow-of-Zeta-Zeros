# Newman Q207-Q208 Adiabatic Bottom-Collar Certificate

Date: 2026-07-25

Status: resumable partial bottom-collar computation.
This is a finite calibration, not an all-j successor theorem
and not a proof of `Lambda<=0` or RH.

## Domain

```text
time collar: [1/1040,1/1035]
x range: [0,245]
panels: half-unit
reference: Q208 bottom phase cells at t=1/1040
```

For `F=16*(1+x^4)*H`, the exact heat equation gives

```text
F_t=-16*(1+x^4)*H_xx,
F_xt=-64*x^3*H_xx-16*(1+x^4)*H_xxx.
```

Each panel proves

```text
delta*sup||(F_t,F_xt)|| < dist(0,C_k),
```

where `C_k` is the stored full Q208 bottom phase cell.
Reverse triangle therefore excludes the origin throughout
the complete collar over that panel.

## Progress

```text
panels=380/490
certified=379
failed=1
maximum transport ratio upper=[1.93177534917108851181666496975308065153941815438368873675273293032583957820968596190324917164165668688059 +/- 4.57e-105]
worst panel=['379/2', '190']
stop reason=failed_panel
```

## Consequence

No collar theorem is promoted from a partial cache.

The already certified derivative-positive strip
`[1/1040,1/5]x[245,246]`, together with the Q207 base, then
gives an independent finite no-contact decomposition of the
Q208 low rectangle.

## Proof Boundary

The complete status certifies only [1/1040,1/1035]x[0,245] by transport from the stored Q208 bottom phase cells. It calibrates one finite successor. No Q208-to-Q209 transport bound, no uniform all-j collar estimate, no all-j right-strip cone, no cofinal theorem, no Lambda<=0, and no RH proof is supplied.
