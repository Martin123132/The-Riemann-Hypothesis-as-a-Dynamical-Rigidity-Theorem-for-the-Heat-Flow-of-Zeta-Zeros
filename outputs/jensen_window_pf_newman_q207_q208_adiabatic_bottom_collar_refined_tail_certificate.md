# Newman Q207-Q208 Refined Adiabatic Bottom Collar

Date: 2026-07-25

Status: rigorous complete hybrid collar certificate.
This is one reverse-endpoint finite successor calibration, not an all-j
theorem and not a proof of `Lambda<=0` or RH.

## Cover

```text
coarse prefix: 379 half-unit panels on [0,189.5]
refined tail: 222 quarter-unit panels on [189.5,245]
time collar: [1/1040,1/1035]
jet scale: ell=1
```

Each refined derivative model is centered on its quarter
panel. Its transport bound is compared with the distance of
the enclosing parent Q208 bottom phase cell from the origin.

## Progress

```text
refined=222/222
refined certified=222
refined failed=0
combined panels=601
maximum ratio upper=[0.720143040423147365587981776315771540377482529162106590967736809795711042778094717184521177685563028967084 +/- 5.35e-106]
worst panel=['189', '379/2']
worst source=coarse_prefix
stop reason=complete
```

## Consequence

The combined cover proves the complete actual bottom collar [1/1040,1/1035]x[0,245] contact-free.

Its reference cells lie on the new Q208 bottom edge. Thus it rigorously
proves the collar and its transport scale, but it is not by itself the
old-edge hypothesis of the forward successor lemma. The separate
forward-successor certificate verifies that hypothesis from Q207 cells.
Together with Q207 and the independently certified
derivative-positive strip `[1/1040,1/5]x[245,246]`, a
complete collar result gives a second finite construction
of the Q208 low rectangle.

## Proof Boundary

The complete status combines a certified coarse prefix and quarter-unit refined tail only for [1/1040,1/1035]x[0,245]. It validates one finite successor. No Q208-to-Q209 transport budget, no uniform all-j collar estimate, no all-j right-strip cone, no cofinal theorem, no Lambda<=0, and no RH proof is supplied.
