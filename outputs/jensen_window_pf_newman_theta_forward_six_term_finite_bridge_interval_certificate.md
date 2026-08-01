# Newman Theta Forward Six-Term Finite-Bridge Interval Certificate

Date: 2026-07-25

Status: complete rigorous finite Q207 theorem.
This is not a cofinal theorem and not a proof of `Lambda<=0` or RH.

## Domain

```text
low-time left strip: [1/1035,1/155] x [38,69]
full-time right strip: [1/1035,1/5] x [69,245]
retained ordinary theta terms: n=1,...,6
```

The completed Q31 theorem supplies the complementary
`[1/155,1/4] x [0,69]` region, while the compact theorem
supplies `x<=38`. Together with the published `Lambda<=1/5`
boundary input, a complete cache promotes the finite rectangle

```text
Q_207=[1/1035,1/4] x [0,245].
```

## Shared Models

Each half-unit x panel uses one 352-bit model centered at
`t=(1/1035+1/5)/2`, shared by every time cell:

```text
time Taylor order: 30
frequency Taylor order: 24
oscillatory moments per panel: 86
maximum absolute moment: 86
cutoff: u=11/5 with analytic forward tail
```

The direct six-term omitted-tail theorem is used without
integration by parts:

```text
|J-J_6^F|<=16*x^4*E0(7),
|J'-(J_6^F)'|<=64*x^3*E0(7)+16*x^4*E1(7).
```

## Progress

Panels: `414/414`.
Certified panels: `414`.
Unresolved panels: `0`.
Initial time cells represented: `7102`.
Certified leaves: `7102`.
Subdivisions: `0`.
Minimum certified ratio: `22234414323684201470099456.0000000000000000000000000000000000000000000000000`.

## Proof Boundary

Completion proves only the finite Q207 diagonal. The surviving
Newman obligation is retained adaptive first-jet separation on
every cofinal transition cell for `x>=245`.
