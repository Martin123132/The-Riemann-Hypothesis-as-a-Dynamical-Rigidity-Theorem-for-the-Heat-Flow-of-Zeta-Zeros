# Carrier-oriented affine A-transition roster

Date: 2026-08-13

Status: rigorous 84-mode canonical tangent-wedge sum; not an exact curved-face
or transformed-amplitude estimate

The exact A-transition roster at `t=10^10` is

```text
target side: 39853..39894, 42 modes,
outer side:  39895..39936, 42 modes.                   (AT1)
```

Every mode uses the common odd-endpoint carrier and `epsilon_A=-1` from
Section 11.431.  The paired endpoint coefficient remains `+1` across the
target edge; only the positive full-line bulk projector coefficient changes.

For `delta=1+rho<0`, the null-coordinate formula is evaluated as two
lower-indented full lines minus monotone right steepest rays.  For
`delta>0`, the first lower interval is evaluated directly on the left
steepest ray and the second upper interval directly on the right ray:

```text
C=(L_left^-(A1,B1;r0)+L_right(A2,B2;r0))/(2pi i).     (AT2)
```

This is the same Abel value as the full-minus-right formula, but its rays
point away from their quadratic saddles.  The superscript minus records that
the direct left ray passes below the pole and already contains the Abel
half-jump; no additional `1/2` is inserted.  Along each ray the modulus is
`exp(-|A|y^2+Ly)` with certified `L<=0`.  With cutoff `20` the
largest omitted ray tail is

```text
[1.214114344290353921994915820649545197093632540426087078e-95 +/- 1.76e-150].               (AT3)
```

The null join remains numerically tight throughout:

```text
min r0=[0.077800061685016209168569504562413303920386863745106864729816157879 +/- 3.92e-67],
max r0=[0.077892933382048249717437668589082867798887176811042743014309990337 +/- 4.60e-67].                         (AT4)
```

After restoring `2(pi/(32t))^(1/4)`, `exp(-i*pi/8)`, the common phase, and
the A endpoint sign, summing in deterministic mode order before taking any
norm gives

```text
target-side affine sum=[-0.102609961802865567357493220639399043949170015220518246614086947480241541069 +/- 4.34e-32],
outer-side affine sum =[-0.0130986049279626132375700037405506226737415385972824762886932843094134995979 +/- 2.65e-32],
complete scalar sum   =[-0.115719253576501078980203829423472127830322948726787662335396020276444324589 +/- 6.99e-32],
complete affine sum   =[-0.115708566730828180595063224379949666622911553817800722902780231789655040667 +/- 6.99e-32],
affine-minus-scalar   =[1.06868456728983851406050435224612074113949089869394326157884867892839225788e-5 +/- 1.40e-31]. (AT5)
```

The termwise absolute affine sum and signed cancellation ratio are

```text
sum |q_A,m^aff|=[0.285165866822512975844202149552109548299974094359266733629940310882136672205 +/- 6.99e-32],
|sum q_A,m^aff|/sum |q_A,m^aff|
 =[0.40575882387370402492582836423632715219624246056767775275063361953 +/- 3.45e-31].       (AT6)
```

Equations (AT5)--(AT6) measure cancellation only inside the canonical affine
tangent-wedge A roster.  They are not a bound for the exact A block or for
`R_Dir`: the curved-face strip and transformed-amplitude remainder remain
unsummed, and the positive full-line projector jump plus the B, zero,
negative, half-current, and outer sectors must stay in the final common
assembly.

Pi provenance: every `pi` in (AT1)--(AT6) comes from equation (9), the exact
bi-Morse/null-coordinate changes, odd-endpoint parity, Gaussian contour
rotation, and the paper normalizer.  No fitted constant is introduced.

Proof boundary: rigorous canonical scalar and affine tangent-wedge values,
ray tails, and signed sums over modes 39853..39936 only.  No exact curved-face
or transformed-amplitude remainder, complete A endpoint theorem, `R_Dir`
estimate, complete `Q_K-T` or `T_upper`, all-height theorem, `Lambda<=0`,
PF-infinity, RH, or prize-level conclusion is proved.
