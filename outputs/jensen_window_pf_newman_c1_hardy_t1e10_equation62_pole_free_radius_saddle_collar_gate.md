# Pole-free radius and five-saddle collar gate

Date: 2026-08-09

Status: exact diagnostic contour-deformation target isolated; not a proof of RH

Lewis (2015), equation (8), uses a circle centred at `t/pi`.  Its left and right real
intercepts are `L=t/pi-R` and `U=2t/pi-L`.  The paper requires the radius to be
adjusted whenever an intercept would hit an odd-integer pole.

At the central output this gate chooses two explicit admissible left
intercepts inside the same gap `(159777,159779)`:

```text
source-tail selector L_s = [159777.08758516150632992061230004331586097236677835194327254377231753769867814640534724306305606880641878319745 +/- 8.92e-106]
diagnostic selector  L_c = [159778.12477143756892831476566535033374263925870443309601192763982602477765211732905705012005715648437862401055 +/- 1.65e-105]
```

The equation-(124) root relation gives tail starts `37946` and `37941`,
respectively.  The first is the published/source endpoint convention; the
second is the introduced diagnostic midpoint convention.  On the far
side both corresponding intercepts lie strictly between the same consecutive
odd poles `6366037945` and `6366037947`.  Both radii are also strictly
smaller than the radial distance to the branch point `z=+a`, with `z=-a`
farther away.  Hence the radial annulus contains neither an odd pole nor either
denominator branch point, and the enclosed alpha-pole roster is unchanged.

Appendix-B equation (B9) has lower limit

```text
u_0(R)=t/(pi R),
```

while its saddle is exactly, by (B15a),

```text
u_sad(N)=1/(1-1/N-2pi N/t).
```

The identity

```text
1/u_sad(N)=1-pi[2N+t/(pi N)]/t
```

shows that saddle crossing is governed by the same equation-(124) threshold.
At all fifteen outputs, `N=37940` remains inside both `D` ranges,
`N=37946` remains outside both, and exactly

```text
N=37941,37942,37943,37944,37945
```

cross the moving lower limit.  Therefore

```text
D(N,R_s)-D(N,R_c)
 = integral_[u_0(R_s)]^[u_0(R_c)] f_D(N,u) du
```

is an exact finite-interval identity for each collar mode, before saddle
approximation.

This does not permit adding the five classical main terms while keeping the
old error unchanged.  Because the two circles enclose the same poles, Cauchy
deformation makes the complete contour invariant.  If `A(R)` denotes the
first theta integral in equation (13), and

```text
S(R)=C(0,R)+D(0,R)+2 sum_(N>=1)(-1)^N[C(N,R)+D(N,R)],
```

then equations (13) and (24) give the exact compensation law

```text
R_s A(R_s)-i t S(R_s)/pi
 =R_c A(R_c)-i t S(R_c)/pi.
```

The diagnostic midpoint is therefore a mathematically defined finite
re-splitting candidate, not a paper-prescribed or proved sharper approximation.
Its transported `C`, `D`, first-integral, endpoint, and saddle remainders would
have to be derived as a new theorem.  This gate proves only the exact geometry;
it supplies no height-uniform remainder and no RH implication.
