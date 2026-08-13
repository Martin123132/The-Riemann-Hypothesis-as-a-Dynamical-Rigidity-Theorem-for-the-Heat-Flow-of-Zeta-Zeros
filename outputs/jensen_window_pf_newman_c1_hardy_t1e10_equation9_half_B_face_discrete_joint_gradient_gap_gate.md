# Discrete joint-gradient gap on the B face

Date: 2026-08-13

Status: exact monotonic reduction plus saved-height interval certificate; not
a proof of the complete B-face estimate

After the sharp interior-saddle allocation, measure the two face directions
by

```text
G(x)=Phi_B'(x)/sqrt(H_B),
q_m(x)=(B*x-2m)/sqrt(2x).                             (DG1)
```

On `0<x<1/2`, both are strictly increasing:

```text
G'(x)=t(1-2x)/[2sqrt(H_B)x^2(1-x)^2]>0,
q_m'(x)=(Bx+2m)/[2sqrt(2)x^(3/2)]>0.                 (DG2)
```

The tangential root is `x_B^-`; the normal root is `2m/B`.  If a crossing
lies to the left of the trace root, the minimum of
`max(|G|,|q_m|)` is the unique balance `-G=q_m` between them.  If it lies to
the right, the balance is `G=-q_m`.  Monotonicity in the integer mode proves
that only the adjacent lattice modes 621 and 622 can minimize the gap.

For mode 621, Arb brackets the balance at

```text
x in [[0.0002425842252836300000000000000000000000000000000000000000000000000000000000000000000000000000000000000 +/- 1.78e-105],
      [0.0002425842252836400000000000000000000000000000000000000000000000000000000000000000000000000000000000000 +/- 2.45e-105]],
xi=[-28.06997072659913464828471057993714771522086259893364736490272005701206856600831677404953933845420740 +/- 1.57e-93],
min max(|G|,|q_621|)>[28.08111806883566024726761215342087595749432171481727969017547623049181821444946510080822214929649263 +/- 4.20e-96]>28.08. (DG3)
```

For mode 622,

```text
x in [[0.0002427575069609500000000000000000000000000000000000000000000000000000000000000000000000000000000000000 +/- 2.50e-105],
      [0.0002427575069609700000000000000000000000000000000000000000000000000000000000000000000000000000000000000 +/- 2.09e-105]],
xi=[22.41971186621432369280077992856595446405681635900997601806842899326080411806448464599902822125488217 +/- 1.57e-93],
min max(|G|,|q_622|)>[22.41260566351748647839147051566440223435497958698472410151239077765238277622543093389317848069314688 +/- 4.07e-96]>22.40. (DG4)
```

Consequently, for every positive integer mode and every B-face point,

```text
max(|G(x)|,|q_m(x)|)>22.40,
G(x)^2+q_m(x)^2>501.                                 (DG5)
```

This is the discrete gain missed by the continuum tangent model.  The B
trace saddle lies between two integer crossings, but in the natural
Gaussian/Fresnel units it is more than 22 units from simultaneous
criticality.  One may therefore use the smooth exact partition

```text
chi_T=G^2/(G^2+q_m^2),
chi_N=q_m^2/(G^2+q_m^2)                              (DG6)
```

and integrate tangentially where `G` dominates and normally where `q_m`
dominates.  No hard crossing collar or singular cotangent split is required
at the theorem level.

The scope is important: (DG5) is a face-residual statement after allocation
of the interior saddle.  For mode 622 the complete triangle still has its
ordinary interior saddle; this certificate does not erase or double-count
that target contribution.

Pi provenance: `pi` comes from the exact equation-(9) face phase and its
Fresnel normal coordinate.  No fitted constant is used.

Proof boundary: exact discrete face noncriticality and a denominator margin
for a two-direction partition only.  No derivative bounds for the partition
or amplitude, completed integration-by-parts estimate, A-fold splice,
complete paired residual, `T_upper`, `Lambda<=0`, PF-infinity, RH, or a
prize-level conclusion is established.
