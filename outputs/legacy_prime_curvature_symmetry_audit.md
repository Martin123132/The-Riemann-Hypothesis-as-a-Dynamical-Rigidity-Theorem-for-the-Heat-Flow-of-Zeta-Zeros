# Legacy Prime-Curvature Symmetry Audit

Date: 2026-07-28

Status: exploratory source reconstruction, artifact control, and
finite spectral audit. This is not a proof of an Xi contact gap,
`Lambda<=0`, PF-infinity, RH, or a Clay-prize result.

```text
work/rh_compute/results/legacy_prime_curvature_symmetry_audit.json
python work/rh_compute/scripts/legacy_prime_curvature_symmetry_audit.py
python work/rh_compute/scripts/check_legacy_prime_curvature_symmetry_audit.py
outputs/legacy_prime_curvature_symmetry_audit_controls.png
```

## Scope And Provenance

The complete plot source was found in the local read-only legacy
Motion-TimeSpace archive. GitHub was not needed, and no external project
file was edited. The eight supplied JPEGs contain 7
unique files; Photos 1 and 2 are byte-identical.

The source uses `N=2000`, every seventh prime starting at 3, and
`threshold=500`. That gives 44 sampled
primes, 1936 ordered pairs, and
488 retained pairs.

## Exact Geometry

The plotted scalar is

```text
E(p,q)=p^2+q^2-floor(sqrt(p^2+q^2))^2.
```

Writing `rho=sqrt(p^2+q^2)` and `m=floor(rho)` gives

```text
E=(rho-m)(rho+m).
```

Thus the raw continuous formula is radial, and `E<=T` selects thin
annular pieces next to the circles `p^2+q^2=m^2`. The visible arcs are
therefore expected before any prime theorem. Swapping `p` and `q`
leaves the formula unchanged, so the exact 44 by 44 defect matrix and
the threshold mask are transpose-symmetric by construction.

The integer `r=m` is not required to be prime. Exact equality is
impossible for odd `p,q`: `p^2+q^2=2 mod 4`, whereas a square is 0 or
1 mod 4. The observed defects occupy only residue classes
[1, 2, 6] modulo 8.

There is genuine arithmetic beneath the picture. On the diagonal,
`E(p,p)=1` is exactly the negative Pell equation

```text
r^2-2p^2=-1.
```

The prime hits below 2000 are `[{'p': 5, 'r': 7}, {'p': 29, 'r': 41}]`,
but the every-seventh-prime source sample misses both. This is a real
quadratic-form thread, not evidence about zeta zeros.

## Rendering Controls

Only points with `E<=500` are passed to cubic `griddata`; values above
the threshold are omitted rather than plotted. Cubic triangulation
breaks the exact transpose symmetry:

```text
raw defect max |E-E^T| = 0
cubic RMS |Z-Z^T|      = 9.18185036
cubic max |Z-Z^T|      = 246.182729
smoothed relative RMS  = 0.0185659845
nearest max |Z-Z^T|    = 0
```

The lower-left density has a concrete boundary cause. Of the
`E<=10` ordered hits, 86
are the `p=3` or `q=3` axes with the deterministic value `E=9`,
a share of 76.106%.

The named "topological skeleton" is the upper 15 percent of a Sobel
gradient magnitude followed by two binary erosions. It is not
skeletonization and computes no topological invariant. The divergence
is the negative discrete Laplacian of the interpolated scalar. The
curl is 1.67e-16 at worst, as expected
when commuting the same roll-based finite differences; `np.roll`
also imposes an unphysical periodic boundary.

## Arithmetic Null

The exact radial and transpose symmetries survive count-matched odd
inputs. A more conservative finite null chooses one prime from each
successive block of seven primes:

| defect threshold | source fraction | null mean | z score |
|---:|---:|---:|---:|
| 10 | 0.058368 | 0.016814 | 2.544 |
| 25 | 0.070248 | 0.030821 | 1.916 |
| 50 | 0.080062 | 0.052725 | 1.238 |
| 100 | 0.098657 | 0.071734 | 1.221 |
| 200 | 0.136364 | 0.128387 | 0.403 |
| 300 | 0.174070 | 0.179514 | -0.323 |
| 500 | 0.252066 | 0.259677 | -0.646 |

These are selected-sample, jointly inspected finite diagnostics. They
do not establish a prime anomaly. In particular, the source pictures
use threshold 500, where the shell geometry dominates.

## Spectral Sign

Every tested raw kernel is real symmetric, but all are indefinite:

| kernel | negative | zero | positive | min eigenvalue | max eigenvalue |
|---|---:|---:|---:|---:|---:|
| defect | 22 | 0 | 22 | -13980.4 | 67907.1 |
| threshold adjacency | 20 | 0 | 24 | -4.61231 | 14.5407 |
| closeness | 21 | 0 | 23 | -2119.64 | 5139.76 |
| exp minus defect over 100 | 20 | 0 | 24 | -4.54758 | 7.84733 |

For example, the defect matrix has
22 negative and
22 positive eigenvalues; the
nonnegative closeness matrix still has
21 negative eigenvalues.
Centering away the constant vector does not fix the sign.

The graph Laplacian `diag(W 1)-W` is positive semidefinite, with
0 negative eigenvalues and
1 zero mode. That positivity is true
for every nonnegative symmetric graph and is created by the
transformation. It becomes relevant to RH only if the actual Xi
contact scalar is proved equal to, or coercively bounded by, that
energy with every signed remainder retained.

## Zero Overlay

The source solves

```text
zeta(1/2+i*t)=0.
```

Here the zero height is `Re(t)`, ranging from
14.1347251417 to 77.1448400689.
The code instead extracts `Im(t)`. At 40 decimal digits those values
span only 7.87e-48;
they are numerical root residuals. Min-max scaling magnifies them by
about 2.54e+50, and changing the
working precision moves a scaled coordinate by as much as
0.816.

Finally, the source sets `rz_x=scaled` and `rz_y=scaled`. The cyan
diagonal is therefore imposed twice: by assigning the same coordinate
to both axes, and by stretching numerical residuals rather than zero
heights. The overlay contains no measured zero/curvature alignment.

## Old Invariance Guard

For any holomorphic `f(s)`,

```text
(partial_sigma^2+partial_t^2)f=0
```

throughout its analytic domain. The legacy note correctly invokes this
at one point, but later asserts that an off-line zero would have a
nonzero Laplacian. That step is incompatible with analyticity. The
proposed two-dimensional Laplacian flow is trivial on holomorphic xi;
similarly, multiplying critical-line xi by a nonzero scalar function
preserves its zeros tautologically.

This is the same logical warning exposed by the pictures: an invariant
can be exact because it was built into the map, without localizing the
unknown zeros.

## Zeta-Aware Replacement

The defensible part of the old idea is a radial prime heat correlation.
With the von Mangoldt weight, define

```text
A(y)=sum Lambda(n) exp(-pi*n^2*y).
```

For `Re(s)>1/2`, termwise Mellin transformation gives

```text
integral A(y)y^(s-1)dy
 =Gamma(s)pi^(-s)[-zeta'(2s)/zeta(2s)].
```

Also

```text
A(y)^2
 =sum Lambda(m)Lambda(n)exp(-pi*y*(m^2+n^2)).
```

This is an exact zeta-aware version of radial prime-pair geometry. The
`pi` is from the standard Gaussian Fourier/theta convention; replacing
it by any `a>0` simply gives `a^(-s)`. It is not inferred from a circle
in the images.

The open thought experiment is precise: construct a completed,
regularized version of this radial kernel and ask whether its
reflection form is exactly a known Weil positivity form or exactly the
current Xi contact scalar. Prime restriction destroys the ordinary
lattice theta functional equation, so positivity cannot be assumed.

## What It Says About The Wall

The current corpus has already converted the useful geometric
intuition into exact mathematics:

- an explicit positive sine-feature kernel `P_(K,R)`;
- a positive edge-feature Gram kernel `G_K`;
- a curvature-weighted energy;
- exact self-adjoint symmetrization.

The diagonal is controlled. The obstruction that remains is the signed
Mobius-labelled off-diagonal correlation; the symmetrized threshold
operator itself has an exact indefinite 2 by 2 minor. On the active Xi
route, the analogous missing result is the contact-conditioned signed
lower bound for `mathcal_C_N`, followed by a one-turn phase budget.

So the old symmetry does explain the repeated wall, but not by giving
the answer: we already have symmetry and even a Gram representation.
What we do not have is a theorem forcing the actual arithmetic
coefficient vector to avoid the negative sector.

## Route Decision

Keep the active endpoint-complete Abel/contact route primary. Retain
the legacy prime-curvature field as a small falsification laboratory
for candidate symmetry-to-energy transformations. Pursue the
von-Mangoldt radial heat bridge only if it yields an exact completed
identity with the actual Xi/contact quantity; otherwise it remains a
separate quadratic-form investigation.

## Boundary

The source reconstruction, radial factorization, modular residues,
transpose identity, interpolation controls, root-coordinate bug,
finite spectra, Pell reduction, and Mellin bridge are exact or
reproducibly finite as stated. No zero/curvature correlation,
reflection-positive prime kernel, Xi contact gap, signed Mobius gain,
one-turn phase theorem, `Lambda<=0`, PF-infinity, RH proof, or
Clay-prize conclusion is obtained.
