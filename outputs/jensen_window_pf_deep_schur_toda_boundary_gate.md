# Jensen-Window PF Deep-Schur Toda And Boundary Gate

Date: 2026-07-16

Status: exact rectangular Toda coordinate, moving-tail zero-boundary
obstruction, strict-Schur/Jensen separation, and a rigorous rejection
of the proposed all-order endpoint hierarchy. The order-ten failures
reject that route, not the Xi/Phi Jensen bridge. This is not a proof
of RH or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_deep_schur_toda_boundary_gate.json
python work/rh_compute/scripts/jensen_window_pf_deep_schur_toda_boundary_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_deep_schur_toda_boundary_gate.py
```

## Rectangular Toda Coordinate

Put

```text
tau_(m,N)=s_((N^m))(h)=Q_(m,N-m+1)(-100)/A_0(-100)^m
```

Desnanot-Jacobi becomes the discrete Toda identity

```text
tau_(m+1,N)*tau_(m-1,N)=tau_(m,N)^2-tau_(m,N-1)*tau_(m,N+1)
```

Consequently, when the lower row is positive,

```text
tau_(m+1,N)>0 iff tau_(m,N)^2>tau_(m,N-1)*tau_(m,N+1), assuming tau_(m-1,N)>0
```

This is exact but not a propagation theorem. The right side is a
difference of positive terms, so ordinary subtraction-free cluster
positivity does not determine its sign. The proposed strict
width-log-concavity input is false at order ten and widths 9 through 12.

## Moving-Tail Boundary Check

A tempting translation is the normalized ordinary tail

```text
b_k^(s)=h_(s+k)/h_s for k>=0 and b_k^(s)=0 for k<0
```

But the ordinary tail resets every negative index to zero. The
shifted deep determinant retains `h_(s+k)` whenever `s+k>=0`.
The two Jacobi-Trudi matrices therefore agree only under the same
deep condition one was trying to escape:

```text
det[b_(mu_i-i+j)] = h_s^(-r)*s_(mu+(s^r))(h) only when mu_r>=r-1, so every tail Jacobi-Trudi index is nonnegative
```

The smallest exact witness is

```text
s_(0,0)(b^(1))-h_1^(-2)*s_(1,1)(h)=h_0*h_2/h_1^2>0
```

The bounded audit checks `251` shapes:
`31` deep translations agree, while all
`220` shallow shapes have a reset
mismatch in the exact rational test specialization. Thus the proposed
moving-tail PF equivalence is false, and a Littlewood-Richardson lift
from those tails is unavailable.

## Strict-Schur Jensen Counterexample

The generic bridge can be rejected even under a much stronger
antecedent. Take

```text
H(z)=exp(z/100)/((1-z)*(1-2*z))
```

This is an Edrei specialization with positive exponential parameter
`1/100` and positive denominator parameters `1,2`. More strongly,
the Plancherel component gives, for every partition `lambda`,

```text
s_lambda[X+Pl_epsilon]>=s_lambda[Pl_epsilon]=f^lambda*epsilon^|lambda|/|lambda|!>0
```

so the sequence is strictly Schur-positive for every partition.
Nevertheless its cubic
Jensen polynomial at shift zero is

```text
6000000+54180000*x+126540900*x^2+90420901*x^3
discriminant=-222484532394597/2000000000000<0
```

A real cubic with negative discriminant has one real zero and one
nonreal conjugate pair. Hence even strict full Schur positivity of
the unweighted sequence does not imply the binomially weighted Jensen
window is hyperbolic. The bridge must use additional Xi/Phi structure.

## Primary-Literature Fit

### lam_postnikov_pylyavskyy_schur_log_concavity

Classification: `formal_identity_support_only`.

Schur log-concavity supports formal rectangular difference identities, but numerical positivity after specialization still needs the relevant Schur values to be positive.

It does not supply: the zeta endpoint rectangle signs.

Source: https://arxiv.org/abs/math/0502446

### wagner_hadamard_products

Classification: `hypotheses_not_met`.

Hadamard products of totally positive Toeplitz matrices are not closed in general; Wagner proves narrower sufficient classes.

It does not supply: a conversion from deep endpoint minors to binomially weighted Jensen minors.

Source: https://doi.org/10.1016/0022-247X(92)90261-B

### angarone_kim_oh_soskin_dual_jacobi_trudi

Classification: `restricted_shape_only`.

The newer preprint leaves the general Hadamard Jacobi-Trudi conjecture open and proves a 3x2-avoiding ribbon-like case.

It does not supply: the open rectangles (N^m), which contain 3x2 blocks for m>=3,N>=2.

Source: https://arxiv.org/abs/2511.08969

### craven_csordas_fox_wright

Classification: `fixed_multiplier_family_support`.

Reciprocal-Gamma arithmetic progressions provide genuine multiplier sequences and explain the fixed-scale model.

It does not supply: closure under the Xi scale mixture or an endpoint-to-Jensen theorem.

Source: https://doi.org/10.1016/j.jmaa.2005.03.058

No direct closing theorem was found in this bounded audited set.
The newer Jacobi-Trudi Hadamard result is especially important to
state accurately: its general conjecture remains open, and the
proved 3x2-avoiding case does not supply a replacement for the
rejected rectangle hierarchy.

## Live Handoffs

The endpoint hierarchy is rejected:

```text
s_((N^10))(h)<0 for N=9,10,11,12
s_((N^m))(h)>0 for every m>=10 and N>=m-1 is false
```

The logically separate Xi/Phi bridge remains open:

```text
prove every binomially weighted Jensen window for the actual Xi/Phi sequence using hypotheses stronger than unweighted Schur positivity and weaker than all-shift signed-Hankel positivity
```

The next route cannot use bare Toda positivity, ordinary shifted-tail
PF, generic Schur positivity, or the false all-order rectangle
antecedent. It must identify a weaker Xi/Phi-specific mechanism for
the binomially weighted Jensen windows.
