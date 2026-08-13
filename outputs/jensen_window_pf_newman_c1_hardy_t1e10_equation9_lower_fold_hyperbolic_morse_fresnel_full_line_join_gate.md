# Full-line hyperbolic Morse--Fresnel prototype join

Date: 2026-08-12

Status: full-line one-mode prototype join certified at three heights; this is
not a proof of 399-event propagation or a complete `T_upper` estimate

Return both the exact retained-Fresnel logistic chart and its Airy model to
the same physical variables `(z,y)`.  Their full real domains then coincide:

```text
Delta I(t)=int_R int_0^64 [A_t(z,y) exp(i Theta_t(z,y))
                           -exp(i Theta_0(z,y))] dy dz. (FLJ1)
```

All normal saddles lie inside `|z|<9`.  Deform each real tail to
`Im z=1`, keeping the vertical connectors at `z=+/-9`.  On every finite
segment, Arb Taylor series integrate the signed difference coefficientwise.
The inner quadratic-exponential integral and its first moment are generated
by the exact recurrence

```text
c_0=1, c_1=i b,
c_n=[i b c_(n-1)+2 i a c_(n-2)]/n,                    (FLJ2)
```

followed by exact termwise integration on `0<=y<=64`.

The signed horizontal pieces are Taylor-integrated through `|Re z|=15`.
Beyond that point, the exact phase obeys a logistic-Gaussian lower bound and
the canonical phase a Gaussian lower bound; their combined absolute tail is
recorded in the machine artifact.  Combining core, connectors, and both
tails gives

```text
lower_face: normalized=[1.6686725118120193656068295240402221679687500000000000000000000000000000000000000e-5 +/- 4.82e-7], physical=[7.5116654514317815483082085847854614257812500000000000000000000000000000000000000e-6 +/- 2.17e-7]
event: normalized=[1.6725211310131271602585911750793457031250000000000000000000000000000000000000000e-5 +/- 4.82e-7], physical=[7.5289903165298710518982261419296264648437500000000000000000000000000000000000000e-6 +/- 2.17e-7]
upper_face: normalized=[1.6763688662990716693457216024398803710937500000000000000000000000000000000000000e-5 +/- 4.82e-7], physical=[7.5463112096940676565282046794891357421875000000000000000000000000000000000000000e-6 +/- 2.17e-7]
```

and hence uniformly at the event and both `pi/16` faces,

```text
|Delta I(t)| < 0.000019,
(sqrt(2)/pi)|Delta I(t)| < 0.0000086.                 (FLJ3)
```

The `pi` normalization is inherited from the Kummer quadratic phase and
integer Fourier--Poisson character.  No fitted geometric constant is used.

## Boundary

This gate certifies one prototype crossing mode at three heights.  It does
not prove uniformity through all 399 turning events, control the remaining
ordinary modes, establish complete `T_upper`, or prove `Lambda<=0`, RH, or a
prize-level conclusion.
