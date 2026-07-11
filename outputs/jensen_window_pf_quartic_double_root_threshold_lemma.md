# Jensen-Window PF Quartic Double-Root Threshold Lemma

Date: 2026-07-10

Status: exact local quartic heat threshold with the global quartic
invariant open. This is not a proof of degree-4 zeta hyperbolicity,
PF-infinity, RH, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_quartic_double_root_threshold_lemma.json
python work/rh_compute/scripts/jensen_window_pf_quartic_double_root_threshold_lemma.py
python work/rh_compute/scripts/check_jensen_window_pf_quartic_double_root_threshold_lemma.py
```

Current result:

```text
validated Jensen-window PF quartic double-root threshold lemma: 10 rows, 0 issues, 5 exact coordinate identities, 1 double-root splitting criterion, 1 branch-aware inward threshold, 1 triple-root equality, 1 tangent factor, 1 explained countermodel, 1 open global-invariant handoff
```

## Double-Root Coordinates

Write a nondegenerate hyperbolic quartic boundary point as

```text
P(w)=(1+a*w)^2*(1+b*w)*(1+c*w),
2*a+b+c=4, p=b*c.
```

Set `A=-3*a^2+8*a+p` and `B=-a^2+2*a+p`. Then

```text
x=(-3*a**2 + 8*a + p)/6
y=18*a*(-a**2 + 2*a + p)/(-3*a**2 + 8*a + p)**2
z=2*p*(-3*a**2 + 8*a + p)/(3*(-a**2 + 2*a + p)**2)
```

Both `A=6*x` and `B=e_3/(2*a)` are positive on the positive-root
boundary stratum.

## Heat Threshold

If `H=P_lambda`, local double-root splitting is real exactly when

```text
H(-1/a)*P''(-1/a)<=0.
```

Exact substitution into the coefficient heat flow gives

```text
U(a,p)=(-a**2 + 2*a + p)*(3*a**2 - 5*a + 5*p)/(6*p**2)
H(-1/a)/r_k=24*(2*k+7)*p^2*(u-U)/(a*A*B).
P''(-1/a)=2*(a-b)*(a-c).
```

Therefore the exact inward condition is

```text
(3*a^2-4*a+p)*(u-U(a,p))<=0
```

If the double root lies outside the two simple roots this requires
`u<=U`; if it lies between them it requires `u>=U`.

## Triple Root

When `3*a^2-4*a+p=0`, one simple root merges with the double root.
A nonzero constant heat perturbation would split the cubic contact into
one real root and a complex pair, so first-order viability requires
`u=U`. At that equality the complete first variation factors as

```text
24*w**2*(a - 1)**2*(a*w + 1)**2*(2*a*k + a - 2*k + 1)/(a - 2)
```

and retains `(1+a*w)^2`. This is tangency, not yet a higher-time
invariance proof.

## Countermodel And Handoff

For the earlier rational obstruction,

```text
curvature=-3749/10000<0
u=104061328/112253067
U=57204145/61133184
U-u=2212607757569/254162496276864>0.
```

The middle-root branch requires `u>=U`, explaining the outward heat
direction exactly. The next task is to express and propagate this
branch-aware threshold as a closed contraction-coordinate invariant,
including the tangent triple-root stratum. Higher degrees remain open.
