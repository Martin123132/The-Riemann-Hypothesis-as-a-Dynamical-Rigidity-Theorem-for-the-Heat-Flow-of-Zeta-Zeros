# Jensen-Window PF Quartic-Quintic Polar-Contact Lemma

Date: 2026-07-10

Status: exact adjacent-degree contact theorem with the infinite-degree
closure open. This is not a proof of PF-infinity, RH, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_quartic_quintic_polar_contact_lemma.json
python work/rh_compute/scripts/jensen_window_pf_quartic_quintic_polar_contact_lemma.py
python work/rh_compute/scripts/check_jensen_window_pf_quartic_quintic_polar_contact_lemma.py
```

Current result:

```text
validated Jensen-window PF quartic-quintic polar-contact lemma: 10 rows, 0 issues, 4 exact polar identities, 1 strict nonroot test, 1 multiplicity rule, 1 double-to-triple theorem, 1 quintic contact factorization, 1 cofactor gate, 1 open all-degree handoff
```

## Polar Identity

The normalized adjacent windows satisfy

```text
P_4=P_5-(w/5)*P_5'
```

If `P_5=product_i(1+alpha_i*w)` with every `alpha_i>0`, then

```text
(5*P_5-w*P_5')/P_5=sum_i 1/(1+alpha_i*w),
d/dw=-sum_i alpha_i/(1+alpha_i*w)^2<0.
```

Thus a polar-derivative root away from the roots of `P_5` is simple.
At a nonzero quintic root, the polar derivative lowers multiplicity by
one. Consequently a double root of `P_4` forces a triple root of any
hyperbolic positive-root quintic extension.

## Threshold Equals Triple Contact

In the quartic double-root coordinates of the threshold lemma,

```text
P_5(-1/a)=-(3*a**4 - 11*a**3 + 2*a**2*p + 10*a**2 - 5*a*p + 6*p**2*u - 5*p**2)/(3*a**2*(-a**2 + 2*a + p))
```

so the shared root condition is exactly `u=U(a,p)`. At this value,

```text
(a*w + 1)**3*(3*a**2*w**2 - 5*a*w**2 - 9*a*w + 5*p*w**2 + 15*w + 3)/3
```

The remaining quadratic factor has discriminant

```text
-5*(-3*a**2 + 14*a + 4*p - 15)/3
```

which is the explicit residual gate for a fully hyperbolic quintic
contact.

## General Handoff

The binomial identity extends the polar relation to every degree:

```text
P_d=P_(d+1)-(w/(d+1))*P_(d+1)'
```

This exposes a plausible higher-minor hierarchy: degree `d+1` controls
boundary multiplicity for degree `d`. It is not yet a proof because no
noncircular infinite-degree closure or uniform compactness theorem has
been supplied for the zeta coefficients.
