# Quartic Outer-Contact Length-14 Survivor Gate

Date: 2026-07-25

Status: exact finite outer-contact length-14 survivor; not Xi,
not a uniform theorem, and not a proof of `Lambda<=0` or RH.

```text
work/rh_compute/results/jensen_window_pf_quartic_outer_contact_length14_survivor_gate.json
python work/rh_compute/scripts/jensen_window_pf_quartic_outer_contact_length14_survivor_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_quartic_outer_contact_length14_survivor_gate.py
```

## Contact

Use the exact normal-form parameters

```text
delta=13/200, q=1/2,
s=T-1/100000,
x_5-U=169/4000000000>0.
```

Thus this is a strict outward left-outer-root quartic contact.
Set

```text
y_6=y_7=...=y_14=99/100.
```

The exact corridor recurrence constructs increasing contractions
through `x_14`. All thirteen contractions are below one and above
their pointwise walls; all twelve scaled-defect, cubic, and
reciprocal-defect steps are strict.

## Signed Layers

Exact rational recurrence proves eleven positive contiguous
order-three gaps and nine positive order-four cap margins. At
512-bit Arb precision, direct arbitrary-column enumeration on
`A_1,...,A_15` gives

```text
order two:  455/455 positive
order three: 1001/1001 positive
order four: 1287/1287 positive
order five: 924/924 positive
order six:  330/330 positive
order seven: 45/45 positive
order eight: 1/1 positive
```

Every Arb lower endpoint is strictly positive and is independently
recomputed from the exact rational inputs by the checker.

## Compatibility

The terminal exact signs are

```text
Delta_14=[9.2006534331808475492048952593230341909823953283692325167253471727672354190166994e-8 +/- 2.98e-88]>0,
Delta_15=[7.8927689829553241289793474975566779179896321629356658743715968124790706817822267e-8 +/- 2.53e-88]>0.
```

Therefore the one-contact interval theorem cannot be promoted to
a uniform all-contact length-14 obstruction under only these
finite scalar and signed-Hankel hypotheses.

## Adjacent-Degree Separation

The same exact coefficients define the normalized adjacent
quintic

```text
P_5(w)=1+5*w+10*x_2*w^2+10*x_2^2*x_3*w^3
       +5*x_2^3*x_3^2*x_4*w^4
       +x_2^4*x_3^3*x_4^2*x_5*w^5.
```

At the quartic double root, exact arithmetic gives

```text
P_5(-1/a)=-2*p^2*(x_5-U)/(a^2*B)
           =-342696672409/5935078556000000000<0,
Disc(P_5)=-1705744471398225051469144773801302220092805577302905409782438390770715908617031/25059098034356086390871870668800000000000000000000000000000000000000000000000000000000000000000000000<0.
```

A real quintic with nonzero negative discriminant has exactly
three distinct real roots and one nonreal conjugate pair. Thus
this survivor is not degree-five Jensen hyperbolic. The finite
signed-Hankel conditions tested here are therefore strictly
weaker than the adjacent-degree condition needed at contact.

This is a finite countermodel gate. It is not the Xi coefficient
sequence, an all-order sign-regular sequence, or a Newman
trajectory. Adjacent-degree closure excludes this witness;
global far-column/all-order structure, theta arithmetic, and
heat compatibility remain candidate ways to establish such
closure for Xi. No
PF-infinity, `Lambda<=0`, RH, or Clay-prize conclusion follows.
