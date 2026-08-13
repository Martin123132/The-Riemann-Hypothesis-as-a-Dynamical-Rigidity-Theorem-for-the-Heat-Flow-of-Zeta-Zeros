# Turning-event characteristic crowding and Fresnel retention

Date: 2026-08-12

Status: exact event-coordinate conversion and a rigorous obstruction to a
full-interior Fresnel replacement are certified; this is not a proof of the
replacement Airy-to-Fresnel-retaining logistic join

For fixed `C=159577`, put

```text
t*=pi*C^2/8,                 beta^3=t*,
d_m=4*beta*(m-C/4)/C,       lambda(t)=(t*-t)/beta.
```

At the event `tau_m=pi*m(C-2m)`, exact algebra gives

```text
lambda(tau_m)=d_m^2,
lambda(t)-d_m^2=(tau_m-t)/beta,                         (FR1)
Delta_m(tau_m)=(4m-C)/C=d_m/beta.                       (FR2)
```

The reduced-saddle Fresnel endpoint coordinate from Section 11.310 is

```text
z_C^*(m,t)=(tau_m-t)/sqrt(pi*(t+2*pi*m^2)).              (FR3)
```

Thus the Airy fold defect and the logistic/Fresnel transition have the same
sign and vanish at exactly the same height.  No characteristic denominator
is introduced.

## Nearest event

For mode `39894`,

```text
tau=t*-pi/8,
d=[-0.01350091479749819087428255402914909480775074027248099775869154129411374374815152748365494323780473255 +/- 4.32e-102],
d^2=pi/(8*beta).
```

On the booked exact-mode buffer `|t-tau|<=pi/16`, the two face coordinates
are

```text
z_lower=[7.833233590235854964207008052715634262226930433444529955784987085268780604113918679193164851923586219e-7 +/- 2.82e-96],
z_upper=[-7.833233590158951976011505453497225468662365699210096724768371282738999790069339783740087041661804873e-7 +/- 2.82e-96].
```

Let `F'(z)=exp(i*pi*z^2/2)`, `F(0)=0`, and
`Q(z)=(1+i)/2-F(z)`.  For real `z`, `|F(z)|<=|z|`.  Replacing `Q(z)` by the
full-interior value `1+i` therefore has relative defect at least

```text
1/2-|z|/sqrt(2)
 >= [0.4999994461067409725981231164675476910875843699215692020276358632233893861574770309402898174952508627 +/- 1.99e-96].          (FR4)
```

So the prototype face is still a half-saddle to better than one part in a
million.  It is not an ordinary full Gaussian face.

## All 399 events

The ordered event modes obey

```text
4*m_n-C=(-1)^(n+1)*(2*n+1).
```

Consequently every atlas event remains in the characteristic band

```text
[6.266567237133170820356317013103392092845460185365058874399192866139857247598338106368712283098441505e-6 +/- 3.43e-106]
 <= |Delta_m(tau_m)| <=
[0.004994454087995137143823984659443403497997831767735951922896156714313466226335875470775863689629457879 +/- 2.85e-103].
```

Even the largest symmetric event-isolation radius anywhere in the selector
cell reaches at most

```text
|z_C^*| <= [0.004975499396900639501153999739442714056481796846052703740201344559530932530466796146709242014373726908 +/- 2.82e-96].
```

Hence a full-interior Fresnel replacement has atlas-wide relative defect at
least

```text
[0.4964817906366619803303478410041072519855686350443438968319742831858731148195121628790933277373226966 +/- 2.00e-96].               (FR5)
```

This closes the proposed direct `exact buffer -> full ordinary-Morse main`
route as inadmissible.  The next valid chart must keep the exact endpoint
current and Fresnel transition inside the universal logistic Morse phase,
then compare that grouped object with the Airy-Fresnel representation before
taking absolute values.

## Phase orientation

Because `C` is odd,

```text
(-1)^m exp(-i*pi*m*C)=1
```

for every integer mode.  Since `C=1 mod 8`, the source-rotated common lower
boundary carrier is also exact:

```text
exp(i[t*-pi/8])=1.
```

No phase was discarded in obtaining (FR1)--(FR5).

## Boundary

This gate proves exact coordinate identities and rules out one nonuniform
approximation strategy.  It does not prove the Airy-to-Fresnel-retaining
logistic join, propagation through the 399 events, complete `T_upper`,
`Lambda<=0`, RH, or a prize-level conclusion.
