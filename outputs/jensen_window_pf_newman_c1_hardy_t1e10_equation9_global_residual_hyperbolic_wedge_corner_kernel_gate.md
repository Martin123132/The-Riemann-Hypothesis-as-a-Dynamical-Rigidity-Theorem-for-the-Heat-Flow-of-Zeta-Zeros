# Hyperbolic wedge and half-boundary corner kernel

Date: 2026-08-13

Status: exact canonical wedge calculus and saved-height corner localization;
not a bound for the curved-face or transformed-amplitude remainder

The half-domain adds a second boundary to the flat-face kernel.  Define

```text
C_rho(a,h)=(2pi)^(-1) Abel-integral_(S>h) exp(-iS^2/2)
             integral_(P<a+rho S)exp(iP^2/2)dP dS,
c=rho^2-1.                                           (WK1)
```

For `c!=0`, differentiation and exact square completion give

```text
partial_a C_rho(a,h)
 =exp[-ia^2/(2c)]/[2sqrt(pi|c|)]
  T_sgn(c)(sqrt(|c|/pi)[h+a rho/c]),                 (WK2)
```

where `T_eps(u)=integral_u^infinity exp(i eps pi v^2/2)dv`.  The other
boundary derivative is

```text
partial_h C_rho(a,h)
 =-exp(-ih^2/2)/(2pi)
   integral_(-infinity)^(a+rho h)exp(iP^2/2)dP,       (WK3)

partial_a partial_h C_rho
 =-exp(i[(a+rho h)^2-h^2]/2)/(2pi).                  (WK4)
```

Thus the complete canonical corner is determined by one Fresnel tail and its
corner density.  Its exact limits are

```text
lim_(h to -infinity) C_rho(a,h)=F_(rho^2-1)(a),
lim_(h to +infinity) C_rho(a,h)=0.                   (WK5)
```

So the flat-face theorem is not a separate approximation: it is the remote
half-boundary limit of the same wedge kernel.

For the exact logistic Morse coordinate,

```text
h=S_H=sqrt(t/2)s(v_H),
v_H=2pi m^2/t,
s(v)=sgn(v-1)sqrt(2[v-1-log v]).                     (WK6)
```

At the B occupancy edge the half boundary is remote:

```text
m=621: h=[-270657.60667401643213927115175830276186660855822192951238949991549341993082160445 +/- 8.46e-76],
m=622: h=[-270598.16637293000615160441678107414880038865958842694839045454806844348519389793 +/- 1.70e-75].     (WK7)
```

At the A face entry it is still remote:

```text
m=39852: h=[-149.72093199180099723352120304044211221886431469263243284538180436597407779303762 +/- 4.65e-78],
m=39853: h=[-146.17478704385563843445143280421652630187952374377725048412292933237891702806400 +/- 3.46e-78]. (WK8)
```

But at the target edge it enters the canonical scale:

```text
m=39894: h=[-0.80838203032890308505596991007222797949032907004523577911376315247172917786964885 +/- 4.73e-81],
m=39895: h=[2.7365176163601948255692296213469496077615048237930160329861756225389997127278104 +/- 9.52e-81]. (WK9)
```

The exact half-boundary event is

```text
m=sqrt(t/(2pi))
 =[39894.2280401432677939946059934381868475858631164934657665925829670657925899301838501252334 +/- 9.28e-87],          (WK10)
```

while the A null-face mode is `A/4=39894.25`; their separation is
`[0.0219598567322060053940065618131524141368835065342334074170329342074100698161498747666092693 +/- 6.36e-93]`.
This quantifies why the A face, its near-null slope, and `x=1/2` cutoff must
be treated as one corner, while the B face may use the full-S profile plus a
remote half-boundary correction.

Pi provenance: `pi` comes from the equation-(9) Fourier/Kummer phase and the
exact Gaussian standardization.  The half-boundary event follows from
`2pi m^2/t=1`; no fitted geometric constant is introduced.

Proof boundary: exact Abel wedge identities and saved-height localization
only.  No curved-face or transformed-amplitude remainder, signed projector
bound, `R_Dir` estimate, complete `Q_K-T` or `T_upper`, all-height theorem,
`Lambda<=0`, PF-infinity, RH, or prize-level conclusion is proved.
