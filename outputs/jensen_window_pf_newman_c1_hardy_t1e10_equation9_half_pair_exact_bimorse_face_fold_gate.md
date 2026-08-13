# Exact bi-Morse face and fold chart

Date: 2026-08-13

Status: exact-lemma certificate; not a proof of the uniform triangle estimate

The separable triangle phase has a second global Morse coordinate.  For the
endpoint phase put

```text
z_D=2m/D,  u=z/z_D,
p=sqrt(u)-1/sqrt(u).                                  (BM1)
```

Then, globally for `z>0`,

```text
phi_(m,D)(z)=pi*m*D+(pi*m*D/2)p^2.                   (BM2)
```

Together with the exact outer coordinate `s`,

```text
psi_m(x)=psi_m(x_m)-(t/4)s^2,                         (BM3)
```

the standardized variables

```text
P=sqrt(pi*m*D)p,   S=sqrt(t/2)s                       (BM4)
```

make the complete phase defect exactly `(P^2-S^2)/2`.  There is no phase
Taylor remainder anywhere in this bi-Morse chart.

The triangle `z<x` becomes the exact curved face

```text
P<P_D(S),
p_D(s)=sqrt[D*x(s)/(2m)]-sqrt[2m/(D*x(s))].           (BM5)
```

At the outer saddle,

```text
p_D(0)=(D-alpha_m)/sqrt(D*alpha_m),                   (BM6)

rho_D^2=[dP_D/dS(0)]^2
 =t(D+alpha_m)^2/(2*pi*m*alpha_m^3).                  (BM7)
```

At an exact endpoint event `t=pi*m(D-2m)`, `alpha_m=D`, and the quadratic
defect induced along the tangent face is

```text
rho_D^2-1=1-4m/D.                                    (BM8)
```

This explains the two qualitatively different endpoints.  At the B crossing,

```text
1-4*621/B=0.99951507304846673086807976150339849067462436219123730751533308176>0.999,
```

so the B face is safely noncharacteristic.  For the A atlas near mode 39894,

```text
1-4*39894/A=0.0000062665672371331708203563170131033920928454601853650588743991928661
            =1/159577.                                   (BM9)
```

The null face occurs exactly at `m=D/4`.  Requiring it to be an endpoint
event also gives

```text
t=pi*D^2/8,  D=sqrt(8t/pi).                           (BM10)
```

Thus the A characteristic fold and the `x=1/2` corner are not neighboring
accidents: they are the same null face in exact bi-Morse coordinates.  The
saved `A` differs slightly from `sqrt(8t/pi)`, producing the small but nonzero
detuning recorded by the interval ledger.

The prospective triangle theorem is now simpler.  Its phase is exact; it
must bound only the transformed amplitude, the curvature of (BM5), and the
finite face/corner truncations.  An ordinary Fresnel face theorem applies at
B.  A null-face fold theorem, compatible with the certified 399-event atlas,
is required at A.

Pi provenance: `pi` is inherited from the exact equation-(9) endpoint and
outer phases.  No fitted or geometric constant is introduced.

Proof boundary: exact bi-Morse coordinates, face detuning, tangent slope, and
fold-defect identity only.  No uniform amplitude/boundary-curvature bound,
complete B face estimate, A-fold splice, complete paired residual,
`T_upper`, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion follows.
