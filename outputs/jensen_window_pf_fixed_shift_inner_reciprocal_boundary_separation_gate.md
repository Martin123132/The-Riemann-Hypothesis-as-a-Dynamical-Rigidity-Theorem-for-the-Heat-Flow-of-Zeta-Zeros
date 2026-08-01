# Jensen-Window PF Fixed-Shift Inner/Reciprocal-Boundary Separation Gate

Date: 2026-07-23

Status: exact rational-inner/reciprocal-boundary separation countermodel
with one open arithmetic target. This is not a proof of RH, PF-infinity,
or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_fixed_shift_inner_reciprocal_boundary_separation_gate.json
python work/rh_compute/scripts/jensen_window_pf_fixed_shift_inner_reciprocal_boundary_separation_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_fixed_shift_inner_reciprocal_boundary_separation_gate.py
```

## Symmetric Four-Zero Model

Fix `omega>0`, `T>0`, and use the centered variable `z=s-1/2`. Define

```text
F(z)
 :=((z-omega)^2+T^2)((z+omega)^2+T^2).              (FSIRBSG.1)
```

This real even quartic obeys

```text
F(-z)=F(z),
F(conj(z))=conj(F(z)),
```

and has the simple zero quartet

```text
z=+/-omega+/-iT.                                     (FSIRBSG.2)
```

It is a finite toy completed function with the functional-equation
symmetries of one off-axis quartet. It makes no assertion that the actual
xi function has such a zero.

## The Fixed-Shift Quotient Is Inner

Form

```text
Q(z):=F(z-omega)/F(z+omega).                         (FSIRBSG.3)
```

Before cancellation,

```text
F(z-omega)
 =((z-2omega)^2+T^2)(z^2+T^2),

F(z+omega)
 =(z^2+T^2)((z+2omega)^2+T^2).
```

The boundary-zero factor cancels, leaving

```text
Q(z)
 =((z-2omega)^2+T^2)/((z+2omega)^2+T^2).             (FSIRBSG.4)
```

Its only poles are `-2omega+/-iT`, strictly in the left half-plane. On the
imaginary axis the numerator and denominator are conjugates, so

```text
|Q(it)|=1.                                            (FSIRBSG.5)
```

More explicitly, (FSIRBSG.4) is the product of the two right-half-plane
Blaschke factors

```text
(z-a)/(z+conj(a)),
a=2omega+iT and 2omega-iT.
```

Thus `Q` is rational inner on `Re(z)>0`.

The causal difference is

```text
Q(z)-1
 =-8omega*z/((z+2omega)^2+T^2).                      (FSIRBSG.6)
```

On the boundary,

```text
|Q(it)-1|^2
 =64omega^2*t^2
  /[(((t-T)^2+4omega^2)((t+T)^2+4omega^2))].
```

Using

```text
t^2/[D_-(t)D_+(t)]
 =t/(4T)[1/D_-(t)-1/D_+(t)],

D_-(t)=(t-T)^2+4omega^2,
D_+(t)=(t+T)^2+4omega^2,
```

and the elementary Cauchy integrals gives

```text
integral_R |Q(it)-1|^2 dt=16*pi*omega,

||Q-1||_(H2)^2=8omega                                  (FSIRBSG.7)
```

under the `1/(2*pi)` Hardy-norm convention. The causal boundary energy is
finite.

## Reciprocal Energy On The Shifted Line Diverges

Now inspect the reciprocal line energy

```text
I_(omega,T)
 :=integral_R dt/[(1+t^2)|F(omega+it)|^2].            (FSIRBSG.8)
```

The shifted line passes through the two zeros

```text
F(omega+iT)=F(omega-iT)=0.                            (FSIRBSG.9)
```

They are simple. Near `t=T`,

```text
1/[(1+t^2)|F(omega+it)|^2]
 ~C_(omega,T)/|t-T|^2,

C_(omega,T)
 =1/[64omega^2*T^2*(omega^2+T^2)*(1+T^2)]>0.          (FSIRBSG.10)
```

Therefore

```text
I_(omega,T)=infinity.                                  (FSIRBSG.11)
```

The same boundary zero is removable in the fixed-shift quotient but is a
pole of the reciprocal shifted-line function.

## Proof-Safety Consequence

This model rejects the generic promotion

```text
fixed-shift innerness or finite limiting causal energy
  =>
finite reciprocal square energy on the same shifted line.              (false)
```

In the RH corpus, the right-hand energy is the spectral form of the
weighted-prefix/logarithmic stable-prefix target. Limiting Jordan/Burnol
innerness cannot supply it at an exceptional fixed shift through a
universal norm multiplier.

The scope is deliberately narrow. The countermodel does not reject the
audited cofinal implication from full determinant or full Burnol control
to RH. One zero quartet can align with one exceptional cancellation shift;
it cannot align with every member of a sequence `omega_j->0`.

The live routes remain:

```text
1. prove sum_N P_(alpha/2,N)/N<infinity by a direct
   Mobius/Mertens mean-square estimate; or

2. prove the full Burnol norm on one cofinal shift sequence,
   which already implies RH through the existing three-gate criterion.
```

Neither route is closed here. The logarithmic tail-energy estimate, RH,
PF-infinity, and `Lambda <= 0` all remain open.
