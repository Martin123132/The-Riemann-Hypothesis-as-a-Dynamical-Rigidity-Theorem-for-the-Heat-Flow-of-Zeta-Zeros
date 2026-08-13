# Incomplete-Gamma Volterra kernel and outer bound

Date: 2026-08-13

Status: interval certificate; not a proof of the endpoint-driver coupling bound

The Volterra transport uses

```text
G_m(z)=integral_z^(1/2) W_t(x)x^(-3/2)
                       exp(-i*pi*m^2/x)dx.             (VK1)
```

Put `w=(1-x)/x` and `W_z=(1-z)/z`.  The Jacobian cancels every
remaining rational power of `x`:

```text
G_m(z)=(-1)^m integral_1^(W_z)
 w^(-1/4+i*t/2)exp(-i*pi*m^2*w)dw.                    (VK2)
```

Thus, with `a=3/4+i*t/2` and `lambda=i*pi*m^2`, this is exactly the
path-defined incomplete-Gamma difference

```text
G_m=(-1)^m lambda^(-a)
 [gamma_path(a,lambda W_z)-gamma_path(a,lambda)].      (VK3)
```

The path is the straight positive-imaginary ray, so (VK3) carries no hidden
branch choice.

For `m>=39895`, the real phase in (VK2) has

```text
phi_m'(w)=t/(2w)-pi*m^2<0,       w>=1.                 (VK4)
```

Write `lambda_m(w)=pi*m^2-t/(2w)`.  Both `lambda_m` and the decay of
`w^(-1/4)` are monotone in the favorable direction.  One integration by
parts and exact total variation therefore give

```text
sup_(0<z<=1/2)|G_m(z)| <=2/[pi*m^2-t/2].               (VK5)
```

Let `nu=sqrt(t/(2*pi))` and `M=39895`.  Summing (VK5), with the decreasing
tail bounded by its integral, yields

```text
sum_(m=M)^infinity sup_z |G_m(z)|
 <=2/[pi(M^2-nu^2)]
   +log[(M+nu)/(M-nu)]/(pi*nu)
 =[0.0001024592139250777024601096162130592540178008608642272364124146630390132760954890507227330639829848688 +/- 1.23e-101] <0.000103.         (VK6)
```

The first mode has spectral gap

```text
pi*M^2-t/2=[193503.5125207458075972594709574470075130898380465352053980618812507859210095476705901020984629221914 +/- 8.42e-92],
```

and individual bound `[1.033572969268750064895190650356156205407155144244537162722061543848830995157131499163462569357709154e-5 +/- 4.50e-102]`.

This proves that the complete outer Volterra-kernel family is summable and
small as an operator kernel.  It does not yet bound the pair residual:
placing an absolute value on the endpoint driver `Delta f_z'` would discard
its `A/B` phase cancellation and produce a useless constant.  The next step
must combine (VK2) with that driver before taking norms, isolating only the
characteristic alignments proved in the transport gate.

Pi provenance: `pi` is inherited from the equation-(9) Fourier/Kummer phase
and its saddle threshold.  No fitted constant is used.

Proof boundary: exact kernel/incomplete-Gamma representation and rigorous
outer-kernel `l1` bound only.  No endpoint-driver convolution bound, complete
outer-pair residual, B crossing, A-fold splice, complete `T_upper`,
height-uniform theorem, `Lambda<=0`, PF-infinity, RH, or prize-level
conclusion is proved.
