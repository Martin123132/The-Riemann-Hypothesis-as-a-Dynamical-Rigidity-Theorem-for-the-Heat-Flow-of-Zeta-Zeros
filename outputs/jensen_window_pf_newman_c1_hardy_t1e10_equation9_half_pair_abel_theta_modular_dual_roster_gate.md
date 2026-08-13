# Abel-theta pair triangle and modular dual roster

Date: 2026-08-13

Status: exact reduction using the Jacobi transformation; not a proof of the
quantitative source-minus-target estimate

For the exact pair triangle put

```text
r(z,x)=1/z-1/x=(x-z)/(xz)>0,
tau_epsilon=r+i epsilon,  epsilon>0,
Theta(tau)=sum_(m in Z) exp(i*pi*tau*m^2).             (AT1)
```

Abel weighting the nonzero symmetric Fourier pairs and summing before either
endpoint is separated gives

```text
S_epsilon=1/(4i*pi) integral_(0<z<x<1/2)
 Delta_g(z) a_t(x)[Theta(tau_epsilon)-1] dx dz,        (AT2)

Delta_g(z)=z^(-1/2){
 e^(i*pi*B^2*z/4)(2+i*pi*B^2*z)
-e^(i*pi*A^2*z/4)(2+i*pi*A^2*z)}.                    (AT3)
```

The endpoint difference in (AT3) is compulsory.  The elementary estimate

```text
|Delta_g(z)| <= (3*pi/2)(B^2+A^2) z^(1/2)             (AT4)
```

makes the grouped corner majorant proportional to `x^(-1/4)`, hence
integrable.  Separating the two `z^(-1/2)` endpoint pieces would lose this
argument.  The established paired `1/m^2` bound then permits
`epsilon -> 0+` by dominated convergence.

For `Im(tau)>0`, apply the standard Jacobi identity

```text
Theta(tau)=(-i*tau)^(-1/2)Theta(-1/tau).              (AT5)
```

The `n`th dual endpoint phase is

```text
Psi_(D,n)=pi*D^2*z/4+(t/2)log((1-x)/x)
            -pi*n^2*x*z/(x-z).                        (AT6)
```

Its interior stationary equations are exactly

```text
alpha_(D,n)=D-2n,
z=x alpha_(D,n)/D,
x(1-x)=2t/[pi alpha_(D,n)^2].                         (AT7)
```

Thus the modular dual index is not an arbitrary new parameter.  On the B
branch,

```text
n=0..2481422  <->  alpha=B,B-2,...,A,                     (AT8)
```

which is the complete contiguous odd source roster of `2481423` values,
in reverse.  Moreover `alpha_(B,L+k)=alpha_(A,k)=A-2k`; the two endpoint
terms label the same continuation below the source interval, although (AT2)
does not make them cancel term by term.

This is both useful and a guard.  The modular transform removes modewise
face poles and recovers the exact arithmetic roster, but it is an involutive
re-expression of the source.  A gain still requires a uniform modular
face/corner estimate that extracts the classical target saddle and bounds
the continuation jointly.

Pi provenance: `pi` is inherited from the equation-(9) Fourier/Kummer phase
and the canonical Jacobi transform.  No geometric fit or arbitrary circle
constant is introduced.

Proof boundary: exact Abel-theta reassembly, corner integrability,
Jacobi-transform reduction, and dual saddle/roster geometry only.  No
quantitative modular remainder, source-minus-target estimate, complete
`T_upper`, height-uniform theorem, `Lambda<=0`, PF-infinity, RH, or
prize-level conclusion is proved.
