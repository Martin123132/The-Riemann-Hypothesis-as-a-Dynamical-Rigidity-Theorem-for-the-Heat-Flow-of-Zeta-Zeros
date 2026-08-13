# Signed B cotangent current on the Gaussian face window

Date: 2026-08-13

Status: cancellation-preserving saved-height interval quadrature for one
boundary current; not a proof of the complete B-face estimate

After modes 621 and 622 are removed exactly, Section 11.365 gives the
analytic first face current

```text
C_hat(x)=-2(2+i*pi*B^2*x)S_hat(Bx/2)/pi^2.           (SC1)
```

All retained modes have the same exact trace phase because

```text
phi_(m,B)(x)+psi_m(x)
 =Phi_B(x)=pi*B^2*x/4+(t/2)log((1-x)/x).             (SC2)
```

Let `x_B` be the lower stationary point of `Phi_B`, let
`H_B=Phi_B''(x_B)`, and put `xi=sqrt(H_B)(x-x_B)`.  The quantity certified
here is the signed equation-(9) projection

```text
E_C=2*(pi/(32t))^(1/4) Re{e^(-i*pi/8)
 integral_(|xi|<=70) W_t(x)e^(i*Phi_B(x))C_hat(x)dx}. (SC3)
```

The apparent poles at `c=621,622` are never evaluated by subtracting large
singular balls.  On the chart centered at an integer `k`, the code combines
the kth pole with cotangent algebraically and evaluates the regular function

```text
pi*cot(pi*delta)-1/delta
 =-2 sum_(n>=1) zeta(2n)delta^(2n-1), delta=c-k.      (SC4)
```

For `|delta|<0.2`, eight terms are used and the omitted series is enclosed
by `4|delta|^17/(1-|delta|^2)`.  Elsewhere the exact cotangent form is used.
The 141 one-xi panels are append-only and fsynced, so the calculation is
resumable.

Arb outward rounding yields

```text
reduced complex integral
 ={'real_ball': '[-0.00464668787188490210768951783559055153056857241361053011322489350395217090472 +/- 4.28e-13]', 'imag_ball': '[0.00355029047728284612004464353856494177174100802490178950534608475529980456157 +/- 4.28e-13]'},

E_C=[4.79191597200657914096854983956099727687375298409419715416966481647269346729e-6 +/- 2.14e-15],
|E_C|<4.792e-6.                                       (SC5)
```

This is below `0.558` of the reference physical target `8.6e-6`.  The
corresponding normalized complex magnitude is

```text
[2.07023354845238357195924230249347534276643844502253962675079876715646604238e-5 +/- 2.13e-15],  (SC6)
```

which is greater than `2.07e-5`.  Thus (SC5) genuinely uses the prescribed
real projection and trace oscillation.  Taking the complex modulus or a
pointwise absolute value would not close.

This gate does not yet add the exact local 621/622 replacements, the next
two face-current terms, the nonlocal Fresnel remainder, or the outside-window
tails.  Their signed sum must be enclosed before (SC5) can be charged to a
complete B-face budget.

Pi provenance: all occurrences come from the equation-(9) face phase,
Fourier cotangent identity, and paper normalization.  The pi/8 rotation is
the exact odd-square half-Kummer reflection.  No fitted constant is used.

Proof boundary: the first pole-subtracted B face current on `|xi|<=70` at
the saved height only.  It proves no complete B estimate, no A-fold splice,
no complete paired residual or `T_upper`, no height-uniform theorem, no
`Lambda<=0`, no PF-infinity statement, no RH, and no prize-level conclusion.
