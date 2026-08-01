# Newman Reciprocal-Saddle Self-Duality Gate

Date: 2026-07-28

Status: exact stationary algebra and theorem-search gate. This is not a
proof of reciprocal block cancellation, `Lambda<=0`, PF-infinity, RH,
or a Clay-prize result.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_first_order_centered_reciprocal_saddle_self_duality_gate.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_reciprocal_saddle_self_duality_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_reciprocal_saddle_self_duality_gate.py
```

## Pi Provenance

The 2*pi in a_omega^2=omega/(2*pi) is fixed by the physical logarithmic oscillation exp(i*omega*log n) written in the standard Poisson phase exp(2*pi*i*f). It is the same Riemann-Siegel saddle normalization already used in a^2=T_0/(2*pi); no fitted circle or new pi is introduced.

## Reciprocal Saddle

```text
For g_(t,sigma)(u)=exp[(t/4)log(u)^2-sigma log(u)] and the physical phase exp[i*omega log(u)], put a_omega^2=omega/(2*pi). In Poisson mode nu, the phase omega log(u)-2*pi*nu*u has the unique stationary point u_nu=a_omega^2/nu, curvature magnitude omega/(2*pi*u_nu^2)=a_omega^2/u_nu^2, and stationary factor u_nu/a_omega.
At u_nu=a_omega^2/nu, omega log(u_nu)-2*pi*nu*u_nu=-omega log(nu)+omega*(2log(a_omega)-1). Thus the physical positive logarithmic phase becomes the conjugate negative logarithmic phase, up to one global unit phase.
```

## Exact Heat Self-Duality

```text
The leading stationary amplitude obeys the exact identity g(u_nu)*(u_nu/a_omega)=g(nu)*(nu/a_omega)^Delta, where Delta=2*sigma-1-t*log(a_omega). Therefore the heat amplitude is exactly reciprocal-self-dual when Delta=0.
For s_*=(1-i*x)/2+t*alpha(s)/2, sigma=Re(s_*)=1/2+t*Re(alpha)/2 and omega=-Im(s_*)=x/2-t*Im(alpha)/2. Hence the actual reciprocal defect is exactly Delta_Xi=t*(Re(alpha)-log(a_omega)).
```

The exact identity is the main new coordinate: the stationary amplitude
is self-dual up to one scalar exponent, rather than merely comparable by
an unrelated upper bound.

## Xi Defect Bound

```text
With T_0=x/2+pi*t/8 and a_0^2=T_0/(2*pi), Im(alpha)=3x/(1+x^2)-(1/2)atan(x), so T_0-omega=(t/2)*[3x/(1+x^2)+(1/2)atan(1/x)] and 0<T_0-omega<7t/(4x). Thus 0<a_0^2-a_omega^2<7t/(8*pi*x).
On L>=50, 0<=tL<=25, the imported 0<=log(a_0)-Re(alpha)<1/(4x), together with omega>x/3, gives 0<log(a_0)-log(a_omega)<21t/(8x^2) and |Delta_Xi|<t/(2x). For a_omega<=nu<=a_omega^2 one has t*log(a_omega)<13, so exp[-13/(2x)]<(nu/a_omega)^Delta_Xi<exp[13/(2x)].
The phase saddle a_omega and the Riemann-Siegel cutoff saddle a_0 differ by less than the displayed O(t/x) square defect, but their floors may still differ on an arbitrarily thin cutoff seam. The adjacent-saddle recurrence must therefore remain in the theorem; the floors may not be identified by size alone.
```

## Critical Block

```text
A block u asymp N^r is sent to the reciprocal scale nu asymp N^(2-r). The current pointwise frontier is exposed at r_*=125662/155153, so its exact reciprocal radius is 2-r_*=184644/155153=1.190076...>1. The obstruction block inside the cutoff is paired naturally with a block outside the cutoff, not with another retained dyadic prefix.
This locates a candidate cancellation mechanism discarded by the pointwise exponent-pair envelope: compose each critical primal block with its reciprocal B-process tail and the Riemann-Siegel endpoint before taking absolute values. The c=2 pointwise exponent deficit 3133668399/48144906818 is the quantitative benchmark: at the active block, the paired estimate must gain a strictly larger power of N over the current pointwise bound.
```

## Hypothetical Stress Test

```text
In the exact hypothetical Delta=0 with corrections and endpoints suppressed, self-duality fixes only the two amplitude moduli. If the global reciprocal phase aligns, A*exp(i*theta)+A*exp(i*theta)=2A*exp(i*theta); if it is opposite, the same pair cancels. Thus reciprocal self-duality alone proves neither cancellation nor a half-plane sign. The physical Xi global phase and endpoint must decide between reinforcement and cancellation.
```

This is the stipulated impossible-or-unphysical calculation used
correctly: it reveals what the exact algebra can and cannot decide before
the physical endpoint is restored.

## Correction Guard

```text
The actual first-order carriers contain (1+d_n), the normalizer phase phi, the C_0+C_1/a endpoint, first x derivatives, and adjacent-cutoff transitions. The leading self-duality identity does not authorize dropping any of them. In particular, the existing recurrence proves that the last saddle and endpoint cancel at leading order and must be composed.
```

## Route Decision

```text
Promote reciprocal primal-tail composition as the next cancellation experiment ahead of another pointwise pair search. First derive an endpoint-complete one-block Poisson/B-process formula with value and first-x-derivative remainders uniform in the saddle fraction. Then test its signed paired main at the exact active radius r_* and at c=2 before attempting a global block sum.
```

## Open Theorem

```text
Prove, or rigorously falsify, an endpoint-complete reciprocal block theorem for the corrected Xi carriers: pair every retained block near r_* with its nu asymp N^(2-r) stationary tail, retain the global normalizer phase, C_0+C_1/a endpoint, d_n corrections, adjacent recurrence, and first x derivative, and obtain a uniform exponent gain strictly exceeding 3133668399/48144906818 at c=2. A successful theorem would lower the cancellation wall toward c=2; it would not cross the separate fixed-c<2 zero-free wall.
```

## Boundary

This artifact proves the continuous stationary-point algebra, exact reciprocal heat-amplitude identity, actual Xi defect formula and bound, critical-radius map, and two nonpromotion guards. It does not prove a discrete Poisson/B-process remainder theorem, reciprocal block cancellation, an improved c_* threshold, the Abel-scalar gap, inner degree closure, contact exclusion, Lambda<=0, PF-infinity, RH, or a Clay-prize conclusion.
