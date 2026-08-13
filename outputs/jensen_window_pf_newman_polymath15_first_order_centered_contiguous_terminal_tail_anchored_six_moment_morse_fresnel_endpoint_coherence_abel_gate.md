# Six-Moment Morse-Fresnel Endpoint-Coherence Abel Gate

Date: 2026-08-02

Status: exact phase and endpoint-coherence route gate with `0 enumerated
physical modes`, `0 endpoint-composed bounds`, `0 signed flow bounds`, and
this is not a proof of RH.

## Physical Collar

Use q=1, L>=50, a^2=exp(L)+1/(32L^2), h=1/a, alpha=xi/(2pi) on T_0-epsilon<=xi<=T_0, and the lower collar ceil(alpha exp(-1/10))<=r<=floor(alpha). Then alpha>2^49>4096.

## Discrete Saddle Phase

G_alpha(r)=alpha[log(alpha/r)-1]+r; because r is integral, e(G_alpha(r))=e(phi_r(alpha/r)).

```text
Delta G=1-alpha log(1+1/r),

Delta^2 G=alpha log((r+1)^2/[r(r+2)])=alpha log(1+1/[r(r+2)]).
```

On the collar and alpha>=4096, 1/(2alpha)<=Delta^2G<=5/(4alpha).

Delta G is strictly increasing and approaches the integer frequency zero at r=floor(alpha); no uniform first-difference gap exists.

## Abel Bound

Let R=floor(alpha), H=ceil(sqrt(alpha)), and split before R-H. On the outer part delta_r=Delta G lies in (-1/9,0).

For z_r=e(G(r)) and A_r=[e(delta_r)-1]^(-1), sum_(m<=r<=n)z_r=A_nz_(n+1)-A_mz_m+sum_(m<r<=n)(A_(r-1)-A_r)z_r.

A(delta)=-1/2-(i/2)cot(pi delta). Since delta_r increases inside (-1/2,0), its total variation telescopes vertically.

Every outer interval ending at n<=R-H has phase sum at most 1/|delta_n|<=(alpha-H+1)/(H-1).

Combining the outer interval with the terminal core proves

```text
Every contiguous interval in the collar has phase sum strictly less than 3sqrt(alpha).
```

and weighted summation by parts gives

```text
For arbitrary complex weights w_r, |sum w_r e(G(r))|<3sqrt(alpha){|w_n|+sum|w_(r+1)-w_r|}.
```

A bounded-variation norm C/alpha therefore gives only 3C/sqrt(alpha), an h-scale estimate rather than h^2.

## Terminal-Core Guard

Put K=floor(sqrt(alpha)/20). For 0<=k<=K, 0<=G(R-k)-G(R)<1/300.

The angular width is below `pi/150<1/40`, so projection onto the endpoint
phase has cosine greater than `99/100`. Together with the checked lower
coefficient floor this gives

```text
Using c_(0,r)(y_1)>=13/(100r), |sum_(k=0)^K c_(0,R-k)(y_1)e(G(R-k))|>1287/[200000sqrt(alpha)]>=1287h/200000.
```

Since a>2^25, the certified terminal subblock exceeds 200000h^2. This is a subblock/Abel nonpromotion guard, not a lower bound for the complete signed transition residual.

## Exact Endpoint Coherence

For s=log(alpha/r), Z(s)=sqrt(2[exp(-s)-1+s]), c_(j,r)(y_1)={A_j(exp(s))-A_j(1)J(exp(-s))}/[rZ(s)].

With kappa=1/(2pi i), exact integration by parts gives Q_(j,r)=kappa c_(j,r)(y_1)-kappa e(alpha log B)c_(j,r)(y_B)+kappa e(phi_r(alpha/r))integral_(y_1)^(y_B)c'_(j,r)(y)e(-y^2/2)dy.

The two phase collapses are

```text
phi_r(alpha/r)-y_1^2/2=-r,

phi_r(alpha/r)-y_B^2/2=alpha log(B)-rB.
```

Thus The lower boundary phase is e(-r)=1 for every integral mode. Because B and r are integral, the upper boundary phase is the common value e(alpha log B).

On R_alpha, |kappa sum c_(0,r)(y_1)|>91/96800. This is one exact coherent component; the upper trace and c' integral may cancel it in the full residual.

Neither saddle-phase Abel cancellation nor Hermitian pairing acts on the two separated coherent endpoint traces. They must remain joined to the physical endpoint E_j, the extracted calB_j endpoint functionals, and the interior c' integral.

## Route Decision

Do not seek the missing h^2 theorem from a coefficient-blind first- or second-difference estimate for the saddle phase. Its endpoint stationary core naturally has h scale, and residual integration by parts removes the saddle phase entirely from both physical boundary traces. Reassemble the complete lower and upper endpoint packages before estimating; apply phase cancellation only to the still-oscillatory interior term.

## Next Target

Derive a single endpoint-composed identity for E_j, the roster boundary traces, calB_j(B)-calB_j(1), and the terminal recurrence. Test exact cancellation of its coherent pieces before bounding the c' interior by grouped saddle-phase methods. Then insert all six composed values into the eight current observations.

## Pi Provenance

The phase e(x)=exp(2pi i x) supplies every discrete phase and kappa=1/(2pi i); A(delta) uses the equivalent identity 1/[e(delta)-1]=-1/2-(i/2)cot(pi delta). The rational sector comparison uses only pi<22/7. No geometric pi is fitted.

## Proof Boundary

This proves exact discrete saddle-phase differences and curvature, an explicit O(sqrt(alpha)) phase partial-sum bound, a certified h-scale terminal subblock guard for the actual lower c coefficient, and exact coherence of both residual endpoint traces. It does not lower-bound the complete signed residual, prove cancellation of the coherent endpoint package, bound either calB_j functional, prove a grouped roster or endpoint-composed h^2 estimate, evaluate the full physical remainder, prove a signed flow or Phi_B bound, exclude contact, establish a retained aggregate or Xi theorem, prove Q209, Lambda<=0, PF-infinity, RH, or a prize-level conclusion.

## Reproduce

```text
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_endpoint_coherence_abel_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_endpoint_coherence_abel_gate.py
```
