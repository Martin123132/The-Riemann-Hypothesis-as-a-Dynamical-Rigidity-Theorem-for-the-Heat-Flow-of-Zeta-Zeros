# Physical P_lin Ideal-Cubic Nilpotent Transport Gate

Date: 2026-08-02

Status: exact ideal four-vector and augmented nilpotent jets, physical saddle matrix, logarithmic mode transport, phase-curvature bridge, and first/second collar variation bounds. This is not a proof of the grouped roster estimate, signed flow, or RH.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_ideal_cubic_nilpotent_transport_gate.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_ideal_cubic_nilpotent_transport_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_ideal_cubic_nilpotent_transport_gate.py
```

## Four-Vector Jet

```text
F_0(x)=(1,-x^2/4-i*u_(N,x)/2,i*x/2,i*x/2)^T.
D_0=[[0,0,0,0],[0,0,i/2,i/2],[i/2,0,0,0],[i/2,0,0,0]].
F_0'=D_0F_0, D_0^3=0, and F_0(x)=exp(xD_0)F_0(0).
```

The row direction (0,0,1,-1) annihilates F_0(x) identically, so only alpha_A+alpha_Q and beta_A+beta_Q survive.

## Joined Augmentation

```text
Y(x)=(F_0(x),xF_0(x))^T.
Dhat=[[D_0,0],[I_4,D_0]], with Dhat^4=0 but Dhat^3!=0.
Y'=Dhat Y and Y(x)=exp(xDhat)Y(0), where the exponential terminates after its cubic term.
bar(P)=(alpha^T,i beta^T)Y; its four coefficients are the Section 11.175 bar(gamma)_0,...,bar(gamma)_3.
```

The ranks of Dhat, Dhat^2, and Dhat^3 are respectively 5, 2, 1.

## Physical Saddle Matrix

On 0<=lambda<=log(B)<ell with ell=log(a)<L, delta_a=ell-Re(alpha_s), and t=1/(2L^2), one has |x+delta_a|=|lambda-Re(alpha_s)|<L. Hence |epsilon_g|<tL/2=1/(4L)<=1/200 and t<=1/5000.

```text
epsilon_g=g_tilde+1/2=t(x+delta_a)/2. On the physical amplitude interval |epsilon_g|<1/200 and t<=1/5000.
M(epsilon_g)=Dhat^2-I_8/12+2epsilon_g Dhat+(epsilon_g^2+t/2)I_8, and L_B[Y]=M(epsilon_g)Y.
At t=epsilon_g=0, M_0=Dhat^2-I_8/12 and det(M_0)=12^(-8), so nilpotence supplies transport but no universal saddle zero.
```

In the induced matrix one-norm, ||Dhat||_1=2. Therefore ||M(epsilon_g)-(Dhat^2-I_8/12)||_1<4/200+1/40000+1/10000=161/8000.

## Exact Mode Transport

```text
d_r=log(1+1/r), x_(r+1)=x_r-d_r, and epsilon_(r+1)=epsilon_r-t*d_r/2.
E(d)=exp(-dDhat)=I-dDhat+d^2Dhat^2/2-d^3Dhat^3/6, so Y_(r+1)=E(d_r)Y_r.
E(d_2)E(d_1)=E(d_1+d_2) exactly.
```

On the terminal collar, d_r=log(1+1/r), Delta G=1-alpha*d_r, and d_r-d_(r+1)=Delta^2G/alpha. Thus 0<d_r-d_(r+1)<=5/(4alpha^2), while d_r<=exp(1/10)/alpha<10/(9alpha).

The exact second field difference is

```text
Delta^2Y_r=[(d_r-d_(r+1))Dhat+(((d_r+d_(r+1))^2)/2-d_r^2)Dhat^2+(2d_r^3-(d_r+d_(r+1))^3)Dhat^3/6]Y_r.
```

and the certified collar bounds are

```text
Using ||Dhat||_1=2, ||Dhat^2||_1=5/2, and ||Dhat^3||_1=3/2, the exact transport gives ||Y_(r+1)-Y_r||_1<(9/(4alpha))||Y_r||_1.
The exact second-transport polynomial and the phase-curvature bounds give ||Y_(r+2)-2Y_(r+1)+Y_r||_1<(9/alpha^2)||Y_r||_1.
```

For the physical saddle field itself, no finite differencing or mode enumeration is needed:

```text
Delta[M_rY_r]=[M(epsilon_r-td_r/2)E(d_r)-M(epsilon_r)]Y_r.
Delta^2[M_rY_r]=[M(epsilon_r-t(d_r+d_(r+1))/2)E(d_r+d_(r+1))-2M(epsilon_r-td_r/2)E(d_r)+M(epsilon_r)]Y_r.
```

## Handoff

The ideal cubic field now has first variation O(alpha^-1) and second variation O(alpha^-2) on exactly the collar where the phase curvature is O(alpha^-1). Insert the exact transported saddle matrix into the weighted Abel identity before bounding the Fresnel integral or endpoint packages.

The exact identity d_r-d_(r+1)=Delta^2G/alpha ties the second variation of the coefficient jet to the already certified phase curvature. This is the structural input needed for a second-order or blockwise Abel calculation; it does not itself perform that calculation.

The invertible core M_0 and arbitrary joined rows prevent a universal zero or sign. Slow mode transport is an Abel input, not a grouped roster bound or signed flow theorem.

## Proof Boundary

This gate proves no numerical retained-observation bound, grouped ideal-cubic c-prime or roster estimate, endpoint cancellation, quadratic residual bound, h^2 flow error, Phi_B upper bound, contact exclusion, retained aggregate or Xi theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion.

## Reproduce

```powershell
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_ideal_cubic_nilpotent_transport_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_ideal_cubic_nilpotent_transport_gate.py
```
