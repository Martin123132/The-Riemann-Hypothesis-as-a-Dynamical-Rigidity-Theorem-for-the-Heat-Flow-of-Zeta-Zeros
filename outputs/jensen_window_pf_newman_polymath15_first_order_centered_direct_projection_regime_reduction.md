# Newman First-Order Centered Direct-Projection Regime Reduction

Date: 2026-07-27

Correction notice (2026-07-31): the displayed endpoint projections that
treat `H_a` as real are superseded by the complex endpoint
source-normalization gate and Formal Core Section 11.146.  The carrier,
derivative, and ordinary/exceptional transport identities remain exact.

Status: exact direct saddle-phase coordinates and an
ordinary/exceptional theorem split. The Xi-specific lower
bounds remain open; this is not a proof of RH.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_first_order_centered_direct_projection_regime_reduction.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_direct_projection_regime_reduction.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_direct_projection_regime_reduction.py
```

## Direct Unit Carriers

The relative coefficient chain has the branch-free split

```text
sigma_*=Re(s_*), omega_*=-Im(s_*); r_n=exp[t*log(n)^2/4-sigma_*log(n)]*|1+d_n|/|1+d_1|>0; zeta_n=phi*exp(i*omega_*log(n))*(1+d_n)/|1+d_n|, |zeta_n|=1; eta*q_n=r_n*zeta_n
```

In the saddle-distance coordinate `u_n=log(a/n)`,

```text
u_n=log(a/n)=-ell_n, delta_a=log(a)-Re(alpha), B_a=1/2-t*delta_a/2, E_a=t*log(a)^2/4-Re(s_*)*log(a); r_n=exp[E_a+B_a*u_n+t*u_n^2/4]*|1+d_n|/|1+d_1|
0<=delta_a<1/(4x), 0<=t<=1/2, hence B_a=1/2-t*delta_a/2>49/100. The amplitudes are explicit and positive; all unresolved oscillation is in the unit carriers zeta_n.
```

The numerical audit at `L=50` gives
`B_a>5.00000000000000000e-01`. Positivity and amplitude
ordering are explicit; the unresolved arithmetic is entirely
in the unit carriers.

The endpoint enters the same frame without an argument branch:

```text
B_0=beta/(|M_t(s)|*|f_1|)>0, J_a=H_(a,x)+mu_a*H_a; Re(eta*r_0)=(-1)^N*B_0*T_0*H_a; Re(eta*r_A)=(-1)^N*B_0*[T_0*Re(J_a)-Im(J_a)]
```

Therefore the two live real observables are

```text
C_n=Re(zeta_n), S_n=Im(zeta_n), s_*'=c+i*b; mathsf_X:=Re(eta*Z_0)=(-1)^N*B_0*T_0*H_a+sum_n r_n*C_n; mathsf_A:=Re(eta*Z_A)=(-1)^N*B_0*[T_0*Re(J_a)-Im(J_a)]+sum_n ell_n*r_n*(-c*C_n+b*S_n)
```

## Derivative Identity

Differentiating `E_[1]=f_1 Z_0` and comparing with the
centered decomposition gives

```text
nu_1=f_(1,x)/f_1=phi'/phi+d_(1,x)/(1+d_1); Z_A=Z_(0,x)+(nu_1-lambda_a)*Z_0-D_(1,x)/f_1
W_0=eta*Z_0=E_[1]/|f_1|, W_A=eta*Z_A=C_a/|f_1|, rho_1=partial_x log|f_1|=Re[d_(1,x)/(1+d_1)]; W_A=W_(0,x)+(rho_1-lambda_a)*W_0-D_(1,x)/|f_1|
For W_0=mathsf_X+i*mathsf_Y and lambda_a=u_a+i*v_a, mathsf_A=partial_x mathsf_X+(rho_1-u_a)*mathsf_X+v_a*mathsf_Y-Re(D_(1,x))/|f_1|
```

The complete nuisance coefficient is tiny:

```text
|rho_1|<8446/x^2, |u_a|<exp(-L), |v_a|<3/x^2, hence |rho_1-lambda_a|<2*exp(-L). Also |D_(1,x)|/|f_1|<2e-7*exp(-5L/4).
```

This is stronger bookkeeping than the previous qualitative
ordinary/exceptional split.

## Ordinary Regime

For `Z_0!=0`,

```text
If Z_0!=0, write h=Z_A/Z_0=p+i*q and R=|Z_0|. Then mathsf_A=p*mathsf_X-q*mathsf_Y, mathsf_X^2+mathsf_Y^2=R^2, and q=-Delta_anchor/R^2.
h=W_(0,x)/W_0+rho_1-lambda_a-D_(1,x)/E_[1]. With radial current r_x=partial_x log|W_0| and branch-free phase current theta_x=Im(W_(0,x)*conj(W_0))/|W_0|^2, p=r_x+rho_1-u_a-Re(D_(1,x)/E_[1]) and q=theta_x-v_a-Im(D_(1,x)/E_[1]).
```

Consequently,

```text
If |mathsf_X|<=delta<tau<=R, |q|>=kappa, and |p|<=K, then |mathsf_A|>=kappa*sqrt(tau^2-delta^2)-K*delta.
On |Z_0|>=tau, |D_(1,x)/E_[1]|<2e-7*exp(-5L/4)/tau, while |rho_1-u_a|<2*exp(-L) and |v_a|<exp(-2L).
Choose tau_L=L*exp(-L). Since |f_1|>1/2, delta_L/tau_L<100000*exp(-L/4)/L<1/100 and sqrt(tau_L^2-delta_L^2)>0.9999*tau_L. Also A_L/tau_L<2*(100000L+1)*exp(-L/4)/L. Thus the ordinary branch closes if 0.9999*kappa_L-K_L*100000*exp(-L/4)/L>2*(100000L+1)*exp(-L/4)/L.
```

The missing input is now explicit: a lower bound for the
phase current, with enough control of the radial current.
Since the imaginary ratio is the normalized
determinant, this does not resurrect the rejected generic
Wronskian route; it identifies exactly what Xi-specific phase
information would have to add.

At `L=50`, the coarse ratios are
`delta_L/tau_L<7.45330634415734220e-03`
and `A_L/tau_L<7.45330783481861103e-01`.

## Exceptional Regime

When `|Z_0|<tau`, no ratio is used:

```text
If |Z_0|<tau, then |mathsf_A-partial_x mathsf_X|<2*exp(-L)*tau+2e-7*exp(-5L/4). Thus |partial_x mathsf_X| exceeding the target by this amount is sufficient.
At Z_0=0, mathsf_A=partial_x mathsf_X-Re(D_(1,x))/|f_1|.
```

At the zero fiber the problem is directional real transversality, not merely simplicity of the complex zero.
The exact guard is

```text
W_*(x,t)=t-x^2/2+i*x satisfies partial_t W_*=-partial_x^2 W_*. At (0,0), W_*=0, |W_*,x|=1, and det D_(x,t)(Re W_*,Im W_*)=-1, but partial_x Re W_*=0. Complex-zero simplicity and a nonzero two-variable Jacobian do not imply the required directional real-slope bound.
```

Thus neither `|Z_(0,x)|>0` nor a nonzero `(x,t)` Jacobian
closes the exceptional branch. The direction selected by the
physical real projection must be controlled.

## Live Theorem

```text
With delta_L=50000*exp(-5L/4)/|f_1| and A_L=(100000L+1)*exp(-5L/4)/|f_1|, choose tau_L=L*exp(-L). On |Z_0|=R>=tau_L prove the direct Xi-specific inequality |Im(Z_A/Z_0)|*sqrt(R^2-delta_L^2)-|Re(Z_A/Z_0)|*delta_L>A_L. Uniform kappa_L,K_L bounds satisfying the displayed coarse condition are one stronger sufficient route, not a required formulation. On |Z_0|<tau_L prove the signed direct bound for partial_x mathsf_X, including the displayed frame and D corrections. Then count mathsf_A-positive crossings below one successor turn.
```

Chart use is restricted by

```text
mathsf_X and mathsf_A inherit the adjacent real-projection bounds, but h=Z_A/Z_0 and its real/imaginary parts do not have a certified complex adjacent-chart bound. Use h only inside one prescribed canonical chart and transfer the final real inequalities.
```

The q<1 obligation remains

```text
The q=2tL^2<1 boundary still requires a separate multiplicity-compatible parabolic/Hermite chart; neither ordinary ratio division nor a uniform zero-fiber slope may be imposed at t=0.
```

This artifact proves no ratio lower bound, exceptional
directional transversality, signed successor count, q<1
chart, contact exclusion, `Lambda<=0`, PF-infinity, RH, or
Clay-prize conclusion.
