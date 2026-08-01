# Newman Centered Carrier-Kernel and Abel-Prefix Reduction

Date: 2026-07-27

Correction notice (2026-07-31): the endpoint Hermitian product must use
`J_a*conj(H_a)` and is superseded in its historical real-`H_a` form by the
complex endpoint source-normalization gate and Formal Core Section 11.146.
The abstract carrier kernels and Abel-prefix identity remain exact.

Status: exact endpoint-complete kernels and one continuous
all-fiber prefix scalar. The Xi-specific lower bound and
signed crossing theorem remain open; this is not a proof of RH.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_first_order_centered_carrier_kernel_abel_prefix_reduction.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_carrier_kernel_abel_prefix_reduction.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_carrier_kernel_abel_prefix_reduction.py
```

## Pi Provenance

```text
pi is the universal constant in the completed-zeta normalization xi(s)=s*(s-1)*pi^(-s/2)*Gamma(s/2)*zeta(s)/2. Equivalently, with theta(y)=sum_(m in Z)exp(-pi*m^2*y), pi^(-s/2)*Gamma(s/2)*zeta(s)=integral_0^infinity[(theta(y)-1)/2]*y^(s/2)*dy/y. For critical height T=x/2 and the heat-shifted height T_0=T+pi*t/8, the Riemann-Siegel saddle gives a^2=T_0/(2*pi)=x/(4*pi)+t/16. Hence at t=0, L=log(x/(4*pi))=log(a^2), and the fixed-t N-cell N<=a<N+1 has endpoints x_N=4*pi*(N^2-t/16), x_(N+1)=4*pi*((N+1)^2-t/16), so Delta x=4*pi*(2*N+1). The later 2*pi phase turn is the period exp(i*(theta+2*pi))=exp(i*theta). The prefix polygon is only the complex partial-sum path and defines none of these appearances of pi.
```

The same universal constant is the circumference-to-diameter
ratio of every Euclidean circle. No particular circle is
selected here, and the complex prefix polygon is not used to
construct or approximate it.

## Anchored Coordinates

```text
z_n=r_n*zeta_n, u_n=log(a/n)>=0, S=sum_(n=1)^N z_n, U=sum_(n=1)^N u_n*z_n, e=eta*r_0=(-1)^N*B_0*(T_0+i)*H_a, g=eta*r_A=(-1)^N*B_0*(T_0+i)*J_a, W_0=e+S, W_A=g+s_*'*U, s_*'=c+i*b
```

## Pairwise Kernels

The bulk radial and tangential numerators are

```text
For theta_nm=arg(z_n)-arg(z_m), P_B=Re(s_*'*U*conj(S))=c*sum_n u_n*r_n^2+sum_(n<m)r_n*r_m*[c*(u_n+u_m)*cos(theta_nm)-b*(u_n-u_m)*sin(theta_nm)]; Q_B=Im(s_*'*U*conj(S))=b*sum_n u_n*r_n^2+sum_(n<m)r_n*r_m*[b*(u_n+u_m)*cos(theta_nm)+c*(u_n-u_m)*sin(theta_nm)]
```

The endpoint-complete form is

```text
P=Re(W_A*conj(W_0))=Re(g*conj(e))+P_B+sum_n Re(g*conj(z_n)+s_*'*u_n*z_n*conj(e)); Q=Im(W_A*conj(W_0))=Im(g*conj(e))+Q_B+sum_n Im(g*conj(z_n)+s_*'*u_n*z_n*conj(e)); g*conj(e)=B_0^2*(T_0^2+1)*H_a*J_a
```

This is the ordinary Wronskian numerator in explicit carrier
coordinates. Its value is not assumed to have one sign.

## Abel Prefix Form

Finite Abel summation gives

```text
F_k=sum_(n=1)^k z_n, h_k=log((k+1)/k)>0, u_n=u_N+sum_(k=n)^(N-1)h_k; U=u_N*S+sum_(k=1)^(N-1)h_k*F_k
V_N=g-s_*'*u_N*e; W_A=s_*'*u_N*W_0+V_N+s_*'*sum_(k=1)^(N-1)h_k*F_k
```

Consequently both Hermitian kernels require only O(N) prefix
projections:

```text
For R=|W_0|, P=c*u_N*R^2+Re(V_N*conj(W_0))+sum_k h_k*[c*Re(F_k*conj(W_0))-b*Im(F_k*conj(W_0))]; Q=b*u_N*R^2+Im(V_N*conj(W_0))+sum_k h_k*[b*Re(F_k*conj(W_0))+c*Im(F_k*conj(W_0))]
```

## Continuous Crossing Scalar

Define the division-free prefix scalar

```text
For W_0=mathsf_X+i*mathsf_Y define mathcal_C_N=Re(V_N)-b*u_N*mathsf_Y+sum_k h_k*[c*Re(F_k)-b*Im(F_k)]. Then mathsf_A=Re(W_A)=mathcal_C_N+c*u_N*mathsf_X. At mathsf_X=0, mathsf_A=mathcal_C_N and Q=-mathsf_Y*mathcal_C_N, including mathsf_Y=0.
```

The terminal real correction on the contact band obeys

```text
On |mathsf_X|<=delta_L, 0<=c=t*D_x/4<exp(-L)/4 and 0<=u_N=log(a/N)<2*exp(-L/2), so |c*u_N*mathsf_X|<25000*exp(-11L/4)/|f_1|<exp(-3L/2)/(4L)*A_L<2e-35*A_L for L>=50.
```

At `L=50` its ratio to the target is below
`1.33931848090403900e-35`.

Thus one sufficient all-fiber theorem is

```text
Let epsilon_term=25000*exp(-11L/4)/|f_1|. On |mathsf_X|<=delta_L, |mathcal_C_N|>A_L+epsilon_term implies |mathsf_A|>A_L. At mathsf_X=0, sign(mathsf_A)=sign(mathcal_C_N) exactly, even when W_0=0.
```

At a real crossing the equality is exact, so the same scalar
also carries the orientation. The ordinary kernel relation is

```text
For R=|W_0|>=tau_L>delta_L, P+i*Q=W_A*conj(W_0) and W_A/W_0=(P+i*Q)/R^2. The previous ordinary inequality is equivalent to |Q|*sqrt(R^2-delta_L^2)-|P|*delta_L>A_L*R^2. At a real crossing Q=-mathsf_Y*mathcal_C_N; the kernel loses the zero fiber only through multiplication by mathsf_Y.
```

The exceptional derivative is

```text
With k_1=rho_1-u_a, partial_x mathsf_X=mathcal_C_N+(c*u_N-k_1)*mathsf_X-v_a*mathsf_Y+Re(D_(1,x))/|f_1|. At W_0=0, partial_x mathsf_X=mathcal_C_N+Re(D_(1,x))/|f_1|.
On |W_0|<tau_L, |partial_x mathsf_X-mathcal_C_N|<epsilon_term+2*exp(-L)*tau_L+2e-7*exp(-5L/4). At W_0=0 only the final D term remains.
```

This is why `mathcal_C_N`, unlike the determinant, retains
the `W_0=0` class.

## Exact Route Guards

The ordinary three-carrier model has

```text
{
  "aligned_Q": "-53/8",
  "amplitudes": "1,4/5,7/10",
  "distances": "3,2,1",
  "interpretation": "Strictly decreasing positive amplitudes and distances allow Q<0, Q>0, and Q=0 with |W_0|^2=61/148>0, even for c>0 and b=-1/2.",
  "kernel_zero": "0",
  "opposed_Q": "7/40",
  "radius_squared_at_zero": "61/148",
  "s_prime": "13/1056-i/2",
  "shared_lower_phase_cos": "-35/37",
  "shared_lower_phase_sin": "-12/37"
}
```

It reaches both current signs and an exact current zero with
`|W_0|^2=61/148`, despite strictly decreasing positive
amplitudes and distances.

The zero-fiber model has

```text
{
  "Re_W_A": "0",
  "U": "7*sqrt(65)*(-8 + I)/325",
  "W_0": "0",
  "W_A": "7*sqrt(65)*I/80",
  "abs_W_A_squared": "637/1280",
  "amplitudes": "1,3/5,2/5",
  "distances": "3,2,1",
  "global_unit_anchor": "(-8+i)/sqrt(65)",
  "interpretation": "At an exact complex-main zero, strict positive amplitude and distance ordering permit a nonzero purely imaginary W_A and zero physical real-direction slope. The Xi anchor and endpoint must enter any successful lower bound.",
  "s_prime": "1/16-i/2"
}
```

It has `W_0=0`, `W_A!=0`, and `Re(W_A)=0`. Therefore carrier
positivity, amplitude decrease, and distance ordering cannot
supply the physical directional slope. The actual Xi anchor
and endpoint phase must be used.

The prefix-sector guard is

```text
The positive h_k make a weighted prefix-sector or aggregate oscillatory estimate a sufficient Xi input, but termwise sign is not inherited from r_(n+1)<r_n, u_(n+1)<u_n, or the common coefficient direction s_*'. The exact countermodels forbid such a generic promotion.
```

The actual Xi relative carriers satisfy

```text
For chi_n=zeta_n*conj(zeta_1), |chi_n|=1, vartheta_(n,x)=Im(chi_(n,x)*conj(chi_n))=(-b)*log(n)+Im[d_(n,x)/(1+d_n)-d_(1,x)/(1+d_1)]. Since -b>=1/2, |d_n|<2189/x<1/2, and |d_(n,x)|<4223/x^2, vartheta_(n,x)>log(2)/2-16892/x^2>1/3 for n>=2. At fixed t a complete cutoff cell has Delta x=4*pi*(2N+1), so chi_2 advances by more than Delta x/3>2*pi and takes every unit-circle value.
The actual Xi relative carriers are monotone rotating, not confined to one fixed sector on a cutoff cell. Therefore a termwise cosine/sine sign or fixed common half-plane cannot be the uniform proof. Relative rotation may instead be used inside an aggregate oscillatory estimate or an explicit phase-cell decomposition.
```

At `L=50` the certified relative-current lower bound is
`3.46573590279972643e-01`.

## Live Theorem

```text
With delta_L=50000*exp(-5L/4)/|f_1|, A_L=(100000L+1)*exp(-5L/4)/|f_1|, and epsilon_term=25000*exp(-11L/4)/|f_1|, prove on q>=1 that |mathsf_X|<=delta_L implies |mathcal_C_N|>A_L+epsilon_term, and prove that the number of crossings mathsf_X=0 with mathcal_C_N>0 is below one successor turn after the prescribed boundary composition. This single division-free prefix theorem includes W_0=0.
```

The q<1 and chart boundaries are

```text
The q=2tL^2<1 boundary remains a separate multiplicity-compatible parabolic/Hermite theorem. The prefix identity remains exact there, but no uniform positive slope is imposed at t=0.
The prefixes F_k and complex kernels are canonical-chart objects. Prove the prefix theorem in one prescribed chart and transfer only mathsf_X and mathsf_A through the certified adjacent real-projection bounds.
```

This artifact proves no Xi prefix lower bound, signed
successor count, q<1 chart, contact exclusion, `Lambda<=0`,
PF-infinity, RH, or Clay-prize conclusion.
