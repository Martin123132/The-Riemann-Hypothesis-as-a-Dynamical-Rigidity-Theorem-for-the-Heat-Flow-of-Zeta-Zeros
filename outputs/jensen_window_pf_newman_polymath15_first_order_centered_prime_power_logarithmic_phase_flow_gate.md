# Prime-Power Logarithmic Phase-Flow Gate

Date: 2026-07-30

Status: exact fixed-ray chain current and a rigorous rejection of
the absolute geometric C1 perturbation route. The replacement
heat-block starlikeness theorem remains open. This gate does not
promote RH or Lambda<=0.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_first_order_centered_prime_power_logarithmic_phase_flow_gate.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_prime_power_logarithmic_phase_flow_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_prime_power_logarithmic_phase_flow_gate.py
```

Current result:

```text
validated prime-power logarithmic phase-flow gate: 20 rows, 7 exact logarithmic-phase identities, 1 dimensionless heat profile, 1 complete dyadic counterexample, 1 rejected absolute C1 target, 1 replacement starlikeness target, 0 joined Abel gaps, 0 successor winding bounds
```

## Exact Chain Flow

On one fixed cutoff chart put `h=log(p)` and

```text
Z_(m,p)=B_m*Q, Q=sum_k A_k, A_k=a_(m,k)*w_(m,p)^k
H=sum_k k*A_k; E=sum_k epsilon_(m*p^k)*A_k
epsilon_n=d_(n,x)/(1+d_n); D_j A_k=x*(epsilon_(m*p^k)-k*h*s_*')*A_k
D_j Q=x*E-x*h*s_*'*H
```

With the positive moving angular clock,

```text
s_*'=c+i*b; tau=-x*h*b>0; rho=c/(-b)>=0
D_j Q=tau*[(i-rho)*H+E/(h*(-b))]
J_ray=|Q|^2+Re(H*conj(Q))-rho*Im(H*conj(Q))+Im(E*conj(Q))/(h*(-b))
D_j arg(wQ)/tau=J_ray/|Q|^2
```

The frozen-circle current is instead

```text
J_ang=|Q|^2+Re(H*conj(Q)); partial_theta arg[z*Q_theta]=J_ang/|Q_theta|^2
```

so a frozen angular theorem cannot be promoted silently to an
`L`-flow theorem. The p-free factor is also moving:

```text
D_j arg(B_m)=x*(beta_t'-b*log(m)); D_j arg(Z_(m,p))/tau=(beta_t'-b*log(m))/(-h*b)+J_ray/|Q|^2-1
```

## Heat Coordinates

After normalizing the correction-free coefficient at `k=0`,

```text
q=exp(-t*h^2/4), v=2*(u-delta_a)/h, u=log(a/(m*p^(M-1))), delta_a=log(a)-Re(alpha)
C_k^0/C_0^0=p^(-k/2)*q^[k*(2M-2-k+v)]
t=0 => q=1 and C_k^0/C_0^0=p^(-k/2)
C_k/C_0=(C_k^0/C_0^0)*(1+d_(m*p^k))/(1+d_m)
```

The exponent contains `k(2M-2-k+v)`. It can be order one or
larger across a long block even when the radial drift per phase
turn is tiny. Heat therefore need not be an absolute C1-small
perturbation of the `t=0` geometric polynomial.

## Complete Counterexample

The obstruction occurs in an actual complete fixed-cutoff chain:

```text
t=1/4, p=2, m=562538277, N=128*m=72004899456, M=8
a=N, L=lambda_N=log(N^2-1/64)>50, tL<25
Delta theta>2*pi*log(2)*(2N+1)>2*pi; there is theta=0 mod 2*pi in the fixed-N cell
At the aligned point 0<=u<log(1+1/N)<1/N
c_k=exp[-k*log(2)/2+(k^2-14k)*log(2)^2/16], 0<=k<=7
```

At the aligned phase, exact rational enclosures prove

```text
Theta_geo(0)>14/5; Theta_heat(0)<11/5; Theta_geo(0)-Theta_heat(0)>3/5
|Theta_actual_flow(0)-Theta_heat(0)|<10^(-8)
|Theta_actual_flow'-Theta_geometric'|>1/2 but mu_8=22/15-sqrt(2)<1/15
```

High-precision values, used only as a readable audit, are

```text
L_left=50.0000000032946126852386460561336836592922647389714028441091
mu_8=0.0524531042935716178649779424569685880969947912897185934899869
Theta_geo(0)=2.88088022903976171546835539087636474523633854204361473984335
Theta_heat(0)=2.09095017702468900136150756487373858868033051740187688727164
difference=0.789930052015072714106847826002626156556008024641737852571706
difference/mu_8=15.0597388401258544811200177324524084374267845177505482520554
```

Thus the old absolute target is false by a wide margin. The
heat deformation is not harmful in this witness: its aligned
speed remains positive. What fails is the choice of norm.

## Replacement Target

The surviving one-sided theorem is

```text
J_heat=|Q_heat|^2+Re(H_heat*conj(Q_heat))>0 on |z|=1
```

Begin with (p,M)=(2,8),(3,4),(5,2), then prove length propagation or an analytic tail theorem.

A quantitative normalized starlike-block theorem supplies both one-sided phase speed and a modulus budget for the d_n C1 perturbation.

Even after that theorem is proved,

```text
Insert any proved long-block estimate into the full p-free, short, singleton, recurrent-endpoint, cutoff, and ray-defect composition
```

The `pi` in the cutoff cell comes from the completed-zeta saddle
normalization. The `2*pi` in the phase sweep is the period of
`exp(i*theta)`. No prime polygon or block defines a new value of
`pi`.

## Boundary

This artifact proves the exact fixed-ray logarithmic derivative of an actual heat-corrected prime-power chain, separates frozen angular, radial, correction, and external-phase currents, gives the exact dimensionless heat profile, and rigorously rejects the previous absolute geometric C1 perturbation target with a complete dyadic cutoff-cell witness. It replaces that failed route by a direct heat-block starlikeness target. It does not prove that target uniformly, a corrected all-family phase margin, a joined p-free or endpoint theorem, the Xi Abel gap, H_j<3*pi/2, contact exclusion, Lambda<=0, PF-infinity, RH, or a Clay-prize conclusion.
