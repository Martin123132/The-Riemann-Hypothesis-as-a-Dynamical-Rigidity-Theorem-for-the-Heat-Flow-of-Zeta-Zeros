# First-Order Target Reconciliation Gate

Date: 2026-07-29

Status: exact target reconciliation and route reduction. This
is not a proof of contact exclusion, Lambda<=0, or RH.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_first_order_target_reconciliation_gate.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_target_reconciliation_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_target_reconciliation_gate.py
```

Current result:

```text
validated Newman first-order target reconciliation gate: 22 rows, 1 superseded handoff, 1 cutoff-uniform first-order remainder, 4 exact energy identities, 2 contact-box coordinates, 1 radial target, 1 box-optimal target, 1 normalized Abel target, 2 nonpromotion guards, 0 pointwise contact exclusions
```

## Supersession

The older handoff used

```text
T_L[J_[0]]>32000000*exp(-3L/2)
```

The later cutoff-uniform theorem instead gives

```text
E_[1]=g_0+sum_(n=1)^N f_n, J_[1]=2*Re(E_[1])
|r_[1]|<100000*exp(-5L/4), |r_[1],x|<200000*L*exp(-5L/4)
```

The first-order corrected main and exp(-5L/4) C1 remainder supersede the older exp(-3L/4) handoff on L>=50, 0<tL<25.

## Cutoff-Local Main

On each critical radius-1/L disk choose one prescribed-N analytic lift E_[1],N; at a cutoff, the adjacent-lift difference is retained in r_[1].

```text
Z_t=J_[1],N+r_[1],N on the real axis
```

At each point, one certified local lift with its matched remainder is enough for pointwise contact exclusion.

Boundary-degree composition must include equality boundaries and every prescribed adjacent chart/homotopy.

## Phase Energy

With the endpoint and every corrected Dirichlet component
retained,

```text
z_0=g_0; z_n=f_n for 1<=n<=N; c_j=Re(z_j), d_j=Re(z_(j,x))
X=sum_(j=0)^N c_j=Re(E_[1]), U=sum_(j=0)^N d_j=Re(E_[1],x)
T_L[J_[1]]=4[X^2+(U/L)^2]
```

For each nonzero component,

```text
For z_j=a_j*exp(i*theta_j)!=0 and z_(j,x)=(u_j+i*v_j)z_j, c_j=a_j*cos(theta_j), d_j=a_j[u_j*cos(theta_j)-v_j*sin(theta_j)].
T_L[J_[1]]/4=[sum_j a_j cos(theta_j)]^2+L^-2[sum_j a_j{u_j cos(theta_j)-v_j sin(theta_j)}]^2
```

The division-free Cartesian formula remains primary when a
component vanishes.

The cross-term structure is

```text
T_L[J_[1]]/4=|sum_j w_j|^2=sum_j|w_j|^2+2sum_(j<k) w_j dot w_k
```

w_1=(A,0), w_2=(-A,0) gives sum_j|w_j|^2=2A^2 but T_L[J_[1]]=0.

This is an algebraic nonpromotion guard, not an actual Xi configuration.

## Contact Box

A contact forces

```text
H=H_x=0 => |X|<delta_0(L) and |U|<gamma_0(L)
```

One radial sufficient target is

```text
T_L[J_[1]]>=50000000000*exp(-5L/2) excludes contact.
new_threshold/old_threshold=(3125/2)*exp(-L)<1 for L>=50
```

But the exact remainder information is rectangular. Define

```text
G_L(X,U)=max(|X|/delta_0(L),|U|/gamma_0(L))
G_L(X,U)>1, equivalently |X|<=delta_0(L) => |U|>gamma_0(L)
```

This excludes exactly the certified axis-aligned contact box; it does not replace it by a larger ellipse.

The radial condition is strictly stronger:

```text
In normalized units delta=1, gamma=2, (X,U)=(0,21/10) has G=21/20>1 but X^2+U^2=441/100<5=delta^2+gamma^2.
```

## Abel Target

The branch-free normalization and terminal shear are

```text
W_0=E_[1]/|f_1|=mathsf_X+i*mathsf_Y
mathsf_X=X/|f_1|
delta_L=50000*exp(-5L/4)/|f_1|, A_L=(100000L+1)*exp(-5L/4)/|f_1|
mathsf_A=mathcal_C_N+c*u_N*mathsf_X
epsilon_term=25000*exp(-11L/4)/|f_1|
```

The proof-facing target is

```text
|mathsf_X|<=delta_L => |mathcal_C_N|>A_L+epsilon_term
```

The Abel target gives |mathsf_A|>A_L, hence |A_a|>(100000L+1)exp(-5L/4); the real-residual bound then gives |U|>100000L exp(-5L/4).

The division-free Abel gap implies the box-optimal signed band and therefore excludes a first-order contact.

## Reduced Domain

For one fixed epsilon, the remaining outer wedge is

```text
0<epsilon<25-c_*
B_epsilon=max(50,L_epsilon)
L>=B_epsilon, q=2tL^2>=1, 0<tL<=c_*+epsilon
```

c_*+epsilon<tL<=25 is closed by the oscillatory-zeta theorem beyond L_epsilon.

tL>=25 is closed by the dominant-saddle theorem.

q<1, the bounded-L shoulder, finite joins, and the multiplicity-compatible inner degree theorem remain separate.

## Route Decision

Retire the old exp(-3L/2) target and do not attack the stronger radial first-order target by diagonal positivity. The exact main-only obligation is rectangle avoidance; the current proof-facing sufficient theorem is its robust division-free Abel form. Restrict that arithmetic work to the oscillatory-spliced low-c q>=1 wedge, and keep q<1/bounded-L degree closure as a separate theorem.

Open handoff:

On one fixed epsilon collar, derive a source-specific lower bound for the endpoint-complete Abel scalar mathcal_C_N in the band |mathsf_X|<=delta_L for L>=B_epsilon, q>=1, and 0<tL<=c_*+epsilon. Preserve W_0=0, the recurrent endpoint, q=1, equality boundaries, and adjacent-cutoff homotopies. Start from the exact prefix sum U=u_N*S+sum_k log((k+1)/k)F_k and test signed cancellation before absolute values.

Pi provenance is unchanged: `pi` is inherited from the
completed-zeta normalization and the Riemann-Siegel saddle.

## Boundary

This artifact reconciles the zeroth- and first-order targets, proves the cutoff-local Cartesian and phase/amplitude energy identities, derives the first-order contact box, radial sufficient threshold, box-optimal signed criterion, and exact Abel-to-band implication, and restricts the outer arithmetic target to the oscillatory-spliced low-c wedge. It does not prove the Xi Abel gap, a one-turn phase budget, the q<1 or bounded-L inner theorem, contact exclusion, Q209, the cofinal descendant theorem, Lambda<=0, PF-infinity, RH, or a Clay-prize conclusion.
