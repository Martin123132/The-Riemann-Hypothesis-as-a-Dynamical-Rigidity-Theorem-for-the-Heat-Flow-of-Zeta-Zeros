# Newman Counterfactual Birth-Signature Atlas

Date: 2026-07-25

Status: exact counterfactual birth-signature atlas and finite-sensor
nonpromotion gate. This is not a proof of `Lambda<=0` or RH.

```text
work/rh_compute/results/jensen_window_pf_newman_counterfactual_birth_signature_atlas.json
python work/rh_compute/scripts/jensen_window_pf_newman_counterfactual_birth_signature_atlas.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_counterfactual_birth_signature_atlas.py
```

Current result:

```text
validated Newman counterfactual birth-signature atlas: 12 rows, 0 issues, 1 symmetric quartet heat collision, 1 exact Jensen threshold, 1 exact Li quartet law, 1 Suzuki shift law, 1 diagonal finite-sensor escape theorem, 2 open global handoffs
```

## Counterfactual Rule

Assume one off-axis quartet provisionally and calculate every
consequence. The construction below is a proof-safety model, not an
assertion that the actual xi function has such a zero.

## Coordinate Map

```text
rho=1/2+epsilon+i*gamma with epsilon,gamma>0
z_rho=2*gamma-2*i*epsilon=x-i*y
s_rho=-z_rho^2=-4*(gamma^2-epsilon^2)+8*i*gamma*epsilon
delta=2*atan(epsilon/gamma)
sin(delta)^2=4*epsilon**2*gamma**2/(epsilon**2 + gamma**2)**2
R=4*(gamma^2+epsilon^2)=x^2+y^2
```

## Symmetric Heat Collision

Use the real even quartet polynomial

```text
P_0(z)=x**4 + 2*x**2*y**2 - 2*x**2*z**2 + y**4 + 2*y**2*z**2 + z**4
P_t(z)=12*t**2 + 4*t*x**2 - 4*t*y**2 - 12*t*z**2 + x**4 + 2*x**2*y**2 - 2*x**2*z**2 + y**4 + 2*y**2*z**2 + z**4
partial_t P_t=-partial_z^2 P_t
Disc_w(P_t)=-16*(-6*t**2 - 2*t*x**2 + 2*t*y**2 + x**2*y**2)
tau=-(x**2 - y**2 - sqrt(x**4 + 4*x**2*y**2 + y**4))/6
tau=(2/3)*(epsilon^2-gamma^2+sqrt(gamma^4+4*gamma^2*epsilon^2+epsilon^4))
P_tau(z)=(z^2-c^2)^2
```

This is an exact four-zero backward-heat birth with the correct
functional-equation symmetries. It is not the Xi heat flow.

At any nondegenerate double-zero boundary the local heat equation gives

```text
H_(tau+h)(c+q)=A*(q^2/2-h)+higher terms
H_x^2-H*H_xx=A^2*(q^2/2+h)+higher terms
L_tau(c)=0 and partial_x L_tau(c)=0
Fourier[K_(1,t)](xi)=L_t(xi/2), so a collision at c creates a zero mode at xi=2*c.
At that time the translates of K_(1,t) are not dense in L1(R).
```

Thus the strict-Laguerre correlation sensor reacts at the collision
itself. It does not need a second index to grow.

## Jensen Sensor

The conjugate squared-zero factor is

```text
G_delta,R(s)=1+2*cos(delta)*s/R+s^2/R^2
A_0=1, A_1=2*cos(delta)/R, A_2=2/R^2
J_d(X)=2*degree*variable*cos(delta)/radius + 1 + 2*variable**2*binomial(degree, 2)/radius**2
Disc(J_d)=4*degree*(-degree*sin(delta)**2 + 1)/radius**2
J_d is hyperbolic iff d*sin(delta)^2<=1
d_first=floor(csc(delta)^2)+1
csc(delta)^2=(epsilon**2 + gamma**2)**2/(4*epsilon**2*gamma**2)
d_first is of order gamma^2/(4*epsilon^2) when gamma/epsilon tends to infinity.
```

The quadratic model exactly saturates the sector cutoff. A high
near-line quartet can therefore remain invisible through an enormous
but finite Jensen degree.

## Li Sensor

For `w=(rho-1)/rho=r*exp(i*theta)`,

```text
w=(2*epsilon + 2*I*gamma - 1)/(2*epsilon + 2*I*gamma + 1)
r^2=(gamma**2 + (1/2 - epsilon)**2)/(gamma**2 + (epsilon + 1/2)**2)
kappa=-log(r)=(-log(epsilon**2 - epsilon + gamma**2 + 1/4) + log(epsilon**2 + epsilon + gamma**2 + 1/4))/2
The four Li multipliers are w, conjugate(w), 1/w, and 1/conjugate(w).
lambda_n^(Q)=4-2*(r^n+r^(-n))*cos(n*theta), where w=r*exp(i*theta) and 0<r<1.
There are infinitely many n with cos(n*theta)>=1/2; along such a subsequence the quartet contribution tends to -infinity.
kappa=-log(r) is asymptotic to epsilon/(gamma^2+1/4) for epsilon/gamma tending to zero.
```

This is the isolated quartet contribution. It does not identify the
first negative complete Xi Li coefficient.

## Suzuki Sensor

```text
For 0<omega<epsilon the shifted quotient has an uncancelled upper-half-plane pole, so D(omega) fails.
At omega=epsilon an equal-height horizontal pair can cancel in the quotient; this is one exceptional shift only.
Every sequence omega_j->0 eventually enters 0<omega_j<epsilon and therefore detects the quartet.
```

One fixed quotient can hide one exact horizontal cancellation. A
cofinal shift sequence cannot hide a quartet with positive offset.

## Diagonal Escape Theorem

Take

```text
gamma_m=m, epsilon_m=1/m, x_m=2m, y_m=2/m
tau_m=-2*(m**4 - sqrt(m**8 + 4*m**4 + 1) - 1)/(3*m**2)
D_m=(m**4 + 1)**2/(4*m**4)
kappa_m=(-log(m**2 + 1/4 - 1/m + m**(-2)) + log(m**2 + 1/4 + 1/m + m**(-2)))/2
For every fixed positive integer n, m^2*lambda_n^(Q) tends to 2*n^2; hence lambda_n^(Q)>0 for all sufficiently large m.
limits={'m2_collision_time': '2', 'jensen_threshold_over_m4': '1/4', 'm3_li_growth_rate': '1', 'li_multiplier': '1', 'fixed_index_li_contribution': '2*li_index**2'}
```

As m tends to infinity, the hypothetical quartet escapes every fixed real band, its heat boundary tends to zero, its first quadratic Jensen failure tends to infinite degree, its Li amplification rate tends to zero, and its Suzuki offset tends to zero.

Consequently, no fixed real band, finite degree, finite Li index,
positive lower time cutoff, or positive fixed Suzuki shift can be
promoted into a global exclusion theorem. At least one coordinate
must be handled cofinally or without a cutoff.

## Route Decision

Retain direct positive-time simplicity or first-jet transversality as the primary route. Positive-boundary attainment makes any Lambda>0 collision finite, and the strict-Laguerre/Wiener sensor reacts at the collision itself without a diverging degree, coefficient index, or shift sequence.

The concrete theorem remains

```text
Prove (H_t(x),H_t'(x))!=(0,0) for every real x and every 0<t<=1/5, equivalently prove L_t(x)>0 or zero-freeness of Fourier[K_(1,t)] throughout that domain.
```

Use the corrected Riemann-Siegel first-jet partition to prove phase-critical-value avoidance in the remaining high-frequency, small-t layer; alternatively prove a boundary-integrable Xi collision energy or a collision-index invariant.

The atlas selects the corrected high-frequency first-jet
phase-critical-value theorem as the primary next attack. A
boundary-integrable collision energy or a collision-index invariant
remains a separate high-risk alternative. Neither is proved here.
