# Geometric Prime-Power Phase-Monotonicity Gate

Date: 2026-07-29

Status: exact correction-free one-block thresholds and complete
short-chain recrossing guard. This is not a proof of the Xi Abel
gap, a successor winding bound, or RH.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_first_order_centered_geometric_prime_power_phase_monotonicity_gate.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_geometric_prime_power_phase_monotonicity_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_geometric_prime_power_phase_monotonicity_gate.py
```

Current result:

```text
validated geometric prime-power phase-monotonicity gate: 18 rows, 1 phase-derivative identity, 3 uniform monotone prime families, 1 complete short-chain recrossing guard, 2 nonpromotion guards, 1 open C1 perturbation target, 0 joined Abel gaps, 0 successor winding bounds
```

## Exact Geometric Block

For the correction-free `t=0` internal chain put

```text
a=p^(-1/2), P_M(z)=z*sum_(k=0)^(M-1)(a*z)^k
P_M(z)=z*[1-(a*z)^M]/(1-a*z)
P_M is nonzero on |z|=1 and has winding one
```

The block is zero-free on the unit circle, but zero-free winding
does not by itself control phase backtracking.

## Phase-Speed Margin

With `Theta_M(theta)=arg P_M(exp(i*theta))`,

```text
K_m(r,theta)=m*(r^2-r*cos(m*theta))/(1-2*r*cos(m*theta)+r^2)
Theta_M'=1+K_M(a^M,theta)-K_1(a,theta)
-m*r/(1-r)<=K_m(r,theta)<=m*r/(1+r)
Theta_M'>=mu_M(a):=1/(1+a)-M*a^M/(1-a^M)
```

mu_M(a)>0 implies strict phase increase by 2*pi and exactly one upward crossing for every constant external phase.

The tail penalty decreases by

```text
a<=M/(M+1) => g_(M+1)(a)<g_M(a)
```

so three uniform families follow:

```text
p>=5: M=1 is exact and M>=2 has mu_M>=mu_2(5^(-1/2))=(3-sqrt(5))/4>0
p=3, M>=4: mu_M>=(2-sqrt(3))/2>0
p=2, M>=8: mu_M>=22/15-sqrt(2)>0
```

These are exact correction-free block theorems.

## Complete Short-Chain Guard

For `a>1/2`, which includes `p=2,3`, choose the constant external
phase `alpha=pi/2`. Then

```text
alpha=pi/2; X_a(theta)=-sin(theta)-a*sin(2theta)=-sin(theta)*(1+2a*cos(theta))
0, pi, theta_a=arccos(-1/(2a)), 2*pi-theta_a
X_a'(0)=-(1+2a)<0, X_a'(pi)=1-2a<0, X_a'(theta_a)=X_a'(2*pi-theta_a)=2a*sin(theta_a)^2>0
two upward and two downward crossings on [0,2*pi)
wind(X_a+i*X_a')=-2
wind(P_2)=1
```

Thus a complete two-level prime-power block may be complex-zero-free
with winding one while its real first jet winds `-2`. This is a
generic external-phase guard, not a claim that the Xi trajectory
attains that phase.

## Surviving Route

The positive margins support the open perturbation target

```text
|Theta_actual'-Theta_geometric'|<mu_M(p^(-1/2))
```

with the heat quadratic, coefficient corrections, and moving base
phase retained. Short dyadic and ternary blocks remain explicit:

```text
Keep p=2 lengths 2..7 and p=3 lengths 2..3 inside an explicit finite block ledger; M=1 is exact
```

Rejoin all p-free bases, singleton chains, the recurrent endpoint, and cutoff homotopies before estimating the ray defect or successor winding.

Open handoff:

Derive a C1 heat-and-correction perturbation bound for the positive-margin long families. Do not sum blockwise winding bounds: insert the long/short decomposition into mathsf_X*(L*D_j mathcal_C_N-mathcal_C_N), preserving external phases and the recurrent endpoint.

The `pi` in the phase period is the standard period of
`exp(i*theta)`. The prime-power chain introduces no new value of
`pi`; the Xi saddle continues to use the completed-zeta constant.

## Boundary

This artifact proves the correction-free geometric block phase formula, elementary monotonicity margin, uniform one-turn families p>=5, p=3 with M>=4, and p=2 with M>=8, and an exact complete two-level recrossing guard for p=2,3. It does not prove C1 stability under the heat quadratic or d_n corrections, a joined p-free or endpoint theorem, the Xi Abel gap, H_j<3*pi/2, contact exclusion, Lambda<=0, PF-infinity, RH, or a Clay-prize conclusion.
