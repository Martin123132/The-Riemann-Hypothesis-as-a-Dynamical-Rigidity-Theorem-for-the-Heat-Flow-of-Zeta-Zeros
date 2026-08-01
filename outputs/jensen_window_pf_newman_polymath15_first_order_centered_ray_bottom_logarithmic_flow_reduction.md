# First-Order Ray-Bottom Logarithmic-Flow Reduction

Date: 2026-07-29

Status: exact fixed-cutoff L-flow and signed phase reduction with
open Xi Abel and joined-ray targets. This is not a proof of contact
exclusion, Lambda<=0, or RH.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_first_order_centered_ray_bottom_logarithmic_flow_reduction.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_ray_bottom_logarithmic_flow_reduction.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_ray_bottom_logarithmic_flow_reduction.py
```

Current result:

```text
validated Newman centered ray-bottom logarithmic-flow reduction: 24 rows, 1 exact cutoff partition, 8 logarithmic-flow identities, 2 phase-flux identities, 1 orientation transfer, 1 cutoff-join contract, 2 nonpromotion guards, 2 open Xi targets, 0 Abel gaps, 0 horizontal phase bounds
```

## Ray And Cutoff Cells

On a generic cofinal ray,

```text
t_j=25/(100+j), x(L)=4*pi*exp(L), a_j(L)^2=exp(L)+t_j/16
q_j(L)=2*t_j*L^2=50*L^2/(100+j)
L_- =max(B_epsilon,sqrt((100+j)/50)), L_+=min(101+j,(c_*+epsilon)*(100+j)/25)
```

The canonical cutoff changes at

```text
lambda_(j,N)=log(N^2-t_j/16)
I_(j,N)=[lambda_(j,N),lambda_(j,N+1)), intersected with [L_-,L_+]
At L=lambda_(j,N), a_j(L)=N and the canonical chart has cutoff N
```

Thus every equality boundary is owned explicitly. The old and new
charts are compared by the certified adjacent homotopy rather than
by identifying their complex Abel coordinates.

## Exact Logarithmic Flow

Inside one fixed-N cell,

```text
D_j=d/dL|_(t=t_j,N fixed)=x*partial_x, x=4*pi*exp(L)>0
u_(N,L)=x/(4*T_0)=exp(L)/[2*(exp(L)+t_j/16)]
widehat_gamma_n=x*gamma_n=x[-s_*'*log(n)+epsilon_n-epsilon_1]
D_j G_k=sum_(n=1)^k widehat_gamma_n*q_n
```

The positive factor `x` preserves the certified carrier signs:

```text
Re(widehat_gamma_n)<-log(2)/(16L^2)<0 for n>=2
Im(widehat_gamma_m-widehat_gamma_n)>=(x/2)*log(m/n)-16892/x>0 for n<m
```

These remain individual-carrier facts. The pair kernels from the
normalized-prefix gate still contain unrestricted sine and cosine
terms, so no aggregate sign is inferred.

## Endpoint-Complete Derivatives

Write `V_tilde_N=r_A-s_*'*u_N*r_0` and
`B_N=V_tilde_N+s_*' sum_k h_k G_k`. Then

```text
D_j Z_0=x*r_(0,x)+sum_(n=1)^N widehat_gamma_n*q_n
D_j V_tilde_N=x*r_(A,x)-x*s_*''*u_N*r_0-s_*'*u_(N,L)*r_0-s_*'*u_N*x*r_(0,x)
D_j B_N=D_j V_tilde_N+x*s_*''*sum_(k=1)^(N-1)h_k*G_k+s_*'*sum_(k=1)^(N-1)h_k*D_j G_k
D_j Z_A=(x*s_*''*u_N+s_*'*u_(N,L))*Z_0+s_*'*u_N*D_j Z_0+D_j B_N
```

The moving unit phase is retained:

```text
eta_L=i*Omega_eta*eta, Omega_eta=x*omega_eta
D_j mathsf_X=Pi_eta(D_j Z_0)-Omega_eta*Im(eta*Z_0)
D_j mathsf_A=Pi_eta(D_j Z_A)-Omega_eta*Im(eta*Z_A)
alpha=c*u_N; D_j mathcal_C_N=D_j mathsf_A-(D_j alpha)*mathsf_X-alpha*D_j mathsf_X
D_j alpha=x*Re(s_*'')*u_N+c*u_(N,L)
mathcal_C_N=Pi_eta(B_N)-b*u_N*mathsf_Y, mathsf_Y=Im(eta*Z_0)
D_j mathcal_C_N=Pi_eta(D_j B_N)-Omega_eta*Im(eta*B_N)-[x*Im(s_*'')*u_N+b*u_(N,L)]*mathsf_Y-b*u_N*[Im(eta*D_j Z_0)+Omega_eta*mathsf_X]
```

This is an exact O(N) derivative with the recurrent endpoint and
all endpoint derivatives still present.

## Signed Ray Flux

Under the still-open pointwise Abel gap, define

```text
Psi_(j,N)(L)=mathsf_X+i*mathcal_C_N/L
D_j arg(Psi)=[L*mathsf_X*D_j mathcal_C_N-mathsf_X*mathcal_C_N-L*mathcal_C_N*D_j mathsf_X]/[L^2*mathsf_X^2+mathcal_C_N^2]
```

Put

```text
E_j=D_j mathsf_X-x*mathcal_C_N
numerator=-L*x*mathcal_C_N^2+mathsf_X*(L*D_j mathcal_C_N-mathcal_C_N)-L*mathcal_C_N*E_j
|E_j|/(x*A_L)<2.03e-14 on |mathsf_X|<=delta_L
```

The negative square is the source-specific crossing core. The
remaining unsigned obstruction is the endpoint-complete term
`mathsf_X*(L*D_j mathcal_C_N-mathcal_C_N)`. At a crossing,

```text
At mathsf_X=0, D_j arg(Psi)=-L*D_j mathsf_X/mathcal_C_N
```

so the open Abel gap and the certified defect ratio imply the same
orientation in `L` as in `x`. Increasing-L and decreasing-L edges
carry opposite intersection signs.

A sufficient, but unproved, pointwise clockwise condition is

```text
|mathsf_X|*|L*D_j mathcal_C_N-mathcal_C_N|+L*|mathcal_C_N|*|E_j|<L*x*mathcal_C_N^2
```

It is recorded as a testable scalar target, not as a theorem.

## Cutoff And Successor Composition

At each cutoff use the joined path

```text
Psi_N ~ (mathsf_X+i*mathsf_A/L)_N ~ (mathsf_X+i*mathsf_A/L)_(N+1) ~ Psi_(N+1)
```

Use the terminal-shear homotopy inside each chart and the certified real first-jet homotopy between charts.

A strong pointwise gap does not bound the number of recrossings, and
ordered individual carriers do not sign the aggregate flux. The
closed successor target therefore remains

```text
After top, bottom, q<1, shoulder, cutoff, chart, and connector-endpoint pieces are joined, prove H_j<3*pi/2
```

Open handoff:

On every nonempty hard interval, work cell by cell with the exact endpoint-complete defect mathsf_X*(L*D_j mathcal_C_N-mathcal_C_N). Test whether Xi multiplicative completion and the recurrent endpoint give a signed or integrable upper bound before taking absolute values. Then compare the two oppositely oriented adjacent rays only inside the fully joined H_j ledger. Preserve q=1, W_0=0, cutoff equalities, and both adjacent charts.

The `pi` in `x=4*pi*exp(L)`, the saddle, and the cutoff cells is
the completed-zeta/Riemann-Siegel constant already recorded in the
source chain. The `2*pi` phase normalization is the period of
`exp(i*theta)`; no circle or fitted polygon is introduced here.

## Boundary

This artifact proves the exact hard-ray cutoff partition, logarithmic chain rule, scaled carrier and prefix currents, endpoint-complete normalized value and slope derivatives, Abel-scalar derivative, ray phase-flux identity, defect decomposition, and orientation transfer. It does not prove the Xi Abel gap, a sign for the aggregate defect, the joined H_j<3*pi/2 bound, q<1 or finite-shoulder closure, contact exclusion, Lambda<=0, PF-infinity, RH, or a Clay-prize conclusion.
