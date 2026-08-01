# Newman Cofinal Boundary-Degree Transfer Target

Date: 2026-07-26

Status: exact cofinal boundary-degree transfer target and
endpoint-uniformity guard. This is not a proof of `Lambda<=0` or RH.

```text
work/rh_compute/results/jensen_window_pf_newman_cofinal_boundary_degree_transfer_target.json
python work/rh_compute/scripts/jensen_window_pf_newman_cofinal_boundary_degree_transfer_target.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_cofinal_boundary_degree_transfer_target.py
```

Current result:

```text
validated Newman cofinal boundary-degree transfer target: 12 rows, 0 issues, 2 orientation identities, 1 vector boundary-Rouche theorem, 1 sign-definite degree composition, 1 Q207 base, 2 explicit C1 proxy budgets, 1 endpoint-uniformity countermodel, 1 open cofinal Xi boundary theorem
```

## Orientation Invariance

Write `H=A*Z`, where `A>0`, and let `ell(t,x)>0`. Then

```text
V_H=(H,H_x/ell)=A*[[1,0],[a/ell,1]]*(Z,Z_x/ell)=A*S_a,ell*V_Z
det(A*S_a,ell)=amplitude**2>0
```

The positive normalizer, derivative shear, and every positive jet scale preserve boundary nonvanishing, local degree, and winding.

At a resolved double contact of a backward-heat solution,

```text
V_F=(F,F_x/ell), where ell(t,x)>0
Use the standard (x,t) orientation, matching the counterclockwise bottom-edge traversal x=0 to x=R.
det D_(x,t)V_F=f_xx**2/scale>0
det D_(t,x)V_F=-f_xx**2/scale<0
At a spatial multiplicity-m contact, ind_(x,t)(V_F)=+floor(m/2)>0; reversing to (t,x) reverses both the index and boundary orientation.
deg(V_F,Omega,0)=wind((F+i*F_x/ell)(partial Omega),0)=+sum_p floor(m_p/2) in the standard (x,t) orientation.
```

Positive normalization and positive scaling therefore do not alter
the sign-definite contact charge.

## Boundary Rouché Transfer

On a compact boundary split the normalized scaled jet as

```text
V_Z=V_P+V_R on partial Omega
V_P=(P,P_x/ell)
V_R=(Z-P,(Z_x-P_x)/ell)
||V_R||_2<||V_P||_2 pointwise on partial Omega
V_s=V_P+s*V_R for 0<=s<=1
||V_s||_2>=||V_P||_2-s*||V_R||_2>0
```

V_Z and V_P are boundary-nonzero and have the same winding.

Combining that homotopy with the same-sign contact charge gives

```text
If the boundary domination holds and wind(V_P)=0, then Omega contains no real multiple-zero contact of Z or H.
The transferred full winding is zero, while every possible interior contact has strictly positive local index in the standard (x,t) orientation.
```

The proxy is allowed to have interior zeros of its own. Only its
boundary winding is used.

## Cofinal Contract

```text
Q_j=[1/(5*j),1/4]x[0,38+j]
Lambda<=0 iff the full Xi first jet is boundary-nonzero with zero winding on every Q_j.
Q_207=[1/1035,1/4]x[0,245] is already contact-free with zero first-jet winding.
For every j>=208, construct a continuous boundary proxy P_j, a positive scale ell_j, and certified errors satisfying strict boundary domination and wind(V_(P_j))=0.
```

The open boundary theorem for every j>=208, together with Q_207 and positive-boundary attainment, implies Lambda<=0 and RH.

This is strictly less demanding than a positive first-jet margin
throughout every two-dimensional rectangle. The arithmetic work is
moved to the bottom and right edges plus one integer winding.

## Ordinary Proxy Budget

On the existing adaptive ordinary-theta domain,

```text
O_(h,t)=J_(N_(F,h),t)^F/(16*x^4*A_t)
E_F0=exp(-3h)/(25*x^(23/4))
E_F1=E_F0*(1/2+(1+4/x)/L)
O_h^2+(O_h'/L)^2>E_F0^2+E_F1^2 only on the selected high-frequency boundary arcs.
O_h^2+(O_h'/L)^2>2*exp(-6h)/(625*x^(23/2))
```

Only selected high-frequency boundary arcs need this inequality.

## Corrected Proxy Budget

On the corrected Riemann--Siegel overlap,

```text
J_hat_(N,t)
r_RS^2+(r_RS'/L)^2<32000000*exp(-3L/2)
J_hat^2+(J_hat'/L)^2>32000000*exp(-3L/2) only on the selected boundary arcs.
```

The ordinary and corrected certified envelopes differ by five exponents; their quantitative targets must not be interchanged.

## Continuity Gate

```text
P_j must be continuous around the complete closed boundary.
Use a fixed retained cutoff on each finite boundary arc whenever the tail theorem remains valid after retaining extra terms.
Alternatively join adjacent cutoff charts through overlap homotopies that retain strict first-jet domination.
A pointwise family with unbridged cutoff jumps has no defined boundary winding and cannot discharge the theorem.
```

A collection of good pointwise estimates is not yet a winding
argument if the retained cutoff jumps without a certified join.

## Endpoint Guard

The exact model

```text
G_(epsilon,t)(x)=x^2+epsilon-2t
partial_t G=-partial_x^2 G
For t>epsilon/2 the zeros are +/-sqrt(2t-epsilon), real and simple.
t=epsilon/2+eta
first-jet norm at x=0: 4*eta**2
4*eta^2 tends to zero as eta tends to zero
```

Contact-free rectangles need not have a time-independent interior first-jet floor; boundary margins may deteriorate with the rectangle.

Thus boundary constants may deteriorate with the cofinal rectangle.
No endpoint-simplicity burden is inserted.

## Route Decision

Use the cofinal boundary-degree contract as the minimal direct collision route. It keeps the sign-definite heat-contact charge and asks arithmetic separation only on one-dimensional edges.

Retain the global corrected phase-critical small-ball theorem as a stronger alternative and the all-degree Jensen/PF route as an independent structural alternative.

The next bounded input is:

```text
A permanent low-time compact base through x=245 and a continuous proxy/winding pilot for the first boundary beyond Q_207 would test the contract without pretending that one finite stage is global.
```

The open cofinal proxy construction is the proof-level obligation.
No finite pilot is to be promoted into that theorem.
