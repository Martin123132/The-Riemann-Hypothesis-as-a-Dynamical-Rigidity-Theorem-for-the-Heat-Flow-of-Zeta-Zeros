# Residual Placement, C1 Transfer, and Cumulative Handoff

Date: 2026-07-31

Status: exact residual-placement and C1 handoff correction. This is not a
proof of an aggregate cross-current bound, Abel gap, `Lambda<=0`, or RH.

## Result

The complete normalized first jet is

```text
V_Z=(Z,partial_x Z/L)=V_J+V_r with V_J=(J_[1],partial_x J_[1]/L) and V_r=(r_[1],partial_x r_[1]/L)
```

The endpoint-terminal edge is a component of `V_J`, not a canonical
component of the global residual.  The earlier request for a complete
`Xi` edge with four residual defects is therefore retired as a proof
obligation.  The retained edge sign and cutoff margins are unchanged.

## Allocation Gauge

For an edge row `v=(c,d)` and an allocated residual row `e=(u,w)`, put
`K(v)=det(v,partial_x v)`.  Exact expansion gives

```text
K(v+lambda e)-K(v)=lambda[c*w_x+u*d_x-d*u_x-w*c_x]+lambda^2[u*w_x-w*u_x]
```

At `lambda=1` this is exactly

```text
K(v+e)-K(v)=c*w_x+u*d_x+u*w_x-d*u_x-w*c_x-w*u_x
```

However, assigning `lambda e` to the edge and `(1-lambda)e` to the
complement leaves the complete Xi jet unchanged.  The individual edge
current varies with `lambda`; the explicit sign-changing witness is

```text
v=(1,0), v_x=(0,-1), e=(0,1), e_x=(-2,0) gives K(v+lambda e)=-1+2lambda^2
```

The four-defect algebra is correct but its proposed Xi-edge meaning
is not invariant.

## Derivative Order

If the second residual coordinate is `w=partial_x r_[1]/L`, then

```text
partial_x[(partial_x r_[1])/L]=partial_x^2 r_[1]/L-partial_x r_[1]/(xL^2)
```

An allocated component current consequently asks for a second x
derivative.  The invariant boundary transfer does not.  It uses

```text
|r_[1]|<B_0=100000exp(-5L/4), |partial_x r_[1]/L|<B_1=200000exp(-5L/4)
M_L(V_J)>1 and M_L(V_r)<1 imply M_L(V_J+sV_r)>0 for 0<=s<=1; no derivative of V_r is used
```

No derivative of `V_r` along the boundary homotopy occurs.

## Cutoff Covariance

On an adjacent overlap,

```text
On an adjacent overlap, V_(J,N)+V_(r,N)=V_(J,N+1)+V_(r,N+1)=V_Z, hence Delta V_r=-Delta V_J
```

The retained edge connector has `Delta_cut<-h/25<0`, while the whole-main
chart jump has `M_L(Delta J)<1/5`.  These are distinct joins with distinct
roles.  Neither requires an Xi-level edge allocation.

## Cumulative Current

Inside the retained main,

```text
Inside V_J only, (mathsf X,mathcal C_N)=(C_edge,D_edge)+sum_(n=1)^(N-1)(c_n,d_n)
```

and exact polarization gives

```text
K(sum_a v_a)=sum_a K(v_a)+sum_(a<b)[det(v_a,partial_x v_b)+det(v_b,partial_x v_a)]
```

The edge theorem and the interior theorem provide the diagonal reserve

```text
R_diag=-K(v_edge)-sum_(n=1)^(N-1)K(v_n)>0; K(V_J)=-R_diag+C_cross
```

but all cross currents remain.  The generic guard

```text
v_1=(cos(theta),-sin(theta)), v_2=(cos(3theta),-sin(3theta)); K(v_1)=-1, K(v_2)=-3, but v_1+v_2=0 at theta=pi/2
```

shows why strict component orientation alone cannot sign the sum.

## Next Theorem

```text
Prove C_cross<R_diag in the required orientation, or prove the division-free Abel gap |mathsf_X|<=delta_L => |mathcal C_N|>A_L+epsilon_term. Then transfer the resulting whole-main boundary theorem to Xi with the C1 homotopy.
```

No cumulative cross-current bound, Abel gap, successor
winding cap, contact exclusion, `Lambda<=0`, PF-infinity, RH, or
prize-level conclusion.
