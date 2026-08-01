# Edge-Anchored Near-Terminal Pair-Current Guard

Date: 2026-07-31

Status: asymptotic one-carrier sign-reversal guard. This is not a proof of a
complete aggregate bound, Abel gap, `Lambda<=0`, or RH.

## Path

```text
Fix p in [-1,1], theta=(1-p)/2, a=N+theta, q=2tL^2=1, and n=N-m with fixed integer m>=1; let N->infinity.
```

The exact saddle expansion gives

```text
2*pi*h^-2[log(1-h(theta+m))-log(1-h theta)]+2*pi*m(h^-1-theta) -> -pi(4m theta+m^2)
z_(N-m)/z_N -> R_m(theta)=(-1)^m*exp(-4*pi*i*m*theta)
```

Hence

```text
Z_m=-Q(p)R_m(theta)=X_m+iY_m, theta=(1-p)/2, k=m+theta
```

## Four Jets

The near-terminal carrier has

```text
c_m/S_a->X_m; a*d_m/S_a->D_m=kY_m/2; a(c_(m,x)-r_Sc_m)/S_a->D_m; a^2(d_(m,x)-r_Sd_m)/S_a->N_m=Y_m/(16*pi)-k^2X_m/4
```

and the retained endpoint-terminal edge has

```text
c_e/S_a->A; a*d_e/S_a->B; a(c_(e,x)-r_Sc_e)/S_a->B; a^2(d_(e,x)-r_Sd_e)/S_a->M
```

## Pair Current

Exact polarization gives

```text
C_(e,m)=A*N_m+X_m*M-2*B*D_m
P_m(p)=(A+X_m)(M+N_m)-(B+D_m)^2=K_edge+K_m+C_(e,m)
```

The 192-bit Arb certificates are

```text
P_3(-5/8)>1/10
P_3(0)<-2
```

Thus the same edge plus `N-3` carrier has both current orientations in the
cofinal `q=1` leading model.  Convergence of all four displayed jets makes
the reversal physical for sufficiently high cells.

## Route Decision

```text
No proof may absorb the interior atoms into the signed edge one at a time by asserting K(v_edge+v_n)<0 or C_(e,n)<-K_edge-K_n uniformly.
```

This rejects termwise edge absorption.  It does not reject the complete
aggregate because the remaining carriers may supply essential signed
cancellation.

```text
The complete Xi-specific signed cross aggregate C_cross<R_diag, or the endpoint-complete division-free Abel gap, remains open. Other carriers may provide essential cancellation.
```

The pi in Q, the saddle phase, and the near-terminal ratio is inherited from the completed-zeta and Riemann-Siegel normalization. It is not fitted from a circle, polygon, or plot.

This artifact does not prove an aggregate cross-current bound, Abel gap,
winding cap, contact exclusion, `Lambda<=0`, PF-infinity, RH, or a prize-level
conclusion.
