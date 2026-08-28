# Endpoint-safe extraction of the crossing-completed B current

Date: 2026-08-13

Status: exact finite-cutoff recombination; not a proof artifact

No bound for the completed exterior current, joined remainder, or RH is
asserted.

For a common finite cutoff and Abel weight `w_(m,epsilon)`, retain

```text
tau_m=1_(622<=m<=39894),
ell_m=1_(m in {621,622}),
sigma_m=(1-tau_m)-1_(q_B,m<0),

N_(M,epsilon)
 =sum_(m=1)^M w_(m,epsilon)(1-ell_m)sigma_m P_bulk,m. (CE1)
```

The trace package by itself is not differentiable across the full released
exterior.  At every nonlocal crossing `q_B,m=0`, `U_B,m` jumps by
`-P_bulk,m`, while the omitted step in (CE1) jumps by `+P_bulk,m`.  This is
not a remote possibility: modes 620 and 623 lie in the fully-on part of the
released left and right exteriors, respectively.  More broadly, the support
with positive cutoff contains the crossing ranges

```text
left:  257..620 (364 crossings),
right: 623..1280605 (1279983 crossings). (CE2)
```

Consequently a classical first- or second-derivative bound for `Btr` alone
cannot be uniform on those exterior intervals.  Distributional jump terms
would have to be retained and cancelled against `R_join`.

The cancellation can instead be made pointwise before differentiation.  From

```text
P_B,m=U_B,m-1_(q_B,m<0)P_bulk,m
```

one obtains the exact mode identity

```text
U_B,m+P_B,-m+sigma_m P_bulk,m
 =P_B,m+P_B,-m+(1-tau_m)P_bulk,m.                    (CE3)
```

The right side is smooth for `x>0`.  Define the crossing-completed trace

```text
Ctr_(M,epsilon)=Btr_(M,epsilon)+N_(M,epsilon).        (CE4)
```

With the same `chi_W`, `E_W=1-chi_W`, and rational C2 cutoff from the
previous gate, replace the old exterior allocation by

```text
R_B,comp=eta_delta E_W Ctr_(M,epsilon),
R_join,comp=G_(M,epsilon)-chi_W Btr_(M,epsilon)-R_B,comp.

G_(M,epsilon)=chi_W Btr_(M,epsilon)
               +R_B,comp+R_join,comp.                (CE5)
```

Equivalently,

```text
R_join,comp=R_join-eta_delta E_W N_(M,epsilon).       (CE6)
```

Because `eta_delta=0` for `x<=delta`, (CE6) changes nothing at the Abel
corner: the endpoint difference, zero mode, negative modes, outer completion,
and half-current remain in the inherited common prescription.  Because
`E_W=0` on the certified window, the negative estimate
`E_Btr,win<-1.3198e-4` is also unchanged.  Away from both protected regions,
the extracted B current is now crossing-complete and classically smooth at
every finite cutoff.

This corrects the next quantitative obligation.  Tangential integration by
parts must be applied to `Ctr`, not to the discontinuous trace package alone.
One must still derive first and second normalized-amplitude bounds uniform in
`M` and `epsilon`, retain the piecewise window-edge boundary terms, and bound
the common-regulator sum with `R_join,comp` below `1.4058e-4`.

Pi provenance: every `pi` remains inherited from the exact equation-(9)
Fresnel current, its bulk completion, and the B phase.  The transfer (CE1)--
(CE6) is algebraic and introduces no geometric or fitted value of `pi`.

Proof boundary: exact finite-cutoff/common-regulator allocation, a proof that
the uncompleted exterior trace has nonlocal jumps, and smooth modewise
crossing completion only.  This gate does not prove a derivative norm,
exterior-current bound, joined-remainder bound, complete `Q_K-T` or `T_upper`,
or height-uniform theorem.  It makes no claim of `Lambda<=0`, PF-infinity,
RH, or a prize-level conclusion.
