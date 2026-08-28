# Half-domain Kummer normalization repair

Date: 2026-08-27

Status: exact normalization repair certified; no numerical enclosure of `Q_K`

The source roster and every joined carrier in the half-Kummer branch live on
`0 <= x <= 1/2`.  The literal full-domain functional printed in Section
11.478 retained the half-domain factor `2 Re` while changing the integration
limit to `1`; read literally, that doubles every full Kummer label.

For an odd label `alpha`, odd-square parity gives

```text
K_alpha(1-x)=exp(i*pi/4) conjugate(K_alpha(x)),

exp(-i*pi/8) integral_0^1 K_alpha(x) dx
 =2 Re[exp(-i*pi/8) integral_0^(1/2) K_alpha(x) dx].   (NR1)
```

Hence the corrected physical functional and full-Kummer label are

```text
c_t=(pi/(32t))^(1/4),

P_t[F]=2c_t Re[exp(-i*pi/8) integral_0^(1/2) W_t(x)F(x) dx],

K_t(alpha)=c_t Re[exp(-i*pi/8) alpha B(a,b)
                   1F1(a;3/2;i*pi*alpha^2/4)].        (NR2)
```

The superseded line had `2c_t` in the last formula.  Its exact ratio to the
corrected label is `2`.  An independent altered-height surrogate verifies the
Euler-Kummer identity, odd-label reflection, and the corrected projection.
The builder discrepancies are respectively
`7.352141691892252849995763137635470382146e-65`,
`1.44976590044783078283848057475599270824e-66`, and
`7.270422655606398277733412262867198128446e-65`.

This repair does not rescale the A, B, Gamma, or source-carrier certificates:
their executable scopes already use the half interval.  It corrects the
source-transform description and the direct full-Kummer formula.  Therefore
the identities

```text
J_Z=Q_K-G-A_transition,
R_nonA=J_Z-E_Btr,win-E_outer-I_(A,42)
```

and their scalar sufficient corridors remain unchanged in the corrected
half-domain normalization.

Proof boundary: exact domain/factor repair and deterministic surrogate checks
only.  No value or interval for `Q_K`, `J_Z`, `Q_K-T`, or the non-A residual,
no all-height theorem, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion
is proved.
