# Newman First-Jet Winding Gate

Date: 2026-07-26

Status: exact all-multiplicity contact index and boundary-flux
reduction. The Xi edge and zero-winding family remains open; this
is not a proof of `Lambda <= 0` or RH.

```text
work/rh_compute/results/jensen_window_pf_newman_first_jet_winding_gate.json
python work/rh_compute/scripts/jensen_window_pf_newman_first_jet_winding_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_first_jet_winding_gate.py
```

## Universal Contact Charge

The domain has the standard (x,t) orientation: the positively oriented half-rectangle runs along its bottom edge from x=0 to x=R. Reversing the coordinate order to (t,x) reverses every local index and the boundary winding.

At a spatial zero of multiplicity `m>=2`, parabolic rescaling gives

```text
P_m(tau,u)=exp(-tau*D_u^2)u^m=(2tau)^(m/2)He_m(u/sqrt(2tau)) for tau>0
Disc_u P_m=(2tau)^(m(m-1)/2)*product_(k=1)^m k^k
For tau>0, P_m has m simple real Hermite roots. For tau<0, its coefficients have one sign in u^2, so it has zero real roots when m is even and only u=0 when m is odd. The contact therefore births floor(m/2) real pairs.
```

The identity `He_m'=m He_(m-1)` and the three-term Hermite
recurrence evaluate the resultant at all roots of `He_(m-1)`,
giving the displayed discriminant product.

At every resolved double contact,

```text
At a nondegenerate common zero of (F,F_x), det D_(x,t)(F,F_x)=F_x*F_xt-F_t*F_xx=F_xx^2>0; in the reversed (t,x) coordinate order the determinant is -F_xx^2<0
```

Degree stability therefore gives

```text
At a real zero of spatial multiplicity m>=2 in any nonzero real-analytic solution F_t=-F_xx, the common-zero point is isolated and its local Brouwer index in the standard (x,t) orientation is +floor(m/2)
If Omega is a compact planar domain whose boundary contains no common zero, wind((F+iF_x)(partial Omega),0)=+sum_(p in Omega) floor(m_p/2)>=0 in the standard (x,t) orientation
Under boundary nonvanishing, wind(F+iF_x)=0 iff Omega contains no real multiple-zero contact
```

For the local degree, perturb inside a fixed parabolic box until
all contacts are double. Each has index `+1` in the standard
`(x,t)` orientation; the real-root jump
is `2floor(m/2)`, so there are exactly `floor(m/2)` pair births.
Homotopy invariance then returns the stated multiple-contact index.

Every interior charge is strictly positive. Hidden contacts cannot
cancel one another in the boundary winding.

| m | real roots before | real roots after | pair births | index | discriminant tau power |
|---:|---:|---:|---:|---:|---:|
| 2 | 0 | 2 | 1 | 1 | 1 |
| 3 | 1 | 3 | 1 | 1 | 3 |
| 4 | 0 | 4 | 2 | 2 | 6 |
| 5 | 1 | 5 | 2 | 2 | 10 |
| 6 | 0 | 6 | 3 | 3 | 15 |
| 7 | 1 | 7 | 3 | 3 | 21 |
| 8 | 0 | 8 | 4 | 4 | 28 |
| 9 | 1 | 9 | 4 | 4 | 36 |
| 10 | 0 | 10 | 5 | 5 | 45 |

## Laguerre Phase Flux

Put `Z=H+iH_x`. Away from common zeros,

```text
Q=H^2+H_x^2
L=H_x^2-H*H_xx
partial_x arg(H+iH_x)=-L/Q
partial_t arg(H+iH_x)=L_x/Q
d arg(H+iH_x)=-(L/Q)dx+(L_x/Q)dt
```

Thus winding is not merely qualitative topology: it is a boundary
integral of the first Laguerre expression and its spatial flux.

## Xi Half-Rectangle Criterion

```text
Let Q_j=[1/(5j),1/4]x[0,38+j]. Then Lambda<=0 iff, for every j, Z=H+iH_x is nonzero on partial Q_j and wind(Z(partial Q_j),0)=0
On x=0, H_t(0)>0. On t=1/4, all zeros are simple because 1/4>1/5>=Lambda. Thus only the bottom edge t=1/(5j), 0<=x<=38+j and the right edge x=38+j, 1/(5j)<=t<=1/4 need new boundary separation
For Q=[delta,T]x[0,R], counterclockwise orientation gives 2pi*wind(Z)=-integral_0^R L_delta/Q_delta dx+integral_delta^T L_x(t,R)/Q(t,R)dt+integral_0^R L_T/Q_T dx
```

The top time `1/4` is strictly above the published upper bound
`Lambda<=1/5`, and the symmetry axis has `H_t(0)>0`. The only
new separation estimates are therefore two one-dimensional
edges. Once those are nonzero, the displayed integer is zero
exactly when the two-dimensional stage has no contact.

## Scope Guard

```text
The universal m=2 and m=4 heat contacts have winding +1 and +2 respectively in the standard (x,t) orientation. Top-edge simplicity, axis positivity, or forward real-rootedness alone does not force zero winding.
```

## Live Handoff

```text
For every diagonal half-rectangle, prove first-jet separation on the bottom and right edges and prove that the displayed phase-flux integer is zero. This is a one-dimensional boundary theorem with sign-definite interior charges; it is not yet supplied by the compact certificate, the dominant ray, root-field balance, or generic heat-flow topology.
```
