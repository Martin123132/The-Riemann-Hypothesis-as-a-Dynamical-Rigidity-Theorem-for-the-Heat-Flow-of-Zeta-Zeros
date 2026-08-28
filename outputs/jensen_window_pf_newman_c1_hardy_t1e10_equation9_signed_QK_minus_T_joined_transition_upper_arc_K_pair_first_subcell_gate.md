# Transition--upper-arc K pair on the first subcell

Date: 2026-08-28

Status: rigorous two-component Hardy-derivative interval; complete `K_T`
remains open

On `I_1=[10^10-10^-4,10^10+10^-4]`, let `K_tr` be the cancellation-stable
transition operator packet and `K_U` the cancellation-stable compact upper-arc
operator packet.  Their common Hardy projection gives

```text
Hardy_t[K_tr]/H=[7.564537099824519827961921691894531250000000000000000000e-5 +/- 1.06e-6],               (TP1)

Hardy_t[K_U]/H =[-7.545256446744565354337042073495922522852197289466857910e-5 +/- 1.87e-11].                      (TP2)
```

The signs are opposite and the centres nearly cancel.  Joining before taking
an absolute value yields

```text
Hardy_t[K_tr+K_U]/H=[1.928065307995447362487961839860872714780271053314209000e-7 +/- 1.07e-6],                 (TP3)

sup_I1 |Hardy_t[K_tr+K_U]/H| < 1.3e-6.             (TP4)
```

For `K_tr`, the operator is pushed through the alternating boundary-jet source
and through the 42-mode Gamma sum per mode.  For `K_U`, the rotated-arc weight
is

```text
delta+i(theta'-log(Y))-H'/H,
```

and the final projection transports the combined phase `theta+phi_U` with
derivative `theta'-log(Y)`.  Thus neither component is obtained by separately
widening its large derivative and phase terms.

The independent checker raises both calculations to 448 bits.  It changes the
transition jet order, derivative slabs, integral panels, and logarithm-series
order, while changing the arc Taylor degree, endpoint cutoff, initial panels,
and all main-panel breakpoints.  Both direct operator packets overlap their
dependency-heavy assembled forms.

Pi provenance: every `pi` comes from the inherited Fresnel phase, exact
integer Fourier spacing, Gamma normalization, Riemann--Siegel phase, or
`p=t/(2*pi)`.  No circle fit, polygon fit, or visual pattern supplies `pi`.

Proof boundary: (TP1)--(TP4) cover only the transition and compact upper arc.
The positive-real tail derivative, lower finite cell, ordinary complement,
and tiny target correction remain outside the pair.  Therefore this does not
prove a complete nonzero-radius `K_T` bound, a wider `Q_K-T` sign interval,
the full event-cell sign theorem, an event-wall handoff, an all-height theorem,
`Lambda<=0`, RH, or a prize-level conclusion.
