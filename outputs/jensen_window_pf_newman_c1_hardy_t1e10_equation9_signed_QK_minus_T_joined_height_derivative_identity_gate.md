# A-free joined height derivative for Q_K-T

Date: 2026-08-28

Status: exact derivative identity and saved-height A-free reassembly
certified; no nonzero-radius sign transport yet

The finite-source transform and the classical target use the same Hardy
projection.  With

```text
D_T(t)=sum_(m=622)^39894 m^(-1/2-it),
Hardy_t[X]=2 Re[exp(i theta(t))X],
```

the exact identities `Q_K=Hardy_t[S_W]/H` and `T=Hardy_t[D_T]` give

```text
Q_K-T=Hardy_t[S_W-H D_T]/H.                         (HD1)
```

This is also obtained from the previously certified route

```text
J_Z=Q_K-G-A_transition,
A_transition=A_endpoint+G_extra,
Q_K-T=J_Z+A_transition+(G-T).                       (HD2)
```

Thus the endpoint coefficient is exactly zero.  The signed B trace and outer
term belong to a valid alternate residual decomposition, but they are not
inputs to (HD1).

Write `h=H'/H`.  Differentiating (HD1) before splitting the joined object
gives

```text
(Q_K-T)'=Hardy_t[K_T]/H,                            (HD3)

K_T=S_W'+(i theta'-h)S_W
    -H[D_T'+i theta'D_T],

D_T'=-i sum_(m=622)^39894 log(m)m^(-1/2-it).        (HD4)
```

The apparent `H'D_T` terms cancel exactly inside (HD3).  The phase and
prefactor derivatives at the saved height are

```text
theta'=[10.5939869317655556783094192036828700648107762715894148006828212957663076699198194637820325 +/- 2.87e-89],
H'/H =[-6.25000000000000000007812500000000000000238281250000000000013525390625000000001233421168368e-32 +/- 4.70e-115].        (HD5)
```

The second value is independently reproduced from the Gamma-duplication
formula.  Its scale is about `6.25e-32`; it was not obtained by subtracting
floating-point approximations.

For the source, differentiation is legal in two exact coordinates:

```text
S_W'=K_0 integral_0^infinity (3pi/4-i log x)
              x^(-s)e^(-pi x^2)G_W(x)dx,            (HD6)

S_W'=E_s integral_0^infinity [pi/2-i log(u/(2pi))]
              u^(-s)e^(iu^2/(4pi))Delta_q(u)du.     (HD7)
```

Near zero the dominating factor is `u^(-1/2)(1+|log u|)`, which is
integrable.  At infinity the fixed finite difference has exponential decay;
the half-line coordinate has Gaussian decay.  These majorants are uniform on
compact subcells of the Section 11.507 event cell.  On the rotated ray
`u=2pi exp(i*pi/4)x`, the weight in (HD7) is exactly the weight in (HD6).
An altered five-label fixture evaluates both derivatives independently and
finds absolute discrepancy
`7.18643290486294994882987682877e-70`.

At `t=10^10`, direct 384-bit summation gives

```text
T =[2.09754436597030700733662243406735756894752882149895737015543412052833543055399311964308855 +/- 3.52e-90],
T'=[1.26431301735607962973890281046974533593549334358889043799637279612599973358058382850295986 +/- 2.33e-90].                 (HD8)
```

The A-free packet assembly

```text
U_unowned+O_join+I_transition-C_G D_(39853..39894)
 +(C_G-H)D_T=S_W-H D_T                              (HD9)
```

produces

```text
Q_K-T=[-0.00322994077760085929185152053833007812500000000000000000000000000000000000000000000000000000 +/- 9.68e-7],   (HD10)
```

overlapping the saved closure and remaining strictly negative.  A second
assembly adds the natural A lift back to `J_arg`, restores the extra Gamma
roster, and reaches the same ball; this checks every cancellation sign.

Half of the saved negativity margin is reserved for transport:

| radius | allowed supremum of `|(Q_K-T)'|` |
|---:|---:|
| `0.01` | `[0.1614486401376120908783585404374205972633377052205903871560019110444647 +/- 3.06e-71]` |
| `0.001` | `[1.614486401376120908783585404374205972633377052205903871560019110444647 +/- 3.06e-70]` |
| `0.0001` | `[16.14486401376120908783585404374205972633377052205903871560019110444647 +/- 3.06e-69]` |
| `0.00001` | `[161.4486401376120908783585404374205972633377052205903871560019110444647 +/- 3.06e-68]` |

No radius is selected by this gate.  The next certificate must bound the
single joined kernel `K_T` on one row of this ladder; separately norming
`Q_K'` and `T'` would discard the cancellation that (HD3) preserves.

Pi provenance: every `pi` in (HD1)--(HD9) comes from the inherited
Riemann--Siegel phase, exact finite Fourier/Gaussian source, Mellin quarter
turn, or Gamma normalization.  No circle, polygon, fitted constant, or
visual pattern supplies `pi`.

Proof boundary: exact A-free object identity, exact joined first derivative,
rigorous saved-height target value/derivative, and a saved-height A-free
negative reassembly only.  No bound for `K_T` on a nonzero-radius interval,
off-height sign, event-wall handoff, all-height theorem, equation-(4) global
infinite-series identity, `Lambda<=0`, PF-infinity, RH, or prize-level
conclusion is proved.
