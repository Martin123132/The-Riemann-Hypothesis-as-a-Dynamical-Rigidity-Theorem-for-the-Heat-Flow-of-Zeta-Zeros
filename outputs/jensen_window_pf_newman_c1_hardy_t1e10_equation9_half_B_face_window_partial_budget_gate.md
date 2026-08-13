# Saved-height B-face window partial budget

Date: 2026-08-13

Status: rigorous join of four already certified window channels; not a
complete B-face estimate

Let `R_3(x)` denote the absolutely summed nonlocal remainder left after the
third B-face integration-by-parts current.  The earlier Fresnel gate proves
uniformly on `|xi|<=70` that

```text
|R_3(x)| <= [7.470760628709106569060729806369077706582819368595599541019052558707095371060331968954272470302183300e-6 +/- 8.41e-100].       (BP1)
```

The physical equation-(9) projection contributes the factor

```text
N_t=2(pi/(32t))^(1/4),
W_t(x)=[x(1-x)]^(-1/4).                              (BP2)
```

This window lies below `x=1/2`.  Therefore `x(1-x)` is increasing and
`W_t` is decreasing, so its maximum is attained at `x_low`, not `x_high`.
Using the certified width and `W_t(x_low)` gives

```text
N_t <= [0.003540217701378688167987246287072571388372104911074868175360479126852313172094241582394783460388360052 +/- 6.88e-103],
sup W_t <= [8.014475644369682631587702791385295193833254535562130905259597073083146217278460161317213354817057757 +/- 4.50e-98],
E_rem <= [1.018469299889615412082094025565786899437817684588824096703618265830466482015461833693955967207514956e-13 +/- 1.38e-107]
      < 1.03e-13.                                    (BP3)
```

This closes the physical normalization of the nonlocal post-third-current
remainder on the Gaussian window.  It is then joined to the signed first
current, the absolute second/third currents, and the normal-safe local
Fresnel remainders without replacing the signed first current by its much
larger complex modulus:

```text
signed first current                 [4.791915972006579140968549839560997276873752984094197154169664816472693467290000000000000000000000000e-6 +/- 2.15e-15]
higher-current absolute allowance    [1.483331100424174054796745933738902400707005723781721238141235663755165833160382838319724960000000000e-6 +/- 1.14e-90]
local normal-safe remainder          [1.947962599224311902349243946192045085592310525934138072572464204156316428956422294583068841143528255e-10 +/- 4.99e-110]
nonlocal post-third remainder        [1.018469299889615412082094025565786899437817684588824096703618265830466482015461833693955967207514956e-13 +/- 1.38e-107]
Fresnel/boundary dictionary          [2.008623530541314263038626341400090849584077132853930785636956278176743558106204164082269591186308018e-20 +/- 2.07e-114]
--------------------------------------------------------------------------
certified partial upper allowance    [6.275441970537625702152377319046551824931421883618776105771895936317957884702708449530947698321351896e-6 +/- 2.15e-15]
                                   < 6.276e-6.        (BP4)
```

Against the working physical target `8.6e-6`, an absolute-value treatment
of every still-open channel therefore has at least

```text
[2.324558029462374297847622680953448175068578116381223894228104063682042115297291550469052301678648104e-6 +/- 2.15e-15] > 2.324e-6      (BP5)
```

available.  This is only accounting headroom.  The exact complementary
outer-safe currents for modes 621 and 622 and both `|xi|>70` grouped trace
tails remain unbounded, so BP4 is not an upper bound for the complete B face.

Pi provenance: `pi` in BP2 is inherited directly from the equation-(9)
normalization.  The Fresnel remainder in BP1 inherits its powers of `pi`
from repeated integration by parts of `exp(i*pi*u^2/2)`.  No fitted or
geometric surrogate for `pi` is introduced.

Proof boundary: four specified saved-height `|xi|<=70` B-face channels
only.  No complementary local current, outside-window tail, complete B
estimate, A-fold splice, complete paired residual, complete `T_upper`,
height-uniform theorem, `Lambda<=0`, PF-infinity, RH, or prize-level
conclusion is proved.
