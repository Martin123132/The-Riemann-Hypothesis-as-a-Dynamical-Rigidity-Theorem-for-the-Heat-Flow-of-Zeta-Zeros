# Joined lower-cell plus ordinary K on the first height subcell

Date: 2026-08-28

Status: rigorous lower-plus-ordinary component certificate; complete `K_T`
is not yet enclosed

On `I_1=[10^10-10^-4,10^10+10^-4]`, retain

```text
O_join=-T_lower-T_upper+Gamma_defect,
T_upper=G_752+F_C+I_L.                              (LOK1)
```

With `L_t=d_t+i theta'-H'/H`, the cancellation-preserving target is

```text
K_(L+O)=L_t[V_L-T_lower-G_752-F_C-I_L+Gamma_defect].
                                                               (LOK2)
```

Every endpoint expansion uses the exact weighted recurrence

```text
T^n[(C-i log y)y^(-1/2)]
 =C T^n[y^(-1/2)]-i T^n[(log y)y^(-1/2)],
C=i theta'-H'/H,       T(f)=(f/(q-y-p/y))'.        (LOK3)
```

For each denominator power, its complete signed complex coefficient is
interval-evaluated on a slab before taking a modulus.  The finite source,
lower complement, and remote upper complement endpoint partials are then
joined as complex balls before any recurrence remainder is added.  The
752-label endpoint-saddle collar is not differentiated after integration:
the weight `i(theta'-log y)-H'/H` is inserted inside its one grouped
geometric-amplitude integral.  Its exact zeros at `L` and `C` survive this
multiplication.

The retained component enclosures are

```text
K_L: real [-5.550521035144591024700789451841786798604412850259680048077e-5 +/- 1.03e-8]
     imag [1.575591933754295400329504103333343296818185409077174744126e-5 +/- 3.61e-8],

K_O: real [0.001131885762403202693009739095911082263411659597604196008272 +/- 2.09e-5]
     imag [-0.0007803223146070577146973592340954016853164112284845942144995 +/- 2.10e-5],

K_(L+O): real [0.001076380552051756782762731201392664395425615469101599207791 +/- 2.09e-5]
         imag [-0.0007645663952695147606940641930620682523482293743938224670582 +/- 2.10e-5].              (LOK4)
```

The explicit sum of all endpoint-recurrence remainders is
`[7.242296992592507531913955875716181508562414791386122354447e-6 +/- 4.20e-64]`.  The Gamma defect is
retained through the exact rationalized expression for `1-C_G`; its Hardy
operator radius remains astronomically small.

After the common Hardy projection and stable positive `H` transport,

```text
Hardy_t[K_O]/H=[0.0007879918030084809288382530212402343750000000000000000000000 +/- 6.22e-5],
Hardy_t[K_(L+O)]/H=[0.0007211635411294992081820964813232421875000000000000000000000 +/- 6.21e-5] > 0. (LOK5)
```

The independent checker raises precision from 384 to 448 bits, changes the
lower finite-cell order, split, and slabs, changes the lower-complement order
`3 -> 4`, changes upper-tail order `8 -> 9`, changes all clustered slab
schemes, raises endpoint powers `16 -> 20`, and changes grouped panels
`192 -> 256`.  It also constructs changed four-label direct quadratures for
both signs of `q-y-p/y` and compares them with the weighted endpoint
recurrence.

Pi provenance: every `pi` comes from the inherited Fresnel phase, exact
half-integer Fourier spacing, root-of-unity endpoint algebra,
Riemann--Siegel phase, Gamma normalization, or `p=t/(2*pi)`.  No fitted
geometric constant supplies `pi`.

Proof boundary: (LOK1)--(LOK5) certify only the joined lower finite-cell plus
ordinary contribution on `I_1`.  They omit the positive-real tail derivative
and tiny target correction and are not yet combined with the certified
transition--upper-arc pair.  Therefore they prove no complete `K_T` bound,
wider `Q_K-T` sign interval, full event-cell theorem, wall handoff,
all-height theorem, `Lambda<=0`, RH, or prize-level conclusion.
