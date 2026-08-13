# Grouped B-endpoint nonstationary enclosure

Date: 2026-08-11

Status: exact upper-endpoint phase reduction and explicit R=9
nonstationary bound validated; not a proof of the lower characteristic fold

After the 84 endpoint characters are reassembled as in the y=64 endpoint
gate, the transformed common phase at `y=y_B` reduces, up to a
`z`-independent constant, to the original source endpoint phase

```text
P_B(z)=t*z/beta+pi*B^2*x(z)/4.                         (BE1)
```

The reduction uses only

```text
x=(1-tanh(z/beta))/2,
sigma=4beta/(pi C),
beta^3=pi C^2/8.                                       (BE2)
```

Thus

```text
P_B'(z)=[t-(pi B^2/8)sech(z/beta)^2]/beta,             (BE3)
P_B''(z)=(pi B^2/(4beta^2))sech(z/beta)^2
          tanh(z/beta).                                (BE4)
```

There is no `B`-endpoint stationary point on `|z|<=9`.  Interval arithmetic
gives

```text
|P_B'| >= [4778016703.758423817360821504318873171196462743725757458275781878497479479133281849495703972286212076 +/- 1.96e-90] >4.77e9. (BE5)
```

The complete grouped endpoint amplitude, including `D_84(y_B)=-84`, obeys

```text
|A_B| <= [0.01762476186986768588593001596347203645035806409020455791291511863741364838917188709219369492669539223 +/- 1.43e-101] <0.0177,    (BE6)
|A_B'| <= [8.197773804310145597463345719806872448546045194434013316143969755202047799076806568830313660286802509e-6 +/- 7.89e-105]. (BE7)
```

One integration by parts on `[-9,9]` therefore gives

```text
|integral_-9^9 A_B(z)exp(iP_B(z))dz|
 <=2Amax/mu+18[A1max/mu+Amax*P2max/mu^2]
 <= [7.408579191603318432795220588787340757158421645674988845899048798314769461562822867182347566354574626e-12 +/- 7.47e-111]
 <7.5e-12.                                             (BE8)
```

This closes the physical `B` endpoint current for the 84 transition modes
on the saved R=9 chart.  It does not estimate the lower `C` endpoint/Fresnel
fold, whose phase is characteristic and must remain grouped.

Pi provenance: the `pi` in (BE1)--(BE4) is inherited from the source Kummer
quadratic and Fourier-Poisson character through the exact scale identities;
no fitted geometric normalization is used.

Proof boundary: a saved-height grouped `B`-endpoint enclosure on `|z|<=9`
only.  No lower characteristic-fold bound, lower-interior join, complete
`T_upper`, height-uniform source theorem, `Lambda<=0`, RH, or prize-level
conclusion is proved.
