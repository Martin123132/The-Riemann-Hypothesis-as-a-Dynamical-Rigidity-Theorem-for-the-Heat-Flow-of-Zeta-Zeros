# Grouped Airy phase and Green packets on the event lattice

Date: 2026-08-13

Status: rigorous finite all-height contiguous-packet and grouped Green-kernel
bound; no completed source or endpoint cancellation

Write the exact Airy event phase relative to event zero as

```text
vartheta(r_j)-vartheta(r_0)=2pi*j*(2q+j)+P_5(j,u)+delta_j(u),
|delta_j(u)|<[4.808566009458293869314943715213386539898294728638859740837964066980273494323764175069919281109437730e-5 +/- 2.76e-105].       (GP1)
```

The integer phase disappears inside every exponential.  At `u=0`, direct
Arb summation of the 399 quintic phasors gives

```text
sum exp(iP_5(j,0))
 =[31.72953963596475630594449123050215762526112312113965227850539315803592845817311810425278260660750704 +/- 1.67e-90]
 +i*[-0.4287931139968062603663683093523764701947280682309090531395764357265126995469250343569574427900822190 +/- 1.60e-90],
absolute value =[31.73243685954278485578872652887189699595862912155425771410357398453719491361038331669753061859562052 +/- 1.69e-90].       (GP2)
```

Enumeration of all `399*400/2` nonempty contiguous blocks proves

```text
max_(a<=b) |sum_(j=a)^b exp(iP_5(j,0))|
 <[52.29819183608537663871395043845224011992798872482164194141560173141320313495068993258744285201994760 +/- 4.78e-99],
candidate block [-30, 31].         (GP3)
```

This is finite exact-roster enumeration with Arb balls, not floating-point
sampling.  Differentiating the explicit polynomial `P_5` in `u` and using
`0<=u<=1/(2C^2)` gives

```text
sum_j sup_u |P_5(j,u)-P_5(j,0)|
 <[0.1958881303902716297327862806097029728909238386483998474824901215757648281734277151469469387002390866 +/- 4.34e-100].                 (GP4)
```

The exponential Lipschitz inequality and (GP1)--(GP4) therefore prove,
uniformly on the complete top corridor,

```text
|sum_(j=-198)^200 exp(i[vartheta(r_j)-vartheta(r_0)])|
 >[31.51736255077477463351737362283849261077351108693858881625514038633417879419460378249205470028217232 +/- 3.10e-99],

sup_(a<=b) |sum_(j=a)^b exp(i[vartheta(r_j)-vartheta(r_0)])|
 <[52.51326614485338686098530334448564450511310675943731083926403532961621925436646946679291876865181334 +/- 2.75e-99]<52.53.           (GP5)
```

The first line is a useful negative result: the cubic packet does not
collapse to an `O(1)` phasor.  The second line still gives substantial,
certified cancellation for every active contiguous roster block.

Using

```text
K(r_j,s)=pi M(-r_j)M(-s)
                 sin(vartheta(r_j)-vartheta(s)),                       (GP6)
```

DLMF 9.8(iii) makes `M(-r_j)` positive and nonincreasing in `j`.
Finite Abel summation, DLMF 9.8.20, and `r,s>=r_min` then give

```text
sup_(contiguous I) |sum_(j in I) K(r_j,s)|
 <[0.02449652675353782363042710617241782344944237247177813976809721694544243513709618715297550801057279523 +/- 2.20e-102]<0.0246.       (GP7)
```

The corresponding termwise pointwise estimate is
`[0.1861265712877299879365740442768799085505258268946836830994676912653377138667936218555350139523900072 +/- 1.47e-101]`; (GP7) improves it by a factor
greater than `[7.598080052750715804234399414038989177536954052001776448349663121461922363863416497279888696355786313 +/- 1.08e-101]` while
retaining phase cancellation.

Equation (GP7) is the kernel packet required after interchanging the finite
event sum with a variation-of-constants integral: at fixed `s`, the active
event indices are contiguous.

Returning from `H(r_m)` to the original transform introduces no modewise
gauge defect.  Indeed, on the exact event lattice,

```text
chi(d_m)=2beta^2 d_m+beta d_m^2
        =2pi m^2-beta^3,
exp[-i chi(d_m)]=exp(i beta^3).                      (GP8)
```

The last factor is common to every integer mode.  Thus the grouped packet
survives the inverse gauge exactly rather than approximately.

The initial-data channel uses

```text
partial_s K(r_j,s)=pi M(-r_j)N(-s)
                    sin(varphi(-s)-vartheta(r_j)).                    (GP9)
```

The event phasor is unchanged.  DLMF 9.8.21 and its signed remainder give

```text
N(-s)^2<=sqrt(s)/pi [1+7/(32s^3)].                  (GP10)
```

The same Abel argument therefore proves

```text
sup_(contiguous I) |sum_(j in I) partial_s K(r_j,s)|
 <[52.77586540497248519805067055145587787983725349165605435454423316906751128487475413576341585572092159 +/- 3.79e-99]<52.8. (GP11)
```

Its termwise counterpart is
`[400.9952502001780176338344059518995644552024262966923969260235076445923438704837765970629546756151092 +/- 1.59e-98]`, again improving by a
factor greater than
`[7.598080052750715804234399414038989177536954052001776448349663121461922363863416497279888693245078356 +/- 3.49e-100]`.

Equations (GP7) and (GP11) do not yet compose the initial data, coherent
endpoint sources, and interior forcing with the completed Poisson remainder.

Pi provenance: the phase integer and `P_5` inherit
`beta^3=pi C^2/8`; (GP6) uses the Airy Wronskian normalization.  The Airy
phase, modulus monotonicity, and modulus remainder are sourced from
`https://dlmf.nist.gov/9.8`.

Machine-audited companion:

```text
outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_detuning_airy_lattice_grouped_green_packet_gate.md
work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_detuning_airy_lattice_grouped_green_packet_gate.json
work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_detuning_airy_lattice_grouped_green_packet_gate.py
work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_detuning_airy_lattice_grouped_green_packet_gate.py
```

No completed endpoint cancellation, grouped 398-mode finite-integral splice,
all-corridor continuation, complete
`Q_K-T` or `T_upper`, `Lambda<=0`, PF-infinity, RH, or prize-level
conclusion is proved.
