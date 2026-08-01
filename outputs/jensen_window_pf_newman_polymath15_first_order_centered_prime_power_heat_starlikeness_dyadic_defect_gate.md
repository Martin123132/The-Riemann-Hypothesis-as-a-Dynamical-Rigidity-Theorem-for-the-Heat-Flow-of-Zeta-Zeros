# Dyadic All-Length Fourier-Defect Localization Gate

Date: 2026-07-30

Status: exact defect-localization theorem for p=2,M>=10.
This does not yet prove dyadic heat-starlikeness or RH.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_first_order_centered_prime_power_heat_starlikeness_dyadic_defect_gate.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_prime_power_heat_starlikeness_dyadic_defect_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_prime_power_heat_starlikeness_dyadic_defect_gate.py
```

## Outcome

```text
p=2,M>=10: D_d>0 for 1<=d<=M-3; D_(M-1)>0
possible negative indices: D_0 and D_(M-2)
```

The old numerical pattern is now an all-length theorem. The
zero-mode defect is genuine, so the remaining step needs a
two-defect positive-kernel or equivalent phase-aware argument.

## Terminal Slope

For the last common group and its two boundary terms,

```text
C(B,d)=B-2(B+1)aq^3+(B+2)a^2q^4+a^2q^(2d+6)[B+2-2(B+3)aq+(B+4)a^2q^(2d+2)]
E_d=1-2aq^3+a^2q^4+a^2q^(2d+6)[1-2aq+a^2q^(2d+2)]
E_d>=(1-a)^2-(2a-1)_+^2/4
E_d>(3-2sqrt(2))/4>1/25
```

Thus `C(B,d)` increases with the actual length parameter `B`.

## Low Offsets

| certificate | lower | shifted Bernstein minimum | degree |
|---|---:|---:|---:|
| d=1 first seven | 1/200 | 4.92967806836108245485e-3 | 76 x 14 |
| d=1, M=10 full | 1/100 | 8.59843886491353405498e-3 | 144 x 16 |
| d=2 first three | 1/25 | 7.20348340752033404922e-3 | 10 x 6 |
| C(15,1) | 1/1 | 3.94658523510873866988e-1 | 12 x 4 |
| C(14,2) | 1/1 | 2.65979032790886194305e-1 | 16 x 4 |

The `d=1` first-seven certificate is used for `M>=11`; the
separate normalized certificate handles `M=10`. The `d=2`
three-group cluster is available for every `M>=10`.

## Middle And Tail

| d | lower | shifted Bernstein minimum | degree |
|---:|---:|---:|---:|
| 3 | 1/25 | 5.98654127689683426806e-2 | 20 x 4 |
| 4 | 1/25 | 1.56346387240449708514e-1 | 24 x 4 |
| 5 | 1/25 | 2.38839087336521810875e-1 | 28 x 4 |
| 6 | 1/25 | 3.13693924044122484606e-1 | 32 x 4 |
| 7 | 1/25 | 3.82620457535385248826e-1 | 36 x 4 |
| 8 | 1/25 | 4.47817532235762768400e-1 | 40 x 4 |
| 9 | 1/25 | 5.10098209099155503344e-1 | 44 x 4 |
| 10 | 1/25 | 5.70189934929834104679e-1 | 48 x 4 |
| 11 | 1/25 | 6.27986499990405112441e-1 | 52 x 4 |
| 12 | 1/25 | 6.84074387808541248425e-1 | 56 x 4 |
| 13 | 1/25 | 7.39446468517832414765e-1 | 60 x 4 |
| 14 | 1/25 | 7.92879956295722161104e-1 | 64 x 4 |

For `d>=15`, no finite scan is used:

```text
C(d+2,d)=dE_d+K_d
delta=(3-2sqrt(2))/4
kappa=5/4-(4/3)sqrt(2)
15delta+kappa=(75-53sqrt(2))/6>0
75^2-2*53^2=7>0
```

## Boundary

The exact common-group decomposition, dyadic cluster-slope bound, low-offset rational certificates, and analytic terminal tail prove D_d>0 for every 1<=d<=M-3 and D_(M-1)>0 when p=2,M>=10. Exactly D_0 and D_(M-2) remain as possible Fourier defects. This gate does not yet prove uniform dyadic heat-starlikeness, an actual-ray theorem, a global join, Lambda<=0, PF-infinity, RH, or a prize-level result.
