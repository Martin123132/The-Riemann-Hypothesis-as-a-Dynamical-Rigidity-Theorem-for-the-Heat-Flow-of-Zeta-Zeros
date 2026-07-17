# Order-Eleven Lambda-Zero Prefix Certificate

Date: 2026-07-16

Status: rigorous certificate for the four shifted-flow complement rows. This is not a proof of RH or `Lambda <= 0`.

```text
python work/rh_compute/scripts/jensen_window_pf_compound_order11_lambda0_prefix_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_compound_order11_lambda0_prefix_certificate.py
```

The source probe parameter `m=10` builds an `11x11` determinant and
uses `sigma=(-1)^(10*11/2)=-1=epsilon_11`. Therefore its shift `n`
is exactly the project coordinate `Q_(11,n)`.

Fresh determinant rebuilds from the rigorous 520-digit coefficient
balls independently overlap all four source rows and prove

```text
Q_(11,n)(0)>0 for every integer 0<=n<=3.
```

| n | rigorous Q_(11,n)(0) enclosure |
|---:|---|
| 0 | `[1.344399299541183706420982662098361242815521675821418372e-330 +/- 5.32e-385]` |
| 1 | `[2.23932940997055059825565619534397912680695701602270644e-358 +/- 6.93e-412]` |
| 2 | `[3.1047661733895941602263053683075325795179262585774414e-386 +/- 8.18e-439]` |
| 3 | `[3.600881351231391643426044110321070168597043930391776e-414 +/- 8.68e-466]` |

This prefix is ready to compose with a future shifted heat ray on
`n>=4`. It does not itself prove that ray, order twelve,
PF-infinity, RH, or `Lambda<=0`.
