# Q207-Q208 Parabolic-Frequency Relative Diagnostic

Date: 2026-07-26

Status: rigorous finite source-reusing relative-transport diagnostic on the
first proved forward successor. This is not an all-stage estimate and not a
proof of `Lambda<=0` or RH.

## Source And Scope

The diagnostic reuses the 525 already certified transport panels covering

```text
[1/1040,1/1035]x[38,245].
```

It performs no new Xi quadrature. The old-edge boxes are for
`F=16(1+x^4)H`, and the stored collar bounds enclose `F_t` and `F_(xt)`.

For a positive scale `s`, put `V_s=(F,sF_x)`. On each panel the stored
componentwise transport bounds give rigorous collar boxes for `F` and
`F_x`, followed by

```text
||partial_t V_s||
 <=sqrt(|F_t|^2+(|s_t||F_x|+s|F_(xt)|)^2),

||V_s||
 >=sqrt(dist(F,0)^2+s_min^2 dist(F_x,0)^2).         (1)
```

Every denominator in (1) is strictly positive.

## Four Scales

```text
unit:       s=1,
frequency:  s=1/L,
parabolic:  s=sqrt(2t),
blended:    s=(L^2+(2t)^(-1))^(-1/2).
```

| scale | max one-step relative cost | max panel | min collar denominator | min panel |
|---|---:|---|---:|---|
| unit | [0.9816633312152211363702219024287642861006044804148409821 +/- 2.91e-56] | [189,379/2] | [4.331444754571908785199800032655506324949273486329727613e-29 +/- 4.02e-85] | [243,973/4] |
| frequency | [2.565613048400747770944442646965811362773969906337190560 +/- 1.89e-56] | [189,379/2] | [3.301092246150930362897456218821063501651460597751952670e-29 +/- 1.24e-84] | [971/4,243] |
| parabolic | [21.43114290405414035253495699725338312077820533849129874 +/- 4.99e-54] | [189,379/2] | [4.287918247712455415051350881557042361542665609843234279e-30 +/- 4.25e-85] | [971/4,243] |
| blended | [21.58229524918177011395351792743838285843242944277557260 +/- 3.77e-54] | [189,379/2] | [4.252195856780897894007187342775197064887727674993792345e-30 +/- 2.78e-85] | [971/4,243] |

The ratio of the blended maximum cost to the parabolic maximum cost is
`[1.007052929738947155933541063048215147959920385265007469 +/- 3.37e-56]`. Its ratio to the
frequency maximum cost is
`[8.412139649287683191042187653834550732903667206940791386 +/- 4.08e-55]`.

## Interpretation

At Q207-Q208, `2tL^2` is small, so the blended scale should be close to the
parabolic scale. This finite audit checks that expectation using rigorous
old-cell and heat-jet enclosures. A finite one-step relative certificate is
not an asymptotic Xi bound; it only tests the conditioning of the proposed
coordinate on the first proved successor.

## Proof Boundary

All 525 finite panel inequalities and all four scale
comparisons are rigorous consequences of stored certified interval data.
They prove no Q209 stage, no uniform old-collar relative theorem, no
ray-aligned all-stage successor, no `Lambda<=0`, no RH, and no Clay-prize
result.

Machine-audited files:

```text
work/rh_compute/results/jensen_window_pf_newman_q207_q208_parabolic_frequency_relative_diagnostic.json
work/rh_compute/scripts/jensen_window_pf_newman_q207_q208_parabolic_frequency_relative_diagnostic.py
work/rh_compute/scripts/check_jensen_window_pf_newman_q207_q208_parabolic_frequency_relative_diagnostic.py
```
