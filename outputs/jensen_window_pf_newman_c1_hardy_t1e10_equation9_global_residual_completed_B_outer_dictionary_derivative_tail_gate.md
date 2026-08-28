# Outer Fresnel-to-boundary dictionary derivative tail

Date: 2026-08-13

Status: uniform analytic dictionary-tail derivative certificate; not a proof
of the complete completed-B estimate

The direct three-term Fresnel expansion and the first three endpoint
boundary currents differ exactly by

```text
D_m=48x^2(B^4x^4+40B^2m^2x^2+80m^4)
    /[pi^4(Bx-2m)^5(Bx+2m)^5].                       (OD1)
```

With `c=Bx/2`, this is

```text
D_m=3c^2(c^4+10m^2c^2+5m^4)
    /[pi^4 B^2(c^2-m^2)^5].                          (OD2)
```

Every monomial in (OD2) is differentiated through order two before norms.
For `m>=B=5122421` and `1e-4<=x<=1/2`, the same
`|c^2-m^2|>=(15/16)m^2` margin and decreasing-series integral test give

```text
sum_(m>=B)|D_m|    <= [8.212668816256876343321415225064040839372463202464211845474144413857878704801571380062727194834765362e-37 +/- 1.55e-136]
sum_(m>=B)|D_m'|   <= [4.344746124475119055731245787343104307014946408689169601032794340390783788984966551696968044838424102e-36 +/- 9.59e-136]
sum_(m>=B)|D_m''|  <= [1.872640700409613832150147173220218014036334238381288114730541910518965433562763374461301616193769872e-35 +/- 6.57e-135].                     (OD3)
```

After the Kummer weight and normalized coordinate,

```text
|w D|          <= [8.212874145810540290961717770075010149662286405265759720632984780337009229980907371424348520381427542e-36 +/- 1.94e-135]
|d_xi(w D)|    <= [7.060897103113070555046108462679378766895409654069237947906255327381557265653729162398176587882636891e-41 +/- 1.59e-138]
|d_xixi(w D)|  <= [3.025488955577184909564061672028677733752852432772210044117984410559645610031339725052634619810541195e-45 +/- 1.36e-142].                        (OD4)
```

Together with the rational-current and exact Fresnel-remainder gates, this
closes the complete analytic outer block `m>=B` through two normalized
derivatives.  The finite block `m<B`, phase-coupled exterior integration,
all indicator/cutoff boundary terms, and the common-regulator joined
remainder remain open.  No complete B estimate, `T_upper`, height-uniform
theorem, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion is proved.

Pi provenance: `pi` is inherited from the exact equation-(9) Fresnel phase
and the endpoint integration-by-parts dictionary; it is not inserted as a
geometric or fitted parameter.
