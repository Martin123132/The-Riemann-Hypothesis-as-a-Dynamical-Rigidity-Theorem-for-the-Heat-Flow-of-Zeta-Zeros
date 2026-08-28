# Sparse source-minus-owned-carriers altered-representation pilot

Date: 2026-08-26

Status: **diagnostic sparse pilot passed; interval physical quadrature remains open**.

## Common mathematical level

This gate evaluates only the pointwise source/ownership block from Core 11.476:

```text
S(x)=sum_(n=0)^L f_x(n)
     -sum_(m=622)^39936 P_m(x)
     -sum_(m=39853)^39936[A_m(x)+A_-m(x)].
```

The saved B-trace, B-outer, and translated A-face artifacts are physical integrals or bounds.
They are deliberately not inserted as pointwise constants.  Their ownership is unchanged.

The primary evaluator joins `P_m+A_m+A_-m` on the 84 A-window modes before
summation.  Selected rows are independently replayed as `-I_m` plus the original
ordinary and A-window endpoint rows.

## Results

- Rational sample rows: `18`.
- Full 39,315-mode altered-form rows: `8`.
- Maximum altered-form relative discrepancy: `5.75208859465524140478933126813786847309831766e-79`.
- Periodic full-source oracle rows: `9`.
- Maximum source-mass-normalized recurrence error: `0.00000000000135315556224451746179354840753824295496510883`.
- Cancellation-ratio range: `0.0431721078725745382502609929896358359529582801` to `0.999998161218292762055707141821646898051060331`.
- Maximum sampled normalized physical-density magnitude: `169212509674.38742390643360708313898107868137`.

| label | x | nearest transition | cancellation ratio | physical density |
|---|---:|---|---:|---:|
| first_B_entry | 1244/5122421 | B_entry 622 (distance 0) | 0.999903179983916320127814855408056007409522872 | 169212509674.38742390643360708313898107868137 |
| dyadic_1_over_4096 | 1/4096 | B_entry 625 (distance 2421/20981436416) | 0.999998161218292762055707141821646898051060331 | -304991574.046641605063574315923188728421282601 |
| decimal_1_over_1000 | 1/1000 | B_entry 2561 (distance 421/5122421000) | 0.99996628942597066912434840773181834193905638 | -6702691.21047449158494721966643996736596525845 |
| first_A_exit | 1244/159577 | A_exit 622 (distance 0) | 0.604312411072885191546975886668233088530442117 | -580752712.470222420368935088532016776625161782 |
| closest_pair_A_630 | 1260/159577 | A_exit 630 (distance 0) | 0.592431972124839288497328694991663413821281107 | -261429100.622114317770295963460320163723300302 |
| closest_pair_B_20223 | 40446/5122421 | B_entry 20223 (distance 0) | 0.674802312466465725302234204195236340512290618 | 108172516.880342445884021090517894510577213272 |
| last_B_entry | 79872/5122421 | B_entry 39936 (distance 0) | 0.0431721078725745382502609929896358359529582801 | 345947.003913671165172650875811347769902893093 |
| dyadic_1_over_64 | 1/64 | A_exit 1247 (distance 39/10212928) | 0.202612669701924852626446121286785805479473318 | 44110.5451361209061497226635721221402828940272 |
| dyadic_1_over_32 | 1/32 | A_exit 2493 (distance 25/5106464) | 0.999985388710887528455598499480754676669522241 | 180198.933756138143639521862667163570068933479 |
| dyadic_1_over_8 | 1/8 | A_exit 9974 (distance 7/1276616) | 0.996668265320141763928626209327849811323252992 | 9941.42511042039754738056433116151978953150978 |
| quarter | 1/4 | A_exit 19947 (distance 1/638308) | 0.891300282873209406386282868210325130965083844 | -18813.1390455771562800424437554558071942810739 |
| third | 1/3 | A_exit 26596 (distance 1/478731) | 0.995611602799610412246465344046648052439659182 | -18279191126.0314309849600228986619477696573852 |
| two_fifths | 2/5 | A_exit 31915 (distance 4/797885) | 0.996950456420534788025334954243482383466288502 | 2680217720.29723527978312600997828481493597027 |
| seven_sixteenths | 7/16 | A_exit 34907 (distance 15/2553232) | 0.962220139666772553164230251904255854169416812 | -53179.5885802587522036691990615589132206738944 |
| last_interior_A_exit | 79704/159577 | A_exit 39852 (distance 0) | 0.996765957548011009354337787489889543241972123 | 12705814.2471627024666361027386803481266123555 |
| first_A_window_exit | 7246/14507 | A_exit 39853 (distance 0) | 0.998057785832153593058409082075934845344346184 | 268790324.625709368563906603465578400635972249 |
| last_lower_A_exit | 79788/159577 | A_exit 39894 (distance 0) | 0.994286871863106373654664734163374225330115176 | 81107788.40987586498693284052523975882082472 |
| midpoint | 1/2 | A_exit 39894 (distance 1/319154) | 0.995163662081068808436607128846196720022683253 | -26299.6738876745155368569394274927196554788197 |

## Decision

The joined pointwise representation is numerically coherent on the sparse roster and is the
selected evaluator for the next bounded integration experiment.  The raw endpoint-plus-Fresnel
form remains a checker only.  This gate does not turn the sparse samples or the asymptotic
Fresnel-tail diagnostics into an interval error theorem.

The next obligation is a physical-transform ownership ledger for the complete source, finite
`P/A` carriers, and unchanged physical channels.  Reuse certified transformed integrals where
they exist; launch event-aware panels only for a term that has no such transform.  Only a joined
interval result below the inherited non-A allowance may be promoted.

## Pi provenance

Every `pi` comes from the original quadratic Kummer character, its integer Fourier transform,
the canonical Fresnel primitive, or the inherited equation-(9) physical normalization.
No fitted or geometric occurrence is introduced.

## Boundary

A deterministic floating sparse-point pilot at t=10^10 for the complete source minus its owned positive full-line P block and 84 paired A atoms, including rational-phase block reseeding, exact transition labels, selected full-roster periodic source checks, and altered endpoint-plus-Fresnel replays. The Fresnel asymptotic stopping diagnostic is not a rigorous remainder enclosure, the O(L) source recurrence has no uniform floating error theorem, and no interval x quadrature is performed. The B-trace, B-outer, and translated A-face physical channels are not reevaluated or mixed into this pointwise gate. No non-A bound, joined R_after_A, R_Dir, Q_K-T, all-height theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.
