# Order-Ten Lambda-Zero Prefix Certificate

Date: 2026-07-16

Status: rigorous finite certificate for the four lambda-zero prefix rows. This is not a proof of RH or `Lambda <= 0`.

```text
python work/rh_compute/scripts/jensen_window_pf_compound_order10_lambda0_prefix_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_compound_order10_lambda0_prefix_certificate.py
```

## Index And Normalization

The coefficient enclosure source uses
`A_k(lambda)=mu_(2k)(lambda) k!/(2k)!`. The project coordinate is

```text
H_(m,n)=det[A_(n+i+j)]_(0<=i,j<m)
epsilon_m=(-1)^binom(m,2)
Q_(m,n)=epsilon_m H_(m,n).
```

The source probe parameter `m=9` builds a `10x10` matrix and uses
`sigma=(-1)^(9*10/2)=-1=epsilon_10`. Thus probe shift `n` is exactly
`Q_(10,n)`, with no rescaling or sign conversion left over.

## Certified Rows

Rebuilding each determinant from the rigorous lambda-zero `A_ball` inputs
at 520 decimal digits, and independently overlapping the pre-existing
signed-Hankel row, proves

```text
Q_(10,n)(0)>0 for every integer 0<=n<=3.
```

| n | rigorous Q_(10,n)(0) enclosure |
|---:|---|
| 0 | `[6.122601728800756973263039882428856342279754341050629910427e-271 +/- 7.74e-329]` |
| 1 | `[4.85410207810852138935088900380689603260497548726993815580e-296 +/- 7.03e-353]` |
| 2 | `[3.2184750025458154961067928863145107326011517462970400960e-321 +/- 5.48e-377]` |
| 3 | `[1.794614559198768285668870688859631526278643197408175440e-346 +/- 1.62e-401]` |

## Boundary

This proves four contiguous signed Hankel rows at lambda=0 only. It does not prove positivity at negative lambda, the order-ten all-shift theorem until the delayed heat tail is composed, PF-infinity, RH, or Lambda<=0.
In particular, this finite certificate alone is not an all-shift
order-ten theorem and is not a proof of PF-infinity, RH, or `Lambda<=0`.
