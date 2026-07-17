# Order-nine localized lower-bridge certificate

Date: 2026-07-14

Status: rigorous first-summand curvature theorem on `1250<=t<=5700`; This is not a proof of RH, `Lambda <= 0`, PF-infinity, or even the
full order-nine entry theorem by itself.

## Certified statement

The hash-bound interval computation proves

`w_1''(t)<=4200/t^2 for every real 1250<=t<=5700`.

It uses exact-t `H2-H22` collars, exact-point `H0-H8` jets, a strict
`t+-7` support collar, positive `B,J,R,S,T,U,V`, a sixth-order point
Taylor polynomial, and a localized seventh-order remainder.

## Coverage

- Root segments: 279
- Accepted blocks: 874
- Accepted widths: `1/2`: 152, `1`: 150, `2`: 140, `4`: 140, `8`: 161, `16`: 131
- Largest scaled curvature upper: 4.19424425522037111044561248832644227467685288026319528153785E+3<4200
- Weakest common coordinate lower: 0.634702287167133801127866466006058615418192489161371203429560 (B)

The adaptive algorithm starts with exact 16-unit rational segments and
bisects failed interval enclosures dyadically. Every accepted block is
checked independently, and the accepted endpoints join exactly from
`t=1250` through `t=5700`.

## Boundary

Separate compact, finite-ray, and asymptotic-ray artifacts cover the
first-summand range above `t=5700`. Their composition and the existing
full-kernel transfer still require an explicit umbrella certificate.
