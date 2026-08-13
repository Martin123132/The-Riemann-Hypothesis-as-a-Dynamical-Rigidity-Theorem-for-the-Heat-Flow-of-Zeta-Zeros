# Hardy block-20 direct MIT=1 kernel gate

Date: 2026-08-06

Status: rigorous finite exact-point audit of all 50 direct kernels; not a proof or outer Hardy theorem

## Direct route

For `MIT=1` the source performs no recurrence and no `q` correction.  It sums

```text
S_hat = sum_(n=0)^104 exp(i*tpm_hat*(a1*n+a2*n^2+a3*n^3)),
```

then applies only the saved conjugation orientation.  All 50 scattered direct
calls are evaluated from their exact binary128 coefficients at 180 and 260
decimal digits.  Every low/high precision pair overlaps.  The saved direct
state was printed with 40 digits after the decimal rather than a hex payload;
each component is therefore enclosed by one full last-decimal unit before
comparison.

The source constant has explicit provenance: `p=4*atan(1)`, `tpp=2*p`, and
`tpm_hat=-tpp`.  The gate separately evaluates the mathematical normalization
`tpm=-2*pi` using Arb's standard pi ball.  No arbitrary circle or polygon is
introduced.

## Finite bounds

```text
maximum exact-source expression versus logged state <= 2.63352964476090825300437347481675978884857497463545E-31
maximum true-2*pi versus source-tpm kernel shift     <= 3.54697478609679005770998174625554860905033211471966E-32
maximum true-2*pi versus logged state gap            <= 2.48161337417395396501533018015609985611061762259977E-31
maximum analytic phase-Lipschitz budget              <= 1.86614931457278521250381799656523055958801692938774E-31
minimum rigorous true-kernel magnitude               >= 4.30886767593764727740389592382671289585277826839587E+0
```

For real phases, `|exp(iu)-exp(iv)|<=|u-v|`; summing this inequality gives the
recorded phase-Lipschitz budget and independently dominates every observed
normalization shift.  The direct source computation is therefore enclosed at
all 50 saved points without borrowing the recursive `q` model.

This does not yet transport coefficient errors from a continuous physical
height interval, bound the outer Hardy representation or the cross-call block
sum, or address the remaining recursive `t5` saddle error.  It proves no
height-uniform theorem, determinant/current sign, `Lambda<=0`, PF-infinity, RH,
or prize-level conclusion.
