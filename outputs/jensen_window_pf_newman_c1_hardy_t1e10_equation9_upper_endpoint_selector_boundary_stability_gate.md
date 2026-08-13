# Upper-endpoint stability at the odd-selector boundary

Date: 2026-08-12

Status: interval certificate for the pinned Hardy source endpoint recurrence;
not a proof of the complete source transition splice or RH

The saved output at `t=10^10` starts Block 0 at the effective odd integer
`159577` and ends Block 35 at

```text
B_saved=5122421.                              (UE1)
```

The endpoint is not chosen by a separate cutoff.  After Block 0 the pinned
source iterates

```text
G_j=floor(X^j ib/2)                       (inner branch),
M_j=the source floor/cap/oddized Gauss-sum length,
B_j=B_(j-1)+2(M_j+1)G_j.                  (UE2)
```

The extracted recurrence reproduces all
`36` printed endpoints,
including (UE1), exactly.

Let

```text
t*=pi*159577^2/8
  =[10000011009.04258833249996651073906089290177178617834443516760678611578387718172618748916453130539395 +/- 1.77e-90].
```

The source `start` routine uses the unsuffixed default-real test `g<3.2`.
Under the admitted build that threshold is

```text
3.200000047683715820312500000000000000000000000000000000000000000000000000000000000000000000000000000.                  (UE3)
```

As the height approaches `t*` from below, this transition advances the
effective start from `159577` to `159579` before
the exact odd-selector boundary.  Its entry height is enclosed by

```text
[10000002368.45285511579661582259313224965527527695873791980491842836971000000000000000000000000000000 +/- 1.01e-24],
```

which is

```text
[8640.589733216703350688145928643246496509219606515362688357746073877181726187489164531305393951659410 +/- 1.01e-24]
```

below `t*`.  At and just above `t*`, the ordinary next-odd rule already gives
the same effective start `159579`.  Consequently no start
index changes at `t*` itself.

Every automatic block count, inner/outer branch, Gauss-sum count, length
floor, cap, and oddization is constant on the Arb-certified interval

```text
t*-8000
 <= t <=
t*+200000.                    (UE4)
```

The conservative minimum raw floor slack is greater than
`4.768371582031250000000e-6`, and the conservative
minimum branch slack is greater than
`[602.7579565380141139030 +/- 4.57e-20]`.  Both the source-literal
and intended-real recurrences produce the same 35 rows.  Their common endpoint
is

```text
B_boundary=5122423.                         (UE5)
```

Thus `B` changed by `+2` during the earlier transition layer, but it does not
change simultaneously with `C -> C+2` at `t*`.  The next analogous transition
begins `[242023.8743084176760132786984325396873464534680666116322316665273561228182738125108354686946060483991 +/- 1.01e-24]`
above `t*`.

For the exact fixed-`B` selector identity at the boundary, using
`B=5122423` gives `2481424`
formal old-selector lattice terms and `2481423`
new-selector terms, again differing by one.  Hence there is no additional
upper strip at the selector height.

Proof boundary: this closes only the question whether the pinned upper
endpoint jumps at the same odd-selector boundary.  It does not yet reassemble
the earlier coupled event consisting of the lower transition term, the direct
roster shift, and the added upper lattice point.  It also does not prove the
399-mode quantitative remainder, complete `T_upper`, `Lambda<=0`, RH, or a
prize-level conclusion.
