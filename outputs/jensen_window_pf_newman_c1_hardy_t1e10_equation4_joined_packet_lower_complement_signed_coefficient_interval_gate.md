# Signed-coefficient interval enclosure of the lower complement

Date: 2026-08-27

Status: the lower complementary tail and complete ordinary packet at the
production height are rigorously enclosed as complex balls; the later full
joined packet remains open; this is not a proof of RH.

The old lower-tail remainder bound expanded every order-eight exact recurrence
coefficient into absolute monomials.  At the limiting endpoint `a`, that
inflates the signed coefficient by factors ranging from
`[1.00000000000000000000000000000000000000000000000000000000000000000000000 +/- 1e-76]` through
`[695088449674257650222135.097467488883418401092982283812563737901644184576 +/- 1.57e-49]`.  The largest losses
occur in denominator powers 13 through 16, exactly where the recurrence
contains strong internal cancellation.

Keep the exact collected coefficient `c_k(y)` for each denominator power and,
on every slab `Y_j`, use the rigorous interval bound

```text
sup_(y in Y_j)|c_k(y)|
 <= abs(c_k(Y_j)).upper.                            (LC1)
```

Only after (LC1) is formed is it multiplied by the positive denominator sum

```text
sum_(n>=0)(delta+n)^(-k)
 <= delta^(-k)+delta^(1-k)/(k-1).                  (LC2)
```

This preserves cancellation internal to each exact coefficient without
assuming cancellation between labels, powers, or slabs.  A rigorous sweep of
the first admissible orders gives

```text
order 2: [2.62777802742344036086419086673865052273724771547004431652154802880403071e-6 +/- 3.94e-78],
order 3: [3.43144948472073618862619562163548514629934634706651380720835718182924818e-7 +/- 2.86e-79],
order 4: [6.02462994098044092064092743214657478036969989005920573370157546178897081e-7 +/- 7.79e-80],
order 5: [1.49132680612525399842894560279749019401125619196776010337220590889835212e-6 +/- 3.89e-78],
order 6: [4.68310873602827019818321894108886990171282796946512884730997519104189044e-6 +/- 2.38e-78],
order 7: [1.77632101539458774486678843030058601650705156131208631290689282301408476e-5 +/- 4.97e-77],
order 8: [7.88577110995210870959603666152270507270145386120493441544005796610750731e-5 +/- 2.44e-77],
order 9: [0.000400749710136320382606525321517097665302117153059001169153781698970809453 +/- 4.14e-77].
```

Order one is excluded because its integrated infinite label sum still has a
denominator-power-one term and is not absolutely summable.  Thus order three
is the certified minimum over orders 2 through 9 under this slab rule.  Using
16384 fourth-power right-endpoint-clustered slabs gives

```text
T_lower endpoint partial =
  [-0.000373272059958648030108458493369361471211836854153018044088714233530343776 +/- 2.82e-76]
 +i[-0.000407687677620205622371538849197569061175897750341206654907622203537301313 +/- 5.92e-76],

T_lower remainder <= [3.43144948472073618862619562163548514629934634706651380720835718182924818e-7 +/- 2.86e-79],

T_lower = [-0.000373272059958648030108458493369361471211836854153018044088714233530343776 +/- 3.44e-7]
         +i[-0.000407687677620205622371538849197569061175897750341206654907622203537301313 +/- 3.44e-7].                    (LC3)
```

The former monomial-triangle remainder bound was
`[25.6198115131043205624420429553941813562250692480943562746261504120688319 +/- 4.17e-71]`; the certified
improvement factor is `[74661776.6841156159246699208345401277173171270899790317993508152656181627 +/- 1.49e-64]`.

The exact complementary identity is

```text
O_join=-T_lower-T_upper
       +(1-C_G)sum_(m=622)^39852 m^(-s).            (LC4)
```

The last packet has log10 absolute upper bound
`[-27287527074.5945097191391910541815681408865708525983385791579325976981661 +/- 2.00e-62]`.  Joining it as a
complex error ball to the certified complete upper complement and (LC3)
before taking a norm gives

```text
O_join = [0.000185680827068395015533301275967734129786506502431127354479477057796135678 +/- 3.44e-7]
        +i[0.000135613304173880147425441503047167047268537480552296189063588367148358990 +/- 3.44e-7],

|O_join| = [0.000229930657610566413495689630508422851562500000000000000000000000000000000 +/- 4.81e-7].             (LC5)
```

The independent checker changes to four recurrence rounds, 20000
fifth-power-clustered slabs, and explicit alternating eta endpoint values.
Its lower-tail and ordinary-packet balls overlap (LC3) and (LC5).  A separate
altered-height negative-denominator direct quadrature overlaps the corresponding
fifth-order recurrence enclosure.

Pi provenance: every `pi` comes from the inherited Fresnel/Fourier phase,
`phi'=2*pi*D`, alternating Fourier characters, or `p=t/(2*pi)`.  The signed
coefficient collection follows exact differentiation; no circle, polygon,
fitted constant, or visual pattern supplies `pi`.

Proof boundary: rigorous signed-coefficient order-three lower-tail enclosure
and the resulting complete production-height ordinary-packet complex ball
only.  No complete later joined packet, `J_Z`, `D_K`, non-A, all-height,
`Lambda<=0`, PF-infinity, RH, or prize-level conclusion is proved.
