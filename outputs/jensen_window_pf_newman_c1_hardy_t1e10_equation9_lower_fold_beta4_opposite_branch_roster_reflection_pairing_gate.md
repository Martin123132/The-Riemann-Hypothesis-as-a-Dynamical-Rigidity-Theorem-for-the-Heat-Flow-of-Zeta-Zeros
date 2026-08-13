# Beta^-4 opposite-branch roster reflection pairing

Date: 2026-08-13
Status: exact finite-roster symmetry; not a proof of the grouped integral

Put `a=beta/C`, `x=a*y`, and write the exact branches from Sections
11.378--11.380 as

```text
C_+(y)=exp(i[xi-pi/4])S_+(X,y),
C_-(y)=conj(C_+(y)).                                  (RP1)
```

For the fixed selector `C=159577`, the 399 modes `39695..40093` have
detuning labels `4m-C=-797,-793,...,795`.  There are 200 negative modes and
199 positive modes.  The branch opposite to the endpoint crossing is `C_+`
on the negative half and `C_-` on the positive half, so its complete roster
sum is

```text
O_C=C_+*sum_(j=0)^199 exp(i*(4j+1)*a*y)+C_-*sum_(j=0)^198 exp(-i*(4j+3)*a*y).                        (RP2)
```

Pair the labels `-(4j+1)` and `+(4j+3)` for `0<=j<=198`.  The only unpaired
old-selector label is `-797`, belonging to mode `39695`.  Exact geometric
summation gives

```text
O_C=C_+*exp(i*399*x)*sin(400*x)/sin(2*x)+C_-*exp(-i*399*x)*sin(398*x)/sin(2*x), x=a*y,                           (RP3)

O_C=exp(i*797*x)C_+ + 2*exp(-i*x)*sin(398*x)/sin(2*x)*Re(exp(i*398*x)C_+).                      (RP4)
```

All quotients in (RP3)--(RP4) are removable: at `sin(2x)=0` their values are
defined by the finite sums in (RP2).  Thus (RP4) is 199 conjugate pairs plus
one edge mode, equivalently one real Dirichlet projection plus that edge,
not 399 unrelated branch errors.

At the shared selector boundary, the already-certified one-cell Fourier block
has modes `39696..40094` but still uses the old strip coordinate `C`.
Its detuning labels are therefore `-793,-789,...,799`, not the adjacent-chart
labels.  Its exact reflected form is

```text
O_F=C_+*exp(i*397*x)*sin(398*x)/sin(2*x)+C_-*exp(-i*401*x)*sin(400*x)/sin(2*x),                     (RP5)

O_F=exp(-i*799*x)C_- + exp(-i*x)*sin(398*x)/sin(2*x)*(exp(i*398*x)C_+ + exp(-i*398*x)C_-).                     (RP6)
```

This block has 199 negative modes, 200 positive modes, and unpaired label
`+799`, mode `40094`.  Set-theoretically, passing from the old lower-edge
roster to this boundary block keeps `39696..40093`, drops `39695`,
and adds `40094`.  This is only a finite mode-set identity; the complete
selector jump also contains the symmetric Poisson complement and endpoint
half-current already certified by the Fourier-completion gate.

The adjacent selector `C+2=159579` has the same mode set but its own event
coordinate, with labels `-795,-791,...,797`.  Its separate mirror formula is

```text
O_(C+2)=C_+'*exp(i*399*x')*sin(398*x')/sin(2*x')+C_-'*exp(-i*399*x')*sin(400*x')/sin(2*x'),                      (RP7)

O_(C+2)=exp(-i*797*x')C_-' + exp(i*x')*sin(398*x')/sin(2*x')*(exp(i*398*x')C_+' + exp(-i*398*x')C_-').                      (RP8)
```

The adjacent unpaired label is `+797`, again mode `40094`.  At `y=0`, the
old lower-edge and shared-boundary projections have opposite one-copy
imaginary imbalances,

```text
O_C(0)=200*C_+ + 199*C_-=399*Re(C_+)+i*Im(C_+),
O_F(0)=199*C_+ + 200*C_-=399*Re(C_+)-i*Im(C_+).                       (RP9)
```

The two canonical charts do not have identical scales.  Direct interval
evaluation gives

```text
a_C=[0.01350091479749819087428255402914909480775074027248099775869154129411374374815152748365494323780473255 +/- 4.32e-102],
a_(C+2)=[0.01350085839504256237300864996144539370895647131059895015102550537810964716026701060879746986015461194 +/- 4.43e-102],
|a_(C+2)-a_C|=[5.640245562850127390406770370109879426896188204760766603591600409658788451687485747337765012060989102e-8 +/- 6.57e-102]<6e-8. (RP10)
```

The boundary Fourier-completion theorem already combines its whole 399 block,
the zero/negative/outer-positive complement, and the endpoint half-current
before absolute values.  Therefore (RP6) is an internal branch organization,
not a license to replace that completed identity by an edge difference.  The
next admissible use is to insert (RP6) into the completed object and determine
which branch projection continues into the ordinary-Morse corridor.  A
termwise absolute sum over 399 modes would erase this exact structure.

This gate proves the three roster reflection algebras only.  It does not bound
the grouped opposite-branch integral, identify the old and adjacent chart
amplitudes, replace the completed Poisson identity by edge cancellation,
close `Q_K-T` or `T_upper`.  No height-uniform theorem, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion is proved.
