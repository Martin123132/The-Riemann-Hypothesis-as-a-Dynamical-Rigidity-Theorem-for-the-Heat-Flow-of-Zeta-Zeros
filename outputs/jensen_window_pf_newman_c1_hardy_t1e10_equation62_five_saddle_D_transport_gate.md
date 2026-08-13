# Five-saddle `D` transport gate

Date: 2026-08-09

Status: finite certified transport measurement; not a proof of RH

The diagnostic pole-free deformation of Section 11.302 moves the B9 lower
limit across exactly `N=37941,...,37945`.  B9 and B15a are sourced to Lewis
(2015).  Direct interval evaluation of the paper's phase
is unstable because two order-`t` terms cancel at each saddle.  Put
`w=sqrt(pc)`, `w_s=a/(4N)`, and `q=w/w_s-1`.  Exact algebra gives

```text
Phi_N(w)-Phi_s
 =t[log(1+q)-q+q^2/(2w_s^2(1+q)^2)].
```

The gate evaluates `log(1+q)-q` by a 40-term analytic series with an explicit
geometric tail on `|q|<0.5`; both real contour endpoints remain inside
`|q|<0.01`.  Independent 70/110/130-digit Arb integration encloses all 75
finite B9 integrals.

At the central output, after applying the exact equation-(13)/(24) prefactor,

```text
Re(transformed five-D collar) = [-0.012021611549728283097321170603987145201879746977306324299057265825477626294385557540534525074687277030995238718 +/- 1.71e-38]
five paper-phase classical terms = [-0.017608042456523167258309400517110434327569506595113977969726826311395121505489265940710310102448076451263440201 +/- 1.63e-101]
D transport minus classical      = [0.0055864309067948841609882299131232891256897596178076536706695604859174952111037084001757850277607994202682014825 +/- 1.71e-38]
```

Thus the exact five-saddle `D` interval supplies about
`[0.68273413012328887096535724473568905075878335602634229512481912424791048560273357449056281664855560546907066455 +/- 9.71e-37]` of the signed classical
boundary correction at the centre, not all of it.  Across all fifteen outputs
the unresolved `D`-versus-classical gap has absolute range

```text
[0.0055479018837189777561396392280461317307061109182269969861789177639912831360766604444424622209766299476261075361 +/- 1.67e-38]
[0.0056248910154013144295084043386481943585794572426663041425318471712176064708979218669032764052256963493152296465 +/- 1.78e-38]
```

and remains strictly above `0.005` everywhere.  The exact Riemann--Siegel theta
versus paper leading phase changes the five-term collar by less than `2e-11`,
so that phase refinement cannot close the gap.

This does not contradict Lemma B1.2: the gate integrates a difference between
two radii through a moving saddle, not a full `D(N,R)` in one fixed lemma
regime.  Its result shows quantitatively why the B18 main term cannot simply be
switched on without transporting the remaining `C`, `D`, first-integral, and
endpoint channels.  Exact contour invariance fixes those channels only in
combination as the negative `D` transport; separate triangle bounds would lose
the required correlation.  No signed bound for the diagnostic midpoint
remainder, no theorem validating that alternative cutoff, no source-aligned
height-uniform hybrid bound, and no RH implication is proved.
