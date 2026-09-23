# The arithmetic energy uncertainty is removed

The source's quadratic energy is now determined, rather than bounded through
an unknown low-frequency spectral measure. The actual common-height ratio
phases supply the missing cancellation. This is an analytic result for the
two original window families, not an extrapolation from the N=1024 example.

The classical exponent pair13/31,16/31 leaves a strict -1/93 power before
epsilon losses. With fixed smooth endpoint cutoffs, the complete ratio
off-diagonal is O(N^(-1/200)). Restoring the exact windows gives the needed
o(log N) result, without claiming that the power rate survives restoration.

The previously unknown Gram measure is therefore ordinary uniform measure
on[0,1]. Both derivative covariance limits follow. This does NOT identify
the nonlinear distribution or justify a Gaussian/random-phase model.

## Numerical consequence

The SAME frozen majorant profile and SAME source normalization give

    0.1512131 < limiting complex energy < 0.1512159.

The former worst-case allowance was approximately0.1536482. The difference
is modest, but the entire unknown energy parameter is now gone. In the
direct-hinge formulation the available signed room grows from roughly
0.002368 to more than0.002976 per log N, an increase of about26%.

The remaining sufficient source targets can be stated without maximizing
over an unknown spectral measure:

* One-mode certificate: limsup M/log N < -0.02125 suffices with b=9/16.
* Richer direct hinge: limsup (J-Q/2)/log N <0.002975 suffices.
* With its existing paid M=16 Fourier truncation, the retained signed
  target can be0.002934, up from the earlier conservative0.00232.

NEITHER signed bound is established. They are alternative sufficient routes,
not measured values or necessary conditions for the actual loss. Keep every
original physical harmonic, weight normalizer and real-source endpoint cost.

## Verification

Three rigorous whole-cell Arb integrations agree:96bits/4096cells,
160bits/32768cells and192bits/262144cells, with a fully paid infinite tail.
Eleven test families pass, including256 exact phase-derivative controls,
26 independent rational profile coefficients and ten elementary-hyperbolic
checks. Four invalid rate certificates and eight altered-result cases are
rejected. A separately coded60/90digit mpmath quadrature agrees, but is
reported only as a diagnostic, not an independent rigorous integration.

The high-resolution computation has128 append-only checkpoints. Clean
extraction and historical integrity receipts accompany the sealed delivery.
No new source-height evaluation, model execution, zero count or repository
publication was needed.

## Next mathematical target

Control the actual signed nonlinear response using the now FIXED quadratic
kernel. The one-mode route needs |C|C with all physical modes0 through3;
the direct-hinge route retains more phase information and has the larger
newly determined signed budget. Do not spend another stage merely optimizing
the old unknown-energy scalar bound: that dependency has been resolved.

This is a reduction of the remaining proof problem, not a proof of RH,
all-height positive hinge control or absolute zero-count closure.
