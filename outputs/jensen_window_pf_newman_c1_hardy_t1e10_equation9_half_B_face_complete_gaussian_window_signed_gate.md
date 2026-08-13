# Complete signed B trace-current package on the Gaussian window

Date: 2026-08-13

Status: rigorous saved-height `|xi|<=70` B trace-current package; not a proof
of the complete endpoint B sector or RH; the nonlocal constant bulk-step
completion and all other global sectors remain open

The nonlocal B modes are partitioned into the signed first pole-subtracted
current, absolute second/third-current allowance, post-third-current
Fresnel remainder, and the exact Fresnel-to-boundary truncation correction.
The local modes 621 and 622 are supplied by the cancellation-safe exact
step-tail theorem.  Thus the complete trace-package upper estimate is

```text
E_Btr,win^+
 = E_nonlocal,0
   + |E_nonlocal,1:2| + |E_nonlocal,rem| + |E_dictionary|
   + E_local,621:622^+ .                              (BW1)
```

The certified terms are

```text
signed nonlocal first       [4.79191597200657914096854983956099727687375298409419715416966481647269346729000000000000000e-6 +/- 2.15e-15]
nonlocal higher allowance   [1.48333110042417405479674593373890240070700572378172123814123566375516583316038283831972496e-6 +/- 1.14e-90]
nonlocal Fresnel remainder  [1.01846929988961541208209402556578689943781768458882409670361826583046648201546183369395597e-13 +/- 2.99e-103]
truncation dictionary       [2.00862353054131426303862634140009084958407713285393078563695627817674355810620416408226959e-20 +/- 1.49e-110]
exact local upper estimate  [-0.000138256127654637100770351436648671729221435261123385096643960000000000000000000000000000000 +/- 1.22e-19]
--------------------------------------------------------------------------
complete trace-package upper [-0.000131980880480359397499389294254019796601012398470818913951995361310102457746940187192698511 +/- 2.15e-15]
                          < -1.3198e-4.               (BW2)
```

There is no double counting in BW2.  The exact local theorem already charges
the full positive normal-safe local Fresnel remainder adversely, so the
same remainder from the earlier partial ledger is not added again.  The
nonlocal Fresnel remainder excludes modes 621 and 622, and the cotangent and
higher-current sums have those modes removed algebraically.

For `tau_m=1_(622<=m<=39894)`, the complete positive B endpoint sector also
contains the step coefficient

```text
sigma_m(x)=(1-tau_m)-1_(q_B,m(x)<0).                  (BW3)
```

BW1 includes `sigma_m P_bulk,m` for the two local crossing modes 621 and
622, because their jumps must remain attached to the exact local tails.  For
every nonlocal mode it contains the sign-adapted trace tail `U_B,m` but not
the constant-on-window term `sigma_m P_bulk,m`.  Those nonlocal steps remain
in the global completion.  Therefore the rigorous negative conclusion is
for the complete trace-current package on the window, not for the complete
endpoint B sector.  The exterior trace package and the omitted nonlocal
steps must be handled in the global grouped remainder.

Pi provenance: all terms retain the equation-(9) normalization and phases;
the `pi/8` projection is the exact odd-square half-Kummer reflection.  No
fitted value of `pi` or empirical phase is introduced.

Proof boundary: the complete saved-height B trace-current package restricted
to `|xi|<=70`, including the exact local crossing steps only.  No bound for
the nonlocal constant bulk-step completion, exterior trace package, complete
endpoint B sector, A-fold splice, complete `Q_K-T` or `T_upper`,
height-uniform theorem, `Lambda<=0`, PF-infinity, RH, or prize-level
conclusion is proved.
