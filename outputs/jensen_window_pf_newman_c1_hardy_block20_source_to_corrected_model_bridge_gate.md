# Hardy block-20 source-to-corrected-model bridge

Date: 2026-08-09
Status: finite_source_to_corrected_model_bridge_closes_at_374_saved_points_and_15_outputs; finite correction theorem only, not a proof of RH

## Signed bridge

The admitted source recurrence and the corrected mathematical model are joined
without adding unrelated absolute maxima.  With `Q_src_drop` denoting the
source nonsaddle correction after deleting the source-only `t2` term,

```text
D_src = P_src-M*C-q_src,
D_1   = D_src-(I69_src-M*C-t5_src),
D_2   = D_1+q_extra = P_src-I69_src-Q_src_drop,
D_corr= D_2+(P_pi-P_src)-(I69_pi-I69_src)
              -(Q_paper-Q_src_drop)
      = P_pi-I69_pi-Q_paper
      = Q_exact-Q_paper.
```

All 374 exact-W1 to `t2`, `t2` to W2--W5, and final endpoint-decomposition
handoffs overlap.  The complete nonsaddle transport also splits exactly into
the independently measured PSI/ERF replacement and a post-special-function
normalization remainder.

```text
maximum exact-W1 replacement             1.25603104855648916744327506343761262107625955214361E-1
maximum t2 deletion correction           7.95774715459476678844418816862571878686315931944912E-2
maximum PSI/ERF q replacement             7.15521887106601530132103123337363958369685514054674E-4
maximum complete nonsaddle transport      7.15521887106601530132103123361733529685432509201032E-4
maximum post-special normalization        1.58106230787332192827157438418188852370088220793027E-9
maximum parent source-to-2*pi transport   2.09595835444283276785241019068565666357858848460625E-31
maximum I69 source-to-2*pi transport      1.17991389847020108804403178845778952105290226953230E-31
maximum source reconstruction roundoff    5.26109688231931655659882087643541965280697931813432E-31
maximum corrected output majorant         2.19652337832631474465367310836148849889919479338472E-3
outputs below 0.005                       15 / 15
```

## Boundary

This proves the signed correction map at the saved finite inputs and connects
its final residual to the independently certified endpoint correction.  It
does not alter the accepted source, enclose the correction map over entire
physical cells, treat other blocks or height, control the outer Hardy
representation, or prove `Lambda<=0`, PF-infinity, RH, or a prize theorem.
