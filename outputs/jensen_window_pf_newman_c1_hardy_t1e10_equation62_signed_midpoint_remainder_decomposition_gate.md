# Signed diagnostic-midpoint decomposition gate

Date: 2026-08-09

Status: finite certified signed decomposition; not a proof of RH

Let `E_src=Q-T_src` be the final equation-(62) residual at the published/source
endpoint, let `T_5=T_mid-T_src` be the exact-theta five-term change to the
introduced diagnostic midpoint, and let `Delta_D` be the certified real
five-saddle B9 transport.  Pure algebra gives

```text
E_mid = E_src-T_5
      = (E_src-Re Delta_D) + (Re Delta_D-T_5).
```

At the central output the two channels and their sum are

```text
E_src-Re Delta_D = [-0.0041182964137970075590400595403186846184790696475438888485910944149858295893902652342229225089884265417546573920 +/- 1.72e-38]
Re Delta_D-T_5   = [0.0055864309068055834077066211896712245656681187696961472778526477290758418536031142659431632725265327575835469950 +/- 1.72e-38]
signed sum        = [0.0014681344930085758486665616493525399471890491221522584292615533140900122642128490317202407635381062158288896030 +/- 3.43e-38]
```

Across all fifteen outputs, the complementary channel is strictly negative in

```text
[-0.0041582987250336848991932763648390901117566832471227153780232306756884478522726059059356545316498168673289933380 +/- 1.79e-38]
[-0.0040782886556580794191937957037847142505253460585377780223481026561352507596445961821299032886613631384361119790 +/- 1.68e-38]
```

while the D-to-boundary channel is strictly positive in

```text
[0.0055479018837298056772630043531304614759725067248252900381656062608263595569703922702772925070750662657478899260 +/- 1.68e-38]
[0.0056248910154118848705125318211745565809326910178214809938758816398844926270310810441530460449477813769642977000 +/- 1.79e-38]
```

Their signed sum reproduces the independently computed diagnostic midpoint residual in

```text
[0.0014665922903781999713192554563354664691760077706987656158526509641960447747584751382173915132979645096353043620 +/- 3.56e-38]
[0.0014696132280717262580692086493457472254471606662875120158175036046911087973257960881473892184137031273117779470 +/- 3.34e-38]
```

and is below `0.005` at every saved height.  In contrast, bounding the two
channels independently gives an absolute triangle between
`[0.0096261905393878850964568000569151757264978527833630680605137089169616103166149884524071957957364294041840019050 +/- 3.34e-38]` and
`[0.0097831897404455697697058081860136466926893742649441963718991123155729404793036869500887005765975982442932910380 +/- 3.56e-38]`, above `0.005` everywhere.
Only about 15 percent of that unsigned size survives the opposite-sign
cancellation.

This finite gate identifies correlation in the diagnostic alternative.
It does not prove that either sign persists for every admissible radius or
height.  More importantly, the `D` split is auxiliary: while the source and
diagnostic-midpoint residue rosters stay fixed, changing a contour selector gives
`dA=-d(Re Delta_D)`, `dB=d(Re Delta_D)`, and hence `d(A+B)=0`.  This does not
make `E_mid` the paper's intrinsic error: the midpoint cutoff itself still
needs an independent alternative-hybrid theorem.  The primary next theorem is
a direct signed bound for the source-aligned equations-(126)--(127) error.
Selector-uniform channel signs matter only if a separate midpoint theorem uses
this contour split.  There is no implication for `Lambda<=0`, PF-infinity, RH,
or the Clay prize.
