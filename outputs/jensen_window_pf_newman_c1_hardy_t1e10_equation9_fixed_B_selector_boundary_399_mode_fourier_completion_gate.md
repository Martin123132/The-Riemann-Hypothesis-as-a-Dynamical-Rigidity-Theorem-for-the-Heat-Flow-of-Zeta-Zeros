# Fixed-B selector-boundary 399-mode Fourier completion

Date: 2026-08-12

Status: exact canonical one-period completion and rigorous 399-mode interval
quadrature validated; this is not a proof of finite-t Airy-model transport,
`T_upper`, `Lambda<=0`, or RH

Hold the theorem endpoint fixed at `B=5122423` and consider the selector
strip

```text
A=159577 <= alpha <= A+2=159579,
t*=pi*A^2/8=[10000011009.04258833249996651073906089290177178617834443516760678611578 +/- 4.07e-60].       (FC1)
```

Use the exact lower-fold scales

```text
beta=t*^(1/3), sigma=4beta/(pi A), h=4beta/A,
alpha=A+sigma*y, 0<=y<=Y=2/sigma.                      (FC2)
```

They satisfy the exact one-period identity

```text
hY=2pi.                                                 (FC3)
```

The joint saddle is `alpha_m=2m+t*/(pi m)`.  Its two roots at
`alpha=A+2` lie in

```text
[39695.01392413987902721080327585248206056114167738689526948765453852349 +/- 5.04e-66],
[40094.48607586012097278919672414751793943885832261310473051234546147651 +/- 4.44e-66].           (FC4)
```

Consequently the strip contains exactly the 399 integer saddle modes
`39696..40094`.  They are the symmetric block

```text
m=39895+k, -199<=k<=199,
d_m=h(k+3/4).                                           (FC5)
```

The closest lower-fold saddle clearance is only
`[3.133303253622098561187145936732340702862585852509149245500574533310525e-6 +/- 3.57e-66]` and the first excluded lower
outer mode misses the upper endpoint by only
`[0.0002802619977327119284544652978964605113994205819372716966872402062649357 +/- 1.19e-66]`.  Its reciprocal canonical
branch margin exceeds `[1272.667086017669035518721940110971297465407710496310742906249304526521 +/- 9.07e-63]`.
Thus a per-mode nonstationary triangle bound is structurally unsuitable.

For each mode define the canonical coefficient

```text
G_m=int_0^Y Ai(-y) exp(i*y^2/(4beta)-i*d_m*y)dy.        (FC6)
```

All 399 coefficients were integrated as independent complex balls in an
append-only cache.  Their grouped value is

```text
2pi sum_(m=39696)^40094 G_m
 =[69.03605876493528260441837879984535362966373755138805650491325413392854 +/- 4.74e-17]
  +i*[-1.704180601685076962575640131572523825798168416766669990981302083618126 +/- 4.74e-17]. (FC7)
```

A second interval calculation using the exact 399-term Dirichlet kernel gives

```text
[69.03605876493528260441837879984535363053241291972437154856666911417188 +/- 1.54e-18]
 +i*[-1.704180601685076962575640131572523825502628752622056715580552249306899 +/- 1.54e-18],       (FC8)
```

and both components overlap.

The key completion is exact.  Put `s=y/Y`.  By (FC3), `G_(39895+k)` is the
`k`th Fourier coefficient of

```text
g(s)=Y Ai(-Ys) exp(i*Y^2*s^2/(4beta)-3pi*i*s/2).       (FC9)
```

The endpoint phase simplifies exactly to `-pi`, so

```text
g(0)=Y Ai(0), g(1)=-Y Ai(-Y).                          (FC10)
```

Dirichlet-Jordan therefore gives the complete symmetric mode sum as
`[g(0)+g(1)]/2`.  The zero, negative, and outer-positive modes, retained as
one cancellation object, are exactly

```text
midpoint-sum_399
 =[-0.1267849362235905264224533026357960975024975008142896709546302321471283 +/- 7.53e-18]
  +i*[0.2712287666795003765320076447956663081493369635859892929772493144837152 +/- 7.53e-18]. (FC11)
```

After the common `2pi` factor this complement is

```text
[-0.7966132484517649224255456604087675297667206530978270539440382689250442 +/- 4.74e-17]
 +i*[1.704180601685076962575640131572523825798168416766669990981302083618126 +/- 4.74e-17].          (FC12)
```

Adding the endpoint half-current, the 399 block, and (FC11) reconstructs
`g(0)` exactly.  No outer mode is bounded separately and no divergent
endpoint current is detached.

This closes the complete canonical Fourier bookkeeping at the selector
boundary.  It does not yet prove that the exact finite-`t` Kummer strip is
within a useful explicit error of this canonical completion.  That
cancellation-preserving finite-t comparison is the next theorem target.

Machine-audited companion:

```text
outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_fixed_B_selector_boundary_399_mode_fourier_completion_gate.md
work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_fixed_B_selector_boundary_399_mode_fourier_completion_gate.json
work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_fixed_B_selector_boundary_399_mode_fourier_completion_cache.jsonl
work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_t1e10_equation9_fixed_B_selector_boundary_399_mode_fourier_completion_gate.py
work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_fixed_B_selector_boundary_399_mode_fourier_completion_gate.py
```

No complete finite-t strip error, `T_upper`, `Lambda<=0`, RH, or prize-level
conclusion is proved.
