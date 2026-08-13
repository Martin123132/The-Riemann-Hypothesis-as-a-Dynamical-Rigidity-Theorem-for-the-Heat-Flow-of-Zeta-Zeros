# Endpoint-safe tangential extraction of the exterior B trace package

Date: 2026-08-13

Status: exact global decomposition and saved-height phase/cutoff certificate;
not a proof of either remaining quantitative bound or RH

Section 11.414 gives `R_end=E_Btr,win+R_group`.  A direct norm of the left
B trace tail is inadmissible because the Abel corner uses the joined endpoint
difference.  Put `delta=10^-4` and define

```text
S(u)=6u^5-15u^4+10u^3,
eta_delta(x)=0                         (x<=delta),
             S((x-delta)/delta)       (delta<x<2delta),
             1                         (x>=2delta).       (TE1)
```

The endpoint values and first two derivatives match, so `eta_delta` is C2.
Moreover

```text
|eta_x| <= 15/(8delta),
|eta_xx|<=10sqrt(3)/(3delta^2).                         (TE2)
```

Let `E_W=1-chi_W`.  At every finite symmetric cutoff define

```text
R_Btr,tan,M = eta_delta E_W Btr_M,
R_join,M    = G_M-chi_W Btr_M-eta_delta E_W Btr_M.     (TE3)
```

Then exactly

```text
G_M=chi_W Btr_M+R_Btr,tan,M+R_join,M.                  (TE4)
```

Because `eta_delta=0` for `x<=delta`, `R_join,M` retains the complete
extracted B trace package there with every omitted nonlocal bulk step, the A
endpoint, zero and negative modes, completion, and endpoint half-current at
the singular corner.  Thus the
previous `Delta_g=g_B-g_A=O(z^(1/2))` majorant and the exact double-Weber
zero/half-current collapse are not broken.

The transition ends at `x=2e-4`, while

```text
x_low=[0.0002424403205207277475284069049234717540085201398212267695529400170214101502687286288150767276271714442 +/- 5.39e-102],
xi(2delta)=[-12435.98320900291736763579846197502347052490460041914213888634868430577822211730004333373035044886384 +/- 1.85e-93]<-12000.       (TE5)
```

Outside the removed window the B trace phase is already tangentially
nonstationary.  Monotonicity gives

```text
-G_B(x_low)=[70.06936513866787354373883150762389179239019416232978219764525104266270345185533994799513262939446318 +/- 1.58e-93]>70.06,
 G_B(x_high)=[69.93077206148275093313326891269402029120313352016808666155479066569716994294537666200026909337132396 +/- 1.57e-93]>69.92,
 |G_B|>69.92 on supp(eta_delta E_W).                   (TE6)
```

In normalized trace coordinate `xi=sqrt(H_B)(x-x_B)`, the cutoff costs only

```text
|eta_xi|  <= [6.435040354771821948161022961107351145071253988958392617976371276200606448308311877126197226975257431e-5 +/- 1.44e-102] <6.5e-5,
|eta_xixi|<= [6.800477029738348063662689015111232456166379608707378349041481564656023317643188134793876792327333210e-9 +/- 3.02e-106] <7e-9. (TE7)
```

Keeping one common Abel regulator therefore leaves the exact target

```text
R_end=E_Btr,win+lim_(eps->0)lim_(M->infinity)
                  [R_Btr,tan(M,eps)+R_join(M,eps)],
limsup_(eps->0)lim_(M->infinity)[R_Btr,tan+R_join]
                  <1.4058e-4  ==> R_end<8.6e-6.        (TE8)
```

No separate `eps->0` limit for the two new summands is asserted here.

The next quantitative step is to bound the first two `xi` derivatives of
the complete summed exterior B trace amplitude before absolute values and apply
tangential integration by parts using (TE6).  The residual `R_join` must
remain in its Abel-theta/double-Weber endpoint grouping.

Pi provenance: `pi` enters only through the exact equation-(9) B trace
phase and its Hessian.  The cutoff is a rational polynomial and introduces
no fitted geometric constant.

Proof boundary: exact endpoint-safe exterior decomposition, C2 cutoff, and
saved-height tangential phase margins only.  Neither a bound for
`R_Btr,tan` or `R_join`, a complete `Q_K-T` or `T_upper` result, nor a
height-uniform theorem is obtained.  These results do not prove
`Lambda<=0`, PF-infinity, RH, or a prize-level conclusion.
