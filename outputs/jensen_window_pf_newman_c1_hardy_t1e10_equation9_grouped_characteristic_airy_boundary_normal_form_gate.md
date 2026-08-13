# Grouped characteristic Airy-boundary normal form

Date: 2026-08-10

Status: exact grouped two-variable normal form and compact-core error budget
validated; not a proof of a complete fold-uniform estimate

The previous lower-endpoint calculation found the correct cubic boundary
phase, but it did not justify treating the remaining grouped current as a
slow scalar amplitude.  Keep the endpoint and Fresnel current together and
sum the finite transition roster before taking absolute values.

Put

```text
eta=pi*C^2/(8t),          beta=(t*eta)^(1/3),
u=2z/beta,                x(z)=[1-tanh(z/beta)]/2,
alpha=C+sigma*y,          sigma=4beta/(pi*C),
d_m=4beta(m-C/4)/C.                                      (GB1)
```

After removing the common constant `pi*C^2/8`, oddness of `C` gives the exact
joint phase

```text
Theta_m(z,y)
 =t[z/beta-eta*tanh(z/beta)]
  -d_m*y-beta*y*tanh(z/beta)+x(z)y^2/(2beta).           (GB2)
```

The transformed Kummer measure, alpha current, and both Jacobians combine to

```text
A(z,y)=sqrt(2)/pi cosh(z/beta)^(-3/2)(1+sigma*y/C).     (GB3)
```

Thus there is no separated endpoint series.  For the 84 modes
`39853..39936`, all mode dependence is the exact finite kernel

```text
D_84(y)=sum_m exp(-i*d_m*y)
 =exp[-i(d_0+83h/2)y] sin(84hy/2)/sin(hy/2),
h=4beta/C.                                              (GB4)
```

At removable zeros of the denominator, (GB4) is interpreted by continuity.
This is the cancellation object that the separated endpoint/Fresnel formula
conceals.

The correct local canonical phase is two-dimensional:

```text
Theta_m^0(z,y)
 =z^3/3-lambda*z-(z+d_m)y+y^2/(4beta),
lambda=(eta-1)t/beta.                                   (GB5)
```

The exact remainder is

```text
R=R_Airy(z)+y[z-beta*tanh(z/beta)]
  +[x(z)-1/2]y^2/(2beta).                               (GB6)
```

Using `|a-tanh(a)|<=|a|^3/3`, `|tanh(a)|<=|a|`, and the
certified Airy remainder, the rectangle `|z|<=4`,
`0<=y<=64` satisfies

```text
|R_Airy| <= [0.000941286188129008107952594802791584350179654934528794474498182145317101861576845180190284947088805975608031353984016027642 +/- 7.06e-123],
|R-R_Airy| <=
  [0.000294151933790315033735185875872370109431142167040248273280681920411594331742764118809464045965251867377509798120005008638 +/- 2.11e-123]
 +[0.000882455801370945101205557627617110328293426501120744819842045761234782995228292356428392137895755602132529394360015025914 +/- 6.12e-123],
|R| <= [0.00211789392329026824289333830628106478790422360268978756762090982696347918854790165542814113094981344511807054646403606219 +/- 1.87e-122] < 0.0022. (GB7)
```

On the same rectangle,

```text
|A/(sqrt(2)/pi)-1|
 <= [6.89418594821050860316841896575867443979239454000581890501598250964674215022103403459681357731059064166038589343761738991e-6 +/- 2.83e-125] < 1e-5. (GB8)
```

The numerical characteristic data are

```text
beta=[2154.43548064036880514538712430952510213643988046170017933872008509078888609877630125920487705916580649509608031828427003 +/- 9.09e-117],
sigma=[0.0171898986102748176776219811452948748606054834923516963032146162647171589391534679396696206285389394191287406205796042223 +/- 9.59e-122],
lambda=[5.10994303949183562830591916748466475096631702036544079267907046952026199304814055685009392261585095593930371371474634517 +/- 1.96e-114],
h=[0.0540036591899927634971302161165963792310029610899239910347661651764549749926061099346197729512189302091177570782326844101 +/- 1.91e-121],
[-2.22765094158720149425662141480960064327887214495936463018410431352876771844500203480306563423778087112610747947709823192 +/- 1.14e-119] <= d_m <=
[2.25465277118219787600518652286789883289437362550432662570148739611699520594130508977037552071339033623066635801621457412 +/- 9.93e-120].
```

The normal core reaches
`[1.37883794862575998593541601365769216701186007544100336514200678866390052833442174326943248068024829891455939980906622941 +/- 2.69e-120]` times the natural Fresnel scale
`sqrt(beta)`.  This is enough to validate the local model but not enough to
discard the complement: for some `(z,m)` the normal saddle lies beyond this
rectangle.  The next obligation is therefore to bound the finite-kernel
canonical integral and make an endpoint/interior partition whose outer piece
contains every displaced normal saddle.  A one-dimensional Airy estimate
with a frozen or slowly varying amplitude is not admissible.

Pi provenance: the `pi` in `sigma`, `h`, and the phase is inherited from the
original Kummer factor `exp(i*pi*x*alpha^2/4)` and the integer
Fourier--Poisson character `exp(-i*pi*m*alpha)`.  The factor `sqrt(2)/pi` in
(GB3) is forced by those coordinates and their two Jacobians.  No circle,
polygon, or fitted geometric normalization is introduced here.

Proof boundary: exact coordinate algebra and a compact numerical error budget
at `t=10^10`.  No canonical-integral bound, outer normal tail, complete
`T_upper` assembly, height-uniform source error, `Lambda<=0`, RH, or
prize-level conclusion is proved.
