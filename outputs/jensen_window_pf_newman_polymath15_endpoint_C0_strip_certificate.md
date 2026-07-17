# Jensen-Window PF Newman Polymath-15 Endpoint C0 Strip Certificate

Date: 2026-07-17

Status: rigorous complex-strip bounds for the endpoint factor. This
is not a proof of `Lambda <= 0` or RH; corrected transversality remains open.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_endpoint_C0_strip_certificate.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_endpoint_C0_strip_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_endpoint_C0_strip_certificate.py
```

## Domain

```text
C0(p)=(exp(pi*i*(p^2/2+3/8))-i*sqrt(2)*cos(pi*p/2))/(2*cos(pi*p))
|Re(p)|<=101/100, |Im(p)|<=1/100
```

The apparent poles at `p=+-1/2` are removable. They are handled by
Cauchy's estimate rather than unstable quotient evaluation.

## Removable Disks

On each radius-`1/4` outer circle,

```text
For z=p0+w, |w|=1/4 and p0=+-1/2, |2*cos(pi*z)|=2*|sin(pi*w)|>=sqrt(2)
sin(pi*|Re w|)>=2*sqrt(2)*|Re w| and sinh(pi*|Im w|)>=2*sqrt(2)*|Im w|
numerator bound = [3.326939835774423430690659078432820102823375035559527541860310199573320 +/- 7.37e-71] < 4
|C0| bound      = [2.352501718475753584400918937474195990201736726366240237577739760870626 +/- 2.24e-71] < 3
```

The gap from outer radius `1/4` to inner radius `3/20` is `1/10`,
so Cauchy gives

```text
|C0'| <= [23.52501718475753584400918937474195990201736726366240237577739760870626 +/- 2.24e-70] < 30
```

## Outside Region

Outside the two inner disks,

```text
Outside the two radius-3/20 disks, the strip width 1/100 forces dist(Re p, Z+1/2)>1/10 and |2*cos(pi*p)|>=2*sin(pi/10)
direct |C0| bound  = [3.958725464747200207594448495304580276282063847079090508337807373303300 +/- 1.73e-70] < 5
direct |C0'| bound = [49.16050018846990534165820020182535524760359511470658899392748619238243 +/- 2.82e-69] < 60
```

Combining the covers proves

```text
On |Re p|<=101/100 and |Im p|<=1/100, |C0(p)|<5 and |C0'(p)|<60
```

## Endpoint Collar

For the critical scale,

```text
For L>=50 and |z-x|<=1/L, the fixed-N map p_N(z)=1-2*(sqrt((z/2+pi*t/8)/(2*pi))-N) stays in this strip
|Delta p| <= [2.777588772992804118932352749217371382079952076041011116709559992176274e-14 +/- 3.73e-84] < 1/100
```

This closes the finite-dimensional derivative input in the endpoint
holomorphic-lift reduction. The next step is to compose it with the
published `e_A+e_B+e_C` bounds and produce explicit normalized `C1`
remainder constants. The corrected first-jet lower bound remains open.
