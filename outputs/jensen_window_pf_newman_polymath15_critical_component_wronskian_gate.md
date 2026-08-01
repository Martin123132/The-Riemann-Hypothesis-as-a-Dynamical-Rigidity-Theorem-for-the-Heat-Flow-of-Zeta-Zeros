# Jensen-Window PF Newman Polymath-15 Critical Component Wronskian Gate

Date: 2026-07-24

Status: exact corrected-component Wronskian reduction, exact route
guard, and finite cancellation diagnostics. The arithmetic small-ball
theorem remains open; this is not a proof of `Lambda <= 0` or RH.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_critical_component_wronskian_gate.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_critical_component_wronskian_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_critical_component_wronskian_gate.py
```

## Pairwise Identity

```text
Write E=sum_(j=0)^N e_j, where e_n=exp(i*beta_t)*exp((t/4)log(n)^2-s_*(x)log(n)) for 1<=n<=N and e_0=-q_(N,t) is the sharp endpoint pseudo-component
For e_j=a_j exp(i*theta_j), write e_j'=(u_j+i*v_j)e_j, u_j=(log a_j)' and v_j=theta_j'
W_E=Im(E'*conj(E))=sum_j a_j^2 v_j+sum_(j<k)a_j*a_k*((u_j-u_k)sin(theta_j-theta_k)+(v_j+v_k)cos(theta_j-theta_k))
```

Equivalently,

```text
For z=(e_0,...,e_N)^T and r_j=u_j+i*v_j, W_E=z^*Kz with K_(jk)=(r_k-conj(r_j))/(2i); K is Hermitian and rank(K)<=2
If m=N+1, the product of the two possible nonzero eigenvalues is lambda_+*lambda_-=-(m*sum_j|r_j|^2-|sum_j r_j|^2)/4=-(m/4)*sum_j|r_j-r_avg|^2; unequal component rates give inertia (1,1,m-2)
```

Thus the unrestricted instantaneous rate form has one positive and
one negative direction whenever the component rates are not all
equal. A positive-semidefinite Gram argument cannot arise from this
matrix alone; arithmetic phase placement or additional Xi structure
must remove its negative direction on the relevant crossing set.

For the actual Dirichlet components,

```text
For 1<=n<=N, u_n=-Re(s_*')log(n) and v_n=beta_t'-Im(s_*')log(n); at t=0, u_n=0 and v_n=beta_0'+log(n)/2
```

At `t=0` the Dirichlet amplitudes have zero radial rate and their phase
speeds increase toward the saddle endpoint. The endpoint correction is
one additional explicit component, so no cross term is discarded.

## Reference Rate

```text
For any real u_* and v_*, set R_*=E'-(u_*+i*v_*)E. Then W_E=v_*|E|^2+Im(R_*conj(E))
At Re(E)=0 with E=iY, Re(E')=Re(R_*)-v_*Y; a nonzero double crossing is equivalent to Re(R_*)=v_*Y
```

This makes endpoint-relative formulations exact, but it does not make
the defect term one-signed.

## Double-Crossing System

```text
Re(E)=sum_j a_j cos(theta_j)=0 and Re(E')=sum_j a_j*(u_j cos(theta_j)-v_j sin(theta_j))=0
```

The two equations require joint cancellation of different trigonometric
observables. Upper bounds for either sum separately do not exclude their
simultaneous smallness.

## Ordered-Speed Countermodel

```text
E_toy(x)=exp(-3ix)-exp(-4ix)+i*exp(-ix)-(i/2)exp(-2ix) has positive component amplitudes and distinct negative phase speeds {-4,-3,-2,-1}, but E_toy(0)=i/2, E_toy'(0)=i, Re(E_toy) has a genuine double zero with second derivative 7, and W_Etoy(0)=0
```

Thus positive component amplitudes, distinct ordered phase speeds, and
a common speed sign are not a transversality theorem.

## Corrected Diagnostics

| label | N | crossing | W_E | phase speed | diagonal | pair radial | pair phase | cancellation |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| ordinary_negative | 6 | 517.2201245916915343451981788255318076677 | -0.80300991218331942075386569183402424 | -0.55142810044966476597942521956113449 | -1.4370300717924581374224329704221245 | -0.00025315178321500870603088573947228110 | 0.63427331139235372537459816432757259 | 5.253481734630345810501480 |
| ordinary_positive | 6 | 519.7486269586060360224240422895051497992 | 0.14735289688380540247293482207795179 | 1.2989222762566859009109017488767268 | -1.4400108434408802277152839682295695 | -0.000059120302277226995896454355426738281 | 1.5874228606269628571841152446629481 | 32.18263398691910292434156 |
| lehmer_left | 33 | 14010.12594382089318324866358493301988117 | 0.059350057399241023592912928392882718 | 0.17890906267654252236906178230689983 | -4.1258998258446731346388431021753458 | -0.000011369180498191117192875631543452886 | 4.1852612524244123493489489061997719 | 771.5103365551727133404041 |
| lehmer_right | 33 | 14010.20090109869997081492153406456328151 | -0.072035416168465488766306515643226443 | -0.15635696774718775311741288197780475 | -4.1259052934109920297648555297577638 | -0.000013386730335568636293413333810094745 | 4.0538832639728621096348424274483475 | 634.1290671153968769113621 |

All component phase speeds in these four rows are negative, the
Dirichlet speeds are strictly increasing, and the endpoint speed is
the largest. Nevertheless the aggregate Wronskian has both signs.
The largest absolute-term cancellation factor is
`771.5103365551727133404041`.
These are finite `t=0` diagnostics, not an asymptotic theorem.

## Live Target

```text
Prove X^2+(U/L)^2>8000000*exp(-3L/2), X=sum_j a_j cos(theta_j), U=sum_j a_j*(u_j cos(theta_j)-v_j sin(theta_j)), on L>=50 and 0<tL<=c_*+o(1)
```

The missing input is a joint arithmetic small-ball or phase-critical-value exclusion for the logarithmic Riemann-Siegel phases; componentwise phase ordering, absolute upper bounds, and a positive-semidefinite argument based only on the instantaneous component-rate matrix cannot supply it.
