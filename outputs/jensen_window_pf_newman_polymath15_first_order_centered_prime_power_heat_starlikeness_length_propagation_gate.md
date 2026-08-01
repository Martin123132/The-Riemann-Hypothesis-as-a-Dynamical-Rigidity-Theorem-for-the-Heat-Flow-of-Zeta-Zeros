# Prime-Power Heat-Starlikeness Length Propagation

Date: 2026-07-30

Status: one exact dyadic next-length certificate and a
uniform all-length theorem for every complete p>=5 chain.
This does not promote RH or Lambda<=0.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_first_order_centered_prime_power_heat_starlikeness_length_propagation_gate.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_prime_power_heat_starlikeness_length_propagation_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_prime_power_heat_starlikeness_length_propagation_gate.py
```

Current result:

```text
validated prime-power heat-starlikeness length propagation gate: 20 rows, 1 append-one identity, 1 exact dyadic next-length certificate, 1 analytic all-length p>=5 family, 2 actual ray-positive propagated families, 2 open small-prime length families, 0 joined Abel gaps, 0 successor winding bounds
```

## Append-One Identity

```text
c_(k,M+1)=q^(2k)c_(k,M), 0<=k<M
d_M=c_(M,M+1)=p^(-M/2)s^M q^(M^2)
Q_(M+1)(z)=Q_M(q^2 z)+d_M z^M
J_(M+1)=J_M(q^2 z)+(M+1)d_M^2+d_M sum_(k=0)^(M-1)(M+k+2)c_(k,M)q^(2k)cos((M-k)theta)
```

The old block is evaluated at `q^2 z`; the final displayed
sum is the entire new cross current. No absolute geometric
`C1` approximation is inserted.

## Dyadic M=9

The exact rational tensor-Bernstein certificate gives

```text
p=2, M=9, q>=47/50, s>=22/25
power degrees=(128, 16, 8)
Bernstein shape=(129, 17, 9)
angular leaves=32
global coefficient lower diagnostic=5.380434103075095297082257837106824887811e-2
J_0>1/20
```

Thus the first unresolved dyadic extension is positive; this
is not yet an induction for `M>=10`.

## All Large-Prime Lengths

Define the Fourier coefficients by

```text
A_0=2sum_(k=0)^(M-1)(k+1)c_k^2; A_d=sum_(k=0)^(M-1-d)(2k+d+2)c_kc_(k+d)
r_k=(s/sqrt(p))q^(2M-3-2k); 0<r_k<=r_(k+1)<=1/sqrt(5)<1/2
```

The ratio sequence is increasing and bounded by
`1/sqrt(5)<1/2`. Direct grouping proves that `A_d` is
decreasing and convex after adjoining `A_M=A_(M+1)=0`.
The only terminal expression is

```text
G_(A,d)=A-2(A+1)aq+(A+2)a^2q^(2d+2)
G_(A,d)>=d+2-2(d+3)q/sqrt(5)+(d+4)q^(2d+2)/5.
```

For `0<=d<=6`, all ordinary degree-`2d+2` Bernstein
coefficients are positive using `2/sqrt(5)<9/10`. For
`d>=7`, dropping the last positive term leaves
`(d-7)/10>=0`.

Fejer's exact identity is

```text
F_d(theta)=1+2sum_(m=1)^d(1-m/(d+1))cos(mtheta)>=0
J_0=(1/2)sum_(d=0)^(M-1)(d+1)(A_d-2A_(d+1)+A_(d+2))F_d
Delta^2 A_0>=14/5-6/sqrt(5)>1/10
p>=5,M>=2: J_0>1/20
```

The circle constant does not enter this argument. The
`cos(d theta)` terms have the ordinary `2*pi` angular period.

## Actual Fixed-Ray Transfer

Because `k log(p)<=log(N)`, the relative coefficient estimate
is uniform in chain length. Geometric sums replace the old
`M<=8` count:

```text
|R_k-1|<r_x=18000/x
sum c_k<20/11; sum k c_k<180/121
|J_ang-J_0|<21r_x=378000/x
|J_ray-J_0|<378006/x+405408/x^2<1/1000
|J_ray-J_0|<3456096/x+1621632/x^2<1/1000
p=2,M=9 and p>=5,M>=2: J_ray>49/1000
```

## Boundary

This gate proves the exact append-one identity, J_0>1/20 and J_ray>49/1000 for the dyadic M=9 family, and J_0>1/20 and J_ray>49/1000 for every complete p>=5 chain of every length M>=2. It does not prove p=2,M>=10, p=3,M>=5, short or singleton chains, p-free or endpoint joins, a Xi Abel gap, a successor winding bound, Lambda<=0, PF-infinity, or RH.
