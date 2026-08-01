# Prime-Power Heat-Starlikeness Base Certificate

Date: 2026-07-30

Status: exact correction-free base-family positivity and a
rigorous transfer to the actual fixed-ray current. This certificate
does not promote RH or Lambda<=0.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_first_order_centered_prime_power_heat_starlikeness_base_certificate.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_prime_power_heat_starlikeness_base_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_prime_power_heat_starlikeness_base_certificate.py
```

Current result:

```text
validated prime-power heat-starlikeness base certificate: 18 rows, 3 exact Bernstein base certificates, 1 analytic p>=5 two-level family, 3 actual ray-positive base families, 0 length-propagation theorems, 0 joined Abel gaps, 0 successor winding bounds
```

## Exact Current

After separating the saddle offset and normalized correction,

```text
c_k=p^(-k/2)*q^[k*(2M-2-k)]*s^k
q=exp(-t*(log p)^2/4); s=exp(-t*(log p)*u/2); 0<=u<log(p)+1/N
Q(z)=sum_(k=0)^(M-1)c_k*z^k
H(z)=sum_(k=0)^(M-1)k*c_k*z^k
J_0=|Q|^2+Re(H*conj(Q))
J_0=sum_k(k+1)c_k^2+sum_(k<l)(k+l+2)c_k*c_l*T_(l-k)(x), x=cos(theta)
```

This identity comes directly from pairing the diagonal and
off-diagonal terms of `|Q|^2+Re(H*conj(Q))`. It replaces
the false absolute geometric-C1 target by a one-sided current.

## Exact Bernstein Certificates

On [0,1]^d, every tensor Bernstein basis function is nonnegative and the basis sums to one; therefore P>=min_I b_I.

Exact Fraction affine transforms, exact Fraction power-to-Bernstein conversion, exact de Casteljau x-splits, and 90-digit outward rational sqrt bounds.

```text
p=2, M=8, q>=47/50, s>=22/25
power degrees=(98, 14, 7), Bernstein shape=(99, 15, 8)
x leaves=('LLL', 'LLR', 'LR', 'RL', 'RRL', 'RRR')
leaf coefficient lower diagnostics=('1.581948618686023320113713421916457116787e-2', '3.899682592748647646504777592590999497967e-3', '6.628727486852878028346800733083400045840e-2', '7.8125e-2', '5.928265078834418285439714128835994766925e-1', '1.249142811973160668888802078601622451776e+0')
global coefficient lower diagnostic=3.899682592748647646504777592590999497967e-3
proved rational current lower=1/300
```

```text
p=3, M=4, q>=17/20, s>=73/100
power degrees=(18, 6, 3), Bernstein shape=(19, 7, 4)
x leaves=('L', 'R')
leaf coefficient lower diagnostics=('4.254623391080937830531237568978401909853e-2', '1.481481481481481481481481481481481481481e-1')
global coefficient lower diagnostic=4.254623391080937830531237568978401909853e-2
proved rational current lower=1/25
```

```text
p=5, M=2, q>=18/25, s>=13/25
power degrees=(2, 2, 1), Bernstein shape=(3, 3, 2)
x leaves=('',)
leaf coefficient lower diagnostics=('5.835921350012618215449579876123425873562e-2',)
global coefficient lower diagnostic=5.835921350012618215449579876123425873562e-2
proved rational current lower=1/20
```

The displayed decimals are outward-rounded lower diagnostics.
The checker decides every inequality using the underlying exact
Fractions and outward rational square-root intervals.

## Parameter Containment

```text
a^2=exp(L)+t/16, N=floor(a), and L>=50 imply N>=2^25
log(2)<7/10; q>=exp[-(log 2)^2/8]>47/50; s>=exp[-log(2)(log(2)+1/N)/4]>22/25
log(3)<11/10; q>=exp[-(log 3)^2/8]>17/20; s>=exp[-log(3)(log(3)+1/N)/4]>73/100
log(5)<161/100; q>=exp[-(log 5)^2/8]>18/25; s>=exp[-log(5)(log(5)+1/N)/4]>13/25
```

The variable `s` contains `u` only. The small positive
`delta_a` factor is deliberately restored as a perturbation;
this avoids silently assuming `u-delta_a>=0`.

## Large Primes

For every `p>=5` with a two-level complete chain,

```text
For M=2, Q=1+r*z with 0<=r=p^(-1/2)*q*s<=1/sqrt(5)
J_0=1+2r^2+3r*cos(theta)
min_theta J_0=(1-r)(1-2r)>=7/5-3/sqrt(5)>1/20
```

Thus the quinary Bernstein box is an exact cross-check, while
the promoted result covers the whole `p>=5, M=2` family.

## Actual Fixed-Ray Transfer

The normalized actual coefficient factor is

```text
R_k=exp[t*k*log(p)*delta_a/2]*(1+d_(m*p^k))/(1+d_m)
|R_k-1|<r_x=18000/x
|J_ang-J_0|<576*r_x+288*r_x^2<900*r_x
J_ray=J_ang-rho*Im(H*conj(Q))+Im(E*conj(Q))/(log(p)*(-b))
|J_ray-J_0|<16200448/x+6486528/x^2<1/1000
```

Consequently,

```text
J_ray>1/300-1/1000=7/3000
J_ray>1/25-1/1000=39/1000
J_ray>1/20-1/1000=49/1000
```

The `pi` entering `x=4*pi*exp(L)` is the usual circle
constant inherited from the completed-zeta saddle
normalization. The angular identity uses the same ordinary
`2*pi` period of `exp(i*theta)`; no new polygonal constant is
inserted.

## Boundary

This certificate proves strict actual fixed-ray phase-current positivity only for the threshold base families (p,M)=(2,8), (3,4), and p>=5 with M=2. It does not prove length propagation, short or singleton chain control, joined p-free base compatibility, recurrent-endpoint composition, a Xi Abel gap, a successor winding bound, or RH.
