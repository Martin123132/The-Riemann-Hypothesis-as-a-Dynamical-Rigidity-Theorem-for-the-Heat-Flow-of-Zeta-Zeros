# Jensen-Window PF Suzuki Fixed-Omega Phase Diagram

Date: 2026-07-23

Status: internally audited theorem-candidate note with exact
countermodel guards and one open arithmetic target. It is not a proof
of the determinant premise, RH, or `Lambda <= 0`.
Independent expert review is required before publication.

```text
work/rh_compute/results/jensen_window_pf_suzuki_fixed_omega_phase_diagram.json
python work/rh_compute/scripts/jensen_window_pf_suzuki_fixed_omega_phase_diagram.py
python work/rh_compute/scripts/check_jensen_window_pf_suzuki_fixed_omega_phase_diagram.py
```

## Fixed-Omega Equivalence

For `omega>0`, let

```text
D(omega): for one (equivalently every) integer nu with nu*omega>1, det(I+/-K_(omega,nu)[t])!=0 for every t>=0
Theta_omega(z)=xi(1/2-omega-i*z)/xi(1/2+omega-i*z).
```

Combining the audited causal-multiplier direction with Suzuki's
fixed-parameter Phragmen-Lindelof and strict-compression arguments gives

```text
D(omega) iff Theta_omega is meromorphic inner in C+ iff m_xi(rho-2*omega)>=m_xi(rho) for every xi-zero rho with Re(rho)>1/2+omega
```

This combined equivalence is a corpus theorem candidate. Its directions
are kept explicit below.

## Determinants To Innerness

All-time determinant nonvanishing is equivalent to strict contraction
of every compact self-adjoint truncation. Compatible compactly supported
forms therefore extend to a bounded full-line Hankel form. Reflection
turns it into a bounded causal translation-invariant convolution
operator, whose multiplier is `Theta_omega^nu` in the high half-plane.
Hence

```text
All finite strict contractions extend to a bounded causal convolution multiplier equal to Theta_omega^nu, hence Theta_omega is meromorphic inner
Theta_omega^nu is inner for one positive integer nu iff Theta_omega is inner iff Theta_omega^mu is inner for every positive integer mu
```

Real-boundary unimodularity is exact:

```text
Theta_omega(u)=conj(xi(1/2+omega-i*u))/xi(1/2+omega-i*u), so |Theta_omega(u)|=1 wherever defined
```

## Innerness To Determinants

Suzuki's earlier fixed-omega argument observes that the quotient is
uniformly bounded in a sufficiently high strip. If it has no poles in
`C+`, Phragmen-Lindelof propagates that bound down to the real boundary.
Thus

```text
For the xi quotient, no poles in C+ plus Suzuki's high-strip bound and Phragmen-Lindelof imply Theta_omega in H-infinity(C+)
```

Once the quotient is inner, its boundary multiplier is an `L2` isometry.
Suzuki's noncompact-support lemma excludes equality after finite
compression, while `nu*omega>1` makes each truncation compact. Therefore

```text
Theta_omega inner gives a full L2 isometry; Suzuki's noncompact-support lemma and compact truncations give ||K_(omega,nu)[t]||<1 for every finite t
```

The growth input is indispensable: `exp(-i*z)` has no poles and has
modulus one on the real line, but grows like `exp(Im(z))` in `C+`.

## Exact Zero-Pair Rule

For `rho=beta+i*gamma`,

```text
rho=beta+i*gamma maps to z_rho=-gamma+i*(beta-1/2-omega); the numerator there is xi(rho-2*omega)
ord_pole(z_rho)=max(0,m_xi(rho)-m_xi(rho-2*omega))
```

Consequently a fixed shift can hide an off-line zero only through an
equal-height horizontal zero pair separated by `2*omega`.

## Phase Diagram

Define

```text
C_xi={(Re(rho)-Re(rho'))/2>0: xi(rho)=xi(rho')=0 and Im(rho)=Im(rho')}; C_xi is countable
delta_xi=sup_{xi(rho)=0}(Re(rho)-1/2) in [0,1/2]
```

Then

```text
[delta_xi,infinity) subset D subset [delta_xi,infinity) union (C_xi intersect (0,delta_xi))
omega not in C_xi => [D(omega) iff xi(s)!=0 for Re(s)>1/2+omega].
```

The possible subthreshold successes are therefore countable and
entirely attributable to exact zero cancellation. More strongly,
isolation of any one off-line zero gives

```text
If RH is false, one off-line zero rho and zero isolation give epsilon_rho>0 such that D(omega) fails for every 0<omega<epsilon_rho
RH iff 0 belongs to closure(D) iff D contains a sequence omega_n->0
```

## Fixed-Shift Countermodels

The cancellation guard is realized exactly by

```text
F(z)=(z^2+1)*(z^2+9)
F(z-i)/F(z+i)=(z-4i)/(z+4i)
```

The quotient is a Blaschke factor although `F` has zeros at
`+/-i` and `+/-3i`. Raising the upper pair's multiplicity breaks the
condition:

```text
F_def(z)=(z^2+1)*(z^2+9)^2
F_def(z-i)/F_def(z+i)=(z-4i)^2*(z+2i)/((z-2i)*(z+4i)^2)
The unmatched multiplicity leaves a pole at z=2i.
```

These polynomial examples test the cancellation logic. They are not
surrogates for Suzuki's decaying arithmetic kernel.

## Jordan-Totient Scalar Target

Suzuki also gives a scalar arithmetic formulation. Put

```text
c_omega(n)=n^omega*product_(p|n)(1-p^(-2*omega))>0
g_omega^<1>(x)=integral_x^1 sqrt(y/x)*g_omega(y)dy/y,
h_omega^<1>(x)=x^(-1)*sum_(n<=x)c_omega(n)*g_omega^<1>(n/x)
```

where `g_omega` is Suzuki's explicit beta-integral weight. Then

```text
Theta_omega inner iff x^(-1/2)*1_(1,infinity)(x)-h_omega^<1>(x) belongs to L2(1,infinity)
If h_omega^<1>(x) has one sign for all sufficiently large x, then Theta_omega is inner
```

Finite samples do not certify it: both conditions concern the full
unbounded half-line.

The sharp scalar obligation is therefore

```text
Prove the L2 residual, or the stronger eventual one-sign condition, for a sequence omega_n->0 by direct arithmetic bounds without assuming a zero-free half-plane
```

Although every `c_omega(n)` is positive, the weight is signed. For
`omega=1/2`, its explicit primitive is negative near zero and positive
elsewhere, so coefficient positivity alone cannot establish the target.

## Sources

- Masatoshi Suzuki, `Hamiltonians arising from L-functions in the Selberg class`: https://arxiv.org/abs/1606.05726
- Masatoshi Suzuki, `A canonical system of differential equations arising from the Riemann zeta-function`: https://arxiv.org/abs/1204.1827
- Masatoshi Suzuki, `On monotonicity of certain weighted summatory functions associated with L-functions`: https://arxiv.org/abs/1204.1823
- Chris Guiver, Hartmut Logemann, and Mark R. Opmeer, causal `L2` multiplier theorem: https://link.springer.com/article/10.1007/s00498-024-00387-4

## Proof Boundary

This artifact identifies the fixed-omega Suzuki determinant property with meromorphic innerness and an exact shifted-zero multiplicity condition, subject to independent review of the combined operator argument. It does not prove that the condition holds at any subcritical omega, does not prove eventual sign or the L2 summatory residual, and does not prove RH or Lambda<=0.
