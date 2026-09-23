# Uniform filtered first moment of the literal rough source

23 September 2026. Private TWO HANDS NETWORK LTD research.
Analytic deduction at the explicitly inherited reciprocal-phase interface.
Not a proof of the favourable signed current or RH. Not machine-formalized.

## 1. Objects, quantifiers and new conclusion

Keep every object in incoming/Nonlinear Increment Separation (full dated
filename in INPUTS.json): fixed selected set S, rho_S=prod_(p in S)(1-1/p),
I0=[1/8,2^24], a(x)>0, frozen c_x, both original windows, d=ceil log(N+1),
A=N+(1+2nu)/4, H comparable to N^(2/3), delta=log(d)/H and

 P_x(s)=sum_(n<=N, (n,Q_S)=1) u_n(x) exp(-is log n),
 u_n=c_x(log(A/n)/d)/(a(x)sqrt(rho_S d n)).

Here s in P_x(s) is a REAL physical height; x remains the original profile.
No source truncation, prime retuning, gate replacement or carrier change occurs.
The fixed compact profile has the inherited uniform coefficient/BV bounds.

Write Kf(t)=int_R k_delta(v)f(t+v)dv for the same positive unit-mass sinc^4
kernel, including its whole infinite tail. Its transform is b(|xi|/delta),
where b(r)=1-6r^2+6r^3 on [0,1/2], 2(1-r)^3 on [1/2,1], and zero beyond 1.
Let

 e(x)=int_0^1 c_x(y)^2/a(x)^2 dy,   e_*=max_(x in I0)e(x)<infinity.

NEW CONCLUSION, under the reciprocal-phase estimate stated in Section 2:

 epsilon_N := max_original_windows sup_(x,t in window) K|P_x|(t) -> 0.       (1)

Meanwhile K|P_x|^2(t)->e(x), uniformly over the same base heights/profiles.
There is no contradiction: the filtered amplitudes are not uniformly
integrable in their squares. We give no effective N threshold or rate in (1).

Unlike the incoming overlap theorem, (1) does not assume the old original-
window L1 theorem or its full-to-rough L2 reassembly. It is proved directly
for the literal S-rough indices using a stronger twisted version of the
incoming pointwise quadratic-energy argument.

## 2. Fixed rational twists of the pointwise filtered energy

For any FIXED coprime positive integers a,b,

 K[exp(is log(a/b)) |P_x(s)|^2](t)
 =sum_(m,n rough <=N) u_m u_n exp(it log(am/(bn)))
                         b(|log(am/(bn))|/delta).                          (2)

This is finite expansion plus the exact full-line transform. The sign and
index convention in (2) follow from P*conjugate(P); no Fourier series of
the nonlinear source is invoked.

The zero-frequency terms satisfy am=bn, hence n=ar, m=br. If (ab,Q_S)>1
there are no such terms. Otherwise they are

 D_(a,b),N(x)=sum_(r<=N/max(a,b), rough) u_(ar)(x)u_(br)(x)
            -> e(x)/sqrt(ab),                                             (3)

uniformly in x. To check (3), insert the coefficients. The two y arguments
differ by fixed log(a)/d and log(b)/d. Bounded profile derivatives replace
their product by c_x(log(A/r)/d)^2 at cost O_(a,b,I0,S)(1/d).
Fixed-period rough harmonic summation replaces sum_(rough r) f(log(A/r)/d)/r
by rho_S times its integral, with an O_(S,I0,a,b)(1) unnormalized error.
The upper truncation shifts the limiting y endpoint by O_(a,b)(1/d).
Dividing by rho_S d gives (3). All constants here may depend on the fixed S.

Off diagonal put h=am-bn!=0. On a dyadic scale D<=m<2D the filter forces
|h|<=C_(a,b) delta D, n comparable to D. A nonempty scale consequently has
D>=c_(a,b)/delta ->infinity. The congruence h=am mod b and both roughness
restrictions split the m sum into a fixed finite union of progressions
modulo b Q_S. For a fixed h the phase is

 t log(am/(am-h)) = -t log(1-h/(am)).

After scaling to a dyadic progression, its normalized finite-order
derivatives approach those of a signed reciprocal phase with errors
O_(a,b,S)(|h|/D+1/D)=O_(a,b,S)(delta). Either sign of h is allowed.
The magnitude of the scale is comparable to t|h|/D, t comparable to N^2.
Cutoffs n<=N, m<=N and filter support only cut partial intervals.

Use the SAME inherited stable partial-progression exponential-sum estimate
with (kappa,lambda)=(13/31,16/31), for every fixed xi>0:

 |partial sum| <= C_(a,b,S,xi) N^xi t^kappa |h|^kappa D^(lambda-2kappa).

This is the sole non-elementary arithmetic interface used here. Its
partial-interval and finite-derivative hypotheses are not replaced by a
sampling claim. The standard framework is ANTEDB, Exponent Pairs,
Definition 5.1 and Lemma 5.2, rechecked 23 September 2026:
https://teorth.github.io/expdb/blueprint/exponent-pairs-chapter.html
The application to these fixed progressions is retained explicitly, as in
proofs/GRAM.md. Constants are NOT asserted uniform in a,b or growing S.

The product u_m u_((am-h)/b) has sup plus total variation O(1/(dD)).
For either sign of h, |log(am/(am-h))| decreases monotonically with m.
Thus the filter multiplier has variation <=1. Abel summation costs no
delta^-1. Summing |h|<=C delta D gives

 O(N^xi t^kappa delta^(kappa+1) D^(lambda-kappa)/d)

per dyadic block. Since lambda-kappa=3/31>0 the scales sum geometrically.
There are no off-band pairs, since their exact multipliers are zero. In total

 offdiag = O_(a,b,S,I0,xi)(N^(-1/93+xi)(log d)^(44/31)/d)=o(1)             (4)

on choosing xi=1/186. In particular, uniformly in the original base heights
and profiles,

 K[exp(is log(a/b))|P|^2](t) -> e(x) r_S(a/b),                            (5)

where r_S=1/sqrt(ab) if (ab,Q_S)=1 and r_S=0 otherwise. For a=b=1,
(5) includes incoming pointwise energy and gives a uniform O(1) bound.

## 3. A finite-prime energy law, not phase independence

Let J be any FIXED finite test-prime set disjoint from S. These primes are
proof instruments only; they are never inserted into the physical multiplier
or removed from P. Put theta_p(s)=-s log p and

 G_J(theta)=prod_(p in J) (1-1/p)/|1-p^(-1/2)exp(i theta_p)|^2.

Use normalized Haar measure on the finite torus for S union J. For every
fixed continuous f on that torus,

 K[f(theta(s)) |P_x(s)|^2](t)
       -> e(x) int f(theta) G_J(theta_J)dtheta,                          (6)
 K[f(theta(s))](t) -> int f(theta)dtheta.                                 (7)

For a character, (6) follows from (5) and the Poisson coefficients
p^(-|k|/2). A nonzero selected-S character has coefficient zero because
the corresponding rational diagonal is forbidden by roughness. Every
nontrivial character in (7) is exactly filtered out once delta is smaller
than its fixed, nonzero rational logarithmic frequency. Trigonometric
approximation extends both formulas uniformly: use positivity and unit mass
for (7), and the uniform K|P|^2 bound for (6). There is no need to bound a
source Fourier l1 norm or use a growing approximation degree.

Equations (6)-(7) concern quadratic energy and finite-prime observations.
They do NOT assert independence of arg(P), any nonlinear Q_l, or the
original Euler multiplier. There is no joint random replacement of P.

## 4. Hellinger bound and the ordered limit

For fixed J, G_J is positive continuous and bounded away from zero.
Positive-kernel Cauchy gives the exact inequality

 K|P| <= [K(|P|^2 G_J^(-1/2)) K(G_J^(1/2))]^(1/2).                       (8)

Both factors are covered by (6)-(7); therefore

 limsup_N epsilon_N <= sqrt(e_*) eta_J,
 eta_J=int sqrt(G_J)dtheta
      =prod_(p in J) sqrt(1-1/p) sum_(k>=0) g_k^2 p^(-k),
 g_k=binom(2k,k)/4^k.                                                     (9)

The constant-term formula in (9) follows by the absolutely convergent
binomial expansion at each fixed p. Set h_k=(1/4)_k/k!. The recurrence gives

 g_(k+1)^2/g_k^2 = (k+1/2)^2/(k+1)^2
 <= (k+1/4)/(k+1)=h_(k+1)/h_k,

because (k+1/4)(k+1)-(k+1/2)^2=k/4>=0. Hence g_k^2<=h_k and

 eta_J <= prod_(p in J)(1-1/p)^(1/4).                                    (10)

Choose J={p<=R : p not in S}, with R fixed BEFORE N tends to infinity.
The finite Euler product includes the harmonic sum H_floor(R), because
every integer <=R has all prime factors <=R. Thus

 eta_J <= [rho_S H_floor(R)]^(-1/4) ->0 as R->infinity.                   (11)

The actual product over removed selected primes up to R is >=rho_S;
this direction gives (11) whether or not R includes every selected prime.
First choose R to make (11) small, then choose N for the finitely many
approximations in (6)-(8). This proves (1), with no J(N), no effective
decay, and no constants uniform in a growing prime set. The very small
fixed rho_S is not suppressed: it can make the onset enormous.

## 5. Original-window overlap with its endpoints retained

Let j(t)dt=dr/Delta on either original window I, zero elsewhere. Since
t=2pi(N+r)^2, j(t)=1/[4pi(N+r)Delta] and
min_I j/max_I j >=1/2 eventually. For every s in I, the difference interval
I-s contains either [0,H/2] or [-H/2,0]. The kernel is even and

 int_(|v|>H/2) k_delta(v)dv <=512/[pi(delta H)^3].                       (12)

Indeed k_delta(v)<=96/[pi delta^3 v^4] and integrate that majorant.
Since delta H=log d ->infinity, eventually int_0^(H/2) k_delta>=1/4.
Consequently, INCLUDING s at either endpoint,

 int_I j(t)k_delta(s-t)dt >= j(s)/8.                                     (13)

Tonelli, positivity, (1) and the original finite profile mass A_I<9 give

 m_N=int (a/4)dx int_I j(s)|P_x(s)|ds
     <=8 int (a/4)dx int_I j(t)K|P_x|(t)dt
     <=8 A_I epsilon_N.                                                   (14)

No exterior mass is subtracted or assumed small. Unlike an approximate
convolution identity, (13) is a one-sided endpoint-inclusive comparison.
It rederives the original ungated first moment without reassembly.

For incoming original-to-shift overlap,

 O_N=int dnu(t,x) |P_x(t)| K|P_x|(t)
     <=epsilon_N m_N <=8 A_I epsilon_N^2.                                (15)

Thus the incoming complete nonlinear additivity error has the stronger
explicit bound h_l epsilon_N m_N, or 8 h_l A_I epsilon_N^2,
with h_1=5/2, h_3=7. It is o(m_N) when m_N>0, but no effective N rate is
claimed. These are comparison inequalities, not a numerical gain at a
specified physical height.

## 6. Two independently shifted heights at the same base

For fixed base t in either original window, draw shifts v,w independently
from the SAME full kernel. Independence here means only product integration;
it is not an assumption about arithmetic phases. Exactly,

 int int k(v)k(w)|P(t+v)||P(t+w)|dvdw = (K|P|(t))^2 <=epsilon_N^2.

The global zero-safe C1,1 estimate for the original odd maps therefore gives

 int int k(v)k(w)|Q_l(P(t+v)-P(t+w))-Q_l(P(t+v))+Q_l(P(t+w))|dvdw
       <=h_l epsilon_N^2.                                                (16)

Each marginal squared amplitude still tends to e(x), not to zero.
This extends the nonlinear splitting to two shifts without requiring their
base points t+v,t+w to lie in the original window. Only the common centre t
is restricted; (6) integrates both full tails before applying (8).

## 7. What this changes in the signed current, and what it does not

Retain I_sym,N, the original D_=,N, global subtraction, every gate harmonic,
both complex multipliers, and both endpoint densities exactly as incoming
Sections 6-7. Its finite representation error now obeys

 |D_filter,N^[d]-I_inc,N|
 <= (4/ZW)[alpha h_1 K_1(0)+beta h_3 K_3(0)] epsilon_N m_N.                (17)

The multiplier tail, exact atom and base-phase linearization are unchanged:

 D_nonzero,N = I_sym,N-D_=,N
              +O_(I0,S)(epsilon_N^2)+O_(I0,S,B)(d^-B)+O(N^-2/3).          (18)

The limiting exact reserve 7e-18 is still used ONCE. There is no second
diagonal subtraction or additional gate normalization. All original
restoration obligations remain. The target is still

 limsup_N max_original_windows [I_sym,N-D_=,N]
                       <0.000920279889999993.                            (19)

Equation (19) is OPEN. Neither (1) nor (16) makes a quadratic response small:
for example a scalar variable equal to L with probability L^-2 and zero
otherwise has mean 1/L but square mean 1. This generic guard is not an
actual-source counterexample; it explains the invalid logical shortcut.
Oddness, external centering and small endpoint overlap do not, on their own,
bound the signed alignment of the external coefficient and the complete
nonlinear increment. No bound on that alignment is supplied here.

The next honest target is that signed, energy-bearing alignment, after
paying the now uniformly controlled additivity error. Any further spectral
or arithmetic cancellation must apply to the joint nonlinear response,
not only the fixed-prime quadratic-energy law (6).
