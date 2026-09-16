# An explicitly evaluated switching correction on the actual height trajectory

**16 September 2026 · v1**  
**Martin Ollett / TWO HANDS NETWORK LTD**  
**Copyright (c) 2026 TWO HANDS NETWORK LTD. All rights reserved.**

## 1. Result and the policy distinction

The incoming package represents the correction caused by an adaptive upper-index
selector as a signed sum of centered primitives at its switches. It bounds every
primitive uniformly, but correctly does not multiply that bound by an unproved
small switch count. Its eight-term finite policy and all historical results are
preserved unchanged.

This continuation introduces an explicitly declared **phase-only reference
selector** alongside that policy. It is not the eight-term argmin, nor asserted
to approximate that argmin with vanishing regret. Its purpose is to calculate a
signed, adaptive ordinary correlation, including its full accumulated switching
correction, using the actual shared-height arithmetic.

For this reference rule, on the declared degree-edge sequence, the new written
argument proves, uniformly over the two original windows,

    <X_N V_(i_phi)>_W -> C_phi,
    C_phi in [1.446119963280524, 1.446120051051692],

and

    Gamma_phi -> -0.01789414360498004329... < 0.

In particular the adaptive ordinary contribution is eventually positive. No
starting cutoff or effective convergence constant is supplied. The switch count
itself is asymptotic to (8/25) N^(2/3) log N, so this result is not obtained by
assuming few switches or vanishing total selection bias.

The finite full-source check uses the SAME original N=64 cells, central gate,
all six coefficients, source-sized scales, degree five, and complete payment
0.271900134725. The reference choice is enclosed on each complete cell, including
both choices on cells containing a possible tie. Every full coefficient remains
in the selection-regret payment. Restored surplus floors are >1.0815 and >1.0512.
Separately paying the ordinary contribution, gate correction, lower-block loss
and regret still gives >1.0553 and >1.0208. These are alternative enclosures of
the same historical surplus, not additional certified heights or counts.

The unbounded GATED minimum remains unresolved. The new theorem concerns the
ordinary product with a specified phase-adaptive index. The reference's gate,
lower-block loss, and selection regret cannot be inferred to have favorable
limits from that ordinary product. An exact transfer back to the original
short-sum selector is given in Section 8.

## 2. Unchanged source and the new reference rule

Retain, for N>=64,

    J=floor((N/64)^(1/3)), M=N-J-1, d=ceil(log(N+1)),
    ell=1/[50(J+1)], a=N+r, t=2*pi*a^2,
    I_(N,nu)=[(1+2nu)/4-ell/2,(1+2nu)/4+ell/2], nu=0,1,
    chi=theta0(t)+pi*(N+nu),
    theta0(t)=2*pi*a^2*log(a)-pi*a^2-pi/8,
    W=(1+cos chi)^2, <f>_W=E_I(Wf)/E_I W.

The original polynomial and positive scales are

    B_(d,j)(z)=sum_(k<=j) binom(j,k)/binom(d,k)*z^(2k)/(4^k*(2k)!),
    s_j=B_(d,j)(d),
    beta_j=Re sum_(n<=M) sqrt(a/n)*exp(i*(chi-t*log n))
                                      *B_(d,j)(log(a/n)+i*pi/4),
    c_N=H_J-E_N, H_J=(2/3)(J+1)-503/1000,
    E_N=(e^2/9)*N^(-(log 16-5/2))+(log N+2)/(4*pi*N),
    Z_j=(beta_j+c_N)/s_j, X=Z_0/sqrt(N), V_j=(2s_d/N)Z_j.

The exact selection surplus is S_N=<X_+ min_j V_j>_W. Its positivity has the
inherited sufficient full-width sign interpretation, not an absolute zero
count. The upper block is j>=ceil(d/4); its complement is the unchanged lower
block. Every field above remains in the finite and analytic accounting.

Choose two upper positions

    j_A(N)=ceil(11*d/32),    j_B(N)=d.

Both belong to the upper block; they are distinct for d>=5. On
N_m=floor(exp(m)), d=m+1. Their limiting spectral parameters in the inherited
relation lambda=s/2+4s^2-2s^3 are respectively s_A=1/4 and s_B=1/2.
Define the FIXED complex pilot constants

    b_A=exp(1/4+i*pi/16), b_B=exp(i*pi/8), z=b_A-b_B != 0.

They are the n=1 constants of the limiting spectral kernels. They are NOT
asserted to equal the exact finite-N first-term amplitudes. This source-derived
reference rule was fixed before the new finite phase calculation:

    i_phi(r)=j_A if Re(exp(i*chi(r))*z)<=0, else j_B.       (P1)

The rule compares two phase projections only. It neither reads full source
values nor inspects whether the final surplus passes. It leaves the physical
coefficient vector, original W, and exact positive-central gate unchanged.
The old eight-term rule is retained separately, not overwritten.

Write z=R exp(i*alpha), R>0. Its first mask is the fixed periodic step function

    I_A(chi)=1_{cos(chi+alpha)<=0}, I_B=1-I_A.

The equality convention affects a zero-measure set, since chi is strictly
increasing on each window.

## 3. A phase-mask mixed-moment theorem

This is the step that replaces a bound on the NUMBER of switches. The analytic
inputs are the explicitly inherited central-energy and upper-kernel estimates
in the included central-energy derivation. In particular, for Lambda=1+log N,

    ||X||_(L2(W))=O(sqrt Lambda),
    |b_(n,j,N)|<=C n^(-beta), |partial_r b_(n,j,N)|<=C n^(-beta)/N,
    ||V_j-V_j^[K]||_2
      <=C[K^(1/2-beta)+N^(2/3-beta)*sqrt Lambda],           (H1)

uniformly for j/d>=1/4 and beta=7/10. Here

    b_(n,j,N)=(2s_d/N)*sqrt(a/n)*B_(d,j)(log(a/n)+i*pi/4)/s_j,
    V_j^[K]=Re(exp(i chi) sum_(n<=K)b_(n,j,N)exp(-it log n))
                                             +(2s_d/N)c_N/s_j.

The full central sum is never shortened. On the degree-edge sequence, for each
fixed n, the inherited kernel limit is

    b_(n,j_A,N) -> b_A*n^(-3/4),
    b_(n,j_B,N) -> b_B*n^(-1),                            (H2)

uniformly over the windows. The uniform envelopes in H1 make the diagonal
series after multiplication by n^(-1/2) summable. Watt's unconditional
mean-square theorem remains an external dependency of the central energy bound;
it is not reproved here or used as a numerical finite-N upper bound.

For a fixed bounded periodic step mask I with finitely many intervals, let

    mu_0(I)=(1/(2*pi)) integral_0^(2*pi) W(chi) I(chi) dchi,
    mu_2(I)=(1/(2*pi)) integral_0^(2*pi) W(chi) I(chi) exp(2i chi) dchi.

For either spectral position (and, by the same proof, any fixed position with
s>1/6), the theorem is

    <I(chi) X V_j>_W ->
      [mu_0(I) Re(b) zeta(1+s)+Re(mu_2(I)b)]/3.          (M1)

Here s is the spectral parameter, either 1/4 or 1/2, and
b=exp(1/2-s+i*pi*s/4). The symbol s_j elsewhere remains the positive
polynomial normalizer. No special-function value in the critical
strip occurs in the limiting constant; zeta(1+s) is an absolutely convergent
positive real series.

### 3.1 Do not approximate a changing indicator by a fixed Fourier polynomial

Set K=min(M,ceil(Lambda^6)). Multiplication by |I|<=1 and H1 bound the first
replacement inside the ordinary mixed moment by

    ||X||_2 ||V_j-V_j^[K]||_2
       =O(Lambda^(-7/10)+N^(-1/30)*Lambda).              (M2)

Let I_H be the Fejer mean of I of degree H=ceil(Lambda^6). Its positivity and
unit mass give 0<=I_H<=1. The fixed step mask's Fourier coefficients satisfy
|I_hat(k)|<=C/(1+|k|). Parseval, including the multiplier
1-|k|/(H+1), proves directly

    mean_circle |I-I_H|^2 <= C/H,
    sum_(|k|<=H) |I_H_hat(k)| <= C(1+log H).             (M3)

Indeed the low-frequency squared error is bounded by
C sum_(1<=k<=H) (k/(H+1))^2/k^2=O(1/H), and the high-frequency square tail
is O(sum_(k>H)k^(-2))=O(1/H). Values at the finitely many jumps do not affect it.

The SAME L2 error is valid, up to an absolute factor, along the actual r-window.
To check this without assuming a sampling distribution, change variable to chi.
Its derivative 4*pi*a*log(a) is positive; its maximum/minimum ratio is bounded,
and the total phase range covers many complete periods. For every nonnegative
periodic function f, bound the inverse Jacobian by its maximum and cover the
range by complete periods. This gives

    E_I f(chi(r)) <= C mean_circle f.

Together with W<=4 and E_I W bounded below, it transfers M3 to L2(W). The bound
is independent of H, so it licenses this increasing degree.

The short upper coefficient has sup norm O(K^(1-beta)). Therefore

    |<X V_j^[K](I-I_H)>_W|
      <=||X||_2 sup|V_j^[K]| ||I-I_H||_2
      =O(Lambda^(1/2) K^(3/10) H^(-1/2))
      =O(Lambda^(-7/10)).                               (M4)

This is the necessary payment for the nonsmooth selector. There is no unjustified
exchange of a discontinuous mask with an unbounded central factor. A fixed-degree
Fourier approximation alone would not yield this estimate.

### 3.2 Every nonconstant harmonic is paid uniformly in the Fourier cutoff

Multiply I_H by W, and expand the two real factors X and V_j^[K]. Every phase has
form h*chi-t*log Q, with Q=n/m or nm, n<=M and m<=K. The shifted harmonic h ranges
over integers of size at most H+4. Crucially Q is still one of these actual integer
products/ratios: the selector depends only on chi and introduces NO new prime
characters or arbitrary rational powers.

The shared-height derivative is exactly

    Phi'=4*pi*a*(h*log(a)-log Q).

For every h outside {0,1}, the relevant nonconstant derivatives are bounded
below by c*N*log N, uniformly even as H grows. For product h>=2 use nm<=N*K,
K polylogarithmic; for ratio harmonics use n/m<=N and m/n<=K. Negative product
harmonics are also separated. The h=0 ratio diagonal n=m and product n=m=1 are
retained, not estimated. The remaining h=0 terms and the slow h=1 terms have
the inherited uniform mixed-sum bounds:

| Class | Sufficient cost, before the Fourier l1 factor |
| --- | --- |
| ratio h=0, n!=m | O(N^(-2/3)K^(4/5)log(2K)+N^(-1/6)K^(3/10)) |
| fast ratio/product harmonics | O(N^(-1/6)K^(3/10)) |
| ratio h=1, m=1 near cutoff | O(N^(-1/6)Lambda) |
| product h=1, nm near N | O(N^(-1/2)K^(4/5)+N^(-1/6)Lambda) |

For the last row the integer progression is retained:

    |N+r-nm| >= (1+|N-nm|)/10,
    sum_(n:N/2<=nm<=2N) min(1,B/(1+|nm-N|))
                      <=2+(2B/m)H_(2N+1), B=O(ell^(-1)).

The amplitudes are at most C n^(-1/2)m^(-beta), with absolute r derivatives
bounded by C/N times that envelope. The phase derivative has at most two
monotone pieces, independently of h. Thus integration by parts has constants
independent of the increasing harmonic cutoff. Summing the entire table against
the Fourier l1 bound O(log H) is still o(1), and is eventually bounded by
O(N^(-1/30)*Lambda). No random phases or diagonal-only lower bound is used.
The source shift terms are O(N^(-1/6)*sqrt Lambda) or smaller and are retained
until their bounded estimates remove them in the limit.

### 3.3 The two exact resonances and spectral passage

After the nonconstant terms have been paid, the ordinary integral has exactly
these leading contributions:

    (1/2) mu_0(I_H) E_I sum_(n<=K) sqrt(a/N)n^(-1/2) Re b_(n,j,N),
    (1/2) Re[mu_2(I_H) E_I sqrt(a/N)b_(1,j,N)].          (M5)

The first comes from equal-index conjugate terms. The second comes from the
n=m=1 nonconjugate term. Division by E_I W->3/2 gives the factor 1/3 in M1.
The mask moments converge to those of I by M3. For the diagonal, fix a finite
F, use H2 on n<=F, and bound the remaining sum uniformly by
C sum_(n>F)n^(-beta-1/2). Then increase F. This justifies the spectral passage
although K grows; it never treats convergence on fixed support as uniform on
all K terms. Equations M2-M5 prove M1 on both actual window families.

The error up to the prelimit diagonal expression is explicitly of the order in
M2. The subsequent spectral convergence has no effective rate supplied by the
inherited concentration theorem. Accordingly only a qualitative o(1) is claimed
for the final limiting constants.

## 4. The complete signed switching correction is explicit

For the half-circle mask I_A, direct elementary integration gives

    mu_0 = 3/4 - 2 cos(alpha)/pi,
    mu_2 = 1/8 + exp(-3i alpha)/(3*pi) - exp(-i alpha)/pi,
    p_A = mu_0/(3/2) = 1/2 - 4 cos(alpha)/(3*pi).         (A1)

For example the mask Fourier integrals are i_0=1/2,
i_1=-exp(-i alpha)/pi, i_2=0, i_3=exp(-3i alpha)/(3*pi), i_4=0.
Multiplying by W=3/2+exp(i chi)+exp(-i chi)
+(exp(2i chi)+exp(-2i chi))/4 proves A1. The other mask is its complement.

The unconditional limiting ordinary correlations are

    C_A=Re(b_A)[zeta(5/4)/2+1/12],
    C_B=Re(b_B)[zeta(3/2)/2+1/12]=c_*.

Subtract the limiting fixed-index mixture p_A C_A+(1-p_A)C_B from the selected
ordinary contribution in M1. ALL equal-index diagonal contributions cancel.
The remaining signed selection correction is

    Gamma_phi,infinity = Re[z*(mu_2-mu_0/6)]/3
      = -[(Re z)^2+4(Im z)^2]/(9*pi*|z|).               (A2)

It is strictly negative. The cancellation is arithmetic: selection changes the
second common-phase harmonic of the n=m=1 product, while the other diagonal
resonances retain their mixture weights. The source terms themselves are not
assigned independent angles.

Consequently

    C_phi = p_A C_A+(1-p_A)C_B+Gamma_phi,infinity.        (A3)

The new rational calculations enclose

    p_A = 0.10513292958332208877...,
    Gamma_phi,infinity = -0.01789414360498004329...,
    C_phi in [1.446119963280524,1.446120051051692].

The signed constant is not claimed to be a finite bound at N=64. In fact both
finite ordinary upper bounds there are below the limiting positive constant.

## 5. Why this resolves a switching sum rather than a switch-count assumption

The inherited primitive identity applies to this rule as well:

    Gamma_phi,N = sum_switches (F_out-F_in).

The rule changes at chi+alpha=pi/2+k*pi. Since chi is strictly increasing,

    number of switches = [chi(r_+)-chi(r_-)]/pi+O(1)
                        =(4+o(1)) N ell log N
                        ~ (8/25) N^(2/3) log N.         (A4)

Thus the number of switches grows substantially, and the signed sum tends to
the nonzero negative constant in A2. We have not bounded it by switch count
times the supremum of a primitive. The mask theorem evaluates the signed sum
through its actual phase harmonics, with the increasing Fourier cutoff paid.
It is consistent with the parent's vanishing primitive bound and its explicit
many-switch counter-control.

This is NOT an evaluation of the original eight-term selector's switching sum.
That selector also depends on the prime phases and its changing short sums.
Those dependencies introduce additional characters; the uniform-harmonic proof
above does not silently cover them.

## 6. Finite source calibration with every reference error paid

The finite policy compares the FIXED constants in P1, choosing indices 2 and 5.
All full coefficients and the original lower block {0,1} remain in acceptance.
For the chosen full value V_i define

    L=min(V_0,V_1), U=min(V_2,V_3,V_4,V_5),
    T=min(L,V_i), r_i=T-min(L,U)>=0.

The exact full surplus satisfies

    S_N=<X V_i>+<(-X)_+V_i>-<X_+(V_i-L)_+>-<X_+ r_i>.   (F1)

For interval evaluation use the equivalent screened regret

    r_i=max(0,max_(j in upper, j!=i) min(L-V_j,V_i-V_j)). (F2)

The self-competitor is exactly zero. This prevents charging its interval width
as a spurious selection error. The original eight-term selector's small-regret
asymptotic is NOT assigned to the new one-term reference. Its finite regret is
recomputed, and its unbounded behavior remains unresolved.

The new elementary calculation encloses Re(exp(i chi)z) on every original closed
cell. It uses a Taylor expression for theta0 about the cell center with remainder
at most [2*pi/(3*a_min)] |u|^3; theta0'''(a)=4*pi/a. This is differentiation of an
explicit phase, not of a source remainder. The two builds use 50 and 80 interval
decimal digits and have identical outward mathematical fields.

On 10 cells in the first window and 11 in the second, the phase decision interval
contains zero. Both reference choices are included on the whole cell. No midpoint
selects one branch there; no tie cell is dropped. There are 375 and 414 cells
where the phase selector is provably different from the old eight-term cell
policy. This is a declared comparison, not a relabeling of that old policy.

The resulting original-W enclosures are:

| Quantity | First window | Second window |
| --- | ---: | ---: |
| ordinary selected product | [1.17142437,1.28623200] | [1.17018480,1.28614960] |
| negative-central gate correction | [-0.04164145,-0.02734795] | [-0.13000271,-0.10441490] |
| lower-block loss | [0.01322192,0.01786276] | [0.00446477,0.00795529] |
| effective reference regret | [0.04243508,0.05652270] | [0.00628626,0.01139376] |
| directly restored surplus | >1.0815 | >1.0512 |
| separately paid F1 surplus | >1.0553 | >1.0208 |

The reference regret is larger than for the old eight-term rule, especially in
the first window. The full arithmetic enclosures pay it. These are certificates
for the same finite surplus, not evidence of new zero-free territory.

The finite covariance for the reference is enclosed by
[-0.03487302,0.00792020] and [-0.05484305,-0.02218303]. Its first-window sign is
UNRESOLVED by these finite ranges. Its negative LIMIT in A2 is not substituted
as a finite allowance or inferred by extrapolating this table.

The independent tuple/Fraction program imports no primary arithmetic or
acceptance code. It reconstructs the complete geometry, signs, gate, lower loss,
regret and all 2,048-cell sums. It shares the certified source and phase images;
it is not an independent full source-function implementation.

## 7. Constants and verification boundary

The primary constants use rational Machin arctangent bounds for pi, a positive
exponential series with its geometric tail, small-argument sine/cosine Taylor
bounds, integer square roots, and 65,536 positive terms for each of zeta(5/4)
and zeta(3/2). For exponent p>1 the unsummed tail lies between
(M+1)^(1-p)/(p-1) and M^(1-p)/(p-1). Every finite term has directed rational
root rounding. No critical-strip zeta or Gamma evaluator is called.

The separate constant calculation uses 32,768 terms, elementary mpmath.iv at
70 digits, and the direct angular moment formula in A1 rather than the reduced
negative expression in A2. Its wider intervals contain the primary ones.
The constants are numerical certification of explicit convergent expressions;
they do not prove the asymptotic transfer by sampling N.

The test suite includes complete replays, exact algebra on all coefficient
configurations in a declared finite control grid, preservation of ambiguous
branches, scale/error/geometry rejection, both precision payloads, and the
asymptotic-versus-finite distinction. Symbolic checks support the angular and
harmonic algebra. Neither finite checks nor source hashes formally verify the
analytic big-O arguments. The written mask theorem is new and awaits separate
review under its listed inherited hypotheses.

## 8. The complete original problem remains attached by an exact transfer

There are two valid ways to use the newly evaluated reference correlation.
For the reference itself, F1 and A3 make the sufficient next target

    limsup (B_phi+T_phi-D_phi) < C_phi,                  (T1)

where B_phi is the lower-block loss, T_phi its actual effective regret, and
D_phi=<(-X)_+ V_(i_phi)>. Neither T_phi=o(1) nor the strict comparison is proved.

Alternatively leave the original growing short-sum selector i_8 entirely in
place. Define the SIGNED policy-difference correlation

    H_N=<X(V_(i_8)-V_(i_phi))>_W.

Then exactly

    S_N=<X V_(i_phi)>+H_N+D_8-B_8-T_8.                  (T2)

This identity neither changes the weight nor drops the lower block. When the
original exact short-argmin policy has its previously paid T_8=o(1), a sufficient
new reference-based target is

    limsup (B_8-D_8-H_N) < C_phi.                       (T3)

This is a different accounting from R-Q<0.5172, not a larger allowance for the
same quantity. No sign for H_N is assumed. In the finite data it is negative in
the first window, [-0.04007987,-0.01944548], but positive in the second,
[0.04552654,0.07046664]. All those signed contributions are retained. Substituting
the independently widened reference pieces still gives original-policy surplus
floors >1.0520 and >1.0079. The older direct bounds remain stronger.

The remaining mathematical task is therefore a shared-height gate/lower-block/
policy-difference inequality, not another fixed-index positive moment or a
small-switch-count assumption. The phase-only result evaluates one full signed
adaptive baseline analytically. It does not complete the original eight-term
selection correction, unbounded gated selection, winding, multiplicity or
absolute zero-count closure.

## 9. Preservation and attribution

All new files are separate stage/date/version deliverables. The 18-member parent
archive is authenticated unchanged. The numerical source and eight-term pilot
are copied byte-for-byte. No laptop access, repository write, desktop change,
bookmark update or historical-pin audit is performed. Only the new reference
phase is evaluated; the complete source-function enclosures are inherited.

[D1] Included unchanged switching-primitive derivation, Sections 2–5: source,
policy correction, switching identity, and current unbounded gap.

[D2] Included unchanged central-energy/joint-spectrum derivation, Sections 4–7:
central energy, coefficient kernel envelopes, spectral convergence and complete
mixed-phase estimates. Watt (2010), DOI 10.1112/jlms/jdq024, remains the explicitly
attributed external unconditional mean-square dependency used there. This note
uses those hypotheses rather than claiming a new proof of Watt's theorem.

[S1] NIST DLMF §1.8, Fourier coefficients and Parseval,
https://dlmf.nist.gov/1.8 . The step-mask Fejer error and its actual-height transfer
are proved in Section 3, not attributed as a source-specific theorem to DLMF.

[S2] NIST DLMF §2.3(i), integration by parts,
https://dlmf.nist.gov/2.3 . The mode ranges, uniformity and near-product payments
needed for this continuation are specified above.

No general reuse permission, publication, novelty priority or RH resolution is
asserted by this private delivery. Earlier mathematical qualifications remain.
