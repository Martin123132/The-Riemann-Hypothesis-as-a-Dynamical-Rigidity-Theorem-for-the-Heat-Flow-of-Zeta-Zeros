# Fixed-Root Shift-Tangency Rigidity Gate

Date: 2026-08-03

Status: exact fixed-root recurrence, EGF classification, finite Hankel obstruction, and fixed-root near-rigidity theorem; drifting-root zeta estimate open; not a proof of RH.

## Adjacent Transfer

For J_A^(d,n)(z)=sum_(j=0)^d binom(d,j)A_(n+j)z^j, (J_A^(d,n))^(s)(r)=d!/(d-s)! sum_(j=0)^(d-s)binom(d-s,j)r^jA_(n+s+j).

Pascal's identity gives J_A^(d+1,n)=J_A^(d,n)+zJ_A^(d,n+1).

Let r!=0. If r has multiplicity at least m in J_A^(d,n) and at least m+1 in J_A^(d+1,n), then J_A^(d,n+1)=(J_A^(d+1,n)-J_A^(d,n))/z has multiplicity at least m at the same r.

The existing polar-contact lemma supplies the m-to-(m+1) lift only when the adjacent degree-(d+1) window is negative-root hyperbolic. The same-root shift transfer is therefore conditional on that higher-degree hypothesis; the recurrence theorem itself merely assumes the common root.

The condition r!=0 is essential. Division by z can lower multiplicity at r=0, so the zero-root channel remains separate.

## Fixed-Root Classification

If one r!=0 is a root of multiplicity at least m in every shifted degree-d Jensen window, then with M=m-1 and p=d-m+1, (1+rE)^p A_(n+M)=0 for every n>=0.

For d=4,m=2, the recurrence derived without relabelling is A_(n+1)+3rA_(n+2)+3r^2A_(n+3)+r^3A_(n+4)=0. Writing A_n through A_(n+3) requires an explicit shift of the recurrence index.

With q=-1/r, every recurrence solution has A_(M+k)=q^kP(k), deg(P)<=p-1=d-m. The lower derivative conditions determine the finite prefix consistently.

For F(z)=sum_(n>=0)A_nz^n/n!, the all-shift root equations are (1+rD)^(d-s)F^(s)=0 for 0<=s<m. In particular F(z)=exp(qz)Q(z), q=-1/r, with deg(Q)<=d-m.

Conversely, F=exp(qz)Q with deg(Q)<=d-m gives A_n=q^nP(n) with the same degree bound, and every J_A^(d,n) has r=-1/q as a root of multiplicity at least m.

The rational OGF sum A_nz^n has a finite pole for a nonzero geometric-polynomial sequence, but the entire zeta generating function is the EGF F=sum A_nz^n/n!. The correct zeta obstruction is the exponential-polynomial classification, not OGF entireness.

## Persistent Rational Witness

A_n=(33n^2+3n+4)/(4*10^n) has the exact common quartic double root r=-10 and satisfies the corrected cubic recurrence at every shift.

sum A_nz^n=(1+7z/10+17z^2/200)/(1-z/10)^3, whose numerator equals 33/2 at z=10.

sum A_nz^n/n!=(1+9z/10+33z^2/400)exp(z/10), which is entire and is exactly the classified exponential-polynomial form.

Every contiguous order-four Hankel determinant of this rank-three sequence vanishes.

The geometric-quadratic class is invariant under A_n'=2(2n+1)A_(n+1), with q'=4q^2, a'=10qa, b'=8qa+6qb, and c'=2q(a+b+c).

This exact witness validates the recurrence and the normalization correction. The separately reported canonical continuation and its A_10 value are not promoted here because its native source data were not supplied.

## Finite Hankel Bridge

Let c_p(r)=(binom(p,0),binom(p,1)r,...,r^p)^T and H_(p+1,M+n)=[A_(M+n+i+j)]_(i,j=0)^p. The vector of p+1 consecutive recurrence defects is H_(p+1,M+n)c_p(r).

If p+1 consecutive fixed-root recurrences are exact, H_(p+1,M+n)c_p(r)=0. Since the first component of c_p(r) is 1, det H_(p+1,M+n)=0.

Let p=d-m+1. If J_A^(d,n) has one nonzero multiplicity-m root and the p consecutive extensions J_A^(d+1,n),...,J_A^(d+1,n+p-1) are negative-root hyperbolic, repeated polar lifting and adjacent transfer give the same multiplicity-m root in shifts n through n+p. Therefore det H_(p+1,m-1+n)=0. At a first loss of global Jensen hyperbolicity, every fixed higher-degree window remains hyperbolic by closure, so a nonzero determinant excludes that degree-d collision.

For d=4,m=2, p=3 and M=1. Four consecutive fixed-root double tangencies force H_(4,n+1)=det[A_(n+1+i+j)]_(i,j=0)^3=0.

The certified anchor theorem H_(4,k)(-100)>0 for every k>=0 excludes four consecutive quartic double tangencies at one fixed nonzero root at lambda=-100.

The certified forward-invariance theorem strengthens the anchor to H_(4,k)(lambda)>0 for every k>=0 and every lambda in [-100,0]. Hence the same fixed-root four-chain is impossible throughout that heat interval.

At any lambda in [-100,0], if a quartic double root at shift n has negative-root hyperbolic quintic extensions at shifts n,n+1,n+2, polar lifting and adjacent transfer create the four forbidden quartic tangencies. The initial tangency and all three higher-degree hypotheses therefore cannot coexist anywhere on the interval.

## Quantitative Near-Rigidity

For epsilon_n=(1+rE)^pA_(M+n), q=-1/r, and B_n=q^(-n)A_(M+n), one has Delta^pB_n=(-1)^p q^(-n)epsilon_n exactly.

Let P_(p-1)(k)=sum_(s=0)^(p-1)binom(k,s)Delta^sB_0. Then B_k-P_(p-1)(k)=sum_(j=0)^(k-p)binom(k-1-j,p-1)Delta^pB_j.

Consequently |B_k-P_(p-1)(k)|<=binom(k,p)max_(0<=j<=k-p)|q|^(-j)|epsilon_j|. For quartic double roots the amplification factor is binom(k,3).

For the p+1 defect vector epsilon=Hc, ||epsilon||_2>=sigma_min(H)||c||_2>=|det H| ||c||_2/||H||_2^p. Hence max_i|epsilon_i| is at least this quantity divided by sqrt(p+1).

If row i is exact at r_i but tested at r_*, its residual is sum_(j=1)^pbinom(p,j)(r_*^j-r_i^j)A_(M+n+i+j). On |r_i|,|r_*|<=R, use |r_*^j-r_i^j|<=jR^(j-1)|r_*-r_i| and compare with the Hankel-margin lower bound.

A zeta-specific near-chain theorem still needs certified singular-value or determinant/norm margins on the relevant heat interval, a compact annulus excluding r=0 and infinity, and control of the adjacent-root drift. Finite closeness to an exponential-polynomial does not contradict entireness by itself.

## Route Decision

Promote the fixed-shift theorem in EGF normalization and use finite Hankel margins before invoking an infinite-chain pole argument.

The existing polar cascade controls a fixed shift while degree grows. This gate controls a fixed degree while shift grows. For a non-exponential-polynomial zeta source, an exact obstruction cannot remain in either fixed-root axis indefinitely.

Build interval lower bounds for sigma_min(H_(4,n)) or a determinant-plus-norm substitute on the heat interval where tangency propagation is needed. Combine those margins with the root-drift residual and the adjacent quintic transfer hypothesis.

This is a Jensen/PF branch. It does not alter the distinct Section 11.193 carrier-plus-near-kernel obligation or license mixing their estimates.

## Pi Provenance

This gate introduces no pi. Every identity is binomial, finite-difference, rational-function, polynomial, differential-operator, or Hankel algebra.

## Proof Boundary

No approximate-chain exclusion on the zeta heat interval, all-degree or all-shift hyperbolicity theorem, PF-infinity, Lambda<=0, RH, or prize-level conclusion is proved. This gate is not a proof of RH.
