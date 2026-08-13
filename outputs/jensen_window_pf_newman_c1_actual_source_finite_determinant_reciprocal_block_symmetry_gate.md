# Newman C1 Finite Determinant Reciprocal-Block Symmetry Gate

Date: 2026-08-05

Status: exact determinant symmetry and reciprocal-block telescope proved; physical signed arithmetic open; not a proof of RH.

## Determinant Chart

Map x=(V,mathcal N,A,Q)^T to M(x)=[[V,A],[Q,mathcal N]]. Then det M(x)=x^TJx=Vmathcal N-AQ.

For B=M(b) and D=M(d), det(B-D)-det(B)=det(D)-tr(adj(B)D)=d^TJd-2b^TJd. No inverse or nonzero determinant is required.

B=M(omega_E+p), D=M(d), and det(B)=Q(omega_E+p)=h^2P_ret at the physical anchor. The completed finite primitive is Q_geom=(2/|nu|^2){det(B-D)-det(B)}.

## Symmetry And Gauges

For real 2x2 L,R with det(L)det(R)=1, the simultaneous action M(x)->LM(x)R preserves Q and every determinant increment. This is the SL(2,R) x SL(2,R) determinant symmetry, enlarged by the determinant-product-one component.

If Delta=det(B)>0 and s=sqrt(Delta), L=sB^(-1), R=I have determinant product one and LBR=sI. In that frozen pointwise gauge, det(B-D)-det(B)=det(D')-s(D'_(11)+D'_(22)).

If Delta=det(B)<0 and s=sqrt(-Delta), take L=sB^(-1) and R=diag(1,-1). Their determinants are both -1, LBR=diag(s,-s), and the increment is det(D')-s(-D'_(11)+D'_(22)).

If det(B)=0, B^(-1) is forbidden. Rank-zero and rank-one bases remain covered by the division-free adjugate formula and must be treated as separate null-cone strata unless h^2P_ret=det(B) is proved nonzero at the calibrated witnesses.

Euclidean coordinate norms are not invariant under determinant-product-one left/right gauges. Any arithmetic estimate must either fix and justify a conditioned gauge or bound the invariant quantities det(D_j) and tr(adj(B_j)D_j), together with their tangent analogues.

## Reciprocal Blocks

Write f=f^(0)+f^(1)+f^(2), where f^(0) uses reciprocal blocks p=q, f^(1) uses |p-q|=1, and f^(2) uses |p-q|>=2. Define d_0=f^(0)-p, d_1=f^(1), d_2=f^(2), so d=f-p=d_0+d_1+d_2.

The p=q blocks contain the reciprocal stationary point u=alpha_P/r, including the two incomplete endpoint diagonals. Their defect is d_0=f^(0)-p and must retain carrier cancellation.

The p=q+1 and p=q-1 blocks have no positive uniform gradient gap and meet the diagonal cells at half-integer boundaries. They form d_1 and must be joined in opposite orientations with half-open transfer.

The |p-q|>=2 blocks form d_2 and obey |alpha_P/u-r|>=alpha_P(|p-q|-1)/{(p+1/2)(q+1/2)} before moduli. Their determinant increment still contains a cross with the updated diagonal-plus-adjacent base and a far self-determinant.

With B_0=b, B_1=b-d_0, B_2=b-d_0-d_1 and Delta_j=Q(B_j-d_j)-Q(B_j), Q(b-d)-Q(b)=Delta_0+Delta_1+Delta_2.

Using the original base b separately in all three increments omits 2d_0^TJd_1+2d_0^TJd_2+2d_1^TJd_2. The diagonal-first order assigns each cross to the later adjacent or far increment.

For I(B,D)=partial_xi{Q(B-D)-Q(B)}=-2B^TJdot(D)-2dot(B)^TJD+2D^TJdot(D), the completed geometric current is the sum I(B_0,d_0)+I(B_1,d_1)+I(B_2,d_2), followed by the fixed factor 2/|nu|^2.

## Exact Rational Audit

The naive original-base split misses `-220150159/259459200` in the stored rational witness.

These independent rational matrices audit universal identities and cross ownership; they are not asserted to be physical Xi observations.

## Tie Rule

At a reciprocal boundary, transfer the complete block mode and its lift. The determinant total is tie invariant; no floor function is differentiated.

## Next Action

Compute the physical diagonal and adjacent block matrices jointly at calibrated B=-49/100 level points, preserving source phases and endpoint halves. Bound the far invariant increment using the certified gradient gap only after its cross with the updated base is retained. Record det(B) to determine whether the nonsingular gauge is available; otherwise stratify the null-cone case. Compare the block telescope against the near/paired-remote chart before selecting an arithmetic lemma.

## Pi Provenance

No new pi is introduced by the determinant symmetry or block telescope. Pi remains inherited from e(x)=exp(2pi i x), alpha_P=xi/(2pi), and the completed-zeta saddle normalization.

## Proof Boundary

This gate proves the determinant-matrix chart, division-free adjugate identity, determinant-product-one symmetry, conditional nonsingular normal forms, singular-base guard, and exact diagonal/adjacent/far primitive and current telescopes. It proves no physical block sign or bound, no nonsingularity of P_ret, no geometric-current inequality, complete-current sign, all-q transport, contact exclusion, Lambda<=0, PF-infinity, RH, or prize-level conclusion.

built Newman C1 finite determinant reciprocal-block symmetry gate: 20 rows, 0 issues, 10 symbolic audits, 6 exact-rational audits, 3 source audits, 3 reciprocal block classes, 1 determinant symmetry, 1 singular-base guard, 1 open arithmetic row
