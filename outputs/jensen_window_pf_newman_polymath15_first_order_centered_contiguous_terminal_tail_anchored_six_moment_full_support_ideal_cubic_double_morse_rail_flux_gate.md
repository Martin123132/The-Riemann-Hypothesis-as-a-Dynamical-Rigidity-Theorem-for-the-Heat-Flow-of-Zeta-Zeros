# Ideal-Cubic Double-Morse Rail-Flux Gate

Date: 2026-08-03

Status: exact double-Morse geometry, curvilinear boundary flux, moving-endpoint shear, tie transport, and ideal rail values proved; signed complete ideal estimate open; not a proof of RH.

## Double-Morse Map

Put rho=pr/alpha_P and v=ur/alpha_P. With eta=sqrt(alpha_P)z(rho) and y=sqrt(alpha_P)z(v), the inverse map is r=alpha_Prho/p and u=pv/rho.

For Psi_p(r,u)=alpha_Plog u-r(u-p), the exact identity is Psi_p=alpha_Plog p+(eta^2-y^2)/2. The sequential positive and negative Morse phases are therefore one two-dimensional saddle chart.

The exact Jacobian is dr du={J(rho)J(v)/rho}deta dy. Hence the transformed amplitude is G_P(eta,y)=exp(S(log(pv/rho)))P(log(pv/rho))J(rho)J(v)/rho.

The physical line u=p is exactly v=rho, hence y=eta and Psi_p=alpha_Plog p. It is the returned-carrier ridge, not an oscillatory error rail.

The reciprocal cell endpoints rho_-=p/(p+1/2) and rho_+=p/(p-1/2) map to fixed eta endpoints. Their inner saddles are u_0=p+1/2 and u_0=p-1/2 respectively.

## Curvilinear Collar

For an internal radius-L collar, U_-=p-L-1/2 and U_+=p+L+1/2. Its Morse image is eta_-<=eta<=eta_+ and y_-(eta)<=y<=y_+(eta), where y_U(eta)=sqrt(alpha_P)z(Urho(eta)/p). This is curvilinear, not Cartesian.

Along a physical rail u=U, y_U'(eta)=(U/p)J(rho)/J(Urho/p)=vJ(rho)/{rho J(v)}>0. The endpoint velocity vanishes nowhere; only its product with eta vanishes at the reciprocal saddle.

For L=1, the cells J_(p-1), J_p, and J_(p+1) tile [p-3/2,p+3/2]. Their two internal Morse interfaces cancel with opposite orientations, leaving only the two outer collar rails.

## Finite Correction Flux

For E(eta,y)=e((eta^2-y^2)/2), div(kappa eta E,-kappa y E)=(eta^2+y^2)E. The full-line first-correction cancellation is the zero-boundary limit of this exact divergence law.

On the true collar image, the paired finite second moments equal kappa times the four-side flux integral int_(boundary Omega)E{eta dy+y deta}. The reciprocal sides contribute kappa eta E dy; a physical rail contributes kappa E{y_U+eta y_U'}deta with its boundary orientation.

The term kappa eta y_U'E is the moving-endpoint shear. It is absent only in an artificial Cartesian Morse rectangle. Thus the four fixed-limit terms in the two one-dimensional moment formulas do not by themselves equal the physical collar rails.

If F(eta)=integral_(y_-(eta))^(y_+(eta))e(-y^2/2)dy, positive-moment integration by parts differentiates F and creates y_+'e(-y_+^2/2)-y_-'e(-y_-^2/2). These are exactly the two shear terms required by the boundary flux.

For a nonconstant coefficient D, div(Dkappa eta E,-Dkappa y E)=D(eta^2+y^2)E+kappa E(eta D_eta-y D_y). The physical ideal-cubic amplitude therefore leaves an interior transport term as well as boundary flux; constant-symbol cancellation cannot be promoted coefficientwise.

## Tie Transport

At the shared reciprocal boundary, Psi_(p+1)-Psi_p=r. When the boundary is an actual integer tie, e(r)=1, so overlapping complete-mode traces use the same Fourier phase.

The radius-one collar changes from [p-3/2,p+3/2] to [p-1/2,p+5/2]. After common-interface cancellation, the only transferred physical strips are the leaving J_(p-1) and entering J_(p+2), exactly as Section 11.189 requires.

The collar-side flux is therefore tie invariant only after the corresponding far-side flux is included. The moving-endpoint shear travels with the same complete physical strips and cannot be assigned to the collar alone.

## Ideal Cubics

For any positive physical rail U, put u_U=log(a/U). Then P_H^0(log U)=i(u_U-u_N)(u_U+u_N)^2/4 and P_T^0(log U)=i(u_U+u_N)(u_U-u_N)^2/4-u_(N,x)(u_U+u_N). These apply at half-integer collar rails as well as integer atoms.

The reciprocal side saddles occur at U=p+1/2 and U=p-1/2, whereas the radius-one physical rails are U=p+3/2 and U=p-3/2. Adjacent-corner transport bridges a full cell; the positive endpoint trace cannot be identified pointwise with the outer collar rail.

At U=N, P_H^0(log N)=0 but P_T^0(log N)=-2u_Nu_(N,x). The exact Hermitian physical endpoint value vanishes; the transpose endpoint remains joined to the conditional terminal quadrature and full complex recurrence.

The identity P_T^0(x)=-P_H^0(-x)+u_(N,x)(x-u_N) does not reverse the curvilinear collar: S(log u), the interval [1,N], reciprocal side traces, and endpoint velocities are not reflection symmetric.

## Handoff

Replace the proposed four fixed-endpoint match by the exact curvilinear boundary flux. Keep reciprocal side integrals, physical outer-rail integrals, and moving-endpoint shear together; internal diagonal/adjacent interfaces cancel only after their oriented union.

Insert the complete transformed ideal amplitudes G_H^0 and G_T^0 into the weighted flux identity, combine the two physical rail traces with the rail-compressed far integration-by-parts functional, and compose both reciprocal side traces with finite-cell tie corrections and the q=N transpose terminal quadrature. Then estimate or obstruct the remaining interior transport eta G_eta-y G_y without blockwise moduli.

## Pi Provenance

The character is e(x)=exp(2pi i x), so kappa=1/(2pi i). The divergence law follows by differentiating e((eta^2-y^2)/2); no geometric constant, circle, polygon, or fitted parameter is inserted.

## Proof Boundary

No signed complete ideal join, weighted rail flux, interior transport, far aggregate, h^2 reserve, completed current, or RH-level theorem is proved. This gate proves no signed physical join, completed-current bound, `Phi_B` theorem, contact exclusion, retained aggregate or Xi theorem, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion.
