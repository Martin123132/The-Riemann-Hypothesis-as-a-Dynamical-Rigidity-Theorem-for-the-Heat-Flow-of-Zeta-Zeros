# Joined Lower/Interior Abel Recombination Gate

Date: 2026-08-02

Status: exact physical kernel expansion and lower/interior/finite-cell recombination proved; signed finite-cell-defect estimate open; not a proof of RH.

## Physical Kernels

B_(H,N)=(R-conjugate(R_N)){conjugate(C_N)Q-conjugate(Q_N)C}+(R_x+conjugate(R_(N,x)))conjugate(C_N)C.

After Q=(R+delta)C+D, B_(H,N)=conjugate(C_N)C{(R-conjugate(R_N))^2+(R-conjugate(R_N))(delta-conjugate(delta))+R_x+conjugate(R_(N,x))}+(R-conjugate(R_N)){conjugate(C_N)D-conjugate(D_N)C}.

B_(T,N)=(R-R_N){C_NQ-Q_NC}+(R_x+R_(N,x))C_NC.

After Q=(R+delta)C+D, B_(T,N)=C_NC{(R-R_N)^2+R_x+R_(N,x)}+(R-R_N){C_ND-D_NC}.

For C=1, D=delta=0, R=ix/2, R_N=-iu_N/2, and R_x=R_(N,x)=-iu_(N,x)/2, B_(H,N)^0=-(x-u_N)^2/4 and B_(T,N)^0=-(x+u_N)^2/4-iu_(N,x).

P_H^0=-i(x+u_N)(x-u_N)^2/4 and P_T^0=-i(x-u_N)(x+u_N)^2/4+u_(N,x)(x-u_N). Both ideal lifts are cubic; the physical correction can restore degrees four and five.

At x=-u_q, P_H^0=i(u_q-u_N)(u_q+u_N)^2/4 and P_T^0=i(u_q+u_N)(u_q-u_N)^2/4-u_(N,x)(u_q+u_N).

Define Delta_H=P_H-P_H^0 and Delta_T=P_T-P_T^0 only after the full terminal conjugations above. Every term in either Delta vanishes under the simultaneous ideal specialization; no modulus or coefficient budget is taken here.

## Lower Endpoint

For P_L=i(lambda-L)P, the exact lower block is -T_1[P_L]-U_1[P_L]=i kappa L P(0)S_1+i kappa^2{[P(0)-L(P'(0)+g_0P(0))]S_2-alpha_P L P(0)S_3}.

Use L=log N in the Hermitian channel and L=log a+u_N in the transpose channel. The logarithms are retained until reverse recombination.

## Exact Recombination

For every outer mode r, -T_(1,r)[P]-U_(1,r)[P]+I_(int,r)[P]=I_P(r)-T_(N,r)[P]-U_(N,r)[P]. Thus the lower S_1/S_2/S_3 traces and the absolute interior reverse exactly to the original mode integral minus its upper boundary.

After symmetric summation, L_N[P]:=-T_1[P]-U_1[P]+I_N[P]=O_N[P]-T_N[P]-U_N[P], and R_N[P]=O_N[P]-T_N[P]. The lower logarithms are representation coordinates, not a separately small error and not a zero.

Let F_N[P]=sum_(r in mathcal R_N)I_P(r), let M_N[P] be the starred physical carrier from full finite Poisson inversion on [1,N], and define D_N[P]=F_N[P]-M_N[P]. Since F_N+O_N=M_N, one has R_N[P]=-D_N[P]-T_N[P] and F_N[P]+R_N[P]=M_N[P]-T_N[P].

The tiled finite-cell formula F_N=tau_band+M_N-E_band gives D_N=tau_band-E_band and R_N=E_band-tau_band-T_N. The finite-cell defect is the global two-sided exterior/tie package, not a new local endpoint estimate.

At a half-open roster tie, F_N gains one complete mode and O_N loses the same mode. F_N+O_N and F_N+R_N are invariant; D_N and R_N change oppositely.

## Relative Channels

For P_H=i(lambda-log N)B_(H,N), T_N[P_H]=0. Hence R_N[P_H]=-D_N[P_H] and F_N[P_H]+R_N[P_H]=M_N[P_H] exactly.

For P_T=i(lambda-log a-u_N)B_(T,N), T_N[P_T]=-2iu_NT_N[B_(T,N)]. Hence R_N[P_T]=-D_N[P_T]+2iu_NT_N[B_(T,N)] and F_N[P_T]+R_N[P_T]=M_N[P_T]+2iu_NT_N[B_(T,N)].

If the joined finite object is the returned carrier M_N rather than the band-mode sum F_N, then M_N[P_H]+R_N[P_H]=M_N[P_H]-D_N[P_H], while the transpose channel also has +2iu_NT_N[B_(T,N)]. Dropping D_N would apply Fourier inversion twice.

M_N[P_H] is the starred q<N sum with multiplier -i(u_q-u_N)B_(H,N)(log q); its q=N atom vanishes. M_N[P_T] is the starred q<=N sum with multiplier -i(u_q+u_N)B_(T,N)(log q).

## Abel Form

For a_q and k_q with A_j=sum_(q<=j)a_q, sum_(q=1)^N a_qk_q=A_Nk_N+sum_(q=1)^(N-1)A_q(k_q-k_(q+1)).

Because the Hermitian relative kernel has k_N=0, its exact Abel form has no terminal boundary: sum_(q=1)^N a_qk_q=sum_(q=1)^(N-1)A_q(k_q-k_(q+1)). This is where phase-difference cancellation can be sought without detaching the finite-cell defect.

The transpose relative kernel generally has k_N=-2iu_NB_(T,N)(log N), so its Abel terminal boundary survives and must be combined with +2iu_NT_N[B_(T,N)] and the pure terminal transpose block.

Abel summation may now be applied to M_N-D_N (and the explicit transpose terminal term), not to M_N alone while silently deleting D_N. No partial-sum estimate or sign is asserted by this gate.

## Handoff

The lower logarithmic traces are no longer an independent target. After exact reverse integration by parts, the live obstruction is the source-specific returned carrier minus the finite-cell defect, plus one explicit transpose terminal term.

Derive a source-specific formula and signed estimate for D_N[P_H] and D_N[P_T] from the global exterior/tie kernel, then Abel-sum M_N-D_N with the ideal cubic kept intact and Delta_H,Delta_T attached only afterward.

No h^2/800 reserve is claimed. The next bound must include both phase families, the correction perturbations, fixed affine row, pure terminal transpose block, and moved-tail defect before the whole-jet C1 transfer.

## Pi Provenance

The constants `kappa=1/(2pi i)` and `alpha_P=xi/(2pi)` and the character `e(x)=exp(2pi i x)` are inherited. No circle, polygon, plotted symmetry, or fitted constant is introduced.

## Proof Boundary

This gate proves exact physical Hermitian and transpose kernel factorizations, their ideal cubic relative lifts, the complete lifted lower endpoint, reverse twofold integration by parts, the distinction between reciprocal-band modes and returned carriers, reduction of the outer remainder to the finite-cell defect, both relative-channel recombinations, and exact Abel identities. It proves no signed finite-cell-defect estimate, completed mixed-current or `Phi_B` bound, contact exclusion, retained aggregate or Xi theorem, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion.
