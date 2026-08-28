# Entire Riemann-auxiliary kernel for the joined equation-(4) packet

Date: 2026-08-27

Status: exact finite-difference reduction; quantitative enclosure remains open.

Define

```text
kappa(z)=[exp(-i*pi*z)-exp(-i*pi*z^2)]/[2i sin(pi*z)]
        =exp[-i*pi(z^2+z)/2]
          sin[pi(z^2-z)/2]/sin(pi*z).
```

The numerator vanishes at every integer and l'Hopital gives

```text
kappa(k)=k-1/2,       k integer.
```

Integer translation leaves the linear quotient fixed and sends the Gaussian quotient to its `N`th csc packet.  Therefore

```text
R_N(z)=z^(-s)[kappa(z)-kappa(z-N)],

g_W(z)=z^(-s)[kappa(z-79788)-kappa(z-2561211)].
```

At `z=0`, the bracket has removable value `2561211-79788=2481423`, exactly reproducing the locally integrable finite source.  No sine pole remains.

For

```text
L_s[h]=integral_0^infinity[Fresnel] z^(-s)h(z)dz
       -integral_C z^(-s)h(z)dz,
```

the complete source is now

```text
S_W=L_s[kappa(z-79788)-kappa(z-2561211)].
```

Thus the fully reassembled target is one entire-kernel finite difference minus the already-owned Gamma and A carriers.  This is the selected exact route for the next endpoint-complete deformation.

The kernel is the trigonometric kernel in Proposition 9 of Arias de Reyna's integral representation of Riemann's auxiliary function: https://arxiv.org/abs/2407.02016.  Only the kernel identity is imported; no auxiliary-function value or bound is transferred to `J_Z`.

Pi provenance: every `pi` is inherited from the Riemann-Siegel Gaussian/sine kernel.  The translation identity uses only integer parity.

## Proof boundary

Exact entire auxiliary-kernel identity, integer removable values, finite-prefix/window finite differences, and common-functional reassembly only. No quantitative contour deformation, Mellin remainder bound, joined-packet, J_Z, or D_K enclosure, and no non-A, all-height, Lambda<=0, PF-infinity, RH, or prize-level conclusion.
