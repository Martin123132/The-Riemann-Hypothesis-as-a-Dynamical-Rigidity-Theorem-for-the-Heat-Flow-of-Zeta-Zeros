# Separable triangular geometry of a signed Fourier pair

Date: 2026-08-13

Status: exact-lemma certificate; not a proof of the uniform triangle estimate

Insert the endpoint-driven solution into its Volterra kernel before taking
absolute values.  For each source endpoint `D` define

```text
phi_(m,D)(z)=pi*m^2/z+pi*D^2*z/4,
psi_m(x)=-pi*m^2/x+(t/2)log((1-x)/x).                 (TG1)
```

Then every symmetric pair has the exact separable representation

```text
K_m+K_-m=sum_(D in {B,A}) eps_D/(2i*pi)
 integral_(0<z<x<1/2)
 z^(-1/2)(2+i*pi*D^2*z)
 x^(-7/4)(1-x)^(-1/4)
 exp[i(phi_(m,D)(z)+psi_m(x))] dx dz,                 (TG2)

eps_B=+1, eps_A=-1.
```

The phase Hessian is diagonal.  Its stationary coordinates are

```text
z_D=2m/D,
x_m=2*pi*m^2/(t+2*pi*m^2),
alpha_m=2m+t/(pi*m).                                  (TG3)
```

Their signed distance to the triangular face is exactly

```text
x_m-z_D=2*pi*m^2(D-alpha_m)/[D(t+2*pi*m^2)].          (TG4)
```

Thus the endpoint characteristic equation `D=alpha_m` is simply the event
where the joint saddle crosses the face `z=x`.  The positive `z` curvature
and negative `x` curvature show that this is the same hyperbolic saddle
geometry as the original portcullis, now with `+m/-m` cancellation already
built in.

At the saved height the exact face ledger is:

```text
B gap at m=621: -0.00000021678693959597129786882090427560385690750922355708063522382819287,
B gap at m=622: 0.000000173393661448749676273577861691477636558473224622493172589709639,

A gap at m=39852: -0.000000005170891355429711968148294120587988200012735094507797265214218598,
A gap at m=39853: 0.0000000079452436065504559092989989394650126082161905296276479849484842658.    (TG5)
```

Consequently:

```text
1..621:       no B joint saddle in the triangle,
622..39852:   B saddle in, A saddle out,
39853..39894: both endpoint saddles in, with opposite signs,
m>=39895:     outer x saddle beyond the half-boundary. (TG6)
```

The ordinary/fold ownership boundary remains `39694|39695`; (TG6) is
stationary geometry and does not override that certified proof allocation.

This unifies the earlier obstructions.  The B tangent profile is the local
half-plane approximation to (TG2) at the `B` face.  The A atlas is the same
face geometry with the opposite endpoint sign.  The `39894|39895` event is
the additional corner `x=1/2`.  They should therefore be estimated with one
triangle theorem and compatible face/corner charts, rather than three
unrelated approximations.

Pi provenance: all `pi` factors come from the exact equation-(9)
Kummer/Fourier phase.  No fitted constant is used.

Proof boundary: exact separable triangular normal form, stationary
coordinates, gap identity, and saved-height face ledger only.  No uniform
triangle estimate, quantitative face/corner chart, complete paired residual,
`T_upper`, height-uniform theorem, `Lambda<=0`, PF-infinity, RH, or
prize-level conclusion is proved.
