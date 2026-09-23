# Reassembly accepted; a concrete shift-cancellation route tested

The incoming Euler reassembly passes audit in its stated analytic scope.
The complete original odd response is now represented by an ordinary finite
Dirichlet polynomial, at o(log N) comparison cost and with no added fixed
charge. The same gate, windows, norm atom, and budget are retained.

## New source-dependent result

A microscopic height shift h=q/d changes the source's logarithmic symbol
by exp(iq y). Its original gated quadratic overlap has the explicit limit

    K_I(q)=int_0^1 [int_I c_x(y)^2/(4a(x)) dx] exp(iq y) dy.

This determines the whole covariance matrix of every fixed finite set of
such shifts. It is a deterministic arithmetic consequence, not a random
phase or Gaussian model.

## The attempted cancellation and its outcome

At q=2pi m, the original gate is preserved up to a vanishing error. The
averaged odd response is shift-invariant at leading order. A natural test
therefore pairs V(t) with V(t+2pi m/d), hoping the latter is nearly -V(t).

The exact limiting antipodal squared defect is

    D_I(2pi m)=2int_0^1 h_I(y)[1+cos(2pi m y)]dy.

The original profile forces D_I0(2pi m)>=6859/168000000 for every integer m.
Consequently the global sharp-stability/Cauchy certificate

    B_I0(m)=(alpha+beta) sqrt(E(I0) D_I0(2pi m))

is always greater than 0.000986957156066316. The required budget is
0.00092128989. No fixed gate-period shift can close the target USING THIS
PARTICULAR CERTIFICATE.

I then strengthened this test rather than stopping at that coarse constant.
Quadratic isotropy gives a sharp averaged gradient bound with
g^2=2(alpha+beta)^2+(alpha-3beta)^2/2. It supplies a smaller secant certificate
which retains the known covariance throughout the connecting segment.

Eight disjoint strips of the actual profile give the stronger uniform defect
floor 0.0001918186133600196. Even this ISOTROPY-AWARE secant/Cauchy certificate
is greater than 0.001330614991372663, still above the original budget.

This is a limitation of both displayed certificates, NOT evidence that the
actual signed response is too large. The new averaged gradient result is
proved using the actual alpha>5beta cone relation and second-energy isotropy.
A signed evaluation of the paired response, a stronger joint correlation
estimate, or the retained sharp-dual residual could still improve it.

## Where to work next

Retain the ordinary finite polynomial and the explicit shift kernel. The
remaining work must control the signed first/third angular response rather
than replace it by this global antipodal distance. The reassembly and quadratic
isotropy remain valid. The favourable signed inequality and RH remain open.

VERIFY_REFINEMENT_final.json records the consolidated accepted checks. The
complete analytic proof and all dependencies are in DERIVATION.md. No actual
source values or new zero counts were computed. Initial failed output export
and its corrected successor are retained. All old evidence is unchanged.
