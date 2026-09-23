# The whole low-label background can now be paid for

The previous stage removed repeated large labels and the all-small
numerator. It correctly left small labels inside mixed terms and the
denominator. This continuation pays for removing that entire low-label
field, rather than treating those interactions as free.

## Main result

For K=floor(exp(sqrt(d))), d=ceil(log(N+1)), and a suitably slower growing
grouping-prime set H_N, the original complete response satisfies

    D_N=J_N^high(H_N,K)+O_(I,S)(sqrt(m_N^max)+d^-1/8).

All five factors in the new numerator have distinct primitive labels
above K. Its two angular modes use one common high-field denominator.
Changing the original denominator is explicitly paid by a quadratic
energy bound and a global homogeneous Lipschitz estimate.

The new remainder tends to zero under the inherited first-moment theorem.
Its rate is weaker than the predecessor's: it pays for a stronger change
of field. The grouping choice depends on the unknown first moment, so
this is not an effective finite-height prescription.

An elementary frequency calculation supplies the finite low-field bound

    E_low <= q_max L_H^2/(Z_W rho_S d)
       [H_K+4K(H_K+1)/T_window] int_I a^-1 dx.

The length payment is K log K rather than the cruder K^2. It also yields
a polynomial-cutoff alternative: for fixed H and 0<delta<2/3, removing
primitive labels up to N^delta costs at most O_(I,S)(L_H^5 sqrt(delta))
in the limit. N must tend to infinity before delta tends to zero.

## The remaining obstacle

The remaining signed current is still unbounded. It contains only five
distinct high-primitive factors, without unpaid small-field feedback.
The whole high-field denominator, both angular modes and global prime
centering must still be retained.

Distinct macroscopic labels can satisfy a*b*(ce)=(ac)*(be), so independent
phase cancellation does not follow. The sufficient target remains a
limsup below 0.00092027989. The complete original budget and all other
proof-wide obligations remain unchanged. RH has not been proved.

## Checks and preservation

4,370 exact assertions passed. Two-precision point, whole-window energy
and scalar checks give 1,080 interval overlaps. The whole-window tests
cover 2,048 complete height cells per precision, but only the declared
finite N64/512 and x1/S257 controls. They do not establish the open sign.

All 8,705 historical pins were unchanged at intake. Final receipts record
the replay and final hash audit. No provider calls, new Forge experiments,
GitHub update or publication. Read DERIVATION.md for the finite payments
and AUDIT.md for the exact limitations.
