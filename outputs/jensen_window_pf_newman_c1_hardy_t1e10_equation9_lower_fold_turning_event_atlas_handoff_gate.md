# Fixed-selector turning-event atlas and exact mode handoff

Date: 2026-08-12

Status: exact fixed-selector event ordering and finite-sum handoff validated;
not a proof of the Airy-to-Morse approximation join or selector jump

For a fixed odd lower endpoint `C`, a mode changes turning ownership when
its joint saddle reaches `alpha=C`:

```text
alpha_m=2m+t/(pi*m)=C,
tau_m=pi*m*(C-2m).                                    (TE1)
```

The distance from the selector top is the exact square

```text
pi*C^2/8-tau_m=pi*(C-4m)^2/8.                         (TE2)
```

Here `C=159577=4q+1`, `q=39894`.  Ordering the odd squares gives

```text
m_(2k)=q-k,       m_(2k+1)=q+k+1,
tau_n=pi*C^2/8-pi[n(n+1)/2+1/8],
tau_n-tau_(n+1)=pi(n+1).                              (TE3)
```

Thus all event heights are distinct and their minimum spacing is

```text
[3.14159265358979323846264338327950288419716939937510582097494459230781640628620899755636517 +/- 1.00e-80].                     (TE4)
```

There are exactly 399 events in the fixed
selector cell.  Immediately above its lower edge their modes form the
contiguous roster `39695..40093`.
At `t=10^10`,

```text
tau_83 > t > tau_84,
transition roster =39853..39936,                       (TE5)

tau_83-t=[57.0578988367820130751563963576447796143211720386268266398370670180952568229788191858721653 +/- 5.00e-81],
t-tau_84=[206.835884064760618955705647837833462658241057508882062322058278735761321305062736697416972 +/- 5.00e-81].     (TE6)
```

The next upward event transfers mode `39936` from the
turning block to the ordinary block.  The next downward event transfers mode
`39852` in the reverse direction.  For every one of
the 399 events, direct set checks certify
the exact finite-sum identity

```text
sum_(m in O) J_m + sum_(m in T) J_m
 =sum_(m in O union {r}) J_m + sum_(m in T minus {r}) J_m. (TE7)
```

No mode is lost or counted twice.  This identity uses the same exact finite
Poisson mode integral `J_m` on both sides; it does not assert that separately
truncated ordinary-Morse and Airy approximations agree.

The odd-selector jump is a different event.  At `C -> C+2`, the old turning
roster is empty near its upper edge, while the new selector begins with the
399-mode roster `39696..40094`
at the removed lower-endpoint strip `[C,C+2]`.  It cannot be replaced by the
one-mode identity (TE7).

Proof boundary: exact event geometry and ownership conservation for exact
finite mode sums only.  No matched Airy/Morse remainder, 399-mode selector
strip reassembly, lower-interior join, complete `T_upper`, `Lambda<=0`, RH,
or prize-level conclusion is proved.
