# Event-parameter contour geometry

Date: 2026-08-12

Status: exact 399-event parameter range and one common enlarged contour
certified; not a proof of a uniform full-line remainder estimate

For event index `0<=n<=398`,

```text
4m_n-C=(-1)^(n+1)(2n+1),
epsilon_n=(4m_n-C)/C,       d_n=beta epsilon_n.       (EPC1)
```

The event square and every `pi/16` height cell give

```text
t*-tau_n=beta d_n^2,
lambda(tau_n+h)=d_n^2-h/beta,
D_(n,tau_n+h)=exp(i h z/beta)D_(n,tau_n).             (EPC2)
```

The full roster has modes `39695..40093` and
`-797<=4m-C<=795`.  Its largest
absolute detuning is

```text
[10.76022909360605812680319556123182856177733999716735521367715841140865376727676740447298976053037185 +/- 4.04e-99].                     (EPC3)
```

The canonical and exact event-cell saddle envelopes are, respectively,

```text
|z_can|<=[13.40830419121808393194148531881860830411441914148959325304795030592823318294307043059252230032032503 +/- 1.64e-99],
|z_exact|<=[13.40839309443151453305515997906104857371606320117144517359263763214173126193026815790329647632081621 +/- 1.82e-94].       (EPC4)
```

Therefore the prototype `R=9` contour is not atlas-uniform: it misses the
outer saddle envelope for 246 events, beginning at
event 153.  A fixed contour with vertical
connectors at `Re z=+/-14`, lift `Im z=1`, and signed horizontal integration
through `|Re z|=21` instead has certified saddle clearances

```text
canonical: [0.5916958087819160680585146811813916958855808585104067469520496940717668170569295694074776996796749685 +/- 2.33e-100],
exact:     [0.5916069055684854669448400209389514262839367988285548264073623678582687380697318420967035236791837864 +/- 1.82e-94].   (EPC5)
```

The minimum factored phase bracket on every lifted Cauchy disk is
`[8.750124927299269612943248932349292328705992463415090357794106692063898131903646497181230031033629023 +/- 8.26e-16]>0`.  At `|Re z|=21`, the
exact logistic-Gaussian and canonical Gaussian exponent margins are
`[40.21715810181322138539558527506100353622849385127164089432000809906730084186458399718123003103362902 +/- 6.24e-99]` and
`[260.8840453824302307883375533311345744804772338745635877471916352771742762042800577369370037957347696 +/- 3.38e-98]`.  The entire contour stays
inside the local hyperbolic chart with `|q|<0.02` and a zero-free cosh
margin.

This establishes common contour geometry only.  It does not show that the
signed exact-minus-canonical integral stays within the prototype numerical
budget at extreme detuning.  It does not propagate all 399 event remainders.
It does not establish complete `T_upper` or prove `Lambda<=0` or RH.  It does
not establish a prize-level conclusion.
