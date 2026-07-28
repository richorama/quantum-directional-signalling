# Quantum Directional Signalling: Research Notes

## Status

Extracted from the `approximate-quantum-coarse-graining` branch of
`richorama/cho` at commit `a81742b`. This project starts from a literature gap
and frozen theorems, not from measured constants.

## Question

For a microscopic unitary channel and a CPTP coarse-graining `B`, quantify the
best autonomous effective channel:

```text
delta(U,B) = inf_E || B Ad_U - E B ||_diamond.
```

The objective is a rigorous relation between autonomy failure, interaction and
recoverable hidden correlations. No Standard Model or cosmological
interpretation is in scope.

## Literature boundary

- Operator entanglement of unitaries: Zanardi, *Phys. Rev. A* 63, 012301
  (2001), arXiv:quant-ph/0010074.
- Causal and semicausal operations: Beckman et al., *Phys. Rev. A* 64, 052309
  (2001), arXiv:quant-ph/0102043.
- Information-disturbance and Stinespring continuity: Kretschmann,
  Schlingemann and Werner, *IEEE Trans. Inf. Theory* 54, 1708 (2008),
  arXiv:quant-ph/0605009.
- Approximate operator-algebra QEC: Beny and Oreshkov, *Phys. Rev. Lett.* 104,
  120501 (2010), arXiv:0907.5391.
- Directional signalling and causal influence: Barsse, Perinotti, Tosini and
  Vaglini, *Phys. Rev. Research* 6, 043305 (2024), arXiv:2309.07771.
- Exact SWAP and CNOT signalling, including parallel and asymptotic uses:
  Barsse, Perinotti, Tosini and Vaglini, *Disentangling signalling and causal
  influence* (2025), arXiv:2505.14120.
- Tensor-product-structure distance: Andreadakis and Zanardi, *Quantum* 9,
  1668 (2025), arXiv:2410.02911.
- Exact quantum model reduction: Grigoletto et al., *Quantum* 9, 1814 (2025),
  arXiv:2412.05102.

The forward perturbative estimate
`delta <= 2 t ||H_interaction||` follows from channel contractivity and is not a
new result. Petz recovery and Stinespring continuity do not by themselves give
the sharp channel-intertwining defect considered here.

The quantity `delta_A(U)` is not a new measure: under the diamond norm and after
exchanging the subsystem labels, it is exactly the directional signalling
measure `S(U)` of Barsse et al., including its optimization over arbitrary correlated and
ancilla-extended inputs. The 2024 paper proves `S(CNOT)<=1`; the 2025 follow-up
proves exact CNOT and SWAP values. AQC1 independently reproduces the CNOT value
and extends it to a continuous locally equivalent Ising family. The novelty
claim is therefore the exact family evaluation, not the definition or endpoint.

## Frozen results

### AQC0 -- normalization audit

Let `T = Tr_B`, `M = Tr_B Ad_U`, and let

```text
D_A(U) = min_E_linear ||M - E T||_F^2.
```

For every bipartite unitary on dimensions `d_A x d_B`,

```text
D_A(U) = d_A^2 d_B [1 - Tr(rho_A(U)^2)],
```

where `rho_A(U)` is the `AA'` reduction of the normalized vectorized operator
`|U>/sqrt(d_A d_B)`. Thus the old raw closure defect is exactly a dimension
factor times linear operator entanglement. Its contraction across changing
dimensions cannot by itself establish RG irrelevance.

Proof: `T T^dagger = d_B I`, so least-squares projection gives
`D=||M||_F^2-||M T^dagger||_F^2/d_B`. Unitary conjugation gives
`||M||_F^2=d_A^2 d_B`; realignment of `U` gives
`||M T^dagger||_F^2=d_A^2 d_B^2 Tr(rho_A^2)`.

### AQC1 -- sharp Ising autonomy defect

For two qubits and

```text
U_theta = exp(-i theta Z x Z) = c I - i s Z x Z,
```

the operational autonomy defect under `B=Tr_B` is

```text
delta_A(U_theta) = |sin(2 theta)| = 2 |c s|.
```

An optimal effective channel is the dephasing channel
`E_theta(X)=c^2 X+s^2 Z X Z`. The residual factors as

```text
Tr_B Ad_U - E_theta Tr_B
  = sin(2 theta) L,
L(X) = -(i/2)[Z, Tr_B(Z_B X)],
```

and `||L||_diamond=1`. The upper bound follows because
`X -> Tr_B(Z_B X)` and `K -> -(i/2)[Z,K]` have completely bounded trace norm at
most one. For the lower bound, take the two physical product inputs
`|+><+|_A x |0><0|_B` and `|+><+|_A x |1><1|_B`. They have the same A marginal,
so every candidate effective channel gives them the same simulated output.
Their true A outputs differ by `2cs Y`, with trace norm `4|cs|`; the triangle
inequality forces one simulation error to be at least `2|cs|`. This matches the
upper bound without a reference ancilla. The exact code additionally certifies
the equivalent signed-operator witness and its flagged-state realization.

At `theta=pi/4`, this unitary is locally equivalent to CZ and CNOT. Hence AQC1
also gives `S(CNOT)=1` in the convention of Barsse et al. Their 2025 follow-up
proves that endpoint independently, so the endpoint is prior art; the
continuous Ising formula remains the candidate new result.

### AQC2 -- SWAP boundary control

For equal dimensions `d` and the SWAP unitary,

```text
delta_A(SWAP) = 2 (1 - 1/d^2).
```

In particular, two-qubit SWAP has defect `3/2`. The optimal effective channel
is completely depolarizing. Indeed, with that choice the residual is
`(id - Depol) Tr_A`, so channel contractivity and the standard covariant-channel
identity

```text
||id_d - Depol_d||_diamond = 2 (1 - 1/d^2)
```

give the upper bound. Conversely, fixing the A input reduces every candidate
effective channel to a replacement state on the relabelled B output. Unitary
twirling shows that the maximally mixed replacement minimizes its distance
from the identity channel. A maximally entangled B-reference input attains the
same value. This is a known channel-discrimination identity applied as a
boundary control, not a new theorem. Barsse et al. (2025) also derive this exact
SWAP signalling value in arbitrary equal dimension.

### AQC3 -- Cartan symmetry reduction

For

```text
U = exp[-i (alpha XX + beta YY + gamma ZZ)]
```

write `U=sum_mu a_mu sigma_mu x sigma_mu`. Coupling B to the maximally mixed
state gives the admissible Pauli channel

```text
F(X) = sum_mu |a_mu|^2 sigma_mu X sigma_mu.
```

The nonnegative weights sum to one, so `F` is CPTP, and
`Tr_B Ad_U(X x I/2)=F(X)` exactly. This does not yet prove that `F` is the
globally optimal effective channel.

The Cartan unitary commutes with every joint Pauli `P x P`. Consequently both
`N=Tr_B Ad_U` and `T=Tr_B` obey the same covariance,

```text
N Ad_(P x P) = Ad_P N,    T Ad_(P x P) = Ad_P T.
```

For any candidate channel `E`, conjugating `E` by `Ad_P` preserves the diamond
defect. Convexity then shows that its Pauli twirl cannot increase the defect.
Therefore at least one optimal `E` is Pauli diagonal. This is a symmetry
reduction from all qubit CPTP maps to the Pauli-channel tetrahedron, not a
solution of the remaining diamond-norm minimization and not a novelty claim by
itself.

### AQC4 -- partial-SWAP scalar theorem

On the equal-angle Cartan line `alpha=beta=gamma=t`, the identity
`XX+YY+ZZ=2 SWAP-I` makes the channel, up to global phase,

```text
U_phi = cos(phi) I - i sin(phi) SWAP,   phi=2t.
```

Its full `V x V` covariance strengthens AQC3: an optimal effective qubit
channel can be chosen depolarizing,

```text
E_lambda(X)=lambda X+(1-lambda)Tr(X)I/2,  -1/3<=lambda<=1.
```

For fixed `lambda`, put

```text
a=(1-lambda)/3,
B=(1+lambda)^2-4 lambda cos^2(phi).
```

The diamond SDP evaluates exactly to

```text
d_phi(lambda) =
  4a                         if B <= 4a^2,
  B/(sqrt(B)-a)              if B > 4a^2.
```

Proof: average the full linear diamond-SDP feasible tuple
`(rho_0,rho_1,X)`, not the nonlinear trace-norm formula, over `V x V`.
Linearity preserves its objective and produces invariant input densities.
Schur-Weyl duality then gives `rho=x P_-+y P_+`, with `x+3y=1`. In the
singlet/triplet basis, the scaled Choi matrix splits into two equivalent
blocks. Four eigenvalues are `(lambda-1)y/2`; the remaining four are two copies
of

```text
[(1-lambda)y +- sqrt((1-lambda)^2 y^2+3xyB)]/2.
```

Maximizing their trace norm over `x` gives the displayed branches. Thus

```text
delta_A(U_phi)=min_{-1/3<=lambda<=1} d_phi(lambda).
```

This scalar minimization has a closed weak-coupling branch:

```text
0 <= sin(phi) <= 1/3:
lambda_*=1,   delta_A(U_phi)=2 sin(phi).
```

For `1/3<sin(phi)<1`, a minimizer lies in `0<lambda_*<1` and satisfies

```text
3[lambda_*+1-2cos^2(phi)] [sqrt(B_*)-2(1-lambda_*)/3] = B_*.
```

It reaches `lambda_*=0` and `delta=3/2` at SWAP. A dedicated search covering
partial-SWAP collision models, covariant channel discrimination, and both
Barsse et al. signalling papers found no treatment of the continuous family or
this two-branch formula. AQC4 therefore appears novel, while its covariance
tools and exact SWAP endpoint are prior art.

### AQC5 -- arbitrary-dimensional partial SWAP

For `U_phi=cos(phi)I-i sin(phi)SWAP` on `C^d x C^d`, full `U(d)` covariance
again reduces the effective channel to

```text
E_lambda(X)=lambda X+(1-lambda)Tr(X)I/d,
-1/(d^2-1) <= lambda <= 1.
```

Let `A=1-lambda`,

```text
B=A^2+d^2 lambda sin^2(phi),
H0=d(d^2-3)/(d^2-1),
H1=2/(d^2-1).
```

For an invariant diamond witness, let `u` be its total weight on the symmetric
subspace.  Write

```text
r_- = d(d-1)/2,  r_+ = d(d+1)/2,
rho = x P_- + y P_+,  r_- x+r_+ y=1,
u=r_+ y.
```

The scaled Choi operator has:

```text
-(1-lambda)y/d, multiplicity d(d-1)(d+2)/2,
-(1-lambda)x/d, multiplicity d(d+1)(d-2)/2,
nu_+, nu_-, each with multiplicity d.
```

Putting `P=(d-1)(d+2)/(d+1)`,
`R=(d+1)(d-2)/(d-1)`, the last two roots obey

```text
nu_+ + nu_- = (1-lambda)[P u+R(1-u)]/d^2,
nu_+ nu_- =
  u(1-u){[(1-lambda)^2/d^2] P R
         -[(1+lambda)^2-4lambda cos^2(phi)]}/d^2.
```

Their signs are opposite.  Summing all eigenvalue magnitudes and setting
`t=2u-1` gives

```text
d_{d,phi}(lambda)
 = max_{-1<=t<=1} (1/d) [
     A(H0+H1 t)
     + sqrt(A^2(H0+H1 t)^2+4B(1-t^2))
   ],
```

where `t=2u-1`.  Put

```text
Q=A^2(H0^2-H1^2)+4B,
C=1-A^2 H1^2/(4B),
D=(A H0+sqrt(Q))/C.
```

If `A H1 D<=4B`, the fixed-channel norm is `D/d`; otherwise its maximum is the
endpoint value `2A(H0+H1)/d`.  This recovers AQC4 when `d=2`.

The outer optimization is convex in `lambda`.  Its left derivative at
`lambda=1` is

```text
sin(phi) - (d^2-3)/(d^2-1).
```

Therefore the exact weak branch is

```text
0 <= sin(phi) <= (d^2-3)/(d^2-1):
lambda_*=1,   delta_A(U_phi)=2 sin(phi).
```

At full SWAP, `lambda_*=0` and
`delta_A=2(1-1/d^2)`, reproducing the known endpoint.  Exact qutrit
certificates give threshold `3/4`, weak value `6/5` at
`(cos phi,sin phi)=(4/5,3/5)`, and SWAP value `16/9`.
For `cos(phi)>0`, the derivative at `lambda=0` is
`-2(d^2-1)cos^2(phi)/d^2<0`; convexity therefore excludes negative
`lambda` from the outer optimum.  Dimension-four certificates give threshold
`13/15`, weak value `8/5` at `(3/5,4/5)`, and SWAP value `15/8`.

### AQC6 -- the XY/iSWAP line

For

```text
U_theta=exp[-i theta(XX+YY)],  0<=theta<=pi/4,
C=cos(2 theta),  S=sin(2 theta),
```

the excitation-zero and excitation-two states are fixed, while the
one-excitation sector rotates by `C I-i S X`.  Joint-Pauli and axial `U(1)`
covariance reduce an optimal effective channel to

```text
E(rho)=q_I rho+q(X rho X+Y rho Y)+q_Z Z rho Z.
```

Equivalently, its Bloch eigenvalues are `(u,u,w)`, with

```text
q_I=(1+2u+w)/4, q=(1-w)/4, q_Z=(1+w-2u)/4,
-1<=w<=1, |u|<=(1+w)/2.
```

For a fixed `(u,w)`, average the complete Watrous-SDP feasible tuple over joint
axial rotations and `X tensor X`.  The remaining odd-sector coherence changes
sign under the antiunitary covariance
`X -> (Z tensor I) X* (Z tensor I)` and may also be averaged away.  Thus an
optimal input density is

```text
rho_A=A(|00><00|+|11><11|)
      +B(|01><01|+|10><10|),  B=1/2-A,  0<=A<=1/2.
```

The scaled residual Choi operator is two unitarily equivalent four-dimensional
blocks.  Each block is the direct spectral contribution `-alpha` together with
the following `M3`.  Put

```text
alpha=A(1-w)/2, gamma^2=AB/2,
g=C-u, h=S,
beta=B[(w-(C^2-S^2))+2iSC]/2.
```

The nontrivial Hermitian block is

```text
M3 = [[alpha, gamma(g+ih),  gamma(-g+ih)],
      [gamma(g-ih), 0,      beta],
      [gamma(-g-ih), beta*, 0]].
```

Consequently the exact fixed-channel norm is

```text
max_{0<=A<=1/2} 2[alpha+||M3||_1].
```

On the closed branches below, `M3` has one positive and two nonpositive
eigenvalues, so this simplifies to

```text
max_{0<=A<=1/2} 4 lambda_max(M3).
```

This gives three regimes.

**I. Weak coupling.**  Let `theta_1=0.258270520262...` be the first positive
root of

```text
2SC=3C^3+C^2-C-1.
```

For `0<=theta<=theta_1`,

```text
q=(C^2,S^2/2,S^2/2,0),
delta_A(U_theta)=2SC=sin(4 theta).
```

The invariant one-excitation density gives the matching universal lower bound.
For the displayed channel, positivity of
`H_1(A)=(SC/2)I-M3(A)` for every `A` gives the upper bound.  The leading
minor is `S(C-2AS)/2`, and the `(2,3)` principal minor is
`A C^2(A-1)(C-1)(C+1)`.  The remaining two-by-two minor is positive because
its residual quadratic has negative leading coefficient and negative
discriminant throughout the weak interval.  Finally

```text
det H_1=A C(C-1)K(A)/4,
K''(A)/2=4C(C-1)^2(C+1)>0,
K(0)=S[2SC-(3C^3+C^2-C-1)],
K(1/2)=-C(C+1)S(C-S).
```

Convexity of `K` and the signs of its endpoints prove `det H_1>=0`.  The first
threshold is exactly where `K(0)` vanishes and the boundary witness at `A=0`
ceases to maximize the norm.

**II. Intermediate coupling.**  For
`theta_1<theta<theta_2`, full SDPs consistently place the optimum on the CPTP
facet `q_Z=0`, or `u=(1+w)/2`, with an interior witness.  This active-facet
selection is a numerical conjecture, not yet a theorem.  The rigorous exact
characterization remains the two-parameter saddle

```text
delta_A(U_theta)
 = min_{(u,w) CPTP} max_{0<=A<=1/2} 2[alpha+||M3||_1].
```

Conditional on the observed `q_Z=0` facet, it reduces to

```text
delta_A(U_theta)
 = min_w max_{0<=A<=1/2} 4 lambda_max(M3)
   subject to u=(1+w)/2.
```

Eliminating `A` and `w` from those stationary equations gives a compact
conditional algebraic answer.  With `t=tan(theta)` and `d=delta_A(U_theta)`,
the physical sheet obeys

```text
Q_+(d,t) =
 d^4(t^6+t^5+3t^4+2t^3+3t^2+t+1)
+d^3(2t^7-46t^6-122t^5-170t^4-138t^3-58t^2-14t+2)
+d^2(-84t^7+24t^6+72t^5+96t^4-68t^3-56t^2-32t)
+d(96t^7+64t^6+160t^5+384t^4+320t^3+128t^2)
+64t^5(2t^2+2t+1) = 0.
```

The degree-eight elimination over `Q(C)` is the field norm
`Q_+(d,t)Q_+(d,-t)`.  Exact Sturm counts show that `Q_+` has exactly one root
`d in (1/2,2)` throughout the intermediate interval, with no collisions or
bracket crossings.  It joins the weak and strong values exactly: the first
threshold is the relevant root of
`t^6-2t^5-3t^4+7t^2+2t-1`, and the second is the root of
`t^3+t^2-1`.

This quartic is rigorous **conditional on the active facet**.  An attempted
KKT/elimination proof found exact inertia and all-witness PSD identities, but
an independent audit found that the load-bearing multiplier sign and physical
sheet tracking still use floating-point branch selection.  Therefore no
analytic active-facet proof, and no unconditional middle formula, is claimed.

**III. Central iSWAP regime.**  Let `theta_2=0.646615513406...` be the root

```text
1+C^2=2C(1+S).
```

Equivalently, `tan(theta_2)` is the positive root of `t^3+t^2-1=0`.  For
`theta_2<=theta<=pi/4`,

```text
u=C(1+S), w=C^2,
q_I=[1+C^2+2C(1+S)]/4,
q_X=q_Y=S^2/4,
q_Z=[1+C^2-2C(1+S)]/4,
delta_A(U_theta)=S+S^2/2.
```

The maximally mixed input density is the matching witness.  At that witness
the largest `M3` eigenvalue is `S(2+S)/8`, with eigenvector proportional to
`(i sqrt(2),1,1)`; its derivatives in both channel parameters vanish, proving
the lower bound by convexity.  For the displayed channel, put
`H_3(A)=[S(2+S)/8]I-M3(A)`.  Its determinant factors as

```text
det H_3 =
 -S^3(4A-1)^2(S+2)
  [A(28S^2+8S-32)+12-13S^2-4S]/512.
```

The bracket is negative at both `A=0` and `A=1/2`.  The two nontrivial
two-by-two principal minors are `S^2G_1/64` and `S^2G_2/64`: `G_1` is a
convex quadratic whose completed-square minimum is positive, while `G_2` is
concave with endpoint values `13S^2+4S-12` and `(S+2)^2`.  Together with
`H_{3,00}=S(S+2-4AS)/8`, all principal minors are nonnegative throughout
regime III.  Hence `H_3(A)>=0` for every `A`, proving the upper bound.  At
iSWAP, `S=1`, the effective channel is completely depolarizing and
`delta_A=3/2`.

The full curve obeys `delta_A(theta)=delta_A(pi/2-theta)` and returns to zero
at `theta=pi/2`, where the unitary is the local gate `Z tensor Z`.  A focused
prior-art search found only exact CNOT and SWAP signalling values, not an
evaluation of the XY line or the iSWAP endpoint.  AQC6 therefore appears
novel, subject to external literature review.  Exact rational controls cover
the weak point `(C,S)=(12/13,5/13)`, the central point `(9/41,40/41)`, and
iSWAP.

## Promotion and kill rules

1. Exact identities must be proved analytically; code supplies finite
   certificates and regression protection.
2. Frobenius norms may diagnose algebra but cannot be advertised as operational.
3. A diamond-norm claim must include a valid upper bound and a matching witness.
4. Known contractivity, continuity or recovery results must be cited rather than
   renamed.
5. Stop if the next result reduces to operator entanglement, a norm
   normalization, or a direct corollary of Stinespring continuity.
6. Prior signalling measures must be named as such; "autonomy" is only the
   coarse-graining interpretation of the same optimization.

## Next theorem

The Ising axis, equal-angle partial-SWAP line, and XY/iSWAP line now have sharp
results.  Do not return directly to the full three-parameter Cartan chamber.
The next tractable target is the two-parameter XXZ surface
`exp[-i(theta(XX+YY)+gamma ZZ)]`, using AQC6's axial covariance and active-facet
phase diagram.  Promote a result only if the additional `gamma` parameter
admits sharp witnesses or a rigorous reduction beyond a generic numerical SDP.
