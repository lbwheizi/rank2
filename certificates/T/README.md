# Exact verification of the T-case calculations

This directory contains the existing coefficient verifier and five additional
GAP verifiers for the revised T classification. The corresponding
arguments occur in `sec:classification-T` and `app:T-calculations` of the
manuscript. References below use stable LaTeX labels.

The revised `eq:T-classification-G2` retains its original five conditions
and adds

```text
Phi(x1*x4^-1,x1*x3^-1,x1*x2^-1) *
Phi(x1*x3^-1,x1*x2^-1,x1*x4^-1) *
Phi(x1*x2^-1,x1*x4^-1,x1*x3^-1) = 1.
```

In the manuscript, `lem:T-pullback-cocycle` identifies this equality with
triviality of the class of the pullback to T. The original five conditions
alone also allow the nontrivial class treated by the additional GAP programs.
Consequently, `epsilon_Phi(W)=1` alone does not justify replacing the
braiding by the constant `(-1)` table.

`verify_T_case.g` is a self-contained GAP-core certificate for the finite
ordered-component calculations.

It was tested with GAP 4.16.0. It loads no GAP package and uses no
floating-point arithmetic.  The coordinate formula for
`epsilon_Phi(W)` is checked by exact Laurent-monomial arithmetic; the
later adjoint calculations use exact coefficient lists in `Z[p]`.

## Run

All six verifiers require GAP 4.16.0 core only. From this directory, run

```sh
gap --bare -A -r -q --quitonbreak --norepl verify_all_T.g
```

This entry point runs `verify_T_case.g` and the five additional verifiers.
To run one verifier, replace `verify_all_T.g` by its basename. Each
verifier has a recorded output file with the same basename and suffix
`.out`. The root `verify_all.sh` runs all 23 repository suites separately
and compares their output with these files.

For an interactive run, start `gap --bare -A -r -q` in this directory,
then enter `Read("verify_all_T.g");` at `gap>`. Individual verifiers can be
loaded in the same way. Success returns to `gap>`; failure calls
`ErrorNoReturn` and stops at `brk>`, where `quit;` returns to the outer
`gap>` prompt. None of the verifiers calls `QUIT_GAP`.

In batch mode, a successful process exits with status zero. The option
`--quitonbreak` makes a failed check exit with nonzero status. The root
runner additionally requires empty standard error and byte-for-byte
agreement with the recorded standard output.

The files `T_model_core.g`, `q8_primitive8_cocycle_data.g`, and
`order2_reflection_certificate_data.g` are loaded by the verifiers that use
them. Keep these files together with the entry points. No external GAP
package or other language runtime is required. Verification does not
rewrite source files, stored data, or recorded outputs.

## Assumptions of `verify_T_case.g`

- The coordinate identity `eq:T-epsilon-explicit` assumes `dim(sigma)=1`
  and `sigma(x1)=-1`. It needs neither `epsilon_Phi(W)=1` nor the new
  parameter equality, and it does not prove either condition.
- The use of the constant `(-1)` table requires `dim(sigma)=1`,
  `sigma(x1)=-1`, `epsilon_Phi(W)=1`, and the new parameter equality. The
  simultaneous choice of tensor bases is supplied by
  `lem:T-compatible-normalization` and `lem:T-tetrahedral-braiding`.
- The higher-adjoint calculations additionally use `dim(rho)=1` and
  `p=rho(x1)*sigma(z)`. The Y2 and Y3 tables, Y4 vanishing, and the mixed
  monodromy eigenvalue are checked modulo `p^2-p+1`. The Y4
  factorization is first checked in the unreduced ring `Z[p]`.

The original verifier neither accepts an arbitrary input cocycle nor checks its
parameter conditions. It verifies the coefficient identities after the
comparison proved in the manuscript. In particular, it does not prove the
cochain equivalence or its compatibility with every tensor degree.

## What `verify_T_case.g` reconstructs

The verifier first reconstructs the tetrahedral rack table from the four
conjugation permutations

```text
Ad(x1)=(2,4,3), Ad(x2)=(1,3,4),
Ad(x3)=(1,4,2), Ad(x4)=(1,2,3).
```

For `eq:T-epsilon-explicit`, it records the three basis definitions

```text
a2=-x4 action a1, a3=-x2 action a1, a4=-x3 action a1
```

and derives the three required actions from the rack and
`eq:YD-projective-action`.  If `a_j=-x_s action a1`, the rack determines the
unique carrier `x_k` having the same action on the degree `x1` as `x_s^2`.
The program therefore constructs the centralizer tail `x_k^-1*x_s^2` and
derives the corresponding quotient of local-cocycle and
projective-representation values.  It
then performs three actual braids beginning with `a2 tensor a3` and compares
the accumulated Laurent monomial with the displayed right-hand side.  This
part assumes neither `epsilon_Phi(W)=1` nor the new parameter equality
and does not use the later constant coefficient table.

The remaining checks do not use the displayed adjoint tables as source data.
They check all 64 instances of the braid equation for the constant `(-1)`
coefficient table in the compatible bases fixed in the manuscript under the
assumptions listed above. They also check that this table is consistent with
`epsilon=1`; they do not derive the value of epsilon for an arbitrary input
cocycle or construct those bases. A central marker `0` is then appended to every
ordered word. Actual adjacent braids in the order

```text
1,2,...,m,m,...,2,1
```

reconstruct the operation h_m. The program evaluates the compatible-basis
coefficient recursion labelled `eq:T-adjoint-recursion`:

```text
phi_m(e_i1...im)
  = e_i1...im - p e_hm(i1,...,im)
    - Prefix(i1 ▷ i2, phi_(m-1)(e_i1,i3,...,im)),

phi_1(e_i) = (1-p)e_i.
```

Polynomials in p are represented directly by integer coefficient lists.
Thus every comparison is exact.

## Identities checked by `verify_T_case.g`

- **`eq:T-epsilon-explicit`:** the three action coefficients are generated
  from the tetrahedral rack, the basis definitions, and
  `eq:YD-projective-action`.  Three successive braidings are then
  accumulated.  The resulting exact Laurent monomial is compared with the
  displayed formula, including its total minus sign, its three local-cocycle
  quotients, and its three `sigma` factors.  The three displayed action
  coefficients are not used as definitions of the generated coefficients.
- **`lem:T-Y2-homogeneous-components` and
  `eq:T-Y2-basis`:** the four displayed basis vectors and all 16 entries
  of the internal complete phi_2 calculation, reduced only modulo
  f(p)=p^2-p+1.
- **`lem:T-Y3-homogeneous-components` and
  `eq:T-y3-generator`:** the displayed generator and all 16 entries of
  the internal complete phi_3 calculation, reduced only modulo f(p).
- **`lem:T-Y4-homogeneous-components`:** the internal exact
  factorization `phi_4(w_1 tensor y)=f(p)Z_1` in the unreduced ring Z[p],
  together with the assertion that the 15 ordered words of the certificate
  target Z_1 are distinct. Only after this exact factorization is checked
  does the verifier reduce modulo f(p).
- **`eq:T-R2-monodromy-cycles` and `eq:T-R2-mixed-monodromy`:** the four
  three-cycles are derived by applying the eight
  adjacent braids to `[1,i,j,k,0]`; they are not used as the definition of
  the monodromy. The program counts six W-W crossings and one r*ell=p
  mixed monodromy, recovers the four displayed coefficient triples of y,
  and checks `M(y)=(1-p)y` modulo f(p).

The basis and generator in `eq:T-Y2-basis` and `eq:T-y3-generator`, and
the coefficient patterns used in the three homogeneous-component lemmas
listed above, are comparison targets.
The complete phi_2 and phi_3 calculations and Z_1 remain internal certificate
data rather than separately labelled formulas in the manuscript. All
left-hand sides are generated afresh from the rack and
`eq:T-adjoint-recursion`.

The original verifier retains its 15 exact checks.

## Additional GAP calculations

The nontrivial class has value `-1` for the new three-factor product, while
the original five conditions hold. The finite model is
`G_*=SL_2(3) x C_6`, with the cocycle and pair in
`eq:T-order-two-finite-cocycle` and `eq:T-order-two-finite-pair`. Its
third adjoint has dimension two. The files below verify the fixed cocycle
and the finite calculations used in the proof that this class gives an
infinite-dimensional rank-two Nichols algebra.

| File | Exact calculation | Corresponding manuscript argument |
| --- | --- | --- |
| `verify_T_cohomology_reduction.g` | The 512 quaternion cocycle entries; all `8^4` cocycle identities; normalization, cyclic invariance and parity; the bar-cycle boundary and evaluation; induced-action exponents; the two possibilities under `q=-1` and `epsilon=1`. | `lem:T-fixed-quaternion-cocycle`, the cocycle evaluation in `lem:T-pullback-cocycle`, and the coefficient calculation in `lem:T-cocycle-two-branches`. |
| `verify_T_sl2_cocycle.g` | The group law on `SL_2(3)`, all sign-cocycle identities, the four supports and induced projective action, diagonal braiding and epsilon, braid relations, and the nontrivial full twist. | The finite model in `eq:T-order-two-finite-cocycle` and `eq:T-order-two-finite-pair`, used by `lem:T-nontrivial-cocycle-infinite`. The cocycle is independent of the central `C_6` coordinate. |
| `verify_T_actual_adjoint.g` | Restricted recursive adjoint dimensions `4,4,2,0`, a nonzero two-coordinate minor for the third adjoint, and direct quantum-symmetrizer dimensions in degrees 2, 3 and 4. | The initial row of `eq:T-order-two-adjoint-dimensions` and the two-dimensional object in `lem:T-order-two-reflection-calculation`. |
| `verify_T_nichols_iterative.g` | Iterated quantum-symmetrizer image dimensions `1,4,8,10,8,4,1,0`, with sum 36. | An auxiliary calculation of `B(W_0)` for the nontrivial finite model. This is not the dimension of `B(V_0 direct_sum W_0)` and is not the ordinary tetrahedral dimension 72. |
| `verify_T_order2_reflections.g` | Four successive recursive adjoint and dual constructions, all five rows of the support/action and adjoint-dimension tables, projective actions, central powers, mixed braid relations, invariant dual evaluation pairings, and the final comparison of all `SL_2(3)` actions. | `lem:T-order-two-reflection-calculation`, `eq:T-order-two-reflected-pairs`, `eq:T-order-two-adjoint-dimensions`, and the finite action comparison used in `lem:T-nontrivial-cocycle-infinite`. |

The complete 512-entry quaternion table is stored in
`q8_primitive8_cocycle_data.g`. The file
`order2_reflection_certificate_data.g` stores all 32 rank records, including
the exact sparse input columns, image columns and independent bases, and
all 27 action records. These are comparison targets: the reflection
verifier reconstructs the corresponding columns and actions and compares
the complete records, without overwriting the data file.

The shared `T_model_core.g` supplies the finite-group operations,
associator and projective-action coefficients, adjacent braids, sparse
linear algebra, and exact coefficient arithmetic used by the finite-model
verifiers. For calculations involving p, GAP uses `p=E(6)`. Since
`p^2-p+1` is irreducible over Q, evaluation at `E(6)` identifies
`Q[p]/(p^2-p+1)` with `Q(E(6))`. Thus the computations are exact and apply
to both roots of the polynomial by the two field embeddings; choosing
`E(6)` is not a numerical specialization test.

These programs verify the stated finite calculations. The homology
computation for T, the reduction of every input cocycle to the fixed
representatives, the cochain equivalence in all tensor degrees, the
existence of reflections under the finite-dimensionality assumption, and
the infinite sequence of real roots remain proofs in the manuscript. In
particular, checking four reflected pairs does not by itself prove that
all subsequent Cartan matrices repeat. That conclusion also uses the
cochain comparison established in `lem:T-nontrivial-cocycle-infinite`.

## Scope of `verify_T_case.g`

This certificate verifies only the finite coefficient calculations just
listed. Its first part verifies the coordinate identity
`eq:T-epsilon-explicit`; it does **not** prove that the resulting scalar is
one. The second part uses the hypotheses of
`lem:T-tetrahedral-braiding`, including the new parameter equality, when
it uses the constant `(-1)` braiding table.  The certificate does not derive the
compatible bases, the action of `x_1^{-1}` on
`(W^*)_{x_1^{-1}}`, or the rigid-dual identity `epsilon_Phi(W^*)=1`.
It also does not prove the normal-form lemma for
arbitrary b_ij, simplicity of the adjoint objects, existence or
involutivity of reflections, standardness of the Cartan graph, the tensor
decomposition indexed by the positive roots, or finite-dimensionality of
a Nichols algebra. Those are mathematical arguments in the manuscript. For
`lem:T-reflection-cocycle-reductions`, the program checks only
`eq:T-R2-monodromy-cycles` and the coefficient calculation in
`eq:T-R2-mixed-monodromy`. It does not check the first-reflection formulas,
the full-twist scalar, `eq:T-R2-central-coefficient`,
`eq:T-R2-cross-coefficient`, or the rigid-dual calculation of `s_2` and
`epsilon_2`. These parts are proved in the manuscript.
