# Exact verifier for the T-case coefficient calculations

`verify_T_case.g` is a self-contained GAP-core certificate for the finite
ordered-component calculations used in the T case.  The corresponding
formulas occur in Section 6 and Appendix D of the accompanying manuscript.
Stable LaTeX labels are used below so that this record remains valid if the
display numbers change.

It was tested with GAP 4.16.0. It loads no GAP package and uses no
floating-point arithmetic.

## Run

From this directory, run

```sh
gap --bare -A -r -q --quitonbreak --norepl verify_T_case.g
```

For an interactive run, first start `gap --bare -A -r -q`, then load the
script at `gap>` with `Read("verify_T_case.g");`.  Success returns to `gap>`;
failure
stops at `brk>`, where `quit;` returns to the outer `gap>` prompt.

In the batch form, the process exits with status zero after printing the
contents of `verify_T_case.out`. Any failed identity writes a diagnostic to
standard error, and `--quitonbreak` makes GAP exit with nonzero status.

## What is reconstructed

The verifier does not use the displayed adjoint tables as source data.
It first reconstructs the tetrahedral rack table from the four conjugation
permutations

```text
Ad(x1)=(2,4,3), Ad(x2)=(1,3,4),
Ad(x3)=(1,4,2), Ad(x4)=(1,2,3).
```

It checks the rack axioms and all 64 instances of the braid equation for
the constant `(-1)` coefficient table in the compatible bases fixed in the
manuscript. It also checks that this table is consistent with the assumed
value `epsilon=1`; it does not derive that value. A central marker `0` is
then appended to every ordered word. Actual adjacent braids in the order

```text
1,2,...,m,m,...,2,1
```

reconstruct the operation h_m. The program evaluates the compatible-basis
coefficient recursion labelled `eq:T-adjoint-recursion` (6.17):

```text
phi_m(e_i1...im)
  = e_i1...im - p e_hm(i1,...,im)
    - Prefix(i1 ▷ i2, phi_(m-1)(e_i1,i3,...,im)),

phi_1(e_i) = (1-p)e_i.
```

Polynomials in p are represented directly by integer coefficient lists.
Thus every comparison is exact.

## Identities checked

- **`eq:T-Y2-basis` and `eq:T-phi2-table` (D.2--D.3):** the four
  displayed basis vectors and all 16 entries of the complete phi_2 table,
  reduced only modulo f(p)=p^2-p+1.
- **`eq:T-y3-generator` (6.18) and `eq:T-phi3-table` (D.4):** the
  displayed generator and all 16 entries of the complete phi_3 table,
  reduced only modulo f(p).
- **`eq:T-Y4-homogeneous-components` and `eq:T-Y4-Z1` (D.5--D.6):**
  the coefficientwise identity `phi_4(w_1 tensor y)=f(p)Z_1` in the
  unreduced ring Z[p], together with the assertion that the 15 displayed
  ordered words of Z_1 are distinct. Only after this exact factorization
  is checked does the verifier reduce modulo f(p).
- **`eq:T-R2-monodromy-cycles` and `eq:T-R2-mixed-monodromy`
  (6.22--6.23):** the four three-cycles are derived by applying the eight
  adjacent braids to `[1,i,j,k,0]`; they are not used as the definition of
  the monodromy. The program counts six W-W crossings and one r*ell=p
  mixed monodromy, recovers the four displayed coefficient triples of y,
  and checks `M(y)=(1-p)y` modulo f(p).

The expected vectors and tables associated with
`eq:T-Y2-basis`--`eq:T-Y4-Z1` are transcribed from the manuscript solely
as comparison targets. All left-hand sides are generated afresh from the
rack and the compatible-basis coefficient recursion.

## Scope

This certificate verifies only the finite coefficient calculations just
listed. It does **not** derive the compatible bases, the assumed identity
`epsilon_Phi(W)=1`, the action of `x_1^{-1}` on
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
