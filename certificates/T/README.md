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
gap --bare -A -r -q --nointeract verify_T_case.g
```

The process exits with status zero after printing the contents of
`verify_T_case.out`. Any failed identity calls `Error` at the first
failing check.

## What is reconstructed

The verifier does not use the displayed adjoint tables as source data.
It first reconstructs the tetrahedral rack table from the four conjugation
permutations

```text
Ad(x1)=(2,4,3), Ad(x2)=(1,3,4),
Ad(x3)=(1,4,2), Ad(x4)=(1,2,3).
```

It checks the rack axioms, epsilon=1, and all 64 instances of the braid
equation for the constant (-1) strict braiding. A central marker `0` is
then appended to every ordered word. Actual adjacent braids in the order

```text
1,2,...,m,m,...,2,1
```

reconstruct the operation h_m. The program evaluates the strict recursion
labelled `eq:T-adjoint-recursion` (D.6):

```text
phi_m(e_i1...im)
  = e_i1...im - p e_hm(i1,...,im)
    - Prefix(i1 ▷ i2, phi_(m-1)(e_i1,i3,...,im)),

phi_1(e_i) = (1-p)e_i.
```

Polynomials in p are represented directly by integer coefficient lists.
Thus every comparison is exact.

## Identities checked

- **`eq:T-Y2-basis` and `eq:T-phi2-table` (D.7--D.8):** the four
  displayed basis vectors and all 16 entries of the complete phi_2 table,
  reduced only modulo f(p)=p^2-p+1.
- **`eq:T-y3-generator` and `eq:T-phi3-table` (D.9--D.10):** the
  displayed generator and all 16 entries of the complete phi_3 table,
  reduced only modulo f(p).
- **`eq:T-Y4-homogeneous-components` and `eq:T-Y4-Z1` (D.11--D.12):**
  the coefficientwise identity `phi_4(w_1 tensor y)=f(p)Z_1` in the
  unreduced ring Z[p], together with the assertion that the 15 displayed
  ordered words of Z_1 are distinct. Only after this exact factorization
  is checked does the verifier reduce modulo f(p).
- **`eq:T-R2-monodromy-cycles` and `eq:T-R2-mixed-monodromy`
  (6.21--6.22):** the four three-cycles are derived by applying the eight
  adjacent braids to `[1,i,j,k,0]`; they are not used as the definition of
  the monodromy. The program counts six W-W crossings and one r*ell=p
  mixed monodromy, recovers the four displayed coefficient triples of y,
  and checks `M(y)=(1-p)y` modulo f(p).

The expected vectors and tables associated with
`eq:T-Y2-basis`--`eq:T-Y4-Z1` are transcribed from the manuscript solely
as comparison targets. All left-hand sides are generated afresh from the
rack and the strict recursion.

## Scope

This certificate verifies only the finite coefficient calculations just
listed. It does **not** prove the theoretical passage from the original
right-associated coquasi category to the canonical strictification, the
normal-form lemma for arbitrary b_ij, simplicity of the adjoint objects,
existence or involutivity of reflections, standardness of the Cartan graph,
root factorization, or finite-dimensionality of a Nichols algebra. Those
are mathematical arguments in the manuscript. In particular, the
cocycle reductions surrounding Lemmas D.1--D.3 and
`lem:T-reflection-cocycle-reductions` are not independently verified here.
