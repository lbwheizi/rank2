# Exact GAP certificates for the coquasi \(\Gamma _3\) calculation

These files verify the integer normalized-bar-complex calculations used
in the three local \(X_2\)-obstruction arguments for
\[
\Gamma _3=
\langle\varepsilon,g,h\mid
\varepsilon^3=1,
g\varepsilon g^{-1}=\varepsilon^{-1},
[h,g]=[h,\varepsilon]=1\rangle.
\]

The scripts use only the GAP core library and exact integer, rational,
or finite-field arithmetic.  No optional GAP package, Python, NumPy,
SciPy, or integer-programming solver is used by the verifier.  The
original Python programs were used to construct and audit candidate
four-chains; the four-chains are hard-coded here and GAP recomputes
their complete integer boundary residuals.

## Files

- `bar_complex.g`: normalized inhomogeneous bar-chain arithmetic.
- `verify_gamma3_32.g`: the \((3,2)\) bar reduction in Lemma E.2
  (`lem:Gamma3-32-bar-reduction`), used in the obstruction argument of
  Proposition 7.15 (`prop:Gamma3-32-X2-obstruction`).
- `verify_gamma3_31_dim2.g`: the \((3,1)_2\) central two-dimensional
  obstruction in Proposition 7.16
  (`prop:Gamma3-dim2-X2-obstruction`).
- `verify_gamma3_31_dim1.g`: the \((3,1)_1\) central one-dimensional
  obstruction in Proposition 7.17
  (`prop:Gamma3-dim1-X2-obstruction`).  It first reruns the shared bar
  certificate from the preceding file and then verifies the three
  scalar-ratio reductions in (7.36)
  (`eq:Gamma3-dim1-m-ratios`) and the branch polynomial identities.
- `verify_s3_cocycle_mod3.g`: an independent normalized
  \(\mathbb F_3\)-valued \(3\)-cocycle on \(S_3\), with kappa exponent
  one.
- `verify_all_gamma3.g`: runs the complete suite.
- `EXPECTED_OUTPUT.txt`: recorded successful output.
- `CHECKSUMS.sha256`: SHA-256 checksums of all preceding deliverables.

## Running the verifier

Run GAP from this directory.  With GAP 4.16.0 the command is

```sh
gap --bare -A -r -q --nointeract verify_all_gamma3.g
```

For an exact comparison with the recorded output:

```sh
gap --bare -A -r -q --nointeract verify_all_gamma3.g > actual_output.txt
diff -u EXPECTED_OUTPUT.txt actual_output.txt
```

The complete verifier exits with status zero only after all four
component success flags have been checked.  Any failed assertion is
written to standard error and terminates GAP with a nonzero status;
the final total-PASS line is then unreachable.

Each component may also be run separately:

```sh
gap --bare -A -r -q --nointeract verify_gamma3_32.g
gap --bare -A -r -q --nointeract verify_gamma3_31_dim2.g
gap --bare -A -r -q --nointeract verify_gamma3_31_dim1.g
gap --bare -A -r -q --nointeract verify_s3_cocycle_mod3.g
```

The one-dimensional script intentionally prints the
\((3,1)_2\) PASS lines first, because it reuses and reruns exactly the
same chains rather than merely asserting that they were checked
elsewhere.

## What is checked exactly

### The \((3,2)\) calculation

`verify_gamma3_32.g` reconstructs from the group normal form all local
chains \(L,U,R\) and the three cycles
\[
\mathsf C_s,qquad \mathsf C_\lambda,qquad \mathsf C_r.
\]
It checks:

1. their normalized \(\Gamma _3\) differentials are exactly zero;
2. their term counts are respectively \(25,44,53\);
3. the boundaries of explicit integral four-chains in
   \(B_4(S_3,\mathbb Z)\) equal their projected chains coefficient by
   coefficient, with zero residual; the four-chains have respectively
   \(5,24,22\) nonzero terms.

These are the bar identities underlying
\(s^3=\kappa_\Phi(\varepsilon)\), the local cubic relation, and the
remaining path-ratio reduction in the \((3,2)\) proof.

### The \((3,1)_2\) calculation

`verify_gamma3_31_dim2.g` reconstructs
\[
\mathsf E,\mathsf A,\mathsf B,\mathsf T_i,\mathsf M_i,
\mathsf\Theta,\mathsf K,
\mathsf C_{02},\mathsf C_{\mathrm{pow}}
\]
and also the auxiliary reductions
\(\mathsf C_{\mathrm{eigen}},\mathsf R,\mathsf R_A\).  It checks all
displayed differentials, verifies that every named difference is a
cycle, and verifies explicit projected integral four-chain boundaries
for \(\mathsf C_{02}\) and \(\mathsf C_{\mathrm{pow}}\).  Their
\(\Gamma _3\)-term counts are \(38,20\), and the two four-chains have
\(8,9\) terms.

### The \((3,1)_1\) calculation

The one-dimensional calculation uses *literally the same* definitions
of
\[
\mathsf E,\mathsf A,\mathsf B,\mathsf T_i,\mathsf M_i,
\mathsf\Theta,\mathsf K,
\mathsf C_{02},\mathsf C_{\mathrm{pow}}
\]
as the two-dimensional calculation.  This is why
`verify_gamma3_31_dim1.g` reruns
`verify_gamma3_31_dim2.g`.  It then checks, over an exact polynomial
ring, the reductions
\[
\frac{m_0}{m_1}=
\kappa\frac{1-t}{p^3t^2},\qquad
\frac{m_0}{m_2}=
-\frac{\kappa}{p^3t(1-t)},\qquad
\frac{m_1}{m_2}=-\frac{t}{(1-t)^2},
\]
after cancellation of the common nonzero \(\mu _0\), together with
the elementary polynomial consequences used for the two branches.

### The \(S_3\) detector

`verify_s3_cocycle_mod3.g` first checks that the recorded triples are
distinct and have nonzero values.  It then evaluates the cochain on
all \(6^3\) triples, reconstructs its actual support, and checks that
this support is exactly the recorded set of size \(32\).  Finally, it
checks all \(6^4=1296\) normalized three-cocycle equations and checks
that
\[
f(\varepsilon,\varepsilon,\varepsilon)
+f(\varepsilon,\varepsilon^2,\varepsilon)=1\in\mathbb F_3.
\]

## What is not proved by these programs

The GAP programs deliberately do **not** prove any of the following:

- that the displayed bar chains are the categorical coefficients of
  the braided adjoint maps (this is the hand derivation in the proof);
- the homology computations
  \(H_3(\Gamma _3,\mathbb Z)\simeq\mathbb Z/3\) and the injectivity of
  \(H_3(\Gamma _3,\mathbb Z)\to H_3(S_3,\mathbb Z)\); the scripts check
  the cycles and the projected four-chain boundaries, while these
  homological lemmas justify the inference back on \(\Gamma _3\);
- the alternating-tricharacter and projective-dimension arguments;
- simplicity or non-simplicity of any Yetter--Drinfeld object;
- the inflation along the epimorphism, the braided monoidal equivalence
  determined by \(J\), or the parameter comparison;
- existence of reflections, standardness of the Cartan graph, root
  factorization, or finite-dimensionality.

Thus the machine-verifiable claim is precisely the exact integer
bar-chain and scalar-polynomial residual computation.  The surrounding
categorical and homological implications remain part of the written
mathematical proof.

## Relation to the construction-stage audit

The GAP suite is a self-contained exact port of the bar-chain checks
used during the construction-stage audit.  The Python finder files are
not part of this deliverable and are not required for reproduction.
The present package removes the Python solver dependency by retaining
only the resulting explicit four-chains and independently recomputing
their boundaries over \(\mathbb Z\).

The one-dimensional and two-dimensional central-support calculations
have identical bar-chain definitions; this equality is stated and
proved in Lemma E.5 (`lem:Gamma3-dim1-chain-identification`) of the
accompanying mathematical text.  The one-dimensional GAP entry point
reruns the shared certificate and separately verifies its own
scalar-ratio consequences.
