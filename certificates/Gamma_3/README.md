# Exact GAP certificates for the coquasi \(\Gamma _3\) calculation

These files verify the coefficient, action, scalar-to-chain, and integer
normalized-bar-complex calculations used in the three local
\(X_2\)-obstruction arguments for
\[
\Gamma _3=
\langle\varepsilon,g,h\mid
\varepsilon^3=1,
g\varepsilon g^{-1}=\varepsilon^{-1},
[h,g]=[h,\varepsilon]=1\rangle.
\]

The scripts use only the GAP core library and exact integer, rational,
cyclotomic, polynomial, or finite-field arithmetic. No optional GAP
package, Python, NumPy, SciPy, or integer-programming solver is used by
the verifier. Atomic values \(\Phi(a,b,c)\) and representation scalars
are kept formal. Local cocycles \(\Phi_x\), the factors \(\Phi^x\) in
`eq:tensor-product`, and the factors \(R(a;x,y)\) are expanded from
their definitions.

The original Python programs were used only to find and audit candidate
four-chains. The four-chains are hard-coded here, and GAP reconstructs
the source residuals and recomputes every integer boundary coefficient.

## Files

- `bar_complex.g`: normalized inhomogeneous bar-chain arithmetic.
- `verify_gamma3_32.g`: the \((3,2)\) second-adjoint coefficient block,
  scalar-to-chain dictionary, and bar reduction in
  `lem:Gamma3-32-second-adjoint-expansion`,
  `eq:Gamma3-32-scalar-chain-dictionary`, and
  `lem:Gamma3-32-bar-reduction`.
- `verify_gamma3_31_dim2.g`: the \((3,1)_2\) coefficient calculation
  in `lem:Gamma3-dim2-epsilon-cubic`,
  `lem:Gamma3-dim2-second-adjoint-expansion`, and
  `lem:Gamma3-dim2-bar-reduction`, including the formulas from
  `eq:Gamma3-dim2-shifted-X1` through
  `eq:Gamma3-dim2-eta-values` and the two comparison-cycle reductions
  leading to `eq:Gamma3-common-three-path-reduction`.
- `verify_gamma3_31_dim1.g`: the \((3,1)_1\) calculation in
  `lem:Gamma3-dim1-three-elementary-tensors`, the chain identification
  in `lem:Gamma3-dim1-chain-identification`, the formulas relating
  \(m_i\) and \(\mu_i\), and the scalar-ratio and branch-polynomial
  reductions. It first reruns the shared \((3,1)_2\) bar certificate.
- `verify_s3_cocycle_mod3.g`: an independent normalized
  \(\mathbb F_3\)-valued \(3\)-cocycle on \(S_3\), with kappa exponent
  one.
- `verify_all_gamma3.g`: runs the complete suite and checks every
  component success flag.
- `EXPECTED_OUTPUT.txt`: recorded successful output.
- `CHECKSUMS.sha256`: SHA-256 checksums of all preceding deliverables.

## Running the verifier

Run GAP from this directory. With GAP 4.16.0 the batch command is

```sh
gap --bare -A -r -q --quitonbreak --norepl verify_all_gamma3.g
```

For an interactive run, start `gap --bare -A -r -q`, then enter

```gap
Read("verify_all_gamma3.g");
```

Success returns to `gap>`. Failure calls `ErrorNoReturn` and stops at
`brk>`, where `quit;` returns to the outer prompt. No script calls
`QUIT_GAP`.

For an exact comparison with the recorded batch output, use

```sh
gap --bare -A -r -q --quitonbreak --norepl verify_all_gamma3.g > actual_output.txt
diff -u EXPECTED_OUTPUT.txt actual_output.txt
```

The batch verifier exits with status zero only after every component
success flag has been checked. A failed assertion is written to
standard error, `--quitonbreak` gives a nonzero status, and the final
total-PASS line is unreachable.

Each component may also be run separately:

```sh
gap --bare -A -r -q --quitonbreak --norepl verify_gamma3_32.g
gap --bare -A -r -q --quitonbreak --norepl verify_gamma3_31_dim2.g
gap --bare -A -r -q --quitonbreak --norepl verify_gamma3_31_dim1.g
gap --bare -A -r -q --quitonbreak --norepl verify_s3_cocycle_mod3.g
```

The one-dimensional entry point intentionally prints the
\((3,1)_2\) PASS lines first because it reruns the shared chains.

## What is checked exactly

### The \((3,2)\) calculation

`verify_gamma3_32.g` starts from the normal form in \(\Gamma _3\), the
transversals defining the induced objects, `eq:YD-projective-action`,
`eq:tensor-product`, the braiding formula, and
`eq:recursive-varphi-2`.
It reconstructs:

1. the first-adjoint coefficients and the relevant induced actions;
2. the central and noncentral second-adjoint paths and their degrees;
3. the cyclic action coefficients \(t_0,t_1,t_2\);
4. the scalar-to-chain dictionary from the original scalar formulas;
5. the evaluations of
   \(\mathsf C_s,\mathsf C_\lambda,\mathsf C_r\).

Three cocycle comparisons which do not cancel formally are certified
by explicit integral four-chains. Their reconstructed residuals have
respectively \(33,34,42\) terms. The program also verifies that the
three comparison cycles have zero normalized \(\Gamma _3\)
differential and term counts \(25,44,53\). Their projections to
\(S_3\) are the boundaries of the recorded four-chains with
respectively \(5,24,22\) terms. Two independent one-factor
perturbations of source expressions are required to leave nonzero
residuals.

These computations cover the coefficient block in
`lem:Gamma3-32-second-adjoint-expansion` and the scalar-chain bridge to
the three identities in `lem:Gamma3-32-bar-reduction`.

### The \((3,1)_2\) calculation

`verify_gamma3_31_dim2.g` retains the original reconstruction of
\[
\mathsf E,\mathsf A,\mathsf B,\mathsf T_i,\mathsf M_i,
\mathsf\Theta,\mathsf K,
\mathsf C_{02},\mathsf C_{\mathrm{pow}}
\]
and of
\(\mathsf C_{\mathrm{eigen}},\mathsf R,\mathsf R_A\). It checks all
displayed differentials and cycles. The projected integral
four-chain boundaries for \(\mathsf C_{02}\) and
\(\mathsf C_{\mathrm{pow}}\) have \(8,9\) terms; the corresponding
\(\Gamma _3\)-chains have \(38,20\) terms.

The added coefficient layer begins with the two rows of the action
table in `eq:Gamma3-dim2-standard-basis`. Under
\(\Delta _3=I\) and \(\rho(g)=-1\), it independently reconstructs:

1. all three shifted \(\varphi _1^\Phi\) formulas;
2. the double-braiding and final \(c_{1,2}^\Phi\) paths in
   \(\varphi _2^\Phi(e_0)\);
3. the common degree of \(e_0,e_1,e_2\);
4. the cyclic action, including its three target indices and three
   coefficients, from the transversals, `eq:YD-projective-action`, and
   `eq:tensor-product`;
5. the evaluations of
   \(\mathsf E,\mathsf R,\mathsf A,\mathsf B,\mathsf T_i,
   \mathsf M_i,\mathsf\Theta,\mathsf K\) and both comparison cycles.

The last exact calculation shows that the quotient between the derived
value of \(m_1\) and \(m_0\kappa^{-1}\) is \(\kappa^3\). Thus the
previously proved relation \(\kappa^3=1\) gives the common reduction.
For the negative control, a new factor is inserted into the source
\(c_{1,2}^\Phi\)-coefficient, the complete last path is recomputed, and
the altered result is required to differ from \(\mathfrak b\).

### The \((3,1)_1\) calculation

The one-dimensional branch uses the same bar-chain definitions as the
two-dimensional branch, so `verify_gamma3_31_dim1.g` first reruns
`verify_gamma3_31_dim2.g`. It then reconstructs its own coefficients
from `eq:YD-projective-action`, `eq:tensor-product`, the braiding
formula, and `eq:recursive-varphi-2`. It checks:

1. the two nontrivial coefficients in
   `lem:Gamma3-dim1-three-elementary-tensors`;
2. all three target indices and coefficients in the cyclic action;
3. the atomic-cocycle evaluations of
   \(\mathsf E,\mathsf A,\mathsf B,\mathsf T_i,\mathsf M_i,
   \mathsf\Theta,\mathsf K\);
4. the three formulas relating \(m_i\) and \(\mu_i\);
5. the two comparison-cycle evaluations and the common-\(\mu\)
   reduction modulo \(\kappa^3=1\).

It finally checks over an exact polynomial ring that
\[
\frac{m_0}{m_1}=
\kappa\frac{1-t}{p^3t^2},\qquad
\frac{m_0}{m_2}=
-\frac{\kappa}{p^3t(1-t)},\qquad
\frac{m_1}{m_2}=-\frac{t}{(1-t)^2},
\]
after cancellation of the common nonzero \(\mu_0\), together with the
branch-polynomial identities. Its negative control perturbs the source
\(c_{1,2}^\Phi\)-coefficient and recomputes the last path.

### The \(S_3\) detector

`verify_s3_cocycle_mod3.g` checks that the recorded triples are distinct
and have nonzero values. It evaluates the cochain on all \(6^3\)
triples, reconstructs its actual support of size \(32\), checks all
\(6^4=1296\) normalized three-cocycle equations, and verifies
\[
f(\varepsilon,\varepsilon,\varepsilon)
+f(\varepsilon,\varepsilon^2,\varepsilon)=1\in\mathbb F_3.
\]

## Written-proof inputs and scope limits

The coefficient verifiers start from `eq:YD-projective-action`,
`eq:tensor-product`, the braiding formula, and
`eq:recursive-varphi-2`. They do not use the desired coefficient
identities as inputs. The following statements are, however,
written-proof inputs or remain outside the programs:

- the homology computations
  \(H_3(\Gamma _3,\mathbb Z)\simeq\mathbb Z/3\) and the injectivity of
  \(H_3(\Gamma _3,\mathbb Z)\to H_3(S_3,\mathbb Z)\); the scripts check
  the cycles and the projected four-chain boundaries, while the
  homological lemma justifies the inference back on \(\Gamma _3\);
- in the \((3,2)\) branch, the previously proved
  \(G\)-linearity of \(\varphi _1^\Phi\), \(\Delta _1=1\), and, only in
  the final scalar dictionary, \(\rho(g)=-1\);
- in the \((3,1)_2\) branch, the normal-form action table,
  \(\Delta _3=I\), \(\rho(g)=-1\), and the homological consequence that
  the verified \(\mathsf R\)-cycle evaluates to one;
- in the \((3,1)_1\) branch,
  `eq:Gamma3-dim1-epsilon-normalization` and
  `eq:Gamma3-dim1-epsilon-square`; the program explicitly uses them as
  rewrite rules and does not claim to prove them;
- `eq:YD-projective-action`, applied to three successive
  \(\varepsilon\)-actions with
  \(\Phi_{\varepsilon g^2h}(\varepsilon,\varepsilon)
  \Phi_{\varepsilon g^2h}(\varepsilon^2,\varepsilon)\), and the earlier
  relation \(\kappa_\Phi(\varepsilon)^3=1\), used in the last common-path
  reduction;
- the alternating-tricharacter and projective-dimension arguments;
- simplicity or non-simplicity of any Yetter--Drinfeld object;
- inflation along the epimorphism, the braided monoidal equivalence
  determined by \(J\), or the parameter comparison;
- existence of reflections, standardness of the Cartan graph, root
  factorization, or finite-dimensionality.

Thus the machine-verifiable claims are precisely the stated exact
coefficient, degree, scalar-to-chain, integer-boundary, and polynomial
residual calculations. The listed structural and homological
implications remain part of the written proof.

## Relation to the construction-stage audit

The GAP suite is a self-contained exact port of the computations used
during the construction-stage audit. The Python finder files are not
part of this deliverable and are not required for reproduction. The
present package retains only the resulting explicit four-chains and
independently reconstructs their source residuals and boundaries.

The one-dimensional and two-dimensional central-support branches have
identical bar-chain definitions. The one-dimensional entry point
therefore reruns the shared certificate, but it reconstructs its own
coefficient paths and independently checks its scalar dictionary and
scalar-ratio consequences.
