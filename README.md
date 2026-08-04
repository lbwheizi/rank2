# Exact certificates for the rank-two coquasi Nichols algebra project

This repository contains computer-verifiable certificates accompanying the
classification of finite-dimensional rank-two Nichols algebras over
coquasi-Hopf algebras.  Every program uses exact integer or finite-ring
arithmetic; no numerical approximation is used.

## Certificate index

| Directory | Main identity or check |
|---|---|
| [`delta_reduction`](certificates/delta_reduction/) | first-adjoint $\Delta$ reduction |
| [`omegaV`](certificates/omegaV/) | $\Omega_V=1$ |
| [`theta_symmetry`](certificates/theta_symmetry/) | $\Theta_\Phi(h,g)=\Theta_\Phi(g,h)$ |
| [`g2_propagation`](certificates/g2_propagation/) | $\mathcal K_3=1$ and $d_4=\lambda^3\kappa_4$ |
| [`g2_pure_parameters/xi_equals_theta_squared`](certificates/g2_pure_parameters/xi_equals_theta_squared/) | $\Xi_\Phi=\Theta_\Phi^2$ |
| [`g2_pure_parameters/theta_r2_invariance`](certificates/g2_pure_parameters/theta_r2_invariance/) | $\Theta_\Phi(hg,h^{-1})=\Theta_\Phi(g,h)$ |
| [`g2_pure_parameters/xi_r2_invariance`](certificates/g2_pure_parameters/xi_r2_invariance/) | auxiliary check $\Xi_\Phi(hg,h^{-1})=\Xi_\Phi(g,h)$ |
| [`g2_pure_parameters/xi_square`](certificates/g2_pure_parameters/xi_square/) | auxiliary check $\Xi_\Phi^2=1$ |
| [`order16_G2`](certificates/order16_G2/) | exhaustive finite-data checks for the order-$16$ $G_2$ example |

Here $\Xi_\Phi$ is the internal certificate name for the third-string scalar
denoted $\mathfrak s_\Phi$ in the paper.  The identity
$\Xi_\Phi=\Theta_\Phi^2$ shows that the oriented $G_2$ parameter locus has
only the three conditions

$$
\Delta=1,\qquad \lambda'=-1,\qquad
\lambda^2=\Theta_\Phi=-1.
$$

## Requirements

- GAP 4.16.0 (the version used for every recorded run);
- GAP core only: no optional package is required;
- a POSIX shell for the convenience script `verify_all.sh`.

## Reproduce all checks

From the repository root, run

```sh
./verify_all.sh
```

Each directory also contains its own `README.md`, recorded output, and
`CHECKSUMS.sha256`.  The convenience script first verifies the checksums,
then runs every GAP program and compares the fresh output byte-for-byte with
the recorded output.

## Scope

The cocycle scripts expand local factors into the free abelian group on
normalized symbols $\Phi(a,b,c)$ and check exact integral combinations of
normalized $3$-cocycle relators.  The order-$16$ verifier separately checks
the finite group, all $16^4$ cocycle identities, projective characters,
the printed twelve-state support graph, its local scalar data, and the root
label arithmetic.  Its README records the precise boundary of the finite
checks; in particular, it does not replace the recursive adjoint calculation
in the paper.
