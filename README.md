# Exact GAP certificates for the rank-two coquasi Nichols algebra project

This repository contains the GAP certificates accompanying the rank-two
classification in the paper.  All calculations use exact integer or
finite-ring arithmetic.  The scripts require only GAP core; no optional GAP
package is loaded.

## Certificate index

| Directory | Identity or finite check |
|---|---|
| [`delta_reduction`](certificates/delta_reduction/) | first-adjoint $\Delta$ reduction |
| [`omegaV`](certificates/omegaV/) | $\Omega_V=1$ |
| [`theta_symmetry`](certificates/theta_symmetry/) | $\Theta_\Phi(h,g)=\Theta_\Phi(g,h)$ |
| [`g2_propagation`](certificates/g2_propagation/) | $\mathcal K_3=1$ and $d_4=\lambda^3\kappa_4$ |
| [`g2_pure_parameters/s_phi_equals_theta_squared`](certificates/g2_pure_parameters/s_phi_equals_theta_squared/) | $s_\Phi=\Theta_\Phi^2$ |
| [`g2_pure_parameters/theta_r2_invariance`](certificates/g2_pure_parameters/theta_r2_invariance/) | $\Theta_\Phi(hg,h^{-1})=\Theta_\Phi(g,h)$ |
| [`g2_pure_parameters/auxiliary/s_phi_r2_invariance`](certificates/g2_pure_parameters/auxiliary/s_phi_r2_invariance/) | auxiliary check $s_\Phi(hg,h^{-1})=s_\Phi(g,h)$ |
| [`g2_pure_parameters/auxiliary/s_phi_square`](certificates/g2_pure_parameters/auxiliary/s_phi_square/) | auxiliary check $s_\Phi^2=1$ |
| [`order16_G2`](certificates/order16_G2/) | finite-data checks for the order-$16$ $G_2$ example |

The identity $s_\Phi=\Theta_\Phi^2$ holds in the abstract $\Gamma_2$
setting, without imposing the $G_2$ parameter conditions.  It is the
calculation used to remove $s_\Phi=1$ from the list of independent
conditions on the oriented $G_2$ parameter locus.  The two directories
under `g2_pure_parameters/auxiliary` contain consistency checks not used as
additional assumptions in the paper.

## Requirements

- GAP 4.16.0 (the version used to produce the recorded output);
- GAP core only;
- a POSIX shell and `sha256sum` for `verify_all.sh`.

## Reproduce all checks

From the repository root, run

```sh
./verify_all.sh
```

The GAP executable may be selected with `GAP_BIN`.  If its root directory is
not detected automatically, set `GAP_ROOT` as well:

```sh
GAP_BIN=/path/to/gap GAP_ROOT=/path/to/gap-root ./verify_all.sh
```

The verifier checks `SHA256SUMS`, checks every directory's
`CHECKSUMS.sha256`, runs each program with
`--bare -A -r -q --nointeract`, and compares the fresh output byte-for-byte
with the recorded output.

## Scope

The cocycle certificates work in the abstract normal form
$\varepsilon^e g^a h^b$ subject only to
$hg=\varepsilon gh$, $\varepsilon^2=1$, and centrality of $\varepsilon$.
They expand products of normalized $3$-cocycle values into exact integral
exponent vectors and reduce them by explicitly listed cocycle relators.

The `order16_G2` script instead checks the explicit finite data of one
order-$16$ quotient: the group and cocycle identities, the stated projective
characters, the twelve printed support states and their transitions, and the
root-label arithmetic.  As its directory README explains, it does not
replace the recursive adjoint calculations proved in the paper.
