# Exact GAP certificates for the rank-two coquasi Nichols algebra project

This repository contains the GAP certificates accompanying the rank-two
classification in the paper.  All calculations use exact integer or
finite-ring arithmetic.  The scripts require only GAP core; no optional GAP
package is loaded.

## Certificate index

The program names reproduce the numbered statement in Appendix B of the
paper.  The LaTeX labels are included so that the correspondence remains
unambiguous if the printed numbers change in a later version.

| Paper reference and LaTeX label | Program | Identity or finite check |
|---|---|---|
| Lemma B.1, `lem:first-coefficient-reduction` | [`lemma_B1_first_coefficient_reduction.g`](certificates/delta_reduction/lemma_B1_first_coefficient_reduction.g) | first-adjoint $\Delta$ reduction |
| Lemma B.4, `lem:Omega-V-one` | [`lemma_B4_Omega_V_equals_one.g`](certificates/omegaV/lemma_B4_Omega_V_equals_one.g) | $\Omega_V=1$ |
| Lemma B.12, `lem:Gamma2-G2-cocycle-reductions` | [`lemma_B12_K3_and_d4_cocycle_reductions.g`](certificates/g2_propagation/lemma_B12_K3_and_d4_cocycle_reductions.g) | $\mathcal K_3=1$ and $d_4=\lambda^3\kappa_4$ |
| Lemma B.14, `lem:Gamma2-third-string-theta-square-cocycle` | [`lemma_B14_s_phi_equals_theta_squared.g`](certificates/g2_pure_parameters/s_phi_equals_theta_squared/lemma_B14_s_phi_equals_theta_squared.g) | $s_\Phi=\Theta_\Phi^2$ |
| Lemma B.15, `lem:Gamma2-R2-Theta-cocycle` | [`lemma_B15_theta_after_second_reflection.g`](certificates/g2_pure_parameters/theta_r2_invariance/lemma_B15_theta_after_second_reflection.g) | $\Theta_\Phi(hg,h^{-1})=\Theta_\Phi(g,h)$ |
| Lemma B.19, `lem:Theta-opposite-symmetry` | [`lemma_B19_theta_symmetry.g`](certificates/theta_symmetry/lemma_B19_theta_symmetry.g) | $\Theta_\Phi(h,g)=\Theta_\Phi(g,h)$ |
| Remark B.21, `rem:order16-G2-certificate` | [`remark_B21_order16_G2_verifier.g`](certificates/order16_G2/remark_B21_order16_G2_verifier.g) | finite-data checks for the order-$16$ $G_2$ example |

Two additional consistency checks, which are not cited as independent
statements in the paper, are stored under
[`g2_pure_parameters/auxiliary`](certificates/g2_pure_parameters/auxiliary/).

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
