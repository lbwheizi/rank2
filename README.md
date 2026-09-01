# Exact GAP computations for rank-two coquasi Nichols algebras

This repository contains the exact GAP computations accompanying the
rank-two classification. The files are grouped by the ambient support group:

```text
certificates/
├── Gamma_2/
├── Gamma_3/
├── Gamma_4/
└── T/
```

All calculations use exact integer or finite-ring arithmetic and require only
GAP core. GAP 4.16.0 was used to produce the recorded output.

Stable LaTeX labels are the primary identifiers for correspondence with the
paper. Some public program filenames retain historical appendix numbers so
that permanent paths do not change; those numbers should not be interpreted
as the current numbering of the paper.

## Gamma_2 index

| Current paper reference and stable label | Program | Verified calculation |
|---|---|---|
| Lemma B.1, `lem:first-coefficient-reduction` | [`lemma_B1_first_coefficient_reduction.g`](certificates/Gamma_2/delta_reduction/lemma_B1_first_coefficient_reduction.g) | first-adjoint Delta reduction |
| Lemma B.3, `lem:Omega-V-one` | [`lemma_B4_Omega_V_equals_one.g`](certificates/Gamma_2/omegaV/lemma_B4_Omega_V_equals_one.g) | `Omega_V=1` |
| Lemma B.10, `lem:Gamma2-G2-cocycle-reductions` | [`lemma_B12_K3_and_d4_cocycle_reductions.g`](certificates/Gamma_2/g2_propagation/lemma_B12_K3_and_d4_cocycle_reductions.g) | `K_3=1` and `d_4=lambda^3 kappa_4` |
| Lemma B.11, `lem:Gamma2-X3-cocycle-factor-identity` | [`lemma_B14_s_phi_equals_theta_squared.g`](certificates/Gamma_2/g2_pure_parameters/s_phi_equals_theta_squared/lemma_B14_s_phi_equals_theta_squared.g) | `s_Phi=Theta_Phi^2` |
| Lemma B.12, `lem:Gamma2-R2-Theta-cocycle` | [`lemma_B15_theta_after_second_reflection.g`](certificates/Gamma_2/g2_pure_parameters/theta_r2_invariance/lemma_B15_theta_after_second_reflection.g) | reflected `Theta` identity |
| Lemma B.14, `lem:Theta-diagonal-braiding-cocycle` | [`theta_diagonal_braiding_cocycle.g`](certificates/Gamma_2/diagonal_braiding_theta/theta_diagonal_braiding_cocycle.g) | verification of the cocycle identity relating `Theta_Phi` to the diagonal braiding coefficients |
| Lemma B.15, `lem:Theta-opposite-symmetry` | [`lemma_B19_theta_symmetry.g`](certificates/Gamma_2/theta_symmetry/lemma_B19_theta_symmetry.g) | opposite symmetry of `Theta_Phi` |

Two auxiliary consistency checks are stored under
[`Gamma_2/g2_pure_parameters/auxiliary`](certificates/Gamma_2/g2_pure_parameters/auxiliary/).

## Gamma_3 index

The flat suite in [`certificates/Gamma_3`](certificates/Gamma_3/) accompanies
the calculations in Appendix E and the three obstruction propositions in the
main text.

| Current paper reference and stable label | Program | Verified calculation |
|---|---|---|
| Lemma E.2, `lem:Gamma3-32-bar-reduction`, and Proposition 7.15 | [`verify_gamma3_32.g`](certificates/Gamma_3/verify_gamma3_32.g) | normalized-bar reductions for the `(3,2)` second-adjoint obstruction |
| Proposition 7.16, `prop:Gamma3-dim2-X2-obstruction` | [`verify_gamma3_31_dim2.g`](certificates/Gamma_3/verify_gamma3_31_dim2.g) | comparison-cycle reductions for the two-dimensional `(3,1)` branch |
| Proposition 7.17, `prop:Gamma3-dim1-X2-obstruction`, including (7.36), `eq:Gamma3-dim1-m-ratios` | [`verify_gamma3_31_dim1.g`](certificates/Gamma_3/verify_gamma3_31_dim1.g) | shared comparison-cycle reductions and the scalar branch identities for the one-dimensional `(3,1)` branch |

The categorical identification in Lemma E.5,
`lem:Gamma3-dim1-chain-identification`, is proved in the paper and is not
delegated to the programs. [`verify_all_gamma3.g`](certificates/Gamma_3/verify_all_gamma3.g)
runs the three branch verifiers together with the finite `S_3` cocycle
detector. Its recorded output is
[`EXPECTED_OUTPUT.txt`](certificates/Gamma_3/EXPECTED_OUTPUT.txt).

## Gamma_4 index

| Current paper reference and stable label | Program | Verified calculation |
|---|---|---|
| Lemma C.3, `lem:Gamma4-first-adjoint`; Lemma C.7, `lem:Gamma4-R1-reflected-parameter-reduction`; (C.11) | [`gamma4_reflected_y2_certificate.g`](certificates/Gamma_4/reflected_y2/gamma4_reflected_y2_certificate.g) | source reconstruction and the seven reflected `Theta_2,j` identities |
| Lemma C.4, `lem:Gamma4-Xi3-cocycle-reduction`; Lemma C.5, `lem:Gamma4-Y3-coefficient-cancellation`; (C.9), `eq:Gamma4-Xi3-direct` | [`gamma4_xi3_certificate.g`](certificates/Gamma_4/xi3/gamma4_xi3_certificate.g) | the `Xi_3` source-to-target reduction and the final coefficient cancellation |
| Lemma C.6, `lem:Gamma4-Y2-support-value`, and (C.10) | [`gamma4_terminal_character_certificate.g`](certificates/Gamma_4/terminal_character/gamma4_terminal_character_certificate.g) | terminal-root inducing-character calculation |
| Lemma C.8, `lem:Gamma4-W-local-Deltas`, and (C.14)--(C.15) | [`gamma4_local_delta_certificate.g`](certificates/Gamma_4/local_delta/gamma4_local_delta_certificate.g) | the corrected third-adjoint recursion and the two local Delta identities for `W` |
| Lemma C.9, `lem:Gamma4-local-Delta-transport`, and (C.16) | [`gamma4_delta_transport_certificate.g`](certificates/Gamma_4/delta_transport/gamma4_delta_transport_certificate.g) | the three displayed identities for quotients involving `Omega_01` and `Omega_10` |

The local-Delta verifier reconstructs the homogeneous objects, twisted tensor
actions, braidings, and the corrected third-adjoint recursion. It obtains all
192 nonzero binomial coordinates before checking the 19 rows used by the exact
Laurent eliminations.

The reflected-`Y_2` verifier reconstructs the packet entries and compares the
source and target cocycle expressions. The assertion that the relevant
homogeneous component of the first adjoint object is a one-dimensional simple
component is the mathematical input proved in Lemma
`lem:Gamma4-first-adjoint`; it is not independently proved by the program.

## T index

The consolidated program [`verify_T_case.g`](certificates/T/verify_T_case.g)
checks the compatible-basis coefficient calculations in Section 6 and
Appendix D, together with two monodromy calculations in Lemma 6.17:

| Current reference | Stable label |
|---|---|
| (6.17) | `eq:T-adjoint-recursion` |
| (6.18) | `eq:T-y3-generator` |
| (D.2) | `eq:T-Y2-basis` |
| (D.3) | `eq:T-phi2-table` |
| (D.4) | `eq:T-phi3-table` |
| (D.5) | `eq:T-Y4-homogeneous-components` |
| (D.6) | `eq:T-Y4-Z1` |
| (6.22) | `eq:T-R2-monodromy-cycles` |
| (6.23) | `eq:T-R2-mixed-monodromy` |

For Lemma 6.17, the program checks only
`eq:T-R2-monodromy-cycles` and `eq:T-R2-mixed-monodromy`. The remaining
reflected-parameter calculations are proved in the manuscript and are not
part of the program.

The program reconstructs the tetrahedral quandle from its four conjugation
permutations, evaluates the compatible-basis coefficient recursion, verifies
the complete `phi_2` and `phi_3` tables and the unreduced fourth-adjoint
factorization, and derives the second-reflection monodromy cycles from
adjacent braids. It assumes the compatibility condition
`epsilon_Phi(W)=1`; it does not derive that condition.

## Reproduce all checks

Requirements:

- GAP 4.16.0;
- GAP core only;
- a POSIX shell and `sha256sum`.

From the repository root, run:

```sh
./verify_all.sh
```

The executable may be selected with `GAP_BIN`. If its root directory is not
detected automatically, set `GAP_ROOT` as well:

```sh
GAP_BIN=/path/to/gap GAP_ROOT=/path/to/gap-root ./verify_all.sh
```

The script checks the root `SHA256SUMS`, every topic-level
`CHECKSUMS.sha256`, runs all 16 suites with
`--bare -A -r -q --nointeract`, and accepts a suite only if GAP exits with
status zero, writes nothing to standard error, and produces standard output
that agrees byte-for-byte with the recorded output. The preliminary version
probe is subject to the same exit-status and standard-error checks and must
report exactly GAP 4.16.0. A successful run ends with:

```text
PASS: all exact certificate suites
```

## Scope

The Gamma_2 cocycle calculations use the abstract normal form subject to the
Gamma_2 relations in the paper. They verify the listed cocycle reductions;
the representation-theoretic and categorical arguments in which those
reductions are used remain in the manuscript.

The Gamma_3 suite checks normalized bar differentials, explicit integral
four-chain residuals after projection to `S_3`, the three second-adjoint
obstruction branches, and the finite `S_3` detector. Inflation along the
epimorphism, the braided monoidal equivalence determined by `J`, and the
parameter comparison remain in the written proof, as specified in the topic
README.

The Gamma_4 programs use abstract normal forms `epsilon^i h^j g^k`, reducing
only `i` modulo four and keeping `j,k` integral. They use neither
quotient-specific orders nor sampled cocycles. Each topic README states which
coefficient or cocycle reduction is verified and which categorical arguments
remain in the paper. In particular, the reflected-`Y_2` program uses the
one-dimensionality and stability of the relevant homogeneous line proved in
the manuscript; it does not establish that categorical statement by checking
a single coordinate.

The T program checks the displayed compatible-basis coefficient recursion,
the two coefficient tables, the fourth-adjoint factorization, and equations
(6.22)--(6.23). It does not construct the compatible bases, derive the action
of `x_1^{-1}` on `(W^*)_{x_1^{-1}}`, or prove the rigid-dual argument yielding
`epsilon_Phi(W^*)=1`. These points and the categorical reflection arguments
remain in the paper.
