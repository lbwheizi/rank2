# Exact GAP certificates for rank-two coquasi Nichols algebras

This repository contains the exact GAP certificates accompanying the
rank-two classification.  Certificates are grouped first by the ambient
support group, so that every calculation can be located from its group case:

```text
certificates/
├── Gamma_2/
├── Gamma_3/
├── Gamma_4/
└── T/
```

All calculations use exact integer or finite-ring arithmetic and require
only GAP core; no optional GAP package is loaded.  `Gamma_3/` and `T/` are
currently index placeholders for future certificates.

## Gamma_2 certificate index

The program names reproduce the numbered statements in Appendix B.  Stable
LaTeX labels are included so that the correspondence remains unambiguous if
printed numbering changes.

| Paper reference and LaTeX label | Program | Identity or finite check |
|---|---|---|
| Lemma B.1, `lem:first-coefficient-reduction` | [`lemma_B1_first_coefficient_reduction.g`](certificates/Gamma_2/delta_reduction/lemma_B1_first_coefficient_reduction.g) | first-adjoint Delta reduction |
| Lemma B.4, `lem:Omega-V-one` | [`lemma_B4_Omega_V_equals_one.g`](certificates/Gamma_2/omegaV/lemma_B4_Omega_V_equals_one.g) | Omega_V = 1 |
| Lemma B.12, `lem:Gamma2-G2-cocycle-reductions` | [`lemma_B12_K3_and_d4_cocycle_reductions.g`](certificates/Gamma_2/g2_propagation/lemma_B12_K3_and_d4_cocycle_reductions.g) | K_3 = 1 and d_4 = lambda^3 kappa_4 |
| Lemma B.14, `lem:Gamma2-third-string-theta-square-cocycle` | [`lemma_B14_s_phi_equals_theta_squared.g`](certificates/Gamma_2/g2_pure_parameters/s_phi_equals_theta_squared/lemma_B14_s_phi_equals_theta_squared.g) | s_Phi = Theta_Phi^2 |
| Lemma B.15, `lem:Gamma2-R2-Theta-cocycle` | [`lemma_B15_theta_after_second_reflection.g`](certificates/Gamma_2/g2_pure_parameters/theta_r2_invariance/lemma_B15_theta_after_second_reflection.g) | reflected Theta invariance |
| Lemma B.19, `lem:Theta-opposite-symmetry` | [`lemma_B19_theta_symmetry.g`](certificates/Gamma_2/theta_symmetry/lemma_B19_theta_symmetry.g) | opposite symmetry of Theta_Phi |
| Remark B.21, `rem:order16-G2-certificate` | [`remark_B21_order16_G2_verifier.g`](certificates/Gamma_2/order16_G2/remark_B21_order16_G2_verifier.g) | finite-data checks for the order-16 G_2 example |

Two additional consistency checks are stored under
[`Gamma_2/g2_pure_parameters/auxiliary`](certificates/Gamma_2/g2_pure_parameters/auxiliary/).

## Gamma_4 certificate index

The following five directories accompany Appendix C.

| Paper reference and LaTeX label | Program | Identity or coefficient check |
|---|---|---|
| Remark C.5, `rem:Gamma4-Xi3-certificate` | [`gamma4_xi3_certificate.g`](certificates/Gamma_4/xi3/gamma4_xi3_certificate.g) | final two-defect cancellation in Xi_3 = 1 |
| Remark C.8, `rem:Gamma4-local-Delta-certificate` | [`gamma4_local_delta_certificate.g`](certificates/Gamma_4/local_delta/gamma4_local_delta_certificate.g) | corrected phi_3 recursion and Y_3 = 0 implies both local Delta identities for W |
| Remark C.10, `rem:Gamma4-Omega-certificate` | [`gamma4_delta_transport_certificate.g`](certificates/Gamma_4/delta_transport/gamma4_delta_transport_certificate.g) | transport of the four local Gamma_2 constants |
| Remark C.12, `rem:Gamma4-terminal-certificate` | [`gamma4_terminal_character_certificate.g`](certificates/Gamma_4/terminal_character/gamma4_terminal_character_certificate.g) | terminal-root self-character cancellation |
| Reflected-Y2 lemma, `lem:Gamma4-R1-reflected-parameter-reduction` | [`gamma4_reflected_y2_certificate.g`](certificates/Gamma_4/reflected_y2/gamma4_reflected_y2_certificate.g) | source reconstruction and seven integral parameter identities proving the reflected second adjoint simple after (R_1) |

The local-Delta verifier reconstructs the homogeneous objects, twisted
tensor actions, braidings, and the corrected third-adjoint recursion inside
GAP.  It obtains all 192 nonzero binomial coordinates before checking the 19
rows used by the exact Laurent eliminations.

## Requirements

- GAP 4.16.0 (the version used to produce the recorded output);
- GAP core only;
- a POSIX shell and `sha256sum` for `verify_all.sh`.

## Reproduce all checks

From the repository root, run:

```sh
./verify_all.sh
```

The GAP executable may be selected with `GAP_BIN`.  If its root directory is
not detected automatically, set `GAP_ROOT` as well:

```sh
GAP_BIN=/path/to/gap GAP_ROOT=/path/to/gap-root ./verify_all.sh
```

The verifier checks the root `SHA256SUMS`, every topic-level
`CHECKSUMS.sha256`, executes each program with
`--bare -A -r -q --nointeract`, and compares fresh output byte-for-byte with
the recorded output.

## Scope

The Gamma_2 cocycle certificates work in the abstract normal form subject to
the Gamma_2 relations stated in the paper.  The order-16 script instead
checks the explicit finite data of that quotient.

The Gamma_4 certificates use abstract normal forms
`epsilon^i h^j g^k`, reducing only `i` modulo four and keeping `j,k` integral.
They use neither quotient-specific orders nor sampled cocycles.  Each topic
README states exactly which coefficient or cocycle reduction is verified and
which categorical arguments remain in the text.
