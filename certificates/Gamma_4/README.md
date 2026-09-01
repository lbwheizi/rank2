# Gamma_4 certificates

This directory contains the exact GAP certificates accompanying Appendix C
of the paper:

- `xi3/` — reconstruction of the corrected third-adjoint common coordinate
  in Lemmas C.4--C.5 and equation (C.9),
  its reduction to the eight-atom Xi_3 target, and the final two-defect
  cancellation Xi_3 = 1;
- `local_delta/` — the corrected third-adjoint recursion and the two local
  Gamma_2 Delta identities for W in Lemma C.8 and equations (C.14)--(C.15);
- `delta_transport/` — the three displayed identities for quotients
  involving \(\Omega_{01}\) and \(\Omega_{10}\) in Lemma C.9 and
  equation (C.16);
- `terminal_character/` — reconstruction of the six action coefficients in
  Lemma C.6 and equation (C.10), and
  the source quotient, its reduction to the 76-atom target, and the
  terminal-root self-character cancellation;
- `reflected_y2/` — reconstruction of the initial and reflected induced
  actions used in Lemmas C.3 and C.7 and equation (C.11), the corrected
  first two adjoint maps and their eight-coordinate \(Y_2\) vectors,
  followed by exact integral
  source-to-target and normalized bar-boundary certificates for the seven
  reflected \(\Theta_{2,j}\) identities.  The first-root character is read
  from one coordinate; stability of that line and the deduction of
  simplicity are supplied by `lem:Gamma4-first-adjoint` in the paper.
  This suite uses abstract \(\Gamma_4\) normal forms
  and the corrected
  `c12 = a (c_WW tensor id) a^-1` convention.

Each topic directory contains its own README and checksum manifest.  From the
repository root, run `./verify_all.sh` to reproduce all checks.
