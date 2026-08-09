# Gamma_4 certificates

This directory contains the exact GAP certificates accompanying Appendix C
of the paper:

- `xi3/` — reconstruction of the corrected third-adjoint common coordinate,
  its reduction to the eight-atom Xi_3 target, and the final two-defect
  cancellation Xi_3 = 1;
- `local_delta/` — the corrected third-adjoint recursion and the two local
  Gamma_2 Delta identities for W;
- `delta_transport/` — transport of the local Gamma_2 constants from W to
  Y_1;
- `terminal_character/` — reconstruction of the six action coefficients and
  the source quotient, its reduction to the 76-atom target, and the
  terminal-root self-character cancellation;
- `reflected_y2/` — an independent reconstruction of the initial and
  reflected induced actions, the corrected first two adjoint maps and their
  eight-coordinate (Y_2) vectors, followed by seven exact integral
  source-to-target and normalized bar-boundary certificates proving
  simplicity of
  \((\operatorname{ad}X_1)^2(V^*)\) after the first reflection.  This suite
  uses abstract \(\Gamma_4\) normal forms and the corrected
  `c12 = a (c_WW tensor id) a^-1` convention.

Each topic directory contains its own README and checksum manifest.  From the
repository root, run `./verify_all.sh` to reproduce all checks.
