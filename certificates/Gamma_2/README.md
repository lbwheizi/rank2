# Gamma_2 certificates

This directory contains the exact GAP certificates for the `Gamma_2` case
and Appendix B of the paper.  The existing certificate suites are grouped by
the calculation they verify:

- `delta_reduction/` — first-adjoint Delta reduction;
- `omegaV/` — the identity Omega_V = 1;
- `g2_propagation/` — the G_2 propagation identities;
- `g2_pure_parameters/` — the pure-parameter identities and auxiliary checks;
- `diagonal_braiding_theta/` — the cocycle identity relating `Theta_Phi` to
  the diagonal braiding coefficients;
- `theta_symmetry/` — opposite symmetry of Theta_Phi.

From the repository root, run `./verify_all.sh` to verify this group together
with the other groups.
