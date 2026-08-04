# GAP certificates for the $G_2$ propagation formulas

This core-only GAP script verifies the two pure normalized-$3$-cocycle
reductions

$$
\mathcal K_3=1,
\qquad
d_4=\lambda^3\kappa_4.
$$

The targets are reconstructed directly from the displayed definitions of
$\mathcal K_3$, $d_4$, $\kappa_3$, and $\kappa_4$.  In the second target the
script uses only the projective-character identity
$\zeta^2=\Phi_g(\varepsilon,\varepsilon)$, exactly as in the paper.

The rebuilt exact certificates are shorter than the lost earlier lists:

- $\mathcal K_3=1$: 13 relators,
  coefficient $\ell^1$-norm 13;
- $d_4=\lambda^3\kappa_4$: 18 relators,
  coefficient $\ell^1$-norm 22.

The previous workspace's 36- and 38-relator lists were not present after
the reset.  These new certificates verify the same two Laurent exponent
targets exactly and do not depend on those missing lists.

## Run

```sh
../../../../tools/bin/gap -A -q --nointeract certificate.g
```

No optional GAP package, floating-point calculation, or finite quotient is
used.  The actual checked output is in `recorded_output.txt`; hashes are in
`CHECKSUMS.sha256`.
