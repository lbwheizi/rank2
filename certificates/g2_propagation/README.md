# GAP certificates for two $G_2$ parameter calculations

The core-only script verifies the normalized-$3$-cocycle reductions

$$
\mathcal K_3=1,
\qquad
d_4=\lambda^3\kappa_4.
$$

The targets are reconstructed from the displayed definitions of
$\mathcal K_3$, $d_4$, $\kappa_3$, and $\kappa_4$.  For the second identity,
the only projective-character relation used is
$\zeta^2=\Phi_g(\varepsilon,\varepsilon)$.

The calculation uses the abstract normal form $\varepsilon^e g^a h^b$ and
only the relations $hg=\varepsilon gh$, $\varepsilon^2=1$, and centrality of
$\varepsilon$.  It neither chooses a finite quotient nor evaluates a
particular cocycle.  The first certificate has 13 relators and coefficient
$\ell^1$-norm 13; the second has 18 relators and coefficient $\ell^1$-norm
22.

## Run

From this directory:

```sh
gap --bare -A -r -q --nointeract certificate.g
```

The checked output is in `recorded_output.txt`; file hashes are in
`CHECKSUMS.sha256`.
