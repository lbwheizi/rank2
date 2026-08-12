# GAP certificate for Lemma B.12

This program corresponds to Lemma B.12 of the paper, with LaTeX label
`lem:Gamma2-R2-Theta-cocycle`.

This self-contained, core-only script verifies

$$
\Theta_\Phi(hg,h^{-1})=\Theta_\Phi(g,h).
$$

Group elements are represented by `[e,a,b]`, meaning
$\varepsilon^e g^a h^b$, and the calculation uses only

$$
hg=\varepsilon gh,
\qquad \varepsilon^2=1,
\qquad \varepsilon\in Z(G).
$$

Every local cocycle and tensor-action factor is expanded as an integer
exponent vector in normalized symbols $\Phi(a,b,c)$.  The script subtracts
the stated product of normalized $3$-cocycle relators and checks that the
integral residual is zero.  It uses no finite quotient, rationalization,
floating-point arithmetic, or optional GAP package.  The certificate has 94
relators and coefficient $\ell^1$-norm 131.

## Run

From this directory:

```sh
gap --bare -A -r -q --nointeract lemma_B15_theta_after_second_reflection.g
```

The checked output is in `lemma_B15_theta_after_second_reflection.out`; file
hashes are in `CHECKSUMS.sha256`.
