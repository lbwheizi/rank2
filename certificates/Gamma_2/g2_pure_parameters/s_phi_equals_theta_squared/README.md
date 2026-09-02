# GAP certificate for Lemma B.11

This program corresponds to Lemma B.11 of the paper, with LaTeX label
`lem:Gamma2-X3-cocycle-factor-identity`.

This self-contained, core-only script verifies

$$
s_\Phi(g,h)=\Theta_\Phi(g,h)^2.
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
floating-point arithmetic, or optional GAP package.  The certificate has 63
relator records, including records with nonunit powers, and total coefficient
$\ell^1$-norm 69.

No $G_2$ parameter condition is used: the identity holds in the abstract
$\Gamma_2$ setting described above.  It is then applied to the $G_2$
parameter locus in the paper.

Consequently, on the $G_2$ parameter locus where $\Theta_\Phi=-1$, the
condition $s_\Phi=1$ follows automatically.

## Run

From this directory:

```sh
gap --bare -A -r -q --quitonbreak --norepl lemma_B14_s_phi_equals_theta_squared.g
```

For an interactive run, first start `gap --bare -A -r -q`, then load the
script at `gap>` with `Read("lemma_B14_s_phi_equals_theta_squared.g");`.
Success returns to `gap>`; failure
stops at `brk>`, where `quit;` returns to the outer `gap>` prompt.

The checked output is in `lemma_B14_s_phi_equals_theta_squared.out`; file
hashes are in `CHECKSUMS.sha256`.
