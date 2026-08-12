# GAP certificate for Lemma B.15

This program corresponds to Lemma B.15 of the paper, with LaTeX label
`lem:Theta-opposite-symmetry`.

This directory gives a self-contained, core-only GAP verification of

$$
\Theta_\Phi(h,g)=\Theta_\Phi(g,h).
$$

Group elements are represented in the abstract normal form
$\varepsilon^e g^a h^b$ by `[e,a,b]`, using only

$$
hg=\varepsilon gh,\qquad \varepsilon^2=1,\qquad
\varepsilon\in Z(G).
$$

The script expands every local cocycle $\Phi_x$, tensor-action factor
$\Phi^x$, and the target scalar into an integer exponent vector in the free
abelian group on normalized values $\Phi(a,b,c)$. It then subtracts the
listed product of normalized $3$-cocycle relators and checks that the exact
integer residual is zero. No rationalization, floating-point arithmetic, or
optional GAP package is used.

The certificate contains **14** relators and has coefficient
$\ell^1$-norm **14**.

## Run

From this directory:

```sh
gap --bare -A -r -q --nointeract lemma_B19_theta_symmetry.g
```

The checked output is saved in `lemma_B19_theta_symmetry.out`. File hashes
are in `CHECKSUMS.sha256`.


## Reduced target used by the certificate

The paper first proves

$$
\Theta_\Phi(a,b)=
\frac{\Phi_a(\varepsilon,a)^2\Phi_a(b,\varepsilon a)}
{\Phi_a(\varepsilon,\varepsilon)\Phi_a(a,b)}.
$$

The 14-relator certificate checks the quotient of this reduced expression
for $(a,b)=(g,h)$ and $(h,g)$. It therefore does not duplicate the earlier
reduction from the long definition of $\Theta_\Phi$.
