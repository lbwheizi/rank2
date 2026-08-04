# GAP certificate: `xi_r2_invariance`

This directory gives a self-contained, core-only GAP verification of

$$
\Xi_\Phi(hg,h^{-1})=\Xi_\Phi(g,h).
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

The certificate contains **101** relators and has coefficient
$\ell^1$-norm **175**.

## Run

From this directory:

```sh
gap -A -q --nointeract certificate.g
```

With the bundled workspace installation, the equivalent command is:

```sh
../../../tools/bin/gap -A -q --nointeract certificate.g
```

The checked output is saved in `recorded_output.txt`. File hashes are in
`CHECKSUMS.sha256`.
