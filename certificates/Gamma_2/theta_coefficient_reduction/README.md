# GAP certificate for the two $\Theta_\Phi$ coefficient reductions

This directory corresponds to the two lemmas with LaTeX labels
`lem:Theta-Phi-reduction` and `lem:opposite-Theta-reduction` in Appendix B.
The program `theta_coefficient_reductions.g` verifies the coefficient
reductions used in the simplicity criteria for $X_2$ and $Y_2$.

For the forward direction, the program constructs the left-hand side of

$$
\frac{
\Phi^h(g,gh)A_1
}{
p^2\Phi^g(g,h)\Phi^h(\varepsilon g,\varepsilon gh)
\Phi_g(h,h)\mu B_1D_1
}
=\Theta_\Phi(g,h)
$$

from the displayed definitions of $A_1$, $B_1$, $D_1$, and
$p=\Phi_g(h,g)^{-1}$.  The target $\Theta_\Phi(g,h)$ is constructed
separately from its displayed definition.  Before imposing any assumption,
the program verifies the exact Laurent identity

$$
\frac{\text{left-hand side}}{\Theta_\Phi(g,h)}
=
\Delta^{-1}
\frac{{\zeta'}^2}{\Phi_h(\varepsilon,\varepsilon)}.
$$

It follows that $\Delta=1$ and
${\zeta'}^2=\Phi_h(\varepsilon,\varepsilon)$ give the asserted reduction.

For the opposite direction, the program independently constructs the
primed definitions $A_1'$, $B_1'$, $D_1'$, and
$p'=\Phi_h(g,h)^{-1}$; it does not obtain them by exchanging variables in
the forward expression.  It verifies

$$
\frac{
\displaystyle
\frac{\Phi^g(h,hg)A_1'}{{p'}^2\Phi^h(h,g)
\Phi^g(\varepsilon h,\varepsilon hg)
\Phi_h(g,g)\mu'B_1'D_1'}
}{\Theta_\Phi(h,g)}
=
{\Delta'}^{-1}
\frac{\zeta^2}{\Phi_g(\varepsilon,\varepsilon)}.
$$

Thus $\Delta'=1$ and
$\zeta^2=\Phi_g(\varepsilon,\varepsilon)$ give the opposite reduction.

## Method and scope

Group elements are stored in the normal form $\varepsilon^e g^a h^b$ and
are multiplied using

$$
hg=\varepsilon gh,\qquad \varepsilon^2=1,
\qquad \varepsilon\in Z(G).
$$

Each character scalar, each local cocycle value $\Phi_x(a,b)$, and each
tensor-action scalar $\Phi^x(a,b)$ is treated as an independent Laurent
generator.  A product is stored as an integral exponent vector.  Hence both
checks are exact: they use neither numerical cocycle samples nor
floating-point arithmetic.  No $3$-cocycle identity is required for these
two substitutions.

As a non-vacuity check, the program deliberately changes one Laurent
exponent in each direction and requires both perturbed comparisons to have a
nonzero residual.

The program does not assume that either coefficient reduction is true.  It
first forms each quotient from the original definitions, then checks that
the quotient differs from the corresponding $\Theta_\Phi$ expression by
exactly the displayed consequence of $\Delta=1$ or $\Delta'=1$ and the
relevant projective-action identity.

## Run

From this directory:

```sh
gap --bare -A -r -q --quitonbreak --norepl theta_coefficient_reductions.g
```

For an interactive run, start `gap --bare -A -r -q` and load the file with
`Read("theta_coefficient_reductions.g");`.  The program contains no
`QUIT_GAP`: success returns to `gap>`, while failure enters `brk>` through
`ErrorNoReturn` and cannot reach the final `PASS` message.

The checked output is archived in `theta_coefficient_reductions.out`. File
hashes are recorded in `CHECKSUMS.sha256`.
