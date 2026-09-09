# GAP certificate for $\Omega_V=1$

The script `lemma_B4_Omega_V_equals_one.g` verifies the identity
$\Omega_V=1$ used in the proof of Lemma B.2 of the current paper,
with LaTeX label `lem:Xi-V-Omega-V`. The proof first establishes
$\Xi_V=\Omega_V/\Phi_g(h,g)$ and then uses this calculation to conclude
$\Xi_V=\Phi_g(h,g)^{-1}$. Thus the program verifies the cocycle calculation
within that proof. The historical script and output filenames are retained;
their numbering does not indicate the current lemma number.

This identity is used in the computation of the second adjoint object in
the $\Gamma_2$ setting. Group elements are written in the normal form
$\varepsilon^e g^a h^b$, subject to

$$
h g=\varepsilon g h,\qquad \varepsilon^2=1,\qquad
\varepsilon\in Z(G).
$$

With the left-left local cocycle conventions used in the accompanying paper,
the scalar checked by the script is

$$
\Omega_V=
\frac{
\zeta^2\Phi(g,\varepsilon g,\varepsilon h)
\Phi_g(g,h)\Phi_g(\varepsilon,h)
}{
\Phi(\varepsilon g,g,\varepsilon h)\Phi^g(g,h)
\Phi_g(\varepsilon,hg)\Phi_g(h,\varepsilon)^2
},
$$

where the one-dimensional projective-action relation
$\zeta^2=\Phi_g(\varepsilon,\varepsilon)$ is substituted before the formal
expansion.

## Method

A formal product of values $\Phi(a,b,c)$ is stored as an exponent vector in
the free abelian group on normalized triples. The script expands $\Omega_V$
and an explicit product of eleven normalized $3$-cocycle relators, then checks
that their exponent vectors agree exactly. It therefore does not select,
enumerate, or numerically evaluate a particular cocycle.

## Scope

The calculation is formal in the abstract $\Gamma_2$ normal form.  It does
not choose a finite quotient or numerically evaluate a cocycle, and it uses
no optional GAP package.

## Run

From this directory:

```sh
gap --bare -A -r -q --quitonbreak --norepl lemma_B4_Omega_V_equals_one.g
```

For an interactive run, first start `gap --bare -A -r -q`, then load the
script at `gap>` with `Read("lemma_B4_Omega_V_equals_one.g");`.  Success
returns to `gap>`; failure
stops at `brk>`, where `quit;` returns to the outer `gap>` prompt.

The checked output is in `lemma_B4_Omega_V_equals_one.out`; file hashes are
in `CHECKSUMS.sha256`.
