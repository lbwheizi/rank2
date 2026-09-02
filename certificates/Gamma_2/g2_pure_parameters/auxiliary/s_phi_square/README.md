# Auxiliary GAP certificate: `s_phi_square`

This self-contained, core-only script verifies

$$
s_\Phi(g,h)^2=1.
$$

This is an auxiliary consequence and is not an additional condition in the
$G_2$ parameter classification.  It is not cited as a separate lemma in the
paper.

Group elements are represented by `[e,a,b]`, meaning
$\varepsilon^e g^a h^b$, and the calculation uses only
$hg=\varepsilon gh$, $\varepsilon^2=1$, and centrality of $\varepsilon$.
The target is reduced exactly by an integral combination of normalized
$3$-cocycle relators.  No finite quotient, rationalization, floating-point
arithmetic, or optional GAP package is used.  The certificate has 37
relators and coefficient $\ell^1$-norm 71.

## Run

From this directory:

```sh
gap --bare -A -r -q --quitonbreak --norepl auxiliary_s_phi_square.g
```

For an interactive run, first start `gap --bare -A -r -q`, then load the
script at `gap>` with `Read("auxiliary_s_phi_square.g");`.  Success returns to
`gap>`; failure
stops at `brk>`, where `quit;` returns to the outer `gap>` prompt.

The checked output is in `auxiliary_s_phi_square.out`; file hashes are in
`CHECKSUMS.sha256`.
