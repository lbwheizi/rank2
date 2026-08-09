# Auxiliary GAP certificate: `s_phi_r2_invariance`

This self-contained, core-only script verifies

$$
s_\Phi(hg,h^{-1})=s_\Phi(g,h).
$$

This is an auxiliary consistency check; the parameter classification uses
the stronger identity $s_\Phi=\Theta_\Phi^2$ together with the corresponding
reflection identity for $\Theta_\Phi$.  It is not cited as a separate lemma
in the paper.

Group elements are represented by `[e,a,b]`, meaning
$\varepsilon^e g^a h^b$, and the calculation uses only
$hg=\varepsilon gh$, $\varepsilon^2=1$, and centrality of $\varepsilon$.
The target is reduced exactly by an integral combination of normalized
$3$-cocycle relators.  No finite quotient, rationalization, floating-point
arithmetic, or optional GAP package is used.  The certificate has 101
relators and coefficient $\ell^1$-norm 175.

## Run

From this directory:

```sh
gap --bare -A -r -q --nointeract auxiliary_s_phi_R2_invariance.g
```

The checked output is in `auxiliary_s_phi_R2_invariance.out`; file hashes are
in `CHECKSUMS.sha256`.
