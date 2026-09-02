# Gamma4 Xi3 cocycle-reduction certificate

`gamma4_xi3_certificate.g` is a core-only GAP 4.16 source-to-target
verifier for the cocycle reduction and coefficient cancellation recorded
in Lemmas C.4 and C.5 and equation (C.9) of the paper
(`lem:Gamma4-Xi3-cocycle-reduction`,
`lem:Gamma4-Y3-coefficient-cancellation`, and
`eq:Gamma4-Xi3-direct`).
It first reconstructs the corrected recursive maps
\(\varphi_1^\Phi,\varphi_2^\Phi,\varphi_3^\Phi\), including
\(c_{12}=a(c_{W,W}\otimes\mathrm{id})a^{-1}\), constructs the eight
nonzero coordinates of the chosen \(Y_2\)-seed, and extracts the two
Laurent monomials in the common zero-based coordinate 45 of \(Y_3\).
The program computes the corresponding GAP position 46 from the four tensor
degrees

\[
(\varepsilon g,\varepsilon g,\varepsilon^2g,\varepsilon^{-1}h)
\]

and checks that their product is \(\varepsilon^2hg^3\).
After the displayed one-dimensional projective-character reduction, it
checks that their signed ratio is exactly the eight-atom target.  With

\[
x=\varepsilon^2g,\qquad b=\varepsilon g,
\]

it checks, in an exact sparse Laurent exponent lattice, that the remaining
cocycle quotient in the second coefficient of the third opposite adjoint
object is

\[
\frac{\mathscr C_\Phi(h,g,x,g)}
     {\mathscr C_\Phi(b,h,x,g)}=1.
\]

It then verifies that target by the two normalized 3-cocycle defects shown
above.  All group elements are handled as abstract normal forms
\(\varepsilon^i h^j g^k\); no quotient-specific order relation is used.
The program therefore checks the full source-to-target bridge for this
chosen coefficient, but it does not verify the simplicity criterion for
\(Y_2\), the other coefficient of \(Y_3\), or any assertion about the
Cartan graph.

Run from this directory with GAP 4.16:

```text
gap --bare -A -r -q --quitonbreak --norepl gamma4_xi3_certificate.g
```

For an interactive run, first start `gap --bare -A -r -q`, then load the
script at `gap>` with `Read("gamma4_xi3_certificate.g");`.  Success returns to
`gap>`; failure
stops at `brk>`, where `quit;` returns to the outer `gap>` prompt.

The checked standard output is stored in
`gamma4_xi3_certificate.out`.  `CHECKSUMS.sha256` records the delivered
files.
