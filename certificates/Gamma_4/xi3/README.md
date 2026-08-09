# Gamma4 Xi3 cocycle-reduction certificate

`gamma4_xi3_certificate.g` is a core-only GAP 4.16 source-to-target
verifier for the coefficient comparison used in Remark C.5 of the paper.
It first reconstructs the corrected recursive maps
\(\varphi_1^\Phi,\varphi_2^\Phi,\varphi_3^\Phi\), including
\(c_{12}=a(c_{W,W}\otimes\mathrm{id})a^{-1}\), constructs the eight
nonzero coordinates of the chosen \(Y_2\)-seed, and extracts the two
Laurent monomials in the common zero-based coordinate 45 of \(Y_3\).
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
gap -A -q gamma4_xi3_certificate.g
```

The checked standard output is stored in
`gamma4_xi3_certificate.out`.  `CHECKSUMS.sha256` records the delivered
files.
