# Gamma4 Xi3 cocycle-reduction certificate

`gamma4_xi3_certificate.g` is a core-only GAP 4.16.0 source-to-target
verifier for the cocycle reduction and coefficient cancellation recorded
in Lemmas C.4 and C.5 and `eq:Gamma4-y3-two-coefficients` of the paper
(`lem:Gamma4-Xi3-cocycle-reduction`,
`lem:Gamma4-Y3-coefficient-cancellation`).
It first reconstructs the corrected recursive maps
\(\varphi_1^\Phi,\varphi_2^\Phi,\varphi_3^\Phi\), including
\(c_{12}=a(c_{W,W}\otimes\mathrm{id})a^{-1}\), constructs the eight
nonzero coordinates of the chosen \(Y_2\)-seed, and verifies both
coefficients in `eq:Gamma4-y3-two-coefficients`.

For the first coefficient, it recovers \(\beta\) from the recursive
\(y_1\)-vector, extracts \(B_1\) from the actual `c12` image, obtains
\(\nu_1\) from the two action columns defining
\(r_1=(\varepsilon^2g)\rhd w\) and \(g\rhd r_1\), and obtains
\(\mu_1\) from a common nonzero coordinate of
\(\varphi_2^\Phi(w\otimes p_1)\) and \(z_1\).  Before using
\(\sigma(g)=-1\), the source calculation gives the stronger exact identity

\[
\beta B_1\nu_1\mu_1=-\sigma(g).
\]

It also constructs the two comparison paths from the actual
\(\varepsilon\)-action on \(w_2\otimes y_1\) and on \(z_1\); both path
residuals are checked to be zero on the common nonzero coordinate, and the
full coordinate supports are checked to agree.  Thus the scalar calculation
in the two comparison identities is a verified output rather than an
assumption.  The one-dimensionality of the relevant homogeneous component,
and hence the fact that the same scalar applies on its other coordinates,
is the mathematical input from the manuscript; the free Laurent source
engine does not independently reprove that consequence of 3-cocycle
coherence.  Substitution of \(\sigma(g)=-1\) then gives
\(\beta B_1\nu_1\mu_1=1\).  As a separate check, the full `phi3` output is
inspected on the corresponding \(w\otimes z_1\) coordinate: its two source
terms have quotient \(\sigma(g)\) and cancel after the same specialization.

For the second coefficient, the verifier extracts the two
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
pair of coefficients, but it does not verify the simplicity criterion for
\(Y_2\), the one-dimensionality argument used to define \(\mu_1\) and
\(\tau_2(\varepsilon)\), or any assertion about the Cartan graph.

Run from this directory with GAP 4.16.0:

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
