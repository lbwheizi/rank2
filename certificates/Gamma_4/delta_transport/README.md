# Gamma4 Delta-transport certificate

`gamma4_delta_transport_certificate.g` is a core-only GAP 4.16 verifier for
the three normalized bar-boundary identities

\[
\frac{\Omega_{01}(g,h)}{\Omega_{10}(g,h)}=1,
\qquad
\frac{\Omega_{01}(hg,h)}{\Omega_{01}(g,h)}=1,
\qquad
\frac{\Omega_{10}(hg,h)}{\Omega_{10}(g,h)}=1.
\]

The script reconstructs the three target Laurent monomials from the displayed
formula for \(\Phi_x\), then checks the embedded integer combinations of
normalized (3)-cocycle boundaries coefficient by coefficient.  The data file
contains respectively 134, 810, and 770 boundary terms.

Group elements are represented by abstract normal forms
\(\varepsilon^i h^j g^k\), where only \(i\) is reduced modulo four.  The
exponents \(j,k\) remain integers.  Consequently the verification does not
use \(h^8=1\), \(g^8=1\), or any other special finite-quotient relation.
No external GAP package, floating-point arithmetic, or cocycle sampling is
used.

Run from this directory with GAP 4.16:

```text
gap -A -q gamma4_delta_transport_certificate.g
```

The checked standard output is stored in
`gamma4_delta_transport_certificate.out`.
It was regenerated in a second independent GAP invocation and compared
byte-for-byte.  `CHECKSUMS.sha256` records the delivered files.
