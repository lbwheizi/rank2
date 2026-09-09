# Gamma4 Delta-transport certificate (auxiliary historical check)

`gamma4_delta_transport_certificate.g` is a core-only GAP 4.16.0 verifier
retained from an earlier version of the paper. The current manuscript no
longer states or uses the Delta-transport lemma or its Omega formulas.
Instead, `lem:Gamma4-global-to-local-Deltas` uses the first reflection and
the first-adjoint simplicity criterion to obtain the local parameter for
`Y1` from that for `W`.

This auxiliary verifier checks the three normalized bar-boundary identities

\[
\frac{\Omega_{01}(g,h)}{\Omega_{10}(g,h)}=1,
\qquad
\frac{\Omega_{01}(hg,h)}{\Omega_{01}(g,h)}=1,
\qquad
\frac{\Omega_{10}(hg,h)}{\Omega_{10}(g,h)}=1.
\]

Here \(\Omega_{01}(x,r)\) and \(\Omega_{10}(x,r)\) denote the two
Laurent monomials constructed by `OmegaPair(x,r)` in the verifier. These
symbols describe the retained calculation and are not notation used in the
current manuscript.

The script expands \(\Phi_x\) in terms of \(\Phi\) to reconstruct the
three target Laurent monomials, then checks the embedded integer combinations
of normalized \(3\)-cocycle boundaries coefficient by coefficient. The data file
contains respectively 134, 810, and 770 boundary terms.

Group elements are represented by abstract normal forms
\(\varepsilon^i h^j g^k\), where only \(i\) is reduced modulo four.  The
exponents \(j,k\) remain integers.  Consequently the verification does not
use \(h^8=1\), \(g^8=1\), or any other special finite-quotient relation.
No external GAP package, floating-point arithmetic, or cocycle sampling is
used.

Run from this directory with GAP 4.16.0:

```text
gap --bare -A -r -q --quitonbreak --norepl gamma4_delta_transport_certificate.g
```

For an interactive run, first start `gap --bare -A -r -q`, then load the
script at `gap>` with `Read("gamma4_delta_transport_certificate.g");`.
Success returns to `gap>`; failure
stops at `brk>`, where `quit;` returns to the outer `gap>` prompt.

The checked standard output is stored in
`gamma4_delta_transport_certificate.out`.
It was regenerated in a second independent GAP invocation and compared
byte-for-byte.  `CHECKSUMS.sha256` records the delivered files.
