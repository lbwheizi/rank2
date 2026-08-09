# Gamma4 Xi3 cocycle-reduction certificate

`gamma4_xi3_certificate.g` is a core-only GAP 4.16 verifier for the
two-defect factorization used in Remark C.5 of the paper.  With

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

The calculation uses only the defining \(\Gamma_4\) relations and normalized
3-cocycle defects.  It does not verify the simplicity criterion for \(Y_2\),
the full coefficient calculation for \(Y_3\), or any assertion about the
Cartan graph.

Run from this directory with GAP 4.16:

```text
gap -A -q gamma4_xi3_certificate.g
```

The checked standard output is stored in
`gamma4_xi3_certificate.out`.  `CHECKSUMS.sha256` records the delivered
files.
