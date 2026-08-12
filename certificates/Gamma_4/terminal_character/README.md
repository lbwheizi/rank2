# Gamma4 terminal-character certificate

`gamma4_terminal_character_certificate.g` is a core-only GAP 4.16 verifier
for the terminal-root self-character cancellation used in Lemma C.9
(`lem:Gamma4-Y2-support-value`) of the paper.  In the abstract normal form
\(\varepsilon^i h^j g^k\), it first
reconstructs the six displayed action coefficients
\(L_{qa},L_{qp},L_{ha},L_{hy},L_{hr},L_{hY}\), then builds
\(K_q,\lambda_1,\lambda_2\) and the source quotient
\(K_q\lambda_2/(\rho(h)\sigma(g)^2\lambda_1)\).  It performs the stated
projective-character eliminations using
\(\rho(h)=\sigma(g)=-1\), \(\Delta_4=1\), and \(\varepsilon^4=1\), and
checks that the result is exactly the recorded 76-atom Laurent target.
Finally it verifies that target as an integral product of 146 normalized
3-cocycle defects.  The certificate coefficient \(\ell^1\)-norm is 245.

The program uses no external GAP package, finite quotient, sampled cocycle,
or floating-point arithmetic.  Its scope now includes the algebraic
source-to-target preprocessing of the six action formulas; it does not
verify the geometric choice of homogeneous vectors or the reflection
argument in the main text.

Run from this directory with GAP 4.16:

```text
gap -A -q gamma4_terminal_character_certificate.g
```

The checked standard output is stored in
`gamma4_terminal_character_certificate.out`.  `CHECKSUMS.sha256` records the
delivered files.
