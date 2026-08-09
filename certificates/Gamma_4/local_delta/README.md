# Gamma4 local-Delta certificate

`gamma4_local_delta_certificate.g` is a core-only GAP 4.16 verifier for the
corrected \(Y_3=0\) coordinate calculation.  Starting from the abstract
\(\Gamma _4\) normal-form multiplication, it reconstructs the homogeneous
objects, twisted tensor actions, braidings, \(\varphi _1^\Phi\),
\(\varphi _2^\Phi\), \(\varphi _3^\Phi\), \(y_2\), and
\(g\rhd y_2\).  It obtains all 192 nonzero third-adjoint coordinate
binomials and independently checks the 19 rows used in the two Laurent
certificates.  It then verifies that these relations imply

\[
L_{01}(W)=L_{10}(W)=1.
\]

The accompanying data file is GAP syntax and contains the expected 19 rows
and the two integer elimination vectors.  Those rows are not merely assumed:
the main verifier reconstructs and compares them individually.  The entire
recursion uses

\[
c_{1,2}^{\Phi}=a\,(c\otimes\operatorname{id})\,a^{-1},
\]

so the input diagonal associator factor is \(\Phi\), and the output factor is
\(\Phi^{-1}\).  This corrects the reversed convention in the earlier Python
prototype.  The exact integer certificate is unchanged after that correction.

The calculation takes place in the free Laurent group on the formal symbols
`R`, `S`, `F`, and `P`.  Group elements are represented as
\(\varepsilon^i h^j g^k\), with only \(i\) reduced modulo four and
\(j,k\in\mathbb Z\).  It therefore uses neither a sampled cocycle nor a
special finite quotient, and it does not assume finite-dimensionality.

Run from this directory with GAP 4.16:

```text
gap -A -q gamma4_local_delta_certificate.g
```

The checked standard output is stored in
`gamma4_local_delta_certificate.out`.
It was regenerated in a second independent GAP invocation and compared
byte-for-byte.  `CHECKSUMS.sha256` records the delivered files.
