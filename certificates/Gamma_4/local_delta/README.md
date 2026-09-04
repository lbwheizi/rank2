# Gamma4 local-Delta certificate

`gamma4_local_delta_certificate.g` uses only core GAP 4.16.0 functions.  It
corresponds to Lemma C.8, `lem:Gamma4-W-local-Deltas`, and the formulas
`eq:Gamma4-local-Delta-01`, `eq:Gamma4-local-Delta-10`, and
`eq:Gamma4-local-action-coefficient-products`.  It
reconstructs the calculation of \(Y_3\) from the \(\Gamma _4\) group law,
`eq:tensor-product`, the braiding formula, and the definitions of
\(\varphi _1^\Phi,\varphi _2^\Phi,\varphi _3^\Phi\).  The program obtains
all 192 nonzero binomial coordinates of \(Y_3\).  It then compares, in the
same quotient orientation, the 19 coordinates used in the two integer
certificates.  Finally, it verifies that these coordinates imply

\[
\Delta_{01}^W=\Delta_{10}^W=1.
\]

The data file contains the expected Laurent quotient for each of the 19
coordinates and the two integer coefficient vectors.  The main program
reconstructs every listed quotient and requires exact equality in the stated
orientation; equality only after inversion is rejected.  It also checks that
each coefficient vector gives the inverse of the corresponding \(\Delta\)
and that all signs cancel.

The calculation uses

\[
c_{1,2}^{\Phi}=a\,(c\otimes\operatorname{id})\,a^{-1},
\]

so the input diagonal associator factor is \(\Phi\), and the output factor is
\(\Phi^{-1}\).

The calculation takes place in the free Laurent group on the formal symbols
`R`, `S`, `F`, and `P`.  Group elements are represented as
\(\varepsilon^i h^j g^k\), with only \(i\) reduced modulo four and
\(j,k\in\mathbb Z\).  It therefore uses neither a sampled cocycle nor a
special finite quotient, and it does not assume finite-dimensionality.

Run from this directory with GAP 4.16.0:

```text
gap --bare -A -r -q --quitonbreak --norepl gamma4_local_delta_certificate.g
```

For an interactive run, first start `gap --bare -A -r -q`, then load the
script at `gap>` with `Read("gamma4_local_delta_certificate.g");`.  Success
returns to `gap>`; failure
stops at `brk>`, where `quit;` returns to the outer `gap>` prompt.

The checked standard output is stored in
`gamma4_local_delta_certificate.out`.
It was regenerated in a second independent GAP invocation and compared
byte-for-byte.  `CHECKSUMS.sha256` records the delivered files.
