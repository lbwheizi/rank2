# Core-only verifier for the order-16 $G_2$ example

This directory contains a self-contained GAP 4.16 verifier for the finite
data displayed in the manuscript's order-16 example. It uses only the GAP
core library and exact exponents modulo four: the exponent `e` represents
the scalar $i^e$. No floating-point arithmetic is used.

## What is checked

The script exhaustively checks:

- the 16 normal forms, closure, identity, inverse formula, all $16^3$
  associativity instances, generation by $g,h$, and the displayed group
  relations;
- normalization of the stated $3$-cocycle and all $16^4=65536$ cocycle
  identities;
- the two stated centralizers, the closed formulas for
  $\Phi_g$ and $\Phi_h$, and every projective-character relation for the
  explicitly defined $\rho$ and $\sigma$;
- the six initial character parameters and the initial values
  $\Delta=\Delta'=1$ and $\Theta_\Phi=-1$;
- for each of the twelve printed states, generation of both centralizers,
  exact reconstruction of the two eight-value projective characters from
  their three printed signature values, and all projective-character
  relations;
- all 24 support transitions obtained from
  $R_1(p,q)=(p^{-1},p^3q)$ and $R_2(p,q)=(qp,q^{-1})$, together with
  involutivity and connectivity of the printed transition table;
- the uniform printed scalars
  $\lambda_j=i$, $\lambda'_j=-1$,
  $\Delta_j=\Delta'_j=1$, and $\Theta_{\Phi,j}=-1$;
- via the displayed induced-action formula, the two self-labels and the
  off-diagonal label product for both two-element support orbits at every
  state; and
- the final arithmetic $16^3 8^3=2^{21}$ from the stated three root factors
  in each orbit.

## Deliberate limits

This is a core-data verifier, not a replacement for the omitted recursive
adjoint calculation. It does **not** derive the terminal reflected-character
signatures from predecessor states, recompute recursive adjoint coefficients,
or prove adjoint simplicity/nonvanishing and terminal vanishing. The script
therefore makes none of those claims. It checks that the printed signatures
are individually valid projective characters and that the printed support
reflection graph is exact and internally closed.

## Run

From the workspace root:

```sh
tools/bin/gap -A -q --nointeract \
  work/gap_supplement/certificates/order16_G2/certificate.g
```

The checked output is in `recorded_output.txt`; SHA-256 hashes are in
`CHECKSUMS.sha256`.
