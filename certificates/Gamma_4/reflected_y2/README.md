# Gamma4 reflected-Y2 parameter certificate

This directory verifies the coordinate calculation in Lemma C.3 and the
seven parameter identities after the first reflection in Lemma C.7 and
equation (C.11) (`lem:Gamma4-R1-reflected-parameter-reduction`):

\[
R_1(V,W)=(V^*,X_1),\qquad
(\varepsilon',h',g')=(\varepsilon^{-1},h^{-1},hg).
\]

Write the eight nonzero coordinates as `E_0,...,E_7`, with coefficients
`calA_0,...,calA_7`.  They are numbered so that conjugation by `h` acts as

```text
(0 1 2 3)(4 5 6 7).
```

If `h.E_j=H_j E_pi(j)`, the seven scalars are

```text
Theta_2_j = H_j calA_j calA_1 / (H_0 calA_0 calA_pi(j)),  1 <= j <= 7.
```

The verifier reconstructs the seven scalars from the initial data and then
performs two exact integral certificate checks.  Together with the thirteen
identities \(\mathcal P_i=1\) proved in the manuscript, these calculations give
\(\Theta^{(1)}_{2,j}=1\) for \(1\leq j\leq7\).  They do not by themselves
prove that the reflected second adjoint object is simple; that conclusion
also uses the mathematical argument in the manuscript.

0. `gamma4_reflected_y2_source.g` builds, from the abstract Gamma4 group
   law and the displayed transversals, the induced actions on the initial
   objects (V,W).  It constructs the displayed first-root vector by
   `phi_1`, reads the scalar action on its line from one nonzero coordinate,
   constructs the rigid dual (V^*), and then builds the reflected induced
   pair

   \[
   (\varepsilon',h',g')=(\varepsilon^{-1},h^{-1},hg).
   \]

   The stability, simplicity, and one-dimensional homogeneous component of
   the first-root object are inputs from
   `lem:Gamma4-first-adjoint` in the manuscript.  In particular, the source
   program does not independently compare both coordinates of the
   first-root vector.  For both the initial and reflected pairs it
   reconstructs `phi_1`, the
   corrected `phi_2`, and all eight nonzero coordinates.  It checks directly
   from the action matrix that the coordinates form the two cycles displayed
   above.  The homogeneous terminal degrees are respectively
   `epsilon*h*g^2` and `epsilon^2*h*g^2`.  The seven scalars in the displayed
   formula are then expanded using only the one-dimensional
   projective-character law.  It also reconstructs the first six packet
   rows

   ```text
   rho(h)/(-1), sigma(g)/(-1), Delta4, Delta4_reflected,
   rho1(h1)/(-1), sigma1(g1)/(-1).
   ```

   GAP checks these six rows, the seven initial Theta rows, and the seven
   reflected target rows atom by atom.  Thus all twenty expanded rows are
   compared with a separate source calculation, subject to the preceding
   first-root-line input from the manuscript.

1. The generated data first records seven base-coordinate ratios
   `R_j=L_0/L_j`.  The verifier applies the exact integral identities

   ```text
   Theta_2_1 = R_1/R_2,  Theta_2_2 = R_1/R_3,
   Theta_2_3 = R_1,      Theta_2_4 = R_1/R_5,
   Theta_2_5 = R_1/R_6,  Theta_2_6 = R_1/R_7,
   Theta_2_7 = R_1/R_4.
   ```

   It applies the same change of basis to the initial rows, reflected target
   rows, pure-cocycle residuals, and bar certificates.  The Stage-A matrix is
   transformed on both its rows and its last seven columns.  All these rows
   come from the recursive `phi_2` construction with the corrected
   right-associated convention

   ```text
   c12 = a (c_WW tensor id) a^-1,
   coefficient = Phi(x,y,l) / Phi(x*y*x^-1,x,l).
   ```

   After the independent source check, the GAP script verifies coefficient
   by coefficient the seven identities

   \[
   \Theta^{(1)}_{2,j}
   \prod_{i=1}^{13}\mathcal P_i^{c_{ji}}
   =\mathcal R_j(\Phi),
   \]

   where the order of the initial scalar rows is

   ```text
   rho(h)/(-1), sigma(g)/(-1), Delta4, Delta4_reflected,
   rho1(h1)/(-1), sigma1(g1)/(-1),
   initial_Theta_2_1, ..., initial_Theta_2_7.
   ```

   The displayed Stage-A matrix has zero columns at `Delta4` and
   `rho1(h1)/(-1)`; the verifier nevertheless retains, reconstructs, and
   audits both packet rows and checks explicitly that those columns vanish.

   The resulting integer matrix `C` is embedded in the verifier.  No rational
   powers are used.  The seven pure-cocycle residual rows have respectively

   ```text
   136, 24, 152, 86, 141, 49, 156
   ```

   nonzero Laurent atoms.

2. For every residual `B_j(Phi)`, the script checks an explicit integral
   normalized bar-boundary certificate.  The seven certificates contain

   ```text
   1561, 276, 1173, 1623, 1095, 516, 1306
   ```

   nonzero 4-simplices.  Their coefficient `l1` norms are

   ```text
   2706, 278, 1952, 2802, 1599, 540, 1999.
   ```

Group elements are represented throughout by abstract normal forms
`epsilon^i h^j g^k`, with only `i` reduced modulo four.  The exponents `j`
and `k` are genuine integers.  Thus the check uses `epsilon^4=1` but no
finite-order relation for `h` or `g`.  The verifier also reports the maximum
absolute `h,g` exponents in certificate entries and adjacent products.

The calculation uses exact integer Laurent arithmetic, core GAP only, and no
floating-point computation, external package, finite cocycle sample, or
finite Gamma4 quotient.  The source engine and the certificate verifier are
both part of the checksum manifest.

Run from this directory with GAP 4.16:

```text
gap --bare -A -r -q --nointeract gamma4_reflected_y2_certificate.g
```

The expected standard output is stored in
`gamma4_reflected_y2_certificate.out`.  `CHECKSUMS.sha256` records the
delivered files.
