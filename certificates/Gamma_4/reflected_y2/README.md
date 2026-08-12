# Gamma4 reflected-Y2 certificate

This directory directly verifies the second-adjoint simplicity statement
after the first reflection used in Lemma C.4
(`lem:Gamma4-R1-reflected-parameter-reduction`):

\[
R_1(V,W)=(V^*,X_1),\qquad
(\varepsilon',h',g')=(\varepsilon^{-1},h^{-1},hg).
\]

The seven reflected line-stability scalars compare the base coordinate `5`
with the coordinate indices `6, 8, 15, 17, 18, 27, 28`.  The verifier proves
that each is one from the initial scalar rows.  It first reconstructs their
source and then performs two exact integral certificate stages.

0. `gamma4_reflected_y2_source.g` independently builds, from the abstract
   Gamma4 group law and the displayed transversals, the induced actions on
   the initial objects (V,W).  It constructs (X_1) by `phi_1`, extracts
   its transported projective character, constructs the rigid dual (V^*),
   and then builds the reflected induced pair

   \[
   (\varepsilon',h',g')=(\varepsilon^{-1},h^{-1},hg).
   \]

   For both the initial and reflected pairs it reconstructs `phi_1`, the
   corrected `phi_2`, and all eight nonzero coordinates.  The supports are
   checked to be exactly

   ```text
   5, 6, 8, 15, 17, 18, 27, 28
   ```

   and the homogeneous terminal degrees are respectively
   `epsilon*h*g^2` and `epsilon^2*h*g^2`.  The seven source ratios are then
   canonically expanded using only the one-dimensional projective-character
   law.  GAP checks them atom by atom against the thirteen initial scalar
   rows and the seven reflected target rows in the data file.  Thus those
   expanded rows are not accepted as unverified source input.

1. The data file contains the 13 fully expanded initial scalar rows and the
   seven fully expanded reflected target rows.  These are the Laurent
   monomial rows obtained directly from the recursive `phi_2` construction,
   with the corrected right-associated convention

   ```text
   c12 = a (c_WW tensor id) a^-1,
   coefficient = Phi(x,y,l) / Phi(x*y*x^-1,x,l).
   ```

   After the independent source check, the GAP script verifies coefficient
   by coefficient the seven identities

   \[
   T_j\prod_{i=0}^{12}P_i^{C_{ji}}=B_j(\Phi),
   \]

   where the order of the initial scalar rows is

   ```text
   R(h)/(-1), S(g)/(-1), Delta4, Delta4_reflected,
   self_Vdual, self_X1, initial_E0, ..., initial_E6.
   ```

   The integer matrix `C` is embedded in the data file.  No rational powers
   are used.  The seven pure-cocycle residual rows have respectively

   ```text
   92, 114, 152, 99, 149, 145, 142
   ```

   nonzero Laurent atoms.

2. For every residual `B_j(Phi)`, the script checks an explicit integral
   normalized bar-boundary certificate.  The seven certificates contain

   ```text
   499, 1085, 1173, 813, 1121, 1137, 1079
   ```

   nonzero 4-simplices.  Their coefficient `l1` norms are

   ```text
   571, 1800, 1952, 1296, 2003, 1862, 1738.
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
gap -A -q gamma4_reflected_y2_certificate.g
```

The expected standard output is stored in
`gamma4_reflected_y2_certificate.out`.  `CHECKSUMS.sha256` records the
delivered files.
