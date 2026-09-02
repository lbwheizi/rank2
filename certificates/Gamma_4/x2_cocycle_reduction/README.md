# Gamma4 X2 cocycle-reduction certificate

`gamma4_x2_cocycle_reduction_certificate.g` is a core-only GAP 4.16
source-to-target verifier for
`lem:Gamma4-X2-cocycle-reduction` in Appendix C of the paper.

The program first constructs the Laurent exponent vectors of

\[
\alpha=\frac{\rho(\varepsilon^{-1})}
              {\Phi_h(g,\varepsilon^{-1})},\qquad
\kappa_X=\frac{\rho(\varepsilon^2)}
 {\Phi_h(\varepsilon^2,h)\Phi_h(g,\varepsilon^2h)}
\]

and of the displayed definition of \(\Xi_X\).  Thus its source is the
quotient \(\Xi_X/\kappa_X\), not the cocycle-defect expression that is to
be proved.  It then applies, as two separately constructed projective
defects, precisely

\[
\rho(\varepsilon^{-1})^2
=\Phi_h(\varepsilon^{-1},\varepsilon^{-1})\rho(\varepsilon^2),
\qquad
\rho(\varepsilon^{-1})\rho(h)
=\Phi_h(\varepsilon^{-1},h)\rho(\varepsilon^{-1}h).
\]

The result is compared atom by atom with the intermediate quotient printed
in the proof.  In particular, the four projective-character atoms in the
source must cancel before the target cocycle defects are constructed.

Next the program treats the three instances of local-cocycle coherence with

\[
(x,y,z,w)=(h,g,\varepsilon^{-1},h),\quad
(h,g,\varepsilon^2,h),\quad
(h,g,\varepsilon^{-1},\varepsilon^{-1}).
\]

It does not merely declare these three quotients equal to one.  For each
instance it expands \(\Phi_x\) into values of \(\Phi\) and independently
checks the four-cocycle-defect certificate obtained from the proof of
`lem:local-cocycle-coherence`.  It then applies the three identities in the
orientations used in Appendix C, expands every remaining \(\Phi_h\) and
\(\Phi^h\), and compares the result with

\[
\frac{
 \mathscr C_\Phi(h,g,h,\varepsilon^{-1})
 \mathscr C_\Phi(\varepsilon^{-1}h,\varepsilon g,
                  \varepsilon^{-1},h)
 \mathscr C_\Phi(\varepsilon^{-1}h,h,g,\varepsilon^{-1})
}{
 \mathscr C_\Phi(h,g,\varepsilon^{-1},h)
 \mathscr C_\Phi(\varepsilon^{-1}h,\varepsilon g,
                  h,\varepsilon^{-1})
 \mathscr C_\Phi(h,\varepsilon^{-1}h,g,\varepsilon^{-1})
}.
\]

The final equality is checked in an exact sparse Laurent exponent lattice.
Group products and conjugations are computed in the abstract normal form
\(\varepsilon^i h^j g^k\), using only the presentation of \(\Gamma_4\).
No finite quotient, sampled cocycle, floating-point arithmetic, or external
GAP package is used.

The certificate verifies only this formal source-to-target cocycle
reduction.  It does not verify the construction or simplicity of \(X_2\),
the finite-dimensionality criterion, or the classification theorem.

Run from this directory with GAP 4.16:

```text
gap --bare -A -r -q --quitonbreak --norepl gamma4_x2_cocycle_reduction_certificate.g
```

For an interactive run, start `gap --bare -A -r -q` and then use

```gap
Read("gamma4_x2_cocycle_reduction_certificate.g");
```

Success returns to `gap>`.  Failure calls `ErrorNoReturn` before the final
`PASS` line and stops at `brk>`; use `quit;` to return to the outer prompt.
The script contains no `QUIT_GAP` call.

The checked standard output is stored in
`gamma4_x2_cocycle_reduction_certificate.out`.
`CHECKSUMS.sha256` records the delivered files.
