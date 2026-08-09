# Gamma4 terminal-character certificate

`gamma4_terminal_character_certificate.g` is a core-only GAP 4.16 verifier
for the terminal-root self-character cancellation used in Remark C.12 of the
paper.  In the abstract normal form \(\varepsilon^i h^j g^k\), it reconstructs
the 76-term Laurent target and verifies its expression as an integral product
of 146 normalized 3-cocycle defects.  The coefficient \(\ell^1\)-norm is 245.

The program uses no external GAP package, finite quotient, sampled cocycle,
or floating-point arithmetic.  Its scope is the final cocycle quotient only;
it does not verify the preceding action formulas or the reflection argument
in the main text.

Run from this directory with GAP 4.16:

```text
gap -A -q gamma4_terminal_character_certificate.g
```

The checked standard output is stored in
`gamma4_terminal_character_certificate.out`.  `CHECKSUMS.sha256` records the
delivered files.
