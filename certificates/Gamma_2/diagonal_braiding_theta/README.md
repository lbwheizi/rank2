# GAP certificate for `lem:Theta-diagonal-braiding-cocycle`

This directory contains the exact GAP certificate for the identity

```text
Theta_Phi(x,y) Phi_x(epsilon,epsilon) Phi_x(x,y)
------------------------------------------------ = 1,
       Phi_x(epsilon,x)^2 Phi_x(y,epsilon*x)
```

where `y*x=epsilon*x*y`, `epsilon^2=1`, and `epsilon` is central.  The
program checks the calculation for the abstract generators `g,h`; the same
calculation applies to every pair `x,y` satisfying these relations.

The script expands the definition of `Theta_Phi`, every local cocycle
`Phi_x`, and every tensor-action factor `Phi^x` into an integer exponent
vector in the normalized symbols `Phi(a,b,c)`.  Independently, it constructs
the product of the twenty normalized 3-cocycle identities displayed in the
paper.  It then checks that the two exponent vectors agree exactly.

The computation uses neither a finite quotient nor sampled cocycles.  It uses
only exact integer arithmetic and the GAP core library.

## Run

With GAP 4.16.0, run from this directory:

```sh
gap --bare -A -r -q --nointeract theta_diagonal_braiding_cocycle.g
```

The expected output is stored in `theta_diagonal_braiding_cocycle.out`.
File hashes are stored in `CHECKSUMS.sha256`.
