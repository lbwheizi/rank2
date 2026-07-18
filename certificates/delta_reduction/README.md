# GAP certificate for the first-adjoint $\Delta$ reduction

The script `delta_certificate.g` verifies the formal $3$-cocycle reduction
used in the simplicity calculation for the first adjoint object $X_1$ in the
$\Gamma_2$ setting. Group elements are written in the normal form
$\varepsilon^e g^a h^b$, subject to

$$
h g=\varepsilon g h,\qquad \varepsilon^2=1,\qquad
\varepsilon\in Z(G).
$$

With the left-left local cocycle conventions used in the accompanying paper,
the script verifies

$$
\frac{
a^2A_1B_1
}{
\lambda'\Phi^h(g,h)\Phi^g(g,h)
}
=\Delta,
$$

where

$$
a=\frac{\zeta}{\Phi_g(h,\varepsilon)},
$$

$$
A_1=
\Phi^h(\varepsilon g,\varepsilon h)
\Phi_g(h,h)\mu
\frac{
\Phi_h(h,g)\lambda'\zeta'
}{
\Phi_h(\varepsilon,h)\Phi_h(g,\varepsilon h)
},
$$

$$
B_1=
\Phi^g(\varepsilon g,\varepsilon h)
\frac{
\Phi_g(g,h)\zeta
}{
\Phi_g(\varepsilon,g)\Phi_g(h,\varepsilon g)
}
\Phi_h(g,g)\mu',
$$

and

$$
\Delta=
\frac{
\zeta\zeta'\mu\mu'
\Phi_h(\varepsilon g,g)\Phi_g(h,h)
}{
\Phi_g(h,\varepsilon)\Phi_h(\varepsilon,g^2)
}.
$$

## Method

After cancelling the one-dimensional projective-character parameters, the
only remaining character factor is replaced using

$$
\zeta^2=\Phi_g(\varepsilon,\varepsilon).
$$

A formal product of values $\Phi(a,b,c)$ is stored as an exponent vector in
the free abelian group on normalized triples. The script expands the quotient
of the first-adjoint coefficient by $\Delta$ and an explicit product of
fifteen normalized $3$-cocycle relators, then checks that their exponent
vectors agree exactly.

The script does not assert that $\Delta$ is identically equal to $1$. The
condition $\Delta=1$ arises from the proportionality condition characterizing
the simplicity of $X_1$.

## Run

From the repository root:

```sh
gap -A -q certificates/delta_reduction/delta_certificate.g
```

Expected output:

```text
Formal Phi terms in the reduced Delta quotient: 27
Formal Phi terms remaining after subtracting the certificate: 0
SUCCESS: the cocycle certificate for the Delta reduction is correct.
Hence the first-adjoint coefficient equals Delta.
```

The saved output of the verified run is in `delta_certificate.out`.

## SHA-256

```text
dcb1fba0925c764d7887f55e9072d3f22f4832e11c8f0e37d22bead7f2c05b3a  delta_certificate.g
3621e395782ac0a83e94d2555818d2eb9f78033940617dc4ee2123d1a0007443  delta_certificate.out
```
