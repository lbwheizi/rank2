# rank2

Computer-verifiable certificates and GAP programs accompanying the study of
finite-dimensional rank-two Nichols algebras in
$\mathcal Z(\operatorname{Vec}_G^\Phi)\simeq {}_G^G\mathcal{YD}^{\Phi}$.

## Certificates

- [$\Omega_V=1$](certificates/omegaV/): a formal normalized $3$-cocycle
  certificate for the scalar appearing in the computation of the second
  adjoint object in the $\Gamma_2$ setting.

## Requirements

- GAP 4.16.0 (the version used for the recorded run)

No additional GAP packages are required by the current certificate.

## Run the current certificate

From the repository root, run

```sh
gap -A -q certificates/omegaV/omegaV_certificate.g
```

Expected output:

```text
Formal Phi terms in the expanded Omega_V: 18
Formal Phi terms remaining after subtracting the certificate: 0
SUCCESS: the cocycle certificate for Omega_V is correct.
Hence Omega_V is a product of cocycle relators and equals 1.
```

## Scope

The script represents every normalized value $\Phi(a,b,c)$ by a formal
generator, expands the local cocycle expressions, and checks coefficientwise
that $\Omega_V$ is the displayed product of normalized $3$-cocycle relators.
It certifies this cocycle identity only; the surrounding mathematical argument
remains in the accompanying paper.

