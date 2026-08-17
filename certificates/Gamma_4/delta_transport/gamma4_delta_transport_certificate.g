#############################################################################
## Gamma4 Omega-transport bar-boundary certificate (core GAP only)
##
## PURPOSE
##   Expand Phi_x in terms of a normalized 3-cocycle Phi and verify exact
##   normalized bar-boundary certificates for
##
##       Omega_01(g,h) / Omega_10(g,h),
##       Omega_01(hg,h) / Omega_01(g,h),
##       Omega_10(hg,h) / Omega_10(g,h).
##
##   A zero residual proves that all four local obstruction scalars agree.
##
## SCOPE
##   The group computations use abstract normal forms epsilon^i h^j g^k,
##   with i modulo 4 and j,k integral.  In particular, no h^8=1, g^8=1,
##   or other special finite-quotient relation is available to this script.
##   All computations are exact integer computations in the free Laurent
##   exponent lattice on normalized symbols Phi(a,b,c).  No external GAP
##   package, floating-point arithmetic, or sampled cocycle is used.
##
## REPRODUCTION (from this directory)
##   gap --bare -A -r -q --nointeract gamma4_delta_transport_certificate.g
#############################################################################

Read("gamma4_delta_transport_certificate_data.g");
SizeScreen([1000, 1000]);

AddExponent := function(v, key, exponent)
    local i;
    if exponent = 0 then
        return;
    fi;
    for i in [1 .. Length(v)] do
        if v[i][1] = key then
            v[i][2] := v[i][2] + exponent;
            if v[i][2] = 0 then
                Remove(v, i);
            fi;
            return;
        fi;
    od;
    Add(v, [key, exponent]);
end;

CanonicalSparse := function(v)
    local ans, term;
    ans := [];
    for term in v do
        AddExponent(ans, term[1], term[2]);
    od;
    Sort(ans, function(left, right)
        return left[1] < right[1];
    end);
    return ans;
end;

AddSparse := function(left, right)
    local ans, term;
    ans := List(left, ShallowCopy);
    for term in right do
        AddExponent(ans, term[1], term[2]);
    od;
    return CanonicalSparse(ans);
end;

ScaleSparse := function(v, scalar)
    return CanonicalSparse(List(v, term -> [term[1], scalar * term[2]]));
end;

IdentityElt := [0, 0, 0];
Epsilon := [1, 0, 0];
H := [0, 1, 0];
G := [0, 0, 1];

## Multiplication in the abstract Gamma4 normal form epsilon^i h^j g^k.
## Only the epsilon exponent is reduced; the h- and g-exponents remain in Z.
Multiply := function(a, b)
    local i, j, k, ii, jj, kk, first;
    i := a[1]; j := a[2]; k := a[3];
    ii := b[1]; jj := b[2]; kk := b[3];
    if (k mod 2) = 0 then
        first := i + ii;
    else
        first := i - ii - jj;
    fi;
    return [first mod 4, j + jj, k + kk];
end;

InverseElt := function(a)
    local i, j, k;
    i := a[1]; j := a[2]; k := a[3];
    if (k mod 2) = 0 then
        return [(-i) mod 4, -j, -k];
    fi;
    return [(i + j) mod 4, -j, -k];
end;

Conjugate := function(a, x)
    return Multiply(Multiply(a, x), InverseElt(a));
end;

PowerElt := function(a, exponent)
    local ans, i;
    if exponent < 0 then
        return PowerElt(InverseElt(a), -exponent);
    fi;
    ans := IdentityElt;
    for i in [1 .. exponent] do
        ans := Multiply(ans, a);
    od;
    return ans;
end;

PhiKey := function(a, b, c)
    return Concatenation("P(", String(a), ",", String(b), ",", String(c), ")");
end;

PhiAtom := function(a, b, c)
    if a = IdentityElt or b = IdentityElt or c = IdentityElt then
        return [];
    fi;
    return [[PhiKey(a, b, c), 1]];
end;

CocycleDefect := function(a, b, c, d)
    local ans;
    ans := [];
    ans := AddSparse(ans, PhiAtom(b, c, d));
    ans := AddSparse(ans, PhiAtom(a, Multiply(b, c), d));
    ans := AddSparse(ans, PhiAtom(a, b, c));
    ans := AddSparse(ans, ScaleSparse(PhiAtom(Multiply(a, b), c, d), -1));
    ans := AddSparse(ans, ScaleSparse(PhiAtom(a, b, Multiply(c, d)), -1));
    return ans;
end;

## Exponent vector of Phi_x(a,b).
PhiX := function(x, a, b)
    local ans;
    ans := [];
    ans := AddSparse(ans, PhiAtom(a, b, x));
    ans := AddSparse(ans, PhiAtom(Conjugate(Multiply(a, b), x), a, b));
    ans := AddSparse(ans, ScaleSparse(PhiAtom(a, Conjugate(b, x), b), -1));
    return ans;
end;

OmegaPair := function(x, r)
    local u, b, square, common, omega01, omega10;
    u := PowerElt(Epsilon, 2);
    b := Multiply(Epsilon, x);
    square := PowerElt(x, 2);
    common := [];
    common := AddSparse(common, PhiX(x, u, u));
    common := AddSparse(common, PhiX(x, u, r));
    common := AddSparse(common, PhiX(x, square, r));
    common := AddSparse(common, ScaleSparse(PhiX(x, x, x), -2));
    common := AddSparse(common, ScaleSparse(PhiX(x, r, u), -1));
    common := AddSparse(common, ScaleSparse(PhiX(x, r, square), -1));

    omega01 := common;
    omega01 := AddSparse(omega01, PhiX(b, Multiply(u, x), x));
    omega01 := AddSparse(omega01, PhiX(x, b, b));
    omega01 := AddSparse(omega01, ScaleSparse(PhiX(x, b, u), -1));
    omega01 := AddSparse(omega01, ScaleSparse(PhiX(b, u, square), -1));

    omega10 := common;
    omega10 := AddSparse(omega10, PhiX(x, Multiply(u, b), b));
    omega10 := AddSparse(omega10, PhiX(b, x, x));
    omega10 := AddSparse(omega10, ScaleSparse(PhiX(b, x, u), -1));
    omega10 := AddSparse(omega10, ScaleSparse(PhiX(x, u, square), -1));
    return [omega01, omega10];
end;

VerifyCertificate := function(entry, target)
    local combination, term, residual, quad, maxEntryH, maxEntryG,
          maxProductH, maxProductG, x, y, product;
    combination := [];
    maxEntryH := 0; maxEntryG := 0;
    maxProductH := 0; maxProductG := 0;
    for term in entry.certificate do
        quad := term[1];
        combination := AddSparse(
            combination,
            ScaleSparse(CocycleDefect(quad[1], quad[2], quad[3], quad[4]), term[2])
        );
        for x in quad do
            maxEntryH := Maximum(maxEntryH, AbsInt(x[2]));
            maxEntryG := Maximum(maxEntryG, AbsInt(x[3]));
        od;
        for x in [1 .. 3] do
            y := quad[x + 1];
            product := Multiply(quad[x], y);
            maxProductH := Maximum(maxProductH, AbsInt(product[2]));
            maxProductG := Maximum(maxProductG, AbsInt(product[3]));
        od;
    od;
    residual := AddSparse(target, ScaleSparse(combination, -1));
    if Length(residual) <> 0 then
        Error("FAIL: nonzero bar-boundary residual for ", entry.name);
    fi;
    Print(entry.name,
          ": certificate terms=", Length(entry.certificate),
          ", coefficient l1=", Sum(entry.certificate, term -> AbsInt(term[2])),
          ", residual atoms=", Length(residual),
          ", max |h,g| in entries=", [maxEntryH, maxEntryG],
          ", max |h,g| in adjacent products=", [maxProductH, maxProductG],
          "\n");
end;

WPair := OmegaPair(G, H);
D := Multiply(H, G);
YPair := OmegaPair(D, H);
Targets := [
    AddSparse(WPair[1], ScaleSparse(WPair[2], -1)),
    AddSparse(YPair[1], ScaleSparse(WPair[1], -1)),
    AddSparse(YPair[2], ScaleSparse(WPair[2], -1))
];

if Length(OmegaCertificates) <> Length(Targets) then
    Error("FAIL: wrong number of embedded certificates");
fi;

VerifyCertificate(OmegaCertificates[1], Targets[1]);
VerifyCertificate(OmegaCertificates[2], Targets[2]);
VerifyCertificate(OmegaCertificates[3], Targets[3]);

Print("Gamma4 Omega-transport bar-boundary certificate\n");
Print("Arithmetic: exact integers in the normalized bar complex\n");
Print("External GAP packages required: none\n");
Print("Group model: epsilon exponent modulo 4; h,g exponents integral\n");
Print("No h- or g-order relation used\n");
Print("PASS: Omega_01(g,h)=Omega_10(g,h)=Omega_01(hg,h)=Omega_10(hg,h)\n");

QUIT_GAP(0);
