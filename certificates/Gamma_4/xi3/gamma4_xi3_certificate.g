#############################################################################
## Gamma4 Xi3 direct-cancellation certificate (core GAP only)
##
## PURPOSE
##   Verify, in the free Laurent group on formal symbols Phi(a,b,c), the exact
##   two-boundary certificate
##
##       Xi3 = D(h,g,ug,g) / D(epsilon*g,h,ug,g) = 1,
##
##   where
##
##       D(a,b,c,d) = Phi(b,c,d) Phi(a,bc,d) Phi(a,b,c)
##                    ----------------------------------- .
##                    Phi(ab,c,d) Phi(a,b,cd)
##
##   For a normalized 3-cocycle, every D(a,b,c,d) is 1.
##
## SCOPE
##   This file checks only the final, corrected direct Xi3 cancellation after
##   the projective-character reduction
##
##       -sigma(ug)/sigma(ut) = Phi(g,ug,g)
##
##   (which uses sigma(g)=-1).  It neither derives the Y2 line-stability
##   equations nor verifies the full Gamma4 classification theorem.  No GAP
##   package, floating-point arithmetic, finite quotient, or sampled cocycle
##   is used: all computations are exact integer computations in a formal
##   sparse Laurent exponent lattice.
##
## GROUP RELATIONS USED BY THIS CERTIFICATE
##   g(ug)=(ug)g=ut, h(ug)=uhg, and (epsilon*g)h=hg=h g.
##
## REPRODUCTION
##   gap -q gamma4_xi3_certificate.g
#############################################################################

## A sparse Laurent monomial is a sorted list [ [key, exponent], ... ].
## Zero exponents are omitted.  Keys are strings such as Phi(h,g,ug).

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

PhiAtom := function(a, b, c)
    return [[Concatenation("Phi(", a, ",", b, ",", c, ")"), 1]];
end;

## Only products needed by the two displayed cocycle defects are accepted.
## An unsupported product is an error rather than an implicit extra relation.
MultiplyWord := function(a, b)
    if a = "g" and b = "ug" then
        return "ut";
    elif a = "ug" and b = "g" then
        return "ut";
    elif a = "h" and b = "ug" then
        return "uhg";
    elif a = "h" and b = "g" then
        return "hg";
    elif a = "epsilon*g" and b = "h" then
        return "hg";
    fi;
    Error("unsupported formal product: ", a, " * ", b);
end;

CocycleDefect := function(a, b, c, d)
    local ans;
    ans := [];
    ans := AddSparse(ans, PhiAtom(b, c, d));
    ans := AddSparse(ans, PhiAtom(a, MultiplyWord(b, c), d));
    ans := AddSparse(ans, PhiAtom(a, b, c));
    ans := AddSparse(ans, ScaleSparse(PhiAtom(MultiplyWord(a, b), c, d), -1));
    ans := AddSparse(ans, ScaleSparse(PhiAtom(a, b, MultiplyWord(c, d)), -1));
    return ans;
end;

## Laurent exponent vector obtained directly from the corrected Xi3 formula,
## after applying the projective-character relation quoted in SCOPE.
XiDirect := [];
XiDirect := AddSparse(XiDirect, PhiAtom("g",         "ug",  "g"));
XiDirect := AddSparse(XiDirect, PhiAtom("h",         "g",   "ug"));
XiDirect := AddSparse(XiDirect, PhiAtom("h",         "ut",  "g"));
XiDirect := AddSparse(XiDirect, PhiAtom("epsilon*g", "h",   "ut"));
XiDirect := AddSparse(XiDirect, ScaleSparse(PhiAtom("h",         "g",   "ut"), -1));
XiDirect := AddSparse(XiDirect, ScaleSparse(PhiAtom("h",         "ug",  "g"),  -1));
XiDirect := AddSparse(XiDirect, ScaleSparse(PhiAtom("epsilon*g", "h",   "ug"), -1));
XiDirect := AddSparse(XiDirect, ScaleSparse(PhiAtom("epsilon*g", "uhg", "g"),  -1));

D1 := CocycleDefect("h",         "g", "ug", "g");
D2 := CocycleDefect("epsilon*g", "h", "ug", "g");
CertificateVector := AddSparse(D1, ScaleSparse(D2, -1));
Residual := AddSparse(XiDirect, ScaleSparse(CertificateVector, -1));

if XiDirect <> CertificateVector then
    Error("FAIL: the direct Xi3 exponent vector is not D1-D2");
fi;
if Length(Residual) <> 0 then
    Error("FAIL: nonzero Laurent residual in the bar-boundary certificate");
fi;

Print("Gamma4 Xi3 direct-cancellation certificate\n");
Print("Arithmetic: exact integers in a formal sparse Laurent lattice\n");
Print("External GAP packages loaded: none\n");
Print("Relation check: (epsilon*g)h = hg = h*g\n");
Print("Direct Xi3 Laurent atoms: ", Length(XiDirect), "\n");
Print("Certificate: +d[h|g|ug|g] - d[epsilon*g|h|ug|g]\n");
Print("Residual Laurent atoms: ", Length(Residual), "\n");
Print("PASS: Xi3 = D(h,g,ug,g)/D(epsilon*g,h,ug,g) = 1\n");

QUIT_GAP(0);
