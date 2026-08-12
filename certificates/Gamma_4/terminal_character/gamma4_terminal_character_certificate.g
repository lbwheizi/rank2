#############################################################################
## Gamma4 terminal-root self-character certificate (core GAP only)
##
## PURPOSE
##   Verify exactly the Laurent identity arising in the terminal-root
##   self-character calculation:
##
##       Q = product_{(a,b,c,d),n in CERTIFICATE}
##               D(a,b,c,d)^n,
##
##   with 146 normalized 3-cocycle defects and coefficient l1 norm 245.
##   Here
##
##       D(a,b,c,d) = Phi(b,c,d) Phi(a,bc,d) Phi(a,b,c)
##                    / (Phi(ab,c,d) Phi(a,b,cd)).
##
##   Thus Q=1 for every normalized 3-cocycle.  The vector Q is the exact
##   target remaining in the terminal-character calculation after the
##   projective-character substitutions rho(h)=sigma(g)=-1, Delta4=1,
##   and epsilon^4=1.  This version does not take that target on trust:
##   it first reconstructs K_q, lambda_1 and lambda_2 from the six displayed
##   action coefficients L_qa,L_qp,L_ha,L_hy,L_hr,L_hY, performs the stated
##   projective-character eliminations, and checks that the result is the
##   recorded 76-atom Laurent target before applying the 146 defects.
##
## SCOPE
##   No external GAP packages, floating-point arithmetic, finite quotient,
##   or sampled cocycle is used.  Group elements are abstract normal forms
##   epsilon^i h^j g^k with i modulo 4 and j,k integral.  The computation is
##   exact in the free Laurent exponent lattice on formal Phi(a,b,c).
##
## REPRODUCTION
##   gap -q gamma4_terminal_character_certificate.g
#############################################################################

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

L1Sparse := function(v)
    return Sum(v, term -> AbsInt(term[2]));
end;

IdentityElt := [ 0, 0, 0 ];

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

Conjugate := function(a, b)
    return Multiply(Multiply(a, b), InverseElt(a));
end;

PowerElt := function(a, n)
    local ans, i;
    if n < 0 then
        return PowerElt(InverseElt(a), -n);
    fi;
    ans := IdentityElt;
    for i in [1 .. n] do
        ans := Multiply(ans, a);
    od;
    return ans;
end;

PhiAtom := function(a, b, c)
    local key;
    if a = IdentityElt or b = IdentityElt or c = IdentityElt then
        return [];
    fi;
    key := Concatenation("P(", String(a), ",", String(b), ",", String(c), ")");
    return [[key, 1]];
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

MultiplyMonomial := function(left, right)
    return AddSparse(left, right);
end;

DivideMonomial := function(left, right)
    return AddSparse(left, ScaleSparse(right, -1));
end;

PowerMonomial := function(value, exponent)
    return ScaleSparse(value, exponent);
end;

CharacterAtom := function(name, a)
    if a = IdentityElt then
        return [];
    fi;
    return [[Concatenation(name, String(a)), 1]];
end;

PhiLocal := function(y, a, b)
    local ans;
    ans := PhiAtom(a, b, y);
    ans := AddSparse(ans, PhiAtom(Conjugate(Multiply(a, b), y), a, b));
    ans := AddSparse(ans, ScaleSparse(PhiAtom(a, Conjugate(b, y), b), -1));
    return ans;
end;

PhiUpper := function(a, x, y)
    local ans;
    ans := PhiAtom(a, x, y);
    ans := AddSparse(ans, PhiAtom(Conjugate(a, x), Conjugate(a, y), a));
    ans := AddSparse(ans, ScaleSparse(PhiAtom(Conjugate(a, x), a, y), -1));
    return ans;
end;

InverseCharacterValue := function(base, generator, value)
    return DivideMonomial(PhiLocal(base, generator, InverseElt(generator)), value);
end;

LettersPower := function(base, generator, value, exponent)
    local letters, i, inverseGenerator, inverseValue;
    letters := [];
    if exponent >= 0 then
        for i in [1 .. exponent] do
            Add(letters, [generator, value]);
        od;
        return letters;
    fi;
    inverseGenerator := InverseElt(generator);
    inverseValue := InverseCharacterValue(base, generator, value);
    for i in [1 .. -exponent] do
        Add(letters, [inverseGenerator, inverseValue]);
    od;
    return letters;
end;

CharacterWord := function(base, letters)
    local current, value, letter;
    current := IdentityElt;
    value := [];
    for letter in letters do
        value := DivideMonomial(
            MultiplyMonomial(value, letter[2]),
            PhiLocal(base, current, letter[1])
        );
        current := Multiply(current, letter[1]);
    od;
    return [value, current];
end;

ExponentOf := function(value, key)
    local term;
    for term in value do
        if term[1] = key then
            return term[2];
        fi;
    od;
    return 0;
end;

RemoveKey := function(value, key)
    return Filtered(value, term -> term[1] <> key);
end;

SubstituteMonomial := function(value, key, replacement)
    local exponent, rest;
    exponent := ExponentOf(value, key);
    if exponent = 0 then
        return value;
    fi;
    rest := RemoveKey(value, key);
    return MultiplyMonomial(rest, PowerMonomial(replacement, exponent));
end;

## Reconstruct the source expression in Lemma C.9.
Epsilon := [1, 0, 0];
H := [0, 1, 0];
G := [0, 0, 1];
R := InverseElt(Epsilon);
U := PowerElt(Epsilon, 2);
T := PowerElt(G, 2);
ZElt := Multiply(R, PowerElt(H, 2));
XElt := Multiply(U, G);
DElt := Multiply(H, G);
QDegree := Multiply(XElt, DElt);
EGElt := Multiply(Epsilon, G);
RGElt := Multiply(R, G);
EDElt := Multiply(Epsilon, DElt);
RDElt := Multiply(R, DElt);
UTElt := Multiply(U, T);
CMinus := Multiply(Multiply(InverseElt(Epsilon), Multiply(QDegree, XElt)), H);
CPlus := Multiply(Multiply(InverseElt(Epsilon), Multiply(H, Multiply(G, H))), Epsilon);

if CMinus <> Multiply(Multiply(U, ZElt), PowerElt(G, 3)) then
    Error("FAIL: c_- normal form is not uzg^3");
fi;
if CPlus <> Multiply(Multiply(U, ZElt), G) then
    Error("FAIL: c_+ normal form is not uzg");
fi;

RhoEpsilon := CharacterAtom("R", Epsilon);
RhoH := CharacterAtom("R", H);
RhoT := CharacterAtom("R", T);
SigmaU := CharacterAtom("S", U);
SigmaZ := CharacterAtom("S", ZElt);
SigmaG := CharacterAtom("S", G);
SignAtom := [["SIGN", 1]];

Rho := function(a)
    local i, j, k, letters, result;
    i := a[1]; j := a[2]; k := a[3];
    if (k mod 2) <> 0 then
        Error("rho requested outside C_G(h): ", a);
    fi;
    letters := [];
    Append(letters, LettersPower(H, Epsilon, RhoEpsilon, i));
    Append(letters, LettersPower(H, H, RhoH, j));
    Append(letters, LettersPower(H, T, RhoT, QuoInt(k, 2)));
    result := CharacterWord(H, letters);
    if result[2] <> a then
        Error("rho word has wrong endpoint: ", a, " versus ", result[2]);
    fi;
    return result[1];
end;

SigmaValue := function(a)
    local i, j, k, m, rhs, aa, letters, result;
    i := a[1]; j := a[2]; k := a[3];
    if (j mod 2) <> 0 then
        Error("sigma requested outside C_G(g): ", a);
    fi;
    m := QuoInt(j, 2);
    rhs := (i + m) mod 4;
    if not rhs in [0, 2] then
        Error("sigma normal form has invalid u exponent: ", a);
    fi;
    aa := QuoInt(rhs, 2);
    letters := [];
    Append(letters, LettersPower(G, U, SigmaU, aa));
    Append(letters, LettersPower(G, ZElt, SigmaZ, m));
    Append(letters, LettersPower(G, G, SigmaG, k));
    result := CharacterWord(G, letters);
    if result[2] <> a then
        Error("sigma word has wrong endpoint: ", a, " versus ", result[2]);
    fi;
    return result[1];
end;

## A_1 and the global first-adjoint scalar Delta_4.
A1 := DivideMonomial(PhiAtom(XElt, EGElt, H), PhiAtom(RGElt, XElt, H));
AlphaCoeff := DivideMonomial(Rho(R), PhiLocal(H, G, R));
BetaCoeff := MultiplyMonomial(
    MultiplyMonomial(AlphaCoeff, PhiLocal(G, Multiply(R, H), H)),
    SigmaValue(ZElt)
);
Delta4 := MultiplyMonomial(
    MultiplyMonomial(BetaCoeff, PhiLocal(H, G, G)),
    MultiplyMonomial(Rho(T), [])
);

## The six displayed action coefficients.
Lqa := DivideMonomial(
    MultiplyMonomial(
        MultiplyMonomial(PhiLocal(EGElt, QDegree, XElt),
                         PhiLocal(G, Multiply(QDegree, XElt), H)),
        SigmaValue(CMinus)
    ),
    PhiLocal(G, Epsilon, CMinus)
);
Lqp := DivideMonomial(
    MultiplyMonomial(
        MultiplyMonomial(PhiUpper(QDegree, XElt, H),
                         PhiLocal(G, QDegree, Epsilon)),
        MultiplyMonomial(SigmaValue(UTElt), Rho(QDegree))
    ),
    PhiLocal(G, H, UTElt)
);
Lha := DivideMonomial(
    MultiplyMonomial(PhiLocal(G, H, Epsilon), SigmaValue(InverseElt(G))),
    MultiplyMonomial(PhiLocal(G, Multiply(XElt, H), InverseElt(G)),
                     PhiLocal(G, XElt, H))
);
Lhy := DivideMonomial(
    MultiplyMonomial(
        MultiplyMonomial(PhiUpper(H, EGElt, H), PhiLocal(G, H, H)),
        MultiplyMonomial(SigmaValue(ZElt), Rho(H))
    ),
    PhiLocal(G, Epsilon, ZElt)
);
Lhr := DivideMonomial(
    MultiplyMonomial(
        MultiplyMonomial(PhiLocal(XElt, H, Multiply(G, H)),
                         PhiLocal(G, Multiply(H, Multiply(G, H)), Epsilon)),
        SigmaValue(CPlus)
    ),
    PhiLocal(G, Epsilon, CPlus)
);
LhY := MultiplyMonomial(
    MultiplyMonomial(
        MultiplyMonomial(BetaCoeff, PhiUpper(XElt, G, Multiply(R, H))),
        PhiLocal(H, XElt, G)
    ),
    MultiplyMonomial(
        MultiplyMonomial(Rho(UTElt), PhiUpper(H, G, H)),
        MultiplyMonomial(SigmaValue(XElt), Rho(H))
    )
);

Kq := MultiplyMonomial(
    MultiplyMonomial(A1, PhiUpper(QDegree, RGElt, EDElt)),
    MultiplyMonomial(Lqa, Lqp)
);
Lambda1 := MultiplyMonomial(PhiUpper(H, EGElt, RDElt),
                            MultiplyMonomial(Lhr, LhY));
Lambda2 := DivideMonomial(
    MultiplyMonomial(PhiUpper(H, XElt, DElt), MultiplyMonomial(Lha, Lhy)),
    A1
);

SourceBeforeElimination := DivideMonomial(
    MultiplyMonomial(Kq, Lambda2),
    MultiplyMonomial(Rho(H), MultiplyMonomial(PowerMonomial(SigmaValue(G), 2),
                                              Lambda1))
);

## rho(h)=sigma(g)=-1.
SourceAfterCharacters := SubstituteMonomial(SourceBeforeElimination,
                                             RhoH[1][1], SignAtom);
SourceAfterCharacters := SubstituteMonomial(SourceAfterCharacters,
                                             SigmaG[1][1], SignAtom);

## Delta_4=1, solved for sigma(z).
DeltaReduced := SubstituteMonomial(Delta4, RhoH[1][1], SignAtom);
DeltaReduced := SubstituteMonomial(DeltaReduced, SigmaG[1][1], SignAtom);
SigmaZExponent := ExponentOf(DeltaReduced, SigmaZ[1][1]);
if not SigmaZExponent in [-1, 1] then
    Error("FAIL: Delta_4 is not linear in sigma(z)");
fi;
DeltaRest := RemoveKey(DeltaReduced, SigmaZ[1][1]);
SigmaZValue := PowerMonomial(DeltaRest, -SigmaZExponent);
SourceAfterCharacters := SubstituteMonomial(SourceAfterCharacters,
                                             SigmaZ[1][1], SigmaZValue);

## epsilon^4=1: rho(epsilon)^4 times the projective correction is one.
EpsilonRelationData := CharacterWord(
    H, LettersPower(H, Epsilon, RhoEpsilon, 4)
);
if EpsilonRelationData[2] <> IdentityElt then
    Error("FAIL: epsilon^4 word does not return to the identity");
fi;
if ExponentOf(EpsilonRelationData[1], RhoEpsilon[1][1]) <> 4 then
    Error("FAIL: epsilon^4 relation has wrong rho(epsilon) exponent");
fi;
EpsilonRest := RemoveKey(EpsilonRelationData[1], RhoEpsilon[1][1]);
if ExponentOf(SourceAfterCharacters, RhoEpsilon[1][1]) <> -4 then
    Error("FAIL: source has wrong residual rho(epsilon) exponent");
fi;
SourceAfterCharacters := RemoveKey(SourceAfterCharacters, RhoEpsilon[1][1]);
SourceTarget := MultiplyMonomial(SourceAfterCharacters, EpsilonRest);

if ForAny(SourceTarget, term -> not StartsWith(term[1], "P(")) then
    Error("FAIL: source-to-target elimination left a character atom");
fi;

## Exact 76-atom Laurent target from the terminal-character reduction.
Target := [];
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 0, 0, 1 ], [ 0, 0, -1 ], [ 0, 0, 1 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 0, 0, 1 ], [ 0, 1, 0 ], [ 2, 0, 1 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 0, 0, 1 ], [ 0, 1, 0 ], [ 3, 0, 0 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 0, 0, 1 ], [ 1, 2, 1 ], [ 0, 0, 1 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 0, 0, 1 ], [ 1, 2, 2 ], [ 0, 0, 1 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 0, 0, 1 ], [ 2, 0, 1 ], [ 0, 0, 1 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 0, 0, 1 ], [ 2, 0, 1 ], [ 3, 1, 0 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 0, 0, 1 ], [ 3, 0, 0 ], [ 0, 1, 0 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 0, 0, 1 ], [ 3, 1, 0 ], [ 0, 1, 0 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 0, 1, 0 ], [ 0, 0, 1 ], [ 0, 1, 0 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 0, 1, 0 ], [ 0, 0, 1 ], [ 2, 0, 2 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 0, 1, 0 ], [ 0, 1, 0 ], [ 0, 0, 1 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 0, 1, 0 ], [ 1, 0, 0 ], [ 0, 0, 1 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 0, 1, 0 ], [ 1, 0, 0 ], [ 0, 1, 0 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 0, 1, 0 ], [ 1, 0, 0 ], [ 1, 0, 0 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 0, 1, 0 ], [ 1, 1, 0 ], [ 0, 0, 2 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 0, 1, 0 ], [ 2, 0, 0 ], [ 0, 0, 2 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 0, 1, 0 ], [ 2, 0, 1 ], [ 0, 0, 1 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 0, 1, 0 ], [ 2, 0, 1 ], [ 0, 1, 1 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 0, 1, 0 ], [ 2, 0, 1 ], [ 1, 0, 0 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 0, 1, 0 ], [ 2, 0, 2 ], [ 0, 0, 1 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 0, 1, 0 ], [ 3, 0, 0 ], [ 1, 0, 0 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 0, 1, 0 ], [ 3, 1, 1 ], [ 2, 0, 1 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 1, 0, 0 ], [ 0, 0, 1 ], [ 1, 2, 1 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 1, 0, 0 ], [ 0, 0, 1 ], [ 1, 2, 3 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 1, 0, 0 ], [ 0, 0, 1 ], [ 3, 2, 0 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 1, 0, 0 ], [ 0, 1, 0 ], [ 1, 0, 0 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 1, 0, 0 ], [ 1, 0, 0 ], [ 0, 1, 0 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 1, 0, 0 ], [ 1, 2, 1 ], [ 0, 0, 1 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 1, 0, 0 ], [ 1, 2, 3 ], [ 0, 0, 1 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 1, 0, 0 ], [ 3, 2, 0 ], [ 0, 0, 1 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 1, 0, 1 ], [ 0, 1, 0 ], [ 1, 1, 2 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 1, 0, 1 ], [ 0, 1, 0 ], [ 2, 0, 2 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 1, 0, 1 ], [ 1, 1, 2 ], [ 0, 1, 0 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 1, 0, 1 ], [ 1, 1, 2 ], [ 1, 0, 0 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 1, 1, 0 ], [ 0, 0, 2 ], [ 0, 1, 0 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 1, 1, 0 ], [ 0, 1, 0 ], [ 0, 0, 2 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 1, 1, 1 ], [ 0, 0, -1 ], [ 0, 0, 1 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 1, 1, 1 ], [ 0, 0, 1 ], [ 0, 0, -1 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 1, 1, 2 ], [ 1, 0, 0 ], [ 0, 0, 1 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 1, 1, 2 ], [ 2, 0, 1 ], [ 0, 1, 0 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 1, 1, 2 ], [ 2, 0, 1 ], [ 1, 0, 0 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 1, 1, 2 ], [ 2, 0, 1 ], [ 1, 0, 1 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 1, 1, 2 ], [ 3, 0, 1 ], [ 1, 1, 1 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 1, 1, 2 ], [ 3, 0, 1 ], [ 2, 0, 1 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 2, 0, 0 ], [ 0, 0, 2 ], [ 0, 1, 0 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 2, 0, 0 ], [ 0, 1, 0 ], [ 0, 0, 2 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 2, 0, 1 ], [ 0, 0, 1 ], [ 0, 1, 0 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 2, 0, 1 ], [ 0, 0, 1 ], [ 3, 1, 0 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 2, 0, 1 ], [ 0, 1, 0 ], [ 0, 0, 1 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 2, 0, 1 ], [ 0, 1, 0 ], [ 0, 1, 0 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 2, 0, 1 ], [ 0, 1, 1 ], [ 0, 1, 0 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 2, 0, 1 ], [ 0, 1, 1 ], [ 1, 1, 2 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 2, 0, 1 ], [ 1, 0, 0 ], [ 1, 2, 1 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 2, 0, 1 ], [ 1, 0, 0 ], [ 1, 2, 3 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 2, 0, 1 ], [ 1, 0, 0 ], [ 3, 2, 0 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 2, 0, 1 ], [ 1, 0, 1 ], [ 0, 1, 0 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 2, 0, 1 ], [ 1, 1, 2 ], [ 1, 1, 1 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 2, 0, 1 ], [ 1, 1, 2 ], [ 2, 0, 1 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 2, 0, 1 ], [ 3, 1, 0 ], [ 0, 0, 1 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 2, 0, 1 ], [ 3, 1, 3 ], [ 0, 1, 0 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 2, 0, 1 ], [ 3, 2, 1 ], [ 1, 0, 0 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 3, 0, 0 ], [ 0, 1, 0 ], [ 1, 0, 0 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 3, 0, 0 ], [ 1, 0, 0 ], [ 0, 1, 0 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 3, 0, 1 ], [ 0, 1, 0 ], [ 0, 1, 1 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 3, 0, 1 ], [ 0, 1, 0 ], [ 1, 0, 0 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 3, 0, 1 ], [ 1, 1, 1 ], [ 0, 0, -1 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 3, 0, 1 ], [ 1, 1, 1 ], [ 0, 1, 0 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 3, 0, 1 ], [ 2, 0, 1 ], [ 0, 1, 0 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 3, 1, 0 ], [ 0, 0, 1 ], [ 3, 0, 0 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 3, 1, 0 ], [ 0, 1, 0 ], [ 0, 0, 1 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 3, 1, 0 ], [ 1, 0, 1 ], [ 0, 1, 0 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 3, 1, 3 ], [ 0, 1, 0 ], [ 0, 0, 1 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 3, 1, 3 ], [ 1, 0, 1 ], [ 0, 1, 0 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 3, 2, 1 ], [ 1, 0, 0 ], [ 0, 0, 1 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 3, 2, 1 ], [ 2, 0, 1 ], [ 1, 0, 0 ]), 1));

## Exact 146-term cocycle-defect certificate.
Certificate := [
    [ [ [ 0, 0, 1 ], [ 0, 1, 0 ], [ 3, 0, 0 ], [ 0, 1, 0 ] ], -2 ],
    [ [ [ 0, 0, 1 ], [ 0, 1, 0 ], [ 3, 0, 0 ], [ 1, 0, 0 ] ], 2 ],
    [ [ [ 0, 0, 1 ], [ 3, 0, 1 ], [ 0, 1, 0 ], [ 1, 0, 1 ] ], 4 ],
    [ [ [ 0, 0, 1 ], [ 3, 0, 1 ], [ 1, 0, 0 ], [ 0, 0, 1 ] ], -4 ],
    [ [ [ 0, 0, 1 ], [ 3, 0, 1 ], [ 1, 0, 0 ], [ 3, 1, 0 ] ], 4 ],
    [ [ [ 0, 0, 1 ], [ 3, 0, 1 ], [ 1, 0, 1 ], [ 3, 1, 0 ] ], -4 ],
    [ [ [ 0, 0, 1 ], [ 3, 1, 0 ], [ 1, 0, 0 ], [ 1, 0, 1 ] ], -1 ],
    [ [ [ 0, 0, 1 ], [ 3, 1, 0 ], [ 1, 0, 0 ], [ 1, 1, 2 ] ], -1 ],
    [ [ [ 0, 0, 1 ], [ 3, 1, 0 ], [ 2, 0, 2 ], [ 0, 1, 0 ] ], 1 ],
    [ [ [ 0, 0, 2 ], [ 0, 1, 0 ], [ 1, 0, 0 ], [ 1, 0, 0 ] ], 1 ],
    [ [ [ 0, 0, 2 ], [ 1, 0, 0 ], [ 1, 0, 0 ], [ 3, 1, 0 ] ], 1 ],
    [ [ [ 0, 0, 2 ], [ 1, 0, 0 ], [ 3, 1, 0 ], [ 1, 0, 0 ] ], -2 ],
    [ [ [ 0, 0, 2 ], [ 2, 0, 0 ], [ 3, 1, 0 ], [ 1, 0, 0 ] ], -1 ],
    [ [ [ 0, 1, 0 ], [ 0, 0, 1 ], [ 0, 1, 0 ], [ 1, 0, 1 ] ], 1 ],
    [ [ [ 0, 1, 0 ], [ 0, 0, 1 ], [ 0, 1, 0 ], [ 2, 0, 1 ] ], -1 ],
    [ [ [ 0, 1, 0 ], [ 0, 0, 1 ], [ 2, 0, 1 ], [ 0, 0, 1 ] ], 3 ],
    [ [ [ 0, 1, 0 ], [ 0, 0, 1 ], [ 2, 0, 1 ], [ 3, 1, 0 ] ], -3 ],
    [ [ [ 0, 1, 0 ], [ 0, 0, 1 ], [ 2, 0, 2 ], [ 3, 1, 0 ] ], 4 ],
    [ [ [ 0, 1, 0 ], [ 0, 0, 1 ], [ 2, 1, 1 ], [ 1, 0, 1 ] ], -4 ],
    [ [ [ 0, 1, 0 ], [ 0, 0, 2 ], [ 1, 0, 0 ], [ 1, 0, 0 ] ], -1 ],
    [ [ [ 0, 1, 0 ], [ 1, 0, 0 ], [ 0, 0, 2 ], [ 0, 1, 0 ] ], 2 ],
    [ [ [ 0, 1, 0 ], [ 1, 0, 0 ], [ 0, 1, 0 ], [ 0, 0, 2 ] ], -4 ],
    [ [ [ 0, 1, 0 ], [ 1, 0, 0 ], [ 1, 0, 0 ], [ 0, 0, 2 ] ], 2 ],
    [ [ [ 0, 1, 0 ], [ 1, 0, 0 ], [ 1, 0, 2 ], [ 3, 1, 0 ] ], 2 ],
    [ [ [ 0, 1, 0 ], [ 1, 0, 2 ], [ 0, 1, 0 ], [ 1, 0, 1 ] ], 4 ],
    [ [ [ 0, 1, 0 ], [ 1, 0, 2 ], [ 1, 0, 0 ], [ 0, 0, 1 ] ], -4 ],
    [ [ [ 0, 1, 0 ], [ 1, 0, 2 ], [ 1, 0, 0 ], [ 3, 1, 0 ] ], 5 ],
    [ [ [ 0, 1, 0 ], [ 1, 0, 2 ], [ 1, 0, 1 ], [ 3, 1, 0 ] ], -4 ],
    [ [ [ 0, 1, 0 ], [ 2, 0, 0 ], [ 0, 0, 2 ], [ 2, 0, 0 ] ], 1 ],
    [ [ [ 0, 1, 0 ], [ 2, 0, 0 ], [ 1, 0, 0 ], [ 1, 0, 0 ] ], -1 ],
    [ [ [ 0, 1, 0 ], [ 2, 0, 0 ], [ 2, 0, 0 ], [ 2, 0, 2 ] ], -1 ],
    [ [ [ 0, 1, 0 ], [ 2, 0, 1 ], [ 0, 0, 1 ], [ 3, 1, 0 ] ], -2 ],
    [ [ [ 0, 1, 0 ], [ 2, 0, 1 ], [ 1, 0, 0 ], [ 3, 1, 1 ] ], 1 ],
    [ [ [ 0, 1, 0 ], [ 2, 0, 2 ], [ 2, 0, 0 ], [ 1, 0, 0 ] ], -1 ],
    [ [ [ 0, 1, 0 ], [ 2, 0, 2 ], [ 3, 0, 0 ], [ 0, 1, 0 ] ], -1 ],
    [ [ [ 0, 1, 0 ], [ 2, 0, 2 ], [ 3, 1, 0 ], [ 1, 0, 0 ] ], -1 ],
    [ [ [ 0, 1, 0 ], [ 3, 0, 0 ], [ 1, 0, 0 ], [ 3, 1, 0 ] ], -2 ],
    [ [ [ 0, 1, 0 ], [ 3, 1, 1 ], [ 2, 0, 1 ], [ 1, 0, 0 ] ], 1 ],
    [ [ [ 0, 1, 1 ], [ 0, 1, 0 ], [ 1, 0, 0 ], [ 0, 0, 1 ] ], -1 ],
    [ [ [ 0, 1, 1 ], [ 0, 1, 0 ], [ 2, 0, 1 ], [ 1, 0, 0 ] ], -2 ],
    [ [ [ 0, 1, 1 ], [ 1, 0, 0 ], [ 0, 1, 0 ], [ 0, 0, 1 ] ], 1 ],
    [ [ [ 0, 1, 1 ], [ 1, 0, 0 ], [ 1, 0, 1 ], [ 0, 1, 0 ] ], -1 ],
    [ [ [ 0, 1, 1 ], [ 1, 0, 0 ], [ 1, 1, 1 ], [ 0, 0, 1 ] ], 1 ],
    [ [ [ 0, 1, 1 ], [ 1, 0, 0 ], [ 3, 1, 0 ], [ 2, 0, 1 ] ], 1 ],
    [ [ [ 0, 1, 1 ], [ 2, 0, 1 ], [ 0, 0, 1 ], [ 0, 1, 0 ] ], 1 ],
    [ [ [ 0, 1, 1 ], [ 2, 0, 1 ], [ 3, 1, 0 ], [ 0, 0, 1 ] ], -1 ],
    [ [ [ 0, 1, 1 ], [ 2, 1, 1 ], [ 1, 0, 0 ], [ 0, 0, 1 ] ], 4 ],
    [ [ [ 0, 1, 1 ], [ 3, 0, 1 ], [ 0, 1, 0 ], [ 1, 0, 0 ] ], -2 ],
    [ [ [ 0, 1, 1 ], [ 3, 0, 1 ], [ 1, 0, 0 ], [ 0, 0, 1 ] ], 2 ],
    [ [ [ 0, 1, 1 ], [ 3, 0, 1 ], [ 1, 0, 0 ], [ 3, 1, 0 ] ], -2 ],
    [ [ [ 0, 1, 1 ], [ 3, 0, 1 ], [ 1, 0, 1 ], [ 3, 1, 0 ] ], 2 ],
    [ [ [ 0, 1, 1 ], [ 3, 0, 1 ], [ 1, 1, 0 ], [ 0, 0, 1 ] ], -2 ],
    [ [ [ 1, 0, 0 ], [ 0, 0, 1 ], [ 0, 1, 0 ], [ 1, 1, 2 ] ], 1 ],
    [ [ [ 1, 0, 0 ], [ 0, 0, 1 ], [ 1, 2, 1 ], [ 0, 0, 1 ] ], -2 ],
    [ [ [ 1, 0, 0 ], [ 0, 0, 1 ], [ 1, 2, 2 ], [ 0, 0, 1 ] ], -1 ],
    [ [ [ 1, 0, 0 ], [ 0, 1, 0 ], [ 0, 0, 1 ], [ 0, 1, 0 ] ], -1 ],
    [ [ [ 1, 0, 0 ], [ 0, 1, 0 ], [ 0, 0, 2 ], [ 1, 0, 0 ] ], 2 ],
    [ [ [ 1, 0, 0 ], [ 0, 1, 0 ], [ 1, 0, 0 ], [ 0, 0, 2 ] ], -2 ],
    [ [ [ 1, 0, 0 ], [ 1, 0, 0 ], [ 0, 0, 2 ], [ 0, 1, 0 ] ], -2 ],
    [ [ [ 1, 0, 0 ], [ 1, 0, 0 ], [ 0, 1, 0 ], [ 0, 0, 2 ] ], 2 ],
    [ [ [ 1, 0, 0 ], [ 1, 0, 0 ], [ 3, 1, 0 ], [ 1, 0, 0 ] ], -1 ],
    [ [ [ 1, 0, 0 ], [ 1, 0, 1 ], [ 1, 0, 0 ], [ 1, 2, 1 ] ], 1 ],
    [ [ [ 1, 0, 0 ], [ 1, 0, 1 ], [ 1, 0, 0 ], [ 3, 1, 0 ] ], -1 ],
    [ [ [ 1, 0, 0 ], [ 1, 0, 1 ], [ 3, 1, 1 ], [ 0, 1, 0 ] ], -1 ],
    [ [ [ 1, 0, 0 ], [ 1, 0, 2 ], [ 3, 1, 0 ], [ 1, 0, 0 ] ], -2 ],
    [ [ [ 1, 0, 0 ], [ 3, 1, 0 ], [ 0, 0, 1 ], [ 1, 1, 1 ] ], -1 ],
    [ [ [ 1, 0, 0 ], [ 3, 1, 0 ], [ 0, 1, 0 ], [ 0, 0, 1 ] ], -1 ],
    [ [ [ 1, 0, 0 ], [ 3, 1, 0 ], [ 1, 0, 0 ], [ 1, 0, 0 ] ], 3 ],
    [ [ [ 1, 0, 0 ], [ 3, 1, 0 ], [ 1, 0, 0 ], [ 3, 1, 1 ] ], -1 ],
    [ [ [ 1, 0, 0 ], [ 3, 1, 0 ], [ 1, 0, 1 ], [ 0, 1, 0 ] ], 1 ],
    [ [ [ 1, 0, 0 ], [ 3, 1, 0 ], [ 2, 0, 0 ], [ 0, 0, 2 ] ], 3 ],
    [ [ [ 1, 0, 0 ], [ 3, 1, 0 ], [ 2, 0, 1 ], [ 0, 0, 1 ] ], -1 ],
    [ [ [ 1, 0, 0 ], [ 3, 1, 0 ], [ 2, 0, 2 ], [ 0, 1, 0 ] ], 1 ],
    [ [ [ 1, 0, 0 ], [ 3, 1, 1 ], [ 2, 0, 1 ], [ 0, 1, 0 ] ], -1 ],
    [ [ [ 1, 0, 1 ], [ 0, 1, 0 ], [ 2, 0, 2 ], [ 3, 1, 0 ] ], -2 ],
    [ [ [ 1, 0, 1 ], [ 1, 0, 0 ], [ 1, 1, 2 ], [ 0, 1, 0 ] ], 1 ],
    [ [ [ 1, 0, 1 ], [ 1, 0, 0 ], [ 1, 2, 1 ], [ 0, 0, 1 ] ], -1 ],
    [ [ [ 1, 0, 1 ], [ 1, 0, 0 ], [ 3, 1, 0 ], [ 2, 0, 2 ] ], 1 ],
    [ [ [ 1, 0, 1 ], [ 1, 1, 2 ], [ 0, 1, 0 ], [ 0, 0, 1 ] ], 1 ],
    [ [ [ 1, 0, 1 ], [ 1, 1, 2 ], [ 1, 0, 0 ], [ 3, 1, 0 ] ], 2 ],
    [ [ [ 1, 0, 1 ], [ 3, 1, 1 ], [ 0, 1, 0 ], [ 0, 0, 1 ] ], 1 ],
    [ [ [ 1, 0, 1 ], [ 3, 1, 1 ], [ 1, 0, 1 ], [ 0, 1, 0 ] ], -1 ],
    [ [ [ 1, 0, 1 ], [ 3, 1, 1 ], [ 2, 0, 1 ], [ 0, 1, 1 ] ], 1 ],
    [ [ [ 1, 0, 1 ], [ 3, 1, 1 ], [ 2, 0, 1 ], [ 1, 0, 0 ] ], -1 ],
    [ [ [ 1, 1, 0 ], [ 0, 0, 1 ], [ 0, 0, -1 ], [ 0, 0, 1 ] ], -3 ],
    [ [ [ 1, 1, 0 ], [ 0, 0, 2 ], [ 1, 0, 0 ], [ 3, 1, 0 ] ], 2 ],
    [ [ [ 1, 1, 1 ], [ 0, 0, 1 ], [ 0, 0, -1 ], [ 0, 0, 1 ] ], 4 ],
    [ [ [ 1, 1, 2 ], [ 0, 1, 0 ], [ 1, 0, 0 ], [ 0, 0, 1 ] ], -2 ],
    [ [ [ 1, 1, 2 ], [ 1, 0, 0 ], [ 0, 0, 1 ], [ 0, 1, 0 ] ], -1 ],
    [ [ [ 1, 1, 2 ], [ 1, 0, 0 ], [ 3, 1, 0 ], [ 0, 0, 1 ] ], 1 ],
    [ [ [ 1, 1, 2 ], [ 2, 0, 1 ], [ 0, 1, 0 ], [ 0, 0, 1 ] ], -2 ],
    [ [ [ 1, 1, 2 ], [ 2, 0, 1 ], [ 1, 0, 0 ], [ 3, 1, 0 ] ], -2 ],
    [ [ [ 1, 1, 2 ], [ 2, 0, 1 ], [ 1, 0, 1 ], [ 0, 1, 0 ] ], 1 ],
    [ [ [ 1, 1, 2 ], [ 3, 0, 1 ], [ 2, 0, 1 ], [ 0, 1, 0 ] ], -1 ],
    [ [ [ 2, 0, 0 ], [ 0, 0, 2 ], [ 0, 1, 0 ], [ 2, 0, 0 ] ], 1 ],
    [ [ [ 2, 0, 0 ], [ 0, 0, 2 ], [ 2, 0, 0 ], [ 0, 1, 0 ] ], -1 ],
    [ [ [ 2, 0, 0 ], [ 0, 1, 0 ], [ 0, 0, 2 ], [ 2, 0, 0 ] ], -1 ],
    [ [ [ 2, 0, 0 ], [ 0, 1, 0 ], [ 1, 0, 0 ], [ 1, 0, 0 ] ], 1 ],
    [ [ [ 2, 0, 0 ], [ 0, 1, 0 ], [ 2, 0, 0 ], [ 2, 0, 2 ] ], 1 ],
    [ [ [ 2, 0, 0 ], [ 1, 0, 0 ], [ 0, 1, 0 ], [ 1, 0, 0 ] ], -1 ],
    [ [ [ 2, 0, 0 ], [ 1, 0, 0 ], [ 1, 0, 0 ], [ 0, 1, 0 ] ], 1 ],
    [ [ [ 2, 0, 0 ], [ 2, 0, 0 ], [ 0, 1, 0 ], [ 2, 0, 2 ] ], -1 ],
    [ [ [ 2, 0, 0 ], [ 2, 0, 0 ], [ 2, 0, 2 ], [ 0, 1, 0 ] ], 1 ],
    [ [ [ 2, 0, 1 ], [ 0, 0, 1 ], [ 0, 1, 0 ], [ 3, 0, 0 ] ], -1 ],
    [ [ [ 2, 0, 1 ], [ 0, 0, 1 ], [ 0, 1, 0 ], [ 3, 1, 0 ] ], -2 ],
    [ [ [ 2, 0, 1 ], [ 0, 0, 1 ], [ 3, 1, 0 ], [ 0, 1, 0 ] ], 1 ],
    [ [ [ 2, 0, 1 ], [ 0, 1, 1 ], [ 0, 1, 0 ], [ 1, 0, 1 ] ], -2 ],
    [ [ [ 2, 0, 1 ], [ 0, 1, 1 ], [ 1, 1, 1 ], [ 0, 0, 1 ] ], -2 ],
    [ [ [ 2, 0, 1 ], [ 1, 0, 0 ], [ 1, 2, 1 ], [ 0, 0, 1 ] ], 2 ],
    [ [ [ 2, 0, 1 ], [ 1, 0, 0 ], [ 1, 2, 2 ], [ 0, 0, 1 ] ], 2 ],
    [ [ [ 2, 0, 1 ], [ 1, 0, 0 ], [ 3, 1, 0 ], [ 0, 1, 0 ] ], 1 ],
    [ [ [ 2, 0, 1 ], [ 1, 0, 0 ], [ 3, 1, 1 ], [ 1, 1, 2 ] ], -1 ],
    [ [ [ 2, 0, 1 ], [ 1, 1, 2 ], [ 2, 0, 1 ], [ 0, 1, 0 ] ], 1 ],
    [ [ [ 2, 0, 1 ], [ 3, 1, 0 ], [ 0, 0, 1 ], [ 3, 0, 0 ] ], 1 ],
    [ [ [ 2, 0, 1 ], [ 3, 1, 0 ], [ 0, 0, 1 ], [ 3, 1, 0 ] ], 1 ],
    [ [ [ 2, 0, 1 ], [ 3, 1, 1 ], [ 3, 1, 0 ], [ 1, 0, 0 ] ], 1 ],
    [ [ [ 2, 0, 1 ], [ 3, 2, 1 ], [ 1, 0, 0 ], [ 0, 0, 1 ] ], -2 ],
    [ [ [ 2, 0, 2 ], [ 0, 1, 0 ], [ 2, 0, 0 ], [ 1, 0, 0 ] ], 1 ],
    [ [ [ 2, 0, 2 ], [ 2, 0, 0 ], [ 0, 1, 0 ], [ 1, 0, 0 ] ], -1 ],
    [ [ [ 2, 0, 2 ], [ 2, 0, 0 ], [ 1, 0, 0 ], [ 0, 1, 0 ] ], 1 ],
    [ [ [ 2, 0, 2 ], [ 3, 1, 0 ], [ 1, 0, 0 ], [ 3, 1, 0 ] ], -2 ],
    [ [ [ 2, 1, 1 ], [ 0, 0, 1 ], [ 3, 0, 0 ], [ 0, 1, 0 ] ], 1 ],
    [ [ [ 3, 0, 0 ], [ 1, 0, 0 ], [ 3, 1, 0 ], [ 1, 0, 0 ] ], 1 ],
    [ [ [ 3, 0, 1 ], [ 0, 1, 0 ], [ 1, 0, 0 ], [ 0, 0, 1 ] ], 4 ],
    [ [ [ 3, 0, 1 ], [ 0, 1, 0 ], [ 1, 0, 0 ], [ 3, 1, 1 ] ], -1 ],
    [ [ [ 3, 0, 1 ], [ 0, 1, 0 ], [ 1, 0, 1 ], [ 0, 1, 0 ] ], 1 ],
    [ [ [ 3, 0, 1 ], [ 1, 0, 0 ], [ 0, 0, 1 ], [ 3, 1, 0 ] ], 1 ],
    [ [ [ 3, 0, 1 ], [ 1, 0, 0 ], [ 0, 0, 1 ], [ 3, 2, 0 ] ], 1 ],
    [ [ [ 3, 0, 1 ], [ 1, 0, 0 ], [ 3, 1, 0 ], [ 0, 1, 1 ] ], -1 ],
    [ [ [ 3, 0, 1 ], [ 1, 0, 0 ], [ 3, 1, 0 ], [ 1, 0, 1 ] ], -1 ],
    [ [ [ 3, 0, 1 ], [ 1, 0, 1 ], [ 3, 1, 0 ], [ 0, 1, 0 ] ], 1 ],
    [ [ [ 3, 0, 1 ], [ 1, 1, 0 ], [ 0, 0, 1 ], [ 0, 0, -1 ] ], 3 ],
    [ [ [ 3, 0, 1 ], [ 1, 1, 0 ], [ 0, 0, 1 ], [ 0, 1, 0 ] ], -1 ],
    [ [ [ 3, 0, 1 ], [ 1, 1, 1 ], [ 0, 0, -1 ], [ 0, 0, 1 ] ], -4 ],
    [ [ [ 3, 1, 0 ], [ 0, 0, 1 ], [ 3, 1, 0 ], [ 2, 0, 1 ] ], -1 ],
    [ [ [ 3, 1, 0 ], [ 1, 0, 0 ], [ 1, 0, 0 ], [ 0, 1, 0 ] ], -3 ],
    [ [ [ 3, 1, 0 ], [ 1, 0, 0 ], [ 1, 0, 1 ], [ 3, 1, 1 ] ], 1 ],
    [ [ [ 3, 1, 0 ], [ 1, 0, 0 ], [ 1, 1, 0 ], [ 0, 0, 2 ] ], -3 ],
    [ [ [ 3, 1, 0 ], [ 1, 0, 0 ], [ 3, 1, 1 ], [ 2, 0, 1 ] ], 1 ],
    [ [ [ 3, 1, 0 ], [ 2, 0, 0 ], [ 0, 0, 2 ], [ 0, 1, 0 ] ], -3 ],
    [ [ [ 3, 1, 0 ], [ 2, 0, 0 ], [ 0, 1, 0 ], [ 0, 0, 2 ] ], 3 ],
    [ [ [ 3, 1, 0 ], [ 2, 0, 1 ], [ 0, 0, 1 ], [ 0, 1, 0 ] ], 1 ],
    [ [ [ 3, 1, 1 ], [ 2, 0, 1 ], [ 0, 1, 0 ], [ 0, 0, 1 ] ], 1 ],
    [ [ [ 3, 1, 1 ], [ 3, 0, 0 ], [ 0, 1, 0 ], [ 1, 0, 0 ] ], -1 ],
    [ [ [ 3, 1, 1 ], [ 3, 0, 0 ], [ 1, 0, 0 ], [ 0, 1, 0 ] ], 1 ],
    [ [ [ 3, 1, 1 ], [ 3, 0, 0 ], [ 1, 0, 0 ], [ 3, 1, 0 ] ], 1 ],
];

if Length(Target) <> 76 then
    Error("FAIL: target must have 76 Laurent atoms");
fi;
if L1Sparse(Target) <> 76 then
    Error("FAIL: target must have l1 norm 76");
fi;
if Length(Certificate) <> 146 then
    Error("FAIL: certificate must contain 146 cocycle defects");
fi;
if Sum(Certificate, entry -> AbsInt(entry[2])) <> 245 then
    Error("FAIL: certificate coefficient l1 norm must be 245");
fi;
if SourceTarget <> Target then
    Error("FAIL: the six action coefficients do not reconstruct the recorded 76-atom target");
fi;

CertificateVector := [];
for entry in Certificate do
    CertificateVector := AddSparse(
        CertificateVector,
        ScaleSparse(CocycleDefect(
            entry[1][1], entry[1][2], entry[1][3], entry[1][4]
        ), entry[2])
    );
od;

Residual := AddSparse(Target, ScaleSparse(CertificateVector, -1));
if Length(Residual) <> 0 then
    Error("FAIL: nonzero Laurent residual; atoms = ", Length(Residual));
fi;

Print("Gamma4 terminal-root self-character certificate\n");
Print("Arithmetic: exact integers in a formal sparse Laurent lattice\n");
Print("External GAP packages loaded: none\n");
Print("Source action coefficients reconstructed: 6\n");
Print("Source scalars reconstructed: K_q, lambda_1, lambda_2\n");
Print("Source-to-recorded-target residual atoms: ",
      Length(AddSparse(SourceTarget, ScaleSparse(Target, -1))), "\n");
Print("Target Laurent atoms: ", Length(Target), "\n");
Print("Target l1 norm: ", L1Sparse(Target), "\n");
Print("Certificate cocycle defects: ", Length(Certificate), "\n");
Print("Certificate coefficient l1 norm: ",
      Sum(Certificate, entry -> AbsInt(entry[2])), "\n");
Print("Residual Laurent atoms: ", Length(Residual), "\n");
Print("PASS: Q is the stated product of normalized 3-cocycle defects, hence Q=1\n");

QUIT_GAP(0);
