#############################################################################
## Gamma4 X2 cocycle-reduction source-to-target certificate (core GAP only)
##
## This program verifies the calculation in
## Lemma `lem:Gamma4-X2-cocycle-reduction` without defining the source
## quotient from the six cocycle defects occurring at the end of the proof.
## It works in the abstract Gamma4 normal form epsilon^i h^j g^k and in an
## exact sparse Laurent exponent lattice.
#############################################################################

Gamma4X2Fail := function(arg)
    CallFuncList(PrintTo, Concatenation(["*errout*"], arg));
    PrintTo("*errout*", "\n");
    ErrorNoReturn("Gamma4 X2 certificate verification failed.");
end;

## Sparse Laurent monomials are lists [ [key, exponent], ... ].  Keys retain
## their structured arguments, so derived Phi_h and Phi^h atoms can be
## expanded independently after the source quotient has been constructed.
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
        return String(left[1]) < String(right[1]);
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

MultiplyMonomial := function(left, right)
    return AddSparse(left, right);
end;

DivideMonomial := function(left, right)
    return AddSparse(left, ScaleSparse(right, -1));
end;

Atom := function(key)
    return [[key, 1]];
end;

## Abstract Gamma4 normal forms epsilon^i h^j g^k.  The multiplication is
## derived only from
##   h*g=epsilon*g*h, g*epsilon=epsilon^-1*g,
##   h*epsilon=epsilon*h, epsilon^4=1.
IdentityElt := [0, 0, 0];

MultiplyElt := function(a, b)
    local first;
    if (a[3] mod 2) = 0 then
        first := a[1] + b[1];
    else
        first := a[1] - b[1] - b[2];
    fi;
    return [first mod 4, a[2] + b[2], a[3] + b[3]];
end;

InverseElt := function(a)
    if (a[3] mod 2) = 0 then
        return [(-a[1]) mod 4, -a[2], -a[3]];
    fi;
    return [(a[1] + a[2]) mod 4, -a[2], -a[3]];
end;

ConjugateElt := function(a, b)
    return MultiplyElt(MultiplyElt(a, b), InverseElt(a));
end;

PowerElt := function(a, n)
    local ans, i;
    if n < 0 then
        return PowerElt(InverseElt(a), -n);
    fi;
    ans := IdentityElt;
    for i in [1 .. n] do
        ans := MultiplyElt(ans, a);
    od;
    return ans;
end;

PhiAtom := function(a, b, c)
    if a = IdentityElt or b = IdentityElt or c = IdentityElt then
        return [];
    fi;
    return Atom(["P", a, b, c]);
end;

PhiLocalAtom := function(base, a, b)
    return Atom(["L", base, a, b]);
end;

PhiUpperAtom := function(base, a, b)
    return Atom(["U", base, a, b]);
end;

CharacterAtom := function(a)
    if a = IdentityElt then
        return [];
    fi;
    return Atom(["R", a]);
end;

## Definitions of Phi_base and Phi^base from Section 2.
ExpandPhiLocal := function(base, a, b)
    local ans;
    ans := PhiAtom(a, b, base);
    ans := AddSparse(ans,
        PhiAtom(ConjugateElt(MultiplyElt(a, b), base), a, b));
    ans := AddSparse(ans,
        ScaleSparse(PhiAtom(a, ConjugateElt(b, base), b), -1));
    return ans;
end;

ExpandPhiUpper := function(base, a, b)
    local ans;
    ans := PhiAtom(base, a, b);
    ans := AddSparse(ans,
        PhiAtom(ConjugateElt(base, a), ConjugateElt(base, b), base));
    ans := AddSparse(ans,
        ScaleSparse(PhiAtom(ConjugateElt(base, a), base, b), -1));
    return ans;
end;

ExpandDerivedAtoms := function(value)
    local ans, term, key, expansion;
    ans := [];
    for term in value do
        key := term[1];
        if key[1] = "P" then
            expansion := Atom(key);
        elif key[1] = "L" then
            expansion := ExpandPhiLocal(key[2], key[3], key[4]);
        elif key[1] = "U" then
            expansion := ExpandPhiUpper(key[2], key[3], key[4]);
        else
            Gamma4X2Fail("uneliminated atom in derived expansion: ", key);
        fi;
        ans := AddSparse(ans, ScaleSparse(expansion, term[2]));
    od;
    return ans;
end;

CocycleDefect := function(a, b, c, d)
    local ans;
    ans := PhiAtom(b, c, d);
    ans := AddSparse(ans, PhiAtom(a, MultiplyElt(b, c), d));
    ans := AddSparse(ans, PhiAtom(a, b, c));
    ans := AddSparse(ans,
        ScaleSparse(PhiAtom(MultiplyElt(a, b), c, d), -1));
    ans := AddSparse(ans,
        ScaleSparse(PhiAtom(a, b, MultiplyElt(c, d)), -1));
    return ans;
end;

## The quotient in the local-cocycle coherence identity.
LocalCoherence := function(x, y, z, w)
    local ans;
    ans := PhiLocalAtom(x, z, w);
    ans := AddSparse(ans, PhiLocalAtom(x, y, MultiplyElt(z, w)));
    ans := AddSparse(ans,
        ScaleSparse(PhiLocalAtom(ConjugateElt(w, x), y, z), -1));
    ans := AddSparse(ans,
        ScaleSparse(PhiLocalAtom(x, MultiplyElt(y, z), w), -1));
    return ans;
end;

## Independent four-defect certificate for local coherence, obtained by
## expanding the two factors A and B in the proof of that lemma.
LocalCoherenceDefects := function(x, y, z, w)
    local v, t, u, ans;
    v := ConjugateElt(w, x);
    t := ConjugateElt(z, v);
    u := ConjugateElt(y, t);
    if MultiplyElt(v, w) <> MultiplyElt(w, x) or
       MultiplyElt(t, z) <> MultiplyElt(z, v) or
       MultiplyElt(u, y) <> MultiplyElt(y, t) then
        Gamma4X2Fail("conjugation identity failed in local coherence");
    fi;
    ans := CocycleDefect(y, z, w, x);
    ans := AddSparse(ans, ScaleSparse(CocycleDefect(y, z, v, w), -1));
    ans := AddSparse(ans, CocycleDefect(y, t, z, w));
    ans := AddSparse(ans, ScaleSparse(CocycleDefect(u, y, z, w), -1));
    return ans;
end;

CheckLocalCoherenceCertificate := function(x, y, z, w)
    local expanded, defects, residual;
    expanded := ExpandDerivedAtoms(LocalCoherence(x, y, z, w));
    defects := LocalCoherenceDefects(x, y, z, w);
    residual := AddSparse(expanded, ScaleSparse(defects, -1));
    if Length(residual) <> 0 then
        Gamma4X2Fail("local-coherence defect certificate has residual ",
                     residual);
    fi;
    return Length(expanded);
end;

CountKind := function(value, kind)
    return Length(Filtered(value, term -> term[1][1] = kind));
end;

## Marked Gamma4 elements.
Epsilon := [1, 0, 0];
H := [0, 1, 0];
G := [0, 0, 1];
R := InverseElt(Epsilon);
U := PowerElt(Epsilon, 2);
RH := MultiplyElt(R, H);
EpsilonG := MultiplyElt(Epsilon, G);
UH := MultiplyElt(U, H);

if R <> [3, 0, 0] or U <> [2, 0, 0] or
   MultiplyElt(H, G) <> MultiplyElt(Epsilon, MultiplyElt(G, H)) or
   MultiplyElt(G, R) <> EpsilonG or
   MultiplyElt(G, U) <> MultiplyElt(U, G) or
   ConjugateElt(H, Epsilon) <> Epsilon then
    Gamma4X2Fail("Gamma4 normal-form regression failed");
fi;

## Construct alpha, kappa_X, and Xi_X independently from their displayed
## definitions.  At this stage the source contains rho-character atoms and
## no cocycle defects have been introduced.
AlphaX := DivideMonomial(CharacterAtom(R), PhiLocalAtom(H, G, R));

KappaX := CharacterAtom(U);
KappaX := DivideMonomial(KappaX, PhiLocalAtom(H, U, H));
KappaX := DivideMonomial(KappaX, PhiLocalAtom(H, G, UH));

XiX := AlphaX;
XiX := MultiplyMonomial(XiX, PhiLocalAtom(H, H, G));
XiX := MultiplyMonomial(XiX, CharacterAtom(RH));
XiX := DivideMonomial(XiX, CharacterAtom(H));
XiX := DivideMonomial(XiX, PhiUpperAtom(H, H, G));
XiX := DivideMonomial(XiX, PhiLocalAtom(H, G, RH));
XiX := MultiplyMonomial(XiX, PhiAtom(H, RH, EpsilonG));
XiX := DivideMonomial(XiX, PhiAtom(RH, H, EpsilonG));

SourceRaw := DivideMonomial(XiX, KappaX);
if CountKind(SourceRaw, "R") <> 4 then
    Gamma4X2Fail("raw Xi_X/kappa_X must contain four rho atoms");
fi;

## For a Phi_h-projective character rho, this is the exponent vector of
## rho(a)rho(b)/(Phi_h(a,b)rho(ab)).
ProjectiveDefect := function(a, b)
    local ans;
    ans := CharacterAtom(a);
    ans := AddSparse(ans, CharacterAtom(b));
    ans := AddSparse(ans, ScaleSparse(PhiLocalAtom(H, a, b), -1));
    ans := AddSparse(ans,
        ScaleSparse(CharacterAtom(MultiplyElt(a, b)), -1));
    return ans;
end;

ProjectiveRR := ProjectiveDefect(R, R);
ProjectiveRH := ProjectiveDefect(R, H);
if MultiplyElt(R, R) <> U or MultiplyElt(R, H) <> RH then
    Gamma4X2Fail("projective-identity products have wrong normal forms");
fi;

## Add the second projective identity and subtract the first.  Both are
## zero for rho.  This eliminates all rho atoms in the raw source.
SourceAfterProjective := AddSparse(SourceRaw, ProjectiveRH);
SourceAfterProjective := AddSparse(SourceAfterProjective,
                                   ScaleSparse(ProjectiveRR, -1));
if CountKind(SourceAfterProjective, "R") <> 0 then
    Gamma4X2Fail("projective substitution did not eliminate rho atoms");
fi;

## Independently encode the quotient displayed immediately after the two
## instances of eq:YD-projective-action and compare it with the transformed
## source.
DisplayedAfterProjective := PhiLocalAtom(H, R, R);
DisplayedAfterProjective := AddSparse(DisplayedAfterProjective,
    PhiLocalAtom(H, U, H));
DisplayedAfterProjective := AddSparse(DisplayedAfterProjective,
    PhiLocalAtom(H, G, UH));
DisplayedAfterProjective := AddSparse(DisplayedAfterProjective,
    PhiAtom(H, RH, EpsilonG));
DisplayedAfterProjective := AddSparse(DisplayedAfterProjective,
    PhiLocalAtom(H, H, G));
DisplayedAfterProjective := AddSparse(DisplayedAfterProjective,
    ScaleSparse(PhiLocalAtom(H, G, R), -1));
DisplayedAfterProjective := AddSparse(DisplayedAfterProjective,
    ScaleSparse(PhiLocalAtom(H, G, RH), -1));
DisplayedAfterProjective := AddSparse(DisplayedAfterProjective,
    ScaleSparse(PhiLocalAtom(H, R, H), -1));
DisplayedAfterProjective := AddSparse(DisplayedAfterProjective,
    ScaleSparse(PhiAtom(RH, H, EpsilonG), -1));
DisplayedAfterProjective := AddSparse(DisplayedAfterProjective,
    ScaleSparse(PhiUpperAtom(H, H, G), -1));

ProjectiveResidual := AddSparse(SourceAfterProjective,
    ScaleSparse(DisplayedAfterProjective, -1));
if Length(ProjectiveResidual) <> 0 then
    Gamma4X2Fail("source does not reduce to the displayed quotient");
fi;

## The three precise local-coherence instances used in the proof.  First
## verify each one independently after expanding Phi_h into Phi-atoms.
Coherence1 := LocalCoherence(H, G, R, H);
Coherence2 := LocalCoherence(H, G, U, H);
Coherence3 := LocalCoherence(H, G, R, R);
CoherenceAtomCounts := [
    CheckLocalCoherenceCertificate(H, G, R, H),
    CheckLocalCoherenceCertificate(H, G, U, H),
    CheckLocalCoherenceCertificate(H, G, R, R)
];

## Apply the three equalities in exactly the orientations used in the proof.
AfterCoherence := AddSparse(SourceAfterProjective, Coherence1);
AfterCoherence := AddSparse(AfterCoherence, ScaleSparse(Coherence2, -1));
AfterCoherence := AddSparse(AfterCoherence, ScaleSparse(Coherence3, -1));
ExpandedAfterCoherence := ExpandDerivedAtoms(AfterCoherence);

## Only now construct the six cocycle defects printed at the end of the
## lemma.  This target is independent of Alpha, KappaX, XiX, and the two
## projective substitutions above.
TargetSixDefects := CocycleDefect(H, G, H, R);
TargetSixDefects := AddSparse(TargetSixDefects,
    CocycleDefect(RH, EpsilonG, R, H));
TargetSixDefects := AddSparse(TargetSixDefects,
    CocycleDefect(RH, H, G, R));
TargetSixDefects := AddSparse(TargetSixDefects,
    ScaleSparse(CocycleDefect(H, G, R, H), -1));
TargetSixDefects := AddSparse(TargetSixDefects,
    ScaleSparse(CocycleDefect(RH, EpsilonG, H, R), -1));
TargetSixDefects := AddSparse(TargetSixDefects,
    ScaleSparse(CocycleDefect(H, RH, G, R), -1));

FinalResidual := AddSparse(ExpandedAfterCoherence,
                           ScaleSparse(TargetSixDefects, -1));
if Length(FinalResidual) <> 0 then
    Gamma4X2Fail("six-defect source-to-target residual is nonzero: ",
                 FinalResidual);
fi;

Print("Gamma4 X2 cocycle-reduction source-to-target certificate\n");
Print("GAP version: ", GAPInfo.Version, "\n");
Print("Arithmetic: exact integers in a formal sparse Laurent lattice\n");
Print("External GAP packages required: none\n");
Print("Raw Xi_X/kappa_X rho atoms: ", CountKind(SourceRaw, "R"), "\n");
Print("Projective multiplication identities used: 2\n");
Print("Projective source-to-displayed residual: ",
      Length(ProjectiveResidual), "\n");
Print("Local-coherence instances expanded and certified: 3\n");
Print("Expanded local-coherence atom counts: ",
      CoherenceAtomCounts, "\n");
Print("Final target cocycle defects: 6\n");
Print("Expanded six-defect Laurent atoms after cancellation: ",
      Length(TargetSixDefects), "\n");
Print("Final source-to-six-defect residual: ", Length(FinalResidual), "\n");
Print("PASS: Xi_X/kappa_X equals the stated six-defect product and hence 1\n");
