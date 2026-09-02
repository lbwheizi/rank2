#############################################################################
## Gamma4 Y3 two-coefficient source-to-target cancellation certificate
## (core GAP only)
##
## PURPOSE
##   Reconstruct beta, B1, nu1, and mu1 from the recursive source and actual
##   epsilon-action coordinates, proving beta*B1*nu1*mu1=-sigma(g) before
##   specializing sigma(g)=-1.  Independently reconstruct the common
##   coordinate used in the Y_3 calculation directly
##   from the corrected recursive operators phi_1, phi_2, phi_3.  The program
##   extracts its two Laurent monomials, converts their ratio with the stated
##   one-dimensional projective-character relation, and verifies that the
##   result is exactly the eight-atom quotient
##
##       Xi3 = D(h,g,ug,g) / D(epsilon*g,h,ug,g) = 1.
##
##   Here D is the normalized 3-cocycle defect.  Thus this file closes both
##   links which are logically separate in the paper:
##
##       corrected phi recursion -> common-coordinate quotient ->
##       eight-atom target -> two cocycle defects -> zero residual.
##
## SCOPE
##   The induced objects use the abstract Gamma4 normal forms
##   epsilon^i h^j g^k and the fixed sections
##
##       r_h=1, r_{epsilon^-1 h}=g,
##       r_g=1, r_{epsilon g}=h, r_{ug}=epsilon,
##       r_{epsilon^-1 g}=epsilon*h.
##
##   No quotient-specific order relation, external GAP package, sampled
##   cocycle, or floating-point arithmetic is used.  Coefficients are exact
##   integers in a formal sparse Laurent lattice.
##
## FAIL-FAST BATCH REPRODUCTION (from this directory)
##   gap --bare -A -r -q --quitonbreak --norepl gamma4_xi3_certificate.g
#############################################################################

Gamma4Fail := function(arg)
    CallFuncList(PrintTo, Concatenation(["*errout*"], arg));
    PrintTo("*errout*", "\n");
    ErrorNoReturn("Certificate verification failed.");
end;

## Sparse Laurent monomials are sorted lists [ [key, exponent], ... ].
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

MultiplyMonomial := function(left, right)
    return AddSparse(left, right);
end;

DivideMonomial := function(left, right)
    return AddSparse(left, ScaleSparse(right, -1));
end;

PowerMonomial := function(value, exponent)
    return ScaleSparse(value, exponent);
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

## Abstract Gamma4 normal forms epsilon^i h^j g^k.
IdentityElt := [0, 0, 0];

MultiplyElt := function(a, b)
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
    local key;
    if a = IdentityElt or b = IdentityElt or c = IdentityElt then
        return [];
    fi;
    key := Concatenation("P(", String(a), ",", String(b), ",", String(c), ")");
    return [[key, 1]];
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
    ans := AddSparse(ans,
        PhiAtom(ConjugateElt(MultiplyElt(a, b), y), a, b));
    ans := AddSparse(ans,
        ScaleSparse(PhiAtom(a, ConjugateElt(b, y), b), -1));
    return ans;
end;

PhiUpper := function(a, x, y)
    local ans;
    ans := PhiAtom(a, x, y);
    ans := AddSparse(ans,
        PhiAtom(ConjugateElt(a, x), ConjugateElt(a, y), a));
    ans := AddSparse(ans,
        ScaleSparse(PhiAtom(ConjugateElt(a, x), a, y), -1));
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

## Sparse Laurent polynomials are lists [ [monomial, integer], ... ].
AddPolyTerm := function(poly, monomial, coefficient)
    local i;
    if coefficient = 0 then
        return;
    fi;
    for i in [1 .. Length(poly)] do
        if poly[i][1] = monomial then
            poly[i][2] := poly[i][2] + coefficient;
            if poly[i][2] = 0 then
                Remove(poly, i);
            fi;
            return;
        fi;
    od;
    Add(poly, [monomial, coefficient]);
end;

CanonicalPoly := function(poly)
    local ans, term;
    ans := [];
    for term in poly do
        AddPolyTerm(ans, term[1], term[2]);
    od;
    Sort(ans, function(left, right)
        return String(left[1]) < String(right[1]);
    end);
    return ans;
end;

MonomialPoly := function(monomial)
    return [[monomial, 1]];
end;

AddPoly := function(left, right, scale)
    local ans, term;
    ans := List(left, term -> [List(term[1], ShallowCopy), term[2]]);
    for term in right do
        AddPolyTerm(ans, term[1], scale * term[2]);
    od;
    return CanonicalPoly(ans);
end;

MultiplyPolyMonomial := function(poly, monomial)
    return CanonicalPoly(List(poly,
        term -> [MultiplyMonomial(term[1], monomial), term[2]]));
end;

ZeroPolyVector := function(n)
    return List([1 .. n], i -> []);
end;

BasisPolyVector := function(n, position)
    local ans;
    ans := ZeroPolyVector(n);
    ans[position] := MonomialPoly([]);
    return ans;
end;

AddPolyVector := function(left, right, scale)
    return List([1 .. Length(left)],
        i -> AddPoly(left[i], right[i], scale));
end;

AddToVectorPosition := function(vector, position, poly)
    vector[position] := AddPoly(vector[position], poly, 1);
end;

PositionDegree := function(degrees, degree)
    local position;
    position := Position(degrees, degree);
    if position = fail then
        Gamma4Fail("degree not found in induced support: ", degree);
    fi;
    return position;
end;

MakeInduced := function(base, name, degrees, representatives)
    return rec(
        kind := "induced",
        base := base,
        name := name,
        degrees := degrees,
        representatives := representatives,
        dim := Length(degrees)
    );
end;

MakeTensor := function(left, right)
    local degrees, x, y;
    degrees := [];
    for x in left.degrees do
        for y in right.degrees do
            Add(degrees, MultiplyElt(x, y));
        od;
    od;
    return rec(
        kind := "tensor",
        left := left,
        right := right,
        degrees := degrees,
        dim := left.dim * right.dim
    );
end;

## Return [row, Laurent monomial] for the unique nonzero entry in an
## action column.  All objects used here have monomial action matrices.
ActionColumn := function(obj, a, column)
    local y, newY, row, repY, repNew, centralizerElt, gamma,
          rightDim, i, j, leftAction, rightAction, correction;
    if obj.kind = "induced" then
        y := obj.degrees[column];
        newY := ConjugateElt(a, y);
        row := PositionDegree(obj.degrees, newY);
        repY := obj.representatives[column];
        repNew := obj.representatives[row];
        centralizerElt := MultiplyElt(
            InverseElt(repNew), MultiplyElt(a, repY));
        gamma := PhiAtom(a, repY, obj.base);
        gamma := AddSparse(gamma, PhiAtom(newY, a, repY));
        gamma := AddSparse(gamma, PhiAtom(repNew, obj.base, centralizerElt));
        gamma := AddSparse(gamma, ScaleSparse(PhiAtom(a, y, repY), -1));
        gamma := AddSparse(gamma,
            ScaleSparse(PhiAtom(repNew, centralizerElt, obj.base), -1));
        gamma := AddSparse(gamma,
            ScaleSparse(PhiAtom(newY, repNew, centralizerElt), -1));
        gamma := MultiplyMonomial(gamma,
            CharacterAtom(obj.name, centralizerElt));
        return [row, gamma];
    fi;

    rightDim := obj.right.dim;
    i := QuoInt(column - 1, rightDim) + 1;
    j := ((column - 1) mod rightDim) + 1;
    leftAction := ActionColumn(obj.left, a, i);
    rightAction := ActionColumn(obj.right, a, j);
    correction := PhiUpper(a, obj.left.degrees[i], obj.right.degrees[j]);
    row := (leftAction[1] - 1) * rightDim + rightAction[1];
    return [row, MultiplyMonomial(correction,
        MultiplyMonomial(leftAction[2], rightAction[2]))];
end;

ApplyActionVector := function(obj, a, vector)
    local ans, column, action;
    ans := ZeroPolyVector(obj.dim);
    for column in [1 .. obj.dim] do
        if Length(vector[column]) > 0 then
            action := ActionColumn(obj, a, column);
            AddToVectorPosition(ans, action[1],
                MultiplyPolyMonomial(vector[column], action[2]));
        fi;
    od;
    return ans;
end;

SignedMonomialRatio := function(numerator, denominator)
    if Length(numerator) <> 1 or Length(denominator) <> 1 then
        Gamma4Fail("FAIL: a signed monomial ratio has a non-monomial entry");
    fi;
    if AbsInt(numerator[1][2]) <> 1 or AbsInt(denominator[1][2]) <> 1 then
        Gamma4Fail("FAIL: a signed monomial ratio has a non-unit coefficient");
    fi;
    return rec(
        monomial := DivideMonomial(numerator[1][1], denominator[1][1]),
        sign := numerator[1][2] * denominator[1][2]
    );
end;

MultiplySignedMonomials := function(left, right)
    return rec(
        monomial := MultiplyMonomial(left.monomial, right.monomial),
        sign := left.sign * right.sign
    );
end;

DivideSignedMonomials := function(left, right)
    return rec(
        monomial := DivideMonomial(left.monomial, right.monomial),
        sign := left.sign * right.sign
    );
end;

SpecializeSignedCharacterMinusOne := function(value, key)
    local exponent, answer;
    exponent := ExponentOf(value.monomial, key);
    answer := rec(
        monomial := AddSparse(value.monomial, [[key, -exponent]]),
        sign := value.sign
    );
    if (exponent mod 2) <> 0 then
        answer.sign := -answer.sign;
    fi;
    return answer;
end;

ApplyBraid := function(left, right, vector)
    local ans, column, i, j, action, row, rightDim;
    ans := ZeroPolyVector(right.dim * left.dim);
    rightDim := right.dim;
    for column in [1 .. Length(vector)] do
        if Length(vector[column]) > 0 then
            i := QuoInt(column - 1, rightDim) + 1;
            j := ((column - 1) mod rightDim) + 1;
            action := ActionColumn(right, left.degrees[i], j);
            row := (action[1] - 1) * left.dim + i;
            AddToVectorPosition(ans, row,
                MultiplyPolyMonomial(vector[column], action[2]));
        fi;
    od;
    return ans;
end;

Phi1Apply := function(acting, target, vector)
    local double;
    double := ApplyBraid(target, acting,
        ApplyBraid(acting, target, vector));
    return AddPolyVector(vector, double, -1);
end;

## Corrected c_12 = a (c_AA tensor id) a^{-1}.  In the convention of the
## paper a^{-1} contributes Phi(x,y,z) on the input and a contributes its
## inverse on the braided output.
C12Apply := function(acting, earlier, vector)
    local ans, column, block, i, remainder, j, k, x, y, z,
          braidAction, outputFirst, row, coefficient;
    ans := ZeroPolyVector(acting.dim * acting.dim * earlier.dim);
    block := acting.dim * earlier.dim;
    for column in [1 .. Length(vector)] do
        if Length(vector[column]) > 0 then
            i := QuoInt(column - 1, block) + 1;
            remainder := (column - 1) mod block;
            j := QuoInt(remainder, earlier.dim) + 1;
            k := (remainder mod earlier.dim) + 1;
            x := acting.degrees[i];
            y := acting.degrees[j];
            z := earlier.degrees[k];
            braidAction := ActionColumn(acting, x, j);
            outputFirst := acting.degrees[braidAction[1]];
            coefficient := MultiplyMonomial(PhiAtom(x, y, z),
                braidAction[2]);
            coefficient := DivideMonomial(coefficient,
                PhiAtom(outputFirst, x, z));
            row := (braidAction[1] - 1) * block
                 + (i - 1) * earlier.dim + k;
            AddToVectorPosition(ans, row,
                MultiplyPolyMonomial(vector[column], coefficient));
        fi;
    od;
    return ans;
end;

ApplyOuterPrevious := function(acting, previous, previousApply, vector)
    local ans, outer, start, slice, image, j;
    ans := ZeroPolyVector(acting.dim * previous.dim);
    for outer in [1 .. acting.dim] do
        start := (outer - 1) * previous.dim;
        slice := vector{[start + 1 .. start + previous.dim]};
        image := previousApply(slice);
        for j in [1 .. previous.dim] do
            ans[start + j] := image[j];
        od;
    od;
    return ans;
end;

HigherPhiApply := function(acting, previous, earlier, previousApply, vector)
    local double, correctedBraid, correction;
    double := ApplyBraid(previous, acting,
        ApplyBraid(acting, previous, vector));
    correctedBraid := C12Apply(acting, earlier, vector);
    correction := ApplyOuterPrevious(acting, previous,
        previousApply, correctedBraid);
    return AddPolyVector(AddPolyVector(vector, double, -1), correction, 1);
end;

## Marked Gamma4 data and induced objects.
Epsilon := [1, 0, 0];
H := [0, 1, 0];
G := [0, 0, 1];
R := InverseElt(Epsilon);
U := PowerElt(Epsilon, 2);
T := PowerElt(G, 2);
RHElt := MultiplyElt(R, H);
EGElt := MultiplyElt(Epsilon, G);
UGElt := MultiplyElt(U, G);
RGElt := MultiplyElt(R, G);
UTElt := MultiplyElt(U, T);

VObj := MakeInduced(H, "R",
    [H, RHElt],
    [IdentityElt, G]);
WObj := MakeInduced(G, "S",
    [G, EGElt, UGElt, RGElt],
    [IdentityElt, H, Epsilon, MultiplyElt(Epsilon, H)]);
D1Obj := MakeTensor(WObj, VObj);
D2Obj := MakeTensor(WObj, D1Obj);

Phi1WVonVector := function(vector)
    return Phi1Apply(WObj, VObj, vector);
end;

Phi2WVonVector := function(vector)
    return HigherPhiApply(WObj, D1Obj, VObj,
        Phi1WVonVector, vector);
end;

## y_1 = phi_1(w_{epsilon g} tensor v_h).
Y1InputPosition := (PositionDegree(WObj.degrees, EGElt) - 1) * VObj.dim
                 + PositionDegree(VObj.degrees, H);
Y1Vector := Phi1WVonVector(BasisPolyVector(D1Obj.dim, Y1InputPosition));

## y_2 = phi_2(w_{ug} tensor y_1).
Y2Input := ZeroPolyVector(D2Obj.dim);
for IndexJ in [1 .. D1Obj.dim] do
    Y2Input[(PositionDegree(WObj.degrees, UGElt) - 1) * D1Obj.dim + IndexJ]
        := Y1Vector[IndexJ];
od;
Y2Vector := Phi2WVonVector(Y2Input);
Y2SupportPositions := Filtered([1 .. D2Obj.dim],
    i -> Length(Y2Vector[i]) > 0);
if Y2SupportPositions <> [6, 7, 9, 16, 18, 19, 28, 29] then
    Gamma4Fail("FAIL: corrected phi_2 has unexpected nonzero coordinates: ",
          Y2SupportPositions);
fi;

#############################################################################
## First coefficient in equation (C.9).
##
## Recover beta, B_1, nu_1, and mu_1 from the recursive vectors, the
## actual c_12 image, and actual action columns.  In particular, neither
## comparison identity at the end of Lemma C.5 is supplied as input.
#############################################################################

WIndexG := PositionDegree(WObj.degrees, G);
WIndexUG := PositionDegree(WObj.degrees, UGElt);
VIndexRH := PositionDegree(VObj.degrees, RHElt);

BetaPosition := (WIndexG - 1) * VObj.dim + VIndexRH;
if Length(Y1Vector[BetaPosition]) <> 1 or
   Y1Vector[BetaPosition][1][2] <> -1 then
    Gamma4Fail("FAIL: cannot recover beta from y_1");
fi;
BetaSigned := rec(
    monomial := Y1Vector[BetaPosition][1][1],
    sign := 1
);

## Extract B_1 from c_12(w_2 tensor (w tensor v_1)).  The c_12 output
## contains the actual scalar of (epsilon^2 g) acting on w; dividing by
## that action scalar leaves B_1.
C12SeedPosition := (WIndexUG - 1) * D1Obj.dim
                 + (WIndexG - 1) * VObj.dim + VIndexRH;
C12Seed := BasisPolyVector(D2Obj.dim, C12SeedPosition);
C12Image := C12Apply(WObj, VObj, C12Seed);
C12Support := Filtered([1 .. D2Obj.dim],
    i -> Length(C12Image[i]) > 0);
R1Action := ActionColumn(WObj, UGElt, WIndexG);
if Length(C12Support) <> 1 or R1Action[1] <> WIndexG then
    Gamma4Fail("FAIL: the c_12 source for B_1 has the wrong output line");
fi;
ExpectedC12Position := (R1Action[1] - 1) * D1Obj.dim
                     + (WIndexUG - 1) * VObj.dim + VIndexRH;
if C12Support[1] <> ExpectedC12Position then
    Gamma4Fail("FAIL: the c_12 source for B_1 has the wrong tensor position");
fi;
B1Signed := SignedMonomialRatio(C12Image[ExpectedC12Position],
    MonomialPoly(R1Action[2]));
B1Direct := DivideMonomial(PhiAtom(UGElt, G, RHElt),
    PhiAtom(ConjugateElt(UGElt, G), UGElt, RHElt));
if B1Signed.sign <> 1 or B1Signed.monomial <> B1Direct then
    Gamma4Fail("FAIL: c_12 does not reconstruct the stated B_1 quotient");
fi;

## r_1=(epsilon^2 g).w and g.r_1=nu_1 w.
GActionOnW := ActionColumn(WObj, G, WIndexG);
if R1Action[1] <> WIndexG or GActionOnW[1] <> WIndexG then
    Gamma4Fail("FAIL: r_1 or g.r_1 leaves the w-line");
fi;
Nu1Signed := rec(
    monomial := MultiplyMonomial(R1Action[2], GActionOnW[2]),
    sign := 1
);

## p_1=phi_1(w_2 tensor v_1), followed by phi_2(w tensor p_1).
P1InputPosition := (WIndexUG - 1) * VObj.dim + VIndexRH;
P1Vector := Phi1WVonVector(BasisPolyVector(D1Obj.dim, P1InputPosition));
WP1Input := ZeroPolyVector(D2Obj.dim);
for IndexJ in [1 .. D1Obj.dim] do
    WP1Input[(WIndexG - 1) * D1Obj.dim + IndexJ] := P1Vector[IndexJ];
od;
Phi2WP1 := Phi2WVonVector(WP1Input);
if Filtered([1 .. D2Obj.dim], i -> Length(Phi2WP1[i]) > 0) <>
   Y2SupportPositions then
    Gamma4Fail("FAIL: phi_2(w tensor p_1) has the wrong coordinate support");
fi;

## Use the nonzero line of multidegree (g,epsilon^2 g,epsilon^-1 h).
## The paper proves that this component of Y_2 is one-dimensional, so this
## actual coordinate determines the scalar mu_1.
Mu1Position := (WIndexG - 1) * D1Obj.dim
             + (WIndexUG - 1) * VObj.dim + VIndexRH;
Mu1Signed := SignedMonomialRatio(Phi2WP1[Mu1Position],
                                 Y2Vector[Mu1Position]);

FirstCoefficientProduct := MultiplySignedMonomials(
    MultiplySignedMonomials(BetaSigned, B1Signed),
    MultiplySignedMonomials(Nu1Signed, Mu1Signed));
SgKey := CharacterAtom("S", G)[1][1];
if FirstCoefficientProduct.sign <> -1 or
   FirstCoefficientProduct.monomial <> CharacterAtom("S", G) then
    Gamma4Fail("FAIL: beta*B1*nu1*mu1 is not -sigma(g)");
fi;
FirstCoefficientSpecialized := SpecializeSignedCharacterMinusOne(
    FirstCoefficientProduct, SgKey);
if FirstCoefficientSpecialized.sign <> 1 or
   Length(FirstCoefficientSpecialized.monomial) <> 0 then
    Gamma4Fail("FAIL: beta*B1*nu1*mu1 is not 1 at sigma(g)=-1");
fi;

## Independently reconstruct the two epsilon-action paths used in the
## manuscript.  The first path acts on w_2 tensor y_1 and then applies
## phi_2; the second acts directly on z_1.
EpsilonSeedAction := ApplyActionVector(D2Obj, Epsilon, Y2Input);
WP1Support := Filtered([1 .. D2Obj.dim], i -> Length(WP1Input[i]) > 0);
if Filtered([1 .. D2Obj.dim], i -> Length(EpsilonSeedAction[i]) > 0) <>
   WP1Support then
    Gamma4Fail("FAIL: epsilon.(w_2 tensor y_1) has the wrong support");
fi;
Chi1Signed := SignedMonomialRatio(EpsilonSeedAction[WP1Support[1]],
                                  WP1Input[WP1Support[1]]);

EpsilonZ1 := ApplyActionVector(D2Obj, Epsilon, Y2Vector);
if Filtered([1 .. D2Obj.dim], i -> Length(EpsilonZ1[i]) > 0) <>
   Y2SupportPositions then
    Gamma4Fail("FAIL: epsilon.z_1 has the wrong coordinate support");
fi;
Tau2EpsilonSigned := SignedMonomialRatio(EpsilonZ1[Mu1Position],
                                         Y2Vector[Mu1Position]);
FirstEpsilonPathResidual := DivideSignedMonomials(
    MultiplySignedMonomials(Chi1Signed, Mu1Signed),
    Tau2EpsilonSigned);
if FirstEpsilonPathResidual.sign <> 1 or
   Length(FirstEpsilonPathResidual.monomial) <> 0 then
    Gamma4Fail("FAIL: the first epsilon-action path has a nonzero residual");
fi;

SecondEpsilonPathResidual := DivideSignedMonomials(
    Chi1Signed,
    MultiplySignedMonomials(Tau2EpsilonSigned,
        MultiplySignedMonomials(
            MultiplySignedMonomials(BetaSigned, B1Signed), Nu1Signed)));
SecondEpsilonPathResidual := SpecializeSignedCharacterMinusOne(
    SecondEpsilonPathResidual, SgKey);
if SecondEpsilonPathResidual.sign <> 1 or
   Length(SecondEpsilonPathResidual.monomial) <> 0 then
    Gamma4Fail("FAIL: the second epsilon-action path has a nonzero residual");
fi;

## y_3 = phi_3(w_g tensor y_2).
Y3Input := ZeroPolyVector(WObj.dim * D2Obj.dim);
for IndexJ in [1 .. D2Obj.dim] do
    Y3Input[(PositionDegree(WObj.degrees, G) - 1) * D2Obj.dim + IndexJ]
        := Y2Vector[IndexJ];
od;
Y3Vector := HigherPhiApply(WObj, D2Obj, D1Obj,
    Phi2WVonVector, Y3Input);

## Check the first coefficient once more inside the full phi_3 output.  The
## selected coordinate is w_g tensored with the same nonzero Y_2 line used
## to determine mu_1.  Its two source terms have quotient sigma(g), hence
## cancel after sigma(g)=-1.
FirstY3Position := (WIndexG - 1) * D2Obj.dim + Mu1Position;
FirstY3Coordinate := Y3Vector[FirstY3Position];
if Length(FirstY3Coordinate) <> 2 or
   Set(List(FirstY3Coordinate, term -> term[2])) <> [-1] then
    Gamma4Fail("FAIL: the first Y_3 coordinate does not have two source terms");
fi;
FirstY3WithSigma := Filtered(FirstY3Coordinate,
    term -> ExponentOf(term[1], SgKey) = 1);
FirstY3WithoutSigma := Filtered(FirstY3Coordinate,
    term -> ExponentOf(term[1], SgKey) = 0);
if Length(FirstY3WithSigma) <> 1 or Length(FirstY3WithoutSigma) <> 1 or
   DivideMonomial(FirstY3WithSigma[1][1],
                  FirstY3WithoutSigma[1][1]) <> CharacterAtom("S", G) then
    Gamma4Fail("FAIL: the two first-coefficient source terms have the wrong quotient");
fi;
FirstY3Specialized := [];
for FirstY3Term in FirstY3Coordinate do
    FirstY3Signed := SpecializeSignedCharacterMinusOne(
        rec(monomial := FirstY3Term[1], sign := FirstY3Term[2]), SgKey);
    AddPolyTerm(FirstY3Specialized, FirstY3Signed.monomial,
                FirstY3Signed.sign);
od;
if Length(FirstY3Specialized) <> 0 then
    Gamma4Fail("FAIL: the first Y_3 coordinate does not cancel at sigma(g)=-1");
fi;

## The common line used for C_- and C_+ in the proof has tensor factors
##
##   (epsilon*g, epsilon*g, epsilon^2*g, epsilon^-1*h).
##
## Compute its ambient position from those four indices instead of taking
## the position from the recorded coefficient calculation.
CommonTensorIndices := [
    PositionDegree(WObj.degrees, EGElt),
    PositionDegree(WObj.degrees, EGElt),
    PositionDegree(WObj.degrees, UGElt),
    PositionDegree(VObj.degrees, RHElt)
];
CommonPosition :=
    (CommonTensorIndices[1] - 1) * D2Obj.dim
  + (CommonTensorIndices[2] - 1) * D1Obj.dim
  + (CommonTensorIndices[3] - 1) * VObj.dim
  + CommonTensorIndices[4];
CommonTensorFactors := [
    WObj.degrees[CommonTensorIndices[1]],
    WObj.degrees[CommonTensorIndices[2]],
    WObj.degrees[CommonTensorIndices[3]],
    VObj.degrees[CommonTensorIndices[4]]
];
if CommonPosition <> 46 or
   CommonTensorFactors <> [EGElt, EGElt, UGElt, RHElt] then
    Gamma4Fail("FAIL: common Y_3 tensor position or multidegree is wrong");
fi;
CommonTotalDegree := IdentityElt;
for CommonFactor in CommonTensorFactors do
    CommonTotalDegree := MultiplyElt(CommonTotalDegree, CommonFactor);
od;
if CommonTotalDegree <> MultiplyElt(MultiplyElt(U, H), PowerElt(G, 3)) then
    Gamma4Fail("FAIL: common Y_3 total degree is not epsilon^2*h*g^3");
fi;

CommonCoordinate := Y3Vector[CommonPosition];
if Length(CommonCoordinate) <> 2 then
    Gamma4Fail("FAIL: common Y_3 coordinate does not have exactly two terms");
fi;
if Set(List(CommonCoordinate, term -> term[2])) <> [-1] then
    Gamma4Fail("FAIL: common-coordinate integer coefficients are not both -1");
fi;

Sg := CharacterAtom("S", G);
Sug := CharacterAtom("S", UGElt);
Sut := CharacterAtom("S", UTElt);
SugKey := Sug[1][1];

## The numerator is intrinsically distinguished by its sigma(ug)-exponent 3;
## the other term has exponent 2.  This is the C_+/C_- orientation used in
## the displayed Xi_3 formula.
if ExponentOf(CommonCoordinate[1][1], SugKey) = 3 then
    NumeratorTerm := CommonCoordinate[1][1];
    DenominatorTerm := CommonCoordinate[2][1];
elif ExponentOf(CommonCoordinate[2][1], SugKey) = 3 then
    NumeratorTerm := CommonCoordinate[2][1];
    DenominatorTerm := CommonCoordinate[1][1];
else
    Gamma4Fail("FAIL: cannot identify the two common-coordinate contributions");
fi;
if ExponentOf(DenominatorTerm, SugKey) <> 2 then
    Gamma4Fail("FAIL: denominator contribution has wrong sigma(ug) exponent");
fi;

## sigma(g)sigma(ug)/(Phi_g(g,ug)sigma(ut))=1.  Dividing by this
## projective defect and inserting sigma(g)=-1 converts the signed common
## coordinate ratio into the statement's Xi_3.
ProjectiveDefect := DivideMonomial(
    MultiplyMonomial(Sg, Sug),
    MultiplyMonomial(PhiLocal(G, G, UGElt), Sut)
);
XiFromCoordinate := DivideMonomial(NumeratorTerm, DenominatorTerm);
XiFromCoordinate := MultiplyMonomial(XiFromCoordinate, Sg);
XiFromCoordinate := DivideMonomial(XiFromCoordinate, ProjectiveDefect);

## Eight-atom target from the displayed formula, with x=ug and b=epsilon*g.
XG := MultiplyElt(UGElt, G);
HX := MultiplyElt(H, UGElt);
XiDirect := PhiAtom(G, UGElt, G);
XiDirect := AddSparse(XiDirect, PhiAtom(H, G, UGElt));
XiDirect := AddSparse(XiDirect, PhiAtom(H, XG, G));
XiDirect := AddSparse(XiDirect, PhiAtom(EGElt, H, XG));
XiDirect := AddSparse(XiDirect, ScaleSparse(PhiAtom(H, G, XG), -1));
XiDirect := AddSparse(XiDirect, ScaleSparse(PhiAtom(H, UGElt, G), -1));
XiDirect := AddSparse(XiDirect, ScaleSparse(PhiAtom(EGElt, H, UGElt), -1));
XiDirect := AddSparse(XiDirect, ScaleSparse(PhiAtom(EGElt, HX, G), -1));

SourceToTargetResidual := AddSparse(XiFromCoordinate,
                                    ScaleSparse(XiDirect, -1));
if Length(SourceToTargetResidual) <> 0 then
    Gamma4Fail("FAIL: recursive common-coordinate ratio is not the eight-atom target");
fi;
if Length(XiDirect) <> 8 then
    Gamma4Fail("FAIL: Xi_3 target must have eight Laurent atoms");
fi;

D1 := CocycleDefect(H, G, UGElt, G);
D2 := CocycleDefect(EGElt, H, UGElt, G);
CertificateVector := AddSparse(D1, ScaleSparse(D2, -1));
Residual := AddSparse(XiDirect, ScaleSparse(CertificateVector, -1));
if Length(Residual) <> 0 then
    Gamma4Fail("FAIL: nonzero Laurent residual in the two-defect certificate");
fi;

Print("Gamma4 Xi3 source-to-target cancellation certificate\n");
Print("Arithmetic: exact integers in a formal sparse Laurent lattice\n");
Print("External GAP packages required: none\n");
Print("Corrected recursive operators reconstructed: phi_1, phi_2, phi_3\n");
Print("Y2 nonzero recursive coordinates: ", Length(Y2SupportPositions), "\n");
Print("First Y3 coefficient: beta, B1, nu1, mu1 reconstructed from source coordinates\n");
Print("Unspecialized product beta*B1*nu1*mu1: -sigma(g)\n");
Print("Epsilon-action comparison paths: residuals=0/0\n");
Print("At sigma(g)=-1: beta*B1*nu1*mu1=1\n");
Print("First Y3 recursive coordinate: terms=2, specialized residual=0\n");
Print("Common Y3 coordinate (zero-based): 45\n");
Print("Common-coordinate Laurent terms: ", Length(CommonCoordinate), "\n");
Print("Source-to-eight-atom-target residual: ",
      Length(SourceToTargetResidual), "\n");
Print("Direct Xi3 Laurent atoms: ", Length(XiDirect), "\n");
Print("Certificate: +d[h|g|ug|g] - d[epsilon*g|h|ug|g]\n");
Print("Final residual Laurent atoms: ", Length(Residual), "\n");
Print("PASS: both Y3 coefficients vanish by independent source calculations\n");
