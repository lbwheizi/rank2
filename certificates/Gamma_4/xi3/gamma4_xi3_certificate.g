#############################################################################
## Gamma4 Xi3 source-to-target and direct-cancellation certificate
## (core GAP only)
##
## PURPOSE
##   Reconstruct the common coordinate used in the Y_3 calculation directly
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
## REPRODUCTION
##   gap --bare -A -r -q --nointeract gamma4_xi3_certificate.g
#############################################################################

Gamma4Fail := function(arg)
    CallFuncList(PrintTo, Concatenation(["*errout*"], arg));
    PrintTo("*errout*", "\n");
    QUIT_GAP(1);
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

## y_3 = phi_3(w_g tensor y_2).
Y3Input := ZeroPolyVector(WObj.dim * D2Obj.dim);
for IndexJ in [1 .. D2Obj.dim] do
    Y3Input[(PositionDegree(WObj.degrees, G) - 1) * D2Obj.dim + IndexJ]
        := Y2Vector[IndexJ];
od;
Y3Vector := HigherPhiApply(WObj, D2Obj, D1Obj,
    Phi2WVonVector, Y3Input);

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
Print("Common Y3 coordinate (zero-based): 45\n");
Print("Common-coordinate Laurent terms: ", Length(CommonCoordinate), "\n");
Print("Source-to-eight-atom-target residual: ",
      Length(SourceToTargetResidual), "\n");
Print("Direct Xi3 Laurent atoms: ", Length(XiDirect), "\n");
Print("Certificate: +d[h|g|ug|g] - d[epsilon*g|h|ug|g]\n");
Print("Final residual Laurent atoms: ", Length(Residual), "\n");
Print("PASS: recursive coordinate quotient equals the stated Xi3 and Xi3=1\n");

QUIT_GAP(0);
