#############################################################################
## Gamma4 z0 factorization and R1-reflected Y2 parameter certificate
## (core GAP only)
##
## The companion source engine first reconstructs all four z0 coordinates
## from phi_1, phi_2, the tensor action, and the corrected c_12.  It checks
## the last-coordinate cancellation at sigma(g)=-1 and the factorization by
## 1-Delta4.  The verifier then checks two exact integral bridges for each
## of the seven
## reflected parameter identities Theta_2_j:
##
##   (A) Theta_2_j + sum_{i=1}^{13} C[j][i] P_i = B_j(Phi)
##       in the free Laurent exponent lattice; and
##   (B) B_j(Phi) is an integral combination of normalized 3-cocycle
##       defects.
##
## Before checking those bridges, the companion source engine independently
## reconstructs the initial induced objects V,W, the first-root character on
## X1, the dual object V*, the reflected induced pair, both phi_1 and phi_2,
## and the two sets of eight nonzero Y2 coordinates.  In each set the
## coordinates are numbered 0,...,7, with h acting as
## (0 1 2 3)(4 5 6 7).  The seven scalars Theta_2_j, 1 <= j <= 7, are
## expanded canonically and checked coefficient by coefficient against the
## supplied packet and target rows.
## The recursive phi_2 construction uses the right-associated convention
##
##   c_12 = a (c_{W,W} tensor id) a^{-1},
##   coefficient = Phi(x,y,l)/Phi(x*y*x^-1,x,l).
##
## Group elements are abstract normal forms epsilon^i h^j g^k: only i is
## reduced modulo four; j,k remain arbitrary integers.  No finite h- or
## g-order relation, external package, floating point, or random sample is
## used.  Thus the expanded target rows are certificate data, not trusted
## source input.
#############################################################################

Gamma4Fail := function(arg)
    CallFuncList(PrintTo, Concatenation(["*errout*"], arg));
    PrintTo("*errout*", "\n");
    ErrorNoReturn("Certificate verification failed.");
end;

Read("gamma4_reflected_y2_certificate_data.g");
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
    local answer, term;
    answer := [];
    for term in v do
        AddExponent(answer, term[1], term[2]);
    od;
    Sort(answer, function(left, right)
        return left[1] < right[1];
    end);
    return answer;
end;

AddSparse := function(left, right)
    local answer, term;
    answer := List(left, ShallowCopy);
    for term in right do
        AddExponent(answer, term[1], term[2]);
    od;
    return CanonicalSparse(answer);
end;

ScaleSparse := function(v, scalar)
    return CanonicalSparse(List(v, term -> [term[1], scalar * term[2]]));
end;

IdentityElt := [0, 0, 0];

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
    local answer;
    answer := [];
    answer := AddSparse(answer, PhiAtom(b, c, d));
    answer := AddSparse(answer, PhiAtom(a, Multiply(b, c), d));
    answer := AddSparse(answer, PhiAtom(a, b, c));
    answer := AddSparse(answer, ScaleSparse(PhiAtom(Multiply(a, b), c, d), -1));
    answer := AddSparse(answer, ScaleSparse(PhiAtom(a, b, Multiply(c, d)), -1));
    return answer;
end;

## This records, at exponent-vector level, the associator quotient used in
## the source construction of the T_j rows.  It is included explicitly to
## make the corrected c_12 orientation auditable.
C12Coefficient := function(x, y, ell)
    return AddSparse(
        PhiAtom(x, y, ell),
        ScaleSparse(PhiAtom(Conjugate(x, y), x, ell), -1)
    );
end;

#############################################################################
## Convert the seven base-coordinate ratios in the generated data to the
## seven scalars used in the paper.
##
## Put R_0=1 and R_j=L_0/L_j for 1<=j<=7, where L_j is the coefficient
## by which h acts on the j-th coordinate line.  For
##
##   pi=(0 1 2 3)(4 5 6 7),
##
## the scalar in the paper is
##
##   Theta_2_j = L_pi(j)/L_1 = R_1/R_pi(j).
##
## Thus, in the Laurent exponent lattice,
##
##   (Theta_2_1,...,Theta_2_7)
##   =(R_1-R_2, R_1-R_3, R_1,
##     R_1-R_5, R_1-R_6, R_1-R_7, R_1-R_4).
##
## The same integral change of basis is applied to the initial packet rows,
## target rows, residual rows, and bar certificates.  The last seven columns
## of the Stage-A matrix are conjugated by this change of basis.
#############################################################################

ReflectedY2ThetaTransform := [
    [ 1, -1,  0,  0,  0,  0,  0 ],
    [ 1,  0, -1,  0,  0,  0,  0 ],
    [ 1,  0,  0,  0,  0,  0,  0 ],
    [ 1,  0,  0,  0, -1,  0,  0 ],
    [ 1,  0,  0,  0,  0, -1,  0 ],
    [ 1,  0,  0,  0,  0,  0, -1 ],
    [ 1,  0,  0, -1,  0,  0,  0 ]
];

ReflectedY2ThetaTransformInverse := [
    [  0,  0, 1,  0,  0,  0,  0 ],
    [ -1,  0, 1,  0,  0,  0,  0 ],
    [  0, -1, 1,  0,  0,  0,  0 ],
    [  0,  0, 1,  0,  0,  0, -1 ],
    [  0,  0, 1, -1,  0,  0,  0 ],
    [  0,  0, 1,  0, -1,  0,  0 ],
    [  0,  0, 1,  0,  0, -1,  0 ]
];

IntegralLinearCombination := function(coefficients, rows)
    local answer, i;
    answer := [];
    for i in [1 .. Length(rows)] do
        answer := AddSparse(answer,
            ScaleSparse(rows[i], coefficients[i]));
    od;
    return answer;
end;

TransformRows := function(matrix, rows)
    return List(matrix,
        coefficients -> IntegralLinearCombination(coefficients, rows));
end;

IntegerMatrixProduct := function(left, right)
    return List(left, row ->
        List([1 .. Length(right[1])], column ->
            Sum([1 .. Length(right)],
                i -> row[i] * right[i][column])));
end;

ReflectedY2PacketNames := Concatenation(
    ReflectedY2BasePacketNames{[1 .. 6]},
    List([1 .. 7], j -> Concatenation("initial_Theta_2_", String(j)))
);
ReflectedY2PacketRows := Concatenation(
    ReflectedY2BasePacketRows{[1 .. 6]},
    TransformRows(ReflectedY2ThetaTransform,
                  ReflectedY2BasePacketRows{[7 .. 13]})
);
ReflectedY2TargetRows := TransformRows(
    ReflectedY2ThetaTransform, ReflectedY2BaseTargetRows);
ReflectedY2ResidualRows := TransformRows(
    ReflectedY2ThetaTransform, ReflectedY2BaseResidualRows);
ReflectedY2Certificates := List([1 .. 7], j -> rec(
    name := Concatenation("Theta_2_", String(j)),
    certificate := IntegralLinearCombination(
        ReflectedY2ThetaTransform[j],
        List(ReflectedY2BaseCertificates, entry -> entry.certificate)
    )
));

ReflectedY2StageACoefficients := List(
    [1 .. 7], i -> Concatenation(
        IntegerMatrixProduct(
            ReflectedY2ThetaTransform,
            List(ReflectedY2BaseStageACoefficients, row -> row{[1 .. 6]}))
        [i],
        IntegerMatrixProduct(
            IntegerMatrixProduct(
                ReflectedY2ThetaTransform,
                List(ReflectedY2BaseStageACoefficients,
                     row -> row{[7 .. 13]})),
            ReflectedY2ThetaTransformInverse
        )[i]
    )
);

ReflectedY2ExpectedStageACoefficients := [
    [  0, 3, 0, 0, 0, -3, -1, 0,  0, 0, 0, 1, 0 ],
    [  0, 0, 0, 0, 0,  0,  0, 0,  0, 0, 0, 0, 0 ],
    [  0, 1, 0, 0, 0, -1,  0, 0, -1, 0, 0, 0, 0 ],
    [  0, 0, 0, 0, 0,  0,  0, 0,  0, 0, 0, 0, 0 ],
    [  2, 1, 0, 0, 0, -1, -1, 0, -1, 0, 0, 1, 1 ],
    [  0, 0, 0, 0, 0,  0,  0, 0,  0, 0, 0, 0, 0 ],
    [  0, 1, 0, 1, 0, -1,  0, 0, -1, 0, 0, 1, 0 ]
];
if ReflectedY2StageACoefficients <>
   ReflectedY2ExpectedStageACoefficients then
    Gamma4Fail("FAIL: incorrect transformed Stage-A coefficient matrix");
fi;
if ForAny(ReflectedY2StageACoefficients,
          row -> row[3] <> 0 or row[5] <> 0) then
    Gamma4Fail("FAIL: the retained Delta4 or rho1(h1) Stage-A column is nonzero");
fi;

## Reconstruct the initial and reflected induced actions, phi_1, phi_2,
## the eight nonzero Y2 coordinates, and Theta_2_j for 1 <= j <= 7.
## This source engine does not read any of the expanded packet/target rows.
Read("gamma4_reflected_y2_source.g");

VerifyZ0SourceReconstruction := function()
    local source;
    source := SRC_Z0Source();
    if Length(source.support) <> 4 or source.termCounts <> [2, 2, 2, 2] then
        Gamma4Fail("FAIL: incomplete z_0 four-coordinate reconstruction");
    fi;
    if Length(source.specializedSupport) <> 2 or not source.factorized then
        Gamma4Fail("FAIL: incomplete z_0 factorization at sigma(g)=-1");
    fi;
    Print("z0 source recursion: multihomogeneous coordinates=4, ",
          "terms per coordinate=", source.termCounts, "\n");
    Print("z0 at sigma(g)=-1: surviving coordinates=2, ",
          "last coefficient=0, factor (1-Delta4) verified\n");
end;

VerifySourceReconstruction := function()
    local source, index, expectedPacket, expectedInitial, expectedReflected;
    source := SRC_ReconstructRows();
    if source.initial.degree <> [1, 1, 2] then
        Gamma4Fail("FAIL: wrong initial terminal degree: ", source.initial.degree);
    fi;
    if source.reflected.degree <> [2, 1, 2] then
        Gamma4Fail("FAIL: wrong reflected terminal degree: ",
              source.reflected.degree);
    fi;
    if source.initial.coordinateLabels <> [0 .. 7] or
       source.reflected.coordinateLabels <> [0 .. 7] then
        Gamma4Fail("FAIL: the Y2 coordinates are not numbered 0,...,7");
    fi;
    for index in [1 .. 6] do
        expectedPacket := CanonicalSparse(ReflectedY2PacketRows[index]);
        if source.packetRows[index] <> expectedPacket then
            Gamma4Fail("FAIL: reconstructed non-Theta packet row at index ",
                       index);
        fi;
    od;
    for index in [1 .. 7] do
        expectedInitial := CanonicalSparse(ReflectedY2PacketRows[index + 6]);
        expectedReflected := CanonicalSparse(ReflectedY2TargetRows[index]);
        if source.initialRows[index] <> expectedInitial then
            Gamma4Fail("FAIL: reconstructed initial Y2 row at index ", index);
        fi;
        if source.reflectedRows[index] <> expectedReflected then
            Gamma4Fail("FAIL: reconstructed reflected Y2 row at index ", index);
        fi;
    od;
    Print("Source recursion: initial/reflected phi2 coordinates=8/8, ",
          "Theta rows checked=7/7, terminal degrees=",
          [source.initial.degree, source.reflected.degree], "\n");
end;

VerifyStageA := function(index)
    local row, i, residual, nonPhi;
    row := ReflectedY2TargetRows[index];
    for i in [1 .. Length(ReflectedY2PacketRows)] do
        row := AddSparse(
            row,
            ScaleSparse(
                ReflectedY2PacketRows[i],
                ReflectedY2StageACoefficients[index][i]
            )
        );
    od;
    residual := AddSparse(
        row,
        ScaleSparse(ReflectedY2ResidualRows[index], -1)
    );
    if Length(residual) <> 0 then
        Gamma4Fail("FAIL: Stage-A source-to-target identity at index ", index);
    fi;
    nonPhi := Filtered(row, term -> PositionSublist(term[1], "P(") <> 1);
    if Length(nonPhi) <> 0 then
        Gamma4Fail("FAIL: Stage-A row still has character/sign atoms at index ", index);
    fi;
    Print("Theta_2_", ReflectedY2ThetaIndices[index],
          ": source atoms=", Length(ReflectedY2TargetRows[index]),
          ", packet residual atoms=", Length(row),
          ", packet coefficient l1=",
          Sum(ReflectedY2StageACoefficients[index], AbsInt),
          "\n");
end;

VerifyBarCertificate := function(index)
    local entry, target, combination, term, quad, residual,
          maxEntryH, maxEntryG, maxProductH, maxProductG,
          x, position, product;
    entry := ReflectedY2Certificates[index];
    target := ReflectedY2ResidualRows[index];
    combination := [];
    maxEntryH := 0; maxEntryG := 0;
    maxProductH := 0; maxProductG := 0;
    for term in entry.certificate do
        if not IsInt(term[2]) then
            Gamma4Fail("FAIL: nonintegral certificate coefficient in ", entry.name);
        fi;
        quad := term[1];
        combination := AddSparse(
            combination,
            ScaleSparse(
                CocycleDefect(quad[1], quad[2], quad[3], quad[4]),
                term[2]
            )
        );
        for x in quad do
            maxEntryH := Maximum(maxEntryH, AbsInt(x[2]));
            maxEntryG := Maximum(maxEntryG, AbsInt(x[3]));
        od;
        for position in [1 .. 3] do
            product := Multiply(quad[position], quad[position + 1]);
            maxProductH := Maximum(maxProductH, AbsInt(product[2]));
            maxProductG := Maximum(maxProductG, AbsInt(product[3]));
        od;
    od;
    residual := AddSparse(target, ScaleSparse(combination, -1));
    if Length(residual) <> 0 then
        Gamma4Fail("FAIL: nonzero abstract bar residual for ", entry.name);
    fi;
    Print(entry.name,
          ": certificate terms=", Length(entry.certificate),
          ", coefficient l1=", Sum(entry.certificate, term -> AbsInt(term[2])),
          ", residual atoms=", Length(residual),
          ", max |h,g| in entries=", [maxEntryH, maxEntryG],
          ", max |h,g| in adjacent products=", [maxProductH, maxProductG],
          "\n");
end;

if Length(ReflectedY2PacketRows) <> 13 then
    Gamma4Fail("FAIL: expected thirteen packet rows");
fi;
if Length(ReflectedY2TargetRows) <> 7 or
   Length(ReflectedY2ResidualRows) <> 7 or
   Length(ReflectedY2Certificates) <> 7 then
    Gamma4Fail("FAIL: expected seven reflected equations");
fi;

## A nontrivial orientation check for the displayed c_12 quotient.
if Length(C12Coefficient([0,0,1], [0,1,0], [0,0,1])) = 0 then
    Gamma4Fail("FAIL: degenerate c_12 orientation check");
fi;

VerifyZ0SourceReconstruction();
VerifySourceReconstruction();

for index in [1 .. 7] do
    VerifyStageA(index);
    VerifyBarCertificate(index);
od;

Print("Gamma4 R1-reflected Y2 parameter certificate\n");
Print("Arithmetic: exact integers in the Laurent lattice and normalized bar complex\n");
Print("c12 convention: Phi(x,y,l)/Phi(x*y*x^-1,x,l)\n");
Print("External GAP packages required: none\n");
Print("Group model: epsilon exponent modulo 4; h,g exponents integral\n");
Print("No h- or g-order relation used\n");
Print("PASS: z0 factorization and seven reflected Theta_2_j identities verified\n");
