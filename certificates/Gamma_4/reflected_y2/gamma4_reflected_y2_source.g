#############################################################################
## Independent source reconstruction for the R1-reflected Gamma4 Y2 rows.
##
## This file is read by gamma4_reflected_y2_certificate.g after the exact
## Laurent helper functions and the delivered data have been loaded.  It
## reconstructs, over abstract Gamma4 normal forms, the induced actions on
## the initial pair (V,W), the scalar formulas for the first-root character
## on X1, the rigid-dual action on V*, and the induced actions on the
## reflected pair (V*,X1).
## It also reconstructs the four z0 coordinates and their factorization at
## sigma(g)=-1.  It then constructs phi_1 and phi_2 with the corrected
## right-associated c_12 convention and returns the seven line-stability ratios
## Theta_2_j, 1 <= j <= 7.  The eight nonzero coordinates are numbered
## 0,...,7 so that conjugation by h acts as (0 1 2 3)(4 5 6 7).
##
## The program reads the scalar for the action on X1 from one nonzero
## coordinate of the displayed first-root vector.  It does not prove that
## the second coordinate gives the same scalar.  Stability and simplicity
## of the one-dimensional homogeneous component are supplied by
## Lemma lem:Gamma4-first-adjoint in the manuscript.
#############################################################################

SRC_One := [0, 0, 0];
SRC_Epsilon := [1, 0, 0];
SRC_U := [2, 0, 0];
SRC_H := [0, 1, 0];
SRC_G := [0, 0, 1];
SRC_T := [0, 0, 2];
SRC_Z := [3, 2, 0];
SRC_D := Multiply(SRC_H, SRC_G);

SRC_KeyLess := function(left, right)
    return String(left) < String(right);
end;
SRC_AddExponent := function(monomial, key, exponent)
    local i;
    if exponent = 0 then
        return;
    fi;
    for i in [1 .. Length(monomial)] do
        if monomial[i][1] = key then
            monomial[i][2] := monomial[i][2] + exponent;
            if monomial[i][2] = 0 then
                Remove(monomial, i);
            fi;
            return;
        fi;
    od;
    Add(monomial, [StructuralCopy(key), exponent]);
end;

SRC_CanonicalMonomial := function(monomial)
    local answer, term;
    answer := [];
    for term in monomial do
        SRC_AddExponent(answer, term[1], term[2]);
    od;
    Sort(answer, function(left, right)
        return SRC_KeyLess(left[1], right[1]);
    end);
    return answer;
end;

SRC_MultiplyMonomials := function(left, right)
    local answer, term;
    answer := List(left, term -> [StructuralCopy(term[1]), term[2]]);
    for term in right do
        SRC_AddExponent(answer, term[1], term[2]);
    od;
    return SRC_CanonicalMonomial(answer);
end;

SRC_ScaleMonomial := function(monomial, exponent)
    return SRC_CanonicalMonomial(List(monomial,
        term -> [StructuralCopy(term[1]), exponent * term[2]]));
end;

SRC_InverseMonomial := function(monomial)
    return SRC_ScaleMonomial(monomial, -1);
end;

SRC_Phi := function(a, b, c)
    if a = SRC_One or b = SRC_One or c = SRC_One then
        return [];
    fi;
    return [[ ["P", StructuralCopy(a), StructuralCopy(b), StructuralCopy(c)], 1 ]];
end;

SRC_Character := function(name, element)
    if element = SRC_One then
        return [];
    fi;
    return [[ [name, StructuralCopy(element)], 1 ]];
end;

SRC_NegativeOne := [[ ["NEG"], 1 ]];

## Phi_x(a,b).
SRC_PhiSubscript := function(x, a, b)
    local answer;
    answer := SRC_Phi(a, b, x);
    answer := SRC_MultiplyMonomials(answer,
        SRC_Phi(Conjugate(Multiply(a, b), x), a, b));
    answer := SRC_MultiplyMonomials(answer,
        SRC_InverseMonomial(SRC_Phi(a, Conjugate(b, x), b)));
    return answer;
end;

## Phi^a(x,y), the tensor-action correction.
SRC_PhiSuperscript := function(a, x, y)
    local answer;
    answer := SRC_Phi(a, x, y);
    answer := SRC_MultiplyMonomials(answer,
        SRC_Phi(Conjugate(a, x), Conjugate(a, y), a));
    answer := SRC_MultiplyMonomials(answer,
        SRC_InverseMonomial(SRC_Phi(Conjugate(a, x), a, y)));
    return answer;
end;

SRC_CocycleDefect := function(a, b, c, d)
    local answer;
    answer := SRC_Phi(b, c, d);
    answer := SRC_MultiplyMonomials(answer,
        SRC_Phi(a, Multiply(b, c), d));
    answer := SRC_MultiplyMonomials(answer, SRC_Phi(a, b, c));
    answer := SRC_MultiplyMonomials(answer,
        SRC_InverseMonomial(SRC_Phi(Multiply(a, b), c, d)));
    answer := SRC_MultiplyMonomials(answer,
        SRC_InverseMonomial(SRC_Phi(a, b, Multiply(c, d))));
    return answer;
end;

#############################################################################
## Sparse Laurent polynomials and matrices.
## A polynomial is [ [monomial, integer coefficient], ... ].  Polynomial
## insertion order is retained, as in the source recursion; monomials
## themselves are canonical exponent rows.
#############################################################################

SRC_PolynomialAdd := function(left, right, scale)
    local answer, term, i, found;
    answer := List(left,
        term -> [List(term[1], item -> [StructuralCopy(item[1]), item[2]]),
                 term[2]]);
    for term in right do
        found := false;
        for i in [1 .. Length(answer)] do
            if answer[i][1] = term[1] then
                answer[i][2] := answer[i][2] + scale * term[2];
                if answer[i][2] = 0 then
                    Remove(answer, i);
                fi;
                found := true;
                break;
            fi;
        od;
        if not found and scale * term[2] <> 0 then
            Add(answer,
                [List(term[1], item -> [StructuralCopy(item[1]), item[2]]),
                 scale * term[2]]);
        fi;
    od;
    return answer;
end;

SRC_PolynomialMultiply := function(left, right)
    local answer, a, b;
    answer := [];
    for a in left do
        for b in right do
            answer := SRC_PolynomialAdd(answer,
                [[SRC_MultiplyMonomials(a[1], b[1]), a[2] * b[2]]], 1);
        od;
    od;
    return answer;
end;

SRC_MonomialPolynomial := function(monomial, coefficient)
    if coefficient = 0 then
        return [];
    fi;
    return [[SRC_CanonicalMonomial(monomial), coefficient]];
end;

SRC_VectorAdd := function(left, right, scale)
    return List([1 .. Length(left)],
        i -> SRC_PolynomialAdd(left[i], right[i], scale));
end;

SRC_ZeroMatrix := function(rows, columns)
    return List([1 .. rows], i -> List([1 .. columns], j -> []));
end;

SRC_IdentityMatrix := function(dimension)
    local matrix, i;
    matrix := SRC_ZeroMatrix(dimension, dimension);
    for i in [1 .. dimension] do
        matrix[i][i] := SRC_MonomialPolynomial([], 1);
    od;
    return matrix;
end;

SRC_DiagonalMatrix := function(monomials)
    local matrix, i;
    matrix := SRC_ZeroMatrix(Length(monomials), Length(monomials));
    for i in [1 .. Length(monomials)] do
        matrix[i][i] := SRC_MonomialPolynomial(monomials[i], 1);
    od;
    return matrix;
end;

SRC_MatrixMultiply := function(left, right)
    local rows, middle, columns, answer, i, k, j, product;
    rows := Length(left);
    middle := Length(right);
    columns := Length(right[1]);
    answer := SRC_ZeroMatrix(rows, columns);
    for i in [1 .. rows] do
        for k in [1 .. middle] do
            if Length(left[i][k]) > 0 then
                for j in [1 .. columns] do
                    if Length(right[k][j]) > 0 then
                        product := SRC_PolynomialMultiply(left[i][k], right[k][j]);
                        answer[i][j] := SRC_PolynomialAdd(answer[i][j], product, 1);
                    fi;
                od;
            fi;
        od;
    od;
    return answer;
end;

SRC_MatrixApply := function(matrix, vector)
    local answer, i, j, product;
    answer := List([1 .. Length(matrix)], i -> []);
    for i in [1 .. Length(matrix)] do
        for j in [1 .. Length(vector)] do
            if Length(matrix[i][j]) > 0 and Length(vector[j]) > 0 then
                product := SRC_PolynomialMultiply(matrix[i][j], vector[j]);
                answer[i] := SRC_PolynomialAdd(answer[i], product, 1);
            fi;
        od;
    od;
    return answer;
end;

SRC_Kronecker := function(left, right)
    local answer, i, j, k, ell, row, column;
    answer := SRC_ZeroMatrix(Length(left) * Length(right),
                             Length(left[1]) * Length(right[1]));
    for i in [1 .. Length(left)] do
        for j in [1 .. Length(left[1])] do
            if Length(left[i][j]) > 0 then
                for k in [1 .. Length(right)] do
                    for ell in [1 .. Length(right[1])] do
                        if Length(right[k][ell]) > 0 then
                            row := (i - 1) * Length(right) + k;
                            column := (j - 1) * Length(right[1]) + ell;
                            answer[row][column] :=
                                SRC_PolynomialMultiply(left[i][j], right[k][ell]);
                        fi;
                    od;
                od;
            fi;
        od;
    od;
    return answer;
end;

SRC_BasisVector := function(dimension, index)
    local vector;
    vector := List([1 .. dimension], i -> []);
    vector[index] := SRC_MonomialPolynomial([], 1);
    return vector;
end;

#############################################################################
## Induced homogeneous objects, tensor products, braiding, and phi_1.
#############################################################################

SRC_MakeObject := function(degrees, builder)
    return rec(
        degrees := List(degrees, StructuralCopy),
        dim := Length(degrees),
        builder := builder,
        cacheActors := [],
        cacheMatrices := []
    );
end;

SRC_ActionMatrix := function(object, actor)
    local position, matrix;
    position := Position(object.cacheActors, actor);
    if position <> fail then
        return object.cacheMatrices[position];
    fi;
    matrix := object.builder(actor);
    Add(object.cacheActors, StructuralCopy(actor));
    Add(object.cacheMatrices, matrix);
    return matrix;
end;

SRC_Representative := function(representatives, degree)
    local matches;
    matches := Filtered(representatives, item -> item[1] = degree);
    if Length(matches) <> 1 then
        Gamma4Fail("missing or repeated transversal representative for degree ", degree);
    fi;
    return matches[1][2];
end;

SRC_InducedObject := function(x, characterBuilder, representatives)
    local degrees, builder;
    degrees := List(representatives, item -> item[1]);
    builder := function(actor)
        local matrix, column, degree, newDegree, row, repDegree, repNew,
              centralizerElement, gamma;
        matrix := SRC_ZeroMatrix(Length(degrees), Length(degrees));
        for column in [1 .. Length(degrees)] do
            degree := degrees[column];
            newDegree := Conjugate(actor, degree);
            row := Position(degrees, newDegree);
            if row = fail then
                Gamma4Fail("induced support is not stable under actor ", actor);
            fi;
            repDegree := SRC_Representative(representatives, degree);
            repNew := SRC_Representative(representatives, newDegree);
            centralizerElement := Multiply(InverseElt(repNew),
                                            Multiply(actor, repDegree));
            gamma := SRC_Phi(actor, repDegree, x);
            gamma := SRC_MultiplyMonomials(gamma,
                SRC_Phi(newDegree, actor, repDegree));
            gamma := SRC_MultiplyMonomials(gamma,
                SRC_Phi(repNew, x, centralizerElement));
            gamma := SRC_MultiplyMonomials(gamma,
                SRC_InverseMonomial(SRC_Phi(actor, degree, repDegree)));
            gamma := SRC_MultiplyMonomials(gamma,
                SRC_InverseMonomial(SRC_Phi(repNew, centralizerElement, x)));
            gamma := SRC_MultiplyMonomials(gamma,
                SRC_InverseMonomial(
                    SRC_Phi(newDegree, repNew, centralizerElement)));
            gamma := SRC_MultiplyMonomials(gamma,
                characterBuilder(centralizerElement));
            matrix[row][column] := SRC_MonomialPolynomial(gamma, 1);
        od;
        return matrix;
    end;
    return SRC_MakeObject(degrees, builder);
end;

SRC_TensorObject := function(left, right)
    local degrees, builder;
    degrees := Concatenation(List(left.degrees,
        x -> List(right.degrees, y -> Multiply(x, y))));
    builder := function(actor)
        local matrix, leftAction, rightAction, i, j, column, correction,
              ii, jj, row, entry;
        matrix := SRC_ZeroMatrix(Length(degrees), Length(degrees));
        leftAction := SRC_ActionMatrix(left, actor);
        rightAction := SRC_ActionMatrix(right, actor);
        for i in [1 .. left.dim] do
            for j in [1 .. right.dim] do
                column := (i - 1) * right.dim + j;
                correction := SRC_PhiSuperscript(actor,
                                                  left.degrees[i],
                                                  right.degrees[j]);
                for ii in [1 .. left.dim] do
                    if Length(leftAction[ii][i]) > 0 then
                        for jj in [1 .. right.dim] do
                            if Length(rightAction[jj][j]) > 0 then
                                row := (ii - 1) * right.dim + jj;
                                entry := SRC_PolynomialMultiply(
                                    SRC_MonomialPolynomial(correction, 1),
                                    SRC_PolynomialMultiply(leftAction[ii][i],
                                                           rightAction[jj][j]));
                                matrix[row][column] := entry;
                            fi;
                        od;
                    fi;
                od;
            od;
        od;
        return matrix;
    end;
    return SRC_MakeObject(degrees, builder);
end;

SRC_Braid := function(left, right)
    local matrix, i, action, j, jj, column, row;
    matrix := SRC_ZeroMatrix(right.dim * left.dim, left.dim * right.dim);
    for i in [1 .. left.dim] do
        action := SRC_ActionMatrix(right, left.degrees[i]);
        for j in [1 .. right.dim] do
            column := (i - 1) * right.dim + j;
            for jj in [1 .. right.dim] do
                if Length(action[jj][j]) > 0 then
                    row := (jj - 1) * left.dim + i;
                    matrix[row][column] := action[jj][j];
                fi;
            od;
        od;
    od;
    return matrix;
end;

SRC_PhiOne := function(left, right)
    local dimension, double, answer, i, j, identityEntry;
    dimension := left.dim * right.dim;
    double := SRC_MatrixMultiply(SRC_Braid(right, left),
                                 SRC_Braid(left, right));
    answer := SRC_ZeroMatrix(dimension, dimension);
    for i in [1 .. dimension] do
        for j in [1 .. dimension] do
            if i = j then
                identityEntry := SRC_MonomialPolynomial([], 1);
            else
                identityEntry := [];
            fi;
            answer[i][j] := SRC_PolynomialAdd(identityEntry, double[i][j], -1);
        od;
    od;
    return answer;
end;

SRC_SignedMonomial := function(polynomial)
    local answer;
    if Length(polynomial) <> 1 then
        Gamma4Fail("expected a single Laurent monomial, got ", Length(polynomial));
    fi;
    if polynomial[1][2] <> 1 and polynomial[1][2] <> -1 then
        Gamma4Fail("expected a unit coefficient");
    fi;
    answer := polynomial[1][1];
    if polynomial[1][2] = -1 then
        answer := SRC_MultiplyMonomials(answer, SRC_NegativeOne);
    fi;
    return answer;
end;

#############################################################################
## Fixed transversals and the first-root/dual characters.
#############################################################################

SRC_OriginalRepresentatives := function()
    local r;
    r := InverseElt(SRC_Epsilon);
    return rec(
        v := [[SRC_H, SRC_One], [Multiply(r, SRC_H), SRC_G]],
        w := [[SRC_G, SRC_One],
              [Multiply(SRC_Epsilon, SRC_G), SRC_H],
              [Multiply(SRC_U, SRC_G), SRC_Epsilon],
              [Multiply(r, SRC_G), Multiply(SRC_Epsilon, SRC_H)]]
    );
end;

SRC_ReflectedRepresentatives := function()
    local ep, hp, gp, rp;
    ep := InverseElt(SRC_Epsilon);
    hp := InverseElt(SRC_H);
    gp := SRC_D;
    rp := SRC_Epsilon;
    return rec(
        v := [[hp, SRC_One], [Multiply(rp, hp), gp]],
        w := [[gp, SRC_One],
              [Multiply(ep, gp), hp],
              [Multiply(SRC_U, gp), ep],
              [Multiply(rp, gp), Multiply(ep, hp)]]
    );
end;

SRC_FirstRootData := fail;
SRC_TauActors := [];
SRC_TauValues := [];

SRC_InitializeFirstRoot := function()
    local reps, v, w, tensor, phiOne, ih, ig, x1, support;
    if SRC_FirstRootData <> fail then
        return;
    fi;
    reps := SRC_OriginalRepresentatives();
    v := SRC_InducedObject(SRC_H,
        c -> SRC_Character("R", c), reps.v);
    w := SRC_InducedObject(SRC_G,
        c -> SRC_Character("S", c), reps.w);
    tensor := SRC_TensorObject(v, w);
    phiOne := SRC_PhiOne(v, w);
    ih := Position(v.degrees, SRC_H);
    ig := Position(w.degrees, SRC_G);
    x1 := SRC_MatrixApply(phiOne,
        SRC_BasisVector(tensor.dim, (ih - 1) * w.dim + ig));
    support := Filtered([1 .. Length(x1)], i -> Length(x1[i]) > 0);
    if Length(support) <> 2 then
        Gamma4Fail("first-root vector must have two nonzero coordinates");
    fi;
    SRC_FirstRootData := rec(tensor := tensor, vector := x1, support := support);
end;

SRC_FirstRootCharacter := function(actor)
    local position, datum, acted, coordinate, value;
    position := Position(SRC_TauActors, actor);
    if position <> fail then
        return SRC_TauValues[position];
    fi;
    SRC_InitializeFirstRoot();
    datum := SRC_FirstRootData;
    acted := SRC_MatrixApply(SRC_ActionMatrix(datum.tensor, actor), datum.vector);
    coordinate := datum.support[1];
    value := SRC_MultiplyMonomials(SRC_SignedMonomial(acted[coordinate]),
        SRC_InverseMonomial(SRC_SignedMonomial(datum.vector[coordinate])));
    Add(SRC_TauActors, StructuralCopy(actor));
    Add(SRC_TauValues, value);
    return value;
end;

SRC_DualCharacter := function(actor)
    local value;
    value := SRC_InverseMonomial(
        SRC_PhiSuperscript(actor, InverseElt(SRC_H), SRC_H));
    return SRC_MultiplyMonomials(value,
        SRC_InverseMonomial(SRC_Character("R", actor)));
end;

SRC_InitialMarking := function()
    local reps;
    reps := SRC_OriginalRepresentatives();
    return rec(
        epsilon := SRC_Epsilon,
        h := SRC_H,
        g := SRC_G,
        v := SRC_InducedObject(SRC_H,
            c -> SRC_Character("R", c), reps.v),
        w := SRC_InducedObject(SRC_G,
            c -> SRC_Character("S", c), reps.w)
    );
end;

SRC_ReflectedMarking := function()
    local ep, hp, gp, reps;
    ep := InverseElt(SRC_Epsilon);
    hp := InverseElt(SRC_H);
    gp := SRC_D;
    reps := SRC_ReflectedRepresentatives();
    return rec(
        epsilon := ep,
        h := hp,
        g := gp,
        v := SRC_InducedObject(hp, SRC_DualCharacter, reps.v),
        w := SRC_InducedObject(gp, SRC_FirstRootCharacter, reps.w)
    );
end;

#############################################################################
## phi_2 and the eight coordinate line-stability ratios.
#############################################################################

SRC_Y2Source := function(marking)
    local epsilon, h, g, u, v, w, d1, f1, wi, vi, y1, d2, double,
          associatorMonomials, associatorIn, associatorOut, c12, w2,
          inputY2, j, term0, term1, term2, y2, hy2, actionMatrix,
          support, nextCoordinate, firstCycle, secondCycle, remaining,
          coordinateOrder, coordinateCoefficients, actionCoefficients,
          permutation, coordinateIndex, degrees, ratios, numerator,
          denominator;
    epsilon := marking.epsilon;
    h := marking.h;
    g := marking.g;
    u := Multiply(epsilon, epsilon);
    v := marking.v;
    w := marking.w;

    d1 := SRC_TensorObject(w, v);
    f1 := SRC_PhiOne(w, v);
    wi := Position(w.degrees, Multiply(epsilon, g));
    vi := Position(v.degrees, h);
    if wi = fail or vi = fail then
        Gamma4Fail("Y1 seed degree is absent from the marked supports");
    fi;
    y1 := SRC_MatrixApply(f1,
        SRC_BasisVector(d1.dim, (wi - 1) * v.dim + vi));

    d2 := SRC_TensorObject(w, d1);
    double := SRC_MatrixMultiply(SRC_Braid(d1, w), SRC_Braid(w, d1));

    ## Corrected convention c12 = a (c_WW tensor id) a^-1:
    ## input factor Phi, output factor Phi^-1.
    associatorMonomials := Concatenation(List(w.degrees,
        x -> Concatenation(List(w.degrees,
            y -> List(v.degrees, z -> SRC_Phi(x, y, z))))));
    associatorIn := SRC_DiagonalMatrix(associatorMonomials);
    associatorOut := SRC_DiagonalMatrix(List(associatorMonomials,
        SRC_InverseMonomial));
    c12 := SRC_MatrixMultiply(associatorOut,
        SRC_MatrixMultiply(
            SRC_Kronecker(SRC_Braid(w, w), SRC_IdentityMatrix(v.dim)),
            associatorIn));

    w2 := Position(w.degrees, Multiply(u, g));
    if w2 = fail then
        Gamma4Fail("Y2 seed degree is absent from the marked support");
    fi;
    inputY2 := List([1 .. d2.dim], i -> []);
    for j in [1 .. Length(y1)] do
        inputY2[(w2 - 1) * d1.dim + j] := y1[j];
    od;
    term0 := inputY2;
    term1 := SRC_MatrixApply(double, inputY2);
    term2 := SRC_MatrixApply(
        SRC_Kronecker(SRC_IdentityMatrix(w.dim), f1),
        SRC_MatrixApply(c12, inputY2));
    y2 := SRC_VectorAdd(SRC_VectorAdd(term0, term1, -1), term2, 1);
    actionMatrix := SRC_ActionMatrix(d2, h);
    hy2 := SRC_MatrixApply(actionMatrix, y2);

    support := Filtered([1 .. Length(y2)], i -> Length(y2[i]) > 0);
    if Length(support) <> 8 then
        Gamma4Fail("Y2 must have exactly eight nonzero coordinates");
    fi;
    if ForAny(support,
              i -> Length(y2[i]) <> 1 or Length(hy2[i]) <> 1) then
        Gamma4Fail("each relevant Y2 and h.Y2 coordinate must be monomial");
    fi;
    degrees := Set(List(support, i -> d2.degrees[i]));
    if Length(degrees) <> 1 then
        Gamma4Fail("the eight Y2 coordinates do not have one homogeneous degree");
    fi;

    ## Recover the two h-orbits directly from the action matrix.  Number the
    ## first orbit 0,1,2,3 and the second orbit 4,5,6,7.  Thus the coordinate
    ## permutation is (0 1 2 3)(4 5 6 7), independently of the ambient
    ## tensor-basis positions used by GAP.
    nextCoordinate := function(coordinate)
        local images;
        images := Filtered([1 .. Length(y2)],
            i -> Length(actionMatrix[i][coordinate]) > 0);
        if Length(images) <> 1 then
            Gamma4Fail("the h-action is not monomial at coordinate ", coordinate);
        fi;
        return images[1];
    end;

    firstCycle := [support[1]];
    for j in [1 .. 3] do
        Add(firstCycle, nextCoordinate(firstCycle[Length(firstCycle)]));
    od;
    if nextCoordinate(firstCycle[4]) <> firstCycle[1] or
       Length(Set(firstCycle)) <> 4 then
        Gamma4Fail("the first coordinate orbit is not a 4-cycle");
    fi;

    remaining := Difference(support, firstCycle);
    if Length(remaining) <> 4 then
        Gamma4Fail("the eight coordinates do not split into two 4-cycles");
    fi;
    secondCycle := [remaining[1]];
    for j in [1 .. 3] do
        Add(secondCycle, nextCoordinate(secondCycle[Length(secondCycle)]));
    od;
    if nextCoordinate(secondCycle[4]) <> secondCycle[1] or
       Set(Concatenation(firstCycle, secondCycle)) <> Set(support) then
        Gamma4Fail("the second coordinate orbit is not a disjoint 4-cycle");
    fi;

    coordinateOrder := Concatenation(firstCycle, secondCycle);
    permutation := [2, 3, 4, 1, 6, 7, 8, 5];
    coordinateCoefficients := List(coordinateOrder,
        coordinate -> SRC_SignedMonomial(y2[coordinate]));
    actionCoefficients := List([1 .. 8], coordinate ->
        SRC_SignedMonomial(
            actionMatrix[
                coordinateOrder[permutation[coordinate]]
            ][
                coordinateOrder[coordinate]
            ]
        ));

    ## If calA_j is the coefficient of coordinate j and
    ## h.E_j=H_j E_pi(j), compute exactly
    ##
    ##   Theta_2_j=H_j calA_j calA_1
    ##             /(H_0 calA_0 calA_pi(j)),  1<=j<=7.
    ratios := [];
    for j in [1 .. 7] do
        coordinateIndex := j + 1;
        numerator := SRC_MultiplyMonomials(
            actionCoefficients[coordinateIndex],
            coordinateCoefficients[coordinateIndex]);
        numerator := SRC_MultiplyMonomials(
            numerator, coordinateCoefficients[2]);
        denominator := SRC_MultiplyMonomials(
            actionCoefficients[1], coordinateCoefficients[1]);
        denominator := SRC_MultiplyMonomials(
            denominator,
            coordinateCoefficients[permutation[coordinateIndex]]);
        Add(ratios, SRC_MultiplyMonomials(numerator,
                                         SRC_InverseMonomial(denominator)));
    od;
    return rec(
        y1 := y1,
        y2 := y2,
        actedY2 := hy2,
        coordinateLabels := [0 .. 7],
        degree := degrees[1],
        ratios := ratios
    );
end;

#############################################################################
## Independent z_0 reconstruction and factorization check.
##
## This starts again from the initial induced pair.  It constructs y_1 by
## phi_1, applies phi_2 to w_g tensor y_1, and only afterwards identifies
## the four multihomogeneous coordinates. 
## The four coefficients used in Lemma
## lem:Gamma4-z0-factorization are outputs of the recursion,
## not source data.
#############################################################################

SRC_SpecializeCharacterMinusOneMonomial := function(monomial, name, element)
    local answer, sign, term;
    answer := [];
    sign := 1;
    for term in monomial do
        if term[1][1] = name and term[1][2] = element then
            if (term[2] mod 2) <> 0 then
                sign := -sign;
            fi;
        else
            Add(answer, [StructuralCopy(term[1]), term[2]]);
        fi;
    od;
    return [SRC_CanonicalMonomial(answer), sign];
end;

SRC_SpecializeCharacterMinusOnePolynomial := function(polynomial, name, element)
    local answer, term, specialized;
    answer := [];
    for term in polynomial do
        specialized := SRC_SpecializeCharacterMinusOneMonomial(
            term[1], name, element);
        answer := SRC_PolynomialAdd(answer,
            [[specialized[1], specialized[2] * term[2]]], 1);
    od;
    return answer;
end;

SRC_SpecializeCharacterMinusOneVector := function(vector, name, element)
    return List(vector,
        polynomial -> SRC_SpecializeCharacterMinusOnePolynomial(
            polynomial, name, element));
end;

SRC_ExpandCharacterPolynomial := function(polynomial)
    local answer, term, expandCharacters;
    answer := [];
    expandCharacters := ValueGlobal("SRC_ExpandCharacters");
    for term in polynomial do
        answer := SRC_PolynomialAdd(answer,
            [[expandCharacters(term[1]), term[2]]], 1);
    od;
    return answer;
end;

SRC_ExpandCharacterVector := function(vector)
    return List(vector, SRC_ExpandCharacterPolynomial);
end;

SRC_Z0Source := function()
    local marking, epsilon, h, g, v, w, d1, f1, wi, vi, y1,
          d2, double, associatorMonomials, associatorIn, associatorOut,
          c12, inputZ0, j, term0, term1, term2, z0, support,
          vInverseH, positionFirst, positionSecond, actedW1,
          positionThird, positionFourth, expectedSupport, betaPosition,
          beta, actionVInverseH, delta, factor, gh, actedW, actedWRow,
          actedWScalar, normalized, expectedCoefficients, sigmaG,
          c12Coefficient, phiGGH, phiRHG, upperOne, upperTwo,
          identityOne, identityTwo, specializedCoefficients,
          specializedExpected, factorizedCoefficients, epsilonH,
          inverseEpsilonG, sourceFourthTerm, expectedFourthTerm,
          coordinateResidual, coordinateCertificate;

    marking := SRC_InitialMarking();
    epsilon := marking.epsilon;
    h := marking.h;
    g := marking.g;
    v := marking.v;
    w := marking.w;

    d1 := SRC_TensorObject(w, v);
    f1 := SRC_PhiOne(w, v);
    wi := Position(w.degrees, Multiply(epsilon, g));
    vi := Position(v.degrees, h);
    y1 := SRC_MatrixApply(f1,
        SRC_BasisVector(d1.dim, (wi - 1) * v.dim + vi));

    d2 := SRC_TensorObject(w, d1);
    double := SRC_MatrixMultiply(SRC_Braid(d1, w), SRC_Braid(w, d1));
    associatorMonomials := Concatenation(List(w.degrees,
        x -> Concatenation(List(w.degrees,
            y -> List(v.degrees, z -> SRC_Phi(x, y, z))))));
    associatorIn := SRC_DiagonalMatrix(associatorMonomials);
    associatorOut := SRC_DiagonalMatrix(List(associatorMonomials,
        SRC_InverseMonomial));
    c12 := SRC_MatrixMultiply(associatorOut,
        SRC_MatrixMultiply(
            SRC_Kronecker(SRC_Braid(w, w), SRC_IdentityMatrix(v.dim)),
            associatorIn));

    inputZ0 := List([1 .. d2.dim], i -> []);
    for j in [1 .. d1.dim] do
        inputZ0[(Position(w.degrees, g) - 1) * d1.dim + j] := y1[j];
    od;
    term0 := inputZ0;
    term1 := SRC_MatrixApply(double, inputZ0);
    term2 := SRC_MatrixApply(
        SRC_Kronecker(SRC_IdentityMatrix(w.dim), f1),
        SRC_MatrixApply(c12, inputZ0));
    z0 := SRC_VectorAdd(SRC_VectorAdd(term0, term1, -1), term2, 1);

    vInverseH := Multiply(InverseElt(epsilon), h);
    positionFirst := (Position(w.degrees, g) - 1) * d1.dim
                   + (Position(w.degrees, Multiply(epsilon, g)) - 1) * v.dim
                   + Position(v.degrees, h);
    positionSecond := (Position(w.degrees, g) - 1) * d1.dim
                    + (Position(w.degrees, g) - 1) * v.dim
                    + Position(v.degrees, vInverseH);

    actedW1 := SRC_ActionMatrix(w, g);
    positionThird := Filtered([1 .. w.dim],
        i -> Length(actedW1[i][Position(w.degrees,
                                      Multiply(epsilon, g))]) > 0);
    if Length(positionThird) <> 1 then
        Gamma4Fail("g acting on w_1 is not monomial");
    fi;
    positionThird := (positionThird[1] - 1) * d1.dim
                   + (Position(w.degrees, g) - 1) * v.dim
                   + Position(v.degrees, h);
    positionFourth := (QuoInt(positionThird - 1, d1.dim)) * d1.dim
                    + (QuoInt(positionThird - 1, d1.dim)) * v.dim
                    + Position(v.degrees, vInverseH);

    expectedSupport := [positionFirst, positionSecond,
                        positionThird, positionFourth];
    support := Filtered([1 .. Length(z0)], i -> Length(z0[i]) > 0);
    if Set(support) <> Set(expectedSupport) then
        Gamma4Fail("z_0 recursion has the wrong four multihomogeneous coordinates");
    fi;
    if ForAny(support, i -> Length(z0[i]) <> 2) then
        Gamma4Fail("a z_0 multihomogeneous coefficient does not have two terms");
    fi;

    ## Extract beta from y_1=w_1 tensor v-beta w tensor v_1, rather than
    ## inserting the displayed formula for beta.
    betaPosition := (Position(w.degrees, g) - 1) * v.dim
                  + Position(v.degrees, vInverseH);
    if Length(y1[betaPosition]) <> 1 or y1[betaPosition][1][2] <> -1 then
        Gamma4Fail("cannot recover beta from the recursive y_1 vector");
    fi;
    beta := y1[betaPosition][1][1];

    ## Recover Phi_h(g,g)rho(g^2) as the actual scalar of g acting on v_1.
    actionVInverseH := SRC_ActionMatrix(v, g);
    if Length(actionVInverseH[Position(v.degrees, h)]
                             [Position(v.degrees, vInverseH)]) <> 1 then
        Gamma4Fail("g acting on v_1 is not monomial");
    fi;
    actionVInverseH := SRC_SignedMonomial(
        actionVInverseH[Position(v.degrees, h)]
                       [Position(v.degrees, vInverseH)]);
    delta := SRC_MultiplyMonomials(beta, actionVInverseH);
    factor := SRC_PolynomialAdd(
        SRC_MonomialPolynomial([], 1),
        SRC_MonomialPolynomial(delta, 1), -1);

    ## Divide the last two standard-basis coordinates by the actual scalar
    ## of (gh) acting on w. 
## The four tensor lines used for the z_0 factorization, rather than
## the transversal basis used internally by the induced module.
    gh := Multiply(g, h);
    actedW := SRC_ActionMatrix(w, gh);
    actedWRow := QuoInt(positionThird - 1, d1.dim) + 1;
    if Length(actedW[actedWRow][Position(w.degrees, g)]) <> 1 then
        Gamma4Fail("(gh) acting on w is not monomial");
    fi;
    actedWScalar := SRC_SignedMonomial(
        actedW[actedWRow][Position(w.degrees, g)]);
    normalized := [
        z0[positionFirst],
        z0[positionSecond],
        SRC_PolynomialMultiply(z0[positionThird],
            SRC_MonomialPolynomial(SRC_InverseMonomial(actedWScalar), 1)),
        SRC_PolynomialMultiply(z0[positionFourth],
            SRC_MonomialPolynomial(
                SRC_ScaleMonomial(actedWScalar, -2), 1))
    ];

    ## Reconstruct the compact four coefficients from the independently
    ## recovered beta, Delta4, tensor-action scalar, and c_12 quotient.
    sigmaG := SRC_Character("S", g);
    c12Coefficient := SRC_MultiplyMonomials(
        SRC_Phi(g, Multiply(epsilon, g), h),
        SRC_InverseMonomial(
            SRC_Phi(Conjugate(g, Multiply(epsilon, g)), g, h)));
    phiGGH := SRC_PhiSubscript(g, g, h);
    phiRHG := SRC_PhiSubscript(g, vInverseH, g);
    upperOne := SRC_PhiSuperscript(g, Multiply(epsilon, g), h);
    upperTwo := SRC_PhiSuperscript(g, g, vInverseH);

    expectedCoefficients := [
        SRC_PolynomialAdd(
            SRC_MonomialPolynomial([], 1),
            SRC_MonomialPolynomial(
                SRC_MultiplyMonomials(sigmaG, delta), 1), 1),
        SRC_PolynomialAdd(
            SRC_MonomialPolynomial(beta, -1),
            SRC_MonomialPolynomial(
                SRC_MultiplyMonomials(beta, sigmaG), -1), 1),
        SRC_PolynomialAdd(
            SRC_MonomialPolynomial(
                SRC_MultiplyMonomials(
                    SRC_MultiplyMonomials(sigmaG, delta), upperTwo), 1),
            SRC_MonomialPolynomial(
                SRC_MultiplyMonomials(c12Coefficient, phiGGH), 1), 1),
        SRC_PolynomialAdd(
            SRC_MonomialPolynomial(
                SRC_MultiplyMonomials(phiGGH, upperOne), -1),
            SRC_MonomialPolynomial(
                SRC_MultiplyMonomials(
                    SRC_MultiplyMonomials(
                        SRC_MultiplyMonomials(phiGGH, c12Coefficient),
                        phiRHG),
                    SRC_InverseMonomial(sigmaG)), -1), 1)
    ];

    normalized := List(normalized, SRC_ExpandCharacterPolynomial);
    expectedCoefficients := List(expectedCoefficients,
        SRC_ExpandCharacterPolynomial);
    if ForAny([1 .. 3],
              j -> Length(SRC_PolynomialAdd(
                  normalized[j], expectedCoefficients[j], -1)) <> 0) then
        Gamma4Fail("one of the first three recursive z_0 coefficients is incorrect");
    fi;
    if Length(normalized[4]) <> 2 or Length(expectedCoefficients[4]) <> 2 or
       normalized[4][1] <> expectedCoefficients[4][1] or
       normalized[4][2][2] <> expectedCoefficients[4][2][2] then
        Gamma4Fail("the fourth recursive z_0 coefficient has the wrong two terms");
    fi;

    ## The second monomial in the fourth coefficient is expressed in the
    ## transversal basis by the source recursion.  Its quotient by the
    ## compact projective-action expression is the following difference of
    ## two normalized 3-cocycle defects.  This bridges the source basis to
  ## the fourth coefficient used for the z_0 factorization without
## assuming that coefficient.
    sourceFourthTerm := normalized[4][2][1];
    expectedFourthTerm := expectedCoefficients[4][2][1];
    coordinateResidual := SRC_MultiplyMonomials(sourceFourthTerm,
        SRC_InverseMonomial(expectedFourthTerm));
    epsilonH := Multiply(epsilon, h);
    inverseEpsilonG := Multiply(InverseElt(epsilon), g);
    coordinateCertificate := SRC_MultiplyMonomials(
        SRC_CocycleDefect(epsilonH, g, Multiply(epsilon, epsilon), g),
        SRC_InverseMonomial(
            SRC_CocycleDefect(inverseEpsilonG, epsilonH,
                              Multiply(epsilon, epsilon), g)));
    if coordinateResidual <> coordinateCertificate then
        Gamma4Fail("the fourth z_0 coordinate has a nonzero source-to-target residual");
    fi;
    normalized := expectedCoefficients;

    ## These two identities follow directly by expanding Phi^g and Phi_g.
    ## They are checked before sigma(g) is specialized.
    identityOne := SRC_MultiplyMonomials(c12Coefficient, phiRHG);
    identityTwo := SRC_MultiplyMonomials(c12Coefficient, phiGGH);
    if upperOne <> identityOne or upperTwo <> identityTwo then
        Gamma4Fail("a cocycle quotient used in the z_0 factorization is incorrect");
    fi;

    specializedCoefficients := List(expectedCoefficients,
        polynomial -> SRC_SpecializeCharacterMinusOnePolynomial(
            polynomial, "S", g));
    specializedExpected := [factor, [],
        SRC_PolynomialMultiply(factor,
            SRC_MonomialPolynomial(identityTwo, 1)), []];
    specializedExpected := List(specializedExpected,
        SRC_ExpandCharacterPolynomial);
    factorizedCoefficients := List(specializedExpected,
        polynomial -> SRC_SpecializeCharacterMinusOnePolynomial(
            polynomial, "S", g));
    if specializedCoefficients[4] <> [] then
        Gamma4Fail("the last z_0 coefficient does not vanish when sigma(g)=-1");
    fi;
    if specializedCoefficients[2] <> [] then
        Gamma4Fail("the w tensor (w tensor v_1) coefficient does not vanish when sigma(g)=-1");
    fi;
    if ForAny([1 .. 4],
              j -> Length(SRC_PolynomialAdd(
                  specializedCoefficients[j],
                  factorizedCoefficients[j], -1)) <> 0) then
        Gamma4Fail("the specialized z_0 coefficients do not have the factor 1-Delta4");
    fi;

    return rec(
        support := support,
        termCounts := List(normalized, Length),
        specializedSupport := Filtered([1 .. 4],
            i -> Length(specializedCoefficients[i]) > 0),
        factorized := true
    );
end;

#############################################################################
## Canonical expansion of one-dimensional projective-character values.
#############################################################################

SRC_AppendFactors := function(factors, generator, exponent)
    local sign, i;
    if exponent >= 0 then sign := 1; else sign := -1; fi;
    for i in [1 .. AbsInt(exponent)] do
        Add(factors, [generator, sign]);
    od;
end;

SRC_CharacterWordValue := function(x, name, factors)
    local value, current, datum, generator, sign, factor, generatorValue;
    value := [];
    current := SRC_One;
    for datum in factors do
        generator := datum[1];
        sign := datum[2];
        if sign = 1 then
            factor := generator;
            generatorValue := SRC_Character(name, generator);
        else
            factor := InverseElt(generator);
            generatorValue := SRC_MultiplyMonomials(
                SRC_PhiSubscript(x, generator, factor),
                SRC_InverseMonomial(SRC_Character(name, generator)));
        fi;
        value := SRC_MultiplyMonomials(value, generatorValue);
        value := SRC_MultiplyMonomials(value,
            SRC_InverseMonomial(SRC_PhiSubscript(x, current, factor)));
        current := Multiply(current, factor);
    od;
    return rec(element := current, value := value);
end;

SRC_CanonicalCharacterValue := function(name, element)
    local e, j, k, factors, half, b, a, datum;
    e := element[1]; j := element[2]; k := element[3];
    factors := [];
    if name = "R" then
        if k mod 2 <> 0 then
            Gamma4Fail("R-value requested outside C_G(h): ", element);
        fi;
        SRC_AppendFactors(factors, SRC_Epsilon, e);
        SRC_AppendFactors(factors, SRC_H, j);
        half := QuoInt(k, 2);
        SRC_AppendFactors(factors, SRC_T, half);
        datum := SRC_CharacterWordValue(SRC_H, name, factors);
    elif name = "S" then
        if j mod 2 <> 0 then
            Gamma4Fail("S-value requested outside C_G(g): ", element);
        fi;
        b := QuoInt(j, 2);
        if (e + b) mod 2 <> 0 then
            Gamma4Fail("invalid epsilon exponent in C_G(g): ", element);
        fi;
        a := QuoInt(e + b, 2) mod 2;
        SRC_AppendFactors(factors, SRC_U, a);
        SRC_AppendFactors(factors, SRC_Z, b);
        SRC_AppendFactors(factors, SRC_G, k);
        datum := SRC_CharacterWordValue(SRC_G, name, factors);
    else
        Gamma4Fail("unknown character name ", name);
    fi;
    if datum.element <> element then
        Gamma4Fail("canonical projective-character word has wrong product for ",
              name, " at ", element, ": ", datum.element);
    fi;
    return datum.value;
end;

SRC_ExpandCharacters := function(monomial)
    local answer, term, key, exponent, atom;
    answer := [];
    for term in monomial do
        key := term[1]; exponent := term[2];
        if key[1] = "R" or key[1] = "S" then
            atom := SRC_CanonicalCharacterValue(key[1], key[2]);
        else
            atom := [[StructuralCopy(key), 1]];
        fi;
        answer := SRC_MultiplyMonomials(answer,
                                        SRC_ScaleMonomial(atom, exponent));
    od;
    return answer;
end;

SRC_AtomString := function(key)
    if key[1] = "P" then
        return PhiKey(key[2], key[3], key[4]);
    elif key[1] = "R" or key[1] = "S" then
        return Concatenation(key[1], "(", String(key[2]), ")");
    elif key[1] = "NEG" then
        return "NEG";
    fi;
    Gamma4Fail("unknown source atom ", key);
end;

SRC_StringRow := function(monomial)
    local answer, term;
    answer := [];
    for term in monomial do
        AddExponent(answer, SRC_AtomString(term[1]), term[2]);
    od;
    return CanonicalSparse(answer);
end;

SRC_Delta4Value := function(epsilon, h, g, rhoValue, sigmaValue)
    local inverseEpsilon, answer;
    inverseEpsilon := InverseElt(epsilon);
    answer := SRC_PhiSubscript(h, g, Multiply(epsilon, g));
    answer := SRC_MultiplyMonomials(answer,
        SRC_PhiSubscript(g, Multiply(inverseEpsilon, h), h));
    answer := SRC_MultiplyMonomials(answer,
        rhoValue(Multiply(inverseEpsilon, Multiply(g, g))));
    answer := SRC_MultiplyMonomials(answer,
        sigmaValue(Multiply(inverseEpsilon, Multiply(h, h))));
    return answer;
end;

SRC_ReconstructPacketRows := function()
    local ep, hp, gp, rows, value;
    ep := InverseElt(SRC_Epsilon);
    hp := InverseElt(SRC_H);
    gp := SRC_D;
    rows := [];

    value := SRC_MultiplyMonomials(
        SRC_Character("R", SRC_H), SRC_InverseMonomial(SRC_NegativeOne));
    Add(rows, SRC_StringRow(value));
    value := SRC_MultiplyMonomials(
        SRC_Character("S", SRC_G), SRC_InverseMonomial(SRC_NegativeOne));
    Add(rows, SRC_StringRow(value));

    value := SRC_Delta4Value(SRC_Epsilon, SRC_H, SRC_G,
        c -> SRC_Character("R", c), c -> SRC_Character("S", c));
    Add(rows, SRC_StringRow(SRC_ExpandCharacters(value)));
    value := SRC_Delta4Value(ep, hp, gp,
        SRC_DualCharacter, SRC_FirstRootCharacter);
    Add(rows, SRC_StringRow(SRC_ExpandCharacters(value)));

    value := SRC_MultiplyMonomials(
        SRC_DualCharacter(hp), SRC_InverseMonomial(SRC_NegativeOne));
    Add(rows, SRC_StringRow(SRC_ExpandCharacters(value)));
    value := SRC_MultiplyMonomials(
        SRC_FirstRootCharacter(gp), SRC_InverseMonomial(SRC_NegativeOne));
    Add(rows, SRC_StringRow(SRC_ExpandCharacters(value)));
    return rows;
end;

SRC_ReconstructRows := function()
    local initial, reflected;
    initial := SRC_Y2Source(SRC_InitialMarking());
    reflected := SRC_Y2Source(SRC_ReflectedMarking());
    return rec(
        initial := initial,
        reflected := reflected,
        packetRows := SRC_ReconstructPacketRows(),
        initialRows := List(initial.ratios,
            monomial -> SRC_StringRow(SRC_ExpandCharacters(monomial))),
        reflectedRows := List(reflected.ratios,
            monomial -> SRC_StringRow(SRC_ExpandCharacters(monomial)))
    );
end;
