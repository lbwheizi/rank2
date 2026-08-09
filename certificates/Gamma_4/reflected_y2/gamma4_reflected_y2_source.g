#############################################################################
## Independent source reconstruction for the R1-reflected Gamma4 Y2 rows.
##
## This file is read by gamma4_reflected_y2_certificate.g after the exact
## Laurent helper functions and the delivered data have been loaded.  It
## reconstructs, over abstract Gamma4 normal forms, the induced actions on
## the initial pair (V,W), the first-root character on X1, the rigid-dual
## action on V*, and the induced actions on the reflected pair (V*,X1).
## It then constructs phi_1 and phi_2 with the corrected right-associated
## c_12 convention and returns the seven line-stability ratios based at the
## first of the eight nonzero Y2 coordinates.
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
        Error("missing or repeated transversal representative for degree ", degree);
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
                Error("induced support is not stable under actor ", actor);
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
        Error("expected a single Laurent monomial, got ", Length(polynomial));
    fi;
    if polynomial[1][2] <> 1 and polynomial[1][2] <> -1 then
        Error("expected a unit coefficient");
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
        Error("first-root vector must have two nonzero coordinates");
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
    local ep, hp, gp, reps, rhoDual;
    ep := InverseElt(SRC_Epsilon);
    hp := InverseElt(SRC_H);
    gp := SRC_D;
    reps := SRC_ReflectedRepresentatives();
    rhoDual := function(c)
        local value;
        value := SRC_InverseMonomial(SRC_PhiSuperscript(c, hp, SRC_H));
        return SRC_MultiplyMonomials(value,
            SRC_InverseMonomial(SRC_Character("R", c)));
    end;
    return rec(
        epsilon := ep,
        h := hp,
        g := gp,
        v := SRC_InducedObject(hp, rhoDual, reps.v),
        w := SRC_InducedObject(gp, SRC_FirstRootCharacter, reps.w)
    );
end;

#############################################################################
## phi_2 and the eight coordinate line-stability ratios.
#############################################################################

SRC_Y2Source := function(marking)
    local epsilon, h, g, u, v, w, d1, f1, wi, vi, y1, d2, double,
          associatorMonomials, associatorIn, associatorOut, c12, w2,
          inputY2, j, term0, term1, term2, y2, hy2, support,
          supportZeroBased, degrees, base, target, ratios, numerator,
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
        Error("Y1 seed degree is absent from the marked supports");
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
        Error("Y2 seed degree is absent from the marked support");
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
    hy2 := SRC_MatrixApply(SRC_ActionMatrix(d2, h), y2);

    support := Filtered([1 .. Length(y2)], i -> Length(y2[i]) > 0);
    supportZeroBased := List(support, i -> i - 1);
    if supportZeroBased <> [5, 6, 8, 15, 17, 18, 27, 28] then
        Error("unexpected Y2 coordinate support: ", supportZeroBased);
    fi;
    if ForAny(support,
              i -> Length(y2[i]) <> 1 or Length(hy2[i]) <> 1) then
        Error("each relevant Y2 and h.Y2 coordinate must be monomial");
    fi;
    degrees := Set(List(support, i -> d2.degrees[i]));
    if Length(degrees) <> 1 then
        Error("the eight Y2 coordinates do not have one homogeneous degree");
    fi;

    base := support[1];
    ratios := [];
    for target in support{[2 .. 8]} do
        numerator := SRC_MultiplyMonomials(
            SRC_SignedMonomial(hy2[base]), SRC_SignedMonomial(y2[target]));
        denominator := SRC_MultiplyMonomials(
            SRC_SignedMonomial(hy2[target]), SRC_SignedMonomial(y2[base]));
        Add(ratios, SRC_MultiplyMonomials(numerator,
                                         SRC_InverseMonomial(denominator)));
    od;
    return rec(
        y1 := y1,
        y2 := y2,
        actedY2 := hy2,
        support := supportZeroBased,
        degree := degrees[1],
        ratios := ratios
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
            Error("R-value requested outside C_G(h): ", element);
        fi;
        SRC_AppendFactors(factors, SRC_Epsilon, e);
        SRC_AppendFactors(factors, SRC_H, j);
        half := QuoInt(k, 2);
        SRC_AppendFactors(factors, SRC_T, half);
        datum := SRC_CharacterWordValue(SRC_H, name, factors);
    elif name = "S" then
        if j mod 2 <> 0 then
            Error("S-value requested outside C_G(g): ", element);
        fi;
        b := QuoInt(j, 2);
        if (e + b) mod 2 <> 0 then
            Error("invalid epsilon exponent in C_G(g): ", element);
        fi;
        a := QuoInt(e + b, 2) mod 2;
        SRC_AppendFactors(factors, SRC_U, a);
        SRC_AppendFactors(factors, SRC_Z, b);
        SRC_AppendFactors(factors, SRC_G, k);
        datum := SRC_CharacterWordValue(SRC_G, name, factors);
    else
        Error("unknown character name ", name);
    fi;
    if datum.element <> element then
        Error("canonical projective-character word has wrong product for ",
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
    Error("unknown source atom ", key);
end;

SRC_StringRow := function(monomial)
    local answer, term;
    answer := [];
    for term in monomial do
        AddExponent(answer, SRC_AtomString(term[1]), term[2]);
    od;
    return CanonicalSparse(answer);
end;

SRC_ReconstructRows := function()
    local initial, reflected;
    initial := SRC_Y2Source(SRC_InitialMarking());
    reflected := SRC_Y2Source(SRC_ReflectedMarking());
    return rec(
        initial := initial,
        reflected := reflected,
        initialRows := List(initial.ratios,
            monomial -> SRC_StringRow(SRC_ExpandCharacters(monomial))),
        reflectedRows := List(reflected.ratios,
            monomial -> SRC_StringRow(SRC_ExpandCharacters(monomial)))
    );
end;
