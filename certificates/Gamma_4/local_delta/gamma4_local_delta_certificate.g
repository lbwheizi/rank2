#############################################################################
## Gamma4 local-Delta Laurent certificate (core GAP only)
##
## PURPOSE
##   Verify the two exact Laurent identities extracted from Y3=0:
##
##       (Delta_01^W)^(-1) = product E_l^(m_l),
##       (Delta_10^W)^(-1) = product E_l^(n_l),
##
##   where every displayed Y3 coordinate equation is E_l=-1 and the sums
##   of the two sets of certificate exponents are zero.  Hence
##   Delta_01^W=Delta_10^W=1.
##
## ASSOCIATOR CONVENTION
##   For right-associated tensors the corrected operator is
##
##       c_12^Phi = a ((c tensor id)) a^(-1),
##
##   so the diagonal factor before c tensor id is Phi and the factor after
##   it is Phi^(-1).  Both the independent recursion and the expected
##   relation snapshot use precisely this convention.
##
## SCOPE
##   The script independently reconstructs V,W, the twisted tensor actions,
##   braidings, phi_1, phi_2, phi_3, y_2 and g acting on y_2 from the Gamma4
##   normal-form group law and the formal coefficients R,S,F,P.  It obtains
##   all 192 nonzero Y3 coordinate binomials, checks the 19 rows used by the
##   certificate in the same quotient orientation, and then verifies the exact
##   Laurent elimination.  It does not numerically sample a cocycle and does
##   not assume finite-dimensionality.  The data file contains only the 19
##   expected rows and two integer certificate vectors; it is not trusted
##   until the independent reconstruction checks it.  No external GAP
##   package is used.
##
## REPRODUCTION (from this directory)
##   gap --bare -A -r -q --nointeract gamma4_local_delta_certificate.g
#############################################################################

Gamma4Fail := function(arg)
    CallFuncList(PrintTo, Concatenation(["*errout*"], arg));
    PrintTo("*errout*", "\n");
    QUIT_GAP(1);
end;

Read("gamma4_local_delta_certificate_data.g");
SizeScreen([1000, 1000]);

if C12AssocIn <> "Phi" or C12AssocOut <> "Phi^-1" then
    Gamma4Fail("FAIL: the relation snapshot has the wrong c12 associator direction");
fi;

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

#############################################################################
## Independent reconstruction of the relevant Y3 coordinate binomials.
##
## Group elements are abstract normal forms epsilon^i h^j g^k.  Only the
## epsilon exponent is reduced modulo four; h- and g-exponents are integral.
#############################################################################

G4One := [0, 0, 0];
G4Epsilon := [1, 0, 0];
G4H := [0, 1, 0];
G4G := [0, 0, 1];

G4Multiply := function(a, b)
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

G4Inverse := function(a)
    local i, j, k;
    i := a[1]; j := a[2]; k := a[3];
    if (k mod 2) = 0 then
        return [(-i) mod 4, -j, -k];
    fi;
    return [(i + j) mod 4, -j, -k];
end;

G4Conjugate := function(a, x)
    return G4Multiply(G4Multiply(a, x), G4Inverse(a));
end;

G4Power := function(a, exponent)
    local ans, i;
    if exponent < 0 then
        return G4Power(G4Inverse(a), -exponent);
    fi;
    ans := G4One;
    for i in [1 .. exponent] do
        ans := G4Multiply(ans, a);
    od;
    return ans;
end;

## Keys deliberately reproduce the readable Python-tuple notation used in
## the stored relation snapshot; the reconstruction itself is pure GAP.
G4EltKey := function(x)
    return Concatenation("(", String(x[1]), ", ", String(x[2]), ", ",
                         String(x[3]), ")");
end;

G4VariableKey := function(kind, arguments)
    return Concatenation("('", kind, "', ",
                         JoinStringsWithSeparator(List(arguments, G4EltKey), ", "),
                         ")");
end;

G4Variable := function(kind, arguments)
    return [[G4VariableKey(kind, arguments), 1]];
end;

G4MonomialMultiply := function(left, right)
    return AddSparse(left, right);
end;

G4MonomialInverse := function(monomial)
    return ScaleSparse(monomial, -1);
end;

G4ActionVariable := function(name, actor, degree)
    return G4Variable(name, [actor, degree]);
end;

G4FVariable := function(actor, x, y)
    return G4Variable("F", [actor, x, y]);
end;

G4PVariable := function(x, y, z)
    return G4Variable("P", [x, y, z]);
end;

#############################################################################
## Sparse Laurent polynomials and matrices.
## A polynomial is [ [monomial, coefficient], ... ].
#############################################################################

G4PolynomialAdd := function(left, right, scale)
    local ans, term, i, found;
    ans := List(left, term -> [List(term[1], ShallowCopy), term[2]]);
    for term in right do
        found := false;
        for i in [1 .. Length(ans)] do
            if ans[i][1] = term[1] then
                ans[i][2] := ans[i][2] + scale * term[2];
                if ans[i][2] = 0 then
                    Remove(ans, i);
                fi;
                found := true;
                break;
            fi;
        od;
        if not found and scale * term[2] <> 0 then
            Add(ans, [List(term[1], ShallowCopy), scale * term[2]]);
        fi;
    od;
    Sort(ans, function(a, b)
        return String(a[1]) < String(b[1]);
    end);
    return ans;
end;

G4PolynomialMultiply := function(left, right)
    local ans, a, b, product;
    ans := [];
    for a in left do
        for b in right do
            product := [[G4MonomialMultiply(a[1], b[1]), a[2] * b[2]]];
            ans := G4PolynomialAdd(ans, product, 1);
        od;
    od;
    return ans;
end;

G4MonomialPolynomial := function(monomial, coefficient)
    if coefficient = 0 then
        return [];
    fi;
    return [[CanonicalSparse(monomial), coefficient]];
end;

G4ZeroMatrix := function(rows, columns)
    return List([1 .. rows], i -> List([1 .. columns], j -> []));
end;

G4IdentityMatrix := function(dimension)
    local matrix, i;
    matrix := G4ZeroMatrix(dimension, dimension);
    for i in [1 .. dimension] do
        matrix[i][i] := G4MonomialPolynomial([], 1);
    od;
    return matrix;
end;

G4DiagonalMatrix := function(monomials)
    local matrix, i;
    matrix := G4ZeroMatrix(Length(monomials), Length(monomials));
    for i in [1 .. Length(monomials)] do
        matrix[i][i] := G4MonomialPolynomial(monomials[i], 1);
    od;
    return matrix;
end;

G4MatrixMultiply := function(left, right)
    local rows, middle, columns, ans, i, k, j, product;
    rows := Length(left);
    middle := Length(right);
    columns := Length(right[1]);
    ans := G4ZeroMatrix(rows, columns);
    for i in [1 .. rows] do
        for k in [1 .. middle] do
            if Length(left[i][k]) > 0 then
                for j in [1 .. columns] do
                    if Length(right[k][j]) > 0 then
                        product := G4PolynomialMultiply(left[i][k], right[k][j]);
                        ans[i][j] := G4PolynomialAdd(ans[i][j], product, 1);
                    fi;
                od;
            fi;
        od;
    od;
    return ans;
end;

G4MatrixApply := function(matrix, vector)
    local ans, i, j, product;
    ans := List([1 .. Length(matrix)], i -> []);
    for i in [1 .. Length(matrix)] do
        for j in [1 .. Length(vector)] do
            if Length(matrix[i][j]) > 0 and Length(vector[j]) > 0 then
                product := G4PolynomialMultiply(matrix[i][j], vector[j]);
                ans[i] := G4PolynomialAdd(ans[i], product, 1);
            fi;
        od;
    od;
    return ans;
end;

G4Kronecker := function(left, right)
    local ans, i, j, k, ell, row, column;
    ans := G4ZeroMatrix(Length(left) * Length(right),
                        Length(left[1]) * Length(right[1]));
    for i in [1 .. Length(left)] do
        for j in [1 .. Length(left[1])] do
            if Length(left[i][j]) > 0 then
                for k in [1 .. Length(right)] do
                    for ell in [1 .. Length(right[1])] do
                        if Length(right[k][ell]) > 0 then
                            row := (i - 1) * Length(right) + k;
                            column := (j - 1) * Length(right[1]) + ell;
                            ans[row][column] :=
                                G4PolynomialMultiply(left[i][j], right[k][ell]);
                        fi;
                    od;
                od;
            fi;
        od;
    od;
    return ans;
end;

G4BasisVector := function(dimension, index)
    local vector;
    vector := List([1 .. dimension], i -> []);
    vector[index] := G4MonomialPolynomial([], 1);
    return vector;
end;

#############################################################################
## Homogeneous objects, tensor action, braiding, and adjoint recursion.
#############################################################################

G4MakeObject := function(degrees, builder)
    return rec(
        degrees := degrees,
        dim := Length(degrees),
        builder := builder,
        cacheActors := [],
        cacheMatrices := []
    );
end;

G4ActionMatrix := function(object, actor)
    local position, matrix;
    position := Position(object.cacheActors, actor);
    if position <> fail then
        return object.cacheMatrices[position];
    fi;
    matrix := object.builder(actor);
    Add(object.cacheActors, ShallowCopy(actor));
    Add(object.cacheMatrices, matrix);
    return matrix;
end;

G4HomogeneousObject := function(degrees, name)
    local builder;
    builder := function(actor)
        local matrix, column, degree, target, row;
        matrix := G4ZeroMatrix(Length(degrees), Length(degrees));
        for column in [1 .. Length(degrees)] do
            degree := degrees[column];
            target := G4Conjugate(actor, degree);
            row := Position(degrees, target);
            if row = fail then
                Gamma4Fail("support is not stable under a required actor");
            fi;
            matrix[row][column] := G4MonomialPolynomial(
                G4ActionVariable(name, actor, degree), 1
            );
        od;
        return matrix;
    end;
    return G4MakeObject(degrees, builder);
end;

G4TensorObject := function(left, right)
    local degrees, builder;
    degrees := Concatenation(List(left.degrees,
        x -> List(right.degrees, y -> G4Multiply(x, y))));
    builder := function(actor)
        local matrix, leftAction, rightAction, i, j, column, correction,
              ii, jj, row, entry;
        matrix := G4ZeroMatrix(Length(degrees), Length(degrees));
        leftAction := G4ActionMatrix(left, actor);
        rightAction := G4ActionMatrix(right, actor);
        for i in [1 .. left.dim] do
            for j in [1 .. right.dim] do
                column := (i - 1) * right.dim + j;
                correction := G4FVariable(actor, left.degrees[i], right.degrees[j]);
                for ii in [1 .. left.dim] do
                    if Length(leftAction[ii][i]) > 0 then
                        for jj in [1 .. right.dim] do
                            if Length(rightAction[jj][j]) > 0 then
                                row := (ii - 1) * right.dim + jj;
                                entry := G4PolynomialMultiply(
                                    G4MonomialPolynomial(correction, 1),
                                    G4PolynomialMultiply(leftAction[ii][i],
                                                         rightAction[jj][j])
                                );
                                matrix[row][column] := entry;
                            fi;
                        od;
                    fi;
                od;
            od;
        od;
        return matrix;
    end;
    return G4MakeObject(degrees, builder);
end;

G4Braid := function(left, right)
    local matrix, i, action, j, jj, column, row;
    matrix := G4ZeroMatrix(right.dim * left.dim, left.dim * right.dim);
    for i in [1 .. left.dim] do
        action := G4ActionMatrix(right, left.degrees[i]);
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

G4PhiOne := function(left, right)
    local dimension, double, ans, i, j, identityEntry;
    dimension := left.dim * right.dim;
    double := G4MatrixMultiply(G4Braid(right, left), G4Braid(left, right));
    ans := G4ZeroMatrix(dimension, dimension);
    for i in [1 .. dimension] do
        for j in [1 .. dimension] do
            if i = j then
                identityEntry := G4MonomialPolynomial([], 1);
            else
                identityEntry := [];
            fi;
            ans[i][j] := G4PolynomialAdd(identityEntry, double[i][j], -1);
        od;
    od;
    return ans;
end;

G4AdjointOperator := function(acting, previous, earlier, previousPhi)
    local current, double, assocMonomials, assocIn, assocOut, c12, correction,
          dimension, phi, i, j, identityEntry;
    current := G4TensorObject(acting, previous);
    double := G4MatrixMultiply(G4Braid(previous, acting),
                               G4Braid(acting, previous));

    ## Corrected convention:
    ## c12 = a (c tensor id) a^(-1).  Thus the input factor is Phi and
    ## the output factor is Phi^(-1).
    assocMonomials := Concatenation(List(acting.degrees,
        x -> Concatenation(List(acting.degrees,
            y -> List(earlier.degrees, z -> G4PVariable(x, y, z))))));
    assocIn := G4DiagonalMatrix(assocMonomials);
    assocOut := G4DiagonalMatrix(List(assocMonomials,
                                      monomial -> G4MonomialInverse(monomial)));
    c12 := G4MatrixMultiply(
        assocOut,
        G4MatrixMultiply(
            G4Kronecker(G4Braid(acting, acting),
                        G4IdentityMatrix(earlier.dim)),
            assocIn
        )
    );
    correction := G4MatrixMultiply(
        G4Kronecker(G4IdentityMatrix(acting.dim), previousPhi), c12
    );

    dimension := current.dim;
    phi := G4ZeroMatrix(dimension, dimension);
    for i in [1 .. dimension] do
        for j in [1 .. dimension] do
            if i = j then
                identityEntry := G4MonomialPolynomial([], 1);
            else
                identityEntry := [];
            fi;
            phi[i][j] := G4PolynomialAdd(
                G4PolynomialAdd(identityEntry, double[i][j], -1),
                correction[i][j], 1
            );
        od;
    od;
    return [current, phi];
end;

G4ThirdAdjointRelations := function()
    local u, r, vDegrees, wDegrees, v, w, d1, f1, wi, vi, y1,
          pair, d2, f2, w2, inputY2, j, y2, d3, f3, gy2,
          vectors, relations, vectorData, tag, vector, outer, inputY3,
          output, coordinate, polynomial, first, second, ratio, sign;
    u := G4Power(G4Epsilon, 2);
    r := G4Power(G4Epsilon, -1);
    vDegrees := [G4H, G4Multiply(r, G4H)];
    wDegrees := [
        G4G,
        G4Multiply(G4Epsilon, G4G),
        G4Multiply(u, G4G),
        G4Multiply(G4Power(G4Epsilon, 3), G4G)
    ];
    v := G4HomogeneousObject(vDegrees, "R");
    w := G4HomogeneousObject(wDegrees, "S");

    d1 := G4TensorObject(w, v);
    f1 := G4PhiOne(w, v);
    wi := Position(w.degrees, G4Multiply(G4Epsilon, G4G));
    vi := Position(v.degrees, G4H);
    y1 := G4MatrixApply(f1,
        G4BasisVector(d1.dim, (wi - 1) * v.dim + vi));

    pair := G4AdjointOperator(w, d1, v, f1);
    d2 := pair[1]; f2 := pair[2];
    w2 := Position(w.degrees, G4Multiply(u, G4G));
    inputY2 := List([1 .. d2.dim], i -> []);
    for j in [1 .. Length(y1)] do
        inputY2[(w2 - 1) * d1.dim + j] := y1[j];
    od;
    y2 := G4MatrixApply(f2, inputY2);

    pair := G4AdjointOperator(w, d2, d1, f2);
    d3 := pair[1]; f3 := pair[2];
    gy2 := G4MatrixApply(G4ActionMatrix(d2, G4G), y2);
    vectors := [rec(tag := "Y3:y2", value := y2),
                rec(tag := "Y3:gy2", value := gy2)];

    relations := [];
    for vectorData in vectors do
        tag := vectorData.tag;
        vector := vectorData.value;
        for outer in [1 .. w.dim] do
            inputY3 := List([1 .. d3.dim], i -> []);
            for j in [1 .. Length(vector)] do
                inputY3[(outer - 1) * d2.dim + j] := vector[j];
            od;
            output := G4MatrixApply(f3, inputY3);
            for coordinate in [1 .. Length(output)] do
                polynomial := output[coordinate];
                if Length(polynomial) > 0 then
                    if Length(polynomial) <> 2 then
                        Gamma4Fail("a reconstructed Y3 coordinate is not binomial");
                    fi;
                    first := polynomial[1];
                    second := polynomial[2];
                    if AbsInt(first[2]) <> 1 or AbsInt(second[2]) <> 1 then
                        Gamma4Fail("a reconstructed Y3 coefficient is not a unit sign");
                    fi;
                    ratio := G4MonomialMultiply(
                        first[1], G4MonomialInverse(second[1])
                    );
                    if first[2] = -second[2] then
                        sign := 0;
                    else
                        sign := 1;
                    fi;
                    Add(relations, rec(
                        label := Concatenation(tag, ":w", String(outer - 1),
                                               ":c", String(coordinate - 1)),
                        sign := sign,
                        ratio := ratio
                    ));
                fi;
            od;
        od;
    od;
    return relations;
end;

VerifyReconstructedRows := function()
    local actual, expected, matches, row;
    actual := G4ThirdAdjointRelations();
    if Length(actual) <> 192 then
        Gamma4Fail("FAIL: expected 192 nonzero reconstructed Y3 binomials, got ",
              Length(actual));
    fi;
    for expected in LocalRelations do
        matches := Filtered(actual, row -> row.label = expected.label);
        if Length(matches) <> 1 then
            Gamma4Fail("FAIL: reconstructed coordinate label missing or repeated: ",
                  expected.label);
        fi;
        row := matches[1];
        if row.sign <> expected.sign then
            Gamma4Fail("FAIL: reconstructed binomial sign differs at ", expected.label);
        fi;
        if row.ratio <> CanonicalSparse(expected.terms) then
            Gamma4Fail("FAIL: reconstructed Laurent row differs at ", expected.label);
        fi;
    od;
    Print("Actual phi3 recursion: nonzero binomial coordinates=", Length(actual),
          ", certificate rows checked in the same orientation=",
          Length(LocalRelations), "\n");
end;

RelationRow := function(relation)
    local ans;
    ans := CanonicalSparse(relation.terms);
    AddExponent(ans, "SIGN", -relation.sign);
    return CanonicalSparse(ans);
end;

RelationByLabel := function(label)
    local matches;
    matches := Filtered(LocalRelations, relation -> relation.label = label);
    if Length(matches) <> 1 then
        Gamma4Fail("relation label is missing or repeated: ", label);
    fi;
    return matches[1];
end;

VerifyCertificate := function(targetName, target, certificate)
    local combination, item, relation, residual, signExponent;
    combination := [];
    signExponent := 0;
    for item in certificate do
        relation := RelationByLabel(item[1]);
        combination := AddSparse(
            combination,
            ScaleSparse(RelationRow(relation), item[2])
        );
        signExponent := signExponent + item[2] * relation.sign;
    od;
    ## The printed target is Delta_ij^W, whereas the certificate is the
    ## exponent vector of (Delta_ij^W)^(-1).  Thus target + combination
    ## must vanish.
    residual := AddSparse(CanonicalSparse(target), combination);
    if Length(residual) <> 0 then
        Gamma4Fail("FAIL: nonzero Laurent residual for ", targetName);
    fi;
    if signExponent <> 0 then
        Gamma4Fail("FAIL: the signs do not cancel for ", targetName);
    fi;
    Print(targetName, ": relations=", Length(certificate),
          ", coefficient sum=", Sum(certificate, item -> item[2]),
          ", sign exponent=", signExponent,
          ", residual atoms=", Length(residual), "\n");
end;

VerifyReconstructedRows();
VerifyCertificate("(Delta_01^W)^(-1)", W01Target, W01Certificate);
VerifyCertificate("(Delta_10^W)^(-1)", W10Target, W10Certificate);

Print("Gamma4 local-Delta Laurent certificate\n");
Print("Arithmetic: exact integers in the free Laurent exponent lattice\n");
Print("External GAP packages required: none\n");
Print("c12 associator direction: input Phi, output Phi^-1\n");
Print("Corrected Y3 certificate rows independently checked: ",
      Length(LocalRelations), "\n");
Print("PASS: Y3=0 forces Delta_01^W=Delta_10^W=1\n");

QUIT_GAP(0);
