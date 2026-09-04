Read("bar_complex.g");

# Paper correspondence:
#   eq:Gamma3-dim2-shifted-X1 through
#   eq:Gamma3-dim2-three-path-vector and eq:Gamma3-dim2-eta-values;
#   eq:Gamma3-dim2-E-chain, eq:Gamma3-dim2-R-reduction,
#   eq:Gamma3-dim2-A-chain through eq:Gamma3-dim2-M-chains, and
#   eq:Gamma3-common-three-path-reduction.

VerifyGamma331Dim2 := function()
    local one, e, e2, g, h, G3Mul, G3Inv, G3Conj,
          Bar, AddC, SubC, Sc, SumC, L, U,
          a, abar, s0, s1, s2, z, Ecoeff, A, B,
          T0, T1, T2, M0, M1, M2, Theta, K,
          Ceigen, Rred, RA, C02, Cpow, expected,
          SOne, S3Mul, Project, Proj, SBar,
          x, x2, y, xy, x2y, D02, Dpow;

    one := [0, 0, 0];
    e   := [1, 0, 0];
    e2  := [2, 0, 0];
    g   := [0, 1, 0];
    h   := [0, 0, 1];

    G3Mul := function(u, v)
        local sign;
        if u[2] mod 2 = 0 then sign := 1; else sign := -1; fi;
        return [(u[1] + sign * v[1]) mod 3,
                u[2] + v[2], u[3] + v[3]];
    end;
    G3Inv := function(u)
        if u[2] mod 2 = 0 then
            return [(-u[1]) mod 3, -u[2], -u[3]];
        else
            return [u[1] mod 3, -u[2], -u[3]];
        fi;
    end;
    G3Conj := function(u, v)
        return G3Mul(G3Mul(u, v), G3Inv(u));
    end;

    Bar  := function(xs) return BC_Bar(one, xs); end;
    AddC := function(u, v) return BC_Add(u, v, one); end;
    SubC := function(u, v) return BC_Sub(u, v, one); end;
    Sc   := function(n, u) return BC_Scale(n, u, one); end;
    SumC := function(xs) return BC_Sum(xs, one); end;
    L := function(x0, u, v)
        return SumC([
            Bar([u, v, x0]),
            Bar([G3Conj(G3Mul(u, v), x0), u, v]),
            Sc(-1, Bar([u, G3Conj(v, x0), v]))
        ]);
    end;
    U := function(u, x0, y0)
        return SumC([
            Bar([u, x0, y0]),
            Bar([G3Conj(u, x0), G3Conj(u, y0), u]),
            Sc(-1, Bar([G3Conj(u, x0), u, y0]))
        ]);
    end;

    a    := G3Mul(e, g);
    abar := G3Mul(e2, g);
    s0   := G3Mul(g, h);
    s1   := G3Conj(e, s0);
    s2   := G3Conj(e2, s0);
    z    := G3Mul(e, G3Mul(G3Mul(g, g), h));

    Ecoeff := SumC([L(h, e, g), Sc(-1, L(h, g, e2)),
                    Sc(-1, L(h, e, e))]);
    A := SubC(L(g, s1, e2), L(s0, e, g));
    B := SumC([Bar([a, g, h]), Sc(-1, Bar([abar, a, h])),
               Sc(2, Ecoeff), L(h, e, e), Sc(-1, L(g, e, g)),
               Sc(-1, U(e2, g, h))]);

    T0 := AddC(U(e, a, s0), L(g, e, e2));
    T1 := AddC(U(e, g, s1), L(s0, e, e));
    T2 := SumC([U(e, abar, s2), L(g, e, e), L(s0, e, e2)]);

    M0 := SubC(T0, A);
    M1 := SubC(AddC(A, T1), B);
    M2 := AddC(B, T2);
    Theta := AddC(L(z, e, e), L(z, e2, e));
    K := AddC(Bar([e, e, e]), Bar([e, e2, e]));

    Ceigen := SumC([Sc(3, Ecoeff), L(h, e, e), L(h, e2, e)]);
    Rred := SumC([L(g, h, e2), Sc(-1, L(g, e2, h)),
                  Sc(-2, L(h, e, g)), L(h, g, e2), L(h, e, e)]);
    RA := SumC([U(a, g, h), L(s0, e, g), Sc(-1, L(g, e, g)),
                Sc(-1, L(h, e, g)), Sc(-1, L(g, g, h)),
                Sc(-1, U(e, g, h))]);
    C02 := SumC([M0, Sc(-1, M2), Sc(-1, K)]);
    Cpow := SumC([Sc(3, M0), Sc(-1, Theta), K]);

    expected := SubC(Bar([e, z]), Bar([z, e]));
    BC_AssertEqual("Gamma31dim2 dM0", BC_Boundary(M0, G3Mul, one), expected, one);
    BC_AssertEqual("Gamma31dim2 dM1", BC_Boundary(M1, G3Mul, one), expected, one);
    BC_AssertEqual("Gamma31dim2 dM2", BC_Boundary(M2, G3Mul, one), expected, one);
    BC_AssertEqual("Gamma31dim2 dTheta", BC_Boundary(Theta, G3Mul, one),
                   Sc(3, expected), one);
    BC_AssertZero("Gamma31dim2 Ceigen cycle", BC_Boundary(Ceigen, G3Mul, one), one);
    BC_AssertZero("Gamma31dim2 R cycle", BC_Boundary(Rred, G3Mul, one), one);
    BC_AssertZero("Gamma31dim2 RA cycle", BC_Boundary(RA, G3Mul, one), one);
    BC_AssertZero("Gamma31dim2 C02 cycle", BC_Boundary(C02, G3Mul, one), one);
    BC_AssertZero("Gamma31dim2 Cpow cycle", BC_Boundary(Cpow, G3Mul, one), one);
    if Length(C02) <> 38 or Length(Cpow) <> 20 then
        G3_Fail("Gamma31dim2 chain term counts differ from the recorded certificate");
    fi;

    SOne := [0, 0, 0];
    S3Mul := function(u, v)
        local sign;
        if u[2] mod 2 = 0 then sign := 1; else sign := -1; fi;
        return [(u[1] + sign * v[1]) mod 3,
                (u[2] + v[2]) mod 2, 0];
    end;
    Project := function(u) return [u[1], u[2] mod 2, 0]; end;
    Proj := function(chain)
        return BC_Project(chain, Project, one, SOne);
    end;
    SBar := function(xs) return BC_Bar(SOne, xs); end;

    x   := Project(e);
    x2  := Project(e2);
    y   := Project(g);
    xy  := S3Mul(x, y);
    x2y := S3Mul(x2, y);

    D02 := BC_Sum([
        SBar([y,x2y,x,x]),
        SBar([x,x2y,x,y]),
        SBar([xy,x,xy,x]),
        BC_Scale(-1,SBar([xy,xy,x,x]),SOne),
        BC_Scale(-1,SBar([xy,xy,x2,x]),SOne),
        BC_Scale(-1,SBar([x2y,x,x,y]),SOne),
        BC_Scale(-1,SBar([x2y,x,xy,x]),SOne),
        SBar([x2y,xy,x2,x])
    ], SOne);

    Dpow := BC_Sum([
        SBar([y,y,y,y]),
        SBar([y,x2y,x,y]),
        BC_Scale(-2,SBar([y,x2y,xy,y]),SOne),
        BC_Scale( 3,SBar([y,x2y,xy,x2]),SOne),
        BC_Scale(-3,SBar([y,x2y,x2,y]),SOne),
        SBar([x,y,y,y]),
        SBar([x,x2y,x,y]),
        SBar([x,x2y,xy,y]),
        SBar([x2y,xy,y,y])
    ], SOne);

    BC_AssertZero("Gamma31dim2 projected Ceigen", Proj(Ceigen), SOne);
    BC_AssertZero("Gamma31dim2 projected R", Proj(Rred), SOne);
    BC_AssertZero("Gamma31dim2 projected RA", Proj(RA), SOne);
    BC_AssertEqual("Gamma31dim2 projected C02 boundary",
        BC_Boundary(D02, S3Mul, SOne), Proj(C02), SOne);
    BC_AssertEqual("Gamma31dim2 projected Cpow boundary",
        BC_Boundary(Dpow, S3Mul, SOne), Proj(Cpow), SOne);

    Print("[PASS] (3,1)_2: all Gamma_3 differentials and reductions agree.\n");
    Print("[PASS] (3,1)_2: explicit S_3 four-chain residuals are zero.\n");
    Print("       C02/Cpow Gamma terms: 38, 20; S3 four-chain terms: 8, 9.\n");
    return true;
end;

G3_GAMMA331_DIM2_OK := VerifyGamma331Dim2();


# The preceding part proves the bar-chain differentials and the two explicit
# S_3 boundary certificates.  The following independent symbolic calculation
# checks the two interfaces which were formerly left to hand calculation:
#
#   (i) the three summands in the recursive formula for varphi_2 on e_0,
#       including the shifted varphi_1 calculation and the cyclic epsilon
#       action on e_0,e_1,e_2;
#
#   (ii) the evaluation of E,R,A,B,T_i,M_i,Theta,K and the two comparison
#        cycles by a normalized multiplicative 3-cocycle.
#
# Atomic values Phi(a,b,c) are independent variables.  Phi_x(a,b) and
# Phi^a(x,y) are expanded from their definitions, rather than entered as
# independent names.  Thus the target coefficient and chain-evaluation
# formulas are not supplied as hypotheses.

VerifyGamma331Dim2CoefficientAndBridge := function()
    local one, e, e2, g, h, gh, eg, e2g, z,
          G3Mul, G3Inv, G3Conj, AssertGroup,
          coefficientField, ring, vars, onePolynomial,
          symbolKeys, nextSymbol, Q, Atom, RatZero, RatOne,
          RatAdd, RatNeg, RatMul, RatProduct, RatInv, RatDiv,
          RatPow, RatEq, AssertRat, RejectRat, Phi, PhiSub, PhiSup,
          Bar, Sc, SumC, L, U, EvalChain,
          VecAdd, VecScale, AssertVec,
          rhoG, rhoH, b, omega, aEpsilon, rRaw, rReduced, rChain,
          phiHEG, phiHEE, shiftedScale, x12Coordinates,
          phiW, phiWbar, shiftedDifference, C0,
          frakA, frakB, rhoGH, actionX1, actionV2,
          doubleBraiding, c12Coefficient, lastPath, x2Paths,
          TransportStep, TensorStep, cyclic0, cyclic1, cyclic2,
          eta0, eta1, eta2, eta0Target, eta1Target, eta2Target,
          m0, m1, m2,
          degreeV0, degreeV1, degreeV2,
          degreeX10, degreeX11, degreeX12, commonDegree,
          Ecoeff, Achain, Bchain, T0chain, T1chain, T2chain,
          M0chain, M1chain, M2chain, ThetaChain, Kchain,
          C02chain, CpowChain, theta, kappa,
          formalM0, formalK, formalTheta, derivedM2, derivedM1,
          expectedM1, commonQuotient, extraFactor, perturbedC12,
          perturbedLastPath;

    #########################################################################
    # Gamma_3 and exact rational functions
    #########################################################################

    one := [0, 0, 0];
    e   := [1, 0, 0];
    e2  := [2, 0, 0];
    g   := [0, 1, 0];
    h   := [0, 0, 1];

    G3Mul := function(x, y)
        local sign;
        if x[2] mod 2 = 0 then sign := 1; else sign := -1; fi;
        return [(x[1] + sign * y[1]) mod 3,
                x[2] + y[2], x[3] + y[3]];
    end;
    G3Inv := function(x)
        if x[2] mod 2 = 0 then
            return [(-x[1]) mod 3, -x[2], -x[3]];
        fi;
        return [x[1] mod 3, -x[2], -x[3]];
    end;
    G3Conj := function(a, x)
        return G3Mul(G3Mul(a, x), G3Inv(a));
    end;
    AssertGroup := function(label, left, right)
        if left <> right then
            G3_Fail(label, ": group elements differ: ", left, " versus ", right);
        fi;
    end;

    gh  := G3Mul(g, h);
    eg  := G3Mul(e, g);
    e2g := G3Mul(e2, g);
    z   := G3Mul(e, G3Mul(G3Mul(g, g), h));

    # Coefficients contain a primitive third root of unity.  All cocycle
    # values and representation scalars remain algebraically independent.
    coefficientField := CF(3);
    ring := PolynomialRing(coefficientField, 180);
    vars := IndeterminatesOfPolynomialRing(ring);
    onePolynomial := One(ring);
    symbolKeys := [];
    nextSymbol := 0;

    Q := function(c)
        return [c * onePolynomial, onePolynomial];
    end;
    RatZero := Q(0);
    RatOne := Q(1);
    Atom := function(key)
        local pos;
        pos := Position(symbolKeys, key);
        if pos = fail then
            nextSymbol := nextSymbol + 1;
            if nextSymbol > Length(vars) then
                G3_Fail("not enough symbolic variables");
            fi;
            Add(symbolKeys, key);
            pos := nextSymbol;
        fi;
        return [vars[pos], onePolynomial];
    end;
    RatAdd := function(x, y)
        return [x[1] * y[2] + y[1] * x[2], x[2] * y[2]];
    end;
    RatNeg := function(x)
        return [-x[1], x[2]];
    end;
    RatMul := function(x, y)
        return [x[1] * y[1], x[2] * y[2]];
    end;
    RatProduct := function(list)
        local out, x;
        out := RatOne;
        for x in list do out := RatMul(out, x); od;
        return out;
    end;
    RatInv := function(x)
        return [x[2], x[1]];
    end;
    RatDiv := function(x, y)
        return RatMul(x, RatInv(y));
    end;
    RatPow := function(x, n)
        local out, i;
        if n < 0 then return RatPow(RatInv(x), -n); fi;
        out := RatOne;
        if n = 0 then return out; fi;
        for i in [1..n] do out := RatMul(out, x); od;
        return out;
    end;
    RatEq := function(x, y)
        return x[1] * y[2] = y[1] * x[2];
    end;
    AssertRat := function(label, left, right)
        if not RatEq(left, right) then
            G3_Fail(label, ": nonzero exact rational-function residual");
        fi;
    end;
    RejectRat := function(label, left, right)
        if RatEq(left, right) then
            G3_Fail(label, ": a deliberately perturbed source was accepted");
        fi;
    end;

    Phi := function(a, c, d)
        if a = one or c = one or d = one then return RatOne; fi;
        return Atom(["Phi", a, c, d]);
    end;
    PhiSub := function(x, a, c)
        return RatDiv(
            RatMul(Phi(a, c, x), Phi(G3Conj(G3Mul(a, c), x), a, c)),
            Phi(a, G3Conj(c, x), c));
    end;
    PhiSup := function(a, x, y)
        return RatDiv(
            RatMul(Phi(a, x, y), Phi(G3Conj(a, x), G3Conj(a, y), a)),
            Phi(G3Conj(a, x), a, y));
    end;

    #########################################################################
    # Independent evaluation of normalized bar chains
    #########################################################################

    Bar := function(entries) return BC_Bar(one, entries); end;
    Sc := function(n, chain) return BC_Scale(n, chain, one); end;
    SumC := function(chains) return BC_Sum(chains, one); end;
    L := function(x, a, c)
        return SumC([
            Bar([a, c, x]),
            Bar([G3Conj(G3Mul(a, c), x), a, c]),
            Sc(-1, Bar([a, G3Conj(c, x), c]))
        ]);
    end;
    U := function(a, x, y)
        return SumC([
            Bar([a, x, y]),
            Bar([G3Conj(a, x), G3Conj(a, y), a]),
            Sc(-1, Bar([G3Conj(a, x), a, y]))
        ]);
    end;
    EvalChain := function(chain)
        local out, term;
        out := RatOne;
        for term in chain do
            out := RatMul(out,
                RatPow(Phi(term[1][1], term[1][2], term[1][3]), term[2]));
        od;
        return out;
    end;

    VecAdd := function(left, right)
        local out, i;
        out := [];
        for i in [1..Length(left)] do
            Add(out, RatAdd(left[i], right[i]));
        od;
        return out;
    end;
    VecScale := function(scalar, vector)
        return List(vector, x -> RatMul(scalar, x));
    end;
    AssertVec := function(label, left, right)
        local i;
        if Length(left) <> Length(right) then
            G3_Fail(label, ": vector lengths differ");
        fi;
        for i in [1..Length(left)] do
            AssertRat(Concatenation(label, " coordinate ", String(i)),
                      left[i], right[i]);
        od;
    end;

    #########################################################################
    # Coefficient calculation from the action table and the recursion
    #########################################################################

    rhoG := Q(-1);
    rhoH := Atom(["representation", "rho(h)"]);
    # Delta_3=1 is rho(h)^2 b=1.
    b := RatInv(RatPow(rhoH, 2));
    omega := Q(E(3));

    phiHEG := PhiSub(h, e, g);
    phiHEE := PhiSub(h, e, e);
    aEpsilon := RatDiv(phiHEG,
        RatMul(PhiSub(h, g, e2), phiHEE));

    # The raw local-cocycle quotient is reduced only after checking that its
    # quotient by the claimed reduced form is exactly <Phi,R>.
    rRaw := RatDiv(PhiSub(g, h, e2), PhiSub(g, e2, h));
    rReduced := RatMul(aEpsilon, phiHEG);
    rChain := SumC([
        L(g, h, e2), Sc(-1, L(g, e2, h)),
        Sc(-2, L(h, e, g)), L(h, g, e2), L(h, e, e)
    ]);
    AssertRat("R-chain evaluation",
              EvalChain(rChain), RatDiv(rRaw, rReduced));

    C0 := RatDiv(RatMul(phiHEE, RatPow(aEpsilon, 2)),
                 PhiSup(e2, g, h));
    shiftedScale := RatDiv(PhiSup(e2, g, h),
                           RatMul(RatPow(aEpsilon, 2), phiHEE));
    x12Coordinates := [
        RatMul(shiftedScale, RatPow(omega, 2)),
        RatNeg(RatProduct([shiftedScale, rhoH, omega]))
    ];

    # These two vectors are obtained directly from the two rows of the
    # action table for w and wbar, after the verified R-chain reduction.
    phiW := [
        RatOne,
        RatNeg(RatDiv(
            RatProduct([rhoH, rReduced, RatPow(omega, 2)]),
            RatMul(aEpsilon, phiHEG)))
    ];
    phiWbar := [
        RatNeg(RatDiv(
            RatProduct([rhoH, rReduced, b, omega]),
            RatMul(aEpsilon, phiHEG))),
        RatOne
    ];
    shiftedDifference := VecAdd(phiW,
        VecScale(RatNeg(rhoH), phiWbar));

    AssertVec("shifted varphi_1 on w", phiW,
        VecScale(RatMul(omega, C0), x12Coordinates));
    AssertVec("shifted varphi_1 on wbar", phiWbar,
        VecScale(RatNeg(RatDiv(RatMul(RatPow(omega, 2), C0), rhoH)),
                 x12Coordinates));
    AssertVec("shifted varphi_1 on w-rho(h)wbar", shiftedDifference,
        VecScale(RatNeg(C0), x12Coordinates));

    frakA := RatDiv(PhiSub(g, G3Conj(e, gh), e2),
                    PhiSub(gh, e, g));
    frakB := RatDiv(
        RatProduct([
            Phi(eg, g, h), phiHEE, RatPow(aEpsilon, 2)
        ]),
        RatProduct([
            Phi(e2g, eg, h), PhiSub(g, e, g), PhiSup(e2, g, h)
        ]));

    # Source coefficients for the two braidings.  The equality
    # Phi^g(g,h)=Phi_g(g,h) is checked from the atomic definitions.
    AssertRat("commuting local/tensor cocycle",
              PhiSup(g, g, h), PhiSub(g, g, h));
    rhoGH := RatDiv(RatMul(rhoG, rhoH), PhiSub(g, g, h));
    actionX1 := RatNeg(RatDiv(
        RatProduct([PhiSup(g, g, h), rhoG]),
        RatMul(rhoH, PhiSub(gh, e, g))));
    actionV2 := RatMul(PhiSub(g, G3Conj(e, gh), e2), rhoGH);
    doubleBraiding := RatMul(actionX1, actionV2);
    AssertRat("minus double braiding path",
              RatNeg(doubleBraiding), frakA);

    c12Coefficient := RatDiv(
        RatMul(rhoG, Phi(eg, g, h)),
        RatProduct([Phi(e2g, eg, h), PhiSub(g, e, g)]));
    lastPath := RatMul(c12Coefficient, RatNeg(C0));
    AssertRat("last recursive path", lastPath, frakB);

    x2Paths := VecAdd(
        VecAdd([RatOne, RatZero, RatZero],
               [RatZero, RatNeg(doubleBraiding), RatZero]),
        [RatZero, RatZero, lastPath]);
    AssertVec("three-path vector", x2Paths, [RatOne, frakA, frakB]);

    degreeV0 := g;
    degreeV1 := G3Conj(e, g);
    degreeV2 := G3Conj(e2, g);
    degreeX10 := gh;
    degreeX11 := G3Conj(e, gh);
    degreeX12 := G3Conj(e2, gh);
    commonDegree := z;
    AssertGroup("degree e_0", G3Mul(degreeV2, degreeX10), commonDegree);
    AssertGroup("degree e_1", G3Mul(degreeV0, degreeX11), commonDegree);
    AssertGroup("degree e_2", G3Mul(degreeV1, degreeX12), commonDegree);

    # Reconstruct epsilon on a transported vector e^i action u from
    # eq:YD-projective-action, and combine the two factors using
    # eq:tensor-product.  The output is [new transversal index,scalar].
    # No eta_i formula is used by this routine.
    TransportStep := function(baseDegree, i)
        if i = 0 then return [1, RatOne]; fi;
        if i = 1 then return [2, PhiSub(baseDegree, e, e)]; fi;
        if i = 2 then return [0, PhiSub(baseDegree, e, e2)]; fi;
        G3_Fail("transport index outside {0,1,2}");
    end;
    TensorStep := function(i, j)
        local vStep, xStep, degreeVi, degreeXj;
        vStep := TransportStep(g, i);
        xStep := TransportStep(gh, j);
        degreeVi := G3Conj([i mod 3, 0, 0], g);
        degreeXj := G3Conj([j mod 3, 0, 0], gh);
        return [vStep[1], xStep[1],
                RatProduct([
                    PhiSup(e, degreeVi, degreeXj), vStep[2], xStep[2]
                ])];
    end;

    cyclic0 := TensorStep(2, 0);
    cyclic1 := TensorStep(0, 1);
    cyclic2 := TensorStep(1, 2);
    if cyclic0{[1,2]} <> [0,1] or cyclic1{[1,2]} <> [1,2]
       or cyclic2{[1,2]} <> [2,0] then
        G3_Fail("epsilon does not cycle e_0,e_1,e_2 with the recorded indices");
    fi;
    eta0 := cyclic0[3];
    eta1 := cyclic1[3];
    eta2 := cyclic2[3];

    # Compare the independently generated coefficients with the three
    # displayed right-hand sides in eq:Gamma3-dim2-eta-values.
    eta0Target := RatMul(PhiSup(e, eg, gh), PhiSub(g, e, e2));
    eta1Target := RatMul(PhiSup(e, g, G3Conj(e, gh)),
                         PhiSub(gh, e, e));
    eta2Target := RatProduct([
        PhiSup(e, e2g, G3Conj(e2, gh)),
        PhiSub(g, e, e), PhiSub(gh, e, e2)
    ]);
    AssertRat("eta_0 from transported actions", eta0, eta0Target);
    AssertRat("eta_1 from transported actions", eta1, eta1Target);
    AssertRat("eta_2 from transported actions", eta2, eta2Target);
    m0 := RatDiv(eta0, frakA);
    m1 := RatDiv(RatMul(frakA, eta1), frakB);
    m2 := RatMul(frakB, eta2);

    #########################################################################
    # Scalar--chain dictionary and the comparison-cycle reduction
    #########################################################################

    Ecoeff := SumC([L(h, e, g), Sc(-1, L(h, g, e2)),
                    Sc(-1, L(h, e, e))]);
    Achain := SumC([L(g, G3Conj(e, gh), e2),
                    Sc(-1, L(gh, e, g))]);
    Bchain := SumC([
        Bar([eg, g, h]), Sc(-1, Bar([e2g, eg, h])),
        Sc(2, Ecoeff), L(h, e, e), Sc(-1, L(g, e, g)),
        Sc(-1, U(e2, g, h))
    ]);
    T0chain := SumC([U(e, eg, gh), L(g, e, e2)]);
    T1chain := SumC([U(e, g, G3Conj(e, gh)), L(gh, e, e)]);
    T2chain := SumC([
        U(e, e2g, G3Conj(e2, gh)),
        L(g, e, e), L(gh, e, e2)
    ]);
    M0chain := SumC([T0chain, Sc(-1, Achain)]);
    M1chain := SumC([Achain, T1chain, Sc(-1, Bchain)]);
    M2chain := SumC([Bchain, T2chain]);
    ThetaChain := SumC([L(z, e, e), L(z, e2, e)]);
    Kchain := SumC([Bar([e, e, e]), Bar([e, e2, e])]);
    C02chain := SumC([M0chain, Sc(-1, M2chain), Sc(-1, Kchain)]);
    CpowChain := SumC([Sc(3, M0chain), Sc(-1, ThetaChain), Kchain]);
    theta := RatMul(PhiSub(z, e, e), PhiSub(z, e2, e));
    kappa := RatMul(Phi(e, e, e), Phi(e, e2, e));

    AssertRat("<Phi,E>", EvalChain(Ecoeff), aEpsilon);
    AssertRat("<Phi,A>", EvalChain(Achain), frakA);
    AssertRat("<Phi,B>", EvalChain(Bchain), frakB);
    AssertRat("<Phi,T_0>", EvalChain(T0chain), eta0);
    AssertRat("<Phi,T_1>", EvalChain(T1chain), eta1);
    AssertRat("<Phi,T_2>", EvalChain(T2chain), eta2);
    AssertRat("<Phi,M_0>", EvalChain(M0chain), m0);
    AssertRat("<Phi,M_1>", EvalChain(M1chain), m1);
    AssertRat("<Phi,M_2>", EvalChain(M2chain), m2);
    AssertRat("<Phi,Theta>", EvalChain(ThetaChain), theta);
    AssertRat("<Phi,K>", EvalChain(Kchain), kappa);
    AssertRat("C_02 evaluation",
        EvalChain(C02chain), RatDiv(RatDiv(m0, m2), kappa));
    AssertRat("C_pow evaluation",
        EvalChain(CpowChain),
        RatMul(RatDiv(RatPow(m0, 3), theta), kappa));

    # Algebra after the two comparison cycles.  The inputs here are precisely
    # m_0/m_2=kappa, m_0^3=theta/kappa, and eq:YD-projective-action applied
    # three times, which gives m_0*m_1*m_2=theta.  The common reduction is
    # computed, not entered.
    formalM0 := Atom(["comparison", "m_0"]);
    formalK := Atom(["comparison", "kappa"]);
    formalTheta := RatMul(RatPow(formalM0, 3), formalK);
    derivedM2 := RatDiv(formalM0, formalK);
    derivedM1 := RatDiv(formalTheta, RatMul(formalM0, derivedM2));
    expectedM1 := RatDiv(formalM0, formalK);
    commonQuotient := RatDiv(derivedM1, expectedM1);
    AssertRat("common reduction before kappa^3=1",
              commonQuotient, RatPow(formalK, 3));

    # A deliberate perturbation is inserted in the source c_12 path and the
    # full last-path construction is repeated.  The target must be rejected.
    extraFactor := Atom(["negative test", "extra c_12 factor"]);
    perturbedC12 := RatMul(c12Coefficient, extraFactor);
    perturbedLastPath := RatMul(perturbedC12, RatNeg(C0));
    RejectRat("perturbed last recursive path", perturbedLastPath, frakB);

    Print("[PASS] (3,1)_2 coefficients: shifted varphi_1 and three paths agree.\n");
    Print("[PASS] (3,1)_2 action: transported epsilon coefficients agree.\n");
    Print("[PASS] (3,1)_2 chain bridge: E/R/A/B/T/M/Theta/K agree.\n");
    Print("[PASS] (3,1)_2 comparisons: the remaining quotient is kappa^3.\n");
    Print("[PASS] (3,1)_2 negative test: a perturbed c_12 path is rejected.\n");
    Print("       Atomic symbolic generators used: ", Length(symbolKeys), ".\n");
    return true;
end;

if not G3_GAMMA331_DIM2_OK then
    G3_Fail("The (3,1)_2 bar verifier did not report success");
fi;
G3_GAMMA331_DIM2_OK := VerifyGamma331Dim2CoefficientAndBridge();
