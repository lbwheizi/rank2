# The one-dimensional central-support calculation uses exactly the same
# normalized bar chains E,A,B,T0,T1,T2,M0,M1,M2,Theta,K,C02,Cpow as the
# (3,1)_2 calculation.  Reading the shared verifier first rechecks those
# chains and their explicit S3 four-chain boundaries.
# Paper correspondence:
#   lem:Gamma3-dim1-three-elementary-tensors;
#   lem:Gamma3-dim1-chain-identification and
#   eq:Gamma3-dim1-chain-evaluations;
#   eq:Gamma3-dim1-m-mu and eq:Gamma3-dim1-m-ratios.
G3_Fail := function(arg)
    CallFuncList(PrintTo,
        Concatenation(["*errout*"], arg, ["\n"]));
    ErrorNoReturn("Certificate verification failed.");
end;

G3_GAMMA331_DIM2_OK := false;
G3_GAMMA331_DIM1_OK := false;
Read("verify_gamma3_31_dim2.g");
if not G3_GAMMA331_DIM2_OK then
    G3_Fail("The shared (3,1)_2 verifier did not report success");
fi;


# This branch has the same normalized bar chains as the two-dimensional
# branch, but different coefficients in the three recursive paths.  The
# following calculation reconstructs those coefficients from the projective
# action, tensor action, the two braidings, and the defining recursion for
# varphi_2.  It also expands every chain evaluation into atomic values
# Phi(a,b,c), and only then derives the m_i--mu_i dictionary.

VerifyGamma331Dim1CoefficientAndBridge := function()
    local one, e, e2, g, h, gh, eg, e2g, z,
          G3Mul, G3Inv, G3Conj, AssertGroup,
          ring, vars, onePolynomial, symbolKeys, nextSymbol,
          Q, Atom, RatZero, RatOne, RatAdd, RatNeg, RatSub,
          RatMul, RatProduct, RatInv, RatDiv, RatPow, RatEq,
          AssertRat, RejectRat,
          Phi, PhiSub, PhiSup, Bar, Sc, SumC, L, U, EvalChain,
          VecAdd, AssertVec,
          rhoG, rhoH, tauG, tauEpsilon, tauEpsilon2,
          phiHEG, phiHEE, aEpsilon, frakA, frakB,
          rho1G, rhoGH, actionX1, actionV2, doubleBraiding,
          c12Coefficient, shiftedX1Scale, lastPath, Acoeff, Bcoeff,
          x2Paths, TransportStep, TensorStep, cyclic0, cyclic1, cyclic2,
          eta0, eta1, eta2, eta0Target, eta1Target, eta2Target,
          degreeV0, degreeV1, degreeV2,
          degreeX10, degreeX11, degreeX12, commonDegree,
          Ecoeff, Achain, Bchain, T0chain, T1chain, T2chain,
          M0chain, M1chain, M2chain, ThetaChain, Kchain,
          C02chain, CpowChain, ActionChain,
          mu0, mu1, mu2, m0, m1, m2, theta, kappa,
          formalMu0, formalK, formalTheta, derivedMu2, derivedMu1,
          expectedMu1, commonQuotient, extraFactor, perturbedC12,
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

    ring := PolynomialRing(Rationals, 180);
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
    RatNeg := function(x) return [-x[1], x[2]]; end;
    RatSub := function(x, y) return RatAdd(x, RatNeg(y)); end;
    RatMul := function(x, y) return [x[1] * y[1], x[2] * y[2]]; end;
    RatProduct := function(list)
        local out, x;
        out := RatOne;
        for x in list do out := RatMul(out, x); od;
        return out;
    end;
    RatInv := function(x) return [x[2], x[1]]; end;
    RatDiv := function(x, y) return RatMul(x, RatInv(y)); end;
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
    # Coefficients from the two braidings and the varphi_2 recursion
    #########################################################################

    rhoG := Atom(["representation", "rho(g)"]);
    rhoH := Atom(["representation", "rho(h)"]);
    tauG := Atom(["representation", "tau(g)"]);
    phiHEG := PhiSub(h, e, g);
    phiHEE := PhiSub(h, e, e);
    aEpsilon := RatDiv(phiHEG,
        RatMul(PhiSub(h, g, e2), phiHEE));

    # The two already-proved written identities
    # eq:Gamma3-dim1-epsilon-normalization and
    # eq:Gamma3-dim1-epsilon-square are used here as rewrite rules.  They are
    # inputs to this coefficient certificate, not conclusions of it.
    tauEpsilon := RatInv(aEpsilon);
    tauEpsilon2 := RatDiv(RatPow(tauEpsilon, 2), phiHEE);

    frakA := RatDiv(PhiSub(g, G3Conj(e, gh), e2),
                    PhiSub(gh, e, g));
    frakB := RatDiv(Phi(eg, g, h),
        RatProduct([
            Phi(e2g, eg, h), PhiSub(g, e, g),
            PhiSup(e2, g, h), tauEpsilon2
        ]));

    AssertRat("commuting local/tensor cocycle",
              PhiSup(g, g, h), PhiSub(g, g, h));
    rho1G := RatProduct([PhiSup(g, g, h), rhoG, tauG]);
    rhoGH := RatDiv(RatMul(rhoG, rhoH), PhiSub(g, g, h));
    actionX1 := RatDiv(rho1G, PhiSub(gh, e, g));
    actionV2 := RatMul(PhiSub(g, G3Conj(e, gh), e2), rhoGH);
    doubleBraiding := RatMul(actionX1, actionV2);
    Acoeff := RatNeg(doubleBraiding);
    AssertRat("double-braiding recursive coefficient", Acoeff,
        RatNeg(RatProduct([
            RatPow(rhoG, 2), rhoH, tauG, frakA
        ])));

    c12Coefficient := RatDiv(
        RatMul(rhoG, Phi(eg, g, h)),
        RatProduct([Phi(e2g, eg, h), PhiSub(g, e, g)]));
    shiftedX1Scale := RatDiv(
        RatSub(RatOne, RatMul(rhoH, tauG)),
        RatMul(PhiSup(e2, g, h), tauEpsilon2));
    lastPath := RatMul(c12Coefficient, shiftedX1Scale);
    Bcoeff := RatProduct([
        rhoG, RatSub(RatOne, RatMul(rhoH, tauG)), frakB
    ]);
    AssertRat("last recursive coefficient", lastPath, Bcoeff);

    x2Paths := VecAdd(
        VecAdd([RatOne, RatZero, RatZero],
               [RatZero, Acoeff, RatZero]),
        [RatZero, RatZero, lastPath]);
    AssertVec("three elementary tensors", x2Paths,
              [RatOne, Acoeff, Bcoeff]);

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

    # Generate the cyclic action from the chosen transversals.  For a vector
    # e^i action u, the projective-action law supplies both the new index and
    # its scalar; the target eta_i formulas are not inputs to this routine.
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

    #########################################################################
    # Scalar--chain dictionary and the common reduction
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
    ActionChain := SumC([
        T0chain, T1chain, T2chain, Sc(-1, ThetaChain)
    ]);

    mu0 := RatDiv(eta0, frakA);
    mu1 := RatDiv(RatMul(frakA, eta1), frakB);
    mu2 := RatMul(frakB, eta2);
    m0 := RatDiv(eta0, Acoeff);
    m1 := RatDiv(RatMul(Acoeff, eta1), Bcoeff);
    m2 := RatMul(Bcoeff, eta2);
    theta := RatMul(PhiSub(z, e, e), PhiSub(z, e2, e));
    kappa := RatMul(Phi(e, e, e), Phi(e, e2, e));

    AssertRat("<Phi,E>", EvalChain(Ecoeff), aEpsilon);
    AssertRat("<Phi,A>", EvalChain(Achain), frakA);
    AssertRat("<Phi,B>", EvalChain(Bchain), frakB);
    AssertRat("<Phi,T_0>", EvalChain(T0chain), eta0);
    AssertRat("<Phi,T_1>", EvalChain(T1chain), eta1);
    AssertRat("<Phi,T_2>", EvalChain(T2chain), eta2);
    AssertRat("<Phi,M_0>", EvalChain(M0chain), mu0);
    AssertRat("<Phi,M_1>", EvalChain(M1chain), mu1);
    AssertRat("<Phi,M_2>", EvalChain(M2chain), mu2);
    AssertRat("<Phi,Theta>", EvalChain(ThetaChain), theta);
    AssertRat("<Phi,K>", EvalChain(Kchain), kappa);
    AssertRat("C_02 evaluation",
        EvalChain(C02chain), RatDiv(RatDiv(mu0, mu2), kappa));
    AssertRat("C_pow evaluation",
        EvalChain(CpowChain),
        RatMul(RatDiv(RatPow(mu0, 3), theta), kappa));
    AssertRat("three-action chain evaluation", EvalChain(ActionChain),
        RatDiv(RatProduct([eta0, eta1, eta2]), theta));

    AssertRat("m_0--mu_0 dictionary", m0,
        RatNeg(RatDiv(mu0, RatProduct([
            RatPow(rhoG, 2), rhoH, tauG
        ]))));
    AssertRat("m_1--mu_1 dictionary", m1,
        RatNeg(RatDiv(RatProduct([rhoG, rhoH, tauG, mu1]),
                         RatSub(RatOne, RatMul(rhoH, tauG)))));
    AssertRat("m_2--mu_2 dictionary", m2,
        RatProduct([
            rhoG, RatSub(RatOne, RatMul(rhoH, tauG)), mu2
        ]));
    AssertRat("product cancellation in the m--mu dictionary",
              RatProduct([m0, m1, m2]),
              RatProduct([mu0, mu1, mu2]));

    # The bar verifier proves that C_02 and C_pow are cycles whose S_3
    # projections are explicit boundaries.  Together with the injectivity
    # lemma proved in the manuscript, this makes them boundaries over
    # Gamma_3.  Their evaluations then give mu_2=mu_0/kappa and
    # theta=kappa*mu_0^3.  Together with the source projective-action law
    # mu_0*mu_1*mu_2=theta, compute mu_1 and compare it with mu_0/kappa.
    formalMu0 := Atom(["comparison", "mu_0"]);
    formalK := Atom(["comparison", "kappa"]);
    formalTheta := RatMul(RatPow(formalMu0, 3), formalK);
    derivedMu2 := RatDiv(formalMu0, formalK);
    derivedMu1 := RatDiv(formalTheta,
                         RatMul(formalMu0, derivedMu2));
    expectedMu1 := RatDiv(formalMu0, formalK);
    commonQuotient := RatDiv(derivedMu1, expectedMu1);
    AssertRat("common mu reduction before kappa^3=1",
              commonQuotient, RatPow(formalK, 3));

    extraFactor := Atom(["negative test", "extra c_12 factor"]);
    perturbedC12 := RatMul(c12Coefficient, extraFactor);
    perturbedLastPath := RatMul(perturbedC12, shiftedX1Scale);
    RejectRat("perturbed last recursive path", perturbedLastPath, Bcoeff);

    Print("[PASS] (3,1)_1 coefficients: all three recursive paths agree.\n");
    Print("[PASS] (3,1)_1 action: transported epsilon coefficients agree.\n");
    Print("[PASS] (3,1)_1 chain bridge: E/A/B/T/M/Theta/K agree.\n");
    Print("[PASS] (3,1)_1 dictionary: all three m_i--mu_i formulas agree.\n");
    Print("[PASS] (3,1)_1 comparisons: common mu quotient is kappa^3.\n");
    Print("[PASS] (3,1)_1 negative test: a perturbed c_12 path is rejected.\n");
    Print("       Atomic symbolic generators used: ", Length(symbolKeys), ".\n");
    return true;
end;

G3_GAMMA331_DIM1_OK := VerifyGamma331Dim1CoefficientAndBridge();

VerifyGamma331Dim1ScalarReduction := function()
    local ring, vars, p, t, k, RatDiv, RatEq,
          m0, m1, m2, r01, r02, r12,
          rhs01, rhs02, rhs12;

    ring := PolynomialRing(Rationals, 3);
    vars := IndeterminatesOfPolynomialRing(ring);
    p := vars[1];
    t := vars[2];
    k := vars[3];

    # Rational functions are represented as [numerator,denominator].
    RatDiv := function(x, y)
        return [x[1]*y[2], x[2]*y[1]];
    end;
    RatEq := function(x, y)
        return x[1]*y[2] = y[1]*x[2];
    end;

    # We cancel the common nonzero factor mu_0 and use
    # mu_1=mu_2=mu_0*kappa^{-1}.
    m0 := [-1, p^2*t];
    m1 := [-p*t, k*(1-t)];
    m2 := [p*(1-t), k];

    r01 := RatDiv(m0,m1);
    r02 := RatDiv(m0,m2);
    r12 := RatDiv(m1,m2);
    rhs01 := [k*(1-t), p^3*t^2];
    rhs02 := [-k, p^3*t*(1-t)];
    rhs12 := [-t, (1-t)^2];

    if not RatEq(r01,rhs01) then G3_Fail("(3,1)_1 m0/m1 reduction"); fi;
    if not RatEq(r02,rhs02) then G3_Fail("(3,1)_1 m0/m2 reduction"); fi;
    if not RatEq(r12,rhs12) then G3_Fail("(3,1)_1 m1/m2 reduction"); fi;

    # Exact polynomial consequences used after m0=m1=m2.
    if t*(1-t)-1 <> -(1-t+t^2) then
        G3_Fail("(3,1)_1 t(1-t)=1 reduction");
    fi;
    if t^3+1 <> (t+1)*(t^2-t+1) then
        G3_Fail("(3,1)_1 t^3=-1 reduction");
    fi;
    if (p*t)^3-1 <> (p*t-1)*((p*t)^2+p*t+1) then
        G3_Fail("(3,1)_1 pt=1 cubic reduction");
    fi;

    Print("[PASS] (3,1)_1: shared C02/Cpow bar certificates were rechecked.\n");
    Print("[PASS] (3,1)_1: all three m_i ratios and branch polynomials agree exactly.\n");
    return true;
end;

if not G3_GAMMA331_DIM1_OK then
    G3_Fail("The (3,1)_1 coefficient/chain verifier did not report success");
fi;
G3_GAMMA331_DIM1_OK := VerifyGamma331Dim1ScalarReduction();
