Read("bar_complex.g");

VerifyGamma332 := function()
    local one, e, e2, g, h, G3Mul, G3Inv, G3Conj,
          Bar, AddC, SubC, Sc, SumC, L, U, R,
          p, q, y0, y1, A, B0, S, E, K, C1, C2,
          T0, T1, T2, Cs, Clambda, Cr,
          SOne, S3Mul, Project, Proj, SBar,
          sg, se, se2, seg, se2g, Qs, Qlambda, Qr;

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
        else
            return [x[1] mod 3, -x[2], -x[3]];
        fi;
    end;
    G3Conj := function(a, x)
        return G3Mul(G3Mul(a, x), G3Inv(a));
    end;

    Bar  := function(xs) return BC_Bar(one, xs); end;
    AddC := function(x, y) return BC_Add(x, y, one); end;
    SubC := function(x, y) return BC_Sub(x, y, one); end;
    Sc   := function(n, x) return BC_Scale(n, x, one); end;
    SumC := function(xs) return BC_Sum(xs, one); end;

    L := function(x, a, b)
        return SumC([
            Bar([a, b, x]),
            Bar([G3Conj(G3Mul(a, b), x), a, b]),
            Sc(-1, Bar([a, G3Conj(b, x), b]))
        ]);
    end;
    U := function(a, x, y)
        return SumC([
            Bar([a, x, y]),
            Bar([G3Conj(a, x), G3Conj(a, y), a]),
            Sc(-1, Bar([G3Conj(a, x), a, y]))
        ]);
    end;
    R := function(a, x, y)
        return SubC(Bar([a, x, y]), Bar([G3Conj(a, x), a, y]));
    end;

    p  := [g, G3Mul(g, e), G3Mul(g, e2)];
    q  := List(p, x -> G3Mul(x, h));
    y0 := G3Mul(e, h);
    y1 := G3Mul(e2, h);

    A := SumC([
        L(g, G3Mul(e2, h), e2),
        Sc(-1, L(g, e, h)),
        Sc(-1, L(y0, g, e2))
    ]);
    B0 := SumC([
        U(g, p[3], y0), L(g, g, e2), Sc(-1, L(g, e, g))
    ]);
    S := SumC([
        B0, Sc(-1, A), Sc(-1, L(g, g, h)), L(y0, e, e)
    ]);
    E := AddC(L(y0, e, e), L(y0, e2, e));
    K := AddC(Bar([e, e, e]), Bar([e, e2, e]));

    C1 := SumC([
        R(g, p[3], y0), L(g, g, e2), Sc(-1, L(g, e, g)),
        Sc(-1, U(e, p[3], y0)), Sc(-1, L(g, e, e2))
    ]);
    C2 := SumC([
        R(g, p[2], y1), L(g, g, e), Sc(-1, L(g, e2, g)),
        Sc(-1, U(e2, p[2], y1)), Sc(-1, L(g, e2, e)),
        Sc(-1, L(y0, e2, g)), L(y0, g, e)
    ]);

    T0 := U(e, p[1], q[1]);
    T1 := SumC([U(e, p[2], q[2]), L(g, e, e),
                L(G3Mul(g, h), e, e)]);
    T2 := SumC([U(e, p[3], q[3]), L(g, e, e2),
                L(G3Mul(g, h), e, e2)]);

    Cs      := SumC([Sc(3, S), Sc(-2, E), Sc(-1, K)]);
    Clambda := SumC([Sc(3, AddC(AddC(C1, T1), Sc(-1, C2))),
                     Sc(-1, T0), Sc(-1, T1), Sc(-1, T2)]);
    Cr      := SumC([Sc(4, C1), Sc(2, T1), Sc(-2, C2),
                     Sc(-2, T0), Sc(-1, S)]);

    BC_AssertZero("Gamma32 Cs differential", BC_Boundary(Cs, G3Mul, one), one);
    BC_AssertZero("Gamma32 Clambda differential", BC_Boundary(Clambda, G3Mul, one), one);
    BC_AssertZero("Gamma32 Cr differential", BC_Boundary(Cr, G3Mul, one), one);
    if Length(Cs) <> 25 or Length(Clambda) <> 44 or Length(Cr) <> 53 then
        G3_Fail("Gamma32 chain term counts differ from the recorded certificate");
    fi;

    SOne := [0, 0, 0];
    S3Mul := function(x, y)
        local sign;
        if x[2] mod 2 = 0 then sign := 1; else sign := -1; fi;
        return [(x[1] + sign * y[1]) mod 3,
                (x[2] + y[2]) mod 2, 0];
    end;
    Project := function(x) return [x[1], x[2] mod 2, 0]; end;
    Proj := function(chain)
        return BC_Project(chain, Project, one, SOne);
    end;
    SBar := function(xs) return BC_Bar(SOne, xs); end;

    sg   := [0, 1, 0];
    se   := [1, 0, 0];
    se2  := [2, 0, 0];
    seg  := [1, 1, 0];
    se2g := [2, 1, 0];

    Qs := BC_Sum([
        BC_Scale(-3, SBar([seg,se,se2,se]), SOne),
        BC_Scale( 3, SBar([se2,seg,se,se2]), SOne),
        BC_Scale(-3, SBar([se2g,se2,seg,se]), SOne),
        BC_Scale( 3, SBar([se2g,se2,seg,se2]), SOne),
        BC_Scale(-3, SBar([se2g,se2,se2,sg]), SOne)
    ], SOne);

    Qlambda := BC_Sum([
        BC_Scale(-2,SBar([sg,sg,se,se2]),SOne),
        BC_Scale(-1,SBar([sg,se,sg,se2]),SOne),
        BC_Scale( 2,SBar([sg,se,se,sg]),SOne),
        BC_Scale(-1,SBar([sg,se,se2,sg]),SOne),
        BC_Scale(-1,SBar([sg,se,se2,se2g]),SOne),
        BC_Scale(-1,SBar([sg,se,se2g,se]),SOne),
        BC_Scale( 2,SBar([sg,seg,se,se]),SOne),
        BC_Scale(-1,SBar([sg,se2,se2,se2]),SOne),
        BC_Scale( 1,SBar([sg,se2g,se,se]),SOne),
        BC_Scale( 1,SBar([se,sg,se2,se2]),SOne),
        BC_Scale( 1,SBar([se,seg,se2,se2g]),SOne),
        BC_Scale( 2,SBar([se,se2,se,sg]),SOne),
        BC_Scale(-2,SBar([se,se2,se2g,se]),SOne),
        BC_Scale( 3,SBar([se,se2,se2g,se2]),SOne),
        BC_Scale(-1,SBar([se,se2g,se,se]),SOne),
        BC_Scale(-1,SBar([se,se2g,se2,sg]),SOne),
        BC_Scale(-1,SBar([seg,sg,se,se]),SOne),
        BC_Scale( 2,SBar([seg,sg,se2,se]),SOne),
        BC_Scale( 2,SBar([seg,se,se,sg]),SOne),
        BC_Scale(-2,SBar([seg,se,se2g,se]),SOne),
        BC_Scale( 1,SBar([seg,se2,sg,se]),SOne),
        BC_Scale(-1,SBar([se2,se2,sg,se]),SOne),
        BC_Scale( 1,SBar([se2,se2,se2,sg]),SOne),
        BC_Scale(-2,SBar([se2g,sg,se,se]),SOne)
    ], SOne);

    Qr := BC_Sum([
        BC_Scale(-2,SBar([sg,sg,se,se2]),SOne),
        BC_Scale(-1,SBar([sg,se,sg,se2]),SOne),
        BC_Scale( 1,SBar([sg,se,se,sg]),SOne),
        BC_Scale(-1,SBar([sg,se,se2,seg]),SOne),
        BC_Scale( 2,SBar([sg,seg,se,se]),SOne),
        BC_Scale( 1,SBar([se,se,sg,se2]),SOne),
        BC_Scale(-1,SBar([se,se,se,sg]),SOne),
        BC_Scale(-1,SBar([se,se,seg,se]),SOne),
        BC_Scale( 1,SBar([se,se,seg,se2]),SOne),
        BC_Scale(-1,SBar([se,se,se2,sg]),SOne),
        BC_Scale( 1,SBar([se,se2,se,sg]),SOne),
        BC_Scale(-2,SBar([se,se2,se2g,se]),SOne),
        BC_Scale( 2,SBar([se,se2,se2g,se2]),SOne),
        BC_Scale(-1,SBar([se,se2g,se,se]),SOne),
        BC_Scale(-2,SBar([se,se2g,se2,sg]),SOne),
        BC_Scale( 1,SBar([seg,sg,se2,se]),SOne),
        BC_Scale( 2,SBar([seg,se,se,sg]),SOne),
        BC_Scale(-2,SBar([seg,se,se2g,se]),SOne),
        BC_Scale( 1,SBar([seg,seg,se,se2]),SOne),
        BC_Scale(-2,SBar([se2g,sg,se,se]),SOne),
        BC_Scale( 1,SBar([se2g,se,se,se2]),SOne),
        BC_Scale(-1,SBar([se2g,se2,seg,se]),SOne)
    ], SOne);

    BC_AssertEqual("Gamma32 projected Cs boundary",
        BC_Boundary(Qs, S3Mul, SOne), Proj(Cs), SOne);
    BC_AssertEqual("Gamma32 projected Clambda boundary",
        BC_Boundary(Qlambda, S3Mul, SOne), Proj(Clambda), SOne);
    BC_AssertEqual("Gamma32 projected Cr boundary",
        BC_Boundary(Qr, S3Mul, SOne), Proj(Cr), SOne);

    Print("[PASS] (3,2): Cs, Clambda, Cr have zero Gamma_3 differential.\n");
    Print("[PASS] (3,2): explicit S_3 four-chain residuals are zero.\n");
    Print("       Gamma terms: 25, 44, 53; S3 four-chain terms: 5, 24, 22.\n");
    return true;
end;

G3_GAMMA332_OK := VerifyGamma332();

# The preceding certificate starts with the bar chains printed in Appendix E.
# The following independent layer starts one step earlier.  It expands the
# local cocycles into values of the normalized 3-cocycle, reconstructs the
# relevant induced actions and braided-adjoint coefficients, and only then
# compares the resulting Laurent monomials with the bar-chain evaluations.

VerifyGamma332CoefficientAndDictionary := function()
    local one, e, e2, g, h, G3Mul, G3Inv, G3Conj,
          Bar, AddC, SubC, Sc, SumC, L, U, R,
          NormalizeAtoms, Monomial, OneM, NegOneM, AtomM, PhiM,
          MulM, InvM, PowM, NegM, ChainM, SameUnsignedM,
          AssertM, ReplaceAtom, RawLocal, RawTensor, RawR,
          NormalizeActionM, AssertCocycleM, RejectPerturbation,
          RhoM, SigmaM, TauM, ei, VAction, WAction, XAction,
          ScaleVec, AssertVec, Phi1VW, TensorActionVW,
          p, q, y0, y1, q0, g2h, eg2h,
          rhoG, rhoH, rhoGH, sigmaE, sigmaE2, sigmaG2,
          aRaw, b0Raw, b1Raw, chiRaw, ell0Raw, ell1Raw,
          k1Raw, k2Raw, c1Raw, c2Raw, tRaw,
          av, aw, ax, coeff, doubleCoeffE1, u0Vec, u2Vec,
          phi1E1Vec, phi1v0w0Vec, phi1v0w1Vec,
          A, B0, S, E, K, C1, C2, T0, T1, T2,
          Cs, Clambda, Cr, Eraw, Normalize32, Assert32,
          sRaw, lambdaRaw, rRaw, kappaRaw, expected,
          i, centralC1, centralC2, centralSource, centralTarget,
          sCentral, noncentralThird,
          Dphi10, DphiE1, Dphi01;

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
        else
            return [x[1] mod 3, -x[2], -x[3]];
        fi;
    end;
    G3Conj := function(a, x)
        return G3Mul(G3Mul(a, x), G3Inv(a));
    end;

    Bar  := function(xs) return BC_Bar(one, xs); end;
    AddC := function(x, y) return BC_Add(x, y, one); end;
    SubC := function(x, y) return BC_Sub(x, y, one); end;
    Sc   := function(n, x) return BC_Scale(n, x, one); end;
    SumC := function(xs) return BC_Sum(xs, one); end;
    L := function(x, a, b)
        return SumC([
            Bar([a, b, x]),
            Bar([G3Conj(G3Mul(a, b), x), a, b]),
            Sc(-1, Bar([a, G3Conj(b, x), b]))
        ]);
    end;
    U := function(a, x, y)
        return SumC([
            Bar([a, x, y]),
            Bar([G3Conj(a, x), G3Conj(a, y), a]),
            Sc(-1, Bar([G3Conj(a, x), a, y]))
        ]);
    end;
    R := function(a, x, y)
        return SubC(Bar([a, x, y]), Bar([G3Conj(a, x), a, y]));
    end;

    # A Laurent monomial consists of a sign, an exponent vector of raw
    # Phi-values (stored as a normalized bar three-chain), and exponents of
    # representation scalars.  There is no substitution by a target formula.
    NormalizeAtoms := function(atoms)
        local out, term, pos;
        out := [];
        for term in atoms do
            if term[2] <> 0 then
                pos := PositionProperty(out, x -> x[1] = term[1]);
                if pos = fail then
                    Add(out, [term[1], term[2]]);
                else
                    out[pos][2] := out[pos][2] + term[2];
                fi;
            fi;
        od;
        out := Filtered(out, x -> x[2] <> 0);
        Sort(out, function(x, y) return x[1] < y[1]; end);
        return out;
    end;
    Monomial := function(sign, phi, atoms)
        if sign <> 1 and sign <> -1 then
            G3_Fail("Gamma32 Laurent monomial has a non-sign coefficient");
        fi;
        return rec(sign := sign,
                   phi := BC_Normalize(phi, one),
                   atoms := NormalizeAtoms(atoms));
    end;
    OneM := Monomial(1, [], []);
    NegOneM := Monomial(-1, [], []);
    AtomM := function(name) return Monomial(1, [], [[name, 1]]); end;
    PhiM := function(a, b, c) return Monomial(1, Bar([a, b, c]), []); end;
    MulM := function(x, y)
        return Monomial(x.sign * y.sign,
                        AddC(x.phi, y.phi),
                        Concatenation(x.atoms, y.atoms));
    end;
    InvM := function(x)
        return Monomial(x.sign, Sc(-1, x.phi),
                        List(x.atoms, a -> [a[1], -a[2]]));
    end;
    PowM := function(x, n)
        if n = 0 then return OneM; fi;
        return Monomial(x.sign^n, Sc(n, x.phi),
                        List(x.atoms, a -> [a[1], n * a[2]]));
    end;
    NegM := function(x) return MulM(NegOneM, x); end;
    ChainM := function(chain) return Monomial(1, chain, []); end;
    SameUnsignedM := function(x, y)
        return x.phi = y.phi and x.atoms = y.atoms;
    end;
    AssertM := function(label, left, right)
        local residual;
        residual := MulM(left, InvM(right));
        if residual.sign <> 1 or Length(residual.phi) <> 0 or
           Length(residual.atoms) <> 0 then
            G3_Fail(label, ": unequal Laurent monomials; residual = ", residual);
        fi;
    end;
    ReplaceAtom := function(m, name, replacement)
        local pos, exponent, remaining;
        pos := PositionProperty(m.atoms, x -> x[1] = name);
        if pos = fail then return m; fi;
        exponent := m.atoms[pos][2];
        remaining := Filtered(m.atoms, x -> x[1] <> name);
        return MulM(Monomial(m.sign, m.phi, remaining),
                    PowM(replacement, exponent));
    end;

    # These are literal expansions of Phi_x, Phi^a, and R(a;x,y).
    RawLocal := function(x, a, b)
        return MulM(MulM(PhiM(a, b, x),
                         PhiM(G3Conj(G3Mul(a, b), x), a, b)),
                    InvM(PhiM(a, G3Conj(b, x), b)));
    end;
    RawTensor := function(a, x, y)
        return MulM(MulM(PhiM(a, x, y),
                         PhiM(G3Conj(a, x), G3Conj(a, y), a)),
                    InvM(PhiM(G3Conj(a, x), a, y)));
    end;
    RawR := function(a, x, y)
        return MulM(PhiM(a, x, y),
                    InvM(PhiM(G3Conj(a, x), a, y)));
    end;

    rhoG := AtomM("rho(g)");
    rhoH := AtomM("rho(h)");
    rhoGH := AtomM("rho(gh)");
    sigmaE := AtomM("sigma(epsilon)");
    sigmaE2 := AtomM("sigma(epsilon^2)");
    sigmaG2 := AtomM("sigma(g^2)");

    RhoM := function(c)
        if c = one then return OneM;
        elif c = g then return rhoG;
        elif c = h then return rhoH;
        elif c = G3Mul(g, h) then return rhoGH;
        else G3_Fail("Gamma32 unexpected rho argument: ", c);
        fi;
    end;
    SigmaM := function(c)
        if c = one then return OneM;
        elif c = e then return sigmaE;
        elif c = e2 then return sigmaE2;
        elif c = G3Mul(g, g) then return sigmaG2;
        elif c = G3Mul(e2, G3Mul(g, g)) then
            # sigma(epsilon^2 g^2), obtained from
            # eq:YD-projective-action applied to
            # sigma(epsilon^2)sigma(g^2).
            return MulM(MulM(sigmaE2, sigmaG2),
                        InvM(RawLocal(y0, e2, G3Mul(g, g))));
        else G3_Fail("Gamma32 unexpected sigma argument: ", c);
        fi;
    end;

    p := [g, G3Mul(g, e), G3Mul(g, e2)];
    q := List(p, x -> G3Mul(x, h));
    y0 := G3Mul(e, h);
    y1 := G3Mul(e2, h);
    q0 := G3Mul(g, h);
    g2h := G3Mul(G3Mul(g, g), h);
    eg2h := G3Mul(e, g2h);
    ei := function(i) return [i mod 3, 0, 0]; end;

    # Canonical consequences of two applications of eq:YD-projective-action:
    # sigma(epsilon)^2 and sigma(epsilon^2)sigma(epsilon).  This
    # normalization is used for direct vector comparisons, without the
    # later assumption rho(g)=-1.
    Eraw := MulM(RawLocal(y0, e, e), RawLocal(y0, e2, e));
    NormalizeActionM := function(m)
        local out, expansion, pos, exponent, remainder, quotient;
        expansion := MulM(PowM(sigmaE, 2),
                          InvM(RawLocal(y0, e, e)));
        out := ReplaceAtom(m, "sigma(epsilon^2)", expansion);
        pos := PositionProperty(out.atoms,
                                x -> x[1] = "sigma(epsilon)");
        if pos <> fail then
            exponent := out.atoms[pos][2];
            remainder := exponent mod 3;
            quotient := (exponent - remainder) / 3;
            out := ReplaceAtom(out, "sigma(epsilon)", OneM);
            out := MulM(out, PowM(sigmaE, remainder));
            out := MulM(out, PowM(Eraw, quotient));
        fi;
        return out;
    end;
    AssertCocycleM := function(label, left, right, fourChain)
        local residual;
        residual := MulM(NormalizeActionM(left),
                         InvM(NormalizeActionM(right)));
        if residual.sign <> 1 or Length(residual.atoms) <> 0 then
            G3_Fail(label, ": non-cocycle residual = ", residual);
        fi;
        BC_AssertZero(Concatenation(label, " cycle check"),
                      BC_Boundary(residual.phi, G3Mul, one), one);
        BC_AssertEqual(label, residual.phi,
                       BC_Boundary(fourChain, G3Mul, one), one);
    end;
    RejectPerturbation := function(label, left, right, normalize)
        local altered, residual, boundary;
        altered := MulM(left, PhiM(e, e, e));
        residual := MulM(normalize(altered), InvM(normalize(right)));
        if residual.sign <> 1 or Length(residual.atoms) <> 0 then
            G3_Fail(label, ": malformed negative-control residual");
        fi;
        boundary := BC_Boundary(residual.phi, G3Mul, one);
        if Length(boundary) = 0 then
            G3_Fail(label, ": one-factor perturbation was not detected");
        fi;
    end;

    # Induced action on v_i = epsilon^i action v_0.  The returned first
    # component is the new basis index, and the second is its coefficient.
    VAction := function(x, i)
        local product, j, c;
        product := G3Mul(x, ei(i));
        j := product[1] mod 3;
        c := G3Mul(G3Inv(ei(j)), product);
        return rec(index := j,
                   coeff := MulM(MulM(RawLocal(g, x, ei(i)), RhoM(c)),
                                 InvM(RawLocal(g, ei(j), c))));
    end;
    # Induced action on (w_0,w_1), with transversals (1,g).
    WAction := function(x, i)
        local ti, product, j, tj, c;
        if i = 0 then ti := one; else ti := g; fi;
        product := G3Mul(x, ti);
        j := product[2] mod 2;
        if j = 0 then tj := one; else tj := g; fi;
        c := G3Mul(G3Inv(tj), product);
        return rec(index := j,
                   coeff := MulM(MulM(RawLocal(y0, x, ti), SigmaM(c)),
                                 InvM(RawLocal(y0, tj, c))));
    end;

    # The character of the degree-q_0 line of X_1 is only needed at g.
    # chiRaw is assigned below, before XAction is called.
    TauM := function(c)
        if c = one then return OneM;
        elif c = g then return chiRaw;
        else G3_Fail("Gamma32 unexpected X_1 character argument: ", c);
        fi;
    end;
    XAction := function(x, i)
        local product, j, c;
        product := G3Mul(x, ei(i));
        j := product[1] mod 3;
        c := G3Mul(G3Inv(ei(j)), product);
        return rec(index := j,
                   coeff := MulM(MulM(RawLocal(q0, x, ei(i)), TauM(c)),
                                 InvM(RawLocal(q0, ei(j), c))));
    end;

    ScaleVec := function(scalar, vector)
        return List(vector, x -> [x[1], MulM(scalar, x[2])]);
    end;
    AssertVec := function(label, left, right)
        local residual, collected, term, normalized, unsigned, pos;
        residual := Concatenation(left, ScaleVec(NegOneM, right));
        collected := [];
        for term in residual do
            normalized := NormalizeActionM(term[2]);
            unsigned := Monomial(1, normalized.phi, normalized.atoms);
            pos := PositionProperty(collected,
                x -> x[1] = term[1] and SameUnsignedM(x[2], unsigned));
            if pos = fail then
                Add(collected, [term[1], unsigned, normalized.sign]);
            else
                collected[pos][3] := collected[pos][3] + normalized.sign;
            fi;
        od;
        collected := Filtered(collected, x -> x[3] <> 0);
        if Length(collected) <> 0 then
            G3_Fail(label, ": nonzero formal vector residual = ", collected);
        fi;
    end;
    Phi1VW := function(i, j)
        local actionW, actionV;
        actionW := WAction(p[i+1], j);
        actionV := VAction([y0, y1][actionW.index+1], i);
        return [
            [[i, j], OneM],
            [[actionV.index, actionW.index],
             NegM(MulM(actionW.coeff, actionV.coeff))]
        ];
    end;
    TensorActionVW := function(x, vector)
        local answer, term, i0, j0, actionV, actionW, factor;
        answer := [];
        for term in vector do
            i0 := term[1][1];
            j0 := term[1][2];
            actionV := VAction(x, i0);
            actionW := WAction(x, j0);
            factor := MulM(MulM(RawTensor(x, p[i0+1], [y0,y1][j0+1]),
                                 actionV.coeff), actionW.coeff);
            Add(answer, [[actionV.index, actionW.index],
                         MulM(term[2], factor)]);
        od;
        return answer;
    end;

    # First adjoint and induced-action coefficients, reconstructed from the
    # two braidings, eq:YD-projective-action, and eq:tensor-product.
    aw := WAction(p[3], 0);
    av := VAction(y1, 2);
    aRaw := MulM(aw.coeff, av.coeff);
    expected := MulM(MulM(RawLocal(g, y1, e2), rhoH), sigmaE2);
    expected := MulM(expected,
        InvM(MulM(RawLocal(g, e, h), RawLocal(y0, g, e2))));
    AssertM("Gamma32 coefficient a from the two braidings", aRaw, expected);

    av := VAction(g, 2);
    aw := WAction(g, 0);
    b0Raw := MulM(MulM(RawTensor(g, p[3], y0), av.coeff), aw.coeff);
    expected := MulM(MulM(RawTensor(g, p[3], y0),
                           MulM(RawLocal(g, g, e2), rhoG)),
                     InvM(RawLocal(g, e, g)));
    AssertM("Gamma32 coefficient b0 from g action E0", b0Raw, expected);
    AssertVec("Gamma32 g action E0=b0 E1",
              TensorActionVW(g, [[[2,0], OneM]]),
              [[[1,1], b0Raw]]);

    av := VAction(g, 1);
    aw := WAction(g, 1);
    b1Raw := MulM(MulM(RawTensor(g, p[2], y1), av.coeff), aw.coeff);
    expected := MulM(MulM(MulM(RawTensor(g, p[2], y1),
                                MulM(RawLocal(g, g, e), rhoG)),
                          InvM(RawLocal(g, e2, g))),
                     MulM(RawLocal(y0, g, g), sigmaG2));
    AssertM("Gamma32 coefficient b1 from g action E1", b1Raw, expected);
    AssertVec("Gamma32 g action E1=b1 E0",
              TensorActionVW(g, [[[1,1], OneM]]),
              [[[2,0], b1Raw]]);

    chiRaw := NegM(MulM(b0Raw, InvM(aRaw)));

    av := VAction(e, 2);
    aw := WAction(e, 0);
    ell0Raw := MulM(MulM(RawTensor(e, p[3], y0), av.coeff), aw.coeff);
    expected := MulM(MulM(RawTensor(e, p[3], y0),
                           RawLocal(g, e, e2)), sigmaE);
    AssertM("Gamma32 ell0 from epsilon action E0", ell0Raw, expected);
    AssertVec("Gamma32 epsilon action E0=ell0(v0 tensor w0)",
              TensorActionVW(e, [[[2,0], OneM]]),
              [[[0,0], ell0Raw]]);

    av := VAction(e2, 1);
    aw := WAction(e2, 1);
    ell1Raw := MulM(MulM(RawTensor(e2, p[2], y1), av.coeff), aw.coeff);
    expected := MulM(MulM(MulM(RawTensor(e2, p[2], y1),
                                RawLocal(g, e2, e)),
                          MulM(RawLocal(y0, e2, g), sigmaE)),
                     InvM(RawLocal(y0, g, e)));
    AssertM("Gamma32 ell1 from epsilon^2 action E1", ell1Raw, expected);
    AssertVec("Gamma32 epsilon^2 action E1=ell1(v0 tensor w1)",
              TensorActionVW(e2, [[[1,1], OneM]]),
              [[[0,1], ell1Raw]]);

    k2Raw := VAction(g, 2).coeff;
    expected := MulM(MulM(RawLocal(g, g, e2), rhoG),
                     InvM(RawLocal(g, e, g)));
    AssertM("Gamma32 k2 from g action v2", k2Raw, expected);
    k1Raw := VAction(g, 1).coeff;
    expected := MulM(MulM(RawLocal(g, g, e), rhoG),
                     InvM(RawLocal(g, e2, g)));
    AssertM("Gamma32 k1 from g action v1", k1Raw, expected);

    u0Vec := Phi1VW(2, 0);
    AssertVec("Gamma32 u0 = E0-a E1", u0Vec,
              [[[2,0], OneM], [[1,1], NegM(aRaw)]]);
    phi1v0w0Vec := Phi1VW(0, 0);
    aw := WAction(g, 0);
    av := VAction(y1, 0);
    coeff := MulM(aw.coeff, av.coeff);
    aw := WAction(e, 1);
    av := VAction(e, 1);
    expected := MulM(MulM(RawTensor(e, p[2], y1), aw.coeff), av.coeff);
    AssertVec("Gamma32 direct phi1(v0 tensor w0) support",
              phi1v0w0Vec,
              [[[0,0], OneM], [[2,1], NegM(coeff)]]);
    Dphi10 := SumC([
        Sc(-1,Bar([g,e,e,h])),
        Sc(-1,Bar([G3Mul(e,g),e,G3Mul(e,h),e])),
        Bar([e,G3Mul(e2,h),g,e]),
        Bar([g,e,e,G3Mul(e,h)]),
        Bar([G3Mul(e2,h),g,e2,e2]),
        Bar([e,G3Mul(e2,g),e,G3Mul(e,h)]),
        Sc(-1,Bar([e,e,g,h])),
        Bar([G3Mul(e,h),e,e,e2]),
        Bar([g,e2,G3Mul(e,h),e]),
        Bar([G3Mul(e2,h),e,e2,g]),
        Sc(-1,Bar([g,G3Mul(e,h),e2,e2])),
        Sc(-1,Bar([G3Mul(e2,h),e,g,e])),
        Sc(-1,Bar([e,G3Mul(e2,h),e2,g])),
        Bar([g,e,G3Mul(e2,h),e2]),
        Sc(-1,Bar([G3Mul(e,g),e,e,e2])),
        Bar([e,e,h,g]),
        Bar([e,g,e2,h]),
        Sc(-1,Bar([e,g,e2,e2])),
        Sc(-1,Bar([e,g,G3Mul(e,h),e])),
        Bar([G3Mul(e2,h),g,e2,e]),
        Sc(-1,Bar([g,G3Mul(e,h),e2,e])),
        Sc(-1,Bar([g,e2,e,G3Mul(e,h)]))
    ]);
    AssertCocycleM("Gamma32 phi1(v0 tensor w0) second coordinate",
                    coeff,
                    MulM(MulM(aRaw, expected), InvM(ell0Raw)),
                    Dphi10);

    # Compute phi1(E1) directly, rather than inserting its value.  The
    # coefficient of E0 in its double-braiding term is then compared with
    # a*b1/b0.  Under the single stated assumption
    # Delta_1=a^2*b1/b0=1, this is a^{-1}, so the checked vector is
    # -a^{-1}u0.
    phi1E1Vec := Phi1VW(1, 1);
    aw := WAction(p[2], 1);
    av := VAction(y0, 1);
    doubleCoeffE1 := MulM(aw.coeff, av.coeff);
    AssertVec("Gamma32 direct phi1(E1) expansion", phi1E1Vec,
              [[[1,1], OneM], [[2,0], NegM(doubleCoeffE1)]]);
    DphiE1 := SumC([
        Bar([e2,g,y1,g]),
        Bar([p[2],e2,h,g]),
        Sc(-1,Bar([e2,e2,h,g])),
        Bar([p[2],e2,p[3],e2]),
        Sc(-1,Bar([g,e2,g,y1])),
        Bar([y1,g,e2,e2]),
        Bar([g,e2,g,e]),
        Bar([y0,e2,g,g]),
        Sc(-1,Bar([p[2],e2,h,e2])),
        Bar([e2,e2,g,h]),
        Sc(-1,Bar([p[2],e2,e2,g])),
        Sc(-1,Bar([g,e2,e2,g])),
        Bar([p[2],e2,e2,h]),
        Sc(-1,Bar([e2,p[3],e2,h])),
        Bar([e2,h,e2,g]),
        Sc(-1,Bar([p[2],e2,p[3],y0])),
        Sc(-1,Bar([e2,g,g,y0])),
        Sc(-1,Bar([g,y0,e2,e2])),
        Sc(-1,Bar([e2,h,g,e])),
        Bar([g,e2,y0,g]),
        Bar([e2,p[3],y0,e]),
        Sc(-1,Bar([e2,y0,g,g])),
        Bar([y0,e2,e2,g]),
        Sc(-1,Bar([y0,e2,g,e]))
    ]);
    AssertCocycleM("Gamma32 double-braiding coefficient on E1",
                    doubleCoeffE1,
                    MulM(MulM(aRaw, b1Raw), InvM(b0Raw)),
                    DphiE1);

    u2Vec := TensorActionVW(e2, u0Vec);
    phi1v0w1Vec := Phi1VW(0, 1);
    av := VAction(e2, 2);
    aw := WAction(e2, 0);
    coeff := MulM(MulM(RawTensor(e2, p[3], y0), av.coeff), aw.coeff);
    aw := WAction(g, 1);
    av := VAction(y0, 0);
    expected := MulM(aw.coeff, av.coeff);
    AssertVec("Gamma32 direct phi1(v0 tensor w1) support",
              phi1v0w1Vec,
              [[[0,1], OneM], [[1,0], NegM(expected)]]);
    AssertVec("Gamma32 direct u2 expansion", u2Vec,
              [[[1,0], coeff],
               [[0,1], NegM(MulM(aRaw, ell1Raw))]]);
    Dphi01 := SumC([
        Bar([e2,g,y1,g]),
        Sc(-1,Bar([e2,e,e,y0])),
        Bar([y0,e2,e,e]),
        Bar([e2,p[2],e,y1]),
        Bar([y0,e2,g,g]),
        Bar([e2,e,h,g]),
        Sc(-1,Bar([g,e2,h,e2])),
        Bar([p[2],e,h,e2]),
        Bar([e2,h,e2,g]),
        Sc(-1,Bar([e2,g,g,y0])),
        Sc(-1,Bar([g,e2,e,y1])),
        Bar([p[2],e2,e2,y0]),
        Sc(-1,Bar([e2,p[3],e2,y0])),
        Sc(-1,Bar([e2,h,g,e])),
        Bar([e,h,e,e]),
        Bar([e2,p[2],e2,h]),
        Bar([e2,p[3],y0,e]),
        Sc(-1,Bar([e2,y0,g,g])),
        Bar([e2,g,e,y0]),
        Sc(-1,Bar([e2,e,g,h])),
        Sc(-1,Bar([e2,p[2],e2,y0])),
        Bar([y0,e2,e2,g]),
        Sc(-1,Bar([y0,e2,g,e])),
        Sc(-1,Bar([e2,e,h,e]))
    ]);
    AssertCocycleM("Gamma32 epsilon^2 equivariance coordinate",
                    MulM(doubleCoeffE1, coeff),
                    MulM(ell1Raw, expected), Dphi01);
    # Substituting Delta_1=1 in the two checked vector identities gives
    # phi1(v0 tensor w1)=-(a ell1)^{-1}u2.

    # Degrees and the complete central coefficient calculation for varphi_2.
    for i in [1..3] do
        if G3Mul(p[i], q[i]) <> g2h then
            G3_Fail("Gamma32 P_i degree mismatch at i = ", i-1);
        fi;
    od;
    if G3Mul(p[3], q[1]) <> eg2h or G3Mul(p[1], q[2]) <> eg2h then
        G3_Fail("Gamma32 noncentral tensor degrees do not agree");
    fi;

    ax := XAction(g, 0);
    av := VAction(q0, 0);
    if ax.index <> 0 or av.index <> 0 then
        G3_Fail("Gamma32 central double braiding has the wrong target");
    fi;
    AssertM("Gamma32 g action on u0", ax.coeff, chiRaw);
    AssertM("Gamma32 gh action on v0", av.coeff, rhoGH);
    AssertM("Gamma32 central double-braiding correction",
            NegM(MulM(ax.coeff, av.coeff)),
            MulM(MulM(b0Raw, rhoGH), InvM(aRaw)));

    # Reconstruct the complete third recursive summand from the two terms
    # E0-aE1 of u0.  The source coefficients use the induced actions and
    # the three phi1 identities checked above.  The targets use the printed
    # definitions of c1 and c2.
    av := VAction(g, 2);
    if av.index <> 1 then
        G3_Fail("Gamma32 E0 route in c12 has the wrong outer V index");
    fi;
    centralC1 := MulM(MulM(RawR(g, p[3], y0), av.coeff),
                      InvM(ell0Raw));
    av := VAction(g, 1);
    if av.index <> 2 then
        G3_Fail("Gamma32 E1 route in c12 has the wrong outer V index");
    fi;
    centralC2 := MulM(MulM(MulM(NegM(aRaw), RawR(g, p[2], y1)),
                            av.coeff),
                      NegM(InvM(MulM(aRaw, ell1Raw))));
    c1Raw := MulM(
        MulM(RawR(g, p[3], y0),
             MulM(MulM(RawLocal(g, g, e2), rhoG),
                  InvM(RawLocal(g, e, g)))),
        InvM(MulM(MulM(RawTensor(e, p[3], y0),
                         RawLocal(g, e, e2)), sigmaE)));
    c2Raw := MulM(
        MulM(RawR(g, p[2], y1),
             MulM(MulM(RawLocal(g, g, e), rhoG),
                  InvM(RawLocal(g, e2, g)))),
        InvM(MulM(MulM(MulM(RawTensor(e2, p[2], y1),
                              RawLocal(g, e2, e)),
                        MulM(RawLocal(y0, e2, g), sigmaE)),
                  InvM(RawLocal(y0, g, e)))));
    centralSource := [[1, centralC1], [2, centralC2]];
    centralTarget := [[1, c1Raw], [2, c2Raw]];
    AssertVec("Gamma32 complete third varphi2 summand on P0",
              centralSource, centralTarget);
    RejectPerturbation("Gamma32 coefficient negative control",
                       centralC1, c1Raw, NormalizeActionM);

    sCentral := MulM(MulM(b0Raw, rhoGH), InvM(aRaw));
    centralSource := Concatenation(
        [[0, OneM], [0, NegM(MulM(chiRaw, rhoGH))]], centralSource);
    centralTarget := Concatenation(
        [[0, OneM], [0, sCentral]], centralTarget);
    AssertVec("Gamma32 complete central varphi2 vector",
              centralSource, centralTarget);

    if G3Conj(p[3], q[1]) <> q[2] or G3Conj(q[2], p[3]) <> p[1] or
       G3Conj(p[3], p[2]) <> p[1] or G3Mul(p[3], y1) <> q[2] then
        G3_Fail("Gamma32 noncentral conjugation targets are incorrect");
    fi;
    AssertM("Gamma32 R(p2;p2,y0)=1", RawR(p[3], p[3], y0), OneM);
    AssertM("Gamma32 p2 action on v2", VAction(p[3], 2).coeff, rhoG);
    noncentralThird := MulM(RawR(p[3], p[3], y0),
                            VAction(p[3], 2).coeff);
    AssertM("Gamma32 noncentral input coordinate in third summand",
            noncentralThird, rhoG);
    ax := XAction(p[3], 0);
    av := VAction(q[2], 2);
    if ax.index <> 1 or av.index <> 0 or VAction(p[3], 1).index <> 0 then
        G3_Fail("Gamma32 noncentral off-coordinate target is incorrect");
    fi;
    # The double braiding and the E1 route therefore lie in
    # C(v0 tensor u1); only id and the checked E0 route contribute to the
    # v2 tensor u0 coordinate, with coefficient 1+rho(g).

    tRaw := [];
    for i in [0..2] do
        av := VAction(e, i);
        ax := XAction(e, i);
        if av.index <> (i+1) mod 3 or ax.index <> (i+1) mod 3 then
            G3_Fail("Gamma32 epsilon action has the wrong P_i target");
        fi;
        coeff := MulM(MulM(RawTensor(e, p[i+1], q[i+1]), av.coeff),
                      ax.coeff);
        Add(tRaw, coeff);
    od;
    AssertM("Gamma32 t0", tRaw[1], RawTensor(e, p[1], q[1]));
    AssertM("Gamma32 t1", tRaw[2],
        MulM(MulM(RawTensor(e, p[2], q[2]), RawLocal(g, e, e)),
             RawLocal(q0, e, e)));
    AssertM("Gamma32 t2", tRaw[3],
        MulM(MulM(RawTensor(e, p[3], q[3]), RawLocal(g, e, e2)),
             RawLocal(q0, e, e2)));

    # Reconstruct the scalar-to-chain dictionary from the original scalar
    # definitions.  L, U, R above are defined independently of RawLocal,
    # RawTensor, RawR, so these checks also audit every Phi-expansion.
    A := SumC([
        L(g, G3Mul(e2, h), e2),
        Sc(-1, L(g, e, h)),
        Sc(-1, L(y0, g, e2))
    ]);
    B0 := SumC([
        U(g, p[3], y0), L(g, g, e2), Sc(-1, L(g, e, g))
    ]);
    S := SumC([
        B0, Sc(-1, A), Sc(-1, L(g, g, h)), L(y0, e, e)
    ]);
    E := AddC(L(y0, e, e), L(y0, e2, e));
    K := AddC(Bar([e, e, e]), Bar([e, e2, e]));
    C1 := SumC([
        R(g, p[3], y0), L(g, g, e2), Sc(-1, L(g, e, g)),
        Sc(-1, U(e, p[3], y0)), Sc(-1, L(g, e, e2))
    ]);
    C2 := SumC([
        R(g, p[2], y1), L(g, g, e), Sc(-1, L(g, e2, g)),
        Sc(-1, U(e2, p[2], y1)), Sc(-1, L(g, e2, e)),
        Sc(-1, L(y0, e2, g)), L(y0, g, e)
    ]);
    T0 := U(e, p[1], q[1]);
    T1 := SumC([U(e, p[2], q[2]), L(g, e, e), L(q0, e, e)]);
    T2 := SumC([U(e, p[3], q[3]), L(g, e, e2), L(q0, e, e2)]);

    AssertM("Gamma32 raw expansion of a", aRaw,
        MulM(MulM(ChainM(A), rhoH), sigmaE2));
    AssertM("Gamma32 raw expansion of b0", b0Raw,
        MulM(ChainM(B0), rhoG));
    AssertM("Gamma32 raw expansion of ell0", ell0Raw,
        MulM(ChainM(AddC(U(e, p[3], y0), L(g, e, e2))), sigmaE));
    AssertM("Gamma32 raw expansion of k2", k2Raw,
        MulM(ChainM(SubC(L(g, g, e2), L(g, e, g))), rhoG));
    AssertM("Gamma32 raw expansion of ell1", ell1Raw,
        MulM(ChainM(SumC([U(e2, p[2], y1), L(g, e2, e),
                          L(y0, e2, g), Sc(-1, L(y0, g, e))])),
             sigmaE));
    AssertM("Gamma32 raw expansion of k1", k1Raw,
        MulM(ChainM(SubC(L(g, g, e), L(g, e2, g))), rhoG));
    AssertM("Gamma32 raw expansion of c1", c1Raw,
        MulM(MulM(RawR(g, p[3], y0), k2Raw), InvM(ell0Raw)));
    AssertM("Gamma32 raw expansion of c2", c2Raw,
        MulM(MulM(RawR(g, p[2], y1), k1Raw), InvM(ell1Raw)));
    AssertM("Gamma32 raw expansion of t0", tRaw[1], ChainM(T0));
    AssertM("Gamma32 raw expansion of t1", tRaw[2], ChainM(T1));
    AssertM("Gamma32 raw expansion of t2", tRaw[3], ChainM(T2));

    # Character relations are consequences of eq:YD-projective-action:
    # rho(gh)=rho(g)rho(h)/Phi_g(g,h),
    # sigma(e)^2=Phi_y0(e,e)sigma(e^2), and
    # sigma(e)^3=Phi_y0(e,e)Phi_y0(e^2,e).
    Eraw := MulM(RawLocal(y0, e, e), RawLocal(y0, e2, e));
    AssertM("Gamma32 raw expansion of the sigma-cubic chain", Eraw, ChainM(E));
    Normalize32 := function(m)
        local out, pos, exponent, remainder, quotient, sigmaE2Expansion,
              rhoGHExpansion;
        out := m;
        rhoGHExpansion := MulM(MulM(rhoG, rhoH),
                               InvM(RawLocal(g, g, h)));
        sigmaE2Expansion := MulM(PowM(sigmaE, 2),
                                 InvM(RawLocal(y0, e, e)));
        out := ReplaceAtom(out, "rho(gh)", rhoGHExpansion);
        out := ReplaceAtom(out, "sigma(epsilon^2)", sigmaE2Expansion);
        # This normalization is used only in the rho(g)=-1 part of the
        # appendix.
        out := ReplaceAtom(out, "rho(g)", NegOneM);
        pos := PositionProperty(out.atoms, x -> x[1] = "sigma(epsilon)");
        if pos <> fail then
            exponent := out.atoms[pos][2];
            remainder := exponent mod 3;
            quotient := (exponent - remainder) / 3;
            out := ReplaceAtom(out, "sigma(epsilon)", OneM);
            out := MulM(out, PowM(sigmaE, remainder));
            out := MulM(out, PowM(Eraw, quotient));
        fi;
        return out;
    end;
    Assert32 := function(label, left, right)
        AssertM(label, Normalize32(left), Normalize32(right));
    end;

    sRaw := MulM(MulM(b0Raw, rhoGH), InvM(aRaw));
    Assert32("Gamma32 dictionary for s", sRaw,
             MulM(PowM(sigmaE, -2), ChainM(S)));
    Assert32("Gamma32 dictionary for sigma(epsilon)^3",
             PowM(sigmaE, 3), ChainM(E));
    Assert32("Gamma32 dictionary for c1", c1Raw,
             MulM(MulM(NegOneM, PowM(sigmaE, -1)), ChainM(C1)));
    Assert32("Gamma32 dictionary for c2", c2Raw,
             MulM(MulM(NegOneM, PowM(sigmaE, -1)), ChainM(C2)));
    Assert32("Gamma32 dictionary for t0", tRaw[1], ChainM(T0));
    Assert32("Gamma32 dictionary for t1", tRaw[2], ChainM(T1));
    Assert32("Gamma32 dictionary for t2", tRaw[3], ChainM(T2));

    Cs := SumC([Sc(3, S), Sc(-2, E), Sc(-1, K)]);
    Clambda := SumC([Sc(3, AddC(AddC(C1, T1), Sc(-1, C2))),
                     Sc(-1, T0), Sc(-1, T1), Sc(-1, T2)]);
    Cr := SumC([Sc(4, C1), Sc(2, T1), Sc(-2, C2),
                Sc(-2, T0), Sc(-1, S)]);
    lambdaRaw := MulM(MulM(c1Raw, tRaw[2]), InvM(c2Raw));
    rRaw := MulM(MulM(PowM(c1Raw, 2), tRaw[2]),
                 InvM(MulM(c2Raw, tRaw[1])));
    kappaRaw := ChainM(K);
    Assert32("Gamma32 comparison-cycle evaluation s^3/kappa",
             MulM(PowM(sRaw, 3), InvM(kappaRaw)), ChainM(Cs));
    RejectPerturbation("Gamma32 dictionary negative control",
                       MulM(PowM(sRaw, 3), InvM(kappaRaw)),
                       ChainM(Cs), Normalize32);
    Assert32("Gamma32 comparison-cycle evaluation lambda^3/(t0t1t2)",
             MulM(PowM(lambdaRaw, 3),
                  InvM(MulM(MulM(tRaw[1], tRaw[2]), tRaw[3]))),
             ChainM(Clambda));
    Assert32("Gamma32 comparison-cycle evaluation r^2/s",
             MulM(PowM(rRaw, 2), InvM(sRaw)), ChainM(Cr));

    Print("[PASS] (3,2): adjoint coefficients and degrees were reconstructed.\n");
    Print("[PASS] (3,2): scalar-chain dictionary and cycle evaluations agree.\n");
    Print("[PASS] (3,2): two one-factor perturbations were rejected.\n");
    return true;
end;

G3_GAMMA332_COEFFICIENT_DICTIONARY_OK :=
    VerifyGamma332CoefficientAndDictionary();
