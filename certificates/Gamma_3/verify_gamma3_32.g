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
