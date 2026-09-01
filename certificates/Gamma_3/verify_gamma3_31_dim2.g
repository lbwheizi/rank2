Read("bar_complex.g");

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
