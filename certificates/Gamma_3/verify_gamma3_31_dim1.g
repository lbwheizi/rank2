# The one-dimensional central-support calculation uses exactly the same
# normalized bar chains E,A,B,T0,T1,T2,M0,M1,M2,Theta,K,C02,Cpow as the
# (3,1)_2 calculation.  Reading the shared verifier first rechecks those
# chains and their explicit S3 four-chain boundaries.
Read("verify_gamma3_31_dim2.g");

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

    if not RatEq(r01,rhs01) then Error("(3,1)_1 m0/m1 reduction"); fi;
    if not RatEq(r02,rhs02) then Error("(3,1)_1 m0/m2 reduction"); fi;
    if not RatEq(r12,rhs12) then Error("(3,1)_1 m1/m2 reduction"); fi;

    # Exact polynomial consequences used after m0=m1=m2.
    if t*(1-t)-1 <> -(1-t+t^2) then
        Error("(3,1)_1 t(1-t)=1 reduction");
    fi;
    if t^3+1 <> (t+1)*(t^2-t+1) then
        Error("(3,1)_1 t^3=-1 reduction");
    fi;
    if (p*t)^3-1 <> (p*t-1)*((p*t)^2+p*t+1) then
        Error("(3,1)_1 pt=1 cubic reduction");
    fi;

    Print("[PASS] (3,1)_1: shared C02/Cpow bar certificates were rechecked.\n");
    Print("[PASS] (3,1)_1: all three m_i ratios and branch polynomials agree exactly.\n");
end;

VerifyGamma331Dim1ScalarReduction();
