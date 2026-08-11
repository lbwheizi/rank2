# Exact core-only check of the F_3-valued normalized S_3 3-cocycle
# used as a detector for the kappa class.

VerifyS3CocycleMod3 := function()
    local one, e, e2, g, eg, e2g, elts, mul, support, value,
          a, b, c, d, delta, kappa;

    one := [0,0,0];
    e   := [1,0,0];
    e2  := [2,0,0];
    g   := [0,1,0];
    eg  := [1,1,0];
    e2g := [2,1,0];
    elts := [one,e,e2,g,eg,e2g];

    mul := function(x, y)
        local sign;
        if x[2] mod 2 = 0 then sign := 1; else sign := -1; fi;
        return [(x[1] + sign*y[1]) mod 3,
                (x[2] + y[2]) mod 2, 0];
    end;

    # Entries are [ [a,b,c], f(a,b,c) ] with values in F_3.
    support := [
        [[e,e,e],1],       [[e,e,e2],2],      [[e,e,g],1],       [[e,e,eg],2],
        [[e,g,e],2],       [[e,g,e2],1],      [[e,g,g],2],       [[e,g,eg],1],
        [[e2,e2,e],2],     [[e2,e2,e2],1],    [[e2,e2,g],2],     [[e2,e2,eg],1],
        [[e2,eg,e],1],     [[e2,eg,e2],2],    [[e2,eg,g],1],     [[e2,eg,eg],2],
        [[g,e2,e],2],      [[g,e2,e2],1],     [[g,e2,g],2],      [[g,e2,eg],1],
        [[g,eg,e],1],      [[g,eg,e2],2],     [[g,eg,g],1],      [[g,eg,eg],2],
        [[eg,e,e],1],      [[eg,e,e2],2],     [[eg,e,g],1],      [[eg,e,eg],2],
        [[eg,g,e],2],      [[eg,g,e2],1],     [[eg,g,g],2],      [[eg,g,eg],1]
    ];

    value := function(x, y, z)
        local pos;
        if x=one or y=one or z=one then return 0; fi;
        pos := PositionProperty(support, t -> t[1]=[x,y,z]);
        if pos=fail then return 0; fi;
        return support[pos][2];
    end;

    if Length(support) <> 32 then Error("S3 cocycle support count"); fi;
    for a in elts do
        for b in elts do
            for c in elts do
                for d in elts do
                    delta := value(b,c,d)
                             - value(mul(a,b),c,d)
                             + value(a,mul(b,c),d)
                             - value(a,b,mul(c,d))
                             + value(a,b,c);
                    if delta mod 3 <> 0 then
                        Error("S3 cocycle equation failed at ", [a,b,c,d]);
                    fi;
                od;
            od;
        od;
    od;

    kappa := (value(e,e,e) + value(e,e2,e)) mod 3;
    if kappa <> 1 then Error("S3 detector has wrong kappa exponent"); fi;

    Print("[PASS] S3 detector: all 1296 normalized 3-cocycle equations hold.\n");
    Print("[PASS] S3 detector: support = 32 and kappa exponent = 1 in F_3.\n");
end;

VerifyS3CocycleMod3();
