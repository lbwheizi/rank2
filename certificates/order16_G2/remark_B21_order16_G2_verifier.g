# Paper correspondence: Remark B.21, label rem:order16-G2-certificate.
# File: remark_B21_order16_G2_verifier.g
# Core-only exact verifier for the explicit order-16 Gamma_2 example.
#
# Scalars are represented by exponents modulo 4: e means i^e.  The
# displayed 3-cocycle has values +/-1, represented by exponents 0 and 2.
# No optional GAP package and no floating-point arithmetic is used.

CertCheckCount := 0;;
CertRequire := function(condition,message)
    if not condition then
        Error(message);
    fi;
    CertCheckCount := CertCheckCount+1;
end;;

CertOne := [0,0,0];;
CertEps := [1,0,0];;
CertG := [0,1,0];;
CertH := [0,0,1];;
CertElements := Cartesian([0..1],[0..1],[0..3]);;

CertMul := function(x,y)
    return [ (x[1]+y[1]+(x[3] mod 2)*y[2]) mod 2,
             (x[2]+y[2]) mod 2,
             (x[3]+y[3]) mod 4 ];
end;;

CertInv := function(x)
    return [ (-x[1]+(x[3] mod 2)*x[2]) mod 2,
             (-x[2]) mod 2,
             (-x[3]) mod 4 ];
end;;

CertProd := function(list)
    local out,x;
    out := CertOne;
    for x in list do
        out := CertMul(out,x);
    od;
    return out;
end;;

CertConj := function(a,x)
    return CertProd([a,x,CertInv(a)]);
end;;

CertCentralizer := function(x)
    return Filtered(CertElements,
        a -> CertMul(a,x)=CertMul(x,a));
end;;

CertGeneratedSet := function(generators)
    local result,changed,old,a,b,z;
    result := [CertOne];
    changed := true;
    while changed do
        changed := false;
        old := ShallowCopy(result);
        for a in old do
            for b in Concatenation(old,generators) do
                for z in [CertMul(a,b),CertMul(b,a)] do
                    if not z in result then
                        Add(result,z);
                        changed := true;
                    fi;
                od;
            od;
        od;
    od;
    return Set(result);
end;;

# omega from equation (order16-G2-cocycle), evaluated in F_2.
CertOmega := function(x,y,z)
    local e,a,b,c,d,u,t,barb,bard,bart,hatt;
    e := x[1]; a := x[2]; b := x[3];
    c := y[2]; d := y[3];
    u := z[2]; t := z[3];
    barb := b mod 2;
    bard := d mod 2;
    bart := t mod 2;
    hatt := QuoInt(t-bart,2) mod 2;
    return (c*(a*u+bart*(e+a*barb+a*bard)+a*hatt)) mod 2;
end;;

# Exponents modulo 4 for Phi_x(a,b) and Phi^x(a,b).
CertLocalBit := function(x,a,b)
    return (CertOmega(a,b,x)
            +CertOmega(CertConj(CertMul(a,b),x),a,b)
            -CertOmega(a,CertConj(b,x),b)) mod 2;
end;;

CertLocalExp := function(x,a,b)
    return (2*CertLocalBit(x,a,b)) mod 4;
end;;

CertUpperBit := function(x,a,b)
    return (CertOmega(x,a,b)
            +CertOmega(CertConj(x,a),CertConj(x,b),x)
            -CertOmega(CertConj(x,a),x,b)) mod 2;
end;;

CertUpperExp := function(x,a,b)
    return (2*CertUpperBit(x,a,b)) mod 4;
end;;

CertMapPosition := function(map,x)
    return Position(List(map,pair -> pair[1]),x);
end;;

CertMapValue := function(map,x)
    local position;
    position := CertMapPosition(map,x);
    CertRequire(position<>fail,
        Concatenation("character value missing at ",String(x)));
    return map[position][2];
end;;

CertAssignMapValue := function(map,x,value)
    local position;
    value := value mod 4;
    position := CertMapPosition(map,x);
    if position=fail then
        Add(map,[x,value]);
        return true;
    fi;
    CertRequire(map[position][2]=value,
        Concatenation("inconsistent character value at ",String(x)));
    return false;
end;;

# Extend the three displayed generator values to the whole centralizer,
# using chi(a)chi(b)=Phi_s(a,b)chi(ab), and check every pair exactly.
CertBuildCharacter := function(s,generators,exponents)
    local map,i,changed,known,pairA,pairB,a,b,ab,value,centralizer;
    map := [[CertOne,0]];
    for i in [1..Length(generators)] do
        CertAssignMapValue(map,generators[i],exponents[i]);
    od;
    changed := true;
    while changed do
        changed := false;
        known := ShallowCopy(map);
        for pairA in known do
            for pairB in known do
                a := pairA[1]; b := pairB[1];
                ab := CertMul(a,b);
                value := (pairA[2]+pairB[2]-CertLocalExp(s,a,b)) mod 4;
                if CertAssignMapValue(map,ab,value) then
                    changed := true;
                fi;
            od;
        od;
    od;
    centralizer := Set(CertCentralizer(s));
    CertRequire(Set(List(map,pair -> pair[1]))=centralizer,
        Concatenation("signature does not generate C_G(",String(s),")"));
    for a in centralizer do
        for b in centralizer do
            CertRequire(
                (CertMapValue(map,a)+CertMapValue(map,b)
                 -CertLocalExp(s,a,b)-CertMapValue(map,CertMul(a,b))) mod 4=0,
                "projective-character relation failed");
        od;
    od;
    return map;
end;;

CertRhoInitialExp := function(x)
    CertRequire((x[3] mod 2)=0,"rho evaluated outside C_G(g)");
    return x[2] mod 4;
end;;

CertSigmaInitialExp := function(x)
    CertRequire(x[2]=0,"sigma evaluated outside C_G(h)");
    return (2*(x[1]+(x[3] mod 2))) mod 4;
end;;

CertThetaExp := function(g,h)
    local e,eg,eh,gh,egh,g2,out;
    e := CertEps;
    eg := CertMul(e,g); eh := CertMul(e,h);
    gh := CertMul(g,h); egh := CertProd([e,g,h]);
    g2 := CertMul(g,g);
    out := 0;
    out := out+CertLocalExp(h,e,e);
    out := out+CertLocalExp(h,eg,g);
    out := out+CertUpperExp(h,g,gh);
    out := out+CertUpperExp(h,eg,eh);
    out := out+2*CertLocalExp(g,h,g);
    out := out-CertLocalExp(g,h,e);
    out := out-CertLocalExp(h,e,g2);
    out := out-CertLocalExp(h,e,h);
    out := out-CertLocalExp(h,g,eh);
    out := out-CertUpperExp(g,g,h);
    out := out+CertLocalExp(h,h,g);
    out := out+CertLocalExp(g,e,g);
    out := out+CertLocalExp(g,h,eg);
    out := out-CertUpperExp(h,eg,egh);
    out := out-CertUpperExp(g,eg,eh);
    out := out-CertUpperExp(h,eg,h);
    out := out-CertLocalExp(g,g,h);
    out := out-CertLocalExp(h,g,g);
    return out mod 4;
end;;

CertDeltaExp := function(p,q,r,s)
    local out,ep,p2;
    ep := CertMul(CertEps,p); p2 := CertMul(p,p);
    out := r[1]+s[1]+r[3]+s[3];
    out := out+CertLocalExp(q,ep,p);
    out := out+CertLocalExp(p,q,q);
    out := out-CertLocalExp(p,q,CertEps);
    out := out-CertLocalExp(q,CertEps,p2);
    return out mod 4;
end;;

CertDeltaPrimeExp := function(p,q,r,s)
    local out,eq,q2;
    eq := CertMul(CertEps,q); q2 := CertMul(q,q);
    out := r[1]+s[1]+r[3]+s[3];
    out := out+CertLocalExp(p,eq,q);
    out := out+CertLocalExp(q,p,p);
    out := out-CertLocalExp(q,p,CertEps);
    out := out-CertLocalExp(p,CertEps,q2);
    return out mod 4;
end;;

CertR1Support := function(p,q)
    return [CertInv(p),CertProd([p,p,p,q])];
end;;

CertR2Support := function(p,q)
    return [CertMul(q,p),CertInv(q)];
end;;

CertRepValue := function(representatives,x)
    local position;
    position := Position(List(representatives,pair -> pair[1]),x);
    CertRequire(position<>fail,
        Concatenation("missing conjugacy representative for ",String(x)));
    return representatives[position][2];
end;;

# Equation (order16-induced-action), in i-exponents.
CertInducedActionExp := function(s,character,representatives,A,z)
    local target,rz,rt,c;
    target := CertConj(A,z);
    rz := CertRepValue(representatives,z);
    rt := CertRepValue(representatives,target);
    c := CertProd([CertInv(rt),A,rz]);
    CertRequire(c in CertCentralizer(s),
        "induced-action decomposition did not land in the centralizer");
    return (CertLocalExp(s,A,rz)-CertLocalExp(s,rt,c)
            +CertMapValue(character,c)) mod 4;
end;;

CertRootBraidingMatrix := function(s,conjugator,character)
    local es,support,reps,row,A,z,matrix;
    es := CertMul(CertEps,s);
    CertRequire(CertConj(conjugator,s)=es,
        "listed opposite support element is not the required conjugator");
    support := [s,es];
    reps := [[s,CertOne],[es,conjugator]];
    matrix := [];
    for A in support do
        row := [];
        for z in support do
            Add(row,CertInducedActionExp(s,character,reps,A,z));
        od;
        Add(matrix,row);
    od;
    return matrix;
end;;

# The twelve rows printed in the manuscript.  R1/R2 use GAP's 1-based
# indices.  The script checks support transformations and internal
# projective-character consistency; it does not re-run recursive adjoints.
CertStates := [
 rec(p:=[0,1,0],q:=[0,0,1],r:=[0,1,0],s:=[2,2,0],R1:=2,R2:=3),
 rec(p:=[0,1,0],q:=[0,1,1],r:=[0,1,2],s:=[0,2,0],R1:=1,R2:=4),
 rec(p:=[1,1,1],q:=[0,0,3],r:=[2,1,0],s:=[2,2,2],R1:=5,R2:=1),
 rec(p:=[1,0,1],q:=[1,1,3],r:=[0,1,2],s:=[2,2,0],R1:=6,R2:=2),
 rec(p:=[0,1,3],q:=[0,1,2],r:=[0,1,0],s:=[0,2,2],R1:=3,R2:=7),
 rec(p:=[1,0,3],q:=[1,1,2],r:=[0,1,0],s:=[2,2,2],R1:=4,R2:=8),
 rec(p:=[0,0,1],q:=[0,1,2],r:=[2,1,0],s:=[0,2,0],R1:=9,R2:=5),
 rec(p:=[0,1,1],q:=[1,1,2],r:=[0,1,0],s:=[2,2,2],R1:=10,R2:=6),
 rec(p:=[0,0,3],q:=[1,1,1],r:=[2,1,0],s:=[2,2,2],R1:=7,R2:=11),
 rec(p:=[1,1,3],q:=[1,0,1],r:=[2,1,2],s:=[0,2,0],R1:=8,R2:=12),
 rec(p:=[1,1,0],q:=[0,1,3],r:=[2,1,2],s:=[0,2,0],R1:=12,R2:=9),
 rec(p:=[1,1,0],q:=[1,0,3],r:=[2,1,2],s:=[0,2,0],R1:=11,R2:=10)
];;

# -------------------------------------------------------------------------
# Exhaustive group checks.
# -------------------------------------------------------------------------
CertRequire(Length(CertElements)=16,"the normal-form set does not have order 16");
CertRequire(Length(Set(CertElements))=16,"normal forms are not unique");
for x in CertElements do
    CertRequire(CertMul(CertOne,x)=x and CertMul(x,CertOne)=x,
        "identity law failed");
    CertRequire(CertMul(x,CertInv(x))=CertOne
                and CertMul(CertInv(x),x)=CertOne,
        "inverse law failed");
    for y in CertElements do
        CertRequire(CertMul(x,y) in CertElements,"closure failed");
        for z in CertElements do
            CertRequire(CertMul(CertMul(x,y),z)=CertMul(x,CertMul(y,z)),
                "associativity failed");
        od;
    od;
od;
CertRequire(CertGeneratedSet([CertG,CertH])=Set(CertElements),
    "g and h do not generate the displayed group");
CertRequire(CertMul(CertH,CertG)=CertProd([CertEps,CertG,CertH]),
    "hg=epsilon gh failed");
CertRequire(CertProd([CertEps,CertEps])=CertOne
            and CertProd([CertG,CertG])=CertOne
            and CertProd([CertH,CertH,CertH,CertH])=CertOne,
    "displayed power relations failed");
CertRequire(ForAll(CertElements,x -> CertMul(CertEps,x)=CertMul(x,CertEps)),
    "epsilon is not central");

# -------------------------------------------------------------------------
# Exhaustive normalized 3-cocycle check.
# -------------------------------------------------------------------------
for x in CertElements do
    for y in CertElements do
        CertRequire(CertOmega(CertOne,x,y)=0
                    and CertOmega(x,CertOne,y)=0
                    and CertOmega(x,y,CertOne)=0,
            "3-cocycle normalization failed");
    od;
od;
for a in CertElements do
    for b in CertElements do
        for c in CertElements do
            for d in CertElements do
                CertRequire((CertOmega(b,c,d)
                    +CertOmega(a,CertMul(b,c),d)
                    +CertOmega(a,b,c)
                    -CertOmega(CertMul(a,b),c,d)
                    -CertOmega(a,b,CertMul(c,d))) mod 2=0,
                    "normalized 3-cocycle identity failed");
            od;
        od;
    od;
od;

# -------------------------------------------------------------------------
# Initial centralizers, local cocycles, and explicit characters.
# -------------------------------------------------------------------------
CertCG := CertCentralizer(CertG);;
CertCH := CertCentralizer(CertH);;
CertRequire(CertCG=Filtered(CertElements,x -> (x[3] mod 2)=0),
    "the displayed formula for C_G(g) failed");
CertRequire(CertCH=Filtered(CertElements,x -> x[2]=0),
    "the displayed formula for C_G(h) failed");
for a in CertCG do
    for b in CertCG do
        CertRequire(CertLocalExp(CertG,a,b)=(2*a[2]*b[2]) mod 4,
            "the closed formula for Phi_g failed");
        CertRequire((CertRhoInitialExp(a)+CertRhoInitialExp(b)
                     -CertLocalExp(CertG,a,b)
                     -CertRhoInitialExp(CertMul(a,b))) mod 4=0,
            "rho projective-character relation failed");
    od;
od;
for a in CertCH do
    for b in CertCH do
        CertRequire(CertLocalExp(CertH,a,b)=0,
            "the closed formula for Phi_h failed");
        CertRequire((CertSigmaInitialExp(a)+CertSigmaInitialExp(b)
                     -CertLocalExp(CertH,a,b)
                     -CertSigmaInitialExp(CertMul(a,b))) mod 4=0,
            "sigma projective-character relation failed");
    od;
od;
CertRequire([CertRhoInitialExp(CertEps),CertRhoInitialExp(CertG),
             CertRhoInitialExp(CertMul(CertH,CertH))]=[0,1,0],
    "initial (zeta,lambda,mu) is not (1,i,1)");
CertRequire([CertSigmaInitialExp(CertEps),CertSigmaInitialExp(CertH),
             CertSigmaInitialExp(CertMul(CertG,CertG))]=[2,2,0],
    "initial (zeta',lambda',mu') is not (-1,-1,1)");

# -------------------------------------------------------------------------
# Twelve-state table: reconstruct each character from its signature, check
# all projective relations, support reflections, local scalars, and labels.
# -------------------------------------------------------------------------
for j in [1..Length(CertStates)] do
    state := CertStates[j];
    p2 := CertMul(state.p,state.p);
    q2 := CertMul(state.q,state.q);
    CertRequire(CertGeneratedSet([CertEps,state.p,q2])
                =Set(CertCentralizer(state.p)),
        "rho centralizer generators in a state are incorrect");
    CertRequire(CertGeneratedSet([CertEps,state.q,p2])
                =Set(CertCentralizer(state.q)),
        "sigma centralizer generators in a state are incorrect");
    state.rhoCharacter := CertBuildCharacter(state.p,
        [CertEps,state.p,q2],state.r);
    state.sigmaCharacter := CertBuildCharacter(state.q,
        [CertEps,state.q,p2],state.s);

    r1support := CertR1Support(state.p,state.q);
    r2support := CertR2Support(state.p,state.q);
    CertRequire(r1support=[CertStates[state.R1].p,CertStates[state.R1].q],
        "an R1 support transition in the table failed");
    CertRequire(r2support=[CertStates[state.R2].p,CertStates[state.R2].q],
        "an R2 support transition in the table failed");

    CertRequire(state.r[2]=1,"a listed lambda_j is not i");
    CertRequire(state.s[2]=2,"a listed lambda'_j is not -1");
    CertRequire(CertDeltaExp(state.p,state.q,state.r,state.s)=0,
        "a listed Delta_j is not 1");
    CertRequire(CertDeltaPrimeExp(state.p,state.q,state.r,state.s)=0,
        "a listed Delta'_j is not 1");
    CertRequire(CertThetaExp(state.p,state.q)=2,
        "a listed Theta_Phi,j is not -1");

    firstMatrix := CertRootBraidingMatrix(
        state.p,state.q,state.rhoCharacter);
    secondMatrix := CertRootBraidingMatrix(
        state.q,state.p,state.sigmaCharacter);
    CertRequire(firstMatrix[1][1]=1 and firstMatrix[2][2]=1
                and (firstMatrix[1][2]+firstMatrix[2][1]) mod 4=0,
        "first-orbit root diagonal labels failed");
    CertRequire(secondMatrix[1][1]=2 and secondMatrix[2][2]=2
                and (secondMatrix[1][2]+secondMatrix[2][1]) mod 4=2,
        "second-orbit root diagonal labels failed");
od;

# Row zero reconstructs the two explicitly defined initial characters.
for x in CertCG do
    CertRequire(CertMapValue(CertStates[1].rhoCharacter,x)
                =CertRhoInitialExp(x),
        "state P0 rho does not equal the explicit initial rho");
od;
for x in CertCH do
    CertRequire(CertMapValue(CertStates[1].sigmaCharacter,x)
                =CertSigmaInitialExp(x),
        "state P0 sigma does not equal the explicit initial sigma");
od;

# Both displayed transitions are involutions, and all twelve rows are
# reachable from P0 under them.
for j in [1..Length(CertStates)] do
    CertRequire(CertStates[CertStates[j].R1].R1=j,
        "R1 is not involutive on the displayed table");
    CertRequire(CertStates[CertStates[j].R2].R2=j,
        "R2 is not involutive on the displayed table");
od;
reachable := [1];;
changed := true;;
while changed do
    changed := false;
    for j in ShallowCopy(reachable) do
        for k in [CertStates[j].R1,CertStates[j].R2] do
            if not k in reachable then
                Add(reachable,k);
                changed := true;
            fi;
        od;
    od;
od;
CertRequire(Set(reachable)=[1..12],
    "the displayed reflection table is not connected from P0");

# Root-factor label arithmetic recorded in the text.
CertFirstRootFactorDimension := 4^2;;
CertSecondRootFactorDimension := 2^3;;
CertDimensionProduct := CertFirstRootFactorDimension^3
                        *CertSecondRootFactorDimension^3;;
CertRequire(CertFirstRootFactorDimension=16,
    "first root-factor dimension arithmetic failed");
CertRequire(CertSecondRootFactorDimension=8,
    "second root-factor dimension arithmetic failed");
CertRequire(CertDimensionProduct=2^21,
    "G2 root-factor dimension product failed");

Print("GAP version: ",GAPInfo.Version,"\n");
Print("Scope: explicit order-16 data; core GAP; exact i-exponents mod 4.\n");
Print("Group: 16 normal forms; 4096 associativity triples; inverses verified.\n");
Print("Cocycle: normalized; 65536 four-tuples verified.\n");
Print("Initial characters: rho/sigma projective relations verified.\n");
Print("Initial centralizers: 8 elements each.\n");
Print("State table: 12 rows; 24 reconstructed projective characters.\n");
Print("Reflection supports: 24 edges; involutions and connectivity verified.\n");
Print("Uniform state scalars: lambda=i, lambda'=-1, Delta=Delta'=1, Theta=-1.\n");
Print("First-orbit root labels: (i,i), cross product 1.\n");
Print("Second-orbit root labels: (-1,-1), cross product -1.\n");
Print("Root-factor arithmetic: 16^3 * 8^3 = ",CertDimensionProduct," = 2^21.\n");
Print("NOT CHECKED: recursive adjoint coefficients or terminal characters.\n");
Print("NOT CHECKED: adjoint simplicity or vanishing.\n");
Print("Exact assertions passed: ",CertCheckCount,"\n");
Print("PASS: all implemented finite checks succeeded.\n");

quit;
