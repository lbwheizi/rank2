# Exact cocycle and induced-action calculations for the revised T case.
# The infinite T coordinates and formal q, eta, zeta_8 exponents are retained.
# No finite sample of input cocycles is used for the general reduction.

Read("q8_primitive8_cocycle_data.g");;
TCo := rec();;
TCo.check := function(condition, message)
    if not condition then
        PrintTo("*errout*", "FAIL: ", message, "\n");
        ErrorNoReturn("T cocycle verification failed.");
    fi;
end;;
TCo.Q := Tuples([0,1],3);;
TCo.oneq := [0,0,0];;
TCo.qm := function(u,v)
    return [(u[1]+v[1]+u[2]*v[2]+u[3]*v[3]+u[3]*v[2]) mod 2,
            (u[2]+v[2]) mod 2, (u[3]+v[3]) mod 2];
end;;
TCo.A := function(q,r)
    local a,b,n,t;
    a:=q[2]; b:=q[3];
    for n in [1..(r mod 3)] do
        t:=a; a:=b; b:=(t+b) mod 2;
    od;
    return [q[1],a,b];
end;;
TCo.qinv := q -> First(TCo.Q,v -> TCo.qm(q,v)=TCo.oneq);;
TCo.F := function(a,b,c)
    return TQuaternionF[Position(TCo.Q,a)][Position(TCo.Q,b)][Position(TCo.Q,c)];
end;;
TCo.E := function(a,b,c)
    local n,result;
    result:=0;
    for n in [1..3] do
        result:=(result+a[2]*b[2]*c[3]) mod 2;
        a:=TCo.A(a,1); b:=TCo.A(b,1); c:=TCo.A(c,1);
    od;
    return result;
end;;
TCo.check(Length(TQuaternionF)=8 and
    ForAll(TQuaternionF,a -> Length(a)=8 and
        ForAll(a,b -> Length(b)=8 and ForAll(b,v -> IsInt(v) and v in [0..7]))),
    "the fixed cocycle has exactly 512 entries in Z/8");
TCo.ok:=true;;
for TCo_a in TCo.Q do
    for TCo_b in TCo.Q do
        TCo.ok:=TCo.ok and
            TCo.A(TCo.qm(TCo_a,TCo_b),1)=TCo.qm(TCo.A(TCo_a,1),TCo.A(TCo_b,1)) and
            TCo.F(TCo.oneq,TCo_a,TCo_b)=0 and
            TCo.F(TCo_a,TCo.oneq,TCo_b)=0 and TCo.F(TCo_a,TCo_b,TCo.oneq)=0;
        for TCo_c in TCo.Q do
            TCo.ok:=TCo.ok and
                TCo.qm(TCo.qm(TCo_a,TCo_b),TCo_c)=TCo.qm(TCo_a,TCo.qm(TCo_b,TCo_c)) and
                TCo.F(TCo.A(TCo_a,1),TCo.A(TCo_b,1),TCo.A(TCo_c,1))=TCo.F(TCo_a,TCo_b,TCo_c) and
                TCo.F(TCo_a,TCo_b,TCo_c) mod 2=TCo.E(TCo_a,TCo_b,TCo_c);
            for TCo_d in TCo.Q do
                TCo.check((TCo.F(TCo_b,TCo_c,TCo_d)-TCo.F(TCo.qm(TCo_a,TCo_b),TCo_c,TCo_d)
                    +TCo.F(TCo_a,TCo.qm(TCo_b,TCo_c),TCo_d)
                    -TCo.F(TCo_a,TCo_b,TCo.qm(TCo_c,TCo_d))+TCo.F(TCo_a,TCo_b,TCo_c)) mod 8=0,
                    "a fixed quaternion cocycle identity");
            od;
        od;
    od;
od;
TCo.check(TCo.ok,"quaternion law, normalization, cyclic invariance and parity");
TCo.i:=[0,1,0];; TCo.j:=[0,0,1];; TCo.k:=[0,1,1];;
TCo.cycle:=[[TCo.i,TCo.j,TCo.k],[TCo.j,TCo.k,TCo.i],[TCo.k,TCo.i,TCo.j]];;
TCo.boundary:=[];;
TCo.addBoundary:=function(key,c)
    local pos;
    pos:=PositionProperty(TCo.boundary,x -> x[1]=key);
    if pos=fail then Add(TCo.boundary,[key,c]);
    else TCo.boundary[pos][2]:=TCo.boundary[pos][2]+c; fi;
end;;
for TCo_term in TCo.cycle do
    TCo_a:=TCo_term[1]; TCo_b:=TCo_term[2]; TCo_c:=TCo_term[3];
    TCo.addBoundary([TCo_b,TCo_c],1);
    TCo.addBoundary([TCo.qm(TCo_a,TCo_b),TCo_c],-1);
    TCo.addBoundary([TCo_a,TCo.qm(TCo_b,TCo_c)],1);
    TCo.addBoundary([TCo_a,TCo_b],-1);
od;
TCo.check(ForAll(TCo.boundary,x -> x[2]=0),"the three-term bar chain has zero boundary");
TCo.values:=List(TCo.cycle,t -> TCo.F(t[1],t[2],t[3]));;
TCo.check(TCo.values=[3,3,3] and Sum(TCo.values) mod 8=1,"bar-cycle values and evaluation");
Print("[OK] Quaternion law, 512 entries, 4096 cocycle identities, and bar cycle.\n");

# T = (Q8 semidirect Z) x Z, with unrestricted integer coordinates.
TCo.one:=[TCo.oneq,0,0];;
TCo.mul:=function(g,h)
    return [TCo.qm(g[1],TCo.A(h[1],g[2])),g[2]+h[2],g[3]+h[3]];
end;;
TCo.inv:=g -> [TCo.A(TCo.qinv(g[1]),-g[2]),-g[2],-g[3]];;
TCo.conj:={g,h} -> TCo.mul(TCo.mul(g,h),TCo.inv(g));;
TCo.phi:={g,h,k} -> TCo.F(g[1],TCo.A(h[1],g[2]),TCo.A(k[1],g[2]+h[2]));;
TCo.alpha:=function(x,g,h)
    return (TCo.phi(g,h,x)+TCo.phi(TCo.conj(TCo.mul(g,h),x),g,h)
            -TCo.phi(g,TCo.conj(h,x),h)) mod 8;
end;;
TCo.x1:=[TCo.oneq,1,0];; TCo.z:=[TCo.oneq,0,1];;
TCo.ts:=[TCo.one,[TCo.i,0,0],[TCo.k,0,0],[TCo.j,0,0]];;
TCo.xs:=List(TCo.ts,t -> TCo.conj(t,TCo.x1));;
TCo.check(TCo.xs=[TCo.x1,[[1,1,1],1,0],[[1,0,1],1,0],[[1,1,0],1,0]],"four support lifts");
TCo.check(TCo.mul(TCo.mul(TCo.xs[1],TCo.xs[2]),TCo.xs[1])=
    TCo.mul(TCo.mul(TCo.xs[2],TCo.xs[1]),TCo.xs[2]),"the support braid relation");
TCo.check(TCo.mul(TCo.mul(TCo.xs[1],TCo.xs[1]),TCo.xs[1])=
    TCo.mul(TCo.mul(TCo.xs[2],TCo.xs[2]),TCo.xs[2]),"the support cube relation");
TCo.C:=[];; TCo.S:=[];; TCo.rack:=[];;
for TCo_g in TCo.xs do
    TCo.cr:=[]; TCo.sr:=[]; TCo.rr:=[];
    for TCo_n in [1..4] do
        TCo_r:=Position(TCo.xs,TCo.conj(TCo_g,TCo.xs[TCo_n]));
        TCo_c:=TCo.mul(TCo.inv(TCo.ts[TCo_r]),TCo.mul(TCo_g,TCo.ts[TCo_n]));
        TCo.check(TCo_c[1]{[2,3]}=[0,0] and TCo_c{[2,3]}=[1,0],"centralizer tail");
        Add(TCo.cr,(TCo.alpha(TCo.x1,TCo_g,TCo.ts[TCo_n])-
            TCo.alpha(TCo.x1,TCo.ts[TCo_r],TCo_c)) mod 8);
        Add(TCo.sr,TCo_c[1][1]); Add(TCo.rr,TCo_r-1);
    od;
    Add(TCo.C,TCo.cr); Add(TCo.S,TCo.sr); Add(TCo.rack,TCo.rr);
od;
TCo.check(TCo.C=[[0,0,0,0],[4,0,7,1],[4,1,0,7],[4,7,1,0]],"induced-action cocycle exponents");
TCo.check(TCo.S=[[0,0,0,0],[1,0,1,0],[1,0,0,1],[1,1,0,0]],"induced-action eta exponents");
Print("[OK] Induced-action exponent matrices over the infinite group T.\n");

# Each triple is the exponent of q, eta, and zeta_8^h, respectively.
TCo.add:={u,v} -> [u[1]+v[1],(u[2]+v[2]) mod 2,(u[3]+v[3]) mod 8];;
TCo.braid:=function(word,r)
    local a,b,out,tail,n,phase,after;
    a:=word[r+1]; b:=word[r+2]; out:=TCo.rack[a+1][b+1]; tail:=TCo.one;
    if r+3<=Length(word) then
        for n in [r+3..Length(word)] do tail:=TCo.mul(tail,TCo.xs[word[n]+1]); od;
    fi;
    phase:=[1,TCo.S[a+1][b+1],(TCo.C[a+1][b+1]+TCo.phi(TCo.xs[a+1],TCo.xs[b+1],tail)
        -TCo.phi(TCo.xs[out+1],TCo.xs[a+1],tail)) mod 8];
    after:=ShallowCopy(word); after[r+1]:=out; after[r+2]:=a;
    return [after,phase];
end;;
TCo.crossings:=function(word,sequence)
    local phase,r,b;
    phase:=[0,0,0];
    for r in sequence do b:=TCo.braid(word,r); word:=b[1]; phase:=TCo.add(phase,b[2]); od;
    return [word,phase];
end;;
TCo.result:=TCo.crossings([1,2],[0,0,0]);;
TCo.check(TCo.result[1]=[1,2] and TCo.add([1,0,0],TCo.result[2])=[4,1,5],"formal epsilon exponents");
TCo.result:=TCo.crossings([0,1,2],[0,1,0,1,0,1]);;
TCo.check(TCo.result=[[0,1,2],[6,1,6]],"formal full-twist exponents");
Print("[OK] epsilon exponents [4,1,5]; full-twist exponents [6,1,6].\n");
TCo.solutions:=[];;
for TCo_h in [0..7] do
    for TCo_s in [0..1] do
        if (5*TCo_h+4*TCo_s) mod 8=0 then Add(TCo.solutions,[TCo_h,(-1)^TCo_s]); fi;
    od;
od;
TCo.check(TCo.solutions=[[0,1],[4,-1]],"the exact two-branch equation");
Print("[OK] q=-1 and epsilon=1 give exactly (h,eta)=(0,1),(4,-1).\n");
TLastCheckPassed:=true;;
Print("PASS: verify_T_cohomology_reduction.g\n");
