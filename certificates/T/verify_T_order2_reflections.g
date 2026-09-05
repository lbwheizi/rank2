# Exact four-reflection calculation for the order-two T cocycle.
# This is the GAP translation of order2_reflections.py.  All sparse columns
# and actions are recomputed before comparison with the frozen certificate.
# Run from this directory or the repository root.  No packages are required.
if IsExistingFile("T_model_core.g") then Read("T_model_core.g");
else Read("certificates/T/T_model_core.g"); fi;
TLastCheckPassed:=false;;
TRefl := rec();;
TRefl.check := TNew.check;;
TRefl.tphi := function(g,x,y)
    return (TNew.phi(g,x,y)+TNew.phi(TNew.conj(g,x),TNew.conj(g,y),g)
            +TNew.phi(TNew.conj(g,x),g,y)) mod 2;
end;;
TRefl.wordact := function(g,word)
    local out,e,n,i,a,tail,k;
    out:=[]; e:=0;
    for n in [1..Length(word)] do
        i:=word[n]; a:=TNew.act(g,i); Add(out,a[1]); e:=e+a[2];
        tail:=TNew.one;
        for k in word{[n+1..Length(word)]} do tail:=TNew.mul(tail,TNew.xs[k+1]); od;
        e:=e+TRefl.tphi(g,TNew.xs[i+1],tail);
    od;
    return [out,e mod 2];
end;;
TRefl.Y := List([0..3],i->[[[i],1]]);;
for TReflLoopM in [2,3] do
    TRefl.Y:=TNew.independent(List(TNew.prefix(TRefl.Y),TNew.apply_recursive));
od;
TRefl.check(Length(TRefl.Y)=2,"initial third adjoint has dimension two");;
TRefl.rows:=Union(List(TRefl.Y,v->List(v,e->e[1])));;
TRefl.minorFound:=false;;
for TReflLoopRow0 in TRefl.rows do
    TRefl.r0:=TReflLoopRow0;
    for TReflLoopRow1 in TRefl.rows do
        TRefl.r1:=TReflLoopRow1;
        TRefl.aa:=TNew.get(TRefl.Y[1],TRefl.r0);
        TRefl.bb:=TNew.get(TRefl.Y[2],TRefl.r0);
        TRefl.cc:=TNew.get(TRefl.Y[1],TRefl.r1);
        TRefl.dd:=TNew.get(TRefl.Y[2],TRefl.r1);
        TRefl.det:=TRefl.aa*TRefl.dd-TRefl.bb*TRefl.cc;
        if TRefl.det<>0 then TRefl.minorFound:=true; break; fi;
    od;
    if TRefl.minorFound then break; fi;
od;
TRefl.check(TRefl.minorFound,"independent rows for third adjoint");;
TRefl.coords := function(v)
    local x,y,out,test,k,term;
    x:=TNew.get(v,TRefl.r0); y:=TNew.get(v,TRefl.r1);
    out:=[(TRefl.dd*x-TRefl.bb*y)/TRefl.det,
          (TRefl.aa*y-TRefl.cc*x)/TRefl.det];
    test:=[];
    for k in [1,2] do
        for term in TRefl.Y[k] do TNew.accumulate(test,term[1],out[k]*term[2]); od;
    od;
    TRefl.check(test=v,"third-adjoint coordinate reconstruction");
    return out;
end;;
TRefl.uCache:=rec();;
TRefl.uact := function(g,i)
    local key,v,term,a,out;
    key:=String([g,i]);
    if IsBound(TRefl.uCache.(key)) then return TRefl.uCache.(key); fi;
    v:=[];
    for term in TRefl.Y[i+1] do
        a:=TRefl.wordact(g,term[1]);
        TNew.accumulate(v,a[1],(-1)^a[2]*term[2]);
    od;
    out:=TRefl.coords(v); TRefl.uCache.(key):=out; return out;
end;;
TRefl.dact := function(g,i)
    local a;
    a:=TNew.act(g,i);
    return [a[1],(a[2]+TRefl.tphi(g,TNew.inv(TNew.xs[i+1]),TNew.xs[i+1])) mod 2];
end;;
TRefl.verifyInitial := function()
    local g,h,i,j,k,a,b,c,left,u,v;
    for g in TNew.G do for h in TNew.G do
        for i in [0..3] do
            a:=TRefl.dact(h,i); b:=TRefl.dact(g,a[1]); c:=TRefl.dact(TNew.mul(g,h),i);
            TRefl.check(b[1]=c[1] and (a[2]+b[2]) mod 2=
                (TNew.alpha(TNew.inv(TNew.xs[i+1]),g,h)+c[2]) mod 2,
                "dual D projective action");
        od;
        for i in [0,1] do
            left:=[0,0]; u:=TRefl.uact(h,i);
            for j in [0,1] do
                v:=TRefl.uact(g,j);
                for k in [0,1] do left[k+1]:=left[k+1]+u[j+1]*v[k+1]; od;
            od;
            TRefl.check(left=TRefl.uact(TNew.mul(g,h),i),"U group action");
        od;
    od; od;
    Print("U and dual D projective actions verified\n");
end;;
TRefl.verifyInitial();;
TRefl.extraDegree:=rec();; TRefl.extraAct:=rec();;
TRefl.rankCertificate:=[];; TRefl.atomCache:=rec();;
TRefl.degree := function(atom)
    if IsBound(TRefl.extraDegree.(atom[1])) then
        return TRefl.extraDegree.(atom[1])[atom[2]+1];
    fi;
    if atom[1]="U" then return [TNew.one,1]; fi;
    TRefl.check(atom[1]="D","registered atom type");
    return [TNew.inv(TNew.xs[atom[2]+1]),0];
end;;
TRefl.atomact := function(g,t,atom)
    local key,typ,i,out,j,c,a,u;
    key:=String([g,t,atom]);
    if IsBound(TRefl.atomCache.(key)) then return TRefl.atomCache.(key); fi;
    typ:=atom[1]; i:=atom[2];
    if IsBound(TRefl.extraAct.(typ)) then out:=TRefl.extraAct.(typ)(g,t,i);
    elif typ="U" then
        out:=[]; u:=TRefl.uact(g,i);
        for j in [0,1] do
            c:=u[j+1]*TNew.p^(2*t);
            if c<>0 then Add(out,[[typ,j],c]); fi;
        od;
    else
        a:=TRefl.dact(g,i); out:=[[[typ,a[1]],(-1)^a[2]*TNew.p^(-t)]];
    fi;
    TRefl.atomCache.(key):=out; return out;
end;;
TRefl.crossCache:=rec();;
TRefl.cross := function(word,r)
    local key,da,db,tail,atom,term,e,w,out;
    key:=String([word,r]);
    if IsBound(TRefl.crossCache.(key)) then return TRefl.crossCache.(key); fi;
    da:=TRefl.degree(word[r+1]); db:=TRefl.degree(word[r+2]); tail:=TNew.one;
    for atom in word{[r+3..Length(word)]} do tail:=TNew.mul(tail,TRefl.degree(atom)[1]); od;
    out:=[];
    for term in TRefl.atomact(da[1],da[2],word[r+2]) do
        e:=(TNew.phi(da[1],db[1],tail)+TNew.phi(TRefl.degree(term[1])[1],da[1],tail)) mod 2;
        w:=Concatenation(word{[1..r]},[term[1],word[r+1]],word{[r+3..Length(word)]});
        TNew.accumulate(out,w,(-1)^e*term[2]);
    od;
    TRefl.crossCache.(key):=out; return out;
end;;
TRefl.applycross := function(v,r)
    local out,term,after;
    out:=[];
    for term in v do for after in TRefl.cross(term[1],r) do
        TNew.accumulate(out,after[1],term[2]*after[2]);
    od; od;
    return out;
end;;
TRefl.recCache:=rec();;
TRefl.recursive := function(word)
    local key,m,out,double,r,term,a,b;
    key:=String(word);
    if IsBound(TRefl.recCache.(key)) then return TRefl.recCache.(key); fi;
    m:=Length(word)-1; out:=[[word,1]]; double:=[[word,1]];
    for r in Concatenation([0..m-1],[m-1,m-2..0]) do double:=TRefl.applycross(double,r); od;
    for term in double do TNew.accumulate(out,term[1],-term[2]); od;
    if m>1 then
        for a in TRefl.cross(word,0) do
            for b in TRefl.recursive(a[1]{[2..Length(a[1])]}) do
                TNew.accumulate(out,Concatenation([a[1][1]],b[1]),a[2]*b[2]);
            od;
        od;
    fi;
    TRefl.recCache.(key):=out; return out;
end;;
TRefl.applyrec := function(v)
    local out,a,b;
    out:=[];
    for a in v do for b in TRefl.recursive(a[1]) do
        TNew.accumulate(out,b[1],a[2]*b[2]);
    od; od;
    return out;
end;;
TRefl.prefixed := function(A,Y)
    local out,a,v;
    out:=[];
    for a in A do for v in Y do
        Add(out,List(v,e->[Concatenation([a],e[1]),e[2]]));
    od; od;
    return out;
end;;
TRefl.adjoints := function(A,B,maxm,expected)
    local Y,result,m,columns,images;
    Y:=List(B,b->[[[b],1]]); result:=[];
    for m in [1..maxm] do
        columns:=TRefl.prefixed(A,Y); images:=List(columns,TRefl.applyrec);
        Y:=TNew.independent(images);
        Add(TRefl.rankCertificate,rec(A:=A[1][1],B:=B[1][1],m:=m,rank:=Length(Y),
            input_columns:=columns,image_columns:=images,independent_basis:=Y));
        Add(result,Length(Y));
        Print("ADJOINT ",A[1][1]," ",m," ",B[1][1]," dimension ",Length(Y),
              " terms ",Sum(Y,Length),"\n");
        if IsEmpty(Y) then break; fi;
    od;
    TRefl.check(result=expected,Concatenation("adjoint ranks ",A[1][1]," on ",B[1][1]));
    return result;
end;;
TRefl.adjointbasis := function(A,B,m)
    local Y,n;
    Y:=List(B,b->[[[b],1]]);
    for n in [1..m] do Y:=TNew.independent(List(TRefl.prefixed(A,Y),TRefl.applyrec)); od;
    return Y;
end;;
TRefl.tensordegree := function(word)
    local h,t,atom,d;
    h:=TNew.one; t:=0;
    for atom in word do d:=TRefl.degree(atom); h:=TNew.mul(h,d[1]); t:=t+d[2]; od;
    return [h,t mod 6];
end;;
TRefl.tensoract := function(g,t,word)
    local out,n,atom,tail,e,after,a,b;
    out:=[[[],1]];
    for n in [1..Length(word)] do
        atom:=word[n]; tail:=TRefl.tensordegree(word{[n+1..Length(word)]})[1];
        e:=TRefl.tphi(g,TRefl.degree(atom)[1],tail); after:=[];
        for a in out do for b in TRefl.atomact(g,t,atom) do
            TNew.accumulate(after,Concatenation(a[1],[b[1]]),(-1)^e*a[2]*b[2]);
        od; od;
        out:=after;
    od;
    return out;
end;;
TRefl.subobject := function(name,basis)
    local degs,pivots,b,gd,cache;
    degs:=[]; pivots:=[];
    for b in basis do
        gd:=Set(List(b,e->TRefl.tensordegree(e[1])));
        TRefl.check(Length(gd)=1,Concatenation(name," homogeneous basis"));
        Add(degs,gd[1]); Add(pivots,b[1][1]);
    od;
    TRefl.check(Length(Set(degs))=Length(degs),Concatenation(name," distinct degrees"));
    TRefl.extraDegree.(name):=degs; cache:=rec();
    TRefl.extraAct.(name):=function(g,t,i)
        local key,v,a,b,j,coeff,out;
        key:=String([g,t,i]);
        if IsBound(cache.(key)) then return cache.(key); fi;
        v:=[];
        for a in basis[i+1] do for b in TRefl.tensoract(g,t,a[1]) do
            TNew.accumulate(v,b[1],a[2]*b[2]);
        od; od;
        j:=Position(degs,[TNew.conj(g,degs[i+1][1]),degs[i+1][2]]);
        TRefl.check(j<>fail,Concatenation(name," conjugate degree exists"));
        coeff:=TNew.get(v,pivots[j])/TNew.get(basis[j],pivots[j]);
        TRefl.check(v=TNew.scale(basis[j],coeff),Concatenation(name," action on subobject"));
        out:=[[[name,j-1],coeff]]; cache.(key):=out; return out;
    end;
    return List([0..Length(basis)-1],i->[name,i]);
end;;
TRefl.dualsimple := function(name,old)
    local cache;
    TRefl.extraDegree.(name):=List(old,a->[TNew.inv(TRefl.degree(a)[1]),(-TRefl.degree(a)[2]) mod 6]);
    cache:=rec();
    TRefl.extraAct.(name):=function(g,t,i)
        local key,out,j,c,h,e;
        key:=String([g,t,i]);
        if IsBound(cache.(key)) then return cache.(key); fi;
        out:=[];
        for j in [0..Length(old)-1] do
            c:=TNew.get(TRefl.atomact(TNew.inv(g),(-t) mod 6,old[j+1]),old[i+1]);
            if c<>0 then
                h:=TRefl.degree(old[j+1])[1];
                e:=(TNew.alpha(h,g,TNew.inv(g))+TRefl.tphi(g,
                    TNew.inv(TRefl.degree(old[i+1])[1]),TRefl.degree(old[i+1])[1])) mod 2;
                Add(out,[[name,j],(-1)^e*c]);
            fi;
        od;
        cache.(key):=out; return out;
    end;
    return List([0..Length(old)-1],i->[name,i]);
end;;
TRefl.verifyBraid := function(atoms)
    local word,a,b,r;
    for word in Cartesian(atoms,atoms,atoms) do
        a:=[[word,1]]; b:=[[word,1]];
        for r in [0,1,0] do a:=TRefl.applycross(a,r); od;
        for r in [1,0,1] do b:=TRefl.applycross(b,r); od;
        TRefl.check(a=b,"mixed braid relation");
    od;
end;;
TRefl.U:=List([0,1],i->["U",i]);;
TRefl.D:=List([0..3],i->["D",i]);;
TRefl.verifyBraid(Concatenation(TRefl.U,TRefl.D));;
Print("Initial mixed braid relations verified\n");
TRefl.adjoints(TRefl.U,TRefl.D,5,[4,4,0]);;
TRefl.adjoints(TRefl.D,TRefl.U,5,[4,4,1,0]);;
TRefl.E:=TRefl.subobject("E",TRefl.adjointbasis(TRefl.U,TRefl.D,2));;
TRefl.F:=TRefl.dualsimple("F",TRefl.U);;
TRefl.adjoints(TRefl.F,TRefl.E,5,[4,4,0]);;
TRefl.adjoints(TRefl.E,TRefl.F,5,[4,4,1,0]);;
TRefl.K:=TRefl.subobject("K",TRefl.adjointbasis(TRefl.E,TRefl.F,3));;
TRefl.L:=TRefl.dualsimple("L",TRefl.E);;
TRefl.adjoints(TRefl.K,TRefl.L,4,[4,0]);;
TRefl.adjoints(TRefl.L,TRefl.K,5,[4,4,2,0]);;
TRefl.N:=TRefl.subobject("N",TRefl.adjointbasis(TRefl.K,TRefl.L,1));;
TRefl.M:=TRefl.dualsimple("M",TRefl.K);;
TRefl.adjoints(TRefl.M,TRefl.N,4,[4,0]);;
TRefl.adjoints(TRefl.N,TRefl.M,5,[4,4,2,0]);;
TRefl.verifyFourth := function()
    local nindex,scalings,known,i,s,g,a,j,b,target;
    nindex:=List(TRefl.N,a->Position(TNew.xs,TRefl.degree(a)[1])-1);
    scalings:=[1]; known:=[0];
    while Length(known)<4 do
        for i in ShallowCopy(known) do
            s:=scalings[i+1];
            for g in TNew.G do
                a:=TRefl.atomact(g,0,TRefl.N[i+1])[1]; j:=a[1][2];
                b:=TNew.act(g,nindex[i+1]);
                TRefl.check(b[1]=nindex[j+1],"fourth pair rack basis matching");
                target:=s*(-1)^b[2]/a[2];
                if j in known then TRefl.check(scalings[j+1]=target,"fourth pair consistent scalings");
                else scalings[j+1]:=target; Add(known,j); fi;
            od;
        od;
    od;
    for i in [0..3] do for g in TNew.G do
        a:=TRefl.atomact(g,0,TRefl.N[i+1])[1]; j:=a[1][2]; b:=TNew.act(g,nindex[i+1]);
        TRefl.check(b[1]=nindex[j+1] and a[2]*scalings[j+1]=scalings[i+1]*(-1)^b[2],
                    "fourth pair H-action comparison");
    od; od;
    TRefl.check(ForAll(TNew.G,g->TRefl.atomact(g,0,TRefl.M[1])=[[TRefl.M[1],1]]),
                "fourth central summand trivial H-action");
    Print("Fourth pair H-actions equivalent to original\n");
    Print("N indices: ",nindex,"\n");
    Print("N scalings: ",scalings,"\n");
end;;
TRefl.verifyFourth();;
TRefl.composedaction := function(g,t,h,s,a)
    local out,b,c;
    out:=[];
    for b in TRefl.atomact(h,s,a) do for c in TRefl.atomact(g,t,b[1]) do
        TNew.accumulate(out,c[1],b[2]*c[2]);
    od; od;
    return out;
end;;
TRefl.verifyReflectedActions := function()
    local atoms,g,h,a,left,right,pair,old,dual,t,i,j,ev,factor,b,c,expected;
    for atoms in [TRefl.E,TRefl.F,TRefl.K,TRefl.L,TRefl.M,TRefl.N] do
        for g in TNew.G do for h in TNew.G do for a in atoms do
            left:=TRefl.composedaction(g,0,h,0,a);
            right:=TNew.scale(TRefl.atomact(TNew.mul(g,h),0,a),(-1)^TNew.alpha(TRefl.degree(a)[1],g,h));
            TRefl.check(left=right,"reflected projective action");
        od; od; od;
        for g in TNew.G do for a in atoms do
            TRefl.check(TRefl.composedaction(g,0,TNew.one,1,a)=
                TRefl.composedaction(TNew.one,1,g,0,a),"central action commutes");
        od; od;
        for a in atoms do TRefl.check(TRefl.atomact(TNew.one,6,a)=[[a,1]],"central action sixth power"); od;
    od;
    Print("All reflected projective actions verified, including central powers\n");
    for pair in [[TRefl.F,TRefl.E],[TRefl.K,TRefl.L],[TRefl.M,TRefl.N]] do
        TRefl.verifyBraid(Concatenation(pair));
    od;
    Print("All reflected mixed braid relations verified\n");
    for pair in [[TRefl.U,TRefl.F],[TRefl.E,TRefl.L],[TRefl.K,TRefl.M]] do
        old:=pair[1]; dual:=pair[2];
        for g in TNew.G do for t in [0,1] do for i in [1..Length(old)] do for j in [1..Length(old)] do
            ev:=0; factor:=(-1)^TRefl.tphi(g,TRefl.degree(dual[i])[1],TRefl.degree(old[j])[1]);
            for b in TRefl.atomact(g,t,dual[i]) do for c in TRefl.atomact(g,t,old[j]) do
                if b[1][2]=c[1][2] then ev:=ev+factor*b[2]*c[2]; fi;
            od; od;
            expected:=0; if i=j then expected:=1; fi;
            TRefl.check(ev=expected,"dual evaluation pairing");
        od; od; od; od;
    od;
    Print("All dual evaluation pairings verified\n");
end;;
TRefl.verifyReflectedActions();;
TRefl.extraDegree.V:=[[TNew.one,1]];;
TRefl.extraDegree.W:=List(TNew.xs,x->[x,0]);;
TRefl.extraAct.V:=function(g,t,i) return [[["V",0],TNew.p^(-t)]]; end;;
TRefl.extraAct.W:=function(g,t,i)
    local a; a:=TNew.act(g,i); return [[["W",a[1]],(-1)^a[2]*TNew.p^t]];
end;;
TRefl.VV:=[["V",0]];; TRefl.WW:=List([0..3],i->["W",i]);;
TRefl.adjoints(TRefl.VV,TRefl.WW,2,[4,0]);;
TRefl.adjoints(TRefl.WW,TRefl.VV,4,[4,4,2,0]);;
TRefl.verifyRows := function()
    local pairs,centralT,centralZ,tetraT,tetraZ,orientation,certificate,rr,aa,bb,hs,a;
    pairs:=[[TRefl.VV,TRefl.WW],[TRefl.U,TRefl.D],[TRefl.F,TRefl.E],
            [TRefl.K,TRefl.L],[TRefl.M,TRefl.N]];
    centralT:=[1,1,5,5,1]; centralZ:=[-1,2,-2,1,-1];
    tetraT:=[0,0,2,4,3]; tetraZ:=[1,-1,3,3,4]; orientation:=[1,-1,-1,1,1];
    certificate:=[];
    for rr in [1..5] do
        aa:=pairs[rr][1]; bb:=pairs[rr][2];
        TRefl.check(Length(aa)=[1,2,2,1,1][rr] and Length(bb)=4,"five table rows: dimensions");
        TRefl.check(ForAll(aa,a->TRefl.degree(a)=[TNew.one,centralT[rr]]),"five table rows: central degree");
        hs:=TNew.xs; if orientation[rr]=-1 then hs:=List(TNew.xs,TNew.inv); fi;
        TRefl.check(Set(List(bb,TRefl.degree))=Set(List(hs,x->[x,tetraT[rr]])),"five table rows: tetrahedral degrees");
        TRefl.check(ForAll(aa,a->TRefl.atomact(TNew.one,1,a)=[[a,TNew.p^centralZ[rr]]]),"five table rows: central eigenvalues");
        TRefl.check(ForAll(bb,a->TRefl.atomact(TNew.one,1,a)=[[a,TNew.p^tetraZ[rr]]]),"five table rows: tetrahedral central eigenvalues");
        TRefl.check(ForAll(bb,a->TRefl.cross([a,a],0)=[[[a,a],-1]]),"five table rows: diagonal self-braiding");
        for a in Concatenation(aa,bb) do
            Add(certificate,rec(pair:=rr-1,atom:=a,degree:=TRefl.degree(a),
                H_actions:=List(TNew.G,g->[g,TRefl.atomact(g,0,a)]),
                z_action:=TRefl.atomact(TNew.one,1,a)));
        od;
    od;
    return certificate;
end;;
TRefl.actionCertificate:=TRefl.verifyRows();;
if IsExistingFile("order2_reflection_certificate_data.g") then
    Read("order2_reflection_certificate_data.g");
else Read("certificates/T/order2_reflection_certificate_data.g"); fi;
TRefl.check(TOrder2Expected.coefficient_field="Q[p]/(p^2-p+1)","certificate coefficient field");;
TRefl.check(Length(TRefl.rankCertificate)=Length(TOrder2Expected.ranks),"number of rank records");;
for TReflLoopRecord in [1..Length(TRefl.rankCertificate)] do
    TRefl.check(TRefl.rankCertificate[TReflLoopRecord]=TOrder2Expected.ranks[TReflLoopRecord],
        Concatenation("full sparse rank certificate record ",String(TReflLoopRecord)));
od;
TRefl.check(TRefl.actionCertificate=TOrder2Expected.actions,"full action certificate");;
Print("All five rows asserted; 32 complete rank records and 27 action records match\n");
TLastCheckPassed:=true;;
Print("PASS: verify_T_order2_reflections.g\n");
