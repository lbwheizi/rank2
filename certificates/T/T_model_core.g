# Exact operations for the order-two cocycle model in Appendix D.
# Group entries and word letters retain the zero-based conventions of the paper.
# Coefficients belong to Q(E(6)).  Sparse vectors are sorted [key, coefficient] lists.
# This file has no package dependency and prints nothing when loaded.

if not IsBound(TNew) then
TNew := rec();;
TNew.p := E(6);;
TNew.one := [0,0,0,0];;
TNew.G := Cartesian([0,1],[0,1],[0,1,2],[0,1]);;
TNew.check := function(test, message)
    if not test then
        PrintTo("*errout*", "FAIL: ", message, "\n");
        ErrorNoReturn(message);
    fi;
end;;
TNew.check(TNew.p^2-TNew.p+1=0 and not IsRat(TNew.p),
           "E(6) gives the exact quadratic coefficient field");;
TNew.A := function(u, r)
    local a,b,t,c;
    a:=u[1]; b:=u[2];
    for t in [1..(r mod 3)] do c:=b; b:=(a+b) mod 2; a:=c; od;
    return [a,b];
end;;
TNew.mul := function(g,h)
    local v;
    v:=TNew.A(h{[1,2]},g[3]);
    return [(g[1]+v[1]) mod 2,(g[2]+v[2]) mod 2,
            (g[3]+h[3]) mod 3,
            (g[4]+h[4]+g[1]*v[1]+g[2]*v[2]+g[2]*v[1]) mod 2];
end;;
TNew.inverses := List(TNew.G,g->First(TNew.G,h->TNew.mul(g,h)=TNew.one));;
TNew.inv := g->TNew.inverses[Position(TNew.G,g)];;
TNew.conj := function(g,h) return TNew.mul(TNew.mul(g,h),TNew.inv(g)); end;;
TNew.parity := function(u,v,w)
    local s,t,a,b,c;
    s:=0;
    for t in [0..2] do
        a:=TNew.A(u,t); b:=TNew.A(v,t); c:=TNew.A(w,t);
        s:=s+a[1]*b[1]*c[2];
    od;
    return s mod 2;
end;;
TNew.phi := function(g,h,k)
    return TNew.parity(g{[1,2]},TNew.A(h{[1,2]},g[3]),
                      TNew.A(k{[1,2]},g[3]+h[3]));
end;;
TNew.alpha := function(x,g,h)
    return (TNew.phi(g,h,x)+TNew.phi(TNew.conj(TNew.mul(g,h),x),g,h)
           +TNew.phi(g,TNew.conj(h,x),h)) mod 2;
end;;
TNew.xs := [[0,0,1,1]];;
Add(TNew.xs,TNew.conj([1,0,0,0],TNew.xs[1]));;
Add(TNew.xs,TNew.conj(TNew.xs[2],TNew.xs[1]));;
Add(TNew.xs,TNew.conj(TNew.xs[1],TNew.xs[2]));;
TNew.ts := List(TNew.xs,x->First(TNew.G,t->TNew.conj(t,TNew.xs[1])=x));;
TNew.act := function(g,i)
    local j,c;
    j:=Position(TNew.xs,TNew.conj(g,TNew.xs[i+1]))-1;
    c:=TNew.mul(TNew.inv(TNew.ts[j+1]),TNew.mul(g,TNew.ts[i+1]));
    TNew.check(c[1]=0 and c[2]=0,"induced-action centralizer");
    return [j,(TNew.alpha(TNew.xs[1],g,TNew.ts[i+1])
           +TNew.alpha(TNew.xs[1],TNew.ts[j+1],c)+c[4]) mod 2];
end;;
TNew.rack := List(TNew.xs,g->List([0..3],i->TNew.act(g,i)[1]));;
TNew.beta := List(TNew.xs,g->List([0..3],i->TNew.act(g,i)[2]));;
TNew.braid := function(word,r)
    local i,j,k,tail,a,e,out;
    i:=word[r+1]; j:=word[r+2]; k:=TNew.rack[i+1][j+1];
    tail:=TNew.one;
    for a in word{[r+3..Length(word)]} do tail:=TNew.mul(tail,TNew.xs[a+1]); od;
    e:=(TNew.beta[i+1][j+1]+TNew.phi(TNew.xs[i+1],TNew.xs[j+1],tail)
       +TNew.phi(TNew.xs[k+1],TNew.xs[i+1],tail)) mod 2;
    out:=ShallowCopy(word); out[r+1]:=k; out[r+2]:=i;
    return [out,e];
end;;
TNew.word_action := function(word,bs)
    local r,out,e;
    e:=0;
    for r in bs do out:=TNew.braid(word,r); word:=out[1]; e:=(e+out[2]) mod 2; od;
    return [word,e];
end;;

# Return the first position whose key is at least the requested key.
TNew.position := function(v,key)
    local lo,hi,mid;
    lo:=1; hi:=Length(v)+1;
    while lo<hi do
        mid:=QuoInt(lo+hi,2);
        if v[mid][1]<key then lo:=mid+1; else hi:=mid; fi;
    od;
    return lo;
end;;
TNew.accumulate := function(v,key,c)
    local pos;
    if c=0 then return; fi;
    pos:=TNew.position(v,key);
    if pos<=Length(v) and v[pos][1]=key then
        v[pos][2]:=v[pos][2]+c;
        if v[pos][2]=0 then Remove(v,pos); fi;
    else Add(v,[ShallowCopy(key),c],pos);
    fi;
end;;
TNew.get := function(v,key)
    local pos;
    pos:=TNew.position(v,key);
    if pos<=Length(v) and v[pos][1]=key then return v[pos][2]; fi;
    return 0;
end;;
TNew.scale := function(v,c)
    if c=0 then return []; fi;
    return List(v,e->[e[1],e[2]*c]);
end;;
TNew.independent := function(columns)
    local echelon,basis,source,v,pivot,pos,c,term;
    echelon:=[]; basis:=[];
    for source in columns do
        v:=List(source,e->[e[1],e[2]]);
        while not IsEmpty(v) do
            pivot:=v[1][1]; pos:=TNew.position(echelon,pivot);
            if pos>Length(echelon) or echelon[pos][1]<>pivot then
                Add(echelon,[pivot,TNew.scale(v,1/v[1][2])],pos);
                Add(basis,source);
                break;
            fi;
            c:=v[1][2];
            for term in echelon[pos][2] do TNew.accumulate(v,term[1],-c*term[2]); od;
        od;
    od;
    return basis;
end;;
TNew.recCache := rec();;
TNew.recursive := function(word)
    local key,m,result,out,crossed,term;
    key:=String(word);
    if IsBound(TNew.recCache.(key)) then return TNew.recCache.(key); fi;
    m:=Length(word);
    if m=1 then result:=[[word,1-TNew.p]];
    else
        result:=[[word,1]];
        out:=TNew.word_action(word,Concatenation([0..m-2],[m-2,m-3..0]));
        TNew.accumulate(result,out[1],-(-1)^out[2]*TNew.p);
        crossed:=TNew.braid(word,0);
        for term in TNew.recursive(crossed[1]{[2..m]}) do
            TNew.accumulate(result,Concatenation([crossed[1][1]],term[1]),
                            (-1)^crossed[2]*term[2]);
        od;
    fi;
    TNew.recCache.(key):=result;
    return result;
end;;
TNew.apply_recursive := function(v)
    local out,term,after;
    out:=[];
    for term in v do
        for after in TNew.recursive(term[1]) do
            TNew.accumulate(out,after[1],term[2]*after[2]);
        od;
    od;
    return out;
end;;
TNew.prefix := function(columns)
    local out,i,v;
    out:=[];
    for i in [0..3] do
        for v in columns do Add(out,List(v,e->[Concatenation([i],e[1]),e[2]])); od;
    od;
    return out;
end;;
TNew.sym_insert := function(v)
    local out,term,word,c,r,after;
    out:=[];
    for term in v do
        word:=term[1]; c:=term[2]; TNew.accumulate(out,word,c);
        for r in [0..Length(word)-2] do
            after:=TNew.braid(word,r); word:=after[1]; c:=(-1)^after[2]*c;
            TNew.accumulate(out,word,c);
        od;
    od;
    return out;
end;;
fi;
