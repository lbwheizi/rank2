# Verify the finite group, normalized sign cocycle, induced action and braiding.
# GAP 4.16 core only.  See T_model_core.g for zero-based conventions.
if IsExistingFile("T_model_core.g") then Read("T_model_core.g");
else Read("certificates/T/T_model_core.g"); fi;
TLastCheckPassed:=false;;

TNew.verify_sl2_cocycle := function()
    local G,one,xs,mt,pt,n,i,j,k,l,g,h,a,b,c,generated,enlarged,word;
    G:=TNew.G; one:=TNew.one; xs:=TNew.xs; n:=Length(G);
    TNew.check(Length(Set(xs))=4,"four distinct support elements");
    TNew.check(Set(List(G,g->TNew.conj(g,xs[1])))=Set(xs),"support conjugacy class");
    TNew.check(Filtered(G,g->TNew.conj(g,xs[1])=xs[1])
              =Filtered(G,g->g[1]=0 and g[2]=0),"centralizer description");
    TNew.check(TNew.mul(TNew.mul(xs[1],xs[2]),xs[1])
              =TNew.mul(TNew.mul(xs[2],xs[1]),xs[2]),"T braid relation");
    TNew.check(TNew.mul(TNew.mul(xs[1],xs[1]),xs[1])
              =TNew.mul(TNew.mul(xs[2],xs[2]),xs[2]),"T equal cubes");
    TNew.check(TNew.mul(TNew.mul(xs[1],xs[2]),xs[3])=one,"central product");
    generated:=[one];
    repeat
        enlarged:=Set(Concatenation(generated,
                      Concatenation(List(generated,g->List(xs,x->TNew.mul(g,x))))));
        if enlarged=generated then break; fi;
        generated:=enlarged;
    until false;
    TNew.check(generated=G,"support generates the group");
    mt:=List(G,g->List(G,h->Position(G,TNew.mul(g,h))));
    pt:=List(G,g->List(G,h->List(G,k->TNew.phi(g,h,k))));
    for i in [1..n] do for j in [1..n] do for k in [1..n] do
        TNew.check(mt[mt[i][j]][k]=mt[i][mt[j][k]],"group associativity");
    od; od; od;
    for g in G do for h in G do
        TNew.check(TNew.phi(one,g,h)=0 and TNew.phi(g,one,h)=0
                  and TNew.phi(g,h,one)=0,"cocycle normalization");
    od; od;
    for i in [1..n] do for j in [1..n] do for k in [1..n] do for l in [1..n] do
        TNew.check((pt[j][k][l]+pt[i][mt[j][k]][l]+pt[i][j][k]
                  +pt[mt[i][j]][k][l]+pt[i][j][mt[k][l]]) mod 2=0,
                  "3-cocycle identity");
    od; od; od; od;
    Print("PASS: group associativity (24^3), normalization, cocycle identity (24^4)\n");
    for g in G do for h in G do for i in [0..3] do
        a:=TNew.act(h,i); b:=TNew.act(g,a[1]); c:=TNew.act(TNew.mul(g,h),i);
        TNew.check(b[1]=c[1] and (a[2]+b[2]) mod 2
                  =(TNew.alpha(xs[i+1],g,h)+c[2]) mod 2,"projective action");
    od; od; od;
    TNew.check(ForAll([1..4],i->TNew.beta[i][i]=1),"self-braiding -1");
    TNew.check((1+TNew.beta[2][3]+TNew.beta[4][2]+TNew.beta[3][4]) mod 2=0,
              "epsilon equals 1");
    TNew.check(TNew.word_action([0,1,2],[0,1,0,1,0,1])=[[0,1,2],1],
              "nontrivial full twist");
    for word in Cartesian([0..3],[0..3],[0..3]) do
        TNew.check(TNew.word_action(word,[0,1,0])=TNew.word_action(word,[1,0,1]),
                  "braid relation on W^3");
    od;
    Print("PASS: projective action (24^2*4), q=-1, epsilon=1\n");
    Print("PASS: nontrivial full twist and braid relation on W^3\n");
    TLastCheckPassed:=true;
    Print("PASS: verify_T_sl2_cocycle.g\n");
end;;
TNew.verify_sl2_cocycle();
