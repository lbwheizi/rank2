# Actual recursive adjoints and full quantum symmetrizers in degrees 2--4.
# Arithmetic is exact in Q(E(6)); no package or external program is used.
if IsExistingFile("T_model_core.g") then Read("T_model_core.g");
else Read("certificates/T/T_model_core.g"); fi;
TLastCheckPassed:=false;;

TNew.verify_actual_adjoint := function()
    local Y,m,columns,dims,u,v,minor,determinant,word,perms,perm,labels,
          bs,j,target,label,tmp,lifts,result,out,basis;
    Y:=List([0..3],i->[[[i],1]]); dims:=[];
    for m in [2..4] do
        columns:=List(TNew.prefix(Y),TNew.apply_recursive);
        Y:=TNew.independent(columns); Add(dims,Length(Y));
        TNew.check(Length(Y)=[4,2,0][m-1],Concatenation("recursive Y",String(m)," dimension"));
        if m=3 then
            u:=TNew.apply_recursive(List(TNew.recursive([1,2]),
                           e->[Concatenation([0],e[1]),e[2]]));
            v:=TNew.apply_recursive(List(TNew.recursive([0,3]),
                           e->[Concatenation([1],e[1]),e[2]]));
            minor:=List([[0,1,2],[3,0,2]],word->[TNew.get(u,word),TNew.get(v,word)]);
            determinant:=minor[1][1]*minor[2][2]-minor[1][2]*minor[2][1];
            TNew.check(minor=[[2-TNew.p,-TNew.p],[-1,-1]],"Y3 specified two-coordinate minor");
            TNew.check(determinant=-2,"Y3 minor determinant -2");
        fi;
    od;
    Print("PASS: Y2,Y3,Y4 dimensions = ",dims,"; specified Y3 minor det=-2\n");
    dims:=[];
    for m in [2..4] do
        lifts:=[];
        for perm in PermutationsList([0..m-1]) do
            labels:=[0..m-1]; bs:=[];
            for target in [1..m] do
                label:=perm[target]; j:=Position(labels,label);
                while j>target do
                    Add(bs,j-2); tmp:=labels[j]; labels[j]:=labels[j-1];
                    labels[j-1]:=tmp; j:=j-1;
                od;
            od;
            Add(lifts,bs);
        od;
        TNew.check(Length(lifts)=Factorial(m),"complete permutation list");
        columns:=[];
        for word in Cartesian(List([1..m],i->[0..3])) do
            result:=[];
            for bs in lifts do
                out:=TNew.word_action(word,bs);
                TNew.accumulate(result,out[1],(-1)^out[2]);
            od;
            Add(columns,result);
        od;
        TNew.check(Length(columns)=4^m,"all symmetrizer columns");
        basis:=TNew.independent(columns); Add(dims,Length(basis));
        TNew.check(Length(basis)=[8,10,8][m-1],"full symmetrizer rank");
    od;
    Print("PASS: full quantum symmetrizer ranks in degrees 2,3,4 = ",dims,"\n");
    TLastCheckPassed:=true;
    Print("PASS: verify_T_actual_adjoint.g\n");
end;;
TNew.verify_actual_adjoint();
