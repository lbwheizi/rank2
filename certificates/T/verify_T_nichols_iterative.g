# Iterated quantum symmetrizers through the first zero homogeneous component.
# Arithmetic is exact in Q(E(6)); no package or external program is used.
if IsExistingFile("T_model_core.g") then Read("T_model_core.g");
else Read("certificates/T/T_model_core.g"); fi;
TLastCheckPassed:=false;;

TNew.verify_nichols_iterative := function()
    local image,dims,m,columns;
    image:=List([0..3],i->[[[i],1]]); dims:=[1,4];
    for m in [2..7] do
        columns:=List(TNew.prefix(image),TNew.sym_insert);
        image:=TNew.independent(columns); Add(dims,Length(image));
        if IsEmpty(image) then break; fi;
    od;
    TNew.check(dims=[1,4,8,10,8,4,1,0],"Nichols homogeneous dimensions through degree 7");
    TNew.check(Sum(dims)=36,"Nichols total dimension 36");
    Print("PASS: Nichols Hilbert coefficients = ",dims,"\n");
    Print("PASS: Nichols total dimension = ",Sum(dims),"\n");
    TLastCheckPassed:=true;
    Print("PASS: verify_T_nichols_iterative.g\n");
end;;
TNew.verify_nichols_iterative();
