# Paper correspondence:
#   lem:Theta-Phi-reduction
#   lem:opposite-Theta-reduction
#
# Exact Laurent-monomial verification of the two coefficient reductions
# used in the simplicity criteria for X_2 and Y_2 in the Gamma_2 case.
#
# The script constructs the two coefficient quotients from the definitions
# of A_1, B_1, D_1, p and of A_1', B_1', D_1', p'.  It constructs the two
# target Theta expressions separately from their displayed definitions.
# It then checks the unconditional formal identities
#
#   coefficient(g,h) / Theta_Phi(g,h)
#       = Delta^(-1) zeta'^2 / Phi_h(epsilon,epsilon),
#
#   coefficient(h,g) / Theta_Phi(h,g)
#       = Delta'^(-1) zeta^2 / Phi_g(epsilon,epsilon).
#
# Therefore Delta=1 and zeta'^2=Phi_h(epsilon,epsilon) give the first
# reduction, while Delta'=1 and zeta^2=Phi_g(epsilon,epsilon) give the
# opposite reduction.
#
# No cocycle value is sampled.  Every scalar, local cocycle Phi_x(a,b),
# and tensor-action scalar Phi^x(a,b) is an independent Laurent generator.


#############################################################################
# The group Gamma_2 in normal form
#############################################################################

one := [0,0,0];;
eps := [1,0,0];;
gg  := [0,1,0];;
hh  := [0,0,1];;

TCRMulElt := function(x,y)
    return [(x[1]+y[1]+x[3]*y[2]) mod 2,
            x[2]+y[2],x[3]+y[3]];
end;;

eg  := TCRMulElt(eps,gg);;
eh  := TCRMulElt(eps,hh);;
g2  := TCRMulElt(gg,gg);;
h2  := TCRMulElt(hh,hh);;
gh  := TCRMulElt(gg,hh);;
hg  := TCRMulElt(hh,gg);;
egh := TCRMulElt(eps,gh);;
ehg := TCRMulElt(eps,hg);;


#############################################################################
# Exact Laurent exponent vectors
#############################################################################

# A Laurent monomial is represented by a list [generator,exponent].
# Generators have one of the following forms:
#   ["scalar",name],
#   ["local",x,a,b] for Phi_x(a,b),
#   ["upper",x,a,b] for Phi^x(a,b).

AddTerm := function(vec,generator,coefficient)
    local generators,pos,newcoefficient;
    if coefficient = 0 then
        return;
    fi;
    generators := List(vec,pair -> pair[1]);
    pos := Position(generators,generator);
    if pos = fail then
        Add(vec,[generator,coefficient]);
    else
        newcoefficient := vec[pos][2]+coefficient;
        if newcoefficient = 0 then
            Remove(vec,pos);
        else
            vec[pos][2] := newcoefficient;
        fi;
    fi;
end;;

AddScalar := function(vec,name,coefficient)
    AddTerm(vec,["scalar",name],coefficient);
end;;

AddLocal := function(vec,x,a,b,coefficient)
    AddTerm(vec,["local",x,a,b],coefficient);
end;;

AddUpper := function(vec,x,a,b,coefficient)
    AddTerm(vec,["upper",x,a,b],coefficient);
end;;

TCRAddVector := function(target,source,coefficient)
    local pair;
    for pair in source do
        AddTerm(target,pair[1],coefficient*pair[2]);
    od;
end;;

CopyVector := function(vec)
    return List(vec,pair -> [pair[1],pair[2]]);
end;;

L1Norm := function(vec)
    return Sum(vec,pair -> AbsInt(pair[2]));
end;;


#############################################################################
# Forward definitions: (g,h)
#############################################################################

ForwardA1 := function()
    local out;
    out := [];
    AddUpper(out,hh,eg,eh,1);
    AddLocal(out,gg,hh,hh,1);
    AddScalar(out,"mu",1);
    AddLocal(out,hh,hh,gg,1);
    AddScalar(out,"lambdaPrime",1);
    AddScalar(out,"zetaPrime",1);
    AddLocal(out,hh,eps,hh,-1);
    AddLocal(out,hh,gg,eh,-1);
    return out;
end;;

ForwardB1 := function()
    local out;
    out := [];
    AddUpper(out,gg,eg,eh,1);
    AddLocal(out,gg,gg,hh,1);
    AddScalar(out,"zeta",1);
    AddLocal(out,gg,eps,gg,-1);
    AddLocal(out,gg,hh,eg,-1);
    AddLocal(out,hh,gg,gg,1);
    AddScalar(out,"muPrime",1);
    return out;
end;;

ForwardD1 := function()
    local out;
    out := [];
    AddUpper(out,hh,eg,hh,1);
    AddLocal(out,gg,hh,hh,1);
    AddScalar(out,"mu",1);
    AddScalar(out,"lambdaPrime",1);
    return out;
end;;

ForwardP := function()
    local out;
    out := [];
    AddLocal(out,gg,hh,gg,-1);
    return out;
end;;

ForwardCoefficient := function()
    local out;
    out := [];
    AddUpper(out,hh,gg,gh,1);
    TCRAddVector(out,ForwardA1(),1);
    TCRAddVector(out,ForwardP(),-2);
    AddUpper(out,gg,gg,hh,-1);
    AddUpper(out,hh,eg,egh,-1);
    AddLocal(out,gg,hh,hh,-1);
    AddScalar(out,"mu",-1);
    TCRAddVector(out,ForwardB1(),-1);
    TCRAddVector(out,ForwardD1(),-1);
    return out;
end;;

ForwardTheta := function()
    local out;
    out := [];
    AddLocal(out,hh,eps,eps,1);
    AddLocal(out,hh,eg,gg,1);
    AddUpper(out,hh,gg,gh,1);
    AddUpper(out,hh,eg,eh,1);
    AddLocal(out,gg,hh,gg,2);
    AddLocal(out,gg,hh,eps,-1);
    AddLocal(out,hh,eps,g2,-1);
    AddLocal(out,hh,eps,hh,-1);
    AddLocal(out,hh,gg,eh,-1);
    AddUpper(out,gg,gg,hh,-1);
    AddLocal(out,hh,hh,gg,1);
    AddLocal(out,gg,eps,gg,1);
    AddLocal(out,gg,hh,eg,1);
    AddUpper(out,hh,eg,egh,-1);
    AddUpper(out,gg,eg,eh,-1);
    AddUpper(out,hh,eg,hh,-1);
    AddLocal(out,gg,gg,hh,-1);
    AddLocal(out,hh,gg,gg,-1);
    return out;
end;;

ForwardDelta := function()
    local out;
    out := [];
    AddScalar(out,"zeta",1);
    AddScalar(out,"zetaPrime",1);
    AddScalar(out,"mu",1);
    AddScalar(out,"muPrime",1);
    AddLocal(out,hh,eg,gg,1);
    AddLocal(out,gg,hh,hh,1);
    AddLocal(out,gg,hh,eps,-1);
    AddLocal(out,hh,eps,g2,-1);
    return out;
end;;

ZetaPrimeSquareRelation := function()
    local out;
    out := [];
    AddScalar(out,"zetaPrime",2);
    AddLocal(out,hh,eps,eps,-1);
    return out;
end;;


#############################################################################
# Opposite definitions: (h,g)
#############################################################################

OppositeA1 := function()
    local out;
    out := [];
    AddUpper(out,gg,eh,eg,1);
    AddLocal(out,hh,gg,gg,1);
    AddScalar(out,"muPrime",1);
    AddLocal(out,gg,gg,hh,1);
    AddScalar(out,"lambda",1);
    AddScalar(out,"zeta",1);
    AddLocal(out,gg,eps,gg,-1);
    AddLocal(out,gg,hh,eg,-1);
    return out;
end;;

OppositeB1 := function()
    local out;
    out := [];
    AddUpper(out,hh,eh,eg,1);
    AddLocal(out,hh,hh,gg,1);
    AddScalar(out,"zetaPrime",1);
    AddLocal(out,hh,eps,hh,-1);
    AddLocal(out,hh,gg,eh,-1);
    AddLocal(out,gg,hh,hh,1);
    AddScalar(out,"mu",1);
    return out;
end;;

OppositeD1 := function()
    local out;
    out := [];
    AddUpper(out,gg,eh,gg,1);
    AddLocal(out,hh,gg,gg,1);
    AddScalar(out,"muPrime",1);
    AddScalar(out,"lambda",1);
    return out;
end;;

OppositeP := function()
    local out;
    out := [];
    AddLocal(out,hh,gg,hh,-1);
    return out;
end;;

OppositeCoefficient := function()
    local out;
    out := [];
    AddUpper(out,gg,hh,hg,1);
    TCRAddVector(out,OppositeA1(),1);
    TCRAddVector(out,OppositeP(),-2);
    AddUpper(out,hh,hh,gg,-1);
    AddUpper(out,gg,eh,ehg,-1);
    AddLocal(out,hh,gg,gg,-1);
    AddScalar(out,"muPrime",-1);
    TCRAddVector(out,OppositeB1(),-1);
    TCRAddVector(out,OppositeD1(),-1);
    return out;
end;;

OppositeTheta := function()
    local out;
    out := [];
    AddLocal(out,gg,eps,eps,1);
    AddLocal(out,gg,eh,hh,1);
    AddUpper(out,gg,hh,hg,1);
    AddUpper(out,gg,eh,eg,1);
    AddLocal(out,hh,gg,hh,2);
    AddLocal(out,hh,gg,eps,-1);
    AddLocal(out,gg,eps,h2,-1);
    AddLocal(out,gg,eps,gg,-1);
    AddLocal(out,gg,hh,eg,-1);
    AddUpper(out,hh,hh,gg,-1);
    AddLocal(out,gg,gg,hh,1);
    AddLocal(out,hh,eps,hh,1);
    AddLocal(out,hh,gg,eh,1);
    AddUpper(out,gg,eh,ehg,-1);
    AddUpper(out,hh,eh,eg,-1);
    AddUpper(out,gg,eh,gg,-1);
    AddLocal(out,hh,hh,gg,-1);
    AddLocal(out,gg,hh,hh,-1);
    return out;
end;;

OppositeDelta := function()
    local out;
    out := [];
    AddScalar(out,"zeta",1);
    AddScalar(out,"zetaPrime",1);
    AddScalar(out,"mu",1);
    AddScalar(out,"muPrime",1);
    AddLocal(out,gg,eh,hh,1);
    AddLocal(out,hh,gg,gg,1);
    AddLocal(out,hh,gg,eps,-1);
    AddLocal(out,gg,eps,h2,-1);
    return out;
end;;

ZetaSquareRelation := function()
    local out;
    out := [];
    AddScalar(out,"zeta",2);
    AddLocal(out,gg,eps,eps,-1);
    return out;
end;;


#############################################################################
# Verification
#############################################################################

if hg <> egh or gh <> ehg then
    ErrorNoReturn("The Gamma_2 normal-form multiplication failed.");
fi;

forwardQuotient := ForwardCoefficient();;
TCRAddVector(forwardQuotient,ForwardTheta(),-1);;
forwardExpected := [];;
TCRAddVector(forwardExpected,ForwardDelta(),-1);;
TCRAddVector(forwardExpected,ZetaPrimeSquareRelation(),1);;
forwardResidual := CopyVector(forwardQuotient);;
TCRAddVector(forwardResidual,forwardExpected,-1);;

oppositeQuotient := OppositeCoefficient();;
TCRAddVector(oppositeQuotient,OppositeTheta(),-1);;
oppositeExpected := [];;
TCRAddVector(oppositeExpected,OppositeDelta(),-1);;
TCRAddVector(oppositeExpected,ZetaSquareRelation(),1);;
oppositeResidual := CopyVector(oppositeQuotient);;
TCRAddVector(oppositeResidual,oppositeExpected,-1);;

# Non-vacuity controls: change one exponent in each source coefficient and
# rerun the same quotient and residual comparisons from that perturbed source.
forwardPerturbedQuotient := ForwardCoefficient();;
AddLocal(forwardPerturbedQuotient,gg,hh,gg,1);;
TCRAddVector(forwardPerturbedQuotient,ForwardTheta(),-1);;
forwardPerturbedResidual := CopyVector(forwardPerturbedQuotient);;
TCRAddVector(forwardPerturbedResidual,forwardExpected,-1);;

oppositePerturbedQuotient := OppositeCoefficient();;
AddLocal(oppositePerturbedQuotient,hh,gg,hh,1);;
TCRAddVector(oppositePerturbedQuotient,OppositeTheta(),-1);;
oppositePerturbedResidual := CopyVector(oppositePerturbedQuotient);;
TCRAddVector(oppositePerturbedResidual,oppositeExpected,-1);;

Print("GAP version: ",GAPInfo.Version,"\n");
Print("Arithmetic: exact integral Laurent exponent vectors\n");
Print("Forward quotient support: ",Length(forwardQuotient),"\n");
Print("Forward quotient L1: ",L1Norm(forwardQuotient),"\n");
Print("Forward exact residual support: ",Length(forwardResidual),"\n");
Print("Opposite quotient support: ",Length(oppositeQuotient),"\n");
Print("Opposite quotient L1: ",L1Norm(oppositeQuotient),"\n");
Print("Opposite exact residual support: ",Length(oppositeResidual),"\n");
Print("Forward perturbed residual support: ",
      Length(forwardPerturbedResidual),"\n");
Print("Opposite perturbed residual support: ",
      Length(oppositePerturbedResidual),"\n");

if Length(forwardResidual) <> 0 then
    Print("FAIL: forward nonzero exact residual: ",forwardResidual,"\n");
    ErrorNoReturn("Forward coefficient reduction failed.");
fi;

if Length(oppositeResidual) <> 0 then
    Print("FAIL: opposite nonzero exact residual: ",oppositeResidual,"\n");
    ErrorNoReturn("Opposite coefficient reduction failed.");
fi;

if Length(forwardPerturbedResidual) = 0
   or Length(oppositePerturbedResidual) = 0 then
    ErrorNoReturn("A non-vacuity control was not detected.");
fi;

Print("PASS: both exact Laurent identities hold.\n");
Print("PASS: both deliberate one-factor perturbations are detected.\n");
Print("Hence both conditional coefficient reductions are verified.\n");
