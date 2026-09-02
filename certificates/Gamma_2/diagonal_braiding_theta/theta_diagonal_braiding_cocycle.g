# Paper correspondence: label lem:Theta-diagonal-braiding-cocycle.
# File: theta_diagonal_braiding_cocycle.g
# Exact integer certificate; GAP core only.
#
# An element [e,a,b] denotes epsilon^e g^a h^b in the abstract group
# h*g=epsilon*g*h, epsilon^2=1, epsilon central.  Only e is reduced.

one := [0,0,0];;
eps := [1,0,0];;
gg := [0,1,0];;
hh := [0,0,1];;

CertMulElt := function(x,y)
    return [ (x[1]+y[1]+x[3]*y[2]) mod 2,
             x[2]+y[2], x[3]+y[3] ];
end;;

CertInvElt := function(x)
    return [ (-x[1]+x[3]*x[2]) mod 2, -x[2], -x[3] ];
end;;

CertConjElt := function(a,x)
    return CertMulElt(CertMulElt(a,x),CertInvElt(a));
end;;

CertAddTerm := function(vec,triple,coefficient)
    local keys,pos,newcoefficient;
    if coefficient = 0 or ForAny(triple,x -> x = one) then
        return;
    fi;
    keys := List(vec,pair -> pair[1]);
    pos := Position(keys,triple);
    if pos = fail then
        Add(vec,[triple,coefficient]);
    else
        newcoefficient := vec[pos][2]+coefficient;
        if newcoefficient = 0 then
            Remove(vec,pos);
        else
            vec[pos][2] := newcoefficient;
        fi;
    fi;
end;;

CertAddPhi := function(vec,a,b,c,coefficient)
    CertAddTerm(vec,[a,b,c],coefficient);
end;;

CertAddVector := function(target,source,coefficient)
    local pair;
    for pair in source do
        CertAddTerm(target,pair[1],coefficient*pair[2]);
    od;
end;;

# Exponent vector of the local cocycle Phi_x(a,b).
CertAddLocal := function(vec,x,a,b,coefficient)
    CertAddPhi(vec,a,b,x,coefficient);
    CertAddPhi(vec,CertConjElt(CertMulElt(a,b),x),a,b,coefficient);
    CertAddPhi(vec,a,CertConjElt(b,x),b,-coefficient);
end;;

# Exponent vector of the tensor-action scalar Phi^x(a,b).
CertAddUpper := function(vec,x,a,b,coefficient)
    CertAddPhi(vec,x,a,b,coefficient);
    CertAddPhi(vec,CertConjElt(x,a),CertConjElt(x,b),x,coefficient);
    CertAddPhi(vec,CertConjElt(x,a),x,b,-coefficient);
end;;

# Exponent vector of the normalized 3-cocycle identity
# C_Phi(a,b,c,d).
CertAddCocycleRelator := function(vec,a,b,c,d,coefficient)
    CertAddPhi(vec,b,c,d,coefficient);
    CertAddPhi(vec,a,CertMulElt(b,c),d,coefficient);
    CertAddPhi(vec,a,b,c,coefficient);
    CertAddPhi(vec,CertMulElt(a,b),c,d,-coefficient);
    CertAddPhi(vec,a,b,CertMulElt(c,d),-coefficient);
end;;

# Exponent vector of Theta_Phi(g,h), expanded from its definition in
# eq:Theta-Phi-definition.
CertThetaVector := function(g,h)
    local out,eg,eh,gh,egh,g2;
    out := [];
    eg := CertMulElt(eps,g);
    eh := CertMulElt(eps,h);
    gh := CertMulElt(g,h);
    egh := CertMulElt(eps,gh);
    g2 := CertMulElt(g,g);

    CertAddLocal(out,h,eps,eps,1);
    CertAddLocal(out,h,eg,g,1);
    CertAddUpper(out,h,g,gh,1);
    CertAddUpper(out,h,eg,eh,1);
    CertAddLocal(out,g,h,g,2);

    CertAddLocal(out,g,h,eps,-1);
    CertAddLocal(out,h,eps,g2,-1);
    CertAddLocal(out,h,eps,h,-1);
    CertAddLocal(out,h,g,eh,-1);
    CertAddUpper(out,g,g,h,-1);

    CertAddLocal(out,h,h,g,1);
    CertAddLocal(out,g,eps,g,1);
    CertAddLocal(out,g,h,eg,1);

    CertAddUpper(out,h,eg,egh,-1);
    CertAddUpper(out,g,eg,eh,-1);
    CertAddUpper(out,h,eg,h,-1);
    CertAddLocal(out,g,g,h,-1);
    CertAddLocal(out,h,g,g,-1);
    return out;
end;;

CertCopyVector := function(vec)
    return List(vec,pair -> [pair[1],pair[2]]);
end;;

CertL1Norm := function(vec)
    return Sum(vec,pair -> AbsInt(pair[2]));
end;;

# Target:
# Theta_Phi(g,h) Phi_g(epsilon,epsilon) Phi_g(g,h)
# -------------------------------------------------
#       Phi_g(epsilon,g)^2 Phi_g(h,epsilon*g)
target := CertThetaVector(gg,hh);;
eg := CertMulElt(eps,gg);;
CertAddLocal(target,gg,eps,eps,1);;
CertAddLocal(target,gg,gg,hh,1);;
CertAddLocal(target,gg,eps,gg,-2);;
CertAddLocal(target,gg,hh,eg,-1);;

eh := CertMulElt(eps,hh);;
gh := CertMulElt(gg,hh);;

# The twenty factors displayed in the paper.  Each record represents
# C_Phi(a,b,c,d)^coefficient.
CertificateTerms := [
    rec(coefficient :=  1, quadruple := [eps,eps,eps,hh]),
    rec(coefficient := -1, quadruple := [eps,eps,eps,eg]),
    rec(coefficient := -1, quadruple := [eps,gg,gg,hh]),
    rec(coefficient :=  1, quadruple := [eps,gg,hh,eps]),
    rec(coefficient :=  1, quadruple := [eps,gg,eh,gg]),
    rec(coefficient := -1, quadruple := [eps,hh,eps,eg]),
    rec(coefficient := -1, quadruple := [eps,hh,gg,gg]),
    rec(coefficient := -1, quadruple := [eps,eh,gg,eps]),
    rec(coefficient :=  1, quadruple := [gg,eps,eps,eg]),
    rec(coefficient := -1, quadruple := [gg,hh,eps,gh]),
    rec(coefficient := -1, quadruple := [gg,eh,gg,hh]),
    rec(coefficient :=  1, quadruple := [hh,eps,eps,eh]),
    rec(coefficient := -1, quadruple := [hh,eps,gg,eps]),
    rec(coefficient :=  1, quadruple := [hh,eps,gg,gg]),
    rec(coefficient := -1, quadruple := [hh,eps,gg,hh]),
    rec(coefficient :=  1, quadruple := [hh,eg,eps,gh]),
    rec(coefficient :=  1, quadruple := [eg,eps,eg,eh]),
    rec(coefficient :=  1, quadruple := [eg,hh,gg,hh]),
    rec(coefficient := -1, quadruple := [eh,eps,eg,eh]),
    rec(coefficient :=  1, quadruple := [eh,eps,eh,gg])
];;

certificate := [];;
for term in CertificateTerms do
    CertAddCocycleRelator(certificate,
        term.quadruple[1],term.quadruple[2],
        term.quadruple[3],term.quadruple[4],term.coefficient);
od;

difference := CertCopyVector(target);;
CertAddVector(difference,certificate,-1);;

Print("GAP version: ",GAPInfo.Version,"\n");
Print("Identity: Theta diagonal-braiding cocycle reduction\n");
Print("Target support: ",Length(target),"\n");
Print("Target L1: ",CertL1Norm(target),"\n");
Print("Certificate relators: ",Length(CertificateTerms),"\n");
Print("Certificate coefficient L1: ",Sum(CertificateTerms,
      term -> AbsInt(term.coefficient)),"\n");
Print("Expanded certificate support: ",Length(certificate),"\n");
Print("Expanded certificate L1: ",CertL1Norm(certificate),"\n");
Print("Exact residual support: ",Length(difference),"\n");

if Length(difference) = 0 then
    Print("PASS: exact integral residual is zero.\n");
else
    Print("FAIL: nonzero exact residual: ",difference,"\n");
    ErrorNoReturn("Certificate verification failed.");
fi;
