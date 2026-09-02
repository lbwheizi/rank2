# Auxiliary consistency check; not cited as a separate lemma in the paper.
# File: auxiliary_s_phi_R2_invariance.g
# This file is generated from an exact integer certificate.
# It uses only the GAP core library; no optional package is loaded.
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

CertAddLocal := function(vec,x,a,b,coefficient)
    CertAddPhi(vec,a,b,x,coefficient);
    CertAddPhi(vec,CertConjElt(CertMulElt(a,b),x),a,b,coefficient);
    CertAddPhi(vec,a,CertConjElt(b,x),b,-coefficient);
end;;

CertAddUpper := function(vec,x,a,b,coefficient)
    CertAddPhi(vec,x,a,b,coefficient);
    CertAddPhi(vec,CertConjElt(x,a),CertConjElt(x,b),x,coefficient);
    CertAddPhi(vec,CertConjElt(x,a),x,b,-coefficient);
end;;

CertAddCocycleRelator := function(vec,a,b,c,d,coefficient)
    CertAddPhi(vec,b,c,d,coefficient);
    CertAddPhi(vec,a,CertMulElt(b,c),d,coefficient);
    CertAddPhi(vec,a,b,c,coefficient);
    CertAddPhi(vec,CertMulElt(a,b),c,d,-coefficient);
    CertAddPhi(vec,a,b,CertMulElt(c,d),-coefficient);
end;;

CertDeltaCocVector := function(g,h)
    local out,eg,g2;
    out := [];
    eg := CertMulElt(eps,g);
    g2 := CertMulElt(g,g);
    CertAddLocal(out,h,eg,g,1);
    CertAddLocal(out,g,h,h,1);
    CertAddLocal(out,g,h,eps,-1);
    CertAddLocal(out,h,eps,g2,-1);
    return out;
end;;

# CertSPhiVector is the exponent vector of s_Phi in the abstract
# Gamma2 normal form.
CertSPhiVector := function(g,h)
    local out,e,eg,gh,egh,g2,q,eq;
    out := [];
    e := eps;
    eg := CertMulElt(e,g);
    gh := CertMulElt(g,h);
    egh := CertMulElt(e,gh);
    g2 := CertMulElt(g,g);
    q := CertMulElt(g2,h);
    eq := CertMulElt(e,q);

    # Twice the pure-cocycle part of kappa_3.
    CertAddPhi(out,g,eg,egh,2);
    CertAddLocal(out,g,g,h,2);
    CertAddLocal(out,g,h,g,-2);
    CertAddPhi(out,eg,g,egh,-2);
    CertAddLocal(out,g,e,g,-2);
    CertAddLocal(out,g,h,eg,-2);
    CertAddUpper(out,g,g,gh,-2);

    # B_h B_g/(A_h A_g), with sigma_2(h) and lambda cancelled.
    CertAddUpper(out,h,eg,eq,1);
    CertAddUpper(out,g,eg,eq,1);
    CertAddLocal(out,g,h,h,1);
    CertAddLocal(out,q,h,g,1);
    CertAddLocal(out,q,eg,h,-1);
    CertAddLocal(out,q,g,e,-1);
    CertAddUpper(out,h,g,q,-1);
    CertAddUpper(out,g,g,q,-1);

    # The induced-character substitutions for sigma_2(epsilon),
    # g acting on hv, and sigma_2(g^2).
    CertAddUpper(out,e,g,gh,1);
    CertAddUpper(out,e,g,h,1);
    CertAddLocal(out,g,g,h,1);
    CertAddLocal(out,g,e,g,-1);
    CertAddLocal(out,g,h,eg,-1);
    CertAddLocal(out,q,g,g,1);
    CertAddUpper(out,g2,g,gh,1);
    CertAddUpper(out,g2,g,h,1);
    CertAddLocal(out,g,g,g,-2);

    # The remaining factor in the displayed definition of s_Phi.
    CertAddLocal(out,g,e,e,2);
    CertAddVector(out,CertDeltaCocVector(g,h),-1);
    return out;
end;;

CertCopyVector := function(vec)
    return List(vec,pair -> [pair[1],pair[2]]);
end;;

CertL1Norm := function(vec)
    return Sum(vec,pair -> AbsInt(pair[2]));
end;;

hg := CertMulElt(hh,gg);;
hinv := CertInvElt(hh);;
target := CertSPhiVector(hg,hinv);;
CertAddVector(target,CertSPhiVector(gg,hh),-1);;
identityText := "s_Phi(hg,h^-1) = s_Phi(g,h)";;

CertificateTerms := [
    rec(coefficient := -1, quadruple := [[0, 0, -1], [0, 1, 1], [1, 1, 0], [1, 1, 1]]),
    rec(coefficient := 1, quadruple := [[0, 0, -1], [1, 1, 1], [1, 0, 0], [1, 1, 0]]),
    rec(coefficient := 1, quadruple := [[0, 0, 1], [1, 1, 0], [0, 1, 1], [0, 1, 0]]),
    rec(coefficient := 1, quadruple := [[0, 0, 1], [1, 1, 0], [1, 0, 0], [1, 1, 1]]),
    rec(coefficient := -1, quadruple := [[0, 1, 0], [0, 0, 1], [0, 1, 0], [1, 1, 1]]),
    rec(coefficient := 1, quadruple := [[0, 1, 0], [0, 1, 0], [0, 1, 0], [0, 0, 1]]),
    rec(coefficient := 1, quadruple := [[0, 1, 0], [0, 1, 0], [0, 1, 0], [0, 1, 1]]),
    rec(coefficient := 1, quadruple := [[0, 1, 0], [1, 0, 0], [0, 2, 1], [0, 1, 0]]),
    rec(coefficient := -1, quadruple := [[0, 1, 0], [1, 0, 0], [1, 1, 1], [0, 2, 0]]),
    rec(coefficient := -1, quadruple := [[0, 1, 0], [1, 0, 0], [1, 1, 1], [1, 1, 1]]),
    rec(coefficient := -3, quadruple := [[0, 1, 0], [1, 1, 1], [1, 0, 0], [1, 1, 1]]),
    rec(coefficient := -1, quadruple := [[0, 1, 0], [1, 1, 1], [1, 1, 0], [0, 0, 1]]),
    rec(coefficient := -4, quadruple := [[0, 1, 1], [0, 0, -1], [0, 1, 1], [1, 1, 1]]),
    rec(coefficient := 1, quadruple := [[0, 1, 1], [0, 0, -1], [1, 0, 0], [1, 2, 2]]),
    rec(coefficient := 3, quadruple := [[0, 1, 1], [0, 0, -1], [1, 1, 1], [0, 1, 1]]),
    rec(coefficient := -1, quadruple := [[0, 1, 1], [1, 0, 0], [0, 0, -1], [1, 2, 2]]),
    rec(coefficient := -3, quadruple := [[0, 1, 1], [1, 0, 0], [1, 1, 0], [0, 1, 1]]),
    rec(coefficient := -1, quadruple := [[1, 0, 0], [0, 0, 1], [0, 1, 0], [1, 0, 0]]),
    rec(coefficient := -1, quadruple := [[1, 0, 0], [0, 0, 1], [1, 1, 0], [1, 0, 0]]),
    rec(coefficient := -4, quadruple := [[1, 0, 0], [0, 1, 1], [1, 0, 0], [1, 2, 1]]),
    rec(coefficient := 4, quadruple := [[1, 0, 0], [0, 1, 1], [1, 1, 0], [1, 1, 1]]),
    rec(coefficient := 1, quadruple := [[1, 0, 0], [1, 0, 0], [0, 1, 1], [1, 0, 0]]),
    rec(coefficient := -1, quadruple := [[1, 0, 0], [1, 0, 0], [1, 1, 0], [0, 0, 1]]),
    rec(coefficient := 4, quadruple := [[1, 0, 0], [1, 0, 0], [1, 1, 0], [1, 0, 0]]),
    rec(coefficient := -3, quadruple := [[1, 0, 0], [1, 0, 0], [1, 1, 0], [1, 1, 1]]),
    rec(coefficient := 3, quadruple := [[1, 0, 0], [1, 0, 0], [1, 1, 1], [0, 1, 0]]),
    rec(coefficient := 1, quadruple := [[1, 0, 0], [1, 0, 0], [1, 1, 1], [1, 0, 0]]),
    rec(coefficient := 2, quadruple := [[1, 0, 0], [1, 1, 0], [0, 0, 1], [0, 1, 0]]),
    rec(coefficient := 1, quadruple := [[1, 0, 0], [1, 1, 0], [0, 0, 1], [1, 0, 0]]),
    rec(coefficient := 1, quadruple := [[1, 0, 0], [1, 1, 0], [0, 1, 1], [0, 1, 0]]),
    rec(coefficient := -2, quadruple := [[1, 0, 0], [1, 1, 0], [0, 1, 1], [1, 0, 0]]),
    rec(coefficient := -1, quadruple := [[1, 0, 0], [1, 1, 0], [1, 0, 0], [0, 0, 1]]),
    rec(coefficient := -1, quadruple := [[1, 0, 0], [1, 1, 0], [1, 0, 0], [0, 1, 0]]),
    rec(coefficient := -1, quadruple := [[1, 0, 0], [1, 1, 0], [1, 0, 0], [0, 1, 1]]),
    rec(coefficient := 2, quadruple := [[1, 0, 0], [1, 1, 0], [1, 0, 0], [1, 0, 0]]),
    rec(coefficient := -1, quadruple := [[1, 0, 0], [1, 1, 0], [1, 0, 0], [1, 1, 1]]),
    rec(coefficient := -1, quadruple := [[1, 0, 0], [1, 1, 0], [1, 1, 0], [0, 0, 1]]),
    rec(coefficient := -4, quadruple := [[1, 0, 0], [1, 1, 0], [1, 1, 1], [1, 0, 0]]),
    rec(coefficient := -1, quadruple := [[1, 0, 0], [1, 1, 1], [0, 0, -1], [0, 1, 1]]),
    rec(coefficient := 5, quadruple := [[1, 0, 0], [1, 1, 1], [0, 0, -1], [1, 1, 1]]),
    rec(coefficient := -1, quadruple := [[1, 0, 0], [1, 1, 1], [0, 1, 0], [0, 1, 0]]),
    rec(coefficient := 3, quadruple := [[1, 0, 0], [1, 1, 1], [0, 1, 0], [1, 0, 0]]),
    rec(coefficient := -4, quadruple := [[1, 0, 0], [1, 1, 1], [0, 1, 0], [1, 1, 1]]),
    rec(coefficient := -1, quadruple := [[1, 0, 0], [1, 1, 1], [1, 0, -1], [1, 1, 1]]),
    rec(coefficient := -1, quadruple := [[1, 0, 0], [1, 1, 1], [1, 0, 0], [0, 0, -1]]),
    rec(coefficient := -1, quadruple := [[1, 0, 0], [1, 1, 1], [1, 0, 0], [0, 1, 0]]),
    rec(coefficient := -4, quadruple := [[1, 0, 0], [1, 1, 1], [1, 0, 0], [1, 1, 0]]),
    rec(coefficient := 2, quadruple := [[1, 0, 0], [1, 1, 1], [1, 1, 0], [1, 0, 0]]),
    rec(coefficient := 1, quadruple := [[1, 0, 0], [1, 1, 1], [1, 1, 0], [1, 1, 1]]),
    rec(coefficient := 1, quadruple := [[1, 0, 0], [1, 1, 1], [1, 1, 1], [0, 0, -1]]),
    rec(coefficient := -1, quadruple := [[1, 0, 1], [0, 1, 0], [1, 0, 0], [1, 1, 0]]),
    rec(coefficient := -1, quadruple := [[1, 0, 1], [1, 1, 0], [1, 0, 0], [0, 1, 0]]),
    rec(coefficient := 2, quadruple := [[1, 1, 0], [0, 0, 1], [0, 1, 0], [1, 0, 0]]),
    rec(coefficient := -3, quadruple := [[1, 1, 0], [0, 0, 1], [0, 1, 0], [1, 1, 0]]),
    rec(coefficient := -1, quadruple := [[1, 1, 0], [0, 0, 1], [1, 0, 0], [0, 1, 0]]),
    rec(coefficient := -1, quadruple := [[1, 1, 0], [0, 0, 1], [1, 0, 0], [0, 2, 0]]),
    rec(coefficient := 4, quadruple := [[1, 1, 0], [0, 0, 1], [1, 1, 0], [0, 1, 0]]),
    rec(coefficient := -1, quadruple := [[1, 1, 0], [0, 1, 1], [0, 1, 0], [1, 0, 0]]),
    rec(coefficient := 1, quadruple := [[1, 1, 0], [1, 0, 0], [0, 0, 1], [0, 2, 0]]),
    rec(coefficient := 1, quadruple := [[1, 1, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]]),
    rec(coefficient := 1, quadruple := [[1, 1, 0], [1, 0, 0], [0, 1, 1], [0, 0, -1]]),
    rec(coefficient := -1, quadruple := [[1, 1, 0], [1, 0, 0], [0, 1, 1], [0, 1, 0]]),
    rec(coefficient := 3, quadruple := [[1, 1, 0], [1, 0, 0], [0, 1, 1], [1, 0, 0]]),
    rec(coefficient := 1, quadruple := [[1, 1, 0], [1, 0, 0], [0, 2, 1], [1, 0, 0]]),
    rec(coefficient := 1, quadruple := [[1, 1, 0], [1, 0, 0], [1, 0, 0], [0, 1, 0]]),
    rec(coefficient := -2, quadruple := [[1, 1, 0], [1, 0, 0], [1, 0, 0], [1, 1, 0]]),
    rec(coefficient := -1, quadruple := [[1, 1, 0], [1, 0, 0], [1, 0, 0], [1, 1, 1]]),
    rec(coefficient := -2, quadruple := [[1, 1, 0], [1, 0, 0], [1, 1, 0], [0, 0, 1]]),
    rec(coefficient := 2, quadruple := [[1, 1, 0], [1, 0, 0], [1, 1, 0], [1, 1, 1]]),
    rec(coefficient := 1, quadruple := [[1, 1, 0], [1, 0, 0], [1, 1, 0], [1, 2, 1]]),
    rec(coefficient := -2, quadruple := [[1, 1, 0], [1, 0, 0], [1, 1, 1], [0, 1, 0]]),
    rec(coefficient := 2, quadruple := [[1, 1, 0], [1, 0, 0], [1, 1, 1], [1, 0, 0]]),
    rec(coefficient := -1, quadruple := [[1, 1, 0], [1, 0, 0], [1, 1, 1], [1, 1, 0]]),
    rec(coefficient := -1, quadruple := [[1, 1, 0], [1, 0, 0], [1, 2, 1], [0, 1, 0]]),
    rec(coefficient := -1, quadruple := [[1, 1, 0], [1, 1, 1], [0, 0, -1], [1, 1, 1]]),
    rec(coefficient := 1, quadruple := [[1, 1, 0], [1, 1, 1], [0, 1, 0], [0, 1, 0]]),
    rec(coefficient := -3, quadruple := [[1, 1, 0], [1, 1, 1], [1, 0, 0], [0, 1, 0]]),
    rec(coefficient := 1, quadruple := [[1, 1, 1], [0, 0, -1], [1, 0, 0], [1, 1, 1]]),
    rec(coefficient := 1, quadruple := [[1, 1, 1], [0, 0, -1], [1, 1, 1], [0, 1, 0]]),
    rec(coefficient := -2, quadruple := [[1, 1, 1], [0, 0, -1], [1, 1, 1], [1, 0, 0]]),
    rec(coefficient := 1, quadruple := [[1, 1, 1], [0, 1, 0], [0, 1, 1], [0, 0, -1]]),
    rec(coefficient := 3, quadruple := [[1, 1, 1], [0, 1, 0], [1, 0, 0], [0, 1, 0]]),
    rec(coefficient := 1, quadruple := [[1, 1, 1], [0, 1, 0], [1, 1, 1], [1, 0, 0]]),
    rec(coefficient := 1, quadruple := [[1, 1, 1], [1, 0, 0], [0, 1, 0], [0, 1, 0]]),
    rec(coefficient := -4, quadruple := [[1, 1, 1], [1, 0, 0], [0, 1, 0], [1, 0, 0]]),
    rec(coefficient := 5, quadruple := [[1, 1, 1], [1, 0, 0], [0, 1, 0], [1, 1, 1]]),
    rec(coefficient := -1, quadruple := [[1, 1, 1], [1, 0, 0], [0, 1, 1], [0, 1, 0]]),
    rec(coefficient := 1, quadruple := [[1, 1, 1], [1, 0, 0], [1, 0, 0], [0, 1, 0]]),
    rec(coefficient := 2, quadruple := [[1, 1, 1], [1, 0, 0], [1, 0, 0], [1, 1, 0]]),
    rec(coefficient := -1, quadruple := [[1, 1, 1], [1, 0, 0], [1, 0, 0], [1, 1, 1]]),
    rec(coefficient := -2, quadruple := [[1, 1, 1], [1, 0, 0], [1, 1, 0], [1, 0, 0]]),
    rec(coefficient := -2, quadruple := [[1, 1, 1], [1, 0, 0], [1, 1, 0], [1, 1, 1]]),
    rec(coefficient := -3, quadruple := [[1, 1, 1], [1, 0, 0], [1, 1, 1], [0, 0, -1]]),
    rec(coefficient := 2, quadruple := [[1, 1, 1], [1, 0, 0], [1, 1, 1], [0, 1, 0]]),
    rec(coefficient := 1, quadruple := [[1, 1, 1], [1, 0, 0], [1, 1, 1], [0, 2, 1]]),
    rec(coefficient := -1, quadruple := [[1, 1, 1], [1, 0, 0], [1, 2, 1], [1, 1, 1]]),
    rec(coefficient := 1, quadruple := [[1, 1, 1], [1, 1, 0], [0, 0, 1], [0, 1, 0]]),
    rec(coefficient := 3, quadruple := [[1, 1, 1], [1, 1, 0], [1, 0, 0], [1, 1, 1]]),
    rec(coefficient := -1, quadruple := [[1, 1, 1], [1, 1, 0], [1, 1, 1], [1, 1, 1]]),
    rec(coefficient := -1, quadruple := [[1, 1, 1], [1, 1, 1], [1, 1, 1], [0, 0, -1]]),
    rec(coefficient := -1, quadruple := [[1, 1, 1], [1, 1, 1], [1, 1, 1], [1, 1, 0]])
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
Print("Identity: ",identityText,"\n");
Print("Target support: ",Length(target),"\n");
Print("Target L1: ",CertL1Norm(target),"\n");
Print("Certificate relators: ",Length(CertificateTerms),"\n");
Print("Certificate L1: ",Sum(CertificateTerms,
      term -> AbsInt(term.coefficient)),"\n");
Print("Exact residual support: ",Length(difference),"\n");

if Length(difference) = 0 then
    Print("PASS: exact integral residual is zero.\n");
else
    Print("FAIL: nonzero exact residual: ",difference,"\n");
    ErrorNoReturn("Certificate verification failed.");
fi;
