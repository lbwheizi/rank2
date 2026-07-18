# delta_certificate.g
#
# Formal verification of the cocycle reduction appearing in the
# simplicity criterion for the first adjoint object X_1.
#
# The script verifies
#
#       a^2 A_1 B_1
#   -------------------------  = Delta,
#   lambda' Phi^h(g,h) Phi^g(g,h)
#
# after cancellation of the one-dimensional projective-character
# parameters and substitution of
#
#       zeta^2 = Phi_g(epsilon,epsilon).
#
# It does not assert that Delta is identically 1.  The condition Delta=1
# comes from the proportionality condition used to characterize the
# simplicity of X_1.
#
# Group elements are represented in the normal form
#     epsilon^e * g^a * h^b
# by the GAP list [e,a,b], where e is reduced modulo 2.  Multiplication
# follows from
#     h*g = epsilon*g*h,  epsilon^2 = 1,  epsilon central.
#
# A formal product of values Phi(a,b,c) is stored as an exponent vector
# in the free abelian group on normalized triples [a|b|c].  Thus the
# calculation is independent of any particular normalized 3-cocycle.


#############################################################################
# The group Gamma_2 in normal form
#############################################################################

one := [ 0, 0, 0 ];
eps := [ 1, 0, 0 ];
gg  := [ 0, 1, 0 ];
hh  := [ 0, 0, 1 ];

MulElt := function(x,y)
    local e,a,b,f,c,d;
    e := x[1];  a := x[2];  b := x[3];
    f := y[1];  c := y[2];  d := y[3];
    return [ (e+f+b*c) mod 2, a+c, b+d ];
end;

InvElt := function(x)
    local e,a,b;
    e := x[1];  a := x[2];  b := x[3];
    return [ (-e+b*a) mod 2, -a, -b ];
end;

ConjElt := function(a,x)
    return MulElt(MulElt(a,x),InvElt(a));
end;


#############################################################################
# Formal exponent vectors for products of Phi(a,b,c)
#############################################################################

# A vector is a list of pairs [triple,coefficient], where
# triple = [a,b,c].  Terms containing the identity are omitted because
# Phi is normalized.

AddTerm := function(vec,triple,coefficient)
    local keys,pos;

    if coefficient = 0 then
        return;
    fi;

    if ForAny(triple,x -> x = one) then
        return;
    fi;

    keys := List(vec,pair -> pair[1]);
    pos := Position(keys,triple);

    if pos = fail then
        Add(vec,[triple,coefficient]);
    else
        vec[pos][2] := vec[pos][2] + coefficient;
    fi;
end;

AddPhi := function(vec,a,b,c,coefficient)
    AddTerm(vec,[a,b,c],coefficient);
end;

AddFormalVector := function(target,source,coefficient)
    local pair;
    for pair in source do
        AddTerm(target,pair[1],coefficient*pair[2]);
    od;
end;

NonzeroPart := function(vec)
    return Filtered(vec,pair -> pair[2] <> 0);
end;


#############################################################################
# Expansions of Phi_x(A,B), Phi^A(x,y), and the cocycle relator
#############################################################################

# Adds coefficient times the exponent vector of
#
# Phi_x(A,B)
#   = Phi(A,B,x) Phi(ABxB^-1A^-1,A,B)
#       / Phi(A,BxB^-1,B).

AddLocal := function(vec,x,A,B,coefficient)
    AddPhi(vec,A,B,x,coefficient);
    AddPhi(vec,ConjElt(MulElt(A,B),x),A,B,coefficient);
    AddPhi(vec,A,ConjElt(B,x),B,-coefficient);
end;

# Adds coefficient times the exponent vector of
#
# Phi^A(x,y)
#   = Phi(A,x,y) Phi(AxA^-1,AyA^-1,A)
#       / Phi(AxA^-1,A,y).

AddUpper := function(vec,A,x,y,coefficient)
    AddPhi(vec,A,x,y,coefficient);
    AddPhi(vec,ConjElt(A,x),ConjElt(A,y),A,coefficient);
    AddPhi(vec,ConjElt(A,x),A,y,-coefficient);
end;

# Adds coefficient times the exponent vector of the normalized
# 3-cocycle relator
#
# C(a,b,c,d)
#   = Phi(b,c,d) Phi(a,bc,d) Phi(a,b,c)
#       / (Phi(ab,c,d) Phi(a,b,cd)).
#
# Every such relator equals 1 for a normalized 3-cocycle.

AddCocycleRelator := function(vec,a,b,c,d,coefficient)
    AddPhi(vec,b,c,d,coefficient);
    AddPhi(vec,a,MulElt(b,c),d,coefficient);
    AddPhi(vec,a,b,c,coefficient);
    AddPhi(vec,MulElt(a,b),c,d,-coefficient);
    AddPhi(vec,a,b,MulElt(c,d),-coefficient);
end;


#############################################################################
# Formal expansion of the quotient of the coefficient by Delta
#############################################################################

eg := MulElt(eps,gg);
eh := MulElt(eps,hh);
g2 := MulElt(gg,gg);

quotient := [];

# After cancelling zeta, zeta', mu, mu', and lambda', the remaining
# character factor is zeta^2.  Substitute
# zeta^2 = Phi_g(epsilon,epsilon).
AddLocal(quotient,gg,eps,eps,1);

# a^2, where a = zeta/Phi_g(h,epsilon).
AddLocal(quotient,gg,hh,eps,-2);

# A_1 without its projective-character factors:
# Phi^h(eg,eh) Phi_g(h,h) Phi_h(h,g)
# / (Phi_h(epsilon,h) Phi_h(g,eh)).
AddUpper(quotient,hh,eg,eh,1);
AddLocal(quotient,gg,hh,hh,1);
AddLocal(quotient,hh,hh,gg,1);
AddLocal(quotient,hh,eps,hh,-1);
AddLocal(quotient,hh,gg,eh,-1);

# B_1 without its projective-character factors:
# Phi^g(eg,eh) Phi_g(g,h) Phi_h(g,g)
# / (Phi_g(epsilon,g) Phi_g(h,eg)).
AddUpper(quotient,gg,eg,eh,1);
AddLocal(quotient,gg,gg,hh,1);
AddLocal(quotient,hh,gg,gg,1);
AddLocal(quotient,gg,eps,gg,-1);
AddLocal(quotient,gg,hh,eg,-1);

# Divide by Phi^h(g,h) Phi^g(g,h).
AddUpper(quotient,hh,gg,hh,-1);
AddUpper(quotient,gg,gg,hh,-1);

# Divide by the cocycle part of
#
# Delta = zeta zeta' mu mu'
#         Phi_h(eg,g) Phi_g(h,h)
#         / (Phi_g(h,epsilon) Phi_h(epsilon,g^2)).
AddLocal(quotient,hh,eg,gg,-1);
AddLocal(quotient,gg,hh,hh,-1);
AddLocal(quotient,gg,hh,eps,1);
AddLocal(quotient,hh,eps,g2,1);


#############################################################################
# Product of cocycle relators claimed to equal the quotient
#############################################################################

certificate := [];

AddCocycleRelator(certificate,eg,  gg,  hh,  eps,  1);
AddCocycleRelator(certificate,eg,  eps, eps, gg,   1);
AddCocycleRelator(certificate,eg,  hh,  eg,  eps, -1);
AddCocycleRelator(certificate,gg,  eh,  eg,  eps,  1);
AddCocycleRelator(certificate,gg,  eg,  hh,  eps, -1);
AddCocycleRelator(certificate,gg,  hh,  eps, eg,  -1);
AddCocycleRelator(certificate,eps, eg,  eps, eps, -1);
AddCocycleRelator(certificate,eps, gg,  eh,  gg,  -1);
AddCocycleRelator(certificate,eps, gg,  gg,  hh,   1);
AddCocycleRelator(certificate,eps, hh,  gg,  gg,   1);
AddCocycleRelator(certificate,hh,  eg,  eps, eg,   1);
AddCocycleRelator(certificate,hh,  eps, eg,  eps,  1);
AddCocycleRelator(certificate,hh,  eps, gg,  eh,   1);
AddCocycleRelator(certificate,hh,  eps, gg,  gg,  -1);
AddCocycleRelator(certificate,hh,  eps, hh,  gg,  -1);


#############################################################################
# Verification
#############################################################################

difference := List(quotient,pair -> [pair[1],pair[2]]);
AddFormalVector(difference,certificate,-1);
difference := NonzeroPart(difference);

Print("Formal Phi terms in the reduced Delta quotient: ",
      Length(NonzeroPart(quotient)),"\n");
Print("Formal Phi terms remaining after subtracting the certificate: ",
      Length(difference),"\n");

if Length(difference) = 0 then
    Print("SUCCESS: the cocycle certificate for the Delta reduction is correct.\n");
    Print("Hence the first-adjoint coefficient equals Delta.\n");
else
    Print("FAILURE: the two formal exponent vectors are different.\n");
    Print("Nonzero difference: ",difference,"\n");
    Error("The Delta certificate did not verify");
fi;

quit;
