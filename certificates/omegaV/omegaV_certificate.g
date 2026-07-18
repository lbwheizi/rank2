# omegaV_certificate.g
#
# Formal verification of the cocycle certificate for Omega_V = 1.
#
# Group elements are represented in the normal form
#     epsilon^e * g^a * h^b
# by the GAP list [e,a,b], where e is reduced modulo 2.
# The multiplication rule follows from
#     h*g = epsilon*g*h,  epsilon^2 = 1,  epsilon central.
#
# A formal product of values Phi(a,b,c) is stored as an exponent
# vector in the free abelian group on normalized triples [a|b|c].
# Thus this script does not choose or enumerate a particular cocycle.


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

# A vector is a list of pairs [ triple, coefficient ], where
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
# Formal expansion of Omega_V
#############################################################################

eg  := MulElt(eps,gg);
eh  := MulElt(eps,hh);
gh  := MulElt(gg,hh);
hg  := MulElt(hh,gg);
egh := MulElt(eps,gh);

omega := [];

# zeta^2 = Phi_g(epsilon,epsilon)
AddLocal(omega,gg,eps,eps,1);

# Remaining factors in
#
# Omega_V = zeta^2 Phi(g,eg,eh) Phi_g(g,h) Phi_g(epsilon,h)
#           -------------------------------------------------- .
#           Phi(eg,g,eh) Phi^g(g,h) Phi_g(epsilon,hg)
#                          Phi_g(h,epsilon)^2

AddPhi(omega,gg,eg,eh,1);
AddLocal(omega,gg,gg,hh,1);
AddLocal(omega,gg,eps,hh,1);
AddPhi(omega,eg,gg,eh,-1);
AddUpper(omega,gg,gg,hh,-1);
AddLocal(omega,gg,eps,hg,-1);
AddLocal(omega,gg,hh,eps,-2);


#############################################################################
# Product of cocycle relators claimed to equal Omega_V
#############################################################################

certificate := [];

AddCocycleRelator(certificate,eg, eps,eg, eh, -1);
AddCocycleRelator(certificate,eg, eps,eg, hh,  1);
AddCocycleRelator(certificate,eg, eps,hh, eg,  1);
AddCocycleRelator(certificate,eg, hh, eps,eg, -1);
AddCocycleRelator(certificate,gg, hh, eps,gg, -1);
AddCocycleRelator(certificate,eps,eg, hh, gg, -1);
AddCocycleRelator(certificate,eps,eg, hh, eps,-1);
AddCocycleRelator(certificate,eps,hh, gg, eps, 1);
AddCocycleRelator(certificate,hh, gg, eps,eps, 1);
AddCocycleRelator(certificate,hh, eps,gg, eps,-1);
AddCocycleRelator(certificate,egh,eps,eps,gg,  1);


#############################################################################
# Verification
#############################################################################

difference := List(omega,pair -> [pair[1],pair[2]]);
AddFormalVector(difference,certificate,-1);
difference := NonzeroPart(difference);

Print("Formal Phi terms in the expanded Omega_V: ",
      Length(NonzeroPart(omega)),"\n");
Print("Formal Phi terms remaining after subtracting the certificate: ",
      Length(difference),"\n");

if Length(difference) = 0 then
    Print("SUCCESS: the cocycle certificate for Omega_V is correct.\n");
    Print("Hence Omega_V is a product of cocycle relators and equals 1.\n");
else
    Print("FAILURE: the two formal exponent vectors are different.\n");
    Print("Nonzero difference: ",difference,"\n");
    Error("The Omega_V certificate did not verify");
fi;

quit;
