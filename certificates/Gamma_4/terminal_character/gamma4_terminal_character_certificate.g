#############################################################################
## Gamma4 terminal-root self-character certificate (core GAP only)
##
## PURPOSE
##   Verify exactly the Laurent identity arising in the terminal-root
##   self-character calculation:
##
##       Q = product_{(a,b,c,d),n in CERTIFICATE}
##               D(a,b,c,d)^n,
##
##   with 146 normalized 3-cocycle defects and coefficient l1 norm 245.
##   Here
##
##       D(a,b,c,d) = Phi(b,c,d) Phi(a,bc,d) Phi(a,b,c)
##                    / (Phi(ab,c,d) Phi(a,b,cd)).
##
##   Thus Q=1 for every normalized 3-cocycle.  The vector Q is the exact
##   target remaining in the terminal-character calculation after the
##   projective-character substitutions rho(h)=sigma(g)=-1, Delta4=1,
##   and epsilon^4=1.
##
## SCOPE
##   No external GAP packages, floating-point arithmetic, finite quotient,
##   or sampled cocycle is used.  Group elements are abstract normal forms
##   epsilon^i h^j g^k with i modulo 4 and j,k integral.  The computation is
##   exact in the free Laurent exponent lattice on formal Phi(a,b,c).
##
## REPRODUCTION
##   gap -q gamma4_terminal_character_certificate.g
#############################################################################

AddExponent := function(v, key, exponent)
    local i;
    if exponent = 0 then
        return;
    fi;
    for i in [1 .. Length(v)] do
        if v[i][1] = key then
            v[i][2] := v[i][2] + exponent;
            if v[i][2] = 0 then
                Remove(v, i);
            fi;
            return;
        fi;
    od;
    Add(v, [key, exponent]);
end;

CanonicalSparse := function(v)
    local ans, term;
    ans := [];
    for term in v do
        AddExponent(ans, term[1], term[2]);
    od;
    Sort(ans, function(left, right)
        return left[1] < right[1];
    end);
    return ans;
end;

AddSparse := function(left, right)
    local ans, term;
    ans := List(left, ShallowCopy);
    for term in right do
        AddExponent(ans, term[1], term[2]);
    od;
    return CanonicalSparse(ans);
end;

ScaleSparse := function(v, scalar)
    return CanonicalSparse(List(v, term -> [term[1], scalar * term[2]]));
end;

L1Sparse := function(v)
    return Sum(v, term -> AbsInt(term[2]));
end;

IdentityElt := [ 0, 0, 0 ];

Multiply := function(a, b)
    local i, j, k, ii, jj, kk, first;
    i := a[1]; j := a[2]; k := a[3];
    ii := b[1]; jj := b[2]; kk := b[3];
    if (k mod 2) = 0 then
        first := i + ii;
    else
        first := i - ii - jj;
    fi;
    return [first mod 4, j + jj, k + kk];
end;

PhiAtom := function(a, b, c)
    local key;
    if a = IdentityElt or b = IdentityElt or c = IdentityElt then
        return [];
    fi;
    key := Concatenation("P(", String(a), ",", String(b), ",", String(c), ")");
    return [[key, 1]];
end;

CocycleDefect := function(a, b, c, d)
    local ans;
    ans := [];
    ans := AddSparse(ans, PhiAtom(b, c, d));
    ans := AddSparse(ans, PhiAtom(a, Multiply(b, c), d));
    ans := AddSparse(ans, PhiAtom(a, b, c));
    ans := AddSparse(ans, ScaleSparse(PhiAtom(Multiply(a, b), c, d), -1));
    ans := AddSparse(ans, ScaleSparse(PhiAtom(a, b, Multiply(c, d)), -1));
    return ans;
end;

## Exact 76-atom Laurent target from the terminal-character reduction.
Target := [];
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 0, 0, 1 ], [ 0, 0, -1 ], [ 0, 0, 1 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 0, 0, 1 ], [ 0, 1, 0 ], [ 2, 0, 1 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 0, 0, 1 ], [ 0, 1, 0 ], [ 3, 0, 0 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 0, 0, 1 ], [ 1, 2, 1 ], [ 0, 0, 1 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 0, 0, 1 ], [ 1, 2, 2 ], [ 0, 0, 1 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 0, 0, 1 ], [ 2, 0, 1 ], [ 0, 0, 1 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 0, 0, 1 ], [ 2, 0, 1 ], [ 3, 1, 0 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 0, 0, 1 ], [ 3, 0, 0 ], [ 0, 1, 0 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 0, 0, 1 ], [ 3, 1, 0 ], [ 0, 1, 0 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 0, 1, 0 ], [ 0, 0, 1 ], [ 0, 1, 0 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 0, 1, 0 ], [ 0, 0, 1 ], [ 2, 0, 2 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 0, 1, 0 ], [ 0, 1, 0 ], [ 0, 0, 1 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 0, 1, 0 ], [ 1, 0, 0 ], [ 0, 0, 1 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 0, 1, 0 ], [ 1, 0, 0 ], [ 0, 1, 0 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 0, 1, 0 ], [ 1, 0, 0 ], [ 1, 0, 0 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 0, 1, 0 ], [ 1, 1, 0 ], [ 0, 0, 2 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 0, 1, 0 ], [ 2, 0, 0 ], [ 0, 0, 2 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 0, 1, 0 ], [ 2, 0, 1 ], [ 0, 0, 1 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 0, 1, 0 ], [ 2, 0, 1 ], [ 0, 1, 1 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 0, 1, 0 ], [ 2, 0, 1 ], [ 1, 0, 0 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 0, 1, 0 ], [ 2, 0, 2 ], [ 0, 0, 1 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 0, 1, 0 ], [ 3, 0, 0 ], [ 1, 0, 0 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 0, 1, 0 ], [ 3, 1, 1 ], [ 2, 0, 1 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 1, 0, 0 ], [ 0, 0, 1 ], [ 1, 2, 1 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 1, 0, 0 ], [ 0, 0, 1 ], [ 1, 2, 3 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 1, 0, 0 ], [ 0, 0, 1 ], [ 3, 2, 0 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 1, 0, 0 ], [ 0, 1, 0 ], [ 1, 0, 0 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 1, 0, 0 ], [ 1, 0, 0 ], [ 0, 1, 0 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 1, 0, 0 ], [ 1, 2, 1 ], [ 0, 0, 1 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 1, 0, 0 ], [ 1, 2, 3 ], [ 0, 0, 1 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 1, 0, 0 ], [ 3, 2, 0 ], [ 0, 0, 1 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 1, 0, 1 ], [ 0, 1, 0 ], [ 1, 1, 2 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 1, 0, 1 ], [ 0, 1, 0 ], [ 2, 0, 2 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 1, 0, 1 ], [ 1, 1, 2 ], [ 0, 1, 0 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 1, 0, 1 ], [ 1, 1, 2 ], [ 1, 0, 0 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 1, 1, 0 ], [ 0, 0, 2 ], [ 0, 1, 0 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 1, 1, 0 ], [ 0, 1, 0 ], [ 0, 0, 2 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 1, 1, 1 ], [ 0, 0, -1 ], [ 0, 0, 1 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 1, 1, 1 ], [ 0, 0, 1 ], [ 0, 0, -1 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 1, 1, 2 ], [ 1, 0, 0 ], [ 0, 0, 1 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 1, 1, 2 ], [ 2, 0, 1 ], [ 0, 1, 0 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 1, 1, 2 ], [ 2, 0, 1 ], [ 1, 0, 0 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 1, 1, 2 ], [ 2, 0, 1 ], [ 1, 0, 1 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 1, 1, 2 ], [ 3, 0, 1 ], [ 1, 1, 1 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 1, 1, 2 ], [ 3, 0, 1 ], [ 2, 0, 1 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 2, 0, 0 ], [ 0, 0, 2 ], [ 0, 1, 0 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 2, 0, 0 ], [ 0, 1, 0 ], [ 0, 0, 2 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 2, 0, 1 ], [ 0, 0, 1 ], [ 0, 1, 0 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 2, 0, 1 ], [ 0, 0, 1 ], [ 3, 1, 0 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 2, 0, 1 ], [ 0, 1, 0 ], [ 0, 0, 1 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 2, 0, 1 ], [ 0, 1, 0 ], [ 0, 1, 0 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 2, 0, 1 ], [ 0, 1, 1 ], [ 0, 1, 0 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 2, 0, 1 ], [ 0, 1, 1 ], [ 1, 1, 2 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 2, 0, 1 ], [ 1, 0, 0 ], [ 1, 2, 1 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 2, 0, 1 ], [ 1, 0, 0 ], [ 1, 2, 3 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 2, 0, 1 ], [ 1, 0, 0 ], [ 3, 2, 0 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 2, 0, 1 ], [ 1, 0, 1 ], [ 0, 1, 0 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 2, 0, 1 ], [ 1, 1, 2 ], [ 1, 1, 1 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 2, 0, 1 ], [ 1, 1, 2 ], [ 2, 0, 1 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 2, 0, 1 ], [ 3, 1, 0 ], [ 0, 0, 1 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 2, 0, 1 ], [ 3, 1, 3 ], [ 0, 1, 0 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 2, 0, 1 ], [ 3, 2, 1 ], [ 1, 0, 0 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 3, 0, 0 ], [ 0, 1, 0 ], [ 1, 0, 0 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 3, 0, 0 ], [ 1, 0, 0 ], [ 0, 1, 0 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 3, 0, 1 ], [ 0, 1, 0 ], [ 0, 1, 1 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 3, 0, 1 ], [ 0, 1, 0 ], [ 1, 0, 0 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 3, 0, 1 ], [ 1, 1, 1 ], [ 0, 0, -1 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 3, 0, 1 ], [ 1, 1, 1 ], [ 0, 1, 0 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 3, 0, 1 ], [ 2, 0, 1 ], [ 0, 1, 0 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 3, 1, 0 ], [ 0, 0, 1 ], [ 3, 0, 0 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 3, 1, 0 ], [ 0, 1, 0 ], [ 0, 0, 1 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 3, 1, 0 ], [ 1, 0, 1 ], [ 0, 1, 0 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 3, 1, 3 ], [ 0, 1, 0 ], [ 0, 0, 1 ]), 1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 3, 1, 3 ], [ 1, 0, 1 ], [ 0, 1, 0 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 3, 2, 1 ], [ 1, 0, 0 ], [ 0, 0, 1 ]), -1));
Target := AddSparse(Target, ScaleSparse(PhiAtom([ 3, 2, 1 ], [ 2, 0, 1 ], [ 1, 0, 0 ]), 1));

## Exact 146-term cocycle-defect certificate.
Certificate := [
    [ [ [ 0, 0, 1 ], [ 0, 1, 0 ], [ 3, 0, 0 ], [ 0, 1, 0 ] ], -2 ],
    [ [ [ 0, 0, 1 ], [ 0, 1, 0 ], [ 3, 0, 0 ], [ 1, 0, 0 ] ], 2 ],
    [ [ [ 0, 0, 1 ], [ 3, 0, 1 ], [ 0, 1, 0 ], [ 1, 0, 1 ] ], 4 ],
    [ [ [ 0, 0, 1 ], [ 3, 0, 1 ], [ 1, 0, 0 ], [ 0, 0, 1 ] ], -4 ],
    [ [ [ 0, 0, 1 ], [ 3, 0, 1 ], [ 1, 0, 0 ], [ 3, 1, 0 ] ], 4 ],
    [ [ [ 0, 0, 1 ], [ 3, 0, 1 ], [ 1, 0, 1 ], [ 3, 1, 0 ] ], -4 ],
    [ [ [ 0, 0, 1 ], [ 3, 1, 0 ], [ 1, 0, 0 ], [ 1, 0, 1 ] ], -1 ],
    [ [ [ 0, 0, 1 ], [ 3, 1, 0 ], [ 1, 0, 0 ], [ 1, 1, 2 ] ], -1 ],
    [ [ [ 0, 0, 1 ], [ 3, 1, 0 ], [ 2, 0, 2 ], [ 0, 1, 0 ] ], 1 ],
    [ [ [ 0, 0, 2 ], [ 0, 1, 0 ], [ 1, 0, 0 ], [ 1, 0, 0 ] ], 1 ],
    [ [ [ 0, 0, 2 ], [ 1, 0, 0 ], [ 1, 0, 0 ], [ 3, 1, 0 ] ], 1 ],
    [ [ [ 0, 0, 2 ], [ 1, 0, 0 ], [ 3, 1, 0 ], [ 1, 0, 0 ] ], -2 ],
    [ [ [ 0, 0, 2 ], [ 2, 0, 0 ], [ 3, 1, 0 ], [ 1, 0, 0 ] ], -1 ],
    [ [ [ 0, 1, 0 ], [ 0, 0, 1 ], [ 0, 1, 0 ], [ 1, 0, 1 ] ], 1 ],
    [ [ [ 0, 1, 0 ], [ 0, 0, 1 ], [ 0, 1, 0 ], [ 2, 0, 1 ] ], -1 ],
    [ [ [ 0, 1, 0 ], [ 0, 0, 1 ], [ 2, 0, 1 ], [ 0, 0, 1 ] ], 3 ],
    [ [ [ 0, 1, 0 ], [ 0, 0, 1 ], [ 2, 0, 1 ], [ 3, 1, 0 ] ], -3 ],
    [ [ [ 0, 1, 0 ], [ 0, 0, 1 ], [ 2, 0, 2 ], [ 3, 1, 0 ] ], 4 ],
    [ [ [ 0, 1, 0 ], [ 0, 0, 1 ], [ 2, 1, 1 ], [ 1, 0, 1 ] ], -4 ],
    [ [ [ 0, 1, 0 ], [ 0, 0, 2 ], [ 1, 0, 0 ], [ 1, 0, 0 ] ], -1 ],
    [ [ [ 0, 1, 0 ], [ 1, 0, 0 ], [ 0, 0, 2 ], [ 0, 1, 0 ] ], 2 ],
    [ [ [ 0, 1, 0 ], [ 1, 0, 0 ], [ 0, 1, 0 ], [ 0, 0, 2 ] ], -4 ],
    [ [ [ 0, 1, 0 ], [ 1, 0, 0 ], [ 1, 0, 0 ], [ 0, 0, 2 ] ], 2 ],
    [ [ [ 0, 1, 0 ], [ 1, 0, 0 ], [ 1, 0, 2 ], [ 3, 1, 0 ] ], 2 ],
    [ [ [ 0, 1, 0 ], [ 1, 0, 2 ], [ 0, 1, 0 ], [ 1, 0, 1 ] ], 4 ],
    [ [ [ 0, 1, 0 ], [ 1, 0, 2 ], [ 1, 0, 0 ], [ 0, 0, 1 ] ], -4 ],
    [ [ [ 0, 1, 0 ], [ 1, 0, 2 ], [ 1, 0, 0 ], [ 3, 1, 0 ] ], 5 ],
    [ [ [ 0, 1, 0 ], [ 1, 0, 2 ], [ 1, 0, 1 ], [ 3, 1, 0 ] ], -4 ],
    [ [ [ 0, 1, 0 ], [ 2, 0, 0 ], [ 0, 0, 2 ], [ 2, 0, 0 ] ], 1 ],
    [ [ [ 0, 1, 0 ], [ 2, 0, 0 ], [ 1, 0, 0 ], [ 1, 0, 0 ] ], -1 ],
    [ [ [ 0, 1, 0 ], [ 2, 0, 0 ], [ 2, 0, 0 ], [ 2, 0, 2 ] ], -1 ],
    [ [ [ 0, 1, 0 ], [ 2, 0, 1 ], [ 0, 0, 1 ], [ 3, 1, 0 ] ], -2 ],
    [ [ [ 0, 1, 0 ], [ 2, 0, 1 ], [ 1, 0, 0 ], [ 3, 1, 1 ] ], 1 ],
    [ [ [ 0, 1, 0 ], [ 2, 0, 2 ], [ 2, 0, 0 ], [ 1, 0, 0 ] ], -1 ],
    [ [ [ 0, 1, 0 ], [ 2, 0, 2 ], [ 3, 0, 0 ], [ 0, 1, 0 ] ], -1 ],
    [ [ [ 0, 1, 0 ], [ 2, 0, 2 ], [ 3, 1, 0 ], [ 1, 0, 0 ] ], -1 ],
    [ [ [ 0, 1, 0 ], [ 3, 0, 0 ], [ 1, 0, 0 ], [ 3, 1, 0 ] ], -2 ],
    [ [ [ 0, 1, 0 ], [ 3, 1, 1 ], [ 2, 0, 1 ], [ 1, 0, 0 ] ], 1 ],
    [ [ [ 0, 1, 1 ], [ 0, 1, 0 ], [ 1, 0, 0 ], [ 0, 0, 1 ] ], -1 ],
    [ [ [ 0, 1, 1 ], [ 0, 1, 0 ], [ 2, 0, 1 ], [ 1, 0, 0 ] ], -2 ],
    [ [ [ 0, 1, 1 ], [ 1, 0, 0 ], [ 0, 1, 0 ], [ 0, 0, 1 ] ], 1 ],
    [ [ [ 0, 1, 1 ], [ 1, 0, 0 ], [ 1, 0, 1 ], [ 0, 1, 0 ] ], -1 ],
    [ [ [ 0, 1, 1 ], [ 1, 0, 0 ], [ 1, 1, 1 ], [ 0, 0, 1 ] ], 1 ],
    [ [ [ 0, 1, 1 ], [ 1, 0, 0 ], [ 3, 1, 0 ], [ 2, 0, 1 ] ], 1 ],
    [ [ [ 0, 1, 1 ], [ 2, 0, 1 ], [ 0, 0, 1 ], [ 0, 1, 0 ] ], 1 ],
    [ [ [ 0, 1, 1 ], [ 2, 0, 1 ], [ 3, 1, 0 ], [ 0, 0, 1 ] ], -1 ],
    [ [ [ 0, 1, 1 ], [ 2, 1, 1 ], [ 1, 0, 0 ], [ 0, 0, 1 ] ], 4 ],
    [ [ [ 0, 1, 1 ], [ 3, 0, 1 ], [ 0, 1, 0 ], [ 1, 0, 0 ] ], -2 ],
    [ [ [ 0, 1, 1 ], [ 3, 0, 1 ], [ 1, 0, 0 ], [ 0, 0, 1 ] ], 2 ],
    [ [ [ 0, 1, 1 ], [ 3, 0, 1 ], [ 1, 0, 0 ], [ 3, 1, 0 ] ], -2 ],
    [ [ [ 0, 1, 1 ], [ 3, 0, 1 ], [ 1, 0, 1 ], [ 3, 1, 0 ] ], 2 ],
    [ [ [ 0, 1, 1 ], [ 3, 0, 1 ], [ 1, 1, 0 ], [ 0, 0, 1 ] ], -2 ],
    [ [ [ 1, 0, 0 ], [ 0, 0, 1 ], [ 0, 1, 0 ], [ 1, 1, 2 ] ], 1 ],
    [ [ [ 1, 0, 0 ], [ 0, 0, 1 ], [ 1, 2, 1 ], [ 0, 0, 1 ] ], -2 ],
    [ [ [ 1, 0, 0 ], [ 0, 0, 1 ], [ 1, 2, 2 ], [ 0, 0, 1 ] ], -1 ],
    [ [ [ 1, 0, 0 ], [ 0, 1, 0 ], [ 0, 0, 1 ], [ 0, 1, 0 ] ], -1 ],
    [ [ [ 1, 0, 0 ], [ 0, 1, 0 ], [ 0, 0, 2 ], [ 1, 0, 0 ] ], 2 ],
    [ [ [ 1, 0, 0 ], [ 0, 1, 0 ], [ 1, 0, 0 ], [ 0, 0, 2 ] ], -2 ],
    [ [ [ 1, 0, 0 ], [ 1, 0, 0 ], [ 0, 0, 2 ], [ 0, 1, 0 ] ], -2 ],
    [ [ [ 1, 0, 0 ], [ 1, 0, 0 ], [ 0, 1, 0 ], [ 0, 0, 2 ] ], 2 ],
    [ [ [ 1, 0, 0 ], [ 1, 0, 0 ], [ 3, 1, 0 ], [ 1, 0, 0 ] ], -1 ],
    [ [ [ 1, 0, 0 ], [ 1, 0, 1 ], [ 1, 0, 0 ], [ 1, 2, 1 ] ], 1 ],
    [ [ [ 1, 0, 0 ], [ 1, 0, 1 ], [ 1, 0, 0 ], [ 3, 1, 0 ] ], -1 ],
    [ [ [ 1, 0, 0 ], [ 1, 0, 1 ], [ 3, 1, 1 ], [ 0, 1, 0 ] ], -1 ],
    [ [ [ 1, 0, 0 ], [ 1, 0, 2 ], [ 3, 1, 0 ], [ 1, 0, 0 ] ], -2 ],
    [ [ [ 1, 0, 0 ], [ 3, 1, 0 ], [ 0, 0, 1 ], [ 1, 1, 1 ] ], -1 ],
    [ [ [ 1, 0, 0 ], [ 3, 1, 0 ], [ 0, 1, 0 ], [ 0, 0, 1 ] ], -1 ],
    [ [ [ 1, 0, 0 ], [ 3, 1, 0 ], [ 1, 0, 0 ], [ 1, 0, 0 ] ], 3 ],
    [ [ [ 1, 0, 0 ], [ 3, 1, 0 ], [ 1, 0, 0 ], [ 3, 1, 1 ] ], -1 ],
    [ [ [ 1, 0, 0 ], [ 3, 1, 0 ], [ 1, 0, 1 ], [ 0, 1, 0 ] ], 1 ],
    [ [ [ 1, 0, 0 ], [ 3, 1, 0 ], [ 2, 0, 0 ], [ 0, 0, 2 ] ], 3 ],
    [ [ [ 1, 0, 0 ], [ 3, 1, 0 ], [ 2, 0, 1 ], [ 0, 0, 1 ] ], -1 ],
    [ [ [ 1, 0, 0 ], [ 3, 1, 0 ], [ 2, 0, 2 ], [ 0, 1, 0 ] ], 1 ],
    [ [ [ 1, 0, 0 ], [ 3, 1, 1 ], [ 2, 0, 1 ], [ 0, 1, 0 ] ], -1 ],
    [ [ [ 1, 0, 1 ], [ 0, 1, 0 ], [ 2, 0, 2 ], [ 3, 1, 0 ] ], -2 ],
    [ [ [ 1, 0, 1 ], [ 1, 0, 0 ], [ 1, 1, 2 ], [ 0, 1, 0 ] ], 1 ],
    [ [ [ 1, 0, 1 ], [ 1, 0, 0 ], [ 1, 2, 1 ], [ 0, 0, 1 ] ], -1 ],
    [ [ [ 1, 0, 1 ], [ 1, 0, 0 ], [ 3, 1, 0 ], [ 2, 0, 2 ] ], 1 ],
    [ [ [ 1, 0, 1 ], [ 1, 1, 2 ], [ 0, 1, 0 ], [ 0, 0, 1 ] ], 1 ],
    [ [ [ 1, 0, 1 ], [ 1, 1, 2 ], [ 1, 0, 0 ], [ 3, 1, 0 ] ], 2 ],
    [ [ [ 1, 0, 1 ], [ 3, 1, 1 ], [ 0, 1, 0 ], [ 0, 0, 1 ] ], 1 ],
    [ [ [ 1, 0, 1 ], [ 3, 1, 1 ], [ 1, 0, 1 ], [ 0, 1, 0 ] ], -1 ],
    [ [ [ 1, 0, 1 ], [ 3, 1, 1 ], [ 2, 0, 1 ], [ 0, 1, 1 ] ], 1 ],
    [ [ [ 1, 0, 1 ], [ 3, 1, 1 ], [ 2, 0, 1 ], [ 1, 0, 0 ] ], -1 ],
    [ [ [ 1, 1, 0 ], [ 0, 0, 1 ], [ 0, 0, -1 ], [ 0, 0, 1 ] ], -3 ],
    [ [ [ 1, 1, 0 ], [ 0, 0, 2 ], [ 1, 0, 0 ], [ 3, 1, 0 ] ], 2 ],
    [ [ [ 1, 1, 1 ], [ 0, 0, 1 ], [ 0, 0, -1 ], [ 0, 0, 1 ] ], 4 ],
    [ [ [ 1, 1, 2 ], [ 0, 1, 0 ], [ 1, 0, 0 ], [ 0, 0, 1 ] ], -2 ],
    [ [ [ 1, 1, 2 ], [ 1, 0, 0 ], [ 0, 0, 1 ], [ 0, 1, 0 ] ], -1 ],
    [ [ [ 1, 1, 2 ], [ 1, 0, 0 ], [ 3, 1, 0 ], [ 0, 0, 1 ] ], 1 ],
    [ [ [ 1, 1, 2 ], [ 2, 0, 1 ], [ 0, 1, 0 ], [ 0, 0, 1 ] ], -2 ],
    [ [ [ 1, 1, 2 ], [ 2, 0, 1 ], [ 1, 0, 0 ], [ 3, 1, 0 ] ], -2 ],
    [ [ [ 1, 1, 2 ], [ 2, 0, 1 ], [ 1, 0, 1 ], [ 0, 1, 0 ] ], 1 ],
    [ [ [ 1, 1, 2 ], [ 3, 0, 1 ], [ 2, 0, 1 ], [ 0, 1, 0 ] ], -1 ],
    [ [ [ 2, 0, 0 ], [ 0, 0, 2 ], [ 0, 1, 0 ], [ 2, 0, 0 ] ], 1 ],
    [ [ [ 2, 0, 0 ], [ 0, 0, 2 ], [ 2, 0, 0 ], [ 0, 1, 0 ] ], -1 ],
    [ [ [ 2, 0, 0 ], [ 0, 1, 0 ], [ 0, 0, 2 ], [ 2, 0, 0 ] ], -1 ],
    [ [ [ 2, 0, 0 ], [ 0, 1, 0 ], [ 1, 0, 0 ], [ 1, 0, 0 ] ], 1 ],
    [ [ [ 2, 0, 0 ], [ 0, 1, 0 ], [ 2, 0, 0 ], [ 2, 0, 2 ] ], 1 ],
    [ [ [ 2, 0, 0 ], [ 1, 0, 0 ], [ 0, 1, 0 ], [ 1, 0, 0 ] ], -1 ],
    [ [ [ 2, 0, 0 ], [ 1, 0, 0 ], [ 1, 0, 0 ], [ 0, 1, 0 ] ], 1 ],
    [ [ [ 2, 0, 0 ], [ 2, 0, 0 ], [ 0, 1, 0 ], [ 2, 0, 2 ] ], -1 ],
    [ [ [ 2, 0, 0 ], [ 2, 0, 0 ], [ 2, 0, 2 ], [ 0, 1, 0 ] ], 1 ],
    [ [ [ 2, 0, 1 ], [ 0, 0, 1 ], [ 0, 1, 0 ], [ 3, 0, 0 ] ], -1 ],
    [ [ [ 2, 0, 1 ], [ 0, 0, 1 ], [ 0, 1, 0 ], [ 3, 1, 0 ] ], -2 ],
    [ [ [ 2, 0, 1 ], [ 0, 0, 1 ], [ 3, 1, 0 ], [ 0, 1, 0 ] ], 1 ],
    [ [ [ 2, 0, 1 ], [ 0, 1, 1 ], [ 0, 1, 0 ], [ 1, 0, 1 ] ], -2 ],
    [ [ [ 2, 0, 1 ], [ 0, 1, 1 ], [ 1, 1, 1 ], [ 0, 0, 1 ] ], -2 ],
    [ [ [ 2, 0, 1 ], [ 1, 0, 0 ], [ 1, 2, 1 ], [ 0, 0, 1 ] ], 2 ],
    [ [ [ 2, 0, 1 ], [ 1, 0, 0 ], [ 1, 2, 2 ], [ 0, 0, 1 ] ], 2 ],
    [ [ [ 2, 0, 1 ], [ 1, 0, 0 ], [ 3, 1, 0 ], [ 0, 1, 0 ] ], 1 ],
    [ [ [ 2, 0, 1 ], [ 1, 0, 0 ], [ 3, 1, 1 ], [ 1, 1, 2 ] ], -1 ],
    [ [ [ 2, 0, 1 ], [ 1, 1, 2 ], [ 2, 0, 1 ], [ 0, 1, 0 ] ], 1 ],
    [ [ [ 2, 0, 1 ], [ 3, 1, 0 ], [ 0, 0, 1 ], [ 3, 0, 0 ] ], 1 ],
    [ [ [ 2, 0, 1 ], [ 3, 1, 0 ], [ 0, 0, 1 ], [ 3, 1, 0 ] ], 1 ],
    [ [ [ 2, 0, 1 ], [ 3, 1, 1 ], [ 3, 1, 0 ], [ 1, 0, 0 ] ], 1 ],
    [ [ [ 2, 0, 1 ], [ 3, 2, 1 ], [ 1, 0, 0 ], [ 0, 0, 1 ] ], -2 ],
    [ [ [ 2, 0, 2 ], [ 0, 1, 0 ], [ 2, 0, 0 ], [ 1, 0, 0 ] ], 1 ],
    [ [ [ 2, 0, 2 ], [ 2, 0, 0 ], [ 0, 1, 0 ], [ 1, 0, 0 ] ], -1 ],
    [ [ [ 2, 0, 2 ], [ 2, 0, 0 ], [ 1, 0, 0 ], [ 0, 1, 0 ] ], 1 ],
    [ [ [ 2, 0, 2 ], [ 3, 1, 0 ], [ 1, 0, 0 ], [ 3, 1, 0 ] ], -2 ],
    [ [ [ 2, 1, 1 ], [ 0, 0, 1 ], [ 3, 0, 0 ], [ 0, 1, 0 ] ], 1 ],
    [ [ [ 3, 0, 0 ], [ 1, 0, 0 ], [ 3, 1, 0 ], [ 1, 0, 0 ] ], 1 ],
    [ [ [ 3, 0, 1 ], [ 0, 1, 0 ], [ 1, 0, 0 ], [ 0, 0, 1 ] ], 4 ],
    [ [ [ 3, 0, 1 ], [ 0, 1, 0 ], [ 1, 0, 0 ], [ 3, 1, 1 ] ], -1 ],
    [ [ [ 3, 0, 1 ], [ 0, 1, 0 ], [ 1, 0, 1 ], [ 0, 1, 0 ] ], 1 ],
    [ [ [ 3, 0, 1 ], [ 1, 0, 0 ], [ 0, 0, 1 ], [ 3, 1, 0 ] ], 1 ],
    [ [ [ 3, 0, 1 ], [ 1, 0, 0 ], [ 0, 0, 1 ], [ 3, 2, 0 ] ], 1 ],
    [ [ [ 3, 0, 1 ], [ 1, 0, 0 ], [ 3, 1, 0 ], [ 0, 1, 1 ] ], -1 ],
    [ [ [ 3, 0, 1 ], [ 1, 0, 0 ], [ 3, 1, 0 ], [ 1, 0, 1 ] ], -1 ],
    [ [ [ 3, 0, 1 ], [ 1, 0, 1 ], [ 3, 1, 0 ], [ 0, 1, 0 ] ], 1 ],
    [ [ [ 3, 0, 1 ], [ 1, 1, 0 ], [ 0, 0, 1 ], [ 0, 0, -1 ] ], 3 ],
    [ [ [ 3, 0, 1 ], [ 1, 1, 0 ], [ 0, 0, 1 ], [ 0, 1, 0 ] ], -1 ],
    [ [ [ 3, 0, 1 ], [ 1, 1, 1 ], [ 0, 0, -1 ], [ 0, 0, 1 ] ], -4 ],
    [ [ [ 3, 1, 0 ], [ 0, 0, 1 ], [ 3, 1, 0 ], [ 2, 0, 1 ] ], -1 ],
    [ [ [ 3, 1, 0 ], [ 1, 0, 0 ], [ 1, 0, 0 ], [ 0, 1, 0 ] ], -3 ],
    [ [ [ 3, 1, 0 ], [ 1, 0, 0 ], [ 1, 0, 1 ], [ 3, 1, 1 ] ], 1 ],
    [ [ [ 3, 1, 0 ], [ 1, 0, 0 ], [ 1, 1, 0 ], [ 0, 0, 2 ] ], -3 ],
    [ [ [ 3, 1, 0 ], [ 1, 0, 0 ], [ 3, 1, 1 ], [ 2, 0, 1 ] ], 1 ],
    [ [ [ 3, 1, 0 ], [ 2, 0, 0 ], [ 0, 0, 2 ], [ 0, 1, 0 ] ], -3 ],
    [ [ [ 3, 1, 0 ], [ 2, 0, 0 ], [ 0, 1, 0 ], [ 0, 0, 2 ] ], 3 ],
    [ [ [ 3, 1, 0 ], [ 2, 0, 1 ], [ 0, 0, 1 ], [ 0, 1, 0 ] ], 1 ],
    [ [ [ 3, 1, 1 ], [ 2, 0, 1 ], [ 0, 1, 0 ], [ 0, 0, 1 ] ], 1 ],
    [ [ [ 3, 1, 1 ], [ 3, 0, 0 ], [ 0, 1, 0 ], [ 1, 0, 0 ] ], -1 ],
    [ [ [ 3, 1, 1 ], [ 3, 0, 0 ], [ 1, 0, 0 ], [ 0, 1, 0 ] ], 1 ],
    [ [ [ 3, 1, 1 ], [ 3, 0, 0 ], [ 1, 0, 0 ], [ 3, 1, 0 ] ], 1 ],
];

if Length(Target) <> 76 then
    Error("FAIL: target must have 76 Laurent atoms");
fi;
if L1Sparse(Target) <> 76 then
    Error("FAIL: target must have l1 norm 76");
fi;
if Length(Certificate) <> 146 then
    Error("FAIL: certificate must contain 146 cocycle defects");
fi;
if Sum(Certificate, entry -> AbsInt(entry[2])) <> 245 then
    Error("FAIL: certificate coefficient l1 norm must be 245");
fi;

CertificateVector := [];
for entry in Certificate do
    CertificateVector := AddSparse(
        CertificateVector,
        ScaleSparse(CocycleDefect(
            entry[1][1], entry[1][2], entry[1][3], entry[1][4]
        ), entry[2])
    );
od;

Residual := AddSparse(Target, ScaleSparse(CertificateVector, -1));
if Length(Residual) <> 0 then
    Error("FAIL: nonzero Laurent residual; atoms = ", Length(Residual));
fi;

Print("Gamma4 terminal-root self-character certificate\n");
Print("Arithmetic: exact integers in a formal sparse Laurent lattice\n");
Print("External GAP packages loaded: none\n");
Print("Target Laurent atoms: ", Length(Target), "\n");
Print("Target l1 norm: ", L1Sparse(Target), "\n");
Print("Certificate cocycle defects: ", Length(Certificate), "\n");
Print("Certificate coefficient l1 norm: ",
      Sum(Certificate, entry -> AbsInt(entry[2])), "\n");
Print("Residual Laurent atoms: ", Length(Residual), "\n");
Print("PASS: Q is the stated product of normalized 3-cocycle defects, hence Q=1\n");

QUIT_GAP(0);
