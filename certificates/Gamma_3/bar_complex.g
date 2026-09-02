# Core-only exact normalized bar-complex utilities.
# A chain is a list [ [bar, integerCoefficient], ... ].

if not IsBound(G3_Fail) then
    G3_Fail := function(arg)
        CallFuncList(PrintTo,
            Concatenation(["*errout*"], arg, ["\n"]));
        ErrorNoReturn("Certificate verification failed.");
    end;
fi;

BC_Slice := function(list, first, last)
    if first > last then
        return [];
    fi;
    return list{[first..last]};
end;
BC_Normalize := function(chain, one)
    local out, term, pos, entry;
    out := [];
    for term in chain do
        if term[2] <> 0 and ForAll(term[1], x -> x <> one) then
            pos := PositionProperty(out, y -> y[1] = term[1]);
            if pos = fail then
                Add(out, [ShallowCopy(term[1]), term[2]]);
            else
                out[pos][2] := out[pos][2] + term[2];
            fi;
        fi;
    od;
    out := Filtered(out, x -> x[2] <> 0);
    Sort(out, function(x, y) return x[1] < y[1]; end);
    return out;
end;

BC_Bar := function(one, entries)
    if ForAny(entries, x -> x = one) then
        return [];
    fi;
    return [[ShallowCopy(entries), 1]];
end;

BC_Add := function(left, right, one)
    return BC_Normalize(Concatenation(left, right), one);
end;

BC_Scale := function(scalar, chain, one)
    return BC_Normalize(List(chain, x -> [x[1], scalar * x[2]]), one);
end;

BC_Sub := function(left, right, one)
    return BC_Add(left, BC_Scale(-1, right, one), one);
end;

BC_Sum := function(chains, one)
    local answer, chain;
    answer := [];
    for chain in chains do
        answer := BC_Add(answer, chain, one);
    od;
    return answer;
end;

BC_Boundary := function(chain, mul, one)
    local answer, term, bar, coeff, n, i, newbar, sign;
    answer := [];
    for term in chain do
        bar := term[1];
        coeff := term[2];
        n := Length(bar);

        newbar := BC_Slice(bar, 2, n);
        answer := BC_Add(answer, BC_Scale(coeff, BC_Bar(one, newbar), one), one);

        for i in [1..n-1] do
            newbar := Concatenation(
                BC_Slice(bar, 1, i-1),
                [mul(bar[i], bar[i+1])],
                BC_Slice(bar, i+2, n)
            );
            sign := (-1)^i;
            answer := BC_Add(answer,
                BC_Scale(coeff * sign, BC_Bar(one, newbar), one), one);
        od;

        newbar := BC_Slice(bar, 1, n-1);
        sign := (-1)^n;
        answer := BC_Add(answer,
            BC_Scale(coeff * sign, BC_Bar(one, newbar), one), one);
    od;
    return BC_Normalize(answer, one);
end;

BC_AssertZero := function(label, chain, one)
    local normalized;
    normalized := BC_Normalize(chain, one);
    if Length(normalized) <> 0 then
        G3_Fail(label, ": nonzero residual: ", normalized);
    fi;
end;

BC_AssertEqual := function(label, left, right, one)
    local residual;
    residual := BC_Sub(left, right, one);
    BC_AssertZero(label, residual, one);
end;

BC_Project := function(chain, project, oneSource, oneTarget)
    local answer, term, imagebar;
    answer := [];
    for term in chain do
        imagebar := List(term[1], project);
        answer := BC_Add(answer,
            BC_Scale(term[2], BC_Bar(oneTarget, imagebar), oneTarget),
            oneTarget);
    od;
    return answer;
end;
