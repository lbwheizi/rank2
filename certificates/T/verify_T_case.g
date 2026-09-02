#############################################################################
#
# verify_T_case.g
#
# Exact GAP-core verifier for the finite coefficient calculations in
# Section 6 and Appendix D of the accompanying manuscript.
#
# No external GAP package and no floating-point arithmetic are used.
# A polynomial a_0 + a_1*p + ... + a_d*p^d is stored as
# [ a_0, a_1, ..., a_d ].
#
#############################################################################

Print("T-case exact core verifier (Appendix D)\n");
Print("Arithmetic: exact coefficient lists in Z[p].\n");
Print("Dependencies: GAP core only; no packages or floating point.\n");

T_CHECKS := 0;;

T_Fail := function(arg)
    PrintTo("*errout*", "FAIL: ", Concatenation(arg), "\n");
    ErrorNoReturn("Certificate verification failed.");
end;;

T_Check := function(condition, message)
    if not condition then
        T_Fail(message);
    fi;
    T_CHECKS := T_CHECKS + 1;
    Print("[OK] ", message, "\n");
end;;

#############################################################################
# Exact polynomial arithmetic in Z[p].
#############################################################################

T_PolyNormalize := function(a)
    local b;
    b := ShallowCopy(a);
    while Length(b) > 1 and b[Length(b)] = 0 do
        Remove(b);
    od;
    if Length(b) = 0 then
        return [0];
    fi;
    return b;
end;;

T_PolyIsZero := function(a)
    local b;
    b := T_PolyNormalize(a);
    return Length(b) = 1 and b[1] = 0;
end;;

T_PolyAdd := function(a, b)
    local c, i, n;
    n := Maximum(Length(a), Length(b));
    c := List([1 .. n], i -> 0);
    for i in [1 .. Length(a)] do
        c[i] := c[i] + a[i];
    od;
    for i in [1 .. Length(b)] do
        c[i] := c[i] + b[i];
    od;
    return T_PolyNormalize(c);
end;;

T_PolyNeg := function(a)
    return List(a, x -> -x);
end;;

T_PolySub := function(a, b)
    return T_PolyAdd(a, T_PolyNeg(b));
end;;

T_PolyMul := function(a, b)
    local c, i, j;
    c := List([1 .. Length(a) + Length(b) - 1], i -> 0);
    for i in [1 .. Length(a)] do
        for j in [1 .. Length(b)] do
            c[i + j - 1] := c[i + j - 1] + a[i] * b[j];
        od;
    od;
    return T_PolyNormalize(c);
end;;

# Exact remainder modulo f(p)=p^2-p+1.  The replacement
# p^d = p^(d-1)-p^(d-2) is applied from the largest degree downwards.
T_PolyReduceF := function(a)
    local b, c, n;
    b := T_PolyNormalize(a);
    while Length(b) > 2 do
        c := b[Length(b)];
        Remove(b);
        n := Length(b);
        b[n] := b[n] + c;
        b[n - 1] := b[n - 1] - c;
    od;
    return T_PolyNormalize(b);
end;;

T_ZERO := [0];;
T_ONE := [1];;
T_MINUS_ONE := [-1];;
T_P := [0, 1];;
T_MINUS_P := [0, -1];;
T_P_MINUS_ONE := [-1, 1];;
T_ONE_MINUS_P := [1, -1];;
T_F := [1, -1, 1];;

#############################################################################
# Dense vectors in the ordered word basis.
#
# A word (i_1,...,i_m), 1 <= i_k <= 4, is encoded lexicographically in
# base four.  At m=4 the ambient vector has only 256 coordinates.
#############################################################################

T_WordIndex := function(word)
    local idx, x;
    idx := 1;
    for x in word do
        idx := 4 * (idx - 1) + x;
    od;
    return idx;
end;;

T_IndexWord := function(m, idx)
    local k, n, word;
    n := idx - 1;
    word := List([1 .. m], k -> 1);
    for k in [m, m - 1 .. 1] do
        word[k] := (n mod 4) + 1;
        n := QuoInt(n, 4);
    od;
    return word;
end;;

T_VecZero := function(m)
    return rec(
        m := m,
        coeffs := List([1 .. 4 ^ m], i -> [0])
    );
end;;

T_VecFromTerms := function(m, terms)
    local idx, term, v;
    v := T_VecZero(m);
    for term in terms do
        if Length(term[1]) <> m then
            T_Fail("A displayed word has the wrong length");
        fi;
        idx := T_WordIndex(term[1]);
        v.coeffs[idx] := T_PolyAdd(v.coeffs[idx], term[2]);
    od;
    return v;
end;;

T_VecAdd := function(a, b)
    local i, out;
    if a.m <> b.m then
        T_Fail("Vector-length mismatch in T_VecAdd");
    fi;
    out := T_VecZero(a.m);
    for i in [1 .. Length(a.coeffs)] do
        out.coeffs[i] := T_PolyAdd(a.coeffs[i], b.coeffs[i]);
    od;
    return out;
end;;

T_VecNeg := function(a)
    local i, out;
    out := T_VecZero(a.m);
    for i in [1 .. Length(a.coeffs)] do
        out.coeffs[i] := T_PolyNeg(a.coeffs[i]);
    od;
    return out;
end;;

T_VecScale := function(a, scalar)
    local i, out;
    out := T_VecZero(a.m);
    for i in [1 .. Length(a.coeffs)] do
        out.coeffs[i] := T_PolyMul(a.coeffs[i], scalar);
    od;
    return out;
end;;

T_VecReduceF := function(a)
    local i, out;
    out := T_VecZero(a.m);
    for i in [1 .. Length(a.coeffs)] do
        out.coeffs[i] := T_PolyReduceF(a.coeffs[i]);
    od;
    return out;
end;;

T_VecEqual := function(a, b)
    if a.m <> b.m then
        return false;
    fi;
    return a.coeffs = b.coeffs;
end;;

T_VecIsZero := function(a)
    return ForAll(a.coeffs, T_PolyIsZero);
end;;

T_VecCoefficient := function(a, word)
    return a.coeffs[T_WordIndex(word)];
end;;

T_VecSupport := function(a)
    local i, support;
    support := [];
    for i in [1 .. Length(a.coeffs)] do
        if not T_PolyIsZero(a.coeffs[i]) then
            Add(support, T_IndexWord(a.m, i));
        fi;
    od;
    return support;
end;;

T_VecPrefix := function(label, a)
    local i, out, word;
    out := T_VecZero(a.m + 1);
    for i in [1 .. Length(a.coeffs)] do
        if not T_PolyIsZero(a.coeffs[i]) then
            word := Concatenation([label], T_IndexWord(a.m, i));
            out.coeffs[T_WordIndex(word)] := a.coeffs[i];
        fi;
    od;
    return out;
end;;

#############################################################################
# The tetrahedral rack and the constant (-1) coefficient table in the
# compatible bases fixed in the manuscript.
#
# The rows are reconstructed from
# Ad(x_1)=(2,4,3), Ad(x_2)=(1,3,4),
# Ad(x_3)=(1,4,2), Ad(x_4)=(1,2,3).
#############################################################################

T_AD := [(2,4,3), (1,3,4), (1,4,2), (1,2,3)];;
T_RACK := List([1 .. 4], function(i)
    return List([1 .. 4], j -> j ^ T_AD[i]);
end);;

T_EXPECTED_RACK := [
    [1, 4, 2, 3],
    [3, 2, 4, 1],
    [4, 1, 3, 2],
    [2, 3, 1, 4]
];;

T_RackOp := function(i, j)
    return T_RACK[i][j];
end;;

T_BRAID_COEFFS := List([1 .. 4],
    i -> List([1 .. 4], j -> -1));;
T_ok := T_RACK = T_EXPECTED_RACK;;
for T_i in [1 .. 4] do
    T_ok := T_ok and T_RackOp(T_i, T_i) = T_i;
    T_ok := T_ok and Set(T_RACK[T_i]) = [1 .. 4];
    for T_j in [1 .. 4] do
        for T_k in [1 .. 4] do
            T_ok := T_ok and
                T_RackOp(T_i, T_RackOp(T_j, T_k)) =
                T_RackOp(T_RackOp(T_i, T_j), T_RackOp(T_i, T_k));
            T_ok := T_ok and
                T_BRAID_COEFFS[T_i][T_k] *
                T_BRAID_COEFFS[T_RackOp(T_i,T_j)][T_RackOp(T_i,T_k)] =
                T_BRAID_COEFFS[T_j][T_k] *
                T_BRAID_COEFFS[T_i][T_RackOp(T_j,T_k)];
        od;
    od;
od;
# For these compatible bases and the assumed value epsilon=1, one has
# s=-1 and epsilon=s*b_23*b_42*b_34.
T_ok := T_ok and
    (-1) * T_BRAID_COEFFS[2][3] * T_BRAID_COEFFS[4][2] *
    T_BRAID_COEFFS[3][4] = 1;;
T_Check(T_ok,
    "lem:T-tetrahedral-braiding: tetrahedral rack table and all 64 constant-(-1) braid instances");;

#############################################################################
# Adjacent label braids and the h_m braid loop.
# The marker 0 denotes the central V-strand.
#############################################################################

T_BraidLabels := function(state, position)
    local a, b, out;
    out := ShallowCopy(state);
    a := out[position];
    b := out[position + 1];
    if a = 0 or b = 0 then
        out[position] := b;
        out[position + 1] := a;
    else
        out[position] := T_RackOp(a, b);
        out[position + 1] := a;
    fi;
    return out;
end;;

T_HWord := function(word)
    local k, m, state;
    m := Length(word);
    state := Concatenation(ShallowCopy(word), [0]);
    for k in [1 .. m] do
        state := T_BraidLabels(state, k);
    od;
    for k in [m, m - 1 .. 1] do
        state := T_BraidLabels(state, k);
    od;
    if state[m + 1] <> 0 then
        T_Fail("The h_m braid loop did not return V to the last position");
    fi;
    return state{[1 .. m]};
end;;

T_ok := true;;
for T_i in [1 .. 4] do
    for T_j in [1 .. 4] do
        T_ok := T_ok and T_HWord([T_i, T_j]) =
            [T_RackOp(T_RackOp(T_i, T_j), T_i), T_RackOp(T_i, T_j)];
    od;
od;
T_Check(T_ok,
    "eq:T-adjoint-recursion: displayed h_m braid-loop direction (all h_2 words)");;

#############################################################################
# Compatible-basis coefficient recursion eq:T-adjoint-recursion (6.17).
#############################################################################

T_PhiWord := fail;;
T_PhiWord := function(word)
    local hword, inner, lower, m, out, prefix;
    m := Length(word);
    if m = 1 then
        return T_VecFromTerms(1, [[word, T_ONE_MINUS_P]]);
    fi;

    out := T_VecFromTerms(m, [
        [word, T_ONE],
        [T_HWord(word), T_MINUS_P]
    ]);

    inner := [word[1]];
    if m >= 3 then
        Append(inner, word{[3 .. m]});
    fi;
    lower := T_PhiWord(inner);
    prefix := T_RackOp(word[1], word[2]);
    out := T_VecAdd(out, T_VecNeg(T_VecPrefix(prefix, lower)));
    return out;
end;;

T_PhiOnPrefixedVector := function(prefix, a)
    local i, image, out, word;
    out := T_VecZero(a.m + 1);
    for i in [1 .. Length(a.coeffs)] do
        if not T_PolyIsZero(a.coeffs[i]) then
            word := Concatenation([prefix], T_IndexWord(a.m, i));
            image := T_PhiWord(word);
            out := T_VecAdd(out, T_VecScale(image, a.coeffs[i]));
        fi;
    od;
    return out;
end;;

#############################################################################
# eq:T-Y2-basis and eq:T-phi2-table (D.2--D.3):
# u_i and the complete phi_2 table, modulo f(p).
#############################################################################

T_u := [];;
T_u[1] := T_VecFromTerms(2, [
    [[2,3], T_ONE], [[3,4], T_MINUS_P], [[4,2], T_P_MINUS_ONE]
]);;
T_u[2] := T_VecFromTerms(2, [
    [[1,4], T_ONE], [[4,3], T_MINUS_P], [[3,1], T_P_MINUS_ONE]
]);;
T_u[3] := T_VecFromTerms(2, [
    [[1,2], T_ONE], [[2,4], T_MINUS_P], [[4,1], T_P_MINUS_ONE]
]);;
T_u[4] := T_VecFromTerms(2, [
    [[1,3], T_ONE], [[3,2], T_MINUS_P], [[2,1], T_P_MINUS_ONE]
]);;

# A non-fail cell [j,a] means a*u_j.
T_PHI2_SPEC := [
    [fail,             [3,T_ONE],         [4,T_ONE],         [2,T_ONE]],
    [[4,T_MINUS_P],    fail,              [1,T_ONE],         [3,T_P_MINUS_ONE]],
    [[2,T_MINUS_P],    [4,T_P_MINUS_ONE], fail,              [1,T_P_MINUS_ONE]],
    [[3,T_MINUS_P],    [1,T_MINUS_P],     [2,T_P_MINUS_ONE], fail]
];;

T_ok := true;;
for T_i in [1 .. 4] do
    for T_j in [1 .. 4] do
        T_actual := T_VecReduceF(T_PhiWord([T_i, T_j]));;
        T_spec := T_PHI2_SPEC[T_i][T_j];;
        if T_spec = fail then
            T_expected := T_VecZero(2);;
        else
            T_expected := T_VecReduceF(
                T_VecScale(T_u[T_spec[1]], T_spec[2]));;
        fi;
        T_ok := T_ok and T_VecEqual(T_actual, T_expected);
    od;
od;
T_Check(T_ok,
    "eq:T-phi2-table: complete phi_2 table modulo p^2-p+1");;

#############################################################################
# eq:T-y3-generator (6.18) and eq:T-phi3-table (D.4):
# y and the complete phi_3 table, modulo f(p).
#############################################################################

T_y := T_VecFromTerms(3, [
    [[1,2,3], T_ONE],
    [[1,4,2], T_P_MINUS_ONE],
    [[4,1,3], T_P_MINUS_ONE],
    [[4,3,2], T_ONE],
    [[4,2,1], T_MINUS_P],
    [[1,3,4], T_MINUS_P],
    [[2,1,4], T_ONE],
    [[2,4,3], T_MINUS_P],
    [[2,3,1], T_P_MINUS_ONE],
    [[3,1,2], T_MINUS_P],
    [[3,2,4], T_P_MINUS_ONE],
    [[3,4,1], T_ONE]
]);;

T_PHI3_DIAGONAL := [T_P, T_P, T_MINUS_ONE, T_ONE_MINUS_P];;
T_ok := true;;
for T_i in [1 .. 4] do
    for T_j in [1 .. 4] do
        T_actual := T_VecReduceF(T_PhiOnPrefixedVector(T_i, T_u[T_j]));;
        if T_i = T_j then
            T_expected := T_VecReduceF(
                T_VecScale(T_y, T_PHI3_DIAGONAL[T_i]));;
        else
            T_expected := T_VecZero(3);;
        fi;
        T_ok := T_ok and T_VecEqual(T_actual, T_expected);
    od;
od;
T_Check(T_ok,
    "eq:T-phi3-table: complete phi_3 table modulo p^2-p+1");;

#############################################################################
# eq:T-Y4-homogeneous-components and eq:T-Y4-Z1 (D.5--D.6):
# exact, unreduced identity phi_4(w_1 tensor y)=f(p) Z_1.
#############################################################################

T_Z1 := T_VecFromTerms(4, [
    [[1,4,3,2], T_ONE],
    [[1,4,2,1], T_MINUS_ONE],
    [[1,3,1,2], T_MINUS_ONE],
    [[1,3,2,4], T_MINUS_ONE],
    [[1,3,4,1], [2]],
    [[3,2,1,2], T_ONE],
    [[3,1,3,2], T_MINUS_ONE],
    [[3,1,2,1], T_ONE],
    [[3,2,4,1], T_MINUS_ONE],
    [[1,2,1,4], T_ONE],
    [[1,2,3,1], T_MINUS_ONE],
    [[2,4,1,4], T_MINUS_ONE],
    [[2,1,2,4], T_ONE],
    [[2,1,4,1], T_MINUS_ONE],
    [[2,4,3,1], T_ONE]
]);;

T_actual4 := T_PhiOnPrefixedVector(1, T_y);;
T_expected4 := T_VecScale(T_Z1, T_F);;
T_Check(Length(T_VecSupport(T_Z1)) = 15,
    "eq:T-Y4-Z1: fifteen pairwise distinct ordered words");;
T_Check(T_VecEqual(T_actual4, T_expected4),
    "eq:T-Y4-homogeneous-components and eq:T-Y4-Z1: exact factorization in unreduced Z[p]");;
T_Check(T_VecIsZero(T_VecReduceF(T_actual4)),
    "eq:T-Y4-homogeneous-components: vanishing after exact factorization is reduced modulo f");;

#############################################################################
# eq:T-R2-monodromy-cycles and eq:T-R2-mixed-monodromy
# (6.22--6.23): derive the four monodromy cycles from adjacent braids.
#############################################################################

# signature = [sign, exponent of r, exponent of ell].
T_BraidWithSignature := function(state, position, signature)
    local a, b, out, sig;
    out := ShallowCopy(state);
    sig := ShallowCopy(signature);
    a := out[position];
    b := out[position + 1];
    if a = 0 then
        sig[3] := sig[3] + 1;
        out[position] := b;
        out[position + 1] := 0;
    elif b = 0 then
        sig[2] := sig[2] + 1;
        out[position] := 0;
        out[position + 1] := a;
    else
        sig[1] := -sig[1];
        out[position] := T_RackOp(a, b);
        out[position + 1] := a;
    fi;
    return rec(state := out, signature := sig);
end;;

T_MonodromyOuterOne := function(triple)
    local k, result, signature, state;
    state := Concatenation([1], ShallowCopy(triple), [0]);
    signature := [1, 0, 0];
    for k in [1,2,3,4,4,3,2,1] do
        result := T_BraidWithSignature(state, k, signature);
        state := result.state;
        signature := result.signature;
    od;
    if state[1] <> 1 or state[5] <> 0 then
        T_Fail("The R2 monodromy braid did not close");
    fi;
    return rec(triple := state{[2 .. 4]}, signature := signature);
end;;

T_CYCLE_STARTS := [[1,2,3], [4,1,3], [4,3,2], [4,2,1]];;
T_EXPECTED_CYCLES := [
    [[1,2,3], [1,4,2], [1,3,4]],
    [[4,1,3], [3,1,2], [2,1,4]],
    [[4,3,2], [3,2,4], [2,4,3]],
    [[4,2,1], [3,4,1], [2,3,1]]
];;

T_cycles := [];;
T_ok := true;;
for T_start in T_CYCLE_STARTS do
    T_current := ShallowCopy(T_start);;
    T_cycle := [];;
    for T_n in [1 .. 3] do
        Add(T_cycle, ShallowCopy(T_current));
        T_result := T_MonodromyOuterOne(T_current);;
        # Six W-W crossings give sign +1, and the two mixed crossings
        # contribute exactly r*ell=p.
        T_ok := T_ok and T_result.signature = [1,1,1];
        T_current := T_result.triple;
    od;
    T_ok := T_ok and T_current = T_start;
    Add(T_cycles, T_cycle);
od;
T_ok := T_ok and T_cycles = T_EXPECTED_CYCLES;
T_ok := T_ok and Set(Concatenation(T_cycles)) = Set(T_VecSupport(T_y));
T_Check(T_ok,
    "eq:T-R2-monodromy-cycles: four cycles derived from the eight adjacent braids");;

T_EXPECTED_TRIPLES := [
    [T_ONE, T_P_MINUS_ONE, T_MINUS_P],
    [T_P_MINUS_ONE, T_MINUS_P, T_ONE],
    [T_ONE, T_P_MINUS_ONE, T_MINUS_P],
    [T_MINUS_P, T_ONE, T_P_MINUS_ONE]
];;

T_coeff_triples := List(T_cycles, function(cycle)
    return List(cycle, word -> T_VecCoefficient(T_y, word));
end);;
T_Check(T_coeff_triples = T_EXPECTED_TRIPLES,
    "eq:T-y3-generator and eq:T-R2-monodromy-cycles: four displayed coefficient triples of y");;

T_ok := true;;
for T_coeffs in T_coeff_triples do
    # If e_1 -> p e_2 -> p e_3 -> p e_1, the coefficient vector is
    # sent to p*(c_3,c_1,c_2).
    T_image_coeffs := [
        T_PolyMul(T_P, T_coeffs[3]),
        T_PolyMul(T_P, T_coeffs[1]),
        T_PolyMul(T_P, T_coeffs[2])
    ];;
    T_target_coeffs := List(T_coeffs,
        c -> T_PolyMul(T_ONE_MINUS_P, c));;
    T_image_coeffs := List(T_image_coeffs, T_PolyReduceF);;
    T_target_coeffs := List(T_target_coeffs, T_PolyReduceF);;
    T_ok := T_ok and T_image_coeffs = T_target_coeffs;
od;
T_Check(T_ok,
    "eq:T-R2-mixed-monodromy: eigenvalue (1-p) on all four cycles");;

Print("PASS: ", T_CHECKS, " exact checks completed.\n");
