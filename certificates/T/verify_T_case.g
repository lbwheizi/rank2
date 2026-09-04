#############################################################################
#
# verify_T_case.g
#
# Exact GAP-core verifier for the finite coefficient calculations in
# Section 6 and Appendix D of the accompanying manuscript.
#
# No external GAP package and no floating-point arithmetic are used.
# Laurent monomials in named local-cocycle and projective-representation
# values are stored by exact exponent lists.  A polynomial
# a_0 + a_1*p + ... + a_d*p^d is stored as [ a_0, a_1, ..., a_d ].
#
#############################################################################

Print("T-case exact core verifier (Appendix D)\n");
Print("Arithmetic: exact Laurent monomials and coefficient lists in Z[p].\n");
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
# Exact Laurent monomials.
#
# A monomial is [coefficient, factors], where factors is a list of pairs
# [name, exponent].  Normalization collects equal names, deletes exponent
# zero, and sorts the remaining pairs.  All calculations below use only
# coefficients +1 and -1, but inversion is kept exact over the rationals.
#############################################################################

T_LaurentNormalize := function(m)
    local factor, factors, out, pos;
    factors := [];
    for factor in m[2] do
        pos := PositionProperty(factors, x -> x[1] = factor[1]);
        if pos = fail then
            Add(factors, [factor[1], factor[2]]);
        else
            factors[pos][2] := factors[pos][2] + factor[2];
        fi;
    od;
    out := Filtered(factors, x -> x[2] <> 0);
    Sort(out, function(a, b) return a[1] < b[1]; end);
    return [m[1], out];
end;;

T_LaurentOne := [1, []];;

T_LaurentConstant := function(c)
    return [c, []];
end;;

T_LaurentGenerator := function(name)
    return [1, [[name, 1]]];
end;;

T_LaurentMul := function(a, b)
    return T_LaurentNormalize([
        a[1] * b[1],
        Concatenation(a[2], b[2])
    ]);
end;;

T_LaurentInv := function(a)
    return T_LaurentNormalize([
        1 / a[1],
        List(a[2], x -> [x[1], -x[2]])
    ]);
end;;

T_LaurentEqual := function(a, b)
    return T_LaurentNormalize(a) = T_LaurentNormalize(b);
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

T_ok := T_RACK = T_EXPECTED_RACK;;
for T_i in [1 .. 4] do
    T_ok := T_ok and T_RackOp(T_i, T_i) = T_i;
    T_ok := T_ok and Set(T_RACK[T_i]) = [1 .. 4];
    for T_j in [1 .. 4] do
        for T_k in [1 .. 4] do
            T_ok := T_ok and
                T_RackOp(T_i, T_RackOp(T_j, T_k)) =
                T_RackOp(T_RackOp(T_i, T_j), T_RackOp(T_i, T_k));
        od;
    od;
od;
T_Check(T_ok,
    "lem:T-group-data: tetrahedral rack reconstructed from the four conjugation permutations");;

#############################################################################
# Part I: eq:T-epsilon-explicit.
#
# This block does not assume epsilon_Phi(W)=1 and does not use the constant
# (-1) coefficient table below.  It starts from
#
#   a_2=-x_4 action a_1,  a_3=-x_2 action a_1,
#   a_4=-x_3 action a_1,
#
# and applies eq:YD-projective-action.  For a basis vector
# a_j=-x_s action a_1, the rack finds the unique carrier x_k for which
# x_k and x_s^2 have the same action on the degree x_1.  Consequently
# c=x_k^-1*x_s^2 lies in C_G(x_1), and
#
# x_s action a_j
#   = Phi_x1(x_s,x_s) / Phi_x1(x_k,c) * sigma(c) * a_t,
#
# where a_t=-x_k action a_1.  The two minus signs are produced from the
# two basis definitions; they are not inserted into the three targets.
#############################################################################

Print("\nPart I: reconstruct eq:T-epsilon-explicit without assuming epsilon_Phi(W)=1.\n");

# The carrier and definition sign for a_1,a_2,a_3,a_4.  The carrier of
# a_1 is recorded only to make the degree list uniform.
T_BASIS_CARRIER := [1, 4, 2, 3];;
T_BASIS_DEFINITION_SIGN := [1, -1, -1, -1];;

T_XName := function(i)
    return Concatenation("x", String(i));
end;;

T_TailName := function(k, s)
    return Concatenation(T_XName(k), "^-1*", T_XName(s), "^2");
end;;

T_LocalPhiName := function(left, right)
    return Concatenation("Phi_x1(", left, ",", right, ")");
end;;

T_SigmaName := function(word)
    return Concatenation("sigma(", word, ")");
end;;

T_BasisDegree := function(i)
    if i = 1 then
        return 1;
    fi;
    return T_RackOp(T_BASIS_CARRIER[i], 1);
end;;

T_DeriveSquareAction := function(s)
    local candidates, coefficient, j, k, tail, target, targetDegree;
    j := Position(T_BASIS_CARRIER, s);
    if j = fail or j = 1 then
        T_Fail("No nontrivial basis definition has the requested carrier");
    fi;

    targetDegree := T_RackOp(s, T_RackOp(s, 1));
    candidates := Filtered([2 .. 4],
        k -> T_RackOp(k, 1) = targetDegree);
    if Length(candidates) <> 1 then
        T_Fail("The rack did not determine a unique carrier for x_s^2 action a_1");
    fi;
    k := candidates[1];
    target := Position(T_BASIS_CARRIER, k);
    if target = fail then
        T_Fail("The derived carrier has no compatible basis vector");
    fi;

    # Since x_k and x_s^2 carry degree x_1 to the same degree, the
    # formally defined tail x_k^-1*x_s^2 fixes x_1.
    tail := T_TailName(k, s);
    coefficient := T_LaurentConstant(
        T_BASIS_DEFINITION_SIGN[j] /
        T_BASIS_DEFINITION_SIGN[target]);
    coefficient := T_LaurentMul(coefficient,
        T_LaurentGenerator(T_LocalPhiName(T_XName(s), T_XName(s))));
    coefficient := T_LaurentMul(coefficient,
        T_LaurentInv(T_LaurentGenerator(
            T_LocalPhiName(T_XName(k), tail))));
    coefficient := T_LaurentMul(coefficient,
        T_LaurentGenerator(T_SigmaName(tail)));

    return rec(
        actor := s,
        source := j,
        target := target,
        tail := tail,
        scalar := coefficient
    );
end;;

T_epsilon_actions := List([2, 4, 3], T_DeriveSquareAction);;
T_Check(
    List(T_epsilon_actions,
        a -> [a.actor, a.source, a.target, a.tail]) = [
        [2, 3, 4, "x3^-1*x2^2"],
        [4, 2, 3, "x2^-1*x4^2"],
        [3, 4, 2, "x4^-1*x3^2"]
    ],
    "eq:T-epsilon-explicit: three coefficients from eq:YD-projective-action derived from the rack and basis definitions");;

# These are comparison targets only.  They are not consulted by
# T_DeriveSquareAction or by the braid accumulator below.
T_epsilon_action_tex_targets := [
    T_LaurentNormalize([1, [
        ["Phi_x1(x2,x2)", 1],
        ["Phi_x1(x3,x3^-1*x2^2)", -1],
        ["sigma(x3^-1*x2^2)", 1]
    ]]),
    T_LaurentNormalize([1, [
        ["Phi_x1(x4,x4)", 1],
        ["Phi_x1(x2,x2^-1*x4^2)", -1],
        ["sigma(x2^-1*x4^2)", 1]
    ]]),
    T_LaurentNormalize([1, [
        ["Phi_x1(x3,x3)", 1],
        ["Phi_x1(x4,x4^-1*x3^2)", -1],
        ["sigma(x4^-1*x3^2)", 1]
    ]])
];;
T_Check(
    List(T_epsilon_actions, a -> a.scalar) =
        T_epsilon_action_tex_targets,
    "eq:T-epsilon-explicit: all three generated action coefficients match their displayed formulas");;

T_BraidBasisPair := function(pair, scalar)
    local action, leftDegree, rightCarrier;
    leftDegree := T_BasisDegree(pair[1]);
    rightCarrier := T_BASIS_CARRIER[pair[2]];
    if leftDegree <> rightCarrier then
        T_Fail("The right basis vector is not defined by the degree acting in this braid");
    fi;
    action := T_DeriveSquareAction(leftDegree);
    if action.source <> pair[2] then
        T_Fail("The action derived from eq:YD-projective-action has the wrong source basis vector");
    fi;
    return rec(
        pair := [action.target, pair[1]],
        scalar := T_LaurentMul(scalar, action.scalar)
    );
end;;

T_epsilon_pair := [2, 3];;
T_epsilon_scalar := T_LaurentOne;;
T_epsilon_states := [ShallowCopy(T_epsilon_pair)];;
for T_n in [1 .. 3] do
    T_result := T_BraidBasisPair(T_epsilon_pair, T_epsilon_scalar);;
    T_epsilon_pair := T_result.pair;;
    T_epsilon_scalar := T_result.scalar;;
    Add(T_epsilon_states, ShallowCopy(T_epsilon_pair));
od;
T_Check(T_epsilon_states = [[2,3], [4,2], [3,4], [2,3]],
    "eq:T-epsilon-explicit: three actual braids close on a_2 tensor a_3");;

# From c^3=sigma(x_1)^-1*epsilon_Phi(W), one has
# epsilon_Phi(W)=sigma(x_1)*c^3.  Here sigma(x_1)=-1 is an explicit
# hypothesis of the displayed coordinate calculation.
T_epsilon_derived := T_LaurentMul(
    T_LaurentConstant(-1), T_epsilon_scalar);;

# Independent transcription of the right-hand side of
# eq:T-epsilon-explicit.  In particular, its sign and all three sigma
# factors are comparison data, whereas the left-hand side above was
# generated by the rack, the basis definitions, and three braid steps.
T_epsilon_tex_rhs := T_LaurentNormalize([-1, [
    ["Phi_x1(x2,x2)", 1],
    ["Phi_x1(x4,x4)", 1],
    ["Phi_x1(x3,x3)", 1],
    ["Phi_x1(x3,x3^-1*x2^2)", -1],
    ["Phi_x1(x2,x2^-1*x4^2)", -1],
    ["Phi_x1(x4,x4^-1*x3^2)", -1],
    ["sigma(x3^-1*x2^2)", 1],
    ["sigma(x2^-1*x4^2)", 1],
    ["sigma(x4^-1*x3^2)", 1]
]]);;
T_Check(T_LaurentEqual(T_epsilon_derived, T_epsilon_tex_rhs),
    "eq:T-epsilon-explicit: exact Laurent identity, including the total minus sign and three sigma factors");;

#############################################################################
# Part II: the pre-existing checks below use the compatible constant (-1)
# braiding table and therefore assume epsilon_Phi(W)=1.
#############################################################################

Print("\nPart II: existing coefficient checks under epsilon_Phi(W)=1.\n");

T_BRAID_COEFFS := List([1 .. 4],
    i -> List([1 .. 4], j -> -1));;
T_ok := true;;
for T_i in [1 .. 4] do
    for T_j in [1 .. 4] do
        for T_k in [1 .. 4] do
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
    "lem:T-tetrahedral-braiding: all 64 constant-(-1) braid instances and assumed epsilon_Phi(W)=1");;

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
# Compatible-basis coefficient recursion eq:T-adjoint-recursion.
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
# Lemma D.1, lem:T-Y2-homogeneous-components, and eq:T-Y2-basis:
# u_i and the complete internal phi_2 calculation, modulo f(p).
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
    "lem:T-Y2-homogeneous-components: complete phi_2 calculation modulo p^2-p+1");;

#############################################################################
# Lemma D.2, lem:T-Y3-homogeneous-components, and eq:T-y3-generator:
# y and the complete internal phi_3 calculation, modulo f(p).
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
    "lem:T-Y3-homogeneous-components: complete phi_3 calculation modulo p^2-p+1");;

#############################################################################
# Lemma D.3, lem:T-Y4-homogeneous-components:
# exact, unreduced internal identity phi_4(w_1 tensor y)=f(p) Z_1.
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
    "lem:T-Y4-homogeneous-components: fifteen pairwise distinct words in the internal Z_1 target");;
T_Check(T_VecEqual(T_actual4, T_expected4),
    "lem:T-Y4-homogeneous-components: exact factorization in unreduced Z[p]");;
T_Check(T_VecIsZero(T_VecReduceF(T_actual4)),
    "lem:T-Y4-homogeneous-components: vanishing after exact factorization is reduced modulo f");;

#############################################################################
# eq:T-R2-monodromy-cycles and eq:T-R2-mixed-monodromy
# derive the four monodromy cycles from adjacent braids.
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
