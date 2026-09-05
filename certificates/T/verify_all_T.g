# Run the existing coefficient calculation and all five additional T suites.
# Start in the directory containing these .g files. GAP core only.
# In an interactive session a successful Read returns to the GAP prompt.

for T_GAP_SUITE in ["verify_T_case.g", "verify_T_cohomology_reduction.g",
        "verify_T_sl2_cocycle.g", "verify_T_actual_adjoint.g",
        "verify_T_nichols_iterative.g", "verify_T_order2_reflections.g"] do
    TLastCheckPassed:=false;
    Read(T_GAP_SUITE);
    if not TLastCheckPassed then
        PrintTo("*errout*", "FAIL: incomplete T suite: ", T_GAP_SUITE, "\n");
        ErrorNoReturn("T verification did not complete.");
    fi;
od;
Print("PASS: all six T GAP suites\n");
