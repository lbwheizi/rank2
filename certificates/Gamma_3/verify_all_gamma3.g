G3_Fail := function(arg)
    CallFuncList(PrintTo,
        Concatenation(["*errout*"], arg, ["\n"]));
    ErrorNoReturn("Certificate verification failed.");
end;

Print("Gamma_3 exact certificate suite (GAP core only)\n");
Print("================================================\n");

G3_GAMMA332_OK := false;
G3_GAMMA332_COEFFICIENT_DICTIONARY_OK := false;
G3_GAMMA331_DIM2_OK := false;
G3_GAMMA331_DIM1_OK := false;
G3_S3_OK := false;

Read("verify_gamma3_32.g");
if not G3_GAMMA332_OK then
    G3_Fail("The (3,2) verifier did not report success");
fi;
if not G3_GAMMA332_COEFFICIENT_DICTIONARY_OK then
    G3_Fail("The (3,2) coefficient/dictionary verifier did not report success");
fi;

Read("verify_gamma3_31_dim1.g");
if not G3_GAMMA331_DIM2_OK then
    G3_Fail("The (3,1)_2 verifier did not report success");
fi;
if not G3_GAMMA331_DIM1_OK then
    G3_Fail("The (3,1)_1 verifier did not report success");
fi;

Read("verify_s3_cocycle_mod3.g");
if not G3_S3_OK then
    G3_Fail("The S3 verifier did not report success");
fi;

if not (G3_GAMMA332_OK and G3_GAMMA332_COEFFICIENT_DICTIONARY_OK
        and G3_GAMMA331_DIM2_OK
        and G3_GAMMA331_DIM1_OK and G3_S3_OK) then
    G3_Fail("The Gamma_3 certificate suite is incomplete");
fi;

Print("================================================\n");
Print("ALL GAMMA_3 CERTIFICATE CHECKS PASSED.\n");
