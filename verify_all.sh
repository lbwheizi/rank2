#!/bin/sh
set -eu

REPOSITORY_ROOT=$(CDPATH= cd -P "$(dirname "$0")" && pwd)
cd "$REPOSITORY_ROOT"

GAP_BIN_INPUT=${GAP_BIN:-gap}
GAP_ROOT_INPUT=${GAP_ROOT:-}

if ! GAP_EXECUTABLE=$(command -v "$GAP_BIN_INPUT"); then
    printf '%s\n' "ERROR: GAP executable not found: $GAP_BIN_INPUT" >&2
    exit 1
fi

case "$GAP_EXECUTABLE" in
    /*) ;;
    *)
        GAP_EXECUTABLE=$(CDPATH= cd -P "$(dirname "$GAP_EXECUTABLE")" && pwd)/$(basename "$GAP_EXECUTABLE")
        ;;
esac

if [ -n "$GAP_ROOT_INPUT" ]; then
    if [ ! -d "$GAP_ROOT_INPUT" ]; then
        printf '%s\n' "ERROR: GAP_ROOT is not a directory: $GAP_ROOT_INPUT" >&2
        exit 1
    fi
    GAP_ROOT_DIRECTORY=$(CDPATH= cd -P "$GAP_ROOT_INPUT" && pwd)
else
    GAP_ROOT_DIRECTORY=
fi

run_gap() {
    script=$1
    if [ -n "$GAP_ROOT_DIRECTORY" ]; then
        "$GAP_EXECUTABLE" --bare -A -r -q --nointeract \
            -l "$GAP_ROOT_DIRECTORY" "$script"
    else
        "$GAP_EXECUTABLE" --bare -A -r -q --nointeract "$script"
    fi
}

run_certificate() {
    directory=$1
    script=$2
    recorded=$3
    temporary=$(mktemp "${TMPDIR:-/tmp}/rank2-gap-output.XXXXXX")
    trap 'rm -f "$temporary"' EXIT HUP INT TERM

    (cd "$directory" && sha256sum -c CHECKSUMS.sha256)
    (cd "$directory" && run_gap "$script") > "$temporary"
    if ! cmp "$temporary" "$directory/$recorded"; then
        printf '%s\n' "ERROR: output differs in $directory" >&2
        exit 1
    fi

    rm -f "$temporary"
    trap - EXIT HUP INT TERM
    printf '%s\n' "PASS: $directory"
}

sha256sum -c SHA256SUMS

run_certificate certificates/Gamma_2/delta_reduction lemma_B1_first_coefficient_reduction.g lemma_B1_first_coefficient_reduction.out
run_certificate certificates/Gamma_2/omegaV lemma_B4_Omega_V_equals_one.g lemma_B4_Omega_V_equals_one.out
run_certificate certificates/Gamma_2/g2_propagation lemma_B12_K3_and_d4_cocycle_reductions.g lemma_B12_K3_and_d4_cocycle_reductions.out
run_certificate certificates/Gamma_2/g2_pure_parameters/s_phi_equals_theta_squared lemma_B14_s_phi_equals_theta_squared.g lemma_B14_s_phi_equals_theta_squared.out
run_certificate certificates/Gamma_2/g2_pure_parameters/theta_r2_invariance lemma_B15_theta_after_second_reflection.g lemma_B15_theta_after_second_reflection.out
run_certificate certificates/Gamma_2/theta_symmetry lemma_B19_theta_symmetry.g lemma_B19_theta_symmetry.out
run_certificate certificates/Gamma_2/order16_G2 remark_B21_order16_G2_verifier.g remark_B21_order16_G2_verifier.out
run_certificate certificates/Gamma_2/g2_pure_parameters/auxiliary/s_phi_r2_invariance auxiliary_s_phi_R2_invariance.g auxiliary_s_phi_R2_invariance.out
run_certificate certificates/Gamma_2/g2_pure_parameters/auxiliary/s_phi_square auxiliary_s_phi_square.g auxiliary_s_phi_square.out

run_certificate certificates/Gamma_3 verify_all_gamma3.g EXPECTED_OUTPUT.txt

run_certificate certificates/Gamma_4/xi3 gamma4_xi3_certificate.g gamma4_xi3_certificate.out
run_certificate certificates/Gamma_4/local_delta gamma4_local_delta_certificate.g gamma4_local_delta_certificate.out
run_certificate certificates/Gamma_4/delta_transport gamma4_delta_transport_certificate.g gamma4_delta_transport_certificate.out
run_certificate certificates/Gamma_4/terminal_character gamma4_terminal_character_certificate.g gamma4_terminal_character_certificate.out
run_certificate certificates/Gamma_4/reflected_y2 gamma4_reflected_y2_certificate.g gamma4_reflected_y2_certificate.out

run_certificate certificates/T verify_T_case.g verify_T_case.out

printf '%s\n' 'PASS: all exact certificate suites'
