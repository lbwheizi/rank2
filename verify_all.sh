#!/bin/sh
set -eu

GAP_BIN="${GAP_BIN:-gap}"

run_certificate() {
    directory="$1"
    script="$2"
    recorded="$3"
    temporary="$(mktemp)"
    trap 'rm -f "$temporary"' EXIT HUP INT TERM

    (cd "$directory" && sha256sum -c CHECKSUMS.sha256)
    (cd "$directory" && "$GAP_BIN" -A -q --nointeract "$script") > "$temporary"
    cmp "$temporary" "$directory/$recorded"
    rm -f "$temporary"
    trap - EXIT HUP INT TERM
    printf '%s\n' "PASS: $directory"
}

run_certificate certificates/delta_reduction delta_certificate.g delta_certificate.out
run_certificate certificates/omegaV omegaV_certificate.g omegaV_certificate.out
run_certificate certificates/theta_symmetry certificate.g recorded_output.txt
run_certificate certificates/g2_propagation certificate.g recorded_output.txt
run_certificate certificates/g2_pure_parameters/xi_equals_theta_squared certificate.g recorded_output.txt
run_certificate certificates/g2_pure_parameters/theta_r2_invariance certificate.g recorded_output.txt
run_certificate certificates/g2_pure_parameters/xi_r2_invariance certificate.g recorded_output.txt
run_certificate certificates/g2_pure_parameters/xi_square certificate.g recorded_output.txt
run_certificate certificates/order16_G2 certificate.g recorded_output.txt

printf '%s\n' 'PASS: all exact certificate suites'
