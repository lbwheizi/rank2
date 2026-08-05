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

run_certificate certificates/delta_reduction delta_certificate.g delta_certificate.out
run_certificate certificates/omegaV omegaV_certificate.g omegaV_certificate.out
run_certificate certificates/theta_symmetry certificate.g recorded_output.txt
run_certificate certificates/g2_propagation certificate.g recorded_output.txt
run_certificate certificates/g2_pure_parameters/s_phi_equals_theta_squared certificate.g recorded_output.txt
run_certificate certificates/g2_pure_parameters/theta_r2_invariance certificate.g recorded_output.txt
run_certificate certificates/g2_pure_parameters/auxiliary/s_phi_r2_invariance certificate.g recorded_output.txt
run_certificate certificates/g2_pure_parameters/auxiliary/s_phi_square certificate.g recorded_output.txt
run_certificate certificates/order16_G2 certificate.g recorded_output.txt

printf '%s\n' 'PASS: all exact certificate suites'
