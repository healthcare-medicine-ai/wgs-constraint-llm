#!/bin/bash
# Submit every fidelity gate and report one verdict.
#
# Until this existed the gates were run by hand, one at a time, which meant
# "the analysis is validated" was a claim about nine separate log files someone
# had read at nine different times. It also let a gate sit red without anyone
# noticing. This submits all of them and prints a single table.
#
#   bash jobs/validate_all.sh            # submit everything, print job ids
#   bash jobs/validate_all.sh --status   # summarise the last submitted run
#
# Ordering: gate_01_definitive consumes the _REGATE input that gate_01_refactor
# produces, so it waits on it. Everything else is independent and runs in
# parallel, bounded by the partition.
#
# Gates that demand bit-identity pin --constraint=CPU_MNF:INTEL, because np.log
# differs in the last bit between AVX-512 and AVX2 and the reference artifacts
# were built on Intel. That pin lives in the individual job scripts.
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1
STATE=logs/validate_all.jobs

if [ "${1:-}" = "--status" ]; then
    [ -f "$STATE" ] || { echo "no run recorded in $STATE"; exit 1; }
    printf "%-26s %-10s %-12s %-10s %s\n" GATE JOBID STATE EXIT ELAPSED
    fail=0
    while read -r name jobid; do
        line=$(sacct -j "$jobid" --format=State%14,ExitCode%8,Elapsed%10 -n 2>/dev/null | head -1)
        st=$(echo "$line" | awk '{print $1}')
        ex=$(echo "$line" | awk '{print $2}')
        el=$(echo "$line" | awk '{print $3}')
        printf "%-26s %-10s %-12s %-10s %s\n" "$name" "$jobid" "${st:-?}" "${ex:-?}" "${el:-?}"
        case "$st" in
            COMPLETED) ;;
            RUNNING|PENDING) fail=2 ;;
            *) fail=1 ;;
        esac
    done < "$STATE"
    echo
    case $fail in
        0) echo "ALL GATES PASS" ;;
        2) echo "still running -- rerun --status later" ;;
        *) echo "AT LEAST ONE GATE FAILED" ;;
    esac
    exit $fail
fi

mkdir -p logs
: > "$STATE"
submit() {                      # submit <name> <script> [dependency-jobid]
    local name=$1 script=$2 dep=${3:-}
    local args=(--parsable)
    [ -n "$dep" ] && args+=(--dependency=afterok:"$dep")
    local jobid
    jobid=$(sbatch "${args[@]}" "$script")
    echo "$name $jobid" >> "$STATE"
    printf "  %-26s %s%s\n" "$name" "$jobid" "${dep:+  (after $dep)}" >&2
    echo "$jobid"
}

echo "submitting fidelity gates:"
REFACTOR=$(submit gate_01_refactor    jobs/gate_01_refactor.sbatch    | tail -1)
submit gate_01_definitive jobs/gate_01_definitive.sbatch "$REFACTOR" >/dev/null
submit gate_01_fixed_path jobs/gate_01_fixed_path.sbatch >/dev/null
submit gate_07_ablations  jobs/gate_07_ablations.sbatch  >/dev/null
submit gate_08_figure4    jobs/gate_08_figure4.sbatch    >/dev/null
submit gate_09_figure3    jobs/gate_09_figure3.sbatch    >/dev/null
submit gate_10_figure2    jobs/gate_10_figure2.sbatch    >/dev/null
submit gate_11_scz        jobs/gate_11_scz.sbatch        >/dev/null
submit gate_12_figure1    jobs/gate_12_figure1.sbatch    >/dev/null

echo
echo "recorded in $STATE"
echo "check with: bash jobs/validate_all.sh --status"
