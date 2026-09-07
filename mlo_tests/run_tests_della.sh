#!/bin/bash
# Verification suite for pair_style lj/cut/mlo on Della.
# Usage: bash mlo_tests/run_tests_della.sh
# Run from the lammps-stable_22Jul2025 directory on the head node.
# Uses at most two local MPI ranks; no Slurm allocation is needed.
set -o pipefail

module purge || exit 1
module load gcc-toolset/14 || exit 1
module load aocc/5.0.0 || exit 1
module load aocl/aocc/5.0.0 || exit 1
module load openmpi/aocc-5.0.0/4.1.6 || exit 1

# Keep each MPI rank single-threaded so tests use at most two CPUs.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export BLIS_NUM_THREADS=1

LMP_BIN=${LMP_BIN:-$HOME/.local/bin/lmp_d9_double_aocc}
[ -x "$LMP_BIN" ] || { echo "no binary at $LMP_BIN"; exit 1; }
T=mlo_tests
pass=0; fail=0
ok()   { printf '  \033[32mPASS\033[0m  %s\n' "$1"; pass=$((pass+1)); }
bad()  { printf '  \033[31mFAIL\033[0m  %s   %s\n' "$1" "$2"; fail=$((fail+1)); }

# Return only validated test records. Show the full MPI/LAMMPS output on
# failure so launch errors cannot disappear into an empty grep result.
run_data() { # $1 ranks  $2 record type  $3 input file  [LAMMPS arguments]
  local np="$1" record="$2" input="$3" output status
  shift 3
  local command=(mpirun --host localhost:2 -np "$np" "$LMP_BIN" -log none -in "$input" "$@")
  if output=$("${command[@]}" 2>&1); then
    :
  else
    status=$?
    printf 'ERROR: MPI/LAMMPS exited with status %s.\nCommand:' "$status" >&2
    printf ' %q' "${command[@]}" >&2
    printf '\n%s\n' "$output" >&2
    return "$status"
  fi

  # These inputs print 14 DATA rows, 11 FD rows, or one S0 row.
  # Reject missing, incomplete, or nonnumeric results before comparison.
  if ! printf '%s\n' "$output" | awk -v record="$record" '
    function number(s) {
      return s ~ /^[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][+-]?[0-9]+)?$/
    }
    function field(s, key, parts) {
      return split(s, parts, "=") == 2 && parts[1] == key && number(parts[2])
    }
    BEGIN { expected = (record == "DATA" ? 14 : (record == "FD" ? 11 : 1)) }
    $1 == record {
      rows++
      if (record == "DATA") {
        if (NF != 6) invalid = 1
        for (i = 2; i <= NF; i++) if (!number($i)) invalid = 1
      } else if (record == "FD") {
        if (NF != 7 || !field($2, "r") || !field($3, "evdwl") ||
            !field($4, "maxerr_x") || !field($5, "maxerr_y") ||
            !field($6, "maxerr_z") || !field($7, "max_fz")) invalid = 1
      } else if (record == "S0") {
        if (NF != 5 || !field($2, "evdwl") || !field($3, "sfx") ||
            !field($4, "sfy") || !field($5, "sfz")) invalid = 1
      } else invalid = 1
    }
    END { exit (invalid || rows != expected) }
  '; then
    printf 'ERROR: Missing, incomplete, or invalid %s data.\nCommand:' "$record" >&2
    printf ' %q' "${command[@]}" >&2
    printf '\n%s\n' "$output" >&2
    return 1
  fi
  printf '%s\n' "$output" | awk -v record="$record" '$1 == record'
}

# ---- 1. collapse onto stock lj/cut in the two exact limits ----------------
# lambda_max = 0 must be exact WCA; lambda_max = 1 with Theta = 1 must be plain lj/cut.
cmp_limit () { # $1 label  $2 mlo-lambda  $3 zpos  $4 reference-style
  local a b d
  a=$(run_data 2 DATA "$T/collapse.in" -var style mlo -var lam "$2" -var zpos "$3") || exit 1
  b=$(run_data 2 DATA "$T/collapse.in" -var style "$4" -var lam 1.0 -var zpos "$3") || exit 1
  d=$(paste <(echo "$a") <(echo "$b") | awk '{if ($2 != $8) invalid=1; du=$3-$9; dfx=$4-$10;
        if(du<0)du=-du; if(dfx<0)dfx=-dfx; m=(du>dfx?du:dfx); if(m>w)w=m} END{if(invalid) exit 1; printf "%.3e", w+0}') || { bad "$1" "mismatched sample positions"; return; }
  [ "$(echo "$d" | awk '{print ($1==0)?1:0}')" = 1 ] && ok "$1 (exact)" || bad "$1" "worst diff $d"
}
cmp_limit "lambda_max=0 reproduces WCA bitwise"          0.0 0.0 wca
cmp_limit "lambda_max=1, Theta=1 reproduces lj/cut bitwise" 1.0 3.0 full

# ---- 2. finite-difference gradient test ----------------------------------
# MUST run at alpha ~ 5. At alpha = 50 the z-forces are numerically zero away from the
# crossover and the test is vacuous -- it would pass with lambda_max dropped from them.
out=$(run_data 2 FD "$T/numdiff.in") || exit 1
worst=$(echo "$out" | awk '{for(i=1;i<=NF;i++){split($i,p,"=");
          if(p[1]=="maxerr_z"&&p[2]+0>w)w=p[2]+0}} END{printf "%.3e", w}')
fzmax=$(echo "$out" | awk '{for(i=1;i<=NF;i++){split($i,p,"=");
          if(p[1]=="max_fz"&&p[2]+0>w)w=p[2]+0}} END{printf "%.3e", w}')
if awk -v worst="$worst" -v fzmax="$fzmax" 'BEGIN{exit !(worst < 1e-6 && fzmax > 0.1)}'; then
  ok "z-force = -dU/dz to $worst  (max|fz| = $fzmax, so not vacuous)"
else bad "finite-difference z-force" "maxerr_z=$worst max|fz|=$fzmax"; fi

# ---- 3. asymmetry invariant ----------------------------------------------
# x,y are cyclic so their force sums must vanish; z is not, so its sum must NOT.
s=$(run_data 2 S0 "$T/mpi_consistency.in" -var nt on) || exit 1
sfx=$(echo "$s"|sed 's/.*sfx=\([^ ]*\).*/\1/'); sfz=$(echo "$s"|sed 's/.*sfz=\([^ ]*\).*/\1/')
if awk -v x="$sfx" -v z="$sfz" 'BEGIN{if(x<0)x=-x; if(z<0)z=-z; exit !(x<1e-12 && z>1e-6)}'; then
  ok "sum fx = $sfx (zero), sum fz = $sfz (nonzero, as the model requires)"
else bad "force-sum invariants" "sfx=$sfx sfz=$sfz"; fi

# ---- 4. MPI / newton consistency -----------------------------------------
# Exercises reverse communication of the ASYMMETRIC z-forces.
ref=""; consistent=1
for nt in on off; do for np in 1 2; do
  s=$(run_data "$np" S0 "$T/mpi_consistency.in" -var nt "$nt") || exit 1
  v=$(echo "$s" | sed 's/.*sfz=\([^ ]*\).*/\1/')
  [ -z "$ref" ] && ref="$v"
  awk -v v="$v" -v ref="$ref" 'BEGIN{d=v-ref; if(d<0)d=-d;
      scale=ref; if(scale<0)scale=-scale; exit !(scale>0 && d/scale < 1e-12)}' || consistent=0
done; done
[ $consistent = 1 ] && ok "sum fz identical across np=1,2 x newton on/off" \
                    || bad "MPI/newton consistency" "sum fz varies"

printf '\n  %d passed, %d failed\n' "$pass" "$fail"
if [ "$fail" -eq 0 ]; then
  module purge || exit 1
fi
exit $((fail > 0))
