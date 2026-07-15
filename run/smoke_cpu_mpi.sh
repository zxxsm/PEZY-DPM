#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$ROOT_DIR"

MPI_RANKS=${MPI_RANKS:-2}
CASE=${1:-both}

if [[ ! -x ./dpm.mpi ]]; then
  echo "Missing ./dpm.mpi. Build it with:" >&2
  echo "  make clean MODE=mpi && make mpi PRECISION=float" >&2
  exit 1
fi

mkdir -p data/output

run_case() {
  local name=$1
  local input=$2
  local output="data/output/cpu_mpi_${name}_smoke.out"

  echo "Running CPU MPI ${name} smoke test with ${MPI_RANKS} rank(s)..."
  mpirun -n "$MPI_RANKS" ./dpm.mpi < "$input" > "$output"

  grep -A1 "No of histories simulated:" "$output" | tail -n 1 \
    | grep -Eq '(^|[[:space:]])10000([[:space:]]|$)'
  grep -q "Have a nice day." "$output"
  echo "PASS: $output"
}

case "$CASE" in
  electron)
    run_case electron data/input/electron_20MeV_smoke1e4.in
    ;;
  photon)
    run_case photon data/input/photon_6MeV_smoke1e4.in
    ;;
  both)
    run_case electron data/input/electron_20MeV_smoke1e4.in
    run_case photon data/input/photon_6MeV_smoke1e4.in
    ;;
  *)
    echo "Usage: $0 [electron|photon|both]" >&2
    exit 2
    ;;
esac
