#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$ROOT_DIR"

MPI_RANKS=${MPI_RANKS:-1}
CASE=${1:-both}

if [[ ! -x ./dpm.sc_mpi || ! -f kernel/kernel.pz ]]; then
  echo "Missing ./dpm.sc_mpi or kernel/kernel.pz. Build them with:" >&2
  echo "  make clean MODE=mpi_sc3s" >&2
  echo "  make mpi_sc3s PRECISION=float KERNEL_VERSION=4th PZC_TARGET_ARCH=sc3s" >&2
  exit 1
fi

mkdir -p data/output

run_case() {
  local name=$1
  local input=$2
  local output="data/output/pezy_sc3s_${name}_smoke.out"

  echo "Running MPI + PEZY-SC3s ${name} smoke test with ${MPI_RANKS} rank(s)..."
  mpirun -n "$MPI_RANKS" ./dpm.sc_mpi < "$input" > "$output"

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
