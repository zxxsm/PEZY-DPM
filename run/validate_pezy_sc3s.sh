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
  local histories=$2
  local input=$3
  local reference=$4
  local threshold=$5
  local output="data/output/pezy_sc3s_${name}_validation.out"

  echo "Running the manuscript ${name} validation case with ${MPI_RANKS} PEZY-SC3s device(s)..."
  mpirun -n "$MPI_RANKS" ./dpm.sc_mpi < "$input" > "$output"

  grep -A1 "No of histories simulated:" "$output" | tail -n 1 \
    | grep -Eq "(^|[[:space:]])${histories}([[:space:]]|$)"
  grep -q "Have a nice day." "$output"
  python3 tools/validate_cax.py \
    --reference "$reference" \
    --candidate "$output" \
    --mean-threshold "$threshold"
}

case "$CASE" in
  electron)
    run_case electron 10000000 \
      data/input/electron_20MeV_validation1e7.in \
      data/reference/electron_20MeV_original_dpm_cax.csv 0.0950
    ;;
  photon)
    run_case photon 50000000 \
      data/input/photon_6MeV_validation5e7.in \
      data/reference/photon_6MeV_original_dpm_cax.csv 0.2712
    ;;
  both)
    "$0" electron
    "$0" photon
    ;;
  *)
    echo "Usage: $0 [electron|photon|both]" >&2
    exit 2
    ;;
esac
