# Test report on server100

Test directory on server100:

```text
/vol1/home/zhuxiaoxiong/PEZY_DPM_clean_codex_20260617/run_20260617_1120
```

Environment detected on server100:

- Host: `server100`
- GNU C++: `/usr/bin/g++`
- GNU Fortran: `/usr/bin/gfortran`
- MPI: `/usr/local/mpich/bin`
- PZSDK: `/opt/pezy/pzsdk`

Smoke-test input:

- Generated from `data/input/dpm.in`.
- Number of histories reduced to `10000` for fast compile-and-run validation.

## Results

| Test step | Result |
| --- | --- |
| `make openmp PRECISION=float` | PASS |
| `./dpm.omp` OpenMP float smoke run | PASS |
| `make openmp PRECISION=double` | PASS |
| `./dpm.omp` OpenMP double smoke run | PASS |
| `make mpi PRECISION=float` | PASS |
| `mpirun -n 2 ./dpm.mpi` MPI float smoke run | PASS |
| `make mpi PRECISION=double` | PASS |
| `mpirun -n 2 ./dpm.mpi` MPI double smoke run | PASS |
| `make mpi_sc3s PRECISION=float KERNEL_VERSION=8th PZC_TARGET_ARCH=sc3s` | PASS |
| `mpirun -n 1 ./dpm.sc_mpi` SC3s float 8th smoke run | PASS |
| `make mpi_sc3s PRECISION=float KERNEL_VERSION=4th PZC_TARGET_ARCH=sc3s` | PASS |
| `mpirun -n 1 ./dpm.sc_mpi` SC3s float 4th smoke run | PASS |
| `make mpi_sc3s PRECISION=double KERNEL_VERSION=8th PZC_TARGET_ARCH=sc3s` | PASS after fixing `atan2` literal type ambiguity |
| `mpirun -n 1 ./dpm.sc_mpi` SC3s double 8th smoke run | PASS |
| `make mpi_sc3s PRECISION=double KERNEL_VERSION=4th PZC_TARGET_ARCH=sc3s` | PASS after fixing `atan2` literal type ambiguity |
| `mpirun -n 1 ./dpm.sc_mpi` SC3s double 4th smoke run | PASS |

The fix applied after the first test pass was limited to the PZC kernels: three calls to `atan2(1.0E0f, ...)` in each of `pzc/kernel_8th.pzc` and `pzc/kernel_4th.pzc` were changed to `atan2((UREAL)1.0, ...)`, so the literal type follows the selected `PRECISION`.

Full logs and smoke outputs were synchronized locally to:

```text
clean_test_logs_outputs_100_20260617/
```

## Retest after source-tree reorganization

After moving C/C++ and Fortran source files into `src/`, the Makefile was updated to use `vpath` and relative include paths for `compile/cpu.mk` and `compile/default_pzcl_host.mk`. The reorganized tree was then retested on server100 in:

```text
/vol1/home/zhuxiaoxiong/PEZY_DPM_clean_codex_20260617/run_20260617_src
```

All build and smoke-run combinations passed again:

| Test step | Result |
| --- | --- |
| `make openmp PRECISION=float` | PASS |
| `./dpm.omp` OpenMP float smoke run | PASS |
| `make openmp PRECISION=double` | PASS |
| `./dpm.omp` OpenMP double smoke run | PASS |
| `make mpi PRECISION=float` | PASS |
| `mpirun -n 2 ./dpm.mpi` MPI float smoke run | PASS |
| `make mpi PRECISION=double` | PASS |
| `mpirun -n 2 ./dpm.mpi` MPI double smoke run | PASS |
| `make mpi_sc3s PRECISION=float KERNEL_VERSION=8th PZC_TARGET_ARCH=sc3s` | PASS |
| `mpirun -n 1 ./dpm.sc_mpi` SC3s float 8th smoke run | PASS |
| `make mpi_sc3s PRECISION=float KERNEL_VERSION=4th PZC_TARGET_ARCH=sc3s` | PASS |
| `mpirun -n 1 ./dpm.sc_mpi` SC3s float 4th smoke run | PASS |
| `make mpi_sc3s PRECISION=double KERNEL_VERSION=8th PZC_TARGET_ARCH=sc3s` | PASS |
| `mpirun -n 1 ./dpm.sc_mpi` SC3s double 8th smoke run | PASS |
| `make mpi_sc3s PRECISION=double KERNEL_VERSION=4th PZC_TARGET_ARCH=sc3s` | PASS |
| `mpirun -n 1 ./dpm.sc_mpi` SC3s double 4th smoke run | PASS |

The second-round logs and smoke outputs were synchronized locally to:

```text
src_reorg_test_logs_outputs_100_20260617/
```

## MPI + two PEZY-SC3s smoke test

After adding the named benchmark inputs `electron_20MeV_2p5e8.in` and
`photon_6MeV_2p5e8.in`, the package was retested on server100 with two MPI ranks
and two PEZY-SC3s devices. The test directory was:

```text
/vol1/home/zhuxiaoxiong/PEZY_DPM_github_codex_20260618_mpi2/pkg
```

The test used the final 4-thread PEZY kernel variant:

```bash
make mpi_sc3s PRECISION=float KERNEL_VERSION=4th PZC_TARGET_ARCH=sc3s
```

For a fast smoke run, the version-controlled inputs
`data/input/electron_20MeV_smoke1e4.in` and
`data/input/photon_6MeV_smoke1e4.in` were used directly. Each input specifies
`10000` source histories. Both runs used:

```bash
mpirun -n 2 ./dpm.sc_mpi
```

| Test step | Result |
| --- | --- |
| `make mpi_sc3s PRECISION=float KERNEL_VERSION=4th PZC_TARGET_ARCH=sc3s` | PASS |
| `mpirun -n 2 ./dpm.sc_mpi < data/input/electron_20MeV_smoke1e4.in` | PASS |
| `mpirun -n 2 ./dpm.sc_mpi < data/input/photon_6MeV_smoke1e4.in` | PASS |

The electron and photon smoke runs both reported `No of histories simulated:
10000` and completed with the normal DPM termination message.

The same named inputs can now be reproduced with the packaged scripts:

```bash
MPI_RANKS=2 bash run/smoke_cpu_mpi.sh
MPI_RANKS=2 bash run/smoke_pezy_sc3s.sh
```

See `docs/USER_GUIDE.md` for build prerequisites, input and output definitions,
single-device commands, and basic result checks.
