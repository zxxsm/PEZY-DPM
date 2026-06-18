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
