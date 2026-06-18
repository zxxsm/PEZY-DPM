# PEZY-DPM

PEZY-DPM is a PEZY-SC3s implementation of the Dose Planning Method (DPM) Monte Carlo dose calculation code for radiotherapy-class photon and electron transport. This repository accompanies the following manuscript, which is currently under revision:

> DPM Monte Carlo dose calculation on the PEZY-SC3s MIMD processor: porting and optimization

The code keeps the original DPM history-based transport structure and adds CPU MPI/OpenMP builds, PEZY-SC3s PZCL kernels, precision-controlled host/kernel types, MPI plus PEZY-SC3s execution, sparse dose scoring, dynamic scheduling, and 4-thread/8-thread PEZY kernel variants.

## License and use

PEZY-DPM additions are released under the GNU Affero General Public License v3.0 only (AGPL-3.0-only). See [`LICENSE`](LICENSE).

Commercial users who wish to use PEZY-DPM without complying with the AGPL-3.0 terms, for example in closed-source products, proprietary services, consulting workflows, or other deployments where the corresponding source code cannot be released under AGPL-3.0, should contact the copyright holders for a separate commercial license.

This software is a research code and is not approved for clinical treatment planning, medical diagnosis, or patient-specific clinical use.

Third-party components and data derived from the original DPM distribution retain their original notices. See [`NOTICE.md`](NOTICE.md) and the original DPM README preserved as [`docs/readme.txt`](docs/readme.txt).

## Repository layout

- `src/`: C/C++ and Fortran source code.
  - `Cfun.cpp`, `Cfun.hpp`, `Cglbvar.hpp`, `Cdglbvar.hpp`: C/C++ host code for CPU, MPI, OpenMP, and PEZY-SC3s execution.
  - `SC3s.cpp`: PEZY-SC3s host-side code using PZCL.
  - `dpm.F`, `time.f`, `getnam.f`, `libmath.f`, `libpenmath.f`, `libeloss.f`, `libphoton.f`, `libgeom.f`: retained Fortran DPM source files required by the current build.
- `pzc/`: PEZY kernel source. The package provides `kernel_4th.pzc` and `kernel_8th.pzc`.
- `kernel/`: PEZY kernel build wrapper.
- `compile/`: common Makefile fragments.
- `run/`: example run scripts.
- `data/input/`: runtime input files, including `command.in`, `electron_20MeV_2p5e8.in`, and `photon_6MeV_2p5e8.in`.
- `data/pre/`: preprocessed DPM physics tables, including electron and photon preprocessing tables such as `pre4elec.*` and `pre4phot.*`.
- `data/vox/`: voxel and geometry files.
- `data/seeds/`: pre-generated random seed files.
- `data/output/`: output directory used by example runs, logs, and simulation outputs.
- `data/benchmarks/`, `data/pendat/`: benchmark and material data inherited from the DPM workflow.
- `docs/readme.txt`: original DPM README and copyright notice.
- `TESTED_ON_100.md`: build and smoke-test report on server100.

## Build requirements

### CPU MPI/OpenMP

- GNU C++ compiler.
- GNU Fortran compiler.
- MPI compiler wrappers for MPI builds, such as `mpicxx` and `mpif90`.

### PEZY-SC3s

- PEZY-SC3s hardware.
- PEZY PZSDK/PZCL environment.
- MPI compiler wrappers.

Example environment on server100:

```bash
export PATH=/usr/local/mpich/bin:/opt/pezy/pzsdk/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/mpich/lib:/opt/pezy/pzsdk/lib:$LD_LIBRARY_PATH
export PZSDK_PATH=/opt/pezy/pzsdk
export PZCLKernelTimeout=1000000000
```

## Precision control

Use `PRECISION=float` or `PRECISION=double`.

- `PRECISION=float` is the default.
- `PRECISION=double` defines `DPM_USE_DOUBLE`; the C++ CPU code and PEZY PZC kernels use `double` for `UREAL`.
- The Fortran main program interface keeps the original DPM structure to avoid changing COMMON-block ABI assumptions.

## Build commands

### MPI

```bash
make clean MODE=mpi
make mpi PRECISION=float
make mpi PRECISION=double
```

The target binary is `dpm.mpi`.

### OpenMP

```bash
make clean MODE=openmp
make openmp PRECISION=float
make openmp PRECISION=double
```

The target binary is `dpm.omp`.

### MPI + PEZY-SC3s

```bash
make clean MODE=mpi_sc3s
make mpi_sc3s PRECISION=float KERNEL_VERSION=8th PZC_TARGET_ARCH=sc3s
make mpi_sc3s PRECISION=float KERNEL_VERSION=4th PZC_TARGET_ARCH=sc3s
make mpi_sc3s PRECISION=double KERNEL_VERSION=8th PZC_TARGET_ARCH=sc3s
make mpi_sc3s PRECISION=double KERNEL_VERSION=4th PZC_TARGET_ARCH=sc3s
```

The target binary is `dpm.sc_mpi`, and the PEZY binary is built as `kernel/kernel.pz`.

## Run examples

```bash
./dpm.omp < data/input/electron_20MeV_2p5e8.in > data/output/omp_electron.out
mpirun -n 64 ./dpm.mpi < data/input/photon_6MeV_2p5e8.in > data/output/mpi64_photon.out
mpirun -n 2 ./dpm.sc_mpi < data/input/photon_6MeV_2p5e8.in > data/output/sc3s_mpi2_photon.out
```

Additional examples are provided in `run/`.

## Validation status

The package was tested on server100 with GNU compilers, MPICH, and PEZY PZSDK. The following combinations passed compilation and short smoke runs:

- OpenMP: `float`, `double`
- MPI: `float`, `double`
- MPI + PEZY-SC3s: `float/double` with `KERNEL_VERSION=4th/8th`

See [`TESTED_ON_100.md`](TESTED_ON_100.md) for the detailed test report.

## Citation

If you use this software in academic work, please cite the associated manuscript when it becomes available and acknowledge the original DPM authors. A provisional citation metadata file is provided in [`CITATION.cff`](CITATION.cff) and should be updated after publication.
