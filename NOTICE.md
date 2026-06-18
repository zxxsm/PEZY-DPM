# Notices and third-party components

This package contains PEZY-DPM additions developed for the manuscript "DPM Monte Carlo dose calculation on the PEZY-SC3s MIMD processor: porting and optimization", together with third-party software components and data files required by the DPM workflow.

## Original DPM

The PEZY-DPM implementation is derived from the Dose Planning Method (DPM) Monte Carlo dose calculation code.

The original DPM README included in this repository as `docs/readme.txt` states:

- Copyright (c) 2000,2001 Polytechnical University of Catalonia.
- Copyright (c) 2000,2001 The University of Michigan.
- Copyright (c) 2000,2001 The University of Barcelona.

The README grants permission to use, copy, modify, distribute, and sell the software and documentation for any purpose without fee, provided that the copyright notice and permission notice appear in copies and supporting documentation. The software is provided "as is" without express or implied warranty.

Users should consult `docs/readme.txt` for the exact original notice.

## Random seed files

Pre-generated random seed files are included in `data/seeds/`. These files are provided as benchmark/runtime data for the DPM workflow. The seed-generation source utility is not included in this repository. If regenerated seed sets are needed, use a compatible RANECU sequence-splitting generator and document the modulus, multiplier, number of generated seeds, and substream spacing.

## PEZY-DPM additions

The PEZY-DPM additions include:

- PEZY-SC3s/PZCL host-device execution support.
- MPI plus PEZY-SC3s execution support.
- Mixed-precision transport changes.
- PEZY-SC3s kernel-side data layout, sparse dose-scoring buffers, atomic dynamic scheduling, and active-thread mapping.
- Build and run scripts added for the manuscript experiments.

These additions are released under the GNU Affero General Public License v3.0 only (AGPL-3.0-only), as provided in `LICENSE`.

Commercial users who wish to use PEZY-DPM without complying with the AGPL-3.0 terms, including users who need closed-source or proprietary redistribution, service deployment, consulting workflow integration, or product integration, should contact the copyright holders for a separate commercial license.

This software is a research code and is not approved for clinical treatment planning, medical diagnosis, or patient-specific clinical use.
