<---- README for DPM programs  -------------------------------------->
      This file was created on May 2, 2001.

CONTENTS

 1) Copyright
 2) Acknowledgement
 3) General
 4) How to install DPM
 5) How to run DPM
 6) Programs
 7) Data files
 8) Parameters and user-defined routines
 9) Sample problem
10) Benchmarks and timings
11) Changes in v1.1



<---- 1. COPYRIGHT  ------------------------------------------------->

Copyright (c) 2000,2001
Polytechnical University of Catalonia

Permission to use, copy, modify, distribute and sell this software and
its documentation for any purpose is hereby granted without fee,
provided that the above copyright notice appear in all copies and that
both that copyright notice and this permission notice appear in
supporting documentation.  The Polytechnical University of Catalonia
makes no representations about the suitability of this software for
any purpose. It is provided "as is" without express or implied
warranty.

Copyright (c) 2000,2001
The University of Michigan

Permission to use, copy, modify, distribute and sell this software and
its documentation for any purpose is hereby granted without fee,
provided that the above copyright notice appear in all copies and that
both that copyright notice and this permission notice appear in
supporting documentation.  The University of Michigan makes no
representations about the suitability of this software for any
purpose. It is provided "as is" without express or implied warranty.

Copyright (c) 2000,2001
The University of Barcelona

Permission to use, copy, modify, distribute and sell this software and
its documentation for any purpose is hereby granted without fee,
provided that the above copyright notice appear in all copies and that
both that copyright notice and this permission notice appear in
supporting documentation.  The University of Barcelona makes no
representations about the suitability of this software for any
purpose. It is provided "as is" without express or implied warranty.



<---- 2. ACKNOWLEDGEMENT  ------------------------------------------->

We (Scott Wilderman, Josep Sempau and Alex Bielajew) would like to
thank Drs. Dick Fraass, Kwok Lam, Dan McShan and Randy Ten Haken, from
the Radiation Oncology department of the University of Michigan, Dr Ed
Larsen from the Department of Nuclear Engineering and Radiological
Sciences, and Dr. F. Salvat, from the ECM Department of the
Universitat de Barcelona, for fruitful discussions. One of the authors
(JS) gratefully acknowledges the financial support of the Direccio
General de Recerca de la Generalitat de Catalunya (Spain) grant no.
1997BEAI400256 as well as the financial support from the Spanish Fondo
de Investigacion Sanitaria under contract no. 98/0047-01. Financial
support for this work has been provided by ADAC Laboratories
(Milpitas, California). In addition, financial support for this work
has been provided under the auspices of the U.S. Department of Energy
by the Lawrence Livermore National Laboratory under contract number
W-7405-ENG-48.



<---- 3. GENERAL ---------------------------------------------------->

DPM is a computer code that simulates the transport of high energy
electrons and photons using Monte Carlo (MC) methods.

DPM has been designed to deal with radiotherapy-class problems. More
specifically, it is optimized for the calculation of the dose
distribution delivered by high energy (~1 MeV up to ~20 MeV) electron
and photon beams to patients, whose geometries are defined in terms of
a (usually) large number of small volume elements, or voxels. Some
approximations in the description of the physical processes involved
assume that voxels are made of low atomic number materials, such as
water, lung, bone, titanium, etc.

There are 4 programs in the DPM suite. The first one is a slightly
modified version of the program MATERIAL that comes with the general
purpose MC code PENELOPE. MATERIAL extracts raw physics data from
PENELOPE's library of constants, and produces tables and cross
sections for given materials. The second program is a pre-processing
program, PREDPM, that uses PENELOPE data (generated with MATERIAL) and
produces a series of input files for a DPM run. The third program,
GENVOXEL, takes a PENELOPE-style surface/solid geometry file and
creates a voxel-by-voxel list of volume densities for DPM. The fourth
program is the Monte Carlo program itself, DPM.

It is assumed that the user is fairly familiar with Monte Carlo
electron-gamma transport programs. Concepts such as cutoff energies
for secondary particle producing processes and absorption energies are
referenced but not explained.

At some points it is also assumed that the user is running on a Unix
environment. This assumption affects things such as the names of the
executables (e.g. 'dpm.x' as opposite to 'dpm.exe') and the commands
in the make-all file (e.g. 'mv' instead of 'move').

Except when otherwise noted, DPM handles energies in eV and distances
in cm.



<---- 4. HOW TO INSTALL DPM ----------------------------------------->

After downloading and unzipping you should have the following files
under your DPM directory,

readme.txt   ->  This file
*.f          ->  FORTRAN 77 sources
make-all     ->  Combined make file (compile+link) for the 4 programs
*.in         ->  Input files for the sample problem
sample.geo   ->  PENELOPE-style geometry file for the sample problem
sample.out   ->  Results of the sample problem
_pendat/*    ->  Material data
benchmarks/* ->  ICCR2000 benchmark figures and photon spectrum

Now follow these steps:

1) Have a look at 'make-all'. It is provided as a way to show what
   FORTRAN files need to be linked together. 'make-all' is prepared
   for a Linux environment using the GNU g77 compiler with no
   optimization options. Adapt it to your system, invoking your
   FORTRAN 77 compiler and using the best optimization options
   available.

2) DPM has been prepared to make use of some time functions in order
   to provide information regarding simulation efficiency after each
   run. However, depending on your operative system and compiler,
   these system functions may not work on your computer or be invoked
   in a way different from the one used by DPM. Therefore, to avoid
   complications in the compilation and link process, all calls to
   time functions are NOT activated by default.

   If you have a GNU compiler and linker you may be interested in
   turning time functions ON at this point. To do so, edit the file
   'time.f' and follow the instructions provided next to each
   occurrence of the string 'opt-Time ON'. If you don't have the GNU
   compiler but you know how to invoke the corresponding functions in
   your system, then you can modify time.f so that each routine
   provides the information described in the corresponding header.

3) Execute

   % make-all

   from the DPM directory. Notice that on an Unix system you may need
   first to change the permissions with

   % chmod u+x make-all

   This step will create four executable files, namely dpm.x,
   predpm.x, genvoxel.x and _pendat/material.x.

You are now ready to go. Proceed to the next section.



<---- 5. HOW TO RUN DPM --------------------------------------------->

Once DPM has been 'installed', follow these steps:

1) Create your DPM geometry file (here named *.vox). In section "Data
   files" the format of a .vox file is presented.

   When the geometry to be described is made of one or several
   homogeneous bodies limited by quadric surfaces (planes, spheres,
   cylinders, etc), the auxiliary program GENVOXEL can be employed. In
   this case, a PENELOPE-style geometry file (*.geo) must be prepared
   beforehand. GENVOXEL will read the .geo file and create the
   corresponding .vox file compatible with DPM for you. GENVOXEL also
   needs an input file (e.g. genvoxel.in). To run,

   % genvoxel.x < genvoxel.in

2) Create your input file for PREDPM (here named predpm.in)

   Two examples of input files, predpm4elec.in and predpm4phot.in, are
   provided. In these two files the values assigned to the simulation
   parameters shigh, slow and ecross are appropriate for the
   simulation of electron and photon therapeutic beams respectively.

   Make sure you have one DPM material file for each material referred
   to in predpm.in. DPM material files are stored in _pendat/. Some of
   the most frequently used materials have already been included in
   the distribution (e.g. Water). Note that the input files provided
   as an example are prepared to generate five materials, namely,
   water, lung (ICRP definition), cortical bone (ICRP), aluminum and
   titanium. If, however, a material is needed and not present, it can
   be easily prepared by running the interactive program MATERIAL,
   from the PENELOPE suite. Notice that this latter program must be
   located in the _pendat/ directory ('make-all' automatically places
   it there).

3) Run predpm.x < predpm.in

   This will create a set of files (with the common prefix specified
   in predpm.in) to be read by dpm. This step needs to be done only
   once for a given set of materials.

4) Create your input file for DPM.

   Material and geometry files are specified by name in the input
   file. An example is provided for the sample problem described below
   (dpm.in). Recall that all the information regarding timings (such
   as "allotted time" in dpm.in) will have no effect if you did not
   activate the routines in time.f (see section "How to install").

5) Run dpm.x < dpm.in > dpm.out

   The result of the simulation will be written to dpm.out after
   completion.

See section "Sample problem" for more details.



<---- 6. PROGRAMS --------------------------------------------------->

The following FORTRAN 77 files are included,

dpm.f                - DPM source
genvoxel.f           - GENVOXEL source
getnam.f             - program for getting file name prefixes
inscat.f             - routines for sampling Compton angles exactly,
                       taken from PENELOPE
libelastic.f         - library of elastic and multiple elastic scattering
                       programs
libeloss.f           - library of energy loss and inelastic scattering routines
libgeom.f            - library of simple voxelized geometry programs
libmat.f             - library of routines for working with compounds, etc.
libmath.f            - library of useful math routines
libpenmath.f         - library of useful math routines adapted from PENELOPE
libphoton.f          - library of photon routines, some used by DPM, some by
                       PREDPM, and some for reference
material.f           - PENELOPE2000 material processing program
pbrem.f              - routines to set up accurate bremsstrahlung sampling in
                       PREDPM
penelope.f           - PENELOPE2000 Monte Carlo and pre-processing library
pengeom2.f           - PENELOPE2000 geometry routines
predpm.f             - PREDPM source
time.f               - Time routines

Note that compiling some of the subprograms which borrow heavily from
PENELOPE with explicit type checking invoked may generate WARNING
messages.


1) DPM

DPM is optimized for radiotherapy-class problems as described in
"DPM, a fast, accurate Monte Carlo code optimized for photon and
electron radiotherapy treatment planning dose calculations", Physics
in Medicine and Biology 45 (2000) p 2263-91. DPM 1.1 introduces some
slight modifications in the physics described in this paper.

It reads its input data from stdin and writes to stdout, both of which
are normally redirected to a file, like in

% dpm.x < dpm.in > dpm.out

The file dpm.in provided for the sample problem should be self-
explanatory.

In addition to the input file, DPM also reads a file called
'command.in'. This file contains instructions to be performed by DPM
at execution time, typically to query about the results obtained up to
that moment (which are written to 'progress.out') or to change the
settings regarding the number of histories or the simulation time
requested. In particular, it can be used to stop the simulation at
once. The command file is a self-explanatory plain text file that can
be edited by the user while DPM is running.


2) GENVOXEL

GENVOXEL reads a PENELOPE formatted geometry file and converts it into
a plain text 'voxel' file in the format required by DPM.

It reads its input data from stdin, which is normally redirected to an
input file, like in

genvoxel.x < genvoxel.in

The file genvoxel.in provided for the sample problem should be self-
explanatory.

It is necessary to have a PENELOPE geometry file prepared in advance.
Refer to PENELOPE documentation for details. An example for the sample
problem is given in sample.geo. PENELOPE is freely available from NEA.
See http://www.nea.fr


3) MATERIAL

MATERIAL is an interactive program that constructs PENELOPE data files
for given materials. To run,

% material.x

Input material is queried from the command line, and can either be
specified in terms of Z and weight fractions, or can be specified from
a default list of PENELOPE materials. The list is found in
_pendat/compdata.tab.

The output file is loaded into _pendat, with a name specified by the
user. Note that this file name must match the material names in the
PREDPM input file (and recall that file names are case sensitive in
Unix!). During the course of execution of MATERIAL, the user will be
asked to input "ENERGY AND OSCILLATOR STRENGTH OF THE PLASMON".  It is
recommended that the default values given in the program output be
used for all non-insulating materials. For insulators, these
quantities should be set equal to 0.


4) PREDPM

PREDPM prepares material data to be read by DPM.

Some local modifications to PREDPM may be required.  In particular
there may be some machine dependencies in the subroutines GETNAM and
GETNAM2, which are found in files inscat.f and getnam.f.  These
routines are used to read in variable length names of materials and
DPM input file prefixes.

PREDPM reads data files that are generated by the PENELOPE pre-
processor MATERIAL (see above).

To run,

% predpm.x < predpm.in

Problem parameters (cutoffs, materials, etc) are specified in the
input file, as described below. A series of output files (DPM data
inputs) are generated, all having the same prefix, as specified in the
input file.




<---- 7. DATA FILES ------------------------------------------------->

Generated by PREDPM
-------------------

PREDPM generates a set of data files to be read by DPM. See section
"Programs" for more details.


The DPM Geometry file
---------------------

Its structure is as follows:

FILE BEGINS HERE *********************************************
 HEADER section: DPM geometry file     --> This is a comment line
                                       --> This blank line MUST be here
 Nxvox,Nyvox,Nzvox:
     61     61    150                  --> No of voxels in {x,y,z} directions
 dxvox,dyvox,dzvox (cm):
  5.000E-01  5.000E-01  2.000E-01      --> Voxel size {dx,dy,dz} (cm)
 Mat#[nx,ny+,nz++] -- dens(g/cm^3):
 1 1.000                               --> 1st column is material no. (following
 1 1.000                                   the order in which they are introduced
 1 1.000                                   in the PREDPM input file). 2nd column is
 1 1.000                                   density.
 1 1.000
 1 1.000                                   Note that the 1st
 1 1.000                                   row refers to the voxel with one corner
 1 1.000                                   in the (0,0,0) position and the center
 1 1.000                                   at (dx/2,dy/2,dz/2); the 2nd row refers
 1 1.000                                   to the voxel centered at (dx/2,dy/2,3dz/2)
 4 2.699                                   and so on. In the example shown, row #151
 4 2.699                                   (there are 150 voxels in the z-direction)
 4 2.699                                   would refer to the voxel centered at
 4 2.699                                   (dx/2,3dy/2,dz/2).
 4 2.699
 2 2.600E-01
 2 2.600E-01
 2 2.600E-01
 2 2.600E-01
 2 2.600E-01
  ...etc

END OF FILE *************************************************


In _pendat/.
------------

penepN.tab     -> electron and positron elemental raw data files for
                  MATERIAL
penphN.tab     -> photon elemental raw data files for MATERIAL
compdata.tab   -> composition and physical parameters for several hundred
                  materials
idlist.tab     -> list of the PENELOPE ID numbers of the materials in
                  compdata.tab
A150           -> MATERIAL output files for common media, provided for
AdiposeTissue     convenience
etc.



<---- 8. PARAMETERS AND USER-DEFINED ROUTINES ----------------------->

PARAMETERS

The operation of DPM is controlled by the parameters shown in the
table below.

Name    File       Value recommended         Comment
------------------------------------------------------------------
Emin    predpm.in  199.0e3 eV     Minimum e- energy in look-up tables

Emax    predpm.in  21.0e6  eV     Maximum e- energy in look-up tables

Eminph  predpm.in  49.0e3  eV     Minimum photon energy in look-up tables

Wcc     predpm.in  200.0e3 eV     Cutoff energy for delta ray production

Wcb     predpm.in  50.0e3  eV     Cutoff energy for bremsstrahlung production

Eabs    dpm.in     = wcc          Absorption energy for e-

Eabsph  dpm.in     = wcb          Absorption energy for photons

shigh   predpm.in  5.0 cm for e-  Step length for energies above 'ecross'
                   1.0 cm for photons

slow    predpm.in  1.0 cm for e-  Step length for energies below 'ecross'
                   0.5 cm for photons

ecross  predpm.in  12.0e6 eV e-   See 'shigh' and 'slow'
                   5.0e6 eV photons


Notice that shigh, slow and ecross should be given different values
depending on whether one wants to simulate electron or photon beams.
Thus, PREDPM must be run twice for each set of materials, changing the
prefix of the generated data files to avoid confusion. Two different
PREDPM input files, named predpm4elec.in (used for the sample problem)
and predpm4phot.in have been provided with the recommended values.

Note also that all the names of the materials in the PREDPM input file
have to match exactly file names created by MATERIAL in _pendat (i.e.
Water, LungICRU, CorticalBone, etc).


SOURCE MODELS

DPM sets the initial state of primary particles in subroutine SOURCE.
As it stands, SOURCE generates monodirectional, monoenergetic
electrons (or photons, depending on the corresponding input in dpm.in)
that impinge on the phantom surface through a square field centered on
the x=y=0 face of the geometry. The side of the field is also defined
in dpm.in.

However, the user may want to define different types of initial
states, e.g. to simulate primary photons with an energy distribution
that reproduces the spectrum produced by an accelerator head. It is
virtually impossible to write a SOURCE subroutine that is able to cope
with all the possibilities that one may think of. Instead, DPM
provides full flexibility by allowing users to write their own SOURCE.
To help in this task, some (commented out) examples of how to code
other initial states are also provided in dpm.f.

In any case, SOURCE must always set the following variables, all
contained in common /dpmpart/,

energy         -> kinetic energy of the particle (eV)
{vx,vy,vz}     -> direction of flight, normalized to 1
{x,y,z}        -> position coordinates (cm)
ptype          -> 0 for photons, -1 for e- and +1 for e+

After setting {x,y,z}, subroutine WHERE must be called. Notice that
the simulated geometry is assumed to lie in the region {x,y,z} >= 0.


SCORING

Some users may find useful to modify subroutine SCORE in order to keep
track of quantities of interest other than the dose (e.g. particle
fluence) or to define different regions of interest (RoIs). Notice
that, in the default version, SCORE scores only energy depositions
occurring inside a user-defined RoI.


REPORT

Similarly to SOURCE, subroutine REPORT can be modified to suit the
users needs. In the present version, REPORT prints out the dose
delivered to the voxels inside some region of interest defined by the
user in the DPM input file.



<---- 9. SAMPLE PROBLEM --------------------------------------------->

Follow these steps to run the sample problem,

1) Follow the steps described in section "How to install DPM".

2) From the DPM directory:

% cd _pendat
% material.x  (prepare water --ID # 277 in PENELOPE's list)
% material.x  (prepare aluminum --ID # 13 in PENELOPE's list)
% cd ..
% genvoxel.x < genvoxel.in
% predpm.x < predpm4elec.in
% dpm.x < dpm.in > dpm.out


MATERIAL will create files for water and aluminum which should be
identical, within machine precision, to 'Water' and 'Al'. When
prompted, provide names for output files diffent from 'Water' and
'Al', so that you do not overwrite the files included in the
distribution. This step is done only to exercise the use of MATERIAL.
Once they have been proven identical, you can delete your water and
aluminum files and use the originals instead. Please, notice that file
names 'Water' and 'Al' are case sensitive. These file names must match
the ones in 'predpm4elec.in'. During the course of execution of
MATERIAL, the user will be asked to input "ENERGY AND OSCILLATOR
STRENGTH OF THE PLASMON".  It is recommended that the default values
given in the program output be used for all non-insulating materials.
For insulators, these quantities should be set to 0.

GENVOXEL will create a file 'sample.vox'.

PREDPM will create a series of input files prefixed 'pre4elec'.

DPM should give results (in dpm.out) similar to those found in
sample.out, within statistical uncertainties. Note that if a round-off
difference between the user's machine and the DPM benchmark platform
leads to a different branching in a single IF test in the sample
problem, the subsequent random walks simulated on the two machines
could be completely different.  Thus, as with all Monte Carlo programs
being run on different platforms, exactly matching results should not
be expected.

The report printed corresponds to the central axis depth dose produced
by a 20 MeV parallel e- beam (1.5 x 1.5 cm^2) on a water/Al/water
(2/3/30 cm) phantom. The simulation of the sample problem should take
less than 7 minutes on an Intel Pentium III/733 MHz computer.


<---- 10. BENCHMARKS AND TIMINGS ------------------------------------>

The results of some benchmarks performed with DPM 1.1 can be found
in benchmarks/.

The following files are included:

iccr-e.ps      -> ICCR benchmark for a 20 MeV e- beam
iccr-ph.ps     -> ICCR benchmark for 18 MeV photon beam
water10m.ps    -> homogeneous water, 10 MeV e- beam
water05m.ps    -> homogeneous water,  5 MeV e- beam
ti20m.ps       -> homogeneous titanium, 20 MeV e- beam
ti-ph.ps       -> homogeneous titanium, 18 MeV photon beam

A detailed description of the ICCR2000 benchmark can be found in
http://www.irs.inms.nrc.ca/inms/irs/papers/iccr00/iccr00.html (as of
the time of writing this file). The spectrum used for the simulation
of the photon beam cases can be obtained from the same web page. It is
also available in the file benchmarks/ICCR2000-18Mphot.spc

All benchmarks refer to the central axis depth dose in a 30.5 x 30.5 x
30 cm phantom filled with 5 x 5 x 2 mm^3 voxels. The source position
and field size (and energy spectrum in the photon case --- e- are
monoenergetic) are the same as in the ICCR benchmark.

In all figures,
- Dashed lines represent depth doses calculated with PENELOPE
- Solid green lines are dose differences (DPM-PENELOPE) relative
  to the PENELOPE dose maximum (in percentage).
- Solid black lines (when shown) are dose differences (EGS4-PENELOPE)
  relative to the PENELOPE dose maximum (in percentage).


The table below shows the timings (in seconds) for each case.

CASE        HISTORIES(*10^6)   SIGMA(%)   TIME(s)
--------    ----------------   --------   -------
iccr-e        1.5               0.249      264

iccr-ph       20                0.387      902

water10m      1.5               0.264      100

water05m      1.5               0.260       51

ti20m         1.5               0.260      422

ti-ph         20                0.394     3260


SIGMA refers to the relative standard deviation averaged over all
voxels (on the central axis) having a dose larger than half the dose
maximum.

These timings were obtained with a PC (AMD K7 processor at 1200 MHz
and 128 MB of RAM) and using the optimization options
g77 *.f -O3 -ffast-math -malign-double -mcpu=i686 -march=i686

It is important to emphasize that DPM does not take advantage of
neighbor voxels being made of the same material. Therefore, very
similar timings can be expected for highly heterogeneous situations,
such as those encountered in some cases in the clinical practice.



<---- 11. CHANGES IN DPM 1.1  --------------------------------------->

Some of the changes introduced in the current version affect the
interface with the user, whereas some other do not. Among the ones
that do affect the interface, the most relevant are

1) The format of the geometry file has been slightly changed. To
   convert an old .vox file to the new format, simply delete the
   following two lines from the header:

        Using predpm set no:
        0

  so that the first 3 lines of the file look like:

        HEADER section: comments go here
                <blank line>
        Nxvox,Nyvox,Nzvox:

2) In the new DPM, two different step lengths (shigh & slow) and one
   crossover energy (Ecross) are defined (in the input file). 'shigh'
   operates at energies above Ecross, whereas slow is active below.

   Due to the different approach adopted for determining distances
   between elastic events, the same step length parameter will not
   produce the same actual steps in the new and old DPM versions,
   being somewhat shorter in DPM 1.1.

3) All data files generated by PREDPM have to be regenerated with the
   new version. The old format is not compatible with DPM 1.1.

4) DPM calls routines contained in time.f which are compatible with
   the GNU FORTRAN 77 compiler, but may not be supported by your
   system. To deactivate these calls, just edit time.f search for the
   string "opt-Time OFF" and follow the instructions provided.

5) Subroutine REPORT prints out the dose delivered to a region of
   interest defined by the user in the DPM input file.

Other changes:

6) The PENELOPE bremsstrahlung DCS has been activated by default.
   Thus, the routines lambre() and stpbre() have been substituted by
   plambr() and pstpbr() respectively.

7) A more accurate method to calculate the rate at which the
   bremsstrahlung fuel is burnt (instead of assuming that it is
   constant with E) has been introduced by default. This implies
   frequent calls to routine ilabip(), which interpolates the inverse
   bremsstrahlung MFP as a function of material and energy, thus
   incrementing the time spent to simulate each history. However, the
   improvement observed in some extreme cases (i.e. 20 MeV e- on Ti)
   provides the motivation.

8) The accuracy with which elastic fuel is burnt has been improved (in
   routine flight()). With the new procedure, the distance to the next
   scattering event is much less dependent on factors such as voxel
   size. The cost of this improvement is more calls to routine
   scpwip(), which interpolates the inverse 1st transport MFP
   (scattering power).

9) To somewhat compensate for the former changes, the following
   routines have been optimized:

   rstpip (stopping power)
   ilabip (bremsstrahlung mean free path)
   scpwip (scattering power)

   Now linear interpolation is used instead of 3splines.

10) The routines electr() and flight() have been optimized.

11) Small changes (mainly aesthetic or to improve speed) to other
   routines in dpm.f and its libraries (such as rotate) have been
   introduced. libelastic.f and libmat.f have been reviewed.

12) The way in which the spreading parameter and the q surface are
   obtained by PREDPM has been changed. Goudsmit&Saunderson theory is
   now used instead of Lewis'.

   The former produces a much faster predpm and gives practically the
   same dose distributions because in both cases the relevant quantity
   is the accumulated scattering strength (defined as the integral of
   the inverse 1st transport MFP (lambda_1^{-1}) over the step
   length), which determines the first term in the respective series.

13) An option that allows using a different q surface for each
   material has been introduced (NOT active by default). This implies
   that predpm.f calculates now the spreading parameter and the q
   surface for every material in the input file. To activate this
   option, search for the string "opt-qSingleMat" and follow the
   instructions provided.


BUGS fixed.

14) A bug in rvoxg() has been corrected. The density was set to the
   value assigned to the reference material, irrespective of the
   values read from the .vox geometry file.

15) A factor Z(Z+1) was incorrectly plugged into the bremsstrahlung
   routines. It should be Z^2.


**** END OF FILE *****************************************************


