# Himalaya

[![Build Status](https://travis-ci.org/Himalaya-Library/Himalaya.svg?branch=master)](https://travis-ci.org/Himalaya-Library/Himalaya)

Himalaya calculates three-loop corrections of order O((αt + αb)*αs^2)
to the MSSM CP-even Higgs mass matrix and to the quartic Higgs
coupling λ in the DR'-bar scheme using the results of:

* R. V. Harlander, P. Kant, L. Mihaila and M. Steinhauser, *Higgs
  boson mass in supersymmetry to three loops*, [*Phys. Rev. Lett.*
  **100** (2008) 191602](https://doi.org/10.1103/PhysRevLett.100.191602),
  [[0803.0672](https://arxiv.org/abs/0803.0672)],

* P. Kant, R. V. Harlander, L. Mihaila and M. Steinhauser, *Light MSSM
  Higgs boson mass to three-loop accuracy*, [*JHEP* **08** (2010)
  104](https://doi.org/10.1007/JHEP08(2010)104),
  [[1005.5709](https://arxiv.org/abs/1005.5709)].

Please refer to these papers as well as

* R. V. Harlander, J. Klappert and A. Voigt, *Higgs mass prediction in
  the MSSM at three-loop level in a pure DR context*, [*Eur.Phys.J.*
  **C77** (2017) no.12, 814](https://doi.org/10.1140/epjc/s10052-017-5368-6),
  [[1708.05720](https://arxiv.org/abs/1708.05720)],

* R. V. Harlander, J. Klappert, A. D. Ochoa Franco and A. Voigt,
  *The light CP-even MSSM Higgs mass resummed to fourth logarithmic order*,
  [*Eur.Phys.J.* **C78** (2018) no.10, 874](https://doi.org/10.1140/epjc/s10052-018-6351-6),
  [[1807.03509](https://arxiv.org/abs/1807.03509)]   

when using Himalaya.

## Requirements

Himalaya requires:
* C++ compiler
* FORTRAN compiler
* [CMake](https://cmake.org/) >= 3.1
* [Eigen3](http://eigen.tuxfamily.org)

## Building Himalaya

Himalaya uses CMake to generate files for build automation.  To build
Himalaya one should first create a separate `build` directory inside
Himalaya's top directory.  Afterwards, `cmake` should be called:

```
cd $HIMALAY_PATH
mkdir build
cd build
cmake ..
```

After calling `cmake` the build directory contains all required build
files. Assuming that GNU make is used, one can start the build by running

```
make
```

By default the code is compiled with optimizations.

## Running Himalaya

When the build is complete, the libraries `libDSZ` and `libHimalaya`
have been created.  The latter must be linked to user-written programs
to call the routines of Himalaya.  The library `libDSZ` is optional
and has to be linked in addition, if the program does not already
incorporate the associated FORTRAN code of P. Slavich for the 2-loop
corrections [[hep-ph/0105096](https://arxiv.org/abs/hep-ph/0105096),
[hep-ph/0112177](https://arxiv.org/abs/hep-ph/0112177),
[hep-ph/0212132](https://arxiv.org/abs/hep-ph/0212132),
[hep-ph/0305127](https://arxiv.org/abs/hep-ph/0305127)].

### C++ interface

We present a brief step by step guide how to run Himalaya at the C++
level and obtain the three-loop corrections to the CP-even Higgs mass
matrix in the MSSM in the DR'-bar scheme and the quartic Higgs
coupling of the Standard Model in the MS-bar scheme.

First one has to include the header

```cpp
#include "HierarchyCalculator.hpp"
```

in the C++ file. The DR'-bar parameters which define the MSSM
parameter must be stored in a `Parameters` object.  Here, an example
for the SPS1a benchmark point is given:

```cpp
himalaya::Parameters pars;           // DR'-bar parameters struct
pars.scale = 4.67491329E+02;         // renormalization scale
pars.mu = 3.52600579E+02;            // mu parameter
pars.g3 = 1.09949966E+00;            // gauge coupling g3 SU(3)
pars.vd = 2.49832484E+01;            // VEV of down Higgs doublet
pars.vu = 2.43650538E+02;            // VEV of up Higgs doublet
pars.mq2 << 2.99793928E+05, 0, 0,
            0, 2.99792102E+05, 0,
            0, 0, 2.49327504E+05;    // soft-breaking squared left-handed squark mass parameters
pars.md2 << 2.78275669E+05, 0, 0,
            0, 2.78273780E+05, 0,
            0, 0, 2.74928741E+05;    // soft-breaking squared right-handed down-squark mass parameters
pars.mu2 << 2.80477426E+05, 0, 0,
            0, 2.80475621E+05, 0,
            0, 0, 1.80478484E+05;    // soft-breaking squared right-handed up-squark mass parameters
pars.Ad << 0, 0, 0,
           0, 0, 0,
           0, 0, -784.3356416708631; // trilinear down type squark-Higgs coupling matrix
pars.Au << 0, 0, 0,
           0, 0, 0,
           0, 0, -527.8746242245387; // trilinear up type squark-Higgs coupling matrix
pars.MA = 3.92960954E+02;            // Mass of the CP-odd Higgs
pars.MG = 5.88220143E+02;            // Mass of the Gluino
pars.MW = 8.04136643E+01;            // Mass of the W boson
pars.MZ = 9.06817306E+01;            // Mass of the Z boson 
pars.Mt = 1.52117491E+02;            // Mass of the top quark
pars.Mb = 2.42010269E+00;            // Mass of the bottom quark
```

Afterwards one can create a `HierarchyCalculator` object for the
chosen parameter set:

```cpp
himalaya::HierarchyCalculator hc(pars);
```

To calculate the DR'-bar loop corrections one needs to call:

```cpp
// the boolean argument switches between corrections proportional to αt (false) or αb (true)
himalaya::HierarchyObject ho = hc.calculateDMh3L(false);
```

All information which has been gathered during the calculation will be
stored in the returned `HierarchyObject` and can be accessed by member
functions.

To obtain the three-loop correction to the Higgs mass matrix one needs
to call:

```cpp
// returns a 2x2 matrix with the αt*αs^2 correction for the given parameter point
auto dMh3L = ho.getDMh(3);
```

The returned matrix should be added to the two-loop mass matrix
**before** diagonalization.

To obtain the three-loop correction to the quartic Higgs coupling λ of
the Standard Model in the DR'-bar scheme in the convention of
[[1407.4081](https://arxiv.org/abs/1407.4081)] one needs to call

```cpp
double delta_lambda_3L = ho.getDLambda(3);
```

The function `getDLambda` returns the loop corrections in
the DR'-bar scheme.  The three-loop shift to the MS-bar scheme,
ocurring when the one- and two-loop corrections are expressed in terms
of the Standard Model MS-bar strong gauge and top Yukawa couplings can
be obtained by calling `ho.getDLambdaDRbarPrimeToMSbarShift(3)`.

An uncertainty estimate of the calculated three-loop λ can be obtained
by calling

```cpp
double delta_lambda_3L_uncertainty = ho.getDLambdaUncertainty(3);
```

The function `getDLambdaUncertainty` returns an uncertainty
estimate by taking into account higher order Xt^n terms missing in
some hierarcy expansions.

A full and detailed example can be found in `source/example.cpp`.

### Mathematica interface

Since version 2.0.0 Himalaya can be run from within Mathematica using
the LibraryLink interface.  To load Himalaya into Mathematica, first,
the file `source/LibraryLink/Himalaya_LibraryLink.m` must be loaded,
which defines the Himalaya's Mathematica functions.  Assuming the
current directory is in the `build/` sub-directory of Himalaya, this
can be done with:

```.m
Get[FileNameJoin[{"..", "source", "LibraryLink", "Himalaya_LibraryLink.m"}]];
```

Afterwards, the LibraryLink `Himalaya_LibraryLink.so` must be loaded
using the `InitializeHimalaya[]` function:

```.m
InitializeHimalaya["Himalaya_LibraryLink.so"];
```

Now, the `HimalayaCalculateDMh3L[]` function is available, which
calculates the loop corrections available by Himalaya.

**Example**:

```.m
Get[FileNameJoin[{"..", "source", "LibraryLink", "Himalaya_LibraryLink.m"}]];

InitializeHimalaya["Himalaya_LibraryLink.so"];

MS = 2000;
TB = 20;
Xt = Sqrt[6] MS;

result = HimalayaCalculateDMh3L[
   settings -> {
       bottom -> False,
       verbose -> True
   },
   parameters -> {
       scale -> MS,
       mu -> MS,
       g1 -> 0.46,
       g2 -> 0.65,
       g3 -> 1.166,
       vd -> 246*Cos[ArcTan[TB]],
       vu -> 246*Sin[ArcTan[TB]],
       mq2 -> MS^2 IdentityMatrix[3],
       md2 -> MS^2 IdentityMatrix[3],
       mu2 -> MS^2 IdentityMatrix[3],
       ml2 -> MS^2 IdentityMatrix[3],
       me2 -> MS^2 IdentityMatrix[3],
       Au -> {{0,0,0},
              {0,0,0},
              {0,0, Xt + MS/TB }},
       Ad -> 0 IdentityMatrix[3],
       Ae -> 0 IdentityMatrix[3],
       Yu -> {{0,0,0},
              {0,0,0},
              {0,0, 0.862 }},
       Yd -> {{0,0,0},
              {0,0,0},
              {0,0, 0.133 }},
       Ye -> {{0,0,0},
              {0,0,0},
              {0,0, 0.101 }},
       MA -> MS,
       M1 -> MS,
       M2 -> MS,
       M3 -> MS
   }
]
```

The `result` contains a list with replacement rules for all loop
corrections:

```.m
{ hierarchyID -> 1, hierarchyName -> h32q2g, 
  MstopMDRPrime -> {1807.42, 2176.21}, 
  Mh2 -> {
    {{3.99005*10^6, -199916.}, {-199916., 18267.1}},
    {{-639.597, 38.1108}, {38.1108, 10354.7}},
    {{-2.06776, 47.4491}, {47.4491, 1872.69}},
    {{-4.18629, 26.4403}, {26.4403, 695.96}}
  },
  Mh2ShiftDRbarPrimeToMDRPrime -> {
    {{0., 0.}, {0., 0.}},
    {{0., 0.}, {0., 0.}},
    {{0., 0.}, {0., 0.}},
    {{-0.0655621, 6.45183}, {6.45183, -47.4757}}
  },
  Mh2ShiftDRbarPrimeToH3m -> {
    {{0., 0.}, {0., 0.}},
    {{0., 0.}, {0.,  0.}},
    {{0., 0.}, {0., 0.}},
    {{-1.53177, 1.95573}, {1.95573, 7.44666}}
  },
  expansionUncertainty -> {0., 0., 0.29937, 0.0298342},
  Mh2EFT -> {8230.07, 10331.1, 1881.15, 695.229},
  lambda -> {0.135998, 0.0625136, 0.00149099, 0.000315613},
  lambdaUncertainty -> {0., 0., 0., 0.00100837},
  lambdaShiftDRbarPrimeToMSbar -> {0., 0., 7.28983*10^-6, -0.000771028} }
```

See `?HimalayaCalculateDMh3L` for a detailed documentation of the
input and output.

## Code Documentation

Doxygen can be used to generate code documentation.
To generate the documentation, run

```
make doc
```

The generated documentation can be found in `doc/html/index.html`.
