# Himalaya

![](https://img.shields.io/github/v/release/Himalaya-Library/Himalaya)
[![Build Status](https://travis-ci.org/Himalaya-Library/Himalaya.svg?branch=master)](https://travis-ci.org/Himalaya-Library/Himalaya)
[![Build Status](https://github.com/Himalaya-Library/Himalaya/workflows/test/badge.svg)](https://github.com/Himalaya-Library/Himalaya/actions)

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

* R. V. Harlander, J. Klappert, and A. Voigt,
  *The light CP-even MSSM Higgs mass including N3LO+N3LL QCD corrections*,
  [*Eur.Phys.J.* **C80** (2020) no.3, 186](https://doi.org/10.1140/epjc/s10052-020-7747-7),
  [[1910.03595](https://arxiv.org/abs/1910.03595)]   

when using Himalaya.

## Requirements

Himalaya requires:
* C++ compiler
* FORTRAN compiler
* [CMake](https://cmake.org/) >= 3.1
* [Eigen3](http://eigen.tuxfamily.org)

## Building Himalaya

### Installation of dependencies

The required Eigen library can be installed using the package manager
of your Linux distribution.  On Debian/Ubuntu, for example, one may
run:

```sh
sudo apt-get install libeigen3-dev
```

Alternatively, the [Conan](https://conan.io/) package manager can be
used to install the dependencies:

```sh
mkdir -p build
cd build
conan install ..
```

### Compilation of Himalaya

Himalaya uses CMake to generate files for build automation.  To build
Himalaya one should first create a separate `build` directory inside
Himalaya's top directory.  Afterwards, `cmake` should be called:

```sh
mkdir -p build
cd build
cmake ..
```

`cmake` will search for required dependencies on your system, see the
previous section. If you want to use the dependencies installed by
Conan (see above), you need to run instead:

```sh
cmake .. -DCMAKE_TOOLCHAIN_FILE=./conan_paths.cmake
```

After calling `cmake` the build directory contains all required build
files. Assuming that GNU make is used, one can start the build by running

```sh
make
```

By default the example executable `example` is created from
[source/example.cpp](source/example.cpp), which prints all loop
corrections calculated by Himalaya for a given MSSM parameter point:

```sh
./example
```

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
#include "himalaya/HierarchyCalculator.hpp"
```

in the C++ file. The MSSM DR'-bar parameters must be stored in a
`Parameters` object.  Here, an example for the SPS1a benchmark point
is given:

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
pars.Mtau = 1.777;                   // Mass of the tau lepton
```

Afterwards one can create a `HierarchyCalculator` object for the
chosen parameter point:

```cpp
himalaya::HierarchyCalculator hc(pars);
```

To calculate the loop corrections in the DR'-bar scheme one needs to
call:

```cpp
// the boolean argument switches between corrections proportional to αt (false) or αb (true)
himalaya::HierarchyObject ho = hc.calculateDMh3L(false);
```

All information which has been gathered during the calculation will be
stored in the returned `HierarchyObject` and can be accessed by member
functions.

To extract the three-loop correction to the Higgs mass matrix one
needs to call:

```cpp
// returns a 2x2 matrix with the αt*αs^2 correction for the given parameter point
auto dMh3L = ho.getDMh(3);
```

To extract the three-loop correction to the quartic Higgs coupling λ of
the Standard Model in the DR'-bar scheme in the convention of
[[1407.4081](https://arxiv.org/abs/1407.4081)] one needs to call

```cpp
double delta_lambda_3L = ho.getDLambda(3);
```

The three-loop shift to the MS-bar scheme, ocurring when the one- and
two-loop corrections are expressed in terms of the Standard Model
MS-bar strong gauge and top Yukawa couplings can be obtained by
calling `ho.getDLambdaDRbarPrimeToMSbarShift(3)`.

An uncertainty estimate of the calculated three-loop λ due to the
truncation of the mass hierarchy expansions can be obtained by calling

```cpp
double delta_lambda_3L_uncertainty = ho.getDLambdaUncertainty(3);
```

A full and detailed example can be found in
[source/example.cpp](source/example.cpp).

**Example**:

```cpp
#include "himalaya/HierarchyCalculator.hpp"
#include <iostream>
#include <cmath>

himalaya::Parameters setup_point(double MS, double tb, double xt)
{
   himalaya::Parameters pars;

   const double MS2 = MS*MS;
   const double Xt = xt*MS;
   const double beta = std::atan(tb);
   pars.scale = MS;
   pars.mu = MS;
   pars.g1 = 0.46;
   pars.g2 = 0.65;
   pars.g3 = 1.166;
   pars.vd = 246*std::cos(beta);
   pars.vu = 246*std::sin(beta);
   pars.mq2 << MS2, 0, 0,
               0, MS2, 0,
               0, 0, MS2;
   pars.md2 << MS2, 0, 0,
               0, MS2, 0,
               0, 0, MS2;
   pars.mu2 << MS2, 0, 0,
               0, MS2, 0,
               0, 0, MS2;
   pars.ml2 << MS2, 0, 0,
               0, MS2, 0,
               0, 0, MS2;
   pars.me2 << MS2, 0, 0,
               0, MS2, 0,
               0, 0, MS2;
   pars.Au << 0, 0, 0,
              0, 0, 0,
              0, 0, Xt + pars.mu/tb;
   pars.Ad << 0, 0, 0, 0, 0, 0, 0, 0, 0;
   pars.Ae << 0, 0, 0, 0, 0, 0, 0, 0, 0;
   pars.Yu << 0, 0, 0, 0, 0, 0, 0, 0, 0.862;
   pars.Yd << 0, 0, 0, 0 ,0 ,0 ,0 ,0, 0.133;
   pars.Ye << 0, 0, 0, 0, 0, 0, 0, 0, 0.101;
   pars.MA = MS;
   pars.M1 = MS;
   pars.M2 = MS;
   pars.MG = MS;

   pars.validate(true);

   return pars;
}

int main()
{
   const auto point = setup_point(2000., 20., std::sqrt(6.));
   himalaya::HierarchyCalculator hc(point);

   try {
      // calculate the 3-loop corrections O(α_t*α_s^2)
      const auto ho = hc.calculateDMh3L(false);

      // extract 2x2 matrix with three-loop O(αt*αs^2) corrections
      const auto dMh3L = ho.getDMh(3);
      // extract three-loop O(αt*αs^2) correction to λ (DR'-bar scheme)
      const double delta_lambda_3L_DR = ho.getDLambda(3);
      // extract uncertainty estimate
      const double delta_lambda_3L_uncertainty = ho.getDLambdaUncertainty(3);
      // convert to MS-bar scheme
      const double delta_lambda_3L_MS =
         delta_lambda_3L_DR + ho.getDLambdaDRbarPrimeToMSbarShift(3);

      std::cout << "Δλ(3-loop,DR') = " << delta_lambda_3L_DR
                << " +- " << delta_lambda_3L_uncertainty << '\n'
                << "Δλ(3-loop,MS) = " << delta_lambda_3L_MS
                << " +- " << delta_lambda_3L_uncertainty << '\n';
   } catch (const std::exception& e) {
      std::cerr << e.what() << '\n';
   }

   return 0;
}
```

Output:

```
Himalaya info: Δλ(3-loop,DR') = 0.000315613 +- 0.00203118
Himalaya info: Δλ(3-loop,MS) = -0.000455415 +- 0.00203118
```

### Mathematica interface

Since version 2.0.0 Himalaya can be run from within Mathematica using
the LibraryLink interface.  To load Himalaya into a Mathematica
session, first, the file
[source/LibraryLink/Himalaya_LibraryLink.m](source/LibraryLink/Himalaya_LibraryLink.m)
must be loaded, which defines the Himalaya's Mathematica interface
functions.  Assuming the current directory is the `build/`
sub-directory of Himalaya, loading `Himalaya_LibraryLink.m` may be
done by calling

```.m
Get[FileNameJoin[{"..", "source", "LibraryLink", "Himalaya_LibraryLink.m"}]];
```

Afterwards, the LibraryLink `Himalaya_LibraryLink.so` must be loaded
using the `InitializeHimalaya[]` function:

```.m
InitializeHimalaya[FileNameJoin[{".", "Himalaya_LibraryLink.so"}]];
```

After the initialization, the function `HimalayaCalculateDMh3L[]` is
available, which calculates the loop corrections implemented in
Himalaya.  A full and detailed example can be found in
[source/example.m](source/example.m).

**Example**:

```.m
Get[FileNameJoin[{"..", "source", "LibraryLink", "Himalaya_LibraryLink.m"}]];

InitializeHimalaya[FileNameJoin[{".", "Himalaya_LibraryLink.so"}]];

MS = 2000;
TB = 20;
Xt = Sqrt[6] MS;

result = HimalayaCalculateDMh3L[
   settings -> {
       bottom -> False,
       loopOrder -> 3,
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

  Mh2EFTAt -> {8230.07, 10337.5, 830.688, 695.611},
  Mh2FO -> {8229.9, 8148.24, 819.187, 696.836},
  Mh2FOAt -> {0., 10331.8, 820.621, 696.836},
  lambda -> {0.135998, 0.0625136, 0.00149099, 0.000315613},
  lambdaUncertainty -> {0., 0., 0., 0.00203118},
  lambdaShiftDRbarPrimeToMSbar -> {0., 0., 7.28983 10^-6  , -0.000771028} }
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
