# Himalaya

[![Build Status](https://travis-ci.org/Himalaya-Library/Himalaya.svg?branch=master)](https://travis-ci.org/Himalaya-Library/Himalaya)

Himalaya can calculate corrections of order O((alpha_t +
alpha_b)*alpha_s^2) to the CP-even Higgs mass matrix in the DR'-bar
scheme using the results of:

* R. V. Harlander, P. Kant, L. Mihaila and M. Steinhauser, *Higgs
  boson mass in supersymmetry to three loops*, [*Phys. Rev. Lett.*
  **100** (2008)
  191602](https://doi.org/10.1103/PhysRevLett.100.191602),
  [[0803.0672](https://arxiv.org/abs/0803.0672)],

* P. Kant, R. V. Harlander, L. Mihaila and M. Steinhauser, *Light MSSM
  Higgs boson mass to three-loop accuracy*, [*JHEP* **08** (2010)
  104](https://doi.org/10.1007/JHEP08(2010)104),
  [[1005.5709](https://arxiv.org/abs/1005.5709)].

Please refer to these papers as well as

* R. V. Harlander, J. Klappert and A. Voigt, *Higgs mass prediction in
  the MSSM at three-loop level in a pure DR context*, [*Eur.Phys.J.* **C77** (2017) no.12, 814]
  (https://doi.org/10.1140/epjc/s10052-017-5368-6),
  [[1708.05720](https://arxiv.org/abs/1708.05720)]

when using Himalaya.

## Requirements
The program requires:
* CMake >= 3.0
* Eigen3

## Installation
CMake is used to generate build files.
To create these build files, one should create a separate build directory inside Himalaya's top directory.
Run:
```
cd $HIMALAY_PATH
mkdir build
cd build
cmake ..
```

After calling `cmake` the build directory contains all required build files. Assuming that Makefiles are used, one can now run:
```
make
```
By default the code is compiled optimized.

## Running the code
After the compilation the static libraries `libDSZ.a` and `libHimalaya.a` have been created. The latter should be linked to the program. `libDSZ.a` is optional and has to be linked, if the program does not incorporate the associated Fortran code of G. Degrassi, P. Slavich and F. Zwirner ([arXiv:hep-ph/0105096](https://arxiv.org/abs/hep-ph/0105096)).

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
for the SPS2 benchmark point is given:

```cpp
himalaya::Parameters pars;                      // DR'-bar parameters struct
pars.scale = 1.11090135E+03;                    // renormalization scale
pars.mu = 3.73337018E+02;                       // mu parameter
pars.g3 = 1.06187116E+00;                       // gauge coupling g3 SU(3)
pars.vd = 2.51008404E+01;                       // VEV of down Higgs doublet
pars.vu = 2.41869332E+02;                       // VEV of up Higgs doublet
pars.mq2 << 2.36646981E+06, 0, 0,
            0, 2.36644973E+06, 0,
            0, 0, 1.63230152E+06;               // soft-breaking squared left-handed squark mass parameters
pars.md2 << 2.35612778E+06, 0, 0,
            0, 2.35610884E+06, 0,
            0, 0, 2.31917415E+06;               // soft-breaking squared right-handed down-squark mass parameters
pars.mu2 << 2.35685097E+06, 0, 0,
            0, 2.35682945E+06, 0,
            0, 0, 9.05923409E+05;               // soft-breaking squared right-handed up-squark mass parameters
pars.Ab = -784.3356416708631;                   // trilinear sbottom-Higgs coupling
pars.At = -527.8746242245387;                   // trilinear stop-Higgs coupling
pars.MA = 1.48446235E+03;                       // Mass of the CP-odd Higgs
pars.MG = 6.69045022E+02;                       // Mass of the Gluino
pars.MW = 8.04001915E+01;                       // Mass of the W boson
pars.MZ = 8.97608307E+01;                       // Mass of the Z boson 
pars.Mt = 1.47685846E+02;                       // Mass of the top quark
pars.Mb = 2.38918959E+00;                       // Mass of the bottom quark
pars.MSt << 9.57566721E+02, 1.28878643E+03;	// Masses of the stop quarks
pars.MSb << 1.27884964E+03, 1.52314587E+03;	// Masses of the sbottom quarks
pars.s2t = sin(2*asin(1.13197339E-01));         // 2 times the sine of the stop mixing angle
pars.s2b = sin(2*asin(-9.99883015E-01));        // 2 times the sine of the sbottom mixing angle
```

The input values of `MSt`, `MSb`, `s2t` and `s2b` are optional. If
they are not provided, they will get calculated internally.

Afterwards one can create a `HierarchyCalculator` object for the
chosen parameter set:

```cpp
himalaya::HierarchyCalculator hc(pars);
```

To calculate the DR'-bar loop corrections one needs to call:

```cpp
// the boolean argument switches between corrections proportional to alpha_t (false) or alpha_b (true)
himalaya::HierarchyObject ho = hc.calculateDMh3L(false);
```

All information which has been gathered during the calculation will be
stored in the returned `HierarchyObject` and can be accessed by member
functions.

To obtain the three-loop correction to the Higgs mass matrix one needs
to call:

```cpp
// returns a 2x2 matrix with the alpha_t*alpha_s^2 correction for the given parameter point
auto dMh3L = ho.getDMh(3);
```

The returned matrix should be added to the two-loop mass matrix
**before** diagonalization.

To obtain the three-loop correction to the quartic Higgs coupling λ of
the Standard Model in the DR'-bar scheme in the convention of
[[arXiv:1407.4081](https://arxiv.org/abs/1407.4081)] one needs to call

```cpp
double delta_lambda_3L = ho.getDeltaLambdaEFT();
```

The function `getDeltaLambdaEFT` returns the three-loop correction in
the DR'-bar scheme.  The three-loop shift to the MS-bar scheme,
ocurring when the one- and two-loop corrections are expressed in terms
of the Standard Model MS-bar strong gauge and top Yukawa couplings can
be obtained by calling `ho.getDRbarPrimeToMSbarShiftEFT()`.

An uncertainty estimate of the calculated three-loop λ can be obtained
by calling

```cpp
double delta_lambda_3L_uncertainty = ho.getDeltaLambdaUncertaintyEFT();
```

The function `getDeltaLambdaUncertainty` returns an uncertainty
estimate by taking into account logarithmic higher order Xt^n terms
missing in some hierarcy expansions.

A full and detailed example can be found in `source/example.cpp`.

### Mathematica interface

Since version 2.0 Himalaya can be run from within Mathematica using
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
  MstopMDRPrime -> {1807.421718155926, 2176.2140733830706}, 
  Mh2 -> {
    {{3.9900456677813977*^6, -199915.84939351628}, 
     {-199915.84939351628, 18267.112558603498}}, 
    {{-639.5967086048405, 38.1108098347363}, 
     {38.1108098347363, 10354.694221748936}}, 
    {{-2.067756466175066, 47.44906816768487}, 
     {47.44906816768487, 1872.6875317790061}}, 
    {{-4.186285565859146, 26.440305208276936}, 
     {26.440305208276936, 695.9601307370843}}
  }, 
  Mh2ShiftDRbarPrimeToMDRPrime -> {
    {{0., 0.}, {0., 0.}},
    {{0., 0.}, {0., 0.}},
    {{0., 0.}, {0., 0.}},
    {{-0.06556209414766911, 6.451828360749449},
     {6.451828360749449, -47.47567025934484}}
  },
  Mh2ShiftDRbarPrimeToH3m -> {
    {{0., 0.}, {0., 0.}},
    {{0., 0.}, {0., 0.}},
    {{0., 0.}, {0., 0.}},
    {{-1.5317742384518245, 1.9557280682757228}, 
     {1.9557280682757228, 7.446661229467915}}
  },
  expansionUncertainty -> {0., 0., 0.29936990221850124, 0.029834152694507742}, 
  lambda -> {0.13599819257964818, 0.06251358668613667, 0.0014909887921878413, 0.00021462203207113314}, 
  lambdaUncertainty -> {0., 0., 0., 0.0010928804743410628}, 
  lambdaShiftDRbarPrimeToMSbar -> {0., 0., 7.2898259577239126*^-6, -0.000771028340668861}, 
  deltaLambda3LoopH3m -> {0.000226477371058642, 0.0012078034614150654}, 
  deltaLambda3LoopH3mShiftDRbarPrimeToMSbar -> -0.0006679606925823672, 
  deltaLambda3LoopEFT -> {0.00021462203207113314, 0.0010928804743410628}, 
  deltaLambda3LoopEFTShiftDRbarPrimeToMSbar -> -0.000771028340668861 }
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
