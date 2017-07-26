# Himalaya

Himalaya can calculate corrections of the order O((alpha_b + alpha_t)*alpha_s^2 ) to the CP-even Higgs mass matrix in the DR-bar scheme using the results of Kant et al. (arXiv:0803.0672, arXiv:1005.5709).

## Requirements
The program requires:
* CMake >= 3.0
* Eigen3

## Installation
CMake is used to generate build files.
To create these build files, you have to create a separate build directory.
Run:
```
mkdir build
cd build
cmake PATH_TO_HIMALAYA_DIR
```
By default the code is compiled optimized.

After calling `cmake` the build directory contains all required build files. Assuming that Makefiles are used, you can now run:
```
make
```

## Running the code
After the compilation the static libraries `libDSZ.a` and `libHimalaya.a` have been created. The latter should be linked to your program. `libDSZ.a` is optional and has to be linked, if your program does not incorporate the associated Fortran code of arXiv:hep-ph/0105096.

### Example
We present a brief step by step guide how to run Himalaya and obtain the three-loop results. First you have to include the headers
```cpp
#include "HierarchyCalculator.hpp"
#include "Himalaya_interface.hpp"
#include "HierarchyObject.hpp"
```
to your file. In the next step you have to initialize all needed parameters in the `Parameters` `struct`. Note that the input has to be provided in the **DR-bar** scheme. Here, an example for a SPS2 benchmark point is given:
```cpp
himalaya::Parameters pars;                      // DR-bar parameters struct
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
pars.Ab = -784.3356416708631;                   // soft-breaking trilinear term Ab
pars.At = -527.8746242245387;                   // soft-breaking trilinear term At
pars.MA = 1.48446235E+03;                       // Mass of the A boson
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
The input values of `MSt`, `MSb`, `s2t` and `s2b` are optional. If they are not provided, they will get calculated internally.

Now you can create a `HierarchyCalculator` object with the given `struct`:
```
himalaya::HierarchyCalculator hc(pars);
```
To calculate the results you just have to call:
```cpp
himalaya::HierarchyObject ho = hc.calculateDMh3L(false, 0); // the bool argument switches between corrections proportional to alpha_t (false) or alpha_b (true). The integer is a flag to choose among the DR- (0) or the MDR-scheme (1)
```
All information which has been gathered during the calculation will be stored in a `HierarchyObject` and can be accessed by member functions. To obtain the 3-loop correction to the Higgs mass matrix you have to call:
```cpp
auto dMh3L = ho.getDMh(3); // this returns a 2x2 matrix which contains the alpha_t*alpha_s^2 corrections for the given parameter point
```
The returned matrix should be added to your two-loop results **before** diagonalization.

A full and detailed example can be found in `source/example.cpp` with its executable `example`.

## Code Documentation
Doxygen can be used to generate code documentation. Go to the `doc` directory
and run
```
doxygen himalaya.conf
```
to generate `html/index.html` and a LaTeX version.
