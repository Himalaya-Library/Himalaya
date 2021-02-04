Himalaya 4.1.0
==============

Changes
-------

 * [commit 9a6b60b]: Performance improvement by ~25% by replacing the
   old FORTRAN implementation of the complex dilogarithm (dating back
   to 20.07.83, written by Wolfgang Hollik and Ansgar Denner) by a
   more performant one from the [polylogarithm package
   6.0.0](https://github.com/Expander/polylogarithm).

 * [commit 1251eba]: Performance improvement by ~10% by replacing the
   complex logarithm by a faster implementation.

 * [commit 313d06d]: Require a C++-14 compatible compiler.

 * The source code directory has been slightly restructured to match
   the namespace hierarchy.

Fixed bugs
----------

 * [commit 2b7fdee]: In the calculation of the average light squark
   mass the left-handed sup was counted twice, while the right-handed
   sup was not counted.


Himalaya 4.0.0
==============

Changes
-------

 * The public Himalaya headers have been moved to the
   `include/himalaya/` directory. When installed, the public headers
   are now located in a `himalaya/` sub-directory of include
   directory. As a consequence, the headers have to be included via

       #include "himalaya/HierarchyCalculator.hpp"

 * [commit 03022da]: Faster implementation of 2-loop function
   phi(x,y,z) [Davydychev and Tausk, Nucl. Phys. B397 (1993) 23] for
   degenerate mass spectra.

 * [commit eb52c7b]: Faster implementation of the dilogarithm function.

 * [commit d957fd7]: Update `FindMathematica.cmake` to version 3.3.0.


Himalaya 3.0.1
==============

Changes
-------

 * The threshold correction loop functions and the dilogarithm
   function have been replaced by cleaner and faster implementations.

Fixed bugs
----------

 * [commit 4f063f1]: Fix linking error on MacOS.


Himalaya 3.0.0
==============

New features
------------

 * The 1-, 2- and 3-loop corrections of O(αt), O(αt\*αs + αt^2) and
   O(αt\*αs^2) to the squared light CP-even Higgs mass are now
   available in the EFT limit v^2 << MSUSY^2.

   At the C++ level one can use the following new accessors of the
   `HierarchyObject` class:

    * `getDMh2EFTAt(int loops)`: Returns loop corrections in the EFT
      limit v^2 << MSUSY^2.
    * `getDMh2FOAt(int loops)`: Returns loop corrections in the full
      MSSM.

   At the Mathematica level the following new symbols are provided
   in the output of `HimalayaCalculateDMh3L[]`:

    * `Mh2EFTAt`: Loop corrections in the EFT limit v^2 << MSUSY^2.
    * `Mh2FOAt`: Loop corrections in the full MSSM.

 * The full 1-loop, 2-loop O((αt+ab)\*αs + (αt+αb)^2 + ab\*aτ + aτ^2)
   and 3-loop O(αt\*αs^2) to the squared light CP-even Higgs mass are
   now available.

   At the C++ level one can use the new accessor `getDMh2FO(int
   loops)` of the `HierarchyObject` class.

   At the Mathematica level the new symbol `Mh2FO` is provided in the
   output of `HimalayaCalculateDMh3L[]`.

Changes
-------

 * The getter to access the loop corrections to the light CP-even
   Higgs mass in the EFT limit has been renamed `getDMh2EFT` ->
   `getDMh2EFTAt`. The new name reflects that fact that only
   αt-enhanced contributions are returned.

 * The Mathematica symbol which contains the loop corrections to the
   light CP-even Higgs mass in the EFT limit has been renamed `Mh2EFT`
   -> `Mh2EFTAt`. The new name reflects that fact that only
   αt-enhanced contributions are included.


Himalaya 2.0.2
==============

Changes
-------

 * [commit 462c43c]: An exception is thrown if the (input) gluino mass
   parameter is negative.

   The implemented 2- and 3-loop contributions usually assume that
   M3 > 0. For this reason the gluino phase factor is absent from the
   expressions (i.e. has always been set to 1) and cannot be adjusted
   in case M3 < 0 to make the gluino mass positive.

Fixed bugs
----------

 * [commit d7bffa7]: Update `FindMathematica.cmake` to version 3.2.5,
   which fixes a build error due to missing libuuid.


Himalaya 2.0.1
==============

Changes
-------

 * The estimate of the non-logarithmic contributions to the quartic
   Higgs coupling λ due to the hierarchy expansion has been changed:
   Xt^5 and Xt^6 terms are no longer taken into account in the
   uncertainty estimate, because they do not appear at 3-loop
   O(αt*αs^2).  Only the Xt^4 term is used in the uncertainty
   estimate, if it is absent from the selected hierarchy.  As a
   result, the uncertainty estimate is now symmetric in Xt.


Himalaya 2.0.0
==============

New features
------------

 * The 3-loop O(αt*αs^2) contribution to the CP-even Higgs mass matrix
   in the MSSM is now calulated in the DR'-bar scheme by default.
   Shifts to the MDR'-bar and H3m schemes are provided separately.

 * The 3-loop O(αt*αs^2) contribution to quartic Higgs coupling λ of
   the Standard Model is calculated in the DR'-bar scheme.  A shift to
   the MS-bar scheme is provided separately.

 * A Mathematica interface for Himalaya has been added.  See README.md
   for a usage example.


Changes
-------

 * The default renormalization scheme for the input parameters has
   been changed to DR'-bar.

 * The names of the following member functions of the
   `HierarchyObject` class have been changed:

    * `getExpUncertainty()` -> `getDMhExpUncertainty()`
    * `getDRToMDRShift()` -> `getDMhDRbarPrimeToMDRbarPrimeShift()`

Himalaya 1.0.1
==============

Improved debug output.


Himalaya 1.0
============

This is the initial release of the Himalaya library.  The Himalaya
library calculates the 3-loop O(at*as^2) contribution to the CP-even
Higgs mass matrix in the MSSM in the DR-bar or MDR-bar scheme, given a
set of MSSM DR-bar input parameters.
