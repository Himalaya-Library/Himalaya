Himalaya 3.0.0
==============

New features
------------

 * The 1-, 2- and 3-loop corrections of O(αt), O(αt\*αs) and
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
