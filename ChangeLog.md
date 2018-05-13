Himalaya 2.0.0
==============

New features
------------

 * The 3-loop O(at*as^2) contribution to the CP-even Higgs mass matrix
   in the MSSM can now also be calulated in the DR'-bar scheme, given
   a set of MSSM DR'-bar input parameters.  This is the new default.

 * The 3-loop O(at*as^2) contribution to quartic Higgs coupling λ of
   the Standard Model can be calculated in the DR'-bar and MS-bar
   scheme.


Changes
-------

 * The default renormalization scheme has been changed to DR'-bar.
   The used renormalization scheme is specified by the second
   (optional) parameter of the `calculateDMh3L()` function.


Himalaya 1.0.1
==============

Improved debug output.


Himalaya 1.0
============

This is the initial release of the Himalaya library.  The Himalaya
library calculates the 3-loop O(at*as^2) contribution to the CP-even
Higgs mass matrix in the MSSM in the DR-bar or MDR-bar scheme, given a
set of MSSM DR-bar input parameters.
