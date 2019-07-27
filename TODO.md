# Todo list for Himalaya 2.1.0

## Interfaces

* decide on C++ interface for the EFT/FO calculation

* decide on Mathematica interface for the EFT/FO calculation

* Our `example.cpp` should only include the one Himalaya header:

      #include "HierarchyCalculator.hpp"

  This should fixed when the decision on the C++ interface is made.

* update references for 2-loop corrections in the output `operator<<`
  of the `HierarchyCalculator`.  See the references in `README.md`.
