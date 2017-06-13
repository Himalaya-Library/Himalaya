#Readme of Himalaya 1.0.0

#Installation
To build Himalaya
	
  1) Check your c++ compiler. Mac users should use g++ and edit the corresponding line in CMakeLists.txt

  2) Just call
    mkdir build
    cd build
    cmake ../
    make
	   
  3) Link the static library libhimalaya.a to your project

#Usage
The Himalaya library provides these functions:

  calculateDMh3L
  compareHierarchies
  calculateHierarchy
  calcDRbarToMDRbarShift
  getMt41L
  getMt42L
  shiftMst1ToMDR
  shiftMst2ToMDR
  getExpansionUncertainty
  checkTerms

They are documented in the doc/ directory. To generate the documentation type

     cd doc/
     doxygen himalaya.conf
     