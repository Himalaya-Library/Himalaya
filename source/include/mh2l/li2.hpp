// ====================================================================
// This file is part of Himalaya.
//
// Himalaya is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================

#include <complex>

namespace himalaya {

/// real dilogarithm from FORTRAN module
double li2(double);

/// complex dilogarithm from FORTRAN module
std::complex<double> li2(const std::complex<double>&);

} // namespace himalaya
