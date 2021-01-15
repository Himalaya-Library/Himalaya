// ====================================================================
// This file is part of Himalaya.
//
// Himalaya is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================

#include "li2.hpp"
#include <complex>

/// complex dilogarithm from FORTRAN module
extern "C" int li2c_(double*, double*, double*, double*);

namespace himalaya {

/// real dilogarithm from FORTRAN module
double li2(double x)
{
   const std::complex<double> z(x, 0.);
   return std::real(li2(z));
}

/// complex dilogarithm from FORTRAN module
std::complex<double> li2(const std::complex<double>& z)
{
   double re_in = z.real();
   double im_in = z.imag();
   double re_out = 0.;
   double im_out = 0.;

   li2c_(&re_in, &im_in, &re_out, &im_out);

   return { re_out, im_out };
}

} // namespace himalaya
