#include "Himalaya_interface.hpp"
#include "Logger.hpp"
#include "Utils.hpp"
#include <Eigen/Eigenvalues>
#include <cmath>
#include <iostream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace himalaya {

namespace {

int sign(double x) noexcept
{
   return x >= 0. ? 1 : -1;
}

double sqr(double x) noexcept
{
   return x*x;
}

/// sorts two eigenvalues
void sort_ew(V2& ew, double& theta) noexcept
{
   if (ew(0) > ew(1)) {
      std::swap(ew(0), ew(1));
      theta *= -1;
   }
}

/// calculates mass eigenvalues (in GeV) and sin(2*theta)
std::pair<V2,double> calculate_MSf_s2f(const RM22& M)
{
   const Eigen::SelfAdjointEigenSolver<RM22> es(M);
   RM22 ev = es.eigenvectors().real();
   V2 ew = es.eigenvalues();

   if (ew.minCoeff() < 0.) {
      throw std::runtime_error(
         "DR sfermion masses are tachyonic: mst1^2 = " + std::to_string(ew(0))
         + ", mst2^2 = " + std::to_string(ew(1)));
   }

   ew = ew.unaryExpr([](double x){ return std::sqrt(x); });

   // extract mixing angle in the convention of hep-ph/0105096
   double theta = 0.;

   if (sign(ev(0,0)) == sign(ev(1,1))) {
      theta = std::acos(std::abs(ev(0,0)));
   } else {
      theta = std::acos(std::abs(ev(0,1)));
      ev.col(0).swap(ev.col(1));
      std::swap(ew(0), ew(1));
   }

   theta = sign(M(0,1) / (M(0,0) - M(1,1))) * std::abs(theta);

   sort_ew(ew, theta);

   const double s2t = std::sin(2 * theta);

   return std::make_pair(ew, s2t);
}

} // anonymous namespace

double Parameters::calculateMsq2() const
{
   using std::sqrt;

   const double beta = std::atan(vu / vd);
   const double cos_2beta = std::cos(2 * beta);
   const double sw2 = 1 - pow2(MW / MZ);
   const double msq =
      (+ 2 * sqrt(mq2(0, 0)) + sqrt(mu2(0, 0)) + sqrt(md2(0, 0))  // sup and sdown
       + 2 * sqrt(mq2(1, 1)) + sqrt(mu2(1, 1)) + sqrt(md2(1, 1))  // scharm and sstrange
       // sbottom
       + sqrt(mq2(2, 2) + pow2(Mb) - (1 / 2. - 1 / 3. * sw2) * pow2(MZ) * cos_2beta)
       + sqrt(md2(2, 2) + pow2(Mb) - 1 / 3. * sw2 * pow2(MZ) * cos_2beta)) / 10.;

   return pow2(msq);
}

/**
 * 	Checks if the stop/sbottom masses and mixing angles are provided. Otherwise calculate them.
 * 	Checks if the stop/sbottom masses are ordered in the right way. If these masses are wrongly ordered
 * 	the right ordering will be introduced.
 * 	Checks if the stops/sbottom masses are degenerated and introduce a small shift to the 1st stop/sbottom mass in this case.
 * 	@param verbose a bool which suppresses the information of the calculation if set to flase
 */
void Parameters::validate(bool verbose)
{
   // check if stop/sbottom masses and/or mixing angles are nan. If so, calculate these quantities.
   if (std::isnan(MSt(0)) || std::isnan(MSt(1)) || std::isnan(s2t)) {
      const double tan_beta = vu / vd;
      const double beta = std::atan(tan_beta);
      const double cos_2beta = std::cos(2 * beta);
      const double Xt = Mt * (At - mu / tan_beta);
      const double sw2 = 1 - MW * MW / MZ / MZ;
      RM22 stopMatrix;
      stopMatrix << mq2(2, 2) + sqr(Mt) + (1/2. - 2/3. * sw2) * sqr(MZ) * cos_2beta, Xt,
	Xt, mu2(2, 2) + sqr(Mt) + 2 / 3. * sw2 * sqr(MZ) * cos_2beta;

      std::tie(MSt, s2t) = calculate_MSf_s2f(stopMatrix);

      if (verbose) {
         INFO_MSG("Stop masses or mixing angle not provided. Calculated values:\n" <<
                  "\tstop masses: " << MSt(0) << " GeV, " << MSt(1) << " GeV,\n" <<
                  "\tmixing angle sin(2*theta): " << s2t);
      }
   }

   if (std::isnan(MSb(0)) || std::isnan(MSb(1)) || std::isnan(s2b)) {
      const double tan_beta = vu / vd;
      const double beta = std::atan(tan_beta);
      const double cos_2beta = std::cos(2 * beta);
      const double Xb = Mb * (Ab - mu * tan_beta);
      const double sw2 = 1 - MW * MW / MZ / MZ;
      RM22 sbottomMatrix;
      sbottomMatrix << mq2(2, 2) + sqr(Mb) - (1/2. - 1/3. * sw2) * sqr(MZ) * cos_2beta, Xb,
         Xb, md2(2, 2) + sqr(Mb) - 1/3. * sw2 * sqr(MZ) * cos_2beta;

      std::tie(MSb, s2b) = calculate_MSf_s2f(sbottomMatrix);

      if (verbose) {
         INFO_MSG("Sbottom masses or mixing angle not provided. Calculated values:\n" <<
                  "\tsbottom masses: " << MSb(0) << " GeV, " << MSb(1) << " GeV,\n" <<
                  "\tmixing angle sin(2*theta): " << s2b << ".");
      }
   }

   // sort stops/sbottoms
   sort_ew(MSt, s2t);
   sort_ew(MSb, s2b);

   // check if the stop/sbottom masses are degenerated. If this is the
   // case one could get spurious poles in Pietro's code. To avoid
   // this numerical issue we shift the stop/bottom 1 mass by a
   // relative (but small) value.
   const double eps = 1e-8;

   if (std::abs(MSt(0) - MSt(1)) < eps) {
      MSt(0) = MSt(1) / (1. + eps);
   }

   if (std::abs(MSb(0) - MSb(1)) < eps) {
      MSb(0) = MSb(0) / (1. + eps);
   }

   if (std::abs(mq2(2, 2) - mu2(2, 2)) < eps) {
      mq2(2,2) = mu2(2,2) / (1. + eps);
   }
}

std::ostream& operator<<(std::ostream& ostr, const Parameters& pars)
{
   ostr << "Himalaya parameters\n"
        << "===================\n"
        << "scale = " << pars.scale << '\n'
        << "mu    = " << pars.mu << '\n'
        << "g3    = " << pars.g3 << '\n'
        << "vd    = " << pars.vd << '\n'
        << "vu    = " << pars.vu << '\n'
        << "mq2   = {{" << pars.mq2(0,0) << ", " << pars.mq2(0,1) << ", " << pars.mq2(0,2) << "}, "
                    "{" << pars.mq2(1,0) << ", " << pars.mq2(1,1) << ", " << pars.mq2(1,2) << "}, "
                    "{" << pars.mq2(2,0) << ", " << pars.mq2(2,1) << ", " << pars.mq2(2,2) << "}}\n"
        << "md2   = {{" << pars.md2(0,0) << ", " << pars.md2(0,1) << ", " << pars.md2(0,2) << "}, "
                    "{" << pars.md2(1,0) << ", " << pars.md2(1,1) << ", " << pars.md2(1,2) << "}, "
                    "{" << pars.md2(2,0) << ", " << pars.md2(2,1) << ", " << pars.md2(2,2) << "}}\n"
        << "mu2   = {{" << pars.mu2(0,0) << ", " << pars.mu2(0,1) << ", " << pars.mu2(0,2) << "}, "
                    "{" << pars.mu2(1,0) << ", " << pars.mu2(1,1) << ", " << pars.mu2(1,2) << "}, "
                    "{" << pars.mu2(2,0) << ", " << pars.mu2(2,1) << ", " << pars.mu2(2,2) << "}}\n"
        << "At    = " << pars.At << '\n'
        << "Ab    = " << pars.Ab << '\n'
        << "MG    = " << pars.MG << '\n'
        << "MW    = " << pars.MW << '\n'
        << "MZ    = " << pars.MZ << '\n'
        << "Mt    = " << pars.Mt << '\n'
        << "Mb    = " << pars.Mb << '\n'
        << "MA    = " << pars.MA << '\n'
        << "MSt   = {" << pars.MSt(0) << ", " << pars.MSt(1) << "}\n"
        << "MSb   = {" << pars.MSb(0) << ", " << pars.MSb(1) << "}\n"
        << "s2t   = " << pars.s2t << '\n'
        << "s2b   = " << pars.s2b << '\n';

   return ostr;
}

} // namespace himalaya
