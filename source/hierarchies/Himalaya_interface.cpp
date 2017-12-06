#include "Himalaya_interface.hpp"
#include <iostream>
#include <vector>
#include <Eigen/Eigenvalues>

namespace himalaya {

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
      const double beta = atan(vu / vd);
      const double Xt = Mt * (At - mu * 1 / tan(beta));
      const double sw2 = 1 - MW * MW / MZ / MZ;
      Eigen::Matrix2d stopMatrix;
      stopMatrix << mq2(2, 2) + Mt * Mt + (1/2. - 2/3. * sw2) * MZ * MZ * cos(2 * beta), Xt,
         Xt, mu2(2, 2) + Mt * Mt + 2 / 3. * sw2 * MZ * MZ * cos(2 * beta);
      // solve eigenvalues and sort them
      const Eigen::EigenSolver<Eigen::Matrix2d> es(stopMatrix);
      std::vector<double> sortedEigenvalues = {sqrt(std::real(es.eigenvalues()(0))), sqrt(std::real(es.eigenvalues()(1)))};
      std::sort(sortedEigenvalues.begin(), sortedEigenvalues.end());
      // set stop masses
      MSt << sortedEigenvalues.at(0), sortedEigenvalues.at(1);

      // extract mixing angle
      const double theta = [&] {
         const double delta1 =
            std::abs(acos(std::real(es.eigenvectors().col(0)(0))) +
                     asin(std::real(es.eigenvectors().col(0)(1))));
         const double delta2 =
            std::abs(acos(std::real(es.eigenvectors().col(1)(0))) +
                     asin(std::real(es.eigenvectors().col(1)(1))));

         if (delta1 < delta2) {
            return acos(std::real(es.eigenvectors().col(0)(0)));
         } else {
            return -acos(std::real(es.eigenvectors().col(1)(0)));
         }
      }();

      s2t = sin(2 * theta);

      if (verbose) {
         std::cout << "\033[1;34m Info:\033[0m Stop masses or mixing angle not provided. Calculated values:\n" <<
            "\tstop masses: " << MSt(0) << " GeV, " << MSt(1) << " GeV,\n" <<
            "\tmixing angle: " << theta << ".\n";
      }
   }

   if (std::isnan(MSb(0)) || std::isnan(MSb(1)) || std::isnan(s2b)) {
      const double beta = atan(vu / vd);
      const double Xb = Mb * (Ab - mu * tan(beta));
      const double sw2 = 1 - MW * MW / MZ / MZ;
      Eigen::Matrix2d sbottomMatrix;
      sbottomMatrix << mq2(2, 2) + Mb * Mb - (1/2. - 1/3. * sw2) * MZ * MZ * cos(2 * beta), Xb,
         Xb, md2(2, 2) + Mb * Mb - 1/3. * sw2 * MZ * MZ * cos(2 * beta);
      // solve eigenvalues and sort them
      const Eigen::EigenSolver<Eigen::Matrix2d> es(sbottomMatrix);
      std::vector<double> sortedEigenvalues = {sqrt(std::real(es.eigenvalues()(0))), sqrt(std::real(es.eigenvalues()(1)))};
      std::sort(sortedEigenvalues.begin(), sortedEigenvalues.end());
      // set sbottom masses
      MSb << sortedEigenvalues.at(0), sortedEigenvalues.at(1);

      // extract mixing angle
      const double theta = [&] {
         const double delta1 =
            std::abs(acos(std::real(es.eigenvectors().col(0)(0))) +
                     asin(std::real(es.eigenvectors().col(0)(1))));
         const double delta2 =
            std::abs(acos(std::real(es.eigenvectors().col(1)(0))) +
                     asin(std::real(es.eigenvectors().col(1)(1))));

         if (delta1 < delta2) {
            return acos(std::real(es.eigenvectors().col(0)(0)));
         } else {
            return -acos(std::real(es.eigenvectors().col(1)(0)));
         }
      }();

      s2b = sin(2 * theta);

      if (verbose) {
         std::cout << "\033[1;34m Info:\033[0m Sbottom masses or mixing angle not provided. Calculated values:\n" <<
            "\tsbottom masses: " << MSb(0) << " GeV, " << MSb(1) << " GeV,\n" <<
            "\tmixing angle: " << theta << ".\n";
      }
   }

   // check the ordering of the stop/sbottom quarks
   if (MSt(0) > MSt(1)) {
      std::swap(MSt(0), MSt(1));
      s2t *= -1;
   }

   if (MSb(0) > MSb(1)) {
      std::swap(MSb(0), MSb(1));
      s2b *= -1;
   }

   // check if the stop/sbottom masses are degenerated. If this is the case one could get spurious poles
   // in Pietro's code. To avoid this numerical issue we shift the stop/bottom 1 mass by a relative (but small)
   // value.
   if (std::abs(MSt(0) - MSt(1)) < 1.0E-5) {
      MSt(0) = MSt(1) / (1. + 1.0E-5);
   }

   if (std::abs(MSb(0) - MSb(1)) < 1.0E-5) {
      MSb(0) = MSb(0) / (1. + 1.0E-5);
   }

   if (std::abs(mq2(2, 2) - mu2(2, 2)) < 1.0E-5) {
      mq2(2,2) = mu2(2,2) / (1. + 1.0E-5);
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
        << "mq2   = {{" << pars.mq2.row(0) << "}, {" << pars.mq2.row(1) << "}, {" << pars.mq2.row(2) << "}}\n"
        << "md2   = {{" << pars.md2.row(0) << "}, {" << pars.md2.row(1) << "}, {" << pars.md2.row(2) << "}}\n"
        << "mu2   = {{" << pars.mu2.row(0) << "}, {" << pars.mu2.row(1) << "}, {" << pars.mu2.row(2) << "}}\n"
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
