// ====================================================================
// This file is part of Himalaya.
//
// Himalaya is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================

#include "himalaya/Himalaya_interface.hpp"

#include "himalaya/misc/Constants.hpp"
#include "himalaya/misc/Logger.hpp"
#include "himalaya/misc/Powers.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>

#include <Eigen/Eigenvalues>

/**
 * @file Himalaya_interface.cpp
 *
 * @brief Definition of the helper functions for the MSSM input
 * parameters
 */

namespace himalaya {

namespace {

int sign(double x) noexcept
{
   return x >= 0. ? 1 : -1;
}

double calc_cw2(double mW, double mZ) noexcept
{
   const double cw2 = pow2(mW/mZ);
   return std::isfinite(cw2) ? cw2 : 1.0;
}

double calc_sw2(double mW, double mZ) noexcept
{
   return 1.0 - calc_cw2(mW, mZ);
}

double calc_v2(double vu, double vd) noexcept
{
   return pow2(vu) + pow2(vd);
}

double calc_cos2beta(double vu, double vd) noexcept
{
   return (pow2(vd) - pow2(vu))/calc_v2(vu, vd);
}

/// sorts two eigenvalues
void sort_ew(V2& ew, double& theta) noexcept
{
   if (ew(0) > ew(1)) {
      std::swap(ew(0), ew(1));
      theta += M_PI/2.;
   }
}

/// calculates mass eigenvalues (in GeV), sin(2*theta) and theta
std::tuple<V2,double,double> calculate_MSf_s2f(const RM22& M)
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

   return std::make_tuple(ew, std::sin(2 * theta), theta);
}

RM33 h_svd(const RM33& M)
{
   Eigen::JacobiSVD<RM33> svd(M);
   return svd.singularValues().reverse().asDiagonal();
}

} // anonymous namespace

double Parameters::calculateMsq2() const
{
   const double cos_2beta = calc_cos2beta(vu, vd);
   const double sw2 = calc_sw2(MW, MZ);
   const double msq = std::pow(
        pow2(mq2(0,0))
       *pow2(mq2(1,1))
       *mq2(0,0)
       *md2(0,0)
       *mu2(1,1)
       *md2(1,1)
       *(mq2(2, 2) + pow2(Mb) - (0.5 - sw2/3) * pow2(MZ) * cos_2beta)
       *(md2(2, 2) + pow2(Mb) - sw2/3 * pow2(MZ) * cos_2beta),
       0.05);

   return pow2(msq);
}

/**
 * Checks if the stop/sbottom masses and mixing angles are provided. Otherwise calculate them.
 * Checks if the stop/sbottom masses are ordered in the right way. If these masses are wrongly ordered
 * the right ordering will be introduced.
 * Checks if the stops/sbottom masses are degenerated and introduce a small shift to the 1st stop/sbottom mass in this case.
 * @param verbose a bool which suppresses the information of the calculation if set to flase
 */
void Parameters::validate(bool verbose)
{
   // check if soft-breaking parameters are greater than zero
   if (mq2.minCoeff() < 0. || md2.minCoeff() < 0. || mu2.minCoeff() < 0. ||
       ml2.minCoeff() < 0. || me2.minCoeff() < 0.) {
      throw std::runtime_error(
         "Soft-breaking squared sfermion mass parameters "
         "must be greater than zero!");
   }

   // The implemented 2- and 3-loop contributions assumed that M3 > 0.
   // For this reason the gluino phase factor is absent from the
   // expressions (i.e. has always been set to 1) and cannot be
   // adjusted in case M3 < 0.
   if (MG < 0.) {
      throw std::runtime_error("Gluino mass parameter must not be negative!");
   }

   // force gluino mass to be positive
   MG = std::abs(MG);

   // diagonalize all yukawa matrices
   Yu = h_svd(Yu);
   Yd = h_svd(Yd);
   Ye = h_svd(Ye);

   if (std::isnan(vu) || vu <= 0.) {
      throw std::runtime_error("Invalid value of vu given!");
   }

   if (std::isnan(vd) || vd <= 0.) {
      throw std::runtime_error("Invalid value of vd given!");
   }

   // calculate all other masses
   if (std::isnan(MW)) {
      MW = 0.5*g2*std::sqrt(calc_v2(vu, vd));
   }
   if (std::isnan(MZ)) {
      MZ = 0.5*std::sqrt((0.6*pow2(g1) + pow2(g2))*calc_v2(vu, vd));
   }
   if (std::isnan(Mt)) {
      Mt = inv_sqrt2*Yu(2,2)*vu;
   }
   if (std::isnan(Mb)) {
      Mb = inv_sqrt2*Yd(2,2)*vd;
   }
   if (std::isnan(Mtau)) {
      Mtau = inv_sqrt2*Ye(2,2)*vd;
   }

   // check if stop/sbottom masses and/or mixing angles are nan. If so, calculate these quantities.
   if (std::isnan(MSt(0)) || std::isnan(MSt(1)) || std::isnan(s2t) || std::isnan(theta_t)) {
      const double tan_beta = vu / vd;
      const double cos_2beta = calc_cos2beta(vu, vd);
      const double Xt = Mt * (Au(2,2) - mu / tan_beta);
      const double sw2 = calc_sw2(MW, MZ);
      RM22 stopMatrix;
      stopMatrix << mq2(2, 2) + pow2(Mt) + (0.5 - 2*sw2/3) * pow2(MZ) * cos_2beta, Xt,
         Xt, mu2(2, 2) + pow2(Mt) + 2*sw2/3 * pow2(MZ) * cos_2beta;

      std::tie(MSt, s2t, theta_t) = calculate_MSf_s2f(stopMatrix);

      if (verbose) {
         INFO_MSG("Stop masses or mixing angle not provided. Calculated values:\n" <<
                  "\tstop masses: " << MSt(0) << " GeV, " << MSt(1) << " GeV,\n" <<
                  "\tmixing angle sin(2*theta_t) = " << s2t << ", theta_t = " << theta_t);
      }
   }

   if (std::isnan(MSb(0)) || std::isnan(MSb(1)) || std::isnan(s2b) || std::isnan(theta_b)) {
      const double tan_beta = vu / vd;
      const double cos_2beta = calc_cos2beta(vu, vd);
      const double Xb = Mb * (Ad(2,2) - mu * tan_beta);
      const double sw2 = calc_sw2(MW, MZ);
      RM22 sbottomMatrix;
      sbottomMatrix << mq2(2, 2) + pow2(Mb) - (0.5 - sw2/3) * pow2(MZ) * cos_2beta, Xb,
         Xb, md2(2, 2) + pow2(Mb) - sw2/3 * pow2(MZ) * cos_2beta;

      std::tie(MSb, s2b, theta_b) = calculate_MSf_s2f(sbottomMatrix);

      if (verbose) {
         INFO_MSG("Sbottom masses or mixing angle not provided. Calculated values:\n" <<
                  "\tsbottom masses: " << MSb(0) << " GeV, " << MSb(1) << " GeV,\n" <<
                  "\tmixing angle sin(2*theta_b) = " << s2b << ", theta_b = " << theta_b);
      }
   }

   if (std::isnan(MStau(0)) || std::isnan(MStau(1)) || std::isnan(s2tau) || std::isnan(theta_tau)) {
      const double tan_beta = vu / vd;
      const double cos_2beta = calc_cos2beta(vu, vd);
      const double Xtau = Mtau * (Ae(2,2) - mu * tan_beta);
      const double sw2 = calc_sw2(MW, MZ);
      RM22 stauMatrix;
      stauMatrix << ml2(2, 2) + pow2(Mtau) - (0.5 - sw2) * pow2(MZ) * cos_2beta, Xtau,
         Xtau, me2(2, 2) + pow2(Mtau) - sw2 * pow2(MZ) * cos_2beta;

      std::tie(MStau, s2tau, theta_tau) = calculate_MSf_s2f(stauMatrix);

      if (verbose) {
         INFO_MSG("Stau masses or mixing angle not provided. Calculated values:\n" <<
                  "\tstau masses: " << MStau(0) << " GeV, " << MStau(1) << " GeV,\n" <<
                  "\tmixing angle sin(2*theta_tau) = " << s2tau << ", theta_tau = " << theta_tau);
      }
   }

   // sort stops/sbottoms/staus
   sort_ew(MSt, s2t);
   sort_ew(MSb, s2b);
   sort_ew(MStau, s2tau);

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
        << "g1    = " << pars.g1 << '\n'
        << "g2    = " << pars.g2 << '\n'
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
        << "ml2   = {{" << pars.ml2(0,0) << ", " << pars.ml2(0,1) << ", " << pars.ml2(0,2) << "}, "
                    "{" << pars.ml2(1,0) << ", " << pars.ml2(1,1) << ", " << pars.ml2(1,2) << "}, "
                    "{" << pars.ml2(2,0) << ", " << pars.ml2(2,1) << ", " << pars.ml2(2,2) << "}}\n"
        << "me2   = {{" << pars.me2(0,0) << ", " << pars.me2(0,1) << ", " << pars.me2(0,2) << "}, "
                    "{" << pars.me2(1,0) << ", " << pars.me2(1,1) << ", " << pars.me2(1,2) << "}, "
                    "{" << pars.me2(2,0) << ", " << pars.me2(2,1) << ", " << pars.me2(2,2) << "}}\n"
        << "Au    = {{" << pars.Au(0,0) << ", " << pars.Au(0,1) << ", " << pars.Au(0,2) << "}, "
                    "{" << pars.Au(1,0) << ", " << pars.Au(1,1) << ", " << pars.Au(1,2) << "}, "
                    "{" << pars.Au(2,0) << ", " << pars.Au(2,1) << ", " << pars.Au(2,2) << "}}\n"
        << "Ad    = {{" << pars.Ad(0,0) << ", " << pars.Ad(0,1) << ", " << pars.Ad(0,2) << "}, "
                    "{" << pars.Ad(1,0) << ", " << pars.Ad(1,1) << ", " << pars.Ad(1,2) << "}, "
                    "{" << pars.Ad(2,0) << ", " << pars.Ad(2,1) << ", " << pars.Ad(2,2) << "}}\n"
        << "Ae    = {{" << pars.Ae(0,0) << ", " << pars.Ae(0,1) << ", " << pars.Ae(0,2) << "}, "
                    "{" << pars.Ae(1,0) << ", " << pars.Ae(1,1) << ", " << pars.Ae(1,2) << "}, "
                    "{" << pars.Ae(2,0) << ", " << pars.Ae(2,1) << ", " << pars.Ae(2,2) << "}}\n"
        << "Yu    = {{" << pars.Yu(0,0) << ", " << pars.Yu(0,1) << ", " << pars.Yu(0,2) << "}, "
                    "{" << pars.Yu(1,0) << ", " << pars.Yu(1,1) << ", " << pars.Yu(1,2) << "}, "
                    "{" << pars.Yu(2,0) << ", " << pars.Yu(2,1) << ", " << pars.Yu(2,2) << "}}\n"
        << "Yd    = {{" << pars.Yd(0,0) << ", " << pars.Yd(0,1) << ", " << pars.Yd(0,2) << "}, "
                    "{" << pars.Yd(1,0) << ", " << pars.Yd(1,1) << ", " << pars.Yd(1,2) << "}, "
                    "{" << pars.Yd(2,0) << ", " << pars.Yd(2,1) << ", " << pars.Yd(2,2) << "}}\n"
        << "Ye    = {{" << pars.Ye(0,0) << ", " << pars.Ye(0,1) << ", " << pars.Ye(0,2) << "}, "
                    "{" << pars.Ye(1,0) << ", " << pars.Ye(1,1) << ", " << pars.Ye(1,2) << "}, "
                    "{" << pars.Ye(2,0) << ", " << pars.Ye(2,1) << ", " << pars.Ye(2,2) << "}}\n"
        << "M1    = " << pars.M1 << '\n'
        << "M2    = " << pars.M2 << '\n'
        << "MG    = " << pars.MG << '\n'
        << "MW    = " << pars.MW << '\n'
        << "MZ    = " << pars.MZ << '\n'
        << "Mt    = " << pars.Mt << '\n'
        << "Mb    = " << pars.Mb << '\n'
        << "Mtau  = " << pars.Mtau << '\n'
        << "MA    = " << pars.MA << '\n'
        << "MSt   = {" << pars.MSt(0) << ", " << pars.MSt(1) << "}\n"
        << "MSb   = {" << pars.MSb(0) << ", " << pars.MSb(1) << "}\n"
        << "MStau = {" << pars.MStau(0) << ", " << pars.MStau(1) << "}\n"
        << "s2t   = " << pars.s2t << '\n'
        << "s2b   = " << pars.s2b << '\n'
        << "s2tau = " << pars.s2tau << '\n';

   return ostr;
}

} // namespace himalaya
