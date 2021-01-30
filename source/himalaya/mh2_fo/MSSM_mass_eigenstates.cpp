// ====================================================================
// This file is part of Himalaya.
//
// Himalaya is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================

#include "himalaya/mh2_fo/MSSM_mass_eigenstates.hpp"
#include "himalaya/mh2_fo/Linalg.hpp"
#include "himalaya/mh2_fo/PV.hpp"
#include "himalaya/mh2_fo/Sum.hpp"
#include "himalaya/mh2l/DSZHiggs.hpp"
#include "himalaya/misc/Constants.hpp"
#include "himalaya/misc/CouplingOrders.hpp"
#include "himalaya/misc/Logger.hpp"
#include "himalaya/misc/Numerics.hpp"
#include "himalaya/misc/Powers.hpp"
#include <cmath>
#include <complex>
#include <iostream>
#include <limits>
#include <tuple>
#include <type_traits>

/**
 * @file MSSM_mass_eigenstates.cpp
 *
 * @brief Contains the implementation of the \a MSSM_spectrum and \a
 * MSSM_mass_eigenstates classes.
 */

namespace himalaya {
namespace mh2_fo {

namespace {

template <typename T>
T constexpr sqr(T x) noexcept { return x*x; }

} // anonymous namespace

/**
 * Constructor
 *
 * @param pars_ MSSM DR' parameters
 * @param only_at if true, only alpha_t-enhanced contributions are calculated
 */
MSSM_mass_eigenstates::MSSM_mass_eigenstates(const Parameters& pars_, bool only_at)
   : pars(pars_)
   , masses(pars_)
   , gaugeless(make_gaugeless(pars_))
{
   orders.fill(1);

   const double eps = 1e-10;

   // TODO(voigt) check if they are still consistent
   if (std::abs(pars_.Mt) < eps) {
      orders.at(CouplingOrders::G32YT4) = 0;
      orders.at(CouplingOrders::YT6) = 0;
      orders.at(CouplingOrders::YB6) = 0;
   }

   if (std::abs(pars_.Mb) < eps) {
      orders.at(CouplingOrders::G32YB4) = 0;
      orders.at(CouplingOrders::YT6) = 0;
      orders.at(CouplingOrders::YB6) = 0;
      orders.at(CouplingOrders::YTAU2YB4) = 0;
      orders.at(CouplingOrders::YTAU4YB2) = 0;
   }

   if (std::abs(pars_.Mtau) < eps) {
      orders.at(CouplingOrders::YTAU6) = 0;
      orders.at(CouplingOrders::YTAU2YB4) = 0;
      orders.at(CouplingOrders::YTAU4YB2) = 0;
   }

   if (only_at) {
      orders.at(CouplingOrders::YB6) = 0;
      orders.at(CouplingOrders::YTAU6) = 0;
      orders.at(CouplingOrders::YTAU2YB4) = 0;
      orders.at(CouplingOrders::YTAU4YB2) = 0;
      orders.at(CouplingOrders::G32YB4) = 0;

      // re-calculate DR' masses for g1 = g2 = yb = ytau = 0
      pars = make_gaugeless(pars_);
      masses = make_3rd_gen(pars);
      masses.M2hh(0) = 0.;
   } else {
     orders.at(CouplingOrders::ONLY_AT_AS) = 0;
   }
}

/**
 * Transform parameter point to be gaugeless, g1 = g2 = 0.
 */
Parameters MSSM_mass_eigenstates::make_gaugeless(const Parameters& pars)
{
   auto gl = pars;
   gl.g1 = 0.;
   gl.g2 = 0.;
   return gl;
}

/// sets yb = ytau = 0
Parameters MSSM_mass_eigenstates::make_3rd_gen(const Parameters& pars)
{
   const double eps = 1e-6;

   auto gen3 = pars;

   gen3.Mb = NaN;
   gen3.Mtau = NaN;
   gen3.MSb.setConstant(NaN);
   gen3.MStau.setConstant(NaN);
   gen3.s2b = NaN;
   gen3.s2tau = NaN;
   gen3.theta_b = NaN;
   gen3.theta_tau = NaN;

   for (int i = 0; i < 3; i++) {
      for (int k = 0; k < 3; k++) {
         if (i == 2 && k == 2) {
            gen3.Yd(i,k) = eps;
            gen3.Ye(i,k) = eps;
         } else {
            gen3.Yd(i,k) = 0.;
            gen3.Ye(i,k) = 0.;
            gen3.Yu(i,k) = 0.;
         }
      }
   }

   gen3.validate(false);

   return gen3;
}

/**
 * Returns the tree-level squared Higgs masses.
 * @return squared Higgs masses
 */
V2 MSSM_mass_eigenstates::calculate_Mh2_tree() const
{
   const auto m0 = get_mass_matrix_hh();
   const auto mh2_tree = fs_diagonalize_hermitian_perturbatively(m0);

   return std::get<0>(mh2_tree);
}

/**
 * Returns the tree-level, 1- and 2-loop contribution to the squared
 * light CP-even Higgs pole mass.  The function does not include
 * implicit or explicit higher orders.
 *
 * @return Tree-level, 1- and 2-loop contribution to light CP-even Higgs mass
 */
std::tuple<double,double,double> MSSM_mass_eigenstates::calculate_Mh2() const
{
   const auto p2    = calculate_Mh2_tree()(0);
   const auto m0    = get_mass_matrix_hh();
   const auto m1    = delta_mh2_1loop(p2);
   const auto m2    = delta_mh2_2loop();

   if (diagonalization == Diagonalization::pert) {
      const auto m0_gl = get_mass_matrix_hh_gaugeless();
      const auto m1_gl = delta_mh2_1loop_gaugeless();

      // tree-level and 1-loop with electroweak gauge couplings
      const auto Mh2 = fs_diagonalize_hermitian_perturbatively(m0, m1);

      // 2-loop in gaugless limit (p = g1 = g2 = 0)
      const auto Mh2_gl = fs_diagonalize_hermitian_perturbatively(m0_gl, m1_gl, m2);

      return std::make_tuple(std::get<0>(Mh2)(0),
                             std::get<1>(Mh2)(0),
                             std::get<2>(Mh2_gl)(0));
   }

   // numeric diagonalization
   RM22 ZH;
   MSSM_spectrum::A2 M2hh_0L, M2hh_1L, M2hh_2L;
   const RM22& MHH_0L = m0;
   const RM22 MHH_1L = m0 + m1;
   const RM22 MHH_2L = m0 + m1 + m2;

   fs_diagonalize_hermitian(MHH_0L, M2hh_0L, ZH);
   fs_diagonalize_hermitian(MHH_1L, M2hh_1L, ZH);
   fs_diagonalize_hermitian(MHH_2L, M2hh_2L, ZH);

   return std::make_tuple(M2hh_0L(0), (M2hh_1L - M2hh_0L)(0), (M2hh_2L - M2hh_1L)(0));
}

/**
 * Higgs 1-loop DR' contribution for arbitrary momentum.
 * The function returns 1/(4Pi)^2 (-selfenergy + tadpole).
 *
 * @param p2 squared momentum
 *
 * @return 1-loop contribution
 */
RM22 MSSM_mass_eigenstates::delta_mh2_1loop(double p2) const
{
   const auto g1     = pars.g1;
   const auto g2     = pars.g2;
   const auto yt     = pars.Yu(2,2);
   const auto yb     = pars.Yd(2,2);
   const auto ytau   = pars.Ye(2,2);
   const auto vu     = pars.vu;
   const auto vd     = pars.vd;
   const auto mu     = pars.mu;
   const auto At     = pars.Au(2,2);
   const auto Ab     = pars.Ad(2,2);
   const auto Atau   = pars.Ae(2,2);
   const auto M2VWm  = masses.M2VWm;
   const auto M2VZ   = masses.M2VZ;
   const auto MFt    = masses.MFt;
   const auto MFb    = masses.MFb;
   const auto MFtau  = masses.MFtau;
   const auto M2SveL = masses.M2SveL;
   const auto M2SvmL = masses.M2SvmL;
   const auto M2SvtL = masses.M2SvtL;
   const auto M2Su   = masses.M2Su;
   const auto M2Sd   = masses.M2Sd;
   const auto M2Sc   = masses.M2Sc;
   const auto M2Ss   = masses.M2Ss;
   const auto M2St   = masses.M2St;
   const auto M2Sb   = masses.M2Sb;
   const auto M2Se   = masses.M2Se;
   const auto M2Sm   = masses.M2Sm;
   const auto M2Stau = masses.M2Stau;
   const auto M2hh   = masses.M2hh;
   const auto M2Ah   = masses.M2Ah;
   const auto M2Hpm  = masses.M2Hpm;
   const auto MChi   = masses.MChi;
   const auto MCha   = masses.MCha;
   const auto ZU     = masses.ZU;
   const auto ZD     = masses.ZD;
   const auto ZC     = masses.ZC;
   const auto ZS     = masses.ZS;
   const auto ZT     = masses.ZT;
   const auto ZB     = masses.ZB;
   const auto ZE     = masses.ZE;
   const auto ZM     = masses.ZM;
   const auto ZTau   = masses.ZTau;
   const auto ZH     = masses.ZH;
   const auto ZA     = masses.ZA;
   const auto ZP     = masses.ZP;
   const auto ZN     = masses.ZN;
   const auto UM     = masses.UM;
   const auto UP     = masses.UP;

   double se11{0.}, se12{0.}, se22{0.};

   se11 += -(A0(M2VWm)*sqr(g2))/2.;
   se11 += (A0(M2VZ)*(-3*sqr(g1) - 5*sqr(g2)))/20.;
   se11 += -(B0(p2,M2SveL,M2SveL)*sqr(vd*(3*sqr(g1) + 5*sqr(g2))))/400.;
   se11 += -(B0(p2,M2SvmL,M2SvmL)*sqr(vd*(3*sqr(g1) + 5*sqr(g2))))/400.;
   se11 += -(B0(p2,M2SvtL,M2SvtL)*sqr(vd*(3*sqr(g1) + 5*sqr(g2))))/400.;
   se11 += (-7*B0(p2,M2VWm,M2VWm)*pow4(g2)*sqr(vd))/8.;
   se11 += (-7*B0(p2,M2VZ,M2VZ)*sqr(vd*(3*sqr(g1) + 5*sqr(g2))))/400.;
   se11 += 3*B0(p2,sqr(MFb),sqr(MFb))*(-p2 + 4*sqr(MFb))*sqr(yb);
   se11 += B0(p2,sqr(MFtau),sqr(MFtau))*(-p2 + 4*sqr(MFtau))*sqr(ytau);
   se11 += SUM(gI1,0,1,SUM(gI2,0,1,-(B0(p2,M2Ah(gI1),M2Ah(gI2))*sqr(vd*(3*sqr(g1) + 5*sqr(g2))*(ZA(gI1,0)*ZA(gI2,0) - ZA(gI1,1)*ZA(gI2,1))))/800.));
   se11 += SUM(gI1,0,1,SUM(gI2,0,1,(-3*B0(p2,M2Sc(gI1),M2Sc(gI2))*sqr(vd*((sqr(g1) - 5*sqr(g2))*ZC(gI1,0)*ZC(gI2,0) - 4*sqr(g1)*ZC(gI1,1)*ZC(gI2,1))))/400.));
   se11 += SUM(gI1,0,1,SUM(gI2,0,1,(-3*B0(p2,M2Sd(gI1),M2Sd(gI2))*sqr(vd*((sqr(g1) + 5*sqr(g2))*ZD(gI1,0)*ZD(gI2,0) + 2*sqr(g1)*ZD(gI1,1)*ZD(gI2,1))))/400.));
   se11 += SUM(gI1,0,1,SUM(gI2,0,1,-(B0(p2,M2Se(gI1),M2Se(gI2))*sqr(vd*((3*sqr(g1) - 5*sqr(g2))*ZE(gI1,0)*ZE(gI2,0) - 6*sqr(g1)*ZE(gI1,1)*ZE(gI2,1))))/400.));
   se11 += SUM(gI1,0,1,SUM(gI2,0,1,-(B0(p2,M2Sm(gI1),M2Sm(gI2))*sqr(vd*((3*sqr(g1) - 5*sqr(g2))*ZM(gI1,0)*ZM(gI2,0) - 6*sqr(g1)*ZM(gI1,1)*ZM(gI2,1))))/400.));
   se11 += SUM(gI1,0,1,SUM(gI2,0,1,(-3*B0(p2,M2Ss(gI1),M2Ss(gI2))*sqr(vd*((sqr(g1) + 5*sqr(g2))*ZS(gI1,0)*ZS(gI2,0) + 2*sqr(g1)*ZS(gI1,1)*ZS(gI2,1))))/400.));
   se11 += SUM(gI1,0,1,SUM(gI2,0,1,(-3*B0(p2,M2Su(gI1),M2Su(gI2))*sqr(vd*((sqr(g1) - 5*sqr(g2))*ZU(gI1,0)*ZU(gI2,0) - 4*sqr(g1)*ZU(gI1,1)*ZU(gI2,1))))/400.));
   se11 += SUM(gI1,0,1,SUM(gI2,0,1,(sqr(g2)*(A0(sqr(MCha(gI1))) + A0(sqr(MCha(gI2))) + B0(p2,sqr(MCha(gI1)),sqr(MCha(gI2)))*(-p2 + sqr(MCha(gI1)) + sqr(MCha(gI2))))*(sqr(UM(gI2,1)*UP(gI1,0)) + sqr(UM(gI1,1)*UP(gI2,0))))/2.) + MCha(gI1)*(SUM(gI2,0,1,2*B0(p2,sqr(MCha(gI1)),sqr(MCha(gI2)))*MCha(gI2)*sqr(g2)*UM(gI1,1)*UM(gI2,1)*UP(gI1,0)*UP(gI2,0)) - (2*g2*sqrt2*A0(sqr(MCha(gI1)))*UM(gI1,1)*UP(gI1,0))/vd));
   se11 += SUM(gI1,0,1,SUM(gI2,0,1,(-3*B0(p2,M2Sb(gI1),M2Sb(gI2))*sqr(ZB(gI1,0)*(vd*(sqr(g1) + 5*sqr(g2) - 20*sqr(yb))*ZB(gI2,0) - 10*Ab*sqrt2*yb*ZB(gI2,1)) + 2*ZB(gI1,1)*(-5*Ab*sqrt2*yb*ZB(gI2,0) + vd*(sqr(g1) - 10*sqr(yb))*ZB(gI2,1))))/400.) + (3*Ab*sqrt2*yb*A0(M2Sb(gI1))*ZB(gI1,0)*ZB(gI1,1))/vd);
   se11 += SUM(gI1,0,1,-(SUM(gI2,0,1,(vd*B0(p2,M2hh(gI1),M2hh(gI2))*sqr((3*sqr(g1) + 5*sqr(g2))*(ZH(gI1,1)*(vu*ZH(gI2,0) + vd*ZH(gI2,1)) + ZH(gI1,0)*(-3*vd*ZH(gI2,0) + vu*ZH(gI2,1)))))/40.) + vu*A0(M2hh(gI1))*(3*sqr(g1) + 5*sqr(g2))*ZH(gI1,0)*ZH(gI1,1))/(20.*vd));
   se11 += SUM(gI1,0,1,SUM(gI2,0,1,-(B0(p2,M2Hpm(gI1),M2Hpm(gI2))*sqr(ZP(gI1,0)*(vd*(3*sqr(g1) + 5*sqr(g2))*ZP(gI2,0) + 5*vu*sqr(g2)*ZP(gI2,1)) + ZP(gI1,1)*(5*vu*sqr(g2)*ZP(gI2,0) + vd*(-3*sqr(g1) + 5*sqr(g2))*ZP(gI2,1))))/400.) + (vu*A0(M2Hpm(gI1))*sqr(g2)*ZP(gI1,0)*ZP(gI1,1))/(2.*vd));
   se11 += SUM(gI1,0,1,(-3*(SUM(gI2,0,1,(vd*B0(p2,M2St(gI1),M2St(gI2))*sqr(ZT(gI1,0)*(vd*(sqr(g1) - 5*sqr(g2))*ZT(gI2,0) + 10*mu*sqrt2*yt*ZT(gI2,1)) + 2*ZT(gI1,1)*(5*mu*sqrt2*yt*ZT(gI2,0) - 2*vd*sqr(g1)*ZT(gI2,1))))/400.) + mu*sqrt2*yt*A0(M2St(gI1))*ZT(gI1,0)*ZT(gI1,1)))/vd);
   se11 += SUM(gI1,0,1,SUM(gI2,0,1,-(B0(p2,M2Stau(gI1),M2Stau(gI2))*sqr(ZTau(gI1,0)*(vd*(3*sqr(g1) - 5*sqr(g2) + 20*sqr(ytau))*ZTau(gI2,0) + 10*Atau*sqrt2*ytau*ZTau(gI2,1)) + 2*ZTau(gI1,1)*(5*Atau*sqrt2*ytau*ZTau(gI2,0) + vd*(-3*sqr(g1) + 10*sqr(ytau))*ZTau(gI2,1))))/400.) + (Atau*sqrt2*ytau*A0(M2Stau(gI1))*ZTau(gI1,0)*ZTau(gI1,1))/vd);
   se11 += SUM(gI1,0,3,SUM(gI2,0,3,((A0(sqr(MChi(gI1))) + A0(sqr(MChi(gI2))) + B0(p2,sqr(MChi(gI1)),sqr(MChi(gI2)))*(-p2 + 2*MChi(gI1)*MChi(gI2) + sqr(MChi(gI1)) + sqr(MChi(gI2))))*sqr(ZN(gI1,2)*(g1*sqrt15*ZN(gI2,0) - 5*g2*ZN(gI2,1)) + (g1*sqrt15*ZN(gI1,0) - 5*g2*ZN(gI1,1))*ZN(gI2,2)))/100.) + (2*A0(sqr(MChi(gI1)))*MChi(gI1)*(g1*sqrt15*ZN(gI1,0) - 5*g2*ZN(gI1,1))*ZN(gI1,2))/(5.*vd));
   se11 += SUM(gI2,0,1,((2*A0(M2VZ) - A0(M2Ah(gI2)) + B0(p2,M2VZ,M2Ah(gI2))*(-M2VZ + 2*p2 + 2*M2Ah(gI2)))*(3*sqr(g1) + 5*sqr(g2))*sqr(ZA(gI2,0)))/20.);
   se11 += SUM(gI2,0,1,((2*A0(M2VWm) - A0(M2Hpm(gI2)) + B0(p2,M2VWm,M2Hpm(gI2))*(-M2VWm + 2*p2 + 2*M2Hpm(gI2)))*sqr(g2*ZP(gI2,0)))/2.);

   se12 += (vd*vu*B0(p2,M2SveL,M2SveL)*sqr(3*sqr(g1) + 5*sqr(g2)))/400.;
   se12 += (vd*vu*B0(p2,M2SvmL,M2SvmL)*sqr(3*sqr(g1) + 5*sqr(g2)))/400.;
   se12 += (vd*vu*B0(p2,M2SvtL,M2SvtL)*sqr(3*sqr(g1) + 5*sqr(g2)))/400.;
   se12 += (-7*vd*vu*B0(p2,M2VWm,M2VWm)*pow4(g2))/8.;
   se12 += (-7*vd*vu*B0(p2,M2VZ,M2VZ)*sqr(3*sqr(g1) + 5*sqr(g2)))/400.;
   se12 += SUM(gI1,0,1,SUM(gI2,0,1,(sqr(g2)*(2*B0(p2,sqr(MCha(gI1)),sqr(MCha(gI2)))*MCha(gI1)*MCha(gI2)*(UM(gI1,1)*UM(gI2,0)*UP(gI1,1)*UP(gI2,0) + UM(gI1,0)*UM(gI2,1)*UP(gI1,0)*UP(gI2,1)) + (A0(sqr(MCha(gI1))) + A0(sqr(MCha(gI2))) + B0(p2,sqr(MCha(gI1)),sqr(MCha(gI2)))*(-p2 + sqr(MCha(gI1)) + sqr(MCha(gI2))))*(UM(gI2,0)*UM(gI2,1)*UP(gI1,0)*UP(gI1,1) + UM(gI1,0)*UM(gI1,1)*UP(gI2,0)*UP(gI2,1))))/2.));
   se12 += SUM(gI1,0,1,SUM(gI2,0,1,(vd*vu*B0(p2,M2Ah(gI1),M2Ah(gI2))*sqr((3*sqr(g1) + 5*sqr(g2))*(ZA(gI1,0)*ZA(gI2,0) - ZA(gI1,1)*ZA(gI2,1))))/800.));
   se12 += SUM(gI1,0,1,SUM(gI2,0,1,(3*B0(p2,M2Sb(gI1),M2Sb(gI2))*(ZB(gI1,0)*(vu*(sqr(g1) + 5*sqr(g2))*ZB(gI2,0) - 10*mu*sqrt2*yb*ZB(gI2,1)) + 2*ZB(gI1,1)*(-5*mu*sqrt2*yb*ZB(gI2,0) + vu*sqr(g1)*ZB(gI2,1)))*(ZB(gI1,0)*(vd*(sqr(g1) + 5*sqr(g2) - 20*sqr(yb))*ZB(gI2,0) - 10*Ab*sqrt2*yb*ZB(gI2,1)) + 2*ZB(gI1,1)*(-5*Ab*sqrt2*yb*ZB(gI2,0) + vd*(sqr(g1) - 10*sqr(yb))*ZB(gI2,1))))/400.));
   se12 += SUM(gI1,0,1,SUM(gI2,0,1,(3*vd*vu*B0(p2,M2Sc(gI1),M2Sc(gI2))*sqr((sqr(g1) - 5*sqr(g2))*ZC(gI1,0)*ZC(gI2,0) - 4*sqr(g1)*ZC(gI1,1)*ZC(gI2,1)))/400.));
   se12 += SUM(gI1,0,1,SUM(gI2,0,1,(3*vd*vu*B0(p2,M2Sd(gI1),M2Sd(gI2))*sqr((sqr(g1) + 5*sqr(g2))*ZD(gI1,0)*ZD(gI2,0) + 2*sqr(g1)*ZD(gI1,1)*ZD(gI2,1)))/400.));
   se12 += SUM(gI1,0,1,SUM(gI2,0,1,(vd*vu*B0(p2,M2Se(gI1),M2Se(gI2))*sqr((3*sqr(g1) - 5*sqr(g2))*ZE(gI1,0)*ZE(gI2,0) - 6*sqr(g1)*ZE(gI1,1)*ZE(gI2,1)))/400.));
   se12 += SUM(gI1,0,1,SUM(gI2,0,1,(vd*vu*B0(p2,M2Sm(gI1),M2Sm(gI2))*sqr((3*sqr(g1) - 5*sqr(g2))*ZM(gI1,0)*ZM(gI2,0) - 6*sqr(g1)*ZM(gI1,1)*ZM(gI2,1)))/400.));
   se12 += SUM(gI1,0,1,SUM(gI2,0,1,(3*vd*vu*B0(p2,M2Ss(gI1),M2Ss(gI2))*sqr((sqr(g1) + 5*sqr(g2))*ZS(gI1,0)*ZS(gI2,0) + 2*sqr(g1)*ZS(gI1,1)*ZS(gI2,1)))/400.));
   se12 += SUM(gI1,0,1,SUM(gI2,0,1,(3*B0(p2,M2St(gI1),M2St(gI2))*(ZT(gI1,0)*(vd*(sqr(g1) - 5*sqr(g2))*ZT(gI2,0) + 10*mu*sqrt2*yt*ZT(gI2,1)) + 2*ZT(gI1,1)*(5*mu*sqrt2*yt*ZT(gI2,0) - 2*vd*sqr(g1)*ZT(gI2,1)))*(ZT(gI1,0)*(vu*(sqr(g1) - 5*sqr(g2) + 20*sqr(yt))*ZT(gI2,0) + 10*At*sqrt2*yt*ZT(gI2,1)) + 2*ZT(gI1,1)*(5*At*sqrt2*yt*ZT(gI2,0) - 2*vu*(sqr(g1) - 5*sqr(yt))*ZT(gI2,1))))/400.));
   se12 += SUM(gI1,0,1,SUM(gI2,0,1,(B0(p2,M2Stau(gI1),M2Stau(gI2))*(ZTau(gI1,0)*(vu*(3*sqr(g1) - 5*sqr(g2))*ZTau(gI2,0) + 10*mu*sqrt2*ytau*ZTau(gI2,1)) + 2*ZTau(gI1,1)*(5*mu*sqrt2*ytau*ZTau(gI2,0) - 3*vu*sqr(g1)*ZTau(gI2,1)))*(ZTau(gI1,0)*(vd*(3*sqr(g1) - 5*sqr(g2) + 20*sqr(ytau))*ZTau(gI2,0) + 10*Atau*sqrt2*ytau*ZTau(gI2,1)) + 2*ZTau(gI1,1)*(5*Atau*sqrt2*ytau*ZTau(gI2,0) + vd*(-3*sqr(g1) + 10*sqr(ytau))*ZTau(gI2,1))))/400.));
   se12 += SUM(gI1,0,1,SUM(gI2,0,1,(3*vd*vu*B0(p2,M2Su(gI1),M2Su(gI2))*sqr((sqr(g1) - 5*sqr(g2))*ZU(gI1,0)*ZU(gI2,0) - 4*sqr(g1)*ZU(gI1,1)*ZU(gI2,1)))/400.));
   se12 += SUM(gI1,0,1,SUM(gI2,0,1,-(B0(p2,M2hh(gI1),M2hh(gI2))*sqr(3*sqr(g1) + 5*sqr(g2))*(ZH(gI1,0)*(vu*ZH(gI2,0) + vd*ZH(gI2,1)) + ZH(gI1,1)*(vd*ZH(gI2,0) - 3*vu*ZH(gI2,1)))*(ZH(gI1,1)*(vu*ZH(gI2,0) + vd*ZH(gI2,1)) + ZH(gI1,0)*(-3*vd*ZH(gI2,0) + vu*ZH(gI2,1))))/800.) + (A0(M2hh(gI1))*(3*sqr(g1) + 5*sqr(g2))*ZH(gI1,0)*ZH(gI1,1))/20.);
   se12 += SUM(gI1,0,1,SUM(gI2,0,1,-(B0(p2,M2Hpm(gI1),M2Hpm(gI2))*(ZP(gI1,0)*(vd*(3*sqr(g1) + 5*sqr(g2))*ZP(gI2,0) + 5*vu*sqr(g2)*ZP(gI2,1)) + ZP(gI1,1)*(5*vu*sqr(g2)*ZP(gI2,0) + vd*(-3*sqr(g1) + 5*sqr(g2))*ZP(gI2,1)))*(ZP(gI1,0)*((-3*vu*sqr(g1) + 5*vu*sqr(g2))*ZP(gI2,0) + 5*vd*sqr(g2)*ZP(gI2,1)) + ZP(gI1,1)*(5*vd*sqr(g2)*ZP(gI2,0) + vu*(3*sqr(g1) + 5*sqr(g2))*ZP(gI2,1))))/400.) - (A0(M2Hpm(gI1))*sqr(g2)*ZP(gI1,0)*ZP(gI1,1))/2.);
   se12 += SUM(gI1,0,3,SUM(gI2,0,3,-((A0(sqr(MChi(gI1))) + A0(sqr(MChi(gI2))) + B0(p2,sqr(MChi(gI1)),sqr(MChi(gI2)))*(-p2 + 2*MChi(gI1)*MChi(gI2) + sqr(MChi(gI1)) + sqr(MChi(gI2))))*(ZN(gI1,2)*(g1*sqrt15*ZN(gI2,0) - 5*g2*ZN(gI2,1)) + (g1*sqrt15*ZN(gI1,0) - 5*g2*ZN(gI1,1))*ZN(gI2,2))*(ZN(gI1,3)*(g1*sqrt15*ZN(gI2,0) - 5*g2*ZN(gI2,1)) + (g1*sqrt15*ZN(gI1,0) - 5*g2*ZN(gI1,1))*ZN(gI2,3)))/100.));
   se12 += SUM(gI2,0,1,-((2*A0(M2VZ) - A0(M2Ah(gI2)) + B0(p2,M2VZ,M2Ah(gI2))*(-M2VZ + 2*p2 + 2*M2Ah(gI2)))*(3*sqr(g1) + 5*sqr(g2))*ZA(gI2,0)*ZA(gI2,1))/20.);
   se12 += SUM(gI2,0,1,-((2*A0(M2VWm) - A0(M2Hpm(gI2)) + B0(p2,M2VWm,M2Hpm(gI2))*(-M2VWm + 2*p2 + 2*M2Hpm(gI2)))*sqr(g2)*ZP(gI2,0)*ZP(gI2,1))/2.);

   se22 += -(A0(M2VWm)*sqr(g2))/2.;
   se22 += (A0(M2VZ)*(-3*sqr(g1) - 5*sqr(g2)))/20.;
   se22 += -(B0(p2,M2SveL,M2SveL)*sqr(vu*(3*sqr(g1) + 5*sqr(g2))))/400.;
   se22 += -(B0(p2,M2SvmL,M2SvmL)*sqr(vu*(3*sqr(g1) + 5*sqr(g2))))/400.;
   se22 += -(B0(p2,M2SvtL,M2SvtL)*sqr(vu*(3*sqr(g1) + 5*sqr(g2))))/400.;
   se22 += (-7*B0(p2,M2VWm,M2VWm)*pow4(g2)*sqr(vu))/8.;
   se22 += (-7*B0(p2,M2VZ,M2VZ)*sqr(vu*(3*sqr(g1) + 5*sqr(g2))))/400.;
   se22 += 3*B0(p2,sqr(MFt),sqr(MFt))*(-p2 + 4*sqr(MFt))*sqr(yt);
   se22 += SUM(gI1,0,1,SUM(gI2,0,1,-(B0(p2,M2Ah(gI1),M2Ah(gI2))*sqr(vu*(3*sqr(g1) + 5*sqr(g2))*(ZA(gI1,0)*ZA(gI2,0) - ZA(gI1,1)*ZA(gI2,1))))/800.));
   se22 += SUM(gI1,0,1,SUM(gI2,0,1,(-3*B0(p2,M2Sc(gI1),M2Sc(gI2))*sqr(vu*((sqr(g1) - 5*sqr(g2))*ZC(gI1,0)*ZC(gI2,0) - 4*sqr(g1)*ZC(gI1,1)*ZC(gI2,1))))/400.));
   se22 += SUM(gI1,0,1,SUM(gI2,0,1,(-3*B0(p2,M2Sd(gI1),M2Sd(gI2))*sqr(vu*((sqr(g1) + 5*sqr(g2))*ZD(gI1,0)*ZD(gI2,0) + 2*sqr(g1)*ZD(gI1,1)*ZD(gI2,1))))/400.));
   se22 += SUM(gI1,0,1,SUM(gI2,0,1,-(B0(p2,M2Se(gI1),M2Se(gI2))*sqr(vu*((3*sqr(g1) - 5*sqr(g2))*ZE(gI1,0)*ZE(gI2,0) - 6*sqr(g1)*ZE(gI1,1)*ZE(gI2,1))))/400.));
   se22 += SUM(gI1,0,1,SUM(gI2,0,1,-(B0(p2,M2Sm(gI1),M2Sm(gI2))*sqr(vu*((3*sqr(g1) - 5*sqr(g2))*ZM(gI1,0)*ZM(gI2,0) - 6*sqr(g1)*ZM(gI1,1)*ZM(gI2,1))))/400.));
   se22 += SUM(gI1,0,1,SUM(gI2,0,1,(-3*B0(p2,M2Ss(gI1),M2Ss(gI2))*sqr(vu*((sqr(g1) + 5*sqr(g2))*ZS(gI1,0)*ZS(gI2,0) + 2*sqr(g1)*ZS(gI1,1)*ZS(gI2,1))))/400.));
   se22 += SUM(gI1,0,1,SUM(gI2,0,1,(-3*B0(p2,M2Su(gI1),M2Su(gI2))*sqr(vu*((sqr(g1) - 5*sqr(g2))*ZU(gI1,0)*ZU(gI2,0) - 4*sqr(g1)*ZU(gI1,1)*ZU(gI2,1))))/400.));
   se22 += SUM(gI1,0,1,SUM(gI2,0,1,(sqr(g2)*(A0(sqr(MCha(gI1))) + A0(sqr(MCha(gI2))) + B0(p2,sqr(MCha(gI1)),sqr(MCha(gI2)))*(-p2 + sqr(MCha(gI1)) + sqr(MCha(gI2))))*(sqr(UM(gI2,0)*UP(gI1,1)) + sqr(UM(gI1,0)*UP(gI2,1))))/2.) + MCha(gI1)*(SUM(gI2,0,1,2*B0(p2,sqr(MCha(gI1)),sqr(MCha(gI2)))*MCha(gI2)*sqr(g2)*UM(gI1,0)*UM(gI2,0)*UP(gI1,1)*UP(gI2,1)) - (2*g2*sqrt2*A0(sqr(MCha(gI1)))*UM(gI1,0)*UP(gI1,1))/vu));
   se22 += SUM(gI1,0,1,(-3*(SUM(gI2,0,1,(vu*B0(p2,M2Sb(gI1),M2Sb(gI2))*sqr(ZB(gI1,0)*(vu*(sqr(g1) + 5*sqr(g2))*ZB(gI2,0) - 10*mu*sqrt2*yb*ZB(gI2,1)) + 2*ZB(gI1,1)*(-5*mu*sqrt2*yb*ZB(gI2,0) + vu*sqr(g1)*ZB(gI2,1))))/400.) + mu*sqrt2*yb*A0(M2Sb(gI1))*ZB(gI1,0)*ZB(gI1,1)))/vu);
   se22 += SUM(gI1,0,1,-(SUM(gI2,0,1,(vu*B0(p2,M2hh(gI1),M2hh(gI2))*sqr((3*sqr(g1) + 5*sqr(g2))*(ZH(gI1,0)*(vu*ZH(gI2,0) + vd*ZH(gI2,1)) + ZH(gI1,1)*(vd*ZH(gI2,0) - 3*vu*ZH(gI2,1)))))/40.) + vd*A0(M2hh(gI1))*(3*sqr(g1) + 5*sqr(g2))*ZH(gI1,0)*ZH(gI1,1))/(20.*vu));
   se22 += SUM(gI1,0,1,SUM(gI2,0,1,-(B0(p2,M2Hpm(gI1),M2Hpm(gI2))*sqr(ZP(gI1,0)*(vu*(3*sqr(g1) - 5*sqr(g2))*ZP(gI2,0) - 5*vd*sqr(g2)*ZP(gI2,1)) - ZP(gI1,1)*(5*vd*sqr(g2)*ZP(gI2,0) + vu*(3*sqr(g1) + 5*sqr(g2))*ZP(gI2,1))))/400.) + (vd*A0(M2Hpm(gI1))*sqr(g2)*ZP(gI1,0)*ZP(gI1,1))/(2.*vu));
   se22 += SUM(gI1,0,1,SUM(gI2,0,1,(-3*B0(p2,M2St(gI1),M2St(gI2))*sqr(ZT(gI1,0)*(vu*(sqr(g1) - 5*sqr(g2) + 20*sqr(yt))*ZT(gI2,0) + 10*At*sqrt2*yt*ZT(gI2,1)) + 2*ZT(gI1,1)*(5*At*sqrt2*yt*ZT(gI2,0) - 2*vu*(sqr(g1) - 5*sqr(yt))*ZT(gI2,1))))/400.) + (3*At*sqrt2*yt*A0(M2St(gI1))*ZT(gI1,0)*ZT(gI1,1))/vu);
   se22 += SUM(gI1,0,1,-((SUM(gI2,0,1,(vu*B0(p2,M2Stau(gI1),M2Stau(gI2))*sqr(ZTau(gI1,0)*(vu*(3*sqr(g1) - 5*sqr(g2))*ZTau(gI2,0) + 10*mu*sqrt2*ytau*ZTau(gI2,1)) + 2*ZTau(gI1,1)*(5*mu*sqrt2*ytau*ZTau(gI2,0) - 3*vu*sqr(g1)*ZTau(gI2,1))))/400.) + mu*sqrt2*ytau*A0(M2Stau(gI1))*ZTau(gI1,0)*ZTau(gI1,1))/vu));
   se22 += SUM(gI1,0,3,SUM(gI2,0,3,((A0(sqr(MChi(gI1))) + A0(sqr(MChi(gI2))) + B0(p2,sqr(MChi(gI1)),sqr(MChi(gI2)))*(-p2 + 2*MChi(gI1)*MChi(gI2) + sqr(MChi(gI1)) + sqr(MChi(gI2))))*sqr(ZN(gI1,3)*(g1*sqrt15*ZN(gI2,0) - 5*g2*ZN(gI2,1)) + (g1*sqrt15*ZN(gI1,0) - 5*g2*ZN(gI1,1))*ZN(gI2,3)))/100.) - (2*A0(sqr(MChi(gI1)))*MChi(gI1)*(g1*sqrt15*ZN(gI1,0) - 5*g2*ZN(gI1,1))*ZN(gI1,3))/(5.*vu));
   se22 += SUM(gI2,0,1,((2*A0(M2VZ) - A0(M2Ah(gI2)) + B0(p2,M2VZ,M2Ah(gI2))*(-M2VZ + 2*p2 + 2*M2Ah(gI2)))*(3*sqr(g1) + 5*sqr(g2))*sqr(ZA(gI2,1)))/20.);
   se22 += SUM(gI2,0,1,((2*A0(M2VWm) - A0(M2Hpm(gI2)) + B0(p2,M2VWm,M2Hpm(gI2))*(-M2VWm + 2*p2 + 2*M2Hpm(gI2)))*sqr(g2*ZP(gI2,1)))/2.);

   RM22 se(RM22::Zero());
   se << se11, se12, se12, se22;

   return se * oneLoop;
}

/**
 * Higgs 1-loop DR' contribution in the gaugeless limit (p = g1 = g2 =
 * 0).  Note, that p = 0 in the gaugeless limit, because in the MSSM p
 * = mh = 0 vanishes when g1 = g2 = 0.
 *
 * The function returns 1/(4Pi)^2 (-selfenergy + tadpole).
 *
 * @return 1-loop contribution for p = g1 = g2 = 0
 */
RM22 MSSM_mass_eigenstates::delta_mh2_1loop_gaugeless() const
{
   const auto yt     = pars.Yu(2,2);
   const auto yb     = orders.at(CouplingOrders::ONLY_AT_AS) == 1 ? 0. : pars.Yd(2,2);
   const auto ytau   = orders.at(CouplingOrders::ONLY_AT_AS) == 1 ? 0. : pars.Ye(2,2);
   const auto vu     = pars.vu;
   const auto vd     = pars.vd;
   const auto mu     = pars.mu;
   const auto At     = pars.Au(2,2);
   const auto Ab     = pars.Ad(2,2);
   const auto Atau   = pars.Ae(2,2);
   const auto MFt    = gaugeless.MFt;
   const auto MFb    = gaugeless.MFb;
   const auto MFtau  = gaugeless.MFtau;
   const auto M2St   = gaugeless.M2St;
   const auto M2Sb   = gaugeless.M2Sb;
   const auto M2Stau = gaugeless.M2Stau;
   const auto ZT     = gaugeless.ZT;
   const auto ZB     = gaugeless.ZB;
   const auto ZTau   = gaugeless.ZTau;

   double se11{0.}, se12{0.}, se22{0.};

   se11 += 12*B0(0,sqr(MFb),sqr(MFb))*sqr(MFb*yb);
   se11 += 4*B0(0,sqr(MFtau),sqr(MFtau))*sqr(MFtau*ytau);
   se11 += (3*Ab*sqrt2*yb*A0(M2Sb(0))*ZB(0,0)*ZB(0,1))/vd;
   se11 += -6*B0(0,M2Sb(0),M2Sb(0))*sqr(yb*(MFb*(sqr(ZB(0,0)) + sqr(ZB(0,1))) + Ab*ZB(0,0)*ZB(0,1)));
   se11 += (3*Ab*sqrt2*yb*A0(M2Sb(1))*ZB(1,0)*ZB(1,1))/vd;
   se11 += -3*B0(0,M2Sb(0),M2Sb(1))*sqr(yb*(2*MFb*ZB(0,0)*ZB(1,0) + Ab*ZB(0,1)*ZB(1,0) + Ab*ZB(0,0)*ZB(1,1) + 2*MFb*ZB(0,1)*ZB(1,1)));
   se11 += -6*B0(0,M2Sb(1),M2Sb(1))*sqr(yb*(MFb*(sqr(ZB(1,0)) + sqr(ZB(1,1))) + Ab*ZB(1,0)*ZB(1,1)));
   se11 += (-3*mu*sqrt2*yt*A0(M2St(0))*ZT(0,0)*ZT(0,1))/vd;
   se11 += -6*B0(0,M2St(0),M2St(0))*sqr(mu*yt*ZT(0,0)*ZT(0,1));
   se11 += (-3*mu*sqrt2*yt*A0(M2St(1))*ZT(1,0)*ZT(1,1))/vd;
   se11 += -6*B0(0,M2St(1),M2St(1))*sqr(mu*yt*ZT(1,0)*ZT(1,1));
   se11 += -3*B0(0,M2St(0),M2St(1))*sqr(mu*yt*(ZT(0,1)*ZT(1,0) + ZT(0,0)*ZT(1,1)));
   se11 += (Atau*sqrt2*ytau*A0(M2Stau(0))*ZTau(0,0)*ZTau(0,1))/vd;
   se11 += -2*B0(0,M2Stau(0),M2Stau(0))*sqr(ytau*(MFtau*(sqr(ZTau(0,0)) + sqr(ZTau(0,1))) + Atau*ZTau(0,0)*ZTau(0,1)));
   se11 += (Atau*sqrt2*ytau*A0(M2Stau(1))*ZTau(1,0)*ZTau(1,1))/vd;
   se11 += -(B0(0,M2Stau(0),M2Stau(1))*sqr(ytau*(2*MFtau*ZTau(0,0)*ZTau(1,0) + Atau*ZTau(0,1)*ZTau(1,0) + Atau*ZTau(0,0)*ZTau(1,1) + 2*MFtau*ZTau(0,1)*ZTau(1,1))));
   se11 += -2*B0(0,M2Stau(1),M2Stau(1))*sqr(ytau*(MFtau*(sqr(ZTau(1,0)) + sqr(ZTau(1,1))) + Atau*ZTau(1,0)*ZTau(1,1)));

   se12 += 6*mu*B0(0,M2Sb(0),M2Sb(0))*sqr(yb)*ZB(0,0)*ZB(0,1)*(MFb*(sqr(ZB(0,0)) + sqr(ZB(0,1))) + Ab*ZB(0,0)*ZB(0,1));
   se12 += 3*mu*B0(0,M2Sb(0),M2Sb(1))*sqr(yb)*(ZB(0,1)*ZB(1,0) + ZB(0,0)*ZB(1,1))*(2*MFb*ZB(0,0)*ZB(1,0) + Ab*ZB(0,1)*ZB(1,0) + Ab*ZB(0,0)*ZB(1,1) + 2*MFb*ZB(0,1)*ZB(1,1));
   se12 += 6*mu*B0(0,M2Sb(1),M2Sb(1))*sqr(yb)*ZB(1,0)*ZB(1,1)*(MFb*(sqr(ZB(1,0)) + sqr(ZB(1,1))) + Ab*ZB(1,0)*ZB(1,1));
   se12 += 6*mu*B0(0,M2St(0),M2St(0))*sqr(yt)*ZT(0,0)*ZT(0,1)*(MFt*(sqr(ZT(0,0)) + sqr(ZT(0,1))) + At*ZT(0,0)*ZT(0,1));
   se12 += 3*mu*B0(0,M2St(0),M2St(1))*sqr(yt)*(ZT(0,1)*ZT(1,0) + ZT(0,0)*ZT(1,1))*(2*MFt*ZT(0,0)*ZT(1,0) + At*ZT(0,1)*ZT(1,0) + At*ZT(0,0)*ZT(1,1) + 2*MFt*ZT(0,1)*ZT(1,1));
   se12 += 6*mu*B0(0,M2St(1),M2St(1))*sqr(yt)*ZT(1,0)*ZT(1,1)*(MFt*(sqr(ZT(1,0)) + sqr(ZT(1,1))) + At*ZT(1,0)*ZT(1,1));
   se12 += 2*mu*B0(0,M2Stau(0),M2Stau(0))*sqr(ytau)*ZTau(0,0)*ZTau(0,1)*(MFtau*(sqr(ZTau(0,0)) + sqr(ZTau(0,1))) + Atau*ZTau(0,0)*ZTau(0,1));
   se12 += mu*B0(0,M2Stau(0),M2Stau(1))*sqr(ytau)*(ZTau(0,1)*ZTau(1,0) + ZTau(0,0)*ZTau(1,1))*(2*MFtau*ZTau(0,0)*ZTau(1,0) + Atau*ZTau(0,1)*ZTau(1,0) + Atau*ZTau(0,0)*ZTau(1,1) + 2*MFtau*ZTau(0,1)*ZTau(1,1));
   se12 += 2*mu*B0(0,M2Stau(1),M2Stau(1))*sqr(ytau)*ZTau(1,0)*ZTau(1,1)*(MFtau*(sqr(ZTau(1,0)) + sqr(ZTau(1,1))) + Atau*ZTau(1,0)*ZTau(1,1));

   se22 += 12*B0(0,sqr(MFt),sqr(MFt))*sqr(MFt*yt);
   se22 += (-3*mu*sqrt2*yb*A0(M2Sb(0))*ZB(0,0)*ZB(0,1))/vu;
   se22 += -6*B0(0,M2Sb(0),M2Sb(0))*sqr(mu*yb*ZB(0,0)*ZB(0,1));
   se22 += (-3*mu*sqrt2*yb*A0(M2Sb(1))*ZB(1,0)*ZB(1,1))/vu;
   se22 += -6*B0(0,M2Sb(1),M2Sb(1))*sqr(mu*yb*ZB(1,0)*ZB(1,1));
   se22 += -3*B0(0,M2Sb(0),M2Sb(1))*sqr(mu*yb*(ZB(0,1)*ZB(1,0) + ZB(0,0)*ZB(1,1)));
   se22 += (3*At*sqrt2*yt*A0(M2St(0))*ZT(0,0)*ZT(0,1))/vu;
   se22 += -6*B0(0,M2St(0),M2St(0))*sqr(yt*(MFt*(sqr(ZT(0,0)) + sqr(ZT(0,1))) + At*ZT(0,0)*ZT(0,1)));
   se22 += (3*At*sqrt2*yt*A0(M2St(1))*ZT(1,0)*ZT(1,1))/vu;
   se22 += -3*B0(0,M2St(0),M2St(1))*sqr(yt*(2*MFt*ZT(0,0)*ZT(1,0) + At*ZT(0,1)*ZT(1,0) + At*ZT(0,0)*ZT(1,1) + 2*MFt*ZT(0,1)*ZT(1,1)));
   se22 += -6*B0(0,M2St(1),M2St(1))*sqr(yt*(MFt*(sqr(ZT(1,0)) + sqr(ZT(1,1))) + At*ZT(1,0)*ZT(1,1)));
   se22 += -((mu*sqrt2*ytau*A0(M2Stau(0))*ZTau(0,0)*ZTau(0,1))/vu);
   se22 += -2*B0(0,M2Stau(0),M2Stau(0))*sqr(mu*ytau*ZTau(0,0)*ZTau(0,1));
   se22 += -((mu*sqrt2*ytau*A0(M2Stau(1))*ZTau(1,0)*ZTau(1,1))/vu);
   se22 += -2*B0(0,M2Stau(1),M2Stau(1))*sqr(mu*ytau*ZTau(1,0)*ZTau(1,1));
   se22 += -(B0(0,M2Stau(0),M2Stau(1))*sqr(mu*ytau*(ZTau(0,1)*ZTau(1,0) + ZTau(0,0)*ZTau(1,1))));

   RM22 se;
   se << se11, se12, se12, se22;

   return se * oneLoop;
}

/**
 * Derivative of Higgs 1-loop DR' contribution w.r.t. p^2, in the
 * gaugeless limit (p = g1 = g2 = 0).  Note, that p = 0 in the
 * gaugeless limit, because in the MSSM p = mh = 0 vanishes when g1 =
 * g2 = 0.
 *
 * The function returns 1/(4Pi)^2 d/dp^2 (-selfenergy + tadpole).
 *
 * @return derivative of 1-loop contribution for p = g1 = g2 = 0
 */
RM22 MSSM_mass_eigenstates::delta_mh2_1loop_gaugeless_deriv() const
{
   const auto yt     = pars.Yu(2,2);
   const auto yb     = orders.at(CouplingOrders::ONLY_AT_AS) == 1 ? 0. : pars.Yd(2,2);
   const auto ytau   = orders.at(CouplingOrders::ONLY_AT_AS) == 1 ? 0. : pars.Ye(2,2);
   const auto mu     = pars.mu;
   const auto At     = pars.Au(2,2);
   const auto Ab     = pars.Ad(2,2);
   const auto Atau   = pars.Ae(2,2);
   const auto MFt    = gaugeless.MFt;
   const auto MFb    = gaugeless.MFb;
   const auto MFtau  = gaugeless.MFtau;
   const auto M2St   = gaugeless.M2St;
   const auto M2Sb   = gaugeless.M2Sb;
   const auto M2Stau = gaugeless.M2Stau;
   const auto ZT     = gaugeless.ZT;
   const auto ZB     = gaugeless.ZB;
   const auto ZTau   = gaugeless.ZTau;

   double se11{0.}, se12{0.}, se22{0.};

   se11 += -3*B0(0,sqr(MFb),sqr(MFb))*sqr(yb);
   se11 += -(B0(0,sqr(MFtau),sqr(MFtau))*sqr(ytau));
   se11 += 12*D1B0(sqr(MFb),sqr(MFb))*sqr(MFb*yb);
   se11 += 4*D1B0(sqr(MFtau),sqr(MFtau))*sqr(MFtau*ytau);
   se11 += -6*D1B0(M2Sb(0),M2Sb(0))*sqr(yb*(MFb*(sqr(ZB(0,0)) + sqr(ZB(0,1))) + Ab*ZB(0,0)*ZB(0,1)));
   se11 += -3*D1B0(M2Sb(0),M2Sb(1))*sqr(yb*(2*MFb*ZB(0,0)*ZB(1,0) + Ab*ZB(0,1)*ZB(1,0) + Ab*ZB(0,0)*ZB(1,1) + 2*MFb*ZB(0,1)*ZB(1,1)));
   se11 += -6*D1B0(M2Sb(1),M2Sb(1))*sqr(yb*(MFb*(sqr(ZB(1,0)) + sqr(ZB(1,1))) + Ab*ZB(1,0)*ZB(1,1)));
   se11 += -6*D1B0(M2St(0),M2St(0))*sqr(mu*yt*ZT(0,0)*ZT(0,1));
   se11 += -6*D1B0(M2St(1),M2St(1))*sqr(mu*yt*ZT(1,0)*ZT(1,1));
   se11 += -3*D1B0(M2St(0),M2St(1))*sqr(mu*yt*(ZT(0,1)*ZT(1,0) + ZT(0,0)*ZT(1,1)));
   se11 += -2*D1B0(M2Stau(0),M2Stau(0))*sqr(ytau*(MFtau*(sqr(ZTau(0,0)) + sqr(ZTau(0,1))) + Atau*ZTau(0,0)*ZTau(0,1)));
   se11 += -(D1B0(M2Stau(0),M2Stau(1))*sqr(ytau*(2*MFtau*ZTau(0,0)*ZTau(1,0) + Atau*ZTau(0,1)*ZTau(1,0) + Atau*ZTau(0,0)*ZTau(1,1) + 2*MFtau*ZTau(0,1)*ZTau(1,1))));
   se11 += -2*D1B0(M2Stau(1),M2Stau(1))*sqr(ytau*(MFtau*(sqr(ZTau(1,0)) + sqr(ZTau(1,1))) + Atau*ZTau(1,0)*ZTau(1,1)));

   se12 += 6*mu*D1B0(M2Sb(0),M2Sb(0))*sqr(yb)*ZB(0,0)*ZB(0,1)*(MFb*(sqr(ZB(0,0)) + sqr(ZB(0,1))) + Ab*ZB(0,0)*ZB(0,1));
   se12 += 3*mu*D1B0(M2Sb(0),M2Sb(1))*sqr(yb)*(ZB(0,1)*ZB(1,0) + ZB(0,0)*ZB(1,1))*(2*MFb*ZB(0,0)*ZB(1,0) + Ab*ZB(0,1)*ZB(1,0) + Ab*ZB(0,0)*ZB(1,1) + 2*MFb*ZB(0,1)*ZB(1,1));
   se12 += 6*mu*D1B0(M2Sb(1),M2Sb(1))*sqr(yb)*ZB(1,0)*ZB(1,1)*(MFb*(sqr(ZB(1,0)) + sqr(ZB(1,1))) + Ab*ZB(1,0)*ZB(1,1));
   se12 += 6*mu*D1B0(M2St(0),M2St(0))*sqr(yt)*ZT(0,0)*ZT(0,1)*(MFt*(sqr(ZT(0,0)) + sqr(ZT(0,1))) + At*ZT(0,0)*ZT(0,1));
   se12 += 3*mu*D1B0(M2St(0),M2St(1))*sqr(yt)*(ZT(0,1)*ZT(1,0) + ZT(0,0)*ZT(1,1))*(2*MFt*ZT(0,0)*ZT(1,0) + At*ZT(0,1)*ZT(1,0) + At*ZT(0,0)*ZT(1,1) + 2*MFt*ZT(0,1)*ZT(1,1));
   se12 += 6*mu*D1B0(M2St(1),M2St(1))*sqr(yt)*ZT(1,0)*ZT(1,1)*(MFt*(sqr(ZT(1,0)) + sqr(ZT(1,1))) + At*ZT(1,0)*ZT(1,1));
   se12 += 2*mu*D1B0(M2Stau(0),M2Stau(0))*sqr(ytau)*ZTau(0,0)*ZTau(0,1)*(MFtau*(sqr(ZTau(0,0)) + sqr(ZTau(0,1))) + Atau*ZTau(0,0)*ZTau(0,1));
   se12 += mu*D1B0(M2Stau(0),M2Stau(1))*sqr(ytau)*(ZTau(0,1)*ZTau(1,0) + ZTau(0,0)*ZTau(1,1))*(2*MFtau*ZTau(0,0)*ZTau(1,0) + Atau*ZTau(0,1)*ZTau(1,0) + Atau*ZTau(0,0)*ZTau(1,1) + 2*MFtau*ZTau(0,1)*ZTau(1,1));
   se12 += 2*mu*D1B0(M2Stau(1),M2Stau(1))*sqr(ytau)*ZTau(1,0)*ZTau(1,1)*(MFtau*(sqr(ZTau(1,0)) + sqr(ZTau(1,1))) + Atau*ZTau(1,0)*ZTau(1,1));

   se22 += -3*B0(0,sqr(MFt),sqr(MFt))*sqr(yt);
   se22 += 12*D1B0(sqr(MFt),sqr(MFt))*sqr(MFt*yt);
   se22 += -6*D1B0(M2Sb(0),M2Sb(0))*sqr(mu*yb*ZB(0,0)*ZB(0,1));
   se22 += -6*D1B0(M2Sb(1),M2Sb(1))*sqr(mu*yb*ZB(1,0)*ZB(1,1));
   se22 += -3*D1B0(M2Sb(0),M2Sb(1))*sqr(mu*yb*(ZB(0,1)*ZB(1,0) + ZB(0,0)*ZB(1,1)));
   se22 += -6*D1B0(M2St(0),M2St(0))*sqr(yt*(MFt*(sqr(ZT(0,0)) + sqr(ZT(0,1))) + At*ZT(0,0)*ZT(0,1)));
   se22 += -3*D1B0(M2St(0),M2St(1))*sqr(yt*(2*MFt*ZT(0,0)*ZT(1,0) + At*ZT(0,1)*ZT(1,0) + At*ZT(0,0)*ZT(1,1) + 2*MFt*ZT(0,1)*ZT(1,1)));
   se22 += -6*D1B0(M2St(1),M2St(1))*sqr(yt*(MFt*(sqr(ZT(1,0)) + sqr(ZT(1,1))) + At*ZT(1,0)*ZT(1,1)));
   se22 += -2*D1B0(M2Stau(0),M2Stau(0))*sqr(mu*ytau*ZTau(0,0)*ZTau(0,1));
   se22 += -2*D1B0(M2Stau(1),M2Stau(1))*sqr(mu*ytau*ZTau(1,0)*ZTau(1,1));
   se22 += -(D1B0(M2Stau(0),M2Stau(1))*sqr(mu*ytau*(ZTau(0,1)*ZTau(1,0) + ZTau(0,0)*ZTau(1,1))));

   RM22 se;
   se << se11, se12, se12, se22;

   return se * oneLoop;
}

/**
 * CP-even Higgs 2-loop DR' contribution in the gaugeless limit (p =
 * g1 = g2 = 0).
 *
 * @note The 2-loop contribution to the heavy CP-even Higgs mass
 * assumes that the tree-level CP-even Higgs mass matrix is expressed
 * in terms of the CP-odd Higgs pole mass, see [hep-ph/0105096]
 * Eqs.(22)-(23).  If the tree-level CP-even Higgs mass matrix is
 * expressed in terms of the running CP-odd Higgs mass, then a
 * corresponding 2-loop contribution to the CP-odd Higgs mass must be
 * included in this function.  See the implementation in
 * SOFTSUSY/FlexibleSUSY for an example.
 *
 * @return 2-loop contribution for p = g1 = g2 = 0
 */
RM22 MSSM_mass_eigenstates::delta_mh2_2loop() const
{
   using namespace himalaya::mh2l;

   const auto g3 = pars.g3;
   const auto mt2 = pow2(gaugeless.MFt);
   const auto mb2 = pow2(gaugeless.MFb);
   const auto mtau2 = pow2(gaugeless.MFtau);
   const auto mg = pars.MG;
   const auto mst12 = pow2(pars.MSt(0));
   const auto mst22 = pow2(pars.MSt(1));
   const auto msb12 = pow2(pars.MSb(0));
   const auto msb22 = pow2(pars.MSb(1));
   const auto mstau12 = pow2(pars.MStau(0));
   const auto mstau22 = pow2(pars.MStau(1));
   const auto msv2 = gaugeless.M2SvtL;
   const auto mA2 = pow2(pars.MA);
   const auto sxt = std::sin(pars.theta_t);
   const auto cxt = std::cos(pars.theta_t);
   const auto sxb = std::sin(pars.theta_b);
   const auto cxb = std::cos(pars.theta_b);
   const auto sxtau = std::sin(pars.theta_tau);
   const auto cxtau = std::cos(pars.theta_tau);
   const auto scale2 = pow2(pars.scale);
   const auto mu = -pars.mu;
   const auto tanb = pars.vu/pars.vd;
   const auto cotb = 1./tanb;
   const auto vev2 = pow2(pars.vu) + pow2(pars.vd);
   const auto include_heavy_higgs = 1;

   RM22 dmh(RM22::Zero());

   // 2-loop contribution from momentum iteration
   switch (mom_it) {
   case Momentum_iteration::pert:
      dmh += delta_mh2_2loop_mom_it_pert();
      break;
   case Momentum_iteration::num:
      dmh += delta_mh2_2loop_mom_it_num(mom_it_precision_goal, mom_it_max_iterations);
      break;
   default:
      break;
   }

   if (orders.at(CouplingOrders::G32YT4)) {
      dmh += delta_mh2_2loop_at_as(
         mt2, mg, mst12, mst22, sxt, cxt, scale2, mu, tanb, vev2, g3,
         include_heavy_higgs);
   }

   if (orders.at(CouplingOrders::G32YB4)) {
      dmh += delta_mh2_2loop_ab_as(
         mb2, mg, msb12, msb22, sxb, cxb, scale2, mu, cotb, vev2, g3,
         include_heavy_higgs);
   }

   if (orders.at(CouplingOrders::YT6)) {
      dmh += delta_mh2_2loop_at_at(
         mt2, mb2, mA2, mst12, mst22, msb12, msb22,
         sxt, cxt, sxb, cxb, scale2, mu, tanb, vev2,
         include_heavy_higgs, orders.at(CouplingOrders::ONLY_AT_AS));
   }

   if (orders.at(CouplingOrders::YTAU6)) {
      dmh += delta_mh2_2loop_atau_atau(
         mtau2, mA2, msv2, mstau12, mstau22,
         sxtau, cxtau, scale2, mu, tanb, vev2,
         include_heavy_higgs);
   }

   if (orders.at(CouplingOrders::YTAU2YB4) || orders.at(CouplingOrders::YTAU4YB2)) {
      dmh += delta_mh2_2loop_ab_atau(
         mtau2, mb2, mstau12, mstau22, msb12, msb22,
         sxtau, cxtau, sxb, cxb, scale2, mu, tanb, vev2,
         include_heavy_higgs);
   }

   return dmh;
}

/**
 * Returns Higgs 2-loop contributions from momentum iteration of
 * 1-loop self-energy in the gaugeless limit p^2 = g1 = g2 = 0.
 *
 * @return 2-loop contribution from momentum iteration
 */
RM22 MSSM_mass_eigenstates::delta_mh2_2loop_mom_it_pert() const
{
   // tree-level Higgs mass matrix in gaugeless limit
   const auto DMH_0L = get_mass_matrix_hh_gaugeless();
   // 1-loop Higgs mass matrix in gaugeless limit
   const auto DMH_1L = delta_mh2_1loop_gaugeless();

   const auto dmh2 = fs_diagonalize_hermitian_perturbatively(DMH_0L, DMH_1L);

   // 1-loop contribution to (squared) Higgs mass eigenvalues.
   const auto dmh2_1L_gl = std::get<1>(dmh2)(0);

   return delta_mh2_1loop_gaugeless_deriv()*dmh2_1L_gl;
}

/**
 * Returns Higgs 2-loop (and higher) contributions from momentum
 * iteration of 1-loop self-energy.
 *
 * @param precision_goal precision goal (fraction, between 0 and 1)
 * @param max_iterations maximum number of iterations
 *
 * @return 2-loop (and higher) contribution from momentum iteration
 */
RM22 MSSM_mass_eigenstates::delta_mh2_2loop_mom_it_num(
   double precision_goal, int max_iterations) const
{
   const auto DMH_0L = get_mass_matrix_hh();
   const auto mh2_0L = masses.M2hh(0);
   auto p2 = mh2_0L;
   int n = 0;
   bool has_converged = false;
   RM22 DMH(RM22::Zero());

   while (!has_converged && n < max_iterations) {
      DMH = DMH_0L + delta_mh2_1loop(p2);

      MSSM_spectrum::A2 M2hh;
      RM22 ZH;
      fs_diagonalize_hermitian(DMH, M2hh, ZH);

      has_converged = is_equal_rel(M2hh(0), p2, precision_goal);
      p2 = M2hh(0);
      n++;
   };

   if (!has_converged) {
      WARNING_MSG("Momentum iteration has not converged after "
                  << n << " iterations.");
   }

   return DMH - DMH_0L - delta_mh2_1loop(mh2_0L);
}

RM22 MSSM_mass_eigenstates::get_mass_matrix_hh() const
{
   return masses.get_mass_matrix_hh();
}

RM22 MSSM_mass_eigenstates::get_mass_matrix_hh_gaugeless() const
{
   return gaugeless.get_mass_matrix_hh();
}

void MSSM_mass_eigenstates::set_correction(CouplingOrders::CouplingOrders order, int flag)
{
   if (flag < 0 || flag > 1)
      INFO_MSG("You can only enable (1) or disable (0) corrections!");

   orders.at(order) = flag;
}

void MSSM_mass_eigenstates::set_mom_it(
   Momentum_iteration mi,
   double mom_it_precision_goal_,
   int mom_it_max_iterations_)
{
   mom_it = mi;
   mom_it_precision_goal = mom_it_precision_goal_;
   mom_it_max_iterations = mom_it_max_iterations_;
}

void MSSM_mass_eigenstates::set_diagonalization(Diagonalization diag)
{
   diagonalization = diag;
}

double MSSM_mass_eigenstates::A0(double m2) const
{
   return a0(m2, sqr(pars.scale));
}

double MSSM_mass_eigenstates::B0(double p2, double m12, double m22) const
{
   return b0(p2, m12, m22, sqr(pars.scale));
}

double MSSM_mass_eigenstates::D1B0(double m12, double m22) const
{
   return d1_b0(m12, m22);
}

std::ostream& operator<<(std::ostream& ostr, const MSSM_mass_eigenstates& me)
{
   ostr << "====================================\n"
           "MSSM_mass_eigenstates\n"
           "====================================\n";
   ostr << me.pars;
   ostr << "------------------------------------\n"
           "DR' masses and mixings\n"
           "------------------------------------\n"
        << me.masses;
   ostr << "------------------------------------\n"
           "DR' masses and mixings (g1 = g2 = 0)\n"
           "------------------------------------\n"
        << me.gaugeless;

   return ostr;
}

} // namespace mh2_fo
} // namespace himalaya
