#include "doctest.h"
#include "himalaya/mh2_eft/Mh2EFTCalculator.hpp"
#include "himalaya/mh2_fo/MSSM_mass_eigenstates.hpp"
#include "himalaya/misc/CouplingOrders.hpp"
#include "himalaya/mh2_fo/Linalg.hpp"

#define CHECK_CLOSE(a,b,eps) CHECK((a) == doctest::Approx(b).epsilon(eps))

namespace {

double pow2(double x) { return x*x; }

/**
 * Calculates the 1st derivative of \f$f(x)\f$ up to order \a Order
 * using the central finite difference.  This function calls \f$f\f$
 * 2 * (Order + 1) times.
 *
 * @param f function
 * @param x point at which derivative is to be calculated
 * @param eps measure for step size \f$h\f$
 * @tparam Order order of accuracy (0, 1, 2, 3)
 *
 * @return derivative
 */
template <int Order, class F, class A>
auto derivative_central(const F& f, A x, A eps = std::numeric_limits<A>::epsilon())
   -> decltype(f(x))
{
   static_assert(Order <= 3, "1st central derivative with order > 3 not implemented");

   using return_type = decltype(f(x));

   // coefficients from Math. Comp. 51 (1988), 699-706, Table 1
   // DOI: http://dx.doi.org/10.1090/S0025-5718-1988-0935077-0
   static const std::vector<std::vector<double> > coeffs = {
      {0.5},
      {2./3., -1./12.},
      {3./4., -3./20., 1./60.},
      {4./5., -1./5., 4./105., -1./280.}
   };

   const A h = std::fabs(x) < eps ? eps : std::sqrt(eps) * x;
   return_type result = 0;

   for (int i = 0; i < Order + 1; i++) {
      const double coeff = coeffs[Order][i];
      const A step = (i + 1) * h;
      result += coeff * (f(x + step) - f(x - step));
   }

   return result / h;
}

himalaya::Parameters make_point(double MS = 20000, double xt = 1, double tb = 20)
{
   const double beta = std::atan(tb);
   const double v = 245.;
   const double msq = MS;
   const double msu = MS;
   const double msd = MS;
   const double mg  = MS;
   const double msl = MS;
   const double mse = MS;

   const double Xt = xt * std::sqrt(msq * msu);

   himalaya::Parameters pars;
   pars.scale = MS;
   pars.mu    = MS;
   pars.g1    = 0.448494;
   pars.g2    = 0.607902;
   pars.g3    = 1.05733;
   pars.vu    = v*std::sin(beta);
   pars.vd    = v*std::cos(beta);
   pars.mq2   << pow2(msq), 0, 0, 0, pow2(msq), 0, 0, 0, pow2(msq);
   pars.md2   << pow2(msd), 0, 0, 0, pow2(msd), 0, 0, 0, pow2(msd);
   pars.mu2   << pow2(msu), 0, 0, 0, pow2(msu), 0, 0, 0, pow2(msu);
   pars.ml2   << pow2(msl), 0, 0, 0, pow2(msl), 0, 0, 0, pow2(msl);
   pars.me2   << pow2(mse), 0, 0, 0, pow2(mse), 0, 0, 0, pow2(mse);
   pars.Au    << 0, 0, 0, 0, 0, 0, 0, 0, Xt + pars.mu/tb;
   pars.Ae    << 0, 0, 0, 0, 0, 0, 0, 0, 0;
   pars.Ad    << 0, 0, 0, 0, 0, 0, 0, 0, 0;
   pars.Yu    << 0, 0, 0, 0, 0, 0, 0, 0, 0.8;
   pars.Ye    << 0, 0, 0, 0, 0, 0, 0, 0, 0.1;
   pars.Yd    << 0, 0, 0, 0, 0, 0, 0, 0, 0.3;
   pars.M1    = MS;
   pars.M2    = MS;
   pars.MG    = mg;
   // pars.MW    = 74.597;
   // pars.MZ    = 85.7704;
   // pars.Mt    = 144.337;
   // pars.Mb    = 2.37054;
   // pars.Mtau  = 1.2;
   pars.MA    = MS;

   pars.validate(false);

   return pars;
}

himalaya::Parameters make_gaugeless(const himalaya::Parameters& pars)
{
   auto gl = pars;
   gl.g1 = 0;
   gl.g2 = 0;
   gl.MW = himalaya::NaN;
   gl.MZ = himalaya::NaN;
   gl.validate(false);
   return gl;
}

/// calculates Mh^2 in the EFT at tree level
double calc_Mh2_EFT_0L(const himalaya::Parameters& pars)
{
   const double beta = std::atan(pars.vu/pars.vd);
   const double cos_2beta = std::cos(2*beta);

   return pow2(pars.MZ * cos_2beta);
}

/// calculates Mh^2 in the EFT at 1-loop level
double calc_Mh2_EFT_1L(const himalaya::Parameters& pars)
{
   himalaya::mh2_eft::Mh2EFTCalculator mhc(pars);

   return mhc.getDeltaMh2EFT1Loop(1,1);
}

/// calculates Mh^2 in the EFT at 1-loop level (only alpha_t contribution)
double calc_Mh2_EFT_1L_at(const himalaya::Parameters& pars)
{
   auto p = make_gaugeless(pars);
   p.Mb = 0.;
   p.Mtau = 0.;

   himalaya::mh2_eft::Mh2EFTCalculator mhc(p);

   return mhc.getDeltaMh2EFT1Loop(1,1);
}

/// calculates Mh^2 in the EFT at 2-loop level
double calc_Mh2_EFT_2L(const himalaya::Parameters& pars)
{
   himalaya::mh2_eft::Mh2EFTCalculator mhc(pars);

   return mhc.getDeltaMh2EFT2Loop(1,1);
}

} // anonymous namespace

TEST_CASE("test_FO_1loop_gaugeless")
{
   using namespace himalaya::mh2_fo;

   const double eps = 1e-10;
   const auto p    = make_point();
   const auto p_gl = make_gaugeless(p);
   const MSSM_mass_eigenstates me(p);
   const MSSM_mass_eigenstates me_gl(p_gl); // gaugeless
   const auto DMh2_1 = me.delta_mh2_1loop_gaugeless();
   const auto DMh2_2 = me_gl.delta_mh2_1loop_gaugeless();
   const auto DMh2_3 = me_gl.delta_mh2_1loop(0);

   CHECK_CLOSE(DMh2_1(0,0), DMh2_2(0,0), eps);
   CHECK_CLOSE(DMh2_1(0,1), DMh2_2(0,1), eps);
   CHECK_CLOSE(DMh2_1(1,0), DMh2_2(1,0), eps);
   CHECK_CLOSE(DMh2_1(1,1), DMh2_2(1,1), eps);
   CHECK_CLOSE(DMh2_2(0,0), DMh2_3(0,0), 1e-7);
   CHECK_CLOSE(DMh2_2(0,1), DMh2_3(0,1), eps);
   CHECK_CLOSE(DMh2_2(1,0), DMh2_3(1,0), eps);
   CHECK_CLOSE(DMh2_2(1,1), DMh2_3(1,1), 1e-7);
}

TEST_CASE("test_EFT_vs_FO_gaugeless")
{
   using namespace himalaya::mh2_fo;
   using namespace himalaya::mh2_eft;

   const auto p = make_gaugeless(make_point());

   const auto DMh2_EFT_0L    = calc_Mh2_EFT_0L(p);
   const auto DMh2_EFT_1L    = calc_Mh2_EFT_1L(p);
   const auto DMh2_EFT_2L    = calc_Mh2_EFT_2L(p);
   const auto DMh2_EFT_1L_at = calc_Mh2_EFT_1L_at(p);
   const auto DMh2_EFT_2L_at = calc_Mh2_EFT_2L(p);

   MSSM_mass_eigenstates me_full(p, false);
   MSSM_mass_eigenstates me_at(p, true);

   const auto DMh2_FO       = me_full.calculate_Mh2();
   const auto DMh2_FO_at    = me_at.calculate_Mh2();

   const auto DMh2_FO_0L    = std::get<0>(DMh2_FO);
   const auto DMh2_FO_1L    = std::get<1>(DMh2_FO);
   const auto DMh2_FO_2L    = std::get<2>(DMh2_FO);
   const auto DMh2_FO_0L_at = std::get<0>(DMh2_FO_at);
   const auto DMh2_FO_1L_at = std::get<1>(DMh2_FO_at);
   const auto DMh2_FO_2L_at = std::get<2>(DMh2_FO_at);

   INFO("DMh2_FO_0L = " << DMh2_FO_0L);
   INFO("DMh2_FO_1L = " << DMh2_FO_1L);
   INFO("DMh2_FO_2L = " << DMh2_FO_2L);

   INFO("DMh2_FO_0L = " << DMh2_FO_0L_at << " O(alpha_t*...)");
   INFO("DMh2_FO_1L = " << DMh2_FO_1L_at << " O(alpha_t*...)");
   INFO("DMh2_FO_2L = " << DMh2_FO_2L_at << " O(alpha_t*...)");

   INFO("DMh2_EFT_0L  = " << DMh2_EFT_0L);
   INFO("DMh2_EFT_1L  = " << DMh2_EFT_1L);
   INFO("DMh2_EFT_2L  = " << DMh2_EFT_2L);

   INFO("DMh2_EFT_1L  = " << DMh2_EFT_1L_at << " O(alpha_t*...)");
   INFO("DMh2_EFT_2L  = " << DMh2_EFT_2L_at << " O(alpha_t*...)");

   CHECK_CLOSE(DMh2_EFT_0L   , DMh2_FO_0L, 1e-6);

   // this test may be enabled once 1L O(αb + ατ) expressions re-enabled
   // CHECK_CLOSE(DMh2_EFT_1L   , DMh2_FO_1L, 1e-5);

   CHECK_CLOSE(0.            , DMh2_FO_0L_at, 1e-7);
   CHECK_CLOSE(DMh2_EFT_1L_at, DMh2_FO_1L_at, 1e-5);
   CHECK_CLOSE(DMh2_EFT_2L_at, DMh2_FO_2L_at, 1e-5);
}

TEST_CASE("test_FO_1loop_derivative")
{
   using namespace himalaya::mh2_fo;

   const auto p = make_gaugeless(make_point());
   const auto p2 = 1e-5;
   const auto eps = 1e-1;

   const MSSM_mass_eigenstates me(p);
   const auto Mh2_1L_deriv = me.delta_mh2_1loop_gaugeless_deriv()(0);

   const auto deriv = [&me] (double p2) { return me.delta_mh2_1loop(p2)(0); };
   const auto Mh2_1L_deriv_num = derivative_central<3>(deriv, p2, eps);

   CHECK_CLOSE(Mh2_1L_deriv, Mh2_1L_deriv_num, 1e-3);
}

TEST_CASE("test_FO_2loop_momentum_iteration_pert")
{
   using namespace himalaya;
   using namespace himalaya::mh2_fo;
   using A2 = Eigen::Array<double,2,1>;

   const auto p = make_gaugeless(make_point());
   MSSM_mass_eigenstates me(p);

   // disable 2-loop corrections
   me.set_correction(CouplingOrders::G32YT4, 0);
   me.set_correction(CouplingOrders::G32YB4, 0);
   me.set_correction(CouplingOrders::YT6, 0);
   me.set_correction(CouplingOrders::YTAU6, 0);
   me.set_correction(CouplingOrders::YTAU2YB4, 0);
   me.set_correction(CouplingOrders::YTAU4YB2, 0);

   // calculates Mh^2 as a function of p^2
   const auto Mh2_1L_p2 = [&me] (double p2) {
      const auto mm_0L = me.get_mass_matrix_hh();
      const auto mm_1L = me.delta_mh2_1loop(p2);
      const RM22 mm = mm_0L + mm_1L;

      A2 M2hh;
      RM22 ZH;
      // diagonalize 1-loop mass matrix for given p^2
      fs_diagonalize_hermitian(mm, M2hh, ZH);

      return M2hh(0);
   };

   const auto Mh2_full    = me.calculate_Mh2();
   const auto Mh2_full_0L = std::get<0>(Mh2_full);
   const auto Mh2_full_1L = Mh2_full_0L + std::get<1>(Mh2_full);
   const auto Mh2_full_2L = Mh2_full_1L + std::get<2>(Mh2_full);

   const auto Mh2_1L = Mh2_full_1L;
   const auto Mh2_2L_mom_it = Mh2_1L_p2(Mh2_1L_p2(0));
   const auto Mh2_2L = Mh2_full_2L;

   INFO("Mh2_1L = " << Mh2_1L);
   INFO("Mh2_2L = " << Mh2_2L);
   INFO("Mh2_2L_mom_it = " << Mh2_2L_mom_it);

   // check that (analytic and numeric) momentum iteration goes into
   // the same direction and is ~ 1% close to each other
   CHECK(std::abs(Mh2_2L - Mh2_2L_mom_it) < 0.01*std::abs(Mh2_2L - Mh2_1L));
}

TEST_CASE("test_FO_2loop_momentum_iteration_num")
{
   using namespace himalaya;
   using namespace himalaya::mh2_fo;

   const auto p = make_gaugeless(make_point());
   MSSM_mass_eigenstates me(p);

   // disable 2-loop corrections
   me.set_correction(CouplingOrders::G32YT4, 0);
   me.set_correction(CouplingOrders::G32YB4, 0);
   me.set_correction(CouplingOrders::YT6, 0);
   me.set_correction(CouplingOrders::YTAU6, 0);
   me.set_correction(CouplingOrders::YTAU2YB4, 0);
   me.set_correction(CouplingOrders::YTAU4YB2, 0);

   // calculate perturbatively
   me.set_mom_it(Momentum_iteration::pert);

   const auto Mh2_pert    = me.calculate_Mh2();
   const auto Mh2_pert_0L = std::get<0>(Mh2_pert);
   const auto Mh2_pert_1L = std::get<1>(Mh2_pert);
   const auto Mh2_pert_2L = std::get<2>(Mh2_pert);

   INFO("Mh2_pert_0L = " << Mh2_pert_0L);
   INFO("Mh2_pert_1L = " << Mh2_pert_1L);
   INFO("Mh2_pert_2L = " << Mh2_pert_2L);

   // calculate numerically
   me.set_mom_it(Momentum_iteration::num);

   const auto Mh2_num    = me.calculate_Mh2();
   const auto Mh2_num_0L = std::get<0>(Mh2_num);
   const auto Mh2_num_1L = std::get<1>(Mh2_num);
   const auto Mh2_num_2L = std::get<2>(Mh2_num);

   INFO("Mh2_num_0L = " << Mh2_num_0L);
   INFO("Mh2_num_1L = " << Mh2_num_1L);
   INFO("Mh2_num_2L = " << Mh2_num_2L);

   // check that (analytic and numeric) momentum iteration goes into
   // the same direction and is ~ 10% close to each other
   CHECK(std::abs(Mh2_pert_2L - Mh2_num_2L) < 0.1*std::abs(Mh2_pert_2L));
}
