#define BOOST_TEST_MODULE test_lambda_normalization

#include <boost/test/unit_test.hpp>
#include "HierarchyCalculator.hpp"
#include "Mh2EFTCalculator.hpp"

namespace {

const double Pi = 3.141592653589793;

double pow2(double x) { return x*x; }
double pow4(double x) { return x*x*x*x; }
double pow6(double x) { return x*x*x*x*x*x; }

himalaya::Parameters make_point()
{
   const double MS = 100000.;
   const double xt = 2.;
   const double tb = 5.;
   const double beta = std::atan(tb);
   const double v = 245.;

   const double msq2 = MS/2;
   const double msu2 = MS*2;
   const double mg   = MS*3;

   const double Xt = xt * std::sqrt(std::sqrt(msq2 * msu2));

   himalaya::Parameters pars;
   pars.scale = MS;
   pars.mu    = MS;
   pars.g3    = 1.05733;
   pars.vu    = v*std::sin(beta);
   pars.vd    = v*std::cos(beta);
   pars.mq2   << pow2(msq2), 0, 0, 0, pow2(msq2), 0, 0, 0, pow2(msq2);
   pars.md2   << pow2(MS), 0, 0, 0, pow2(MS), 0, 0, 0, pow2(MS);
   pars.mu2   << pow2(msu2), 0, 0, 0, pow2(msu2), 0, 0, 0, pow2(msu2);
   pars.At    = Xt + pars.mu/tb;
   pars.Ab    = 0;
   pars.MG    = mg;
   pars.MW    = 74.597;
   pars.MZ    = 85.7704;
   pars.Mt    = 144.337;
   pars.Mb    = 2.37054;
   pars.MA    = MS;

   return pars;
}

// calculates y_t in the MSSM
double calc_yt(const himalaya::Parameters& pars)
{
   return std::sqrt(2.)*pars.Mt/pars.vu;
}

// calculates y_t in the SM
double calc_gt(const himalaya::Parameters& pars)
{
   const double v = std::sqrt(pow2(pars.vu) + pow2(pars.vd));
   return std::sqrt(2.)*pars.Mt/v;
}

// calculates a_t
double calc_at(const himalaya::Parameters& pars)
{
   const double yt = calc_yt(pars);
   const double tb = pars.vu/pars.vd;
   const double beta = std::atan(tb);

   return pow2(yt*std::sin(beta))/(4*Pi);
}

/// calculates lightest mass eigenvalue of given mass matrix
double calc_Mh2(const Eigen::Matrix2d& Mh_mat)
{
   const Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> es(Mh_mat);
   return es.eigenvalues()(0);
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
   const double at = calc_at(pars);
   const double tb = pars.vu/pars.vd;
   const double mt = pars.Mt;
   const double mQ32 = pars.mq2(2,2);
   const double mU32 = pars.mu2(2,2);
   const double Xt = pars.At - pars.mu/tb;
   const double MR2 = pow2(pars.scale);

   himalaya::mh2_eft::Mh2EFTCalculator mhc;

   return mhc.Mh2_EFT_1loop(at, mt, mQ32, mU32, Xt, MR2);
}

/// calculates Mh^2 in the EFT at 2-loop level
double calc_Mh2_EFT_2L(const himalaya::Parameters& pars)
{
   const double at = calc_at(pars);
   const double tb = pars.vu/pars.vd;
   const double mt = pars.Mt;
   const double mQ32 = pars.mq2(2,2);
   const double mU32 = pars.mu2(2,2);
   const double Xt = pars.At - pars.mu/tb;
   const double MR2 = pow2(pars.scale);
   const double g3 = pars.g3;
   const double m3 = pars.MG;

   himalaya::mh2_eft::Mh2EFTCalculator mhc;

   return mhc.Mh2_EFT_2loop(at, mt, mQ32, mU32, Xt, MR2, g3, m3);
}

/// calculates Mh^2 in the EFT at 3-loop level
double calc_Mh2_EFT_3L(const himalaya::Parameters& pars,
                       double zeta_lambda_3L)
{
   const double at = calc_at(pars);
   const double tb = pars.vu/pars.vd;
   const double mt = pars.Mt;
   const double mQ32 = pars.mq2(2,2);
   const double mU32 = pars.mu2(2,2);
   const double Xt = pars.At - pars.mu/tb;
   const double MR2 = pow2(pars.scale);
   const double g3 = pars.g3;
   const double m3 = pars.MG;
   const double msq2 = std::pow(
      pars.mq2(0,0) * pars.mq2(1,1) *
      pars.mu2(0,0) * pars.mu2(1,1) *
      pars.md2(0,0) * pars.md2(1,1), 1./6.);
   const double gt  = calc_gt(pars);
   const double v2  = pow2(pars.vu) + pow2(pars.vd);
   const double v   = std::sqrt(v2);   // ~ 245 GeV
   const double vDO = v/std::sqrt(2.); // ~ 173 GeV

   himalaya::mh2_eft::Mh2EFTCalculator mhc;

   const double DMh2_EFT_3L_logs =
      mhc.Mh2_EFT_3loop(at, mt, mQ32, mU32, Xt, MR2, g3, m3, msq2);

   const double as = pow2(g3)/(4*pow2(Pi));
   // prefactor from paper
   const double pref = pow4(mt)/pow2(4*Pi*vDO)*pow2(as);
   // prefactor from HSSUSY
   const double pref2 = 8*pow4(gt*g3)/pow6(4*Pi) * v2;
   const double DMh2_EFT_3L_const = pref * zeta_lambda_3L;

   BOOST_CHECK_CLOSE_FRACTION(pref, pref2, 1e-10);

   return DMh2_EFT_3L_logs + DMh2_EFT_3L_const;
}

} // anonymous namespace

BOOST_AUTO_TEST_CASE(test_lambda_normalization)
{
   using namespace himalaya;

   const auto pars = make_point();
   auto hc = HierarchyCalculator(pars);
   const auto ho = hc.calculateDMh3L(false);

   const double zeta_Himalaya = ho.getZetaHimalaya();
   const double zeta_EFT      = ho.getZetaEFT();

   const auto DMh_0L = ho.getDMh(0);
   const auto DMh_1L = ho.getDMh(1);
   const auto DMh_2L = ho.getDMh(2);
   const auto DMh_3L = ho.getDMh(3);

   const auto Mh2_0L = calc_Mh2(DMh_0L);
   const auto Mh2_1L = calc_Mh2(DMh_0L + DMh_1L);
   const auto Mh2_2L = calc_Mh2(DMh_0L + DMh_1L + DMh_2L);
   const auto Mh2_3L = calc_Mh2(DMh_0L + DMh_1L + DMh_2L + DMh_3L);

   BOOST_TEST_MESSAGE("Mh(0L) = " << std::sqrt(Mh2_0L));
   BOOST_TEST_MESSAGE("Mh(1L) = " << std::sqrt(Mh2_1L));
   BOOST_TEST_MESSAGE("Mh(2L) = " << std::sqrt(Mh2_2L));
   BOOST_TEST_MESSAGE("Mh(3L) = " << std::sqrt(Mh2_3L));

   BOOST_CHECK_CLOSE_FRACTION(Mh2_0L, calc_Mh2_EFT_0L(pars), 1e-5);

   BOOST_CHECK_CLOSE_FRACTION(Mh2_1L, calc_Mh2_EFT_0L(pars)
                                    + calc_Mh2_EFT_1L(pars), 1e-6);

   BOOST_CHECK_CLOSE_FRACTION(Mh2_2L, calc_Mh2_EFT_0L(pars)
                                    + calc_Mh2_EFT_1L(pars)
                                    + calc_Mh2_EFT_2L(pars), 1e-6);

   BOOST_CHECK_CLOSE_FRACTION(Mh2_3L, calc_Mh2_EFT_0L(pars)
                                    + calc_Mh2_EFT_1L(pars)
                                    + calc_Mh2_EFT_2L(pars)
                                    + calc_Mh2_EFT_3L(pars,zeta_Himalaya),
                              0.002);

   BOOST_CHECK_CLOSE_FRACTION(Mh2_3L, calc_Mh2_EFT_0L(pars)
                                    + calc_Mh2_EFT_1L(pars)
                                    + calc_Mh2_EFT_2L(pars)
                                    + calc_Mh2_EFT_3L(pars,zeta_EFT),
                              0.03);

}
