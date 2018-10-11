// ====================================================================
// This file is part of Himalaya.
//
// Himalaya is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================

#include "MSSM_mass_eigenstates.hpp"
#include "linalg2.hpp"
#include "pv.hpp"
#include "Logger.hpp"

namespace himalaya {
namespace mh1l {

namespace {

const double sqrt2 = 1.414213562373095;
const double inv_sqrt2 = 0.7071067811865475; // 1/sqrt2
const double one_loop = 0.006332573977646111; // 1/(4Pi)^2

double sqr(double x) { return x*x; }

/**
 * Normalize each element of the given real matrix to be within the
 * interval [min, max].  Values < min are set to min.  Values > max
 * are set to max.
 *
 * @param m matrix
 * @param min minimum
 * @param max maximum
 */
template <int M, int N>
void normalize_to_interval(Eigen::Matrix<double,M,N>& m, double min = -1., double max = 1.)
{
   auto data = m.data();
   const auto size = m.size();

   for (int i = 0; i < size; i++) {
      if (data[i] < min)
         data[i] = min;
      else if (data[i] > max)
         data[i] = max;
   }
}

/**
 * Normalize each element of the given complex matrix to have a
 * magnitude within the interval [0, max].  If the magnitude of a
 * matrix element is > max, then the magnitude is set to max.  The
 * phase angles are not modified.
 *
 * @param m matrix
 * @param max_mag maximum magnitude
 */
template <int M, int N>
void normalize_to_interval(Eigen::Matrix<std::complex<double>,M,N>& m, double max_mag = 1.)
{
   auto data = m.data();
   const auto size = m.size();

   for (int i = 0; i < size; i++) {
      if (std::abs(data[i]) > max_mag)
         data[i] = std::polar(max_mag, std::arg(data[i]));
   }
}

} // anonymous namespace

MSSM_mass_eigenstates::MSSM_mass_eigenstates(const Parameters& pars_)
   : pars(pars_)
{
   calculate_parameters();
}

void MSSM_mass_eigenstates::calculate_parameters()
{
   calculate_MFt();
   calculate_MFb();
   calculate_MFtau();
   calculate_MSt();
   calculate_MSb();
   calculate_MStau();
}

void MSSM_mass_eigenstates::calculate_MFt()
{
   MFt = pars.Yu(2,2) * pars.vu * inv_sqrt2;
}

void MSSM_mass_eigenstates::calculate_MFb()
{
   MFb = pars.Yd(2,2) * pars.vd * inv_sqrt2;
}

void MSSM_mass_eigenstates::calculate_MFtau()
{
   MFtau = pars.Ye(2,2) * pars.vd * inv_sqrt2;
}

RM22 MSSM_mass_eigenstates::get_mass_matrix_St() const
{
   const auto g1 = pars.g1;
   const auto g2 = pars.g2;
   const auto yt = pars.Yu(2,2);
   const auto vu = pars.vu;
   const auto vd = pars.vd;
   const auto mu = pars.mu;
   const auto Tyt = pars.Au(2,2) * yt;
   const auto mq2 = pars.mq2(2,2);
   const auto mu2 = pars.mu2(2,2);

   RM22 mass_matrix_St;

   mass_matrix_St(0,0) = mq2 - 0.025*sqr(g1)*sqr(vd) + 0.125*sqr(g2)*sqr(
      vd) + 0.5*sqr(yt)*sqr(vu) + 0.025*sqr(g1)*sqr(vu) - 0.125*sqr(g2)
      *sqr(vu);
   mass_matrix_St(0,1) = inv_sqrt2*vu*Tyt - inv_sqrt2*vd*yt*mu;
   mass_matrix_St(1,0) = mass_matrix_St(0,1);
   mass_matrix_St(1,1) = mu2 + 0.1*sqr(g1)*sqr(vd) + 0.5*sqr(yt)*
      sqr(vu) - 0.1*sqr(g1)*sqr(vu);

   return mass_matrix_St;
}

void MSSM_mass_eigenstates::calculate_MSt()
{
   const auto mass_matrix_St = get_mass_matrix_St();
   flexiblesusy::fs_diagonalize_hermitian(mass_matrix_St, M2St, ZT);
   normalize_to_interval(ZT);
}

RM22 MSSM_mass_eigenstates::get_mass_matrix_Sb() const
{
   const auto g1 = pars.g1;
   const auto g2 = pars.g2;
   const auto yb = pars.Yd(2,2);
   const auto vu = pars.vu;
   const auto vd = pars.vd;
   const auto mu = pars.mu;
   const auto Tyb = pars.Ad(2,2) * yb;
   const auto mq2 = pars.mq2(2,2);
   const auto md2 = pars.md2(2,2);

   RM22 mass_matrix_Sb;

   mass_matrix_Sb(0,0) = mq2 + 0.5*sqr(yb)*sqr(vd) - 0.025*sqr(g1)
      *sqr(vd) - 0.125*sqr(g2)*sqr(vd) + 0.025*sqr(g1)*sqr(vu) + 0.125*sqr(g2)*
      sqr(vu);
   mass_matrix_Sb(0,1) = inv_sqrt2*vd*Tyb - inv_sqrt2*vu*yb*mu;
   mass_matrix_Sb(1,0) = mass_matrix_Sb(0,1);
   mass_matrix_Sb(1,1) = md2 + 0.5*sqr(yb)*sqr(vd) - 0.05*sqr(g1)*
      sqr(vd) + 0.05*sqr(g1)*sqr(vu);

   return mass_matrix_Sb;
}

void MSSM_mass_eigenstates::calculate_MSb()
{
   const auto mass_matrix_Sb = get_mass_matrix_Sb();
   flexiblesusy::fs_diagonalize_hermitian(mass_matrix_Sb, M2Sb, ZB);
   normalize_to_interval(ZB);
}

RM22 MSSM_mass_eigenstates::get_mass_matrix_Stau() const
{
   const auto g1 = pars.g1;
   const auto g2 = pars.g2;
   const auto ye = pars.Ye(2,2);
   const auto vu = pars.vu;
   const auto vd = pars.vd;
   const auto mu = pars.mu;
   const auto Tye = pars.Ae(2,2) * ye;
   const auto ml2 = pars.ml2(2,2);
   const auto me2 = pars.me2(2,2);

   RM22 mass_matrix_Stau;

   mass_matrix_Stau(0,0) = ml2 + 0.5*sqr(ye)*sqr(vd) + 0.075*sqr(
      g1)*sqr(vd) - 0.125*sqr(g2)*sqr(vd) - 0.075*sqr(g1)*sqr(vu) + 0.125*sqr(
      g2)*sqr(vu);
   mass_matrix_Stau(0,1) = inv_sqrt2*vd*Tye - inv_sqrt2*vu*ye*mu;
   mass_matrix_Stau(1,0) = mass_matrix_Stau(0,1);
   mass_matrix_Stau(1,1) = me2 + 0.5*sqr(ye)*sqr(vd) - 0.15*sqr(g1
      )*sqr(vd) + 0.15*sqr(g1)*sqr(vu);

   return mass_matrix_Stau;
}

void MSSM_mass_eigenstates::calculate_MStau()
{
   const auto mass_matrix_Stau = get_mass_matrix_Stau();
   flexiblesusy::fs_diagonalize_hermitian(mass_matrix_Stau, M2Stau, ZTau);
   normalize_to_interval(ZTau);
}

/**
 * Higgs 1-loop DR' contribution for arbitrary momentum.
 *
 * @param p2 squared momentum
 *
 * @return 1-loop contribution
 */
RM22 MSSM_mass_eigenstates::delta_mh2_1loop(double p2) const
{
   RM22 se(RM22::Zero());

   throw "self_energy_h_1loop is not implemented";

   return se;
}

/**
 * Higgs 1-loop DR' contribution for p = 0 in the gaugeless limit (g1
 * = g2 = 0).
 *
 * @return 1-loop contribution for p = g1 = g2 = 0
 */
RM22 MSSM_mass_eigenstates::delta_mh2_1loop_effpot_gaugeless() const
{
   const auto yt = pars.Yu(2,2);
   const auto yb = pars.Yd(2,2);
   const auto ytau = pars.Ye(2,2);
   const auto vu = pars.vu;
   const auto vd = pars.vd;
   const auto mu = pars.mu;
   const auto At = pars.Au(2,2);
   const auto Ab = pars.Ad(2,2);
   const auto Atau = pars.Ae(2,2);

   double se11{0.}, se12{0.}, se22{0.};

   se11 += 12*B0(0,sqr(MFb),sqr(MFb))*sqr(MFb)*sqr(yb);
   se11 += 4*B0(0,sqr(MFtau),sqr(MFtau))*sqr(MFtau)*sqr(ytau);
   se11 += -3*B0(0,M2Sb(0),M2Sb(0))*sqr(yb)*sqr(vd*yb*(sqr(ZB(0,0)) + sqr(ZB(0,1))) + Ab*sqrt2*ZB(0,0)*ZB(0,1));
   se11 += -3*B0(0,M2Sb(1),M2Sb(1))*sqr(yb)*sqr(vd*yb*(sqr(ZB(1,0)) + sqr(ZB(1,1))) + Ab*sqrt2*ZB(1,0)*ZB(1,1));
   se11 += (-3*B0(0,M2Sb(0),M2Sb(1))*sqr(yb)*sqr(Ab*sqrt2*(ZB(0,1)*ZB(1,0) + ZB(0,0)*ZB(1,1)) + 2*vd*yb*(ZB(0,0)*ZB(1,0) + ZB(0,1)*ZB(1,1))))/4.;
   se11 += (-3*B0(0,M2Sb(1),M2Sb(0))*sqr(yb)*sqr(Ab*sqrt2*(ZB(0,1)*ZB(1,0) + ZB(0,0)*ZB(1,1)) + 2*vd*yb*(ZB(0,0)*ZB(1,0) + ZB(0,1)*ZB(1,1))))/4.;
   se11 += -6*B0(0,M2St(0),M2St(0))*sqr(mu)*sqr(yt)*sqr(ZT(0,0))*sqr(ZT(0,1));
   se11 += -6*B0(0,M2St(1),M2St(1))*sqr(mu)*sqr(yt)*sqr(ZT(1,0))*sqr(ZT(1,1));
   se11 += (-3*B0(0,M2St(0),M2St(1))*sqr(mu)*sqr(yt)*sqr(ZT(0,1)*ZT(1,0) + ZT(0,0)*ZT(1,1)))/2.;
   se11 += (-3*B0(0,M2St(1),M2St(0))*sqr(mu)*sqr(yt)*sqr(ZT(0,1)*ZT(1,0) + ZT(0,0)*ZT(1,1)))/2.;
   se11 += -(B0(0,M2Stau(0),M2Stau(0))*sqr(ytau)*sqr(vd*ytau*(sqr(ZTau(0,0)) + sqr(ZTau(0,1))) + Atau*sqrt2*ZTau(0,0)*ZTau(0,1)));
   se11 += -(B0(0,M2Stau(1),M2Stau(1))*sqr(ytau)*sqr(vd*ytau*(sqr(ZTau(1,0)) + sqr(ZTau(1,1))) + Atau*sqrt2*ZTau(1,0)*ZTau(1,1)));
   se11 += -(B0(0,M2Stau(0),M2Stau(1))*sqr(ytau)*sqr(Atau*sqrt2*(ZTau(0,1)*ZTau(1,0) + ZTau(0,0)*ZTau(1,1)) + 2*vd*ytau*(ZTau(0,0)*ZTau(1,0) + ZTau(0,1)*ZTau(1,1))))/4.;
   se11 += -(B0(0,M2Stau(1),M2Stau(0))*sqr(ytau)*sqr(Atau*sqrt2*(ZTau(0,1)*ZTau(1,0) + ZTau(0,0)*ZTau(1,1)) + 2*vd*ytau*(ZTau(0,0)*ZTau(1,0) + ZTau(0,1)*ZTau(1,1))))/4.;
   se11 += (3*Ab*sqrt2*yb*A0(M2Sb(0))*ZB(0,0)*ZB(0,1))/vd;
   se11 += (3*Ab*sqrt2*yb*A0(M2Sb(1))*ZB(1,0)*ZB(1,1))/vd;
   se11 += (-3*mu*sqrt2*yt*A0(M2St(0))*ZT(0,0)*ZT(0,1))/vd;
   se11 += (-3*mu*sqrt2*yt*A0(M2St(1))*ZT(1,0)*ZT(1,1))/vd;
   se11 += (Atau*sqrt2*ytau*A0(M2Stau(0))*ZTau(0,0)*ZTau(0,1))/vd;
   se11 += (Atau*sqrt2*ytau*A0(M2Stau(1))*ZTau(1,0)*ZTau(1,1))/vd;

   se12 += 3*mu*B0(0,M2Sb(0),M2Sb(0))*sqr(yb)*ZB(0,0)*ZB(0,1)*(sqrt2*vd*yb*(sqr(ZB(0,0)) + sqr(ZB(0,1))) + 2*Ab*ZB(0,0)*ZB(0,1));
   se12 += 3*mu*B0(0,M2Sb(1),M2Sb(1))*sqr(yb)*ZB(1,0)*ZB(1,1)*(sqrt2*vd*yb*(sqr(ZB(1,0)) + sqr(ZB(1,1))) + 2*Ab*ZB(1,0)*ZB(1,1));
   se12 += (3*mu*B0(0,M2Sb(0),M2Sb(1))*sqr(yb)*(ZB(0,1)*ZB(1,0) + ZB(0,0)*ZB(1,1))*(Ab*(ZB(0,1)*ZB(1,0) + ZB(0,0)*ZB(1,1)) + sqrt2*vd*yb*(ZB(0,0)*ZB(1,0) + ZB(0,1)*ZB(1,1))))/2.;
   se12 += (3*mu*B0(0,M2Sb(1),M2Sb(0))*sqr(yb)*(ZB(0,1)*ZB(1,0) + ZB(0,0)*ZB(1,1))*(Ab*(ZB(0,1)*ZB(1,0) + ZB(0,0)*ZB(1,1)) + sqrt2*vd*yb*(ZB(0,0)*ZB(1,0) + ZB(0,1)*ZB(1,1))))/2.;
   se12 += 3*mu*B0(0,M2St(0),M2St(0))*sqr(yt)*ZT(0,0)*ZT(0,1)*(sqrt2*vu*yt*(sqr(ZT(0,0)) + sqr(ZT(0,1))) + 2*At*ZT(0,0)*ZT(0,1));
   se12 += 3*mu*B0(0,M2St(1),M2St(1))*sqr(yt)*ZT(1,0)*ZT(1,1)*(sqrt2*vu*yt*(sqr(ZT(1,0)) + sqr(ZT(1,1))) + 2*At*ZT(1,0)*ZT(1,1));
   se12 += (3*mu*B0(0,M2St(0),M2St(1))*sqr(yt)*(ZT(0,1)*ZT(1,0) + ZT(0,0)*ZT(1,1))*(At*(ZT(0,1)*ZT(1,0) + ZT(0,0)*ZT(1,1)) + sqrt2*vu*yt*(ZT(0,0)*ZT(1,0) + ZT(0,1)*ZT(1,1))))/2.;
   se12 += (3*mu*B0(0,M2St(1),M2St(0))*sqr(yt)*(ZT(0,1)*ZT(1,0) + ZT(0,0)*ZT(1,1))*(At*(ZT(0,1)*ZT(1,0) + ZT(0,0)*ZT(1,1)) + sqrt2*vu*yt*(ZT(0,0)*ZT(1,0) + ZT(0,1)*ZT(1,1))))/2.;
   se12 += mu*sqrt2*B0(0,M2Stau(0),M2Stau(0))*sqr(ytau)*ZTau(0,0)*ZTau(0,1)*(vd*ytau*(sqr(ZTau(0,0)) + sqr(ZTau(0,1))) + Atau*sqrt2*ZTau(0,0)*ZTau(0,1));
   se12 += mu*sqrt2*B0(0,M2Stau(1),M2Stau(1))*sqr(ytau)*ZTau(1,0)*ZTau(1,1)*(vd*ytau*(sqr(ZTau(1,0)) + sqr(ZTau(1,1))) + Atau*sqrt2*ZTau(1,0)*ZTau(1,1));
   se12 += (mu*B0(0,M2Stau(0),M2Stau(1))*sqr(ytau)*(ZTau(0,1)*ZTau(1,0) + ZTau(0,0)*ZTau(1,1))*(Atau*sqrt2*(ZTau(0,1)*ZTau(1,0) + ZTau(0,0)*ZTau(1,1)) + 2*vd*ytau*(ZTau(0,0)*ZTau(1,0) + ZTau(0,1)*ZTau(1,1))))/(2.*sqrt2);
   se12 += (mu*B0(0,M2Stau(1),M2Stau(0))*sqr(ytau)*(ZTau(0,1)*ZTau(1,0) + ZTau(0,0)*ZTau(1,1))*(Atau*sqrt2*(ZTau(0,1)*ZTau(1,0) + ZTau(0,0)*ZTau(1,1)) + 2*vd*ytau*(ZTau(0,0)*ZTau(1,0) + ZTau(0,1)*ZTau(1,1))))/(2.*sqrt2);

   se22 += 12*B0(0,sqr(MFt),sqr(MFt))*sqr(MFt)*sqr(yt);
   se22 += -6*B0(0,M2Sb(0),M2Sb(0))*sqr(mu)*sqr(yb)*sqr(ZB(0,0))*sqr(ZB(0,1));
   se22 += -6*B0(0,M2Sb(1),M2Sb(1))*sqr(mu)*sqr(yb)*sqr(ZB(1,0))*sqr(ZB(1,1));
   se22 += (-3*B0(0,M2Sb(0),M2Sb(1))*sqr(mu)*sqr(yb)*sqr(ZB(0,1)*ZB(1,0) + ZB(0,0)*ZB(1,1)))/2.;
   se22 += (-3*B0(0,M2Sb(1),M2Sb(0))*sqr(mu)*sqr(yb)*sqr(ZB(0,1)*ZB(1,0) + ZB(0,0)*ZB(1,1)))/2.;
   se22 += -3*B0(0,M2St(0),M2St(0))*sqr(yt)*sqr(vu*yt*(sqr(ZT(0,0)) + sqr(ZT(0,1))) + At*sqrt2*ZT(0,0)*ZT(0,1));
   se22 += -3*B0(0,M2St(1),M2St(1))*sqr(yt)*sqr(vu*yt*(sqr(ZT(1,0)) + sqr(ZT(1,1))) + At*sqrt2*ZT(1,0)*ZT(1,1));
   se22 += (-3*B0(0,M2St(0),M2St(1))*sqr(yt)*sqr(At*sqrt2*(ZT(0,1)*ZT(1,0) + ZT(0,0)*ZT(1,1)) + 2*vu*yt*(ZT(0,0)*ZT(1,0) + ZT(0,1)*ZT(1,1))))/4.;
   se22 += (-3*B0(0,M2St(1),M2St(0))*sqr(yt)*sqr(At*sqrt2*(ZT(0,1)*ZT(1,0) + ZT(0,0)*ZT(1,1)) + 2*vu*yt*(ZT(0,0)*ZT(1,0) + ZT(0,1)*ZT(1,1))))/4.;
   se22 += -2*B0(0,M2Stau(0),M2Stau(0))*sqr(mu)*sqr(ytau)*sqr(ZTau(0,0))*sqr(ZTau(0,1));
   se22 += -2*B0(0,M2Stau(1),M2Stau(1))*sqr(mu)*sqr(ytau)*sqr(ZTau(1,0))*sqr(ZTau(1,1));
   se22 += -(B0(0,M2Stau(0),M2Stau(1))*sqr(mu)*sqr(ytau)*sqr(ZTau(0,1)*ZTau(1,0) + ZTau(0,0)*ZTau(1,1)))/2.;
   se22 += -(B0(0,M2Stau(1),M2Stau(0))*sqr(mu)*sqr(ytau)*sqr(ZTau(0,1)*ZTau(1,0) + ZTau(0,0)*ZTau(1,1)))/2.;
   se22 += (-3*mu*sqrt2*yb*A0(M2Sb(0))*ZB(0,0)*ZB(0,1))/vu;
   se22 += (-3*mu*sqrt2*yb*A0(M2Sb(1))*ZB(1,0)*ZB(1,1))/vu;
   se22 += (3*At*sqrt2*yt*A0(M2St(0))*ZT(0,0)*ZT(0,1))/vu;
   se22 += (3*At*sqrt2*yt*A0(M2St(1))*ZT(1,0)*ZT(1,1))/vu;
   se22 += -((mu*sqrt2*ytau*A0(M2Stau(0))*ZTau(0,0)*ZTau(0,1))/vu);
   se22 += -((mu*sqrt2*ytau*A0(M2Stau(1))*ZTau(1,0)*ZTau(1,1))/vu);

   RM22 se;
   se << se11, se12, se12, se22;

   return se * one_loop;
}

double MSSM_mass_eigenstates::A0(double m2) const
{
   return a0(m2, sqr(pars.scale));
}

double MSSM_mass_eigenstates::B0(double p2, double m12, double m22) const
{
   return b0(p2, m12, m22, sqr(pars.scale));
}

} // namespace mh1l
} // namespace himalaya
