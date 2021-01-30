// ====================================================================
// This file is part of Himalaya.
//
// Himalaya is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================

#include "himalaya/mh2_fo/MSSM_spectrum.hpp"
#include "himalaya/mh2_fo/Linalg.hpp"
#include "himalaya/misc/Constants.hpp"
#include "himalaya/misc/Numerics.hpp"
#include "himalaya/misc/Powers.hpp"

/**
 * @file MSSM_spectrum.cpp
 *
 * @brief Contains the implementation of the \a MSSM_spectrum.
 */

namespace himalaya {
namespace mh2_fo {

namespace {

template <typename T>
T constexpr sqr(T x) noexcept { return x*x; }

#define DEFINE_COMMUTATIVE_OPERATOR_COMPLEX_INT(op)                     \
   template <typename T>                                                \
   std::complex<T> operator op(const std::complex<T>& lhs, int rhs)     \
   {                                                                    \
      return lhs op static_cast<T>(rhs);                                \
   }                                                                    \
                                                                        \
   template <typename T>                                                \
   std::complex<T> operator op(int lhs, const std::complex<T>& rhs)     \
   {                                                                    \
      return static_cast<T>(lhs) op rhs;                                \
   }

DEFINE_COMMUTATIVE_OPERATOR_COMPLEX_INT(*)
DEFINE_COMMUTATIVE_OPERATOR_COMPLEX_INT(/)
DEFINE_COMMUTATIVE_OPERATOR_COMPLEX_INT(+)
DEFINE_COMMUTATIVE_OPERATOR_COMPLEX_INT(-)

/**
 * Converts the given vector of masses and the corresponding (complex)
 * mixing matrix to SLHA convention: Matrix rows with non-zero
 * imaginary parts are multiplied by i and the corresponding mass
 * eigenvalue is multiplied by -1.  As a result the mixing matrix will
 * be real and the mass eigenvalues might be positive or negative.  It
 * is assumed that these mixings result from diagonalizing a symmetric
 * fermion mass matrix in the convention of Haber and Kane,
 * Phys. Rept. 117 (1985) 75-263.  This conversion makes sense only if
 * the original symmetric mass matrix is real-valued.
 *
 * @param m vector of masses
 * @param z mixing matrix
 */
template<int N>
void convert_symmetric_fermion_mixings_to_slha(
   Eigen::Array<double, N, 1>& m,
   Eigen::Matrix<std::complex<double>, N, N>& z)
{
   for (int i = 0; i < N; i++) {
      // check if i'th row contains non-zero imaginary parts
      if (!is_zero(z.row(i).imag().cwiseAbs().maxCoeff())) {
         z.row(i) *= std::complex<double>(0.0,1.0);
         m(i) *= -1;
      }
   }
}

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

/**
 * Fills lower triangle of symmetric matrix from values in upper
 * triangle.
 *
 * @param m matrix
 */
template <typename Derived>
void symmetrize(Eigen::MatrixBase<Derived>& m)
{
   static_assert(Eigen::MatrixBase<Derived>::RowsAtCompileTime ==
                 Eigen::MatrixBase<Derived>::ColsAtCompileTime,
                 "Symmetrize is only defined for squared matrices");

   for (int i = 0; i < Eigen::MatrixBase<Derived>::RowsAtCompileTime; i++)
      for (int k = 0; k < i; k++)
         m(i,k) = m(k,i);
}

int sign(double x)
{
   return x < 0 ? -1 : 1;
}

double signed_sqrt(double x)
{
   return sign(x)*std::sqrt(std::abs(x));
}

} // anonymous namespace

MSSM_spectrum::MSSM_spectrum(const Parameters& pars_)
   : pars(pars_)
{
   calculate_spectrum();
}

void MSSM_spectrum::calculate_spectrum()
{
   calculate_MVZ();
   calculate_MVWm();
   calculate_MFt();
   calculate_MFb();
   calculate_MFtau();
   calculate_MSveL();
   calculate_MSvmL();
   calculate_MSvtL();
   calculate_MSu();
   calculate_MSd();
   calculate_MSc();
   calculate_MSs();
   calculate_MSt();
   calculate_MSb();
   calculate_MSe();
   calculate_MSm();
   calculate_MStau();
   calculate_Mhh();
   calculate_MAh();
   calculate_MHpm();
   calculate_MChi();
   calculate_MCha();
}

void MSSM_spectrum::calculate_MVWm()
{
   const auto g2 = pars.g2;
   const auto vu = pars.vu;
   const auto vd = pars.vd;

   M2VWm = 0.25*sqr(g2)*(sqr(vd) + sqr(vu));
}

void MSSM_spectrum::calculate_MVZ()
{
   const auto gY = pars.g1 * sqrt35;
   const auto g2 = pars.g2;
   const auto vu = pars.vu;
   const auto vd = pars.vd;

   M2VZ = 0.25*(sqr(gY) + sqr(g2))*(sqr(vd) + sqr(vu));
}

void MSSM_spectrum::calculate_MFt()
{
   MFt = pars.Yu(2,2) * pars.vu * inv_sqrt2;
}

void MSSM_spectrum::calculate_MFb()
{
   MFb = pars.Yd(2,2) * pars.vd * inv_sqrt2;
}

void MSSM_spectrum::calculate_MFtau()
{
   MFtau = pars.Ye(2,2) * pars.vd * inv_sqrt2;
}

void MSSM_spectrum::calculate_MSveL()
{
   const auto g1 = pars.g1;
   const auto g2 = pars.g2;
   const auto vu = pars.vu;
   const auto vd = pars.vd;
   const auto ml2 = pars.ml2(0,0);

   M2SveL = 0.125*(8*ml2 - 0.6*sqr(g1)*(-sqr(vd)
      + sqr(vu)) - sqr(g2)*(-sqr(vd) + sqr(vu)));
}

void MSSM_spectrum::calculate_MSvmL()
{
   const auto g1 = pars.g1;
   const auto g2 = pars.g2;
   const auto vu = pars.vu;
   const auto vd = pars.vd;
   const auto ml2 = pars.ml2(1,1);

   M2SvmL = 0.125*(8*ml2 - 0.6*sqr(g1)*(-sqr(vd)
      + sqr(vu)) - sqr(g2)*(-sqr(vd) + sqr(vu)));
}

void MSSM_spectrum::calculate_MSvtL()
{
   const auto g1 = pars.g1;
   const auto g2 = pars.g2;
   const auto vu = pars.vu;
   const auto vd = pars.vd;
   const auto ml2 = pars.ml2(2,2);

   M2SvtL = 0.125*(8*ml2 - 0.6*sqr(g1)*(-sqr(vd)
      + sqr(vu)) - sqr(g2)*(-sqr(vd) + sqr(vu)));
}

RM22 MSSM_spectrum::get_mass_matrix_Su() const
{
   const auto g1 = pars.g1;
   const auto g2 = pars.g2;
   const auto yu = 0.; // pars.Yu(0,0);
   const auto vu = pars.vu;
   const auto vd = pars.vd;
   const auto mu = pars.mu;
   const auto Tyu = pars.Au(0,0) * yu;
   const auto mq2 = pars.mq2(0,0);
   const auto mu2 = pars.mu2(0,0);

   RM22 mass_matrix_Su;

   mass_matrix_Su(0,0) = mq2 - 0.025*sqr(g1)*sqr(vd) + 0.125*sqr(g2)*sqr(
      vd) + 0.5*sqr(yu)*sqr(vu) + 0.025*sqr(g1)*sqr(vu) - 0.125*sqr(g2)
      *sqr(vu);
   mass_matrix_Su(0,1) = inv_sqrt2*vu*Tyu - inv_sqrt2*vd*yu*mu;
   mass_matrix_Su(1,0) = mass_matrix_Su(0,1);
   mass_matrix_Su(1,1) = mu2 + 0.1*sqr(g1)*sqr(vd) + 0.5*sqr(yu)*
      sqr(vu) - 0.1*sqr(g1)*sqr(vu);

   return mass_matrix_Su;
}

void MSSM_spectrum::calculate_MSu()
{
   const auto mass_matrix_Su = get_mass_matrix_Su();
   fs_diagonalize_hermitian(mass_matrix_Su, M2Su, ZU);
   normalize_to_interval(ZU);
}

RM22 MSSM_spectrum::get_mass_matrix_Sd() const
{
   const auto g1 = pars.g1;
   const auto g2 = pars.g2;
   const auto yd = 0.; // pars.Yd(0,0);
   const auto vu = pars.vu;
   const auto vd = pars.vd;
   const auto mu = pars.mu;
   const auto Tyd = pars.Ad(0,0) * yd;
   const auto mq2 = pars.mq2(0,0);
   const auto md2 = pars.md2(0,0);

   RM22 mass_matrix_Sd;

   mass_matrix_Sd(0,0) = mq2 + 0.5*sqr(yd)*sqr(vd) - 0.025*sqr(g1)
      *sqr(vd) - 0.125*sqr(g2)*sqr(vd) + 0.025*sqr(g1)*sqr(vu) + 0.125*sqr(g2)*
      sqr(vu);
   mass_matrix_Sd(0,1) = inv_sqrt2*vd*Tyd - inv_sqrt2*vu*yd*mu;
   mass_matrix_Sd(1,0) = mass_matrix_Sd(0,1);
   mass_matrix_Sd(1,1) = md2 + 0.5*sqr(yd)*sqr(vd) - 0.05*sqr(g1)*
      sqr(vd) + 0.05*sqr(g1)*sqr(vu);

   return mass_matrix_Sd;
}

void MSSM_spectrum::calculate_MSd()
{
   const auto mass_matrix_Sd = get_mass_matrix_Sd();
   fs_diagonalize_hermitian(mass_matrix_Sd, M2Sd, ZD);
   normalize_to_interval(ZD);
}

RM22 MSSM_spectrum::get_mass_matrix_Sc() const
{
   const auto g1 = pars.g1;
   const auto g2 = pars.g2;
   const auto yc = 0.; // pars.Yu(1,1);
   const auto vu = pars.vu;
   const auto vd = pars.vd;
   const auto mu = pars.mu;
   const auto Tyc = pars.Au(1,1) * yc;
   const auto mq2 = pars.mq2(1,1);
   const auto mu2 = pars.mu2(1,1);

   RM22 mass_matrix_Sc;

   mass_matrix_Sc(0,0) = mq2 - 0.025*sqr(g1)*sqr(vd) + 0.125*sqr(g2)*sqr(
      vd) + 0.5*sqr(yc)*sqr(vu) + 0.025*sqr(g1)*sqr(vu) - 0.125*sqr(g2)
      *sqr(vu);
   mass_matrix_Sc(0,1) = inv_sqrt2*vu*Tyc - inv_sqrt2*vd*yc*mu;
   mass_matrix_Sc(1,0) = mass_matrix_Sc(0,1);
   mass_matrix_Sc(1,1) = mu2 + 0.1*sqr(g1)*sqr(vd) + 0.5*sqr(yc)*
      sqr(vu) - 0.1*sqr(g1)*sqr(vu);

   return mass_matrix_Sc;
}

void MSSM_spectrum::calculate_MSc()
{
   const auto mass_matrix_Sc = get_mass_matrix_Sc();
   fs_diagonalize_hermitian(mass_matrix_Sc, M2Sc, ZC);
   normalize_to_interval(ZC);
}

RM22 MSSM_spectrum::get_mass_matrix_Ss() const
{
   const auto g1 = pars.g1;
   const auto g2 = pars.g2;
   const auto ys = 0.; // pars.Yd(1,1);
   const auto vu = pars.vu;
   const auto vd = pars.vd;
   const auto mu = pars.mu;
   const auto Tys = pars.Ad(1,1) * ys;
   const auto mq2 = pars.mq2(1,1);
   const auto md2 = pars.md2(1,1);

   RM22 mass_matrix_Ss;

   mass_matrix_Ss(0,0) = mq2 + 0.5*sqr(ys)*sqr(vd) - 0.025*sqr(g1)
      *sqr(vd) - 0.125*sqr(g2)*sqr(vd) + 0.025*sqr(g1)*sqr(vu) + 0.125*sqr(g2)*
      sqr(vu);
   mass_matrix_Ss(0,1) = inv_sqrt2*vd*Tys - inv_sqrt2*vu*ys*mu;
   mass_matrix_Ss(1,0) = mass_matrix_Ss(0,1);
   mass_matrix_Ss(1,1) = md2 + 0.5*sqr(ys)*sqr(vd) - 0.05*sqr(g1)*
      sqr(vd) + 0.05*sqr(g1)*sqr(vu);

   return mass_matrix_Ss;
}

void MSSM_spectrum::calculate_MSs()
{
   const auto mass_matrix_Ss = get_mass_matrix_Ss();
   fs_diagonalize_hermitian(mass_matrix_Ss, M2Ss, ZS);
   normalize_to_interval(ZS);
}

RM22 MSSM_spectrum::get_mass_matrix_St() const
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

void MSSM_spectrum::calculate_MSt()
{
   const auto mass_matrix_St = get_mass_matrix_St();
   fs_diagonalize_hermitian(mass_matrix_St, M2St, ZT);
   normalize_to_interval(ZT);
}

RM22 MSSM_spectrum::get_mass_matrix_Sb() const
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

void MSSM_spectrum::calculate_MSb()
{
   const auto mass_matrix_Sb = get_mass_matrix_Sb();
   fs_diagonalize_hermitian(mass_matrix_Sb, M2Sb, ZB);
   normalize_to_interval(ZB);
}

RM22 MSSM_spectrum::get_mass_matrix_Se() const
{
   const auto g1 = pars.g1;
   const auto g2 = pars.g2;
   const auto ye = 0.; // pars.Ye(0,0);
   const auto vu = pars.vu;
   const auto vd = pars.vd;
   const auto mu = pars.mu;
   const auto Tye = pars.Ae(0,0) * ye;
   const auto ml2 = pars.ml2(0,0);
   const auto me2 = pars.me2(0,0);

   RM22 mass_matrix_Se;

   mass_matrix_Se(0,0) = ml2 + 0.5*sqr(ye)*sqr(vd) + 0.075*sqr(
      g1)*sqr(vd) - 0.125*sqr(g2)*sqr(vd) - 0.075*sqr(g1)*sqr(vu) + 0.125*sqr(
      g2)*sqr(vu);
   mass_matrix_Se(0,1) = inv_sqrt2*vd*Tye - inv_sqrt2*vu*ye*mu;
   mass_matrix_Se(1,0) = mass_matrix_Se(0,1);
   mass_matrix_Se(1,1) = me2 + 0.5*sqr(ye)*sqr(vd) - 0.15*sqr(g1
      )*sqr(vd) + 0.15*sqr(g1)*sqr(vu);

   return mass_matrix_Se;
}

void MSSM_spectrum::calculate_MSe()
{
   const auto mass_matrix_Se = get_mass_matrix_Se();
   fs_diagonalize_hermitian(mass_matrix_Se, M2Se, ZE);
   normalize_to_interval(ZE);
}

RM22 MSSM_spectrum::get_mass_matrix_Sm() const
{
   const auto g1 = pars.g1;
   const auto g2 = pars.g2;
   const auto ym = 0.; // pars.Ye(1,1);
   const auto vu = pars.vu;
   const auto vd = pars.vd;
   const auto mu = pars.mu;
   const auto Tym = pars.Ae(1,1) * ym;
   const auto ml2 = pars.ml2(1,1);
   const auto me2 = pars.me2(1,1);

   RM22 mass_matrix_Sm;

   mass_matrix_Sm(0,0) = ml2 + 0.5*sqr(ym)*sqr(vd) + 0.075*sqr(
      g1)*sqr(vd) - 0.125*sqr(g2)*sqr(vd) - 0.075*sqr(g1)*sqr(vu) + 0.125*sqr(
      g2)*sqr(vu);
   mass_matrix_Sm(0,1) = inv_sqrt2*vd*Tym - inv_sqrt2*vu*ym*mu;
   mass_matrix_Sm(1,0) = mass_matrix_Sm(0,1);
   mass_matrix_Sm(1,1) = me2 + 0.5*sqr(ym)*sqr(vd) - 0.15*sqr(g1
      )*sqr(vd) + 0.15*sqr(g1)*sqr(vu);

   return mass_matrix_Sm;
}

void MSSM_spectrum::calculate_MSm()
{
   const auto mass_matrix_Sm = get_mass_matrix_Sm();
   fs_diagonalize_hermitian(mass_matrix_Sm, M2Sm, ZM);
   normalize_to_interval(ZM);
}

RM22 MSSM_spectrum::get_mass_matrix_Stau() const
{
   const auto g1 = pars.g1;
   const auto g2 = pars.g2;
   const auto ytau = pars.Ye(2,2);
   const auto vu = pars.vu;
   const auto vd = pars.vd;
   const auto mu = pars.mu;
   const auto Tytau = pars.Ae(2,2) * ytau;
   const auto ml2 = pars.ml2(2,2);
   const auto me2 = pars.me2(2,2);

   RM22 mass_matrix_Stau;

   mass_matrix_Stau(0,0) = ml2 + 0.5*sqr(ytau)*sqr(vd) + 0.075*sqr(
      g1)*sqr(vd) - 0.125*sqr(g2)*sqr(vd) - 0.075*sqr(g1)*sqr(vu) + 0.125*sqr(
      g2)*sqr(vu);
   mass_matrix_Stau(0,1) = inv_sqrt2*vd*Tytau - inv_sqrt2*vu*ytau*mu;
   mass_matrix_Stau(1,0) = mass_matrix_Stau(0,1);
   mass_matrix_Stau(1,1) = me2 + 0.5*sqr(ytau)*sqr(vd) - 0.15*sqr(g1
      )*sqr(vd) + 0.15*sqr(g1)*sqr(vu);

   return mass_matrix_Stau;
}

void MSSM_spectrum::calculate_MStau()
{
   const auto mass_matrix_Stau = get_mass_matrix_Stau();
   fs_diagonalize_hermitian(mass_matrix_Stau, M2Stau, ZTau);
   normalize_to_interval(ZTau);
}

RM22 MSSM_spectrum::get_mass_matrix_hh() const
{
   const auto g1 = pars.g1;
   const auto g2 = pars.g2;
   const auto vu = pars.vu;
   const auto vd = pars.vd;
   const auto mA2 = sqr(pars.MA);
   const auto v2 = sqr(vu) + sqr(vd);
   const auto Bmu = mA2*vu*vd/v2;

   RM22 mass_matrix_hh;

   mass_matrix_hh(0,0) = (Bmu*vu)/vd + ((3*sqr(g1) + 5*sqr(g2))*sqr(vd))/20.;
   mass_matrix_hh(0,1) = -Bmu - (vd*vu*(3*sqr(g1) + 5*sqr(g2)))/20.;
   mass_matrix_hh(1,0) = mass_matrix_hh(0,1);
   mass_matrix_hh(1,1) = (Bmu*vd)/vu + ((3*sqr(g1) + 5*sqr(g2))*sqr(vu))/20.;

   return mass_matrix_hh;
}

void MSSM_spectrum::calculate_Mhh()
{
   const auto mass_matrix_hh = get_mass_matrix_hh();
   fs_diagonalize_hermitian(mass_matrix_hh, M2hh, ZH);
   normalize_to_interval(ZH);
}

RM22 MSSM_spectrum::get_mass_matrix_Ah() const
{
   const auto g1 = pars.g1;
   const auto g2 = pars.g2;
   const auto vu = pars.vu;
   const auto vd = pars.vd;
   const auto mA2 = sqr(pars.MA);
   const auto v2 = sqr(vu) + sqr(vd);
   const auto Bmu = mA2*vu*vd/v2;
   const auto xi = 1.;

   RM22 mass_matrix_Ah;

   mass_matrix_Ah(0,0) = (20*Bmu*vu + pow3(vd)*xi*(3*sqr(g1) + 5*sqr(g2)))/(20.*vd);
   mass_matrix_Ah(0,1) = Bmu - (vd*vu*xi*(3*sqr(g1) + 5*sqr(g2)))/20.;
   mass_matrix_Ah(1,0) = mass_matrix_Ah(0,1);
   mass_matrix_Ah(1,1) = (20*Bmu*vd + pow3(vu)*xi*(3*sqr(g1) + 5*sqr(g2)))/(20.*vu);

   return mass_matrix_Ah;
}

void MSSM_spectrum::calculate_MAh()
{
   const auto mass_matrix_Ah = get_mass_matrix_Ah();
   fs_diagonalize_hermitian(mass_matrix_Ah, M2Ah, ZA);
   normalize_to_interval(ZA);
}

RM22 MSSM_spectrum::get_mass_matrix_Hpm() const
{
   const auto g2 = pars.g2;
   const auto vu = pars.vu;
   const auto vd = pars.vd;
   const auto mA2 = sqr(pars.MA);
   const auto v2 = sqr(vu) + sqr(vd);
   const auto Bmu = mA2*vu*vd/v2;
   const auto xi = 1.;

   RM22 mass_matrix_Hpm;

   mass_matrix_Hpm(0,0) = (4*Bmu*vu + vd*sqr(g2)*(xi*sqr(vd) + sqr(vu)))/(4.*vd);
   mass_matrix_Hpm(0,1) = Bmu - (vd*vu*(-1 + xi)*sqr(g2))/4.;
   mass_matrix_Hpm(1,0) = mass_matrix_Hpm(0,1);
   mass_matrix_Hpm(1,1) = (4*Bmu*vd + vu*sqr(g2)*(sqr(vd) + xi*sqr(vu)))/(4.*vu);

   return mass_matrix_Hpm;
}

void MSSM_spectrum::calculate_MHpm()
{
   const auto mass_matrix_Hpm = get_mass_matrix_Hpm();
   fs_diagonalize_hermitian(mass_matrix_Hpm, M2Hpm, ZP);
   normalize_to_interval(ZP);
}

RM44 MSSM_spectrum::get_mass_matrix_Chi() const
{
   const auto g1 = pars.g1;
   const auto g2 = pars.g2;
   const auto vu = pars.vu;
   const auto vd = pars.vd;
   const auto mu = pars.mu;
   const auto M1 = pars.M1;
   const auto M2 = pars.M2;

   RM44 mass_matrix_Chi;

   mass_matrix_Chi(0,0) = M1;
   mass_matrix_Chi(0,1) = 0;
   mass_matrix_Chi(0,2) = -0.3872983346207417*g1*vd;
   mass_matrix_Chi(0,3) = 0.3872983346207417*g1*vu;
   mass_matrix_Chi(1,1) = M2;
   mass_matrix_Chi(1,2) = 0.5*g2*vd;
   mass_matrix_Chi(1,3) = -0.5*g2*vu;
   mass_matrix_Chi(2,2) = 0;
   mass_matrix_Chi(2,3) = -mu;
   mass_matrix_Chi(3,3) = 0;

   symmetrize(mass_matrix_Chi);

   return mass_matrix_Chi;
}

void MSSM_spectrum::calculate_MChi()
{
   CM44 ZN_tmp(CM44::Zero());

   const auto mass_matrix_Chi = get_mass_matrix_Chi();
   fs_diagonalize_symmetric(mass_matrix_Chi, MChi, ZN_tmp);
   normalize_to_interval(ZN_tmp);

   // convert to SLHA convention to avoid imaginary parts
   convert_symmetric_fermion_mixings_to_slha(MChi, ZN_tmp);
   ZN = ZN_tmp.real();
}

RM22 MSSM_spectrum::get_mass_matrix_Cha() const
{
   const auto g2 = pars.g2;
   const auto vu = pars.vu;
   const auto vd = pars.vd;
   const auto mu = pars.mu;
   const auto M2 = pars.M2;

   RM22 mass_matrix_Cha;

   mass_matrix_Cha(0,0) = M2;
   mass_matrix_Cha(0,1) = inv_sqrt2*g2*vu;
   mass_matrix_Cha(1,0) = inv_sqrt2*g2*vd;
   mass_matrix_Cha(1,1) = mu;

   return mass_matrix_Cha;
}

void MSSM_spectrum::calculate_MCha()
{
   const auto mass_matrix_Cha = get_mass_matrix_Cha();
   fs_svd(mass_matrix_Cha, MCha, UM, UP);
}

std::ostream& operator<<(std::ostream& ostr, const MSSM_spectrum& spec)
{
   const auto ssqrt = [] (double x) { return signed_sqrt(x); };

   ostr << "MFb = " << spec.MFb << '\n';
   ostr << "MFt = " << spec.MFt << '\n';
   ostr << "MFtau = " << spec.MFtau << '\n';
   ostr << "MSveL = " << signed_sqrt(spec.M2SveL) << '\n';
   ostr << "MSvmL = " << signed_sqrt(spec.M2SvmL) << '\n';
   ostr << "MSvtL = " << signed_sqrt(spec.M2SvtL) << '\n';
   ostr << "MSd = " << spec.M2Sd.transpose().unaryExpr(ssqrt) << '\n';
   ostr << "MSu = " << spec.M2Su.transpose().unaryExpr(ssqrt) << '\n';
   ostr << "MSe = " << spec.M2Se.transpose().unaryExpr(ssqrt) << '\n';
   ostr << "MSm = " << spec.M2Sm.transpose().unaryExpr(ssqrt) << '\n';
   ostr << "MStau = " << spec.M2Stau.transpose().unaryExpr(ssqrt) << '\n';
   ostr << "MSs = " << spec.M2Ss.transpose().unaryExpr(ssqrt) << '\n';
   ostr << "MSc = " << spec.M2Sc.transpose().unaryExpr(ssqrt) << '\n';
   ostr << "MSb = " << spec.M2Sb.transpose().unaryExpr(ssqrt) << '\n';
   ostr << "MSt = " << spec.M2St.transpose().unaryExpr(ssqrt) << '\n';
   ostr << "Mhh = " << spec.M2hh.transpose().unaryExpr(ssqrt) << '\n';
   ostr << "MAh = " << spec.M2Ah.transpose().unaryExpr(ssqrt) << '\n';
   ostr << "MHpm = " << spec.M2Hpm.transpose().unaryExpr(ssqrt) << '\n';
   ostr << "MChi = " << spec.MChi.transpose() << '\n';
   ostr << "MCha = " << spec.MCha.transpose() << '\n';
   ostr << "MVWm = " << signed_sqrt(spec.M2VWm) << '\n';
   ostr << "MVZ = " << signed_sqrt(spec.M2VZ) << '\n';
   ostr << "ZD = " << spec.ZD << '\n';
   ostr << "ZU = " << spec.ZU << '\n';
   ostr << "ZE = " << spec.ZE << '\n';
   ostr << "ZM = " << spec.ZM << '\n';
   ostr << "ZTau = " << spec.ZTau << '\n';
   ostr << "ZS = " << spec.ZS << '\n';
   ostr << "ZC = " << spec.ZC << '\n';
   ostr << "ZB = " << spec.ZB << '\n';
   ostr << "ZT = " << spec.ZT << '\n';
   ostr << "ZH = " << spec.ZH << '\n';
   ostr << "ZA = " << spec.ZA << '\n';
   ostr << "ZP = " << spec.ZP << '\n';
   ostr << "ZN = " << spec.ZN << '\n';
   ostr << "UM = " << spec.UM << '\n';
   ostr << "UP = " << spec.UP << '\n';

   return ostr;
}

} // namespace mh2_fo
} // namespace himalaya
