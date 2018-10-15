// ====================================================================
// This file is part of Himalaya.
//
// Himalaya is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================

#include "MSSM_mass_eigenstates.hpp"
#include "linalg2.hpp"
#include "pv.hpp"
#include "sum.hpp"
#include "Logger.hpp"
#include <cmath>

namespace himalaya {
namespace mh1l {

namespace {

const double sqrt2 = 1.414213562373095;
const double sqrt15 = 3.872983346207417; // sqrt(15)
const double sqrt35 = 0.7745966692414834; // sqrt(3/5)
const double inv_sqrt2 = 0.7071067811865475; // 1/sqrt2
const double one_loop = 0.006332573977646111; // 1/(4Pi)^2

template <typename T> T sqr(T x) { return x*x; }
double pow2(double x) { return x*x; }
double pow3(double x) { return x*x*x; }
double pow4(double x) { return x*x*x*x; }
double pow6(double x) { return x*x*x*x*x*x; }

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

template <typename T>
typename std::enable_if<!std::is_unsigned<T>::value, bool>::type
is_zero(T a, T prec = std::numeric_limits<T>::epsilon()) noexcept
{
   return std::abs(a) <= prec;
}

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

void MSSM_spectrum::calculate_spectrum(const Parameters& pars)
{
   calculate_MVZ(pars);
   calculate_MVWm(pars);
   calculate_MFt(pars);
   calculate_MFb(pars);
   calculate_MFtau(pars);
   calculate_MSveL(pars);
   calculate_MSvmL(pars);
   calculate_MSvtL(pars);
   calculate_MSu(pars);
   calculate_MSd(pars);
   calculate_MSc(pars);
   calculate_MSs(pars);
   calculate_MSt(pars);
   calculate_MSb(pars);
   calculate_MSe(pars);
   calculate_MSm(pars);
   calculate_MStau(pars);
   calculate_Mhh(pars);
   calculate_MAh(pars);
   calculate_MHpm(pars);
   calculate_MChi(pars);
   calculate_MCha(pars);
}

void MSSM_spectrum::calculate_MVWm(const Parameters& pars)
{
   const auto g2 = pars.g2;
   const auto vu = pars.vu;
   const auto vd = pars.vd;

   M2VWm = 0.25*sqr(g2)*(sqr(vd) + sqr(vu));
}

void MSSM_spectrum::calculate_MVZ(const Parameters& pars)
{
   const auto gY = pars.g1 * sqrt35;
   const auto g2 = pars.g2;
   const auto vu = pars.vu;
   const auto vd = pars.vd;

   M2VZ = 0.25*(sqr(gY) + sqr(g2))*(sqr(vd) + sqr(vu));
}

void MSSM_spectrum::calculate_MFt(const Parameters& pars)
{
   MFt = pars.Yu(2,2) * pars.vu * inv_sqrt2;
}

void MSSM_spectrum::calculate_MFb(const Parameters& pars)
{
   MFb = pars.Yd(2,2) * pars.vd * inv_sqrt2;
}

void MSSM_spectrum::calculate_MFtau(const Parameters& pars)
{
   MFtau = pars.Ye(2,2) * pars.vd * inv_sqrt2;
}

void MSSM_spectrum::calculate_MSveL(const Parameters& pars)
{
   const auto g1 = pars.g1;
   const auto g2 = pars.g2;
   const auto vu = pars.vu;
   const auto vd = pars.vd;
   const auto ml2 = pars.ml2(0,0);

   M2SveL = 0.125*(8*ml2 - 0.6*sqr(g1)*(-sqr(vd)
      + sqr(vu)) - sqr(g2)*(-sqr(vd) + sqr(vu)));
}

void MSSM_spectrum::calculate_MSvmL(const Parameters& pars)
{
   const auto g1 = pars.g1;
   const auto g2 = pars.g2;
   const auto vu = pars.vu;
   const auto vd = pars.vd;
   const auto ml2 = pars.ml2(1,1);

   M2SvmL = 0.125*(8*ml2 - 0.6*sqr(g1)*(-sqr(vd)
      + sqr(vu)) - sqr(g2)*(-sqr(vd) + sqr(vu)));
}

void MSSM_spectrum::calculate_MSvtL(const Parameters& pars)
{
   const auto g1 = pars.g1;
   const auto g2 = pars.g2;
   const auto vu = pars.vu;
   const auto vd = pars.vd;
   const auto ml2 = pars.ml2(2,2);

   M2SvtL = 0.125*(8*ml2 - 0.6*sqr(g1)*(-sqr(vd)
      + sqr(vu)) - sqr(g2)*(-sqr(vd) + sqr(vu)));
}

RM22 MSSM_spectrum::get_mass_matrix_Su(const Parameters& pars) const
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

void MSSM_spectrum::calculate_MSu(const Parameters& pars)
{
   const auto mass_matrix_Su = get_mass_matrix_Su(pars);
   flexiblesusy::fs_diagonalize_hermitian(mass_matrix_Su, M2Su, ZU);
   normalize_to_interval(ZU);
}

RM22 MSSM_spectrum::get_mass_matrix_Sd(const Parameters& pars) const
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

void MSSM_spectrum::calculate_MSd(const Parameters& pars)
{
   const auto mass_matrix_Sd = get_mass_matrix_Sd(pars);
   flexiblesusy::fs_diagonalize_hermitian(mass_matrix_Sd, M2Sd, ZD);
   normalize_to_interval(ZD);
}

RM22 MSSM_spectrum::get_mass_matrix_Sc(const Parameters& pars) const
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

void MSSM_spectrum::calculate_MSc(const Parameters& pars)
{
   const auto mass_matrix_Sc = get_mass_matrix_Sc(pars);
   flexiblesusy::fs_diagonalize_hermitian(mass_matrix_Sc, M2Sc, ZC);
   normalize_to_interval(ZC);
}

RM22 MSSM_spectrum::get_mass_matrix_Ss(const Parameters& pars) const
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

void MSSM_spectrum::calculate_MSs(const Parameters& pars)
{
   const auto mass_matrix_Ss = get_mass_matrix_Ss(pars);
   flexiblesusy::fs_diagonalize_hermitian(mass_matrix_Ss, M2Ss, ZS);
   normalize_to_interval(ZS);
}

RM22 MSSM_spectrum::get_mass_matrix_St(const Parameters& pars) const
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

void MSSM_spectrum::calculate_MSt(const Parameters& pars)
{
   const auto mass_matrix_St = get_mass_matrix_St(pars);
   flexiblesusy::fs_diagonalize_hermitian(mass_matrix_St, M2St, ZT);
   normalize_to_interval(ZT);
}

RM22 MSSM_spectrum::get_mass_matrix_Sb(const Parameters& pars) const
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

void MSSM_spectrum::calculate_MSb(const Parameters& pars)
{
   const auto mass_matrix_Sb = get_mass_matrix_Sb(pars);
   flexiblesusy::fs_diagonalize_hermitian(mass_matrix_Sb, M2Sb, ZB);
   normalize_to_interval(ZB);
}

RM22 MSSM_spectrum::get_mass_matrix_Se(const Parameters& pars) const
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

void MSSM_spectrum::calculate_MSe(const Parameters& pars)
{
   const auto mass_matrix_Se = get_mass_matrix_Se(pars);
   flexiblesusy::fs_diagonalize_hermitian(mass_matrix_Se, M2Se, ZE);
   normalize_to_interval(ZE);
}

RM22 MSSM_spectrum::get_mass_matrix_Sm(const Parameters& pars) const
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

void MSSM_spectrum::calculate_MSm(const Parameters& pars)
{
   const auto mass_matrix_Sm = get_mass_matrix_Sm(pars);
   flexiblesusy::fs_diagonalize_hermitian(mass_matrix_Sm, M2Sm, ZM);
   normalize_to_interval(ZM);
}

RM22 MSSM_spectrum::get_mass_matrix_Stau(const Parameters& pars) const
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

void MSSM_spectrum::calculate_MStau(const Parameters& pars)
{
   const auto mass_matrix_Stau = get_mass_matrix_Stau(pars);
   flexiblesusy::fs_diagonalize_hermitian(mass_matrix_Stau, M2Stau, ZTau);
   normalize_to_interval(ZTau);
}

RM22 MSSM_spectrum::get_mass_matrix_hh(const Parameters& pars) const
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

void MSSM_spectrum::calculate_Mhh(const Parameters& pars)
{
   const auto mass_matrix_hh = get_mass_matrix_hh(pars);
   flexiblesusy::fs_diagonalize_hermitian(mass_matrix_hh, M2hh, ZH);
   normalize_to_interval(ZH);
}

RM22 MSSM_spectrum::get_mass_matrix_Ah(const Parameters& pars) const
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

void MSSM_spectrum::calculate_MAh(const Parameters& pars)
{
   const auto mass_matrix_Ah = get_mass_matrix_Ah(pars);
   flexiblesusy::fs_diagonalize_hermitian(mass_matrix_Ah, M2Ah, ZA);
   normalize_to_interval(ZA);
}

RM22 MSSM_spectrum::get_mass_matrix_Hpm(const Parameters& pars) const
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

void MSSM_spectrum::calculate_MHpm(const Parameters& pars)
{
   const auto mass_matrix_Hpm = get_mass_matrix_Hpm(pars);
   flexiblesusy::fs_diagonalize_hermitian(mass_matrix_Hpm, M2Hpm, ZP);
   normalize_to_interval(ZP);
}

RM44 MSSM_spectrum::get_mass_matrix_Chi(const Parameters& pars) const
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

void MSSM_spectrum::calculate_MChi(const Parameters& pars)
{
   CM44 ZN_tmp(CM44::Zero());

   const auto mass_matrix_Chi = get_mass_matrix_Chi(pars);
   flexiblesusy::fs_diagonalize_symmetric(mass_matrix_Chi, MChi, ZN_tmp);
   normalize_to_interval(ZN_tmp);

   // convert to SLHA convention to avoid imaginary parts
   convert_symmetric_fermion_mixings_to_slha(MChi, ZN_tmp);
   ZN = ZN_tmp.real();
}

RM22 MSSM_spectrum::get_mass_matrix_Cha(const Parameters& pars) const
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

void MSSM_spectrum::calculate_MCha(const Parameters& pars)
{
   const auto mass_matrix_Cha = get_mass_matrix_Cha(pars);
   flexiblesusy::fs_svd(mass_matrix_Cha, MCha, UM, UP);
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

/* ************************************************************ */

MSSM_mass_eigenstates::MSSM_mass_eigenstates(const Parameters& pars_)
   : pars(pars_)
{
   calculate_parameters();
}

void MSSM_mass_eigenstates::calculate_parameters()
{
   masses.calculate_spectrum(pars);
   gaugeless.calculate_spectrum(make_gaugeless(pars));
}

Parameters MSSM_mass_eigenstates::make_gaugeless(const Parameters& pars) const
{
   auto gl = pars;
   gl.g1 = 0.;
   gl.g2 = 0.;
   return gl;
}

/**
 * Returns the Higgs pole masses at a given loop level.  The function
 * does not include implicit or explicit higher orders.
 *
 * @param loops number of loops
 *
 * @return Higgs pole masses
 */
V2 MSSM_mass_eigenstates::calculate_Mh2(int loops) const
{
   using std::sqrt;

   V2 mh2(V2::Zero());

   if (loops >= 0) {
      const auto m0 = masses.get_mass_matrix_hh(pars);
      const auto a11 = m0(0,0), a12 = m0(0,1), a22 = m0(1,1);
      const auto c1 = a11 + a22;
      const auto c2 = sqrt(sqr(a11) + 4*sqr(a12) - 2*a11*a22 + sqr(a22));

      mh2(0) += 0.5*(c1 - c2);
      mh2(1) += 0.5*(c1 + c2);

      if (loops > 0) {
         const auto p2 = mh2(0);
         const auto m1 = delta_mh2_1loop(p2);
         const auto b11 = m1(0,0), b12 = m1(0,1), b22 = m1(1,1);
         const auto c3 = b11 + b22;
         const auto c4 = (a11*b11 - a22*b11 + 4*a12*b12 - a11*b22 + a22*b22)/c2;

         mh2(0) += 0.5*(c3 - c4);
         mh2(1) += 0.5*(c3 + c4);

         if (loops > 1) {
            // 2-loop contribution from momentum iteration
            const RM22 m2_mom = delta_mh2_1loop_gaugeless() * delta_mh2_1loop_gaugeless_deriv();
            const auto d11 = m2_mom(0,0), d12 = m2_mom(0,1), d22 = m2_mom(1,1);
            const auto c5 = 0.5*c2*(
               - pow2(c4)/pow6(c2)
               + (pow2(b11) + 4*pow2(b12) - 2*b11*b22 + pow2(b22) + 2*a11*d11
                  - 2*a22*d11 + 8*a12*d12 - 2*a11*d22 + 2*a22*d22)/pow4(c2));

            mh2(0) += 0.5*(d11 + d22 - c5);
            mh2(1) += 0.5*(d11 + d22 + c5);
         }
      }
   }

   return mh2;
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

   se11 += (A0(M2VZ)*(-3*sqr(g1) - 5*sqr(g2)))/20.;
   se11 += -(A0(M2VWm)*sqr(g2))/2.;
   se11 += (-7*B0(p2,M2VWm,M2VWm)*pow4(g2)*sqr(vd))/8.;
   se11 += -3*B0(p2,sqr(MFb),sqr(MFb))*(p2 - 4*sqr(MFb))*sqr(yb);
   se11 += -(B0(p2,sqr(MFtau),sqr(MFtau))*(p2 - 4*sqr(MFtau))*sqr(ytau));
   se11 += -(B0(p2,M2SveL,M2SveL)*sqr(vd*(3*sqr(g1) + 5*sqr(g2))))/400.;
   se11 += -(B0(p2,M2SvmL,M2SvmL)*sqr(vd*(3*sqr(g1) + 5*sqr(g2))))/400.;
   se11 += -(B0(p2,M2SvtL,M2SvtL)*sqr(vd*(3*sqr(g1) + 5*sqr(g2))))/400.;
   se11 += (-7*B0(p2,M2VZ,M2VZ)*sqr(vd*(3*sqr(g1) + 5*sqr(g2))))/400.;
   se11 += SUM(gI1,0,1,-(A0(M2Ah(gI1))*(3*sqr(g1) + 5*sqr(g2))*(sqr(ZA(gI1,0)) - sqr(ZA(gI1,1))))/20.)/2.;
   se11 += -SUM(gI1,0,1,-(vd*A0(M2Ah(gI1))*(3*sqr(g1) + 5*sqr(g2))*(sqr(ZA(gI1,0)) - sqr(ZA(gI1,1))))/20.)/(2.*vd);
   se11 += 3*SUM(gI1,0,1,(A0(M2Sb(gI1))*((sqr(g1) + 5*sqr(g2) - 20*sqr(yb))*sqr(ZB(gI1,0)) + 2*(sqr(g1) - 10*sqr(yb))*sqr(ZB(gI1,1))))/20.);
   se11 += 3*SUM(gI1,0,1,(A0(M2Sc(gI1))*((sqr(g1) - 5*sqr(g2))*sqr(ZC(gI1,0)) - 4*sqr(g1*ZC(gI1,1))))/20.);
   se11 += (-3*SUM(gI1,0,1,(vd*A0(M2Sc(gI1))*((sqr(g1) - 5*sqr(g2))*sqr(ZC(gI1,0)) - 4*sqr(g1*ZC(gI1,1))))/20.))/vd;
   se11 += 3*SUM(gI1,0,1,(A0(M2Sd(gI1))*((sqr(g1) + 5*sqr(g2))*sqr(ZD(gI1,0)) + 2*sqr(g1*ZD(gI1,1))))/20.);
   se11 += (-3*SUM(gI1,0,1,(vd*A0(M2Sd(gI1))*((sqr(g1) + 5*sqr(g2))*sqr(ZD(gI1,0)) + 2*sqr(g1*ZD(gI1,1))))/20.))/vd;
   se11 += SUM(gI1,0,1,-(A0(M2Se(gI1))*((3*sqr(g1) - 5*sqr(g2))*sqr(ZE(gI1,0)) - 6*sqr(g1*ZE(gI1,1))))/20.);
   se11 += -(SUM(gI1,0,1,-(vd*A0(M2Se(gI1))*((3*sqr(g1) - 5*sqr(g2))*sqr(ZE(gI1,0)) - 6*sqr(g1*ZE(gI1,1))))/20.)/vd);
   se11 += SUM(gI1,0,1,-(A0(M2hh(gI1))*(3*sqr(g1) + 5*sqr(g2))*(3*sqr(ZH(gI1,0)) - sqr(ZH(gI1,1))))/20.)/2.;
   se11 += SUM(gI1,0,1,-(A0(M2Sm(gI1))*((3*sqr(g1) - 5*sqr(g2))*sqr(ZM(gI1,0)) - 6*sqr(g1*ZM(gI1,1))))/20.);
   se11 += -(SUM(gI1,0,1,-(vd*A0(M2Sm(gI1))*((3*sqr(g1) - 5*sqr(g2))*sqr(ZM(gI1,0)) - 6*sqr(g1*ZM(gI1,1))))/20.)/vd);
   se11 += SUM(gI1,0,1,-(A0(M2Hpm(gI1))*((3*sqr(g1) + 5*sqr(g2))*sqr(ZP(gI1,0)) + (-3*sqr(g1) + 5*sqr(g2))*sqr(ZP(gI1,1))))/20.);
   se11 += 3*SUM(gI1,0,1,(A0(M2Ss(gI1))*((sqr(g1) + 5*sqr(g2))*sqr(ZS(gI1,0)) + 2*sqr(g1*ZS(gI1,1))))/20.);
   se11 += (-3*SUM(gI1,0,1,(vd*A0(M2Ss(gI1))*((sqr(g1) + 5*sqr(g2))*sqr(ZS(gI1,0)) + 2*sqr(g1*ZS(gI1,1))))/20.))/vd;
   se11 += 3*SUM(gI1,0,1,(A0(M2St(gI1))*((sqr(g1) - 5*sqr(g2))*sqr(ZT(gI1,0)) - 4*sqr(g1*ZT(gI1,1))))/20.);
   se11 += SUM(gI1,0,1,-(A0(M2Stau(gI1))*((3*sqr(g1) - 5*sqr(g2) + 20*sqr(ytau))*sqr(ZTau(gI1,0)) + 2*(-3*sqr(g1) + 10*sqr(ytau))*sqr(ZTau(gI1,1))))/20.);
   se11 += 3*SUM(gI1,0,1,(A0(M2Su(gI1))*((sqr(g1) - 5*sqr(g2))*sqr(ZU(gI1,0)) - 4*sqr(g1*ZU(gI1,1))))/20.);
   se11 += (-3*SUM(gI1,0,1,(vd*A0(M2Su(gI1))*((sqr(g1) - 5*sqr(g2))*sqr(ZU(gI1,0)) - 4*sqr(g1*ZU(gI1,1))))/20.))/vd;
   se11 += -SUM(gI1,0,1,SUM(gI2,0,1,-(sqr(g2)*(A0(sqr(MCha(gI1))) + A0(sqr(MCha(gI2))) + B0(p2,sqr(MCha(gI1)),sqr(MCha(gI2)))*(-p2 + sqr(MCha(gI1)) + sqr(MCha(gI2))))*(sqr(UM(gI2,1)*UP(gI1,0)) + sqr(UM(gI1,1)*UP(gI2,0))))/2.));
   se11 += -SUM(gI1,0,1,SUM(gI2,0,1,(B0(p2,M2Ah(gI1),M2Ah(gI2))*sqr(vd*(3*sqr(g1) + 5*sqr(g2))*(ZA(gI1,0)*ZA(gI2,0) - ZA(gI1,1)*ZA(gI2,1))))/400.))/2.;
   se11 += -3*SUM(gI1,0,1,SUM(gI2,0,1,(B0(p2,M2Sb(gI1),M2Sb(gI2))*sqr(ZB(gI1,0)*(vd*(sqr(g1) + 5*sqr(g2) - 20*sqr(yb))*ZB(gI2,0) - 10*Ab*sqrt2*yb*ZB(gI2,1)) + 2*ZB(gI1,1)*(-5*Ab*sqrt2*yb*ZB(gI2,0) + vd*(sqr(g1) - 10*sqr(yb))*ZB(gI2,1))))/400.));
   se11 += -3*SUM(gI1,0,1,SUM(gI2,0,1,(B0(p2,M2Sc(gI1),M2Sc(gI2))*sqr(vd*((sqr(g1) - 5*sqr(g2))*ZC(gI1,0)*ZC(gI2,0) - 4*sqr(g1)*ZC(gI1,1)*ZC(gI2,1))))/400.));
   se11 += -3*SUM(gI1,0,1,SUM(gI2,0,1,(B0(p2,M2Sd(gI1),M2Sd(gI2))*sqr(vd*((sqr(g1) + 5*sqr(g2))*ZD(gI1,0)*ZD(gI2,0) + 2*sqr(g1)*ZD(gI1,1)*ZD(gI2,1))))/400.));
   se11 += -SUM(gI1,0,1,SUM(gI2,0,1,(B0(p2,M2Se(gI1),M2Se(gI2))*sqr(vd*((3*sqr(g1) - 5*sqr(g2))*ZE(gI1,0)*ZE(gI2,0) - 6*sqr(g1)*ZE(gI1,1)*ZE(gI2,1))))/400.));
   se11 += -SUM(gI1,0,1,SUM(gI2,0,1,(B0(p2,M2hh(gI1),M2hh(gI2))*sqr((3*sqr(g1) + 5*sqr(g2))*(ZH(gI1,1)*(vu*ZH(gI2,0) + vd*ZH(gI2,1)) + ZH(gI1,0)*(-3*vd*ZH(gI2,0) + vu*ZH(gI2,1)))))/400.))/2.;
   se11 += -SUM(gI1,0,1,SUM(gI2,0,1,(B0(p2,M2Sm(gI1),M2Sm(gI2))*sqr(vd*((3*sqr(g1) - 5*sqr(g2))*ZM(gI1,0)*ZM(gI2,0) - 6*sqr(g1)*ZM(gI1,1)*ZM(gI2,1))))/400.));
   se11 += -SUM(gI1,0,1,SUM(gI2,0,1,(B0(p2,M2Hpm(gI1),M2Hpm(gI2))*sqr(ZP(gI1,0)*(vd*(3*sqr(g1) + 5*sqr(g2))*ZP(gI2,0) + 5*vu*sqr(g2)*ZP(gI2,1)) + ZP(gI1,1)*(5*vu*sqr(g2)*ZP(gI2,0) + vd*(-3*sqr(g1) + 5*sqr(g2))*ZP(gI2,1))))/400.));
   se11 += -3*SUM(gI1,0,1,SUM(gI2,0,1,(B0(p2,M2Ss(gI1),M2Ss(gI2))*sqr(vd*((sqr(g1) + 5*sqr(g2))*ZS(gI1,0)*ZS(gI2,0) + 2*sqr(g1)*ZS(gI1,1)*ZS(gI2,1))))/400.));
   se11 += -3*SUM(gI1,0,1,SUM(gI2,0,1,(B0(p2,M2St(gI1),M2St(gI2))*sqr(ZT(gI1,0)*(vd*(sqr(g1) - 5*sqr(g2))*ZT(gI2,0) + 10*mu*sqrt2*yt*ZT(gI2,1)) + 2*ZT(gI1,1)*(5*mu*sqrt2*yt*ZT(gI2,0) - 2*vd*sqr(g1)*ZT(gI2,1))))/400.));
   se11 += -SUM(gI1,0,1,SUM(gI2,0,1,(B0(p2,M2Stau(gI1),M2Stau(gI2))*sqr(ZTau(gI1,0)*(vd*(3*sqr(g1) - 5*sqr(g2) + 20*sqr(ytau))*ZTau(gI2,0) + 10*Atau*sqrt2*ytau*ZTau(gI2,1)) + 2*ZTau(gI1,1)*(5*Atau*sqrt2*ytau*ZTau(gI2,0) + vd*(-3*sqr(g1) + 10*sqr(ytau))*ZTau(gI2,1))))/400.));
   se11 += -3*SUM(gI1,0,1,SUM(gI2,0,1,(B0(p2,M2Su(gI1),M2Su(gI2))*sqr(vd*((sqr(g1) - 5*sqr(g2))*ZU(gI1,0)*ZU(gI2,0) - 4*sqr(g1)*ZU(gI1,1)*ZU(gI2,1))))/400.));
   se11 += 2*SUM(gI1,0,1,MCha(gI1)*SUM(gI2,0,1,B0(p2,sqr(MCha(gI1)),sqr(MCha(gI2)))*MCha(gI2)*sqr(g2)*UM(gI1,1)*UM(gI2,1)*UP(gI1,0)*UP(gI2,0)));
   se11 += (2*SUM(gI1,0,1,-(g2*sqrt2*A0(sqr(MCha(gI1)))*MCha(gI1)*UM(gI1,1)*UP(gI1,0))))/vd;
   se11 += (-3*SUM(gI1,0,1,(A0(M2Sb(gI1))*(vd*(sqr(g1) + 5*sqr(g2) - 20*sqr(yb))*sqr(ZB(gI1,0)) + 2*vd*(sqr(g1) - 10*sqr(yb))*sqr(ZB(gI1,1)) - 20*Ab*sqrt2*yb*ZB(gI1,0)*ZB(gI1,1)))/20.))/vd;
   se11 += -SUM(gI1,0,1,-(A0(M2hh(gI1))*(3*sqr(g1) + 5*sqr(g2))*(3*vd*sqr(ZH(gI1,0)) - vd*sqr(ZH(gI1,1)) - 2*vu*ZH(gI1,0)*ZH(gI1,1)))/20.)/(2.*vd);
   se11 += -(SUM(gI1,0,1,-(A0(M2Hpm(gI1))*(vd*(3*sqr(g1) + 5*sqr(g2))*sqr(ZP(gI1,0)) + vd*(-3*sqr(g1) + 5*sqr(g2))*sqr(ZP(gI1,1)) + 10*vu*sqr(g2)*ZP(gI1,0)*ZP(gI1,1)))/20.)/vd);
   se11 += (-3*SUM(gI1,0,1,(A0(M2St(gI1))*(vd*(sqr(g1) - 5*sqr(g2))*sqr(ZT(gI1,0)) - 4*vd*sqr(g1*ZT(gI1,1)) + 20*mu*sqrt2*yt*ZT(gI1,0)*ZT(gI1,1)))/20.))/vd;
   se11 += -(SUM(gI1,0,1,-(A0(M2Stau(gI1))*(vd*(3*sqr(g1) - 5*sqr(g2) + 20*sqr(ytau))*sqr(ZTau(gI1,0)) + 2*vd*(-3*sqr(g1) + 10*sqr(ytau))*sqr(ZTau(gI1,1)) + 20*Atau*sqrt2*ytau*ZTau(gI1,0)*ZTau(gI1,1)))/20.)/vd);
   se11 += SUM(gI1,0,3,MChi(gI1)*SUM(gI2,0,3,(B0(p2,sqr(MChi(gI1)),sqr(MChi(gI2)))*MChi(gI2)*sqr(ZN(gI1,2)*(g1*sqrt15*ZN(gI2,0) - 5*g2*ZN(gI2,1)) + (g1*sqrt15*ZN(gI1,0) - 5*g2*ZN(gI1,1))*ZN(gI2,2)))/50.));
   se11 += -SUM(gI1,0,3,SUM(gI2,0,3,-((A0(sqr(MChi(gI1))) + A0(sqr(MChi(gI2))) + B0(p2,sqr(MChi(gI1)),sqr(MChi(gI2)))*(-p2 + sqr(MChi(gI1)) + sqr(MChi(gI2))))*sqr(ZN(gI1,2)*(g1*sqrt15*ZN(gI2,0) - 5*g2*ZN(gI2,1)) + (g1*sqrt15*ZN(gI1,0) - 5*g2*ZN(gI1,1))*ZN(gI2,2)))/50.))/2.;
   se11 += SUM(gI1,0,3,(2*A0(sqr(MChi(gI1)))*MChi(gI1)*(g1*sqrt15*ZN(gI1,0) - 5*g2*ZN(gI1,1))*ZN(gI1,2))/5.)/vd;
   se11 += -SUM(gI2,0,1,-((2*A0(M2VZ) - A0(M2Ah(gI2)) + B0(p2,M2VZ,M2Ah(gI2))*(-M2VZ + 2*p2 + 2*M2Ah(gI2)))*(3*sqr(g1) + 5*sqr(g2))*sqr(ZA(gI2,0)))/20.);
   se11 += -2*SUM(gI2,0,1,-((2*A0(M2VWm) - A0(M2Hpm(gI2)) + B0(p2,M2VWm,M2Hpm(gI2))*(-M2VWm + 2*p2 + 2*M2Hpm(gI2)))*sqr(g2*ZP(gI2,0)))/4.);

   se12 += (-7*vd*vu*B0(p2,M2VWm,M2VWm)*pow4(g2))/8.;
   se12 += (vd*vu*B0(p2,M2SveL,M2SveL)*sqr(3*sqr(g1) + 5*sqr(g2)))/400.;
   se12 += (vd*vu*B0(p2,M2SvmL,M2SvmL)*sqr(3*sqr(g1) + 5*sqr(g2)))/400.;
   se12 += (vd*vu*B0(p2,M2SvtL,M2SvtL)*sqr(3*sqr(g1) + 5*sqr(g2)))/400.;
   se12 += (-7*vd*vu*B0(p2,M2VZ,M2VZ)*sqr(3*sqr(g1) + 5*sqr(g2)))/400.;
   se12 += -SUM(gI1,0,1,SUM(gI2,0,1,-(vd*vu*B0(p2,M2Ah(gI1),M2Ah(gI2))*sqr((3*sqr(g1) + 5*sqr(g2))*(ZA(gI1,0)*ZA(gI2,0) - ZA(gI1,1)*ZA(gI2,1))))/400.))/2.;
   se12 += -3*SUM(gI1,0,1,SUM(gI2,0,1,-(vd*vu*B0(p2,M2Sc(gI1),M2Sc(gI2))*sqr((sqr(g1) - 5*sqr(g2))*ZC(gI1,0)*ZC(gI2,0) - 4*sqr(g1)*ZC(gI1,1)*ZC(gI2,1)))/400.));
   se12 += -3*SUM(gI1,0,1,SUM(gI2,0,1,-(vd*vu*B0(p2,M2Sd(gI1),M2Sd(gI2))*sqr((sqr(g1) + 5*sqr(g2))*ZD(gI1,0)*ZD(gI2,0) + 2*sqr(g1)*ZD(gI1,1)*ZD(gI2,1)))/400.));
   se12 += -SUM(gI1,0,1,SUM(gI2,0,1,-(vd*vu*B0(p2,M2Se(gI1),M2Se(gI2))*sqr((3*sqr(g1) - 5*sqr(g2))*ZE(gI1,0)*ZE(gI2,0) - 6*sqr(g1)*ZE(gI1,1)*ZE(gI2,1)))/400.));
   se12 += -SUM(gI1,0,1,SUM(gI2,0,1,-(vd*vu*B0(p2,M2Sm(gI1),M2Sm(gI2))*sqr((3*sqr(g1) - 5*sqr(g2))*ZM(gI1,0)*ZM(gI2,0) - 6*sqr(g1)*ZM(gI1,1)*ZM(gI2,1)))/400.));
   se12 += -3*SUM(gI1,0,1,SUM(gI2,0,1,-(vd*vu*B0(p2,M2Ss(gI1),M2Ss(gI2))*sqr((sqr(g1) + 5*sqr(g2))*ZS(gI1,0)*ZS(gI2,0) + 2*sqr(g1)*ZS(gI1,1)*ZS(gI2,1)))/400.));
   se12 += -3*SUM(gI1,0,1,SUM(gI2,0,1,-(vd*vu*B0(p2,M2Su(gI1),M2Su(gI2))*sqr((sqr(g1) - 5*sqr(g2))*ZU(gI1,0)*ZU(gI2,0) - 4*sqr(g1)*ZU(gI1,1)*ZU(gI2,1)))/400.));
   se12 += 2*SUM(gI1,0,1,MCha(gI1)*SUM(gI2,0,1,(B0(p2,sqr(MCha(gI1)),sqr(MCha(gI2)))*MCha(gI2)*sqr(g2)*(UM(gI1,1)*UM(gI2,0)*UP(gI1,1)*UP(gI2,0) + UM(gI1,0)*UM(gI2,1)*UP(gI1,0)*UP(gI2,1)))/2.));
   se12 += -SUM(gI1,0,1,SUM(gI2,0,1,-(sqr(g2)*(A0(sqr(MCha(gI1))) + A0(sqr(MCha(gI2))) + B0(p2,sqr(MCha(gI1)),sqr(MCha(gI2)))*(-p2 + sqr(MCha(gI1)) + sqr(MCha(gI2))))*(UM(gI2,0)*UM(gI2,1)*UP(gI1,0)*UP(gI1,1) + UM(gI1,0)*UM(gI1,1)*UP(gI2,0)*UP(gI2,1)))/2.));
   se12 += -3*SUM(gI1,0,1,SUM(gI2,0,1,-(B0(p2,M2Sb(gI1),M2Sb(gI2))*(ZB(gI1,0)*(vu*(sqr(g1) + 5*sqr(g2))*ZB(gI2,0) - 10*mu*sqrt2*yb*ZB(gI2,1)) + 2*ZB(gI1,1)*(-5*mu*sqrt2*yb*ZB(gI2,0) + vu*sqr(g1)*ZB(gI2,1)))*(ZB(gI1,0)*(vd*(sqr(g1) + 5*sqr(g2) - 20*sqr(yb))*ZB(gI2,0) - 10*Ab*sqrt2*yb*ZB(gI2,1)) + 2*ZB(gI1,1)*(-5*Ab*sqrt2*yb*ZB(gI2,0) + vd*(sqr(g1) - 10*sqr(yb))*ZB(gI2,1))))/400.));
   se12 += -SUM(gI1,0,1,SUM(gI2,0,1,(B0(p2,M2hh(gI1),M2hh(gI2))*sqr(3*sqr(g1) + 5*sqr(g2))*(ZH(gI1,0)*(vu*ZH(gI2,0) + vd*ZH(gI2,1)) + ZH(gI1,1)*(vd*ZH(gI2,0) - 3*vu*ZH(gI2,1)))*(ZH(gI1,1)*(vu*ZH(gI2,0) + vd*ZH(gI2,1)) + ZH(gI1,0)*(-3*vd*ZH(gI2,0) + vu*ZH(gI2,1))))/400.))/2.;
   se12 += -SUM(gI1,0,1,SUM(gI2,0,1,(B0(p2,M2Hpm(gI1),M2Hpm(gI2))*(ZP(gI1,0)*(vd*(3*sqr(g1) + 5*sqr(g2))*ZP(gI2,0) + 5*vu*sqr(g2)*ZP(gI2,1)) + ZP(gI1,1)*(5*vu*sqr(g2)*ZP(gI2,0) + vd*(-3*sqr(g1) + 5*sqr(g2))*ZP(gI2,1)))*(ZP(gI1,0)*((-3*vu*sqr(g1) + 5*vu*sqr(g2))*ZP(gI2,0) + 5*vd*sqr(g2)*ZP(gI2,1)) + ZP(gI1,1)*(5*vd*sqr(g2)*ZP(gI2,0) + vu*(3*sqr(g1) + 5*sqr(g2))*ZP(gI2,1))))/400.));
   se12 += -3*SUM(gI1,0,1,SUM(gI2,0,1,-(B0(p2,M2St(gI1),M2St(gI2))*(ZT(gI1,0)*(vd*(sqr(g1) - 5*sqr(g2))*ZT(gI2,0) + 10*mu*sqrt2*yt*ZT(gI2,1)) + 2*ZT(gI1,1)*(5*mu*sqrt2*yt*ZT(gI2,0) - 2*vd*sqr(g1)*ZT(gI2,1)))*(ZT(gI1,0)*(vu*(sqr(g1) - 5*sqr(g2) + 20*sqr(yt))*ZT(gI2,0) + 10*At*sqrt2*yt*ZT(gI2,1)) + 2*ZT(gI1,1)*(5*At*sqrt2*yt*ZT(gI2,0) - 2*vu*(sqr(g1) - 5*sqr(yt))*ZT(gI2,1))))/400.));
   se12 += -SUM(gI1,0,1,SUM(gI2,0,1,-(B0(p2,M2Stau(gI1),M2Stau(gI2))*(ZTau(gI1,0)*(vu*(3*sqr(g1) - 5*sqr(g2))*ZTau(gI2,0) + 10*mu*sqrt2*ytau*ZTau(gI2,1)) + 2*ZTau(gI1,1)*(5*mu*sqrt2*ytau*ZTau(gI2,0) - 3*vu*sqr(g1)*ZTau(gI2,1)))*(ZTau(gI1,0)*(vd*(3*sqr(g1) - 5*sqr(g2) + 20*sqr(ytau))*ZTau(gI2,0) + 10*Atau*sqrt2*ytau*ZTau(gI2,1)) + 2*ZTau(gI1,1)*(5*Atau*sqrt2*ytau*ZTau(gI2,0) + vd*(-3*sqr(g1) + 10*sqr(ytau))*ZTau(gI2,1))))/400.));
   se12 += SUM(gI1,0,1,(A0(M2hh(gI1))*(3*sqr(g1) + 5*sqr(g2))*ZH(gI1,0)*ZH(gI1,1))/10.)/2.;
   se12 += SUM(gI1,0,1,-(A0(M2Hpm(gI1))*sqr(g2)*ZP(gI1,0)*ZP(gI1,1))/2.);
   se12 += SUM(gI1,0,3,MChi(gI1)*SUM(gI2,0,3,-(B0(p2,sqr(MChi(gI1)),sqr(MChi(gI2)))*MChi(gI2)*(ZN(gI1,2)*(g1*sqrt15*ZN(gI2,0) - 5*g2*ZN(gI2,1)) + (g1*sqrt15*ZN(gI1,0) - 5*g2*ZN(gI1,1))*ZN(gI2,2))*(ZN(gI1,3)*(g1*sqrt15*ZN(gI2,0) - 5*g2*ZN(gI2,1)) + (g1*sqrt15*ZN(gI1,0) - 5*g2*ZN(gI1,1))*ZN(gI2,3)))/50.));
   se12 += -SUM(gI1,0,3,SUM(gI2,0,3,((A0(sqr(MChi(gI1))) + A0(sqr(MChi(gI2))) + B0(p2,sqr(MChi(gI1)),sqr(MChi(gI2)))*(-p2 + sqr(MChi(gI1)) + sqr(MChi(gI2))))*(ZN(gI1,2)*(g1*sqrt15*ZN(gI2,0) - 5*g2*ZN(gI2,1)) + (g1*sqrt15*ZN(gI1,0) - 5*g2*ZN(gI1,1))*ZN(gI2,2))*(ZN(gI1,3)*(g1*sqrt15*ZN(gI2,0) - 5*g2*ZN(gI2,1)) + (g1*sqrt15*ZN(gI1,0) - 5*g2*ZN(gI1,1))*ZN(gI2,3)))/50.))/2.;
   se12 += -SUM(gI2,0,1,((2*A0(M2VZ) - A0(M2Ah(gI2)) + B0(p2,M2VZ,M2Ah(gI2))*(-M2VZ + 2*p2 + 2*M2Ah(gI2)))*(3*sqr(g1) + 5*sqr(g2))*ZA(gI2,0)*ZA(gI2,1))/20.);
   se12 += -2*SUM(gI2,0,1,((2*A0(M2VWm) - A0(M2Hpm(gI2)) + B0(p2,M2VWm,M2Hpm(gI2))*(-M2VWm + 2*p2 + 2*M2Hpm(gI2)))*sqr(g2)*ZP(gI2,0)*ZP(gI2,1))/4.);

   se22 += (A0(M2VZ)*(-3*sqr(g1) - 5*sqr(g2)))/20.;
   se22 += -(A0(M2VWm)*sqr(g2))/2.;
   se22 += (-7*B0(p2,M2VWm,M2VWm)*pow4(g2)*sqr(vu))/8.;
   se22 += -3*B0(p2,sqr(MFt),sqr(MFt))*(p2 - 4*sqr(MFt))*sqr(yt);
   se22 += -(B0(p2,M2SveL,M2SveL)*sqr(vu*(3*sqr(g1) + 5*sqr(g2))))/400.;
   se22 += -(B0(p2,M2SvmL,M2SvmL)*sqr(vu*(3*sqr(g1) + 5*sqr(g2))))/400.;
   se22 += -(B0(p2,M2SvtL,M2SvtL)*sqr(vu*(3*sqr(g1) + 5*sqr(g2))))/400.;
   se22 += (-7*B0(p2,M2VZ,M2VZ)*sqr(vu*(3*sqr(g1) + 5*sqr(g2))))/400.;
   se22 += SUM(gI1,0,1,(A0(M2Ah(gI1))*(3*sqr(g1) + 5*sqr(g2))*(sqr(ZA(gI1,0)) - sqr(ZA(gI1,1))))/20.)/2.;
   se22 += -SUM(gI1,0,1,(vu*A0(M2Ah(gI1))*(3*sqr(g1) + 5*sqr(g2))*(sqr(ZA(gI1,0)) - sqr(ZA(gI1,1))))/20.)/(2.*vu);
   se22 += 3*SUM(gI1,0,1,-(A0(M2Sb(gI1))*((sqr(g1) + 5*sqr(g2))*sqr(ZB(gI1,0)) + 2*sqr(g1*ZB(gI1,1))))/20.);
   se22 += 3*SUM(gI1,0,1,-(A0(M2Sc(gI1))*((sqr(g1) - 5*sqr(g2))*sqr(ZC(gI1,0)) - 4*sqr(g1*ZC(gI1,1))))/20.);
   se22 += (-3*SUM(gI1,0,1,-(vu*A0(M2Sc(gI1))*((sqr(g1) - 5*sqr(g2))*sqr(ZC(gI1,0)) - 4*sqr(g1*ZC(gI1,1))))/20.))/vu;
   se22 += 3*SUM(gI1,0,1,-(A0(M2Sd(gI1))*((sqr(g1) + 5*sqr(g2))*sqr(ZD(gI1,0)) + 2*sqr(g1*ZD(gI1,1))))/20.);
   se22 += (-3*SUM(gI1,0,1,-(vu*A0(M2Sd(gI1))*((sqr(g1) + 5*sqr(g2))*sqr(ZD(gI1,0)) + 2*sqr(g1*ZD(gI1,1))))/20.))/vu;
   se22 += SUM(gI1,0,1,(A0(M2Se(gI1))*((3*sqr(g1) - 5*sqr(g2))*sqr(ZE(gI1,0)) - 6*sqr(g1*ZE(gI1,1))))/20.);
   se22 += -(SUM(gI1,0,1,(vu*A0(M2Se(gI1))*((3*sqr(g1) - 5*sqr(g2))*sqr(ZE(gI1,0)) - 6*sqr(g1*ZE(gI1,1))))/20.)/vu);
   se22 += SUM(gI1,0,1,(A0(M2hh(gI1))*(3*sqr(g1) + 5*sqr(g2))*(sqr(ZH(gI1,0)) - 3*sqr(ZH(gI1,1))))/20.)/2.;
   se22 += SUM(gI1,0,1,(A0(M2Sm(gI1))*((3*sqr(g1) - 5*sqr(g2))*sqr(ZM(gI1,0)) - 6*sqr(g1*ZM(gI1,1))))/20.);
   se22 += -(SUM(gI1,0,1,(vu*A0(M2Sm(gI1))*((3*sqr(g1) - 5*sqr(g2))*sqr(ZM(gI1,0)) - 6*sqr(g1*ZM(gI1,1))))/20.)/vu);
   se22 += SUM(gI1,0,1,(A0(M2Hpm(gI1))*((3*sqr(g1) - 5*sqr(g2))*sqr(ZP(gI1,0)) - (3*sqr(g1) + 5*sqr(g2))*sqr(ZP(gI1,1))))/20.);
   se22 += 3*SUM(gI1,0,1,-(A0(M2Ss(gI1))*((sqr(g1) + 5*sqr(g2))*sqr(ZS(gI1,0)) + 2*sqr(g1*ZS(gI1,1))))/20.);
   se22 += (-3*SUM(gI1,0,1,-(vu*A0(M2Ss(gI1))*((sqr(g1) + 5*sqr(g2))*sqr(ZS(gI1,0)) + 2*sqr(g1*ZS(gI1,1))))/20.))/vu;
   se22 += 3*SUM(gI1,0,1,-(A0(M2St(gI1))*((sqr(g1) - 5*sqr(g2) + 20*sqr(yt))*sqr(ZT(gI1,0)) - 4*(sqr(g1) - 5*sqr(yt))*sqr(ZT(gI1,1))))/20.);
   se22 += SUM(gI1,0,1,(A0(M2Stau(gI1))*((3*sqr(g1) - 5*sqr(g2))*sqr(ZTau(gI1,0)) - 6*sqr(g1*ZTau(gI1,1))))/20.);
   se22 += 3*SUM(gI1,0,1,-(A0(M2Su(gI1))*((sqr(g1) - 5*sqr(g2))*sqr(ZU(gI1,0)) - 4*sqr(g1*ZU(gI1,1))))/20.);
   se22 += (-3*SUM(gI1,0,1,-(vu*A0(M2Su(gI1))*((sqr(g1) - 5*sqr(g2))*sqr(ZU(gI1,0)) - 4*sqr(g1*ZU(gI1,1))))/20.))/vu;
   se22 += -SUM(gI1,0,1,SUM(gI2,0,1,-(sqr(g2)*(A0(sqr(MCha(gI1))) + A0(sqr(MCha(gI2))) + B0(p2,sqr(MCha(gI1)),sqr(MCha(gI2)))*(-p2 + sqr(MCha(gI1)) + sqr(MCha(gI2))))*(sqr(UM(gI2,0)*UP(gI1,1)) + sqr(UM(gI1,0)*UP(gI2,1))))/2.));
   se22 += -SUM(gI1,0,1,SUM(gI2,0,1,(B0(p2,M2Ah(gI1),M2Ah(gI2))*sqr(vu*(3*sqr(g1) + 5*sqr(g2))*(ZA(gI1,0)*ZA(gI2,0) - ZA(gI1,1)*ZA(gI2,1))))/400.))/2.;
   se22 += -3*SUM(gI1,0,1,SUM(gI2,0,1,(B0(p2,M2Sb(gI1),M2Sb(gI2))*sqr(ZB(gI1,0)*(vu*(sqr(g1) + 5*sqr(g2))*ZB(gI2,0) - 10*mu*sqrt2*yb*ZB(gI2,1)) + 2*ZB(gI1,1)*(-5*mu*sqrt2*yb*ZB(gI2,0) + vu*sqr(g1)*ZB(gI2,1))))/400.));
   se22 += -3*SUM(gI1,0,1,SUM(gI2,0,1,(B0(p2,M2Sc(gI1),M2Sc(gI2))*sqr(vu*((sqr(g1) - 5*sqr(g2))*ZC(gI1,0)*ZC(gI2,0) - 4*sqr(g1)*ZC(gI1,1)*ZC(gI2,1))))/400.));
   se22 += -3*SUM(gI1,0,1,SUM(gI2,0,1,(B0(p2,M2Sd(gI1),M2Sd(gI2))*sqr(vu*((sqr(g1) + 5*sqr(g2))*ZD(gI1,0)*ZD(gI2,0) + 2*sqr(g1)*ZD(gI1,1)*ZD(gI2,1))))/400.));
   se22 += -SUM(gI1,0,1,SUM(gI2,0,1,(B0(p2,M2Se(gI1),M2Se(gI2))*sqr(vu*((3*sqr(g1) - 5*sqr(g2))*ZE(gI1,0)*ZE(gI2,0) - 6*sqr(g1)*ZE(gI1,1)*ZE(gI2,1))))/400.));
   se22 += -SUM(gI1,0,1,SUM(gI2,0,1,(B0(p2,M2hh(gI1),M2hh(gI2))*sqr((3*sqr(g1) + 5*sqr(g2))*(ZH(gI1,0)*(vu*ZH(gI2,0) + vd*ZH(gI2,1)) + ZH(gI1,1)*(vd*ZH(gI2,0) - 3*vu*ZH(gI2,1)))))/400.))/2.;
   se22 += -SUM(gI1,0,1,SUM(gI2,0,1,(B0(p2,M2Sm(gI1),M2Sm(gI2))*sqr(vu*((3*sqr(g1) - 5*sqr(g2))*ZM(gI1,0)*ZM(gI2,0) - 6*sqr(g1)*ZM(gI1,1)*ZM(gI2,1))))/400.));
   se22 += -SUM(gI1,0,1,SUM(gI2,0,1,(B0(p2,M2Hpm(gI1),M2Hpm(gI2))*sqr(ZP(gI1,0)*(vu*(3*sqr(g1) - 5*sqr(g2))*ZP(gI2,0) - 5*vd*sqr(g2)*ZP(gI2,1)) - ZP(gI1,1)*(5*vd*sqr(g2)*ZP(gI2,0) + vu*(3*sqr(g1) + 5*sqr(g2))*ZP(gI2,1))))/400.));
   se22 += -3*SUM(gI1,0,1,SUM(gI2,0,1,(B0(p2,M2Ss(gI1),M2Ss(gI2))*sqr(vu*((sqr(g1) + 5*sqr(g2))*ZS(gI1,0)*ZS(gI2,0) + 2*sqr(g1)*ZS(gI1,1)*ZS(gI2,1))))/400.));
   se22 += -3*SUM(gI1,0,1,SUM(gI2,0,1,(B0(p2,M2St(gI1),M2St(gI2))*sqr(ZT(gI1,0)*(vu*(sqr(g1) - 5*sqr(g2) + 20*sqr(yt))*ZT(gI2,0) + 10*At*sqrt2*yt*ZT(gI2,1)) + 2*ZT(gI1,1)*(5*At*sqrt2*yt*ZT(gI2,0) - 2*vu*(sqr(g1) - 5*sqr(yt))*ZT(gI2,1))))/400.));
   se22 += -SUM(gI1,0,1,SUM(gI2,0,1,(B0(p2,M2Stau(gI1),M2Stau(gI2))*sqr(ZTau(gI1,0)*(vu*(3*sqr(g1) - 5*sqr(g2))*ZTau(gI2,0) + 10*mu*sqrt2*ytau*ZTau(gI2,1)) + 2*ZTau(gI1,1)*(5*mu*sqrt2*ytau*ZTau(gI2,0) - 3*vu*sqr(g1)*ZTau(gI2,1))))/400.));
   se22 += -3*SUM(gI1,0,1,SUM(gI2,0,1,(B0(p2,M2Su(gI1),M2Su(gI2))*sqr(vu*((sqr(g1) - 5*sqr(g2))*ZU(gI1,0)*ZU(gI2,0) - 4*sqr(g1)*ZU(gI1,1)*ZU(gI2,1))))/400.));
   se22 += 2*SUM(gI1,0,1,MCha(gI1)*SUM(gI2,0,1,B0(p2,sqr(MCha(gI1)),sqr(MCha(gI2)))*MCha(gI2)*sqr(g2)*UM(gI1,0)*UM(gI2,0)*UP(gI1,1)*UP(gI2,1)));
   se22 += (2*SUM(gI1,0,1,-(g2*sqrt2*A0(sqr(MCha(gI1)))*MCha(gI1)*UM(gI1,0)*UP(gI1,1))))/vu;
   se22 += (-3*SUM(gI1,0,1,-(A0(M2Sb(gI1))*(vu*(sqr(g1) + 5*sqr(g2))*sqr(ZB(gI1,0)) + 2*vu*sqr(g1*ZB(gI1,1)) - 20*mu*sqrt2*yb*ZB(gI1,0)*ZB(gI1,1)))/20.))/vu;
   se22 += -SUM(gI1,0,1,(A0(M2hh(gI1))*(3*sqr(g1) + 5*sqr(g2))*(vu*sqr(ZH(gI1,0)) - 3*vu*sqr(ZH(gI1,1)) + 2*vd*ZH(gI1,0)*ZH(gI1,1)))/20.)/(2.*vu);
   se22 += -(SUM(gI1,0,1,(A0(M2Hpm(gI1))*(vu*(3*sqr(g1) - 5*sqr(g2))*sqr(ZP(gI1,0)) - vu*(3*sqr(g1) + 5*sqr(g2))*sqr(ZP(gI1,1)) - 10*vd*sqr(g2)*ZP(gI1,0)*ZP(gI1,1)))/20.)/vu);
   se22 += (-3*SUM(gI1,0,1,-(A0(M2St(gI1))*(vu*(sqr(g1) - 5*sqr(g2) + 20*sqr(yt))*sqr(ZT(gI1,0)) - 4*vu*(sqr(g1) - 5*sqr(yt))*sqr(ZT(gI1,1)) + 20*At*sqrt2*yt*ZT(gI1,0)*ZT(gI1,1)))/20.))/vu;
   se22 += -(SUM(gI1,0,1,(A0(M2Stau(gI1))*(vu*(3*sqr(g1) - 5*sqr(g2))*sqr(ZTau(gI1,0)) - 6*vu*sqr(g1*ZTau(gI1,1)) + 20*mu*sqrt2*ytau*ZTau(gI1,0)*ZTau(gI1,1)))/20.)/vu);
   se22 += SUM(gI1,0,3,MChi(gI1)*SUM(gI2,0,3,(B0(p2,sqr(MChi(gI1)),sqr(MChi(gI2)))*MChi(gI2)*sqr(ZN(gI1,3)*(g1*sqrt15*ZN(gI2,0) - 5*g2*ZN(gI2,1)) + (g1*sqrt15*ZN(gI1,0) - 5*g2*ZN(gI1,1))*ZN(gI2,3)))/50.));
   se22 += -SUM(gI1,0,3,SUM(gI2,0,3,-((A0(sqr(MChi(gI1))) + A0(sqr(MChi(gI2))) + B0(p2,sqr(MChi(gI1)),sqr(MChi(gI2)))*(-p2 + sqr(MChi(gI1)) + sqr(MChi(gI2))))*sqr(ZN(gI1,3)*(g1*sqrt15*ZN(gI2,0) - 5*g2*ZN(gI2,1)) + (g1*sqrt15*ZN(gI1,0) - 5*g2*ZN(gI1,1))*ZN(gI2,3)))/50.))/2.;
   se22 += SUM(gI1,0,3,(-2*A0(sqr(MChi(gI1)))*MChi(gI1)*(g1*sqrt15*ZN(gI1,0) - 5*g2*ZN(gI1,1))*ZN(gI1,3))/5.)/vu;
   se22 += -SUM(gI2,0,1,-((2*A0(M2VZ) - A0(M2Ah(gI2)) + B0(p2,M2VZ,M2Ah(gI2))*(-M2VZ + 2*p2 + 2*M2Ah(gI2)))*(3*sqr(g1) + 5*sqr(g2))*sqr(ZA(gI2,1)))/20.);
   se22 += -2*SUM(gI2,0,1,-((2*A0(M2VWm) - A0(M2Hpm(gI2)) + B0(p2,M2VWm,M2Hpm(gI2))*(-M2VWm + 2*p2 + 2*M2Hpm(gI2)))*sqr(g2*ZP(gI2,1)))/4.);

   RM22 se(RM22::Zero());
   se << se11, se12, se12, se22;

   return se * one_loop;
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
   const auto yb     = pars.Yd(2,2);
   const auto ytau   = pars.Ye(2,2);
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
   const auto yb     = pars.Yd(2,2);
   const auto ytau   = pars.Ye(2,2);
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

   se11 += -3*B0(0,sqr(MFb),sqr(MFb))*sqr(yb);
   se11 += 12*D1B0(sqr(MFb),sqr(MFb))*sqr(MFb*yb);
   se11 += -(B0(0,sqr(MFtau),sqr(MFtau))*sqr(ytau));
   se11 += 4*D1B0(sqr(MFtau),sqr(MFtau))*sqr(MFtau*ytau);
   se11 += -3*D1B0(M2Sb(0),M2Sb(0))*sqr(yb*(vd*yb*(sqr(ZB(0,0)) + sqr(ZB(0,1))) + Ab*sqrt2*ZB(0,0)*ZB(0,1)));
   se11 += -3*D1B0(M2Sb(1),M2Sb(1))*sqr(yb*(vd*yb*(sqr(ZB(1,0)) + sqr(ZB(1,1))) + Ab*sqrt2*ZB(1,0)*ZB(1,1)));
   se11 += (-3*D1B0(M2Sb(0),M2Sb(1))*sqr(yb*(Ab*sqrt2*(ZB(0,1)*ZB(1,0) + ZB(0,0)*ZB(1,1)) + 2*vd*yb*(ZB(0,0)*ZB(1,0) + ZB(0,1)*ZB(1,1)))))/4.;
   se11 += (-3*D1B0(M2Sb(1),M2Sb(0))*sqr(yb*(Ab*sqrt2*(ZB(0,1)*ZB(1,0) + ZB(0,0)*ZB(1,1)) + 2*vd*yb*(ZB(0,0)*ZB(1,0) + ZB(0,1)*ZB(1,1)))))/4.;
   se11 += -6*D1B0(M2St(0),M2St(0))*sqr(mu*yt*ZT(0,0)*ZT(0,1));
   se11 += -6*D1B0(M2St(1),M2St(1))*sqr(mu*yt*ZT(1,0)*ZT(1,1));
   se11 += (-3*D1B0(M2St(0),M2St(1))*sqr(mu*yt*(ZT(0,1)*ZT(1,0) + ZT(0,0)*ZT(1,1))))/2.;
   se11 += (-3*D1B0(M2St(1),M2St(0))*sqr(mu*yt*(ZT(0,1)*ZT(1,0) + ZT(0,0)*ZT(1,1))))/2.;
   se11 += -(D1B0(M2Stau(0),M2Stau(0))*sqr(ytau*(vd*ytau*(sqr(ZTau(0,0)) + sqr(ZTau(0,1))) + Atau*sqrt2*ZTau(0,0)*ZTau(0,1))));
   se11 += -(D1B0(M2Stau(1),M2Stau(1))*sqr(ytau*(vd*ytau*(sqr(ZTau(1,0)) + sqr(ZTau(1,1))) + Atau*sqrt2*ZTau(1,0)*ZTau(1,1))));
   se11 += -(D1B0(M2Stau(0),M2Stau(1))*sqr(ytau*(Atau*sqrt2*(ZTau(0,1)*ZTau(1,0) + ZTau(0,0)*ZTau(1,1)) + 2*vd*ytau*(ZTau(0,0)*ZTau(1,0) + ZTau(0,1)*ZTau(1,1)))))/4.;
   se11 += -(D1B0(M2Stau(1),M2Stau(0))*sqr(ytau*(Atau*sqrt2*(ZTau(0,1)*ZTau(1,0) + ZTau(0,0)*ZTau(1,1)) + 2*vd*ytau*(ZTau(0,0)*ZTau(1,0) + ZTau(0,1)*ZTau(1,1)))))/4.;

   se12 += 3*mu*D1B0(M2Sb(0),M2Sb(0))*sqr(yb)*ZB(0,0)*ZB(0,1)*(sqrt2*vd*yb*(sqr(ZB(0,0)) + sqr(ZB(0,1))) + 2*Ab*ZB(0,0)*ZB(0,1));
   se12 += 3*mu*D1B0(M2Sb(1),M2Sb(1))*sqr(yb)*ZB(1,0)*ZB(1,1)*(sqrt2*vd*yb*(sqr(ZB(1,0)) + sqr(ZB(1,1))) + 2*Ab*ZB(1,0)*ZB(1,1));
   se12 += (3*mu*D1B0(M2Sb(0),M2Sb(1))*sqr(yb)*(ZB(0,1)*ZB(1,0) + ZB(0,0)*ZB(1,1))*(Ab*(ZB(0,1)*ZB(1,0) + ZB(0,0)*ZB(1,1)) + sqrt2*vd*yb*(ZB(0,0)*ZB(1,0) + ZB(0,1)*ZB(1,1))))/2.;
   se12 += (3*mu*D1B0(M2Sb(1),M2Sb(0))*sqr(yb)*(ZB(0,1)*ZB(1,0) + ZB(0,0)*ZB(1,1))*(Ab*(ZB(0,1)*ZB(1,0) + ZB(0,0)*ZB(1,1)) + sqrt2*vd*yb*(ZB(0,0)*ZB(1,0) + ZB(0,1)*ZB(1,1))))/2.;
   se12 += 3*mu*D1B0(M2St(0),M2St(0))*sqr(yt)*ZT(0,0)*ZT(0,1)*(sqrt2*vu*yt*(sqr(ZT(0,0)) + sqr(ZT(0,1))) + 2*At*ZT(0,0)*ZT(0,1));
   se12 += 3*mu*D1B0(M2St(1),M2St(1))*sqr(yt)*ZT(1,0)*ZT(1,1)*(sqrt2*vu*yt*(sqr(ZT(1,0)) + sqr(ZT(1,1))) + 2*At*ZT(1,0)*ZT(1,1));
   se12 += (3*mu*D1B0(M2St(0),M2St(1))*sqr(yt)*(ZT(0,1)*ZT(1,0) + ZT(0,0)*ZT(1,1))*(At*(ZT(0,1)*ZT(1,0) + ZT(0,0)*ZT(1,1)) + sqrt2*vu*yt*(ZT(0,0)*ZT(1,0) + ZT(0,1)*ZT(1,1))))/2.;
   se12 += (3*mu*D1B0(M2St(1),M2St(0))*sqr(yt)*(ZT(0,1)*ZT(1,0) + ZT(0,0)*ZT(1,1))*(At*(ZT(0,1)*ZT(1,0) + ZT(0,0)*ZT(1,1)) + sqrt2*vu*yt*(ZT(0,0)*ZT(1,0) + ZT(0,1)*ZT(1,1))))/2.;
   se12 += mu*sqrt2*D1B0(M2Stau(0),M2Stau(0))*sqr(ytau)*ZTau(0,0)*ZTau(0,1)*(vd*ytau*(sqr(ZTau(0,0)) + sqr(ZTau(0,1))) + Atau*sqrt2*ZTau(0,0)*ZTau(0,1));
   se12 += mu*sqrt2*D1B0(M2Stau(1),M2Stau(1))*sqr(ytau)*ZTau(1,0)*ZTau(1,1)*(vd*ytau*(sqr(ZTau(1,0)) + sqr(ZTau(1,1))) + Atau*sqrt2*ZTau(1,0)*ZTau(1,1));
   se12 += (mu*sqrt2*D1B0(M2Stau(0),M2Stau(1))*sqr(ytau)*(ZTau(0,1)*ZTau(1,0) + ZTau(0,0)*ZTau(1,1))*(Atau*sqrt2*(ZTau(0,1)*ZTau(1,0) + ZTau(0,0)*ZTau(1,1)) + 2*vd*ytau*(ZTau(0,0)*ZTau(1,0) + ZTau(0,1)*ZTau(1,1))))/4.;
   se12 += (mu*sqrt2*D1B0(M2Stau(1),M2Stau(0))*sqr(ytau)*(ZTau(0,1)*ZTau(1,0) + ZTau(0,0)*ZTau(1,1))*(Atau*sqrt2*(ZTau(0,1)*ZTau(1,0) + ZTau(0,0)*ZTau(1,1)) + 2*vd*ytau*(ZTau(0,0)*ZTau(1,0) + ZTau(0,1)*ZTau(1,1))))/4.;

   se22 += -3*B0(0,sqr(MFt),sqr(MFt))*sqr(yt);
   se22 += 12*D1B0(sqr(MFt),sqr(MFt))*sqr(MFt*yt);
   se22 += -6*D1B0(M2Sb(0),M2Sb(0))*sqr(mu*yb*ZB(0,0)*ZB(0,1));
   se22 += -6*D1B0(M2Sb(1),M2Sb(1))*sqr(mu*yb*ZB(1,0)*ZB(1,1));
   se22 += (-3*D1B0(M2Sb(0),M2Sb(1))*sqr(mu*yb*(ZB(0,1)*ZB(1,0) + ZB(0,0)*ZB(1,1))))/2.;
   se22 += (-3*D1B0(M2Sb(1),M2Sb(0))*sqr(mu*yb*(ZB(0,1)*ZB(1,0) + ZB(0,0)*ZB(1,1))))/2.;
   se22 += -3*D1B0(M2St(0),M2St(0))*sqr(yt*(vu*yt*(sqr(ZT(0,0)) + sqr(ZT(0,1))) + At*sqrt2*ZT(0,0)*ZT(0,1)));
   se22 += -3*D1B0(M2St(1),M2St(1))*sqr(yt*(vu*yt*(sqr(ZT(1,0)) + sqr(ZT(1,1))) + At*sqrt2*ZT(1,0)*ZT(1,1)));
   se22 += (-3*D1B0(M2St(0),M2St(1))*sqr(yt*(At*sqrt2*(ZT(0,1)*ZT(1,0) + ZT(0,0)*ZT(1,1)) + 2*vu*yt*(ZT(0,0)*ZT(1,0) + ZT(0,1)*ZT(1,1)))))/4.;
   se22 += (-3*D1B0(M2St(1),M2St(0))*sqr(yt*(At*sqrt2*(ZT(0,1)*ZT(1,0) + ZT(0,0)*ZT(1,1)) + 2*vu*yt*(ZT(0,0)*ZT(1,0) + ZT(0,1)*ZT(1,1)))))/4.;
   se22 += -2*D1B0(M2Stau(0),M2Stau(0))*sqr(mu*ytau*ZTau(0,0)*ZTau(0,1));
   se22 += -2*D1B0(M2Stau(1),M2Stau(1))*sqr(mu*ytau*ZTau(1,0)*ZTau(1,1));
   se22 += -(D1B0(M2Stau(0),M2Stau(1))*sqr(mu*ytau*(ZTau(0,1)*ZTau(1,0) + ZTau(0,0)*ZTau(1,1))))/2.;
   se22 += -(D1B0(M2Stau(1),M2Stau(0))*sqr(mu*ytau*(ZTau(0,1)*ZTau(1,0) + ZTau(0,0)*ZTau(1,1))))/2.;

   RM22 se;
   se << se11, se12, se12, se22;

   return se * one_loop;
}

RM22 MSSM_mass_eigenstates::get_mass_matrix_hh() const
{
   return masses.get_mass_matrix_hh(pars);
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

} // namespace mh1l
} // namespace himalaya
