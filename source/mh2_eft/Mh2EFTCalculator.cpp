#define Pi M_PI

#include "Mh2EFTCalculator.hpp"

template <typename T> T pow2(T x)  { return x*x; }
template <typename T> T pow3(T x)  { return x*x*x; }
template <typename T> T pow4(T x)  { return x*x*x*x; }
template <typename T> T pow5(T x)  { return x*x*x*x*x; }
template <typename T> T pow6(T x)  { return x*x*x*x*x*x; }
template <typename T> T pow7(T x)  { return x*x*x*x*x*x*x; }
template <typename T> T pow8(T x)  { return x*x*x*x*x*x*x*x; }
template <typename T> T pow9(T x)  { return x*x*x*x*x*x*x*x*x; }
template <typename T> T pow11(T x) { return pow2(x)*pow9(x); }
template <typename T> T pow12(T x) { return pow2(pow6(x)); }
template <typename T> T pow13(T x) { return pow4(x)*pow9(x); }
template <typename T> T pow14(T x) { return pow2(pow7(x)); }
template <typename T> T pow15(T x) { return pow6(x)*pow9(x); }
template <typename T> T pow16(T x) { return pow2(pow8(x)); }
template <typename T> T pow18(T x) { return pow2(pow9(x)); }

namespace himalaya {
namespace mh2_eft {

   const double Mh2EFTCalculator::zt2 = 1.6449340668482264364724151666460;
   const double Mh2EFTCalculator::zt3 = 1.2020569031595942853997381615114;
   const double Mh2EFTCalculator::log2 = std::log(2.);
}
}
   
//himalaya::mh2_eft::Mh2EFTCalculator::Mh2EFTCalculator(const Parameters& p){
//}


   /**
    * 	fin[] function from arXiv:hep-ph/0507139 .
    *
    * 	@param m12 squared mass \f$m_1^2\f$
    * 	@param m22 squared mass \f$m_2^2\f$
    * 	@param MR2 squared renormalization scale
    *
    * 	@return fin(m12, m22)
    */
   double himalaya::mh2_eft::Mh2EFTCalculator::fin(double m12, double m22, double MR2)
   {
      using std::log;
      using gm2calc::dilog;
      return (6*(m12*log(m12/MR2) + m22*log(m22/MR2)) +
         (-m12 - m22)*(7 + zt2) +
         (m12 - m22)*(2*dilog(1 - m12/m22) +
            pow2(log(m12/m22))/2.) +
         ((m12 + m22)*pow2(log(m12/m22)))/2. -
         2*(m12*pow2(log(m12/MR2)) + m22*pow2(log(m22/MR2))))/2.;
   }

/// 1-loop coefficient O(at*log^0)
double himalaya::mh2_eft::Mh2EFTCalculator::coeff_as_0_log_0(double mQ32, double mU32, double Xt, double MR2)
{
   using std::log;

   const double result = 
      (12*pow4(Xt))/pow2(mQ32 - mU32) + log(mQ32/MR2)*(6 + (12*pow2(Xt
         ))/(mQ32 - mU32) + ((-6*mQ32)/pow3(mQ32 - mU32) - (6*mU32)/pow3(mQ32 -
         mU32))*pow4(Xt)) + log(mU32/MR2)*(6 - (12*pow2(Xt))/(mQ32 - mU32) + (
         (6*mQ32)/pow3(mQ32 - mU32) + (6*mU32)/pow3(mQ32 - mU32))*pow4(Xt));

   return result;
}

/// 1-loop coefficient O(at*log^1)
double himalaya::mh2_eft::Mh2EFTCalculator::coeff_as_0_log_1()
{
   return 12;
}

/// 2-loop coefficient O(at*as*log^0)
double himalaya::mh2_eft::Mh2EFTCalculator::coeff_as_1_log_0(double mQ32, double mU32, double Xt, double m3, double MR2)
{
   using std::log;
   using gm2calc::dilog;
   const double m32 = pow2(m3);

   const double result = 
      (16*(-(mQ32*mU32) + pow2(m32)))/((m32 - mQ32)*(m32 - mU32)) - (
         64*m32*pow2(Xt))/(mQ32*mU32) + (256*Xt*(-(mU32*dilog(1 - m32/mQ32)) +
         m32*(dilog(1 - m32/mQ32) - dilog(1 - m32/mU32)) + mQ32*dilog(1 -
         m32/mU32))*pow3(m3))/((m32 - mQ32)*(m32 - mU32)*(mQ32 - mU32)) + (128*
         m3*(2*m32*(dilog(1 - m32/mQ32) - dilog(1 - m32/mU32)) + mQ32*(-4 -
         dilog(1 - m32/mQ32) + dilog(1 - m32/mU32)) + mU32*(4 - dilog(1 -
         m32/mQ32) + dilog(1 - m32/mU32)))*pow3(Xt))/pow3(mQ32 - mU32) + (32*(2
         *m32*mQ32*mU32*dilog(1 - m32/mQ32) - 2*m32*mQ32*mU32*dilog(1 -
         m32/mU32) + m32*pow2(mQ32) + 6*mU32*pow2(mQ32) - mU32*dilog(1 -
         m32/mQ32)*pow2(mQ32) + mU32*dilog(1 - m32/mU32)*pow2(mQ32) - m32*pow2(
         mU32) - 6*mQ32*pow2(mU32) - mQ32*dilog(1 - m32/mQ32)*pow2(mU32) + mQ32
         *dilog(1 - m32/mU32)*pow2(mU32))*pow4(Xt))/(mQ32*mU32*pow3(mQ32 - mU32
         )) + log(mU32/MR2)*((32*(-2*m32*mU32 + 3*pow2(m32) + pow2(mU32)))/pow2
         (m32 - mU32) + (64*(3*m32 - 2*mU32)*pow2(Xt))/((m32 - mU32)*(-mQ32 +
         mU32)) - (64*Xt*(-(m3*mU32) + 2*pow3(m3)))/((m32 - mU32)*(-mQ32 + mU32
         )) + (128*(m3*mQ32 + 3*m3*mU32)*pow3(Xt))/pow3(-mQ32 + mU32) - (32*(2*
         m32 + 3*mQ32 + 7*mU32)*pow4(Xt))/pow3(-mQ32 + mU32)) + pow2(log(
         mQ32/MR2))*((-16*(-2*m32*mQ32 + 3*pow2(m32) + pow2(mQ32)))/pow2(m32 -
         mQ32) - (32*(3*mQ32 - mU32)*pow2(Xt))/pow2(mQ32 - mU32) + (128*Xt*pow3
         (m3))/((m32 - mQ32)*(mQ32 - mU32)) + (64*m3*(2*m32 - mQ32 - mU32)*pow3
         (Xt))/pow3(mQ32 - mU32) + (32*(m32*mQ32 - m32*mU32 + 2*mQ32*mU32 + 2*
         pow2(mQ32))*pow4(Xt))/pow4(mQ32 - mU32)) + pow2(log(mU32/MR2))*((-16*(
         -2*m32*mU32 + 3*pow2(m32) + pow2(mU32)))/pow2(m32 - mU32) - (32*(-mQ32
         + 3*mU32)*pow2(Xt))/pow2(-mQ32 + mU32) + (128*Xt*pow3(m3))/((m32 -
         mU32)*(-mQ32 + mU32)) + (64*m3*(2*m32 - mQ32 - mU32)*pow3(Xt))/pow3(
         -mQ32 + mU32) + (32*(-(m32*mQ32) + m32*mU32 + 2*mQ32*mU32 + 2*pow2(
         mU32))*pow4(Xt))/pow4(-mQ32 + mU32)) + log(mQ32/MR2)*((32*(-2*m32*mQ32
         + 3*pow2(m32) + pow2(mQ32)))/pow2(m32 - mQ32) + (64*(3*m32 - 2*mQ32)*
         pow2(Xt))/((m32 - mQ32)*(mQ32 - mU32)) - (64*Xt*(-(m3*mQ32) + 2*pow3(
         m3)))/((m32 - mQ32)*(mQ32 - mU32)) + (128*(3*m3*mQ32 + m3*mU32)*pow3(
         Xt))/pow3(mQ32 - mU32) - (32*(2*m32 + 7*mQ32 + 3*mU32)*pow4(Xt))/pow3(
         mQ32 - mU32) + log(mU32/MR2)*((64*(mQ32 + mU32)*pow2(Xt))/pow2(mQ32 -
         mU32) - (64*(2*mQ32*mU32 + pow2(mQ32) + pow2(mU32))*pow4(Xt))/pow4(
         mQ32 - mU32))) + (16*(9*m32*(mQ32 + mU32)*pow2(mQ32)*pow2(mU32) - 3*
         pow3(mQ32)*pow3(mU32) + pow3(m32)*(mU32*(19 + 8*dilog(1 - m32/mU32))*
         pow2(mQ32) + mQ32*(19 + 8*dilog(1 - m32/mQ32))*pow2(mU32) + 2*pow3(
         mQ32) + 2*pow3(mU32)) - 4*pow2(m32)*(6*pow2(mQ32)*pow2(mU32) + mU32*(2
         + dilog(1 - m32/mU32))*pow3(mQ32) + mQ32*(2 + dilog(1 - m32/mQ32))*
         pow3(mU32)) - (mQ32*mU32*(13 + 4*dilog(1 - m32/mQ32) + 4*dilog(1 -
         m32/mU32)) + 4*pow2(mQ32) + 4*pow2(mU32))*pow4(m32) + 2*(mQ32 + mU32)*
         pow5(m32)))/(mQ32*mU32*pow2(m32 - mQ32)*pow2(m32 - mU32)) + log(
         m32/MR2)*((64*(-m32 + mQ32 + mU32)*pow2(m32)*pow2(Xt))/(mQ32*(-m32 +
         mQ32)*(m32 - mU32)*mU32) + (64*Xt*pow3(m3))/((m32 - mQ32)*(m32 - mU32)
         ) + log(mQ32/MR2)*((32*pow2(m32))/pow2(m32 - mQ32) - (128*Xt*pow3(m3))
         /((m32 - mQ32)*(mQ32 - mU32)) - (256*pow3(m3)*pow3(Xt))/pow3(mQ32 -
         mU32)) + log(mU32/MR2)*((32*pow2(m32))/pow2(m32 - mU32) - (128*Xt*pow3
         (m3))/((m32 - mU32)*(-mQ32 + mU32)) - (256*pow3(m3)*pow3(Xt))/pow3(
         -mQ32 + mU32)) - (32*(m32*mQ32 + m32*mU32)*pow4(Xt))/(mQ32*mU32*pow2(
         mQ32 - mU32)) - (32*(-2*mU32*pow2(mQ32)*pow3(m32) - 2*mQ32*pow2(mU32)*
         pow3(m32) + mU32*pow2(m32)*pow3(mQ32) + pow3(m32)*pow3(mQ32) + mQ32*
         pow2(m32)*pow3(mU32) + pow3(m32)*pow3(mU32) + 2*mQ32*mU32*pow4(m32) -
         2*pow2(mQ32)*pow4(m32) - 2*pow2(mU32)*pow4(m32) + mQ32*pow5(m32) +
         mU32*pow5(m32)))/(mQ32*mU32*pow2(m32 - mQ32)*pow2(m32 - mU32)));

   return result;
}

/// 2-loop coefficient O(at*as*log^1)
double himalaya::mh2_eft::Mh2EFTCalculator::coeff_as_1_log_1(double mQ32, double mU32, double Xt, double m3, double MR2)
{
   using std::log;
   const double m32 = pow2(m3);

   const double result = 
      (-32*(-(m32*mQ32) - m32*mU32 + 2*pow2(m32)))/((m32 - mQ32)*(m32
         - mU32)) + log(mQ32/MR2)*((128*m3*mQ32*Xt)/((m32 - mQ32)*(mQ32 - mU32)
         ) - (64*m32*mQ32)/pow2(m32 - mQ32) + (32*pow2(mQ32))/pow2(m32 - mQ32))
         + log(mU32/MR2)*((128*m3*mU32*Xt)/((m32 - mU32)*(-mQ32 + mU32)) - (64
         *m32*mU32)/pow2(m32 - mU32) + (32*pow2(mU32))/pow2(m32 - mU32)) + log(
         m32/MR2)*((-128*Xt*pow3(m3))/((m32 - mQ32)*(m32 - mU32)) + (32*(pow2(
         m32)*pow2(mQ32) + pow2(m32)*pow2(mU32) - 2*mQ32*pow3(m32) - 2*mU32*
         pow3(m32) + 2*pow4(m32)))/(pow2(m32 - mQ32)*pow2(m32 - mU32)));

   return result;
}

/// 2-loop coefficient O(at*as*log^2)
double himalaya::mh2_eft::Mh2EFTCalculator::coeff_as_1_log_2()
{
   return 96;
}

/// 3-loop coefficient O(at*as^2*log^0)
double himalaya::mh2_eft::Mh2EFTCalculator::coeff_as_2_log_0(double mQ32, double mU32, double Xt, double m3, double msq2, double MR2)
{
   using std::log;
   using gm2calc::dilog;
   const double m32 = pow2(m3);

   const double deltalambda3Loop = 0.;
   const double SMConst3Loop = 0.;

   const double result =
      deltalambda3Loop + SMConst3Loop - (176*(-(mQ32*mU32) + pow2(m32)
         ))/(3.*(m32 - mQ32)*(m32 - mU32)) - (64*m32*(3*m32*mQ32 + 3*m32*mU32 +
         5*mQ32*mU32 - 11*pow2(m32))*pow2(Xt))/(3.*(m32 - mQ32)*mQ32*(m32 -
         mU32)*mU32) - (48*pow2(-(mQ32*mU32) + pow2(m32)))/(pow2(m32 - mQ32)*
         pow2(m32 - mU32)) + (128*m3*(2*m32*(dilog(1 - m32/mQ32) - dilog(1 -
         m32/mU32)) + mQ32*(-4 - dilog(1 - m32/mQ32) + dilog(1 - m32/mU32)) +
         mU32*(4 - dilog(1 - m32/mQ32) + dilog(1 - m32/mU32)))*(3*m32*mQ32 + 3*
         m32*mU32 + 5*mQ32*mU32 - 11*pow2(m32))*pow3(Xt))/(3.*(m32 - mQ32)*(m32
         - mU32)*pow3(mQ32 - mU32)) + Xt*((256*(-(mU32*dilog(1 - m32/mQ32)) +
         m32*(dilog(1 - m32/mQ32) - dilog(1 - m32/mU32)) + mQ32*dilog(1 -
         m32/mU32))*(3*m32*mQ32 + 3*m32*mU32 + 5*mQ32*mU32 - 11*pow2(m32))*pow3
         (m3))/(3.*(mQ32 - mU32)*pow2(m32 - mQ32)*pow2(m32 - mU32)) + 24*((64*(
         m3*fin(mQ32,m32,MR2) - m3*fin(mQ32,mU32,MR2) + m3*fin(mU32,m32,MR2) +
         3*pow3(m3) + zt2*pow3(m3)))/(9.*(m32 - mQ32)*(m32 - mU32)) - (4*((20*(
         m3*mQ32 - m3*msq2)*fin(mQ32,msq2,MR2))/((mQ32 - mU32)*pow2(m32 - mQ32)
         ) - (2*m3*(mQ32 - mU32)*fin(mU32,m32,MR2))/((m32 - mU32)*pow2(m32 -
         mQ32)) + (2*m3*(mQ32 - mU32)*fin(mQ32,m32,MR2))/((m32 - mQ32)*pow2(m32
         - mU32)) + (20*(-(m3*msq2) + m3*mU32)*fin(mU32,msq2,MR2))/((-mQ32 +
         mU32)*pow2(m32 - mU32)) + (2*(6*m3*mQ32 + 60*m3*msq2 + 6*m3*mU32 + m3*
         mQ32*zt2 + 10*m3*msq2*zt2 + m3*mU32*zt2 + 24*pow3(m3)))/((m32 - mQ32)*
         (m32 - mU32)) - (20*fin(m32,msq2,MR2)*(m3*mQ32*msq2 - 2*m3*mQ32*mU32 +
         m3*msq2*mU32 + mQ32*pow3(m3) - 2*msq2*pow3(m3) + mU32*pow3(m3)))/(
         pow2(m32 - mQ32)*pow2(m32 - mU32)) + 3*((-2*m3*fin(mQ32,mU32,MR2))/((
         m32 - mQ32)*(m32 - mU32)) + (2*fin(mQ32,m32,MR2)*(m3*mQ32 - 3*m3*mU32
         + 2*pow3(m3)))/((m32 - mQ32)*(m32 - mU32)*(mQ32 - mU32)) - (2*fin(mU32
         ,m32,MR2)*(-3*m3*mQ32 + m3*mU32 + 2*pow3(m3)))/((m32 - mQ32)*(m32 -
         mU32)*(mQ32 - mU32)) + (2*(11*pow3(m3) + 3*zt2*pow3(m3)))/((m32 - mQ32
         )*(m32 - mU32))) + (2*fin(mQ32,mU32,MR2)*(m3*pow2(mQ32) + m3*pow2(mU32
         ) - 2*mQ32*pow3(m3) - 2*mU32*pow3(m3) + 2*pow5(m3)))/(pow2(m32 - mQ32)
         *pow2(m32 - mU32))))/3.)) + (16*(3*m32*mQ32 + 3*m32*mU32 + 5*mQ32*mU32
         - 11*pow2(m32))*(9*m32*(mQ32 + mU32)*pow2(mQ32)*pow2(mU32) - 3*pow3(
         mQ32)*pow3(mU32) + pow3(m32)*(mU32*(19 + 8*dilog(1 - m32/mU32))*pow2(
         mQ32) + mQ32*(19 + 8*dilog(1 - m32/mQ32))*pow2(mU32) + 2*pow3(mQ32) +
         2*pow3(mU32)) - 4*pow2(m32)*(6*pow2(mQ32)*pow2(mU32) + mU32*(2 + dilog
         (1 - m32/mU32))*pow3(mQ32) + mQ32*(2 + dilog(1 - m32/mQ32))*pow3(mU32)
         ) - (mQ32*mU32*(13 + 4*dilog(1 - m32/mQ32) + 4*dilog(1 - m32/mU32)) +
         4*pow2(mQ32) + 4*pow2(mU32))*pow4(m32) + 2*(mQ32 + mU32)*pow5(m32)))/(
         3.*mQ32*mU32*pow3(m32 - mQ32)*pow3(m32 - mU32)) + 24*(
         -4.111111111111111 - (2*mQ32)/(3.*(-m32 + mQ32)) - (2*mU32)/(3.*(-m32
         + mU32)) + pow2(-2 + (2*mQ32)/(3.*(-m32 + mQ32)) + (2*mU32)/(3.*(-m32
         + mU32))) - (16*(-29.625 - (21*mQ32)/(4.*(-m32 + mQ32)) - (56*mQ32)/(
         mQ32 - mU32) + (56*mU32)/(mQ32 - mU32) + (111*mQ32)/(4.*(-m32 + mU32))
         + (45*mU32)/(2.*(-m32 + mU32)) + zt2/2. - (5*mQ32*zt2)/(2.*(-m32 +
         mQ32)) - (8*mQ32*zt2)/(mQ32 - mU32) + (8*mU32*zt2)/(mQ32 - mU32) + (4*
         mQ32*zt2)/(-m32 + mU32) + (3*mU32*zt2)/(2.*(-m32 + mU32)) + (2*(mQ32 +
         mU32)*fin(mQ32,mU32,MR2))/((m32 - mQ32)*(m32 - mU32)) + (56*pow2(mQ32
         ))/((-m32 + mQ32)*(mQ32 - mU32)) - (111*pow2(mQ32))/(4.*(-m32 + mQ32)*
         (-m32 + mU32)) + (8*zt2*pow2(mQ32))/((-m32 + mQ32)*(mQ32 - mU32)) - (4
         *zt2*pow2(mQ32))/((-m32 + mQ32)*(-m32 + mU32)) + (23*pow2(mQ32))/(8.*
         pow2(-m32 + mQ32)) + (3*zt2*pow2(mQ32))/(2.*pow2(-m32 + mQ32)) - (56*
         pow2(mU32))/((mQ32 - mU32)*(-m32 + mU32)) - (8*zt2*pow2(mU32))/((mQ32
         - mU32)*(-m32 + mU32)) + (23*pow2(mU32))/(8.*pow2(-m32 + mU32)) + (3*
         zt2*pow2(mU32))/(2.*pow2(-m32 + mU32)) + fin(mQ32,m32,MR2)*(-5/(2.*(
         -m32 + mQ32)) - 8/(mQ32 - mU32) + (8*mQ32)/((-m32 + mQ32)*(mQ32 - mU32
         )) + 2/(-m32 + mU32) - (4*mQ32)/((-m32 + mQ32)*(-m32 + mU32)) - (3*
         pow2(mQ32))/pow3(-m32 + mQ32)) - (21*pow3(mQ32))/pow3(-m32 + mQ32) - (
         3*zt2*pow3(mQ32))/pow3(-m32 + mQ32) + fin(mU32,m32,MR2)*(-2/(-m32 +
         mQ32) + 8/(mQ32 - mU32) + 3/(2.*(-m32 + mU32)) - (4*mQ32)/((-m32 +
         mQ32)*(-m32 + mU32)) - (8*mU32)/((mQ32 - mU32)*(-m32 + mU32)) - (3*
         pow2(mU32))/pow3(-m32 + mU32)) - (21*pow3(mU32))/pow3(-m32 + mU32) - (
         3*zt2*pow3(mU32))/pow3(-m32 + mU32)))/9. - (4*(-16.333333333333332 + (
         31*mQ32)/(-m32 + mQ32) + (45*msq2)/(-m32 + mQ32) + (9*mU32)/(2.*(-m32
         + mQ32)) + (9*mQ32)/(2.*(-m32 + mU32)) + (45*msq2)/(-m32 + mU32) + (31
         *mU32)/(-m32 + mU32) - 3*zt2 + (mQ32*zt2)/(4.*(-m32 + mQ32)) + (15*
         msq2*zt2)/(2.*(-m32 + mQ32)) + (3*mU32*zt2)/(4.*(-m32 + mQ32)) + (3*
         mQ32*zt2)/(4.*(-m32 + mU32)) + (15*msq2*zt2)/(2.*(-m32 + mU32)) + (
         mU32*zt2)/(4.*(-m32 + mU32)) - (27275*mQ32*msq2)/(432.*pow2(-m32 +
         mQ32)) - (2825*mQ32*mU32)/(432.*pow2(-m32 + mQ32)) - (65*mQ32*msq2*zt2
         )/(6.*pow2(-m32 + mQ32)) - (mQ32*mU32*zt2)/pow2(-m32 + mQ32) - (23*
         pow2(mQ32))/pow2(-m32 + mQ32) + (1355*msq2*pow2(mQ32))/(432.*(mQ32 -
         msq2)*pow2(-m32 + mQ32)) + (233*mU32*pow2(mQ32))/(432.*(mQ32 - mU32)*
         pow2(-m32 + mQ32)) + (5*msq2*zt2*pow2(mQ32))/(6.*(mQ32 - msq2)*pow2(
         -m32 + mQ32)) - (515*mQ32*pow2(msq2))/(144.*(mQ32 - msq2)*pow2(-m32 +
         mQ32)) - (5*mQ32*zt2*pow2(msq2))/(3.*(mQ32 - msq2)*pow2(-m32 + mQ32))
         + (95*pow2(mQ32)*pow2(msq2))/(216.*pow2(-m32 + mQ32)*pow2(mQ32 - msq2)
         ) + (5*zt2*pow2(mQ32)*pow2(msq2))/(6.*pow2(-m32 + mQ32)*pow2(mQ32 -
         msq2)) - (233*mQ32*pow2(mU32))/(288.*(mQ32 - mU32)*pow2(-m32 + mQ32))
         + (233*pow2(mQ32)*pow2(mU32))/(864.*pow2(-m32 + mQ32)*pow2(mQ32 - mU32
         )) + 3*(16.09722222222222 - (22*mQ32)/(-m32 + mQ32) - (14*mQ32)/(-m32
         + mU32) - (36*mU32)/(-m32 + mU32) + zt2 - (3*mQ32*zt2)/(2.*(-m32 +
         mQ32)) - (2*mQ32*zt2)/(-m32 + mU32) - (7*mU32*zt2)/(2.*(-m32 + mU32))
         - ((mQ32 + mU32)*fin(mQ32,mU32,MR2))/((m32 - mQ32)*(m32 - mU32)) + (14
         *pow2(mQ32))/((-m32 + mQ32)*(-m32 + mU32)) + (2*zt2*pow2(mQ32))/((-m32
         + mQ32)*(-m32 + mU32)) + (fin(mQ32,m32,MR2)*(3*m32*mQ32 - 3*mQ32*mU32
         + pow2(m32) - pow2(mQ32)))/((m32 - mU32)*pow2(-m32 + mQ32)) + (31*
         pow2(mQ32))/pow2(-m32 + mQ32) + (3*zt2*pow2(mQ32))/pow2(-m32 + mQ32) +
         (fin(mU32,m32,MR2)*(3*m32*mU32 - 3*mQ32*mU32 + pow2(m32) - pow2(mU32)
         ))/((m32 - mQ32)*pow2(m32 - mU32)) + (31*pow2(mU32))/pow2(-m32 + mU32)
         + (3*zt2*pow2(mU32))/pow2(-m32 + mU32)) - (6*mQ32*mU32)/pow2(-m32 +
         mU32) - (27275*msq2*mU32)/(432.*pow2(-m32 + mU32)) - (mQ32*mU32*zt2)
         /pow2(-m32 + mU32) - (65*msq2*mU32*zt2)/(6.*pow2(-m32 + mU32)) - (515*
         mU32*pow2(msq2))/(144.*(-msq2 + mU32)*pow2(-m32 + mU32)) - (5*mU32*zt2
         *pow2(msq2))/(3.*(-msq2 + mU32)*pow2(-m32 + mU32)) - (23*pow2(mU32))
         /pow2(-m32 + mU32) + (1355*msq2*pow2(mU32))/(432.*(-msq2 + mU32)*pow2(
         -m32 + mU32)) + (5*msq2*zt2*pow2(mU32))/(6.*(-msq2 + mU32)*pow2(-m32 +
         mU32)) + (95*pow2(msq2)*pow2(mU32))/(216.*pow2(-m32 + mU32)*pow2(
         -msq2 + mU32)) + (5*zt2*pow2(msq2)*pow2(mU32))/(6.*pow2(-m32 + mU32)*
         pow2(-msq2 + mU32)) - (5*fin(mQ32,msq2,MR2)*(-3*m32*mQ32 + 3*m32*msq2
         + mQ32*msq2 - 2*pow2(m32) + pow2(mQ32)))/(2.*pow3(-m32 + mQ32)) + (fin
         (mU32,m32,MR2)*(6*m32*mQ32*mU32 - 8*mQ32*pow2(m32) + 2*mU32*pow2(m32)
         + 3*m32*pow2(mQ32) - 3*m32*pow2(mU32) - mQ32*pow2(mU32) + 2*pow3(m32)
         - pow3(mQ32)))/(4.*(m32 - mU32)*pow3(-m32 + mQ32)) - (95*mQ32*pow3(
         msq2))/(216.*pow2(-m32 + mQ32)*pow2(mQ32 - msq2)) - (5*mQ32*zt2*pow3(
         msq2))/(6.*pow2(-m32 + mQ32)*pow2(mQ32 - msq2)) - (95*mU32*pow3(msq2))
         /(216.*pow2(-m32 + mU32)*pow2(-msq2 + mU32)) - (5*mU32*zt2*pow3(msq2))
         /(6.*pow2(-m32 + mU32)*pow2(-msq2 + mU32)) + (271*pow3(mU32))/(864.*(
         mQ32 - mU32)*pow2(-m32 + mQ32)) + (zt2*pow3(mU32))/(12.*(mQ32 - mU32)*
         pow2(-m32 + mQ32)) - (7*mQ32*pow3(mU32))/(12.*pow2(-m32 + mQ32)*pow2(
         mQ32 - mU32)) - (mQ32*zt2*pow3(mU32))/(12.*pow2(-m32 + mQ32)*pow2(mQ32
         - mU32)) + (fin(mQ32,m32,MR2)*(-6*m32*mQ32*mU32 - 2*mQ32*pow2(m32) +
         8*mU32*pow2(m32) + 3*m32*pow2(mQ32) + mU32*pow2(mQ32) - 3*m32*pow2(
         mU32) - 2*pow3(m32) + pow3(mU32)))/(4.*(m32 - mQ32)*pow3(m32 - mU32))
         + fin(m32,msq2,MR2)*(5/(2.*(-m32 + mQ32)) + 5/(2.*(-m32 + mU32)) + (15
         *mQ32)/(2.*pow2(-m32 + mQ32)) - (15*msq2)/(2.*pow2(-m32 + mQ32)) - (15
         *msq2)/(2.*pow2(-m32 + mU32)) + (15*mU32)/(2.*pow2(-m32 + mU32)) + (10
         *mQ32*msq2)/pow3(-m32 + mQ32) - (10*pow2(mQ32))/pow3(-m32 + mQ32) + (
         10*msq2*mU32)/pow3(-m32 + mU32) - (10*pow2(mU32))/pow3(-m32 + mU32)) -
         (5*fin(mU32,msq2,MR2)*(3*m32*msq2 - 3*m32*mU32 + msq2*mU32 - 2*pow2(
         m32) + pow2(mU32)))/(2.*pow3(-m32 + mU32)) + (271*pow4(mU32))/(864.*
         pow2(-m32 + mQ32)*pow2(mQ32 - mU32)) + (zt2*pow4(mU32))/(12.*pow2(-m32
         + mQ32)*pow2(mQ32 - mU32)) + (fin(mQ32,mU32,MR2)*(15*mU32*pow2(m32)*
         pow2(mQ32) + 15*mQ32*pow2(m32)*pow2(mU32) - 6*m32*pow2(mQ32)*pow2(mU32
         ) - 20*mQ32*mU32*pow3(m32) + 14*pow2(mQ32)*pow3(m32) + 14*pow2(mU32)*
         pow3(m32) - 6*m32*mU32*pow3(mQ32) - 11*pow2(m32)*pow3(mQ32) + pow2(
         mU32)*pow3(mQ32) - 6*m32*mQ32*pow3(mU32) - 11*pow2(m32)*pow3(mU32) +
         pow2(mQ32)*pow3(mU32) - 6*mQ32*pow4(m32) - 6*mU32*pow4(m32) + 3*m32*
         pow4(mQ32) + mU32*pow4(mQ32) + 3*m32*pow4(mU32) + mQ32*pow4(mU32) + 4*
         pow5(m32)))/(4.*pow3(-m32 + mQ32)*pow3(m32 - mU32))))/3.) + pow4(Xt)*(
         (32*(3*m32*mQ32 + 3*m32*mU32 + 5*mQ32*mU32 - 11*pow2(m32))*((2*mQ32*
         mU32*(-dilog(1 - m32/mQ32) + dilog(1 - m32/mU32)) - pow2(mQ32) + pow2(
         mU32))*pow3(m32) + (-5 + dilog(1 - m32/mQ32) - dilog(1 - m32/mU32))*
         pow2(mU32)*pow3(mQ32) + pow2(m32)*(3*mU32*(-2 + dilog(1 - m32/mQ32) -
         dilog(1 - m32/mU32))*pow2(mQ32) + 3*mQ32*(2 + dilog(1 - m32/mQ32) -
         dilog(1 - m32/mU32))*pow2(mU32) + pow3(mQ32) - pow3(mU32)) + (5 +
         dilog(1 - m32/mQ32) - dilog(1 - m32/mU32))*pow2(mQ32)*pow3(mU32) + m32
         *(4*(-dilog(1 - m32/mQ32) + dilog(1 - m32/mU32))*pow2(mQ32)*pow2(mU32)
         + mU32*(5 - dilog(1 - m32/mQ32) + dilog(1 - m32/mU32))*pow3(mQ32) +
         mQ32*(-5 - dilog(1 - m32/mQ32) + dilog(1 - m32/mU32))*pow3(mU32))))/(
         3.*(m32 - mQ32)*mQ32*(-m32 + mQ32)*mU32*pow2(m32 - mU32)*pow3(mQ32 -
         mU32)) + (12*((8*(-(mQ32*mU32) + pow2(m32)))/(3.*(m32 - mQ32)*(m32 -
         mU32)) + (40*pow2(-(mQ32*mU32) + pow2(m32)))/(9.*pow2(m32 - mQ32)*pow2
         (m32 - mU32)) - 4*(-4.111111111111111 - (2*mQ32)/(3.*(-m32 + mQ32)) -
         (2*mU32)/(3.*(-m32 + mU32)) + pow2(-2 + (2*mQ32)/(3.*(-m32 + mQ32)) +
         (2*mU32)/(3.*(-m32 + mU32))) - (16*(-29.625 - (21*mQ32)/(4.*(-m32 +
         mQ32)) - (56*mQ32)/(mQ32 - mU32) + (56*mU32)/(mQ32 - mU32) + (111*mQ32
         )/(4.*(-m32 + mU32)) + (45*mU32)/(2.*(-m32 + mU32)) + zt2/2. - (5*mQ32
         *zt2)/(2.*(-m32 + mQ32)) - (8*mQ32*zt2)/(mQ32 - mU32) + (8*mU32*zt2)/(
         mQ32 - mU32) + (4*mQ32*zt2)/(-m32 + mU32) + (3*mU32*zt2)/(2.*(-m32 +
         mU32)) + (2*(mQ32 + mU32)*fin(mQ32,mU32,MR2))/((m32 - mQ32)*(m32 -
         mU32)) + (56*pow2(mQ32))/((-m32 + mQ32)*(mQ32 - mU32)) - (111*pow2(
         mQ32))/(4.*(-m32 + mQ32)*(-m32 + mU32)) + (8*zt2*pow2(mQ32))/((-m32 +
         mQ32)*(mQ32 - mU32)) - (4*zt2*pow2(mQ32))/((-m32 + mQ32)*(-m32 + mU32)
         ) + (23*pow2(mQ32))/(8.*pow2(-m32 + mQ32)) + (3*zt2*pow2(mQ32))/(2.*
         pow2(-m32 + mQ32)) - (56*pow2(mU32))/((mQ32 - mU32)*(-m32 + mU32)) - (
         8*zt2*pow2(mU32))/((mQ32 - mU32)*(-m32 + mU32)) + (23*pow2(mU32))/(8.*
         pow2(-m32 + mU32)) + (3*zt2*pow2(mU32))/(2.*pow2(-m32 + mU32)) + fin(
         mQ32,m32,MR2)*(-5/(2.*(-m32 + mQ32)) - 8/(mQ32 - mU32) + (8*mQ32)/((
         -m32 + mQ32)*(mQ32 - mU32)) + 2/(-m32 + mU32) - (4*mQ32)/((-m32 + mQ32
         )*(-m32 + mU32)) - (3*pow2(mQ32))/pow3(-m32 + mQ32)) - (21*pow3(mQ32))
         /pow3(-m32 + mQ32) - (3*zt2*pow3(mQ32))/pow3(-m32 + mQ32) + fin(mU32,
         m32,MR2)*(-2/(-m32 + mQ32) + 8/(mQ32 - mU32) + 3/(2.*(-m32 + mU32)) -
         (4*mQ32)/((-m32 + mQ32)*(-m32 + mU32)) - (8*mU32)/((mQ32 - mU32)*(-m32
         + mU32)) - (3*pow2(mU32))/pow3(-m32 + mU32)) - (21*pow3(mU32))/pow3(
         -m32 + mU32) - (3*zt2*pow3(mU32))/pow3(-m32 + mU32)))/9. - (4*(
         -16.333333333333332 + (31*mQ32)/(-m32 + mQ32) + (45*msq2)/(-m32 + mQ32
         ) + (9*mU32)/(2.*(-m32 + mQ32)) + (9*mQ32)/(2.*(-m32 + mU32)) + (45*
         msq2)/(-m32 + mU32) + (31*mU32)/(-m32 + mU32) - 3*zt2 + (mQ32*zt2)/(4.
         *(-m32 + mQ32)) + (15*msq2*zt2)/(2.*(-m32 + mQ32)) + (3*mU32*zt2)/(4.*
         (-m32 + mQ32)) + (3*mQ32*zt2)/(4.*(-m32 + mU32)) + (15*msq2*zt2)/(2.*(
         -m32 + mU32)) + (mU32*zt2)/(4.*(-m32 + mU32)) - (27275*mQ32*msq2)/(
         432.*pow2(-m32 + mQ32)) - (2825*mQ32*mU32)/(432.*pow2(-m32 + mQ32)) -
         (65*mQ32*msq2*zt2)/(6.*pow2(-m32 + mQ32)) - (mQ32*mU32*zt2)/pow2(-m32
         + mQ32) - (23*pow2(mQ32))/pow2(-m32 + mQ32) + (1355*msq2*pow2(mQ32))/(
         432.*(mQ32 - msq2)*pow2(-m32 + mQ32)) + (233*mU32*pow2(mQ32))/(432.*(
         mQ32 - mU32)*pow2(-m32 + mQ32)) + (5*msq2*zt2*pow2(mQ32))/(6.*(mQ32 -
         msq2)*pow2(-m32 + mQ32)) - (515*mQ32*pow2(msq2))/(144.*(mQ32 - msq2)*
         pow2(-m32 + mQ32)) - (5*mQ32*zt2*pow2(msq2))/(3.*(mQ32 - msq2)*pow2(
         -m32 + mQ32)) + (95*pow2(mQ32)*pow2(msq2))/(216.*pow2(-m32 + mQ32)*
         pow2(mQ32 - msq2)) + (5*zt2*pow2(mQ32)*pow2(msq2))/(6.*pow2(-m32 +
         mQ32)*pow2(mQ32 - msq2)) - (233*mQ32*pow2(mU32))/(288.*(mQ32 - mU32)*
         pow2(-m32 + mQ32)) + (233*pow2(mQ32)*pow2(mU32))/(864.*pow2(-m32 +
         mQ32)*pow2(mQ32 - mU32)) + 3*(16.09722222222222 - (22*mQ32)/(-m32 +
         mQ32) - (14*mQ32)/(-m32 + mU32) - (36*mU32)/(-m32 + mU32) + zt2 - (3*
         mQ32*zt2)/(2.*(-m32 + mQ32)) - (2*mQ32*zt2)/(-m32 + mU32) - (7*mU32*
         zt2)/(2.*(-m32 + mU32)) - ((mQ32 + mU32)*fin(mQ32,mU32,MR2))/((m32 -
         mQ32)*(m32 - mU32)) + (14*pow2(mQ32))/((-m32 + mQ32)*(-m32 + mU32)) +
         (2*zt2*pow2(mQ32))/((-m32 + mQ32)*(-m32 + mU32)) + (fin(mQ32,m32,MR2)*
         (3*m32*mQ32 - 3*mQ32*mU32 + pow2(m32) - pow2(mQ32)))/((m32 - mU32)*
         pow2(-m32 + mQ32)) + (31*pow2(mQ32))/pow2(-m32 + mQ32) + (3*zt2*pow2(
         mQ32))/pow2(-m32 + mQ32) + (fin(mU32,m32,MR2)*(3*m32*mU32 - 3*mQ32*
         mU32 + pow2(m32) - pow2(mU32)))/((m32 - mQ32)*pow2(m32 - mU32)) + (31*
         pow2(mU32))/pow2(-m32 + mU32) + (3*zt2*pow2(mU32))/pow2(-m32 + mU32))
         - (6*mQ32*mU32)/pow2(-m32 + mU32) - (27275*msq2*mU32)/(432.*pow2(-m32
         + mU32)) - (mQ32*mU32*zt2)/pow2(-m32 + mU32) - (65*msq2*mU32*zt2)/(6.*
         pow2(-m32 + mU32)) - (515*mU32*pow2(msq2))/(144.*(-msq2 + mU32)*pow2(
         -m32 + mU32)) - (5*mU32*zt2*pow2(msq2))/(3.*(-msq2 + mU32)*pow2(-m32 +
         mU32)) - (23*pow2(mU32))/pow2(-m32 + mU32) + (1355*msq2*pow2(mU32))/(
         432.*(-msq2 + mU32)*pow2(-m32 + mU32)) + (5*msq2*zt2*pow2(mU32))/(6.*(
         -msq2 + mU32)*pow2(-m32 + mU32)) + (95*pow2(msq2)*pow2(mU32))/(216.*
         pow2(-m32 + mU32)*pow2(-msq2 + mU32)) + (5*zt2*pow2(msq2)*pow2(mU32))/
         (6.*pow2(-m32 + mU32)*pow2(-msq2 + mU32)) - (5*fin(mQ32,msq2,MR2)*(-3*
         m32*mQ32 + 3*m32*msq2 + mQ32*msq2 - 2*pow2(m32) + pow2(mQ32)))/(2.*
         pow3(-m32 + mQ32)) + (fin(mU32,m32,MR2)*(6*m32*mQ32*mU32 - 8*mQ32*pow2
         (m32) + 2*mU32*pow2(m32) + 3*m32*pow2(mQ32) - 3*m32*pow2(mU32) - mQ32*
         pow2(mU32) + 2*pow3(m32) - pow3(mQ32)))/(4.*(m32 - mU32)*pow3(-m32 +
         mQ32)) - (95*mQ32*pow3(msq2))/(216.*pow2(-m32 + mQ32)*pow2(mQ32 - msq2
         )) - (5*mQ32*zt2*pow3(msq2))/(6.*pow2(-m32 + mQ32)*pow2(mQ32 - msq2))
         - (95*mU32*pow3(msq2))/(216.*pow2(-m32 + mU32)*pow2(-msq2 + mU32)) - (
         5*mU32*zt2*pow3(msq2))/(6.*pow2(-m32 + mU32)*pow2(-msq2 + mU32)) + (
         271*pow3(mU32))/(864.*(mQ32 - mU32)*pow2(-m32 + mQ32)) + (zt2*pow3(
         mU32))/(12.*(mQ32 - mU32)*pow2(-m32 + mQ32)) - (7*mQ32*pow3(mU32))/(
         12.*pow2(-m32 + mQ32)*pow2(mQ32 - mU32)) - (mQ32*zt2*pow3(mU32))/(12.*
         pow2(-m32 + mQ32)*pow2(mQ32 - mU32)) + (fin(mQ32,m32,MR2)*(-6*m32*mQ32
         *mU32 - 2*mQ32*pow2(m32) + 8*mU32*pow2(m32) + 3*m32*pow2(mQ32) + mU32*
         pow2(mQ32) - 3*m32*pow2(mU32) - 2*pow3(m32) + pow3(mU32)))/(4.*(m32 -
         mQ32)*pow3(m32 - mU32)) + fin(m32,msq2,MR2)*(5/(2.*(-m32 + mQ32)) + 5/
         (2.*(-m32 + mU32)) + (15*mQ32)/(2.*pow2(-m32 + mQ32)) - (15*msq2)/(2.*
         pow2(-m32 + mQ32)) - (15*msq2)/(2.*pow2(-m32 + mU32)) + (15*mU32)/(2.*
         pow2(-m32 + mU32)) + (10*mQ32*msq2)/pow3(-m32 + mQ32) - (10*pow2(mQ32)
         )/pow3(-m32 + mQ32) + (10*msq2*mU32)/pow3(-m32 + mU32) - (10*pow2(mU32
         ))/pow3(-m32 + mU32)) - (5*fin(mU32,msq2,MR2)*(3*m32*msq2 - 3*m32*mU32
         + msq2*mU32 - 2*pow2(m32) + pow2(mU32)))/(2.*pow3(-m32 + mU32)) + (
         271*pow4(mU32))/(864.*pow2(-m32 + mQ32)*pow2(mQ32 - mU32)) + (zt2*pow4
         (mU32))/(12.*pow2(-m32 + mQ32)*pow2(mQ32 - mU32)) + (fin(mQ32,mU32,MR2
         )*(15*mU32*pow2(m32)*pow2(mQ32) + 15*mQ32*pow2(m32)*pow2(mU32) - 6*m32
         *pow2(mQ32)*pow2(mU32) - 20*mQ32*mU32*pow3(m32) + 14*pow2(mQ32)*pow3(
         m32) + 14*pow2(mU32)*pow3(m32) - 6*m32*mU32*pow3(mQ32) - 11*pow2(m32)*
         pow3(mQ32) + pow2(mU32)*pow3(mQ32) - 6*m32*mQ32*pow3(mU32) - 11*pow2(
         m32)*pow3(mU32) + pow2(mQ32)*pow3(mU32) - 6*mQ32*pow4(m32) - 6*mU32*
         pow4(m32) + 3*m32*pow4(mQ32) + mU32*pow4(mQ32) + 3*m32*pow4(mU32) +
         mQ32*pow4(mU32) + 4*pow5(m32)))/(4.*pow3(-m32 + mQ32)*pow3(m32 - mU32)
         )))/3.)))/pow2(mQ32 - mU32)) - (48*((64*(m3*fin(mQ32,m32,MR2) - m3*fin
         (mQ32,mU32,MR2) + m3*fin(mU32,m32,MR2) + 3*pow3(m3) + zt2*pow3(m3)))/(
         9.*(m32 - mQ32)*(m32 - mU32)) - (4*((20*(m3*mQ32 - m3*msq2)*fin(mQ32,
         msq2,MR2))/((mQ32 - mU32)*pow2(m32 - mQ32)) - (2*m3*(mQ32 - mU32)*fin(
         mU32,m32,MR2))/((m32 - mU32)*pow2(m32 - mQ32)) + (2*m3*(mQ32 - mU32)*
         fin(mQ32,m32,MR2))/((m32 - mQ32)*pow2(m32 - mU32)) + (20*(-(m3*msq2) +
         m3*mU32)*fin(mU32,msq2,MR2))/((-mQ32 + mU32)*pow2(m32 - mU32)) + (2*(
         6*m3*mQ32 + 60*m3*msq2 + 6*m3*mU32 + m3*mQ32*zt2 + 10*m3*msq2*zt2 + m3
         *mU32*zt2 + 24*pow3(m3)))/((m32 - mQ32)*(m32 - mU32)) - (20*fin(m32,
         msq2,MR2)*(m3*mQ32*msq2 - 2*m3*mQ32*mU32 + m3*msq2*mU32 + mQ32*pow3(m3
         ) - 2*msq2*pow3(m3) + mU32*pow3(m3)))/(pow2(m32 - mQ32)*pow2(m32 -
         mU32)) + 3*((-2*m3*fin(mQ32,mU32,MR2))/((m32 - mQ32)*(m32 - mU32)) + (
         2*fin(mQ32,m32,MR2)*(m3*mQ32 - 3*m3*mU32 + 2*pow3(m3)))/((m32 - mQ32)*
         (m32 - mU32)*(mQ32 - mU32)) - (2*fin(mU32,m32,MR2)*(-3*m3*mQ32 + m3*
         mU32 + 2*pow3(m3)))/((m32 - mQ32)*(m32 - mU32)*(mQ32 - mU32)) + (2*(11
         *pow3(m3) + 3*zt2*pow3(m3)))/((m32 - mQ32)*(m32 - mU32))) + (2*fin(
         mQ32,mU32,MR2)*(m3*pow2(mQ32) + m3*pow2(mU32) - 2*mQ32*pow3(m3) - 2*
         mU32*pow3(m3) + 2*pow5(m3)))/(pow2(m32 - mQ32)*pow2(m32 - mU32))))/3.)
         *pow5(Xt))/pow2(mQ32 - mU32) + log(msq2/MR2)*((-96*m32*pow2(Xt))/(mQ32
         *mU32) + Xt*((2560*m3*msq2)/((m32 - mQ32)*(m32 - mU32)) + (384*(-(mU32
         *dilog(1 - m32/mQ32)) + m32*(dilog(1 - m32/mQ32) - dilog(1 - m32/mU32)
         ) + mQ32*dilog(1 - m32/mU32))*pow3(m3))/((m32 - mQ32)*(m32 - mU32)*(
         mQ32 - mU32))) + (192*m3*(2*m32*(dilog(1 - m32/mQ32) - dilog(1 -
         m32/mU32)) + mQ32*(-4 - dilog(1 - m32/mQ32) + dilog(1 - m32/mU32)) +
         mU32*(4 - dilog(1 - m32/mQ32) + dilog(1 - m32/mU32)))*pow3(Xt))/pow3(
         mQ32 - mU32) + ((48*((2*mQ32*mU32*(-dilog(1 - m32/mQ32) + dilog(1 -
         m32/mU32)) - pow2(mQ32) + pow2(mU32))*pow3(m32) + (-5 + dilog(1 -
         m32/mQ32) - dilog(1 - m32/mU32))*pow2(mU32)*pow3(mQ32) + pow2(m32)*(3*
         mU32*(-2 + dilog(1 - m32/mQ32) - dilog(1 - m32/mU32))*pow2(mQ32) + 3*
         mQ32*(2 + dilog(1 - m32/mQ32) - dilog(1 - m32/mU32))*pow2(mU32) + pow3
         (mQ32) - pow3(mU32)) + (5 + dilog(1 - m32/mQ32) - dilog(1 - m32/mU32))
         *pow2(mQ32)*pow3(mU32) + m32*(4*(-dilog(1 - m32/mQ32) + dilog(1 -
         m32/mU32))*pow2(mQ32)*pow2(mU32) + mU32*(5 - dilog(1 - m32/mQ32) +
         dilog(1 - m32/mU32))*pow3(mQ32) + mQ32*(-5 - dilog(1 - m32/mQ32) +
         dilog(1 - m32/mU32))*pow3(mU32))))/(mQ32*(-m32 + mQ32)*(m32 - mU32)*
         mU32*pow3(mQ32 - mU32)) + (16*(-480*m32*mQ32*msq2*mU32 - 600*mQ32*msq2
         *pow2(m32) + 480*mQ32*mU32*pow2(m32) - 600*msq2*mU32*pow2(m32) + 360*
         m32*msq2*pow2(mQ32) - 239*m32*mU32*pow2(mQ32) + 120*msq2*mU32*pow2(
         mQ32) + 120*pow2(m32)*pow2(mQ32) - 239*m32*mQ32*pow2(mU32) + 360*m32*
         msq2*pow2(mU32) + 120*mQ32*msq2*pow2(mU32) + 120*pow2(m32)*pow2(mU32)
         + 119*pow2(mQ32)*pow2(mU32) - 241*mQ32*pow3(m32) + 720*msq2*pow3(m32)
         - 241*mU32*pow3(m32) + 121*pow4(m32)))/(3.*pow2(m32 - mQ32)*pow2(m32 -
         mU32)*pow2(mQ32 - mU32)))*pow4(Xt) + (16*(300*msq2*mU32*pow2(m32)*
         pow2(mQ32) + 300*mQ32*msq2*pow2(m32)*pow2(mU32) + 240*m32*msq2*pow2(
         mQ32)*pow2(mU32) - 348*pow2(m32)*pow2(mQ32)*pow2(mU32) - 360*mQ32*msq2
         *mU32*pow3(m32) + 206*mU32*pow2(mQ32)*pow3(m32) + 36*mU32*dilog(1 -
         m32/mU32)*pow2(mQ32)*pow3(m32) + 206*mQ32*pow2(mU32)*pow3(m32) + 36*
         mQ32*dilog(1 - m32/mQ32)*pow2(mU32)*pow3(m32) - 180*m32*msq2*mU32*pow3
         (mQ32) - 96*mU32*pow2(m32)*pow3(mQ32) - 18*mU32*dilog(1 - m32/mU32)*
         pow2(m32)*pow3(mQ32) + 160*m32*pow2(mU32)*pow3(mQ32) - 60*msq2*pow2(
         mU32)*pow3(mQ32) + 9*pow3(m32)*pow3(mQ32) - 180*m32*mQ32*msq2*pow3(
         mU32) - 96*mQ32*pow2(m32)*pow3(mU32) - 18*mQ32*dilog(1 - m32/mQ32)*
         pow2(m32)*pow3(mU32) + 160*m32*pow2(mQ32)*pow3(mU32) - 60*msq2*pow2(
         mQ32)*pow3(mU32) + 9*pow3(m32)*pow3(mU32) - 73*pow3(mQ32)*pow3(mU32) -
         119*mQ32*mU32*pow4(m32) - 18*mQ32*mU32*dilog(1 - m32/mQ32)*pow4(m32)
         - 18*mQ32*mU32*dilog(1 - m32/mU32)*pow4(m32) - 18*pow2(mQ32)*pow4(m32)
         - 18*pow2(mU32)*pow4(m32) + 9*mQ32*pow5(m32) + 9*mU32*pow5(m32)))/(3.
         *mQ32*mU32*pow2(m32 - mQ32)*pow2(m32 - mU32)) - (5120*m3*msq2*pow5(Xt)
         )/((m32 - mQ32)*(m32 - mU32)*pow2(mQ32 - mU32))) + pow2(log(msq2/MR2))
         *((-320*m3*msq2*Xt)/((m32 - mQ32)*(m32 - mU32)) - (40*(4*m32*mQ32*msq2
         *mU32 + 5*mQ32*msq2*pow2(m32) + 8*mQ32*mU32*pow2(m32) + 5*msq2*mU32*
         pow2(m32) - 3*m32*msq2*pow2(mQ32) - 4*m32*mU32*pow2(mQ32) - msq2*mU32*
         pow2(mQ32) + 2*pow2(m32)*pow2(mQ32) - 4*m32*mQ32*pow2(mU32) - 3*m32*
         msq2*pow2(mU32) - mQ32*msq2*pow2(mU32) + 2*pow2(m32)*pow2(mU32) + 2*
         pow2(mQ32)*pow2(mU32) - 4*mQ32*pow3(m32) - 6*msq2*pow3(m32) - 4*mU32*
         pow3(m32) + 2*pow4(m32)))/(pow2(m32 - mQ32)*pow2(m32 - mU32)) + (80*(4
         *m32*mQ32*msq2*mU32 + 5*mQ32*msq2*pow2(m32) + 8*mQ32*mU32*pow2(m32) +
         5*msq2*mU32*pow2(m32) - 3*m32*msq2*pow2(mQ32) - 4*m32*mU32*pow2(mQ32)
         - msq2*mU32*pow2(mQ32) + 2*pow2(m32)*pow2(mQ32) - 4*m32*mQ32*pow2(mU32
         ) - 3*m32*msq2*pow2(mU32) - mQ32*msq2*pow2(mU32) + 2*pow2(m32)*pow2(
         mU32) + 2*pow2(mQ32)*pow2(mU32) - 4*mQ32*pow3(m32) - 6*msq2*pow3(m32)
         - 4*mU32*pow3(m32) + 2*pow4(m32))*pow4(Xt))/(pow2(m32 - mQ32)*pow2(m32
         - mU32)*pow2(mQ32 - mU32)) + (640*m3*msq2*pow5(Xt))/((m32 - mQ32)*(
         m32 - mU32)*pow2(mQ32 - mU32))) + log(mU32/MR2)*((-256*m32*mU32)/(3.*
         pow2(m32 - mU32)) + (128*pow2(mU32))/(3.*pow2(m32 - mU32)) + (8*(11*
         m32*mQ32*mU32 - 13*mU32*pow2(m32) + 6*m32*pow2(mU32) - 5*mQ32*pow2(
         mU32) + pow3(m32)))/(3.*(m32 - mQ32)*pow2(m32 - mU32)) - (96*(-(mQ32*
         mU32) + pow2(m32))*(2*m32*mU32 - pow2(mU32)))/((m32 - mQ32)*pow3(m32 -
         mU32)) + (16*(3*m32*mQ32 + 3*m32*mU32 + 5*mQ32*mU32 - 11*pow2(m32))*(
         -((6*mQ32 + 7*mU32)*pow2(m32)) - 2*mQ32*pow2(mU32) + m32*(5*mQ32*mU32
         + 3*pow2(mU32)) + 7*pow3(m32)))/(3.*pow2(m32 - mQ32)*pow3(m32 - mU32))
         + pow3(Xt)*((-2048*pow3(m3))/(3.*mQ32*(m32 - mU32)*(-mQ32 + mU32)) -
         (128*m3*(mQ32 + 3*mU32)*(3*m32*mQ32 + 3*m32*mU32 + 5*mQ32*mU32 - 11*
         pow2(m32)))/(3.*(m32 - mQ32)*(m32 - mU32)*pow3(mQ32 - mU32)) + (64*m3*
         (2*m32*(dilog(1 - m32/mQ32) - dilog(1 - m32/mU32)) + mQ32*(-4 - dilog(
         1 - m32/mQ32) + dilog(1 - m32/mU32)) + mU32*(4 - dilog(1 - m32/mQ32) +
         dilog(1 - m32/mU32)))*(-34*m32*mU32 + pow2(m32) + 17*pow2(mU32)))/(3.
         *pow2(m32 - mU32)*pow3(mQ32 - mU32)) + (48*((64*(m3*fin(mQ32,m32,MR2)
         - m3*fin(mQ32,mU32,MR2) + m3*fin(mU32,m32,MR2) + 3*pow3(m3) + zt2*pow3
         (m3)))/(9.*(m32 - mQ32)*(m32 - mU32)) - (4*((20*(m3*mQ32 - m3*msq2)*
         fin(mQ32,msq2,MR2))/((mQ32 - mU32)*pow2(m32 - mQ32)) - (2*m3*(mQ32 -
         mU32)*fin(mU32,m32,MR2))/((m32 - mU32)*pow2(m32 - mQ32)) + (2*m3*(mQ32
         - mU32)*fin(mQ32,m32,MR2))/((m32 - mQ32)*pow2(m32 - mU32)) + (20*(-(
         m3*msq2) + m3*mU32)*fin(mU32,msq2,MR2))/((-mQ32 + mU32)*pow2(m32 -
         mU32)) + (2*(6*m3*mQ32 + 60*m3*msq2 + 6*m3*mU32 + m3*mQ32*zt2 + 10*m3*
         msq2*zt2 + m3*mU32*zt2 + 24*pow3(m3)))/((m32 - mQ32)*(m32 - mU32)) - (
         20*fin(m32,msq2,MR2)*(m3*mQ32*msq2 - 2*m3*mQ32*mU32 + m3*msq2*mU32 +
         mQ32*pow3(m3) - 2*msq2*pow3(m3) + mU32*pow3(m3)))/(pow2(m32 - mQ32)*
         pow2(m32 - mU32)) + 3*((-2*m3*fin(mQ32,mU32,MR2))/((m32 - mQ32)*(m32 -
         mU32)) + (2*fin(mQ32,m32,MR2)*(m3*mQ32 - 3*m3*mU32 + 2*pow3(m3)))/((
         m32 - mQ32)*(m32 - mU32)*(mQ32 - mU32)) - (2*fin(mU32,m32,MR2)*(-3*m3*
         mQ32 + m3*mU32 + 2*pow3(m3)))/((m32 - mQ32)*(m32 - mU32)*(mQ32 - mU32)
         ) + (2*(11*pow3(m3) + 3*zt2*pow3(m3)))/((m32 - mQ32)*(m32 - mU32))) +
         (2*fin(mQ32,mU32,MR2)*(m3*pow2(mQ32) + m3*pow2(mU32) - 2*mQ32*pow3(m3)
         - 2*mU32*pow3(m3) + 2*pow5(m3)))/(pow2(m32 - mQ32)*pow2(m32 - mU32)))
         )/3.))/(mQ32 - mU32)) + (8*(-34*m32*mU32 + pow2(m32) + 17*pow2(mU32))*
         (9*m32*(mQ32 + mU32)*pow2(mQ32)*pow2(mU32) - 3*pow3(mQ32)*pow3(mU32) +
         pow3(m32)*(mU32*(19 + 8*dilog(1 - m32/mU32))*pow2(mQ32) + mQ32*(19 +
         8*dilog(1 - m32/mQ32))*pow2(mU32) + 2*pow3(mQ32) + 2*pow3(mU32)) - 4*
         pow2(m32)*(6*pow2(mQ32)*pow2(mU32) + mU32*(2 + dilog(1 - m32/mU32))*
         pow3(mQ32) + mQ32*(2 + dilog(1 - m32/mQ32))*pow3(mU32)) - (mQ32*mU32*(
         13 + 4*dilog(1 - m32/mQ32) + 4*dilog(1 - m32/mU32)) + 4*pow2(mQ32) + 4
         *pow2(mU32))*pow4(m32) + 2*(mQ32 + mU32)*pow5(m32)))/(3.*mQ32*mU32*
         pow2(m32 - mQ32)*pow4(m32 - mU32)) + 6*((8*(-(mQ32*mU32) + pow2(m32)))
         /(3.*(m32 - mQ32)*(m32 - mU32)) + (40*pow2(-(mQ32*mU32) + pow2(m32)))/
         (9.*pow2(m32 - mQ32)*pow2(m32 - mU32)) - 4*(-4.111111111111111 - (2*
         mQ32)/(3.*(-m32 + mQ32)) - (2*mU32)/(3.*(-m32 + mU32)) + pow2(-2 + (2*
         mQ32)/(3.*(-m32 + mQ32)) + (2*mU32)/(3.*(-m32 + mU32))) - (16*(-29.625
         - (21*mQ32)/(4.*(-m32 + mQ32)) - (56*mQ32)/(mQ32 - mU32) + (56*mU32)/
         (mQ32 - mU32) + (111*mQ32)/(4.*(-m32 + mU32)) + (45*mU32)/(2.*(-m32 +
         mU32)) + zt2/2. - (5*mQ32*zt2)/(2.*(-m32 + mQ32)) - (8*mQ32*zt2)/(mQ32
         - mU32) + (8*mU32*zt2)/(mQ32 - mU32) + (4*mQ32*zt2)/(-m32 + mU32) + (
         3*mU32*zt2)/(2.*(-m32 + mU32)) + (2*(mQ32 + mU32)*fin(mQ32,mU32,MR2))/
         ((m32 - mQ32)*(m32 - mU32)) + (56*pow2(mQ32))/((-m32 + mQ32)*(mQ32 -
         mU32)) - (111*pow2(mQ32))/(4.*(-m32 + mQ32)*(-m32 + mU32)) + (8*zt2*
         pow2(mQ32))/((-m32 + mQ32)*(mQ32 - mU32)) - (4*zt2*pow2(mQ32))/((-m32
         + mQ32)*(-m32 + mU32)) + (23*pow2(mQ32))/(8.*pow2(-m32 + mQ32)) + (3*
         zt2*pow2(mQ32))/(2.*pow2(-m32 + mQ32)) - (56*pow2(mU32))/((mQ32 - mU32
         )*(-m32 + mU32)) - (8*zt2*pow2(mU32))/((mQ32 - mU32)*(-m32 + mU32)) +
         (23*pow2(mU32))/(8.*pow2(-m32 + mU32)) + (3*zt2*pow2(mU32))/(2.*pow2(
         -m32 + mU32)) + fin(mQ32,m32,MR2)*(-5/(2.*(-m32 + mQ32)) - 8/(mQ32 -
         mU32) + (8*mQ32)/((-m32 + mQ32)*(mQ32 - mU32)) + 2/(-m32 + mU32) - (4*
         mQ32)/((-m32 + mQ32)*(-m32 + mU32)) - (3*pow2(mQ32))/pow3(-m32 + mQ32)
         ) - (21*pow3(mQ32))/pow3(-m32 + mQ32) - (3*zt2*pow3(mQ32))/pow3(-m32 +
         mQ32) + fin(mU32,m32,MR2)*(-2/(-m32 + mQ32) + 8/(mQ32 - mU32) + 3/(2.
         *(-m32 + mU32)) - (4*mQ32)/((-m32 + mQ32)*(-m32 + mU32)) - (8*mU32)/((
         mQ32 - mU32)*(-m32 + mU32)) - (3*pow2(mU32))/pow3(-m32 + mU32)) - (21*
         pow3(mU32))/pow3(-m32 + mU32) - (3*zt2*pow3(mU32))/pow3(-m32 + mU32)))
         /9. - (4*(-16.333333333333332 + (31*mQ32)/(-m32 + mQ32) + (45*msq2)/(
         -m32 + mQ32) + (9*mU32)/(2.*(-m32 + mQ32)) + (9*mQ32)/(2.*(-m32 + mU32
         )) + (45*msq2)/(-m32 + mU32) + (31*mU32)/(-m32 + mU32) - 3*zt2 + (mQ32
         *zt2)/(4.*(-m32 + mQ32)) + (15*msq2*zt2)/(2.*(-m32 + mQ32)) + (3*mU32*
         zt2)/(4.*(-m32 + mQ32)) + (3*mQ32*zt2)/(4.*(-m32 + mU32)) + (15*msq2*
         zt2)/(2.*(-m32 + mU32)) + (mU32*zt2)/(4.*(-m32 + mU32)) - (27275*mQ32*
         msq2)/(432.*pow2(-m32 + mQ32)) - (2825*mQ32*mU32)/(432.*pow2(-m32 +
         mQ32)) - (65*mQ32*msq2*zt2)/(6.*pow2(-m32 + mQ32)) - (mQ32*mU32*zt2)
         /pow2(-m32 + mQ32) - (23*pow2(mQ32))/pow2(-m32 + mQ32) + (1355*msq2*
         pow2(mQ32))/(432.*(mQ32 - msq2)*pow2(-m32 + mQ32)) + (233*mU32*pow2(
         mQ32))/(432.*(mQ32 - mU32)*pow2(-m32 + mQ32)) + (5*msq2*zt2*pow2(mQ32)
         )/(6.*(mQ32 - msq2)*pow2(-m32 + mQ32)) - (515*mQ32*pow2(msq2))/(144.*(
         mQ32 - msq2)*pow2(-m32 + mQ32)) - (5*mQ32*zt2*pow2(msq2))/(3.*(mQ32 -
         msq2)*pow2(-m32 + mQ32)) + (95*pow2(mQ32)*pow2(msq2))/(216.*pow2(-m32
         + mQ32)*pow2(mQ32 - msq2)) + (5*zt2*pow2(mQ32)*pow2(msq2))/(6.*pow2(
         -m32 + mQ32)*pow2(mQ32 - msq2)) - (233*mQ32*pow2(mU32))/(288.*(mQ32 -
         mU32)*pow2(-m32 + mQ32)) + (233*pow2(mQ32)*pow2(mU32))/(864.*pow2(-m32
         + mQ32)*pow2(mQ32 - mU32)) + 3*(16.09722222222222 - (22*mQ32)/(-m32 +
         mQ32) - (14*mQ32)/(-m32 + mU32) - (36*mU32)/(-m32 + mU32) + zt2 - (3*
         mQ32*zt2)/(2.*(-m32 + mQ32)) - (2*mQ32*zt2)/(-m32 + mU32) - (7*mU32*
         zt2)/(2.*(-m32 + mU32)) - ((mQ32 + mU32)*fin(mQ32,mU32,MR2))/((m32 -
         mQ32)*(m32 - mU32)) + (14*pow2(mQ32))/((-m32 + mQ32)*(-m32 + mU32)) +
         (2*zt2*pow2(mQ32))/((-m32 + mQ32)*(-m32 + mU32)) + (fin(mQ32,m32,MR2)*
         (3*m32*mQ32 - 3*mQ32*mU32 + pow2(m32) - pow2(mQ32)))/((m32 - mU32)*
         pow2(-m32 + mQ32)) + (31*pow2(mQ32))/pow2(-m32 + mQ32) + (3*zt2*pow2(
         mQ32))/pow2(-m32 + mQ32) + (fin(mU32,m32,MR2)*(3*m32*mU32 - 3*mQ32*
         mU32 + pow2(m32) - pow2(mU32)))/((m32 - mQ32)*pow2(m32 - mU32)) + (31*
         pow2(mU32))/pow2(-m32 + mU32) + (3*zt2*pow2(mU32))/pow2(-m32 + mU32))
         - (6*mQ32*mU32)/pow2(-m32 + mU32) - (27275*msq2*mU32)/(432.*pow2(-m32
         + mU32)) - (mQ32*mU32*zt2)/pow2(-m32 + mU32) - (65*msq2*mU32*zt2)/(6.*
         pow2(-m32 + mU32)) - (515*mU32*pow2(msq2))/(144.*(-msq2 + mU32)*pow2(
         -m32 + mU32)) - (5*mU32*zt2*pow2(msq2))/(3.*(-msq2 + mU32)*pow2(-m32 +
         mU32)) - (23*pow2(mU32))/pow2(-m32 + mU32) + (1355*msq2*pow2(mU32))/(
         432.*(-msq2 + mU32)*pow2(-m32 + mU32)) + (5*msq2*zt2*pow2(mU32))/(6.*(
         -msq2 + mU32)*pow2(-m32 + mU32)) + (95*pow2(msq2)*pow2(mU32))/(216.*
         pow2(-m32 + mU32)*pow2(-msq2 + mU32)) + (5*zt2*pow2(msq2)*pow2(mU32))/
         (6.*pow2(-m32 + mU32)*pow2(-msq2 + mU32)) - (5*fin(mQ32,msq2,MR2)*(-3*
         m32*mQ32 + 3*m32*msq2 + mQ32*msq2 - 2*pow2(m32) + pow2(mQ32)))/(2.*
         pow3(-m32 + mQ32)) + (fin(mU32,m32,MR2)*(6*m32*mQ32*mU32 - 8*mQ32*pow2
         (m32) + 2*mU32*pow2(m32) + 3*m32*pow2(mQ32) - 3*m32*pow2(mU32) - mQ32*
         pow2(mU32) + 2*pow3(m32) - pow3(mQ32)))/(4.*(m32 - mU32)*pow3(-m32 +
         mQ32)) - (95*mQ32*pow3(msq2))/(216.*pow2(-m32 + mQ32)*pow2(mQ32 - msq2
         )) - (5*mQ32*zt2*pow3(msq2))/(6.*pow2(-m32 + mQ32)*pow2(mQ32 - msq2))
         - (95*mU32*pow3(msq2))/(216.*pow2(-m32 + mU32)*pow2(-msq2 + mU32)) - (
         5*mU32*zt2*pow3(msq2))/(6.*pow2(-m32 + mU32)*pow2(-msq2 + mU32)) + (
         271*pow3(mU32))/(864.*(mQ32 - mU32)*pow2(-m32 + mQ32)) + (zt2*pow3(
         mU32))/(12.*(mQ32 - mU32)*pow2(-m32 + mQ32)) - (7*mQ32*pow3(mU32))/(
         12.*pow2(-m32 + mQ32)*pow2(mQ32 - mU32)) - (mQ32*zt2*pow3(mU32))/(12.*
         pow2(-m32 + mQ32)*pow2(mQ32 - mU32)) + (fin(mQ32,m32,MR2)*(-6*m32*mQ32
         *mU32 - 2*mQ32*pow2(m32) + 8*mU32*pow2(m32) + 3*m32*pow2(mQ32) + mU32*
         pow2(mQ32) - 3*m32*pow2(mU32) - 2*pow3(m32) + pow3(mU32)))/(4.*(m32 -
         mQ32)*pow3(m32 - mU32)) + fin(m32,msq2,MR2)*(5/(2.*(-m32 + mQ32)) + 5/
         (2.*(-m32 + mU32)) + (15*mQ32)/(2.*pow2(-m32 + mQ32)) - (15*msq2)/(2.*
         pow2(-m32 + mQ32)) - (15*msq2)/(2.*pow2(-m32 + mU32)) + (15*mU32)/(2.*
         pow2(-m32 + mU32)) + (10*mQ32*msq2)/pow3(-m32 + mQ32) - (10*pow2(mQ32)
         )/pow3(-m32 + mQ32) + (10*msq2*mU32)/pow3(-m32 + mU32) - (10*pow2(mU32
         ))/pow3(-m32 + mU32)) - (5*fin(mU32,msq2,MR2)*(3*m32*msq2 - 3*m32*mU32
         + msq2*mU32 - 2*pow2(m32) + pow2(mU32)))/(2.*pow3(-m32 + mU32)) + (
         271*pow4(mU32))/(864.*pow2(-m32 + mQ32)*pow2(mQ32 - mU32)) + (zt2*pow4
         (mU32))/(12.*pow2(-m32 + mQ32)*pow2(mQ32 - mU32)) + (fin(mQ32,mU32,MR2
         )*(15*mU32*pow2(m32)*pow2(mQ32) + 15*mQ32*pow2(m32)*pow2(mU32) - 6*m32
         *pow2(mQ32)*pow2(mU32) - 20*mQ32*mU32*pow3(m32) + 14*pow2(mQ32)*pow3(
         m32) + 14*pow2(mU32)*pow3(m32) - 6*m32*mU32*pow3(mQ32) - 11*pow2(m32)*
         pow3(mQ32) + pow2(mU32)*pow3(mQ32) - 6*m32*mQ32*pow3(mU32) - 11*pow2(
         m32)*pow3(mU32) + pow2(mQ32)*pow3(mU32) - 6*mQ32*pow4(m32) - 6*mU32*
         pow4(m32) + 3*m32*pow4(mQ32) + mU32*pow4(mQ32) + 3*m32*pow4(mU32) +
         mQ32*pow4(mU32) + 4*pow5(m32)))/(4.*pow3(-m32 + mQ32)*pow3(m32 - mU32)
         )))/3.)) + pow2(Xt)*((8192*mU32*(-(mU32*dilog(1 - m32/mQ32)) + m32*(
         dilog(1 - m32/mQ32) - dilog(1 - m32/mU32)) + mQ32*dilog(1 - m32/mU32))
         *pow2(m32))/(3.*(m32 - mQ32)*(mQ32 - mU32)*(-mQ32 + mU32)*pow2(m32 -
         mU32)) - (32*(3*m32*mQ32 + 3*m32*mU32 + 5*mQ32*mU32 - 11*pow2(m32))*(3
         *mQ32*mU32 - 2*m32*(3*mQ32 + 2*mU32) + 7*pow2(m32)))/(3.*(mQ32 - mU32)
         *pow2(m32 - mQ32)*pow2(m32 - mU32)) - (32*m32*(-34*m32*mU32 + pow2(m32
         ) + 17*pow2(mU32)))/(3.*mQ32*mU32*pow2(m32 - mU32)) - (12*((8*(-(mQ32*
         mU32) + pow2(m32)))/(3.*(m32 - mQ32)*(m32 - mU32)) + (40*pow2(-(mQ32*
         mU32) + pow2(m32)))/(9.*pow2(m32 - mQ32)*pow2(m32 - mU32)) - 4*(
         -4.111111111111111 - (2*mQ32)/(3.*(-m32 + mQ32)) - (2*mU32)/(3.*(-m32
         + mU32)) + pow2(-2 + (2*mQ32)/(3.*(-m32 + mQ32)) + (2*mU32)/(3.*(-m32
         + mU32))) - (16*(-29.625 - (21*mQ32)/(4.*(-m32 + mQ32)) - (56*mQ32)/(
         mQ32 - mU32) + (56*mU32)/(mQ32 - mU32) + (111*mQ32)/(4.*(-m32 + mU32))
         + (45*mU32)/(2.*(-m32 + mU32)) + zt2/2. - (5*mQ32*zt2)/(2.*(-m32 +
         mQ32)) - (8*mQ32*zt2)/(mQ32 - mU32) + (8*mU32*zt2)/(mQ32 - mU32) + (4*
         mQ32*zt2)/(-m32 + mU32) + (3*mU32*zt2)/(2.*(-m32 + mU32)) + (2*(mQ32 +
         mU32)*fin(mQ32,mU32,MR2))/((m32 - mQ32)*(m32 - mU32)) + (56*pow2(mQ32
         ))/((-m32 + mQ32)*(mQ32 - mU32)) - (111*pow2(mQ32))/(4.*(-m32 + mQ32)*
         (-m32 + mU32)) + (8*zt2*pow2(mQ32))/((-m32 + mQ32)*(mQ32 - mU32)) - (4
         *zt2*pow2(mQ32))/((-m32 + mQ32)*(-m32 + mU32)) + (23*pow2(mQ32))/(8.*
         pow2(-m32 + mQ32)) + (3*zt2*pow2(mQ32))/(2.*pow2(-m32 + mQ32)) - (56*
         pow2(mU32))/((mQ32 - mU32)*(-m32 + mU32)) - (8*zt2*pow2(mU32))/((mQ32
         - mU32)*(-m32 + mU32)) + (23*pow2(mU32))/(8.*pow2(-m32 + mU32)) + (3*
         zt2*pow2(mU32))/(2.*pow2(-m32 + mU32)) + fin(mQ32,m32,MR2)*(-5/(2.*(
         -m32 + mQ32)) - 8/(mQ32 - mU32) + (8*mQ32)/((-m32 + mQ32)*(mQ32 - mU32
         )) + 2/(-m32 + mU32) - (4*mQ32)/((-m32 + mQ32)*(-m32 + mU32)) - (3*
         pow2(mQ32))/pow3(-m32 + mQ32)) - (21*pow3(mQ32))/pow3(-m32 + mQ32) - (
         3*zt2*pow3(mQ32))/pow3(-m32 + mQ32) + fin(mU32,m32,MR2)*(-2/(-m32 +
         mQ32) + 8/(mQ32 - mU32) + 3/(2.*(-m32 + mU32)) - (4*mQ32)/((-m32 +
         mQ32)*(-m32 + mU32)) - (8*mU32)/((mQ32 - mU32)*(-m32 + mU32)) - (3*
         pow2(mU32))/pow3(-m32 + mU32)) - (21*pow3(mU32))/pow3(-m32 + mU32) - (
         3*zt2*pow3(mU32))/pow3(-m32 + mU32)))/9. - (4*(-16.333333333333332 + (
         31*mQ32)/(-m32 + mQ32) + (45*msq2)/(-m32 + mQ32) + (9*mU32)/(2.*(-m32
         + mQ32)) + (9*mQ32)/(2.*(-m32 + mU32)) + (45*msq2)/(-m32 + mU32) + (31
         *mU32)/(-m32 + mU32) - 3*zt2 + (mQ32*zt2)/(4.*(-m32 + mQ32)) + (15*
         msq2*zt2)/(2.*(-m32 + mQ32)) + (3*mU32*zt2)/(4.*(-m32 + mQ32)) + (3*
         mQ32*zt2)/(4.*(-m32 + mU32)) + (15*msq2*zt2)/(2.*(-m32 + mU32)) + (
         mU32*zt2)/(4.*(-m32 + mU32)) - (27275*mQ32*msq2)/(432.*pow2(-m32 +
         mQ32)) - (2825*mQ32*mU32)/(432.*pow2(-m32 + mQ32)) - (65*mQ32*msq2*zt2
         )/(6.*pow2(-m32 + mQ32)) - (mQ32*mU32*zt2)/pow2(-m32 + mQ32) - (23*
         pow2(mQ32))/pow2(-m32 + mQ32) + (1355*msq2*pow2(mQ32))/(432.*(mQ32 -
         msq2)*pow2(-m32 + mQ32)) + (233*mU32*pow2(mQ32))/(432.*(mQ32 - mU32)*
         pow2(-m32 + mQ32)) + (5*msq2*zt2*pow2(mQ32))/(6.*(mQ32 - msq2)*pow2(
         -m32 + mQ32)) - (515*mQ32*pow2(msq2))/(144.*(mQ32 - msq2)*pow2(-m32 +
         mQ32)) - (5*mQ32*zt2*pow2(msq2))/(3.*(mQ32 - msq2)*pow2(-m32 + mQ32))
         + (95*pow2(mQ32)*pow2(msq2))/(216.*pow2(-m32 + mQ32)*pow2(mQ32 - msq2)
         ) + (5*zt2*pow2(mQ32)*pow2(msq2))/(6.*pow2(-m32 + mQ32)*pow2(mQ32 -
         msq2)) - (233*mQ32*pow2(mU32))/(288.*(mQ32 - mU32)*pow2(-m32 + mQ32))
         + (233*pow2(mQ32)*pow2(mU32))/(864.*pow2(-m32 + mQ32)*pow2(mQ32 - mU32
         )) + 3*(16.09722222222222 - (22*mQ32)/(-m32 + mQ32) - (14*mQ32)/(-m32
         + mU32) - (36*mU32)/(-m32 + mU32) + zt2 - (3*mQ32*zt2)/(2.*(-m32 +
         mQ32)) - (2*mQ32*zt2)/(-m32 + mU32) - (7*mU32*zt2)/(2.*(-m32 + mU32))
         - ((mQ32 + mU32)*fin(mQ32,mU32,MR2))/((m32 - mQ32)*(m32 - mU32)) + (14
         *pow2(mQ32))/((-m32 + mQ32)*(-m32 + mU32)) + (2*zt2*pow2(mQ32))/((-m32
         + mQ32)*(-m32 + mU32)) + (fin(mQ32,m32,MR2)*(3*m32*mQ32 - 3*mQ32*mU32
         + pow2(m32) - pow2(mQ32)))/((m32 - mU32)*pow2(-m32 + mQ32)) + (31*
         pow2(mQ32))/pow2(-m32 + mQ32) + (3*zt2*pow2(mQ32))/pow2(-m32 + mQ32) +
         (fin(mU32,m32,MR2)*(3*m32*mU32 - 3*mQ32*mU32 + pow2(m32) - pow2(mU32)
         ))/((m32 - mQ32)*pow2(m32 - mU32)) + (31*pow2(mU32))/pow2(-m32 + mU32)
         + (3*zt2*pow2(mU32))/pow2(-m32 + mU32)) - (6*mQ32*mU32)/pow2(-m32 +
         mU32) - (27275*msq2*mU32)/(432.*pow2(-m32 + mU32)) - (mQ32*mU32*zt2)
         /pow2(-m32 + mU32) - (65*msq2*mU32*zt2)/(6.*pow2(-m32 + mU32)) - (515*
         mU32*pow2(msq2))/(144.*(-msq2 + mU32)*pow2(-m32 + mU32)) - (5*mU32*zt2
         *pow2(msq2))/(3.*(-msq2 + mU32)*pow2(-m32 + mU32)) - (23*pow2(mU32))
         /pow2(-m32 + mU32) + (1355*msq2*pow2(mU32))/(432.*(-msq2 + mU32)*pow2(
         -m32 + mU32)) + (5*msq2*zt2*pow2(mU32))/(6.*(-msq2 + mU32)*pow2(-m32 +
         mU32)) + (95*pow2(msq2)*pow2(mU32))/(216.*pow2(-m32 + mU32)*pow2(
         -msq2 + mU32)) + (5*zt2*pow2(msq2)*pow2(mU32))/(6.*pow2(-m32 + mU32)*
         pow2(-msq2 + mU32)) - (5*fin(mQ32,msq2,MR2)*(-3*m32*mQ32 + 3*m32*msq2
         + mQ32*msq2 - 2*pow2(m32) + pow2(mQ32)))/(2.*pow3(-m32 + mQ32)) + (fin
         (mU32,m32,MR2)*(6*m32*mQ32*mU32 - 8*mQ32*pow2(m32) + 2*mU32*pow2(m32)
         + 3*m32*pow2(mQ32) - 3*m32*pow2(mU32) - mQ32*pow2(mU32) + 2*pow3(m32)
         - pow3(mQ32)))/(4.*(m32 - mU32)*pow3(-m32 + mQ32)) - (95*mQ32*pow3(
         msq2))/(216.*pow2(-m32 + mQ32)*pow2(mQ32 - msq2)) - (5*mQ32*zt2*pow3(
         msq2))/(6.*pow2(-m32 + mQ32)*pow2(mQ32 - msq2)) - (95*mU32*pow3(msq2))
         /(216.*pow2(-m32 + mU32)*pow2(-msq2 + mU32)) - (5*mU32*zt2*pow3(msq2))
         /(6.*pow2(-m32 + mU32)*pow2(-msq2 + mU32)) + (271*pow3(mU32))/(864.*(
         mQ32 - mU32)*pow2(-m32 + mQ32)) + (zt2*pow3(mU32))/(12.*(mQ32 - mU32)*
         pow2(-m32 + mQ32)) - (7*mQ32*pow3(mU32))/(12.*pow2(-m32 + mQ32)*pow2(
         mQ32 - mU32)) - (mQ32*zt2*pow3(mU32))/(12.*pow2(-m32 + mQ32)*pow2(mQ32
         - mU32)) + (fin(mQ32,m32,MR2)*(-6*m32*mQ32*mU32 - 2*mQ32*pow2(m32) +
         8*mU32*pow2(m32) + 3*m32*pow2(mQ32) + mU32*pow2(mQ32) - 3*m32*pow2(
         mU32) - 2*pow3(m32) + pow3(mU32)))/(4.*(m32 - mQ32)*pow3(m32 - mU32))
         + fin(m32,msq2,MR2)*(5/(2.*(-m32 + mQ32)) + 5/(2.*(-m32 + mU32)) + (15
         *mQ32)/(2.*pow2(-m32 + mQ32)) - (15*msq2)/(2.*pow2(-m32 + mQ32)) - (15
         *msq2)/(2.*pow2(-m32 + mU32)) + (15*mU32)/(2.*pow2(-m32 + mU32)) + (10
         *mQ32*msq2)/pow3(-m32 + mQ32) - (10*pow2(mQ32))/pow3(-m32 + mQ32) + (
         10*msq2*mU32)/pow3(-m32 + mU32) - (10*pow2(mU32))/pow3(-m32 + mU32)) -
         (5*fin(mU32,msq2,MR2)*(3*m32*msq2 - 3*m32*mU32 + msq2*mU32 - 2*pow2(
         m32) + pow2(mU32)))/(2.*pow3(-m32 + mU32)) + (271*pow4(mU32))/(864.*
         pow2(-m32 + mQ32)*pow2(mQ32 - mU32)) + (zt2*pow4(mU32))/(12.*pow2(-m32
         + mQ32)*pow2(mQ32 - mU32)) + (fin(mQ32,mU32,MR2)*(15*mU32*pow2(m32)*
         pow2(mQ32) + 15*mQ32*pow2(m32)*pow2(mU32) - 6*m32*pow2(mQ32)*pow2(mU32
         ) - 20*mQ32*mU32*pow3(m32) + 14*pow2(mQ32)*pow3(m32) + 14*pow2(mU32)*
         pow3(m32) - 6*m32*mU32*pow3(mQ32) - 11*pow2(m32)*pow3(mQ32) + pow2(
         mU32)*pow3(mQ32) - 6*m32*mQ32*pow3(mU32) - 11*pow2(m32)*pow3(mU32) +
         pow2(mQ32)*pow3(mU32) - 6*mQ32*pow4(m32) - 6*mU32*pow4(m32) + 3*m32*
         pow4(mQ32) + mU32*pow4(mQ32) + 3*m32*pow4(mU32) + mQ32*pow4(mU32) + 4*
         pow5(m32)))/(4.*pow3(-m32 + mQ32)*pow3(m32 - mU32))))/3.)))/(mQ32 -
         mU32)) + (8*(-360*msq2*mU32*pow2(m32)*pow2(mQ32) + 420*mQ32*msq2*pow2(
         m32)*pow2(mU32) - 300*m32*msq2*pow2(mQ32)*pow2(mU32) + 777*pow2(m32)*
         pow2(mQ32)*pow2(mU32) + 180*mQ32*msq2*mU32*pow3(m32) - 70*mU32*pow2(
         mQ32)*pow3(m32) + 245*mQ32*pow2(mU32)*pow3(m32) - 180*msq2*pow2(mU32)*
         pow3(m32) + 180*m32*msq2*mU32*pow3(mQ32) + 34*mU32*pow2(m32)*pow3(mQ32
         ) - 535*m32*pow2(mU32)*pow3(mQ32) + 60*msq2*pow2(mU32)*pow3(mQ32) + 53
         *pow3(m32)*pow3(mQ32) + 120*m32*mQ32*msq2*pow3(mU32) - 1256*mQ32*pow2(
         m32)*pow3(mU32) - 60*msq2*pow2(m32)*pow3(mU32) + 264*m32*pow2(mQ32)*
         pow3(mU32) - 60*msq2*pow2(mQ32)*pow3(mU32) + 540*pow3(m32)*pow3(mU32)
         + 164*pow3(mQ32)*pow3(mU32) - 106*mQ32*mU32*pow4(m32) - 105*pow2(mQ32)
         *pow4(m32) - 301*pow2(mU32)*pow4(m32) + 18*m32*mU32*pow4(mQ32) + 6*
         pow2(mU32)*pow4(mQ32) + 417*m32*mQ32*pow4(mU32) - 67*pow2(m32)*pow4(
         mU32) - 158*pow2(mQ32)*pow4(mU32) + 52*mQ32*pow5(m32) + 76*mU32*pow5(
         m32) - 36*m32*pow5(mU32) - 12*mQ32*pow5(mU32)))/(3.*(mQ32 - mU32)*pow2
         (-m32 + mQ32)*pow3(m32 - mU32)) + pow4(Xt)*((4096*m32*mU32*(2*m32*(
         dilog(1 - m32/mQ32) - dilog(1 - m32/mU32)) + mQ32*(-4 - dilog(1 -
         m32/mQ32) + dilog(1 - m32/mU32)) + mU32*(4 - dilog(1 - m32/mQ32) +
         dilog(1 - m32/mU32))))/(3.*(m32 - mU32)*(-mQ32 + mU32)*pow3(mQ32 -
         mU32)) + (16*(-34*m32*mU32 + pow2(m32) + 17*pow2(mU32))*((2*mQ32*mU32*
         (-dilog(1 - m32/mQ32) + dilog(1 - m32/mU32)) - pow2(mQ32) + pow2(mU32)
         )*pow3(m32) + (-5 + dilog(1 - m32/mQ32) - dilog(1 - m32/mU32))*pow2(
         mU32)*pow3(mQ32) + pow2(m32)*(3*mU32*(-2 + dilog(1 - m32/mQ32) - dilog
         (1 - m32/mU32))*pow2(mQ32) + 3*mQ32*(2 + dilog(1 - m32/mQ32) - dilog(1
         - m32/mU32))*pow2(mU32) + pow3(mQ32) - pow3(mU32)) + (5 + dilog(1 -
         m32/mQ32) - dilog(1 - m32/mU32))*pow2(mQ32)*pow3(mU32) + m32*(4*(
         -dilog(1 - m32/mQ32) + dilog(1 - m32/mU32))*pow2(mQ32)*pow2(mU32) +
         mU32*(5 - dilog(1 - m32/mQ32) + dilog(1 - m32/mU32))*pow3(mQ32) + mQ32
         *(-5 - dilog(1 - m32/mQ32) + dilog(1 - m32/mU32))*pow3(mU32))))/(3.*
         mQ32*(-m32 + mQ32)*mU32*pow3(m32 - mU32)*pow3(mQ32 - mU32)) + (16*(3*
         m32*mQ32 + 3*m32*mU32 + 5*mQ32*mU32 - 11*pow2(m32))*(-3*mQ32*(mQ32 + 5
         *mU32)*pow2(mU32) - pow2(m32)*(15*mQ32*mU32 + 6*pow2(mQ32) + 29*pow2(
         mU32)) + (3*mQ32 + 7*mU32)*pow3(m32) + m32*(7*mU32*pow2(mQ32) + 31*
         mQ32*pow2(mU32) + 16*pow3(mU32)) + 4*pow4(m32)))/(3.*pow2(m32 - mQ32)*
         pow3(m32 - mU32)*pow3(mQ32 - mU32)) + 2*((3*(mQ32 + mU32)*((8*(-(mQ32*
         mU32) + pow2(m32)))/(3.*(m32 - mQ32)*(m32 - mU32)) + (40*pow2(-(mQ32*
         mU32) + pow2(m32)))/(9.*pow2(m32 - mQ32)*pow2(m32 - mU32)) - 4*(
         -4.111111111111111 - (2*mQ32)/(3.*(-m32 + mQ32)) - (2*mU32)/(3.*(-m32
         + mU32)) + pow2(-2 + (2*mQ32)/(3.*(-m32 + mQ32)) + (2*mU32)/(3.*(-m32
         + mU32))) - (16*(-29.625 - (21*mQ32)/(4.*(-m32 + mQ32)) - (56*mQ32)/(
         mQ32 - mU32) + (56*mU32)/(mQ32 - mU32) + (111*mQ32)/(4.*(-m32 + mU32))
         + (45*mU32)/(2.*(-m32 + mU32)) + zt2/2. - (5*mQ32*zt2)/(2.*(-m32 +
         mQ32)) - (8*mQ32*zt2)/(mQ32 - mU32) + (8*mU32*zt2)/(mQ32 - mU32) + (4*
         mQ32*zt2)/(-m32 + mU32) + (3*mU32*zt2)/(2.*(-m32 + mU32)) + (2*(mQ32 +
         mU32)*fin(mQ32,mU32,MR2))/((m32 - mQ32)*(m32 - mU32)) + (56*pow2(mQ32
         ))/((-m32 + mQ32)*(mQ32 - mU32)) - (111*pow2(mQ32))/(4.*(-m32 + mQ32)*
         (-m32 + mU32)) + (8*zt2*pow2(mQ32))/((-m32 + mQ32)*(mQ32 - mU32)) - (4
         *zt2*pow2(mQ32))/((-m32 + mQ32)*(-m32 + mU32)) + (23*pow2(mQ32))/(8.*
         pow2(-m32 + mQ32)) + (3*zt2*pow2(mQ32))/(2.*pow2(-m32 + mQ32)) - (56*
         pow2(mU32))/((mQ32 - mU32)*(-m32 + mU32)) - (8*zt2*pow2(mU32))/((mQ32
         - mU32)*(-m32 + mU32)) + (23*pow2(mU32))/(8.*pow2(-m32 + mU32)) + (3*
         zt2*pow2(mU32))/(2.*pow2(-m32 + mU32)) + fin(mQ32,m32,MR2)*(-5/(2.*(
         -m32 + mQ32)) - 8/(mQ32 - mU32) + (8*mQ32)/((-m32 + mQ32)*(mQ32 - mU32
         )) + 2/(-m32 + mU32) - (4*mQ32)/((-m32 + mQ32)*(-m32 + mU32)) - (3*
         pow2(mQ32))/pow3(-m32 + mQ32)) - (21*pow3(mQ32))/pow3(-m32 + mQ32) - (
         3*zt2*pow3(mQ32))/pow3(-m32 + mQ32) + fin(mU32,m32,MR2)*(-2/(-m32 +
         mQ32) + 8/(mQ32 - mU32) + 3/(2.*(-m32 + mU32)) - (4*mQ32)/((-m32 +
         mQ32)*(-m32 + mU32)) - (8*mU32)/((mQ32 - mU32)*(-m32 + mU32)) - (3*
         pow2(mU32))/pow3(-m32 + mU32)) - (21*pow3(mU32))/pow3(-m32 + mU32) - (
         3*zt2*pow3(mU32))/pow3(-m32 + mU32)))/9. - (4*(-16.333333333333332 + (
         31*mQ32)/(-m32 + mQ32) + (45*msq2)/(-m32 + mQ32) + (9*mU32)/(2.*(-m32
         + mQ32)) + (9*mQ32)/(2.*(-m32 + mU32)) + (45*msq2)/(-m32 + mU32) + (31
         *mU32)/(-m32 + mU32) - 3*zt2 + (mQ32*zt2)/(4.*(-m32 + mQ32)) + (15*
         msq2*zt2)/(2.*(-m32 + mQ32)) + (3*mU32*zt2)/(4.*(-m32 + mQ32)) + (3*
         mQ32*zt2)/(4.*(-m32 + mU32)) + (15*msq2*zt2)/(2.*(-m32 + mU32)) + (
         mU32*zt2)/(4.*(-m32 + mU32)) - (27275*mQ32*msq2)/(432.*pow2(-m32 +
         mQ32)) - (2825*mQ32*mU32)/(432.*pow2(-m32 + mQ32)) - (65*mQ32*msq2*zt2
         )/(6.*pow2(-m32 + mQ32)) - (mQ32*mU32*zt2)/pow2(-m32 + mQ32) - (23*
         pow2(mQ32))/pow2(-m32 + mQ32) + (1355*msq2*pow2(mQ32))/(432.*(mQ32 -
         msq2)*pow2(-m32 + mQ32)) + (233*mU32*pow2(mQ32))/(432.*(mQ32 - mU32)*
         pow2(-m32 + mQ32)) + (5*msq2*zt2*pow2(mQ32))/(6.*(mQ32 - msq2)*pow2(
         -m32 + mQ32)) - (515*mQ32*pow2(msq2))/(144.*(mQ32 - msq2)*pow2(-m32 +
         mQ32)) - (5*mQ32*zt2*pow2(msq2))/(3.*(mQ32 - msq2)*pow2(-m32 + mQ32))
         + (95*pow2(mQ32)*pow2(msq2))/(216.*pow2(-m32 + mQ32)*pow2(mQ32 - msq2)
         ) + (5*zt2*pow2(mQ32)*pow2(msq2))/(6.*pow2(-m32 + mQ32)*pow2(mQ32 -
         msq2)) - (233*mQ32*pow2(mU32))/(288.*(mQ32 - mU32)*pow2(-m32 + mQ32))
         + (233*pow2(mQ32)*pow2(mU32))/(864.*pow2(-m32 + mQ32)*pow2(mQ32 - mU32
         )) + 3*(16.09722222222222 - (22*mQ32)/(-m32 + mQ32) - (14*mQ32)/(-m32
         + mU32) - (36*mU32)/(-m32 + mU32) + zt2 - (3*mQ32*zt2)/(2.*(-m32 +
         mQ32)) - (2*mQ32*zt2)/(-m32 + mU32) - (7*mU32*zt2)/(2.*(-m32 + mU32))
         - ((mQ32 + mU32)*fin(mQ32,mU32,MR2))/((m32 - mQ32)*(m32 - mU32)) + (14
         *pow2(mQ32))/((-m32 + mQ32)*(-m32 + mU32)) + (2*zt2*pow2(mQ32))/((-m32
         + mQ32)*(-m32 + mU32)) + (fin(mQ32,m32,MR2)*(3*m32*mQ32 - 3*mQ32*mU32
         + pow2(m32) - pow2(mQ32)))/((m32 - mU32)*pow2(-m32 + mQ32)) + (31*
         pow2(mQ32))/pow2(-m32 + mQ32) + (3*zt2*pow2(mQ32))/pow2(-m32 + mQ32) +
         (fin(mU32,m32,MR2)*(3*m32*mU32 - 3*mQ32*mU32 + pow2(m32) - pow2(mU32)
         ))/((m32 - mQ32)*pow2(m32 - mU32)) + (31*pow2(mU32))/pow2(-m32 + mU32)
         + (3*zt2*pow2(mU32))/pow2(-m32 + mU32)) - (6*mQ32*mU32)/pow2(-m32 +
         mU32) - (27275*msq2*mU32)/(432.*pow2(-m32 + mU32)) - (mQ32*mU32*zt2)
         /pow2(-m32 + mU32) - (65*msq2*mU32*zt2)/(6.*pow2(-m32 + mU32)) - (515*
         mU32*pow2(msq2))/(144.*(-msq2 + mU32)*pow2(-m32 + mU32)) - (5*mU32*zt2
         *pow2(msq2))/(3.*(-msq2 + mU32)*pow2(-m32 + mU32)) - (23*pow2(mU32))
         /pow2(-m32 + mU32) + (1355*msq2*pow2(mU32))/(432.*(-msq2 + mU32)*pow2(
         -m32 + mU32)) + (5*msq2*zt2*pow2(mU32))/(6.*(-msq2 + mU32)*pow2(-m32 +
         mU32)) + (95*pow2(msq2)*pow2(mU32))/(216.*pow2(-m32 + mU32)*pow2(
         -msq2 + mU32)) + (5*zt2*pow2(msq2)*pow2(mU32))/(6.*pow2(-m32 + mU32)*
         pow2(-msq2 + mU32)) - (5*fin(mQ32,msq2,MR2)*(-3*m32*mQ32 + 3*m32*msq2
         + mQ32*msq2 - 2*pow2(m32) + pow2(mQ32)))/(2.*pow3(-m32 + mQ32)) + (fin
         (mU32,m32,MR2)*(6*m32*mQ32*mU32 - 8*mQ32*pow2(m32) + 2*mU32*pow2(m32)
         + 3*m32*pow2(mQ32) - 3*m32*pow2(mU32) - mQ32*pow2(mU32) + 2*pow3(m32)
         - pow3(mQ32)))/(4.*(m32 - mU32)*pow3(-m32 + mQ32)) - (95*mQ32*pow3(
         msq2))/(216.*pow2(-m32 + mQ32)*pow2(mQ32 - msq2)) - (5*mQ32*zt2*pow3(
         msq2))/(6.*pow2(-m32 + mQ32)*pow2(mQ32 - msq2)) - (95*mU32*pow3(msq2))
         /(216.*pow2(-m32 + mU32)*pow2(-msq2 + mU32)) - (5*mU32*zt2*pow3(msq2))
         /(6.*pow2(-m32 + mU32)*pow2(-msq2 + mU32)) + (271*pow3(mU32))/(864.*(
         mQ32 - mU32)*pow2(-m32 + mQ32)) + (zt2*pow3(mU32))/(12.*(mQ32 - mU32)*
         pow2(-m32 + mQ32)) - (7*mQ32*pow3(mU32))/(12.*pow2(-m32 + mQ32)*pow2(
         mQ32 - mU32)) - (mQ32*zt2*pow3(mU32))/(12.*pow2(-m32 + mQ32)*pow2(mQ32
         - mU32)) + (fin(mQ32,m32,MR2)*(-6*m32*mQ32*mU32 - 2*mQ32*pow2(m32) +
         8*mU32*pow2(m32) + 3*m32*pow2(mQ32) + mU32*pow2(mQ32) - 3*m32*pow2(
         mU32) - 2*pow3(m32) + pow3(mU32)))/(4.*(m32 - mQ32)*pow3(m32 - mU32))
         + fin(m32,msq2,MR2)*(5/(2.*(-m32 + mQ32)) + 5/(2.*(-m32 + mU32)) + (15
         *mQ32)/(2.*pow2(-m32 + mQ32)) - (15*msq2)/(2.*pow2(-m32 + mQ32)) - (15
         *msq2)/(2.*pow2(-m32 + mU32)) + (15*mU32)/(2.*pow2(-m32 + mU32)) + (10
         *mQ32*msq2)/pow3(-m32 + mQ32) - (10*pow2(mQ32))/pow3(-m32 + mQ32) + (
         10*msq2*mU32)/pow3(-m32 + mU32) - (10*pow2(mU32))/pow3(-m32 + mU32)) -
         (5*fin(mU32,msq2,MR2)*(3*m32*msq2 - 3*m32*mU32 + msq2*mU32 - 2*pow2(
         m32) + pow2(mU32)))/(2.*pow3(-m32 + mU32)) + (271*pow4(mU32))/(864.*
         pow2(-m32 + mQ32)*pow2(mQ32 - mU32)) + (zt2*pow4(mU32))/(12.*pow2(-m32
         + mQ32)*pow2(mQ32 - mU32)) + (fin(mQ32,mU32,MR2)*(15*mU32*pow2(m32)*
         pow2(mQ32) + 15*mQ32*pow2(m32)*pow2(mU32) - 6*m32*pow2(mQ32)*pow2(mU32
         ) - 20*mQ32*mU32*pow3(m32) + 14*pow2(mQ32)*pow3(m32) + 14*pow2(mU32)*
         pow3(m32) - 6*m32*mU32*pow3(mQ32) - 11*pow2(m32)*pow3(mQ32) + pow2(
         mU32)*pow3(mQ32) - 6*m32*mQ32*pow3(mU32) - 11*pow2(m32)*pow3(mU32) +
         pow2(mQ32)*pow3(mU32) - 6*mQ32*pow4(m32) - 6*mU32*pow4(m32) + 3*m32*
         pow4(mQ32) + mU32*pow4(mQ32) + 3*m32*pow4(mU32) + mQ32*pow4(mU32) + 4*
         pow5(m32)))/(4.*pow3(-m32 + mQ32)*pow3(m32 - mU32))))/3.)))/pow3(mQ32
         - mU32) - (8*(-360*msq2*mU32*pow2(m32)*pow2(mQ32) + 420*mQ32*msq2*pow2
         (m32)*pow2(mU32) - 300*m32*msq2*pow2(mQ32)*pow2(mU32) + 773*pow2(m32)*
         pow2(mQ32)*pow2(mU32) + 180*mQ32*msq2*mU32*pow3(m32) - 5*mU32*pow2(
         mQ32)*pow3(m32) + 219*mQ32*pow2(mU32)*pow3(m32) - 180*msq2*pow2(mU32)*
         pow3(m32) + 180*m32*msq2*mU32*pow3(mQ32) + 23*mU32*pow2(m32)*pow3(mQ32
         ) - 559*m32*pow2(mU32)*pow3(mQ32) + 60*msq2*pow2(mU32)*pow3(mQ32) + 53
         *pow3(m32)*pow3(mQ32) + 120*m32*mQ32*msq2*pow3(mU32) - 1247*mQ32*pow2(
         m32)*pow3(mU32) - 60*msq2*pow2(m32)*pow3(mU32) + 279*m32*pow2(mQ32)*
         pow3(mU32) - 60*msq2*pow2(mQ32)*pow3(mU32) + 501*pow3(m32)*pow3(mU32)
         + 179*pow3(mQ32)*pow3(mU32) - 159*mQ32*mU32*pow4(m32) - 106*pow2(mQ32)
         *pow4(m32) - 247*pow2(mU32)*pow4(m32) + 18*m32*mU32*pow4(mQ32) + 6*
         pow2(mU32)*pow4(mQ32) + 426*m32*mQ32*pow4(mU32) - 61*pow2(m32)*pow4(
         mU32) - 173*pow2(mQ32)*pow4(mU32) + 53*mQ32*pow5(m32) + 75*mU32*pow5(
         m32) - 36*m32*pow5(mU32) - 12*mQ32*pow5(mU32)))/(3.*pow2(m32 - mQ32)*
         pow3(m32 - mU32)*pow3(mQ32 - mU32)))) + pow2(log(msq2/MR2))*((320*m3*
         msq2*Xt)/((m32 - mQ32)*(m32 - mU32)) - (640*m3*msq2*pow3(Xt))/((m32 -
         mQ32)*(m32 - mU32)*(mQ32 - mU32)) - (40*(-4*m32*mQ32*msq2*mU32 - 5*
         mQ32*msq2*pow2(m32) - 8*mQ32*mU32*pow2(m32) - 5*msq2*mU32*pow2(m32) +
         3*m32*msq2*pow2(mQ32) + 4*m32*mU32*pow2(mQ32) + msq2*mU32*pow2(mQ32) -
         2*pow2(m32)*pow2(mQ32) + 4*m32*mQ32*pow2(mU32) + 3*m32*msq2*pow2(mU32
         ) + mQ32*msq2*pow2(mU32) - 2*pow2(m32)*pow2(mU32) - 2*pow2(mQ32)*pow2(
         mU32) + 4*mQ32*pow3(m32) + 6*msq2*pow3(m32) + 4*mU32*pow3(m32) - 2*
         pow4(m32)))/(pow2(-m32 + mQ32)*pow2(m32 - mU32)) + (80*pow2(Xt)*(-4*
         m32*mQ32*msq2*mU32 - 5*mQ32*msq2*pow2(m32) - 8*mQ32*mU32*pow2(m32) - 5
         *msq2*mU32*pow2(m32) + 3*m32*msq2*pow2(mQ32) + 4*m32*mU32*pow2(mQ32) +
         msq2*mU32*pow2(mQ32) - 2*pow2(m32)*pow2(mQ32) + 4*m32*mQ32*pow2(mU32)
         + 3*m32*msq2*pow2(mU32) + mQ32*msq2*pow2(mU32) - 2*pow2(m32)*pow2(
         mU32) - 2*pow2(mQ32)*pow2(mU32) + 4*mQ32*pow3(m32) + 6*msq2*pow3(m32)
         + 4*mU32*pow3(m32) - 2*pow4(m32)))/((mQ32 - mU32)*pow2(-m32 + mQ32)*
         pow2(m32 - mU32)) - (40*(mQ32 + mU32)*(-4*m32*mQ32*msq2*mU32 - 5*mQ32*
         msq2*pow2(m32) - 8*mQ32*mU32*pow2(m32) - 5*msq2*mU32*pow2(m32) + 3*m32
         *msq2*pow2(mQ32) + 4*m32*mU32*pow2(mQ32) + msq2*mU32*pow2(mQ32) - 2*
         pow2(m32)*pow2(mQ32) + 4*m32*mQ32*pow2(mU32) + 3*m32*msq2*pow2(mU32) +
         mQ32*msq2*pow2(mU32) - 2*pow2(m32)*pow2(mU32) - 2*pow2(mQ32)*pow2(
         mU32) + 4*mQ32*pow3(m32) + 6*msq2*pow3(m32) + 4*mU32*pow3(m32) - 2*
         pow4(m32))*pow4(Xt))/(pow2(-m32 + mQ32)*pow2(m32 - mU32)*pow3(mQ32 -
         mU32)) + (320*m3*msq2*(mQ32 + mU32)*pow5(Xt))/((m32 - mQ32)*(m32 -
         mU32)*pow3(mQ32 - mU32))) + log(msq2/MR2)*(Xt*((-2560*m3*msq2)/((m32 -
         mQ32)*(m32 - mU32)) + (192*m3)/(mQ32 - mU32) - (96*m3*mU32)/((m32 -
         mU32)*(-mQ32 + mU32)) + (320*m3*(m32 + 6*msq2 - mU32)*mU32)/(3.*(-mQ32
         + mU32)*pow2(m32 - mU32))) + ((5120*m3*msq2)/((m32 - mQ32)*(m32 -
         mU32)*(mQ32 - mU32)) - (192*m3*(mQ32 + 3*mU32))/pow3(mQ32 - mU32))*
         pow3(Xt) - (32*pow2(Xt)*(-240*m32*mQ32*msq2*mU32 - 300*mQ32*msq2*pow2(
         m32) + 330*mQ32*mU32*pow2(m32) - 300*msq2*mU32*pow2(m32) + 180*m32*
         msq2*pow2(mQ32) - 160*m32*mU32*pow2(mQ32) + 60*msq2*mU32*pow2(mQ32) +
         87*pow2(m32)*pow2(mQ32) - 151*m32*mQ32*pow2(mU32) + 180*m32*msq2*pow2(
         mU32) + 60*mQ32*msq2*pow2(mU32) + 78*pow2(m32)*pow2(mU32) + 73*pow2(
         mQ32)*pow2(mU32) - 179*mQ32*pow3(m32) + 360*msq2*pow3(m32) - 170*mU32*
         pow3(m32) + 92*pow4(m32)))/(3.*(mQ32 - mU32)*pow2(m32 - mQ32)*pow2(m32
         - mU32)) + ((24*(-3*mQ32*(mQ32 + 5*mU32)*pow2(mU32) - pow2(m32)*(15*
         mQ32*mU32 + 6*pow2(mQ32) + 29*pow2(mU32)) + (3*mQ32 + 7*mU32)*pow3(m32
         ) + m32*(7*mU32*pow2(mQ32) + 31*mQ32*pow2(mU32) + 16*pow3(mU32)) + 4*
         pow4(m32)))/((m32 - mQ32)*pow2(m32 - mU32)*pow3(mQ32 - mU32)) + 2*((8*
         (90*m32*msq2*mU32 + 2*mU32*pow2(m32) - 3*m32*pow2(mU32) + 30*msq2*pow2
         (mU32) + pow3(mU32)))/(3.*pow2(mQ32 - mU32)*pow3(m32 - mU32)) + (4*(
         mQ32 + mU32)*(-480*m32*mQ32*msq2*mU32 - 600*mQ32*msq2*pow2(m32) + 480*
         mQ32*mU32*pow2(m32) - 600*msq2*mU32*pow2(m32) + 360*m32*msq2*pow2(mQ32
         ) - 239*m32*mU32*pow2(mQ32) + 120*msq2*mU32*pow2(mQ32) + 120*pow2(m32)
         *pow2(mQ32) - 239*m32*mQ32*pow2(mU32) + 360*m32*msq2*pow2(mU32) + 120*
         mQ32*msq2*pow2(mU32) + 120*pow2(m32)*pow2(mU32) + 119*pow2(mQ32)*pow2(
         mU32) - 241*mQ32*pow3(m32) + 720*msq2*pow3(m32) - 241*mU32*pow3(m32) +
         121*pow4(m32)))/(3.*pow2(m32 - mQ32)*pow2(m32 - mU32)*pow3(mQ32 -
         mU32))))*pow4(Xt) + (16*(150*mQ32*msq2*mU32*pow2(m32) - 165*m32*msq2*
         mU32*pow2(mQ32) + 180*msq2*pow2(m32)*pow2(mQ32) - 230*mU32*pow2(m32)*
         pow2(mQ32) + 330*m32*mQ32*msq2*pow2(mU32) - 439*mQ32*pow2(m32)*pow2(
         mU32) + 465*msq2*pow2(m32)*pow2(mU32) + 212*m32*pow2(mQ32)*pow2(mU32)
         - 75*msq2*pow2(mQ32)*pow2(mU32) - 300*mQ32*msq2*pow3(m32) + 475*mQ32*
         mU32*pow3(m32) - 705*msq2*mU32*pow3(m32) + 87*pow2(mQ32)*pow3(m32) +
         227*pow2(mU32)*pow3(m32) + 143*m32*mQ32*pow3(mU32) - 180*m32*msq2*pow3
         (mU32) - 60*mQ32*msq2*pow3(mU32) - 74*pow2(m32)*pow3(mU32) - 69*pow2(
         mQ32)*pow3(mU32) - 179*mQ32*pow4(m32) + 360*msq2*pow4(m32) - 245*mU32*
         pow4(m32) + 92*pow5(m32)))/(3.*pow2(m32 - mQ32)*pow3(m32 - mU32)) - (
         640*m3*(-12*m32*mQ32*msq2 - m32*mQ32*mU32 - 6*m32*msq2*mU32 + 6*mQ32*
         msq2*mU32 + mU32*pow2(m32) - m32*pow2(mU32) + mQ32*pow2(mU32) + 12*
         msq2*pow2(mU32))*pow5(Xt))/(3.*(m32 - mQ32)*pow2(m32 - mU32)*pow3(
         -mQ32 + mU32))) + Xt*((704*m3*mU32)/(3.*(m32 - mU32)*(-mQ32 + mU32)) +
         (128*m3*(3*m32*mQ32 + 3*m32*mU32 + 5*mQ32*mU32 - 11*pow2(m32)))/(3.*(
         m32 - mQ32)*(m32 - mU32)*(mQ32 - mU32)) - (384*m3*mU32*(-(mQ32*mU32) +
         pow2(m32)))/((m32 - mQ32)*(mQ32 - mU32)*pow2(m32 - mU32)) + (128*(-(
         mU32*dilog(1 - m32/mQ32)) + m32*(dilog(1 - m32/mQ32) - dilog(1 -
         m32/mU32)) + mQ32*dilog(1 - m32/mU32))*(-34*m32*mU32 + pow2(m32) + 17*
         pow2(mU32))*pow3(m3))/(3.*(m32 - mQ32)*(mQ32 - mU32)*pow3(m32 - mU32))
         - 24*((64*(m3*fin(mQ32,m32,MR2) - m3*fin(mQ32,mU32,MR2) + m3*fin(mU32
         ,m32,MR2) + 3*pow3(m3) + zt2*pow3(m3)))/(9.*(m32 - mQ32)*(m32 - mU32))
         - (4*((20*(m3*mQ32 - m3*msq2)*fin(mQ32,msq2,MR2))/((mQ32 - mU32)*pow2
         (m32 - mQ32)) - (2*m3*(mQ32 - mU32)*fin(mU32,m32,MR2))/((m32 - mU32)*
         pow2(m32 - mQ32)) + (2*m3*(mQ32 - mU32)*fin(mQ32,m32,MR2))/((m32 -
         mQ32)*pow2(m32 - mU32)) + (20*(-(m3*msq2) + m3*mU32)*fin(mU32,msq2,MR2
         ))/((-mQ32 + mU32)*pow2(m32 - mU32)) + (2*(6*m3*mQ32 + 60*m3*msq2 + 6*
         m3*mU32 + m3*mQ32*zt2 + 10*m3*msq2*zt2 + m3*mU32*zt2 + 24*pow3(m3)))/(
         (m32 - mQ32)*(m32 - mU32)) - (20*fin(m32,msq2,MR2)*(m3*mQ32*msq2 - 2*
         m3*mQ32*mU32 + m3*msq2*mU32 + mQ32*pow3(m3) - 2*msq2*pow3(m3) + mU32*
         pow3(m3)))/(pow2(m32 - mQ32)*pow2(m32 - mU32)) + 3*((-2*m3*fin(mQ32,
         mU32,MR2))/((m32 - mQ32)*(m32 - mU32)) + (2*fin(mQ32,m32,MR2)*(m3*mQ32
         - 3*m3*mU32 + 2*pow3(m3)))/((m32 - mQ32)*(m32 - mU32)*(mQ32 - mU32))
         - (2*fin(mU32,m32,MR2)*(-3*m3*mQ32 + m3*mU32 + 2*pow3(m3)))/((m32 -
         mQ32)*(m32 - mU32)*(mQ32 - mU32)) + (2*(11*pow3(m3) + 3*zt2*pow3(m3)))
         /((m32 - mQ32)*(m32 - mU32))) + (2*fin(mQ32,mU32,MR2)*(m3*pow2(mQ32) +
         m3*pow2(mU32) - 2*mQ32*pow3(m3) - 2*mU32*pow3(m3) + 2*pow5(m3)))/(
         pow2(m32 - mQ32)*pow2(m32 - mU32))))/3.) + (512*m3*(9*m32*(mQ32 + mU32
         )*pow2(mQ32)*pow2(mU32) - 3*pow3(mQ32)*pow3(mU32) + pow3(m32)*(mU32*(
         19 + 8*dilog(1 - m32/mU32))*pow2(mQ32) + mQ32*(19 + 8*dilog(1 -
         m32/mQ32))*pow2(mU32) + 2*pow3(mQ32) + 2*pow3(mU32)) - 4*pow2(m32)*(6*
         pow2(mQ32)*pow2(mU32) + mU32*(2 + dilog(1 - m32/mU32))*pow3(mQ32) +
         mQ32*(2 + dilog(1 - m32/mQ32))*pow3(mU32)) - (mQ32*mU32*(13 + 4*dilog(
         1 - m32/mQ32) + 4*dilog(1 - m32/mU32)) + 4*pow2(mQ32) + 4*pow2(mU32))*
         pow4(m32) + 2*(mQ32 + mU32)*pow5(m32)))/(3.*mQ32*(-mQ32 + mU32)*pow2(
         m32 - mQ32)*pow3(m32 - mU32)) - (64*(-60*m3*mQ32*msq2*mU32 - 6*m3*mU32
         *pow2(mQ32) + 89*m3*mQ32*pow2(mU32) + 19*mQ32*mU32*pow3(m3) + 60*msq2*
         mU32*pow3(m3) - 115*pow2(mU32)*pow3(m3) + 12*m3*pow3(mU32) - 16*mQ32*
         pow5(m3) + mU32*pow5(m3) + 16*pow7(m3)))/(3.*(m32 - mQ32)*(-mQ32 +
         mU32)*pow2(m32 - mU32))) + pow5(Xt)*((1024*m3*((2*mQ32*mU32*(-dilog(1
         - m32/mQ32) + dilog(1 - m32/mU32)) - pow2(mQ32) + pow2(mU32))*pow3(m32
         ) + (-5 + dilog(1 - m32/mQ32) - dilog(1 - m32/mU32))*pow2(mU32)*pow3(
         mQ32) + pow2(m32)*(3*mU32*(-2 + dilog(1 - m32/mQ32) - dilog(1 -
         m32/mU32))*pow2(mQ32) + 3*mQ32*(2 + dilog(1 - m32/mQ32) - dilog(1 -
         m32/mU32))*pow2(mU32) + pow3(mQ32) - pow3(mU32)) + (5 + dilog(1 -
         m32/mQ32) - dilog(1 - m32/mU32))*pow2(mQ32)*pow3(mU32) + m32*(4*(
         -dilog(1 - m32/mQ32) + dilog(1 - m32/mU32))*pow2(mQ32)*pow2(mU32) +
         mU32*(5 - dilog(1 - m32/mQ32) + dilog(1 - m32/mU32))*pow3(mQ32) + mQ32
         *(-5 - dilog(1 - m32/mQ32) + dilog(1 - m32/mU32))*pow3(mU32))))/(3.*
         mQ32*(-m32 + mQ32)*(-mQ32 + mU32)*pow2(m32 - mU32)*pow3(mQ32 - mU32))
         - (128*m3*mU32*(3*m32*mQ32 + 3*m32*mU32 + 5*mQ32*mU32 - 11*pow2(m32)))
         /(3.*(m32 - mQ32)*pow2(m32 - mU32)*pow3(-mQ32 + mU32)) + 2*((-12*(mQ32
         + mU32)*((64*(m3*fin(mQ32,m32,MR2) - m3*fin(mQ32,mU32,MR2) + m3*fin(
         mU32,m32,MR2) + 3*pow3(m3) + zt2*pow3(m3)))/(9.*(m32 - mQ32)*(m32 -
         mU32)) - (4*((20*(m3*mQ32 - m3*msq2)*fin(mQ32,msq2,MR2))/((mQ32 - mU32
         )*pow2(m32 - mQ32)) - (2*m3*(mQ32 - mU32)*fin(mU32,m32,MR2))/((m32 -
         mU32)*pow2(m32 - mQ32)) + (2*m3*(mQ32 - mU32)*fin(mQ32,m32,MR2))/((m32
         - mQ32)*pow2(m32 - mU32)) + (20*(-(m3*msq2) + m3*mU32)*fin(mU32,msq2,
         MR2))/((-mQ32 + mU32)*pow2(m32 - mU32)) + (2*(6*m3*mQ32 + 60*m3*msq2 +
         6*m3*mU32 + m3*mQ32*zt2 + 10*m3*msq2*zt2 + m3*mU32*zt2 + 24*pow3(m3))
         )/((m32 - mQ32)*(m32 - mU32)) - (20*fin(m32,msq2,MR2)*(m3*mQ32*msq2 -
         2*m3*mQ32*mU32 + m3*msq2*mU32 + mQ32*pow3(m3) - 2*msq2*pow3(m3) + mU32
         *pow3(m3)))/(pow2(m32 - mQ32)*pow2(m32 - mU32)) + 3*((-2*m3*fin(mQ32,
         mU32,MR2))/((m32 - mQ32)*(m32 - mU32)) + (2*fin(mQ32,m32,MR2)*(m3*mQ32
         - 3*m3*mU32 + 2*pow3(m3)))/((m32 - mQ32)*(m32 - mU32)*(mQ32 - mU32))
         - (2*fin(mU32,m32,MR2)*(-3*m3*mQ32 + m3*mU32 + 2*pow3(m3)))/((m32 -
         mQ32)*(m32 - mU32)*(mQ32 - mU32)) + (2*(11*pow3(m3) + 3*zt2*pow3(m3)))
         /((m32 - mQ32)*(m32 - mU32))) + (2*fin(mQ32,mU32,MR2)*(m3*pow2(mQ32) +
         m3*pow2(mU32) - 2*mQ32*pow3(m3) - 2*mU32*pow3(m3) + 2*pow5(m3)))/(
         pow2(m32 - mQ32)*pow2(m32 - mU32))))/3.))/pow3(mQ32 - mU32) - (128*(
         -30*m3*mQ32*msq2*mU32 - 3*m3*mU32*pow2(mQ32) + 48*m3*mQ32*pow2(mU32) +
         11*mQ32*mU32*pow3(m3) + 30*msq2*mU32*pow3(m3) - 56*pow2(mU32)*pow3(m3
         ) + 6*m3*pow3(mU32) - 8*mQ32*pow5(m3) - 6*mU32*pow5(m3) + 8*pow7(m3)))
         /(3.*(m32 - mQ32)*pow2(m32 - mU32)*pow3(mQ32 - mU32))))) + pow2(log(
         mU32/MR2))*((-16*pow2(Xt)*(360*msq2*mU32*pow2(m32)*pow2(mQ32) - 420*
         mQ32*msq2*pow2(m32)*pow2(mU32) + 300*m32*msq2*pow2(mQ32)*pow2(mU32) -
         1447*pow2(m32)*pow2(mQ32)*pow2(mU32) - 180*mQ32*msq2*mU32*pow3(m32) +
         794*mU32*pow2(mQ32)*pow3(m32) + 243*mQ32*pow2(mU32)*pow3(m32) + 180*
         msq2*pow2(mU32)*pow3(m32) - 180*m32*msq2*mU32*pow3(mQ32) - 244*mU32*
         pow2(m32)*pow3(mQ32) + 735*m32*pow2(mU32)*pow3(mQ32) - 60*msq2*pow2(
         mU32)*pow3(mQ32) - 41*pow3(m32)*pow3(mQ32) - 120*m32*mQ32*msq2*pow3(
         mU32) + 1418*mQ32*pow2(m32)*pow3(mU32) + 60*msq2*pow2(m32)*pow3(mU32)
         - 172*m32*pow2(mQ32)*pow3(mU32) + 60*msq2*pow2(mQ32)*pow3(mU32) - 748*
         pow3(m32)*pow3(mU32) - 210*pow3(mQ32)*pow3(mU32) - 592*mQ32*mU32*pow4(
         m32) + 65*pow2(mQ32)*pow4(m32) + 283*pow2(mU32)*pow4(m32) - 18*m32*
         mU32*pow4(mQ32) - 6*pow2(mU32)*pow4(mQ32) - 529*m32*mQ32*pow4(mU32) +
         153*pow2(m32)*pow4(mU32) + 184*pow2(mQ32)*pow4(mU32) - 24*mQ32*pow5(
         m32) + 108*mU32*pow5(m32) + 36*m32*pow5(mU32) + 12*mQ32*pow5(mU32)))/(
         3.*pow2(m32 - mQ32)*pow2(-mQ32 + mU32)*pow3(m32 - mU32)) + log(
         msq2/MR2)*((-8*(-90*m32*msq2*mU32 - 65*mU32*pow2(m32) + 57*m32*pow2(
         mU32) - 30*msq2*pow2(mU32) + 27*pow3(m32) - 19*pow3(mU32)))/(3.*pow3(
         m32 - mU32)) + pow2(Xt)*((48*((mQ32 - 3*mU32)*pow2(m32) + 2*(mQ32 - 2*
         mU32)*pow2(mU32) + m32*(-4*mQ32*mU32 + 8*pow2(mU32))))/(pow2(m32 -
         mU32)*pow2(mQ32 - mU32)) - (16*(90*m32*msq2*mU32 + 2*mU32*pow2(m32) -
         3*m32*pow2(mU32) + 30*msq2*pow2(mU32) + pow3(mU32)))/(3.*(mQ32 - mU32)
         *pow3(m32 - mU32))) + ((-64*m3*(m32 + 60*msq2 - mU32)*mU32)/(3.*pow2(
         m32 - mU32)*pow2(-mQ32 + mU32)) + (96*m3*(3*mQ32*mU32 - m32*(mQ32 + 3*
         mU32) + 2*pow2(m32) - pow2(mU32)))/((m32 - mU32)*pow3(-mQ32 + mU32)))*
         pow3(Xt) + ((8*(mQ32 + mU32)*(90*m32*msq2*mU32 + 2*mU32*pow2(m32) - 3*
         m32*pow2(mU32) + 30*msq2*pow2(mU32) + pow3(mU32)))/(3.*pow3(m32 - mU32
         )*pow3(mQ32 - mU32)) + (24*(8*mQ32*mU32*pow2(m32) + 2*m32*mU32*(-5*
         mQ32*mU32 + pow2(mQ32) - 4*pow2(mU32)) - pow2(mQ32)*pow2(mU32) - 2*(
         mQ32 - mU32)*pow3(m32) + 4*mQ32*pow3(mU32) + 5*pow4(mU32)))/(pow2(m32
         - mU32)*pow4(mQ32 - mU32)))*pow4(Xt) + (64*Xt*(-30*m3*msq2*mU32 + 5*m3
         *pow2(mU32) - 14*mU32*pow3(m3) + 9*pow5(m3)))/(3.*(-mQ32 + mU32)*pow2(
         m32 - mU32)) + (320*m3*(m32 + 6*msq2 - mU32)*mU32*(mQ32 + mU32)*pow5(
         Xt))/(3.*pow2(m32 - mU32)*pow4(-mQ32 + mU32))) + (4*(150*msq2*pow2(m32
         )*pow2(mQ32)*pow2(mU32) - 900*msq2*mU32*pow2(mQ32)*pow3(m32) + 600*
         mQ32*msq2*pow2(mU32)*pow3(m32) + 2355*pow2(mQ32)*pow2(mU32)*pow3(m32)
         + 450*msq2*mU32*pow2(m32)*pow3(mQ32) - 300*m32*msq2*pow2(mU32)*pow3(
         mQ32) - 1673*pow2(m32)*pow2(mU32)*pow3(mQ32) + 375*mU32*pow3(m32)*pow3
         (mQ32) - 750*mQ32*msq2*pow2(m32)*pow3(mU32) + 600*m32*msq2*pow2(mQ32)*
         pow3(mU32) - 1055*pow2(m32)*pow2(mQ32)*pow3(mU32) - 4711*mQ32*pow3(m32
         )*pow3(mU32) + 300*msq2*pow3(m32)*pow3(mU32) + 1757*m32*pow3(mQ32)*
         pow3(mU32) - 150*msq2*pow3(mQ32)*pow3(mU32) + 450*mQ32*msq2*mU32*pow4(
         m32) - 920*mU32*pow2(mQ32)*pow4(m32) + 1752*mQ32*pow2(mU32)*pow4(m32)
         - 450*msq2*pow2(mU32)*pow4(m32) + 64*pow3(mQ32)*pow4(m32) + 2944*pow3(
         mU32)*pow4(m32) + 45*mU32*pow2(m32)*pow4(mQ32) - 30*m32*pow2(mU32)*
         pow4(mQ32) - 15*pow3(mU32)*pow4(mQ32) - 300*m32*mQ32*msq2*pow4(mU32) +
         4459*mQ32*pow2(m32)*pow4(mU32) + 150*msq2*pow2(m32)*pow4(mU32) - 1133
         *m32*pow2(mQ32)*pow4(mU32) + 150*msq2*pow2(mQ32)*pow4(mU32) - 1859*
         pow3(m32)*pow4(mU32) - 387*pow3(mQ32)*pow4(mU32) - 26*mQ32*mU32*pow5(
         m32) - 30*pow2(mQ32)*pow5(m32) - 1864*pow2(mU32)*pow5(m32) - 1059*m32*
         mQ32*pow5(mU32) + 144*pow2(m32)*pow5(mU32) + 375*pow2(mQ32)*pow5(mU32)
         - 34*mQ32*pow6(m32) + 418*mU32*pow6(m32) + 81*m32*pow6(mU32) + 27*
         mQ32*pow6(mU32)))/(3.*(-mQ32 + mU32)*pow2(m32 - mQ32)*pow4(m32 - mU32)
         ) - (2560*m32*pow2(mU32)*pow6(Xt))/(3.*pow2(m32 - mU32)*pow4(-mQ32 +
         mU32)) - (8*pow4(Xt)*(-420*msq2*pow2(mQ32)*pow2(mU32)*pow3(m32) + 300*
         msq2*pow2(m32)*pow2(mU32)*pow3(mQ32) - 180*msq2*mU32*pow3(m32)*pow3(
         mQ32) + 3248*pow2(mU32)*pow3(m32)*pow3(mQ32) - 60*msq2*pow2(m32)*pow2(
         mQ32)*pow3(mU32) + 420*mQ32*msq2*pow3(m32)*pow3(mU32) + 8872*pow2(mQ32
         )*pow3(m32)*pow3(mU32) - 60*m32*msq2*pow3(mQ32)*pow3(mU32) - 4220*pow2
         (m32)*pow3(mQ32)*pow3(mU32) + 90*msq2*mU32*pow2(mQ32)*pow4(m32) + 180*
         mQ32*msq2*pow2(mU32)*pow4(m32) - 4816*pow2(mQ32)*pow2(mU32)*pow4(m32)
         - 796*mU32*pow3(mQ32)*pow4(m32) - 9028*mQ32*pow3(mU32)*pow4(m32) - 270
         *msq2*pow3(mU32)*pow4(m32) + 90*msq2*mU32*pow2(m32)*pow4(mQ32) - 60*
         m32*msq2*pow2(mU32)*pow4(mQ32) - 803*pow2(m32)*pow2(mU32)*pow4(mQ32) +
         200*mU32*pow3(m32)*pow4(mQ32) + 782*m32*pow3(mU32)*pow4(mQ32) - 30*
         msq2*pow3(mU32)*pow4(mQ32) + 41*pow4(m32)*pow4(mQ32) - 420*mQ32*msq2*
         pow2(m32)*pow4(mU32) + 300*m32*msq2*pow2(mQ32)*pow4(mU32) - 5272*pow2(
         m32)*pow2(mQ32)*pow4(mU32) + 5088*mQ32*pow3(m32)*pow4(mU32) + 180*msq2
         *pow3(m32)*pow4(mU32) + 2117*m32*pow3(mQ32)*pow4(mU32) - 60*msq2*pow3(
         mQ32)*pow4(mU32) - 1961*pow4(m32)*pow4(mU32) - 212*pow4(mQ32)*pow4(
         mU32) + 699*mU32*pow2(mQ32)*pow5(m32) + 4021*mQ32*pow2(mU32)*pow5(m32)
         - 97*pow3(mQ32)*pow5(m32) + 3537*pow3(mU32)*pow5(m32) + 9*mU32*pow2(
         m32)*pow5(mQ32) - 6*m32*pow2(mU32)*pow5(mQ32) - 3*pow3(mU32)*pow5(mQ32
         ) - 180*m32*mQ32*msq2*pow5(mU32) + 1095*mQ32*pow2(m32)*pow5(mU32) + 90
         *msq2*pow2(m32)*pow5(mU32) + 57*m32*pow2(mQ32)*pow5(mU32) + 90*msq2*
         pow2(mQ32)*pow5(mU32) - 448*pow3(m32)*pow5(mU32) - 284*pow3(mQ32)*pow5
         (mU32) - 196*mQ32*mU32*pow6(m32) + 104*pow2(mQ32)*pow6(m32) - 1524*
         pow2(mU32)*pow6(m32) - 979*m32*mQ32*pow6(mU32) + 311*pow2(m32)*pow6(
         mU32) + 404*pow2(mQ32)*pow6(mU32) - 48*mQ32*pow7(m32) + 48*mU32*pow7(
         m32) + 45*m32*pow7(mU32) + 15*mQ32*pow7(mU32)))/(3.*pow2(m32 - mQ32)*
         pow4(m32 - mU32)*pow4(-mQ32 + mU32)) + (128*pow3(Xt)*(-60*m3*msq2*mU32
         *pow2(mQ32) + 60*m3*mQ32*msq2*pow2(mU32) + 127*m3*pow2(mQ32)*pow2(mU32
         ) + 60*mQ32*msq2*mU32*pow3(m3) - 45*mU32*pow2(mQ32)*pow3(m3) - 173*
         mQ32*pow2(mU32)*pow3(m3) - 60*msq2*pow2(mU32)*pow3(m3) - 6*m3*mU32*
         pow3(mQ32) - 80*m3*mQ32*pow3(mU32) + 120*pow3(m3)*pow3(mU32) - 12*m3*
         pow4(mU32) + 96*mQ32*mU32*pow5(m3) - 14*pow2(mQ32)*pow5(m3) + 6*pow2(
         mU32)*pow5(m3) + 7*mQ32*pow7(m3) - 37*mU32*pow7(m3) + 11*pow9(m3)))/(
         3.*(m32 - mQ32)*pow2(m32 - mU32)*pow3(mQ32 - mU32)) - (32*Xt*(150*m3*
         msq2*pow2(mQ32)*pow2(mU32) - 150*msq2*mU32*pow2(mQ32)*pow3(m3) + 283*
         pow2(mQ32)*pow2(mU32)*pow3(m3) + 15*m3*pow2(mU32)*pow3(mQ32) - 15*mU32
         *pow3(m3)*pow3(mQ32) - 150*m3*mQ32*msq2*pow3(mU32) - 265*m3*pow2(mQ32)
         *pow3(mU32) + 47*mQ32*pow3(m3)*pow3(mU32) + 150*msq2*pow3(m3)*pow3(
         mU32) + 239*m3*mQ32*pow4(mU32) - 363*pow3(m3)*pow4(mU32) + 150*mQ32*
         msq2*mU32*pow5(m3) - 80*mU32*pow2(mQ32)*pow5(m3) - 332*mQ32*pow2(mU32)
         *pow5(m3) - 150*msq2*pow2(mU32)*pow5(m3) + 460*pow3(mU32)*pow5(m3) +
         27*m3*pow5(mU32) + 220*mQ32*mU32*pow7(m3) - 18*pow2(mQ32)*pow7(m3) -
         218*pow2(mU32)*pow7(m3) - 14*mQ32*pow9(m3) + 14*mU32*pow9(m3)))/(3.*(
         m32 - mQ32)*pow2(-mQ32 + mU32)*pow3(m32 - mU32)) - (64*pow5(Xt)*(30*m3
         *msq2*pow2(mQ32)*pow2(mU32) - 30*msq2*mU32*pow2(mQ32)*pow3(m3) - 120*
         mQ32*msq2*pow2(mU32)*pow3(m3) + 91*pow2(mQ32)*pow2(mU32)*pow3(m3) + 3*
         m3*pow2(mU32)*pow3(mQ32) - 3*mU32*pow3(m3)*pow3(mQ32) + 90*m3*mQ32*
         msq2*pow3(mU32) - 71*m3*pow2(mQ32)*pow3(mU32) + 463*mQ32*pow3(m3)*pow3
         (mU32) - 90*msq2*pow3(m3)*pow3(mU32) - 271*m3*mQ32*pow4(mU32) + 315*
         pow3(m3)*pow4(mU32) + 30*mQ32*msq2*mU32*pow5(m3) - 20*mU32*pow2(mQ32)*
         pow5(m3) - 168*mQ32*pow2(mU32)*pow5(m3) + 90*msq2*pow2(mU32)*pow5(m3)
         - 418*pow3(mU32)*pow5(m3) - 15*m3*pow5(mU32) - 8*mQ32*mU32*pow7(m3) -
         16*pow2(mQ32)*pow7(m3) + 54*pow2(mU32)*pow7(m3) + 16*mQ32*pow9(m3) +
         48*mU32*pow9(m3)))/(3.*(m32 - mQ32)*pow3(m32 - mU32)*pow4(mQ32 - mU32)
         )) + pow3(log(mQ32/MR2))*((8*pow2(Xt)*(-30*msq2*pow2(m32)*pow2(mQ32)*
         pow2(mU32) - 120*msq2*mU32*pow2(mQ32)*pow3(m32) + 180*mQ32*msq2*pow2(
         mU32)*pow3(m32) - 1405*pow2(mQ32)*pow2(mU32)*pow3(m32) + 150*msq2*mU32
         *pow2(m32)*pow3(mQ32) - 120*m32*msq2*pow2(mU32)*pow3(mQ32) + 501*pow2(
         m32)*pow2(mU32)*pow3(mQ32) - 60*msq2*pow3(m32)*pow3(mQ32) + 333*mU32*
         pow3(m32)*pow3(mQ32) - 90*mQ32*msq2*pow2(m32)*pow3(mU32) + 60*m32*msq2
         *pow2(mQ32)*pow3(mU32) + 415*pow2(m32)*pow2(mQ32)*pow3(mU32) - 77*mQ32
         *pow3(m32)*pow3(mU32) - 415*m32*pow3(mQ32)*pow3(mU32) + 30*msq2*pow3(
         mQ32)*pow3(mU32) - 90*mQ32*msq2*mU32*pow4(m32) + 90*msq2*pow2(mQ32)*
         pow4(m32) + 1478*mU32*pow2(mQ32)*pow4(m32) + 562*mQ32*pow2(mU32)*pow4(
         m32) - 302*pow3(mQ32)*pow4(m32) - 2*pow3(mU32)*pow4(m32) + 60*m32*msq2
         *mU32*pow4(mQ32) - 30*msq2*pow2(m32)*pow4(mQ32) - 1285*mU32*pow2(m32)*
         pow4(mQ32) + 473*m32*pow2(mU32)*pow4(mQ32) - 30*msq2*pow2(mU32)*pow4(
         mQ32) + 605*pow3(m32)*pow4(mQ32) + 87*pow3(mU32)*pow4(mQ32) - 9*mQ32*
         pow2(m32)*pow4(mU32) + 6*m32*pow2(mQ32)*pow4(mU32) + 3*pow3(mQ32)*pow4
         (mU32) - 842*mQ32*mU32*pow5(m32) - 554*pow2(mQ32)*pow5(m32) + 4*pow2(
         mU32)*pow5(m32) + 345*m32*mU32*pow5(mQ32) - 126*pow2(m32)*pow5(mQ32) -
         159*pow2(mU32)*pow5(mQ32) + 378*mQ32*pow6(m32) - 2*mU32*pow6(m32) - 9
         *m32*pow6(mQ32) - 3*mU32*pow6(mQ32)))/(3.*pow2(m32 - mU32)*pow2(mQ32 -
         mU32)*pow4(m32 - mQ32)) - (4*(30*msq2*pow2(m32)*pow2(mQ32)*pow2(mU32)
         + 120*msq2*mU32*pow2(mQ32)*pow3(m32) - 180*mQ32*msq2*pow2(mU32)*pow3(
         m32) + 381*pow2(mQ32)*pow2(mU32)*pow3(m32) - 150*msq2*mU32*pow2(m32)*
         pow3(mQ32) + 120*m32*msq2*pow2(mU32)*pow3(mQ32) + 35*pow2(m32)*pow2(
         mU32)*pow3(mQ32) + 60*msq2*pow3(m32)*pow3(mQ32) - 1405*mU32*pow3(m32)*
         pow3(mQ32) + 90*mQ32*msq2*pow2(m32)*pow3(mU32) - 60*m32*msq2*pow2(mQ32
         )*pow3(mU32) - 487*pow2(m32)*pow2(mQ32)*pow3(mU32) + 221*mQ32*pow3(m32
         )*pow3(mU32) + 415*m32*pow3(mQ32)*pow3(mU32) - 30*msq2*pow3(mQ32)*pow3
         (mU32) + 90*mQ32*msq2*mU32*pow4(m32) - 90*msq2*pow2(mQ32)*pow4(m32) +
         786*mU32*pow2(mQ32)*pow4(m32) - 338*mQ32*pow2(mU32)*pow4(m32) + 838*
         pow3(mQ32)*pow4(m32) - 6*pow3(mU32)*pow4(m32) - 60*m32*msq2*mU32*pow4(
         mQ32) + 30*msq2*pow2(m32)*pow4(mQ32) + 1029*mU32*pow2(m32)*pow4(mQ32)
         - 345*m32*pow2(mU32)*pow4(mQ32) + 30*msq2*pow2(mU32)*pow4(mQ32) - 477*
         pow3(m32)*pow4(mQ32) - 87*pow3(mU32)*pow4(mQ32) + 9*mQ32*pow2(m32)*
         pow4(mU32) - 6*m32*pow2(mQ32)*pow4(mU32) - 3*pow3(mQ32)*pow4(mU32) -
         38*mQ32*mU32*pow5(m32) - 614*pow2(mQ32)*pow5(m32) + 12*pow2(mU32)*pow5
         (m32) - 201*m32*mU32*pow5(mQ32) + 54*pow2(m32)*pow5(mQ32) + 87*pow2(
         mU32)*pow5(mQ32) + 134*mQ32*pow6(m32) - 6*mU32*pow6(m32) + 9*m32*pow6(
         mQ32) + 3*mU32*pow6(mQ32)))/(3.*(mQ32 - mU32)*pow2(m32 - mU32)*pow4(
         m32 - mQ32)) + (1280*m32*(mQ32 + mU32)*pow2(mQ32)*pow6(Xt))/(3.*pow2(
         m32 - mQ32)*pow5(mQ32 - mU32)) + (4*pow4(Xt)*(-60*msq2*pow2(mQ32)*pow2
         (mU32)*pow3(m32) - 120*msq2*pow2(m32)*pow2(mU32)*pow3(mQ32) + 180*msq2
         *mU32*pow3(m32)*pow3(mQ32) + 6656*pow2(mU32)*pow3(m32)*pow3(mQ32) +
         120*msq2*pow2(m32)*pow2(mQ32)*pow3(mU32) - 180*mQ32*msq2*pow3(m32)*
         pow3(mU32) + 1578*pow2(mQ32)*pow3(m32)*pow3(mU32) + 60*m32*msq2*pow3(
         mQ32)*pow3(mU32) - 1956*pow2(m32)*pow3(mQ32)*pow3(mU32) + 90*mQ32*msq2
         *pow2(mU32)*pow4(m32) - 6012*pow2(mQ32)*pow2(mU32)*pow4(m32) - 90*msq2
         *pow3(mQ32)*pow4(m32) - 8648*mU32*pow3(mQ32)*pow4(m32) - 232*mQ32*pow3
         (mU32)*pow4(m32) - 120*msq2*mU32*pow2(m32)*pow4(mQ32) + 60*m32*msq2*
         pow2(mU32)*pow4(mQ32) - 996*pow2(m32)*pow2(mU32)*pow4(mQ32) + 60*msq2*
         pow3(m32)*pow4(mQ32) + 798*mU32*pow3(m32)*pow4(mQ32) + 502*m32*pow3(
         mU32)*pow4(mQ32) - 250*pow4(m32)*pow4(mQ32) + 90*mQ32*msq2*pow2(m32)*
         pow4(mU32) - 60*m32*msq2*pow2(mQ32)*pow4(mU32) - 222*pow2(m32)*pow2(
         mQ32)*pow4(mU32) - 3*mQ32*pow3(m32)*pow4(mU32) + 265*m32*pow3(mQ32)*
         pow4(mU32) - 30*msq2*pow3(mQ32)*pow4(mU32) + 6*pow4(m32)*pow4(mU32) -
         54*pow4(mQ32)*pow4(mU32) + 8316*mU32*pow2(mQ32)*pow5(m32) + 1454*mQ32*
         pow2(mU32)*pow5(m32) + 3626*pow3(mQ32)*pow5(m32) - 20*pow3(mU32)*pow5(
         m32) - 60*m32*msq2*mU32*pow5(mQ32) + 30*msq2*pow2(m32)*pow5(mQ32) +
         2203*mU32*pow2(m32)*pow5(mQ32) - 1322*m32*pow2(mU32)*pow5(mQ32) + 30*
         msq2*pow2(mU32)*pow5(mQ32) - 965*pow3(m32)*pow5(mQ32) + 144*pow3(mU32)
         *pow5(mQ32) + 9*mQ32*pow2(m32)*pow5(mU32) - 6*m32*pow2(mQ32)*pow5(mU32
         ) - 3*pow3(mQ32)*pow5(mU32) - 2272*mQ32*mU32*pow6(m32) - 3606*pow2(
         mQ32)*pow6(m32) + 22*pow2(mU32)*pow6(m32) - 408*m32*mU32*pow6(mQ32) +
         162*pow2(m32)*pow6(mQ32) + 198*pow2(mU32)*pow6(mQ32) + 1032*mQ32*pow7(
         m32) - 8*mU32*pow7(m32) + 9*m32*pow7(mQ32) + 3*mU32*pow7(mQ32)))/(3.*
         pow2(m32 - mU32)*pow4(m32 - mQ32)*pow4(mQ32 - mU32)) + (32*Xt*(30*m3*
         msq2*pow2(mQ32)*pow2(mU32) - 30*mQ32*msq2*pow2(mU32)*pow3(m3) + 131*
         pow2(mQ32)*pow2(mU32)*pow3(m3) - 30*m3*msq2*mU32*pow3(mQ32) - 70*m3*
         pow2(mU32)*pow3(mQ32) + 30*msq2*pow3(m3)*pow3(mQ32) - 84*mU32*pow3(m3)
         *pow3(mQ32) + 3*m3*pow2(mQ32)*pow3(mU32) - 3*mQ32*pow3(m3)*pow3(mU32)
         + 80*m3*mU32*pow4(mQ32) - 92*pow3(m3)*pow4(mQ32) + 30*mQ32*msq2*mU32*
         pow5(m3) - 30*msq2*pow2(mQ32)*pow5(m3) - 11*mU32*pow2(mQ32)*pow5(m3) -
         113*mQ32*pow2(mU32)*pow5(m3) + 172*pow3(mQ32)*pow5(m3) + 3*m3*pow5(
         mQ32) + 115*mQ32*mU32*pow7(m3) - 135*pow2(mQ32)*pow7(m3) + 4*pow2(mU32
         )*pow7(m3) + 4*mQ32*pow9(m3) - 4*mU32*pow9(m3)))/(3.*(m32 - mU32)*pow2
         (mQ32 - mU32)*pow3(m32 - mQ32)) + (64*pow3(Xt)*(2*pow11(m3) + 30*m3*
         msq2*pow2(mQ32)*pow2(mU32) - 30*mQ32*msq2*pow2(mU32)*pow3(m3) + 122*
         pow2(mQ32)*pow2(mU32)*pow3(m3) - 30*m3*msq2*mU32*pow3(mQ32) - 79*m3*
         pow2(mU32)*pow3(mQ32) + 30*msq2*pow3(m3)*pow3(mQ32) - 58*mU32*pow3(m3)
         *pow3(mQ32) + 3*m3*pow2(mQ32)*pow3(mU32) - 3*mQ32*pow3(m3)*pow3(mU32)
         + 103*m3*mU32*pow4(mQ32) - 115*pow3(m3)*pow4(mQ32) + 30*mQ32*msq2*mU32
         *pow5(m3) - 30*msq2*pow2(mQ32)*pow5(m3) - 147*mU32*pow2(mQ32)*pow5(m3)
         - 28*mQ32*pow2(mU32)*pow5(m3) + 155*pow3(mQ32)*pow5(m3) + 3*m3*pow5(
         mQ32) + 73*mQ32*mU32*pow7(m3) + 10*pow2(mQ32)*pow7(m3) + pow2(mU32)*
         pow7(m3) - 39*mQ32*pow9(m3) - 3*mU32*pow9(m3)))/(3.*(m32 - mU32)*pow3(
         m32 - mQ32)*pow3(mQ32 - mU32)) + (32*pow5(Xt)*(30*msq2*pow2(mQ32)*pow2
         (mU32)*pow3(m3) - 30*msq2*mU32*pow3(m3)*pow3(mQ32) + 129*pow2(mU32)*
         pow3(m3)*pow3(mQ32) - 30*m3*msq2*pow2(mQ32)*pow3(mU32) + 30*mQ32*msq2*
         pow3(m3)*pow3(mU32) - 60*pow2(mQ32)*pow3(m3)*pow3(mU32) + 51*m3*pow3(
         mQ32)*pow3(mU32) + 30*m3*msq2*mU32*pow4(mQ32) - 74*m3*pow2(mU32)*pow4(
         mQ32) - 30*msq2*pow3(m3)*pow4(mQ32) + 268*mU32*pow3(m3)*pow4(mQ32) - 3
         *m3*pow2(mQ32)*pow4(mU32) + 3*mQ32*pow3(m3)*pow4(mU32) - 30*mQ32*msq2*
         pow2(mU32)*pow5(m3) - 72*pow2(mQ32)*pow2(mU32)*pow5(m3) + 30*msq2*pow3
         (mQ32)*pow5(m3) - 201*mU32*pow3(mQ32)*pow5(m3) - 7*mQ32*pow3(mU32)*
         pow5(m3) - 200*pow4(mQ32)*pow5(m3) - 131*m3*mU32*pow5(mQ32) + 140*pow3
         (m3)*pow5(mQ32) - 3*m3*pow6(mQ32) + 112*mU32*pow2(mQ32)*pow7(m3) + 33*
         mQ32*pow2(mU32)*pow7(m3) + 15*pow3(mQ32)*pow7(m3) - 32*mQ32*mU32*pow9(
         m3) + 32*pow2(mQ32)*pow9(m3)))/(3.*(m32 - mU32)*pow3(m32 - mQ32)*pow5(
         mQ32 - mU32))) + pow3(log(mU32/MR2))*((8*pow2(Xt)*(-30*msq2*pow2(m32)*
         pow2(mQ32)*pow2(mU32) + 180*msq2*mU32*pow2(mQ32)*pow3(m32) - 120*mQ32*
         msq2*pow2(mU32)*pow3(m32) - 1405*pow2(mQ32)*pow2(mU32)*pow3(m32) - 90*
         msq2*mU32*pow2(m32)*pow3(mQ32) + 60*m32*msq2*pow2(mU32)*pow3(mQ32) +
         403*pow2(m32)*pow2(mU32)*pow3(mQ32) - 69*mU32*pow3(m32)*pow3(mQ32) +
         150*mQ32*msq2*pow2(m32)*pow3(mU32) - 120*m32*msq2*pow2(mQ32)*pow3(mU32
         ) + 521*pow2(m32)*pow2(mQ32)*pow3(mU32) + 269*mQ32*pow3(m32)*pow3(mU32
         ) - 60*msq2*pow3(m32)*pow3(mU32) - 407*m32*pow3(mQ32)*pow3(mU32) + 30*
         msq2*pow3(mQ32)*pow3(mU32) - 90*mQ32*msq2*mU32*pow4(m32) + 552*mU32*
         pow2(mQ32)*pow4(m32) + 1514*mQ32*pow2(mU32)*pow4(m32) + 90*msq2*pow2(
         mU32)*pow4(m32) - 4*pow3(mQ32)*pow4(m32) - 266*pow3(mU32)*pow4(m32) -
         9*mU32*pow2(m32)*pow4(mQ32) + 6*m32*pow2(mU32)*pow4(mQ32) + 3*pow3(
         mU32)*pow4(mQ32) + 60*m32*mQ32*msq2*pow4(mU32) - 1239*mQ32*pow2(m32)*
         pow4(mU32) - 30*msq2*pow2(m32)*pow4(mU32) + 453*m32*pow2(mQ32)*pow4(
         mU32) - 30*msq2*pow2(mQ32)*pow4(mU32) + 581*pow3(m32)*pow4(mU32) + 85*
         pow3(mQ32)*pow4(mU32) - 846*mQ32*mU32*pow5(m32) + 8*pow2(mQ32)*pow5(
         m32) - 578*pow2(mU32)*pow5(m32) + 333*m32*mQ32*pow5(mU32) - 120*pow2(
         m32)*pow5(mU32) - 153*pow2(mQ32)*pow5(mU32) - 4*mQ32*pow6(m32) + 384*
         mU32*pow6(m32) - 9*m32*pow6(mU32) - 3*mQ32*pow6(mU32)))/(3.*pow2(m32 -
         mQ32)*pow2(-mQ32 + mU32)*pow4(m32 - mU32)) - (4*(30*msq2*pow2(m32)*
         pow2(mQ32)*pow2(mU32) - 180*msq2*mU32*pow2(mQ32)*pow3(m32) + 120*mQ32*
         msq2*pow2(mU32)*pow3(m32) + 365*pow2(mQ32)*pow2(mU32)*pow3(m32) + 90*
         msq2*mU32*pow2(m32)*pow3(mQ32) - 60*m32*msq2*pow2(mU32)*pow3(mQ32) -
         471*pow2(m32)*pow2(mU32)*pow3(mQ32) + 205*mU32*pow3(m32)*pow3(mQ32) -
         150*mQ32*msq2*pow2(m32)*pow3(mU32) + 120*m32*msq2*pow2(mQ32)*pow3(mU32
         ) + 35*pow2(m32)*pow2(mQ32)*pow3(mU32) - 1381*mQ32*pow3(m32)*pow3(mU32
         ) + 60*msq2*pow3(m32)*pow3(mU32) + 407*m32*pow3(mQ32)*pow3(mU32) - 30*
         msq2*pow3(mQ32)*pow3(mU32) + 90*mQ32*msq2*mU32*pow4(m32) - 312*mU32*
         pow2(mQ32)*pow4(m32) + 770*mQ32*pow2(mU32)*pow4(m32) - 90*msq2*pow2(
         mU32)*pow4(m32) + 822*pow3(mU32)*pow4(m32) + 9*mU32*pow2(m32)*pow4(
         mQ32) - 6*m32*pow2(mU32)*pow4(mQ32) - 3*pow3(mU32)*pow4(mQ32) - 60*m32
         *mQ32*msq2*pow4(mU32) + 1015*mQ32*pow2(m32)*pow4(mU32) + 30*msq2*pow2(
         m32)*pow4(mU32) - 341*m32*pow2(mQ32)*pow4(mU32) + 30*msq2*pow2(mQ32)*
         pow4(mU32) - 469*pow3(m32)*pow4(mU32) - 85*pow3(mQ32)*pow4(mU32) - 42*
         mQ32*mU32*pow5(m32) - 598*pow2(mU32)*pow5(m32) - 197*m32*mQ32*pow5(
         mU32) + 52*pow2(m32)*pow5(mU32) + 85*pow2(mQ32)*pow5(mU32) + 128*mU32*
         pow6(m32) + 9*m32*pow6(mU32) + 3*mQ32*pow6(mU32)))/(3.*(-mQ32 + mU32)*
         pow2(m32 - mQ32)*pow4(m32 - mU32)) + (1280*m32*(mQ32 + mU32)*pow2(mU32
         )*pow6(Xt))/(3.*pow2(m32 - mU32)*pow5(-mQ32 + mU32)) + (4*pow4(Xt)*(
         -60*msq2*pow2(mQ32)*pow2(mU32)*pow3(m32) + 120*msq2*pow2(m32)*pow2(
         mU32)*pow3(mQ32) - 180*msq2*mU32*pow3(m32)*pow3(mQ32) + 1634*pow2(mU32
         )*pow3(m32)*pow3(mQ32) - 120*msq2*pow2(m32)*pow2(mQ32)*pow3(mU32) +
         180*mQ32*msq2*pow3(m32)*pow3(mU32) + 6792*pow2(mQ32)*pow3(m32)*pow3(
         mU32) + 60*m32*msq2*pow3(mQ32)*pow3(mU32) - 2020*pow2(m32)*pow3(mQ32)*
         pow3(mU32) + 90*msq2*mU32*pow2(mQ32)*pow4(m32) - 6116*pow2(mQ32)*pow2(
         mU32)*pow4(m32) - 256*mU32*pow3(mQ32)*pow4(m32) - 8728*mQ32*pow3(mU32)
         *pow4(m32) - 90*msq2*pow3(mU32)*pow4(m32) + 90*msq2*mU32*pow2(m32)*
         pow4(mQ32) - 60*m32*msq2*pow2(mU32)*pow4(mQ32) - 222*pow2(m32)*pow2(
         mU32)*pow4(mQ32) - 3*mU32*pow3(m32)*pow4(mQ32) + 265*m32*pow3(mU32)*
         pow4(mQ32) - 30*msq2*pow3(mU32)*pow4(mQ32) + 6*pow4(m32)*pow4(mQ32) -
         120*mQ32*msq2*pow2(m32)*pow4(mU32) + 60*m32*msq2*pow2(mQ32)*pow4(mU32)
         - 1100*pow2(m32)*pow2(mQ32)*pow4(mU32) + 898*mQ32*pow3(m32)*pow4(mU32
         ) + 60*msq2*pow3(m32)*pow4(mU32) + 538*m32*pow3(mQ32)*pow4(mU32) - 282
         *pow4(m32)*pow4(mU32) - 54*pow4(mQ32)*pow4(mU32) + 1498*mU32*pow2(mQ32
         )*pow5(m32) + 8356*mQ32*pow2(mU32)*pow5(m32) - 16*pow3(mQ32)*pow5(m32)
         + 3634*pow3(mU32)*pow5(m32) + 9*mU32*pow2(m32)*pow5(mQ32) - 6*m32*
         pow2(mU32)*pow5(mQ32) - 3*pow3(mU32)*pow5(mQ32) - 60*m32*mQ32*msq2*
         pow5(mU32) + 2139*mQ32*pow2(m32)*pow5(mU32) + 30*msq2*pow2(m32)*pow5(
         mU32) - 1278*m32*pow2(mQ32)*pow5(mU32) + 30*msq2*pow2(mQ32)*pow5(mU32)
         - 937*pow3(m32)*pow5(mU32) + 136*pow3(mQ32)*pow5(mU32) - 2288*mQ32*
         mU32*pow6(m32) + 14*pow2(mQ32)*pow6(m32) - 3598*pow2(mU32)*pow6(m32) -
         392*m32*mQ32*pow6(mU32) + 154*pow2(m32)*pow6(mU32) + 190*pow2(mQ32)*
         pow6(mU32) - 4*mQ32*pow7(m32) + 1028*mU32*pow7(m32) + 9*m32*pow7(mU32)
         + 3*mQ32*pow7(mU32)))/(3.*pow2(m32 - mQ32)*pow4(m32 - mU32)*pow4(
         -mQ32 + mU32)) + (32*pow3(Xt)*(2*pow11(m3) + 60*m3*msq2*pow2(mQ32)*
         pow2(mU32) - 60*msq2*mU32*pow2(mQ32)*pow3(m3) + 241*pow2(mQ32)*pow2(
         mU32)*pow3(m3) + 6*m3*pow2(mU32)*pow3(mQ32) - 6*mU32*pow3(m3)*pow3(
         mQ32) - 60*m3*mQ32*msq2*pow3(mU32) - 157*m3*pow2(mQ32)*pow3(mU32) -
         122*mQ32*pow3(m3)*pow3(mU32) + 60*msq2*pow3(m3)*pow3(mU32) + 207*m3*
         mQ32*pow4(mU32) - 231*pow3(m3)*pow4(mU32) + 60*mQ32*msq2*mU32*pow5(m3)
         - 53*mU32*pow2(mQ32)*pow5(m3) - 282*mQ32*pow2(mU32)*pow5(m3) - 60*
         msq2*pow2(mU32)*pow5(m3) + 315*pow3(mU32)*pow5(m3) + 6*m3*pow5(mU32) +
         136*mQ32*mU32*pow7(m3) + pow2(mQ32)*pow7(m3) + 11*pow2(mU32)*pow7(m3)
         - 3*mQ32*pow9(m3) - 71*mU32*pow9(m3)))/(3.*(m32 - mQ32)*pow3(m32 -
         mU32)*pow3(-mQ32 + mU32)) + (32*Xt*(30*m3*msq2*pow2(mQ32)*pow2(mU32) -
         30*msq2*mU32*pow2(mQ32)*pow3(m3) + 129*pow2(mQ32)*pow2(mU32)*pow3(m3)
         + 3*m3*pow2(mU32)*pow3(mQ32) - 3*mU32*pow3(m3)*pow3(mQ32) - 30*m3*
         mQ32*msq2*pow3(mU32) - 70*m3*pow2(mQ32)*pow3(mU32) - 82*mQ32*pow3(m3)*
         pow3(mU32) + 30*msq2*pow3(m3)*pow3(mU32) + 80*m3*mQ32*pow4(mU32) - 92*
         pow3(m3)*pow4(mU32) + 30*mQ32*msq2*mU32*pow5(m3) - 109*mU32*pow2(mQ32)
         *pow5(m3) - 13*mQ32*pow2(mU32)*pow5(m3) - 30*msq2*pow2(mU32)*pow5(m3)
         + 170*pow3(mU32)*pow5(m3) + 3*m3*pow5(mU32) + 113*mQ32*mU32*pow7(m3) +
         2*pow2(mQ32)*pow7(m3) - 131*pow2(mU32)*pow7(m3) - 2*mQ32*pow9(m3) + 2
         *mU32*pow9(m3)))/(3.*(m32 - mQ32)*pow2(-mQ32 + mU32)*pow3(m32 - mU32))
         + (32*pow5(Xt)*(30*msq2*pow2(mQ32)*pow2(mU32)*pow3(m3) - 30*m3*msq2*
         pow2(mU32)*pow3(mQ32) + 30*msq2*mU32*pow3(m3)*pow3(mQ32) - 60*pow2(
         mU32)*pow3(m3)*pow3(mQ32) - 30*mQ32*msq2*pow3(m3)*pow3(mU32) + 129*
         pow2(mQ32)*pow3(m3)*pow3(mU32) + 51*m3*pow3(mQ32)*pow3(mU32) - 3*m3*
         pow2(mU32)*pow4(mQ32) + 3*mU32*pow3(m3)*pow4(mQ32) + 30*m3*mQ32*msq2*
         pow4(mU32) - 74*m3*pow2(mQ32)*pow4(mU32) + 268*mQ32*pow3(m3)*pow4(mU32
         ) - 30*msq2*pow3(m3)*pow4(mU32) - 30*msq2*mU32*pow2(mQ32)*pow5(m3) -
         72*pow2(mQ32)*pow2(mU32)*pow5(m3) - 7*mU32*pow3(mQ32)*pow5(m3) - 201*
         mQ32*pow3(mU32)*pow5(m3) + 30*msq2*pow3(mU32)*pow5(m3) - 200*pow4(mU32
         )*pow5(m3) - 131*m3*mQ32*pow5(mU32) + 140*pow3(m3)*pow5(mU32) - 3*m3*
         pow6(mU32) + 33*mU32*pow2(mQ32)*pow7(m3) + 112*mQ32*pow2(mU32)*pow7(m3
         ) + 15*pow3(mU32)*pow7(m3) - 32*mQ32*mU32*pow9(m3) + 32*pow2(mU32)*
         pow9(m3)))/(3.*(m32 - mQ32)*pow3(m32 - mU32)*pow5(-mQ32 + mU32))) +
         log(mQ32/MR2)*((-128*(2*m32*mQ32 - pow2(mQ32)))/(3.*pow2(m32 - mQ32))
         + (16*(5*m32*mQ32*mU32 - 7*mQ32*pow2(m32) + 3*m32*pow2(mQ32) - 2*mU32*
         pow2(mQ32) + pow3(m32)))/(3.*(m32 - mU32)*pow2(m32 - mQ32)) - (96*(-(
         mQ32*mU32) + pow2(m32))*(2*m32*mQ32 - pow2(mQ32)))/((m32 - mU32)*pow3(
         m32 - mQ32)) + (16*(3*m32*mQ32 + 3*m32*mU32 + 5*mQ32*mU32 - 11*pow2(
         m32))*(-((7*mQ32 + 6*mU32)*pow2(m32)) - 2*mU32*pow2(mQ32) + m32*(5*
         mQ32*mU32 + 3*pow2(mQ32)) + 7*pow3(m32)))/(3.*pow2(m32 - mU32)*pow3(
         m32 - mQ32)) + pow3(Xt)*((-2048*pow3(m3))/(3.*(m32 - mQ32)*(mQ32 -
         mU32)*mU32) + (128*m3*(3*mQ32 + mU32)*(3*m32*mQ32 + 3*m32*mU32 + 5*
         mQ32*mU32 - 11*pow2(m32)))/(3.*(m32 - mQ32)*(m32 - mU32)*pow3(mQ32 -
         mU32)) + (128*m3*(2*m32*(dilog(1 - m32/mQ32) - dilog(1 - m32/mU32)) +
         mQ32*(-4 - dilog(1 - m32/mQ32) + dilog(1 - m32/mU32)) + mU32*(4 -
         dilog(1 - m32/mQ32) + dilog(1 - m32/mU32)))*(-18*m32*mQ32 + pow2(m32)
         + 9*pow2(mQ32)))/(3.*pow2(m32 - mQ32)*pow3(mQ32 - mU32)) - (48*((64*(
         m3*fin(mQ32,m32,MR2) - m3*fin(mQ32,mU32,MR2) + m3*fin(mU32,m32,MR2) +
         3*pow3(m3) + zt2*pow3(m3)))/(9.*(m32 - mQ32)*(m32 - mU32)) - (4*((20*(
         m3*mQ32 - m3*msq2)*fin(mQ32,msq2,MR2))/((mQ32 - mU32)*pow2(m32 - mQ32)
         ) - (2*m3*(mQ32 - mU32)*fin(mU32,m32,MR2))/((m32 - mU32)*pow2(m32 -
         mQ32)) + (2*m3*(mQ32 - mU32)*fin(mQ32,m32,MR2))/((m32 - mQ32)*pow2(m32
         - mU32)) + (20*(-(m3*msq2) + m3*mU32)*fin(mU32,msq2,MR2))/((-mQ32 +
         mU32)*pow2(m32 - mU32)) + (2*(6*m3*mQ32 + 60*m3*msq2 + 6*m3*mU32 + m3*
         mQ32*zt2 + 10*m3*msq2*zt2 + m3*mU32*zt2 + 24*pow3(m3)))/((m32 - mQ32)*
         (m32 - mU32)) - (20*fin(m32,msq2,MR2)*(m3*mQ32*msq2 - 2*m3*mQ32*mU32 +
         m3*msq2*mU32 + mQ32*pow3(m3) - 2*msq2*pow3(m3) + mU32*pow3(m3)))/(
         pow2(m32 - mQ32)*pow2(m32 - mU32)) + 3*((-2*m3*fin(mQ32,mU32,MR2))/((
         m32 - mQ32)*(m32 - mU32)) + (2*fin(mQ32,m32,MR2)*(m3*mQ32 - 3*m3*mU32
         + 2*pow3(m3)))/((m32 - mQ32)*(m32 - mU32)*(mQ32 - mU32)) - (2*fin(mU32
         ,m32,MR2)*(-3*m3*mQ32 + m3*mU32 + 2*pow3(m3)))/((m32 - mQ32)*(m32 -
         mU32)*(mQ32 - mU32)) + (2*(11*pow3(m3) + 3*zt2*pow3(m3)))/((m32 - mQ32
         )*(m32 - mU32))) + (2*fin(mQ32,mU32,MR2)*(m3*pow2(mQ32) + m3*pow2(mU32
         ) - 2*mQ32*pow3(m3) - 2*mU32*pow3(m3) + 2*pow5(m3)))/(pow2(m32 - mQ32)
         *pow2(m32 - mU32))))/3.))/(mQ32 - mU32)) + (16*(-18*m32*mQ32 + pow2(
         m32) + 9*pow2(mQ32))*(9*m32*(mQ32 + mU32)*pow2(mQ32)*pow2(mU32) - 3*
         pow3(mQ32)*pow3(mU32) + pow3(m32)*(mU32*(19 + 8*dilog(1 - m32/mU32))*
         pow2(mQ32) + mQ32*(19 + 8*dilog(1 - m32/mQ32))*pow2(mU32) + 2*pow3(
         mQ32) + 2*pow3(mU32)) - 4*pow2(m32)*(6*pow2(mQ32)*pow2(mU32) + mU32*(2
         + dilog(1 - m32/mU32))*pow3(mQ32) + mQ32*(2 + dilog(1 - m32/mQ32))*
         pow3(mU32)) - (mQ32*mU32*(13 + 4*dilog(1 - m32/mQ32) + 4*dilog(1 -
         m32/mU32)) + 4*pow2(mQ32) + 4*pow2(mU32))*pow4(m32) + 2*(mQ32 + mU32)*
         pow5(m32)))/(3.*mQ32*mU32*pow2(m32 - mU32)*pow4(m32 - mQ32)) + 6*((8*(
         -(mQ32*mU32) + pow2(m32)))/(3.*(m32 - mQ32)*(m32 - mU32)) + (40*pow2(-
         (mQ32*mU32) + pow2(m32)))/(9.*pow2(m32 - mQ32)*pow2(m32 - mU32)) - 4*(
         -4.111111111111111 - (2*mQ32)/(3.*(-m32 + mQ32)) - (2*mU32)/(3.*(-m32
         + mU32)) + pow2(-2 + (2*mQ32)/(3.*(-m32 + mQ32)) + (2*mU32)/(3.*(-m32
         + mU32))) - (16*(-29.625 - (21*mQ32)/(4.*(-m32 + mQ32)) - (56*mQ32)/(
         mQ32 - mU32) + (56*mU32)/(mQ32 - mU32) + (111*mQ32)/(4.*(-m32 + mU32))
         + (45*mU32)/(2.*(-m32 + mU32)) + zt2/2. - (5*mQ32*zt2)/(2.*(-m32 +
         mQ32)) - (8*mQ32*zt2)/(mQ32 - mU32) + (8*mU32*zt2)/(mQ32 - mU32) + (4*
         mQ32*zt2)/(-m32 + mU32) + (3*mU32*zt2)/(2.*(-m32 + mU32)) + (2*(mQ32 +
         mU32)*fin(mQ32,mU32,MR2))/((m32 - mQ32)*(m32 - mU32)) + (56*pow2(mQ32
         ))/((-m32 + mQ32)*(mQ32 - mU32)) - (111*pow2(mQ32))/(4.*(-m32 + mQ32)*
         (-m32 + mU32)) + (8*zt2*pow2(mQ32))/((-m32 + mQ32)*(mQ32 - mU32)) - (4
         *zt2*pow2(mQ32))/((-m32 + mQ32)*(-m32 + mU32)) + (23*pow2(mQ32))/(8.*
         pow2(-m32 + mQ32)) + (3*zt2*pow2(mQ32))/(2.*pow2(-m32 + mQ32)) - (56*
         pow2(mU32))/((mQ32 - mU32)*(-m32 + mU32)) - (8*zt2*pow2(mU32))/((mQ32
         - mU32)*(-m32 + mU32)) + (23*pow2(mU32))/(8.*pow2(-m32 + mU32)) + (3*
         zt2*pow2(mU32))/(2.*pow2(-m32 + mU32)) + fin(mQ32,m32,MR2)*(-5/(2.*(
         -m32 + mQ32)) - 8/(mQ32 - mU32) + (8*mQ32)/((-m32 + mQ32)*(mQ32 - mU32
         )) + 2/(-m32 + mU32) - (4*mQ32)/((-m32 + mQ32)*(-m32 + mU32)) - (3*
         pow2(mQ32))/pow3(-m32 + mQ32)) - (21*pow3(mQ32))/pow3(-m32 + mQ32) - (
         3*zt2*pow3(mQ32))/pow3(-m32 + mQ32) + fin(mU32,m32,MR2)*(-2/(-m32 +
         mQ32) + 8/(mQ32 - mU32) + 3/(2.*(-m32 + mU32)) - (4*mQ32)/((-m32 +
         mQ32)*(-m32 + mU32)) - (8*mU32)/((mQ32 - mU32)*(-m32 + mU32)) - (3*
         pow2(mU32))/pow3(-m32 + mU32)) - (21*pow3(mU32))/pow3(-m32 + mU32) - (
         3*zt2*pow3(mU32))/pow3(-m32 + mU32)))/9. - (4*(-16.333333333333332 + (
         31*mQ32)/(-m32 + mQ32) + (45*msq2)/(-m32 + mQ32) + (9*mU32)/(2.*(-m32
         + mQ32)) + (9*mQ32)/(2.*(-m32 + mU32)) + (45*msq2)/(-m32 + mU32) + (31
         *mU32)/(-m32 + mU32) - 3*zt2 + (mQ32*zt2)/(4.*(-m32 + mQ32)) + (15*
         msq2*zt2)/(2.*(-m32 + mQ32)) + (3*mU32*zt2)/(4.*(-m32 + mQ32)) + (3*
         mQ32*zt2)/(4.*(-m32 + mU32)) + (15*msq2*zt2)/(2.*(-m32 + mU32)) + (
         mU32*zt2)/(4.*(-m32 + mU32)) - (27275*mQ32*msq2)/(432.*pow2(-m32 +
         mQ32)) - (2825*mQ32*mU32)/(432.*pow2(-m32 + mQ32)) - (65*mQ32*msq2*zt2
         )/(6.*pow2(-m32 + mQ32)) - (mQ32*mU32*zt2)/pow2(-m32 + mQ32) - (23*
         pow2(mQ32))/pow2(-m32 + mQ32) + (1355*msq2*pow2(mQ32))/(432.*(mQ32 -
         msq2)*pow2(-m32 + mQ32)) + (233*mU32*pow2(mQ32))/(432.*(mQ32 - mU32)*
         pow2(-m32 + mQ32)) + (5*msq2*zt2*pow2(mQ32))/(6.*(mQ32 - msq2)*pow2(
         -m32 + mQ32)) - (515*mQ32*pow2(msq2))/(144.*(mQ32 - msq2)*pow2(-m32 +
         mQ32)) - (5*mQ32*zt2*pow2(msq2))/(3.*(mQ32 - msq2)*pow2(-m32 + mQ32))
         + (95*pow2(mQ32)*pow2(msq2))/(216.*pow2(-m32 + mQ32)*pow2(mQ32 - msq2)
         ) + (5*zt2*pow2(mQ32)*pow2(msq2))/(6.*pow2(-m32 + mQ32)*pow2(mQ32 -
         msq2)) - (233*mQ32*pow2(mU32))/(288.*(mQ32 - mU32)*pow2(-m32 + mQ32))
         + (233*pow2(mQ32)*pow2(mU32))/(864.*pow2(-m32 + mQ32)*pow2(mQ32 - mU32
         )) + 3*(16.09722222222222 - (22*mQ32)/(-m32 + mQ32) - (14*mQ32)/(-m32
         + mU32) - (36*mU32)/(-m32 + mU32) + zt2 - (3*mQ32*zt2)/(2.*(-m32 +
         mQ32)) - (2*mQ32*zt2)/(-m32 + mU32) - (7*mU32*zt2)/(2.*(-m32 + mU32))
         - ((mQ32 + mU32)*fin(mQ32,mU32,MR2))/((m32 - mQ32)*(m32 - mU32)) + (14
         *pow2(mQ32))/((-m32 + mQ32)*(-m32 + mU32)) + (2*zt2*pow2(mQ32))/((-m32
         + mQ32)*(-m32 + mU32)) + (fin(mQ32,m32,MR2)*(3*m32*mQ32 - 3*mQ32*mU32
         + pow2(m32) - pow2(mQ32)))/((m32 - mU32)*pow2(-m32 + mQ32)) + (31*
         pow2(mQ32))/pow2(-m32 + mQ32) + (3*zt2*pow2(mQ32))/pow2(-m32 + mQ32) +
         (fin(mU32,m32,MR2)*(3*m32*mU32 - 3*mQ32*mU32 + pow2(m32) - pow2(mU32)
         ))/((m32 - mQ32)*pow2(m32 - mU32)) + (31*pow2(mU32))/pow2(-m32 + mU32)
         + (3*zt2*pow2(mU32))/pow2(-m32 + mU32)) - (6*mQ32*mU32)/pow2(-m32 +
         mU32) - (27275*msq2*mU32)/(432.*pow2(-m32 + mU32)) - (mQ32*mU32*zt2)
         /pow2(-m32 + mU32) - (65*msq2*mU32*zt2)/(6.*pow2(-m32 + mU32)) - (515*
         mU32*pow2(msq2))/(144.*(-msq2 + mU32)*pow2(-m32 + mU32)) - (5*mU32*zt2
         *pow2(msq2))/(3.*(-msq2 + mU32)*pow2(-m32 + mU32)) - (23*pow2(mU32))
         /pow2(-m32 + mU32) + (1355*msq2*pow2(mU32))/(432.*(-msq2 + mU32)*pow2(
         -m32 + mU32)) + (5*msq2*zt2*pow2(mU32))/(6.*(-msq2 + mU32)*pow2(-m32 +
         mU32)) + (95*pow2(msq2)*pow2(mU32))/(216.*pow2(-m32 + mU32)*pow2(
         -msq2 + mU32)) + (5*zt2*pow2(msq2)*pow2(mU32))/(6.*pow2(-m32 + mU32)*
         pow2(-msq2 + mU32)) - (5*fin(mQ32,msq2,MR2)*(-3*m32*mQ32 + 3*m32*msq2
         + mQ32*msq2 - 2*pow2(m32) + pow2(mQ32)))/(2.*pow3(-m32 + mQ32)) + (fin
         (mU32,m32,MR2)*(6*m32*mQ32*mU32 - 8*mQ32*pow2(m32) + 2*mU32*pow2(m32)
         + 3*m32*pow2(mQ32) - 3*m32*pow2(mU32) - mQ32*pow2(mU32) + 2*pow3(m32)
         - pow3(mQ32)))/(4.*(m32 - mU32)*pow3(-m32 + mQ32)) - (95*mQ32*pow3(
         msq2))/(216.*pow2(-m32 + mQ32)*pow2(mQ32 - msq2)) - (5*mQ32*zt2*pow3(
         msq2))/(6.*pow2(-m32 + mQ32)*pow2(mQ32 - msq2)) - (95*mU32*pow3(msq2))
         /(216.*pow2(-m32 + mU32)*pow2(-msq2 + mU32)) - (5*mU32*zt2*pow3(msq2))
         /(6.*pow2(-m32 + mU32)*pow2(-msq2 + mU32)) + (271*pow3(mU32))/(864.*(
         mQ32 - mU32)*pow2(-m32 + mQ32)) + (zt2*pow3(mU32))/(12.*(mQ32 - mU32)*
         pow2(-m32 + mQ32)) - (7*mQ32*pow3(mU32))/(12.*pow2(-m32 + mQ32)*pow2(
         mQ32 - mU32)) - (mQ32*zt2*pow3(mU32))/(12.*pow2(-m32 + mQ32)*pow2(mQ32
         - mU32)) + (fin(mQ32,m32,MR2)*(-6*m32*mQ32*mU32 - 2*mQ32*pow2(m32) +
         8*mU32*pow2(m32) + 3*m32*pow2(mQ32) + mU32*pow2(mQ32) - 3*m32*pow2(
         mU32) - 2*pow3(m32) + pow3(mU32)))/(4.*(m32 - mQ32)*pow3(m32 - mU32))
         + fin(m32,msq2,MR2)*(5/(2.*(-m32 + mQ32)) + 5/(2.*(-m32 + mU32)) + (15
         *mQ32)/(2.*pow2(-m32 + mQ32)) - (15*msq2)/(2.*pow2(-m32 + mQ32)) - (15
         *msq2)/(2.*pow2(-m32 + mU32)) + (15*mU32)/(2.*pow2(-m32 + mU32)) + (10
         *mQ32*msq2)/pow3(-m32 + mQ32) - (10*pow2(mQ32))/pow3(-m32 + mQ32) + (
         10*msq2*mU32)/pow3(-m32 + mU32) - (10*pow2(mU32))/pow3(-m32 + mU32)) -
         (5*fin(mU32,msq2,MR2)*(3*m32*msq2 - 3*m32*mU32 + msq2*mU32 - 2*pow2(
         m32) + pow2(mU32)))/(2.*pow3(-m32 + mU32)) + (271*pow4(mU32))/(864.*
         pow2(-m32 + mQ32)*pow2(mQ32 - mU32)) + (zt2*pow4(mU32))/(12.*pow2(-m32
         + mQ32)*pow2(mQ32 - mU32)) + (fin(mQ32,mU32,MR2)*(15*mU32*pow2(m32)*
         pow2(mQ32) + 15*mQ32*pow2(m32)*pow2(mU32) - 6*m32*pow2(mQ32)*pow2(mU32
         ) - 20*mQ32*mU32*pow3(m32) + 14*pow2(mQ32)*pow3(m32) + 14*pow2(mU32)*
         pow3(m32) - 6*m32*mU32*pow3(mQ32) - 11*pow2(m32)*pow3(mQ32) + pow2(
         mU32)*pow3(mQ32) - 6*m32*mQ32*pow3(mU32) - 11*pow2(m32)*pow3(mU32) +
         pow2(mQ32)*pow3(mU32) - 6*mQ32*pow4(m32) - 6*mU32*pow4(m32) + 3*m32*
         pow4(mQ32) + mU32*pow4(mQ32) + 3*m32*pow4(mU32) + mQ32*pow4(mU32) + 4*
         pow5(m32)))/(4.*pow3(-m32 + mQ32)*pow3(m32 - mU32))))/3.)) + pow2(Xt)*
         ((-64*m32*(-18*m32*mQ32 + pow2(m32) + 9*pow2(mQ32)))/(3.*mQ32*mU32*
         pow2(m32 - mQ32)) + (32*(3*m32*mQ32 + 3*m32*mU32 + 5*mQ32*mU32 - 11*
         pow2(m32))*(3*mQ32*mU32 - 2*m32*(2*mQ32 + 3*mU32) + 7*pow2(m32)))/(3.*
         (mQ32 - mU32)*pow2(m32 - mQ32)*pow2(m32 - mU32)) + (8192*mQ32*(-(mU32*
         dilog(1 - m32/mQ32)) + m32*(dilog(1 - m32/mQ32) - dilog(1 - m32/mU32))
         + mQ32*dilog(1 - m32/mU32))*pow2(m32))/(3.*(m32 - mU32)*pow2(m32 -
         mQ32)*pow2(mQ32 - mU32)) + (12*((8*(-(mQ32*mU32) + pow2(m32)))/(3.*(
         m32 - mQ32)*(m32 - mU32)) + (40*pow2(-(mQ32*mU32) + pow2(m32)))/(9.*
         pow2(m32 - mQ32)*pow2(m32 - mU32)) - 4*(-4.111111111111111 - (2*mQ32)/
         (3.*(-m32 + mQ32)) - (2*mU32)/(3.*(-m32 + mU32)) + pow2(-2 + (2*mQ32)/
         (3.*(-m32 + mQ32)) + (2*mU32)/(3.*(-m32 + mU32))) - (16*(-29.625 - (21
         *mQ32)/(4.*(-m32 + mQ32)) - (56*mQ32)/(mQ32 - mU32) + (56*mU32)/(mQ32
         - mU32) + (111*mQ32)/(4.*(-m32 + mU32)) + (45*mU32)/(2.*(-m32 + mU32))
         + zt2/2. - (5*mQ32*zt2)/(2.*(-m32 + mQ32)) - (8*mQ32*zt2)/(mQ32 -
         mU32) + (8*mU32*zt2)/(mQ32 - mU32) + (4*mQ32*zt2)/(-m32 + mU32) + (3*
         mU32*zt2)/(2.*(-m32 + mU32)) + (2*(mQ32 + mU32)*fin(mQ32,mU32,MR2))/((
         m32 - mQ32)*(m32 - mU32)) + (56*pow2(mQ32))/((-m32 + mQ32)*(mQ32 -
         mU32)) - (111*pow2(mQ32))/(4.*(-m32 + mQ32)*(-m32 + mU32)) + (8*zt2*
         pow2(mQ32))/((-m32 + mQ32)*(mQ32 - mU32)) - (4*zt2*pow2(mQ32))/((-m32
         + mQ32)*(-m32 + mU32)) + (23*pow2(mQ32))/(8.*pow2(-m32 + mQ32)) + (3*
         zt2*pow2(mQ32))/(2.*pow2(-m32 + mQ32)) - (56*pow2(mU32))/((mQ32 - mU32
         )*(-m32 + mU32)) - (8*zt2*pow2(mU32))/((mQ32 - mU32)*(-m32 + mU32)) +
         (23*pow2(mU32))/(8.*pow2(-m32 + mU32)) + (3*zt2*pow2(mU32))/(2.*pow2(
         -m32 + mU32)) + fin(mQ32,m32,MR2)*(-5/(2.*(-m32 + mQ32)) - 8/(mQ32 -
         mU32) + (8*mQ32)/((-m32 + mQ32)*(mQ32 - mU32)) + 2/(-m32 + mU32) - (4*
         mQ32)/((-m32 + mQ32)*(-m32 + mU32)) - (3*pow2(mQ32))/pow3(-m32 + mQ32)
         ) - (21*pow3(mQ32))/pow3(-m32 + mQ32) - (3*zt2*pow3(mQ32))/pow3(-m32 +
         mQ32) + fin(mU32,m32,MR2)*(-2/(-m32 + mQ32) + 8/(mQ32 - mU32) + 3/(2.
         *(-m32 + mU32)) - (4*mQ32)/((-m32 + mQ32)*(-m32 + mU32)) - (8*mU32)/((
         mQ32 - mU32)*(-m32 + mU32)) - (3*pow2(mU32))/pow3(-m32 + mU32)) - (21*
         pow3(mU32))/pow3(-m32 + mU32) - (3*zt2*pow3(mU32))/pow3(-m32 + mU32)))
         /9. - (4*(-16.333333333333332 + (31*mQ32)/(-m32 + mQ32) + (45*msq2)/(
         -m32 + mQ32) + (9*mU32)/(2.*(-m32 + mQ32)) + (9*mQ32)/(2.*(-m32 + mU32
         )) + (45*msq2)/(-m32 + mU32) + (31*mU32)/(-m32 + mU32) - 3*zt2 + (mQ32
         *zt2)/(4.*(-m32 + mQ32)) + (15*msq2*zt2)/(2.*(-m32 + mQ32)) + (3*mU32*
         zt2)/(4.*(-m32 + mQ32)) + (3*mQ32*zt2)/(4.*(-m32 + mU32)) + (15*msq2*
         zt2)/(2.*(-m32 + mU32)) + (mU32*zt2)/(4.*(-m32 + mU32)) - (27275*mQ32*
         msq2)/(432.*pow2(-m32 + mQ32)) - (2825*mQ32*mU32)/(432.*pow2(-m32 +
         mQ32)) - (65*mQ32*msq2*zt2)/(6.*pow2(-m32 + mQ32)) - (mQ32*mU32*zt2)
         /pow2(-m32 + mQ32) - (23*pow2(mQ32))/pow2(-m32 + mQ32) + (1355*msq2*
         pow2(mQ32))/(432.*(mQ32 - msq2)*pow2(-m32 + mQ32)) + (233*mU32*pow2(
         mQ32))/(432.*(mQ32 - mU32)*pow2(-m32 + mQ32)) + (5*msq2*zt2*pow2(mQ32)
         )/(6.*(mQ32 - msq2)*pow2(-m32 + mQ32)) - (515*mQ32*pow2(msq2))/(144.*(
         mQ32 - msq2)*pow2(-m32 + mQ32)) - (5*mQ32*zt2*pow2(msq2))/(3.*(mQ32 -
         msq2)*pow2(-m32 + mQ32)) + (95*pow2(mQ32)*pow2(msq2))/(216.*pow2(-m32
         + mQ32)*pow2(mQ32 - msq2)) + (5*zt2*pow2(mQ32)*pow2(msq2))/(6.*pow2(
         -m32 + mQ32)*pow2(mQ32 - msq2)) - (233*mQ32*pow2(mU32))/(288.*(mQ32 -
         mU32)*pow2(-m32 + mQ32)) + (233*pow2(mQ32)*pow2(mU32))/(864.*pow2(-m32
         + mQ32)*pow2(mQ32 - mU32)) + 3*(16.09722222222222 - (22*mQ32)/(-m32 +
         mQ32) - (14*mQ32)/(-m32 + mU32) - (36*mU32)/(-m32 + mU32) + zt2 - (3*
         mQ32*zt2)/(2.*(-m32 + mQ32)) - (2*mQ32*zt2)/(-m32 + mU32) - (7*mU32*
         zt2)/(2.*(-m32 + mU32)) - ((mQ32 + mU32)*fin(mQ32,mU32,MR2))/((m32 -
         mQ32)*(m32 - mU32)) + (14*pow2(mQ32))/((-m32 + mQ32)*(-m32 + mU32)) +
         (2*zt2*pow2(mQ32))/((-m32 + mQ32)*(-m32 + mU32)) + (fin(mQ32,m32,MR2)*
         (3*m32*mQ32 - 3*mQ32*mU32 + pow2(m32) - pow2(mQ32)))/((m32 - mU32)*
         pow2(-m32 + mQ32)) + (31*pow2(mQ32))/pow2(-m32 + mQ32) + (3*zt2*pow2(
         mQ32))/pow2(-m32 + mQ32) + (fin(mU32,m32,MR2)*(3*m32*mU32 - 3*mQ32*
         mU32 + pow2(m32) - pow2(mU32)))/((m32 - mQ32)*pow2(m32 - mU32)) + (31*
         pow2(mU32))/pow2(-m32 + mU32) + (3*zt2*pow2(mU32))/pow2(-m32 + mU32))
         - (6*mQ32*mU32)/pow2(-m32 + mU32) - (27275*msq2*mU32)/(432.*pow2(-m32
         + mU32)) - (mQ32*mU32*zt2)/pow2(-m32 + mU32) - (65*msq2*mU32*zt2)/(6.*
         pow2(-m32 + mU32)) - (515*mU32*pow2(msq2))/(144.*(-msq2 + mU32)*pow2(
         -m32 + mU32)) - (5*mU32*zt2*pow2(msq2))/(3.*(-msq2 + mU32)*pow2(-m32 +
         mU32)) - (23*pow2(mU32))/pow2(-m32 + mU32) + (1355*msq2*pow2(mU32))/(
         432.*(-msq2 + mU32)*pow2(-m32 + mU32)) + (5*msq2*zt2*pow2(mU32))/(6.*(
         -msq2 + mU32)*pow2(-m32 + mU32)) + (95*pow2(msq2)*pow2(mU32))/(216.*
         pow2(-m32 + mU32)*pow2(-msq2 + mU32)) + (5*zt2*pow2(msq2)*pow2(mU32))/
         (6.*pow2(-m32 + mU32)*pow2(-msq2 + mU32)) - (5*fin(mQ32,msq2,MR2)*(-3*
         m32*mQ32 + 3*m32*msq2 + mQ32*msq2 - 2*pow2(m32) + pow2(mQ32)))/(2.*
         pow3(-m32 + mQ32)) + (fin(mU32,m32,MR2)*(6*m32*mQ32*mU32 - 8*mQ32*pow2
         (m32) + 2*mU32*pow2(m32) + 3*m32*pow2(mQ32) - 3*m32*pow2(mU32) - mQ32*
         pow2(mU32) + 2*pow3(m32) - pow3(mQ32)))/(4.*(m32 - mU32)*pow3(-m32 +
         mQ32)) - (95*mQ32*pow3(msq2))/(216.*pow2(-m32 + mQ32)*pow2(mQ32 - msq2
         )) - (5*mQ32*zt2*pow3(msq2))/(6.*pow2(-m32 + mQ32)*pow2(mQ32 - msq2))
         - (95*mU32*pow3(msq2))/(216.*pow2(-m32 + mU32)*pow2(-msq2 + mU32)) - (
         5*mU32*zt2*pow3(msq2))/(6.*pow2(-m32 + mU32)*pow2(-msq2 + mU32)) + (
         271*pow3(mU32))/(864.*(mQ32 - mU32)*pow2(-m32 + mQ32)) + (zt2*pow3(
         mU32))/(12.*(mQ32 - mU32)*pow2(-m32 + mQ32)) - (7*mQ32*pow3(mU32))/(
         12.*pow2(-m32 + mQ32)*pow2(mQ32 - mU32)) - (mQ32*zt2*pow3(mU32))/(12.*
         pow2(-m32 + mQ32)*pow2(mQ32 - mU32)) + (fin(mQ32,m32,MR2)*(-6*m32*mQ32
         *mU32 - 2*mQ32*pow2(m32) + 8*mU32*pow2(m32) + 3*m32*pow2(mQ32) + mU32*
         pow2(mQ32) - 3*m32*pow2(mU32) - 2*pow3(m32) + pow3(mU32)))/(4.*(m32 -
         mQ32)*pow3(m32 - mU32)) + fin(m32,msq2,MR2)*(5/(2.*(-m32 + mQ32)) + 5/
         (2.*(-m32 + mU32)) + (15*mQ32)/(2.*pow2(-m32 + mQ32)) - (15*msq2)/(2.*
         pow2(-m32 + mQ32)) - (15*msq2)/(2.*pow2(-m32 + mU32)) + (15*mU32)/(2.*
         pow2(-m32 + mU32)) + (10*mQ32*msq2)/pow3(-m32 + mQ32) - (10*pow2(mQ32)
         )/pow3(-m32 + mQ32) + (10*msq2*mU32)/pow3(-m32 + mU32) - (10*pow2(mU32
         ))/pow3(-m32 + mU32)) - (5*fin(mU32,msq2,MR2)*(3*m32*msq2 - 3*m32*mU32
         + msq2*mU32 - 2*pow2(m32) + pow2(mU32)))/(2.*pow3(-m32 + mU32)) + (
         271*pow4(mU32))/(864.*pow2(-m32 + mQ32)*pow2(mQ32 - mU32)) + (zt2*pow4
         (mU32))/(12.*pow2(-m32 + mQ32)*pow2(mQ32 - mU32)) + (fin(mQ32,mU32,MR2
         )*(15*mU32*pow2(m32)*pow2(mQ32) + 15*mQ32*pow2(m32)*pow2(mU32) - 6*m32
         *pow2(mQ32)*pow2(mU32) - 20*mQ32*mU32*pow3(m32) + 14*pow2(mQ32)*pow3(
         m32) + 14*pow2(mU32)*pow3(m32) - 6*m32*mU32*pow3(mQ32) - 11*pow2(m32)*
         pow3(mQ32) + pow2(mU32)*pow3(mQ32) - 6*m32*mQ32*pow3(mU32) - 11*pow2(
         m32)*pow3(mU32) + pow2(mQ32)*pow3(mU32) - 6*mQ32*pow4(m32) - 6*mU32*
         pow4(m32) + 3*m32*pow4(mQ32) + mU32*pow4(mQ32) + 3*m32*pow4(mU32) +
         mQ32*pow4(mU32) + 4*pow5(m32)))/(4.*pow3(-m32 + mQ32)*pow3(m32 - mU32)
         )))/3.)))/(mQ32 - mU32)) + (8*(-420*msq2*mU32*pow2(m32)*pow2(mQ32) +
         360*mQ32*msq2*pow2(m32)*pow2(mU32) + 300*m32*msq2*pow2(mQ32)*pow2(mU32
         ) - 747*pow2(m32)*pow2(mQ32)*pow2(mU32) - 180*mQ32*msq2*mU32*pow3(m32)
         + 180*msq2*pow2(mQ32)*pow3(m32) - 215*mU32*pow2(mQ32)*pow3(m32) + 20*
         mQ32*pow2(mU32)*pow3(m32) - 120*m32*msq2*mU32*pow3(mQ32) + 60*msq2*
         pow2(m32)*pow3(mQ32) + 1206*mU32*pow2(m32)*pow3(mQ32) - 254*m32*pow2(
         mU32)*pow3(mQ32) + 60*msq2*pow2(mU32)*pow3(mQ32) - 510*pow3(m32)*pow3(
         mQ32) - 180*m32*mQ32*msq2*pow3(mU32) - 4*mQ32*pow2(m32)*pow3(mU32) +
         505*m32*pow2(mQ32)*pow3(mU32) - 60*msq2*pow2(mQ32)*pow3(mU32) - 63*
         pow3(m32)*pow3(mU32) - 154*pow3(mQ32)*pow3(mU32) + 116*mQ32*mU32*pow4(
         m32) + 271*pow2(mQ32)*pow4(m32) + 125*pow2(mU32)*pow4(m32) - 397*m32*
         mU32*pow4(mQ32) + 57*pow2(m32)*pow4(mQ32) + 148*pow2(mU32)*pow4(mQ32)
         - 18*m32*mQ32*pow4(mU32) - 6*pow2(mQ32)*pow4(mU32) - 66*mQ32*pow5(m32)
         - 62*mU32*pow5(m32) + 36*m32*pow5(mQ32) + 12*mU32*pow5(mQ32)))/(3.*(
         mQ32 - mU32)*pow2(m32 - mU32)*pow3(m32 - mQ32)) + pow4(Xt)*((32*(-18*
         m32*mQ32 + pow2(m32) + 9*pow2(mQ32))*((2*mQ32*mU32*(-dilog(1 -
         m32/mQ32) + dilog(1 - m32/mU32)) - pow2(mQ32) + pow2(mU32))*pow3(m32)
         + (-5 + dilog(1 - m32/mQ32) - dilog(1 - m32/mU32))*pow2(mU32)*pow3(
         mQ32) + pow2(m32)*(3*mU32*(-2 + dilog(1 - m32/mQ32) - dilog(1 -
         m32/mU32))*pow2(mQ32) + 3*mQ32*(2 + dilog(1 - m32/mQ32) - dilog(1 -
         m32/mU32))*pow2(mU32) + pow3(mQ32) - pow3(mU32)) + (5 + dilog(1 -
         m32/mQ32) - dilog(1 - m32/mU32))*pow2(mQ32)*pow3(mU32) + m32*(4*(
         -dilog(1 - m32/mQ32) + dilog(1 - m32/mU32))*pow2(mQ32)*pow2(mU32) +
         mU32*(5 - dilog(1 - m32/mQ32) + dilog(1 - m32/mU32))*pow3(mQ32) + mQ32
         *(-5 - dilog(1 - m32/mQ32) + dilog(1 - m32/mU32))*pow3(mU32))))/(3.*
         mQ32*(-m32 + mQ32)*(m32 - mU32)*mU32*pow2(m32 - mQ32)*pow3(mQ32 - mU32
         )) - (16*(3*m32*mQ32 + 3*m32*mU32 + 5*mQ32*mU32 - 11*pow2(m32))*(-3*
         mU32*(5*mQ32 + mU32)*pow2(mQ32) - pow2(m32)*(15*mQ32*mU32 + 29*pow2(
         mQ32) + 6*pow2(mU32)) + (7*mQ32 + 3*mU32)*pow3(m32) + m32*(31*mU32*
         pow2(mQ32) + 7*mQ32*pow2(mU32) + 16*pow3(mQ32)) + 4*pow4(m32)))/(3.*
         pow2(m32 - mU32)*pow3(m32 - mQ32)*pow3(mQ32 - mU32)) + (4096*m32*mQ32*
         (2*m32*(dilog(1 - m32/mQ32) - dilog(1 - m32/mU32)) + mQ32*(-4 - dilog(
         1 - m32/mQ32) + dilog(1 - m32/mU32)) + mU32*(4 - dilog(1 - m32/mQ32) +
         dilog(1 - m32/mU32))))/(3.*(m32 - mQ32)*pow4(mQ32 - mU32)) + 2*((-3*(
         mQ32 + mU32)*((8*(-(mQ32*mU32) + pow2(m32)))/(3.*(m32 - mQ32)*(m32 -
         mU32)) + (40*pow2(-(mQ32*mU32) + pow2(m32)))/(9.*pow2(m32 - mQ32)*pow2
         (m32 - mU32)) - 4*(-4.111111111111111 - (2*mQ32)/(3.*(-m32 + mQ32)) -
         (2*mU32)/(3.*(-m32 + mU32)) + pow2(-2 + (2*mQ32)/(3.*(-m32 + mQ32)) +
         (2*mU32)/(3.*(-m32 + mU32))) - (16*(-29.625 - (21*mQ32)/(4.*(-m32 +
         mQ32)) - (56*mQ32)/(mQ32 - mU32) + (56*mU32)/(mQ32 - mU32) + (111*mQ32
         )/(4.*(-m32 + mU32)) + (45*mU32)/(2.*(-m32 + mU32)) + zt2/2. - (5*mQ32
         *zt2)/(2.*(-m32 + mQ32)) - (8*mQ32*zt2)/(mQ32 - mU32) + (8*mU32*zt2)/(
         mQ32 - mU32) + (4*mQ32*zt2)/(-m32 + mU32) + (3*mU32*zt2)/(2.*(-m32 +
         mU32)) + (2*(mQ32 + mU32)*fin(mQ32,mU32,MR2))/((m32 - mQ32)*(m32 -
         mU32)) + (56*pow2(mQ32))/((-m32 + mQ32)*(mQ32 - mU32)) - (111*pow2(
         mQ32))/(4.*(-m32 + mQ32)*(-m32 + mU32)) + (8*zt2*pow2(mQ32))/((-m32 +
         mQ32)*(mQ32 - mU32)) - (4*zt2*pow2(mQ32))/((-m32 + mQ32)*(-m32 + mU32)
         ) + (23*pow2(mQ32))/(8.*pow2(-m32 + mQ32)) + (3*zt2*pow2(mQ32))/(2.*
         pow2(-m32 + mQ32)) - (56*pow2(mU32))/((mQ32 - mU32)*(-m32 + mU32)) - (
         8*zt2*pow2(mU32))/((mQ32 - mU32)*(-m32 + mU32)) + (23*pow2(mU32))/(8.*
         pow2(-m32 + mU32)) + (3*zt2*pow2(mU32))/(2.*pow2(-m32 + mU32)) + fin(
         mQ32,m32,MR2)*(-5/(2.*(-m32 + mQ32)) - 8/(mQ32 - mU32) + (8*mQ32)/((
         -m32 + mQ32)*(mQ32 - mU32)) + 2/(-m32 + mU32) - (4*mQ32)/((-m32 + mQ32
         )*(-m32 + mU32)) - (3*pow2(mQ32))/pow3(-m32 + mQ32)) - (21*pow3(mQ32))
         /pow3(-m32 + mQ32) - (3*zt2*pow3(mQ32))/pow3(-m32 + mQ32) + fin(mU32,
         m32,MR2)*(-2/(-m32 + mQ32) + 8/(mQ32 - mU32) + 3/(2.*(-m32 + mU32)) -
         (4*mQ32)/((-m32 + mQ32)*(-m32 + mU32)) - (8*mU32)/((mQ32 - mU32)*(-m32
         + mU32)) - (3*pow2(mU32))/pow3(-m32 + mU32)) - (21*pow3(mU32))/pow3(
         -m32 + mU32) - (3*zt2*pow3(mU32))/pow3(-m32 + mU32)))/9. - (4*(
         -16.333333333333332 + (31*mQ32)/(-m32 + mQ32) + (45*msq2)/(-m32 + mQ32
         ) + (9*mU32)/(2.*(-m32 + mQ32)) + (9*mQ32)/(2.*(-m32 + mU32)) + (45*
         msq2)/(-m32 + mU32) + (31*mU32)/(-m32 + mU32) - 3*zt2 + (mQ32*zt2)/(4.
         *(-m32 + mQ32)) + (15*msq2*zt2)/(2.*(-m32 + mQ32)) + (3*mU32*zt2)/(4.*
         (-m32 + mQ32)) + (3*mQ32*zt2)/(4.*(-m32 + mU32)) + (15*msq2*zt2)/(2.*(
         -m32 + mU32)) + (mU32*zt2)/(4.*(-m32 + mU32)) - (27275*mQ32*msq2)/(
         432.*pow2(-m32 + mQ32)) - (2825*mQ32*mU32)/(432.*pow2(-m32 + mQ32)) -
         (65*mQ32*msq2*zt2)/(6.*pow2(-m32 + mQ32)) - (mQ32*mU32*zt2)/pow2(-m32
         + mQ32) - (23*pow2(mQ32))/pow2(-m32 + mQ32) + (1355*msq2*pow2(mQ32))/(
         432.*(mQ32 - msq2)*pow2(-m32 + mQ32)) + (233*mU32*pow2(mQ32))/(432.*(
         mQ32 - mU32)*pow2(-m32 + mQ32)) + (5*msq2*zt2*pow2(mQ32))/(6.*(mQ32 -
         msq2)*pow2(-m32 + mQ32)) - (515*mQ32*pow2(msq2))/(144.*(mQ32 - msq2)*
         pow2(-m32 + mQ32)) - (5*mQ32*zt2*pow2(msq2))/(3.*(mQ32 - msq2)*pow2(
         -m32 + mQ32)) + (95*pow2(mQ32)*pow2(msq2))/(216.*pow2(-m32 + mQ32)*
         pow2(mQ32 - msq2)) + (5*zt2*pow2(mQ32)*pow2(msq2))/(6.*pow2(-m32 +
         mQ32)*pow2(mQ32 - msq2)) - (233*mQ32*pow2(mU32))/(288.*(mQ32 - mU32)*
         pow2(-m32 + mQ32)) + (233*pow2(mQ32)*pow2(mU32))/(864.*pow2(-m32 +
         mQ32)*pow2(mQ32 - mU32)) + 3*(16.09722222222222 - (22*mQ32)/(-m32 +
         mQ32) - (14*mQ32)/(-m32 + mU32) - (36*mU32)/(-m32 + mU32) + zt2 - (3*
         mQ32*zt2)/(2.*(-m32 + mQ32)) - (2*mQ32*zt2)/(-m32 + mU32) - (7*mU32*
         zt2)/(2.*(-m32 + mU32)) - ((mQ32 + mU32)*fin(mQ32,mU32,MR2))/((m32 -
         mQ32)*(m32 - mU32)) + (14*pow2(mQ32))/((-m32 + mQ32)*(-m32 + mU32)) +
         (2*zt2*pow2(mQ32))/((-m32 + mQ32)*(-m32 + mU32)) + (fin(mQ32,m32,MR2)*
         (3*m32*mQ32 - 3*mQ32*mU32 + pow2(m32) - pow2(mQ32)))/((m32 - mU32)*
         pow2(-m32 + mQ32)) + (31*pow2(mQ32))/pow2(-m32 + mQ32) + (3*zt2*pow2(
         mQ32))/pow2(-m32 + mQ32) + (fin(mU32,m32,MR2)*(3*m32*mU32 - 3*mQ32*
         mU32 + pow2(m32) - pow2(mU32)))/((m32 - mQ32)*pow2(m32 - mU32)) + (31*
         pow2(mU32))/pow2(-m32 + mU32) + (3*zt2*pow2(mU32))/pow2(-m32 + mU32))
         - (6*mQ32*mU32)/pow2(-m32 + mU32) - (27275*msq2*mU32)/(432.*pow2(-m32
         + mU32)) - (mQ32*mU32*zt2)/pow2(-m32 + mU32) - (65*msq2*mU32*zt2)/(6.*
         pow2(-m32 + mU32)) - (515*mU32*pow2(msq2))/(144.*(-msq2 + mU32)*pow2(
         -m32 + mU32)) - (5*mU32*zt2*pow2(msq2))/(3.*(-msq2 + mU32)*pow2(-m32 +
         mU32)) - (23*pow2(mU32))/pow2(-m32 + mU32) + (1355*msq2*pow2(mU32))/(
         432.*(-msq2 + mU32)*pow2(-m32 + mU32)) + (5*msq2*zt2*pow2(mU32))/(6.*(
         -msq2 + mU32)*pow2(-m32 + mU32)) + (95*pow2(msq2)*pow2(mU32))/(216.*
         pow2(-m32 + mU32)*pow2(-msq2 + mU32)) + (5*zt2*pow2(msq2)*pow2(mU32))/
         (6.*pow2(-m32 + mU32)*pow2(-msq2 + mU32)) - (5*fin(mQ32,msq2,MR2)*(-3*
         m32*mQ32 + 3*m32*msq2 + mQ32*msq2 - 2*pow2(m32) + pow2(mQ32)))/(2.*
         pow3(-m32 + mQ32)) + (fin(mU32,m32,MR2)*(6*m32*mQ32*mU32 - 8*mQ32*pow2
         (m32) + 2*mU32*pow2(m32) + 3*m32*pow2(mQ32) - 3*m32*pow2(mU32) - mQ32*
         pow2(mU32) + 2*pow3(m32) - pow3(mQ32)))/(4.*(m32 - mU32)*pow3(-m32 +
         mQ32)) - (95*mQ32*pow3(msq2))/(216.*pow2(-m32 + mQ32)*pow2(mQ32 - msq2
         )) - (5*mQ32*zt2*pow3(msq2))/(6.*pow2(-m32 + mQ32)*pow2(mQ32 - msq2))
         - (95*mU32*pow3(msq2))/(216.*pow2(-m32 + mU32)*pow2(-msq2 + mU32)) - (
         5*mU32*zt2*pow3(msq2))/(6.*pow2(-m32 + mU32)*pow2(-msq2 + mU32)) + (
         271*pow3(mU32))/(864.*(mQ32 - mU32)*pow2(-m32 + mQ32)) + (zt2*pow3(
         mU32))/(12.*(mQ32 - mU32)*pow2(-m32 + mQ32)) - (7*mQ32*pow3(mU32))/(
         12.*pow2(-m32 + mQ32)*pow2(mQ32 - mU32)) - (mQ32*zt2*pow3(mU32))/(12.*
         pow2(-m32 + mQ32)*pow2(mQ32 - mU32)) + (fin(mQ32,m32,MR2)*(-6*m32*mQ32
         *mU32 - 2*mQ32*pow2(m32) + 8*mU32*pow2(m32) + 3*m32*pow2(mQ32) + mU32*
         pow2(mQ32) - 3*m32*pow2(mU32) - 2*pow3(m32) + pow3(mU32)))/(4.*(m32 -
         mQ32)*pow3(m32 - mU32)) + fin(m32,msq2,MR2)*(5/(2.*(-m32 + mQ32)) + 5/
         (2.*(-m32 + mU32)) + (15*mQ32)/(2.*pow2(-m32 + mQ32)) - (15*msq2)/(2.*
         pow2(-m32 + mQ32)) - (15*msq2)/(2.*pow2(-m32 + mU32)) + (15*mU32)/(2.*
         pow2(-m32 + mU32)) + (10*mQ32*msq2)/pow3(-m32 + mQ32) - (10*pow2(mQ32)
         )/pow3(-m32 + mQ32) + (10*msq2*mU32)/pow3(-m32 + mU32) - (10*pow2(mU32
         ))/pow3(-m32 + mU32)) - (5*fin(mU32,msq2,MR2)*(3*m32*msq2 - 3*m32*mU32
         + msq2*mU32 - 2*pow2(m32) + pow2(mU32)))/(2.*pow3(-m32 + mU32)) + (
         271*pow4(mU32))/(864.*pow2(-m32 + mQ32)*pow2(mQ32 - mU32)) + (zt2*pow4
         (mU32))/(12.*pow2(-m32 + mQ32)*pow2(mQ32 - mU32)) + (fin(mQ32,mU32,MR2
         )*(15*mU32*pow2(m32)*pow2(mQ32) + 15*mQ32*pow2(m32)*pow2(mU32) - 6*m32
         *pow2(mQ32)*pow2(mU32) - 20*mQ32*mU32*pow3(m32) + 14*pow2(mQ32)*pow3(
         m32) + 14*pow2(mU32)*pow3(m32) - 6*m32*mU32*pow3(mQ32) - 11*pow2(m32)*
         pow3(mQ32) + pow2(mU32)*pow3(mQ32) - 6*m32*mQ32*pow3(mU32) - 11*pow2(
         m32)*pow3(mU32) + pow2(mQ32)*pow3(mU32) - 6*mQ32*pow4(m32) - 6*mU32*
         pow4(m32) + 3*m32*pow4(mQ32) + mU32*pow4(mQ32) + 3*m32*pow4(mU32) +
         mQ32*pow4(mU32) + 4*pow5(m32)))/(4.*pow3(-m32 + mQ32)*pow3(m32 - mU32)
         )))/3.)))/pow3(mQ32 - mU32) + (8*(420*msq2*mU32*pow2(m32)*pow2(mQ32) -
         360*mQ32*msq2*pow2(m32)*pow2(mU32) - 300*m32*msq2*pow2(mQ32)*pow2(
         mU32) + 743*pow2(m32)*pow2(mQ32)*pow2(mU32) + 180*mQ32*msq2*mU32*pow3(
         m32) - 180*msq2*pow2(mQ32)*pow3(m32) + 189*mU32*pow2(mQ32)*pow3(m32) +
         46*mQ32*pow2(mU32)*pow3(m32) + 120*m32*msq2*mU32*pow3(mQ32) - 60*msq2
         *pow2(m32)*pow3(mQ32) - 1198*mU32*pow2(m32)*pow3(mQ32) + 270*m32*pow2(
         mU32)*pow3(mQ32) - 60*msq2*pow2(mU32)*pow3(mQ32) + 470*pow3(m32)*pow3(
         mQ32) + 180*m32*mQ32*msq2*pow3(mU32) - 6*mQ32*pow2(m32)*pow3(mU32) -
         531*m32*pow2(mQ32)*pow3(mU32) + 60*msq2*pow2(mQ32)*pow3(mU32) + 63*
         pow3(m32)*pow3(mU32) + 170*pow3(mQ32)*pow3(mU32) - 170*mQ32*mU32*pow4(
         m32) - 215*pow2(mQ32)*pow4(m32) - 127*pow2(mU32)*pow4(m32) + 407*m32*
         mU32*pow4(mQ32) - 51*pow2(m32)*pow4(mQ32) - 164*pow2(mU32)*pow4(mQ32)
         + 18*m32*mQ32*pow4(mU32) + 6*pow2(mQ32)*pow4(mU32) + 64*mQ32*pow5(m32)
         + 64*mU32*pow5(m32) - 36*m32*pow5(mQ32) - 12*mU32*pow5(mQ32)))/(3.*
         pow2(m32 - mU32)*pow3(m32 - mQ32)*pow3(mQ32 - mU32)))) + pow2(log(
         msq2/MR2))*((320*m3*msq2*Xt)/((m32 - mQ32)*(m32 - mU32)) + (640*m3*
         msq2*pow3(Xt))/((m32 - mQ32)*(m32 - mU32)*(mQ32 - mU32)) - (40*(-4*m32
         *mQ32*msq2*mU32 - 5*mQ32*msq2*pow2(m32) - 8*mQ32*mU32*pow2(m32) - 5*
         msq2*mU32*pow2(m32) + 3*m32*msq2*pow2(mQ32) + 4*m32*mU32*pow2(mQ32) +
         msq2*mU32*pow2(mQ32) - 2*pow2(m32)*pow2(mQ32) + 4*m32*mQ32*pow2(mU32)
         + 3*m32*msq2*pow2(mU32) + mQ32*msq2*pow2(mU32) - 2*pow2(m32)*pow2(mU32
         ) - 2*pow2(mQ32)*pow2(mU32) + 4*mQ32*pow3(m32) + 6*msq2*pow3(m32) + 4*
         mU32*pow3(m32) - 2*pow4(m32)))/(pow2(-m32 + mQ32)*pow2(m32 - mU32)) -
         (80*pow2(Xt)*(-4*m32*mQ32*msq2*mU32 - 5*mQ32*msq2*pow2(m32) - 8*mQ32*
         mU32*pow2(m32) - 5*msq2*mU32*pow2(m32) + 3*m32*msq2*pow2(mQ32) + 4*m32
         *mU32*pow2(mQ32) + msq2*mU32*pow2(mQ32) - 2*pow2(m32)*pow2(mQ32) + 4*
         m32*mQ32*pow2(mU32) + 3*m32*msq2*pow2(mU32) + mQ32*msq2*pow2(mU32) - 2
         *pow2(m32)*pow2(mU32) - 2*pow2(mQ32)*pow2(mU32) + 4*mQ32*pow3(m32) + 6
         *msq2*pow3(m32) + 4*mU32*pow3(m32) - 2*pow4(m32)))/((mQ32 - mU32)*pow2
         (-m32 + mQ32)*pow2(m32 - mU32)) + (40*(mQ32 + mU32)*(-4*m32*mQ32*msq2*
         mU32 - 5*mQ32*msq2*pow2(m32) - 8*mQ32*mU32*pow2(m32) - 5*msq2*mU32*
         pow2(m32) + 3*m32*msq2*pow2(mQ32) + 4*m32*mU32*pow2(mQ32) + msq2*mU32*
         pow2(mQ32) - 2*pow2(m32)*pow2(mQ32) + 4*m32*mQ32*pow2(mU32) + 3*m32*
         msq2*pow2(mU32) + mQ32*msq2*pow2(mU32) - 2*pow2(m32)*pow2(mU32) - 2*
         pow2(mQ32)*pow2(mU32) + 4*mQ32*pow3(m32) + 6*msq2*pow3(m32) + 4*mU32*
         pow3(m32) - 2*pow4(m32))*pow4(Xt))/(pow2(-m32 + mQ32)*pow2(m32 - mU32)
         *pow3(mQ32 - mU32)) - (320*m3*msq2*(mQ32 + mU32)*pow5(Xt))/((m32 -
         mQ32)*(m32 - mU32)*pow3(mQ32 - mU32))) + log(msq2/MR2)*(Xt*((-2560*m3*
         msq2)/((m32 - mQ32)*(m32 - mU32)) - (192*m3)/(mQ32 - mU32) - (96*m3*
         mQ32)/((m32 - mQ32)*(mQ32 - mU32)) + (320*m3*mQ32*(m32 - mQ32 + 6*msq2
         ))/(3.*(mQ32 - mU32)*pow2(m32 - mQ32))) + ((-5120*m3*msq2)/((m32 -
         mQ32)*(m32 - mU32)*(mQ32 - mU32)) + (192*m3*(3*mQ32 + mU32))/pow3(mQ32
         - mU32))*pow3(Xt) + (32*pow2(Xt)*(-240*m32*mQ32*msq2*mU32 - 300*mQ32*
         msq2*pow2(m32) + 330*mQ32*mU32*pow2(m32) - 300*msq2*mU32*pow2(m32) +
         180*m32*msq2*pow2(mQ32) - 151*m32*mU32*pow2(mQ32) + 60*msq2*mU32*pow2(
         mQ32) + 78*pow2(m32)*pow2(mQ32) - 160*m32*mQ32*pow2(mU32) + 180*m32*
         msq2*pow2(mU32) + 60*mQ32*msq2*pow2(mU32) + 87*pow2(m32)*pow2(mU32) +
         73*pow2(mQ32)*pow2(mU32) - 170*mQ32*pow3(m32) + 360*msq2*pow3(m32) -
         179*mU32*pow3(m32) + 92*pow4(m32)))/(3.*(mQ32 - mU32)*pow2(m32 - mQ32)
         *pow2(m32 - mU32)) + ((-24*(-3*mU32*(5*mQ32 + mU32)*pow2(mQ32) - pow2(
         m32)*(15*mQ32*mU32 + 29*pow2(mQ32) + 6*pow2(mU32)) + (7*mQ32 + 3*mU32)
         *pow3(m32) + m32*(31*mU32*pow2(mQ32) + 7*mQ32*pow2(mU32) + 16*pow3(
         mQ32)) + 4*pow4(m32)))/((m32 - mU32)*pow2(m32 - mQ32)*pow3(mQ32 - mU32
         )) + 2*((8*(90*m32*mQ32*msq2 + 2*mQ32*pow2(m32) - 3*m32*pow2(mQ32) +
         30*msq2*pow2(mQ32) + pow3(mQ32)))/(3.*pow2(mQ32 - mU32)*pow3(m32 -
         mQ32)) - (4*(mQ32 + mU32)*(-480*m32*mQ32*msq2*mU32 - 600*mQ32*msq2*
         pow2(m32) + 480*mQ32*mU32*pow2(m32) - 600*msq2*mU32*pow2(m32) + 360*
         m32*msq2*pow2(mQ32) - 239*m32*mU32*pow2(mQ32) + 120*msq2*mU32*pow2(
         mQ32) + 120*pow2(m32)*pow2(mQ32) - 239*m32*mQ32*pow2(mU32) + 360*m32*
         msq2*pow2(mU32) + 120*mQ32*msq2*pow2(mU32) + 120*pow2(m32)*pow2(mU32)
         + 119*pow2(mQ32)*pow2(mU32) - 241*mQ32*pow3(m32) + 720*msq2*pow3(m32)
         - 241*mU32*pow3(m32) + 121*pow4(m32)))/(3.*pow2(m32 - mQ32)*pow2(m32 -
         mU32)*pow3(mQ32 - mU32))))*pow4(Xt) + (16*(150*mQ32*msq2*mU32*pow2(
         m32) + 330*m32*msq2*mU32*pow2(mQ32) + 465*msq2*pow2(m32)*pow2(mQ32) -
         439*mU32*pow2(m32)*pow2(mQ32) - 165*m32*mQ32*msq2*pow2(mU32) - 230*
         mQ32*pow2(m32)*pow2(mU32) + 180*msq2*pow2(m32)*pow2(mU32) + 212*m32*
         pow2(mQ32)*pow2(mU32) - 75*msq2*pow2(mQ32)*pow2(mU32) - 705*mQ32*msq2*
         pow3(m32) + 475*mQ32*mU32*pow3(m32) - 300*msq2*mU32*pow3(m32) + 227*
         pow2(mQ32)*pow3(m32) + 87*pow2(mU32)*pow3(m32) - 180*m32*msq2*pow3(
         mQ32) + 143*m32*mU32*pow3(mQ32) - 60*msq2*mU32*pow3(mQ32) - 74*pow2(
         m32)*pow3(mQ32) - 69*pow2(mU32)*pow3(mQ32) - 245*mQ32*pow4(m32) + 360*
         msq2*pow4(m32) - 179*mU32*pow4(m32) + 92*pow5(m32)))/(3.*pow2(m32 -
         mU32)*pow3(m32 - mQ32)) - (640*m3*(-6*m32*mQ32*msq2 - m32*mQ32*mU32 -
         12*m32*msq2*mU32 + 6*mQ32*msq2*mU32 + mQ32*pow2(m32) - m32*pow2(mQ32)
         + 12*msq2*pow2(mQ32) + mU32*pow2(mQ32))*pow5(Xt))/(3.*(m32 - mU32)*
         pow2(m32 - mQ32)*pow3(mQ32 - mU32))) + Xt*((704*m3*mQ32)/(3.*(m32 -
         mQ32)*(mQ32 - mU32)) - (128*m3*(3*m32*mQ32 + 3*m32*mU32 + 5*mQ32*mU32
         - 11*pow2(m32)))/(3.*(m32 - mQ32)*(m32 - mU32)*(mQ32 - mU32)) + (384*
         m3*mQ32*(-(mQ32*mU32) + pow2(m32)))/((m32 - mU32)*(mQ32 - mU32)*pow2(
         m32 - mQ32)) + (256*(-(mU32*dilog(1 - m32/mQ32)) + m32*(dilog(1 -
         m32/mQ32) - dilog(1 - m32/mU32)) + mQ32*dilog(1 - m32/mU32))*(-18*m32*
         mQ32 + pow2(m32) + 9*pow2(mQ32))*pow3(m3))/(3.*(m32 - mU32)*(mQ32 -
         mU32)*pow3(m32 - mQ32)) - 24*((64*(m3*fin(mQ32,m32,MR2) - m3*fin(mQ32,
         mU32,MR2) + m3*fin(mU32,m32,MR2) + 3*pow3(m3) + zt2*pow3(m3)))/(9.*(
         m32 - mQ32)*(m32 - mU32)) - (4*((20*(m3*mQ32 - m3*msq2)*fin(mQ32,msq2,
         MR2))/((mQ32 - mU32)*pow2(m32 - mQ32)) - (2*m3*(mQ32 - mU32)*fin(mU32,
         m32,MR2))/((m32 - mU32)*pow2(m32 - mQ32)) + (2*m3*(mQ32 - mU32)*fin(
         mQ32,m32,MR2))/((m32 - mQ32)*pow2(m32 - mU32)) + (20*(-(m3*msq2) + m3*
         mU32)*fin(mU32,msq2,MR2))/((-mQ32 + mU32)*pow2(m32 - mU32)) + (2*(6*m3
         *mQ32 + 60*m3*msq2 + 6*m3*mU32 + m3*mQ32*zt2 + 10*m3*msq2*zt2 + m3*
         mU32*zt2 + 24*pow3(m3)))/((m32 - mQ32)*(m32 - mU32)) - (20*fin(m32,
         msq2,MR2)*(m3*mQ32*msq2 - 2*m3*mQ32*mU32 + m3*msq2*mU32 + mQ32*pow3(m3
         ) - 2*msq2*pow3(m3) + mU32*pow3(m3)))/(pow2(m32 - mQ32)*pow2(m32 -
         mU32)) + 3*((-2*m3*fin(mQ32,mU32,MR2))/((m32 - mQ32)*(m32 - mU32)) + (
         2*fin(mQ32,m32,MR2)*(m3*mQ32 - 3*m3*mU32 + 2*pow3(m3)))/((m32 - mQ32)*
         (m32 - mU32)*(mQ32 - mU32)) - (2*fin(mU32,m32,MR2)*(-3*m3*mQ32 + m3*
         mU32 + 2*pow3(m3)))/((m32 - mQ32)*(m32 - mU32)*(mQ32 - mU32)) + (2*(11
         *pow3(m3) + 3*zt2*pow3(m3)))/((m32 - mQ32)*(m32 - mU32))) + (2*fin(
         mQ32,mU32,MR2)*(m3*pow2(mQ32) + m3*pow2(mU32) - 2*mQ32*pow3(m3) - 2*
         mU32*pow3(m3) + 2*pow5(m3)))/(pow2(m32 - mQ32)*pow2(m32 - mU32))))/3.)
         + (512*m3*(9*m32*(mQ32 + mU32)*pow2(mQ32)*pow2(mU32) - 3*pow3(mQ32)*
         pow3(mU32) + pow3(m32)*(mU32*(19 + 8*dilog(1 - m32/mU32))*pow2(mQ32) +
         mQ32*(19 + 8*dilog(1 - m32/mQ32))*pow2(mU32) + 2*pow3(mQ32) + 2*pow3(
         mU32)) - 4*pow2(m32)*(6*pow2(mQ32)*pow2(mU32) + mU32*(2 + dilog(1 -
         m32/mU32))*pow3(mQ32) + mQ32*(2 + dilog(1 - m32/mQ32))*pow3(mU32)) - (
         mQ32*mU32*(13 + 4*dilog(1 - m32/mQ32) + 4*dilog(1 - m32/mU32)) + 4*
         pow2(mQ32) + 4*pow2(mU32))*pow4(m32) + 2*(mQ32 + mU32)*pow5(m32)))/(3.
         *(mQ32 - mU32)*mU32*pow2(m32 - mU32)*pow3(m32 - mQ32)) - (64*(-60*m3*
         mQ32*msq2*mU32 + 89*m3*mU32*pow2(mQ32) - 6*m3*mQ32*pow2(mU32) + 60*
         mQ32*msq2*pow3(m3) + 19*mQ32*mU32*pow3(m3) - 115*pow2(mQ32)*pow3(m3) +
         12*m3*pow3(mQ32) + mQ32*pow5(m3) - 16*mU32*pow5(m3) + 16*pow7(m3)))/(
         3.*(m32 - mU32)*(mQ32 - mU32)*pow2(m32 - mQ32))) + pow5(Xt)*((-128*m3*
         mQ32*(3*m32*mQ32 + 3*m32*mU32 + 5*mQ32*mU32 - 11*pow2(m32)))/(3.*(m32
         - mU32)*pow2(m32 - mQ32)*pow3(mQ32 - mU32)) + (1024*m3*((2*mQ32*mU32*(
         -dilog(1 - m32/mQ32) + dilog(1 - m32/mU32)) - pow2(mQ32) + pow2(mU32))
         *pow3(m32) + (-5 + dilog(1 - m32/mQ32) - dilog(1 - m32/mU32))*pow2(
         mU32)*pow3(mQ32) + pow2(m32)*(3*mU32*(-2 + dilog(1 - m32/mQ32) - dilog
         (1 - m32/mU32))*pow2(mQ32) + 3*mQ32*(2 + dilog(1 - m32/mQ32) - dilog(1
         - m32/mU32))*pow2(mU32) + pow3(mQ32) - pow3(mU32)) + (5 + dilog(1 -
         m32/mQ32) - dilog(1 - m32/mU32))*pow2(mQ32)*pow3(mU32) + m32*(4*(
         -dilog(1 - m32/mQ32) + dilog(1 - m32/mU32))*pow2(mQ32)*pow2(mU32) +
         mU32*(5 - dilog(1 - m32/mQ32) + dilog(1 - m32/mU32))*pow3(mQ32) + mQ32
         *(-5 - dilog(1 - m32/mQ32) + dilog(1 - m32/mU32))*pow3(mU32))))/(3.*(
         m32 - mQ32)*(-m32 + mQ32)*(m32 - mU32)*mU32*pow4(mQ32 - mU32)) + 2*((
         12*(mQ32 + mU32)*((64*(m3*fin(mQ32,m32,MR2) - m3*fin(mQ32,mU32,MR2) +
         m3*fin(mU32,m32,MR2) + 3*pow3(m3) + zt2*pow3(m3)))/(9.*(m32 - mQ32)*(
         m32 - mU32)) - (4*((20*(m3*mQ32 - m3*msq2)*fin(mQ32,msq2,MR2))/((mQ32
         - mU32)*pow2(m32 - mQ32)) - (2*m3*(mQ32 - mU32)*fin(mU32,m32,MR2))/((
         m32 - mU32)*pow2(m32 - mQ32)) + (2*m3*(mQ32 - mU32)*fin(mQ32,m32,MR2))
         /((m32 - mQ32)*pow2(m32 - mU32)) + (20*(-(m3*msq2) + m3*mU32)*fin(mU32
         ,msq2,MR2))/((-mQ32 + mU32)*pow2(m32 - mU32)) + (2*(6*m3*mQ32 + 60*m3*
         msq2 + 6*m3*mU32 + m3*mQ32*zt2 + 10*m3*msq2*zt2 + m3*mU32*zt2 + 24*
         pow3(m3)))/((m32 - mQ32)*(m32 - mU32)) - (20*fin(m32,msq2,MR2)*(m3*
         mQ32*msq2 - 2*m3*mQ32*mU32 + m3*msq2*mU32 + mQ32*pow3(m3) - 2*msq2*
         pow3(m3) + mU32*pow3(m3)))/(pow2(m32 - mQ32)*pow2(m32 - mU32)) + 3*((
         -2*m3*fin(mQ32,mU32,MR2))/((m32 - mQ32)*(m32 - mU32)) + (2*fin(mQ32,
         m32,MR2)*(m3*mQ32 - 3*m3*mU32 + 2*pow3(m3)))/((m32 - mQ32)*(m32 - mU32
         )*(mQ32 - mU32)) - (2*fin(mU32,m32,MR2)*(-3*m3*mQ32 + m3*mU32 + 2*pow3
         (m3)))/((m32 - mQ32)*(m32 - mU32)*(mQ32 - mU32)) + (2*(11*pow3(m3) + 3
         *zt2*pow3(m3)))/((m32 - mQ32)*(m32 - mU32))) + (2*fin(mQ32,mU32,MR2)*(
         m3*pow2(mQ32) + m3*pow2(mU32) - 2*mQ32*pow3(m3) - 2*mU32*pow3(m3) + 2*
         pow5(m3)))/(pow2(m32 - mQ32)*pow2(m32 - mU32))))/3.))/pow3(mQ32 - mU32
         ) + (128*(-30*m3*mQ32*msq2*mU32 + 48*m3*mU32*pow2(mQ32) - 3*m3*mQ32*
         pow2(mU32) + 30*mQ32*msq2*pow3(m3) + 11*mQ32*mU32*pow3(m3) - 56*pow2(
         mQ32)*pow3(m3) + 6*m3*pow3(mQ32) - 6*mQ32*pow5(m3) - 8*mU32*pow5(m3) +
         8*pow7(m3)))/(3.*(m32 - mU32)*pow2(m32 - mQ32)*pow3(mQ32 - mU32)))) +
         log(mU32/MR2)*((8*(-360*msq2*mU32*pow2(m32)*pow2(mQ32) - 360*mQ32*
         msq2*pow2(m32)*pow2(mU32) - 360*m32*msq2*pow2(mQ32)*pow2(mU32) + 3816*
         pow2(m32)*pow2(mQ32)*pow2(mU32) + 1080*mQ32*msq2*mU32*pow3(m32) - 60*
         msq2*pow2(mQ32)*pow3(m32) - 2597*mU32*pow2(mQ32)*pow3(m32) - 2585*mQ32
         *pow2(mU32)*pow3(m32) - 60*msq2*pow2(mU32)*pow3(m32) + 180*m32*msq2*
         mU32*pow3(mQ32) + 654*mU32*pow2(m32)*pow3(mQ32) - 1229*m32*pow2(mU32)*
         pow3(mQ32) + 60*msq2*pow2(mU32)*pow3(mQ32) - 115*pow3(m32)*pow3(mQ32)
         + 180*m32*mQ32*msq2*pow3(mU32) + 650*mQ32*pow2(m32)*pow3(mU32) - 1229*
         m32*pow2(mQ32)*pow3(mU32) + 60*msq2*pow2(mQ32)*pow3(mU32) - 111*pow3(
         m32)*pow3(mU32) + 382*pow3(mQ32)*pow3(mU32) - 180*mQ32*msq2*pow4(m32)
         + 1626*mQ32*mU32*pow4(m32) - 180*msq2*mU32*pow4(m32) + 697*pow2(mQ32)*
         pow4(m32) + 685*pow2(mU32)*pow4(m32) + 51*m32*mU32*pow4(mQ32) - 36*
         pow2(m32)*pow4(mQ32) + 21*pow2(mU32)*pow4(mQ32) + 51*m32*mQ32*pow4(
         mU32) - 36*pow2(m32)*pow4(mU32) + 21*pow2(mQ32)*pow4(mU32) - 378*mQ32*
         pow5(m32) - 370*mU32*pow5(m32) + 32*pow6(m32)))/(3.*pow3(m32 - mQ32)*
         pow3(m32 - mU32)) + (32*pow2(Xt)*(720*msq2*pow2(m32)*pow2(mQ32)*pow2(
         mU32) - 30*msq2*mU32*pow2(mQ32)*pow3(m32) - 30*mQ32*msq2*pow2(mU32)*
         pow3(m32) - 1817*pow2(mQ32)*pow2(mU32)*pow3(m32) - 360*msq2*mU32*pow2(
         m32)*pow3(mQ32) - 90*m32*msq2*pow2(mU32)*pow3(mQ32) + 381*pow2(m32)*
         pow2(mU32)*pow3(mQ32) + 30*msq2*pow3(m32)*pow3(mQ32) + 589*mU32*pow3(
         m32)*pow3(mQ32) - 360*mQ32*msq2*pow2(m32)*pow3(mU32) - 90*m32*msq2*
         pow2(mQ32)*pow3(mU32) + 349*pow2(m32)*pow2(mQ32)*pow3(mU32) + 629*mQ32
         *pow3(m32)*pow3(mU32) + 30*msq2*pow3(m32)*pow3(mU32) - 107*m32*pow3(
         mQ32)*pow3(mU32) - 60*msq2*pow3(mQ32)*pow3(mU32) - 180*mQ32*msq2*mU32*
         pow4(m32) + 90*msq2*pow2(mQ32)*pow4(m32) + 668*mU32*pow2(mQ32)*pow4(
         m32) + 640*mQ32*pow2(mU32)*pow4(m32) + 90*msq2*pow2(mU32)*pow4(m32) -
         401*pow3(mQ32)*pow4(m32) - 415*pow3(mU32)*pow4(m32) + 90*m32*msq2*mU32
         *pow4(mQ32) - 289*mU32*pow2(m32)*pow4(mQ32) + 12*m32*pow2(mU32)*pow4(
         mQ32) + 30*msq2*pow2(mU32)*pow4(mQ32) + 113*pow3(m32)*pow4(mQ32) + 16*
         pow3(mU32)*pow4(mQ32) + 90*m32*mQ32*msq2*pow4(mU32) - 305*mQ32*pow2(
         m32)*pow4(mU32) + 29*m32*pow2(mQ32)*pow4(mU32) + 30*msq2*pow2(mQ32)*
         pow4(mU32) + 118*pow3(m32)*pow4(mU32) + 10*pow3(mQ32)*pow4(mU32) - 593
         *mQ32*mU32*pow5(m32) + 126*pow2(mQ32)*pow5(m32) + 139*pow2(mU32)*pow5(
         m32) - 3*m32*mU32*pow5(mQ32) + 18*pow2(m32)*pow5(mQ32) - 3*pow2(mU32)*
         pow5(mQ32) - 3*m32*mQ32*pow5(mU32) + 18*pow2(m32)*pow5(mU32) - 3*pow2(
         mQ32)*pow5(mU32) + 44*mQ32*pow6(m32) + 40*mU32*pow6(m32)))/(3.*pow2(
         mQ32 - mU32)*pow3(m32 - mQ32)*pow3(m32 - mU32)) + (5120*m32*mQ32*mU32*
         pow6(Xt))/(3.*(m32 - mQ32)*(m32 - mU32)*pow4(mQ32 - mU32)) + log(
         msq2/MR2)*((24*(4*m32*mQ32*mU32*(mQ32 + mU32) - 2*pow2(mQ32)*pow2(mU32
         ) - pow2(m32)*(8*mQ32*mU32 + pow2(mQ32) + pow2(mU32)) + 2*(mQ32 + mU32
         )*pow3(m32)))/(pow2(m32 - mQ32)*pow2(m32 - mU32)) + 2*((4*(90*m32*mQ32
         *msq2 + 2*mQ32*pow2(m32) - 3*m32*pow2(mQ32) + 30*msq2*pow2(mQ32) +
         pow3(mQ32)))/(3.*pow3(m32 - mQ32)) + (4*(90*m32*msq2*mU32 + 2*mU32*
         pow2(m32) - 3*m32*pow2(mU32) + 30*msq2*pow2(mU32) + pow3(mU32)))/(3.*
         pow3(m32 - mU32))) + (2*((32*m3*mQ32*(m32 - mQ32 + 60*msq2))/(3.*pow2(
         m32 - mQ32)*pow2(mQ32 - mU32)) + (32*m3*(m32 + 60*msq2 - mU32)*mU32)/(
         3.*pow2(m32 - mU32)*pow2(-mQ32 + mU32))) + (192*(-2*m3*mQ32*mU32 + (
         mQ32 + mU32)*pow3(m3)))/((m32 - mQ32)*(m32 - mU32)*pow2(mQ32 - mU32)))
         *pow3(Xt) + pow2(Xt)*(2*((-8*(90*m32*mQ32*msq2 + 2*mQ32*pow2(m32) - 3*
         m32*pow2(mQ32) + 30*msq2*pow2(mQ32) + pow3(mQ32)))/(3.*(mQ32 - mU32)*
         pow3(m32 - mQ32)) + (8*(90*m32*msq2*mU32 + 2*mU32*pow2(m32) - 3*m32*
         pow2(mU32) + 30*msq2*pow2(mU32) + pow3(mU32)))/(3.*(mQ32 - mU32)*pow3(
         m32 - mU32))) + (48*(2*(mQ32 + mU32)*pow2(mQ32)*pow2(mU32) - 4*m32*
         mQ32*mU32*pow2(mQ32 + mU32) - 2*(2*mQ32*mU32 + 3*pow2(mQ32) + 3*pow2(
         mU32))*pow3(m32) + 3*pow2(m32)*pow3(mQ32 + mU32) + 2*(mQ32 + mU32)*
         pow4(m32)))/(pow2(m32 - mQ32)*pow2(m32 - mU32)*pow2(mQ32 - mU32))) + (
         2*((4*(mQ32 + mU32)*(90*m32*mQ32*msq2 + 2*mQ32*pow2(m32) - 3*m32*pow2(
         mQ32) + 30*msq2*pow2(mQ32) + pow3(mQ32)))/(3.*pow3(m32 - mQ32)*pow3(
         mQ32 - mU32)) - (4*(mQ32 + mU32)*(90*m32*msq2*mU32 + 2*mU32*pow2(m32)
         - 3*m32*pow2(mU32) + 30*msq2*pow2(mU32) + pow3(mU32)))/(3.*pow3(m32 -
         mU32)*pow3(mQ32 - mU32))) - (24*(mQ32 + mU32)*(4*(mQ32 + mU32)*pow2(
         mQ32)*pow2(mU32) - 8*m32*mQ32*mU32*pow2(mQ32 + mU32) - 2*(6*mQ32*mU32
         + 5*pow2(mQ32) + 5*pow2(mU32))*pow3(m32) + pow2(m32)*(19*mU32*pow2(
         mQ32) + 19*mQ32*pow2(mU32) + 5*pow3(mQ32) + 5*pow3(mU32)) + 4*(mQ32 +
         mU32)*pow4(m32)))/(pow2(m32 - mQ32)*pow2(m32 - mU32)*pow4(mQ32 - mU32)
         ))*pow4(Xt) + ((-96*m3*(mQ32 + mU32)*(-2*mQ32*mU32 + m32*(mQ32 + mU32)
         ))/((m32 - mQ32)*(m32 - mU32)*pow4(mQ32 - mU32)) + 2*((-16*m3*mQ32*(
         m32 - mQ32 + 60*msq2)*(mQ32 + mU32))/(3.*pow2(m32 - mQ32)*pow4(mQ32 -
         mU32)) - (16*m3*(m32 + 60*msq2 - mU32)*mU32*(mQ32 + mU32))/(3.*pow2(
         m32 - mU32)*pow4(-mQ32 + mU32))))*pow5(Xt) - (320*Xt*(-6*m3*mQ32*msq2*
         mU32 + mQ32*mU32*pow3(m3) - mQ32*pow5(m3) + 6*msq2*pow5(m3) - mU32*
         pow5(m3) + pow7(m3)))/(3.*pow2(m32 - mQ32)*pow2(m32 - mU32))) + (16*
         pow4(Xt)*(60*msq2*pow2(mQ32)*pow2(mU32)*pow3(m32) - 360*msq2*pow2(m32)
         *pow2(mU32)*pow3(mQ32) + 7752*pow2(mU32)*pow3(m32)*pow3(mQ32) - 360*
         msq2*pow2(m32)*pow2(mQ32)*pow3(mU32) + 7728*pow2(mQ32)*pow3(m32)*pow3(
         mU32) + 180*m32*msq2*pow3(mQ32)*pow3(mU32) - 4272*pow2(m32)*pow3(mQ32)
         *pow3(mU32) + 90*msq2*mU32*pow2(mQ32)*pow4(m32) + 90*mQ32*msq2*pow2(
         mU32)*pow4(m32) - 9160*pow2(mQ32)*pow2(mU32)*pow4(m32) - 90*msq2*pow3(
         mQ32)*pow4(m32) - 3974*mU32*pow3(mQ32)*pow4(m32) - 3942*mQ32*pow3(mU32
         )*pow4(m32) - 90*msq2*pow3(mU32)*pow4(m32) + 360*msq2*mU32*pow2(m32)*
         pow4(mQ32) - 2394*pow2(m32)*pow2(mU32)*pow4(mQ32) - 30*msq2*pow3(m32)*
         pow4(mQ32) + 859*mU32*pow3(m32)*pow4(mQ32) + 703*m32*pow3(mU32)*pow4(
         mQ32) + 30*msq2*pow3(mU32)*pow4(mQ32) + 256*pow4(m32)*pow4(mQ32) + 360
         *mQ32*msq2*pow2(m32)*pow4(mU32) - 2382*pow2(m32)*pow2(mQ32)*pow4(mU32)
         + 847*mQ32*pow3(m32)*pow4(mU32) - 30*msq2*pow3(m32)*pow4(mU32) + 699*
         m32*pow3(mQ32)*pow4(mU32) + 30*msq2*pow3(mQ32)*pow4(mU32) + 260*pow4(
         m32)*pow4(mU32) + 132*pow4(mQ32)*pow4(mU32) + 3901*mU32*pow2(mQ32)*
         pow5(m32) + 3877*mQ32*pow2(mU32)*pow5(m32) + 197*pow3(mQ32)*pow5(m32)
         + 185*pow3(mU32)*pow5(m32) - 90*m32*msq2*mU32*pow5(mQ32) + 102*mU32*
         pow2(m32)*pow5(mQ32) + 313*m32*pow2(mU32)*pow5(mQ32) - 30*msq2*pow2(
         mU32)*pow5(mQ32) - 113*pow3(m32)*pow5(mQ32) - 106*pow3(mU32)*pow5(mQ32
         ) - 90*m32*mQ32*msq2*pow5(mU32) + 102*mQ32*pow2(m32)*pow5(mU32) + 313*
         m32*pow2(mQ32)*pow5(mU32) - 30*msq2*pow2(mQ32)*pow5(mU32) - 113*pow3(
         m32)*pow5(mU32) - 106*pow3(mQ32)*pow5(mU32) - 1180*mQ32*mU32*pow6(m32)
         - 224*pow2(mQ32)*pow6(m32) - 212*pow2(mU32)*pow6(m32) - 6*m32*mU32*
         pow6(mQ32) - 18*pow2(m32)*pow6(mQ32) - 6*m32*mQ32*pow6(mU32) - 18*pow2
         (m32)*pow6(mU32) + 2*mQ32*pow7(m32) - 2*mU32*pow7(m32)))/(3.*pow3(m32
         - mQ32)*pow3(m32 - mU32)*pow4(-mQ32 + mU32)) - (64*pow3(Xt)*(-480*msq2
         *mU32*pow2(mQ32)*pow3(m3) + 480*mQ32*msq2*pow2(mU32)*pow3(m3) - 8*pow2
         (mQ32)*pow2(mU32)*pow3(m3) + 120*m3*msq2*mU32*pow3(mQ32) - 481*m3*pow2
         (mU32)*pow3(mQ32) + 612*mU32*pow3(m3)*pow3(mQ32) - 120*m3*mQ32*msq2*
         pow3(mU32) + 485*m3*pow2(mQ32)*pow3(mU32) - 620*mQ32*pow3(m3)*pow3(
         mU32) - 12*m3*mU32*pow4(mQ32) + 24*pow3(m3)*pow4(mQ32) + 12*m3*mQ32*
         pow4(mU32) - 24*pow3(m3)*pow4(mU32) + 120*msq2*pow2(mQ32)*pow5(m3) -
         57*mU32*pow2(mQ32)*pow5(m3) + 77*mQ32*pow2(mU32)*pow5(m3) - 120*msq2*
         pow2(mU32)*pow5(m3) - 279*pow3(mQ32)*pow5(m3) + 283*pow3(mU32)*pow5(m3
         ) - 8*mQ32*mU32*pow7(m3) + 16*pow2(mQ32)*pow7(m3) - 24*pow2(mU32)*pow7
         (m3) + 63*mQ32*pow9(m3) - 59*mU32*pow9(m3)))/(3.*pow2(m32 - mQ32)*pow2
         (m32 - mU32)*pow3(mQ32 - mU32)) + (32*Xt*(240*m3*msq2*pow2(mQ32)*pow2(
         mU32) - 1011*pow2(mQ32)*pow2(mU32)*pow3(m3) - 120*m3*msq2*mU32*pow3(
         mQ32) + 59*m3*pow2(mU32)*pow3(mQ32) + 431*mU32*pow3(m3)*pow3(mQ32) -
         120*m3*mQ32*msq2*pow3(mU32) + 57*m3*pow2(mQ32)*pow3(mU32) + 436*mQ32*
         pow3(m3)*pow3(mU32) - 42*m3*mU32*pow4(mQ32) + 24*pow3(m3)*pow4(mQ32) -
         42*m3*mQ32*pow4(mU32) + 24*pow3(m3)*pow4(mU32) - 240*mQ32*msq2*mU32*
         pow5(m3) + 120*msq2*pow2(mQ32)*pow5(m3) + 378*mU32*pow2(mQ32)*pow5(m3)
         + 369*mQ32*pow2(mU32)*pow5(m3) + 120*msq2*pow2(mU32)*pow5(m3) - 324*
         pow3(mQ32)*pow5(m3) - 327*pow3(mU32)*pow5(m3) - 251*mQ32*mU32*pow7(m3)
         + 106*pow2(mQ32)*pow7(m3) + 113*pow2(mU32)*pow7(m3) + 2*mQ32*pow9(m3)
         - 2*mU32*pow9(m3)))/(3.*pow2(m32 - mQ32)*pow2(m32 - mU32)*pow2(mQ32 -
         mU32)) + (128*pow5(Xt)*(60*m3*msq2*pow2(mQ32)*pow2(mU32) - 120*msq2*
         mU32*pow2(mQ32)*pow3(m3) - 120*mQ32*msq2*pow2(mU32)*pow3(m3) + 474*
         pow2(mQ32)*pow2(mU32)*pow3(m3) + 30*m3*msq2*mU32*pow3(mQ32) - 177*m3*
         pow2(mU32)*pow3(mQ32) + 190*mU32*pow3(m3)*pow3(mQ32) + 30*m3*mQ32*msq2
         *pow3(mU32) - 177*m3*pow2(mQ32)*pow3(mU32) + 190*mQ32*pow3(m3)*pow3(
         mU32) + 6*pow3(m3)*pow4(mQ32) + 6*pow3(m3)*pow4(mU32) + 60*mQ32*msq2*
         mU32*pow5(m3) + 30*msq2*pow2(mQ32)*pow5(m3) - 230*mU32*pow2(mQ32)*pow5
         (m3) - 230*mQ32*pow2(mU32)*pow5(m3) + 30*msq2*pow2(mU32)*pow5(m3) - 73
         *pow3(mQ32)*pow5(m3) - 73*pow3(mU32)*pow5(m3) + 40*mQ32*mU32*pow7(m3)
         - 5*pow2(mQ32)*pow7(m3) - 5*pow2(mU32)*pow7(m3) + 32*mQ32*pow9(m3) +
         32*mU32*pow9(m3)))/(3.*pow2(m32 - mQ32)*pow2(m32 - mU32)*pow4(mQ32 -
         mU32))) + pow2(log(mU32/MR2))*((1280*m32*(mQ32 + mU32)*(2*m32*mQ32*
         mU32 + m32*pow2(mU32) - 3*mQ32*pow2(mU32))*pow6(Xt))/(3.*(m32 - mQ32)*
         pow2(m32 - mU32)*pow5(mQ32 - mU32)) - (4*(90*msq2*pow2(mQ32)*pow2(mU32
         )*pow3(m32) + 90*msq2*pow2(m32)*pow2(mU32)*pow3(mQ32) - 270*msq2*mU32*
         pow3(m32)*pow3(mQ32) + 1452*pow2(mU32)*pow3(m32)*pow3(mQ32) - 270*msq2
         *pow2(m32)*pow2(mQ32)*pow3(mU32) + 210*mQ32*msq2*pow3(m32)*pow3(mU32)
         - 2028*pow2(mQ32)*pow3(m32)*pow3(mU32) + 150*m32*msq2*pow3(mQ32)*pow3(
         mU32) - 712*pow2(m32)*pow3(mQ32)*pow3(mU32) + 270*msq2*mU32*pow2(mQ32)
         *pow4(m32) - 210*mQ32*msq2*pow2(mU32)*pow4(m32) + 639*pow2(mQ32)*pow2(
         mU32)*pow4(m32) - 849*mU32*pow3(mQ32)*pow4(m32) + 2517*mQ32*pow3(mU32)
         *pow4(m32) - 60*msq2*pow3(mU32)*pow4(m32) + 90*msq2*mU32*pow2(m32)*
         pow4(mQ32) - 60*m32*msq2*pow2(mU32)*pow4(mQ32) - 701*pow2(m32)*pow2(
         mU32)*pow4(mQ32) + 306*mU32*pow3(m32)*pow4(mQ32) + 616*m32*pow3(mU32)*
         pow4(mQ32) - 30*msq2*pow3(mU32)*pow4(mQ32) - 104*pow4(m32)*pow4(mQ32)
         + 90*mQ32*msq2*pow2(m32)*pow4(mU32) - 90*m32*msq2*pow2(mQ32)*pow4(mU32
         ) + 1930*pow2(m32)*pow2(mQ32)*pow4(mU32) - 1626*mQ32*pow3(m32)*pow4(
         mU32) - 30*msq2*pow3(m32)*pow4(mU32) - 296*m32*pow3(mQ32)*pow4(mU32) +
         30*msq2*pow3(mQ32)*pow4(mU32) + 357*pow4(m32)*pow4(mU32) - 125*pow4(
         mQ32)*pow4(mU32) - 90*mQ32*msq2*mU32*pow5(m32) + 70*mU32*pow2(mQ32)*
         pow5(m32) - 1678*mQ32*pow2(mU32)*pow5(m32) + 90*msq2*pow2(mU32)*pow5(
         m32) + 312*pow3(mQ32)*pow5(m32) - 624*pow3(mU32)*pow5(m32) + 27*mU32*
         pow2(m32)*pow5(mQ32) - 18*m32*pow2(mU32)*pow5(mQ32) - 9*pow3(mU32)*
         pow5(mQ32) + 233*mQ32*pow2(m32)*pow5(mU32) - 454*m32*pow2(mQ32)*pow5(
         mU32) - 24*pow3(m32)*pow5(mU32) + 125*pow3(mQ32)*pow5(mU32) + 556*mQ32
         *mU32*pow6(m32) - 214*pow2(mQ32)*pow6(m32) + 426*pow2(mU32)*pow6(m32)
         + 24*m32*mQ32*pow6(mU32) - 9*pow2(m32)*pow6(mU32) + 9*pow2(mQ32)*pow6(
         mU32) + 6*mQ32*pow7(m32) - 134*mU32*pow7(m32)))/(3.*(mQ32 - mU32)*pow3
         (m32 - mQ32)*pow4(m32 - mU32)) - (8*pow2(Xt)*(90*msq2*pow2(mQ32)*pow2(
         mU32)*pow3(m32) + 90*msq2*pow2(m32)*pow2(mU32)*pow3(mQ32) - 270*msq2*
         mU32*pow3(m32)*pow3(mQ32) + 716*pow2(mU32)*pow3(m32)*pow3(mQ32) - 270*
         msq2*pow2(m32)*pow2(mQ32)*pow3(mU32) + 210*mQ32*msq2*pow3(m32)*pow3(
         mU32) + 1036*pow2(mQ32)*pow3(m32)*pow3(mU32) + 150*m32*msq2*pow3(mQ32)
         *pow3(mU32) - 48*pow2(m32)*pow3(mQ32)*pow3(mU32) + 270*msq2*mU32*pow2(
         mQ32)*pow4(m32) - 210*mQ32*msq2*pow2(mU32)*pow4(m32) - 3269*pow2(mQ32)
         *pow2(mU32)*pow4(m32) - 801*mU32*pow3(mQ32)*pow4(m32) + 753*mQ32*pow3(
         mU32)*pow4(m32) - 60*msq2*pow3(mU32)*pow4(m32) + 90*msq2*mU32*pow2(m32
         )*pow4(mQ32) - 60*m32*msq2*pow2(mU32)*pow4(mQ32) + 523*pow2(m32)*pow2(
         mU32)*pow4(mQ32) - 330*mU32*pow3(m32)*pow4(mQ32) - 340*m32*pow3(mU32)*
         pow4(mQ32) - 30*msq2*pow3(mU32)*pow4(mQ32) + 48*pow4(m32)*pow4(mQ32) +
         90*mQ32*msq2*pow2(m32)*pow4(mU32) - 90*m32*msq2*pow2(mQ32)*pow4(mU32)
         + 166*pow2(m32)*pow2(mQ32)*pow4(mU32) - 1030*mQ32*pow3(m32)*pow4(mU32
         ) - 30*msq2*pow3(m32)*pow4(mU32) - 176*m32*pow3(mQ32)*pow4(mU32) + 30*
         msq2*pow3(mQ32)*pow4(mU32) + 709*pow4(m32)*pow4(mU32) + 91*pow4(mQ32)*
         pow4(mU32) - 90*mQ32*msq2*mU32*pow5(m32) + 2830*mU32*pow2(mQ32)*pow5(
         m32) + 1262*mQ32*pow2(mU32)*pow5(m32) + 90*msq2*pow2(mU32)*pow5(m32) +
         368*pow3(mQ32)*pow5(m32) - 1164*pow3(mU32)*pow5(m32) - 9*mU32*pow2(
         m32)*pow5(mQ32) + 6*m32*pow2(mU32)*pow5(mQ32) + 3*pow3(mU32)*pow5(mQ32
         ) + 113*mQ32*pow2(m32)*pow5(mU32) + 106*m32*pow2(mQ32)*pow5(mU32) - 72
         *pow3(m32)*pow5(mU32) - 27*pow3(mQ32)*pow5(mU32) - 1552*mQ32*mU32*pow6
         (m32) - 914*pow2(mQ32)*pow6(m32) + 642*pow2(mU32)*pow6(m32) - 12*m32*
         mQ32*pow6(mU32) - 9*pow2(m32)*pow6(mU32) - 3*pow2(mQ32)*pow6(mU32) +
         498*mQ32*pow7(m32) - 114*mU32*pow7(m32)))/(3.*pow2(mQ32 - mU32)*pow3(
         m32 - mQ32)*pow4(m32 - mU32)) - (4*pow4(Xt)*(180*msq2*pow2(mU32)*pow3(
         m32)*pow3(mQ32) - 300*msq2*pow2(mQ32)*pow3(m32)*pow3(mU32) + 180*msq2*
         pow2(m32)*pow3(mQ32)*pow3(mU32) - 14888*pow3(m32)*pow3(mQ32)*pow3(mU32
         ) - 60*msq2*pow2(mQ32)*pow2(mU32)*pow4(m32) - 270*msq2*mU32*pow3(mQ32)
         *pow4(m32) + 20822*pow2(mU32)*pow3(mQ32)*pow4(m32) + 270*mQ32*msq2*
         pow3(mU32)*pow4(m32) + 9860*pow2(mQ32)*pow3(mU32)*pow4(m32) - 180*msq2
         *pow2(m32)*pow2(mU32)*pow4(mQ32) + 270*msq2*mU32*pow3(m32)*pow4(mQ32)
         - 5650*pow2(mU32)*pow3(m32)*pow4(mQ32) - 90*m32*msq2*pow3(mU32)*pow4(
         mQ32) + 5189*pow2(m32)*pow3(mU32)*pow4(mQ32) + 2577*mU32*pow4(m32)*
         pow4(mQ32) + 180*msq2*pow2(m32)*pow2(mQ32)*pow4(mU32) - 180*mQ32*msq2*
         pow3(m32)*pow4(mU32) + 5666*pow2(mQ32)*pow3(m32)*pow4(mU32) - 60*m32*
         msq2*pow3(mQ32)*pow4(mU32) + 930*pow2(m32)*pow3(mQ32)*pow4(mU32) -
         7470*mQ32*pow4(m32)*pow4(mU32) + 60*msq2*pow4(m32)*pow4(mU32) - 1476*
         m32*pow4(mQ32)*pow4(mU32) + 90*msq2*mU32*pow2(mQ32)*pow5(m32) - 20588*
         pow2(mQ32)*pow2(mU32)*pow5(m32) - 11094*mU32*pow3(mQ32)*pow5(m32) +
         558*mQ32*pow3(mU32)*pow5(m32) - 90*msq2*pow3(mU32)*pow5(m32) - 536*
         pow4(mQ32)*pow5(m32) + 2476*pow4(mU32)*pow5(m32) - 90*msq2*mU32*pow2(
         m32)*pow5(mQ32) + 60*m32*msq2*pow2(mU32)*pow5(mQ32) - 642*pow2(m32)*
         pow2(mU32)*pow5(mQ32) + 330*mU32*pow3(m32)*pow5(mQ32) + 462*m32*pow3(
         mU32)*pow5(mQ32) + 30*msq2*pow3(mU32)*pow5(mQ32) - 16*pow4(m32)*pow5(
         mQ32) - 126*pow4(mU32)*pow5(mQ32) - 90*mQ32*msq2*pow2(m32)*pow5(mU32)
         + 90*m32*msq2*pow2(mQ32)*pow5(mU32) - 4751*pow2(m32)*pow2(mQ32)*pow5(
         mU32) + 4422*mQ32*pow3(m32)*pow5(mU32) + 30*msq2*pow3(m32)*pow5(mU32)
         + 2030*m32*pow3(mQ32)*pow5(mU32) - 30*msq2*pow3(mQ32)*pow5(mU32) -
         1453*pow4(m32)*pow5(mU32) - 128*pow4(mQ32)*pow5(mU32) + 12522*mU32*
         pow2(mQ32)*pow6(m32) + 6014*mQ32*pow2(mU32)*pow6(m32) + 2218*pow3(mQ32
         )*pow6(m32) - 1298*pow3(mU32)*pow6(m32) + 9*mU32*pow2(m32)*pow6(mQ32)
         - 6*m32*pow2(mU32)*pow6(mQ32) - 3*pow3(mU32)*pow6(mQ32) - 232*mQ32*
         pow2(m32)*pow6(mU32) + 2*m32*pow2(mQ32)*pow6(mU32) + 136*pow3(m32)*
         pow6(mU32) - 2*pow3(mQ32)*pow6(mU32) - 4360*mQ32*mU32*pow7(m32) - 2698
         *pow2(mQ32)*pow7(m32) + 146*pow2(mU32)*pow7(m32) + 12*m32*mQ32*pow7(
         mU32) + 9*pow2(m32)*pow7(mU32) + 3*pow2(mQ32)*pow7(mU32) + 1032*mQ32*
         pow8(m32) - 8*mU32*pow8(m32)))/(3.*pow3(m32 - mQ32)*pow4(m32 - mU32)*
         pow4(mQ32 - mU32)) - (32*Xt*(4*mQ32*pow11(m3) - 4*mU32*pow11(m3) - 30*
         msq2*pow2(mQ32)*pow2(mU32)*pow3(m3) + 30*m3*msq2*pow2(mU32)*pow3(mQ32)
         - 30*msq2*mU32*pow3(m3)*pow3(mQ32) + 57*pow2(mU32)*pow3(m3)*pow3(mQ32
         ) - 30*m3*msq2*pow2(mQ32)*pow3(mU32) + 60*mQ32*msq2*pow3(m3)*pow3(mU32
         ) + 162*pow2(mQ32)*pow3(m3)*pow3(mU32) - 60*m3*pow3(mQ32)*pow3(mU32) +
         9*m3*pow2(mU32)*pow4(mQ32) - 9*mU32*pow3(m3)*pow4(mQ32) + 26*m3*pow2(
         mQ32)*pow4(mU32) - 143*mQ32*pow3(m3)*pow4(mU32) + 60*msq2*mU32*pow2(
         mQ32)*pow5(m3) - 30*mQ32*msq2*pow2(mU32)*pow5(m3) - 329*pow2(mQ32)*
         pow2(mU32)*pow5(m3) + 32*mU32*pow3(mQ32)*pow5(m3) + 124*mQ32*pow3(mU32
         )*pow5(m3) - 30*msq2*pow3(mU32)*pow5(m3) + 77*pow4(mU32)*pow5(m3) + 9*
         m3*mQ32*pow5(mU32) - 3*pow3(m3)*pow5(mU32) - 30*mQ32*msq2*mU32*pow7(m3
         ) + 116*mU32*pow2(mQ32)*pow7(m3) + 71*mQ32*pow2(mU32)*pow7(m3) + 30*
         msq2*pow2(mU32)*pow7(m3) - 13*pow3(mQ32)*pow7(m3) - 110*pow3(mU32)*
         pow7(m3) - 17*mQ32*mU32*pow9(m3) - 23*pow2(mQ32)*pow9(m3) + 24*pow2(
         mU32)*pow9(m3)))/(3.*pow2(m32 - mQ32)*pow2(mQ32 - mU32)*pow3(m32 -
         mU32)) - (64*pow3(Xt)*(-37*mQ32*pow11(m3) - 7*mU32*pow11(m3) + 2*pow13
         (m3) - 30*msq2*pow2(mQ32)*pow2(mU32)*pow3(m3) + 30*m3*msq2*pow2(mU32)*
         pow3(mQ32) - 30*msq2*mU32*pow3(m3)*pow3(mQ32) + 8*pow2(mU32)*pow3(m3)*
         pow3(mQ32) - 30*m3*msq2*pow2(mQ32)*pow3(mU32) + 60*mQ32*msq2*pow3(m3)*
         pow3(mU32) - 441*pow2(mQ32)*pow3(m3)*pow3(mU32) - 7*m3*pow3(mQ32)*pow3
         (mU32) - 3*m3*pow2(mU32)*pow4(mQ32) + 3*mU32*pow3(m3)*pow4(mQ32) + 175
         *m3*pow2(mQ32)*pow4(mU32) - 251*mQ32*pow3(m3)*pow4(mU32) + 60*msq2*
         mU32*pow2(mQ32)*pow5(m3) - 30*mQ32*msq2*pow2(mU32)*pow5(m3) + 436*pow2
         (mQ32)*pow2(mU32)*pow5(m3) + 7*mU32*pow3(mQ32)*pow5(m3) + 567*mQ32*
         pow3(mU32)*pow5(m3) - 30*msq2*pow3(mU32)*pow5(m3) + 108*pow4(mU32)*
         pow5(m3) - 3*m3*mQ32*pow5(mU32) - 3*pow3(m3)*pow5(mU32) - 30*mQ32*msq2
         *mU32*pow7(m3) - 237*mU32*pow2(mQ32)*pow7(m3) - 472*mQ32*pow2(mU32)*
         pow7(m3) + 30*msq2*pow2(mU32)*pow7(m3) + 8*pow3(mQ32)*pow7(m3) - 171*
         pow3(mU32)*pow7(m3) + 244*mQ32*mU32*pow9(m3) + 19*pow2(mQ32)*pow9(m3)
         + 55*pow2(mU32)*pow9(m3)))/(3.*pow2(m32 - mQ32)*pow3(m32 - mU32)*pow3(
         mQ32 - mU32)) - (32*pow5(Xt)*(-32*mQ32*mU32*pow11(m3) + 32*pow11(m3)*
         pow2(mQ32) + 60*msq2*pow2(mU32)*pow3(m3)*pow3(mQ32) - 30*msq2*pow2(
         mQ32)*pow3(m3)*pow3(mU32) + 679*pow3(m3)*pow3(mQ32)*pow3(mU32) - 30*m3
         *msq2*pow2(mU32)*pow4(mQ32) + 30*msq2*mU32*pow3(m3)*pow4(mQ32) + 74*
         pow2(mU32)*pow3(m3)*pow4(mQ32) - 29*m3*pow3(mU32)*pow4(mQ32) + 30*m3*
         msq2*pow2(mQ32)*pow4(mU32) - 60*mQ32*msq2*pow3(m3)*pow4(mU32) + 851*
         pow2(mQ32)*pow3(m3)*pow4(mU32) - 246*m3*pow3(mQ32)*pow4(mU32) - 30*
         msq2*pow2(mQ32)*pow2(mU32)*pow5(m3) - 60*msq2*mU32*pow3(mQ32)*pow5(m3)
         - 713*pow2(mU32)*pow3(mQ32)*pow5(m3) + 60*mQ32*msq2*pow3(mU32)*pow5(
         m3) - 1177*pow2(mQ32)*pow3(mU32)*pow5(m3) - 60*mU32*pow4(mQ32)*pow5(m3
         ) - 791*mQ32*pow4(mU32)*pow5(m3) + 30*msq2*pow4(mU32)*pow5(m3) + 3*m3*
         pow2(mU32)*pow5(mQ32) - 3*mU32*pow3(m3)*pow5(mQ32) - 211*m3*pow2(mQ32)
         *pow5(mU32) + 316*mQ32*pow3(m3)*pow5(mU32) - 139*pow5(m3)*pow5(mU32) +
         3*m3*mQ32*pow6(mU32) + 3*pow3(m3)*pow6(mU32) + 30*msq2*mU32*pow2(mQ32
         )*pow7(m3) + 731*pow2(mQ32)*pow2(mU32)*pow7(m3) + 343*mU32*pow3(mQ32)*
         pow7(m3) + 617*mQ32*pow3(mU32)*pow7(m3) - 30*msq2*pow3(mU32)*pow7(m3)
         - pow4(mQ32)*pow7(m3) + 230*pow4(mU32)*pow7(m3) - 226*mU32*pow2(mQ32)*
         pow9(m3) - 145*mQ32*pow2(mU32)*pow9(m3) - 31*pow3(mQ32)*pow9(m3) - 78*
         pow3(mU32)*pow9(m3)))/(3.*pow2(m32 - mQ32)*pow3(m32 - mU32)*pow5(mQ32
         - mU32)))) + pow2(log(mQ32/MR2))*((-16*pow2(Xt)*(-420*msq2*mU32*pow2(
         m32)*pow2(mQ32) + 360*mQ32*msq2*pow2(m32)*pow2(mU32) + 300*m32*msq2*
         pow2(mQ32)*pow2(mU32) - 1429*pow2(m32)*pow2(mQ32)*pow2(mU32) - 180*
         mQ32*msq2*mU32*pow3(m32) + 180*msq2*pow2(mQ32)*pow3(m32) + 255*mU32*
         pow2(mQ32)*pow3(m32) + 770*mQ32*pow2(mU32)*pow3(m32) - 120*m32*msq2*
         mU32*pow3(mQ32) + 60*msq2*pow2(m32)*pow3(mQ32) + 1392*mU32*pow2(m32)*
         pow3(mQ32) - 168*m32*pow2(mU32)*pow3(mQ32) + 60*msq2*pow2(mU32)*pow3(
         mQ32) - 732*pow3(m32)*pow3(mQ32) - 180*m32*mQ32*msq2*pow3(mU32) - 230*
         mQ32*pow2(m32)*pow3(mU32) + 719*m32*pow2(mQ32)*pow3(mU32) - 60*msq2*
         pow2(mQ32)*pow3(mU32) - 45*pow3(m32)*pow3(mU32) - 204*pow3(mQ32)*pow3(
         mU32) - 586*mQ32*mU32*pow4(m32) + 269*pow2(mQ32)*pow4(m32) + 73*pow2(
         mU32)*pow4(m32) - 517*m32*mU32*pow4(mQ32) + 147*pow2(m32)*pow4(mQ32) +
         178*pow2(mU32)*pow4(mQ32) - 18*m32*mQ32*pow4(mU32) - 6*pow2(mQ32)*
         pow4(mU32) + 112*mQ32*pow5(m32) - 28*mU32*pow5(m32) + 36*m32*pow5(mQ32
         ) + 12*mU32*pow5(mQ32)))/(3.*pow2(m32 - mU32)*pow2(mQ32 - mU32)*pow3(
         m32 - mQ32)) + log(msq2/MR2)*((-8*(-90*m32*mQ32*msq2 - 65*mQ32*pow2(
         m32) + 57*m32*pow2(mQ32) - 30*msq2*pow2(mQ32) + 27*pow3(m32) - 19*pow3
         (mQ32)))/(3.*pow3(m32 - mQ32)) + pow2(Xt)*((16*(90*m32*mQ32*msq2 + 2*
         mQ32*pow2(m32) - 3*m32*pow2(mQ32) + 30*msq2*pow2(mQ32) + pow3(mQ32)))/
         (3.*(mQ32 - mU32)*pow3(m32 - mQ32)) - (48*((3*mQ32 - mU32)*pow2(m32) +
         m32*(4*mQ32*mU32 - 8*pow2(mQ32)) - 2*mU32*pow2(mQ32) + 4*pow3(mQ32)))
         /(pow2(m32 - mQ32)*pow2(mQ32 - mU32))) + ((-64*m3*mQ32*(m32 - mQ32 +
         60*msq2))/(3.*pow2(m32 - mQ32)*pow2(mQ32 - mU32)) + (96*m3*(3*mQ32*
         mU32 - m32*(3*mQ32 + mU32) + 2*pow2(m32) - pow2(mQ32)))/((m32 - mQ32)*
         pow3(mQ32 - mU32)))*pow3(Xt) + ((-8*(mQ32 + mU32)*(90*m32*mQ32*msq2 +
         2*mQ32*pow2(m32) - 3*m32*pow2(mQ32) + 30*msq2*pow2(mQ32) + pow3(mQ32))
         )/(3.*pow3(m32 - mQ32)*pow3(mQ32 - mU32)) + (24*(8*mQ32*mU32*pow2(m32)
         - pow2(mQ32)*pow2(mU32) + 2*(mQ32 - mU32)*pow3(m32) + 4*mU32*pow3(
         mQ32) - 2*m32*(5*mU32*pow2(mQ32) - mQ32*pow2(mU32) + 4*pow3(mQ32)) + 5
         *pow4(mQ32)))/(pow2(m32 - mQ32)*pow4(mQ32 - mU32)))*pow4(Xt) + (64*Xt*
         (-30*m3*mQ32*msq2 + 5*m3*pow2(mQ32) - 14*mQ32*pow3(m3) + 9*pow5(m3)))/
         (3.*(mQ32 - mU32)*pow2(m32 - mQ32)) + (320*m3*mQ32*(m32 - mQ32 + 6*
         msq2)*(mQ32 + mU32)*pow5(Xt))/(3.*pow2(m32 - mQ32)*pow4(mQ32 - mU32)))
         + (4*(150*msq2*pow2(m32)*pow2(mQ32)*pow2(mU32) + 600*msq2*mU32*pow2(
         mQ32)*pow3(m32) - 900*mQ32*msq2*pow2(mU32)*pow3(m32) + 2227*pow2(mQ32)
         *pow2(mU32)*pow3(m32) - 750*msq2*mU32*pow2(m32)*pow3(mQ32) + 600*m32*
         msq2*pow2(mU32)*pow3(mQ32) - 1015*pow2(m32)*pow2(mU32)*pow3(mQ32) +
         300*msq2*pow3(m32)*pow3(mQ32) - 4599*mU32*pow3(m32)*pow3(mQ32) + 450*
         mQ32*msq2*pow2(m32)*pow3(mU32) - 300*m32*msq2*pow2(mQ32)*pow3(mU32) -
         1585*pow2(m32)*pow2(mQ32)*pow3(mU32) + 327*mQ32*pow3(m32)*pow3(mU32) +
         1693*m32*pow3(mQ32)*pow3(mU32) - 150*msq2*pow3(mQ32)*pow3(mU32) + 450
         *mQ32*msq2*mU32*pow4(m32) - 450*msq2*pow2(mQ32)*pow4(m32) + 1744*mU32*
         pow2(mQ32)*pow4(m32) - 832*mQ32*pow2(mU32)*pow4(m32) + 2856*pow3(mQ32)
         *pow4(m32) + 72*pow3(mU32)*pow4(m32) - 300*m32*msq2*mU32*pow4(mQ32) +
         150*msq2*pow2(m32)*pow4(mQ32) + 4347*mU32*pow2(m32)*pow4(mQ32) - 1101*
         m32*pow2(mU32)*pow4(mQ32) + 150*msq2*pow2(mU32)*pow4(mQ32) - 1795*pow3
         (m32)*pow4(mQ32) - 371*pow3(mU32)*pow4(mQ32) + 45*mQ32*pow2(m32)*pow4(
         mU32) - 30*m32*pow2(mQ32)*pow4(mU32) - 15*pow3(mQ32)*pow4(mU32) - 58*
         mQ32*mU32*pow5(m32) - 1816*pow2(mQ32)*pow5(m32) - 46*pow2(mU32)*pow5(
         m32) - 1027*m32*mU32*pow5(mQ32) + 128*pow2(m32)*pow5(mQ32) + 359*pow2(
         mU32)*pow5(mQ32) + 410*mQ32*pow6(m32) - 26*mU32*pow6(m32) + 81*m32*
         pow6(mQ32) + 27*mU32*pow6(mQ32)))/(3.*(mQ32 - mU32)*pow2(m32 - mU32)*
         pow4(m32 - mQ32)) - (2560*m32*pow2(mQ32)*pow6(Xt))/(3.*pow2(m32 - mQ32
         )*pow4(mQ32 - mU32)) - (8*pow4(Xt)*(-420*msq2*pow2(mQ32)*pow2(mU32)*
         pow3(m32) - 60*msq2*pow2(m32)*pow2(mU32)*pow3(mQ32) + 420*msq2*mU32*
         pow3(m32)*pow3(mQ32) + 8928*pow2(mU32)*pow3(m32)*pow3(mQ32) + 300*msq2
         *pow2(m32)*pow2(mQ32)*pow3(mU32) - 180*mQ32*msq2*pow3(m32)*pow3(mU32)
         + 3208*pow2(mQ32)*pow3(m32)*pow3(mU32) - 60*m32*msq2*pow3(mQ32)*pow3(
         mU32) - 4220*pow2(m32)*pow3(mQ32)*pow3(mU32) + 180*msq2*mU32*pow2(mQ32
         )*pow4(m32) + 90*mQ32*msq2*pow2(mU32)*pow4(m32) - 4820*pow2(mQ32)*pow2
         (mU32)*pow4(m32) - 270*msq2*pow3(mQ32)*pow4(m32) - 9076*mU32*pow3(mQ32
         )*pow4(m32) - 756*mQ32*pow3(mU32)*pow4(m32) - 420*msq2*mU32*pow2(m32)*
         pow4(mQ32) + 300*m32*msq2*pow2(mU32)*pow4(mQ32) - 5316*pow2(m32)*pow2(
         mU32)*pow4(mQ32) + 180*msq2*pow3(m32)*pow4(mQ32) + 5100*mU32*pow3(m32)
         *pow4(mQ32) + 2137*m32*pow3(mU32)*pow4(mQ32) - 60*msq2*pow3(mU32)*pow4
         (mQ32) - 1953*pow4(m32)*pow4(mQ32) + 90*mQ32*msq2*pow2(m32)*pow4(mU32)
         - 60*m32*msq2*pow2(mQ32)*pow4(mU32) - 779*pow2(m32)*pow2(mQ32)*pow4(
         mU32) + 184*mQ32*pow3(m32)*pow4(mU32) + 766*m32*pow3(mQ32)*pow4(mU32)
         - 30*msq2*pow3(mQ32)*pow4(mU32) + 45*pow4(m32)*pow4(mU32) - 208*pow4(
         mQ32)*pow4(mU32) + 4053*mU32*pow2(mQ32)*pow5(m32) + 671*mQ32*pow2(mU32
         )*pow5(m32) + 3545*pow3(mQ32)*pow5(m32) - 109*pow3(mU32)*pow5(m32) -
         180*m32*msq2*mU32*pow5(mQ32) + 90*msq2*pow2(m32)*pow5(mQ32) + 1111*
         mU32*pow2(m32)*pow5(mQ32) + 61*m32*pow2(mU32)*pow5(mQ32) + 90*msq2*
         pow2(mU32)*pow5(mQ32) - 460*pow3(m32)*pow5(mQ32) - 292*pow3(mU32)*pow5
         (mQ32) + 9*mQ32*pow2(m32)*pow5(mU32) - 6*m32*pow2(mQ32)*pow5(mU32) - 3
         *pow3(mQ32)*pow5(mU32) - 196*mQ32*mU32*pow6(m32) - 1536*pow2(mQ32)*
         pow6(m32) + 116*pow2(mU32)*pow6(m32) - 987*m32*mU32*pow6(mQ32) + 315*
         pow2(m32)*pow6(mQ32) + 408*pow2(mU32)*pow6(mQ32) + 52*mQ32*pow7(m32) -
         52*mU32*pow7(m32) + 45*m32*pow7(mQ32) + 15*mU32*pow7(mQ32)))/(3.*pow2
         (m32 - mU32)*pow4(m32 - mQ32)*pow4(mQ32 - mU32)) - (64*pow3(Xt)*(120*
         m3*msq2*mU32*pow2(mQ32) - 120*m3*mQ32*msq2*pow2(mU32) + 255*m3*pow2(
         mQ32)*pow2(mU32) + 120*mQ32*msq2*mU32*pow3(m3) - 120*msq2*pow2(mQ32)*
         pow3(m3) - 353*mU32*pow2(mQ32)*pow3(m3) - 92*mQ32*pow2(mU32)*pow3(m3)
         - 157*m3*mU32*pow3(mQ32) + 237*pow3(m3)*pow3(mQ32) - 12*m3*mQ32*pow3(
         mU32) - 24*m3*pow4(mQ32) + 197*mQ32*mU32*pow5(m3) + 18*pow2(mQ32)*pow5
         (m3) - 27*pow2(mU32)*pow5(m3) - 77*mQ32*pow7(m3) + 13*mU32*pow7(m3) +
         22*pow9(m3)))/(3.*(m32 - mU32)*pow2(m32 - mQ32)*pow3(mQ32 - mU32)) - (
         32*Xt*(150*m3*msq2*pow2(mQ32)*pow2(mU32) - 150*mQ32*msq2*pow2(mU32)*
         pow3(m3) + 287*pow2(mQ32)*pow2(mU32)*pow3(m3) - 150*m3*msq2*mU32*pow3(
         mQ32) - 266*m3*pow2(mU32)*pow3(mQ32) + 150*msq2*pow3(m3)*pow3(mQ32) +
         44*mU32*pow3(m3)*pow3(mQ32) + 15*m3*pow2(mQ32)*pow3(mU32) - 15*mQ32*
         pow3(m3)*pow3(mU32) + 240*m3*mU32*pow4(mQ32) - 364*pow3(m3)*pow4(mQ32)
         + 150*mQ32*msq2*mU32*pow5(m3) - 150*msq2*pow2(mQ32)*pow5(m3) - 331*
         mU32*pow2(mQ32)*pow5(m3) - 85*mQ32*pow2(mU32)*pow5(m3) + 464*pow3(mQ32
         )*pow5(m3) + 27*m3*pow5(mQ32) + 223*mQ32*mU32*pow7(m3) - 223*pow2(mQ32
         )*pow7(m3) - 16*pow2(mU32)*pow7(m3) + 16*mQ32*pow9(m3) - 16*mU32*pow9(
         m3)))/(3.*(m32 - mU32)*pow2(mQ32 - mU32)*pow3(m32 - mQ32)) - (64*pow5(
         Xt)*(30*m3*msq2*pow2(mQ32)*pow2(mU32) - 120*msq2*mU32*pow2(mQ32)*pow3(
         m3) - 30*mQ32*msq2*pow2(mU32)*pow3(m3) + 91*pow2(mQ32)*pow2(mU32)*pow3
         (m3) + 90*m3*msq2*mU32*pow3(mQ32) - 71*m3*pow2(mU32)*pow3(mQ32) - 90*
         msq2*pow3(m3)*pow3(mQ32) + 463*mU32*pow3(m3)*pow3(mQ32) + 3*m3*pow2(
         mQ32)*pow3(mU32) - 3*mQ32*pow3(m3)*pow3(mU32) - 271*m3*mU32*pow4(mQ32)
         + 315*pow3(m3)*pow4(mQ32) + 30*mQ32*msq2*mU32*pow5(m3) + 90*msq2*pow2
         (mQ32)*pow5(m3) - 168*mU32*pow2(mQ32)*pow5(m3) - 20*mQ32*pow2(mU32)*
         pow5(m3) - 418*pow3(mQ32)*pow5(m3) - 15*m3*pow5(mQ32) - 8*mQ32*mU32*
         pow7(m3) + 54*pow2(mQ32)*pow7(m3) - 16*pow2(mU32)*pow7(m3) + 48*mQ32*
         pow9(m3) + 16*mU32*pow9(m3)))/(3.*(m32 - mU32)*pow3(m32 - mQ32)*pow4(
         mQ32 - mU32)) + log(mU32/MR2)*((-1280*m32*(mQ32 + mU32)*(2*m32*mQ32*
         mU32 + m32*pow2(mQ32) - 3*mU32*pow2(mQ32))*pow6(Xt))/(3.*(m32 - mU32)*
         pow2(m32 - mQ32)*pow5(mQ32 - mU32)) - (4*(-90*msq2*pow2(mQ32)*pow2(
         mU32)*pow3(m32) + 270*msq2*pow2(m32)*pow2(mU32)*pow3(mQ32) - 210*msq2*
         mU32*pow3(m32)*pow3(mQ32) + 2004*pow2(mU32)*pow3(m32)*pow3(mQ32) - 90*
         msq2*pow2(m32)*pow2(mQ32)*pow3(mU32) + 270*mQ32*msq2*pow3(m32)*pow3(
         mU32) - 1420*pow2(mQ32)*pow3(m32)*pow3(mU32) - 150*m32*msq2*pow3(mQ32)
         *pow3(mU32) + 704*pow2(m32)*pow3(mQ32)*pow3(mU32) + 210*msq2*mU32*pow2
         (mQ32)*pow4(m32) - 270*mQ32*msq2*pow2(mU32)*pow4(m32) - 639*pow2(mQ32)
         *pow2(mU32)*pow4(m32) + 60*msq2*pow3(mQ32)*pow4(m32) - 2477*mU32*pow3(
         mQ32)*pow4(m32) + 807*mQ32*pow3(mU32)*pow4(m32) - 90*msq2*mU32*pow2(
         m32)*pow4(mQ32) + 90*m32*msq2*pow2(mU32)*pow4(mQ32) - 1912*pow2(m32)*
         pow2(mU32)*pow4(mQ32) + 30*msq2*pow3(m32)*pow4(mQ32) + 1604*mU32*pow3(
         m32)*pow4(mQ32) + 294*m32*pow3(mU32)*pow4(mQ32) - 30*msq2*pow3(mU32)*
         pow4(mQ32) - 349*pow4(m32)*pow4(mQ32) - 90*mQ32*msq2*pow2(m32)*pow4(
         mU32) + 60*m32*msq2*pow2(mQ32)*pow4(mU32) + 685*pow2(m32)*pow2(mQ32)*
         pow4(mU32) - 290*mQ32*pow3(m32)*pow4(mU32) - 608*m32*pow3(mQ32)*pow4(
         mU32) + 30*msq2*pow3(mQ32)*pow4(mU32) + 98*pow4(m32)*pow4(mU32) + 123*
         pow4(mQ32)*pow4(mU32) + 90*mQ32*msq2*mU32*pow5(m32) - 90*msq2*pow2(
         mQ32)*pow5(m32) + 1646*mU32*pow2(mQ32)*pow5(m32) - 40*mQ32*pow2(mU32)*
         pow5(m32) + 608*pow3(mQ32)*pow5(m32) - 294*pow3(mU32)*pow5(m32) - 227*
         mU32*pow2(m32)*pow5(mQ32) + 448*m32*pow2(mU32)*pow5(mQ32) + 22*pow3(
         m32)*pow5(mQ32) - 123*pow3(mU32)*pow5(mQ32) - 27*mQ32*pow2(m32)*pow5(
         mU32) + 18*m32*pow2(mQ32)*pow5(mU32) + 9*pow3(mQ32)*pow5(mU32) - 554*
         mQ32*mU32*pow6(m32) - 410*pow2(mQ32)*pow6(m32) + 196*pow2(mU32)*pow6(
         m32) - 24*m32*mU32*pow6(mQ32) + 9*pow2(m32)*pow6(mQ32) - 9*pow2(mU32)*
         pow6(mQ32) + 128*mQ32*pow7(m32)))/(3.*(mQ32 - mU32)*pow3(m32 - mU32)*
         pow4(m32 - mQ32)) + (8*pow2(Xt)*(-90*msq2*pow2(mQ32)*pow2(mU32)*pow3(
         m32) + 270*msq2*pow2(m32)*pow2(mU32)*pow3(mQ32) - 210*msq2*mU32*pow3(
         m32)*pow3(mQ32) - 832*pow2(mU32)*pow3(m32)*pow3(mQ32) - 90*msq2*pow2(
         m32)*pow2(mQ32)*pow3(mU32) + 270*mQ32*msq2*pow3(m32)*pow3(mU32) - 640*
         pow2(mQ32)*pow3(m32)*pow3(mU32) - 150*m32*msq2*pow3(mQ32)*pow3(mU32) -
         36*pow2(m32)*pow3(mQ32)*pow3(mU32) + 210*msq2*mU32*pow2(mQ32)*pow4(
         m32) - 270*mQ32*msq2*pow2(mU32)*pow4(m32) + 3113*pow2(mQ32)*pow2(mU32)
         *pow4(m32) + 60*msq2*pow3(mQ32)*pow4(m32) - 941*mU32*pow3(mQ32)*pow4(
         m32) + 767*mQ32*pow3(mU32)*pow4(m32) - 90*msq2*mU32*pow2(m32)*pow4(
         mQ32) + 90*m32*msq2*pow2(mU32)*pow4(mQ32) - 292*pow2(m32)*pow2(mU32)*
         pow4(mQ32) + 30*msq2*pow3(m32)*pow4(mQ32) + 1152*mU32*pow3(m32)*pow4(
         mQ32) + 222*m32*pow3(mU32)*pow4(mQ32) - 30*msq2*pow3(mU32)*pow4(mQ32)
         - 749*pow4(m32)*pow4(mQ32) - 90*mQ32*msq2*pow2(m32)*pow4(mU32) + 60*
         m32*msq2*pow2(mQ32)*pow4(mU32) - 535*pow2(m32)*pow2(mQ32)*pow4(mU32) +
         338*mQ32*pow3(m32)*pow4(mU32) + 348*m32*pow3(mQ32)*pow4(mU32) + 30*
         msq2*pow3(mQ32)*pow4(mU32) - 50*pow4(m32)*pow4(mU32) - 93*pow4(mQ32)*
         pow4(mU32) + 90*mQ32*msq2*mU32*pow5(m32) - 90*msq2*pow2(mQ32)*pow5(m32
         ) - 1130*mU32*pow2(mQ32)*pow5(m32) - 2776*mQ32*pow2(mU32)*pow5(m32) +
         1224*pow3(mQ32)*pow5(m32) - 362*pow3(mU32)*pow5(m32) - 143*mU32*pow2(
         m32)*pow5(mQ32) - 76*m32*pow2(mU32)*pow5(mQ32) + 82*pow3(m32)*pow5(
         mQ32) + 17*pow3(mU32)*pow5(mQ32) + 9*mQ32*pow2(m32)*pow5(mU32) - 6*m32
         *pow2(mQ32)*pow5(mU32) - 3*pow3(mQ32)*pow5(mU32) + 1514*mQ32*mU32*pow6
         (m32) - 682*pow2(mQ32)*pow6(m32) + 908*pow2(mU32)*pow6(m32) + 12*m32*
         mU32*pow6(mQ32) + 9*pow2(m32)*pow6(mQ32) + 3*pow2(mU32)*pow6(mQ32) +
         124*mQ32*pow7(m32) - 496*mU32*pow7(m32)))/(3.*pow2(mQ32 - mU32)*pow3(
         m32 - mU32)*pow4(m32 - mQ32)) + (4*pow4(Xt)*(300*msq2*pow2(mU32)*pow3(
         m32)*pow3(mQ32) - 180*msq2*pow2(mQ32)*pow3(m32)*pow3(mU32) - 180*msq2*
         pow2(m32)*pow3(mQ32)*pow3(mU32) + 14272*pow3(m32)*pow3(mQ32)*pow3(mU32
         ) + 60*msq2*pow2(mQ32)*pow2(mU32)*pow4(m32) - 270*msq2*mU32*pow3(mQ32)
         *pow4(m32) - 9228*pow2(mU32)*pow3(mQ32)*pow4(m32) + 270*mQ32*msq2*pow3
         (mU32)*pow4(m32) - 20318*pow2(mQ32)*pow3(mU32)*pow4(m32) - 180*msq2*
         pow2(m32)*pow2(mU32)*pow4(mQ32) + 180*msq2*mU32*pow3(m32)*pow4(mQ32) -
         6214*pow2(mU32)*pow3(m32)*pow4(mQ32) + 60*m32*msq2*pow3(mU32)*pow4(
         mQ32) - 526*pow2(m32)*pow3(mU32)*pow4(mQ32) - 60*msq2*pow4(m32)*pow4(
         mQ32) + 7810*mU32*pow4(m32)*pow4(mQ32) + 180*msq2*pow2(m32)*pow2(mQ32)
         *pow4(mU32) - 270*mQ32*msq2*pow3(m32)*pow4(mU32) + 5386*pow2(mQ32)*
         pow3(m32)*pow4(mU32) + 90*m32*msq2*pow3(mQ32)*pow4(mU32) - 4933*pow2(
         m32)*pow3(mQ32)*pow4(mU32) - 2441*mQ32*pow4(m32)*pow4(mU32) + 1352*m32
         *pow4(mQ32)*pow4(mU32) - 90*mQ32*msq2*pow2(mU32)*pow5(m32) + 20180*
         pow2(mQ32)*pow2(mU32)*pow5(m32) + 90*msq2*pow3(mQ32)*pow5(m32) - 838*
         mU32*pow3(mQ32)*pow5(m32) + 10882*mQ32*pow3(mU32)*pow5(m32) - 2556*
         pow4(mQ32)*pow5(m32) + 508*pow4(mU32)*pow5(m32) + 90*msq2*mU32*pow2(
         m32)*pow5(mQ32) - 90*m32*msq2*pow2(mU32)*pow5(mQ32) + 5003*pow2(m32)*
         pow2(mU32)*pow5(mQ32) - 30*msq2*pow3(m32)*pow5(mQ32) - 4626*mU32*pow3(
         m32)*pow5(mQ32) - 2162*m32*pow3(mU32)*pow5(mQ32) + 30*msq2*pow3(mU32)*
         pow5(mQ32) + 1513*pow4(m32)*pow5(mQ32) + 152*pow4(mU32)*pow5(mQ32) +
         90*mQ32*msq2*pow2(m32)*pow5(mU32) - 60*m32*msq2*pow2(mQ32)*pow5(mU32)
         + 690*pow2(m32)*pow2(mQ32)*pow5(mU32) - 362*mQ32*pow3(m32)*pow5(mU32)
         - 494*m32*pow3(mQ32)*pow5(mU32) - 30*msq2*pow3(mQ32)*pow5(mU32) + 24*
         pow4(m32)*pow5(mU32) + 134*pow4(mQ32)*pow5(mU32) - 5894*mU32*pow2(mQ32
         )*pow6(m32) - 12382*mQ32*pow2(mU32)*pow6(m32) + 1338*pow3(mQ32)*pow6(
         m32) - 2182*pow3(mU32)*pow6(m32) + 280*mU32*pow2(m32)*pow6(mQ32) - 50*
         m32*pow2(mU32)*pow6(mQ32) - 152*pow3(m32)*pow6(mQ32) + 18*pow3(mU32)*
         pow6(mQ32) - 9*mQ32*pow2(m32)*pow6(mU32) + 6*m32*pow2(mQ32)*pow6(mU32)
         + 3*pow3(mQ32)*pow6(mU32) + 4332*mQ32*mU32*pow7(m32) - 146*pow2(mQ32)
         *pow7(m32) + 2678*pow2(mU32)*pow7(m32) - 12*m32*mU32*pow7(mQ32) - 9*
         pow2(m32)*pow7(mQ32) - 3*pow2(mU32)*pow7(mQ32) + 4*mQ32*pow8(m32) -
         1028*mU32*pow8(m32)))/(3.*pow3(m32 - mU32)*pow4(m32 - mQ32)*pow4(mQ32
         - mU32)) + (32*Xt*(2*mQ32*pow11(m3) - 2*mU32*pow11(m3) + 30*msq2*pow2(
         mQ32)*pow2(mU32)*pow3(m3) + 30*m3*msq2*pow2(mU32)*pow3(mQ32) - 60*msq2
         *mU32*pow3(m3)*pow3(mQ32) - 164*pow2(mU32)*pow3(m3)*pow3(mQ32) - 30*m3
         *msq2*pow2(mQ32)*pow3(mU32) + 30*mQ32*msq2*pow3(m3)*pow3(mU32) - 55*
         pow2(mQ32)*pow3(m3)*pow3(mU32) + 60*m3*pow3(mQ32)*pow3(mU32) - 26*m3*
         pow2(mU32)*pow4(mQ32) + 143*mU32*pow3(m3)*pow4(mQ32) - 9*m3*pow2(mQ32)
         *pow4(mU32) + 9*mQ32*pow3(m3)*pow4(mU32) + 30*msq2*mU32*pow2(mQ32)*
         pow5(m3) - 60*mQ32*msq2*pow2(mU32)*pow5(m3) + 329*pow2(mQ32)*pow2(mU32
         )*pow5(m3) + 30*msq2*pow3(mQ32)*pow5(m3) - 120*mU32*pow3(mQ32)*pow5(m3
         ) - 36*mQ32*pow3(mU32)*pow5(m3) - 77*pow4(mQ32)*pow5(m3) - 9*m3*mU32*
         pow5(mQ32) + 3*pow3(m3)*pow5(mQ32) + 30*mQ32*msq2*mU32*pow7(m3) - 30*
         msq2*pow2(mQ32)*pow7(m3) - 77*mU32*pow2(mQ32)*pow7(m3) - 110*mQ32*pow2
         (mU32)*pow7(m3) + 108*pow3(mQ32)*pow7(m3) + 15*pow3(mU32)*pow7(m3) +
         17*mQ32*mU32*pow9(m3) - 20*pow2(mQ32)*pow9(m3) + 19*pow2(mU32)*pow9(m3
         )))/(3.*pow2(m32 - mU32)*pow2(mQ32 - mU32)*pow3(m32 - mQ32)) + (32*
         pow3(Xt)*(-7*mQ32*pow11(m3) - 69*mU32*pow11(m3) + 2*pow13(m3) - 60*
         msq2*pow2(mQ32)*pow2(mU32)*pow3(m3) - 60*m3*msq2*pow2(mU32)*pow3(mQ32)
         + 120*msq2*mU32*pow3(m3)*pow3(mQ32) - 875*pow2(mU32)*pow3(m3)*pow3(
         mQ32) + 60*m3*msq2*pow2(mQ32)*pow3(mU32) - 60*mQ32*msq2*pow3(m3)*pow3(
         mU32) + 19*pow2(mQ32)*pow3(m3)*pow3(mU32) - 15*m3*pow3(mQ32)*pow3(mU32
         ) + 349*m3*pow2(mU32)*pow4(mQ32) - 500*mU32*pow3(m3)*pow4(mQ32) - 6*m3
         *pow2(mQ32)*pow4(mU32) + 6*mQ32*pow3(m3)*pow4(mU32) - 60*msq2*mU32*
         pow2(mQ32)*pow5(m3) + 120*mQ32*msq2*pow2(mU32)*pow5(m3) + 857*pow2(
         mQ32)*pow2(mU32)*pow5(m3) - 60*msq2*pow3(mQ32)*pow5(m3) + 1123*mU32*
         pow3(mQ32)*pow5(m3) + 11*mQ32*pow3(mU32)*pow5(m3) + 215*pow4(mQ32)*
         pow5(m3) - 6*m3*mU32*pow5(mQ32) - 6*pow3(m3)*pow5(mQ32) - 60*mQ32*msq2
         *mU32*pow7(m3) + 60*msq2*pow2(mQ32)*pow7(m3) - 923*mU32*pow2(mQ32)*
         pow7(m3) - 461*mQ32*pow2(mU32)*pow7(m3) - 337*pow3(mQ32)*pow7(m3) + 17
         *pow3(mU32)*pow7(m3) + 471*mQ32*mU32*pow9(m3) + 101*pow2(mQ32)*pow9(m3
         ) + 34*pow2(mU32)*pow9(m3)))/(3.*pow2(m32 - mU32)*pow3(m32 - mQ32)*
         pow3(mQ32 - mU32)) - (32*pow5(Xt)*(32*mQ32*mU32*pow11(m3) - 32*pow11(
         m3)*pow2(mU32) + 30*msq2*pow2(mU32)*pow3(m3)*pow3(mQ32) - 60*msq2*pow2
         (mQ32)*pow3(m3)*pow3(mU32) - 679*pow3(m3)*pow3(mQ32)*pow3(mU32) - 30*
         m3*msq2*pow2(mU32)*pow4(mQ32) + 60*msq2*mU32*pow3(m3)*pow4(mQ32) - 851
         *pow2(mU32)*pow3(m3)*pow4(mQ32) + 246*m3*pow3(mU32)*pow4(mQ32) + 30*m3
         *msq2*pow2(mQ32)*pow4(mU32) - 30*mQ32*msq2*pow3(m3)*pow4(mU32) - 74*
         pow2(mQ32)*pow3(m3)*pow4(mU32) + 29*m3*pow3(mQ32)*pow4(mU32) + 30*msq2
         *pow2(mQ32)*pow2(mU32)*pow5(m3) - 60*msq2*mU32*pow3(mQ32)*pow5(m3) +
         1177*pow2(mU32)*pow3(mQ32)*pow5(m3) + 60*mQ32*msq2*pow3(mU32)*pow5(m3)
         + 713*pow2(mQ32)*pow3(mU32)*pow5(m3) - 30*msq2*pow4(mQ32)*pow5(m3) +
         791*mU32*pow4(mQ32)*pow5(m3) + 60*mQ32*pow4(mU32)*pow5(m3) + 211*m3*
         pow2(mU32)*pow5(mQ32) - 316*mU32*pow3(m3)*pow5(mQ32) + 139*pow5(m3)*
         pow5(mQ32) - 3*m3*pow2(mQ32)*pow5(mU32) + 3*mQ32*pow3(m3)*pow5(mU32) -
         3*m3*mU32*pow6(mQ32) - 3*pow3(m3)*pow6(mQ32) - 30*mQ32*msq2*pow2(mU32
         )*pow7(m3) - 731*pow2(mQ32)*pow2(mU32)*pow7(m3) + 30*msq2*pow3(mQ32)*
         pow7(m3) - 617*mU32*pow3(mQ32)*pow7(m3) - 343*mQ32*pow3(mU32)*pow7(m3)
         - 230*pow4(mQ32)*pow7(m3) + pow4(mU32)*pow7(m3) + 145*mU32*pow2(mQ32)
         *pow9(m3) + 226*mQ32*pow2(mU32)*pow9(m3) + 78*pow3(mQ32)*pow9(m3) + 31
         *pow3(mU32)*pow9(m3)))/(3.*pow2(m32 - mU32)*pow3(m32 - mQ32)*pow5(mQ32
         - mU32)))) + pow2(log(m32/MR2))*(pow2(Xt)*((-1792*pow3(m32))/(3.*pow2
         (m32 - mQ32)*pow2(m32 - mU32)) + (64*(-m32 + mQ32 + mU32)*pow2(m32)*(2
         + (8*(pow2(m32)*pow2(mQ32) + pow2(m32)*pow2(mU32) - 2*mQ32*pow3(m32)
         - 2*mU32*pow3(m32) + 2*pow4(m32)))/(3.*pow2(m32 - mQ32)*pow2(m32 -
         mU32))))/(mQ32*(-m32 + mQ32)*(m32 - mU32)*mU32)) - (2560*pow3(m32)*
         pow6(Xt))/(3.*pow2(m32 - mQ32)*pow2(m32 - mU32)*pow2(mQ32 - mU32)) - (
         2048*(m32 - mQ32 - mU32)*pow3(Xt)*pow7(m3))/(3.*mQ32*mU32*pow2(m32 -
         mQ32)*pow2(m32 - mU32)) + (64*pow5(Xt)*(34*mQ32*mU32*pow11(m3) + 16*
         mQ32*pow13(m3) + 16*mU32*pow13(m3) - 32*pow11(m3)*pow2(mQ32) - 32*
         pow11(m3)*pow2(mU32) - 30*msq2*pow2(mU32)*pow3(m3)*pow3(mQ32) - 30*
         msq2*pow2(mQ32)*pow3(m3)*pow3(mU32) + 260*pow3(m3)*pow3(mQ32)*pow3(
         mU32) - 3*pow2(mU32)*pow3(m3)*pow4(mQ32) - 3*pow2(mQ32)*pow3(m3)*pow4(
         mU32) + 120*msq2*pow2(mQ32)*pow2(mU32)*pow5(m3) + 30*msq2*mU32*pow3(
         mQ32)*pow5(m3) - 401*pow2(mU32)*pow3(mQ32)*pow5(m3) + 30*mQ32*msq2*
         pow3(mU32)*pow5(m3) - 401*pow2(mQ32)*pow3(mU32)*pow5(m3) + 3*mU32*pow4
         (mQ32)*pow5(m3) + 3*mQ32*pow4(mU32)*pow5(m3) - 90*msq2*mU32*pow2(mQ32)
         *pow7(m3) - 90*mQ32*msq2*pow2(mU32)*pow7(m3) + 582*pow2(mQ32)*pow2(
         mU32)*pow7(m3) + 141*mU32*pow3(mQ32)*pow7(m3) + 141*mQ32*pow3(mU32)*
         pow7(m3) + 60*mQ32*msq2*mU32*pow9(m3) - 178*mU32*pow2(mQ32)*pow9(m3) -
         178*mQ32*pow2(mU32)*pow9(m3) + 16*pow3(mQ32)*pow9(m3) + 16*pow3(mU32)
         *pow9(m3)))/(3.*mQ32*mU32*pow2(mQ32 - mU32)*pow3(m32 - mQ32)*pow3(m32
         - mU32)) + (32*Xt*(42*mQ32*mU32*pow11(m3) + 32*mQ32*pow13(m3) + 32*
         mU32*pow13(m3) - 64*pow11(m3)*pow2(mQ32) - 64*pow11(m3)*pow2(mU32) +
         30*msq2*pow2(mU32)*pow3(m3)*pow3(mQ32) + 30*msq2*pow2(mQ32)*pow3(m3)*
         pow3(mU32) - 248*pow3(m3)*pow3(mQ32)*pow3(mU32) + 3*pow2(mU32)*pow3(m3
         )*pow4(mQ32) + 3*pow2(mQ32)*pow3(m3)*pow4(mU32) - 120*msq2*pow2(mQ32)*
         pow2(mU32)*pow5(m3) - 30*msq2*mU32*pow3(mQ32)*pow5(m3) + 393*pow2(mU32
         )*pow3(mQ32)*pow5(m3) - 30*mQ32*msq2*pow3(mU32)*pow5(m3) + 393*pow2(
         mQ32)*pow3(mU32)*pow5(m3) - 3*mU32*pow4(mQ32)*pow5(m3) - 3*mQ32*pow4(
         mU32)*pow5(m3) + 90*msq2*mU32*pow2(mQ32)*pow7(m3) + 90*mQ32*msq2*pow2(
         mU32)*pow7(m3) - 598*pow2(mQ32)*pow2(mU32)*pow7(m3) - 97*mU32*pow3(
         mQ32)*pow7(m3) - 97*mQ32*pow3(mU32)*pow7(m3) - 60*mQ32*msq2*mU32*pow9(
         m3) + 106*mU32*pow2(mQ32)*pow9(m3) + 106*mQ32*pow2(mU32)*pow9(m3) + 32
         *pow3(mQ32)*pow9(m3) + 32*pow3(mU32)*pow9(m3)))/(3.*mQ32*mU32*pow3(m32
         - mQ32)*pow3(m32 - mU32)) - (8*pow4(Xt)*(-360*msq2*pow3(m32)*pow3(
         mQ32)*pow3(mU32) - 240*msq2*pow2(mU32)*pow3(mQ32)*pow4(m32) - 240*msq2
         *pow2(mQ32)*pow3(mU32)*pow4(m32) + 7264*pow3(mQ32)*pow3(mU32)*pow4(m32
         ) + 240*msq2*pow2(mU32)*pow3(m32)*pow4(mQ32) + 120*msq2*pow2(m32)*pow3
         (mU32)*pow4(mQ32) - 1628*pow3(m32)*pow3(mU32)*pow4(mQ32) - 360*msq2*
         mU32*pow4(m32)*pow4(mQ32) + 2380*pow2(mU32)*pow4(m32)*pow4(mQ32) + 240
         *msq2*pow2(mQ32)*pow3(m32)*pow4(mU32) + 120*msq2*pow2(m32)*pow3(mQ32)*
         pow4(mU32) - 1628*pow3(m32)*pow3(mQ32)*pow4(mU32) - 360*mQ32*msq2*pow4
         (m32)*pow4(mU32) + 2380*pow2(mQ32)*pow4(m32)*pow4(mU32) - 768*pow2(m32
         )*pow4(mQ32)*pow4(mU32) + 480*msq2*pow2(mQ32)*pow2(mU32)*pow5(m32) +
         510*msq2*mU32*pow3(mQ32)*pow5(m32) - 7189*pow2(mU32)*pow3(mQ32)*pow5(
         m32) + 510*mQ32*msq2*pow3(mU32)*pow5(m32) - 7189*pow2(mQ32)*pow3(mU32)
         *pow5(m32) - 77*mU32*pow4(mQ32)*pow5(m32) - 77*mQ32*pow4(mU32)*pow5(
         m32) - 60*msq2*pow2(m32)*pow2(mU32)*pow5(mQ32) + 90*msq2*mU32*pow3(m32
         )*pow5(mQ32) - 323*pow2(mU32)*pow3(m32)*pow5(mQ32) - 30*m32*msq2*pow3(
         mU32)*pow5(mQ32) + 148*pow2(m32)*pow3(mU32)*pow5(mQ32) - 86*mU32*pow4(
         m32)*pow5(mQ32) + 285*m32*pow4(mU32)*pow5(mQ32) + 56*pow5(m32)*pow5(
         mQ32) - 60*msq2*pow2(m32)*pow2(mQ32)*pow5(mU32) + 90*mQ32*msq2*pow3(
         m32)*pow5(mU32) - 323*pow2(mQ32)*pow3(m32)*pow5(mU32) - 30*m32*msq2*
         pow3(mQ32)*pow5(mU32) + 148*pow2(m32)*pow3(mQ32)*pow5(mU32) - 86*mQ32*
         pow4(m32)*pow5(mU32) + 285*m32*pow4(mQ32)*pow5(mU32) + 56*pow5(m32)*
         pow5(mU32) - 72*pow5(mQ32)*pow5(mU32) - 420*msq2*mU32*pow2(mQ32)*pow6(
         m32) - 420*mQ32*msq2*pow2(mU32)*pow6(m32) + 6516*pow2(mQ32)*pow2(mU32)
         *pow6(m32) + 1088*mU32*pow3(mQ32)*pow6(m32) + 1088*mQ32*pow3(mU32)*
         pow6(m32) - 224*pow4(mQ32)*pow6(m32) - 224*pow4(mU32)*pow6(m32) - 6*
         pow2(m32)*pow2(mU32)*pow6(mQ32) + 9*mU32*pow3(m32)*pow6(mQ32) - 3*m32*
         pow3(mU32)*pow6(mQ32) - 6*pow2(m32)*pow2(mQ32)*pow6(mU32) + 9*mQ32*
         pow3(m32)*pow6(mU32) - 3*m32*pow3(mQ32)*pow6(mU32) + 180*mQ32*msq2*
         mU32*pow7(m32) - 1122*mU32*pow2(mQ32)*pow7(m32) - 1122*mQ32*pow2(mU32)
         *pow7(m32) + 368*pow3(mQ32)*pow7(m32) + 368*pow3(mU32)*pow7(m32) + 108
         *mQ32*mU32*pow8(m32) - 288*pow2(mQ32)*pow8(m32) - 288*pow2(mU32)*pow8(
         m32) + 88*mQ32*pow9(m32) + 88*mU32*pow9(m32)))/(3.*mQ32*mU32*pow2(mQ32
         - mU32)*pow4(m32 - mQ32)*pow4(m32 - mU32)) - (4*(360*msq2*pow3(m32)*
         pow3(mQ32)*pow3(mU32) + 240*msq2*pow2(mU32)*pow3(mQ32)*pow4(m32) + 240
         *msq2*pow2(mQ32)*pow3(mU32)*pow4(m32) - 7936*pow3(mQ32)*pow3(mU32)*
         pow4(m32) - 240*msq2*pow2(mU32)*pow3(m32)*pow4(mQ32) - 120*msq2*pow2(
         m32)*pow3(mU32)*pow4(mQ32) + 1996*pow3(m32)*pow3(mU32)*pow4(mQ32) +
         360*msq2*mU32*pow4(m32)*pow4(mQ32) - 3052*pow2(mU32)*pow4(m32)*pow4(
         mQ32) - 240*msq2*pow2(mQ32)*pow3(m32)*pow4(mU32) - 120*msq2*pow2(m32)*
         pow3(mQ32)*pow4(mU32) + 1996*pow3(m32)*pow3(mQ32)*pow4(mU32) + 360*
         mQ32*msq2*pow4(m32)*pow4(mU32) - 3052*pow2(mQ32)*pow4(m32)*pow4(mU32)
         + 576*pow2(m32)*pow4(mQ32)*pow4(mU32) - 480*msq2*pow2(mQ32)*pow2(mU32)
         *pow5(m32) - 510*msq2*mU32*pow3(mQ32)*pow5(m32) + 8341*pow2(mU32)*pow3
         (mQ32)*pow5(m32) - 510*mQ32*msq2*pow3(mU32)*pow5(m32) + 8341*pow2(mQ32
         )*pow3(mU32)*pow5(m32) + 645*mU32*pow4(mQ32)*pow5(m32) + 645*mQ32*pow4
         (mU32)*pow5(m32) + 60*msq2*pow2(m32)*pow2(mU32)*pow5(mQ32) - 90*msq2*
         mU32*pow3(m32)*pow5(mQ32) + 451*pow2(mU32)*pow3(m32)*pow5(mQ32) + 30*
         m32*msq2*pow3(mU32)*pow5(mQ32) - 196*pow2(m32)*pow3(mU32)*pow5(mQ32) -
         42*mU32*pow4(m32)*pow5(mQ32) - 261*m32*pow4(mU32)*pow5(mQ32) + 112*
         pow5(m32)*pow5(mQ32) + 60*msq2*pow2(m32)*pow2(mQ32)*pow5(mU32) - 90*
         mQ32*msq2*pow3(m32)*pow5(mU32) + 451*pow2(mQ32)*pow3(m32)*pow5(mU32) +
         30*m32*msq2*pow3(mQ32)*pow5(mU32) - 196*pow2(m32)*pow3(mQ32)*pow5(
         mU32) - 42*mQ32*pow4(m32)*pow5(mU32) - 261*m32*pow4(mQ32)*pow5(mU32) +
         112*pow5(m32)*pow5(mU32) + 72*pow5(mQ32)*pow5(mU32) + 420*msq2*mU32*
         pow2(mQ32)*pow6(m32) + 420*mQ32*msq2*pow2(mU32)*pow6(m32) - 7924*pow2(
         mQ32)*pow2(mU32)*pow6(m32) - 2128*mU32*pow3(mQ32)*pow6(m32) - 2128*
         mQ32*pow3(mU32)*pow6(m32) - 448*pow4(mQ32)*pow6(m32) - 448*pow4(mU32)*
         pow6(m32) + 6*pow2(m32)*pow2(mU32)*pow6(mQ32) - 9*mU32*pow3(m32)*pow6(
         mQ32) + 3*m32*pow3(mU32)*pow6(mQ32) + 6*pow2(m32)*pow2(mQ32)*pow6(mU32
         ) - 9*mQ32*pow3(m32)*pow6(mU32) + 3*m32*pow3(mQ32)*pow6(mU32) - 180*
         mQ32*msq2*mU32*pow7(m32) + 2210*mU32*pow2(mQ32)*pow7(m32) + 2210*mQ32*
         pow2(mU32)*pow7(m32) + 736*pow3(mQ32)*pow7(m32) + 736*pow3(mU32)*pow7(
         m32) - 716*mQ32*mU32*pow8(m32) - 576*pow2(mQ32)*pow8(m32) - 576*pow2(
         mU32)*pow8(m32) + 176*mQ32*pow9(m32) + 176*mU32*pow9(m32)))/(3.*mQ32*
         mU32*pow4(m32 - mQ32)*pow4(m32 - mU32)) + log(mU32/MR2)*((-1280*(mQ32
         + mU32)*pow3(m32)*pow6(Xt))/(3.*pow2(m32 - mQ32)*pow2(m32 - mU32)*pow3
         (mQ32 - mU32)) + (8*pow2(Xt)*(-360*msq2*pow2(mQ32)*pow2(mU32)*pow3(m32
         ) + 120*msq2*pow2(m32)*pow2(mU32)*pow3(mQ32) + 240*msq2*mU32*pow3(m32)
         *pow3(mQ32) - 1836*pow2(mU32)*pow3(m32)*pow3(mQ32) + 120*msq2*pow2(m32
         )*pow2(mQ32)*pow3(mU32) + 240*mQ32*msq2*pow3(m32)*pow3(mU32) - 2028*
         pow2(mQ32)*pow3(m32)*pow3(mU32) - 576*pow2(m32)*pow3(mQ32)*pow3(mU32)
         - 240*msq2*mU32*pow2(mQ32)*pow4(m32) - 240*mQ32*msq2*pow2(mU32)*pow4(
         m32) + 7136*pow2(mQ32)*pow2(mU32)*pow4(m32) - 360*msq2*pow3(mQ32)*pow4
         (m32) + 2540*mU32*pow3(mQ32)*pow4(m32) + 2924*mQ32*pow3(mU32)*pow4(m32
         ) - 360*msq2*pow3(mU32)*pow4(m32) - 60*msq2*mU32*pow2(m32)*pow4(mQ32)
         - 30*m32*msq2*pow2(mU32)*pow4(mQ32) + 244*pow2(m32)*pow2(mU32)*pow4(
         mQ32) + 90*msq2*pow3(m32)*pow4(mQ32) - 499*mU32*pow3(m32)*pow4(mQ32) +
         261*m32*pow3(mU32)*pow4(mQ32) + 74*pow4(m32)*pow4(mQ32) - 60*mQ32*
         msq2*pow2(m32)*pow4(mU32) - 30*m32*msq2*pow2(mQ32)*pow4(mU32) + 244*
         pow2(m32)*pow2(mQ32)*pow4(mU32) - 499*mQ32*pow3(m32)*pow4(mU32) + 90*
         msq2*pow3(m32)*pow4(mU32) + 261*m32*pow3(mQ32)*pow4(mU32) + 74*pow4(
         m32)*pow4(mU32) - 72*pow4(mQ32)*pow4(mU32) + 480*mQ32*msq2*mU32*pow5(
         m32) + 510*msq2*pow2(mQ32)*pow5(m32) - 6485*mU32*pow2(mQ32)*pow5(m32)
         - 7061*mQ32*pow2(mU32)*pow5(m32) + 510*msq2*pow2(mU32)*pow5(m32) - 421
         *pow3(mQ32)*pow5(m32) - 613*pow3(mU32)*pow5(m32) - 6*mU32*pow2(m32)*
         pow5(mQ32) - 3*m32*pow2(mU32)*pow5(mQ32) + 9*pow3(m32)*pow5(mQ32) - 6*
         mQ32*pow2(m32)*pow5(mU32) - 3*m32*pow2(mQ32)*pow5(mU32) + 9*pow3(m32)*
         pow5(mU32) - 420*mQ32*msq2*pow6(m32) + 5620*mQ32*mU32*pow6(m32) - 420*
         msq2*mU32*pow6(m32) + 1184*pow2(mQ32)*pow6(m32) + 1568*pow2(mU32)*pow6
         (m32) - 1010*mQ32*pow7(m32) + 180*msq2*pow7(m32) - 1202*mU32*pow7(m32)
         + 172*pow8(m32)))/(3.*(mQ32 - mU32)*pow4(m32 - mQ32)*pow4(m32 - mU32)
         ) - (4*(-360*msq2*pow2(mQ32)*pow2(mU32)*pow3(m32) + 120*msq2*pow2(m32)
         *pow2(mU32)*pow3(mQ32) + 240*msq2*mU32*pow3(m32)*pow3(mQ32) - 1996*
         pow2(mU32)*pow3(m32)*pow3(mQ32) + 120*msq2*pow2(m32)*pow2(mQ32)*pow3(
         mU32) + 240*mQ32*msq2*pow3(m32)*pow3(mU32) - 2188*pow2(mQ32)*pow3(m32)
         *pow3(mU32) - 576*pow2(m32)*pow3(mQ32)*pow3(mU32) - 240*msq2*mU32*pow2
         (mQ32)*pow4(m32) - 240*mQ32*msq2*pow2(mU32)*pow4(m32) + 8320*pow2(mQ32
         )*pow2(mU32)*pow4(m32) - 360*msq2*pow3(mQ32)*pow4(m32) + 2860*mU32*
         pow3(mQ32)*pow4(m32) + 3244*mQ32*pow3(mU32)*pow4(m32) - 360*msq2*pow3(
         mU32)*pow4(m32) - 60*msq2*mU32*pow2(m32)*pow4(mQ32) - 30*m32*msq2*pow2
         (mU32)*pow4(mQ32) + 196*pow2(m32)*pow2(mU32)*pow4(mQ32) + 90*msq2*pow3
         (m32)*pow4(mQ32) - 403*mU32*pow3(m32)*pow4(mQ32) + 261*m32*pow3(mU32)*
         pow4(mQ32) - 38*pow4(m32)*pow4(mQ32) - 60*mQ32*msq2*pow2(m32)*pow4(
         mU32) - 30*m32*msq2*pow2(mQ32)*pow4(mU32) + 244*pow2(m32)*pow2(mQ32)*
         pow4(mU32) - 499*mQ32*pow3(m32)*pow4(mU32) + 90*msq2*pow3(m32)*pow4(
         mU32) + 261*m32*pow3(mQ32)*pow4(mU32) + 74*pow4(m32)*pow4(mU32) - 72*
         pow4(mQ32)*pow4(mU32) + 480*mQ32*msq2*mU32*pow5(m32) + 510*msq2*pow2(
         mQ32)*pow5(m32) - 8373*mU32*pow2(mQ32)*pow5(m32) - 8757*mQ32*pow2(mU32
         )*pow5(m32) + 510*msq2*pow2(mU32)*pow5(m32) - 325*pow3(mQ32)*pow5(m32)
         - 773*pow3(mU32)*pow5(m32) - 6*mU32*pow2(m32)*pow5(mQ32) - 3*m32*pow2
         (mU32)*pow5(mQ32) + 9*pow3(m32)*pow5(mQ32) - 6*mQ32*pow2(m32)*pow5(
         mU32) - 3*m32*pow2(mQ32)*pow5(mU32) + 9*pow3(m32)*pow5(mU32) - 420*
         mQ32*msq2*pow6(m32) + 8052*mQ32*mU32*pow6(m32) - 420*msq2*mU32*pow6(
         m32) + 1664*pow2(mQ32)*pow6(m32) + 2288*pow2(mU32)*pow6(m32) - 1810*
         mQ32*pow7(m32) + 180*msq2*pow7(m32) - 2162*mU32*pow7(m32) + 508*pow8(
         m32)))/(3.*pow4(m32 - mQ32)*pow4(m32 - mU32)) + (32*pow5(Xt)*(98*mQ32*
         pow11(m3) + 98*mU32*pow11(m3) - 60*msq2*pow2(mQ32)*pow2(mU32)*pow3(m3)
         - 30*msq2*mU32*pow3(m3)*pow3(mQ32) + 257*pow2(mU32)*pow3(m3)*pow3(
         mQ32) - 30*mQ32*msq2*pow3(m3)*pow3(mU32) + 257*pow2(mQ32)*pow3(m3)*
         pow3(mU32) - 3*mU32*pow3(m3)*pow4(mQ32) - 3*mQ32*pow3(m3)*pow4(mU32) +
         150*msq2*mU32*pow2(mQ32)*pow5(m3) + 150*mQ32*msq2*pow2(mU32)*pow5(m3)
         - 834*pow2(mQ32)*pow2(mU32)*pow5(m3) + 30*msq2*pow3(mQ32)*pow5(m3) -
         414*mU32*pow3(mQ32)*pow5(m3) - 414*mQ32*pow3(mU32)*pow5(m3) + 30*msq2*
         pow3(mU32)*pow5(m3) + 3*pow4(mQ32)*pow5(m3) + 3*pow4(mU32)*pow5(m3) -
         180*mQ32*msq2*mU32*pow7(m3) - 90*msq2*pow2(mQ32)*pow7(m3) + 819*mU32*
         pow2(mQ32)*pow7(m3) + 819*mQ32*pow2(mU32)*pow7(m3) - 90*msq2*pow2(mU32
         )*pow7(m3) + 173*pow3(mQ32)*pow7(m3) + 173*pow3(mU32)*pow7(m3) + 60*
         mQ32*msq2*pow9(m3) - 516*mQ32*mU32*pow9(m3) + 60*msq2*mU32*pow9(m3) -
         258*pow2(mQ32)*pow9(m3) - 258*pow2(mU32)*pow9(m3)))/(3.*pow3(m32 -
         mQ32)*pow3(m32 - mU32)*pow3(mQ32 - mU32)) + (32*Xt*(-134*mQ32*pow11(m3
         ) - 178*mU32*pow11(m3) + 88*pow13(m3) - 30*msq2*mU32*pow3(m3)*pow3(
         mQ32) + 239*pow2(mU32)*pow3(m3)*pow3(mQ32) + 30*mQ32*msq2*pow3(m3)*
         pow3(mU32) - 263*pow2(mQ32)*pow3(m3)*pow3(mU32) - 3*mU32*pow3(m3)*pow4
         (mQ32) + 3*mQ32*pow3(m3)*pow4(mU32) + 90*msq2*mU32*pow2(mQ32)*pow5(m3)
         - 90*mQ32*msq2*pow2(mU32)*pow5(m3) + 72*pow2(mQ32)*pow2(mU32)*pow5(m3
         ) + 30*msq2*pow3(mQ32)*pow5(m3) - 372*mU32*pow3(mQ32)*pow5(m3) + 420*
         mQ32*pow3(mU32)*pow5(m3) - 30*msq2*pow3(mU32)*pow5(m3) + 3*pow4(mQ32)*
         pow5(m3) - 3*pow4(mU32)*pow5(m3) - 90*msq2*pow2(mQ32)*pow7(m3) + 361*
         mU32*pow2(mQ32)*pow7(m3) - 577*mQ32*pow2(mU32)*pow7(m3) + 90*msq2*pow2
         (mU32)*pow7(m3) + 85*pow3(mQ32)*pow7(m3) - 173*pow3(mU32)*pow7(m3) +
         60*mQ32*msq2*pow9(m3) + 144*mQ32*mU32*pow9(m3) - 60*msq2*mU32*pow9(m3)
         - 26*pow2(mQ32)*pow9(m3) + 314*pow2(mU32)*pow9(m3)))/(3.*(mQ32 - mU32
         )*pow3(m32 - mQ32)*pow3(m32 - mU32)) + pow3(Xt)*((128*(4*mQ32*mU32 - 2
         *m32*(mQ32 + mU32) + 2*pow2(m32) - pow2(mQ32) - pow2(mU32))*pow3(m3)*(
         2 + (8*(pow2(m32)*pow2(mQ32) + pow2(m32)*pow2(mU32) - 2*mQ32*pow3(m32)
         - 2*mU32*pow3(m32) + 2*pow4(m32)))/(3.*pow2(m32 - mQ32)*pow2(m32 -
         mU32))))/((m32 - mQ32)*(m32 - mU32)*pow3(mQ32 - mU32)) - (1024*(-2*m32
         *(mQ32 + mU32) + 2*pow2(m32) + pow2(mQ32) + pow2(mU32))*pow7(m3))/(3.*
         (mQ32 - mU32)*pow3(m32 - mQ32)*pow3(m32 - mU32)) - (64*(22*pow11(m3) -
         30*msq2*mU32*pow2(mQ32)*pow3(m3) - 30*mQ32*msq2*pow2(mU32)*pow3(m3) +
         248*pow2(mQ32)*pow2(mU32)*pow3(m3) - 3*mU32*pow3(m3)*pow3(mQ32) - 3*
         mQ32*pow3(m3)*pow3(mU32) + 120*mQ32*msq2*mU32*pow5(m3) + 30*msq2*pow2(
         mQ32)*pow5(m3) - 393*mU32*pow2(mQ32)*pow5(m3) - 393*mQ32*pow2(mU32)*
         pow5(m3) + 30*msq2*pow2(mU32)*pow5(m3) + 3*pow3(mQ32)*pow5(m3) + 3*
         pow3(mU32)*pow5(m3) - 90*mQ32*msq2*pow7(m3) + 598*mQ32*mU32*pow7(m3) -
         90*msq2*mU32*pow7(m3) + 129*pow2(mQ32)*pow7(m3) + 129*pow2(mU32)*pow7
         (m3) - 170*mQ32*pow9(m3) + 60*msq2*pow9(m3) - 170*mU32*pow9(m3)))/(3.*
         (mQ32 - mU32)*pow3(m32 - mQ32)*pow3(m32 - mU32))) - (4*pow4(Xt)*(-120*
         msq2*pow2(mU32)*pow3(m32)*pow3(mQ32) - 120*msq2*pow2(mQ32)*pow3(m32)*
         pow3(mU32) + 240*msq2*pow2(m32)*pow3(mQ32)*pow3(mU32) - 1048*pow3(m32)
         *pow3(mQ32)*pow3(mU32) - 480*msq2*pow2(mQ32)*pow2(mU32)*pow4(m32) -
         600*msq2*mU32*pow3(mQ32)*pow4(m32) + 4492*pow2(mU32)*pow3(mQ32)*pow4(
         m32) - 600*mQ32*msq2*pow3(mU32)*pow4(m32) + 4492*pow2(mQ32)*pow3(mU32)
         *pow4(m32) + 60*msq2*pow2(m32)*pow2(mU32)*pow4(mQ32) + 330*msq2*mU32*
         pow3(m32)*pow4(mQ32) - 3327*pow2(mU32)*pow3(m32)*pow4(mQ32) - 30*m32*
         msq2*pow3(mU32)*pow4(mQ32) - 332*pow2(m32)*pow3(mU32)*pow4(mQ32) - 360
         *msq2*pow4(m32)*pow4(mQ32) + 4598*mU32*pow4(m32)*pow4(mQ32) + 60*msq2*
         pow2(m32)*pow2(mQ32)*pow4(mU32) + 330*mQ32*msq2*pow3(m32)*pow4(mU32) -
         3327*pow2(mQ32)*pow3(m32)*pow4(mU32) - 30*m32*msq2*pow3(mQ32)*pow4(
         mU32) - 332*pow2(m32)*pow3(mQ32)*pow4(mU32) + 4598*mQ32*pow4(m32)*pow4
         (mU32) - 360*msq2*pow4(m32)*pow4(mU32) + 522*m32*pow4(mQ32)*pow4(mU32)
         + 990*msq2*mU32*pow2(mQ32)*pow5(m32) + 990*mQ32*msq2*pow2(mU32)*pow5(
         m32) + 4118*pow2(mQ32)*pow2(mU32)*pow5(m32) + 510*msq2*pow3(mQ32)*pow5
         (m32) - 4986*mU32*pow3(mQ32)*pow5(m32) - 4986*mQ32*pow3(mU32)*pow5(m32
         ) + 510*msq2*pow3(mU32)*pow5(m32) - 1413*pow4(mQ32)*pow5(m32) - 1413*
         pow4(mU32)*pow5(m32) - 60*msq2*mU32*pow2(m32)*pow5(mQ32) - 30*m32*msq2
         *pow2(mU32)*pow5(mQ32) + 238*pow2(m32)*pow2(mU32)*pow5(mQ32) + 90*msq2
         *pow3(m32)*pow5(mQ32) - 490*mU32*pow3(m32)*pow5(mQ32) + 258*m32*pow3(
         mU32)*pow5(mQ32) + 74*pow4(m32)*pow5(mQ32) - 72*pow4(mU32)*pow5(mQ32)
         - 60*mQ32*msq2*pow2(m32)*pow5(mU32) - 30*m32*msq2*pow2(mQ32)*pow5(mU32
         ) + 238*pow2(m32)*pow2(mQ32)*pow5(mU32) - 490*mQ32*pow3(m32)*pow5(mU32
         ) + 90*msq2*pow3(m32)*pow5(mU32) + 258*m32*pow3(mQ32)*pow5(mU32) + 74*
         pow4(m32)*pow5(mU32) - 72*pow4(mQ32)*pow5(mU32) - 840*mQ32*msq2*mU32*
         pow6(m32) - 420*msq2*pow2(mQ32)*pow6(m32) - 8620*mU32*pow2(mQ32)*pow6(
         m32) - 8620*mQ32*pow2(mU32)*pow6(m32) - 420*msq2*pow2(mU32)*pow6(m32)
         + 1632*pow3(mQ32)*pow6(m32) + 1632*pow3(mU32)*pow6(m32) - 6*mU32*pow2(
         m32)*pow6(mQ32) - 3*m32*pow2(mU32)*pow6(mQ32) + 9*pow3(m32)*pow6(mQ32)
         - 6*mQ32*pow2(m32)*pow6(mU32) - 3*m32*pow2(mQ32)*pow6(mU32) + 9*pow3(
         m32)*pow6(mU32) + 180*mQ32*msq2*pow7(m32) + 14940*mQ32*mU32*pow7(m32)
         + 180*msq2*mU32*pow7(m32) + 3118*pow2(mQ32)*pow7(m32) + 3118*pow2(mU32
         )*pow7(m32) - 5460*mQ32*pow8(m32) - 5460*mU32*pow8(m32) + 2048*pow9(
         m32)))/(3.*pow3(mQ32 - mU32)*pow4(m32 - mQ32)*pow4(m32 - mU32))) + log
         (mQ32/MR2)*((1280*(mQ32 + mU32)*pow3(m32)*pow6(Xt))/(3.*pow2(m32 -
         mQ32)*pow2(m32 - mU32)*pow3(mQ32 - mU32)) - (8*pow2(Xt)*(-360*msq2*
         pow2(mQ32)*pow2(mU32)*pow3(m32) + 120*msq2*pow2(m32)*pow2(mU32)*pow3(
         mQ32) + 240*msq2*mU32*pow3(m32)*pow3(mQ32) - 2028*pow2(mU32)*pow3(m32)
         *pow3(mQ32) + 120*msq2*pow2(m32)*pow2(mQ32)*pow3(mU32) + 240*mQ32*msq2
         *pow3(m32)*pow3(mU32) - 1836*pow2(mQ32)*pow3(m32)*pow3(mU32) - 576*
         pow2(m32)*pow3(mQ32)*pow3(mU32) - 240*msq2*mU32*pow2(mQ32)*pow4(m32) -
         240*mQ32*msq2*pow2(mU32)*pow4(m32) + 7136*pow2(mQ32)*pow2(mU32)*pow4(
         m32) - 360*msq2*pow3(mQ32)*pow4(m32) + 2924*mU32*pow3(mQ32)*pow4(m32)
         + 2540*mQ32*pow3(mU32)*pow4(m32) - 360*msq2*pow3(mU32)*pow4(m32) - 60*
         msq2*mU32*pow2(m32)*pow4(mQ32) - 30*m32*msq2*pow2(mU32)*pow4(mQ32) +
         244*pow2(m32)*pow2(mU32)*pow4(mQ32) + 90*msq2*pow3(m32)*pow4(mQ32) -
         499*mU32*pow3(m32)*pow4(mQ32) + 261*m32*pow3(mU32)*pow4(mQ32) + 74*
         pow4(m32)*pow4(mQ32) - 60*mQ32*msq2*pow2(m32)*pow4(mU32) - 30*m32*msq2
         *pow2(mQ32)*pow4(mU32) + 244*pow2(m32)*pow2(mQ32)*pow4(mU32) - 499*
         mQ32*pow3(m32)*pow4(mU32) + 90*msq2*pow3(m32)*pow4(mU32) + 261*m32*
         pow3(mQ32)*pow4(mU32) + 74*pow4(m32)*pow4(mU32) - 72*pow4(mQ32)*pow4(
         mU32) + 480*mQ32*msq2*mU32*pow5(m32) + 510*msq2*pow2(mQ32)*pow5(m32) -
         7061*mU32*pow2(mQ32)*pow5(m32) - 6485*mQ32*pow2(mU32)*pow5(m32) + 510
         *msq2*pow2(mU32)*pow5(m32) - 613*pow3(mQ32)*pow5(m32) - 421*pow3(mU32)
         *pow5(m32) - 6*mU32*pow2(m32)*pow5(mQ32) - 3*m32*pow2(mU32)*pow5(mQ32)
         + 9*pow3(m32)*pow5(mQ32) - 6*mQ32*pow2(m32)*pow5(mU32) - 3*m32*pow2(
         mQ32)*pow5(mU32) + 9*pow3(m32)*pow5(mU32) - 420*mQ32*msq2*pow6(m32) +
         5620*mQ32*mU32*pow6(m32) - 420*msq2*mU32*pow6(m32) + 1568*pow2(mQ32)*
         pow6(m32) + 1184*pow2(mU32)*pow6(m32) - 1202*mQ32*pow7(m32) + 180*msq2
         *pow7(m32) - 1010*mU32*pow7(m32) + 172*pow8(m32)))/(3.*(mQ32 - mU32)*
         pow4(m32 - mQ32)*pow4(m32 - mU32)) - (4*(-360*msq2*pow2(mQ32)*pow2(
         mU32)*pow3(m32) + 120*msq2*pow2(m32)*pow2(mU32)*pow3(mQ32) + 240*msq2*
         mU32*pow3(m32)*pow3(mQ32) - 2188*pow2(mU32)*pow3(m32)*pow3(mQ32) + 120
         *msq2*pow2(m32)*pow2(mQ32)*pow3(mU32) + 240*mQ32*msq2*pow3(m32)*pow3(
         mU32) - 1996*pow2(mQ32)*pow3(m32)*pow3(mU32) - 576*pow2(m32)*pow3(mQ32
         )*pow3(mU32) - 240*msq2*mU32*pow2(mQ32)*pow4(m32) - 240*mQ32*msq2*pow2
         (mU32)*pow4(m32) + 8320*pow2(mQ32)*pow2(mU32)*pow4(m32) - 360*msq2*
         pow3(mQ32)*pow4(m32) + 3244*mU32*pow3(mQ32)*pow4(m32) + 2860*mQ32*pow3
         (mU32)*pow4(m32) - 360*msq2*pow3(mU32)*pow4(m32) - 60*msq2*mU32*pow2(
         m32)*pow4(mQ32) - 30*m32*msq2*pow2(mU32)*pow4(mQ32) + 244*pow2(m32)*
         pow2(mU32)*pow4(mQ32) + 90*msq2*pow3(m32)*pow4(mQ32) - 499*mU32*pow3(
         m32)*pow4(mQ32) + 261*m32*pow3(mU32)*pow4(mQ32) + 74*pow4(m32)*pow4(
         mQ32) - 60*mQ32*msq2*pow2(m32)*pow4(mU32) - 30*m32*msq2*pow2(mQ32)*
         pow4(mU32) + 196*pow2(m32)*pow2(mQ32)*pow4(mU32) - 403*mQ32*pow3(m32)*
         pow4(mU32) + 90*msq2*pow3(m32)*pow4(mU32) + 261*m32*pow3(mQ32)*pow4(
         mU32) - 38*pow4(m32)*pow4(mU32) - 72*pow4(mQ32)*pow4(mU32) + 480*mQ32*
         msq2*mU32*pow5(m32) + 510*msq2*pow2(mQ32)*pow5(m32) - 8757*mU32*pow2(
         mQ32)*pow5(m32) - 8373*mQ32*pow2(mU32)*pow5(m32) + 510*msq2*pow2(mU32)
         *pow5(m32) - 773*pow3(mQ32)*pow5(m32) - 325*pow3(mU32)*pow5(m32) - 6*
         mU32*pow2(m32)*pow5(mQ32) - 3*m32*pow2(mU32)*pow5(mQ32) + 9*pow3(m32)*
         pow5(mQ32) - 6*mQ32*pow2(m32)*pow5(mU32) - 3*m32*pow2(mQ32)*pow5(mU32)
         + 9*pow3(m32)*pow5(mU32) - 420*mQ32*msq2*pow6(m32) + 8052*mQ32*mU32*
         pow6(m32) - 420*msq2*mU32*pow6(m32) + 2288*pow2(mQ32)*pow6(m32) + 1664
         *pow2(mU32)*pow6(m32) - 2162*mQ32*pow7(m32) + 180*msq2*pow7(m32) -
         1810*mU32*pow7(m32) + 508*pow8(m32)))/(3.*pow4(m32 - mQ32)*pow4(m32 -
         mU32)) - (32*pow5(Xt)*(98*mQ32*pow11(m3) + 98*mU32*pow11(m3) - 60*msq2
         *pow2(mQ32)*pow2(mU32)*pow3(m3) - 30*msq2*mU32*pow3(m3)*pow3(mQ32) +
         257*pow2(mU32)*pow3(m3)*pow3(mQ32) - 30*mQ32*msq2*pow3(m3)*pow3(mU32)
         + 257*pow2(mQ32)*pow3(m3)*pow3(mU32) - 3*mU32*pow3(m3)*pow4(mQ32) - 3*
         mQ32*pow3(m3)*pow4(mU32) + 150*msq2*mU32*pow2(mQ32)*pow5(m3) + 150*
         mQ32*msq2*pow2(mU32)*pow5(m3) - 834*pow2(mQ32)*pow2(mU32)*pow5(m3) +
         30*msq2*pow3(mQ32)*pow5(m3) - 414*mU32*pow3(mQ32)*pow5(m3) - 414*mQ32*
         pow3(mU32)*pow5(m3) + 30*msq2*pow3(mU32)*pow5(m3) + 3*pow4(mQ32)*pow5(
         m3) + 3*pow4(mU32)*pow5(m3) - 180*mQ32*msq2*mU32*pow7(m3) - 90*msq2*
         pow2(mQ32)*pow7(m3) + 819*mU32*pow2(mQ32)*pow7(m3) + 819*mQ32*pow2(
         mU32)*pow7(m3) - 90*msq2*pow2(mU32)*pow7(m3) + 173*pow3(mQ32)*pow7(m3)
         + 173*pow3(mU32)*pow7(m3) + 60*mQ32*msq2*pow9(m3) - 516*mQ32*mU32*
         pow9(m3) + 60*msq2*mU32*pow9(m3) - 258*pow2(mQ32)*pow9(m3) - 258*pow2(
         mU32)*pow9(m3)))/(3.*pow3(m32 - mQ32)*pow3(m32 - mU32)*pow3(mQ32 -
         mU32)) - (32*Xt*(-178*mQ32*pow11(m3) - 134*mU32*pow11(m3) + 88*pow13(
         m3) + 30*msq2*mU32*pow3(m3)*pow3(mQ32) - 263*pow2(mU32)*pow3(m3)*pow3(
         mQ32) - 30*mQ32*msq2*pow3(m3)*pow3(mU32) + 239*pow2(mQ32)*pow3(m3)*
         pow3(mU32) + 3*mU32*pow3(m3)*pow4(mQ32) - 3*mQ32*pow3(m3)*pow4(mU32) -
         90*msq2*mU32*pow2(mQ32)*pow5(m3) + 90*mQ32*msq2*pow2(mU32)*pow5(m3) +
         72*pow2(mQ32)*pow2(mU32)*pow5(m3) - 30*msq2*pow3(mQ32)*pow5(m3) + 420
         *mU32*pow3(mQ32)*pow5(m3) - 372*mQ32*pow3(mU32)*pow5(m3) + 30*msq2*
         pow3(mU32)*pow5(m3) - 3*pow4(mQ32)*pow5(m3) + 3*pow4(mU32)*pow5(m3) +
         90*msq2*pow2(mQ32)*pow7(m3) - 577*mU32*pow2(mQ32)*pow7(m3) + 361*mQ32*
         pow2(mU32)*pow7(m3) - 90*msq2*pow2(mU32)*pow7(m3) - 173*pow3(mQ32)*
         pow7(m3) + 85*pow3(mU32)*pow7(m3) - 60*mQ32*msq2*pow9(m3) + 144*mQ32*
         mU32*pow9(m3) + 60*msq2*mU32*pow9(m3) + 314*pow2(mQ32)*pow9(m3) - 26*
         pow2(mU32)*pow9(m3)))/(3.*(mQ32 - mU32)*pow3(m32 - mQ32)*pow3(m32 -
         mU32)) + pow3(Xt)*((-128*(4*mQ32*mU32 - 2*m32*(mQ32 + mU32) + 2*pow2(
         m32) - pow2(mQ32) - pow2(mU32))*pow3(m3)*(2 + (8*(pow2(m32)*pow2(mQ32)
         + pow2(m32)*pow2(mU32) - 2*mQ32*pow3(m32) - 2*mU32*pow3(m32) + 2*pow4
         (m32)))/(3.*pow2(m32 - mQ32)*pow2(m32 - mU32))))/((m32 - mQ32)*(m32 -
         mU32)*pow3(mQ32 - mU32)) + (1024*(-2*m32*(mQ32 + mU32) + 2*pow2(m32) +
         pow2(mQ32) + pow2(mU32))*pow7(m3))/(3.*(mQ32 - mU32)*pow3(m32 - mQ32)
         *pow3(m32 - mU32)) + (64*(22*pow11(m3) - 30*msq2*mU32*pow2(mQ32)*pow3(
         m3) - 30*mQ32*msq2*pow2(mU32)*pow3(m3) + 248*pow2(mQ32)*pow2(mU32)*
         pow3(m3) - 3*mU32*pow3(m3)*pow3(mQ32) - 3*mQ32*pow3(m3)*pow3(mU32) +
         120*mQ32*msq2*mU32*pow5(m3) + 30*msq2*pow2(mQ32)*pow5(m3) - 393*mU32*
         pow2(mQ32)*pow5(m3) - 393*mQ32*pow2(mU32)*pow5(m3) + 30*msq2*pow2(mU32
         )*pow5(m3) + 3*pow3(mQ32)*pow5(m3) + 3*pow3(mU32)*pow5(m3) - 90*mQ32*
         msq2*pow7(m3) + 598*mQ32*mU32*pow7(m3) - 90*msq2*mU32*pow7(m3) + 129*
         pow2(mQ32)*pow7(m3) + 129*pow2(mU32)*pow7(m3) - 170*mQ32*pow9(m3) + 60
         *msq2*pow9(m3) - 170*mU32*pow9(m3)))/(3.*(mQ32 - mU32)*pow3(m32 - mQ32
         )*pow3(m32 - mU32))) + (4*pow4(Xt)*(-120*msq2*pow2(mU32)*pow3(m32)*
         pow3(mQ32) - 120*msq2*pow2(mQ32)*pow3(m32)*pow3(mU32) + 240*msq2*pow2(
         m32)*pow3(mQ32)*pow3(mU32) - 1048*pow3(m32)*pow3(mQ32)*pow3(mU32) -
         480*msq2*pow2(mQ32)*pow2(mU32)*pow4(m32) - 600*msq2*mU32*pow3(mQ32)*
         pow4(m32) + 4492*pow2(mU32)*pow3(mQ32)*pow4(m32) - 600*mQ32*msq2*pow3(
         mU32)*pow4(m32) + 4492*pow2(mQ32)*pow3(mU32)*pow4(m32) + 60*msq2*pow2(
         m32)*pow2(mU32)*pow4(mQ32) + 330*msq2*mU32*pow3(m32)*pow4(mQ32) - 3327
         *pow2(mU32)*pow3(m32)*pow4(mQ32) - 30*m32*msq2*pow3(mU32)*pow4(mQ32) -
         332*pow2(m32)*pow3(mU32)*pow4(mQ32) - 360*msq2*pow4(m32)*pow4(mQ32) +
         4598*mU32*pow4(m32)*pow4(mQ32) + 60*msq2*pow2(m32)*pow2(mQ32)*pow4(
         mU32) + 330*mQ32*msq2*pow3(m32)*pow4(mU32) - 3327*pow2(mQ32)*pow3(m32)
         *pow4(mU32) - 30*m32*msq2*pow3(mQ32)*pow4(mU32) - 332*pow2(m32)*pow3(
         mQ32)*pow4(mU32) + 4598*mQ32*pow4(m32)*pow4(mU32) - 360*msq2*pow4(m32)
         *pow4(mU32) + 522*m32*pow4(mQ32)*pow4(mU32) + 990*msq2*mU32*pow2(mQ32)
         *pow5(m32) + 990*mQ32*msq2*pow2(mU32)*pow5(m32) + 4118*pow2(mQ32)*pow2
         (mU32)*pow5(m32) + 510*msq2*pow3(mQ32)*pow5(m32) - 4986*mU32*pow3(mQ32
         )*pow5(m32) - 4986*mQ32*pow3(mU32)*pow5(m32) + 510*msq2*pow3(mU32)*
         pow5(m32) - 1413*pow4(mQ32)*pow5(m32) - 1413*pow4(mU32)*pow5(m32) - 60
         *msq2*mU32*pow2(m32)*pow5(mQ32) - 30*m32*msq2*pow2(mU32)*pow5(mQ32) +
         238*pow2(m32)*pow2(mU32)*pow5(mQ32) + 90*msq2*pow3(m32)*pow5(mQ32) -
         490*mU32*pow3(m32)*pow5(mQ32) + 258*m32*pow3(mU32)*pow5(mQ32) + 74*
         pow4(m32)*pow5(mQ32) - 72*pow4(mU32)*pow5(mQ32) - 60*mQ32*msq2*pow2(
         m32)*pow5(mU32) - 30*m32*msq2*pow2(mQ32)*pow5(mU32) + 238*pow2(m32)*
         pow2(mQ32)*pow5(mU32) - 490*mQ32*pow3(m32)*pow5(mU32) + 90*msq2*pow3(
         m32)*pow5(mU32) + 258*m32*pow3(mQ32)*pow5(mU32) + 74*pow4(m32)*pow5(
         mU32) - 72*pow4(mQ32)*pow5(mU32) - 840*mQ32*msq2*mU32*pow6(m32) - 420*
         msq2*pow2(mQ32)*pow6(m32) - 8620*mU32*pow2(mQ32)*pow6(m32) - 8620*mQ32
         *pow2(mU32)*pow6(m32) - 420*msq2*pow2(mU32)*pow6(m32) + 1632*pow3(mQ32
         )*pow6(m32) + 1632*pow3(mU32)*pow6(m32) - 6*mU32*pow2(m32)*pow6(mQ32)
         - 3*m32*pow2(mU32)*pow6(mQ32) + 9*pow3(m32)*pow6(mQ32) - 6*mQ32*pow2(
         m32)*pow6(mU32) - 3*m32*pow2(mQ32)*pow6(mU32) + 9*pow3(m32)*pow6(mU32)
         + 180*mQ32*msq2*pow7(m32) + 14940*mQ32*mU32*pow7(m32) + 180*msq2*mU32
         *pow7(m32) + 3118*pow2(mQ32)*pow7(m32) + 3118*pow2(mU32)*pow7(m32) -
         5460*mQ32*pow8(m32) - 5460*mU32*pow8(m32) + 2048*pow9(m32)))/(3.*pow3(
         mQ32 - mU32)*pow4(m32 - mQ32)*pow4(m32 - mU32)))) + log(m32/MR2)*((-16
         *(3*m32*mQ32 + 3*m32*mU32 + 5*mQ32*mU32 - 11*pow2(m32))*pow2(m32)*(2*
         m32*(mQ32 + mU32)*pow2(mQ32 - mU32) + pow2(m32)*(2*mQ32*mU32 - 4*pow2(
         mQ32) - 4*pow2(mU32)) + mQ32*mU32*(pow2(mQ32) + pow2(mU32)) + 2*(mQ32
         + mU32)*pow3(m32)))/(3.*mQ32*mU32*pow3(m32 - mQ32)*pow3(m32 - mU32)) +
         (128*(pow2(m32)*pow2(mQ32) + pow2(m32)*pow2(mU32) - 2*mQ32*pow3(m32)
         - 2*mU32*pow3(m32) + 2*pow4(m32)))/(3.*pow2(m32 - mQ32)*pow2(m32 -
         mU32)) + (96*(-(mQ32*mU32) + pow2(m32))*(pow2(m32)*pow2(mQ32) + pow2(
         m32)*pow2(mU32) - 2*mQ32*pow3(m32) - 2*mU32*pow3(m32) + 2*pow4(m32)))/
         (pow3(m32 - mQ32)*pow3(m32 - mU32)) + (16*(2*m32*mU32*pow2(mQ32) +
         pow2(m32)*pow2(mQ32) + 2*m32*mQ32*pow2(mU32) + pow2(m32)*pow2(mU32) -
         2*pow2(mQ32)*pow2(mU32) - 4*mQ32*pow3(m32) - 4*mU32*pow3(m32) + 4*pow4
         (m32)))/(pow2(m32 - mQ32)*pow2(m32 - mU32)) + pow2(Xt)*((-64*(m32 -
         mQ32 - mU32)*pow2(m32)*(-3*m32*mQ32 - 3*m32*mU32 - 5*mQ32*mU32 + 11*
         pow2(m32)))/(3.*mQ32*mU32*pow2(m32 - mQ32)*pow2(m32 - mU32)) - (8192*(
         -(mU32*dilog(1 - m32/mQ32)) + m32*(dilog(1 - m32/mQ32) - dilog(1 -
         m32/mU32)) + mQ32*dilog(1 - m32/mU32))*pow3(m32))/(3.*(mQ32 - mU32)*
         pow2(m32 - mQ32)*pow2(m32 - mU32)) - (64*m32*(2 + (8*(pow2(m32)*pow2(
         mQ32) + pow2(m32)*pow2(mU32) - 2*mQ32*pow3(m32) - 2*mU32*pow3(m32) + 2
         *pow4(m32)))/(3.*pow2(m32 - mQ32)*pow2(m32 - mU32))))/(mQ32*mU32)) +
         pow3(Xt)*((128*m3*(2*m32*(dilog(1 - m32/mQ32) - dilog(1 - m32/mU32)) +
         mQ32*(-4 - dilog(1 - m32/mQ32) + dilog(1 - m32/mU32)) + mU32*(4 -
         dilog(1 - m32/mQ32) + dilog(1 - m32/mU32)))*(2 + (8*(pow2(m32)*pow2(
         mQ32) + pow2(m32)*pow2(mU32) - 2*mQ32*pow3(m32) - 2*mU32*pow3(m32) + 2
         *pow4(m32)))/(3.*pow2(m32 - mQ32)*pow2(m32 - mU32))))/pow3(mQ32 - mU32
         ) + (2048*pow5(m3))/(3.*(m32 - mQ32)*mQ32*(m32 - mU32)*mU32)) + (16*(2
         + (8*(pow2(m32)*pow2(mQ32) + pow2(m32)*pow2(mU32) - 2*mQ32*pow3(m32)
         - 2*mU32*pow3(m32) + 2*pow4(m32)))/(3.*pow2(m32 - mQ32)*pow2(m32 -
         mU32)))*(9*m32*(mQ32 + mU32)*pow2(mQ32)*pow2(mU32) - 3*pow3(mQ32)*pow3
         (mU32) + pow3(m32)*(mU32*(19 + 8*dilog(1 - m32/mU32))*pow2(mQ32) +
         mQ32*(19 + 8*dilog(1 - m32/mQ32))*pow2(mU32) + 2*pow3(mQ32) + 2*pow3(
         mU32)) - 4*pow2(m32)*(6*pow2(mQ32)*pow2(mU32) + mU32*(2 + dilog(1 -
         m32/mU32))*pow3(mQ32) + mQ32*(2 + dilog(1 - m32/mQ32))*pow3(mU32)) - (
         mQ32*mU32*(13 + 4*dilog(1 - m32/mQ32) + 4*dilog(1 - m32/mU32)) + 4*
         pow2(mQ32) + 4*pow2(mU32))*pow4(m32) + 2*(mQ32 + mU32)*pow5(m32)))/(
         mQ32*mU32*pow2(m32 - mQ32)*pow2(m32 - mU32)) - (16*(90*msq2*mU32*pow2(
         m32)*pow2(mQ32) + 90*mQ32*msq2*pow2(m32)*pow2(mU32) + 1358*pow2(m32)*
         pow2(mQ32)*pow2(mU32) - 180*mQ32*msq2*mU32*pow3(m32) + 270*msq2*pow2(
         mQ32)*pow3(m32) - 2829*mU32*pow2(mQ32)*pow3(m32) - 2829*mQ32*pow2(mU32
         )*pow3(m32) + 270*msq2*pow2(mU32)*pow3(m32) - 30*m32*msq2*mU32*pow3(
         mQ32) - 90*msq2*pow2(m32)*pow3(mQ32) + 610*mU32*pow2(m32)*pow3(mQ32) -
         213*m32*pow2(mU32)*pow3(mQ32) - 135*pow3(m32)*pow3(mQ32) - 30*m32*
         mQ32*msq2*pow3(mU32) + 610*mQ32*pow2(m32)*pow3(mU32) - 90*msq2*pow2(
         m32)*pow3(mU32) - 213*m32*pow2(mQ32)*pow3(mU32) - 135*pow3(m32)*pow3(
         mU32) + 24*pow3(mQ32)*pow3(mU32) - 240*mQ32*msq2*pow4(m32) + 4588*mQ32
         *mU32*pow4(m32) - 240*msq2*mU32*pow4(m32) + 898*pow2(mQ32)*pow4(m32) +
         898*pow2(mU32)*pow4(m32) - 3*m32*mU32*pow4(mQ32) - 9*pow2(m32)*pow4(
         mQ32) - 3*m32*mQ32*pow4(mU32) - 9*pow2(m32)*pow4(mU32) - 1604*mQ32*
         pow5(m32) + 180*msq2*pow5(m32) - 1604*mU32*pow5(m32) + 600*pow6(m32)))
         /(3.*pow3(m32 - mQ32)*pow3(m32 - mU32)) + pow4(Xt)*((-4096*(2*m32*(
         dilog(1 - m32/mQ32) - dilog(1 - m32/mU32)) + mQ32*(-4 - dilog(1 -
         m32/mQ32) + dilog(1 - m32/mU32)) + mU32*(4 - dilog(1 - m32/mQ32) +
         dilog(1 - m32/mU32)))*pow2(m32))/(3.*(m32 - mQ32)*(m32 - mU32)*pow3(
         mQ32 - mU32)) - (32*m32*(3*m32*mQ32 + 3*m32*mU32 + 5*mQ32*mU32 - 11*
         pow2(m32))*((mQ32 + mU32)*pow2(mQ32)*pow2(mU32) - m32*mQ32*mU32*(4*
         mQ32*mU32 + pow2(mQ32) + pow2(mU32)) - 2*(mQ32*mU32 + pow2(mQ32) +
         pow2(mU32))*pow3(m32) + pow2(m32)*pow3(mQ32 + mU32) + (mQ32 + mU32)*
         pow4(m32)))/(3.*mQ32*mU32*pow2(mQ32 - mU32)*pow3(m32 - mQ32)*pow3(m32
         - mU32)) + (32*((2*mQ32*mU32*(-dilog(1 - m32/mQ32) + dilog(1 -
         m32/mU32)) - pow2(mQ32) + pow2(mU32))*pow3(m32) + (-5 + dilog(1 -
         m32/mQ32) - dilog(1 - m32/mU32))*pow2(mU32)*pow3(mQ32) + pow2(m32)*(3*
         mU32*(-2 + dilog(1 - m32/mQ32) - dilog(1 - m32/mU32))*pow2(mQ32) + 3*
         mQ32*(2 + dilog(1 - m32/mQ32) - dilog(1 - m32/mU32))*pow2(mU32) + pow3
         (mQ32) - pow3(mU32)) + (5 + dilog(1 - m32/mQ32) - dilog(1 - m32/mU32))
         *pow2(mQ32)*pow3(mU32) + m32*(4*(-dilog(1 - m32/mQ32) + dilog(1 -
         m32/mU32))*pow2(mQ32)*pow2(mU32) + mU32*(5 - dilog(1 - m32/mQ32) +
         dilog(1 - m32/mU32))*pow3(mQ32) + mQ32*(-5 - dilog(1 - m32/mQ32) +
         dilog(1 - m32/mU32))*pow3(mU32)))*(2 + (8*(pow2(m32)*pow2(mQ32) + pow2
         (m32)*pow2(mU32) - 2*mQ32*pow3(m32) - 2*mU32*pow3(m32) + 2*pow4(m32)))
         /(3.*pow2(m32 - mQ32)*pow2(m32 - mU32))))/(mQ32*(-m32 + mQ32)*(m32 -
         mU32)*mU32*pow3(mQ32 - mU32)) + (32*(90*msq2*mU32*pow2(m32)*pow2(mQ32)
         + 90*mQ32*msq2*pow2(m32)*pow2(mU32) + 1376*pow2(m32)*pow2(mQ32)*pow2(
         mU32) - 180*mQ32*msq2*mU32*pow3(m32) + 270*msq2*pow2(mQ32)*pow3(m32) -
         2840*mU32*pow2(mQ32)*pow3(m32) - 2840*mQ32*pow2(mU32)*pow3(m32) + 270
         *msq2*pow2(mU32)*pow3(m32) - 30*m32*msq2*mU32*pow3(mQ32) - 90*msq2*
         pow2(m32)*pow3(mQ32) + 623*mU32*pow2(m32)*pow3(mQ32) - 225*m32*pow2(
         mU32)*pow3(mQ32) - 132*pow3(m32)*pow3(mQ32) - 30*m32*mQ32*msq2*pow3(
         mU32) + 623*mQ32*pow2(m32)*pow3(mU32) - 90*msq2*pow2(m32)*pow3(mU32) -
         225*m32*pow2(mQ32)*pow3(mU32) - 132*pow3(m32)*pow3(mU32) + 30*pow3(
         mQ32)*pow3(mU32) - 240*mQ32*msq2*pow4(m32) + 4572*mQ32*mU32*pow4(m32)
         - 240*msq2*mU32*pow4(m32) + 873*pow2(mQ32)*pow4(m32) + 873*pow2(mU32)*
         pow4(m32) - 3*m32*mU32*pow4(mQ32) - 9*pow2(m32)*pow4(mQ32) - 3*m32*
         mQ32*pow4(mU32) - 9*pow2(m32)*pow4(mU32) - 1560*mQ32*pow5(m32) + 180*
         msq2*pow5(m32) - 1560*mU32*pow5(m32) + 568*pow6(m32)))/(3.*pow2(mQ32 -
         mU32)*pow3(m32 - mQ32)*pow3(m32 - mU32))) + log(msq2/MR2)*((96*(-m32
         + mQ32 + mU32)*pow2(m32)*pow2(Xt))/(mQ32*(-m32 + mQ32)*(m32 - mU32)*
         mU32) - (24*pow2(m32)*(2*m32*(mQ32 + mU32)*pow2(mQ32 - mU32) + pow2(
         m32)*(2*mQ32*mU32 - 4*pow2(mQ32) - 4*pow2(mU32)) + mQ32*mU32*(pow2(
         mQ32) + pow2(mU32)) + 2*(mQ32 + mU32)*pow3(m32)))/(mQ32*mU32*pow2(m32
         - mQ32)*pow2(m32 - mU32)) - (24*(pow2(m32)*pow2(mQ32) + pow2(m32)*pow2
         (mU32) - 2*mQ32*pow3(m32) - 2*mU32*pow3(m32) + 2*pow4(m32)))/(pow2(m32
         - mQ32)*pow2(m32 - mU32)) + pow4(Xt)*((-48*m32*((mQ32 + mU32)*pow2(
         mQ32)*pow2(mU32) - m32*mQ32*mU32*(4*mQ32*mU32 + pow2(mQ32) + pow2(mU32
         )) - 2*(mQ32*mU32 + pow2(mQ32) + pow2(mU32))*pow3(m32) + pow2(m32)*
         pow3(mQ32 + mU32) + (mQ32 + mU32)*pow4(m32)))/(mQ32*mU32*pow2(m32 -
         mQ32)*pow2(m32 - mU32)*pow2(mQ32 - mU32)) + (12*((4*(pow2(m32)*pow2(
         mQ32) + pow2(m32)*pow2(mU32) - 2*mQ32*pow3(m32) - 2*mU32*pow3(m32) + 2
         *pow4(m32)))/(pow2(m32 - mQ32)*pow2(m32 - mU32)) - 4*((10*(pow2(m32)*
         pow2(mQ32) + pow2(m32)*pow2(mU32) - 2*mQ32*pow3(m32) - 2*mU32*pow3(m32
         ) + 2*pow4(m32)))/(9.*pow2(-m32 + mQ32)*pow2(m32 - mU32)) + (10*msq2*(
         -3*mU32*pow2(m32)*pow2(mQ32) - 3*mQ32*pow2(m32)*pow2(mU32) + 6*mQ32*
         mU32*pow3(m32) - 9*pow2(mQ32)*pow3(m32) - 9*pow2(mU32)*pow3(m32) + m32
         *mU32*pow3(mQ32) + 3*pow2(m32)*pow3(mQ32) + m32*mQ32*pow3(mU32) + 3*
         pow2(m32)*pow3(mU32) + 8*mQ32*pow4(m32) + 8*mU32*pow4(m32) - 6*pow5(
         m32)))/(3.*pow3(-m32 + mQ32)*pow3(m32 - mU32)))))/pow2(mQ32 - mU32)) +
         24*((10*(pow2(m32)*pow2(mQ32) + pow2(m32)*pow2(mU32) - 2*mQ32*pow3(
         m32) - 2*mU32*pow3(m32) + 2*pow4(m32)))/(9.*pow2(-m32 + mQ32)*pow2(m32
         - mU32)) + (10*msq2*(-3*mU32*pow2(m32)*pow2(mQ32) - 3*mQ32*pow2(m32)*
         pow2(mU32) + 6*mQ32*mU32*pow3(m32) - 9*pow2(mQ32)*pow3(m32) - 9*pow2(
         mU32)*pow3(m32) + m32*mU32*pow3(mQ32) + 3*pow2(m32)*pow3(mQ32) + m32*
         mQ32*pow3(mU32) + 3*pow2(m32)*pow3(mU32) + 8*mQ32*pow4(m32) + 8*mU32*
         pow4(m32) - 6*pow5(m32)))/(3.*pow3(-m32 + mQ32)*pow3(m32 - mU32))) + (
         640*pow5(Xt)*(-6*mQ32*msq2*pow3(m3) + mQ32*mU32*pow3(m3) - 6*msq2*mU32
         *pow3(m3) - mQ32*pow5(m3) + 12*msq2*pow5(m3) - mU32*pow5(m3) + pow7(m3
         )))/(3.*pow2(m32 - mQ32)*pow2(m32 - mU32)*pow2(mQ32 - mU32)) - (32*Xt*
         (-60*mQ32*msq2*pow3(m3) + mQ32*mU32*pow3(m3) - 60*msq2*mU32*pow3(m3) -
         mQ32*pow5(m3) + 120*msq2*pow5(m3) - mU32*pow5(m3) + pow7(m3)))/(3.*
         pow2(m32 - mQ32)*pow2(m32 - mU32))) + pow5(Xt)*((128*(3*m32*mQ32 + 3*
         m32*mU32 + 5*mQ32*mU32 - 11*pow2(m32))*pow3(m3))/(3.*pow2(m32 - mQ32)*
         pow2(m32 - mU32)*pow2(mQ32 - mU32)) - (1024*pow3(m3)*((2*mQ32*mU32*(
         -dilog(1 - m32/mQ32) + dilog(1 - m32/mU32)) - pow2(mQ32) + pow2(mU32))
         *pow3(m32) + (-5 + dilog(1 - m32/mQ32) - dilog(1 - m32/mU32))*pow2(
         mU32)*pow3(mQ32) + pow2(m32)*(3*mU32*(-2 + dilog(1 - m32/mQ32) - dilog
         (1 - m32/mU32))*pow2(mQ32) + 3*mQ32*(2 + dilog(1 - m32/mQ32) - dilog(1
         - m32/mU32))*pow2(mU32) + pow3(mQ32) - pow3(mU32)) + (5 + dilog(1 -
         m32/mQ32) - dilog(1 - m32/mU32))*pow2(mQ32)*pow3(mU32) + m32*(4*(
         -dilog(1 - m32/mQ32) + dilog(1 - m32/mU32))*pow2(mQ32)*pow2(mU32) +
         mU32*(5 - dilog(1 - m32/mQ32) + dilog(1 - m32/mU32))*pow3(mQ32) + mQ32
         *(-5 - dilog(1 - m32/mQ32) + dilog(1 - m32/mU32))*pow3(mU32))))/(3.*(
         m32 - mQ32)*mQ32*(-m32 + mQ32)*mU32*pow2(m32 - mU32)*pow3(mQ32 - mU32)
         ) - (256*(-10*mQ32*msq2*pow3(m3) + 54*mQ32*mU32*pow3(m3) - 10*msq2*
         mU32*pow3(m3) - pow2(mQ32)*pow3(m3) - pow2(mU32)*pow3(m3) - 37*mQ32*
         pow5(m3) + 20*msq2*pow5(m3) - 37*mU32*pow5(m3) + 22*pow7(m3)))/(pow2(
         m32 - mQ32)*pow2(m32 - mU32)*pow2(mQ32 - mU32))) + Xt*((-704*pow3(m3))
         /(3.*(m32 - mQ32)*(m32 - mU32)) - (384*(-(mQ32*mU32) + pow2(m32))*pow3
         (m3))/(pow2(m32 - mQ32)*pow2(m32 - mU32)) + (256*(-(mU32*dilog(1 -
         m32/mQ32)) + m32*(dilog(1 - m32/mQ32) - dilog(1 - m32/mU32)) + mQ32*
         dilog(1 - m32/mU32))*pow3(m3)*(2 + (8*(pow2(m32)*pow2(mQ32) + pow2(m32
         )*pow2(mU32) - 2*mQ32*pow3(m32) - 2*mU32*pow3(m32) + 2*pow4(m32)))/(3.
         *pow2(m32 - mQ32)*pow2(m32 - mU32))))/((m32 - mQ32)*(m32 - mU32)*(mQ32
         - mU32)) - (512*pow3(m3)*(9*m32*(mQ32 + mU32)*pow2(mQ32)*pow2(mU32) -
         3*pow3(mQ32)*pow3(mU32) + pow3(m32)*(mU32*(19 + 8*dilog(1 - m32/mU32)
         )*pow2(mQ32) + mQ32*(19 + 8*dilog(1 - m32/mQ32))*pow2(mU32) + 2*pow3(
         mQ32) + 2*pow3(mU32)) - 4*pow2(m32)*(6*pow2(mQ32)*pow2(mU32) + mU32*(2
         + dilog(1 - m32/mU32))*pow3(mQ32) + mQ32*(2 + dilog(1 - m32/mQ32))*
         pow3(mU32)) - (mQ32*mU32*(13 + 4*dilog(1 - m32/mQ32) + 4*dilog(1 -
         m32/mU32)) + 4*pow2(mQ32) + 4*pow2(mU32))*pow4(m32) + 2*(mQ32 + mU32)*
         pow5(m32)))/(3.*mQ32*mU32*pow3(m32 - mQ32)*pow3(m32 - mU32)) + (64*(
         -60*mQ32*msq2*pow3(m3) + 317*mQ32*mU32*pow3(m3) - 60*msq2*mU32*pow3(m3
         ) - 6*pow2(mQ32)*pow3(m3) - 6*pow2(mU32)*pow3(m3) - 225*mQ32*pow5(m3)
         + 120*msq2*pow5(m3) - 225*mU32*pow5(m3) + 145*pow7(m3)))/(3.*pow2(m32
         - mQ32)*pow2(m32 - mU32))) + pow2(log(mQ32/MR2))*((-2560*mQ32*(mQ32 +
         mU32)*pow2(m32)*pow6(Xt))/(3.*(m32 - mU32)*pow2(m32 - mQ32)*pow4(mQ32
         - mU32)) - (8*(189*pow2(mQ32)*pow2(mU32)*pow3(m32) - 369*pow2(m32)*
         pow2(mU32)*pow3(mQ32) + 531*mU32*pow3(m32)*pow3(mQ32) - 59*pow2(m32)*
         pow2(mQ32)*pow3(mU32) - 75*mQ32*pow3(m32)*pow3(mU32) + 119*m32*pow3(
         mQ32)*pow3(mU32) - 427*mU32*pow2(mQ32)*pow4(m32) + 221*mQ32*pow2(mU32)
         *pow4(m32) - 209*pow3(mQ32)*pow4(m32) - 17*pow3(mU32)*pow4(m32) - 125*
         mU32*pow2(m32)*pow4(mQ32) + 76*m32*pow2(mU32)*pow4(mQ32) + 25*pow3(m32
         )*pow4(mQ32) - 24*pow3(mU32)*pow4(mQ32) - 50*mQ32*mU32*pow5(m32) + 249
         *pow2(mQ32)*pow5(m32) + 51*pow2(mU32)*pow5(m32) + 3*m32*mU32*pow5(mQ32
         ) + 9*pow2(m32)*pow5(mQ32) - 84*mQ32*pow6(m32) - 100*mU32*pow6(m32) +
         66*pow7(m32)))/(3.*pow3(m32 - mU32)*pow4(m32 - mQ32)) - (16*pow2(Xt)*(
         1246*pow2(mU32)*pow3(m32)*pow3(mQ32) - 180*pow2(mQ32)*pow3(m32)*pow3(
         mU32) + 26*pow2(m32)*pow3(mQ32)*pow3(mU32) - 956*pow2(mQ32)*pow2(mU32)
         *pow4(m32) - 2702*mU32*pow3(mQ32)*pow4(m32) + 246*mQ32*pow3(mU32)*pow4
         (m32) - 692*pow2(m32)*pow2(mU32)*pow4(mQ32) + 1242*mU32*pow3(m32)*pow4
         (mQ32) + 139*m32*pow3(mU32)*pow4(mQ32) - 593*pow4(m32)*pow4(mQ32) + 71
         *pow2(m32)*pow2(mQ32)*pow4(mU32) + 51*mQ32*pow3(m32)*pow4(mU32) - 119*
         m32*pow3(mQ32)*pow4(mU32) - 35*pow4(m32)*pow4(mU32) + 24*pow4(mQ32)*
         pow4(mU32) + 2641*mU32*pow2(mQ32)*pow5(m32) - 215*mQ32*pow2(mU32)*pow5
         (m32) + 1429*pow3(mQ32)*pow5(m32) - 151*pow3(mU32)*pow5(m32) - 238*
         mU32*pow2(m32)*pow5(mQ32) + 145*m32*pow2(mU32)*pow5(mQ32) + 81*pow3(
         m32)*pow5(mQ32) - 48*pow3(mU32)*pow5(mQ32) - 676*mQ32*mU32*pow6(m32) -
         1516*pow2(mQ32)*pow6(m32) + 424*pow2(mU32)*pow6(m32) + 3*m32*mU32*
         pow6(mQ32) + 9*pow2(m32)*pow6(mQ32) + 582*mQ32*pow7(m32) - 238*mU32*
         pow7(m32)))/(3.*pow2(mQ32 - mU32)*pow3(m32 - mU32)*pow4(m32 - mQ32)) -
         (8*pow4(Xt)*(-4362*pow3(m32)*pow3(mQ32)*pow3(mU32) + 15882*pow2(mU32)
         *pow3(mQ32)*pow4(m32) + 26*pow2(mQ32)*pow3(mU32)*pow4(m32) - 7300*pow2
         (mU32)*pow3(m32)*pow4(mQ32) + 3118*pow2(m32)*pow3(mU32)*pow4(mQ32) +
         5951*mU32*pow4(m32)*pow4(mQ32) + 1513*pow2(mQ32)*pow3(m32)*pow4(mU32)
         - 673*pow2(m32)*pow3(mQ32)*pow4(mU32) - 747*mQ32*pow4(m32)*pow4(mU32)
         - 176*m32*pow4(mQ32)*pow4(mU32) - 10006*pow2(mQ32)*pow2(mU32)*pow5(m32
         ) - 15702*mU32*pow3(mQ32)*pow5(m32) + 2222*mQ32*pow3(mU32)*pow5(m32) -
         1581*pow4(mQ32)*pow5(m32) + 11*pow4(mU32)*pow5(m32) + 562*pow2(m32)*
         pow2(mU32)*pow5(mQ32) - 395*mU32*pow3(m32)*pow5(mQ32) - 380*m32*pow3(
         mU32)*pow5(mQ32) + 105*pow4(m32)*pow5(mQ32) + 48*pow4(mU32)*pow5(mQ32)
         + 17*pow2(m32)*pow2(mQ32)*pow5(mU32) - 131*mQ32*pow3(m32)*pow5(mU32)
         + 71*m32*pow3(mQ32)*pow5(mU32) + 63*pow4(m32)*pow5(mU32) - 12*pow4(
         mQ32)*pow5(mU32) + 13092*mU32*pow2(mQ32)*pow6(m32) + 44*mQ32*pow2(mU32
         )*pow6(m32) + 4844*pow3(mQ32)*pow6(m32) - 668*pow3(mU32)*pow6(m32) +
         281*mU32*pow2(m32)*pow6(mQ32) - 184*m32*pow2(mU32)*pow6(mQ32) - 109*
         pow3(m32)*pow6(mQ32) + 60*pow3(mU32)*pow6(mQ32) - 2824*mQ32*mU32*pow7(
         m32) - 4690*pow2(mQ32)*pow7(m32) + 1018*pow2(mU32)*pow7(m32) - 3*m32*
         mU32*pow7(mQ32) - 9*pow2(m32)*pow7(mQ32) + 1448*mQ32*pow8(m32) - 424*
         mU32*pow8(m32)))/(3.*pow3(m32 - mU32)*pow4(m32 - mQ32)*pow4(mQ32 -
         mU32)) + (32*Xt*(52*pow11(m3) - 45*pow2(mQ32)*pow2(mU32)*pow3(m3) - 12
         *m3*pow2(mU32)*pow3(mQ32) - 25*mU32*pow3(m3)*pow3(mQ32) + 6*pow3(m3)*
         pow4(mQ32) + 221*mU32*pow2(mQ32)*pow5(m3) + 84*mQ32*pow2(mU32)*pow5(m3
         ) + 3*pow3(mQ32)*pow5(m3) - 299*mQ32*mU32*pow7(m3) - 94*pow2(mQ32)*
         pow7(m3) + 69*pow2(mU32)*pow7(m3) + 129*mQ32*pow9(m3) - 89*mU32*pow9(
         m3)))/(3.*(mQ32 - mU32)*pow2(m32 - mU32)*pow3(m32 - mQ32)) + (32*pow5(
         Xt)*(64*mU32*pow11(m3) + 122*pow2(mU32)*pow3(m3)*pow3(mQ32) + 49*pow2(
         mQ32)*pow3(m3)*pow3(mU32) + 12*m3*pow3(mQ32)*pow3(mU32) + 12*m3*pow2(
         mU32)*pow4(mQ32) + 67*mU32*pow3(m3)*pow4(mQ32) - 401*pow2(mQ32)*pow2(
         mU32)*pow5(m3) - 328*mU32*pow3(mQ32)*pow5(m3) - 60*mQ32*pow3(mU32)*
         pow5(m3) - 51*pow4(mQ32)*pow5(m3) - 6*pow3(m3)*pow5(mQ32) + 413*mU32*
         pow2(mQ32)*pow7(m3) + 330*mQ32*pow2(mU32)*pow7(m3) + 178*pow3(mQ32)*
         pow7(m3) + 31*pow3(mU32)*pow7(m3) - 248*mQ32*mU32*pow9(m3) - 89*pow2(
         mQ32)*pow9(m3) - 95*pow2(mU32)*pow9(m3)))/(3.*pow2(m32 - mU32)*pow3(
         m32 - mQ32)*pow4(mQ32 - mU32)) + (64*pow3(Xt)*(-78*mQ32*pow11(m3) - 38
         *mU32*pow11(m3) + 40*pow13(m3) + 18*pow2(mU32)*pow3(m3)*pow3(mQ32) +
         15*pow2(mQ32)*pow3(m3)*pow3(mU32) + 18*m3*pow3(mQ32)*pow3(mU32) - 6*m3
         *pow2(mU32)*pow4(mQ32) - 75*mU32*pow3(m3)*pow4(mQ32) - 243*pow2(mQ32)*
         pow2(mU32)*pow5(m3) + 228*mU32*pow3(mQ32)*pow5(m3) - 2*mQ32*pow3(mU32)
         *pow5(m3) + 49*pow4(mQ32)*pow5(m3) + 6*pow3(m3)*pow5(mQ32) + 41*mU32*
         pow2(mQ32)*pow7(m3) + 130*mQ32*pow2(mU32)*pow7(m3) - 212*pow3(mQ32)*
         pow7(m3) + pow3(mU32)*pow7(m3) - 60*mQ32*mU32*pow9(m3) + 163*pow2(mQ32
         )*pow9(m3) + 5*pow2(mU32)*pow9(m3)))/(3.*pow2(m32 - mU32)*pow3(m32 -
         mQ32)*pow3(mQ32 - mU32))) + log(mU32/MR2)*((-5120*mU32*pow2(m32)*pow6(
         Xt))/(3.*(m32 - mQ32)*pow2(m32 - mU32)*pow3(mQ32 - mU32)) + log(
         msq2/MR2)*((-16*(45*msq2*mU32*pow2(m32)*pow2(mQ32) + 45*mQ32*msq2*pow2
         (m32)*pow2(mU32) - 90*mQ32*msq2*mU32*pow3(m32) + 135*msq2*pow2(mQ32)*
         pow3(m32) + 12*mU32*pow2(mQ32)*pow3(m32) - 15*mQ32*pow2(mU32)*pow3(m32
         ) + 135*msq2*pow2(mU32)*pow3(m32) - 15*m32*msq2*mU32*pow3(mQ32) - 45*
         msq2*pow2(m32)*pow3(mQ32) - 4*mU32*pow2(m32)*pow3(mQ32) + 4*pow3(m32)*
         pow3(mQ32) - 15*m32*mQ32*msq2*pow3(mU32) + 5*mQ32*pow2(m32)*pow3(mU32)
         - 45*msq2*pow2(m32)*pow3(mU32) - 5*pow3(m32)*pow3(mU32) - 120*mQ32*
         msq2*pow4(m32) + 3*mQ32*mU32*pow4(m32) - 120*msq2*mU32*pow4(m32) - 12*
         pow2(mQ32)*pow4(m32) + 15*pow2(mU32)*pow4(m32) + 7*mQ32*pow5(m32) + 90
         *msq2*pow5(m32) - 11*mU32*pow5(m32) + pow6(m32)))/(3.*pow3(m32 - mQ32)
         *pow3(m32 - mU32)) + (160*pow2(Xt)*(9*msq2*mU32*pow2(m32)*pow2(mQ32) +
         9*mQ32*msq2*pow2(m32)*pow2(mU32) - 18*mQ32*msq2*mU32*pow3(m32) + 27*
         msq2*pow2(mQ32)*pow3(m32) - 3*mU32*pow2(mQ32)*pow3(m32) - 3*mQ32*pow2(
         mU32)*pow3(m32) + 27*msq2*pow2(mU32)*pow3(m32) - 3*m32*msq2*mU32*pow3(
         mQ32) - 9*msq2*pow2(m32)*pow3(mQ32) + mU32*pow2(m32)*pow3(mQ32) - pow3
         (m32)*pow3(mQ32) - 3*m32*mQ32*msq2*pow3(mU32) + mQ32*pow2(m32)*pow3(
         mU32) - 9*msq2*pow2(m32)*pow3(mU32) - pow3(m32)*pow3(mU32) - 24*mQ32*
         msq2*pow4(m32) + 6*mQ32*mU32*pow4(m32) - 24*msq2*mU32*pow4(m32) + 3*
         pow2(mQ32)*pow4(m32) + 3*pow2(mU32)*pow4(m32) - 4*mQ32*pow5(m32) + 18*
         msq2*pow5(m32) - 4*mU32*pow5(m32) + 2*pow6(m32)))/(3.*(mQ32 - mU32)*
         pow3(m32 - mQ32)*pow3(m32 - mU32)) - (80*(mQ32 + mU32)*pow4(Xt)*(9*
         msq2*mU32*pow2(m32)*pow2(mQ32) + 9*mQ32*msq2*pow2(m32)*pow2(mU32) - 18
         *mQ32*msq2*mU32*pow3(m32) + 27*msq2*pow2(mQ32)*pow3(m32) - 3*mU32*pow2
         (mQ32)*pow3(m32) - 3*mQ32*pow2(mU32)*pow3(m32) + 27*msq2*pow2(mU32)*
         pow3(m32) - 3*m32*msq2*mU32*pow3(mQ32) - 9*msq2*pow2(m32)*pow3(mQ32) +
         mU32*pow2(m32)*pow3(mQ32) - pow3(m32)*pow3(mQ32) - 3*m32*mQ32*msq2*
         pow3(mU32) + mQ32*pow2(m32)*pow3(mU32) - 9*msq2*pow2(m32)*pow3(mU32) -
         pow3(m32)*pow3(mU32) - 24*mQ32*msq2*pow4(m32) + 6*mQ32*mU32*pow4(m32)
         - 24*msq2*mU32*pow4(m32) + 3*pow2(mQ32)*pow4(m32) + 3*pow2(mU32)*pow4
         (m32) - 4*mQ32*pow5(m32) + 18*msq2*pow5(m32) - 4*mU32*pow5(m32) + 2*
         pow6(m32)))/(3.*pow3(m32 - mQ32)*pow3(m32 - mU32)*pow3(mQ32 - mU32)) +
         (320*(mQ32 + mU32)*pow5(Xt)*(-6*mQ32*msq2*pow3(m3) + mQ32*mU32*pow3(
         m3) - 6*msq2*mU32*pow3(m3) - mQ32*pow5(m3) + 12*msq2*pow5(m3) - mU32*
         pow5(m3) + pow7(m3)))/(3.*pow2(m32 - mQ32)*pow2(m32 - mU32)*pow3(mQ32
         - mU32)) + Xt*((96*(2*m32 - mQ32 - mU32)*pow3(m3))/((m32 - mQ32)*(m32
         - mU32)*(mQ32 - mU32)) + (32*(-60*mQ32*msq2*pow3(m3) + mQ32*mU32*pow3(
         m3) - 60*msq2*mU32*pow3(m3) - mQ32*pow5(m3) + 120*msq2*pow5(m3) - mU32
         *pow5(m3) + pow7(m3)))/(3.*pow2(m32 - mQ32)*pow2(m32 - mU32))) + pow3(
         Xt)*((192*(4*mQ32*mU32 - 2*m32*(mQ32 + mU32) + 2*pow2(m32) - pow2(mQ32
         ) - pow2(mU32))*pow3(m3))/((m32 - mQ32)*(m32 - mU32)*pow3(mQ32 - mU32)
         ) - (64*(-60*mQ32*msq2*pow3(m3) + mQ32*mU32*pow3(m3) - 60*msq2*mU32*
         pow3(m3) - mQ32*pow5(m3) + 120*msq2*pow5(m3) - mU32*pow5(m3) + pow7(m3
         )))/(3.*(mQ32 - mU32)*pow2(m32 - mQ32)*pow2(m32 - mU32)))) + (32*pow2(
         Xt)*(180*msq2*pow2(mQ32)*pow2(mU32)*pow3(m32) - 90*msq2*pow2(m32)*pow2
         (mU32)*pow3(mQ32) - 270*msq2*mU32*pow3(m32)*pow3(mQ32) + 3396*pow2(
         mU32)*pow3(m32)*pow3(mQ32) - 90*msq2*pow2(m32)*pow2(mQ32)*pow3(mU32) -
         270*mQ32*msq2*pow3(m32)*pow3(mU32) + 3162*pow2(mQ32)*pow3(m32)*pow3(
         mU32) - 1626*pow2(m32)*pow3(mQ32)*pow3(mU32) + 240*msq2*mU32*pow2(mQ32
         )*pow4(m32) + 240*mQ32*msq2*pow2(mU32)*pow4(m32) - 5392*pow2(mQ32)*
         pow2(mU32)*pow4(m32) - 1386*mU32*pow3(mQ32)*pow4(m32) - 1014*mQ32*pow3
         (mU32)*pow4(m32) + 90*msq2*mU32*pow2(m32)*pow4(mQ32) + 30*m32*msq2*
         pow2(mU32)*pow4(mQ32) - 749*pow2(m32)*pow2(mU32)*pow4(mQ32) + 253*mU32
         *pow3(m32)*pow4(mQ32) + 297*m32*pow3(mU32)*pow4(mQ32) - pow4(m32)*pow4
         (mQ32) + 90*mQ32*msq2*pow2(m32)*pow4(mU32) + 30*m32*msq2*pow2(mQ32)*
         pow4(mU32) - 691*pow2(m32)*pow2(mQ32)*pow4(mU32) + 157*mQ32*pow3(m32)*
         pow4(mU32) + 285*m32*pow3(mQ32)*pow4(mU32) + 17*pow4(m32)*pow4(mU32) -
         48*pow4(mQ32)*pow4(mU32) - 180*mQ32*msq2*mU32*pow5(m32) + 2273*mU32*
         pow2(mQ32)*pow5(m32) + 1963*mQ32*pow2(mU32)*pow5(m32) + 3*pow3(mQ32)*
         pow5(m32) - 51*pow3(mU32)*pow5(m32) + 9*mU32*pow2(m32)*pow5(mQ32) + 3*
         m32*pow2(mU32)*pow5(mQ32) + 9*mQ32*pow2(m32)*pow5(mU32) + 3*m32*pow2(
         mQ32)*pow5(mU32) - 904*mQ32*mU32*pow6(m32) - 3*pow2(mQ32)*pow6(m32) +
         35*pow2(mU32)*pow6(m32) + mQ32*pow7(m32) - mU32*pow7(m32)))/(3.*mQ32*(
         mQ32 - mU32)*mU32*pow3(m32 - mQ32)*pow3(m32 - mU32)) - (8*(360*msq2*
         pow2(mU32)*pow3(m32)*pow3(mQ32) - 540*msq2*pow2(mQ32)*pow3(m32)*pow3(
         mU32) + 180*msq2*pow2(m32)*pow3(mQ32)*pow3(mU32) - 8861*pow3(m32)*pow3
         (mQ32)*pow3(mU32) - 120*msq2*pow2(mQ32)*pow2(mU32)*pow4(m32) - 540*
         msq2*mU32*pow3(mQ32)*pow4(m32) + 7583*pow2(mU32)*pow3(mQ32)*pow4(m32)
         - 1020*mQ32*msq2*pow3(mU32)*pow4(m32) + 15379*pow2(mQ32)*pow3(mU32)*
         pow4(m32) - 120*msq2*pow2(m32)*pow2(mU32)*pow4(mQ32) + 180*msq2*mU32*
         pow3(m32)*pow4(mQ32) - 1545*pow2(mU32)*pow3(m32)*pow4(mQ32) - 60*m32*
         msq2*pow3(mU32)*pow4(mQ32) + 1811*pow2(m32)*pow3(mU32)*pow4(mQ32) +
         429*mU32*pow4(m32)*pow4(mQ32) + 240*msq2*pow2(m32)*pow2(mQ32)*pow4(
         mU32) + 720*mQ32*msq2*pow3(m32)*pow4(mU32) - 7381*pow2(mQ32)*pow3(m32)
         *pow4(mU32) + 3655*pow2(m32)*pow3(mQ32)*pow4(mU32) + 2335*mQ32*pow4(
         m32)*pow4(mU32) - 641*m32*pow4(mQ32)*pow4(mU32) + 480*msq2*mU32*pow2(
         mQ32)*pow5(m32) + 840*mQ32*msq2*pow2(mU32)*pow5(m32) - 12656*pow2(mQ32
         )*pow2(mU32)*pow5(m32) - 2253*mU32*pow3(mQ32)*pow5(m32) - 5423*mQ32*
         pow3(mU32)*pow5(m32) - 2*pow4(mQ32)*pow5(m32) - 136*pow4(mU32)*pow5(
         m32) - 12*pow2(m32)*pow2(mU32)*pow5(mQ32) + 18*mU32*pow3(m32)*pow5(
         mQ32) - 6*m32*pow3(mU32)*pow5(mQ32) - 180*mQ32*msq2*pow2(m32)*pow5(
         mU32) - 60*m32*msq2*pow2(mQ32)*pow5(mU32) + 1341*pow2(m32)*pow2(mQ32)*
         pow5(mU32) - 289*mQ32*pow3(m32)*pow5(mU32) - 538*m32*pow3(mQ32)*pow5(
         mU32) + 34*pow4(m32)*pow5(mU32) + 84*pow4(mQ32)*pow5(mU32) - 360*mQ32*
         msq2*mU32*pow6(m32) + 3740*mU32*pow2(mQ32)*pow6(m32) + 4714*mQ32*pow2(
         mU32)*pow6(m32) + 6*pow3(mQ32)*pow6(m32) + 172*pow3(mU32)*pow6(m32) -
         27*mQ32*pow2(m32)*pow6(mU32) - 9*m32*pow2(mQ32)*pow6(mU32) - 1448*mQ32
         *mU32*pow7(m32) - 6*pow2(mQ32)*pow7(m32) - 72*pow2(mU32)*pow7(m32) + 2
         *mQ32*pow8(m32) + 2*mU32*pow8(m32)))/(3.*mQ32*mU32*pow3(m32 - mQ32)*
         pow4(m32 - mU32)) - (16*pow4(Xt)*(-90*msq2*pow3(m32)*pow3(mQ32)*pow3(
         mU32) - 330*msq2*pow2(mU32)*pow3(mQ32)*pow4(m32) - 570*msq2*pow2(mQ32)
         *pow3(mU32)*pow4(m32) + 16894*pow3(mQ32)*pow3(mU32)*pow4(m32) + 270*
         msq2*pow2(mU32)*pow3(m32)*pow4(mQ32) + 30*msq2*pow2(m32)*pow3(mU32)*
         pow4(mQ32) - 6970*pow3(m32)*pow3(mU32)*pow4(mQ32) - 270*msq2*mU32*pow4
         (m32)*pow4(mQ32) + 5964*pow2(mU32)*pow4(m32)*pow4(mQ32) + 90*msq2*pow2
         (mQ32)*pow3(m32)*pow4(mU32) + 210*msq2*pow2(m32)*pow3(mQ32)*pow4(mU32)
         - 12199*pow3(m32)*pow3(mQ32)*pow4(mU32) - 510*mQ32*msq2*pow4(m32)*
         pow4(mU32) + 15469*pow2(mQ32)*pow4(m32)*pow4(mU32) - 30*m32*msq2*pow4(
         mQ32)*pow4(mU32) + 3405*pow2(m32)*pow4(mQ32)*pow4(mU32) + 660*msq2*
         pow2(mQ32)*pow2(mU32)*pow5(m32) + 240*msq2*mU32*pow3(mQ32)*pow5(m32) -
         10834*pow2(mU32)*pow3(mQ32)*pow5(m32) + 420*mQ32*msq2*pow3(mU32)*pow5
         (m32) - 15617*pow2(mQ32)*pow3(mU32)*pow5(m32) - 1653*mU32*pow4(mQ32)*
         pow5(m32) - 5675*mQ32*pow4(mU32)*pow5(m32) - 60*msq2*pow2(m32)*pow2(
         mU32)*pow5(mQ32) + 90*msq2*mU32*pow3(m32)*pow5(mQ32) - 1187*pow2(mU32)
         *pow3(m32)*pow5(mQ32) - 30*m32*msq2*pow3(mU32)*pow5(mQ32) + 1160*pow2(
         m32)*pow3(mU32)*pow5(mQ32) + 318*mU32*pow4(m32)*pow5(mQ32) - 330*m32*
         pow4(mU32)*pow5(mQ32) - pow5(m32)*pow5(mQ32) + 30*msq2*pow2(m32)*pow2(
         mQ32)*pow5(mU32) + 360*mQ32*msq2*pow3(m32)*pow5(mU32) - 6641*pow2(mQ32
         )*pow3(m32)*pow5(mU32) - 30*m32*msq2*pow3(mQ32)*pow5(mU32) + 4122*pow2
         (m32)*pow3(mQ32)*pow5(mU32) + 2428*mQ32*pow4(m32)*pow5(mU32) - 901*m32
         *pow4(mQ32)*pow5(mU32) + 68*pow5(m32)*pow5(mU32) + 48*pow5(mQ32)*pow5(
         mU32) - 180*msq2*mU32*pow2(mQ32)*pow6(m32) - 180*mQ32*msq2*pow2(mU32)*
         pow6(m32) + 6667*pow2(mQ32)*pow2(mU32)*pow6(m32) + 2556*mU32*pow3(mQ32
         )*pow6(m32) + 4964*mQ32*pow3(mU32)*pow6(m32) + 3*pow4(mQ32)*pow6(m32)
         - 86*pow4(mU32)*pow6(m32) - 6*pow2(m32)*pow2(mU32)*pow6(mQ32) + 9*mU32
         *pow3(m32)*pow6(mQ32) - 3*m32*pow3(mU32)*pow6(mQ32) - 90*mQ32*msq2*
         pow2(m32)*pow6(mU32) - 30*m32*msq2*pow2(mQ32)*pow6(mU32) + 993*pow2(
         m32)*pow2(mQ32)*pow6(mU32) - 268*mQ32*pow3(m32)*pow6(mU32) - 488*m32*
         pow3(mQ32)*pow6(mU32) - 17*pow4(m32)*pow6(mU32) + 120*pow4(mQ32)*pow6(
         mU32) - 898*mU32*pow2(mQ32)*pow7(m32) - 1335*mQ32*pow2(mU32)*pow7(m32)
         - 3*pow3(mQ32)*pow7(m32) + 36*pow3(mU32)*pow7(m32) - 18*mQ32*pow2(m32
         )*pow7(mU32) - 6*m32*pow2(mQ32)*pow7(mU32) - 88*mQ32*mU32*pow8(m32) +
         pow2(mQ32)*pow8(m32) - pow2(mU32)*pow8(m32)))/(3.*mQ32*mU32*pow3(m32 -
         mQ32)*pow3(mQ32 - mU32)*pow4(m32 - mU32)) + (64*Xt*(22*mQ32*pow11(m3)
         + 16*mU32*pow11(m3) - 60*msq2*mU32*pow3(m3)*pow3(mQ32) + 326*pow2(
         mU32)*pow3(m3)*pow3(mQ32) + 60*mQ32*msq2*pow3(m3)*pow3(mU32) - 309*
         pow2(mQ32)*pow3(m3)*pow3(mU32) - 12*m3*pow3(mQ32)*pow3(mU32) - 6*mU32*
         pow3(m3)*pow4(mQ32) + 9*mQ32*pow3(m3)*pow4(mU32) + 120*msq2*mU32*pow2(
         mQ32)*pow5(m3) - 180*mQ32*msq2*pow2(mU32)*pow5(m3) + 355*pow2(mQ32)*
         pow2(mU32)*pow5(m3) + 60*msq2*pow3(mQ32)*pow5(m3) - 529*mU32*pow3(mQ32
         )*pow5(m3) + 222*mQ32*pow3(mU32)*pow5(m3) + 6*pow4(mQ32)*pow5(m3) +
         120*mQ32*msq2*mU32*pow7(m3) - 120*msq2*pow2(mQ32)*pow7(m3) + 21*mU32*
         pow2(mQ32)*pow7(m3) - 414*mQ32*pow2(mU32)*pow7(m3) + 295*pow3(mQ32)*
         pow7(m3) + 16*pow3(mU32)*pow7(m3) + 241*mQ32*mU32*pow9(m3) - 227*pow2(
         mQ32)*pow9(m3) - 32*pow2(mU32)*pow9(m3)))/(3.*mQ32*(mQ32 - mU32)*pow2(
         m32 - mQ32)*pow3(m32 - mU32)) + (128*pow5(Xt)*(8*mQ32*pow11(m3) + 8*
         mU32*pow11(m3) - 60*msq2*pow2(mQ32)*pow2(mU32)*pow3(m3) - 30*msq2*mU32
         *pow3(m3)*pow3(mQ32) + 209*pow2(mU32)*pow3(m3)*pow3(mQ32) - 30*mQ32*
         msq2*pow3(m3)*pow3(mU32) + 221*pow2(mQ32)*pow3(m3)*pow3(mU32) + 6*m3*
         pow3(mQ32)*pow3(mU32) - 3*mU32*pow3(m3)*pow4(mQ32) - 6*mQ32*pow3(m3)*
         pow4(mU32) + 120*msq2*mU32*pow2(mQ32)*pow5(m3) + 90*mQ32*msq2*pow2(
         mU32)*pow5(m3) - 620*pow2(mQ32)*pow2(mU32)*pow5(m3) + 30*msq2*pow3(
         mQ32)*pow5(m3) - 359*mU32*pow3(mQ32)*pow5(m3) - 175*mQ32*pow3(mU32)*
         pow5(m3) + 3*pow4(mQ32)*pow5(m3) - 60*mQ32*msq2*mU32*pow7(m3) - 60*
         msq2*pow2(mQ32)*pow7(m3) + 498*mU32*pow2(mQ32)*pow7(m3) + 365*mQ32*
         pow2(mU32)*pow7(m3) + 160*pow3(mQ32)*pow7(m3) + 8*pow3(mU32)*pow7(m3)
         - 176*mQ32*mU32*pow9(m3) - 131*pow2(mQ32)*pow9(m3) - 16*pow2(mU32)*
         pow9(m3)))/(3.*mQ32*pow2(m32 - mQ32)*pow3(m32 - mU32)*pow3(mQ32 - mU32
         )) - (128*pow3(Xt)*(22*mQ32*pow11(m3) - 60*msq2*pow2(mQ32)*pow2(mU32)*
         pow3(m3) - 60*msq2*mU32*pow3(m3)*pow3(mQ32) + 640*pow2(mU32)*pow3(m3)*
         pow3(mQ32) + 60*mQ32*msq2*pow3(m3)*pow3(mU32) - 391*pow2(mQ32)*pow3(m3
         )*pow3(mU32) + 18*m3*pow3(mQ32)*pow3(mU32) + 6*m3*pow2(mU32)*pow4(mQ32
         ) + 60*msq2*pow3(m3)*pow4(mQ32) - 367*mU32*pow3(m3)*pow4(mQ32) + 6*
         mQ32*pow3(m3)*pow4(mU32) + 240*msq2*mU32*pow2(mQ32)*pow5(m3) - 120*
         mQ32*msq2*pow2(mU32)*pow5(m3) - 167*pow2(mQ32)*pow2(mU32)*pow5(m3) -
         120*msq2*pow3(mQ32)*pow5(m3) - 235*mU32*pow3(mQ32)*pow5(m3) + 315*mQ32
         *pow3(mU32)*pow5(m3) + 303*pow4(mQ32)*pow5(m3) + 6*pow3(m3)*pow5(mQ32)
         + 358*mU32*pow2(mQ32)*pow7(m3) - 277*mQ32*pow2(mU32)*pow7(m3) - 253*
         pow3(mQ32)*pow7(m3) - 16*pow3(mU32)*pow7(m3) + 6*mQ32*mU32*pow9(m3) +
         10*pow2(mQ32)*pow9(m3) + 16*pow2(mU32)*pow9(m3)))/(3.*mQ32*pow2(m32 -
         mQ32)*pow2(m32 - mU32)*pow3(mQ32 - mU32))) + pow2(log(mU32/MR2))*((
         -2560*mU32*(mQ32 + mU32)*pow2(m32)*pow6(Xt))/(3.*(m32 - mQ32)*pow2(m32
         - mU32)*pow4(mQ32 - mU32)) - (8*(195*pow2(mQ32)*pow2(mU32)*pow3(m32)
         - 61*pow2(m32)*pow2(mU32)*pow3(mQ32) - 71*mU32*pow3(m32)*pow3(mQ32) -
         369*pow2(m32)*pow2(mQ32)*pow3(mU32) + 531*mQ32*pow3(m32)*pow3(mU32) +
         119*m32*pow3(mQ32)*pow3(mU32) + 209*mU32*pow2(mQ32)*pow4(m32) - 433*
         mQ32*pow2(mU32)*pow4(m32) - 19*pow3(mQ32)*pow4(m32) - 209*pow3(mU32)*
         pow4(m32) - 125*mQ32*pow2(m32)*pow4(mU32) + 76*m32*pow2(mQ32)*pow4(
         mU32) + 25*pow3(m32)*pow4(mU32) - 24*pow3(mQ32)*pow4(mU32) - 38*mQ32*
         mU32*pow5(m32) + 57*pow2(mQ32)*pow5(m32) + 251*pow2(mU32)*pow5(m32) +
         3*m32*mQ32*pow5(mU32) + 9*pow2(m32)*pow5(mU32) - 106*mQ32*pow6(m32) -
         88*mU32*pow6(m32) + 68*pow7(m32)))/(3.*pow3(m32 - mQ32)*pow4(m32 -
         mU32)) + (16*pow2(Xt)*(180*pow2(mU32)*pow3(m32)*pow3(mQ32) - 1246*pow2
         (mQ32)*pow3(m32)*pow3(mU32) - 26*pow2(m32)*pow3(mQ32)*pow3(mU32) + 956
         *pow2(mQ32)*pow2(mU32)*pow4(m32) - 246*mU32*pow3(mQ32)*pow4(m32) +
         2702*mQ32*pow3(mU32)*pow4(m32) - 71*pow2(m32)*pow2(mU32)*pow4(mQ32) -
         51*mU32*pow3(m32)*pow4(mQ32) + 119*m32*pow3(mU32)*pow4(mQ32) + 35*pow4
         (m32)*pow4(mQ32) + 692*pow2(m32)*pow2(mQ32)*pow4(mU32) - 1242*mQ32*
         pow3(m32)*pow4(mU32) - 139*m32*pow3(mQ32)*pow4(mU32) + 593*pow4(m32)*
         pow4(mU32) - 24*pow4(mQ32)*pow4(mU32) + 215*mU32*pow2(mQ32)*pow5(m32)
         - 2641*mQ32*pow2(mU32)*pow5(m32) + 151*pow3(mQ32)*pow5(m32) - 1429*
         pow3(mU32)*pow5(m32) + 238*mQ32*pow2(m32)*pow5(mU32) - 145*m32*pow2(
         mQ32)*pow5(mU32) - 81*pow3(m32)*pow5(mU32) + 48*pow3(mQ32)*pow5(mU32)
         + 676*mQ32*mU32*pow6(m32) - 424*pow2(mQ32)*pow6(m32) + 1516*pow2(mU32)
         *pow6(m32) - 3*m32*mQ32*pow6(mU32) - 9*pow2(m32)*pow6(mU32) + 238*mQ32
         *pow7(m32) - 582*mU32*pow7(m32)))/(3.*pow2(mQ32 - mU32)*pow3(m32 -
         mQ32)*pow4(m32 - mU32)) + (8*pow4(Xt)*(4362*pow3(m32)*pow3(mQ32)*pow3(
         mU32) - 26*pow2(mU32)*pow3(mQ32)*pow4(m32) - 15882*pow2(mQ32)*pow3(
         mU32)*pow4(m32) - 1513*pow2(mU32)*pow3(m32)*pow4(mQ32) + 673*pow2(m32)
         *pow3(mU32)*pow4(mQ32) + 747*mU32*pow4(m32)*pow4(mQ32) + 7300*pow2(
         mQ32)*pow3(m32)*pow4(mU32) - 3118*pow2(m32)*pow3(mQ32)*pow4(mU32) -
         5951*mQ32*pow4(m32)*pow4(mU32) + 176*m32*pow4(mQ32)*pow4(mU32) + 10006
         *pow2(mQ32)*pow2(mU32)*pow5(m32) - 2222*mU32*pow3(mQ32)*pow5(m32) +
         15702*mQ32*pow3(mU32)*pow5(m32) - 11*pow4(mQ32)*pow5(m32) + 1581*pow4(
         mU32)*pow5(m32) - 17*pow2(m32)*pow2(mU32)*pow5(mQ32) + 131*mU32*pow3(
         m32)*pow5(mQ32) - 71*m32*pow3(mU32)*pow5(mQ32) - 63*pow4(m32)*pow5(
         mQ32) + 12*pow4(mU32)*pow5(mQ32) - 562*pow2(m32)*pow2(mQ32)*pow5(mU32)
         + 395*mQ32*pow3(m32)*pow5(mU32) + 380*m32*pow3(mQ32)*pow5(mU32) - 105
         *pow4(m32)*pow5(mU32) - 48*pow4(mQ32)*pow5(mU32) - 44*mU32*pow2(mQ32)*
         pow6(m32) - 13092*mQ32*pow2(mU32)*pow6(m32) + 668*pow3(mQ32)*pow6(m32)
         - 4844*pow3(mU32)*pow6(m32) - 281*mQ32*pow2(m32)*pow6(mU32) + 184*m32
         *pow2(mQ32)*pow6(mU32) + 109*pow3(m32)*pow6(mU32) - 60*pow3(mQ32)*pow6
         (mU32) + 2824*mQ32*mU32*pow7(m32) - 1018*pow2(mQ32)*pow7(m32) + 4690*
         pow2(mU32)*pow7(m32) + 3*m32*mQ32*pow7(mU32) + 9*pow2(m32)*pow7(mU32)
         + 424*mQ32*pow8(m32) - 1448*mU32*pow8(m32)))/(3.*pow3(m32 - mQ32)*pow4
         (m32 - mU32)*pow4(mQ32 - mU32)) - (32*Xt*(54*pow11(m3) - 43*pow2(mQ32)
         *pow2(mU32)*pow3(m3) - 12*m3*pow2(mQ32)*pow3(mU32) - 25*mQ32*pow3(m3)*
         pow3(mU32) + 6*pow3(m3)*pow4(mU32) + 80*mU32*pow2(mQ32)*pow5(m3) + 217
         *mQ32*pow2(mU32)*pow5(m3) + 3*pow3(mU32)*pow5(m3) - 291*mQ32*mU32*pow7
         (m3) + 71*pow2(mQ32)*pow7(m3) - 92*pow2(mU32)*pow7(m3) - 93*mQ32*pow9(
         m3) + 125*mU32*pow9(m3)))/(3.*(mQ32 - mU32)*pow2(m32 - mQ32)*pow3(m32
         - mU32)) + (32*pow5(Xt)*(64*mQ32*pow11(m3) + 49*pow2(mU32)*pow3(m3)*
         pow3(mQ32) + 122*pow2(mQ32)*pow3(m3)*pow3(mU32) + 12*m3*pow3(mQ32)*
         pow3(mU32) + 12*m3*pow2(mQ32)*pow4(mU32) + 67*mQ32*pow3(m3)*pow4(mU32)
         - 401*pow2(mQ32)*pow2(mU32)*pow5(m3) - 60*mU32*pow3(mQ32)*pow5(m3) -
         328*mQ32*pow3(mU32)*pow5(m3) - 51*pow4(mU32)*pow5(m3) - 6*pow3(m3)*
         pow5(mU32) + 330*mU32*pow2(mQ32)*pow7(m3) + 413*mQ32*pow2(mU32)*pow7(
         m3) + 31*pow3(mQ32)*pow7(m3) + 178*pow3(mU32)*pow7(m3) - 248*mQ32*mU32
         *pow9(m3) - 95*pow2(mQ32)*pow9(m3) - 89*pow2(mU32)*pow9(m3)))/(3.*pow2
         (m32 - mQ32)*pow3(m32 - mU32)*pow4(mQ32 - mU32)) - (64*pow3(Xt)*(-42*
         mQ32*pow11(m3) - 84*mU32*pow11(m3) + 42*pow13(m3) + 15*pow2(mU32)*pow3
         (m3)*pow3(mQ32) + 16*pow2(mQ32)*pow3(m3)*pow3(mU32) + 18*m3*pow3(mQ32)
         *pow3(mU32) - 6*m3*pow2(mQ32)*pow4(mU32) - 75*mQ32*pow3(m3)*pow4(mU32)
         - 237*pow2(mQ32)*pow2(mU32)*pow5(m3) - 2*mU32*pow3(mQ32)*pow5(m3) +
         232*mQ32*pow3(mU32)*pow5(m3) + 49*pow4(mU32)*pow5(m3) + 6*pow3(m3)*
         pow5(mU32) + 124*mU32*pow2(mQ32)*pow7(m3) + 29*mQ32*pow2(mU32)*pow7(m3
         ) + pow3(mQ32)*pow7(m3) - 214*pow3(mU32)*pow7(m3) - 48*mQ32*mU32*pow9(
         m3) + 7*pow2(mQ32)*pow9(m3) + 169*pow2(mU32)*pow9(m3)))/(3.*pow2(m32 -
         mQ32)*pow3(m32 - mU32)*pow3(mQ32 - mU32))) + log(mQ32/MR2)*((5120*
         mQ32*pow2(m32)*pow6(Xt))/(3.*(m32 - mU32)*pow2(m32 - mQ32)*pow3(mQ32 -
         mU32)) + log(msq2/MR2)*((-16*(45*msq2*mU32*pow2(m32)*pow2(mQ32) + 45*
         mQ32*msq2*pow2(m32)*pow2(mU32) - 90*mQ32*msq2*mU32*pow3(m32) + 135*
         msq2*pow2(mQ32)*pow3(m32) - 15*mU32*pow2(mQ32)*pow3(m32) + 12*mQ32*
         pow2(mU32)*pow3(m32) + 135*msq2*pow2(mU32)*pow3(m32) - 15*m32*msq2*
         mU32*pow3(mQ32) - 45*msq2*pow2(m32)*pow3(mQ32) + 5*mU32*pow2(m32)*pow3
         (mQ32) - 5*pow3(m32)*pow3(mQ32) - 15*m32*mQ32*msq2*pow3(mU32) - 4*mQ32
         *pow2(m32)*pow3(mU32) - 45*msq2*pow2(m32)*pow3(mU32) + 4*pow3(m32)*
         pow3(mU32) - 120*mQ32*msq2*pow4(m32) + 3*mQ32*mU32*pow4(m32) - 120*
         msq2*mU32*pow4(m32) + 15*pow2(mQ32)*pow4(m32) - 12*pow2(mU32)*pow4(m32
         ) - 11*mQ32*pow5(m32) + 90*msq2*pow5(m32) + 7*mU32*pow5(m32) + pow6(
         m32)))/(3.*pow3(m32 - mQ32)*pow3(m32 - mU32)) - (160*pow2(Xt)*(9*msq2*
         mU32*pow2(m32)*pow2(mQ32) + 9*mQ32*msq2*pow2(m32)*pow2(mU32) - 18*mQ32
         *msq2*mU32*pow3(m32) + 27*msq2*pow2(mQ32)*pow3(m32) - 3*mU32*pow2(mQ32
         )*pow3(m32) - 3*mQ32*pow2(mU32)*pow3(m32) + 27*msq2*pow2(mU32)*pow3(
         m32) - 3*m32*msq2*mU32*pow3(mQ32) - 9*msq2*pow2(m32)*pow3(mQ32) + mU32
         *pow2(m32)*pow3(mQ32) - pow3(m32)*pow3(mQ32) - 3*m32*mQ32*msq2*pow3(
         mU32) + mQ32*pow2(m32)*pow3(mU32) - 9*msq2*pow2(m32)*pow3(mU32) - pow3
         (m32)*pow3(mU32) - 24*mQ32*msq2*pow4(m32) + 6*mQ32*mU32*pow4(m32) - 24
         *msq2*mU32*pow4(m32) + 3*pow2(mQ32)*pow4(m32) + 3*pow2(mU32)*pow4(m32)
         - 4*mQ32*pow5(m32) + 18*msq2*pow5(m32) - 4*mU32*pow5(m32) + 2*pow6(
         m32)))/(3.*(mQ32 - mU32)*pow3(m32 - mQ32)*pow3(m32 - mU32)) + (80*(
         mQ32 + mU32)*pow4(Xt)*(9*msq2*mU32*pow2(m32)*pow2(mQ32) + 9*mQ32*msq2*
         pow2(m32)*pow2(mU32) - 18*mQ32*msq2*mU32*pow3(m32) + 27*msq2*pow2(mQ32
         )*pow3(m32) - 3*mU32*pow2(mQ32)*pow3(m32) - 3*mQ32*pow2(mU32)*pow3(m32
         ) + 27*msq2*pow2(mU32)*pow3(m32) - 3*m32*msq2*mU32*pow3(mQ32) - 9*msq2
         *pow2(m32)*pow3(mQ32) + mU32*pow2(m32)*pow3(mQ32) - pow3(m32)*pow3(
         mQ32) - 3*m32*mQ32*msq2*pow3(mU32) + mQ32*pow2(m32)*pow3(mU32) - 9*
         msq2*pow2(m32)*pow3(mU32) - pow3(m32)*pow3(mU32) - 24*mQ32*msq2*pow4(
         m32) + 6*mQ32*mU32*pow4(m32) - 24*msq2*mU32*pow4(m32) + 3*pow2(mQ32)*
         pow4(m32) + 3*pow2(mU32)*pow4(m32) - 4*mQ32*pow5(m32) + 18*msq2*pow5(
         m32) - 4*mU32*pow5(m32) + 2*pow6(m32)))/(3.*pow3(m32 - mQ32)*pow3(m32
         - mU32)*pow3(mQ32 - mU32)) - (320*(mQ32 + mU32)*pow5(Xt)*(-6*mQ32*msq2
         *pow3(m3) + mQ32*mU32*pow3(m3) - 6*msq2*mU32*pow3(m3) - mQ32*pow5(m3)
         + 12*msq2*pow5(m3) - mU32*pow5(m3) + pow7(m3)))/(3.*pow2(m32 - mQ32)*
         pow2(m32 - mU32)*pow3(mQ32 - mU32)) + Xt*((96*(-2*m32 + mQ32 + mU32)*
         pow3(m3))/((m32 - mQ32)*(m32 - mU32)*(mQ32 - mU32)) + (32*(-60*mQ32*
         msq2*pow3(m3) + mQ32*mU32*pow3(m3) - 60*msq2*mU32*pow3(m3) - mQ32*pow5
         (m3) + 120*msq2*pow5(m3) - mU32*pow5(m3) + pow7(m3)))/(3.*pow2(m32 -
         mQ32)*pow2(m32 - mU32))) + pow3(Xt)*((-192*(4*mQ32*mU32 - 2*m32*(mQ32
         + mU32) + 2*pow2(m32) - pow2(mQ32) - pow2(mU32))*pow3(m3))/((m32 -
         mQ32)*(m32 - mU32)*pow3(mQ32 - mU32)) + (64*(-60*mQ32*msq2*pow3(m3) +
         mQ32*mU32*pow3(m3) - 60*msq2*mU32*pow3(m3) - mQ32*pow5(m3) + 120*msq2*
         pow5(m3) - mU32*pow5(m3) + pow7(m3)))/(3.*(mQ32 - mU32)*pow2(m32 -
         mQ32)*pow2(m32 - mU32)))) + (32*pow2(Xt)*(-180*msq2*pow2(mQ32)*pow2(
         mU32)*pow3(m32) + 90*msq2*pow2(m32)*pow2(mU32)*pow3(mQ32) + 270*msq2*
         mU32*pow3(m32)*pow3(mQ32) - 3159*pow2(mU32)*pow3(m32)*pow3(mQ32) + 90*
         msq2*pow2(m32)*pow2(mQ32)*pow3(mU32) + 270*mQ32*msq2*pow3(m32)*pow3(
         mU32) - 3399*pow2(mQ32)*pow3(m32)*pow3(mU32) + 1626*pow2(m32)*pow3(
         mQ32)*pow3(mU32) - 240*msq2*mU32*pow2(mQ32)*pow4(m32) - 240*mQ32*msq2*
         pow2(mU32)*pow4(m32) + 5392*pow2(mQ32)*pow2(mU32)*pow4(m32) + 1008*
         mU32*pow3(mQ32)*pow4(m32) + 1392*mQ32*pow3(mU32)*pow4(m32) - 90*msq2*
         mU32*pow2(m32)*pow4(mQ32) - 30*m32*msq2*pow2(mU32)*pow4(mQ32) + 690*
         pow2(m32)*pow2(mU32)*pow4(mQ32) - 155*mU32*pow3(m32)*pow4(mQ32) - 285*
         m32*pow3(mU32)*pow4(mQ32) - 18*pow4(m32)*pow4(mQ32) - 90*mQ32*msq2*
         pow2(m32)*pow4(mU32) - 30*m32*msq2*pow2(mQ32)*pow4(mU32) + 750*pow2(
         m32)*pow2(mQ32)*pow4(mU32) - 255*mQ32*pow3(m32)*pow4(mU32) - 297*m32*
         pow3(mQ32)*pow4(mU32) + 2*pow4(m32)*pow4(mU32) + 48*pow4(mQ32)*pow4(
         mU32) + 180*mQ32*msq2*mU32*pow5(m32) - 1958*mU32*pow2(mQ32)*pow5(m32)
         - 2278*mQ32*pow2(mU32)*pow5(m32) + 54*pow3(mQ32)*pow5(m32) - 6*pow3(
         mU32)*pow5(m32) - 9*mU32*pow2(m32)*pow5(mQ32) - 3*m32*pow2(mU32)*pow5(
         mQ32) - 9*mQ32*pow2(m32)*pow5(mU32) - 3*m32*pow2(mQ32)*pow5(mU32) +
         904*mQ32*mU32*pow6(m32) - 38*pow2(mQ32)*pow6(m32) + 6*pow2(mU32)*pow6(
         m32) + 2*mQ32*pow7(m32) - 2*mU32*pow7(m32)))/(3.*mQ32*(mQ32 - mU32)*
         mU32*pow3(m32 - mQ32)*pow3(m32 - mU32)) - (8*(-540*msq2*pow2(mU32)*
         pow3(m32)*pow3(mQ32) + 360*msq2*pow2(mQ32)*pow3(m32)*pow3(mU32) + 180*
         msq2*pow2(m32)*pow3(mQ32)*pow3(mU32) - 8855*pow3(m32)*pow3(mQ32)*pow3(
         mU32) - 120*msq2*pow2(mQ32)*pow2(mU32)*pow4(m32) - 1020*msq2*mU32*pow3
         (mQ32)*pow4(m32) + 15361*pow2(mU32)*pow3(mQ32)*pow4(m32) - 540*mQ32*
         msq2*pow3(mU32)*pow4(m32) + 7577*pow2(mQ32)*pow3(mU32)*pow4(m32) + 240
         *msq2*pow2(m32)*pow2(mU32)*pow4(mQ32) + 720*msq2*mU32*pow3(m32)*pow4(
         mQ32) - 7373*pow2(mU32)*pow3(m32)*pow4(mQ32) + 3655*pow2(m32)*pow3(
         mU32)*pow4(mQ32) + 2335*mU32*pow4(m32)*pow4(mQ32) - 120*msq2*pow2(m32)
         *pow2(mQ32)*pow4(mU32) + 180*mQ32*msq2*pow3(m32)*pow4(mU32) - 1543*
         pow2(mQ32)*pow3(m32)*pow4(mU32) - 60*m32*msq2*pow3(mQ32)*pow4(mU32) +
         1809*pow2(m32)*pow3(mQ32)*pow4(mU32) + 431*mQ32*pow4(m32)*pow4(mU32) -
         641*m32*pow4(mQ32)*pow4(mU32) + 840*msq2*mU32*pow2(mQ32)*pow5(m32) +
         480*mQ32*msq2*pow2(mU32)*pow5(m32) - 12642*pow2(mQ32)*pow2(mU32)*pow5(
         m32) - 5421*mU32*pow3(mQ32)*pow5(m32) - 2259*mQ32*pow3(mU32)*pow5(m32)
         - 144*pow4(mQ32)*pow5(m32) - 4*pow4(mU32)*pow5(m32) - 180*msq2*mU32*
         pow2(m32)*pow5(mQ32) - 60*m32*msq2*pow2(mU32)*pow5(mQ32) + 1339*pow2(
         m32)*pow2(mU32)*pow5(mQ32) - 289*mU32*pow3(m32)*pow5(mQ32) - 538*m32*
         pow3(mU32)*pow5(mQ32) + 36*pow4(m32)*pow5(mQ32) + 84*pow4(mU32)*pow5(
         mQ32) - 12*pow2(m32)*pow2(mQ32)*pow5(mU32) + 18*mQ32*pow3(m32)*pow5(
         mU32) - 6*m32*pow3(mQ32)*pow5(mU32) - 360*mQ32*msq2*mU32*pow6(m32) +
         4712*mU32*pow2(mQ32)*pow6(m32) + 3744*mQ32*pow2(mU32)*pow6(m32) + 184*
         pow3(mQ32)*pow6(m32) + 12*pow3(mU32)*pow6(m32) - 27*mU32*pow2(m32)*
         pow6(mQ32) - 9*m32*pow2(mU32)*pow6(mQ32) - 1450*mQ32*mU32*pow7(m32) -
         80*pow2(mQ32)*pow7(m32) - 12*pow2(mU32)*pow7(m32) + 4*mQ32*pow8(m32) +
         4*mU32*pow8(m32)))/(3.*mQ32*mU32*pow3(m32 - mU32)*pow4(m32 - mQ32)) -
         (16*pow4(Xt)*(90*msq2*pow3(m32)*pow3(mQ32)*pow3(mU32) + 570*msq2*pow2
         (mU32)*pow3(mQ32)*pow4(m32) + 330*msq2*pow2(mQ32)*pow3(mU32)*pow4(m32)
         - 16902*pow3(mQ32)*pow3(mU32)*pow4(m32) - 90*msq2*pow2(mU32)*pow3(m32
         )*pow4(mQ32) - 210*msq2*pow2(m32)*pow3(mU32)*pow4(mQ32) + 12196*pow3(
         m32)*pow3(mU32)*pow4(mQ32) + 510*msq2*mU32*pow4(m32)*pow4(mQ32) -
         15452*pow2(mU32)*pow4(m32)*pow4(mQ32) - 270*msq2*pow2(mQ32)*pow3(m32)*
         pow4(mU32) - 30*msq2*pow2(m32)*pow3(mQ32)*pow4(mU32) + 6982*pow3(m32)*
         pow3(mQ32)*pow4(mU32) + 270*mQ32*msq2*pow4(m32)*pow4(mU32) - 5982*pow2
         (mQ32)*pow4(m32)*pow4(mU32) + 30*m32*msq2*pow4(mQ32)*pow4(mU32) - 3408
         *pow2(m32)*pow4(mQ32)*pow4(mU32) - 660*msq2*pow2(mQ32)*pow2(mU32)*pow5
         (m32) - 420*msq2*mU32*pow3(mQ32)*pow5(m32) + 15609*pow2(mU32)*pow3(
         mQ32)*pow5(m32) - 240*mQ32*msq2*pow3(mU32)*pow5(m32) + 10851*pow2(mQ32
         )*pow3(mU32)*pow5(m32) + 5657*mU32*pow4(mQ32)*pow5(m32) + 1665*mQ32*
         pow4(mU32)*pow5(m32) - 30*msq2*pow2(m32)*pow2(mU32)*pow5(mQ32) - 360*
         msq2*mU32*pow3(m32)*pow5(mQ32) + 6629*pow2(mU32)*pow3(m32)*pow5(mQ32)
         + 30*m32*msq2*pow3(mU32)*pow5(mQ32) - 4118*pow2(m32)*pow3(mU32)*pow5(
         mQ32) - 2416*mU32*pow4(m32)*pow5(mQ32) + 901*m32*pow4(mU32)*pow5(mQ32)
         - 72*pow5(m32)*pow5(mQ32) + 60*msq2*pow2(m32)*pow2(mQ32)*pow5(mU32) -
         90*mQ32*msq2*pow3(m32)*pow5(mU32) + 1193*pow2(mQ32)*pow3(m32)*pow5(
         mU32) + 30*m32*msq2*pow3(mQ32)*pow5(mU32) - 1164*pow2(m32)*pow3(mQ32)*
         pow5(mU32) - 322*mQ32*pow4(m32)*pow5(mU32) + 331*m32*pow4(mQ32)*pow5(
         mU32) + 2*pow5(m32)*pow5(mU32) - 48*pow5(mQ32)*pow5(mU32) + 180*msq2*
         mU32*pow2(mQ32)*pow6(m32) + 180*mQ32*msq2*pow2(mU32)*pow6(m32) - 6670*
         pow2(mQ32)*pow2(mU32)*pow6(m32) - 4952*mU32*pow3(mQ32)*pow6(m32) -
         2568*mQ32*pow3(mU32)*pow6(m32) + 92*pow4(mQ32)*pow6(m32) - 6*pow4(mU32
         )*pow6(m32) + 90*msq2*mU32*pow2(m32)*pow6(mQ32) + 30*m32*msq2*pow2(
         mU32)*pow6(mQ32) - 990*pow2(m32)*pow2(mU32)*pow6(mQ32) + 265*mU32*pow3
         (m32)*pow6(mQ32) + 487*m32*pow3(mU32)*pow6(mQ32) + 18*pow4(m32)*pow6(
         mQ32) - 120*pow4(mU32)*pow6(mQ32) + 6*pow2(m32)*pow2(mQ32)*pow6(mU32)
         - 9*mQ32*pow3(m32)*pow6(mU32) + 3*m32*pow3(mQ32)*pow6(mU32) + 1332*
         mU32*pow2(mQ32)*pow7(m32) + 902*mQ32*pow2(mU32)*pow7(m32) - 40*pow3(
         mQ32)*pow7(m32) + 6*pow3(mU32)*pow7(m32) + 18*mU32*pow2(m32)*pow7(mQ32
         ) + 6*m32*pow2(mU32)*pow7(mQ32) + 88*mQ32*mU32*pow8(m32) + 2*pow2(mQ32
         )*pow8(m32) - 2*pow2(mU32)*pow8(m32)))/(3.*mQ32*mU32*pow3(m32 - mU32)*
         pow3(mQ32 - mU32)*pow4(m32 - mQ32)) - (32*Xt*(32*mQ32*pow11(m3) + 44*
         mU32*pow11(m3) + 120*msq2*mU32*pow3(m3)*pow3(mQ32) - 617*pow2(mU32)*
         pow3(m3)*pow3(mQ32) - 120*mQ32*msq2*pow3(m3)*pow3(mU32) + 651*pow2(
         mQ32)*pow3(m3)*pow3(mU32) - 24*m3*pow3(mQ32)*pow3(mU32) + 18*mU32*pow3
         (m3)*pow4(mQ32) - 12*mQ32*pow3(m3)*pow4(mU32) - 360*msq2*mU32*pow2(
         mQ32)*pow5(m3) + 240*mQ32*msq2*pow2(mU32)*pow5(m3) + 709*pow2(mQ32)*
         pow2(mU32)*pow5(m3) + 443*mU32*pow3(mQ32)*pow5(m3) - 1056*mQ32*pow3(
         mU32)*pow5(m3) + 120*msq2*pow3(mU32)*pow5(m3) + 12*pow4(mU32)*pow5(m3)
         + 240*mQ32*msq2*mU32*pow7(m3) - 826*mU32*pow2(mQ32)*pow7(m3) + 41*
         mQ32*pow2(mU32)*pow7(m3) - 240*msq2*pow2(mU32)*pow7(m3) + 32*pow3(mQ32
         )*pow7(m3) + 589*pow3(mU32)*pow7(m3) + 481*mQ32*mU32*pow9(m3) - 64*
         pow2(mQ32)*pow9(m3) - 453*pow2(mU32)*pow9(m3)))/(3.*(mQ32 - mU32)*mU32
         *pow2(m32 - mU32)*pow3(m32 - mQ32)) - (128*pow5(Xt)*(8*mQ32*pow11(m3)
         + 8*mU32*pow11(m3) - 60*msq2*pow2(mQ32)*pow2(mU32)*pow3(m3) - 30*msq2*
         mU32*pow3(m3)*pow3(mQ32) + 221*pow2(mU32)*pow3(m3)*pow3(mQ32) - 30*
         mQ32*msq2*pow3(m3)*pow3(mU32) + 209*pow2(mQ32)*pow3(m3)*pow3(mU32) + 6
         *m3*pow3(mQ32)*pow3(mU32) - 6*mU32*pow3(m3)*pow4(mQ32) - 3*mQ32*pow3(
         m3)*pow4(mU32) + 90*msq2*mU32*pow2(mQ32)*pow5(m3) + 120*mQ32*msq2*pow2
         (mU32)*pow5(m3) - 620*pow2(mQ32)*pow2(mU32)*pow5(m3) - 175*mU32*pow3(
         mQ32)*pow5(m3) - 359*mQ32*pow3(mU32)*pow5(m3) + 30*msq2*pow3(mU32)*
         pow5(m3) + 3*pow4(mU32)*pow5(m3) - 60*mQ32*msq2*mU32*pow7(m3) + 365*
         mU32*pow2(mQ32)*pow7(m3) + 498*mQ32*pow2(mU32)*pow7(m3) - 60*msq2*pow2
         (mU32)*pow7(m3) + 8*pow3(mQ32)*pow7(m3) + 160*pow3(mU32)*pow7(m3) -
         176*mQ32*mU32*pow9(m3) - 16*pow2(mQ32)*pow9(m3) - 131*pow2(mU32)*pow9(
         m3)))/(3.*mU32*pow2(m32 - mU32)*pow3(m32 - mQ32)*pow3(mQ32 - mU32)) +
         (128*pow3(Xt)*(22*mU32*pow11(m3) - 60*msq2*pow2(mQ32)*pow2(mU32)*pow3(
         m3) + 60*msq2*mU32*pow3(m3)*pow3(mQ32) - 391*pow2(mU32)*pow3(m3)*pow3(
         mQ32) - 60*mQ32*msq2*pow3(m3)*pow3(mU32) + 640*pow2(mQ32)*pow3(m3)*
         pow3(mU32) + 18*m3*pow3(mQ32)*pow3(mU32) + 6*mU32*pow3(m3)*pow4(mQ32)
         + 6*m3*pow2(mQ32)*pow4(mU32) - 367*mQ32*pow3(m3)*pow4(mU32) + 60*msq2*
         pow3(m3)*pow4(mU32) - 120*msq2*mU32*pow2(mQ32)*pow5(m3) + 240*mQ32*
         msq2*pow2(mU32)*pow5(m3) - 167*pow2(mQ32)*pow2(mU32)*pow5(m3) + 315*
         mU32*pow3(mQ32)*pow5(m3) - 235*mQ32*pow3(mU32)*pow5(m3) - 120*msq2*
         pow3(mU32)*pow5(m3) + 303*pow4(mU32)*pow5(m3) + 6*pow3(m3)*pow5(mU32)
         - 277*mU32*pow2(mQ32)*pow7(m3) + 358*mQ32*pow2(mU32)*pow7(m3) - 16*
         pow3(mQ32)*pow7(m3) - 253*pow3(mU32)*pow7(m3) + 6*mQ32*mU32*pow9(m3) +
         16*pow2(mQ32)*pow9(m3) + 10*pow2(mU32)*pow9(m3)))/(3.*mU32*pow2(m32 -
         mQ32)*pow2(m32 - mU32)*pow3(mQ32 - mU32)) + log(mU32/MR2)*((2560*(
         mQ32 + mU32)*(m32*mQ32 + m32*mU32 - 2*mQ32*mU32)*pow2(m32)*pow6(Xt))/(
         3.*pow2(m32 - mQ32)*pow2(m32 - mU32)*pow4(mQ32 - mU32)) + (8*(348*pow2
         (mU32)*pow3(m32)*pow3(mQ32) + 356*pow2(mQ32)*pow3(m32)*pow3(mU32) -
         592*pow2(m32)*pow3(mQ32)*pow3(mU32) + 292*pow2(mQ32)*pow2(mU32)*pow4(
         m32) + 388*mU32*pow3(mQ32)*pow4(m32) + 372*mQ32*pow3(mU32)*pow4(m32) -
         60*pow2(m32)*pow2(mU32)*pow4(mQ32) - 133*mU32*pow3(m32)*pow4(mQ32) +
         123*m32*pow3(mU32)*pow4(mQ32) + 102*pow4(m32)*pow4(mQ32) - 62*pow2(m32
         )*pow2(mQ32)*pow4(mU32) - 129*mQ32*pow3(m32)*pow4(mU32) + 123*m32*pow3
         (mQ32)*pow4(mU32) + 100*pow4(m32)*pow4(mU32) - 24*pow4(mQ32)*pow4(mU32
         ) - 1059*mU32*pow2(mQ32)*pow5(m32) - 1043*mQ32*pow2(mU32)*pow5(m32) -
         299*pow3(mQ32)*pow5(m32) - 291*pow3(mU32)*pow5(m32) + 6*mU32*pow2(m32)
         *pow5(mQ32) + 3*m32*pow2(mU32)*pow5(mQ32) - 9*pow3(m32)*pow5(mQ32) + 6
         *mQ32*pow2(m32)*pow5(mU32) + 3*m32*pow2(mQ32)*pow5(mU32) - 9*pow3(m32)
         *pow5(mU32) + 1196*mQ32*mU32*pow6(m32) + 518*pow2(mQ32)*pow6(m32) +
         508*pow2(mU32)*pow6(m32) - 434*mQ32*pow7(m32) - 430*mU32*pow7(m32) +
         130*pow8(m32)))/(3.*pow4(m32 - mQ32)*pow4(m32 - mU32)) + (16*pow2(Xt)*
         (-2056*pow3(m32)*pow3(mQ32)*pow3(mU32) + 4340*pow2(mU32)*pow3(mQ32)*
         pow4(m32) + 4340*pow2(mQ32)*pow3(mU32)*pow4(m32) - 1875*pow2(mU32)*
         pow3(m32)*pow4(mQ32) + 686*pow2(m32)*pow3(mU32)*pow4(mQ32) + 2024*mU32
         *pow4(m32)*pow4(mQ32) - 1875*pow2(mQ32)*pow3(m32)*pow4(mU32) + 686*
         pow2(m32)*pow3(mQ32)*pow4(mU32) + 2024*mQ32*pow4(m32)*pow4(mU32) - 230
         *m32*pow4(mQ32)*pow4(mU32) - 6042*pow2(mQ32)*pow2(mU32)*pow5(m32) -
         4054*mU32*pow3(mQ32)*pow5(m32) - 4054*mQ32*pow3(mU32)*pow5(m32) - 669*
         pow4(mQ32)*pow5(m32) - 669*pow4(mU32)*pow5(m32) + 312*pow2(m32)*pow2(
         mU32)*pow5(mQ32) - 370*mU32*pow3(m32)*pow5(mQ32) - 74*m32*pow3(mU32)*
         pow5(mQ32) + 116*pow4(m32)*pow5(mQ32) + 24*pow4(mU32)*pow5(mQ32) + 312
         *pow2(m32)*pow2(mQ32)*pow5(mU32) - 370*mQ32*pow3(m32)*pow5(mU32) - 74*
         m32*pow3(mQ32)*pow5(mU32) + 116*pow4(m32)*pow5(mU32) + 24*pow4(mQ32)*
         pow5(mU32) + 4234*mU32*pow2(mQ32)*pow6(m32) + 4234*mQ32*pow2(mU32)*
         pow6(m32) + 1238*pow3(mQ32)*pow6(m32) + 1238*pow3(mU32)*pow6(m32) - 6*
         mU32*pow2(m32)*pow6(mQ32) - 3*m32*pow2(mU32)*pow6(mQ32) + 9*pow3(m32)*
         pow6(mQ32) - 6*mQ32*pow2(m32)*pow6(mU32) - 3*m32*pow2(mQ32)*pow6(mU32)
         + 9*pow3(m32)*pow6(mU32) - 2132*mQ32*mU32*pow7(m32) - 1046*pow2(mQ32)
         *pow7(m32) - 1046*pow2(mU32)*pow7(m32) + 344*mQ32*pow8(m32) + 344*mU32
         *pow8(m32)))/(3.*pow2(mQ32 - mU32)*pow4(m32 - mQ32)*pow4(m32 - mU32))
         + (64*Xt*(-23*mQ32*pow11(m3) + 18*mU32*pow11(m3) + pow13(m3) - 44*pow2
         (mU32)*pow3(m3)*pow3(mQ32) + 43*pow2(mQ32)*pow3(m3)*pow3(mU32) - 3*
         mU32*pow3(m3)*pow4(mQ32) + 3*mQ32*pow3(m3)*pow4(mU32) + 3*pow2(mQ32)*
         pow2(mU32)*pow5(m3) + 52*mU32*pow3(mQ32)*pow5(m3) - 50*mQ32*pow3(mU32)
         *pow5(m3) + 3*pow4(mQ32)*pow5(m3) - 3*pow4(mU32)*pow5(m3) - 25*mU32*
         pow2(mQ32)*pow7(m3) + 16*mQ32*pow2(mU32)*pow7(m3) - 24*pow3(mQ32)*pow7
         (m3) + 23*pow3(mU32)*pow7(m3) + 6*mQ32*mU32*pow9(m3) + 27*pow2(mQ32)*
         pow9(m3) - 23*pow2(mU32)*pow9(m3)))/(3.*(mQ32 - mU32)*pow3(m32 - mQ32)
         *pow3(m32 - mU32)) + (128*pow3(Xt)*(9*mQ32*mU32*pow11(m3) - 3*mQ32*
         pow13(m3) - 3*mU32*pow13(m3) + pow15(m3) - 57*pow11(m3)*pow2(mQ32) +
         63*pow11(m3)*pow2(mU32) + pow3(m3)*pow3(mQ32)*pow3(mU32) - 42*pow2(
         mU32)*pow3(m3)*pow4(mQ32) - 12*m3*pow3(mU32)*pow4(mQ32) + 42*pow2(mQ32
         )*pow3(m3)*pow4(mU32) + 12*m3*pow3(mQ32)*pow4(mU32) + 231*pow2(mU32)*
         pow3(mQ32)*pow5(m3) - 237*pow2(mQ32)*pow3(mU32)*pow5(m3) + 63*mU32*
         pow4(mQ32)*pow5(m3) - 63*mQ32*pow4(mU32)*pow5(m3) + 3*mU32*pow3(m3)*
         pow5(mQ32) - 3*pow5(m3)*pow5(mQ32) - 3*mQ32*pow3(m3)*pow5(mU32) + 3*
         pow5(m3)*pow5(mU32) + 9*pow2(mQ32)*pow2(mU32)*pow7(m3) - 283*mU32*pow3
         (mQ32)*pow7(m3) + 289*mQ32*pow3(mU32)*pow7(m3) - 25*pow4(mQ32)*pow7(m3
         ) + 25*pow4(mU32)*pow7(m3) + 147*mU32*pow2(mQ32)*pow9(m3) - 165*mQ32*
         pow2(mU32)*pow9(m3) + 103*pow3(mQ32)*pow9(m3) - 105*pow3(mU32)*pow9(m3
         )))/(3.*pow3(m32 - mQ32)*pow3(m32 - mU32)*pow3(mQ32 - mU32)) - (64*
         pow5(Xt)*(-248*mQ32*mU32*pow11(m3) + 32*mQ32*pow13(m3) + 32*mU32*pow13
         (m3) - 124*pow11(m3)*pow2(mQ32) - 124*pow11(m3)*pow2(mU32) - 110*pow3(
         m3)*pow3(mQ32)*pow3(mU32) - 52*pow2(mU32)*pow3(m3)*pow4(mQ32) - 12*m3*
         pow3(mU32)*pow4(mQ32) - 52*pow2(mQ32)*pow3(m3)*pow4(mU32) - 12*m3*pow3
         (mQ32)*pow4(mU32) + 450*pow2(mU32)*pow3(mQ32)*pow5(m3) + 450*pow2(mQ32
         )*pow3(mU32)*pow5(m3) + 89*mU32*pow4(mQ32)*pow5(m3) + 89*mQ32*pow4(
         mU32)*pow5(m3) + 3*mU32*pow3(m3)*pow5(mQ32) - 3*pow5(m3)*pow5(mQ32) +
         3*mQ32*pow3(m3)*pow5(mU32) - 3*pow5(m3)*pow5(mU32) - 814*pow2(mQ32)*
         pow2(mU32)*pow7(m3) - 448*mU32*pow3(mQ32)*pow7(m3) - 448*mQ32*pow3(
         mU32)*pow7(m3) - 41*pow4(mQ32)*pow7(m3) - 41*pow4(mU32)*pow7(m3) + 540
         *mU32*pow2(mQ32)*pow9(m3) + 540*mQ32*pow2(mU32)*pow9(m3) + 152*pow3(
         mQ32)*pow9(m3) + 152*pow3(mU32)*pow9(m3)))/(3.*pow3(m32 - mQ32)*pow3(
         m32 - mU32)*pow4(mQ32 - mU32)) + (8*pow4(Xt)*(-40488*pow3(mQ32)*pow3(
         mU32)*pow4(m32) + 14107*pow3(m32)*pow3(mU32)*pow4(mQ32) - 11764*pow2(
         mU32)*pow4(m32)*pow4(mQ32) + 14107*pow3(m32)*pow3(mQ32)*pow4(mU32) -
         11764*pow2(mQ32)*pow4(m32)*pow4(mU32) - 6588*pow2(m32)*pow4(mQ32)*pow4
         (mU32) + 41616*pow2(mU32)*pow3(mQ32)*pow5(m32) + 41616*pow2(mQ32)*pow3
         (mU32)*pow5(m32) + 4563*mU32*pow4(mQ32)*pow5(m32) + 4563*mQ32*pow4(
         mU32)*pow5(m32) - 539*pow2(mU32)*pow3(m32)*pow5(mQ32) - 198*pow2(m32)*
         pow3(mU32)*pow5(mQ32) + 116*mU32*pow4(m32)*pow5(mQ32) + 592*m32*pow4(
         mU32)*pow5(mQ32) + 157*pow5(m32)*pow5(mQ32) - 539*pow2(mQ32)*pow3(m32)
         *pow5(mU32) - 198*pow2(m32)*pow3(mQ32)*pow5(mU32) + 116*mQ32*pow4(m32)
         *pow5(mU32) + 592*m32*pow4(mQ32)*pow5(mU32) + 157*pow5(m32)*pow5(mU32)
         - 96*pow5(mQ32)*pow5(mU32) - 46196*pow2(mQ32)*pow2(mU32)*pow6(m32) -
         18368*mU32*pow3(mQ32)*pow6(m32) - 18368*mQ32*pow3(mU32)*pow6(m32) -
         902*pow4(mQ32)*pow6(m32) - 902*pow4(mU32)*pow6(m32) - 482*pow2(m32)*
         pow2(mU32)*pow6(mQ32) + 521*mU32*pow3(m32)*pow6(mQ32) + 173*m32*pow3(
         mU32)*pow6(mQ32) - 172*pow4(m32)*pow6(mQ32) - 48*pow4(mU32)*pow6(mQ32)
         - 482*pow2(m32)*pow2(mQ32)*pow6(mU32) + 521*mQ32*pow3(m32)*pow6(mU32)
         + 173*m32*pow3(mQ32)*pow6(mU32) - 172*pow4(m32)*pow6(mU32) - 48*pow4(
         mQ32)*pow6(mU32) + 20650*mU32*pow2(mQ32)*pow7(m32) + 20650*mQ32*pow2(
         mU32)*pow7(m32) + 3158*pow3(mQ32)*pow7(m32) + 3158*pow3(mU32)*pow7(m32
         ) + 6*mU32*pow2(m32)*pow7(mQ32) + 3*m32*pow2(mU32)*pow7(mQ32) - 9*pow3
         (m32)*pow7(mQ32) + 6*mQ32*pow2(m32)*pow7(mU32) + 3*m32*pow2(mQ32)*pow7
         (mU32) - 9*pow3(m32)*pow7(mU32) - 8544*mQ32*mU32*pow8(m32) - 3248*pow2
         (mQ32)*pow8(m32) - 3248*pow2(mU32)*pow8(m32) + 1024*mQ32*pow9(m32) +
         1024*mU32*pow9(m32)))/(3.*pow4(m32 - mQ32)*pow4(m32 - mU32)*pow4(mQ32
         - mU32)))));

   return result;
}

/// 3-loop coefficient O(at*as^2*log^1)
double himalaya::mh2_eft::Mh2EFTCalculator::coeff_as_2_log_1(double mQ32, double mU32, double Xt, double m3, double msq2, double MR2)
{
   using std::log;
   const double m32 = pow2(m3);

   const double result = 
      -653.3333333333334 - 128*zt3 + (1120*(-(mQ32*mU32) + pow2(m32)))
         /(3.*(m32 - mQ32)*(m32 - mU32)) + (160*pow2(-(mQ32*mU32) + pow2(m32)))
         /(3.*pow2(m32 - mQ32)*pow2(m32 - mU32)) + pow2(log(msq2/MR2))*((640*m3
         *msq2*Xt)/((m32 - mQ32)*(m32 - mU32)) - (80*(-4*m32*mQ32*msq2*mU32 - 5
         *mQ32*msq2*pow2(m32) - 8*mQ32*mU32*pow2(m32) - 5*msq2*mU32*pow2(m32) +
         3*m32*msq2*pow2(mQ32) + 4*m32*mU32*pow2(mQ32) + msq2*mU32*pow2(mQ32)
         - 2*pow2(m32)*pow2(mQ32) + 4*m32*mQ32*pow2(mU32) + 3*m32*msq2*pow2(
         mU32) + mQ32*msq2*pow2(mU32) - 2*pow2(m32)*pow2(mU32) - 2*pow2(mQ32)*
         pow2(mU32) + 4*mQ32*pow3(m32) + 6*msq2*pow3(m32) + 4*mU32*pow3(m32) -
         2*pow4(m32)))/(pow2(-m32 + mQ32)*pow2(m32 - mU32))) + log(msq2/MR2)*((
         -5120*m3*msq2*Xt)/((m32 - mQ32)*(m32 - mU32)) + (16*(-480*m32*mQ32*
         msq2*mU32 - 600*mQ32*msq2*pow2(m32) + 444*mQ32*mU32*pow2(m32) - 600*
         msq2*mU32*pow2(m32) + 360*m32*msq2*pow2(mQ32) - 221*m32*mU32*pow2(mQ32
         ) + 120*msq2*mU32*pow2(mQ32) + 111*pow2(m32)*pow2(mQ32) - 221*m32*mQ32
         *pow2(mU32) + 360*m32*msq2*pow2(mU32) + 120*mQ32*msq2*pow2(mU32) + 111
         *pow2(m32)*pow2(mU32) + 110*pow2(mQ32)*pow2(mU32) - 223*mQ32*pow3(m32)
         + 720*msq2*pow3(m32) - 223*mU32*pow3(m32) + 112*pow4(m32)))/(3.*pow2(
         m32 - mQ32)*pow2(m32 - mU32))) - 48*Xt*((64*(m3*fin(mQ32,m32,MR2) - m3
         *fin(mQ32,mU32,MR2) + m3*fin(mU32,m32,MR2) + 3*pow3(m3) + zt2*pow3(m3)
         ))/(9.*(m32 - mQ32)*(m32 - mU32)) - (4*((20*(m3*mQ32 - m3*msq2)*fin(
         mQ32,msq2,MR2))/((mQ32 - mU32)*pow2(m32 - mQ32)) - (2*m3*(mQ32 - mU32)
         *fin(mU32,m32,MR2))/((m32 - mU32)*pow2(m32 - mQ32)) + (2*m3*(mQ32 -
         mU32)*fin(mQ32,m32,MR2))/((m32 - mQ32)*pow2(m32 - mU32)) + (20*(-(m3*
         msq2) + m3*mU32)*fin(mU32,msq2,MR2))/((-mQ32 + mU32)*pow2(m32 - mU32))
         + (2*(6*m3*mQ32 + 60*m3*msq2 + 6*m3*mU32 + m3*mQ32*zt2 + 10*m3*msq2*
         zt2 + m3*mU32*zt2 + 24*pow3(m3)))/((m32 - mQ32)*(m32 - mU32)) - (20*
         fin(m32,msq2,MR2)*(m3*mQ32*msq2 - 2*m3*mQ32*mU32 + m3*msq2*mU32 + mQ32
         *pow3(m3) - 2*msq2*pow3(m3) + mU32*pow3(m3)))/(pow2(m32 - mQ32)*pow2(
         m32 - mU32)) + 3*((-2*m3*fin(mQ32,mU32,MR2))/((m32 - mQ32)*(m32 - mU32
         )) + (2*fin(mQ32,m32,MR2)*(m3*mQ32 - 3*m3*mU32 + 2*pow3(m3)))/((m32 -
         mQ32)*(m32 - mU32)*(mQ32 - mU32)) - (2*fin(mU32,m32,MR2)*(-3*m3*mQ32 +
         m3*mU32 + 2*pow3(m3)))/((m32 - mQ32)*(m32 - mU32)*(mQ32 - mU32)) + (2
         *(11*pow3(m3) + 3*zt2*pow3(m3)))/((m32 - mQ32)*(m32 - mU32))) + (2*fin
         (mQ32,mU32,MR2)*(m3*pow2(mQ32) + m3*pow2(mU32) - 2*mQ32*pow3(m3) - 2*
         mU32*pow3(m3) + 2*pow5(m3)))/(pow2(m32 - mQ32)*pow2(m32 - mU32))))/3.)
         - 48*(-4.111111111111111 - (2*mQ32)/(3.*(-m32 + mQ32)) - (2*mU32)/(3.
         *(-m32 + mU32)) + pow2(-2 + (2*mQ32)/(3.*(-m32 + mQ32)) + (2*mU32)/(3.
         *(-m32 + mU32))) - (16*(-29.625 - (21*mQ32)/(4.*(-m32 + mQ32)) - (56*
         mQ32)/(mQ32 - mU32) + (56*mU32)/(mQ32 - mU32) + (111*mQ32)/(4.*(-m32 +
         mU32)) + (45*mU32)/(2.*(-m32 + mU32)) + zt2/2. - (5*mQ32*zt2)/(2.*(
         -m32 + mQ32)) - (8*mQ32*zt2)/(mQ32 - mU32) + (8*mU32*zt2)/(mQ32 - mU32
         ) + (4*mQ32*zt2)/(-m32 + mU32) + (3*mU32*zt2)/(2.*(-m32 + mU32)) + (2*
         (mQ32 + mU32)*fin(mQ32,mU32,MR2))/((m32 - mQ32)*(m32 - mU32)) + (56*
         pow2(mQ32))/((-m32 + mQ32)*(mQ32 - mU32)) - (111*pow2(mQ32))/(4.*(-m32
         + mQ32)*(-m32 + mU32)) + (8*zt2*pow2(mQ32))/((-m32 + mQ32)*(mQ32 -
         mU32)) - (4*zt2*pow2(mQ32))/((-m32 + mQ32)*(-m32 + mU32)) + (23*pow2(
         mQ32))/(8.*pow2(-m32 + mQ32)) + (3*zt2*pow2(mQ32))/(2.*pow2(-m32 +
         mQ32)) - (56*pow2(mU32))/((mQ32 - mU32)*(-m32 + mU32)) - (8*zt2*pow2(
         mU32))/((mQ32 - mU32)*(-m32 + mU32)) + (23*pow2(mU32))/(8.*pow2(-m32 +
         mU32)) + (3*zt2*pow2(mU32))/(2.*pow2(-m32 + mU32)) + fin(mQ32,m32,MR2
         )*(-5/(2.*(-m32 + mQ32)) - 8/(mQ32 - mU32) + (8*mQ32)/((-m32 + mQ32)*(
         mQ32 - mU32)) + 2/(-m32 + mU32) - (4*mQ32)/((-m32 + mQ32)*(-m32 + mU32
         )) - (3*pow2(mQ32))/pow3(-m32 + mQ32)) - (21*pow3(mQ32))/pow3(-m32 +
         mQ32) - (3*zt2*pow3(mQ32))/pow3(-m32 + mQ32) + fin(mU32,m32,MR2)*(-2/(
         -m32 + mQ32) + 8/(mQ32 - mU32) + 3/(2.*(-m32 + mU32)) - (4*mQ32)/((
         -m32 + mQ32)*(-m32 + mU32)) - (8*mU32)/((mQ32 - mU32)*(-m32 + mU32)) -
         (3*pow2(mU32))/pow3(-m32 + mU32)) - (21*pow3(mU32))/pow3(-m32 + mU32)
         - (3*zt2*pow3(mU32))/pow3(-m32 + mU32)))/9. - (4*(-16.333333333333332
         + (31*mQ32)/(-m32 + mQ32) + (45*msq2)/(-m32 + mQ32) + (9*mU32)/(2.*(
         -m32 + mQ32)) + (9*mQ32)/(2.*(-m32 + mU32)) + (45*msq2)/(-m32 + mU32)
         + (31*mU32)/(-m32 + mU32) - 3*zt2 + (mQ32*zt2)/(4.*(-m32 + mQ32)) + (
         15*msq2*zt2)/(2.*(-m32 + mQ32)) + (3*mU32*zt2)/(4.*(-m32 + mQ32)) + (3
         *mQ32*zt2)/(4.*(-m32 + mU32)) + (15*msq2*zt2)/(2.*(-m32 + mU32)) + (
         mU32*zt2)/(4.*(-m32 + mU32)) - (27275*mQ32*msq2)/(432.*pow2(-m32 +
         mQ32)) - (2825*mQ32*mU32)/(432.*pow2(-m32 + mQ32)) - (65*mQ32*msq2*zt2
         )/(6.*pow2(-m32 + mQ32)) - (mQ32*mU32*zt2)/pow2(-m32 + mQ32) - (23*
         pow2(mQ32))/pow2(-m32 + mQ32) + (1355*msq2*pow2(mQ32))/(432.*(mQ32 -
         msq2)*pow2(-m32 + mQ32)) + (233*mU32*pow2(mQ32))/(432.*(mQ32 - mU32)*
         pow2(-m32 + mQ32)) + (5*msq2*zt2*pow2(mQ32))/(6.*(mQ32 - msq2)*pow2(
         -m32 + mQ32)) - (515*mQ32*pow2(msq2))/(144.*(mQ32 - msq2)*pow2(-m32 +
         mQ32)) - (5*mQ32*zt2*pow2(msq2))/(3.*(mQ32 - msq2)*pow2(-m32 + mQ32))
         + (95*pow2(mQ32)*pow2(msq2))/(216.*pow2(-m32 + mQ32)*pow2(mQ32 - msq2)
         ) + (5*zt2*pow2(mQ32)*pow2(msq2))/(6.*pow2(-m32 + mQ32)*pow2(mQ32 -
         msq2)) - (233*mQ32*pow2(mU32))/(288.*(mQ32 - mU32)*pow2(-m32 + mQ32))
         + (233*pow2(mQ32)*pow2(mU32))/(864.*pow2(-m32 + mQ32)*pow2(mQ32 - mU32
         )) + 3*(16.09722222222222 - (22*mQ32)/(-m32 + mQ32) - (14*mQ32)/(-m32
         + mU32) - (36*mU32)/(-m32 + mU32) + zt2 - (3*mQ32*zt2)/(2.*(-m32 +
         mQ32)) - (2*mQ32*zt2)/(-m32 + mU32) - (7*mU32*zt2)/(2.*(-m32 + mU32))
         - ((mQ32 + mU32)*fin(mQ32,mU32,MR2))/((m32 - mQ32)*(m32 - mU32)) + (14
         *pow2(mQ32))/((-m32 + mQ32)*(-m32 + mU32)) + (2*zt2*pow2(mQ32))/((-m32
         + mQ32)*(-m32 + mU32)) + (fin(mQ32,m32,MR2)*(3*m32*mQ32 - 3*mQ32*mU32
         + pow2(m32) - pow2(mQ32)))/((m32 - mU32)*pow2(-m32 + mQ32)) + (31*
         pow2(mQ32))/pow2(-m32 + mQ32) + (3*zt2*pow2(mQ32))/pow2(-m32 + mQ32) +
         (fin(mU32,m32,MR2)*(3*m32*mU32 - 3*mQ32*mU32 + pow2(m32) - pow2(mU32)
         ))/((m32 - mQ32)*pow2(m32 - mU32)) + (31*pow2(mU32))/pow2(-m32 + mU32)
         + (3*zt2*pow2(mU32))/pow2(-m32 + mU32)) - (6*mQ32*mU32)/pow2(-m32 +
         mU32) - (27275*msq2*mU32)/(432.*pow2(-m32 + mU32)) - (mQ32*mU32*zt2)
         /pow2(-m32 + mU32) - (65*msq2*mU32*zt2)/(6.*pow2(-m32 + mU32)) - (515*
         mU32*pow2(msq2))/(144.*(-msq2 + mU32)*pow2(-m32 + mU32)) - (5*mU32*zt2
         *pow2(msq2))/(3.*(-msq2 + mU32)*pow2(-m32 + mU32)) - (23*pow2(mU32))
         /pow2(-m32 + mU32) + (1355*msq2*pow2(mU32))/(432.*(-msq2 + mU32)*pow2(
         -m32 + mU32)) + (5*msq2*zt2*pow2(mU32))/(6.*(-msq2 + mU32)*pow2(-m32 +
         mU32)) + (95*pow2(msq2)*pow2(mU32))/(216.*pow2(-m32 + mU32)*pow2(
         -msq2 + mU32)) + (5*zt2*pow2(msq2)*pow2(mU32))/(6.*pow2(-m32 + mU32)*
         pow2(-msq2 + mU32)) - (5*fin(mQ32,msq2,MR2)*(-3*m32*mQ32 + 3*m32*msq2
         + mQ32*msq2 - 2*pow2(m32) + pow2(mQ32)))/(2.*pow3(-m32 + mQ32)) + (fin
         (mU32,m32,MR2)*(6*m32*mQ32*mU32 - 8*mQ32*pow2(m32) + 2*mU32*pow2(m32)
         + 3*m32*pow2(mQ32) - 3*m32*pow2(mU32) - mQ32*pow2(mU32) + 2*pow3(m32)
         - pow3(mQ32)))/(4.*(m32 - mU32)*pow3(-m32 + mQ32)) - (95*mQ32*pow3(
         msq2))/(216.*pow2(-m32 + mQ32)*pow2(mQ32 - msq2)) - (5*mQ32*zt2*pow3(
         msq2))/(6.*pow2(-m32 + mQ32)*pow2(mQ32 - msq2)) - (95*mU32*pow3(msq2))
         /(216.*pow2(-m32 + mU32)*pow2(-msq2 + mU32)) - (5*mU32*zt2*pow3(msq2))
         /(6.*pow2(-m32 + mU32)*pow2(-msq2 + mU32)) + (271*pow3(mU32))/(864.*(
         mQ32 - mU32)*pow2(-m32 + mQ32)) + (zt2*pow3(mU32))/(12.*(mQ32 - mU32)*
         pow2(-m32 + mQ32)) - (7*mQ32*pow3(mU32))/(12.*pow2(-m32 + mQ32)*pow2(
         mQ32 - mU32)) - (mQ32*zt2*pow3(mU32))/(12.*pow2(-m32 + mQ32)*pow2(mQ32
         - mU32)) + (fin(mQ32,m32,MR2)*(-6*m32*mQ32*mU32 - 2*mQ32*pow2(m32) +
         8*mU32*pow2(m32) + 3*m32*pow2(mQ32) + mU32*pow2(mQ32) - 3*m32*pow2(
         mU32) - 2*pow3(m32) + pow3(mU32)))/(4.*(m32 - mQ32)*pow3(m32 - mU32))
         + fin(m32,msq2,MR2)*(5/(2.*(-m32 + mQ32)) + 5/(2.*(-m32 + mU32)) + (15
         *mQ32)/(2.*pow2(-m32 + mQ32)) - (15*msq2)/(2.*pow2(-m32 + mQ32)) - (15
         *msq2)/(2.*pow2(-m32 + mU32)) + (15*mU32)/(2.*pow2(-m32 + mU32)) + (10
         *mQ32*msq2)/pow3(-m32 + mQ32) - (10*pow2(mQ32))/pow3(-m32 + mQ32) + (
         10*msq2*mU32)/pow3(-m32 + mU32) - (10*pow2(mU32))/pow3(-m32 + mU32)) -
         (5*fin(mU32,msq2,MR2)*(3*m32*msq2 - 3*m32*mU32 + msq2*mU32 - 2*pow2(
         m32) + pow2(mU32)))/(2.*pow3(-m32 + mU32)) + (271*pow4(mU32))/(864.*
         pow2(-m32 + mQ32)*pow2(mQ32 - mU32)) + (zt2*pow4(mU32))/(12.*pow2(-m32
         + mQ32)*pow2(mQ32 - mU32)) + (fin(mQ32,mU32,MR2)*(15*mU32*pow2(m32)*
         pow2(mQ32) + 15*mQ32*pow2(m32)*pow2(mU32) - 6*m32*pow2(mQ32)*pow2(mU32
         ) - 20*mQ32*mU32*pow3(m32) + 14*pow2(mQ32)*pow3(m32) + 14*pow2(mU32)*
         pow3(m32) - 6*m32*mU32*pow3(mQ32) - 11*pow2(m32)*pow3(mQ32) + pow2(
         mU32)*pow3(mQ32) - 6*m32*mQ32*pow3(mU32) - 11*pow2(m32)*pow3(mU32) +
         pow2(mQ32)*pow3(mU32) - 6*mQ32*pow4(m32) - 6*mU32*pow4(m32) + 3*m32*
         pow4(mQ32) + mU32*pow4(mQ32) + 3*m32*pow4(mU32) + mQ32*pow4(mU32) + 4*
         pow5(m32)))/(4.*pow3(-m32 + mQ32)*pow3(m32 - mU32))))/3.) + log(
         mU32/MR2)*(log(msq2/MR2)*((-64*m3*(m32 + 60*msq2 - mU32)*mU32*Xt)/(3.*
         (-mQ32 + mU32)*pow2(m32 - mU32)) + (16*(90*m32*msq2*mU32 + 2*mU32*pow2
         (m32) - 3*m32*pow2(mU32) + 30*msq2*pow2(mU32) + pow3(mU32)))/(3.*pow3(
         m32 - mU32))) + (32*(-180*msq2*mU32*pow2(m32)*pow2(mQ32) + 210*mQ32*
         msq2*pow2(m32)*pow2(mU32) - 150*m32*msq2*pow2(mQ32)*pow2(mU32) + 257*
         pow2(m32)*pow2(mQ32)*pow2(mU32) + 90*mQ32*msq2*mU32*pow3(m32) + 128*
         mU32*pow2(mQ32)*pow3(m32) + 76*mQ32*pow2(mU32)*pow3(m32) - 90*msq2*
         pow2(mU32)*pow3(m32) + 90*m32*msq2*mU32*pow3(mQ32) - 54*mU32*pow2(m32)
         *pow3(mQ32) - 182*m32*pow2(mU32)*pow3(mQ32) + 30*msq2*pow2(mU32)*pow3(
         mQ32) + 27*pow3(m32)*pow3(mQ32) + 60*m32*mQ32*msq2*pow3(mU32) - 461*
         mQ32*pow2(m32)*pow3(mU32) - 30*msq2*pow2(m32)*pow3(mU32) + 107*m32*
         pow2(mQ32)*pow3(mU32) - 30*msq2*pow2(mQ32)*pow3(mU32) + 153*pow3(m32)*
         pow3(mU32) + 57*pow3(mQ32)*pow3(mU32) - 144*mQ32*mU32*pow4(m32) - 54*
         pow2(mQ32)*pow4(m32) - 58*pow2(mU32)*pow4(m32) + 9*m32*mU32*pow4(mQ32)
         + 3*pow2(mU32)*pow4(mQ32) + 148*m32*mQ32*pow4(mU32) + 2*pow2(m32)*
         pow4(mU32) - 54*pow2(mQ32)*pow4(mU32) + 27*mQ32*pow5(m32) + 37*mU32*
         pow5(m32) - 18*m32*pow5(mU32) - 6*mQ32*pow5(mU32)))/(3.*(-mQ32 + mU32)
         *pow2(m32 - mQ32)*pow3(m32 - mU32)) - (256*Xt*(-30*m3*mQ32*msq2*mU32 -
         3*m3*mU32*pow2(mQ32) + 32*m3*mQ32*pow2(mU32) + 27*mQ32*mU32*pow3(m3)
         + 30*msq2*mU32*pow3(m3) - 40*pow2(mU32)*pow3(m3) + 6*m3*pow3(mU32) - 8
         *mQ32*pow5(m3) - 22*mU32*pow5(m3) + 8*pow7(m3)))/(3.*(m32 - mQ32)*(
         mQ32 - mU32)*pow2(m32 - mU32))) + pow2(log(mQ32/MR2))*((-32*(2*m32*
         mQ32 - pow2(mQ32)))/(3.*pow2(m32 - mQ32)) + (512*m32*pow2(mQ32)*pow2(
         Xt))/(pow2(m32 - mQ32)*pow2(mQ32 - mU32)) + (160*pow2(-2*m32*mQ32 +
         pow2(mQ32)))/(3.*pow4(m32 - mQ32)) - 48*((-2*m32*mQ32 + pow2(mQ32))/(
         9.*pow2(m32 - mQ32)) + (4*pow2(2*m32*mQ32 - pow2(mQ32)))/(9.*pow4(m32
         - mQ32)) - (16*((-25*mQ32)/(4.*(-m32 + mQ32)) - (4*mQ32)/(mQ32 - mU32)
         + (4*pow2(mQ32))/((-m32 + mQ32)*(mQ32 - mU32)) + (15*pow2(mQ32))/(2.*
         pow2(-m32 + mQ32)) - (6*pow3(mQ32))/pow3(-m32 + mQ32) + (3*pow4(mQ32))
         /(8.*pow4(-m32 + mQ32))))/9. - (4*((3*(-6*mQ32*pow2(m32) + 5*m32*pow2(
         mQ32) - 7*pow3(mQ32)))/(4.*pow3(m32 - mQ32)) + (-60*mQ32*msq2*mU32*
         pow2(m32) - 20*m32*msq2*mU32*pow2(mQ32) + 10*msq2*pow2(m32)*pow2(mQ32)
         + 42*mU32*pow2(m32)*pow2(mQ32) + 30*m32*mQ32*msq2*pow2(mU32) - 31*
         mQ32*pow2(m32)*pow2(mU32) - 21*m32*pow2(mQ32)*pow2(mU32) + 10*msq2*
         pow2(mQ32)*pow2(mU32) + 30*mQ32*msq2*pow3(m32) + 52*mQ32*mU32*pow3(m32
         ) - 10*pow2(mQ32)*pow3(m32) + 2*pow2(mU32)*pow3(m32) - 15*m32*mU32*
         pow3(mQ32) - 3*pow2(m32)*pow3(mQ32) + 6*pow2(mU32)*pow3(mQ32) + 3*m32*
         mQ32*pow3(mU32) + pow2(mQ32)*pow3(mU32) - 28*mQ32*pow4(m32) - 4*mU32*
         pow4(m32) + 3*m32*pow4(mQ32) + mU32*pow4(mQ32) + 2*pow5(m32))/(8.*pow2
         (m32 - mU32)*pow3(m32 - mQ32))))/3.) - (64*Xt*(-30*m3*msq2*pow2(mQ32)*
         pow2(mU32) + 30*mQ32*msq2*pow2(mU32)*pow3(m3) + 5*pow2(mQ32)*pow2(mU32
         )*pow3(m3) + 30*m3*msq2*mU32*pow3(mQ32) + 20*m3*pow2(mU32)*pow3(mQ32)
         - 30*msq2*pow3(m3)*pow3(mQ32) - 2*mU32*pow3(m3)*pow3(mQ32) - 3*m3*pow2
         (mQ32)*pow3(mU32) + 3*mQ32*pow3(m3)*pow3(mU32) - 30*m3*mU32*pow4(mQ32)
         + 42*pow3(m3)*pow4(mQ32) - 30*mQ32*msq2*mU32*pow5(m3) + 30*msq2*pow2(
         mQ32)*pow5(m3) - 3*mU32*pow2(mQ32)*pow5(m3) - 9*mQ32*pow2(mU32)*pow5(
         m3) - 36*pow3(mQ32)*pow5(m3) - 3*m3*pow5(mQ32) + 3*mQ32*mU32*pow7(m3)
         + 13*pow2(mQ32)*pow7(m3)))/(3.*(m32 - mU32)*pow2(mQ32 - mU32)*pow3(m32
         - mQ32))) + pow2(log(mU32/MR2))*((-16*(2*m32*mU32 - pow2(mU32)))/(3.*
         pow2(m32 - mU32)) + (512*m32*pow2(mU32)*pow2(Xt))/(pow2(m32 - mU32)*
         pow2(-mQ32 + mU32)) + (160*pow2(-2*m32*mU32 + pow2(mU32)))/(3.*pow4(
         m32 - mU32)) - 48*((-2*m32*mU32 + pow2(mU32))/(9.*pow2(m32 - mU32)) +
         (4*pow2(2*m32*mU32 - pow2(mU32)))/(9.*pow4(m32 - mU32)) - (16*((4*mU32
         )/(mQ32 - mU32) - (25*mU32)/(4.*(-m32 + mU32)) - (4*pow2(mU32))/((mQ32
         - mU32)*(-m32 + mU32)) + (15*pow2(mU32))/(2.*pow2(-m32 + mU32)) - (6*
         pow3(mU32))/pow3(-m32 + mU32) + (3*pow4(mU32))/(8.*pow4(-m32 + mU32)))
         )/9. - (4*((3*(-6*mU32*pow2(m32) + 5*m32*pow2(mU32) - 7*pow3(mU32)))/(
         4.*pow3(m32 - mU32)) + (-60*mQ32*msq2*mU32*pow2(m32) + 30*m32*msq2*
         mU32*pow2(mQ32) - 31*mU32*pow2(m32)*pow2(mQ32) - 20*m32*mQ32*msq2*pow2
         (mU32) + 42*mQ32*pow2(m32)*pow2(mU32) + 10*msq2*pow2(m32)*pow2(mU32) -
         21*m32*pow2(mQ32)*pow2(mU32) + 10*msq2*pow2(mQ32)*pow2(mU32) + 52*
         mQ32*mU32*pow3(m32) + 30*msq2*mU32*pow3(m32) + 2*pow2(mQ32)*pow3(m32)
         - 10*pow2(mU32)*pow3(m32) + 3*m32*mU32*pow3(mQ32) + pow2(mU32)*pow3(
         mQ32) - 15*m32*mQ32*pow3(mU32) - 3*pow2(m32)*pow3(mU32) + 6*pow2(mQ32)
         *pow3(mU32) - 4*mQ32*pow4(m32) - 28*mU32*pow4(m32) + 3*m32*pow4(mU32)
         + mQ32*pow4(mU32) + 2*pow5(m32))/(8.*pow2(-m32 + mQ32)*pow3(m32 - mU32
         ))))/3.) - (64*Xt*(-30*m3*msq2*pow2(mQ32)*pow2(mU32) + 30*msq2*mU32*
         pow2(mQ32)*pow3(m3) + 3*pow2(mQ32)*pow2(mU32)*pow3(m3) - 3*m3*pow2(
         mU32)*pow3(mQ32) + 3*mU32*pow3(m3)*pow3(mQ32) + 30*m3*mQ32*msq2*pow3(
         mU32) + 21*m3*pow2(mQ32)*pow3(mU32) - mQ32*pow3(m3)*pow3(mU32) - 30*
         msq2*pow3(m3)*pow3(mU32) - 31*m3*mQ32*pow4(mU32) + 43*pow3(m3)*pow4(
         mU32) - 30*mQ32*msq2*mU32*pow5(m3) - 8*mU32*pow2(mQ32)*pow5(m3) - 2*
         mQ32*pow2(mU32)*pow5(m3) + 30*msq2*pow2(mU32)*pow5(m3) - 38*pow3(mU32)
         *pow5(m3) - 3*m3*pow5(mU32) + 2*mQ32*mU32*pow7(m3) + 14*pow2(mU32)*
         pow7(m3)))/(3.*(m32 - mQ32)*pow2(-mQ32 + mU32)*pow3(m32 - mU32))) +
         log(mQ32/MR2)*(log(msq2/MR2)*((-64*m3*mQ32*(m32 - mQ32 + 60*msq2)*Xt)/
         (3.*(mQ32 - mU32)*pow2(m32 - mQ32)) + (16*(90*m32*mQ32*msq2 + 2*mQ32*
         pow2(m32) - 3*m32*pow2(mQ32) + 30*msq2*pow2(mQ32) + pow3(mQ32)))/(3.*
         pow3(m32 - mQ32))) + (16*(420*msq2*mU32*pow2(m32)*pow2(mQ32) - 360*
         mQ32*msq2*pow2(m32)*pow2(mU32) - 300*m32*msq2*pow2(mQ32)*pow2(mU32) +
         481*pow2(m32)*pow2(mQ32)*pow2(mU32) + 180*mQ32*msq2*mU32*pow3(m32) -
         180*msq2*pow2(mQ32)*pow3(m32) + 119*mU32*pow2(mQ32)*pow3(m32) + 312*
         mQ32*pow2(mU32)*pow3(m32) + 120*m32*msq2*mU32*pow3(mQ32) - 60*msq2*
         pow2(m32)*pow3(mQ32) - 868*mU32*pow2(m32)*pow3(mQ32) + 204*m32*pow2(
         mU32)*pow3(mQ32) - 60*msq2*pow2(mU32)*pow3(mQ32) + 272*pow3(m32)*pow3(
         mQ32) + 180*m32*mQ32*msq2*pow3(mU32) - 140*mQ32*pow2(m32)*pow3(mU32) -
         333*m32*pow2(mQ32)*pow3(mU32) + 60*msq2*pow2(mQ32)*pow3(mU32) + 65*
         pow3(m32)*pow3(mU32) + 104*pow3(mQ32)*pow3(mU32) - 300*mQ32*mU32*pow4(
         m32) - 81*pow2(mQ32)*pow4(m32) - 131*pow2(mU32)*pow4(m32) + 275*m32*
         mU32*pow4(mQ32) + 15*pow2(m32)*pow4(mQ32) - 98*pow2(mU32)*pow4(mQ32) +
         18*m32*mQ32*pow4(mU32) + 6*pow2(mQ32)*pow4(mU32) + 62*mQ32*pow5(m32)
         + 66*mU32*pow5(m32) - 36*m32*pow5(mQ32) - 12*mU32*pow5(mQ32)))/(3.*(
         mQ32 - mU32)*pow2(m32 - mU32)*pow3(m32 - mQ32)) + (256*Xt*(-30*m3*mQ32
         *msq2*mU32 + 32*m3*mU32*pow2(mQ32) - 3*m3*mQ32*pow2(mU32) + 30*mQ32*
         msq2*pow3(m3) + 27*mQ32*mU32*pow3(m3) - 40*pow2(mQ32)*pow3(m3) + 6*m3*
         pow3(mQ32) - 22*mQ32*pow5(m3) - 8*mU32*pow5(m3) + 8*pow7(m3)))/(3.*(
         m32 - mU32)*(mQ32 - mU32)*pow2(m32 - mQ32)) + log(mU32/MR2)*((-1024*
         m32*mQ32*mU32*pow2(Xt))/((m32 - mQ32)*(m32 - mU32)*pow2(mQ32 - mU32))
         - (16*(-119*pow2(m32)*pow2(mQ32)*pow2(mU32) + 118*mU32*pow2(mQ32)*pow3
         (m32) + 121*mQ32*pow2(mU32)*pow3(m32) - 57*mU32*pow2(m32)*pow3(mQ32) +
         34*m32*pow2(mU32)*pow3(mQ32) - 58*mQ32*pow2(m32)*pow3(mU32) + 34*m32*
         pow2(mQ32)*pow3(mU32) + pow3(m32)*pow3(mU32) - 13*pow3(mQ32)*pow3(mU32
         ) - 84*mQ32*mU32*pow4(m32) - 3*pow2(mU32)*pow4(m32) + 9*m32*mU32*pow4(
         mQ32) + 3*pow2(mU32)*pow4(mQ32) + 9*m32*mQ32*pow4(mU32) + 3*pow2(mQ32)
         *pow4(mU32) + 2*mU32*pow5(m32)))/(3.*pow3(m32 - mQ32)*pow3(m32 - mU32)
         ) + (64*Xt*(55*pow2(mQ32)*pow2(mU32)*pow3(m3) - 21*m3*pow2(mU32)*pow3(
         mQ32) + 19*mU32*pow3(m3)*pow3(mQ32) - 23*m3*pow2(mQ32)*pow3(mU32) + 22
         *mQ32*pow3(m3)*pow3(mU32) + 6*m3*mU32*pow4(mQ32) + 6*m3*mQ32*pow4(mU32
         ) - 46*mU32*pow2(mQ32)*pow5(m3) - 49*mQ32*pow2(mU32)*pow5(m3) - pow3(
         mU32)*pow5(m3) + 31*mQ32*mU32*pow7(m3) + pow2(mU32)*pow7(m3)))/(3.*
         pow2(m32 - mQ32)*pow2(m32 - mU32)*pow2(mQ32 - mU32)))) + pow2(log(
         m32/MR2))*((512*pow2(Xt)*pow3(m32))/(pow2(m32 - mQ32)*pow2(m32 - mU32)
         ) - (8*(-360*msq2*pow2(mQ32)*pow2(mU32)*pow3(m32) + 120*msq2*pow2(m32)
         *pow2(mU32)*pow3(mQ32) + 240*msq2*mU32*pow3(m32)*pow3(mQ32) - 2092*
         pow2(mU32)*pow3(m32)*pow3(mQ32) + 120*msq2*pow2(m32)*pow2(mQ32)*pow3(
         mU32) + 240*mQ32*msq2*pow3(m32)*pow3(mU32) - 2092*pow2(mQ32)*pow3(m32)
         *pow3(mU32) - 576*pow2(m32)*pow3(mQ32)*pow3(mU32) - 240*msq2*mU32*pow2
         (mQ32)*pow4(m32) - 240*mQ32*msq2*pow2(mU32)*pow4(m32) + 8320*pow2(mQ32
         )*pow2(mU32)*pow4(m32) - 360*msq2*pow3(mQ32)*pow4(m32) + 3052*mU32*
         pow3(mQ32)*pow4(m32) + 3052*mQ32*pow3(mU32)*pow4(m32) - 360*msq2*pow3(
         mU32)*pow4(m32) - 60*msq2*mU32*pow2(m32)*pow4(mQ32) - 30*m32*msq2*pow2
         (mU32)*pow4(mQ32) + 220*pow2(m32)*pow2(mU32)*pow4(mQ32) + 90*msq2*pow3
         (m32)*pow4(mQ32) - 451*mU32*pow3(m32)*pow4(mQ32) + 261*m32*pow3(mU32)*
         pow4(mQ32) + 18*pow4(m32)*pow4(mQ32) - 60*mQ32*msq2*pow2(m32)*pow4(
         mU32) - 30*m32*msq2*pow2(mQ32)*pow4(mU32) + 220*pow2(m32)*pow2(mQ32)*
         pow4(mU32) - 451*mQ32*pow3(m32)*pow4(mU32) + 90*msq2*pow3(m32)*pow4(
         mU32) + 261*m32*pow3(mQ32)*pow4(mU32) + 18*pow4(m32)*pow4(mU32) - 72*
         pow4(mQ32)*pow4(mU32) + 480*mQ32*msq2*mU32*pow5(m32) + 510*msq2*pow2(
         mQ32)*pow5(m32) - 8565*mU32*pow2(mQ32)*pow5(m32) - 8565*mQ32*pow2(mU32
         )*pow5(m32) + 510*msq2*pow2(mU32)*pow5(m32) - 549*pow3(mQ32)*pow5(m32)
         - 549*pow3(mU32)*pow5(m32) - 6*mU32*pow2(m32)*pow5(mQ32) - 3*m32*pow2
         (mU32)*pow5(mQ32) + 9*pow3(m32)*pow5(mQ32) - 6*mQ32*pow2(m32)*pow5(
         mU32) - 3*m32*pow2(mQ32)*pow5(mU32) + 9*pow3(m32)*pow5(mU32) - 420*
         mQ32*msq2*pow6(m32) + 8052*mQ32*mU32*pow6(m32) - 420*msq2*mU32*pow6(
         m32) + 1976*pow2(mQ32)*pow6(m32) + 1976*pow2(mU32)*pow6(m32) - 1986*
         mQ32*pow7(m32) + 180*msq2*pow7(m32) - 1986*mU32*pow7(m32) + 508*pow8(
         m32)))/(3.*pow4(m32 - mQ32)*pow4(m32 - mU32)) + (64*Xt*(22*pow11(m3) -
         30*msq2*mU32*pow2(mQ32)*pow3(m3) - 30*mQ32*msq2*pow2(mU32)*pow3(m3) +
         248*pow2(mQ32)*pow2(mU32)*pow3(m3) - 3*mU32*pow3(m3)*pow3(mQ32) - 3*
         mQ32*pow3(m3)*pow3(mU32) + 120*mQ32*msq2*mU32*pow5(m3) + 30*msq2*pow2(
         mQ32)*pow5(m3) - 393*mU32*pow2(mQ32)*pow5(m3) - 393*mQ32*pow2(mU32)*
         pow5(m3) + 30*msq2*pow2(mU32)*pow5(m3) + 3*pow3(mQ32)*pow5(m3) + 3*
         pow3(mU32)*pow5(m3) - 90*mQ32*msq2*pow7(m3) + 598*mQ32*mU32*pow7(m3) -
         90*msq2*mU32*pow7(m3) + 129*pow2(mQ32)*pow7(m3) + 129*pow2(mU32)*pow7
         (m3) - 170*mQ32*pow9(m3) + 60*msq2*pow9(m3) - 170*mU32*pow9(m3)))/(3.*
         pow3(m32 - mQ32)*pow3(m32 - mU32))) + log(m32/MR2)*((32*(90*msq2*mU32*
         pow2(m32)*pow2(mQ32) + 90*mQ32*msq2*pow2(m32)*pow2(mU32) + 1322*pow2(
         m32)*pow2(mQ32)*pow2(mU32) - 180*mQ32*msq2*mU32*pow3(m32) + 270*msq2*
         pow2(mQ32)*pow3(m32) - 2690*mU32*pow2(mQ32)*pow3(m32) - 2690*mQ32*pow2
         (mU32)*pow3(m32) + 270*msq2*pow2(mU32)*pow3(m32) - 30*m32*msq2*mU32*
         pow3(mQ32) - 90*msq2*pow2(m32)*pow3(mQ32) + 573*mU32*pow2(m32)*pow3(
         mQ32) - 207*m32*pow2(mU32)*pow3(mQ32) - 94*pow3(m32)*pow3(mQ32) - 30*
         m32*mQ32*msq2*pow3(mU32) + 573*mQ32*pow2(m32)*pow3(mU32) - 90*msq2*
         pow2(m32)*pow3(mU32) - 207*m32*pow2(mQ32)*pow3(mU32) - 94*pow3(m32)*
         pow3(mU32) + 24*pow3(mQ32)*pow3(mU32) - 240*mQ32*msq2*pow4(m32) + 4326
         *mQ32*mU32*pow4(m32) - 240*msq2*mU32*pow4(m32) + 759*pow2(mQ32)*pow4(
         m32) + 759*pow2(mU32)*pow4(m32) - 3*m32*mU32*pow4(mQ32) - 9*pow2(m32)*
         pow4(mQ32) - 3*m32*mQ32*pow4(mU32) - 9*pow2(m32)*pow4(mU32) - 1414*
         mQ32*pow5(m32) + 180*msq2*pow5(m32) - 1414*mU32*pow5(m32) + 498*pow6(
         m32)))/(3.*pow3(m32 - mQ32)*pow3(m32 - mU32)) - (256*Xt*(-30*mQ32*msq2
         *pow3(m3) + 146*mQ32*mU32*pow3(m3) - 30*msq2*mU32*pow3(m3) - 3*pow2(
         mQ32)*pow3(m3) - 3*pow2(mU32)*pow3(m3) - 95*mQ32*pow5(m3) + 60*msq2*
         pow5(m3) - 95*mU32*pow5(m3) + 50*pow7(m3)))/(3.*pow2(m32 - mQ32)*pow2(
         m32 - mU32)) + log(msq2/MR2)*((48*(pow2(m32)*pow2(mQ32) + pow2(m32)*
         pow2(mU32) - 2*mQ32*pow3(m32) - 2*mU32*pow3(m32) + 2*pow4(m32)))/(pow2
         (m32 - mQ32)*pow2(m32 - mU32)) - 48*((10*(pow2(m32)*pow2(mQ32) + pow2(
         m32)*pow2(mU32) - 2*mQ32*pow3(m32) - 2*mU32*pow3(m32) + 2*pow4(m32)))/
         (9.*pow2(-m32 + mQ32)*pow2(m32 - mU32)) + (10*msq2*(-3*mU32*pow2(m32)*
         pow2(mQ32) - 3*mQ32*pow2(m32)*pow2(mU32) + 6*mQ32*mU32*pow3(m32) - 9*
         pow2(mQ32)*pow3(m32) - 9*pow2(mU32)*pow3(m32) + m32*mU32*pow3(mQ32) +
         3*pow2(m32)*pow3(mQ32) + m32*mQ32*pow3(mU32) + 3*pow2(m32)*pow3(mU32)
         + 8*mQ32*pow4(m32) + 8*mU32*pow4(m32) - 6*pow5(m32)))/(3.*pow3(-m32 +
         mQ32)*pow3(m32 - mU32))) + (64*Xt*(-60*mQ32*msq2*pow3(m3) + mQ32*mU32*
         pow3(m3) - 60*msq2*mU32*pow3(m3) - mQ32*pow5(m3) + 120*msq2*pow5(m3) -
         mU32*pow5(m3) + pow7(m3)))/(3.*pow2(m32 - mQ32)*pow2(m32 - mU32))) +
         log(mU32/MR2)*((1024*mU32*pow2(m32)*pow2(Xt))/((m32 - mQ32)*(mQ32 -
         mU32)*pow2(m32 - mU32)) + (16*(318*pow2(mQ32)*pow2(mU32)*pow3(m32) -
         110*pow2(m32)*pow2(mU32)*pow3(mQ32) + 221*mU32*pow3(m32)*pow3(mQ32) +
         81*pow2(m32)*pow2(mQ32)*pow3(mU32) - 47*mQ32*pow3(m32)*pow3(mU32) - 23
         *m32*pow3(mQ32)*pow3(mU32) - 659*mU32*pow2(mQ32)*pow4(m32) - 374*mQ32*
         pow2(mU32)*pow4(m32) - 64*pow3(mQ32)*pow4(m32) - 83*pow3(mU32)*pow4(
         m32) + 4*mQ32*pow2(m32)*pow4(mU32) - 4*m32*pow2(mQ32)*pow4(mU32) + 48*
         pow3(m32)*pow4(mU32) + 684*mQ32*mU32*pow5(m32) + 192*pow2(mQ32)*pow5(
         m32) + 214*pow2(mU32)*pow5(m32) - 3*m32*mQ32*pow5(mU32) - 9*pow2(m32)*
         pow5(mU32) - 192*mQ32*pow6(m32) - 258*mU32*pow6(m32) + 64*pow7(m32)))/
         (3.*pow3(m32 - mQ32)*pow4(m32 - mU32)) + (128*Xt*(16*pow11(m3) + 53*
         pow2(mQ32)*pow2(mU32)*pow3(m3) - 3*pow3(m3)*pow4(mU32) - 85*mU32*pow2(
         mQ32)*pow5(m3) - 98*mQ32*pow2(mU32)*pow5(m3) + pow3(mU32)*pow5(m3) +
         162*mQ32*mU32*pow7(m3) + 16*pow2(mQ32)*pow7(m3) + 52*pow2(mU32)*pow7(
         m3) - 32*mQ32*pow9(m3) - 82*mU32*pow9(m3)))/(3.*(mQ32 - mU32)*pow2(m32
         - mQ32)*pow3(m32 - mU32))) + log(mQ32/MR2)*((-1024*mQ32*pow2(m32)*
         pow2(Xt))/((m32 - mU32)*(mQ32 - mU32)*pow2(m32 - mQ32)) + (16*(321*
         pow2(mQ32)*pow2(mU32)*pow3(m32) + 81*pow2(m32)*pow2(mU32)*pow3(mQ32) -
         43*mU32*pow3(m32)*pow3(mQ32) - 111*pow2(m32)*pow2(mQ32)*pow3(mU32) +
         223*mQ32*pow3(m32)*pow3(mU32) - 23*m32*pow3(mQ32)*pow3(mU32) - 383*
         mU32*pow2(mQ32)*pow4(m32) - 665*mQ32*pow2(mU32)*pow4(m32) - 87*pow3(
         mQ32)*pow4(m32) - 65*pow3(mU32)*pow4(m32) + 3*mU32*pow2(m32)*pow4(mQ32
         ) - 4*m32*pow2(mU32)*pow4(mQ32) + 49*pow3(m32)*pow4(mQ32) + 694*mQ32*
         mU32*pow5(m32) + 221*pow2(mQ32)*pow5(m32) + 195*pow2(mU32)*pow5(m32) -
         3*m32*mU32*pow5(mQ32) - 9*pow2(m32)*pow5(mQ32) - 264*mQ32*pow6(m32) -
         196*mU32*pow6(m32) + 66*pow7(m32)))/(3.*pow3(m32 - mU32)*pow4(m32 -
         mQ32)) - (64*Xt*(32*pow11(m3) + 107*pow2(mQ32)*pow2(mU32)*pow3(m3) -
         mU32*pow3(m3)*pow3(mQ32) - 6*pow3(m3)*pow4(mQ32) - 195*mU32*pow2(mQ32)
         *pow5(m3) - 172*mQ32*pow2(mU32)*pow5(m3) + 3*pow3(mQ32)*pow5(m3) + 325
         *mQ32*mU32*pow7(m3) + 102*pow2(mQ32)*pow7(m3) + 33*pow2(mU32)*pow7(m3)
         - 163*mQ32*pow9(m3) - 65*mU32*pow9(m3)))/(3.*(mQ32 - mU32)*pow2(m32 -
         mU32)*pow3(m32 - mQ32))));

   return result;
}

/// 3-loop coefficient O(at*as^2*log^2)
double himalaya::mh2_eft::Mh2EFTCalculator::coeff_as_2_log_2(double mQ32, double mU32, double Xt, double m3, double msq2, double MR2)
{
   using std::log;
   const double m32 = pow2(m3);

   const double result = 
      144*log(msq2/MR2) - (64*(m32*mQ32 + m32*mU32 - 5*mQ32*mU32 + 3*
         pow2(m32)))/((m32 - mQ32)*(m32 - mU32)) + log(mQ32/MR2)*((1024*m3*mQ32
         *Xt)/((m32 - mQ32)*(mQ32 - mU32)) + (32*(-18*m32*mQ32 + pow2(m32) + 9*
         pow2(mQ32)))/pow2(m32 - mQ32)) + log(mU32/MR2)*((1024*m3*mU32*Xt)/((
         m32 - mU32)*(-mQ32 + mU32)) + (16*(-34*m32*mU32 + pow2(m32) + 17*pow2(
         mU32)))/pow2(m32 - mU32)) + log(m32/MR2)*(192 - (1024*Xt*pow3(m3))/((
         m32 - mQ32)*(m32 - mU32)) + (256*(pow2(m32)*pow2(mQ32) + pow2(m32)*
         pow2(mU32) - 2*mQ32*pow3(m32) - 2*mU32*pow3(m32) + 2*pow4(m32)))/(pow2
         (m32 - mQ32)*pow2(m32 - mU32)));

   return result;
}

/// 3-loop coefficient O(at*as^2*log^3)
double himalaya::mh2_eft::Mh2EFTCalculator::coeff_as_2_log_3()
{
   return 736;
}

double himalaya::mh2_eft::Mh2EFTCalculator::Mh2_EFT_1loop(
   double yt, double tb, double mt, double mQ32, double mU32,
   double Xt, double MR2)
{
   const double mt2 = pow2(mt);
   const double L = std::log(MR2/mt2);
   const double at = pow2(yt)/(4*Pi);
   const double sb2 = pow2(tb)/(1 + pow2(tb));
   const double oneLoop = 1./(4*Pi);

   return oneLoop * mt2 * sb2 * at * (
      coeff_as_0_log_0(mQ32, mU32, Xt, MR2) +
      coeff_as_0_log_1() * L
   );
}

double himalaya::mh2_eft::Mh2EFTCalculator::Mh2_EFT_2loop(
   double yt, double tb, double mt, double mQ32, double mU32,
   double Xt, double MR2, double g3, double m3)
{
   const double mt2 = pow2(mt);
   const double L = std::log(MR2/mt2);
   const double at = pow2(yt)/(4*Pi);
   const double as = pow2(g3)/(4*Pi);
   const double sb2 = pow2(tb)/(1 + pow2(tb));
   const double twoLoop = 1./pow2(4*Pi);

   return twoLoop * mt2 * sb2 * at * as * (
      coeff_as_1_log_0(mQ32, mU32, Xt, m3, MR2) +
      coeff_as_1_log_1(mQ32, mU32, Xt, m3, MR2) * L +
      coeff_as_1_log_2() * pow2(L)
   );
}

double himalaya::mh2_eft::Mh2EFTCalculator::Mh2_EFT_3loop(
   double yt, double tb, double mt, double mQ32, double mU32,
   double Xt, double MR2, double g3, double m3, double msq2)
{
   const double mt2 = pow2(mt);
   const double L = std::log(MR2/mt2);
   const double at = pow2(yt)/(4*Pi);
   const double as = pow2(g3)/(4*Pi);
   const double sb2 = pow2(tb)/(1 + pow2(tb));
   const double threeLoop = 1./pow3(4*Pi);

   return threeLoop * mt2 * sb2 * at * pow2(as) * (
      coeff_as_2_log_0(mQ32, mU32, Xt, m3, msq2, MR2) +
      coeff_as_2_log_1(mQ32, mU32, Xt, m3, msq2, MR2) * L +
      coeff_as_2_log_2(mQ32, mU32, Xt, m3, msq2, MR2) * pow2(L) +
      coeff_as_2_log_3() * pow3(L)
   );
}

/**
 * 	A function which maps a boolean to a string.
 * 	@param tf a boolean.
 * 	@return A string which is 'true' if tf is true or 'false' if tf is false.
 */
std::string himalaya::mh2_eft::Mh2EFTCalculator::tf(const bool tf){
   return tf ? "true" : "false";
}

/**
 * 	Function to check terms
 */
void himalaya::mh2_eft::Mh2EFTCalculator::checkTerms(double mQ32, double mU32, double Xt, double MR2, double m3, double msq2){
   bool checkPassed = false;
   double dum = coeff_as_0_log_0(mQ32, mU32, Xt, MR2);
   checkPassed = dum - 182.7043188198589 < 1e-03;
   std::cout << "Check as^0 log^0 passed: " << tf(checkPassed) << ".\n";
   dum = coeff_as_0_log_1();
   checkPassed = dum - 12. < 1e-03;
   std::cout << "Check as^0 log^1 passed: " << tf(checkPassed) << ".\n";
   dum = coeff_as_1_log_0(mQ32, mU32, Xt, m3, MR2);
   checkPassed = dum - 136.3533935546875 < 1e-03;
   std::cout << "Check as^1 log^0 passed: " << tf(checkPassed) << ".\n";
   dum = coeff_as_1_log_1(mQ32, mU32, Xt, m3, MR2);
   checkPassed = dum - (-126.2886909716763) < 1e-03;
   std::cout << "Check as^1 log^1 passed: " << tf(checkPassed) << ".\n";
   dum = coeff_as_1_log_2();
   checkPassed = dum - 96. < 1e-03;
   std::cout << "Check as^1 log^2 passed: " << tf(checkPassed) << ".\n";
   dum = coeff_as_2_log_0(mQ32, mU32, Xt, m3, msq2, MR2);
   checkPassed = dum - (-9.11632e17) < 1e-03;
   std::cout << "Check as^2 log^0 passed: " << tf(checkPassed) << ".\n";
   std::cout << "is: "<< dum << " should be: " << -9.11632e17 << "\n";
   dum = coeff_as_2_log_1(mQ32, mU32, Xt, m3, msq2, MR2);
   checkPassed = dum - (-1.12658e7) < 1e-03;
   std::cout << "Check as^2 log^1 passed: " << tf(checkPassed) << ".\n";
   std::cout << "is: "<< dum << " should be: " << -1.12658e7 << "\n";
   dum = coeff_as_2_log_2(mQ32, mU32, Xt, m3, msq2, MR2);
   checkPassed = dum - (-679.9017579145730) < 1e-03;
   std::cout << "Check as^2 log^2 passed: " << tf(checkPassed) << ".\n";
   dum = coeff_as_2_log_3();
   checkPassed = dum - 736. < 1e-03;
   std::cout << "Check as^2 log^3 passed: " << tf(checkPassed) << ".\n";
}
