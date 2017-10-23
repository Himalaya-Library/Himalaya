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
      -64 + (32*m32)/mQ32 + (32*m32)/mU32 - (64*m32*pow2(Xt))/(mQ32*
         mU32) - (512*m3*pow3(Xt))/pow2(mQ32 - mU32) + (192/pow2(mQ32 - mU32) +
         (32*m32)/(mQ32*pow2(mQ32 - mU32)) + (32*m32)/(mU32*pow2(mQ32 - mU32))
         )*pow4(Xt) + dilog(1 - m32/mQ32)*((-64*pow2(m32))/pow2(-m32 + mQ32) +
         (256*Xt*pow3(m3))/((m32 - mQ32)*(mQ32 - mU32)) + (128*m3*(2*m32 - mQ32
         - mU32)*pow3(Xt))/pow3(mQ32 - mU32) - (32*(-2*m32 + mQ32 + mU32)*pow4
         (Xt))/pow3(mQ32 - mU32)) + dilog(1 - m32/mU32)*((-64*pow2(m32))/pow2(
         m32 - mU32) + (256*Xt*pow3(m3))/((m32 - mU32)*(-mQ32 + mU32)) + (128*
         m3*(-2*m32 + mQ32 + mU32)*pow3(Xt))/pow3(mQ32 - mU32) + (32*(-2*m32 +
         mQ32 + mU32)*pow4(Xt))/pow3(mQ32 - mU32)) + log(mU32/MR2)*((32*(-2*m32
         *mU32 + 3*pow2(m32) + pow2(mU32)))/pow2(m32 - mU32) + (64*(3*m32 - 2*
         mU32)*pow2(Xt))/((m32 - mU32)*(-mQ32 + mU32)) - (64*Xt*(-(m3*mU32) + 2
         *pow3(m3)))/((m32 - mU32)*(-mQ32 + mU32)) + (128*(m3*mQ32 + 3*m3*mU32)
         *pow3(Xt))/pow3(-mQ32 + mU32) - (32*(2*m32 + 3*mQ32 + 7*mU32)*pow4(Xt)
         )/pow3(-mQ32 + mU32)) + pow2(log(mQ32/MR2))*((-16*(-2*m32*mQ32 + 3*
         pow2(m32) + pow2(mQ32)))/pow2(m32 - mQ32) - (32*(3*mQ32 - mU32)*pow2(
         Xt))/pow2(mQ32 - mU32) + (128*Xt*pow3(m3))/((m32 - mQ32)*(mQ32 - mU32)
         ) + (64*m3*(2*m32 - mQ32 - mU32)*pow3(Xt))/pow3(mQ32 - mU32) + (32*(
         m32*mQ32 - m32*mU32 + 2*mQ32*mU32 + 2*pow2(mQ32))*pow4(Xt))/pow4(mQ32
         - mU32)) + pow2(log(mU32/MR2))*((-16*(-2*m32*mU32 + 3*pow2(m32) + pow2
         (mU32)))/pow2(m32 - mU32) - (32*(-mQ32 + 3*mU32)*pow2(Xt))/pow2(-mQ32
         + mU32) + (128*Xt*pow3(m3))/((m32 - mU32)*(-mQ32 + mU32)) + (64*m3*(2*
         m32 - mQ32 - mU32)*pow3(Xt))/pow3(-mQ32 + mU32) + (32*(-(m32*mQ32) +
         m32*mU32 + 2*mQ32*mU32 + 2*pow2(mU32))*pow4(Xt))/pow4(-mQ32 + mU32)) +
         log(mQ32/MR2)*((32*(-2*m32*mQ32 + 3*pow2(m32) + pow2(mQ32)))/pow2(m32
         - mQ32) + (64*(3*m32 - 2*mQ32)*pow2(Xt))/((m32 - mQ32)*(mQ32 - mU32))
         - (64*Xt*(-(m3*mQ32) + 2*pow3(m3)))/((m32 - mQ32)*(mQ32 - mU32)) + (
         128*(3*m3*mQ32 + m3*mU32)*pow3(Xt))/pow3(mQ32 - mU32) - (32*(2*m32 + 7
         *mQ32 + 3*mU32)*pow4(Xt))/pow3(mQ32 - mU32) + log(mU32/MR2)*((64*(mQ32
         + mU32)*pow2(Xt))/pow2(mQ32 - mU32) - (64*(2*mQ32*mU32 + pow2(mQ32) +
         pow2(mU32))*pow4(Xt))/pow4(mQ32 - mU32))) + log(m32/MR2)*((64*(-m32 +
         mQ32 + mU32)*pow2(m32)*pow2(Xt))/(mQ32*(-m32 + mQ32)*(m32 - mU32)*
         mU32) + (64*Xt*pow3(m3))/((m32 - mQ32)*(m32 - mU32)) + log(mQ32/MR2)*(
         (32*pow2(m32))/pow2(m32 - mQ32) - (128*Xt*pow3(m3))/((m32 - mQ32)*(
         mQ32 - mU32)) - (256*pow3(m3)*pow3(Xt))/pow3(mQ32 - mU32)) + log(
         mU32/MR2)*((32*pow2(m32))/pow2(m32 - mU32) - (128*Xt*pow3(m3))/((m32 -
         mU32)*(-mQ32 + mU32)) - (256*pow3(m3)*pow3(Xt))/pow3(-mQ32 + mU32)) -
         (32*(m32*mQ32 + m32*mU32)*pow4(Xt))/(mQ32*mU32*pow2(mQ32 - mU32)) - (
         32*(-2*mU32*pow2(mQ32)*pow3(m32) - 2*mQ32*pow2(mU32)*pow3(m32) + mU32*
         pow2(m32)*pow3(mQ32) + pow3(m32)*pow3(mQ32) + mQ32*pow2(m32)*pow3(mU32
         ) + pow3(m32)*pow3(mU32) + 2*mQ32*mU32*pow4(m32) - 2*pow2(mQ32)*pow4(
         m32) - 2*pow2(mU32)*pow4(m32) + mQ32*pow5(m32) + mU32*pow5(m32)))/(
         mQ32*mU32*pow2(m32 - mQ32)*pow2(m32 - mU32)));

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
   const double dlatas2 = 0.; // 3L threshold correction \[CapitalDelta]\[Lambda]^(3L) O(at*as^2)

   const double result =
      -772.8212616582708 + dlatas2 + (128*mQ32*mU32)/(3.*(m32 - mQ32)*
         (m32 - mU32)) + 320*zt3 - (128*pow2(m32))/(3.*(m32 - mQ32)*(m32 - mU32
         )) - (16*(-(mQ32*mU32) + pow2(m32)))/((m32 - mQ32)*(m32 - mU32)) - (64
         *m32*(3*m32*mQ32 + 3*m32*mU32 + 5*mQ32*mU32 - 11*pow2(m32))*pow2(Xt))/
         (3.*(m32 - mQ32)*mQ32*(m32 - mU32)*mU32) + (2816*pow2(zt2))/15. + (256
         *(6*zt2 - pow2(log(2)))*pow2(log(2)))/9. - (48*pow2(-(mQ32*mU32) +
         pow2(m32)))/(pow2(m32 - mQ32)*pow2(m32 - mU32)) - (64*Xt*(-3*m3*mQ32 -
         30*m3*msq2 - 3*m3*mU32 + 14*pow3(m3)))/(3.*(m32 - mQ32)*(m32 - mU32))
         - (16*(3*m32*mQ32 + 3*m32*mU32 + 5*mQ32*mU32 - 11*pow2(m32))*(9*mQ32*
         mU32*pow2(m32) - 6*m32*mU32*pow2(mQ32) + 2*pow2(m32)*pow2(mQ32) - 6*
         m32*mQ32*pow2(mU32) + 2*pow2(m32)*pow2(mU32) + 3*pow2(mQ32)*pow2(mU32)
         - 2*mQ32*pow3(m32) - 2*mU32*pow3(m32)))/(3.*mQ32*mU32*pow2(-m32 +
         mQ32)*pow2(m32 - mU32)) - (512*m3*(3*m32*mQ32 + 3*m32*mU32 + 5*mQ32*
         mU32 - 11*pow2(m32))*pow3(Xt))/(3.*(m32 - mQ32)*(m32 - mU32)*pow2(mQ32
         - mU32)) + (4*(240*m32*mQ32*msq2*mU32 + 300*mQ32*msq2*pow2(m32) - 112
         *mQ32*mU32*pow2(m32) + 300*msq2*mU32*pow2(m32) - 180*m32*msq2*pow2(
         mQ32) + 208*m32*mU32*pow2(mQ32) - 60*msq2*mU32*pow2(mQ32) + 93*pow2(
         m32)*pow2(mQ32) + 208*m32*mQ32*pow2(mU32) - 180*m32*msq2*pow2(mU32) -
         60*mQ32*msq2*pow2(mU32) + 93*pow2(m32)*pow2(mU32) - 169*pow2(mQ32)*
         pow2(mU32) - 282*mQ32*pow3(m32) - 360*msq2*pow3(m32) - 282*mU32*pow3(
         m32) - 18*m32*pow3(mQ32) - 6*mU32*pow3(mQ32) - 18*m32*pow3(mU32) - 6*
         mQ32*pow3(mU32) + 291*pow4(m32)))/(3.*pow2(m32 - mQ32)*pow2(m32 - mU32
         )) + (8*pow4(Xt)*(-300*msq2*mU32*pow2(m32)*pow2(mQ32) - 300*mQ32*msq2*
         pow2(m32)*pow2(mU32) - 240*m32*msq2*pow2(mQ32)*pow2(mU32) - 168*pow2(
         m32)*pow2(mQ32)*pow2(mU32) + 360*mQ32*msq2*mU32*pow3(m32) + 558*mU32*
         pow2(mQ32)*pow3(m32) + 558*mQ32*pow2(mU32)*pow3(m32) + 180*m32*msq2*
         mU32*pow3(mQ32) - 173*mU32*pow2(m32)*pow3(mQ32) - 236*m32*pow2(mU32)*
         pow3(mQ32) + 60*msq2*pow2(mU32)*pow3(mQ32) - 12*pow3(m32)*pow3(mQ32) +
         180*m32*mQ32*msq2*pow3(mU32) - 173*mQ32*pow2(m32)*pow3(mU32) - 236*
         m32*pow2(mQ32)*pow3(mU32) + 60*msq2*pow2(mQ32)*pow3(mU32) - 12*pow3(
         m32)*pow3(mU32) + 277*pow3(mQ32)*pow3(mU32) - 455*mQ32*mU32*pow4(m32)
         + 56*pow2(mQ32)*pow4(m32) + 56*pow2(mU32)*pow4(m32) + 18*m32*mU32*pow4
         (mQ32) + 6*pow2(mU32)*pow4(mQ32) + 18*m32*mQ32*pow4(mU32) + 6*pow2(
         mQ32)*pow4(mU32) - 44*mQ32*pow5(m32) - 44*mU32*pow5(m32)))/(3.*mQ32*
         mU32*pow2(-m32 + mQ32)*pow2(m32 - mU32)*pow2(mQ32 - mU32)) + (128*(-3*
         m3*mQ32 - 30*m3*msq2 - 3*m3*mU32 + 14*pow3(m3))*pow5(Xt))/(3.*(m32 -
         mQ32)*(m32 - mU32)*pow2(mQ32 - mU32)) + log(msq2/MR2)*((640*m3*msq2*Xt
         )/((m32 - mQ32)*(m32 - mU32)) - (96*m32*pow2(Xt))/(mQ32*mU32) - (768*
         m3*pow3(Xt))/pow2(mQ32 - mU32) + dilog(1 - m32/mQ32)*((-96*pow2(m32))
         /pow2(-m32 + mQ32) + (384*Xt*pow3(m3))/((m32 - mQ32)*(mQ32 - mU32)) +
         (192*m3*(2*m32 - mQ32 - mU32)*pow3(Xt))/pow3(mQ32 - mU32) - (48*(-2*
         m32 + mQ32 + mU32)*pow4(Xt))/pow3(mQ32 - mU32)) + dilog(1 - m32/mU32)*
         ((-96*pow2(m32))/pow2(m32 - mU32) + (384*Xt*pow3(m3))/((m32 - mU32)*(
         -mQ32 + mU32)) + (192*m3*(-2*m32 + mQ32 + mU32)*pow3(Xt))/pow3(mQ32 -
         mU32) + (48*(-2*m32 + mQ32 + mU32)*pow4(Xt))/pow3(mQ32 - mU32)) + (16*
         (75*msq2*mU32*pow2(m32)*pow2(mQ32) + 75*mQ32*msq2*pow2(m32)*pow2(mU32)
         + 60*m32*msq2*pow2(mQ32)*pow2(mU32) - 328*pow2(m32)*pow2(mQ32)*pow2(
         mU32) - 90*mQ32*msq2*mU32*pow3(m32) + 196*mU32*pow2(mQ32)*pow3(m32) +
         196*mQ32*pow2(mU32)*pow3(m32) - 45*m32*msq2*mU32*pow3(mQ32) - 91*mU32*
         pow2(m32)*pow3(mQ32) + 150*m32*pow2(mU32)*pow3(mQ32) - 15*msq2*pow2(
         mU32)*pow3(mQ32) + 9*pow3(m32)*pow3(mQ32) - 45*m32*mQ32*msq2*pow3(mU32
         ) - 91*mQ32*pow2(m32)*pow3(mU32) + 150*m32*pow2(mQ32)*pow3(mU32) - 15*
         msq2*pow2(mQ32)*pow3(mU32) + 9*pow3(m32)*pow3(mU32) - 68*pow3(mQ32)*
         pow3(mU32) - 114*mQ32*mU32*pow4(m32) - 18*pow2(mQ32)*pow4(m32) - 18*
         pow2(mU32)*pow4(m32) + 9*mQ32*pow5(m32) + 9*mU32*pow5(m32)))/(3.*mQ32*
         mU32*pow2(m32 - mQ32)*pow2(m32 - mU32)) + (16*pow4(Xt)*(-150*msq2*mU32
         *pow2(m32)*pow2(mQ32) - 150*mQ32*msq2*pow2(m32)*pow2(mU32) - 120*m32*
         msq2*pow2(mQ32)*pow2(mU32) + 620*pow2(m32)*pow2(mQ32)*pow2(mU32) + 180
         *mQ32*msq2*mU32*pow3(m32) - 293*mU32*pow2(mQ32)*pow3(m32) - 293*mQ32*
         pow2(mU32)*pow3(m32) + 90*m32*msq2*mU32*pow3(mQ32) + 146*mU32*pow2(m32
         )*pow3(mQ32) - 309*m32*pow2(mU32)*pow3(mQ32) + 30*msq2*pow2(mU32)*pow3
         (mQ32) + 9*pow3(m32)*pow3(mQ32) + 90*m32*mQ32*msq2*pow3(mU32) + 146*
         mQ32*pow2(m32)*pow3(mU32) - 309*m32*pow2(mQ32)*pow3(mU32) + 30*msq2*
         pow2(mQ32)*pow3(mU32) + 9*pow3(m32)*pow3(mU32) + 154*pow3(mQ32)*pow3(
         mU32) + 138*mQ32*mU32*pow4(m32) - 18*pow2(mQ32)*pow4(m32) - 18*pow2(
         mU32)*pow4(m32) + 9*mQ32*pow5(m32) + 9*mU32*pow5(m32)))/(3.*mQ32*mU32*
         pow2(-m32 + mQ32)*pow2(m32 - mU32)*pow2(mQ32 - mU32)) - (1280*m3*msq2*
         pow5(Xt))/((m32 - mQ32)*(m32 - mU32)*pow2(mQ32 - mU32))) + pow2(log(
         msq2/MR2))*(24*Xt*((40*m3*pow2(mQ32 - msq2))/(3.*(mQ32 - mU32)*pow2(
         m32 - mQ32)) + (40*m3*pow2(-msq2 + mU32))/(3.*(-mQ32 + mU32)*pow2(m32
         - mU32))) + 24*((-5*(2*mQ32 - 2*msq2)*(-3*m32*mQ32 + 3*m32*msq2 + mQ32
         *msq2 - 2*pow2(m32) + pow2(mQ32)))/(6.*pow3(-m32 + mQ32)) - (5*(-2*
         msq2 + 2*mU32)*(3*m32*msq2 - 3*m32*mU32 + msq2*mU32 - 2*pow2(m32) +
         pow2(mU32)))/(6.*pow3(-m32 + mU32))) - (48*((-5*(2*mQ32 - 2*msq2)*(-3*
         m32*mQ32 + 3*m32*msq2 + mQ32*msq2 - 2*pow2(m32) + pow2(mQ32)))/(6.*
         pow3(-m32 + mQ32)) - (5*(-2*msq2 + 2*mU32)*(3*m32*msq2 - 3*m32*mU32 +
         msq2*mU32 - 2*pow2(m32) + pow2(mU32)))/(6.*pow3(-m32 + mU32)))*pow4(Xt
         ))/pow2(mQ32 - mU32) - (48*((40*m3*pow2(mQ32 - msq2))/(3.*(mQ32 - mU32
         )*pow2(m32 - mQ32)) + (40*m3*pow2(-msq2 + mU32))/(3.*(-mQ32 + mU32)*
         pow2(m32 - mU32)))*pow5(Xt))/pow2(mQ32 - mU32)) + dilog(1 - m32/msq2)*
         ((640*(m32 - msq2)*Xt*(m3*mQ32*msq2 - 2*m3*mQ32*mU32 + m3*msq2*mU32 +
         mQ32*pow3(m3) - 2*msq2*pow3(m3) + mU32*pow3(m3)))/(pow2(m32 - mQ32)*
         pow2(m32 - mU32)) + (80*(m32 - msq2)*(-6*mQ32*msq2*mU32*pow2(m32) + 3*
         m32*msq2*mU32*pow2(mQ32) + 9*msq2*pow2(m32)*pow2(mQ32) - 15*mU32*pow2(
         m32)*pow2(mQ32) + 3*m32*mQ32*msq2*pow2(mU32) - 15*mQ32*pow2(m32)*pow2(
         mU32) + 9*msq2*pow2(m32)*pow2(mU32) - 8*mQ32*msq2*pow3(m32) + 30*mQ32*
         mU32*pow3(m32) - 8*msq2*mU32*pow3(m32) + 3*pow2(mQ32)*pow3(m32) + 3*
         pow2(mU32)*pow3(m32) - 3*m32*msq2*pow3(mQ32) + 5*m32*mU32*pow3(mQ32) -
         msq2*mU32*pow3(mQ32) - pow2(m32)*pow3(mQ32) + 5*m32*mQ32*pow3(mU32) -
         3*m32*msq2*pow3(mU32) - mQ32*msq2*pow3(mU32) - pow2(m32)*pow3(mU32) -
         8*mQ32*pow4(m32) + 6*msq2*pow4(m32) - 8*mU32*pow4(m32) + 2*pow5(m32))
         )/(pow3(m32 - mQ32)*pow3(m32 - mU32)) - (160*(m32 - msq2)*pow4(Xt)*(-6
         *mQ32*msq2*mU32*pow2(m32) + 3*m32*msq2*mU32*pow2(mQ32) + 9*msq2*pow2(
         m32)*pow2(mQ32) - 15*mU32*pow2(m32)*pow2(mQ32) + 3*m32*mQ32*msq2*pow2(
         mU32) - 15*mQ32*pow2(m32)*pow2(mU32) + 9*msq2*pow2(m32)*pow2(mU32) - 8
         *mQ32*msq2*pow3(m32) + 30*mQ32*mU32*pow3(m32) - 8*msq2*mU32*pow3(m32)
         + 3*pow2(mQ32)*pow3(m32) + 3*pow2(mU32)*pow3(m32) - 3*m32*msq2*pow3(
         mQ32) + 5*m32*mU32*pow3(mQ32) - msq2*mU32*pow3(mQ32) - pow2(m32)*pow3(
         mQ32) + 5*m32*mQ32*pow3(mU32) - 3*m32*msq2*pow3(mU32) - mQ32*msq2*pow3
         (mU32) - pow2(m32)*pow3(mU32) - 8*mQ32*pow4(m32) + 6*msq2*pow4(m32) -
         8*mU32*pow4(m32) + 2*pow5(m32)))/(pow2(mQ32 - mU32)*pow3(m32 - mQ32)*
         pow3(m32 - mU32)) - (1280*(m32 - msq2)*(m3*mQ32*msq2 - 2*m3*mQ32*mU32
         + m3*msq2*mU32 + mQ32*pow3(m3) - 2*msq2*pow3(m3) + mU32*pow3(m3))*pow5
         (Xt))/(pow2(m32 - mQ32)*pow2(m32 - mU32)*pow2(mQ32 - mU32))) + dilog(1
         - msq2/mQ32)*((640*m3*Xt*pow2(mQ32 - msq2))/((mQ32 - mU32)*pow2(m32 -
         mQ32)) - (80*(mQ32 - msq2)*(3*m32*mQ32 - 3*m32*msq2 - mQ32*msq2 + 2*
         pow2(m32) - pow2(mQ32)))/pow3(m32 - mQ32) + (160*(mQ32 - msq2)*(3*m32*
         mQ32 - 3*m32*msq2 - mQ32*msq2 + 2*pow2(m32) - pow2(mQ32))*pow4(Xt))/(
         pow2(mQ32 - mU32)*pow3(m32 - mQ32)) - (1280*m3*pow2(mQ32 - msq2)*pow5(
         Xt))/(pow2(m32 - mQ32)*pow3(mQ32 - mU32))) + dilog(1 - msq2/mU32)*((
         640*m3*Xt*pow2(-msq2 + mU32))/((-mQ32 + mU32)*pow2(m32 - mU32)) - (80*
         (-msq2 + mU32)*(-3*m32*msq2 + 3*m32*mU32 - msq2*mU32 + 2*pow2(m32) -
         pow2(mU32)))/pow3(m32 - mU32) + (160*(-msq2 + mU32)*(-3*m32*msq2 + 3*
         m32*mU32 - msq2*mU32 + 2*pow2(m32) - pow2(mU32))*pow4(Xt))/(pow2(-mQ32
         + mU32)*pow3(m32 - mU32)) - (1280*m3*pow2(-msq2 + mU32)*pow5(Xt))/(
         pow2(m32 - mU32)*pow3(-mQ32 + mU32))) + dilog(1 - mQ32/mU32)*((-32*Xt*
         (-8*m3*mU32*pow2(mQ32) + 8*m3*mQ32*pow2(mU32) - 10*pow2(mQ32)*pow3(m3)
         + 10*pow2(mU32)*pow3(m3) + 6*m3*pow3(mQ32) - 6*m3*pow3(mU32) + 10*
         mQ32*pow5(m3) - 10*mU32*pow5(m3)))/(3.*pow2(m32 - mQ32)*pow2(m32 -
         mU32)) - (16*pow4(Xt)*(65*mU32*pow2(m32)*pow2(mQ32) + 65*mQ32*pow2(m32
         )*pow2(mU32) - 34*m32*pow2(mQ32)*pow2(mU32) - 76*mQ32*mU32*pow3(m32) +
         34*pow2(mQ32)*pow3(m32) + 34*pow2(mU32)*pow3(m32) - 26*m32*mU32*pow3(
         mQ32) - 29*pow2(m32)*pow3(mQ32) + 7*pow2(mU32)*pow3(mQ32) - 26*m32*
         mQ32*pow3(mU32) - 29*pow2(m32)*pow3(mU32) + 7*pow2(mQ32)*pow3(mU32) -
         14*mQ32*pow4(m32) - 14*mU32*pow4(m32) + 9*m32*pow4(mQ32) + 3*mU32*pow4
         (mQ32) + 9*m32*pow4(mU32) + 3*mQ32*pow4(mU32) + 12*pow5(m32)))/(3.*(
         mQ32 - mU32)*pow3(m32 - mQ32)*pow3(m32 - mU32)) + (4*(220*mU32*pow2(
         mQ32)*pow3(m32) - 220*mQ32*pow2(mU32)*pow3(m32) - 188*mU32*pow2(m32)*
         pow3(mQ32) + 16*m32*pow2(mU32)*pow3(mQ32) - 68*pow3(m32)*pow3(mQ32) +
         188*mQ32*pow2(m32)*pow3(mU32) - 16*m32*pow2(mQ32)*pow3(mU32) + 68*pow3
         (m32)*pow3(mU32) + 28*pow2(mQ32)*pow4(m32) - 28*pow2(mU32)*pow4(m32) +
         70*m32*mU32*pow4(mQ32) + 58*pow2(m32)*pow4(mQ32) - 8*pow2(mU32)*pow4(
         mQ32) - 70*m32*mQ32*pow4(mU32) - 58*pow2(m32)*pow4(mU32) + 8*pow2(mQ32
         )*pow4(mU32) - 24*mQ32*pow5(m32) + 24*mU32*pow5(m32) - 18*m32*pow5(
         mQ32) - 6*mU32*pow5(mQ32) + 18*m32*pow5(mU32) + 6*mQ32*pow5(mU32)))/(
         3.*pow3(-m32 + mQ32)*pow3(m32 - mU32)) + (128*(-(m3*mQ32*mU32) + 3*m3*
         pow2(mQ32) + 3*m3*pow2(mU32) - 5*mQ32*pow3(m3) - 5*mU32*pow3(m3) + 5*
         pow5(m3))*pow5(Xt))/(3.*(mQ32 - mU32)*pow2(m32 - mQ32)*pow2(m32 - mU32
         ))) + dilog(1 - m32/mQ32)*((128*m3*(2*m32 - mQ32 - mU32)*(3*m32*mQ32 +
         3*m32*mU32 + 5*mQ32*mU32 - 11*pow2(m32))*pow3(Xt))/(3.*(m32 - mQ32)*(
         m32 - mU32)*pow3(mQ32 - mU32)) + Xt*((256*(3*m32*mQ32 + 3*m32*mU32 + 5
         *mQ32*mU32 - 11*pow2(m32))*pow3(m3))/(3.*(m32 - mU32)*(mQ32 - mU32)*
         pow2(m32 - mQ32)) + (64*(-7*m3*mQ32*mU32 + 3*m3*pow2(mQ32) + 22*m3*
         pow2(mU32) + mQ32*pow3(m3) - 37*mU32*pow3(m3) + 18*pow5(m3)))/(3.*(
         -mQ32 + mU32)*pow2(m32 - mU32))) - (128*(-7*m3*mQ32*mU32 + 3*m3*pow2(
         mQ32) + 22*m3*pow2(mU32) + mQ32*pow3(m3) - 37*mU32*pow3(m3) + 18*pow5(
         m3))*pow5(Xt))/(3.*pow2(m32 - mU32)*pow3(-mQ32 + mU32)) - (16*pow4(Xt)
         *(411*pow2(m32)*pow2(mQ32)*pow2(mU32) - 313*mU32*pow2(mQ32)*pow3(m32)
         - 1151*mQ32*pow2(mU32)*pow3(m32) - 149*mU32*pow2(m32)*pow3(mQ32) + 55*
         m32*pow2(mU32)*pow3(mQ32) + 3*pow3(m32)*pow3(mQ32) + 409*mQ32*pow2(m32
         )*pow3(mU32) - 129*m32*pow2(mQ32)*pow3(mU32) - 59*pow3(m32)*pow3(mU32)
         - 29*pow3(mQ32)*pow3(mU32) + 1084*mQ32*mU32*pow4(m32) + 138*pow2(mQ32
         )*pow4(m32) + 398*pow2(mU32)*pow4(m32) + 41*m32*mU32*pow4(mQ32) + 20*
         pow2(m32)*pow4(mQ32) - pow2(mU32)*pow4(mQ32) - 30*m32*mQ32*pow4(mU32)
         - 31*pow2(m32)*pow4(mU32) + 13*pow2(mQ32)*pow4(mU32) - 372*mQ32*pow5(
         m32) - 468*mU32*pow5(m32) - 9*m32*pow5(mQ32) - 3*mU32*pow5(mQ32) + 172
         *pow6(m32)))/(3.*pow2(m32 - mQ32)*pow3(m32 - mU32)*pow3(mQ32 - mU32))
         + (8*(1498*pow2(mQ32)*pow2(mU32)*pow3(m32) - 376*pow2(m32)*pow2(mU32)*
         pow3(mQ32) + 100*mU32*pow3(m32)*pow3(mQ32) - 598*pow2(m32)*pow2(mQ32)*
         pow3(mU32) + 324*mQ32*pow3(m32)*pow3(mU32) + 148*m32*pow3(mQ32)*pow3(
         mU32) - 1145*mU32*pow2(mQ32)*pow4(m32) - 1177*mQ32*pow2(mU32)*pow4(m32
         ) - 89*pow3(mQ32)*pow4(m32) - 149*pow3(mU32)*pow4(m32) + 192*mU32*pow2
         (m32)*pow4(mQ32) - 42*m32*pow2(mU32)*pow4(mQ32) + 11*pow3(m32)*pow4(
         mQ32) + 19*pow3(mU32)*pow4(mQ32) + 43*mQ32*pow2(m32)*pow4(mU32) + 57*
         m32*pow2(mQ32)*pow4(mU32) - 13*pow3(m32)*pow4(mU32) - 23*pow3(mQ32)*
         pow4(mU32) + 1072*mQ32*mU32*pow5(m32) + 368*pow2(mQ32)*pow5(m32) + 480
         *pow2(mU32)*pow5(m32) - 44*m32*mU32*pow5(mQ32) - 29*pow2(m32)*pow5(
         mQ32) + pow2(mU32)*pow5(mQ32) - 334*mQ32*pow6(m32) - 434*mU32*pow6(m32
         ) + 9*m32*pow6(mQ32) + 3*mU32*pow6(mQ32) + 128*pow7(m32)))/(3.*(mQ32 -
         mU32)*pow3(m32 - mQ32)*pow3(m32 - mU32))) + dilog(1 - m32/mU32)*((128
         *m3*(-2*m32 + mQ32 + mU32)*(3*m32*mQ32 + 3*m32*mU32 + 5*mQ32*mU32 - 11
         *pow2(m32))*pow3(Xt))/(3.*(m32 - mQ32)*(m32 - mU32)*pow3(mQ32 - mU32))
         + Xt*((256*(3*m32*mQ32 + 3*m32*mU32 + 5*mQ32*mU32 - 11*pow2(m32))*
         pow3(m3))/(3.*(m32 - mQ32)*(-mQ32 + mU32)*pow2(m32 - mU32)) + (64*(-7*
         m3*mQ32*mU32 + 22*m3*pow2(mQ32) + 3*m3*pow2(mU32) - 37*mQ32*pow3(m3) +
         mU32*pow3(m3) + 18*pow5(m3)))/(3.*(mQ32 - mU32)*pow2(m32 - mQ32))) -
         (128*(-7*m3*mQ32*mU32 + 22*m3*pow2(mQ32) + 3*m3*pow2(mU32) - 37*mQ32*
         pow3(m3) + mU32*pow3(m3) + 18*pow5(m3))*pow5(Xt))/(3.*pow2(m32 - mQ32)
         *pow3(mQ32 - mU32)) + (16*pow4(Xt)*(411*pow2(m32)*pow2(mQ32)*pow2(mU32
         ) - 1151*mU32*pow2(mQ32)*pow3(m32) - 313*mQ32*pow2(mU32)*pow3(m32) +
         409*mU32*pow2(m32)*pow3(mQ32) - 129*m32*pow2(mU32)*pow3(mQ32) - 59*
         pow3(m32)*pow3(mQ32) - 149*mQ32*pow2(m32)*pow3(mU32) + 55*m32*pow2(
         mQ32)*pow3(mU32) + 3*pow3(m32)*pow3(mU32) - 29*pow3(mQ32)*pow3(mU32) +
         1084*mQ32*mU32*pow4(m32) + 398*pow2(mQ32)*pow4(m32) + 138*pow2(mU32)*
         pow4(m32) - 30*m32*mU32*pow4(mQ32) - 31*pow2(m32)*pow4(mQ32) + 13*pow2
         (mU32)*pow4(mQ32) + 41*m32*mQ32*pow4(mU32) + 20*pow2(m32)*pow4(mU32) -
         pow2(mQ32)*pow4(mU32) - 468*mQ32*pow5(m32) - 372*mU32*pow5(m32) - 9*
         m32*pow5(mU32) - 3*mQ32*pow5(mU32) + 172*pow6(m32)))/(3.*pow2(m32 -
         mU32)*pow3(m32 - mQ32)*pow3(mQ32 - mU32)) - (8*(1498*pow2(mQ32)*pow2(
         mU32)*pow3(m32) - 598*pow2(m32)*pow2(mU32)*pow3(mQ32) + 324*mU32*pow3(
         m32)*pow3(mQ32) - 376*pow2(m32)*pow2(mQ32)*pow3(mU32) + 100*mQ32*pow3(
         m32)*pow3(mU32) + 148*m32*pow3(mQ32)*pow3(mU32) - 1177*mU32*pow2(mQ32)
         *pow4(m32) - 1145*mQ32*pow2(mU32)*pow4(m32) - 149*pow3(mQ32)*pow4(m32)
         - 89*pow3(mU32)*pow4(m32) + 43*mU32*pow2(m32)*pow4(mQ32) + 57*m32*
         pow2(mU32)*pow4(mQ32) - 13*pow3(m32)*pow4(mQ32) - 23*pow3(mU32)*pow4(
         mQ32) + 192*mQ32*pow2(m32)*pow4(mU32) - 42*m32*pow2(mQ32)*pow4(mU32) +
         11*pow3(m32)*pow4(mU32) + 19*pow3(mQ32)*pow4(mU32) + 1072*mQ32*mU32*
         pow5(m32) + 480*pow2(mQ32)*pow5(m32) + 368*pow2(mU32)*pow5(m32) - 44*
         m32*mQ32*pow5(mU32) - 29*pow2(m32)*pow5(mU32) + pow2(mQ32)*pow5(mU32)
         - 434*mQ32*pow6(m32) - 334*mU32*pow6(m32) + 9*m32*pow6(mU32) + 3*mQ32*
         pow6(mU32) + 128*pow7(m32)))/(3.*(mQ32 - mU32)*pow3(m32 - mQ32)*pow3(
         m32 - mU32))) + pow3(log(mU32/MR2))*((-4*(-120*mQ32*msq2*mU32*pow2(m32
         ) - 33*mU32*pow2(m32)*pow2(mQ32) - 60*m32*mQ32*mU32*pow2(msq2) + 90*
         mQ32*pow2(m32)*pow2(msq2) - 90*mU32*pow2(m32)*pow2(msq2) + 180*m32*
         mQ32*msq2*pow2(mU32) - 429*mQ32*pow2(m32)*pow2(mU32) + 120*msq2*pow2(
         m32)*pow2(mU32) + 36*m32*pow2(mQ32)*pow2(mU32) + 60*m32*pow2(msq2)*
         pow2(mU32) - 30*mQ32*pow2(msq2)*pow2(mU32) - 60*mQ32*msq2*pow3(m32) +
         276*mQ32*mU32*pow3(m32) + 60*msq2*mU32*pow3(m32) - 2*pow2(mQ32)*pow3(
         m32) + 494*pow2(mU32)*pow3(m32) - 6*m32*mU32*pow3(mQ32) + 9*pow2(m32)*
         pow3(mQ32) - 3*pow2(mU32)*pow3(mQ32) + 234*m32*mQ32*pow3(mU32) - 180*
         m32*msq2*pow3(mU32) - 59*pow2(m32)*pow3(mU32) - pow2(mQ32)*pow3(mU32)
         + 30*pow2(msq2)*pow3(mU32) + 38*mQ32*pow4(m32) - 550*mU32*pow4(m32) -
         136*m32*pow4(mU32) - 63*mQ32*pow4(mU32) + 128*pow5(m32) + 67*pow5(mU32
         )))/(3.*(-mQ32 + mU32)*pow4(m32 - mU32)) - (8*pow2(Xt)*(-120*mQ32*msq2
         *mU32*pow2(m32) - 33*mU32*pow2(m32)*pow2(mQ32) - 60*m32*mQ32*mU32*pow2
         (msq2) + 90*mQ32*pow2(m32)*pow2(msq2) - 90*mU32*pow2(m32)*pow2(msq2) +
         180*m32*mQ32*msq2*pow2(mU32) - 361*mQ32*pow2(m32)*pow2(mU32) + 120*
         msq2*pow2(m32)*pow2(mU32) + 36*m32*pow2(mQ32)*pow2(mU32) + 60*m32*pow2
         (msq2)*pow2(mU32) - 30*mQ32*pow2(msq2)*pow2(mU32) - 60*mQ32*msq2*pow3(
         m32) + 140*mQ32*mU32*pow3(m32) + 60*msq2*mU32*pow3(m32) - 2*pow2(mQ32)
         *pow3(m32) + 1670*pow2(mU32)*pow3(m32) - 6*m32*mU32*pow3(mQ32) + 9*
         pow2(m32)*pow3(mQ32) - 3*pow2(mU32)*pow3(mQ32) + 234*m32*mQ32*pow3(
         mU32) - 180*m32*msq2*pow3(mU32) - 615*pow2(m32)*pow3(mU32) - pow2(mQ32
         )*pow3(mU32) + 30*pow2(msq2)*pow3(mU32) + 42*mQ32*pow4(m32) - 1062*
         mU32*pow4(m32) - 248*m32*pow4(mU32) - 63*mQ32*pow4(mU32) + 128*pow5(
         m32) + 135*pow5(mU32)))/(3.*pow2(-mQ32 + mU32)*pow4(m32 - mU32)) + (4*
         pow4(Xt)*(-120*msq2*mU32*pow2(m32)*pow2(mQ32) - 60*m32*mU32*pow2(mQ32)
         *pow2(msq2) + 90*pow2(m32)*pow2(mQ32)*pow2(msq2) + 180*m32*msq2*pow2(
         mQ32)*pow2(mU32) - 222*pow2(m32)*pow2(mQ32)*pow2(mU32) - 90*pow2(m32)*
         pow2(msq2)*pow2(mU32) - 30*pow2(mQ32)*pow2(msq2)*pow2(mU32) - 60*msq2*
         pow2(mQ32)*pow3(m32) + 66*mU32*pow2(mQ32)*pow3(m32) + 2314*mQ32*pow2(
         mU32)*pow3(m32) + 60*msq2*pow2(mU32)*pow3(m32) - 24*mU32*pow2(m32)*
         pow3(mQ32) + 30*m32*pow2(mU32)*pow3(mQ32) - 2*pow3(m32)*pow3(mQ32) -
         2344*mQ32*pow2(m32)*pow3(mU32) + 120*msq2*pow2(m32)*pow3(mU32) + 134*
         m32*pow2(mQ32)*pow3(mU32) + 60*m32*pow2(msq2)*pow3(mU32) + 4726*pow3(
         m32)*pow3(mU32) - 4*pow3(mQ32)*pow3(mU32) - 872*mQ32*mU32*pow4(m32) +
         44*pow2(mQ32)*pow4(m32) - 4276*pow2(mU32)*pow4(m32) - 6*m32*mU32*pow4(
         mQ32) + 9*pow2(m32)*pow4(mQ32) - 3*pow2(mU32)*pow4(mQ32) + 638*m32*
         mQ32*pow4(mU32) - 180*m32*msq2*pow4(mU32) - 1163*pow2(m32)*pow4(mU32)
         - 30*pow2(mQ32)*pow4(mU32) + 30*pow2(msq2)*pow4(mU32) + 124*mQ32*pow5(
         m32) + 1156*mU32*pow5(m32) - 604*m32*pow5(mU32) + 140*mQ32*pow5(mU32)
         + 169*pow6(mU32)))/(3.*pow4(m32 - mU32)*pow4(-mQ32 + mU32)) + (1280*
         m32*(mQ32 + mU32)*pow2(mU32)*pow6(Xt))/(3.*pow2(m32 - mU32)*pow5(-mQ32
         + mU32)) + (32*Xt*(-30*m3*mQ32*mU32*pow2(msq2) + 60*m3*mQ32*msq2*pow2
         (mU32) + 10*m3*pow2(mQ32)*pow2(mU32) + 30*m3*pow2(msq2)*pow2(mU32) -
         60*mQ32*msq2*mU32*pow3(m3) - 11*mU32*pow2(mQ32)*pow3(m3) + 30*mQ32*
         pow2(msq2)*pow3(m3) - 30*mU32*pow2(msq2)*pow3(m3) - 149*mQ32*pow2(mU32
         )*pow3(m3) + 60*msq2*pow2(mU32)*pow3(m3) - 3*m3*mU32*pow3(mQ32) + 3*
         pow3(m3)*pow3(mQ32) + 32*m3*mQ32*pow3(mU32) - 60*m3*msq2*pow3(mU32) +
         189*pow3(m3)*pow3(mU32) - 55*m3*pow4(mU32) + 185*mQ32*mU32*pow5(m3) +
         pow2(mQ32)*pow5(m3) - 202*pow2(mU32)*pow5(m3) - 20*mQ32*pow7(m3) + 20*
         mU32*pow7(m3)))/(3.*pow2(-mQ32 + mU32)*pow3(m32 - mU32)) + (32*pow5(Xt
         )*(30*m3*mU32*pow2(mQ32)*pow2(msq2) - 60*m3*msq2*pow2(mQ32)*pow2(mU32)
         + 60*msq2*mU32*pow2(mQ32)*pow3(m3) - 30*pow2(mQ32)*pow2(msq2)*pow3(m3
         ) + 94*pow2(mQ32)*pow2(mU32)*pow3(m3) + 30*pow2(msq2)*pow2(mU32)*pow3(
         m3) - 7*m3*pow2(mU32)*pow3(mQ32) + 8*mU32*pow3(m3)*pow3(mQ32) - 26*m3*
         pow2(mQ32)*pow3(mU32) - 30*m3*pow2(msq2)*pow3(mU32) - 200*mQ32*pow3(m3
         )*pow3(mU32) - 60*msq2*pow3(m3)*pow3(mU32) + 3*m3*mU32*pow4(mQ32) - 3*
         pow3(m3)*pow4(mQ32) + 87*m3*mQ32*pow4(mU32) + 60*m3*msq2*pow4(mU32) -
         219*pow3(m3)*pow4(mU32) - 70*mU32*pow2(mQ32)*pow5(m3) + 145*mQ32*pow2(
         mU32)*pow5(m3) - pow3(mQ32)*pow5(m3) + 86*pow3(mU32)*pow5(m3) + 103*m3
         *pow5(mU32) - 32*mQ32*mU32*pow7(m3) + 18*pow2(mQ32)*pow7(m3) + 14*pow2
         (mU32)*pow7(m3)))/(3.*pow3(m32 - mU32)*pow5(-mQ32 + mU32)) + (32*pow3(
         Xt)*(-60*m3*mQ32*mU32*pow2(msq2) + 120*m3*mQ32*msq2*pow2(mU32) + 20*m3
         *pow2(mQ32)*pow2(mU32) + 60*m3*pow2(msq2)*pow2(mU32) - 120*mQ32*msq2*
         mU32*pow3(m3) - 22*mU32*pow2(mQ32)*pow3(m3) + 60*mQ32*pow2(msq2)*pow3(
         m3) - 60*mU32*pow2(msq2)*pow3(m3) - 281*mQ32*pow2(mU32)*pow3(m3) + 120
         *msq2*pow2(mU32)*pow3(m3) - 6*m3*mU32*pow3(mQ32) + 6*pow3(m3)*pow3(
         mQ32) + 81*m3*mQ32*pow3(mU32) - 120*m3*msq2*pow3(mU32) + 353*pow3(m3)*
         pow3(mU32) - 157*m3*pow4(mU32) + 205*mQ32*mU32*pow5(m3) + 2*pow2(mQ32)
         *pow5(m3) - 131*pow2(mU32)*pow5(m3) - 37*mQ32*pow7(m3) - 35*mU32*pow7(
         m3) + 2*pow9(m3)))/(3.*pow3(m32 - mU32)*pow3(-mQ32 + mU32))) + pow2(
         log(mU32/MR2))*((-16*pow2(Xt)*(-180*msq2*mU32*pow2(m32)*pow2(mQ32) +
         210*mQ32*msq2*pow2(m32)*pow2(mU32) - 150*m32*msq2*pow2(mQ32)*pow2(mU32
         ) - 1612*pow2(m32)*pow2(mQ32)*pow2(mU32) + 90*mQ32*msq2*mU32*pow3(m32)
         + 1346*mU32*pow2(mQ32)*pow3(m32) + 1473*mQ32*pow2(mU32)*pow3(m32) -
         90*msq2*pow2(mU32)*pow3(m32) + 90*m32*msq2*mU32*pow3(mQ32) - 373*mU32*
         pow2(m32)*pow3(mQ32) + 471*m32*pow2(mU32)*pow3(mQ32) + 30*msq2*pow2(
         mU32)*pow3(mQ32) - 41*pow3(m32)*pow3(mQ32) + 60*m32*mQ32*msq2*pow3(
         mU32) + 137*mQ32*pow2(m32)*pow3(mU32) - 30*msq2*pow2(m32)*pow3(mU32) +
         338*m32*pow2(mQ32)*pow3(mU32) - 30*msq2*pow2(mQ32)*pow3(mU32) - 226*
         pow3(m32)*pow3(mU32) - 177*pow3(mQ32)*pow3(mU32) - 1462*mQ32*mU32*pow4
         (m32) + 65*pow2(mQ32)*pow4(m32) - 383*pow2(mU32)*pow4(m32) + 9*m32*
         mU32*pow4(mQ32) + 3*pow2(mU32)*pow4(mQ32) - 391*m32*mQ32*pow4(mU32) +
         192*pow2(m32)*pow4(mU32) + 151*pow2(mQ32)*pow4(mU32) - 24*mQ32*pow5(
         m32) + 492*mU32*pow5(m32) + 9*m32*pow5(mU32) + 3*mQ32*pow5(mU32)))/(3.
         *pow2(m32 - mQ32)*pow2(-mQ32 + mU32)*pow3(m32 - mU32)) + log(msq2/MR2)
         *((-8*(-180*m32*msq2*mU32 - 60*msq2*pow2(m32) - 65*mU32*pow2(m32) + 90
         *m32*pow2(msq2) + 30*mU32*pow2(msq2) + 57*m32*pow2(mU32) + 27*pow3(m32
         ) - 19*pow3(mU32)))/(3.*pow3(m32 - mU32)) + pow2(Xt)*((48*((mQ32 - 3*
         mU32)*pow2(m32) + 2*(mQ32 - 2*mU32)*pow2(mU32) + m32*(-4*mQ32*mU32 + 8
         *pow2(mU32))))/(pow2(m32 - mU32)*pow2(mQ32 - mU32)) - (16*(180*m32*
         msq2*mU32 + 60*msq2*pow2(m32) + 2*mU32*pow2(m32) - 90*m32*pow2(msq2) -
         30*mU32*pow2(msq2) - 3*m32*pow2(mU32) + pow3(mU32)))/(3.*(mQ32 - mU32
         )*pow3(m32 - mU32))) + ((-64*(120*m3*msq2*mU32 - 60*m3*pow2(msq2) - m3
         *pow2(mU32) + mU32*pow3(m3)))/(3.*pow2(m32 - mU32)*pow2(-mQ32 + mU32))
         + (96*m3*(3*mQ32*mU32 - m32*(mQ32 + 3*mU32) + 2*pow2(m32) - pow2(mU32
         )))/((m32 - mU32)*pow3(-mQ32 + mU32)))*pow3(Xt) + ((8*(mQ32 + mU32)*(
         180*m32*msq2*mU32 + 60*msq2*pow2(m32) + 2*mU32*pow2(m32) - 90*m32*pow2
         (msq2) - 30*mU32*pow2(msq2) - 3*m32*pow2(mU32) + pow3(mU32)))/(3.*pow3
         (m32 - mU32)*pow3(mQ32 - mU32)) + (24*(8*mQ32*mU32*pow2(m32) + 2*m32*
         mU32*(-5*mQ32*mU32 + pow2(mQ32) - 4*pow2(mU32)) - pow2(mQ32)*pow2(mU32
         ) - 2*(mQ32 - mU32)*pow3(m32) + 4*mQ32*pow3(mU32) + 5*pow4(mU32)))/(
         pow2(m32 - mU32)*pow4(mQ32 - mU32)))*pow4(Xt) + (64*Xt*(-60*m3*msq2*
         mU32 + 30*m3*pow2(msq2) + 5*m3*pow2(mU32) - 14*mU32*pow3(m3) + 9*pow5(
         m3)))/(3.*(-mQ32 + mU32)*pow2(m32 - mU32)) + (320*(mQ32 + mU32)*(12*m3
         *msq2*mU32 - 6*m3*pow2(msq2) - m3*pow2(mU32) + mU32*pow3(m3))*pow5(Xt)
         )/(3.*pow2(m32 - mU32)*pow4(-mQ32 + mU32))) - (2560*m32*pow2(mU32)*
         pow6(Xt))/(3.*pow2(m32 - mU32)*pow4(-mQ32 + mU32)) + (4*(30*mU32*pow2(
         m32)*pow2(mQ32)*pow2(msq2) - 300*msq2*pow2(m32)*pow2(mQ32)*pow2(mU32)
         - 150*mQ32*pow2(m32)*pow2(msq2)*pow2(mU32) + 120*m32*pow2(mQ32)*pow2(
         msq2)*pow2(mU32) + 660*msq2*mU32*pow2(mQ32)*pow3(m32) + 120*mQ32*mU32*
         pow2(msq2)*pow3(m32) - 180*pow2(mQ32)*pow2(msq2)*pow3(m32) - 300*mQ32*
         msq2*pow2(mU32)*pow3(m32) + 4812*pow2(mQ32)*pow2(mU32)*pow3(m32) + 60*
         pow2(msq2)*pow2(mU32)*pow3(m32) - 300*msq2*mU32*pow2(m32)*pow3(mQ32) -
         60*m32*mU32*pow2(msq2)*pow3(mQ32) + 90*pow2(m32)*pow2(msq2)*pow3(mQ32
         ) + 300*m32*msq2*pow2(mU32)*pow3(mQ32) - 1436*pow2(m32)*pow2(mU32)*
         pow3(mQ32) - 30*pow2(msq2)*pow2(mU32)*pow3(mQ32) - 60*msq2*pow3(m32)*
         pow3(mQ32) + 764*mU32*pow3(m32)*pow3(mQ32) + 660*mQ32*msq2*pow2(m32)*
         pow3(mU32) - 420*m32*msq2*pow2(mQ32)*pow3(mU32) - 2968*pow2(m32)*pow2(
         mQ32)*pow3(mU32) - 60*m32*mQ32*pow2(msq2)*pow3(mU32) + 30*pow2(m32)*
         pow2(msq2)*pow3(mU32) + 30*pow2(mQ32)*pow2(msq2)*pow3(mU32) + 2044*
         mQ32*pow3(m32)*pow3(mU32) - 300*msq2*pow3(m32)*pow3(mU32) + 992*m32*
         pow3(mQ32)*pow3(mU32) + 60*msq2*pow3(mQ32)*pow3(mU32) - 420*mQ32*msq2*
         mU32*pow4(m32) + 120*msq2*pow2(mQ32)*pow4(m32) - 2847*mU32*pow2(mQ32)*
         pow4(m32) + 90*mQ32*pow2(msq2)*pow4(m32) - 90*mU32*pow2(msq2)*pow4(m32
         ) - 4635*mQ32*pow2(mU32)*pow4(m32) + 300*msq2*pow2(mU32)*pow4(m32) +
         115*pow3(mQ32)*pow4(m32) - 313*pow3(mU32)*pow4(m32) - 39*mU32*pow2(m32
         )*pow4(mQ32) + 54*m32*pow2(mU32)*pow4(mQ32) - 20*pow3(m32)*pow4(mQ32)
         + 5*pow3(mU32)*pow4(mQ32) + 120*m32*mQ32*msq2*pow4(mU32) + 815*mQ32*
         pow2(m32)*pow4(mU32) - 60*msq2*pow2(m32)*pow4(mU32) + 284*m32*pow2(
         mQ32)*pow4(mU32) - 60*msq2*pow2(mQ32)*pow4(mU32) - 560*pow3(m32)*pow4(
         mU32) - 299*pow3(mQ32)*pow4(mU32) - 60*mQ32*msq2*pow5(m32) + 3132*mQ32
         *mU32*pow5(m32) + 60*msq2*mU32*pow5(m32) + 20*pow2(mQ32)*pow5(m32) +
         1328*pow2(mU32)*pow5(m32) - 6*m32*mU32*pow5(mQ32) + 9*pow2(m32)*pow5(
         mQ32) - 3*pow2(mU32)*pow5(mQ32) - 702*m32*mQ32*pow5(mU32) + 291*pow2(
         m32)*pow5(mU32) + 291*pow2(mQ32)*pow5(mU32) - 252*mQ32*pow6(m32) -
         1028*mU32*pow6(m32) + 18*m32*pow6(mU32) + 6*mQ32*pow6(mU32) + 128*pow7
         (m32)))/(3.*(-mQ32 + mU32)*pow2(m32 - mQ32)*pow4(m32 - mU32)) - (8*
         pow4(Xt)*(180*pow2(m32)*pow2(mQ32)*pow2(msq2)*pow2(mU32) - 300*mU32*
         pow2(mQ32)*pow2(msq2)*pow3(m32) + 420*msq2*pow2(mQ32)*pow2(mU32)*pow3(
         m32) + 60*mQ32*pow2(msq2)*pow2(mU32)*pow3(m32) + 60*mU32*pow2(m32)*
         pow2(msq2)*pow3(mQ32) - 180*m32*pow2(msq2)*pow2(mU32)*pow3(mQ32) - 180
         *msq2*mU32*pow3(m32)*pow3(mQ32) + 180*pow2(msq2)*pow3(m32)*pow3(mQ32)
         + 3202*pow2(mU32)*pow3(m32)*pow3(mQ32) - 480*msq2*pow2(m32)*pow2(mQ32)
         *pow3(mU32) - 180*mQ32*pow2(m32)*pow2(msq2)*pow3(mU32) + 180*m32*pow2(
         mQ32)*pow2(msq2)*pow3(mU32) - 60*mQ32*msq2*pow3(m32)*pow3(mU32) +
         11390*pow2(mQ32)*pow3(m32)*pow3(mU32) + 60*pow2(msq2)*pow3(m32)*pow3(
         mU32) + 300*m32*msq2*pow3(mQ32)*pow3(mU32) - 4230*pow2(m32)*pow3(mQ32)
         *pow3(mU32) - 60*pow2(msq2)*pow3(mQ32)*pow3(mU32) + 270*msq2*mU32*pow2
         (mQ32)*pow4(m32) + 180*mQ32*mU32*pow2(msq2)*pow4(m32) - 90*pow2(mQ32)*
         pow2(msq2)*pow4(m32) - 360*mQ32*msq2*pow2(mU32)*pow4(m32) - 6104*pow2(
         mQ32)*pow2(mU32)*pow4(m32) - 90*pow2(msq2)*pow2(mU32)*pow4(m32) - 120*
         msq2*pow3(mQ32)*pow4(m32) - 474*mU32*pow3(mQ32)*pow4(m32) - 13622*mQ32
         *pow3(mU32)*pow4(m32) + 210*msq2*pow3(mU32)*pow4(m32) + 30*msq2*mU32*
         pow2(m32)*pow4(mQ32) + 60*m32*mU32*pow2(msq2)*pow4(mQ32) - 90*pow2(m32
         )*pow2(msq2)*pow4(mQ32) - 120*m32*msq2*pow2(mU32)*pow4(mQ32) - 692*
         pow2(m32)*pow2(mU32)*pow4(mQ32) + 30*pow2(msq2)*pow2(mU32)*pow4(mQ32)
         + 60*msq2*pow3(m32)*pow4(mQ32) + 178*mU32*pow3(m32)*pow4(mQ32) + 722*
         m32*pow3(mU32)*pow4(mQ32) + 30*msq2*pow3(mU32)*pow4(mQ32) - 10*pow4(
         m32)*pow4(mQ32) + 480*mQ32*msq2*pow2(m32)*pow4(mU32) - 240*m32*msq2*
         pow2(mQ32)*pow4(mU32) - 7123*pow2(m32)*pow2(mQ32)*pow4(mU32) - 60*m32*
         mQ32*pow2(msq2)*pow4(mU32) + 30*pow2(m32)*pow2(msq2)*pow4(mU32) + 30*
         pow2(mQ32)*pow2(msq2)*pow4(mU32) + 9482*mQ32*pow3(m32)*pow4(mU32) -
         240*msq2*pow3(m32)*pow4(mU32) + 1861*m32*pow3(mQ32)*pow4(mU32) - 4030*
         pow4(m32)*pow4(mU32) - 190*pow4(mQ32)*pow4(mU32) - 120*mQ32*msq2*mU32*
         pow5(m32) + 60*msq2*pow2(mQ32)*pow5(m32) + 201*mU32*pow2(mQ32)*pow5(
         m32) + 6267*mQ32*pow2(mU32)*pow5(m32) + 60*msq2*pow2(mU32)*pow5(m32) -
         147*pow3(mQ32)*pow5(m32) + 5679*pow3(mU32)*pow5(m32) + 21*mU32*pow2(
         m32)*pow5(mQ32) - 42*m32*pow2(mU32)*pow5(mQ32) + 20*pow3(m32)*pow5(
         mQ32) + pow3(mU32)*pow5(mQ32) + 60*m32*mQ32*msq2*pow5(mU32) - 1079*
         mQ32*pow2(m32)*pow5(mU32) - 30*msq2*pow2(m32)*pow5(mU32) + 1015*m32*
         pow2(mQ32)*pow5(mU32) - 30*msq2*pow2(mQ32)*pow5(mU32) + 368*pow3(m32)*
         pow5(mU32) - 244*pow3(mQ32)*pow5(mU32) - 120*mQ32*mU32*pow6(m32) + 322
         *pow2(mQ32)*pow6(m32) - 2586*pow2(mU32)*pow6(m32) + 6*m32*mU32*pow6(
         mQ32) - 9*pow2(m32)*pow6(mQ32) + 3*pow2(mU32)*pow6(mQ32) - 787*m32*
         mQ32*pow6(mU32) + 392*pow2(m32)*pow6(mU32) + 347*pow2(mQ32)*pow6(mU32)
         - 176*mQ32*pow7(m32) + 176*mU32*pow7(m32) + 9*m32*pow7(mU32) + 3*mQ32
         *pow7(mU32)))/(3.*pow2(m32 - mQ32)*pow4(m32 - mU32)*pow4(-mQ32 + mU32)
         ) + (128*pow3(Xt)*(30*m3*msq2*mU32*pow2(mQ32) - 30*m3*mQ32*msq2*pow2(
         mU32) + 82*m3*pow2(mQ32)*pow2(mU32) - 30*mQ32*msq2*mU32*pow3(m3) - 117
         *mU32*pow2(mQ32)*pow3(m3) - 47*mQ32*pow2(mU32)*pow3(m3) + 30*msq2*pow2
         (mU32)*pow3(m3) + 3*m3*mU32*pow3(mQ32) - 53*m3*mQ32*pow3(mU32) + 66*
         pow3(m3)*pow3(mU32) - 3*m3*pow4(mU32) + 150*mQ32*mU32*pow5(m3) - 14*
         pow2(mQ32)*pow5(m3) - 48*pow2(mU32)*pow5(m3) + 7*mQ32*pow7(m3) - 37*
         mU32*pow7(m3) + 11*pow9(m3)))/(3.*(m32 - mQ32)*pow2(m32 - mU32)*pow3(
         mQ32 - mU32)) - (64*pow5(Xt)*(-30*m3*mU32*pow2(mQ32)*pow2(msq2) + 30*
         m3*msq2*pow2(mQ32)*pow2(mU32) + 30*m3*mQ32*pow2(msq2)*pow2(mU32) - 30*
         msq2*mU32*pow2(mQ32)*pow3(m3) + 30*pow2(mQ32)*pow2(msq2)*pow3(m3) + 60
         *mQ32*msq2*pow2(mU32)*pow3(m3) + 106*pow2(mQ32)*pow2(mU32)*pow3(m3) -
         30*pow2(msq2)*pow2(mU32)*pow3(m3) + 7*m3*pow2(mU32)*pow3(mQ32) - 5*
         mU32*pow3(m3)*pow3(mQ32) - 90*m3*mQ32*msq2*pow3(mU32) - 82*m3*pow2(
         mQ32)*pow3(mU32) + 511*mQ32*pow3(m3)*pow3(mU32) + 90*msq2*pow3(m3)*
         pow3(mU32) - 3*m3*mU32*pow4(mQ32) + 3*pow3(m3)*pow4(mQ32) - 201*m3*
         mQ32*pow4(mU32) + 215*pow3(m3)*pow4(mU32) + 30*mQ32*msq2*mU32*pow5(m3)
         - 5*mU32*pow2(mQ32)*pow5(m3) - 30*mQ32*pow2(msq2)*pow5(m3) + 30*mU32*
         pow2(msq2)*pow5(m3) - 306*mQ32*pow2(mU32)*pow5(m3) - 90*msq2*pow2(mU32
         )*pow5(m3) - 2*pow3(mQ32)*pow5(m3) - 437*pow3(mU32)*pow5(m3) - 3*m3*
         pow5(mU32) - 6*mQ32*mU32*pow7(m3) - 35*pow2(mQ32)*pow7(m3) + 179*pow2(
         mU32)*pow7(m3) + 34*mQ32*pow9(m3) + 30*mU32*pow9(m3)))/(3.*(m32 - mQ32
         )*pow3(m32 - mU32)*pow4(mQ32 - mU32)) - (32*Xt*(30*m3*mU32*pow2(mQ32)*
         pow2(msq2) - 120*m3*msq2*pow2(mQ32)*pow2(mU32) - 30*m3*mQ32*pow2(msq2)
         *pow2(mU32) + 120*msq2*mU32*pow2(mQ32)*pow3(m3) - 30*pow2(mQ32)*pow2(
         msq2)*pow3(m3) + 367*pow2(mQ32)*pow2(mU32)*pow3(m3) + 30*pow2(msq2)*
         pow2(mU32)*pow3(m3) - 16*m3*pow2(mU32)*pow3(mQ32) + 14*mU32*pow3(m3)*
         pow3(mQ32) + 120*m3*mQ32*msq2*pow3(mU32) - 137*m3*pow2(mQ32)*pow3(mU32
         ) - 226*mQ32*pow3(m3)*pow3(mU32) - 120*msq2*pow3(m3)*pow3(mU32) + 3*m3
         *mU32*pow4(mQ32) - 3*pow3(m3)*pow4(mQ32) + 160*m3*mQ32*pow4(mU32) -
         200*pow3(m3)*pow4(mU32) - 120*mQ32*msq2*mU32*pow5(m3) - 311*mU32*pow2(
         mQ32)*pow5(m3) + 30*mQ32*pow2(msq2)*pow5(m3) - 30*mU32*pow2(msq2)*pow5
         (m3) - 122*mQ32*pow2(mU32)*pow5(m3) + 120*msq2*pow2(mU32)*pow5(m3) + 2
         *pow3(mQ32)*pow5(m3) + 479*pow3(mU32)*pow5(m3) + 6*m3*pow5(mU32) + 380
         *mQ32*mU32*pow7(m3) + pow2(mQ32)*pow7(m3) - 397*pow2(mU32)*pow7(m3) -
         32*mQ32*pow9(m3) + 32*mU32*pow9(m3)))/(3.*(m32 - mQ32)*pow2(-mQ32 +
         mU32)*pow3(m32 - mU32))) + pow3(log(mQ32/MR2))*((1280*m32*(mQ32 + mU32
         )*pow2(mQ32)*pow6(Xt))/(3.*pow2(m32 - mQ32)*pow5(mQ32 - mU32)) - (4*(
         270*pow2(m32)*pow2(mQ32)*pow2(msq2)*pow2(mU32) - 210*mU32*pow2(mQ32)*
         pow2(msq2)*pow3(m32) - 180*msq2*pow2(mQ32)*pow2(mU32)*pow3(m32) - 90*
         mQ32*pow2(msq2)*pow2(mU32)*pow3(m32) - 90*mU32*pow2(m32)*pow2(msq2)*
         pow3(mQ32) - 540*msq2*pow2(m32)*pow2(mU32)*pow3(mQ32) + 90*m32*pow2(
         msq2)*pow2(mU32)*pow3(mQ32) + 540*msq2*mU32*pow3(m32)*pow3(mQ32) + 30*
         pow2(msq2)*pow3(m32)*pow3(mQ32) - 1170*pow2(mU32)*pow3(m32)*pow3(mQ32)
         + 420*msq2*pow2(m32)*pow2(mQ32)*pow3(mU32) - 90*mQ32*pow2(m32)*pow2(
         msq2)*pow3(mU32) - 150*m32*pow2(mQ32)*pow2(msq2)*pow3(mU32) - 420*mQ32
         *msq2*pow3(m32)*pow3(mU32) - 1874*pow2(mQ32)*pow3(m32)*pow3(mU32) +
         270*pow2(msq2)*pow3(m32)*pow3(mU32) + 180*m32*msq2*pow3(mQ32)*pow3(
         mU32) + 882*pow2(m32)*pow3(mQ32)*pow3(mU32) - 30*pow2(msq2)*pow3(mQ32)
         *pow3(mU32) - 180*msq2*mU32*pow2(mQ32)*pow4(m32) + 210*mQ32*mU32*pow2(
         msq2)*pow4(m32) + 60*pow2(mQ32)*pow2(msq2)*pow4(m32) + 540*mQ32*msq2*
         pow2(mU32)*pow4(m32) + 3011*pow2(mQ32)*pow2(mU32)*pow4(m32) - 270*pow2
         (msq2)*pow2(mU32)*pow4(m32) - 180*msq2*pow3(mQ32)*pow4(m32) + 213*mU32
         *pow3(mQ32)*pow4(m32) + 1383*mQ32*pow3(mU32)*pow4(m32) - 180*msq2*pow3
         (mU32)*pow4(m32) - 116*pow2(m32)*pow2(mU32)*pow4(mQ32) + 634*mU32*pow3
         (m32)*pow4(mQ32) - 71*m32*pow3(mU32)*pow4(mQ32) - 96*pow4(m32)*pow4(
         mQ32) + 120*mQ32*msq2*pow2(m32)*pow4(mU32) - 180*m32*msq2*pow2(mQ32)*
         pow4(mU32) + 466*pow2(m32)*pow2(mQ32)*pow4(mU32) + 60*m32*mQ32*pow2(
         msq2)*pow4(mU32) - 90*pow2(m32)*pow2(msq2)*pow4(mU32) + 30*pow2(mQ32)*
         pow2(msq2)*pow4(mU32) - 284*mQ32*pow3(m32)*pow4(mU32) + 60*msq2*pow3(
         m32)*pow4(mU32) - 276*m32*pow3(mQ32)*pow4(mU32) - 31*pow4(m32)*pow4(
         mU32) + 69*pow4(mQ32)*pow4(mU32) - 300*mQ32*msq2*mU32*pow5(m32) + 120*
         msq2*pow2(mQ32)*pow5(m32) - 1841*mU32*pow2(mQ32)*pow5(m32) - 90*mQ32*
         pow2(msq2)*pow5(m32) + 90*mU32*pow2(msq2)*pow5(m32) - 2539*mQ32*pow2(
         mU32)*pow5(m32) + 180*msq2*pow2(mU32)*pow5(m32) - 69*pow3(mQ32)*pow5(
         m32) - 31*pow3(mU32)*pow5(m32) - 374*mU32*pow2(m32)*pow5(mQ32) + 181*
         m32*pow2(mU32)*pow5(mQ32) + 6*pow3(m32)*pow5(mQ32) - 65*pow3(mU32)*
         pow5(mQ32) + 60*mQ32*msq2*pow6(m32) + 1900*mQ32*mU32*pow6(m32) - 60*
         msq2*mU32*pow6(m32) + 490*pow2(mQ32)*pow6(m32) + 298*pow2(mU32)*pow6(
         m32) + 47*m32*mU32*pow6(mQ32) + 38*pow2(m32)*pow6(mQ32) - pow2(mU32)*
         pow6(mQ32) - 544*mQ32*pow7(m32) - 352*mU32*pow7(m32) - 9*m32*pow7(mQ32
         ) - 3*mU32*pow7(mQ32) + 128*pow8(m32)))/(3.*(mQ32 - mU32)*pow3(m32 -
         mU32)*pow4(m32 - mQ32)) - (8*pow2(Xt)*(270*pow2(m32)*pow2(mQ32)*pow2(
         msq2)*pow2(mU32) - 210*mU32*pow2(mQ32)*pow2(msq2)*pow3(m32) - 180*msq2
         *pow2(mQ32)*pow2(mU32)*pow3(m32) - 90*mQ32*pow2(msq2)*pow2(mU32)*pow3(
         m32) - 90*mU32*pow2(m32)*pow2(msq2)*pow3(mQ32) - 540*msq2*pow2(m32)*
         pow2(mU32)*pow3(mQ32) + 90*m32*pow2(msq2)*pow2(mU32)*pow3(mQ32) + 540*
         msq2*mU32*pow3(m32)*pow3(mQ32) + 30*pow2(msq2)*pow3(m32)*pow3(mQ32) -
         2778*pow2(mU32)*pow3(m32)*pow3(mQ32) + 420*msq2*pow2(m32)*pow2(mQ32)*
         pow3(mU32) - 90*mQ32*pow2(m32)*pow2(msq2)*pow3(mU32) - 150*m32*pow2(
         mQ32)*pow2(msq2)*pow3(mU32) - 420*mQ32*msq2*pow3(m32)*pow3(mU32) -
         2826*pow2(mQ32)*pow3(m32)*pow3(mU32) + 270*pow2(msq2)*pow3(m32)*pow3(
         mU32) + 180*m32*msq2*pow3(mQ32)*pow3(mU32) + 1418*pow2(m32)*pow3(mQ32)
         *pow3(mU32) - 30*pow2(msq2)*pow3(mQ32)*pow3(mU32) - 180*msq2*mU32*pow2
         (mQ32)*pow4(m32) + 210*mQ32*mU32*pow2(msq2)*pow4(m32) + 60*pow2(mQ32)*
         pow2(msq2)*pow4(m32) + 540*mQ32*msq2*pow2(mU32)*pow4(m32) + 6299*pow2(
         mQ32)*pow2(mU32)*pow4(m32) - 270*pow2(msq2)*pow2(mU32)*pow4(m32) - 180
         *msq2*pow3(mQ32)*pow4(m32) + 1821*mU32*pow3(mQ32)*pow4(m32) + 1463*
         mQ32*pow3(mU32)*pow4(m32) - 180*msq2*pow3(mU32)*pow4(m32) - 500*pow2(
         m32)*pow2(mU32)*pow4(mQ32) + 1018*mU32*pow3(m32)*pow4(mQ32) + 57*m32*
         pow3(mU32)*pow4(mQ32) - 224*pow4(m32)*pow4(mQ32) + 120*mQ32*msq2*pow2(
         m32)*pow4(mU32) - 180*m32*msq2*pow2(mQ32)*pow4(mU32) + 394*pow2(m32)*
         pow2(mQ32)*pow4(mU32) + 60*m32*mQ32*pow2(msq2)*pow4(mU32) - 90*pow2(
         m32)*pow2(msq2)*pow4(mU32) + 30*pow2(mQ32)*pow2(msq2)*pow4(mU32) - 140
         *mQ32*pow3(m32)*pow4(mU32) + 60*msq2*pow3(m32)*pow4(mU32) - 276*m32*
         pow3(mQ32)*pow4(mU32) - 39*pow4(m32)*pow4(mU32) + 69*pow4(mQ32)*pow4(
         mU32) - 300*mQ32*msq2*mU32*pow5(m32) + 120*msq2*pow2(mQ32)*pow5(m32) -
         5273*mU32*pow2(mQ32)*pow5(m32) - 90*mQ32*pow2(msq2)*pow5(m32) + 90*
         mU32*pow2(msq2)*pow5(m32) - 3643*mQ32*pow2(mU32)*pow5(m32) + 180*msq2*
         pow2(mU32)*pow5(m32) - 605*pow3(mQ32)*pow5(m32) - 7*pow3(mU32)*pow5(
         m32) - 590*mU32*pow2(m32)*pow5(mQ32) + 397*m32*pow2(mU32)*pow5(mQ32) +
         78*pow3(m32)*pow5(mQ32) - 137*pow3(mU32)*pow5(mQ32) + 60*mQ32*msq2*
         pow6(m32) + 3292*mQ32*mU32*pow6(m32) - 60*msq2*mU32*pow6(m32) + 1658*
         pow2(mQ32)*pow6(m32) + 274*pow2(mU32)*pow6(m32) + 47*m32*mU32*pow6(
         mQ32) + 38*pow2(m32)*pow6(mQ32) - pow2(mU32)*pow6(mQ32) - 1056*mQ32*
         pow7(m32) - 344*mU32*pow7(m32) - 9*m32*pow7(mQ32) - 3*mU32*pow7(mQ32)
         + 128*pow8(m32)))/(3.*pow2(mQ32 - mU32)*pow3(m32 - mU32)*pow4(m32 -
         mQ32)) + (4*pow4(Xt)*(-300*pow2(mQ32)*pow2(msq2)*pow2(mU32)*pow3(m32)
         + 180*pow2(m32)*pow2(msq2)*pow2(mU32)*pow3(mQ32) - 180*mU32*pow2(msq2)
         *pow3(m32)*pow3(mQ32) + 360*msq2*pow2(mU32)*pow3(m32)*pow3(mQ32) + 180
         *pow2(m32)*pow2(mQ32)*pow2(msq2)*pow3(mU32) - 600*msq2*pow2(mQ32)*pow3
         (m32)*pow3(mU32) + 180*mQ32*pow2(msq2)*pow3(m32)*pow3(mU32) - 120*msq2
         *pow2(m32)*pow3(mQ32)*pow3(mU32) - 60*m32*pow2(msq2)*pow3(mQ32)*pow3(
         mU32) - 12228*pow3(m32)*pow3(mQ32)*pow3(mU32) + 270*mU32*pow2(mQ32)*
         pow2(msq2)*pow4(m32) + 360*msq2*pow2(mQ32)*pow2(mU32)*pow4(m32) - 60*
         mQ32*pow2(msq2)*pow2(mU32)*pow4(m32) - 360*msq2*mU32*pow3(mQ32)*pow4(
         m32) + 60*pow2(msq2)*pow3(mQ32)*pow4(m32) + 21176*pow2(mU32)*pow3(mQ32
         )*pow4(m32) + 360*mQ32*msq2*pow3(mU32)*pow4(m32) + 11830*pow2(mQ32)*
         pow3(mU32)*pow4(m32) - 270*pow2(msq2)*pow3(mU32)*pow4(m32) - 90*mU32*
         pow2(m32)*pow2(msq2)*pow4(mQ32) - 540*msq2*pow2(m32)*pow2(mU32)*pow4(
         mQ32) + 90*m32*pow2(msq2)*pow2(mU32)*pow4(mQ32) + 540*msq2*mU32*pow3(
         m32)*pow4(mQ32) + 30*pow2(msq2)*pow3(m32)*pow4(mQ32) - 5276*pow2(mU32)
         *pow3(m32)*pow4(mQ32) + 180*m32*msq2*pow3(mU32)*pow4(mQ32) + 3258*pow2
         (m32)*pow3(mU32)*pow4(mQ32) - 30*pow2(msq2)*pow3(mU32)*pow4(mQ32) -
         180*msq2*pow4(m32)*pow4(mQ32) + 3885*mU32*pow4(m32)*pow4(mQ32) + 540*
         msq2*pow2(m32)*pow2(mQ32)*pow4(mU32) - 180*mQ32*pow2(m32)*pow2(msq2)*
         pow4(mU32) - 90*m32*pow2(mQ32)*pow2(msq2)*pow4(mU32) - 360*mQ32*msq2*
         pow3(m32)*pow4(mU32) - 2878*pow2(mQ32)*pow3(m32)*pow4(mU32) + 270*pow2
         (msq2)*pow3(m32)*pow4(mU32) + 2708*pow2(m32)*pow3(mQ32)*pow4(mU32) +
         1016*mQ32*pow4(m32)*pow4(mU32) - 180*msq2*pow4(m32)*pow4(mU32) - 743*
         m32*pow4(mQ32)*pow4(mU32) - 180*msq2*mU32*pow2(mQ32)*pow5(m32) - 90*
         pow2(mQ32)*pow2(msq2)*pow5(m32) - 120*mQ32*msq2*pow2(mU32)*pow5(m32) -
         19808*pow2(mQ32)*pow2(mU32)*pow5(m32) + 90*pow2(msq2)*pow2(mU32)*pow5
         (m32) + 120*msq2*pow3(mQ32)*pow5(m32) - 16422*mU32*pow3(mQ32)*pow5(m32
         ) - 3938*mQ32*pow3(mU32)*pow5(m32) + 180*msq2*pow3(mU32)*pow5(m32) -
         1157*pow4(mQ32)*pow5(m32) + 13*pow4(mU32)*pow5(m32) - 2386*pow2(m32)*
         pow2(mU32)*pow5(mQ32) + 2248*mU32*pow3(m32)*pow5(mQ32) + 1030*m32*pow3
         (mU32)*pow5(mQ32) - 584*pow4(m32)*pow5(mQ32) - 140*pow4(mU32)*pow5(
         mQ32) + 120*mQ32*msq2*pow2(m32)*pow5(mU32) - 180*m32*msq2*pow2(mQ32)*
         pow5(mU32) + 210*pow2(m32)*pow2(mQ32)*pow5(mU32) + 60*m32*mQ32*pow2(
         msq2)*pow5(mU32) - 90*pow2(m32)*pow2(msq2)*pow5(mU32) + 30*pow2(mQ32)*
         pow2(msq2)*pow5(mU32) - 60*mQ32*pow3(m32)*pow5(mU32) + 60*msq2*pow3(
         m32)*pow5(mU32) - 132*m32*pow3(mQ32)*pow5(mU32) - 43*pow4(m32)*pow5(
         mU32) + 33*pow4(mQ32)*pow5(mU32) + 60*msq2*pow2(mQ32)*pow6(m32) +
         15098*mU32*pow2(mQ32)*pow6(m32) + 6078*mQ32*pow2(mU32)*pow6(m32) - 60*
         msq2*pow2(mU32)*pow6(m32) + 4730*pow3(mQ32)*pow6(m32) + 238*pow3(mU32)
         *pow6(m32) - 660*mU32*pow2(m32)*pow6(mQ32) + 552*m32*pow2(mU32)*pow6(
         mQ32) + 114*pow3(m32)*pow6(mQ32) - 174*pow3(mU32)*pow6(mQ32) - 4328*
         mQ32*mU32*pow7(m32) - 4284*pow2(mQ32)*pow7(m32) - 316*pow2(mU32)*pow7(
         m32) + 38*m32*mU32*pow7(mQ32) + 38*pow2(m32)*pow7(mQ32) - 4*pow2(mU32)
         *pow7(mQ32) + 1160*mQ32*pow8(m32) + 120*mU32*pow8(m32) - 9*m32*pow8(
         mQ32) - 3*mU32*pow8(mQ32)))/(3.*pow3(m32 - mU32)*pow4(m32 - mQ32)*pow4
         (mQ32 - mU32)) + (64*pow3(Xt)*(-21*mQ32*pow11(m3) - 23*mU32*pow11(m3)
         + 2*pow13(m3) + 30*m3*pow2(mQ32)*pow2(msq2)*pow2(mU32) - 60*mU32*pow2(
         mQ32)*pow2(msq2)*pow3(m3) - 60*msq2*pow2(mQ32)*pow2(mU32)*pow3(m3) +
         30*mQ32*pow2(msq2)*pow2(mU32)*pow3(m3) - 60*m3*msq2*pow2(mU32)*pow3(
         mQ32) + 120*msq2*mU32*pow3(m3)*pow3(mQ32) + 107*pow2(mU32)*pow3(m3)*
         pow3(mQ32) + 60*m3*msq2*pow2(mQ32)*pow3(mU32) - 30*m3*mQ32*pow2(msq2)*
         pow3(mU32) - 60*mQ32*msq2*pow3(m3)*pow3(mU32) - 167*pow2(mQ32)*pow3(m3
         )*pow3(mU32) + 30*pow2(msq2)*pow3(m3)*pow3(mU32) + 56*m3*pow3(mQ32)*
         pow3(mU32) - 93*m3*pow2(mU32)*pow4(mQ32) + 136*mU32*pow3(m3)*pow4(mQ32
         ) - 60*msq2*mU32*pow2(mQ32)*pow5(m3) + 30*mQ32*mU32*pow2(msq2)*pow5(m3
         ) + 30*pow2(mQ32)*pow2(msq2)*pow5(m3) + 120*mQ32*msq2*pow2(mU32)*pow5(
         m3) + 218*pow2(mQ32)*pow2(mU32)*pow5(m3) - 60*pow2(msq2)*pow2(mU32)*
         pow5(m3) - 60*msq2*pow3(mQ32)*pow5(m3) - 282*mU32*pow3(mQ32)*pow5(m3)
         + 118*mQ32*pow3(mU32)*pow5(m3) - 88*pow4(mQ32)*pow5(m3) + 10*m3*mU32*
         pow5(mQ32) + 8*pow3(m3)*pow5(mQ32) - 3*m3*pow6(mQ32) - 60*mQ32*msq2*
         mU32*pow7(m3) + 60*msq2*pow2(mQ32)*pow7(m3) - 35*mU32*pow2(mQ32)*pow7(
         m3) - 30*mQ32*pow2(msq2)*pow7(m3) + 30*mU32*pow2(msq2)*pow7(m3) - 225*
         mQ32*pow2(mU32)*pow7(m3) + 179*pow3(mQ32)*pow7(m3) - 23*pow3(mU32)*
         pow7(m3) + 146*mQ32*mU32*pow9(m3) - 61*pow2(mQ32)*pow9(m3) + 41*pow2(
         mU32)*pow9(m3)))/(3.*pow2(m32 - mU32)*pow3(m32 - mQ32)*pow3(mQ32 -
         mU32)) + (32*Xt*(22*mQ32*pow11(m3) - 22*mU32*pow11(m3) + 30*m3*pow2(
         mQ32)*pow2(msq2)*pow2(mU32) - 60*mU32*pow2(mQ32)*pow2(msq2)*pow3(m3) -
         60*msq2*pow2(mQ32)*pow2(mU32)*pow3(m3) + 30*mQ32*pow2(msq2)*pow2(mU32
         )*pow3(m3) - 60*m3*msq2*pow2(mU32)*pow3(mQ32) + 120*msq2*mU32*pow3(m3)
         *pow3(mQ32) + 142*pow2(mU32)*pow3(m3)*pow3(mQ32) + 60*m3*msq2*pow2(
         mQ32)*pow3(mU32) - 30*m3*mQ32*pow2(msq2)*pow3(mU32) - 60*mQ32*msq2*
         pow3(m3)*pow3(mU32) - 176*pow2(mQ32)*pow3(m3)*pow3(mU32) + 30*pow2(
         msq2)*pow3(m3)*pow3(mU32) + 47*m3*pow3(mQ32)*pow3(mU32) - 70*m3*pow2(
         mU32)*pow4(mQ32) + 90*mU32*pow3(m3)*pow4(mQ32) - 60*msq2*mU32*pow2(
         mQ32)*pow5(m3) + 30*mQ32*mU32*pow2(msq2)*pow5(m3) + 30*pow2(mQ32)*pow2
         (msq2)*pow5(m3) + 120*mQ32*msq2*pow2(mU32)*pow5(m3) + 91*pow2(mQ32)*
         pow2(mU32)*pow5(m3) - 60*pow2(msq2)*pow2(mU32)*pow5(m3) - 60*msq2*pow3
         (mQ32)*pow5(m3) - 325*mU32*pow3(mQ32)*pow5(m3) + 203*mQ32*pow3(mU32)*
         pow5(m3) - 65*pow4(mQ32)*pow5(m3) + 10*m3*mU32*pow5(mQ32) + 8*pow3(m3)
         *pow5(mQ32) - 3*m3*pow6(mQ32) - 60*mQ32*msq2*mU32*pow7(m3) + 60*msq2*
         pow2(mQ32)*pow7(m3) + 246*mU32*pow2(mQ32)*pow7(m3) - 30*mQ32*pow2(msq2
         )*pow7(m3) + 30*mU32*pow2(msq2)*pow7(m3) - 352*mQ32*pow2(mU32)*pow7(m3
         ) + 196*pow3(mQ32)*pow7(m3) - 26*pow3(mU32)*pow7(m3) + 145*mQ32*mU32*
         pow9(m3) - 206*pow2(mQ32)*pow9(m3) + 45*pow2(mU32)*pow9(m3)))/(3.*pow2
         (m32 - mU32)*pow2(mQ32 - mU32)*pow3(m32 - mQ32)) + (32*pow5(Xt)*(-32*
         mQ32*mU32*pow11(m3) + 14*pow11(m3)*pow2(mQ32) + 18*pow11(m3)*pow2(mU32
         ) + 30*pow2(mQ32)*pow2(msq2)*pow2(mU32)*pow3(m3) - 30*m3*pow2(msq2)*
         pow2(mU32)*pow3(mQ32) + 60*mU32*pow2(msq2)*pow3(m3)*pow3(mQ32) - 60*
         msq2*pow2(mU32)*pow3(m3)*pow3(mQ32) + 120*msq2*pow2(mQ32)*pow3(m3)*
         pow3(mU32) - 60*mQ32*pow2(msq2)*pow3(m3)*pow3(mU32) - 158*pow3(m3)*
         pow3(mQ32)*pow3(mU32) + 60*m3*msq2*pow2(mU32)*pow4(mQ32) - 120*msq2*
         mU32*pow3(m3)*pow4(mQ32) - 388*pow2(mU32)*pow3(m3)*pow4(mQ32) + 87*m3*
         pow3(mU32)*pow4(mQ32) - 60*m3*msq2*pow2(mQ32)*pow4(mU32) + 30*m3*mQ32*
         pow2(msq2)*pow4(mU32) + 60*mQ32*msq2*pow3(m3)*pow4(mU32) + 108*pow2(
         mQ32)*pow3(m3)*pow4(mU32) - 30*pow2(msq2)*pow3(m3)*pow4(mU32) - 31*m3*
         pow3(mQ32)*pow4(mU32) - 60*mU32*pow2(mQ32)*pow2(msq2)*pow5(m3) - 60*
         msq2*pow2(mQ32)*pow2(mU32)*pow5(m3) + 30*mQ32*pow2(msq2)*pow2(mU32)*
         pow5(m3) + 120*msq2*mU32*pow3(mQ32)*pow5(m3) - 30*pow2(msq2)*pow3(mQ32
         )*pow5(m3) + 450*pow2(mU32)*pow3(mQ32)*pow5(m3) - 120*mQ32*msq2*pow3(
         mU32)*pow5(m3) - 30*pow2(mQ32)*pow3(mU32)*pow5(m3) + 60*pow2(msq2)*
         pow3(mU32)*pow5(m3) + 60*msq2*pow4(mQ32)*pow5(m3) + 510*mU32*pow4(mQ32
         )*pow5(m3) - 83*mQ32*pow4(mU32)*pow5(m3) + 108*m3*pow2(mU32)*pow5(mQ32
         ) - 194*mU32*pow3(m3)*pow5(mQ32) + 113*pow5(m3)*pow5(mQ32) - 7*m3*mU32
         *pow6(mQ32) - 8*pow3(m3)*pow6(mQ32) + 30*pow2(mQ32)*pow2(msq2)*pow7(m3
         ) + 60*mQ32*msq2*pow2(mU32)*pow7(m3) - 182*pow2(mQ32)*pow2(mU32)*pow7(
         m3) - 30*pow2(msq2)*pow2(mU32)*pow7(m3) - 60*msq2*pow3(mQ32)*pow7(m3)
         - 362*mU32*pow3(mQ32)*pow7(m3) + 106*mQ32*pow3(mU32)*pow7(m3) - 224*
         pow4(mQ32)*pow7(m3) + 22*pow4(mU32)*pow7(m3) + 3*m3*pow7(mQ32) + 117*
         mU32*pow2(mQ32)*pow9(m3) - 6*mQ32*pow2(mU32)*pow9(m3) + 86*pow3(mQ32)*
         pow9(m3) - 37*pow3(mU32)*pow9(m3)))/(3.*pow2(m32 - mU32)*pow3(m32 -
         mQ32)*pow5(mQ32 - mU32))) + log(mU32/MR2)*(pow3(Xt)*((-128*(-3*m3*mQ32
         - 30*m3*msq2 - 3*m3*mU32 + 14*pow3(m3)))/(3.*(m32 - mQ32)*(m32 - mU32
         )*(mQ32 - mU32)) - (128*m3*(mQ32 + 3*mU32)*(3*m32*mQ32 + 3*m32*mU32 +
         5*mQ32*mU32 - 11*pow2(m32)))/(3.*(m32 - mQ32)*(m32 - mU32)*pow3(mQ32 -
         mU32)) + (256*(-17*m3*mQ32*pow2(mU32) + 26*mQ32*mU32*pow3(m3) + 8*
         pow2(mU32)*pow3(m3) + 7*mQ32*pow5(m3) - 8*mU32*pow5(m3)))/(3.*mQ32*
         pow2(m32 - mU32)*pow2(mQ32 - mU32))) - (8*pow2(Xt)*(-300*msq2*mU32*
         pow2(m32)*pow2(mQ32) - 300*mQ32*msq2*pow2(m32)*pow2(mU32) - 240*m32*
         msq2*pow2(mQ32)*pow2(mU32) - 40*pow2(m32)*pow2(mQ32)*pow2(mU32) + 360*
         mQ32*msq2*mU32*pow3(m32) + 886*mU32*pow2(mQ32)*pow3(m32) + 326*mQ32*
         pow2(mU32)*pow3(m32) + 180*m32*msq2*mU32*pow3(mQ32) - 301*mU32*pow2(
         m32)*pow3(mQ32) - 212*m32*pow2(mU32)*pow3(mQ32) + 60*msq2*pow2(mU32)*
         pow3(mQ32) + 4*pow3(m32)*pow3(mQ32) + 180*m32*mQ32*msq2*pow3(mU32) - 5
         *mQ32*pow2(m32)*pow3(mU32) - 308*m32*pow2(mQ32)*pow3(mU32) + 60*msq2*
         pow2(mQ32)*pow3(mU32) - 68*pow3(m32)*pow3(mU32) + 237*pow3(mQ32)*pow3(
         mU32) - 695*mQ32*mU32*pow4(m32) - 8*pow2(mQ32)*pow4(m32) + 136*pow2(
         mU32)*pow4(m32) + 18*m32*mU32*pow4(mQ32) + 6*pow2(mU32)*pow4(mQ32) +
         18*m32*mQ32*pow4(mU32) + 6*pow2(mQ32)*pow4(mU32) + 4*mQ32*pow5(m32) -
         4*mU32*pow5(m32)))/(3.*mQ32*(mQ32 - mU32)*mU32*pow2(m32 - mQ32)*pow2(
         m32 - mU32)) + pow2(log(msq2/MR2))*(-24*Xt*((40*m3*pow2(mQ32 - msq2))/
         (3.*(mQ32 - mU32)*pow2(m32 - mQ32)) + (40*m3*pow2(-msq2 + mU32))/(3.*(
         -mQ32 + mU32)*pow2(m32 - mU32))) - 24*((-5*(2*mQ32 - 2*msq2)*(-3*m32*
         mQ32 + 3*m32*msq2 + mQ32*msq2 - 2*pow2(m32) + pow2(mQ32)))/(6.*pow3(
         -m32 + mQ32)) - (5*(-2*msq2 + 2*mU32)*(3*m32*msq2 - 3*m32*mU32 + msq2*
         mU32 - 2*pow2(m32) + pow2(mU32)))/(6.*pow3(-m32 + mU32))) + (48*pow2(
         Xt)*((-5*(2*mQ32 - 2*msq2)*(-3*m32*mQ32 + 3*m32*msq2 + mQ32*msq2 - 2*
         pow2(m32) + pow2(mQ32)))/(6.*pow3(-m32 + mQ32)) - (5*(-2*msq2 + 2*mU32
         )*(3*m32*msq2 - 3*m32*mU32 + msq2*mU32 - 2*pow2(m32) + pow2(mU32)))/(
         6.*pow3(-m32 + mU32))))/(mQ32 - mU32) + (48*((40*m3*pow2(mQ32 - msq2))
         /(3.*(mQ32 - mU32)*pow2(m32 - mQ32)) + (40*m3*pow2(-msq2 + mU32))/(3.*
         (-mQ32 + mU32)*pow2(m32 - mU32)))*pow3(Xt))/(mQ32 - mU32) - (24*(mQ32
         + mU32)*((-5*(2*mQ32 - 2*msq2)*(-3*m32*mQ32 + 3*m32*msq2 + mQ32*msq2 -
         2*pow2(m32) + pow2(mQ32)))/(6.*pow3(-m32 + mQ32)) - (5*(-2*msq2 + 2*
         mU32)*(3*m32*msq2 - 3*m32*mU32 + msq2*mU32 - 2*pow2(m32) + pow2(mU32))
         )/(6.*pow3(-m32 + mU32)))*pow4(Xt))/pow3(mQ32 - mU32) - (24*(mQ32 +
         mU32)*((40*m3*pow2(mQ32 - msq2))/(3.*(mQ32 - mU32)*pow2(m32 - mQ32)) +
         (40*m3*pow2(-msq2 + mU32))/(3.*(-mQ32 + mU32)*pow2(m32 - mU32)))*pow5
         (Xt))/pow3(mQ32 - mU32)) + dilog(1 - m32/msq2)*((-640*(m32 - msq2)*Xt*
         (m3*mQ32*msq2 - 2*m3*mQ32*mU32 + m3*msq2*mU32 + mQ32*pow3(m3) - 2*msq2
         *pow3(m3) + mU32*pow3(m3)))/(pow2(m32 - mQ32)*pow2(m32 - mU32)) - 12*(
         2*m32 - 2*msq2)*(-10/(3.*(-m32 + mQ32)) - 10/(3.*(-m32 + mU32)) - (10*
         mQ32)/pow2(-m32 + mQ32) + (10*msq2)/pow2(-m32 + mQ32) + (10*msq2)/pow2
         (-m32 + mU32) - (10*mU32)/pow2(-m32 + mU32) - (40*mQ32*msq2)/(3.*pow3(
         -m32 + mQ32)) + (40*pow2(mQ32))/(3.*pow3(-m32 + mQ32)) - (40*msq2*mU32
         )/(3.*pow3(-m32 + mU32)) + (40*pow2(mU32))/(3.*pow3(-m32 + mU32))) + (
         24*(2*m32 - 2*msq2)*pow2(Xt)*(-10/(3.*(-m32 + mQ32)) - 10/(3.*(-m32 +
         mU32)) - (10*mQ32)/pow2(-m32 + mQ32) + (10*msq2)/pow2(-m32 + mQ32) + (
         10*msq2)/pow2(-m32 + mU32) - (10*mU32)/pow2(-m32 + mU32) - (40*mQ32*
         msq2)/(3.*pow3(-m32 + mQ32)) + (40*pow2(mQ32))/(3.*pow3(-m32 + mQ32))
         - (40*msq2*mU32)/(3.*pow3(-m32 + mU32)) + (40*pow2(mU32))/(3.*pow3(
         -m32 + mU32))))/(mQ32 - mU32) + (1280*(m32 - msq2)*(m3*mQ32*msq2 - 2*
         m3*mQ32*mU32 + m3*msq2*mU32 + mQ32*pow3(m3) - 2*msq2*pow3(m3) + mU32*
         pow3(m3))*pow3(Xt))/((mQ32 - mU32)*pow2(m32 - mQ32)*pow2(m32 - mU32))
         - (12*(2*m32 - 2*msq2)*(mQ32 + mU32)*(-10/(3.*(-m32 + mQ32)) - 10/(3.*
         (-m32 + mU32)) - (10*mQ32)/pow2(-m32 + mQ32) + (10*msq2)/pow2(-m32 +
         mQ32) + (10*msq2)/pow2(-m32 + mU32) - (10*mU32)/pow2(-m32 + mU32) - (
         40*mQ32*msq2)/(3.*pow3(-m32 + mQ32)) + (40*pow2(mQ32))/(3.*pow3(-m32 +
         mQ32)) - (40*msq2*mU32)/(3.*pow3(-m32 + mU32)) + (40*pow2(mU32))/(3.*
         pow3(-m32 + mU32)))*pow4(Xt))/pow3(mQ32 - mU32) - (640*(m32 - msq2)*(
         mQ32 + mU32)*(m3*mQ32*msq2 - 2*m3*mQ32*mU32 + m3*msq2*mU32 + mQ32*pow3
         (m3) - 2*msq2*pow3(m3) + mU32*pow3(m3))*pow5(Xt))/(pow2(m32 - mQ32)*
         pow2(m32 - mU32)*pow3(mQ32 - mU32))) + log(msq2/MR2)*(Xt*((-640*m3*
         msq2)/((m32 - mQ32)*(m32 - mU32)) + (192*m3)/(mQ32 - mU32) - (96*m3*
         mU32)/((m32 - mU32)*(-mQ32 + mU32)) + (320*(12*m3*msq2*mU32 - 6*m3*
         pow2(msq2) - m3*pow2(mU32) + mU32*pow3(m3)))/(3.*(-mQ32 + mU32)*pow2(
         m32 - mU32))) + ((1280*m3*msq2)/((m32 - mQ32)*(m32 - mU32)*(mQ32 -
         mU32)) - (192*m3*(mQ32 + 3*mU32))/pow3(mQ32 - mU32))*pow3(Xt) - (32*
         pow2(Xt)*(-60*m32*mQ32*msq2*mU32 - 75*mQ32*msq2*pow2(m32) + 310*mQ32*
         mU32*pow2(m32) - 75*msq2*mU32*pow2(m32) + 45*m32*msq2*pow2(mQ32) - 150
         *m32*mU32*pow2(mQ32) + 15*msq2*mU32*pow2(mQ32) + 82*pow2(m32)*pow2(
         mQ32) - 141*m32*mQ32*pow2(mU32) + 45*m32*msq2*pow2(mU32) + 15*mQ32*
         msq2*pow2(mU32) + 73*pow2(m32)*pow2(mU32) + 68*pow2(mQ32)*pow2(mU32) -
         169*mQ32*pow3(m32) + 90*msq2*pow3(m32) - 160*mU32*pow3(m32) + 87*pow4
         (m32)))/(3.*(mQ32 - mU32)*pow2(m32 - mQ32)*pow2(m32 - mU32)) + ((16*(
         180*m32*msq2*mU32 + 60*msq2*pow2(m32) + 2*mU32*pow2(m32) - 90*m32*pow2
         (msq2) - 30*mU32*pow2(msq2) - 3*m32*pow2(mU32) + pow3(mU32)))/(3.*pow2
         (mQ32 - mU32)*pow3(m32 - mU32)) + (24*(-3*mQ32*(mQ32 + 5*mU32)*pow2(
         mU32) - pow2(m32)*(15*mQ32*mU32 + 6*pow2(mQ32) + 29*pow2(mU32)) + (3*
         mQ32 + 7*mU32)*pow3(m32) + m32*(7*mU32*pow2(mQ32) + 31*mQ32*pow2(mU32)
         + 16*pow3(mU32)) + 4*pow4(m32)))/((m32 - mQ32)*pow2(m32 - mU32)*pow3(
         mQ32 - mU32)) + (8*(mQ32 + mU32)*(-120*m32*mQ32*msq2*mU32 - 150*mQ32*
         msq2*pow2(m32) + 440*mQ32*mU32*pow2(m32) - 150*msq2*mU32*pow2(m32) +
         90*m32*msq2*pow2(mQ32) - 219*m32*mU32*pow2(mQ32) + 30*msq2*mU32*pow2(
         mQ32) + 110*pow2(m32)*pow2(mQ32) - 219*m32*mQ32*pow2(mU32) + 90*m32*
         msq2*pow2(mU32) + 30*mQ32*msq2*pow2(mU32) + 110*pow2(m32)*pow2(mU32) +
         109*pow2(mQ32)*pow2(mU32) - 221*mQ32*pow3(m32) + 180*msq2*pow3(m32) -
         221*mU32*pow3(m32) + 111*pow4(m32)))/(3.*pow2(m32 - mQ32)*pow2(m32 -
         mU32)*pow3(mQ32 - mU32)))*pow4(Xt) + (16*(195*mQ32*msq2*mU32*pow2(m32)
         - 120*m32*msq2*mU32*pow2(mQ32) + 15*msq2*pow2(m32)*pow2(mQ32) - 215*
         mU32*pow2(m32)*pow2(mQ32) - 30*m32*mQ32*mU32*pow2(msq2) - 90*mQ32*pow2
         (m32)*pow2(msq2) + 15*mU32*pow2(m32)*pow2(msq2) + 45*m32*pow2(mQ32)*
         pow2(msq2) + 15*mU32*pow2(mQ32)*pow2(msq2) + 75*m32*mQ32*msq2*pow2(
         mU32) - 409*mQ32*pow2(m32)*pow2(mU32) + 120*msq2*pow2(m32)*pow2(mU32)
         + 197*m32*pow2(mQ32)*pow2(mU32) - 15*msq2*pow2(mQ32)*pow2(mU32) - 15*
         mQ32*msq2*pow3(m32) + 445*mQ32*mU32*pow3(m32) - 255*msq2*mU32*pow3(m32
         ) + 82*pow2(mQ32)*pow3(m32) + 45*pow2(msq2)*pow3(m32) + 212*pow2(mU32)
         *pow3(m32) + 133*m32*mQ32*pow3(mU32) - 45*m32*msq2*pow3(mU32) - 15*
         mQ32*msq2*pow3(mU32) - 69*pow2(m32)*pow3(mU32) - 64*pow2(mQ32)*pow3(
         mU32) - 169*mQ32*pow4(m32) + 60*msq2*pow4(m32) - 230*mU32*pow4(m32) +
         87*pow5(m32)))/(3.*pow2(m32 - mQ32)*pow3(m32 - mU32)) + ((-640*m3*msq2
         *(mQ32 + mU32))/((m32 - mQ32)*(m32 - mU32)*pow3(mQ32 - mU32)) - (192*
         m3*mU32)/((m32 - mU32)*pow3(-mQ32 + mU32)) - (64*(120*m3*msq2*mU32 -
         60*m3*pow2(msq2) - m3*pow2(mU32) + mU32*pow3(m3)))/(3.*pow2(m32 - mU32
         )*pow3(-mQ32 + mU32)))*pow5(Xt)) + dilog(1 - msq2/mQ32)*((-640*m3*Xt*
         pow2(mQ32 - msq2))/((mQ32 - mU32)*pow2(m32 - mQ32)) + (40*(2*mQ32 - 2*
         msq2)*(-3*m32*mQ32 + 3*m32*msq2 + mQ32*msq2 - 2*pow2(m32) + pow2(mQ32)
         ))/pow3(-m32 + mQ32) - (80*(2*mQ32 - 2*msq2)*(-3*m32*mQ32 + 3*m32*msq2
         + mQ32*msq2 - 2*pow2(m32) + pow2(mQ32))*pow2(Xt))/((mQ32 - mU32)*pow3
         (-m32 + mQ32)) + (1280*m3*pow2(mQ32 - msq2)*pow3(Xt))/(pow2(m32 - mQ32
         )*pow2(mQ32 - mU32)) + (40*(2*mQ32 - 2*msq2)*(mQ32 + mU32)*(-3*m32*
         mQ32 + 3*m32*msq2 + mQ32*msq2 - 2*pow2(m32) + pow2(mQ32))*pow4(Xt))/(
         pow3(-m32 + mQ32)*pow3(mQ32 - mU32)) - (640*m3*(mQ32 + mU32)*pow2(mQ32
         - msq2)*pow5(Xt))/(pow2(m32 - mQ32)*pow4(mQ32 - mU32))) + dilog(1 -
         msq2/mU32)*((-640*m3*Xt*pow2(-msq2 + mU32))/((-mQ32 + mU32)*pow2(m32 -
         mU32)) + (40*(-2*msq2 + 2*mU32)*(3*m32*msq2 - 3*m32*mU32 + msq2*mU32
         - 2*pow2(m32) + pow2(mU32)))/pow3(-m32 + mU32) - (160*(-msq2 + mU32)*(
         3*m32*msq2 - 3*m32*mU32 + msq2*mU32 - 2*pow2(m32) + pow2(mU32))*pow2(
         Xt))/((mQ32 - mU32)*pow3(-m32 + mU32)) - (1280*m3*pow2(-msq2 + mU32)*
         pow3(Xt))/(pow2(m32 - mU32)*pow2(-mQ32 + mU32)) + (80*(mQ32 + mU32)*(
         -msq2 + mU32)*(3*m32*msq2 - 3*m32*mU32 + msq2*mU32 - 2*pow2(m32) +
         pow2(mU32))*pow4(Xt))/(pow3(mQ32 - mU32)*pow3(-m32 + mU32)) + (640*m3*
         (mQ32 + mU32)*pow2(-msq2 + mU32)*pow5(Xt))/(pow2(m32 - mU32)*pow4(
         -mQ32 + mU32))) + dilog(1 - mQ32/mU32)*((-128*pow3(Xt)*(-(m3*mQ32*mU32
         ) + 3*m3*pow2(mQ32) + 3*m3*pow2(mU32) - 5*mQ32*pow3(m3) - 5*mU32*pow3(
         m3) + 5*pow5(m3)))/(3.*pow2(m32 - mQ32)*pow2(m32 - mU32)) + (32*Xt*(-8
         *m3*mU32*pow2(mQ32) + 8*m3*mQ32*pow2(mU32) - 10*pow2(mQ32)*pow3(m3) +
         10*pow2(mU32)*pow3(m3) + 6*m3*pow3(mQ32) - 6*m3*pow3(mU32) + 10*mQ32*
         pow5(m3) - 10*mU32*pow5(m3)))/(3.*pow2(m32 - mQ32)*pow2(m32 - mU32)) +
         (16*pow2(Xt)*(65*mU32*pow2(m32)*pow2(mQ32) + 65*mQ32*pow2(m32)*pow2(
         mU32) - 34*m32*pow2(mQ32)*pow2(mU32) - 76*mQ32*mU32*pow3(m32) + 34*
         pow2(mQ32)*pow3(m32) + 34*pow2(mU32)*pow3(m32) - 26*m32*mU32*pow3(mQ32
         ) - 29*pow2(m32)*pow3(mQ32) + 7*pow2(mU32)*pow3(mQ32) - 26*m32*mQ32*
         pow3(mU32) - 29*pow2(m32)*pow3(mU32) + 7*pow2(mQ32)*pow3(mU32) - 14*
         mQ32*pow4(m32) - 14*mU32*pow4(m32) + 9*m32*pow4(mQ32) + 3*mU32*pow4(
         mQ32) + 9*m32*pow4(mU32) + 3*mQ32*pow4(mU32) + 12*pow5(m32)))/(3.*pow3
         (m32 - mQ32)*pow3(m32 - mU32)) - (8*(mQ32 + mU32)*pow4(Xt)*(65*mU32*
         pow2(m32)*pow2(mQ32) + 65*mQ32*pow2(m32)*pow2(mU32) - 34*m32*pow2(mQ32
         )*pow2(mU32) - 76*mQ32*mU32*pow3(m32) + 34*pow2(mQ32)*pow3(m32) + 34*
         pow2(mU32)*pow3(m32) - 26*m32*mU32*pow3(mQ32) - 29*pow2(m32)*pow3(mQ32
         ) + 7*pow2(mU32)*pow3(mQ32) - 26*m32*mQ32*pow3(mU32) - 29*pow2(m32)*
         pow3(mU32) + 7*pow2(mQ32)*pow3(mU32) - 14*mQ32*pow4(m32) - 14*mU32*
         pow4(m32) + 9*m32*pow4(mQ32) + 3*mU32*pow4(mQ32) + 9*m32*pow4(mU32) +
         3*mQ32*pow4(mU32) + 12*pow5(m32)))/(3.*pow2(mQ32 - mU32)*pow3(m32 -
         mQ32)*pow3(m32 - mU32)) - (4*(220*mU32*pow2(mQ32)*pow3(m32) - 220*mQ32
         *pow2(mU32)*pow3(m32) - 188*mU32*pow2(m32)*pow3(mQ32) + 16*m32*pow2(
         mU32)*pow3(mQ32) - 68*pow3(m32)*pow3(mQ32) + 188*mQ32*pow2(m32)*pow3(
         mU32) - 16*m32*pow2(mQ32)*pow3(mU32) + 68*pow3(m32)*pow3(mU32) + 28*
         pow2(mQ32)*pow4(m32) - 28*pow2(mU32)*pow4(m32) + 70*m32*mU32*pow4(mQ32
         ) + 58*pow2(m32)*pow4(mQ32) - 8*pow2(mU32)*pow4(mQ32) - 70*m32*mQ32*
         pow4(mU32) - 58*pow2(m32)*pow4(mU32) + 8*pow2(mQ32)*pow4(mU32) - 24*
         mQ32*pow5(m32) + 24*mU32*pow5(m32) - 18*m32*pow5(mQ32) - 6*mU32*pow5(
         mQ32) + 18*m32*pow5(mU32) + 6*mQ32*pow5(mU32)))/(3.*pow3(-m32 + mQ32)*
         pow3(m32 - mU32)) + (64*(mQ32 + mU32)*(-(m3*mQ32*mU32) + 3*m3*pow2(
         mQ32) + 3*m3*pow2(mU32) - 5*mQ32*pow3(m3) - 5*mU32*pow3(m3) + 5*pow5(
         m3))*pow5(Xt))/(3.*pow2(m32 - mQ32)*pow2(m32 - mU32)*pow2(mQ32 - mU32)
         )) + (4*(540*msq2*pow2(mQ32)*pow2(mU32)*pow3(m32) - 240*msq2*pow2(m32)
         *pow2(mU32)*pow3(mQ32) + 300*msq2*mU32*pow3(m32)*pow3(mQ32) + 2217*
         pow2(mU32)*pow3(m32)*pow3(mQ32) - 840*mQ32*msq2*pow3(m32)*pow3(mU32) +
         1975*pow2(mQ32)*pow3(m32)*pow3(mU32) - 720*m32*msq2*pow3(mQ32)*pow3(
         mU32) - 2359*pow2(m32)*pow3(mQ32)*pow3(mU32) - 360*msq2*mU32*pow2(mQ32
         )*pow4(m32) + 360*mQ32*msq2*pow2(mU32)*pow4(m32) - 2127*pow2(mQ32)*
         pow2(mU32)*pow4(m32) - 704*mU32*pow3(mQ32)*pow4(m32) + 583*mQ32*pow3(
         mU32)*pow4(m32) - 180*msq2*mU32*pow2(m32)*pow4(mQ32) + 300*m32*msq2*
         pow2(mU32)*pow4(mQ32) - 541*pow2(m32)*pow2(mU32)*pow4(mQ32) + 207*mU32
         *pow3(m32)*pow4(mQ32) + 505*m32*pow3(mU32)*pow4(mQ32) + 120*msq2*pow3(
         mU32)*pow4(mQ32) - 4*pow4(m32)*pow4(mQ32) + 420*mQ32*msq2*pow2(m32)*
         pow4(mU32) + 600*m32*msq2*pow2(mQ32)*pow4(mU32) + 283*pow2(m32)*pow2(
         mQ32)*pow4(mU32) - 1259*mQ32*pow3(m32)*pow4(mU32) + 571*m32*pow3(mQ32)
         *pow4(mU32) - 60*msq2*pow3(mQ32)*pow4(mU32) + 204*pow4(m32)*pow4(mU32)
         - 183*pow4(mQ32)*pow4(mU32) + 611*mU32*pow2(mQ32)*pow5(m32) + 33*mQ32
         *pow2(mU32)*pow5(m32) + 8*pow3(mQ32)*pow5(m32) - 140*pow3(mU32)*pow5(
         m32) - 18*mU32*pow2(m32)*pow5(mQ32) + 30*m32*pow2(mU32)*pow5(mQ32) +
         12*pow3(mU32)*pow5(mQ32) - 180*m32*mQ32*msq2*pow5(mU32) + 587*mQ32*
         pow2(m32)*pow5(mU32) - 594*m32*pow2(mQ32)*pow5(mU32) - 60*msq2*pow2(
         mQ32)*pow5(mU32) - 68*pow3(m32)*pow5(mU32) + 171*pow3(mQ32)*pow5(mU32)
         - 4*pow2(mQ32)*pow6(m32) + 4*pow2(mU32)*pow6(m32)))/(3.*mQ32*mU32*(
         -mQ32 + mU32)*pow2(m32 - mQ32)*pow3(m32 - mU32)) + (4*pow4(Xt)*(600*
         msq2*pow2(mQ32)*pow2(mU32)*pow3(m32) + 480*msq2*pow2(m32)*pow2(mU32)*
         pow3(mQ32) + 300*msq2*mU32*pow3(m32)*pow3(mQ32) - 6913*pow2(mU32)*pow3
         (m32)*pow3(mQ32) - 1380*msq2*pow2(m32)*pow2(mQ32)*pow3(mU32) + 1020*
         mQ32*msq2*pow3(m32)*pow3(mU32) - 17667*pow2(mQ32)*pow3(m32)*pow3(mU32)
         + 420*m32*msq2*pow3(mQ32)*pow3(mU32) + 13435*pow2(m32)*pow3(mQ32)*
         pow3(mU32) - 360*msq2*mU32*pow2(mQ32)*pow4(m32) - 360*mQ32*msq2*pow2(
         mU32)*pow4(m32) + 10101*pow2(mQ32)*pow2(mU32)*pow4(m32) - 1222*mU32*
         pow3(mQ32)*pow4(m32) + 7283*mQ32*pow3(mU32)*pow4(m32) - 180*msq2*mU32*
         pow2(m32)*pow4(mQ32) - 240*m32*msq2*pow2(mU32)*pow4(mQ32) + 1361*pow2(
         m32)*pow2(mU32)*pow4(mQ32) + 493*mU32*pow3(m32)*pow4(mQ32) - 2665*m32*
         pow3(mU32)*pow4(mQ32) - 60*msq2*pow3(mU32)*pow4(mQ32) - 4*pow4(m32)*
         pow4(mQ32) - 360*mQ32*msq2*pow2(m32)*pow4(mU32) - 360*m32*msq2*pow2(
         mQ32)*pow4(mU32) + 6073*pow2(m32)*pow2(mQ32)*pow4(mU32) - 1881*mQ32*
         pow3(m32)*pow4(mU32) - 5659*m32*pow3(mQ32)*pow4(mU32) + 240*msq2*pow3(
         mQ32)*pow4(mU32) + 204*pow4(m32)*pow4(mU32) + 1167*pow4(mQ32)*pow4(
         mU32) + 667*mU32*pow2(mQ32)*pow5(m32) - 4793*mQ32*pow2(mU32)*pow5(m32)
         + 8*pow3(mQ32)*pow5(m32) - 140*pow3(mU32)*pow5(m32) - 18*mU32*pow2(
         m32)*pow5(mQ32) - 24*m32*pow2(mU32)*pow5(mQ32) - 6*pow3(mU32)*pow5(
         mQ32) + 180*m32*mQ32*msq2*pow5(mU32) - 975*mQ32*pow2(m32)*pow5(mU32) +
         1508*m32*pow2(mQ32)*pow5(mU32) + 60*msq2*pow2(mQ32)*pow5(mU32) - 68*
         pow3(m32)*pow5(mU32) - 417*pow3(mQ32)*pow5(mU32) + 176*mQ32*mU32*pow6(
         m32) - 4*pow2(mQ32)*pow6(m32) + 4*pow2(mU32)*pow6(m32) - 18*m32*mQ32*
         pow6(mU32) - 6*pow2(mQ32)*pow6(mU32)))/(3.*mQ32*mU32*pow2(m32 - mQ32)*
         pow3(m32 - mU32)*pow3(-mQ32 + mU32)) + (64*Xt*(-60*m3*msq2*mU32*pow2(
         mQ32) + 30*m3*mQ32*msq2*pow2(mU32) - 74*m3*pow2(mQ32)*pow2(mU32) + 30*
         msq2*pow2(mQ32)*pow3(m3) + 100*mU32*pow2(mQ32)*pow3(m3) + 87*mQ32*pow2
         (mU32)*pow3(m3) - 6*m3*mU32*pow3(mQ32) + 3*pow3(m3)*pow3(mQ32) - 112*
         mQ32*mU32*pow5(m3) - 20*pow2(mQ32)*pow5(m3) - 16*pow2(mU32)*pow5(m3) +
         22*mQ32*pow7(m3) + 16*mU32*pow7(m3)))/(3.*(m32 - mQ32)*mQ32*(-mQ32 +
         mU32)*pow2(m32 - mU32)) - (64*pow5(Xt)*(30*m3*msq2*mU32*pow2(mQ32) -
         30*m3*mQ32*msq2*pow2(mU32) + 184*m3*pow2(mQ32)*pow2(mU32) - 30*mQ32*
         msq2*mU32*pow3(m3) + 30*msq2*pow2(mQ32)*pow3(m3) - 166*mU32*pow2(mQ32)
         *pow3(m3) - 185*mQ32*pow2(mU32)*pow3(m3) + 3*m3*mU32*pow3(mQ32) + 3*
         pow3(m3)*pow3(mQ32) + 3*m3*mQ32*pow3(mU32) + 172*mQ32*mU32*pow5(m3) -
         62*pow2(mQ32)*pow5(m3) - 16*pow2(mU32)*pow5(m3) + 48*mQ32*pow7(m3) +
         16*mU32*pow7(m3)))/(3.*(m32 - mQ32)*mQ32*pow2(m32 - mU32)*pow3(mQ32 -
         mU32)) + dilog(1 - m32/mU32)*((-32*pow2(m32)*(-34*m32*mU32 + pow2(m32)
         + 17*pow2(mU32)))/(3.*pow4(m32 - mU32)) + pow3(Xt)*((64*m3*(-2*m32 +
         mQ32 + mU32)*(-34*m32*mU32 + pow2(m32) + 17*pow2(mU32)))/(3.*pow2(m32
         - mU32)*pow3(mQ32 - mU32)) + (128*(-7*m3*mQ32*mU32 + 22*m3*pow2(mQ32)
         + 3*m3*pow2(mU32) - 37*mQ32*pow3(m3) + mU32*pow3(m3) + 18*pow5(m3)))/(
         3.*pow2(m32 - mQ32)*pow2(mQ32 - mU32))) + ((1024*m3*(2*m32 - mQ32 -
         mU32)*mU32)/(3.*(m32 - mU32)*pow4(-mQ32 + mU32)) - (64*(mQ32 + mU32)*(
         -7*m3*mQ32*mU32 + 22*m3*pow2(mQ32) + 3*m3*pow2(mU32) - 37*mQ32*pow3(m3
         ) + mU32*pow3(m3) + 18*pow5(m3)))/(3.*pow2(m32 - mQ32)*pow4(mQ32 -
         mU32)))*pow5(Xt) - (8*(-417*pow2(m32)*pow2(mQ32)*pow2(mU32) + 1025*
         mU32*pow2(mQ32)*pow3(m32) + 251*mQ32*pow2(mU32)*pow3(m32) - 391*mU32*
         pow2(m32)*pow3(mQ32) + 167*m32*pow2(mU32)*pow3(mQ32) + 13*pow3(m32)*
         pow3(mQ32) + 151*mQ32*pow2(m32)*pow3(mU32) - 41*m32*pow2(mQ32)*pow3(
         mU32) - 9*pow3(m32)*pow3(mU32) + 19*pow3(mQ32)*pow3(mU32) - 902*mQ32*
         mU32*pow4(m32) - 280*pow2(mQ32)*pow4(m32) - 98*pow2(mU32)*pow4(m32) +
         34*m32*mU32*pow4(mQ32) + 37*pow2(m32)*pow4(mQ32) - 23*pow2(mU32)*pow4(
         mQ32) - 41*m32*mQ32*pow4(mU32) - 20*pow2(m32)*pow4(mU32) + pow2(mQ32)*
         pow4(mU32) + 346*mQ32*pow5(m32) + 294*mU32*pow5(m32) + 9*m32*pow5(mU32
         ) + 3*mQ32*pow5(mU32) - 128*pow6(m32)))/(3.*(mQ32 - mU32)*pow2(m32 -
         mU32)*pow3(m32 - mQ32)) - (16*pow2(Xt)*(417*pow2(m32)*pow2(mQ32)*pow2(
         mU32) - 2561*mU32*pow2(mQ32)*pow3(m32) - 251*mQ32*pow2(mU32)*pow3(m32)
         + 903*mU32*pow2(m32)*pow3(mQ32) - 167*m32*pow2(mU32)*pow3(mQ32) - 13*
         pow3(m32)*pow3(mQ32) - 151*mQ32*pow2(m32)*pow3(mU32) + 41*m32*pow2(
         mQ32)*pow3(mU32) + 9*pow3(m32)*pow3(mU32) - 19*pow3(mQ32)*pow3(mU32) +
         2438*mQ32*mU32*pow4(m32) + 280*pow2(mQ32)*pow4(m32) + 98*pow2(mU32)*
         pow4(m32) - 34*m32*mU32*pow4(mQ32) - 37*pow2(m32)*pow4(mQ32) + 23*pow2
         (mU32)*pow4(mQ32) + 41*m32*mQ32*pow4(mU32) + 20*pow2(m32)*pow4(mU32) -
         pow2(mQ32)*pow4(mU32) - 346*mQ32*pow5(m32) - 806*mU32*pow5(m32) - 9*
         m32*pow5(mU32) - 3*mQ32*pow5(mU32) + 128*pow6(m32)))/(3.*pow2(m32 -
         mU32)*pow2(-mQ32 + mU32)*pow3(m32 - mQ32)) + (8*pow4(Xt)*(-7596*pow2(
         mQ32)*pow2(mU32)*pow3(m32) + 3712*pow2(m32)*pow2(mU32)*pow3(mQ32) -
         2990*mU32*pow3(m32)*pow3(mQ32) + 2210*pow2(m32)*pow2(mQ32)*pow3(mU32)
         - 2186*mQ32*pow3(m32)*pow3(mU32) - 774*m32*pow3(mQ32)*pow3(mU32) +
         5326*mU32*pow2(mQ32)*pow4(m32) + 6466*mQ32*pow2(mU32)*pow4(m32) + 262*
         pow3(mQ32)*pow4(m32) + 746*pow3(mU32)*pow4(m32) + 526*mU32*pow2(m32)*
         pow4(mQ32) - 543*m32*pow2(mU32)*pow4(mQ32) - 3*pow3(m32)*pow4(mQ32) +
         4*pow3(mU32)*pow4(mQ32) - 29*mQ32*pow2(m32)*pow4(mU32) - 20*m32*pow2(
         mQ32)*pow4(mU32) - 25*pow3(m32)*pow4(mU32) + 14*pow3(mQ32)*pow4(mU32)
         - 4100*mQ32*mU32*pow5(m32) - 332*pow2(mQ32)*pow5(m32) - 1968*pow2(mU32
         )*pow5(m32) + 34*m32*mU32*pow5(mQ32) - 39*pow2(m32)*pow5(mQ32) - 11*
         pow2(mU32)*pow5(mQ32) + 32*m32*mQ32*pow5(mU32) + 20*pow2(m32)*pow5(
         mU32) - 4*pow2(mQ32)*pow5(mU32) + 124*mQ32*pow6(m32) + 1156*mU32*pow6(
         m32) - 9*m32*pow6(mU32) - 3*mQ32*pow6(mU32)))/(3.*pow2(m32 - mU32)*
         pow3(m32 - mQ32)*pow4(-mQ32 + mU32)) + Xt*((-64*(-7*m3*mQ32*mU32 + 22*
         m3*pow2(mQ32) + 3*m3*pow2(mU32) - 37*mQ32*pow3(m3) + mU32*pow3(m3) +
         18*pow5(m3)))/(3.*(mQ32 - mU32)*pow2(m32 - mQ32)) + (128*(17*pow2(mU32
         )*pow3(m3) - 50*mU32*pow5(m3) + pow7(m3)))/(3.*(-mQ32 + mU32)*pow3(m32
         - mU32)))) + dilog(1 - m32/mQ32)*((-64*(4*m3*mU32*pow2(mQ32) + m3*
         mQ32*pow2(mU32) + 20*mQ32*mU32*pow3(m3) - pow2(mQ32)*pow3(m3) - 11*
         pow2(mU32)*pow3(m3) - 3*m3*pow3(mQ32) - 6*m3*pow3(mU32) - 18*mQ32*pow5
         (m3) + 14*mU32*pow5(m3))*pow5(Xt))/(3.*pow2(m32 - mU32)*pow4(-mQ32 +
         mU32)) + (16*pow2(Xt)*(417*pow2(m32)*pow2(mQ32)*pow2(mU32) - 251*mU32*
         pow2(mQ32)*pow3(m32) - 2049*mQ32*pow2(mU32)*pow3(m32) - 151*mU32*pow2(
         m32)*pow3(mQ32) + 41*m32*pow2(mU32)*pow3(mQ32) + 9*pow3(m32)*pow3(mQ32
         ) + 903*mQ32*pow2(m32)*pow3(mU32) - 167*m32*pow2(mQ32)*pow3(mU32) -
         525*pow3(m32)*pow3(mU32) - 19*pow3(mQ32)*pow3(mU32) + 1414*mQ32*mU32*
         pow4(m32) + 98*pow2(mQ32)*pow4(m32) + 1304*pow2(mU32)*pow4(m32) + 41*
         m32*mU32*pow4(mQ32) + 20*pow2(m32)*pow4(mQ32) - pow2(mU32)*pow4(mQ32)
         - 34*m32*mQ32*pow4(mU32) - 37*pow2(m32)*pow4(mU32) + 23*pow2(mQ32)*
         pow4(mU32) - 294*mQ32*pow5(m32) - 858*mU32*pow5(m32) - 9*m32*pow5(mQ32
         ) - 3*mU32*pow5(mQ32) + 128*pow6(m32)))/(3.*pow2(m32 - mQ32)*pow2(mQ32
         - mU32)*pow3(m32 - mU32)) - (8*(417*pow2(m32)*pow2(mQ32)*pow2(mU32) -
         251*mU32*pow2(mQ32)*pow3(m32) - 821*mQ32*pow2(mU32)*pow3(m32) - 151*
         mU32*pow2(m32)*pow3(mQ32) + 41*m32*pow2(mU32)*pow3(mQ32) + 9*pow3(m32)
         *pow3(mQ32) + 323*mQ32*pow2(m32)*pow3(mU32) - 167*m32*pow2(mQ32)*pow3(
         mU32) - 217*pow3(m32)*pow3(mU32) - 19*pow3(mQ32)*pow3(mU32) + 762*mQ32
         *mU32*pow4(m32) + 98*pow2(mQ32)*pow4(m32) + 420*pow2(mU32)*pow4(m32) +
         41*m32*mU32*pow4(mQ32) + 20*pow2(m32)*pow4(mQ32) - pow2(mU32)*pow4(
         mQ32) - 34*m32*mQ32*pow4(mU32) + 31*pow2(m32)*pow4(mU32) + 23*pow2(
         mQ32)*pow4(mU32) - 290*mQ32*pow5(m32) - 350*mU32*pow5(m32) - 9*m32*
         pow5(mQ32) - 3*mU32*pow5(mQ32) + 128*pow6(m32)))/(3.*(mQ32 - mU32)*
         pow2(m32 - mQ32)*pow3(m32 - mU32)) - (8*pow4(Xt)*(-5516*pow2(mQ32)*
         pow2(mU32)*pow3(m32) + 882*pow2(m32)*pow2(mU32)*pow3(mQ32) - 474*mU32*
         pow3(m32)*pow3(mQ32) + 3984*pow2(m32)*pow2(mQ32)*pow3(mU32) - 6126*
         mQ32*pow3(m32)*pow3(mU32) - 502*m32*pow3(mQ32)*pow3(mU32) + 2702*mU32*
         pow2(mQ32)*pow4(m32) + 7406*mQ32*pow2(mU32)*pow4(m32) + 90*pow3(mQ32)*
         pow4(m32) + 2602*pow3(mU32)*pow4(m32) - 201*mU32*pow2(m32)*pow4(mQ32)
         + 184*m32*pow2(mU32)*pow4(mQ32) + 11*pow3(m32)*pow4(mQ32) - 54*pow3(
         mU32)*pow4(mQ32) + 1718*mQ32*pow2(m32)*pow4(mU32) - 883*m32*pow2(mQ32)
         *pow4(mU32) - 695*pow3(m32)*pow4(mU32) + 4*pow3(mQ32)*pow4(mU32) -
         3068*mQ32*mU32*pow5(m32) - 284*pow2(mQ32)*pow5(m32) - 3048*pow2(mU32)*
         pow5(m32) + 32*m32*mU32*pow5(mQ32) + 20*pow2(m32)*pow5(mQ32) - 4*pow2(
         mU32)*pow5(mQ32) - 102*m32*mQ32*pow5(mU32) - 3*pow2(m32)*pow5(mU32) +
         57*pow2(mQ32)*pow5(mU32) + 124*mQ32*pow6(m32) + 1156*mU32*pow6(m32) -
         9*m32*pow6(mQ32) - 3*mU32*pow6(mQ32)))/(3.*pow2(m32 - mQ32)*pow3(m32 -
         mU32)*pow4(-mQ32 + mU32)) - (64*pow3(Xt)*(20*m3*mU32*pow2(mQ32) - 75*
         m3*mQ32*pow2(mU32) + 110*mQ32*mU32*pow3(m3) - 2*pow2(mQ32)*pow3(m3) -
         6*pow2(mU32)*pow3(m3) - 6*m3*pow3(mQ32) + 27*m3*pow3(mU32) - 37*mQ32*
         pow5(m3) - 33*mU32*pow5(m3) + 2*pow7(m3)))/(3.*pow2(m32 - mU32)*pow3(
         -mQ32 + mU32)) - (64*Xt*(22*m3*pow2(mQ32)*pow2(mU32) - 23*mU32*pow2(
         mQ32)*pow3(m3) - 78*mQ32*pow2(mU32)*pow3(m3) - 7*m3*mU32*pow3(mQ32) -
         5*pow3(m3)*pow3(mQ32) + 3*m3*pow4(mQ32) + 135*mQ32*mU32*pow5(m3) + 19*
         pow2(mQ32)*pow5(m3) + 24*pow2(mU32)*pow5(m3) - 37*mQ32*pow7(m3) - 73*
         mU32*pow7(m3) + 20*pow9(m3)))/(3.*(-mQ32 + mU32)*pow2(m32 - mQ32)*pow2
         (m32 - mU32)))) + log(mQ32/MR2)*(pow3(Xt)*((128*(-3*m3*mQ32 - 30*m3*
         msq2 - 3*m3*mU32 + 14*pow3(m3)))/(3.*(m32 - mQ32)*(m32 - mU32)*(mQ32 -
         mU32)) + (128*m3*(3*mQ32 + mU32)*(3*m32*mQ32 + 3*m32*mU32 + 5*mQ32*
         mU32 - 11*pow2(m32)))/(3.*(m32 - mQ32)*(m32 - mU32)*pow3(mQ32 - mU32))
         - (512*(9*m3*mU32*pow2(mQ32) - 14*mQ32*mU32*pow3(m3) - 4*pow2(mQ32)*
         pow3(m3) + 4*mQ32*pow5(m3) - 3*mU32*pow5(m3)))/(3.*mU32*pow2(m32 -
         mQ32)*pow2(mQ32 - mU32))) - (8*pow2(Xt)*(300*msq2*mU32*pow2(m32)*pow2(
         mQ32) + 300*mQ32*msq2*pow2(m32)*pow2(mU32) + 240*m32*msq2*pow2(mQ32)*
         pow2(mU32) + 40*pow2(m32)*pow2(mQ32)*pow2(mU32) - 360*mQ32*msq2*mU32*
         pow3(m32) - 314*mU32*pow2(mQ32)*pow3(m32) - 898*mQ32*pow2(mU32)*pow3(
         m32) - 180*m32*msq2*mU32*pow3(mQ32) - 3*mU32*pow2(m32)*pow3(mQ32) +
         312*m32*pow2(mU32)*pow3(mQ32) - 60*msq2*pow2(mU32)*pow3(mQ32) + 72*
         pow3(m32)*pow3(mQ32) - 180*m32*mQ32*msq2*pow3(mU32) + 309*mQ32*pow2(
         m32)*pow3(mU32) + 208*m32*pow2(mQ32)*pow3(mU32) - 60*msq2*pow2(mQ32)*
         pow3(mU32) - 8*pow3(m32)*pow3(mU32) - 237*pow3(mQ32)*pow3(mU32) + 695*
         mQ32*mU32*pow4(m32) - 144*pow2(mQ32)*pow4(m32) + 16*pow2(mU32)*pow4(
         m32) - 18*m32*mU32*pow4(mQ32) - 6*pow2(mU32)*pow4(mQ32) - 18*m32*mQ32*
         pow4(mU32) - 6*pow2(mQ32)*pow4(mU32) + 8*mQ32*pow5(m32) - 8*mU32*pow5(
         m32)))/(3.*mQ32*(mQ32 - mU32)*mU32*pow2(m32 - mQ32)*pow2(m32 - mU32))
         + log(msq2/MR2)*(Xt*((-640*m3*msq2)/((m32 - mQ32)*(m32 - mU32)) - (192
         *m3)/(mQ32 - mU32) - (96*m3*mQ32)/((m32 - mQ32)*(mQ32 - mU32)) + (320*
         (12*m3*mQ32*msq2 - m3*pow2(mQ32) - 6*m3*pow2(msq2) + mQ32*pow3(m3)))/(
         3.*(mQ32 - mU32)*pow2(m32 - mQ32))) + ((-1280*m3*msq2)/((m32 - mQ32)*(
         m32 - mU32)*(mQ32 - mU32)) + (192*m3*(3*mQ32 + mU32))/pow3(mQ32 - mU32
         ))*pow3(Xt) + (32*pow2(Xt)*(-60*m32*mQ32*msq2*mU32 - 75*mQ32*msq2*pow2
         (m32) + 310*mQ32*mU32*pow2(m32) - 75*msq2*mU32*pow2(m32) + 45*m32*msq2
         *pow2(mQ32) - 141*m32*mU32*pow2(mQ32) + 15*msq2*mU32*pow2(mQ32) + 73*
         pow2(m32)*pow2(mQ32) - 150*m32*mQ32*pow2(mU32) + 45*m32*msq2*pow2(mU32
         ) + 15*mQ32*msq2*pow2(mU32) + 82*pow2(m32)*pow2(mU32) + 68*pow2(mQ32)*
         pow2(mU32) - 160*mQ32*pow3(m32) + 90*msq2*pow3(m32) - 169*mU32*pow3(
         m32) + 87*pow4(m32)))/(3.*(mQ32 - mU32)*pow2(m32 - mQ32)*pow2(m32 -
         mU32)) + ((16*(180*m32*mQ32*msq2 + 2*mQ32*pow2(m32) + 60*msq2*pow2(m32
         ) - 3*m32*pow2(mQ32) - 90*m32*pow2(msq2) - 30*mQ32*pow2(msq2) + pow3(
         mQ32)))/(3.*pow2(mQ32 - mU32)*pow3(m32 - mQ32)) - (24*(-3*mU32*(5*mQ32
         + mU32)*pow2(mQ32) - pow2(m32)*(15*mQ32*mU32 + 29*pow2(mQ32) + 6*pow2
         (mU32)) + (7*mQ32 + 3*mU32)*pow3(m32) + m32*(31*mU32*pow2(mQ32) + 7*
         mQ32*pow2(mU32) + 16*pow3(mQ32)) + 4*pow4(m32)))/((m32 - mU32)*pow2(
         m32 - mQ32)*pow3(mQ32 - mU32)) - (8*(mQ32 + mU32)*(-120*m32*mQ32*msq2*
         mU32 - 150*mQ32*msq2*pow2(m32) + 440*mQ32*mU32*pow2(m32) - 150*msq2*
         mU32*pow2(m32) + 90*m32*msq2*pow2(mQ32) - 219*m32*mU32*pow2(mQ32) + 30
         *msq2*mU32*pow2(mQ32) + 110*pow2(m32)*pow2(mQ32) - 219*m32*mQ32*pow2(
         mU32) + 90*m32*msq2*pow2(mU32) + 30*mQ32*msq2*pow2(mU32) + 110*pow2(
         m32)*pow2(mU32) + 109*pow2(mQ32)*pow2(mU32) - 221*mQ32*pow3(m32) + 180
         *msq2*pow3(m32) - 221*mU32*pow3(m32) + 111*pow4(m32)))/(3.*pow2(m32 -
         mQ32)*pow2(m32 - mU32)*pow3(mQ32 - mU32)))*pow4(Xt) + (16*(195*mQ32*
         msq2*mU32*pow2(m32) + 75*m32*msq2*mU32*pow2(mQ32) + 120*msq2*pow2(m32)
         *pow2(mQ32) - 409*mU32*pow2(m32)*pow2(mQ32) - 30*m32*mQ32*mU32*pow2(
         msq2) + 15*mQ32*pow2(m32)*pow2(msq2) - 90*mU32*pow2(m32)*pow2(msq2) -
         120*m32*mQ32*msq2*pow2(mU32) - 215*mQ32*pow2(m32)*pow2(mU32) + 15*msq2
         *pow2(m32)*pow2(mU32) + 197*m32*pow2(mQ32)*pow2(mU32) - 15*msq2*pow2(
         mQ32)*pow2(mU32) + 45*m32*pow2(msq2)*pow2(mU32) + 15*mQ32*pow2(msq2)*
         pow2(mU32) - 255*mQ32*msq2*pow3(m32) + 445*mQ32*mU32*pow3(m32) - 15*
         msq2*mU32*pow3(m32) + 212*pow2(mQ32)*pow3(m32) + 45*pow2(msq2)*pow3(
         m32) + 82*pow2(mU32)*pow3(m32) - 45*m32*msq2*pow3(mQ32) + 133*m32*mU32
         *pow3(mQ32) - 15*msq2*mU32*pow3(mQ32) - 69*pow2(m32)*pow3(mQ32) - 64*
         pow2(mU32)*pow3(mQ32) - 230*mQ32*pow4(m32) + 60*msq2*pow4(m32) - 169*
         mU32*pow4(m32) + 87*pow5(m32)))/(3.*pow2(m32 - mU32)*pow3(m32 - mQ32))
         + ((-192*m3*mQ32)/((m32 - mQ32)*pow3(mQ32 - mU32)) + (640*m3*msq2*(
         mQ32 + mU32))/((m32 - mQ32)*(m32 - mU32)*pow3(mQ32 - mU32)) - (64*(120
         *m3*mQ32*msq2 - m3*pow2(mQ32) - 60*m3*pow2(msq2) + mQ32*pow3(m3)))/(3.
         *pow2(m32 - mQ32)*pow3(mQ32 - mU32)))*pow5(Xt)) + pow2(log(msq2/MR2))*
         (-24*Xt*((40*m3*pow2(mQ32 - msq2))/(3.*(mQ32 - mU32)*pow2(m32 - mQ32))
         + (40*m3*pow2(-msq2 + mU32))/(3.*(-mQ32 + mU32)*pow2(m32 - mU32))) -
         24*((-5*(2*mQ32 - 2*msq2)*(-3*m32*mQ32 + 3*m32*msq2 + mQ32*msq2 - 2*
         pow2(m32) + pow2(mQ32)))/(6.*pow3(-m32 + mQ32)) - (5*(-2*msq2 + 2*mU32
         )*(3*m32*msq2 - 3*m32*mU32 + msq2*mU32 - 2*pow2(m32) + pow2(mU32)))/(
         6.*pow3(-m32 + mU32))) - (48*pow2(Xt)*((-5*(2*mQ32 - 2*msq2)*(-3*m32*
         mQ32 + 3*m32*msq2 + mQ32*msq2 - 2*pow2(m32) + pow2(mQ32)))/(6.*pow3(
         -m32 + mQ32)) - (5*(-2*msq2 + 2*mU32)*(3*m32*msq2 - 3*m32*mU32 + msq2*
         mU32 - 2*pow2(m32) + pow2(mU32)))/(6.*pow3(-m32 + mU32))))/(mQ32 -
         mU32) - (48*((40*m3*pow2(mQ32 - msq2))/(3.*(mQ32 - mU32)*pow2(m32 -
         mQ32)) + (40*m3*pow2(-msq2 + mU32))/(3.*(-mQ32 + mU32)*pow2(m32 - mU32
         )))*pow3(Xt))/(mQ32 - mU32) + (24*(mQ32 + mU32)*((-5*(2*mQ32 - 2*msq2)
         *(-3*m32*mQ32 + 3*m32*msq2 + mQ32*msq2 - 2*pow2(m32) + pow2(mQ32)))/(
         6.*pow3(-m32 + mQ32)) - (5*(-2*msq2 + 2*mU32)*(3*m32*msq2 - 3*m32*mU32
         + msq2*mU32 - 2*pow2(m32) + pow2(mU32)))/(6.*pow3(-m32 + mU32)))*pow4
         (Xt))/pow3(mQ32 - mU32) + (24*(mQ32 + mU32)*((40*m3*pow2(mQ32 - msq2))
         /(3.*(mQ32 - mU32)*pow2(m32 - mQ32)) + (40*m3*pow2(-msq2 + mU32))/(3.*
         (-mQ32 + mU32)*pow2(m32 - mU32)))*pow5(Xt))/pow3(mQ32 - mU32)) + dilog
         (1 - m32/msq2)*((-640*(m32 - msq2)*Xt*(m3*mQ32*msq2 - 2*m3*mQ32*mU32 +
         m3*msq2*mU32 + mQ32*pow3(m3) - 2*msq2*pow3(m3) + mU32*pow3(m3)))/(
         pow2(m32 - mQ32)*pow2(m32 - mU32)) - 12*(2*m32 - 2*msq2)*(-10/(3.*(
         -m32 + mQ32)) - 10/(3.*(-m32 + mU32)) - (10*mQ32)/pow2(-m32 + mQ32) +
         (10*msq2)/pow2(-m32 + mQ32) + (10*msq2)/pow2(-m32 + mU32) - (10*mU32)
         /pow2(-m32 + mU32) - (40*mQ32*msq2)/(3.*pow3(-m32 + mQ32)) + (40*pow2(
         mQ32))/(3.*pow3(-m32 + mQ32)) - (40*msq2*mU32)/(3.*pow3(-m32 + mU32))
         + (40*pow2(mU32))/(3.*pow3(-m32 + mU32))) - (24*(2*m32 - 2*msq2)*pow2(
         Xt)*(-10/(3.*(-m32 + mQ32)) - 10/(3.*(-m32 + mU32)) - (10*mQ32)/pow2(
         -m32 + mQ32) + (10*msq2)/pow2(-m32 + mQ32) + (10*msq2)/pow2(-m32 +
         mU32) - (10*mU32)/pow2(-m32 + mU32) - (40*mQ32*msq2)/(3.*pow3(-m32 +
         mQ32)) + (40*pow2(mQ32))/(3.*pow3(-m32 + mQ32)) - (40*msq2*mU32)/(3.*
         pow3(-m32 + mU32)) + (40*pow2(mU32))/(3.*pow3(-m32 + mU32))))/(mQ32 -
         mU32) - (1280*(m32 - msq2)*(m3*mQ32*msq2 - 2*m3*mQ32*mU32 + m3*msq2*
         mU32 + mQ32*pow3(m3) - 2*msq2*pow3(m3) + mU32*pow3(m3))*pow3(Xt))/((
         mQ32 - mU32)*pow2(m32 - mQ32)*pow2(m32 - mU32)) + (12*(2*m32 - 2*msq2)
         *(mQ32 + mU32)*(-10/(3.*(-m32 + mQ32)) - 10/(3.*(-m32 + mU32)) - (10*
         mQ32)/pow2(-m32 + mQ32) + (10*msq2)/pow2(-m32 + mQ32) + (10*msq2)/pow2
         (-m32 + mU32) - (10*mU32)/pow2(-m32 + mU32) - (40*mQ32*msq2)/(3.*pow3(
         -m32 + mQ32)) + (40*pow2(mQ32))/(3.*pow3(-m32 + mQ32)) - (40*msq2*mU32
         )/(3.*pow3(-m32 + mU32)) + (40*pow2(mU32))/(3.*pow3(-m32 + mU32)))*
         pow4(Xt))/pow3(mQ32 - mU32) + (640*(m32 - msq2)*(mQ32 + mU32)*(m3*mQ32
         *msq2 - 2*m3*mQ32*mU32 + m3*msq2*mU32 + mQ32*pow3(m3) - 2*msq2*pow3(m3
         ) + mU32*pow3(m3))*pow5(Xt))/(pow2(m32 - mQ32)*pow2(m32 - mU32)*pow3(
         mQ32 - mU32))) + dilog(1 - msq2/mQ32)*((-640*m3*Xt*pow2(mQ32 - msq2))/
         ((mQ32 - mU32)*pow2(m32 - mQ32)) + (40*(2*mQ32 - 2*msq2)*(-3*m32*mQ32
         + 3*m32*msq2 + mQ32*msq2 - 2*pow2(m32) + pow2(mQ32)))/pow3(-m32 + mQ32
         ) + (80*(2*mQ32 - 2*msq2)*(-3*m32*mQ32 + 3*m32*msq2 + mQ32*msq2 - 2*
         pow2(m32) + pow2(mQ32))*pow2(Xt))/((mQ32 - mU32)*pow3(-m32 + mQ32)) -
         (1280*m3*pow2(mQ32 - msq2)*pow3(Xt))/(pow2(m32 - mQ32)*pow2(mQ32 -
         mU32)) - (40*(2*mQ32 - 2*msq2)*(mQ32 + mU32)*(-3*m32*mQ32 + 3*m32*msq2
         + mQ32*msq2 - 2*pow2(m32) + pow2(mQ32))*pow4(Xt))/(pow3(-m32 + mQ32)*
         pow3(mQ32 - mU32)) + (640*m3*(mQ32 + mU32)*pow2(mQ32 - msq2)*pow5(Xt))
         /(pow2(m32 - mQ32)*pow4(mQ32 - mU32))) + dilog(1 - msq2/mU32)*((-640*
         m3*Xt*pow2(-msq2 + mU32))/((-mQ32 + mU32)*pow2(m32 - mU32)) + (40*(-2*
         msq2 + 2*mU32)*(3*m32*msq2 - 3*m32*mU32 + msq2*mU32 - 2*pow2(m32) +
         pow2(mU32)))/pow3(-m32 + mU32) + (160*(-msq2 + mU32)*(3*m32*msq2 - 3*
         m32*mU32 + msq2*mU32 - 2*pow2(m32) + pow2(mU32))*pow2(Xt))/((mQ32 -
         mU32)*pow3(-m32 + mU32)) + (1280*m3*pow2(-msq2 + mU32)*pow3(Xt))/(pow2
         (m32 - mU32)*pow2(-mQ32 + mU32)) - (80*(mQ32 + mU32)*(-msq2 + mU32)*(3
         *m32*msq2 - 3*m32*mU32 + msq2*mU32 - 2*pow2(m32) + pow2(mU32))*pow4(Xt
         ))/(pow3(mQ32 - mU32)*pow3(-m32 + mU32)) - (640*m3*(mQ32 + mU32)*pow2(
         -msq2 + mU32)*pow5(Xt))/(pow2(m32 - mU32)*pow4(-mQ32 + mU32))) + dilog
         (1 - mQ32/mU32)*((128*pow3(Xt)*(-(m3*mQ32*mU32) + 3*m3*pow2(mQ32) + 3*
         m3*pow2(mU32) - 5*mQ32*pow3(m3) - 5*mU32*pow3(m3) + 5*pow5(m3)))/(3.*
         pow2(m32 - mQ32)*pow2(m32 - mU32)) + (32*Xt*(-8*m3*mU32*pow2(mQ32) + 8
         *m3*mQ32*pow2(mU32) - 10*pow2(mQ32)*pow3(m3) + 10*pow2(mU32)*pow3(m3)
         + 6*m3*pow3(mQ32) - 6*m3*pow3(mU32) + 10*mQ32*pow5(m3) - 10*mU32*pow5(
         m3)))/(3.*pow2(m32 - mQ32)*pow2(m32 - mU32)) - (16*pow2(Xt)*(65*mU32*
         pow2(m32)*pow2(mQ32) + 65*mQ32*pow2(m32)*pow2(mU32) - 34*m32*pow2(mQ32
         )*pow2(mU32) - 76*mQ32*mU32*pow3(m32) + 34*pow2(mQ32)*pow3(m32) + 34*
         pow2(mU32)*pow3(m32) - 26*m32*mU32*pow3(mQ32) - 29*pow2(m32)*pow3(mQ32
         ) + 7*pow2(mU32)*pow3(mQ32) - 26*m32*mQ32*pow3(mU32) - 29*pow2(m32)*
         pow3(mU32) + 7*pow2(mQ32)*pow3(mU32) - 14*mQ32*pow4(m32) - 14*mU32*
         pow4(m32) + 9*m32*pow4(mQ32) + 3*mU32*pow4(mQ32) + 9*m32*pow4(mU32) +
         3*mQ32*pow4(mU32) + 12*pow5(m32)))/(3.*pow3(m32 - mQ32)*pow3(m32 -
         mU32)) + (8*(mQ32 + mU32)*pow4(Xt)*(65*mU32*pow2(m32)*pow2(mQ32) + 65*
         mQ32*pow2(m32)*pow2(mU32) - 34*m32*pow2(mQ32)*pow2(mU32) - 76*mQ32*
         mU32*pow3(m32) + 34*pow2(mQ32)*pow3(m32) + 34*pow2(mU32)*pow3(m32) -
         26*m32*mU32*pow3(mQ32) - 29*pow2(m32)*pow3(mQ32) + 7*pow2(mU32)*pow3(
         mQ32) - 26*m32*mQ32*pow3(mU32) - 29*pow2(m32)*pow3(mU32) + 7*pow2(mQ32
         )*pow3(mU32) - 14*mQ32*pow4(m32) - 14*mU32*pow4(m32) + 9*m32*pow4(mQ32
         ) + 3*mU32*pow4(mQ32) + 9*m32*pow4(mU32) + 3*mQ32*pow4(mU32) + 12*pow5
         (m32)))/(3.*pow2(mQ32 - mU32)*pow3(m32 - mQ32)*pow3(m32 - mU32)) - (4*
         (220*mU32*pow2(mQ32)*pow3(m32) - 220*mQ32*pow2(mU32)*pow3(m32) - 188*
         mU32*pow2(m32)*pow3(mQ32) + 16*m32*pow2(mU32)*pow3(mQ32) - 68*pow3(m32
         )*pow3(mQ32) + 188*mQ32*pow2(m32)*pow3(mU32) - 16*m32*pow2(mQ32)*pow3(
         mU32) + 68*pow3(m32)*pow3(mU32) + 28*pow2(mQ32)*pow4(m32) - 28*pow2(
         mU32)*pow4(m32) + 70*m32*mU32*pow4(mQ32) + 58*pow2(m32)*pow4(mQ32) - 8
         *pow2(mU32)*pow4(mQ32) - 70*m32*mQ32*pow4(mU32) - 58*pow2(m32)*pow4(
         mU32) + 8*pow2(mQ32)*pow4(mU32) - 24*mQ32*pow5(m32) + 24*mU32*pow5(m32
         ) - 18*m32*pow5(mQ32) - 6*mU32*pow5(mQ32) + 18*m32*pow5(mU32) + 6*mQ32
         *pow5(mU32)))/(3.*pow3(-m32 + mQ32)*pow3(m32 - mU32)) - (64*(mQ32 +
         mU32)*(-(m3*mQ32*mU32) + 3*m3*pow2(mQ32) + 3*m3*pow2(mU32) - 5*mQ32*
         pow3(m3) - 5*mU32*pow3(m3) + 5*pow5(m3))*pow5(Xt))/(3.*pow2(m32 - mQ32
         )*pow2(m32 - mU32)*pow2(mQ32 - mU32))) + (4*(540*msq2*pow2(mQ32)*pow2(
         mU32)*pow3(m32) - 840*msq2*mU32*pow3(m32)*pow3(mQ32) + 1943*pow2(mU32)
         *pow3(m32)*pow3(mQ32) - 240*msq2*pow2(m32)*pow2(mQ32)*pow3(mU32) + 300
         *mQ32*msq2*pow3(m32)*pow3(mU32) + 2281*pow2(mQ32)*pow3(m32)*pow3(mU32)
         - 720*m32*msq2*pow3(mQ32)*pow3(mU32) - 2391*pow2(m32)*pow3(mQ32)*pow3
         (mU32) + 360*msq2*mU32*pow2(mQ32)*pow4(m32) - 360*mQ32*msq2*pow2(mU32)
         *pow4(m32) - 2143*pow2(mQ32)*pow2(mU32)*pow4(m32) + 631*mU32*pow3(mQ32
         )*pow4(m32) - 744*mQ32*pow3(mU32)*pow4(m32) + 420*msq2*mU32*pow2(m32)*
         pow4(mQ32) + 600*m32*msq2*pow2(mU32)*pow4(mQ32) + 335*pow2(m32)*pow2(
         mU32)*pow4(mQ32) - 1307*mU32*pow3(m32)*pow4(mQ32) + 563*m32*pow3(mU32)
         *pow4(mQ32) - 60*msq2*pow3(mU32)*pow4(mQ32) + 216*pow4(m32)*pow4(mQ32)
         - 180*mQ32*msq2*pow2(m32)*pow4(mU32) + 300*m32*msq2*pow2(mQ32)*pow4(
         mU32) - 577*pow2(m32)*pow2(mQ32)*pow4(mU32) + 227*mQ32*pow3(m32)*pow4(
         mU32) + 533*m32*pow3(mQ32)*pow4(mU32) + 120*msq2*pow3(mQ32)*pow4(mU32)
         - 8*pow4(m32)*pow4(mU32) - 191*pow4(mQ32)*pow4(mU32) + 17*mU32*pow2(
         mQ32)*pow5(m32) + 631*mQ32*pow2(mU32)*pow5(m32) - 152*pow3(mQ32)*pow5(
         m32) + 16*pow3(mU32)*pow5(m32) - 180*m32*msq2*mU32*pow5(mQ32) + 603*
         mU32*pow2(m32)*pow5(mQ32) - 614*m32*pow2(mU32)*pow5(mQ32) - 60*msq2*
         pow2(mU32)*pow5(mQ32) - 72*pow3(m32)*pow5(mQ32) + 179*pow3(mU32)*pow5(
         mQ32) - 18*mQ32*pow2(m32)*pow5(mU32) + 30*m32*pow2(mQ32)*pow5(mU32) +
         12*pow3(mQ32)*pow5(mU32) + 8*pow2(mQ32)*pow6(m32) - 8*pow2(mU32)*pow6(
         m32)))/(3.*mQ32*(mQ32 - mU32)*mU32*pow2(m32 - mU32)*pow3(m32 - mQ32))
         + (4*pow4(Xt)*(600*msq2*pow2(mQ32)*pow2(mU32)*pow3(m32) - 1380*msq2*
         pow2(m32)*pow2(mU32)*pow3(mQ32) + 1020*msq2*mU32*pow3(m32)*pow3(mQ32)
         - 17603*pow2(mU32)*pow3(m32)*pow3(mQ32) + 480*msq2*pow2(m32)*pow2(mQ32
         )*pow3(mU32) + 300*mQ32*msq2*pow3(m32)*pow3(mU32) - 7009*pow2(mQ32)*
         pow3(m32)*pow3(mU32) + 420*m32*msq2*pow3(mQ32)*pow3(mU32) + 13499*pow2
         (m32)*pow3(mQ32)*pow3(mU32) - 360*msq2*mU32*pow2(mQ32)*pow4(m32) - 360
         *mQ32*msq2*pow2(mU32)*pow4(m32) + 10117*pow2(mQ32)*pow2(mU32)*pow4(m32
         ) + 7235*mU32*pow3(mQ32)*pow4(m32) - 1198*mQ32*pow3(mU32)*pow4(m32) -
         360*msq2*mU32*pow2(m32)*pow4(mQ32) - 360*m32*msq2*pow2(mU32)*pow4(mQ32
         ) + 5965*pow2(m32)*pow2(mU32)*pow4(mQ32) - 1833*mU32*pow3(m32)*pow4(
         mQ32) - 5635*m32*pow3(mU32)*pow4(mQ32) + 240*msq2*pow3(mU32)*pow4(mQ32
         ) + 216*pow4(m32)*pow4(mQ32) - 180*mQ32*msq2*pow2(m32)*pow4(mU32) -
         240*m32*msq2*pow2(mQ32)*pow4(mU32) + 1421*pow2(m32)*pow2(mQ32)*pow4(
         mU32) + 481*mQ32*pow3(m32)*pow4(mU32) - 2733*m32*pow3(mQ32)*pow4(mU32)
         - 60*msq2*pow3(mQ32)*pow4(mU32) - 8*pow4(m32)*pow4(mU32) + 1191*pow4(
         mQ32)*pow4(mU32) - 4777*mU32*pow2(mQ32)*pow5(m32) + 655*mQ32*pow2(mU32
         )*pow5(m32) - 152*pow3(mQ32)*pow5(m32) + 16*pow3(mU32)*pow5(m32) + 180
         *m32*msq2*mU32*pow5(mQ32) - 991*mU32*pow2(m32)*pow5(mQ32) + 1552*m32*
         pow2(mU32)*pow5(mQ32) + 60*msq2*pow2(mU32)*pow5(mQ32) - 72*pow3(m32)*
         pow5(mQ32) - 441*pow3(mU32)*pow5(mQ32) - 18*mQ32*pow2(m32)*pow5(mU32)
         - 24*m32*pow2(mQ32)*pow5(mU32) - 6*pow3(mQ32)*pow5(mU32) + 176*mQ32*
         mU32*pow6(m32) + 8*pow2(mQ32)*pow6(m32) - 8*pow2(mU32)*pow6(m32) - 18*
         m32*mU32*pow6(mQ32) - 6*pow2(mU32)*pow6(mQ32)))/(3.*mQ32*mU32*pow2(m32
         - mU32)*pow3(m32 - mQ32)*pow3(mQ32 - mU32)) + (64*Xt*(30*m3*msq2*mU32
         *pow2(mQ32) - 60*m3*mQ32*msq2*pow2(mU32) - 74*m3*pow2(mQ32)*pow2(mU32)
         + 87*mU32*pow2(mQ32)*pow3(m3) + 100*mQ32*pow2(mU32)*pow3(m3) + 30*
         msq2*pow2(mU32)*pow3(m3) - 6*m3*mQ32*pow3(mU32) + 3*pow3(m3)*pow3(mU32
         ) - 112*mQ32*mU32*pow5(m3) - 16*pow2(mQ32)*pow5(m3) - 20*pow2(mU32)*
         pow5(m3) + 16*mQ32*pow7(m3) + 22*mU32*pow7(m3)))/(3.*(m32 - mU32)*(
         mQ32 - mU32)*mU32*pow2(m32 - mQ32)) + (64*pow5(Xt)*(-30*m3*msq2*mU32*
         pow2(mQ32) + 30*m3*mQ32*msq2*pow2(mU32) + 184*m3*pow2(mQ32)*pow2(mU32)
         - 30*mQ32*msq2*mU32*pow3(m3) - 185*mU32*pow2(mQ32)*pow3(m3) - 166*
         mQ32*pow2(mU32)*pow3(m3) + 30*msq2*pow2(mU32)*pow3(m3) + 3*m3*mU32*
         pow3(mQ32) + 3*m3*mQ32*pow3(mU32) + 3*pow3(m3)*pow3(mU32) + 172*mQ32*
         mU32*pow5(m3) - 16*pow2(mQ32)*pow5(m3) - 62*pow2(mU32)*pow5(m3) + 16*
         mQ32*pow7(m3) + 48*mU32*pow7(m3)))/(3.*(m32 - mU32)*mU32*pow2(m32 -
         mQ32)*pow3(mQ32 - mU32)) + dilog(1 - m32/mQ32)*((-64*pow2(m32)*(-18*
         m32*mQ32 + pow2(m32) + 9*pow2(mQ32)))/(3.*pow4(m32 - mQ32)) + pow3(Xt)
         *((128*m3*(2*m32 - mQ32 - mU32)*(-18*m32*mQ32 + pow2(m32) + 9*pow2(
         mQ32)))/(3.*pow2(m32 - mQ32)*pow3(mQ32 - mU32)) + (128*(-7*m3*mQ32*
         mU32 + 3*m3*pow2(mQ32) + 22*m3*pow2(mU32) + mQ32*pow3(m3) - 37*mU32*
         pow3(m3) + 18*pow5(m3)))/(3.*pow2(m32 - mU32)*pow2(-mQ32 + mU32))) + (
         (-1024*m3*mQ32*(-2*m32 + mQ32 + mU32))/(3.*(m32 - mQ32)*pow4(mQ32 -
         mU32)) - (64*(mQ32 + mU32)*(-7*m3*mQ32*mU32 + 3*m3*pow2(mQ32) + 22*m3*
         pow2(mU32) + mQ32*pow3(m3) - 37*mU32*pow3(m3) + 18*pow5(m3)))/(3.*pow2
         (m32 - mU32)*pow4(-mQ32 + mU32)))*pow5(Xt) + (8*(-417*pow2(m32)*pow2(
         mQ32)*pow2(mU32) + 251*mU32*pow2(mQ32)*pow3(m32) + 1025*mQ32*pow2(mU32
         )*pow3(m32) + 151*mU32*pow2(m32)*pow3(mQ32) - 41*m32*pow2(mU32)*pow3(
         mQ32) - 9*pow3(m32)*pow3(mQ32) - 391*mQ32*pow2(m32)*pow3(mU32) + 167*
         m32*pow2(mQ32)*pow3(mU32) + 13*pow3(m32)*pow3(mU32) + 19*pow3(mQ32)*
         pow3(mU32) - 902*mQ32*mU32*pow4(m32) - 98*pow2(mQ32)*pow4(m32) - 280*
         pow2(mU32)*pow4(m32) - 41*m32*mU32*pow4(mQ32) - 20*pow2(m32)*pow4(mQ32
         ) + pow2(mU32)*pow4(mQ32) + 34*m32*mQ32*pow4(mU32) + 37*pow2(m32)*pow4
         (mU32) - 23*pow2(mQ32)*pow4(mU32) + 294*mQ32*pow5(m32) + 346*mU32*pow5
         (m32) + 9*m32*pow5(mQ32) + 3*mU32*pow5(mQ32) - 128*pow6(m32)))/(3.*(
         mQ32 - mU32)*pow2(m32 - mQ32)*pow3(m32 - mU32)) - (16*pow2(Xt)*(417*
         pow2(m32)*pow2(mQ32)*pow2(mU32) - 251*mU32*pow2(mQ32)*pow3(m32) - 2561
         *mQ32*pow2(mU32)*pow3(m32) - 151*mU32*pow2(m32)*pow3(mQ32) + 41*m32*
         pow2(mU32)*pow3(mQ32) + 9*pow3(m32)*pow3(mQ32) + 903*mQ32*pow2(m32)*
         pow3(mU32) - 167*m32*pow2(mQ32)*pow3(mU32) - 13*pow3(m32)*pow3(mU32) -
         19*pow3(mQ32)*pow3(mU32) + 2438*mQ32*mU32*pow4(m32) + 98*pow2(mQ32)*
         pow4(m32) + 280*pow2(mU32)*pow4(m32) + 41*m32*mU32*pow4(mQ32) + 20*
         pow2(m32)*pow4(mQ32) - pow2(mU32)*pow4(mQ32) - 34*m32*mQ32*pow4(mU32)
         - 37*pow2(m32)*pow4(mU32) + 23*pow2(mQ32)*pow4(mU32) - 806*mQ32*pow5(
         m32) - 346*mU32*pow5(m32) - 9*m32*pow5(mQ32) - 3*mU32*pow5(mQ32) + 128
         *pow6(m32)))/(3.*pow2(m32 - mQ32)*pow2(mQ32 - mU32)*pow3(m32 - mU32))
         + (8*pow4(Xt)*(-7612*pow2(mQ32)*pow2(mU32)*pow3(m32) + 2234*pow2(m32)*
         pow2(mU32)*pow3(mQ32) - 2210*mU32*pow3(m32)*pow3(mQ32) + 3704*pow2(m32
         )*pow2(mQ32)*pow3(mU32) - 2958*mQ32*pow3(m32)*pow3(mU32) - 782*m32*
         pow3(mQ32)*pow3(mU32) + 6492*mU32*pow2(mQ32)*pow4(m32) + 5310*mQ32*
         pow2(mU32)*pow4(m32) + 754*pow3(mQ32)*pow4(m32) + 244*pow3(mU32)*pow4(
         m32) - 23*mU32*pow2(m32)*pow4(mQ32) - 26*m32*pow2(mU32)*pow4(mQ32) -
         27*pow3(m32)*pow4(mQ32) + 16*pow3(mU32)*pow4(mQ32) + 506*mQ32*pow2(m32
         )*pow4(mU32) - 533*m32*pow2(mQ32)*pow4(mU32) + 7*pow3(m32)*pow4(mU32)
         + 4*pow3(mQ32)*pow4(mU32) - 4104*mQ32*mU32*pow5(m32) - 1978*pow2(mQ32)
         *pow5(m32) - 318*pow2(mU32)*pow5(m32) + 32*m32*mU32*pow5(mQ32) + 20*
         pow2(m32)*pow5(mQ32) - 4*pow2(mU32)*pow5(mQ32) + 38*m32*mQ32*pow5(mU32
         ) - 41*pow2(m32)*pow5(mU32) - 13*pow2(mQ32)*pow5(mU32) + 1160*mQ32*
         pow6(m32) + 120*mU32*pow6(m32) - 9*m32*pow6(mQ32) - 3*mU32*pow6(mQ32))
         )/(3.*pow2(m32 - mQ32)*pow3(m32 - mU32)*pow4(mQ32 - mU32)) + Xt*((-64*
         (-7*m3*mQ32*mU32 + 3*m3*pow2(mQ32) + 22*m3*pow2(mU32) + mQ32*pow3(m3)
         - 37*mU32*pow3(m3) + 18*pow5(m3)))/(3.*(-mQ32 + mU32)*pow2(m32 - mU32)
         ) + (256*(9*pow2(mQ32)*pow3(m3) - 26*mQ32*pow5(m3) + pow7(m3)))/(3.*(
         mQ32 - mU32)*pow3(m32 - mQ32)))) + dilog(1 - m32/mU32)*((-64*(m3*mU32*
         pow2(mQ32) + 4*m3*mQ32*pow2(mU32) + 20*mQ32*mU32*pow3(m3) - 11*pow2(
         mQ32)*pow3(m3) - pow2(mU32)*pow3(m3) - 6*m3*pow3(mQ32) - 3*m3*pow3(
         mU32) + 14*mQ32*pow5(m3) - 18*mU32*pow5(m3))*pow5(Xt))/(3.*pow2(m32 -
         mQ32)*pow4(mQ32 - mU32)) + (16*pow2(Xt)*(417*pow2(m32)*pow2(mQ32)*pow2
         (mU32) - 2049*mU32*pow2(mQ32)*pow3(m32) - 251*mQ32*pow2(mU32)*pow3(m32
         ) + 903*mU32*pow2(m32)*pow3(mQ32) - 167*m32*pow2(mU32)*pow3(mQ32) -
         525*pow3(m32)*pow3(mQ32) - 151*mQ32*pow2(m32)*pow3(mU32) + 41*m32*pow2
         (mQ32)*pow3(mU32) + 9*pow3(m32)*pow3(mU32) - 19*pow3(mQ32)*pow3(mU32)
         + 1414*mQ32*mU32*pow4(m32) + 1304*pow2(mQ32)*pow4(m32) + 98*pow2(mU32)
         *pow4(m32) - 34*m32*mU32*pow4(mQ32) - 37*pow2(m32)*pow4(mQ32) + 23*
         pow2(mU32)*pow4(mQ32) + 41*m32*mQ32*pow4(mU32) + 20*pow2(m32)*pow4(
         mU32) - pow2(mQ32)*pow4(mU32) - 858*mQ32*pow5(m32) - 294*mU32*pow5(m32
         ) - 9*m32*pow5(mU32) - 3*mQ32*pow5(mU32) + 128*pow6(m32)))/(3.*pow2(
         m32 - mU32)*pow2(mQ32 - mU32)*pow3(m32 - mQ32)) + (8*(417*pow2(m32)*
         pow2(mQ32)*pow2(mU32) - 809*mU32*pow2(mQ32)*pow3(m32) - 251*mQ32*pow2(
         mU32)*pow3(m32) + 319*mU32*pow2(m32)*pow3(mQ32) - 167*m32*pow2(mU32)*
         pow3(mQ32) - 229*pow3(m32)*pow3(mQ32) - 151*mQ32*pow2(m32)*pow3(mU32)
         + 41*m32*pow2(mQ32)*pow3(mU32) + 9*pow3(m32)*pow3(mU32) - 19*pow3(mQ32
         )*pow3(mU32) + 750*mQ32*mU32*pow4(m32) + 432*pow2(mQ32)*pow4(m32) + 98
         *pow2(mU32)*pow4(m32) - 34*m32*mU32*pow4(mQ32) + 35*pow2(m32)*pow4(
         mQ32) + 23*pow2(mU32)*pow4(mQ32) + 41*m32*mQ32*pow4(mU32) + 20*pow2(
         m32)*pow4(mU32) - pow2(mQ32)*pow4(mU32) - 354*mQ32*pow5(m32) - 286*
         mU32*pow5(m32) - 9*m32*pow5(mU32) - 3*mQ32*pow5(mU32) + 128*pow6(m32))
         )/(3.*(mQ32 - mU32)*pow2(m32 - mU32)*pow3(m32 - mQ32)) - (8*pow4(Xt)*(
         -5500*pow2(mQ32)*pow2(mU32)*pow3(m32) + 3992*pow2(m32)*pow2(mU32)*pow3
         (mQ32) - 6158*mU32*pow3(m32)*pow3(mQ32) + 858*pow2(m32)*pow2(mQ32)*
         pow3(mU32) - 450*mQ32*pow3(m32)*pow3(mU32) - 494*m32*pow3(mQ32)*pow3(
         mU32) + 7422*mU32*pow2(mQ32)*pow4(m32) + 2676*mQ32*pow2(mU32)*pow4(m32
         ) + 2620*pow3(mQ32)*pow4(m32) + 82*pow3(mU32)*pow4(m32) + 1738*mU32*
         pow2(m32)*pow4(mQ32) - 893*m32*pow2(mU32)*pow4(mQ32) - 705*pow3(m32)*
         pow4(mQ32) + 4*pow3(mU32)*pow4(mQ32) - 207*mQ32*pow2(m32)*pow4(mU32) +
         190*m32*pow2(mQ32)*pow4(mU32) + 13*pow3(m32)*pow4(mU32) - 56*pow3(
         mQ32)*pow4(mU32) - 3064*mQ32*mU32*pow5(m32) - 3062*pow2(mQ32)*pow5(m32
         ) - 274*pow2(mU32)*pow5(m32) - 106*m32*mU32*pow5(mQ32) - pow2(m32)*
         pow5(mQ32) + 59*pow2(mU32)*pow5(mQ32) + 32*m32*mQ32*pow5(mU32) + 20*
         pow2(m32)*pow5(mU32) - 4*pow2(mQ32)*pow5(mU32) + 1160*mQ32*pow6(m32) +
         120*mU32*pow6(m32) - 9*m32*pow6(mU32) - 3*mQ32*pow6(mU32)))/(3.*pow2(
         m32 - mU32)*pow3(m32 - mQ32)*pow4(mQ32 - mU32)) - (128*pow3(Xt)*(-38*
         m3*mU32*pow2(mQ32) + 10*m3*mQ32*pow2(mU32) + 56*mQ32*mU32*pow3(m3) -
         pow2(mQ32)*pow3(m3) - pow2(mU32)*pow3(m3) + 13*m3*pow3(mQ32) - 3*m3*
         pow3(mU32) - 19*mQ32*pow5(m3) - 19*mU32*pow5(m3) + 2*pow7(m3)))/(3.*
         pow2(m32 - mQ32)*pow3(mQ32 - mU32)) - (64*Xt*(22*m3*pow2(mQ32)*pow2(
         mU32) - 80*mU32*pow2(mQ32)*pow3(m3) - 23*mQ32*pow2(mU32)*pow3(m3) - 7*
         m3*mQ32*pow3(mU32) - 5*pow3(m3)*pow3(mU32) + 3*m3*pow4(mU32) + 139*
         mQ32*mU32*pow5(m3) + 26*pow2(mQ32)*pow5(m3) + 19*pow2(mU32)*pow5(m3) -
         77*mQ32*pow7(m3) - 39*mU32*pow7(m3) + 22*pow9(m3)))/(3.*(mQ32 - mU32)
         *pow2(m32 - mQ32)*pow2(m32 - mU32))) + log(mU32/MR2)*((-8*(-180*msq2*
         mU32*pow2(m32)*pow2(mQ32) - 180*mQ32*msq2*pow2(m32)*pow2(mU32) - 180*
         m32*msq2*pow2(mQ32)*pow2(mU32) - 2035*pow2(m32)*pow2(mQ32)*pow2(mU32)
         + 540*mQ32*msq2*mU32*pow3(m32) - 30*msq2*pow2(mQ32)*pow3(m32) + 1528*
         mU32*pow2(mQ32)*pow3(m32) + 1626*mQ32*pow2(mU32)*pow3(m32) - 30*msq2*
         pow2(mU32)*pow3(m32) + 90*m32*msq2*mU32*pow3(mQ32) - 703*mU32*pow2(m32
         )*pow3(mQ32) + 826*m32*pow2(mU32)*pow3(mQ32) + 30*msq2*pow2(mU32)*pow3
         (mQ32) + 198*pow3(m32)*pow3(mQ32) + 90*m32*mQ32*msq2*pow3(mU32) - 793*
         mQ32*pow2(m32)*pow3(mU32) + 834*m32*pow2(mQ32)*pow3(mU32) + 30*msq2*
         pow2(mQ32)*pow3(mU32) + 160*pow3(m32)*pow3(mU32) - 301*pow3(mQ32)*pow3
         (mU32) - 90*mQ32*msq2*pow4(m32) - 734*mQ32*mU32*pow4(m32) - 90*msq2*
         mU32*pow4(m32) - 219*pow2(mQ32)*pow4(m32) - 193*pow2(mU32)*pow4(m32) -
         32*m32*mU32*pow4(mQ32) - 20*pow2(m32)*pow4(mQ32) + 4*pow2(mU32)*pow4(
         mQ32) + 3*m32*mQ32*pow4(mU32) + 9*pow2(m32)*pow4(mU32) - 246*mQ32*pow5
         (m32) - 266*mU32*pow5(m32) + 9*m32*pow5(mQ32) + 3*mU32*pow5(mQ32) +
         342*pow6(m32)))/(3.*pow3(m32 - mQ32)*pow3(m32 - mU32)) + (16*pow2(Xt)*
         (-720*msq2*pow2(m32)*pow2(mQ32)*pow2(mU32) + 30*msq2*mU32*pow2(mQ32)*
         pow3(m32) + 30*mQ32*msq2*pow2(mU32)*pow3(m32) - 6424*pow2(mQ32)*pow2(
         mU32)*pow3(m32) + 360*msq2*mU32*pow2(m32)*pow3(mQ32) + 90*m32*msq2*
         pow2(mU32)*pow3(mQ32) + 2394*pow2(m32)*pow2(mU32)*pow3(mQ32) - 30*msq2
         *pow3(m32)*pow3(mQ32) - 1226*mU32*pow3(m32)*pow3(mQ32) + 360*mQ32*msq2
         *pow2(m32)*pow3(mU32) + 90*m32*msq2*pow2(mQ32)*pow3(mU32) + 2450*pow2(
         m32)*pow2(mQ32)*pow3(mU32) - 1306*mQ32*pow3(m32)*pow3(mU32) - 30*msq2*
         pow3(m32)*pow3(mU32) - 1168*m32*pow3(mQ32)*pow3(mU32) + 60*msq2*pow3(
         mQ32)*pow3(mU32) + 180*mQ32*msq2*mU32*pow4(m32) - 90*msq2*pow2(mQ32)*
         pow4(m32) + 4594*mU32*pow2(mQ32)*pow4(m32) + 4658*mQ32*pow2(mU32)*pow4
         (m32) - 90*msq2*pow2(mU32)*pow4(m32) - 310*pow3(mQ32)*pow4(m32) - 278*
         pow3(mU32)*pow4(m32) - 90*m32*msq2*mU32*pow4(mQ32) - 353*mU32*pow2(m32
         )*pow4(mQ32) + 156*m32*pow2(mU32)*pow4(mQ32) - 30*msq2*pow2(mU32)*pow4
         (mQ32) + 275*pow3(m32)*pow4(mQ32) + 22*pow3(mU32)*pow4(mQ32) - 90*m32*
         mQ32*msq2*pow4(mU32) - 325*mQ32*pow2(m32)*pow4(mU32) + 130*m32*pow2(
         mQ32)*pow4(mU32) - 30*msq2*pow2(mQ32)*pow4(mU32) + 265*pow3(m32)*pow4(
         mU32) + 30*pow3(mQ32)*pow4(mU32) - 3694*mQ32*mU32*pow5(m32) - 384*pow2
         (mQ32)*pow5(m32) - 418*pow2(mU32)*pow5(m32) - 15*m32*mU32*pow5(mQ32) +
         9*pow2(m32)*pow5(mQ32) - 6*pow2(mU32)*pow5(mQ32) - 15*m32*mQ32*pow5(
         mU32) + 9*pow2(m32)*pow5(mU32) - 6*pow2(mQ32)*pow5(mU32) + 462*mQ32*
         pow6(m32) + 474*mU32*pow6(m32)))/(3.*pow2(mQ32 - mU32)*pow3(m32 - mQ32
         )*pow3(m32 - mU32)) + (5120*m32*mQ32*mU32*pow6(Xt))/(3.*(m32 - mQ32)*(
         m32 - mU32)*pow4(mQ32 - mU32)) + log(msq2/MR2)*((24*(4*m32*mQ32*mU32*(
         mQ32 + mU32) - 2*pow2(mQ32)*pow2(mU32) - pow2(m32)*(8*mQ32*mU32 + pow2
         (mQ32) + pow2(mU32)) + 2*(mQ32 + mU32)*pow3(m32)))/(pow2(m32 - mQ32)*
         pow2(m32 - mU32)) + (8*(180*m32*mQ32*msq2 + 2*mQ32*pow2(m32) + 60*msq2
         *pow2(m32) - 3*m32*pow2(mQ32) - 90*m32*pow2(msq2) - 30*mQ32*pow2(msq2)
         + pow3(mQ32)))/(3.*pow3(m32 - mQ32)) + (8*(180*m32*msq2*mU32 + 60*
         msq2*pow2(m32) + 2*mU32*pow2(m32) - 90*m32*pow2(msq2) - 30*mU32*pow2(
         msq2) - 3*m32*pow2(mU32) + pow3(mU32)))/(3.*pow3(m32 - mU32)) + ((64*(
         120*m3*mQ32*msq2 - m3*pow2(mQ32) - 60*m3*pow2(msq2) + mQ32*pow3(m3)))/
         (3.*pow2(m32 - mQ32)*pow2(mQ32 - mU32)) + (64*(120*m3*msq2*mU32 - 60*
         m3*pow2(msq2) - m3*pow2(mU32) + mU32*pow3(m3)))/(3.*pow2(m32 - mU32)*
         pow2(-mQ32 + mU32)) + (192*(-2*m3*mQ32*mU32 + (mQ32 + mU32)*pow3(m3)))
         /((m32 - mQ32)*(m32 - mU32)*pow2(mQ32 - mU32)))*pow3(Xt) + pow2(Xt)*((
         -16*(180*m32*mQ32*msq2 + 2*mQ32*pow2(m32) + 60*msq2*pow2(m32) - 3*m32*
         pow2(mQ32) - 90*m32*pow2(msq2) - 30*mQ32*pow2(msq2) + pow3(mQ32)))/(3.
         *(mQ32 - mU32)*pow3(m32 - mQ32)) + (16*(180*m32*msq2*mU32 + 60*msq2*
         pow2(m32) + 2*mU32*pow2(m32) - 90*m32*pow2(msq2) - 30*mU32*pow2(msq2)
         - 3*m32*pow2(mU32) + pow3(mU32)))/(3.*(mQ32 - mU32)*pow3(m32 - mU32))
         + (48*(2*(mQ32 + mU32)*pow2(mQ32)*pow2(mU32) - 4*m32*mQ32*mU32*pow2(
         mQ32 + mU32) - 2*(2*mQ32*mU32 + 3*pow2(mQ32) + 3*pow2(mU32))*pow3(m32)
         + 3*pow2(m32)*pow3(mQ32 + mU32) + 2*(mQ32 + mU32)*pow4(m32)))/(pow2(
         m32 - mQ32)*pow2(m32 - mU32)*pow2(mQ32 - mU32))) + ((8*(mQ32 + mU32)*(
         180*m32*mQ32*msq2 + 2*mQ32*pow2(m32) + 60*msq2*pow2(m32) - 3*m32*pow2(
         mQ32) - 90*m32*pow2(msq2) - 30*mQ32*pow2(msq2) + pow3(mQ32)))/(3.*pow3
         (m32 - mQ32)*pow3(mQ32 - mU32)) - (8*(mQ32 + mU32)*(180*m32*msq2*mU32
         + 60*msq2*pow2(m32) + 2*mU32*pow2(m32) - 90*m32*pow2(msq2) - 30*mU32*
         pow2(msq2) - 3*m32*pow2(mU32) + pow3(mU32)))/(3.*pow3(m32 - mU32)*pow3
         (mQ32 - mU32)) - (24*(mQ32 + mU32)*(4*(mQ32 + mU32)*pow2(mQ32)*pow2(
         mU32) - 8*m32*mQ32*mU32*pow2(mQ32 + mU32) - 2*(6*mQ32*mU32 + 5*pow2(
         mQ32) + 5*pow2(mU32))*pow3(m32) + pow2(m32)*(19*mU32*pow2(mQ32) + 19*
         mQ32*pow2(mU32) + 5*pow3(mQ32) + 5*pow3(mU32)) + 4*(mQ32 + mU32)*pow4(
         m32)))/(pow2(m32 - mQ32)*pow2(m32 - mU32)*pow4(mQ32 - mU32)))*pow4(Xt)
         + ((-96*m3*(mQ32 + mU32)*(-2*mQ32*mU32 + m32*(mQ32 + mU32)))/((m32 -
         mQ32)*(m32 - mU32)*pow4(mQ32 - mU32)) - (32*(mQ32 + mU32)*(120*m3*mQ32
         *msq2 - m3*pow2(mQ32) - 60*m3*pow2(msq2) + mQ32*pow3(m3)))/(3.*pow2(
         m32 - mQ32)*pow4(mQ32 - mU32)) - (32*(mQ32 + mU32)*(120*m3*msq2*mU32 -
         60*m3*pow2(msq2) - m3*pow2(mU32) + mU32*pow3(m3)))/(3.*pow2(m32 -
         mU32)*pow4(-mQ32 + mU32)))*pow5(Xt) - (320*Xt*(-12*m3*mQ32*msq2*mU32 +
         6*m3*mQ32*pow2(msq2) + 6*m3*mU32*pow2(msq2) + mQ32*mU32*pow3(m3) - 12
         *pow2(msq2)*pow3(m3) - mQ32*pow5(m3) + 12*msq2*pow5(m3) - mU32*pow5(m3
         ) + pow7(m3)))/(3.*pow2(m32 - mQ32)*pow2(m32 - mU32))) + (8*pow4(Xt)*(
         -60*msq2*pow2(mQ32)*pow2(mU32)*pow3(m32) + 360*msq2*pow2(m32)*pow2(
         mU32)*pow3(mQ32) + 21138*pow2(mU32)*pow3(m32)*pow3(mQ32) + 360*msq2*
         pow2(m32)*pow2(mQ32)*pow3(mU32) + 20522*pow2(mQ32)*pow3(m32)*pow3(mU32
         ) - 180*m32*msq2*pow3(mQ32)*pow3(mU32) - 12116*pow2(m32)*pow3(mQ32)*
         pow3(mU32) - 90*msq2*mU32*pow2(mQ32)*pow4(m32) - 90*mQ32*msq2*pow2(
         mU32)*pow4(m32) - 24928*pow2(mQ32)*pow2(mU32)*pow4(m32) + 90*msq2*pow3
         (mQ32)*pow4(m32) - 11670*mU32*pow3(mQ32)*pow4(m32) - 11842*mQ32*pow3(
         mU32)*pow4(m32) + 90*msq2*pow3(mU32)*pow4(m32) - 360*msq2*mU32*pow2(
         m32)*pow4(mQ32) - 6833*pow2(m32)*pow2(mU32)*pow4(mQ32) + 30*msq2*pow3(
         m32)*pow4(mQ32) + 3785*mU32*pow3(m32)*pow4(mQ32) + 2260*m32*pow3(mU32)
         *pow4(mQ32) - 30*msq2*pow3(mU32)*pow4(mQ32) - 8*pow4(m32)*pow4(mQ32) -
         360*mQ32*msq2*pow2(m32)*pow4(mU32) - 6555*pow2(m32)*pow2(mQ32)*pow4(
         mU32) + 4297*mQ32*pow3(m32)*pow4(mU32) + 30*msq2*pow3(m32)*pow4(mU32)
         + 2334*m32*pow3(mQ32)*pow4(mU32) - 30*msq2*pow3(mQ32)*pow4(mU32) - 32*
         pow4(m32)*pow4(mU32) + 256*pow4(mQ32)*pow4(mU32) + 10898*mU32*pow2(
         mQ32)*pow5(m32) + 10982*mQ32*pow2(mU32)*pow5(m32) + 1054*pow3(mQ32)*
         pow5(m32) + 1066*pow3(mU32)*pow5(m32) + 90*m32*msq2*mU32*pow5(mQ32) +
         252*mU32*pow2(m32)*pow5(mQ32) + 557*m32*pow2(mU32)*pow5(mQ32) + 30*
         msq2*pow2(mU32)*pow5(mQ32) - 207*pow3(m32)*pow5(mQ32) - 210*pow3(mU32)
         *pow5(mQ32) + 90*m32*mQ32*msq2*pow5(mU32) - 112*mQ32*pow2(m32)*pow5(
         mU32) + 475*m32*pow2(mQ32)*pow5(mU32) + 30*msq2*pow2(mQ32)*pow5(mU32)
         - 255*pow3(m32)*pow5(mU32) - 220*pow3(mQ32)*pow5(mU32) - 3128*mQ32*
         mU32*pow6(m32) - 822*pow2(mQ32)*pow6(m32) - 818*pow2(mU32)*pow6(m32) -
         91*m32*mU32*pow6(mQ32) - 67*pow2(m32)*pow6(mQ32) + 2*pow2(mU32)*pow6(
         mQ32) + 15*m32*mQ32*pow6(mU32) - 9*pow2(m32)*pow6(mU32) + 6*pow2(mQ32)
         *pow6(mU32) + 4*mQ32*pow7(m32) - 4*mU32*pow7(m32) + 18*m32*pow7(mQ32)
         + 6*mU32*pow7(mQ32)))/(3.*pow3(m32 - mQ32)*pow3(m32 - mU32)*pow4(-mQ32
         + mU32)) - (64*pow3(Xt)*(240*msq2*mU32*pow2(mQ32)*pow3(m3) - 240*mQ32
         *msq2*pow2(mU32)*pow3(m3) - 8*pow2(mQ32)*pow2(mU32)*pow3(m3) - 60*m3*
         msq2*mU32*pow3(mQ32) - 337*m3*pow2(mU32)*pow3(mQ32) + 612*mU32*pow3(m3
         )*pow3(mQ32) + 60*m3*mQ32*msq2*pow3(mU32) + 341*m3*pow2(mQ32)*pow3(
         mU32) - 620*mQ32*pow3(m3)*pow3(mU32) - 12*m3*mU32*pow4(mQ32) + 6*pow3(
         m3)*pow4(mQ32) + 12*m3*mQ32*pow4(mU32) - 6*pow3(m3)*pow4(mU32) - 60*
         msq2*pow2(mQ32)*pow5(m3) - 669*mU32*pow2(mQ32)*pow5(m3) + 689*mQ32*
         pow2(mU32)*pow5(m3) + 60*msq2*pow2(mU32)*pow5(m3) - 171*pow3(mQ32)*
         pow5(m3) + 175*pow3(mU32)*pow5(m3) - 8*mQ32*mU32*pow7(m3) + 124*pow2(
         mQ32)*pow7(m3) - 132*pow2(mU32)*pow7(m3) + 63*mQ32*pow9(m3) - 59*mU32*
         pow9(m3)))/(3.*pow2(m32 - mQ32)*pow2(m32 - mU32)*pow3(mQ32 - mU32)) +
         (32*Xt*(-120*m3*msq2*pow2(mQ32)*pow2(mU32) - 317*pow2(mQ32)*pow2(mU32)
         *pow3(m3) + 60*m3*msq2*mU32*pow3(mQ32) + 39*m3*pow2(mU32)*pow3(mQ32) +
         117*mU32*pow3(m3)*pow3(mQ32) + 60*m3*mQ32*msq2*pow3(mU32) + 7*m3*pow2
         (mQ32)*pow3(mU32) + 102*mQ32*pow3(m3)*pow3(mU32) - 20*m3*mU32*pow4(
         mQ32) - 4*pow3(m3)*pow4(mQ32) + 6*pow3(m3)*pow4(mU32) + 120*mQ32*msq2*
         mU32*pow5(m3) - 60*msq2*pow2(mQ32)*pow5(m3) + 250*mU32*pow2(mQ32)*pow5
         (m3) + 271*mQ32*pow2(mU32)*pow5(m3) - 60*msq2*pow2(mU32)*pow5(m3) -
         206*pow3(mQ32)*pow5(m3) - 219*pow3(mU32)*pow5(m3) + 6*m3*pow5(mQ32) -
         467*mQ32*mU32*pow7(m3) + 214*pow2(mQ32)*pow7(m3) + 221*pow2(mU32)*pow7
         (m3) + 2*mQ32*pow9(m3) - 2*mU32*pow9(m3)))/(3.*pow2(m32 - mQ32)*pow2(
         m32 - mU32)*pow2(mQ32 - mU32)) + (64*pow5(Xt)*(-60*m3*msq2*pow2(mQ32)*
         pow2(mU32) + 120*msq2*mU32*pow2(mQ32)*pow3(m3) + 120*mQ32*msq2*pow2(
         mU32)*pow3(m3) + 920*pow2(mQ32)*pow2(mU32)*pow3(m3) - 30*m3*msq2*mU32*
         pow3(mQ32) - 298*m3*pow2(mU32)*pow3(mQ32) + 352*mU32*pow3(m3)*pow3(
         mQ32) - 30*m3*mQ32*msq2*pow3(mU32) - 268*m3*pow2(mQ32)*pow3(mU32) +
         372*mQ32*pow3(m3)*pow3(mU32) + 14*m3*mU32*pow4(mQ32) + 13*pow3(m3)*
         pow4(mQ32) - 6*m3*mQ32*pow4(mU32) + 3*pow3(m3)*pow4(mU32) - 60*mQ32*
         msq2*mU32*pow5(m3) - 30*msq2*pow2(mQ32)*pow5(m3) - 638*mU32*pow2(mQ32)
         *pow5(m3) - 668*mQ32*pow2(mU32)*pow5(m3) - 30*msq2*pow2(mU32)*pow5(m3)
         - 102*pow3(mQ32)*pow5(m3) - 92*pow3(mU32)*pow5(m3) - 6*m3*pow5(mQ32)
         + 188*mQ32*mU32*pow7(m3) + 44*pow2(mQ32)*pow7(m3) + 44*pow2(mU32)*pow7
         (m3) + 64*mQ32*pow9(m3) + 64*mU32*pow9(m3)))/(3.*pow2(m32 - mQ32)*pow2
         (m32 - mU32)*pow4(mQ32 - mU32))) + pow2(log(mU32/MR2))*((-16*(-18*m32*
         mQ32 + pow2(m32) + 9*pow2(mQ32))*(-4*m32*mU32 + 3*pow2(m32) + 2*pow2(
         mU32)))/(3.*pow2(m32 - mQ32)*pow2(m32 - mU32)) + (8*(-34*m32*mU32 +
         pow2(m32) + 17*pow2(mU32))*(4*m32*mQ32*mU32*(mQ32 + mU32) - 2*pow2(
         mQ32)*pow2(mU32) - pow2(m32)*(8*mQ32*mU32 + pow2(mQ32) + pow2(mU32)) +
         2*(mQ32 + mU32)*pow3(m32)))/(3.*pow2(m32 - mQ32)*pow4(m32 - mU32)) +
         (8*(184*pow2(m32)*pow2(mQ32)*pow2(mU32) - 194*mU32*pow2(mQ32)*pow3(m32
         ) - 87*mQ32*pow2(mU32)*pow3(m32) + 122*mU32*pow2(m32)*pow3(mQ32) - 68*
         m32*pow2(mU32)*pow3(mQ32) + 34*pow3(m32)*pow3(mQ32) + 29*mQ32*pow2(m32
         )*pow3(mU32) - 60*m32*pow2(mQ32)*pow3(mU32) - pow3(m32)*pow3(mU32) +
         20*pow3(mQ32)*pow3(mU32) + 70*mQ32*mU32*pow4(m32) - 14*pow2(mQ32)*pow4
         (m32) + 3*pow2(mU32)*pow4(m32) - 35*m32*mU32*pow4(mQ32) - 29*pow2(m32)
         *pow4(mQ32) + 4*pow2(mU32)*pow4(mQ32) + 12*mQ32*pow5(m32) - 2*mU32*
         pow5(m32) + 9*m32*pow5(mQ32) + 3*mU32*pow5(mQ32)))/(3.*pow3(m32 - mQ32
         )*pow3(m32 - mU32)) - (4*(-120*mQ32*msq2*mU32*pow2(m32) - 33*mU32*pow2
         (m32)*pow2(mQ32) - 60*m32*mQ32*mU32*pow2(msq2) + 90*mQ32*pow2(m32)*
         pow2(msq2) - 90*mU32*pow2(m32)*pow2(msq2) + 180*m32*mQ32*msq2*pow2(
         mU32) - 51*mQ32*pow2(m32)*pow2(mU32) + 120*msq2*pow2(m32)*pow2(mU32) +
         36*m32*pow2(mQ32)*pow2(mU32) + 60*m32*pow2(msq2)*pow2(mU32) - 30*mQ32
         *pow2(msq2)*pow2(mU32) - 60*mQ32*msq2*pow3(m32) + 64*mQ32*mU32*pow3(
         m32) + 60*msq2*mU32*pow3(m32) - 2*pow2(mQ32)*pow3(m32) + 706*pow2(mU32
         )*pow3(m32) - 6*m32*mU32*pow3(mQ32) + 9*pow2(m32)*pow3(mQ32) - 3*pow2(
         mU32)*pow3(mQ32) - 38*m32*mQ32*pow3(mU32) - 180*m32*msq2*pow3(mU32) -
         437*pow2(m32)*pow3(mU32) - pow2(mQ32)*pow3(mU32) + 30*pow2(msq2)*pow3(
         mU32) + 44*mQ32*pow4(m32) - 556*mU32*pow4(m32) + 136*m32*pow4(mU32) +
         5*mQ32*pow4(mU32) + 128*pow5(m32) - pow5(mU32)))/(3.*(-mQ32 + mU32)*
         pow4(m32 - mU32)) + (1280*m32*(mQ32 + mU32)*(2*m32*mQ32*mU32 + m32*
         pow2(mU32) - 3*mQ32*pow2(mU32))*pow6(Xt))/(3.*(m32 - mQ32)*pow2(m32 -
         mU32)*pow5(mQ32 - mU32)) + (8*pow2(Xt)*(270*pow2(m32)*pow2(mQ32)*pow2(
         msq2)*pow2(mU32) - 90*mU32*pow2(mQ32)*pow2(msq2)*pow3(m32) - 180*msq2*
         pow2(mQ32)*pow2(mU32)*pow3(m32) - 210*mQ32*pow2(msq2)*pow2(mU32)*pow3(
         m32) - 90*mU32*pow2(m32)*pow2(msq2)*pow3(mQ32) + 420*msq2*pow2(m32)*
         pow2(mU32)*pow3(mQ32) - 150*m32*pow2(msq2)*pow2(mU32)*pow3(mQ32) - 420
         *msq2*mU32*pow3(m32)*pow3(mQ32) + 270*pow2(msq2)*pow3(m32)*pow3(mQ32)
         - 1992*pow2(mU32)*pow3(m32)*pow3(mQ32) - 540*msq2*pow2(m32)*pow2(mQ32)
         *pow3(mU32) - 90*mQ32*pow2(m32)*pow2(msq2)*pow3(mU32) + 90*m32*pow2(
         mQ32)*pow2(msq2)*pow3(mU32) + 540*mQ32*msq2*pow3(m32)*pow3(mU32) -
         2924*pow2(mQ32)*pow3(m32)*pow3(mU32) + 30*pow2(msq2)*pow3(m32)*pow3(
         mU32) + 180*m32*msq2*pow3(mQ32)*pow3(mU32) + 424*pow2(m32)*pow3(mQ32)*
         pow3(mU32) - 30*pow2(msq2)*pow3(mQ32)*pow3(mU32) + 540*msq2*mU32*pow2(
         mQ32)*pow4(m32) + 210*mQ32*mU32*pow2(msq2)*pow4(m32) - 270*pow2(mQ32)*
         pow2(msq2)*pow4(m32) - 180*mQ32*msq2*pow2(mU32)*pow4(m32) + 6259*pow2(
         mQ32)*pow2(mU32)*pow4(m32) + 60*pow2(msq2)*pow2(mU32)*pow4(m32) - 180*
         msq2*pow3(mQ32)*pow4(m32) + 1947*mU32*pow3(mQ32)*pow4(m32) + 1957*mQ32
         *pow3(mU32)*pow4(m32) - 180*msq2*pow3(mU32)*pow4(m32) + 120*msq2*mU32*
         pow2(m32)*pow4(mQ32) + 60*m32*mU32*pow2(msq2)*pow4(mQ32) - 90*pow2(m32
         )*pow2(msq2)*pow4(mQ32) - 180*m32*msq2*pow2(mU32)*pow4(mQ32) - 238*
         pow2(m32)*pow2(mU32)*pow4(mQ32) + 30*pow2(msq2)*pow2(mU32)*pow4(mQ32)
         + 60*msq2*pow3(m32)*pow4(mQ32) + 49*mU32*pow3(m32)*pow4(mQ32) + 497*
         m32*pow3(mU32)*pow4(mQ32) - 187*pow4(m32)*pow4(mQ32) + 713*pow2(m32)*
         pow2(mQ32)*pow4(mU32) - 235*mQ32*pow3(m32)*pow4(mU32) - 109*m32*pow3(
         mQ32)*pow4(mU32) - 376*pow4(m32)*pow4(mU32) - 113*pow4(mQ32)*pow4(mU32
         ) - 300*mQ32*msq2*mU32*pow5(m32) + 180*msq2*pow2(mQ32)*pow5(m32) -
         5047*mU32*pow2(mQ32)*pow5(m32) + 90*mQ32*pow2(msq2)*pow5(m32) - 90*
         mU32*pow2(msq2)*pow5(m32) - 4593*mQ32*pow2(mU32)*pow5(m32) + 120*msq2*
         pow2(mU32)*pow5(m32) - 339*pow3(mQ32)*pow5(m32) + 283*pow3(mU32)*pow5(
         m32) + 45*mU32*pow2(m32)*pow5(mQ32) - 135*m32*pow2(mU32)*pow5(mQ32) +
         87*pow3(m32)*pow5(mQ32) + 3*pow3(mU32)*pow5(mQ32) + 11*mQ32*pow2(m32)*
         pow5(mU32) - 111*m32*pow2(mQ32)*pow5(mU32) + 87*pow3(m32)*pow5(mU32) +
         37*pow3(mQ32)*pow5(mU32) - 60*mQ32*msq2*pow6(m32) + 3672*mQ32*mU32*
         pow6(m32) + 60*msq2*mU32*pow6(m32) + 1158*pow2(mQ32)*pow6(m32) + 450*
         pow2(mU32)*pow6(m32) + 18*m32*mU32*pow6(mQ32) - 27*pow2(m32)*pow6(mQ32
         ) + 9*pow2(mU32)*pow6(mQ32) - 844*mQ32*pow7(m32) - 564*mU32*pow7(m32)
         + 128*pow8(m32)))/(3.*pow2(mQ32 - mU32)*pow3(m32 - mQ32)*pow4(m32 -
         mU32)) - (4*pow4(Xt)*(-300*pow2(mQ32)*pow2(msq2)*pow2(mU32)*pow3(m32)
         + 180*pow2(m32)*pow2(msq2)*pow2(mU32)*pow3(mQ32) + 180*mU32*pow2(msq2)
         *pow3(m32)*pow3(mQ32) - 600*msq2*pow2(mU32)*pow3(m32)*pow3(mQ32) + 180
         *pow2(m32)*pow2(mQ32)*pow2(msq2)*pow3(mU32) + 360*msq2*pow2(mQ32)*pow3
         (m32)*pow3(mU32) - 180*mQ32*pow2(msq2)*pow3(m32)*pow3(mU32) - 120*msq2
         *pow2(m32)*pow3(mQ32)*pow3(mU32) - 60*m32*pow2(msq2)*pow3(mQ32)*pow3(
         mU32) - 18052*pow3(m32)*pow3(mQ32)*pow3(mU32) - 60*mU32*pow2(mQ32)*
         pow2(msq2)*pow4(m32) + 360*msq2*pow2(mQ32)*pow2(mU32)*pow4(m32) + 270*
         mQ32*pow2(msq2)*pow2(mU32)*pow4(m32) + 360*msq2*mU32*pow3(mQ32)*pow4(
         m32) - 270*pow2(msq2)*pow3(mQ32)*pow4(m32) + 24958*pow2(mU32)*pow3(
         mQ32)*pow4(m32) - 360*mQ32*msq2*pow3(mU32)*pow4(m32) + 15560*pow2(mQ32
         )*pow3(mU32)*pow4(m32) + 60*pow2(msq2)*pow3(mU32)*pow4(m32) - 180*mU32
         *pow2(m32)*pow2(msq2)*pow4(mQ32) + 540*msq2*pow2(m32)*pow2(mU32)*pow4(
         mQ32) - 90*m32*pow2(msq2)*pow2(mU32)*pow4(mQ32) - 360*msq2*mU32*pow3(
         m32)*pow4(mQ32) + 270*pow2(msq2)*pow3(m32)*pow4(mQ32) - 7207*pow2(mU32
         )*pow3(m32)*pow4(mQ32) + 5850*pow2(m32)*pow3(mU32)*pow4(mQ32) - 180*
         msq2*pow4(m32)*pow4(mQ32) + 3584*mU32*pow4(m32)*pow4(mQ32) - 540*msq2*
         pow2(m32)*pow2(mQ32)*pow4(mU32) - 90*mQ32*pow2(m32)*pow2(msq2)*pow4(
         mU32) + 90*m32*pow2(mQ32)*pow2(msq2)*pow4(mU32) + 540*mQ32*msq2*pow3(
         m32)*pow4(mU32) + 2513*pow2(mQ32)*pow3(m32)*pow4(mU32) + 30*pow2(msq2)
         *pow3(m32)*pow4(mU32) + 180*m32*msq2*pow3(mQ32)*pow4(mU32) + 2185*pow2
         (m32)*pow3(mQ32)*pow4(mU32) - 30*pow2(msq2)*pow3(mQ32)*pow4(mU32) -
         4427*mQ32*pow4(m32)*pow4(mU32) - 180*msq2*pow4(m32)*pow4(mU32) - 1604*
         m32*pow4(mQ32)*pow4(mU32) - 120*msq2*mU32*pow2(mQ32)*pow5(m32) + 90*
         pow2(mQ32)*pow2(msq2)*pow5(m32) - 180*mQ32*msq2*pow2(mU32)*pow5(m32) -
         26136*pow2(mQ32)*pow2(mU32)*pow5(m32) - 90*pow2(msq2)*pow2(mU32)*pow5
         (m32) + 180*msq2*pow3(mQ32)*pow5(m32) - 13282*mU32*pow3(mQ32)*pow5(m32
         ) - 3654*mQ32*pow3(mU32)*pow5(m32) + 120*msq2*pow3(mU32)*pow5(m32) -
         507*pow4(mQ32)*pow5(m32) + 1595*pow4(mU32)*pow5(m32) + 120*msq2*mU32*
         pow2(m32)*pow5(mQ32) + 60*m32*mU32*pow2(msq2)*pow5(mQ32) - 90*pow2(m32
         )*pow2(msq2)*pow5(mQ32) - 180*m32*msq2*pow2(mU32)*pow5(mQ32) - 321*
         pow2(m32)*pow2(mU32)*pow5(mQ32) + 30*pow2(msq2)*pow2(mU32)*pow5(mQ32)
         + 60*msq2*pow3(m32)*pow5(mQ32) + 136*mU32*pow3(m32)*pow5(mQ32) + 490*
         m32*pow3(mU32)*pow5(mQ32) - 155*pow4(m32)*pow5(mQ32) - 142*pow4(mU32)*
         pow5(mQ32) - 3748*pow2(m32)*pow2(mQ32)*pow5(mU32) + 3172*mQ32*pow3(m32
         )*pow5(mU32) + 1740*m32*pow3(mQ32)*pow5(mU32) - 1120*pow4(m32)*pow5(
         mU32) - 140*pow4(mQ32)*pow5(mU32) - 60*msq2*pow2(mQ32)*pow6(m32) +
         14886*mU32*pow2(mQ32)*pow6(m32) + 9226*mQ32*pow2(mU32)*pow6(m32) + 60*
         msq2*pow2(mU32)*pow6(m32) + 2462*pow3(mQ32)*pow6(m32) - 206*pow3(mU32)
         *pow6(m32) + 18*mU32*pow2(m32)*pow6(mQ32) - 117*m32*pow2(mU32)*pow6(
         mQ32) + 87*pow3(m32)*pow6(mQ32) + 12*pow3(mU32)*pow6(mQ32) - 117*mQ32*
         pow2(m32)*pow6(mU32) - 15*m32*pow2(mQ32)*pow6(mU32) + 151*pow3(m32)*
         pow6(mU32) + 5*pow3(mQ32)*pow6(mU32) - 5384*mQ32*mU32*pow7(m32) - 3044
         *pow2(mQ32)*pow7(m32) - 532*pow2(mU32)*pow7(m32) + 18*m32*mU32*pow7(
         mQ32) - 27*pow2(m32)*pow7(mQ32) + 9*pow2(mU32)*pow7(mQ32) + 1160*mQ32*
         pow8(m32) + 120*mU32*pow8(m32)))/(3.*pow3(m32 - mQ32)*pow4(m32 - mU32)
         *pow4(mQ32 - mU32)) - (32*Xt*(22*mQ32*pow11(m3) - 22*mU32*pow11(m3) -
         30*m3*pow2(mQ32)*pow2(msq2)*pow2(mU32) - 30*mU32*pow2(mQ32)*pow2(msq2)
         *pow3(m3) + 60*msq2*pow2(mQ32)*pow2(mU32)*pow3(m3) + 60*mQ32*pow2(msq2
         )*pow2(mU32)*pow3(m3) + 30*m3*mU32*pow2(msq2)*pow3(mQ32) - 60*m3*msq2*
         pow2(mU32)*pow3(mQ32) + 60*msq2*mU32*pow3(m3)*pow3(mQ32) - 30*pow2(
         msq2)*pow3(m3)*pow3(mQ32) + 106*pow2(mU32)*pow3(m3)*pow3(mQ32) + 60*m3
         *msq2*pow2(mQ32)*pow3(mU32) - 120*mQ32*msq2*pow3(m3)*pow3(mU32) + 25*
         pow2(mQ32)*pow3(m3)*pow3(mU32) - 38*m3*pow3(mQ32)*pow3(mU32) + 10*m3*
         pow2(mU32)*pow4(mQ32) - 5*mU32*pow3(m3)*pow4(mQ32) + 15*m3*pow2(mQ32)*
         pow4(mU32) - 65*mQ32*pow3(m3)*pow4(mU32) - 120*msq2*mU32*pow2(mQ32)*
         pow5(m3) - 30*mQ32*mU32*pow2(msq2)*pow5(m3) + 60*pow2(mQ32)*pow2(msq2)
         *pow5(m3) + 60*mQ32*msq2*pow2(mU32)*pow5(m3) - 263*pow2(mQ32)*pow2(
         mU32)*pow5(m3) - 30*pow2(msq2)*pow2(mU32)*pow5(m3) - 66*mU32*pow3(mQ32
         )*pow5(m3) + 198*mQ32*pow3(mU32)*pow5(m3) + 60*msq2*pow3(mU32)*pow5(m3
         ) - 5*pow4(mQ32)*pow5(m3) + 40*pow4(mU32)*pow5(m3) - 3*m3*mU32*pow5(
         mQ32) + 3*pow3(m3)*pow5(mQ32) + 60*mQ32*msq2*mU32*pow7(m3) + 235*mU32*
         pow2(mQ32)*pow7(m3) - 30*mQ32*pow2(msq2)*pow7(m3) + 30*mU32*pow2(msq2)
         *pow7(m3) - 56*mQ32*pow2(mU32)*pow7(m3) - 60*msq2*pow2(mU32)*pow7(m3)
         + 14*pow3(mQ32)*pow7(m3) - 129*pow3(mU32)*pow7(m3) - 51*mQ32*mU32*pow9
         (m3) - 60*pow2(mQ32)*pow9(m3) + 95*pow2(mU32)*pow9(m3)))/(3.*pow2(m32
         - mQ32)*pow2(mQ32 - mU32)*pow3(m32 - mU32)) - (64*pow3(Xt)*(-19*mQ32*
         pow11(m3) - 25*mU32*pow11(m3) + 2*pow13(m3) - 30*m3*pow2(mQ32)*pow2(
         msq2)*pow2(mU32) - 30*mU32*pow2(mQ32)*pow2(msq2)*pow3(m3) + 60*msq2*
         pow2(mQ32)*pow2(mU32)*pow3(m3) + 60*mQ32*pow2(msq2)*pow2(mU32)*pow3(m3
         ) + 30*m3*mU32*pow2(msq2)*pow3(mQ32) - 60*m3*msq2*pow2(mU32)*pow3(mQ32
         ) + 60*msq2*mU32*pow3(m3)*pow3(mQ32) - 30*pow2(msq2)*pow3(m3)*pow3(
         mQ32) + 45*pow2(mU32)*pow3(m3)*pow3(mQ32) + 60*m3*msq2*pow2(mQ32)*pow3
         (mU32) - 120*mQ32*msq2*pow3(m3)*pow3(mU32) - 530*pow2(mQ32)*pow3(m3)*
         pow3(mU32) + 47*m3*pow3(mQ32)*pow3(mU32) - 30*m3*pow2(mU32)*pow4(mQ32)
         + 15*mU32*pow3(m3)*pow4(mQ32) + 136*m3*pow2(mQ32)*pow4(mU32) - 205*
         mQ32*pow3(m3)*pow4(mU32) - 120*msq2*mU32*pow2(mQ32)*pow5(m3) - 30*mQ32
         *mU32*pow2(msq2)*pow5(m3) + 60*pow2(mQ32)*pow2(msq2)*pow5(m3) + 60*
         mQ32*msq2*pow2(mU32)*pow5(m3) + 442*pow2(mQ32)*pow2(mU32)*pow5(m3) -
         30*pow2(msq2)*pow2(mU32)*pow5(m3) - 91*mU32*pow3(mQ32)*pow5(m3) + 681*
         mQ32*pow3(mU32)*pow5(m3) + 60*msq2*pow3(mU32)*pow5(m3) + 15*pow4(mQ32)
         *pow5(m3) + 71*pow4(mU32)*pow5(m3) + 9*m3*mU32*pow5(mQ32) - 9*pow3(m3)
         *pow5(mQ32) + 60*mQ32*msq2*mU32*pow7(m3) - 78*mU32*pow2(mQ32)*pow7(m3)
         - 30*mQ32*pow2(msq2)*pow7(m3) + 30*mU32*pow2(msq2)*pow7(m3) - 619*
         mQ32*pow2(mU32)*pow7(m3) - 60*msq2*pow2(mU32)*pow7(m3) + 15*pow3(mQ32)
         *pow7(m3) - 190*pow3(mU32)*pow7(m3) + 210*mQ32*mU32*pow9(m3) - 18*pow2
         (mQ32)*pow9(m3) + 126*pow2(mU32)*pow9(m3)))/(3.*pow2(m32 - mQ32)*pow3(
         m32 - mU32)*pow3(mQ32 - mU32)) - (32*pow5(Xt)*(-32*mQ32*mU32*pow11(m3)
         + 14*pow11(m3)*pow2(mQ32) + 18*pow11(m3)*pow2(mU32) - 30*pow2(mQ32)*
         pow2(msq2)*pow2(mU32)*pow3(m3) + 60*mU32*pow2(msq2)*pow3(m3)*pow3(mQ32
         ) - 120*msq2*pow2(mU32)*pow3(m3)*pow3(mQ32) + 30*m3*pow2(mQ32)*pow2(
         msq2)*pow3(mU32) + 60*msq2*pow2(mQ32)*pow3(m3)*pow3(mU32) - 60*mQ32*
         pow2(msq2)*pow3(m3)*pow3(mU32) + 731*pow3(m3)*pow3(mQ32)*pow3(mU32) -
         30*m3*mU32*pow2(msq2)*pow4(mQ32) + 60*m3*msq2*pow2(mU32)*pow4(mQ32) -
         60*msq2*mU32*pow3(m3)*pow4(mQ32) + 30*pow2(msq2)*pow3(m3)*pow4(mQ32) +
         25*pow2(mU32)*pow3(m3)*pow4(mQ32) - 56*m3*pow3(mU32)*pow4(mQ32) - 60*
         m3*msq2*pow2(mQ32)*pow4(mU32) + 120*mQ32*msq2*pow3(m3)*pow4(mU32) +
         894*pow2(mQ32)*pow3(m3)*pow4(mU32) - 261*m3*pow3(mQ32)*pow4(mU32) - 30
         *mU32*pow2(mQ32)*pow2(msq2)*pow5(m3) + 60*msq2*pow2(mQ32)*pow2(mU32)*
         pow5(m3) + 60*mQ32*pow2(msq2)*pow2(mU32)*pow5(m3) + 120*msq2*mU32*pow3
         (mQ32)*pow5(m3) - 60*pow2(msq2)*pow3(mQ32)*pow5(m3) - 621*pow2(mU32)*
         pow3(mQ32)*pow5(m3) - 120*mQ32*msq2*pow3(mU32)*pow5(m3) - 1297*pow2(
         mQ32)*pow3(mU32)*pow5(m3) + 30*pow2(msq2)*pow3(mU32)*pow5(m3) + 23*
         mU32*pow4(mQ32)*pow5(m3) - 868*mQ32*pow4(mU32)*pow5(m3) - 60*msq2*pow4
         (mU32)*pow5(m3) + 21*m3*pow2(mU32)*pow5(mQ32) - 6*mU32*pow3(m3)*pow5(
         mQ32) - 15*pow5(m3)*pow5(mQ32) - 175*m3*pow2(mQ32)*pow5(mU32) + 267*
         mQ32*pow3(m3)*pow5(mU32) - 102*pow5(m3)*pow5(mU32) - 9*m3*mU32*pow6(
         mQ32) + 9*pow3(m3)*pow6(mQ32) - 60*msq2*mU32*pow2(mQ32)*pow7(m3) + 30*
         pow2(mQ32)*pow2(msq2)*pow7(m3) + 719*pow2(mQ32)*pow2(mU32)*pow7(m3) -
         30*pow2(msq2)*pow2(mU32)*pow7(m3) + 177*mU32*pow3(mQ32)*pow7(m3) + 783
         *mQ32*pow3(mU32)*pow7(m3) + 60*msq2*pow3(mU32)*pow7(m3) - 8*pow4(mQ32)
         *pow7(m3) + 249*pow4(mU32)*pow7(m3) - 155*mU32*pow2(mQ32)*pow9(m3) -
         182*mQ32*pow2(mU32)*pow9(m3) + 6*pow3(mQ32)*pow9(m3) - 149*pow3(mU32)*
         pow9(m3)))/(3.*pow2(m32 - mQ32)*pow3(m32 - mU32)*pow5(mQ32 - mU32))))
         + pow2(log(mQ32/MR2))*((-16*pow2(Xt)*(210*msq2*mU32*pow2(m32)*pow2(
         mQ32) - 180*mQ32*msq2*pow2(m32)*pow2(mU32) - 150*m32*msq2*pow2(mQ32)*
         pow2(mU32) - 1624*pow2(m32)*pow2(mQ32)*pow2(mU32) + 90*mQ32*msq2*mU32*
         pow3(m32) - 90*msq2*pow2(mQ32)*pow3(m32) + 1455*mU32*pow2(mQ32)*pow3(
         m32) + 1372*mQ32*pow2(mU32)*pow3(m32) + 60*m32*msq2*mU32*pow3(mQ32) -
         30*msq2*pow2(m32)*pow3(mQ32) + 161*mU32*pow2(m32)*pow3(mQ32) + 332*m32
         *pow2(mU32)*pow3(mQ32) - 30*msq2*pow2(mU32)*pow3(mQ32) - 240*pow3(m32)
         *pow3(mQ32) + 90*m32*mQ32*msq2*pow3(mU32) - 389*mQ32*pow2(m32)*pow3(
         mU32) + 485*m32*pow2(mQ32)*pow3(mU32) + 30*msq2*pow2(mQ32)*pow3(mU32)
         - 35*pow3(m32)*pow3(mU32) - 181*pow3(mQ32)*pow3(mU32) - 1466*mQ32*mU32
         *pow4(m32) - 367*pow2(mQ32)*pow4(m32) + 53*pow2(mU32)*pow4(m32) - 399*
         m32*mU32*pow4(mQ32) + 196*pow2(m32)*pow4(mQ32) + 155*pow2(mU32)*pow4(
         mQ32) + 9*m32*mQ32*pow4(mU32) + 3*pow2(mQ32)*pow4(mU32) + 486*mQ32*
         pow5(m32) - 18*mU32*pow5(m32) + 9*m32*pow5(mQ32) + 3*mU32*pow5(mQ32)))
         /(3.*pow2(m32 - mU32)*pow2(mQ32 - mU32)*pow3(m32 - mQ32)) + log(
         msq2/MR2)*((-8*(-180*m32*mQ32*msq2 - 65*mQ32*pow2(m32) - 60*msq2*pow2(
         m32) + 57*m32*pow2(mQ32) + 90*m32*pow2(msq2) + 30*mQ32*pow2(msq2) + 27
         *pow3(m32) - 19*pow3(mQ32)))/(3.*pow3(m32 - mQ32)) + pow2(Xt)*((16*(
         180*m32*mQ32*msq2 + 2*mQ32*pow2(m32) + 60*msq2*pow2(m32) - 3*m32*pow2(
         mQ32) - 90*m32*pow2(msq2) - 30*mQ32*pow2(msq2) + pow3(mQ32)))/(3.*(
         mQ32 - mU32)*pow3(m32 - mQ32)) - (48*((3*mQ32 - mU32)*pow2(m32) + m32*
         (4*mQ32*mU32 - 8*pow2(mQ32)) - 2*mU32*pow2(mQ32) + 4*pow3(mQ32)))/(
         pow2(m32 - mQ32)*pow2(mQ32 - mU32))) + ((-64*(120*m3*mQ32*msq2 - m3*
         pow2(mQ32) - 60*m3*pow2(msq2) + mQ32*pow3(m3)))/(3.*pow2(m32 - mQ32)*
         pow2(mQ32 - mU32)) + (96*m3*(3*mQ32*mU32 - m32*(3*mQ32 + mU32) + 2*
         pow2(m32) - pow2(mQ32)))/((m32 - mQ32)*pow3(mQ32 - mU32)))*pow3(Xt) +
         ((-8*(mQ32 + mU32)*(180*m32*mQ32*msq2 + 2*mQ32*pow2(m32) + 60*msq2*
         pow2(m32) - 3*m32*pow2(mQ32) - 90*m32*pow2(msq2) - 30*mQ32*pow2(msq2)
         + pow3(mQ32)))/(3.*pow3(m32 - mQ32)*pow3(mQ32 - mU32)) + (24*(8*mQ32*
         mU32*pow2(m32) - pow2(mQ32)*pow2(mU32) + 2*(mQ32 - mU32)*pow3(m32) + 4
         *mU32*pow3(mQ32) - 2*m32*(5*mU32*pow2(mQ32) - mQ32*pow2(mU32) + 4*pow3
         (mQ32)) + 5*pow4(mQ32)))/(pow2(m32 - mQ32)*pow4(mQ32 - mU32)))*pow4(Xt
         ) + (64*Xt*(-60*m3*mQ32*msq2 + 5*m3*pow2(mQ32) + 30*m3*pow2(msq2) - 14
         *mQ32*pow3(m3) + 9*pow5(m3)))/(3.*(mQ32 - mU32)*pow2(m32 - mQ32)) + (
         320*(mQ32 + mU32)*(12*m3*mQ32*msq2 - m3*pow2(mQ32) - 6*m3*pow2(msq2) +
         mQ32*pow3(m3))*pow5(Xt))/(3.*pow2(m32 - mQ32)*pow4(mQ32 - mU32))) - (
         2560*m32*pow2(mQ32)*pow6(Xt))/(3.*pow2(m32 - mQ32)*pow4(mQ32 - mU32))
         + (4*(270*pow2(m32)*pow2(mQ32)*pow2(msq2)*pow2(mU32) - 210*mU32*pow2(
         mQ32)*pow2(msq2)*pow3(m32) - 90*mQ32*pow2(msq2)*pow2(mU32)*pow3(m32) -
         90*mU32*pow2(m32)*pow2(msq2)*pow3(mQ32) - 1080*msq2*pow2(m32)*pow2(
         mU32)*pow3(mQ32) + 90*m32*pow2(msq2)*pow2(mU32)*pow3(mQ32) + 960*msq2*
         mU32*pow3(m32)*pow3(mQ32) + 30*pow2(msq2)*pow3(m32)*pow3(mQ32) - 5278*
         pow2(mU32)*pow3(m32)*pow3(mQ32) + 600*msq2*pow2(m32)*pow2(mQ32)*pow3(
         mU32) - 90*mQ32*pow2(m32)*pow2(msq2)*pow3(mU32) - 150*m32*pow2(mQ32)*
         pow2(msq2)*pow3(mU32) - 960*mQ32*msq2*pow3(m32)*pow3(mU32) - 6262*pow2
         (mQ32)*pow3(m32)*pow3(mU32) + 270*pow2(msq2)*pow3(m32)*pow3(mU32) +
         480*m32*msq2*pow3(mQ32)*pow3(mU32) + 4086*pow2(m32)*pow3(mQ32)*pow3(
         mU32) - 30*pow2(msq2)*pow3(mQ32)*pow3(mU32) - 600*msq2*mU32*pow2(mQ32)
         *pow4(m32) + 210*mQ32*mU32*pow2(msq2)*pow4(m32) + 60*pow2(mQ32)*pow2(
         msq2)*pow4(m32) + 1080*mQ32*msq2*pow2(mU32)*pow4(m32) + 9653*pow2(mQ32
         )*pow2(mU32)*pow4(m32) - 270*pow2(msq2)*pow2(mU32)*pow4(m32) - 300*
         msq2*pow3(mQ32)*pow4(m32) + 2119*mU32*pow3(mQ32)*pow4(m32) + 3565*mQ32
         *pow3(mU32)*pow4(m32) - 180*msq2*pow3(mU32)*pow4(m32) + 180*msq2*mU32*
         pow2(m32)*pow4(mQ32) - 180*m32*msq2*pow2(mU32)*pow4(mQ32) - 446*pow2(
         m32)*pow2(mU32)*pow4(mQ32) - 60*msq2*pow3(m32)*pow4(mQ32) + 1686*mU32*
         pow3(m32)*pow4(mQ32) - 599*m32*pow3(mU32)*pow4(mQ32) + 60*msq2*pow3(
         mU32)*pow4(mQ32) - 528*pow4(m32)*pow4(mQ32) + 300*mQ32*msq2*pow2(m32)*
         pow4(mU32) - 300*m32*msq2*pow2(mQ32)*pow4(mU32) + 1426*pow2(m32)*pow2(
         mQ32)*pow4(mU32) + 60*m32*mQ32*pow2(msq2)*pow4(mU32) - 90*pow2(m32)*
         pow2(msq2)*pow4(mU32) + 30*pow2(mQ32)*pow2(msq2)*pow4(mU32) - 746*mQ32
         *pow3(m32)*pow4(mU32) + 60*msq2*pow3(m32)*pow4(mU32) - 1034*m32*pow3(
         mQ32)*pow4(mU32) - 60*msq2*pow3(mQ32)*pow4(mU32) - 89*pow4(m32)*pow4(
         mU32) + 307*pow4(mQ32)*pow4(mU32) - 480*mQ32*msq2*mU32*pow5(m32) + 300
         *msq2*pow2(mQ32)*pow5(m32) - 5861*mU32*pow2(mQ32)*pow5(m32) - 90*mQ32*
         pow2(msq2)*pow5(m32) + 90*mU32*pow2(msq2)*pow5(m32) - 6037*mQ32*pow2(
         mU32)*pow5(m32) + 180*msq2*pow2(mU32)*pow5(m32) - 307*pow3(mQ32)*pow5(
         m32) + 45*pow3(mU32)*pow5(m32) - 1172*mU32*pow2(m32)*pow5(mQ32) + 979*
         m32*pow2(mU32)*pow5(mQ32) + 232*pow3(m32)*pow5(mQ32) - 291*pow3(mU32)*
         pow5(mQ32) + 18*mQ32*pow2(m32)*pow5(mU32) - 12*m32*pow2(mQ32)*pow5(
         mU32) - 6*pow3(mQ32)*pow5(mU32) + 60*mQ32*msq2*pow6(m32) + 4132*mQ32*
         mU32*pow6(m32) - 60*msq2*mU32*pow6(m32) + 1308*pow2(mQ32)*pow6(m32) +
         320*pow2(mU32)*pow6(m32) + 35*m32*mU32*pow6(mQ32) + 56*pow2(m32)*pow6(
         mQ32) - 7*pow2(mU32)*pow6(mQ32) - 1016*mQ32*pow7(m32) - 392*mU32*pow7(
         m32) - 9*m32*pow7(mQ32) - 3*mU32*pow7(mQ32) + 128*pow8(m32)))/(3.*(
         mQ32 - mU32)*pow3(m32 - mU32)*pow4(m32 - mQ32)) - (8*pow4(Xt)*(120*
         pow2(mQ32)*pow2(msq2)*pow2(mU32)*pow3(m32) + 360*pow2(m32)*pow2(msq2)*
         pow2(mU32)*pow3(mQ32) - 240*mU32*pow2(msq2)*pow3(m32)*pow3(mQ32) - 420
         *msq2*pow2(mU32)*pow3(m32)*pow3(mQ32) - 360*pow2(m32)*pow2(mQ32)*pow2(
         msq2)*pow3(mU32) - 420*msq2*pow2(mQ32)*pow3(m32)*pow3(mU32) + 360*mQ32
         *pow2(msq2)*pow3(m32)*pow3(mU32) + 780*msq2*pow2(m32)*pow3(mQ32)*pow3(
         mU32) - 240*m32*pow2(msq2)*pow3(mQ32)*pow3(mU32) - 15392*pow3(m32)*
         pow3(mQ32)*pow3(mU32) + 150*mU32*pow2(mQ32)*pow2(msq2)*pow4(m32) + 780
         *msq2*pow2(mQ32)*pow2(mU32)*pow4(m32) - 480*mQ32*pow2(msq2)*pow2(mU32)
         *pow4(m32) - 270*msq2*mU32*pow3(mQ32)*pow4(m32) + 60*pow2(msq2)*pow3(
         mQ32)*pow4(m32) + 25400*pow2(mU32)*pow3(mQ32)*pow4(m32) - 450*mQ32*
         msq2*pow3(mU32)*pow4(m32) + 9104*pow2(mQ32)*pow3(mU32)*pow4(m32) + 270
         *pow2(msq2)*pow3(mU32)*pow4(m32) - 90*mU32*pow2(m32)*pow2(msq2)*pow4(
         mQ32) - 720*msq2*pow2(m32)*pow2(mU32)*pow4(mQ32) + 90*m32*pow2(msq2)*
         pow2(mU32)*pow4(mQ32) + 720*msq2*mU32*pow3(m32)*pow4(mQ32) + 30*pow2(
         msq2)*pow3(m32)*pow4(mQ32) - 17072*pow2(mU32)*pow3(m32)*pow4(mQ32) +
         240*m32*msq2*pow3(mU32)*pow4(mQ32) + 9007*pow2(m32)*pow3(mU32)*pow4(
         mQ32) - 30*pow2(msq2)*pow3(mU32)*pow4(mQ32) - 240*msq2*pow4(m32)*pow4(
         mQ32) + 13130*mU32*pow4(m32)*pow4(mQ32) - 120*msq2*pow2(m32)*pow2(mQ32
         )*pow4(mU32) + 210*m32*pow2(mQ32)*pow2(msq2)*pow4(mU32) + 210*mQ32*
         msq2*pow3(m32)*pow4(mU32) - 3971*pow2(mQ32)*pow3(m32)*pow4(mU32) - 270
         *pow2(msq2)*pow3(m32)*pow4(mU32) - 270*m32*msq2*pow3(mQ32)*pow4(mU32)
         + 4850*pow2(m32)*pow3(mQ32)*pow4(mU32) + 60*pow2(msq2)*pow3(mQ32)*pow4
         (mU32) + 880*mQ32*pow4(m32)*pow4(mU32) + 180*msq2*pow4(m32)*pow4(mU32)
         - 2108*m32*pow4(mQ32)*pow4(mU32) - 420*msq2*mU32*pow2(mQ32)*pow5(m32)
         + 180*mQ32*mU32*pow2(msq2)*pow5(m32) - 90*pow2(mQ32)*pow2(msq2)*pow5(
         m32) + 390*mQ32*msq2*pow2(mU32)*pow5(m32) - 12473*pow2(mQ32)*pow2(mU32
         )*pow5(m32) - 90*pow2(msq2)*pow2(mU32)*pow5(m32) + 210*msq2*pow3(mQ32)
         *pow5(m32) - 19173*mU32*pow3(mQ32)*pow5(m32) - 743*mQ32*pow3(mU32)*
         pow5(m32) - 180*msq2*pow3(mU32)*pow5(m32) - 3988*pow4(mQ32)*pow5(m32)
         + 137*pow4(mU32)*pow5(m32) + 90*msq2*mU32*pow2(m32)*pow5(mQ32) - 90*
         m32*msq2*pow2(mU32)*pow5(mQ32) + 2250*pow2(m32)*pow2(mU32)*pow5(mQ32)
         - 30*msq2*pow3(m32)*pow5(mQ32) - 969*mU32*pow3(m32)*pow5(mQ32) - 1225*
         m32*pow3(mU32)*pow5(mQ32) + 30*msq2*pow3(mU32)*pow5(mQ32) + 364*pow4(
         m32)*pow5(mQ32) + 252*pow4(mU32)*pow5(mQ32) - 30*mQ32*msq2*pow2(m32)*
         pow5(mU32) + 120*m32*msq2*pow2(mQ32)*pow5(mU32) + 767*pow2(m32)*pow2(
         mQ32)*pow5(mU32) - 60*m32*mQ32*pow2(msq2)*pow5(mU32) + 90*pow2(m32)*
         pow2(msq2)*pow5(mU32) - 30*pow2(mQ32)*pow2(msq2)*pow5(mU32) - 299*mQ32
         *pow3(m32)*pow5(mU32) - 60*msq2*pow3(m32)*pow5(mU32) - 669*m32*pow3(
         mQ32)*pow5(mU32) - 30*msq2*pow3(mQ32)*pow5(mU32) + 2*pow4(m32)*pow5(
         mU32) + 191*pow4(mQ32)*pow5(mU32) - 120*mQ32*msq2*mU32*pow6(m32) + 60*
         msq2*pow2(mQ32)*pow6(m32) + 8831*mU32*pow2(mQ32)*pow6(m32) + 369*mQ32*
         pow2(mU32)*pow6(m32) + 60*msq2*pow2(mU32)*pow6(m32) + 5659*pow3(mQ32)*
         pow6(m32) - 475*pow3(mU32)*pow6(m32) - 1426*mU32*pow2(m32)*pow6(mQ32)
         + 1103*m32*pow2(mU32)*pow6(mQ32) + 343*pow3(m32)*pow6(mQ32) - 356*pow3
         (mU32)*pow6(mQ32) + 9*mQ32*pow2(m32)*pow6(mU32) - 6*m32*pow2(mQ32)*
         pow6(mU32) - 3*pow3(mQ32)*pow6(mU32) - 300*mQ32*mU32*pow7(m32) - 2588*
         pow2(mQ32)*pow7(m32) + 504*pow2(mU32)*pow7(m32) + 50*m32*mU32*pow7(
         mQ32) + 47*pow2(m32)*pow7(mQ32) - pow2(mU32)*pow7(mQ32) + 180*mQ32*
         pow8(m32) - 180*mU32*pow8(m32) - 9*m32*pow8(mQ32) - 3*mU32*pow8(mQ32))
         )/(3.*pow3(m32 - mU32)*pow4(m32 - mQ32)*pow4(mQ32 - mU32)) - (64*pow3(
         Xt)*(-60*m3*msq2*mU32*pow2(mQ32) + 60*m3*mQ32*msq2*pow2(mU32) + 165*m3
         *pow2(mQ32)*pow2(mU32) - 60*mQ32*msq2*mU32*pow3(m3) + 60*msq2*pow2(
         mQ32)*pow3(m3) - 101*mU32*pow2(mQ32)*pow3(m3) - 236*mQ32*pow2(mU32)*
         pow3(m3) - 103*m3*mU32*pow3(mQ32) + 129*pow3(m3)*pow3(mQ32) + 6*m3*
         mQ32*pow3(mU32) - 6*m3*pow4(mQ32) + 305*mQ32*mU32*pow5(m3) - 90*pow2(
         mQ32)*pow5(m3) - 27*pow2(mU32)*pow5(m3) - 77*mQ32*pow7(m3) + 13*mU32*
         pow7(m3) + 22*pow9(m3)))/(3.*(m32 - mU32)*pow2(m32 - mQ32)*pow3(mQ32 -
         mU32)) - (64*pow5(Xt)*(30*mQ32*pow11(m3) + 34*mU32*pow11(m3) - 30*m3*
         pow2(mQ32)*pow2(msq2)*pow2(mU32) + 60*mU32*pow2(mQ32)*pow2(msq2)*pow3(
         m3) - 30*msq2*pow2(mQ32)*pow2(mU32)*pow3(m3) - 30*mQ32*pow2(msq2)*pow2
         (mU32)*pow3(m3) + 90*m3*msq2*pow2(mU32)*pow3(mQ32) - 180*msq2*mU32*
         pow3(m3)*pow3(mQ32) - 608*pow2(mU32)*pow3(m3)*pow3(mQ32) - 30*m3*msq2*
         pow2(mQ32)*pow3(mU32) + 30*m3*mQ32*pow2(msq2)*pow3(mU32) + 30*mQ32*
         msq2*pow3(m3)*pow3(mU32) - 94*pow2(mQ32)*pow3(m3)*pow3(mU32) - 30*pow2
         (msq2)*pow3(m3)*pow3(mU32) + 67*m3*pow3(mQ32)*pow3(mU32) + 216*m3*pow2
         (mU32)*pow4(mQ32) - 396*mU32*pow3(m3)*pow4(mQ32) + 3*m3*pow2(mQ32)*
         pow4(mU32) - 3*mQ32*pow3(m3)*pow4(mU32) + 150*msq2*mU32*pow2(mQ32)*
         pow5(m3) - 30*mQ32*mU32*pow2(msq2)*pow5(m3) - 30*pow2(mQ32)*pow2(msq2)
         *pow5(m3) - 60*mQ32*msq2*pow2(mU32)*pow5(m3) + 427*pow2(mQ32)*pow2(
         mU32)*pow5(m3) + 60*pow2(msq2)*pow2(mU32)*pow5(m3) + 90*msq2*pow3(mQ32
         )*pow5(m3) + 923*mU32*pow3(mQ32)*pow5(m3) + 5*mQ32*pow3(mU32)*pow5(m3)
         + 225*pow4(mQ32)*pow5(m3) - 7*m3*mU32*pow5(mQ32) - 11*pow3(m3)*pow5(
         mQ32) + 3*m3*pow6(mQ32) + 30*mQ32*msq2*mU32*pow7(m3) - 90*msq2*pow2(
         mQ32)*pow7(m3) - 470*mU32*pow2(mQ32)*pow7(m3) + 30*mQ32*pow2(msq2)*
         pow7(m3) - 30*mU32*pow2(msq2)*pow7(m3) - 14*mQ32*pow2(mU32)*pow7(m3) -
         442*pow3(mQ32)*pow7(m3) + 38*pow3(mU32)*pow7(m3) - 36*mQ32*mU32*pow9(
         m3) + 179*pow2(mQ32)*pow9(m3) - 69*pow2(mU32)*pow9(m3)))/(3.*pow2(m32
         - mU32)*pow3(m32 - mQ32)*pow4(mQ32 - mU32)) - (32*Xt*(34*mQ32*pow11(m3
         ) - 34*mU32*pow11(m3) + 30*m3*pow2(mQ32)*pow2(msq2)*pow2(mU32) - 60*
         mU32*pow2(mQ32)*pow2(msq2)*pow3(m3) - 120*msq2*pow2(mQ32)*pow2(mU32)*
         pow3(m3) + 30*mQ32*pow2(msq2)*pow2(mU32)*pow3(m3) - 120*m3*msq2*pow2(
         mU32)*pow3(mQ32) + 240*msq2*mU32*pow3(m3)*pow3(mQ32) + 106*pow2(mU32)*
         pow3(m3)*pow3(mQ32) + 120*m3*msq2*pow2(mQ32)*pow3(mU32) - 30*m3*mQ32*
         pow2(msq2)*pow3(mU32) - 120*mQ32*msq2*pow3(m3)*pow3(mU32) - 392*pow2(
         mQ32)*pow3(m3)*pow3(mU32) + 30*pow2(msq2)*pow3(m3)*pow3(mU32) + 153*m3
         *pow3(mQ32)*pow3(mU32) - 176*m3*pow2(mU32)*pow4(mQ32) + 342*mU32*pow3(
         m3)*pow4(mQ32) + 6*m3*pow2(mQ32)*pow4(mU32) - 6*mQ32*pow3(m3)*pow4(
         mU32) - 120*msq2*mU32*pow2(mQ32)*pow5(m3) + 30*mQ32*mU32*pow2(msq2)*
         pow5(m3) + 30*pow2(mQ32)*pow2(msq2)*pow5(m3) + 240*mQ32*msq2*pow2(mU32
         )*pow5(m3) + 477*pow2(mQ32)*pow2(mU32)*pow5(m3) - 60*pow2(msq2)*pow2(
         mU32)*pow5(m3) - 120*msq2*pow3(mQ32)*pow5(m3) - 687*mU32*pow3(mQ32)*
         pow5(m3) + 325*mQ32*pow3(mU32)*pow5(m3) - 211*pow4(mQ32)*pow5(m3) + 4*
         m3*mU32*pow5(mQ32) + 14*pow3(m3)*pow5(mQ32) - 3*m3*pow6(mQ32) - 120*
         mQ32*msq2*mU32*pow7(m3) + 120*msq2*pow2(mQ32)*pow7(m3) + 266*mU32*pow2
         (mQ32)*pow7(m3) - 30*mQ32*pow2(msq2)*pow7(m3) + 30*mU32*pow2(msq2)*
         pow7(m3) - 684*mQ32*pow2(mU32)*pow7(m3) + 488*pow3(mQ32)*pow7(m3) - 6*
         pow3(mU32)*pow7(m3) + 349*mQ32*mU32*pow9(m3) - 402*pow2(mQ32)*pow9(m3)
         + 37*pow2(mU32)*pow9(m3)))/(3.*pow2(m32 - mU32)*pow2(mQ32 - mU32)*
         pow3(m32 - mQ32)) + log(mU32/MR2)*((-1280*m32*(mQ32 + mU32)*(2*m32*
         mQ32*mU32 + m32*pow2(mQ32) - 3*mU32*pow2(mQ32))*pow6(Xt))/(3.*(m32 -
         mU32)*pow2(m32 - mQ32)*pow5(mQ32 - mU32)) + (8*pow2(Xt)*(270*pow2(m32)
         *pow2(mQ32)*pow2(msq2)*pow2(mU32) - 210*mU32*pow2(mQ32)*pow2(msq2)*
         pow3(m32) - 180*msq2*pow2(mQ32)*pow2(mU32)*pow3(m32) - 90*mQ32*pow2(
         msq2)*pow2(mU32)*pow3(m32) - 90*mU32*pow2(m32)*pow2(msq2)*pow3(mQ32) -
         540*msq2*pow2(m32)*pow2(mU32)*pow3(mQ32) + 90*m32*pow2(msq2)*pow2(
         mU32)*pow3(mQ32) + 540*msq2*mU32*pow3(m32)*pow3(mQ32) + 30*pow2(msq2)*
         pow3(m32)*pow3(mQ32) - 3662*pow2(mU32)*pow3(m32)*pow3(mQ32) + 420*msq2
         *pow2(m32)*pow2(mQ32)*pow3(mU32) - 90*mQ32*pow2(m32)*pow2(msq2)*pow3(
         mU32) - 150*m32*pow2(mQ32)*pow2(msq2)*pow3(mU32) - 420*mQ32*msq2*pow3(
         m32)*pow3(mU32) - 1766*pow2(mQ32)*pow3(m32)*pow3(mU32) + 270*pow2(msq2
         )*pow3(m32)*pow3(mU32) + 180*m32*msq2*pow3(mQ32)*pow3(mU32) + 670*pow2
         (m32)*pow3(mQ32)*pow3(mU32) - 30*pow2(msq2)*pow3(mQ32)*pow3(mU32) -
         180*msq2*mU32*pow2(mQ32)*pow4(m32) + 210*mQ32*mU32*pow2(msq2)*pow4(m32
         ) + 60*pow2(mQ32)*pow2(msq2)*pow4(m32) + 540*mQ32*msq2*pow2(mU32)*pow4
         (m32) + 6721*pow2(mQ32)*pow2(mU32)*pow4(m32) - 270*pow2(msq2)*pow2(
         mU32)*pow4(m32) - 180*msq2*pow3(mQ32)*pow4(m32) + 1295*mU32*pow3(mQ32)
         *pow4(m32) + 1523*mQ32*pow3(mU32)*pow4(m32) - 180*msq2*pow3(mU32)*pow4
         (m32) + 950*pow2(m32)*pow2(mU32)*pow4(mQ32) + 688*mU32*pow3(m32)*pow4(
         mQ32) - 123*m32*pow3(mU32)*pow4(mQ32) - 272*pow4(m32)*pow4(mQ32) + 120
         *mQ32*msq2*pow2(m32)*pow4(mU32) - 180*m32*msq2*pow2(mQ32)*pow4(mU32) -
         538*pow2(m32)*pow2(mQ32)*pow4(mU32) + 60*m32*mQ32*pow2(msq2)*pow4(
         mU32) - 90*pow2(m32)*pow2(msq2)*pow4(mU32) + 30*pow2(mQ32)*pow2(msq2)*
         pow4(mU32) + 324*mQ32*pow3(m32)*pow4(mU32) + 60*msq2*pow3(m32)*pow4(
         mU32) + 412*m32*pow3(mQ32)*pow4(mU32) - 87*pow4(m32)*pow4(mU32) - 103*
         pow4(mQ32)*pow4(mU32) - 300*mQ32*msq2*mU32*pow5(m32) + 120*msq2*pow2(
         mQ32)*pow5(m32) - 4347*mU32*pow2(mQ32)*pow5(m32) - 90*mQ32*pow2(msq2)*
         pow5(m32) + 90*mU32*pow2(msq2)*pow5(m32) - 4987*mQ32*pow2(mU32)*pow5(
         m32) + 180*msq2*pow2(mU32)*pow5(m32) + 265*pow3(mQ32)*pow5(m32) - 375*
         pow3(mU32)*pow5(m32) - 520*mU32*pow2(m32)*pow5(mQ32) - 159*m32*pow2(
         mU32)*pow5(mQ32) - 92*pow3(m32)*pow5(mQ32) + 39*pow3(mU32)*pow5(mQ32)
         + 60*mQ32*msq2*pow6(m32) + 3562*mQ32*mU32*pow6(m32) - 60*msq2*mU32*
         pow6(m32) + 446*pow2(mQ32)*pow6(m32) + 1188*pow2(mU32)*pow6(m32) + 141
         *m32*mU32*pow6(mQ32) + 114*pow2(m32)*pow6(mQ32) - 3*pow2(mU32)*pow6(
         mQ32) - 554*mQ32*pow7(m32) - 842*mU32*pow7(m32) - 27*m32*pow7(mQ32) -
         9*mU32*pow7(mQ32) + 128*pow8(m32)))/(3.*pow2(mQ32 - mU32)*pow3(m32 -
         mU32)*pow4(m32 - mQ32)) - (4*(270*pow2(m32)*pow2(mQ32)*pow2(msq2)*pow2
         (mU32) - 210*mU32*pow2(mQ32)*pow2(msq2)*pow3(m32) - 180*msq2*pow2(mQ32
         )*pow2(mU32)*pow3(m32) - 90*mQ32*pow2(msq2)*pow2(mU32)*pow3(m32) - 90*
         mU32*pow2(m32)*pow2(msq2)*pow3(mQ32) - 540*msq2*pow2(m32)*pow2(mU32)*
         pow3(mQ32) + 90*m32*pow2(msq2)*pow2(mU32)*pow3(mQ32) + 540*msq2*mU32*
         pow3(m32)*pow3(mQ32) + 30*pow2(msq2)*pow3(m32)*pow3(mQ32) - 386*pow2(
         mU32)*pow3(m32)*pow3(mQ32) + 420*msq2*pow2(m32)*pow2(mQ32)*pow3(mU32)
         - 90*mQ32*pow2(m32)*pow2(msq2)*pow3(mU32) - 150*m32*pow2(mQ32)*pow2(
         msq2)*pow3(mU32) - 420*mQ32*msq2*pow3(m32)*pow3(mU32) - 2306*pow2(mQ32
         )*pow3(m32)*pow3(mU32) + 270*pow2(msq2)*pow3(m32)*pow3(mU32) + 180*m32
         *msq2*pow3(mQ32)*pow3(mU32) + 1002*pow2(m32)*pow3(mQ32)*pow3(mU32) -
         30*pow2(msq2)*pow3(mQ32)*pow3(mU32) - 180*msq2*mU32*pow2(mQ32)*pow4(
         m32) + 210*mQ32*mU32*pow2(msq2)*pow4(m32) + 60*pow2(mQ32)*pow2(msq2)*
         pow4(m32) + 540*mQ32*msq2*pow2(mU32)*pow4(m32) + 2585*pow2(mQ32)*pow2(
         mU32)*pow4(m32) - 270*pow2(msq2)*pow2(mU32)*pow4(m32) - 180*msq2*pow3(
         mQ32)*pow4(m32) + 199*mU32*pow3(mQ32)*pow4(m32) + 1699*mQ32*pow3(mU32)
         *pow4(m32) - 180*msq2*pow3(mU32)*pow4(m32) - 638*pow2(m32)*pow2(mU32)*
         pow4(mQ32) + 324*mU32*pow3(m32)*pow4(mQ32) - 19*m32*pow3(mU32)*pow4(
         mQ32) - 64*pow4(m32)*pow4(mQ32) + 120*mQ32*msq2*pow2(m32)*pow4(mU32) -
         180*m32*msq2*pow2(mQ32)*pow4(mU32) + 658*pow2(m32)*pow2(mQ32)*pow4(
         mU32) + 60*m32*mQ32*pow2(msq2)*pow4(mU32) - 90*pow2(m32)*pow2(msq2)*
         pow4(mU32) + 30*pow2(mQ32)*pow2(msq2)*pow4(mU32) - 420*mQ32*pow3(m32)*
         pow4(mU32) + 60*msq2*pow3(m32)*pow4(mU32) - 388*m32*pow3(mQ32)*pow4(
         mU32) + 61*pow4(m32)*pow4(mU32) + 97*pow4(mQ32)*pow4(mU32) - 300*mQ32*
         msq2*mU32*pow5(m32) + 120*msq2*pow2(mQ32)*pow5(m32) - 1619*mU32*pow2(
         mQ32)*pow5(m32) - 90*mQ32*pow2(msq2)*pow5(m32) + 90*mU32*pow2(msq2)*
         pow5(m32) - 2307*mQ32*pow2(mU32)*pow5(m32) + 180*msq2*pow2(mU32)*pow5(
         m32) - 247*pow3(mQ32)*pow5(m32) - 307*pow3(mU32)*pow5(m32) - 88*mU32*
         pow2(m32)*pow5(mQ32) + 317*m32*pow2(mU32)*pow5(mQ32) + 100*pow3(m32)*
         pow5(mQ32) - 101*pow3(mU32)*pow5(mQ32) + 60*mQ32*msq2*pow6(m32) + 1542
         *mQ32*mU32*pow6(m32) - 60*msq2*mU32*pow6(m32) + 670*pow2(mQ32)*pow6(
         m32) + 476*pow2(mU32)*pow6(m32) - 47*m32*mU32*pow6(mQ32) - 38*pow2(m32
         )*pow6(mQ32) + pow2(mU32)*pow6(mQ32) - 550*mQ32*pow7(m32) - 346*mU32*
         pow7(m32) + 9*m32*pow7(mQ32) + 3*mU32*pow7(mQ32) + 128*pow8(m32)))/(3.
         *(mQ32 - mU32)*pow3(m32 - mU32)*pow4(m32 - mQ32)) - (4*pow4(Xt)*(-300*
         pow2(mQ32)*pow2(msq2)*pow2(mU32)*pow3(m32) + 180*pow2(m32)*pow2(msq2)*
         pow2(mU32)*pow3(mQ32) - 180*mU32*pow2(msq2)*pow3(m32)*pow3(mQ32) + 360
         *msq2*pow2(mU32)*pow3(m32)*pow3(mQ32) + 180*pow2(m32)*pow2(mQ32)*pow2(
         msq2)*pow3(mU32) - 600*msq2*pow2(mQ32)*pow3(m32)*pow3(mU32) + 180*mQ32
         *pow2(msq2)*pow3(m32)*pow3(mU32) - 120*msq2*pow2(m32)*pow3(mQ32)*pow3(
         mU32) - 60*m32*pow2(msq2)*pow3(mQ32)*pow3(mU32) - 18228*pow3(m32)*pow3
         (mQ32)*pow3(mU32) + 270*mU32*pow2(mQ32)*pow2(msq2)*pow4(m32) + 360*
         msq2*pow2(mQ32)*pow2(mU32)*pow4(m32) - 60*mQ32*pow2(msq2)*pow2(mU32)*
         pow4(m32) - 360*msq2*mU32*pow3(mQ32)*pow4(m32) + 60*pow2(msq2)*pow3(
         mQ32)*pow4(m32) + 15072*pow2(mU32)*pow3(mQ32)*pow4(m32) + 360*mQ32*
         msq2*pow3(mU32)*pow4(m32) + 24682*pow2(mQ32)*pow3(mU32)*pow4(m32) -
         270*pow2(msq2)*pow3(mU32)*pow4(m32) - 90*mU32*pow2(m32)*pow2(msq2)*
         pow4(mQ32) - 540*msq2*pow2(m32)*pow2(mU32)*pow4(mQ32) + 90*m32*pow2(
         msq2)*pow2(mU32)*pow4(mQ32) + 540*msq2*mU32*pow3(m32)*pow4(mQ32) + 30*
         pow2(msq2)*pow3(m32)*pow4(mQ32) + 2920*pow2(mU32)*pow3(m32)*pow4(mQ32)
         + 180*m32*msq2*pow3(mU32)*pow4(mQ32) + 2474*pow2(m32)*pow3(mU32)*pow4
         (mQ32) - 30*pow2(msq2)*pow3(mU32)*pow4(mQ32) - 180*msq2*pow4(m32)*pow4
         (mQ32) - 5097*mU32*pow4(m32)*pow4(mQ32) + 540*msq2*pow2(m32)*pow2(mQ32
         )*pow4(mU32) - 180*mQ32*pow2(m32)*pow2(msq2)*pow4(mU32) - 90*m32*pow2(
         mQ32)*pow2(msq2)*pow4(mU32) - 360*mQ32*msq2*pow3(m32)*pow4(mU32) -
         6526*pow2(mQ32)*pow3(m32)*pow4(mU32) + 270*pow2(msq2)*pow3(m32)*pow4(
         mU32) + 5636*pow2(m32)*pow3(mQ32)*pow4(mU32) + 3160*mQ32*pow4(m32)*
         pow4(mU32) - 180*msq2*pow4(m32)*pow4(mU32) - 1633*m32*pow4(mQ32)*pow4(
         mU32) - 180*msq2*mU32*pow2(mQ32)*pow5(m32) - 90*pow2(mQ32)*pow2(msq2)*
         pow5(m32) - 120*mQ32*msq2*pow2(mU32)*pow5(m32) - 25608*pow2(mQ32)*pow2
         (mU32)*pow5(m32) + 90*pow2(msq2)*pow2(mU32)*pow5(m32) + 120*msq2*pow3(
         mQ32)*pow5(m32) - 3338*mU32*pow3(mQ32)*pow5(m32) - 13106*mQ32*pow3(
         mU32)*pow5(m32) + 180*msq2*pow3(mU32)*pow5(m32) + 1597*pow4(mQ32)*pow5
         (m32) - 521*pow4(mU32)*pow5(m32) - 4138*pow2(m32)*pow2(mU32)*pow5(mQ32
         ) + 3988*mU32*pow3(m32)*pow5(mQ32) + 1734*m32*pow3(mU32)*pow5(mQ32) -
         1036*pow4(m32)*pow5(mQ32) - 140*pow4(mU32)*pow5(mQ32) + 120*mQ32*msq2*
         pow2(m32)*pow5(mU32) - 180*m32*msq2*pow2(mQ32)*pow5(mU32) - 702*pow2(
         m32)*pow2(mQ32)*pow5(mU32) + 60*m32*mQ32*pow2(msq2)*pow5(mU32) - 90*
         pow2(m32)*pow2(msq2)*pow5(mU32) + 30*pow2(mQ32)*pow2(msq2)*pow5(mU32)
         + 348*mQ32*pow3(m32)*pow5(mU32) + 60*msq2*pow3(m32)*pow5(mU32) + 564*
         m32*pow3(mQ32)*pow5(mU32) - 61*pow4(m32)*pow5(mU32) - 141*pow4(mQ32)*
         pow5(mU32) + 60*msq2*pow2(mQ32)*pow6(m32) + 9070*mU32*pow2(mQ32)*pow6(
         m32) + 14710*mQ32*pow2(mU32)*pow6(m32) - 60*msq2*pow2(mU32)*pow6(m32)
         - 210*pow3(mQ32)*pow6(m32) + 2462*pow3(mU32)*pow6(m32) - 552*mU32*pow2
         (m32)*pow6(mQ32) + 96*m32*pow2(mU32)*pow6(mQ32) - 22*pow3(m32)*pow6(
         mQ32) - 2*pow3(mU32)*pow6(mQ32) - 5356*mQ32*mU32*pow7(m32) - 532*pow2(
         mQ32)*pow7(m32) - 3024*pow2(mU32)*pow7(m32) + 114*m32*mU32*pow7(mQ32)
         + 114*pow2(m32)*pow7(mQ32) - 12*pow2(mU32)*pow7(mQ32) + 124*mQ32*pow8(
         m32) + 1156*mU32*pow8(m32) - 27*m32*pow8(mQ32) - 9*mU32*pow8(mQ32)))/(
         3.*pow3(m32 - mU32)*pow4(m32 - mQ32)*pow4(mQ32 - mU32)) + (32*pow3(Xt)
         *(-43*mQ32*pow11(m3) - 33*mU32*pow11(m3) + 2*pow13(m3) - 60*m3*pow2(
         mQ32)*pow2(msq2)*pow2(mU32) + 120*mU32*pow2(mQ32)*pow2(msq2)*pow3(m3)
         + 120*msq2*pow2(mQ32)*pow2(mU32)*pow3(m3) - 60*mQ32*pow2(msq2)*pow2(
         mU32)*pow3(m3) + 120*m3*msq2*pow2(mU32)*pow3(mQ32) - 240*msq2*mU32*
         pow3(m3)*pow3(mQ32) - 1143*pow2(mU32)*pow3(m3)*pow3(mQ32) - 120*m3*
         msq2*pow2(mQ32)*pow3(mU32) + 60*m3*mQ32*pow2(msq2)*pow3(mU32) + 120*
         mQ32*msq2*pow3(m3)*pow3(mU32) + 123*pow2(mQ32)*pow3(m3)*pow3(mU32) -
         60*pow2(msq2)*pow3(m3)*pow3(mU32) + 3*m3*pow3(mQ32)*pow3(mU32) + 361*
         m3*pow2(mU32)*pow4(mQ32) - 288*mU32*pow3(m3)*pow4(mQ32) + 120*msq2*
         mU32*pow2(mQ32)*pow5(m3) - 60*mQ32*mU32*pow2(msq2)*pow5(m3) - 60*pow2(
         mQ32)*pow2(msq2)*pow5(m3) - 240*mQ32*msq2*pow2(mU32)*pow5(m3) + 959*
         pow2(mQ32)*pow2(mU32)*pow5(m3) + 120*pow2(msq2)*pow2(mU32)*pow5(m3) +
         120*msq2*pow3(mQ32)*pow5(m3) + 1201*mU32*pow3(mQ32)*pow5(m3) - 155*
         mQ32*pow3(mU32)*pow5(m3) + 201*pow4(mQ32)*pow5(m3) - 60*m3*mU32*pow5(
         mQ32) - 48*pow3(m3)*pow5(mQ32) + 18*m3*pow6(mQ32) + 120*mQ32*msq2*mU32
         *pow7(m3) - 120*msq2*pow2(mQ32)*pow7(m3) - 1127*mU32*pow2(mQ32)*pow7(
         m3) + 60*mQ32*pow2(msq2)*pow7(m3) - 60*mU32*pow2(msq2)*pow7(m3) - 233*
         mQ32*pow2(mU32)*pow7(m3) - 405*pow3(mQ32)*pow7(m3) + 61*pow3(mU32)*
         pow7(m3) + 403*mQ32*mU32*pow9(m3) + 243*pow2(mQ32)*pow9(m3) - 40*pow2(
         mU32)*pow9(m3)))/(3.*pow2(m32 - mU32)*pow3(m32 - mQ32)*pow3(mQ32 -
         mU32)) + (32*Xt*(20*mQ32*pow11(m3) - 20*mU32*pow11(m3) + 30*m3*pow2(
         mQ32)*pow2(msq2)*pow2(mU32) - 60*mU32*pow2(mQ32)*pow2(msq2)*pow3(m3) -
         60*msq2*pow2(mQ32)*pow2(mU32)*pow3(m3) + 30*mQ32*pow2(msq2)*pow2(mU32
         )*pow3(m3) - 60*m3*msq2*pow2(mU32)*pow3(mQ32) + 120*msq2*mU32*pow3(m3)
         *pow3(mQ32) - 42*pow2(mU32)*pow3(m3)*pow3(mQ32) + 60*m3*msq2*pow2(mQ32
         )*pow3(mU32) - 30*m3*mQ32*pow2(msq2)*pow3(mU32) - 60*mQ32*msq2*pow3(m3
         )*pow3(mU32) - 99*pow2(mQ32)*pow3(m3)*pow3(mU32) + 30*pow2(msq2)*pow3(
         m3)*pow3(mU32) + 23*m3*pow3(mQ32)*pow3(mU32) + 85*mU32*pow3(m3)*pow4(
         mQ32) - 60*msq2*mU32*pow2(mQ32)*pow5(m3) + 30*mQ32*mU32*pow2(msq2)*
         pow5(m3) + 30*pow2(mQ32)*pow2(msq2)*pow5(m3) + 120*mQ32*msq2*pow2(mU32
         )*pow5(m3) + 278*pow2(mQ32)*pow2(mU32)*pow5(m3) - 60*pow2(msq2)*pow2(
         mU32)*pow5(m3) - 60*msq2*pow3(mQ32)*pow5(m3) - 219*mU32*pow3(mQ32)*
         pow5(m3) + 67*mQ32*pow3(mU32)*pow5(m3) - 30*pow4(mQ32)*pow5(m3) - 10*
         m3*mU32*pow5(mQ32) - 8*pow3(m3)*pow5(mQ32) + 3*m3*pow6(mQ32) - 60*mQ32
         *msq2*mU32*pow7(m3) + 60*msq2*pow2(mQ32)*pow7(m3) + 65*mU32*pow2(mQ32)
         *pow7(m3) - 30*mQ32*pow2(msq2)*pow7(m3) + 30*mU32*pow2(msq2)*pow7(m3)
         - 244*mQ32*pow2(mU32)*pow7(m3) + 122*pow3(mQ32)*pow7(m3) - 7*pow3(mU32
         )*pow7(m3) + 51*mQ32*mU32*pow9(m3) - 91*pow2(mQ32)*pow9(m3) + 56*pow2(
         mU32)*pow9(m3)))/(3.*pow2(m32 - mU32)*pow2(mQ32 - mU32)*pow3(m32 -
         mQ32)) + (32*pow5(Xt)*(-32*mQ32*mU32*pow11(m3) + 18*pow11(m3)*pow2(
         mQ32) + 14*pow11(m3)*pow2(mU32) - 30*pow2(mQ32)*pow2(msq2)*pow2(mU32)*
         pow3(m3) + 30*m3*pow2(msq2)*pow2(mU32)*pow3(mQ32) - 60*mU32*pow2(msq2)
         *pow3(m3)*pow3(mQ32) + 60*msq2*pow2(mU32)*pow3(m3)*pow3(mQ32) - 120*
         msq2*pow2(mQ32)*pow3(m3)*pow3(mU32) + 60*mQ32*pow2(msq2)*pow3(m3)*pow3
         (mU32) + 761*pow3(m3)*pow3(mQ32)*pow3(mU32) - 60*m3*msq2*pow2(mU32)*
         pow4(mQ32) + 120*msq2*mU32*pow3(m3)*pow4(mQ32) + 879*pow2(mU32)*pow3(
         m3)*pow4(mQ32) - 261*m3*pow3(mU32)*pow4(mQ32) + 60*m3*msq2*pow2(mQ32)*
         pow4(mU32) - 30*m3*mQ32*pow2(msq2)*pow4(mU32) - 60*mQ32*msq2*pow3(m3)*
         pow4(mU32) + 25*pow2(mQ32)*pow3(m3)*pow4(mU32) + 30*pow2(msq2)*pow3(m3
         )*pow4(mU32) - 41*m3*pow3(mQ32)*pow4(mU32) + 60*mU32*pow2(mQ32)*pow2(
         msq2)*pow5(m3) + 60*msq2*pow2(mQ32)*pow2(mU32)*pow5(m3) - 30*mQ32*pow2
         (msq2)*pow2(mU32)*pow5(m3) - 120*msq2*mU32*pow3(mQ32)*pow5(m3) + 30*
         pow2(msq2)*pow3(mQ32)*pow5(m3) - 1267*pow2(mU32)*pow3(mQ32)*pow5(m3) +
         120*mQ32*msq2*pow3(mU32)*pow5(m3) - 681*pow2(mQ32)*pow3(mU32)*pow5(m3
         ) - 60*pow2(msq2)*pow3(mU32)*pow5(m3) - 60*msq2*pow4(mQ32)*pow5(m3) -
         823*mU32*pow4(mQ32)*pow5(m3) + 23*mQ32*pow4(mU32)*pow5(m3) - 190*m3*
         pow2(mU32)*pow5(mQ32) + 231*mU32*pow3(m3)*pow5(mQ32) - 132*pow5(m3)*
         pow5(mQ32) + 21*m3*mU32*pow6(mQ32) + 24*pow3(m3)*pow6(mQ32) - 30*pow2(
         mQ32)*pow2(msq2)*pow7(m3) - 60*mQ32*msq2*pow2(mU32)*pow7(m3) + 719*
         pow2(mQ32)*pow2(mU32)*pow7(m3) + 30*pow2(msq2)*pow2(mU32)*pow7(m3) +
         60*msq2*pow3(mQ32)*pow7(m3) + 753*mU32*pow3(mQ32)*pow7(m3) + 207*mQ32*
         pow3(mU32)*pow7(m3) + 264*pow4(mQ32)*pow7(m3) - 23*pow4(mU32)*pow7(m3)
         - 9*m3*pow7(mQ32) - 182*mU32*pow2(mQ32)*pow9(m3) - 155*mQ32*pow2(mU32
         )*pow9(m3) - 149*pow3(mQ32)*pow9(m3) + 6*pow3(mU32)*pow9(m3)))/(3.*
         pow2(m32 - mU32)*pow3(m32 - mQ32)*pow5(mQ32 - mU32)))) + pow2(log(
         m32/MR2))*(pow2(Xt)*((-1792*pow3(m32))/(3.*pow2(m32 - mQ32)*pow2(m32 -
         mU32)) + (64*(-m32 + mQ32 + mU32)*pow2(m32)*(2 + (8*(pow2(m32)*pow2(
         mQ32) + pow2(m32)*pow2(mU32) - 2*mQ32*pow3(m32) - 2*mU32*pow3(m32) + 2
         *pow4(m32)))/(3.*pow2(m32 - mQ32)*pow2(m32 - mU32))))/(mQ32*(-m32 +
         mQ32)*(m32 - mU32)*mU32)) - (2560*pow3(m32)*pow6(Xt))/(3.*pow2(m32 -
         mQ32)*pow2(m32 - mU32)*pow2(mQ32 - mU32)) - (2048*(m32 - mQ32 - mU32)*
         pow3(Xt)*pow7(m3))/(3.*mQ32*mU32*pow2(m32 - mQ32)*pow2(m32 - mU32)) +
         (128*pow5(Xt)*(7*mQ32*mU32*pow11(m3) + 8*mQ32*pow13(m3) + 8*mU32*pow13
         (m3) - 16*pow11(m3)*pow2(mQ32) - 16*pow11(m3)*pow2(mU32) + 87*pow3(m3)
         *pow3(mQ32)*pow3(mU32) - 134*pow2(mU32)*pow3(mQ32)*pow5(m3) - 134*pow2
         (mQ32)*pow3(mU32)*pow5(m3) + 188*pow2(mQ32)*pow2(mU32)*pow7(m3) + 47*
         mU32*pow3(mQ32)*pow7(m3) + 47*mQ32*pow3(mU32)*pow7(m3) - 54*mU32*pow2(
         mQ32)*pow9(m3) - 54*mQ32*pow2(mU32)*pow9(m3) + 8*pow3(mQ32)*pow9(m3) +
         8*pow3(mU32)*pow9(m3)))/(3.*mQ32*mU32*pow2(mQ32 - mU32)*pow3(m32 -
         mQ32)*pow3(m32 - mU32)) + (64*Xt*(31*mQ32*mU32*pow11(m3) + 16*mQ32*
         pow13(m3) + 16*mU32*pow13(m3) - 32*pow11(m3)*pow2(mQ32) - 32*pow11(m3)
         *pow2(mU32) - 81*pow3(m3)*pow3(mQ32)*pow3(mU32) + 130*pow2(mU32)*pow3(
         mQ32)*pow5(m3) + 130*pow2(mQ32)*pow3(mU32)*pow5(m3) - 196*pow2(mQ32)*
         pow2(mU32)*pow7(m3) - 25*mU32*pow3(mQ32)*pow7(m3) - 25*mQ32*pow3(mU32)
         *pow7(m3) + 18*mU32*pow2(mQ32)*pow9(m3) + 18*mQ32*pow2(mU32)*pow9(m3)
         + 16*pow3(mQ32)*pow9(m3) + 16*pow3(mU32)*pow9(m3)))/(3.*mQ32*mU32*pow3
         (m32 - mQ32)*pow3(m32 - mU32)) - (16*pow4(Xt)*(1712*pow3(mQ32)*pow3(
         mU32)*pow4(m32) - 88*pow3(m32)*pow3(mU32)*pow4(mQ32) + 572*pow2(mU32)*
         pow4(m32)*pow4(mQ32) - 88*pow3(m32)*pow3(mQ32)*pow4(mU32) + 572*pow2(
         mQ32)*pow4(m32)*pow4(mU32) - 548*pow2(m32)*pow4(mQ32)*pow4(mU32) -
         2012*pow2(mU32)*pow3(mQ32)*pow5(m32) - 2012*pow2(mQ32)*pow3(mU32)*pow5
         (m32) + 4*mU32*pow4(mQ32)*pow5(m32) + 4*mQ32*pow4(mU32)*pow5(m32) - 70
         *pow2(mU32)*pow3(m32)*pow5(mQ32) - 49*pow2(m32)*pow3(mU32)*pow5(mQ32)
         - 25*mU32*pow4(m32)*pow5(mQ32) + 156*m32*pow4(mU32)*pow5(mQ32) + 28*
         pow5(m32)*pow5(mQ32) - 70*pow2(mQ32)*pow3(m32)*pow5(mU32) - 49*pow2(
         m32)*pow3(mQ32)*pow5(mU32) - 25*mQ32*pow4(m32)*pow5(mU32) + 156*m32*
         pow4(mQ32)*pow5(mU32) + 28*pow5(m32)*pow5(mU32) - 36*pow5(mQ32)*pow5(
         mU32) + 1932*pow2(mQ32)*pow2(mU32)*pow6(m32) + 277*mU32*pow3(mQ32)*
         pow6(m32) + 277*mQ32*pow3(mU32)*pow6(m32) - 112*pow4(mQ32)*pow6(m32) -
         112*pow4(mU32)*pow6(m32) - 294*mU32*pow2(mQ32)*pow7(m32) - 294*mQ32*
         pow2(mU32)*pow7(m32) + 184*pow3(mQ32)*pow7(m32) + 184*pow3(mU32)*pow7(
         m32) - 2*mQ32*mU32*pow8(m32) - 144*pow2(mQ32)*pow8(m32) - 144*pow2(
         mU32)*pow8(m32) + 44*mQ32*pow9(m32) + 44*mU32*pow9(m32)))/(3.*mQ32*
         mU32*pow2(mQ32 - mU32)*pow4(m32 - mQ32)*pow4(m32 - mU32)) - (8*(-2048*
         pow3(mQ32)*pow3(mU32)*pow4(m32) + 272*pow3(m32)*pow3(mU32)*pow4(mQ32)
         - 908*pow2(mU32)*pow4(m32)*pow4(mQ32) + 272*pow3(m32)*pow3(mQ32)*pow4(
         mU32) - 908*pow2(mQ32)*pow4(m32)*pow4(mU32) + 452*pow2(m32)*pow4(mQ32)
         *pow4(mU32) + 2588*pow2(mU32)*pow3(mQ32)*pow5(m32) + 2588*pow2(mQ32)*
         pow3(mU32)*pow5(m32) + 280*mU32*pow4(mQ32)*pow5(m32) + 280*mQ32*pow4(
         mU32)*pow5(m32) + 134*pow2(mU32)*pow3(m32)*pow5(mQ32) + 25*pow2(m32)*
         pow3(mU32)*pow5(mQ32) - 39*mU32*pow4(m32)*pow5(mQ32) - 144*m32*pow4(
         mU32)*pow5(mQ32) + 56*pow5(m32)*pow5(mQ32) + 134*pow2(mQ32)*pow3(m32)*
         pow5(mU32) + 25*pow2(m32)*pow3(mQ32)*pow5(mU32) - 39*mQ32*pow4(m32)*
         pow5(mU32) - 144*m32*pow4(mQ32)*pow5(mU32) + 56*pow5(m32)*pow5(mU32) +
         36*pow5(mQ32)*pow5(mU32) - 2636*pow2(mQ32)*pow2(mU32)*pow6(m32) - 797
         *mU32*pow3(mQ32)*pow6(m32) - 797*mQ32*pow3(mU32)*pow6(m32) - 224*pow4(
         mQ32)*pow6(m32) - 224*pow4(mU32)*pow6(m32) + 838*mU32*pow2(mQ32)*pow7(
         m32) + 838*mQ32*pow2(mU32)*pow7(m32) + 368*pow3(mQ32)*pow7(m32) + 368*
         pow3(mU32)*pow7(m32) - 302*mQ32*mU32*pow8(m32) - 288*pow2(mQ32)*pow8(
         m32) - 288*pow2(mU32)*pow8(m32) + 88*mQ32*pow9(m32) + 88*mU32*pow9(m32
         )))/(3.*mQ32*mU32*pow4(m32 - mQ32)*pow4(m32 - mU32)) + log(mU32/MR2)*(
         (-1280*(mQ32 + mU32)*pow3(m32)*pow6(Xt))/(3.*pow2(m32 - mQ32)*pow2(m32
         - mU32)*pow3(mQ32 - mU32)) + (16*pow2(Xt)*(-192*pow2(mU32)*pow3(m32)*
         pow3(mQ32) - 288*pow2(mQ32)*pow3(m32)*pow3(mU32) - 452*pow2(m32)*pow3(
         mQ32)*pow3(mU32) + 1648*pow2(mQ32)*pow2(mU32)*pow4(m32) + 652*mU32*
         pow3(mQ32)*pow4(m32) + 844*mQ32*pow3(mU32)*pow4(m32) - pow2(m32)*pow2(
         mU32)*pow4(mQ32) - 158*mU32*pow3(m32)*pow4(mQ32) + 144*m32*pow3(mU32)*
         pow4(mQ32) + 55*pow4(m32)*pow4(mQ32) - pow2(m32)*pow2(mQ32)*pow4(mU32)
         - 158*mQ32*pow3(m32)*pow4(mU32) + 144*m32*pow3(mQ32)*pow4(mU32) + 55*
         pow4(m32)*pow4(mU32) - 36*pow4(mQ32)*pow4(mU32) - 1660*mU32*pow2(mQ32)
         *pow5(m32) - 1948*mQ32*pow2(mU32)*pow5(m32) - 168*pow3(mQ32)*pow5(m32)
         - 264*pow3(mU32)*pow5(m32) + 1484*mQ32*mU32*pow6(m32) + 325*pow2(mQ32
         )*pow6(m32) + 517*pow2(mU32)*pow6(m32) - 238*mQ32*pow7(m32) - 334*mU32
         *pow7(m32) + 30*pow8(m32)))/(3.*(mQ32 - mU32)*pow4(m32 - mQ32)*pow4(
         m32 - mU32)) - (8*(-272*pow2(mU32)*pow3(m32)*pow3(mQ32) - 368*pow2(
         mQ32)*pow3(m32)*pow3(mU32) - 452*pow2(m32)*pow3(mQ32)*pow3(mU32) +
         2240*pow2(mQ32)*pow2(mU32)*pow4(m32) + 812*mU32*pow3(mQ32)*pow4(m32) +
         1004*mQ32*pow3(mU32)*pow4(m32) - 25*pow2(m32)*pow2(mU32)*pow4(mQ32) -
         110*mU32*pow3(m32)*pow4(mQ32) + 144*m32*pow3(mU32)*pow4(mQ32) - pow4(
         m32)*pow4(mQ32) - pow2(m32)*pow2(mQ32)*pow4(mU32) - 158*mQ32*pow3(m32)
         *pow4(mU32) + 144*m32*pow3(mQ32)*pow4(mU32) + 55*pow4(m32)*pow4(mU32)
         - 36*pow4(mQ32)*pow4(mU32) - 2604*mU32*pow2(mQ32)*pow5(m32) - 2796*
         mQ32*pow2(mU32)*pow5(m32) - 120*pow3(mQ32)*pow5(m32) - 344*pow3(mU32)*
         pow5(m32) + 2700*mQ32*mU32*pow6(m32) + 565*pow2(mQ32)*pow6(m32) + 877*
         pow2(mU32)*pow6(m32) - 638*mQ32*pow7(m32) - 814*mU32*pow7(m32) + 198*
         pow8(m32)))/(3.*pow4(m32 - mQ32)*pow4(m32 - mU32)) + (64*pow5(Xt)*(39*
         mQ32*pow11(m3) + 39*mU32*pow11(m3) + 87*pow2(mU32)*pow3(m3)*pow3(mQ32)
         + 87*pow2(mQ32)*pow3(m3)*pow3(mU32) - 284*pow2(mQ32)*pow2(mU32)*pow5(
         m3) - 142*mU32*pow3(mQ32)*pow5(m3) - 142*mQ32*pow3(mU32)*pow5(m3) +
         283*mU32*pow2(mQ32)*pow7(m3) + 283*mQ32*pow2(mU32)*pow7(m3) + 63*pow3(
         mQ32)*pow7(m3) + 63*pow3(mU32)*pow7(m3) - 188*mQ32*mU32*pow9(m3) - 94*
         pow2(mQ32)*pow9(m3) - 94*pow2(mU32)*pow9(m3)))/(3.*pow3(m32 - mQ32)*
         pow3(m32 - mU32)*pow3(mQ32 - mU32)) + (64*Xt*(-77*mQ32*pow11(m3) - 79*
         mU32*pow11(m3) + 44*pow13(m3) + 75*pow2(mU32)*pow3(m3)*pow3(mQ32) - 87
         *pow2(mQ32)*pow3(m3)*pow3(mU32) + 36*pow2(mQ32)*pow2(mU32)*pow5(m3) -
         118*mU32*pow3(mQ32)*pow5(m3) + 142*mQ32*pow3(mU32)*pow5(m3) + 101*mU32
         *pow2(mQ32)*pow7(m3) - 209*mQ32*pow2(mU32)*pow7(m3) + 19*pow3(mQ32)*
         pow7(m3) - 63*pow3(mU32)*pow7(m3) + 72*mQ32*mU32*pow9(m3) + 22*pow2(
         mQ32)*pow9(m3) + 122*pow2(mU32)*pow9(m3)))/(3.*(mQ32 - mU32)*pow3(m32
         - mQ32)*pow3(m32 - mU32)) + (128*pow3(Xt)*(282*mQ32*mU32*pow11(m3) -
         100*mQ32*pow13(m3) - 100*mU32*pow13(m3) + 44*pow15(m3) + 45*pow11(m3)*
         pow2(mQ32) + 45*pow11(m3)*pow2(mU32) + 186*pow3(m3)*pow3(mQ32)*pow3(
         mU32) - 87*pow2(mU32)*pow3(m3)*pow4(mQ32) - 87*pow2(mQ32)*pow3(m3)*
         pow4(mU32) - 178*pow2(mU32)*pow3(mQ32)*pow5(m3) - 178*pow2(mQ32)*pow3(
         mU32)*pow5(m3) + 142*mU32*pow4(mQ32)*pow5(m3) + 142*mQ32*pow4(mU32)*
         pow5(m3) + 422*pow2(mQ32)*pow2(mU32)*pow7(m3) - 42*mU32*pow3(mQ32)*
         pow7(m3) - 42*mQ32*pow3(mU32)*pow7(m3) - 63*pow4(mQ32)*pow7(m3) - 63*
         pow4(mU32)*pow7(m3) - 250*mU32*pow2(mQ32)*pow9(m3) - 250*mQ32*pow2(
         mU32)*pow9(m3) + 66*pow3(mQ32)*pow9(m3) + 66*pow3(mU32)*pow9(m3)))/(3.
         *pow3(m32 - mQ32)*pow3(m32 - mU32)*pow3(mQ32 - mU32)) - (8*pow4(Xt)*(
         928*pow3(m32)*pow3(mQ32)*pow3(mU32) - 292*pow2(mU32)*pow3(mQ32)*pow4(
         m32) - 292*pow2(mQ32)*pow3(mU32)*pow4(m32) - 846*pow2(mU32)*pow3(m32)*
         pow4(mQ32) - 453*pow2(m32)*pow3(mU32)*pow4(mQ32) + 1699*mU32*pow4(m32)
         *pow4(mQ32) - 846*pow2(mQ32)*pow3(m32)*pow4(mU32) - 453*pow2(m32)*pow3
         (mQ32)*pow4(mU32) + 1699*mQ32*pow4(m32)*pow4(mU32) + 288*m32*pow4(mQ32
         )*pow4(mU32) + 5224*pow2(mQ32)*pow2(mU32)*pow5(m32) - 868*mU32*pow3(
         mQ32)*pow5(m32) - 868*mQ32*pow3(mU32)*pow5(m32) - 664*pow4(mQ32)*pow5(
         m32) - 664*pow4(mU32)*pow5(m32) - pow2(m32)*pow2(mU32)*pow5(mQ32) -
         158*mU32*pow3(m32)*pow5(mQ32) + 144*m32*pow3(mU32)*pow5(mQ32) + 55*
         pow4(m32)*pow5(mQ32) - 36*pow4(mU32)*pow5(mQ32) - pow2(m32)*pow2(mQ32)
         *pow5(mU32) - 158*mQ32*pow3(m32)*pow5(mU32) + 144*m32*pow3(mQ32)*pow5(
         mU32) + 55*pow4(m32)*pow5(mU32) - 36*pow4(mQ32)*pow5(mU32) - 5903*mU32
         *pow2(mQ32)*pow6(m32) - 5903*mQ32*pow2(mU32)*pow6(m32) + 549*pow3(mQ32
         )*pow6(m32) + 549*pow3(mU32)*pow6(m32) + 8004*mQ32*mU32*pow7(m32) +
         1826*pow2(mQ32)*pow7(m32) + 1826*pow2(mU32)*pow7(m32) - 2786*mQ32*pow8
         (m32) - 2786*mU32*pow8(m32) + 1024*pow9(m32)))/(3.*pow3(mQ32 - mU32)*
         pow4(m32 - mQ32)*pow4(m32 - mU32))) + log(mQ32/MR2)*((1280*(mQ32 +
         mU32)*pow3(m32)*pow6(Xt))/(3.*pow2(m32 - mQ32)*pow2(m32 - mU32)*pow3(
         mQ32 - mU32)) - (16*pow2(Xt)*(-288*pow2(mU32)*pow3(m32)*pow3(mQ32) -
         192*pow2(mQ32)*pow3(m32)*pow3(mU32) - 452*pow2(m32)*pow3(mQ32)*pow3(
         mU32) + 1648*pow2(mQ32)*pow2(mU32)*pow4(m32) + 844*mU32*pow3(mQ32)*
         pow4(m32) + 652*mQ32*pow3(mU32)*pow4(m32) - pow2(m32)*pow2(mU32)*pow4(
         mQ32) - 158*mU32*pow3(m32)*pow4(mQ32) + 144*m32*pow3(mU32)*pow4(mQ32)
         + 55*pow4(m32)*pow4(mQ32) - pow2(m32)*pow2(mQ32)*pow4(mU32) - 158*mQ32
         *pow3(m32)*pow4(mU32) + 144*m32*pow3(mQ32)*pow4(mU32) + 55*pow4(m32)*
         pow4(mU32) - 36*pow4(mQ32)*pow4(mU32) - 1948*mU32*pow2(mQ32)*pow5(m32)
         - 1660*mQ32*pow2(mU32)*pow5(m32) - 264*pow3(mQ32)*pow5(m32) - 168*
         pow3(mU32)*pow5(m32) + 1484*mQ32*mU32*pow6(m32) + 517*pow2(mQ32)*pow6(
         m32) + 325*pow2(mU32)*pow6(m32) - 334*mQ32*pow7(m32) - 238*mU32*pow7(
         m32) + 30*pow8(m32)))/(3.*(mQ32 - mU32)*pow4(m32 - mQ32)*pow4(m32 -
         mU32)) - (8*(-368*pow2(mU32)*pow3(m32)*pow3(mQ32) - 272*pow2(mQ32)*
         pow3(m32)*pow3(mU32) - 452*pow2(m32)*pow3(mQ32)*pow3(mU32) + 2240*pow2
         (mQ32)*pow2(mU32)*pow4(m32) + 1004*mU32*pow3(mQ32)*pow4(m32) + 812*
         mQ32*pow3(mU32)*pow4(m32) - pow2(m32)*pow2(mU32)*pow4(mQ32) - 158*mU32
         *pow3(m32)*pow4(mQ32) + 144*m32*pow3(mU32)*pow4(mQ32) + 55*pow4(m32)*
         pow4(mQ32) - 25*pow2(m32)*pow2(mQ32)*pow4(mU32) - 110*mQ32*pow3(m32)*
         pow4(mU32) + 144*m32*pow3(mQ32)*pow4(mU32) - pow4(m32)*pow4(mU32) - 36
         *pow4(mQ32)*pow4(mU32) - 2796*mU32*pow2(mQ32)*pow5(m32) - 2604*mQ32*
         pow2(mU32)*pow5(m32) - 344*pow3(mQ32)*pow5(m32) - 120*pow3(mU32)*pow5(
         m32) + 2700*mQ32*mU32*pow6(m32) + 877*pow2(mQ32)*pow6(m32) + 565*pow2(
         mU32)*pow6(m32) - 814*mQ32*pow7(m32) - 638*mU32*pow7(m32) + 198*pow8(
         m32)))/(3.*pow4(m32 - mQ32)*pow4(m32 - mU32)) - (64*pow5(Xt)*(39*mQ32*
         pow11(m3) + 39*mU32*pow11(m3) + 87*pow2(mU32)*pow3(m3)*pow3(mQ32) + 87
         *pow2(mQ32)*pow3(m3)*pow3(mU32) - 284*pow2(mQ32)*pow2(mU32)*pow5(m3) -
         142*mU32*pow3(mQ32)*pow5(m3) - 142*mQ32*pow3(mU32)*pow5(m3) + 283*
         mU32*pow2(mQ32)*pow7(m3) + 283*mQ32*pow2(mU32)*pow7(m3) + 63*pow3(mQ32
         )*pow7(m3) + 63*pow3(mU32)*pow7(m3) - 188*mQ32*mU32*pow9(m3) - 94*pow2
         (mQ32)*pow9(m3) - 94*pow2(mU32)*pow9(m3)))/(3.*pow3(m32 - mQ32)*pow3(
         m32 - mU32)*pow3(mQ32 - mU32)) - (64*Xt*(-79*mQ32*pow11(m3) - 77*mU32*
         pow11(m3) + 44*pow13(m3) - 87*pow2(mU32)*pow3(m3)*pow3(mQ32) + 75*pow2
         (mQ32)*pow3(m3)*pow3(mU32) + 36*pow2(mQ32)*pow2(mU32)*pow5(m3) + 142*
         mU32*pow3(mQ32)*pow5(m3) - 118*mQ32*pow3(mU32)*pow5(m3) - 209*mU32*
         pow2(mQ32)*pow7(m3) + 101*mQ32*pow2(mU32)*pow7(m3) - 63*pow3(mQ32)*
         pow7(m3) + 19*pow3(mU32)*pow7(m3) + 72*mQ32*mU32*pow9(m3) + 122*pow2(
         mQ32)*pow9(m3) + 22*pow2(mU32)*pow9(m3)))/(3.*(mQ32 - mU32)*pow3(m32 -
         mQ32)*pow3(m32 - mU32)) - (128*pow3(Xt)*(282*mQ32*mU32*pow11(m3) -
         100*mQ32*pow13(m3) - 100*mU32*pow13(m3) + 44*pow15(m3) + 45*pow11(m3)*
         pow2(mQ32) + 45*pow11(m3)*pow2(mU32) + 186*pow3(m3)*pow3(mQ32)*pow3(
         mU32) - 87*pow2(mU32)*pow3(m3)*pow4(mQ32) - 87*pow2(mQ32)*pow3(m3)*
         pow4(mU32) - 178*pow2(mU32)*pow3(mQ32)*pow5(m3) - 178*pow2(mQ32)*pow3(
         mU32)*pow5(m3) + 142*mU32*pow4(mQ32)*pow5(m3) + 142*mQ32*pow4(mU32)*
         pow5(m3) + 422*pow2(mQ32)*pow2(mU32)*pow7(m3) - 42*mU32*pow3(mQ32)*
         pow7(m3) - 42*mQ32*pow3(mU32)*pow7(m3) - 63*pow4(mQ32)*pow7(m3) - 63*
         pow4(mU32)*pow7(m3) - 250*mU32*pow2(mQ32)*pow9(m3) - 250*mQ32*pow2(
         mU32)*pow9(m3) + 66*pow3(mQ32)*pow9(m3) + 66*pow3(mU32)*pow9(m3)))/(3.
         *pow3(m32 - mQ32)*pow3(m32 - mU32)*pow3(mQ32 - mU32)) + (8*pow4(Xt)*(
         928*pow3(m32)*pow3(mQ32)*pow3(mU32) - 292*pow2(mU32)*pow3(mQ32)*pow4(
         m32) - 292*pow2(mQ32)*pow3(mU32)*pow4(m32) - 846*pow2(mU32)*pow3(m32)*
         pow4(mQ32) - 453*pow2(m32)*pow3(mU32)*pow4(mQ32) + 1699*mU32*pow4(m32)
         *pow4(mQ32) - 846*pow2(mQ32)*pow3(m32)*pow4(mU32) - 453*pow2(m32)*pow3
         (mQ32)*pow4(mU32) + 1699*mQ32*pow4(m32)*pow4(mU32) + 288*m32*pow4(mQ32
         )*pow4(mU32) + 5224*pow2(mQ32)*pow2(mU32)*pow5(m32) - 868*mU32*pow3(
         mQ32)*pow5(m32) - 868*mQ32*pow3(mU32)*pow5(m32) - 664*pow4(mQ32)*pow5(
         m32) - 664*pow4(mU32)*pow5(m32) - pow2(m32)*pow2(mU32)*pow5(mQ32) -
         158*mU32*pow3(m32)*pow5(mQ32) + 144*m32*pow3(mU32)*pow5(mQ32) + 55*
         pow4(m32)*pow5(mQ32) - 36*pow4(mU32)*pow5(mQ32) - pow2(m32)*pow2(mQ32)
         *pow5(mU32) - 158*mQ32*pow3(m32)*pow5(mU32) + 144*m32*pow3(mQ32)*pow5(
         mU32) + 55*pow4(m32)*pow5(mU32) - 36*pow4(mQ32)*pow5(mU32) - 5903*mU32
         *pow2(mQ32)*pow6(m32) - 5903*mQ32*pow2(mU32)*pow6(m32) + 549*pow3(mQ32
         )*pow6(m32) + 549*pow3(mU32)*pow6(m32) + 8004*mQ32*mU32*pow7(m32) +
         1826*pow2(mQ32)*pow7(m32) + 1826*pow2(mU32)*pow7(m32) - 2786*mQ32*pow8
         (m32) - 2786*mU32*pow8(m32) + 1024*pow9(m32)))/(3.*pow3(mQ32 - mU32)*
         pow4(m32 - mQ32)*pow4(m32 - mU32)))) + log(m32/MR2)*(pow3(Xt)*((-512*
         m3*(2 + (8*(pow2(m32)*pow2(mQ32) + pow2(m32)*pow2(mU32) - 2*mQ32*pow3(
         m32) - 2*mU32*pow3(m32) + 2*pow4(m32)))/(3.*pow2(m32 - mQ32)*pow2(m32
         - mU32))))/pow2(mQ32 - mU32) + (2048*pow5(m3))/(3.*(m32 - mQ32)*mQ32*(
         m32 - mU32)*mU32)) - (64*pow2(Xt)*(-7*mU32*pow2(m32)*pow2(mQ32) - 7*
         mQ32*pow2(m32)*pow2(mU32) + 6*m32*pow2(mQ32)*pow2(mU32) + 25*mQ32*mU32
         *pow3(m32) + 17*pow2(mQ32)*pow3(m32) + 17*pow2(mU32)*pow3(m32) - 42*
         mQ32*pow4(m32) - 42*mU32*pow4(m32) + 33*pow5(m32)))/(3.*mQ32*mU32*pow2
         (m32 - mQ32)*pow2(m32 - mU32)) + (16*pow4(Xt)*(180*msq2*pow2(mQ32)*
         pow2(mU32)*pow3(m32) - 90*msq2*pow2(m32)*pow2(mU32)*pow3(mQ32) - 270*
         msq2*mU32*pow3(m32)*pow3(mQ32) - 5645*pow2(mU32)*pow3(m32)*pow3(mQ32)
         - 90*msq2*pow2(m32)*pow2(mQ32)*pow3(mU32) - 270*mQ32*msq2*pow3(m32)*
         pow3(mU32) - 5645*pow2(mQ32)*pow3(m32)*pow3(mU32) + 3522*pow2(m32)*
         pow3(mQ32)*pow3(mU32) + 240*msq2*mU32*pow2(mQ32)*pow4(m32) + 240*mQ32*
         msq2*pow2(mU32)*pow4(m32) + 8964*pow2(mQ32)*pow2(mU32)*pow4(m32) +
         2756*mU32*pow3(mQ32)*pow4(m32) + 2756*mQ32*pow3(mU32)*pow4(m32) + 90*
         msq2*mU32*pow2(m32)*pow4(mQ32) + 30*m32*msq2*pow2(mU32)*pow4(mQ32) +
         850*pow2(m32)*pow2(mU32)*pow4(mQ32) - 465*mU32*pow3(m32)*pow4(mQ32) -
         559*m32*pow3(mU32)*pow4(mQ32) - 34*pow4(m32)*pow4(mQ32) + 90*mQ32*msq2
         *pow2(m32)*pow4(mU32) + 30*m32*msq2*pow2(mQ32)*pow4(mU32) + 850*pow2(
         m32)*pow2(mQ32)*pow4(mU32) - 465*mQ32*pow3(m32)*pow4(mU32) - 559*m32*
         pow3(mQ32)*pow4(mU32) - 34*pow4(m32)*pow4(mU32) + 120*pow4(mQ32)*pow4(
         mU32) - 180*mQ32*msq2*mU32*pow5(m32) - 4190*mU32*pow2(mQ32)*pow5(m32)
         - 4190*mQ32*pow2(mU32)*pow5(m32) + 118*pow3(mQ32)*pow5(m32) + 118*pow3
         (mU32)*pow5(m32) + 9*mU32*pow2(m32)*pow5(mQ32) + 3*m32*pow2(mU32)*pow5
         (mQ32) + 9*mQ32*pow2(m32)*pow5(mU32) + 3*m32*pow2(mQ32)*pow5(mU32) +
         1876*mQ32*mU32*pow6(m32) - 150*pow2(mQ32)*pow6(m32) - 150*pow2(mU32)*
         pow6(m32) + 66*mQ32*pow7(m32) + 66*mU32*pow7(m32)))/(3.*mQ32*mU32*pow2
         (mQ32 - mU32)*pow3(m32 - mQ32)*pow3(m32 - mU32)) + (8*(-180*msq2*pow2(
         mQ32)*pow2(mU32)*pow3(m32) + 90*msq2*pow2(m32)*pow2(mU32)*pow3(mQ32) +
         270*msq2*mU32*pow3(m32)*pow3(mQ32) + 3785*pow2(mU32)*pow3(m32)*pow3(
         mQ32) + 90*msq2*pow2(m32)*pow2(mQ32)*pow3(mU32) + 270*mQ32*msq2*pow3(
         m32)*pow3(mU32) + 3785*pow2(mQ32)*pow3(m32)*pow3(mU32) - 2470*pow2(m32
         )*pow3(mQ32)*pow3(mU32) - 240*msq2*mU32*pow2(mQ32)*pow4(m32) - 240*
         mQ32*msq2*pow2(mU32)*pow4(m32) - 5432*pow2(mQ32)*pow2(mU32)*pow4(m32)
         - 2016*mU32*pow3(mQ32)*pow4(m32) - 2016*mQ32*pow3(mU32)*pow4(m32) - 90
         *msq2*mU32*pow2(m32)*pow4(mQ32) - 30*m32*msq2*pow2(mU32)*pow4(mQ32) -
         860*pow2(m32)*pow2(mU32)*pow4(mQ32) + 519*mU32*pow3(m32)*pow4(mQ32) +
         513*m32*pow3(mU32)*pow4(mQ32) - 68*pow4(m32)*pow4(mQ32) - 90*mQ32*msq2
         *pow2(m32)*pow4(mU32) - 30*m32*msq2*pow2(mQ32)*pow4(mU32) - 860*pow2(
         m32)*pow2(mQ32)*pow4(mU32) + 519*mQ32*pow3(m32)*pow4(mU32) + 513*m32*
         pow3(mQ32)*pow4(mU32) - 68*pow4(m32)*pow4(mU32) - 96*pow4(mQ32)*pow4(
         mU32) + 180*mQ32*msq2*mU32*pow5(m32) + 2738*mU32*pow2(mQ32)*pow5(m32)
         + 2738*mQ32*pow2(mU32)*pow5(m32) + 236*pow3(mQ32)*pow5(m32) + 236*pow3
         (mU32)*pow5(m32) - 9*mU32*pow2(m32)*pow5(mQ32) - 3*m32*pow2(mU32)*pow5
         (mQ32) - 9*mQ32*pow2(m32)*pow5(mU32) - 3*m32*pow2(mQ32)*pow5(mU32) -
         1336*mQ32*mU32*pow6(m32) - 300*pow2(mQ32)*pow6(m32) - 300*pow2(mU32)*
         pow6(m32) + 132*mQ32*pow7(m32) + 132*mU32*pow7(m32)))/(3.*mQ32*mU32*
         pow3(m32 - mQ32)*pow3(m32 - mU32)) + log(msq2/MR2)*((96*(-m32 + mQ32 +
         mU32)*pow2(m32)*pow2(Xt))/(mQ32*(-m32 + mQ32)*(m32 - mU32)*mU32) + (
         640*pow5(Xt)*(-11*mQ32*mU32*pow3(m3) + 5*mQ32*pow5(m3) + 5*mU32*pow5(
         m3) + pow7(m3)))/(3.*pow2(m32 - mQ32)*pow2(m32 - mU32)*pow2(mQ32 -
         mU32)) - (32*Xt*(-119*mQ32*mU32*pow3(m3) + 59*mQ32*pow5(m3) + 59*mU32*
         pow5(m3) + pow7(m3)))/(3.*pow2(m32 - mQ32)*pow2(m32 - mU32)) - (16*
         pow4(Xt)*(528*pow2(mU32)*pow3(m32)*pow3(mQ32) + 528*pow2(mQ32)*pow3(
         m32)*pow3(mU32) - 54*pow2(m32)*pow3(mQ32)*pow3(mU32) - 1002*pow2(mQ32)
         *pow2(mU32)*pow4(m32) - 150*mU32*pow3(mQ32)*pow4(m32) - 150*mQ32*pow3(
         mU32)*pow4(m32) - 167*pow2(m32)*pow2(mU32)*pow4(mQ32) + 47*mU32*pow3(
         m32)*pow4(mQ32) + 9*m32*pow3(mU32)*pow4(mQ32) - 9*pow4(m32)*pow4(mQ32)
         - 167*pow2(m32)*pow2(mQ32)*pow4(mU32) + 47*mQ32*pow3(m32)*pow4(mU32)
         + 9*m32*pow3(mQ32)*pow4(mU32) - 9*pow4(m32)*pow4(mU32) + 308*mU32*pow2
         (mQ32)*pow5(m32) + 308*mQ32*pow2(mU32)*pow5(m32) + 27*pow3(mQ32)*pow5(
         m32) + 27*pow3(mU32)*pow5(m32) - 94*mQ32*mU32*pow6(m32) - 27*pow2(mQ32
         )*pow6(m32) - 27*pow2(mU32)*pow6(m32) + 9*mQ32*pow7(m32) + 9*mU32*pow7
         (m32)))/(3.*mQ32*mU32*pow2(mQ32 - mU32)*pow3(m32 - mQ32)*pow3(m32 -
         mU32)) - (16*(-237*pow2(mU32)*pow3(m32)*pow3(mQ32) - 237*pow2(mQ32)*
         pow3(m32)*pow3(mU32) + 474*pow2(mQ32)*pow2(mU32)*pow4(m32) + 30*mU32*
         pow3(mQ32)*pow4(m32) + 30*mQ32*pow3(mU32)*pow4(m32) + 79*pow2(m32)*
         pow2(mU32)*pow4(mQ32) - 10*mU32*pow3(m32)*pow4(mQ32) - 9*pow4(m32)*
         pow4(mQ32) + 79*pow2(m32)*pow2(mQ32)*pow4(mU32) - 10*mQ32*pow3(m32)*
         pow4(mU32) - 9*pow4(m32)*pow4(mU32) - 109*mU32*pow2(mQ32)*pow5(m32) -
         109*mQ32*pow2(mU32)*pow5(m32) + 27*pow3(mQ32)*pow5(m32) + 27*pow3(mU32
         )*pow5(m32) + 20*mQ32*mU32*pow6(m32) - 27*pow2(mQ32)*pow6(m32) - 27*
         pow2(mU32)*pow6(m32) + 9*mQ32*pow7(m32) + 9*mU32*pow7(m32)))/(3.*mQ32*
         mU32*pow3(m32 - mQ32)*pow3(m32 - mU32))) - (128*pow5(Xt)*(30*msq2*mU32
         *pow2(mQ32)*pow3(m3) + 30*mQ32*msq2*pow2(mU32)*pow3(m3) + 101*pow2(
         mQ32)*pow2(mU32)*pow3(m3) + 3*mU32*pow3(m3)*pow3(mQ32) + 3*mQ32*pow3(
         m3)*pow3(mU32) - 60*mQ32*msq2*mU32*pow5(m3) - 115*mU32*pow2(mQ32)*pow5
         (m3) - 115*mQ32*pow2(mU32)*pow5(m3) + 123*mQ32*mU32*pow7(m3) - 8*pow2(
         mQ32)*pow7(m3) - 8*pow2(mU32)*pow7(m3) + 8*mQ32*pow9(m3) + 8*mU32*pow9
         (m3)))/(3.*mQ32*mU32*pow2(m32 - mQ32)*pow2(m32 - mU32)*pow2(mQ32 -
         mU32)) - (64*Xt*(-30*msq2*mU32*pow2(mQ32)*pow3(m3) - 30*mQ32*msq2*pow2
         (mU32)*pow3(m3) - 90*pow2(mQ32)*pow2(mU32)*pow3(m3) - 3*mU32*pow3(m3)*
         pow3(mQ32) - 3*mQ32*pow3(m3)*pow3(mU32) + 60*mQ32*msq2*mU32*pow5(m3) +
         112*mU32*pow2(mQ32)*pow5(m3) + 112*mQ32*pow2(mU32)*pow5(m3) - 128*
         mQ32*mU32*pow7(m3) - 16*pow2(mQ32)*pow7(m3) - 16*pow2(mU32)*pow7(m3) +
         16*mQ32*pow9(m3) + 16*mU32*pow9(m3)))/(3.*mQ32*mU32*pow2(m32 - mQ32)*
         pow2(m32 - mU32)) + dilog(1 - m32/mQ32)*((-8192*pow2(Xt)*pow3(m32))/(
         3.*(m32 - mU32)*(mQ32 - mU32)*pow2(m32 - mQ32)) - (64*pow2(m32)*(2 + (
         8*(pow2(m32)*pow2(mQ32) + pow2(m32)*pow2(mU32) - 2*mQ32*pow3(m32) - 2*
         mU32*pow3(m32) + 2*pow4(m32)))/(3.*pow2(m32 - mQ32)*pow2(m32 - mU32)))
         )/pow2(-m32 + mQ32) + (128*m3*(2*m32 - mQ32 - mU32)*pow3(Xt)*(2 + (8*(
         pow2(m32)*pow2(mQ32) + pow2(m32)*pow2(mU32) - 2*mQ32*pow3(m32) - 2*
         mU32*pow3(m32) + 2*pow4(m32)))/(3.*pow2(m32 - mQ32)*pow2(m32 - mU32)))
         )/pow3(mQ32 - mU32) - (64*(2*m32 - mQ32 - mU32)*(52*mQ32*mU32*pow2(m32
         ) + 6*m32*mU32*pow2(mQ32) - 7*pow2(m32)*pow2(mQ32) + 6*m32*mQ32*pow2(
         mU32) - 7*pow2(m32)*pow2(mU32) - 3*pow2(mQ32)*pow2(mU32) - 50*mQ32*
         pow3(m32) - 50*mU32*pow3(m32) + 53*pow4(m32))*pow4(Xt))/(3.*pow2(m32 -
         mQ32)*pow2(m32 - mU32)*pow3(mQ32 - mU32)) + (1024*(-2*m32 + mQ32 +
         mU32)*pow3(m3)*pow5(Xt))/(3.*(m32 - mQ32)*(m32 - mU32)*pow3(mQ32 -
         mU32)) + (512*Xt*(11*pow11(m3) + 3*pow2(mQ32)*pow2(mU32)*pow3(m3) - 6*
         mU32*pow2(mQ32)*pow5(m3) - 6*mQ32*pow2(mU32)*pow5(m3) + 8*mQ32*mU32*
         pow7(m3) + 7*pow2(mQ32)*pow7(m3) + 11*pow2(mU32)*pow7(m3) - 10*mQ32*
         pow9(m3) - 18*mU32*pow9(m3)))/(3.*(mQ32 - mU32)*pow2(m32 - mU32)*pow3(
         m32 - mQ32))) + dilog(1 - m32/mU32)*((-8192*pow2(Xt)*pow3(m32))/(3.*(
         m32 - mQ32)*(-mQ32 + mU32)*pow2(m32 - mU32)) - (64*pow2(m32)*(2 + (8*(
         pow2(m32)*pow2(mQ32) + pow2(m32)*pow2(mU32) - 2*mQ32*pow3(m32) - 2*
         mU32*pow3(m32) + 2*pow4(m32)))/(3.*pow2(m32 - mQ32)*pow2(m32 - mU32)))
         )/pow2(m32 - mU32) + (128*m3*(-2*m32 + mQ32 + mU32)*pow3(Xt)*(2 + (8*(
         pow2(m32)*pow2(mQ32) + pow2(m32)*pow2(mU32) - 2*mQ32*pow3(m32) - 2*
         mU32*pow3(m32) + 2*pow4(m32)))/(3.*pow2(m32 - mQ32)*pow2(m32 - mU32)))
         )/pow3(mQ32 - mU32) + (64*(2*m32 - mQ32 - mU32)*(52*mQ32*mU32*pow2(m32
         ) + 6*m32*mU32*pow2(mQ32) - 7*pow2(m32)*pow2(mQ32) + 6*m32*mQ32*pow2(
         mU32) - 7*pow2(m32)*pow2(mU32) - 3*pow2(mQ32)*pow2(mU32) - 50*mQ32*
         pow3(m32) - 50*mU32*pow3(m32) + 53*pow4(m32))*pow4(Xt))/(3.*pow2(m32 -
         mQ32)*pow2(m32 - mU32)*pow3(mQ32 - mU32)) - (1024*(-2*m32 + mQ32 +
         mU32)*pow3(m3)*pow5(Xt))/(3.*(m32 - mQ32)*(m32 - mU32)*pow3(mQ32 -
         mU32)) - (512*Xt*(11*pow11(m3) + 3*pow2(mQ32)*pow2(mU32)*pow3(m3) - 6*
         mU32*pow2(mQ32)*pow5(m3) - 6*mQ32*pow2(mU32)*pow5(m3) + 8*mQ32*mU32*
         pow7(m3) + 11*pow2(mQ32)*pow7(m3) + 7*pow2(mU32)*pow7(m3) - 18*mQ32*
         pow9(m3) - 10*mU32*pow9(m3)))/(3.*(mQ32 - mU32)*pow2(m32 - mQ32)*pow3(
         m32 - mU32))) + pow2(log(mQ32/MR2))*((-2560*mQ32*(mQ32 + mU32)*pow2(
         m32)*pow6(Xt))/(3.*(m32 - mU32)*pow2(m32 - mQ32)*pow4(mQ32 - mU32)) +
         (16*(-17*pow2(mU32)*pow3(m32)*pow3(mQ32) - 147*pow2(mQ32)*pow3(m32)*
         pow3(mU32) - 81*pow2(m32)*pow3(mQ32)*pow3(mU32) + 397*pow2(mQ32)*pow2(
         mU32)*pow4(m32) + 159*mU32*pow3(mQ32)*pow4(m32) + 321*mQ32*pow3(mU32)*
         pow4(m32) + 101*pow2(m32)*pow2(mU32)*pow4(mQ32) - 157*mU32*pow3(m32)*
         pow4(mQ32) - 12*m32*pow3(mU32)*pow4(mQ32) + 110*pow4(m32)*pow4(mQ32) -
         pow2(m32)*pow2(mQ32)*pow4(mU32) - 36*mQ32*pow3(m32)*pow4(mU32) + 48*
         m32*pow3(mQ32)*pow4(mU32) - 27*pow4(m32)*pow4(mU32) - 12*pow4(mQ32)*
         pow4(mU32) - 427*mU32*pow2(mQ32)*pow5(m32) - 703*mQ32*pow2(mU32)*pow5(
         m32) - 169*pow3(mQ32)*pow5(m32) + 19*pow3(mU32)*pow5(m32) + 45*mU32*
         pow2(m32)*pow5(mQ32) - 36*m32*pow2(mU32)*pow5(mQ32) - 27*pow3(m32)*
         pow5(mQ32) + 12*pow3(mU32)*pow5(mQ32) + 632*mQ32*mU32*pow6(m32) + 238*
         pow2(mQ32)*pow6(m32) + 90*pow2(mU32)*pow6(m32) - 244*mQ32*pow7(m32) -
         140*mU32*pow7(m32) + 64*pow8(m32)))/(3.*(mQ32 - mU32)*pow3(m32 - mU32)
         *pow4(m32 - mQ32)) + (32*pow2(Xt)*(-811*pow2(mU32)*pow3(m32)*pow3(mQ32
         ) - 189*pow2(mQ32)*pow3(m32)*pow3(mU32) + 61*pow2(m32)*pow3(mQ32)*pow3
         (mU32) + 1199*pow2(mQ32)*pow2(mU32)*pow4(m32) + 1401*mU32*pow3(mQ32)*
         pow4(m32) + 79*mQ32*pow3(mU32)*pow4(m32) + 325*pow2(m32)*pow2(mU32)*
         pow4(mQ32) - 525*mU32*pow3(m32)*pow4(mQ32) - 60*m32*pow3(mU32)*pow4(
         mQ32) + 302*pow4(m32)*pow4(mQ32) - 7*pow2(m32)*pow2(mQ32)*pow4(mU32) -
         24*mQ32*pow3(m32)*pow4(mU32) + 48*m32*pow3(mQ32)*pow4(mU32) - pow4(
         m32)*pow4(mU32) - 12*pow4(mQ32)*pow4(mU32) - 1897*mU32*pow2(mQ32)*pow5
         (m32) - 545*mQ32*pow2(mU32)*pow5(m32) - 759*pow3(mQ32)*pow5(m32) + 69*
         pow3(mU32)*pow5(m32) + 97*mU32*pow2(m32)*pow5(mQ32) - 72*m32*pow2(mU32
         )*pow5(mQ32) - 55*pow3(m32)*pow5(mQ32) + 24*pow3(mU32)*pow5(mQ32) +
         962*mQ32*mU32*pow6(m32) + 954*pow2(mQ32)*pow6(m32) - 72*pow2(mU32)*
         pow6(m32) - 502*mQ32*pow7(m32) - 54*mU32*pow7(m32) + 64*pow8(m32)))/(
         3.*pow2(mQ32 - mU32)*pow3(m32 - mU32)*pow4(m32 - mQ32)) - (16*pow4(Xt)
         *(-2648*pow3(m32)*pow3(mQ32)*pow3(mU32) + 8712*pow2(mU32)*pow3(mQ32)*
         pow4(m32) + 936*pow2(mQ32)*pow3(mU32)*pow4(m32) - 3742*pow2(mU32)*pow3
         (m32)*pow4(mQ32) + 1612*pow2(m32)*pow3(mU32)*pow4(mQ32) + 3031*mU32*
         pow4(m32)*pow4(mQ32) + 479*pow2(mQ32)*pow3(m32)*pow4(mU32) - 234*pow2(
         m32)*pow3(mQ32)*pow4(mU32) - 190*mQ32*pow4(m32)*pow4(mU32) - 90*m32*
         pow4(mQ32)*pow4(mU32) - 6232*pow2(mQ32)*pow2(mU32)*pow5(m32) - 8472*
         mU32*pow3(mQ32)*pow5(m32) + 452*mQ32*pow3(mU32)*pow5(m32) - 835*pow4(
         mQ32)*pow5(m32) - pow4(mU32)*pow5(m32) + 238*pow2(m32)*pow2(mU32)*pow5
         (mQ32) - 116*mU32*pow3(m32)*pow5(mQ32) - 180*m32*pow3(mU32)*pow5(mQ32)
         + 58*pow4(m32)*pow5(mQ32) + 24*pow4(mU32)*pow5(mQ32) + 37*pow2(m32)*
         pow2(mQ32)*pow5(mU32) - 64*mQ32*pow3(m32)*pow5(mU32) + 24*m32*pow3(
         mQ32)*pow5(mU32) + 13*pow4(m32)*pow5(mU32) - 6*pow4(mQ32)*pow5(mU32) +
         7366*mU32*pow2(mQ32)*pow6(m32) + 786*mQ32*pow2(mU32)*pow6(m32) + 2618
         *pow3(mQ32)*pow6(m32) - 194*pow3(mU32)*pow6(m32) + 123*mU32*pow2(m32)*
         pow6(mQ32) - 90*m32*pow2(mU32)*pow6(mQ32) - 69*pow3(m32)*pow6(mQ32) +
         30*pow3(mU32)*pow6(mQ32) - 1796*mQ32*mU32*pow7(m32) - 2556*pow2(mQ32)*
         pow7(m32) + 336*pow2(mU32)*pow7(m32) + 788*mQ32*pow8(m32) - 148*mU32*
         pow8(m32)))/(3.*pow3(m32 - mU32)*pow4(m32 - mQ32)*pow4(mQ32 - mU32)) +
         (32*Xt*(16*pow11(m3) - 89*pow2(mQ32)*pow2(mU32)*pow3(m3) - 12*m3*pow2
         (mU32)*pow3(mQ32) - 11*mU32*pow3(m3)*pow3(mQ32) + 267*mU32*pow2(mQ32)*
         pow5(m3) + 172*mQ32*pow2(mU32)*pow5(m3) + 13*pow3(mQ32)*pow5(m3) - 433
         *mQ32*mU32*pow7(m3) - 132*pow2(mQ32)*pow7(m3) + 25*pow2(mU32)*pow7(m3)
         + 199*mQ32*pow9(m3) - 15*mU32*pow9(m3)))/(3.*(mQ32 - mU32)*pow2(m32 -
         mU32)*pow3(m32 - mQ32)) + (32*pow5(Xt)*(36*mQ32*pow11(m3) + 100*mU32*
         pow11(m3) + 152*pow2(mU32)*pow3(m3)*pow3(mQ32) + 93*pow2(mQ32)*pow3(m3
         )*pow3(mU32) + 12*m3*pow3(mQ32)*pow3(mU32) + 12*m3*pow2(mU32)*pow4(
         mQ32) + 59*mU32*pow3(m3)*pow4(mQ32) - 535*pow2(mQ32)*pow2(mU32)*pow5(
         m3) - 384*mU32*pow3(mQ32)*pow5(m3) - 148*mQ32*pow3(mU32)*pow5(m3) - 61
         *pow4(mQ32)*pow5(m3) + 585*mU32*pow2(mQ32)*pow7(m3) + 508*mQ32*pow2(
         mU32)*pow7(m3) + 216*pow3(mQ32)*pow7(m3) + 75*pow3(mU32)*pow7(m3) -
         392*mQ32*mU32*pow9(m3) - 159*pow2(mQ32)*pow9(m3) - 169*pow2(mU32)*pow9
         (m3)))/(3.*pow2(m32 - mU32)*pow3(m32 - mQ32)*pow4(mQ32 - mU32)) + (64*
         pow3(Xt)*(-114*mQ32*pow11(m3) - 2*mU32*pow11(m3) + 40*pow13(m3) - 40*
         pow2(mU32)*pow3(m3)*pow3(mQ32) + 59*pow2(mQ32)*pow3(m3)*pow3(mU32) +
         18*m3*pow3(mQ32)*pow3(mU32) - 6*m3*pow2(mU32)*pow4(mQ32) - 55*mU32*
         pow3(m3)*pow4(mQ32) - 201*pow2(mQ32)*pow2(mU32)*pow5(m3) + 264*mU32*
         pow3(mQ32)*pow5(m3) - 90*mQ32*pow3(mU32)*pow5(m3) + 59*pow4(mQ32)*pow5
         (m3) - 55*mU32*pow2(mQ32)*pow7(m3) + 220*mQ32*pow2(mU32)*pow7(m3) -
         250*pow3(mQ32)*pow7(m3) + 45*pow3(mU32)*pow7(m3) - 56*mQ32*mU32*pow9(
         m3) + 233*pow2(mQ32)*pow9(m3) - 69*pow2(mU32)*pow9(m3)))/(3.*pow2(m32
         - mU32)*pow3(m32 - mQ32)*pow3(mQ32 - mU32))) + pow2(log(mU32/MR2))*((
         -2560*mU32*(mQ32 + mU32)*pow2(m32)*pow6(Xt))/(3.*(m32 - mQ32)*pow2(m32
         - mU32)*pow4(mQ32 - mU32)) + (32*pow2(Xt)*(-189*pow2(mU32)*pow3(m32)*
         pow3(mQ32) - 811*pow2(mQ32)*pow3(m32)*pow3(mU32) + 61*pow2(m32)*pow3(
         mQ32)*pow3(mU32) + 1199*pow2(mQ32)*pow2(mU32)*pow4(m32) + 79*mU32*pow3
         (mQ32)*pow4(m32) + 1401*mQ32*pow3(mU32)*pow4(m32) - 7*pow2(m32)*pow2(
         mU32)*pow4(mQ32) - 24*mU32*pow3(m32)*pow4(mQ32) + 48*m32*pow3(mU32)*
         pow4(mQ32) - pow4(m32)*pow4(mQ32) + 325*pow2(m32)*pow2(mQ32)*pow4(mU32
         ) - 525*mQ32*pow3(m32)*pow4(mU32) - 60*m32*pow3(mQ32)*pow4(mU32) + 302
         *pow4(m32)*pow4(mU32) - 12*pow4(mQ32)*pow4(mU32) - 545*mU32*pow2(mQ32)
         *pow5(m32) - 1897*mQ32*pow2(mU32)*pow5(m32) + 69*pow3(mQ32)*pow5(m32)
         - 759*pow3(mU32)*pow5(m32) + 97*mQ32*pow2(m32)*pow5(mU32) - 72*m32*
         pow2(mQ32)*pow5(mU32) - 55*pow3(m32)*pow5(mU32) + 24*pow3(mQ32)*pow5(
         mU32) + 962*mQ32*mU32*pow6(m32) - 72*pow2(mQ32)*pow6(m32) + 954*pow2(
         mU32)*pow6(m32) - 54*mQ32*pow7(m32) - 502*mU32*pow7(m32) + 64*pow8(m32
         )))/(3.*pow2(mQ32 - mU32)*pow3(m32 - mQ32)*pow4(m32 - mU32)) - (16*(
         -146*pow2(mU32)*pow3(m32)*pow3(mQ32) - 20*pow2(mQ32)*pow3(m32)*pow3(
         mU32) - 80*pow2(m32)*pow3(mQ32)*pow3(mU32) + 400*pow2(mQ32)*pow2(mU32)
         *pow4(m32) + 316*mU32*pow3(mQ32)*pow4(m32) + 162*mQ32*pow3(mU32)*pow4(
         m32) - 2*pow2(m32)*pow2(mU32)*pow4(mQ32) - 34*mU32*pow3(m32)*pow4(mQ32
         ) + 48*m32*pow3(mU32)*pow4(mQ32) - 28*pow4(m32)*pow4(mQ32) + 101*pow2(
         m32)*pow2(mQ32)*pow4(mU32) - 157*mQ32*pow3(m32)*pow4(mU32) - 12*m32*
         pow3(mQ32)*pow4(mU32) + 110*pow4(m32)*pow4(mU32) - 12*pow4(mQ32)*pow4(
         mU32) - 700*mU32*pow2(mQ32)*pow5(m32) - 432*mQ32*pow2(mU32)*pow5(m32)
         + 22*pow3(mQ32)*pow5(m32) - 170*pow3(mU32)*pow5(m32) + 45*mQ32*pow2(
         m32)*pow5(mU32) - 36*m32*pow2(mQ32)*pow5(mU32) - 27*pow3(m32)*pow5(
         mU32) + 12*pow3(mQ32)*pow5(mU32) + 633*mQ32*mU32*pow6(m32) + 87*pow2(
         mQ32)*pow6(m32) + 240*pow2(mU32)*pow6(m32) - 139*mQ32*pow7(m32) - 245*
         mU32*pow7(m32) + 64*pow8(m32)))/(3.*(mQ32 - mU32)*pow3(m32 - mQ32)*
         pow4(m32 - mU32)) + (16*pow4(Xt)*(2648*pow3(m32)*pow3(mQ32)*pow3(mU32)
         - 936*pow2(mU32)*pow3(mQ32)*pow4(m32) - 8712*pow2(mQ32)*pow3(mU32)*
         pow4(m32) - 479*pow2(mU32)*pow3(m32)*pow4(mQ32) + 234*pow2(m32)*pow3(
         mU32)*pow4(mQ32) + 190*mU32*pow4(m32)*pow4(mQ32) + 3742*pow2(mQ32)*
         pow3(m32)*pow4(mU32) - 1612*pow2(m32)*pow3(mQ32)*pow4(mU32) - 3031*
         mQ32*pow4(m32)*pow4(mU32) + 90*m32*pow4(mQ32)*pow4(mU32) + 6232*pow2(
         mQ32)*pow2(mU32)*pow5(m32) - 452*mU32*pow3(mQ32)*pow5(m32) + 8472*mQ32
         *pow3(mU32)*pow5(m32) + pow4(mQ32)*pow5(m32) + 835*pow4(mU32)*pow5(m32
         ) - 37*pow2(m32)*pow2(mU32)*pow5(mQ32) + 64*mU32*pow3(m32)*pow5(mQ32)
         - 24*m32*pow3(mU32)*pow5(mQ32) - 13*pow4(m32)*pow5(mQ32) + 6*pow4(mU32
         )*pow5(mQ32) - 238*pow2(m32)*pow2(mQ32)*pow5(mU32) + 116*mQ32*pow3(m32
         )*pow5(mU32) + 180*m32*pow3(mQ32)*pow5(mU32) - 58*pow4(m32)*pow5(mU32)
         - 24*pow4(mQ32)*pow5(mU32) - 786*mU32*pow2(mQ32)*pow6(m32) - 7366*
         mQ32*pow2(mU32)*pow6(m32) + 194*pow3(mQ32)*pow6(m32) - 2618*pow3(mU32)
         *pow6(m32) - 123*mQ32*pow2(m32)*pow6(mU32) + 90*m32*pow2(mQ32)*pow6(
         mU32) + 69*pow3(m32)*pow6(mU32) - 30*pow3(mQ32)*pow6(mU32) + 1796*mQ32
         *mU32*pow7(m32) - 336*pow2(mQ32)*pow7(m32) + 2556*pow2(mU32)*pow7(m32)
         + 148*mQ32*pow8(m32) - 788*mU32*pow8(m32)))/(3.*pow3(m32 - mQ32)*pow4
         (m32 - mU32)*pow4(mQ32 - mU32)) - (32*Xt*(18*pow11(m3) - 87*pow2(mQ32)
         *pow2(mU32)*pow3(m3) - 12*m3*pow2(mQ32)*pow3(mU32) - 11*mQ32*pow3(m3)*
         pow3(mU32) + 168*mU32*pow2(mQ32)*pow5(m3) + 263*mQ32*pow2(mU32)*pow5(
         m3) + 13*pow3(mU32)*pow5(m3) - 425*mQ32*mU32*pow7(m3) + 27*pow2(mQ32)*
         pow7(m3) - 130*pow2(mU32)*pow7(m3) - 19*mQ32*pow9(m3) + 195*mU32*pow9(
         m3)))/(3.*(mQ32 - mU32)*pow2(m32 - mQ32)*pow3(m32 - mU32)) + (32*pow5(
         Xt)*(100*mQ32*pow11(m3) + 36*mU32*pow11(m3) + 93*pow2(mU32)*pow3(m3)*
         pow3(mQ32) + 152*pow2(mQ32)*pow3(m3)*pow3(mU32) + 12*m3*pow3(mQ32)*
         pow3(mU32) + 12*m3*pow2(mQ32)*pow4(mU32) + 59*mQ32*pow3(m3)*pow4(mU32)
         - 535*pow2(mQ32)*pow2(mU32)*pow5(m3) - 148*mU32*pow3(mQ32)*pow5(m3) -
         384*mQ32*pow3(mU32)*pow5(m3) - 61*pow4(mU32)*pow5(m3) + 508*mU32*pow2
         (mQ32)*pow7(m3) + 585*mQ32*pow2(mU32)*pow7(m3) + 75*pow3(mQ32)*pow7(m3
         ) + 216*pow3(mU32)*pow7(m3) - 392*mQ32*mU32*pow9(m3) - 169*pow2(mQ32)*
         pow9(m3) - 159*pow2(mU32)*pow9(m3)))/(3.*pow2(m32 - mQ32)*pow3(m32 -
         mU32)*pow4(mQ32 - mU32)) - (64*pow3(Xt)*(-6*mQ32*pow11(m3) - 120*mU32*
         pow11(m3) + 42*pow13(m3) + 59*pow2(mU32)*pow3(m3)*pow3(mQ32) - 42*pow2
         (mQ32)*pow3(m3)*pow3(mU32) + 18*m3*pow3(mQ32)*pow3(mU32) - 6*m3*pow2(
         mQ32)*pow4(mU32) - 55*mQ32*pow3(m3)*pow4(mU32) - 195*pow2(mQ32)*pow2(
         mU32)*pow5(m3) - 90*mU32*pow3(mQ32)*pow5(m3) + 268*mQ32*pow3(mU32)*
         pow5(m3) + 59*pow4(mU32)*pow5(m3) + 214*mU32*pow2(mQ32)*pow7(m3) - 67*
         mQ32*pow2(mU32)*pow7(m3) + 45*pow3(mQ32)*pow7(m3) - 252*pow3(mU32)*
         pow7(m3) - 44*mQ32*mU32*pow9(m3) - 67*pow2(mQ32)*pow9(m3) + 239*pow2(
         mU32)*pow9(m3)))/(3.*pow2(m32 - mQ32)*pow3(m32 - mU32)*pow3(mQ32 -
         mU32))) + log(mU32/MR2)*((-5120*mU32*pow2(m32)*pow6(Xt))/(3.*(m32 -
         mQ32)*pow2(m32 - mU32)*pow3(mQ32 - mU32)) + (16*pow2(Xt)*(-180*msq2*
         pow2(mQ32)*pow2(mU32)*pow3(m32) + 90*msq2*pow2(m32)*pow2(mU32)*pow3(
         mQ32) + 270*msq2*mU32*pow3(m32)*pow3(mQ32) + 3987*pow2(mU32)*pow3(m32)
         *pow3(mQ32) + 90*msq2*pow2(m32)*pow2(mQ32)*pow3(mU32) + 270*mQ32*msq2*
         pow3(m32)*pow3(mU32) + 3519*pow2(mQ32)*pow3(m32)*pow3(mU32) - 2430*
         pow2(m32)*pow3(mQ32)*pow3(mU32) - 240*msq2*mU32*pow2(mQ32)*pow4(m32) -
         240*mQ32*msq2*pow2(mU32)*pow4(m32) - 5696*pow2(mQ32)*pow2(mU32)*pow4(
         m32) - 2436*mU32*pow3(mQ32)*pow4(m32) - 1692*mQ32*pow3(mU32)*pow4(m32)
         - 90*msq2*mU32*pow2(m32)*pow4(mQ32) - 30*m32*msq2*pow2(mU32)*pow4(
         mQ32) - 850*pow2(m32)*pow2(mU32)*pow4(mQ32) + 587*mU32*pow3(m32)*pow4(
         mQ32) + 513*m32*pow3(mU32)*pow4(mQ32) - 2*pow4(m32)*pow4(mQ32) - 90*
         mQ32*msq2*pow2(m32)*pow4(mU32) - 30*m32*msq2*pow2(mQ32)*pow4(mU32) -
         734*pow2(m32)*pow2(mQ32)*pow4(mU32) + 395*mQ32*pow3(m32)*pow4(mU32) +
         489*m32*pow3(mQ32)*pow4(mU32) + 34*pow4(m32)*pow4(mU32) - 96*pow4(mQ32
         )*pow4(mU32) + 180*mQ32*msq2*mU32*pow5(m32) + 3280*mU32*pow2(mQ32)*
         pow5(m32) + 2660*mQ32*pow2(mU32)*pow5(m32) + 6*pow3(mQ32)*pow5(m32) -
         102*pow3(mU32)*pow5(m32) - 9*mU32*pow2(m32)*pow5(mQ32) - 3*m32*pow2(
         mU32)*pow5(mQ32) - 9*mQ32*pow2(m32)*pow5(mU32) - 3*m32*pow2(mQ32)*pow5
         (mU32) - 1472*mQ32*mU32*pow6(m32) - 6*pow2(mQ32)*pow6(m32) + 70*pow2(
         mU32)*pow6(m32) + 2*mQ32*pow7(m32) - 2*mU32*pow7(m32)))/(3.*mQ32*(mQ32
         - mU32)*mU32*pow3(m32 - mQ32)*pow3(m32 - mU32)) - (8*(450*msq2*pow3(
         m32)*pow3(mQ32)*pow3(mU32) - 210*msq2*pow2(mU32)*pow3(mQ32)*pow4(m32)
         + 450*msq2*pow2(mQ32)*pow3(mU32)*pow4(m32) + 1602*pow3(mQ32)*pow3(mU32
         )*pow4(m32) - 90*msq2*pow2(mU32)*pow3(m32)*pow4(mQ32) - 150*msq2*pow2(
         m32)*pow3(mU32)*pow4(mQ32) - 3698*pow3(m32)*pow3(mU32)*pow4(mQ32) +
         270*msq2*mU32*pow4(m32)*pow4(mQ32) + 3528*pow2(mU32)*pow4(m32)*pow4(
         mQ32) - 630*msq2*pow2(mQ32)*pow3(m32)*pow4(mU32) - 30*msq2*pow2(m32)*
         pow3(mQ32)*pow4(mU32) + 1682*pow3(m32)*pow3(mQ32)*pow4(mU32) - 510*
         mQ32*msq2*pow4(m32)*pow4(mU32) - 5506*pow2(mQ32)*pow4(m32)*pow4(mU32)
         - 30*m32*msq2*pow4(mQ32)*pow4(mU32) + 1522*pow2(m32)*pow4(mQ32)*pow4(
         mU32) - 180*msq2*pow2(mQ32)*pow2(mU32)*pow5(m32) - 240*msq2*mU32*pow3(
         mQ32)*pow5(m32) - 3080*pow2(mU32)*pow3(mQ32)*pow5(m32) + 420*mQ32*msq2
         *pow3(mU32)*pow5(m32) + 3634*pow2(mQ32)*pow3(mU32)*pow5(m32) - 1902*
         mU32*pow4(mQ32)*pow5(m32) + 3774*mQ32*pow4(mU32)*pow5(m32) + 60*msq2*
         pow2(m32)*pow2(mU32)*pow5(mQ32) - 90*msq2*mU32*pow3(m32)*pow5(mQ32) -
         972*pow2(mU32)*pow3(m32)*pow5(mQ32) + 30*m32*msq2*pow3(mU32)*pow5(mQ32
         ) + 1019*pow2(m32)*pow3(mU32)*pow5(mQ32) + 547*mU32*pow4(m32)*pow5(
         mQ32) - 540*m32*pow4(mU32)*pow5(mQ32) - 2*pow5(m32)*pow5(mQ32) + 210*
         msq2*pow2(m32)*pow2(mQ32)*pow5(mU32) + 360*mQ32*msq2*pow3(m32)*pow5(
         mU32) + 3339*pow2(mQ32)*pow3(m32)*pow5(mU32) + 30*m32*msq2*pow3(mQ32)*
         pow5(mU32) - 2026*pow2(m32)*pow3(mQ32)*pow5(mU32) - 2057*mQ32*pow4(m32
         )*pow5(mU32) + 84*m32*pow4(mQ32)*pow5(mU32) + 136*pow5(m32)*pow5(mU32)
         + 84*pow5(mQ32)*pow5(mU32) + 180*msq2*mU32*pow2(mQ32)*pow6(m32) - 180
         *mQ32*msq2*pow2(mU32)*pow6(m32) - 610*pow2(mQ32)*pow2(mU32)*pow6(m32)
         + 2188*mU32*pow3(mQ32)*pow6(m32) - 3332*mQ32*pow3(mU32)*pow6(m32) + 6*
         pow4(mQ32)*pow6(m32) - 172*pow4(mU32)*pow6(m32) + 6*pow2(m32)*pow2(
         mU32)*pow6(mQ32) - 9*mU32*pow3(m32)*pow6(mQ32) + 3*m32*pow3(mU32)*pow6
         (mQ32) - 90*mQ32*msq2*pow2(m32)*pow6(mU32) - 30*m32*msq2*pow2(mQ32)*
         pow6(mU32) - 640*pow2(m32)*pow2(mQ32)*pow6(mU32) + 426*mQ32*pow3(m32)*
         pow6(mU32) + 456*m32*pow3(mQ32)*pow6(mU32) - 34*pow4(m32)*pow6(mU32) -
         84*pow4(mQ32)*pow6(mU32) - 760*mU32*pow2(mQ32)*pow7(m32) + 1462*mQ32*
         pow2(mU32)*pow7(m32) - 6*pow3(mQ32)*pow7(m32) + 72*pow3(mU32)*pow7(m32
         ) - 9*mQ32*pow2(m32)*pow7(mU32) - 3*m32*pow2(mQ32)*pow7(mU32) - 128*
         mQ32*mU32*pow8(m32) + 2*pow2(mQ32)*pow8(m32) - 2*pow2(mU32)*pow8(m32))
         )/(3.*mQ32*(mQ32 - mU32)*mU32*pow3(m32 - mQ32)*pow4(m32 - mU32)) - (8*
         pow4(Xt)*(90*msq2*pow3(m32)*pow3(mQ32)*pow3(mU32) + 330*msq2*pow2(mU32
         )*pow3(mQ32)*pow4(m32) + 570*msq2*pow2(mQ32)*pow3(mU32)*pow4(m32) +
         25638*pow3(mQ32)*pow3(mU32)*pow4(m32) - 270*msq2*pow2(mU32)*pow3(m32)*
         pow4(mQ32) - 30*msq2*pow2(m32)*pow3(mU32)*pow4(mQ32) - 10862*pow3(m32)
         *pow3(mU32)*pow4(mQ32) + 270*msq2*mU32*pow4(m32)*pow4(mQ32) + 9676*
         pow2(mU32)*pow4(m32)*pow4(mQ32) - 90*msq2*pow2(mQ32)*pow3(m32)*pow4(
         mU32) - 210*msq2*pow2(m32)*pow3(mQ32)*pow4(mU32) - 18070*pow3(m32)*
         pow3(mQ32)*pow4(mU32) + 510*mQ32*msq2*pow4(m32)*pow4(mU32) + 22990*
         pow2(mQ32)*pow4(m32)*pow4(mU32) + 30*m32*msq2*pow4(mQ32)*pow4(mU32) +
         5474*pow2(m32)*pow4(mQ32)*pow4(mU32) - 660*msq2*pow2(mQ32)*pow2(mU32)*
         pow5(m32) - 240*msq2*mU32*pow3(mQ32)*pow5(m32) - 17588*pow2(mU32)*pow3
         (mQ32)*pow5(m32) - 420*mQ32*msq2*pow3(mU32)*pow5(m32) - 25584*pow2(
         mQ32)*pow3(mU32)*pow5(m32) - 2996*mU32*pow4(mQ32)*pow5(m32) - 9926*
         mQ32*pow4(mU32)*pow5(m32) + 60*msq2*pow2(m32)*pow2(mU32)*pow5(mQ32) -
         90*msq2*mU32*pow3(m32)*pow5(mQ32) - 1828*pow2(mU32)*pow3(m32)*pow5(
         mQ32) + 30*m32*msq2*pow3(mU32)*pow5(mQ32) + 1723*pow2(m32)*pow3(mU32)*
         pow5(mQ32) + 643*mU32*pow4(m32)*pow5(mQ32) - 616*m32*pow4(mU32)*pow5(
         mQ32) - 2*pow5(m32)*pow5(mQ32) - 30*msq2*pow2(m32)*pow2(mQ32)*pow5(
         mU32) - 360*mQ32*msq2*pow3(m32)*pow5(mU32) - 9553*pow2(mQ32)*pow3(m32)
         *pow5(mU32) + 30*m32*msq2*pow3(mQ32)*pow5(mU32) + 6600*pow2(m32)*pow3(
         mQ32)*pow5(mU32) + 4623*mQ32*pow4(m32)*pow5(mU32) - 1602*m32*pow4(mQ32
         )*pow5(mU32) + 136*pow5(m32)*pow5(mU32) + 96*pow5(mQ32)*pow5(mU32) +
         180*msq2*mU32*pow2(mQ32)*pow6(m32) + 180*mQ32*msq2*pow2(mU32)*pow6(m32
         ) + 12962*pow2(mQ32)*pow2(mU32)*pow6(m32) + 4406*mU32*pow3(mQ32)*pow6(
         m32) + 9110*mQ32*pow3(mU32)*pow6(m32) + 6*pow4(mQ32)*pow6(m32) - 172*
         pow4(mU32)*pow6(m32) + 6*pow2(m32)*pow2(mU32)*pow6(mQ32) - 9*mU32*pow3
         (m32)*pow6(mQ32) + 3*m32*pow3(mU32)*pow6(mQ32) + 90*mQ32*msq2*pow2(m32
         )*pow6(mU32) + 30*m32*msq2*pow2(mQ32)*pow6(mU32) + 1268*pow2(m32)*pow2
         (mQ32)*pow6(mU32) - 702*mQ32*pow3(m32)*pow6(mU32) - 884*m32*pow3(mQ32)
         *pow6(mU32) - 34*pow4(m32)*pow6(mU32) + 240*pow4(mQ32)*pow6(mU32) -
         2152*mU32*pow2(mQ32)*pow7(m32) - 3178*mQ32*pow2(mU32)*pow7(m32) - 6*
         pow3(mQ32)*pow7(m32) + 72*pow3(mU32)*pow7(m32) + 9*mQ32*pow2(m32)*pow7
         (mU32) + 3*m32*pow2(mQ32)*pow7(mU32) + 80*mQ32*mU32*pow8(m32) + 2*pow2
         (mQ32)*pow8(m32) - 2*pow2(mU32)*pow8(m32)))/(3.*mQ32*mU32*pow3(m32 -
         mQ32)*pow3(mQ32 - mU32)*pow4(m32 - mU32)) + (64*pow5(Xt)*(52*mQ32*
         pow11(m3) + 16*mU32*pow11(m3) + 60*msq2*pow2(mQ32)*pow2(mU32)*pow3(m3)
         + 30*msq2*mU32*pow3(m3)*pow3(mQ32) + 213*pow2(mU32)*pow3(m3)*pow3(
         mQ32) + 30*mQ32*msq2*pow3(m3)*pow3(mU32) + 179*pow2(mQ32)*pow3(m3)*
         pow3(mU32) + 12*m3*pow3(mQ32)*pow3(mU32) + 3*mU32*pow3(m3)*pow4(mQ32)
         + 3*mQ32*pow3(m3)*pow4(mU32) - 120*msq2*mU32*pow2(mQ32)*pow5(m3) - 90*
         mQ32*msq2*pow2(mU32)*pow5(m3) - 737*pow2(mQ32)*pow2(mU32)*pow5(m3) -
         30*msq2*pow3(mQ32)*pow5(m3) - 407*mU32*pow3(mQ32)*pow5(m3) - 219*mQ32*
         pow3(mU32)*pow5(m3) - 3*pow4(mQ32)*pow5(m3) + 60*mQ32*msq2*mU32*pow7(
         m3) + 60*msq2*pow2(mQ32)*pow7(m3) + 770*mU32*pow2(mQ32)*pow7(m3) + 558
         *mQ32*pow2(mU32)*pow7(m3) + 214*pow3(mQ32)*pow7(m3) + 16*pow3(mU32)*
         pow7(m3) - 362*mQ32*mU32*pow9(m3) - 276*pow2(mQ32)*pow9(m3) - 32*pow2(
         mU32)*pow9(m3)))/(3.*mQ32*pow2(m32 - mQ32)*pow3(m32 - mU32)*pow3(mQ32
         - mU32)) + (64*Xt*(4*mQ32*pow11(m3) + 16*mU32*pow11(m3) + 30*msq2*mU32
         *pow3(m3)*pow3(mQ32) + 37*pow2(mU32)*pow3(m3)*pow3(mQ32) - 30*mQ32*
         msq2*pow3(m3)*pow3(mU32) - 35*pow2(mQ32)*pow3(m3)*pow3(mU32) - 12*m3*
         pow3(mQ32)*pow3(mU32) + 3*mU32*pow3(m3)*pow4(mQ32) - 3*mQ32*pow3(m3)*
         pow4(mU32) - 60*msq2*mU32*pow2(mQ32)*pow5(m3) + 90*mQ32*msq2*pow2(mU32
         )*pow5(m3) + 111*pow2(mQ32)*pow2(mU32)*pow5(m3) - 30*msq2*pow3(mQ32)*
         pow5(m3) - 68*mU32*pow3(mQ32)*pow5(m3) + 86*mQ32*pow3(mU32)*pow5(m3) -
         3*pow4(mQ32)*pow5(m3) - 60*mQ32*msq2*mU32*pow7(m3) + 60*msq2*pow2(
         mQ32)*pow7(m3) - 106*mU32*pow2(mQ32)*pow7(m3) - 223*mQ32*pow2(mU32)*
         pow7(m3) + 123*pow3(mQ32)*pow7(m3) + 16*pow3(mU32)*pow7(m3) + 216*mQ32
         *mU32*pow9(m3) - 130*pow2(mQ32)*pow9(m3) - 32*pow2(mU32)*pow9(m3)))/(
         3.*mQ32*(mQ32 - mU32)*pow2(m32 - mQ32)*pow3(m32 - mU32)) - (128*pow3(
         Xt)*(22*mQ32*pow11(m3) + 30*msq2*pow2(mQ32)*pow2(mU32)*pow3(m3) + 30*
         msq2*mU32*pow3(m3)*pow3(mQ32) + 106*pow2(mU32)*pow3(m3)*pow3(mQ32) -
         30*mQ32*msq2*pow3(m3)*pow3(mU32) - 115*pow2(mQ32)*pow3(m3)*pow3(mU32)
         + 18*m3*pow3(mQ32)*pow3(mU32) + 6*m3*pow2(mU32)*pow4(mQ32) - 30*msq2*
         pow3(m3)*pow4(mQ32) - 91*mU32*pow3(m3)*pow4(mQ32) - 3*mQ32*pow3(m3)*
         pow4(mU32) - 120*msq2*mU32*pow2(mQ32)*pow5(m3) + 60*mQ32*msq2*pow2(
         mU32)*pow5(m3) - 17*pow2(mQ32)*pow2(mU32)*pow5(m3) + 60*msq2*pow3(mQ32
         )*pow5(m3) - 85*mU32*pow3(mQ32)*pow5(m3) + 165*mQ32*pow3(mU32)*pow5(m3
         ) + 153*pow4(mQ32)*pow5(m3) - 3*pow3(m3)*pow5(mQ32) + 238*mU32*pow2(
         mQ32)*pow7(m3) - 217*mQ32*pow2(mU32)*pow7(m3) - 193*pow3(mQ32)*pow7(m3
         ) - 16*pow3(mU32)*pow7(m3) + 6*mQ32*mU32*pow9(m3) + 10*pow2(mQ32)*pow9
         (m3) + 16*pow2(mU32)*pow9(m3)))/(3.*mQ32*pow2(m32 - mQ32)*pow2(m32 -
         mU32)*pow3(mQ32 - mU32)) + log(msq2/MR2)*((-320*pow2(Xt)*(-21*mU32*
         pow2(mQ32)*pow3(m32) - 21*mQ32*pow2(mU32)*pow3(m32) + 7*mU32*pow2(m32)
         *pow3(mQ32) - pow3(m32)*pow3(mQ32) + 7*mQ32*pow2(m32)*pow3(mU32) -
         pow3(m32)*pow3(mU32) + 42*mQ32*mU32*pow4(m32) + 3*pow2(mQ32)*pow4(m32)
         + 3*pow2(mU32)*pow4(m32) - 10*mQ32*pow5(m32) - 10*mU32*pow5(m32) + 2*
         pow6(m32)))/(3.*(mQ32 - mU32)*pow3(m32 - mQ32)*pow3(m32 - mU32)) + (
         160*(mQ32 + mU32)*pow4(Xt)*(-21*mU32*pow2(mQ32)*pow3(m32) - 21*mQ32*
         pow2(mU32)*pow3(m32) + 7*mU32*pow2(m32)*pow3(mQ32) - pow3(m32)*pow3(
         mQ32) + 7*mQ32*pow2(m32)*pow3(mU32) - pow3(m32)*pow3(mU32) + 42*mQ32*
         mU32*pow4(m32) + 3*pow2(mQ32)*pow4(m32) + 3*pow2(mU32)*pow4(m32) - 10*
         mQ32*pow5(m32) - 10*mU32*pow5(m32) + 2*pow6(m32)))/(3.*pow3(m32 - mQ32
         )*pow3(m32 - mU32)*pow3(mQ32 - mU32)) + (16*(-237*mU32*pow2(mQ32)*pow3
         (m32) - 210*mQ32*pow2(mU32)*pow3(m32) + 79*mU32*pow2(m32)*pow3(mQ32) -
         19*pow3(m32)*pow3(mQ32) + 70*mQ32*pow2(m32)*pow3(mU32) - 10*pow3(m32)
         *pow3(mU32) + 447*mQ32*mU32*pow4(m32) + 57*pow2(mQ32)*pow4(m32) + 30*
         pow2(mU32)*pow4(m32) - 127*mQ32*pow5(m32) - 109*mU32*pow5(m32) + 29*
         pow6(m32)))/(3.*pow3(m32 - mQ32)*pow3(m32 - mU32)) + (320*(mQ32 + mU32
         )*pow5(Xt)*(-11*mQ32*mU32*pow3(m3) + 5*mQ32*pow5(m3) + 5*mU32*pow5(m3)
         + pow7(m3)))/(3.*pow2(m32 - mQ32)*pow2(m32 - mU32)*pow3(mQ32 - mU32))
         + pow3(Xt)*((192*(4*mQ32*mU32 - 2*m32*(mQ32 + mU32) + 2*pow2(m32) -
         pow2(mQ32) - pow2(mU32))*pow3(m3))/((m32 - mQ32)*(m32 - mU32)*pow3(
         mQ32 - mU32)) - (64*(-119*mQ32*mU32*pow3(m3) + 59*mQ32*pow5(m3) + 59*
         mU32*pow5(m3) + pow7(m3)))/(3.*(mQ32 - mU32)*pow2(m32 - mQ32)*pow2(m32
         - mU32))) + (64*Xt*(-64*mU32*pow2(mQ32)*pow3(m3) + 55*mQ32*pow2(mU32)
         *pow3(m3) + 18*mQ32*mU32*pow5(m3) + 34*pow2(mQ32)*pow5(m3) - 25*pow2(
         mU32)*pow5(m3) - 13*mQ32*pow7(m3) - 14*mU32*pow7(m3) + 9*pow9(m3)))/(
         3.*(mQ32 - mU32)*pow2(m32 - mQ32)*pow2(m32 - mU32)))) + log(mQ32/MR2)*
         ((5120*mQ32*pow2(m32)*pow6(Xt))/(3.*(m32 - mU32)*pow2(m32 - mQ32)*pow3
         (mQ32 - mU32)) + log(msq2/MR2)*((320*pow2(Xt)*(-21*mU32*pow2(mQ32)*
         pow3(m32) - 21*mQ32*pow2(mU32)*pow3(m32) + 7*mU32*pow2(m32)*pow3(mQ32)
         - pow3(m32)*pow3(mQ32) + 7*mQ32*pow2(m32)*pow3(mU32) - pow3(m32)*pow3
         (mU32) + 42*mQ32*mU32*pow4(m32) + 3*pow2(mQ32)*pow4(m32) + 3*pow2(mU32
         )*pow4(m32) - 10*mQ32*pow5(m32) - 10*mU32*pow5(m32) + 2*pow6(m32)))/(
         3.*(mQ32 - mU32)*pow3(m32 - mQ32)*pow3(m32 - mU32)) - (160*(mQ32 +
         mU32)*pow4(Xt)*(-21*mU32*pow2(mQ32)*pow3(m32) - 21*mQ32*pow2(mU32)*
         pow3(m32) + 7*mU32*pow2(m32)*pow3(mQ32) - pow3(m32)*pow3(mQ32) + 7*
         mQ32*pow2(m32)*pow3(mU32) - pow3(m32)*pow3(mU32) + 42*mQ32*mU32*pow4(
         m32) + 3*pow2(mQ32)*pow4(m32) + 3*pow2(mU32)*pow4(m32) - 10*mQ32*pow5(
         m32) - 10*mU32*pow5(m32) + 2*pow6(m32)))/(3.*pow3(m32 - mQ32)*pow3(m32
         - mU32)*pow3(mQ32 - mU32)) + (16*(-210*mU32*pow2(mQ32)*pow3(m32) -
         237*mQ32*pow2(mU32)*pow3(m32) + 70*mU32*pow2(m32)*pow3(mQ32) - 10*pow3
         (m32)*pow3(mQ32) + 79*mQ32*pow2(m32)*pow3(mU32) - 19*pow3(m32)*pow3(
         mU32) + 447*mQ32*mU32*pow4(m32) + 30*pow2(mQ32)*pow4(m32) + 57*pow2(
         mU32)*pow4(m32) - 109*mQ32*pow5(m32) - 127*mU32*pow5(m32) + 29*pow6(
         m32)))/(3.*pow3(m32 - mQ32)*pow3(m32 - mU32)) - (320*(mQ32 + mU32)*
         pow5(Xt)*(-11*mQ32*mU32*pow3(m3) + 5*mQ32*pow5(m3) + 5*mU32*pow5(m3) +
         pow7(m3)))/(3.*pow2(m32 - mQ32)*pow2(m32 - mU32)*pow3(mQ32 - mU32)) +
         Xt*((96*(-2*m32 + mQ32 + mU32)*pow3(m3))/((m32 - mQ32)*(m32 - mU32)*(
         mQ32 - mU32)) + (32*(-119*mQ32*mU32*pow3(m3) + 59*mQ32*pow5(m3) + 59*
         mU32*pow5(m3) + pow7(m3)))/(3.*pow2(m32 - mQ32)*pow2(m32 - mU32))) +
         pow3(Xt)*((-192*(4*mQ32*mU32 - 2*m32*(mQ32 + mU32) + 2*pow2(m32) -
         pow2(mQ32) - pow2(mU32))*pow3(m3))/((m32 - mQ32)*(m32 - mU32)*pow3(
         mQ32 - mU32)) + (64*(-119*mQ32*mU32*pow3(m3) + 59*mQ32*pow5(m3) + 59*
         mU32*pow5(m3) + pow7(m3)))/(3.*(mQ32 - mU32)*pow2(m32 - mQ32)*pow2(m32
         - mU32)))) + (16*pow2(Xt)*(180*msq2*pow2(mQ32)*pow2(mU32)*pow3(m32) -
         90*msq2*pow2(m32)*pow2(mU32)*pow3(mQ32) - 270*msq2*mU32*pow3(m32)*
         pow3(mQ32) - 3513*pow2(mU32)*pow3(m32)*pow3(mQ32) - 90*msq2*pow2(m32)*
         pow2(mQ32)*pow3(mU32) - 270*mQ32*msq2*pow3(m32)*pow3(mU32) - 3993*pow2
         (mQ32)*pow3(m32)*pow3(mU32) + 2430*pow2(m32)*pow3(mQ32)*pow3(mU32) +
         240*msq2*mU32*pow2(mQ32)*pow4(m32) + 240*mQ32*msq2*pow2(mU32)*pow4(m32
         ) + 5696*pow2(mQ32)*pow2(mU32)*pow4(m32) + 1680*mU32*pow3(mQ32)*pow4(
         m32) + 2448*mQ32*pow3(mU32)*pow4(m32) + 90*msq2*mU32*pow2(m32)*pow4(
         mQ32) + 30*m32*msq2*pow2(mU32)*pow4(mQ32) + 732*pow2(m32)*pow2(mU32)*
         pow4(mQ32) - 391*mU32*pow3(m32)*pow4(mQ32) - 489*m32*pow3(mU32)*pow4(
         mQ32) - 36*pow4(m32)*pow4(mQ32) + 90*mQ32*msq2*pow2(m32)*pow4(mU32) +
         30*m32*msq2*pow2(mQ32)*pow4(mU32) + 852*pow2(m32)*pow2(mQ32)*pow4(mU32
         ) - 591*mQ32*pow3(m32)*pow4(mU32) - 513*m32*pow3(mQ32)*pow4(mU32) + 4*
         pow4(m32)*pow4(mU32) + 96*pow4(mQ32)*pow4(mU32) - 180*mQ32*msq2*mU32*
         pow5(m32) - 2650*mU32*pow2(mQ32)*pow5(m32) - 3290*mQ32*pow2(mU32)*pow5
         (m32) + 108*pow3(mQ32)*pow5(m32) - 12*pow3(mU32)*pow5(m32) + 9*mU32*
         pow2(m32)*pow5(mQ32) + 3*m32*pow2(mU32)*pow5(mQ32) + 9*mQ32*pow2(m32)*
         pow5(mU32) + 3*m32*pow2(mQ32)*pow5(mU32) + 1472*mQ32*mU32*pow6(m32) -
         76*pow2(mQ32)*pow6(m32) + 12*pow2(mU32)*pow6(m32) + 4*mQ32*pow7(m32) -
         4*mU32*pow7(m32)))/(3.*mQ32*(mQ32 - mU32)*mU32*pow3(m32 - mQ32)*pow3(
         m32 - mU32)) - (8*pow4(Xt)*(-90*msq2*pow3(m32)*pow3(mQ32)*pow3(mU32) -
         570*msq2*pow2(mU32)*pow3(mQ32)*pow4(m32) - 330*msq2*pow2(mQ32)*pow3(
         mU32)*pow4(m32) - 25654*pow3(mQ32)*pow3(mU32)*pow4(m32) + 90*msq2*pow2
         (mU32)*pow3(m32)*pow4(mQ32) + 210*msq2*pow2(m32)*pow3(mU32)*pow4(mQ32)
         + 18064*pow3(m32)*pow3(mU32)*pow4(mQ32) - 510*msq2*mU32*pow4(m32)*
         pow4(mQ32) - 22956*pow2(mU32)*pow4(m32)*pow4(mQ32) + 270*msq2*pow2(
         mQ32)*pow3(m32)*pow4(mU32) + 30*msq2*pow2(m32)*pow3(mQ32)*pow4(mU32) +
         10886*pow3(m32)*pow3(mQ32)*pow4(mU32) - 270*mQ32*msq2*pow4(m32)*pow4(
         mU32) - 9712*pow2(mQ32)*pow4(m32)*pow4(mU32) - 30*m32*msq2*pow4(mQ32)*
         pow4(mU32) - 5480*pow2(m32)*pow4(mQ32)*pow4(mU32) + 660*msq2*pow2(mQ32
         )*pow2(mU32)*pow5(m32) + 420*msq2*mU32*pow3(mQ32)*pow5(m32) + 25568*
         pow2(mU32)*pow3(mQ32)*pow5(m32) + 240*mQ32*msq2*pow3(mU32)*pow5(m32) +
         17622*pow2(mQ32)*pow3(mU32)*pow5(m32) + 9890*mU32*pow4(mQ32)*pow5(m32
         ) + 3020*mQ32*pow4(mU32)*pow5(m32) + 30*msq2*pow2(m32)*pow2(mU32)*pow5
         (mQ32) + 360*msq2*mU32*pow3(m32)*pow5(mQ32) + 9529*pow2(mU32)*pow3(m32
         )*pow5(mQ32) - 30*m32*msq2*pow3(mU32)*pow5(mQ32) - 6592*pow2(m32)*pow3
         (mU32)*pow5(mQ32) - 4599*mU32*pow4(m32)*pow5(mQ32) + 1602*m32*pow4(
         mU32)*pow5(mQ32) - 144*pow5(m32)*pow5(mQ32) - 60*msq2*pow2(m32)*pow2(
         mQ32)*pow5(mU32) + 90*mQ32*msq2*pow3(m32)*pow5(mU32) + 1840*pow2(mQ32)
         *pow3(m32)*pow5(mU32) - 30*m32*msq2*pow3(mQ32)*pow5(mU32) - 1731*pow2(
         m32)*pow3(mQ32)*pow5(mU32) - 651*mQ32*pow4(m32)*pow5(mU32) + 618*m32*
         pow4(mQ32)*pow5(mU32) + 4*pow5(m32)*pow5(mU32) - 96*pow5(mQ32)*pow5(
         mU32) - 180*msq2*mU32*pow2(mQ32)*pow6(m32) - 180*mQ32*msq2*pow2(mU32)*
         pow6(m32) - 12968*pow2(mQ32)*pow2(mU32)*pow6(m32) - 9086*mU32*pow3(
         mQ32)*pow6(m32) - 4430*mQ32*pow3(mU32)*pow6(m32) + 184*pow4(mQ32)*pow6
         (m32) - 12*pow4(mU32)*pow6(m32) - 90*msq2*mU32*pow2(m32)*pow6(mQ32) -
         30*m32*msq2*pow2(mU32)*pow6(mQ32) - 1262*pow2(m32)*pow2(mU32)*pow6(
         mQ32) + 696*mU32*pow3(m32)*pow6(mQ32) + 882*m32*pow3(mU32)*pow6(mQ32)
         + 36*pow4(m32)*pow6(mQ32) - 240*pow4(mU32)*pow6(mQ32) - 6*pow2(m32)*
         pow2(mQ32)*pow6(mU32) + 9*mQ32*pow3(m32)*pow6(mU32) - 3*m32*pow3(mQ32)
         *pow6(mU32) + 3172*mU32*pow2(mQ32)*pow7(m32) + 2160*mQ32*pow2(mU32)*
         pow7(m32) - 80*pow3(mQ32)*pow7(m32) + 12*pow3(mU32)*pow7(m32) - 9*mU32
         *pow2(m32)*pow7(mQ32) - 3*m32*pow2(mU32)*pow7(mQ32) - 80*mQ32*mU32*
         pow8(m32) + 4*pow2(mQ32)*pow8(m32) - 4*pow2(mU32)*pow8(m32)))/(3.*mQ32
         *mU32*pow3(m32 - mU32)*pow3(mQ32 - mU32)*pow4(m32 - mQ32)) - (8*(-450*
         msq2*pow3(m32)*pow3(mQ32)*pow3(mU32) - 450*msq2*pow2(mU32)*pow3(mQ32)*
         pow4(m32) + 210*msq2*pow2(mQ32)*pow3(mU32)*pow4(m32) - 1590*pow3(mQ32)
         *pow3(mU32)*pow4(m32) + 630*msq2*pow2(mU32)*pow3(m32)*pow4(mQ32) + 30*
         msq2*pow2(m32)*pow3(mU32)*pow4(mQ32) - 1684*pow3(m32)*pow3(mU32)*pow4(
         mQ32) + 510*msq2*mU32*pow4(m32)*pow4(mQ32) + 5488*pow2(mU32)*pow4(m32)
         *pow4(mQ32) + 90*msq2*pow2(mQ32)*pow3(m32)*pow4(mU32) + 150*msq2*pow2(
         m32)*pow3(mQ32)*pow4(mU32) + 3694*pow3(m32)*pow3(mQ32)*pow4(mU32) -
         270*mQ32*msq2*pow4(m32)*pow4(mU32) - 3520*pow2(mQ32)*pow4(m32)*pow4(
         mU32) + 30*m32*msq2*pow4(mQ32)*pow4(mU32) - 1524*pow2(m32)*pow4(mQ32)*
         pow4(mU32) + 180*msq2*pow2(mQ32)*pow2(mU32)*pow5(m32) - 420*msq2*mU32*
         pow3(mQ32)*pow5(m32) - 3622*pow2(mU32)*pow3(mQ32)*pow5(m32) + 240*mQ32
         *msq2*pow3(mU32)*pow5(m32) + 3060*pow2(mQ32)*pow3(mU32)*pow5(m32) -
         3764*mU32*pow4(mQ32)*pow5(m32) + 1906*mQ32*pow4(mU32)*pow5(m32) - 210*
         msq2*pow2(m32)*pow2(mU32)*pow5(mQ32) - 360*msq2*mU32*pow3(m32)*pow5(
         mQ32) - 3331*pow2(mU32)*pow3(m32)*pow5(mQ32) - 30*m32*msq2*pow3(mU32)*
         pow5(mQ32) + 2028*pow2(m32)*pow3(mU32)*pow5(mQ32) + 2055*mU32*pow4(m32
         )*pow5(mQ32) - 84*m32*pow4(mU32)*pow5(mQ32) - 144*pow5(m32)*pow5(mQ32)
         - 60*msq2*pow2(m32)*pow2(mQ32)*pow5(mU32) + 90*mQ32*msq2*pow3(m32)*
         pow5(mU32) + 970*pow2(mQ32)*pow3(m32)*pow5(mU32) - 30*m32*msq2*pow3(
         mQ32)*pow5(mU32) - 1017*pow2(m32)*pow3(mQ32)*pow5(mU32) - 549*mQ32*
         pow4(m32)*pow5(mU32) + 540*m32*pow4(mQ32)*pow5(mU32) + 4*pow5(m32)*
         pow5(mU32) - 84*pow5(mQ32)*pow5(mU32) + 180*msq2*mU32*pow2(mQ32)*pow6(
         m32) - 180*mQ32*msq2*pow2(mU32)*pow6(m32) + 616*pow2(mQ32)*pow2(mU32)*
         pow6(m32) + 3318*mU32*pow3(mQ32)*pow6(m32) - 2186*mQ32*pow3(mU32)*pow6
         (m32) + 184*pow4(mQ32)*pow6(m32) - 12*pow4(mU32)*pow6(m32) + 90*msq2*
         mU32*pow2(m32)*pow6(mQ32) + 30*m32*msq2*pow2(mU32)*pow6(mQ32) + 638*
         pow2(m32)*pow2(mU32)*pow6(mQ32) - 426*mU32*pow3(m32)*pow6(mQ32) - 456*
         m32*pow3(mU32)*pow6(mQ32) + 36*pow4(m32)*pow6(mQ32) + 84*pow4(mU32)*
         pow6(mQ32) - 6*pow2(m32)*pow2(mQ32)*pow6(mU32) + 9*mQ32*pow3(m32)*pow6
         (mU32) - 3*m32*pow3(mQ32)*pow6(mU32) - 1456*mU32*pow2(mQ32)*pow7(m32)
         + 756*mQ32*pow2(mU32)*pow7(m32) - 80*pow3(mQ32)*pow7(m32) + 12*pow3(
         mU32)*pow7(m32) + 9*mU32*pow2(m32)*pow7(mQ32) + 3*m32*pow2(mU32)*pow7(
         mQ32) + 128*mQ32*mU32*pow8(m32) + 4*pow2(mQ32)*pow8(m32) - 4*pow2(mU32
         )*pow8(m32)))/(3.*mQ32*(mQ32 - mU32)*mU32*pow3(m32 - mU32)*pow4(m32 -
         mQ32)) - (64*pow5(Xt)*(16*mQ32*pow11(m3) + 52*mU32*pow11(m3) + 60*msq2
         *pow2(mQ32)*pow2(mU32)*pow3(m3) + 30*msq2*mU32*pow3(m3)*pow3(mQ32) +
         179*pow2(mU32)*pow3(m3)*pow3(mQ32) + 30*mQ32*msq2*pow3(m3)*pow3(mU32)
         + 213*pow2(mQ32)*pow3(m3)*pow3(mU32) + 12*m3*pow3(mQ32)*pow3(mU32) + 3
         *mU32*pow3(m3)*pow4(mQ32) + 3*mQ32*pow3(m3)*pow4(mU32) - 90*msq2*mU32*
         pow2(mQ32)*pow5(m3) - 120*mQ32*msq2*pow2(mU32)*pow5(m3) - 737*pow2(
         mQ32)*pow2(mU32)*pow5(m3) - 219*mU32*pow3(mQ32)*pow5(m3) - 407*mQ32*
         pow3(mU32)*pow5(m3) - 30*msq2*pow3(mU32)*pow5(m3) - 3*pow4(mU32)*pow5(
         m3) + 60*mQ32*msq2*mU32*pow7(m3) + 558*mU32*pow2(mQ32)*pow7(m3) + 770*
         mQ32*pow2(mU32)*pow7(m3) + 60*msq2*pow2(mU32)*pow7(m3) + 16*pow3(mQ32)
         *pow7(m3) + 214*pow3(mU32)*pow7(m3) - 362*mQ32*mU32*pow9(m3) - 32*pow2
         (mQ32)*pow9(m3) - 276*pow2(mU32)*pow9(m3)))/(3.*mU32*pow2(m32 - mU32)*
         pow3(m32 - mQ32)*pow3(mQ32 - mU32)) - (32*Xt*(32*mQ32*pow11(m3) + 8*
         mU32*pow11(m3) - 60*msq2*mU32*pow3(m3)*pow3(mQ32) - 69*pow2(mU32)*pow3
         (m3)*pow3(mQ32) + 60*mQ32*msq2*pow3(m3)*pow3(mU32) + 73*pow2(mQ32)*
         pow3(m3)*pow3(mU32) - 24*m3*pow3(mQ32)*pow3(mU32) - 6*mU32*pow3(m3)*
         pow4(mQ32) + 6*mQ32*pow3(m3)*pow4(mU32) + 180*msq2*mU32*pow2(mQ32)*
         pow5(m3) - 120*mQ32*msq2*pow2(mU32)*pow5(m3) + 221*pow2(mQ32)*pow2(
         mU32)*pow5(m3) + 171*mU32*pow3(mQ32)*pow5(m3) - 134*mQ32*pow3(mU32)*
         pow5(m3) - 60*msq2*pow3(mU32)*pow5(m3) - 6*pow4(mU32)*pow5(m3) - 120*
         mQ32*msq2*mU32*pow7(m3) - 444*mU32*pow2(mQ32)*pow7(m3) - 213*mQ32*pow2
         (mU32)*pow7(m3) + 120*msq2*pow2(mU32)*pow7(m3) + 32*pow3(mQ32)*pow7(m3
         ) + 245*pow3(mU32)*pow7(m3) + 431*mQ32*mU32*pow9(m3) - 64*pow2(mQ32)*
         pow9(m3) - 259*pow2(mU32)*pow9(m3)))/(3.*(mQ32 - mU32)*mU32*pow2(m32 -
         mU32)*pow3(m32 - mQ32)) + (128*pow3(Xt)*(22*mU32*pow11(m3) + 30*msq2*
         pow2(mQ32)*pow2(mU32)*pow3(m3) - 30*msq2*mU32*pow3(m3)*pow3(mQ32) -
         115*pow2(mU32)*pow3(m3)*pow3(mQ32) + 30*mQ32*msq2*pow3(m3)*pow3(mU32)
         + 106*pow2(mQ32)*pow3(m3)*pow3(mU32) + 18*m3*pow3(mQ32)*pow3(mU32) - 3
         *mU32*pow3(m3)*pow4(mQ32) + 6*m3*pow2(mQ32)*pow4(mU32) - 91*mQ32*pow3(
         m3)*pow4(mU32) - 30*msq2*pow3(m3)*pow4(mU32) + 60*msq2*mU32*pow2(mQ32)
         *pow5(m3) - 120*mQ32*msq2*pow2(mU32)*pow5(m3) - 17*pow2(mQ32)*pow2(
         mU32)*pow5(m3) + 165*mU32*pow3(mQ32)*pow5(m3) - 85*mQ32*pow3(mU32)*
         pow5(m3) + 60*msq2*pow3(mU32)*pow5(m3) + 153*pow4(mU32)*pow5(m3) - 3*
         pow3(m3)*pow5(mU32) - 217*mU32*pow2(mQ32)*pow7(m3) + 238*mQ32*pow2(
         mU32)*pow7(m3) - 16*pow3(mQ32)*pow7(m3) - 193*pow3(mU32)*pow7(m3) + 6*
         mQ32*mU32*pow9(m3) + 16*pow2(mQ32)*pow9(m3) + 10*pow2(mU32)*pow9(m3)))
         /(3.*mU32*pow2(m32 - mQ32)*pow2(m32 - mU32)*pow3(mQ32 - mU32)) + log(
         mU32/MR2)*((2560*(mQ32 + mU32)*(m32*mQ32 + m32*mU32 - 2*mQ32*mU32)*
         pow2(m32)*pow6(Xt))/(3.*pow2(m32 - mQ32)*pow2(m32 - mU32)*pow4(mQ32 -
         mU32)) + (16*(-252*pow2(mU32)*pow3(m32)*pow3(mQ32) - 248*pow2(mQ32)*
         pow3(m32)*pow3(mU32) - 132*pow2(m32)*pow3(mQ32)*pow3(mU32) + 1166*pow2
         (mQ32)*pow2(mU32)*pow4(m32) + 452*mU32*pow3(mQ32)*pow4(m32) + 444*mQ32
         *pow3(mU32)*pow4(m32) + 18*pow2(m32)*pow2(mU32)*pow4(mQ32) - 68*mU32*
         pow3(m32)*pow4(mQ32) + 48*m32*pow3(mU32)*pow4(mQ32) + 18*pow4(m32)*
         pow4(mQ32) + 17*pow2(m32)*pow2(mQ32)*pow4(mU32) - 66*mQ32*pow3(m32)*
         pow4(mU32) + 48*m32*pow3(mQ32)*pow4(mU32) + 17*pow4(m32)*pow4(mU32) -
         12*pow4(mQ32)*pow4(mU32) - 1272*mU32*pow2(mQ32)*pow5(m32) - 1264*mQ32*
         pow2(mU32)*pow5(m32) - 132*pow3(mQ32)*pow5(m32) - 128*pow3(mU32)*pow5(
         m32) + 1204*mQ32*mU32*pow6(m32) + 361*pow2(mQ32)*pow6(m32) + 356*pow2(
         mU32)*pow6(m32) - 334*mQ32*pow7(m32) - 332*mU32*pow7(m32) + 91*pow8(
         m32)))/(3.*pow4(m32 - mQ32)*pow4(m32 - mU32)) + (64*Xt*(-43*mQ32*pow11
         (m3) + 38*mU32*pow11(m3) + pow13(m3) - 73*pow2(mU32)*pow3(m3)*pow3(
         mQ32) + 72*pow2(mQ32)*pow3(m3)*pow3(mU32) + 3*pow2(mQ32)*pow2(mU32)*
         pow5(m3) + 98*mU32*pow3(mQ32)*pow5(m3) - 96*mQ32*pow3(mU32)*pow5(m3) -
         94*mU32*pow2(mQ32)*pow7(m3) + 85*mQ32*pow2(mU32)*pow7(m3) - 41*pow3(
         mQ32)*pow7(m3) + 40*pow3(mU32)*pow7(m3) + 6*mQ32*mU32*pow9(m3) + 67*
         pow2(mQ32)*pow9(m3) - 63*pow2(mU32)*pow9(m3)))/(3.*(mQ32 - mU32)*pow3(
         m32 - mQ32)*pow3(m32 - mU32)) + (128*pow3(Xt)*(9*mQ32*mU32*pow11(m3) +
         33*mQ32*pow13(m3) - 39*mU32*pow13(m3) + pow15(m3) - 147*pow11(m3)*
         pow2(mQ32) + 153*pow11(m3)*pow2(mU32) + pow3(m3)*pow3(mQ32)*pow3(mU32)
         - 54*pow2(mU32)*pow3(m3)*pow4(mQ32) - 12*m3*pow3(mU32)*pow4(mQ32) +
         54*pow2(mQ32)*pow3(m3)*pow4(mU32) + 12*m3*pow3(mQ32)*pow4(mU32) + 279*
         pow2(mU32)*pow3(mQ32)*pow5(m3) - 285*pow2(mQ32)*pow3(mU32)*pow5(m3) +
         102*mU32*pow4(mQ32)*pow5(m3) - 102*mQ32*pow4(mU32)*pow5(m3) + 9*pow2(
         mQ32)*pow2(mU32)*pow7(m3) - 409*mU32*pow3(mQ32)*pow7(m3) + 415*mQ32*
         pow3(mU32)*pow7(m3) - 52*pow4(mQ32)*pow7(m3) + 52*pow4(mU32)*pow7(m3)
         + 273*mU32*pow2(mQ32)*pow9(m3) - 291*mQ32*pow2(mU32)*pow9(m3) + 181*
         pow3(mQ32)*pow9(m3) - 183*pow3(mU32)*pow9(m3)))/(3.*pow3(m32 - mQ32)*
         pow3(m32 - mU32)*pow3(mQ32 - mU32)) - (128*pow5(Xt)*(-214*mQ32*mU32*
         pow11(m3) + 34*mQ32*pow13(m3) + 34*mU32*pow13(m3) - 107*pow11(m3)*pow2
         (mQ32) - 107*pow11(m3)*pow2(mU32) - 70*pow3(m3)*pow3(mQ32)*pow3(mU32)
         - 35*pow2(mU32)*pow3(m3)*pow4(mQ32) - 6*m3*pow3(mU32)*pow4(mQ32) - 35*
         pow2(mQ32)*pow3(m3)*pow4(mU32) - 6*m3*pow3(mQ32)*pow4(mU32) + 291*pow2
         (mU32)*pow3(mQ32)*pow5(m3) + 291*pow2(mQ32)*pow3(mU32)*pow5(m3) + 67*
         mU32*pow4(mQ32)*pow5(m3) + 67*mQ32*pow4(mU32)*pow5(m3) - 560*pow2(mQ32
         )*pow2(mU32)*pow7(m3) - 314*mU32*pow3(mQ32)*pow7(m3) - 314*mQ32*pow3(
         mU32)*pow7(m3) - 34*pow4(mQ32)*pow7(m3) - 34*pow4(mU32)*pow7(m3) + 411
         *mU32*pow2(mQ32)*pow9(m3) + 411*mQ32*pow2(mU32)*pow9(m3) + 115*pow3(
         mQ32)*pow9(m3) + 115*pow3(mU32)*pow9(m3)))/(3.*pow3(m32 - mQ32)*pow3(
         m32 - mU32)*pow4(mQ32 - mU32)) - (64*pow2(Xt)*(776*pow3(m32)*pow3(mQ32
         )*pow3(mU32) - 1704*pow2(mU32)*pow3(mQ32)*pow4(m32) - 1704*pow2(mQ32)*
         pow3(mU32)*pow4(m32) + 564*pow2(mU32)*pow3(m32)*pow4(mQ32) - 199*pow2(
         m32)*pow3(mU32)*pow4(mQ32) - 561*mU32*pow4(m32)*pow4(mQ32) + 564*pow2(
         mQ32)*pow3(m32)*pow4(mU32) - 199*pow2(m32)*pow3(mQ32)*pow4(mU32) - 561
         *mQ32*pow4(m32)*pow4(mU32) + 48*m32*pow4(mQ32)*pow4(mU32) + 2808*pow2(
         mQ32)*pow2(mU32)*pow5(m32) + 1488*mU32*pow3(mQ32)*pow5(m32) + 1488*
         mQ32*pow3(mU32)*pow5(m32) + 164*pow4(mQ32)*pow5(m32) + 164*pow4(mU32)*
         pow5(m32) - 81*pow2(m32)*pow2(mU32)*pow5(mQ32) + 88*mU32*pow3(m32)*
         pow5(mQ32) + 24*m32*pow3(mU32)*pow5(mQ32) - 27*pow4(m32)*pow5(mQ32) -
         6*pow4(mU32)*pow5(mQ32) - 81*pow2(m32)*pow2(mQ32)*pow5(mU32) + 88*mQ32
         *pow3(m32)*pow5(mU32) + 24*m32*pow3(mQ32)*pow5(mU32) - 27*pow4(m32)*
         pow5(mU32) - 6*pow4(mQ32)*pow5(mU32) - 2083*mU32*pow2(mQ32)*pow6(m32)
         - 2083*mQ32*pow2(mU32)*pow6(m32) - 405*pow3(mQ32)*pow6(m32) - 405*pow3
         (mU32)*pow6(m32) + 1368*mQ32*mU32*pow7(m32) + 516*pow2(mQ32)*pow7(m32)
         + 516*pow2(mU32)*pow7(m32) - 310*mQ32*pow8(m32) - 310*mU32*pow8(m32)
         + 64*pow9(m32)))/(3.*pow2(mQ32 - mU32)*pow4(m32 - mQ32)*pow4(m32 -
         mU32)) + (32*pow4(Xt)*(-11360*pow3(mQ32)*pow3(mU32)*pow4(m32) + 3884*
         pow3(m32)*pow3(mU32)*pow4(mQ32) - 3615*pow2(mU32)*pow4(m32)*pow4(mQ32)
         + 3884*pow3(m32)*pow3(mQ32)*pow4(mU32) - 3615*pow2(mQ32)*pow4(m32)*
         pow4(mU32) - 1702*pow2(m32)*pow4(mQ32)*pow4(mU32) + 12176*pow2(mU32)*
         pow3(mQ32)*pow5(m32) + 12176*pow2(mQ32)*pow3(mU32)*pow5(m32) + 1612*
         mU32*pow4(mQ32)*pow5(m32) + 1612*mQ32*pow4(mU32)*pow5(m32) - 44*pow2(
         mU32)*pow3(m32)*pow5(mQ32) - 80*pow2(m32)*pow3(mU32)*pow5(mQ32) - 24*
         mU32*pow4(m32)*pow5(mQ32) + 144*m32*pow4(mU32)*pow5(mQ32) + 36*pow5(
         m32)*pow5(mQ32) - 44*pow2(mQ32)*pow3(m32)*pow5(mU32) - 80*pow2(m32)*
         pow3(mQ32)*pow5(mU32) - 24*mQ32*pow4(m32)*pow5(mU32) + 144*m32*pow4(
         mQ32)*pow5(mU32) + 36*pow5(m32)*pow5(mU32) - 24*pow5(mQ32)*pow5(mU32)
         - 13598*pow2(mQ32)*pow2(mU32)*pow6(m32) - 5712*mU32*pow3(mQ32)*pow6(
         m32) - 5712*mQ32*pow3(mU32)*pow6(m32) - 321*pow4(mQ32)*pow6(m32) - 321
         *pow4(mU32)*pow6(m32) - 125*pow2(m32)*pow2(mU32)*pow6(mQ32) + 128*mU32
         *pow3(m32)*pow6(mQ32) + 48*m32*pow3(mU32)*pow6(mQ32) - 41*pow4(m32)*
         pow6(mQ32) - 12*pow4(mU32)*pow6(mQ32) - 125*pow2(m32)*pow2(mQ32)*pow6(
         mU32) + 128*mQ32*pow3(m32)*pow6(mU32) + 48*m32*pow3(mQ32)*pow6(mU32) -
         41*pow4(m32)*pow6(mU32) - 12*pow4(mQ32)*pow6(mU32) + 6252*mU32*pow2(
         mQ32)*pow7(m32) + 6252*mQ32*pow2(mU32)*pow7(m32) + 1044*pow3(mQ32)*
         pow7(m32) + 1044*pow3(mU32)*pow7(m32) - 2584*mQ32*mU32*pow8(m32) -
         1036*pow2(mQ32)*pow8(m32) - 1036*pow2(mU32)*pow8(m32) + 320*mQ32*pow9(
         m32) + 320*mU32*pow9(m32)))/(3.*pow4(m32 - mQ32)*pow4(m32 - mU32)*pow4
         (mQ32 - mU32)))));

   return result;
}

/// 3-loop coefficient O(at*as^2*log^1)
double himalaya::mh2_eft::Mh2EFTCalculator::coeff_as_2_log_1(double mQ32, double mU32, double Xt, double m3, double msq2, double MR2)
{
   using std::log;
   using gm2calc::dilog;
   const double m32 = pow2(m3);

   const double result =
      -128*zt3 + Xt*((-128*m3*mQ32)/((m32 - mQ32)*(m32 - mU32)) - (
         1280*m3*msq2)/((m32 - mQ32)*(m32 - mU32)) - (128*m3*mU32)/((m32 - mQ32
         )*(m32 - mU32)) + (1792*pow3(m3))/(3.*(m32 - mQ32)*(m32 - mU32))) +
         dilog(1 - msq2/mQ32)*((-1280*m3*Xt*pow2(mQ32 - msq2))/((mQ32 - mU32)*
         pow2(m32 - mQ32)) + (80*(2*mQ32 - 2*msq2)*(-3*m32*mQ32 + 3*m32*msq2 +
         mQ32*msq2 - 2*pow2(m32) + pow2(mQ32)))/pow3(-m32 + mQ32)) + dilog(1 -
         m32/msq2)*((-1280*(m32 - msq2)*Xt*(m3*mQ32*msq2 - 2*m3*mQ32*mU32 + m3*
         msq2*mU32 + mQ32*pow3(m3) - 2*msq2*pow3(m3) + mU32*pow3(m3)))/(pow2(
         m32 - mQ32)*pow2(m32 - mU32)) - 24*(2*m32 - 2*msq2)*(-10/(3.*(-m32 +
         mQ32)) - 10/(3.*(-m32 + mU32)) - (10*mQ32)/pow2(-m32 + mQ32) + (10*
         msq2)/pow2(-m32 + mQ32) + (10*msq2)/pow2(-m32 + mU32) - (10*mU32)/pow2
         (-m32 + mU32) - (40*mQ32*msq2)/(3.*pow3(-m32 + mQ32)) + (40*pow2(mQ32)
         )/(3.*pow3(-m32 + mQ32)) - (40*msq2*mU32)/(3.*pow3(-m32 + mU32)) + (40
         *pow2(mU32))/(3.*pow3(-m32 + mU32)))) + pow2(log(msq2/MR2))*(-48*Xt*((
         40*m3*pow2(mQ32 - msq2))/(3.*(mQ32 - mU32)*pow2(m32 - mQ32)) + (40*m3*
         pow2(-msq2 + mU32))/(3.*(-mQ32 + mU32)*pow2(m32 - mU32))) - 48*((-5*(2
         *mQ32 - 2*msq2)*(-3*m32*mQ32 + 3*m32*msq2 + mQ32*msq2 - 2*pow2(m32) +
         pow2(mQ32)))/(6.*pow3(-m32 + mQ32)) - (5*(-2*msq2 + 2*mU32)*(3*m32*
         msq2 - 3*m32*mU32 + msq2*mU32 - 2*pow2(m32) + pow2(mU32)))/(6.*pow3(
         -m32 + mU32)))) + dilog(1 - msq2/mU32)*((-1280*m3*Xt*pow2(-msq2 + mU32
         ))/((-mQ32 + mU32)*pow2(m32 - mU32)) + (80*(-2*msq2 + 2*mU32)*(3*m32*
         msq2 - 3*m32*mU32 + msq2*mU32 - 2*pow2(m32) + pow2(mU32)))/pow3(-m32 +
         mU32)) - (16*(120*m32*mQ32*msq2*mU32 + 150*mQ32*msq2*pow2(m32) + 454*
         mQ32*mU32*pow2(m32) + 150*msq2*mU32*pow2(m32) - 90*m32*msq2*pow2(mQ32)
         - 211*m32*mU32*pow2(mQ32) - 30*msq2*mU32*pow2(mQ32) + 169*pow2(m32)*
         pow2(mQ32) - 211*m32*mQ32*pow2(mU32) - 90*m32*msq2*pow2(mU32) - 30*
         mQ32*msq2*pow2(mU32) + 169*pow2(m32)*pow2(mU32) + 98*pow2(mQ32)*pow2(
         mU32) - 316*mQ32*pow3(m32) - 180*msq2*pow3(m32) - 316*mU32*pow3(m32) -
         9*m32*pow3(mQ32) - 3*mU32*pow3(mQ32) - 9*m32*pow3(mU32) - 3*mQ32*pow3
         (mU32) + 188*pow4(m32)))/(3.*pow2(m32 - mQ32)*pow2(m32 - mU32)) + log(
         msq2/MR2)*((-1280*m3*msq2*Xt)/((m32 - mQ32)*(m32 - mU32)) + (16*(-120*
         m32*mQ32*msq2*mU32 - 150*mQ32*msq2*pow2(m32) + 404*mQ32*mU32*pow2(m32)
         - 150*msq2*mU32*pow2(m32) + 90*m32*msq2*pow2(mQ32) - 201*m32*mU32*
         pow2(mQ32) + 30*msq2*mU32*pow2(mQ32) + 101*pow2(m32)*pow2(mQ32) - 201*
         m32*mQ32*pow2(mU32) + 90*m32*msq2*pow2(mU32) + 30*mQ32*msq2*pow2(mU32)
         + 101*pow2(m32)*pow2(mU32) + 100*pow2(mQ32)*pow2(mU32) - 203*mQ32*
         pow3(m32) + 180*msq2*pow3(m32) - 203*mU32*pow3(m32) + 102*pow4(m32)))/
         (3.*pow2(m32 - mQ32)*pow2(m32 - mU32))) + dilog(1 - mQ32/mU32)*((64*Xt
         *(-8*m3*mU32*pow2(mQ32) + 8*m3*mQ32*pow2(mU32) - 10*pow2(mQ32)*pow3(m3
         ) + 10*pow2(mU32)*pow3(m3) + 6*m3*pow3(mQ32) - 6*m3*pow3(mU32) + 10*
         mQ32*pow5(m3) - 10*mU32*pow5(m3)))/(3.*pow2(m32 - mQ32)*pow2(m32 -
         mU32)) - (8*(220*mU32*pow2(mQ32)*pow3(m32) - 220*mQ32*pow2(mU32)*pow3(
         m32) - 188*mU32*pow2(m32)*pow3(mQ32) + 16*m32*pow2(mU32)*pow3(mQ32) -
         68*pow3(m32)*pow3(mQ32) + 188*mQ32*pow2(m32)*pow3(mU32) - 16*m32*pow2(
         mQ32)*pow3(mU32) + 68*pow3(m32)*pow3(mU32) + 28*pow2(mQ32)*pow4(m32) -
         28*pow2(mU32)*pow4(m32) + 70*m32*mU32*pow4(mQ32) + 58*pow2(m32)*pow4(
         mQ32) - 8*pow2(mU32)*pow4(mQ32) - 70*m32*mQ32*pow4(mU32) - 58*pow2(m32
         )*pow4(mU32) + 8*pow2(mQ32)*pow4(mU32) - 24*mQ32*pow5(m32) + 24*mU32*
         pow5(m32) - 18*m32*pow5(mQ32) - 6*mU32*pow5(mQ32) + 18*m32*pow5(mU32)
         + 6*mQ32*pow5(mU32)))/(3.*pow3(-m32 + mQ32)*pow3(m32 - mU32))) + dilog
         (1 - m32/mQ32)*((-128*Xt*(-7*m3*mQ32*mU32 + 3*m3*pow2(mQ32) + 22*m3*
         pow2(mU32) + mQ32*pow3(m3) - 37*mU32*pow3(m3) + 18*pow5(m3)))/(3.*(
         -mQ32 + mU32)*pow2(m32 - mU32)) + (16*(-417*pow2(m32)*pow2(mQ32)*pow2(
         mU32) + 251*mU32*pow2(mQ32)*pow3(m32) + 1025*mQ32*pow2(mU32)*pow3(m32)
         + 151*mU32*pow2(m32)*pow3(mQ32) - 41*m32*pow2(mU32)*pow3(mQ32) - 9*
         pow3(m32)*pow3(mQ32) - 391*mQ32*pow2(m32)*pow3(mU32) + 167*m32*pow2(
         mQ32)*pow3(mU32) + 13*pow3(m32)*pow3(mU32) + 19*pow3(mQ32)*pow3(mU32)
         - 902*mQ32*mU32*pow4(m32) - 98*pow2(mQ32)*pow4(m32) - 280*pow2(mU32)*
         pow4(m32) - 41*m32*mU32*pow4(mQ32) - 20*pow2(m32)*pow4(mQ32) + pow2(
         mU32)*pow4(mQ32) + 34*m32*mQ32*pow4(mU32) + 37*pow2(m32)*pow4(mU32) -
         23*pow2(mQ32)*pow4(mU32) + 294*mQ32*pow5(m32) + 346*mU32*pow5(m32) + 9
         *m32*pow5(mQ32) + 3*mU32*pow5(mQ32) - 128*pow6(m32)))/(3.*(mQ32 - mU32
         )*pow2(m32 - mQ32)*pow3(m32 - mU32))) + dilog(1 - m32/mU32)*((-128*Xt*
         (-7*m3*mQ32*mU32 + 22*m3*pow2(mQ32) + 3*m3*pow2(mU32) - 37*mQ32*pow3(
         m3) + mU32*pow3(m3) + 18*pow5(m3)))/(3.*(mQ32 - mU32)*pow2(m32 - mQ32)
         ) - (16*(-417*pow2(m32)*pow2(mQ32)*pow2(mU32) + 1025*mU32*pow2(mQ32)*
         pow3(m32) + 251*mQ32*pow2(mU32)*pow3(m32) - 391*mU32*pow2(m32)*pow3(
         mQ32) + 167*m32*pow2(mU32)*pow3(mQ32) + 13*pow3(m32)*pow3(mQ32) + 151*
         mQ32*pow2(m32)*pow3(mU32) - 41*m32*pow2(mQ32)*pow3(mU32) - 9*pow3(m32)
         *pow3(mU32) + 19*pow3(mQ32)*pow3(mU32) - 902*mQ32*mU32*pow4(m32) - 280
         *pow2(mQ32)*pow4(m32) - 98*pow2(mU32)*pow4(m32) + 34*m32*mU32*pow4(
         mQ32) + 37*pow2(m32)*pow4(mQ32) - 23*pow2(mU32)*pow4(mQ32) - 41*m32*
         mQ32*pow4(mU32) - 20*pow2(m32)*pow4(mU32) + pow2(mQ32)*pow4(mU32) +
         346*mQ32*pow5(m32) + 294*mU32*pow5(m32) + 9*m32*pow5(mU32) + 3*mQ32*
         pow5(mU32) - 128*pow6(m32)))/(3.*(mQ32 - mU32)*pow2(m32 - mU32)*pow3(
         m32 - mQ32))) + log(mU32/MR2)*(log(msq2/MR2)*((-64*Xt*(120*m3*msq2*
         mU32 - 60*m3*pow2(msq2) - m3*pow2(mU32) + mU32*pow3(m3)))/(3.*(-mQ32 +
         mU32)*pow2(m32 - mU32)) + (16*(180*m32*msq2*mU32 + 60*msq2*pow2(m32)
         + 2*mU32*pow2(m32) - 90*m32*pow2(msq2) - 30*mU32*pow2(msq2) - 3*m32*
         pow2(mU32) + pow3(mU32)))/(3.*pow3(m32 - mU32))) + (16*(180*msq2*mU32*
         pow2(m32)*pow2(mQ32) - 210*mQ32*msq2*pow2(m32)*pow2(mU32) + 150*m32*
         msq2*pow2(mQ32)*pow2(mU32) + 679*pow2(m32)*pow2(mQ32)*pow2(mU32) - 90*
         mQ32*msq2*mU32*pow3(m32) - 296*mU32*pow2(mQ32)*pow3(m32) - 1078*mQ32*
         pow2(mU32)*pow3(m32) + 90*msq2*pow2(mU32)*pow3(m32) - 90*m32*msq2*mU32
         *pow3(mQ32) + 21*mU32*pow2(m32)*pow3(mQ32) - 100*m32*pow2(mU32)*pow3(
         mQ32) - 30*msq2*pow2(mU32)*pow3(mQ32) + 54*pow3(m32)*pow3(mQ32) - 60*
         m32*mQ32*msq2*pow3(mU32) + 359*mQ32*pow2(m32)*pow3(mU32) + 30*msq2*
         pow2(m32)*pow3(mU32) - 296*m32*pow2(mQ32)*pow3(mU32) + 30*msq2*pow2(
         mQ32)*pow3(mU32) - 216*pow3(m32)*pow3(mU32) + 81*pow3(mQ32)*pow3(mU32)
         + 582*mQ32*mU32*pow4(m32) - 108*pow2(mQ32)*pow4(m32) + 550*pow2(mU32)
         *pow4(m32) - 9*m32*mU32*pow4(mQ32) - 3*pow2(mU32)*pow4(mQ32) + 158*m32
         *mQ32*pow4(mU32) - 35*pow2(m32)*pow4(mU32) - 75*pow2(mQ32)*pow4(mU32)
         + 54*mQ32*pow5(m32) - 310*mU32*pow5(m32) - 9*m32*pow5(mU32) - 3*mQ32*
         pow5(mU32)))/(3.*(-mQ32 + mU32)*pow2(m32 - mQ32)*pow3(m32 - mU32)) - (
         128*Xt*(30*m3*mQ32*msq2*mU32 + 3*m3*mU32*pow2(mQ32) + 28*m3*mQ32*pow2(
         mU32) - 18*mQ32*mU32*pow3(m3) - 30*msq2*mU32*pow3(m3) - 26*pow2(mU32)*
         pow3(m3) + 3*m3*pow3(mU32) - 16*mQ32*pow5(m3) + 10*mU32*pow5(m3) + 16*
         pow7(m3)))/(3.*(m32 - mQ32)*(mQ32 - mU32)*pow2(m32 - mU32))) + pow2(
         log(mU32/MR2))*((512*m32*pow2(mU32)*pow2(Xt))/(pow2(m32 - mU32)*pow2(
         -mQ32 + mU32)) - (8*(-120*mQ32*msq2*mU32*pow2(m32) - 33*mU32*pow2(m32)
         *pow2(mQ32) - 60*m32*mQ32*mU32*pow2(msq2) + 90*mQ32*pow2(m32)*pow2(
         msq2) - 90*mU32*pow2(m32)*pow2(msq2) + 180*m32*mQ32*msq2*pow2(mU32) -
         51*mQ32*pow2(m32)*pow2(mU32) + 120*msq2*pow2(m32)*pow2(mU32) + 36*m32*
         pow2(mQ32)*pow2(mU32) + 60*m32*pow2(msq2)*pow2(mU32) - 30*mQ32*pow2(
         msq2)*pow2(mU32) - 60*mQ32*msq2*pow3(m32) + 64*mQ32*mU32*pow3(m32) +
         60*msq2*mU32*pow3(m32) - 2*pow2(mQ32)*pow3(m32) + 706*pow2(mU32)*pow3(
         m32) - 6*m32*mU32*pow3(mQ32) + 9*pow2(m32)*pow3(mQ32) - 3*pow2(mU32)*
         pow3(mQ32) - 38*m32*mQ32*pow3(mU32) - 180*m32*msq2*pow3(mU32) - 437*
         pow2(m32)*pow3(mU32) - pow2(mQ32)*pow3(mU32) + 30*pow2(msq2)*pow3(mU32
         ) + 44*mQ32*pow4(m32) - 556*mU32*pow4(m32) + 136*m32*pow4(mU32) + 5*
         mQ32*pow4(mU32) + 128*pow5(m32) - pow5(mU32)))/(3.*(-mQ32 + mU32)*pow4
         (m32 - mU32)) + (64*Xt*(-30*m3*mQ32*mU32*pow2(msq2) + 60*m3*mQ32*msq2*
         pow2(mU32) + 10*m3*pow2(mQ32)*pow2(mU32) + 30*m3*pow2(msq2)*pow2(mU32)
         - 60*mQ32*msq2*mU32*pow3(m3) - 11*mU32*pow2(mQ32)*pow3(m3) + 30*mQ32*
         pow2(msq2)*pow3(m3) - 30*mU32*pow2(msq2)*pow3(m3) - 17*mQ32*pow2(mU32)
         *pow3(m3) + 60*msq2*pow2(mU32)*pow3(m3) - 3*m3*mU32*pow3(mQ32) + 3*
         pow3(m3)*pow3(mQ32) - 17*m3*mQ32*pow3(mU32) - 60*m3*msq2*pow3(mU32) +
         57*pow3(m3)*pow3(mU32) - 6*m3*pow4(mU32) + 68*mQ32*mU32*pow5(m3) +
         pow2(mQ32)*pow5(m3) - 85*pow2(mU32)*pow5(m3) - 18*mQ32*pow7(m3) + 18*
         mU32*pow7(m3)))/(3.*pow2(-mQ32 + mU32)*pow3(m32 - mU32))) + log(
         mQ32/MR2)*(log(msq2/MR2)*((-64*Xt*(120*m3*mQ32*msq2 - m3*pow2(mQ32) -
         60*m3*pow2(msq2) + mQ32*pow3(m3)))/(3.*(mQ32 - mU32)*pow2(m32 - mQ32))
         + (16*(180*m32*mQ32*msq2 + 2*mQ32*pow2(m32) + 60*msq2*pow2(m32) - 3*
         m32*pow2(mQ32) - 90*m32*pow2(msq2) - 30*mQ32*pow2(msq2) + pow3(mQ32)))
         /(3.*pow3(m32 - mQ32))) + (16*(-210*msq2*mU32*pow2(m32)*pow2(mQ32) +
         180*mQ32*msq2*pow2(m32)*pow2(mU32) + 150*m32*msq2*pow2(mQ32)*pow2(mU32
         ) + 676*pow2(m32)*pow2(mQ32)*pow2(mU32) - 90*mQ32*msq2*mU32*pow3(m32)
         + 90*msq2*pow2(mQ32)*pow3(m32) - 1081*mU32*pow2(mQ32)*pow3(m32) - 290*
         mQ32*pow2(mU32)*pow3(m32) - 60*m32*msq2*mU32*pow3(mQ32) + 30*msq2*pow2
         (m32)*pow3(mQ32) + 363*mU32*pow2(m32)*pow3(mQ32) - 296*m32*pow2(mU32)*
         pow3(mQ32) + 30*msq2*pow2(mU32)*pow3(mQ32) - 220*pow3(m32)*pow3(mQ32)
         - 90*m32*mQ32*msq2*pow3(mU32) + 19*mQ32*pow2(m32)*pow3(mU32) - 99*m32*
         pow2(mQ32)*pow3(mU32) - 30*msq2*pow2(mQ32)*pow3(mU32) + 55*pow3(m32)*
         pow3(mU32) + 81*pow3(mQ32)*pow3(mU32) + 580*mQ32*mU32*pow4(m32) + 555*
         pow2(mQ32)*pow4(m32) - 111*pow2(mU32)*pow4(m32) + 157*m32*mU32*pow4(
         mQ32) - 34*pow2(m32)*pow4(mQ32) - 75*pow2(mU32)*pow4(mQ32) - 9*m32*
         mQ32*pow4(mU32) - 3*pow2(mQ32)*pow4(mU32) - 312*mQ32*pow5(m32) + 56*
         mU32*pow5(m32) - 9*m32*pow5(mQ32) - 3*mU32*pow5(mQ32)))/(3.*(mQ32 -
         mU32)*pow2(m32 - mU32)*pow3(m32 - mQ32)) + (128*Xt*(30*m3*mQ32*msq2*
         mU32 + 28*m3*mU32*pow2(mQ32) + 3*m3*mQ32*pow2(mU32) - 30*mQ32*msq2*
         pow3(m3) - 18*mQ32*mU32*pow3(m3) - 26*pow2(mQ32)*pow3(m3) + 3*m3*pow3(
         mQ32) + 10*mQ32*pow5(m3) - 16*mU32*pow5(m3) + 16*pow7(m3)))/(3.*(m32 -
         mU32)*(mQ32 - mU32)*pow2(m32 - mQ32)) + log(mU32/MR2)*((-1024*m32*
         mQ32*mU32*pow2(Xt))/((m32 - mQ32)*(m32 - mU32)*pow2(mQ32 - mU32)) + (
         16*(184*pow2(m32)*pow2(mQ32)*pow2(mU32) - 194*mU32*pow2(mQ32)*pow3(m32
         ) - 87*mQ32*pow2(mU32)*pow3(m32) + 122*mU32*pow2(m32)*pow3(mQ32) - 68*
         m32*pow2(mU32)*pow3(mQ32) + 34*pow3(m32)*pow3(mQ32) + 29*mQ32*pow2(m32
         )*pow3(mU32) - 60*m32*pow2(mQ32)*pow3(mU32) - pow3(m32)*pow3(mU32) +
         20*pow3(mQ32)*pow3(mU32) + 70*mQ32*mU32*pow4(m32) - 14*pow2(mQ32)*pow4
         (m32) + 3*pow2(mU32)*pow4(m32) - 35*m32*mU32*pow4(mQ32) - 29*pow2(m32)
         *pow4(mQ32) + 4*pow2(mU32)*pow4(mQ32) + 12*mQ32*pow5(m32) - 2*mU32*
         pow5(m32) + 9*m32*pow5(mQ32) + 3*mU32*pow5(mQ32)))/(3.*pow3(m32 - mQ32
         )*pow3(m32 - mU32)) + (64*Xt*(45*pow2(mQ32)*pow2(mU32)*pow3(m3) - 37*
         m3*pow2(mU32)*pow3(mQ32) + 9*mU32*pow3(m3)*pow3(mQ32) - 9*m3*pow2(mQ32
         )*pow3(mU32) + 32*mQ32*pow3(m3)*pow3(mU32) + 20*m3*mU32*pow4(mQ32) +
         10*pow3(m3)*pow4(mQ32) - 26*mU32*pow2(mQ32)*pow5(m3) - 59*mQ32*pow2(
         mU32)*pow5(m3) - 10*pow3(mQ32)*pow5(m3) - pow3(mU32)*pow5(m3) - 6*m3*
         pow5(mQ32) + 31*mQ32*mU32*pow7(m3) + pow2(mU32)*pow7(m3)))/(3.*pow2(
         m32 - mQ32)*pow2(m32 - mU32)*pow2(mQ32 - mU32)))) + pow2(log(m32/MR2))
         *((512*pow2(Xt)*pow3(m32))/(pow2(m32 - mQ32)*pow2(m32 - mU32)) - (16*(
         -320*pow2(mU32)*pow3(m32)*pow3(mQ32) - 320*pow2(mQ32)*pow3(m32)*pow3(
         mU32) - 452*pow2(m32)*pow3(mQ32)*pow3(mU32) + 2240*pow2(mQ32)*pow2(
         mU32)*pow4(m32) + 908*mU32*pow3(mQ32)*pow4(m32) + 908*mQ32*pow3(mU32)*
         pow4(m32) - 13*pow2(m32)*pow2(mU32)*pow4(mQ32) - 134*mU32*pow3(m32)*
         pow4(mQ32) + 144*m32*pow3(mU32)*pow4(mQ32) + 27*pow4(m32)*pow4(mQ32) -
         13*pow2(m32)*pow2(mQ32)*pow4(mU32) - 134*mQ32*pow3(m32)*pow4(mU32) +
         144*m32*pow3(mQ32)*pow4(mU32) + 27*pow4(m32)*pow4(mU32) - 36*pow4(mQ32
         )*pow4(mU32) - 2700*mU32*pow2(mQ32)*pow5(m32) - 2700*mQ32*pow2(mU32)*
         pow5(m32) - 232*pow3(mQ32)*pow5(m32) - 232*pow3(mU32)*pow5(m32) + 2700
         *mQ32*mU32*pow6(m32) + 721*pow2(mQ32)*pow6(m32) + 721*pow2(mU32)*pow6(
         m32) - 726*mQ32*pow7(m32) - 726*mU32*pow7(m32) + 198*pow8(m32)))/(3.*
         pow4(m32 - mQ32)*pow4(m32 - mU32)) + (128*Xt*(pow11(m3) + 81*pow2(mQ32
         )*pow2(mU32)*pow3(m3) - 130*mU32*pow2(mQ32)*pow5(m3) - 130*mQ32*pow2(
         mU32)*pow5(m3) + 196*mQ32*mU32*pow7(m3) + 41*pow2(mQ32)*pow7(m3) + 41*
         pow2(mU32)*pow7(m3) - 50*mQ32*pow9(m3) - 50*mU32*pow9(m3)))/(3.*pow3(
         m32 - mQ32)*pow3(m32 - mU32))) + pow2(log(mQ32/MR2))*((512*m32*pow2(
         mQ32)*pow2(Xt))/(pow2(m32 - mQ32)*pow2(mQ32 - mU32)) - (8*(270*pow2(
         m32)*pow2(mQ32)*pow2(msq2)*pow2(mU32) - 210*mU32*pow2(mQ32)*pow2(msq2)
         *pow3(m32) - 180*msq2*pow2(mQ32)*pow2(mU32)*pow3(m32) - 90*mQ32*pow2(
         msq2)*pow2(mU32)*pow3(m32) - 90*mU32*pow2(m32)*pow2(msq2)*pow3(mQ32) -
         540*msq2*pow2(m32)*pow2(mU32)*pow3(mQ32) + 90*m32*pow2(msq2)*pow2(
         mU32)*pow3(mQ32) + 540*msq2*mU32*pow3(m32)*pow3(mQ32) + 30*pow2(msq2)*
         pow3(m32)*pow3(mQ32) - 1518*pow2(mU32)*pow3(m32)*pow3(mQ32) + 420*msq2
         *pow2(m32)*pow2(mQ32)*pow3(mU32) - 90*mQ32*pow2(m32)*pow2(msq2)*pow3(
         mU32) - 150*m32*pow2(mQ32)*pow2(msq2)*pow3(mU32) - 420*mQ32*msq2*pow3(
         m32)*pow3(mU32) - 894*pow2(mQ32)*pow3(m32)*pow3(mU32) + 270*pow2(msq2)
         *pow3(m32)*pow3(mU32) + 180*m32*msq2*pow3(mQ32)*pow3(mU32) + 422*pow2(
         m32)*pow3(mQ32)*pow3(mU32) - 30*pow2(msq2)*pow3(mQ32)*pow3(mU32) - 180
         *msq2*mU32*pow2(mQ32)*pow4(m32) + 210*mQ32*mU32*pow2(msq2)*pow4(m32) +
         60*pow2(mQ32)*pow2(msq2)*pow4(m32) + 540*mQ32*msq2*pow2(mU32)*pow4(
         m32) + 2495*pow2(mQ32)*pow2(mU32)*pow4(m32) - 270*pow2(msq2)*pow2(mU32
         )*pow4(m32) - 180*msq2*pow3(mQ32)*pow4(m32) + 1137*mU32*pow3(mQ32)*
         pow4(m32) + 699*mQ32*pow3(mU32)*pow4(m32) - 180*msq2*pow3(mU32)*pow4(
         m32) + 532*pow2(m32)*pow2(mU32)*pow4(mQ32) - 158*mU32*pow3(m32)*pow4(
         mQ32) - 143*m32*pow3(mU32)*pow4(mQ32) + 192*pow4(m32)*pow4(mQ32) + 120
         *mQ32*msq2*pow2(m32)*pow4(mU32) - 180*m32*msq2*pow2(mQ32)*pow4(mU32) +
         62*pow2(m32)*pow2(mQ32)*pow4(mU32) + 60*m32*mQ32*pow2(msq2)*pow4(mU32
         ) - 90*pow2(m32)*pow2(msq2)*pow4(mU32) + 30*pow2(mQ32)*pow2(msq2)*pow4
         (mU32) - 52*mQ32*pow3(m32)*pow4(mU32) + 60*msq2*pow3(m32)*pow4(mU32) +
         12*m32*pow3(mQ32)*pow4(mU32) - 43*pow4(m32)*pow4(mU32) - 3*pow4(mQ32)
         *pow4(mU32) - 300*mQ32*msq2*mU32*pow5(m32) + 120*msq2*pow2(mQ32)*pow5(
         m32) - 2133*mU32*pow2(mQ32)*pow5(m32) - 90*mQ32*pow2(msq2)*pow5(m32) +
         90*mU32*pow2(msq2)*pow5(m32) - 1879*mQ32*pow2(mU32)*pow5(m32) + 180*
         msq2*pow2(mU32)*pow5(m32) - 473*pow3(mQ32)*pow5(m32) + 5*pow3(mU32)*
         pow5(m32) - 158*mU32*pow2(m32)*pow5(mQ32) - 35*m32*pow2(mU32)*pow5(
         mQ32) - 66*pow3(m32)*pow5(mQ32) + 7*pow3(mU32)*pow5(mQ32) + 60*mQ32*
         msq2*pow6(m32) + 1704*mQ32*mU32*pow6(m32) - 60*msq2*mU32*pow6(m32) +
         722*pow2(mQ32)*pow6(m32) + 262*pow2(mU32)*pow6(m32) + 47*m32*mU32*pow6
         (mQ32) + 38*pow2(m32)*pow6(mQ32) - pow2(mU32)*pow6(mQ32) - 556*mQ32*
         pow7(m32) - 340*mU32*pow7(m32) - 9*m32*pow7(mQ32) - 3*mU32*pow7(mQ32)
         + 128*pow8(m32)))/(3.*(mQ32 - mU32)*pow3(m32 - mU32)*pow4(m32 - mQ32))
         + (64*Xt*(18*mQ32*pow11(m3) - 18*mU32*pow11(m3) + 30*m3*pow2(mQ32)*
         pow2(msq2)*pow2(mU32) - 60*mU32*pow2(mQ32)*pow2(msq2)*pow3(m3) - 60*
         msq2*pow2(mQ32)*pow2(mU32)*pow3(m3) + 30*mQ32*pow2(msq2)*pow2(mU32)*
         pow3(m3) - 60*m3*msq2*pow2(mU32)*pow3(mQ32) + 120*msq2*mU32*pow3(m3)*
         pow3(mQ32) + 106*pow2(mU32)*pow3(m3)*pow3(mQ32) + 60*m3*msq2*pow2(mQ32
         )*pow3(mU32) - 30*m3*mQ32*pow2(msq2)*pow3(mU32) - 60*mQ32*msq2*pow3(m3
         )*pow3(mU32) - 40*pow2(mQ32)*pow3(m3)*pow3(mU32) + 30*pow2(msq2)*pow3(
         m3)*pow3(mU32) - 3*m3*pow3(mQ32)*pow3(mU32) - 20*m3*pow2(mU32)*pow4(
         mQ32) - 10*mU32*pow3(m3)*pow4(mQ32) - 60*msq2*mU32*pow2(mQ32)*pow5(m3)
         + 30*mQ32*mU32*pow2(msq2)*pow5(m3) + 30*pow2(mQ32)*pow2(msq2)*pow5(m3
         ) + 120*mQ32*msq2*pow2(mU32)*pow5(m3) - 59*pow2(mQ32)*pow2(mU32)*pow5(
         m3) - 60*pow2(msq2)*pow2(mU32)*pow5(m3) - 60*msq2*pow3(mQ32)*pow5(m3)
         - 103*mU32*pow3(mQ32)*pow5(m3) + 81*mQ32*pow3(mU32)*pow5(m3) - 15*pow4
         (mQ32)*pow5(m3) + 10*m3*mU32*pow5(mQ32) + 8*pow3(m3)*pow5(mQ32) - 3*m3
         *pow6(mQ32) - 60*mQ32*msq2*mU32*pow7(m3) + 60*msq2*pow2(mQ32)*pow7(m3)
         + 138*mU32*pow2(mQ32)*pow7(m3) - 30*mQ32*pow2(msq2)*pow7(m3) + 30*
         mU32*pow2(msq2)*pow7(m3) - 112*mQ32*pow2(mU32)*pow7(m3) + 60*pow3(mQ32
         )*pow7(m3) - 22*pow3(mU32)*pow7(m3) + 31*mQ32*mU32*pow9(m3) - 84*pow2(
         mQ32)*pow9(m3) + 37*pow2(mU32)*pow9(m3)))/(3.*pow2(m32 - mU32)*pow2(
         mQ32 - mU32)*pow3(m32 - mQ32))) + log(m32/MR2)*((16*(-90*msq2*mU32*
         pow2(m32)*pow2(mQ32) - 90*mQ32*msq2*pow2(m32)*pow2(mU32) + 1822*pow2(
         m32)*pow2(mQ32)*pow2(mU32) + 180*mQ32*msq2*mU32*pow3(m32) - 270*msq2*
         pow2(mQ32)*pow3(m32) - 2575*mU32*pow2(mQ32)*pow3(m32) - 2575*mQ32*pow2
         (mU32)*pow3(m32) - 270*msq2*pow2(mU32)*pow3(m32) + 30*m32*msq2*mU32*
         pow3(mQ32) + 90*msq2*pow2(m32)*pow3(mQ32) + 498*mU32*pow2(m32)*pow3(
         mQ32) - 333*m32*pow2(mU32)*pow3(mQ32) - 269*pow3(m32)*pow3(mQ32) + 30*
         m32*mQ32*msq2*pow3(mU32) + 498*mQ32*pow2(m32)*pow3(mU32) + 90*msq2*
         pow2(m32)*pow3(mU32) - 333*m32*pow2(mQ32)*pow3(mU32) - 269*pow3(m32)*
         pow3(mU32) + 48*pow3(mQ32)*pow3(mU32) + 240*mQ32*msq2*pow4(m32) + 3564
         *mQ32*mU32*pow4(m32) + 240*msq2*mU32*pow4(m32) + 1182*pow2(mQ32)*pow4(
         m32) + 1182*pow2(mU32)*pow4(m32) + 3*m32*mU32*pow4(mQ32) + 9*pow2(m32)
         *pow4(mQ32) + 3*m32*mQ32*pow4(mU32) + 9*pow2(m32)*pow4(mU32) - 1562*
         mQ32*pow5(m32) - 180*msq2*pow5(m32) - 1562*mU32*pow5(m32) + 660*pow6(
         m32)))/(3.*pow3(m32 - mQ32)*pow3(m32 - mU32)) - (128*Xt*(30*mQ32*msq2*
         pow3(m3) + 34*mQ32*mU32*pow3(m3) + 30*msq2*mU32*pow3(m3) + 3*pow2(mQ32
         )*pow3(m3) + 3*pow2(mU32)*pow3(m3) - 40*mQ32*pow5(m3) - 60*msq2*pow5(
         m3) - 40*mU32*pow5(m3) + 40*pow7(m3)))/(3.*pow2(m32 - mQ32)*pow2(m32 -
         mU32)) + log(msq2/MR2)*((16*(-447*mU32*pow2(mQ32)*pow3(m32) - 447*
         mQ32*pow2(mU32)*pow3(m32) + 149*mU32*pow2(m32)*pow3(mQ32) - 29*pow3(
         m32)*pow3(mQ32) + 149*mQ32*pow2(m32)*pow3(mU32) - 29*pow3(m32)*pow3(
         mU32) + 894*mQ32*mU32*pow4(m32) + 87*pow2(mQ32)*pow4(m32) + 87*pow2(
         mU32)*pow4(m32) - 236*mQ32*pow5(m32) - 236*mU32*pow5(m32) + 58*pow6(
         m32)))/(3.*pow3(m32 - mQ32)*pow3(m32 - mU32)) + (64*Xt*(-119*mQ32*mU32
         *pow3(m3) + 59*mQ32*pow5(m3) + 59*mU32*pow5(m3) + pow7(m3)))/(3.*pow2(
         m32 - mQ32)*pow2(m32 - mU32))) + log(mQ32/MR2)*((-1024*mQ32*pow2(m32)*
         pow2(Xt))/((m32 - mU32)*(mQ32 - mU32)*pow2(m32 - mQ32)) + (64*(-3*pow2
         (mU32)*pow3(m32)*pow3(mQ32) - 164*pow2(mQ32)*pow3(m32)*pow3(mU32) - 11
         *pow2(m32)*pow3(mQ32)*pow3(mU32) + 290*pow2(mQ32)*pow2(mU32)*pow4(m32)
         - 49*mU32*pow3(mQ32)*pow4(m32) + 251*mQ32*pow3(mU32)*pow4(m32) + 9*
         pow2(m32)*pow2(mU32)*pow4(mQ32) + 25*mU32*pow3(m32)*pow4(mQ32) - 19*
         pow4(m32)*pow4(mQ32) + 42*pow2(m32)*pow2(mQ32)*pow4(mU32) - 55*mQ32*
         pow3(m32)*pow4(mU32) + 7*pow4(m32)*pow4(mU32) - 170*mU32*pow2(mQ32)*
         pow5(m32) - 451*mQ32*pow2(mU32)*pow5(m32) + 33*pow3(mQ32)*pow5(m32) -
         52*pow3(mU32)*pow5(m32) - 8*mU32*pow2(m32)*pow5(mQ32) + 5*pow3(m32)*
         pow5(mQ32) + 329*mQ32*mU32*pow6(m32) + 32*pow2(mQ32)*pow6(m32) + 119*
         pow2(mU32)*pow6(m32) - 89*mQ32*pow7(m32) - 103*mU32*pow7(m32) + 32*
         pow8(m32)))/(3.*(mQ32 - mU32)*pow3(m32 - mU32)*pow4(m32 - mQ32)) - (64
         *Xt*(68*pow11(m3) + 151*pow2(mQ32)*pow2(mU32)*pow3(m3) - 15*mU32*pow3(
         m3)*pow3(mQ32) - 241*mU32*pow2(mQ32)*pow5(m3) - 260*mQ32*pow2(mU32)*
         pow5(m3) - 7*pow3(mQ32)*pow5(m3) + 459*mQ32*mU32*pow7(m3) + 140*pow2(
         mQ32)*pow7(m3) + 77*pow2(mU32)*pow7(m3) - 233*mQ32*pow9(m3) - 139*mU32
         *pow9(m3)))/(3.*(mQ32 - mU32)*pow2(m32 - mU32)*pow3(m32 - mQ32))) +
         log(mU32/MR2)*((1024*mU32*pow2(m32)*pow2(Xt))/((m32 - mQ32)*(mQ32 -
         mU32)*pow2(m32 - mU32)) - (16*(-655*pow2(mU32)*pow3(m32)*pow3(mQ32) -
         11*pow2(mQ32)*pow3(m32)*pow3(mU32) - 43*pow2(m32)*pow3(mQ32)*pow3(mU32
         ) + 1157*pow2(mQ32)*pow2(mU32)*pow4(m32) + 999*mU32*pow3(mQ32)*pow4(
         m32) - 191*mQ32*pow3(mU32)*pow4(m32) + 167*pow2(m32)*pow2(mU32)*pow4(
         mQ32) - 218*mU32*pow3(m32)*pow4(mQ32) + 27*pow4(m32)*pow4(mQ32) + 35*
         pow2(m32)*pow2(mQ32)*pow4(mU32) + 97*mQ32*pow3(m32)*pow4(mU32) - 72*
         pow4(m32)*pow4(mU32) - 1797*mU32*pow2(mQ32)*pow5(m32) - 683*mQ32*pow2(
         mU32)*pow5(m32) - 205*pow3(mQ32)*pow5(m32) + 125*pow3(mU32)*pow5(m32)
         - 31*mQ32*pow2(m32)*pow5(mU32) + 19*pow3(m32)*pow5(mU32) + 1314*mQ32*
         mU32*pow6(m32) + 472*pow2(mQ32)*pow6(m32) + 134*pow2(mU32)*pow6(m32) -
         410*mQ32*pow7(m32) - 358*mU32*pow7(m32) + 128*pow8(m32)))/(3.*(mQ32 -
         mU32)*pow3(m32 - mQ32)*pow4(m32 - mU32)) + (128*Xt*(34*pow11(m3) + 75
         *pow2(mQ32)*pow2(mU32)*pow3(m3) - 7*mQ32*pow3(m3)*pow3(mU32) - 129*
         mU32*pow2(mQ32)*pow5(m3) - 121*mQ32*pow2(mU32)*pow5(m3) - 4*pow3(mU32)
         *pow5(m3) + 229*mQ32*mU32*pow7(m3) + 38*pow2(mQ32)*pow7(m3) + 71*pow2(
         mU32)*pow7(m3) - 69*mQ32*pow9(m3) - 117*mU32*pow9(m3)))/(3.*(mQ32 -
         mU32)*pow2(m32 - mQ32)*pow3(m32 - mU32))));

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

/// 1-loop correction O(at), at = yt^2 Sin[beta]^2/(4 Pi)
double himalaya::mh2_eft::Mh2EFTCalculator::Mh2_EFT_1loop(
   double at, double mt, double mQ32, double mU32,
   double Xt, double MR2)
{
   const double mt2 = pow2(mt);
   const double L = std::log(MR2/mt2);
   const double oneLoop = 1./(4*Pi);

   return oneLoop * mt2 * at * (
      coeff_as_0_log_0(mQ32, mU32, Xt, MR2) +
      coeff_as_0_log_1() * L
   );
}

/// 2-loop correction O(at*as)
/// at = yt^2 Sin[beta]^2/(4 Pi), as = g3^2/(4 Pi), 
double himalaya::mh2_eft::Mh2EFTCalculator::Mh2_EFT_2loop(
   double at, double mt, double mQ32, double mU32,
   double Xt, double MR2, double g3, double m3)
{
   const double mt2 = pow2(mt);
   const double L = std::log(MR2/mt2);
   const double as = pow2(g3)/(4*Pi);
   const double twoLoop = 1./pow2(4*Pi);

   return twoLoop * mt2 * at * as * (
      coeff_as_1_log_0(mQ32, mU32, Xt, m3, MR2) +
      coeff_as_1_log_1(mQ32, mU32, Xt, m3, MR2) * L +
      coeff_as_1_log_2() * pow2(L)
   );
}

/// 3-loop correction O(at*as^2), without \[CapitalDelta]\[Lambda]^(3L),
/// at = yt^2 Sin[beta]^2/(4 Pi), as = g3^2/(4 Pi), 
double himalaya::mh2_eft::Mh2EFTCalculator::Mh2_EFT_3loop(
   double at, double mt, double mQ32, double mU32,
   double Xt, double MR2, double g3, double m3, double msq2)
{
   const double mt2 = pow2(mt);
   const double L = std::log(MR2/mt2);
   const double as = pow2(g3)/(4*Pi);
   const double threeLoop = 1./pow3(4*Pi);

   return threeLoop * mt2 * at * pow2(as) * (
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
   checkPassed = dum - 65.75234703 < 1e-03;
   std::cout << "Check as^0 log^0 passed: " << tf(checkPassed) << ".\n";
   dum = coeff_as_0_log_1();
   checkPassed = dum - 12. < 1e-03;
   std::cout << "Check as^0 log^1 passed: " << tf(checkPassed) << ".\n";
   dum = coeff_as_1_log_0(mQ32, mU32, Xt, m3, MR2);
   checkPassed = dum - 977.1683309 < 1e-03;
   std::cout << "Check as^1 log^0 passed: " << tf(checkPassed) << ".\n";
   dum = coeff_as_1_log_1(mQ32, mU32, Xt, m3, MR2);
   checkPassed = dum - 186.2109432 < 1e-03;
   std::cout << "Check as^1 log^1 passed: " << tf(checkPassed) << ".\n";
   dum = coeff_as_1_log_2();
   checkPassed = dum - 96. < 1e-03;
   std::cout << "Check as^1 log^2 passed: " << tf(checkPassed) << ".\n";
   dum = coeff_as_2_log_0(mQ32, mU32, Xt, m3, msq2, MR2);
   checkPassed = dum - (-13930.41545) < 1e-03;
   std::cout << "Check as^2 log^0 passed: " << tf(checkPassed) << ".\n";
   dum = coeff_as_2_log_1(mQ32, mU32, Xt, m3, msq2, MR2);
   checkPassed = dum - (-5580.985505) < 1e-03;
   std::cout << "Check as^2 log^1 passed: " << tf(checkPassed) << ".\n";
   dum = coeff_as_2_log_2(mQ32, mU32, Xt, m3, msq2, MR2);
   checkPassed = dum - 2929.488094 < 1e-03;
   std::cout << "Check as^2 log^2 passed: " << tf(checkPassed) << ".\n";
   dum = coeff_as_2_log_3();
   checkPassed = dum - 736. < 1e-03;
   std::cout << "Check as^2 log^3 passed: " << tf(checkPassed) << ".\n";
}
