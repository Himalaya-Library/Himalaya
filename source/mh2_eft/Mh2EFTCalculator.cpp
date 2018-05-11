#include "Mh2EFTCalculator.hpp"
#include "ThresholdCalculator.hpp"
#include "Hierarchies.hpp"
#include <iostream>

namespace himalaya {
namespace mh2_eft {
namespace {
   const double zt3 = 1.2020569031595942853997381615114;
   const double Pi  = 3.1415926535897932384626433832795;
   const double log2 = std::log(2.);
   template <typename T> T pow2(T x)  { return x*x; }
   template <typename T> T pow3(T x)  { return x*x*x; }
   template <typename T> T pow4(T x)  { return x*x*x*x; }
   template <typename T> T pow5(T x)  { return x*x*x*x*x; }
   template <typename T> T pow6(T x)  { return x*x*x*x*x*x; }
   template <typename T> T pow7(T x)  { return x*x*x*x*x*x*x; }
   template <typename T> T pow8(T x)  { return x*x*x*x*x*x*x*x; }
   template <typename T> T pow9(T x)  { return x*x*x*x*x*x*x*x*x; }
   template <typename T> T power10(T x) { return pow9(x)*x; }
   template <typename T> T pow11(T x) { return pow2(x)*pow9(x); }
   template <typename T> T pow12(T x) { return pow2(pow6(x)); }
   template <typename T> T pow13(T x) { return pow4(x)*pow9(x); }
   template <typename T> T pow14(T x) { return pow2(pow7(x)); }
   template <typename T> T pow15(T x) { return pow6(x)*pow9(x); }
   template <typename T> T pow16(T x) { return pow2(pow8(x)); }
   template <typename T> T pow17(T x) { return pow16(x)*x; }
   template <typename T> T pow18(T x) { return pow2(pow9(x)); }
   template <typename T> T pow19(T x) { return pow18(x)*x; }
   template <typename T> T pow20(T x) { return pow19(x)*x; }

} // anonymous namespace
} // namespace mh2_eft
} // namespace himalaya

/**
 *	Constructor
 * 	@param p a HimalayaInterface struct
 * 	@param msq2 the averaged squark mass of the first two generations squared
 * 	@param verbose a bool enable the output of the parameter validation. Enabled by default
 */
himalaya::mh2_eft::Mh2EFTCalculator::Mh2EFTCalculator(
   const himalaya::Parameters& p_, double msq2_, bool verbose)
   : p(p_), msq2(msq2_)
{
   p.validate(verbose);

   if (!std::isfinite(msq2_))
      msq2 = p.calculateMsq2();
}

/**
 * 	Returns the 1-loop EFT contribution to the Higgs mass
 * 	@param omitSMLogs an integer flag to remove all Log(mu^2/mt^2) terms
 * 	@param omitMSSMLogs an integer flag to remove all Log(mu^2/Mx^2) terms
 */
double himalaya::mh2_eft::Mh2EFTCalculator::getDeltaMh2EFT1Loop(int omitSMLogs, int omitMSSMLogs){
   ThresholdCalculator thresholdCalculator(p, msq2);
   
   using std::log;
   const double lmMt = omitSMLogs * log(pow2(p.scale / p.Mt));
   
   return 12 * lmMt + 
      thresholdCalculator.getThresholdCorrection(ThresholdVariables::LAMBDA_AT,
	 RenSchemes::DRBARPRIME, omitMSSMLogs);
}

/**
 * 	Returns the 2-loop EFT contribution to the Higgs mass
 * 	@param omitSMLogs an integer flag to remove all Log(mu^2/mt^2) terms
 * 	@param omitMSSMLogs an integer flag to remove all Log(mu^2/Mx^2) terms
 */
double himalaya::mh2_eft::Mh2EFTCalculator::getDeltaMh2EFT2Loop(int omitSMLogs, int omitMSSMLogs){
   ThresholdCalculator thresholdCalculator(p, msq2);
   
   using std::log;
   const double lmMt = omitSMLogs * log(pow2(p.scale / p.Mt));
   const double dytas = thresholdCalculator.getThresholdCorrection(
      ThresholdVariables::YT_AS, RenSchemes::DRBARPRIME, omitMSSMLogs);

   return 96 * pow2(lmMt) + (-32 + 48 * dytas) * lmMt - 24 * dytas
      + thresholdCalculator.getThresholdCorrection(ThresholdVariables::LAMBDA_AT_AS,
	 RenSchemes::DRBARPRIME, omitMSSMLogs);
}

/**
 * 	Returns the 3-loop EFT contribution to the Higgs mass
 * 	@param omitSMLogs an integer flag to remove all Log(mu^2/mt^2) terms
 * 	@param omitMSSMLogs an integer flag to remove all Log(mu^2/Mx^2) terms
 * 	@param xtOrder an integer to subtract all Xt contributions starting at Xt^4 up to Xt^6 (e.g. xtOrder = 4 will subtract O(Xt^5 + Xt^6) terms). By default no Xt terms are subtracted. SM logs are set to 0 by default.
 */
double himalaya::mh2_eft::Mh2EFTCalculator::getDeltaMh2EFT3Loop(int omitSMLogs, int omitMSSMLogs, int xtOrder){
   ThresholdCalculator thresholdCalculator(p, msq2);
   
   using std::log;
   
   const double catas2 = 248.1215180432007;
   const double lmMt = omitSMLogs * log(pow2(p.scale / p.Mt));
   // threshold correction of yt_as DRbar'
   const double dytas = thresholdCalculator.getThresholdCorrection(
      ThresholdVariables::YT_AS, RenSchemes::DRBARPRIME, omitMSSMLogs);
   // threshold correction of yt_as2 DRbar'
   const double dytas2 = thresholdCalculator.getThresholdCorrection(
      ThresholdVariables::YT_AS2, RenSchemes::DRBARPRIME, omitMSSMLogs);
   // threshold correction of g3_as DRbar'
   const double dg3as = thresholdCalculator.getThresholdCorrection(
      ThresholdVariables::G3_AS, RenSchemes::DRBARPRIME, omitMSSMLogs);
   double xtSubtraction;
   switch(xtOrder){
      case 5:
	 xtSubtraction = thresholdCalculator.getXtTerms(6, omitMSSMLogs);
	 break;
      case 4:
	 xtSubtraction = thresholdCalculator.getXtTerms(6, omitMSSMLogs)
	    + thresholdCalculator.getXtTerms(5, omitMSSMLogs);
	 break;
      case 3:
	 xtSubtraction = thresholdCalculator.getXtTerms(6, omitMSSMLogs)
	    + thresholdCalculator.getXtTerms(5, omitMSSMLogs)
	    + thresholdCalculator.getXtTerms(4, omitMSSMLogs);
	 break;
      default:
	 xtSubtraction = 0;
	 break;
   }

   return 736 * pow3(lmMt) + (160 + 192 * dg3as + 384 * dytas) * pow2(lmMt)
      + (-128 * zt3 - 2056 / 3. + -64 * dg3as - 512 * dytas + 72 * pow2(dytas)
      + 48 * dytas2) * lmMt + 64 * dytas - 84 * pow2(dytas) - 24 * dytas2
      + catas2
      + thresholdCalculator.getThresholdCorrection(ThresholdVariables::LAMBDA_AT_AS2,
	 RenSchemes::DRBARPRIME, omitMSSMLogs) - xtSubtraction;
}

/**
 *   Returns the matching relation of zeta_lambda^(2) for the degenerate mass case
 *   @param scale the renormalization scale
 *   @param mst1 the mass of the light stop quark
 *   @return zeta_lambda^(2)
 */
double himalaya::mh2_eft::Mh2EFTCalculator::getZetaDegenerate(double scale, double mst1, double Xt, int omitlogs) const{
   using std::log;
   
   const double li4 = 0.5174790616738993; // Li4(1/2)
   
   const double LS = omitlogs*log(pow2(scale / p.MSt(0))) + log(pow2(p.MSt(0) / mst1));
   
   const double zt2lam = (29365 + 23040*li4 - 49320*zt3 + 160*LS*(-226 + 27*zt3) + 26520*pow2(LS) -
        80*Xt*(823 - 100*LS - 477*zt3 + 366*pow2(LS)) + 30*(-67 - 220*LS - 990*
        zt3 + 84*pow2(LS))*pow2(Xt) - 960*pow2(Pi)*pow2(log2) - 10080*pow3(LS)
        + 20*(3568 + 20*LS - 2259*zt3 + 108*pow2(LS))*pow3(Xt) - 176*pow4(Pi) +
        960*pow4(log2))/540.;
   return zt2lam;
}

/**
 * 	A function which maps a boolean to a string.
 * 	@param tf a boolean.
 * 	@return A string which is 'true' if tf is true or 'false' if tf is false.
 */
std::string himalaya::mh2_eft::Mh2EFTCalculator::tf(const bool tf){
   return tf ? "\033[1;32mtrue\033[0m" : "\033[1;31mfalse\033[0m";
}
