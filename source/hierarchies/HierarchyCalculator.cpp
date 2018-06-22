// ====================================================================
// This file is part of Himalaya.
//
// Himalaya is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================

#include "HierarchyCalculator.hpp"
#include "Mh2EFTCalculator.hpp"
#include "H3.hpp"
#include "H32q2g.hpp"
#include "H3q22g.hpp"
#include "H4.hpp"
#include "H5.hpp"
#include "H5g1.hpp"
#include "H6.hpp"
#include "H6b.hpp"
#include "H6b2qg2.hpp"
#include "H6bq22g.hpp"
#include "H6bq2g2.hpp"
#include "H6g2.hpp"
#include "H9.hpp"
#include "H9q2.hpp"
#include "Constants.hpp"
#include "Logger.hpp"
#include "Utils.hpp"
#include "ThresholdCalculator.hpp"
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <type_traits>

extern "C" void dszhiggs_(double *t, double *mg, double *T1, double *T2, double *st, double *ct, double *q, double *mu, double *tanb,
      double *v2, double *gs, int *OS, double *S11, double *S22, double *S12);

static bool isInfoPrinted; /**< If this bool is true, than no info will be printed in further runs */

/**
 * 	Define static variables
 */
namespace himalaya {
namespace {

/** The hierarchy map which maps all hierarchies to their mother hierarchies */
const std::map<int, int> hierarchyMap = {
   { Hierarchies::h3     , Hierarchies::h3  },
   { Hierarchies::h32q2g , Hierarchies::h3  },
   { Hierarchies::h3q22g , Hierarchies::h3  },
   { Hierarchies::h4     , Hierarchies::h4  },
   { Hierarchies::h5     , Hierarchies::h5  },
   { Hierarchies::h5g1   , Hierarchies::h5  },
   { Hierarchies::h6     , Hierarchies::h6  },
   { Hierarchies::h6g2   , Hierarchies::h6  },
   { Hierarchies::h6b    , Hierarchies::h6b },
   { Hierarchies::h6b2qg2, Hierarchies::h6b },
   { Hierarchies::h6bq22g, Hierarchies::h6b },
   { Hierarchies::h6bq2g2, Hierarchies::h6b },
   { Hierarchies::h9     , Hierarchies::h9  },
   { Hierarchies::h9q2   , Hierarchies::h9  }
};

} // anonymous namespace
} // namespace himalaya

/**
 * 	Constructor 
 * 	@param p_ a HimalayaInterface struct
 * 	@param verbose_ suppress informative output during the calculation, if set to false
 */
himalaya::HierarchyCalculator::HierarchyCalculator(const Parameters& p_,
						   const bool verbose_)
   : p(p_)
   , verbose(verbose_)
{
   if(!isInfoPrinted && verbose){
      printInfo();
      isInfoPrinted = true;
   }

   p.validate(verbose);

   // init common variables
   init();
}

/**
 * 	Initializes all common variables.
 */
void himalaya::HierarchyCalculator::init(){
   // fill flag list
   flagMap.clear();
   for (int i = ExpansionDepth::FIRST; i < ExpansionDepth::NUMBER_OF_EXPANSIONS; i++) {
      flagMap.emplace(i, 1u);
   }

   // beta
   const double beta = atan(p.vu / p.vd);

   // Al4p
   Al4p = pow2(p.g3 / (4 * Pi));

   // MGl
   Mgl = p.MG;

   // Msq, checked
   Msq = std::sqrt(std::abs(p.calculateMsq2()));

   // lmMsq, checked
   lmMsq = log(pow2(p.scale / Msq));

   // lmMgl, checked
   lmMgl = log(pow2(p.scale / Mgl));

   // prefactor, GF = 1/(sqrt(2) * (vu^2 + vd^2)) (here, GF is calculated in the DRbar scheme, checked)
   prefac = (3. / (sqrt(2) * (pow2(p.vu) + pow2(p.vd)) * sqrt(2) * pow2(Pi) * pow2(sin(beta))));
}

/**
 * 	Calculates the 3-loop mass matrix and other information of the hierarchy selection process.
 * 	@param isAlphab a bool which determines if the returned object is proportinal to alpha_b.
 * 	@return A HierarchyObject which holds all information of the calculation.
 */
himalaya::HierarchyObject himalaya::HierarchyCalculator::calculateDMh3L(bool isAlphab){
   HierarchyObject ho (isAlphab);
   
   if (isAlphab)
      INFO_MSG("3-loop threshold correction Δλ not available for O(ab*as^2)!");
   
   const int mdrFlag = 0;
   
   // set mdrFlag
   ho.setMDRFlag(mdrFlag);

   // compare hierarchies and get the best fitting hierarchy
   compareHierarchies(ho);
   
   // set Xt order truncation for EFT contribution to be consistent with H3m
   int xtOrder = 4;
   const int suitableHierarchy = ho.getSuitableHierarchy();
   if(suitableHierarchy == himalaya::Hierarchies::h3
      || suitableHierarchy == himalaya::Hierarchies::h32q2g 
      || suitableHierarchy == himalaya::Hierarchies::h3q22g
      || suitableHierarchy == himalaya::Hierarchies::h9
      || suitableHierarchy == himalaya::Hierarchies::h9q2) xtOrder = 3;
   
   // calculate the 3-loop Higgs mass matrix for the obtained hierachy in the (M)DRbar' scheme
   ho.setDMh(3, calculateHierarchy(ho, 0, 0, 1) + shiftH3mToDRbarPrime(ho));
   
   // set the alpha_x contributions
   ho.setDMh(1, getMt41L(ho, mdrFlag, mdrFlag));
   
   // set the alpha_x*alpha_s contributions
   ho.setDMh(2, getMt42L(ho, mdrFlag, mdrFlag));
   
   // estimate the uncertainty of the expansion at 3-loop level
   ho.setDMhExpUncertainty(3, getExpansionUncertainty(ho,
						   ho.getDMh(0) + ho.getDMh(1) + ho.getDMh(2), 0, 0, 1));
   
   // set the uncertainty of the expansion at 1-loop level to 0 by default, 
   // if the user needs this value getExpansionUncertainty should be called
   ho.setDMhExpUncertainty(1, 0.);

   // calculate delta_lambda
   // create a modified parameters struct and construct Mh2EFTCalculator and ThresholdCalculator
   auto p_mass_ES = p;
   p_mass_ES.mu2(2,2) = pow2(p.MSt(0));
   p_mass_ES.mq2(2,2) = pow2(p.MSt(1));
   himalaya::mh2_eft::Mh2EFTCalculator mh2EFTCalculator(p_mass_ES);
   himalaya::ThresholdCalculator tc (p_mass_ES);
   
   // to obtain delta_lambda one has to divide the difference of the two calculations by v^2
   const double v2 = pow2(p.vu) + pow2(p.vd);
   
   const double gt = sqrt(2)*p.Mt/std::sqrt(v2);
   
   const double pref = 1./pow6(4*Pi) * pow2(p.Mt * gt * pow2(p.g3));
   
   // calculate the (non-)logarithmic part of Mh2 without delta_lambda_3L
   // the first line is equivalent to  64 * dytas - 84 * pow2(dytas) - 24 * dytas2 + catas2 including all log(mu^2/Mst1^2)
   // which are also included in the H3m result, and the second line omits these log(mu^2/Mst1^2) terms since they originate
   // from SM contributions. The non-logarithmic part contains only Xt orders up to O(Xt^2) which are also included in H3m
   // so no subtraction is needed. Checked.
   const double subtractionTermH3m = mh2EFTCalculator.getDeltaMh2EFT3Loop(0,1,0);
   const double subtractionTermEFT = mh2EFTCalculator.getDeltaMh2EFT3Loop(0,0,0);

   // calculate the EFT logs. In the first call we calculate the full reconstructed contribution to delta_lambda_3L
   // including all logarithmic contributions. In the second line we subtract all non-logarithmic contributions
   // to isolate the logarithmic ones. Checked.
   const double eftLogs = pref*(
      tc.getThresholdCorrection(
	    ThresholdVariables::LAMBDA_AT_AS2, RenSchemes::DRBARPRIME, 1)
      - tc.getThresholdCorrection(
	    ThresholdVariables::LAMBDA_AT_AS2, RenSchemes::DRBARPRIME, 0));

   // calculate the non-logarithmic part of delta_lambda at 3L
   const double deltaLambda3LNonLog = pref*(ho.getDLambdaNonLog()
      + shiftH3mToDRbarPrimeMh2(ho,0)) - subtractionTermEFT;

   // calculate delta_lambda_H3m
   ho.setDLambdaH3m((pref*(ho.getDLambdaH3m()
      + shiftH3mToDRbarPrimeMh2(ho,1)) - subtractionTermH3m)/v2);
   
   // caluclate delta_lambda_EFT
   ho.setDLambdaEFT((deltaLambda3LNonLog + eftLogs)/v2);

   // save the non-logarithmic part of delta_lambda 3L
   ho.setDLambdaNonLog(deltaLambda3LNonLog/v2);
   
   // calculate DR' -> MS shift for delta_lambda 3L
   ho.setDLambdaH3mDRbarPrimeToMSbarShift(pref*tc.getDRbarPrimeToMSbarShift(xtOrder,1,1)/v2);
   // this shift generates Xt^5*Log(mu) terms for the EFT expression
   ho.setDLambdaEFTDRbarPrimeToMSbarShift(pref*tc.getDRbarPrimeToMSbarShift(xtOrder,1,0)/v2);
   
   // set the uncertainty of delta_lambda due to missing Xt terms
   const int xt4Flag = xtOrder == 3 ? 1 : 0;
   ho.setDLambdaH3mXtUncertainty(pref*(xt4Flag*tc.getDRbarPrimeToMSbarXtTerms(tc.getLimit(), 4, 0)
      + tc.getDRbarPrimeToMSbarXtTerms(tc.getLimit(), 5, 0)
      + tc.getDRbarPrimeToMSbarXtTerms(tc.getLimit(), 6, 0))/v2);
   ho.setDLambdaEFTXtUncertainty(pref*(xt4Flag*tc.getDRbarPrimeToMSbarXtTerms(tc.getLimit(), 4, 1)
      + tc.getDRbarPrimeToMSbarXtTerms(tc.getLimit(), 5, 1)
      + tc.getDRbarPrimeToMSbarXtTerms(tc.getLimit(), 6, 0))/v2);
   
   // set uncertainty of delta_lambda
   ho.setDLambdaExpUncertainty(ho.getDMhExpUncertainty(3)/sqrt(v2));
   
   // calculate shifts needed to convert DR' to other renormalization schemes,
   // here one needs it with the minus sign to convert DR -> H3m
   ho.setDMhDRbarPrimeToH3mShift(-shiftH3mToDRbarPrime(ho));
   
   // set flags to omit all but O(at*as^n)
   mh2EFTCalculator.setCorrectionFlag(himalaya::EFTOrders::G12G22, 0);
   mh2EFTCalculator.setCorrectionFlag(himalaya::EFTOrders::G12YB2, 0);
   mh2EFTCalculator.setCorrectionFlag(himalaya::EFTOrders::G14, 0);
   mh2EFTCalculator.setCorrectionFlag(himalaya::EFTOrders::G24, 0);
   mh2EFTCalculator.setCorrectionFlag(himalaya::EFTOrders::G12YB2, 0);
   mh2EFTCalculator.setCorrectionFlag(himalaya::EFTOrders::G22YB2, 0);
   mh2EFTCalculator.setCorrectionFlag(himalaya::EFTOrders::YB4, 0);
   mh2EFTCalculator.setCorrectionFlag(himalaya::EFTOrders::G12YTAU2, 0);
   mh2EFTCalculator.setCorrectionFlag(himalaya::EFTOrders::G22YTAU2, 0);
   mh2EFTCalculator.setCorrectionFlag(himalaya::EFTOrders::YTAU4, 0);
   mh2EFTCalculator.setCorrectionFlag(himalaya::EFTOrders::G12YT2, 0);
   mh2EFTCalculator.setCorrectionFlag(himalaya::EFTOrders::G22YT2, 0);
   mh2EFTCalculator.setCorrectionFlag(himalaya::EFTOrders::G32YB4, 0);
   mh2EFTCalculator.setCorrectionFlag(himalaya::EFTOrders::YB6, 0);
   mh2EFTCalculator.setCorrectionFlag(himalaya::EFTOrders::YT6, 0);
   mh2EFTCalculator.setCorrectionFlag(himalaya::EFTOrders::YTAU2YB4, 0);
   mh2EFTCalculator.setCorrectionFlag(himalaya::EFTOrders::YTAU6, 0);
   mh2EFTCalculator.setCorrectionFlag(himalaya::EFTOrders::YT2YB4, 0);
   mh2EFTCalculator.setCorrectionFlag(himalaya::EFTOrders::YB2YT4, 0);
   // mh_eft^2
   const double mh2_eft = mh2EFTCalculator.getDeltaMh2EFT0Loop();
   // 1-Loop prefactor at
   const double pref_1L = 1./pow2(4*Pi) * pow2(p.Mt * gt);
   // 2-Loop prefactor at*as
   const double pref_2L = 1./pow4(4*Pi) * pow2(p.Mt * gt * p.g3);
   // fill in results of EFT calculation
   ho.setDMh2EFT(0, mh2_eft);
   ho.setDMh2EFT(1, mh2EFTCalculator.getDeltaMh2EFT1Loop(1, 1));
   ho.setDMh2EFT(2, mh2EFTCalculator.getDeltaMh2EFT2Loop(1, 1));
   ho.setDMh2EFT(3, mh2EFTCalculator.getDeltaMh2EFT3Loop(1, 1, 0) 
      + ho.getDLambdaEFT()*v2);
   ho.setDLambda(0, mh2_eft/v2);
   ho.setDLambda(1, pref_1L*(tc.getThresholdCorrection(
      ThresholdVariables::LAMBDA_AT,RenSchemes::DRBARPRIME, 1))/v2);
   ho.setDLambda(2, pref_2L*(tc.getThresholdCorrection(
      ThresholdVariables::LAMBDA_AT_AS, RenSchemes::DRBARPRIME, 1))/v2);
   ho.setDLambda(3, ho.getDLambdaEFT());
   ho.setDLambdaDRbarPrimeToMSbarShift(0, 0.);
   ho.setDLambdaDRbarPrimeToMSbarShift(1, 0.);
   ho.setDLambdaDRbarPrimeToMSbarShift(2, pref_2L*(-4*ho.getDLambda(1)
      *tc.getThresholdCorrection(ThresholdVariables::YT_AS,
				 RenSchemes::DRBARPRIME, 1))/v2);
   ho.setDLambdaDRbarPrimeToMSbarShift(3, ho.getDLambdaEFTDRbarPrimeToMSbarShift());
   //mh2EFTCalculator.getDeltaLambdaDegenerate(p.scale, p.MSt(1), p.Au(2,2) - p.mu*p.vd/p.vu,1);
   auto ho_mdr = ho;
   ho_mdr.setMDRFlag(1);
   // calculate the DR to MDR shift with the obtained hierarchy
   ho_mdr.setDMhDRbarPrimeToMDRbarPrimeShift(calcDRbarToMDRbarShift(ho_mdr, true, true));
   ho_mdr.setDMh(3, calculateHierarchy(ho_mdr, 0, 0, 1) + shiftH3mToDRbarPrime(ho_mdr));
   Eigen::Vector2d mdrMasses;
   mdrMasses(0) = ho_mdr.getMDRMasses()(0);
   mdrMasses(1) = ho_mdr.getMDRMasses()(1);
   ho.setMDRMasses(mdrMasses);
   ho.setDMhDRbarPrimeToMDRbarPrimeShift(ho_mdr.getDMhDRbarPrimeToMDRbarPrimeShift() 
      + ho_mdr.getDMh(3) - ho.getDMh(3));
   
   return ho;
}

/**
 * 	Compares deviation of all hierarchies with the exact two-loop result and returns the hierarchy which minimizes the error.
 * 	@param ho a HierarchyObject with constant isAlphab.
 * 	@return An integer which is identified with the suitable hierarchy.
 */
int himalaya::HierarchyCalculator::compareHierarchies(himalaya::HierarchyObject& ho){
   // set flags to truncate the expansion
   flagMap.at(ExpansionDepth::xx) = 0;
   flagMap.at(ExpansionDepth::xxMst) = 0;
   double error = -1.;
   int suitableHierarchy = -1;
   
   // sine of 2 times beta
   const double s2b = sin(2*atan(p.vu/p.vd));
   const double tbeta = p.vu/p.vd;
   
   // tree level Higgs mass matrix
   Eigen::Matrix2d treelvl;
   treelvl (0,0) = s2b/2.*(pow2(p.MZ) / tbeta + pow2(p.MA) * tbeta);
   treelvl (1,0) = s2b/2.*(-pow2(p.MZ) - pow2(p.MA));
   treelvl (0,1) = treelvl (1,0);
   treelvl (1,1) = s2b/2.*(pow2(p.MZ) * tbeta + pow2(p.MA) / tbeta);
   
   ho.setDMh(0, treelvl);

   // compare the exact higgs mass at 2-loop level with the expanded expressions to find a suitable hierarchy
   for (int hierarchy = Hierarchies::FIRST; hierarchy < Hierarchies::NUMBER_OF_HIERARCHIES; hierarchy++) {
      // first, check if the hierarchy is suitable to the mass spectrum
      ho.setSuitableHierarchy(hierarchy);
      
      if(isHierarchySuitable(ho)){
	 // calculate the exact 1-loop result (only alpha_t/b)
	 const Eigen::Matrix2d Mt41L = getMt41L(ho, ho.getMDRFlag(), 0);
	 
	 // call the routine of Pietro Slavich to get the alpha_s alpha_t/b corrections with the MDRbar masses
	 const Eigen::Matrix2d Mt42L = getMt42L(ho, ho.getMDRFlag(), 0);
	 
	 // Note: spurious poles are handled by the validate method
	 // of the Himalaya_Interface struct
	 
	 //calculate the exact Higgs mass at 2-loop (only up to alpha_s alpha_t/b)
	 const Eigen::EigenSolver<Eigen::Matrix2d> es2L (treelvl + Mt41L + Mt42L);
	 const double Mh2l = sortEigenvalues(es2L).at(0);

	 // calculate the expanded 2-loop expression with the specific hierarchy
	 const Eigen::EigenSolver<Eigen::Matrix2d> esExpanded (treelvl + Mt41L 
	    + calculateHierarchy(ho, 0, 1, 0));
	 
	 // calculate the higgs mass in the given mass hierarchy and compare the result to estimate the error
	 const double Mh2LExpanded = sortEigenvalues(esExpanded).at(0);

	 // estimate the error
	 const double twoLoopError = std::abs((Mh2l - Mh2LExpanded));

	 // estimate the uncertainty of the expansion at 2L
	 const double expUncertainty2L = getExpansionUncertainty(ho, treelvl 
	    + Mt41L, 0, 1, 0);

	 // estimate the uncertainty of the expansion at 3L
	 const double expUncertainty3L = getExpansionUncertainty(ho, treelvl
	    + Mt41L + Mt42L, 0, 0, 1); 

	 // add these errors to include the error of the expansion in the comparison
	 const double currError = sqrt(pow2(twoLoopError) 
	    + pow2(expUncertainty2L) + pow2(expUncertainty3L));

	 // if the error is negative, it is the first iteration and there is no hierarchy which fits better
	 if(error < 0){
	    error = currError;
	    suitableHierarchy = hierarchy;
	    ho.setAbsDiff2L(twoLoopError);
	    ho.setRelDiff2L(twoLoopError/Mh2l);
	    ho.setDMhExpUncertainty(2, expUncertainty2L);
	    ho.setDMhExpUncertainty(3, expUncertainty3L);
	 }
	 // compare the current error with the last error and choose the hierarchy which fits best (lowest error)
	 else if(currError < error){
	    error = currError;
	    suitableHierarchy = hierarchy;
	    ho.setAbsDiff2L(twoLoopError);
	    ho.setRelDiff2L(twoLoopError/Mh2l);
	    ho.setDMhExpUncertainty(2, expUncertainty2L);
	    ho.setDMhExpUncertainty(3, expUncertainty3L);
	 }
      }
   }
   ho.setSuitableHierarchy(suitableHierarchy);
   // reset the flags
   flagMap.at(ExpansionDepth::xx) = 1;
   flagMap.at(ExpansionDepth::xxMst) = 1;
   return suitableHierarchy;
}

/**
 * 	Checks if a hierarchy is suitable to the given mass spectrum.
 * 	@param ho a HierarchyObject with constant isAlphab and a hierarchy candidate.
 * 	@returns A bool if the hierarchy candidate is suitable to the given mass spectrum.
 */
bool himalaya::HierarchyCalculator::isHierarchySuitable(const himalaya::HierarchyObject& ho){
   double Mst1, Mst2;
   if(!ho.getIsAlphab()){
      Mst1 = p.MSt(0);
      Mst2 = p.MSt(1);
   }
   else{
      Mst1 = p.MSb(0);
      Mst2 = p.MSb(1);
   }
   
   // check if the squark mass is the heaviest mass
   const double delta = Msq*0.30;	// allow for an offset of 30%
   
   if(std::abs(Msq - Mst2) > delta) return false;
   if(std::abs(Msq - p.MG) > delta) return false;
   
   switch (ho.getSuitableHierarchy()){
      case Hierarchies::h3:
	 return Mgl > Mst2;
      case Hierarchies::h32q2g:
	 return (Mst2 >= Msq) && (Mst2 > Mgl);
      case Hierarchies::h3q22g:
	 return (Msq > Mst2) && (Mst2 > Mgl);
      case Hierarchies::h4:
	 return (Mst1 < Msq) && (Mst1 >= Mgl);
      case Hierarchies::h5:
	 return (Mst2 - Mst1 > 0.1*Mst1) && ((Mgl - Mst1) < std::abs(Mgl - Mst2)) && (Mst2 < Msq) && (Mst1 >= Mgl);
      case Hierarchies::h5g1:
	 return (Mst2 - Mst1 > 0.1*Mst1) && ((Mgl - Mst1) < std::abs(Mgl - Mst2)) && (Mst2 < Msq) && (Mgl > Mst1);
      case Hierarchies::h6:
	 return (Mst2 - Mst1 > 0.1*Mst1) && ((Mst2 - Mgl) < std::abs(Mgl - Mst1)) && (Mst2 < Msq) && (Mst2 >= Mgl);
      case Hierarchies::h6g2:
	 return (Mst2 - Mst1 > 0.1*Mst1) && ((Mst2 - Mgl) < std::abs(Mgl - Mst1)) && (Mst2 < Msq) && (Mgl > Mst2);
      case Hierarchies::h6b:
	 return (Mst2 - Mst1 > 0.1*Mst1) && ((Mst2 - Mgl) < std::abs(Mgl - Mst1)) && (Mst2 >= Msq) && (Mst2 >= Mgl);
      case Hierarchies::h6b2qg2:
	 return (Mst2 - Mst1 > 0.1*Mst1) && ((Mst2 - Mgl) < std::abs(Mgl - Mst1)) && (Mst2 >= Msq) && (Mgl > Mst2);
      case Hierarchies::h6bq22g:
	 return (Mst2 - Mst1 > 0.1*Mst1) && ((Mst2 - Mgl) < std::abs(Mgl - Mst1)) && (Msq > Mst2) && (Mst2 >= Mgl);
      case Hierarchies::h6bq2g2:
	 return (Mst2 - Mst1 > 0.1*Mst1) && ((Mst2 - Mgl) < std::abs(Mgl - Mst1)) && (Msq > Mst2) && (Mgl > Mst2);
      case Hierarchies::h9:
	 return (Mst2 >= Msq) && ((Mst2 - Mst1) < (Mst1 - Mgl));
      case Hierarchies::h9q2:
	 return (Msq > Mst2) && ((Mst1 - Mst1) < (Mst1 - Mgl));
   }
   return false;
}

// TODO(avoigt): if one is interested in the expansion at one- and two-loop choose a unified choice for the MDR scheme
/**
 * 	Calculates the hierarchy contributions for a specific hierarchy at a specific loop order.
 * 	@param ho a HierarchyObject with constant isAlphab.
 * 	@param oneLoopFlagIn an integer flag which is 0 or 1 in order to add or omit the expanded one-loop results to the returned value, respectivley.
 * 	@param twoLoopFlagIn an integer flag which is 0 or 1 in order to add or omit the expanded two-loop results to the returned value, respectivley.
 * 	@param threeLoopFlagIn an integer flag which is 0 or 1 in order to add or omit the expanded three-loop results to the returned value, respectivley.
 * 	@throws runtime_error Throws a runtime_error if the tree-level is requested in terms of hierarchies.
 * 	@return The loop corrected Higgs mass matrix which contains the expanded corrections at the given order.
 */
Eigen::Matrix2d himalaya::HierarchyCalculator::calculateHierarchy(himalaya::HierarchyObject& ho, const int oneLoopFlagIn,
								  const int twoLoopFlagIn, const int threeLoopFlagIn) {
   // get the hierarchy
   const int hierarchy = ho.getSuitableHierarchy();

   // the hierarchy files containing 1-, 2- and 3-loop terms (alpha_s^0 alpha_t/b, alpha_s alpha_t/b, alpha_s^2 alpha_t/b)
   double sigS1Full = 0., sigS2Full = 0., sigS12Full = 0.;

   // common variables
   double At, Mt, s2t, Mst1 = 0., Mst2 = 0.;
   if (!ho.getIsAlphab()) {
      At = p.Au(2,2);
      Mt = p.Mt;
      s2t = p.s2t;
   }
   else {
      At = p.Ad(2,2);
      Mt = p.Mb;
      s2t = p.s2b;
   }

   const double beta = atan(p.vu / p.vd);
   const double lmMt = log(pow2(p.scale / Mt));

   // this loop is needed to calculate the suitable mass shift order by order
   for(int currentLoopOrder = 1; currentLoopOrder <= 3; currentLoopOrder ++){
      bool runThisOrder;
      double curSig1 = 0., curSig2 = 0., curSig12 = 0.;
      int oneLoopFlag = 0, twoLoopFlag = 0, threeLoopFlag = 0;
      switch (currentLoopOrder){
	 case 1:
	    oneLoopFlag = 1;
	    runThisOrder = oneLoopFlag == oneLoopFlagIn;
	 break;
	 case 2:
	    twoLoopFlag = 1;
	    runThisOrder = twoLoopFlag == twoLoopFlagIn;
	 break;
	 case 3:
	    threeLoopFlag = 1;
	    runThisOrder = threeLoopFlag == threeLoopFlagIn;
	 break;
      }
      if(runThisOrder){
	 // set the Msx masses according to MDRFlag
	 if(oneLoopFlag == 1){
	    Mst1 = shiftMst1ToMDR(ho, 0, 0);
	    Mst2 = shiftMst2ToMDR(ho, 0, 0);
	 }
	 else if(twoLoopFlag == 1){
	    Mst1 = shiftMst1ToMDR(ho, ho.getMDRFlag(), 0);
	    Mst2 = shiftMst2ToMDR(ho, ho.getMDRFlag(), 0);
	 }
	 else if(threeLoopFlag == 1){
	    Mst1 = shiftMst1ToMDR(ho, ho.getMDRFlag(), ho.getMDRFlag());
	    Mst2 = shiftMst2ToMDR(ho, ho.getMDRFlag(), ho.getMDRFlag());
	 }
	 else{
	    throw std::runtime_error("There are no tree-level hierarchies included!");
	 }
	 // select the suitable hierarchy for the specific hierarchy and set variables
	 switch(getCorrectHierarchy(hierarchy)){
	    case Hierarchies::h3:{
	       const double Dmglst1 = Mgl - Mst1;
	       const double Dmsqst1 = pow2(Msq) - pow2(Mst1);
	       const double Dmst12 = pow2(Mst1) - pow2(Mst2);
	       const double lmMst1 = log(pow2(p.scale / Mst1));
	       switch(hierarchy){
		  case Hierarchies::h3:{
		     const H3 hierarchy3(flagMap, Al4p, beta,
			Dmglst1, Dmst12, Dmsqst1, lmMt, lmMst1,
			Mgl, Mt, Mst1, Mst2, Msq, p.mu,
			s2t, 
			ho.getMDRFlag(), oneLoopFlag, twoLoopFlag, threeLoopFlag);
		     curSig1 = hierarchy3.getS1();
		     curSig2 = hierarchy3.getS2();
		     curSig12 = hierarchy3.getS12();
		     if(oneLoopFlagIn == 0 && twoLoopFlagIn == 0 && threeLoopFlagIn == 1){
                        const double c = hierarchy3.calc_coef_at_as2_no_sm_logs_log0();
			ho.setDLambdaH3m(c
			   + lmMst1 * hierarchy3.calc_coef_at_as2_no_sm_logs_log1()
			   + pow2(lmMst1) * hierarchy3.calc_coef_at_as2_no_sm_logs_log2()
			   + pow3(lmMst1) * hierarchy3.calc_coef_at_as2_no_sm_logs_log3());
                        ho.setDLambdaNonLog(c);
		     }
		  }
		  break;
		  case Hierarchies::h32q2g:{
		     const H32q2g hierarchy32q2g(flagMap, Al4p, beta,
			Dmglst1, Dmst12, Dmsqst1, lmMt, lmMst1,
			Mt, Mst1, Mst2, p.mu,
			s2t,
			ho.getMDRFlag(), oneLoopFlag, twoLoopFlag, threeLoopFlag);
		     curSig1 = hierarchy32q2g.getS1();
		     curSig2 = hierarchy32q2g.getS2();
		     curSig12 = hierarchy32q2g.getS12();
		     if(oneLoopFlagIn == 0 && twoLoopFlagIn == 0 && threeLoopFlagIn == 1){
                        const double c = hierarchy32q2g.calc_coef_at_as2_no_sm_logs_log0();
			ho.setDLambdaH3m(c
			   + lmMst1 * hierarchy32q2g.calc_coef_at_as2_no_sm_logs_log1()
			   + pow2(lmMst1) * hierarchy32q2g.calc_coef_at_as2_no_sm_logs_log2()
			   + pow3(lmMst1) * hierarchy32q2g.calc_coef_at_as2_no_sm_logs_log3());
                        ho.setDLambdaNonLog(c);
		     }
		  }
		  break;
		  case Hierarchies::h3q22g:{
		     const H3q22g hierarchy3q22g(flagMap, Al4p, beta,
			Dmglst1, Dmst12, Dmsqst1, lmMt, lmMst1,
			Mt, Mst1, Mst2, Msq, p.mu,
			s2t,
			ho.getMDRFlag(), oneLoopFlag, twoLoopFlag, threeLoopFlag);
		     curSig1 = hierarchy3q22g.getS1();
		     curSig2 = hierarchy3q22g.getS2();
		     curSig12 = hierarchy3q22g.getS12();
		     if(oneLoopFlagIn == 0 && twoLoopFlagIn == 0 && threeLoopFlagIn == 1){
                        const double c = hierarchy3q22g.calc_coef_at_as2_no_sm_logs_log0();
			ho.setDLambdaH3m(c
			   + lmMst1 * hierarchy3q22g.calc_coef_at_as2_no_sm_logs_log1()
			   + pow2(lmMst1) * hierarchy3q22g.calc_coef_at_as2_no_sm_logs_log2()
			   + pow3(lmMst1) * hierarchy3q22g.calc_coef_at_as2_no_sm_logs_log3());
                        ho.setDLambdaNonLog(c);
		     }
		  }
		  break;
	       }
	    }
	    break;
	    case Hierarchies::h4:{
	       const double Msusy = (Mst1 + Mst2 + Mgl) / 3.;
	       const double lmMsusy = log(pow2(p.scale / Msusy));
	       const double lmMst1 = log(pow2(p.scale / Mst1));
	       const H4 hierarchy4(flagMap, Al4p, At, beta,
		  lmMt, lmMsq, lmMsusy, Mt, Msusy, Msq,
		  ho.getMDRFlag(), oneLoopFlag, twoLoopFlag, threeLoopFlag);
	       curSig1 = hierarchy4.getS1();
	       curSig2 = hierarchy4.getS2();
	       curSig12 = hierarchy4.getS12();
	       if(oneLoopFlagIn == 0 && twoLoopFlagIn == 0 && threeLoopFlagIn == 1){
                  const double c = hierarchy4.calc_coef_at_as2_no_sm_logs_log0();
		  ho.setDLambdaH3m(c
		     + lmMst1 * hierarchy4.calc_coef_at_as2_no_sm_logs_log1()
		     + pow2(lmMst1) * hierarchy4.calc_coef_at_as2_no_sm_logs_log2()
		     + pow3(lmMst1) * hierarchy4.calc_coef_at_as2_no_sm_logs_log3());
                  ho.setDLambdaNonLog(c);
	       }
	    }
	    break;
	    case Hierarchies::h5:{
	       const double Dmglst1 = Mgl - Mst1;
	       const double lmMst1 = log(pow2(p.scale / Mst1));
	       const double lmMst2 = log(pow2(p.scale / Mst2));
	       switch(hierarchy){
		  case Hierarchies::h5:{
		     const H5 hierarchy5(flagMap, Al4p, beta, Dmglst1,
			lmMt, lmMst1, lmMst2, lmMsq, Mt, Mst1,
			Mst2, Msq, p.mu,
			s2t,
			ho.getMDRFlag(), oneLoopFlag, twoLoopFlag, threeLoopFlag);
		     curSig1 = hierarchy5.getS1();
		     curSig2 = hierarchy5.getS2();
		     curSig12 = hierarchy5.getS12();
		     if(oneLoopFlagIn == 0 && twoLoopFlagIn == 0 && threeLoopFlagIn == 1){
                        const double c = hierarchy5.calc_coef_at_as2_no_sm_logs_log0();
			ho.setDLambdaH3m(c
			   + lmMst1 * hierarchy5.calc_coef_at_as2_no_sm_logs_log1()
			   + pow2(lmMst1) * hierarchy5.calc_coef_at_as2_no_sm_logs_log2()
			   + pow3(lmMst1) * hierarchy5.calc_coef_at_as2_no_sm_logs_log3());
                        ho.setDLambdaNonLog(c);
		     }
		  }
		  break;
		  case Hierarchies::h5g1:{
		     const H5g1 hierarchy5g1(flagMap, Al4p, beta, Dmglst1,
			lmMt, lmMst1, lmMst2, lmMsq, Mgl, Mt, Mst1,
			Mst2, Msq, p.mu,
			s2t,
			ho.getMDRFlag(), oneLoopFlag, twoLoopFlag, threeLoopFlag);
		     curSig1 = hierarchy5g1.getS1();
		     curSig2 = hierarchy5g1.getS2();
		     curSig12 = hierarchy5g1.getS12();
		     if(oneLoopFlagIn == 0 && twoLoopFlagIn == 0 && threeLoopFlagIn == 1){
                        const double c = hierarchy5g1.calc_coef_at_as2_no_sm_logs_log0();
			ho.setDLambdaH3m(c
			   + lmMst1 * hierarchy5g1.calc_coef_at_as2_no_sm_logs_log1()
			   + pow2(lmMst1) * hierarchy5g1.calc_coef_at_as2_no_sm_logs_log2()
			   + pow3(lmMst1) * hierarchy5g1.calc_coef_at_as2_no_sm_logs_log3());
                        ho.setDLambdaNonLog(c);
		     }
		  }
		  break;
	       }
	    }
	    break;
	    case Hierarchies::h6:{
	       const double Dmglst2 = Mgl - Mst2;
	       const double lmMst1 = log(pow2(p.scale / Mst1));
	       const double lmMst2 = log(pow2(p.scale / Mst2));
	       switch(hierarchy){
		  case Hierarchies::h6:{
		     const H6 hierarchy6(flagMap, Al4p, beta, Dmglst2,
			lmMt, lmMst1, lmMst2, lmMsq,
			Mt, Mst1, Mst2, Msq, p.mu,
			s2t,
			ho.getMDRFlag(), oneLoopFlag, twoLoopFlag, threeLoopFlag);
		     curSig1 = hierarchy6.getS1();
		     curSig2 = hierarchy6.getS2();
		     curSig12 = hierarchy6.getS12();
		     if(oneLoopFlagIn == 0 && twoLoopFlagIn == 0 && threeLoopFlagIn == 1){
                        const double c = hierarchy6.calc_coef_at_as2_no_sm_logs_log0();
			ho.setDLambdaH3m(c
			   + lmMst1 * hierarchy6.calc_coef_at_as2_no_sm_logs_log1()
			   + pow2(lmMst1) * hierarchy6.calc_coef_at_as2_no_sm_logs_log2()
			   + pow3(lmMst1) * hierarchy6.calc_coef_at_as2_no_sm_logs_log3());
                        ho.setDLambdaNonLog(c);
		     };
		  }
		  break;
		  case Hierarchies::h6g2:{
		     const H6g2 hierarchy6g2(flagMap, Al4p, beta, Dmglst2,
			lmMt, lmMst1, lmMst2, lmMsq,
			Mgl, Mt, Mst1, Mst2, Msq, p.mu,
			s2t,
			ho.getMDRFlag(), oneLoopFlag, twoLoopFlag, threeLoopFlag);
		     curSig1 = hierarchy6g2.getS1();
		     curSig2 = hierarchy6g2.getS2();
		     curSig12 = hierarchy6g2.getS12();
		     if(oneLoopFlagIn == 0 && twoLoopFlagIn == 0 && threeLoopFlagIn == 1){
                        const double c = hierarchy6g2.calc_coef_at_as2_no_sm_logs_log0();
			ho.setDLambdaH3m(c
			   + lmMst1 * hierarchy6g2.calc_coef_at_as2_no_sm_logs_log1()
			   + pow2(lmMst1) * hierarchy6g2.calc_coef_at_as2_no_sm_logs_log2()
			   + pow3(lmMst1) * hierarchy6g2.calc_coef_at_as2_no_sm_logs_log3());
                        ho.setDLambdaNonLog(c);
		     }
		  }
		  break;
	       }
	    }
	    break;
	    case Hierarchies::h6b:{
	       const double Dmglst2 = Mgl - Mst2;
	       const double Dmsqst2 = Msq - Mst2;
	       const double lmMst1 = log(pow2(p.scale / Mst1));
	       const double lmMst2 = log(pow2(p.scale / Mst2));
	       switch(hierarchy){
		  case Hierarchies::h6b:{
		     const H6b hierarchy6b(flagMap, Al4p, beta, Dmglst2,
			Dmsqst2, lmMt, lmMst1, lmMst2,
			Mt, Mst1, Mst2, p.mu,
			s2t,
			ho.getMDRFlag(), oneLoopFlag, twoLoopFlag, threeLoopFlag);
		     curSig1 = hierarchy6b.getS1();
		     curSig2 = hierarchy6b.getS2();
		     curSig12 = hierarchy6b.getS12();
		     if(oneLoopFlagIn == 0 && twoLoopFlagIn == 0 && threeLoopFlagIn == 1){
                        const double c = hierarchy6b.calc_coef_at_as2_no_sm_logs_log0();
			ho.setDLambdaH3m(c
			   + lmMst1 * hierarchy6b.calc_coef_at_as2_no_sm_logs_log1()
			   + pow2(lmMst1) * hierarchy6b.calc_coef_at_as2_no_sm_logs_log2()
			   + pow3(lmMst1) * hierarchy6b.calc_coef_at_as2_no_sm_logs_log3());
                        ho.setDLambdaNonLog(c);
		     }
		  }
		  break;
		  case Hierarchies::h6b2qg2:{
		     const H6b2qg2 hierarchy6b2qg2(flagMap, Al4p, beta, Dmglst2,
			Dmsqst2, lmMt, lmMst1, lmMst2,
			Mgl, Mt, Mst1, Mst2, p.mu,
			s2t,
			ho.getMDRFlag(), oneLoopFlag, twoLoopFlag, threeLoopFlag);
		     curSig1 = hierarchy6b2qg2.getS1();
		     curSig2 = hierarchy6b2qg2.getS2();
		     curSig12 = hierarchy6b2qg2.getS12();
		     if(oneLoopFlagIn == 0 && twoLoopFlagIn == 0 && threeLoopFlagIn == 1){
                        const double c = hierarchy6b2qg2.calc_coef_at_as2_no_sm_logs_log0();
			ho.setDLambdaH3m(c
			   + lmMst1 * hierarchy6b2qg2.calc_coef_at_as2_no_sm_logs_log1()
			   + pow2(lmMst1) * hierarchy6b2qg2.calc_coef_at_as2_no_sm_logs_log2()
			   + pow3(lmMst1) * hierarchy6b2qg2.calc_coef_at_as2_no_sm_logs_log3());
                        ho.setDLambdaNonLog(c);
		     }
		  }
		  break;
		  case Hierarchies::h6bq22g:{
		     const H6bq22g hierarchy6bq22g(flagMap, Al4p, beta, Dmglst2,
			Dmsqst2, lmMt, lmMst1, lmMst2,
			Mt, Mst1, Mst2, Msq, p.mu,
			s2t,
			ho.getMDRFlag(), oneLoopFlag, twoLoopFlag, threeLoopFlag);
		     curSig1 = hierarchy6bq22g.getS1();
		     curSig2 = hierarchy6bq22g.getS2();
		     curSig12 = hierarchy6bq22g.getS12();
		     if(oneLoopFlagIn == 0 && twoLoopFlagIn == 0 && threeLoopFlagIn == 1){
                        const double c = hierarchy6bq22g.calc_coef_at_as2_no_sm_logs_log0();
			ho.setDLambdaH3m(c
			   + lmMst1 * hierarchy6bq22g.calc_coef_at_as2_no_sm_logs_log1()
			   + pow2(lmMst1) * hierarchy6bq22g.calc_coef_at_as2_no_sm_logs_log2()
			   + pow3(lmMst1) * hierarchy6bq22g.calc_coef_at_as2_no_sm_logs_log3());
                        ho.setDLambdaNonLog(c);
		     }
		  }
		  break;
		  case Hierarchies::h6bq2g2:{
		     const H6bq2g2 hierarchy6bq2g2(flagMap, Al4p, beta, Dmglst2,
			Dmsqst2, lmMt, lmMst1, lmMst2,
			Mgl, Mt, Mst1,Mst2, Msq, p.mu,
			s2t,
			ho.getMDRFlag(), oneLoopFlag, twoLoopFlag, threeLoopFlag);
		     curSig1 = hierarchy6bq2g2.getS1();
		     curSig2 = hierarchy6bq2g2.getS2();
		     curSig12 = hierarchy6bq2g2.getS12();
		     if(oneLoopFlagIn == 0 && twoLoopFlagIn == 0 && threeLoopFlagIn == 1){
                        const double c = hierarchy6bq2g2.calc_coef_at_as2_no_sm_logs_log0();
			ho.setDLambdaH3m(c
			   + lmMst1 * hierarchy6bq2g2.calc_coef_at_as2_no_sm_logs_log1()
			   + pow2(lmMst1) * hierarchy6bq2g2.calc_coef_at_as2_no_sm_logs_log2()
			   + pow3(lmMst1) * hierarchy6bq2g2.calc_coef_at_as2_no_sm_logs_log3());
                        ho.setDLambdaNonLog(c);
		     }
		  }
		  break;
	       }
	    }
	    break;
	    case Hierarchies::h9:{
	       const double lmMst1 = log(pow2(p.scale / Mst1));
	       const double Dmst12 = pow2(Mst1) - pow2(Mst2);
	       const double Dmsqst1 = pow2(Msq) - pow2(Mst1);
	       switch(hierarchy){
		  case Hierarchies::h9:{
		     const H9 hierarchy9(flagMap, Al4p, beta, Dmst12, Dmsqst1,
			lmMt, lmMgl, lmMst1,
			Mgl, Mt, Mst1, Mst2, p.mu,
			s2t,
			ho.getMDRFlag(), oneLoopFlag, twoLoopFlag, threeLoopFlag);
		     curSig1 = hierarchy9.getS1();
		     curSig2 = hierarchy9.getS2();
		     curSig12 = hierarchy9.getS12();
		     if(oneLoopFlagIn == 0 && twoLoopFlagIn == 0 && threeLoopFlagIn == 1){
                        const double c = hierarchy9.calc_coef_at_as2_no_sm_logs_log0();
			ho.setDLambdaH3m(c
			   + lmMst1 * hierarchy9.calc_coef_at_as2_no_sm_logs_log1()
			   + pow2(lmMst1) * hierarchy9.calc_coef_at_as2_no_sm_logs_log2()
			   + pow3(lmMst1) * hierarchy9.calc_coef_at_as2_no_sm_logs_log3());
                        ho.setDLambdaNonLog(c);
		     }
		  }
		  break;
		  case Hierarchies::h9q2:{
		     const H9q2 hierarchy9q2(flagMap, Al4p, beta, Dmst12, Dmsqst1,
			lmMt, lmMgl, lmMst1,
			Mgl, Mt, Mst1, Mst2, Msq, p.mu,
			s2t,
			ho.getMDRFlag(), oneLoopFlag, twoLoopFlag, threeLoopFlag);
		     curSig1 = hierarchy9q2.getS1();
		     curSig2 = hierarchy9q2.getS2();
		     curSig12 = hierarchy9q2.getS12();
		     if(oneLoopFlagIn == 0 && twoLoopFlagIn == 0 && threeLoopFlagIn == 1){
                        const double c = hierarchy9q2.calc_coef_at_as2_no_sm_logs_log0();
			ho.setDLambdaH3m(c
			   + lmMst1 * hierarchy9q2.calc_coef_at_as2_no_sm_logs_log1()
			   + pow2(lmMst1) * hierarchy9q2.calc_coef_at_as2_no_sm_logs_log2()
			   + pow3(lmMst1) * hierarchy9q2.calc_coef_at_as2_no_sm_logs_log3());
                        ho.setDLambdaNonLog(c);
		     }
		  }
		  break;
	       }
	    }
	    break;
	 }
      }
      sigS1Full += curSig1;
      sigS2Full += curSig2;
      sigS12Full += curSig12;
   }
   
   // add the MDR masses to the hierarchy object only if a 3-loop calculation has to be done, otherwise let the user decide
   if(oneLoopFlagIn == 0 && twoLoopFlagIn == 0 && threeLoopFlagIn == 1){
      Eigen::Vector2d mdrMasses;
      mdrMasses(0) = Mst1;
      mdrMasses(1) = Mst2;
      ho.setMDRMasses(mdrMasses);
   }

   Eigen::Matrix2d higgsMassMatrix;
   higgsMassMatrix(0, 0) = prefac * sigS1Full;
   higgsMassMatrix(0, 1) = prefac * sigS12Full;
   higgsMassMatrix(1, 0) = higgsMassMatrix(0, 1);
   higgsMassMatrix(1, 1) = prefac * sigS2Full;
   return higgsMassMatrix;
}

/**
 * 	Shifts Msx1 according to the hierarchy to the MDR scheme.
 * 	@param ho a HierarchyObject with constant isAlphab.
 * 	@param oneLoopFlag an integer flag which is 0 or 1 in order to shift the order O(alpha_s).
 * 	@param twoLoopFlag an integer flag which is 0 or 1 in order to shift the order O(alpha_s^2).
 * 	@return A double which is the MDR sx_1 mass.
 */
double himalaya::HierarchyCalculator::shiftMst1ToMDR(const himalaya::HierarchyObject& ho,
						     const unsigned int oneLoopFlag,
						     const unsigned int twoLoopFlag) {
   double Mst1mod = 0., Mst1, Mst2;
   if(!ho.getIsAlphab()){
      Mst1 = p.MSt(0);
      Mst2 = p.MSt(1);
   }
   else{
      Mst1 = p.MSb(0);
      Mst2 = p.MSb(1);
   }
   const double lmMst2 = log(pow2(p.scale) / pow2(Mst2));
   const double Dmglst2 = Mgl - Mst2;
   const double mdr2mst1ka = (-8. * twoLoopFlag * pow2(Al4p)
      * (10 * pow2(Msq) * (-1 + 2 * lmMsq + 2 * z2) + pow2(Mst2)
      * (-1 + 2 * lmMst2 + 2 * z2))) / (3. * pow2(Mst1));
   switch (getCorrectHierarchy(ho.getSuitableHierarchy())) {
   case Hierarchies::h3:
      Mst1mod = (1 + mdr2mst1ka);
      break;
   case Hierarchies::h4:
      Mst1mod = (1 + mdr2mst1ka);
      break;
   case Hierarchies::h5:
      Mst1mod = (1 + mdr2mst1ka);
      break;
   case Hierarchies::h6:
      Mst1mod = (144 * oneLoopFlag * Al4p * (1 + lmMgl) * pow2(Mgl) * pow4(Msq) + 27 * (1 + mdr2mst1ka) * pow4(Msq) * pow2(Mst1) +
         twoLoopFlag * pow2(Al4p) * Mgl * (-5 * (67 + 84 * lmMgl - 84 * lmMsq) * pow5(Mgl) - 40 * (43 + 30 * lmMgl - 30 * lmMsq) * pow3(Mgl) * pow2(Msq) +
            288 * Dmglst2 * pow4(Msq) * (1 - 2 * z2) + 12 * Mgl * pow4(Msq) * (79 + 144 * pow2(lmMgl) - 150 * lmMsq +
               90 * pow2(lmMsq) - 90 * lmMgl * (-3 + 2 * lmMsq) + 208 * z2))) / (27. * pow4(Msq) * pow2(Mst1));
      break;
   case Hierarchies::h6b:
      Mst1mod = (48 * oneLoopFlag * Al4p * (1 + lmMgl) * pow2(Mgl) + 9 * (1 + mdr2mst1ka) * pow2(Mst1) +
         8 * twoLoopFlag * pow2(Al4p) * (-135 * pow2(Msq) + 12 * Dmglst2 * Mgl * (1 - 22 * z2) +
            pow2(Mgl) * (77 + 135 * lmMgl + 72 * pow2(lmMgl) - 75 * lmMsq -
               90 * lmMgl * lmMsq + 45 * pow2(lmMsq) + 104 * z2))) / (9. * pow2(Mst1));
      break;
   case Hierarchies::h9:
      Mst1mod = (1 + mdr2mst1ka);
      break;
   }
   return Mst1 * sqrt(Mst1mod);
}

/**
 * 	Shifts Msx2 according to the hierarchy to the MDR scheme.
 * 	@param ho a HierarchyObject with constant isAlphab.
 * 	@param oneLoopFlag an integer flag which is 0 or 1 in order to shift the order O(alpha_s).
 * 	@param twoLoopFlag an integer flag which is 0 or 1 in order to shift the order O(alpha_s^2).
 * 	@return A double which is the MDR sx_2 mass.
 */
double himalaya::HierarchyCalculator::shiftMst2ToMDR(const himalaya::HierarchyObject& ho,
						     const unsigned int oneLoopFlag,
						     const unsigned int twoLoopFlag) {
   double Mst2mod = 0., Mst2;
   if(!ho.getIsAlphab()){
      Mst2 = p.MSt(1);
   }
   else{
      Mst2 = p.MSb(1);
   }
   const double Dmglst2 = Mgl - Mst2;
   const double mdr2mst2ka = (-80. * twoLoopFlag * pow2(Al4p) 
      * pow2(Msq) * (-1 + 2 * lmMsq + 2 * z2)) / (3. * pow2(Mst2));
   switch (getCorrectHierarchy(ho.getSuitableHierarchy())) {
   case Hierarchies::h3:
      Mst2mod = (1 + mdr2mst2ka);
      break;
   case Hierarchies::h4:
      Mst2mod = (1 + mdr2mst2ka);
      break;
   case Hierarchies::h5:
      Mst2mod = (1 + mdr2mst2ka);
      break;
   case Hierarchies::h6:
      Mst2mod = (144 * oneLoopFlag * Al4p * (1 + lmMgl) * pow2(Mgl) * pow4(Msq) + 27 * (1 + mdr2mst2ka) * pow4(Msq) * pow2(Mst2) +
         twoLoopFlag * pow2(Al4p) * Mgl * (-5 * (67 + 84 * lmMgl - 84 * lmMsq) * pow5(Mgl) - 40 * (43 + 30 * lmMgl - 30 * lmMsq) * pow3(Mgl) * pow2(Msq) +
            288 * Dmglst2 * pow4(Msq) * (1 - 2 * z2) + 12 * Mgl * pow4(Msq) * (79 + 144 * pow2(lmMgl) - 150 * lmMsq +
               90 * pow2(lmMsq) - 90 * lmMgl * (-3 + 2 * lmMsq) + 208 * z2))) / (27. * pow4(Msq) * pow2(Mst2));
      break;
   case Hierarchies::h6b:
      Mst2mod = (48 * oneLoopFlag * Al4p * (1 + lmMgl) * pow2(Mgl) + 9 * (1 + mdr2mst2ka) * pow2(Mst2) +
         8 * twoLoopFlag * pow2(Al4p) * (-135 * pow2(Msq) + 12 * Dmglst2 * Mgl * (1 - 22 * z2) +
            pow2(Mgl) * (77 + 135 * lmMgl + 72 * pow2(lmMgl) - 75 * lmMsq -
               90 * lmMgl * lmMsq + 45 * pow2(lmMsq) + 104 * z2))) / (9. * pow2(Mst2));
      break;
   case Hierarchies::h9:
      Mst2mod = (1 + mdr2mst2ka);
      break;
   }
   return Mst2 * sqrt(Mst2mod);
}

/**
 * 	Shifts the H3m renormalization scheme to DR' scheme. This shift has to be added to the H3m result!
 * 	@param ho a HierarchyObject with constant isAlphab.
 * 	@return A matrix which shifts the H3m scheme to the DR' scheme at three-loop level
 *
 */
Eigen::Matrix2d himalaya::HierarchyCalculator::shiftH3mToDRbarPrime(const himalaya::HierarchyObject& ho){
   Eigen::Matrix2d shift;
   
   // truncate shift at O(Xt^2) to be consistent with H3m result
   int truncateXt = 1;
   const int suitableHierarchy = ho.getSuitableHierarchy();
   if(suitableHierarchy == himalaya::Hierarchies::h3
      || suitableHierarchy == himalaya::Hierarchies::h32q2g 
      || suitableHierarchy == himalaya::Hierarchies::h3q22g
      || suitableHierarchy == himalaya::Hierarchies::h9
      || suitableHierarchy == himalaya::Hierarchies::h9q2) truncateXt = 0;
   
   // pre-factor of shift -> checked normalization against H3m normalization and they coincide
   const double k = 1 / (16 * pow2(Pi));
   const double yt = sqrt(2) * p.Mt / p.vu;
   const double prefac = pow4(p.g3) * pow3(k) * pow2(p.Mt * yt);
   
   // tanbeta
   const double Tbeta = p.vu / p.vd;
   
   // stop masses
   const double Mst1 = shiftMst1ToMDR(ho, ho.getMDRFlag(), ho.getMDRFlag());
   const double Mst2 = shiftMst2ToMDR(ho, ho.getMDRFlag(), ho.getMDRFlag());
   double Xt = p.Au(2,2) - p.mu / Tbeta;
   // Hierarchy h4 only covers O(Xt^0)
   if(suitableHierarchy == himalaya::Hierarchies::h4) Xt = 0;
   
   // threshold for degenerate squark mass case is 1% of the stop mass
   const double eps = Mst1 * 0.01;
   
   // squared masses
   const double Mst12 = pow2(Mst1);
   const double Mgl2 = pow2(p.MG);
   const double Msq2 = pow2(Msq);
   const double scale2 = pow2(p.scale);
   const double Xt2 = pow2(Xt);
   const double Mst22 = pow2(Mst2);
   const double Dmst12 = Mst12 - Mst22;
   
   // logarithms
   const double lmMst1 = log(scale2 / Mst12);
   const double lmMst2 = log(scale2 / Mst22);
   const double lmMsq = log(scale2 / Msq2);
   const double lmMgl = log(scale2 / Mgl2);
   
   // degenerate mass case flag
   bool isDegen = false;
   
   // check for degenerate squark masses
   if(std::abs(Mst1 - Mst2) < eps){
      const double Mst2shift = Mst1 + sqrt(std::abs(Dmst12))/2.;
      const double lmMst2shift = log(scale2 / pow2(Mst2shift));
      // limit
      const double lim = (32*Xt2*(-3*(1 + lmMgl)*pow2(Mgl) + 5*(1 + lmMsq)*pow2(Msq) + (1 +
        lmMst1)*pow2(Mst1))*pow2(p.mu))/(3.*pow6(Mst1));
      // exact result
      const double exact = (16*Xt2*(-6*(1 + lmMgl)*pow2(Mgl) + 10*(1 + lmMsq)*pow2(Msq) + (1 +
        lmMst1)*pow2(Mst1) + (1 + lmMst2)*pow2(Mst2))*pow2(p.mu)*(-4*log(Mst1/
        Mst2)*pow2(Mst1)*pow2(Mst2) + pow4(Mst1) - pow4(Mst2)))/(pow2(Mst1)*
        pow2(Mst2)*pow3(pow2(Mst1) - pow2(Mst2)));
      // exact result with shifted stop_2 mass
      const double exactShifted = (16*Xt2*(-6*(1 + lmMgl)*pow2(Mgl) + 10*(1 + lmMsq)*pow2(Msq) + (1 +
        lmMst1)*pow2(Mst1) + (1 + lmMst2shift)*pow2(Mst2shift))*pow2(p.mu)*(-4*log(Mst1/
        Mst2shift)*pow2(Mst1)*pow2(Mst2shift) + pow4(Mst1) - pow4(Mst2shift)))/(pow2(Mst1)*
        pow2(Mst2shift)*pow3(pow2(Mst1) - pow2(Mst2shift)));
      
      isDegen = std::abs(exactShifted - lim) >= std::abs(exact - lim)
	 || !std::isfinite(exact) || !std::isfinite(exactShifted);
   }

   if(isDegen){
      // matrix elements
      shift(0, 0) = (32*Xt2*(-3*(1 + lmMgl)*pow2(Mgl) + 5*(1 + lmMsq)*pow2(Msq) + (1 +
        lmMst1)*pow2(Mst1))*pow2(p.mu))/(3.*pow6(Mst1));
      shift(1, 0) = (-32*p.mu*(-3*(1 + lmMgl)*pow2(Mgl) + 5*(1 + lmMsq)*pow2(Msq) + (1 +
        lmMst1)*pow2(Mst1))*(p.mu*Xt2 - 3*Tbeta*Xt*pow2(Mst1) + Tbeta*pow3(Xt)))/
        (3.*Tbeta*pow6(Mst1));
      shift(0, 1) = shift(1,0);
      shift(1, 1) = (32*(-3*(1 + lmMgl)*pow2(Mgl) + 5*(1 + lmMsq)*pow2(Msq) + (1 + lmMst1)*
        pow2(Mst1))*(-6*Tbeta*(p.mu*Xt + Tbeta*Xt2)*pow2(Mst1) + Xt2*pow2(p.mu) +
        truncateXt*pow2(Tbeta)*pow2(Xt2) + 2*p.mu*Tbeta*pow3(Xt) + 6*pow2(Tbeta)*
        pow4(Mst1)))/(3.*pow2(Tbeta)*pow6(Mst1));
   }
   else{
      // matrix elements
      shift(0, 0) = (16*Xt2*(-6*(1 + lmMgl)*pow2(Mgl) + 10*(1 + lmMsq)*pow2(Msq) + (1 +
        lmMst1)*pow2(Mst1) + (1 + lmMst2)*pow2(Mst2))*pow2(p.mu)*(-4*log(Mst1/
        Mst2)*pow2(Mst1)*pow2(Mst2) + pow4(Mst1) - pow4(Mst2)))/(pow2(Mst1)*
        pow2(Mst2)*pow3(pow2(Mst1) - pow2(Mst2)));
      shift(1, 0) = (16*p.mu*(-6*(1 + lmMgl)*pow2(Mgl) + 10*(1 + lmMsq)*pow2(Msq) + (1 +
        lmMst1)*pow2(Mst1) + (1 + lmMst2)*pow2(Mst2))*(4*log(Mst1/Mst2)*pow2(
        Mst1)*pow2(Mst2)*(p.mu*Xt2 + Tbeta*pow3(Xt)) + p.mu*Xt2*(-pow4(Mst1) +
        pow4(Mst2)) + Tbeta*Xt*(pow3(pow2(Mst1) - pow2(Mst2)) + pow2(Xt)*(-
        pow4(Mst1) + pow4(Mst2)))))/(Tbeta*pow2(Mst1)*pow2(Mst2)*pow3(pow2(
        Mst1) - pow2(Mst2)));
      shift(0, 1) = shift(1,0);
      shift(1, 1) = (16*(-6*(1 + lmMgl)*pow2(Mgl) + 10*(1 + lmMsq)*pow2(Msq) + (1 + lmMst1)*
        pow2(Mst1) + (1 + lmMst2)*pow2(Mst2))*(-4*log(Mst1/Mst2)*pow2(Mst1)*
        pow2(Mst2)*(Xt2*pow2(p.mu) + truncateXt*pow2(Tbeta)*pow2(Xt2) + 2*p.mu*
        Tbeta*pow3(Xt)) + Xt2*(-2*pow2(Tbeta)*pow3(pow2(Mst1) - pow2(Mst2)) +
        pow2(p.mu)*(pow4(Mst1) - pow4(Mst2))) + Tbeta*(-2*p.mu*Xt*pow3(pow2(Mst1) -
        pow2(Mst2)) + Tbeta*(pow2(Mst1) + pow2(Mst2))*pow3(pow2(Mst1) - pow2(
        Mst2)) + 2*p.mu*pow3(Xt)*(pow4(Mst1) - pow4(Mst2))) + truncateXt*pow2(
        Tbeta)*pow2(Xt2)*(pow4(Mst1) - pow4(Mst2))))/(pow2(Mst1)*pow2(Mst2)*
        pow2(Tbeta)*pow3(pow2(Mst1) - pow2(Mst2)));
   }

   return prefac * shift;
}

/**
 * 	Shifts the H3m renormalization scheme to DR' scheme. This shift has to be added to the H3m result!
 * 	Note: This shift is WITHOUT the three-loop pre-factor g3^4*k^3*Mt^2*yt^2*Sin[beta]^2 with k = 1 / (16 Pi^2)
 * 	@param ho a HierarchyObject with constant isAlphab.
 * 	@return A double which shifts the H3m scheme to the DR' scheme at three-loop level
 *
 */
double himalaya::HierarchyCalculator::shiftH3mToDRbarPrimeMh2(const himalaya::HierarchyObject& ho, int omitLogs){
   double shift;
   
   // truncate shift at O(Xt^2) to be consistent with H3m result
   int truncateXt = 1;
   const int suitableHierarchy = ho.getSuitableHierarchy();
   if(suitableHierarchy == himalaya::Hierarchies::h3
      || suitableHierarchy == himalaya::Hierarchies::h32q2g 
      || suitableHierarchy == himalaya::Hierarchies::h3q22g
      || suitableHierarchy == himalaya::Hierarchies::h9
      || suitableHierarchy == himalaya::Hierarchies::h9q2) truncateXt = 0;
   
   // tanbeta
   const double Tbeta = p.vu / p.vd;
   
   // stop masses
   const double Mst1 = p.MSt(0);
   const double Mst2 = p.MSt(1);
   
   double Xt = p.Au(2,2) - p.mu / Tbeta;
   // Hierarchy h4 only covers O(Xt^0)
   if(suitableHierarchy == himalaya::Hierarchies::h4) Xt = 0;

   // threshold for degenerate squark mass case is 1% of the stop mass
   const double eps = Mst1 * 0.01;
   
   // squared masses
   const double Mst12 = pow2(Mst1);
   const double Mgl2 = pow2(p.MG);
   const double Msq2 = pow2(Msq);
   const double scale2 = pow2(p.scale);
   const double Xt2 = pow2(Xt);
   const double Mst22 = pow2(Mst2);
   const double Dmst12 = Mst12 - Mst22;

   // logarithms
   const double lmMst1 = omitLogs * log(scale2 / Mst12);
   const double lmMst2 = omitLogs * log(scale2 / Mst12) - log(Mst22 / Mst12);
   const double lmMsq = omitLogs * log(scale2 / Mst12) - log(Msq2 / Mst12);
   const double lmMgl = omitLogs * log(scale2 / Mst12) - log(Mgl2 / Mst12);
   
   // degenerate mass case flag
   bool isDegen = false;
   
   // check for degenerate squark masses
   if(std::abs(Mst1 - Mst2) < eps){
      const double Mst2shift = Mst1 + sqrt(std::abs(Dmst12))/2.;
      const double lmMst2shift = log(scale2 / pow2(Mst2shift));
      // limit
      const double lim = (32*(-3*(1 + lmMgl)*pow2(Mgl) + 5*(1 + lmMsq)*pow2(Msq) + (1 + lmMst1)*
        pow2(Mst1))*(-6*Xt2*pow2(Mst1) + truncateXt*pow2(Xt2) + 6*pow4(Mst1)))/
        (3.*pow6(Mst1));
      // exact result
      const double exact =  (16*(-6*(1 + lmMgl)*pow2(Mgl) + 10*(1 + lmMsq)*pow2(Msq) + (1 + lmMst1)*
        pow2(Mst1) + (1 + lmMst2)*pow2(Mst2))*(4*truncateXt*log(Mst2/Mst1)*
        pow2(Mst1)*pow2(Mst2)*pow2(Xt2) - 2*Xt2*pow3(pow2(Mst1) - pow2(Mst2)) +
        (pow2(Mst1) + pow2(Mst2))*pow3(pow2(Mst1) - pow2(Mst2)) + truncateXt*
        pow2(Xt2)*(pow4(Mst1) - pow4(Mst2))))/(pow2(Mst1)*pow2(Mst2)*pow3(pow2(
        Mst1) - pow2(Mst2)));
      // exact result with shifted stop_2 mass
      const double exactShifted =  (16*(-6*(1 + lmMgl)*pow2(Mgl) + 10*(1 + lmMsq)*pow2(Msq) + (1 + lmMst1)*
        pow2(Mst1) + (1 + lmMst2shift)*pow2(Mst2shift))*(4*truncateXt*log(Mst2shift/Mst1)*
        pow2(Mst1)*pow2(Mst2shift)*pow2(Xt2) - 2*Xt2*pow3(pow2(Mst1) - pow2(Mst2shift)) +
        (pow2(Mst1) + pow2(Mst2shift))*pow3(pow2(Mst1) - pow2(Mst2shift)) + truncateXt*
        pow2(Xt2)*(pow4(Mst1) - pow4(Mst2shift))))/(pow2(Mst1)*pow2(Mst2shift)*pow3(pow2(
        Mst1) - pow2(Mst2shift)));

      isDegen = std::abs(exactShifted - lim) >= std::abs(exact - lim)
         || !std::isfinite(exact) || !std::isfinite(exactShifted);
   }
   
   if(isDegen){
      shift = (32*(-3*(1 + lmMgl)*pow2(Mgl) + 5*(1 + lmMsq)*pow2(Msq) + (1 + lmMst1)*
        pow2(Mst1))*(-6*Xt2*pow2(Mst1) + truncateXt*pow2(Xt2) + 6*pow4(Mst1)))/
        (3.*pow6(Mst1));
   }
   else{
      shift = (16*(-6*(1 + lmMgl)*pow2(Mgl) + 10*(1 + lmMsq)*pow2(Msq) + (1 + lmMst1)*
        pow2(Mst1) + (1 + lmMst2)*pow2(Mst2))*(4*truncateXt*log(Mst2/Mst1)*
        pow2(Mst1)*pow2(Mst2)*pow2(Xt2) - 2*Xt2*pow3(pow2(Mst1) - pow2(Mst2)) +
        (pow2(Mst1) + pow2(Mst2))*pow3(pow2(Mst1) - pow2(Mst2)) + truncateXt*
        pow2(Xt2)*(pow4(Mst1) - pow4(Mst2))))/(pow2(Mst1)*pow2(Mst2)*pow3(pow2(
        Mst1) - pow2(Mst2)));
   }

   return shift;
}


/**
 * 	Sorts the eigenvalues of a 2x2 matrix.
 * 	@param es the EigenSolver object corresponding to the matrix whose eigenvalues should be sorted.
 * 	@return A sorted vector with the lowest eigenvalue at position 0.
 */
std::vector<double> himalaya::HierarchyCalculator::sortEigenvalues(const Eigen::EigenSolver<Eigen::Matrix2d>& es){
  std::vector<double> sortedEigenvalues = {sqrt(std::real(es.eigenvalues()(0))), sqrt(std::real(es.eigenvalues()(1)))};
  std::sort(sortedEigenvalues.begin(), sortedEigenvalues.end());
  return sortedEigenvalues;
}

/**
 * 	Calculates the loop corrected Higgs mass matrix at the order O(alpha_x). Here, x can be t or b.
 * 	@param ho a HierarchyObject with constant isAlphab.
 * 	@param shiftOneLoop An integer flag which is 0 or 1 in order to shift the one-loop terms to the MDR scheme.
 * 	@param shiftTwoLoop An integer flag which is 0 or 1 in order to shift the two-loop terms to the MDR scheme.
 * 	@return The loop corrected Higgs mass matrix at the order O(alpha_x).
 */
Eigen::Matrix2d himalaya::HierarchyCalculator::getMt41L(const himalaya::HierarchyObject& ho, const unsigned int shiftOneLoop, const unsigned int shiftTwoLoop){
   Eigen::Matrix2d Mt41L;
   const double GF = 1/(sqrt(2) * (pow2(p.vu) + pow2(p.vd)));
   const double beta = atan(p.vu/p.vd);
   const double Mst1 = shiftMst1ToMDR(ho, shiftOneLoop, shiftTwoLoop);
   const double Mst2 = shiftMst2ToMDR(ho, shiftOneLoop, shiftTwoLoop);
   const double sbeta = sin(beta);
   const double cbeta = cos(beta);
   double Mt;
   double s2t;
   if(!ho.getIsAlphab()){
      s2t = p.s2t;
      Mt = p.Mt;
   }
   else{
      s2t = p.s2b;
      Mt = p.Mb;
   }
   
   Mt41L(0, 0) = (-3 * GF * pow2(Mt) * pow2(p.mu) * pow2(1 / sbeta) *
      (-pow2(Mst1) + pow2(Mst2) + pow2(Mst1) * log(Mst1) +
      pow2(Mst2) * log(Mst1) - pow2(Mst1) * log(Mst2) -
      pow2(Mst2) * log(Mst2)) * pow2(s2t)) /
      (4. * sqrt(2) * (pow2(Mst1) - pow2(Mst2)) * pow2(Pi));
   
   Mt41L(0, 1) = (3 * GF * pow2(1 / sbeta) *
      (-(pow3(Mt) * p.mu * (log(Mst1) - log(Mst2)) * s2t) / 2. +
      (pow2(Mt) * pow2(p.mu) * 1 / tan(beta) *
      (-pow2(Mst1) + pow2(Mst2) + pow2(Mst1) * log(Mst1) +
      pow2(Mst2) * log(Mst1) - pow2(Mst1) * log(Mst2) -
      pow2(Mst2) * log(Mst2)) * pow2(s2t)) /
      (4. * (pow2(Mst1) - pow2(Mst2))) +
      (Mt * p.mu * (-pow2(Mst1) + pow2(Mst2) + pow2(Mst1) * log(Mst1) +
      pow2(Mst2) * log(Mst1) - pow2(Mst1) * log(Mst2) -
      pow2(Mst2) * log(Mst2)) * pow3(s2t)) / 8.)) /
      (sqrt(2) * pow2(Pi));
   
   Mt41L (1,0) = Mt41L(0,1);
   
   Mt41L(1, 1) = (3 * GF * pow2(1 / sbeta) *
      (pow4(Mt) * (log(Mst1) + log(Mst2) - 2 * log(Mt)) +
      pow3(Mt) * p.mu * 1 / tan(beta) * (log(Mst1) - log(Mst2)) * s2t +
      (pow2(Mt) * pow2(1 / sbeta) *
      (pow2(Mst1) * pow2(p.mu) * pow2(cbeta) -
      pow2(Mst2) * pow2(p.mu) * pow2(cbeta) -
      pow2(Mst1) * pow2(p.mu) * pow2(cbeta) * log(Mst1) -
      pow2(Mst2) * pow2(p.mu) * pow2(cbeta) * log(Mst1) +
      pow2(Mst1) * pow2(p.mu) * pow2(cbeta) * log(Mst2) +
      pow2(Mst2) * pow2(p.mu) * pow2(cbeta) * log(Mst2) +
      2 * pow4(Mst1) * log(Mst1) * pow2(sbeta) -
      4 * pow2(Mst1) * pow2(Mst2) * log(Mst1) * pow2(sbeta) +
      2 * pow4(Mst2) * log(Mst1) * pow2(sbeta) -
      2 * pow4(Mst1) * log(Mst2) * pow2(sbeta) +
      4 * pow2(Mst1) * pow2(Mst2) * log(Mst2) * pow2(sbeta) -
      2 * pow4(Mst2) * log(Mst2) * pow2(sbeta)) * pow2(s2t)) /
      (4. * (pow2(Mst1) - pow2(Mst2))) -
      (Mt * p.mu * 1 / tan(beta) * (-pow2(Mst1) + pow2(Mst2) +
      pow2(Mst1) * log(Mst1) + pow2(Mst2) * log(Mst1) -
      pow2(Mst1) * log(Mst2) - pow2(Mst2) * log(Mst2)) *
      pow3(s2t)) / 4. -
      ((pow2(Mst1) - pow2(Mst2)) *
      (-pow2(Mst1) + pow2(Mst2) + pow2(Mst1) * log(Mst1) +
      pow2(Mst2) * log(Mst1) - pow2(Mst1) * log(Mst2) -
      pow2(Mst2) * log(Mst2)) * pow4(s2t)) / 16.)) /
      (sqrt(2) * pow2(Pi));
   
    return Mt41L;
}

/**
 * 	Calculates the loop corrected Higgs mass matrix at the order O(alpha_x*alpha_s). Here, x can be t or b.
 * 	@param ho a HierarchyObject with constant isAlphab.
 * 	@param shiftOneLoop An integer flag which is 0 or 1 in order to shift the one-loop terms to the MDR scheme.
 * 	@param shiftTwoLoop An integer flag which is 0 or 1 in order to shift the two-loop terms to the MDR scheme.
 * 	@return The loop corrected Higgs mass matrix at the order O(alpha_x*alpha_s).
 */
Eigen::Matrix2d himalaya::HierarchyCalculator::getMt42L(const himalaya::HierarchyObject& ho,
							const unsigned int shiftOneLoop,
							const unsigned int shiftTwoLoop){
   Eigen::Matrix2d Mt42L;
   double S11, S12, S22;
   double Mt2;
   double MG = p.MG;
   double st;
   double ct;
   double Mst12 = pow2(shiftMst1ToMDR(ho, shiftOneLoop, shiftTwoLoop));
   double Mst22 = pow2(shiftMst2ToMDR(ho, shiftOneLoop, shiftTwoLoop));
   if(!ho.getIsAlphab()){
      const double theta = asin(p.s2t)/2.;
      Mt2 = pow2(p.Mt);
      st = sin(theta);
      ct = cos(theta);
   }
   else{
      const double theta = asin(p.s2b)/2.;
      Mt2 = pow2(p.Mb);
      st = sin(theta);
      ct = cos(theta);
   }
   double scale2 = pow2(p.scale);
   // note the sign difference in mu
   double mu = - p.mu;
   double tanb = p.vu/p.vd;
   double v2 = pow2(p.vu) + pow2(p.vd);
   double gs = p.g3;
   int os = 0;
   dszhiggs_(&Mt2, &MG, &Mst12, &Mst22, &st, &ct, &scale2, &mu, &tanb, &v2, &gs, &os, &S11, &S22, &S12);
   Mt42L(0, 0) = S11;
   Mt42L(1, 0) = S12;
   Mt42L(0, 1) = S12;
   Mt42L(1, 1) = S22;
   return Mt42L;
}

/**
 * 	Calculates the contribution to the order (alpha_x) and (alpha_s alpha_x) as the difference
 * 	of the Higgs mass matrices of the MDR and DR scheme. Here, x can be t or b.
 * 	@param ho a HierarchyObject with constant isAlphab.
 * 	@param shiftOneLoop a bool to shift the terms at one-loop level.
 * 	@param shiftTwoLoop a bool to shift the terms at two-loop level.
 * 	@return The loop corrected Higgs mass matrix difference of the MDR and DR scheme at the given order.
 */
Eigen::Matrix2d himalaya::HierarchyCalculator::calcDRbarToMDRbarShift(const himalaya::HierarchyObject& ho,
								      const bool shiftOneLoop,
								      const bool shiftTwoLoop){
   if(shiftOneLoop && shiftTwoLoop){
      return getMt41L(ho, 1, 1) + getMt42L(ho, 1, 1) - getMt41L(ho, 0, 0) - getMt42L(ho, 0, 0);
   }
   else if(shiftOneLoop){
      return getMt41L(ho, 1, 1) - getMt41L(ho, 0, 0);
   }
   else if(shiftTwoLoop){
      return getMt42L(ho, 1, 1) - getMt42L(ho, 0, 0);
   }
   else{
      return Eigen::Matrix2d::Zero();
   }
}


/**
 * 	Estimates the uncertainty of the expansion at a given order.
 * 	@param ho a HierarchyObject with constant isAlphab.
 * 	@param massMatrix the CP-even Higgs mass matrix without the corrections whose uncertainty should be estimated.
 * 	@param oneLoopFlag an integer flag which is 0 or 1 in order to estimate the uncertainty of the one-loop expansion terms.
 * 	@param twoLoopFlag an integer flag which is 0 or 1 in order to estimate the uncertainty of the two-loop expansion terms.
 * 	@param threeLoopFlag an integer flag which is 0 or 1 in order to estimte the uncertainty of the three-loop expansion terms.
 * 	@return A double which is the estimated uncertainty.
 */
double himalaya::HierarchyCalculator::getExpansionUncertainty(himalaya::HierarchyObject& ho,
							      const Eigen::Matrix2d& massMatrix,
							      const unsigned int oneLoopFlag,
							      const unsigned int twoLoopFlag,
							      const unsigned int threeLoopFlag){
   double Mh;
   double Mhcut;
   std::vector<double> errors;
   // reset flags
   flagMap.at(ExpansionDepth::xxMst) = 1;
   Eigen::EigenSolver<Eigen::Matrix2d> es;
   switch (getCorrectHierarchy(ho.getSuitableHierarchy())) {
   case Hierarchies::h3:
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mh = sortEigenvalues(es).at(0);
      // truncate the expansion at all variables with one order lower than the expansion depth and evaluate the expansion uncertainty 
      flagMap.at(ExpansionDepth::xxDmglst1) = 0;
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mhcut = sortEigenvalues(es).at(0);
      errors.push_back(std::abs(Mh - Mhcut));
      flagMap.at(ExpansionDepth::xxDmglst1) = 1;
      flagMap.at(ExpansionDepth::xxDmsqst1) = 0;
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mhcut = sortEigenvalues(es).at(0);
      errors.push_back(std::abs(Mh - Mhcut));
      flagMap.at(ExpansionDepth::xxDmsqst1) = 1;
      flagMap.at(ExpansionDepth::xxDmst12) = 0;
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mhcut = sortEigenvalues(es).at(0);
      errors.push_back(std::abs(Mh - Mhcut));
      flagMap.at(ExpansionDepth::xxDmst12) = 1;
      break;
   case Hierarchies::h4:
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mh = sortEigenvalues(es).at(0);
      flagMap.at(ExpansionDepth::xxAt) = 0;
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mhcut = sortEigenvalues(es).at(0);
      errors.push_back(std::abs(Mh - Mhcut));
      flagMap.at(ExpansionDepth::xxAt) = 1;
      flagMap.at(ExpansionDepth::xxlmMsusy) = 0;
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mhcut = sortEigenvalues(es).at(0);
      errors.push_back(std::abs(Mh - Mhcut));
      flagMap.at(ExpansionDepth::xxlmMsusy) = 1;
      flagMap.at(ExpansionDepth::xxMsq) = 0;
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mhcut = sortEigenvalues(es).at(0);
      errors.push_back(std::abs(Mh - Mhcut));
      flagMap.at(ExpansionDepth::xxMsq) = 1;
      flagMap.at(ExpansionDepth::xxMsusy) = 0;
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mhcut = sortEigenvalues(es).at(0);
      errors.push_back(std::abs(Mh - Mhcut));
      flagMap.at(ExpansionDepth::xxMsusy) = 1;
      break;
   case Hierarchies::h5:
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mh = sortEigenvalues(es).at(0);
      flagMap.at(ExpansionDepth::xxDmglst1) = 0;
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mhcut = sortEigenvalues(es).at(0);
      errors.push_back(std::abs(Mh - Mhcut));
      flagMap.at(ExpansionDepth::xxDmglst1) = 1;
      flagMap.at(ExpansionDepth::xxMsq) = 0;
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mhcut = sortEigenvalues(es).at(0);
      errors.push_back(std::abs(Mh - Mhcut));
      flagMap.at(ExpansionDepth::xxMsq) = 1;
      break;
   case Hierarchies::h6:
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mh = sortEigenvalues(es).at(0);
      flagMap.at(ExpansionDepth::xxDmglst2) = 0;
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mhcut = sortEigenvalues(es).at(0);
      errors.push_back(std::abs(Mh - Mhcut));
      flagMap.at(ExpansionDepth::xxDmglst2) = 1;
      flagMap.at(ExpansionDepth::xxMsq) = 0;
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mhcut = sortEigenvalues(es).at(0);
      errors.push_back(std::abs(Mh - Mhcut));
      flagMap.at(ExpansionDepth::xxMsq) = 1;
      break;
   case Hierarchies::h6b:
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mh = sortEigenvalues(es).at(0);
      flagMap.at(ExpansionDepth::xxDmglst2) = 0;
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mhcut = sortEigenvalues(es).at(0);
      errors.push_back(std::abs(Mh - Mhcut));
      flagMap.at(ExpansionDepth::xxDmglst2) = 1;
      flagMap.at(ExpansionDepth::xxDmsqst2) = 0;
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mhcut = sortEigenvalues(es).at(0);
      errors.push_back(std::abs(Mh - Mhcut));
      flagMap.at(ExpansionDepth::xxDmsqst2) = 1;
      break;
   case Hierarchies::h9:
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mh = sortEigenvalues(es).at(0);
      flagMap.at(ExpansionDepth::xxDmsqst1) = 0;
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mhcut = sortEigenvalues(es).at(0);
      errors.push_back(std::abs(Mh - Mhcut));
      flagMap.at(ExpansionDepth::xxDmsqst1) = 1;
      flagMap.at(ExpansionDepth::xxDmst12) = 0;
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mhcut = sortEigenvalues(es).at(0);
      errors.push_back(std::abs(Mh - Mhcut));
      flagMap.at(ExpansionDepth::xxDmst12) = 1;
      flagMap.at(ExpansionDepth::xxMgl) = 0;
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mhcut = sortEigenvalues(es).at(0);
      errors.push_back(std::abs(Mh - Mhcut));
      flagMap.at(ExpansionDepth::xxMgl) = 1;
      break;
   }
   // evalue the sqrt of the squared errors
   const double squaredErrorSum =
      std::accumulate(errors.cbegin(), errors.cend(), 0.,
                      [](double l, double r) { return l + r*r; });

   // set the expansion depth for the next comparison
   flagMap.at(ExpansionDepth::xxMst) = 0;
   flagMap.at(ExpansionDepth::xx) = 0;

   return std::sqrt(squaredErrorSum);
}

/**
 * 	Maps a hierarchy to it's mother hierarchy.
 * 	@param hierarchy the key to a hierarchy.
 * 	@throws runtime_error Throws a runtime_error if the given hierarchy is not included.
 * 	@returns The key of the mother hierarchy.
 */
int himalaya::HierarchyCalculator::getCorrectHierarchy(const int hierarchy){
   if(hierarchy < 0 || hierarchy > 13){
      if(hierarchy == -1){
	 throw std::runtime_error("\033[1;31m Error: No suitable hierarchy found!\033[0m");
      } else{
	 throw std::runtime_error("\033[1;31m Error: Hierarchy " + std::to_string(hierarchy) + " not included!\033[0m"); 
      }
   }
   return hierarchyMap.at(hierarchy);
}

/**
 * 	Prints out some information about Himalaya.
 */
void himalaya::HierarchyCalculator::printInfo(){
   std::cout << "............................................................................\n";
   std::cout << "Himalaya " << Himalaya_VERSION_MAJOR << "." << Himalaya_VERSION_MINOR << "." << Himalaya_VERSION_RELEASE << "\tѧѦѧ \n";
   std::cout << "Uses code by P. Slavich et al. (2-loop αt*αs) [hep-ph/0105096].\n";
   std::cout << "Uses the 3-loop αt*αs^2 contributions of Kant et al. [arXiv:1005.5709],\n";
   std::cout << "  2-loop αt*αs contributions of Bagnaschi et.al. [arXiv:1407.4081],\n";
   std::cout << "  2-loop αs^2 contributions Bednyakov et.al. [hep-ph/0210258,hep-ph/0507139]\n";
   std::cout << "............................................................................\n";
}
