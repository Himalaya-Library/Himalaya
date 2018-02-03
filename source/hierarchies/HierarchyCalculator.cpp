#include "HierarchyCalculator.hpp"
#include "Hierarchies.hpp"
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
#include "Utils.hpp"
#include <iomanip>
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
 * 	@param verbose a bool which suppresses the information of the calculation if set to flase
 */
himalaya::HierarchyCalculator::HierarchyCalculator(const Parameters& p_, const bool verbose_)
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

   //sw2
   const double sw2 = 1 - pow2(p.MW / p.MZ);

   // Al4p
   Al4p = pow2(p.g3 / (4 * Pi));

   // MGl
   Mgl = p.MG;

   // Msq, checked
   Msq = (2 * sqrt(p.mq2(0, 0)) + sqrt(p.mu2(0, 0)) + sqrt(p.md2(0, 0))	// sup and sdown
      + 2 * sqrt(p.mq2(1, 1)) + sqrt(p.mu2(1, 1)) + sqrt(p.md2(1, 1))	// scharm and sstrange
      // sbottom
      + sqrt(p.mq2(2, 2) + pow2(p.Mb) - (1 / 2. - 1 / 3. * sw2) * pow2(p.MZ) * cos(2 * beta))
      + sqrt(p.md2(2, 2) + pow2(p.Mb) - 1 / 3. * sw2 * pow2(p.MZ) * cos(2 * beta))) / 10.;
   
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
 *        @param mdrFlag an integer to choose among DR and MDR scheme. Set to DR by default.
 * 	@return A HierarchyObject which holds all information of the calculation.
 */
himalaya::HierarchyObject himalaya::HierarchyCalculator::calculateDMh3L(bool isAlphab, const int mdrFlag){
   HierarchyObject ho (isAlphab);
   
   // set mdrFlag
   ho.setMDRFlag(mdrFlag);
   
   // compare hierarchies and get the best fitting hierarchy
   compareHierarchies(ho);
   
   // calculate the DR to MDR shift with the obtained hierarchy
   if(mdrFlag == 1){
      ho.setDRToMDRShift(calcDRbarToMDRbarShift(ho, true, true));
   }
   else{
      ho.setDRToMDRShift(calcDRbarToMDRbarShift(ho, false, false));
   }
   
   // calculate the 3-loop Higgs mass matrix for the obtained hierarhy
   ho.setDMh(3, calculateHierarchy(ho, 0, 0, 1));
   
   // set the alpha_x contributions
   ho.setDMh(1, getMt41L(ho, mdrFlag, mdrFlag));
   
   // set the alpha_x*alpha_s contributions
   ho.setDMh(2, getMt42L(ho, mdrFlag, mdrFlag));
   
   // estimate the uncertainty of the expansion at 3-loop level
   ho.setExpUncertainty(3, getExpansionUncertainty(ho,
						   ho.getDMh(0) + ho.getDMh(1) + ho.getDMh(2), 0, 0, 1));
   
   // set the uncertainty of the expansion at 1-loop level to 0 by convention, if the user needs this value getExpansionUncertainty should be called
   ho.setExpUncertainty(1, 0.);

   double Xt = p.At - p.mu * p.vd / p.vu;
   double mQ3 = sqrt(p.mq2(2,2));
   double mU3 = sqrt(p.mu2(2,2));
   double m3 = p.MG;
   double msq = Msq;
   
   // take care of degenerated masses
   checkForDegenerateCase(mQ3, mU3, m3, msq, Xt);
   
   // calc the delta zeta
   himalaya::mh2_eft::Mh2EFTCalculator mh2EFTCalculator;
   //const double lmMQ3 = log(pow2(p.scale / mQ3));
   const double lmMst1 = log(pow2(p.scale / p.MSt(0)));
   const double zeta_3L_const =
      mh2EFTCalculator.coeff_as2_susy_log02(mQ3, mU3, m3, msq, p.MSt(0), Xt);
   const double prefactor = 1./16.;

   // Factorization: Himalaya_const + Log(mu^2/M_X^2) * Himalaya_coeff_log^1 + Log(mu^2/M_X^2)^2 Himalaya_coeff_log^2 
   //	+ Log(mu^2/M_X^2)^3 Himalaya_coeff_log^3 - EFT_const_w/o_dlatas2_and_Log(M_X^2/M_Y^2)
   // the factor 1/16 is a partial loop factor which gets factorized into zeta
   // M_X is a susy mass
   ho.setZetaHimalaya(prefactor * (ho.getZetaHimalaya() - zeta_3L_const));
   
   // add the EFT logs and subtract constant part twice to avoid double counting
   // Factorization: Himalaya_const - EFT_const_w/o_dlatas2_and_Log(M_X^2/M_Y^2)
   //	+ Log(mu^2/mst1^2)^1 EFT_coeff_log^1  + Log(mu^2/mst1^2)^2 EFT_coeff_log^2 + Log(mu^2/mst1^2)^3 EFT_coeff_log^3
   ho.setZetaEFT(prefactor * (ho.getZetaEFT() - zeta_3L_const
      + lmMst1 * mh2EFTCalculator.coeff_as2_susy_log12(mQ3, mU3, m3, msq, p.MSt(0), Xt)
      + pow2(lmMst1) * mh2EFTCalculator.coeff_as2_susy_log22(mQ3, mU3, m3, msq, p.MSt(0), Xt)
      + pow3(lmMst1) * mh2EFTCalculator.coeff_as2_susy_log32()));

   // Set zeta for degenerated mass case, where MS = mQ3, the argument Xt is here suppressed by the SUSY scale (named xt)
   ho.setZetaDegenerated(mh2EFTCalculator.getZetaDegenerated(p.scale, mQ3, Xt / mQ3));
   /*   std::cout << "xt " << Xt/mQ3 << " do const " << zeta_3L_const << " logs " <<
      lmMst1 * mh2EFTCalculator.coeff_as2_susy_log12(mQ3, mU3, m3, msq, p.MSt(0), Xt)
     + pow2(lmMst1) * mh2EFTCalculator.coeff_as2_susy_log22(mQ3, mU3, m3, msq, p.MSt(0), Xt)
     + pow3(lmMst1) * mh2EFTCalculator.coeff_as2_susy_log32() << "\n";*/

   // set zeta (non-logarithmic part)
   ho.setZetaConst(prefactor * (ho.getZetaConst() - zeta_3L_const));

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
	 
	 //DEPRECATED calc 1-loop shift for DRbar -> MDRbar
	 //calc difference of Mt41L or Mt41L in the MDRbar scheme directly
	 //it seems that in H3m the sign of the function getShift is wrong (Mt4LDRbar - Mt4LMDRbar)????
	 //to be consistent in the MDRbar-scheme we should subtract Mt4LDRbar, shouldn't we?
	 //Eigen::Matrix2d shift = getShift(hierarchyMap.at(hierarchy), isAlphab);
	 
	 //calculate the exact Higgs mass at 2-loop (only up to alpha_s alpha_t/b)
	 const Eigen::EigenSolver<Eigen::Matrix2d> es2L (treelvl + Mt41L + Mt42L);
	 const double Mh2l = sortEigenvalues(es2L).at(0);

	 // calculate the expanded 2-loop expression with the specific hierarchy
	 const Eigen::EigenSolver<Eigen::Matrix2d> esExpanded (treelvl + Mt41L + calculateHierarchy(ho, 0, 1, 0));
	 
	 // calculate the higgs mass in the given mass hierarchy and compare the result to estimate the error
	 const double Mh2LExpanded = sortEigenvalues(esExpanded).at(0);

	 // estimate the error
	 const double twoLoopError = fabs((Mh2l - Mh2LExpanded));

	 // estimate the uncertainty of the expansion
	 const double expUncertainty = getExpansionUncertainty(ho, treelvl + Mt41L, 0, 1, 0);
	 
	 // add these errors to include the error of the expansion in the comparison
	 const double currError = sqrt(pow2(twoLoopError) + pow2(expUncertainty));

	 // if the error is negative, it is the first iteration and there is no hierarchy which fits better
	 if(error < 0){
	    error = currError;
	    suitableHierarchy = hierarchy;
	    ho.setAbsDiff2L(twoLoopError);
	    ho.setRelDiff2L(twoLoopError/Mh2l);
	    ho.setExpUncertainty(2, expUncertainty);
	 }
	 // compare the current error with the last error and choose the hierarchy which fits best (lowest error)
	 else if(currError < error){
	    error = currError;
	    suitableHierarchy = hierarchy;
	    ho.setAbsDiff2L(twoLoopError);
	    ho.setRelDiff2L(twoLoopError/Mh2l);
	    ho.setExpUncertainty(2, expUncertainty);
	 }
      }
   }
   ho.setSuitableHierarchy(suitableHierarchy);
   // reset the flags
   flagMap.at(ExpansionDepth::xx) = 1;
   flagMap.at(ExpansionDepth::xxMst) = 1;
   return suitableHierarchy;
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
      At = p.At;
      Mt = p.Mt;
      s2t = p.s2t;
   }
   else {
      At = p.Ab;
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
			ho.setZetaHimalaya(c
			   + lmMst1 * hierarchy3.calc_coef_at_as2_no_sm_logs_log1()
			   + pow2(lmMst1) * hierarchy3.calc_coef_at_as2_no_sm_logs_log2()
			   + pow3(lmMst1) * hierarchy3.calc_coef_at_as2_no_sm_logs_log3());
			ho.setZetaEFT(c);
                        ho.setZetaConst(c);
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
			ho.setZetaHimalaya(c
			   + lmMst1 * hierarchy32q2g.calc_coef_at_as2_no_sm_logs_log1()
			   + pow2(lmMst1) * hierarchy32q2g.calc_coef_at_as2_no_sm_logs_log2()
			   + pow3(lmMst1) * hierarchy32q2g.calc_coef_at_as2_no_sm_logs_log3());
			ho.setZetaEFT(c);
                        ho.setZetaConst(c);
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
			ho.setZetaHimalaya(c
			   + lmMst1 * hierarchy3q22g.calc_coef_at_as2_no_sm_logs_log1()
			   + pow2(lmMst1) * hierarchy3q22g.calc_coef_at_as2_no_sm_logs_log2()
			   + pow3(lmMst1) * hierarchy3q22g.calc_coef_at_as2_no_sm_logs_log3());
			ho.setZetaEFT(c);
                        ho.setZetaConst(c);
		     }
		  }
		  break;
	       }
	    }
	    break;
	    case Hierarchies::h4:{
	       const double Msusy = (Mst1 + Mst2 + Mgl) / 3.;
	       const double lmMsusy = log(pow2(p.scale / Msusy));
	       const H4 hierarchy4(flagMap, Al4p, At, beta,
		  lmMt, lmMsq, lmMsusy, Mt, Msusy, Msq,
		  ho.getMDRFlag(), oneLoopFlag, twoLoopFlag, threeLoopFlag);
	       curSig1 = hierarchy4.getS1();
	       curSig2 = hierarchy4.getS2();
	       curSig12 = hierarchy4.getS12();
	       if(oneLoopFlagIn == 0 && twoLoopFlagIn == 0 && threeLoopFlagIn == 1){
                  const double c = hierarchy4.calc_coef_at_as2_no_sm_logs_log0();
		  ho.setZetaHimalaya(c
		     + lmMsusy * hierarchy4.calc_coef_at_as2_no_sm_logs_log1()
		     + pow2(lmMsusy) * hierarchy4.calc_coef_at_as2_no_sm_logs_log2()
		     + pow3(lmMsusy) * hierarchy4.calc_coef_at_as2_no_sm_logs_log3());
		  ho.setZetaEFT(c);
                  ho.setZetaConst(c);
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
			ho.setZetaHimalaya(c
			   + lmMst1 * hierarchy5.calc_coef_at_as2_no_sm_logs_log1()
			   + pow2(lmMst1) * hierarchy5.calc_coef_at_as2_no_sm_logs_log2()
			   + pow3(lmMst1) * hierarchy5.calc_coef_at_as2_no_sm_logs_log3());
			ho.setZetaEFT(c);
                        ho.setZetaConst(c);
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
			ho.setZetaHimalaya(c
			   + lmMst1 * hierarchy5g1.calc_coef_at_as2_no_sm_logs_log1()
			   + pow2(lmMst1) * hierarchy5g1.calc_coef_at_as2_no_sm_logs_log2()
			   + pow3(lmMst1) * hierarchy5g1.calc_coef_at_as2_no_sm_logs_log3());
			ho.setZetaEFT(c);
                        ho.setZetaConst(c);
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
			ho.setZetaHimalaya(c
			   + lmMst1 * hierarchy6.calc_coef_at_as2_no_sm_logs_log1()
			   + pow2(lmMst1) * hierarchy6.calc_coef_at_as2_no_sm_logs_log2()
			   + pow3(lmMst1) * hierarchy6.calc_coef_at_as2_no_sm_logs_log3());
			ho.setZetaEFT(c);
                        ho.setZetaConst(c);
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
			ho.setZetaHimalaya(c
			   + lmMst1 * hierarchy6g2.calc_coef_at_as2_no_sm_logs_log1()
			   + pow2(lmMst1) * hierarchy6g2.calc_coef_at_as2_no_sm_logs_log2()
			   + pow3(lmMst1) * hierarchy6g2.calc_coef_at_as2_no_sm_logs_log3());
			ho.setZetaEFT(c);
                        ho.setZetaConst(c);
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
			ho.setZetaHimalaya(c
			   + lmMst1 * hierarchy6b.calc_coef_at_as2_no_sm_logs_log1()
			   + pow2(lmMst1) * hierarchy6b.calc_coef_at_as2_no_sm_logs_log2()
			   + pow3(lmMst1) * hierarchy6b.calc_coef_at_as2_no_sm_logs_log3());
			ho.setZetaEFT(c);
                        ho.setZetaConst(c);
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
			ho.setZetaHimalaya(c
			   + lmMst1 * hierarchy6b2qg2.calc_coef_at_as2_no_sm_logs_log1()
			   + pow2(lmMst1) * hierarchy6b2qg2.calc_coef_at_as2_no_sm_logs_log2()
			   + pow3(lmMst1) * hierarchy6b2qg2.calc_coef_at_as2_no_sm_logs_log3());
			ho.setZetaEFT(c);
                        ho.setZetaConst(c);
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
			ho.setZetaHimalaya(c
			   + lmMst1 * hierarchy6bq22g.calc_coef_at_as2_no_sm_logs_log1()
			   + pow2(lmMst1) * hierarchy6bq22g.calc_coef_at_as2_no_sm_logs_log2()
			   + pow3(lmMst1) * hierarchy6bq22g.calc_coef_at_as2_no_sm_logs_log3());
			ho.setZetaEFT(c);
                        ho.setZetaConst(c);
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
			ho.setZetaHimalaya(c
			   + lmMst1 * hierarchy6bq2g2.calc_coef_at_as2_no_sm_logs_log1()
			   + pow2(lmMst1) * hierarchy6bq2g2.calc_coef_at_as2_no_sm_logs_log2()
			   + pow3(lmMst1) * hierarchy6bq2g2.calc_coef_at_as2_no_sm_logs_log3());
			ho.setZetaEFT(c);
                        ho.setZetaConst(c);
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
			ho.setZetaHimalaya(c
			   + lmMst1 * hierarchy9.calc_coef_at_as2_no_sm_logs_log1()
			   + pow2(lmMst1) * hierarchy9.calc_coef_at_as2_no_sm_logs_log2()
			   + pow3(lmMst1) * hierarchy9.calc_coef_at_as2_no_sm_logs_log3());
			ho.setZetaEFT(c);
                        ho.setZetaConst(c);
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
			ho.setZetaHimalaya(c
			   + lmMst1 * hierarchy9q2.calc_coef_at_as2_no_sm_logs_log1()
			   + pow2(lmMst1) * hierarchy9q2.calc_coef_at_as2_no_sm_logs_log2()
			   + pow3(lmMst1) * hierarchy9q2.calc_coef_at_as2_no_sm_logs_log3());
			ho.setZetaEFT(c);
                        ho.setZetaConst(c);
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
      Eigen::Matrix<double, 2, 1> mdrMasses;
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

/**
 * 	Shifts Msx1 according to the hierarchy to the MDR scheme.
 * 	@param ho a HierarchyObject with constant isAlphab.
 * 	@param oneLoopFlag an integer flag which is 0 or 1 in order to shift the order O(alpha_s).
 * 	@param twoLoopFlag an integer flag which is 0 or 1 in order to shift the order O(alpha_s^2).
 * 	@return A double which is the MDR sx_1 mass.
 */
double himalaya::HierarchyCalculator::shiftMst1ToMDR(const himalaya::HierarchyObject& ho, const unsigned int oneLoopFlag, const unsigned int twoLoopFlag) {
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
   const double mdr2mst1ka = (-8. * twoLoopFlag * pow2(Al4p) * (10 * pow2(Msq) * (-1 + 2 * lmMsq + 2 * z2) + pow2(Mst2) * (-1 + 2 * lmMst2 + 2 * z2))) / (3. * pow2(Mst1));
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
double himalaya::HierarchyCalculator::shiftMst2ToMDR(const himalaya::HierarchyObject& ho, const unsigned int oneLoopFlag, const unsigned int twoLoopFlag) {
   double Mst2mod = 0., Mst2;
   if(!ho.getIsAlphab()){
      Mst2 = p.MSt(1);
   }
   else{
      Mst2 = p.MSb(1);
   }
   const double Dmglst2 = Mgl - Mst2;
   const double mdr2mst2ka = (-80. * twoLoopFlag * pow2(Al4p) * pow2(Msq) * (-1 + 2 * lmMsq + 2 * z2)) / (3. * pow2(Mst2));
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
   Mt41L (0,0) = (-3*GF*pow2(Mt)*pow2(p.mu)*pow2(1/sin(beta))*
		  (-pow2(Mst1) + pow2(Mst2) + pow2(Mst1)*log(Mst1) + 
		   pow2(Mst2)*log(Mst1) - pow2(Mst1)*log(Mst2) - 
		   pow2(Mst2)*log(Mst2))*pow2(s2t))/
     (4.*sqrt(2)*(pow2(Mst1) - pow2(Mst2))*pow2(Pi));
   Mt41L (0,1) = (3*GF*pow2(1/sin(beta))*
		  (-(pow3(Mt)*p.mu*(log(Mst1) - log(Mst2))*s2t)/2. + 
		   (pow2(Mt)*pow2(p.mu)*1/tan(beta)*
		    (-pow2(Mst1) + pow2(Mst2) + pow2(Mst1)*log(Mst1) + 
		     pow2(Mst2)*log(Mst1) - pow2(Mst1)*log(Mst2) - 
		     pow2(Mst2)*log(Mst2))*pow2(s2t))/
		   (4.*(pow2(Mst1) - pow2(Mst2))) + 
		   (Mt*p.mu*(-pow2(Mst1) + pow2(Mst2) + pow2(Mst1)*log(Mst1) + 
			     pow2(Mst2)*log(Mst1) - pow2(Mst1)*log(Mst2) - 
			     pow2(Mst2)*log(Mst2))*pow3(s2t))/8.))/
     (sqrt(2)*pow2(Pi));
   Mt41L (1,0) = Mt41L(0,1);
   Mt41L (1,1) =  (3*GF*pow2(1/sin(beta))*
		   (pow4(Mt)*(log(Mst1) + log(Mst2) - 2*log(Mt)) + 
		    pow3(Mt)*p.mu*1/tan(beta)*(log(Mst1) - log(Mst2))*s2t + 
		    (pow2(Mt)*pow2(1/sin(beta))*
		     (pow2(Mst1)*pow2(p.mu)*pow2(cos(beta)) - 
		      pow2(Mst2)*pow2(p.mu)*pow2(cos(beta)) - 
		      pow2(Mst1)*pow2(p.mu)*pow2(cos(beta))*log(Mst1) - 
		      pow2(Mst2)*pow2(p.mu)*pow2(cos(beta))*log(Mst1) + 
		      pow2(Mst1)*pow2(p.mu)*pow2(cos(beta))*log(Mst2) + 
		      pow2(Mst2)*pow2(p.mu)*pow2(cos(beta))*log(Mst2) + 
		      2*pow4(Mst1)*log(Mst1)*pow2(sin(beta)) - 
		      4*pow2(Mst1)*pow2(Mst2)*log(Mst1)*pow2(sin(beta)) + 
		      2*pow4(Mst2)*log(Mst1)*pow2(sin(beta)) - 
		      2*pow4(Mst1)*log(Mst2)*pow2(sin(beta)) + 
		      4*pow2(Mst1)*pow2(Mst2)*log(Mst2)*pow2(sin(beta)) - 
		      2*pow4(Mst2)*log(Mst2)*pow2(sin(beta)))*pow2(s2t))/
		    (4.*(pow2(Mst1) - pow2(Mst2))) - 
		    (Mt*p.mu*1/tan(beta)*(-pow2(Mst1) + pow2(Mst2) + 
					  pow2(Mst1)*log(Mst1) + pow2(Mst2)*log(Mst1) - 
					  pow2(Mst1)*log(Mst2) - pow2(Mst2)*log(Mst2))*
		     pow3(s2t))/4. - 
		    ((pow2(Mst1) - pow2(Mst2))*
		     (-pow2(Mst1) + pow2(Mst2) + pow2(Mst1)*log(Mst1) + 
		      pow2(Mst2)*log(Mst1) - pow2(Mst1)*log(Mst2) - 
		      pow2(Mst2)*log(Mst2))*pow4(s2t))/16.))/
     (sqrt(2)*pow2(Pi));
    return Mt41L;
}

/**
 * 	@deprecated
 * 	Shifts the 1-loop terms to the MDRbar scheme.
 * 	@param ho a HierarchyObject with constant isAlphab and a hierarchy candidate.
 * 	@return A matrix which corresponds to the difference of the DR and MDR scheme.
 */
Eigen::Matrix2d himalaya::HierarchyCalculator::getShift(const himalaya::HierarchyObject& ho){
   Eigen::Matrix2d shift;
   const double GF = 1/(sqrt(2) * (pow2(p.vu) + pow2(p.vd)));
   double Mst1;
   double Mst2;
   double Mt;
   double s2t;
   const double beta = atan(p.vu/p.vd);
   double deltamst1, deltamst2;
   if(!ho.getIsAlphab()){
      Mst1 = p.MSt(0);
      Mst2 = p.MSt(1);
      s2t = p.s2t;
      Mt = p.Mt;
      deltamst1 = -(Mst1 - shiftMst1ToMDR(ho, 1, 0));
      deltamst2 = -(Mst2 - shiftMst2ToMDR(ho, 1, 0));
   }
   else{
      Mst1 = p.MSb(0);
      Mst2 = p.MSb(1);
      s2t = p.s2b;
      Mt = p.Mb;
      deltamst1 = -(Mst1 - shiftMst1ToMDR(ho, 1, 0));
      deltamst2 = -(Mst2 - shiftMst2ToMDR(ho, 1, 0));
   }
   shift(0, 0) = (3 * GF * (deltamst2 * Mst1 - deltamst1 * Mst2) * pow2(Mt) * pow2(p.mu) * pow2(1 / Pi) *
      pow2(1/sin(beta)) * pow2(1 / (pow2(Mst1) - pow2(Mst2))) * pow2(s2t) *
      (4 * (log(Mst1) - log(Mst2)) * pow2(Mst1) * pow2(Mst2) - pow4(Mst1) +
	 pow4(Mst2))) / (4. * sqrt(2) * Mst1 * Mst2);
   shift(0, 1) = (-3 * GF * Mt * p.mu * pow2(1 / Pi) * pow2(1/sin(beta)) *
      pow2(1 / (pow2(Mst1) - pow2(Mst2))) * s2t *
      (-(pow2(pow2(Mst1) - pow2(Mst2)) *
	    (4 * (-(deltamst2 * Mst1) + deltamst1 * Mst2) * pow2(Mt) +
	       (-2 * Mst1 * Mst2 * (deltamst1 * Mst1 + deltamst2 * Mst2) *
		  (log(Mst1) - log(Mst2)) +
		  (deltamst2 * Mst1 + deltamst1 * Mst2) * pow2(Mst1) -
		  (deltamst2 * Mst1 + deltamst1 * Mst2) * pow2(Mst2)) * pow2(s2t))) 
	 + 2 * (deltamst2 * Mst1 - deltamst1 * Mst2) * Mt * p.mu * 1/tan(beta) *
	 (4 * (log(Mst1) - log(Mst2)) * pow2(Mst1) * pow2(Mst2) - pow4(Mst1) +
	    pow4(Mst2)) * s2t)) / (8. * sqrt(2) * Mst1 * Mst2);
   shift(1, 0) = shift(0, 1);
   shift(1, 1) = (3 * GF * pow2(1 / Pi) * pow2(1/sin(beta)) *
      ((Mt * p.mu * 1/tan(beta) * (-(deltamst1 * Mst1) + deltamst2 * Mst2 +
	       2 * deltamst1 * Mst1 * log(Mst1) + 2 * deltamst2 * Mst2 * log(Mst1) -
	       2 * deltamst1 * Mst1 * log(Mst2) - 2 * deltamst2 * Mst2 * log(Mst2) -
	       (deltamst2 * pow2(Mst1)) / Mst2 + (deltamst1 * pow2(Mst2)) / Mst1) *
	    pow3(s2t)) / 4. + (pow2(Mt) * pow2(1 / (pow2(Mst1) - pow2(Mst2))) *
	    pow2(s2t) * (2 * deltamst2 * pow7(Mst1) -
	       2 * deltamst1 * pow6(Mst1) * Mst2 -
	       2 * deltamst2 * Mst1 * pow6(Mst2) + 2 * deltamst1 * pow7(Mst2) -
	       4 * deltamst1 * pow6(Mst1) * Mst2 * log(Mst1) +
	       4 * deltamst2 * Mst1 * pow6(Mst2) * log(Mst1) +
	       4 * deltamst1 * pow6(Mst1) * Mst2 * log(Mst2) -
	       4 * deltamst2 * Mst1 * pow6(Mst2) * log(Mst2) -
	       deltamst2 * pow5(Mst1) * pow2(p.mu) * pow2(1/tan(beta)) -
	       deltamst1 * pow5(Mst2) * pow2(p.mu) * pow2(1/tan(beta)) +
	       2 * deltamst2 * pow2(Mst2) *
	       (pow5(Mst1) * (-3 + 2 * log(Mst1) - 2 * log(Mst2)) +
		  2 * (log(Mst1) - log(Mst2)) * pow2(p.mu) * pow2(1/tan(beta)) *
		  pow3(Mst1)) - 2 * deltamst1 * pow2(Mst1) *
	       (pow5(Mst2) * (3 + 2 * log(Mst1) - 2 * log(Mst2)) +
		  2 * (log(Mst1) - log(Mst2)) * pow2(p.mu) * pow2(1/tan(beta)) *
		  pow3(Mst2)) + deltamst1 * Mst2 * pow2(p.mu) * pow2(1/tan(beta)) *
	       pow4(Mst1) + 6 * deltamst1 * pow3(Mst2) * pow4(Mst1) +
	       8 * deltamst1 * log(Mst1) * pow3(Mst2) * pow4(Mst1) -
	       8 * deltamst1 * log(Mst2) * pow3(Mst2) * pow4(Mst1) +
	       deltamst2 * Mst1 * pow2(p.mu) * pow2(1/tan(beta)) * pow4(Mst2) +
	       6 * deltamst2 * pow3(Mst1) * pow4(Mst2) -
	       8 * deltamst2 * log(Mst1) * pow3(Mst1) * pow4(Mst2) +
	       8 * deltamst2 * log(Mst2) * pow3(Mst1) * pow4(Mst2))) / (4. * Mst1 * Mst2) -
	 ((deltamst2 * Mst1 + deltamst1 * Mst2) * pow4(Mt)) / (Mst1 * Mst2) -
	 ((deltamst2 * pow5(Mst1) + deltamst1 * pow5(Mst2) -
	    4 * deltamst2 * pow2(Mst2) * pow3(Mst1) -
	    4 * deltamst1 * pow2(Mst1) * pow3(Mst2) +
	    3 * deltamst1 * Mst2 * pow4(Mst1) -
	    4 * deltamst1 * Mst2 * log(Mst1) * pow4(Mst1) +
	    4 * deltamst1 * Mst2 * log(Mst2) * pow4(Mst1) +
	    3 * deltamst2 * Mst1 * pow4(Mst2) +
	    4 * deltamst2 * Mst1 * log(Mst1) * pow4(Mst2) -
	    4 * deltamst2 * Mst1 * log(Mst2) * pow4(Mst2)) * pow4(s2t)) /
	 (16. * Mst1 * Mst2) + (-(deltamst1 / Mst1) + deltamst2 / Mst2) * p.mu *
	 1/tan(beta) * pow3(Mt) * s2t)) / sqrt(2);
   return shift;
}


/**
 * 	Calculates the loop corrected Higgs mass matrix at the order O(alpha_x*alpha_s). Here, x can be t or b.
 * 	@param ho a HierarchyObject with constant isAlphab.
 * 	@param shiftOneLoop An integer flag which is 0 or 1 in order to shift the one-loop terms to the MDR scheme.
 * 	@param shiftTwoLoop An integer flag which is 0 or 1 in order to shift the two-loop terms to the MDR scheme.
 * 	@return The loop corrected Higgs mass matrix at the order O(alpha_x*alpha_s).
 */
Eigen::Matrix2d himalaya::HierarchyCalculator::getMt42L(const himalaya::HierarchyObject& ho, const unsigned int shiftOneLoop, const unsigned int shiftTwoLoop){
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
Eigen::Matrix2d himalaya::HierarchyCalculator::calcDRbarToMDRbarShift(const himalaya::HierarchyObject& ho, const bool shiftOneLoop, const bool shiftTwoLoop){
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
double himalaya::HierarchyCalculator::getExpansionUncertainty(himalaya::HierarchyObject& ho, const Eigen::Matrix2d& massMatrix, const unsigned int oneLoopFlag, const unsigned int twoLoopFlag, const unsigned int threeLoopFlag){
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
      errors.push_back(fabs(Mh - Mhcut));
      flagMap.at(ExpansionDepth::xxDmglst1) = 1;
      flagMap.at(ExpansionDepth::xxDmsqst1) = 0;
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mhcut = sortEigenvalues(es).at(0);
      errors.push_back(fabs(Mh - Mhcut));
      flagMap.at(ExpansionDepth::xxDmsqst1) = 1;
      flagMap.at(ExpansionDepth::xxDmst12) = 0;
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mhcut = sortEigenvalues(es).at(0);
      errors.push_back(fabs(Mh - Mhcut));
      flagMap.at(ExpansionDepth::xxDmst12) = 1;
      break;
   case Hierarchies::h4:
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mh = sortEigenvalues(es).at(0);
      flagMap.at(ExpansionDepth::xxAt) = 0;
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mhcut = sortEigenvalues(es).at(0);
      errors.push_back(fabs(Mh - Mhcut));
      flagMap.at(ExpansionDepth::xxAt) = 1;
      flagMap.at(ExpansionDepth::xxlmMsusy) = 0;
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mhcut = sortEigenvalues(es).at(0);
      errors.push_back(fabs(Mh - Mhcut));
      flagMap.at(ExpansionDepth::xxlmMsusy) = 1;
      flagMap.at(ExpansionDepth::xxMsq) = 0;
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mhcut = sortEigenvalues(es).at(0);
      errors.push_back(fabs(Mh - Mhcut));
      flagMap.at(ExpansionDepth::xxMsq) = 1;
      flagMap.at(ExpansionDepth::xxMsusy) = 0;
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mhcut = sortEigenvalues(es).at(0);
      errors.push_back(fabs(Mh - Mhcut));
      flagMap.at(ExpansionDepth::xxMsusy) = 1;
      break;
   case Hierarchies::h5:
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mh = sortEigenvalues(es).at(0);
      flagMap.at(ExpansionDepth::xxDmglst1) = 0;
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mhcut = sortEigenvalues(es).at(0);
      errors.push_back(fabs(Mh - Mhcut));
      flagMap.at(ExpansionDepth::xxDmglst1) = 1;
      flagMap.at(ExpansionDepth::xxMsq) = 0;
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mhcut = sortEigenvalues(es).at(0);
      errors.push_back(fabs(Mh - Mhcut));
      flagMap.at(ExpansionDepth::xxMsq) = 1;
      break;
   case Hierarchies::h6:
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mh = sortEigenvalues(es).at(0);
      flagMap.at(ExpansionDepth::xxDmglst2) = 0;
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mhcut = sortEigenvalues(es).at(0);
      errors.push_back(fabs(Mh - Mhcut));
      flagMap.at(ExpansionDepth::xxDmglst2) = 1;
      flagMap.at(ExpansionDepth::xxMsq) = 0;
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mhcut = sortEigenvalues(es).at(0);
      errors.push_back(fabs(Mh - Mhcut));
      flagMap.at(ExpansionDepth::xxMsq) = 1;
      break;
   case Hierarchies::h6b:
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mh = sortEigenvalues(es).at(0);
      flagMap.at(ExpansionDepth::xxDmglst2) = 0;
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mhcut = sortEigenvalues(es).at(0);
      errors.push_back(fabs(Mh - Mhcut));
      flagMap.at(ExpansionDepth::xxDmglst2) = 1;
      flagMap.at(ExpansionDepth::xxDmsqst2) = 0;
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mhcut = sortEigenvalues(es).at(0);
      errors.push_back(fabs(Mh - Mhcut));
      flagMap.at(ExpansionDepth::xxDmsqst2) = 1;
      break;
   case Hierarchies::h9:
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mh = sortEigenvalues(es).at(0);
      flagMap.at(ExpansionDepth::xxDmsqst1) = 0;
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mhcut = sortEigenvalues(es).at(0);
      errors.push_back(fabs(Mh - Mhcut));
      flagMap.at(ExpansionDepth::xxDmsqst1) = 1;
      flagMap.at(ExpansionDepth::xxDmst12) = 0;
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mhcut = sortEigenvalues(es).at(0);
      errors.push_back(fabs(Mh - Mhcut));
      flagMap.at(ExpansionDepth::xxDmst12) = 1;
      flagMap.at(ExpansionDepth::xxMgl) = 0;
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mhcut = sortEigenvalues(es).at(0);
      errors.push_back(fabs(Mh - Mhcut));
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
 *      A function to check for a degenerated mass case and shifting the masses to avoid huge contributions
 *      @param mQ3 the left-handed third generation soft-breaking parameter
 *      @param mU3 the right-handed third generation soft-breaking parameters
 *      @param m3 the gluino mass
 *      @param msq the average squark mass
 *      @param Xt the stop mixing parameters
 */
void himalaya::HierarchyCalculator::checkForDegenerateCase(double &mQ3, double &mU3, double &m3, double &msq, double &Xt){
   const double eps = 0.00;
   bool changedMQ3 = false;
   bool changedM3 = false;
   bool changedMsq = false;
   bool changedXt = false;
   if(abs(mU3 - m3) < m3 * eps){ 
      m3 = m3 + m3 * eps;
      changedM3 = true;
   }
   if(abs(mU3 - msq) < msq * eps){
      msq = msq + msq * eps;
      changedMsq = true;
   }
   if(abs(mU3 - Xt) < abs(Xt) * eps){
      Xt = Xt + Xt * eps;
      changedXt = true;
   }
   if(abs(mQ3 - m3) < m3 * eps){
      m3 = m3 + m3 * eps;
      changedM3 = true;
   }
   if(abs(mQ3 - msq) < msq * eps){
      msq = msq + msq * eps;
      changedMsq = true;
   }
   if(abs(mQ3 - Xt) < abs(Xt) * eps){
      Xt = Xt + Xt * eps;
      changedXt = true;
   }
   if(abs(mQ3 - mU3) < mQ3 * eps){
      mQ3 = mQ3 + mQ3 * eps;
      changedMQ3 = true;
   }
   if(abs(m3 - msq) < msq * eps){
      msq = msq + msq * eps;
      changedMsq = true;
   }
   if(abs(m3 - Xt) < abs(Xt) * eps){
      Xt = Xt + Xt * eps;
      changedXt = true;
   }
   if(abs(msq - Xt) < abs(Xt) * eps){
      Xt = Xt + Xt * eps;
      changedXt = true;
   }
   if(verbose && changedMQ3) std::cout << "\033[1;34mHimalaya info:\033[0m Changed mQ3 to " << mQ3 << " GeV to avoid poles in the EFT result!\n";
   if(verbose && changedM3) std::cout << "\033[1;34mHimalaya info:\033[0m Changed MG to " << m3 << " GeV to avoid poles in the EFT result!\n";
   if(verbose && changedMsq) std::cout << "\033[1;34mHimalaya info:\033[0m Changed Msq to " << msq << " GeV to avoid poles in the EFT result!\n";
   if(verbose && changedXt) std::cout << "\033[1;34mHimalaya info:\033[0m Changed Xt to " << Xt << " GeV to avoid poles in the EFT result!\n";
}

/**
 * 	Maps a hierarchy to it's mother hierarchy.
 * 	@param hierarchy the key to a hierarchy.
 * 	@throws runtime_error Throws a runtime_error if the given hierarchy is not included.
 * 	@returns The key of the mother hierarchy.
 */
int himalaya::HierarchyCalculator::getCorrectHierarchy(const int hierarchy){
   if(hierarchy < 0 || hierarchy > 13){
      throw std::runtime_error("\033[1;31m Error: Hierarchy " + std::to_string(hierarchy) + " not included!\033[0m");
   }
   return hierarchyMap.at(hierarchy);
}

/**
 * 	Prints out some information about Himalaya.
 */
void himalaya::HierarchyCalculator::printInfo(){
   std::cout << "......................................................................" << "\n";
   std::cout << "Himalaya " << Himalaya_VERSION_MAJOR << "." << Himalaya_VERSION_MINOR << "." << Himalaya_VERSION_RELEASE << "\tѧѦѧ \n";
   std::cout << "Uses code by: P. Slavich et al. (2-loop at*as) [hep-ph/0105096]" << "\n";
   std::cout << "Uses the 3-loop at*as^2 contributions of Kant et al. [arXiv:1005.5709]" << "\n";
   std::cout << "......................................................................" << "\n";
}
