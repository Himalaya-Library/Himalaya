#include "H4.hpp"
#include "HierarchyCalculator.hpp"
#include "Constants.hpp"
#include "Utils.hpp"
#include <cmath>
#include <type_traits>

/**
 * 	Constructor
 * 	@param flagMap the flagMap for the truncation of expansion variables
 * 	@param Al4p a double alpha_s/4/Pi
 * 	@param At a double tri-linear breaking term
 * 	@param beta a double which is the mixing angle beta
 * 	@param lmMt a double log((<renormalization scale> / Mt)^2)
 * 	@param lmMsq a double log((<renormalization scale> / Msq)^2)
 * 	@param lmMst1 a double log((<renormalization scale> / Mst1)^2)
 * 	@param Mt a double top/bottom quark mass
 * 	@param Msusy a double (Mst1 + Mst2 + Mgl) / 3.
 * 	@param Msq a double the average squark mass w/o the top squark
 * 	@param mdrFlag an int 0 for DR and 1 for MDR scheme
 * 	@param oneLoopFlag an int flag to consider the one-loop expansion terms
 * 	@param twoLoopFlag an int flag to consider the two-loop expansion terms
 * 	@param threeLoopFlag an int flag to consider the three-loop expansion terms
 */
himalaya::H4::H4(const std::map<unsigned int, unsigned int>& flagMap, double Al4p, double At, double beta,
		 double lmMt, double lmMsq, double lmMsusy, double Mt, double Msusy, double Msq,
		 int mdrFlag, int oneLoopFlag, int twoLoopFlag, int threeLoopFlag){
   // abbrev for cos(beta) and sin(beta)
   Cbeta = cos(beta);
   Sbeta = sin(beta);
   this -> At = At;
   this -> lmMt = lmMt;
   this -> lmMsq = lmMsq;
   this -> lmMsusy = lmMsusy;
   this -> Mt = Mt;
   this -> Msusy = Msusy;
   this -> Msq = Msq;
   // mdr flags, indicates if one wants to shift the dr stop mass to the mdr stop mass
   shiftst1 = mdrFlag;
   shiftst2 = mdrFlag;
   shiftst3 = mdrFlag;
   // expansion flags
   xAt = flagMap.at(HierarchyCalculator::xxAt);
   xMsq = flagMap.at(HierarchyCalculator::xxMsq);
   xlmMsusy = flagMap.at(HierarchyCalculator::xxlmMsusy);
   xMsusy = flagMap.at(HierarchyCalculator::xxMsusy);
   
   s1 = 
   #include "../hierarchies/h4/sigS1Full.inc"
   ;
   s2 = 
   #include "../hierarchies/h4/sigS2Full.inc"
   ;
   s12 = 
   #include "../hierarchies/h4/sigS12Full.inc"
   ;
}

/**
 * 	@return The diagonal (1, 1) matrix element of the Higgs mass matrix as a double for the hierarchy 'H4'
 */
double himalaya::H4::getS1() const {
   return s1;
}

/**
 * 	@return The diagonal (2, 2) matrix element of the Higgs mass matrix as a double for the hierarchy 'H4'
 */
double himalaya::H4::getS2() const {
   return s2;
}

/**
 * 	@return The off-diagonal (1, 2) = (2, 1) matrix element of the Higgs mass matrix as a double for the hierarchy 'H4'
 */
double himalaya::H4::getS12() const {
   return s12;
}

/**
 * 	@return returns the susy log^0 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double himalaya::H4::calc_coef_at_as2_no_sm_logs_log0() const {

   const double result =
      (-(pow2(Sbeta)*pow4(Mt)*(-591666768*(-10589 + 7500*z2)*pow4(Msusy)*pow6(
        Msq) - 1724976*(-2819419 + 1800750*z2)*pow4(Msq)*pow6(Msusy) -
        14791669200*(-691 + 270*z2 - 6*z3)*pow2(Msusy)*pow8(Msq) +
        221875038000*pow2(Msusy)*pow3(log(pow2(Msq)/pow2(Msusy)))*pow8(Msq) -
        665500*(-6262157 + 4000752*z2)*pow2(Msq)*pow8(Msusy) - 96049800*pow2(
        Msusy)*pow2(log(pow2(Msq)/pow2(Msusy)))*(14586*pow4(Msq)*pow4(Msusy) +
        20328*pow2(Msusy)*pow6(Msq) + 12760*pow2(Msq)*pow6(Msusy) + 23100*pow8(
        Msq) + 11865*pow8(Msusy)) + 1331250228000*(-1 + 2*z2)*power10(Msq) -
        55440*log(pow2(Msq)/pow2(Msusy))*(-28388052*pow4(Msusy)*pow6(Msq) -
        51750369*pow4(Msq)*pow6(Msusy) - 16008300*(-5 + 3*z2)*pow2(Msusy)*pow8(
        Msq) - 58536775*pow2(Msq)*pow8(Msusy) + 48024900*power10(Msq) -
        63123270*power10(Msusy)) - 5145*(-742606013 + 474368400*z2)*power10(
        Msusy)))/(4.992188355e10*pow2(Msusy)*pow8(Msq)))/pow4(Mt)/pow2(Sbeta)*
        12.;

   return result;
}

/**
 * 	@return returns the susy log^1 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double himalaya::H4::calc_coef_at_as2_no_sm_logs_log1() const {

   const double result =
      (-(pow2(Sbeta)*pow4(Mt)*(2160*pow4(Msusy)*pow6(Msq) + 540*pow4(Msq)*pow6(
        Msusy) - 16*(-173 + 135*z2 + 54*z3)*pow2(Msusy)*pow8(Msq) + 240*pow2(
        Msq)*pow8(Msusy) + 180*log(pow2(Msq)/pow2(Msusy))*pow2(Msusy)*(14*pow4(
        Msq)*pow4(Msusy) + 20*pow2(Msusy)*pow6(Msq) + 12*pow2(Msq)*pow6(Msusy)
        + 28*pow8(Msq) + 11*pow8(Msusy)) + 4320*power10(Msq) + 135*power10(
        Msusy)))/(81.*pow2(Msusy)*pow8(Msq)))/pow4(Mt)/pow2(Sbeta)*12.; 

   return result;
}

/**
 * 	@return returns the susy log^2 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double himalaya::H4::calc_coef_at_as2_no_sm_logs_log2() const {

   const double result =
      ((8*(221 + 45*log(pow2(Msq)/pow2(Msusy)))*pow2(Sbeta)*pow4(Mt))/27.)/
      pow4(Mt)/pow2(Sbeta)*12.; 

   return result;
}

/**
 * 	@return returns the susy log^3 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double himalaya::H4::calc_coef_at_as2_no_sm_logs_log3() const {

   const double result =
      ((-224*pow2(Sbeta)*pow4(Mt))/9.)/pow4(Mt)/pow2(Sbeta)*12.; 

   return result;
}


