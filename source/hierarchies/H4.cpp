// ====================================================================
// This file is part of Himalaya.
//
// Himalaya is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================

#include "H4.hpp"
#include "Constants.hpp"
#include "HimalayaFlags.hpp"
#include "Powers.hpp"
#include <cmath>

namespace himalaya {
namespace hierarchies {

/**
 * Constructor
 * @param expansionDepth the flagMap for the truncation of expansion variables
 * @param Al4p a double alpha_s/4/Pi
 * @param At a double tri-linear breaking term
 * @param beta a double which is the mixing angle beta
 * @param lmMt a double log((<renormalization scale> / Mt)^2)
 * @param lmMsq a double log((<renormalization scale> / Msq)^2)
 * @param lmMsusy a double log((<renormalization scale> / Msusy)^2)
 * @param Mt a double top/bottom quark mass
 * @param Msusy a double (Mst1 + Mst2 + Mgl) / 3.
 * @param Msq a double the average squark mass w/o the top squark
 * @param mdrFlag an int 0 for DR and 1 for MDR scheme
 * @param oneLoopFlag an int flag to consider the one-loop expansion terms
 * @param twoLoopFlag an int flag to consider the two-loop expansion terms
 * @param threeLoopFlag an int flag to consider the three-loop expansion terms
 */
H4::H4(const ExpansionFlags_t& expansionDepth, double Al4p, double At, double beta,
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
   this -> oneLoopFlag = oneLoopFlag;
   this -> twoLoopFlag = twoLoopFlag;
   this -> threeLoopFlag = threeLoopFlag;
   this -> Al4p = Al4p;
   // expansion flags
   xAt = expansionDepth.at(ExpansionDepth::At);
   xMsq = expansionDepth.at(ExpansionDepth::Msq);
   xlmMsusy = expansionDepth.at(ExpansionDepth::lmMsusy);
   xMsusy = expansionDepth.at(ExpansionDepth::Msusy);
}

/**
 * @return The diagonal (1, 1) matrix element of the Higgs mass matrix as a double for the hierarchy 'H4'
 */
double H4::getS1() const {
   return (2*threeLoopFlag*xAt*pow2(Al4p)*pow2(At)*(349 - 56*lmMsusy + 24*lmMt -
        282*z3 - 32*pow2(lmMsusy))*pow4(Mt)*(pow12(Sbeta) + pow2(Cbeta)*(pow2(
        Sbeta) + pow4(Sbeta) + pow6(Sbeta) + pow8(Sbeta) + power10(Sbeta))))/(
        27.*pow2(Cbeta)*pow2(Msusy));
}

/**
 * @return The diagonal (2, 2) matrix element of the Higgs mass matrix as a double for the hierarchy 'H4'
 */
double H4::getS2() const {
   return -(pow4(Mt)*((lmMsusy - lmMt)*oneLoopFlag + (8*Al4p*twoLoopFlag*(At*(1 +
        lmMsusy + lmMt) + Msusy*(1 + lmMsusy - lmMt + 2*lmMsusy*lmMt + pow2(
        lmMsusy) - 3*pow2(lmMt))))/(3.*Msusy) + threeLoopFlag*pow2(Al4p)*(
        204.74074074074073 + (800*lmMsq)/9. + (4*(1 - 2*lmMsusy)*shiftst3)/3. -
        (400*pow2(lmMsq))/9. + (8*lmMsusy*(-554 + 270*lmMsq + 135*pow2(lmMsq)))
        /81. - (1288*pow2(lmMsusy))/27. - (8*lmMt*(125 - 71*lmMsusy + 30*lmMsq*
        (-7 + 3*lmMsusy) + 39*pow2(lmMsusy)))/27. + (8*(-18 + 5*lmMsq + 23*
        lmMsusy)*pow2(lmMt))/3. + (8*(10589 + lmMsusy*(4910 - 3750*lmMt) -
        2250*lmMt + 10*lmMsq*(-266 + 285*lmMsusy + 375*lmMt) - 3300*pow2(lmMsq)
        + 450*pow2(lmMsusy))*pow2(Msusy))/(675.*pow2(Msq)) + (2*(-636*At*Msusy*
        z3 + xAt*pow2(At)*(-349 + 56*lmMsusy - 24*lmMt + 282*z3 + 32*pow2(
        lmMsusy)) - 180*(-1 + 2*lmMsq)*(-2 + shiftst1 + shiftst2)*xMsq*pow2(
        Msq) - 24*(-1 + 6*lmMsusy - 6*lmMt)*z3*pow2(Msusy)))/(27.*pow2(Msusy))
        - (40*pow3(lmMsq))/9. + 16*xlmMsusy*pow3(lmMsusy) - (184*pow3(lmMt))/3.
         - (4*(-85750*At*(-33 + 2*lmMsusy*(-11 + 6*lmMt) - 2*lmMsq*(-11 + 6*
        lmMsusy + 6*lmMt) + 12*pow2(lmMsq)) + Msusy*(-5638838 + 385875*lmMt -
        70*lmMsq*(-47521 + 20685*lmMsusy + 25725*lmMt) + 35*lmMsusy*(-106067 +
        51450*lmMt) + 1624350*pow2(lmMsq) - 176400*pow2(lmMsusy)))*pow3(Msusy))
        /(231525.*pow4(Msq)) - (10*(-222264*At*(-17 + 2*lmMsusy*(-7 + 3*lmMt) -
        2*lmMsq*(-7 + 3*lmMsusy + 3*lmMt) + 6*pow2(lmMsq)) + Msusy*(-6262157 +
        222264*lmMt + 252*lmMsusy*(-20233 + 7938*lmMt) - 252*lmMsq*(-19351 +
        6678*lmMsusy + 7938*lmMt) + 1841616*pow2(lmMsq) - 158760*pow2(lmMsusy))
        )*pow5(Msusy))/(750141.*pow6(Msq)) + (xMsusy*(76.53372960293683 - (
        687056*lmMsq)/9801. + (71.76726864605652 + (700*lmMsq)/33.)*lmMsusy + (
        5*(-3 + 44*lmMsq - 44*lmMsusy)*lmMt)/9. - (2260*pow2(lmMsq))/99. + (
        160*pow2(lmMsusy))/99.)*pow8(Msusy) + (2*At*(240*(-9 + lmMsq*(4 - 3*
        lmMsusy - 3*lmMt) + lmMsusy*(-4 + 3*lmMt) + 3*pow2(lmMsq))*pow2(Msusy)*
        pow6(Msq) + 8*(104 + 157*lmMsusy + 84*lmMt + 24*lmMsusy*lmMt - 30*
        lmMsq*(5 + 6*lmMsusy + 3*lmMt) + 135*pow2(lmMsq) + 228*pow2(lmMsusy) +
        108*pow2(lmMt))*pow8(Msq) + 5*(-419 + 36*lmMsusy*(-11 + 4*lmMt) - 36*
        lmMsq*(-11 + 4*lmMsusy + 4*lmMt) + 144*pow2(lmMsq))*pow8(Msusy)))/(81.*
        Msusy) - (8*z2*(100*pow4(Msusy)*pow6(Msq) + 70*pow4(Msq)*pow6(Msusy) +
        3*(30 + 20*lmMsq - 10*(lmMsusy + lmMt) + shiftst3)*pow2(Msusy)*pow8(
        Msq) - 20*At*(2*pow3(Msusy)*(pow2(Msusy)*pow4(Msq) + pow2(Msq)*pow4(
        Msusy) + pow6(Msq) + pow6(Msusy)) + 3*Msusy*pow8(Msq)) + 60*pow2(Msq)*
        pow8(Msusy) + 30*(-2 + shiftst1 + shiftst2)*xMsq*power10(Msq) + 55*
        xMsusy*power10(Msusy)))/(9.*pow2(Msusy)))/pow8(Msq))));
}

/**
 * @return The off-diagonal (1, 2) = (2, 1) matrix element of the Higgs mass matrix as a double for the hierarchy 'H4'
 */
double H4::getS12() const {
   return (Al4p*At*Sbeta*pow4(Mt)*(108*(1 + lmMsusy + lmMt)*Msusy*twoLoopFlag*pow8(
        Msq) - Al4p*threeLoopFlag*(60*(33 - 22*lmMsq + 2*lmMsusy*(11 - 6*lmMt)
        + 12*lmMsq*(lmMsusy + lmMt) - 24*z2 - 12*pow2(lmMsq))*pow4(Msq)*pow5(
        Msusy) + 240*(9 - 4*lmMsq + lmMsusy*(4 - 3*lmMt) + 3*lmMsq*(lmMsusy +
        lmMt) - 6*z2 - 3*pow2(lmMsq))*pow3(Msusy)*pow6(Msq) + 120*(17 - 14*
        lmMsq + 2*lmMsusy*(7 - 3*lmMt) + 6*lmMsq*(lmMsusy + lmMt) - 12*z2 - 6*
        pow2(lmMsq))*pow2(Msq)*pow7(Msusy) - 2*(3*At*xAt*(-349 + 56*lmMsusy -
        24*lmMt + 282*z3 + 32*pow2(lmMsusy)) + 2*Msusy*(208 + 314*lmMsusy + 24*
        (7 + 2*lmMsusy)*lmMt - 60*lmMsq*(5 + 6*lmMsusy + 3*lmMt) + 540*z2 -
        477*z3 + 270*pow2(lmMsq) + 456*pow2(lmMsusy) + 216*pow2(lmMt)))*pow8(
        Msq) + 5*(419 + 36*lmMsusy*(11 - 4*lmMt) + 36*lmMsq*(-11 + 4*lmMsusy +
        4*lmMt) - 288*z2 - 144*pow2(lmMsq))*pow9(Msusy))))/(81.*Cbeta*pow2(
        Msusy)*pow8(Msq));
}

/**
 * @return returns the susy log^0 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double H4::calc_coef_at_as2_no_sm_logs_log0() const {

   const double result =
      (-(pow4(Mt)*(-591666768*(-10589 + 7500*z2)*pow4(Msusy)*pow6(
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
        Msusy)))/(4.992188355e10*pow2(Msusy)*pow8(Msq)))/pow4(Mt)*
        12.;

   return result;
}

/**
 * @return returns the susy log^1 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double H4::calc_coef_at_as2_no_sm_logs_log1() const {

   const double result =
      (-(pow4(Mt)*(2160*pow4(Msusy)*pow6(Msq) + 540*pow4(Msq)*pow6(
        Msusy) - 16*(-173 + 135*z2 + 54*z3)*pow2(Msusy)*pow8(Msq) + 240*pow2(
        Msq)*pow8(Msusy) + 180*log(pow2(Msq)/pow2(Msusy))*pow2(Msusy)*(14*pow4(
        Msq)*pow4(Msusy) + 20*pow2(Msusy)*pow6(Msq) + 12*pow2(Msq)*pow6(Msusy)
        + 28*pow8(Msq) + 11*pow8(Msusy)) + 4320*power10(Msq) + 135*power10(
        Msusy)))/(81.*pow2(Msusy)*pow8(Msq)))/pow4(Mt)*12.;

   return result;
}

/**
 * @return returns the susy log^2 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double H4::calc_coef_at_as2_no_sm_logs_log2() const {

   const double result =
      ((8*(221 + 45*log(pow2(Msq)/pow2(Msusy)))*pow4(Mt))/27.)/
      pow4(Mt)*12.;

   return result;
}

/**
 * @return returns the susy log^3 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double H4::calc_coef_at_as2_no_sm_logs_log3() const {

   const double result =
      ((-224*pow4(Mt))/9.)/pow4(Mt)*12.;

   return result;
}

} // namespace hierarchies
} // namespace himalaya
