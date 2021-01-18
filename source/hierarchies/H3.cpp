// ====================================================================
// This file is part of Himalaya.
//
// Himalaya is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================

#include "H3.hpp"
#include "Constants.hpp"
#include "HimalayaFlags.hpp"
#include "Powers.hpp"
#include <cmath>

namespace himalaya {
namespace hierarchies {

/**
 * Constuctor
 * @param expansionDepth the flagMap for the truncation of expansion variables
 * @param Al4p a double alpha_s/4/Pi
 * @param beta a double which is the mixing angle beta
 * @param Dmglst1 a double Mgl - Mst1
 * @param Dmst12 a double Mst1^2 - Mst2^2
 * @param Dmsqst1 a double Msq^2 - Mst1^2
 * @param lmMt a double log((renormalization scale / Mt)^2)
 * @param lmMst1 a double log((renormalization scale / Mst1)^2)
 * @param Mgl a double gluino mass
 * @param Mt a double top/bottom quark mass
 * @param Mst1 a double stop 1 mass
 * @param Mst2 a double stop 2 mass
 * @param MuSUSY a double mu parameter
 * @param s2t a double 2 times the sine of the stop/sbottom quark mixing angle
 * @param mdrFlag an int 0 for DR and 1 for MDR scheme
 * @param oneLoopFlag an int flag to consider the one-loop expansion terms
 * @param twoLoopFlag an int flag to consider the two-loop expansion terms
 * @param threeLoopFlag an int flag to consider the three-loop expansion terms
 */
H3::H3 (const ExpansionFlags_t& expansionDepth, double Al4p, double beta,
                   double Dmglst1, double Dmst12, double Dmsqst1, double lmMt, double lmMst1,
                   double Mgl, double Mt, double Mst1, double Mst2, double Msq, double MuSUSY,
                   double s2t,
                   int mdrFlag, int oneLoopFlag, int twoLoopFlag, int threeLoopFlag) {
   // abbrev for tan(beta) and sin(beta)
   Tbeta = tan(beta);
   Sbeta = sin(beta);
   this -> Dmglst1 = Dmglst1;
   this -> Dmst12 = Dmst12;
   this -> Dmsqst1 = Dmsqst1;
   this -> lmMst1 = lmMst1;
   this -> Mgl = Mgl;
   this -> Mt = Mt;
   this -> Mst1 = Mst1;
   this -> Mst2 = Mst2;
   this -> Msq = Msq;
   this -> MuSUSY = MuSUSY;
   this -> s2t = s2t;
   this -> lmMt = lmMt;
   this -> oneLoopFlag = oneLoopFlag;
   this -> twoLoopFlag = twoLoopFlag;
   this -> threeLoopFlag = threeLoopFlag;
   this -> Al4p = Al4p;
   // mdr flags, indicates if one wants to shift the dr stop mass to the mdr stop mass
   shiftst1 = mdrFlag;
   shiftst2 = mdrFlag;
   shiftst3 = mdrFlag;
   // expansion flags
   xDmglst1 = expansionDepth.at(ExpansionDepth::Dmglst1);
   xDmst12 = expansionDepth.at(ExpansionDepth::Dmglst1);
   xDmsqst1 = expansionDepth.at(ExpansionDepth::Dmsqst1);
}

/**
 * @return The diagonal (1, 1) matrix element of the Higgs mass matrix as a double for the hierarchy 'H3'
 */
double H3::getS1() const {
   return -(pow2(Mt)*pow2(MuSUSY)*(-4*xDmst12*pow3(Dmst12)*(4050*oneLoopFlag*pow2(
        s2t) + (8*Al4p*(450*Mgl*Mst1*s2t*twoLoopFlag*(4*Dmglst1*Mgl*((5 - 6*
        lmMst1)*Mt + 4*(1 - 3*lmMst1)*Mst1*s2t) - 4*(Mt + 6*lmMst1*Mt + 3*(-5 +
        6*lmMst1)*Mst1*s2t)*pow2(Dmglst1) + (-4*(5 + 6*lmMst1)*Mt + (1 + 6*
        lmMst1)*Mst1*s2t)*pow2(Mgl)) + (xDmglst1*(-360*Mst1*s2t*(3*(9 + 10*
        lmMst1)*Mt + 10*(-13 + 12*lmMst1)*Mst1*s2t)*twoLoopFlag*pow2(Msq) +
        Al4p*threeLoopFlag*(30000*Dmsqst1*Mst1*Mt*s2t + pow2(Msq)*(2*Mst1*Mt*
        s2t*(37824007 + 770520*lmMst1 - 131400*pow2(lmMst1)) + (59957863 +
        480000*lmMst1 - 26880*lmMt - 230400*pow2(lmMst1))*pow2(Mt) + 15*(-
        3044017 - 27472*lmMst1 + 48480*pow2(lmMst1))*pow2(Mst1)*pow2(s2t))))*
        pow3(Dmglst1))/pow2(Msq)))/(pow2(Mst1)*pow3(Mgl)) + threeLoopFlag*pow2(
        Al4p)*((Mt*(Mt*(72*pow2(Dmglst1)*(3891491 + 27200*lmMst1 - 960*lmMt -
        19200*pow2(lmMst1)) + 200*Dmglst1*Mgl*(403559 + 384*(lmMst1 + lmMt) -
        4608*pow2(lmMst1)) + 15*(-1763661 + 47104*lmMst1 - 5120*lmMt + 24576*
        pow2(lmMst1))*pow2(Mgl)) + (240*Mst1*s2t*(5*pow2(Mgl)*(840*Dmsqst1 + (-
        36863 + 80*lmMst1 + 552*pow2(lmMst1))*pow2(Msq)) + 10*pow2(Dmglst1)*(
        100*Dmsqst1 + (-32829 + 1852*lmMst1 + 660*pow2(lmMst1))*pow2(Msq)) +
        Dmglst1*Mgl*(1000*Dmsqst1 + (-1282471 + 7264*lmMst1 + 18120*pow2(
        lmMst1))*pow2(Msq))))/pow2(Msq)))/(pow2(Mgl)*pow2(Mst1)) + 15*(350605 +
        4320*shiftst1 + 2880*shiftst2 + 8352*shiftst3 - 96*lmMst1*(-115 + 90*
        shiftst1 + 60*shiftst2 + 54*shiftst3) - 2160*pow2(lmMst1) + (40*
        Dmglst1*(-84209 - 1264*lmMst1 + 240*pow2(lmMst1)))/Mgl + (8*pow2(
        Dmglst1)*(-1732531 - 16896*lmMst1 + 24840*pow2(lmMst1)))/pow2(Mgl) + (
        2400*Dmsqst1*(7 - 24*lmMst1*(-1 + shiftst2) + 36*shiftst2))/pow2(Msq))*
        pow2(s2t) + (12000*s2t*xDmsqst1*pow2(Dmsqst1)*(20*Dmglst1*Mt*(Dmglst1*
        Mgl + xDmglst1*pow2(Dmglst1) + pow2(Mgl)) + 3*(28*Mt + Mst1*s2t*(7 + 6*
        shiftst1 - 24*lmMst1*(-1 + shiftst2) + 30*shiftst2))*pow3(Mgl)))/(Mst1*
        pow3(Mgl)*pow4(Msq)))) + 8*pow2(Mst2)*((-1800*Al4p*Dmst12*s2t*
        twoLoopFlag*(pow2(Mgl)*(Dmst12*(4*Mt - 3*Mst1*s2t) + 12*(1 + 2*lmMst1)*
        Mt*pow2(Mst2)) + pow2(Dmglst1)*(6*Dmst12*(4*Mt + (-5 + 6*lmMst1)*Mst1*
        s2t) + 4*(-11 + 6*lmMst1)*Mt*pow2(Mst2)) + 4*Dmglst1*Mgl*(2*Dmst12*(-1
        + 3*lmMst1)*Mst1*s2t + (-5 + 6*lmMst1)*Mt*pow2(Mst2))))/(Mst1*pow2(Mgl)
        ) + 2025*oneLoopFlag*pow2(Dmst12)*pow2(s2t) - (2*Al4p*xDmglst1*pow3(
        Dmglst1)*(720*Dmst12*Mst1*s2t*twoLoopFlag*pow2(Msq)*(52*Dmst12*Mt + 5*
        Dmst12*(-13 + 12*lmMst1)*Mst1*s2t + (-77 + 30*lmMst1)*Mt*pow2(Mst2)) +
        Al4p*threeLoopFlag*(2*Dmst12*Mt*(30000*Dmsqst1*Mst1*s2t + (4*Mst1*s2t*(
        31025111 + 290880*lmMst1 - 251100*pow2(lmMst1)) + Mt*(59957863 +
        480000*lmMst1 - 26880*lmMt - 230400*pow2(lmMst1)))*pow2(Msq))*pow2(
        Mst2) + pow2(Dmst12)*(-60000*Dmsqst1*Mst1*Mt*s2t + pow2(Msq)*(2*Mst1*
        Mt*s2t*(-99874229 - 1352280*lmMst1 + 633600*pow2(lmMst1)) + 2*(-
        59957863 - 480000*lmMst1 + 26880*lmMt + 230400*pow2(lmMst1))*pow2(Mt) +
        15*(3044017 + 27472*lmMst1 - 48480*pow2(lmMst1))*pow2(Mst1)*pow2(s2t)))
        - 24*(-3877891 - 46400*lmMst1 + 960*lmMt + 19200*pow2(lmMst1))*pow2(
        Msq)*pow2(Mt)*pow4(Mst2))))/(pow2(Msq)*pow2(Mst1)*pow3(Mgl))) + (
        threeLoopFlag*pow2(Al4p)*(-480*Mgl*pow2(Mst1)*pow2(s2t)*((Dmst12*pow2(
        Msq)*pow2(Mst2)*(300*Dmsqst1*pow2(Mgl)*(Dmst12*(-7 - 18*shiftst1 - 18*
        shiftst2 + 12*lmMst1*(-2 + shiftst1 + shiftst2)) - 12*(-3 + 2*lmMst1)*(
        shiftst1 - shiftst2)*pow2(Mst2)) + pow2(Msq)*(Dmst12*pow2(Dmglst1)*(
        1732531 + 16896*lmMst1 - 24840*pow2(lmMst1)) + 5*Dmglst1*Dmst12*Mgl*(
        84209 + 1264*lmMst1 - 240*pow2(lmMst1)) + 10*pow2(Mgl)*(Dmst12*(1429 -
        180*shiftst1 - 180*shiftst2 + lmMst1*(-454 + 360*shiftst1 + 360*
        shiftst2) + 24*pow2(lmMst1)) - 360*(-1 + 2*lmMst1)*(shiftst1 -
        shiftst2)*pow2(Mst2)))))/2. + pow2(Mgl)*(90*Dmst12*shiftst3*pow2(Mst2)*
        (Dmst12*(-7 + 6*lmMst1) + 2*(1 - 2*lmMst1)*pow2(Mst2))*pow4(Msq) + 36*
        z2*(-2*xDmst12*pow3(Dmst12)*(100*Dmsqst1*shiftst2*(Dmsqst1*xDmsqst1 +
        pow2(Msq)) + (15*shiftst1 + 10*shiftst2 + 9*shiftst3)*pow4(Msq)) +
        Dmst12*pow2(Mst2)*(5*shiftst3*(3*Dmst12 - 2*pow2(Mst2))*pow4(Msq) - 50*
        (-(Dmst12*(shiftst1 + shiftst2)) + 2*(shiftst1 - shiftst2)*pow2(Mst2))*
        (xDmsqst1*pow2(Dmsqst1) + Dmsqst1*pow2(Msq) + pow4(Msq)))))) + pow2(
        Mst2)*(24000*Dmst12*Mst1*s2t*xDmsqst1*pow2(Dmsqst1)*(40*Dmglst1*Mt*(
        Dmglst1*Mgl + xDmglst1*pow2(Dmglst1) + pow2(Mgl))*(Dmst12 - pow2(Mst2))
        - 3*(-(Dmst12*(56*Mt + Mst1*s2t*(7 + 30*shiftst1 + 6*shiftst2 - 12*
        lmMst1*(-2 + shiftst1 + shiftst2)))) + 8*(7*Mt - 3*(-2 + lmMst1)*Mst1*
        s2t*(shiftst1 - shiftst2))*pow2(Mst2))*pow3(Mgl)) + 2*Mgl*Mt*pow2(Msq)*
        (240*Dmst12*Mst1*s2t*(5*pow2(Mgl)*(1680*Dmsqst1*(Dmst12 - pow2(Mst2)) +
        pow2(Msq)*(3*Dmst12*(-10473 + 40*lmMst1 + 256*pow2(lmMst1)) - 8*(1361 +
        10*lmMst1 + 54*pow2(lmMst1))*pow2(Mst2))) + Dmglst1*Mgl*(2000*Dmsqst1*(
        Dmst12 - pow2(Mst2)) + pow2(Msq)*(Dmst12*(-949861 + 1944*lmMst1 +
        11520*pow2(lmMst1)) + 20*(-33261 + 532*lmMst1 + 660*pow2(lmMst1))*pow2(
        Mst2))) + pow2(Dmglst1)*(2000*Dmsqst1*(Dmst12 - pow2(Mst2)) + pow2(Msq)
        *(Dmst12*(958501 + 24456*lmMst1 - 11520*pow2(lmMst1)) + 2*(-1286791 -
        5936*lmMst1 + 18120*pow2(lmMst1))*pow2(Mst2)))) + Mt*pow2(Msq)*(16*
        pow2(Dmglst1)*(-9*Dmst12*(-3891491 - 27200*lmMst1 + 960*lmMt + 19200*
        pow2(lmMst1))*(Dmst12 - pow2(Mst2)) - 50*(345581 + 4896*lmMst1 + 96*
        lmMt - 3456*pow2(lmMst1))*pow4(Mst2)) + 400*Dmglst1*Mgl*(Dmst12*(403559
        + 384*(lmMst1 + lmMt) - 4608*pow2(lmMst1))*(Dmst12 - pow2(Mst2)) - 24*(
        9631 + 16*lmMst1 + 48*lmMt - 192*pow2(lmMst1))*pow4(Mst2)) - 15*pow2(
        Mgl)*(pow2(Dmst12)*(852541 + 9216*lmMst1 - 10240*lmMt + 6144*pow2(
        lmMst1)) + 160*Dmst12*(11389 - 704*lmMst1 + 192*lmMt - 384*pow2(lmMst1)
        )*pow2(Mst2) + 1920*(349 - 56*lmMst1 + 24*lmMt - 32*pow2(lmMst1))*pow4(
        Mst2))))) - 225*z3*(-2*xDmst12*pow3(Dmst12)*((-16*Mgl*Mst1*Mt*s2t*(
        142987*Dmglst1*Mgl + 37582*pow2(Dmglst1) + 20297*pow2(Mgl)) + 3*Mgl*(
        197112*Dmglst1*Mgl + 687960*pow2(Dmglst1) - 65963*pow2(Mgl))*pow2(Mt) +
        8*xDmglst1*(557078*Mst1*Mt*s2t + 442053*pow2(Mt) - 349745*pow2(Mst1)*
        pow2(s2t))*pow3(Dmglst1))*pow4(Msq) - Mgl*pow2(Mst1)*pow2(s2t)*(-10080*
        Dmsqst1*pow2(Mgl)*(Dmsqst1*xDmsqst1 + pow2(Msq)) + (403880*Dmglst1*Mgl
        + 1600920*pow2(Dmglst1) - 37669*pow2(Mgl))*pow4(Msq))) + pow2(Mst2)*(-
        8*(Mgl*pow2(Dmst12)*pow2(Mst1)*pow2(s2t)*(-1260*Dmsqst1*pow2(Mgl)*(
        Dmsqst1*xDmsqst1 + pow2(Msq)) + (50485*Dmglst1*Mgl + 200115*pow2(
        Dmglst1) + 2574*pow2(Mgl))*pow4(Msq)) + xDmglst1*pow3(Dmglst1)*pow4(
        Msq)*(18*Dmst12*Mt*(49117*Mt + 102024*Mst1*s2t)*pow2(Mst2) + pow2(
        Dmst12)*(-1475294*Mst1*Mt*s2t - 884106*pow2(Mt) + 349745*pow2(Mst1)*
        pow2(s2t)) + 687960*pow2(Mt)*pow4(Mst2))) + Mgl*Mt*pow4(Msq)*(16*
        Dmst12*Mst1*s2t*(5*Dmst12*(-21081*Dmglst1*Mgl + 21081*pow2(Dmglst1) -
        3457*pow2(Mgl)) - 2*(37582*Dmglst1*Mgl + 142987*pow2(Dmglst1) + 3012*
        pow2(Mgl))*pow2(Mst2)) + 3*Mt*(48*Dmglst1*Mgl*(8213*pow2(Dmst12) -
        8213*Dmst12*pow2(Mst2) - 4664*pow4(Mst2)) + 432*pow2(Dmglst1)*(3185*
        pow2(Dmst12) - 3185*Dmst12*pow2(Mst2) - 1566*pow4(Mst2)) - pow2(Mgl)*(
        31963*pow2(Dmst12) + 68000*Dmst12*pow2(Mst2) + 24064*pow4(Mst2))))))))/
        (pow2(Mst1)*pow3(Mgl)*pow4(Msq))))/(777600.*pow6(Mst2));
}

/**
 * @return The diagonal (2, 2) matrix element of the Higgs mass matrix as a double for the hierarchy 'H3'
 */
double H3::getS2() const {
   return (pow2(Mst2)*(2*Al4p*Mgl*Mt*pow2(Msq)*(-70560*Mst1*Mt*twoLoopFlag*pow2(
        Msq)*(2*MuSUSY*pow2(Sbeta)*(-3*pow2(Dmglst1)*(2*Dmst12*Mt*((477 - 330*
        lmMst1 + 30*lmMt)*Mt + 50*(13 - 18*lmMst1)*Mst1*s2t)*pow2(Mst2) - pow2(
        Dmst12)*(70*Mst1*Mt*s2t + 8*(92 - 30*lmMst1 + 5*lmMt)*pow2(Mt) + 25*(11
        - 6*lmMst1)*pow2(Mst1)*pow2(s2t)) + 400*(4 - 3*lmMst1)*pow2(Mt)*pow4(
        Mst2)) + 25*pow2(Mgl)*(4*Dmst12*Mt*(2*(5 + 6*lmMst1 + 3*lmMt)*Mt - 9*(1
        + 2*lmMst1)*Mst1*s2t)*pow2(Mst2) - pow2(Dmst12)*(24*(1 - 3*lmMst1)*
        Mst1*Mt*s2t + 2*(5 + 6*(lmMst1 + lmMt))*pow2(Mt) + 9*(1 + 2*lmMst1)*
        pow2(Mst1)*pow2(s2t)) + 72*(1 + lmMst1 + lmMt)*pow2(Mt)*pow4(Mst2)) -
        Dmglst1*Mgl*(-300*Dmst12*Mt*((-3 + 6*lmMst1)*Mt + 2*(-1 + 6*lmMst1)*
        Mst1*s2t)*pow2(Mst2) + 3*pow2(Dmst12)*(-100*Mst1*Mt*s2t + 2*(-41 + 90*
        lmMst1 + 10*lmMt)*pow2(Mt) + 25*(-5 + 6*lmMst1)*pow2(Mst1)*pow2(s2t)) -
        200*(-7 + 15*lmMst1 + 3*lmMt)*pow2(Mt)*pow4(Mst2))) + Tbeta*(-25*
        Dmst12*s2t*(4*Dmglst1*Mgl*((-5 + 6*lmMst1)*Mt*pow2(Mst2)*pow2(MuSUSY) -
        2*Dmst12*Mst1*s2t*((-1 + 3*lmMst1)*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 3*
        (-1 + 6*lmMst1)*pow2(Mst1)*pow2(Sbeta))) + 2*pow2(Dmglst1)*(2*(-11 + 6*
        lmMst1)*Mt*pow2(Mst2)*pow2(MuSUSY) + 3*Dmst12*(4*Mt*pow2(MuSUSY) - (-5
        + 6*lmMst1)*Mst1*s2t*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 2*(13 - 18*
        lmMst1)*s2t*pow2(Sbeta)*pow3(Mst1))) + pow2(Mgl)*(12*(1 + 2*lmMst1)*Mt*
        pow2(Mst2)*pow2(MuSUSY) + Dmst12*(4*Mt*pow2(MuSUSY) + 3*Mst1*s2t*pow2(
        MuSUSY)*(-1 + pow2(Sbeta)) + 36*(1 + 2*lmMst1)*s2t*pow2(Sbeta)*pow3(
        Mst1)))) + pow2(Sbeta)*(4*Dmst12*Mt*s2t*(-25*pow2(Mgl)*(Dmst12*(4*(1 +
        3*lmMst1 + 6*lmMt)*pow2(Mst1) - pow2(MuSUSY)) - 3*pow2(Mst2)*(6*(1 + 2*
        (lmMst1 + lmMt))*pow2(Mst1) + (1 + 2*lmMst1)*pow2(MuSUSY))) + pow2(
        Dmglst1)*(-6*Dmst12*((-169 + 135*lmMst1 + 15*lmMt)*pow2(Mst1) - 25*
        pow2(MuSUSY)) + 25*pow2(Mst2)*(24*(-4 + 3*lmMst1)*pow2(Mst1) + (-11 +
        6*lmMst1)*pow2(MuSUSY))) + 25*Dmglst1*Mgl*(-4*Dmst12*(-4 + 6*lmMst1 +
        3*lmMt)*pow2(Mst1) + pow2(Mst2)*(2*(-17 + 30*lmMst1 + 6*lmMt)*pow2(
        Mst1) + (-5 + 6*lmMst1)*pow2(MuSUSY)))) - 2*Mst1*pow2(Mt)*(4*Dmglst1*
        Mgl*((11 + 60*lmMst1 - 60*lmMt)*pow2(Dmst12) + 25*Dmst12*(1 + 30*lmMst1
        + 6*lmMt)*pow2(Mst2) + 100*(1 + 12*lmMst1 + 6*lmMt)*pow4(Mst2)) + 2*
        pow2(Dmglst1)*((-83 + 120*lmMst1 - 120*lmMt)*pow2(Dmst12) + Dmst12*(-
        1381 + 2490*lmMst1 + 210*lmMt)*pow2(Mst2) + 200*(-11 + 21*lmMst1 + 6*
        lmMt)*pow4(Mst2)) + 25*pow2(Mgl)*((5 + 42*lmMst1 + 30*lmMt)*pow2(
        Dmst12) - 8*Dmst12*(4 + 3*lmMst1 + 6*lmMt)*pow2(Mst2) - 72*(1 + lmMst1
        - lmMt + 2*lmMst1*lmMt + pow2(lmMst1) - 3*pow2(lmMt))*pow4(Mst2)))))) +
        Al4p*threeLoopFlag*(784*Mst1*MuSUSY*pow2(Sbeta)*(225*Mt*pow2(Dmst12)*(
        pow2(Dmglst1)*(1000*Dmsqst1 + (1286791 + 5936*lmMst1 - 18120*pow2(
        lmMst1))*pow2(Msq)) + 10*Dmglst1*Mgl*(100*Dmsqst1 + (33261 - 532*lmMst1
        - 660*pow2(lmMst1))*pow2(Msq)) + 20*pow2(Mgl)*(210*Dmsqst1 + (1361 +
        10*lmMst1 + 54*pow2(lmMst1))*pow2(Msq)))*pow2(Mst1)*pow2(s2t) + pow3(
        Mt)*(2*Dmglst1*(1200*Dmglst1*Dmsqst1*(Dmst12*(557 + 120*lmMst1 - 120*
        lmMt)*(Dmst12 - pow2(Mst2)) - 16*(4 + 15*lmMst1 - 15*lmMt)*pow4(Mst2))
        + 1200*Dmsqst1*Mgl*(Dmst12*(557 + 120*lmMst1 - 120*lmMt)*(Dmst12 -
        pow2(Mst2)) - 50*(26 + 3*lmMst1 - 3*lmMt)*pow4(Mst2)) + Dmglst1*pow2(
        Msq)*(-3*pow2(Dmst12)*(-129193181 - 401100*lmMt + 100*lmMst1*(-7351 +
        1164*lmMt) + 1336800*pow2(lmMst1) + 28800*pow2(lmMt)) + 2*Dmst12*(-
        327941741 - 686700*lmMt + 1080*lmMst1*(-44 + 445*lmMt) + 3531600*pow2(
        lmMst1) + 64800*pow2(lmMt))*pow2(Mst2) + 80*(-3777727 - 24255*lmMt +
        486*lmMst1*(118 + 45*lmMt) + 52380*pow2(lmMst1))*pow4(Mst2)) - Mgl*
        pow2(Msq)*(-3*pow2(Dmst12)*(-39474953 + lmMst1*(52780 - 87600*lmMt) -
        3300*lmMt + 319200*pow2(lmMst1) + 14400*pow2(lmMt)) + 40*Dmst12*(
        3746977 + 4005*lmMt - 18*lmMst1*(2711 + 1215*lmMt) - 52380*pow2(lmMst1)
        )*pow2(Mst2) + 6000*(9961 + 66*lmMt - 10*lmMst1*(67 + 24*lmMt) + 42*
        pow2(lmMst1) + 72*pow2(lmMt))*pow4(Mst2))) + 15*pow2(Mgl)*(4000*
        Dmsqst1*(Dmst12*(65 - 6*lmMst1 + 6*lmMt)*(Dmst12 - pow2(Mst2)) + 6*(-20
        + 3*lmMst1 - 3*lmMt)*pow4(Mst2)) - pow2(Msq)*(-(pow2(Dmst12)*(-2284899
        + 49840*lmMt - 32*lmMst1*(-1793 + 555*lmMt) + 87360*pow2(lmMst1) +
        28800*pow2(lmMt))) + 200*Dmst12*(11697 + 448*lmMst1 + 330*lmMt - 372*
        lmMst1*lmMt + 408*pow2(lmMst1) + 288*pow2(lmMt))*pow2(Mst2) + 1600*(434
        - 83*lmMst1 + 174*lmMt - 66*lmMst1*lmMt + 183*pow2(lmMst1) + 108*pow2(
        lmMt))*pow4(Mst2)))) + 450*Mst1*s2t*(Dmst12*(pow2(Dmglst1)*(19160*
        Dmsqst1*(Dmst12 - pow2(Mst2)) + pow2(Msq)*(Dmst12*(-586073 - 9268*
        lmMst1 - 448*lmMt + 19200*pow2(lmMst1)) + 2*(13289 - 916*lmMst1 + 80*
        lmMt - 60*pow2(lmMst1))*pow2(Mst2))) + 5*pow2(Mgl)*(5*Dmsqst1*(161*
        Dmst12 + 24*(1 + 12*lmMst1)*pow2(Mst2)) + pow2(Msq)*(Dmst12*(2071 +
        296*lmMst1 - 96*lmMt - 600*pow2(lmMst1)) + 8*(631 - lmMst1 + 36*lmMt +
        21*pow2(lmMst1))*pow2(Mst2))) + 4*Dmglst1*Mgl*(4700*Dmsqst1*(Dmst12 -
        pow2(Mst2)) - pow2(Msq)*(Dmst12*(17539 + 574*lmMst1 + 160*lmMt - 1920*
        pow2(lmMst1)) + 5*(-4539 + 596*lmMst1 - 48*lmMt + 516*pow2(lmMst1))*
        pow2(Mst2))))*pow2(Mt) + 45*pow2(Mgl)*(shiftst3*pow2(Msq)*(-8*Dmst12*(-
        3 + 2*lmMst1)*pow2(Mst2)*pow2(Mt) + pow2(Dmst12)*(16*(-2 + lmMst1)*
        pow2(Mt) + (1 - 2*lmMst1)*pow2(Mst1)*pow2(s2t)) + 8*(-1 + 2*lmMst1)*
        pow2(Mt)*pow4(Mst2)) + 10*(Dmsqst1*(-3 + 2*lmMst1) + (-1 + 2*lmMst1)*
        pow2(Msq))*(-8*Dmst12*shiftst2*pow2(Mst2)*pow2(Mt) + (-shiftst1 +
        shiftst2)*pow2(Dmst12)*pow2(Mst1)*pow2(s2t) + 8*(shiftst1 - shiftst2)*
        pow2(Mt)*pow4(Mst2))))) + Mt*Tbeta*(245*(pow2(MuSUSY)*(8*pow2(Dmglst1)*
        (6*Dmst12*Mt*(10000*Dmsqst1*Mst1*s2t + (3*Mt*(3891491 + 27200*lmMst1 -
        960*lmMt - 19200*pow2(lmMst1)) + 10*Mst1*s2t*(1286791 + 5936*lmMst1 -
        18120*pow2(lmMst1)))*pow2(Msq))*pow2(Mst2) - 3*pow2(Dmst12)*(20000*
        Dmsqst1*Mst1*Mt*s2t + pow2(Msq)*(10*Mst1*Mt*s2t*(958501 + 24456*lmMst1
        - 11520*pow2(lmMst1)) + 6*(3891491 + 27200*lmMst1 - 960*lmMt - 19200*
        pow2(lmMst1))*pow2(Mt) + 5*(-1732531 - 16896*lmMst1 + 24840*pow2(
        lmMst1))*pow2(Mst1)*pow2(s2t))) + 100*(345581 + 4896*lmMst1 + 96*lmMt -
        3456*pow2(lmMst1))*pow2(Msq)*pow2(Mt)*pow4(Mst2)) + 40*Dmglst1*Mgl*(10*
        Dmst12*Mt*(1200*Dmsqst1*Mst1*s2t + (Mt*(403559 + 384*(lmMst1 + lmMt) -
        4608*pow2(lmMst1)) - 12*Mst1*s2t*(-33261 + 532*lmMst1 + 660*pow2(
        lmMst1)))*pow2(Msq))*pow2(Mst2) - pow2(Dmst12)*(12000*Dmsqst1*Mst1*Mt*
        s2t + pow2(Msq)*(6*Mst1*Mt*s2t*(-949861 + 1944*lmMst1 + 11520*pow2(
        lmMst1)) + 10*(403559 + 384*(lmMst1 + lmMt) - 4608*pow2(lmMst1))*pow2(
        Mt) + 15*(-84209 - 1264*lmMst1 + 240*pow2(lmMst1))*pow2(Mst1)*pow2(s2t)
        )) + 240*(9631 + 16*lmMst1 + 48*lmMt - 192*pow2(lmMst1))*pow2(Msq)*
        pow2(Mt)*pow4(Mst2)) + 15*pow2(Mgl)*(160*Dmst12*Mt*(840*Dmsqst1*Mst1*
        s2t + (Mt*(11389 - 704*lmMst1 + 192*lmMt - 384*pow2(lmMst1)) + 4*Mst1*
        s2t*(1361 + 10*lmMst1 + 54*pow2(lmMst1)))*pow2(Msq))*pow2(Mst2) + pow2(
        Dmst12)*(-2400*Dmsqst1*Mst1*s2t*(56*Mt + (7 + 24*lmMst1)*Mst1*s2t) +
        pow2(Msq)*(-240*Mst1*Mt*s2t*(-10473 + 40*lmMst1 + 256*pow2(lmMst1)) + (
        852541 + 9216*lmMst1 - 10240*lmMt + 6144*pow2(lmMst1))*pow2(Mt) + 80*(
        1429 - 454*lmMst1 + 24*pow2(lmMst1))*pow2(Mst1)*pow2(s2t))) + 1920*(349
        - 56*lmMst1 + 24*lmMt - 32*pow2(lmMst1))*pow2(Msq)*pow2(Mt)*pow4(Mst2))
        ) + 21600*pow2(Mgl)*pow2(Mst1)*(-(shiftst3*pow2(Msq)*(-2*Dmst12*pow2(
        Mst2)*((-1 + 2*lmMst1)*pow2(MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta)) + 12*(
        (3 - 2*lmMst1)*pow2(Mt) + (-1 + 2*lmMst1)*pow2(Mst1)*pow2(s2t))*pow2(
        Sbeta)) + pow2(Dmst12)*((-7 + 6*lmMst1)*pow2(MuSUSY)*pow2(s2t)*(-1 +
        pow2(Sbeta)) + 6*(-8*(-2 + lmMst1)*pow2(Mt) + (-15 + 14*lmMst1)*pow2(
        Mst1)*pow2(s2t))*pow2(Sbeta)) + 24*(1 - 2*lmMst1)*pow2(Mt)*pow2(Sbeta)*
        pow4(Mst2))) - 10*(Dmsqst1*(3 - 2*lmMst1) + (1 - 2*lmMst1)*pow2(Msq))*(
        pow2(Dmst12)*pow2(s2t)*(-((shiftst1 + shiftst2)*pow2(MuSUSY)*(-1 +
        pow2(Sbeta))) + 6*(-3*shiftst1 + shiftst2)*pow2(Mst1)*pow2(Sbeta)) + 2*
        Dmst12*pow2(Mst2)*((shiftst1 - shiftst2)*pow2(MuSUSY)*pow2(s2t)*(-1 +
        pow2(Sbeta)) + 12*(shiftst2*pow2(Mt) + (shiftst1 - shiftst2)*pow2(Mst1)
        *pow2(s2t))*pow2(Sbeta)) + 24*(shiftst1 + shiftst2)*pow2(Mt)*pow2(
        Sbeta)*pow4(Mst2)))) + pow2(Sbeta)*(784*Dmst12*Mst1*Mt*s2t*(375*pow2(
        Mgl)*(-80*Dmsqst1*(Dmst12*((-98 + 24*lmMst1 - 24*lmMt)*pow2(Mst1) - 21*
        pow2(MuSUSY)) + 3*pow2(Mst2)*((74 - 12*lmMst1 + 12*lmMt)*pow2(Mst1) +
        7*pow2(MuSUSY))) + pow2(Msq)*(-8*pow2(Mst2)*(8*(347 - 50*lmMst1 + 66*
        lmMt - 66*lmMst1*lmMt + 183*pow2(lmMst1) + 108*pow2(lmMt))*pow2(Mst1) +
        (1361 + 10*lmMst1 + 54*pow2(lmMst1))*pow2(MuSUSY)) + Dmst12*(16*(-4378
        - 517*lmMst1 + 243*lmMt - 78*lmMst1*lmMt + 528*pow2(lmMst1) + 288*pow2(
        lmMt))*pow2(Mst1) + 3*(-10473 + 40*lmMst1 + 256*pow2(lmMst1))*pow2(
        MuSUSY)))) - pow2(Dmglst1)*(1200*Dmsqst1*(Dmst12*((866 - 240*lmMst1 +
        240*lmMt)*pow2(Mst1) - 125*pow2(MuSUSY)) + pow2(Mst2)*(16*(23 + 30*
        lmMst1 - 30*lmMt)*pow2(Mst1) + 125*pow2(MuSUSY))) + pow2(Msq)*(10*pow2(
        Mst2)*(8*(7531199 + 48510*lmMt - 486*lmMst1*(191 + 90*lmMt) - 104760*
        pow2(lmMst1))*pow2(Mst1) + 15*(1286791 + 5936*lmMst1 - 18120*pow2(
        lmMst1))*pow2(MuSUSY)) + Dmst12*(-4*(-176974411 + 218700*lmMt - 540*
        lmMst1*(3971 + 730*lmMt) + 1436400*pow2(lmMst1) + 64800*pow2(lmMt))*
        pow2(Mst1) + 75*(-958501 - 24456*lmMst1 + 11520*pow2(lmMst1))*pow2(
        MuSUSY)))) - 5*Dmglst1*Mgl*(240*Dmsqst1*(25*pow2(Mst2)*(2*(55 + 6*
        lmMst1 - 6*lmMt)*pow2(Mst1) + 5*pow2(MuSUSY)) - Dmst12*(4*(379 + 15*
        lmMst1 - 15*lmMt)*pow2(Mst1) + 125*pow2(MuSUSY))) + pow2(Msq)*(300*
        pow2(Mst2)*(16*(4964 - 3*lmMt - 5*lmMst1*(55 + 24*lmMt) + 21*pow2(
        lmMst1) + 36*pow2(lmMt))*pow2(Mst1) + (33261 - 532*lmMst1 - 660*pow2(
        lmMst1))*pow2(MuSUSY)) + Dmst12*(8*(4511549 + 9810*lmMt + 6*lmMst1*(
        14879 + 4710*lmMt) - 117360*pow2(lmMst1) - 21600*pow2(lmMt))*pow2(Mst1)
        - 15*(-949861 + 1944*lmMst1 + 11520*pow2(lmMst1))*pow2(MuSUSY))))) -
        29400*pow2(Dmst12)*pow2(Mst1)*(pow2(Dmglst1)*(114960*Dmsqst1*pow2(Mst1)
        + pow2(Msq)*(12*(-13249 + 916*lmMst1 - 80*lmMt + 60*pow2(lmMst1))*pow2(
        Mst1) + (1732531 + 16896*lmMst1 - 24840*pow2(lmMst1))*pow2(MuSUSY))) +
        5*Dmglst1*Mgl*(22560*Dmsqst1*pow2(Mst1) + pow2(Msq)*(24*(-4515 + 596*
        lmMst1 - 48*lmMt + 516*pow2(lmMst1))*pow2(Mst1) + (84209 + 1264*lmMst1
        - 240*pow2(lmMst1))*pow2(MuSUSY))) - 10*pow2(Mgl)*(30*Dmsqst1*(12*(1 +
        12*lmMst1)*pow2(Mst1) + (7 + 24*lmMst1)*pow2(MuSUSY)) + pow2(Msq)*(24*(
        613 - lmMst1 + 36*lmMt + 21*pow2(lmMst1))*pow2(Mst1) + (-1429 + 454*
        lmMst1 - 24*pow2(lmMst1))*pow2(MuSUSY))))*pow2(s2t) - pow2(Mt)*(16*
        pow2(Mst1)*(pow2(Dmst12)*(980*pow2(Mgl)*(15*Dmsqst1*(4561 + 510*lmMst1
        - 510*lmMt) + (93973 - 61305*lmMt + 54*lmMst1*(2337 + 1255*lmMt) -
        48420*pow2(lmMst1) - 58050*pow2(lmMt))*pow2(Msq)) + 98*Dmglst1*Mgl*(
        10800*Dmsqst1*(423 - 20*lmMst1 + 20*lmMt) + (-22574599 - 21060*lmMt +
        60*lmMst1*(6677 + 8880*lmMt) + 1598400*pow2(lmMst1) + 172800*pow2(lmMt)
        )*pow2(Msq)) + 3*pow2(Dmglst1)*(11760*Dmsqst1*(13147 - 90*lmMst1 + 90*
        lmMt) + (-5259064991 - 32822300*lmMt + 60*lmMst1*(-104761 + 256200*
        lmMt) + 148327200*pow2(lmMst1) + 5644800*pow2(lmMt))*pow2(Msq))) - 196*
        Dmst12*(10*Dmglst1*Mgl*(540*Dmsqst1*(423 - 20*lmMst1 + 20*lmMt) + (-
        5136871 - 29790*lmMt + 12*lmMst1*(8942 + 4305*lmMt) + 86040*pow2(
        lmMst1) + 21600*pow2(lmMt))*pow2(Msq)) + pow2(Dmglst1)*(180*Dmsqst1*(
        13147 - 90*lmMst1 + 90*lmMt) + (-125277461 - 1327590*lmMt + 60*lmMst1*(
        -2303 + 5955*lmMt) + 153000*pow2(lmMst1) + 151200*pow2(lmMt))*pow2(Msq)
        ) - 750*pow2(Mgl)*(5*Dmsqst1*(79 + 204*lmMst1 + 12*lmMt) + (-4*lmMst1*(
        131 + 264*lmMt) + 264*pow2(lmMst1) + 21*(783 + 14*lmMt + 30*pow2(lmMt))
        )*pow2(Msq)))*pow2(Mst2) + 7840*(-15*Dmsqst1*(50*Dmglst1*(146 - 15*
        lmMst1 + 15*lmMt)*Mgl + (8014 - 435*lmMst1 + 435*lmMt)*pow2(Dmglst1) +
        225*(14 + 6*lmMt - 6*lmMst1*(3 + lmMt) + 3*(pow2(lmMst1) + pow2(lmMt)))
        *pow2(Mgl)) + pow2(Msq)*(-75*Dmglst1*Mgl*(-15571 + 852*lmMt + lmMst1*(-
        508 + 48*lmMt) + 978*pow2(lmMst1) + 612*pow2(lmMt)) - pow2(Dmglst1)*(-
        3223748 + 2205*lmMt + 6*lmMst1*(-16808 + 4605*lmMt) + 79695*pow2(
        lmMst1) + 33750*pow2(lmMt)) + 75*pow2(Mgl)*(623 - 555*lmMt - 3*(221 +
        129*lmMt)*pow2(lmMst1) - 486*pow2(lmMt) + lmMst1*(1066 + 843*lmMt +
        756*pow2(lmMt)) + 252*pow3(lmMst1) - 621*pow3(lmMt))))*pow4(Mst2)) -
        245*pow2(Msq)*pow2(MuSUSY)*(16*pow2(Dmglst1)*(-9*Dmst12*(-3891491 -
        27200*lmMst1 + 960*lmMt + 19200*pow2(lmMst1))*(Dmst12 - pow2(Mst2)) -
        50*(345581 + 4896*lmMst1 + 96*lmMt - 3456*pow2(lmMst1))*pow4(Mst2)) +
        400*Dmglst1*Mgl*(Dmst12*(403559 + 384*(lmMst1 + lmMt) - 4608*pow2(
        lmMst1))*(Dmst12 - pow2(Mst2)) - 24*(9631 + 16*lmMst1 + 48*lmMt - 192*
        pow2(lmMst1))*pow4(Mst2)) - 15*pow2(Mgl)*(pow2(Dmst12)*(852541 + 9216*
        lmMst1 - 10240*lmMt + 6144*pow2(lmMst1)) + 160*Dmst12*(11389 - 704*
        lmMst1 + 192*lmMt - 384*pow2(lmMst1))*pow2(Mst2) + 1920*(349 - 56*
        lmMst1 + 24*lmMt - 32*pow2(lmMst1))*pow4(Mst2)))))))) + pow2(Mt)*(
        3969000*oneLoopFlag*pow2(Mst1)*pow3(Mgl)*pow4(Msq)*(s2t*pow2(Dmst12)*(-
        12*Mt*MuSUSY*pow2(Sbeta) + s2t*Tbeta*(pow2(MuSUSY)*(-1 + pow2(Sbeta)) +
        12*pow2(Mst1)*pow2(Sbeta))) + 12*Mt*pow2(Sbeta)*(2*Dmst12*MuSUSY*s2t*
        pow2(Mst2) + Mt*Tbeta*(pow2(Dmst12) - 2*Dmst12*pow2(Mst2) + 4*(-lmMst1
        + lmMt)*pow4(Mst2)))) + 16*Al4p*xDmglst1*pow2(Msq)*pow3(Dmglst1)*(
        17640*Mst1*twoLoopFlag*pow2(Msq)*(-(pow2(Dmst12)*(8*((476 - 90*lmMst1 +
        15*lmMt)*MuSUSY - 3*(51 + 10*lmMst1 - 10*lmMt)*Mst1*Tbeta)*pow2(Mt)*
        pow2(Sbeta) + 5*Mst1*pow2(s2t)*(10*(-13 + 12*lmMst1)*Tbeta*pow2(MuSUSY)
        *(-1 + pow2(Sbeta)) + 3*Mst1*((77 - 30*lmMst1)*MuSUSY + 4*(-71 + 60*
        lmMst1)*Mst1*Tbeta)*pow2(Sbeta)) + 4*Mt*s2t*(130*Tbeta*pow2(MuSUSY)*(-1
        + pow2(Sbeta)) + 3*Mst1*(120*MuSUSY + (277 - 155*lmMst1 + 5*lmMt)*Mst1*
        Tbeta)*pow2(Sbeta)))) + 2*Dmst12*Mt*pow2(Mst2)*(Mt*((2231 - 990*lmMst1
        + 90*lmMt)*MuSUSY + 12*(-301 + 290*lmMst1 + 10*lmMt)*Mst1*Tbeta)*pow2(
        Sbeta) + s2t*(-5*(-77 + 30*lmMst1)*Tbeta*pow2(MuSUSY)*(-1 + pow2(Sbeta)
        ) + (60*(71 - 60*lmMst1)*Mst1*MuSUSY + 8*(481 - 240*lmMst1 + 15*lmMt)*
        Tbeta*pow2(Mst1))*pow2(Sbeta))) + 8*((977 - 480*lmMst1 + 30*lmMt)*
        MuSUSY + 3*(-509 + 510*lmMst1 + 90*lmMt)*Mst1*Tbeta)*pow2(Mt)*pow2(
        Sbeta)*pow4(Mst2)) + Al4p*threeLoopFlag*(4*Mst1*MuSUSY*pow2(Sbeta)*(
        2940*Dmsqst1*(-20*Dmst12*Mt*((557 + 120*lmMst1 - 120*lmMt)*Mt + 3660*
        Mst1*s2t)*pow2(Mst2) + 5*pow2(Dmst12)*(14640*Mst1*Mt*s2t + 4*(557 +
        120*lmMst1 - 120*lmMt)*pow2(Mt) + 375*pow2(Mst1)*pow2(s2t)) + 24*(463 -
        135*lmMst1 + 135*lmMt)*pow2(Mt)*pow4(Mst2)) - pow2(Msq)*(2*Dmst12*Mt*(
        7350*Mst1*s2t*(548999 + 10980*lmMst1 + 288*lmMt - 19080*pow2(lmMst1)) +
        Mt*(49723877243 + 60*lmMst1*(4936063 - 389970*lmMt) + 57352680*lmMt -
        342543600*pow2(lmMst1) - 3175200*pow2(lmMt)))*pow2(Mst2) - 15*pow2(
        Dmst12)*(490*Mst1*Mt*s2t*(-611423 + 9984*lmMst1 - 768*lmMt + 23040*
        pow2(lmMst1)) + (5753390765 + 580*lmMst1*(79969 - 1932*lmMt) + 7091364*
        lmMt - 35700000*pow2(lmMst1) - 282240*pow2(lmMt))*pow2(Mt) + 49*(
        31025111 + 290880*lmMst1 - 251100*pow2(lmMst1))*pow2(Mst1)*pow2(s2t)) +
        392*(122282257 + 60*lmMst1*(8318 - 3885*lmMt) + 479550*lmMt - 1351800*
        pow2(lmMst1) - 21600*pow2(lmMt))*pow2(Mt)*pow4(Mst2))) + Tbeta*(245*
        pow2(MuSUSY)*(2*Dmst12*Mt*(30000*Dmsqst1*Mst1*s2t + (4*Mst1*s2t*(
        31025111 + 290880*lmMst1 - 251100*pow2(lmMst1)) + Mt*(59957863 +
        480000*lmMst1 - 26880*lmMt - 230400*pow2(lmMst1)))*pow2(Msq))*pow2(
        Mst2) + pow2(Dmst12)*(-60000*Dmsqst1*Mst1*Mt*s2t + pow2(Msq)*(2*Mst1*
        Mt*s2t*(-99874229 - 1352280*lmMst1 + 633600*pow2(lmMst1)) + 2*(-
        59957863 - 480000*lmMst1 + 26880*lmMt + 230400*pow2(lmMst1))*pow2(Mt) +
        15*(3044017 + 27472*lmMst1 - 48480*pow2(lmMst1))*pow2(Mst1)*pow2(s2t)))
        - 24*(-3877891 - 46400*lmMst1 + 960*lmMt + 19200*pow2(lmMst1))*pow2(
        Msq)*pow2(Mt)*pow4(Mst2)) + pow2(Sbeta)*(2*Dmst12*Mst1*Mt*s2t*(11760*
        Dmsqst1*(4*Dmst12*(-2729 + 105*lmMst1 - 105*lmMt)*pow2(Mst1) + 6*(791 -
        270*lmMst1 + 270*lmMt)*pow2(Mst1)*pow2(Mst2) + 625*Dmst12*pow2(MuSUSY)
        - 625*pow2(Mst2)*pow2(MuSUSY)) - pow2(Msq)*(196*pow2(Mst2)*(8*(61021241
        + 15*lmMst1*(20521 - 7770*lmMt) + 250575*lmMt - 675900*pow2(lmMst1) -
        10800*pow2(lmMt))*pow2(Mst1) + 5*(31025111 + 290880*lmMst1 - 251100*
        pow2(lmMst1))*pow2(MuSUSY)) + Dmst12*(4*(25774874431 - 37697520*lmMt +
        600*lmMst1*(311999 + 37149*lmMt) - 77590800*pow2(lmMst1) + 1058400*
        pow2(lmMt))*pow2(Mst1) + 245*(-99874229 - 1352280*lmMst1 + 633600*pow2(
        lmMst1))*pow2(MuSUSY)))) - 3675*pow2(Dmst12)*pow2(Mst1)*(8*(14640*
        Dmsqst1 + (548855 + 10980*lmMst1 + 288*lmMt - 19080*pow2(lmMst1))*pow2(
        Msq))*pow2(Mst1) + (3044017 + 27472*lmMst1 - 48480*pow2(lmMst1))*pow2(
        Msq)*pow2(MuSUSY))*pow2(s2t) - 2*pow2(Mt)*(47040*Dmsqst1*pow2(Mst1)*(3*
        Dmst12*(3401 + 105*lmMst1 - 105*lmMt)*(Dmst12 - pow2(Mst2)) + (-19241 +
        420*lmMst1 - 420*lmMt)*pow4(Mst2)) + pow2(Msq)*(Dmst12*pow2(Mst2)*(16*(
        2524164367 + 8198205*lmMst1 + 19269705*lmMt - 2463300*lmMst1*lmMt +
        27997200*pow2(lmMst1) - 1058400*pow2(lmMt))*pow2(Mst1) - 245*(-59957863
        - 480000*lmMst1 + 26880*lmMt + 230400*pow2(lmMst1))*pow2(MuSUSY)) +
        pow2(Dmst12)*(4*(-7672052891 - 10289580*lmMt + 900*lmMst1*(15649 +
        11284*lmMt) + 98506800*pow2(lmMst1) + 4233600*pow2(lmMt))*pow2(Mst1) +
        245*(-59957863 - 480000*lmMst1 + 26880*lmMt + 230400*pow2(lmMst1))*
        pow2(MuSUSY)) - 588*(4*(-21126629 - 218510*lmMt + 20*lmMst1*(-28958 +
        5055*lmMt) + 194100*pow2(lmMst1) + 91800*pow2(lmMt))*pow2(Mst1) + 5*(-
        3877891 - 46400*lmMst1 + 960*lmMt + 19200*pow2(lmMst1))*pow2(MuSUSY))*
        pow4(Mst2)))))))) + 23520*Mst1*threeLoopFlag*xDmsqst1*pow2(Al4p)*pow2(
        Dmsqst1)*(4500*Mst1*Tbeta*pow2(Mt)*pow3(Mgl)*(24*shiftst2*(pow2(Dmst12)
        + 2*(-2 + lmMst1)*pow2(Mst2)*(Dmst12 + pow2(Mst2)))*pow2(Mt)*pow2(
        Sbeta) + Dmst12*shiftst2*pow2(s2t)*(-(Dmst12*(-1 + 2*lmMst1)*pow2(
        MuSUSY)*(-1 + pow2(Sbeta))) + 12*Dmst12*(-2 + lmMst1)*pow2(Mst1)*pow2(
        Sbeta) - 4*(-2 + lmMst1)*pow2(Mst2)*(pow2(MuSUSY)*(-1 + pow2(Sbeta)) +
        12*pow2(Mst1)*pow2(Sbeta))) - Dmst12*shiftst1*pow2(s2t)*(Dmst12*(-5 +
        2*lmMst1)*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 36*Dmst12*(-2 + lmMst1)*
        pow2(Mst1)*pow2(Sbeta) - 4*(-2 + lmMst1)*pow2(Mst2)*(pow2(MuSUSY)*(-1 +
        pow2(Sbeta)) + 12*pow2(Mst1)*pow2(Sbeta))) + 48*(-2 + lmMst1)*shiftst1*
        pow2(Mt)*pow2(Sbeta)*pow4(Mst2)) + 8*xDmglst1*pow2(Mt)*pow3(Dmglst1)*(-
        2*Dmst12*Mt*pow2(Mst2)*(2*Mt*(5*(557 + 120*lmMst1 - 120*lmMt)*MuSUSY -
        6*(3401 + 105*lmMst1 - 105*lmMt)*Mst1*Tbeta)*pow2(Sbeta) + s2t*(625*
        Tbeta*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + Mst1*(36600*MuSUSY + (-7121 +
        870*lmMst1 - 870*lmMt)*Mst1*Tbeta)*pow2(Sbeta))) + pow2(Dmst12)*((4*(5*
        (557 + 120*lmMst1 - 120*lmMt)*MuSUSY - 6*(3401 + 105*lmMst1 - 105*lmMt)
        *Mst1*Tbeta)*pow2(Mt) + 75*(25*MuSUSY - 488*Mst1*Tbeta)*pow2(Mst1)*
        pow2(s2t))*pow2(Sbeta) + 2*Mt*s2t*(625*Tbeta*pow2(MuSUSY)*(-1 + pow2(
        Sbeta)) + Mst1*(36600*MuSUSY + (-13291 - 330*lmMst1 + 330*lmMt)*Mst1*
        Tbeta)*pow2(Sbeta))) + 4*((3778 - 435*lmMst1 + 435*lmMt)*MuSUSY + 25*(
        1225 - 39*lmMst1 + 39*lmMt)*Mst1*Tbeta)*pow2(Mt)*pow2(Sbeta)*pow4(Mst2)
        ) + Mgl*(5*Tbeta*pow2(Mt)*(5*Dmst12*s2t*(-30*Dmst12*(7 + 24*lmMst1)*
        Mst1*s2t*pow2(Mgl)*pow2(MuSUSY) - 80*Dmst12*Mt*(5*Dmglst1*Mgl + 5*pow2(
        Dmglst1) + 21*pow2(Mgl))*pow2(MuSUSY) + 80*Mt*(5*Dmglst1*Mgl + 5*pow2(
        Dmglst1) + 21*pow2(Mgl))*pow2(Mst2)*pow2(MuSUSY) - 3*Dmst12*Mst1*s2t*((
        3760*Dmglst1*Mgl + 3832*pow2(Dmglst1) - 5*(229 + 288*lmMst1)*pow2(Mgl))
        *pow2(Mst1) - 10*(7 + 24*lmMst1)*pow2(Mgl)*pow2(MuSUSY))*pow2(Sbeta)) +
        4*Mt*pow2(Sbeta)*(4*Dmst12*s2t*pow2(Mst1)*(50*pow2(Mgl)*(Dmst12*(23 -
        6*lmMst1 + 6*lmMt) + (-85 + 12*lmMst1 - 12*lmMt)*pow2(Mst2)) + Dmglst1*
        Mgl*(Dmst12*(1041 - 90*lmMst1 + 90*lmMt) - 25*(91 + 6*lmMst1 - 6*lmMt)*
        pow2(Mst2)) + pow2(Dmglst1)*(9*Dmst12*(-149 + 10*lmMst1 - 10*lmMt) + (
        107 - 330*lmMst1 + 330*lmMt)*pow2(Mst2))) + 100*Dmst12*s2t*(5*Dmglst1*
        Mgl + 5*pow2(Dmglst1) + 21*pow2(Mgl))*(Dmst12 - pow2(Mst2))*pow2(
        MuSUSY) - (4*Mst1*Mt*pow2(Dmglst1)*(3*Dmst12*(13147 - 90*lmMst1 + 90*
        lmMt)*(Dmst12 - pow2(Mst2)) + (-64033 + 3360*lmMst1 - 3360*lmMt)*pow4(
        Mst2)))/5. - 8*Dmglst1*Mgl*Mst1*Mt*(9*Dmst12*(423 - 20*lmMst1 + 20*
        lmMt)*(Dmst12 - pow2(Mst2)) + 25*(-226 + 21*lmMst1 - 21*lmMt)*pow4(
        Mst2)) + Mst1*Mt*pow2(Mgl)*(6*(419 - 60*lmMst1 + 60*lmMt)*pow2(Dmst12)
        - 225*Dmst12*(49 + 46*lmMst1 + 2*lmMt)*pow2(Mst2) + 50*(205 + 174*lmMt
        - 6*lmMst1*(101 + 18*lmMt) + 54*(pow2(lmMst1) + pow2(lmMt)))*pow4(Mst2)
        ))) + 50*MuSUSY*pow2(Sbeta)*(60*pow2(Dmst12)*(5*Dmglst1*Mgl + 5*pow2(
        Dmglst1) + 21*pow2(Mgl))*pow2(Mst1)*pow2(Mt)*pow2(s2t) + Mst1*(3*
        Dmst12*s2t*(8*Dmglst1*(479*Dmglst1 + 470*Mgl)*(Dmst12 - pow2(Mst2)) -
        5*pow2(Mgl)*(44*Dmst12 - (229 + 288*lmMst1)*pow2(Mst2)))*pow3(Mt) +
        540*Mt*s2t*pow2(Mgl)*(-8*Dmst12*(-2 + lmMst1)*shiftst2*pow2(Mst2)*pow2(
        Mt) + pow2(Dmst12)*(-4*shiftst2*pow2(Mt) - (-2 + lmMst1)*(shiftst1 -
        shiftst2)*pow2(Mst1)*pow2(s2t)) + 8*(-2 + lmMst1)*(shiftst1 - shiftst2)
        *pow2(Mt)*pow4(Mst2))) + (16*(25*pow2(Mgl)*(Dmst12*(65 - 6*lmMst1 + 6*
        lmMt)*(Dmst12 - pow2(Mst2)) + (-91 + 12*lmMst1 - 12*lmMt)*pow4(Mst2)) +
        Dmglst1*Mgl*(Dmst12*(557 + 120*lmMst1 - 120*lmMt)*(Dmst12 - pow2(Mst2))
        - 25*(44 + 3*lmMst1 - 3*lmMt)*pow4(Mst2)) + pow2(Dmglst1)*(Dmst12*(557
        + 120*lmMst1 - 120*lmMt)*(Dmst12 - pow2(Mst2)) + (136 - 165*lmMst1 +
        165*lmMt)*pow4(Mst2)))*pow4(Mt))/5.)))) + 11025*threeLoopFlag*pow2(
        Al4p)*(384*z2*pow2(Mst1)*pow3(Mgl)*(pow2(Mst2)*(-50*Mt*xDmsqst1*pow2(
        Dmsqst1)*(-(pow2(Dmst12)*pow2(s2t)*(-(Mt*(shiftst1 + shiftst2)*Tbeta*
        pow2(MuSUSY)*(-1 + pow2(Sbeta))) - 3*(MuSUSY*s2t*(shiftst1 - shiftst2)
        + 2*Mt*(3*shiftst1 - shiftst2)*Tbeta)*pow2(Mst1)*pow2(Sbeta))) + 2*
        Dmst12*Mt*pow2(Mst2)*(12*Mt*MuSUSY*s2t*shiftst2*pow2(Sbeta) - 12*
        shiftst2*Tbeta*pow2(Mt)*pow2(Sbeta) - (shiftst1 - shiftst2)*Tbeta*pow2(
        s2t)*(pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 12*pow2(Mst1)*pow2(Sbeta))) -
        24*(MuSUSY*s2t*(shiftst1 - shiftst2) + Mt*(shiftst1 + shiftst2)*Tbeta)*
        pow2(Mt)*pow2(Sbeta)*pow4(Mst2)) + pow2(Msq)*(50*Tbeta*(Dmsqst1 + pow2(
        Msq))*pow2(Mt)*(-(pow2(Dmst12)*pow2(s2t)*((shiftst1 + shiftst2)*pow2(
        MuSUSY)*(-1 + pow2(Sbeta)) + 6*(3*shiftst1 - shiftst2)*pow2(Mst1)*pow2(
        Sbeta))) + 2*Dmst12*pow2(Mst2)*((shiftst1 - shiftst2)*pow2(MuSUSY)*
        pow2(s2t)*(-1 + pow2(Sbeta)) + 12*(shiftst2*pow2(Mt) + (shiftst1 -
        shiftst2)*pow2(Mst1)*pow2(s2t))*pow2(Sbeta)) + 24*(shiftst1 + shiftst2)
        *pow2(Mt)*pow2(Sbeta)*pow4(Mst2)) + 15*MuSUSY*pow2(Sbeta)*(10*(Dmsqst1
        + pow2(Msq))*(-8*Dmst12*s2t*shiftst2*pow2(Mst2)*pow3(Mt) + Mt*(-
        shiftst1 + shiftst2)*pow2(Dmst12)*pow2(Mst1)*pow3(s2t) + 8*s2t*(
        shiftst1 - shiftst2)*pow3(Mt)*pow4(Mst2)) + shiftst3*pow2(Msq)*(-(Mt*
        pow2(Dmst12)*pow2(Mst1)*pow3(s2t)) + 8*s2t*pow3(Mt)*(pow2(Dmst12) -
        Dmst12*pow2(Mst2) + pow4(Mst2))))) + 5*shiftst3*Tbeta*pow4(Msq)*(-(
        Dmst12*pow2(Mt)*pow2(s2t)*((3*Dmst12 - 2*pow2(Mst2))*pow2(MuSUSY)*(-1 +
        pow2(Sbeta)) + 6*pow2(Mst1)*(7*Dmst12 - 4*pow2(Mst2))*pow2(Sbeta))) +
        24*pow2(Sbeta)*(pow2(Dmst12) - Dmst12*pow2(Mst2) + pow4(Mst2))*pow4(Mt)
        )) - xDmst12*pow3(Dmst12)*(50*Dmsqst1*(Dmsqst1*xDmsqst1 + pow2(Msq))*(-
        4*shiftst2*pow2(Mt)*(Tbeta*pow2(MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta)) +
        6*Mt*(MuSUSY*s2t - Mt*Tbeta)*pow2(Sbeta)) + pow2(Sbeta)*(2*Mt*(MuSUSY*
        s2t*(-2*shiftst1 + shiftst2) + 2*Mt*(-4*shiftst1 + shiftst2)*Tbeta)*
        pow2(Mst1)*pow2(s2t) + (shiftst1 - shiftst2)*Tbeta*pow4(Mst1)*pow4(s2t)
        )) + pow4(Msq)*(2*pow2(Mt)*(-((15*shiftst1 + 10*shiftst2 + 9*shiftst3)*
        Tbeta*pow2(MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta))) + 60*Mt*shiftst3*(
        MuSUSY*s2t + Mt*Tbeta)*pow2(Sbeta)) + pow2(Sbeta)*(-5*Mt*(MuSUSY*s2t*(
        40*shiftst1 - 20*shiftst2 + 7*shiftst3) + 2*Mt*(80*shiftst1 - 20*
        shiftst2 + 29*shiftst3)*Tbeta)*pow2(Mst1)*pow2(s2t) + 5*(10*shiftst1 -
        10*shiftst2 + shiftst3)*Tbeta*pow4(Mst1)*pow4(s2t))))) - 5*Mt*z3*(Mt*
        pow2(Mst2)*(8*Dmst12*Mgl*Mst1*s2t*Tbeta*pow2(Msq)*(Dmst12*Mst1*s2t*(12*
        (420*Dmsqst1*pow2(Mgl) + (2898*Dmglst1*Mgl - 3295*pow2(Dmglst1) + 660*
        pow2(Mgl))*pow2(Msq))*pow2(Mst1) + (1260*Dmsqst1*pow2(Mgl) - (50485*
        Dmglst1*Mgl + 200115*pow2(Dmglst1) + 2574*pow2(Mgl))*pow2(Msq))*pow2(
        MuSUSY)) - 2*Mt*pow2(Msq)*(12*pow2(Mst1)*(Dmst12*(22100*Dmglst1*Mgl +
        86801*pow2(Dmglst1) + 3290*pow2(Mgl)) + 8*(1963*Dmglst1*Mgl + 9451*
        pow2(Dmglst1) + 106*pow2(Mgl))*pow2(Mst2)) - (5*Dmst12*(-21081*Dmglst1*
        Mgl + 21081*pow2(Dmglst1) - 3457*pow2(Mgl)) - 2*(37582*Dmglst1*Mgl +
        142987*pow2(Dmglst1) + 3012*pow2(Mgl))*pow2(Mst2))*pow2(MuSUSY)))*pow2(
        Sbeta) - 5040*xDmsqst1*pow2(Dmsqst1)*pow2(Mst1)*pow3(Mgl)*(6*Dmst12*Mt*
        (-5*MuSUSY*s2t + 6*Mt*Tbeta)*pow2(Mst2)*pow2(Sbeta) + Tbeta*pow2(
        Dmst12)*pow2(s2t)*(-2*pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 15*pow2(Mst1)*
        pow2(Sbeta)) + 24*Tbeta*pow2(Mt)*pow2(Sbeta)*pow4(Mst2)) - 8*xDmglst1*
        pow3(Dmglst1)*pow4(Msq)*(6*Dmst12*Mt*pow2(Mst2)*(147351*Mt*Tbeta*pow2(
        MuSUSY)*(-1 + pow2(Sbeta)) + 48*(3977*MuSUSY*s2t + 7908*Mt*Tbeta)*pow2(
        Mst1)*pow2(Sbeta) + 4*Mst1*MuSUSY*(76518*MuSUSY*s2t*Tbeta*(-1 + pow2(
        Sbeta)) + 500581*Mt*pow2(Sbeta)) + 966992*s2t*Tbeta*pow2(Sbeta)*pow3(
        Mst1)) + pow2(Dmst12)*(-2*Mst1*Mt*MuSUSY*(737647*MuSUSY*s2t*Tbeta*(-1 +
        pow2(Sbeta)) + 5205992*Mt*pow2(Sbeta)) + pow2(Mst1)*(349745*Tbeta*pow2(
        MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta)) - 24*Mt*(-22939*MuSUSY*s2t +
        77043*Mt*Tbeta)*pow2(Sbeta)) - 24*s2t*(114777*MuSUSY*s2t - 258833*Mt*
        Tbeta)*pow2(Sbeta)*pow3(Mst1) + Tbeta*(-884106*pow2(Mt)*pow2(MuSUSY)*(-
        1 + pow2(Sbeta)) + 572688*pow2(s2t)*pow2(Sbeta)*pow4(Mst1))) + 24*pow2(
        Mt)*(28665*Tbeta*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 4*Mst1*(60437*MuSUSY
        + 28466*Mst1*Tbeta)*pow2(Sbeta))*pow4(Mst2)) + Mgl*pow2(Msq)*(-3*Tbeta*
        pow2(Mt)*pow2(Sbeta)*(20160*Dmsqst1*Dmst12*pow2(Mgl)*pow2(Mst1)*(Dmst12
        + 2*pow2(Mst2)) + pow2(Msq)*(pow2(Mgl)*(160*Dmst12*pow2(Mst2)*(2052*
        pow2(Mst1) + 425*pow2(MuSUSY)) + pow2(Dmst12)*(29472*pow2(Mst1) +
        31963*pow2(MuSUSY)) + 512*((-214 + 24*lmMst1 - 24*lmMt)*pow2(Mst1) +
        47*pow2(MuSUSY))*pow4(Mst2)) - 48*Dmglst1*Mgl*(pow2(Dmst12)*(6970*pow2(
        Mst1) + 8213*pow2(MuSUSY)) - Dmst12*pow2(Mst2)*(30912*pow2(Mst1) +
        8213*pow2(MuSUSY)) - 8*(3356*pow2(Mst1) + 583*pow2(MuSUSY))*pow4(Mst2))
        - 48*pow2(Dmglst1)*(7*pow2(Dmst12)*(7481*pow2(Mst1) + 4095*pow2(MuSUSY)
        ) - 5*Dmst12*pow2(Mst2)*(15248*pow2(Mst1) + 5733*pow2(MuSUSY)) - 2*(
        37624*pow2(Mst1) + 7047*pow2(MuSUSY))*pow4(Mst2)))) - Tbeta*pow2(
        MuSUSY)*(8*pow2(Dmst12)*(1260*Dmsqst1*pow2(Mgl) - (50485*Dmglst1*Mgl +
        200115*pow2(Dmglst1) + 2574*pow2(Mgl))*pow2(Msq))*pow2(Mst1)*pow2(s2t)
        + Mt*pow2(Msq)*(16*Dmst12*Mst1*s2t*(5*Dmst12*(-21081*Dmglst1*Mgl +
        21081*pow2(Dmglst1) - 3457*pow2(Mgl)) - 2*(37582*Dmglst1*Mgl + 142987*
        pow2(Dmglst1) + 3012*pow2(Mgl))*pow2(Mst2)) + 3*Mt*(48*Dmglst1*Mgl*(
        8213*pow2(Dmst12) - 8213*Dmst12*pow2(Mst2) - 4664*pow4(Mst2)) + 432*
        pow2(Dmglst1)*(3185*pow2(Dmst12) - 3185*Dmst12*pow2(Mst2) - 1566*pow4(
        Mst2)) - pow2(Mgl)*(31963*pow2(Dmst12) + 68000*Dmst12*pow2(Mst2) +
        24064*pow4(Mst2))))) + 16*Mst1*MuSUSY*pow2(Sbeta)*(6*Dmst12*Mst1*Mt*
        s2t*(105*Dmsqst1*pow2(Mgl)*(7*Dmst12 + 8*pow2(Mst2)) + pow2(Msq)*(3*
        pow2(Mgl)*(403*Dmst12 + 440*pow2(Mst2)) - 2*pow2(Dmglst1)*(32498*Dmst12
        + 3295*pow2(Mst2)) + Dmglst1*(-7642*Dmst12*Mgl + 5796*Mgl*pow2(Mst2))))
        + pow2(Msq)*(3*pow2(Dmst12)*(37582*Dmglst1*Mgl + 142987*pow2(Dmglst1) +
        3012*pow2(Mgl))*pow2(Mst1)*pow2(s2t) + pow2(Mt)*(4*pow2(Dmglst1)*(
        286982*pow2(Dmst12) - 487227*Dmst12*pow2(Mst2) - 226824*pow4(Mst2)) -
        pow2(Mgl)*(51181*pow2(Dmst12) + 49656*Dmst12*pow2(Mst2) + 10176*pow4(
        Mst2)) - 4*Dmglst1*Mgl*(86833*pow2(Dmst12) + 113412*Dmst12*pow2(Mst2) +
        47112*pow4(Mst2))))))) + 2*xDmst12*pow3(Dmst12)*(pow3(Mgl)*(-2520*
        Dmsqst1*pow2(Mst1)*(Dmsqst1*xDmsqst1*(Mt*Tbeta*pow2(s2t)*(4*pow2(
        MuSUSY)*(-1 + pow2(Sbeta)) + 15*pow2(Mst1)*pow2(Sbeta)) - 36*Tbeta*
        pow2(Sbeta)*pow3(Mt) + 2*MuSUSY*pow2(Sbeta)*(15*s2t*pow2(Mt) + pow2(
        Mst1)*pow3(s2t))) + pow2(Msq)*(Mt*Tbeta*pow2(s2t)*(4*pow2(MuSUSY)*(-1 +
        pow2(Sbeta)) + pow2(Mst1)*pow2(Sbeta)) - 48*Tbeta*pow2(Sbeta)*pow3(Mt)
        + 2*MuSUSY*pow2(Sbeta)*(22*s2t*pow2(Mt) + pow2(Mst1)*pow3(s2t)))) +
        pow4(Msq)*(16*Mst1*MuSUSY*pow2(Mt)*(20297*MuSUSY*s2t*Tbeta*(-1 + pow2(
        Sbeta)) + 76009*Mt*pow2(Sbeta)) - Mt*pow2(Mst1)*(37669*Tbeta*pow2(
        MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta)) - 4*Mt*(-32783*MuSUSY*s2t + 33933*
        Mt*Tbeta)*pow2(Sbeta)) + 4*s2t*(Mt*(33783*MuSUSY*s2t - 23402*Mt*Tbeta)
        + 18*Mst1*s2t*(143*MuSUSY*s2t - 37*Mt*Tbeta))*pow2(Sbeta)*pow3(Mst1) +
        Tbeta*(197889*pow2(MuSUSY)*(-1 + pow2(Sbeta))*pow3(Mt) + 24096*pow2(
        Sbeta)*pow3(s2t)*pow5(Mst1)))) + 4*Dmglst1*pow4(Msq)*(pow2(Mgl)*(4*
        Mst1*MuSUSY*pow2(Mt)*(142987*MuSUSY*s2t*Tbeta*(-1 + pow2(Sbeta)) +
        574156*Mt*pow2(Sbeta)) + 2*Mt*pow2(Mst1)*(50485*Tbeta*pow2(MuSUSY)*
        pow2(s2t)*(-1 + pow2(Sbeta)) + 12*Mt*(4744*MuSUSY*s2t + 12729*Mt*Tbeta)
        *pow2(Sbeta)) + s2t*(Mt*(90723*MuSUSY*s2t - 164264*Mt*Tbeta) + Mst1*
        s2t*(50485*MuSUSY*s2t - 80628*Mt*Tbeta))*pow2(Sbeta)*pow3(Mst1) + 86*
        Tbeta*(-1719*pow2(MuSUSY)*(-1 + pow2(Sbeta))*pow3(Mt) + 874*pow2(Sbeta)
        *pow3(s2t)*pow5(Mst1))) + Dmglst1*Mgl*(-8*Mst1*MuSUSY*pow2(Mt)*(-18791*
        MuSUSY*s2t*Tbeta*(-1 + pow2(Sbeta)) + 86737*Mt*pow2(Sbeta)) + 6*Mt*
        pow2(Mst1)*(66705*Tbeta*pow2(MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta)) + (
        273164*Mt*MuSUSY*s2t - 85482*Tbeta*pow2(Mt))*pow2(Sbeta)) + s2t*(-3*Mt*
        s2t*(391379*MuSUSY + 116812*Mst1*Tbeta) + 4379080*Tbeta*pow2(Mt) +
        200115*Mst1*MuSUSY*pow2(s2t))*pow2(Sbeta)*pow3(Mst1) + Tbeta*(-515970*
        pow2(MuSUSY)*(-1 + pow2(Sbeta))*pow3(Mt) + 285974*pow2(Sbeta)*pow3(s2t)
        *pow5(Mst1))) + xDmglst1*pow2(Dmglst1)*(-4*Mst1*MuSUSY*pow2(Mt)*(
        278539*MuSUSY*s2t*Tbeta*(-1 + pow2(Sbeta)) + 2202506*Mt*pow2(Sbeta)) +
        2*Mt*pow2(Mst1)*(349745*Tbeta*pow2(MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta))
        - 24*Mt*(-46801*MuSUSY*s2t + 29595*Mt*Tbeta)*pow2(Sbeta)) + s2t*(-
        4967589*Mt*MuSUSY*s2t + 297420*Mst1*Mt*s2t*Tbeta + 16623976*Tbeta*pow2(
        Mt) + 349745*Mst1*MuSUSY*pow2(s2t))*pow2(Sbeta)*pow3(Mst1) + Tbeta*(-
        884106*pow2(MuSUSY)*(-1 + pow2(Sbeta))*pow3(Mt) + 918216*pow2(Sbeta)*
        pow3(s2t)*pow5(Mst1))))))) + 4*xDmst12*pow3(Dmst12)*(5880*Mst1*
        threeLoopFlag*xDmsqst1*pow2(Al4p)*pow2(Dmsqst1)*(8*Dmglst1*Mt*(5*pow2(
        Mgl)*(-2*s2t*pow2(Mt)*(125*Tbeta*pow2(MuSUSY)*(-1 + pow2(Sbeta)) +
        Mst1*(7050*MuSUSY + (-193 - 330*lmMst1 + 330*lmMt)*Mst1*Tbeta)*pow2(
        Sbeta)) + pow2(Sbeta)*(150*Mt*(-5*MuSUSY + 94*Mst1*Tbeta)*pow2(Mst1)*
        pow2(s2t) - 4*((557 + 120*lmMst1 - 120*lmMt)*MuSUSY + 9*(-423 + 20*
        lmMst1 - 20*lmMt)*Mst1*Tbeta)*pow3(Mt) + 125*Tbeta*pow3(s2t)*pow4(Mst1)
        )) + xDmglst1*pow2(Dmglst1)*(-2*s2t*pow2(Mt)*(625*Tbeta*pow2(MuSUSY)*(-
        1 + pow2(Sbeta)) + (36600*Mst1*MuSUSY - 3*(6487 + 510*lmMst1 - 510*
        lmMt)*Tbeta*pow2(Mst1))*pow2(Sbeta)) + pow2(Sbeta)*(150*Mt*(-25*MuSUSY
        + 488*Mst1*Tbeta)*pow2(Mst1)*pow2(s2t) - 4*(5*(557 + 120*lmMst1 - 120*
        lmMt)*MuSUSY - 6*(3401 + 105*lmMst1 - 105*lmMt)*Mst1*Tbeta)*pow3(Mt) +
        625*Tbeta*pow3(s2t)*pow4(Mst1))) + Dmglst1*Mgl*(-50*s2t*pow2(Mt)*(25*
        Tbeta*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + Mst1*(1437*MuSUSY - 5*(103 + 6*
        lmMst1 - 6*lmMt)*Mst1*Tbeta)*pow2(Sbeta)) + pow2(Sbeta)*(150*Mt*(-25*
        MuSUSY + 479*Mst1*Tbeta)*pow2(Mst1)*pow2(s2t) + (-20*(557 + 120*lmMst1
        - 120*lmMt)*MuSUSY + 6*(13147 - 90*lmMst1 + 90*lmMt)*Mst1*Tbeta)*pow3(
        Mt) + 625*Tbeta*pow3(s2t)*pow4(Mst1)))) + 5*pow3(Mgl)*(-75*Mst1*pow2(
        Mt)*pow2(s2t)*(4*(7 + 6*shiftst1 - 24*lmMst1*(-1 + shiftst2) + 30*
        shiftst2)*Tbeta*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 3*Mst1*(112*MuSUSY +
        Mst1*(91 + 240*shiftst1 + 32*lmMst1*(3 - 4*shiftst1 + shiftst2))*Tbeta)
        *pow2(Sbeta)) - 150*s2t*(56*Tbeta*pow2(MuSUSY)*(-1 + pow2(Sbeta)) +
        Mst1*(3*MuSUSY*(47 + 96*(lmMst1 + shiftst2 - lmMst1*shiftst2)) - 208*
        Mst1*Tbeta)*pow2(Sbeta))*pow3(Mt) + pow2(Sbeta)*(-150*Mt*(MuSUSY*(7 +
        108*shiftst1 - 72*shiftst2 + 24*lmMst1*(1 - 2*shiftst1 + shiftst2)) -
        28*Mst1*Tbeta)*pow3(Mst1)*pow3(s2t) - 4*(200*(65 - 6*lmMst1 + 6*lmMt)*
        MuSUSY - 3*Mst1*(1999 - 90*lmMt + 90*lmMst1*(41 - 40*shiftst2) + 3600*
        shiftst2)*Tbeta)*pow4(Mt) - 1800*(-2 + lmMst1)*(shiftst1 - shiftst2)*
        Tbeta*pow4(s2t)*pow5(Mst1)))) + Al4p*Mgl*pow2(Msq)*(-35280*Mst1*Mt*
        twoLoopFlag*pow2(Msq)*(pow2(Mgl)*(25*Mst1*Mt*pow2(s2t)*((1 + 6*lmMst1)*
        Tbeta*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 12*Mst1*(MuSUSY + 3*lmMst1*
        MuSUSY + (1 + 12*lmMst1)*Mst1*Tbeta)*pow2(Sbeta)) + 100*s2t*pow2(Mt)*(-
        ((5 + 6*lmMst1)*Tbeta*pow2(MuSUSY)*(-1 + pow2(Sbeta))) + 2*Mst1*(-12*(-
        1 + lmMst1)*MuSUSY + (1 + 3*lmMst1 + 9*lmMt)*Mst1*Tbeta)*pow2(Sbeta)) +
        pow2(Sbeta)*(-4*(50*(5 + 6*lmMst1)*MuSUSY + (137 - 330*lmMst1 - 270*
        lmMt)*Mst1*Tbeta)*pow3(Mt) + 75*(MuSUSY - 2*(1 + 2*lmMst1)*Mst1*Tbeta)*
        pow3(Mst1)*pow3(s2t))) + 2*Dmglst1*(Mgl*(-25*Mst1*Mt*pow2(s2t)*(8*(-1 +
        3*lmMst1)*Tbeta*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 3*Mst1*((5 - 6*
        lmMst1)*MuSUSY + 6*(-1 + 4*lmMst1)*Mst1*Tbeta)*pow2(Sbeta)) - 2*s2t*
        pow2(Mt)*(30*lmMst1*(5*Tbeta*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + Mst1*(
        60*MuSUSY - 11*Mst1*Tbeta)*pow2(Sbeta)) + Tbeta*(-125*pow2(MuSUSY)*(-1
        + pow2(Sbeta)) + 2*(131 - 135*lmMt)*pow2(Mst1)*pow2(Sbeta))) + pow2(
        Sbeta)*(4*(6*(17 - 30*lmMst1 + 5*lmMt)*MuSUSY + (47 + 870*lmMst1 + 30*
        lmMt)*Mst1*Tbeta)*pow3(Mt) + 25*((4 - 12*lmMst1)*MuSUSY + (5 - 6*
        lmMst1)*Mst1*Tbeta)*pow3(Mst1)*pow3(s2t))) + Dmglst1*(-15*Mst1*Mt*pow2(
        s2t)*(10*(-5 + 6*lmMst1)*Tbeta*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + Mst1*(
        (85 - 30*lmMst1)*MuSUSY + (-137 + 180*lmMst1)*Mst1*Tbeta)*pow2(Sbeta))
        + 10*s2t*pow2(Mt)*(-5*(1 + 6*lmMst1)*Tbeta*pow2(MuSUSY)*(-1 + pow2(
        Sbeta)) + 6*Mst1*((58 - 90*lmMst1)*MuSUSY + (2 + 15*lmMst1 + 5*lmMt)*
        Mst1*Tbeta)*pow2(Sbeta)) + pow2(Sbeta)*(-2*(3*(259 + 90*lmMst1 + 10*
        lmMt)*MuSUSY + (1547 - 2730*lmMst1 + 30*lmMt)*Mst1*Tbeta)*pow3(Mt) +
        25*(3*(5 - 6*lmMst1)*MuSUSY + (11 - 6*lmMst1)*Mst1*Tbeta)*pow3(Mst1)*
        pow3(s2t))))) - Al4p*threeLoopFlag*(196*Mst1*Mt*MuSUSY*pow2(Sbeta)*(
        Dmglst1*Mgl*(2400*Dmsqst1*Mt*(7050*Mst1*Mt*s2t + 2*(557 + 120*lmMst1 -
        120*lmMt)*pow2(Mt) + 375*pow2(Mst1)*pow2(s2t)) - pow2(Msq)*(3600*Mst1*
        s2t*(12383 + 4128*lmMst1 + 80*lmMt - 1260*pow2(lmMst1))*pow2(Mt) + 225*
        Mt*(284641 + 8696*lmMst1 + 1680*pow2(lmMst1))*pow2(Mst1)*pow2(s2t) + 8*
        (193364399 + 90000*lmMt - 300*lmMst1*(3781 + 582*lmMt) - 2005200*pow2(
        lmMst1) - 43200*pow2(lmMt))*pow3(Mt) + 375*(84209 + 1264*lmMst1 - 240*
        pow2(lmMst1))*pow3(Mst1)*pow3(s2t))) + 15*pow2(Mgl)*(500*Dmsqst1*(6*
        Mst1*s2t*(173 - 144*lmMst1*(-1 + shiftst2) + 216*shiftst2)*pow2(Mt) +
        504*Mt*pow2(Mst1)*pow2(s2t) + 16*(65 - 6*lmMst1 + 6*lmMt)*pow3(Mt) + 3*
        (7 + 72*shiftst1 - 36*shiftst2 + 24*lmMst1*(1 - 2*shiftst1 + shiftst2))
        *pow3(Mst1)*pow3(s2t)) - pow2(Msq)*(15*Mst1*s2t*(-102747 + 640*lmMt +
        6720*shiftst3 - 32*lmMst1*(331 + 90*shiftst3) + 13888*pow2(lmMst1))*
        pow2(Mt) - 75*Mt*(-20531 + 200*lmMst1 + 1200*pow2(lmMst1))*pow2(Mst1)*
        pow2(s2t) - 4*(-3454599 + 16840*lmMt + 48*lmMst1*(262 + 405*lmMt) +
        46560*pow2(lmMst1))*pow3(Mt) + 50*(1429 - 720*shiftst1 + 360*shiftst2 -
        234*shiftst3 + 2*lmMst1*(-227 + 720*shiftst1 - 360*shiftst2 + 126*
        shiftst3) + 24*pow2(lmMst1))*pow3(Mst1)*pow3(s2t))) + pow2(Dmglst1)*(
        2400*Dmsqst1*Mt*(7185*Mst1*Mt*s2t + 2*(557 + 120*lmMst1 - 120*lmMt)*
        pow2(Mt) + 375*pow2(Mst1)*pow2(s2t)) + pow2(Msq)*(-7200*Mst1*s2t*(
        143196 + 2546*lmMst1 + 92*lmMt - 4785*pow2(lmMst1))*pow2(Mt) + 225*Mt*(
        3532083 + 36328*lmMst1 - 47760*pow2(lmMst1))*pow2(Mst1)*pow2(s2t) + 16*
        (29818901 + 258300*lmMt + 30*lmMst1*(35963 + 2190*lmMt) - 239400*pow2(
        lmMst1) - 10800*pow2(lmMt))*pow3(Mt) + 75*(-1732531 - 16896*lmMst1 +
        24840*pow2(lmMst1))*pow3(Mst1)*pow3(s2t)))) + Tbeta*(-245*(24*pow2(
        Dmglst1)*(10000*Dmsqst1*Mst1*s2t + (3*Mt*(3891491 + 27200*lmMst1 - 960*
        lmMt - 19200*pow2(lmMst1)) + 100*Mst1*s2t*(-32829 + 1852*lmMst1 + 660*
        pow2(lmMst1)))*pow2(Msq)) + 40*Dmglst1*Mgl*(6000*Dmsqst1*Mst1*s2t + (5*
        Mt*(403559 + 384*(lmMst1 + lmMt) - 4608*pow2(lmMst1)) + 6*Mst1*s2t*(-
        1282471 + 7264*lmMst1 + 18120*pow2(lmMst1)))*pow2(Msq)) + 15*pow2(Mgl)*
        (67200*Dmsqst1*Mst1*s2t + (80*Mst1*s2t*(-36863 + 80*lmMst1 + 552*pow2(
        lmMst1)) + Mt*(-1763661 + 47104*lmMst1 - 5120*lmMt + 24576*pow2(lmMst1)
        ))*pow2(Msq)))*pow2(MuSUSY)*pow3(Mt) + 3675*pow2(Mst1)*(-(pow2(Mt)*
        pow2(s2t)*((2400*Dmsqst1*(7 + 24*lmMst1)*pow2(Mgl) + (40*Dmglst1*Mgl*(-
        84209 - 1264*lmMst1 + 240*pow2(lmMst1)) + 8*pow2(Dmglst1)*(-1732531 -
        16896*lmMst1 + 24840*pow2(lmMst1)) + 5*(70121 + 2208*lmMst1 - 432*pow2(
        lmMst1))*pow2(Mgl))*pow2(Msq))*pow2(MuSUSY) - 1440*shiftst1*pow2(Mgl)*(
        80*Dmsqst1*(3 - 2*lmMst1)*pow2(Mst1) + (1 - 2*lmMst1)*pow2(Msq)*(80*
        pow2(Mst1) + 3*pow2(MuSUSY)))*pow2(Sbeta))) + 144*pow2(Mgl)*(-(
        shiftst3*pow2(Msq)*(2*pow2(Mt)*pow2(s2t)*((-29 + 18*lmMst1)*pow2(
        MuSUSY)*(-1 + pow2(Sbeta)) + 5*(-83 + 58*lmMst1)*pow2(Mst1)*pow2(Sbeta)
        ) - 80*(-7 + 3*lmMst1)*pow2(Sbeta)*pow4(Mt) + 5*(1 - 2*lmMst1)*pow2(
        Sbeta)*pow4(Mst1)*pow4(s2t))) - 10*(3*(1 - 2*lmMst1)*shiftst1*pow2(Msq)
        *pow2(Mt)*pow2(MuSUSY)*pow2(s2t) - (1 - 2*lmMst1)*shiftst2*pow2(Msq)*
        pow2(s2t)*(2*pow2(Mt)*(pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 10*pow2(Mst1)*
        pow2(Sbeta)) + 5*pow2(s2t)*pow2(Sbeta)*pow4(Mst1)) + 5*shiftst1*(
        Dmsqst1*(3 - 2*lmMst1) + (1 - 2*lmMst1)*pow2(Msq))*pow2(Sbeta)*pow4(
        Mst1)*pow4(s2t) + 5*Dmsqst1*(3 - 2*lmMst1)*shiftst2*(4*pow2(Mt)*pow2(
        s2t)*(-(pow2(MuSUSY)*(-1 + pow2(Sbeta))) + pow2(Mst1)*pow2(Sbeta)) +
        pow2(Sbeta)*(24*pow4(Mt) - pow4(Mst1)*pow4(s2t)))))) + Mt*pow2(Sbeta)*(
        -392*Mst1*s2t*pow2(Mt)*(15*pow2(Mgl)*(2000*Dmsqst1*(2*(13 + 6*lmMst1 -
        6*lmMt)*pow2(Mst1) - 21*pow2(MuSUSY)) - pow2(Msq)*((558619 + 76160*lmMt
        - 224*lmMst1*(1219 + 60*lmMt) + 123840*pow2(lmMst1) + 86400*pow2(lmMt))
        *pow2(Mst1) + 50*(-36863 + 80*lmMst1 + 552*pow2(lmMst1))*pow2(MuSUSY)))
        + 10*pow2(Dmglst1)*(3000*Dmsqst1*(84*pow2(Mst1) - 5*pow2(MuSUSY)) +
        pow2(Msq)*((148185343 + 1333716*lmMst1 + 170460*lmMt + 87840*lmMst1*
        lmMt - 1376640*pow2(lmMst1) - 43200*pow2(lmMt))*pow2(Mst1) - 150*(-
        32829 + 1852*lmMst1 + 660*pow2(lmMst1))*pow2(MuSUSY))) - 2*Dmglst1*Mgl*
        (600*Dmsqst1*(6*(47 - 30*lmMst1 + 30*lmMt)*pow2(Mst1) + 125*pow2(
        MuSUSY)) + pow2(Msq)*((28188929 - 143100*lmMt - 3780*lmMst1*(549 + 80*
        lmMt) + 1389600*pow2(lmMst1) + 388800*pow2(lmMt))*pow2(Mst1) + 75*(-
        1282471 + 7264*lmMst1 + 18120*pow2(lmMst1))*pow2(MuSUSY)))) - 3675*Mt*
        pow2(Mst1)*(8*pow2(Dmglst1)*(114960*Dmsqst1*pow2(Mst1) + pow2(Msq)*(3*(
        -612347 - 7436*lmMst1 - 608*lmMt + 19320*pow2(lmMst1))*pow2(Mst1) + (
        1732531 + 16896*lmMst1 - 24840*pow2(lmMst1))*pow2(MuSUSY))) - 5*pow2(
        Mgl)*(120*Dmsqst1*((-137 + 288*lmMst1)*pow2(Mst1) + 4*(7 + 24*lmMst1)*
        pow2(MuSUSY)) + pow2(Msq)*(24*(2785 - 304*lmMst1 + 384*lmMt + 768*pow2(
        lmMst1))*pow2(Mst1) + (70121 + 2208*lmMst1 - 432*pow2(lmMst1))*pow2(
        MuSUSY))) + 8*Dmglst1*Mgl*(112800*Dmsqst1*pow2(Mst1) + pow2(Msq)*(24*(-
        20017 + 1203*lmMst1 - 200*lmMt + 2250*pow2(lmMst1))*pow2(Mst1) + 5*(
        84209 + 1264*lmMst1 - 240*pow2(lmMst1))*pow2(MuSUSY))))*pow2(s2t) - (
        40*pow2(Dmglst1)*(7056*Dmsqst1*(13147 - 90*lmMst1 + 90*lmMt)*pow2(Mst1)
        + 2*(-700000759 + 6327384*lmMt + 12*lmMst1*(-88589 + 185010*lmMt) +
        85997520*pow2(lmMst1) + 423360*pow2(lmMt))*pow2(Msq)*pow2(Mst1) + 441*(
        -3891491 - 27200*lmMst1 + 960*lmMt + 19200*pow2(lmMst1))*pow2(Msq)*
        pow2(MuSUSY)) + 49*pow2(Mgl)*(9600*Dmsqst1*(3268 + 2805*lmMst1 - 105*
        lmMt)*pow2(Mst1) + pow2(Msq)*((83430364 - 8607840*lmMt + 480*lmMst1*(
        36107 + 13380*lmMt) - 9273600*pow2(lmMst1) - 6652800*pow2(lmMt))*pow2(
        Mst1) + 75*(1763661 - 47104*lmMst1 + 5120*lmMt - 24576*pow2(lmMst1))*
        pow2(MuSUSY))) + 392*Dmglst1*Mgl*(21600*Dmsqst1*(423 - 20*lmMst1 + 20*
        lmMt)*pow2(Mst1) + pow2(Msq)*(12*(9598037 + 92280*lmMt + 20*lmMst1*(-
        11207 + 270*lmMt) + 246000*pow2(lmMst1) - 14400*pow2(lmMt))*pow2(Mst1)
        - 125*(403559 + 384*(lmMst1 + lmMt) - 4608*pow2(lmMst1))*pow2(MuSUSY)))
        )*pow3(Mt) - 29400*(pow2(Dmglst1)*(1000*Dmsqst1 + (1286791 + 5936*
        lmMst1 - 18120*pow2(lmMst1))*pow2(Msq)) + 10*Dmglst1*Mgl*(100*Dmsqst1 +
        (33261 - 532*lmMst1 - 660*pow2(lmMst1))*pow2(Msq)) + 20*pow2(Mgl)*(210*
        Dmsqst1 + (1361 + 10*lmMst1 + 54*pow2(lmMst1))*pow2(Msq)))*pow3(s2t)*
        pow5(Mst1))))) + Mt*(992250*oneLoopFlag*pow2(Mst1)*pow3(Mgl)*(-(Mt*
        Tbeta*pow2(s2t)*(pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 18*pow2(Mst1)*pow2(
        Sbeta))) - pow2(Sbeta)*(-8*MuSUSY*s2t*pow2(Mt) + 8*Tbeta*pow3(Mt) +
        MuSUSY*pow2(Mst1)*pow3(s2t)))*pow4(Msq) + 4*Al4p*xDmglst1*pow2(Msq)*
        pow3(Dmglst1)*(-17640*Mst1*twoLoopFlag*pow2(Msq)*(-5*Mst1*Mt*pow2(s2t)*
        (20*(-13 + 12*lmMst1)*Tbeta*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 3*Mst1*((
        129 - 30*lmMst1)*MuSUSY + 4*(-83 + 60*lmMst1)*Mst1*Tbeta)*pow2(Sbeta))
        + 2*s2t*pow2(Mt)*(-15*(9 + 10*lmMst1)*Tbeta*pow2(MuSUSY)*(-1 + pow2(
        Sbeta)) + 2*Mst1*(30*(47 - 60*lmMst1)*MuSUSY + (106 + 285*lmMst1 + 15*
        lmMt)*Mst1*Tbeta)*pow2(Sbeta)) + pow2(Sbeta)*(-2*((1577 + 270*lmMst1 +
        30*lmMt)*MuSUSY + 12*(199 - 310*lmMst1 + 10*lmMt)*Mst1*Tbeta)*pow3(Mt)
        + 5*((130 - 120*lmMst1)*MuSUSY + (77 - 30*lmMst1)*Mst1*Tbeta)*pow3(
        Mst1)*pow3(s2t))) + Al4p*threeLoopFlag*(11760*Dmsqst1*Mst1*(-2*s2t*
        pow2(Mt)*(625*Tbeta*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + (36600*Mst1*
        MuSUSY - 2*(8543 + 390*lmMst1 - 390*lmMt)*Tbeta*pow2(Mst1))*pow2(Sbeta)
        ) + pow2(Sbeta)*(150*Mt*(-25*MuSUSY + 488*Mst1*Tbeta)*pow2(Mst1)*pow2(
        s2t) - 4*(5*(557 + 120*lmMst1 - 120*lmMt)*MuSUSY - 6*(3401 + 105*lmMst1
        - 105*lmMt)*Mst1*Tbeta)*pow3(Mt) + 625*Tbeta*pow3(s2t)*pow4(Mst1))) +
        pow2(Msq)*(-4*Mst1*MuSUSY*pow2(Mt)*(-245*MuSUSY*s2t*Tbeta*(-37824007 -
        770520*lmMst1 + 131400*pow2(lmMst1))*(-1 + pow2(Sbeta)) + 8*Mt*(
        9144246058 + 12254445*lmMt + 90*lmMst1*(1109907 + 18305*lmMt) -
        48239100*pow2(lmMst1) - 264600*pow2(lmMt))*pow2(Sbeta)) + 6*Mt*pow2(
        Mst1)*(-1225*Tbeta*(-3044017 - 27472*lmMst1 + 48480*pow2(lmMst1))*pow2(
        MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta)) + 8*Mt*(2450*MuSUSY*s2t*(580211 +
        498*lmMst1 + 528*lmMt - 21060*pow2(lmMst1)) + Mt*Tbeta*(-874574719 +
        9416610*lmMt + 10*lmMst1*(1016017 + 174300*lmMt) + 51500400*pow2(
        lmMst1) + 705600*pow2(lmMt)))*pow2(Sbeta)) + s2t*(3675*Mst1*s2t*(
        MuSUSY*s2t*(3044017 + 27472*lmMst1 - 48480*pow2(lmMst1)) + 4*Mt*Tbeta*(
        486671 + 31944*lmMst1 - 192*lmMt - 15120*pow2(lmMst1))) + Mt*(735*
        MuSUSY*s2t*(-223974673 - 2515800*lmMst1 + 1638000*pow2(lmMst1)) + 4*Mt*
        Tbeta*(137797425107 + 35209020*lmMt + 300*lmMst1*(3595111 + 92568*lmMt)
        - 690681600*pow2(lmMst1) - 2116800*pow2(lmMt))))*pow2(Sbeta)*pow3(Mst1)
        + 490*Tbeta*((-59957863 - 480000*lmMst1 + 26880*lmMt + 230400*pow2(
        lmMst1))*pow2(MuSUSY)*(-1 + pow2(Sbeta))*pow3(Mt) + 2*(31025111 +
        290880*lmMst1 - 251100*pow2(lmMst1))*pow2(Sbeta)*pow3(s2t)*pow5(Mst1)))
        )))))/(1.90512e8*Tbeta*pow2(Mst1)*pow2(Sbeta)*pow3(Mgl)*pow4(Msq)*pow6(
        Mst2));
}

/**
 * @return The off-diagonal (1, 2) = (2, 1) matrix element of the Higgs mass matrix as a double for the hierarchy 'H3'
 */
double H3::getS12() const {
   return -(MuSUSY*(Mt*(1984500*Dmst12*oneLoopFlag*s2t*pow2(Mst1)*pow3(Mgl)*pow4(
        Msq)*(-2*Dmst12*Mt*(MuSUSY*s2t + 6*Mt*Tbeta)*pow2(Mst2) + pow2(Dmst12)*
        (2*Mt*MuSUSY*s2t + 8*Tbeta*pow2(Mt) - Tbeta*pow2(Mst1)*pow2(s2t)) + 24*
        Tbeta*pow2(Mt)*pow4(Mst2)) + 8*Al4p*xDmglst1*pow2(Msq)*pow3(Dmglst1)*(
        17640*Mst1*twoLoopFlag*pow2(Msq)*(Mt*pow2(Dmst12)*pow2(Mst2)*(80*Mt*
        s2t*(13*MuSUSY - 18*Mst1*Tbeta) + 8*(-476 + 90*lmMst1 - 15*lmMt)*Tbeta*
        pow2(Mt) + 5*Mst1*(20*(-13 + 12*lmMst1)*MuSUSY + 3*(-77 + 30*lmMst1)*
        Mst1*Tbeta)*pow2(s2t)) + pow3(Dmst12)*(-60*s2t*((9 + 10*lmMst1)*MuSUSY
        + 2*(47 - 60*lmMst1)*Mst1*Tbeta)*pow2(Mt) - 5*Mst1*Mt*(40*(-13 + 12*
        lmMst1)*MuSUSY + 9*(-43 + 10*lmMst1)*Mst1*Tbeta)*pow2(s2t) + 2*(1577 +
        270*lmMst1 + 30*lmMt)*Tbeta*pow3(Mt) + 50*(-13 + 12*lmMst1)*Tbeta*pow3(
        Mst1)*pow3(s2t)) - 2*Dmst12*(10*(77 - 30*lmMst1)*MuSUSY*s2t + (-2231 +
        990*lmMst1 - 90*lmMt)*Mt*Tbeta + 60*(-71 + 60*lmMst1)*Mst1*s2t*Tbeta)*
        pow2(Mt)*pow4(Mst2) + 8*(977 - 480*lmMst1 + 30*lmMt)*Tbeta*pow3(Mt)*
        pow6(Mst2)) + Al4p*threeLoopFlag*(490*Mt*MuSUSY*(-(pow2(Dmst12)*pow2(
        Mst2)*(60000*Dmsqst1*Mst1*Mt*s2t + pow2(Msq)*(2*Mst1*Mt*s2t*(99874229 +
        1352280*lmMst1 - 633600*pow2(lmMst1)) + (119915726 + 960000*lmMst1 -
        53760*lmMt - 460800*pow2(lmMst1))*pow2(Mt) + 15*(-3044017 - 27472*
        lmMst1 + 48480*pow2(lmMst1))*pow2(Mst1)*pow2(s2t)))) + 2*((30000*
        Dmsqst1*Mst1*Mt*s2t + pow2(Msq)*(2*Mst1*Mt*s2t*(37824007 + 770520*
        lmMst1 - 131400*pow2(lmMst1)) + (59957863 + 480000*lmMst1 - 26880*lmMt
        - 230400*pow2(lmMst1))*pow2(Mt) + 15*(-3044017 - 27472*lmMst1 + 48480*
        pow2(lmMst1))*pow2(Mst1)*pow2(s2t)))*pow3(Dmst12) + Dmst12*Mt*(30000*
        Dmsqst1*Mst1*s2t + (4*Mst1*s2t*(31025111 + 290880*lmMst1 - 251100*pow2(
        lmMst1)) + Mt*(59957863 + 480000*lmMst1 - 26880*lmMt - 230400*pow2(
        lmMst1)))*pow2(Msq))*pow4(Mst2)) + 24*(3877891 + 46400*lmMst1 - 960*
        lmMt - 19200*pow2(lmMst1))*pow2(Msq)*pow2(Mt)*pow6(Mst2)) - Mst1*Tbeta*
        (735*Mt*pow2(Dmst12)*pow2(Mst1)*(30000*Dmsqst1*(2*Dmst12 - pow2(Mst2))
        + pow2(Msq)*(Dmst12*(223974673 + 2515800*lmMst1 - 1638000*pow2(lmMst1))
        + 4*(-31025111 - 290880*lmMst1 + 251100*pow2(lmMst1))*pow2(Mst2)))*
        pow2(s2t) + 3675*(-3044017 - 27472*lmMst1 + 48480*pow2(lmMst1))*pow2(
        Msq)*pow3(Dmst12)*pow3(Mst1)*pow3(s2t) + 29400*Dmst12*Mst1*s2t*pow2(Mt)
        *(29280*Dmsqst1*(pow2(Dmst12) - Dmst12*pow2(Mst2) + pow4(Mst2)) + pow2(
        Msq)*(4*pow2(Dmst12)*(-580211 - 498*lmMst1 - 528*lmMt + 21060*pow2(
        lmMst1)) + Dmst12*(611423 - 9984*lmMst1 + 768*lmMt - 23040*pow2(lmMst1)
        )*pow2(Mst2) + 2*(548999 + 10980*lmMst1 + 288*lmMt - 19080*pow2(lmMst1)
        )*pow4(Mst2))) + 4*pow3(Mt)*(11760*Dmsqst1*(5*Dmst12*(557 + 120*lmMst1
        - 120*lmMt)*(pow2(Dmst12) - Dmst12*pow2(Mst2) + pow4(Mst2)) + 6*(-463 +
        135*lmMst1 - 135*lmMt)*pow6(Mst2)) + pow2(Msq)*(-15*pow2(Dmst12)*(
        5753390765 + 580*lmMst1*(79969 - 1932*lmMt) + 7091364*lmMt - 35700000*
        pow2(lmMst1) - 282240*pow2(lmMt))*pow2(Mst2) + 8*(9144246058 +
        12254445*lmMt + 90*lmMst1*(1109907 + 18305*lmMt) - 48239100*pow2(
        lmMst1) - 264600*pow2(lmMt))*pow3(Dmst12) + 2*Dmst12*(49723877243 + 60*
        lmMst1*(4936063 - 389970*lmMt) + 57352680*lmMt - 342543600*pow2(lmMst1)
        - 3175200*pow2(lmMt))*pow4(Mst2) + 392*(122282257 + 60*lmMst1*(8318 -
        3885*lmMt) + 479550*lmMt - 1351800*pow2(lmMst1) - 21600*pow2(lmMt))*
        pow6(Mst2))))))) + 1470*Al4p*Mgl*pow2(Msq)*(48*Mst1*Mt*twoLoopFlag*
        pow2(Msq)*(-25*s2t*pow3(Dmst12)*(8*MuSUSY*(Dmglst1*(-5 + 6*lmMst1)*Mgl
        + (1 + 6*lmMst1)*pow2(Dmglst1) + (5 + 6*lmMst1)*pow2(Mgl))*pow2(Mt) +
        Tbeta*(8*Dmglst1*(1 - 3*lmMst1)*Mgl + 6*(5 - 6*lmMst1)*pow2(Dmglst1) +
        3*pow2(Mgl))*pow2(s2t)*pow3(Mst1)) + Mt*(50*Dmst12*MuSUSY*s2t*(Mst1*
        s2t*pow2(Dmst12)*(16*Dmglst1*(1 - 3*lmMst1)*Mgl + (60 - 72*lmMst1)*
        pow2(Dmglst1) + (1 + 6*lmMst1)*pow2(Mgl)) + 4*Dmst12*Mt*(6*pow2(
        Dmglst1) + pow2(Mgl))*pow2(Mst2) - Dmst12*Mst1*s2t*(8*Dmglst1*(1 - 3*
        lmMst1)*Mgl + 6*(5 - 6*lmMst1)*pow2(Dmglst1) + 3*pow2(Mgl))*pow2(Mst2)
        - 4*Mt*(Dmglst1*(5 - 6*lmMst1)*Mgl + (11 - 6*lmMst1)*pow2(Dmglst1) - 3*
        (1 + 2*lmMst1)*pow2(Mgl))*pow4(Mst2)) + 2*Tbeta*(3*pow2(Dmglst1)*(-(
        pow2(Dmst12)*pow2(Mst2)*(70*Mst1*Mt*s2t + 8*(92 - 30*lmMst1 + 5*lmMt)*
        pow2(Mt) + 25*(11 - 6*lmMst1)*pow2(Mst1)*pow2(s2t))) + (40*(-29 + 45*
        lmMst1)*Mst1*Mt*s2t + 2*(259 + 90*lmMst1 + 10*lmMt)*pow2(Mt) + 25*(17 -
        6*lmMst1)*pow2(Mst1)*pow2(s2t))*pow3(Dmst12) + 2*Dmst12*Mt*((477 - 330*
        lmMst1 + 30*lmMt)*Mt + 50*(13 - 18*lmMst1)*Mst1*s2t)*pow4(Mst2) + 400*(
        4 - 3*lmMst1)*pow2(Mt)*pow6(Mst2)) - 25*pow2(Mgl)*(-(pow2(Dmst12)*pow2(
        Mst2)*(24*(1 - 3*lmMst1)*Mst1*Mt*s2t + 2*(5 + 6*lmMst1 + 6*lmMt)*pow2(
        Mt) + 9*(1 + 2*lmMst1)*pow2(Mst1)*pow2(s2t))) + (-48*(-1 + lmMst1)*
        Mst1*Mt*s2t - 4*(5 + 6*lmMst1)*pow2(Mt) + 6*(1 + 3*lmMst1)*pow2(Mst1)*
        pow2(s2t))*pow3(Dmst12) + 4*Dmst12*Mt*(2*(5 + 6*lmMst1 + 3*lmMt)*Mt -
        9*(1 + 2*lmMst1)*Mst1*s2t)*pow4(Mst2) + 72*(1 + lmMst1 + lmMt)*pow2(Mt)
        *pow6(Mst2)) + Dmglst1*Mgl*(-3*pow2(Dmst12)*(pow2(Mst2)*(100*Mst1*Mt*
        s2t - 2*(-41 + 90*lmMst1 + 10*lmMt)*pow2(Mt) + 25*(5 - 6*lmMst1)*pow2(
        Mst1)*pow2(s2t)) + Dmst12*(-1200*lmMst1*Mst1*Mt*s2t + 8*(17 - 30*lmMst1
        + 5*lmMt)*pow2(Mt) + 25*(-5 + 6*lmMst1)*pow2(Mst1)*pow2(s2t))) + 300*
        Dmst12*Mt*((3 - 6*lmMst1)*Mt + 2*(1 - 6*lmMst1)*Mst1*s2t)*pow4(Mst2) -
        200*(-7 + 15*lmMst1 + 3*lmMt)*pow2(Mt)*pow6(Mst2))))) + Al4p*
        threeLoopFlag*((MuSUSY*pow2(Mt)*(8*pow2(Dmglst1)*(-3*pow2(Dmst12)*pow2(
        Mst2)*(20000*Dmsqst1*Mst1*Mt*s2t + pow2(Msq)*(10*Mst1*Mt*s2t*(958501 +
        24456*lmMst1 - 11520*pow2(lmMst1)) + 6*(3891491 + 27200*lmMst1 - 960*
        lmMt - 19200*pow2(lmMst1))*pow2(Mt) + 5*(-1732531 - 16896*lmMst1 +
        24840*pow2(lmMst1))*pow2(Mst1)*pow2(s2t))) + 6*((10000*Dmsqst1*Mst1*Mt*
        s2t + pow2(Msq)*(100*Mst1*Mt*s2t*(-32829 + 1852*lmMst1 + 660*pow2(
        lmMst1)) + 3*(3891491 + 27200*lmMst1 - 960*lmMt - 19200*pow2(lmMst1))*
        pow2(Mt) + 5*(-1732531 - 16896*lmMst1 + 24840*pow2(lmMst1))*pow2(Mst1)*
        pow2(s2t)))*pow3(Dmst12) + Dmst12*Mt*(10000*Dmsqst1*Mst1*s2t + (3*Mt*(
        3891491 + 27200*lmMst1 - 960*lmMt - 19200*pow2(lmMst1)) + 10*Mst1*s2t*(
        1286791 + 5936*lmMst1 - 18120*pow2(lmMst1)))*pow2(Msq))*pow4(Mst2)) +
        100*(345581 + 4896*lmMst1 + 96*lmMt - 3456*pow2(lmMst1))*pow2(Msq)*
        pow2(Mt)*pow6(Mst2)) + 40*Dmglst1*Mgl*(-(pow2(Dmst12)*pow2(Mst2)*(
        12000*Dmsqst1*Mst1*Mt*s2t + pow2(Msq)*(6*Mst1*Mt*s2t*(-949861 + 1944*
        lmMst1 + 11520*pow2(lmMst1)) + 10*(403559 + 384*(lmMst1 + lmMt) - 4608*
        pow2(lmMst1))*pow2(Mt) + 15*(-84209 - 1264*lmMst1 + 240*pow2(lmMst1))*
        pow2(Mst1)*pow2(s2t)))) + 2*(6000*Dmsqst1*Mst1*Mt*s2t + pow2(Msq)*(6*
        Mst1*Mt*s2t*(-1282471 + 7264*lmMst1 + 18120*pow2(lmMst1)) + 5*(403559 +
        384*(lmMst1 + lmMt) - 4608*pow2(lmMst1))*pow2(Mt) + 15*(-84209 - 1264*
        lmMst1 + 240*pow2(lmMst1))*pow2(Mst1)*pow2(s2t)))*pow3(Dmst12) + 10*
        Dmst12*Mt*(1200*Dmsqst1*Mst1*s2t + (Mt*(403559 + 384*(lmMst1 + lmMt) -
        4608*pow2(lmMst1)) - 12*Mst1*s2t*(-33261 + 532*lmMst1 + 660*pow2(
        lmMst1)))*pow2(Msq))*pow4(Mst2) + 240*(9631 + 16*lmMst1 + 48*lmMt -
        192*pow2(lmMst1))*pow2(Msq)*pow2(Mt)*pow6(Mst2)) + 15*pow2(Mgl)*(-(
        pow2(Dmst12)*pow2(Mst2)*(2400*Dmsqst1*Mst1*s2t*(56*Mt + Mst1*s2t*(7 -
        12*lmMst1*(-2 + shiftst1 + shiftst2) + 18*(shiftst1 + shiftst2))) -
        pow2(Msq)*(-240*Mst1*Mt*s2t*(-10473 + 40*lmMst1 + 256*pow2(lmMst1)) + (
        852541 + 9216*lmMst1 - 10240*lmMt + 6144*pow2(lmMst1))*pow2(Mt) + 80*(
        1429 - 454*lmMst1 - 180*(shiftst1 + shiftst2) + 360*lmMst1*(shiftst1 +
        shiftst2) - 126*shiftst3 + 108*lmMst1*shiftst3 + 24*pow2(lmMst1))*pow2(
        Mst1)*pow2(s2t)))) + 2*(-2400*Dmsqst1*Mst1*s2t*(-28*Mt + Mst1*s2t*(-7 +
        24*lmMst1*(-1 + shiftst2) - 36*shiftst2)) + pow2(Msq)*(80*Mst1*Mt*s2t*(
        -36863 + 80*lmMst1 + 552*pow2(lmMst1)) + (-1763661 + 47104*lmMst1 -
        5120*lmMt + 24576*pow2(lmMst1))*pow2(Mt) + (350605 + 4320*shiftst1 +
        2880*shiftst2 + 8352*shiftst3 - 96*lmMst1*(-115 + 90*shiftst1 + 60*
        shiftst2 + 54*shiftst3) - 2160*pow2(lmMst1))*pow2(Mst1)*pow2(s2t)))*
        pow3(Dmst12) + 160*Dmst12*(60*Dmsqst1*Mst1*s2t*(14*Mt + 3*(3 - 2*
        lmMst1)*Mst1*s2t*(shiftst1 - shiftst2)) + pow2(Msq)*(4*Mst1*Mt*s2t*(
        1361 + 10*lmMst1 + 54*pow2(lmMst1)) + (11389 - 704*lmMst1 + 192*lmMt -
        384*pow2(lmMst1))*pow2(Mt) + 18*(1 - 2*lmMst1)*(10*shiftst1 - 10*
        shiftst2 + shiftst3)*pow2(Mst1)*pow2(s2t)))*pow4(Mst2) + 1920*(349 -
        56*lmMst1 + 24*lmMt - 32*pow2(lmMst1))*pow2(Msq)*pow2(Mt)*pow6(Mst2))))
        /3. - 4*Mst1*Tbeta*(15*Dmst12*Mst1*s2t*pow2(Mt)*(Dmst12*Mst1*s2t*(5*
        pow2(Mgl)*(1680*Dmsqst1*(2*Dmst12 - pow2(Mst2)) + pow2(Msq)*(Dmst12*(-
        20531 + 200*lmMst1 + 1200*pow2(lmMst1)) - 8*(1361 + 10*lmMst1 + 54*
        pow2(lmMst1))*pow2(Mst2))) + Dmglst1*Mgl*(2000*Dmsqst1*(2*Dmst12 -
        pow2(Mst2)) - pow2(Msq)*(Dmst12*(284641 + 8696*lmMst1 + 1680*pow2(
        lmMst1)) - 20*(-33261 + 532*lmMst1 + 660*pow2(lmMst1))*pow2(Mst2))) +
        pow2(Dmglst1)*(2000*Dmsqst1*(2*Dmst12 - pow2(Mst2)) + pow2(Msq)*(
        Dmst12*(3532083 + 36328*lmMst1 - 47760*pow2(lmMst1)) + 2*(-1286791 -
        5936*lmMst1 + 18120*pow2(lmMst1))*pow2(Mst2)))) + Mt*(pow2(Mgl)*(100*
        Dmsqst1*((346 + 288*lmMst1)*pow2(Dmst12) - 161*Dmst12*pow2(Mst2) - 24*(
        1 + 12*lmMst1)*pow4(Mst2)) + pow2(Msq)*(pow2(Dmst12)*(102747 + 10592*
        lmMst1 - 640*lmMt - 13888*pow2(lmMst1)) + 20*Dmst12*(-2071 - 296*lmMst1
        + 96*lmMt + 600*pow2(lmMst1))*pow2(Mst2) - 160*(631 - lmMst1 + 36*lmMt
        + 21*pow2(lmMst1))*pow4(Mst2))) + 4*pow2(Dmglst1)*(19160*Dmsqst1*(pow2(
        Dmst12) - Dmst12*pow2(Mst2) + pow4(Mst2)) + pow2(Msq)*(8*pow2(Dmst12)*(
        -143196 - 2546*lmMst1 - 92*lmMt + 4785*pow2(lmMst1)) + Dmst12*(586073 +
        9268*lmMst1 + 448*lmMt - 19200*pow2(lmMst1))*pow2(Mst2) + 2*(-13289 +
        916*lmMst1 - 80*lmMt + 60*pow2(lmMst1))*pow4(Mst2))) + 16*Dmglst1*Mgl*(
        4700*Dmsqst1*(pow2(Dmst12) - Dmst12*pow2(Mst2) + pow4(Mst2)) + pow2(
        Msq)*(pow2(Dmst12)*(-12383 - 4128*lmMst1 - 80*lmMt + 1260*pow2(lmMst1))
        + Dmst12*(17539 + 574*lmMst1 + 160*lmMt - 1920*pow2(lmMst1))*pow2(Mst2)
        + 5*(-4539 + 596*lmMst1 - 48*lmMt + 516*pow2(lmMst1))*pow4(Mst2))))) +
        900*Mst1*Mt*s2t*pow2(Mgl)*(-10*shiftst2*pow2(Dmst12)*(Dmsqst1*(3 - 2*
        lmMst1) + (1 - 2*lmMst1)*pow2(Msq))*pow2(Mst1)*(2*Dmst12 - 3*pow2(Mst2)
        )*pow2(s2t) + shiftst3*pow2(Dmst12)*pow2(Msq)*pow2(Mst1)*(Dmst12*(13 -
        14*lmMst1) + 3*(-1 + 2*lmMst1)*pow2(Mst2))*pow2(s2t) - 240*shiftst2*
        pow2(Mt)*((1 - 2*lmMst1)*pow2(Msq)*(Dmst12 + pow2(Mst2))*pow4(Mst2) +
        Dmsqst1*(-3 + 2*lmMst1)*(pow3(Dmst12) - Dmst12*pow4(Mst2) - pow6(Mst2))
        ) - 8*shiftst3*pow2(Msq)*pow2(Mt)*(6*(-2 + lmMst1)*pow2(Dmst12)*pow2(
        Mst2) + (14 - 6*lmMst1)*pow3(Dmst12) + 3*Dmst12*(3 - 2*lmMst1)*pow4(
        Mst2) + 3*(-1 + 2*lmMst1)*pow6(Mst2))) + Mt*(50*(30*Dmsqst1*(7 + 24*
        lmMst1)*pow2(Mgl) - ((pow2(Dmglst1)*(1732531 + 16896*lmMst1 - 24840*
        pow2(lmMst1)) + 5*Dmglst1*Mgl*(84209 + 1264*lmMst1 - 240*pow2(lmMst1))
        + 10*(1429 - 454*lmMst1 + 24*pow2(lmMst1))*pow2(Mgl))*pow2(Msq))/10.)*
        pow3(Dmst12)*pow3(Mst1)*pow3(s2t) + 9000*Mst1*s2t*shiftst1*pow2(Mgl)*(
        Dmsqst1*(3 - 2*lmMst1) + (1 - 2*lmMst1)*pow2(Msq))*(pow2(Dmst12)*pow2(
        Mst1)*(4*Dmst12 - 3*pow2(Mst2))*pow2(s2t) + 24*pow2(Mt)*pow6(Mst2))) +
        (2*pow4(Mt)*(2*Dmglst1*(1200*Dmglst1*Dmsqst1*(Dmst12*(557 + 120*lmMst1
        - 120*lmMt)*(pow2(Dmst12) - Dmst12*pow2(Mst2) + pow4(Mst2)) + 16*(4 +
        15*lmMst1 - 15*lmMt)*pow6(Mst2)) + 1200*Dmsqst1*Mgl*(Dmst12*(557 + 120*
        lmMst1 - 120*lmMt)*(pow2(Dmst12) - Dmst12*pow2(Mst2) + pow4(Mst2)) +
        50*(26 + 3*lmMst1 - 3*lmMt)*pow6(Mst2)) + Dmglst1*pow2(Msq)*(3*pow2(
        Dmst12)*(-129193181 - 401100*lmMt + 100*lmMst1*(-7351 + 1164*lmMt) +
        1336800*pow2(lmMst1) + 28800*pow2(lmMt))*pow2(Mst2) + 4*(29818901 +
        258300*lmMt + 30*lmMst1*(35963 + 2190*lmMt) - 239400*pow2(lmMst1) -
        10800*pow2(lmMt))*pow3(Dmst12) + 2*Dmst12*(327941741 + 1080*lmMst1*(44
        - 445*lmMt) + 686700*lmMt - 3531600*pow2(lmMst1) - 64800*pow2(lmMt))*
        pow4(Mst2) + 80*(3777727 + 24255*lmMt - 486*lmMst1*(118 + 45*lmMt) -
        52380*pow2(lmMst1))*pow6(Mst2)) - Mgl*pow2(Msq)*(3*pow2(Dmst12)*(-
        39474953 + lmMst1*(52780 - 87600*lmMt) - 3300*lmMt + 319200*pow2(
        lmMst1) + 14400*pow2(lmMt))*pow2(Mst2) - 2*(-193364399 - 90000*lmMt +
        300*lmMst1*(3781 + 582*lmMt) + 2005200*pow2(lmMst1) + 43200*pow2(lmMt))
        *pow3(Dmst12) + 40*Dmst12*(-3746977 - 4005*lmMt + 18*lmMst1*(2711 +
        1215*lmMt) + 52380*pow2(lmMst1))*pow4(Mst2) - 6000*(9961 + 66*lmMt -
        10*lmMst1*(67 + 24*lmMt) + 42*pow2(lmMst1) + 72*pow2(lmMt))*pow6(Mst2))
        ) + 15*pow2(Mgl)*(4000*Dmsqst1*(Dmst12*(65 - 6*lmMst1 + 6*lmMt)*(pow2(
        Dmst12) - Dmst12*pow2(Mst2) + pow4(Mst2)) + 6*(20 - 3*lmMst1 + 3*lmMt)*
        pow6(Mst2)) - pow2(Msq)*(pow2(Dmst12)*(-2284899 + 49840*lmMt - 32*
        lmMst1*(-1793 + 555*lmMt) + 87360*pow2(lmMst1) + 28800*pow2(lmMt))*
        pow2(Mst2) - 2*(-3454599 + 16840*lmMt + 48*lmMst1*(262 + 405*lmMt) +
        46560*pow2(lmMst1))*pow3(Dmst12) - 200*Dmst12*(11697 + 448*lmMst1 +
        330*lmMt - 372*lmMst1*lmMt + 408*pow2(lmMst1) + 288*pow2(lmMt))*pow4(
        Mst2) - 1600*(434 - 83*lmMst1 + 174*lmMt - 66*lmMst1*lmMt + 183*pow2(
        lmMst1) + 108*pow2(lmMt))*pow6(Mst2)))))/15.))) + 735*threeLoopFlag*
        pow2(Al4p)*(-2880*z2*pow2(Mst1)*pow3(Mgl)*(-50*Mt*s2t*xDmsqst1*pow2(
        Dmsqst1)*(s2t*pow2(Dmst12)*(2*Mt*MuSUSY*(shiftst1 + shiftst2) + 3*s2t*(
        -shiftst1 + shiftst2)*Tbeta*pow2(Mst1))*pow2(Mst2) + (-8*Mt*MuSUSY*s2t*
        shiftst2 + 24*shiftst2*Tbeta*pow2(Mt) + 2*(2*shiftst1 - shiftst2)*
        Tbeta*pow2(Mst1)*pow2(s2t))*pow3(Dmst12) - 4*Dmst12*Mt*(MuSUSY*s2t*(
        shiftst1 - shiftst2) + 6*Mt*shiftst2*Tbeta)*pow4(Mst2) + 24*(shiftst1 -
        shiftst2)*Tbeta*pow2(Mt)*pow6(Mst2)) + pow2(Msq)*(2*Dmst12*MuSUSY*pow2(
        Mt)*pow2(s2t)*(2*pow2(Dmst12)*(100*Dmsqst1*shiftst2 + (15*shiftst1 +
        10*shiftst2 + 9*shiftst3)*pow2(Msq)) - 5*Dmst12*(10*Dmsqst1*(shiftst1 +
        shiftst2) + (10*(shiftst1 + shiftst2) + 3*shiftst3)*pow2(Msq))*pow2(
        Mst2) + 10*(10*Dmsqst1*(shiftst1 - shiftst2) + (10*shiftst1 - 10*
        shiftst2 + shiftst3)*pow2(Msq))*pow4(Mst2)) + 5*shiftst3*Tbeta*pow2(
        Msq)*(-(Mt*pow2(Dmst12)*pow2(Mst1)*(7*Dmst12 - 3*pow2(Mst2))*pow3(s2t))
        + 24*s2t*pow3(Mt)*(-(pow2(Dmst12)*pow2(Mst2)) + pow3(Dmst12) + Dmst12*
        pow4(Mst2) - pow6(Mst2))) - 50*Tbeta*(shiftst2*(-(Mt*pow2(Dmst12)*(
        Dmsqst1 + pow2(Msq))*pow2(Mst1)*(2*Dmst12 - 3*pow2(Mst2))*pow3(s2t)) +
        24*s2t*pow3(Mt)*(Dmsqst1*pow3(Dmst12) - (Dmsqst1 + pow2(Msq))*(Dmst12 +
        pow2(Mst2))*pow4(Mst2))) + shiftst1*(Dmsqst1 + pow2(Msq))*(Mt*pow2(
        Dmst12)*pow2(Mst1)*(4*Dmst12 - 3*pow2(Mst2))*pow3(s2t) + 24*s2t*pow3(
        Mt)*pow6(Mst2))))) + 32*Mst1*xDmsqst1*pow2(Dmsqst1)*(4*xDmglst1*pow2(
        Mt)*pow3(Dmglst1)*(2500*Dmst12*Mt*MuSUSY*s2t*(pow2(Dmst12) - Dmst12*
        pow2(Mst2) + pow4(Mst2)) - Tbeta*(1875*pow2(Dmst12)*pow2(Mst1)*(2*
        Dmst12 - pow2(Mst2))*pow2(s2t) + 73200*Dmst12*Mst1*Mt*s2t*(pow2(Dmst12)
        - Dmst12*pow2(Mst2) + pow4(Mst2)) + 4*pow2(Mt)*(5*Dmst12*(557 + 120*
        lmMst1 - 120*lmMt)*(pow2(Dmst12) - Dmst12*pow2(Mst2) + pow4(Mst2)) + (-
        3778 + 435*lmMst1 - 435*lmMt)*pow6(Mst2)))) + Mgl*(250*Dmst12*MuSUSY*
        s2t*pow2(Mt)*(40*Dmglst1*(Dmglst1 + Mgl)*Mt*(pow2(Dmst12) - Dmst12*
        pow2(Mst2) + pow4(Mst2)) + 3*pow2(Mgl)*(2*(28*Mt + Mst1*s2t*(7 + 6*
        shiftst1 - 24*lmMst1*(-1 + shiftst2) + 30*shiftst2))*pow2(Dmst12) -
        Dmst12*(56*Mt + Mst1*s2t*(7 + 30*shiftst1 + 6*shiftst2 - 12*lmMst1*(-2
        + shiftst1 + shiftst2)))*pow2(Mst2) + 8*(7*Mt - 3*(-2 + lmMst1)*Mst1*
        s2t*(shiftst1 - shiftst2))*pow4(Mst2))) - 25*Tbeta*(60*pow2(Dmst12)*(5*
        Dmglst1*Mgl + 5*pow2(Dmglst1) + 21*pow2(Mgl))*pow2(Mst1)*(2*Dmst12 -
        pow2(Mst2))*pow2(Mt)*pow2(s2t) + 3*Dmst12*Mst1*s2t*pow3(Mt)*(8*Dmglst1*
        (479*Dmglst1 + 470*Mgl)*(pow2(Dmst12) - Dmst12*pow2(Mst2) + pow4(Mst2))
        + 5*pow2(Mgl)*(3*(47 + 96*lmMst1)*pow2(Dmst12) + 44*Dmst12*pow2(Mst2) -
        (229 + 288*lmMst1)*pow4(Mst2))) + 15*Mst1*Mt*s2t*pow2(Mgl)*(36*pow2(
        Dmst12)*pow2(Mst2)*(4*shiftst2*pow2(Mt) + (-2 + lmMst1)*(shiftst1 -
        shiftst2)*pow2(Mst1)*pow2(s2t)) + (-288*(-1 + lmMst1)*shiftst2*pow2(Mt)
        + (7 + 108*shiftst1 - 72*shiftst2 + 24*lmMst1*(1 - 2*shiftst1 +
        shiftst2))*pow2(Mst1)*pow2(s2t))*pow3(Dmst12) + 288*Dmst12*(-2 +
        lmMst1)*shiftst2*pow2(Mt)*pow4(Mst2) - 288*(-2 + lmMst1)*(shiftst1 -
        shiftst2)*pow2(Mt)*pow6(Mst2)) + (16*pow4(Mt)*(pow2(Dmglst1)*(Dmst12*(
        557 + 120*lmMst1 - 120*lmMt)*(pow2(Dmst12) - Dmst12*pow2(Mst2) + pow4(
        Mst2)) + (-136 + 165*lmMst1 - 165*lmMt)*pow6(Mst2)) + Dmglst1*Mgl*(
        Dmst12*(557 + 120*lmMst1 - 120*lmMt)*(pow2(Dmst12) - Dmst12*pow2(Mst2)
        + pow4(Mst2)) + 25*(44 + 3*lmMst1 - 3*lmMt)*pow6(Mst2)) + 25*pow2(Mgl)*
        (Dmst12*(65 - 6*lmMst1 + 6*lmMt)*(pow2(Dmst12) - Dmst12*pow2(Mst2) +
        pow4(Mst2)) + (91 - 12*lmMst1 + 12*lmMt)*pow6(Mst2))))/5.))) - 75*Mt*
        z3*(-5040*Dmst12*s2t*xDmsqst1*pow2(Dmsqst1)*pow2(Mst1)*pow3(Mgl)*(2*
        Dmst12*Mt*MuSUSY*s2t*pow2(Mst2) + pow2(Dmst12)*(-4*Mt*MuSUSY*s2t + 15*
        Tbeta*pow2(Mt) + Tbeta*pow2(Mst1)*pow2(s2t)) - 15*Tbeta*pow2(Mt)*pow4(
        Mst2)) + 4*xDmglst1*pow3(Dmglst1)*pow4(Msq)*(-2*Mt*pow2(Dmst12)*pow2(
        Mst2)*(2*Mst1*Mt*s2t*(737647*MuSUSY + 137634*Mst1*Tbeta) + (884106*
        MuSUSY - 5205992*Mst1*Tbeta)*pow2(Mt) - (349745*MuSUSY + 1377324*Mst1*
        Tbeta)*pow2(Mst1)*pow2(s2t)) + pow3(Dmst12)*(8*Mst1*s2t*(278539*MuSUSY
        + 280806*Mst1*Tbeta)*pow2(Mt) - 11*Mt*(127180*MuSUSY + 451599*Mst1*
        Tbeta)*pow2(Mst1)*pow2(s2t) + 4*(442053*MuSUSY - 2202506*Mst1*Tbeta)*
        pow3(Mt) + 349745*Tbeta*pow3(s2t)*pow4(Mst1)) + 12*Dmst12*(Mt*(147351*
        MuSUSY - 1001162*Mst1*Tbeta) + 24*Mst1*s2t*(12753*MuSUSY - 3977*Mst1*
        Tbeta))*pow2(Mt)*pow4(Mst2) + 624*(2205*MuSUSY - 9298*Mst1*Tbeta)*pow3(
        Mt)*pow6(Mst2)) + Mgl*pow2(Msq)*(Mt*MuSUSY*(2*pow2(Dmst12)*pow2(Mst1)*(
        5040*Dmsqst1*pow2(Mgl)*(2*Dmst12 - pow2(Mst2)) - 20*Dmglst1*(40023*
        Dmglst1 + 10097*Mgl)*pow2(Msq)*(2*Dmst12 - pow2(Mst2)) + pow2(Mgl)*
        pow2(Msq)*(37669*Dmst12 + 10296*pow2(Mst2)))*pow2(s2t) + Mt*pow2(Msq)*(
        -16*Dmst12*Mst1*s2t*(2*pow2(Dmst12)*(142987*Dmglst1*Mgl + 37582*pow2(
        Dmglst1) + 20297*pow2(Mgl)) + 5*Dmst12*(-21081*Dmglst1*Mgl + 21081*
        pow2(Dmglst1) - 3457*pow2(Mgl))*pow2(Mst2) - 2*(37582*Dmglst1*Mgl +
        142987*pow2(Dmglst1) + 3012*pow2(Mgl))*pow4(Mst2)) + 3*Mt*(-(pow2(Mgl)*
        (-31963*pow2(Dmst12)*pow2(Mst2) + 131926*pow3(Dmst12) - 68000*Dmst12*
        pow4(Mst2) - 24064*pow6(Mst2))) + 432*pow2(Dmglst1)*(-3185*pow2(Dmst12)
        *pow2(Mst2) + 3185*Dmst12*(pow2(Dmst12) + pow4(Mst2)) + 1566*pow6(Mst2)
        ) + 48*Dmglst1*Mgl*(-8213*pow2(Dmst12)*pow2(Mst2) + 8213*Dmst12*(pow2(
        Dmst12) + pow4(Mst2)) + 4664*pow6(Mst2))))) + 4*Mst1*Tbeta*(-((1260*
        Dmsqst1*pow2(Mgl) - (50485*Dmglst1*Mgl + 200115*pow2(Dmglst1) + 2574*
        pow2(Mgl))*pow2(Msq))*pow3(Dmst12)*pow3(Mst1)*pow3(s2t)) - Dmst12*Mst1*
        s2t*pow2(Mt)*(1260*Dmsqst1*pow2(Mgl)*(22*pow2(Dmst12) - 7*Dmst12*pow2(
        Mst2) - 8*pow4(Mst2)) + pow2(Msq)*(pow2(Mgl)*(32783*pow2(Dmst12) -
        14508*Dmst12*pow2(Mst2) - 15840*pow4(Mst2)) - 24*Dmglst1*(Dmglst1*(
        68291*pow2(Dmst12) - 32498*Dmst12*pow2(Mst2) - 3295*pow4(Mst2)) + Mgl*(
        4744*pow2(Dmst12) - 3821*Dmst12*pow2(Mst2) + 2898*pow4(Mst2))))) +
        pow2(Msq)*(-3*Mt*pow2(Dmst12)*pow2(Mst1)*(Dmst12*(-30241*Dmglst1*Mgl +
        391379*pow2(Dmglst1) - 11261*pow2(Mgl)) - 2*(37582*Dmglst1*Mgl +
        142987*pow2(Dmglst1) + 3012*pow2(Mgl))*pow2(Mst2))*pow2(s2t) + 2*pow3(
        Mt)*(4*Dmglst1*Mgl*(-86833*pow2(Dmst12)*pow2(Mst2) + 287078*pow3(
        Dmst12) - 113412*Dmst12*pow4(Mst2) - 47112*pow6(Mst2)) + pow2(Mgl)*(-
        51181*pow2(Dmst12)*pow2(Mst2) + 152018*pow3(Dmst12) - 49656*Dmst12*
        pow4(Mst2) - 10176*pow6(Mst2)) - 4*pow2(Dmglst1)*(-286982*pow2(Dmst12)*
        pow2(Mst2) + 86737*pow3(Dmst12) + 487227*Dmst12*pow4(Mst2) + 226824*
        pow6(Mst2))))))))))/(1.90512e8*Tbeta*pow2(Mst1)*pow3(Mgl)*pow4(Msq)*
        pow6(Mst2));
}

/**
 * @return returns the susy log^0 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double H3::calc_coef_at_as2_no_sm_logs_log0() const {

   const double result =
         (49*pow3(Mgl)*(-150*Dmsqst1*pow2(Msq)*(-4*Mt*pow2(Dmst12)*pow2(Mst2)*(
        19600*Mst1*Mt*s2t + (-9122 + 14175*z3)*pow2(Mt) + 450*(2 - 21*z3)*pow2(
        Mst1)*pow2(s2t)) + pow3(Dmst12)*(-20800*Mst1*s2t*pow2(Mt) - 75*Mt*(274
        + 63*z3)*pow2(Mst1)*pow2(s2t) + 16*(-6536 + 14175*z3)*pow3(Mt) - 8400*
        pow3(Mst1)*pow3(s2t)) - 200*Dmst12*(-888*Mst1*s2t + Mt*(-158 + 567*z3))
        *pow2(Mt)*pow4(Mst2) - 201600*pow3(Mt)*pow6(Mst2)) - 150*pow2(Dmsqst1)*
        (-(Mt*pow2(Dmst12)*pow2(Mst2)*(36800*Mst1*Mt*s2t + 20112*pow2(Mt) + 75*
        (458 - 945*z3)*pow2(Mst1)*pow2(s2t))) + 3*pow3(Dmst12)*(-20800*Mst1*
        s2t*pow2(Mt) - 525*Mt*(-26 + 45*z3)*pow2(Mst1)*pow2(s2t) + 4*(-3998 +
        14175*z3)*pow3(Mt) - 2800*pow3(Mst1)*pow3(s2t)) - 100*Dmst12*(-1360*
        Mst1*s2t + 63*Mt*(-14 + 27*z3))*pow2(Mt)*pow4(Mst2) - 200*(410 + 567*
        z3)*pow3(Mt)*pow6(Mst2)) - pow4(Msq)*(-20*Mt*pow2(Dmst12)*pow2(Mst2)*(
        300*Mst1*Mt*s2t*(-17512 + 14805*z3) + (-375892 + 621675*z3)*pow2(Mt) -
        900*(-1226 + 495*z3)*pow2(Mst1)*pow2(s2t)) + pow3(Dmst12)*(-30*Mst1*
        s2t*(-1117238 + 877575*z3)*pow2(Mt) - 2250*Mt*(-5570 + 333*z3)*pow2(
        Mst1)*pow2(s2t) + (-41715182 + 38174625*z3)*pow3(Mt) + 3000*(-2722 +
        2259*z3)*pow3(Mst1)*pow3(s2t)) - 6000*Dmst12*(81*Mt*(-406 + 285*z3) +
        8*Mst1*s2t*(-694 + 477*z3))*pow2(Mt)*pow4(Mst2) + 48000*(623 + 963*z3)*
        pow3(Mt)*pow6(Mst2))) + 196*Dmglst1*pow2(Mgl)*(600*pow2(Dmsqst1)*(-6*
        Mt*pow2(Dmst12)*pow2(Mst2)*(-347*Mst1*Mt*s2t + 2538*pow2(Mt) + 1175*
        pow2(Mst1)*pow2(s2t)) + pow3(Dmst12)*(386*Mst1*s2t*pow2(Mt) + 14100*Mt*
        pow2(Mst1)*pow2(s2t) + 15228*pow3(Mt) + 125*pow3(Mst1)*pow3(s2t)) + 2*
        Dmst12*(7614*Mt - 2275*Mst1*s2t)*pow2(Mt)*pow4(Mst2) + 22600*pow3(Mt)*
        pow6(Mst2)) + 600*Dmsqst1*pow2(Msq)*(-2*Mt*pow2(Dmst12)*pow2(Mst2)*(-
        1516*Mst1*Mt*s2t + 7614*pow2(Mt) + 3525*pow2(Mst1)*pow2(s2t)) + pow3(
        Dmst12)*(-564*Mst1*s2t*pow2(Mt) + 14100*Mt*pow2(Mst1)*pow2(s2t) +
        15228*pow3(Mt) + 125*pow3(Mst1)*pow3(s2t)) + 4*Dmst12*(3807*Mt - 1375*
        Mst1*s2t)*pow2(Mt)*pow4(Mst2) + 29200*pow3(Mt)*pow6(Mst2)) + pow4(Msq)*
        (Mt*pow2(Dmst12)*pow2(Mst2)*(40*Mst1*Mt*s2t*(-4511549 + 3729375*z3) + (
        45149198 - 35285625*z3)*pow2(Mt) - 47250*(-430 + 207*z3)*pow2(Mst1)*
        pow2(s2t)) + pow3(Dmst12)*(2*Mst1*s2t*(-28188929 + 23099625*z3)*pow2(
        Mt) + 225*Mt*(-160136 + 100785*z3)*pow2(Mst1)*pow2(s2t) + (115176444 -
        85920750*z3)*pow3(Mt) + 1125*(22174 - 18791*z3)*pow3(Mst1)*pow3(s2t)) +
        40*Dmst12*(150*Mst1*s2t*(-19856 + 17667*z3) + Mt*(-5136871 + 3912300*
        z3))*pow2(Mt)*pow4(Mst2) + 6000*(-31142 + 22653*z3)*pow3(Mt)*pow6(Mst2)
        )) + 2*Mgl*pow2(Dmglst1)*(11760*pow2(Dmsqst1)*(-3*Mt*pow2(Dmst12)*pow2(
        Mst2)*(4470*Mst1*Mt*s2t + 26294*pow2(Mt) + 11975*pow2(Mst1)*pow2(s2t))
        + pow3(Dmst12)*(25750*Mst1*s2t*pow2(Mt) + 71850*Mt*pow2(Mst1)*pow2(s2t)
        + 78882*pow3(Mt) + 625*pow3(Mst1)*pow3(s2t)) + 2*Dmst12*(39441*Mt +
        535*Mst1*s2t)*pow2(Mt)*pow4(Mst2) + 128066*pow3(Mt)*pow6(Mst2)) +
        11760*Dmsqst1*pow2(Msq)*(-(Mt*pow2(Dmst12)*pow2(Mst2)*(8660*Mst1*Mt*s2t
        + 78882*pow2(Mt) + 35925*pow2(Mst1)*pow2(s2t))) + pow3(Dmst12)*(21000*
        Mst1*s2t*pow2(Mt) + 71850*Mt*pow2(Mst1)*pow2(s2t) + 78882*pow3(Mt) +
        625*pow3(Mst1)*pow3(s2t)) + 2*Dmst12*(39441*Mt - 1840*Mst1*s2t)*pow2(
        Mt)*pow4(Mst2) + 160280*pow3(Mt)*pow6(Mst2)) + pow4(Msq)*(Mt*pow2(
        Dmst12)*pow2(Mst2)*(196*Mst1*Mt*s2t*(-353948822 + 292953375*z3) + (
        31554389946 - 25980577875*z3)*pow2(Mt) + 22050*(26498 + 49425*z3)*pow2(
        Mst1)*pow2(s2t)) + 5*pow3(Dmst12)*(-196*Mst1*s2t*(-148185343 +
        123161625*z3)*pow2(Mt) + 4410*Mt*(-612347 + 438045*z3)*pow2(Mst1)*pow2(
        s2t) + (-2800003036 + 2827317150*z3)*pow3(Mt) - 735*(-2573582 +
        2144805*z3)*pow3(Mst1)*pow3(s2t)) + 392*Dmst12*(260*Mst1*s2t*(-579323 +
        490725*z3) + Mt*(-125277461 + 96491250*z3))*pow2(Mt)*pow4(Mst2) + 3920*
        (-12894992 + 9523575*z3)*pow3(Mt)*pow6(Mst2))) + 4*pow3(Dmglst1)*(5880*
        pow2(Dmsqst1)*(-2*Mt*pow2(Dmst12)*pow2(Mst2)*(13291*Mst1*Mt*s2t +
        40812*pow2(Mt) + 18300*pow2(Mst1)*pow2(s2t)) + pow3(Dmst12)*(38922*
        Mst1*s2t*pow2(Mt) + 73200*Mt*pow2(Mst1)*pow2(s2t) + 81624*pow3(Mt) +
        625*pow3(Mst1)*pow3(s2t)) + 2*Dmst12*(40812*Mt + 7121*Mst1*s2t)*pow2(
        Mt)*pow4(Mst2) + 122500*pow3(Mt)*pow6(Mst2)) + 5880*Dmsqst1*pow2(Msq)*(
        -8*Mt*pow2(Dmst12)*pow2(Mst2)*(2729*Mst1*Mt*s2t + 10203*pow2(Mt) +
        4575*pow2(Mst1)*pow2(s2t)) + pow3(Dmst12)*(34172*Mst1*s2t*pow2(Mt) +
        73200*Mt*pow2(Mst1)*pow2(s2t) + 81624*pow3(Mt) + 625*pow3(Mst1)*pow3(
        s2t)) + 12*Dmst12*(6802*Mt + 791*Mst1*s2t)*pow2(Mt)*pow4(Mst2) +
        153928*pow3(Mt)*pow6(Mst2)) + pow4(Msq)*(2*Mt*pow2(Dmst12)*pow2(Mst2)*(
        Mst1*Mt*s2t*(-51549748862 + 42804507375*z3) + (15344105782 -
        12740986125*z3)*pow2(Mt) + 36750*(-109771 + 107379*z3)*pow2(Mst1)*pow2(
        s2t)) + pow3(Dmst12)*(2*Mst1*s2t*(137797425107 - 114549584625*z3)*pow2(
        Mt) - 3675*Mt*(-973342 + 1115325*z3)*pow2(Mst1)*pow2(s2t) + 12*(-
        1749149438 + 1631424375*z3)*pow3(Mt) - 6370*(-2386547 + 1986525*z3)*
        pow3(Mst1)*pow3(s2t)) + 8*Dmst12*(49*Mst1*s2t*(-244084964 + 203974875*
        z3) + Mt*(-5048328734 + 3923356500*z3))*pow2(Mt)*pow4(Mst2) + 2352*(-
        21126629 + 16012125*z3)*pow3(Mt)*pow6(Mst2))))/(1.9845e6*pow3(Mgl)*
        pow3(Mt)*pow4(Msq)*pow6(Mst2));

   return result;
}

/**
 * @return returns the susy log^1 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double H3::calc_coef_at_as2_no_sm_logs_log1() const {

   const double result =
         (4*(2*pow3(Dmglst1)*(2940*pow2(Dmsqst1)*pow2(Mt)*(-((42*Mt + 11*Mst1*s2t)
        *pow2(Dmst12)*pow2(Mst2)) + (42*Mt + 51*Mst1*s2t)*pow3(Dmst12) +
        Dmst12*(42*Mt - 29*Mst1*s2t)*pow4(Mst2) - 65*Mt*pow6(Mst2)) + 5880*
        Dmsqst1*pow2(Msq)*pow2(Mt)*(7*(-3*Mt + Mst1*s2t)*pow2(Dmst12)*pow2(
        Mst2) + (21*Mt + 13*Mst1*s2t)*pow3(Dmst12) + 3*Dmst12*(7*Mt - 9*Mst1*
        s2t)*pow4(Mst2) - 28*Mt*pow6(Mst2)) + pow4(Msq)*(-10*Mt*pow2(Dmst12)*
        pow2(Mst2)*(623998*Mst1*Mt*s2t + 46947*pow2(Mt) + 134505*pow2(Mst1)*
        pow2(s2t)) + pow3(Dmst12)*(17975555*Mst1*s2t*pow2(Mt) + 1956570*Mt*
        pow2(Mst1)*pow2(s2t) + 2032034*pow3(Mt) + 1187760*pow3(Mst1)*pow3(s2t))
        - 2*Dmst12*(546547*Mt + 2011058*Mst1*s2t)*pow2(Mt)*pow4(Mst2) -
        11351536*pow3(Mt)*pow6(Mst2))) - 98*Dmglst1*pow2(Mgl)*(600*Dmsqst1*
        pow2(Msq)*pow2(Mt)*(-((6*Mt + Mst1*s2t)*pow2(Dmst12)*pow2(Mst2)) + (6*
        Mt - 3*Mst1*s2t)*pow3(Dmst12) + Dmst12*(6*Mt + 5*Mst1*s2t)*pow4(Mst2) +
        25*Mt*pow6(Mst2)) + 300*pow2(Dmsqst1)*pow2(Mt)*(3*(-4*Mt + Mst1*s2t)*
        pow2(Dmst12)*pow2(Mst2) + (12*Mt - 11*Mst1*s2t)*pow3(Dmst12) + Dmst12*(
        12*Mt + 5*Mst1*s2t)*pow4(Mst2) + 35*Mt*pow6(Mst2)) + pow4(Msq)*(Mt*
        pow2(Dmst12)*pow2(Mst2)*(29758*Mst1*Mt*s2t + 6677*pow2(Mt) + 22350*
        pow2(Mst1)*pow2(s2t)) + pow3(Dmst12)*(-34587*Mst1*s2t*pow2(Mt) - 18045*
        Mt*pow2(Mst1)*pow2(s2t) + 22414*pow3(Mt) + 3325*pow3(Mst1)*pow3(s2t)) -
        8*Dmst12*(4471*Mt + 6875*Mst1*s2t)*pow2(Mt)*pow4(Mst2) + 50800*pow3(Mt)
        *pow6(Mst2))) - Mgl*pow2(Dmglst1)*(5880*pow2(Dmsqst1)*pow2(Mt)*(-3*(3*
        Mt + 5*Mst1*s2t)*pow2(Dmst12)*pow2(Mst2) + (9*Mt - 25*Mst1*s2t)*pow3(
        Dmst12) + Dmst12*(9*Mt + 55*Mst1*s2t)*pow4(Mst2) + 112*Mt*pow6(Mst2)) +
        5880*Dmsqst1*pow2(Msq)*pow2(Mt)*(-((9*Mt + 40*Mst1*s2t)*pow2(Dmst12)*
        pow2(Mst2)) + 9*Mt*pow3(Dmst12) + Dmst12*(9*Mt + 80*Mst1*s2t)*pow4(
        Mst2) + 145*Mt*pow6(Mst2)) + pow4(Msq)*(3*Mt*pow2(Dmst12)*pow2(Mst2)*(
        2334948*Mst1*Mt*s2t - 104761*pow2(Mt) + 112210*pow2(Mst1)*pow2(s2t)) +
        pow3(Dmst12)*(-10892014*Mst1*s2t*pow2(Mt) + 1366365*Mt*pow2(Mst1)*pow2(
        s2t) + 177178*pow3(Mt) - 363580*pow3(Mst1)*pow3(s2t)) + 196*Dmst12*(
        2303*Mt - 30942*Mst1*s2t)*pow2(Mt)*pow4(Mst2) + 13177472*pow3(Mt)*pow6(
        Mst2))) - 49*pow3(Mgl)*(-150*Dmsqst1*Mt*pow2(Msq)*(pow2(Dmst12)*pow2(
        Mst2)*(-80*Mst1*Mt*s2t - 17*pow2(Mt) + 180*pow2(Mst1)*pow2(s2t)) + 2*(
        20*Mst1*Mt*s2t + 187*pow2(Mt) - 90*pow2(Mst1)*pow2(s2t))*pow3(Dmst12) -
        20*Dmst12*Mt*(17*Mt - 6*Mst1*s2t)*pow4(Mst2) - 1080*pow2(Mt)*pow6(Mst2)
        ) - 150*Mt*pow2(Dmsqst1)*(-4*pow2(Dmst12)*pow2(Mst2)*(10*Mst1*Mt*s2t +
        3*pow2(Mt) - 45*pow2(Mst1)*pow2(s2t)) + 9*(41*pow2(Mt) - 20*pow2(Mst1)*
        pow2(s2t))*pow3(Dmst12) + 5*Dmst12*Mt*(-69*Mt + 16*Mst1*s2t)*pow4(Mst2)
        - 1010*pow2(Mt)*pow6(Mst2)) - pow4(Msq)*(-2*Mt*pow2(Dmst12)*pow2(Mst2)*
        (25850*Mst1*Mt*s2t + 21033*pow2(Mt) + 75*pow2(Mst1)*pow2(s2t)) + pow3(
        Dmst12)*(68264*Mst1*s2t*pow2(Mt) + 5700*Mt*pow2(Mst1)*pow2(s2t) +
        36107*pow3(Mt) + 250*pow3(Mst1)*pow3(s2t)) + 200*Dmst12*(131*Mt + 100*
        Mst1*s2t)*pow2(Mt)*pow4(Mst2) + 400*(-533 + 54*z3)*pow3(Mt)*pow6(Mst2))
        )))/(33075.*pow3(Mgl)*pow3(Mt)*pow4(Msq)*pow6(Mst2));

   return result;
}

/**
 * @return returns the susy log^2 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double H3::calc_coef_at_as2_no_sm_logs_log2() const {

   const double result =
         (8*(7*Dmglst1*pow2(Mgl)*pow4(Msq)*(-2*Mt*pow2(Dmst12)*pow2(Mst2)*(-1304*
        Mst1*Mt*s2t + 888*pow2(Mt) + 645*pow2(Mst1)*pow2(s2t)) + pow3(Dmst12)*(
        -1544*Mst1*s2t*pow2(Mt) + 2250*Mt*pow2(Mst1)*pow2(s2t) + 1640*pow3(Mt)
        - 275*pow3(Mst1)*pow3(s2t)) + 8*Dmst12*(239*Mt - 35*Mst1*s2t)*pow2(Mt)*
        pow4(Mst2) + 6520*pow3(Mt)*pow6(Mst2)) + pow3(Dmglst1)*pow4(Msq)*(-4*
        Mt*pow2(Dmst12)*pow2(Mst2)*(-6158*Mst1*Mt*s2t + 7818*pow2(Mt) - 5565*
        pow2(Mst1)*pow2(s2t)) + pow3(Dmst12)*(-109632*Mst1*s2t*pow2(Mt) - 8820*
        Mt*pow2(Mst1)*pow2(s2t) + 98096*pow3(Mt) - 9765*pow3(Mst1)*pow3(s2t)) +
        16*Dmst12*(-2222*Mt + 5257*Mst1*s2t)*pow2(Mt)*pow4(Mst2) + 36232*pow3(
        Mt)*pow6(Mst2)) + Mgl*pow2(Dmglst1)*pow4(Msq)*(-3*Mt*pow2(Dmst12)*pow2(
        Mst2)*(-7448*Mst1*Mt*s2t + 11772*pow2(Mt) + 35*pow2(Mst1)*pow2(s2t)) +
        pow3(Dmst12)*(-53536*Mst1*s2t*pow2(Mt) + 16905*Mt*pow2(Mst1)*pow2(s2t)
        + 68252*pow3(Mt) - 5285*pow3(Mst1)*pow3(s2t)) + 28*Dmst12*(85*Mt +
        1164*Mst1*s2t)*pow2(Mt)*pow4(Mst2) + 49588*pow3(Mt)*pow6(Mst2)) + 7*
        pow3(Mgl)*(Mt*pow2(Dmst12)*pow2(Mst2)*(1760*Mst1*Mt*s2t + 538*pow2(Mt)
        + 105*pow2(Mst1)*pow2(s2t))*pow4(Msq) - pow3(Dmst12)*(1032*Mst1*s2t*
        pow2(Mt) + 480*Mt*pow2(Mst1)*pow2(s2t) + 644*pow3(Mt) - 45*pow3(Mst1)*
        pow3(s2t))*pow4(Msq) - 40*Dmst12*(11*Mt + 61*Mst1*s2t)*pow2(Mt)*pow4(
        Msq)*pow4(Mst2) + 10*pow3(Mt)*(45*pow2(Dmsqst1) + 90*Dmsqst1*pow2(Msq)
        + 442*pow4(Msq))*pow6(Mst2))))/(315.*pow3(Mgl)*pow3(Mt)*pow4(Msq)*pow6(
        Mst2));

   return result;
}

/**
 * @return returns the susy log^3 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double H3::calc_coef_at_as2_no_sm_logs_log3() const {

   const double result = -298.6666666666667;

   return result;
}

} // namespace hierarchies
} // namespace himalaya
