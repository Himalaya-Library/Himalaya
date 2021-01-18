// ====================================================================
// This file is part of Himalaya.
//
// Himalaya is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================

#include "H3q22g.hpp"
#include "Constants.hpp"
#include "HierarchyFlags.hpp"
#include "Powers.hpp"
#include <cmath>

namespace himalaya {
namespace hierarchies {

/**
 * Constructor
 * @param expansionDepth the flagMap for the truncation of expansion variables
 * @param Al4p a double alpha_s/4/Pi
 * @param beta a double which is the mixing angle beta
 * @param Dmglst1 a double Mgl - Mst1
 * @param Dmst12 a double Mst1^2 - Mst2^2
 * @param Dmsqst1 a double Msq^2 - Mst1^2
 * @param lmMt a double log((<renormalization scale> / Mt)^2)
 * @param lmMst1 a double log((<renormalization scale> / Mst1)^2)
 * @param Mt a double top/bottom quark mass
 * @param Mst1 a double stop 1 mass
 * @param Mst2 a double stop 2 mass
 * @param Msq a double average squark mass w/o the stop quarks
 * @param MuSUSY a double mu parameter
 * @param s2t a double 2 times the sine of the stop/sbottom quark mixing angle
 * @param mdrFlag an int 0 for DR and 1 for MDR scheme
 * @param oneLoopFlag an int flag to consider the one-loop expansion terms
 * @param twoLoopFlag an int flag to consider the two-loop expansion terms
 * @param threeLoopFlag an int flag to consider the three-loop expansion terms
 */
H3q22g::H3q22g(const ExpansionFlags_t& expansionDepth, double Al4p, double beta,
                 double Dmglst1, double Dmst12, double Dmsqst1, double lmMt, double lmMst1,
                 double Mt, double Mst1, double Mst2, double Msq, double MuSUSY,
                 double s2t,
                 int mdrFlag, int oneLoopFlag, int twoLoopFlag, int threeLoopFlag){
   // abbrev for tan(beta) and sin(beta)
   Tbeta = tan(beta);
   Sbeta = sin(beta);
   this -> Dmglst1 = Dmglst1;
   this -> Dmst12 = Dmst12;
   this -> Dmsqst1 = Dmsqst1;
   this -> lmMst1 = lmMst1;
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
 * @return The diagonal (1, 1) matrix element of the Higgs mass matrix as a double for the hierarchy 'H3q22g'
 */
double H3q22g::getS1() const {
   return -(pow2(Mt)*pow2(MuSUSY)*(-4*Mst1*xDmst12*pow3(Dmst12)*(12000*Dmsqst1*s2t*
        threeLoopFlag*pow2(Al4p)*(20*Dmglst1*Mt*(Dmsqst1*xDmsqst1 + pow2(Msq))
        + 3*Mst1*(Dmsqst1*(28*Mt + Mst1*s2t*(7 + 6*shiftst1 - 24*lmMst1*(-1 +
        shiftst2) + 30*shiftst2))*xDmsqst1 + (28*Mt + Mst1*s2t*(7 - 24*lmMst1*(
        -1 + shiftst2) + 36*shiftst2))*pow2(Msq)))*pow2(Mst1) + pow4(Msq)*(
        Al4p*(8*Dmglst1*Mt*(4*Dmglst1*(720*Dmglst1*s2t*twoLoopFlag*xDmglst1 +
        Al4p*threeLoopFlag*(Dmglst1*s2t*xDmglst1*(14217821 + 161940*lmMst1 -
        28800*pow2(lmMst1)) + Mt*(6233611 + 58800*lmMst1 - 4560*lmMt - 14400*
        pow2(lmMst1)))) + 5*Mst1*(-2160*Dmglst1*s2t*twoLoopFlag + Al4p*
        threeLoopFlag*(6*Dmglst1*s2t*(954181 + 11256*lmMst1 - 11520*pow2(
        lmMst1)) + 5*Mt*(403559 + 384*(lmMst1 + lmMt) - 4608*pow2(lmMst1))))) -
        15*pow2(Mst1)*(960*Dmglst1*s2t*((-5 + 6*lmMst1)*Mt + Dmglst1*(-11 + 6*
        lmMst1)*s2t)*twoLoopFlag + Al4p*threeLoopFlag*(16*Dmglst1*s2t*(Mt*(
        1282471 - 7264*lmMst1 - 18120*pow2(lmMst1)) + Dmglst1*s2t*(655743 +
        5288*lmMst1 - 11820*pow2(lmMst1))) + (1763661 - 47104*lmMst1 + 5120*
        lmMt - 24576*pow2(lmMst1))*pow2(Mt))) - 600*s2t*(24*((5 + 6*lmMst1)*Mt
        + 4*Dmglst1*(-1 + 3*lmMst1)*s2t)*twoLoopFlag + Al4p*threeLoopFlag*(
        Dmglst1*s2t*(84209 + 1264*lmMst1 - 240*pow2(lmMst1)) - 2*Mt*(-36863 +
        80*lmMst1 + 552*pow2(lmMst1))))*pow3(Mst1)) + 15*(270*oneLoopFlag +
        Al4p*(240*(twoLoopFlag + 6*lmMst1*twoLoopFlag) + Al4p*threeLoopFlag*(
        350605 + 4320*shiftst1 + 2880*shiftst2 + 8352*shiftst3 - 96*lmMst1*(-
        115 + 90*shiftst1 + 60*shiftst2 + 54*shiftst3) - 2160*pow2(lmMst1))))*
        pow2(s2t)*pow4(Mst1))) - threeLoopFlag*pow2(Al4p)*(480*pow2(s2t)*pow3(
        Mst1)*((Dmst12*pow2(Msq)*pow2(Mst2)*(2*Dmst12*pow2(Dmglst1)*(655743 +
        5288*lmMst1 - 11820*pow2(lmMst1))*pow2(Msq) + 5*Dmglst1*Dmst12*Mst1*(
        84209 + 1264*lmMst1 - 240*pow2(lmMst1))*pow2(Msq) + 10*pow2(Mst1)*(30*
        Dmsqst1*(Dmst12*(-7 - 18*shiftst1 - 18*shiftst2 + 12*lmMst1*(-2 +
        shiftst1 + shiftst2)) - 12*(-3 + 2*lmMst1)*(shiftst1 - shiftst2)*pow2(
        Mst2)) + pow2(Msq)*(Dmst12*(1429 - 180*shiftst1 - 180*shiftst2 +
        lmMst1*(-454 + 360*shiftst1 + 360*shiftst2) + 24*pow2(lmMst1)) - 360*(-
        1 + 2*lmMst1)*(shiftst1 - shiftst2)*pow2(Mst2)))))/2. + pow2(Mst1)*(90*
        Dmst12*shiftst3*pow2(Mst2)*(Dmst12*(-7 + 6*lmMst1) + 2*(1 - 2*lmMst1)*
        pow2(Mst2))*pow4(Msq) + 36*z2*(-2*xDmst12*pow3(Dmst12)*(100*Dmsqst1*
        shiftst2*(Dmsqst1*xDmsqst1 + pow2(Msq)) + (15*shiftst1 + 10*shiftst2 +
        9*shiftst3)*pow4(Msq)) + Dmst12*pow2(Mst2)*(5*shiftst3*(3*Dmst12 - 2*
        pow2(Mst2))*pow4(Msq) - 50*(-(Dmst12*(shiftst1 + shiftst2)) + 2*(
        shiftst1 - shiftst2)*pow2(Mst2))*(xDmsqst1*pow2(Dmsqst1) + Dmsqst1*
        pow2(Msq) + pow4(Msq)))))) + 2*Mst1*pow2(Mst2)*(240*Dmst12*Mst1*s2t*(-
        50*Mst1*xDmsqst1*pow2(Dmsqst1)*(40*Dmglst1*Mt*(Dmst12 - pow2(Mst2)) +
        3*Mst1*(Dmst12*(56*Mt + Mst1*s2t*(7 + 30*shiftst1 + 6*shiftst2 - 12*
        lmMst1*(-2 + shiftst1 + shiftst2))) + 8*(-7*Mt + 3*(-2 + lmMst1)*Mst1*
        s2t*(shiftst1 - shiftst2))*pow2(Mst2))) - Mt*pow2(Msq)*(2*pow2(Dmglst1)
        *(954181 + 11256*lmMst1 - 11520*pow2(lmMst1))*pow2(Msq)*(Dmst12 - pow2(
        Mst2)) + 5*pow2(Mst1)*(1680*Dmsqst1*(Dmst12 - pow2(Mst2)) + pow2(Msq)*(
        3*Dmst12*(-10473 + 40*lmMst1 + 256*pow2(lmMst1)) - 8*(1361 + 10*lmMst1
        + 54*pow2(lmMst1))*pow2(Mst2))) + Dmglst1*Mst1*(2000*Dmsqst1*(Dmst12 -
        pow2(Mst2)) + pow2(Msq)*(Dmst12*(-949861 + 1944*lmMst1 + 11520*pow2(
        lmMst1)) + 20*(-33261 + 532*lmMst1 + 660*pow2(lmMst1))*pow2(Mst2))))) +
        pow2(Mt)*pow4(Msq)*(-400*Dmglst1*Mst1*(Dmst12*(403559 + 384*(lmMst1 +
        lmMt) - 4608*pow2(lmMst1))*(Dmst12 - pow2(Mst2)) - 24*(9631 + 16*lmMst1
        + 48*lmMt - 192*pow2(lmMst1))*pow4(Mst2)) + 15*pow2(Mst1)*(pow2(Dmst12)
        *(852541 + 9216*lmMst1 - 10240*lmMt + 6144*pow2(lmMst1)) + 160*Dmst12*(
        11389 - 704*lmMst1 + 192*lmMt - 384*pow2(lmMst1))*pow2(Mst2) + 1920*(
        349 - 56*lmMst1 + 24*lmMt - 32*pow2(lmMst1))*pow4(Mst2)) - 32*pow2(
        Dmglst1)*(-2*Dmst12*(-6233611 - 58800*lmMst1 + 4560*lmMt + 14400*pow2(
        lmMst1))*(Dmst12 - pow2(Mst2)) + 25*(-230009 - 4704*lmMst1 + 480*lmMt +
        1152*pow2(lmMst1))*pow4(Mst2)))) + 225*z3*(-2*Mst1*xDmst12*pow3(Dmst12)
        *(Mt*(-8*Dmglst1*Mst1*(-73917*Mt + 285974*Mst1*s2t) + 48*(30678*Mt +
        35135*Mst1*s2t)*pow2(Dmglst1) - (197889*Mt + 324752*Mst1*s2t)*pow2(
        Mst1) + 3371456*s2t*xDmglst1*pow3(Dmglst1))*pow4(Msq) - pow2(Mst1)*
        pow2(s2t)*(-10080*Dmsqst1*(Dmsqst1*xDmsqst1 + pow2(Msq))*pow2(Mst1) + (
        403880*Dmglst1*Mst1 + 1197040*pow2(Dmglst1) - 37669*pow2(Mst1))*pow4(
        Msq))) + pow2(Mst2)*(-8*pow2(Dmst12)*pow2(s2t)*pow3(Mst1)*(-1260*
        Dmsqst1*(Dmsqst1*xDmsqst1 + pow2(Msq))*pow2(Mst1) + (50485*Dmglst1*Mst1
        + 149630*pow2(Dmglst1) + 2574*pow2(Mst1))*pow4(Msq)) + Mt*pow4(Msq)*(
        16*Dmst12*s2t*pow2(Mst1)*(210810*pow2(Dmglst1)*(Dmst12 - pow2(Mst2)) -
        pow2(Mst1)*(17285*Dmst12 + 6024*pow2(Mst2)) - Dmglst1*Mst1*(105405*
        Dmst12 + 75164*pow2(Mst2))) - 1408*xDmglst1*pow3(Dmglst1)*(-4789*
        Dmst12*Mst1*s2t*(Dmst12 - pow2(Mst2)) + 1503*Mt*pow4(Mst2)) + 3*Mst1*
        Mt*(96*pow2(Dmglst1)*(10226*pow2(Dmst12) - 10226*Dmst12*pow2(Mst2) -
        4715*pow4(Mst2)) + 48*Dmglst1*Mst1*(8213*pow2(Dmst12) - 8213*Dmst12*
        pow2(Mst2) - 4664*pow4(Mst2)) - pow2(Mst1)*(31963*pow2(Dmst12) + 68000*
        Dmst12*pow2(Mst2) + 24064*pow4(Mst2))))))) + 8*pow2(Mst2)*pow4(Msq)*(
        1800*Al4p*Dmst12*s2t*twoLoopFlag*pow2(Mst1)*(Dmst12*Mst1*s2t*(8*
        Dmglst1*(1 - 3*lmMst1)*Mst1 + (22 - 12*lmMst1)*pow2(Dmglst1) + 3*pow2(
        Mst1)) + 4*Mt*(-6*pow2(Dmglst1)*(Dmst12 - pow2(Mst2)) + Dmglst1*(5 - 6*
        lmMst1)*Mst1*pow2(Mst2) - pow2(Mst1)*(Dmst12 + 3*(1 + 2*lmMst1)*pow2(
        Mst2)))) - 16*Al4p*Mt*xDmglst1*pow3(Dmglst1)*(-720*Dmst12*Mst1*s2t*
        twoLoopFlag*(Dmst12 - pow2(Mst2)) + Al4p*threeLoopFlag*(Dmst12*Mst1*
        s2t*(-14217821 - 161940*lmMst1 + 28800*pow2(lmMst1))*(Dmst12 - pow2(
        Mst2)) + 2*(2219399 + 9600*lmMst1 + 960*lmMt)*Mt*pow4(Mst2))) + 2025*
        oneLoopFlag*pow2(Dmst12)*pow2(s2t)*pow5(Mst1))))/(777600.*pow4(Msq)*
        pow5(Mst1)*pow6(Mst2));
}

/**
 * @return The diagonal (2, 2) matrix element of the Higgs mass matrix as a double for the hierarchy 'H3q22g'
 */
double H3q22g::getS2() const {
   return (pow2(Mst2)*(Mst1*threeLoopFlag*pow2(Al4p)*(23520*Mt*xDmsqst1*pow2(
        Dmsqst1)*(64*xDmglst1*pow2(Mt)*pow2(Sbeta)*pow3(Dmglst1)*(Dmst12*(1331
        - 420*lmMst1 + 420*lmMt)*s2t*Tbeta*pow2(Mst1)*(Dmst12 - pow2(Mst2)) -
        Mt*((1541 - 420*lmMst1 + 420*lmMt)*MuSUSY + (2579 + 120*lmMst1 - 120*
        lmMt)*Mst1*Tbeta)*pow4(Mst2)) - 5*pow3(Mst1)*(50*Dmst12*Mt*pow2(Mst2)*(
        2*(8*(65 - 6*lmMst1 + 6*lmMt)*MuSUSY + 9*Mst1*(49 + 2*lmMt + lmMst1*(46
        - 48*shiftst2) + 96*shiftst2)*Tbeta)*pow2(Mt)*pow2(Sbeta) - Mt*s2t*(-
        168*Tbeta*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + Mst1*(MuSUSY*(687 - 864*
        lmMst1*(-1 + shiftst2) + 1728*shiftst2) + 16*(-85 + 12*lmMst1 - 12*
        lmMt)*Mst1*Tbeta)*pow2(Sbeta)) - 72*(-2 + lmMst1)*Mst1*(shiftst1 -
        shiftst2)*Tbeta*pow2(s2t)*(pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 12*pow2(
        Mst1)*pow2(Sbeta))) - pow2(Dmst12)*(-200*s2t*pow2(Mt)*(-42*Tbeta*pow2(
        MuSUSY)*(-1 + pow2(Sbeta)) + Mst1*(3*MuSUSY*(11 + 36*shiftst2) + 4*(-23
        + 6*lmMst1 - 6*lmMt)*Mst1*Tbeta)*pow2(Sbeta)) + 75*Mst1*Mt*pow2(s2t)*(-
        2*(7 + 30*shiftst1 + 6*shiftst2 - 12*lmMst1*(-2 + shiftst1 + shiftst2))
        *Tbeta*pow2(MuSUSY)*(1 - pow2(Sbeta)) + Mst1*(168*MuSUSY + Mst1*(229 +
        864*shiftst1 - 288*shiftst2 + 144*lmMst1*(2 - 3*shiftst1 + shiftst2))*
        Tbeta)*pow2(Sbeta)) + pow2(Sbeta)*(8*(100*(65 - 6*lmMst1 + 6*lmMt)*
        MuSUSY + 3*Mst1*(419 - 60*lmMst1 + 60*lmMt + 900*shiftst2)*Tbeta)*pow3(
        Mt) - 5400*(-2 + lmMst1)*MuSUSY*(shiftst1 - shiftst2)*pow3(Mst1)*pow3(
        s2t))) + 200*(MuSUSY*((364 - 48*lmMst1 + 48*lmMt)*Mt - 216*(-2 +
        lmMst1)*Mst1*s2t*(shiftst1 - shiftst2)) - Mst1*Mt*Tbeta*(205 + 174*lmMt
        - 432*(shiftst1 + shiftst2) - 6*lmMst1*(101 + 18*lmMt - 36*(shiftst1 +
        shiftst2)) + 54*(pow2(lmMst1) + pow2(lmMt))))*pow2(Mt)*pow2(Sbeta)*
        pow4(Mst2)) + 8*Dmglst1*Mst1*Mt*(3*Dmglst1*pow2(Sbeta)*(Dmst12*Mst1*(2*
        Mt*(-225*MuSUSY*s2t + (457 + 510*lmMst1 - 510*lmMt)*Mt*Tbeta + 10*(397
        - 30*lmMst1 + 30*lmMt)*Mst1*s2t*Tbeta)*pow2(Mst2) - Dmst12*(-10*Mt*s2t*
        (45*MuSUSY + 2*(-397 + 30*lmMst1 - 30*lmMt)*Mst1*Tbeta) + 2*(457 + 510*
        lmMst1 - 510*lmMt)*Tbeta*pow2(Mt) + 225*Tbeta*pow2(Mst1)*pow2(s2t))) +
        2*(20*(206 - 15*lmMst1 + 15*lmMt)*MuSUSY + 9*(279 + 70*lmMst1 - 70*
        lmMt)*Mst1*Tbeta)*pow2(Mt)*pow4(Mst2)) - 5*Mst1*(2*Dmst12*Mt*pow2(Mst2)
        *(2*Mt*((557 + 120*lmMst1 - 120*lmMt)*MuSUSY + 9*(-423 + 20*lmMst1 -
        20*lmMt)*Mst1*Tbeta)*pow2(Sbeta) + 25*s2t*(5*Tbeta*pow2(MuSUSY)*(-1 +
        pow2(Sbeta)) + Mst1*(282*MuSUSY + (91 + 6*lmMst1 - 6*lmMt)*Mst1*Tbeta)*
        pow2(Sbeta))) - pow2(Dmst12)*((4*((557 + 120*lmMst1 - 120*lmMt)*MuSUSY
        + 9*(-423 + 20*lmMst1 - 20*lmMt)*Mst1*Tbeta)*pow2(Mt) + 75*(5*MuSUSY -
        94*Mst1*Tbeta)*pow2(Mst1)*pow2(s2t))*pow2(Sbeta) + 2*Mt*s2t*(125*Tbeta*
        pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 3*Mst1*(2350*MuSUSY + (347 - 30*
        lmMst1 + 30*lmMt)*Mst1*Tbeta)*pow2(Sbeta))) + 100*((44 + 3*lmMst1 - 3*
        lmMt)*MuSUSY + (-226 + 21*lmMst1 - 21*lmMt)*Mst1*Tbeta)*pow2(Mt)*pow2(
        Sbeta)*pow4(Mst2)))) + 2*Mt*pow2(Msq)*(784*Mst1*MuSUSY*pow2(Sbeta)*(
        450*Dmst12*Mst1*s2t*(5*pow2(Mst1)*(5*Dmsqst1*(161*Dmst12 + 24*(1 + 12*
        lmMst1)*pow2(Mst2)) + pow2(Msq)*(Dmst12*(2071 + 296*lmMst1 - 96*lmMt -
        600*pow2(lmMst1)) + 8*(631 - lmMst1 + 36*lmMt + 21*pow2(lmMst1))*pow2(
        Mst2))) + 4*Dmglst1*Mst1*(4700*Dmsqst1*(Dmst12 - pow2(Mst2)) - pow2(
        Msq)*(Dmst12*(17539 + 574*lmMst1 + 160*lmMt - 1920*pow2(lmMst1)) + 5*(-
        4539 + 596*lmMst1 - 48*lmMt + 516*pow2(lmMst1))*pow2(Mst2))) + pow2(
        Dmglst1)*(360*Dmsqst1*(Dmst12 - pow2(Mst2)) + pow2(Msq)*(Dmst12*(-
        515917 - 6972*lmMst1 + 192*lmMt + 11520*pow2(lmMst1)) + 2*(-32101 +
        5044*lmMst1 - 400*lmMt + 5100*pow2(lmMst1))*pow2(Mst2))))*pow2(Mt) +
        225*Mt*pow2(Dmst12)*pow2(Mst1)*(pow2(Dmglst1)*(954181 + 11256*lmMst1 -
        11520*pow2(lmMst1))*pow2(Msq) + 10*Dmglst1*Mst1*(100*Dmsqst1 + (33261 -
        532*lmMst1 - 660*pow2(lmMst1))*pow2(Msq)) + 20*(210*Dmsqst1 + (1361 +
        10*lmMst1 + 54*pow2(lmMst1))*pow2(Msq))*pow2(Mst1))*pow2(s2t) + 20250*
        s2t*pow3(Mst1)*(shiftst3*pow2(Msq)*(-8*Dmst12*(-3 + 2*lmMst1)*pow2(
        Mst2)*pow2(Mt) + pow2(Dmst12)*(16*(-2 + lmMst1)*pow2(Mt) + (1 - 2*
        lmMst1)*pow2(Mst1)*pow2(s2t)) + 8*(-1 + 2*lmMst1)*pow2(Mt)*pow4(Mst2))
        + 10*(Dmsqst1*(-3 + 2*lmMst1) + (-1 + 2*lmMst1)*pow2(Msq))*(-8*Dmst12*
        shiftst2*pow2(Mst2)*pow2(Mt) + (-shiftst1 + shiftst2)*pow2(Dmst12)*
        pow2(Mst1)*pow2(s2t) + 8*(shiftst1 - shiftst2)*pow2(Mt)*pow4(Mst2))) +
        pow3(Mt)*(4*pow2(Dmglst1)*(3*Dmst12*(84334067 + 120*lmMst1*(2843 - 120*
        lmMt) + 202200*lmMt - 828000*pow2(lmMst1) - 21600*pow2(lmMt))*pow2(Msq)
        *(Dmst12 - pow2(Mst2)) + 40*(90*Dmsqst1*(206 - 15*lmMst1 + 15*lmMt) + (
        -3030652 - 19305*lmMt + 6*lmMst1*(1183 + 645*lmMt) + 55530*pow2(lmMst1)
        + 5400*pow2(lmMt))*pow2(Msq))*pow4(Mst2)) + 2*Dmglst1*Mst1*(1200*
        Dmsqst1*(Dmst12*(557 + 120*lmMst1 - 120*lmMt)*(Dmst12 - pow2(Mst2)) -
        50*(26 + 3*lmMst1 - 3*lmMt)*pow4(Mst2)) - pow2(Msq)*(-3*pow2(Dmst12)*(-
        39474953 + lmMst1*(52780 - 87600*lmMt) - 3300*lmMt + 319200*pow2(
        lmMst1) + 14400*pow2(lmMt)) + 40*Dmst12*(3746977 + 4005*lmMt - 18*
        lmMst1*(2711 + 1215*lmMt) - 52380*pow2(lmMst1))*pow2(Mst2) + 6000*(9961
        + 66*lmMt - 10*lmMst1*(67 + 24*lmMt) + 42*pow2(lmMst1) + 72*pow2(lmMt))
        *pow4(Mst2))) + 15*pow2(Mst1)*(4000*Dmsqst1*(Dmst12*(65 - 6*lmMst1 + 6*
        lmMt)*(Dmst12 - pow2(Mst2)) + 6*(-20 + 3*lmMst1 - 3*lmMt)*pow4(Mst2)) -
        pow2(Msq)*(-(pow2(Dmst12)*(-2284899 + 49840*lmMt - 32*lmMst1*(-1793 +
        555*lmMt) + 87360*pow2(lmMst1) + 28800*pow2(lmMt))) + 200*Dmst12*(11697
        + 448*lmMst1 + 330*lmMt - 372*lmMst1*lmMt + 408*pow2(lmMst1) + 288*
        pow2(lmMt))*pow2(Mst2) + 1600*(434 - 83*lmMst1 + 174*lmMt - 66*lmMst1*
        lmMt + 183*pow2(lmMst1) + 108*pow2(lmMt))*pow4(Mst2))))) + Mt*Tbeta*(
        245*(pow2(MuSUSY)*(16*pow2(Dmglst1)*pow2(Msq)*(2*Dmst12*Mt*(15*Mst1*
        s2t*(954181 + 11256*lmMst1 - 11520*pow2(lmMst1)) - 2*Mt*(-6233611 -
        58800*lmMst1 + 4560*lmMt + 14400*pow2(lmMst1)))*pow2(Mst2) + pow2(
        Dmst12)*(30*Mst1*Mt*s2t*(-954181 - 11256*lmMst1 + 11520*pow2(lmMst1)) +
        4*(-6233611 - 58800*lmMst1 + 4560*lmMt + 14400*pow2(lmMst1))*pow2(Mt) +
        15*(655743 + 5288*lmMst1 - 11820*pow2(lmMst1))*pow2(Mst1)*pow2(s2t)) -
        50*(-230009 - 4704*lmMst1 + 480*lmMt + 1152*pow2(lmMst1))*pow2(Mt)*
        pow4(Mst2)) + 40*Dmglst1*Mst1*(10*Dmst12*Mt*(1200*Dmsqst1*Mst1*s2t + (
        Mt*(403559 + 384*(lmMst1 + lmMt) - 4608*pow2(lmMst1)) - 12*Mst1*s2t*(-
        33261 + 532*lmMst1 + 660*pow2(lmMst1)))*pow2(Msq))*pow2(Mst2) - pow2(
        Dmst12)*(12000*Dmsqst1*Mst1*Mt*s2t + pow2(Msq)*(6*Mst1*Mt*s2t*(-949861
        + 1944*lmMst1 + 11520*pow2(lmMst1)) + 10*(403559 + 384*(lmMst1 + lmMt)
        - 4608*pow2(lmMst1))*pow2(Mt) + 15*(-84209 - 1264*lmMst1 + 240*pow2(
        lmMst1))*pow2(Mst1)*pow2(s2t))) + 240*(9631 + 16*lmMst1 + 48*lmMt -
        192*pow2(lmMst1))*pow2(Msq)*pow2(Mt)*pow4(Mst2)) + 15*pow2(Mst1)*(160*
        Dmst12*Mt*(840*Dmsqst1*Mst1*s2t + (Mt*(11389 - 704*lmMst1 + 192*lmMt -
        384*pow2(lmMst1)) + 4*Mst1*s2t*(1361 + 10*lmMst1 + 54*pow2(lmMst1)))*
        pow2(Msq))*pow2(Mst2) + pow2(Dmst12)*(-2400*Dmsqst1*Mst1*s2t*(56*Mt + (
        7 + 24*lmMst1)*Mst1*s2t) + pow2(Msq)*(-240*Mst1*Mt*s2t*(-10473 + 40*
        lmMst1 + 256*pow2(lmMst1)) + (852541 + 9216*lmMst1 - 10240*lmMt + 6144*
        pow2(lmMst1))*pow2(Mt) + 80*(1429 - 454*lmMst1 + 24*pow2(lmMst1))*pow2(
        Mst1)*pow2(s2t))) + 1920*(349 - 56*lmMst1 + 24*lmMt - 32*pow2(lmMst1))*
        pow2(Msq)*pow2(Mt)*pow4(Mst2))) + 21600*pow4(Mst1)*(-(shiftst3*pow2(
        Msq)*(-2*Dmst12*pow2(Mst2)*((-1 + 2*lmMst1)*pow2(MuSUSY)*pow2(s2t)*(-1
        + pow2(Sbeta)) + 12*((3 - 2*lmMst1)*pow2(Mt) + (-1 + 2*lmMst1)*pow2(
        Mst1)*pow2(s2t))*pow2(Sbeta)) + pow2(Dmst12)*((-7 + 6*lmMst1)*pow2(
        MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta)) + 6*(-8*(-2 + lmMst1)*pow2(Mt) + (
        -15 + 14*lmMst1)*pow2(Mst1)*pow2(s2t))*pow2(Sbeta)) + 24*(1 - 2*lmMst1)
        *pow2(Mt)*pow2(Sbeta)*pow4(Mst2))) - 10*(Dmsqst1*(3 - 2*lmMst1) + (1 -
        2*lmMst1)*pow2(Msq))*(pow2(Dmst12)*pow2(s2t)*(-((shiftst1 + shiftst2)*
        pow2(MuSUSY)*(-1 + pow2(Sbeta))) + 6*(-3*shiftst1 + shiftst2)*pow2(
        Mst1)*pow2(Sbeta)) + 2*Dmst12*pow2(Mst2)*((shiftst1 - shiftst2)*pow2(
        MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta)) + 12*(shiftst2*pow2(Mt) + (
        shiftst1 - shiftst2)*pow2(Mst1)*pow2(s2t))*pow2(Sbeta)) + 24*(shiftst1
        + shiftst2)*pow2(Mt)*pow2(Sbeta)*pow4(Mst2)))) + pow2(Sbeta)*(784*
        Dmst12*Mst1*Mt*s2t*(375*pow2(Mst1)*(-80*Dmsqst1*(Dmst12*((-98 + 24*
        lmMst1 - 24*lmMt)*pow2(Mst1) - 21*pow2(MuSUSY)) + 3*pow2(Mst2)*((74 -
        12*lmMst1 + 12*lmMt)*pow2(Mst1) + 7*pow2(MuSUSY))) + pow2(Msq)*(-8*
        pow2(Mst2)*(8*(347 - 50*lmMst1 + 66*lmMt - 66*lmMst1*lmMt + 183*pow2(
        lmMst1) + 108*pow2(lmMt))*pow2(Mst1) + (1361 + 10*lmMst1 + 54*pow2(
        lmMst1))*pow2(MuSUSY)) + Dmst12*(16*(-4378 - 517*lmMst1 + 243*lmMt -
        78*lmMst1*lmMt + 528*pow2(lmMst1) + 288*pow2(lmMt))*pow2(Mst1) + 3*(-
        10473 + 40*lmMst1 + 256*pow2(lmMst1))*pow2(MuSUSY)))) - 2*pow2(Dmglst1)
        *(3600*Dmsqst1*(397 - 30*lmMst1 + 30*lmMt)*pow2(Mst1)*(Dmst12 - pow2(
        Mst2)) + pow2(Msq)*(5*pow2(Mst2)*(-8*(-6041999 - 49410*lmMt + 6*lmMst1*
        (1721 + 1290*lmMt) + 111060*pow2(lmMst1) + 10800*pow2(lmMt))*pow2(Mst1)
        + 15*(954181 + 11256*lmMst1 - 11520*pow2(lmMst1))*pow2(MuSUSY)) +
        Dmst12*((263717842 - 633600*lmMt + 2400*lmMst1*(1043 + 93*lmMt) -
        525600*pow2(lmMst1) + 302400*pow2(lmMt))*pow2(Mst1) + 75*(-954181 -
        11256*lmMst1 + 11520*pow2(lmMst1))*pow2(MuSUSY)))) - 5*Dmglst1*Mst1*(
        240*Dmsqst1*(25*pow2(Mst2)*(2*(55 + 6*lmMst1 - 6*lmMt)*pow2(Mst1) + 5*
        pow2(MuSUSY)) - Dmst12*(4*(379 + 15*lmMst1 - 15*lmMt)*pow2(Mst1) + 125*
        pow2(MuSUSY))) + pow2(Msq)*(300*pow2(Mst2)*(16*(4964 - 3*lmMt - 5*
        lmMst1*(55 + 24*lmMt) + 21*pow2(lmMst1) + 36*pow2(lmMt))*pow2(Mst1) + (
        33261 - 532*lmMst1 - 660*pow2(lmMst1))*pow2(MuSUSY)) + Dmst12*(8*(
        4511549 + 9810*lmMt + 6*lmMst1*(14879 + 4710*lmMt) - 117360*pow2(
        lmMst1) - 21600*pow2(lmMt))*pow2(Mst1) - 15*(-949861 + 1944*lmMst1 +
        11520*pow2(lmMst1))*pow2(MuSUSY))))) - 29400*pow2(Dmst12)*pow2(Mst1)*
        pow2(s2t)*(2*pow2(Dmglst1)*(1080*Dmsqst1*pow2(Mst1) + pow2(Msq)*(6*(
        31901 - 5044*lmMst1 + 400*lmMt - 5100*pow2(lmMst1))*pow2(Mst1) + (
        655743 + 5288*lmMst1 - 11820*pow2(lmMst1))*pow2(MuSUSY))) - 10*pow2(
        Mst1)*(30*Dmsqst1*(12*(1 + 12*lmMst1)*pow2(Mst1) + (7 + 24*lmMst1)*
        pow2(MuSUSY)) + pow2(Msq)*(24*(613 - lmMst1 + 36*lmMt + 21*pow2(lmMst1)
        )*pow2(Mst1) + (-1429 + 454*lmMst1 - 24*pow2(lmMst1))*pow2(MuSUSY))) +
        5*Dmglst1*(Mst1*pow2(Msq)*(24*(-4515 + 596*lmMst1 - 48*lmMt + 516*pow2(
        lmMst1))*pow2(Mst1) + (84209 + 1264*lmMst1 - 240*pow2(lmMst1))*pow2(
        MuSUSY)) + 22560*Dmsqst1*pow3(Mst1))) - pow2(Mt)*(9408000*(10*Dmglst1*
        Dmsqst1*(-146 + 15*lmMst1 - 15*lmMt) - 45*Dmsqst1*Mst1*(14 + 6*lmMt -
        6*lmMst1*(3 + lmMt) + 3*pow2(lmMst1) + 3*pow2(lmMt)) - Dmglst1*(-15571
        + 852*lmMt + lmMst1*(-508 + 48*lmMt) + 978*pow2(lmMst1) + 612*pow2(
        lmMt))*pow2(Msq) - Mst1*pow2(Msq)*(-623 + 555*lmMt + (663 + 387*lmMt)*
        pow2(lmMst1) + 486*pow2(lmMt) - lmMst1*(1066 + 843*lmMt + 756*pow2(
        lmMt)) - 252*pow3(lmMst1) + 621*pow3(lmMt)))*pow3(Mst1)*pow4(Mst2) +
        16*pow2(Mst1)*(pow2(Dmst12)*(pow2(Dmglst1)*(35280*Dmsqst1*(457 + 510*
        lmMst1 - 510*lmMt) + (-13564884271 - 96403020*lmMt - 60*lmMst1*(968629
        + 101640*lmMt) + 288338400*pow2(lmMst1))*pow2(Msq)) + 98*Dmglst1*Mst1*(
        10800*Dmsqst1*(423 - 20*lmMst1 + 20*lmMt) + (-22574599 - 21060*lmMt +
        60*lmMst1*(6677 + 8880*lmMt) + 1598400*pow2(lmMst1) + 172800*pow2(lmMt)
        )*pow2(Msq)) + 980*(15*Dmsqst1*(4561 + 510*lmMst1 - 510*lmMt) + (93973
        - 61305*lmMt + 54*lmMst1*(2337 + 1255*lmMt) - 48420*pow2(lmMst1) -
        58050*pow2(lmMt))*pow2(Msq))*pow2(Mst1)) - 196*Dmst12*(10*Dmglst1*Mst1*
        (540*Dmsqst1*(423 - 20*lmMst1 + 20*lmMt) + (-5136871 - 29790*lmMt + 12*
        lmMst1*(8942 + 4305*lmMt) + 86040*pow2(lmMst1) + 21600*pow2(lmMt))*
        pow2(Msq)) + pow2(Dmglst1)*(180*Dmsqst1*(457 + 510*lmMst1 - 510*lmMt) -
        (73908751 + 1029690*lmMt + 540*lmMst1*(2243 + 295*lmMt) + 707400*pow2(
        lmMst1) + 64800*pow2(lmMt))*pow2(Msq)) - 750*(5*Dmsqst1*(79 + 204*
        lmMst1 + 12*lmMt) + (-4*lmMst1*(131 + 264*lmMt) + 264*pow2(lmMst1) +
        21*(783 + 14*lmMt + 30*pow2(lmMt)))*pow2(Msq))*pow2(Mst1))*pow2(Mst2) -
        7840*pow2(Dmglst1)*(315*Dmsqst1*(34 + 15*lmMst1 - 15*lmMt) + (-2055923
        - 61695*lmMt + 54*lmMst1*(-1162 + 445*lmMt) + 6345*pow2(lmMst1) -
        12150*pow2(lmMt))*pow2(Msq))*pow4(Mst2)) + 245*pow2(Msq)*pow2(MuSUSY)*(
        -400*Dmglst1*Mst1*(Dmst12*(403559 + 384*(lmMst1 + lmMt) - 4608*pow2(
        lmMst1))*(Dmst12 - pow2(Mst2)) - 24*(9631 + 16*lmMst1 + 48*lmMt - 192*
        pow2(lmMst1))*pow4(Mst2)) + 15*pow2(Mst1)*(pow2(Dmst12)*(852541 + 9216*
        lmMst1 - 10240*lmMt + 6144*pow2(lmMst1)) + 160*Dmst12*(11389 - 704*
        lmMst1 + 192*lmMt - 384*pow2(lmMst1))*pow2(Mst2) + 1920*(349 - 56*
        lmMst1 + 24*lmMt - 32*pow2(lmMst1))*pow4(Mst2)) - 32*pow2(Dmglst1)*(-2*
        Dmst12*(-6233611 - 58800*lmMst1 + 4560*lmMt + 14400*pow2(lmMst1))*(
        Dmst12 - pow2(Mst2)) + 25*(-230009 - 4704*lmMst1 + 480*lmMt + 1152*
        pow2(lmMst1))*pow4(Mst2)))))))) + pow2(Mt)*(Al4p*(-64*xDmglst1*pow2(
        Msq)*pow3(Dmglst1)*(8820*Mst1*twoLoopFlag*pow2(Msq)*(Dmst12*Mt*pow2(
        Mst2)*(40*s2t*Tbeta*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + Mt*((181 - 90*
        lmMst1 + 90*lmMt)*MuSUSY + 800*Mst1*Tbeta)*pow2(Sbeta) - 6*Mst1*s2t*(
        110*MuSUSY + (-17 + 30*lmMst1 - 30*lmMt)*Mst1*Tbeta)*pow2(Sbeta)) -
        pow2(Dmst12)*((((181 - 90*lmMst1 + 90*lmMt)*MuSUSY + 800*Mst1*Tbeta)*
        pow2(Mt) + 30*(2*MuSUSY + 11*Mst1*Tbeta)*pow2(Mst1)*pow2(s2t))*pow2(
        Sbeta) - 2*Mt*s2t*(-20*Tbeta*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + Mst1*(
        330*MuSUSY + (17 + 45*lmMst1 - 45*lmMt)*Mst1*Tbeta)*pow2(Sbeta))) + 4*(
        (48 - 45*lmMst1 + 45*lmMt)*MuSUSY + (377 - 30*lmMst1 + 30*lmMt)*Mst1*
        Tbeta)*pow2(Mt)*pow2(Sbeta)*pow4(Mst2)) + Al4p*threeLoopFlag*(2*Dmst12*
        Mst1*Mt*pow2(Mst2)*(11760*Dmsqst1*(1331 - 420*lmMst1 + 420*lmMt)*s2t*
        Tbeta*pow2(Mst1)*pow2(Sbeta) + pow2(Msq)*(5*Mt*(4*Mst1*Tbeta*(16826654
        - 1921395*lmMt - 3*lmMst1*(555463 + 2520*lmMt) + 4241160*pow2(lmMst1))
        + 3*MuSUSY*(1417174939 - 401268*lmMt + 4*lmMst1*(4061413 + 37800*lmMt)
        - 3185280*pow2(lmMst1) + 211680*pow2(lmMt)))*pow2(Sbeta) + 49*s2t*(-5*
        Tbeta*(-14217821 - 161940*lmMst1 + 28800*pow2(lmMst1))*pow2(MuSUSY)*(-1
        + pow2(Sbeta)) + 2*Mst1*(75*MuSUSY*(520781 + 17172*lmMst1 - 192*lmMt -
        11520*pow2(lmMst1)) + 4*Mst1*Tbeta*(27088246 + 5775*lmMt + 45*lmMst1*(
        12571 + 270*lmMt) - 136350*pow2(lmMst1) + 16200*pow2(lmMt)))*pow2(
        Sbeta)))) - pow2(Dmst12)*(Mst1*pow2(Msq)*((10*(4*Mst1*Tbeta*(16826654 -
        1921395*lmMt - 3*lmMst1*(555463 + 2520*lmMt) + 4241160*pow2(lmMst1)) +
        3*MuSUSY*(1417174939 - 401268*lmMt + 4*lmMst1*(4061413 + 37800*lmMt) -
        3185280*pow2(lmMst1) + 211680*pow2(lmMt)))*pow2(Mt) + 735*(MuSUSY*(
        14217821 + 161940*lmMst1 - 28800*pow2(lmMst1)) + 10*Mst1*Tbeta*(-520877
        - 17172*lmMst1 + 192*lmMt + 11520*pow2(lmMst1)))*pow2(Mst1)*pow2(s2t))*
        pow2(Sbeta) + 2*Mt*s2t*(-245*Tbeta*(-14217821 - 161940*lmMst1 + 28800*
        pow2(lmMst1))*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + Mst1*(-7350*MuSUSY*(-
        520781 - 17172*lmMst1 + 192*lmMt + 11520*pow2(lmMst1)) + Mst1*Tbeta*(-
        10642041163 + 11458020*lmMt + 60*lmMst1*(-346639 + 41580*lmMt) -
        5670000*pow2(lmMst1) + 3175200*pow2(lmMt)))*pow2(Sbeta))) + 23520*
        Dmsqst1*(1331 - 420*lmMst1 + 420*lmMt)*Mt*s2t*Tbeta*pow2(Sbeta)*pow3(
        Mst1)) + 196*pow2(Mt)*(120*Dmsqst1*Mst1*((1541 - 420*lmMst1 + 420*lmMt)
        *MuSUSY + (2579 + 120*lmMst1 - 120*lmMt)*Mst1*Tbeta)*pow2(Sbeta) +
        pow2(Msq)*(5*(2219399 + 9600*lmMst1 + 960*lmMt)*Tbeta*pow2(MuSUSY)*(-1
        + pow2(Sbeta)) + 2*Mst1*(Mst1*Tbeta*(10583177 + 60630*lmMt + 540*
        lmMst1*(188 + 395*lmMt) + 278100*pow2(lmMst1) - 59400*pow2(lmMt)) +
        MuSUSY*(54198467 + 43950*lmMt + 180*lmMst1*(6353 + 135*lmMt) - 272700*
        pow2(lmMst1) + 32400*pow2(lmMt)))*pow2(Sbeta)))*pow4(Mst2))) + 141120*
        twoLoopFlag*pow2(Mst1)*pow4(Msq)*(2*MuSUSY*pow2(Sbeta)*(-25*pow2(Mst1)*
        (4*Dmst12*Mt*(2*(5 + 6*lmMst1 + 3*lmMt)*Mt - 9*(1 + 2*lmMst1)*Mst1*s2t)
        *pow2(Mst2) - pow2(Dmst12)*(24*(1 - 3*lmMst1)*Mst1*Mt*s2t + 2*(5 + 6*(
        lmMst1 + lmMt))*pow2(Mt) + 9*(1 + 2*lmMst1)*pow2(Mst1)*pow2(s2t)) + 72*
        (1 + lmMst1 + lmMt)*pow2(Mt)*pow4(Mst2)) + 2*pow2(Dmglst1)*(3*Dmst12*
        Mt*((327 - 30*lmMst1 + 30*lmMt)*Mt + 50*(11 - 6*lmMst1)*Mst1*s2t)*pow2(
        Mst2) + 9*pow2(Dmst12)*(5*Mst1*Mt*s2t + (-109 + 10*lmMst1 - 10*lmMt)*
        pow2(Mt) - 25*pow2(Mst1)*pow2(s2t)) + 100*(17 - 3*lmMst1 + 3*lmMt)*
        pow2(Mt)*pow4(Mst2)) + Dmglst1*Mst1*(-300*Dmst12*Mt*((-3 + 6*lmMst1)*Mt
        + 2*(-1 + 6*lmMst1)*Mst1*s2t)*pow2(Mst2) + 3*pow2(Dmst12)*(-100*Mst1*
        Mt*s2t + 2*(-41 + 90*lmMst1 + 10*lmMt)*pow2(Mt) + 25*(-5 + 6*lmMst1)*
        pow2(Mst1)*pow2(s2t)) - 200*(-7 + 15*lmMst1 + 3*lmMt)*pow2(Mt)*pow4(
        Mst2))) + Tbeta*(25*Dmst12*s2t*(2*pow2(Dmglst1)*(12*Dmst12*Mt*pow2(
        MuSUSY) - 12*Mt*pow2(Mst2)*pow2(MuSUSY) - Dmst12*(-11 + 6*lmMst1)*Mst1*
        s2t*(pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 6*pow2(Mst1)*pow2(Sbeta))) - 4*
        Dmglst1*Mst1*((5 - 6*lmMst1)*Mt*pow2(Mst2)*pow2(MuSUSY) + 2*Dmst12*
        Mst1*s2t*((-1 + 3*lmMst1)*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 3*(-1 + 6*
        lmMst1)*pow2(Mst1)*pow2(Sbeta))) + pow2(Mst1)*(12*(1 + 2*lmMst1)*Mt*
        pow2(Mst2)*pow2(MuSUSY) + Dmst12*(4*Mt*pow2(MuSUSY) + 3*Mst1*s2t*pow2(
        MuSUSY)*(-1 + pow2(Sbeta)) + 36*(1 + 2*lmMst1)*s2t*pow2(Sbeta)*pow3(
        Mst1)))) + pow2(Sbeta)*(-4*Dmst12*Mt*s2t*(2*pow2(Dmglst1)*(-25*pow2(
        Mst2)*((31 - 6*lmMst1 + 6*lmMt)*pow2(Mst1) + 3*pow2(MuSUSY)) + Dmst12*(
        (307 - 105*lmMst1 + 105*lmMt)*pow2(Mst1) + 75*pow2(MuSUSY))) - 25*pow2(
        Mst1)*(Dmst12*(4*(1 + 3*lmMst1 + 6*lmMt)*pow2(Mst1) - pow2(MuSUSY)) -
        3*pow2(Mst2)*(6*(1 + 2*(lmMst1 + lmMt))*pow2(Mst1) + (1 + 2*lmMst1)*
        pow2(MuSUSY))) + 25*Dmglst1*(Mst1*pow2(Mst2)*(2*(-17 + 30*lmMst1 + 6*
        lmMt)*pow2(Mst1) + (-5 + 6*lmMst1)*pow2(MuSUSY)) - 4*Dmst12*(-4 + 6*
        lmMst1 + 3*lmMt)*pow3(Mst1))) - 2*Mst1*pow2(Mt)*(6*pow2(Dmglst1)*(35*
        pow2(Dmst12) + 3*Dmst12*(159 - 110*lmMst1 + 10*lmMt)*pow2(Mst2) + 200*(
        4 - 3*lmMst1)*pow4(Mst2)) - 4*Dmglst1*Mst1*((11 + 60*lmMst1 - 60*lmMt)*
        pow2(Dmst12) + 25*Dmst12*(1 + 30*lmMst1 + 6*lmMt)*pow2(Mst2) + 100*(1 +
        12*lmMst1 + 6*lmMt)*pow4(Mst2)) - 25*pow2(Mst1)*((5 + 42*lmMst1 + 30*
        lmMt)*pow2(Dmst12) - 8*Dmst12*(4 + 3*lmMst1 + 6*lmMt)*pow2(Mst2) - 72*(
        1 + lmMst1 - lmMt + 2*lmMst1*lmMt + pow2(lmMst1) - 3*pow2(lmMt))*pow4(
        Mst2))))))) + 3969000*oneLoopFlag*pow4(Msq)*(s2t*pow2(Dmst12)*(-12*Mt*
        MuSUSY*pow2(Sbeta) + s2t*Tbeta*(pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 12*
        pow2(Mst1)*pow2(Sbeta))) + 12*Mt*pow2(Sbeta)*(2*Dmst12*MuSUSY*s2t*pow2(
        Mst2) + Mt*Tbeta*(pow2(Dmst12) - 2*Dmst12*pow2(Mst2) + 4*(-lmMst1 +
        lmMt)*pow4(Mst2))))*pow5(Mst1))) + 4*Mst1*xDmst12*pow3(Dmst12)*(Mt*(
        992250*oneLoopFlag*(-(Mt*Tbeta*pow2(s2t)*(pow2(MuSUSY)*(-1 + pow2(
        Sbeta)) + 18*pow2(Mst1)*pow2(Sbeta))) - pow2(Sbeta)*(-8*MuSUSY*s2t*
        pow2(Mt) + 8*Tbeta*pow3(Mt) + MuSUSY*pow2(Mst1)*pow3(s2t)))*pow4(Msq)*
        pow4(Mst1) + Al4p*(-35280*Mst1*twoLoopFlag*(pow2(Mst1)*(25*Mst1*Mt*
        pow2(s2t)*((1 + 6*lmMst1)*Tbeta*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 12*
        Mst1*(MuSUSY + 3*lmMst1*MuSUSY + (1 + 12*lmMst1)*Mst1*Tbeta)*pow2(
        Sbeta)) + 100*s2t*pow2(Mt)*(-((5 + 6*lmMst1)*Tbeta*pow2(MuSUSY)*(-1 +
        pow2(Sbeta))) + 2*Mst1*(-12*(-1 + lmMst1)*MuSUSY + (1 + 3*lmMst1 + 9*
        lmMt)*Mst1*Tbeta)*pow2(Sbeta)) + pow2(Sbeta)*(-4*(50*(5 + 6*lmMst1)*
        MuSUSY + (137 - 330*lmMst1 - 270*lmMt)*Mst1*Tbeta)*pow3(Mt) + 75*(
        MuSUSY - 2*(1 + 2*lmMst1)*Mst1*Tbeta)*pow3(Mst1)*pow3(s2t))) + 2*
        Dmglst1*(Dmglst1*(-5*Mst1*Mt*pow2(s2t)*(10*(-11 + 6*lmMst1)*Tbeta*pow2(
        MuSUSY)*(-1 + pow2(Sbeta)) + 3*Mst1*(60*MuSUSY + (-107 + 60*lmMst1)*
        Mst1*Tbeta)*pow2(Sbeta)) + 4*s2t*pow2(Mt)*(-75*Tbeta*pow2(MuSUSY)*(-1 +
        pow2(Sbeta)) + Mst1*((870 - 450*lmMst1)*MuSUSY + (161 + 60*lmMst1 - 60*
        lmMt)*Mst1*Tbeta)*pow2(Sbeta)) + pow2(Sbeta)*(-6*((327 - 30*lmMst1 +
        30*lmMt)*MuSUSY + (547 - 330*lmMst1 + 30*lmMt)*Mst1*Tbeta)*pow3(Mt) +
        25*((11 - 6*lmMst1)*MuSUSY + 6*Mst1*Tbeta)*pow3(Mst1)*pow3(s2t))) +
        Mst1*(-25*Mst1*Mt*pow2(s2t)*(8*(-1 + 3*lmMst1)*Tbeta*pow2(MuSUSY)*(-1 +
        pow2(Sbeta)) + 3*Mst1*((5 - 6*lmMst1)*MuSUSY + 6*(-1 + 4*lmMst1)*Mst1*
        Tbeta)*pow2(Sbeta)) - 2*s2t*pow2(Mt)*(30*lmMst1*(5*Tbeta*pow2(MuSUSY)*(
        -1 + pow2(Sbeta)) + Mst1*(60*MuSUSY - 11*Mst1*Tbeta)*pow2(Sbeta)) +
        Tbeta*(-125*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 2*(131 - 135*lmMt)*pow2(
        Mst1)*pow2(Sbeta))) + pow2(Sbeta)*(4*(6*(17 - 30*lmMst1 + 5*lmMt)*
        MuSUSY + (47 + 870*lmMst1 + 30*lmMt)*Mst1*Tbeta)*pow3(Mt) + 25*((4 -
        12*lmMst1)*MuSUSY + (5 - 6*lmMst1)*Mst1*Tbeta)*pow3(Mst1)*pow3(s2t)))))
        *pow4(Msq) - 16*xDmglst1*pow2(Msq)*pow3(Dmglst1)*(-8820*twoLoopFlag*
        pow2(Msq)*(10*s2t*pow2(Mt)*(-4*Tbeta*pow2(MuSUSY)*(-1 + pow2(Sbeta)) +
        Mst1*(66*MuSUSY + 17*Mst1*Tbeta)*pow2(Sbeta)) + pow2(Sbeta)*(-60*Mt*(2*
        MuSUSY + 11*Mst1*Tbeta)*pow2(Mst1)*pow2(s2t) + ((-181 + 90*lmMst1 - 90*
        lmMt)*MuSUSY - 800*Mst1*Tbeta)*pow3(Mt) + 20*Tbeta*pow3(s2t)*pow4(Mst1)
        )) + Al4p*threeLoopFlag*(23520*Dmsqst1*(1331 - 420*lmMst1 + 420*lmMt)*
        s2t*Tbeta*pow2(Mst1)*pow2(Mt)*pow2(Sbeta) + 2*s2t*pow2(Msq)*pow2(Mt)*(-
        245*Tbeta*(-14217821 - 161940*lmMst1 + 28800*pow2(lmMst1))*pow2(MuSUSY)
        *(-1 + pow2(Sbeta)) + 2*Mst1*(3675*MuSUSY*(520781 + 17172*lmMst1 - 192*
        lmMt - 11520*pow2(lmMst1)) + Mst1*Tbeta*(-15951337379 + 10326120*lmMt +
        120*lmMst1*(-1097288 + 945*lmMt) + 21054600*pow2(lmMst1)))*pow2(Sbeta))
        - 5*pow2(Msq)*pow2(Sbeta)*(-294*Mt*(MuSUSY*(14217821 + 161940*lmMst1 -
        28800*pow2(lmMst1)) + 10*Mst1*Tbeta*(-520877 - 17172*lmMst1 + 192*lmMt
        + 11520*pow2(lmMst1)))*pow2(Mst1)*pow2(s2t) + 2*(4*Mst1*Tbeta*(-
        16826654 + 1921395*lmMt + 3*lmMst1*(555463 + 2520*lmMt) - 4241160*pow2(
        lmMst1)) + 3*MuSUSY*(-1417174939 + 401268*lmMt - 4*lmMst1*(4061413 +
        37800*lmMt) + 3185280*pow2(lmMst1) - 211680*pow2(lmMt)))*pow3(Mt) + 49*
        Tbeta*(14217821 + 161940*lmMst1 - 28800*pow2(lmMst1))*pow3(s2t)*pow4(
        Mst1)))))) + threeLoopFlag*pow2(Al4p)*(-(pow2(Msq)*(196*Mst1*Mt*MuSUSY*
        pow2(Sbeta)*(Dmglst1*Mst1*(2400*Dmsqst1*Mt*(7050*Mst1*Mt*s2t + 2*(557 +
        120*lmMst1 - 120*lmMt)*pow2(Mt) + 375*pow2(Mst1)*pow2(s2t)) - pow2(Msq)
        *(3600*Mst1*s2t*(12383 + 4128*lmMst1 + 80*lmMt - 1260*pow2(lmMst1))*
        pow2(Mt) + 225*Mt*(284641 + 8696*lmMst1 + 1680*pow2(lmMst1))*pow2(Mst1)
        *pow2(s2t) + 8*(193364399 + 90000*lmMt - 300*lmMst1*(3781 + 582*lmMt) -
        2005200*pow2(lmMst1) - 43200*pow2(lmMt))*pow3(Mt) + 375*(84209 + 1264*
        lmMst1 - 240*pow2(lmMst1))*pow3(Mst1)*pow3(s2t))) + 15*pow2(Mst1)*(500*
        Dmsqst1*(6*Mst1*s2t*(173 - 144*lmMst1*(-1 + shiftst2) + 216*shiftst2)*
        pow2(Mt) + 504*Mt*pow2(Mst1)*pow2(s2t) + 16*(65 - 6*lmMst1 + 6*lmMt)*
        pow3(Mt) + 3*(7 + 72*shiftst1 - 36*shiftst2 + 24*lmMst1*(1 - 2*shiftst1
        + shiftst2))*pow3(Mst1)*pow3(s2t)) - pow2(Msq)*(15*Mst1*s2t*(-102747 +
        640*lmMt + 6720*shiftst3 - 32*lmMst1*(331 + 90*shiftst3) + 13888*pow2(
        lmMst1))*pow2(Mt) - 75*Mt*(-20531 + 200*lmMst1 + 1200*pow2(lmMst1))*
        pow2(Mst1)*pow2(s2t) - 4*(-3454599 + 16840*lmMt + 48*lmMst1*(262 + 405*
        lmMt) + 46560*pow2(lmMst1))*pow3(Mt) + 50*(1429 - 720*shiftst1 + 360*
        shiftst2 - 234*shiftst3 + 2*lmMst1*(-227 + 720*shiftst1 - 360*shiftst2
        + 126*shiftst3) + 24*pow2(lmMst1))*pow3(Mst1)*pow3(s2t))) + 6*pow2(
        Dmglst1)*(54000*Dmsqst1*Mst1*s2t*pow2(Mt) + pow2(Msq)*(-600*Mst1*s2t*(
        274009 + 964*lmMst1 + 104*lmMt - 8310*pow2(lmMst1))*pow2(Mt) + 150*Mt*(
        954181 + 11256*lmMst1 - 11520*pow2(lmMst1))*pow2(Mst1)*pow2(s2t) + 4*(
        84334067 + 120*lmMst1*(2843 - 120*lmMt) + 202200*lmMt - 828000*pow2(
        lmMst1) - 21600*pow2(lmMt))*pow3(Mt) + 25*(-655743 - 5288*lmMst1 +
        11820*pow2(lmMst1))*pow3(Mst1)*pow3(s2t)))) + Tbeta*(3675*pow2(Mst1)*
        pow2(Mt)*pow2(s2t)*((8*Dmglst1*(2*Dmglst1*(655743 + 5288*lmMst1 -
        11820*pow2(lmMst1)) + 5*Mst1*(84209 + 1264*lmMst1 - 240*pow2(lmMst1)))*
        pow2(Msq) - 5*(480*Dmsqst1*(7 + 24*lmMst1) + (70121 + 2208*lmMst1 -
        432*pow2(lmMst1))*pow2(Msq))*pow2(Mst1))*pow2(MuSUSY) + 1440*shiftst1*
        pow2(Mst1)*(80*Dmsqst1*(3 - 2*lmMst1)*pow2(Mst1) + (1 - 2*lmMst1)*pow2(
        Msq)*(80*pow2(Mst1) + 3*pow2(MuSUSY)))*pow2(Sbeta)) - 245*(16*pow2(
        Dmglst1)*(2*Mt*(6233611 + 58800*lmMst1 - 4560*lmMt - 14400*pow2(lmMst1)
        ) + 15*Mst1*s2t*(954181 + 11256*lmMst1 - 11520*pow2(lmMst1)))*pow2(Msq)
        + 40*Dmglst1*Mst1*(6000*Dmsqst1*Mst1*s2t + (5*Mt*(403559 + 384*(lmMst1
        + lmMt) - 4608*pow2(lmMst1)) + 6*Mst1*s2t*(-1282471 + 7264*lmMst1 +
        18120*pow2(lmMst1)))*pow2(Msq)) + 15*(67200*Dmsqst1*Mst1*s2t + (80*
        Mst1*s2t*(-36863 + 80*lmMst1 + 552*pow2(lmMst1)) + Mt*(-1763661 +
        47104*lmMst1 - 5120*lmMt + 24576*pow2(lmMst1)))*pow2(Msq))*pow2(Mst1))*
        pow2(MuSUSY)*pow3(Mt) + 529200*pow4(Mst1)*(-(shiftst3*pow2(Msq)*(2*
        pow2(Mt)*pow2(s2t)*((-29 + 18*lmMst1)*pow2(MuSUSY)*(-1 + pow2(Sbeta)) +
        5*(-83 + 58*lmMst1)*pow2(Mst1)*pow2(Sbeta)) - 80*(-7 + 3*lmMst1)*pow2(
        Sbeta)*pow4(Mt) + 5*(1 - 2*lmMst1)*pow2(Sbeta)*pow4(Mst1)*pow4(s2t))) -
        10*(3*(1 - 2*lmMst1)*shiftst1*pow2(Msq)*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)
        - (1 - 2*lmMst1)*shiftst2*pow2(Msq)*pow2(s2t)*(2*pow2(Mt)*(pow2(MuSUSY)
        *(-1 + pow2(Sbeta)) - 10*pow2(Mst1)*pow2(Sbeta)) + 5*pow2(s2t)*pow2(
        Sbeta)*pow4(Mst1)) + 5*shiftst1*(Dmsqst1*(3 - 2*lmMst1) + (1 - 2*
        lmMst1)*pow2(Msq))*pow2(Sbeta)*pow4(Mst1)*pow4(s2t) + 5*Dmsqst1*(3 - 2*
        lmMst1)*shiftst2*(4*pow2(Mt)*pow2(s2t)*(-(pow2(MuSUSY)*(-1 + pow2(
        Sbeta))) + pow2(Mst1)*pow2(Sbeta)) + pow2(Sbeta)*(24*pow4(Mt) - pow4(
        Mst1)*pow4(s2t))))) + Mt*pow2(Sbeta)*(-392*Mst1*s2t*pow2(Mt)*(15*pow2(
        Mst1)*(2000*Dmsqst1*(2*(13 + 6*lmMst1 - 6*lmMt)*pow2(Mst1) - 21*pow2(
        MuSUSY)) - pow2(Msq)*((558619 + 76160*lmMt - 224*lmMst1*(1219 + 60*
        lmMt) + 123840*pow2(lmMst1) + 86400*pow2(lmMt))*pow2(Mst1) + 50*(-36863
        + 80*lmMst1 + 552*pow2(lmMst1))*pow2(MuSUSY))) + 2*pow2(Dmglst1)*(3600*
        Dmsqst1*(397 - 30*lmMst1 + 30*lmMt)*pow2(Mst1) + pow2(Msq)*(4*(
        192278911 + 177300*lmMt + 60*lmMst1*(19139 + 570*lmMt) - 1373400*pow2(
        lmMst1) + 43200*pow2(lmMt))*pow2(Mst1) + 75*(-954181 - 11256*lmMst1 +
        11520*pow2(lmMst1))*pow2(MuSUSY))) - 2*Dmglst1*Mst1*(600*Dmsqst1*(6*(47
        - 30*lmMst1 + 30*lmMt)*pow2(Mst1) + 125*pow2(MuSUSY)) + pow2(Msq)*((
        28188929 - 143100*lmMt - 3780*lmMst1*(549 + 80*lmMt) + 1389600*pow2(
        lmMst1) + 388800*pow2(lmMt))*pow2(Mst1) + 75*(-1282471 + 7264*lmMst1 +
        18120*pow2(lmMst1))*pow2(MuSUSY)))) - 3675*Mt*pow2(Mst1)*pow2(s2t)*(8*
        pow2(Dmglst1)*(2160*Dmsqst1*pow2(Mst1) + pow2(Msq)*(3*(-452211 - 17060*
        lmMst1 + 992*lmMt + 1320*pow2(lmMst1))*pow2(Mst1) + 2*(655743 + 5288*
        lmMst1 - 11820*pow2(lmMst1))*pow2(MuSUSY))) - 5*pow2(Mst1)*(120*
        Dmsqst1*((-137 + 288*lmMst1)*pow2(Mst1) + 4*(7 + 24*lmMst1)*pow2(
        MuSUSY)) + pow2(Msq)*(24*(2785 - 304*lmMst1 + 384*lmMt + 768*pow2(
        lmMst1))*pow2(Mst1) + (70121 + 2208*lmMst1 - 432*pow2(lmMst1))*pow2(
        MuSUSY))) + 8*Dmglst1*(Mst1*pow2(Msq)*(24*(-20017 + 1203*lmMst1 - 200*
        lmMt + 2250*pow2(lmMst1))*pow2(Mst1) + 5*(84209 + 1264*lmMst1 - 240*
        pow2(lmMst1))*pow2(MuSUSY)) + 112800*Dmsqst1*pow3(Mst1))) - pow3(Mt)*(
        16*pow2(Dmglst1)*(17640*Dmsqst1*(457 + 510*lmMst1 - 510*lmMt)*pow2(
        Mst1) + pow2(Msq)*((-6321826673 + 60581820*lmMst1 + 4506600*lmMt +
        9513000*lmMst1*lmMt + 357663600*pow2(lmMst1) + 6350400*pow2(lmMt))*
        pow2(Mst1) + 490*(-6233611 - 58800*lmMst1 + 4560*lmMt + 14400*pow2(
        lmMst1))*pow2(MuSUSY))) + 392*Dmglst1*(Mst1*pow2(Msq)*(12*(9598037 +
        92280*lmMt + 20*lmMst1*(-11207 + 270*lmMt) + 246000*pow2(lmMst1) -
        14400*pow2(lmMt))*pow2(Mst1) - 125*(403559 + 384*(lmMst1 + lmMt) -
        4608*pow2(lmMst1))*pow2(MuSUSY)) + 21600*Dmsqst1*(423 - 20*lmMst1 + 20*
        lmMt)*pow3(Mst1)) + 49*(pow2(Msq)*pow2(Mst1)*((83430364 - 8607840*lmMt
        + 480*lmMst1*(36107 + 13380*lmMt) - 9273600*pow2(lmMst1) - 6652800*
        pow2(lmMt))*pow2(Mst1) + 75*(1763661 - 47104*lmMst1 + 5120*lmMt -
        24576*pow2(lmMst1))*pow2(MuSUSY)) + 9600*Dmsqst1*(3268 + 2805*lmMst1 -
        105*lmMt)*pow4(Mst1))) - 29400*(pow2(Dmglst1)*(954181 + 11256*lmMst1 -
        11520*pow2(lmMst1))*pow2(Msq) + 10*Dmglst1*Mst1*(100*Dmsqst1 + (33261 -
        532*lmMst1 - 660*pow2(lmMst1))*pow2(Msq)) + 20*(210*Dmsqst1 + (1361 +
        10*lmMst1 + 54*pow2(lmMst1))*pow2(Msq))*pow2(Mst1))*pow3(s2t)*pow5(
        Mst1))))) + 5880*xDmsqst1*pow2(Dmsqst1)*pow2(Mst1)*(16*pow2(Dmglst1)*
        pow2(Mt)*(Mt*s2t*(-675*MuSUSY - 30*(-397 + 30*lmMst1 - 30*lmMt)*Mst1*
        Tbeta + 4*Dmglst1*(-1331 + 420*lmMst1 - 420*lmMt)*Tbeta*xDmglst1) + 3*(
        457 + 510*lmMst1 - 510*lmMt)*Tbeta*pow2(Mt) + 675*Tbeta*pow2(Mst1)*
        pow2(s2t))*pow2(Sbeta) + 40*Dmglst1*Mt*(-2*s2t*pow2(Mt)*(125*Tbeta*
        pow2(MuSUSY)*(-1 + pow2(Sbeta)) + Mst1*(7050*MuSUSY + (-193 - 330*
        lmMst1 + 330*lmMt)*Mst1*Tbeta)*pow2(Sbeta)) + pow2(Sbeta)*(150*Mt*(-5*
        MuSUSY + 94*Mst1*Tbeta)*pow2(Mst1)*pow2(s2t) - 4*((557 + 120*lmMst1 -
        120*lmMt)*MuSUSY + 9*(-423 + 20*lmMst1 - 20*lmMt)*Mst1*Tbeta)*pow3(Mt)
        + 125*Tbeta*pow3(s2t)*pow4(Mst1))) + 5*Mst1*(-75*Mst1*pow2(Mt)*pow2(
        s2t)*(4*(7 + 6*shiftst1 - 24*lmMst1*(-1 + shiftst2) + 30*shiftst2)*
        Tbeta*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 3*Mst1*(112*MuSUSY + Mst1*(91 +
        240*shiftst1 + 32*lmMst1*(3 - 4*shiftst1 + shiftst2))*Tbeta)*pow2(
        Sbeta)) - 150*s2t*(56*Tbeta*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + Mst1*(3*
        MuSUSY*(47 + 96*(lmMst1 + shiftst2 - lmMst1*shiftst2)) - 208*Mst1*
        Tbeta)*pow2(Sbeta))*pow3(Mt) + pow2(Sbeta)*(-150*Mt*(MuSUSY*(7 + 108*
        shiftst1 - 72*shiftst2 + 24*lmMst1*(1 - 2*shiftst1 + shiftst2)) - 28*
        Mst1*Tbeta)*pow3(Mst1)*pow3(s2t) - 4*(200*(65 - 6*lmMst1 + 6*lmMt)*
        MuSUSY - 3*Mst1*(1999 - 90*lmMt + 90*lmMst1*(41 - 40*shiftst2) + 3600*
        shiftst2)*Tbeta)*pow4(Mt) - 1800*(-2 + lmMst1)*(shiftst1 - shiftst2)*
        Tbeta*pow4(s2t)*pow5(Mst1)))))) + 11025*threeLoopFlag*pow2(Al4p)*(384*
        z2*(pow2(Mst2)*(-50*Mt*xDmsqst1*pow2(Dmsqst1)*(-(pow2(Dmst12)*pow2(s2t)
        *(-(Mt*(shiftst1 + shiftst2)*Tbeta*pow2(MuSUSY)*(-1 + pow2(Sbeta))) -
        3*(MuSUSY*s2t*(shiftst1 - shiftst2) + 2*Mt*(3*shiftst1 - shiftst2)*
        Tbeta)*pow2(Mst1)*pow2(Sbeta))) + 2*Dmst12*Mt*pow2(Mst2)*(12*Mt*MuSUSY*
        s2t*shiftst2*pow2(Sbeta) - 12*shiftst2*Tbeta*pow2(Mt)*pow2(Sbeta) - (
        shiftst1 - shiftst2)*Tbeta*pow2(s2t)*(pow2(MuSUSY)*(-1 + pow2(Sbeta)) +
        12*pow2(Mst1)*pow2(Sbeta))) - 24*(MuSUSY*s2t*(shiftst1 - shiftst2) +
        Mt*(shiftst1 + shiftst2)*Tbeta)*pow2(Mt)*pow2(Sbeta)*pow4(Mst2)) +
        pow2(Msq)*(50*Tbeta*(Dmsqst1 + pow2(Msq))*pow2(Mt)*(-(pow2(Dmst12)*
        pow2(s2t)*((shiftst1 + shiftst2)*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 6*(
        3*shiftst1 - shiftst2)*pow2(Mst1)*pow2(Sbeta))) + 2*Dmst12*pow2(Mst2)*(
        (shiftst1 - shiftst2)*pow2(MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta)) + 12*(
        shiftst2*pow2(Mt) + (shiftst1 - shiftst2)*pow2(Mst1)*pow2(s2t))*pow2(
        Sbeta)) + 24*(shiftst1 + shiftst2)*pow2(Mt)*pow2(Sbeta)*pow4(Mst2)) +
        15*MuSUSY*pow2(Sbeta)*(10*(Dmsqst1 + pow2(Msq))*(-8*Dmst12*s2t*
        shiftst2*pow2(Mst2)*pow3(Mt) + Mt*(-shiftst1 + shiftst2)*pow2(Dmst12)*
        pow2(Mst1)*pow3(s2t) + 8*s2t*(shiftst1 - shiftst2)*pow3(Mt)*pow4(Mst2))
        + shiftst3*pow2(Msq)*(-(Mt*pow2(Dmst12)*pow2(Mst1)*pow3(s2t)) + 8*s2t*
        pow3(Mt)*(pow2(Dmst12) - Dmst12*pow2(Mst2) + pow4(Mst2))))) + 5*
        shiftst3*Tbeta*pow4(Msq)*(-(Dmst12*pow2(Mt)*pow2(s2t)*((3*Dmst12 - 2*
        pow2(Mst2))*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 6*pow2(Mst1)*(7*Dmst12 -
        4*pow2(Mst2))*pow2(Sbeta))) + 24*pow2(Sbeta)*(pow2(Dmst12) - Dmst12*
        pow2(Mst2) + pow4(Mst2))*pow4(Mt))) - xDmst12*pow3(Dmst12)*(50*Dmsqst1*
        (Dmsqst1*xDmsqst1 + pow2(Msq))*(-4*shiftst2*pow2(Mt)*(Tbeta*pow2(
        MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta)) + 6*Mt*(MuSUSY*s2t - Mt*Tbeta)*
        pow2(Sbeta)) + pow2(Sbeta)*(2*Mt*(MuSUSY*s2t*(-2*shiftst1 + shiftst2) +
        2*Mt*(-4*shiftst1 + shiftst2)*Tbeta)*pow2(Mst1)*pow2(s2t) + (shiftst1 -
        shiftst2)*Tbeta*pow4(Mst1)*pow4(s2t))) + pow4(Msq)*(2*pow2(Mt)*(-((15*
        shiftst1 + 10*shiftst2 + 9*shiftst3)*Tbeta*pow2(MuSUSY)*pow2(s2t)*(-1 +
        pow2(Sbeta))) + 60*Mt*shiftst3*(MuSUSY*s2t + Mt*Tbeta)*pow2(Sbeta)) +
        pow2(Sbeta)*(-5*Mt*(MuSUSY*s2t*(40*shiftst1 - 20*shiftst2 + 7*shiftst3)
        + 2*Mt*(80*shiftst1 - 20*shiftst2 + 29*shiftst3)*Tbeta)*pow2(Mst1)*
        pow2(s2t) + 5*(10*shiftst1 - 10*shiftst2 + shiftst3)*Tbeta*pow4(Mst1)*
        pow4(s2t)))))*pow5(Mst1) - 5*Mt*z3*(Mt*pow2(Mst2)*(8*Dmst12*s2t*Tbeta*
        pow2(Msq)*pow2(Mst1)*pow2(Sbeta)*(-2*Mt*pow2(Msq)*(pow2(Mst1)*(12*
        Dmst12*(22100*Dmglst1*Mst1 + 64701*pow2(Dmglst1) + 3290*pow2(Mst1)) +
        96*(1963*Dmglst1*Mst1 + 7488*pow2(Dmglst1) + 106*pow2(Mst1))*pow2(Mst2)
        ) + (-210810*pow2(Dmglst1)*(Dmst12 - pow2(Mst2)) + pow2(Mst1)*(17285*
        Dmst12 + 6024*pow2(Mst2)) + Dmglst1*Mst1*(105405*Dmst12 + 75164*pow2(
        Mst2)))*pow2(MuSUSY)) + Dmst12*Mst1*s2t*(18*pow2(Mst1)*(11*pow2(Msq)*(
        40*pow2(Mst1) - 13*pow2(MuSUSY)) + 70*Dmsqst1*(4*pow2(Mst1) + pow2(
        MuSUSY))) - Dmglst1*pow2(Msq)*(74316*Dmglst1*pow2(Mst1) + 149630*
        Dmglst1*pow2(MuSUSY) + 50485*Mst1*pow2(MuSUSY) - 34776*pow3(Mst1)))) -
        64*xDmglst1*pow3(Dmglst1)*pow4(Msq)*(Dmst12*Mst1*(-(Dmst12*((9*(71189*
        MuSUSY + 1240*Mst1*Tbeta)*pow2(Mt) + 33*(4789*MuSUSY - 1738*Mst1*Tbeta)
        *pow2(Mst1)*pow2(s2t))*pow2(Sbeta) + Mt*s2t*(105358*Tbeta*pow2(MuSUSY)*
        (-1 + pow2(Sbeta)) + 3*Mst1*(38236*MuSUSY - 107331*Mst1*Tbeta)*pow2(
        Sbeta)))) + Mt*pow2(Mst2)*(9*Mt*(71189*MuSUSY + 1240*Mst1*Tbeta)*pow2(
        Sbeta) + 2*s2t*(52679*Tbeta*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 6*Mst1*(
        9559*MuSUSY + 26559*Mst1*Tbeta)*pow2(Sbeta)))) + 6*pow2(Mt)*(5511*
        Tbeta*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 2*Mst1*(26559*MuSUSY + 5282*
        Mst1*Tbeta)*pow2(Sbeta))*pow4(Mst2)) + Mst1*pow2(Msq)*(-3*Tbeta*pow2(
        Mt)*pow2(Sbeta)*(20160*Dmsqst1*Dmst12*(Dmst12 + 2*pow2(Mst2))*pow4(
        Mst1) + pow2(Msq)*pow2(Mst1)*(160*Dmst12*pow2(Mst2)*(2052*pow2(Mst1) +
        425*pow2(MuSUSY)) + pow2(Dmst12)*(29472*pow2(Mst1) + 31963*pow2(MuSUSY)
        ) + 512*((-214 + 24*lmMst1 - 24*lmMt)*pow2(Mst1) + 47*pow2(MuSUSY))*
        pow4(Mst2)) + 48*Dmglst1*pow2(Msq)*(Mst1*(-(pow2(Dmst12)*(6970*pow2(
        Mst1) + 8213*pow2(MuSUSY))) + Dmst12*pow2(Mst2)*(30912*pow2(Mst1) +
        8213*pow2(MuSUSY)) + 8*(3356*pow2(Mst1) + 583*pow2(MuSUSY))*pow4(Mst2))
        + Dmglst1*(4*Dmst12*pow2(Mst2)*(11332*pow2(Mst1) + 5113*pow2(MuSUSY)) -
        pow2(Dmst12)*(45397*pow2(Mst1) + 20452*pow2(MuSUSY)) + 10*(4840*pow2(
        Mst1) + 943*pow2(MuSUSY))*pow4(Mst2)))) - Tbeta*pow2(MuSUSY)*(-8*pow2(
        Dmst12)*pow2(Mst1)*(5*Dmglst1*(29926*Dmglst1 + 10097*Mst1)*pow2(Msq) -
        18*(70*Dmsqst1 - 143*pow2(Msq))*pow2(Mst1))*pow2(s2t) + Mt*pow2(Msq)*(
        16*Dmst12*Mst1*s2t*(210810*pow2(Dmglst1)*(Dmst12 - pow2(Mst2)) - pow2(
        Mst1)*(17285*Dmst12 + 6024*pow2(Mst2)) - Dmglst1*Mst1*(105405*Dmst12 +
        75164*pow2(Mst2))) + 3*Mt*(96*pow2(Dmglst1)*(10226*pow2(Dmst12) -
        10226*Dmst12*pow2(Mst2) - 4715*pow4(Mst2)) + 48*Dmglst1*Mst1*(8213*
        pow2(Dmst12) - 8213*Dmst12*pow2(Mst2) - 4664*pow4(Mst2)) - pow2(Mst1)*(
        31963*pow2(Dmst12) + 68000*Dmst12*pow2(Mst2) + 24064*pow4(Mst2))))) +
        16*Mst1*MuSUSY*pow2(Sbeta)*(-6*Dmst12*Mst1*Mt*s2t*(2*Dmglst1*pow2(Msq)*
        (28677*Dmglst1*Dmst12 + 3821*Dmst12*Mst1 + 6193*Dmglst1*pow2(Mst2) -
        2898*Mst1*pow2(Mst2)) - 3*pow2(Mst1)*(35*Dmsqst1*(7*Dmst12 + 8*pow2(
        Mst2)) + pow2(Msq)*(403*Dmst12 + 440*pow2(Mst2)))) + pow2(Msq)*(3*pow2(
        Dmst12)*pow2(Mst1)*(37582*Dmglst1*Mst1 + 105405*pow2(Dmglst1) + 3012*
        pow2(Mst1))*pow2(s2t) + pow2(Mt)*(1404*pow2(Dmglst1)*(1065*pow2(Dmst12)
        - 1065*Dmst12*pow2(Mst2) - 512*pow4(Mst2)) - pow2(Mst1)*(51181*pow2(
        Dmst12) + 49656*Dmst12*pow2(Mst2) + 10176*pow4(Mst2)) - 4*Dmglst1*Mst1*
        (86833*pow2(Dmst12) + 113412*Dmst12*pow2(Mst2) + 47112*pow4(Mst2))))))
        - 5040*xDmsqst1*pow2(Dmsqst1)*(6*Dmst12*Mt*(-5*MuSUSY*s2t + 6*Mt*Tbeta)
        *pow2(Mst2)*pow2(Sbeta) + Tbeta*pow2(Dmst12)*pow2(s2t)*(-2*pow2(MuSUSY)
        *(-1 + pow2(Sbeta)) - 15*pow2(Mst1)*pow2(Sbeta)) + 24*Tbeta*pow2(Mt)*
        pow2(Sbeta)*pow4(Mst2))*pow5(Mst1)) - 2*Mst1*xDmst12*pow3(Dmst12)*(
        2520*Dmsqst1*(Dmsqst1*xDmsqst1*(Mt*Tbeta*pow2(s2t)*(4*pow2(MuSUSY)*(-1
        + pow2(Sbeta)) + 15*pow2(Mst1)*pow2(Sbeta)) - 36*Tbeta*pow2(Sbeta)*
        pow3(Mt) + 2*MuSUSY*pow2(Sbeta)*(15*s2t*pow2(Mt) + pow2(Mst1)*pow3(s2t)
        )) + pow2(Msq)*(Mt*Tbeta*pow2(s2t)*(4*pow2(MuSUSY)*(-1 + pow2(Sbeta)) +
        pow2(Mst1)*pow2(Sbeta)) - 48*Tbeta*pow2(Sbeta)*pow3(Mt) + 2*MuSUSY*
        pow2(Sbeta)*(22*s2t*pow2(Mt) + pow2(Mst1)*pow3(s2t))))*pow4(Mst1) +
        pow4(Msq)*(-32*xDmglst1*pow3(Dmglst1)*(-2*s2t*pow2(Mt)*(52679*Tbeta*
        pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 237*Mst1*(-242*MuSUSY + 2031*Mst1*
        Tbeta)*pow2(Sbeta)) + pow2(Sbeta)*(66*Mt*(-4789*MuSUSY + 1738*Mst1*
        Tbeta)*pow2(Mst1)*pow2(s2t) - 9*(71189*MuSUSY + 1240*Mst1*Tbeta)*pow3(
        Mt) + 52679*Tbeta*pow3(s2t)*pow4(Mst1))) - 4*Dmglst1*Mst1*(4*Mst1*
        MuSUSY*pow2(Mt)*(142987*MuSUSY*s2t*Tbeta*(-1 + pow2(Sbeta)) + 574156*
        Mt*pow2(Sbeta)) + 2*Mt*pow2(Mst1)*(50485*Tbeta*pow2(MuSUSY)*pow2(s2t)*(
        -1 + pow2(Sbeta)) + 12*Mt*(4744*MuSUSY*s2t + 12729*Mt*Tbeta)*pow2(
        Sbeta)) + s2t*(Mt*(90723*MuSUSY*s2t - 164264*Mt*Tbeta) + Mst1*s2t*(
        50485*MuSUSY*s2t - 80628*Mt*Tbeta))*pow2(Sbeta)*pow3(Mst1) + 86*Tbeta*(
        -1719*pow2(MuSUSY)*(-1 + pow2(Sbeta))*pow3(Mt) + 874*pow2(Sbeta)*pow3(
        s2t)*pow5(Mst1))) - pow2(Mst1)*(16*Mst1*MuSUSY*pow2(Mt)*(20297*MuSUSY*
        s2t*Tbeta*(-1 + pow2(Sbeta)) + 76009*Mt*pow2(Sbeta)) - Mt*pow2(Mst1)*(
        37669*Tbeta*pow2(MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta)) - 4*Mt*(-32783*
        MuSUSY*s2t + 33933*Mt*Tbeta)*pow2(Sbeta)) + 4*s2t*(Mt*(33783*MuSUSY*s2t
        - 23402*Mt*Tbeta) + 18*Mst1*s2t*(143*MuSUSY*s2t - 37*Mt*Tbeta))*pow2(
        Sbeta)*pow3(Mst1) + Tbeta*(197889*pow2(MuSUSY)*(-1 + pow2(Sbeta))*pow3(
        Mt) + 24096*pow2(Sbeta)*pow3(s2t)*pow5(Mst1))) - 8*pow2(Dmglst1)*(-30*
        Mst1*MuSUSY*pow2(Mt)*(7027*MuSUSY*s2t*Tbeta*(-1 + pow2(Sbeta)) + 49842*
        Mt*pow2(Sbeta)) + 2*Mt*pow2(Mst1)*(74815*Tbeta*pow2(MuSUSY)*pow2(s2t)*(
        -1 + pow2(Sbeta)) + 3*Mt*(127094*MuSUSY*s2t - 68199*Mt*Tbeta)*pow2(
        Sbeta)) + s2t*(-6*Mt*s2t*(105405*MuSUSY + 22484*Mst1*Tbeta) + 2271672*
        Tbeta*pow2(Mt) + 74815*Mst1*MuSUSY*pow2(s2t))*pow2(Sbeta)*pow3(Mst1) +
        Tbeta*(-184068*pow2(MuSUSY)*(-1 + pow2(Sbeta))*pow3(Mt) + 105405*pow2(
        Sbeta)*pow3(s2t)*pow5(Mst1))))))))/(1.90512e8*Tbeta*pow2(Sbeta)*pow4(
        Msq)*pow5(Mst1)*pow6(Mst2));
}

/**
 * @return The off-diagonal (1, 2) = (2, 1) matrix element of the Higgs mass matrix as a double for the hierarchy 'H3q22g'
 */
double H3q22g::getS12() const {
   return -(MuSUSY*(40500*Dmst12*Mt*oneLoopFlag*s2t*pow4(Msq)*(-2*Dmst12*Mt*(
        MuSUSY*s2t + 6*Mt*Tbeta)*pow2(Mst2) + pow2(Dmst12)*(2*Mt*MuSUSY*s2t +
        8*Tbeta*pow2(Mt) - Tbeta*pow2(Mst1)*pow2(s2t)) + 24*Tbeta*pow2(Mt)*
        pow4(Mst2))*pow5(Mst1) + 32*Al4p*pow2(Msq)*(4*xDmglst1*pow3(Dmglst1)*(
        45*Mst1*twoLoopFlag*pow2(Msq)*pow2(Mt)*(80*Dmst12*Mt*MuSUSY*s2t*(pow2(
        Dmst12) - Dmst12*pow2(Mst2) + pow4(Mst2)) - Tbeta*(pow2(Dmst12)*pow2(
        Mst2)*(660*Mst1*Mt*s2t + (-181 + 90*lmMst1 - 90*lmMt)*pow2(Mt) - 60*
        pow2(Mst1)*pow2(s2t)) + (-660*Mst1*Mt*s2t + (181 - 90*lmMst1 + 90*lmMt)
        *pow2(Mt) + 120*pow2(Mst1)*pow2(s2t))*pow3(Dmst12) + Dmst12*Mt*((181 -
        90*lmMst1 + 90*lmMt)*Mt - 660*Mst1*s2t)*pow4(Mst2) + 12*(16 - 15*lmMst1
        + 15*lmMt)*pow2(Mt)*pow6(Mst2))) + Al4p*threeLoopFlag*(5*MuSUSY*pow2(
        Msq)*pow3(Mt)*(Dmst12*Mst1*s2t*(14217821 + 161940*lmMst1 - 28800*pow2(
        lmMst1))*(pow2(Dmst12) - Dmst12*pow2(Mst2) + pow4(Mst2)) + 2*(2219399 +
        9600*lmMst1 + 960*lmMt)*Mt*pow6(Mst2)) - (Mst1*Tbeta*pow2(Mt)*(pow2(
        Msq)*(-15*pow2(Dmst12)*pow2(Mst2)*(980*Mst1*Mt*s2t*(520781 + 17172*
        lmMst1 - 192*lmMt - 11520*pow2(lmMst1)) + (2834349878 - 802536*lmMt +
        8*lmMst1*(4061413 + 37800*lmMt) - 6370560*pow2(lmMst1) + 423360*pow2(
        lmMt))*pow2(Mt) + 49*(14217821 + 161940*lmMst1 - 28800*pow2(lmMst1))*
        pow2(Mst1)*pow2(s2t)) + 30*((490*Mst1*Mt*s2t*(520781 + 17172*lmMst1 -
        192*lmMt - 11520*pow2(lmMst1)) + (1417174939 - 401268*lmMt + 4*lmMst1*(
        4061413 + 37800*lmMt) - 3185280*pow2(lmMst1) + 211680*pow2(lmMt))*pow2(
        Mt) + 49*(14217821 + 161940*lmMst1 - 28800*pow2(lmMst1))*pow2(Mst1)*
        pow2(s2t))*pow3(Dmst12) + Dmst12*Mt*(490*Mst1*s2t*(520781 + 17172*
        lmMst1 - 192*lmMt - 11520*pow2(lmMst1)) + Mt*(1417174939 - 401268*lmMt
        + 4*lmMst1*(4061413 + 37800*lmMt) - 3185280*pow2(lmMst1) + 211680*pow2(
        lmMt)))*pow4(Mst2))) + 392*(60*Dmsqst1*(1541 - 420*lmMst1 + 420*lmMt) +
        (54198467 + 43950*lmMt + 180*lmMst1*(6353 + 135*lmMt) - 272700*pow2(
        lmMst1) + 32400*pow2(lmMt))*pow2(Msq))*pow2(Mt)*pow6(Mst2)))/196.)) +
        45*Mt*twoLoopFlag*pow2(Msq)*pow2(Mst1)*(-25*s2t*pow3(Dmst12)*(8*(5 + 6*
        lmMst1)*MuSUSY*pow2(Mst1)*pow2(Mt) + 2*pow2(Dmglst1)*(24*MuSUSY*pow2(
        Mt) + (11 - 6*lmMst1)*Tbeta*pow2(s2t)*pow3(Mst1)) - 8*Dmglst1*((5 - 6*
        lmMst1)*Mst1*MuSUSY*pow2(Mt) + (-1 + 3*lmMst1)*Tbeta*pow2(s2t)*pow4(
        Mst1)) + 3*Tbeta*pow2(s2t)*pow5(Mst1)) + Mt*(50*Dmst12*MuSUSY*s2t*(
        Mst1*s2t*pow2(Dmst12)*(16*Dmglst1*(1 - 3*lmMst1)*Mst1 + (44 - 24*
        lmMst1)*pow2(Dmglst1) + (1 + 6*lmMst1)*pow2(Mst1)) + 4*Dmst12*Mt*(6*
        pow2(Dmglst1) + pow2(Mst1))*pow2(Mst2) - Dmst12*Mst1*s2t*(8*Dmglst1*(1
        - 3*lmMst1)*Mst1 + (22 - 12*lmMst1)*pow2(Dmglst1) + 3*pow2(Mst1))*pow2(
        Mst2) - 4*Mt*(Dmglst1*(5 - 6*lmMst1)*Mst1 + 6*pow2(Dmglst1) - 3*(1 + 2*
        lmMst1)*pow2(Mst1))*pow4(Mst2)) + 2*Tbeta*(-25*pow2(Mst1)*(-(pow2(
        Dmst12)*pow2(Mst2)*(24*(1 - 3*lmMst1)*Mst1*Mt*s2t + 2*(5 + 6*lmMst1 +
        6*lmMt)*pow2(Mt) + 9*(1 + 2*lmMst1)*pow2(Mst1)*pow2(s2t))) + (-48*(-1 +
        lmMst1)*Mst1*Mt*s2t - 4*(5 + 6*lmMst1)*pow2(Mt) + 6*(1 + 3*lmMst1)*
        pow2(Mst1)*pow2(s2t))*pow3(Dmst12) + 4*Dmst12*Mt*(2*(5 + 6*lmMst1 + 3*
        lmMt)*Mt - 9*(1 + 2*lmMst1)*Mst1*s2t)*pow4(Mst2) + 72*(1 + lmMst1 +
        lmMt)*pow2(Mt)*pow6(Mst2)) + 2*pow2(Dmglst1)*(9*pow2(Dmst12)*pow2(Mst2)
        *(5*Mst1*Mt*s2t + (-109 + 10*lmMst1 - 10*lmMt)*pow2(Mt) - 25*pow2(Mst1)
        *pow2(s2t)) + (60*(-29 + 15*lmMst1)*Mst1*Mt*s2t + (981 - 90*lmMst1 +
        90*lmMt)*pow2(Mt) + 450*pow2(Mst1)*pow2(s2t))*pow3(Dmst12) + 3*Dmst12*
        Mt*((327 - 30*lmMst1 + 30*lmMt)*Mt + 50*(11 - 6*lmMst1)*Mst1*s2t)*pow4(
        Mst2) + 100*(17 - 3*lmMst1 + 3*lmMt)*pow2(Mt)*pow6(Mst2)) + Dmglst1*
        Mst1*(-3*pow2(Dmst12)*(pow2(Mst2)*(100*Mst1*Mt*s2t - 2*(-41 + 90*lmMst1
        + 10*lmMt)*pow2(Mt) + 25*(5 - 6*lmMst1)*pow2(Mst1)*pow2(s2t)) + Dmst12*
        (-1200*lmMst1*Mst1*Mt*s2t + 8*(17 - 30*lmMst1 + 5*lmMt)*pow2(Mt) + 25*(
        -5 + 6*lmMst1)*pow2(Mst1)*pow2(s2t))) + 300*Dmst12*Mt*((3 - 6*lmMst1)*
        Mt + 2*(1 - 6*lmMst1)*Mst1*s2t)*pow4(Mst2) - 200*(-7 + 15*lmMst1 + 3*
        lmMt)*pow2(Mt)*pow6(Mst2)))))) + threeLoopFlag*pow2(Al4p)*(4*Mst1*pow2(
        Msq)*((5*MuSUSY*pow2(Mt)*(16*pow2(Dmglst1)*pow2(Msq)*(pow2(Dmst12)*
        pow2(Mst2)*(30*Mst1*Mt*s2t*(-954181 - 11256*lmMst1 + 11520*pow2(lmMst1)
        ) + 4*(-6233611 - 58800*lmMst1 + 4560*lmMt + 14400*pow2(lmMst1))*pow2(
        Mt) + 15*(655743 + 5288*lmMst1 - 11820*pow2(lmMst1))*pow2(Mst1)*pow2(
        s2t)) - 2*(15*Mst1*Mt*s2t*(-954181 - 11256*lmMst1 + 11520*pow2(lmMst1))
        + 2*(-6233611 - 58800*lmMst1 + 4560*lmMt + 14400*pow2(lmMst1))*pow2(Mt)
        + 15*(655743 + 5288*lmMst1 - 11820*pow2(lmMst1))*pow2(Mst1)*pow2(s2t))*
        pow3(Dmst12) + 2*Dmst12*Mt*(15*Mst1*s2t*(954181 + 11256*lmMst1 - 11520*
        pow2(lmMst1)) - 2*Mt*(-6233611 - 58800*lmMst1 + 4560*lmMt + 14400*pow2(
        lmMst1)))*pow4(Mst2) - 50*(-230009 - 4704*lmMst1 + 480*lmMt + 1152*
        pow2(lmMst1))*pow2(Mt)*pow6(Mst2)) + 40*Dmglst1*Mst1*(-(pow2(Dmst12)*
        pow2(Mst2)*(12000*Dmsqst1*Mst1*Mt*s2t + pow2(Msq)*(6*Mst1*Mt*s2t*(-
        949861 + 1944*lmMst1 + 11520*pow2(lmMst1)) + 10*(403559 + 384*(lmMst1 +
        lmMt) - 4608*pow2(lmMst1))*pow2(Mt) + 15*(-84209 - 1264*lmMst1 + 240*
        pow2(lmMst1))*pow2(Mst1)*pow2(s2t)))) + 2*(6000*Dmsqst1*Mst1*Mt*s2t +
        pow2(Msq)*(6*Mst1*Mt*s2t*(-1282471 + 7264*lmMst1 + 18120*pow2(lmMst1))
        + 5*(403559 + 384*(lmMst1 + lmMt) - 4608*pow2(lmMst1))*pow2(Mt) + 15*(-
        84209 - 1264*lmMst1 + 240*pow2(lmMst1))*pow2(Mst1)*pow2(s2t)))*pow3(
        Dmst12) + 10*Dmst12*Mt*(1200*Dmsqst1*Mst1*s2t + (Mt*(403559 + 384*(
        lmMst1 + lmMt) - 4608*pow2(lmMst1)) - 12*Mst1*s2t*(-33261 + 532*lmMst1
        + 660*pow2(lmMst1)))*pow2(Msq))*pow4(Mst2) + 240*(9631 + 16*lmMst1 +
        48*lmMt - 192*pow2(lmMst1))*pow2(Msq)*pow2(Mt)*pow6(Mst2)) + 15*pow2(
        Mst1)*(-(pow2(Dmst12)*pow2(Mst2)*(2400*Dmsqst1*Mst1*s2t*(56*Mt + Mst1*
        s2t*(7 - 12*lmMst1*(-2 + shiftst1 + shiftst2) + 18*(shiftst1 +
        shiftst2))) - pow2(Msq)*(-240*Mst1*Mt*s2t*(-10473 + 40*lmMst1 + 256*
        pow2(lmMst1)) + (852541 + 9216*lmMst1 - 10240*lmMt + 6144*pow2(lmMst1))
        *pow2(Mt) + 80*(1429 - 454*lmMst1 - 180*(shiftst1 + shiftst2) + 360*
        lmMst1*(shiftst1 + shiftst2) - 126*shiftst3 + 108*lmMst1*shiftst3 + 24*
        pow2(lmMst1))*pow2(Mst1)*pow2(s2t)))) + 2*(-2400*Dmsqst1*Mst1*s2t*(-28*
        Mt + Mst1*s2t*(-7 + 24*lmMst1*(-1 + shiftst2) - 36*shiftst2)) + pow2(
        Msq)*(80*Mst1*Mt*s2t*(-36863 + 80*lmMst1 + 552*pow2(lmMst1)) + (-
        1763661 + 47104*lmMst1 - 5120*lmMt + 24576*pow2(lmMst1))*pow2(Mt) + (
        350605 + 4320*shiftst1 + 2880*shiftst2 + 8352*shiftst3 - 96*lmMst1*(-
        115 + 90*shiftst1 + 60*shiftst2 + 54*shiftst3) - 2160*pow2(lmMst1))*
        pow2(Mst1)*pow2(s2t)))*pow3(Dmst12) + 160*Dmst12*(60*Dmsqst1*Mst1*s2t*(
        14*Mt + 3*(3 - 2*lmMst1)*Mst1*s2t*(shiftst1 - shiftst2)) + pow2(Msq)*(
        4*Mst1*Mt*s2t*(1361 + 10*lmMst1 + 54*pow2(lmMst1)) + (11389 - 704*
        lmMst1 + 192*lmMt - 384*pow2(lmMst1))*pow2(Mt) + 18*(1 - 2*lmMst1)*(10*
        shiftst1 - 10*shiftst2 + shiftst3)*pow2(Mst1)*pow2(s2t)))*pow4(Mst2) +
        1920*(349 - 56*lmMst1 + 24*lmMt - 32*pow2(lmMst1))*pow2(Msq)*pow2(Mt)*
        pow6(Mst2))))/2. - 2*Mst1*Tbeta*(225*Dmst12*Mst1*s2t*pow2(Mt)*(Dmst12*
        Mst1*s2t*(2*pow2(Dmglst1)*(954181 + 11256*lmMst1 - 11520*pow2(lmMst1))*
        pow2(Msq)*(2*Dmst12 - pow2(Mst2)) + 5*pow2(Mst1)*(1680*Dmsqst1*(2*
        Dmst12 - pow2(Mst2)) + pow2(Msq)*(Dmst12*(-20531 + 200*lmMst1 + 1200*
        pow2(lmMst1)) - 8*(1361 + 10*lmMst1 + 54*pow2(lmMst1))*pow2(Mst2))) +
        Dmglst1*Mst1*(2000*Dmsqst1*(2*Dmst12 - pow2(Mst2)) - pow2(Msq)*(Dmst12*
        (284641 + 8696*lmMst1 + 1680*pow2(lmMst1)) - 20*(-33261 + 532*lmMst1 +
        660*pow2(lmMst1))*pow2(Mst2)))) + Mt*(4*pow2(Dmglst1)*(360*Dmsqst1*(
        pow2(Dmst12) - Dmst12*pow2(Mst2) + pow4(Mst2)) + pow2(Msq)*(4*pow2(
        Dmst12)*(-274009 - 964*lmMst1 - 104*lmMt + 8310*pow2(lmMst1)) + Dmst12*
        (515917 + 6972*lmMst1 - 192*lmMt - 11520*pow2(lmMst1))*pow2(Mst2) + 2*(
        32101 - 5044*lmMst1 + 400*lmMt - 5100*pow2(lmMst1))*pow4(Mst2))) +
        pow2(Mst1)*(100*Dmsqst1*((346 + 288*lmMst1)*pow2(Dmst12) - 161*Dmst12*
        pow2(Mst2) - 24*(1 + 12*lmMst1)*pow4(Mst2)) + pow2(Msq)*(pow2(Dmst12)*(
        102747 + 10592*lmMst1 - 640*lmMt - 13888*pow2(lmMst1)) + 20*Dmst12*(-
        2071 - 296*lmMst1 + 96*lmMt + 600*pow2(lmMst1))*pow2(Mst2) - 160*(631 -
        lmMst1 + 36*lmMt + 21*pow2(lmMst1))*pow4(Mst2))) + 16*Dmglst1*Mst1*(
        4700*Dmsqst1*(pow2(Dmst12) - Dmst12*pow2(Mst2) + pow4(Mst2)) + pow2(
        Msq)*(pow2(Dmst12)*(-12383 - 4128*lmMst1 - 80*lmMt + 1260*pow2(lmMst1))
        + Dmst12*(17539 + 574*lmMst1 + 160*lmMt - 1920*pow2(lmMst1))*pow2(Mst2)
        + 5*(-4539 + 596*lmMst1 - 48*lmMt + 516*pow2(lmMst1))*pow4(Mst2))))) +
        2*pow4(Mt)*(-3*pow2(Dmst12)*(4*pow2(Dmglst1)*(84334067 + 120*lmMst1*(
        2843 - 120*lmMt) + 202200*lmMt - 828000*pow2(lmMst1) - 21600*pow2(lmMt)
        )*pow2(Msq) + 2*Dmglst1*Mst1*(400*Dmsqst1*(557 + 120*lmMst1 - 120*lmMt)
        + (-39474953 + lmMst1*(52780 - 87600*lmMt) - 3300*lmMt + 319200*pow2(
        lmMst1) + 14400*pow2(lmMt))*pow2(Msq)) + 5*(4000*Dmsqst1*(65 - 6*lmMst1
        + 6*lmMt) + (-2284899 + 49840*lmMt - 32*lmMst1*(-1793 + 555*lmMt) +
        87360*pow2(lmMst1) + 28800*pow2(lmMt))*pow2(Msq))*pow2(Mst1))*pow2(
        Mst2) + 2*(6*pow2(Dmglst1)*(84334067 + 120*lmMst1*(2843 - 120*lmMt) +
        202200*lmMt - 828000*pow2(lmMst1) - 21600*pow2(lmMt))*pow2(Msq) + 2*
        Dmglst1*Mst1*(600*Dmsqst1*(557 + 120*lmMst1 - 120*lmMt) + (-193364399 -
        90000*lmMt + 300*lmMst1*(3781 + 582*lmMt) + 2005200*pow2(lmMst1) +
        43200*pow2(lmMt))*pow2(Msq)) + 15*(2000*Dmsqst1*(65 - 6*lmMst1 + 6*
        lmMt) + (-3454599 + 16840*lmMt + 48*lmMst1*(262 + 405*lmMt) + 46560*
        pow2(lmMst1))*pow2(Msq))*pow2(Mst1))*pow3(Dmst12) + 40*Dmst12*((3*pow2(
        Dmglst1)*(84334067 + 120*lmMst1*(2843 - 120*lmMt) + 202200*lmMt -
        828000*pow2(lmMst1) - 21600*pow2(lmMt))*pow2(Msq))/10. + 2*Dmglst1*
        Mst1*(30*Dmsqst1*(557 + 120*lmMst1 - 120*lmMt) + (3746977 + 4005*lmMt -
        18*lmMst1*(2711 + 1215*lmMt) - 52380*pow2(lmMst1))*pow2(Msq)) + 75*(20*
        Dmsqst1*(65 - 6*lmMst1 + 6*lmMt) + (11697 + 448*lmMst1 + 330*lmMt -
        372*lmMst1*lmMt + 408*pow2(lmMst1) + 288*pow2(lmMt))*pow2(Msq))*pow2(
        Mst1))*pow4(Mst2) + 160*(pow2(Dmglst1)*(90*Dmsqst1*(-206 + 15*lmMst1 -
        15*lmMt) + (3030652 + 19305*lmMt - 6*lmMst1*(1183 + 645*lmMt) - 55530*
        pow2(lmMst1) - 5400*pow2(lmMt))*pow2(Msq)) + 75*Dmglst1*Mst1*(10*
        Dmsqst1*(26 + 3*lmMst1 - 3*lmMt) + (9961 + 66*lmMt - 10*lmMst1*(67 +
        24*lmMt) + 42*pow2(lmMst1) + 72*pow2(lmMt))*pow2(Msq)) + 150*(15*
        Dmsqst1*(20 - 3*lmMst1 + 3*lmMt) + (434 - 83*lmMst1 + 174*lmMt - 66*
        lmMst1*lmMt + 183*pow2(lmMst1) + 108*pow2(lmMt))*pow2(Msq))*pow2(Mst1))
        *pow6(Mst2)) + 750*pow3(Mst1)*(-36*shiftst2*(5*Mt*pow2(Dmst12)*(
        Dmsqst1*(3 - 2*lmMst1) + (1 - 2*lmMst1)*pow2(Msq))*pow2(Mst1)*(2*Dmst12
        - 3*pow2(Mst2))*pow3(s2t) + 120*s2t*pow3(Mt)*((1 - 2*lmMst1)*pow2(Msq)*
        (Dmst12 + pow2(Mst2))*pow4(Mst2) + Dmsqst1*(-3 + 2*lmMst1)*(pow3(
        Dmst12) - Dmst12*pow4(Mst2) - pow6(Mst2)))) + Mt*s2t*(-(((pow2(Dmglst1)
        *(655743 + 5288*lmMst1 - 11820*pow2(lmMst1))*pow2(Msq))/5. + (Dmglst1*
        Mst1*(84209 + 1264*lmMst1 - 240*pow2(lmMst1))*pow2(Msq))/2. - 30*
        Dmsqst1*(7 + 24*lmMst1)*pow2(Mst1) + (1429 - 454*lmMst1 + 24*pow2(
        lmMst1))*pow2(Msq)*pow2(Mst1))*pow2(s2t)*pow3(Dmst12)) + 180*shiftst1*(
        Dmsqst1*(3 - 2*lmMst1) + (1 - 2*lmMst1)*pow2(Msq))*(pow2(Dmst12)*pow2(
        Mst1)*(4*Dmst12 - 3*pow2(Mst2))*pow2(s2t) + 24*pow2(Mt)*pow6(Mst2)) -
        18*shiftst3*pow2(Msq)*(pow2(Dmst12)*pow2(Mst1)*(Dmst12*(-13 + 14*
        lmMst1) + 3*(1 - 2*lmMst1)*pow2(Mst2))*pow2(s2t) + 8*pow2(Mt)*(6*(-2 +
        lmMst1)*pow2(Dmst12)*pow2(Mst2) + (14 - 6*lmMst1)*pow3(Dmst12) + 3*
        Dmst12*(3 - 2*lmMst1)*pow4(Mst2) + 3*(-1 + 2*lmMst1)*pow6(Mst2))))))) +
        75*(-576*z2*pow5(Mst1)*(-50*Mt*s2t*xDmsqst1*pow2(Dmsqst1)*(s2t*pow2(
        Dmst12)*(2*Mt*MuSUSY*(shiftst1 + shiftst2) + 3*s2t*(-shiftst1 +
        shiftst2)*Tbeta*pow2(Mst1))*pow2(Mst2) + (-8*Mt*MuSUSY*s2t*shiftst2 +
        24*shiftst2*Tbeta*pow2(Mt) + 2*(2*shiftst1 - shiftst2)*Tbeta*pow2(Mst1)
        *pow2(s2t))*pow3(Dmst12) - 4*Dmst12*Mt*(MuSUSY*s2t*(shiftst1 -
        shiftst2) + 6*Mt*shiftst2*Tbeta)*pow4(Mst2) + 24*(shiftst1 - shiftst2)*
        Tbeta*pow2(Mt)*pow6(Mst2)) + pow2(Msq)*(2*Dmst12*MuSUSY*pow2(Mt)*pow2(
        s2t)*(2*pow2(Dmst12)*(100*Dmsqst1*shiftst2 + (15*shiftst1 + 10*shiftst2
        + 9*shiftst3)*pow2(Msq)) - 5*Dmst12*(10*Dmsqst1*(shiftst1 + shiftst2) +
        (10*(shiftst1 + shiftst2) + 3*shiftst3)*pow2(Msq))*pow2(Mst2) + 10*(10*
        Dmsqst1*(shiftst1 - shiftst2) + (10*shiftst1 - 10*shiftst2 + shiftst3)*
        pow2(Msq))*pow4(Mst2)) + 5*shiftst3*Tbeta*pow2(Msq)*(-(Mt*pow2(Dmst12)*
        pow2(Mst1)*(7*Dmst12 - 3*pow2(Mst2))*pow3(s2t)) + 24*s2t*pow3(Mt)*(-(
        pow2(Dmst12)*pow2(Mst2)) + pow3(Dmst12) + Dmst12*pow4(Mst2) - pow6(
        Mst2))) - 50*Tbeta*(shiftst2*(-(Mt*pow2(Dmst12)*(Dmsqst1 + pow2(Msq))*
        pow2(Mst1)*(2*Dmst12 - 3*pow2(Mst2))*pow3(s2t)) + 24*s2t*pow3(Mt)*(
        Dmsqst1*pow3(Dmst12) - (Dmsqst1 + pow2(Msq))*(Dmst12 + pow2(Mst2))*
        pow4(Mst2))) + shiftst1*(Dmsqst1 + pow2(Msq))*(Mt*pow2(Dmst12)*pow2(
        Mst1)*(4*Dmst12 - 3*pow2(Mst2))*pow3(s2t) + 24*s2t*pow3(Mt)*pow6(Mst2))
        ))) + Mt*((32*Mst1*xDmsqst1*pow2(Dmsqst1)*(-20*Dmglst1*Mt*pow2(Mst1)*(-
        (pow2(Dmst12)*pow2(Mst2)*(-100*Mt*s2t*(5*MuSUSY - 141*Mst1*Tbeta) + 4*(
        557 + 120*lmMst1 - 120*lmMt)*Tbeta*pow2(Mt) + 375*Tbeta*pow2(Mst1)*
        pow2(s2t))) + 2*(-50*Mt*s2t*(5*MuSUSY - 141*Mst1*Tbeta) + 2*(557 + 120*
        lmMst1 - 120*lmMt)*Tbeta*pow2(Mt) + 375*Tbeta*pow2(Mst1)*pow2(s2t))*
        pow3(Dmst12) + 4*Dmst12*Mt*(-125*MuSUSY*s2t + (557 + 120*lmMst1 - 120*
        lmMt)*Mt*Tbeta + 3525*Mst1*s2t*Tbeta)*pow4(Mst2) + 100*(44 + 3*lmMst1 -
        3*lmMt)*Tbeta*pow2(Mt)*pow6(Mst2)) - 125*pow3(Mst1)*(-2*pow2(Dmst12)*
        pow2(Mst2)*(-6*s2t*(28*MuSUSY + Mst1*(11 + 36*shiftst2)*Tbeta)*pow2(Mt)
        - 3*Mst1*Mt*(MuSUSY*(7 + 30*shiftst1 + 6*shiftst2 - 12*lmMst1*(-2 +
        shiftst1 + shiftst2)) - 42*Mst1*Tbeta)*pow2(s2t) + Tbeta*(8*(65 - 6*
        lmMst1 + 6*lmMt)*pow3(Mt) - 54*(-2 + lmMst1)*(shiftst1 - shiftst2)*
        pow3(Mst1)*pow3(s2t))) + pow3(Dmst12)*(-(s2t*(336*MuSUSY - 9*Mst1*(47 +
        96*(lmMst1 + shiftst2 - lmMst1*shiftst2))*Tbeta)*pow2(Mt)) - 12*Mst1*
        Mt*(MuSUSY*(7 + 6*shiftst1 - 24*lmMst1*(-1 + shiftst2) + 30*shiftst2) -
        42*Mst1*Tbeta)*pow2(s2t) + Tbeta*(16*(65 - 6*lmMst1 + 6*lmMt)*pow3(Mt)
        + 3*(7 + 108*shiftst1 - 72*shiftst2 + 24*lmMst1*(1 - 2*shiftst1 +
        shiftst2))*pow3(Mst1)*pow3(s2t))) - Dmst12*Mt*(3*Mt*s2t*(112*MuSUSY +
        Mst1*(229 - 288*lmMst1*(-1 + shiftst2) + 576*shiftst2)*Tbeta) + 16*(-65
        + 6*lmMst1 - 6*lmMt)*Tbeta*pow2(Mt) - 144*(-2 + lmMst1)*Mst1*MuSUSY*(
        shiftst1 - shiftst2)*pow2(s2t))*pow4(Mst2) + 16*((91 - 12*lmMst1 + 12*
        lmMt)*Mt - 54*(-2 + lmMst1)*Mst1*s2t*(shiftst1 - shiftst2))*Tbeta*pow2(
        Mt)*pow6(Mst2)) + 8*Tbeta*pow2(Dmglst1)*pow2(Mt)*(4*Dmglst1*(-1541 +
        420*lmMst1 - 420*lmMt)*Mt*xDmglst1*pow6(Mst2) + 15*Mst1*(-45*Dmst12*
        Mst1*s2t*(pow2(Dmst12) - Dmst12*pow2(Mst2) + pow4(Mst2)) + 4*(206 - 15*
        lmMst1 + 15*lmMt)*Mt*pow6(Mst2)))))/5. - 15*z3*(-5040*Dmst12*s2t*
        xDmsqst1*pow2(Dmsqst1)*(2*Dmst12*Mt*MuSUSY*s2t*pow2(Mst2) + pow2(
        Dmst12)*(-4*Mt*MuSUSY*s2t + 15*Tbeta*pow2(Mt) + Tbeta*pow2(Mst1)*pow2(
        s2t)) - 15*Tbeta*pow2(Mt)*pow4(Mst2))*pow5(Mst1) + 32*Mt*xDmglst1*pow3(
        Dmglst1)*pow4(Msq)*(44*Mt*MuSUSY*(4789*Dmst12*Mst1*s2t*(pow2(Dmst12) -
        Dmst12*pow2(Mst2) + pow4(Mst2)) + 1503*Mt*pow6(Mst2)) - 3*Mst1*Tbeta*(
        52679*pow2(Dmst12)*pow2(Mst1)*(2*Dmst12 - pow2(Mst2))*pow2(s2t) +
        38236*Dmst12*Mst1*Mt*s2t*(pow2(Dmst12) - Dmst12*pow2(Mst2) + pow4(Mst2)
        ) + 3*pow2(Mt)*(-71189*pow2(Dmst12)*pow2(Mst2) + 71189*Dmst12*(pow2(
        Dmst12) + pow4(Mst2)) + 35412*pow6(Mst2)))) + Mst1*pow2(Msq)*(Mt*
        MuSUSY*(-2*pow2(Dmst12)*pow2(Mst1)*(20*Dmglst1*(29926*Dmglst1 + 10097*
        Mst1)*pow2(Msq)*(2*Dmst12 - pow2(Mst2)) - pow2(Mst1)*(5040*Dmsqst1*(2*
        Dmst12 - pow2(Mst2)) + pow2(Msq)*(37669*Dmst12 + 10296*pow2(Mst2))))*
        pow2(s2t) + Mt*pow2(Msq)*(16*Dmst12*Mst1*s2t*(210810*pow2(Dmglst1)*(
        pow2(Dmst12) - Dmst12*pow2(Mst2) + pow4(Mst2)) + pow2(Mst1)*(-40594*
        pow2(Dmst12) + 17285*Dmst12*pow2(Mst2) + 6024*pow4(Mst2)) + Dmglst1*
        Mst1*(-285974*pow2(Dmst12) + 105405*Dmst12*pow2(Mst2) + 75164*pow4(
        Mst2))) + 3*Mt*(-(pow2(Mst1)*(-31963*pow2(Dmst12)*pow2(Mst2) + 131926*
        pow3(Dmst12) - 68000*Dmst12*pow4(Mst2) - 24064*pow6(Mst2))) + 48*
        Dmglst1*Mst1*(-8213*pow2(Dmst12)*pow2(Mst2) + 8213*Dmst12*(pow2(Dmst12)
        + pow4(Mst2)) + 4664*pow6(Mst2)) + 96*pow2(Dmglst1)*(-10226*pow2(
        Dmst12)*pow2(Mst2) + 10226*Dmst12*(pow2(Dmst12) + pow4(Mst2)) + 4715*
        pow6(Mst2))))) - 4*Mst1*Tbeta*(-((5*Dmglst1*(29926*Dmglst1 + 10097*
        Mst1)*pow2(Msq) - 18*(70*Dmsqst1 - 143*pow2(Msq))*pow2(Mst1))*pow3(
        Dmst12)*pow3(Mst1)*pow3(s2t)) - Dmst12*Mst1*s2t*pow2(Mt)*(-(pow2(Mst1)*
        (pow2(Msq)*(32783*pow2(Dmst12) - 14508*Dmst12*pow2(Mst2) - 15840*pow4(
        Mst2)) + 1260*Dmsqst1*(22*pow2(Dmst12) - 7*Dmst12*pow2(Mst2) - 8*pow4(
        Mst2)))) + 24*Dmglst1*pow2(Msq)*(11*Dmglst1*(5777*pow2(Dmst12) - 2607*
        Dmst12*pow2(Mst2) - 563*pow4(Mst2)) + Mst1*(4744*pow2(Dmst12) - 3821*
        Dmst12*pow2(Mst2) + 2898*pow4(Mst2)))) + pow2(Msq)*(3*Mt*pow2(Dmst12)*
        pow2(Mst1)*(210810*pow2(Dmglst1)*(2*Dmst12 - pow2(Mst2)) - pow2(Mst1)*(
        11261*Dmst12 + 6024*pow2(Mst2)) - Dmglst1*Mst1*(30241*Dmst12 + 75164*
        pow2(Mst2)))*pow2(s2t) + 2*pow3(Mt)*(1404*pow2(Dmglst1)*(-1065*pow2(
        Dmst12)*pow2(Mst2) + 1065*pow3(Dmst12) + 1065*Dmst12*pow4(Mst2) + 512*
        pow6(Mst2)) + pow2(Mst1)*(51181*pow2(Dmst12)*pow2(Mst2) - 152018*pow3(
        Dmst12) + 49656*Dmst12*pow4(Mst2) + 10176*pow6(Mst2)) + 4*Dmglst1*Mst1*
        (86833*pow2(Dmst12)*pow2(Mst2) - 287078*pow3(Dmst12) + 113412*Dmst12*
        pow4(Mst2) + 47112*pow6(Mst2))))))))))))/(3.888e6*Tbeta*pow4(Msq)*pow5(
        Mst1)*pow6(Mst2));
}

/**
 * @return returns the susy log^0 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double H3q22g::calc_coef_at_as2_no_sm_logs_log0() const {

   const double result =
      (8*pow3(Dmglst1)*(Mt*pow2(Dmst12)*pow2(Mst2)*(31305120*Mst1*Mt*s2t*pow2(
        Dmsqst1) + 31305120*Dmsqst1*Mst1*Mt*s2t*pow2(Msq) + (Mst1*Mt*s2t*(-
        21284082326 + 17749864125*z3) + (673066160 - 615195000*z3)*pow2(Mt) +
        7350*(-520877 + 430155*z3)*pow2(Mst1)*pow2(s2t))*pow4(Msq)) + pow3(
        Dmst12)*(-31305120*Mst1*s2t*pow2(Dmsqst1)*pow2(Mt) - 31305120*Dmsqst1*
        Mst1*s2t*pow2(Msq)*pow2(Mt) + (2*Mst1*s2t*(31902674758 - 26534253375*
        z3)*pow2(Mt) - 14700*Mt*(-520877 + 430155*z3)*pow2(Mst1)*pow2(s2t) +
        40*(-16826654 + 15379875*z3)*pow3(Mt) + 245*(14217821 - 11852775*z3)*
        pow3(Mst1)*pow3(s2t))*pow4(Msq)) + 4*Dmst12*pow2(Mt)*(-7826280*Mst1*
        s2t*pow2(Dmsqst1) - 7826280*Dmsqst1*Mst1*s2t*pow2(Msq) + (10*Mt*(-
        16826654 + 15379875*z3) + 49*Mst1*s2t*(-108352984 + 89636625*z3))*pow4(
        Msq))*pow4(Mst2) - 392*pow3(Mt)*(154740*pow2(Dmsqst1) + 154740*Dmsqst1*
        pow2(Msq) + (10583177 - 8913375*z3)*pow4(Msq))*pow6(Mst2)) + 49*pow3(
        Mst1)*(-150*Dmsqst1*pow2(Msq)*(-4*Mt*pow2(Dmst12)*pow2(Mst2)*(19600*
        Mst1*Mt*s2t + (-9122 + 14175*z3)*pow2(Mt) + 450*(2 - 21*z3)*pow2(Mst1)*
        pow2(s2t)) + pow3(Dmst12)*(-20800*Mst1*s2t*pow2(Mt) - 75*Mt*(274 + 63*
        z3)*pow2(Mst1)*pow2(s2t) + 16*(-6536 + 14175*z3)*pow3(Mt) - 8400*pow3(
        Mst1)*pow3(s2t)) - 200*Dmst12*(-888*Mst1*s2t + Mt*(-158 + 567*z3))*
        pow2(Mt)*pow4(Mst2) - 201600*pow3(Mt)*pow6(Mst2)) - 150*pow2(Dmsqst1)*(
        -(Mt*pow2(Dmst12)*pow2(Mst2)*(36800*Mst1*Mt*s2t + 20112*pow2(Mt) + 75*(
        458 - 945*z3)*pow2(Mst1)*pow2(s2t))) + 3*pow3(Dmst12)*(-20800*Mst1*s2t*
        pow2(Mt) - 525*Mt*(-26 + 45*z3)*pow2(Mst1)*pow2(s2t) + 4*(-3998 +
        14175*z3)*pow3(Mt) - 2800*pow3(Mst1)*pow3(s2t)) - 100*Dmst12*(-1360*
        Mst1*s2t + 63*Mt*(-14 + 27*z3))*pow2(Mt)*pow4(Mst2) - 200*(410 + 567*
        z3)*pow3(Mt)*pow6(Mst2)) - pow4(Msq)*(-20*Mt*pow2(Dmst12)*pow2(Mst2)*(
        300*Mst1*Mt*s2t*(-17512 + 14805*z3) + (-375892 + 621675*z3)*pow2(Mt) -
        900*(-1226 + 495*z3)*pow2(Mst1)*pow2(s2t)) + pow3(Dmst12)*(-30*Mst1*
        s2t*(-1117238 + 877575*z3)*pow2(Mt) - 2250*Mt*(-5570 + 333*z3)*pow2(
        Mst1)*pow2(s2t) + (-41715182 + 38174625*z3)*pow3(Mt) + 3000*(-2722 +
        2259*z3)*pow3(Mst1)*pow3(s2t)) - 6000*Dmst12*(81*Mt*(-406 + 285*z3) +
        8*Mst1*s2t*(-694 + 477*z3))*pow2(Mt)*pow4(Mst2) + 48000*(623 + 963*z3)*
        pow3(Mt)*pow6(Mst2))) + 196*Dmglst1*pow2(Mst1)*(600*pow2(Dmsqst1)*(-6*
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
        )) + 2*Mst1*pow2(Dmglst1)*(35280*Dmsqst1*Mt*pow2(Msq)*(-(pow2(Dmst12)*
        pow2(Mst2)*(7940*Mst1*Mt*s2t + 914*pow2(Mt) + 225*pow2(Mst1)*pow2(s2t))
        ) + (7940*Mst1*Mt*s2t + 914*pow2(Mt) + 450*pow2(Mst1)*pow2(s2t))*pow3(
        Dmst12) + 2*Dmst12*Mt*(457*Mt + 3970*Mst1*s2t)*pow4(Mst2) + 4760*pow2(
        Mt)*pow6(Mst2)) + 35280*Mt*pow2(Dmsqst1)*(-(pow2(Dmst12)*pow2(Mst2)*(
        7940*Mst1*Mt*s2t + 914*pow2(Mt) + 225*pow2(Mst1)*pow2(s2t))) + (7940*
        Mst1*Mt*s2t + 914*pow2(Mt) + 450*pow2(Mst1)*pow2(s2t))*pow3(Dmst12) +
        2*Dmst12*Mt*(457*Mt + 3970*Mst1*s2t)*pow4(Mst2) + 5022*pow2(Mt)*pow6(
        Mst2)) + pow4(Msq)*(Mt*pow2(Dmst12)*pow2(Mst2)*(196*Mst1*Mt*s2t*(-
        263717842 + 218365875*z3) + (27129768542 - 22522586625*z3)*pow2(Mt) +
        22050*(-63802 + 92895*z3)*pow2(Mst1)*pow2(s2t)) + pow3(Dmst12)*(-392*
        Mst1*s2t*(-384557822 + 319453875*z3)*pow2(Mt) + 66150*Mt*(-150737 +
        112420*z3)*pow2(Mst1)*pow2(s2t) + (-25287306692 + 22556819250*z3)*pow3(
        Mt) + 3675*(1908362 - 1581075*z3)*pow3(Mst1)*pow3(s2t)) + 392*Dmst12*(
        20*Mst1*s2t*(-6041999 + 5054400*z3) + Mt*(-73908751 + 57368250*z3))*
        pow2(Mt)*pow4(Mst2) + 3920*(-8223692 + 6125625*z3)*pow3(Mt)*pow6(Mst2))
        ))/(1.9845e6*pow3(Mst1)*pow3(Mt)*pow4(Msq)*pow6(Mst2));

   return result;
}

/**
 * @return returns the susy log^1 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double H3q22g::calc_coef_at_as2_no_sm_logs_log1() const {

   const double result =
      (4*(2*pow3(Dmglst1)*(-2*Mt*pow2(Dmst12)*pow2(Mst2)*(82320*Mst1*Mt*s2t*
        pow2(Dmsqst1) + 82320*Dmsqst1*Mst1*Mt*s2t*pow2(Msq) + (346639*Mst1*Mt*
        s2t + 555463*pow2(Mt) + 1051785*pow2(Mst1)*pow2(s2t))*pow4(Msq)) +
        pow3(Dmst12)*(164640*Mst1*s2t*pow2(Dmsqst1)*pow2(Mt) + 164640*Dmsqst1*
        Mst1*s2t*pow2(Msq)*pow2(Mt) + (8778304*Mst1*s2t*pow2(Mt) + 4207140*Mt*
        pow2(Mst1)*pow2(s2t) + 1110926*pow3(Mt) + 661255*pow3(Mst1)*pow3(s2t))*
        pow4(Msq)) + 2*Dmst12*pow2(Mt)*(82320*Mst1*s2t*pow2(Dmsqst1) + 82320*
        Dmsqst1*Mst1*s2t*pow2(Msq) + (555463*Mt - 3695874*Mst1*s2t)*pow4(Msq))*
        pow4(Mst2) - 4704*pow3(Mt)*(10*pow2(Dmsqst1) + 10*Dmsqst1*pow2(Msq) +
        141*pow4(Msq))*pow6(Mst2)) + Mst1*pow2(Dmglst1)*(17640*pow2(Dmsqst1)*
        pow2(Mt)*((-17*Mt + 10*Mst1*s2t)*pow2(Dmst12)*pow2(Mst2) + (17*Mt - 10*
        Mst1*s2t)*pow3(Dmst12) + Dmst12*(17*Mt - 10*Mst1*s2t)*pow4(Mst2) + 21*
        Mt*pow6(Mst2)) + 17640*Dmsqst1*pow2(Msq)*pow2(Mt)*((-17*Mt + 10*Mst1*
        s2t)*pow2(Dmst12)*pow2(Mst2) + (17*Mt - 10*Mst1*s2t)*pow3(Dmst12) +
        Dmst12*(17*Mt - 10*Mst1*s2t)*pow4(Mst2) + 35*Mt*pow6(Mst2)) + pow4(Msq)
        *(Mt*pow2(Dmst12)*pow2(Mst2)*(-4088560*Mst1*Mt*s2t + 968629*pow2(Mt) +
        1853670*pow2(Mst1)*pow2(s2t)) + pow3(Dmst12)*(7502488*Mst1*s2t*pow2(Mt)
        - 3134775*Mt*pow2(Mst1)*pow2(s2t) + 2019394*pow3(Mt) + 689430*pow3(
        Mst1)*pow3(s2t)) - 196*Dmst12*(20187*Mt - 3442*Mst1*s2t)*pow2(Mt)*pow4(
        Mst2) - 8199072*pow3(Mt)*pow6(Mst2))) - 98*Dmglst1*pow2(Mst1)*(600*
        Dmsqst1*pow2(Msq)*pow2(Mt)*(-((6*Mt + Mst1*s2t)*pow2(Dmst12)*pow2(Mst2)
        ) + (6*Mt - 3*Mst1*s2t)*pow3(Dmst12) + Dmst12*(6*Mt + 5*Mst1*s2t)*pow4(
        Mst2) + 25*Mt*pow6(Mst2)) + 300*pow2(Dmsqst1)*pow2(Mt)*(3*(-4*Mt +
        Mst1*s2t)*pow2(Dmst12)*pow2(Mst2) + (12*Mt - 11*Mst1*s2t)*pow3(Dmst12)
        + Dmst12*(12*Mt + 5*Mst1*s2t)*pow4(Mst2) + 35*Mt*pow6(Mst2)) + pow4(
        Msq)*(Mt*pow2(Dmst12)*pow2(Mst2)*(29758*Mst1*Mt*s2t + 6677*pow2(Mt) +
        22350*pow2(Mst1)*pow2(s2t)) + pow3(Dmst12)*(-34587*Mst1*s2t*pow2(Mt) -
        18045*Mt*pow2(Mst1)*pow2(s2t) + 22414*pow3(Mt) + 3325*pow3(Mst1)*pow3(
        s2t)) - 8*Dmst12*(4471*Mt + 6875*Mst1*s2t)*pow2(Mt)*pow4(Mst2) + 50800*
        pow3(Mt)*pow6(Mst2))) - 49*pow3(Mst1)*(-150*Dmsqst1*Mt*pow2(Msq)*(pow2(
        Dmst12)*pow2(Mst2)*(-80*Mst1*Mt*s2t - 17*pow2(Mt) + 180*pow2(Mst1)*
        pow2(s2t)) + 2*(20*Mst1*Mt*s2t + 187*pow2(Mt) - 90*pow2(Mst1)*pow2(s2t)
        )*pow3(Dmst12) - 20*Dmst12*Mt*(17*Mt - 6*Mst1*s2t)*pow4(Mst2) - 1080*
        pow2(Mt)*pow6(Mst2)) - 150*Mt*pow2(Dmsqst1)*(-4*pow2(Dmst12)*pow2(Mst2)
        *(10*Mst1*Mt*s2t + 3*pow2(Mt) - 45*pow2(Mst1)*pow2(s2t)) + 9*(41*pow2(
        Mt) - 20*pow2(Mst1)*pow2(s2t))*pow3(Dmst12) + 5*Dmst12*Mt*(-69*Mt + 16*
        Mst1*s2t)*pow4(Mst2) - 1010*pow2(Mt)*pow6(Mst2)) - pow4(Msq)*(-2*Mt*
        pow2(Dmst12)*pow2(Mst2)*(25850*Mst1*Mt*s2t + 21033*pow2(Mt) + 75*pow2(
        Mst1)*pow2(s2t)) + pow3(Dmst12)*(68264*Mst1*s2t*pow2(Mt) + 5700*Mt*
        pow2(Mst1)*pow2(s2t) + 36107*pow3(Mt) + 250*pow3(Mst1)*pow3(s2t)) +
        200*Dmst12*(131*Mt + 100*Mst1*s2t)*pow2(Mt)*pow4(Mst2) + 400*(-533 +
        54*z3)*pow3(Mt)*pow6(Mst2)))))/(33075.*pow3(Mst1)*pow3(Mt)*pow4(Msq)*
        pow6(Mst2));

   return result;
}

/**
 * @return returns the susy log^2 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double H3q22g::calc_coef_at_as2_no_sm_logs_log2() const {

   const double result =
      (8*(-7*Dmglst1*pow2(Mst1)*pow4(Msq)*(2*Mt*pow2(Dmst12)*pow2(Mst2)*(-1304*
        Mst1*Mt*s2t + 888*pow2(Mt) + 645*pow2(Mst1)*pow2(s2t)) + pow3(Dmst12)*(
        1544*Mst1*s2t*pow2(Mt) - 2250*Mt*pow2(Mst1)*pow2(s2t) - 1640*pow3(Mt) +
        275*pow3(Mst1)*pow3(s2t)) - 8*Dmst12*(239*Mt - 35*Mst1*s2t)*pow2(Mt)*
        pow4(Mst2) - 6520*pow3(Mt)*pow6(Mst2)) - 8*pow3(Dmglst1)*pow4(Msq)*(-3*
        Mt*pow2(Dmst12)*pow2(Mst2)*(-75*Mst1*Mt*s2t + 1122*pow2(Mt) + 560*pow2(
        Mst1)*pow2(s2t)) + pow3(Dmst12)*(1671*Mst1*s2t*pow2(Mt) + 3360*Mt*pow2(
        Mst1)*pow2(s2t) + 3366*pow3(Mt) + 140*pow3(Mst1)*pow3(s2t)) + 3*Dmst12*
        (1122*Mt - 707*Mst1*s2t)*pow2(Mt)*pow4(Mst2) + 2163*pow3(Mt)*pow6(Mst2)
        ) + Mst1*pow2(Dmglst1)*pow4(Msq)*(Mt*pow2(Dmst12)*pow2(Mst2)*(4088*
        Mst1*Mt*s2t - 22884*pow2(Mt) + 8925*pow2(Mst1)*pow2(s2t)) + pow3(
        Dmst12)*(-42728*Mst1*s2t*pow2(Mt) + 1155*Mt*pow2(Mst1)*pow2(s2t) +
        56772*pow3(Mt) - 3360*pow3(Mst1)*pow3(s2t)) - 28*Dmst12*(393*Mt - 1234*
        Mst1*s2t)*pow2(Mt)*pow4(Mst2) + 3948*pow3(Mt)*pow6(Mst2)) + 7*pow3(
        Mst1)*(Mt*pow2(Dmst12)*pow2(Mst2)*(1760*Mst1*Mt*s2t + 538*pow2(Mt) +
        105*pow2(Mst1)*pow2(s2t))*pow4(Msq) - pow3(Dmst12)*(1032*Mst1*s2t*pow2(
        Mt) + 480*Mt*pow2(Mst1)*pow2(s2t) + 644*pow3(Mt) - 45*pow3(Mst1)*pow3(
        s2t))*pow4(Msq) - 40*Dmst12*(11*Mt + 61*Mst1*s2t)*pow2(Mt)*pow4(Msq)*
        pow4(Mst2) + 10*pow3(Mt)*(45*pow2(Dmsqst1) + 90*Dmsqst1*pow2(Msq) +
        442*pow4(Msq))*pow6(Mst2))))/(315.*pow3(Mst1)*pow3(Mt)*pow4(Msq)*pow6(
        Mst2));

   return result;
}

/**
 * @return returns the susy log^3 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double H3q22g::calc_coef_at_as2_no_sm_logs_log3() const {

   const double result =
      -298.6666666666667;

   return result;
}

} // namespace hierarchies
} // namespace himalaya
