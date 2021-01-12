// ====================================================================
// This file is part of Himalaya.
//
// Himalaya is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================

#include "H6b.hpp"
#include "enums.hpp"
#include "Constants.hpp"
#include "power.hpp"
#include <cmath>

namespace himalaya{
namespace hierarchies{

/**
 * Constructor
 * @param expansionDepth the flagMap for the truncation of expansion variables
 * @param Al4p a double alpha_s/4/Pi
 * @param beta a double which is the mixing angle beta
 * @param Dmglst2 a double Mgl - Mst2
 * @param Dmsqst2 a double Msq - Mst2
 * @param lmMt a double log((<renormalization scale> / Mt)^2)
 * @param lmMst1 a double log((<renormalization scale> / Mst1)^2)
 * @param lmMst2 a double log((<renormalization scale> / Mst2)^2)
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
H6b::H6b(const ExpansionFlags_t& expansionDepth, double Al4p, double beta, double Dmglst2,
                 double Dmsqst2, double lmMt, double lmMst1, double lmMst2,
                 double Mt, double Mst1, double Mst2, double MuSUSY,
                 double s2t,
                 int mdrFlag, int oneLoopFlag, int twoLoopFlag, int threeLoopFlag){
   // abbrev for tan(beta) and sin(beta)
   Tbeta = tan(beta);
   Sbeta = sin(beta);
   this -> Dmglst2 = Dmglst2;
   this -> Dmsqst2 = Dmsqst2;
   this -> lmMt = lmMt;
   this -> lmMst1 = lmMst1;
   this -> lmMst2 = lmMst2;
   this -> Mt = Mt;
   this -> Mst1 = Mst1;
   this -> Mst2 = Mst2;
   this -> MuSUSY = MuSUSY;
   this -> s2t = s2t;
   this -> oneLoopFlag = oneLoopFlag;
   this -> twoLoopFlag = twoLoopFlag;
   this -> threeLoopFlag = threeLoopFlag;
   this -> Al4p = Al4p;
   // mdr flags, indicates if one wants to shift the dr stop mass to the mdr stop mass
   shiftst1 = mdrFlag;
   shiftst2 = mdrFlag;
   shiftst3 = mdrFlag;
   xDR2DRMOD = mdrFlag;
   // expansion flags
   xDmglst2 = expansionDepth.at(ExpansionDepth::Dmglst2);
   xDmsqst2 = expansionDepth.at(ExpansionDepth::Dmsqst2);
   xMst = expansionDepth.at(ExpansionDepth::Mst);
}

/**
 * @return The diagonal (1, 1) matrix element of the Higgs mass matrix as a double for the hierarchy 'H6b'
 */
double H6b::getS1() const {
   return (Tbeta*(324*z2*pow2(Mst1)*(6*Al4p*pow2(Mt)*pow2(MuSUSY)*(-6*s2t*
        twoLoopFlag*pow2(Mst1)*(4*Mt*pow3(Mst2)*((9*Dmglst2 + Mst2)*pow2(Mst1)*
        pow2(Mst2) + (13*Dmglst2 + Mst2)*pow4(Mst1) + 5*Dmglst2*pow4(Mst2) +
        pow5(Mst2)) - s2t*pow4(Mst2)*(4*Dmglst2*(pow2(Mst1)*pow2(Mst2) + pow4(
        Mst1) + pow4(Mst2)) + pow5(Mst2)) + 2*xMst*(-2*Dmglst2*Mst2*(-17*Mt +
        Mst2*s2t) + (48*Mt - Mst2*s2t)*xDmglst2*pow2(Dmglst2) + 2*Mt*pow2(Mst2)
        )*pow6(Mst1)) + (xDmglst2*pow2(Dmglst2)*pow2(Mst2)*(3240*s2t*
        twoLoopFlag*pow2(Mst1)*((-24*Mt + Mst2*s2t)*pow2(Mst1)*pow2(Mst2) + (-
        36*Mt + Mst2*s2t)*pow4(Mst1) + (-12*Mt + Mst2*s2t)*pow4(Mst2)) + Al4p*
        Mst2*threeLoopFlag*(-3*pow2(Mst1)*pow2(Mst2)*(8*(46987 + 15390*lmMst1 +
        4050*lmMst2)*Mst2*Mt*s2t + 192*(17 + 1320*lmMst2)*pow2(Mt) + (-40001 -
        52560*lmMst1 + 37080*lmMst2)*pow2(Mst2)*pow2(s2t)) - 2*(5*(-184517 +
        74952*lmMst1 + 18360*lmMst2)*Mst2*Mt*s2t + 6*(12761 + 5400*lmMst1 +
        178920*lmMst2)*pow2(Mt) + (261247 - 163890*lmMst1 + 140670*lmMst2)*
        pow2(Mst2)*pow2(s2t))*pow4(Mst1) + 28080*pow2(s2t)*pow6(Mst2))))/270.)
        + (Mt*threeLoopFlag*pow2(Al4p)*pow3(Mst2)*(450*s2t*xDmsqst2*pow2(
        Dmsqst2)*(Mt*pow2(MuSUSY)*(6*Mst2*(-144*(3 + lmMst1 - lmMst2)*Mt + (209
        + 18*lmMst1 - 18*lmMst2)*Mst2*s2t)*pow2(Mst1) + 8*(421 + 54*lmMst1 -
        54*lmMst2)*s2t*pow4(Mst1) + 108*s2t*pow4(Mst2)) + 27*s2t*((51 - 6*
        lmMst1 + 6*lmMst2)*Mt + 2*(11 + 10*lmMst1 - 10*lmMst2)*Mst2*s2t)*(-1 +
        pow2(Sbeta))*pow2(Sbeta)*pow6(Mst1)) + Mt*(2*pow2(MuSUSY)*(27*pow2(
        Mst1)*pow2(Mst2)*(960*Mst2*s2t*(8*(2 - 3*lmMst2)*Mt + (1 - 10*lmMst1 +
        12*lmMst2)*Mst2*s2t)*xDmglst2*xDR2DRMOD*pow2(Dmglst2) - Dmglst2*Mt*s2t*
        (-14400*Dmsqst2*(2 + 3*lmMst1 - 3*lmMst2) + (173947 - 25080*lmMst1 +
        7680*xDR2DRMOD + 360*lmMst2*(191 + 64*xDR2DRMOD))*pow2(Mst2)) - 7680*
        Dmglst2*(5 + 6*lmMst2)*Mst2*pow2(Mt) + 240*Dmglst2*(7 - 14*lmMst1 + 30*
        lmMst2)*xDR2DRMOD*pow2(s2t)*pow3(Mst2) - 5*Mst2*(Mt*s2t*(2880*Dmsqst2*(
        2 + lmMst1 - lmMst2) + 3*(2269 + 664*lmMst1 + 512*xDR2DRMOD + 8*lmMst2*
        (-7 + 64*xDR2DRMOD))*pow2(Mst2)) + 64*(53 + 24*lmMst2)*Mst2*pow2(Mt) +
        48*(5 - 26*lmMst1 + 18*lmMst2)*xDR2DRMOD*pow2(s2t)*pow3(Mst2))) + 36*(-
        (Dmglst2*(Mt*s2t*(10800*Dmsqst2*(1 - 9*lmMst1 + 9*lmMst2) + (364291 -
        88560*lmMst1 + 11520*xDR2DRMOD + 1080*lmMst2*(137 + 32*xDR2DRMOD))*
        pow2(Mst2)) + 30*(3643 - 120*lmMst1 + 3192*lmMst2)*Mst2*pow2(Mt) + 360*
        (7*lmMst1 - 15*lmMst2)*xDR2DRMOD*pow2(s2t)*pow3(Mst2))) + 5*Mst2*(288*
        s2t*(-8*(-2 + 3*lmMst2)*Mt + (-2 - 5*lmMst1 + 6*lmMst2)*Mst2*s2t)*
        xDmglst2*xDR2DRMOD*pow2(Dmglst2) + Mt*s2t*(-2160*Dmsqst2*(5 + 3*lmMst1
        - 3*lmMst2) - 2*(3937 + 2988*lmMst1 + 1152*xDR2DRMOD + 72*lmMst2*(-31 +
        16*xDR2DRMOD))*pow2(Mst2)) - 3*(-353 + 72*lmMst1 + 696*lmMst2)*Mst2*
        pow2(Mt) + 72*(4 + 13*lmMst1 - 9*lmMst2)*xDR2DRMOD*pow2(s2t)*pow3(Mst2)
        ))*pow4(Mst1) + 2*(-(Mt*s2t*(3*Mst2*(266863 + 396360*lmMst1 + 103680*
        xDR2DRMOD + 3240*lmMst2*(-107 + 32*xDR2DRMOD)) + Dmglst2*(15057833 -
        4014360*lmMst1 + 311040*xDR2DRMOD + 3240*lmMst2*(1717 + 288*xDR2DRMOD))
        )) + 45*(38401 + 1080*lmMst1 - 7992*lmMst2)*pow2(Mt) + 6480*Mst2*(
        Dmglst2*(-7*lmMst1 + 15*lmMst2) + (4 + 13*lmMst1 - 9*lmMst2)*Mst2)*
        xDR2DRMOD*pow2(s2t))*pow6(Mst1) + 6480*xDR2DRMOD*(7*Dmglst2*Mst2 + 20*
        xDmglst2*pow2(Dmglst2) - 13*pow2(Mst2))*pow2(s2t)*pow6(Mst2)) + pow2(
        s2t)*(300*Dmsqst2*(pow2(MuSUSY)*(-27*(72*Dmglst2*(2 + lmMst1 - lmMst2)
        + Mst2*(53 + 12*(-1 + 2*lmMst1 - 2*lmMst2)*shiftst1))*pow2(Mst1)*pow3(
        Mst2) - 9*Mst2*(72*Dmglst2*(3 + 8*lmMst1 - 8*lmMst2) + Mst2*(1097 + 72*
        (lmMst1 - lmMst2)*(-1 + shiftst1)))*pow4(Mst1) + 324*(-2*Dmglst2 + Mst2
        + Mst2*shiftst1)*pow5(Mst2)) + (-((20701 + 648*(lmMst1 - lmMst2)*(-1 +
        shiftst1))*pow2(MuSUSY)) + 162*Dmglst2*(-93 + 10*lmMst1 - 10*lmMst2)*
        Mst2*(-1 + pow2(Sbeta))*pow2(Sbeta))*pow6(Mst1)) + Mst2*pow2(MuSUSY)*(
        4*Dmglst2*(162*(6803 - 1810*lmMst1 + 2670*lmMst2)*pow2(Mst2)*pow4(Mst1)
        - 135*(648*lmMst1 - 7*(233 + 240*lmMst2))*pow2(Mst1)*pow4(Mst2) + (
        4256978 - 615600*lmMst1 + 754920*lmMst2)*pow6(Mst1) - 22680*pow6(Mst2))
        + 3*Mst2*(-6*(533629 - 1080*shiftst3 + 180*lmMst1*(-5 + 60*shiftst1 +
        12*shiftst3) - 180*lmMst2*(-1 + 60*shiftst1 + 12*shiftst3))*pow2(Mst2)*
        pow4(Mst1) - 90*(5597 - 360*shiftst1 - 72*shiftst3 - 24*lmMst2*(46 +
        30*shiftst1 + 3*shiftst3) + 12*lmMst1*(47 + 60*shiftst1 + 6*shiftst3))*
        pow2(Mst1)*pow4(Mst2) - (6649153 - 6480*shiftst3 + 540*lmMst1*(-151 +
        120*shiftst1 + 36*shiftst3) - 540*lmMst2*(-143 + 120*shiftst1 + 36*
        shiftst3))*pow6(Mst1) + 1080*(19 + 30*shiftst1 + 3*shiftst3)*pow6(Mst2)
        ))))))/270.) - (Mt*threeLoopFlag*pow2(Al4p)*pow3(Mst2)*(-75*z3*pow4(
        Mst1)*(45*s2t*xDmsqst2*pow2(Dmsqst2)*(Mt*(3456*Mst2*Mt - 30208*s2t*
        pow2(Mst1) - 13209*s2t*pow2(Mst2))*pow2(MuSUSY) - 216*s2t*(-3*Mt + 10*
        Mst2*s2t)*(-1 + pow2(Sbeta))*pow2(Sbeta)*pow4(Mst1)) + Mt*(-12*
        xDmglst2*pow2(Dmglst2)*pow2(MuSUSY)*(8*pow2(Mst1)*(-161918*Mst2*Mt*s2t
        + 77334*pow2(Mt) - 114571*pow2(Mst2)*pow2(s2t)) - 3*pow2(Mst2)*(262184*
        Mst2*Mt*s2t - 87552*pow2(Mt) + 122917*pow2(Mst2)*pow2(s2t))) + 8*pow2(
        s2t)*(-24300*Dmglst2*Dmsqst2*Mst2*(-1 + pow2(Sbeta))*pow2(Sbeta)*pow4(
        Mst1) + pow2(MuSUSY)*(30*Dmsqst2*(9*Mst2*(288*Dmglst2 + 35*Mst2)*pow2(
        Mst1) + 27*(36*Dmglst2 - 23*Mst2)*pow3(Mst2) + 1771*pow4(Mst1)) + Mst2*
        (Dmglst2*(1052676*pow2(Mst1)*pow2(Mst2) + 1412464*pow4(Mst1) + 475605*
        pow4(Mst2)) + 6*Mst2*(3*(35719 + 108*lmMst1 - 108*lmMst2)*pow2(Mst1)*
        pow2(Mst2) + 2*(91963 + 162*lmMst1 - 162*lmMst2)*pow4(Mst1) + 27*(1319
        + 6*lmMst1 - 6*lmMst2)*pow4(Mst2))))))) + 18*s2t*xDmsqst2*pow2(Dmsqst2)
        *pow2(Mst1)*(-300*Mt*(2160*(3 + 2*lmMst1 - 2*lmMst2)*Mst2*Mt*pow2(Mst1)
        + 270*s2t*shiftst1*(2*(lmMst1 - lmMst2)*pow2(Mst1) - pow2(Mst2))*(pow2(
        Mst1) + pow2(Mst2)) + 40*s2t*T1ep*pow2(Mst1)*(14*pow2(Mst1) + 3*pow2(
        Mst2)))*pow2(MuSUSY) + s2t*(10125*pow2(Sbeta)*(2*(11 + 39*lmMst1 - 39*
        lmMst2)*Mst2*s2t*(-1 + pow2(Sbeta)) + (167 - 5*lmMst1 + 5*lmMst2)*Mt*
        pow2(Sbeta))*pow6(Mst1) + Mt*(pow2(MuSUSY)*(375*(1925 + 64*OepS2 -
        27864*S2 - 72*(lmMst1 - lmMst2)*(-1 + 18*S2))*pow2(Mst1)*pow2(Mst2) +
        2*(7*(-46499 + 8000*OepS2 - 1728000*S2) + 90*lmMst2*(-2137 + 12600*S2)
        - 90*lmMst1*(-2137 + 1080*lmMst2 + 12600*S2) + 48600*(pow2(lmMst1) +
        pow2(lmMst2)))*pow4(Mst1) + 162000*pow4(Mst2)) + 10125*(-167 + 5*lmMst1
        - 5*lmMst2)*pow2(Sbeta)*pow6(Mst1)))) + Mt*(50*pow2(s2t)*(4860*Dmglst2*
        Dmsqst2*(-1081 + 165*lmMst1 - 165*lmMst2)*Mst2*(-1 + pow2(Sbeta))*pow2(
        Sbeta)*pow8(Mst1) - (pow2(MuSUSY)*(120*Dmsqst2*(1157625*(288*Dmglst2*(4
        + lmMst1 - lmMst2) - Mst2*(1141 - 32*OepS2 - 4860*S2 + 24*lmMst1*(14 +
        27*S2 + 6*lmMst2*(-1 + shiftst1) - 9*shiftst1) - 24*lmMst2*(11 + 27*S2
        - 6*shiftst1) + 108*shiftst1 - 144*(-1 + shiftst1)*pow2(lmMst2)))*pow3(
        Mst2)*pow4(Mst1) + 41674500*(12*Dmglst2 + Mst2 + 2*lmMst2*Mst2*(-1 +
        shiftst1) - 3*Mst2*shiftst1)*pow2(Mst1)*pow5(Mst2) - 3087*Mst2*(13500*
        Dmglst2*(23 - 42*lmMst1 + 42*lmMst2) + Mst2*(817051 - 52000*OepS2 -
        229500*S2 - 540*lmMst2*(311 + 1950*S2 - 150*shiftst1) + 540*lmMst1*(311
        + 1950*S2 - 150*shiftst1 + 10*lmMst2*(1 + 10*shiftst1)) - 29700*pow2(
        lmMst1) + 2700*(9 - 20*shiftst1)*pow2(lmMst2)))*pow6(Mst1) - 5*(
        593331163 - 60642400*OepS2 + 1143733500*S2 + 1260*lmMst2*(8263 -
        974610*S2 + 39690*shiftst1) - 1260*lmMst1*(8263 - 974610*S2 + 39690*
        shiftst1 - 3780*lmMst2*(19 + 7*shiftst1)) - 61916400*pow2(lmMst1) -
        4762800*(6 + 7*shiftst1)*pow2(lmMst2))*pow8(Mst1)) - Mst2*(4*Dmglst2*(-
        2315250*(14267 - 432*B4 + 576*D3 - 360*DN + 4752*lmMst1 - 224*OepS2 +
        324*(-1677 + 14*lmMst1 - 14*lmMst2)*S2 - 1404*pow2(lmMst1) + 72*lmMst2*
        (-63 - 232*lmMst1 + 16*pow2(lmMst1)) - 72*(-281 + 123*lmMst1)*pow2(
        lmMst2) + 7704*pow3(lmMst2))*pow4(Mst1)*pow4(Mst2) + 24696*pow2(Mst2)*(
        -2695042 + 40500*B4 - 54000*D3 + 33750*DN + 326895*lmMst1 + 294105*
        lmMst2 + 1935450*lmMst1*lmMst2 + 258000*OepS2 + 105590250*S2 - 5224500*
        lmMst1*S2 + 5224500*lmMst2*S2 + 324900*pow2(lmMst1) - 938250*lmMst2*
        pow2(lmMst1) - 2260350*pow2(lmMst2) + 2524500*lmMst1*pow2(lmMst2) -
        11250*pow3(lmMst1) - 1575000*pow3(lmMst2))*pow6(Mst1) - 83349000*(-20 +
        2*(99 + 16*lmMst1)*lmMst2 + 123*pow2(lmMst2))*pow2(Mst1)*pow6(Mst2) + (
        -219487237477 + 1000188000*B4 - 1333584000*D3 + 833490000*DN +
        25338520620*lmMst1 - 6418297620*lmMst2 + 57402853200*lmMst1*lmMst2 +
        23746576000*OepS2 + 3479135436000*S2 - 480868164000*lmMst1*S2 +
        480868164000*lmMst2*S2 + 6305153400*pow2(lmMst1) - 49175910000*lmMst2*
        pow2(lmMst1) - 63708006600*pow2(lmMst2) + 109687284000*lmMst1*pow2(
        lmMst2) + 1278018000*pow3(lmMst1) - 61789392000*pow3(lmMst2))*pow8(
        Mst1) + 2667168000*lmMst2*(1 + lmMst2)*pow8(Mst2)) + 3*Mst2*(1543500*(
        17297 + 1440*B4 - 144*D3 + 72*DN + 20064*lmMst1 - 2208*OepS2 + 972*(307
        + 46*lmMst1 - 46*lmMst2)*S2 + 216*(-1 + lmMst1 - lmMst2)*(-1 + 2*
        lmMst2)*shiftst3 + 1908*pow2(lmMst1) - 12*lmMst2*(3508 - 1944*lmMst1 +
        117*pow2(lmMst1)) + 1080*shiftst1*(1 + lmMst1*(-2 + 4*lmMst2) - 4*pow2(
        lmMst2)) + 72*(-470 + 147*lmMst1)*pow2(lmMst2) + 576*pow3(lmMst1) -
        9756*pow3(lmMst2))*pow4(Mst1)*pow4(Mst2) + 20580*pow2(Mst2)*(10552777 +
        205200*B4 - 5400*DN + 4145280*lmMst1 - 845600*OepS2 + 8100*(11243 +
        2114*lmMst1 - 2114*lmMst2)*S2 + 162000*(lmMst1 - lmMst2)*(-1 + 2*
        lmMst2)*shiftst1 - 1181700*pow2(lmMst1) - 240*lmMst2*(16192 - 26430*
        lmMst1 + 3465*pow2(lmMst1)) + 16200*shiftst3*(1 + lmMst1*(-2 + 4*
        lmMst2) - 4*pow2(lmMst2)) + 900*(-5735 + 3072*lmMst1)*pow2(lmMst2) -
        55800*pow3(lmMst1) - 1877400*pow3(lmMst2))*pow6(Mst1) - 55566000*(103 +
        32*lmMst1*(1 + lmMst2) - 30*shiftst1 - 3*shiftst3 + 6*lmMst2*(31 + 10*
        shiftst1 + shiftst3) + 91*pow2(lmMst2))*pow2(Mst1)*pow6(Mst2) + (
        314145861731 + 4223016000*B4 - 111132000*DN + 7960477140*lmMst1 -
        36322328000*OepS2 + 18522000*(194357 + 39711*lmMst1 - 39711*lmMst2)*S2
        + 3333960000*(lmMst1 - lmMst2)*(-1 + 2*lmMst2)*shiftst1 + 333396000*(-1
        + 3*lmMst1 - 3*lmMst2)*(-1 + 2*lmMst2)*shiftst3 - 105268993200*pow2(
        lmMst1) - 1260*lmMst2*(4972789 - 230465340*lmMst1 + 36250200*pow2(
        lmMst1)) + 3175200*(-58301 + 37135*lmMst1)*pow2(lmMst2) - 2444904000*
        pow3(lmMst1) - 69790896000*pow3(lmMst2))*pow8(Mst1) + 889056000*pow2(1
        + lmMst2)*pow8(Mst2)))))/171500.) - 10*pow2(MuSUSY)*(30*pow4(Mst1)*((-
        4*Mt*s2t*(-3*Mst2*(3764*T1ep + 99632*z3 + 11764*z4 - 20505*pow2(z2)) +
        Dmglst2*(65684*T1ep - 1582324*z3 - 2636*z4 + 259215*pow2(z2))) - pow2(
        s2t)*(13260*Dmsqst2*(4*T1ep + 2*z4 + 3*pow2(z2)) + Mst2*(-8*Dmglst2*(
        17308*T1ep + 11327*z4 + 15897*pow2(z2)) + 3*Mst2*(52948*T1ep + 27446*z4
        + 35823*pow2(z2)))) - 3*pow2(Mt)*(-383185 + 2592*B4 - 2592*D3 + 1296*DN
        + 187704*lmMst1 - 374328*lmMst2 + 136080*lmMst1*lmMst2 - 10368*lmMt +
        20736*lmMst1*lmMt - 20736*lmMst2*lmMt + 6912*lmMst1*lmMst2*lmMt +
        17440*OepS2 + 36936*S2 - 353160*lmMst1*S2 + 353160*lmMst2*S2 - 26160*
        T1ep - 18024*z3 - 24744*z4 + 7992*pow2(lmMst1) + 5616*lmMst2*pow2(
        lmMst1) - 3456*lmMt*pow2(lmMst1) - 185544*pow2(lmMst2) + 53136*lmMst1*
        pow2(lmMst2) - 3456*lmMt*pow2(lmMst2) - 19620*pow2(z2) - 720*pow3(
        lmMst1) - 58032*pow3(lmMst2)))*pow4(Mst1) - 18*pow2(Mst2)*(-12960*
        Dmsqst2*Mst2*Mt*s2t*z3 - 30*xDmsqst2*pow2(Dmsqst2)*pow2(s2t)*(2*z4 + 3*
        pow2(z2)) - xDmglst2*pow2(Dmglst2)*(648*z4*pow2(Mt) + 8*Mst2*Mt*s2t*(
        109*z4 - 1902*pow2(z2)) + pow2(Mst2)*pow2(s2t)*(3172*z4 + 627*pow2(z2))
        ) + 6*pow2(Mst2)*(15*Dmsqst2*pow2(s2t)*(4*T1ep + 2*z4 + 3*pow2(z2)) +
        4*pow2(Mt)*(-436 + 6*B4 - 6*D3 + 3*DN + 48*lmMst1 - 408*lmMst2 + 96*
        lmMst1*lmMst2 - 24*lmMt + 972*S2 - 230*z3 - 27*z4 - 192*pow2(lmMst2) +
        48*lmMst1*pow2(lmMst2) - 48*pow3(lmMst2))) - 3*Mt*s2t*(84*T1ep + 11930*
        z3 + 510*z4 - 1665*pow2(z2))*pow3(Mst2) - 3*Dmglst2*(-1440*Dmsqst2*Mt*
        s2t*z3 + Mt*s2t*pow2(Mst2)*(-140*T1ep + 30274*z3 + 542*z4 - 5289*pow2(
        z2)) - 16*Mst2*pow2(Mt)*(-540 + 6*B4 - 6*D3 + 3*DN - 48*lmMst1 - 432*
        lmMst2 + 648*S2 + 83*z3 - 27*z4 - 96*pow2(lmMst2) + 48*lmMst1*pow2(
        lmMst2) - 48*pow3(lmMst2)) + 2*pow2(s2t)*(28*T1ep + 212*z4 + 237*pow2(
        z2))*pow3(Mst2)) + 9*pow2(s2t)*(92*T1ep + 40*z4 - 3*pow2(z2))*pow4(
        Mst2)) - 6*pow2(Mst1)*(-2*xDmglst2*pow2(Dmglst2)*(Mst2*Mt*s2t*(8090*z4
        - 37437*pow2(z2)) + pow2(Mst2)*pow2(s2t)*(7489*z4 - 795*pow2(z2)) + 6*
        pow2(Mt)*(338*z4 + 21*pow2(z2))) - 6*Dmglst2*(-10800*Dmsqst2*Mt*s2t*z3
        - 2*Mt*s2t*pow2(Mst2)*(1012*T1ep - 47456*z3 - 709*z4 + 8535*pow2(z2)) +
        Mst2*pow2(Mt)*(28405 - 288*B4 + 288*D3 - 144*DN + 5472*lmMst1 + 25632*
        lmMst2 - 576*lmMst1*lmMst2 - 1152*lmMt + 576*lmMst1*lmMt - 576*lmMst2*
        lmMt + 224*OepS2 - 21060*S2 - 4536*lmMst1*S2 + 4536*lmMst2*S2 - 336*
        T1ep - 3000*z3 + 1128*z4 + 5184*pow2(lmMst2) - 4608*lmMst1*pow2(lmMst2)
        - 252*pow2(z2) + 4608*pow3(lmMst2)) + 6*pow2(s2t)*(172*T1ep + 185*z4 +
        237*pow2(z2))*pow3(Mst2)) + 3*(-30240*Dmsqst2*Mst2*Mt*s2t*z3 - 140*
        xDmsqst2*pow2(Dmsqst2)*pow2(s2t)*(2*z4 + 3*pow2(z2)) + 3*pow2(Mst2)*(
        130*Dmsqst2*pow2(s2t)*(4*T1ep + 2*z4 + 3*pow2(z2)) + pow2(Mt)*(-10667 +
        96*B4 - 96*D3 + 48*DN + 3072*lmMst1 - 9600*lmMst2 + 3456*lmMst1*lmMst2
        - 384*lmMt + 384*lmMst1*lmMt - 384*lmMst2*lmMt + 224*OepS2 + 13932*S2 -
        4536*lmMst1*S2 + 4536*lmMst2*S2 - 336*T1ep - 2616*z3 - 600*z4 - 4992*
        pow2(lmMst2) + 1536*lmMst1*pow2(lmMst2) - 252*pow2(z2) - 1536*pow3(
        lmMst2))) - 4*Mt*s2t*(272*T1ep + 14636*z3 + 1135*z4 - 2388*pow2(z2))*
        pow3(Mst2) + pow2(s2t)*(4228*T1ep + 2276*z4 + 2523*pow2(z2))*pow4(Mst2)
        ))) + s2t*(Mt*pow2(Mst1)*(Dmglst2*(36*pow2(Mst2)*(66761 + 301320*B4 -
        4860*DN - 205380*lmMst1 + 40480*OepS2 - 216*(2489 + 3795*lmMst1 - 3795*
        lmMst2)*S2 + 23760*pow2(lmMst1) + 180*lmMst2*(4993 - 1956*lmMst1 + 48*
        pow2(lmMst1)) - 1080*(-482 + 331*lmMst1)*pow2(lmMst2) + 348840*pow3(
        lmMst2))*pow4(Mst1) + 2332800*Dmsqst2*((-2 - 5*lmMst1 + 5*lmMst2)*pow2(
        Mst1)*pow2(Mst2) + (8 - 15*lmMst1 + 15*lmMst2)*pow4(Mst1)) + 27*pow2(
        Mst1)*(23917 + 188640*B4 - 3600*DN - 37440*lmMst1 + 5600*OepS2 - 324*(-
        453 + 350*lmMst1 - 350*lmMst2)*S2 + 11520*pow2(lmMst1) - 2880*lmMst2*(-
        237 + 55*lmMst1 + 4*pow2(lmMst1)) - 1440*(-280 + 121*lmMst1)*pow2(
        lmMst2) + 185760*pow3(lmMst2))*pow4(Mst2) - 10*(2773621 - 1660176*B4 +
        25272*DN + 2004408*lmMst1 - 525472*OepS2 + 108*(123113 + 98526*lmMst1 -
        98526*lmMst2)*S2 + 3888*pow2(lmMst1) - 144*lmMst2*(36802 - 11421*lmMst1
        + 1728*pow2(lmMst1)) + 167184*(-14 + 15*lmMst1)*pow2(lmMst2) - 31104*
        pow3(lmMst1) - 2227824*pow3(lmMst2))*pow6(Mst1) + 622080*(-1 + 2*lmMst2
        + 3*pow2(lmMst2))*pow6(Mst2)) + 15*Mst2*(24*pow2(Mst2)*(75569 + 13716*
        B4 - 54*DN - 33426*lmMst1 - 1088*OepS2 + 162*(169 + 136*lmMst1 - 136*
        lmMst2)*S2 - 2376*pow2(lmMst1) + 54*lmMst2*(1427 - 1012*lmMst1 + 16*
        pow2(lmMst1)) - 108*(-642 + 203*lmMst1)*pow2(lmMst2) + 21060*pow3(
        lmMst2))*pow4(Mst1) + 155520*Dmsqst2*((1 + lmMst1 - lmMst2)*pow2(Mst1)*
        pow2(Mst2) + (1 + 3*lmMst1 - 3*lmMst2)*pow4(Mst1)) + 27*pow2(Mst1)*(
        28683 + 5280*B4 - 48*DN - 5952*lmMst1 - 224*OepS2 + 324*(-1 + 14*lmMst1
        - 14*lmMst2)*S2 - 768*pow2(lmMst1) - 192*lmMst2*(-214 + 73*lmMst1 + 4*
        pow2(lmMst1)) - 96*(-268 + 57*lmMst1)*pow2(lmMst2) + 6240*pow3(lmMst2))
        *pow4(Mst2) + 2*(1702429 + 257904*B4 - 648*DN - 748656*lmMst1 - 30112*
        OepS2 + 108*(9185 + 5646*lmMst1 - 5646*lmMst2)*S2 + 41904*pow2(lmMst1)
        + 216*lmMst2*(5971 - 6106*lmMst1 + 576*pow2(lmMst1)) - 41904*(-34 + 15*
        lmMst1)*pow2(lmMst2) - 3456*pow3(lmMst1) + 507600*pow3(lmMst2))*pow6(
        Mst1) + 41472*pow2(1 + lmMst2)*pow6(Mst2))) + 4860*xDR2DRMOD*(64*Mt*
        pow2(Mst1)*(2*pow2(Mst2)*((1 + lmMst2)*Mst2*(1 - 10*lmMst2 + 4*lmMst1*(
        2 + lmMst2) - 4*pow2(lmMst2)) + 2*Dmglst2*(8 + 13*lmMst2 - 8*pow2(
        lmMst2) + lmMst1*(-5 + 5*lmMst2 + 6*pow2(lmMst2)) - 6*pow3(lmMst2)))*
        pow4(Mst1) + (1 + lmMst2)*(-2*(1 + 5*lmMst2 - lmMst1*(3 + 2*lmMst2) +
        2*pow2(lmMst2))*pow2(Mst1)*pow5(Mst2) + Mst2*(7 - 32*lmMst2 + 4*lmMst1*
        (7 + 3*lmMst2) - 12*pow2(lmMst2))*pow6(Mst1)) + Dmglst2*((49 + 103*
        lmMst2 - 36*(1 + lmMst2)*pow2(lmMst2) + 4*lmMst1*(-11 + 6*lmMst2 + 9*
        pow2(lmMst2)))*pow6(Mst1) + 2*(pow2(Mst1)*(8 + 7*lmMst2 - 11*pow2(
        lmMst2) + lmMst1*(-3 + 5*lmMst2 + 6*pow2(lmMst2)) - 6*pow3(lmMst2))*
        pow4(Mst2) + (1 - 2*lmMst2 - 3*pow2(lmMst2))*pow6(Mst2))) - 2*pow2(1 +
        lmMst2)*pow7(Mst2)) + s2t*((75*Dmsqst2*(-1 + 2*lmMst1 - 2*lmMst2) -
        pow2(Mst2)*(189 + 726*lmMst2 + 32*(1 + lmMst2)*pow2(lmMst1) + 707*pow2(
        lmMst2) - 2*lmMst1*(253 + 332*lmMst2 + 123*pow2(lmMst2)) + 214*pow3(
        lmMst2)) - 2*pow2(Mst1)*(32 + 285*lmMst2 + 144*(1 + lmMst2)*pow2(
        lmMst1) + 444*pow2(lmMst2) - lmMst1*(253 + 588*lmMst2 + 379*pow2(
        lmMst2)) + 235*pow3(lmMst2)))*pow4(Mst1)*pow4(Mst2) + 2*pow2(Mst2)*(75*
        Dmsqst2*(lmMst1 - lmMst2) - pow2(Mst1)*(40 + 277*lmMst2 + 272*(1 +
        lmMst2)*pow2(lmMst1) + 556*pow2(lmMst2) - lmMst1*(237 + 828*lmMst2 +
        635*pow2(lmMst2)) + 363*pow3(lmMst2)))*pow6(Mst1) + pow2(Mst1)*(135*
        shiftst1*(Dmsqst2 + pow2(Mst2))*((pow2(Mst1) + pow2(Mst2))*pow4(Mst2) -
        2*(lmMst1 - lmMst2)*pow2(Mst1)*(pow2(Mst1)*pow2(Mst2) + pow4(Mst1) +
        pow4(Mst2))) + 15*Dmsqst2*(-4*Dmglst2*(lmMst1 - lmMst2)*Mst2*pow4(Mst1)
        + 10*(lmMst1 - lmMst2)*pow6(Mst1) - 5*pow6(Mst2)) + 2*Dmglst2*Mst2*(15*
        Dmsqst2*pow4(Mst2) + pow2(Mst1)*(15*Dmsqst2*(1 - 2*lmMst1 + 2*lmMst2)*
        pow2(Mst2) - (60 + 206*lmMst2 + 32*lmMst2*pow2(lmMst1) + lmMst1*(8 -
        460*lmMst2 - 246*pow2(lmMst2)) + 519*pow2(lmMst2) + 214*pow3(lmMst2))*
        pow4(Mst2)) - 2*(pow2(Mst2)*(48 + 4*lmMst2*(31 + 36*pow2(lmMst1)) +
        278*pow2(lmMst2) - lmMst1*(44 + 278*lmMst2 + 379*pow2(lmMst2)) + 235*
        pow3(lmMst2))*pow4(Mst1) + (48 + 4*lmMst2*(45 + 68*pow2(lmMst1)) + 310*
        pow2(lmMst2) - lmMst1*(92 + 310*lmMst2 + 635*pow2(lmMst2)) + 363*pow3(
        lmMst2))*pow6(Mst1)) - (-20 + 2*(99 + 16*lmMst1)*lmMst2 + 123*pow2(
        lmMst2))*pow6(Mst2)) - (205 + 252*lmMst2 + 32*lmMst1*(1 + lmMst2) + 91*
        pow2(lmMst2))*pow8(Mst2)) + 16*(1 + lmMst2)*(4*Dmglst2*lmMst2 + Mst2 +
        lmMst2*Mst2)*pow9(Mst2))))))))/25.) + pow2(Mt)*(3*Tbeta*pow2(Mst2)*
        pow2(MuSUSY)*(729*oneLoopFlag*pow2(s2t)*pow3(Mst2)*pow4(Mst1)*(-2*(
        lmMst1 - lmMst2)*pow2(Mst1)*(pow2(Mst1) + pow2(Mst2)) + (2 - lmMst1 +
        lmMst2)*pow4(Mst2)) + Al4p*(8*xDmglst2*pow2(Dmglst2)*(-243*s2t*
        twoLoopFlag*pow2(Mst1)*(2*(-8*(10 + lmMst1 - lmMst2)*Mt + Mst2*s2t*(8 -
        9*lmMst1 + 9*lmMst2 + 2*lmMst1*lmMst2 - 2*pow2(lmMst2)))*pow2(Mst2)*
        pow4(Mst1) + 2*(-4*(8 + lmMst1 - lmMst2)*Mt + Mst2*s2t*(8 - 5*lmMst1 +
        4*lmMst2 + 2*lmMst1*lmMst2 - 2*pow2(lmMst2)))*pow2(Mst1)*pow4(Mst2) + (
        -268*Mt + Mst2*s2t*(17 + 4*lmMst1*(-7 + lmMst2) + 28*lmMst2 - 4*pow2(
        lmMst2)))*pow6(Mst1) - 2*(-2 + lmMst2)*s2t*pow7(Mst2)) + (Al4p*Mst2*
        threeLoopFlag*(75*pow2(Mst2)*(-8*Mst2*Mt*s2t*(334138 - 61560*B4 + 1620*
        DN - 236520*lmMst1 - 180*(-823 + 438*lmMst1)*lmMst2 + 2240*OepS2 - 81*(
        -1373 + 560*lmMst1 - 560*lmMst2)*S2 - 3360*T1ep - 4320*pow2(lmMst1) +
        1080*(29 + 48*lmMst1)*pow2(lmMst2) - 51840*pow3(lmMst2)) + pow2(Mst2)*
        pow2(s2t)*(13169 - 41040*B4 + 43200*D3 - 22680*DN + 282960*lmMst1 +
        1120*OepS2 - 324*(65819 + 70*lmMst1 - 70*lmMst2)*S2 - 1680*T1ep -
        13500*pow2(lmMst1) + 720*lmMst2*(-775 + 12*lmMst1 + 24*pow2(lmMst1)) -
        1080*(-2 + 123*lmMst1)*pow2(lmMst2) + 115560*pow3(lmMst2)) + 96*pow2(
        Mt)*(33934 - 90*B4 + 90*D3 - 45*DN + 120*(163 + 24*lmMst1)*lmMst2 -
        120*lmMt - 10206*S2 - 720*(2 + lmMst1)*pow2(lmMst2) + 720*(lmMst1 +
        pow3(lmMst2))))*pow4(Mst1) - 40500*s2t*(Mst2*s2t*(812 - 32*lmMst1*(-2 +
        lmMst2) + 38*lmMst2 - 251*pow2(lmMst2)) + 128*Mt*(3 + 4*lmMst2 - 3*
        pow2(lmMst2)))*pow2(Mst1)*pow5(Mst2) + 2*(-2*pow2(Mst2)*pow2(s2t)*(
        2011073 + 1417500*B4 - 1458000*D3 + 749250*DN + 934245*lmMst1 -
        1178000*OepS2 + 1350*(620417 + 17670*lmMst1 - 17670*lmMst2)*S2 +
        1767000*T1ep + 3150900*pow2(lmMst1) - 45*lmMst2*(-124139 + 189090*
        lmMst1 + 71550*pow2(lmMst1)) + 4050*(1323 + 1970*lmMst1)*pow2(lmMst2) +
        101250*pow3(lmMst1) - 4860000*pow3(lmMst2)) - 125*Mst2*Mt*s2t*(996211 -
        295488*B4 + 7776*DN - 1030176*lmMst1 - 144*(-1883 + 3618*lmMst1)*lmMst2
        + 98336*OepS2 - 756*(259 + 2634*lmMst1 - 2634*lmMst2)*S2 - 147504*T1ep
        + 49248*pow2(lmMst1) + 5184*(67 + 48*lmMst1)*pow2(lmMst2) - 248832*
        pow3(lmMst2)) + 150*pow2(Mt)*(2199511 - 4320*B4 + 4320*D3 - 2160*DN +
        140160*lmMst1 + 960*(1426 + 303*lmMst1)*lmMst2 - 2880*(-16 + 5*lmMst1 -
        5*lmMst2)*lmMt - 1120*OepS2 + 324*(-3411 + 70*lmMst1 - 70*lmMst2)*S2 +
        1680*T1ep - 2880*(77 + 24*lmMst1)*pow2(lmMst2) + 69120*pow3(lmMst2)))*
        pow6(Mst1) + 1296000*(2 + lmMst2 - 3*pow2(lmMst2))*pow2(s2t)*pow8(Mst2)
        ))/500.) + Mst2*s2t*(972*twoLoopFlag*pow2(Mst1)*(4*(4*Mt*(5 + 6*lmMst2
        - lmMst1*(4 + 3*lmMst2) + 3*pow2(lmMst2)) + Mst2*s2t*(-1 + 13*lmMst2 -
        lmMst1*(13 + 8*lmMst2) + pow2(lmMst1) + 7*pow2(lmMst2)))*pow3(Mst2)*
        pow4(Mst1) + 2*(-(Mst2*s2t*(-14 - 20*lmMst2 + 2*lmMst1*(5 + 3*lmMst2) +
        pow2(lmMst1) - 7*pow2(lmMst2))) + 8*Mt*(4 + 3*lmMst2 - lmMst1*(1 +
        lmMst2) + pow2(lmMst2)))*pow2(Mst1)*pow5(Mst2) + 4*Dmglst2*(Mst2*s2t*(-
        5 + 8*lmMst2 - 4*lmMst1*(2 + lmMst2) + 4*pow2(lmMst2)) + Mt*(65 +
        lmMst1*(34 - 20*lmMst2) - 26*lmMst2 + 20*pow2(lmMst2)))*pow6(Mst1) +
        Mst2*(Mst2*s2t*(-1 + 50*lmMst2 - 2*lmMst1*(25 + 32*lmMst2) + 20*pow2(
        lmMst1) + 44*pow2(lmMst2)) + Mt*(84 + 152*lmMst2 - 40*lmMst1*(3 + 2*
        lmMst2) + 80*pow2(lmMst2)))*pow6(Mst1) + 8*Dmglst2*((Mst2*s2t*(-2 + 3*
        lmMst2 - lmMst1*(3 + 2*lmMst2) + 2*pow2(lmMst2)) + Mt*(22 + lmMst1*(8 -
        6*lmMst2) - 4*lmMst2 + 6*pow2(lmMst2)))*pow2(Mst2)*pow4(Mst1) + (-(
        Mst2*s2t*(1 + lmMst1 - 2*lmMst2 + 2*lmMst1*lmMst2 - 2*pow2(lmMst2))) +
        2*Mt*(6 + lmMst1 + lmMst2 - lmMst1*lmMst2 + pow2(lmMst2)))*pow2(Mst1)*
        pow4(Mst2) + lmMst2*s2t*pow7(Mst2)) + 4*(1 + lmMst2)*s2t*pow8(Mst2)) -
        108*xDR2DRMOD*(9*s2t*pow2(Mst1)*(-(Al4p*threeLoopFlag*xDmsqst2*pow2(
        Dmsqst2)*(10*lmMst1*pow2(Mst1)*(14*Dmglst2*Mst2 + 19*xDmglst2*pow2(
        Dmglst2) + pow2(Mst1) + pow2(Mst2)) - 5*((14*Dmglst2*Mst2 + 19*
        xDmglst2*pow2(Dmglst2) + pow2(Mst2))*(pow2(Mst1) + pow2(Mst2)) + 2*
        lmMst2*pow2(Mst1)*(14*Dmglst2*Mst2 + 19*xDmglst2*pow2(Dmglst2) + pow2(
        Mst1) + pow2(Mst2))))) + 4*twoLoopFlag*(2*Dmglst2*lmMst2*Mst2 + (-2 +
        lmMst2)*xDmglst2*pow2(Dmglst2) + (1 + lmMst2)*pow2(Mst2))*((pow2(Mst1)
        + pow2(Mst2))*pow4(Mst2) - 2*(lmMst1 - lmMst2)*pow2(Mst1)*(pow2(Mst1)*
        pow2(Mst2) + pow4(Mst1) + pow4(Mst2)))) + 2*Al4p*threeLoopFlag*
        xDmglst2*pow2(Dmglst2)*(pow2(Mst2)*(585*Dmsqst2*(-1 + 2*lmMst1 - 2*
        lmMst2)*s2t - 64*Mst2*Mt*(5 + 149*lmMst2 + 12*pow2(lmMst2) + 6*lmMst1*(
        -7 - 8*lmMst2 + 6*pow2(lmMst2)) - 36*pow3(lmMst2)) + s2t*pow2(Mst2)*(-
        1540 - 2506*lmMst2 + 96*(-2 + lmMst2)*pow2(lmMst1) + lmMst1*(2760 + 36*
        lmMst2 - 738*pow2(lmMst2)) + 141*pow2(lmMst2) + 642*pow3(lmMst2)))*
        pow4(Mst1) + 3*pow2(Mst1)*(-195*Dmsqst2*s2t + 128*Mst2*Mt*(-3 - 4*
        lmMst2 + 3*pow2(lmMst2)) + s2t*(-380 + 32*lmMst1*(-2 + lmMst2) - 38*
        lmMst2 + 251*pow2(lmMst2))*pow2(Mst2))*pow4(Mst2) + 2*(585*Dmsqst2*(
        lmMst1 - lmMst2)*s2t + 32*Mst2*Mt*(85 - 215*lmMst2 + lmMst1*(66 + 90*
        lmMst2 - 72*pow2(lmMst2)) - 54*pow2(lmMst2) + 72*pow3(lmMst2)) + 3*s2t*
        pow2(Mst2)*(-336 - 716*lmMst2 + 144*(-2 + lmMst2)*pow2(lmMst1) +
        lmMst1*(668 + 534*lmMst2 - 379*pow2(lmMst2)) - 246*pow2(lmMst2) + 235*
        pow3(lmMst2)))*pow6(Mst1) + 96*s2t*(2 + lmMst2 - 3*pow2(lmMst2))*pow8(
        Mst2)))))) - 162*s2t*xMst*(27*(lmMst1 - lmMst2)*oneLoopFlag*s2t*Tbeta*
        pow2(MuSUSY)*pow3(Mst2) - 2*Al4p*twoLoopFlag*(4*Mst2*Tbeta*(-18*(-2 +
        lmMst2)*(-lmMst1 + lmMst2)*s2t*xDmglst2*xDR2DRMOD*pow2(Dmglst2) +
        Dmglst2*Mt*(785 + 6*lmMst1*(85 - 42*lmMst2) - 438*lmMst2 + 252*pow2(
        lmMst2)) + Mst2*Mt*(193 + 474*lmMst2 - 6*lmMst1*(67 + 42*lmMst2) + 252*
        pow2(lmMst2)) - Dmglst2*Mst2*s2t*(49 - 84*lmMst2 + lmMst1*(84 - 36*
        lmMst2*(-1 + xDR2DRMOD)) + 36*(-1 + xDR2DRMOD)*pow2(lmMst2)) - s2t*(1 +
        3*lmMst2*(-37 + 6*xDR2DRMOD) - 3*lmMst1*(-37 + 6*lmMst2*(-12 +
        xDR2DRMOD) + 6*xDR2DRMOD) - 81*pow2(lmMst1) + 9*(-15 + 2*xDR2DRMOD)*
        pow2(lmMst2))*pow2(Mst2))*pow2(MuSUSY) + xDmglst2*pow2(Dmglst2)*(48*(
        143 - 18*lmMst1 + 18*lmMst2)*Mt*Tbeta*pow2(MuSUSY) - 2*Mst2*s2t*Tbeta*(
        157 + 348*lmMst2 + 12*lmMst1*(-29 + 3*lmMst2) - 36*pow2(lmMst2))*pow2(
        MuSUSY) + 60*(43 - 60*lmMst1 + 60*lmMst2)*Mst2*Mt*MuSUSY*pow2(Sbeta) -
        15*(-43 + 60*lmMst1 - 60*lmMst2)*s2t*Tbeta*(-1 + pow2(Sbeta))*pow2(
        Sbeta)*pow3(Mst2))))*power10(Mst1)))/(17496.*Tbeta*pow4(Mst1)*pow9(
        Mst2));
}

/**
 * @return The diagonal (2, 2) matrix element of the Higgs mass matrix as a double for the hierarchy 'H6b'
 */
double H6b::getS2() const {
   return -(oneLoopFlag*((4*Mt*MuSUSY*s2t*(4*(lmMst1 - lmMst2)*pow2(Mt) + ((-2 -
        lmMst1 + lmMst2)*pow2(Mst1) + (2 - lmMst1 + lmMst2)*pow2(Mst2))*pow2(
        s2t)))/Tbeta + 4*pow2(Mt)*pow2(s2t)*(2*(lmMst1 - lmMst2)*(pow2(Mst1) -
        pow2(Mst2)) + pow2(MuSUSY)*(2 - lmMst1 + lmMst2 - (2*(lmMst1 - lmMst2)*
        pow2(Mst1)*(pow2(Mst1) + pow2(Mst2)))/pow4(Mst2))) - (4*pow2(Mt)*pow2(
        MuSUSY)*pow2(s2t)*(2 - lmMst1 + lmMst2 - (2*(lmMst1 - lmMst2)*pow2(
        Mst1)*(pow2(Mst1) + pow2(Mst2)))/pow4(Mst2)))/pow2(Sbeta) + 16*(lmMst1
        + lmMst2 - 2*lmMt)*pow4(Mt) - (-4*pow2(Mst1)*pow2(Mst2) + (2 + lmMst1 -
        lmMst2)*pow4(Mst1) + (2 - lmMst1 + lmMst2)*pow4(Mst2))*pow4(s2t)))/32.
         - (threeLoopFlag*pow2(Al4p)*(12*xDmglst2*pow2(Dmglst2)*(-2*pow2(Mst1)*
        (-8*pow2(Mst2)*pow2(Mt)*(-(Tbeta*pow2(MuSUSY)*pow2(s2t)*(-1 + pow2(
        Sbeta))*(-229142*z3 - 7489*z4 + 795*pow2(z2))) + pow2(Sbeta)*(-72*Mt*
        MuSUSY*s2t*(2534*z3 + 33*z4 + 90*pow2(z2)) - 2*Tbeta*pow2(Mt)*(-619510*
        z3 + 11338*z4 + 17007*pow2(z2)))) - 4*s2t*pow2(Mt)*pow2(Sbeta)*(3*
        MuSUSY*s2t*(127198*z3 + 6782*z4 - 14613*pow2(z2)) + 16*Mt*Tbeta*(-
        61328*z3 + 2114*z4 + 3171*pow2(z2)))*pow3(Mst2) - 8*Mst2*MuSUSY*(-2*Mt*
        pow2(Sbeta)*(-342281*z3 + 8792*z4 + 13188*pow2(z2)) - MuSUSY*s2t*Tbeta*
        (-1 + pow2(Sbeta))*(-323836*z3 - 8090*z4 + 37437*pow2(z2)))*pow3(Mt) -
        2*Mt*pow2(s2t)*pow2(Sbeta)*(MuSUSY*s2t*(547817*z3 + 10924*z4 - 6942*
        pow2(z2)) + 9*Mt*Tbeta*(-15053*z3 + 792*z4 + 1188*pow2(z2)))*pow4(Mst2)
        + Tbeta*(-48*pow2(MuSUSY)*(-1 + pow2(Sbeta))*(-25778*z3 + 338*z4 + 21*
        pow2(z2))*pow4(Mt) + pow2(Sbeta)*(Mst2*s2t*(89533*z3 - 4054*z4 - 5352*
        pow2(z2)) + 28*Mt*(-9920*z3 + 782*z4 + 1173*pow2(z2)))*pow3(s2t)*pow5(
        Mst2))) + 3*pow2(Mst2)*(-4*pow2(Mst2)*pow2(Mt)*(-(Tbeta*pow2(MuSUSY)*
        pow2(s2t)*(-1 + pow2(Sbeta))*(122917*z3 + 6344*z4 + 1254*pow2(z2))) +
        2*Mt*pow2(Sbeta)*(7*Mt*Tbeta*(-29483*z3 + 152*z4 + 228*pow2(z2)) + 3*
        MuSUSY*s2t*(55597*z3 - 264*z4 + 252*pow2(z2)))) + 4*Mt*s2t*pow2(Sbeta)*
        (4*Mt*(3*MuSUSY*s2t*(32773*z3 + 218*z4 - 3804*pow2(z2)) + 2*Mt*Tbeta*(-
        32323*z3 + 112*z4 + 168*pow2(z2))) + Mst2*s2t*(3*Mt*Tbeta*(55597*z3 -
        264*z4 + 252*pow2(z2)) + MuSUSY*s2t*(122917*z3 + 6344*z4 + 1254*pow2(
        z2))))*pow3(Mst2) + 32*Mst2*MuSUSY*(Mt*pow2(Sbeta)*(32323*z3 - 112*z4 -
        168*pow2(z2)) - MuSUSY*s2t*Tbeta*(-1 + pow2(Sbeta))*(-32773*z3 - 218*z4
        + 3804*pow2(z2)))*pow3(Mt) + Tbeta*(-576*(608*z3 - 9*z4)*pow2(MuSUSY)*(
        -1 + pow2(Sbeta))*pow4(Mt) + pow2(Sbeta)*(-(Mst2*s2t*(122917*z3 + 6344*
        z4 + 1254*pow2(z2))) + 16*Mt*(-32773*z3 - 218*z4 + 3804*pow2(z2)))*
        pow3(s2t)*pow5(Mst2)))) - 8*Dmglst2*(27*pow2(Mst2)*(4*Mst2*pow2(Mt)*(8*
        Tbeta*(-135*Dmsqst2*z3*pow2(MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta)) + 2*
        pow2(Mt)*((83*z3 - 27*z4)*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 900*
        Dmsqst2*z3*pow2(Sbeta))) + Mst2*(2160*Dmsqst2*MuSUSY*z3*pow2(s2t)*pow2(
        Sbeta) + 2*MuSUSY*pow2(Mt)*pow2(Sbeta)*(18*z3 + 506*z4 - 105*pow2(z2))
        + Mt*s2t*Tbeta*(-17280*Dmsqst2*z3*pow2(Sbeta) + pow2(MuSUSY)*(-1 +
        pow2(Sbeta))*(-30274*z3 - 542*z4 + 5289*pow2(z2))))) + 5760*Dmsqst2*
        MuSUSY*z3*(MuSUSY*s2t*Tbeta*(-1 + pow2(Sbeta)) + 12*Mt*pow2(Sbeta))*
        pow3(Mt) - 4*Mt*pow3(Mst2)*(Mt*Tbeta*pow2(MuSUSY)*pow2(s2t)*(-1 + pow2(
        Sbeta))*(17615*z3 + 424*z4 + 474*pow2(z2)) + pow2(Sbeta)*(8*MuSUSY*s2t*
        (-439*z3 + 108*z4)*pow2(Mt) + 2*Tbeta*(9017*z3 + 56*z4 + 84*pow2(z2))*
        pow3(Mt) + 1080*Dmsqst2*MuSUSY*z3*pow3(s2t))) + s2t*pow2(Sbeta)*(2*s2t*
        pow2(Mt)*(8*Mst2*Tbeta*(-439*z3 + 108*z4) + 3*MuSUSY*(-30274*z3 - 542*
        z4 + 5289*pow2(z2))) - 2*Mt*pow2(s2t)*(1440*Dmsqst2*Tbeta*z3 + 2*Mst2*
        MuSUSY*(17615*z3 + 424*z4 + 474*pow2(z2)) + Tbeta*pow2(Mst2)*(-30274*z3
        - 542*z4 + 5289*pow2(z2))) + 8*Tbeta*(-18*z3 - 506*z4 + 105*pow2(z2))*
        pow3(Mt) + Mst2*Tbeta*(1080*Dmsqst2*z3 + pow2(Mst2)*(17615*z3 + 424*z4
        + 474*pow2(z2)))*pow3(s2t))*pow4(Mst2)) - pow4(Mst1)*(8*s2t*(-(Tbeta*
        pow2(MuSUSY)*(-1 + pow2(Sbeta))*(-4*(395581*z3 + 659*z4) + 259215*pow2(
        z2))) + pow2(Sbeta)*(855360*Dmsqst2*Tbeta*z3 - 36*Mst2*MuSUSY*(2806*z3
        + 1514*z4 + 2757*pow2(z2)) - 4*Tbeta*pow2(Mst2)*(-653582*z3 + 25730*z4
        + 38595*pow2(z2))))*pow3(Mt) + pow2(Sbeta)*(2*Mt*(8*MuSUSY*(89947*z3 +
        6332*z4 + 9498*pow2(z2)) + Mst2*Tbeta*(-565214*z3 + 31142*z4 + 46713*
        pow2(z2)))*pow3(Mst2)*pow3(s2t) + 32*(2*Mst2*Tbeta*(-591386*z3 + 5582*
        z4 + 8373*pow2(z2)) + MuSUSY*(-834482*z3 + 28538*z4 + 48639*pow2(z2)))*
        pow4(Mt) + Tbeta*(48600*Dmsqst2*z3 + pow2(Mst2)*(217283*z3 - 16796*z4 -
        25194*pow2(z2)))*pow3(Mst2)*pow4(s2t)) - 4*pow2(Mt)*pow2(s2t)*(pow2(
        Sbeta)*(233280*Dmsqst2*MuSUSY*z3 + 3*MuSUSY*pow2(Mst2)*(-728116*z3 +
        10126*z4 + 105585*pow2(z2)) - 36*Tbeta*(-734*z3 + 1538*z4 + 2307*pow2(
        z2))*pow3(Mst2)) - 4*Mst2*Tbeta*(pow2(MuSUSY)*(-1 + pow2(Sbeta))*(
        353116*z3 + 11327*z4 + 15897*pow2(z2)) - 6075*Dmsqst2*z3*(-2 + pow2(
        Sbeta))*pow4(Sbeta)))) + 18*pow2(Mst1)*(8*pow2(Mst2)*pow2(Mt)*(6480*
        Dmsqst2*MuSUSY*z3*pow2(s2t)*pow2(Sbeta) - 24*MuSUSY*pow2(Mt)*pow2(
        Sbeta)*(-1675*z3 + 26*z4 + 93*pow2(z2)) + Mt*s2t*Tbeta*(-30240*Dmsqst2*
        z3*pow2(Sbeta) + pow2(MuSUSY)*(-1 + pow2(Sbeta))*(-47456*z3 - 709*z4 +
        8535*pow2(z2)))) + 43200*Dmsqst2*MuSUSY*z3*(MuSUSY*s2t*Tbeta*(-1 +
        pow2(Sbeta)) + 8*Mt*pow2(Sbeta))*pow3(Mt) + Mt*pow3(Mst2)*(12*s2t*pow2(
        Mt)*pow2(Sbeta)*(8*MuSUSY*(590*z3 - 4*z4 + 75*pow2(z2)) + Mst2*Tbeta*(-
        26782*z3 + 922*z4 + 1383*pow2(z2))) + 3*Mt*MuSUSY*pow2(s2t)*(-8*MuSUSY*
        Tbeta*(-1 + pow2(Sbeta))*(9747*z3 + 185*z4 + 237*pow2(z2)) + Mst2*pow2(
        Sbeta)*(-99002*z3 - 1210*z4 + 18273*pow2(z2))) - 16*Tbeta*pow2(Sbeta)*(
        -28846*z3 + 322*z4 + 483*pow2(z2))*pow3(Mt) - 2160*Dmsqst2*(5*MuSUSY +
        6*Mst2*Tbeta)*z3*pow2(Sbeta)*pow3(s2t)) + 6*pow2(s2t)*pow2(Sbeta)*(180*
        Dmsqst2*Tbeta*z3*pow2(s2t) - 4*Tbeta*pow2(Mt)*(741*z3 + 100*z4 + 150*
        pow2(z2)) - Mt*MuSUSY*s2t*(21373*z3 + 316*z4 + 474*pow2(z2)))*pow5(
        Mst2) + Tbeta*(48*Mst2*pow2(Mt)*(-360*Dmsqst2*z3*pow2(MuSUSY)*pow2(s2t)
        *(-1 + pow2(Sbeta)) + pow2(Mt)*(6480*Dmsqst2*z3*pow2(Sbeta) + pow2(
        MuSUSY)*(-1 + pow2(Sbeta))*(250*z3 - 94*z4 + 21*pow2(z2)))) + pow2(
        Sbeta)*(3*Mst2*s2t*(1879*z3 - 54*z4) + Mt*(8180*z3 - 416*z4 - 2406*
        pow2(z2)))*pow3(s2t)*pow6(Mst2)))) + 3*(-(pow4(Mst1)*(-32*s2t*pow2(Mt)*
        pow2(Sbeta)*(3*MuSUSY*s2t*(11816*z3 + 4954*z4 - 6177*pow2(z2)) + 8*Mt*
        Tbeta*(-57758*z3 + 1370*z4 + 2055*pow2(z2)))*pow3(Mst2) + 4*pow2(s2t)*
        pow2(Sbeta)*(4*Tbeta*pow2(Mt)*(-65326*z3 + 13714*z4 + 20571*pow2(z2)) +
        4*Mt*s2t*(Mst2*Tbeta*(-44630*z3 + 878*z4 + 1317*pow2(z2)) + 7*MuSUSY*(-
        43868*z3 + 1970*z4 + 2955*pow2(z2))) - Tbeta*pow2(s2t)*(40*Dmsqst2*(-
        260*z3 + 14*z4 + 21*pow2(z2)) + pow2(Mst2)*(4*(-5225 + 162*lmMst1 -
        162*lmMst2)*z3 + 2294*z4 + 5385*pow2(z2))))*pow4(Mst2) - pow2(Mst2)*(
        16*Tbeta*pow2(Mt)*pow2(s2t)*(-160*Dmsqst2*pow2(Sbeta)*(2*(-7*z3 + z4) +
        3*pow2(z2)) + pow2(MuSUSY)*(1 - pow2(Sbeta))*(-8*(91963 + 162*lmMst1 -
        162*lmMst2)*z3 + 27446*z4 + 35823*pow2(z2))) + pow2(Sbeta)*(64*MuSUSY*
        s2t*(-27086*z3 + 12062*z4 + 18093*pow2(z2))*pow3(Mt) - 33280*Dmsqst2*
        Mt*MuSUSY*(2*(-7*z3 + z4) + 3*pow2(z2))*pow3(s2t) + 1024*Tbeta*(24241*
        z3 + 146*z4 + 219*pow2(z2))*pow4(Mt) - 5*Tbeta*xDmsqst2*pow2(Dmsqst2)*(
        2453*z3 + 448*z4 + 672*pow2(z2))*pow4(s2t))) - 64*Mst2*Mt*(-(s2t*Tbeta*
        pow2(Mt)*(25920*Dmsqst2*z3*pow2(Sbeta) + pow2(MuSUSY)*(-1 + pow2(Sbeta)
        )*(-99632*z3 - 11764*z4 + 20505*pow2(z2)))) - 4*MuSUSY*pow2(Sbeta)*(-
        107072*z3 + 3326*z4 + 3045*pow2(z2))*pow3(Mt) + 405*Dmsqst2*z3*pow2(
        s2t)*pow2(Sbeta)*(96*Mt*MuSUSY - 5*Dmsqst2*s2t*Tbeta*xDmsqst2*(-2 - 2*
        pow2(Sbeta) + pow4(Sbeta)))) + 32*Mt*(-10*Dmsqst2*Mt*(-(Tbeta*pow2(
        MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta))*(-1771*z3 + 442*z4 + 663*pow2(z2))
        ) + 2*Mt*pow2(Sbeta)*(3888*Mt*Tbeta*z3 + MuSUSY*s2t*(-1930*z3 + 106*z4
        + 159*pow2(z2)))) - 6*Tbeta*pow2(MuSUSY)*(-1 + pow2(Sbeta))*(1502*z3 +
        2062*z4 + 1635*pow2(z2))*pow3(Mt) - 5*xDmsqst2*pow2(Dmsqst2)*pow2(s2t)*
        pow2(Sbeta)*(2*MuSUSY*s2t*(-3034*z3 + 94*z4 + 141*pow2(z2)) + 3*Mt*
        Tbeta*(20*z4 + 30*pow2(z2) + z3*(-140 - 162*pow2(Sbeta) + 81*pow4(
        Sbeta))))))) + 6*pow2(Mst1)*(34560*Dmsqst2*Mst2*z3*pow2(Mt)*(MuSUSY*(
        16*pow2(Mt) + 3*Dmsqst2*xDmsqst2*pow2(s2t))*pow2(Sbeta) - 2*Mt*s2t*
        Tbeta*(-7*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 4*Dmsqst2*xDmsqst2*pow2(
        Sbeta))) + 320*xDmsqst2*pow2(Dmsqst2)*pow2(Mt)*(1080*Tbeta*z3*pow2(Mt)*
        pow2(Sbeta) - 8*Mt*MuSUSY*s2t*pow2(Sbeta)*(2*(-7*z3 + z4) + 3*pow2(z2))
        + Tbeta*pow2(MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta))*(-944*z3 + 14*z4 +
        21*pow2(z2))) - 2*Mt*pow2(Mst2)*(-5*xDmsqst2*pow2(Dmsqst2)*pow2(s2t)*
        pow2(Sbeta)*(2*Mt*Tbeta*(7*z3 + 80*z4 + 120*pow2(z2)) + MuSUSY*s2t*(-
        16999*z3 + 352*z4 + 528*pow2(z2))) - 240*Dmsqst2*Mt*(-(Tbeta*pow2(
        MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta))*(-35*z3 + 26*z4 + 39*pow2(z2))) +
        2*Mt*pow2(Sbeta)*(288*Mt*Tbeta*z3 + MuSUSY*s2t*(-202*z3 + 10*z4 + 15*
        pow2(z2)))) - 288*Tbeta*pow2(MuSUSY)*(-1 + pow2(Sbeta))*(218*z3 + 50*z4
        + 21*pow2(z2))*pow3(Mt)) + 64*Mt*pow3(Mst2)*(-270*Dmsqst2*(-24*Mt*
        MuSUSY + Dmsqst2*s2t*Tbeta*xDmsqst2)*z3*pow2(s2t)*pow2(Sbeta) - s2t*
        Tbeta*pow2(Mt)*(4320*Dmsqst2*z3*pow2(Sbeta) + pow2(MuSUSY)*(-1 + pow2(
        Sbeta))*(-14636*z3 - 1135*z4 + 2388*pow2(z2))) - 4*MuSUSY*pow2(Sbeta)*(
        -8219*z3 + 326*z4 + 165*pow2(z2))*pow3(Mt)) - pow4(Mst2)*(-16*Tbeta*
        pow2(Mt)*pow2(s2t)*(-60*Dmsqst2*pow2(Sbeta)*(2*(-7*z3 + z4) + 3*pow2(
        z2)) - pow2(MuSUSY)*(-1 + pow2(Sbeta))*((-71438 - 216*lmMst1 + 216*
        lmMst2)*z3 + 2276*z4 + 2523*pow2(z2))) + pow2(Sbeta)*(-16*MuSUSY*s2t*(
        3718*z3 + 3470*z4 + 5205*pow2(z2))*pow3(Mt) + 960*Dmsqst2*Mt*MuSUSY*(-
        52*z3 + 10*z4 + 15*pow2(z2))*pow3(s2t) - 32*Tbeta*(69050*z3 + 418*z4 +
        627*pow2(z2))*pow4(Mt) + 5*Tbeta*xDmsqst2*pow2(Dmsqst2)*(-1895*z3 +
        128*z4 + 192*pow2(z2))*pow4(s2t))) + 8*s2t*pow2(Sbeta)*(-4*Mt*pow2(s2t)
        *(1080*Dmsqst2*Tbeta*z3 + Mst2*MuSUSY*((-23848 - 54*lmMst1 + 54*lmMst2)
        *z3 + 958*z4 + 1275*pow2(z2))) + s2t*pow2(Mt)*(MuSUSY*(68262*z3 + 9030*
        z4 - 13671*pow2(z2)) - 2*Mst2*Tbeta*(4*(-55 + 54*lmMst1 - 54*lmMst2)*z3
        + 826*z4 + 1887*pow2(z2))) + 4*Tbeta*(-48454*z3 + 754*z4 + 1131*pow2(
        z2))*pow3(Mt) + 15*Dmsqst2*Mst2*Tbeta*(-173*z3 + 14*z4 + 21*pow2(z2))*
        pow3(s2t))*pow5(Mst2) + 4*Tbeta*pow2(Sbeta)*(-4*Mt*(-6518*z3 + 740*z4 +
        219*pow2(z2)) + Mst2*s2t*(-23954*z3 + 1556*z4 + 2577*pow2(z2)))*pow3(
        s2t)*pow7(Mst2)) + 9*Mst2*(11520*Dmsqst2*z3*pow2(Mt)*(2*Dmsqst2*Mt*
        MuSUSY*xDmsqst2*(MuSUSY*s2t*Tbeta*(-1 + pow2(Sbeta)) + 4*Mt*pow2(Sbeta)
        ) + pow2(Mst2)*(MuSUSY*(16*pow2(Mt) + 3*Dmsqst2*xDmsqst2*pow2(s2t))*
        pow2(Sbeta) - 4*Mt*s2t*Tbeta*(-3*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 2*
        Dmsqst2*xDmsqst2*pow2(Sbeta)))) + 20*Mst2*xDmsqst2*pow2(Dmsqst2)*pow2(
        Mt)*(4*Mt*pow2(Sbeta)*(675*Mt*Tbeta*z3 + MuSUSY*s2t*(301*z3 - 16*z4 -
        24*pow2(z2))) + Tbeta*pow2(MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta))*(-4403*
        z3 + 32*z4 + 48*pow2(z2))) - 4*Mt*pow3(Mst2)*(-5*xDmsqst2*pow2(Dmsqst2)
        *pow2(s2t)*pow2(Sbeta)*(2*Mt*Tbeta*(-301*z3 + 16*z4 + 24*pow2(z2)) +
        MuSUSY*s2t*(-4403*z3 + 32*z4 + 48*pow2(z2))) - 240*Dmsqst2*Mt*(-(Tbeta*
        pow2(MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta))*(23*z3 + 2*z4 + 3*pow2(z2)))
        + 2*Mt*pow2(Sbeta)*(48*Mt*Tbeta*z3 + MuSUSY*s2t*(-58*z3 + 2*z4 + 3*
        pow2(z2)))) - 64*Tbeta*(230*z3 + 27*z4)*pow2(MuSUSY)*(-1 + pow2(Sbeta))
        *pow3(Mt)) - (96*Tbeta*pow2(Mt)*pow2(s2t)*(pow2(MuSUSY)*(1 - pow2(
        Sbeta))*(2*(1319 + 6*lmMst1 - 6*lmMst2)*z3 - 40*z4 + 3*pow2(z2)) + 10*
        Dmsqst2*pow2(Sbeta)*(-58*z3 + 2*z4 + 3*pow2(z2))) + pow2(Sbeta)*(-96*
        MuSUSY*s2t*(6*(77 - 8*lmMst1 + 8*lmMst2)*z3 + 202*z4 + 159*pow2(z2))*
        pow3(Mt) + 960*Dmsqst2*Mt*MuSUSY*(23*z3 + 2*z4 + 3*pow2(z2))*pow3(s2t)
        + Tbeta*(-192*(8*(514 - 3*(lmMst1 + lmMst2) + 6*lmMt)*z3 + 14*z4 + 21*
        pow2(z2))*pow4(Mt) + 5*xDmsqst2*pow2(Dmsqst2)*(-4403*z3 + 32*z4 + 48*
        pow2(z2))*pow4(s2t))))*pow5(Mst2) + Mt*(32*(-360*Dmsqst2*(-18*Mt*MuSUSY
        + Dmsqst2*s2t*Tbeta*xDmsqst2)*z3*pow2(s2t)*pow2(Sbeta) - 5*s2t*Tbeta*
        pow2(Mt)*(1152*Dmsqst2*z3*pow2(Sbeta) + pow2(MuSUSY)*(-1 + pow2(Sbeta))
        *(-2386*z3 - 102*z4 + 333*pow2(z2))) - 6*MuSUSY*pow2(Sbeta)*(-1922*z3 +
        206*z4 + 21*pow2(z2))*pow3(Mt))*pow4(Mst2) - 48*s2t*pow2(Sbeta)*(1440*
        Dmsqst2*Tbeta*z3*pow2(s2t) - 4*Tbeta*pow2(Mt)*(-1922*z3 + 206*z4 + 21*
        pow2(z2)) + 5*Mt*MuSUSY*s2t*(-2386*z3 - 102*z4 + 333*pow2(z2)))*pow6(
        Mst2)) + 48*pow2(s2t)*pow2(Sbeta)*(2*Mt*MuSUSY*s2t*(2*(1319 + 6*lmMst1
        - 6*lmMst2)*z3 - 40*z4 + 3*pow2(z2)) + Tbeta*(5*Dmsqst2*pow2(s2t)*(23*
        z3 + 2*z4 + 3*pow2(z2)) - pow2(Mt)*(6*(77 - 8*lmMst1 + 8*lmMst2)*z3 +
        202*z4 + 159*pow2(z2))))*pow7(Mst2) + 8*Tbeta*pow2(Sbeta)*(-3*Mst2*s2t*
        (2*(1319 + 6*lmMst1 - 6*lmMst2)*z3 - 40*z4 + 3*pow2(z2)) + 10*Mt*(-
        2386*z3 - 102*z4 + 333*pow2(z2)))*pow3(s2t)*pow8(Mst2)))))/(23328.*
        Tbeta*pow2(Sbeta)*pow6(Mst2)) - Al4p*((twoLoopFlag*((-72*Mt*pow3(s2t)*(
        -(Mst2*(4*(3 + lmMst1*(2 + lmMst2) - pow2(lmMst2))*pow2(Mst1)*pow2(
        Mst2) + (3 + 2*lmMst1 - 2*lmMst2)*pow4(Mst1) + 4*(-4 + lmMst1 - 3*
        lmMst2 + lmMst1*lmMst2 - pow2(lmMst2))*pow4(Mst2))) + Dmglst2*(-4*(1 +
        lmMst1*(-2 + lmMst2) + 4*lmMst2 - pow2(lmMst2))*pow2(Mst1)*pow2(Mst2) +
        (1 + 6*lmMst1 - 6*lmMst2)*pow4(Mst1) + 4*(6 + lmMst1 + lmMst2 - lmMst1*
        lmMst2 + pow2(lmMst2))*pow4(Mst2))))/pow2(Mst2) - (144*s2t*pow3(Mt)*(4*
        pow2(Mst1)*((4 - 5*lmMst2 + lmMst1*(5 + 2*lmMst2 - 2*lmMt) + 2*lmMst2*
        lmMt - 2*pow2(lmMst2))*pow2(Mst2) + (-5 - 6*lmMst2 + lmMst1*(4 + 3*
        lmMst2) - 3*pow2(lmMst2))*pow2(MuSUSY))*pow3(Mst2) + Mst2*(4*(5 - 7*
        lmMst2 + lmMst1*(7 + 2*lmMst2 - 2*lmMt) + 2*lmMst2*lmMt - 2*pow2(
        lmMst2))*pow2(Mst2) + (-21 - 38*lmMst2 + 10*lmMst1*(3 + 2*lmMst2) - 20*
        pow2(lmMst2))*pow2(MuSUSY))*pow4(Mst1) + 4*(-4 + lmMst1 - 3*lmMst2 +
        lmMst1*lmMst2 - pow2(lmMst2))*pow2(MuSUSY)*pow5(Mst2) + Dmglst2*(4*
        pow2(Mst1)*pow2(Mst2)*((6 + 13*lmMst2 - 4*lmMt - 2*lmMst2*lmMt +
        lmMst1*(-9 - 2*lmMst2 + 2*lmMt) + 2*pow2(lmMst2))*pow2(Mst2) + (-11 +
        2*lmMst2 + lmMst1*(-4 + 3*lmMst2) - 3*pow2(lmMst2))*pow2(MuSUSY)) - (4*
        (1 - 29*lmMst2 + lmMst1*(25 + 6*lmMst2 - 6*lmMt) + 4*lmMt + 6*lmMst2*
        lmMt - 6*pow2(lmMst2))*pow2(Mst2) + (65 + lmMst1*(34 - 20*lmMst2) - 26*
        lmMst2 + 20*pow2(lmMst2))*pow2(MuSUSY))*pow4(Mst1) - 4*(6 + lmMst1 +
        lmMst2 - lmMst1*lmMst2 + pow2(lmMst2))*pow2(MuSUSY)*pow4(Mst2) + 8*(5 +
        lmMst1*(-1 + lmMst2) + lmMst2 - pow2(lmMst2))*pow6(Mst2)) + 4*(3 - 4*
        lmMst2 + 2*(lmMst1 + lmMst1*lmMst2 + lmMt) - 2*pow2(lmMst2))*pow7(Mst2)
        ))/pow6(Mst2) + (32*pow4(Mt)*(2*Dmglst2*((53 - 18*lmMst1 + 24*lmMst2 -
        24*lmMt)*pow2(Mst1)*pow4(Mst2) + 18*((5 + lmMst2*(11 - 2*lmMt) - 2*
        lmMst1*(4 + lmMst2 - lmMt) - 3*lmMt + 2*pow2(lmMst2))*pow2(Mst2)*pow4(
        Mst1) + (3 - 6*lmMst2*(-5 + lmMt) - 5*lmMt + lmMst1*(-25 - 6*lmMst2 +
        6*lmMt) + 6*pow2(lmMst2))*pow6(Mst1)) - 18*lmMst2*pow6(Mst2)) + 9*(2*(2
        + 2*lmMst1*(3 + lmMst2 - lmMt) + lmMt + lmMst2*(-7 + 2*lmMt) - 2*pow2(
        lmMst2))*pow3(Mst2)*pow4(Mst1) + (3 + 2*lmMst1*(3 + lmMst2) + lmMt +
        lmMst2*(-5 + 4*lmMt) + pow2(lmMst1) - pow2(lmMst2) - 6*pow2(lmMt))*
        pow2(Mst1)*pow5(Mst2) + Mst2*(9 + lmMst1*(22 + 6*lmMst2 - 6*lmMt) + 6*
        lmMst2*(-4 + lmMt) + 2*lmMt - 6*pow2(lmMst2))*pow6(Mst1) - 2*(1 +
        lmMst2)*pow7(Mst2))))/(pow2(Mst1)*pow5(Mst2)) + (36*Mt*MuSUSY*(16*pow3(
        Mt)*(2*(4 + lmMst1*(3 + 2*lmMst2 - lmMt) + lmMst2*(-4 + lmMt) + lmMt -
        2*pow2(lmMst2))*pow2(Mst1)*pow3(Mst2) + Mst2*(13 + 2*lmMst1*(6 + 3*
        lmMst2 - 2*lmMt) + 2*lmMt + 2*lmMst2*(-7 + 2*lmMt) - 6*pow2(lmMst2))*
        pow4(Mst1) + Dmglst2*(11 - 6*lmMst1*lmMst2 - 8*lmMst2*(-5 + lmMt) + 8*
        lmMst1*(-4 + lmMt) - 8*lmMt + 6*pow2(lmMst2))*pow4(Mst1) + 2*Dmglst2*((
        7 + 7*lmMst2 + lmMst1*(-5 + lmMt) - 2*lmMt - lmMst2*lmMt)*pow2(Mst1)*
        pow2(Mst2) + (5 + lmMst1*(-1 + lmMst2) + lmMst2 - pow2(lmMst2))*pow4(
        Mst2)) + 2*(2 + lmMst1 + (-2 + lmMst1)*lmMst2 + lmMt - pow2(lmMst2))*
        pow5(Mst2)) + 6*Mt*pow2(Mst2)*pow2(s2t)*(4*(1 + 3*lmMst2 - lmMst1*(3 +
        2*lmMst2) + 2*pow2(lmMst2))*pow2(Mst1)*pow3(Mst2) + Mst2*(1 + 14*lmMst2
        - 2*lmMst1*(7 + 4*lmMst2) + 8*pow2(lmMst2))*pow4(Mst1) + Dmglst2*(4*(5
        + lmMst1*(3 - 2*lmMst2) - 3*lmMst2 + 2*pow2(lmMst2))*pow2(Mst1)*pow2(
        Mst2) + (21 + lmMst1*(18 - 8*lmMst2) - 18*lmMst2 + 8*pow2(lmMst2))*
        pow4(Mst1) + 4*(6 + lmMst1 + lmMst2 - lmMst1*lmMst2 + pow2(lmMst2))*
        pow4(Mst2)) + 4*(4 + 3*lmMst2 - lmMst1*(1 + lmMst2) + pow2(lmMst2))*
        pow5(Mst2)) + (8*Mst2*s2t*pow2(Mt)*(Dmglst2*(4*(1 - lmMst1 + lmMst2)*
        pow2(Mst2)*pow4(Mst1) + 2*(1 + 2*lmMst2)*pow2(Mst1)*pow4(Mst2) + (4 -
        8*lmMst1 + 8*lmMst2)*pow6(Mst1) - 4*lmMst2*pow6(Mst2)) - Mst2*(-2*
        lmMst1*((1 + 2*lmMst2)*(pow2(Mst1) + pow2(Mst2))*pow4(Mst1) + (3 +
        lmMst2)*pow2(Mst1)*pow4(Mst2)) + pow2(lmMst1)*(2*pow2(Mst2)*pow4(Mst1)
        - pow2(Mst1)*pow4(Mst2) + 2*pow6(Mst1)) + pow2(lmMst2)*(2*pow2(Mst2)*
        pow4(Mst1) + 3*pow2(Mst1)*pow4(Mst2) + 2*pow6(Mst1)) + 2*pow6(Mst2) +
        2*lmMst2*(pow2(Mst2)*pow4(Mst1) + 2*pow2(Mst1)*pow4(Mst2) + pow6(Mst1)
        + pow6(Mst2)))))/pow2(Mst1) + (pow3(Mst2)*pow3(s2t)*(2*(-16 + 6*lmMst2
        - 2*lmMst1*(8 + 5*lmMst2) + 3*pow2(lmMst1) + 7*pow2(lmMst2))*pow3(Mst2)
        *pow4(Mst1) - 2*(-12 - 18*lmMst2 + 2*lmMst1*(5 + 3*lmMst2) + pow2(
        lmMst1) - 7*pow2(lmMst2))*pow2(Mst1)*pow5(Mst2) + Mst2*(3 + lmMst1*(2 -
        32*lmMst2) - 2*lmMst2 + 16*pow2(lmMst1) + 16*pow2(lmMst2))*pow6(Mst1) -
        4*Dmglst2*(2*(1 + 2*lmMst1 - lmMst2)*pow2(Mst2)*pow4(Mst1) + 2*(1 +
        lmMst1 - lmMst2 + 2*lmMst1*lmMst2 - 2*pow2(lmMst2))*pow2(Mst1)*pow4(
        Mst2) + (1 + 2*lmMst1 - 2*lmMst2)*pow6(Mst1) - 2*lmMst2*pow6(Mst2)) +
        4*(1 + lmMst2)*pow7(Mst2)))/pow2(Mst1)))/(Tbeta*pow6(Mst2)) + 36*pow2(
        Mt)*pow2(s2t)*(4*(-2*lmMst1*(-2 + lmMst2) - lmMst2*(2 + lmMst2) + 3*
        pow2(lmMst1))*pow2(Mst1) - 4*(2 - 2*lmMst2 + 2*lmMst1*(3 + lmMst2) +
        pow2(lmMst1) - 3*pow2(lmMst2))*pow2(Mst2) + (8*(2*Dmglst2*lmMst2 + Mst2
        + lmMst2*Mst2)*pow3(Mst2))/pow2(Mst1) - (8*Dmglst2*(pow2(Mst1)*pow2(
        Mst2) - 2*lmMst1*pow2(Mst1)*(pow2(Mst1) + pow2(Mst2)) + 2*lmMst2*pow4(
        Mst1) + pow4(Mst2) + 4*lmMst2*pow4(Mst2)))/pow3(Mst2) + (pow2(MuSUSY)*(
        4*(-1 - 13*lmMst1 + 13*lmMst2 - 8*lmMst1*lmMst2 + pow2(lmMst1) + 7*
        pow2(lmMst2))*pow3(Mst2)*pow4(Mst1) - 2*(-14 + 10*lmMst1 - 20*lmMst2 +
        6*lmMst1*lmMst2 + pow2(lmMst1) - 7*pow2(lmMst2))*pow2(Mst1)*pow5(Mst2)
        + Mst2*(-1 - 50*lmMst1 + 50*lmMst2 - 64*lmMst1*lmMst2 + 20*pow2(lmMst1)
        + 44*pow2(lmMst2))*pow6(Mst1) - 4*Dmglst2*(2*(2 - 3*lmMst2 + lmMst1*(3
        + 2*lmMst2) - 2*pow2(lmMst2))*pow2(Mst2)*pow4(Mst1) + 2*(1 + lmMst1 -
        2*lmMst2 + 2*lmMst1*lmMst2 - 2*pow2(lmMst2))*pow2(Mst1)*pow4(Mst2) + (5
        - 8*lmMst2 + 4*lmMst1*(2 + lmMst2) - 4*pow2(lmMst2))*pow6(Mst1) - 2*
        lmMst2*pow6(Mst2)) + 4*(1 + lmMst2)*pow7(Mst2)))/(pow2(Mst1)*pow5(Mst2)
        )) + ((-9*pow4(s2t)*(4*(-14 - 3*lmMst1 - 6*lmMst2 - 2*lmMst1*lmMst2 +
        2*pow2(lmMst1))*pow3(Mst2)*pow4(Mst1) - 2*(-10 + 10*lmMst1 - 16*lmMst2
        + 6*lmMst1*lmMst2 + pow2(lmMst1) - 7*pow2(lmMst2))*pow2(Mst1)*pow5(
        Mst2) + Mst2*(35 + 34*lmMst1 - 14*lmMst2 - 12*lmMst1*lmMst2 + 10*pow2(
        lmMst1) + 2*pow2(lmMst2))*pow6(Mst1) + 4*Dmglst2*(-2*(lmMst1 - 2*
        lmMst1*lmMst2 + 2*pow2(lmMst2))*pow2(Mst2)*pow4(Mst1) - 2*(1 + lmMst1 +
        2*lmMst1*lmMst2 - 2*pow2(lmMst2))*pow2(Mst1)*pow4(Mst2) + (1 + 2*
        lmMst1)*pow6(Mst1) + 2*lmMst2*pow6(Mst2)) + 4*(1 + lmMst2)*pow7(Mst2)))
        /Mst2 - (36*s2t*pow2(Mt)*pow2(MuSUSY)*(4*(4*Mt*(5 + 6*lmMst2 - lmMst1*(
        4 + 3*lmMst2) + 3*pow2(lmMst2)) + Mst2*s2t*(-1 + 13*lmMst2 - lmMst1*(13
        + 8*lmMst2) + pow2(lmMst1) + 7*pow2(lmMst2)))*pow3(Mst2)*pow4(Mst1) +
        2*(-(Mst2*s2t*(-14 - 20*lmMst2 + 2*lmMst1*(5 + 3*lmMst2) + pow2(lmMst1)
        - 7*pow2(lmMst2))) + 8*Mt*(4 + 3*lmMst2 - lmMst1*(1 + lmMst2) + pow2(
        lmMst2)))*pow2(Mst1)*pow5(Mst2) + 4*Dmglst2*(Mst2*s2t*(-5 + 8*lmMst2 -
        4*lmMst1*(2 + lmMst2) + 4*pow2(lmMst2)) + Mt*(65 + lmMst1*(34 - 20*
        lmMst2) - 26*lmMst2 + 20*pow2(lmMst2)))*pow6(Mst1) + Mst2*(Mst2*s2t*(-1
        + 50*lmMst2 - 2*lmMst1*(25 + 32*lmMst2) + 20*pow2(lmMst1) + 44*pow2(
        lmMst2)) + Mt*(84 + 152*lmMst2 - 40*lmMst1*(3 + 2*lmMst2) + 80*pow2(
        lmMst2)))*pow6(Mst1) + 8*Dmglst2*((Mst2*s2t*(-2 + 3*lmMst2 - lmMst1*(3
        + 2*lmMst2) + 2*pow2(lmMst2)) + Mt*(22 + lmMst1*(8 - 6*lmMst2) - 4*
        lmMst2 + 6*pow2(lmMst2)))*pow2(Mst2)*pow4(Mst1) + (-(Mst2*s2t*(1 +
        lmMst1 - 2*lmMst2 + 2*lmMst1*lmMst2 - 2*pow2(lmMst2))) + 2*Mt*(6 +
        lmMst1 + lmMst2 - lmMst1*lmMst2 + pow2(lmMst2)))*pow2(Mst1)*pow4(Mst2)
        + lmMst2*s2t*pow7(Mst2)) + 4*(1 + lmMst2)*s2t*pow8(Mst2)))/(pow2(Sbeta)
        *pow6(Mst2)))/pow2(Mst1)))/216. + (xDR2DRMOD*((3375*Al4p*threeLoopFlag*
        xDmsqst2*pow2(Dmsqst2)*pow2(Mst1)*(-(pow2(s2t)*pow4(Mst1)*(-8*(lmMst1 -
        lmMst2)*Tbeta*pow2(Mt)*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 4*Mt*(-(
        MuSUSY*s2t) + 2*Mt*Tbeta)*pow2(Mst2)*pow2(Sbeta) + (1 - 2*lmMst1 + 2*
        lmMst2)*Tbeta*pow2(s2t)*pow2(Sbeta)*pow4(Mst2))) + pow4(Mst2)*(-4*
        Tbeta*pow2(Mt)*pow2(s2t)*(pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 2*pow2(
        Mst2)*pow2(Sbeta)) + pow2(Sbeta)*(16*MuSUSY*s2t*pow3(Mt) - 4*Mt*MuSUSY*
        pow2(Mst2)*pow3(s2t) + 16*Tbeta*pow4(Mt) + Tbeta*pow4(Mst2)*pow4(s2t)))
        - pow2(Mst1)*pow2(Mst2)*(-4*Tbeta*pow2(Mt)*pow2(s2t)*((-1 + 2*lmMst1 -
        2*lmMst2)*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 4*pow2(Mst2)*pow2(Sbeta)) +
        pow2(Sbeta)*(16*MuSUSY*s2t*pow3(Mt) + 8*(-lmMst1 + lmMst2)*Mt*MuSUSY*
        pow2(Mst2)*pow3(s2t) - 16*Tbeta*pow4(Mt) + (1 + 2*lmMst1 - 2*lmMst2)*
        Tbeta*pow4(Mst2)*pow4(s2t))) + 19*xDmglst2*pow2(Dmglst2)*(pow2(s2t)*(4*
        Mt*MuSUSY*s2t - 8*Tbeta*pow2(Mt) + (-1 + 2*lmMst1 - 2*lmMst2)*Tbeta*
        pow2(Mst2)*pow2(s2t))*pow2(Sbeta)*pow4(Mst1) + pow2(Mst2)*(-4*Tbeta*
        pow2(Mt)*pow2(s2t)*(pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 2*pow2(Mst2)*
        pow2(Sbeta)) + pow2(Sbeta)*(16*MuSUSY*s2t*pow3(Mt) - 4*Mt*MuSUSY*pow2(
        Mst2)*pow3(s2t) + 16*Tbeta*pow4(Mt) + Tbeta*pow4(Mst2)*pow4(s2t))) -
        pow2(Mst1)*(-4*Tbeta*pow2(Mt)*pow2(s2t)*((-1 + 2*lmMst1 - 2*lmMst2)*
        pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 4*pow2(Mst2)*pow2(Sbeta)) + pow2(
        Sbeta)*(16*MuSUSY*s2t*pow3(Mt) + 8*(-lmMst1 + lmMst2)*Mt*MuSUSY*pow2(
        Mst2)*pow3(s2t) - 16*Tbeta*pow4(Mt) + (1 + 2*lmMst1 - 2*lmMst2)*Tbeta*
        pow4(Mst2)*pow4(s2t)))) + 14*Dmglst2*Mst2*(pow2(s2t)*(4*Mt*MuSUSY*s2t -
        8*Tbeta*pow2(Mt) + Tbeta*(pow2(Mst1) + (-1 + 2*lmMst1 - 2*lmMst2)*pow2(
        Mst2))*pow2(s2t))*pow2(Sbeta)*pow4(Mst1) + pow2(Mst2)*(-4*Tbeta*pow2(
        Mt)*pow2(s2t)*(pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 2*pow2(Mst2)*pow2(
        Sbeta)) + pow2(Sbeta)*(16*MuSUSY*s2t*pow3(Mt) - 4*Mt*MuSUSY*pow2(Mst2)*
        pow3(s2t) + 16*Tbeta*pow4(Mt) + Tbeta*pow4(Mst2)*pow4(s2t))) - pow2(
        Mst1)*(-4*Tbeta*pow2(Mt)*pow2(s2t)*((-1 + 2*lmMst1 - 2*lmMst2)*pow2(
        MuSUSY)*(-1 + pow2(Sbeta)) + 4*pow2(Mst2)*pow2(Sbeta)) + pow2(Sbeta)*(
        16*MuSUSY*s2t*pow3(Mt) + 8*(-lmMst1 + lmMst2)*Mt*MuSUSY*pow2(Mst2)*
        pow3(s2t) - 16*Tbeta*pow4(Mt) + (1 + 2*lmMst1 - 2*lmMst2)*Tbeta*pow4(
        Mst2)*pow4(s2t)))) + Tbeta*pow2(Mst2)*pow2(Sbeta)*pow4(s2t)*pow6(Mst1))
        )/(Tbeta*pow2(Sbeta)) + (2700*Mst2*(2*Dmglst2*lmMst2 + Mst2 + lmMst2*
        Mst2)*twoLoopFlag*pow2(Mst1)*(-(pow2(Mst2)*pow2(s2t)*pow4(Mst1)*(-8*(
        lmMst1 - lmMst2)*Tbeta*pow2(Mt)*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 4*Mt*
        (-(MuSUSY*s2t) + 2*Mt*Tbeta)*pow2(Mst2)*pow2(Sbeta) + (1 - 2*lmMst1 +
        2*lmMst2)*Tbeta*pow2(s2t)*pow2(Sbeta)*pow4(Mst2))) - pow2(Mst1)*pow4(
        Mst2)*(-4*Tbeta*pow2(Mt)*pow2(s2t)*((-1 + 2*lmMst1 - 2*lmMst2)*pow2(
        MuSUSY)*(-1 + pow2(Sbeta)) + 4*pow2(Mst2)*pow2(Sbeta)) + pow2(Sbeta)*(
        16*MuSUSY*s2t*pow3(Mt) + 8*(-lmMst1 + lmMst2)*Mt*MuSUSY*pow2(Mst2)*
        pow3(s2t) - 16*Tbeta*pow4(Mt) + (1 + 2*lmMst1 - 2*lmMst2)*Tbeta*pow4(
        Mst2)*pow4(s2t))) + Tbeta*pow2(s2t)*(8*(lmMst1 - lmMst2)*pow2(Mt)*pow2(
        MuSUSY)*(-1 + pow2(Sbeta)) + pow2(s2t)*pow2(Sbeta)*pow4(Mst2))*pow6(
        Mst1) + (-4*Tbeta*pow2(Mt)*pow2(s2t)*(pow2(MuSUSY)*(-1 + pow2(Sbeta)) +
        2*pow2(Mst2)*pow2(Sbeta)) + pow2(Sbeta)*(16*MuSUSY*s2t*pow3(Mt) - 4*Mt*
        MuSUSY*pow2(Mst2)*pow3(s2t) + 16*Tbeta*pow4(Mt) + Tbeta*pow4(Mst2)*
        pow4(s2t)))*pow6(Mst2)))/(Tbeta*pow2(Sbeta)) + 2*xDmglst2*pow2(Dmglst2)
        *((-1350*(2 - lmMst2)*twoLoopFlag*pow2(Mst1)*(-(pow2(Mst2)*pow2(s2t)*
        pow4(Mst1)*(-8*(lmMst1 - lmMst2)*Tbeta*pow2(Mt)*pow2(MuSUSY)*(-1 +
        pow2(Sbeta)) + 4*Mt*(-(MuSUSY*s2t) + 2*Mt*Tbeta)*pow2(Mst2)*pow2(Sbeta)
        + (1 - 2*lmMst1 + 2*lmMst2)*Tbeta*pow2(s2t)*pow2(Sbeta)*pow4(Mst2))) -
        pow2(Mst1)*pow4(Mst2)*(-4*Tbeta*pow2(Mt)*pow2(s2t)*((-1 + 2*lmMst1 - 2*
        lmMst2)*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 4*pow2(Mst2)*pow2(Sbeta)) +
        pow2(Sbeta)*(16*MuSUSY*s2t*pow3(Mt) + 8*(-lmMst1 + lmMst2)*Mt*MuSUSY*
        pow2(Mst2)*pow3(s2t) - 16*Tbeta*pow4(Mt) + (1 + 2*lmMst1 - 2*lmMst2)*
        Tbeta*pow4(Mst2)*pow4(s2t))) + Tbeta*pow2(s2t)*(8*(lmMst1 - lmMst2)*
        pow2(Mt)*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + pow2(s2t)*pow2(Sbeta)*pow4(
        Mst2))*pow6(Mst1) + (-4*Tbeta*pow2(Mt)*pow2(s2t)*(pow2(MuSUSY)*(-1 +
        pow2(Sbeta)) + 2*pow2(Mst2)*pow2(Sbeta)) + pow2(Sbeta)*(16*MuSUSY*s2t*
        pow3(Mt) - 4*Mt*MuSUSY*pow2(Mst2)*pow3(s2t) + 16*Tbeta*pow4(Mt) +
        Tbeta*pow4(Mst2)*pow4(s2t)))*pow6(Mst2)))/(Tbeta*pow2(Sbeta)) + Al4p*
        threeLoopFlag*(6400*Mst2*s2t*pow2(Mst1)*pow3(Mt)*(-3*(6*pow2(Mst2)*(116
        + lmMst2*(181 - 30*lmMt) - 36*lmMt - 4*(-14 + lmMt)*pow2(lmMst2) - 2*
        lmMst1*(42 + lmMst2*(21 - 2*lmMt) - 8*lmMt + 2*pow2(lmMst2)) + 4*pow3(
        lmMst2)) + pow2(MuSUSY)*(85 - 215*lmMst2 + lmMst1*(66 + 90*lmMst2 - 72*
        pow2(lmMst2)) - 54*pow2(lmMst2) + 72*pow3(lmMst2)))*pow4(Mst1) + 18*(3
        + 4*lmMst2 - 3*pow2(lmMst2))*(2*pow2(Mst2) + pow2(MuSUSY))*pow4(Mst2) -
        pow2(Mst1)*(-3*pow2(Mst2)*pow2(MuSUSY)*(5 + 149*lmMst2 + 12*pow2(
        lmMst2) + 6*lmMst1*(-7 - 8*lmMst2 + 6*pow2(lmMst2)) - 36*pow3(lmMst2))
        + (553 + 194*lmMst2 + 18*lmMst1*(-1 + 4*lmMst2 - 2*lmMt) - 30*lmMt +
        48*lmMst2*lmMt - 192*pow2(lmMst2))*pow4(Mst2))) + 9600*Mt*pow2(Mst1)*
        pow3(s2t)*(pow2(Mst1)*pow2(Mst2)*(31 - 101*lmMst2 + 6*lmMst1*(7 + 8*
        lmMst2 - 6*pow2(lmMst2)) - 48*pow2(lmMst2) + 36*pow3(lmMst2)) + (77 +
        59*lmMst2 - 6*lmMst1*(3 + lmMst2) - 12*pow2(lmMst2))*pow4(Mst1) + 6*(-3
        - 4*lmMst2 + 3*pow2(lmMst2))*pow4(Mst2))*pow5(Mst2) + 16*pow2(Mst2)*
        pow4(Mt)*(-43875*Dmsqst2*pow2(Mst1)*(pow2(Mst1) + pow2(Mst2)) + (394748
        - 14400*lmMst1*(8 + 2*lmMst2 - lmMt) - 102480*lmMt - 6*lmMst2*(-50363 +
        4680*lmMt) + 41355*pow2(lmMst2))*pow2(Mst2)*pow4(Mst1) + 225*(-188 +
        32*lmMst1*(-2 + lmMst2) - 166*lmMst2 + 59*pow2(lmMst2))*pow2(Mst1)*
        pow4(Mst2) + 7200*(283 + lmMst2*(489 - 93*lmMt) - 94*lmMt + 9*(21 - 2*
        lmMt)*pow2(lmMst2) - 2*lmMst1*(115 - 9*lmMst2*(-8 + lmMt) - 24*lmMt +
        9*pow2(lmMst2)) + 18*pow3(lmMst2))*pow6(Mst1) + 7200*(2 + lmMst2 - 3*
        pow2(lmMst2))*pow6(Mst2)) + (100*Mst2*Mt*MuSUSY*(2*Mst2*(3510*Dmsqst2*
        s2t*pow2(Mt) - 144*Mt*pow2(s2t)*(13 - 125*lmMst2 + 6*lmMst1*(7 + 8*
        lmMst2 - 6*pow2(lmMst2)) - 30*pow2(lmMst2) + 36*pow3(lmMst2))*pow3(
        Mst2) + 64*Mst2*(215 + 37*lmMst2 + 3*(-5 + 8*lmMst2)*lmMt - 18*lmMst1*(
        1 - 2*lmMst2 + lmMt) - 42*pow2(lmMst2))*pow3(Mt) + 3*pow2(Mst2)*(2*s2t*
        pow2(Mt)*(964 + 306*lmMst2 + 96*(-2 + lmMst2)*pow2(lmMst1) - 81*pow2(
        lmMst2) + 96*(lmMst1 + 3*lmMst1*lmMst2 - 2*lmMst1*pow2(lmMst2) + pow3(
        lmMst2))) - 585*Dmsqst2*(lmMst1 - lmMst2)*pow3(s2t) + pow2(Mst2)*(200 +
        1196*lmMst2 - 48*(-2 + lmMst2)*pow2(lmMst1) + 306*pow2(lmMst2) + 3*
        lmMst1*(-492 + 10*lmMst2 + 123*pow2(lmMst2)) - 321*pow3(lmMst2))*pow3(
        s2t)))*pow4(Mst1) - 9*pow2(Mst1)*pow3(Mst2)*(780*Dmsqst2*s2t*pow2(Mt) +
        192*Mt*(-3 - 4*lmMst2 + 3*pow2(lmMst2))*pow2(s2t)*pow3(Mst2) + 256*
        Mst2*(3 + 4*lmMst2 - 3*pow2(lmMst2))*pow3(Mt) + pow2(Mst2)*(4*s2t*(220
        - 32*lmMst1*(-2 + lmMst2) + 70*lmMst2 - 59*pow2(lmMst2))*pow2(Mt) -
        195*Dmsqst2*pow3(s2t)) + (-444 + 32*lmMst1*(-2 + lmMst2) - 70*lmMst2 +
        347*pow2(lmMst2))*pow3(s2t)*pow4(Mst2)) + (1152*Mst2*s2t*pow2(Mt)*(26 +
        29*lmMst2 + (-2 + lmMst2)*pow2(lmMst1) + 4*pow2(lmMst2) - 2*lmMst1*(8 +
        lmMst2 + pow2(lmMst2)) + pow3(lmMst2)) - 1728*Mt*pow2(Mst2)*pow2(s2t)*(
        15 - 11*lmMst2 + lmMst1*(4 + 7*lmMst2 - 6*pow2(lmMst2)) - 7*pow2(
        lmMst2) + 6*pow3(lmMst2)) + 128*(1097 + lmMst2*(1531 - 246*lmMt) - 339*
        lmMt + (444 - 36*lmMt)*pow2(lmMst2) - 18*lmMst1*(39 - 2*lmMst2*(-9 +
        lmMt) - 7*lmMt + 2*pow2(lmMst2)) + 36*pow3(lmMst2))*pow3(Mt) - 3*Mst2*(
        585*Dmsqst2 + pow2(Mst2)*(-476 - 1790*lmMst2 + 768*(-2 + lmMst2)*pow2(
        lmMst1) - 1617*pow2(lmMst2) - 96*lmMst1*(-13 - 33*lmMst2 + 16*pow2(
        lmMst2)) + 768*pow3(lmMst2)))*pow3(s2t))*pow6(Mst1) + 288*s2t*(2 +
        lmMst2 - 3*pow2(lmMst2))*(4*pow2(Mt) - pow2(Mst2)*pow2(s2t))*pow7(Mst2)
        ))/Tbeta - 75*pow4(Mst2)*pow4(s2t)*((-932 + 2182*lmMst2 - 96*(-2 +
        lmMst2)*pow2(lmMst1) + 1653*pow2(lmMst2) + 6*lmMst1*(-524 + 26*lmMst2 +
        123*pow2(lmMst2)) - 642*pow3(lmMst2))*pow4(Mst1)*pow4(Mst2) + 585*
        Dmsqst2*pow2(Mst1)*(-((1 + 2*lmMst1 - 2*lmMst2)*pow2(Mst1)*pow2(Mst2))
        + (-1 + 2*lmMst1 - 2*lmMst2)*pow4(Mst1) + pow4(Mst2)) + pow2(Mst2)*(76
        - 602*lmMst2 - 672*(-2 + lmMst2)*pow2(lmMst1) + 1005*pow2(lmMst2) + 6*
        lmMst1*(284 - 538*lmMst2 + 133*pow2(lmMst2)) - 126*pow3(lmMst2))*pow6(
        Mst1) - 3*(-508 + 32*lmMst1*(-2 + lmMst2) - 102*lmMst2 + 443*pow2(
        lmMst2))*pow2(Mst1)*pow6(Mst2) + 96*(-2 - lmMst2 + 3*pow2(lmMst2))*
        pow8(Mst2)) - (300*s2t*pow2(Mt)*pow2(MuSUSY)*(pow2(Mst2)*(-585*Dmsqst2*
        (-1 + 2*lmMst1 - 2*lmMst2)*s2t + s2t*pow2(Mst2)*(1540 + 2506*lmMst2 -
        96*(-2 + lmMst2)*pow2(lmMst1) - 141*pow2(lmMst2) + 6*lmMst1*(-460 - 6*
        lmMst2 + 123*pow2(lmMst2)) - 642*pow3(lmMst2)) + 64*Mst2*Mt*(5 + 149*
        lmMst2 + 12*pow2(lmMst2) + 6*lmMst1*(-7 - 8*lmMst2 + 6*pow2(lmMst2)) -
        36*pow3(lmMst2)))*pow4(Mst1) - 3*pow2(Mst1)*(-195*Dmsqst2*s2t + 128*
        Mst2*Mt*(-3 - 4*lmMst2 + 3*pow2(lmMst2)) + s2t*(-380 + 32*lmMst1*(-2 +
        lmMst2) - 38*lmMst2 + 251*pow2(lmMst2))*pow2(Mst2))*pow4(Mst2) - 2*(
        585*Dmsqst2*(lmMst1 - lmMst2)*s2t + 32*Mst2*Mt*(85 - 215*lmMst2 +
        lmMst1*(66 + 90*lmMst2 - 72*pow2(lmMst2)) - 54*pow2(lmMst2) + 72*pow3(
        lmMst2)) + 3*s2t*pow2(Mst2)*(-336 - 716*lmMst2 + 144*(-2 + lmMst2)*
        pow2(lmMst1) + lmMst1*(668 + 534*lmMst2 - 379*pow2(lmMst2)) - 246*pow2(
        lmMst2) + 235*pow3(lmMst2)))*pow6(Mst1) + 96*s2t*(-2 - lmMst2 + 3*pow2(
        lmMst2))*pow8(Mst2)))/pow2(Sbeta) - 300*pow2(Mt)*pow2(s2t)*((pow2(
        MuSUSY)*(-1540 - 2506*lmMst2 + 96*(-2 + lmMst2)*pow2(lmMst1) + lmMst1*(
        2760 + 36*lmMst2 - 738*pow2(lmMst2)) + 141*pow2(lmMst2) + 642*pow3(
        lmMst2)) + 4*pow2(Mst2)*(812 + 258*lmMst2 + 48*(-2 + lmMst2)*pow2(
        lmMst1) - 129*pow2(lmMst2) + 48*(lmMst1*(3 + 2*lmMst2 - 2*pow2(lmMst2))
        + pow3(lmMst2))))*pow4(Mst1)*pow4(Mst2) + 2*pow2(Mst2)*((1532 + 2478*
        lmMst2 - 96*lmMst1*(17 + 5*lmMst2) + 465*pow2(lmMst2))*pow2(Mst2) + 3*
        pow2(MuSUSY)*(-336 - 716*lmMst2 + 144*(-2 + lmMst2)*pow2(lmMst1) +
        lmMst1*(668 + 534*lmMst2 - 379*pow2(lmMst2)) - 246*pow2(lmMst2) + 235*
        pow3(lmMst2)))*pow6(Mst1) - 585*Dmsqst2*(pow4(Mst1)*((1 - 2*lmMst1 + 2*
        lmMst2)*pow2(Mst2)*pow2(MuSUSY) - 4*pow4(Mst2)) + pow2(Mst1)*(2*pow2(
        Mst2) + pow2(MuSUSY))*pow4(Mst2) + 2*(pow2(Mst2) + (-lmMst1 + lmMst2)*
        pow2(MuSUSY))*pow6(Mst1)) - 3*pow2(Mst1)*((568 - 64*lmMst1*(-2 +
        lmMst2) + 204*lmMst2 - 310*pow2(lmMst2))*pow2(Mst2) + (380 - 32*lmMst1*
        (-2 + lmMst2) + 38*lmMst2 - 251*pow2(lmMst2))*pow2(MuSUSY))*pow6(Mst2)
        + 96*(2 + lmMst2 - 3*pow2(lmMst2))*(2*pow2(Mst2) + pow2(MuSUSY))*pow8(
        Mst2))))))/(16200.*pow4(Mst1)*pow6(Mst2)) + xDmglst2*pow2(Dmglst2)*(-(
        twoLoopFlag*(72*Mt*pow3(s2t)*pow4(Mst2)*(8*pow2(Mst1)*pow2(Mst2) + (3 -
        6*lmMst1 + 6*lmMst2)*pow4(Mst1) + 2*(8 + lmMst1 - lmMst2)*pow4(Mst2)) +
        (96*Mst2*pow4(Mt)*(3*(37 + 41*lmMst2 - 13*lmMt - 6*lmMst2*lmMt +
        lmMst1*(-28 - 6*lmMst2 + 6*lmMt) + 6*pow2(lmMst2))*pow2(Mst2)*pow4(
        Mst1) + (17 - 3*lmMst1 + 9*lmMst2 - 3*lmMt)*pow2(Mst1)*pow4(Mst2) + 3*(
        75 + 174*lmMst2 - 37*lmMt - 30*lmMst2*lmMt + lmMst1*(-137 - 30*lmMst2 +
        30*lmMt) + 30*pow2(lmMst2))*pow6(Mst1) + 3*(-2 + lmMst2)*pow6(Mst2)))/
        pow2(Mst1) - 16*s2t*pow3(Mt)*(9*(2*(27 + 70*lmMst2 - 4*lmMst1*(14 + 3*
        lmMst2 - 3*lmMt) - 2*(7 + 6*lmMst2)*lmMt + 12*pow2(lmMst2))*pow2(Mst2)
        + 67*pow2(MuSUSY))*pow4(Mst1) + 18*(8 + lmMst1 - lmMst2)*pow2(MuSUSY)*
        pow4(Mst2) + 18*pow2(Mst1)*(2*(10 + lmMst1 - lmMst2)*pow2(Mst2)*pow2(
        MuSUSY) + (19 + 17*lmMst2 - 6*lmMt - 2*lmMst2*lmMt + lmMst1*(-11 - 2*
        lmMst2 + 2*lmMt) + 2*pow2(lmMst2))*pow4(Mst2)) + (19 + 36*lmMst1 - 42*
        lmMst2 + 6*lmMt)*pow6(Mst2)) - (12*Mst2*pow2(Mt)*pow2(s2t)*(-2*pow4(
        Mst1)*(3*(8 - 9*lmMst1 + 9*lmMst2 + 2*lmMst1*lmMst2 - 2*pow2(lmMst2))*
        pow2(Mst2)*pow2(MuSUSY) + (-13 + 18*lmMst1 - 24*lmMst2)*pow4(Mst2)) -
        3*(4*(-4 + 7*lmMst1 - 7*lmMst2)*pow2(Mst2) + (17 + 4*lmMst1*(-7 +
        lmMst2) + 28*lmMst2 - 4*pow2(lmMst2))*pow2(MuSUSY))*pow6(Mst1) + 6*(-2
        + lmMst2)*(2*pow2(Mst2) + pow2(MuSUSY))*pow6(Mst2) - 2*pow2(Mst1)*(3*(8
        - 5*lmMst1 + 4*lmMst2 + 2*lmMst1*lmMst2 - 2*pow2(lmMst2))*pow2(MuSUSY)*
        pow4(Mst2) + (-29 + 12*lmMst2)*pow6(Mst2))))/pow2(Mst1) + (-9*pow4(s2t)
        *pow5(Mst2)*(2*(-6 + lmMst1 - 2*lmMst1*lmMst2 + 2*pow2(lmMst2))*pow2(
        Mst2)*pow4(Mst1) + 2*(4 + 6*lmMst2 + lmMst1*(-5 + 2*lmMst2) - 2*pow2(
        lmMst2))*pow2(Mst1)*pow4(Mst2) + (1 - 2*lmMst1)*pow6(Mst1) - 2*(-2 +
        lmMst2)*pow6(Mst2)) + (36*s2t*pow2(Mt)*pow2(MuSUSY)*(-2*(-8*(10 +
        lmMst1 - lmMst2)*Mt + Mst2*s2t*(8 - 9*lmMst1 + 9*lmMst2 + 2*lmMst1*
        lmMst2 - 2*pow2(lmMst2)))*pow2(Mst2)*pow4(Mst1) - 2*(-4*(8 + lmMst1 -
        lmMst2)*Mt + Mst2*s2t*(8 - 5*lmMst1 + 4*lmMst2 + 2*lmMst1*lmMst2 - 2*
        pow2(lmMst2)))*pow2(Mst1)*pow4(Mst2) + (268*Mt + Mst2*s2t*(-17 - 4*
        lmMst1*(-7 + lmMst2) - 28*lmMst2 + 4*pow2(lmMst2)))*pow6(Mst1) + 2*(-2
        + lmMst2)*s2t*pow7(Mst2)))/pow2(Sbeta) - (4*Mt*MuSUSY*(-2*(pow2(Mst2)*(
        36*(5 - 3*lmMst1 + 3*lmMst2)*Mst2*s2t*pow2(Mt) - 54*(12 + lmMst1 -
        lmMst2)*Mt*pow2(Mst2)*pow2(s2t) + 4*(155 + 3*lmMst2*(41 - 6*lmMt) - 18*
        lmMst1*(4 + lmMst2 - lmMt) - 51*lmMt + 18*pow2(lmMst2))*pow3(Mt) - 9*(
        4*lmMst1 - 5*lmMst2)*pow3(Mst2)*pow3(s2t))*pow4(Mst1) + pow2(Mst1)*(6*(
        17 - 6*lmMst2)*Mst2*s2t*pow2(Mt) - 54*(8 + lmMst1 - lmMst2)*Mt*pow2(
        Mst2)*pow2(s2t) + 4*(11 + 18*lmMst1 - 21*lmMst2 + 3*lmMt)*pow3(Mt) + 9*
        (6 + 5*lmMst2 + lmMst1*(-5 + 2*lmMst2) - 2*pow2(lmMst2))*pow3(Mst2)*
        pow3(s2t))*pow4(Mst2)) - (72*(9 - 10*lmMst1 + 10*lmMst2)*Mst2*s2t*pow2(
        Mt) + 54*(-27 + 4*lmMst1 - 4*lmMst2)*Mt*pow2(Mst2)*pow2(s2t) + 8*(335 +
        699*lmMst2 - 18*lmMst1*(29 + 7*lmMst2 - 7*lmMt) - 3*(59 + 42*lmMst2)*
        lmMt + 126*pow2(lmMst2))*pow3(Mt) + 9*(1 - 10*lmMst1 + 10*lmMst2)*pow3(
        Mst2)*pow3(s2t))*pow6(Mst1) + 18*(-2 + lmMst2)*s2t*(-4*pow2(Mt) + pow2(
        Mst2)*pow2(s2t))*pow7(Mst2)))/Tbeta)/pow2(Mst1)))/(108.*pow7(Mst2)) + (
        Al4p*threeLoopFlag*((175*Mt*pow3(s2t)*(12*pow2(Mst1)*pow2(Mst2)*(282298
        - 61560*B4 + 1620*DN - 236520*lmMst1 - 180*(-439 + 438*lmMst1)*lmMst2 +
        2240*OepS2 - 81*(-1373 + 560*lmMst1 - 560*lmMst2)*S2 - 4320*pow2(
        lmMst1) + 1080*(77 + 48*lmMst1)*pow2(lmMst2) - 51840*pow3(lmMst2)) - (
        2727217 - 437920*OepS2 + 360*lmMst2*(4958 - 24633*S2) + 3648132*S2 +
        360*lmMst1*(-1460 + 1980*lmMst2 + 24633*S2) - 349920*pow2(lmMst1) -
        673920*pow2(lmMst2))*pow4(Mst1) + 103680*(3 + 4*lmMst2 - 3*pow2(lmMst2)
        )*pow4(Mst2)))/(Mst2*pow2(Mst1)) - (50*s2t*pow3(Mt)*((8*pow2(Mst2)*(
        7384247 - 12322800*lmMt + 1183840*OepS2 + 5821092*S2 + 22680*(41 + 4*
        lmMst2 - 4*lmMt)*pow2(lmMst1) + 1890*(-2489 + 48*lmMt)*pow2(lmMst2) +
        315*lmMst1*(1073 + 7734*lmMst2 + 690*lmMt - 76104*S2 + 864*pow2(lmMst2)
        - 864*pow2(lmMt)) + 816480*pow2(lmMt) + 315*lmMst2*(11587 + 966*lmMt +
        76104*S2 + 864*pow2(lmMt)) - 362880*pow3(lmMst2)) + 35*pow2(MuSUSY)*(
        996211 - 295488*B4 + 7776*DN - 1030176*lmMst1 - 144*(-1883 + 3618*
        lmMst1)*lmMst2 + 98336*OepS2 - 756*(259 + 2634*lmMst1 - 2634*lmMst2)*S2
        + 49248*pow2(lmMst1) + 5184*(67 + 48*lmMst1)*pow2(lmMst2) - 248832*
        pow3(lmMst2)))*pow4(Mst1) + 725760*(3 + 4*lmMst2 - 3*pow2(lmMst2))*(2*
        pow2(Mst2) + pow2(MuSUSY))*pow4(Mst2) + 12*pow2(Mst1)*(7*pow2(Mst2)*
        pow2(MuSUSY)*(334138 - 61560*B4 + 1620*DN - 236520*lmMst1 - 180*(-823 +
        438*lmMst1)*lmMst2 + 2240*OepS2 - 81*(-1373 + 560*lmMst1 - 560*lmMst2)*
        S2 - 4320*pow2(lmMst1) + 1080*(29 + 48*lmMst1)*pow2(lmMst2) - 51840*
        pow3(lmMst2)) - 2*(1243547 + 1157940*lmMt - 15680*OepS2 + 1680*lmMst2*(
        -1091 + 42*lmMt - 189*S2) - 953208*S2 + 7560*lmMst1*(135 + 33*lmMst2 -
        8*lmMt + 42*S2) + 30240*pow2(lmMst1) - 486360*pow2(lmMst2) + 15120*
        pow2(lmMt))*pow4(Mst2))))/(pow2(Mst1)*pow5(Mst2)) + 2551500*pow4(s2t)*(
        (-2*OepS2*(1136*pow2(Mst1) + 21*pow2(Mst2)))/729. + (S2*(4*(14023 +
        8520*lmMst1 - 8520*lmMst2)*pow2(Mst1) + 9*(65819 + 70*lmMst1 - 70*
        lmMst2)*pow2(Mst2)))/540. + pow2(Mst1)*(29.427737997256514 - B4/3. + (
        4*D3)/9. - (5*DN)/18. + (270961*lmMst1)/8100. + (653*pow2(lmMst1))/90.
         - (lmMst2*(332311 + 189090*lmMst1 + 57150*pow2(lmMst1)))/
        8100. + (7.95 + (74*lmMst1)/9.)*pow2(lmMst2) + (5*pow3(lmMst1))/18. - (
        13*pow3(lmMst2))/9.) - pow2(Mst2)*(47.56630658436214 - (19*B4)/9. + (
        20*D3)/9. - (7*DN)/6. + (163*lmMst1)/9. - (25*pow2(lmMst1))/36. + (2*
        lmMst2*(-347 - 18*lmMst1 + 12*pow2(lmMst1)))/27. - ((99 + 41*lmMst1)*
        pow2(lmMst2))/6. + (107*pow3(lmMst2))/18.) + ((940 - 32*lmMst1*(-2 +
        lmMst2) + 102*lmMst2 - 443*pow2(lmMst2))*pow4(Mst2))/(36.*pow2(Mst1)) -
        (8*(2 + lmMst2 - 3*pow2(lmMst2))*pow6(Mst2))/(9.*pow4(Mst1))) + (pow2(
        Mt)*pow2(s2t)*(-75*pow4(Mst1)*(-7*pow2(Mst2)*pow2(MuSUSY)*(13169 -
        41040*B4 + 43200*D3 - 22680*DN + 282960*lmMst1 + 1120*OepS2 - 324*(
        65819 + 70*lmMst1 - 70*lmMst2)*S2 - 13500*pow2(lmMst1) + 720*lmMst2*(-
        775 + 12*lmMst1 + 24*pow2(lmMst1)) - 1080*(-2 + 123*lmMst1)*pow2(
        lmMst2) + 115560*pow3(lmMst2)) + 18*(2164753 - 3360*B4 + 3360*D3 -
        1680*DN + 89040*lmMst1 - 6720*lmMt + 7840*OepS2 - 324*(523 + 490*lmMst1
        - 490*lmMst2)*S2 + 17220*pow2(lmMst1) - 140*lmMst2*(-6485 - 1098*lmMst1
        + 96*pow2(lmMst1)) - 420*(129 + 32*lmMst1)*pow2(lmMst2) + 26880*pow3(
        lmMst2))*pow4(Mst2)) - 2*(14*pow2(MuSUSY)*(2011073 + 1417500*B4 -
        1458000*D3 + 749250*DN + 934245*lmMst1 - 1178000*OepS2 + 1350*(620417 +
        17670*lmMst1 - 17670*lmMst2)*S2 + 3150900*pow2(lmMst1) - 45*lmMst2*(-
        124139 + 189090*lmMst1 + 71550*pow2(lmMst1)) + 4050*(1323 + 1970*
        lmMst1)*pow2(lmMst2) + 101250*pow3(lmMst1) - 4860000*pow3(lmMst2)) + 9*
        pow2(Mst2)*(68526791 + 9072000*lmMt + 2772000*OepS2 - 107403300*S2 +
        420*lmMst2*(102713 + 6000*lmMt + 133650*S2) + 6300*(131 - 30*lmMst2)*
        pow2(lmMst1) - 35116200*pow2(lmMst2) - 840*lmMst1*(-47431 - 41010*
        lmMst2 + 3000*lmMt + 66825*S2 + 6975*pow2(lmMst2)) + 63000*pow3(lmMst1)
        + 5985000*pow3(lmMst2)))*pow6(Mst1) + 9072000*(2 + lmMst2 - 3*pow2(
        lmMst2))*(2*pow2(Mst2) + pow2(MuSUSY))*pow6(Mst2) - 283500*pow2(Mst1)*(
        (812 - 32*lmMst1*(-2 + lmMst2) + 38*lmMst2 - 251*pow2(lmMst2))*pow2(
        MuSUSY)*pow4(Mst2) + (1432 - 64*lmMst1*(-2 + lmMst2) + 204*lmMst2 -
        310*pow2(lmMst2))*pow6(Mst2))))/(pow4(Mst1)*pow4(Mst2)) + (2551500*
        MuSUSY*(-(pow2(Mt)*pow2(s2t)*((64*Mst2*(3 + 4*lmMst2 - 3*pow2(lmMst2)))
        /(3.*pow2(Mst1)) + (761.0320987654321 - 152*B4 + 4*DN - 584*lmMst1 + (
        4*(631 - 438*lmMst1)*lmMst2)/9. - (32*pow2(lmMst1))/3. + (8*(53 + 48*
        lmMst1)*pow2(lmMst2))/3. - 128*pow3(lmMst2))/Mst2 + ((56*OepS2*(415*
        pow2(Mst1) + 24*pow2(Mst2)))/243. - (S2*((21422 + 87150*lmMst1 - 87150*
        lmMst2)*pow2(Mst1) + 9*(-1373 + 560*lmMst1 - 560*lmMst2)*pow2(Mst2)))/
        45.)/pow3(Mst2) + (pow2(Mst1)*(199.87633744855967 - 152*B4 + 4*DN - (
        12848*lmMst1)/27. - (8*(293 + 1152*lmMst1)*lmMst2)/27. + (184*pow2(
        lmMst1))/3. + 8*(35 + 16*lmMst1)*pow2(lmMst2) - 128*pow3(lmMst2)))/
        pow3(Mst2))) + Mt*(pow3(s2t)*(92.93189300411522 - (76*B4)/9. + (80*D3)/
        9. - (14*DN)/3. + (196*lmMst1)/3. - (25*pow2(lmMst1))/9. + (2*lmMst2*(-
        1493 - 24*lmMst1 + 48*pow2(lmMst1)))/27. - ((247 + 246*lmMst1)*pow2(
        lmMst2))/9. + (8*OepS2*(21 + (1157*pow2(Mst1))/pow2(Mst2)))/729. + ((-
        876 + 32*lmMst1*(-2 + lmMst2) - 70*lmMst2 + 347*pow2(lmMst2))*pow2(
        Mst2))/(9.*pow2(Mst1)) - (S2*((648463 + 34710*lmMst1 - 34710*lmMst2)*
        pow2(Mst1) + 9*(65819 + 70*lmMst1 - 70*lmMst2)*pow2(Mst2)))/(135.*pow2(
        Mst2)) - (pow2(Mst1)*(24.779058984910836 + (64*B4)/9. - (64*D3)/9. + (
        32*DN)/9. + (138661*lmMst1)/2025. + (159*pow2(lmMst1))/5. - lmMst2*(
        53.5116049382716 + (458*lmMst1)/5. + (286*pow2(lmMst1))/9.) + (2*(1333
        + 1355*lmMst1)*pow2(lmMst2))/45. + (10*pow3(lmMst1))/9. - (266*pow3(
        lmMst2))/9.))/pow2(Mst2) + (214*pow3(lmMst2))/9. + (32*(2 + lmMst2 - 3*
        pow2(lmMst2))*pow4(Mst2))/(9.*pow4(Mst1))) - (4*T1ep*(21*pow2(Mst2)*(
        36*Mst2*s2t*pow2(Mt) - 24*Mt*pow2(Mst2)*pow2(s2t) + 32*pow3(Mt) + pow3(
        Mst2)*pow3(s2t)) + pow2(Mst1)*(4320*Mst2*s2t*pow2(Mt) - 8715*Mt*pow2(
        Mst2)*pow2(s2t) + 17584*pow3(Mt) + 1157*pow3(Mst2)*pow3(s2t))))/(243.*
        pow5(Mst2))) + (4*(-3*(1648637 + 1173060*lmMt - 15680*OepS2 + 1680*
        lmMst2*(-1214 + 42*lmMt - 189*S2) - 953208*S2 + 7560*lmMst1*(131 + 33*
        lmMst2 - 8*lmMt + 42*S2) + 30240*pow2(lmMst1) - 304920*pow2(lmMst2) +
        15120*pow2(lmMt))*pow2(Mst1)*pow2(Mst2) + 2*(-1657412 - 7512750*lmMt +
        615440*OepS2 + 4340358*S2 + 11340*(35 + 4*lmMst2 - 4*lmMt)*pow2(lmMst1)
        + 945*(-1981 + 48*lmMt)*pow2(lmMst2) + 315*lmMst1*(-4223 + 2679*lmMst2
        + 201*lmMt - 39564*S2 + 432*pow2(lmMst2) - 432*pow2(lmMt)) + 385560*
        pow2(lmMt) + 945*lmMst2*(5321 + 193*lmMt + 13188*S2 + 144*pow2(lmMt)) -
        181440*pow3(lmMst2))*pow4(Mst1) - 181440*(-3 - 4*lmMst2 + 3*pow2(
        lmMst2))*pow4(Mst2))*pow4(Mt))/(25515.*pow2(Mst1)*pow5(Mst2)) + (s2t*
        pow3(Mt)*(75*pow2(Mst2)*(2435233 - 3360*B4 + 3360*D3 - 1680*DN +
        115920*lmMst1 - 6720*lmMt + 7840*OepS2 - 324*(523 + 490*lmMst1 - 490*
        lmMst2)*S2 + 17220*pow2(lmMst1) - 140*lmMst2*(-6695 - 1002*lmMst1 + 96*
        pow2(lmMst1)) - 1680*(47 + 8*lmMst1)*pow2(lmMst2) + 26880*pow3(lmMst2))
        *pow4(Mst1) + 31500*(652 - 32*lmMst1*(-2 + lmMst2) + 70*lmMst2 - 59*
        pow2(lmMst2))*pow2(Mst1)*pow4(Mst2) + 2*(127852633 - 126000*B4 +
        126000*D3 - 63000*DN + 23638020*lmMst1 - 252000*(-17 + 5*lmMst1 - 5*
        lmMst2)*lmMt + 1680000*OepS2 - 2700*(22243 + 12600*lmMst1 - 12600*
        lmMst2)*S2 + 1058400*pow2(lmMst1) - 420*lmMst2*(-136544 - 53535*lmMst1
        + 1425*pow2(lmMst1)) - 6300*(3257 + 545*lmMst1)*pow2(lmMst2) + 31500*
        pow3(lmMst1) + 4000500*pow3(lmMst2))*pow6(Mst1) - 1008000*(2 + lmMst2 -
        3*pow2(lmMst2))*pow6(Mst2)))/(70875.*pow4(Mst1)*pow4(Mst2))))/Tbeta + (
        4*pow4(Mt)*(15*pow4(Mst1)*(840*pow2(Mst2)*pow2(MuSUSY)*(33934 - 90*B4 +
        90*D3 - 45*DN + 720*lmMst1 + 120*(163 + 24*lmMst1)*lmMst2 - 120*lmMt -
        10206*S2 - 720*(2 + lmMst1)*pow2(lmMst2) + 720*pow3(lmMst2)) - (
        46632377 + 17278380*lmMt - 744800*OepS2 - 29679210*S2 - 168*lmMst2*(-
        56837 + 7170*lmMt + 89775*S2) + 132300*pow2(lmMst1) + 10427760*pow2(
        lmMst2) + 37800*lmMst1*(66 - 240*lmMst2 + 399*S2 + 32*(lmMt + pow2(
        lmMst2))) - 1587600*pow2(lmMt) - 1209600*pow3(lmMst2))*pow4(Mst2)) + (
        pow2(Mst2)*(648916519 - 1795878000*lmMt + 158732000*OepS2 + 66753450*S2
        + 113400*(1403 + 285*lmMst2 - 225*lmMt)*pow2(lmMst1) + 94500*(-15839 +
        558*lmMt)*pow2(lmMst2) + 3780*lmMst1*(132381 + lmMst2*(272035 - 7200*
        lmMt) - 21415*lmMt - 850350*S2 + 22200*pow2(lmMst2) - 21600*pow2(lmMt))
        + 176904000*pow2(lmMt) + 1260*lmMst2*(-245893 + 169395*lmMt + 2551050*
        S2 + 64800*pow2(lmMt)) - 2268000*pow3(lmMst1) - 113967000*pow3(lmMst2))
        + 525*pow2(MuSUSY)*(2199511 - 4320*B4 + 4320*D3 - 2160*DN + 140160*
        lmMst1 + 960*(1426 + 303*lmMst1)*lmMst2 - 2880*(-16 + 5*lmMst1 - 5*
        lmMst2)*lmMt - 1120*OepS2 + 324*(-3411 + 70*lmMst1 - 70*lmMst2)*S2 -
        2880*(77 + 24*lmMst1)*pow2(lmMst2) + 69120*pow3(lmMst2)))*pow6(Mst1) +
        283500*(620 - 32*lmMst1*(-2 + lmMst2) + 166*lmMst2 - 59*pow2(lmMst2))*
        pow2(Mst1)*pow6(Mst2) - 9072000*(2 + lmMst2 - 3*pow2(lmMst2))*pow8(
        Mst2)))/(pow4(Mst1)*pow6(Mst2)) + (-10500*T1ep*(4*pow2(Mst1)*(-2*pow2(
        Mst2)*pow2(Mt)*pow2(s2t)*(-589*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 891*
        pow2(Mst2)*pow2(Sbeta)) - 14*Mst2*s2t*(439*pow2(MuSUSY)*(-1 + pow2(
        Sbeta)) + 1208*pow2(Mst2)*pow2(Sbeta))*pow3(Mt) + 4*(-21*pow2(MuSUSY)*(
        -1 + pow2(Sbeta)) + 5669*pow2(Mst2)*pow2(Sbeta))*pow4(Mt) + (2737*Mt -
        284*Mst2*s2t)*pow2(Sbeta)*pow3(s2t)*pow5(Mst2)) + 21*pow3(Mst2)*(-4*
        Mst2*pow2(Mt)*pow2(s2t)*(-(pow2(MuSUSY)*(-1 + pow2(Sbeta))) + 18*pow2(
        Mst2)*pow2(Sbeta)) - 64*s2t*(pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 2*pow2(
        Mst2)*pow2(Sbeta))*pow3(Mt) + pow2(Sbeta)*(32*Mt*pow3(s2t)*pow4(Mst2) +
        304*Mst2*pow4(Mt) - pow4(s2t)*pow5(Mst2)))) - (7*pow2(Mt)*pow2(MuSUSY)*
        (75*pow2(Mst2)*(-8*Mst2*Mt*s2t*(334138 - 61560*B4 + 1620*DN - 236520*
        lmMst1 - 180*(-823 + 438*lmMst1)*lmMst2 + 2240*OepS2 - 81*(-1373 + 560*
        lmMst1 - 560*lmMst2)*S2 - 4320*pow2(lmMst1) + 1080*(29 + 48*lmMst1)*
        pow2(lmMst2) - 51840*pow3(lmMst2)) + 96*pow2(Mt)*(33934 - 90*B4 + 90*D3
        - 45*DN + 720*lmMst1 + 120*(163 + 24*lmMst1)*lmMst2 - 120*lmMt - 10206*
        S2 - 720*(2 + lmMst1)*pow2(lmMst2) + 720*pow3(lmMst2)) + pow2(Mst2)*
        pow2(s2t)*(13169 - 41040*B4 + 43200*D3 - 22680*DN + 282960*lmMst1 +
        1120*OepS2 - 324*(65819 + 70*lmMst1 - 70*lmMst2)*S2 - 13500*pow2(
        lmMst1) + 720*lmMst2*(-775 + 12*lmMst1 + 24*pow2(lmMst1)) - 1080*(-2 +
        123*lmMst1)*pow2(lmMst2) + 115560*pow3(lmMst2)))*pow4(Mst1) - 40500*
        s2t*(Mst2*s2t*(812 - 32*lmMst1*(-2 + lmMst2) + 38*lmMst2 - 251*pow2(
        lmMst2)) + 128*Mt*(3 + 4*lmMst2 - 3*pow2(lmMst2)))*pow2(Mst1)*pow5(
        Mst2) + 2*(-2*pow2(Mst2)*pow2(s2t)*(2011073 + 1417500*B4 - 1458000*D3 +
        749250*DN + 934245*lmMst1 - 1178000*OepS2 + 1350*(620417 + 17670*lmMst1
        - 17670*lmMst2)*S2 + 3150900*pow2(lmMst1) - 45*lmMst2*(-124139 +
        189090*lmMst1 + 71550*pow2(lmMst1)) + 4050*(1323 + 1970*lmMst1)*pow2(
        lmMst2) + 101250*pow3(lmMst1) - 4860000*pow3(lmMst2)) - 125*Mst2*Mt*
        s2t*(996211 - 295488*B4 + 7776*DN - 1030176*lmMst1 - 144*(-1883 + 3618*
        lmMst1)*lmMst2 + 98336*OepS2 - 756*(259 + 2634*lmMst1 - 2634*lmMst2)*S2
        + 49248*pow2(lmMst1) + 5184*(67 + 48*lmMst1)*pow2(lmMst2) - 248832*
        pow3(lmMst2)) + 150*pow2(Mt)*(2199511 - 4320*B4 + 4320*D3 - 2160*DN +
        140160*lmMst1 + 960*(1426 + 303*lmMst1)*lmMst2 - 2880*(-16 + 5*lmMst1 -
        5*lmMst2)*lmMt - 1120*OepS2 + 324*(-3411 + 70*lmMst1 - 70*lmMst2)*S2 -
        2880*(77 + 24*lmMst1)*pow2(lmMst2) + 69120*pow3(lmMst2)))*pow6(Mst1) +
        1296000*(2 + lmMst2 - 3*pow2(lmMst2))*pow2(s2t)*pow8(Mst2)))/pow4(Mst1)
        )/(pow2(Sbeta)*pow6(Mst2))))/2.5515e6)) - ((xMst*pow6(Mst1)*((-27*(
        lmMst1 - lmMst2)*oneLoopFlag*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*(-1 +
        pow2(Sbeta))*pow3(Mst2))/pow2(Sbeta) - Al4p*twoLoopFlag*(8*Mst2*s2t*
        pow2(Mt)*(Mst2*(Mt*(-193 - 474*lmMst2 + 6*lmMst1*(67 + 42*lmMst2) -
        252*pow2(lmMst2)) + Dmglst2*s2t*(49 - 84*lmMst2 + 12*lmMst1*(7 + 3*
        lmMst2) - 36*pow2(lmMst2)))*pow2(MuSUSY) + Dmglst2*Mt*(-785 + 438*
        lmMst2 + 6*lmMst1*(-85 + 42*lmMst2) - 252*pow2(lmMst2))*pow2(MuSUSY) +
        s2t*(1 - 111*lmMst2 + 3*lmMst1*(37 + 72*lmMst2) - 81*pow2(lmMst1) -
        135*pow2(lmMst2))*pow2(Mst2)*pow2(MuSUSY) + 4*(9*Dmglst2*(-lmMst1 +
        lmMst2)*s2t + Mt*(49 - 75*lmMst2 + 3*lmMst1*(25 + 6*lmMst2 - 6*lmMt) +
        18*lmMst2*lmMt - 18*pow2(lmMst2)))*pow3(Mst2)) - 32*(Dmglst2*s2t*(95 -
        447*lmMst2 + lmMst1*(411 + 90*lmMst2 - 90*lmMt) + 36*lmMt + 90*lmMst2*
        lmMt - 90*pow2(lmMst2)) + Mt*(65 - 159*lmMst2 + 6*lmMst1*(25 + 6*lmMst2
        - 6*lmMt) + 9*lmMt + 36*lmMst2*lmMt - 36*pow2(lmMst2)))*pow3(Mst2)*
        pow3(Mt) + 2*Mst2*Mt*MuSUSY*((4*Mt*MuSUSY*s2t*(Dmglst2*Mst2*s2t*(-49 +
        84*lmMst2 - 12*lmMst1*(7 + 3*lmMst2) + 36*pow2(lmMst2)) + Dmglst2*Mt*(
        785 - 438*lmMst2 - 6*lmMst1*(-85 + 42*lmMst2) + 252*pow2(lmMst2)) +
        Mst2*Mt*(193 + 474*lmMst2 - 6*lmMst1*(67 + 42*lmMst2) + 252*pow2(
        lmMst2)) + s2t*(-1 + 111*lmMst2 - 3*lmMst1*(37 + 72*lmMst2) + 81*pow2(
        lmMst1) + 135*pow2(lmMst2))*pow2(Mst2)))/pow2(Sbeta) + (8*Dmglst2*(36*(
        -1 + 3*lmMst1 - 3*lmMst2)*Mst2*s2t*pow2(Mt) + 3*Mt*(-50 + 51*lmMst2 +
        3*lmMst1*(-17 + 6*lmMst2) - 18*pow2(lmMst2))*pow2(Mst2)*pow2(s2t) + 4*(
        7 - 381*lmMst2 + lmMst1*(327 + 72*lmMst2 - 81*lmMt) + 54*lmMt + 81*
        lmMst2*lmMt - 72*pow2(lmMst2))*pow3(Mt) + 2*(1 + 3*lmMst1 - 3*lmMst2)*
        pow3(Mst2)*pow3(s2t)) - Mst2*(-144*Mst2*s2t*(-lmMst1 + lmMst2 - 2*
        lmMst1*lmMst2 + pow2(lmMst1) + pow2(lmMst2))*pow2(Mt) + 24*Mt*(1 + 33*
        lmMst2 - 3*lmMst1*(11 + 6*lmMst2) + 18*pow2(lmMst2))*pow2(Mst2)*pow2(
        s2t) + 32*(83 - 96*lmMst2 + 3*lmMst1*(29 + 12*lmMst2 - 9*lmMt) + 9*lmMt
        + 27*lmMst2*lmMt - 36*pow2(lmMst2))*pow3(Mt) + (5 + lmMst1*(6 - 288*
        lmMst2) - 6*lmMst2 + 144*pow2(lmMst1) + 144*pow2(lmMst2))*pow3(Mst2)*
        pow3(s2t)))/Tbeta) + 576*Dmglst2*(4 + 6*lmMst1*(9 + 2*lmMst2 - 2*lmMt)
        + 7*lmMt + lmMst2*(-61 + 12*lmMt) - 12*pow2(lmMst2))*pow2(Mst2)*pow4(
        Mt) + 4*Dmglst2*(11 + 42*lmMst1 - 42*lmMst2)*Mt*pow3(s2t)*pow5(Mst2) +
        (-144*(lmMst1 - lmMst2)*Mst2*xDR2DRMOD*(2*Dmglst2*lmMst2*Mst2 + (-2 +
        lmMst2)*xDmglst2*pow2(Dmglst2) + (1 + lmMst2)*pow2(Mst2))*pow2(Mt)*
        pow2(MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta)) + (xDmglst2*pow2(Dmglst2)*(
        24*s2t*(4*(-143 + 18*lmMst1 - 18*lmMst2)*Tbeta*pow2(MuSUSY)*(-1 + pow2(
        Sbeta)) + Mst2*pow2(Sbeta)*(-12*Mst2*Tbeta*(18 + lmMst2*(169 - 30*lmMt)
        - 3*lmMst1*(49 + 10*lmMst2 - 10*lmMt) - 22*lmMt + 30*pow2(lmMst2)) +
        MuSUSY*(156 + 252*lmMst2 - 5*(43 + 60*lmMst2)*pow2(Sbeta) + 12*lmMst1*(
        -21 + 25*pow2(Sbeta)))))*pow3(Mt) - 2*Mst2*pow2(Mt)*pow2(s2t)*(-2*
        Tbeta*(157 + 348*lmMst2 + 12*lmMst1*(-29 + 3*lmMst2) - 36*pow2(lmMst2))
        *pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 36*(85 - 36*lmMst1 + 36*lmMst2)*
        Mst2*MuSUSY*pow2(Sbeta) + 3*Tbeta*pow2(Mst2)*pow2(Sbeta)*(96 - 430*
        pow2(Sbeta) + 215*pow4(Sbeta) - 12*(lmMst1 - lmMst2)*(22 - 50*pow2(
        Sbeta) + 25*pow4(Sbeta)))) + pow2(Sbeta)*(16*Mt*((1 - 24*lmMst1 + 24*
        lmMst2)*MuSUSY + 6*(1 - 6*lmMst1 + 6*lmMst2)*Mst2*Tbeta)*pow3(Mst2)*
        pow3(s2t) + 32*(9*Mst2*Tbeta*(94 + 475*lmMst2 - 6*lmMst1*(67 + 14*
        lmMst2 - 14*lmMt) - 73*lmMt - 84*lmMst2*lmMt + 84*pow2(lmMst2)) +
        MuSUSY*(398 + 2085*lmMst2 - 18*lmMst1*(95 + 22*lmMst2 - 22*lmMt) - 375*
        lmMt - 396*lmMst2*lmMt + 396*pow2(lmMst2)))*pow4(Mt) + (5 + 6*lmMst1 -
        6*lmMst2)*Tbeta*pow4(s2t)*pow5(Mst2))))/Tbeta)/pow2(Sbeta) + 2*(5 + 6*
        lmMst1 - 6*lmMst2)*(-2*Mt + Dmglst2*s2t)*pow3(s2t)*pow6(Mst2) - (11 +
        6*lmMst1 - 6*lmMst2)*pow4(s2t)*pow7(Mst2))))/108. + (Al4p*z2*(-6*
        xDmglst2*pow2(Dmglst2)*pow2(Mst2)*(22680*twoLoopFlag*pow2(Mst1)*(-4*
        pow2(Mt)*(-8*(7*MuSUSY + 15*Mst2*Tbeta)*pow2(Mt)*pow2(Sbeta) + Mst2*
        MuSUSY*pow2(s2t)*(MuSUSY*Tbeta*(-1 + pow2(Sbeta)) - 18*Mst2*pow2(Sbeta)
        ) + 12*Mt*s2t*Tbeta*(-3*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 4*pow2(Mst2)*
        pow2(Sbeta)))*pow4(Mst1) + s2t*(4*(12*Mt - Mst2*s2t)*Tbeta*pow2(Mt)*
        pow2(MuSUSY)*(-1 + pow2(Sbeta)) + s2t*pow2(Mst2)*pow2(Sbeta)*(-4*Mst2*
        Mt*s2t*(MuSUSY + 6*Mst2*Tbeta) + 72*MuSUSY*pow2(Mt) + Tbeta*pow2(s2t)*
        pow3(Mst2)))*pow4(Mst2) - pow2(Mst1)*pow2(Mst2)*(-4*Mst2*MuSUSY*pow2(
        Mt)*pow2(s2t)*(-(MuSUSY*Tbeta*(-1 + pow2(Sbeta))) + 18*Mst2*pow2(Sbeta)
        ) + 32*s2t*Tbeta*(-3*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + pow2(Mst2)*pow2(
        Sbeta))*pow3(Mt) - 32*(MuSUSY + 3*Mst2*Tbeta)*pow2(Sbeta)*pow4(Mt) +
        Tbeta*pow2(Sbeta)*pow4(s2t)*pow5(Mst2))) + Al4p*Mst2*threeLoopFlag*(4*
        pow4(Mst1)*(2*pow2(Mst2)*pow2(Mt)*(-7*(-261247 + 163890*lmMst1 -
        140670*lmMst2)*Tbeta*pow2(MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta)) + 2*Mt*(
        36*(297316 - 20265*lmMst1 + 231945*lmMst2)*MuSUSY*s2t + (71735177 -
        8595720*lmMst1 + 18272520*lmMst2 - 6940080*lmMt)*Mt*Tbeta)*pow2(Sbeta))
        - s2t*(-21*(-1486429 + 190080*lmMst1 + 43200*lmMst2)*MuSUSY*s2t + 32*(
        4999046 - 615195*lmMst1 + 1302210*lmMst2 - 566055*lmMt)*Mt*Tbeta)*pow2(
        Mt)*pow2(Sbeta)*pow3(Mst2) + 2*Mst2*MuSUSY*(35*(-184517 + 74952*lmMst1
        + 18360*lmMst2)*MuSUSY*s2t*Tbeta*(-1 + pow2(Sbeta)) + 8*(10736701 -
        1553580*lmMst1 + 3303720*lmMst2 - 1508220*lmMt)*Mt*pow2(Sbeta))*pow3(
        Mt) + Mt*(7*(642497 - 170100*lmMst1 + 170100*lmMst2)*MuSUSY*s2t + 90*(-
        186703 + 16632*lmMst1 - 97272*lmMst2)*Mt*Tbeta)*pow2(s2t)*pow2(Sbeta)*
        pow4(Mst2) + 84*(12761 + 5400*lmMst1 + 178920*lmMst2)*Tbeta*pow2(
        MuSUSY)*(-1 + pow2(Sbeta))*pow4(Mt) + 7*((2050273 - 5400*lmMst1 + 5400*
        lmMst2)*Mt + 5*(-36721 + 621*lmMst1 - 2943*lmMst2)*Mst2*s2t)*Tbeta*
        pow2(Sbeta)*pow3(s2t)*pow5(Mst2)) - 3*pow2(Mst1)*pow2(Mst2)*(-4*pow2(
        Mst2)*pow2(Mt)*(-7*(40001 + 52560*lmMst1 - 37080*lmMst2)*Tbeta*pow2(
        MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta)) + 4*Mt*(3*(255749 + 2100*lmMst1 +
        441420*lmMst2)*MuSUSY*s2t + (3572731 - 486360*lmMst1 + 912240*lmMst2 -
        587160*lmMt)*Mt*Tbeta)*pow2(Sbeta)) + 4*Mt*s2t*(Mst2*s2t*(7*(30641 +
        52560*lmMst1 - 37080*lmMst2)*MuSUSY*s2t + 6*(233909 + 2100*lmMst1 +
        441420*lmMst2)*Mt*Tbeta) - 4*Mt*(21*(46987 + 15390*lmMst1 + 4050*
        lmMst2)*MuSUSY*s2t + 4*(-497573 + 107730*lmMst1 - 233100*lmMst2 +
        125370*lmMt)*Mt*Tbeta))*pow2(Sbeta)*pow3(Mst2) - 32*Mst2*MuSUSY*(7*(
        46987 + 15390*lmMst1 + 4050*lmMst2)*MuSUSY*s2t*Tbeta*(-1 + pow2(Sbeta))
        + 4*(217444 - 53865*lmMst1 + 116550*lmMst2 - 62685*lmMt)*Mt*pow2(Sbeta)
        )*pow3(Mt) + Tbeta*(-5376*(17 + 1320*lmMst2)*pow2(MuSUSY)*(-1 + pow2(
        Sbeta))*pow4(Mt) + 7*(16*(46987 + 15390*lmMst1 + 4050*lmMst2)*Mt + (-
        21281 - 52560*lmMst1 + 37080*lmMst2)*Mst2*s2t)*pow2(Sbeta)*pow3(s2t)*
        pow5(Mst2))) + 196560*(-4*Tbeta*pow2(Mt)*pow2(s2t)*(pow2(MuSUSY)*(-1 +
        pow2(Sbeta)) + 2*pow2(Mst2)*pow2(Sbeta)) + pow2(Sbeta)*(16*MuSUSY*s2t*
        pow3(Mt) - 4*Mt*MuSUSY*pow2(Mst2)*pow3(s2t) + 16*Tbeta*pow4(Mt) +
        Tbeta*pow4(Mst2)*pow4(s2t)))*pow6(Mst2))) + 544320*twoLoopFlag*xMst*
        pow2(Mt)*(xDmglst2*pow2(Dmglst2)*(-16*(11*MuSUSY + 21*Mst2*Tbeta)*pow2(
        Mt)*pow2(Sbeta) + Mst2*MuSUSY*pow2(s2t)*(MuSUSY*Tbeta*(-1 + pow2(Sbeta)
        ) - 18*Mst2*pow2(Sbeta)) - 24*Mt*s2t*Tbeta*(2*pow2(MuSUSY)*(-1 + pow2(
        Sbeta)) - 5*pow2(Mst2)*pow2(Sbeta))) + 2*Mt*pow2(Mst2)*(-8*Mt*(2*MuSUSY
        + Mst2*Tbeta)*pow2(Sbeta) + s2t*Tbeta*(-(pow2(MuSUSY)*(-1 + pow2(Sbeta)
        )) + 4*pow2(Mst2)*pow2(Sbeta))) - 2*Dmglst2*Mst2*(-16*(2*MuSUSY + 3*
        Mst2*Tbeta)*pow2(Mt)*pow2(Sbeta) + Mst2*MuSUSY*pow2(s2t)*(-(MuSUSY*
        Tbeta*(-1 + pow2(Sbeta))) + 6*Mst2*pow2(Sbeta)) + Mt*s2t*Tbeta*(17*
        pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 20*pow2(Mst2)*pow2(Sbeta))))*pow8(
        Mst1) - 420*Al4p*threeLoopFlag*pow3(Mst2)*(xDmsqst2*pow2(Dmsqst2)*(4*
        Mt*MuSUSY*pow2(Sbeta)*(243*Mst2*Mt*pow2(Mst1)*(80*(1 + lmMst1 - lmMst2)
        *pow2(Mst1) + 40*(3 + lmMst1 - lmMst2)*pow2(Mst2))*pow2(s2t) - 12960*(-
        2 + lmMst1 - 2*lmMst2 + lmMt)*Mst2*pow2(Mst1)*pow3(Mt) + 30*s2t*pow2(
        Mt)*(627*pow2(Mst1)*pow2(Mst2) + 1348*pow4(Mst1) + 108*pow4(Mst2)) - 5*
        pow3(s2t)*(3*(1057 + 162*lmMst1 - 162*lmMst2)*pow2(Mst2)*pow4(Mst1) +
        9*(191 + 18*lmMst1 - 18*lmMst2)*pow2(Mst1)*pow4(Mst2) + 4180*pow6(Mst1)
        + 162*pow6(Mst2))) + Tbeta*(-25920*(3 + lmMst1 - lmMst2)*Mst2*s2t*pow2(
        Mst1)*pow2(MuSUSY)*pow3(Mt) + 162*Mt*pow2(s2t)*((10*(209 + 18*lmMst1 -
        18*lmMst2)*Mt*pow2(Mst1)*pow2(Mst2)*pow2(MuSUSY))/9. + (40*(421 + 54*
        lmMst1 - 54*lmMst2)*Mt*pow2(MuSUSY)*pow4(Mst1))/27. + 20*Mt*pow2(
        MuSUSY)*pow4(Mst2) + 5*((-51 + 6*lmMst1 - 6*lmMst2)*Mt - 2*(11 + 10*
        lmMst1 - 10*lmMst2)*Mst2*s2t)*(-2 + pow2(Sbeta))*pow4(Sbeta)*pow6(Mst1)
        ) + pow2(Sbeta)*(-25920*Mst2*s2t*pow2(Mst1)*((14 - 4*lmMst1 + 8*lmMst2
        - 4*lmMt)*pow2(Mst1) + (5 - 2*lmMst1 + 4*lmMst2 - 2*lmMt)*pow2(Mst2) +
        (-3 - lmMst1 + lmMst2)*pow2(MuSUSY))*pow3(Mt) + 3240*Mst2*Mt*pow2(Mst1)
        *pow3(s2t)*(4*(1 - lmMst1 + lmMst2)*pow2(Mst1)*pow2(Mst2) + (17 + 10*
        lmMst1 - 10*lmMst2)*pow4(Mst1) - 4*(3 + lmMst1 - lmMst2)*pow4(Mst2)) +
        12960*((13 - 3*lmMst1 + 6*lmMst2 - 3*lmMt)*pow2(Mst1)*pow2(Mst2) + 5*(9
        - 2*lmMst1 + 4*lmMst2 - 2*lmMt)*pow4(Mst1) + pow4(Mst2))*pow4(Mt) + 9*
        pow2(Mst2)*pow4(s2t)*((20*(121 + 27*lmMst1 - 27*lmMst2)*pow2(Mst2)*
        pow4(Mst1))/3. + 5*(173 + 18*lmMst1 - 18*lmMst2)*pow2(Mst1)*pow4(Mst2)
        + (560.5555555555555 - 270*lmMst1 + 270*lmMst2)*pow6(Mst1) + 90*pow6(
        Mst2)) - 60*pow2(Mt)*pow2(s2t)*(721*(pow2(Mst1) + pow2(Mst2))*pow4(
        Mst1) + 519*pow2(Mst1)*pow4(Mst2) + pow2(MuSUSY)*(3*(209 + 18*lmMst1 -
        18*lmMst2)*pow2(Mst1)*pow2(Mst2) + 4*(421 + 54*lmMst1 - 54*lmMst2)*
        pow4(Mst1) + 54*pow4(Mst2)) + 108*pow6(Mst2))))) - 108*xDR2DRMOD*(8*Mt*
        MuSUSY*pow2(Sbeta)*(pow2(Mst1)*pow2(Mst2)*pow2(s2t)*(-(s2t*pow3(Mst2)*(
        (7*Dmglst2 - 13*Mst2)*pow2(Mst1) + 2*Dmglst2*(7*lmMst1 - 15*lmMst2)*
        pow2(Mst2) - 2*(4 + 13*lmMst1 - 9*lmMst2)*pow3(Mst2))) - 48*(Dmglst2 +
        3*Dmglst2*lmMst2 + Mst2 + lmMst2*Mst2)*Mt*(pow2(Mst1)*pow2(Mst2) +
        pow4(Mst1) + pow4(Mst2))) + 64*pow2(Mst1)*pow3(Mt)*(-(Dmglst2*(7 +
        lmMst2)*pow2(Mst1)*pow2(Mst2)) - 2*Dmglst2*(11 + 5*lmMst2)*pow4(Mst1) +
        Dmglst2*(-1 + lmMst2)*pow4(Mst2) + (1 + lmMst2)*Mst2*(3*pow2(Mst1)*
        pow2(Mst2) + 6*pow4(Mst1) + pow4(Mst2))) + 4*(7*Dmglst2 - 13*Mst2)*s2t*
        (pow2(Mst1) - pow2(Mst2))*pow2(Mt)*pow5(Mst2)) + 8*Mt*pow2(Sbeta)*(-80*
        Tbeta*xDmglst2*pow2(Dmglst2)*pow3(Mt) + (7*Dmglst2 - 13*Mst2)*MuSUSY*
        pow3(Mst2)*pow3(s2t))*pow6(Mst2) + Mst2*xDmglst2*pow2(Dmglst2)*(32*Mt*
        MuSUSY*pow2(Sbeta)*(12*(2 - 3*lmMst2)*Mt*pow2(Mst1)*pow2(Mst2)*(pow2(
        Mst1) + pow2(Mst2))*pow2(s2t) + 20*s2t*(pow2(Mst1) - pow2(Mst2))*pow2(
        Mt)*pow3(Mst2) + 16*pow2(Mst1)*((7 + 2*lmMst2)*pow2(Mst1) - pow2(Mst2))
        *pow3(Mt) - pow3(Mst2)*pow3(s2t)*(2*(2 + 5*lmMst1 - 6*lmMst2)*pow2(
        Mst1)*pow2(Mst2) + 5*pow4(Mst1) - 5*pow4(Mst2))) + Tbeta*(-32*s2t*pow2(
        Mt)*pow2(MuSUSY)*((-8*(-2 + 3*lmMst2)*Mt + (1 - 10*lmMst1 + 12*lmMst2)*
        Mst2*s2t)*pow2(Mst1)*pow2(Mst2) - 2*(8*(-2 + 3*lmMst2)*Mt + (2 + 5*
        lmMst1 - 6*lmMst2)*Mst2*s2t)*pow4(Mst1) + 5*s2t*pow5(Mst2)) + pow2(
        Sbeta)*(128*pow2(Mst1)*(pow2(Mst2)*(-(Mst2*Mt) + 4*s2t*pow2(Mst2) + 2*(
        2 - 3*lmMst2)*s2t*pow2(MuSUSY)) - 4*pow2(Mst1)*(-3*(8 + 3*lmMst2)*Mst2*
        Mt + 2*(4 + lmMst2)*s2t*pow2(Mst2) + (-2 + 3*lmMst2)*s2t*pow2(MuSUSY)))
        *pow3(Mt) - 32*Mst2*pow2(Mt)*pow2(s2t)*((2*(2 + 5*lmMst1 - 6*lmMst2)*
        pow2(Mst1) - 5*pow2(Mst2))*(pow2(Mst1) + pow2(Mst2))*pow2(MuSUSY) + 20*
        pow2(Mst1)*pow4(Mst2) - 10*pow2(Mst2)*(pow4(Mst1) + pow4(Mst2))) + (8*(
        9 + 10*lmMst1 - 12*lmMst2)*pow2(Mst1)*pow2(Mst2) + (8 - 80*lmMst1 + 96*
        lmMst2)*pow4(Mst1) - 40*pow4(Mst2))*pow4(s2t)*pow5(Mst2) + 128*(-2 + 3*
        lmMst2)*Mt*pow2(Mst1)*pow3(s2t)*pow6(Mst2)))) + Tbeta*(s2t*pow2(Mt)*
        pow2(MuSUSY)*(256*(Dmglst2 + 3*Dmglst2*lmMst2 + Mst2 + lmMst2*Mst2)*Mt*
        pow2(Mst1)*(2*pow2(Mst1)*pow2(Mst2) + 3*pow4(Mst1) + pow4(Mst2)) + 8*
        Mst2*s2t*(2*(Dmglst2*(7*lmMst1 - 15*lmMst2) + (-4 - 13*lmMst1 + 9*
        lmMst2)*Mst2)*(pow2(Mst1) + pow2(Mst2))*pow4(Mst1) - Dmglst2*((7 - 14*
        lmMst1 + 30*lmMst2)*pow2(Mst1) + 7*pow2(Mst2))*pow4(Mst2) + (5 - 26*
        lmMst1 + 18*lmMst2)*pow2(Mst1)*pow5(Mst2) + 13*pow7(Mst2))) + pow2(
        Sbeta)*(2*(pow2(Mst1) - pow2(Mst2))*pow4(s2t)*(2*(4 + 13*lmMst1 - 9*
        lmMst2)*pow2(Mst1)*pow3(Mst2) + 13*Mst2*pow4(Mst1) + Dmglst2*(2*(-7*
        lmMst1 + 15*lmMst2)*pow2(Mst1)*pow2(Mst2) - 7*pow4(Mst1) + 7*pow4(Mst2)
        ) - 13*pow5(Mst2))*pow5(Mst2) + 32*Mst2*pow4(Mt)*(48*(1 + lmMst2)*Mst2*
        (2*pow2(Mst1) + pow2(Mst2))*pow4(Mst1) + (29 + 16*lmMst2)*pow2(Mst1)*
        pow5(Mst2) - Dmglst2*(96*(2 + lmMst2)*pow2(Mst2)*pow4(Mst1) + 39*pow2(
        Mst1)*pow4(Mst2) + 192*(3 + 2*lmMst2)*pow6(Mst1) + 7*pow6(Mst2)) + 13*
        pow7(Mst2)) - 256*s2t*pow2(Mst1)*pow3(Mt)*(-6*Dmglst2*(5 + 3*lmMst2)*
        pow2(Mst2)*pow4(Mst1) + 4*(-(Dmglst2*(3 + lmMst2)) + (1 + lmMst2)*Mst2)
        *pow2(Mst1)*pow4(Mst2) + (Dmglst2 + 3*Dmglst2*lmMst2 + Mst2 + lmMst2*
        Mst2)*pow2(MuSUSY)*(2*pow2(Mst1)*pow2(Mst2) + 3*pow4(Mst1) + pow4(Mst2)
        ) + 2*Dmglst2*(-1 + lmMst2)*pow6(Mst2) + 2*(1 + lmMst2)*(3*pow3(Mst2)*
        pow4(Mst1) + pow7(Mst2))) + 128*(Dmglst2 + 3*Dmglst2*lmMst2 + Mst2 +
        lmMst2*Mst2)*Mt*pow2(Mst1)*pow3(s2t)*pow8(Mst2) + 8*Mst2*pow2(Mt)*pow2(
        s2t)*(2*(7*Dmglst2 - 13*Mst2)*pow4(Mst1)*pow4(Mst2) + 4*(-7*Dmglst2 +
        13*Mst2)*pow2(Mst1)*pow6(Mst2) + pow2(MuSUSY)*(2*(4 + 13*lmMst1 - 9*
        lmMst2)*Mst2*(pow2(Mst1) + pow2(Mst2))*pow4(Mst1) + 7*Dmglst2*(pow2(
        Mst1) + pow2(Mst2))*pow4(Mst2) + Dmglst2*(-14*lmMst1 + 30*lmMst2)*pow2(
        Mst1)*(pow2(Mst1)*pow2(Mst2) + pow4(Mst1) + pow4(Mst2)) + (-5 + 26*
        lmMst1 - 18*lmMst2)*pow2(Mst1)*pow5(Mst2) - 13*pow7(Mst2)) + 2*(7*
        Dmglst2 - 13*Mst2)*pow8(Mst2)))))) + pow3(Mst2)*(-68040*twoLoopFlag*
        pow2(Mst1)*(pow5(Mst2)*(4*Mst2*MuSUSY*pow2(Mt)*pow2(s2t)*(-(MuSUSY*
        Tbeta*(-1 + pow2(Sbeta))) + 6*Mst2*pow2(Sbeta)) + 16*s2t*Tbeta*(pow2(
        MuSUSY)*(-1 + pow2(Sbeta)) - 4*pow2(Mst2)*pow2(Sbeta))*pow3(Mt) - 4*Mt*
        (MuSUSY + 2*Mst2*Tbeta)*pow2(Sbeta)*pow3(Mst2)*pow3(s2t) + 32*(2*MuSUSY
        + Mst2*Tbeta)*pow2(Sbeta)*pow4(Mt) + Tbeta*pow2(Sbeta)*pow4(s2t)*pow5(
        Mst2)) + Mst2*pow4(Mst1)*(16*s2t*Tbeta*(pow2(MuSUSY)*(-1 + pow2(Sbeta))
        - 4*pow2(Mst2)*pow2(Sbeta))*pow3(Mt) + 96*(2*MuSUSY + Mst2*Tbeta)*pow2(
        Sbeta)*pow4(Mt) + Tbeta*pow2(Sbeta)*pow4(s2t)*pow5(Mst2)) - 2*pow2(
        Mst1)*pow3(Mst2)*(pow2(Sbeta)*(-2*Mt*(MuSUSY + 2*Mst2*Tbeta)*pow3(Mst2)
        *pow3(s2t) - 32*(2*MuSUSY + Mst2*Tbeta)*pow4(Mt)) + Tbeta*(8*s2t*(pow2(
        MuSUSY) + (4*pow2(Mst2) - pow2(MuSUSY))*pow2(Sbeta))*pow3(Mt) + pow2(
        Sbeta)*pow4(s2t)*pow5(Mst2))) + 4*Dmglst2*(4*pow2(Mt)*(-12*(MuSUSY + 2*
        Mst2*Tbeta)*pow2(Mt)*pow2(Sbeta) + Mst2*MuSUSY*pow2(s2t)*(-(MuSUSY*
        Tbeta*(-1 + pow2(Sbeta))) + 6*Mst2*pow2(Sbeta)) + Mt*s2t*Tbeta*(13*
        pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 12*pow2(Mst2)*pow2(Sbeta)))*pow4(
        Mst1) + pow2(Mst1)*pow2(Mst2)*(8*s2t*(3*MuSUSY*s2t + 2*Mt*Tbeta)*pow2(
        Mst2)*pow2(Mt)*pow2(Sbeta) + Tbeta*(-4*Mst2*pow2(Mt)*pow2(MuSUSY)*pow2(
        s2t)*(-1 + pow2(Sbeta)) + 36*s2t*pow2(MuSUSY)*(-1 + pow2(Sbeta))*pow3(
        Mt) + (2*Mt - Mst2*s2t)*pow2(Sbeta)*pow3(s2t)*pow4(Mst2) - 32*Mst2*
        pow2(Sbeta)*pow4(Mt))) + pow4(Mst2)*(2*Mst2*MuSUSY*pow2(Mt)*pow2(s2t)*(
        -2*MuSUSY*Tbeta*(-1 + pow2(Sbeta)) + 15*Mst2*pow2(Sbeta)) + 4*s2t*
        Tbeta*(5*pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 4*pow2(Mst2)*pow2(Sbeta))*
        pow3(Mt) - 2*Mt*(2*MuSUSY + 5*Mst2*Tbeta)*pow2(Sbeta)*pow3(Mst2)*pow3(
        s2t) + 16*MuSUSY*pow2(Sbeta)*pow4(Mt) + Tbeta*pow2(Sbeta)*pow4(s2t)*
        pow5(Mst2)))) + Al4p*threeLoopFlag*(Tbeta*(-56*pow2(Mst1)*pow2(MuSUSY)*
        pow3(Mt)*(-27*pow2(Mst2)*(5*Mst2*(64*(53 + 24*lmMst2)*Mst2*Mt + 2880*
        Dmsqst2*(2 + lmMst1 - lmMst2)*s2t + 3*(2269 + 664*lmMst1 - 56*lmMst2)*
        s2t*pow2(Mst2)) + Dmglst2*(7680*(5 + 6*lmMst2)*Mst2*Mt - 14400*Dmsqst2*
        (2 + 3*lmMst1 - 3*lmMst2)*s2t + (173947 - 25080*lmMst1 + 68760*lmMst2)*
        s2t*pow2(Mst2))) - 36*pow2(Mst1)*(5*Mst2*(3*(-353 + 72*lmMst1 + 696*
        lmMst2)*Mst2*Mt + 2160*Dmsqst2*(5 + 3*lmMst1 - 3*lmMst2)*s2t + 2*(3937
        + 2988*lmMst1 - 2232*lmMst2)*s2t*pow2(Mst2)) + Dmglst2*(30*(3643 - 120*
        lmMst1 + 3192*lmMst2)*Mst2*Mt + 10800*Dmsqst2*(1 - 9*lmMst1 + 9*lmMst2)
        *s2t + (364291 - 88560*lmMst1 + 147960*lmMst2)*s2t*pow2(Mst2))) + (90*(
        38401 + 1080*lmMst1 - 7992*lmMst2)*Mt + 2*Dmglst2*(-15057833 + 4014360*
        lmMst1 - 5563080*lmMst2)*s2t - 6*(266863 + 396360*lmMst1 - 346680*
        lmMst2)*Mst2*s2t)*pow4(Mst1)) - 68040*shiftst3*pow2(Mst2)*(-4*pow2(
        Mst1)*pow2(s2t)*pow4(Mst2)*(2*(-1 + lmMst1 - lmMst2)*pow2(Mt)*(-(pow2(
        MuSUSY)*(-1 + pow2(Sbeta))) + pow2(Mst2)*pow2(Sbeta)) + pow2(s2t)*pow2(
        Sbeta)*pow4(Mst2)) + 8*(-1 + 3*lmMst1 - 3*lmMst2)*pow2(Mt)*pow2(MuSUSY)
        *pow2(s2t)*(-1 + pow2(Sbeta))*pow6(Mst1) + (-4*pow2(Mt)*pow2(s2t)*(
        pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 2*pow2(Mst2)*pow2(Sbeta)) + pow2(
        Sbeta)*(16*pow4(Mt) + pow4(Mst2)*pow4(s2t)))*pow6(Mst2) + pow4(Mst1)*(
        8*(-1 + 2*lmMst1 - 2*lmMst2)*pow2(Mst2)*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)
        *(-1 + pow2(Sbeta)) + (3 + 2*lmMst1 - 2*lmMst2)*pow2(Sbeta)*pow4(s2t)*
        pow6(Mst2))) - 28*pow2(Mt)*pow2(s2t)*(-300*Dmsqst2*(9*Mst2*pow2(MuSUSY)
        *(159*pow2(Mst1)*pow3(Mst2) + (1097 - 72*lmMst1 + 72*lmMst2)*Mst2*pow4(
        Mst1) + 72*Dmglst2*(3*(2 + lmMst1 - lmMst2)*pow2(Mst1)*pow2(Mst2) + (3
        + 8*lmMst1 - 8*lmMst2)*pow4(Mst1) + pow4(Mst2)) - 36*pow5(Mst2)) + ((
        20701 - 648*lmMst1 + 648*lmMst2)*pow2(MuSUSY) + 162*Dmglst2*(-93 + 10*
        lmMst1 - 10*lmMst2)*Mst2*(-2 + pow2(Sbeta))*pow4(Sbeta))*pow6(Mst1)) +
        Mst2*pow2(MuSUSY)*(4*Dmglst2*(162*(6803 - 1810*lmMst1 + 2670*lmMst2)*
        pow2(Mst2)*pow4(Mst1) + 135*(1631 - 648*lmMst1 + 1680*lmMst2)*pow2(
        Mst1)*pow4(Mst2) + (4256978 - 615600*lmMst1 + 754920*lmMst2)*pow6(Mst1)
        - 22680*pow6(Mst2)) - 3*(6*(533629 - 900*lmMst1 + 180*lmMst2)*pow3(
        Mst2)*pow4(Mst1) + 90*(5597 + 564*lmMst1 - 1104*lmMst2)*pow2(Mst1)*
        pow5(Mst2) + (6649153 - 81540*lmMst1 + 77220*lmMst2)*Mst2*pow6(Mst1) -
        20520*pow7(Mst2)))) - 680400*shiftst1*(Dmsqst2 + pow2(Mst2))*(pow4(
        Mst2)*(pow2(s2t)*(-8*pow2(Mt) + (pow2(Mst1) - pow2(Mst2))*pow2(s2t))*
        pow2(Sbeta)*pow4(Mst1) + pow2(Mst1)*(4*pow2(Mt)*pow2(s2t)*(pow2(MuSUSY)
        + (4*pow2(Mst2) - pow2(MuSUSY))*pow2(Sbeta)) + pow2(Sbeta)*(16*pow4(Mt)
        - pow4(Mst2)*pow4(s2t))) + pow2(Mst2)*(-4*pow2(Mt)*pow2(s2t)*(pow2(
        MuSUSY)*(-1 + pow2(Sbeta)) + 2*pow2(Mst2)*pow2(Sbeta)) + pow2(Sbeta)*(
        16*pow4(Mt) + pow4(Mst2)*pow4(s2t)))) + (-2*lmMst1 + 2*lmMst2)*pow2(
        Mst1)*pow2(s2t)*(-4*pow2(Mt)*pow2(MuSUSY)*(-1 + pow2(Sbeta))*(pow4(
        Mst1) + pow4(Mst2)) - pow2(Mst1)*pow2(Mst2)*(4*pow2(Mt)*pow2(MuSUSY)*(-
        1 + pow2(Sbeta)) + pow2(s2t)*pow2(Sbeta)*pow4(Mst2)) + pow2(s2t)*pow2(
        Sbeta)*pow8(Mst2))) + pow2(Sbeta)*(-56*s2t*pow2(Mst1)*pow3(Mt)*(
        Dmglst2*(18*pow2(Mst1)*pow2(Mst2)*(3*(510149 - 55320*lmMst1 + 118200*
        lmMst2 - 32160*lmMt)*pow2(Mst2) + 2*(364291 - 88560*lmMst1 + 147960*
        lmMst2)*pow2(MuSUSY)) + 2*(4*(9908167 - 949320*lmMst1 + 1769040*lmMst2
        - 217080*lmMt)*pow2(Mst2) + (15057833 - 4014360*lmMst1 + 5563080*
        lmMst2)*pow2(MuSUSY))*pow4(Mst1) + 27*(30*(-1672*lmMst1 + 1448*lmMst2 +
        59*(71 - 32*lmMt))*pow2(Mst2) + (173947 - 25080*lmMst1 + 68760*lmMst2)*
        pow2(MuSUSY))*pow4(Mst2) + 388800*Dmsqst2*((-2 - 3*lmMst1 + 3*lmMst2)*
        pow2(Mst2)*pow2(MuSUSY) + pow2(Mst1)*((51 - 14*lmMst1 + 28*lmMst2 - 14*
        lmMt)*pow2(Mst2) + (1 - 9*lmMst1 + 9*lmMst2)*pow2(MuSUSY)) + (87 - 22*
        lmMst1 + 44*lmMst2 - 22*lmMt)*pow4(Mst1) + (19 - 6*lmMst1 + 12*lmMst2 -
        6*lmMt)*pow4(Mst2))) - 3*Mst2*(6*pow2(Mst1)*pow2(Mst2)*((470657 -
        86040*lmMst1 + 178200*lmMst2)*pow2(Mst2) - 20*(3937 + 2988*lmMst1 -
        2232*lmMst2)*pow2(MuSUSY)) + 2*((2299036 - 388800*lmMst1 + 656640*
        lmMst2)*pow2(Mst2) + (-266863 - 396360*lmMst1 + 346680*lmMst2)*pow2(
        MuSUSY))*pow4(Mst1) + 9*(2*(68693 - 9960*lmMst1 + 45480*lmMst2 - 3840*
        lmMt)*pow2(Mst2) - 15*(2269 + 664*lmMst1 - 56*lmMst2)*pow2(MuSUSY))*
        pow4(Mst2) + 129600*Dmsqst2*((-2 - lmMst1 + lmMst2)*pow2(Mst2)*pow2(
        MuSUSY) + pow2(Mst1)*((5 - 2*lmMst1 + 4*lmMst2 - 2*lmMt)*pow2(Mst2) + (
        -5 - 3*lmMst1 + 3*lmMst2)*pow2(MuSUSY)) + (5 - 2*lmMst1 + 4*lmMst2 - 2*
        lmMt)*pow4(Mst1) + (1 - 2*lmMst1 + 4*lmMst2 - 2*lmMt)*pow4(Mst2)))) +
        7*pow3(Mst2)*pow4(s2t)*(-4*Dmglst2*(-108*(-5917 + 1095*lmMst1 + 195*
        lmMst2)*pow4(Mst1)*pow4(Mst2) + (2272991 - 116640*lmMst1 + 116640*
        lmMst2)*pow2(Mst2)*pow6(Mst1) - 135*(648*lmMst1 - 7*(281 + 240*lmMst2))
        *pow2(Mst1)*pow6(Mst2) + 24300*Dmsqst2*(4*(4 - lmMst1 + lmMst2)*pow2(
        Mst2)*pow4(Mst1) + (7 + 10*lmMst1 - 10*lmMst2)*pow6(Mst1) - 2*((4 + 3*
        lmMst1 - 3*lmMst2)*pow2(Mst1)*pow4(Mst2) + pow6(Mst2))) - 22680*pow8(
        Mst2)) + 3*Mst2*(6*(362299 - 17820*lmMst1 + 33300*lmMst2)*pow4(Mst1)*
        pow4(Mst2) - 5*(-149867 + 3996*lmMst1 + 4860*lmMst2)*pow2(Mst2)*pow6(
        Mst1) + 100*Dmsqst2*(9*(743 - 72*lmMst1 + 72*lmMst2)*pow2(Mst2)*pow4(
        Mst1) + 2079*pow2(Mst1)*pow4(Mst2) + (2386 + 648*lmMst1 - 648*lmMst2)*
        pow6(Mst1) - 324*pow6(Mst2)) + 90*(6053 + 564*lmMst1 - 1104*lmMst2)*
        pow2(Mst1)*pow6(Mst2) - 20520*pow8(Mst2))) - 16*pow4(Mt)*(680400*
        Dmsqst2*(-2*((Dmglst2*(23 - 5*lmMst1 + 10*lmMst2 - 5*lmMt) + (-3 +
        lmMst1 - 2*lmMst2 + lmMt)*Mst2)*pow2(Mst1)*pow3(Mst2) + Mst2*(Dmglst2*(
        77 - 18*lmMst1 + 36*lmMst2 - 18*lmMt) + 2*(-3 + lmMst1 - 2*lmMst2 +
        lmMt)*Mst2)*pow4(Mst1)) + (-2*Dmglst2 + Mst2)*pow5(Mst2) + 2*(7 - 3*
        lmMst1 + 6*lmMst2 - 3*lmMt)*pow6(Mst1)) - 4*Dmglst2*Mst2*(27*pow2(Mst1)
        *pow2(Mst2)*((209341 - 26880*lmMst1 + 82320*lmMst2 - 1680*lmMt)*pow2(
        Mst2) - 6720*(5 + 6*lmMst2)*pow2(MuSUSY)) + 9*((3184397 - 476280*lmMst1
        + 1156680*lmMst2 - 161280*lmMt)*pow2(Mst2) + 105*(-3643 + 120*lmMst1 -
        3192*lmMst2)*pow2(MuSUSY))*pow4(Mst1) + 4*(19152737 - 3674160*lmMst1 +
        6872040*lmMst2 - 850500*lmMt)*pow6(Mst1) + 158760*pow6(Mst2)) + 21*(3*
        pow2(Mst2)*((310039 - 91080*lmMst1 + 145080*lmMst2 - 43920*lmMt)*pow2(
        Mst2) + 30*(-353 + 72*lmMst1 + 696*lmMst2)*pow2(MuSUSY))*pow4(Mst1) +
        9*pow2(Mst1)*((48983 - 12480*lmMst1 + 26820*lmMst2 - 12180*lmMt)*pow2(
        Mst2) + 160*(53 + 24*lmMst2)*pow2(MuSUSY))*pow4(Mst2) + (2*(626869 -
        300240*lmMst1 + 395280*lmMst2 - 78840*lmMt)*pow2(Mst2) - 15*(38401 +
        1080*lmMst1 - 7992*lmMst2)*pow2(MuSUSY))*pow6(Mst1) + 20520*pow8(Mst2))
        ) - 28*Mt*pow2(s2t)*(s2t*pow2(Mst1)*pow2(Mst2)*(51*(6167 - 9720*lmMst1
        + 9720*lmMst2)*pow3(Mst2)*pow4(Mst1) + 194400*Dmsqst2*Mst2*(-2*(1 +
        lmMst1 - lmMst2)*pow2(Mst1)*pow2(Mst2) + pow4(Mst1) - 2*(2 + lmMst1 -
        lmMst2)*pow4(Mst2)) + 90*(4673 - 5976*lmMst1 + 8424*lmMst2)*pow2(Mst1)*
        pow5(Mst2) - Dmglst2*((8583283 - 2329560*lmMst1 + 2329560*lmMst2)*pow2(
        Mst2)*pow4(Mst1) + 18*(206741 - 101880*lmMst1 + 89640*lmMst2)*pow2(
        Mst1)*pow4(Mst2) + 194400*Dmsqst2*(2*(5 - 3*lmMst1 + 3*lmMst2)*pow2(
        Mst1)*pow2(Mst2) + 11*pow4(Mst1) - 2*(2 + 3*lmMst1 - 3*lmMst2)*pow4(
        Mst2)) + 27*(173947 - 25080*lmMst1 + 68760*lmMst2)*pow6(Mst2)) - 405*(
        2269 + 664*lmMst1 - 56*lmMst2)*pow7(Mst2)) + Mt*(4*Dmglst2*Mst2*(-162*
        pow2(Mst2)*(20*(505 + 68*lmMst1 + 124*lmMst2)*pow2(Mst2) + (6803 -
        1810*lmMst1 + 2670*lmMst2)*pow2(MuSUSY))*pow4(Mst1) - 27*pow2(Mst1)*(4*
        (4849 - 270*lmMst1 + 8910*lmMst2)*pow2(Mst2) + 5*(1631 - 648*lmMst1 +
        1680*lmMst2)*pow2(MuSUSY))*pow4(Mst2) - 2*(45*(78533 + 6732*lmMst1 +
        180*lmMst2)*pow2(Mst2) + (2128489 - 307800*lmMst1 + 377460*lmMst2)*
        pow2(MuSUSY))*pow6(Mst1) + 48600*Dmsqst2*((8*pow2(Mst2) + (3 + 8*lmMst1
        - 8*lmMst2)*pow2(MuSUSY))*pow4(Mst1) + (2*pow2(Mst2) + pow2(MuSUSY))*
        pow4(Mst2) + pow2(Mst1)*(3*(2 + lmMst1 - lmMst2)*pow2(Mst2)*pow2(
        MuSUSY) + 8*pow4(Mst2)) + 8*pow6(Mst1)) + 22680*(2*pow2(Mst2) + pow2(
        MuSUSY))*pow6(Mst2)) + 3*(6*((309749 + 12240*lmMst1 - 12240*lmMst2)*
        pow2(Mst2) + (533629 - 900*lmMst1 + 180*lmMst2)*pow2(MuSUSY))*pow4(
        Mst1)*pow4(Mst2) + pow2(Mst2)*((3291029 + 210600*lmMst1 - 210600*
        lmMst2)*pow2(Mst2) + (6649153 - 81540*lmMst1 + 77220*lmMst2)*pow2(
        MuSUSY))*pow6(Mst1) + 9*pow2(Mst1)*((5927 + 13560*lmMst1 - 36600*
        lmMst2)*pow2(Mst2) + 10*(5597 + 564*lmMst1 - 1104*lmMst2)*pow2(MuSUSY))
        *pow6(Mst2) + 100*Dmsqst2*(9*pow2(Mst2)*(202*pow2(Mst2) + (1097 - 72*
        lmMst1 + 72*lmMst2)*pow2(MuSUSY))*pow4(Mst1) + (808*pow2(Mst2) + (20701
        - 648*lmMst1 + 648*lmMst2)*pow2(MuSUSY))*pow6(Mst1) - 324*(2*pow2(Mst2)
        + pow2(MuSUSY))*pow6(Mst2) + 27*pow2(Mst1)*(53*pow2(MuSUSY)*pow4(Mst2)
        + 77*pow6(Mst2))) - 20520*(2*pow2(Mst2) + pow2(MuSUSY))*pow8(Mst2))))))
        + 28*Mt*MuSUSY*pow2(Sbeta)*(-3*(3*pow3(Mst2)*(-2*(623599 + 65160*lmMst1
        - 134280*lmMst2)*Mst2*s2t*pow2(Mt) + 15*(11075 + 17928*lmMst1 - 17352*
        lmMst2)*Mt*pow2(Mst2)*pow2(s2t) + 64*(41936 - 7245*lmMst1 + 19665*
        lmMst2 - 720*lmMt)*pow3(Mt) + 4*(224837 - 4680*lmMst1 + 8370*lmMst2 +
        2700*shiftst1)*pow3(Mst2)*pow3(s2t))*pow4(Mst1) + 9*pow2(Mst1)*(-2*
        Mst2*s2t*(1367 + 13560*lmMst1 - 36600*lmMst2 + 7200*shiftst1)*pow2(Mt)
        + 45*(2269 + 664*lmMst1 - 56*lmMst2)*Mt*pow2(Mst2)*pow2(s2t) + 4*(66773
        - 9960*lmMst1 + 45480*lmMst2 - 3840*lmMt)*pow3(Mt) + 10*(5825 + 36*
        shiftst3 + 12*lmMst1*(47 + 60*shiftst1 + 3*shiftst3) - 12*lmMst2*(92 +
        60*shiftst1 + 3*shiftst3))*pow3(Mst2)*pow3(s2t))*pow5(Mst2) + Mst2*(-4*
        (2580913 + 203040*lmMst1 - 306720*lmMst2)*Mst2*s2t*pow2(Mt) + 6*(30643
        + 217080*lmMst1 - 212760*lmMst2)*Mt*pow2(Mst2)*pow2(s2t) + 16*(1077991
        - 184140*lmMst1 + 400140*lmMst2 - 8640*lmMt)*pow3(Mt) + (3447379 -
        76140*lmMst1 + 76140*lmMst2)*pow3(Mst2)*pow3(s2t))*pow6(Mst1) - 100*
        Dmsqst2*(27*pow2(Mst1)*pow3(Mst2)*(106*Mst2*s2t*pow2(Mt) + 48*Mst2*s2t*
        shiftst1*pow2(Mt) - 288*Mt*pow2(Mst2)*pow2(s2t) + 144*lmMst2*Mt*pow2(
        Mst2)*pow2(s2t) - 384*lmMst2*pow3(Mt) + 192*lmMt*pow3(Mt) - 65*pow3(
        Mst2)*pow3(s2t) + 24*lmMst2*shiftst1*pow3(Mst2)*pow3(s2t) + 24*lmMst1*(
        -6*Mt*pow2(Mst2)*pow2(s2t) + 8*pow3(Mt) - shiftst1*pow3(Mst2)*pow3(s2t)
        )) + 18*Mst2*(361*Mst2*s2t*pow2(Mt) - 216*(3 + 2*lmMst1 - 2*lmMst2)*Mt*
        pow2(Mst2)*pow2(s2t) + 576*(-1 + lmMst1 - 2*lmMst2 + lmMt)*pow3(Mt) + (
        -469 + 36*lmMst1 - 36*lmMst2 - 18*shiftst1)*pow3(Mst2)*pow3(s2t))*pow4(
        Mst1) + 2*s2t*(-972*(5 + 4*lmMst1 - 4*lmMst2)*Mst2*Mt*s2t + 4057*pow2(
        Mt) - 5414*pow2(Mst2)*pow2(s2t))*pow6(Mst1) + 324*s2t*(1 + shiftst1)*(-
        4*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow6(Mst2)) - 1080*s2t*(19 + 30*
        shiftst1 + 3*shiftst3)*(-4*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow8(Mst2))
        + Dmglst2*(27*(pow2(Mst2)*(-32*(19579 + 1770*lmMst1 + 12630*lmMst2)*
        Mst2*s2t*pow2(Mt) + (-935323 + 279000*lmMst1 - 385560*lmMst2)*Mt*pow2(
        Mst2)*pow2(s2t) + 32*(67843 - 10050*lmMst1 + 17490*lmMst2 - 7560*lmMt)*
        pow3(Mt) + 4*(32663 - 7620*lmMst1 + 7620*lmMst2)*pow3(Mst2)*pow3(s2t))*
        pow4(Mst1) + pow2(Mst1)*(32*(-4429 + 270*lmMst1 - 8910*lmMst2)*Mst2*
        s2t*pow2(Mt) + 3*(-173947 + 25080*lmMst1 - 68760*lmMst2)*Mt*pow2(Mst2)*
        pow2(s2t) - 60*(-3245 + 1672*lmMst1 - 1448*lmMst2 + 1888*lmMt)*pow3(Mt)
        + 20*(1799 - 648*lmMst1 + 1680*lmMst2)*pow3(Mst2)*pow3(s2t))*pow4(Mst2)
        ) + 2*(-72*(510139 + 44280*lmMst1 + 76680*lmMst2)*Mst2*s2t*pow2(Mt) +
        15*(-1700119 + 484056*lmMst1 - 579960*lmMst2)*Mt*pow2(Mst2)*pow2(s2t) +
        8*(13463149 - 1492020*lmMst1 + 2713500*lmMst2 - 625320*lmMt)*pow3(Mt) +
        8*(788723 - 80595*lmMst1 + 80595*lmMst2)*pow3(Mst2)*pow3(s2t))*pow6(
        Mst1) + 97200*Dmsqst2*(2*pow2(Mst1)*pow2(Mst2)*(20*Mst2*s2t*pow2(Mt) +
        6*(2 + 3*lmMst1 - 3*lmMst2)*Mt*pow2(Mst2)*pow2(s2t) - 8*(-8 + 3*lmMst1
        - 6*lmMst2 + 3*lmMt)*pow3(Mt) + (-5 - 3*lmMst1 + 3*lmMst2)*pow3(Mst2)*
        pow3(s2t)) - 2*(-36*Mst2*s2t*pow2(Mt) + 18*(1 - 2*lmMst1 + 2*lmMst2)*
        Mt*pow2(Mst2)*pow2(s2t) + 80*(-3 + lmMst1 - 2*lmMst2 + lmMt)*pow3(Mt) +
        (-3 + 5*lmMst1 - 5*lmMst2)*pow3(Mst2)*pow3(s2t))*pow4(Mst1) + 8*s2t*
        pow2(Mt)*pow5(Mst2) + (6*(-17 + 12*lmMst1 - 12*lmMst2)*Mt + 13*Mst2*
        s2t)*pow2(s2t)*pow6(Mst1) - 2*pow3(s2t)*pow7(Mst2)) - 90720*(-4*s2t*
        pow2(Mt)*pow7(Mst2) + pow3(s2t)*pow9(Mst2))))))))/(408240.*Tbeta*pow2(
        Mst1)*pow2(Sbeta)))/pow9(Mst2) + pow2(Al4p)*(-(threeLoopFlag*(xDmsqst2*
        pow2(Dmsqst2)*pow2(Mst1)*(16*Mt*MuSUSY*pow2(Sbeta)*((45*Mst2*pow2(Mst1)
        *(16*Mst2*s2t*(88 + 8*OepS2 - 81*S2 - 18*(lmMst1 - lmMst2)*(2 + 9*S2) -
        54*shiftst1 - 12*T1ep)*pow2(Mt) + 2592*(3 + 2*lmMst1 - 2*lmMst2)*Mt*
        pow2(Mst2)*pow2(s2t) + 864*(4 - 8*lmMst1 + 15*lmMst2 - 7*lmMt)*pow3(Mt)
        - (1493 + 64*OepS2 - 27864*S2 - 72*(lmMst1 - lmMst2)*(-1 + 18*S2 + 6*
        shiftst1) - 96*T1ep)*pow3(Mst2)*pow3(s2t)))/4. + (3*s2t*(972000*(1 + 4*
        lmMst1 - 4*lmMst2)*Mst2*Mt*s2t + 8*(187757 + 16000*OepS2 - 540000*S2 +
        90*lmMst2*(-587 + 3600*S2) - 90*lmMst1*(-587 + 630*lmMst2 + 3600*S2) -
        24000*T1ep + 28350*(pow2(lmMst1) + pow2(lmMst2)))*pow2(Mt) + (1372861 -
        88000*OepS2 + 180*lmMst2*(1987 - 9900*S2) + 13743000*S2 + 180*lmMst1*(-
        1987 + 1080*lmMst2 + 9900*S2) + 81000*shiftst1 + 132000*T1ep - 97200*(
        pow2(lmMst1) + pow2(lmMst2)))*pow2(Mst2)*pow2(s2t))*pow4(Mst1))/100. -
        2430*s2t*(2 + shiftst1)*(-4*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow4(Mst2)
        + ((610359812 - 64484000*OepS2 + 315*lmMst2*(540413 - 4145400*S2) +
        8117266500*S2 + 315*lmMst1*(-540413 + 884520*lmMst2 + 4145400*S2) +
        96726000*T1ep - 139311900*(pow2(lmMst1) + pow2(lmMst2)))*pow3(s2t)*
        pow6(Mst1))/17150.) + Tbeta*(-311040*(3 + 2*lmMst1 - 2*lmMst2)*Mst2*
        s2t*pow2(Mst1)*pow2(MuSUSY)*pow3(Mt) - 144*pow2(Mt)*pow2(s2t)*(270*
        shiftst1*(2*(lmMst1 - lmMst2)*pow2(Mst1) - pow2(Mst2))*(pow2(Mst1) +
        pow2(Mst2))*pow2(MuSUSY) - 40*T1ep*pow2(Mst1)*pow2(Sbeta)*(3*pow2(Mst2)
        *(pow2(Mst2) + pow2(MuSUSY)) + pow2(Mst1)*(5*pow2(Mst2) + 14*pow2(
        MuSUSY)) + 5*pow4(Mst1))) + 1944*shiftst1*pow2(Sbeta)*(20*pow2(Mt)*
        pow2(s2t)*((2*(lmMst1 - lmMst2)*pow2(Mst1) - pow2(Mst2))*(pow2(Mst1) +
        pow2(Mst2))*pow2(MuSUSY) + 4*pow2(Mst1)*pow4(Mst2) - 2*pow2(Mst2)*(
        pow4(Mst1) + pow4(Mst2))) + 80*pow2(Mst2)*(pow2(Mst1) + pow2(Mst2))*
        pow4(Mt) + 5*(pow2(Mst1) - pow2(Mst2))*pow2(Mst2)*(2*(lmMst1 - lmMst2)*
        pow2(Mst1)*pow2(Mst2) + pow4(Mst1) - pow4(Mst2))*pow4(s2t)) - 480*T1ep*
        pow2(Mst1)*(12*(14*pow2(Mst1) + 3*pow2(Mst2))*pow2(Mt)*pow2(MuSUSY)*
        pow2(s2t) + pow2(Mst2)*pow2(Sbeta)*(24*pow2(Mst1)*pow2(Mst2) + 14*pow4(
        Mst1) + 9*pow4(Mst2))*pow4(s2t)) + pow2(s2t)*(-4860*Mt*(2*(11 + 39*
        lmMst1 - 39*lmMst2)*Mst2*s2t*(-2 + pow2(Sbeta)) + (167 - 5*lmMst1 + 5*
        lmMst2)*Mt*pow2(Sbeta))*pow4(Sbeta)*pow6(Mst1) + (12*pow2(Mt)*(pow2(
        MuSUSY)*(375*(1925 + 64*OepS2 - 27864*S2 - 72*(lmMst1 - lmMst2)*(-1 +
        18*S2))*pow2(Mst1)*pow2(Mst2) + 2*(7*(-46499 + 8000*OepS2 - 1728000*S2)
        + 90*lmMst2*(-2137 + 12600*S2) - 90*lmMst1*(-2137 + 1080*lmMst2 +
        12600*S2) + 48600*(pow2(lmMst1) + pow2(lmMst2)))*pow4(Mst1) + 162000*
        pow4(Mst2)) + 20250*(167 - 5*lmMst1 + 5*lmMst2)*pow4(Sbeta)*pow6(Mst1))
        )/25.) + pow2(Sbeta)*(-77760*Mst2*s2t*pow2(Mst1)*(4*(11 - 6*lmMst1 +
        12*lmMst2 - 6*lmMt)*pow2(Mst1) + (15 - 16*lmMst1 + 30*lmMst2 - 14*lmMt)
        *pow2(Mst2) - 4*(3 + 2*lmMst1 - 2*lmMst2)*pow2(MuSUSY))*pow3(Mt) +
        19440*Mst2*Mt*pow2(Mst1)*pow3(s2t)*(16*(1 - lmMst1 + lmMst2)*pow2(Mst1)
        *pow2(Mst2) + (29 + 41*lmMst1 - 41*lmMst2)*pow4(Mst1) - 8*(3 + 2*lmMst1
        - 2*lmMst2)*pow4(Mst2)) + (144*(125*(3095 - 1224*lmMst1 + 6*lmMst2*(365
        - 36*lmMt) - 966*lmMt + 108*pow2(lmMst2) + 108*pow2(lmMt))*pow2(Mst1)*
        pow2(Mst2) + 6*(176473 + 30*lmMst1*(-1771 + 90*lmMst2 - 150*lmMt) -
        63000*lmMt + 30*lmMst2*(3871 + 150*lmMt) + 900*pow2(lmMst1) - 3600*
        pow2(lmMst2))*pow4(Mst1) + 54000*pow4(Mst2))*pow4(Mt))/25. - (12*pow2(
        Mt)*pow2(s2t)*(686*(2*(121757 + 10000*OepS2 - 479250*S2 + 90*lmMst2*(-
        887 + 2250*S2) - 90*lmMst1*(-887 + 630*lmMst2 + 2250*S2) + 28350*(pow2(
        lmMst1) + pow2(lmMst2)))*pow2(Mst2) + (7*(-46499 + 8000*OepS2 -
        1728000*S2) + 90*lmMst2*(-2137 + 12600*S2) - 90*lmMst1*(-2137 + 1080*
        lmMst2 + 12600*S2) + 48600*(pow2(lmMst1) + pow2(lmMst2)))*pow2(MuSUSY))
        *pow4(Mst1) + 55566000*(2*pow2(Mst2) + pow2(MuSUSY))*pow4(Mst2) -
        128625*pow2(Mst1)*(-((1925 + 64*OepS2 - 27864*S2 - 72*(lmMst1 - lmMst2)
        *(-1 + 18*S2))*pow2(Mst2)*pow2(MuSUSY)) + 8*(20 - 8*OepS2 + 81*S2 + 18*
        (lmMst1 - lmMst2)*(2 + 9*S2))*pow4(Mst2)) + 250*(397737 + 54880*OepS2 -
        2580732*S2 + 252*lmMst2*(-2333 + 4410*S2) - 252*lmMst1*(-2333 + 1134*
        lmMst2 + 4410*S2) + 142884*(pow2(lmMst1) + pow2(lmMst2)))*pow6(Mst1)))/
        8575. + (pow2(Mst2)*pow4(s2t)*(-4116*(483184 - 16000*OepS2 + 45*lmMst2*
        (1837 - 7200*S2) + 823500*S2 + 45*lmMst1*(-1837 + 1080*lmMst2 + 7200*
        S2) - 24300*(pow2(lmMst1) + pow2(lmMst2)))*pow2(Mst2)*pow4(Mst1) +
        385875*(1061 + 64*OepS2 - 27864*S2 - 72*(lmMst1 - lmMst2)*(-1 + 18*S2))
        *pow2(Mst1)*pow4(Mst2) + 5*(38390869 + 7683200*OepS2 - 418597200*S2 +
        630*lmMst2*(8753 + 246960*S2) - 630*lmMst1*(8753 + 113400*lmMst2 +
        246960*S2) + 35721000*(pow2(lmMst1) + pow2(lmMst2)))*pow6(Mst1) +
        166698000*pow6(Mst2)))/8575.))) - 36*xDR2DRMOD*(4*Mt*MuSUSY*pow2(Sbeta)
        *(2*Dmglst2*(2*pow3(Mst2)*pow4(Mst1)*(-270*Dmsqst2*s2t*pow2(Mt) + 3*
        s2t*pow2(Mst2)*(45*Dmsqst2*(lmMst1 - lmMst2)*pow2(s2t) + 2*pow2(Mt)*(52
        + 370*lmMst2 - 96*lmMst2*pow2(lmMst1) + 81*pow2(lmMst2) + 96*lmMst1*(-1
        + lmMst2 + 2*pow2(lmMst2)) - 96*pow3(lmMst2))) - 432*Mt*pow2(s2t)*(7 +
        9*lmMst2 - 8*pow2(lmMst2) + lmMst1*(-3 + 5*lmMst2 + 6*pow2(lmMst2)) -
        6*pow3(lmMst2))*pow3(Mst2) + 32*Mst2*(119 + lmMst2*(224 - 15*lmMt) -
        51*lmMt + 3*(5 + 6*lmMt)*pow2(lmMst2) + 18*lmMst1*(-4 + lmMt - lmMst2*(
        1 + lmMt) + pow2(lmMst2)) - 18*pow3(lmMst2))*pow3(Mt) + 9*(40 + 4*
        lmMst2 + 16*lmMst2*pow2(lmMst1) + lmMst1*(4 - 246*lmMst2 - 123*pow2(
        lmMst2)) + 198*pow2(lmMst2) + 107*pow3(lmMst2))*pow3(s2t)*pow4(Mst2)) +
        9*pow2(Mst1)*(60*Dmsqst2*s2t*pow2(Mt) + 96*Mt*(-1 + 2*lmMst2 + 3*pow2(
        lmMst2))*pow2(s2t)*pow3(Mst2) - 128*Mst2*(-1 + 2*lmMst2 + 3*pow2(
        lmMst2))*pow3(Mt) - pow2(Mst2)*(4*s2t*(-52 + 102*lmMst2 + 32*lmMst1*
        lmMst2 + 59*pow2(lmMst2))*pow2(Mt) + 15*Dmsqst2*pow3(s2t)) + (-20 + (
        230 + 32*lmMst1)*lmMst2 + 155*pow2(lmMst2))*pow3(s2t)*pow4(Mst2))*pow5(
        Mst2) + pow2(Mst2)*(-1152*Mst2*s2t*(lmMst2*pow2(lmMst1) + lmMst1*(4 +
        2*lmMst2 - 2*pow2(lmMst2)) + (-4 + lmMst2)*pow2(1 + lmMst2))*pow2(Mt) -
        864*Mt*pow2(Mst2)*pow2(s2t)*(8 + 19*lmMst2 - 5*pow2(lmMst2) + lmMst1*(-
        7 + 5*lmMst2 + 6*pow2(lmMst2)) - 6*pow3(lmMst2)) + 64*(263 + lmMst2*(
        962 - 213*lmMt) - 177*lmMt + (393 - 18*lmMt)*pow2(lmMst2) - 18*lmMst1*(
        26 - lmMst2*(-17 + lmMt) - 7*lmMt + pow2(lmMst2)) + 18*pow3(lmMst2))*
        pow3(Mt) + 9*Mst2*(15*Dmsqst2 + pow2(Mst2)*(36 + 42*lmMst2 + 256*
        lmMst2*pow2(lmMst1) + 37*pow2(lmMst2) - 32*lmMst1*(3 + 3*lmMst2 + 16*
        pow2(lmMst2)) + 256*pow3(lmMst2)))*pow3(s2t))*pow6(Mst1) + 16*(-72*
        Mst2*s2t*pow2(Mt)*(-6 - 14*lmMst2 + lmMst2*pow2(lmMst1) + lmMst1*(9 +
        6*lmMst2 - 2*pow2(lmMst2)) - 6*pow2(lmMst2) + pow3(lmMst2)) + 27*Mt*
        pow2(Mst2)*pow2(s2t)*(-17 - 51*lmMst2 + 4*pow2(lmMst2) - 4*lmMst1*(-6 +
        lmMst2 + 3*pow2(lmMst2)) + 12*pow3(lmMst2)) + 4*(290 + lmMst2*(2447 -
        645*lmMt) - 375*lmMt - 3*(-509 + 60*lmMt)*pow2(lmMst2) - 18*lmMst1*(87
        + lmMst2*(71 - 10*lmMt) - 22*lmMt + 10*pow2(lmMst2)) + 180*pow3(lmMst2)
        )*pow3(Mt) + 9*(16*lmMst2*pow2(lmMst1) - 2*lmMst1*(3 + 2*lmMst2 + 16*
        pow2(lmMst2)) + lmMst2*(7 + 4*lmMst2 + 16*pow2(lmMst2)))*pow3(Mst2)*
        pow3(s2t))*pow8(Mst1) - 288*lmMst2*(1 + lmMst2)*s2t*(-4*pow2(Mt) +
        pow2(Mst2)*pow2(s2t))*pow9(Mst2)) + 9*Mst2*(-2*pow3(Mst2)*pow4(Mst1)*(
        30*Dmsqst2*s2t*(-5 + 9*shiftst1)*pow2(Mt) + s2t*pow2(Mst2)*(15*Dmsqst2*
        (lmMst1 - lmMst2)*(5 - 9*shiftst1)*pow2(s2t) - 2*pow2(Mt)*(109 + 76*
        lmMst2 - 135*shiftst1 - 21*pow2(lmMst2) + 64*lmMst1*pow2(1 + lmMst2) -
        32*((1 + lmMst2)*pow2(lmMst1) + pow3(lmMst2)))) + (1 + lmMst2)*(-96*(2*
        lmMst2*(2 + lmMst2) - lmMst1*(3 + 2*lmMst2))*Mt*pow2(s2t)*pow3(Mst2) +
        64*Mst2*(1 + lmMst2*(7 - 2*lmMt) - 2*lmMst1*(2 + lmMst2 - lmMt) - 3*
        lmMt + 2*pow2(lmMst2))*pow3(Mt)) + (8 + 3*lmMst2*(-79 + 45*shiftst1) -
        16*(1 + lmMst2)*pow2(lmMst1) - 308*pow2(lmMst2) + lmMst1*(269 + 348*
        lmMst2 - 135*shiftst1 + 123*pow2(lmMst2)) - 107*pow3(lmMst2))*pow3(s2t)
        *pow4(Mst2)) - pow2(Mst1)*(60*Dmsqst2*s2t*(5 - 9*shiftst1)*pow2(Mt) -
        64*pow2(1 + lmMst2)*(3*Mt*pow2(s2t)*pow3(Mst2) - 4*Mst2*pow3(Mt)) +
        pow2(Mst2)*(4*s2t*(173 + 188*lmMst2 + 32*lmMst1*(1 + lmMst2) - 135*
        shiftst1 + 59*pow2(lmMst2))*pow2(Mt) + 15*Dmsqst2*(-5 + 9*shiftst1)*
        pow3(s2t)) - (221 + 284*lmMst2 + 32*lmMst1*(1 + lmMst2) - 135*shiftst1
        + 107*pow2(lmMst2))*pow3(s2t)*pow4(Mst2))*pow5(Mst2) + pow2(Mst2)*(64*(
        1 + lmMst2)*Mt*(-2*Mst2*Mt*s2t*pow2(1 - lmMst1 + lmMst2) + 2*(3 - 21*
        lmMst2 + 2*lmMst1*(8 + 3*lmMst2 - 3*lmMt) + 5*lmMt + 6*lmMst2*lmMt - 6*
        pow2(lmMst2))*pow2(Mt) + 3*(-2 + 5*lmMst2 - lmMst1*(5 + 2*lmMst2) + 2*
        pow2(lmMst2))*pow2(Mst2)*pow2(s2t)) - Mst2*(15*Dmsqst2*(5 - 9*shiftst1)
        + pow2(Mst2)*(125 + 156*lmMst2 + 512*lmMst1*lmMst2*(1 + lmMst2) - 135*
        shiftst1 - 181*pow2(lmMst2) - 256*((1 + lmMst2)*pow2(lmMst1) + pow3(
        lmMst2))))*pow3(s2t))*pow6(Mst1) + 16*(1 + lmMst2)*((-8*Mst2*s2t*(1 +
        3*lmMst2 - lmMst1*(3 + 2*lmMst2) + pow2(lmMst1) + pow2(lmMst2))*pow2(
        Mt) - 6*Mt*(5 - 12*lmMst2 + 4*lmMst1*(3 + lmMst2) - 4*pow2(lmMst2))*
        pow2(Mst2)*pow2(s2t) + 8*(12 - 45*lmMst2 + 2*lmMst1*(19 + 6*lmMst2 - 6*
        lmMt) + 7*lmMt + 12*lmMst2*lmMt - 12*pow2(lmMst2))*pow3(Mt) + (1 +
        lmMst1*(2 - 32*lmMst2) - 2*lmMst2 + 16*pow2(lmMst1) + 16*pow2(lmMst2))*
        pow3(Mst2)*pow3(s2t))*pow8(Mst1) + (1 + lmMst2)*s2t*(4*pow2(Mt) - pow2(
        Mst2)*pow2(s2t))*pow9(Mst2)))) + Tbeta*(pow2(Mst1)*(-9*pow6(Mst2)*(512*
        Dmglst2*s2t*(-1 + 2*lmMst2 + 3*pow2(lmMst2))*pow2(MuSUSY)*pow3(Mt) + 2*
        Dmglst2*(-20 + (262 + 32*lmMst1)*lmMst2 + 187*pow2(lmMst2))*pow2(Sbeta)
        *pow4(s2t)*pow5(Mst2) + (237 + 316*lmMst2 + 32*lmMst1*(1 + lmMst2) +
        123*pow2(lmMst2))*pow2(Sbeta)*pow4(s2t)*pow6(Mst2)) + 1215*shiftst1*(
        pow2(Sbeta)*pow4(Mst2)*(16*(Dmsqst2*pow2(Mst1)*pow2(Mst2)*pow2(Mt)*
        pow2(s2t) + (Dmsqst2 + pow2(Mst2))*(pow2(Mst1) + pow2(Mst2))*pow4(Mt))
        + (pow2(Mst1) - pow2(Mst2))*(Dmsqst2 + pow2(Mst2))*(2*(lmMst1 - lmMst2)
        *pow2(Mst1)*pow2(Mst2) + pow4(Mst1) - pow4(Mst2))*pow4(s2t)) + pow2(Mt)
        *pow2(s2t)*(4*(Dmsqst2 + pow2(Mst2))*pow2(MuSUSY)*((pow2(Mst1) + pow2(
        Mst2))*pow4(Mst2) - 2*(lmMst1 - lmMst2)*pow2(Mst1)*(pow2(Mst1)*pow2(
        Mst2) + pow4(Mst1) + pow4(Mst2))) - 4*pow2(Sbeta)*((Dmsqst2 + pow2(
        Mst2))*pow2(MuSUSY)*((pow2(Mst1) + pow2(Mst2))*pow4(Mst2) - 2*(lmMst1 -
        lmMst2)*pow2(Mst1)*(pow2(Mst1)*pow2(Mst2) + pow4(Mst1) + pow4(Mst2))) +
        2*(Dmsqst2*pow4(Mst2)*(pow4(Mst1) + pow4(Mst2)) + pow2(pow2(Mst1) -
        pow2(Mst2))*pow6(Mst2)))))) + pow2(Sbeta)*(256*s2t*pow2(Mst1)*pow3(Mt)*
        (9*(1 + lmMst2)*Mst2*(-2*pow2(Mst2)*((3 + 2*lmMst1*(7 + 2*lmMst2 - 2*
        lmMt) + 4*lmMst2*(-4 + lmMt) + 2*lmMt - 4*pow2(lmMst2))*pow2(Mst2) + (1
        - 10*lmMst2 + 4*lmMst1*(2 + lmMst2) - 4*pow2(lmMst2))*pow2(MuSUSY))*
        pow4(Mst1) + pow2(Mst1)*((1 - 2*lmMst1*(5 + 2*lmMst2 - 2*lmMt) - 4*
        lmMst2*(-3 + lmMt) - 6*lmMt + 4*pow2(lmMst2))*pow2(Mst2) + 2*(1 + 5*
        lmMst2 - lmMst1*(3 + 2*lmMst2) + 2*pow2(lmMst2))*pow2(MuSUSY))*pow4(
        Mst2) - (2*(8 - 27*lmMst2 + lmMst1*(25 + 6*lmMst2 - 6*lmMt) + 2*lmMt +
        6*lmMst2*lmMt - 6*pow2(lmMst2))*pow2(Mst2) + (7 - 32*lmMst2 + 4*lmMst1*
        (7 + 3*lmMst2) - 12*pow2(lmMst2))*pow2(MuSUSY))*pow6(Mst1) + 2*(1 +
        lmMst2)*(2*pow2(Mst2) + pow2(MuSUSY))*pow6(Mst2)) - Dmglst2*(pow2(Mst1)
        *(pow2(Mst2)*(253 + 5*lmMst2*(107 - 6*lmMt) - 102*lmMt - 18*lmMst1*(9 +
        lmMst2 - 2*lmMt + 2*lmMst2*lmMt - 2*pow2(lmMst2)) + 12*(10 + 3*lmMt)*
        pow2(lmMst2) - 36*pow3(lmMst2)) + 18*pow2(MuSUSY)*(8 + 7*lmMst2 - 11*
        pow2(lmMst2) + lmMst1*(-3 + 5*lmMst2 + 6*pow2(lmMst2)) - 6*pow3(lmMst2)
        ))*pow4(Mst2) + 9*((49 + 103*lmMst2 - 36*(1 + lmMst2)*pow2(lmMst2) + 4*
        lmMst1*(-11 + 6*lmMst2 + 9*pow2(lmMst2)))*pow2(MuSUSY) + pow2(Mst2)*(28
        + lmMst2*(378 - 96*lmMt) - 44*lmMt + lmMst1*(-274 + 60*lmMt + 18*
        lmMst2*(-13 + 2*lmMt) - 36*pow2(lmMst2)) + (270 - 36*lmMt)*pow2(lmMst2)
        + 36*pow3(lmMst2)))*pow6(Mst1) + 18*(pow2(Mst2)*(2*pow2(MuSUSY)*(8 +
        13*lmMst2 - 8*pow2(lmMst2) + lmMst1*(-5 + 5*lmMst2 + 6*pow2(lmMst2)) -
        6*pow3(lmMst2)) + pow2(Mst2)*(23 + lmMst2*(93 - 22*lmMt) - 14*lmMt - 4*
        (-11 + lmMt)*pow2(lmMst2) - 2*lmMst1*(25 + lmMst2*(17 - 2*lmMt) - 6*
        lmMt + 2*pow2(lmMst2)) + 4*pow3(lmMst2)))*pow4(Mst1) + (1 - 2*lmMst2 -
        3*pow2(lmMst2))*(2*pow2(Mst2) + pow2(MuSUSY))*pow6(Mst2)))) + 9*pow2(
        Mst1)*pow3(Mst2)*pow3(s2t)*(128*Mst2*Mt*((1 + lmMst2)*Mst2*(2*(2 + 2*
        lmMst1 - lmMst2)*pow2(Mst2)*pow4(Mst1) + 2*(1 - 3*lmMst2 + lmMst1*(3 +
        2*lmMst2) - 2*pow2(lmMst2))*pow2(Mst1)*pow4(Mst2) + (1 + 2*lmMst1 - 2*
        lmMst2)*pow6(Mst1) - 2*(1 + lmMst2)*pow6(Mst2)) + Dmglst2*(2*(1 - 4*
        lmMst1 + 10*lmMst2 + 3*pow2(lmMst2))*pow2(Mst2)*pow4(Mst1) + 2*pow2(
        Mst1)*(6 + 11*lmMst2 - 5*pow2(lmMst2) + lmMst1*(-3 + 5*lmMst2 + 6*pow2(
        lmMst2)) - 6*pow3(lmMst2))*pow4(Mst2) + (1 + 13*lmMst2 - 2*lmMst1*(5 +
        3*lmMst2) + 6*pow2(lmMst2))*pow6(Mst1) - 2*(-1 + 2*lmMst2 + 3*pow2(
        lmMst2))*pow6(Mst2))) + s2t*(30*Dmglst2*Dmsqst2*(pow2(Mst1) - pow2(
        Mst2))*(2*(lmMst1 - lmMst2)*pow2(Mst1)*pow2(Mst2) + pow4(Mst1) - pow4(
        Mst2)) - 2*Dmglst2*pow2(Mst1)*pow2(Mst2)*(pow2(Mst1)*pow2(Mst2)*(-44 +
        34*lmMst2 + 224*lmMst2*pow2(lmMst1) - 359*pow2(lmMst2) - 2*lmMst1*(52 -
        198*lmMst2 + 133*pow2(lmMst2)) + 42*pow3(lmMst2)) + (-36 + (70 + 32*
        lmMst1)*lmMst2 + 27*pow2(lmMst2))*pow4(Mst1) + (100 - 222*lmMst2 + 32*
        lmMst2*pow2(lmMst1) + lmMst1*(8 - 524*lmMst2 - 246*pow2(lmMst2)) + 241*
        pow2(lmMst2) + 214*pow3(lmMst2))*pow4(Mst2)) - Mst2*(75*Dmsqst2*(pow2(
        Mst1) - pow2(Mst2))*(2*(lmMst1 - lmMst2)*pow2(Mst1)*pow2(Mst2) + pow4(
        Mst1) - pow4(Mst2)) + (-109 - 630*lmMst2 + 224*(1 + lmMst2)*pow2(
        lmMst1) + lmMst1*(538 + 184*lmMst2 - 266*pow2(lmMst2)) - 435*pow2(
        lmMst2) + 42*pow3(lmMst2))*pow4(Mst1)*pow4(Mst2) + (141 + 140*lmMst2 +
        32*lmMst1*(1 + lmMst2) + 43*pow2(lmMst2))*pow2(Mst2)*pow6(Mst1) + pow2(
        Mst1)*(-237 + 190*lmMst2 + 32*(1 + lmMst2)*pow2(lmMst1) + 509*pow2(
        lmMst2) - 2*lmMst1*(285 + 364*lmMst2 + 123*pow2(lmMst2)) + 214*pow3(
        lmMst2))*pow6(Mst2)))) + 16*Mst2*pow4(Mt)*(Mst2*(-(pow4(Mst1)*(675*
        Dmsqst2*pow2(Mst2) + (2117 + lmMst2*(4796 - 1248*lmMt) - 576*lmMst1*(1
        + lmMst2)*(3 + lmMst2 - lmMt) - 672*lmMt + (3651 - 576*lmMt)*pow2(
        lmMst2) + 576*pow3(lmMst2))*pow4(Mst2))) + 288*(1 + lmMst2)*((5 - 57*
        lmMst2 + 2*lmMst1*(25 + 6*lmMst2 - 6*lmMt) + 7*lmMt + 12*lmMst2*lmMt -
        12*pow2(lmMst2))*pow2(Mst1) + (-2 - 27*lmMst2 + lmMst1*(22 + 6*lmMst2 -
        6*lmMt) + 5*lmMt + 6*lmMst2*lmMt - 6*pow2(lmMst2))*pow2(Mst2))*pow6(
        Mst1) - 9*pow2(Mst1)*(75*Dmsqst2*pow4(Mst2) + (173 + 188*lmMst2 + 32*
        lmMst1*(1 + lmMst2) + 59*pow2(lmMst2))*pow6(Mst2)) + 144*pow2(1 +
        lmMst2)*pow8(Mst2)) + 2*Dmglst2*(pow4(Mst1)*(135*Dmsqst2*pow2(Mst2) + (
        4220 + 6114*lmMst2 - 576*lmMst1*(4 + 2*lmMst2 - lmMt) - 48*(29 + 27*
        lmMst2)*lmMt + 1341*pow2(lmMst2))*pow4(Mst2)) + 9*pow2(Mst1)*(15*
        Dmsqst2*pow4(Mst2) - (-84 + (70 + 32*lmMst1)*lmMst2 + 59*pow2(lmMst2))*
        pow6(Mst2)) + 288*(pow2(Mst2)*(31 + lmMst2*(95 - 23*lmMt) - 16*lmMt + (
        51 - 6*lmMt)*pow2(lmMst2) - 2*lmMst1*(25 + lmMst2*(20 - 3*lmMt) - 6*
        lmMt + 3*pow2(lmMst2)) + 6*pow3(lmMst2))*pow6(Mst1) + (42 + lmMst2*(242
        - 62*lmMt) - 33*lmMt + 6*(29 - 4*lmMt)*pow2(lmMst2) - 2*lmMst1*(81 +
        lmMst2*(74 - 12*lmMt) - 18*lmMt + 12*pow2(lmMst2)) + 24*pow3(lmMst2))*
        pow8(Mst1) + lmMst2*(1 + lmMst2)*pow8(Mst2)))) - 12*pow2(Mt)*pow2(s2t)*
        (225*Dmsqst2*pow2(Mst1)*(-2*pow4(Mst1)*((-lmMst1 + lmMst2)*pow2(Mst2)*
        pow2(MuSUSY) + pow4(Mst2)) + 2*(lmMst1 - lmMst2)*pow2(MuSUSY)*pow6(
        Mst1) - (2*pow2(Mst2) + pow2(MuSUSY))*pow6(Mst2) + pow2(Mst1)*((-1 + 2*
        lmMst1 - 2*lmMst2)*pow2(MuSUSY)*pow4(Mst2) + 4*pow6(Mst2))) + 2*
        Dmglst2*Mst2*(-((3*pow2(MuSUSY)*(60 + 206*lmMst2 + 32*lmMst2*pow2(
        lmMst1) + lmMst1*(8 - 460*lmMst2 - 246*pow2(lmMst2)) + 519*pow2(lmMst2)
        + 214*pow3(lmMst2)) + 4*pow2(Mst2)*(52 + lmMst2*(-338 + 48*pow2(lmMst1)
        ) - 129*pow2(lmMst2) + 48*(lmMst1 - 2*lmMst1*lmMst2*(1 + lmMst2) +
        pow3(lmMst2))))*pow4(Mst1)*pow4(Mst2)) + 2*pow2(Mst2)*((332 + 302*
        lmMst2 - 288*lmMst1*(1 + lmMst2) + 111*pow2(lmMst2))*pow2(Mst2) - 3*
        pow2(MuSUSY)*(48 + 4*lmMst2*(31 + 36*pow2(lmMst1)) + 278*pow2(lmMst2) -
        lmMst1*(44 + 278*lmMst2 + 379*pow2(lmMst2)) + 235*pow3(lmMst2)))*pow6(
        Mst1) + 45*Dmsqst2*(pow4(Mst1)*((1 - 2*lmMst1 + 2*lmMst2)*pow2(Mst2)*
        pow2(MuSUSY) - 4*pow4(Mst2)) + pow2(Mst1)*(2*pow2(Mst2) + pow2(MuSUSY))
        *pow4(Mst2) + 2*(pow2(Mst2) + (-lmMst1 + lmMst2)*pow2(MuSUSY))*pow6(
        Mst1)) + 3*pow2(Mst1)*(-2*(-52 + 2*(67 + 16*lmMst1)*lmMst2 + 91*pow2(
        lmMst2))*pow2(Mst2) - (-20 + 2*(99 + 16*lmMst1)*lmMst2 + 123*pow2(
        lmMst2))*pow2(MuSUSY))*pow6(Mst2) + 6*(32*(2 + 7*lmMst2 - lmMst1*(5 +
        4*lmMst2) + 4*pow2(lmMst2))*pow2(Mst2) - pow2(MuSUSY)*(48 + 4*lmMst2*(
        45 + 68*pow2(lmMst1)) + 310*pow2(lmMst2) - lmMst1*(92 + 310*lmMst2 +
        635*pow2(lmMst2)) + 363*pow3(lmMst2)))*pow8(Mst1) + 96*lmMst2*(1 +
        lmMst2)*(2*pow2(Mst2) + pow2(MuSUSY))*pow8(Mst2)) - 3*pow2(Mst2)*(282*
        pow4(Mst2)*pow6(Mst1) - 564*pow4(Mst1)*pow6(Mst2) + pow4(Mst1)*(32*(1 +
        lmMst2)*pow2(lmMst1)*(pow2(MuSUSY)*(9*pow2(Mst1)*pow2(Mst2) + 17*pow4(
        Mst1) + pow4(Mst2)) + 2*pow6(Mst2)) + 2*pow3(lmMst2)*(pow2(MuSUSY)*(
        235*pow2(Mst1)*pow2(Mst2) + 363*pow4(Mst1) + 107*pow4(Mst2)) + 32*pow6(
        Mst2))) - 2*lmMst1*((253 + 588*lmMst2 + 379*pow2(lmMst2))*pow2(Mst2)*
        pow2(MuSUSY)*pow6(Mst1) - 16*(1 + lmMst2)*pow2(Mst1)*(2*pow2(Mst2) +
        pow2(MuSUSY))*pow6(Mst2) + pow4(Mst1)*((253 + 332*lmMst2 + 123*pow2(
        lmMst2))*pow2(MuSUSY)*pow4(Mst2) + 32*(3 + 5*lmMst2 + 2*pow2(lmMst2))*
        pow6(Mst2)) + (32*(1 + lmMst2)*pow2(Mst2) + (237 + 828*lmMst2 + 635*
        pow2(lmMst2))*pow2(MuSUSY))*pow8(Mst1)) + pow2(MuSUSY)*(189*pow4(Mst1)*
        pow4(Mst2) + 64*pow2(Mst2)*pow6(Mst1) + 205*pow2(Mst1)*pow6(Mst2) + 80*
        pow8(Mst1) - 16*pow8(Mst2)) + 378*pow2(Mst1)*pow8(Mst2) + 2*lmMst2*((
        285*pow2(Mst2)*pow2(MuSUSY) + 172*pow4(Mst2))*pow6(Mst1) - 33*pow4(
        Mst1)*(-11*pow2(MuSUSY)*pow4(Mst2) + 8*pow6(Mst2)) + (32*pow2(Mst2) +
        277*pow2(MuSUSY))*pow8(Mst1) - 16*(2*pow2(Mst2) + pow2(MuSUSY))*pow8(
        Mst2) + 2*pow2(Mst1)*(63*pow2(MuSUSY)*pow6(Mst2) + 110*pow8(Mst2))) -
        pow2(lmMst2)*(-6*(148*pow2(Mst2)*pow2(MuSUSY) + 25*pow4(Mst2))*pow6(
        Mst1) + pow4(Mst1)*(-707*pow2(MuSUSY)*pow4(Mst2) + 76*pow6(Mst2)) - 8*(
        8*pow2(Mst2) + 139*pow2(MuSUSY))*pow8(Mst1) + 16*(2*pow2(Mst2) + pow2(
        MuSUSY))*pow8(Mst2) - pow2(Mst1)*(91*pow2(MuSUSY)*pow6(Mst2) + 150*
        pow8(Mst2))) - 32*power10(Mst2)))) + 36*s2t*pow2(Mt)*pow2(MuSUSY)*(-((
        75*Dmsqst2*(1 - 2*lmMst1 + 2*lmMst2)*s2t + 128*(1 + lmMst2)*Mst2*Mt*(1
        + 5*lmMst2 - lmMst1*(3 + 2*lmMst2) + 2*pow2(lmMst2)) + s2t*pow2(Mst2)*(
        189 + 726*lmMst2 + 32*(1 + lmMst2)*pow2(lmMst1) + 707*pow2(lmMst2) - 2*
        lmMst1*(253 + 332*lmMst2 + 123*pow2(lmMst2)) + 214*pow3(lmMst2)))*pow4(
        Mst1)*pow4(Mst2)) - pow2(Mst1)*(75*Dmsqst2*s2t + 128*Mst2*Mt*pow2(1 +
        lmMst2) + s2t*(205 + 252*lmMst2 + 32*lmMst1*(1 + lmMst2) + 91*pow2(
        lmMst2))*pow2(Mst2))*pow6(Mst2) + 2*(pow2(Mst2)*(75*Dmsqst2*(lmMst1 -
        lmMst2)*s2t + 64*(1 + lmMst2)*Mst2*Mt*(1 - 10*lmMst2 + 4*lmMst1*(2 +
        lmMst2) - 4*pow2(lmMst2)) - s2t*pow2(Mst2)*(32 + 285*lmMst2 + 144*(1 +
        lmMst2)*pow2(lmMst1) + 444*pow2(lmMst2) - lmMst1*(253 + 588*lmMst2 +
        379*pow2(lmMst2)) + 235*pow3(lmMst2)))*pow6(Mst1) + (75*Dmsqst2*(lmMst1
        - lmMst2)*s2t + 32*(1 + lmMst2)*Mst2*Mt*(7 - 32*lmMst2 + 4*lmMst1*(7 +
        3*lmMst2) - 12*pow2(lmMst2)) - s2t*pow2(Mst2)*(40 + 277*lmMst2 + 272*(1
        + lmMst2)*pow2(lmMst1) + 556*pow2(lmMst2) - lmMst1*(237 + 828*lmMst2 +
        635*pow2(lmMst2)) + 363*pow3(lmMst2)))*pow8(Mst1) + Dmglst2*((15*
        Dmsqst2*(1 - 2*lmMst1 + 2*lmMst2)*s2t + 64*Mst2*Mt*(8 + 7*lmMst2 - 11*
        pow2(lmMst2) + lmMst1*(-3 + 5*lmMst2 + 6*pow2(lmMst2)) - 6*pow3(lmMst2)
        ) - s2t*pow2(Mst2)*(60 + 206*lmMst2 + 32*lmMst2*pow2(lmMst1) + lmMst1*(
        8 - 460*lmMst2 - 246*pow2(lmMst2)) + 519*pow2(lmMst2) + 214*pow3(
        lmMst2)))*pow3(Mst2)*pow4(Mst1) + 2*Mst2*(15*Dmsqst2*(-lmMst1 + lmMst2)
        *s2t + 64*Mst2*Mt*(8 + 13*lmMst2 - 8*pow2(lmMst2) + lmMst1*(-5 + 5*
        lmMst2 + 6*pow2(lmMst2)) - 6*pow3(lmMst2)) - s2t*pow2(Mst2)*(48 + 124*
        lmMst2 + 144*lmMst2*pow2(lmMst1) + 278*pow2(lmMst2) - lmMst1*(44 + 278*
        lmMst2 + 379*pow2(lmMst2)) + 235*pow3(lmMst2)))*pow6(Mst1) + 2*(16*Mt*(
        49 + 103*lmMst2 - 36*(1 + lmMst2)*pow2(lmMst2) + 4*lmMst1*(-11 + 6*
        lmMst2 + 9*pow2(lmMst2))) - Mst2*s2t*(48 + 4*lmMst2*(45 + 68*pow2(
        lmMst1)) + 310*pow2(lmMst2) - lmMst1*(92 + 310*lmMst2 + 635*pow2(
        lmMst2)) + 363*pow3(lmMst2)))*pow8(Mst1) + s2t*(pow2(Mst1)*(15*Dmsqst2
        - (-20 + 2*(99 + 16*lmMst1)*lmMst2 + 123*pow2(lmMst2))*pow2(Mst2))*
        pow5(Mst2) + 32*lmMst2*(1 + lmMst2)*pow9(Mst2)))) + 16*s2t*pow2(1 +
        lmMst2)*power10(Mst2))))))/(11664.*Tbeta*pow2(Sbeta)*pow4(Mst1)*pow6(
        Mst2)) + threeLoopFlag*(pow4(s2t)*(Dmsqst2*(3.2221604938271606 +
        lmMst1*(5.188888888888889 - 7*lmMst2) - (317*lmMst2)/90. + (11*pow2(
        lmMst1))/6. + (31*pow2(lmMst2))/6.)*pow2(Mst1) - pow2(Mst2)*((5*
        Dmsqst2*(1213 + 48*lmMst1*(7 - 3*lmMst2) - 408*lmMst2 + 144*pow2(
        lmMst2)))/216. + pow2(Mst1)*(79.01365226337448 - B4/9. + (2*D3)/9. -
        DN/6. + (4372*lmMst1)/405. + (lmMst2*(16051 + 22980*lmMst1 - 5175*pow2(
        lmMst1)))/810. - (1631*pow2(lmMst1))/108. + ((-92 + 327*lmMst1)*pow2(
        lmMst2))/27. - (79*pow3(lmMst1))/54. - (115*pow3(lmMst2))/27.)) +
        Dmglst2*(24.239197530864196 - (2*B4)/3. + (8*D3)/9. - (5*DN)/9. + (22*
        lmMst1)/3. - (13*pow2(lmMst1))/6. + (lmMst2*(-277 - 264*lmMst1 + 16*
        pow2(lmMst1)))/9. + ((142 - 123*lmMst1)*pow2(lmMst2))/9. + (10*Dmsqst2)
        /pow2(Mst1) + (107*pow3(lmMst2))/9.)*pow3(Mst2) + (46.745471895783595 +
        B4 + D3/9. - DN/9. + (34838747*lmMst1)/529200. + (50723*pow2(lmMst1))/
        1890. + (lmMst2*(-23468297 - 17276980*lmMst1 + 3601500*pow2(lmMst1)))/
        529200. + (4*(2941 - 2415*lmMst1)*pow2(lmMst2))/945. + (Dmsqst2*(
        15.13649301931237 + lmMst1*(13.99649785840262 - (106*lmMst2)/21.) - (
        621671*lmMst2)/39690. + (53*pow2(lmMst1))/21. + (53*pow2(lmMst2))/21.))
        /pow2(Mst2) - (10*pow3(lmMst1))/27. + (409*pow3(lmMst2))/108. + (5*
        Dmglst2*Dmsqst2*(76 - 249*lmMst1 + 249*lmMst2))/(54.*pow3(Mst2)))*pow4(
        Mst1) + (Dmglst2*((5*Dmsqst2*(-75 + 26*lmMst1 - 26*lmMst2)*pow2(Mst1))/
        6. + pow2(Mst2)*((20*Dmsqst2*(1 + lmMst1 - lmMst2))/3. - pow2(Mst1)*(
        0.7822304526748971 - (2*B4)/3. + (8*D3)/9. - (5*DN)/9. + (81193*lmMst1)
        /4050. + (137*pow2(lmMst1))/135. - (lmMst2*(81643 + 86970*lmMst1 +
        48150*pow2(lmMst1)))/4050. + ((4969 + 3840*lmMst1)*pow2(lmMst2))/270. -
        (5*pow3(lmMst1))/27. - (58*pow3(lmMst2))/27.)) + (79.58863384550371 + (
        8287903*lmMst2)/1.1907e6 + (4.326984126984127 + (11*lmMst2)/3.)*pow2(
        lmMst1) + lmMst1*(1.2061367262954565 - (101*lmMst2)/315. - (11*pow2(
        lmMst2))/3.) - (51*pow2(lmMst2))/70. - (11*pow3(lmMst1))/9. + (11*pow3(
        lmMst2))/9.)*pow4(Mst1)))/Mst2 - ((25289 + 1440*B4 - 144*D3 + 72*DN +
        22368*lmMst1 + 1908*pow2(lmMst1) - 12*lmMst2*(2296 - 2136*lmMst1 + 117*
        pow2(lmMst1)) + 504*(-53 + 21*lmMst1)*pow2(lmMst2) + (1080*Dmsqst2*(-1
        + 2*lmMst2))/pow2(Mst1) + 576*pow3(lmMst1) - 9756*pow3(lmMst2))*pow4(
        Mst2))/1296. + ((Mst2*(135 + 250*lmMst2 + 32*lmMst1*(1 + lmMst2) + 123*
        pow2(lmMst2)) + Dmglst2*(-40 + (524 + 64*lmMst1)*lmMst2 + 374*pow2(
        lmMst2)))*pow5(Mst2))/(36.*pow2(Mst1)) + (30*Dmsqst2*(9*(56*OepS2 - 27*
        (253 + 42*lmMst1 - 42*lmMst2)*S2)*pow2(Mst1)*pow2(Mst2) + 2*(56*OepS2 -
        27*(163 + 42*lmMst1 - 42*lmMst2)*S2)*pow4(Mst1) + 27*(8*OepS2 - 81*(-15
        + 2*lmMst1 - 2*lmMst2)*S2)*pow4(Mst2)) + Mst2*(9*(5144*OepS2 - 81*(5717
        + 1286*lmMst1 - 1286*lmMst2)*S2)*pow2(Mst1)*pow3(Mst2) + 15*Mst2*(1436*
        OepS2 - 27*(3370 + 1077*lmMst1 - 1077*lmMst2)*S2)*pow4(Mst1) - 2*
        Dmglst2*(972*(16*OepS2 + (275 - 324*lmMst1 + 324*lmMst2)*S2)*pow2(Mst1)
        *pow2(Mst2) + (33592*OepS2 - 27*(51635 + 25194*lmMst1 - 25194*lmMst2)*
        S2)*pow4(Mst1) + 27*(56*OepS2 - 81*(-1677 + 14*lmMst1 - 14*lmMst2)*S2)*
        pow4(Mst2)) + 81*(184*OepS2 - 81*(307 + 46*lmMst1 - 46*lmMst2)*S2)*
        pow5(Mst2)))/(8748.*pow2(Mst2))) - pow2(Mt)*pow2(s2t)*((-5*Dmglst2*
        Dmsqst2*(-1081 + 165*lmMst1 - 165*lmMst2)*(-2 + pow2(Sbeta))*pow2(
        Sbeta)*pow4(Mst1))/(9.*pow5(Mst2)) + (pow2(MuSUSY)*(53.385802469135804
         + (40*B4)/
        9. - (4*D3)/9. + (2*DN)/9. + (1672*lmMst1)/27. + (53*pow2(lmMst1))/9. -
        lmMst2*(129.92592592592592 - 72*lmMst1 + (13*pow2(lmMst1))/3.) + (2*(-
        470 + 147*lmMst1)*pow2(lmMst2))/9. - (2*Dmglst2*Mst2*(-20 + 2*(99 + 16*
        lmMst1)*lmMst2 + 123*pow2(lmMst2)))/(9.*pow2(Mst1)) - ((103 + 186*
        lmMst2 + 32*lmMst1*(1 + lmMst2) + 91*pow2(lmMst2))*pow2(Mst2))/(9.*
        pow2(Mst1)) + (16*pow3(lmMst1))/9. + (2250*Dmsqst2*(1141 + 48*lmMst1*(7
        - 3*lmMst2) - 264*lmMst2 + 144*pow2(lmMst2)) + pow2(Mst1)*(10552777 +
        205200*B4 - 5400*DN + 4145280*lmMst1 - 1181700*pow2(lmMst1) - 240*
        lmMst2*(16192 - 26430*lmMst1 + 3465*pow2(lmMst1)) + 900*(-5735 + 3072*
        lmMst1)*pow2(lmMst2) - 55800*pow3(lmMst1) - 1877400*pow3(lmMst2)))/(
        24300.*pow2(Mst2)) - (271*pow3(lmMst2))/9. - (Dmglst2*(14267 - 432*B4 +
        576*D3 - 360*DN + 4752*lmMst1 - 1404*pow2(lmMst1) + 72*lmMst2*(-63 -
        232*lmMst1 + 16*pow2(lmMst1)) + 72*(281 - 123*lmMst1)*pow2(lmMst2) + (
        6480*Dmsqst2)/pow2(Mst1) + 7704*pow3(lmMst2)))/(162.*Mst2) - (2*
        Dmglst2*(405000*Dmsqst2*(4 + lmMst1 - lmMst2) + pow2(Mst1)*(2695042 -
        40500*B4 + 54000*D3 - 33750*DN - 326895*lmMst1 - 324900*pow2(lmMst1) +
        15*lmMst2*(-19607 - 129030*lmMst1 + 62550*pow2(lmMst1)) - 450*(-5023 +
        5610*lmMst1)*pow2(lmMst2) + 11250*pow3(lmMst1) + 1575000*pow3(lmMst2)))
        )/(30375.*pow3(Mst2)) + (16*(1 + lmMst2)*(4*Dmglst2*lmMst2 + Mst2 +
        lmMst2*Mst2)*pow3(Mst2))/(9.*pow4(Mst1)) + (Dmsqst2*(201.74098765432097
         - (622*lmMst2)/15. + (2*lmMst1*(311 + 10*lmMst2))/15. - (22*pow2(
        lmMst1))/3. + 6*pow2(lmMst2))*pow2(Mst1) + (628.1736268201578 + (76*B4)
        /9. - (2*DN)/9. + (6317839*lmMst1)/396900. - (66307*pow2(lmMst1))/315.
         - lmMst2*(
        12.52907281431091 - (182909*lmMst1)/315. + (274*pow2(lmMst1))/3.) + (2*
        (-58301 + 37135*lmMst1)*pow2(lmMst2))/315. - (44*pow3(lmMst1))/9. - (
        1256*pow3(lmMst2))/9.)*pow4(Mst1))/pow4(Mst2) - (Dmglst2*(
        585.1892843532082 - (8*B4)/3. + (32*D3)/9. - (20*DN)/9. - (20109937*
        lmMst1)/297675. - (15886*pow2(lmMst1))/945. + lmMst2*(
        17.112243218274966 - (144628*lmMst1)/945. + (1180*pow2(lmMst1))/9.) + (
        2*(80257 - 138180*lmMst1)*pow2(lmMst2))/945. - (92*pow3(lmMst1))/27. +
        (4448*pow3(lmMst2))/27.)*pow4(Mst1))/pow5(Mst2) - (Dmsqst2*(10 - 20*
        lmMst2 + (10*Dmglst2*(-23 + 42*lmMst1 - 42*lmMst2)*pow4(Mst1))/pow5(
        Mst2) - (3*(237.28785508324435 + (16526*lmMst2)/3969. + (2*lmMst1*(-
        8263 + 71820*lmMst2))/3969. - (520*pow2(lmMst1))/21. - (80*pow2(lmMst2)
        )/7.)*pow6(Mst1))/pow6(Mst2)))/(3.*pow2(Mst1)) + (-30*Dmsqst2*(9*(104*
        OepS2 + 27*(17 - 78*lmMst1 + 78*lmMst2)*S2)*pow2(Mst1)*pow2(Mst2) + 13*
        (136*OepS2 - 27*(95 + 102*lmMst1 - 102*lmMst2)*S2)*pow4(Mst1) + 27*(8*
        OepS2 - 81*(-15 + 2*lmMst1 - 2*lmMst2)*S2)*pow4(Mst2)) + Mst2*(-9*(
        8456*OepS2 - 81*(11243 + 2114*lmMst1 - 2114*lmMst2)*S2)*pow2(Mst1)*
        pow3(Mst2) + 3*Mst2*(-52948*OepS2 + 27*(194357 + 39711*lmMst1 - 39711*
        lmMst2)*S2)*pow4(Mst1) + 2*Dmglst2*(54*(344*OepS2 + 9*(15643 - 774*
        lmMst1 + 774*lmMst2)*S2)*pow2(Mst1)*pow2(Mst2) + 4*(17308*OepS2 + 27*(
        93919 - 12981*lmMst1 + 12981*lmMst2)*S2)*pow4(Mst1) + 27*(56*OepS2 -
        81*(-1677 + 14*lmMst1 - 14*lmMst2)*S2)*pow4(Mst2)) + 81*(-184*OepS2 +
        81*(307 + 46*lmMst1 - 46*lmMst2)*S2)*pow5(Mst2)))/(2187.*pow6(Mst2))))/
        pow2(Sbeta)) - (pow2(MuSUSY)*pow3(Mt)*(-27*pow2(Mst1)*pow2(Mst2)*(15*
        Mst2*(32*Mst2*Mt*(436 - 6*B4 + 6*D3 - 3*DN - 48*lmMst1 + (408 - 96*
        lmMst1)*lmMst2 + 24*lmMt - 972*S2 - 48*(-4 + lmMst1)*pow2(lmMst2) + 48*
        pow3(lmMst2)) + s2t*(5760*Dmsqst2*(1 + lmMst1 - lmMst2) + pow2(Mst2)*(
        28683 + 5280*B4 - 48*DN - 5952*lmMst1 - 224*OepS2 + 324*(-1 + 14*lmMst1
        - 14*lmMst2)*S2 - 768*pow2(lmMst1) - 192*lmMst2*(-214 + 73*lmMst1 + 4*
        pow2(lmMst1)) - 96*(-268 + 57*lmMst1)*pow2(lmMst2) + 6240*pow3(lmMst2))
        )) + Dmglst2*(2880*Mst2*Mt*(180 - 2*B4 + 2*D3 - DN + 16*lmMst1 + 144*
        lmMst2 - 216*S2 - 16*(-2 + lmMst1)*pow2(lmMst2) + 16*pow3(lmMst2)) +
        s2t*(-86400*Dmsqst2*(2 + 5*lmMst1 - 5*lmMst2) + pow2(Mst2)*(23917 +
        188640*B4 - 3600*DN - 37440*lmMst1 + 5600*OepS2 - 324*(-453 + 350*
        lmMst1 - 350*lmMst2)*S2 + 11520*pow2(lmMst1) - 2880*lmMst2*(-237 + 55*
        lmMst1 + 4*pow2(lmMst1)) - 1440*(-280 + 121*lmMst1)*pow2(lmMst2) +
        185760*pow3(lmMst2))))) - 36*(5*Mst2*(9*Mst2*Mt*(10667 - 96*B4 + 96*D3
        - 48*DN - 3072*lmMst1 - 384*(-25 + 9*lmMst1)*lmMst2 - 384*(-1 + lmMst1
        - lmMst2)*lmMt - 224*OepS2 + 324*(-43 + 14*lmMst1 - 14*lmMst2)*S2 -
        384*(-13 + 4*lmMst1)*pow2(lmMst2) + 1536*pow3(lmMst2)) + s2t*(12960*
        Dmsqst2*(1 + 3*lmMst1 - 3*lmMst2) + 2*pow2(Mst2)*(75569 + 13716*B4 -
        54*DN - 33426*lmMst1 - 1088*OepS2 + 162*(169 + 136*lmMst1 - 136*lmMst2)
        *S2 - 2376*pow2(lmMst1) + 54*lmMst2*(1427 - 1012*lmMst1 + 16*pow2(
        lmMst1)) - 108*(-642 + 203*lmMst1)*pow2(lmMst2) + 21060*pow3(lmMst2))))
        + Dmglst2*(30*Mst2*Mt*(28405 - 288*B4 + 288*D3 - 144*DN + 5472*lmMst1 -
        288*(-89 + 2*lmMst1)*lmMst2 + 576*(-2 + lmMst1 - lmMst2)*lmMt + 224*
        OepS2 - 324*(65 + 14*lmMst1 - 14*lmMst2)*S2 - 576*(-9 + 8*lmMst1)*pow2(
        lmMst2) + 4608*pow3(lmMst2)) + s2t*(64800*Dmsqst2*(8 - 15*lmMst1 + 15*
        lmMst2) + pow2(Mst2)*(66761 + 301320*B4 - 4860*DN - 205380*lmMst1 +
        40480*OepS2 - 216*(2489 + 3795*lmMst1 - 3795*lmMst2)*S2 + 23760*pow2(
        lmMst1) + 180*lmMst2*(4993 - 1956*lmMst1 + 48*pow2(lmMst1)) - 1080*(-
        482 + 331*lmMst1)*pow2(lmMst2) + 348840*pow3(lmMst2)))))*pow4(Mst1) -
        10*(9*Mt*(383185 - 2592*B4 + 2592*D3 - 1296*DN - 187704*lmMst1 - 17440*
        OepS2 + 648*(-57 + 545*lmMst1 - 545*lmMst2)*S2 - 7992*pow2(lmMst1) -
        216*lmMst2*(-1733 + 630*lmMst1 + 26*pow2(lmMst1)) - 216*(-859 + 246*
        lmMst1)*pow2(lmMst2) + 3456*lmMt*(3 + 6*lmMst2 - 2*lmMst1*(3 + lmMst2)
        + pow2(lmMst1) + pow2(lmMst2)) + 720*pow3(lmMst1) + 58032*pow3(lmMst2))
        + s2t*(-(Dmglst2*(2773621 - 1660176*B4 + 25272*DN + 2004408*lmMst1 -
        525472*OepS2 + 108*(123113 + 98526*lmMst1 - 98526*lmMst2)*S2 + 3888*
        pow2(lmMst1) - 144*lmMst2*(36802 - 11421*lmMst1 + 1728*pow2(lmMst1)) +
        167184*(-14 + 15*lmMst1)*pow2(lmMst2) - 31104*pow3(lmMst1) - 2227824*
        pow3(lmMst2))) + 3*Mst2*(1702429 + 257904*B4 - 648*DN - 748656*lmMst1 -
        30112*OepS2 + 108*(9185 + 5646*lmMst1 - 5646*lmMst2)*S2 + 41904*pow2(
        lmMst1) + 216*lmMst2*(5971 - 6106*lmMst1 + 576*pow2(lmMst1)) - 41904*(-
        34 + 15*lmMst1)*pow2(lmMst2) - 3456*pow3(lmMst1) + 507600*pow3(lmMst2))
        ))*pow6(Mst1) - 622080*(1 + lmMst2)*(-Dmglst2 + 3*Dmglst2*lmMst2 + Mst2
        + lmMst2*Mst2)*s2t*pow6(Mst2)))/(43740.*pow2(Mst1)*pow2(Sbeta)*pow6(
        Mst2)) + pow2(Mt)*pow2(s2t)*((5*Dmsqst2*(391 + 32*lmMst1 - 128*lmMst2))
        /18. + pow2(Mst2)*(312.7354938271605 + (8*D3)/9. - (8*DN)/9. + (806*
        lmMst1)/27. + (32*lmMt)/3. + (2*lmMst2*(2267 + 90*lmMst1 - 117*pow2(
        lmMst1)))/27. + (124*pow2(lmMst1))/9. + (4*(133 + 8*lmMst1)*pow2(
        lmMst2))/9. + (20*Dmsqst2*(-1 + 2*lmMst2))/(3.*pow2(Mst1)) + (32*pow3(
        lmMst1))/9. + (14*pow3(lmMst2))/9.) + pow2(Mst1)*(97.4135390946502 - (
        8*D3)/9. + (8*DN)/9. - (88984*lmMst1)/405. - (64*(lmMst1 - lmMst2)*
        lmMt)/3. - (1738*pow2(lmMst1))/27. + (2*lmMst2*(52322 - 4710*lmMst1 +
        315*pow2(lmMst1)))/405. + (4*(479 - 237*lmMst1)*pow2(lmMst2))/27. - (
        260*pow3(lmMst1))/27. + (1166*pow3(lmMst2))/27.) - (Dmglst2*(
        175.16754355781114 - (17578814*lmMst1)/33075. + lmMst2*(
        257.7056386999244 + (105592*lmMst1)/315. - (208*pow2(lmMst1))/9.) - (
        35576*pow2(lmMst1))/315. + (16*(-4376 + 2555*lmMst1)*pow2(lmMst2))/315.
         + (64*lmMt*(2 + 5*lmMst2 - lmMst1*(5 + 2*lmMst2) + pow2(lmMst1) +
        pow2(lmMst2)))/3. + (16*pow3(lmMst1))/27. - (2896*pow3(lmMst2))/27.)*
        pow4(Mst1))/pow3(Mst2) - (Dmsqst2*(17.999506172839506 + lmMst1*(
        32.62222222222222 - 16*lmMst2) - (2068*lmMst2)/45. + 8*pow2(lmMst1) +
        8*pow2(lmMst2))*pow2(Mst1) - (199.98139767323744 - (18614063*lmMst1)/
        66150. + lmMst2*(293.39173091458804 - (25514*lmMst1)/945. - (286*pow2(
        lmMst1))/9.) - (48143*pow2(lmMst1))/945. + (32*lmMt*(-2*lmMst1*(1 +
        lmMst2) + lmMst2*(2 + lmMst2) + pow2(lmMst1)))/3. + (77.94391534391535
         - (2*lmMst1)/
        9.)*pow2(lmMst2) + (190*pow3(lmMst1))/27. + (674*pow3(lmMst2))/27.)*
        pow4(Mst1))/pow2(Mst2) + ((-2*(Mst2*(87 + 154*lmMst2 + 32*lmMst1*(1 +
        lmMst2) + 75*pow2(lmMst2)) + 2*Dmglst2*(-52 + 2*(67 + 16*lmMst1)*lmMst2
        + 91*pow2(lmMst2)))*pow3(Mst2))/9. + (4*Dmglst2*(-33750*Dmsqst2*(23*
        pow2(Mst1) + 18*pow2(Mst2)) + 225*pow2(Mst1)*pow2(Mst2)*(11674 - 120*B4
        + 120*D3 - 60*DN + 690*lmMst1 + 345*pow2(lmMst1) - 5*lmMst2*(-3427 -
        54*lmMst1 + 96*pow2(lmMst1)) + (4515 - 480*lmMst1)*pow2(lmMst2) + 960*
        pow3(lmMst2)) + (1216808 + 1164855*lmMst2 - 162000*(2 + lmMst2)*lmMt +
        225*(-341 + 30*lmMst2)*pow2(lmMst1) + 30*lmMst1*(34484 - 16260*lmMst2 +
        5400*lmMt - 21825*pow2(lmMst2)) + 365400*pow2(lmMst2) - 2250*pow3(
        lmMst1) + 650250*pow3(lmMst2))*pow4(Mst1)))/(30375.*Mst2))/pow2(Mst1) +
        (Dmsqst2*(Mst2*(22430251 + 148309560*lmMst2 + 1260*lmMst1*(-117706 +
        34965*lmMst2) - 22027950*pow2(lmMst1) - 22027950*pow2(lmMst2))*pow4(
        Mst1) - 4630500*Dmglst2*(3*(113 - 6*lmMst1 + 6*lmMst2)*pow2(Mst1)*pow2(
        Mst2) + (328 - 96*lmMst1 + 96*lmMst2)*pow4(Mst1))))/(6.251175e6*pow5(
        Mst2)) + (32*(1 + lmMst2)*(4*Dmglst2*lmMst2 + Mst2 + lmMst2*Mst2)*pow5(
        Mst2))/(9.*pow4(Mst1)) + pow2(MuSUSY)*(53.385802469135804 + (40*B4)/9.
         - (4*D3)/
        9. + (2*DN)/9. + (1672*lmMst1)/27. + (53*pow2(lmMst1))/9. - lmMst2*(
        129.92592592592592 - 72*lmMst1 + (13*pow2(lmMst1))/3.) + (2*(-470 +
        147*lmMst1)*pow2(lmMst2))/9. - (2*Dmglst2*Mst2*(-20 + 2*(99 + 16*
        lmMst1)*lmMst2 + 123*pow2(lmMst2)))/(9.*pow2(Mst1)) - ((103 + 186*
        lmMst2 + 32*lmMst1*(1 + lmMst2) + 91*pow2(lmMst2))*pow2(Mst2))/(9.*
        pow2(Mst1)) + (16*pow3(lmMst1))/9. + (2250*Dmsqst2*(1141 + 48*lmMst1*(7
        - 3*lmMst2) - 264*lmMst2 + 144*pow2(lmMst2)) + pow2(Mst1)*(10552777 +
        205200*B4 - 5400*DN + 4145280*lmMst1 - 1181700*pow2(lmMst1) - 240*
        lmMst2*(16192 - 26430*lmMst1 + 3465*pow2(lmMst1)) + 900*(-5735 + 3072*
        lmMst1)*pow2(lmMst2) - 55800*pow3(lmMst1) - 1877400*pow3(lmMst2)))/(
        24300.*pow2(Mst2)) - (271*pow3(lmMst2))/9. - (Dmglst2*(14267 - 432*B4 +
        576*D3 - 360*DN + 4752*lmMst1 - 1404*pow2(lmMst1) + 72*lmMst2*(-63 -
        232*lmMst1 + 16*pow2(lmMst1)) + 72*(281 - 123*lmMst1)*pow2(lmMst2) + (
        6480*Dmsqst2)/pow2(Mst1) + 7704*pow3(lmMst2)))/(162.*Mst2) - (2*
        Dmglst2*(405000*Dmsqst2*(4 + lmMst1 - lmMst2) + pow2(Mst1)*(2695042 -
        40500*B4 + 54000*D3 - 33750*DN - 326895*lmMst1 - 324900*pow2(lmMst1) +
        15*lmMst2*(-19607 - 129030*lmMst1 + 62550*pow2(lmMst1)) - 450*(-5023 +
        5610*lmMst1)*pow2(lmMst2) + 11250*pow3(lmMst1) + 1575000*pow3(lmMst2)))
        )/(30375.*pow3(Mst2)) + (16*(1 + lmMst2)*(4*Dmglst2*lmMst2 + Mst2 +
        lmMst2*Mst2)*pow3(Mst2))/(9.*pow4(Mst1)) + (Dmsqst2*(201.74098765432097
         - (622*lmMst2)/15. + (2*lmMst1*(311 + 10*lmMst2))/15. - (22*pow2(
        lmMst1))/3. + 6*pow2(lmMst2))*pow2(Mst1) + (628.1736268201578 + (76*B4)
        /9. - (2*DN)/9. + (6317839*lmMst1)/396900. - (66307*pow2(lmMst1))/315.
         - lmMst2*(
        12.52907281431091 - (182909*lmMst1)/315. + (274*pow2(lmMst1))/3.) + (2*
        (-58301 + 37135*lmMst1)*pow2(lmMst2))/315. - (44*pow3(lmMst1))/9. - (
        1256*pow3(lmMst2))/9.)*pow4(Mst1))/pow4(Mst2) - (Dmglst2*(
        585.1892843532082 - (8*B4)/3. + (32*D3)/9. - (20*DN)/9. - (20109937*
        lmMst1)/297675. - (15886*pow2(lmMst1))/945. + lmMst2*(
        17.112243218274966 - (144628*lmMst1)/945. + (1180*pow2(lmMst1))/9.) + (
        169.85608465608465 - (2632*lmMst1)/9.)*pow2(lmMst2) - (92*pow3(lmMst1))
        /27. + (4448*pow3(lmMst2))/27.)*pow4(Mst1))/pow5(Mst2) - (Dmsqst2*(10 -
        20*lmMst2 + (10*Dmglst2*(-23 + 42*lmMst1 - 42*lmMst2)*pow4(Mst1))/pow5(
        Mst2) - (3*(237.28785508324435 + (16526*lmMst2)/3969. + (2*lmMst1*(-
        8263 + 71820*lmMst2))/3969. - (520*pow2(lmMst1))/21. - (80*pow2(lmMst2)
        )/7.)*pow6(Mst1))/pow6(Mst2)))/(3.*pow2(Mst1))) + (-300*Dmsqst2*(9*
        pow2(Mst1)*pow2(Mst2)*(2*(8*OepS2 - 27*(23 + 6*lmMst1 - 6*lmMst2)*S2)*
        pow2(Mst2) + (104*OepS2 + 27*(17 - 78*lmMst1 + 78*lmMst2)*S2)*pow2(
        MuSUSY)) + (4*(16*OepS2 - 27*(35 + 12*lmMst1 - 12*lmMst2)*S2)*pow2(
        Mst2) + 13*(136*OepS2 - 27*(95 + 102*lmMst1 - 102*lmMst2)*S2)*pow2(
        MuSUSY))*pow4(Mst1) + 27*((8*OepS2 - 81*(9 + 2*lmMst1 - 2*lmMst2)*S2)*
        pow2(Mst2) + (8*OepS2 - 81*(-15 + 2*lmMst1 - 2*lmMst2)*S2)*pow2(MuSUSY)
        )*pow4(Mst2)) + Mst2*(-3*Mst2*(6*pow2(Mst1)*pow2(Mst2)*((25160*OepS2 -
        81*(9191 + 6290*lmMst1 - 6290*lmMst2)*S2)*pow2(Mst2) + 5*(8456*OepS2 -
        81*(11243 + 2114*lmMst1 - 2114*lmMst2)*S2)*pow2(MuSUSY)) + ((274280*
        OepS2 - 27*(399127 + 205710*lmMst1 - 205710*lmMst2)*S2)*pow2(Mst2) +
        10*(52948*OepS2 - 27*(194357 + 39711*lmMst1 - 39711*lmMst2)*S2)*pow2(
        MuSUSY))*pow4(Mst1) + 27*((2120*OepS2 - 81*(-141 + 530*lmMst1 - 530*
        lmMst2)*S2)*pow2(Mst2) + 10*(184*OepS2 - 81*(307 + 46*lmMst1 - 46*
        lmMst2)*S2)*pow2(MuSUSY))*pow4(Mst2)) + 4*Dmglst2*(54*pow2(Mst1)*pow2(
        Mst2)*((2000*OepS2 - 162*(31 + 250*lmMst1 - 250*lmMst2)*S2)*pow2(Mst2)
        + 5*(344*OepS2 + 9*(15643 - 774*lmMst1 + 774*lmMst2)*S2)*pow2(MuSUSY))
        + 2*(9*(30760*OepS2 - 27*(28283 + 23070*lmMst1 - 23070*lmMst2)*S2)*
        pow2(Mst2) + 10*(17308*OepS2 + 27*(93919 - 12981*lmMst1 + 12981*lmMst2)
        *S2)*pow2(MuSUSY))*pow4(Mst1) + 135*(56*OepS2 - 81*(-1677 + 14*lmMst1 -
        14*lmMst2)*S2)*pow2(MuSUSY)*pow4(Mst2) - 2768742*S2*pow6(Mst2))))/(
        21870.*pow6(Mst2))) + pow4(Mt)*(480.98395061728394 - (640*lmMst1)/9. -
        (392*pow2(lmMst1))/9. + (4*lmMst2*(-553 - 1224*lmMst1 + 63*pow2(lmMst1)
        ))/81. + (32*(121 - 18*lmMst1)*pow2(lmMst2))/27. + (4*lmMt*(926 - 749*
        lmMst2 + 6*lmMst1*(27 + 32*lmMst2) + 9*pow2(lmMst1) + 57*pow2(lmMst2)))
        /27. - (4*(-12 + lmMst1 + 55*lmMst2)*pow2(lmMt))/3. + (4*(71 + 122*
        lmMst2 + 32*lmMst1*(1 + lmMst2) + 59*pow2(lmMst2))*pow2(Mst2))/(9.*
        pow2(Mst1)) - (68*pow3(lmMst1))/9. + (18000*Dmsqst2*(-103 + 42*lmMst1 +
        78*lmMt - 6*lmMst2*(23 + 3*lmMt) + 9*pow2(lmMst2) + 9*pow2(lmMt)) +
        pow2(Mst1)*(7319011 + 3223800*lmMt - 1800*(391 + 114*lmMst2 - 90*lmMt)*
        pow2(lmMst1) + 3600*(568 + 99*lmMt)*pow2(lmMst2) - 259200*pow2(lmMt) -
        120*lmMst2*(22352 - 2835*lmMt + 4320*pow2(lmMt)) + 120*lmMst1*(4217 +
        1215*lmMt - 15*lmMst2*(871 + 288*lmMt) + 3240*pow2(lmMst2) + 4320*pow2(
        lmMt)) + 14400*pow3(lmMst1) - 198000*pow3(lmMst2)))/(12150.*pow2(Mst2))
        + (4*Dmglst2*(113400*Dmsqst2 + 630*(-84 + (70 + 32*lmMst1)*lmMst2 + 59*
        pow2(lmMst2))*pow2(Mst2) + pow2(Mst1)*(773533 + 131775*lmMt + 35*
        lmMst2*(-6569 + 2220*lmMt) - 4410*pow2(lmMst1) - 199920*pow2(lmMst2) +
        5040*lmMst1*(13 + 28*lmMst2 - 8*lmMt + 8*pow2(lmMst2)) + 63000*pow2(
        lmMt) - 40320*pow3(lmMst2))))/(2835.*Mst2*pow2(Mst1)) + (8*pow3(lmMst2)
        )/9. + (184*pow3(lmMt))/3. + (2*Dmglst2*(52500*Dmsqst2*(2954 - 864*
        lmMst1 + 1713*lmMst2 - 849*lmMt) + pow2(Mst1)*(98884301 - 57865500*lmMt
        + 12600*(1213 + 285*lmMst2 - 225*lmMt)*pow2(lmMst1) + 18900*(-5953 +
        310*lmMt)*pow2(lmMst2) + 13608000*pow2(lmMt) + 1260*(lmMst1*(28194 +
        lmMst2*(61415 - 2400*lmMt) - 1075*lmMt + 13800*pow2(lmMst2) - 7200*
        pow2(lmMt)) + lmMst2*(8581 + 6025*lmMt + 7200*pow2(lmMt))) - 252000*
        pow3(lmMst1) - 20727000*pow3(lmMst2))))/(212625.*pow3(Mst2)) - (64*(1 +
        lmMst2)*(4*Dmglst2*lmMst2 + Mst2 + lmMst2*Mst2)*pow3(Mst2))/(9.*pow4(
        Mst1)) - ((4*Dmsqst2*(91321 + lmMst2*(92010 - 4500*lmMt) - 30*lmMst1*(
        1717 + 120*lmMst2 - 150*lmMt) - 40500*lmMt - 450*pow2(lmMst1) + 4050*
        pow2(lmMst2))*pow2(Mst1))/675. - (1500.0856244066115 + (43574647*
        lmMst1)/99225. - (139432*pow2(lmMst1))/945. - (2*lmMst2*(35585111 -
        2612085*lmMst1 + 2072700*pow2(lmMst1)))/99225. + (2*(34339 + 46620*
        lmMst1)*pow2(lmMst2))/945. + (lmMt*(3155 + 1418*lmMst2 - 2*lmMst1*(513
        + 572*lmMst2) + 312*pow2(lmMst1) + 832*pow2(lmMst2)))/9. + (64*(-1 + 3*
        lmMst1 - 3*lmMst2)*pow2(lmMt))/3. + (64*pow3(lmMst1))/27. - (1600*pow3(
        lmMst2))/27.)*pow4(Mst1))/pow4(Mst2) - (Dmglst2*(2929.938520304849 + (
        55510684*lmMst1)/59535. - (126272*pow2(lmMst1))/189. - (4*lmMst2*(
        42300121 + 12578580*lmMst1 + 2487240*pow2(lmMst1)))/59535. + (32*(10166
        - 693*lmMst1)*pow2(lmMst2))/189. + (8*lmMt*(5695 + 1974*lmMst2 - 12*
        lmMst1*(163 + 47*lmMst2) + 468*pow2(lmMst1) + 96*pow2(lmMst2)))/27. + (
        128*(-5 + 6*lmMst1 - 6*lmMst2)*pow2(lmMt))/3. + (256*pow3(lmMst1))/27.
         + (7424*pow3(lmMst2))/
        27.)*pow4(Mst1))/pow5(Mst2) + Dmsqst2*((40*(1 - 2*lmMst2))/(3.*pow2(
        Mst1)) + (40*Dmglst2*(401 + lmMst2*(274 - 4*lmMt) - 2*lmMst1*(71 + 2*
        lmMst2 - 2*lmMt) - 132*lmMt + 4*pow2(lmMst2))*pow2(Mst1))/(3.*pow5(
        Mst2)) + (8*(-26331136 + 105*lmMst1*(236317 + 35070*lmMst2 - 36750*
        lmMt) + 11962125*lmMt + 210*lmMst2*(-175121 + 18375*lmMt) + 88200*pow2(
        lmMst1) - 3770550*pow2(lmMst2))*pow4(Mst1))/(231525.*pow6(Mst2))) + (-(
        pow2(MuSUSY)*(6*Mst2*pow2(Mst1)*(3*Mst2*(10667 - 96*B4 + 96*D3 - 48*DN
        - 3072*lmMst1 - 384*(-25 + 9*lmMst1)*lmMst2 - 384*(-1 + lmMst1 -
        lmMst2)*lmMt - 384*(-13 + 4*lmMst1)*pow2(lmMst2) + 1536*pow3(lmMst2)) +
        2*Dmglst2*(28405 - 288*B4 + 288*D3 - 144*DN + 5472*lmMst1 - 288*(-89 +
        2*lmMst1)*lmMst2 + 576*(-2 + lmMst1 - lmMst2)*lmMt - 576*(-9 + 8*
        lmMst1)*pow2(lmMst2) + 4608*pow3(lmMst2))) + 144*(6*Dmglst2*(180 - 2*B4
        + 2*D3 - DN + 16*lmMst1 + 144*lmMst2 - 16*(-2 + lmMst1)*pow2(lmMst2) +
        16*pow3(lmMst2)) + Mst2*(436 - 6*B4 + 6*D3 - 3*DN - 48*lmMst1 + (408 -
        96*lmMst1)*lmMst2 + 24*lmMt - 48*(-4 + lmMst1)*pow2(lmMst2) + 48*pow3(
        lmMst2)))*pow3(Mst2) + (383185 - 2592*B4 + 2592*D3 - 1296*DN - 187704*
        lmMst1 - 7992*pow2(lmMst1) - 216*lmMst2*(-1733 + 630*lmMst1 + 26*pow2(
        lmMst1)) - 216*(-859 + 246*lmMst1)*pow2(lmMst2) + 3456*lmMt*(3 + 6*
        lmMst2 - 2*lmMst1*(3 + lmMst2) + pow2(lmMst1) + pow2(lmMst2)) + 720*
        pow3(lmMst1) + 58032*pow3(lmMst2))*pow4(Mst1)))/486. + (16*OepS2*(4*
        Dmglst2*Mst2*(63*pow2(Mst1)*(23*pow2(Mst2) - 3*pow2(MuSUSY)) + 5582*
        pow4(Mst1) + 189*pow4(Mst2)) + 3*((1168*pow2(Mst2) + 1635*pow2(MuSUSY))
        *pow4(Mst1) + pow2(Mst1)*(378*pow2(Mst2)*pow2(MuSUSY) + 627*pow4(Mst2))
        + 189*pow6(Mst2))))/2187. - (S2*(4*Dmglst2*Mst2*(45*pow2(Mst1)*(23*(115
        + 588*lmMst1 - 588*lmMst2)*pow2(Mst2) - 126*(65 + 14*lmMst1 - 14*
        lmMst2)*pow2(MuSUSY)) + (2001242 + 2344440*lmMst1 - 2344440*lmMst2)*
        pow4(Mst1) - 81*(3360*pow2(Mst2)*pow2(MuSUSY) + (1593 - 980*lmMst1 +
        980*lmMst2)*pow4(Mst2))) + 42*((8*(8653 + 4380*lmMst1 - 4380*lmMst2)*
        pow2(Mst2) + 90*(-57 + 545*lmMst1 - 545*lmMst2)*pow2(MuSUSY))*pow4(
        Mst1) + 9*pow2(Mst1)*(90*(-43 + 14*lmMst1 - 14*lmMst2)*pow2(Mst2)*pow2(
        MuSUSY) + (4163 + 2090*lmMst1 - 2090*lmMst2)*pow4(Mst2)) + 81*(-240*
        pow2(MuSUSY)*pow4(Mst2) + (109 + 70*lmMst1 - 70*lmMst2)*pow6(Mst2)))))/
        2835.)/pow6(Mst2)) - (s2t*pow3(Mt)*(18*Dmglst2*pow2(Mst2)*(-3*pow2(
        Mst2)*(322421 + 414720*lmMt - 73760*OepS2 + 564876*S2 - 2880*(21 + 12*
        lmMst2 - 4*lmMt)*pow2(lmMst1) + 480*(1097 - 24*lmMt)*pow2(lmMst2) -
        69120*pow2(lmMt) - 40*lmMst2*(-11408 + 924*lmMt + 37341*S2 + 864*pow2(
        lmMt)) + 40*lmMst1*(-11876 - 8844*lmMst2 - 156*lmMt + 37341*S2 + 288*
        pow2(lmMst2) + 864*pow2(lmMt)) + 23040*pow3(lmMst2)) + 2*pow2(MuSUSY)*(
        66761 + 301320*B4 - 4860*DN - 205380*lmMst1 + 40480*OepS2 - 216*(2489 +
        3795*lmMst1 - 3795*lmMst2)*S2 + 23760*pow2(lmMst1) + 180*lmMst2*(4993 -
        1956*lmMst1 + 48*pow2(lmMst1)) - 1080*(-482 + 331*lmMst1)*pow2(lmMst2)
        + 348840*pow3(lmMst2)))*pow4(Mst1) + 18*(20*pow2(MuSUSY)*(75569 +
        13716*B4 - 54*DN - 33426*lmMst1 - 1088*OepS2 + 162*(169 + 136*lmMst1 -
        136*lmMst2)*S2 - 2376*pow2(lmMst1) + 54*lmMst2*(1427 - 1012*lmMst1 +
        16*pow2(lmMst1)) - 108*(-642 + 203*lmMst1)*pow2(lmMst2) + 21060*pow3(
        lmMst2)) + pow2(Mst2)*(1156193 + 198720*lmMt - 60320*OepS2 + 1414908*S2
        - 3240*lmMst2*(42 + 32*lmMt*(1 + lmMt) + 377*S2) + 8640*(-13 + 4*lmMst2
        + 4*lmMt)*pow2(lmMst1) + 17280*(51 - 2*lmMt)*pow2(lmMst2) - 1080*
        lmMst1*(590 + 720*lmMst2 - 104*lmMt - 1131*S2 + 224*pow2(lmMst2) - 96*
        pow2(lmMt)) + 207360*pow3(lmMst2)))*pow3(Mst2)*pow4(Mst1) - 2*Dmglst2*(
        5*pow2(MuSUSY)*(2773621 - 1660176*B4 + 25272*DN + 2004408*lmMst1 -
        525472*OepS2 + 108*(123113 + 98526*lmMst1 - 98526*lmMst2)*S2 + 3888*
        pow2(lmMst1) - 144*lmMst2*(36802 - 11421*lmMst1 + 1728*pow2(lmMst1)) +
        167184*(-14 + 15*lmMst1)*pow2(lmMst2) - 31104*pow3(lmMst1) - 2227824*
        pow3(lmMst2)) + 4*pow2(Mst2)*(20964397 + 4563540*lmMt - 2058400*OepS2 +
        47545272*S2 - 9720*(266 + 74*lmMst2 - 55*lmMt)*pow2(lmMst1) + 3240*(
        1493 + 39*lmMt)*pow2(lmMst2) + 360*lmMst2*(4724 + 5328*lmMt - 115785*S2
        - 1944*pow2(lmMt)) - 466560*pow2(lmMt) + 180*lmMst1*(-32857 - 11844*
        lmMt - 18*lmMst2*(485 + 204*lmMt) + 231570*S2 + 1674*pow2(lmMst2) +
        3888*pow2(lmMt)) + 9720*pow3(lmMst1) + 408240*pow3(lmMst2)))*pow6(Mst1)
        + 6*Mst2*(4*pow2(Mst2)*(2016907 + 110160*lmMt - 109600*OepS2 + 3520152*
        S2 + 1080*(-158 + 6*lmMst2 + 39*lmMt)*pow2(lmMst1) + 1080*(565 - 3*
        lmMt)*pow2(lmMst2) - 1080*lmMst1*(421 + 60*lmMt + lmMst2*(413 + 36*
        lmMt) - 2055*S2 + 129*pow2(lmMst2) - 72*pow2(lmMt)) - 1080*lmMst2*(167
        - 66*lmMt + 2055*S2 + 72*pow2(lmMt)) + 1080*pow3(lmMst1) + 131760*pow3(
        lmMst2)) + 5*pow2(MuSUSY)*(1702429 + 257904*B4 - 648*DN - 748656*lmMst1
        - 30112*OepS2 + 108*(9185 + 5646*lmMst1 - 5646*lmMst2)*S2 + 41904*pow2(
        lmMst1) + 216*lmMst2*(5971 - 6106*lmMst1 + 576*pow2(lmMst1)) - 41904*(-
        34 + 15*lmMst1)*pow2(lmMst2) - 3456*pow3(lmMst1) + 507600*pow3(lmMst2))
        )*pow6(Mst1) + 622080*Dmglst2*(-1 + 2*lmMst2 + 3*pow2(lmMst2))*(2*pow2(
        Mst2) + pow2(MuSUSY))*pow6(Mst2) - Dmglst2*pow2(Mst1)*(-27*pow2(MuSUSY)
        *(23917 + 188640*B4 - 3600*DN - 37440*lmMst1 + 5600*OepS2 - 324*(-453 +
        350*lmMst1 - 350*lmMst2)*S2 + 11520*pow2(lmMst1) - 2880*lmMst2*(-237 +
        55*lmMst1 + 4*pow2(lmMst1)) - 1440*(-280 + 121*lmMst1)*pow2(lmMst2) +
        185760*pow3(lmMst2))*pow4(Mst2) - 64800*Dmsqst2*(-36*(2 + 5*lmMst1 - 5*
        lmMst2)*pow2(Mst2)*pow2(MuSUSY) + 36*pow2(Mst1)*(8*(7 - 3*lmMst1 + 6*
        lmMst2 - 3*lmMt)*pow2(Mst2) + (8 - 15*lmMst1 + 15*lmMst2)*pow2(MuSUSY))
        + 9*(465 + 374*lmMst2 - 4*lmMst1*(52 + 3*lmMst2 - 3*lmMt) - 2*(83 + 6*
        lmMst2)*lmMt + 12*pow2(lmMst2))*pow4(Mst1) + (578 - 360*lmMst1 + 708*
        lmMst2 - 348*lmMt)*pow4(Mst2)) + 270*(30183 + 45408*lmMt - 1120*OepS2 -
        43740*S2 - 8*lmMst2*(-8398 + 588*lmMt + 2835*S2) + 2304*(-1 + lmMst2)*
        pow2(lmMst1) + 96*(199 + 48*lmMt)*pow2(lmMst2) + 72*lmMst1*(-120 + 64*
        lmMt - 8*lmMst2*(5 + 8*lmMt) + 315*S2 + 164*pow2(lmMst2)) - 14112*pow3(
        lmMst2))*pow6(Mst2)) + 81*Mst2*pow2(Mst1)*(5*pow2(MuSUSY)*(28683 +
        5280*B4 - 48*DN - 5952*lmMst1 - 224*OepS2 + 324*(-1 + 14*lmMst1 - 14*
        lmMst2)*S2 - 768*pow2(lmMst1) - 192*lmMst2*(-214 + 73*lmMst1 + 4*pow2(
        lmMst1)) - 96*(-268 + 57*lmMst1)*pow2(lmMst2) + 6240*pow3(lmMst2))*
        pow4(Mst2) - 7200*Dmsqst2*(-4*(1 + lmMst1 - lmMst2)*pow2(Mst2)*pow2(
        MuSUSY) + 4*pow2(Mst1)*((6 - 4*lmMst1 + 8*lmMst2 - 4*lmMt)*pow2(Mst2) +
        (-1 - 3*lmMst1 + 3*lmMst2)*pow2(MuSUSY)) + (35 + 50*lmMst2 - 4*lmMst1*(
        8 + lmMst2 - lmMt) - 2*(9 + 2*lmMst2)*lmMt + 4*pow2(lmMst2))*pow4(Mst1)
        - 2*(7 + 4*lmMst1 - 10*lmMst2 + 6*lmMt)*pow4(Mst2)) + 2*(60759 - 12480*
        lmMt - 1120*OepS2 - 25596*S2 - 120*lmMst2*(26 + 132*lmMt + 189*S2) -
        3840*(1 + lmMst2)*pow2(lmMst1) + 480*(155 - 16*lmMt)*pow2(lmMst2) -
        120*lmMst1*(216 + 8*lmMst2*(61 - 8*lmMt) - 64*lmMt - 189*S2 + 164*pow2(
        lmMst2)) - 11520*pow2(lmMt) + 23520*pow3(lmMst2))*pow6(Mst2)) + 622080*
        pow2(1 + lmMst2)*(2*pow2(Mst2) + pow2(MuSUSY))*pow7(Mst2)))/(43740.*
        pow2(Mst1)*pow6(Mst2)) + ((-10*shiftst1*(Dmsqst2*(3 - 2*lmMst2) + (1 -
        2*lmMst2)*pow2(Mst2))*(pow4(Mst2)*(pow2(s2t)*(-8*pow2(Mt) + (pow2(Mst1)
        - pow2(Mst2))*pow2(s2t))*pow2(Sbeta)*pow4(Mst1) + pow2(Mst1)*(4*pow2(
        Mt)*pow2(s2t)*(pow2(MuSUSY) + (4*pow2(Mst2) - pow2(MuSUSY))*pow2(Sbeta)
        ) + pow2(Sbeta)*(16*pow4(Mt) - pow4(Mst2)*pow4(s2t))) + pow2(Mst2)*(-4*
        pow2(Mt)*pow2(s2t)*(pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 2*pow2(Mst2)*
        pow2(Sbeta)) + pow2(Sbeta)*(16*pow4(Mt) + pow4(Mst2)*pow4(s2t)))) + (-
        2*lmMst1 + 2*lmMst2)*pow2(Mst1)*pow2(s2t)*(-4*pow2(Mt)*pow2(MuSUSY)*(-1
        + pow2(Sbeta))*(pow4(Mst1) + pow4(Mst2)) - pow2(Mst1)*pow2(Mst2)*(4*
        pow2(Mt)*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + pow2(s2t)*pow2(Sbeta)*pow4(
        Mst2)) + pow2(s2t)*pow2(Sbeta)*pow8(Mst2))))/pow2(Sbeta) + (1 - 2*
        lmMst2)*shiftst3*pow2(Mst2)*((8*(-1 + lmMst1 - lmMst2)*pow2(Mst1)*pow2(
        Mt)*pow2(s2t) - 16*pow4(Mt) - (-4*pow2(Mst1)*pow2(Mst2) + (3 + 2*lmMst1
        - 2*lmMst2)*pow4(Mst1) + pow4(Mst2))*pow4(s2t))*pow6(Mst2) - (4*pow2(
        Mt)*pow2(MuSUSY)*pow2(s2t)*(2*(1 - 2*lmMst1 + 2*lmMst2)*pow2(Mst2)*
        pow4(Mst1) + 2*(1 - lmMst1 + lmMst2)*pow2(Mst1)*pow4(Mst2) + (2 - 6*
        lmMst1 + 6*lmMst2)*pow6(Mst1) + pow6(Mst2)))/pow2(Sbeta) + 4*pow2(Mt)*
        pow2(s2t)*(pow2(MuSUSY)*pow6(Mst2) + 2*(pow2(MuSUSY)*((1 - 2*lmMst1 +
        2*lmMst2)*pow2(Mst2)*pow4(Mst1) + (1 - lmMst1 + lmMst2)*pow2(Mst1)*
        pow4(Mst2) + (1 - 3*lmMst1 + 3*lmMst2)*pow6(Mst1)) + pow8(Mst2)))))/(
        12.*pow2(Mst1)*pow6(Mst2)) + (Mt*pow3(s2t)*(-(Dmglst2*(-18*(27211 +
        36720*B4 + 1080*DN - 298440*lmMst1 + 64160*OepS2 - 108*(14033 + 12030*
        lmMst1 - 12030*lmMst2)*S2 + 12960*pow2(lmMst1) + 360*lmMst2*(-503 -
        636*lmMst1 + 144*pow2(lmMst1)) - 2160*(30 + 89*lmMst1)*pow2(lmMst2) +
        140400*pow3(lmMst2))*pow4(Mst1)*pow4(Mst2) + pow2(Mst2)*(31897243 -
        2491360*OepS2 + 90290268*S2 - 360*lmMst2*(18652 + 140139*S2) + 38880*(
        37 - 40*lmMst2)*pow2(lmMst1) + 3188160*pow2(lmMst2) + 360*lmMst1*(17410
        - 12852*lmMst2 + 140139*S2 + 11232*pow2(lmMst2)) - 311040*pow3(lmMst1)
        - 2177280*pow3(lmMst2))*pow6(Mst1) + pow2(Mst1)*(-583200*Dmsqst2*(4*(12
        - 5*lmMst1 + 5*lmMst2)*pow2(Mst1)*pow2(Mst2) + (40 - 3*lmMst1 + 3*
        lmMst2)*pow4(Mst1) - 4*(2 + 5*lmMst1 - 5*lmMst2)*pow4(Mst2)) - 27*(
        69997 + 188640*B4 - 3600*DN - 37440*lmMst1 + 5600*OepS2 - 324*(-453 +
        350*lmMst1 - 350*lmMst2)*S2 + 11520*pow2(lmMst1) - 2880*lmMst2*(-205 +
        55*lmMst1 + 4*pow2(lmMst1)) - 1440*(-184 + 121*lmMst1)*pow2(lmMst2) +
        185760*pow3(lmMst2))*pow6(Mst2)) - 622080*(-1 + 2*lmMst2 + 3*pow2(
        lmMst2))*pow8(Mst2))) + 15*(-38880*Dmsqst2*Mst2*pow2(Mst1)*(4*(1 -
        lmMst1 + lmMst2)*pow2(Mst1)*pow2(Mst2) + (4 - lmMst1 + lmMst2)*pow4(
        Mst1) - 4*(1 + lmMst1 - lmMst2)*pow4(Mst2)) + 6*(51041 + 7344*B4 + 216*
        DN - 80136*lmMst1 - 2336*OepS2 + 324*(347 + 146*lmMst1 - 146*lmMst2)*S2
        - 2592*pow2(lmMst1) + 216*lmMst2*(-221 - 428*lmMst1 + 48*pow2(lmMst1))
        - 432*(-122 + 89*lmMst1)*pow2(lmMst2) + 28080*pow3(lmMst2))*pow4(Mst1)*
        pow5(Mst2) + (551987 - 14048*OepS2 + 661068*S2 - 216*lmMst2*(46 + 1317*
        S2) + 864*(205 + 216*lmMst2)*pow2(lmMst1) + 216000*pow2(lmMst2) - 216*
        lmMst1*(248 + 1820*lmMst2 - 1317*S2 + 1632*pow2(lmMst2)) - 6912*pow3(
        lmMst1) + 172800*pow3(lmMst2))*pow3(Mst2)*pow6(Mst1) + 27*pow2(Mst1)*(
        25611 + 5280*B4 - 48*DN - 5952*lmMst1 - 224*OepS2 + 324*(-1 + 14*lmMst1
        - 14*lmMst2)*S2 - 768*pow2(lmMst1) - 192*lmMst2*(-182 + 73*lmMst1 + 4*
        pow2(lmMst1)) - 96*(-236 + 57*lmMst1)*pow2(lmMst2) + 6240*pow3(lmMst2))
        *pow7(Mst2) + 41472*pow2(1 + lmMst2)*pow9(Mst2))))/(87480.*pow2(Mst1)*
        pow4(Mst2)) + (-(T1ep*(pow4(Mst1)*(-4*pow2(Mt)*pow2(s2t)*(60*Dmsqst2*(
        221*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 8*pow2(Mst2)*pow2(Sbeta)) + Mst2*
        (39711*Mst2*pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 8*Dmglst2*(4327*pow2(
        MuSUSY)*(-1 + pow2(Sbeta)) + 6921*pow2(Mst2)*pow2(Sbeta)) + 20571*pow2(
        Sbeta)*pow3(Mst2))) + 16*s2t*(-16421*Dmglst2*pow2(MuSUSY)*(-1 + pow2(
        Sbeta)) + 2823*Mst2*pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 51460*Dmglst2*
        pow2(Mst2)*pow2(Sbeta) + 8220*pow2(Sbeta)*pow3(Mst2))*pow3(Mt) + 4*(
        15571*Dmglst2 - 1317*Mst2)*Mt*pow2(Sbeta)*pow3(s2t)*pow4(Mst2) + 16*(
        4905*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 8*Mst2*(2791*Dmglst2 + 438*Mst2)
        *pow2(Sbeta))*pow4(Mt) + (840*Dmsqst2 + Mst2*(-16796*Dmglst2 + 5385*
        Mst2))*pow2(Sbeta)*pow4(Mst2)*pow4(s2t)) - 54*pow4(Mst2)*(-14*Dmglst2*(
        4*Mst2*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta)) - 10*s2t*(
        pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 2*pow2(Mst2)*pow2(Sbeta))*pow3(Mt) +
        5*Mt*pow2(Sbeta)*pow3(s2t)*pow4(Mst2) + 16*Mst2*pow2(Sbeta)*pow4(Mt) -
        pow2(Sbeta)*pow4(s2t)*pow5(Mst2)) - 3*(-40*Dmsqst2*pow2(Mt)*pow2(
        MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta)) + 4*pow2(Mst2)*pow2(Mt)*(-23*pow2(
        MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta)) + 2*(7*pow2(Mt) - 5*Dmsqst2*pow2(
        s2t))*pow2(Sbeta)) + 28*Mst2*s2t*pow2(MuSUSY)*(-1 + pow2(Sbeta))*pow3(
        Mt) + 56*s2t*pow2(Sbeta)*pow3(Mst2)*pow3(Mt) + 2*pow2(s2t)*(-53*pow2(
        Mt) + 5*Dmsqst2*pow2(s2t))*pow2(Sbeta)*pow4(Mst2) - 14*Mt*pow2(Sbeta)*
        pow3(s2t)*pow5(Mst2) + 23*pow2(Sbeta)*pow4(s2t)*pow6(Mst2))) + 18*Mst2*
        pow2(Mst1)*(Mst2*(8*Mst2*s2t*(136*pow2(MuSUSY)*(-1 + pow2(Sbeta)) +
        377*pow2(Mst2)*pow2(Sbeta))*pow3(Mt) - 4*pow2(Mt)*pow2(s2t)*(390*
        Dmsqst2*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 1057*pow2(Mst2)*pow2(MuSUSY)*
        (-1 + pow2(Sbeta)) + 60*Dmsqst2*pow2(Mst2)*pow2(Sbeta) + 629*pow2(
        Sbeta)*pow4(Mst2)) + 8*(126*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 209*pow2(
        Mst2)*pow2(Sbeta))*pow4(Mt) + (210*Dmsqst2 + 643*pow2(Mst2))*pow2(
        Sbeta)*pow4(Mst2)*pow4(s2t) - 292*Mt*pow2(Sbeta)*pow3(s2t)*pow5(Mst2))
        - 4*Dmglst2*(-12*pow2(Mst2)*pow2(Mt)*pow2(s2t)*(43*pow2(MuSUSY)*(-1 +
        pow2(Sbeta)) + 50*pow2(Mst2)*pow2(Sbeta)) + 2*Mst2*s2t*(506*pow2(
        MuSUSY)*(-1 + pow2(Sbeta)) + 1383*pow2(Mst2)*pow2(Sbeta))*pow3(Mt) +
        56*(3*pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 23*pow2(Mst2)*pow2(Sbeta))*
        pow4(Mt) - 401*Mt*pow2(Sbeta)*pow3(s2t)*pow5(Mst2) + 108*pow2(Sbeta)*
        pow4(s2t)*pow6(Mst2)))))/(1458.*pow2(Sbeta)) + (Mt*MuSUSY*(500094000*
        s2t*pow2(Mst1)*pow4(Mst2)*(-10*shiftst1*((1 - 2*lmMst2)*pow2(Mst2)*(4*
        pow2(Mst2)*pow2(Mt) - 2*pow2(Mst1)*(2*pow2(Mt) + (-lmMst1 + lmMst2)*
        pow2(Mst2)*pow2(s2t)) + pow2(s2t)*pow4(Mst1)) + Dmsqst2*(3 - 2*lmMst2)*
        (-4*pow2(Mst1)*pow2(Mt) + 4*pow2(Mst2)*pow2(Mt) + 2*(lmMst1 - lmMst2)*
        pow2(Mst1)*pow2(Mst2)*pow2(s2t) + pow2(s2t)*(pow4(Mst1) - pow4(Mst2))))
        + (1 - 2*lmMst2)*((10*shiftst1 + shiftst3)*pow2(Mst2)*pow2(s2t) -
        shiftst3*(4*pow2(Mt) + (1 + lmMst1 - lmMst2)*pow2(Mst1)*pow2(s2t)))*
        pow4(Mst2)) + 68600*pow2(Mst1)*pow3(Mt)*(-(Dmglst2*(216*pow2(Mst2)*(
        97837 + 71580*lmMt - 9920*OepS2 + 43272*S2 - 360*(23 + 8*lmMst2 - 4*
        lmMt)*pow2(lmMst1) + 720*(97 + 2*lmMt)*pow2(lmMst2) - 8640*pow2(lmMt) -
        30*lmMst2*(-2911 + 396*lmMt + 6696*S2 + 144*pow2(lmMt)) + 10*lmMst1*(-
        6157 + 642*lmMt - 6*lmMst2*(791 + 48*lmMt) + 20088*S2 + 882*pow2(
        lmMst2) + 432*pow2(lmMt)) - 5940*pow3(lmMst2))*pow4(Mst1) - 129600*
        Dmsqst2*((101 - 90*lmMst1 + 177*lmMst2 - 87*lmMt)*pow2(Mst1)*pow2(Mst2)
        + (497 - 306*lmMst1 + 609*lmMst2 - 303*lmMt)*pow4(Mst1)) + 135*pow2(
        Mst1)*(57495 + 45408*lmMt - 1120*OepS2 - 43740*S2 - 8*lmMst2*(-6952 +
        588*lmMt + 2835*S2) + 2304*(-1 + lmMst2)*pow2(lmMst1) + 96*(79 + 48*
        lmMt)*pow2(lmMst2) + 72*lmMst1*(-88 + 64*lmMt - 8*lmMst2*(9 + 8*lmMt) +
        315*S2 + 164*pow2(lmMst2)) - 14112*pow3(lmMst2))*pow4(Mst2) + 20*(
        5659217 + 1592460*lmMt - 518816*OepS2 + 9976392*S2 - 972*(569 + 180*
        lmMst2 - 126*lmMt)*pow2(lmMst1) + 324*(5353 + 126*lmMt)*pow2(lmMst2) -
        186624*pow2(lmMt) + 72*lmMst1*(-27653 - 3015*lmMt - 18*lmMst2*(689 +
        126*lmMt) + 145917*S2 + 2160*pow2(lmMst2) + 2592*pow2(lmMt)) - 36*
        lmMst2*(-39031 - 3204*lmMt + 291834*S2 + 5184*pow2(lmMt)) + 1944*pow3(
        lmMst1) + 17496*pow3(lmMst2))*pow6(Mst1) - 622080*(-1 + 2*lmMst2 + 3*
        pow2(lmMst2))*pow6(Mst2))) + 3*Mst2*(48*pow2(Mst2)*(110219 - 1080*lmMt
        - 4400*OepS2 + 74034*S2 + 540*(-15 + 4*lmMt)*pow2(lmMst1) + 810*(121 -
        8*lmMt)*pow2(lmMst2) - 135*lmMst1*(337 + lmMst2*(588 - 32*lmMt) - 132*
        lmMt - 660*S2 + 194*pow2(lmMst2) - 48*pow2(lmMt)) - 6480*pow2(lmMt) -
        405*lmMst2*(31 + 54*lmMt + 220*S2 + 16*pow2(lmMt)) + 26190*pow3(lmMst2)
        )*pow4(Mst1) + 388800*Dmsqst2*((5 + 2*lmMst1 - 5*lmMst2 + 3*lmMt)*pow2(
        Mst1)*pow2(Mst2) + (1 + 6*lmMst1 - 13*lmMst2 + 7*lmMt)*pow4(Mst1)) +
        27*pow2(Mst1)*(56439 - 24000*lmMt - 1120*OepS2 - 25596*S2 - 360*lmMst2*
        (-12 + 44*lmMt + 63*S2) - 3840*(1 + lmMst2)*pow2(lmMst1) + 480*(163 -
        16*lmMt)*pow2(lmMst2) - 120*lmMst1*(184 + 8*lmMst2*(57 - 8*lmMt) - 64*
        lmMt - 189*S2 + 164*pow2(lmMst2)) - 11520*pow2(lmMt) + 23520*pow3(
        lmMst2))*pow4(Mst2) + 20*(678923 + 19440*lmMt - 32480*OepS2 + 881712*S2
        + 108*(-457 + 12*lmMst2 + 126*lmMt)*pow2(lmMst1) + 540*(661 - 30*lmMt)*
        pow2(lmMst2) - 108*lmMst1*(1841 + lmMst2*(2626 - 24*lmMt) - 420*lmMt -
        6090*S2 + 840*pow2(lmMst2) - 288*pow2(lmMt)) - 15552*pow2(lmMt) - 108*
        lmMst2*(619 + 498*lmMt + 6090*S2 + 288*pow2(lmMt)) + 216*pow3(lmMst1) +
        89208*pow3(lmMst2))*pow6(Mst1) + 207360*pow2(1 + lmMst2)*pow6(Mst2))) -
        4116000*T1ep*pow4(Mst1)*(-3*(27*pow4(Mst2)*(-40*Dmsqst2*s2t*pow2(Mt) -
        21*Mt*pow2(s2t)*pow3(Mst2) + 28*Mst2*pow3(Mt) + pow2(Mst2)*(-106*s2t*
        pow2(Mt) + 20*Dmsqst2*pow3(s2t)) + 46*pow3(s2t)*pow4(Mst2)) + 3*pow2(
        Mst1)*pow2(Mst2)*(-600*Dmsqst2*s2t*pow2(Mt) - 627*Mt*pow2(s2t)*pow3(
        Mst2) + 1760*Mst2*pow3(Mt) + pow2(Mst2)*(-3470*s2t*pow2(Mt) + 600*
        Dmsqst2*pow3(s2t)) + 1700*pow3(s2t)*pow4(Mst2)) + pow4(Mst1)*(-2120*
        Dmsqst2*s2t*pow2(Mt) - 3198*Mt*pow2(s2t)*pow3(Mst2) + 16240*Mst2*pow3(
        Mt) - 4*pow2(Mst2)*(6031*s2t*pow2(Mt) - 520*Dmsqst2*pow3(s2t)) + 6895*
        pow3(s2t)*pow4(Mst2))) + Dmglst2*(27*pow2(Mst1)*pow2(Mst2)*(-800*Mst2*
        s2t*pow2(Mt) - 907*Mt*pow2(Mst2)*pow2(s2t) + 1984*pow3(Mt) + 316*pow3(
        Mst2)*pow3(s2t)) + 2*(-66168*Mst2*s2t*pow2(Mt) - 35601*Mt*pow2(Mst2)*
        pow2(s2t) + 129704*pow3(Mt) + 12664*pow3(Mst2)*pow3(s2t))*pow4(Mst1) +
        189*(20*pow3(Mt)*pow4(Mst2) - 15*Mt*pow2(s2t)*pow6(Mst2) + 4*pow3(s2t)*
        pow7(Mst2)))) + 12*s2t*pow2(Mt)*(20*Dmsqst2*(1157625*(1968*Dmglst2 +
        Mst2*(-1101 + 240*lmMst2 + 32*OepS2 - 2916*S2 + 648*lmMst2*S2 - 24*
        lmMst1*(4 + 27*S2)))*pow3(Mst2)*pow4(Mst1) + 83349000*(12*Dmglst2 +
        Mst2 - 2*lmMst2*Mst2)*pow2(Mst1)*pow5(Mst2) + 3087*Mst2*(27000*Dmglst2*
        (65 - 2*lmMst1 + 2*lmMst2) - Mst2*(339977 - 20000*OepS2 + 1080*lmMst2*(
        89 - 375*S2) + 1714500*S2 + 1080*lmMst1*(-89 + 60*lmMst2 + 375*S2) -
        32400*(pow2(lmMst1) + pow2(lmMst2))))*pow6(Mst1) - (1094369501 -
        72716000*OepS2 + 2520*lmMst2*(235453 - 584325*S2) + 5940931500*S2 +
        2520*lmMst1*(-235453 + 114345*lmMst2 + 584325*S2) - 144074700*(pow2(
        lmMst1) + pow2(lmMst2)))*pow8(Mst1)) + Mst2*(-77175*(979423 + 2880*D3 -
        2880*DN + 73680*lmMst1 + 34560*lmMt - 25440*OepS2 + 972*(-141 + 530*
        lmMst1 - 530*lmMst2)*S2 + 44640*pow2(lmMst1) - 240*lmMst2*(-1901 + 6*
        lmMst1 + 117*pow2(lmMst1)) + 720*(207 + 16*lmMst1)*pow2(lmMst2) +
        11520*pow3(lmMst1) + 5040*pow3(lmMst2))*pow4(Mst1)*pow5(Mst2) - 5145*(
        19425643 + 518400*lmMt - 1388000*OepS2 + 120*lmMst2*(165994 + 8640*lmMt
        - 234225*S2) + 27723060*S2 - 3600*(683 + 96*lmMst2)*pow2(lmMst1) +
        5684400*pow2(lmMst2) - 120*lmMst1*(84094 + 9600*lmMst2 + 8640*lmMt -
        234225*S2 + 12780*pow2(lmMst2)) - 295200*pow3(lmMst1) + 2174400*pow3(
        lmMst2))*pow3(Mst2)*pow6(Mst1) + 55566000*(71 + 122*lmMst2 + 32*lmMst1*
        (1 + lmMst2) + 59*pow2(lmMst2))*pow2(Mst1)*pow7(Mst2) - Mst2*(
        149949681779 + 2667168000*lmMt - 16549064000*OepS2 + 1260*lmMst2*(
        141677449 + 8467200*lmMt - 265967100*S2) + 512266658400*S2 - 264600*(
        90913 + 36750*lmMst2 - 10080*lmMt)*pow2(lmMst1) + 264600*(189227 +
        10080*lmMt)*pow2(lmMst2) - 1260*lmMst1*(99165049 + 8467200*lmMt + 420*
        lmMst2*(28997 + 10080*lmMt) - 265967100*S2 + 6306300*pow2(lmMst2)) +
        240786000*pow3(lmMst1) + 17429202000*pow3(lmMst2))*pow8(Mst1) - 12*
        Dmglst2*(617400*(12454 - 120*B4 + 120*D3 - 60*DN + 690*lmMst1 - 17091*
        S2 + 345*pow2(lmMst1) - 5*lmMst2*(-3121 + 42*lmMst1 + 96*pow2(lmMst1))
        + (3630 - 480*lmMst1)*pow2(lmMst2) + 960*pow3(lmMst2))*pow4(Mst1)*pow4(
        Mst2) + 2744*pow2(Mst2)*(3856958 - 27000*B4 + 27000*D3 - 13500*DN +
        1270770*lmMst1 + 162000*(-2 + lmMst1 - lmMst2)*lmMt + 150000*OepS2 -
        30375*(139 + 100*lmMst1 - 100*lmMst2)*S2 + 900*pow2(lmMst1) - 30*
        lmMst2*(-153166 + 17835*lmMst1 + 3375*pow2(lmMst1)) - 450*(-2627 +
        1695*lmMst1)*pow2(lmMst2) - 2250*pow3(lmMst1) + 866250*pow3(lmMst2))*
        pow6(Mst1) - 9261000*(-52 + 2*(51 + 16*lmMst1)*lmMst2 + 59*pow2(lmMst2)
        )*pow2(Mst1)*pow6(Mst2) + 5*(1297790971 - 14817600*B4 + 14817600*D3 -
        7408800*DN + 3134593140*lmMst1 + 504347200*OepS2 - 740880*(17269 +
        13785*lmMst1 - 13785*lmMst2)*S2 + 426711600*pow2(lmMst1) + 420*lmMst2*(
        2917823 - 3813600*lmMst1 + 97020*pow2(lmMst1)) - 176400*(-8677 + 5439*
        lmMst1)*pow2(lmMst2) - 88905600*lmMt*(4 + 6*lmMst2 - 2*lmMst1*(3 +
        lmMst2) + pow2(lmMst1) + pow2(lmMst2)) - 3704400*pow3(lmMst1) +
        922395600*pow3(lmMst2))*pow8(Mst1) + 296352000*lmMst2*(1 + lmMst2)*
        pow8(Mst2)) - 889056000*pow2(1 + lmMst2)*pow9(Mst2))) - 51450*Mt*pow2(
        Mst1)*pow2(s2t)*(Dmglst2*(9*(195293 + 639360*B4 - 8640*DN - 709200*
        lmMst1 + 145120*OepS2 - 108*(23989 + 27210*lmMst1 - 27210*lmMst2)*S2 +
        60480*pow2(lmMst1) + 720*lmMst2*(2149 - 1296*lmMst1 + 96*pow2(lmMst1))
        - 8640*(-101 + 105*lmMst1)*pow2(lmMst2) + 838080*pow3(lmMst2))*pow4(
        Mst1)*pow4(Mst2) - 2*pow2(Mst2)*(15069803 - 2877120*B4 + 38880*DN +
        6325200*lmMst1 - 1898720*OepS2 + 108*(525961 + 356010*lmMst1 - 356010*
        lmMst2)*S2 + 447120*pow2(lmMst1) - 360*lmMst2*(28667 - 5238*lmMst1 +
        3024*pow2(lmMst1)) + 38880*(-60 + 157*lmMst1)*pow2(lmMst2) - 155520*
        pow3(lmMst1) - 4860000*pow3(lmMst2))*pow6(Mst1) + 583200*Dmsqst2*(40*(1
        - lmMst1 + lmMst2)*pow2(Mst2)*pow4(Mst1) - 4*(2 + 5*lmMst1 - 5*lmMst2)*
        pow2(Mst1)*pow4(Mst2) + (80 - 43*lmMst1 + 43*lmMst2)*pow6(Mst1)) + 27*
        pow2(Mst1)*(46957 + 188640*B4 - 3600*DN - 37440*lmMst1 + 5600*OepS2 -
        324*(-453 + 350*lmMst1 - 350*lmMst2)*S2 + 11520*pow2(lmMst1) - 2880*
        lmMst2*(-221 + 55*lmMst1 + 4*pow2(lmMst1)) - 1440*(-232 + 121*lmMst1)*
        pow2(lmMst2) + 185760*pow3(lmMst2))*pow6(Mst2) + 622080*(-1 + 2*lmMst2
        + 3*pow2(lmMst2))*pow8(Mst2)) + 15*(-38880*Dmsqst2*Mst2*pow2(Mst1)*(8*(
        -lmMst1 + lmMst2)*pow2(Mst1)*pow2(Mst2) + (4 - 9*lmMst1 + 9*lmMst2)*
        pow4(Mst1) - 4*(1 + lmMst1 - lmMst2)*pow4(Mst2)) + 3*(346405 + 62208*B4
        + 246672*lmMst2 - 6688*OepS2 - 324*(-685 + 418*lmMst2)*S2 + 1728*(-7 +
        8*lmMst2)*pow2(lmMst1) + 323136*pow2(lmMst2) - 216*lmMst1*(990 + 1440*
        lmMst2 - 627*S2 + 584*pow2(lmMst2)) + 112320*pow3(lmMst2))*pow4(Mst1)*
        pow5(Mst2) + 2*(795601 + 93312*B4 + 365040*lmMst2 - 17056*OepS2 +
        663444*S2 - 345384*lmMst2*S2 + 432*(163 + 264*lmMst2)*pow2(lmMst1) +
        592704*pow2(lmMst2) - 216*lmMst1*(1609 + 3070*lmMst2 - 1599*S2 + 1692*
        pow2(lmMst2)) - 3456*pow3(lmMst1) + 254880*pow3(lmMst2))*pow3(Mst2)*
        pow6(Mst1) + 27*pow2(Mst1)*(27147 + 5280*B4 - 48*DN - 5952*lmMst1 -
        224*OepS2 + 324*(-1 + 14*lmMst1 - 14*lmMst2)*S2 - 768*pow2(lmMst1) -
        192*lmMst2*(-198 + 73*lmMst1 + 4*pow2(lmMst1)) - 288*(-84 + 19*lmMst1)*
        pow2(lmMst2) + 6240*pow3(lmMst2))*pow7(Mst2) + 41472*pow2(1 + lmMst2)*
        pow9(Mst2))) - Mst2*pow3(s2t)*(120*Dmsqst2*Mst2*pow2(Mst1)*(-12348*(
        97294 - 10000*OepS2 + 398250*S2 - 45*lmMst2*(383 + 4500*S2) + 45*
        lmMst1*(233 + 330*lmMst2 + 4500*S2) - 7425*(pow2(lmMst1) + pow2(lmMst2)
        ))*pow2(Mst2)*pow4(Mst1) - 1157625*(1177 - 32*OepS2 - 4860*S2 - 24*
        lmMst2*(14 + 27*S2) + 24*lmMst1*(14 - 6*lmMst2 + 27*S2) + 144*pow2(
        lmMst2))*pow2(Mst1)*pow4(Mst2) - 2*(222209689 - 71344000*OepS2 + 10080*
        lmMst2*(28298 - 143325*S2) + 3213567000*S2 + 1260*lmMst1*(-226384 +
        172935*lmMst2 + 1146600*S2) - 108949050*(pow2(lmMst1) + pow2(lmMst2)))*
        pow6(Mst1) + 41674500*(1 - 2*lmMst2)*pow6(Mst2)) + 3*pow3(Mst2)*(-
        1543500*(21005 + 1440*B4 - 144*D3 + 72*DN + 21216*lmMst1 - 2208*OepS2 +
        972*(307 + 46*lmMst1 - 46*lmMst2)*S2 + 1908*pow2(lmMst1) - 12*lmMst2*(
        2950 - 2040*lmMst1 + 117*pow2(lmMst1)) + 108*(-283 + 98*lmMst1)*pow2(
        lmMst2) + 576*pow3(lmMst1) - 9756*pow3(lmMst2))*pow4(Mst1)*pow4(Mst2) -
        41160*pow2(Mst2)*(4627751 + 48600*B4 + 5400*D3 - 5400*DN + 1320240*
        lmMst1 - 340000*OepS2 + 81000*(424 + 85*lmMst1 - 85*lmMst2)*S2 -
        662400*pow2(lmMst1) - 30*lmMst2*(12148 - 76560*lmMst1 + 12105*pow2(
        lmMst1)) + 2250*(-583 + 438*lmMst1)*pow2(lmMst2) - 49500*pow3(lmMst1) -
        572850*pow3(lmMst2))*pow6(Mst1) + 55566000*(119 + 218*lmMst2 + 32*
        lmMst1*(1 + lmMst2) + 107*pow2(lmMst2))*pow2(Mst1)*pow6(Mst2) - (
        96969711071 - 18919880000*OepS2 + 1260*lmMst2*(58499851 - 304069500*S2)
        + 1725694740000*S2 - 1058400*(76483 + 26985*lmMst2)*pow2(lmMst1) -
        78893665200*pow2(lmMst2) + 1260*lmMst1*(-61388401 + 126859740*lmMst2 +
        304069500*S2 + 48421800*pow2(lmMst2)) - 1296540000*pow3(lmMst1) -
        31154004000*pow3(lmMst2))*pow8(Mst1) - 889056000*pow2(1 + lmMst2)*pow8(
        Mst2)) - 4*Dmglst2*(-6174*(5430043 - 948000*OepS2 + 60*lmMst2*(8743 -
        319950*S2) - 218605500*S2 + 900*(-859 + 3690*lmMst2)*pow2(lmMst1) +
        1454400*pow2(lmMst2) - 60*lmMst1*(51493 + 24630*lmMst2 - 319950*S2 +
        112950*pow2(lmMst2)) + 45000*pow3(lmMst1) + 3411000*pow3(lmMst2))*pow4(
        Mst2)*pow6(Mst1) - 2315250*(14987 - 432*B4 + 576*D3 - 360*DN + 4752*
        lmMst1 - 224*OepS2 + 324*(-1677 + 14*lmMst1 - 14*lmMst2)*S2 - 1404*
        pow2(lmMst1) + 144*lmMst2*(-81 - 124*lmMst1 + 8*pow2(lmMst1)) - 36*(-
        439 + 246*lmMst1)*pow2(lmMst2) + 7704*pow3(lmMst2))*pow4(Mst1)*pow6(
        Mst2) - 5*pow2(Mst2)*(30586096049 - 3475001600*OepS2 + 1260*lmMst2*(
        2171669 - 55848240*S2) - 174295724400*S2 + 793800*(433 + 6552*lmMst2)*
        pow2(lmMst1) + 1577280600*pow2(lmMst2) - 1260*lmMst1*(2740559 +
        1524600*lmMst2 - 55848240*S2 + 7514640*pow2(lmMst2)) - 311169600*pow3(
        lmMst1) + 4578638400*pow3(lmMst2))*pow8(Mst1) + 138915000*Dmsqst2*(-36*
        (5 + 2*lmMst1 - 2*lmMst2)*pow4(Mst1)*pow4(Mst2) + 9*(55 - 34*lmMst1 +
        34*lmMst2)*pow2(Mst2)*pow6(Mst1) - 108*pow2(Mst1)*pow6(Mst2) + (419 -
        57*lmMst1 + 57*lmMst2)*pow8(Mst1)) - 83349000*(-20 + (230 + 32*lmMst1)*
        lmMst2 + 155*pow2(lmMst2))*pow2(Mst1)*pow8(Mst2) + 2667168000*lmMst2*(1
        + lmMst2)*power10(Mst2)))))/(1.500282e9*Tbeta*pow4(Mst1)))/pow6(Mst2)));
}

/**
 * @return The off-diagonal (1, 2) = (2, 1) matrix element of the Higgs mass matrix as a double for the hierarchy 'H6b'
 */
double H6b::getS12() const {
   return MuSUSY*((Mt*oneLoopFlag*s2t*(4*(lmMst1 - lmMst2)*pow2(Mt) - ((2 + lmMst1
        - lmMst2)*pow2(Mst1) + (-2 + lmMst1 - lmMst2)*pow2(Mst2))*pow2(s2t) - (
        2*Mt*MuSUSY*s2t*(2 - lmMst1 + lmMst2 - (2*(lmMst1 - lmMst2)*pow2(Mst1)*
        (pow2(Mst1) + pow2(Mst2)))/pow4(Mst2)))/Tbeta))/16. + (Al4p*Mt*
        twoLoopFlag*(12*Mt*pow2(s2t)*(2*(Mst2*(4 + 3*lmMst2 - lmMst1*(1 +
        lmMst2) + pow2(lmMst2)) + ((Dmglst2*(5 + lmMst1*(3 - 2*lmMst2) - 3*
        lmMst2 + 2*pow2(lmMst2)) + Mst2*(1 + 3*lmMst2 - lmMst1*(3 + 2*lmMst2) +
        2*pow2(lmMst2)))*pow2(Mst1))/pow2(Mst2)) + ((0.5 + 7*lmMst2 - lmMst1*(7
        + 4*lmMst2) + 4*pow2(lmMst2))*pow4(Mst1))/pow3(Mst2) + Dmglst2*(12 - 2*
        lmMst1*lmMst2 + 2*(lmMst1 + lmMst2 + pow2(lmMst2)) + ((10.5 + lmMst1*(9
        - 4*lmMst2) - 9*lmMst2 + 4*pow2(lmMst2))*pow4(Mst1))/pow4(Mst2))) + 4*
        s2t*pow2(Mt)*(12*lmMst1 + 4*(-2 + lmMst1)*lmMst2 + 2*pow2(lmMst1) - 6*
        pow2(lmMst2) - (4*(1 + lmMst2)*pow2(Mst2))/pow2(Mst1) + (Dmglst2*(4 +
        lmMst2*(8 - (8*pow2(Mst2))/pow2(Mst1))))/Mst2 + (4*(1 - lmMst1 +
        lmMst2)*(2*Dmglst2 + (lmMst1 - lmMst2)*Mst2)*pow2(Mst1))/pow3(Mst2) - (
        4*(Dmglst2*(-2 + 4*lmMst1 - 4*lmMst2) + Mst2*(-lmMst1 + lmMst2 - 2*
        lmMst1*lmMst2 + pow2(lmMst1) + pow2(lmMst2)))*pow4(Mst1))/pow5(Mst2)) -
        (16*pow3(Mt)*(Dmglst2*(2*(-7 + lmMst2*(-7 + lmMt) - lmMst1*(-5 + lmMt)
        + 2*lmMt)*pow2(Mst1)*pow2(Mst2) + (-11 + lmMst1*(32 + 6*lmMst2 - 8*
        lmMt) + 8*lmMst2*(-5 + lmMt) + 8*lmMt - 6*pow2(lmMst2))*pow4(Mst1) - 2*
        (5 + lmMst1*(-1 + lmMst2) + lmMst2 - pow2(lmMst2))*pow4(Mst2)) - Mst2*(
        2*(4 + lmMst1*(3 + 2*lmMst2 - lmMt) + lmMst2*(-4 + lmMt) + lmMt - 2*
        pow2(lmMst2))*pow2(Mst1)*pow2(Mst2) + (13 + 2*lmMst1*(6 + 3*lmMst2 - 2*
        lmMt) + 2*lmMt + 2*lmMst2*(-7 + 2*lmMt) - 6*pow2(lmMst2))*pow4(Mst1) +
        2*(2 + lmMst1 - 2*lmMst2 + lmMst1*lmMst2 + lmMt - pow2(lmMst2))*pow4(
        Mst2))))/pow6(Mst2) + ((pow3(s2t)*(2*(-16 + 6*lmMst2 - 2*lmMst1*(8 + 5*
        lmMst2) + 3*pow2(lmMst1) + 7*pow2(lmMst2))*pow3(Mst2)*pow4(Mst1) - 2*(-
        12 - 18*lmMst2 + 2*lmMst1*(5 + 3*lmMst2) + pow2(lmMst1) - 7*pow2(
        lmMst2))*pow2(Mst1)*pow5(Mst2) + Mst2*(3 + lmMst1*(2 - 32*lmMst2) - 2*
        lmMst2 + 16*pow2(lmMst1) + 16*pow2(lmMst2))*pow6(Mst1) - 4*Dmglst2*(2*(
        1 + 2*lmMst1 - lmMst2)*pow2(Mst2)*pow4(Mst1) + 2*(1 + lmMst1 - lmMst2 +
        2*lmMst1*lmMst2 - 2*pow2(lmMst2))*pow2(Mst1)*pow4(Mst2) + (1 + 2*lmMst1
        - 2*lmMst2)*pow6(Mst1) - 2*lmMst2*pow6(Mst2)) + 4*(1 + lmMst2)*pow7(
        Mst2)))/pow3(Mst2) - (2*Mt*MuSUSY*s2t*(4*(4*Mt*(5 + 6*lmMst2 - lmMst1*(
        4 + 3*lmMst2) + 3*pow2(lmMst2)) + Mst2*s2t*(-1 + 13*lmMst2 - lmMst1*(13
        + 8*lmMst2) + pow2(lmMst1) + 7*pow2(lmMst2)))*pow3(Mst2)*pow4(Mst1) +
        2*(-(Mst2*s2t*(-14 - 20*lmMst2 + 2*lmMst1*(5 + 3*lmMst2) + pow2(lmMst1)
        - 7*pow2(lmMst2))) + 8*Mt*(4 + 3*lmMst2 - lmMst1*(1 + lmMst2) + pow2(
        lmMst2)))*pow2(Mst1)*pow5(Mst2) + 4*Dmglst2*(Mst2*s2t*(-5 + 8*lmMst2 -
        4*lmMst1*(2 + lmMst2) + 4*pow2(lmMst2)) + Mt*(65 + lmMst1*(34 - 20*
        lmMst2) - 26*lmMst2 + 20*pow2(lmMst2)))*pow6(Mst1) + Mst2*(Mst2*s2t*(-1
        + 50*lmMst2 - 2*lmMst1*(25 + 32*lmMst2) + 20*pow2(lmMst1) + 44*pow2(
        lmMst2)) + Mt*(84 + 152*lmMst2 - 40*lmMst1*(3 + 2*lmMst2) + 80*pow2(
        lmMst2)))*pow6(Mst1) + 8*Dmglst2*((Mst2*s2t*(-2 + 3*lmMst2 - lmMst1*(3
        + 2*lmMst2) + 2*pow2(lmMst2)) + Mt*(22 + lmMst1*(8 - 6*lmMst2) - 4*
        lmMst2 + 6*pow2(lmMst2)))*pow2(Mst2)*pow4(Mst1) + (-(Mst2*s2t*(1 +
        lmMst1 - 2*lmMst2 + 2*lmMst1*lmMst2 - 2*pow2(lmMst2))) + 2*Mt*(6 +
        lmMst1 + lmMst2 - lmMst1*lmMst2 + pow2(lmMst2)))*pow2(Mst1)*pow4(Mst2)
        + lmMst2*s2t*pow7(Mst2)) + 4*(1 + lmMst2)*s2t*pow8(Mst2)))/(Tbeta*pow6(
        Mst2)))/pow2(Mst1)))/12. + Al4p*xDmglst2*pow2(Dmglst2)*((Mt*
        twoLoopFlag*(-2*(pow2(Mst2)*(36*(5 - 3*lmMst1 + 3*lmMst2)*Mst2*s2t*
        pow2(Mt) - 54*(12 + lmMst1 - lmMst2)*Mt*pow2(Mst2)*pow2(s2t) + 4*(155 +
        3*lmMst2*(41 - 6*lmMt) - 18*lmMst1*(4 + lmMst2 - lmMt) - 51*lmMt + 18*
        pow2(lmMst2))*pow3(Mt) - 9*(4*lmMst1 - 5*lmMst2)*pow3(Mst2)*pow3(s2t))*
        pow4(Mst1) + pow2(Mst1)*(6*(17 - 6*lmMst2)*Mst2*s2t*pow2(Mt) - 54*(8 +
        lmMst1 - lmMst2)*Mt*pow2(Mst2)*pow2(s2t) + 4*(11 + 18*lmMst1 - 21*
        lmMst2 + 3*lmMt)*pow3(Mt) + 9*(6 + 5*lmMst2 + lmMst1*(-5 + 2*lmMst2) -
        2*pow2(lmMst2))*pow3(Mst2)*pow3(s2t))*pow4(Mst2)) - (72*(9 - 10*lmMst1
        + 10*lmMst2)*Mst2*s2t*pow2(Mt) + 54*(-27 + 4*lmMst1 - 4*lmMst2)*Mt*
        pow2(Mst2)*pow2(s2t) + 8*(335 + 699*lmMst2 - 18*lmMst1*(29 + 7*lmMst2 -
        7*lmMt) - 3*(59 + 42*lmMst2)*lmMt + 126*pow2(lmMst2))*pow3(Mt) + 9*(1 -
        10*lmMst1 + 10*lmMst2)*pow3(Mst2)*pow3(s2t))*pow6(Mst1) + 18*(-2 +
        lmMst2)*s2t*(-4*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow7(Mst2) + (18*Mt*
        MuSUSY*s2t*(2*(-8*(10 + lmMst1 - lmMst2)*Mt + Mst2*s2t*(8 - 9*lmMst1 +
        9*lmMst2 + 2*lmMst1*lmMst2 - 2*pow2(lmMst2)))*pow2(Mst2)*pow4(Mst1) +
        2*(-4*(8 + lmMst1 - lmMst2)*Mt + Mst2*s2t*(8 - 5*lmMst1 + 4*lmMst2 + 2*
        lmMst1*lmMst2 - 2*pow2(lmMst2)))*pow2(Mst1)*pow4(Mst2) + (-268*Mt +
        Mst2*s2t*(17 + 4*lmMst1*(-7 + lmMst2) + 28*lmMst2 - 4*pow2(lmMst2)))*
        pow6(Mst1) - 2*(-2 + lmMst2)*s2t*pow7(Mst2)))/Tbeta))/(54.*pow2(Mst1)*
        pow7(Mst2)) + Al4p*threeLoopFlag*(s2t*pow3(Mt)*((2*(652 - 32*lmMst1*(-2
        + lmMst2) + 70*lmMst2 - 59*pow2(lmMst2)))/(9.*pow2(Mst1)) + (
        1288.4830687830688 - (16*B4)/9. + (16*D3)/9. + (184*lmMst1)/3. - (32*
        lmMt)/9. + (2*lmMst2*(6695 + 1002*lmMst1 - 96*pow2(lmMst1)))/27. + (82*
        pow2(lmMst1))/9. - (8*(DN + (47 + 8*lmMst1)*pow2(lmMst2)))/9. + (128*
        pow3(lmMst2))/9.)/pow2(Mst2) - (64*(2 + lmMst2 - 3*pow2(lmMst2))*pow2(
        Mst2))/(9.*pow4(Mst1)) - S2*((89.65714285714286 + 84*lmMst1 - 84*
        lmMst2)/pow2(Mst2) + ((847.3523809523809 + 480*lmMst1 - 480*lmMst2)*
        pow2(Mst1))/pow4(Mst2)) + ((16*OepS2*(40*pow2(Mst1) + 7*pow2(Mst2)))/
        27. + pow2(Mst1)*(1803.9172204585539 - (16*B4)/9. + (16*D3)/9. - (8*DN)
        /9. + (225124*lmMst1)/675. + (32*(17 - 5*lmMst1 + 5*lmMst2)*lmMt)/9. +
        (4*lmMst2*(136544 + 53535*lmMst1 - 1425*pow2(lmMst1)))/675. + (224*
        pow2(lmMst1))/15. - (4*(3257 + 545*lmMst1)*pow2(lmMst2))/45. + (4*pow3(
        lmMst1))/9. + (508*pow3(lmMst2))/9.))/pow4(Mst2)) + Mt*(pow3(s2t)*(
        46.46594650205761 - (38*B4)/9. + (40*D3)/9. - (7*DN)/3. + (98*lmMst1)/
        3. - (25*pow2(lmMst1))/18. + (lmMst2*(-1493 - 24*lmMst1 + 48*pow2(
        lmMst1)))/27. + (4*OepS2*(21 + (1157*pow2(Mst1))/pow2(Mst2)))/729. + (-
        ((247 + 246*lmMst1)*pow2(lmMst2)) + ((-876 + 32*lmMst1*(-2 + lmMst2) -
        70*lmMst2 + 347*pow2(lmMst2))*pow2(Mst2))/pow2(Mst1))/18. - ((S2*((
        648463 + 34710*lmMst1 - 34710*lmMst2)*pow2(Mst1) + 9*(65819 + 70*lmMst1
        - 70*lmMst2)*pow2(Mst2)))/270. + pow2(Mst1)*(12.389529492455418 + (32*
        B4)/9. - (32*D3)/9. + (16*DN)/9. + (138661*lmMst1)/4050. + (159*pow2(
        lmMst1))/10. - lmMst2*(26.7558024691358 + (229*lmMst1)/5. + (143*pow2(
        lmMst1))/9.) + ((1333 + 1355*lmMst1)*pow2(lmMst2))/45. + (5*pow3(
        lmMst1))/9. - (133*pow3(lmMst2))/9.))/pow2(Mst2) + (107*pow3(lmMst2))/
        9. + (16*(2 + lmMst2 - 3*pow2(lmMst2))*pow4(Mst2))/(9.*pow4(Mst1))) - (
        2*T1ep*(21*pow2(Mst2)*(36*Mst2*s2t*pow2(Mt) - 24*Mt*pow2(Mst2)*pow2(
        s2t) + 32*pow3(Mt) + pow3(Mst2)*pow3(s2t)) + pow2(Mst1)*(4320*Mst2*s2t*
        pow2(Mt) - 8715*Mt*pow2(Mst2)*pow2(s2t) + 17584*pow3(Mt) + 1157*pow3(
        Mst2)*pow3(s2t))))/(243.*pow5(Mst2))) + (2*(-3*(1648637 + 1173060*lmMt
        - 15680*OepS2 + 1680*lmMst2*(-1214 + 42*lmMt - 189*S2) - 953208*S2 +
        7560*lmMst1*(131 + 33*lmMst2 - 8*lmMt + 42*S2) + 30240*pow2(lmMst1) -
        304920*pow2(lmMst2) + 15120*pow2(lmMt))*pow2(Mst1)*pow2(Mst2) + 2*(-
        1657412 - 7512750*lmMt + 615440*OepS2 + 4340358*S2 + 11340*(35 + 4*
        lmMst2 - 4*lmMt)*pow2(lmMst1) + 945*(-1981 + 48*lmMt)*pow2(lmMst2) +
        315*lmMst1*(-4223 + 2679*lmMst2 + 201*lmMt - 39564*S2 + 432*pow2(
        lmMst2) - 432*pow2(lmMt)) + 385560*pow2(lmMt) + 945*lmMst2*(5321 + 193*
        lmMt + 13188*S2 + 144*pow2(lmMt)) - 181440*pow3(lmMst2))*pow4(Mst1) -
        181440*(-3 - 4*lmMst2 + 3*pow2(lmMst2))*pow4(Mst2))*pow4(Mt))/(25515.*
        pow2(Mst1)*pow5(Mst2)) + pow2(Mt)*(-(pow2(s2t)*((32*Mst2*(3 + 4*lmMst2
        - 3*pow2(lmMst2)))/(3.*pow2(Mst1)) + (380.51604938271606 - 76*B4 + 2*DN
        - 292*lmMst1 + (2*(631 - 438*lmMst1)*lmMst2)/9. - (16*pow2(lmMst1))/3.
         + (4*(53 + 48*lmMst1)*pow2(lmMst2))/
        3. - 64*pow3(lmMst2))/Mst2 + ((28*OepS2*(415*pow2(Mst1) + 24*pow2(Mst2)
        ))/243. - (S2*((21422 + 87150*lmMst1 - 87150*lmMst2)*pow2(Mst1) + 9*(-
        1373 + 560*lmMst1 - 560*lmMst2)*pow2(Mst2)))/90. + pow2(Mst1)*(
        99.93816872427983 - 76*B4 + 2*DN - (6424*lmMst1)/27. - (4*(293 + 1152*
        lmMst1)*lmMst2)/27. + (92*pow2(lmMst1))/3. + 4*(35 + 16*lmMst1)*pow2(
        lmMst2) - 64*pow3(lmMst2)))/pow3(Mst2))) - (MuSUSY*(75*pow2(Mst2)*(-8*
        Mst2*Mt*s2t*(334138 - 61560*B4 + 1620*DN - 236520*lmMst1 - 180*(-823 +
        438*lmMst1)*lmMst2 + 2240*OepS2 - 81*(-1373 + 560*lmMst1 - 560*lmMst2)*
        S2 - 3360*T1ep - 4320*pow2(lmMst1) + 1080*(29 + 48*lmMst1)*pow2(lmMst2)
        - 51840*pow3(lmMst2)) + pow2(Mst2)*pow2(s2t)*(13169 - 41040*B4 + 43200*
        D3 - 22680*DN + 282960*lmMst1 + 1120*OepS2 - 324*(65819 + 70*lmMst1 -
        70*lmMst2)*S2 - 1680*T1ep - 13500*pow2(lmMst1) + 720*lmMst2*(-775 + 12*
        lmMst1 + 24*pow2(lmMst1)) - 1080*(-2 + 123*lmMst1)*pow2(lmMst2) +
        115560*pow3(lmMst2)) + 96*pow2(Mt)*(33934 - 90*B4 + 90*D3 - 45*DN +
        120*(163 + 24*lmMst1)*lmMst2 - 120*lmMt - 10206*S2 - 720*(2 + lmMst1)*
        pow2(lmMst2) + 720*(lmMst1 + pow3(lmMst2))))*pow4(Mst1) - 40500*s2t*(
        Mst2*s2t*(812 - 32*lmMst1*(-2 + lmMst2) + 38*lmMst2 - 251*pow2(lmMst2))
        + 128*Mt*(3 + 4*lmMst2 - 3*pow2(lmMst2)))*pow2(Mst1)*pow5(Mst2) + 2*(-
        2*pow2(Mst2)*pow2(s2t)*(2011073 + 1417500*B4 - 1458000*D3 + 749250*DN +
        934245*lmMst1 - 1178000*OepS2 + 1350*(620417 + 17670*lmMst1 - 17670*
        lmMst2)*S2 + 1767000*T1ep + 3150900*pow2(lmMst1) - 45*lmMst2*(-124139 +
        189090*lmMst1 + 71550*pow2(lmMst1)) + 4050*(1323 + 1970*lmMst1)*pow2(
        lmMst2) + 101250*pow3(lmMst1) - 4860000*pow3(lmMst2)) - 125*Mst2*Mt*
        s2t*(996211 - 295488*B4 + 7776*DN - 1030176*lmMst1 - 144*(-1883 + 3618*
        lmMst1)*lmMst2 + 98336*OepS2 - 756*(259 + 2634*lmMst1 - 2634*lmMst2)*S2
        - 147504*T1ep + 49248*pow2(lmMst1) + 5184*(67 + 48*lmMst1)*pow2(lmMst2)
        - 248832*pow3(lmMst2)) + 150*pow2(Mt)*(2199511 - 4320*B4 + 4320*D3 -
        2160*DN + 140160*lmMst1 + 960*(1426 + 303*lmMst1)*lmMst2 - 2880*(-16 +
        5*lmMst1 - 5*lmMst2)*lmMt - 1120*OepS2 + 324*(-3411 + 70*lmMst1 - 70*
        lmMst2)*S2 + 1680*T1ep - 2880*(77 + 24*lmMst1)*pow2(lmMst2) + 69120*
        pow3(lmMst2)))*pow6(Mst1) + 1296000*(2 + lmMst2 - 3*pow2(lmMst2))*pow2(
        s2t)*pow8(Mst2)))/(364500.*Tbeta*pow4(Mst1)*pow6(Mst2))))) - (Al4p*Mt*
        xDR2DRMOD*(54*s2t*pow2(Mst1)*(-5*Al4p*threeLoopFlag*xDmsqst2*pow2(
        Dmsqst2)*((-2*lmMst1 + 2*lmMst2)*s2t*pow2(Mst1)*(2*Mt*MuSUSY*pow2(Mst1)
        + (14*Dmglst2*Mst2 + 19*xDmglst2*pow2(Dmglst2) + pow2(Mst2))*(2*Mt*
        MuSUSY - s2t*Tbeta*pow2(Mst2))) + (14*Dmglst2*Mst2 + 19*xDmglst2*pow2(
        Dmglst2) + pow2(Mst2))*(2*Mt*(MuSUSY*s2t - 2*Mt*Tbeta)*pow2(Mst1) +
        pow2(Mst2)*(2*Mt*MuSUSY*s2t + 4*Tbeta*pow2(Mt) - Tbeta*pow2(Mst2)*pow2(
        s2t)) + Tbeta*pow2(s2t)*pow4(Mst1))) + 4*Mst2*(2*Dmglst2*lmMst2 + Mst2
        + lmMst2*Mst2)*twoLoopFlag*(-((2*Mt*(MuSUSY*s2t - 2*Mt*Tbeta)*pow2(
        Mst1) + pow2(Mst2)*(2*Mt*MuSUSY*s2t + 4*Tbeta*pow2(Mt) - Tbeta*pow2(
        Mst2)*pow2(s2t)) + Tbeta*pow2(s2t)*pow4(Mst1))*pow4(Mst2)) + 2*(lmMst1
        - lmMst2)*s2t*pow2(Mst1)*(2*Mt*MuSUSY*(pow2(Mst1)*pow2(Mst2) + pow4(
        Mst1) + pow4(Mst2)) - s2t*Tbeta*pow6(Mst2)))) + xDmglst2*pow2(Dmglst2)*
        (216*(2 - lmMst2)*s2t*twoLoopFlag*pow2(Mst1)*(Tbeta*(-4*pow2(Mst1)*
        pow2(Mt) + 4*pow2(Mst2)*pow2(Mt) + 2*(lmMst1 - lmMst2)*pow2(Mst1)*pow2(
        Mst2)*pow2(s2t) + pow2(s2t)*(pow4(Mst1) - pow4(Mst2)))*pow4(Mst2) + 2*
        Mt*MuSUSY*s2t*((pow2(Mst1) + pow2(Mst2))*pow4(Mst2) - 2*(lmMst1 -
        lmMst2)*pow2(Mst1)*(pow2(Mst1)*pow2(Mst2) + pow4(Mst1) + pow4(Mst2))))
        + Al4p*threeLoopFlag*(-8*Tbeta*pow2(Mst2)*(3510*Dmsqst2*s2t*pow2(Mt) -
        144*Mt*pow2(s2t)*(13 - 125*lmMst2 + 6*lmMst1*(7 + 8*lmMst2 - 6*pow2(
        lmMst2)) - 30*pow2(lmMst2) + 36*pow3(lmMst2))*pow3(Mst2) + 64*Mst2*(215
        + 37*lmMst2 + 3*(-5 + 8*lmMst2)*lmMt - 18*lmMst1*(1 - 2*lmMst2 + lmMt)
        - 42*pow2(lmMst2))*pow3(Mt) + 3*pow2(Mst2)*(2*s2t*pow2(Mt)*(964 + 306*
        lmMst2 + 96*(-2 + lmMst2)*pow2(lmMst1) - 81*pow2(lmMst2) + 96*(lmMst1 +
        3*lmMst1*lmMst2 - 2*lmMst1*pow2(lmMst2) + pow3(lmMst2))) - 585*Dmsqst2*
        (lmMst1 - lmMst2)*pow3(s2t) + pow2(Mst2)*(200 + 1196*lmMst2 - 48*(-2 +
        lmMst2)*pow2(lmMst1) + 306*pow2(lmMst2) + 3*lmMst1*(-492 + 10*lmMst2 +
        123*pow2(lmMst2)) - 321*pow3(lmMst2))*pow3(s2t)))*pow4(Mst1) - 4*Mst2*
        Tbeta*(-1728*Mt*pow2(Mst2)*pow2(s2t)*(15 - 11*lmMst2 + lmMst1*(4 + 7*
        lmMst2 - 6*pow2(lmMst2)) - 7*pow2(lmMst2) + 6*pow3(lmMst2)) + 128*(1097
        + lmMst2*(1531 - 246*lmMt) - 339*lmMt + (444 - 36*lmMt)*pow2(lmMst2) -
        18*lmMst1*(39 - 2*lmMst2*(-9 + lmMt) - 7*lmMt + 2*pow2(lmMst2)) + 36*
        pow3(lmMst2))*pow3(Mt) + Mst2*(1152*s2t*pow2(Mt)*(26 + 29*lmMst2 + (-2
        + lmMst2)*pow2(lmMst1) + 4*pow2(lmMst2) - 2*lmMst1*(8 + lmMst2 + pow2(
        lmMst2)) + pow3(lmMst2)) - 3*(585*Dmsqst2 + pow2(Mst2)*(-476 - 1790*
        lmMst2 + 768*(-2 + lmMst2)*pow2(lmMst1) - 1617*pow2(lmMst2) - 96*
        lmMst1*(-13 - 33*lmMst2 + 16*pow2(lmMst2)) + 768*pow3(lmMst2)))*pow3(
        s2t)))*pow6(Mst1) - 24*s2t*(48*Tbeta*(2 + lmMst2 - 3*pow2(lmMst2))*(4*
        pow2(Mt) - pow2(Mst2)*pow2(s2t))*pow8(Mst2) + Mt*MuSUSY*(pow2(Mst2)*(
        585*Dmsqst2*(-1 + 2*lmMst1 - 2*lmMst2)*s2t - 64*Mst2*Mt*(5 + 149*lmMst2
        + 12*pow2(lmMst2) + 6*lmMst1*(-7 - 8*lmMst2 + 6*pow2(lmMst2)) - 36*
        pow3(lmMst2)) + s2t*pow2(Mst2)*(-1540 - 2506*lmMst2 + 96*(-2 + lmMst2)*
        pow2(lmMst1) + lmMst1*(2760 + 36*lmMst2 - 738*pow2(lmMst2)) + 141*pow2(
        lmMst2) + 642*pow3(lmMst2)))*pow4(Mst1) + 3*pow2(Mst1)*(-195*Dmsqst2*
        s2t + 128*Mst2*Mt*(-3 - 4*lmMst2 + 3*pow2(lmMst2)) + s2t*(-380 + 32*
        lmMst1*(-2 + lmMst2) - 38*lmMst2 + 251*pow2(lmMst2))*pow2(Mst2))*pow4(
        Mst2) + 2*(585*Dmsqst2*(lmMst1 - lmMst2)*s2t + 32*Mst2*Mt*(85 - 215*
        lmMst2 + lmMst1*(66 + 90*lmMst2 - 72*pow2(lmMst2)) - 54*pow2(lmMst2) +
        72*pow3(lmMst2)) + 3*s2t*pow2(Mst2)*(-336 - 716*lmMst2 + 144*(-2 +
        lmMst2)*pow2(lmMst1) + lmMst1*(668 + 534*lmMst2 - 379*pow2(lmMst2)) -
        246*pow2(lmMst2) + 235*pow3(lmMst2)))*pow6(Mst1) + 96*s2t*(2 + lmMst2 -
        3*pow2(lmMst2))*pow8(Mst2))) + 9*Tbeta*pow2(Mst1)*(3120*Dmsqst2*s2t*
        pow2(Mt)*pow4(Mst2) + 4*(4*s2t*(220 - 32*lmMst1*(-2 + lmMst2) + 70*
        lmMst2 - 59*pow2(lmMst2))*pow2(Mt) - 195*Dmsqst2*pow3(s2t))*pow6(Mst2)
        + 256*(-3 - 4*lmMst2 + 3*pow2(lmMst2))*(-4*pow3(Mt)*pow5(Mst2) + 3*Mt*
        pow2(s2t)*pow7(Mst2)) + pow3(s2t)*(45*pow4(Dmsqst2) + 4*(-444 + 32*
        lmMst1*(-2 + lmMst2) - 70*lmMst2 + 347*pow2(lmMst2))*pow8(Mst2)))))))/(
        648.*Tbeta*pow4(Mst1)*pow6(Mst2))) + (Mt*((xMst*(27*(lmMst1 - lmMst2)*
        Mt*oneLoopFlag*pow2(MuSUSY)*pow2(s2t)*pow3(Mst2) - Al4p*twoLoopFlag*(2*
        xDmglst2*pow2(Dmglst2)*(-(Mst2*Mt*pow2(s2t)*(18*(85 - 36*lmMst1 + 36*
        lmMst2)*Mst2*MuSUSY*Tbeta + (314 + 24*lmMst2*(29 - 6*xDR2DRMOD) - 24*
        lmMst1*(29 + 3*lmMst2*(-1 + xDR2DRMOD) - 6*xDR2DRMOD) + 72*(-1 +
        xDR2DRMOD)*pow2(lmMst2))*pow2(MuSUSY) + 15*(-43 + 60*lmMst1 - 60*
        lmMst2)*pow2(Mst2)*(-1 + pow2(Sbeta))*pow2(Sbeta))) + MuSUSY*(6*s2t*
        pow2(Mt)*(8*(143 - 18*lmMst1 + 18*lmMst2)*MuSUSY + Mst2*Tbeta*(371 +
        lmMst2*(552 - 600*pow2(Sbeta)) - 430*pow2(Sbeta) + 24*lmMst1*(-23 + 25*
        pow2(Sbeta)))) + Tbeta*(8*(398 + lmMst2*(2085 - 396*lmMt) - 18*lmMst1*(
        95 + 22*lmMst2 - 22*lmMt) - 375*lmMt + 396*pow2(lmMst2))*pow3(Mt) + 4*(
        1 - 24*lmMst1 + 24*lmMst2)*pow3(Mst2)*pow3(s2t)))) + Mst2*MuSUSY*(8*
        Dmglst2*(s2t*(36*(-1 + 3*lmMst1 - 3*lmMst2)*Mst2*Tbeta + MuSUSY*(785 +
        6*lmMst1*(85 - 42*lmMst2) - 438*lmMst2 + 252*pow2(lmMst2)))*pow2(Mt) -
        Mst2*Mt*(3*Mst2*Tbeta*(50 + lmMst1*(51 - 18*lmMst2) - 51*lmMst2 + 18*
        pow2(lmMst2)) + MuSUSY*(49 - 84*lmMst2 + lmMst1*(84 - 36*lmMst2*(-1 +
        xDR2DRMOD)) + 36*(-1 + xDR2DRMOD)*pow2(lmMst2)))*pow2(s2t) + Tbeta*(4*(
        7 - 381*lmMst2 + lmMst1*(327 + 72*lmMst2 - 81*lmMt) + 54*lmMt + 81*
        lmMst2*lmMt - 72*pow2(lmMst2))*pow3(Mt) + 2*(1 + 3*lmMst1 - 3*lmMst2)*
        pow3(Mst2)*pow3(s2t))) - Mst2*(-8*s2t*(18*Mst2*Tbeta*(-lmMst1 + lmMst2
        - 2*lmMst1*lmMst2 + pow2(lmMst1) + pow2(lmMst2)) + MuSUSY*(193 + 474*
        lmMst2 - 6*lmMst1*(67 + 42*lmMst2) + 252*pow2(lmMst2)))*pow2(Mt) + 8*
        Mst2*Mt*(3*Mst2*Tbeta*(1 + 33*lmMst2 - 3*lmMst1*(11 + 6*lmMst2) + 18*
        pow2(lmMst2)) + MuSUSY*(1 + 3*lmMst2*(-37 + 6*xDR2DRMOD) - 3*lmMst1*(-
        37 + 6*lmMst2*(-12 + xDR2DRMOD) + 6*xDR2DRMOD) - 81*pow2(lmMst1) + 9*(-
        15 + 2*xDR2DRMOD)*pow2(lmMst2)))*pow2(s2t) + Tbeta*(32*(83 + 3*lmMst1*(
        29 + 12*lmMst2 - 9*lmMt) + 9*lmMt + 3*lmMst2*(-32 + 9*lmMt) - 36*pow2(
        lmMst2))*pow3(Mt) + (5 + lmMst1*(6 - 288*lmMst2) - 6*lmMst2 + 144*(
        pow2(lmMst1) + pow2(lmMst2)))*pow3(Mst2)*pow3(s2t))))))*pow6(Mst1))/
        108. - (Al4p*z2*(MuSUSY*(6*xDmglst2*pow2(Dmglst2)*pow2(Mst2)*(-22680*
        twoLoopFlag*pow2(Mst1)*(-2*Mt*pow2(Mst1)*pow2(Mst2)*(-24*Mt*MuSUSY*s2t
        + 4*Tbeta*pow2(Mt) + Mst2*(MuSUSY + 9*Mst2*Tbeta)*pow2(s2t)) - 2*Mt*(-
        36*Mt*MuSUSY*s2t + 28*Tbeta*pow2(Mt) + Mst2*(MuSUSY + 9*Mst2*Tbeta)*
        pow2(s2t))*pow4(Mst1) + s2t*(-2*Mst2*Mt*s2t*(MuSUSY + 9*Mst2*Tbeta) +
        24*MuSUSY*pow2(Mt) + Tbeta*pow2(s2t)*pow3(Mst2))*pow4(Mst2)) + Al4p*
        Mst2*threeLoopFlag*(-3*pow2(Mst1)*pow2(Mst2)*(4*Mst2*s2t*(28*(46987 +
        15390*lmMst1 + 4050*lmMst2)*MuSUSY - 3*(255749 + 2100*lmMst1 + 441420*
        lmMst2)*Mst2*Tbeta)*pow2(Mt) - 14*Mt*((40001 + 52560*lmMst1 - 37080*
        lmMst2)*MuSUSY + 6*(46987 + 15390*lmMst1 + 4050*lmMst2)*Mst2*Tbeta)*
        pow2(Mst2)*pow2(s2t) + 32*(84*(17 + 1320*lmMst2)*MuSUSY + (-217444 +
        53865*lmMst1 - 116550*lmMst2 + 62685*lmMt)*Mst2*Tbeta)*pow3(Mt) + 7*(
        30641 + 52560*lmMst1 - 37080*lmMst2)*Tbeta*pow3(s2t)*pow4(Mst2)) +
        pow4(Mst1)*(-4*Mst2*s2t*(35*(-184517 + 74952*lmMst1 + 18360*lmMst2)*
        MuSUSY + 36*(-297316 + 20265*lmMst1 - 231945*lmMst2)*Mst2*Tbeta)*pow2(
        Mt) + 7*Mt*(4*(-261247 + 163890*lmMst1 - 140670*lmMst2)*MuSUSY + 3*(-
        1486429 + 190080*lmMst1 + 43200*lmMst2)*Mst2*Tbeta)*pow2(Mst2)*pow2(
        s2t) - 8*(21*(12761 + 5400*lmMst1 + 178920*lmMst2)*MuSUSY + 2*(-
        10736701 + 1553580*lmMst1 - 3303720*lmMst2 + 1508220*lmMt)*Mst2*Tbeta)*
        pow3(Mt) + 7*(642497 - 170100*lmMst1 + 170100*lmMst2)*Tbeta*pow3(s2t)*
        pow4(Mst2)) - 196560*s2t*(-2*Mt*MuSUSY*s2t - 4*Tbeta*pow2(Mt) + Tbeta*
        pow2(Mst2)*pow2(s2t))*pow6(Mst2))) - 68040*twoLoopFlag*pow2(Mst1)*(4*
        Mt*xMst*(xDmglst2*pow2(Dmglst2)*(48*Mt*MuSUSY*s2t - 88*Tbeta*pow2(Mt) -
        Mst2*(MuSUSY + 9*Mst2*Tbeta)*pow2(s2t)) + 2*Mst2*(Mst2*Mt*(MuSUSY*s2t -
        8*Mt*Tbeta) + Dmglst2*(17*Mt*MuSUSY*s2t + 16*Tbeta*pow2(Mt) - Mst2*(
        MuSUSY + 3*Mst2*Tbeta)*pow2(s2t))))*pow6(Mst1) + pow3(Mst2)*(2*Mt*
        MuSUSY*s2t*(4*Mt*((9*Dmglst2 + Mst2)*pow2(Mst1)*pow2(Mst2) + (13*
        Dmglst2 + Mst2)*pow4(Mst1) + 5*Dmglst2*pow4(Mst2) + pow5(Mst2)) - Mst2*
        s2t*(4*Dmglst2*(pow2(Mst1)*pow2(Mst2) + pow4(Mst1) + pow4(Mst2)) +
        pow5(Mst2))) + Tbeta*(16*pow3(Mt)*(-2*pow2(Mst1)*pow3(Mst2) + 3*(
        Dmglst2 - Mst2)*pow4(Mst1) - (Dmglst2 + Mst2)*pow4(Mst2)) - 6*Mt*pow2(
        Mst2)*pow2(s2t)*(Dmglst2*(4*pow2(Mst1)*pow2(Mst2) + 4*pow4(Mst1) + 5*
        pow4(Mst2)) + pow5(Mst2)) + (Mst2*(4*Dmglst2 + Mst2) - pow2(Mst1))*
        pow3(s2t)*pow6(Mst2))))) + 7*Al4p*threeLoopFlag*pow3(Mst2)*(300*
        xDmsqst2*pow2(Dmsqst2)*(3*s2t*(Mt*pow2(MuSUSY)*(6*Mst2*(-144*(3 +
        lmMst1 - lmMst2)*Mt + (209 + 18*lmMst1 - 18*lmMst2)*Mst2*s2t)*pow2(
        Mst1) + 8*(421 + 54*lmMst1 - 54*lmMst2)*s2t*pow4(Mst1) + 108*s2t*pow4(
        Mst2)) + 27*s2t*((51 - 6*lmMst1 + 6*lmMst2)*Mt + 2*(11 + 10*lmMst1 -
        10*lmMst2)*Mst2*s2t)*(-1 + pow2(Sbeta))*pow2(Sbeta)*pow6(Mst1)) -
        MuSUSY*Tbeta*(9*Mst2*pow2(Mst1)*(-418*Mst2*s2t*pow2(Mt) - 216*(3 +
        lmMst1 - lmMst2)*Mt*pow2(Mst2)*pow2(s2t) + 288*(-2 + lmMst1 - 2*lmMst2
        + lmMt)*pow3(Mt) + (191 + 18*lmMst1 - 18*lmMst2)*pow3(Mst2)*pow3(s2t))
        + 3*s2t*(-1296*(1 + lmMst1 - lmMst2)*Mst2*Mt*s2t - 2696*pow2(Mt) + (
        1057 + 162*lmMst1 - 162*lmMst2)*pow2(Mst2)*pow2(s2t))*pow4(Mst1) +
        4180*pow3(s2t)*pow6(Mst1) + 162*(-4*s2t*pow2(Mt)*pow4(Mst2) + pow3(s2t)
        *pow6(Mst2)))) + 2*Mt*(2*Mt*pow2(Mst1)*pow2(MuSUSY)*(-27*pow2(Mst2)*(5*
        Mst2*(64*(53 + 24*lmMst2)*Mst2*Mt + 2880*Dmsqst2*(2 + lmMst1 - lmMst2)*
        s2t + 3*(2269 + 664*lmMst1 - 56*lmMst2)*s2t*pow2(Mst2)) + Dmglst2*(
        7680*(5 + 6*lmMst2)*Mst2*Mt - 14400*Dmsqst2*(2 + 3*lmMst1 - 3*lmMst2)*
        s2t + (173947 - 25080*lmMst1 + 68760*lmMst2)*s2t*pow2(Mst2))) - 36*
        pow2(Mst1)*(5*Mst2*(3*(-353 + 72*lmMst1 + 696*lmMst2)*Mst2*Mt + 2160*
        Dmsqst2*(5 + 3*lmMst1 - 3*lmMst2)*s2t + 2*(3937 + 2988*lmMst1 - 2232*
        lmMst2)*s2t*pow2(Mst2)) + Dmglst2*(30*(3643 - 120*lmMst1 + 3192*lmMst2)
        *Mst2*Mt + 10800*Dmsqst2*(1 - 9*lmMst1 + 9*lmMst2)*s2t + (364291 -
        88560*lmMst1 + 147960*lmMst2)*s2t*pow2(Mst2))) + 2*(45*(38401 + 1080*
        lmMst1 - 7992*lmMst2)*Mt + Dmglst2*(-15057833 + 4014360*lmMst1 -
        5563080*lmMst2)*s2t - 3*(266863 + 396360*lmMst1 - 346680*lmMst2)*Mst2*
        s2t)*pow4(Mst1)) + pow2(s2t)*(300*Dmsqst2*(pow2(MuSUSY)*(-27*(72*
        Dmglst2*(2 + lmMst1 - lmMst2) + Mst2*(53 + 12*(-1 + 2*lmMst1 - 2*
        lmMst2)*shiftst1))*pow2(Mst1)*pow3(Mst2) - 9*Mst2*(72*Dmglst2*(3 + 8*
        lmMst1 - 8*lmMst2) + Mst2*(1097 + 72*(lmMst1 - lmMst2)*(-1 + shiftst1))
        )*pow4(Mst1) + 324*(-2*Dmglst2 + Mst2 + Mst2*shiftst1)*pow5(Mst2)) + (-
        ((20701 + 648*(lmMst1 - lmMst2)*(-1 + shiftst1))*pow2(MuSUSY)) + 162*
        Dmglst2*(-93 + 10*lmMst1 - 10*lmMst2)*Mst2*(-1 + pow2(Sbeta))*pow2(
        Sbeta))*pow6(Mst1)) + Mst2*pow2(MuSUSY)*(4*Dmglst2*(162*(6803 - 1810*
        lmMst1 + 2670*lmMst2)*pow2(Mst2)*pow4(Mst1) - 135*(648*lmMst1 - 7*(233
        + 240*lmMst2))*pow2(Mst1)*pow4(Mst2) + (4256978 - 615600*lmMst1 +
        754920*lmMst2)*pow6(Mst1) - 22680*pow6(Mst2)) + 3*Mst2*(-6*(533629 -
        1080*shiftst3 + 180*lmMst1*(-5 + 60*shiftst1 + 12*shiftst3) - 180*
        lmMst2*(-1 + 60*shiftst1 + 12*shiftst3))*pow2(Mst2)*pow4(Mst1) - 90*(
        5597 - 360*shiftst1 - 72*shiftst3 - 24*lmMst2*(46 + 30*shiftst1 + 3*
        shiftst3) + 12*lmMst1*(47 + 60*shiftst1 + 6*shiftst3))*pow2(Mst1)*pow4(
        Mst2) - (6649153 - 6480*shiftst3 + 540*lmMst1*(-151 + 120*shiftst1 +
        36*shiftst3) - 540*lmMst2*(-143 + 120*shiftst1 + 36*shiftst3))*pow6(
        Mst1) + 1080*(19 + 30*shiftst1 + 3*shiftst3)*pow6(Mst2))))) + MuSUSY*(-
        3240*xDR2DRMOD*(-16*Mst2*xDmglst2*pow2(Dmglst2)*(2*pow2(Mst1)*pow2(
        Mst2)*(-2*s2t*(4*(-2 + 3*lmMst2)*MuSUSY + 5*Mst2*Tbeta)*pow2(Mt) +
        Mst2*Mt*((1 - 10*lmMst1 + 12*lmMst2)*MuSUSY + 6*(-2 + 3*lmMst2)*Mst2*
        Tbeta)*pow2(s2t) + 8*Tbeta*pow3(Mt) + (2 + 5*lmMst1 - 6*lmMst2)*Tbeta*
        pow3(Mst2)*pow3(s2t)) + (32*(2 - 3*lmMst2)*MuSUSY*s2t*pow2(Mt) + 4*
        Mst2*Mt*((-2 - 5*lmMst1 + 6*lmMst2)*MuSUSY + 3*(-2 + 3*lmMst2)*Mst2*
        Tbeta)*pow2(s2t) - 16*(7 + 2*lmMst2)*Tbeta*pow3(Mt) + 5*Tbeta*pow3(
        Mst2)*pow3(s2t))*pow4(Mst1) - 5*s2t*(-2*Mt*MuSUSY*s2t - 4*Tbeta*pow2(
        Mt) + Tbeta*pow2(Mst2)*pow2(s2t))*pow5(Mst2)) + 4*Tbeta*(pow2(Mst1)*(-
        48*(Dmglst2 + 3*Dmglst2*lmMst2 + Mst2 + lmMst2*Mst2)*Mt*pow2(Mst2)*
        pow2(s2t)*(pow2(Mst1)*pow2(Mst2) + pow4(Mst1) + pow4(Mst2)) + 64*pow3(
        Mt)*(-(Dmglst2*(7 + lmMst2)*pow2(Mst1)*pow2(Mst2)) - 2*Dmglst2*(11 + 5*
        lmMst2)*pow4(Mst1) + Dmglst2*(-1 + lmMst2)*pow4(Mst2) + (1 + lmMst2)*
        Mst2*(3*pow2(Mst1)*pow2(Mst2) + 6*pow4(Mst1) + pow4(Mst2)))) + (4*(7*
        Dmglst2 - 13*Mst2)*s2t*(pow2(Mst1) - pow2(Mst2))*pow2(Mt) + pow3(s2t)*(
        2*(4 + 13*lmMst1 - 9*lmMst2)*pow2(Mst1)*pow3(Mst2) + 13*Mst2*pow4(Mst1)
        + Dmglst2*(2*(-7*lmMst1 + 15*lmMst2)*pow2(Mst1)*pow2(Mst2) - 7*pow4(
        Mst1) + 7*pow4(Mst2)) - 13*pow5(Mst2)))*pow5(Mst2)) + Mt*MuSUSY*s2t*(
        256*(Dmglst2 + 3*Dmglst2*lmMst2 + Mst2 + lmMst2*Mst2)*Mt*pow2(Mst1)*(2*
        pow2(Mst1)*pow2(Mst2) + 3*pow4(Mst1) + pow4(Mst2)) + 8*Mst2*s2t*(2*(
        Dmglst2*(7*lmMst1 - 15*lmMst2) + (-4 - 13*lmMst1 + 9*lmMst2)*Mst2)*(
        pow2(Mst1) + pow2(Mst2))*pow4(Mst1) - Dmglst2*((7 - 14*lmMst1 + 30*
        lmMst2)*pow2(Mst1) + 7*pow2(Mst2))*pow4(Mst2) + (5 - 26*lmMst1 + 18*
        lmMst2)*pow2(Mst1)*pow5(Mst2) + 13*pow7(Mst2)))) + 3*Tbeta*(3*pow3(
        Mst2)*(-2*(623599 + 65160*lmMst1 - 134280*lmMst2)*Mst2*s2t*pow2(Mt) +
        15*(11075 + 17928*lmMst1 - 17352*lmMst2)*Mt*pow2(Mst2)*pow2(s2t) + 64*(
        41936 - 7245*lmMst1 + 19665*lmMst2 - 720*lmMt)*pow3(Mt) + 4*(224837 -
        4680*lmMst1 + 8370*lmMst2 + 2700*shiftst1)*pow3(Mst2)*pow3(s2t))*pow4(
        Mst1) + 9*pow2(Mst1)*(-2*Mst2*s2t*(1367 + 13560*lmMst1 - 36600*lmMst2 +
        7200*shiftst1)*pow2(Mt) + 45*(2269 + 664*lmMst1 - 56*lmMst2)*Mt*pow2(
        Mst2)*pow2(s2t) + 4*(66773 - 9960*lmMst1 + 45480*lmMst2 - 3840*lmMt)*
        pow3(Mt) + 10*(5825 + 36*shiftst3 + 12*lmMst1*(47 + 60*shiftst1 + 3*
        shiftst3) - 12*lmMst2*(92 + 60*shiftst1 + 3*shiftst3))*pow3(Mst2)*pow3(
        s2t))*pow5(Mst2) + Mst2*(-4*(2580913 + 203040*lmMst1 - 306720*lmMst2)*
        Mst2*s2t*pow2(Mt) + 6*(30643 + 217080*lmMst1 - 212760*lmMst2)*Mt*pow2(
        Mst2)*pow2(s2t) + 16*(1077991 - 184140*lmMst1 + 400140*lmMst2 - 8640*
        lmMt)*pow3(Mt) + (3447379 - 76140*lmMst1 + 76140*lmMst2)*pow3(Mst2)*
        pow3(s2t))*pow6(Mst1) - 100*Dmsqst2*(27*pow2(Mst1)*pow3(Mst2)*(2*Mst2*
        s2t*(53 + 24*shiftst1)*pow2(Mt) - 288*Mt*pow2(Mst2)*pow2(s2t) + 192*
        lmMt*pow3(Mt) - 65*pow3(Mst2)*pow3(s2t) + 24*lmMst1*(-6*Mt*pow2(Mst2)*
        pow2(s2t) + 8*pow3(Mt) - shiftst1*pow3(Mst2)*pow3(s2t)) - 24*lmMst2*(-
        6*Mt*pow2(Mst2)*pow2(s2t) + 16*pow3(Mt) - shiftst1*pow3(Mst2)*pow3(s2t)
        )) - 18*Mst2*(-361*Mst2*s2t*pow2(Mt) + 216*(3 + 2*lmMst1 - 2*lmMst2)*
        Mt*pow2(Mst2)*pow2(s2t) - 576*(-1 + lmMst1 - 2*lmMst2 + lmMt)*pow3(Mt)
        + (469 - 36*lmMst1 + 36*lmMst2 + 18*shiftst1)*pow3(Mst2)*pow3(s2t))*
        pow4(Mst1) + s2t*(2*(-972*(5 + 4*lmMst1 - 4*lmMst2)*Mst2*Mt*s2t + 4057*
        pow2(Mt) - 5414*pow2(Mst2)*pow2(s2t))*pow6(Mst1) + 324*(1 + shiftst1)*(
        -4*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow6(Mst2))) - 1080*s2t*(19 + 30*
        shiftst1 + 3*shiftst3)*(-4*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow8(Mst2))
        - Dmglst2*Tbeta*(27*(pow2(Mst2)*(-32*(19579 + 1770*lmMst1 + 12630*
        lmMst2)*Mst2*s2t*pow2(Mt) + (-935323 + 279000*lmMst1 - 385560*lmMst2)*
        Mt*pow2(Mst2)*pow2(s2t) + 32*(67843 - 10050*lmMst1 + 17490*lmMst2 -
        7560*lmMt)*pow3(Mt) + 4*(32663 - 7620*lmMst1 + 7620*lmMst2)*pow3(Mst2)*
        pow3(s2t))*pow4(Mst1) + pow2(Mst1)*(32*(-4429 + 270*lmMst1 - 8910*
        lmMst2)*Mst2*s2t*pow2(Mt) + 3*(-173947 + 25080*lmMst1 - 68760*lmMst2)*
        Mt*pow2(Mst2)*pow2(s2t) - 60*(-3245 + 1672*lmMst1 - 1448*lmMst2 + 1888*
        lmMt)*pow3(Mt) + 20*(1799 - 648*lmMst1 + 1680*lmMst2)*pow3(Mst2)*pow3(
        s2t))*pow4(Mst2)) + 2*(-72*(510139 + 44280*lmMst1 + 76680*lmMst2)*Mst2*
        s2t*pow2(Mt) + 15*(-1700119 + 484056*lmMst1 - 579960*lmMst2)*Mt*pow2(
        Mst2)*pow2(s2t) + 8*(13463149 - 1492020*lmMst1 + 2713500*lmMst2 -
        625320*lmMt)*pow3(Mt) + 8*(788723 - 80595*lmMst1 + 80595*lmMst2)*pow3(
        Mst2)*pow3(s2t))*pow6(Mst1) + 97200*Dmsqst2*(2*(pow2(Mst1)*pow2(Mst2)*(
        20*Mst2*s2t*pow2(Mt) + 6*(2 + 3*lmMst1 - 3*lmMst2)*Mt*pow2(Mst2)*pow2(
        s2t) + 8*(8 + 6*lmMst2 - 3*(lmMst1 + lmMt))*pow3(Mt) - (5 + 3*lmMst1 -
        3*lmMst2)*pow3(Mst2)*pow3(s2t)) + (36*Mst2*s2t*pow2(Mt) + 18*(-1 + 2*
        lmMst1 - 2*lmMst2)*Mt*pow2(Mst2)*pow2(s2t) - 80*(-3 + lmMst1 - 2*lmMst2
        + lmMt)*pow3(Mt) + (3 - 5*lmMst1 + 5*lmMst2)*pow3(Mst2)*pow3(s2t))*
        pow4(Mst1)) + 8*s2t*pow2(Mt)*pow5(Mst2) + (6*(-17 + 12*lmMst1 - 12*
        lmMst2)*Mt + 13*Mst2*s2t)*pow2(s2t)*pow6(Mst1) - 2*pow3(s2t)*pow7(Mst2)
        ) - 90720*(-4*s2t*pow2(Mt)*pow7(Mst2) + pow3(s2t)*pow9(Mst2)))))))/(
        204120.*pow2(Mst1))))/(Tbeta*pow9(Mst2)) - (threeLoopFlag*pow2(Al4p)*((
        Mt*((-4*Mt*pow2(MuSUSY)*(-27*pow2(Mst1)*pow2(Mst2)*(15*Mst2*(32*Mst2*(-
        15*Dmsqst2*T1ep*pow2(s2t) + pow2(Mt)*(436 - 6*B4 + 6*D3 - 3*DN - 48*
        lmMst1 + (408 - 96*lmMst1)*lmMst2 + 24*lmMt - 972*S2 - 48*(-4 + lmMst1)
        *pow2(lmMst2) + 48*pow3(lmMst2))) + Mt*s2t*(5760*Dmsqst2*(1 + lmMst1 -
        lmMst2) + pow2(Mst2)*(28683 + 5280*B4 - 48*DN - 5952*lmMst1 - 224*OepS2
        + 324*(-1 + 14*lmMst1 - 14*lmMst2)*S2 + 336*T1ep - 768*pow2(lmMst1) -
        192*lmMst2*(-214 + 73*lmMst1 + 4*pow2(lmMst1)) - 96*(-268 + 57*lmMst1)*
        pow2(lmMst2) + 6240*pow3(lmMst2))) - 1104*T1ep*pow2(s2t)*pow3(Mst2)) +
        Dmglst2*(2880*Mst2*pow2(Mt)*(180 - 2*B4 + 2*D3 - DN + 144*lmMst2 - 216*
        S2 - 16*(-2 + lmMst1)*pow2(lmMst2) + 16*(lmMst1 + pow3(lmMst2))) + Mt*
        s2t*(-86400*Dmsqst2*(2 + 5*lmMst1 - 5*lmMst2) + pow2(Mst2)*(23917 +
        188640*B4 - 3600*DN - 37440*lmMst1 + 5600*OepS2 - 324*(-453 + 350*
        lmMst1 - 350*lmMst2)*S2 - 8400*T1ep + 11520*pow2(lmMst1) - 2880*lmMst2*
        (-237 + 55*lmMst1 + 4*pow2(lmMst1)) - 1440*(-280 + 121*lmMst1)*pow2(
        lmMst2) + 185760*pow3(lmMst2))) + 3360*T1ep*pow2(s2t)*pow3(Mst2))) -
        36*(5*Mst2*(9*Mst2*(-520*Dmsqst2*T1ep*pow2(s2t) + pow2(Mt)*(10667 - 96*
        B4 + 96*D3 - 48*DN - 3072*lmMst1 - 384*(-25 + 9*lmMst1)*lmMst2 - 384*
        lmMst1*lmMt + 384*(1 + lmMst2)*lmMt - 224*OepS2 + 324*(-43 + 14*lmMst1
        - 14*lmMst2)*S2 + 336*T1ep - 384*(-13 + 4*lmMst1)*pow2(lmMst2) + 1536*
        pow3(lmMst2))) + Mt*s2t*(12960*Dmsqst2*(1 + 3*lmMst1 - 3*lmMst2) + 2*
        pow2(Mst2)*(75569 + 13716*B4 - 54*DN - 33426*lmMst1 - 1088*OepS2 + 162*
        (169 + 136*lmMst1 - 136*lmMst2)*S2 + 1632*T1ep - 2376*pow2(lmMst1) +
        54*lmMst2*(1427 - 1012*lmMst1 + 16*pow2(lmMst1)) - 108*(-642 + 203*
        lmMst1)*pow2(lmMst2) + 21060*pow3(lmMst2))) - 12684*T1ep*pow2(s2t)*
        pow3(Mst2)) + Dmglst2*(30*Mst2*pow2(Mt)*(28405 - 288*B4 + 288*D3 - 144*
        DN + 5472*lmMst1 - 288*(-89 + 2*lmMst1)*lmMst2 + 576*(-2 + lmMst1 -
        lmMst2)*lmMt + 224*OepS2 - 324*(65 + 14*lmMst1 - 14*lmMst2)*S2 - 336*
        T1ep - 576*(-9 + 8*lmMst1)*pow2(lmMst2) + 4608*pow3(lmMst2)) + Mt*s2t*(
        64800*Dmsqst2*(8 - 15*lmMst1 + 15*lmMst2) + pow2(Mst2)*(66761 + 301320*
        B4 - 4860*DN - 205380*lmMst1 + 40480*OepS2 - 216*(2489 + 3795*lmMst1 -
        3795*lmMst2)*S2 - 60720*T1ep + 23760*pow2(lmMst1) + 180*lmMst2*(4993 -
        1956*lmMst1 + 48*pow2(lmMst1)) - 1080*(-482 + 331*lmMst1)*pow2(lmMst2)
        + 348840*pow3(lmMst2))) + 30960*T1ep*pow2(s2t)*pow3(Mst2)))*pow4(Mst1)
        - 10*(-12*(13260*Dmsqst2 + Mst2*(-34616*Dmglst2 + 39711*Mst2))*T1ep*
        pow2(s2t) + 9*pow2(Mt)*(383185 - 2592*B4 + 2592*D3 - 1296*DN - 187704*
        lmMst1 - 17440*OepS2 + 648*(-57 + 545*lmMst1 - 545*lmMst2)*S2 + 26160*
        T1ep - 7992*pow2(lmMst1) - 216*lmMst2*(-1733 + 630*lmMst1 + 26*pow2(
        lmMst1)) - 216*(-859 + 246*lmMst1)*pow2(lmMst2) + 3456*lmMt*(3 + 6*
        lmMst2 - 2*lmMst1*(3 + lmMst2) + pow2(lmMst1) + pow2(lmMst2)) + 720*
        pow3(lmMst1) + 58032*pow3(lmMst2)) - Mt*s2t*(Dmglst2*(2773621 -
        1660176*B4 + 25272*DN + 2004408*lmMst1 - 525472*OepS2 + 108*(123113 +
        98526*lmMst1 - 98526*lmMst2)*S2 + 788208*T1ep + 3888*pow2(lmMst1) -
        144*lmMst2*(36802 - 11421*lmMst1 + 1728*pow2(lmMst1)) + 167184*(-14 +
        15*lmMst1)*pow2(lmMst2) - 31104*pow3(lmMst1) - 2227824*pow3(lmMst2)) -
        3*Mst2*(1702429 + 257904*B4 - 648*DN - 748656*lmMst1 - 30112*OepS2 +
        108*(9185 + 5646*lmMst1 - 5646*lmMst2)*S2 + 45168*T1ep + 41904*pow2(
        lmMst1) + 216*lmMst2*(5971 - 6106*lmMst1 + 576*pow2(lmMst1)) - 41904*(-
        34 + 15*lmMst1)*pow2(lmMst2) - 3456*pow3(lmMst1) + 507600*pow3(lmMst2))
        ))*pow6(Mst1) - 622080*(1 + lmMst2)*(-Dmglst2 + 3*Dmglst2*lmMst2 + Mst2
        + lmMst2*Mst2)*Mt*s2t*pow6(Mst2)))/(15.*pow2(Mst1)*pow6(Mst2)) - 11664*
        Mt*pow2(s2t)*((5*Dmglst2*Dmsqst2*(-1081 + 165*lmMst1 - 165*lmMst2)*(-1
        + pow2(Sbeta))*pow2(Sbeta)*pow4(Mst1))/(9.*pow5(Mst2)) + pow2(MuSUSY)*(
        53.385802469135804 + (40*B4)/9. - (4*D3)/9. + (2*DN)/9. + (1672*lmMst1)
        /27. + (53*pow2(lmMst1))/9. - lmMst2*(129.92592592592592 - 72*lmMst1 +
        (13*pow2(lmMst1))/3.) + (2*(-470 + 147*lmMst1)*pow2(lmMst2))/9. - ((103
        + 186*lmMst2 + 32*lmMst1*(1 + lmMst2) + 91*pow2(lmMst2))*pow2(Mst2))/(
        9.*pow2(Mst1)) + (2250*Dmsqst2*(1141 + 48*lmMst1*(7 - 3*lmMst2) - 264*
        lmMst2 + 144*pow2(lmMst2)) + pow2(Mst1)*(10552777 + 205200*B4 - 5400*DN
        + 4145280*lmMst1 - 1181700*pow2(lmMst1) - 240*lmMst2*(16192 - 26430*
        lmMst1 + 3465*pow2(lmMst1)) + 900*(-5735 + 3072*lmMst1)*pow2(lmMst2) -
        55800*pow3(lmMst1) - 1877400*pow3(lmMst2)))/(24300.*pow2(Mst2)) - (271*
        pow3(lmMst2))/9. + (16*(pow3(lmMst1) + ((1 + lmMst2)*(4*Dmglst2*lmMst2
        + Mst2 + lmMst2*Mst2)*pow3(Mst2))/pow4(Mst1)))/9. + (Dmsqst2*(
        201.74098765432097 - (622*lmMst2)/15. + (2*lmMst1*(311 + 10*lmMst2))/
        15. - (22*pow2(lmMst1))/3. + 6*pow2(lmMst2))*pow2(Mst1) + (
        628.1736268201578 + (76*B4)/9. - (2*DN)/9. + (6317839*lmMst1)/396900. -
        (66307*pow2(lmMst1))/315. - lmMst2*(12.52907281431091 - (182909*lmMst1)
        /315. + (274*pow2(lmMst1))/3.) + (2*(-58301 + 37135*lmMst1)*pow2(
        lmMst2))/315. - (44*pow3(lmMst1))/9. - (1256*pow3(lmMst2))/9.)*pow4(
        Mst1))/pow4(Mst2) - Dmglst2*((2*Mst2*(-20 + 2*(99 + 16*lmMst1)*lmMst2 +
        123*pow2(lmMst2)))/(9.*pow2(Mst1)) + (14267 - 432*B4 + 576*D3 - 360*DN
        + 4752*lmMst1 - 1404*pow2(lmMst1) + 72*lmMst2*(-63 - 232*lmMst1 + 16*
        pow2(lmMst1)) + 72*(281 - 123*lmMst1)*pow2(lmMst2) + (6480*Dmsqst2)/
        pow2(Mst1) + 7704*pow3(lmMst2))/(162.*Mst2) + (810000*Dmsqst2*(4 +
        lmMst1 - lmMst2) + 2*pow2(Mst1)*(2695042 - 40500*B4 + 54000*D3 - 33750*
        DN - 326895*lmMst1 - 324900*pow2(lmMst1) + 15*lmMst2*(-19607 - 129030*
        lmMst1 + 62550*pow2(lmMst1)) - 450*(-5023 + 5610*lmMst1)*pow2(lmMst2) +
        11250*pow3(lmMst1) + 1575000*pow3(lmMst2)))/(30375.*pow3(Mst2)) + ((
        585.1892843532082 - (8*B4)/3. + (32*D3)/9. - (20*DN)/9. - (20109937*
        lmMst1)/297675. - (15886*pow2(lmMst1))/945. + lmMst2*(
        17.112243218274966 - (144628*lmMst1)/945. + (1180*pow2(lmMst1))/9.) + (
        2*(80257 - 138180*lmMst1)*pow2(lmMst2))/945. - (92*pow3(lmMst1))/27. +
        (4448*pow3(lmMst2))/27.)*pow4(Mst1))/pow5(Mst2)) - (Dmsqst2*(10 - 20*
        lmMst2 + (10*Dmglst2*(-23 + 42*lmMst1 - 42*lmMst2)*pow4(Mst1))/pow5(
        Mst2) - (3*(237.28785508324435 + (16526*lmMst2)/3969. + (2*lmMst1*(-
        8263 + 71820*lmMst2))/3969. - (520*pow2(lmMst1))/21. - (80*pow2(lmMst2)
        )/7.)*pow6(Mst1))/pow6(Mst2)))/(3.*pow2(Mst1)) + ((10*shiftst1*(
        Dmsqst2*(3 - 2*lmMst2) + (1 - 2*lmMst2)*pow2(Mst2))*((pow2(Mst1) +
        pow2(Mst2))*pow4(Mst2) - 2*(lmMst1 - lmMst2)*pow2(Mst1)*(pow2(Mst1)*
        pow2(Mst2) + pow4(Mst1) + pow4(Mst2))))/(3.*pow2(Mst1)) + (S2*(30*
        Dmsqst2*(9*(-17 + 78*lmMst1 - 78*lmMst2)*pow2(Mst1)*pow2(Mst2) + 13*(95
        + 102*lmMst1 - 102*lmMst2)*pow4(Mst1) + 81*(-15 + 2*lmMst1 - 2*lmMst2)*
        pow4(Mst2)) + Mst2*(3*Mst2*(9*(11243 + 2114*lmMst1 - 2114*lmMst2)*pow2(
        Mst1)*pow2(Mst2) + (194357 + 39711*lmMst1 - 39711*lmMst2)*pow4(Mst1) +
        81*(307 + 46*lmMst1 - 46*lmMst2)*pow4(Mst2)) + Dmglst2*(36*(15643 -
        774*lmMst1 + 774*lmMst2)*pow2(Mst1)*pow2(Mst2) + 8*(93919 - 12981*
        lmMst1 + 12981*lmMst2)*pow4(Mst1) + 162*(1677 - 14*lmMst1 + 14*lmMst2)*
        pow4(Mst2)))))/81. - (4*OepS2*(60*Dmsqst2*(117*pow2(Mst1)*pow2(Mst2) +
        221*pow4(Mst1) + 27*pow4(Mst2)) + Mst2*(19026*pow2(Mst1)*pow3(Mst2) +
        39711*Mst2*pow4(Mst1) - 4*Dmglst2*(2322*pow2(Mst1)*pow2(Mst2) + 8654*
        pow4(Mst1) + 189*pow4(Mst2)) + 3726*pow5(Mst2))))/2187.)/pow6(Mst2) + (
        (1 - 2*lmMst2)*shiftst3*(2*(1 - 2*lmMst1 + 2*lmMst2)*pow2(Mst2)*pow4(
        Mst1) + 2*(1 - lmMst1 + lmMst2)*pow2(Mst1)*pow4(Mst2) + (2 - 6*lmMst1 +
        6*lmMst2)*pow6(Mst1) + pow6(Mst2)))/(3.*pow2(Mst1)*pow4(Mst2)))) + (16*
        z3*(2*MuSUSY*(2*Mt*(6759*MuSUSY + 834482*Dmglst2*Tbeta - 321216*Mst2*
        Tbeta) + s2t*(1582324*Dmglst2*MuSUSY + 298896*Mst2*MuSUSY + 28950*
        Dmsqst2*Tbeta + 50508*Dmglst2*Mst2*Tbeta + 40629*Tbeta*pow2(Mst2)))*
        pow2(Mt) + Mt*pow2(s2t)*(2*Mst2*MuSUSY*(Dmglst2*(706232*MuSUSY -
        546087*Mst2*Tbeta) + 6*Mst2*((91963 + 162*lmMst1 - 162*lmMst2)*MuSUSY -
        4431*Mst2*Tbeta)) + 3645*xDmsqst2*pow2(Dmsqst2)*(-1 + pow2(Sbeta))*
        pow2(Sbeta) + 30*Dmsqst2*(3888*(Dmglst2 - Mst2)*MuSUSY*Tbeta + 1771*
        pow2(MuSUSY) - 810*Dmglst2*Mst2*(-1 + pow2(Sbeta))*pow2(Sbeta))) - ((
        21840*Dmsqst2 + 11*Mst2*(16354*Dmglst2 + 20937*Mst2))*MuSUSY*Tbeta*
        pow2(Mst2) - 30*xDmsqst2*pow2(Dmsqst2)*(1517*MuSUSY*Tbeta - 405*Mst2*(-
        1 + pow2(Sbeta))*pow2(Sbeta)))*pow3(s2t))*pow4(Mst1) + 8*MuSUSY*z4*((4*
        s2t*(2636*Dmglst2*MuSUSY + 35292*Mst2*MuSUSY - 1590*Dmsqst2*Tbeta +
        27252*Dmglst2*Mst2*Tbeta - 18093*Tbeta*pow2(Mst2))*pow2(Mt) - 2*Mt*(
        13260*Dmsqst2*MuSUSY + Mst2*(-45308*Dmglst2*MuSUSY + 41169*Mst2*MuSUSY
        - 15189*Dmglst2*Mst2*Tbeta + 22293*Tbeta*pow2(Mst2)))*pow2(s2t) + 8*(
        9279*MuSUSY - 28538*Dmglst2*Tbeta + 9978*Mst2*Tbeta)*pow3(Mt) + Tbeta*(
        -2820*xDmsqst2*pow2(Dmsqst2) + 6240*Dmsqst2*pow2(Mst2) + (-25328*
        Dmglst2 + 20685*Mst2)*pow3(Mst2))*pow3(s2t))*pow4(Mst1) - 3*pow2(Mst1)*
        (2*xDmglst2*pow2(Dmglst2)*(-4*Mst2*s2t*(4045*MuSUSY + 594*Mst2*Tbeta)*
        pow2(Mt) + Mt*(-14978*MuSUSY + 10173*Mst2*Tbeta)*pow2(Mst2)*pow2(s2t) -
        8*(507*MuSUSY + 2198*Mst2*Tbeta)*pow3(Mt) + 2731*Tbeta*pow3(s2t)*pow4(
        Mst2)) + 3*(-8*(1135*MuSUSY*s2t + 652*Mt*Tbeta)*pow2(Mt)*pow3(Mst2) -
        20*pow2(Mst2)*(-30*Dmsqst2*s2t*Tbeta*pow2(Mt) - 78*Dmsqst2*Mt*MuSUSY*
        pow2(s2t) + 180*MuSUSY*pow3(Mt) - 11*Tbeta*xDmsqst2*pow2(Dmsqst2)*pow3(
        s2t)) + s2t*(-80*Mt*(7*MuSUSY*s2t + 4*Mt*Tbeta)*xDmsqst2*pow2(Dmsqst2)
        + (4552*Mt*MuSUSY*s2t + 3470*Tbeta*pow2(Mt) - 600*Dmsqst2*Tbeta*pow2(
        s2t))*pow4(Mst2)) + Dmglst2*Mst2*(8*Mst2*s2t*(-709*MuSUSY + 24*Mst2*
        Tbeta)*pow2(Mt) + 15*Mt*(-296*MuSUSY + 121*Mst2*Tbeta)*pow2(Mst2)*pow2(
        s2t) - 96*(47*MuSUSY - 26*Mst2*Tbeta)*pow3(Mt) + 948*Tbeta*pow3(s2t)*
        pow4(Mst2)) + (4515*Mt - 1916*Mst2*s2t)*Tbeta*pow2(s2t)*pow5(Mst2))) -
        9*pow2(Mst2)*(4*xDmglst2*pow2(Dmglst2)*(2*Mst2*s2t*(-218*MuSUSY + 99*
        Mst2*Tbeta)*pow2(Mt) + Mt*(-1586*MuSUSY + 327*Mst2*Tbeta)*pow2(Mst2)*
        pow2(s2t) - 4*(81*MuSUSY + 28*Mst2*Tbeta)*pow3(Mt) + 793*Tbeta*pow3(
        s2t)*pow4(Mst2)) + 3*(-40*Mt*s2t*(MuSUSY*s2t + Mt*Tbeta)*xDmsqst2*pow2(
        Dmsqst2) - 4*Dmglst2*Mst2*(Mst2*s2t*(271*MuSUSY - 216*Mst2*Tbeta) + Mt*
        (216*MuSUSY + 253*Mst2*Tbeta))*pow2(Mt) - 12*(85*MuSUSY*s2t + 103*Mt*
        Tbeta)*pow2(Mt)*pow3(Mst2) + Dmglst2*Mt*(-848*MuSUSY + 813*Mst2*Tbeta)*
        pow2(s2t)*pow3(Mst2) - 4*pow2(Mst2)*(-30*Dmsqst2*s2t*Tbeta*pow2(Mt) -
        30*Dmsqst2*Mt*MuSUSY*pow2(s2t) + 108*MuSUSY*pow3(Mt) - 5*Tbeta*
        xDmsqst2*pow2(Dmsqst2)*pow3(s2t)) + 6*s2t*(40*Mt*MuSUSY*s2t + 101*
        Tbeta*pow2(Mt) - 10*Dmsqst2*Tbeta*pow2(s2t))*pow4(Mst2) - 15*(-51*Mt +
        8*Mst2*s2t)*Tbeta*pow2(s2t)*pow5(Mst2) + 424*Dmglst2*Tbeta*pow3(s2t)*
        pow5(Mst2)))) - (16*xDmsqst2*pow2(Dmsqst2)*(-9*Mt*s2t*(2160*(3 + 2*
        lmMst1 - 2*lmMst2)*Mst2*Mt*pow2(Mst1) + 270*s2t*shiftst1*(2*(lmMst1 -
        lmMst2)*pow2(Mst1) - pow2(Mst2))*(pow2(Mst1) + pow2(Mst2)) + 40*s2t*
        T1ep*pow2(Mst1)*(14*pow2(Mst1) + 3*pow2(Mst2)))*pow2(MuSUSY) + MuSUSY*
        Tbeta*((45*Mst2*pow2(Mst1)*(16*Mst2*s2t*(88 + 8*OepS2 - 81*S2 - 18*(
        lmMst1 - lmMst2)*(2 + 9*S2) - 54*shiftst1 - 12*T1ep)*pow2(Mt) + 2592*(3
        + 2*lmMst1 - 2*lmMst2)*Mt*pow2(Mst2)*pow2(s2t) + 864*(4 - 8*lmMst1 +
        15*lmMst2 - 7*lmMt)*pow3(Mt) - (1493 + 64*OepS2 - 27864*S2 - 72*(lmMst1
        - lmMst2)*(-1 + 18*S2 + 6*shiftst1) - 96*T1ep)*pow3(Mst2)*pow3(s2t)))/
        8. + s2t*((3*(972000*(1 + 4*lmMst1 - 4*lmMst2)*Mst2*Mt*s2t + 8*(187757
        + 16000*OepS2 - 540000*S2 + 90*lmMst2*(-587 + 3600*S2) - 90*lmMst1*(-
        587 + 630*lmMst2 + 3600*S2) - 24000*T1ep + 28350*(pow2(lmMst1) + pow2(
        lmMst2)))*pow2(Mt) + (1372861 - 88000*OepS2 + 180*lmMst2*(1987 - 9900*
        S2) + 13743000*S2 + 180*lmMst1*(-1987 + 1080*lmMst2 + 9900*S2) + 81000*
        shiftst1 + 132000*T1ep - 97200*(pow2(lmMst1) + pow2(lmMst2)))*pow2(
        Mst2)*pow2(s2t))*pow4(Mst1))/200. - 1215*(2 + shiftst1)*(-4*pow2(Mt) +
        pow2(Mst2)*pow2(s2t))*pow4(Mst2)) + ((610359812 - 64484000*OepS2 + 315*
        lmMst2*(540413 - 4145400*S2) + 8117266500*S2 + 315*lmMst1*(-540413 +
        884520*lmMst2 + 4145400*S2) + 96726000*T1ep - 139311900*(pow2(lmMst1) +
        pow2(lmMst2)))*pow3(s2t)*pow6(Mst1))/34300.) + (3*pow2(s2t)*(10125*
        pow2(Sbeta)*(2*(11 + 39*lmMst1 - 39*lmMst2)*Mst2*s2t*(-1 + pow2(Sbeta))
        + (167 - 5*lmMst1 + 5*lmMst2)*Mt*pow2(Sbeta))*pow6(Mst1) + Mt*(pow2(
        MuSUSY)*(375*(1925 + 64*OepS2 - 27864*S2 - 72*(lmMst1 - lmMst2)*(-1 +
        18*S2))*pow2(Mst1)*pow2(Mst2) + 2*(7*(-46499 + 8000*OepS2 - 1728000*S2)
        + 90*lmMst2*(-2137 + 12600*S2) - 90*lmMst1*(-2137 + 1080*lmMst2 +
        12600*S2) + 48600*(pow2(lmMst1) + pow2(lmMst2)))*pow4(Mst1) + 162000*
        pow4(Mst2)) + 10125*(-167 + 5*lmMst1 - 5*lmMst2)*pow2(Sbeta)*pow6(Mst1)
        )))/100.))/pow2(Mst1) - MuSUSY*z3*(3*pow2(Mst1)*(384*Mt*(5*s2t*(472*
        MuSUSY*s2t + 28*Mt*Tbeta + 81*Mst2*s2t*Tbeta)*xDmsqst2*pow2(Dmsqst2) +
        540*Dmsqst2*Mst2*(-7*Mt*MuSUSY*s2t + 4*Tbeta*pow2(Mt) + 3*Tbeta*pow2(
        Mst2)*pow2(s2t)) + Mt*(-7318*MuSUSY*s2t + 8219*Mt*Tbeta)*pow3(Mst2)) -
        3*pow2(Mst2)*(96960*Dmsqst2*s2t*Tbeta*pow2(Mt) + 16800*Dmsqst2*Mt*
        MuSUSY*pow2(s2t) + 125568*MuSUSY*pow3(Mt) + 84995*Tbeta*xDmsqst2*pow2(
        Dmsqst2)*pow3(s2t)) - 48*s2t*(2*(35719 + 108*lmMst1 - 108*lmMst2)*Mt*
        MuSUSY*s2t - 1859*Tbeta*pow2(Mt) - 1560*Dmsqst2*Tbeta*pow2(s2t))*pow4(
        Mst2) + 4*xDmglst2*pow2(Dmglst2)*(-32*Mst2*s2t*(80959*MuSUSY + 22806*
        Mst2*Tbeta)*pow2(Mt) + 4*Mt*(-458284*MuSUSY + 190797*Mst2*Tbeta)*pow2(
        Mst2)*pow2(s2t) + 8*(154668*MuSUSY + 342281*Mst2*Tbeta)*pow3(Mt) +
        547817*Tbeta*pow3(s2t)*pow4(Mst2)) + 24*(34131*Mt + 4*(11924 + 27*
        lmMst1 - 27*lmMst2)*Mst2*s2t)*Tbeta*pow2(s2t)*pow5(Mst2) + 24*Dmglst2*(
        43200*Dmsqst2*(MuSUSY*s2t - 4*Mt*Tbeta)*pow2(Mt) + 24*s2t*(-9747*Mt*
        MuSUSY*s2t - 1180*Tbeta*pow2(Mt) + 225*Dmsqst2*Tbeta*pow2(s2t))*pow3(
        Mst2) + Mst2*Mt*(-379648*Mst2*Mt*MuSUSY*s2t + 2400*(5*MuSUSY - 67*Mst2*
        Tbeta)*pow2(Mt) - 3*pow2(s2t)*(2880*Dmsqst2*(2*MuSUSY + 3*Mst2*Tbeta) -
        49501*Tbeta*pow3(Mst2))) + 64119*Tbeta*pow3(s2t)*pow5(Mst2))) + 9*Mst2*
        (4*Dmglst2*Mst2*(Dmglst2*xDmglst2*(-2*Mst2*Mt*s2t*(Mst2*s2t*(122917*
        MuSUSY - 196638*Mst2*Tbeta) + Mt*(262184*MuSUSY + 166791*Mst2*Tbeta)) +
        8*(21888*MuSUSY + 32323*Mst2*Tbeta)*pow3(Mt) + 122917*Tbeta*pow3(s2t)*
        pow4(Mst2)) + 6*(2880*Dmsqst2*(MuSUSY*s2t - 6*Mt*Tbeta)*pow2(Mt) + 2*
        s2t*(-17615*Mt*MuSUSY*s2t - 1756*Tbeta*pow2(Mt) + 540*Dmsqst2*Tbeta*
        pow2(s2t))*pow3(Mst2) + Mst2*Mt*(-60548*Mst2*Mt*MuSUSY*s2t + 4*(664*
        MuSUSY - 9*Mst2*Tbeta)*pow2(Mt) - 3*pow2(s2t)*(720*Dmsqst2*(MuSUSY +
        Mst2*Tbeta) - 15137*Tbeta*pow3(Mst2))) + 17615*Tbeta*pow3(s2t)*pow5(
        Mst2))) + 3*(10*Mt*(1152*Mt*(-(MuSUSY*s2t) + 2*Mt*Tbeta) + 7*Mst2*s2t*(
        629*MuSUSY*s2t + 86*Mt*Tbeta))*xDmsqst2*pow2(Dmsqst2) + 32*Mt*pow2(
        Mst2)*(pow2(Mst2)*(-5965*Mt*MuSUSY*s2t + 2883*Tbeta*pow2(Mt) + 1620*
        Dmsqst2*Tbeta*pow2(s2t)) + 90*Dmsqst2*(-24*Mt*MuSUSY*s2t + 16*Tbeta*
        pow2(Mt) + 3*Dmsqst2*Tbeta*xDmsqst2*pow2(s2t))) - 5*pow3(Mst2)*(5568*
        Dmsqst2*s2t*Tbeta*pow2(Mt) - 2208*Dmsqst2*Mt*MuSUSY*pow2(s2t) + 5888*
        MuSUSY*pow3(Mt) + 4403*Tbeta*xDmsqst2*pow2(Dmsqst2)*pow3(s2t)) - 48*
        s2t*(2*(1319 + 6*lmMst1 - 6*lmMst2)*Mt*MuSUSY*s2t + 3*(-77 + 8*lmMst1 -
        8*lmMst2)*Tbeta*pow2(Mt) + 115*Dmsqst2*Tbeta*pow2(s2t))*pow5(Mst2) +
        24*(5965*Mt + 2*(1319 + 6*lmMst1 - 6*lmMst2)*Mst2*s2t)*Tbeta*pow2(s2t)*
        pow6(Mst2)))))/pow6(Mst2)))/Tbeta + 6*MuSUSY*((2*Mt*pow2(z2)*((-4*s2t*(
        172810*Dmglst2*MuSUSY + 41010*Mst2*MuSUSY + 1590*Dmsqst2*Tbeta - 33084*
        Dmglst2*Mst2*Tbeta + 18093*Tbeta*pow2(Mst2))*pow2(Mt) + 2*Mt*(-13260*
        Dmsqst2*MuSUSY + Mst2*(42392*Dmglst2*MuSUSY - 35823*Mst2*MuSUSY +
        105585*Dmglst2*Mst2*Tbeta + 18531*Tbeta*pow2(Mst2)))*pow2(s2t) + 8*(
        4905*MuSUSY - 32426*Dmglst2*Tbeta + 6090*Mst2*Tbeta)*pow3(Mt) + Tbeta*(
        -2820*xDmsqst2*pow2(Dmsqst2) + 6240*Dmsqst2*pow2(Mst2) + (-25328*
        Dmglst2 + 20685*Mst2)*pow3(Mst2))*pow3(s2t))*pow4(Mst1) - 9*pow2(Mst2)*
        (2*Mst2*xDmglst2*pow2(Dmglst2)*(4*s2t*(2536*MuSUSY - 63*Mst2*Tbeta)*
        pow2(Mt) - 2*Mst2*Mt*(209*MuSUSY + 3804*Mst2*Tbeta)*pow2(s2t) - 224*
        Tbeta*pow3(Mt) + 209*Tbeta*pow3(Mst2)*pow3(s2t)) + 3*Dmglst2*pow2(Mst2)
        *(7052*MuSUSY*s2t*pow2(Mt) - Mst2*Mt*(632*MuSUSY + 5289*Mst2*Tbeta)*
        pow2(s2t) + 140*Tbeta*pow3(Mt) + 316*Tbeta*pow3(Mst2)*pow3(s2t)) + 3*(-
        40*Mt*s2t*(MuSUSY*s2t + Mt*Tbeta)*xDmsqst2*pow2(Dmsqst2) + 20*Dmsqst2*
        s2t*pow2(Mst2)*(6*Mt*MuSUSY*s2t + 6*Tbeta*pow2(Mt) + Dmsqst2*Tbeta*
        xDmsqst2*pow2(s2t)) + 12*(185*MuSUSY*s2t - 7*Mt*Tbeta)*pow2(Mt)*pow3(
        Mst2) - 6*s2t*(2*Mt*MuSUSY*s2t - 53*Tbeta*pow2(Mt) + 10*Dmsqst2*Tbeta*
        pow2(s2t))*pow4(Mst2) - 1665*Mt*Tbeta*pow2(s2t)*pow5(Mst2) + 6*Tbeta*
        pow3(s2t)*pow6(Mst2))) + 3*pow2(Mst1)*(-3*Dmglst2*Mst2*(-80*Mst2*s2t*(-
        569*MuSUSY + 30*Mst2*Tbeta)*pow2(Mt) - 3*Mt*(1264*MuSUSY + 6091*Mst2*
        Tbeta)*pow2(Mst2)*pow2(s2t) + 96*(7*MuSUSY + 62*Mst2*Tbeta)*pow3(Mt) +
        948*Tbeta*pow3(s2t)*pow4(Mst2)) + 2*xDmglst2*pow2(Dmglst2)*(4*Mst2*s2t*
        (-12479*MuSUSY + 1080*Mst2*Tbeta)*pow2(Mt) + Mt*(-1060*MuSUSY + 14613*
        Mst2*Tbeta)*pow2(Mst2)*pow2(s2t) + 56*(3*MuSUSY + 314*Mst2*Tbeta)*pow3(
        Mt) + 1157*Tbeta*pow3(s2t)*pow4(Mst2)) + 3*(80*Mt*s2t*(7*MuSUSY*s2t +
        4*Mt*Tbeta)*xDmsqst2*pow2(Dmsqst2) + 32*(-398*MuSUSY*s2t + 55*Mt*Tbeta)
        *pow2(Mt)*pow3(Mst2) - 4*pow2(Mst2)*(150*Dmsqst2*s2t*Tbeta*pow2(Mt) +
        390*Dmsqst2*Mt*MuSUSY*pow2(s2t) - 252*MuSUSY*pow3(Mt) + 55*Tbeta*
        xDmsqst2*pow2(Dmsqst2)*pow3(s2t)) + 2*s2t*(-1682*Mt*MuSUSY*s2t - 1735*
        Tbeta*pow2(Mt) + 300*Dmsqst2*Tbeta*pow2(s2t))*pow4(Mst2) + 4557*Mt*
        Tbeta*pow2(s2t)*pow5(Mst2) + 1700*Tbeta*pow3(s2t)*pow6(Mst2)))))/(
        Tbeta*pow6(Mst2)) + 1944*(-(s2t*pow3(Mt)*(302.29104938271604 + (8*D3)/
        9. - (8*DN)/9. + (614*lmMst1)/27. + (32*lmMt)/3. + (124*pow2(lmMst1))/
        9. - (2*lmMst2*(-1901 + 6*lmMst1 + 117*pow2(lmMst1)))/27. + (46 + (32*
        lmMst1)/9.)*pow2(lmMst2) + (14*pow3(lmMst2))/9. + (13500*Dmsqst2*(367 +
        32*lmMst1 - 80*lmMst2) + pow2(Mst1)*(19425643 + 518400*lmMt + 240*
        lmMst2*(82997 + 4320*lmMt) - 3600*(683 + 96*lmMst2)*pow2(lmMst1) +
        5684400*pow2(lmMst2) - 240*lmMst1*(42047 + 4800*lmMst2 + 4320*lmMt +
        6390*pow2(lmMst2)) - 295200*pow3(lmMst1) + 2174400*pow3(lmMst2)))/(
        48600.*pow2(Mst2)) + (32*(pow3(lmMst1) + ((1 + lmMst2)*(4*Dmglst2*
        lmMst2 + Mst2 + lmMst2*Mst2)*pow3(Mst2))/pow4(Mst1)))/9. + (Dmsqst2*(
        83.94493827160494 - (356*lmMst1)/15. + (4*(89 + 60*lmMst1)*lmMst2)/15.
         - 8*(pow2(lmMst1) + pow2(lmMst2)))*pow2(Mst1) + (
        599.6859861506036 - (99165049*lmMst1)/198450. + lmMst2*(
        713.9201259763164 - (57994*lmMst1)/945. - (350*pow2(lmMst1))/9.) - (
        90913*pow2(lmMst1))/945. + (200.24021164021164 - (286*lmMst1)/9.)*pow2(
        lmMst2) + (32*lmMt*(1 + 4*lmMst2 - 2*lmMst1*(2 + lmMst2) + pow2(lmMst1)
        + pow2(lmMst2)))/3. + (26*pow3(lmMst1))/27. + (1882*pow3(lmMst2))/27.)*
        pow4(Mst1))/pow4(Mst2) + Dmglst2*((-4*Mst2*(-52 + 2*(51 + 16*lmMst1)*
        lmMst2 + 59*pow2(lmMst2)))/(9.*pow2(Mst1)) + (4*(12454 - 120*B4 + 120*
        D3 - 60*DN + 690*lmMst1 + 345*pow2(lmMst1) - 5*lmMst2*(-3121 + 42*
        lmMst1 + 96*pow2(lmMst1)) + (3630 - 480*lmMst1)*pow2(lmMst2) - (2700*
        Dmsqst2)/pow2(Mst1) + 960*pow3(lmMst2)))/(135.*Mst2) - (8*(691875*
        Dmsqst2 - pow2(Mst1)*(1928479 - 13500*B4 + 13500*D3 - 6750*DN + 635385*
        lmMst1 + 81000*(-2 + lmMst1 - lmMst2)*lmMt + 450*pow2(lmMst1) - 15*
        lmMst2*(-153166 + 17835*lmMst1 + 3375*pow2(lmMst1)) - 225*(-2627 +
        1695*lmMst1)*pow2(lmMst2) - 1125*pow3(lmMst1) + 433125*pow3(lmMst2))))/
        (30375.*pow3(Mst2)) + ((311.41128771790903 - (32*B4)/9. + (32*D3)/9. -
        (16*DN)/9. + (14926634*lmMst1)/19845. + (19352*pow2(lmMst1))/189. + (2*
        lmMst2*(2917823 - 3813600*lmMst1 + 97020*pow2(lmMst1)))/19845. + (8*(
        8677 - 5439*lmMst1)*pow2(lmMst2))/189. - (64*lmMt*(4 + 6*lmMst2 - 2*
        lmMst1*(3 + lmMst2) + pow2(lmMst1) + pow2(lmMst2)))/3. - (8*pow3(
        lmMst1))/9. + (664*pow3(lmMst2))/3.)*pow4(Mst1))/pow5(Mst2)) + ((-2*(71
        + 122*lmMst2 + 32*lmMst1*(1 + lmMst2) + 59*pow2(lmMst2))*pow2(Mst2))/9.
         + (Dmsqst2*(40*lmMst2 - 20*(1 + (Dmglst2*(65 - 2*lmMst1 + 2*lmMst2)*
        pow4(Mst1))/pow5(Mst2)) + (3*(87.53310385647498 - (941812*lmMst1)/
        19845. + (47.45840262030738 + (484*lmMst1)/21.)*lmMst2 - (242*(pow2(
        lmMst1) + pow2(lmMst2)))/21.)*pow6(Mst1))/pow6(Mst2)))/3.)/pow2(Mst1) +
        ((-4*OepS2*(20*Dmsqst2*(45*pow2(Mst1)*pow2(Mst2) + 53*pow4(Mst1) + 27*
        pow4(Mst2)) + Mst2*(5205*pow2(Mst1)*pow3(Mst2) + 12062*Mst2*pow4(Mst1)
        - 24*Dmglst2*(150*pow2(Mst1)*pow2(Mst2) + 919*pow4(Mst1)) + 1431*pow5(
        Mst2))))/729. + (S2*(100*Dmsqst2*(9*(127 + 30*lmMst1 - 30*lmMst2)*pow2(
        Mst1)*pow2(Mst2) + (1283 + 318*lmMst1 - 318*lmMst2)*pow4(Mst1) + 81*(9
        + 2*lmMst1 - 2*lmMst2)*pow4(Mst2)) + Mst2*(9*(17113 + 17350*lmMst1 -
        17350*lmMst2)*pow2(Mst1)*pow3(Mst2) + 4*(138286 + 90465*lmMst1 - 90465*
        lmMst2)*Mst2*pow4(Mst1) - 24*Dmglst2*(45*(139 + 100*lmMst1 - 100*
        lmMst2)*pow2(Mst1)*pow2(Mst2) + (34538 + 27570*lmMst1 - 27570*lmMst2)*
        pow4(Mst1) + 5697*pow4(Mst2)) + 81*(-141 + 530*lmMst1 - 530*lmMst2)*
        pow5(Mst2))))/270.)/pow6(Mst2))) + Mt*(((1 - 2*lmMst2)*s2t*pow2(Mst2)*(
        -2*shiftst3*pow2(Mt) + 5*shiftst1*pow2(Mst2)*pow2(s2t)))/(3.*pow2(Mst1)
        ) + pow3(s2t)*((pow2(Mst2)*(21005 + 1440*B4 - 144*D3 + 72*DN + 21216*
        lmMst1 + 1908*pow2(lmMst1) - 12*lmMst2*(2950 - 2040*lmMst1 + 117*pow2(
        lmMst1)) + 108*(-283 + 98*lmMst1)*pow2(lmMst2) + (1080*Dmsqst2*(-1 + 2*
        lmMst2))/pow2(Mst1) - (108*(-1 + 2*lmMst2)*shiftst3*((-1 - lmMst1 +
        lmMst2)*pow2(Mst1) + pow2(Mst2)))/pow2(Mst1) + 576*pow3(lmMst1) - 9756*
        pow3(lmMst2)))/648. + pow2(Mst1)*(190.4424279835391 + 2*B4 + (2*D3)/9.
         - (2*DN)/
        9. + (22004*lmMst1)/405. - (736*pow2(lmMst1))/27. - (lmMst2*(12148 -
        76560*lmMst1 + 12105*pow2(lmMst1)))/810. + (5*(-583 + 438*lmMst1)*pow2(
        lmMst2))/54. - (55*pow3(lmMst1))/27. - (1273*pow3(lmMst2))/54.) - ((
        Mst2*(119 + 218*lmMst2 + 32*lmMst1*(1 + lmMst2) + 107*pow2(lmMst2)) +
        Dmglst2*(-40 + (460 + 64*lmMst1)*lmMst2 + 310*pow2(lmMst2)))*pow3(Mst2)
        )/(18.*pow2(Mst1)) + ((Dmsqst2*(97294 + 10485*lmMst1 + 45*(-383 + 330*
        lmMst1)*lmMst2 - 7425*(pow2(lmMst1) + pow2(lmMst2)))*pow2(Mst1))/2025.
         + (
        96.95148419197191 + (58499851*lmMst2)/793800. - ((76483 + 26985*lmMst2)
        *pow2(lmMst1))/945. - (149081*pow2(lmMst2))/1890. + lmMst1*(-
        77.33484630889393 + (302047*lmMst2)/1890. + 61*pow2(lmMst2)) - (35*
        pow3(lmMst1))/27. - (841*pow3(lmMst2))/27.)*pow4(Mst1))/pow2(Mst2) -
        Dmglst2*(Mst2*(46.25617283950617 - (4*B4)/3. + (16*D3)/9. - (10*DN)/9.
         + (44*lmMst1)/
        3. - (13*pow2(lmMst1))/3. + (4*lmMst2*(-81 - 124*lmMst1 + 8*pow2(
        lmMst1)))/9. + (48.77777777777778 - (82*lmMst1)/3.)*pow2(lmMst2) + (20*
        Dmsqst2)/pow2(Mst1) + (214*pow3(lmMst2))/9.) + (810000*Dmsqst2*(5 + 2*
        lmMst1 - 2*lmMst2) + pow2(Mst1)*(5430043 + 524580*lmMst2 + 900*(-859 +
        3690*lmMst2)*pow2(lmMst1) + 1454400*pow2(lmMst2) - 60*lmMst1*(51493 +
        24630*lmMst2 + 112950*pow2(lmMst2)) + 45000*pow3(lmMst1) + 3411000*
        pow3(lmMst2)))/(121500.*Mst2) + ((203.8689796251638 + (2171669*lmMst2)/
        119070. + ((433 + 6552*lmMst2)*pow2(lmMst1))/189. + (1987*pow2(lmMst2))
        /189. - (lmMst1*(2740559 + 1524600*lmMst2 + 7514640*pow2(lmMst2)))/
        119070. - (56*pow3(lmMst1))/27. + (824*pow3(lmMst2))/27.)*pow4(Mst1))/
        pow3(Mst2)) + Dmsqst2*(54.49074074074074 + (140*lmMst1)/9. - (20*(7 +
        3*lmMst1)*lmMst2)/9. + (20*pow2(lmMst2))/3. + (5*Dmglst2*(55 - 34*
        lmMst1 + 34*lmMst2)*pow2(Mst1))/(3.*pow3(Mst2)) + (((5*Dmglst2*(419 -
        57*lmMst1 + 57*lmMst2))/27. + Mst2*(17.77343371446168 - (452768*lmMst1)
        /19845. + (22.815217939027463 + (122*lmMst1)/7.)*lmMst2 - (61*(pow2(
        lmMst1) + pow2(lmMst2)))/7.))*pow4(Mst1))/pow5(Mst2)) + (8*(1 + lmMst2)
        *(4*Dmglst2*lmMst2 + Mst2 + lmMst2*Mst2)*pow5(Mst2))/(9.*pow4(Mst1)) +
        (-30*Dmsqst2*(18*(40*OepS2 - 27*(59 + 30*lmMst1 - 30*lmMst2)*S2)*pow2(
        Mst1)*pow2(Mst2) + 4*(208*OepS2 - 27*(347 + 156*lmMst1 - 156*lmMst2)*
        S2)*pow4(Mst1) + 27*(8*OepS2 - 81*(-15 + 2*lmMst1 - 2*lmMst2)*S2)*pow4(
        Mst2)) + Mst2*(54*Dmglst2*(632*OepS2 + 9*(16193 - 1422*lmMst1 + 1422*
        lmMst2)*S2)*pow2(Mst1)*pow2(Mst2) - 180*(340*OepS2 - 81*(424 + 85*
        lmMst1 - 85*lmMst2)*S2)*pow2(Mst1)*pow3(Mst2) - 105*Mst2*(788*OepS2 -
        27*(2662 + 591*lmMst1 - 591*lmMst2)*S2)*pow4(Mst1) + 4*Dmglst2*(25328*
        OepS2 + 27*(47051 - 18996*lmMst1 + 18996*lmMst2)*S2)*pow4(Mst1) + 54*
        Dmglst2*(56*OepS2 - 81*(-1677 + 14*lmMst1 - 14*lmMst2)*S2)*pow4(Mst2) +
        81*(-184*OepS2 + 81*(307 + 46*lmMst1 - 46*lmMst2)*S2)*pow5(Mst2)))/(
        4374.*pow4(Mst2))) + (T1ep*(81*pow4(Mst2)*(-40*Dmsqst2*s2t*pow2(Mt) -
        21*Mt*pow2(s2t)*pow3(Mst2) + 28*Mst2*pow3(Mt) + pow2(Mst2)*(-106*s2t*
        pow2(Mt) + 20*Dmsqst2*pow3(s2t)) + 46*pow3(s2t)*pow4(Mst2)) + 9*pow2(
        Mst1)*pow2(Mst2)*(-600*Dmsqst2*s2t*pow2(Mt) - 627*Mt*pow2(s2t)*pow3(
        Mst2) + 1760*Mst2*pow3(Mt) + pow2(Mst2)*(-3470*s2t*pow2(Mt) + 600*
        Dmsqst2*pow3(s2t)) + 1700*pow3(s2t)*pow4(Mst2)) + 3*pow4(Mst1)*(-2120*
        Dmsqst2*s2t*pow2(Mt) - 3198*Mt*pow2(s2t)*pow3(Mst2) + 16240*Mst2*pow3(
        Mt) - 4*pow2(Mst2)*(6031*s2t*pow2(Mt) - 520*Dmsqst2*pow3(s2t)) + 6895*
        pow3(s2t)*pow4(Mst2)) + Dmglst2*(-27*pow2(Mst1)*pow2(Mst2)*(-800*Mst2*
        s2t*pow2(Mt) - 907*Mt*pow2(Mst2)*pow2(s2t) + 1984*pow3(Mt) + 316*pow3(
        Mst2)*pow3(s2t)) + (132336*Mst2*s2t*pow2(Mt) + 71202*Mt*pow2(Mst2)*
        pow2(s2t) - 259408*pow3(Mt) - 25328*pow3(Mst2)*pow3(s2t))*pow4(Mst1) -
        189*(20*pow3(Mt)*pow4(Mst2) - 15*Mt*pow2(s2t)*pow6(Mst2) + 4*pow3(s2t)*
        pow7(Mst2)))))/(729.*pow6(Mst2))) + ((-5*Mt*s2t*shiftst1*((1 - 2*
        lmMst2)*pow2(Mst2)*(4*pow2(Mst2)*pow2(Mt) - 2*pow2(Mst1)*(2*pow2(Mt) +
        (-lmMst1 + lmMst2)*pow2(Mst2)*pow2(s2t)) + pow2(s2t)*pow4(Mst1)) +
        Dmsqst2*(3 - 2*lmMst2)*(-4*pow2(Mst1)*pow2(Mt) + 4*pow2(Mst2)*pow2(Mt)
        + 2*(lmMst1 - lmMst2)*pow2(Mst1)*pow2(Mst2)*pow2(s2t) + pow2(s2t)*(
        pow4(Mst1) - pow4(Mst2)))))/(3.*pow2(Mst2)) + (pow2(Mt)*(4*pow2(Mt)*(-(
        Dmglst2*(216*pow2(Mst2)*(97837 + 71580*lmMt - 9920*OepS2 + 43272*S2 -
        360*(23 + 8*lmMst2 - 4*lmMt)*pow2(lmMst1) + 720*(97 + 2*lmMt)*pow2(
        lmMst2) - 8640*pow2(lmMt) - 30*lmMst2*(-2911 + 396*lmMt + 6696*S2 +
        144*pow2(lmMt)) + 10*lmMst1*(-6157 + 642*lmMt - 6*lmMst2*(791 + 48*
        lmMt) + 20088*S2 + 882*pow2(lmMst2) + 432*pow2(lmMt)) - 5940*pow3(
        lmMst2))*pow4(Mst1) - 129600*Dmsqst2*((101 - 90*lmMst1 + 177*lmMst2 -
        87*lmMt)*pow2(Mst1)*pow2(Mst2) + (497 - 306*lmMst1 + 609*lmMst2 - 303*
        lmMt)*pow4(Mst1)) + 135*pow2(Mst1)*(57495 + 45408*lmMt - 1120*OepS2 -
        43740*S2 - 8*lmMst2*(-6952 + 588*lmMt + 2835*S2) + 2304*(-1 + lmMst2)*
        pow2(lmMst1) + 96*(79 + 48*lmMt)*pow2(lmMst2) + 72*lmMst1*(-88 + 64*
        lmMt - 8*lmMst2*(9 + 8*lmMt) + 315*S2 + 164*pow2(lmMst2)) - 14112*pow3(
        lmMst2))*pow4(Mst2) + 20*(5659217 + 1592460*lmMt - 518816*OepS2 +
        9976392*S2 - 972*(569 + 180*lmMst2 - 126*lmMt)*pow2(lmMst1) + 324*(5353
        + 126*lmMt)*pow2(lmMst2) - 186624*pow2(lmMt) + 72*lmMst1*(-27653 -
        3015*lmMt - 18*lmMst2*(689 + 126*lmMt) + 145917*S2 + 2160*pow2(lmMst2)
        + 2592*pow2(lmMt)) - 36*lmMst2*(-39031 - 3204*lmMt + 291834*S2 + 5184*
        pow2(lmMt)) + 1944*pow3(lmMst1) + 17496*pow3(lmMst2))*pow6(Mst1) -
        622080*(-1 + 2*lmMst2 + 3*pow2(lmMst2))*pow6(Mst2))) + 3*Mst2*(48*pow2(
        Mst2)*(110219 - 1080*lmMt - 4400*OepS2 + 74034*S2 + 540*(-15 + 4*lmMt)*
        pow2(lmMst1) + 810*(121 - 8*lmMt)*pow2(lmMst2) - 135*lmMst1*(337 +
        lmMst2*(588 - 32*lmMt) - 132*lmMt - 660*S2 + 194*pow2(lmMst2) - 48*
        pow2(lmMt)) - 6480*pow2(lmMt) - 405*lmMst2*(31 + 54*lmMt + 220*S2 + 16*
        pow2(lmMt)) + 26190*pow3(lmMst2))*pow4(Mst1) + 388800*Dmsqst2*((5 + 2*
        lmMst1 - 5*lmMst2 + 3*lmMt)*pow2(Mst1)*pow2(Mst2) + (1 + 6*lmMst1 - 13*
        lmMst2 + 7*lmMt)*pow4(Mst1)) + 27*pow2(Mst1)*(56439 - 24000*lmMt -
        1120*OepS2 - 25596*S2 - 360*lmMst2*(-12 + 44*lmMt + 63*S2) - 3840*(1 +
        lmMst2)*pow2(lmMst1) + 480*(163 - 16*lmMt)*pow2(lmMst2) - 120*lmMst1*(
        184 + 8*lmMst2*(57 - 8*lmMt) - 64*lmMt - 189*S2 + 164*pow2(lmMst2)) -
        11520*pow2(lmMt) + 23520*pow3(lmMst2))*pow4(Mst2) + 20*(678923 + 19440*
        lmMt - 32480*OepS2 + 881712*S2 + 540*(661 - 30*lmMt)*pow2(lmMst2) -
        15552*pow2(lmMt) - 108*((457 - 12*lmMst2 - 126*lmMt)*pow2(lmMst1) +
        lmMst1*(1841 + lmMst2*(2626 - 24*lmMt) - 420*lmMt - 6090*S2 + 840*pow2(
        lmMst2) - 288*pow2(lmMt)) + lmMst2*(619 + 498*lmMt + 6090*S2 + 288*
        pow2(lmMt))) + 216*pow3(lmMst1) + 89208*pow3(lmMst2))*pow6(Mst1) +
        207360*pow2(1 + lmMst2)*pow6(Mst2))) - 3*pow2(s2t)*(Dmglst2*(9*(195293
        + 639360*B4 - 8640*DN - 709200*lmMst1 + 145120*OepS2 - 108*(23989 +
        27210*lmMst1 - 27210*lmMst2)*S2 + 60480*pow2(lmMst1) + 720*lmMst2*(2149
        - 1296*lmMst1 + 96*pow2(lmMst1)) - 8640*(-101 + 105*lmMst1)*pow2(
        lmMst2) + 838080*pow3(lmMst2))*pow4(Mst1)*pow4(Mst2) - 2*pow2(Mst2)*(
        15069803 - 2877120*B4 + 38880*DN + 6325200*lmMst1 - 1898720*OepS2 +
        108*(525961 + 356010*lmMst1 - 356010*lmMst2)*S2 + 447120*pow2(lmMst1) -
        360*lmMst2*(28667 - 5238*lmMst1 + 3024*pow2(lmMst1)) + 38880*(-60 +
        157*lmMst1)*pow2(lmMst2) - 155520*pow3(lmMst1) - 4860000*pow3(lmMst2))*
        pow6(Mst1) + 583200*Dmsqst2*(40*(1 - lmMst1 + lmMst2)*pow2(Mst2)*pow4(
        Mst1) - 4*(2 + 5*lmMst1 - 5*lmMst2)*pow2(Mst1)*pow4(Mst2) + (80 - 43*
        lmMst1 + 43*lmMst2)*pow6(Mst1)) + 27*pow2(Mst1)*(46957 + 188640*B4 -
        3600*DN - 37440*lmMst1 + 5600*OepS2 - 324*(-453 + 350*lmMst1 - 350*
        lmMst2)*S2 + 11520*pow2(lmMst1) - 2880*lmMst2*(-221 + 55*lmMst1 + 4*
        pow2(lmMst1)) - 1440*(-232 + 121*lmMst1)*pow2(lmMst2) + 185760*pow3(
        lmMst2))*pow6(Mst2) + 622080*(-1 + 2*lmMst2 + 3*pow2(lmMst2))*pow8(
        Mst2)) + 15*(3*(346405 + 62208*B4 + 246672*lmMst2 - 6688*OepS2 - 324*(-
        685 + 418*lmMst2)*S2 + 1728*(-7 + 8*lmMst2)*pow2(lmMst1) + 323136*pow2(
        lmMst2) - 216*lmMst1*(990 + 1440*lmMst2 - 627*S2 + 584*pow2(lmMst2)) +
        112320*pow3(lmMst2))*pow4(Mst1)*pow5(Mst2) + 2*(795601 + 93312*B4 +
        365040*lmMst2 - 17056*OepS2 + 663444*S2 - 345384*lmMst2*S2 + 432*(163 +
        264*lmMst2)*pow2(lmMst1) + 592704*pow2(lmMst2) - 216*lmMst1*(1609 +
        3070*lmMst2 - 1599*S2 + 1692*pow2(lmMst2)) - 3456*pow3(lmMst1) +
        254880*pow3(lmMst2))*pow3(Mst2)*pow6(Mst1) + pow2(Mst1)*(-38880*
        Dmsqst2*Mst2*(8*(-lmMst1 + lmMst2)*pow2(Mst1)*pow2(Mst2) + (4 - 9*
        lmMst1 + 9*lmMst2)*pow4(Mst1) - 4*(1 + lmMst1 - lmMst2)*pow4(Mst2)) +
        27*(27147 + 5280*B4 - 48*DN - 5952*lmMst1 - 224*OepS2 + 324*(-1 + 14*
        lmMst1 - 14*lmMst2)*S2 - 768*pow2(lmMst1) - 192*lmMst2*(-198 + 73*
        lmMst1 + 4*pow2(lmMst1)) - 288*(-84 + 19*lmMst1)*pow2(lmMst2) + 6240*
        pow3(lmMst2))*pow7(Mst2)) + 41472*pow2(1 + lmMst2)*pow9(Mst2)))))/(
        174960.*pow6(Mst2)))/pow2(Mst1)) + (3*xDR2DRMOD*((-36*Mt*s2t*pow2(Mst1)
        *(256*Dmglst2*MuSUSY*(-1 + 2*lmMst2 + 3*pow2(lmMst2))*pow2(Mt) + 135*
        shiftst1*Tbeta*pow2(s2t)*pow4(Mst2)))/Tbeta + (pow2(Mst1)*(-4860*
        shiftst1*pow4(Mst2)*(4*s2t*(pow2(Mst1) - pow2(Mst2))*(Dmsqst2 + pow2(
        Mst2))*pow3(Mt) - Mt*pow3(s2t)*(2*Dmsqst2*(lmMst1 - lmMst2)*pow2(Mst1)*
        pow2(Mst2) + (Dmsqst2 + pow2(Mst2))*pow4(Mst1) - (Dmsqst2 + 2*(-lmMst1
        + lmMst2)*pow2(Mst1))*pow4(Mst2))) + 512*pow4(Mt)*(9*(1 + lmMst2)*Mst2*
        ((3 - 21*lmMst2 + 2*lmMst1*(8 + 3*lmMst2 - 3*lmMt) + 5*lmMt + 6*lmMst2*
        lmMt - 6*pow2(lmMst2))*pow2(Mst2)*pow4(Mst1) + (-1 - 7*lmMst2 + 2*
        lmMst1*(2 + lmMst2 - lmMt) + 3*lmMt + 2*lmMst2*lmMt - 2*pow2(lmMst2))*
        pow2(Mst1)*pow4(Mst2) + (12 - 45*lmMst2 + 2*lmMst1*(19 + 6*lmMst2 - 6*
        lmMt) + 7*lmMt + 12*lmMst2*lmMt - 12*pow2(lmMst2))*pow6(Mst1) - 2*(1 +
        lmMst2)*pow6(Mst2)) + Dmglst2*(pow2(Mst2)*(263 + lmMst2*(962 - 213*
        lmMt) - 177*lmMt + (393 - 18*lmMt)*pow2(lmMst2) - 18*lmMst1*(26 -
        lmMst2*(-17 + lmMt) - 7*lmMt + pow2(lmMst2)) + 18*pow3(lmMst2))*pow4(
        Mst1) - pow2(Mst1)*(17*(-7 + 3*lmMt) + lmMst2*(-224 + 15*lmMt) - 3*(5 +
        6*lmMt)*pow2(lmMst2) - 18*lmMst1*(-4 + lmMt - lmMst2*(1 + lmMt) + pow2(
        lmMst2)) + 18*pow3(lmMst2))*pow4(Mst2) + (290 + lmMst2*(2447 - 645*
        lmMt) - 375*lmMt - 3*(-509 + 60*lmMt)*pow2(lmMst2) - 18*lmMst1*(87 +
        lmMst2*(71 - 10*lmMt) - 22*lmMt + 10*pow2(lmMst2)) + 180*pow3(lmMst2))*
        pow6(Mst1) - 18*(-1 + 2*lmMst2 + 3*pow2(lmMst2))*pow6(Mst2))) - 3456*
        pow2(Mst2)*pow2(Mt)*pow2(s2t)*((1 + lmMst2)*Mst2*(2*(2 - 5*lmMst2 +
        lmMst1*(5 + 2*lmMst2) - 2*pow2(lmMst2))*pow2(Mst2)*pow4(Mst1) + 2*(-2*
        lmMst2*(2 + lmMst2) + lmMst1*(3 + 2*lmMst2))*pow2(Mst1)*pow4(Mst2) + (5
        - 12*lmMst2 + 4*lmMst1*(3 + lmMst2) - 4*pow2(lmMst2))*pow6(Mst1) - 2*(1
        + lmMst2)*pow6(Mst2)) + Dmglst2*(2*pow2(Mst2)*(8 + 19*lmMst2 - 5*pow2(
        lmMst2) + lmMst1*(-7 + 5*lmMst2 + 6*pow2(lmMst2)) - 6*pow3(lmMst2))*
        pow4(Mst1) + 2*pow2(Mst1)*(7 + 9*lmMst2 - 8*pow2(lmMst2) + lmMst1*(-3 +
        5*lmMst2 + 6*pow2(lmMst2)) - 6*pow3(lmMst2))*pow4(Mst2) + (17 + 51*
        lmMst2 - 4*pow2(lmMst2) + 4*lmMst1*(-6 + lmMst2 + 3*pow2(lmMst2)) - 12*
        pow3(lmMst2))*pow6(Mst1) - 2*(-1 + 2*lmMst2 + 3*pow2(lmMst2))*pow6(
        Mst2)))) + 9*Mt*pow3(s2t)*(pow2(Mst1)*(45*(2*Dmglst2*Mst2 - 2*lmMst1*
        pow2(Mst1) + 2*lmMst2*pow2(Mst1) + pow2(Mst2))*pow4(Dmsqst2) + 60*
        Dmsqst2*(2*Dmglst2 - 5*Mst2)*pow3(Mst2)*(2*(lmMst1 - lmMst2)*pow2(Mst1)
        *pow2(Mst2) + pow4(Mst1) - pow4(Mst2))) + 4*pow3(Mst2)*(-2*(8 - 237*
        lmMst2 - 16*(1 + lmMst2)*pow2(lmMst1) - 308*pow2(lmMst2) + lmMst1*(269
        + 348*lmMst2 + 123*pow2(lmMst2)) - 107*pow3(lmMst2))*pow4(Mst1)*pow5(
        Mst2) - (125 + lmMst2*(156 + 512*lmMst1*(1 + lmMst2)) - 181*pow2(
        lmMst2) - 256*((1 + lmMst2)*pow2(lmMst1) + pow3(lmMst2)))*pow3(Mst2)*
        pow6(Mst1) + (221 + 284*lmMst2 + 32*lmMst1*(1 + lmMst2) + 107*pow2(
        lmMst2))*pow2(Mst1)*pow7(Mst2) + 16*(1 + lmMst2)*Mst2*(1 + lmMst1*(2 -
        32*lmMst2) - 2*lmMst2 + 16*(pow2(lmMst1) + pow2(lmMst2)))*pow8(Mst1) +
        2*Dmglst2*(4*pow2(Mst1)*pow2(Mst2)*(20*pow2(Mst1)*pow2(Mst2) + 9*pow4(
        Mst1) - 5*pow4(Mst2)) - 2*lmMst1*(-4 + 246*lmMst2 + 123*pow2(lmMst2))*
        pow4(Mst1)*pow4(Mst2) + 214*pow3(lmMst2)*pow4(Mst1)*pow4(Mst2) + 32*
        lmMst2*pow2(lmMst1)*pow4(Mst1)*(8*pow2(Mst1)*pow2(Mst2) + 8*pow4(Mst1)
        + pow4(Mst2)) + 256*(pow2(Mst1) + pow2(Mst2))*pow3(lmMst2)*pow6(Mst1) +
        32*lmMst1*lmMst2*pow2(Mst1)*pow6(Mst2) - 32*lmMst1*((3 + 3*lmMst2 + 16*
        pow2(lmMst2))*pow2(Mst2)*pow6(Mst1) + (3 + 2*lmMst2 + 16*pow2(lmMst2))*
        pow8(Mst1)) + pow2(lmMst2)*(396*pow4(Mst1)*pow4(Mst2) + 37*pow2(Mst2)*
        pow6(Mst1) + 155*pow2(Mst1)*pow6(Mst2) + 64*pow8(Mst1) - 32*pow8(Mst2))
        + 2*lmMst2*(4*pow4(Mst1)*pow4(Mst2) + 21*pow2(Mst2)*pow6(Mst1) + 115*
        pow2(Mst1)*pow6(Mst2) + 56*pow8(Mst1) - 16*pow8(Mst2))) - 16*pow2(1 +
        lmMst2)*pow9(Mst2))) + 6*s2t*pow2(Mt)*(-90*Dmsqst2*Mt*pow2(Mst1)*(3*
        pow3(Dmsqst2) + 4*(2*Dmglst2 - 5*Mst2)*(pow2(Mst1) - pow2(Mst2))*pow3(
        Mst2)) + 8*Mst2*Mt*(3*(109 + 76*lmMst2 - 21*pow2(lmMst2) + 64*lmMst1*
        pow2(1 + lmMst2) - 32*((1 + lmMst2)*pow2(lmMst1) + pow3(lmMst2)))*pow4(
        Mst1)*pow5(Mst2) - 3*(173 + 188*lmMst2 + 32*lmMst1*(1 + lmMst2) + 59*
        pow2(lmMst2))*pow2(Mst1)*pow7(Mst2) - 96*(1 + lmMst2)*Mst2*(pow2(1 -
        lmMst1 + lmMst2)*pow2(Mst2)*pow6(Mst1) + (1 + 3*lmMst2 - lmMst1*(3 + 2*
        lmMst2) + pow2(lmMst1) + pow2(lmMst2))*pow8(Mst1)) + 2*Dmglst2*((52 +
        370*lmMst2 - 96*lmMst2*pow2(lmMst1) + 81*pow2(lmMst2) + 96*lmMst1*(-1 +
        lmMst2 + 2*pow2(lmMst2)) - 96*pow3(lmMst2))*pow4(Mst1)*pow4(Mst2) - 96*
        (lmMst2*pow2(lmMst1) + lmMst1*(4 + 2*lmMst2 - 2*pow2(lmMst2)) + (-4 +
        lmMst2)*pow2(1 + lmMst2))*pow2(Mst2)*pow6(Mst1) - 3*(-52 + 2*(51 + 16*
        lmMst1)*lmMst2 + 59*pow2(lmMst2))*pow2(Mst1)*pow6(Mst2) - 96*(-6 - 14*
        lmMst2 + lmMst2*pow2(lmMst1) + lmMst1*(9 + 6*lmMst2 - 2*pow2(lmMst2)) -
        6*pow2(lmMst2) + pow3(lmMst2))*pow8(Mst1) + 96*lmMst2*(1 + lmMst2)*
        pow8(Mst2)) + 48*pow2(1 + lmMst2)*pow9(Mst2)) + (3*MuSUSY*(-4*(128*(1 +
        lmMst2)*Mst2*Mt*(1 + 5*lmMst2 - lmMst1*(3 + 2*lmMst2) + 2*pow2(lmMst2))
        + s2t*(15*Dmsqst2*(1 - 2*lmMst1 + 2*lmMst2)*(5 - 9*shiftst1) + pow2(
        Mst2)*(189 + 726*lmMst2 - 135*(1 + 2*lmMst2)*shiftst1 + 32*(1 + lmMst2)
        *pow2(lmMst1) + 707*pow2(lmMst2) - 2*lmMst1*(253 + 332*lmMst2 - 135*
        shiftst1 + 123*pow2(lmMst2)) + 214*pow3(lmMst2))))*pow4(Mst1)*pow4(
        Mst2) - pow2(Mst1)*(512*Mt*pow2(1 + lmMst2)*pow7(Mst2) + s2t*(45*pow4(
        Dmsqst2) + 60*Dmsqst2*(5 - 9*shiftst1)*pow6(Mst2) + 4*(205 + 252*lmMst2
        + 32*lmMst1*(1 + lmMst2) - 135*shiftst1 + 91*pow2(lmMst2))*pow8(Mst2)))
        + 8*(pow2(Mst2)*(64*(1 + lmMst2)*Mst2*Mt*(1 - 10*lmMst2 + 4*lmMst1*(2 +
        lmMst2) - 4*pow2(lmMst2)) + s2t*(15*Dmsqst2*(lmMst1 - lmMst2)*(5 - 9*
        shiftst1) - pow2(Mst2)*(32 + 15*lmMst2*(19 - 9*shiftst1) + 144*(1 +
        lmMst2)*pow2(lmMst1) + 444*pow2(lmMst2) - lmMst1*(253 + 588*lmMst2 -
        135*shiftst1 + 379*pow2(lmMst2)) + 235*pow3(lmMst2))))*pow6(Mst1) + (-
        15*Dmsqst2*(lmMst1 - lmMst2)*s2t*(-5 + 9*shiftst1) + 32*(1 + lmMst2)*
        Mst2*Mt*(7 - 32*lmMst2 + 4*lmMst1*(7 + 3*lmMst2) - 12*pow2(lmMst2)) -
        s2t*pow2(Mst2)*(40 + 277*lmMst2 - 135*lmMst2*shiftst1 + 272*(1 +
        lmMst2)*pow2(lmMst1) + 556*pow2(lmMst2) - lmMst1*(237 + 828*lmMst2 -
        135*shiftst1 + 635*pow2(lmMst2)) + 363*pow3(lmMst2)))*pow8(Mst1) +
        Dmglst2*((15*Dmsqst2*(1 - 2*lmMst1 + 2*lmMst2)*s2t + 64*Mst2*Mt*(8 + 7*
        lmMst2 - 11*pow2(lmMst2) + lmMst1*(-3 + 5*lmMst2 + 6*pow2(lmMst2)) - 6*
        pow3(lmMst2)) - s2t*pow2(Mst2)*(60 + 206*lmMst2 + 32*lmMst2*pow2(
        lmMst1) + lmMst1*(8 - 460*lmMst2 - 246*pow2(lmMst2)) + 519*pow2(lmMst2)
        + 214*pow3(lmMst2)))*pow3(Mst2)*pow4(Mst1) + 2*(Mst2*(15*Dmsqst2*(-
        lmMst1 + lmMst2)*s2t + 64*Mst2*Mt*(8 + 13*lmMst2 - 8*pow2(lmMst2) +
        lmMst1*(-5 + 5*lmMst2 + 6*pow2(lmMst2)) - 6*pow3(lmMst2)) - s2t*pow2(
        Mst2)*(48 + 124*lmMst2 + 144*lmMst2*pow2(lmMst1) + 278*pow2(lmMst2) -
        lmMst1*(44 + 278*lmMst2 + 379*pow2(lmMst2)) + 235*pow3(lmMst2)))*pow6(
        Mst1) + (16*Mt*(49 + 103*lmMst2 - 36*(1 + lmMst2)*pow2(lmMst2) + 4*
        lmMst1*(-11 + 6*lmMst2 + 9*pow2(lmMst2))) - Mst2*s2t*(48 + 4*lmMst2*(45
        + 68*pow2(lmMst1)) + 310*pow2(lmMst2) - lmMst1*(92 + 310*lmMst2 + 635*
        pow2(lmMst2)) + 363*pow3(lmMst2)))*pow8(Mst1)) + s2t*(pow2(Mst1)*(15*
        Dmsqst2 - (-20 + 2*(99 + 16*lmMst1)*lmMst2 + 123*pow2(lmMst2))*pow2(
        Mst2))*pow5(Mst2) + 32*lmMst2*(1 + lmMst2)*pow9(Mst2)))) + 64*s2t*pow2(
        1 + lmMst2)*power10(Mst2)))/Tbeta))/pow6(Mst2)))/pow4(Mst1))))/11664.;
}

/**
 * @return returns the susy log^0 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double H6b::calc_coef_at_as2_no_sm_logs_log0() const {

   const double result =
      (30*pow2(Dmsqst2)*pow2(Mst1)*(166698000*(-2 + z2)*pow2(-4*pow2(Mt) +
        pow2(Mst2)*pow2(s2t))*pow4(Mst2) - 385875*pow2(Mst1)*pow2(Mst2)*(-8*
        pow2(Mst2)*pow2(Mt)*pow2(s2t)*(-160 + 64*OepS2 - 648*S2 - 96*T1ep -
        2076*z2 + 903*z3 - 48*z4 - 72*pow2(z2)) + 3456*Mst2*s2t*(-15 + 20*z2 -
        16*z3)*pow3(Mt) + 6912*Mt*(-3 + 3*z2 - z3)*pow3(Mst2)*pow3(s2t) + (
        99040 - 89856*z2 + 32400*z3)*pow4(Mt) + (2122 + 128*OepS2 - 55728*S2 -
        192*T1ep - 4152*z2 + 13209*z3 - 96*z4 - 144*pow2(z2))*pow4(Mst2)*pow4(
        s2t)) + 2058*pow4(Mst1)*(4*pow2(Mst2)*pow2(Mt)*pow2(s2t)*(487028 +
        40000*OepS2 - 1917000*S2 - 60000*T1ep - 1081500*z2 - 2625*z3 - 30000*z4
        - 45000*pow2(z2)) - 2592000*Mst2*s2t*(-11 + 14*z2 - 8*z3)*pow3(Mt) +
        1296000*Mt*(-2 + z2 + z3)*pow3(Mst2)*pow3(s2t) + 288*(-176473 + 202500*
        z2 - 90000*z3)*pow4(Mt) + (1932736 - 64000*OepS2 + 3294000*S2 + 96000*
        T1ep + 726000*z2 - 710625*z3 + 48000*z4 + 72000*pow2(z2))*pow4(Mst2)*
        pow4(s2t)) - 5*pow2(s2t)*(-100018800*Mst2*Mt*s2t*(-23 + 30*z2 - 20*z3)
        + pow2(Mst2)*pow2(s2t)*(76781738 + 15366400*OepS2 - 837194400*S2 -
        23049600*T1ep - 207652200*z2 - 63103425*z3 - 11524800*z4 - 17287200*
        pow2(z2)) - 600*pow2(Mt)*(109760*OepS2 - 3*(1281429 + 1720488*S2 +
        54880*T1ep + 44590*z2 - 303212*z3 + 27440*z4 + 41160*pow2(z2))))*pow6(
        Mst1)) + 120*Dmsqst2*pow2(Mst1)*(3087*Mst2*pow4(Mst1)*(-4500*Dmglst2*(
        8*(113 - 48*z2)*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 192*Mst2*s2t*(51*z2 -
        28*(2 + z3))*pow3(Mt) + 96*Mt*(5*z2 + 3*(-4 + z3))*pow3(Mst2)*pow3(s2t)
        + 48*(-401 + 308*z2 - 144*z3)*pow4(Mt) - 3*(-75 + 32*z2 + 8*z3)*pow4(
        Mst2)*pow4(s2t)) + Mst2*(-8*pow2(Mst2)*pow2(Mt)*pow2(s2t)*(36449 +
        4000*OepS2 - 310500*S2 - 6000*T1ep - 75750*z2 + 21000*z3 - 3000*z4 -
        4500*pow2(z2)) - 864000*Mst2*s2t*(-6 + 5*z2 - 4*z3)*pow3(Mt) - 432000*
        Mt*(1 + z2 - z3)*pow3(Mst2)*pow3(s2t) + 96*(-91321 + 54000*z2 - 36000*
        z3)*pow4(Mt) + (52199 + 28000*OepS2 - 3415500*S2 - 42000*T1ep - 557250*
        z2 + 259500*z3 - 21000*z4 - 31500*pow2(z2))*pow4(Mst2)*pow4(s2t))) +
        385875*pow2(Mst1)*pow3(Mst2)*(-16*Dmglst2*(36*(23 - 24*z2)*pow2(Mst2)*
        pow2(Mt)*pow2(s2t) + 24*Mst2*s2t*(289 - 342*z2 + 216*z3)*pow3(Mt) -
        216*Mt*(-2 + 2*z2 - z3)*pow3(Mst2)*pow3(s2t) + 8*(-1477 + 1242*z2 -
        540*z3)*pow4(Mt) + 27*(-2 + 4*z2 - 3*z3)*pow4(Mst2)*pow4(s2t)) + 3*
        Mst2*(-4*pow2(Mst2)*pow2(Mt)*pow2(s2t)*(32*OepS2 - 3*(391 + 972*S2 +
        16*T1ep + 154*z2 - 232*z3 + 8*z4 + 12*pow2(z2))) - 1152*Mst2*s2t*(7 +
        2*z2 - 8*z3)*pow3(Mt) - 1152*Mt*(-1 + 2*z2 - 3*z3)*pow3(Mst2)*pow3(s2t)
        + 64*(-103 + 108*z2 - 72*z3)*pow4(Mt) + (-1213 + 32*OepS2 + 4860*S2 -
        48*T1ep - 462*z2 - 276*z3 - 24*z4 - 36*pow2(z2))*pow4(Mst2)*pow4(s2t)))
        + 41674500*(Mst2 - 4*Dmglst2*(-3 + z2) + 2*Mst2*z2)*pow2(-4*pow2(Mt) +
        pow2(Mst2)*pow2(s2t))*pow5(Mst2) + 2*(-4*Mst2*pow2(Mt)*pow2(s2t)*(
        1157625*Dmglst2*(-1931 + 1098*z2 + 180*z3) + Mst2*(-22430251 + 5488000*
        OepS2 - 324135000*S2 - 8232000*T1ep - 103929000*z2 + 28812000*z3 -
        4116000*z4 - 6174000*pow2(z2))) + 333396000*s2t*(Dmglst2*(-465 + 348*z2
        - 176*z3) + Mst2*(35 - 20*z2 + 16*z3))*pow3(Mt) + 333396000*Mt*(
        Dmglst2*(20 - 11*z2) + Mst2*(-2 + z2))*pow2(Mst2)*pow3(s2t) + 3456*(-
        6582784 + 2701125*z2 - 2315250*z3)*pow4(Mt) + (4630500*Dmglst2*(38 +
        63*z2 - 90*z3) + Mst2*(378483467 + 9604000*OepS2 - 754771500*S2 -
        14406000*T1ep - 306899250*z2 + 133770000*z3 - 7203000*z4 - 10804500*
        pow2(z2)))*pow3(Mst2)*pow4(s2t))*pow6(Mst1)) + 55566000*Mst2*pow3(log(
        pow2(Mst1)/pow2(Mst2)))*pow4(Mst1)*(4*pow2(Mst1)*pow3(Mst2)*(1166*pow2(
        Mst2)*pow2(Mt)*pow2(s2t) - 2304*Mst2*s2t*pow3(Mt) + 780*Mt*pow3(Mst2)*
        pow3(s2t) - 440*pow4(Mt) + 115*pow4(Mst2)*pow4(s2t)) + Mst2*pow4(Mst1)*
        (2696*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 7808*Mst2*s2t*pow3(Mt) + 3200*Mt*
        pow3(Mst2)*pow3(s2t) - 6400*pow4(Mt) + 409*pow4(Mst2)*pow4(s2t)) - 4*
        Dmglst2*(pow4(Mst1)*(-2896*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 2016*Mst2*
        s2t*pow3(Mt) - 672*Mt*pow3(Mst2)*pow3(s2t) + 7424*pow4(Mt) - 33*pow4(
        Mst2)*pow4(s2t)) - 2*pow2(Mst1)*pow2(Mst2)*(1156*pow2(Mst2)*pow2(Mt)*
        pow2(s2t) + 384*Mst2*s2t*pow3(Mt) + 390*Mt*pow3(Mst2)*pow3(s2t) - 2632*
        pow4(Mt) + 29*pow4(Mst2)*pow4(s2t)) - 3*pow4(Mst2)*(256*pow2(Mst2)*
        pow2(Mt)*pow2(s2t) - 784*Mst2*s2t*pow3(Mt) + 516*Mt*pow3(Mst2)*pow3(
        s2t) - 512*pow4(Mt) + 107*pow4(Mst2)*pow4(s2t))) + 3*(56*pow2(Mst2)*
        pow2(Mt)*pow2(s2t) - 3136*Mst2*s2t*pow3(Mt) + 1040*Mt*pow3(Mst2)*pow3(
        s2t) + 32*pow4(Mt) + 271*pow4(Mst2)*pow4(s2t))*pow5(Mst2) + 6*pow2(
        Dmglst2)*(-512*pow3(Mst2)*pow4(Mt) + 2*Mst2*pow2(Mst1)*(380*pow2(Mst2)*
        pow2(Mt)*pow2(s2t) - 512*Mst2*s2t*pow3(Mt) + 1608*pow4(Mt) + 13*pow4(
        Mst2)*pow4(s2t)) + 256*pow2(Mt)*pow2(s2t)*pow5(Mst2) + 768*Mt*pow3(s2t)
        *pow6(Mst2) + 107*pow4(s2t)*pow7(Mst2))) + Mst2*(4*Dmglst2*pow2(Mst1)*(
        -1764*pow2(Mst2)*pow4(Mst1)*(-112*pow2(Mst2)*pow2(Mt)*pow2(s2t)*(
        1216808 + 150000*OepS2 - 376650*S2 - 225000*T1ep - 3408750*z2 - 833625*
        z3 - 112500*z4 - 168750*pow2(z2)) + 1050*Mst2*s2t*(-322421 + 73760*
        OepS2 - 564876*S2 - 110640*T1ep - 3060894*z2 + 1606920*z3 - 55320*z4 -
        82980*pow2(z2))*pow3(Mt) - 175*Mt*(27211 + 36720*B4 + 1080*DN + 64160*
        OepS2 - 1515564*S2 - 96240*T1ep - 1240446*z2 + 245400*z3 - 12480*z4 -
        72180*pow2(z2))*pow3(Mst2)*pow3(s2t) - 8*(98884301 + 4508000*OepS2 -
        17853750*S2 - 6762000*T1ep - 477659550*z2 + 302883000*z3 - 3381000*z4 -
        5071500*pow2(z2))*pow4(Mt) - 7*(-95041 + 81000*B4 - 108000*D3 + 67500*
        DN - 432000*OepS2 - 7425000*S2 + 648000*T1ep + 5325300*z2 + 4227750*z3
        - 121500*z4)*pow4(Mst2)*pow4(s2t)) - 66150*pow2(Mst1)*pow4(Mst2)*(672*(
        -11674 + 120*B4 - 120*D3 + 60*DN + 17091*S2 + 4849*z2 + 2195*z3 - 540*
        z4)*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 140*Mst2*s2t*(1120*OepS2 + 3*(-
        10061 + 14580*S2 - 560*T1ep - 25134*z2 + 72*z3 + 2024*z4 - 420*pow2(z2)
        ))*pow3(Mt) - 7*Mt*(69997 + 188640*B4 - 3600*DN + 5600*OepS2 + 146772*
        S2 - 8400*T1ep - 1043682*z2 + 1816440*z3 + 32520*z4 - 317340*pow2(z2))*
        pow3(Mst2)*pow3(s2t) - 16*(1547066 + 7840*OepS2 + 258066*S2 - 11760*
        T1ep - 1256046*z2 - 946785*z3 - 5880*z4 - 8820*pow2(z2))*pow4(Mt) + 35*
        (-15707 + 432*B4 - 576*D3 + 360*DN + 224*OepS2 + 543348*S2 - 336*T1ep -
        11802*z2 - 105690*z3 - 2544*z4 - 2844*pow2(z2))*pow4(Mst2)*pow4(s2t)) +
        (72*pow2(Mst2)*pow2(Mt)*pow2(s2t)*(-3650009897 + 2110136000*OepS2 -
        52385772600*S2 - 3165204000*T1ep - 40405228500*z2 + 755286000*z3 -
        1582602000*z4 - 2373903000*pow2(z2)) - 274400*Mst2*s2t*(-20964397 +
        2058400*OepS2 - 47545272*S2 - 3087600*T1ep - 59449002*z2 + 39214920*z3
        - 1543800*z4 - 2315700*pow2(z2))*pow3(Mt) + 17150*Mt*(-31897243 +
        2491360*OepS2 - 90290268*S2 - 3737040*T1ep - 51499698*z2 + 33912840*z3
        - 1868520*z4 - 2802780*pow2(z2))*pow3(Mst2)*pow3(s2t) + 80*(-
        54946675289 + 3063401600*OepS2 - 52952863320*S2 - 4595102400*T1ep -
        225236187120*z2 + 243414477600*z3 - 2297551200*z4 - 3446326800*pow2(z2)
        )*pow4(Mt) + (119405394763 - 11522056000*OepS2 + 478191735000*S2 +
        17283084000*T1ep + 233890773900*z2 - 111792103500*z3 + 8641542000*z4 +
        12962313000*pow2(z2))*pow4(Mst2)*pow4(s2t))*pow6(Mst1) - 333396000*(2*
        Mt + Mst2*s2t)*(14*Mt*(3 + z2) + Mst2*s2t*(5 + 7*z2))*pow2(-2*Mt +
        Mst2*s2t)*pow6(Mst2)) + 588*Mst2*pow2(Dmglst2)*(1134000*pow2(Mst1)*
        pow4(Mst2)*(-8*(-179 + 26*z2)*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 768*Mst2*
        s2t*pow3(Mt) - 192*Mt*pow3(Mst2)*pow3(s2t) + 16*(-155 + 26*z2)*pow4(Mt)
        + (-235 + 26*z2)*pow4(Mst2)*pow4(s2t)) + 15*pow2(Mst2)*pow4(Mst1)*(-
        360*pow2(Mst2)*pow2(Mt)*pow2(s2t)*(-2164753 + 3360*B4 - 3360*D3 + 1680*
        DN - 7840*OepS2 + 169452*S2 + 11760*T1ep + 467818*z2 + 1945895*z3 -
        9240*z4 + 8820*pow2(z2)) + 320*Mst2*s2t*(-1243547 + 15680*OepS2 +
        953208*S2 - 23520*T1ep - 2985438*z2 + 3393915*z3 - 11760*z4 - 17640*
        pow2(z2))*pow3(Mt) + 560*Mt*(-282298 + 61560*B4 - 1620*DN - 2240*OepS2
        - 111213*S2 + 3360*T1ep - 281922*z2 + 983190*z3 + 6540*z4 - 114120*
        pow2(z2))*pow3(Mst2)*pow3(s2t) - 16*(-46632377 + 744800*OepS2 +
        29679210*S2 - 1117200*T1ep - 107181930*z2 + 108350025*z3 - 558600*z4 -
        837900*pow2(z2))*pow4(Mt) - 35*(-924689 + 41040*B4 - 43200*D3 + 22680*
        DN - 1120*OepS2 + 21325356*S2 + 1680*T1ep - 127686*z2 - 3687510*z3 -
        190320*z4 - 37620*pow2(z2))*pow4(Mst2)*pow4(s2t)) - 2*(-36*pow2(Mst2)*
        pow2(Mt)*pow2(s2t)*(68526791 + 2772000*OepS2 - 107403300*S2 - 4158000*
        T1ep - 140027250*z2 + 39514125*z3 - 2079000*z4 - 3118500*pow2(z2)) -
        800*Mst2*s2t*(7384247 + 1183840*OepS2 + 5821092*S2 - 1775760*T1ep -
        59988552*z2 + 25757760*z3 - 887880*z4 - 1331820*pow2(z2))*pow3(Mt) +
        350*Mt*(-2727217 + 437920*OepS2 - 3648132*S2 - 656880*T1ep - 12301638*
        z2 + 4166400*z3 - 328440*z4 - 492660*pow2(z2))*pow3(Mst2)*pow3(s2t) +
        8*(648916519 + 158732000*OepS2 + 66753450*S2 - 238098000*T1ep -
        10760276550*z2 + 6504855000*z3 - 119049000*z4 - 178573500*pow2(z2))*
        pow4(Mt) - 7*(-21452821 + 243000*B4 - 324000*D3 + 202500*DN + 2272000*
        OepS2 - 75724200*S2 - 3408000*T1ep - 55081500*z2 + 67149750*z3 -
        3040500*z4 - 4014000*pow2(z2))*pow4(Mst2)*pow4(s2t))*pow6(Mst1) +
        18144000*pow2(-4*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow6(Mst2)) + 3*Mst2*
        (308700*pow4(Mst1)*pow4(Mst2)*(2*pow2(Mst2)*pow2(Mt)*pow2(s2t)*(1013263
        + 2880*D3 - 2880*DN - 25440*OepS2 - 137052*S2 + 38160*T1ep + 35562*z2 +
        83160*z3 + 36360*z4 + 28620*pow2(z2)) + 24*Mst2*s2t*(-60759 + 1120*
        OepS2 + 25596*S2 - 1680*T1ep - 137386*z2 + 115320*z3 - 12360*z4 - 1260*
        pow2(z2))*pow3(Mt) + 30*Mt*(25611 + 5280*B4 - 48*DN - 224*OepS2 - 324*
        S2 + 336*T1ep - 13614*z2 + 47720*z3 + 2040*z4 - 6660*pow2(z2))*pow3(
        Mst2)*pow3(s2t) + 8*(389597 + 3360*OepS2 - 105948*S2 - 5040*T1ep +
        293898*z2 - 740160*z3 - 2520*z4 - 3780*pow2(z2))*pow4(Mt) - 5*(25289 +
        1440*B4 - 144*D3 + 72*DN - 2208*OepS2 + 298404*S2 + 3312*T1ep + 36318*
        z2 - 94968*z3 + 1440*z4 - 108*pow2(z2))*pow4(Mst2)*pow4(s2t)) + 20580*
        pow2(Mst2)*(4*pow2(Mst2)*pow2(Mt)*pow2(s2t)*(2367149 - 21600*D3 +
        21600*DN - 503200*OepS2 + 14889420*S2 + 754800*T1ep + 9292470*z2 -
        66000*z3 + 247800*z4 + 566100*pow2(z2)) + 40*Mst2*s2t*(-1156193 +
        60320*OepS2 - 1414908*S2 - 90480*T1ep - 2823942*z2 + 2907240*z3 -
        45240*z4 - 67860*pow2(z2))*pow3(Mt) + 100*Mt*(51041 + 7344*B4 + 216*DN
        - 2336*OepS2 + 112428*S2 + 3504*T1ep + 28038*z2 - 78216*z3 + 8880*z4 +
        2628*pow2(z2))*pow3(Mst2)*pow3(s2t) + 8*(7319011 + 167200*OepS2 -
        6744060*S2 - 250800*T1ep + 9301170*z2 - 20715000*z3 - 125400*z4 -
        188100*pow2(z2))*pow4(Mt) + (-7680127 + 10800*B4 - 21600*D3 + 16200*DN
        + 514400*OepS2 - 46307700*S2 - 771600*T1ep - 10868970*z2 + 7186200*z3 -
        466800*z4 - 773100*pow2(z2))*pow4(Mst2)*pow4(s2t))*pow6(Mst1) +
        55566000*pow2(Mst1)*(-8*(87 + 38*z2)*pow2(Mst2)*pow2(Mt)*pow2(s2t) -
        1024*Mst2*s2t*pow3(Mt) + 256*Mt*pow3(Mst2)*pow3(s2t) + 16*(71 + 38*z2)*
        pow4(Mt) + (135 + 38*z2)*pow4(Mst2)*pow4(s2t))*pow6(Mst2) + (-16*pow2(
        Mst2)*pow2(Mt)*pow2(s2t)*(-25002374272 + 4703902000*OepS2 -
        184815757350*S2 - 7055853000*T1ep - 84661721025*z2 + 16805113500*z3 -
        3527926500*z4 - 5291889750*pow2(z2)) + 1097600*Mst2*s2t*(-2016907 +
        109600*OepS2 - 3520152*S2 - 164400*T1ep - 3448554*z2 + 3465480*z3 -
        82200*z4 - 123300*pow2(z2))*pow3(Mt) - 68600*Mt*(-2759935 + 70240*OepS2
        - 3305340*S2 - 105360*T1ep - 629034*z2 + 2677800*z3 - 52680*z4 - 79020*
        pow2(z2))*pow3(Mst2)*pow3(s2t) + 16*(187545955063 + 3204992000*OepS2 -
        128216692800*S2 - 4807488000*T1ep + 129009640200*z2 - 399103824000*z3 -
        2403744000*z4 - 3605616000*pow2(z2))*pow4(Mt) + (93508520089 +
        2000376000*B4 + 222264000*D3 - 222264000*DN + 4925480000*OepS2 -
        312095700000*S2 - 7388220000*T1ep - 77106571500*z2 + 21506100000*z3 -
        2360526000*z4 - 5541165000*pow2(z2))*pow4(Mst2)*pow4(s2t))*pow8(Mst1) +
        7112448000*pow2(Mt)*(-2*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow8(Mst2))) +
        1260*log(pow2(Mst1)/pow2(Mst2))*(30*pow2(Dmsqst2)*pow4(Mst1)*(5*pow2(
        s2t)*(-26460*Mst2*Mt*s2t*(-121 + 60*z2) + 30*(-20869 + 35280*S2 + 5292*
        z2)*pow2(Mt) + (-8753 - 246960*S2 + 79380*z2)*pow2(Mst2)*pow2(s2t))*
        pow4(Mst1) + 294*pow2(Mst1)*(8*(-887 + 2250*S2)*pow2(Mst2)*pow2(Mt)*
        pow2(s2t) - 57600*Mst2*s2t*(-3 + 2*z2)*pow3(Mt) + 7200*Mt*(-2 + z2)*
        pow3(Mst2)*pow3(s2t) + 48*(-3871 + 3000*z2)*pow4(Mt) + (1837 - 7200*S2
        - 900*z2)*pow4(Mst2)*pow4(s2t)) - 14700*pow2(Mst2)*(-24*(2 + 9*S2)*
        pow2(Mst2)*pow2(Mt)*pow2(s2t) + 144*Mst2*s2t*(-15 + 8*z2)*pow3(Mt) -
        144*Mt*(-2 + z2)*pow3(Mst2)*pow3(s2t) + (1460 - 864*z2)*pow4(Mt) + 3*(-
        1 + 18*S2 + 3*z2)*pow4(Mst2)*pow4(s2t))) + 120*Dmsqst2*pow2(Mst1)*(441*
        Mst2*pow4(Mst1)*(Mst2*(-8*(-517 + 450*S2)*pow2(Mst2)*pow2(Mt)*pow2(s2t)
        - 19200*Mst2*s2t*(-2 + z2)*pow3(Mt) + 2400*Mt*(-1 + z2)*pow3(Mst2)*
        pow3(s2t) + 16*(-3067 + 1200*z2)*pow4(Mt) + (-317 + 3150*S2 - 300*z2)*
        pow4(Mst2)*pow4(s2t)) - 150*Dmglst2*(8*pow2(Mst2)*pow2(Mt)*pow2(s2t) -
        128*Mst2*s2t*(-12 + 7*z2)*pow3(Mt) + 16*Mt*(-5 + 3*z2)*pow3(Mst2)*pow3(
        s2t) + 16*(-137 + 72*z2)*pow4(Mt) + (13 - 4*z2)*pow4(Mst2)*pow4(s2t)))
        + 7350*pow2(Mst1)*pow3(Mst2)*(3*Mst2*(-4*(16 + 27*S2)*pow2(Mst2)*pow2(
        Mt)*pow2(s2t) - 96*Mst2*s2t*(-5 + 4*z2)*pow3(Mt) + 48*Mt*(-1 + z2)*
        pow3(Mst2)*pow3(s2t) + 16*(-23 + 12*z2)*pow4(Mt) + (17 + 27*S2)*pow4(
        Mst2)*pow4(s2t)) + Dmglst2*(96*Mst2*s2t*(-59 + 36*z2)*pow3(Mt) - 144*
        Mt*(-5 + 3*z2)*pow3(Mst2)*pow3(s2t) - 8*(-571 + 360*z2)*pow4(Mt) + 18*(
        -2 + 3*z2)*pow4(Mst2)*pow4(s2t))) + (2*Mst2*(8*Mst2*(58853 - 44100*S2)
        - 11025*Dmglst2*(-37 + 60*z2))*pow2(Mt)*pow2(s2t) + 1058400*s2t*(Mst2*(
        25 - 8*z2) + 11*Dmglst2*(-17 + 8*z2))*pow3(Mt) + 264600*(3*Dmglst2 -
        Mst2)*Mt*pow2(Mst2)*pow3(s2t) + 288*(-175121 + 44100*z2)*pow4(Mt) + (-
        11025*Dmglst2*(-83 + 60*z2) + Mst2*(-621671 + 308700*S2 + 132300*z2))*
        pow3(Mst2)*pow4(s2t))*pow6(Mst1) - 66150*pow2(-4*pow2(Mt) + pow2(Mst2)*
        pow2(s2t))*pow6(Mst2)) + Mst2*(588*Mst2*pow2(Dmglst2)*(-450*pow2(Mst1)*
        pow4(Mst2)*(-408*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 2048*Mst2*s2t*pow3(Mt)
        + 512*Mt*pow3(Mst2)*pow3(s2t) + 1328*pow4(Mt) + 51*pow4(Mst2)*pow4(s2t)
        ) + 2*pow2(Mst2)*pow4(Mst1)*(300*(6485 + 1134*S2 - 6306*z2)*pow2(Mst2)*
        pow2(Mt)*pow2(s2t) + 1600*Mst2*s2t*(2182 + 378*S2 - 1665*z2)*pow3(Mt) -
        600*Mt*(439 + 252*S2 + 135*z2)*pow3(Mst2)*pow3(s2t) - 16*(-56837 +
        89775*S2 - 162900*z2)*pow4(Mt) + 75*(-1388 + 63*S2 - 618*z2)*pow4(Mst2)
        *pow4(s2t)) + (24*(102713 + 133650*S2 - 173700*z2)*pow2(Mst2)*pow2(Mt)*
        pow2(s2t) + 400*Mst2*s2t*(11587 + 76104*S2 - 49608*z2)*pow3(Mt) - 200*
        Mt*(-4958 + 24633*S2 - 90*z2)*pow3(Mst2)*pow3(s2t) + (3934288 -
        40816800*S2 + 34804800*z2)*pow4(Mt) + (332311 + 511200*S2 - 49050*z2)*
        pow4(Mst2)*pow4(s2t))*pow6(Mst1) + 7200*pow2(-4*pow2(Mt) + pow2(Mst2)*
        pow2(s2t))*pow6(Mst2)) + 3*Mst2*(9800*pow4(Mst1)*pow4(Mst2)*(-6*(-4534
        + 4293*S2 + 1830*z2 - 72*z3)*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 72*Mst2*
        s2t*(26 + 189*S2 - 758*z2)*pow3(Mt) - 126*Mt*(-208 + 27*S2 - 2*z2)*
        pow3(Mst2)*pow3(s2t) + 8*(-553 + 1701*S2 + 4023*z2 + 108*z3)*pow4(Mt) +
        3*(1148 + 1863*S2 + 276*z2 - 18*z3)*pow4(Mst2)*pow4(s2t)) + 1960*pow2(
        Mst2)*(-4*(-52322 + 84915*S2 + 3060*z2 + 540*z3)*pow2(Mst2)*pow2(Mt)*
        pow2(s2t) + 1080*Mst2*s2t*(42 + 377*S2 - 330*z2)*pow3(Mt) - 180*Mt*(221
        + 219*S2 - 234*z2)*pow3(Mst2)*pow3(s2t) + 8*(-22352 + 28215*S2 + 36270*
        z2)*pow4(Mt) + (-16051 + 86805*S2 - 8325*z2)*pow4(Mst2)*pow4(s2t))*
        pow6(Mst1) + 88200*pow2(Mst1)*(-616*pow2(Mst2)*pow2(Mt)*pow2(s2t) -
        1024*Mst2*s2t*pow3(Mt) + 256*Mt*pow3(Mst2)*pow3(s2t) + 976*pow4(Mt) +
        125*pow4(Mst2)*pow4(s2t))*pow6(Mst2) + (-24*(-19407863 + 50398950*S2 +
        2866500*z2)*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 940800*Mst2*s2t*(167 +
        2055*S2 - 912*z2)*pow3(Mt) - 58800*Mt*(46 + 1317*S2 - 918*z2)*pow3(
        Mst2)*pow3(s2t) + 32*(-35585111 + 25754400*S2 + 32281200*z2)*pow4(Mt) +
        3*(-23468297 + 26386500*S2 + 661500*z2 + 176400*z3)*pow4(Mst2)*pow4(
        s2t))*pow8(Mst1) + 11289600*pow2(Mt)*(-2*pow2(Mt) + pow2(Mst2)*pow2(
        s2t))*pow8(Mst2)) + 4*Dmglst2*(-14700*pow4(Mst1)*pow4(Mst2)*(12*(-3427
        + 1782*z2)*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 4*Mst2*s2t*(-8398 + 2835*S2
        - 3258*z2)*pow3(Mt) - 9*Mt*(1640 + 315*S2 - 1146*z2)*pow3(Mst2)*pow3(
        s2t) + (26276 - 9072*S2 + 28224*z2)*pow4(Mt) + 9*(277 + 63*S2 - 140*z2)
        *pow4(Mst2)*pow4(s2t)) - 294*pow2(Mst2)*(-8*(77657 + 202500*S2 - 55800*
        z2)*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 200*Mst2*s2t*(-11408 + 37341*S2 -
        17730*z2)*pow3(Mt) - 300*Mt*(-503 + 3609*S2 - 1494*z2)*pow3(Mst2)*pow3(
        s2t) - 48*(8581 + 72450*S2 - 137700*z2)*pow4(Mt) + (-81643 + 291600*S2
        + 5850*z2)*pow4(Mst2)*pow4(s2t))*pow6(Mst1) + 132300*pow2(Mst1)*(-536*
        pow2(Mst2)*pow2(Mt)*pow2(s2t) - 512*Mst2*s2t*pow3(Mt) + 128*Mt*pow3(
        Mst2)*pow3(s2t) + 560*pow4(Mt) + 131*pow4(Mst2)*pow4(s2t))*pow6(Mst2) +
        (72*(-4261807 + 33912900*S2 - 73500*z2)*pow2(Mst2)*pow2(Mt)*pow2(s2t) -
        78400*Mst2*s2t*(-4724 + 115785*S2 - 29484*z2)*pow3(Mt) + 4900*Mt*(18652
        + 140139*S2 - 38826*z2)*pow3(Mst2)*pow3(s2t) + 80*(42300121 + 49233240*
        S2 - 64139040*z2)*pow4(Mt) + (8287903 - 185175900*S2 + 9525600*z2)*
        pow4(Mst2)*pow4(s2t))*pow8(Mst1) + 16934400*pow2(Mt)*(-2*pow2(Mt) +
        pow2(Mst2)*pow2(s2t))*pow8(Mst2)))) + 1587600*pow2(log(pow2(Mst1)/pow2(
        Mst2)))*(90*pow2(Dmsqst2)*(-280*pow2(Mst2)*pow4(Mst1)*pow4(Mt) + 14*(
        14*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 32*pow4(Mt) - 3*pow4(Mst2)*pow4(s2t)
        )*pow6(Mst1) + 15*(12*pow2(Mt)*pow2(s2t) - 5*pow2(Mst2)*pow4(s2t))*
        pow8(Mst1)) - 90*Dmsqst2*(-7*Mst2*(-48*pow2(Mt)*pow2(s2t)*pow3(Mst2) +
        320*Dmglst2*pow4(Mt) - 144*Mst2*pow4(Mt) + 31*pow4(s2t)*pow5(Mst2))*
        pow6(Mst1) + 2*(74*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 1120*(3*Dmglst2 -
        Mst2)*s2t*pow3(Mt) + 2736*pow4(Mt) - 53*pow4(Mst2)*pow4(s2t))*pow8(
        Mst1) + 140*pow4(Mst1)*(-4*pow4(Mst2)*pow4(Mt) + pow4(s2t)*pow8(Mst2)))
        + Mst2*(70*pow4(Mst1)*(3192*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 14880*Mst2*
        s2t*pow3(Mt) + 5664*Mt*pow3(Mst2)*pow3(s2t) + 7744*pow4(Mt) + 1113*
        pow4(Mst2)*pow4(s2t))*pow5(Mst2) + 560*pow3(Mst2)*(479*pow2(Mst2)*pow2(
        Mt)*pow2(s2t) - 2448*Mst2*s2t*pow3(Mt) + 366*Mt*pow3(Mst2)*pow3(s2t) +
        1136*pow4(Mt) + 23*pow4(Mst2)*pow4(s2t))*pow6(Mst1) - 7*Mst2*pow2(
        Dmglst2)*(-15*pow2(Mst1)*pow4(Mst2)*(-1240*pow2(Mst2)*pow2(Mt)*pow2(
        s2t) - 3072*Mst2*s2t*pow3(Mt) + 768*Mt*pow3(Mst2)*pow3(s2t) + 944*pow4(
        Mt) + 443*pow4(Mst2)*pow4(s2t)) + 2*pow2(Mst2)*pow4(Mst1)*(7740*pow2(
        Mst2)*pow2(Mt)*pow2(s2t) - 61760*Mst2*s2t*pow3(Mt) + 18480*Mt*pow3(
        Mst2)*pow3(s2t) - 66208*pow4(Mt) + 4455*pow4(Mst2)*pow4(s2t)) + (
        133776*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 398240*Mst2*s2t*pow3(Mt) +
        24960*Mt*pow3(Mst2)*pow3(s2t) - 1267120*pow4(Mt) + 4293*pow4(Mst2)*
        pow4(s2t))*pow6(Mst1) + 1440*pow2(-4*pow2(Mt) + pow2(Mst2)*pow2(s2t))*
        pow6(Mst2)) + 105*pow2(Mst1)*(-600*pow2(Mst2)*pow2(Mt)*pow2(s2t) -
        1024*Mst2*s2t*pow3(Mt) + 256*Mt*pow3(Mst2)*pow3(s2t) + 944*pow4(Mt) +
        123*pow4(Mst2)*pow4(s2t))*pow7(Mst2) + 4*Mst2*(73657*pow2(Mst2)*pow2(
        Mt)*pow2(s2t) - 316400*Mst2*s2t*pow3(Mt) + 35000*Mt*pow3(Mst2)*pow3(
        s2t) + 68678*pow4(Mt) + 11764*pow4(Mst2)*pow4(s2t))*pow8(Mst1) - 2*
        Dmglst2*(-140*pow4(Mst1)*pow4(Mst2)*(1806*pow2(Mst2)*pow2(Mt)*pow2(s2t)
        + 1592*Mst2*s2t*pow3(Mt) + 1104*Mt*pow3(Mst2)*pow3(s2t) - 3808*pow4(Mt)
        + 213*pow4(Mst2)*pow4(s2t)) + 7*pow2(Mst2)*(-12992*pow2(Mst2)*pow2(Mt)*
        pow2(s2t) - 175520*Mst2*s2t*pow3(Mt) + 3600*Mt*pow3(Mst2)*pow3(s2t) +
        285744*pow4(Mt) + 4969*pow4(Mst2)*pow4(s2t))*pow6(Mst1) - 105*pow2(
        Mst1)*(-728*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 1536*Mst2*s2t*pow3(Mt) +
        384*Mt*pow3(Mst2)*pow3(s2t) + 944*pow4(Mt) + 187*pow4(Mst2)*pow4(s2t))*
        pow6(Mst2) + (-420096*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 1672160*Mst2*s2t*
        pow3(Mt) + 68880*Mt*pow3(Mst2)*pow3(s2t) + 3253120*pow4(Mt) + 1377*
        pow4(Mst2)*pow4(s2t))*pow8(Mst1) - 26880*pow2(Mt)*(-2*pow2(Mt) + pow2(
        Mst2)*pow2(s2t))*pow8(Mst2)) + 13440*pow2(Mt)*(-2*pow2(Mt) + pow2(Mst2)
        *pow2(s2t))*pow9(Mst2))))/(5.00094e8*pow4(Mst1)*pow4(Mt)*pow6(Mst2));

   return result;
}

/**
 * @return returns the susy log^1 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double H6b::calc_coef_at_as2_no_sm_logs_log1() const {

   const double result =
      (6000*pow2(Dmsqst2)*(72*(2*Mst2*s2t*(3 - 2*z2) + Mt*(-7 + 5*z2))*pow2(
        Mst1) + (36*Mst2*s2t*(7 - 4*z2) + Mt*(-161 + 108*z2))*pow2(Mst2))*pow3(
        Mt)*pow4(Mst1) + 3000*Dmsqst2*pow2(Mst1)*(9*Mst2*pow4(Mst1)*(-32*
        Dmglst2*(2*Mst2*s2t*(12 - 7*z2) + 3*Mt*(-11 + 6*z2))*pow3(Mt) + Mst2*(
        8*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 64*Mst2*s2t*(-2 + z2)*pow3(Mt) + 16*(
        -9 + 4*z2)*pow4(Mt) + pow4(Mst2)*pow4(s2t))) + pow2(Mst1)*pow3(Mst2)*(-
        8*Dmglst2*(12*Mst2*s2t*(29 - 18*z2) + Mt*(-283 + 180*z2))*pow3(Mt) + 3*
        Mst2*(-48*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 96*Mst2*s2t*(-3 + 2*z2)*pow3(
        Mt) + 32*(-8 + 3*z2)*pow4(Mt) + 3*pow4(Mst2)*pow4(s2t))) + 9*(16*s2t*(-
        83*Dmglst2 + 9*Mst2 + 44*Dmglst2*z2 - 4*Mst2*z2)*pow3(Mt) + 8*(-31 +
        12*z2)*pow4(Mt) - pow4(Mst2)*pow4(s2t))*pow6(Mst1) - 9*pow2(-4*pow2(Mt)
        + pow2(Mst2)*pow2(s2t))*pow6(Mst2)) + 450*Mst2*pow2(log(pow2(Mst1)/
        pow2(Mst2)))*pow4(Mst1)*(8*pow2(Mst1)*pow3(Mst2)*(425*pow2(Mst2)*pow2(
        Mt)*pow2(s2t) - 704*Mst2*s2t*pow3(Mt) + 212*Mt*pow3(Mst2)*pow3(s2t) -
        76*pow4(Mt) + 3*pow4(Mst2)*pow4(s2t)) + Mst2*pow4(Mst1)*(2688*pow2(
        Mst2)*pow2(Mt)*pow2(s2t) - 5056*Mst2*s2t*pow3(Mt) + 1024*Mt*pow3(Mst2)*
        pow3(s2t) - 2848*pow4(Mt) + 41*pow4(Mst2)*pow4(s2t)) - 8*Dmglst2*(8*Mt*
        (-157*Mst2*s2t*pow2(Mt) - 108*Mt*pow2(Mst2)*pow2(s2t) + 398*pow3(Mt) -
        16*pow3(Mst2)*pow3(s2t))*pow4(Mst1) + pow4(Mst2)*(-320*pow2(Mst2)*pow2(
        Mt)*pow2(s2t) + 848*Mst2*s2t*pow3(Mt) - 532*Mt*pow3(Mst2)*pow3(s2t) +
        512*pow4(Mt) - 99*pow4(Mst2)*pow4(s2t)) + pow2(Mst1)*pow2(Mst2)*(-768*
        pow2(Mst2)*pow2(Mt)*pow2(s2t) - 448*Mst2*s2t*pow3(Mt) - 212*Mt*pow3(
        Mst2)*pow3(s2t) + 1896*pow4(Mt) + 35*pow4(Mst2)*pow4(s2t))) + (296*
        pow2(Mst2)*pow2(Mt)*pow2(s2t) - 6784*Mst2*s2t*pow3(Mt) + 2208*Mt*pow3(
        Mst2)*pow3(s2t) - 672*pow4(Mt) + 519*pow4(Mst2)*pow4(s2t))*pow5(Mst2) +
        4*Mst2*pow2(Dmglst2)*(320*pow2(Mt)*pow2(s2t)*pow4(Mst2) - 512*pow2(
        Mst2)*pow4(Mt) + pow2(Mst1)*(768*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 1152*
        Mst2*s2t*pow3(Mt) + 3640*pow4(Mt) - 35*pow4(Mst2)*pow4(s2t)) + 768*Mt*
        pow3(s2t)*pow5(Mst2) + 99*pow4(s2t)*pow6(Mst2))) + Mst2*(4*Mst2*pow2(
        Dmglst2)*(-225*pow2(Mst1)*pow4(Mst2)*(-664*pow2(Mst2)*pow2(Mt)*pow2(
        s2t) - 2048*Mst2*s2t*pow3(Mt) + 512*Mt*pow3(Mst2)*pow3(s2t) + 1840*
        pow4(Mt) + 83*pow4(Mst2)*pow4(s2t)) + 2*pow2(Mst2)*pow4(Mst1)*(-150*(-
        7121 + 6336*z2)*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 400*Mst2*s2t*(-1934 +
        1791*z2)*pow3(Mt) - 300*Mt*(-875 + 648*z2)*pow3(Mst2)*pow3(s2t) + 8*(
        71687 + 76050*z2)*pow4(Mt) + 75*(-205 + 129*z2)*pow4(Mst2)*pow4(s2t)) +
        25*(12*(7903 - 5760*z2)*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 96*Mst2*s2t*(-
        1055 + 2181*z2)*pow3(Mt) + 13992*Mt*pow3(Mst2)*pow3(s2t) + 80*(-605 +
        4608*z2)*pow4(Mt) + 3*(409 - 258*z2)*pow4(Mst2)*pow4(s2t))*pow6(Mst1) +
        3600*pow2(-4*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow6(Mst2)) - 100*
        Dmglst2*(2*pow4(Mst1)*pow4(Mst2)*(12*(-3565 + 1728*z2)*pow2(Mst2)*pow2(
        Mt)*pow2(s2t) + 8*Mst2*s2t*(-3659 + 252*z2)*pow3(Mt) + 72*Mt*(-192 +
        91*z2)*pow3(Mst2)*pow3(s2t) + 44*(427 + 432*z2)*pow4(Mt) + 9*(211 - 86*
        z2)*pow4(Mst2)*pow4(s2t)) + 6*pow2(Mst2)*(4*(-1955 + 1152*z2)*pow2(
        Mst2)*pow2(Mt)*pow2(s2t) - 48*Mst2*s2t*(-13 + 262*z2)*pow3(Mt) + 24*Mt*
        (111 - 17*z2)*pow3(Mst2)*pow3(s2t) + 8*(-1471 + 3240*z2)*pow4(Mt) + 3*(
        -1 + 86*z2)*pow4(Mst2)*pow4(s2t))*pow6(Mst1) - 18*pow2(Mst1)*(-536*
        pow2(Mst2)*pow2(Mt)*pow2(s2t) - 512*Mst2*s2t*pow3(Mt) + 128*Mt*pow3(
        Mst2)*pow3(s2t) + 560*pow4(Mt) + 131*pow4(Mst2)*pow4(s2t))*pow6(Mst2) +
        3*(192*(-77 + 48*z2)*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 48*Mst2*s2t*(-867
        + 1012*z2)*pow3(Mt) - 276*Mt*pow3(Mst2)*pow3(s2t) + 16*(-6445 + 6768*
        z2)*pow4(Mt) - 441*pow4(Mst2)*pow4(s2t))*pow8(Mst1) - 2304*pow2(Mt)*(-
        2*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow8(Mst2)) + 25*Mst2*(4*pow4(Mst1)*
        pow4(Mst2)*(-72*(-445 + 96*z2)*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 144*
        Mst2*s2t*(-121 + 296*z2)*pow3(Mt) - 144*Mt*(-151 + 19*z2)*pow3(Mst2)*
        pow3(s2t) + 8*(-1993 + 2151*z2 + 216*z3)*pow4(Mt) + 81*(8 + 5*z2)*pow4(
        Mst2)*pow4(s2t)) - 36*pow2(Mst2)*(-696*pow2(Mst2)*pow2(Mt)*pow2(s2t) +
        32*Mst2*s2t*(-179 + 128*z2)*pow3(Mt) + 16*Mt*(148 - 17*z2)*pow3(Mst2)*
        pow3(s2t) - 8*(-403 + 300*z2)*pow4(Mt) + (551 + 86*z2)*pow4(Mst2)*pow4(
        s2t))*pow6(Mst1) + 36*pow2(Mst1)*(-744*pow2(Mst2)*pow2(Mt)*pow2(s2t) -
        1024*Mst2*s2t*pow3(Mt) + 256*Mt*pow3(Mst2)*pow3(s2t) + 1232*pow4(Mt) +
        141*pow4(Mst2)*pow4(s2t))*pow6(Mst2) + 9*(864*pow2(Mst2)*pow2(Mt)*pow2(
        s2t) - 512*Mst2*s2t*(-49 + 31*z2)*pow3(Mt) - 784*Mt*pow3(Mst2)*pow3(
        s2t) + 8*(-2503 + 1408*z2)*pow4(Mt) + (1547 + 164*z2)*pow4(Mst2)*pow4(
        s2t))*pow8(Mst1) + 4608*pow2(Mt)*(-2*pow2(Mt) + pow2(Mst2)*pow2(s2t))*
        pow8(Mst2))) + 30*log(pow2(Mst1)/pow2(Mst2))*(7200*pow2(Dmsqst2)*(pow2(
        Mst1) - pow2(Mst2))*pow4(Mst1)*pow4(Mt) - 1800*Dmsqst2*pow4(Mst1)*(8*(
        5*Mt + 6*Dmglst2*s2t - 2*Mst2*s2t)*pow3(Mt)*pow4(Mst1) - 8*pow4(Mst2)*
        pow4(Mt) - Mst2*pow2(Mst1)*(16*Dmglst2*pow4(Mt) - 8*Mst2*pow4(Mt) +
        pow4(s2t)*pow5(Mst2)) + pow4(s2t)*pow8(Mst2)) + Mst2*(2*Mst2*pow2(
        Dmglst2)*(2*pow2(Mst2)*pow4(Mst1)*(3240*pow2(Mst2)*pow2(Mt)*pow2(s2t) +
        45920*Mst2*s2t*pow3(Mt) - 9720*Mt*pow3(Mst2)*pow3(s2t) + 37408*pow4(Mt)
        - 4635*pow4(Mst2)*pow4(s2t)) + 45*pow2(Mst1)*pow4(Mst2)*(-456*pow2(
        Mst2)*pow2(Mt)*pow2(s2t) - 1024*Mst2*s2t*pow3(Mt) + 256*Mt*pow3(Mst2)*
        pow3(s2t) + 400*pow4(Mt) + 153*pow4(Mst2)*pow4(s2t)) + 2*(-34080*pow2(
        Mst2)*pow2(Mt)*pow2(s2t) - 147560*Mst2*s2t*pow3(Mt) - 5880*Mt*pow3(
        Mst2)*pow3(s2t) + 415932*pow4(Mt) + 1005*pow4(Mst2)*pow4(s2t))*pow6(
        Mst1) - 1440*pow2(-4*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow6(Mst2)) - 20*
        Dmglst2*(-4*pow4(Mst1)*pow4(Mst2)*(1860*pow2(Mst2)*pow2(Mt)*pow2(s2t) +
        1472*Mst2*s2t*pow3(Mt) + 774*Mt*pow3(Mst2)*pow3(s2t) - 2464*pow4(Mt) +
        15*pow4(Mst2)*pow4(s2t)) + 2*pow2(Mst2)*(-432*pow2(Mst2)*pow2(Mt)*pow2(
        s2t) - 11656*Mst2*s2t*pow3(Mt) + 996*Mt*pow3(Mst2)*pow3(s2t) + 18748*
        pow4(Mt) + 207*pow4(Mst2)*pow4(s2t))*pow6(Mst1) - 3*pow2(Mst1)*(-856*
        pow2(Mst2)*pow2(Mt)*pow2(s2t) - 1536*Mst2*s2t*pow3(Mt) + 384*Mt*pow3(
        Mst2)*pow3(s2t) + 1200*pow4(Mt) + 203*pow4(Mst2)*pow4(s2t))*pow6(Mst2)
        + 4*(-738*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 10004*Mst2*s2t*pow3(Mt) +
        135*Mt*pow3(Mst2)*pow3(s2t) + 17532*pow4(Mt) + 12*pow4(Mst2)*pow4(s2t))
        *pow8(Mst1) - 768*pow2(Mt)*(-2*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow8(
        Mst2)) + 5*Mst2*(4*pow4(Mst1)*pow4(Mst2)*(3372*pow2(Mst2)*pow2(Mt)*
        pow2(s2t) - 9024*Mst2*s2t*pow3(Mt) + 3912*Mt*pow3(Mst2)*pow3(s2t) +
        6112*pow4(Mt) + 579*pow4(Mst2)*pow4(s2t)) - 24*pow2(Mst2)*(-534*pow2(
        Mst2)*pow2(Mt)*pow2(s2t) + 1824*Mst2*s2t*pow3(Mt) - 60*Mt*pow3(Mst2)*
        pow3(s2t) - 934*pow4(Mt) + 97*pow4(Mst2)*pow4(s2t))*pow6(Mst1) + 6*
        pow2(Mst1)*(-728*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 1024*Mst2*s2t*pow3(Mt)
        + 256*Mt*pow3(Mst2)*pow3(s2t) + 1200*pow4(Mt) + 139*pow4(Mst2)*pow4(
        s2t))*pow6(Mst2) + 3*(4640*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 15296*Mst2*
        s2t*pow3(Mt) + 240*Mt*pow3(Mst2)*pow3(s2t) + 7128*pow4(Mt) - 279*pow4(
        Mst2)*pow4(s2t))*pow8(Mst1) + 768*pow2(Mt)*(-2*pow2(Mt) + pow2(Mst2)*
        pow2(s2t))*pow8(Mst2)))))/(1350.*pow4(Mst1)*pow4(Mt)*pow6(Mst2));

   return result;
}

/**
 * @return returns the susy log^2 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double H6b::calc_coef_at_as2_no_sm_logs_log2() const {

   const double result =
      (Mst2*pow2(Dmglst2)*(pow2(Mst2)*pow4(Mst1)*(33360*pow2(Mst2)*pow2(Mt)*
        pow2(s2t) + 52480*Mst2*s2t*pow3(Mt) + 18896*pow4(Mt) - 10005*pow4(Mst2)
        *pow4(s2t)) + 15*pow2(Mst1)*pow4(Mst2)*(-1496*pow2(Mst2)*pow2(Mt)*pow2(
        s2t) - 3072*Mst2*s2t*pow3(Mt) + 768*Mt*pow3(Mst2)*pow3(s2t) + 1456*
        pow4(Mt) + 475*pow4(Mst2)*pow4(s2t)) + 5*(120*pow2(Mst2)*pow2(Mt)*pow2(
        s2t) - 22656*Mst2*s2t*pow3(Mt) - 2304*Mt*pow3(Mst2)*pow3(s2t) + 52384*
        pow4(Mt) + 879*pow4(Mst2)*pow4(s2t))*pow6(Mst1) - 1440*pow2(-4*pow2(Mt)
        + pow2(Mst2)*pow2(s2t))*pow6(Mst2)) + 30*log(pow2(Mst1)/pow2(Mst2))*
        pow4(Mst1)*(pow2(Mst1)*pow3(Mst2)*(1096*pow2(Mst2)*pow2(Mt)*pow2(s2t) -
        1280*Mst2*s2t*pow3(Mt) + 328*Mt*pow3(Mst2)*pow3(s2t) - 32*pow4(Mt) -
        91*pow4(Mst2)*pow4(s2t)) - Mst2*pow4(Mst1)*(-768*pow2(Mst2)*pow2(Mt)*
        pow2(s2t) + 1216*Mst2*s2t*pow3(Mt) + 400*pow4(Mt) + 41*pow4(Mst2)*pow4(
        s2t)) - 2*Dmglst2*(32*pow2(Mt)*(-57*Mst2*Mt*s2t + 119*pow2(Mt) - 24*
        pow2(Mst2)*pow2(s2t))*pow4(Mst1) + pow4(Mst2)*(-384*pow2(Mst2)*pow2(Mt)
        *pow2(s2t) + 912*Mst2*s2t*pow3(Mt) - 548*Mt*pow3(Mst2)*pow3(s2t) + 512*
        pow4(Mt) - 91*pow4(Mst2)*pow4(s2t)) + pow2(Mst1)*pow2(Mst2)*(-768*pow2(
        Mst2)*pow2(Mt)*pow2(s2t) - 640*Mst2*s2t*pow3(Mt) - 164*Mt*pow3(Mst2)*
        pow3(s2t) + 2016*pow4(Mt) + 91*pow4(Mst2)*pow4(s2t))) + 4*(14*pow2(
        Mst2)*pow2(Mt)*pow2(s2t) - 456*Mst2*s2t*pow3(Mt) + 146*Mt*pow3(Mst2)*
        pow3(s2t) - 166*pow4(Mt) + 33*pow4(Mst2)*pow4(s2t))*pow5(Mst2) + pow2(
        Dmglst2)*(-512*pow3(Mst2)*pow4(Mt) + Mst2*pow2(Mst1)*(768*pow2(Mst2)*
        pow2(Mt)*pow2(s2t) - 1280*Mst2*s2t*pow3(Mt) + 4000*pow4(Mt) - 91*pow4(
        Mst2)*pow4(s2t)) + 384*pow2(Mt)*pow2(s2t)*pow5(Mst2) + 768*Mt*pow3(s2t)
        *pow6(Mst2) + 91*pow4(s2t)*pow7(Mst2))) - 10*Dmglst2*(pow4(Mst1)*pow4(
        Mst2)*(-8208*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 4640*Mst2*s2t*pow3(Mt) -
        1968*Mt*pow3(Mst2)*pow3(s2t) + 4816*pow4(Mt) + 849*pow4(Mst2)*pow4(s2t)
        ) - 3*pow2(Mst2)*(-472*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 2496*Mst2*s2t*
        pow3(Mt) - 1040*Mt*pow3(Mst2)*pow3(s2t) - 3360*pow4(Mt) + 37*pow4(Mst2)
        *pow4(s2t))*pow6(Mst1) - 3*pow2(Mst1)*(-984*pow2(Mst2)*pow2(Mt)*pow2(
        s2t) - 1536*Mst2*s2t*pow3(Mt) + 384*Mt*pow3(Mst2)*pow3(s2t) + 1456*
        pow4(Mt) + 219*pow4(Mst2)*pow4(s2t))*pow6(Mst2) + 3*(-2240*Mst2*s2t*
        pow3(Mt) + 3744*pow4(Mt) - 59*pow4(Mst2)*pow4(s2t))*pow8(Mst1) - 768*
        pow2(Mt)*(-2*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow8(Mst2)) + 5*Mst2*(-
        720*pow2(Dmsqst2)*pow4(Mst1)*pow4(Mt) + 1440*Dmsqst2*pow2(Mst2)*pow4(
        Mst1)*pow4(Mt) + pow4(Mst1)*pow4(Mst2)*(8592*pow2(Mst2)*pow2(Mt)*pow2(
        s2t) - 4800*Mst2*s2t*pow3(Mt) + 3936*Mt*pow3(Mst2)*pow3(s2t) + 4256*
        pow4(Mt) - 69*pow4(Mst2)*pow4(s2t)) - 3*pow2(Mst2)*(600*pow2(Mst2)*
        pow2(Mt)*pow2(s2t) - 128*Mst2*s2t*pow3(Mt) + 1568*Mt*pow3(Mst2)*pow3(
        s2t) + 672*pow4(Mt) + 355*pow4(Mst2)*pow4(s2t))*pow6(Mst1) + 3*pow2(
        Mst1)*(-856*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 1024*Mst2*s2t*pow3(Mt) +
        256*Mt*pow3(Mst2)*pow3(s2t) + 1456*pow4(Mt) + 155*pow4(Mst2)*pow4(s2t))
        *pow6(Mst2) + (384*Mst2*s2t*pow3(Mt) - 2400*pow4(Mt) + 717*pow4(Mst2)*
        pow4(s2t))*pow8(Mst1) + 384*pow2(Mt)*(-2*pow2(Mt) + pow2(Mst2)*pow2(
        s2t))*pow8(Mst2)))/(45.*pow4(Mst1)*pow4(Mt)*pow5(Mst2));

   return result;
}

/**
 * @return returns the susy log^3 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double H6b::calc_coef_at_as2_no_sm_logs_log3() const {

   const double result =
      -298.6666666666667;

   return result;
}

}        // hierarchies
}        // himalaya
