// ====================================================================
// This file is part of Himalaya.
//
// Himalaya is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================

#include "H32q2g.hpp"
#include "enums.hpp"
#include "constants.hpp"
#include "power.hpp"
#include <cmath>

namespace himalaya{
namespace hierarchies{

/**
 * Constuctor
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
 * @param MuSUSY a double mu parameter
 * @param s2t a double 2 times the sine of the stop/sbottom quark mixing angle
 * @param mdrFlag an int 0 for DR and 1 for MDR scheme
 * @param oneLoopFlag an int flag to consider the one-loop expansion terms
 * @param twoLoopFlag an int flag to consider the two-loop expansion terms
 * @param threeLoopFlag an int flag to consider the three-loop expansion terms
 */
H32q2g::H32q2g(const ExpansionFlags_t& expansionDepth, double Al4p, double beta,
                 double Dmglst1, double Dmst12, double Dmsqst1, double lmMt, double lmMst1,
                 double Mt, double Mst1, double Mst2, double MuSUSY,
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
 * @return The diagonal (1, 1) matrix element of the Higgs mass matrix as a double for the hierarchy 'H32q2g'
 */
double H32q2g::getS1() const {
   return -(pow2(Mt)*pow2(MuSUSY)*(360*Dmst12*Mst1*s2t*pow2(Mst2)*(40*Al4p*Mst1*
        twoLoopFlag*(Dmst12*Mst1*s2t*(8*Dmglst1*(1 - 3*lmMst1)*Mst1 + (22 - 12*
        lmMst1)*pow2(Dmglst1) + 3*pow2(Mst1)) + 4*Mt*(-6*pow2(Dmglst1)*(Dmst12
        - pow2(Mst2)) + Dmglst1*(5 - 6*lmMst1)*Mst1*pow2(Mst2) - pow2(Mst1)*(
        Dmst12 + 3*(1 + 2*lmMst1)*pow2(Mst2)))) + 32*Al4p*(Dmst12 - pow2(Mst2))
        *(75*Al4p*s2t*(shiftst1 - shiftst2)*threeLoopFlag*xDmsqst1*pow2(
        Dmsqst1) + 8*Mt*twoLoopFlag*xDmglst1*pow3(Dmglst1)) + 45*Dmst12*
        oneLoopFlag*s2t*pow4(Mst1)) - 4*Mst1*xDmst12*pow3(Dmst12)*(Al4p*(40*
        Mst1*Mt*(-2160*s2t*twoLoopFlag*pow2(Dmglst1) + Al4p*threeLoopFlag*(s2t*
        (25200*Dmsqst1 + pow2(Dmglst1)*(5725086 + 67536*lmMst1 - 69120*pow2(
        lmMst1))) + 5*Dmglst1*Mt*(403559 + 384*(lmMst1 + lmMt) - 4608*pow2(
        lmMst1)))) - 15*pow2(Mst1)*(960*Dmglst1*s2t*((-5 + 6*lmMst1)*Mt +
        Dmglst1*(-11 + 6*lmMst1)*s2t)*twoLoopFlag + Al4p*threeLoopFlag*(16*s2t*
        (150*Dmsqst1*s2t*(-7 + 24*lmMst1*(-1 + shiftst2) - 36*shiftst2) +
        Dmglst1*Mt*(1282471 - 7264*lmMst1 - 18120*pow2(lmMst1)) + s2t*pow2(
        Dmglst1)*(655743 + 5288*lmMst1 - 11820*pow2(lmMst1))) + (1763661 -
        47104*lmMst1 + 5120*lmMt - 24576*pow2(lmMst1))*pow2(Mt))) + 32*(720*Mt*
        s2t*twoLoopFlag*xDmglst1*pow3(Dmglst1) + Al4p*threeLoopFlag*(7500*
        Dmglst1*Dmsqst1*Mt*s2t + pow2(Dmglst1)*(6233611 + 58800*lmMst1 - 4560*
        lmMt - 14400*pow2(lmMst1))*pow2(Mt) + 6750*(shiftst1 - shiftst2)*
        xDmsqst1*pow2(Dmsqst1)*pow2(s2t) + Mt*s2t*xDmglst1*(14217821 + 161940*
        lmMst1 - 28800*pow2(lmMst1))*pow3(Dmglst1))) - 600*s2t*(24*((5 + 6*
        lmMst1)*Mt + 4*Dmglst1*(-1 + 3*lmMst1)*s2t)*twoLoopFlag + Al4p*
        threeLoopFlag*(Dmglst1*s2t*(84209 + 1264*lmMst1 - 240*pow2(lmMst1)) -
        2*Mt*(-36863 + 80*lmMst1 + 552*pow2(lmMst1))))*pow3(Mst1)) + 15*(270*
        oneLoopFlag + Al4p*(240*(twoLoopFlag + 6*lmMst1*twoLoopFlag) + Al4p*
        threeLoopFlag*(350605 + 4320*shiftst1 + 2880*shiftst2 + 8352*shiftst3 -
        96*lmMst1*(-115 + 90*shiftst1 + 60*shiftst2 + 54*shiftst3) - 2160*pow2(
        lmMst1))))*pow2(s2t)*pow4(Mst1)) - threeLoopFlag*pow2(Al4p)*(128*Mt*
        xDmglst1*pow3(Dmglst1)*(Mst1*s2t*(Dmst12*(-14217821 - 161940*lmMst1 +
        11852775*z3 + 28800*pow2(lmMst1))*(Dmst12 - pow2(Mst2))*pow2(Mst2) -
        11852775*xDmst12*z3*pow3(Dmst12)) + Mt*(4438798 + 19200*lmMst1 + 1920*
        lmMt - 3719925*z3)*pow6(Mst2)) + 15*pow2(Mst1)*(-(pow2(Dmst12)*pow2(
        Mst2)*(240*Mt*s2t*(1120*Dmsqst1 + (-20946 + 80*lmMst1 + 17285*z3 + 512*
        pow2(lmMst1))*pow2(Mst1)) - Mst1*(1705082 + 18432*lmMst1 - 20480*lmMt -
        1438335*z3 + 12288*pow2(lmMst1))*pow2(Mt) + 2400*Dmsqst1*Mst1*(14 - 24*
        lmMst1*(-2 + shiftst1 + shiftst2) + (shiftst1 + shiftst2)*(36 - 24*z2)
        - 63*z3)*pow2(s2t) - 80*(2858 - 908*lmMst1 + 720*lmMst1*(shiftst1 +
        shiftst2) - 252*shiftst3 + 216*lmMst1*shiftst3 + 72*(10*shiftst2 + 3*
        shiftst3)*z2 - 360*(shiftst1 + shiftst2 - 2*shiftst1*z2) - 3861*z3 +
        48*pow2(lmMst1))*pow2(s2t)*pow3(Mst1))) + 160*Dmst12*(120*Dmsqst1*s2t*(
        14*Mt + 3*Mst1*s2t*(shiftst1 - shiftst2)*(3 - 2*(lmMst1 + z2))) + 4*Mt*
        s2t*(2722 + 20*lmMst1 - 2259*z3 + 108*pow2(lmMst1))*pow2(Mst1) + Mst1*(
        22778 - 1408*lmMst1 + 384*lmMt - 19125*z3 - 768*pow2(lmMst1))*pow2(Mt)
        + 36*(10*shiftst1 - 10*shiftst2 + shiftst3)*(1 - 2*(lmMst1 + z2))*pow2(
        s2t)*pow3(Mst1))*pow4(Mst2) + Mst1*(6*xDmst12*(5*Mt*(197889*Mt +
        324752*Mst1*s2t)*z3 + (-2400*Dmsqst1*(16*shiftst2*z2 + 21*z3) - (5760*
        shiftst1*z2 + 3840*shiftst2*z2 + 3456*shiftst3*z2 + 188345*z3)*pow2(
        Mst1))*pow2(s2t))*pow3(Dmst12) - 3840*(-349 + 56*lmMst1 - 24*lmMt +
        282*z3 + 32*pow2(lmMst1))*pow2(Mt)*pow6(Mst2))) + 8*Dmglst1*Mst1*(2*
        Dmglst1*pow2(Dmst12)*pow2(Mst2)*(30*Mst1*Mt*s2t*(-1908362 - 22512*
        lmMst1 + 1581075*z3 + 23040*pow2(lmMst1)) + 4*(-12467222 - 117600*
        lmMst1 + 9120*lmMt + 10353825*z3 + 28800*pow2(lmMst1))*pow2(Mt) + 15*(
        1311486 + 10576*lmMst1 - 1122225*z3 - 23640*pow2(lmMst1))*pow2(Mst1)*
        pow2(s2t)) - 900*Dmglst1*xDmst12*z3*(105405*Mst1*Mt*s2t + 92034*pow2(
        Mt) - 74815*pow2(Mst1)*pow2(s2t))*pow3(Dmst12) - 5*pow2(Dmst12)*pow2(
        Mst2)*(Mt*s2t*(24000*Dmsqst1 + 6*(-1899722 + 3888*lmMst1 + 1581075*z3 +
        23040*pow2(lmMst1))*pow2(Mst1)) + 10*Mst1*(807118 + 768*(lmMst1 + lmMt)
        - 665253*z3 - 9216*pow2(lmMst1))*pow2(Mt) + 15*(-168418 - 2528*lmMst1 +
        151455*z3 + 480*pow2(lmMst1))*pow2(s2t)*pow3(Mst1)) + 4*Dmglst1*Dmst12*
        Mt*(15*Mst1*s2t*(1908362 + 22512*lmMst1 - 1581075*z3 - 23040*pow2(
        lmMst1)) - 2*Mt*(-12467222 - 117600*lmMst1 + 9120*lmMt + 10353825*z3 +
        28800*pow2(lmMst1)))*pow4(Mst2) + 50*Dmst12*Mt*(Mst1*Mt*(807118 + 768*(
        lmMst1 + lmMt) - 665253*z3 - 9216*pow2(lmMst1)) + s2t*(2400*Dmsqst1 -
        12*(-66522 + 1064*lmMst1 + 56373*z3 + 1320*pow2(lmMst1))*pow2(Mst1)))*
        pow4(Mst2) - 100*Dmglst1*(-460018 - 9408*lmMst1 + 960*lmMt + 381915*z3
        + 2304*pow2(lmMst1))*pow2(Mt)*pow6(Mst2) - 5*Mst1*(90*xDmst12*z3*(-
        285974*Mst1*Mt*s2t + 73917*pow2(Mt) - 50485*pow2(Mst1)*pow2(s2t))*pow3(
        Dmst12) + 240*(-19262 - 32*lmMst1 - 96*lmMt + 15741*z3 + 384*pow2(
        lmMst1))*pow2(Mt)*pow6(Mst2))))))/(777600.*pow5(Mst1)*pow6(Mst2));
}

/**
 * @return The diagonal (2, 2) matrix element of the Higgs mass matrix as a double for the hierarchy 'H32q2g'
 */
double H32q2g::getS2() const {
   return -(-4*xDmst12*pow2(Mst1)*pow3(Dmst12)*(Mt*(20250*oneLoopFlag*(-(Mt*Tbeta*
        pow2(s2t)*(pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 18*pow2(Mst1)*pow2(Sbeta))
        ) - 8*Tbeta*pow2(Sbeta)*pow3(Mt) + MuSUSY*pow2(Sbeta)*(8*s2t*pow2(Mt) -
        pow2(Mst1)*pow3(s2t)))*pow4(Mst1) + Al4p*((16*xDmglst1*pow3(Dmglst1)*(
        8820*twoLoopFlag*(10*s2t*pow2(Mt)*(-4*Tbeta*pow2(MuSUSY)*(-1 + pow2(
        Sbeta)) + Mst1*(66*MuSUSY + 17*Mst1*Tbeta)*pow2(Sbeta)) + pow2(Sbeta)*(
        -60*Mt*(2*MuSUSY + 11*Mst1*Tbeta)*pow2(Mst1)*pow2(s2t) + ((-181 + 90*
        lmMst1 - 90*lmMt)*MuSUSY - 800*Mst1*Tbeta)*pow3(Mt) + 20*Tbeta*pow3(
        s2t)*pow4(Mst1))) + Al4p*threeLoopFlag*(-2*s2t*pow2(Mt)*(2*Mst1*(3675*
        MuSUSY*(520781 + 17172*lmMst1 - 192*lmMt - 11520*pow2(lmMst1)) + Mst1*
        Tbeta*(-15951337379 + 10326120*lmMt + 120*lmMst1*(-1097288 + 945*lmMt)
        + 21054600*pow2(lmMst1)))*pow2(Sbeta) + 245*Tbeta*(-((-14217821 -
        161940*lmMst1 + 28800*pow2(lmMst1))*pow2(MuSUSY)*(-1 + pow2(Sbeta))) +
        48*Dmsqst1*(1331 - 420*lmMst1 + 420*lmMt)*pow2(Sbeta))) + 5*pow2(Sbeta)
        *(-294*Mt*(MuSUSY*(14217821 + 161940*lmMst1 - 28800*pow2(lmMst1)) + 10*
        Mst1*Tbeta*(-520877 - 17172*lmMst1 + 192*lmMt + 11520*pow2(lmMst1)))*
        pow2(Mst1)*pow2(s2t) + 2*(4*Mst1*Tbeta*(-16826654 + 1921395*lmMt + 3*
        lmMst1*(555463 + 2520*lmMt) - 4241160*pow2(lmMst1)) + 3*MuSUSY*(-
        1417174939 + 401268*lmMt - 4*lmMst1*(4061413 + 37800*lmMt) + 3185280*
        pow2(lmMst1) - 211680*pow2(lmMt)))*pow3(Mt) + 49*Tbeta*(14217821 +
        161940*lmMst1 - 28800*pow2(lmMst1))*pow3(s2t)*pow4(Mst1)))))/49. + 720*
        Mst1*twoLoopFlag*(MuSUSY*pow2(Sbeta)*(25*pow2(Mst1)*(96*(-1 + lmMst1)*
        Mst1*s2t*pow2(Mt) - 12*(1 + 3*lmMst1)*Mt*pow2(Mst1)*pow2(s2t) + 8*(5 +
        6*lmMst1)*pow3(Mt) - 3*pow3(Mst1)*pow3(s2t)) - 2*Dmglst1*Mst1*(-3600*
        lmMst1*Mst1*s2t*pow2(Mt) + 75*(-5 + 6*lmMst1)*Mt*pow2(Mst1)*pow2(s2t) +
        24*(17 - 30*lmMst1 + 5*lmMt)*pow3(Mt) + 100*(1 - 3*lmMst1)*pow3(Mst1)*
        pow3(s2t)) + pow2(Dmglst1)*(240*(-29 + 15*lmMst1)*Mst1*s2t*pow2(Mt) +
        1800*Mt*pow2(Mst1)*pow2(s2t) + (3924 - 360*lmMst1 + 360*lmMt)*pow3(Mt)
        + 50*(-11 + 6*lmMst1)*pow3(Mst1)*pow3(s2t))) + Tbeta*(25*Mt*s2t*(-4*
        Dmglst1*Mst1*((-5 + 6*lmMst1)*Mt + 4*(-1 + 3*lmMst1)*Mst1*s2t) - 4*(6*
        Mt + (-11 + 6*lmMst1)*Mst1*s2t)*pow2(Dmglst1) + (-4*(5 + 6*lmMst1)*Mt +
        (1 + 6*lmMst1)*Mst1*s2t)*pow2(Mst1))*pow2(MuSUSY) + pow2(Sbeta)*(-50*(
        Dmglst1*(5 - 6*lmMst1)*Mst1 + 6*pow2(Dmglst1) - 3*(1 + 2*lmMst1)*pow2(
        Mst1))*pow3(s2t)*pow4(Mst1) - 4*s2t*pow2(Mt)*(2*(161 + 60*lmMst1 - 60*
        lmMt)*pow2(Dmglst1)*pow2(Mst1) - 25*(Dmglst1*(-5 + 6*lmMst1)*Mst1 + 6*
        pow2(Dmglst1) + (5 + 6*lmMst1)*pow2(Mst1))*pow2(MuSUSY) + 2*Dmglst1*(-
        131 + 165*lmMst1 + 135*lmMt)*pow3(Mst1) + 50*(1 + 3*lmMst1 + 9*lmMt)*
        pow4(Mst1)) + Mst1*(4*(-2*Dmglst1*(47 + 870*lmMst1 + 30*lmMt)*Mst1 + (
        1641 - 990*lmMst1 + 90*lmMt)*pow2(Dmglst1) + (137 - 330*lmMst1 - 270*
        lmMt)*pow2(Mst1))*pow3(Mt) - 25*Mt*pow2(s2t)*((1 + 6*lmMst1)*pow2(Mst1)
        *pow2(MuSUSY) + pow2(Dmglst1)*((128.4 - 72*lmMst1)*pow2(Mst1) + 4*(11 -
        6*lmMst1)*pow2(MuSUSY)) + 4*Dmglst1*(4*(1 - 3*lmMst1)*Mst1*pow2(MuSUSY)
        + (9 - 36*lmMst1)*pow3(Mst1)) + 12*(1 + 12*lmMst1)*pow4(Mst1)))))))) +
        threeLoopFlag*pow2(Al4p)*((8*Mt*pow2(Dmglst1)*(2*Mt*pow2(Mst1)*(Mt*(
        44100*MuSUSY*s2t*(274009 + 964*lmMst1 + 104*lmMt - 8310*pow2(lmMst1)) +
        Mt*Tbeta*(-6321826673 + 4506600*lmMt + 60*lmMst1*(1009697 + 158550*
        lmMt) + 357663600*pow2(lmMst1) + 6350400*pow2(lmMt)))*pow2(Sbeta) +
        3675*Tbeta*pow2(s2t)*(-((-655743 - 5288*lmMst1 + 11820*pow2(lmMst1))*
        pow2(MuSUSY)*(-1 + pow2(Sbeta))) + 1080*Dmsqst1*pow2(Sbeta))) + pow2(
        Mt)*(980*(Mt*Tbeta*(-6233611 - 58800*lmMst1 + 4560*lmMt + 14400*pow2(
        lmMst1))*pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 36*Dmsqst1*(225*MuSUSY*s2t +
        (-457 - 510*lmMst1 + 510*lmMt)*Mt*Tbeta)*pow2(Sbeta)) - 294*Mst1*(2*Mt*
        MuSUSY*(84334067 + 120*lmMst1*(2843 - 120*lmMt) + 202200*lmMt - 828000*
        pow2(lmMst1) - 21600*pow2(lmMt))*pow2(Sbeta) - 25*s2t*Tbeta*((-954181 -
        11256*lmMst1 + 11520*pow2(lmMst1))*pow2(MuSUSY)*(-1 + pow2(Sbeta)) +
        48*Dmsqst1*(397 - 30*lmMst1 + 30*lmMt)*pow2(Sbeta)))) + 49*s2t*(450*Mt*
        MuSUSY*s2t*(-954181 - 11256*lmMst1 + 11520*pow2(lmMst1)) + 75*Mst1*s2t*
        (MuSUSY*s2t*(655743 + 5288*lmMst1 - 11820*pow2(lmMst1)) + Mst1*s2t*
        Tbeta*(954181 + 11256*lmMst1 - 11520*pow2(lmMst1)) + 3*Mt*Tbeta*(-
        452211 - 17060*lmMst1 + 992*lmMt + 1320*pow2(lmMst1))) + 8*Tbeta*(
        192278911 + 177300*lmMt + 60*lmMst1*(19139 + 570*lmMt) - 1373400*pow2(
        lmMst1) + 43200*pow2(lmMt))*pow2(Mt))*pow2(Sbeta)*pow3(Mst1)))/49. -
        30000*xDmsqst1*pow2(Dmsqst1)*(-3*Tbeta*pow2(Mt)*pow2(s2t)*(-12*(
        shiftst1 - shiftst2)*pow2(MuSUSY)*(-1 + pow2(Sbeta)) - (205 + 72*(
        shiftst1 + shiftst2))*pow2(Mst1)*pow2(Sbeta)) + pow2(Sbeta)*(-(s2t*(
        MuSUSY*(615 + 432*shiftst2) + 8*(Dmglst1*(19 + 6*lmMst1 - 6*lmMt) + 4*(
        13 - 3*lmMst1 + 3*lmMt)*Mst1)*Tbeta)*pow3(Mt)) - 18*(shiftst1 -
        shiftst2)*pow2(Mst1)*(-6*Mt*MuSUSY + s2t*Tbeta*pow2(Mst1))*pow3(s2t) +
        2*(283 + 6*lmMst1 - 6*lmMt + 216*shiftst2)*Tbeta*pow4(Mt))) + pow2(
        Mst1)*(-(pow2(Mst1)*pow2(Mt)*(-75*Tbeta*(350605 + 4320*shiftst1 + 2880*
        shiftst2 + 8352*shiftst3 - 96*lmMst1*(-115 + 90*shiftst1 + 60*shiftst2
        + 54*shiftst3) - 2160*pow2(lmMst1))*pow2(MuSUSY)*pow2(s2t)*(1 - pow2(
        Sbeta)) + 4*Mt*(-225*MuSUSY*s2t*(-102747 + 640*lmMt + 6720*shiftst3 -
        32*lmMst1*(331 + 90*shiftst3) + 13888*pow2(lmMst1)) + Mt*Tbeta*(-
        20857591 + 2151960*lmMt - 120*lmMst1*(36107 + 13380*lmMt - 5400*
        shiftst3) - 1512000*shiftst3 + 2318400*pow2(lmMst1) + 1663200*pow2(
        lmMt)))*pow2(Sbeta))) + 60*Mt*s2t*(-75*Mt*s2t*(2*Mst1*Tbeta*(2785 +
        384*lmMt + 960*shiftst1 - 240*shiftst2 + 996*shiftst3 - 8*lmMst1*(38 +
        240*shiftst1 - 60*shiftst2 + 87*shiftst3) + 768*pow2(lmMst1)) + MuSUSY*
        (-20531 + 200*lmMst1 + 1200*pow2(lmMst1))) - 2*Tbeta*(558619 + 76160*
        lmMt - 224*lmMst1*(1219 + 60*lmMt) + 123840*pow2(lmMst1) + 86400*pow2(
        lmMt))*pow2(Mt) + 50*Mst1*MuSUSY*(1429 - 720*shiftst1 + 360*shiftst2 -
        234*shiftst3 + 2*lmMst1*(-227 + 720*shiftst1 - 360*shiftst2 + 126*
        shiftst3) + 24*pow2(lmMst1))*pow2(s2t))*pow2(Sbeta)*pow3(Mst1) + 240*
        Mst1*MuSUSY*(-25*MuSUSY*s2t*Tbeta*(-36863 + 80*lmMst1 + 552*pow2(
        lmMst1))*(-1 + pow2(Sbeta)) - Mt*(-3454599 + 16840*lmMt + 48*lmMst1*(
        262 + 405*lmMt) + 46560*pow2(lmMst1))*pow2(Sbeta))*pow3(Mt) + Tbeta*(-
        75*(-1763661 + 47104*lmMst1 - 5120*lmMt + 24576*pow2(lmMst1))*pow2(
        MuSUSY)*(-1 + pow2(Sbeta))*pow4(Mt) + 6000*(9*(1 - 2*lmMst1)*Mst1*s2t*(
        10*shiftst1 - 10*shiftst2 + shiftst3) + 2*Mt*(1361 + 10*lmMst1 + 54*
        pow2(lmMst1)))*pow2(Sbeta)*pow3(s2t)*pow5(Mst1))) + 600*Dmsqst1*Mst1*(-
        75*Mst1*pow2(Mt)*pow2(s2t)*(-4*(-7 + 24*lmMst1*(-1 + shiftst2) - 36*
        shiftst2)*Tbeta*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 336*Mst1*MuSUSY*pow2(
        Sbeta) - (137 - 576*shiftst1 + 144*shiftst2 - 96*lmMst1*(3 - 4*shiftst1
        + shiftst2))*Tbeta*pow2(Mst1)*pow2(Sbeta)) - 100*s2t*(84*Tbeta*pow2(
        MuSUSY)*(-1 + pow2(Sbeta)) + Mst1*(MuSUSY*(519 - 432*lmMst1*(-1 +
        shiftst2) + 648*shiftst2) - 8*(13 + 6*lmMst1 - 6*lmMt)*Mst1*Tbeta)*
        pow2(Sbeta))*pow3(Mt) + pow2(Sbeta)*(-150*Mt*(MuSUSY*(7 + 72*shiftst1 -
        36*shiftst2 + 24*lmMst1*(1 - 2*shiftst1 + shiftst2)) - 28*Mst1*Tbeta)*
        pow3(Mst1)*pow3(s2t) - 16*(50*(65 - 6*lmMst1 + 6*lmMt)*MuSUSY + Mst1*(-
        3268 + 105*lmMt - 4050*shiftst2 + 15*lmMst1*(-187 + 180*shiftst2))*
        Tbeta)*pow4(Mt) + 900*(3 - 2*lmMst1)*(shiftst1 - shiftst2)*Tbeta*pow4(
        s2t)*pow5(Mst1))) + 4*Dmglst1*Mt*(1200*Dmsqst1*(-2*s2t*pow2(Mt)*(125*
        Tbeta*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 6*Mst1*(1175*MuSUSY + (47 - 30*
        lmMst1 + 30*lmMt)*Mst1*Tbeta)*pow2(Sbeta)) + pow2(Sbeta)*(150*Mt*(-5*
        MuSUSY + 94*Mst1*Tbeta)*pow2(Mst1)*pow2(s2t) - 4*((557 + 120*lmMst1 -
        120*lmMt)*MuSUSY + 9*(-423 + 20*lmMst1 - 20*lmMt)*Mst1*Tbeta)*pow3(Mt)
        + 125*Tbeta*pow3(s2t)*pow4(Mst1))) + Mst1*(6*Mt*pow2(Mst1)*(-125*Tbeta*
        (-84209 - 1264*lmMst1 + 240*pow2(lmMst1))*pow2(MuSUSY)*pow2(s2t)*(-1 +
        pow2(Sbeta)) + 4*Mt*(150*MuSUSY*s2t*(12383 + 4128*lmMst1 + 80*lmMt -
        1260*pow2(lmMst1)) + Mt*Tbeta*(9598037 + 92280*lmMt + 20*lmMst1*(-11207
        + 270*lmMt) + 246000*pow2(lmMst1) - 14400*pow2(lmMt)))*pow2(Sbeta)) +
        4*Mst1*MuSUSY*pow2(Mt)*(-75*MuSUSY*s2t*Tbeta*(-1282471 + 7264*lmMst1 +
        18120*pow2(lmMst1))*(-1 + pow2(Sbeta)) - 2*Mt*(-193364399 - 90000*lmMt
        + 300*lmMst1*(3781 + 582*lmMt) + 2005200*pow2(lmMst1) + 43200*pow2(
        lmMt))*pow2(Sbeta)) + s2t*(75*Mst1*s2t*(5*MuSUSY*s2t*(84209 + 1264*
        lmMst1 - 240*pow2(lmMst1)) + 48*Mt*Tbeta*(-20017 + 1203*lmMst1 - 200*
        lmMt + 2250*pow2(lmMst1))) + Mt*(225*MuSUSY*s2t*(284641 + 8696*lmMst1 +
        1680*pow2(lmMst1)) + 4*Mt*Tbeta*(-28188929 + 143100*lmMt + 3780*lmMst1*
        (549 + 80*lmMt) - 1389600*pow2(lmMst1) - 388800*pow2(lmMt))))*pow2(
        Sbeta)*pow3(Mst1) + 250*Tbeta*((403559 + 384*(lmMst1 + lmMt) - 4608*
        pow2(lmMst1))*pow2(MuSUSY)*(1 - pow2(Sbeta))*pow3(Mt) - 6*(-33261 +
        532*lmMst1 + 660*pow2(lmMst1))*pow2(Sbeta)*pow3(s2t)*pow5(Mst1)))))) +
        pow2(Mst2)*(160*threeLoopFlag*pow2(Al4p)*(-75*xDmsqst1*pow2(Dmsqst1)*(
        360*Tbeta*pow2(Mst1)*pow2(Mt)*(Dmst12*pow2(Mst2)*(-((shiftst1 -
        shiftst2)*pow2(MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta))) - 12*(shiftst2*
        pow2(Mt) + (shiftst1 - shiftst2)*pow2(Mst1)*pow2(s2t))*pow2(Sbeta)) +
        pow2(Dmst12)*((shiftst1 - shiftst2)*pow2(MuSUSY)*pow2(s2t)*(-1 + pow2(
        Sbeta)) + 3*(4*shiftst2*pow2(Mt) + (3*shiftst1 - shiftst2)*pow2(Mst1)*
        pow2(s2t))*pow2(Sbeta)) - 12*(shiftst1 + shiftst2)*pow2(Mt)*pow2(Sbeta)
        *pow4(Mst2)) + pow2(Sbeta)*(10*Mt*MuSUSY*(16*Dmglst1*(8 + 3*lmMst1 - 3*
        lmMt)*pow3(Mt)*pow4(Mst2) + Mst1*(3*Mst1*s2t*pow2(Dmst12)*(-((205 +
        144*shiftst2)*pow2(Mt)) + 18*(shiftst1 - shiftst2)*pow2(Mst1)*pow2(s2t)
        ) + pow2(Mt)*(3*Dmst12*Mst1*s2t*(205 + 144*shiftst2)*pow2(Mst2) + 16*((
        29 - 6*lmMst1 + 6*lmMt)*Mt + 27*Mst1*s2t*(-shiftst1 + shiftst2))*pow4(
        Mst2)))) + Tbeta*(3075*pow2(Dmst12)*pow2(Mt)*pow2(s2t)*pow4(Mst1) + 20*
        pow2(Mst1)*pow3(Mt)*(((283 + 6*lmMst1 - 6*lmMt)*Mt - 4*(Dmglst1*(19 +
        6*lmMst1 - 6*lmMt) + 4*(13 - 3*lmMst1 + 3*lmMt)*Mst1)*s2t)*pow2(Dmst12)
        + Dmst12*((-283 - 6*lmMst1 + 6*lmMt)*Mt + 4*(Dmglst1*(19 + 6*lmMst1 -
        6*lmMt) + 4*(13 - 3*lmMst1 + 3*lmMt)*Mst1)*s2t)*pow2(Mst2) - 2*Mt*(299
        + 42*lmMt - 6*lmMst1*(7 + 18*lmMt) + 54*pow2(lmMst1) + 54*pow2(lmMt))*
        pow4(Mst2)) + (48*Dmglst1*(Dmglst1*(131 - 420*lmMst1 + 420*lmMt) + 250*
        (-22 + 3*lmMst1 - 3*lmMt)*Mst1)*pow4(Mst2)*pow4(Mt))/25.))) + Mst1*((
        Mt*MuSUSY*pow2(Sbeta)*(15*(200*pow2(Mst2)*(-6*Dmst12*Mst1*s2t*(631 +
        36*lmMt + 90*shiftst2 + 27*shiftst3 - lmMst1*(1 + 180*shiftst2 + 18*
        shiftst3) + 21*pow2(lmMst1)) + Dmst12*Mt*(11697 + 448*lmMst1 + 330*lmMt
        - 372*lmMst1*lmMt + 408*pow2(lmMst1) + 288*pow2(lmMt)) + 54*(1 - 2*
        lmMst1)*Mst1*s2t*(10*shiftst1 - 10*shiftst2 + shiftst3)*pow2(Mst2) + 8*
        Mt*(434 - 83*lmMst1 + 174*lmMt - 66*lmMst1*lmMt + 183*pow2(lmMst1) +
        108*pow2(lmMt))*pow2(Mst2))*pow2(Mt) + pow2(Dmst12)*(150*Mst1*s2t*(-
        2071 + 96*lmMt + 288*shiftst3 - 8*lmMst1*(37 + 18*shiftst3) + 600*pow2(
        lmMst1))*pow2(Mt) - 300*Mt*(1361 + 10*lmMst1 + 54*pow2(lmMst1))*pow2(
        Mst1)*pow2(s2t) - (-2284899 + 49840*lmMt - 32*lmMst1*(-1793 + 555*lmMt)
        + 87360*pow2(lmMst1) + 28800*pow2(lmMt))*pow3(Mt) + 1350*(-1 + 2*
        lmMst1)*(10*shiftst1 - 10*shiftst2 + shiftst3)*pow3(Mst1)*pow3(s2t)))*
        pow4(Mst1) - Mt*pow2(Dmglst1)*(3*Dmst12*Mst1*(-4*Mt*(Mst1*Mt*(84334067
        + 120*lmMst1*(2843 - 120*lmMt) + 202200*lmMt - 828000*pow2(lmMst1) -
        21600*pow2(lmMt)) + 75*s2t*(180*Dmsqst1 + (32101 - 5044*lmMst1 + 400*
        lmMt - 5100*pow2(lmMst1))*pow2(Mst1)))*pow2(Mst2) + Dmst12*(150*Mt*s2t*
        (360*Dmsqst1 + (-515917 - 6972*lmMst1 + 192*lmMt + 11520*pow2(lmMst1))*
        pow2(Mst1)) + 4*Mst1*(84334067 + 120*lmMst1*(2843 - 120*lmMt) + 202200*
        lmMt - 828000*pow2(lmMst1) - 21600*pow2(lmMt))*pow2(Mt) + 75*(954181 +
        11256*lmMst1 - 11520*pow2(lmMst1))*pow2(s2t)*pow3(Mst1))) + 160*(90*
        Dmsqst1*(206 - 15*lmMst1 + 15*lmMt) + (-3030652 - 19305*lmMt + 6*
        lmMst1*(1183 + 645*lmMt) + 55530*pow2(lmMst1) + 5400*pow2(lmMt))*pow2(
        Mst1))*pow2(Mt)*pow4(Mst2)) + 3750*Dmsqst1*pow2(Mst1)*(-(pow2(Dmst12)*(
        483*Mst1*s2t*pow2(Mt) + 252*Mt*pow2(Mst1)*pow2(s2t) + 16*(65 - 6*lmMst1
        + 6*lmMt)*pow3(Mt) + 54*(3 - 2*lmMst1)*(shiftst1 - shiftst2)*pow3(Mst1)
        *pow3(s2t))) + pow2(Mt)*(8*Dmst12*(2*(65 - 6*lmMst1 + 6*lmMt)*Mt - 9*
        Mst1*s2t*(1 - 12*lmMst1*(-1 + shiftst2) + 18*shiftst2))*pow2(Mst2) +
        48*((40 - 6*lmMst1 + 6*lmMt)*Mt + 9*(3 - 2*lmMst1)*Mst1*s2t*(shiftst1 -
        shiftst2))*pow4(Mst2))) + Dmglst1*Mt*(600*Dmsqst1*Mst1*(4*Dmst12*Mt*((
        557 + 120*lmMst1 - 120*lmMt)*Mt + 3525*Mst1*s2t)*pow2(Mst2) - pow2(
        Dmst12)*(14100*Mst1*Mt*s2t + 4*(557 + 120*lmMst1 - 120*lmMt)*pow2(Mt) +
        375*pow2(Mst1)*pow2(s2t)) + 200*(26 + 3*lmMst1 - 3*lmMt)*pow2(Mt)*pow4(
        Mst2)) + 2*pow3(Mst1)*(20*Dmst12*Mt*(225*Mst1*s2t*(-4539 + 596*lmMst1 -
        48*lmMt + 516*pow2(lmMst1)) - 2*Mt*(-3746977 - 4005*lmMt + 18*lmMst1*(
        2711 + 1215*lmMt) + 52380*pow2(lmMst1)))*pow2(Mst2) + 3*pow2(Dmst12)*(
        300*Mst1*Mt*s2t*(17539 + 574*lmMst1 + 160*lmMt - 1920*pow2(lmMst1)) + (
        39474953 + 3300*lmMt + 20*lmMst1*(-2639 + 4380*lmMt) - 319200*pow2(
        lmMst1) - 14400*pow2(lmMt))*pow2(Mt) + 375*(-33261 + 532*lmMst1 + 660*
        pow2(lmMst1))*pow2(Mst1)*pow2(s2t)) + 6000*(9961 + 66*lmMt - 10*lmMst1*
        (67 + 24*lmMt) + 42*pow2(lmMst1) + 72*pow2(lmMt))*pow2(Mt)*pow4(Mst2)))
        ))/5. + Mst1*Tbeta*(-(pow2(Mt)*pow2(MuSUSY)*(40*Dmglst1*(10*Dmst12*Mt*(
        Mst1*Mt*(403559 + 384*(lmMst1 + lmMt) - 4608*pow2(lmMst1)) + s2t*(1200*
        Dmsqst1 - 12*(-33261 + 532*lmMst1 + 660*pow2(lmMst1))*pow2(Mst1)))*
        pow2(Mst2) - pow2(Dmst12)*(Mt*s2t*(12000*Dmsqst1 + 6*(-949861 + 1944*
        lmMst1 + 11520*pow2(lmMst1))*pow2(Mst1)) + 10*Mst1*(403559 + 384*(
        lmMst1 + lmMt) - 4608*pow2(lmMst1))*pow2(Mt) + 15*(-84209 - 1264*lmMst1
        + 240*pow2(lmMst1))*pow2(s2t)*pow3(Mst1)) + 240*Mst1*(9631 + 16*lmMst1
        + 48*lmMt - 192*pow2(lmMst1))*pow2(Mt)*pow4(Mst2)) + 15*Mst1*(160*
        Dmst12*Mt*(840*Dmsqst1*s2t + Mst1*Mt*(11389 - 704*lmMst1 + 192*lmMt -
        384*pow2(lmMst1)) + 4*s2t*(1361 + 10*lmMst1 + 54*pow2(lmMst1))*pow2(
        Mst1))*pow2(Mst2) + pow2(Dmst12)*(240*Mt*s2t*(-560*Dmsqst1 + (10473 -
        40*lmMst1 - 256*pow2(lmMst1))*pow2(Mst1)) + Mst1*(852541 + 9216*lmMst1
        - 10240*lmMt + 6144*pow2(lmMst1))*pow2(Mt) + 80*Mst1*(-30*Dmsqst1*(7 +
        24*lmMst1) + (1429 - 454*lmMst1 + 24*pow2(lmMst1))*pow2(Mst1))*pow2(
        s2t)) + 1920*Mst1*(349 - 56*lmMst1 + 24*lmMt - 32*pow2(lmMst1))*pow2(
        Mt)*pow4(Mst2)) + 16*pow2(Dmglst1)*(2*Dmst12*Mt*(15*Mst1*s2t*(954181 +
        11256*lmMst1 - 11520*pow2(lmMst1)) - 2*Mt*(-6233611 - 58800*lmMst1 +
        4560*lmMt + 14400*pow2(lmMst1)))*pow2(Mst2) + pow2(Dmst12)*(30*Mst1*Mt*
        s2t*(-954181 - 11256*lmMst1 + 11520*pow2(lmMst1)) + 4*(-6233611 -
        58800*lmMst1 + 4560*lmMt + 14400*pow2(lmMst1))*pow2(Mt) + 15*(655743 +
        5288*lmMst1 - 11820*pow2(lmMst1))*pow2(Mst1)*pow2(s2t)) - 50*(-230009 -
        4704*lmMst1 + 480*lmMt + 1152*pow2(lmMst1))*pow2(Mt)*pow4(Mst2))))/16.
         + pow2(Mst1)*pow2(Mt)*((15*pow2(Dmst12)*pow2(s2t)*pow2(Sbeta)*(2*pow2(
        Dmglst1)*(1080*Dmsqst1 - 6*(-31901 + 5044*lmMst1 - 400*lmMt + 5100*
        pow2(lmMst1))*pow2(Mst1) + (655743 + 5288*lmMst1 - 11820*pow2(lmMst1))*
        pow2(MuSUSY)) + 5*Dmglst1*Mst1*(22560*Dmsqst1 + 24*(-4515 + 596*lmMst1
        - 48*lmMt + 516*pow2(lmMst1))*pow2(Mst1) + (84209 + 1264*lmMst1 - 240*
        pow2(lmMst1))*pow2(MuSUSY)) - 10*((-1429 + 454*lmMst1 - 24*pow2(lmMst1)
        )*pow2(Mst1)*pow2(MuSUSY) + 30*Dmsqst1*(12*(1 + 12*lmMst1)*pow2(Mst1) +
        (7 + 24*lmMst1)*pow2(MuSUSY)) + 24*(613 - lmMst1 + 36*lmMt + 21*pow2(
        lmMst1))*pow4(Mst1))))/2. + 13500*(Dmsqst1*(3 - 2*lmMst1) + (1 - 2*
        lmMst1)*pow2(Mst1))*(pow2(Dmst12)*pow2(s2t)*(-((shiftst1 + shiftst2)*
        pow2(MuSUSY)*(-1 + pow2(Sbeta))) + 6*(-3*shiftst1 + shiftst2)*pow2(
        Mst1)*pow2(Sbeta)) + 2*Dmst12*pow2(Mst2)*((shiftst1 - shiftst2)*pow2(
        MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta)) + 12*(shiftst2*pow2(Mt) + (
        shiftst1 - shiftst2)*pow2(Mst1)*pow2(s2t))*pow2(Sbeta)) + 24*(shiftst1
        + shiftst2)*pow2(Mt)*pow2(Sbeta)*pow4(Mst2))) + (pow2(Sbeta)*pow3(Mt)*(
        784*Dmst12*s2t*(2*Mst1*pow2(Dmglst1)*(3600*Dmsqst1*(397 - 30*lmMst1 +
        30*lmMt)*(Dmst12 - pow2(Mst2)) + 5*pow2(Mst2)*(-8*(-6041999 - 49410*
        lmMt + 6*lmMst1*(1721 + 1290*lmMt) + 111060*pow2(lmMst1) + 10800*pow2(
        lmMt))*pow2(Mst1) + 15*(954181 + 11256*lmMst1 - 11520*pow2(lmMst1))*
        pow2(MuSUSY)) + Dmst12*((263717842 - 633600*lmMt + 2400*lmMst1*(1043 +
        93*lmMt) - 525600*pow2(lmMst1) + 302400*pow2(lmMt))*pow2(Mst1) + 75*(-
        954181 - 11256*lmMst1 + 11520*pow2(lmMst1))*pow2(MuSUSY))) + 5*Dmglst1*
        (240*Dmsqst1*(25*pow2(Mst2)*(2*(55 + 6*lmMst1 - 6*lmMt)*pow2(Mst1) + 5*
        pow2(MuSUSY)) - Dmst12*(4*(379 + 15*lmMst1 - 15*lmMt)*pow2(Mst1) + 125*
        pow2(MuSUSY))) + pow2(Mst1)*(300*pow2(Mst2)*(16*(4964 - 3*lmMt - 5*
        lmMst1*(55 + 24*lmMt) + 21*pow2(lmMst1) + 36*pow2(lmMt))*pow2(Mst1) + (
        33261 - 532*lmMst1 - 660*pow2(lmMst1))*pow2(MuSUSY)) + Dmst12*(8*(
        4511549 + 9810*lmMt + 6*lmMst1*(14879 + 4710*lmMt) - 117360*pow2(
        lmMst1) - 21600*pow2(lmMt))*pow2(Mst1) - 15*(-949861 + 1944*lmMst1 +
        11520*pow2(lmMst1))*pow2(MuSUSY)))) - 375*(-80*Dmsqst1*Mst1*(Dmst12*((-
        98 + 24*lmMst1 - 24*lmMt)*pow2(Mst1) - 21*pow2(MuSUSY)) + 3*pow2(Mst2)*
        ((74 - 12*lmMst1 + 12*lmMt)*pow2(Mst1) + 7*pow2(MuSUSY))) + (-8*pow2(
        Mst2)*(8*(347 - 50*lmMst1 + 66*lmMt - 66*lmMst1*lmMt + 183*pow2(lmMst1)
        + 108*pow2(lmMt))*pow2(Mst1) + (1361 + 10*lmMst1 + 54*pow2(lmMst1))*
        pow2(MuSUSY)) + Dmst12*(16*(-4378 - 517*lmMst1 + 243*lmMt - 78*lmMst1*
        lmMt + 528*pow2(lmMst1) + 288*pow2(lmMt))*pow2(Mst1) + 3*(-10473 + 40*
        lmMst1 + 256*pow2(lmMst1))*pow2(MuSUSY)))*pow3(Mst1))) + Mt*(784*
        Dmglst1*Mst1*(5*Dmst12*pow2(Mst2)*(8*(5136871 + 29790*lmMt - 12*lmMst1*
        (8942 + 4305*lmMt) - 86040*pow2(lmMst1) - 21600*pow2(lmMt))*pow2(Mst1)
        + 25*(403559 + 384*(lmMst1 + lmMt) - 4608*pow2(lmMst1))*pow2(MuSUSY)) -
        pow2(Dmst12)*(2*(22574599 + 21060*lmMt - 60*lmMst1*(6677 + 8880*lmMt) -
        1598400*pow2(lmMst1) - 172800*pow2(lmMt))*pow2(Mst1) + 125*(403559 +
        384*(lmMst1 + lmMt) - 4608*pow2(lmMst1))*pow2(MuSUSY)) + 3000*(-4*(-
        15571 + 852*lmMt + lmMst1*(-508 + 48*lmMt) + 978*pow2(lmMst1) + 612*
        pow2(lmMt))*pow2(Mst1) + (9631 + 16*lmMst1 + 48*lmMt - 192*pow2(lmMst1)
        )*pow2(MuSUSY))*pow4(Mst2) + 2400*Dmsqst1*(9*Dmst12*(423 - 20*lmMst1 +
        20*lmMt)*(Dmst12 - pow2(Mst2)) + 50*(-146 + 15*lmMst1 - 15*lmMt)*pow4(
        Mst2))) + 16*pow2(Dmglst1)*(196*Dmst12*pow2(Mst2)*((73908751 + 1029690*
        lmMt + 540*lmMst1*(2243 + 295*lmMt) + 707400*pow2(lmMst1) + 64800*pow2(
        lmMt))*pow2(Mst1) + 5*(6233611 + 58800*lmMst1 - 4560*lmMt - 14400*pow2(
        lmMst1))*pow2(MuSUSY)) - pow2(Dmst12)*((13564884271 + 96403020*lmMt +
        60*lmMst1*(968629 + 101640*lmMt) - 288338400*pow2(lmMst1))*pow2(Mst1) +
        980*(6233611 + 58800*lmMst1 - 4560*lmMt - 14400*pow2(lmMst1))*pow2(
        MuSUSY)) + 490*(16*(2055923 + 54*lmMst1*(1162 - 445*lmMt) + 61695*lmMt
        - 6345*pow2(lmMst1) + 12150*pow2(lmMt))*pow2(Mst1) + 25*(230009 + 4704*
        lmMst1 - 480*lmMt - 1152*pow2(lmMst1))*pow2(MuSUSY))*pow4(Mst2) +
        35280*Dmsqst1*(Dmst12*(457 + 510*lmMst1 - 510*lmMt)*(Dmst12 - pow2(
        Mst2)) - 70*(34 + 15*lmMst1 - 15*lmMt)*pow4(Mst2))) + 245*pow2(Mst1)*(
        2400*Dmst12*pow2(Mst2)*(4*(-4*lmMst1*(131 + 264*lmMt) + 264*pow2(
        lmMst1) + 21*(783 + 14*lmMt + 30*pow2(lmMt)))*pow2(Mst1) + (11389 -
        704*lmMst1 + 192*lmMt - 384*pow2(lmMst1))*pow2(MuSUSY)) + pow2(Dmst12)*
        (-64*(-93973 + 61305*lmMt - 54*lmMst1*(2337 + 1255*lmMt) + 48420*pow2(
        lmMst1) + 58050*pow2(lmMt))*pow2(Mst1) + 15*(852541 + 9216*lmMst1 -
        10240*lmMt + 6144*pow2(lmMst1))*pow2(MuSUSY)) + 9600*(3*(349 - 56*
        lmMst1 + 24*lmMt - 32*pow2(lmMst1))*pow2(MuSUSY) + 4*pow2(Mst1)*(623 -
        555*lmMt - 3*(221 + 129*lmMt)*pow2(lmMst1) - 486*pow2(lmMt) + lmMst1*(
        1066 + 843*lmMt + 756*pow2(lmMt)) + 252*pow3(lmMst1) - 621*pow3(lmMt)))
        *pow4(Mst2) + 960*Dmsqst1*((4561 + 510*lmMst1 - 510*lmMt)*pow2(Dmst12)
        + 50*Dmst12*(79 + 204*lmMst1 + 12*lmMt)*pow2(Mst2) - 1800*(14 + 6*lmMt
        - 6*lmMst1*(3 + lmMt) + 3*(pow2(lmMst1) + pow2(lmMt)))*pow4(Mst2))))))/
        3920. + 8100*shiftst3*pow4(Mst1)*((Dmst12*pow2(Mt)*pow2(s2t)*(-2*(-1 +
        2*lmMst1)*pow2(Mst2)*(pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 12*pow2(Mst1)*
        pow2(Sbeta)) + Dmst12*((-7 + 6*lmMst1)*pow2(MuSUSY)*(-1 + pow2(Sbeta))
        + 6*(-15 + 14*lmMst1)*pow2(Mst1)*pow2(Sbeta))))/6. + 4*pow2(Sbeta)*(-2*
        (-2 + lmMst1)*pow2(Dmst12) + Dmst12*(-3 + 2*lmMst1)*pow2(Mst2) + (1 -
        2*lmMst1)*pow4(Mst2))*pow4(Mt))))) + pow2(Mt)*(64*Al4p*(2*xDmglst1*
        pow3(Dmglst1)*((Al4p*Dmst12*threeLoopFlag*pow2(Mst1)*(2*Mt*pow2(Mst2)*(
        5*Mt*(4*Mst1*Tbeta*(16826654 - 1921395*lmMt - 3*lmMst1*(555463 + 2520*
        lmMt) + 4241160*pow2(lmMst1)) + 3*MuSUSY*(1417174939 - 401268*lmMt + 4*
        lmMst1*(4061413 + 37800*lmMt) - 3185280*pow2(lmMst1) + 211680*pow2(
        lmMt)))*pow2(Sbeta) + 49*s2t*(2*Mst1*(75*MuSUSY*(520781 + 17172*lmMst1
        - 192*lmMt - 11520*pow2(lmMst1)) + 4*Mst1*Tbeta*(27088246 + 5775*lmMt +
        45*lmMst1*(12571 + 270*lmMt) - 136350*pow2(lmMst1) + 16200*pow2(lmMt)))
        *pow2(Sbeta) + 5*Tbeta*(-((-14217821 - 161940*lmMst1 + 28800*pow2(
        lmMst1))*pow2(MuSUSY)*(-1 + pow2(Sbeta))) + 48*Dmsqst1*(1331 - 420*
        lmMst1 + 420*lmMt)*pow2(Sbeta)))) - Dmst12*((10*(4*Mst1*Tbeta*(16826654
        - 1921395*lmMt - 3*lmMst1*(555463 + 2520*lmMt) + 4241160*pow2(lmMst1))
        + 3*MuSUSY*(1417174939 - 401268*lmMt + 4*lmMst1*(4061413 + 37800*lmMt)
        - 3185280*pow2(lmMst1) + 211680*pow2(lmMt)))*pow2(Mt) + 735*(MuSUSY*(
        14217821 + 161940*lmMst1 - 28800*pow2(lmMst1)) + 10*Mst1*Tbeta*(-520877
        - 17172*lmMst1 + 192*lmMt + 11520*pow2(lmMst1)))*pow2(Mst1)*pow2(s2t))*
        pow2(Sbeta) + 2*Mt*s2t*(Mst1*(-7350*MuSUSY*(-520781 - 17172*lmMst1 +
        192*lmMt + 11520*pow2(lmMst1)) + Mst1*Tbeta*(-10642041163 + 11458020*
        lmMt + 60*lmMst1*(-346639 + 41580*lmMt) - 5670000*pow2(lmMst1) +
        3175200*pow2(lmMt)))*pow2(Sbeta) + 245*Tbeta*(-((-14217821 - 161940*
        lmMst1 + 28800*pow2(lmMst1))*pow2(MuSUSY)*(-1 + pow2(Sbeta))) + 48*
        Dmsqst1*(1331 - 420*lmMst1 + 420*lmMt)*pow2(Sbeta))))))/98. + 2*Al4p*
        threeLoopFlag*pow2(Mt)*(120*Dmsqst1*((1541 - 420*lmMst1 + 420*lmMt)*
        MuSUSY + (2579 + 120*lmMst1 - 120*lmMt)*Mst1*Tbeta)*pow2(Sbeta) + Mst1*
        (5*(2219399 + 9600*lmMst1 + 960*lmMt)*Tbeta*pow2(MuSUSY)*(-1 + pow2(
        Sbeta)) + 2*Mst1*(Mst1*Tbeta*(10583177 + 60630*lmMt + 540*lmMst1*(188 +
        395*lmMt) + 278100*pow2(lmMst1) - 59400*pow2(lmMt)) + MuSUSY*(54198467
        + 43950*lmMt + 180*lmMst1*(6353 + 135*lmMt) - 272700*pow2(lmMst1) +
        32400*pow2(lmMt)))*pow2(Sbeta)))*pow4(Mst2) + 45*twoLoopFlag*pow2(Mst1)
        *(2*Dmst12*Mt*pow2(Mst2)*(40*s2t*Tbeta*pow2(MuSUSY)*(-1 + pow2(Sbeta))
        + Mt*((181 - 90*lmMst1 + 90*lmMt)*MuSUSY + 800*Mst1*Tbeta)*pow2(Sbeta)
        - 6*Mst1*s2t*(110*MuSUSY + (-17 + 30*lmMst1 - 30*lmMt)*Mst1*Tbeta)*
        pow2(Sbeta)) - 2*pow2(Dmst12)*((((181 - 90*lmMst1 + 90*lmMt)*MuSUSY +
        800*Mst1*Tbeta)*pow2(Mt) + 30*(2*MuSUSY + 11*Mst1*Tbeta)*pow2(Mst1)*
        pow2(s2t))*pow2(Sbeta) - 2*Mt*s2t*(-20*Tbeta*pow2(MuSUSY)*(-1 + pow2(
        Sbeta)) + Mst1*(330*MuSUSY + (17 + 45*lmMst1 - 45*lmMt)*Mst1*Tbeta)*
        pow2(Sbeta))) + 8*((48 - 45*lmMst1 + 45*lmMt)*MuSUSY + (377 - 30*lmMst1
        + 30*lmMt)*Mst1*Tbeta)*pow2(Mt)*pow2(Sbeta)*pow4(Mst2))) - 45*
        twoLoopFlag*pow3(Mst1)*(2*MuSUSY*pow2(Sbeta)*(-25*pow2(Mst1)*(4*Dmst12*
        Mt*(2*(5 + 6*lmMst1 + 3*lmMt)*Mt - 9*(1 + 2*lmMst1)*Mst1*s2t)*pow2(
        Mst2) - pow2(Dmst12)*(24*(1 - 3*lmMst1)*Mst1*Mt*s2t + 2*(5 + 6*(lmMst1
        + lmMt))*pow2(Mt) + 9*(1 + 2*lmMst1)*pow2(Mst1)*pow2(s2t)) + 72*(1 +
        lmMst1 + lmMt)*pow2(Mt)*pow4(Mst2)) + 2*pow2(Dmglst1)*(3*Dmst12*Mt*((
        327 - 30*lmMst1 + 30*lmMt)*Mt + 50*(11 - 6*lmMst1)*Mst1*s2t)*pow2(Mst2)
        + 9*pow2(Dmst12)*(5*Mst1*Mt*s2t + (-109 + 10*lmMst1 - 10*lmMt)*pow2(Mt)
        - 25*pow2(Mst1)*pow2(s2t)) + 100*(17 - 3*lmMst1 + 3*lmMt)*pow2(Mt)*
        pow4(Mst2)) + Dmglst1*Mst1*(-300*Dmst12*Mt*((-3 + 6*lmMst1)*Mt + 2*(-1
        + 6*lmMst1)*Mst1*s2t)*pow2(Mst2) + 3*pow2(Dmst12)*(-100*Mst1*Mt*s2t +
        2*(-41 + 90*lmMst1 + 10*lmMt)*pow2(Mt) + 25*(-5 + 6*lmMst1)*pow2(Mst1)*
        pow2(s2t)) - 200*(-7 + 15*lmMst1 + 3*lmMt)*pow2(Mt)*pow4(Mst2))) +
        Tbeta*(25*Dmst12*s2t*(2*pow2(Dmglst1)*(12*Dmst12*Mt*pow2(MuSUSY) - 12*
        Mt*pow2(Mst2)*pow2(MuSUSY) - Dmst12*(-11 + 6*lmMst1)*Mst1*s2t*(pow2(
        MuSUSY)*(-1 + pow2(Sbeta)) + 6*pow2(Mst1)*pow2(Sbeta))) - 4*Dmglst1*
        Mst1*((5 - 6*lmMst1)*Mt*pow2(Mst2)*pow2(MuSUSY) + 2*Dmst12*Mst1*s2t*((-
        1 + 3*lmMst1)*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 3*(-1 + 6*lmMst1)*pow2(
        Mst1)*pow2(Sbeta))) + pow2(Mst1)*(12*(1 + 2*lmMst1)*Mt*pow2(Mst2)*pow2(
        MuSUSY) + Dmst12*(4*Mt*pow2(MuSUSY) + 3*Mst1*s2t*pow2(MuSUSY)*(-1 +
        pow2(Sbeta)) + 36*(1 + 2*lmMst1)*s2t*pow2(Sbeta)*pow3(Mst1)))) + pow2(
        Sbeta)*(-4*Dmst12*Mt*s2t*(2*pow2(Dmglst1)*(-25*pow2(Mst2)*((31 - 6*
        lmMst1 + 6*lmMt)*pow2(Mst1) + 3*pow2(MuSUSY)) + Dmst12*((307 - 105*
        lmMst1 + 105*lmMt)*pow2(Mst1) + 75*pow2(MuSUSY))) - 25*pow2(Mst1)*(
        Dmst12*(4*(1 + 3*lmMst1 + 6*lmMt)*pow2(Mst1) - pow2(MuSUSY)) - 3*pow2(
        Mst2)*(6*(1 + 2*(lmMst1 + lmMt))*pow2(Mst1) + (1 + 2*lmMst1)*pow2(
        MuSUSY))) + 25*Dmglst1*(Mst1*pow2(Mst2)*(2*(-17 + 30*lmMst1 + 6*lmMt)*
        pow2(Mst1) + (-5 + 6*lmMst1)*pow2(MuSUSY)) - 4*Dmst12*(-4 + 6*lmMst1 +
        3*lmMt)*pow3(Mst1))) - 2*Mst1*pow2(Mt)*(6*pow2(Dmglst1)*(35*pow2(
        Dmst12) + 3*Dmst12*(159 - 110*lmMst1 + 10*lmMt)*pow2(Mst2) + 200*(4 -
        3*lmMst1)*pow4(Mst2)) - 4*Dmglst1*Mst1*((11 + 60*lmMst1 - 60*lmMt)*
        pow2(Dmst12) + 25*Dmst12*(1 + 30*lmMst1 + 6*lmMt)*pow2(Mst2) + 100*(1 +
        12*lmMst1 + 6*lmMt)*pow4(Mst2)) - 25*pow2(Mst1)*((5 + 42*lmMst1 + 30*
        lmMt)*pow2(Dmst12) - 8*Dmst12*(4 + 3*lmMst1 + 6*lmMt)*pow2(Mst2) - 72*(
        1 + lmMst1 - lmMt + 2*lmMst1*lmMt + pow2(lmMst1) - 3*pow2(lmMt))*pow4(
        Mst2))))))) + 81000*oneLoopFlag*(12*Dmst12*Mt*MuSUSY*s2t*(Dmst12 - 2*
        pow2(Mst2))*pow2(Sbeta) - Tbeta*pow2(Dmst12)*pow2(s2t)*(pow2(MuSUSY)*(-
        1 + pow2(Sbeta)) + 12*pow2(Mst1)*pow2(Sbeta)) - 12*Tbeta*pow2(Mt)*pow2(
        Sbeta)*(pow2(Dmst12) - 2*Dmst12*pow2(Mst2) + 4*(-lmMst1 + lmMt)*pow4(
        Mst2)))*pow6(Mst1))) - 225*Mst1*threeLoopFlag*pow2(Al4p)*(-5*Mt*z3*(Mt*
        pow2(Mst2)*(-64*xDmglst1*pow3(Dmglst1)*(Dmst12*Mst1*(-(Dmst12*((9*(
        71189*MuSUSY + 1240*Mst1*Tbeta)*pow2(Mt) + 33*(4789*MuSUSY - 1738*Mst1*
        Tbeta)*pow2(Mst1)*pow2(s2t))*pow2(Sbeta) + Mt*s2t*(105358*Tbeta*pow2(
        MuSUSY)*(-1 + pow2(Sbeta)) + 3*Mst1*(38236*MuSUSY - 107331*Mst1*Tbeta)*
        pow2(Sbeta)))) + Mt*pow2(Mst2)*(9*Mt*(71189*MuSUSY + 1240*Mst1*Tbeta)*
        pow2(Sbeta) + 2*s2t*(52679*Tbeta*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 6*
        Mst1*(9559*MuSUSY + 26559*Mst1*Tbeta)*pow2(Sbeta)))) + 6*pow2(Mt)*(
        5511*Tbeta*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 2*Mst1*(26559*MuSUSY +
        5282*Mst1*Tbeta)*pow2(Sbeta))*pow4(Mst2)) - Mst1*Tbeta*pow2(MuSUSY)*(
        16*Dmst12*Mst1*Mt*s2t*(210810*pow2(Dmglst1)*(Dmst12 - pow2(Mst2)) -
        pow2(Mst1)*(17285*Dmst12 + 6024*pow2(Mst2)) - Dmglst1*Mst1*(105405*
        Dmst12 + 75164*pow2(Mst2))) - 8*pow2(Dmst12)*pow2(Mst1)*(-1260*Dmsqst1
        + 50485*Dmglst1*Mst1 + 149630*pow2(Dmglst1) + 2574*pow2(Mst1))*pow2(
        s2t) + 3*pow2(Mt)*(96*pow2(Dmglst1)*(10226*pow2(Dmst12) - 10226*Dmst12*
        pow2(Mst2) - 4715*pow4(Mst2)) + 48*Dmglst1*Mst1*(8213*pow2(Dmst12) -
        8213*Dmst12*pow2(Mst2) - 4664*pow4(Mst2)) - pow2(Mst1)*(31963*pow2(
        Dmst12) + 68000*Dmst12*pow2(Mst2) + 24064*pow4(Mst2)))) + Mst1*Tbeta*
        pow2(Sbeta)*(-16*Dmst12*Mst1*Mt*s2t*(pow2(Mst1)*(12*Dmst12*(22100*
        Dmglst1*Mst1 + 64701*pow2(Dmglst1) + 3290*pow2(Mst1)) + 96*(1963*
        Dmglst1*Mst1 + 7488*pow2(Dmglst1) + 106*pow2(Mst1))*pow2(Mst2)) + (-
        210810*pow2(Dmglst1)*(Dmst12 - pow2(Mst2)) + pow2(Mst1)*(17285*Dmst12 +
        6024*pow2(Mst2)) + Dmglst1*Mst1*(105405*Dmst12 + 75164*pow2(Mst2)))*
        pow2(MuSUSY)) - 8*pow2(Dmst12)*pow2(Mst1)*pow2(s2t)*(2*pow2(Dmglst1)*(
        37158*pow2(Mst1) + 74815*pow2(MuSUSY)) + Dmglst1*(50485*Mst1*pow2(
        MuSUSY) - 34776*pow3(Mst1)) - 18*(-143*pow2(Mst1)*pow2(MuSUSY) + 70*
        Dmsqst1*(4*pow2(Mst1) + pow2(MuSUSY)) + 440*pow4(Mst1))) + 3*pow2(Mt)*(
        -512*(2517*Dmglst1 + 2*(-107 + 12*lmMst1 - 12*lmMt)*Mst1)*pow3(Mst1)*
        pow4(Mst2) + 48*pow2(Mst1)*(pow2(Dmst12)*(-420*Dmsqst1 + 6970*Dmglst1*
        Mst1 + 45397*pow2(Dmglst1) - 614*pow2(Mst1)) - 8*Dmst12*(105*Dmsqst1 +
        3864*Dmglst1*Mst1 + 5666*pow2(Dmglst1) + 855*pow2(Mst1))*pow2(Mst2) -
        48400*pow2(Dmglst1)*pow4(Mst2)) + pow2(MuSUSY)*(96*pow2(Dmglst1)*(
        10226*pow2(Dmst12) - 10226*Dmst12*pow2(Mst2) - 4715*pow4(Mst2)) + 48*
        Dmglst1*Mst1*(8213*pow2(Dmst12) - 8213*Dmst12*pow2(Mst2) - 4664*pow4(
        Mst2)) - pow2(Mst1)*(31963*pow2(Dmst12) + 68000*Dmst12*pow2(Mst2) +
        24064*pow4(Mst2))))) + 16*Mst1*pow2(Sbeta)*(315*xDmsqst1*pow2(Dmsqst1)*
        (7*s2t*pow2(Dmst12)*(-2*Mt*MuSUSY + s2t*Tbeta*pow2(Mst1)) + 14*Dmst12*
        Mt*MuSUSY*s2t*pow2(Mst2) + 12*Tbeta*pow2(Mt)*(pow2(Dmst12) - Dmst12*
        pow2(Mst2) - 2*pow4(Mst2))) + Mst1*MuSUSY*(6*Dmst12*Mst1*Mt*s2t*(105*
        Dmsqst1*(7*Dmst12 + 8*pow2(Mst2)) + 3*pow2(Mst1)*(403*Dmst12 + 440*
        pow2(Mst2)) - 22*pow2(Dmglst1)*(2607*Dmst12 + 563*pow2(Mst2)) +
        Dmglst1*(-7642*Dmst12*Mst1 + 5796*Mst1*pow2(Mst2))) + 3*pow2(Dmst12)*
        pow2(Mst1)*(37582*Dmglst1*Mst1 + 105405*pow2(Dmglst1) + 3012*pow2(Mst1)
        )*pow2(s2t) + pow2(Mt)*(1404*pow2(Dmglst1)*(1065*pow2(Dmst12) - 1065*
        Dmst12*pow2(Mst2) - 512*pow4(Mst2)) - pow2(Mst1)*(51181*pow2(Dmst12) +
        49656*Dmst12*pow2(Mst2) + 10176*pow4(Mst2)) - 4*Dmglst1*Mst1*(86833*
        pow2(Dmst12) + 113412*Dmst12*pow2(Mst2) + 47112*pow4(Mst2)))))) - 2*
        Mst1*xDmst12*pow3(Dmst12)*(pow2(Sbeta)*(5040*Mt*xDmsqst1*pow2(Dmsqst1)*
        (-7*Mt*MuSUSY*s2t + 6*Tbeta*pow2(Mt) + 7*Tbeta*pow2(Mst1)*pow2(s2t)) +
        4*Mst1*MuSUSY*(Mst1*s2t*(27720*Dmsqst1 - 113856*Dmglst1*Mst1 - 1525128*
        pow2(Dmglst1) + 32783*pow2(Mst1))*pow2(Mt) - 3*Mt*pow2(Mst1)*(30241*
        Dmglst1*Mst1 - 421620*pow2(Dmglst1) + 11261*pow2(Mst1))*pow2(s2t) + 4*(
        -574156*Dmglst1*Mst1 + 747630*pow2(Dmglst1) - 76009*pow2(Mst1))*pow3(
        Mt) - (-1260*Dmsqst1 + 50485*Dmglst1*Mst1 + 149630*pow2(Dmglst1) +
        2574*pow2(Mst1))*pow3(Mst1)*pow3(s2t))) - 32*xDmglst1*pow3(Dmglst1)*(-
        2*s2t*pow2(Mt)*(52679*Tbeta*pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 237*Mst1*
        (-242*MuSUSY + 2031*Mst1*Tbeta)*pow2(Sbeta)) + pow2(Sbeta)*(66*Mt*(-
        4789*MuSUSY + 1738*Mst1*Tbeta)*pow2(Mst1)*pow2(s2t) - 9*(71189*MuSUSY +
        1240*Mst1*Tbeta)*pow3(Mt) + 52679*Tbeta*pow3(s2t)*pow4(Mst1))) + Tbeta*
        (8*pow2(Dmglst1)*(2*Mt*pow2(Mst1)*pow2(s2t)*(-74815*pow2(MuSUSY)*(-1 +
        pow2(Sbeta)) + 67452*pow2(Mst1)*pow2(Sbeta)) - 6*Mst1*s2t*pow2(Mt)*(-
        35135*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 378612*pow2(Mst1)*pow2(Sbeta))
        + 18*(10226*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 22733*pow2(Mst1)*pow2(
        Sbeta))*pow3(Mt) - 105405*pow2(Sbeta)*pow3(s2t)*pow5(Mst1)) + pow2(
        Mst1)*(8*Mst1*s2t*pow2(Mt)*(-40594*pow2(MuSUSY)*(-1 + pow2(Sbeta)) +
        11701*pow2(Mst1)*pow2(Sbeta)) - 3*(65963*pow2(MuSUSY)*(-1 + pow2(Sbeta)
        ) + 4*(10080*Dmsqst1 + 11311*pow2(Mst1))*pow2(Sbeta))*pow3(Mt) + Mt*
        pow2(s2t)*(37669*pow2(Mst1)*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 2520*
        Dmsqst1*(4*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + pow2(Mst1)*pow2(Sbeta)) +
        2664*pow2(Sbeta)*pow4(Mst1)) - 24096*pow2(Sbeta)*pow3(s2t)*pow5(Mst1))
        - 8*Dmglst1*Mst1*(Mt*pow2(Mst1)*pow2(s2t)*(50485*pow2(MuSUSY)*(-1 +
        pow2(Sbeta)) - 40314*pow2(Mst1)*pow2(Sbeta)) + pow2(Mt)*(285974*Mst1*
        s2t*pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 82132*s2t*pow2(Sbeta)*pow3(Mst1))
        - 9*(8213*pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 16972*pow2(Mst1)*pow2(
        Sbeta))*pow3(Mt) + 37582*pow2(Sbeta)*pow3(s2t)*pow5(Mst1))))) + 384*z2*
        pow3(Mst1)*(pow2(Mst2)*(15*MuSUSY*pow2(Sbeta)*(10*Mt*s2t*(Dmsqst1 +
        pow2(Mst1))*(-8*Dmst12*shiftst2*pow2(Mst2)*pow2(Mt) + (-shiftst1 +
        shiftst2)*pow2(Dmst12)*pow2(Mst1)*pow2(s2t) + 8*(shiftst1 - shiftst2)*
        pow2(Mt)*pow4(Mst2)) + shiftst3*pow2(Mst1)*(-(Mt*pow2(Dmst12)*pow2(
        Mst1)*pow3(s2t)) + 8*s2t*pow3(Mt)*(pow2(Dmst12) - Dmst12*pow2(Mst2) +
        pow4(Mst2)))) + Tbeta*(50*(Dmsqst1 + pow2(Mst1))*pow2(Mt)*(pow2(Dmst12)
        *pow2(s2t)*(-((shiftst1 + shiftst2)*pow2(MuSUSY)*(-1 + pow2(Sbeta))) +
        6*(-3*shiftst1 + shiftst2)*pow2(Mst1)*pow2(Sbeta)) + 2*Dmst12*pow2(
        Mst2)*((shiftst1 - shiftst2)*pow2(MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta))
        + 12*(shiftst2*pow2(Mt) + (shiftst1 - shiftst2)*pow2(Mst1)*pow2(s2t))*
        pow2(Sbeta)) + 24*(shiftst1 + shiftst2)*pow2(Mt)*pow2(Sbeta)*pow4(Mst2)
        ) + 5*shiftst3*pow2(Mst1)*(-(Dmst12*pow2(Mt)*pow2(s2t)*((3*Dmst12 - 2*
        pow2(Mst2))*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 6*pow2(Mst1)*(7*Dmst12 -
        4*pow2(Mst2))*pow2(Sbeta))) + 24*pow2(Sbeta)*(pow2(Dmst12) - Dmst12*
        pow2(Mst2) + pow4(Mst2))*pow4(Mt)))) - xDmst12*pow3(Dmst12)*(2*pow2(
        Mst1)*pow2(Mt)*(-((15*shiftst1 + 10*shiftst2 + 9*shiftst3)*Tbeta*pow2(
        MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta))) + 60*Mt*shiftst3*(MuSUSY*s2t +
        Mt*Tbeta)*pow2(Sbeta)) + 50*Dmsqst1*(-4*shiftst2*pow2(Mt)*(Tbeta*pow2(
        MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta)) + 6*Mt*(MuSUSY*s2t - Mt*Tbeta)*
        pow2(Sbeta)) + pow2(Sbeta)*(2*Mt*(MuSUSY*s2t*(-2*shiftst1 + shiftst2) +
        2*Mt*(-4*shiftst1 + shiftst2)*Tbeta)*pow2(Mst1)*pow2(s2t) + (shiftst1 -
        shiftst2)*Tbeta*pow4(Mst1)*pow4(s2t))) + pow2(Sbeta)*(-5*Mt*(MuSUSY*
        s2t*(40*shiftst1 - 20*shiftst2 + 7*shiftst3) + 2*Mt*(80*shiftst1 - 20*
        shiftst2 + 29*shiftst3)*Tbeta)*pow2(s2t)*pow4(Mst1) + 5*(10*shiftst1 -
        10*shiftst2 + shiftst3)*Tbeta*pow4(s2t)*pow6(Mst1))))))/(3.888e6*Tbeta*
        pow2(Sbeta)*pow6(Mst1)*pow6(Mst2));
}

/**
 * @return The off-diagonal (1, 2) = (2, 1) matrix element of the Higgs mass matrix as a double for the hierarchy 'H32q2g'
 */
double H32q2g::getS12() const {
   return -(128*Al4p*xDmglst1*pow3(Dmglst1)*(45*MuSUSY*twoLoopFlag*pow2(Mst1)*pow2(
        Mt)*(80*Dmst12*Mt*MuSUSY*s2t*(pow2(Dmst12) - Dmst12*pow2(Mst2) + pow4(
        Mst2)) - Tbeta*(pow2(Dmst12)*pow2(Mst2)*(660*Mst1*Mt*s2t + (-181 + 90*
        lmMst1 - 90*lmMt)*pow2(Mt) - 60*pow2(Mst1)*pow2(s2t)) + (-660*Mst1*Mt*
        s2t + (181 - 90*lmMst1 + 90*lmMt)*pow2(Mt) + 120*pow2(Mst1)*pow2(s2t))*
        pow3(Dmst12) + Dmst12*Mt*((181 - 90*lmMst1 + 90*lmMt)*Mt - 660*Mst1*
        s2t)*pow4(Mst2) + 12*(16 - 15*lmMst1 + 15*lmMt)*pow2(Mt)*pow6(Mst2))) +
        Al4p*threeLoopFlag*(10*MuSUSY*((2219399 + 9600*lmMst1 + 960*lmMt)*Mst1*
        MuSUSY + 12*Dmsqst1*(-1541 + 420*lmMst1 - 420*lmMt)*Tbeta)*pow4(Mt)*
        pow6(Mst2) + (MuSUSY*pow2(Mst1)*pow2(Mt)*(735*Dmst12*Mst1*s2t*Tbeta*(
        Dmst12*Mst1*s2t*(-14217821 - 161940*lmMst1 + 28800*pow2(lmMst1))*(2*
        Dmst12 - pow2(Mst2)) + 20*Mt*(-520781 - 17172*lmMst1 + 192*lmMt +
        11520*pow2(lmMst1))*(pow2(Dmst12) - Dmst12*pow2(Mst2) + pow4(Mst2))) -
        2*Mt*(-5*Dmst12*(98*MuSUSY*s2t*(14217821 + 161940*lmMst1 - 28800*pow2(
        lmMst1)) - 3*Mt*Tbeta*(1417174939 - 401268*lmMt + 4*lmMst1*(4061413 +
        37800*lmMt) - 3185280*pow2(lmMst1) + 211680*pow2(lmMt)))*(pow2(Dmst12)
        - Dmst12*pow2(Mst2) + pow4(Mst2)) + 196*Mt*Tbeta*(54198467 + 43950*lmMt
        + 180*lmMst1*(6353 + 135*lmMt) - 272700*pow2(lmMst1) + 32400*pow2(lmMt)
        )*pow6(Mst2))))/196.)) + MuSUSY*(225*threeLoopFlag*pow2(Al4p)*(-192*z2*
        pow4(Mst1)*(2*Dmst12*MuSUSY*pow2(Mt)*pow2(s2t)*(2*pow2(Dmst12)*(100*
        Dmsqst1*shiftst2 + (15*shiftst1 + 10*shiftst2 + 9*shiftst3)*pow2(Mst1))
        - 5*Dmst12*(10*Dmsqst1*(shiftst1 + shiftst2) + (10*(shiftst1 +
        shiftst2) + 3*shiftst3)*pow2(Mst1))*pow2(Mst2) + 10*(10*Dmsqst1*(
        shiftst1 - shiftst2) + (10*shiftst1 - 10*shiftst2 + shiftst3)*pow2(
        Mst1))*pow4(Mst2)) - 5*Tbeta*(Mt*shiftst3*pow2(Dmst12)*(7*Dmst12 - 3*
        pow2(Mst2))*pow3(s2t)*pow4(Mst1) + 10*shiftst2*(-(Mt*pow2(Dmst12)*pow2(
        Mst1)*(Dmsqst1 + pow2(Mst1))*(2*Dmst12 - 3*pow2(Mst2))*pow3(s2t)) + 24*
        s2t*pow3(Mt)*(Dmsqst1*pow3(Dmst12) - (Dmsqst1 + pow2(Mst1))*(Dmst12 +
        pow2(Mst2))*pow4(Mst2))) - 24*s2t*shiftst3*pow2(Mst1)*pow3(Mt)*(-(pow2(
        Dmst12)*pow2(Mst2)) + pow3(Dmst12) + Dmst12*pow4(Mst2) - pow6(Mst2)) +
        10*shiftst1*(Dmsqst1 + pow2(Mst1))*(Mt*pow2(Dmst12)*pow2(Mst1)*(4*
        Dmst12 - 3*pow2(Mst2))*pow3(s2t) + 24*s2t*pow3(Mt)*pow6(Mst2)))) + Mt*(
        (800*xDmsqst1*pow2(Dmsqst1)*(3*Dmst12*s2t*pow2(Mst1)*(pow2(Dmst12)*(12*
        s2t*(shiftst1 - shiftst2)*(2*Mt*MuSUSY - 3*s2t*Tbeta*pow2(Mst1)) + (205
        + 144*shiftst2)*Tbeta*pow2(Mt)) - Dmst12*pow2(Mst2)*(6*s2t*(shiftst1 -
        shiftst2)*(4*Mt*MuSUSY - 3*s2t*Tbeta*pow2(Mst1)) + (205 + 144*shiftst2)
        *Tbeta*pow2(Mt)) + Mt*(24*MuSUSY*s2t*(shiftst1 - shiftst2) + Mt*(205 +
        144*shiftst2)*Tbeta)*pow4(Mst2)) + 16*(Dmglst1*(8 + 3*lmMst1 - 3*lmMt)*
        Mt + Mst1*((29 - 6*lmMst1 + 6*lmMt)*Mt + 27*Mst1*s2t*(-shiftst1 +
        shiftst2)))*Tbeta*pow2(Mt)*pow6(Mst2)))/3. - 5*Mst1*z3*(35280*Dmst12*
        Mst1*s2t*Tbeta*xDmsqst1*pow2(Dmsqst1)*pow2(Mt)*(pow2(Dmst12) - Dmst12*
        pow2(Mst2) + pow4(Mst2)) + Mst1*Mt*MuSUSY*(-2*pow2(Dmst12)*pow2(Mst1)*(
        20*(-252*Dmsqst1 + 10097*Dmglst1*Mst1 + 29926*pow2(Dmglst1))*(2*Dmst12
        - pow2(Mst2)) - pow2(Mst1)*(37669*Dmst12 + 10296*pow2(Mst2)))*pow2(s2t)
        + 16*Dmst12*Mst1*Mt*s2t*(210810*pow2(Dmglst1)*(pow2(Dmst12) - Dmst12*
        pow2(Mst2) + pow4(Mst2)) + pow2(Mst1)*(-40594*pow2(Dmst12) + 17285*
        Dmst12*pow2(Mst2) + 6024*pow4(Mst2)) + Dmglst1*Mst1*(-285974*pow2(
        Dmst12) + 105405*Dmst12*pow2(Mst2) + 75164*pow4(Mst2))) + 3*pow2(Mt)*(-
        (pow2(Mst1)*(-31963*pow2(Dmst12)*pow2(Mst2) + 131926*pow3(Dmst12) -
        68000*Dmst12*pow4(Mst2) - 24064*pow6(Mst2))) + 48*Dmglst1*Mst1*(-8213*
        pow2(Dmst12)*pow2(Mst2) + 8213*Dmst12*(pow2(Dmst12) + pow4(Mst2)) +
        4664*pow6(Mst2)) + 96*pow2(Dmglst1)*(-10226*pow2(Dmst12)*pow2(Mst2) +
        10226*Dmst12*(pow2(Dmst12) + pow4(Mst2)) + 4715*pow6(Mst2)))) + 32*Mt*
        xDmglst1*pow3(Dmglst1)*(44*Mt*MuSUSY*(4789*Dmst12*Mst1*s2t*(pow2(
        Dmst12) - Dmst12*pow2(Mst2) + pow4(Mst2)) + 1503*Mt*pow6(Mst2)) - 3*
        Mst1*Tbeta*(52679*pow2(Dmst12)*pow2(Mst1)*(2*Dmst12 - pow2(Mst2))*pow2(
        s2t) + 38236*Dmst12*Mst1*Mt*s2t*(pow2(Dmst12) - Dmst12*pow2(Mst2) +
        pow4(Mst2)) + 3*pow2(Mt)*(-71189*pow2(Dmst12)*pow2(Mst2) + 71189*
        Dmst12*(pow2(Dmst12) + pow4(Mst2)) + 35412*pow6(Mst2)))) - 4*Tbeta*
        pow2(Mst1)*(3*Mt*pow2(Dmst12)*pow2(Mst1)*(210810*pow2(Dmglst1)*(2*
        Dmst12 - pow2(Mst2)) - pow2(Mst1)*(11261*Dmst12 + 6024*pow2(Mst2)) -
        Dmglst1*Mst1*(30241*Dmst12 + 75164*pow2(Mst2)))*pow2(s2t) - (-1260*
        Dmsqst1 + 50485*Dmglst1*Mst1 + 149630*pow2(Dmglst1) + 2574*pow2(Mst1))*
        pow3(Dmst12)*pow3(Mst1)*pow3(s2t) - Dmst12*Mst1*s2t*pow2(Mt)*(pow2(
        Dmst12)*(-27720*Dmsqst1 + 113856*Dmglst1*Mst1 + 1525128*pow2(Dmglst1) -
        32783*pow2(Mst1)) - 12*Dmst12*(-735*Dmsqst1 + 7642*Dmglst1*Mst1 +
        57354*pow2(Dmglst1) - 1209*pow2(Mst1))*pow2(Mst2) - 24*(-2898*Dmglst1*
        Mst1 + 6193*pow2(Dmglst1) - 60*(7*Dmsqst1 + 11*pow2(Mst1)))*pow4(Mst2))
        + 2*pow3(Mt)*(1404*pow2(Dmglst1)*(-1065*pow2(Dmst12)*pow2(Mst2) + 1065*
        pow3(Dmst12) + 1065*Dmst12*pow4(Mst2) + 512*pow6(Mst2)) + pow2(Mst1)*(
        51181*pow2(Dmst12)*pow2(Mst2) - 152018*pow3(Dmst12) + 49656*Dmst12*
        pow4(Mst2) + 10176*pow6(Mst2)) + 4*Dmglst1*Mst1*(86833*pow2(Dmst12)*
        pow2(Mst2) - 287078*pow3(Dmst12) + 113412*Dmst12*pow4(Mst2) + 47112*
        pow6(Mst2))))))) + Mt*(40500*Dmst12*oneLoopFlag*s2t*(-2*Dmst12*Mt*(
        MuSUSY*s2t + 6*Mt*Tbeta)*pow2(Mst2) + pow2(Dmst12)*(2*Mt*MuSUSY*s2t +
        8*Tbeta*pow2(Mt) - Tbeta*pow2(Mst1)*pow2(s2t)) + 24*Tbeta*pow2(Mt)*
        pow4(Mst2))*pow6(Mst1) + 1440*Al4p*twoLoopFlag*pow3(Mst1)*(-25*s2t*
        pow3(Dmst12)*(8*(5 + 6*lmMst1)*MuSUSY*pow2(Mst1)*pow2(Mt) + 2*pow2(
        Dmglst1)*(24*MuSUSY*pow2(Mt) + (11 - 6*lmMst1)*Tbeta*pow2(s2t)*pow3(
        Mst1)) - 8*Dmglst1*((5 - 6*lmMst1)*Mst1*MuSUSY*pow2(Mt) + (-1 + 3*
        lmMst1)*Tbeta*pow2(s2t)*pow4(Mst1)) + 3*Tbeta*pow2(s2t)*pow5(Mst1)) +
        Mt*(50*Dmst12*MuSUSY*s2t*(Mst1*s2t*pow2(Dmst12)*(16*Dmglst1*(1 - 3*
        lmMst1)*Mst1 + (44 - 24*lmMst1)*pow2(Dmglst1) + (1 + 6*lmMst1)*pow2(
        Mst1)) + 4*Dmst12*Mt*(6*pow2(Dmglst1) + pow2(Mst1))*pow2(Mst2) -
        Dmst12*Mst1*s2t*(8*Dmglst1*(1 - 3*lmMst1)*Mst1 + (22 - 12*lmMst1)*pow2(
        Dmglst1) + 3*pow2(Mst1))*pow2(Mst2) - 4*Mt*(Dmglst1*(5 - 6*lmMst1)*Mst1
        + 6*pow2(Dmglst1) - 3*(1 + 2*lmMst1)*pow2(Mst1))*pow4(Mst2)) + 2*Tbeta*
        (-25*pow2(Mst1)*(-(pow2(Dmst12)*pow2(Mst2)*(24*(1 - 3*lmMst1)*Mst1*Mt*
        s2t + 2*(5 + 6*lmMst1 + 6*lmMt)*pow2(Mt) + 9*(1 + 2*lmMst1)*pow2(Mst1)*
        pow2(s2t))) + (-48*(-1 + lmMst1)*Mst1*Mt*s2t - 4*(5 + 6*lmMst1)*pow2(
        Mt) + 6*(1 + 3*lmMst1)*pow2(Mst1)*pow2(s2t))*pow3(Dmst12) + 4*Dmst12*
        Mt*(2*(5 + 6*lmMst1 + 3*lmMt)*Mt - 9*(1 + 2*lmMst1)*Mst1*s2t)*pow4(
        Mst2) + 72*(1 + lmMst1 + lmMt)*pow2(Mt)*pow6(Mst2)) + 2*pow2(Dmglst1)*(
        9*pow2(Dmst12)*pow2(Mst2)*(5*Mst1*Mt*s2t + (-109 + 10*lmMst1 - 10*lmMt)
        *pow2(Mt) - 25*pow2(Mst1)*pow2(s2t)) + (60*(-29 + 15*lmMst1)*Mst1*Mt*
        s2t + (981 - 90*lmMst1 + 90*lmMt)*pow2(Mt) + 450*pow2(Mst1)*pow2(s2t))*
        pow3(Dmst12) + 3*Dmst12*Mt*((327 - 30*lmMst1 + 30*lmMt)*Mt + 50*(11 -
        6*lmMst1)*Mst1*s2t)*pow4(Mst2) + 100*(17 - 3*lmMst1 + 3*lmMt)*pow2(Mt)*
        pow6(Mst2)) + Dmglst1*Mst1*(-3*pow2(Dmst12)*(pow2(Mst2)*(100*Mst1*Mt*
        s2t - 2*(-41 + 90*lmMst1 + 10*lmMt)*pow2(Mt) + 25*(5 - 6*lmMst1)*pow2(
        Mst1)*pow2(s2t)) + Dmst12*(-1200*lmMst1*Mst1*Mt*s2t + 8*(17 - 30*lmMst1
        + 5*lmMt)*pow2(Mt) + 25*(-5 + 6*lmMst1)*pow2(Mst1)*pow2(s2t))) + 300*
        Dmst12*Mt*((3 - 6*lmMst1)*Mt + 2*(1 - 6*lmMst1)*Mst1*s2t)*pow4(Mst2) -
        200*(-7 + 15*lmMst1 + 3*lmMt)*pow2(Mt)*pow6(Mst2))))) + 2*Mst1*
        threeLoopFlag*pow2(Al4p)*(5*Mst1*Mt*MuSUSY*(40*Dmglst1*(2*pow3(Dmst12)*
        (Mt*s2t*(6000*Dmsqst1 + 6*(-1282471 + 7264*lmMst1 + 18120*pow2(lmMst1))
        *pow2(Mst1)) + 5*Mst1*(403559 + 384*(lmMst1 + lmMt) - 4608*pow2(lmMst1)
        )*pow2(Mt) + 15*(-84209 - 1264*lmMst1 + 240*pow2(lmMst1))*pow2(s2t)*
        pow3(Mst1)) - pow2(Dmst12)*pow2(Mst2)*(Mt*s2t*(12000*Dmsqst1 + 6*(-
        949861 + 1944*lmMst1 + 11520*pow2(lmMst1))*pow2(Mst1)) + 10*Mst1*(
        403559 + 384*(lmMst1 + lmMt) - 4608*pow2(lmMst1))*pow2(Mt) + 15*(-84209
        - 1264*lmMst1 + 240*pow2(lmMst1))*pow2(s2t)*pow3(Mst1)) + 10*Dmst12*Mt*
        (Mst1*Mt*(403559 + 384*(lmMst1 + lmMt) - 4608*pow2(lmMst1)) + s2t*(
        1200*Dmsqst1 - 12*(-33261 + 532*lmMst1 + 660*pow2(lmMst1))*pow2(Mst1)))
        *pow4(Mst2) + 240*Mst1*(9631 + 16*lmMst1 + 48*lmMt - 192*pow2(lmMst1))*
        pow2(Mt)*pow6(Mst2)) + 15*Mst1*(-(pow2(Dmst12)*pow2(Mst2)*(240*Mt*s2t*(
        560*Dmsqst1 + (-10473 + 40*lmMst1 + 256*pow2(lmMst1))*pow2(Mst1)) -
        Mst1*(852541 + 9216*lmMst1 - 10240*lmMt + 6144*pow2(lmMst1))*pow2(Mt) +
        2400*Dmsqst1*Mst1*(7 - 12*lmMst1*(-2 + shiftst1 + shiftst2) + 18*(
        shiftst1 + shiftst2))*pow2(s2t) - 80*(1429 - 454*lmMst1 - 180*(shiftst1
        + shiftst2) + 360*lmMst1*(shiftst1 + shiftst2) - 126*shiftst3 + 108*
        lmMst1*shiftst3 + 24*pow2(lmMst1))*pow2(s2t)*pow3(Mst1))) - 2*pow3(
        Dmst12)*(80*Mt*s2t*(-840*Dmsqst1 + (36863 - 80*lmMst1 - 552*pow2(
        lmMst1))*pow2(Mst1)) + Mst1*(1763661 - 47104*lmMst1 + 5120*lmMt -
        24576*pow2(lmMst1))*pow2(Mt) + pow2(s2t)*(2400*Dmsqst1*Mst1*(-7 + 24*
        lmMst1*(-1 + shiftst2) - 36*shiftst2) + (-350605 - 4320*shiftst1 -
        2880*shiftst2 - 8352*shiftst3 + 96*lmMst1*(-115 + 90*shiftst1 + 60*
        shiftst2 + 54*shiftst3) + 2160*pow2(lmMst1))*pow3(Mst1))) + 160*Dmst12*
        (Mt*s2t*(840*Dmsqst1 + 4*(1361 + 10*lmMst1 + 54*pow2(lmMst1))*pow2(
        Mst1)) + Mst1*(11389 - 704*lmMst1 + 192*lmMt - 384*pow2(lmMst1))*pow2(
        Mt) + pow2(s2t)*(180*Dmsqst1*(3 - 2*lmMst1)*Mst1*(shiftst1 - shiftst2)
        + 18*(1 - 2*lmMst1)*(10*shiftst1 - 10*shiftst2 + shiftst3)*pow3(Mst1)))
        *pow4(Mst2) + 1920*Mst1*(349 - 56*lmMst1 + 24*lmMt - 32*pow2(lmMst1))*
        pow2(Mt)*pow6(Mst2)) + 16*pow2(Dmglst1)*(pow2(Dmst12)*pow2(Mst2)*(30*
        Mst1*Mt*s2t*(-954181 - 11256*lmMst1 + 11520*pow2(lmMst1)) + 4*(-6233611
        - 58800*lmMst1 + 4560*lmMt + 14400*pow2(lmMst1))*pow2(Mt) + 15*(655743
        + 5288*lmMst1 - 11820*pow2(lmMst1))*pow2(Mst1)*pow2(s2t)) - 2*(15*Mst1*
        Mt*s2t*(-954181 - 11256*lmMst1 + 11520*pow2(lmMst1)) + 2*(-6233611 -
        58800*lmMst1 + 4560*lmMt + 14400*pow2(lmMst1))*pow2(Mt) + 15*(655743 +
        5288*lmMst1 - 11820*pow2(lmMst1))*pow2(Mst1)*pow2(s2t))*pow3(Dmst12) +
        2*Dmst12*Mt*(15*Mst1*s2t*(954181 + 11256*lmMst1 - 11520*pow2(lmMst1)) -
        2*Mt*(-6233611 - 58800*lmMst1 + 4560*lmMt + 14400*pow2(lmMst1)))*pow4(
        Mst2) - 50*(-230009 - 4704*lmMst1 + 480*lmMt + 1152*pow2(lmMst1))*pow2(
        Mt)*pow6(Mst2))) + 4*Tbeta*(2*pow2(Dmglst1)*(3*Dmst12*Mst1*(Dmst12*Mt*
        pow2(Mst2)*(150*Mt*s2t*(360*Dmsqst1 + (-515917 - 6972*lmMst1 + 192*lmMt
        + 11520*pow2(lmMst1))*pow2(Mst1)) + 4*Mst1*(84334067 + 120*lmMst1*(2843
        - 120*lmMt) + 202200*lmMt - 828000*pow2(lmMst1) - 21600*pow2(lmMt))*
        pow2(Mt) + 75*(954181 + 11256*lmMst1 - 11520*pow2(lmMst1))*pow2(s2t)*
        pow3(Mst1)) - pow2(Dmst12)*(600*s2t*(90*Dmsqst1 + (-274009 - 964*lmMst1
        - 104*lmMt + 8310*pow2(lmMst1))*pow2(Mst1))*pow2(Mt) + 150*Mt*(954181 +
        11256*lmMst1 - 11520*pow2(lmMst1))*pow2(s2t)*pow3(Mst1) + 4*Mst1*(
        84334067 + 120*lmMst1*(2843 - 120*lmMt) + 202200*lmMt - 828000*pow2(
        lmMst1) - 21600*pow2(lmMt))*pow3(Mt) + 25*(-655743 - 5288*lmMst1 +
        11820*pow2(lmMst1))*pow3(s2t)*pow4(Mst1)) - 4*(Mst1*Mt*(84334067 + 120*
        lmMst1*(2843 - 120*lmMt) + 202200*lmMt - 828000*pow2(lmMst1) - 21600*
        pow2(lmMt)) + 75*s2t*(180*Dmsqst1 + (32101 - 5044*lmMst1 + 400*lmMt -
        5100*pow2(lmMst1))*pow2(Mst1)))*pow2(Mt)*pow4(Mst2)) + 160*(90*Dmsqst1*
        (206 - 15*lmMst1 + 15*lmMt) + (-3030652 - 19305*lmMt + 6*lmMst1*(1183 +
        645*lmMt) + 55530*pow2(lmMst1) + 5400*pow2(lmMt))*pow2(Mst1))*pow3(Mt)*
        pow6(Mst2)) - Dmglst1*(1200*Dmsqst1*Mst1*Mt*(-(pow2(Dmst12)*pow2(Mst2)*
        (14100*Mst1*Mt*s2t + 4*(557 + 120*lmMst1 - 120*lmMt)*pow2(Mt) + 375*
        pow2(Mst1)*pow2(s2t))) + 2*(7050*Mst1*Mt*s2t + 2*(557 + 120*lmMst1 -
        120*lmMt)*pow2(Mt) + 375*pow2(Mst1)*pow2(s2t))*pow3(Dmst12) + 4*Dmst12*
        Mt*((557 + 120*lmMst1 - 120*lmMt)*Mt + 3525*Mst1*s2t)*pow4(Mst2) + 200*
        (26 + 3*lmMst1 - 3*lmMt)*pow2(Mt)*pow6(Mst2)) + pow3(Mst1)*(12*Mt*pow2(
        Dmst12)*pow2(Mst2)*(300*Mst1*Mt*s2t*(17539 + 574*lmMst1 + 160*lmMt -
        1920*pow2(lmMst1)) + (39474953 + 3300*lmMt + 20*lmMst1*(-2639 + 4380*
        lmMt) - 319200*pow2(lmMst1) - 14400*pow2(lmMt))*pow2(Mt) + 375*(-33261
        + 532*lmMst1 + 660*pow2(lmMst1))*pow2(Mst1)*pow2(s2t)) - pow3(Dmst12)*(
        3600*Mst1*s2t*(12383 + 4128*lmMst1 + 80*lmMt - 1260*pow2(lmMst1))*pow2(
        Mt) + 225*Mt*(284641 + 8696*lmMst1 + 1680*pow2(lmMst1))*pow2(Mst1)*
        pow2(s2t) + 8*(193364399 + 90000*lmMt - 300*lmMst1*(3781 + 582*lmMt) -
        2005200*pow2(lmMst1) - 43200*pow2(lmMt))*pow3(Mt) + 375*(84209 + 1264*
        lmMst1 - 240*pow2(lmMst1))*pow3(Mst1)*pow3(s2t)) + 80*Dmst12*(225*Mst1*
        s2t*(-4539 + 596*lmMst1 - 48*lmMt + 516*pow2(lmMst1)) - 2*Mt*(-3746977
        - 4005*lmMt + 18*lmMst1*(2711 + 1215*lmMt) + 52380*pow2(lmMst1)))*pow2(
        Mt)*pow4(Mst2) + 24000*(9961 + 66*lmMt - 10*lmMst1*(67 + 24*lmMt) + 42*
        pow2(lmMst1) + 72*pow2(lmMt))*pow3(Mt)*pow6(Mst2))) - 15*(pow4(Mst1)*(
        2*pow2(Dmst12)*pow2(Mst2)*(150*Mst1*s2t*(-2071 + 96*lmMt + 288*shiftst3
        - 8*lmMst1*(37 + 18*shiftst3) + 600*pow2(lmMst1))*pow2(Mt) - 300*Mt*(
        1361 + 10*lmMst1 + 54*pow2(lmMst1))*pow2(Mst1)*pow2(s2t) - (-2284899 +
        49840*lmMt - 32*lmMst1*(-1793 + 555*lmMt) + 87360*pow2(lmMst1) + 28800*
        pow2(lmMt))*pow3(Mt) + 1350*(-1 + 2*lmMst1)*(10*shiftst1 - 10*shiftst2
        + shiftst3)*pow3(Mst1)*pow3(s2t)) - pow3(Dmst12)*(15*Mst1*s2t*(-102747
        + 640*lmMt + 6720*shiftst3 - 32*lmMst1*(331 + 90*shiftst3) + 13888*
        pow2(lmMst1))*pow2(Mt) - 75*Mt*(-20531 + 200*lmMst1 + 1200*pow2(lmMst1)
        )*pow2(Mst1)*pow2(s2t) - 4*(-3454599 + 16840*lmMt + 48*lmMst1*(262 +
        405*lmMt) + 46560*pow2(lmMst1))*pow3(Mt) + 50*(1429 - 720*shiftst1 +
        360*shiftst2 - 234*shiftst3 + 2*lmMst1*(-227 + 720*shiftst1 - 360*
        shiftst2 + 126*shiftst3) + 24*pow2(lmMst1))*pow3(Mst1)*pow3(s2t)) +
        400*(-6*Dmst12*Mst1*s2t*(631 + 36*lmMt + 90*shiftst2 + 27*shiftst3 -
        lmMst1*(1 + 180*shiftst2 + 18*shiftst3) + 21*pow2(lmMst1)) + Dmst12*Mt*
        (11697 + 448*lmMst1 + 330*lmMt - 372*lmMst1*lmMt + 408*pow2(lmMst1) +
        288*pow2(lmMt)) + 54*(1 - 2*lmMst1)*Mst1*s2t*(10*shiftst1 - 10*shiftst2
        + shiftst3)*pow2(Mst2) + 8*Mt*(434 - 83*lmMst1 + 174*lmMt - 66*lmMst1*
        lmMt + 183*pow2(lmMst1) + 108*pow2(lmMt))*pow2(Mst2))*pow2(Mt)*pow4(
        Mst2)) + 500*Dmsqst1*pow2(Mst1)*(-(pow2(Dmst12)*pow2(Mst2)*(483*Mst1*
        s2t*pow2(Mt) + 252*Mt*pow2(Mst1)*pow2(s2t) + 16*(65 - 6*lmMst1 + 6*
        lmMt)*pow3(Mt) + 54*(3 - 2*lmMst1)*(shiftst1 - shiftst2)*pow3(Mst1)*
        pow3(s2t))) + pow3(Dmst12)*(6*Mst1*s2t*(173 - 144*lmMst1*(-1 +
        shiftst2) + 216*shiftst2)*pow2(Mt) + 504*Mt*pow2(Mst1)*pow2(s2t) + 16*(
        65 - 6*lmMst1 + 6*lmMt)*pow3(Mt) + 3*(7 + 72*shiftst1 - 36*shiftst2 +
        24*lmMst1*(1 - 2*shiftst1 + shiftst2))*pow3(Mst1)*pow3(s2t)) + pow2(Mt)
        *(8*Dmst12*(2*(65 - 6*lmMst1 + 6*lmMt)*Mt - 9*Mst1*s2t*(1 - 12*lmMst1*(
        -1 + shiftst2) + 18*shiftst2))*pow4(Mst2) + 48*((40 - 6*lmMst1 + 6*
        lmMt)*Mt + 9*(3 - 2*lmMst1)*Mst1*s2t*(shiftst1 - shiftst2))*pow6(Mst2))
        )))))))/(3.888e6*Tbeta*pow6(Mst1)*pow6(Mst2));
}

/**
 * @return returns the susy log^0 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double H32q2g::calc_coef_at_as2_no_sm_logs_log0() const {

   const double result =
      (-8*Mst1*pow3(Dmglst1)*(Mst1*Mt*pow2(Dmst12)*pow2(Mst2)*(-31305120*
        Dmsqst1*Mt*s2t + Mt*s2t*(21284082326 - 17749864125*z3)*pow2(Mst1) + 40*
        Mst1*(-16826654 + 15379875*z3)*pow2(Mt) - 7350*(-520877 + 430155*z3)*
        pow2(s2t)*pow3(Mst1)) + Mst1*pow3(Dmst12)*(31305120*Dmsqst1*s2t*pow2(
        Mt) + 2*s2t*(-31902674758 + 26534253375*z3)*pow2(Mst1)*pow2(Mt) +
        14700*Mt*(-520877 + 430155*z3)*pow2(s2t)*pow3(Mst1) + 40*Mst1*(16826654
        - 15379875*z3)*pow3(Mt) + 245*(-14217821 + 11852775*z3)*pow3(s2t)*pow4(
        Mst1)) - 4*Dmst12*Mst1*(-7826280*Dmsqst1*s2t + 10*Mst1*Mt*(-16826654 +
        15379875*z3) + 49*s2t*(-108352984 + 89636625*z3)*pow2(Mst1))*pow2(Mt)*
        pow4(Mst2) + 392*(154740*Dmsqst1 + (10583177 - 8913375*z3)*pow2(Mst1))*
        pow3(Mt)*pow6(Mst2)) + 49*(3750*Mt*pow2(Dmsqst1)*pow2(Mst1)*(pow2(
        Dmst12)*pow2(Mst2)*(-1664*Mst1*Mt*s2t - 4*(-566 + 567*z3)*pow2(Mt) + 3*
        (410 - 441*z3)*pow2(Mst1)*pow2(s2t)) + 2*(832*Mst1*Mt*s2t + 2*(-566 +
        567*z3)*pow2(Mt) + 3*(-410 + 441*z3)*pow2(Mst1)*pow2(s2t))*pow3(Dmst12)
        + 4*Dmst12*Mt*(416*Mst1*s2t + Mt*(-566 + 567*z3))*pow4(Mst2) + 8*(-598
        + 567*z3)*pow2(Mt)*pow6(Mst2)) + 150*Dmsqst1*pow4(Mst1)*(4*Mt*pow2(
        Dmst12)*pow2(Mst2)*(19600*Mst1*Mt*s2t + (-9122 + 14175*z3)*pow2(Mt) +
        450*(2 - 21*z3)*pow2(Mst1)*pow2(s2t)) + pow3(Dmst12)*(20800*Mst1*s2t*
        pow2(Mt) + 75*Mt*(274 + 63*z3)*pow2(Mst1)*pow2(s2t) - 16*(-6536 +
        14175*z3)*pow3(Mt) + 8400*pow3(Mst1)*pow3(s2t)) + 200*Dmst12*(-888*
        Mst1*s2t + Mt*(-158 + 567*z3))*pow2(Mt)*pow4(Mst2) + 201600*pow3(Mt)*
        pow6(Mst2)) - pow6(Mst1)*(-20*Mt*pow2(Dmst12)*pow2(Mst2)*(300*Mst1*Mt*
        s2t*(-17512 + 14805*z3) + (-375892 + 621675*z3)*pow2(Mt) - 900*(-1226 +
        495*z3)*pow2(Mst1)*pow2(s2t)) + pow3(Dmst12)*(-30*Mst1*s2t*(-1117238 +
        877575*z3)*pow2(Mt) - 2250*Mt*(-5570 + 333*z3)*pow2(Mst1)*pow2(s2t) + (
        -41715182 + 38174625*z3)*pow3(Mt) + 3000*(-2722 + 2259*z3)*pow3(Mst1)*
        pow3(s2t)) - 6000*Dmst12*(81*Mt*(-406 + 285*z3) + 8*Mst1*s2t*(-694 +
        477*z3))*pow2(Mt)*pow4(Mst2) + 48000*(623 + 963*z3)*pow3(Mt)*pow6(Mst2)
        )) + 196*Dmglst1*(30000*Mst1*pow2(Dmsqst1)*pow2(Mt)*(-19*Mst1*s2t*pow2(
        Dmst12)*pow2(Mst2) + 19*Mst1*s2t*pow3(Dmst12) + 19*Dmst12*Mst1*s2t*
        pow4(Mst2) - 132*Mt*pow6(Mst2)) + 600*Dmsqst1*pow3(Mst1)*(-2*Mt*pow2(
        Dmst12)*pow2(Mst2)*(-1516*Mst1*Mt*s2t + 7614*pow2(Mt) + 3525*pow2(Mst1)
        *pow2(s2t)) + pow3(Dmst12)*(-564*Mst1*s2t*pow2(Mt) + 14100*Mt*pow2(
        Mst1)*pow2(s2t) + 15228*pow3(Mt) + 125*pow3(Mst1)*pow3(s2t)) + 4*
        Dmst12*(3807*Mt - 1375*Mst1*s2t)*pow2(Mt)*pow4(Mst2) + 29200*pow3(Mt)*
        pow6(Mst2)) + pow5(Mst1)*(Mt*pow2(Dmst12)*pow2(Mst2)*(40*Mst1*Mt*s2t*(-
        4511549 + 3729375*z3) + (45149198 - 35285625*z3)*pow2(Mt) - 47250*(-430
        + 207*z3)*pow2(Mst1)*pow2(s2t)) + pow3(Dmst12)*(2*Mst1*s2t*(-28188929 +
        23099625*z3)*pow2(Mt) + 225*Mt*(-160136 + 100785*z3)*pow2(Mst1)*pow2(
        s2t) + (115176444 - 85920750*z3)*pow3(Mt) + 1125*(22174 - 18791*z3)*
        pow3(Mst1)*pow3(s2t)) + 40*Dmst12*(150*Mst1*s2t*(-19856 + 17667*z3) +
        Mt*(-5136871 + 3912300*z3))*pow2(Mt)*pow4(Mst2) + 6000*(-31142 + 22653*
        z3)*pow3(Mt)*pow6(Mst2))) + 2*pow2(Dmglst1)*(9243360*pow2(Dmsqst1)*
        pow3(Mt)*pow6(Mst2) + 35280*Dmsqst1*Mt*pow2(Mst1)*(-(pow2(Dmst12)*pow2(
        Mst2)*(7940*Mst1*Mt*s2t + 914*pow2(Mt) + 225*pow2(Mst1)*pow2(s2t))) + (
        7940*Mst1*Mt*s2t + 914*pow2(Mt) + 450*pow2(Mst1)*pow2(s2t))*pow3(
        Dmst12) + 2*Dmst12*Mt*(457*Mt + 3970*Mst1*s2t)*pow4(Mst2) + 4760*pow2(
        Mt)*pow6(Mst2)) + pow4(Mst1)*(Mt*pow2(Dmst12)*pow2(Mst2)*(196*Mst1*Mt*
        s2t*(-263717842 + 218365875*z3) + (27129768542 - 22522586625*z3)*pow2(
        Mt) + 22050*(-63802 + 92895*z3)*pow2(Mst1)*pow2(s2t)) + pow3(Dmst12)*(-
        392*Mst1*s2t*(-384557822 + 319453875*z3)*pow2(Mt) + 66150*Mt*(-150737 +
        112420*z3)*pow2(Mst1)*pow2(s2t) + (-25287306692 + 22556819250*z3)*pow3(
        Mt) + 3675*(1908362 - 1581075*z3)*pow3(Mst1)*pow3(s2t)) + 392*Dmst12*(
        20*Mst1*s2t*(-6041999 + 5054400*z3) + Mt*(-73908751 + 57368250*z3))*
        pow2(Mt)*pow4(Mst2) + 3920*(-8223692 + 6125625*z3)*pow3(Mt)*pow6(Mst2))
        ))/(1.9845e6*pow3(Mt)*pow6(Mst1)*pow6(Mst2));

   return result;
}

/**
 * @return returns the susy log^1 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double H32q2g::calc_coef_at_as2_no_sm_logs_log1() const {

   const double result =
      (4*(2*Mst1*pow3(Dmglst1)*(-2*Mst1*Mt*pow2(Dmst12)*pow2(Mst2)*(82320*
        Dmsqst1*Mt*s2t + 346639*Mt*s2t*pow2(Mst1) + 555463*Mst1*pow2(Mt) +
        1051785*pow2(s2t)*pow3(Mst1)) + Mst1*pow3(Dmst12)*(164640*Dmsqst1*s2t*
        pow2(Mt) + 8778304*s2t*pow2(Mst1)*pow2(Mt) + 4207140*Mt*pow2(s2t)*pow3(
        Mst1) + 1110926*Mst1*pow3(Mt) + 661255*pow3(s2t)*pow4(Mst1)) + 2*
        Dmst12*Mst1*(555463*Mst1*Mt + 82320*Dmsqst1*s2t - 3695874*s2t*pow2(
        Mst1))*pow2(Mt)*pow4(Mst2) - 4704*(10*Dmsqst1 + 141*pow2(Mst1))*pow3(
        Mt)*pow6(Mst2)) + pow2(Dmglst1)*(-246960*pow2(Dmsqst1)*pow3(Mt)*pow6(
        Mst2) - 17640*Dmsqst1*pow2(Mst1)*pow2(Mt)*((17*Mt - 10*Mst1*s2t)*pow2(
        Dmst12)*pow2(Mst2) + (-17*Mt + 10*Mst1*s2t)*pow3(Dmst12) + Dmst12*(-17*
        Mt + 10*Mst1*s2t)*pow4(Mst2) - 35*Mt*pow6(Mst2)) + pow4(Mst1)*(Mt*pow2(
        Dmst12)*pow2(Mst2)*(-4088560*Mst1*Mt*s2t + 968629*pow2(Mt) + 1853670*
        pow2(Mst1)*pow2(s2t)) + pow3(Dmst12)*(7502488*Mst1*s2t*pow2(Mt) -
        3134775*Mt*pow2(Mst1)*pow2(s2t) + 2019394*pow3(Mt) + 689430*pow3(Mst1)*
        pow3(s2t)) - 196*Dmst12*(20187*Mt - 3442*Mst1*s2t)*pow2(Mt)*pow4(Mst2)
        - 8199072*pow3(Mt)*pow6(Mst2))) + 98*Dmglst1*(600*Dmsqst1*pow2(Mt)*
        pow3(Mst1)*((6*Mt + Mst1*s2t)*pow2(Dmst12)*pow2(Mst2) + (-6*Mt + 3*
        Mst1*s2t)*pow3(Dmst12) - Dmst12*(6*Mt + 5*Mst1*s2t)*pow4(Mst2) - 25*Mt*
        pow6(Mst2)) + 1500*Mst1*pow2(Dmsqst1)*pow2(Mt)*(-(Mst1*s2t*pow2(Dmst12)
        *pow2(Mst2)) + Mst1*s2t*pow3(Dmst12) + Dmst12*Mst1*s2t*pow4(Mst2) + 3*
        Mt*pow6(Mst2)) - pow5(Mst1)*(Mt*pow2(Dmst12)*pow2(Mst2)*(29758*Mst1*Mt*
        s2t + 6677*pow2(Mt) + 22350*pow2(Mst1)*pow2(s2t)) + pow3(Dmst12)*(-
        34587*Mst1*s2t*pow2(Mt) - 18045*Mt*pow2(Mst1)*pow2(s2t) + 22414*pow3(
        Mt) + 3325*pow3(Mst1)*pow3(s2t)) - 8*Dmst12*(4471*Mt + 6875*Mst1*s2t)*
        pow2(Mt)*pow4(Mst2) + 50800*pow3(Mt)*pow6(Mst2))) - 49*(750*pow2(
        Dmsqst1)*pow2(Mst1)*pow2(Mt)*(-((Mt + 8*Mst1*s2t)*pow2(Dmst12)*pow2(
        Mst2)) + (Mt + 8*Mst1*s2t)*pow3(Dmst12) + Dmst12*(Mt + 8*Mst1*s2t)*
        pow4(Mst2) - 14*Mt*pow6(Mst2)) + 150*Dmsqst1*Mt*pow4(Mst1)*(pow2(
        Dmst12)*pow2(Mst2)*(80*Mst1*Mt*s2t + 17*pow2(Mt) - 180*pow2(Mst1)*pow2(
        s2t)) - 2*(20*Mst1*Mt*s2t + 187*pow2(Mt) - 90*pow2(Mst1)*pow2(s2t))*
        pow3(Dmst12) + 20*Dmst12*Mt*(17*Mt - 6*Mst1*s2t)*pow4(Mst2) + 1080*
        pow2(Mt)*pow6(Mst2)) - pow6(Mst1)*(-2*Mt*pow2(Dmst12)*pow2(Mst2)*(
        25850*Mst1*Mt*s2t + 21033*pow2(Mt) + 75*pow2(Mst1)*pow2(s2t)) + pow3(
        Dmst12)*(68264*Mst1*s2t*pow2(Mt) + 5700*Mt*pow2(Mst1)*pow2(s2t) +
        36107*pow3(Mt) + 250*pow3(Mst1)*pow3(s2t)) + 200*Dmst12*(131*Mt + 100*
        Mst1*s2t)*pow2(Mt)*pow4(Mst2) + 400*(-533 + 54*z3)*pow3(Mt)*pow6(Mst2))
        )))/(33075.*pow3(Mt)*pow6(Mst1)*pow6(Mst2));

   return result;
}

/**
 * @return returns the susy log^2 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double H32q2g::calc_coef_at_as2_no_sm_logs_log2() const {

   const double result =
      (-8*(7*Dmglst1*pow3(Mst1)*(2*Mt*pow2(Dmst12)*pow2(Mst2)*(-1304*Mst1*Mt*
        s2t + 888*pow2(Mt) + 645*pow2(Mst1)*pow2(s2t)) + pow3(Dmst12)*(1544*
        Mst1*s2t*pow2(Mt) - 2250*Mt*pow2(Mst1)*pow2(s2t) - 1640*pow3(Mt) + 275*
        pow3(Mst1)*pow3(s2t)) - 8*Dmst12*(239*Mt - 35*Mst1*s2t)*pow2(Mt)*pow4(
        Mst2) - 6520*pow3(Mt)*pow6(Mst2)) + pow2(Dmglst1)*pow2(Mst1)*(Mt*pow2(
        Dmst12)*pow2(Mst2)*(-4088*Mst1*Mt*s2t + 22884*pow2(Mt) - 8925*pow2(
        Mst1)*pow2(s2t)) + pow3(Dmst12)*(42728*Mst1*s2t*pow2(Mt) - 1155*Mt*
        pow2(Mst1)*pow2(s2t) - 56772*pow3(Mt) + 3360*pow3(Mst1)*pow3(s2t)) +
        28*Dmst12*(393*Mt - 1234*Mst1*s2t)*pow2(Mt)*pow4(Mst2) - 3948*pow3(Mt)*
        pow6(Mst2)) + 8*Mst1*pow3(Dmglst1)*(-3*Mt*pow2(Dmst12)*pow2(Mst2)*(-75*
        Mst1*Mt*s2t + 1122*pow2(Mt) + 560*pow2(Mst1)*pow2(s2t)) + pow3(Dmst12)*
        (1671*Mst1*s2t*pow2(Mt) + 3360*Mt*pow2(Mst1)*pow2(s2t) + 3366*pow3(Mt)
        + 140*pow3(Mst1)*pow3(s2t)) + 3*Dmst12*(1122*Mt - 707*Mst1*s2t)*pow2(
        Mt)*pow4(Mst2) + 2163*pow3(Mt)*pow6(Mst2)) + 7*(-(Mt*pow2(Dmst12)*pow2(
        Mst2)*(1760*Mst1*Mt*s2t + 538*pow2(Mt) + 105*pow2(Mst1)*pow2(s2t))*
        pow4(Mst1)) + pow3(Dmst12)*(1032*Mst1*s2t*pow2(Mt) + 480*Mt*pow2(Mst1)*
        pow2(s2t) + 644*pow3(Mt) - 45*pow3(Mst1)*pow3(s2t))*pow4(Mst1) + 40*
        Dmst12*(11*Mt + 61*Mst1*s2t)*pow2(Mt)*pow4(Mst1)*pow4(Mst2) + 10*pow3(
        Mt)*(45*pow2(Dmsqst1) - 90*Dmsqst1*pow2(Mst1) - 442*pow4(Mst1))*pow6(
        Mst2))))/(315.*pow3(Mt)*pow4(Mst1)*pow6(Mst2));

   return result;
}

/**
 * @return returns the susy log^3 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double H32q2g::calc_coef_at_as2_no_sm_logs_log3() const {

   const double result =
      -298.6666666666667;

   return result;
}

}        // hierarchies
}        // himalaya
