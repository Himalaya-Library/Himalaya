// ====================================================================
// This file is part of Himalaya.
//
// Himalaya is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================

#include "H5.hpp"
#include "Hierarchies.hpp"
#include "Constants.hpp"
#include "Utils.hpp"
#include <cmath>
#include <type_traits>

namespace himalaya{
namespace hierarchies{

/**
 * Constructor
 * @param flagMap the flagMap for the truncation of expansion variables
 * @param Al4p a double alpha_s/4/Pi
 * @param beta a double which is the mixing angle beta
 * @param Dmglst1 a double Mgl - Mst1
 * @param lmMt a double log((<renormalization scale> / Mt)^2)
 * @param lmMst1 a double log((<renormalization scale> / Mst1)^2)
 * @param lmMst2 a double log((<renormalization scale> / Mst2)^2)
 * @param lmMsq a double log((<renormalization scale> / Msq)^2)
 * @param Mt a double top/bottom quark mass
 * @param Mst1 a double stop 1 mass
 * @param Mst2 a double stop 2 mass
 * @param Msq a double average squark mass w/o the stop quark
 * @param MuSUSY a double mu parameter
 * @param s2t a double 2 times the sine of the stop/sbottom quark mixing angle
 * @param mdrFlag an int 0 for DR and 1 for MDR scheme
 * @param oneLoopFlag an int flag to consider the one-loop expansion terms
 * @param twoLoopFlag an int flag to consider the two-loop expansion terms
 * @param threeLoopFlag an int flag to consider the three-loop expansion terms
 */
H5::H5(const std::map<unsigned int, unsigned int>& flagMap, double Al4p, double beta, double Dmglst1,
                 double lmMt, double lmMst1, double lmMst2, double lmMsq, double Mt, double Mst1,
                 double Mst2, double Msq, double MuSUSY,
                 double s2t,
                 int mdrFlag, int oneLoopFlag, int twoLoopFlag, int threeLoopFlag){
   // abbrev for tan(beta) and sin(beta)
   Tbeta = tan(beta);
   Cbeta = cos(beta);
   Sbeta = sin(beta);
   this -> Dmglst1 = Dmglst1;
   this -> lmMt = lmMt;
   this -> lmMst1 = lmMst1;
   this -> lmMst2 = lmMst2;
   this -> lmMsq = lmMsq;
   this -> Mt = Mt;
   this -> Mst1 = Mst1;
   this -> Mst2 = Mst2;
   this -> Msq = Msq;
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
   // expansion flags
   xDmglst1 = flagMap.at(ExpansionDepth::xxDmglst1);
   xMsq = flagMap.at(ExpansionDepth::xxMsq);
}

/**
 * @return The diagonal (1, 1) matrix element of the Higgs mass matrix as a double for the hierarchy 'H5'
 */
double H5::getS1() const {
   return (pow2(Mt)*pow2(MuSUSY)*(243*oneLoopFlag*pow2(s2t)*(2 - lmMst1 + lmMst2 -
        (2*(lmMst1 - lmMst2)*pow2(Mst1)*(pow2(Mst1) + pow2(Mst2)))/pow4(Mst2))
        + (Al4p*((324*s2t*twoLoopFlag*(8*Dmglst1*Mt*(Mst1*(12 - 2*lmMst1 + 6*
        lmMst2 - pow2(lmMst1) + pow2(lmMst2))*pow2(Mst2) + (10 + 12*lmMst2 + 4*
        lmMst1*(-2 + 3*lmMst2) - 9*pow2(lmMst1) - 3*pow2(lmMst2))*pow3(Mst1)) +
        8*Mst1*Mt*(Mst1*(8 - 2*lmMst1 + 6*lmMst2 - pow2(lmMst1) + pow2(lmMst2))
        *pow2(Mst2) + (6 + 4*lmMst1*(-3 + lmMst2) + 16*lmMst2 - 5*pow2(lmMst1)
        + pow2(lmMst2))*pow3(Mst1)) - 8*Dmglst1*s2t*((-(lmMst2*(1 + lmMst2)) +
        pow2(lmMst1))*pow2(Mst1)*pow2(Mst2) + (2 + 3*lmMst1 - 3*lmMst2 + pow2(
        lmMst1) - pow2(lmMst2))*pow4(Mst1) - (-1 + lmMst1)*pow4(Mst2)) + Mst1*
        s2t*(4*(2 + 11*lmMst2 - 2*lmMst1*(5 + 4*lmMst2) + pow2(lmMst1) + 7*
        pow2(lmMst2))*pow2(Mst1)*pow2(Mst2) + (1 + 46*lmMst2 - 2*lmMst1*(23 +
        32*lmMst2) + 20*pow2(lmMst1) + 44*pow2(lmMst2))*pow4(Mst1) + 4*(7 + 7*
        lmMst2 - 2*lmMst1*(1 + lmMst2) + 2*pow2(lmMst2))*pow4(Mst2))))/Mst1 +
        108*z2*(12*s2t*twoLoopFlag*(4*Mst1*Mt*pow2(Mst2) + 4*Mt*pow3(Mst1) + 4*
        Dmglst1*(-3*Mt*pow2(Mst1) + Mt*pow2(Mst2) + Mst1*s2t*pow2(Mst2) + s2t*
        pow3(Mst1)) - s2t*pow4(Mst2)) + xDmglst1*pow2(Dmglst1)*(24*s2t*
        twoLoopFlag*(-12*Mst1*Mt + s2t*pow2(Mst1) + s2t*pow2(Mst2)) + Al4p*
        threeLoopFlag*(32*(55 + 16*lmMst2)*pow2(Mt) + pow2(s2t)*(-((2043 + 360*
        lmMsq + 1196*lmMst1 - 1900*lmMst2)*pow2(Mst1)) + 4*(390 - 90*lmMsq -
        139*lmMst1 + 315*lmMst2)*pow2(Mst2) + (198*pow4(Mst2))/pow2(Mst1)) +
        s2t*(16*(-1187 + 360*lmMsq + 171*lmMst1 - 747*lmMst2)*Mst1*Mt + (150*
        pow2(Mst1)*pow2(Mst2)*(-8*Mst1*Mt + s2t*pow2(Mst2)))/pow4(Msq) + (60*(-
        16*Mst1*Mt*pow2(Mst2) - 18*s2t*pow2(Mst1)*pow2(Mst2) + 112*Mt*pow3(
        Mst1) + s2t*pow4(Mst2)))/pow2(Msq))))) + (3*xDmglst1*pow2(Dmglst1)*(-
        432*s2t*twoLoopFlag*pow4(Msq)*(-8*Mst1*Mt*pow2(Mst2) - s2t*(4 + lmMst2
        - pow2(lmMst1) + pow2(lmMst2))*pow2(Mst1)*pow2(Mst2) + 4*Mt*(4 + 7*
        lmMst2 - lmMst1*(7 + 6*lmMst2) + 3*(pow2(lmMst1) + pow2(lmMst2)))*pow3(
        Mst1) + s2t*(-4 + 9*lmMst1 - 9*lmMst2 + pow2(lmMst1) - pow2(lmMst2))*
        pow4(Mst1) - (-2 + lmMst1)*s2t*pow4(Mst2)) + Al4p*threeLoopFlag*(64*
        pow2(Mst1)*pow2(Mt)*(5483 - 18*B4 - 9*DN + 418*lmMst1 + 102*pow2(
        lmMst1) + 2*lmMst2*(739 - 42*lmMst1 + 9*pow2(lmMst1)) + 54*(5 + lmMst1)
        *pow2(lmMst2) + 24*lmMt*(17 + 17*lmMst2 - lmMst1*(17 + 6*lmMst2) + 3*
        pow2(lmMst1) + 3*pow2(lmMst2)) - 78*pow3(lmMst1) + 6*pow3(lmMst2))*
        pow4(Msq) - 48*Mst1*Mt*s2t*(((5875 + 1020*lmMsq - 478*lmMst1 - 610*
        lmMst2 - 576*pow2(lmMst1))*pow2(Mst2) + 2*pow2(Mst1)*(6129 + 684*B4 -
        18*DN - 3630*lmMsq + 6*(-251 + 600*lmMsq)*lmMst1 + lmMst2*(6914 - 1122*
        lmMst1 + 720*lmMsq*(-5 + 3*lmMst1) - 1179*pow2(lmMst1)) - 3*(997 + 360*
        lmMsq)*pow2(lmMst1) - 9*(-393 + 120*lmMsq + 125*lmMst1)*pow2(lmMst2) +
        1353*pow3(lmMst1) + 951*pow3(lmMst2)))*pow4(Msq) + 5*(983 - 159*lmMst1
        + 180*lmMsq*(-2 + lmMst1 - lmMst2) + 519*lmMst2 - 90*pow2(lmMst1) + 90*
        pow2(lmMst2))*pow2(Mst2)*pow4(Mst1) + 20*pow2(Msq)*((311 - 66*lmMst1 +
        36*lmMsq*(-2 + lmMst1 - lmMst2) + 138*lmMst2 - 18*pow2(lmMst1) + 18*
        pow2(lmMst2))*pow2(Mst1)*pow2(Mst2) + (29 - 6*lmMst1 + 36*lmMsq*(-2 +
        3*lmMst1 - 3*lmMst2) + 78*lmMst2 + 360*lmMst1*lmMst2 - 234*pow2(lmMst1)
        - 126*pow2(lmMst2))*pow4(Mst1)) - 90*pow2(Mst1)*pow4(Mst2)) + 3*pow2(
        s2t)*((40*(160 - 213*lmMst1 + 18*lmMsq*(8 + lmMst1 - lmMst2) + 69*
        lmMst2 - 9*pow2(lmMst1) + 9*pow2(lmMst2))*pow2(Msq)*pow2(Mst1) + 2*(
        6415 + 60*lmMsq*(37 - 18*lmMst1) - 1314*lmMst2 + 10*lmMst1*(-97 + 66*
        lmMst2) + 540*pow2(lmMsq) - 234*pow2(lmMst1) + 54*pow2(lmMst2))*pow4(
        Msq))*pow4(Mst2) + pow4(Mst1)*(216*(129.29128086419752 - (52*B4)/9. - (
        22*DN)/9. - (380*lmMsq)/3. - lmMst1*(150.47685185185185 - 160*lmMsq +
        10*pow2(lmMsq)) + (5*(-57 + 8*lmMsq)*pow2(lmMst1))/4. + lmMst2*(
        324.6990740740741 - 160*lmMsq - (3383*lmMst1)/18. + 10*pow2(lmMsq) + (
        833*pow2(lmMst1))/9.) - ((-8947 + 360*lmMsq + 4868*lmMst1)*pow2(lmMst2)
        )/36. - (721*pow3(lmMst1))/27. + (1873*pow3(lmMst2))/27.)*pow4(Msq) + (
        10585 - 8070*lmMst1 + 5370*lmMst2 - 900*pow2(lmMst1) + 900*(lmMsq*(3 +
        2*lmMst1 - 2*lmMst2) + pow2(lmMst2)))*pow4(Mst2)) + pow2(Mst1)*(8*pow2(
        Msq)*pow2(Mst2)*(-90*(25 - 11*lmMst1 + 2*lmMsq*(-5 + 9*lmMst1 - 9*
        lmMst2) + 21*lmMst2 - 9*pow2(lmMst1) + 9*pow2(lmMst2))*pow2(Mst1) +
        pow2(Msq)*(265 + 36*B4 + 30*DN - 1935*lmMsq + lmMst1*(2178 + 270*lmMsq
        - 270*pow2(lmMsq)) + 135*pow2(lmMsq) + 3*lmMst2*(353 - 180*lmMsq + 86*
        lmMst1 + 90*pow2(lmMsq) - 79*pow2(lmMst1)) + 30*(-13 + 9*lmMsq)*pow2(
        lmMst1) - 3*(-372 + 90*lmMsq + 49*lmMst1)*pow2(lmMst2) - 65*pow3(
        lmMst1) + 449*pow3(lmMst2))) - 90*pow6(Mst2))))))/(pow2(Mst1)*pow4(Msq)
        )))/pow4(Mst2) + 2*threeLoopFlag*pow2(Al4p)*(972*pow2(s2t)*(
        83.61265432098766 + 4*B4 - (4*DN)/9. - (185*lmMsq)/9. + (25*pow2(lmMsq)
        )/3. - (lmMst1*(3781 - 300*lmMsq + 180*pow2(lmMsq)))/108. + (lmMst2*(
        14065 - 3498*lmMst1 + 60*lmMsq*(-35 + 12*lmMst1) + 180*pow2(lmMsq) -
        18*pow2(lmMst1)))/108. + (6.361111111111111 - (5*lmMsq)/3.)*pow2(
        lmMst1) - ((-2193 + 180*lmMsq + 442*lmMst1)*pow2(lmMst2))/36. +
        Dmglst1*((7 - 2*lmMst1 - 60*lmMsq*(-5 + 6*lmMst1) - 226*lmMst2 + 220*
        lmMst1*lmMst2 + 180*pow2(lmMsq) + 178*pow2(lmMst1) + 18*pow2(lmMst2))/(
        18.*Mst1) + (10*Mst1*(82 - 93*lmMst1 + 6*lmMsq*(4 + 3*lmMst1 - 3*
        lmMst2) + 69*lmMst2 - 9*pow2(lmMst1) + 9*pow2(lmMst2)))/(27.*pow2(Msq))
        ) + ((1 - 2*lmMst2)*pow2(Mst2))/(3.*pow2(Mst1)) + (11*pow3(lmMst1))/18.
         + (71*pow3(lmMst2))/6. + (Mst1*(-(Mst1*(26.575942386831276 - (76*B4)/
        9. + (2*DN)/9. + (35*lmMsq)/3. - 5*pow2(lmMsq) + lmMst1*(
        224.2274074074074 - (380*lmMsq)/9. + (40*pow2(lmMsq))/3.) - (149*pow2(
        lmMst1))/90. - (lmMst2*(372457 - 242670*lmMst1 + 1500*lmMsq*(-47 + 24*
        lmMst1) + 18000*pow2(lmMsq) + 17700*pow2(lmMst1)))/1350. + ((-17709 +
        2400*lmMsq + 6940*lmMst1)*pow2(lmMst2))/90. + (8*pow3(lmMst1))/27. - (
        1736*pow3(lmMst2))/27.)) + (2*Dmglst1*(586 + 36*B4 + 30*DN - 495*lmMsq
        + 135*pow2(lmMsq) - 6*lmMst1*(-73 - 45*lmMsq + 45*pow2(lmMsq)) + 3*
        lmMst2*(229 - 180*lmMsq + 86*lmMst1 + 90*pow2(lmMsq) - 79*pow2(lmMst1))
        + 18*(-43 + 15*lmMsq)*pow2(lmMst1) - 3*(-372 + 90*lmMsq + 49*lmMst1)*
        pow2(lmMst2) - 65*pow3(lmMst1) + 449*pow3(lmMst2)))/27.))/pow2(Mst2) -
        ((Mst1*(265.6519022158834 - (76*B4)/9. + (2*DN)/9. - (25*lmMsq)/2. +
        lmMst1*(426.37458616780043 - (605*lmMsq)/9. + (40*pow2(lmMsq))/3.) - (
        66.24060846560846 - (40*lmMsq)/3.)*pow2(lmMst1) - lmMst2*(
        411.7079195011338 - (1386139*lmMst1)/3780. + (5*lmMsq*(-121 + 96*
        lmMst1))/9. + (40*pow2(lmMsq))/3. + (298*pow2(lmMst1))/3.) - (
        298.6850529100529 - 40*lmMsq - (2174*lmMst1)/9.)*pow2(lmMst2) + (80*
        pow3(lmMst1))/27. - (3920*pow3(lmMst2))/27.) + Dmglst1*(
        112.49099794238683 - (8*B4)/3. - (20*DN)/9. + lmMst1*(
        254.24382716049382 - 120*lmMsq + 20*pow2(lmMsq)) + (129.57407407407408
         - 20*lmMsq)*pow2(lmMst1) - lmMst2*(
        211.57716049382717 - 120*lmMsq - (1439*lmMst1)/27. + 20*pow2(lmMsq) + (
        338*pow2(lmMst1))/9.) + ((-9491 + 1080*lmMsq + 6636*lmMst1)*pow2(
        lmMst2))/54. + (38*pow3(lmMst1))/9. - (806*pow3(lmMst2))/9.))*pow3(
        Mst1))/pow4(Mst2) - (-((836 - 353*lmMst2 + 30*lmMsq*(7 - 4*lmMst1 + 4*
        lmMst2) + lmMst1*(143 + 750*lmMst2) - 315*pow2(lmMst1) - 435*pow2(
        lmMst2))*pow2(Mst1)*pow2(Mst2)) + 100*Dmglst1*(55 + lmMst1 + 6*lmMsq*(-
        5 + 7*lmMst1 - 7*lmMst2) + 29*lmMst2 - 21*pow2(lmMst1) + 21*pow2(
        lmMst2))*pow3(Mst1) + (1525 + 30*lmMsq*(-25 + 44*lmMst1 - 44*lmMst2) +
        1928*lmMst2 - 2*lmMst1*(589 + 900*lmMst2) + 240*pow2(lmMst1) + 1560*
        pow2(lmMst2))*pow4(Mst1) + (939 + 30*lmMsq*(-12 + 5*lmMst1 - 5*lmMst2)
        + 760*lmMst2 - 50*lmMst1*(8 + 3*lmMst2) + 150*pow2(lmMst2))*pow4(Mst2))
        /(135.*pow2(Msq)*pow2(Mst2)) + (-2*(9723887 + 420*lmMsq*(-13109 +
        11760*lmMst1 - 10920*lmMst2) + 13243860*lmMst2 - 164640*lmMst1*(47 +
        30*lmMst2) - 176400*pow2(lmMsq) + 4762800*pow2(lmMst2))*pow2(Mst1)*
        pow2(Mst2) + 514500*Dmglst1*(-18*Mst1*pow2(Mst2) + (485 - 346*lmMst1 +
        12*lmMsq*(1 + 10*lmMst1 - 10*lmMst2) + 334*lmMst2 - 60*pow2(lmMst1) +
        60*pow2(lmMst2))*pow3(Mst1)) + (41220947 - 420*lmMsq*(12479 + 4830*
        lmMst1 - 5670*lmMst2) - 1081710*lmMst2 + 210*lmMst1*(30109 + 123480*
        lmMst2) - 176400*pow2(lmMsq) - 11951100*pow2(lmMst1) - 14156100*pow2(
        lmMst2))*pow4(Mst1) - (7399303 + 420*lmMsq*(-9571 + 2940*lmMst1 - 3780*
        lmMst2) + 6571740*lmMst2 - 82320*lmMst1*(31 + 15*lmMst2) + 176400*pow2(
        lmMsq) + 1411200*pow2(lmMst2))*pow4(Mst2))/(1.11132e7*pow4(Msq)) + (8*
        OepS2*(204*pow2(Mst1)*pow2(Mst2) + 370*pow4(Mst1) + 27*pow4(Mst2)) -
        12*T1ep*(204*pow2(Mst1)*pow2(Mst2) + 370*pow4(Mst1) + 27*pow4(Mst2)) -
        27*S2*(36*(49 + 34*lmMst1 - 34*lmMst2)*pow2(Mst1)*pow2(Mst2) + (4138 +
        2220*lmMst1 - 2220*lmMst2)*pow4(Mst1) + 81*(-11 + 2*lmMst1 - 2*lmMst2)*
        pow4(Mst2)) + (2430*xMsq*(1 - 2*(lmMsq + z2))*pow2(Msq)*(pow2(Mst1)*
        pow2(Mst2) + 2*lmMst1*(-1 + shiftst1)*pow2(Mst1)*pow2(Mst2) - 2*lmMst2*
        (-1 + shiftst1)*pow2(Mst1)*pow2(Mst2) + shiftst2*pow2(Mst1)*(2*pow2(
        Mst1) + pow2(Mst2)) - 2*lmMst1*(1 - 2*shiftst1 + shiftst2)*pow4(Mst1) +
        2*lmMst2*(1 - 2*shiftst1 + shiftst2)*pow4(Mst1) + pow4(Mst2) -
        shiftst1*(2*pow2(Mst1)*pow2(Mst2) + 2*pow4(Mst1) + pow4(Mst2))))/pow2(
        Mst1) + (243*(-1 + 2*lmMst2)*shiftst3*(2*(1 - 2*lmMst1 + 2*lmMst2)*
        pow2(Mst2)*pow4(Mst1) + 2*(1 - lmMst1 + lmMst2)*pow2(Mst1)*pow4(Mst2) +
        (2 - 6*lmMst1 + 6*lmMst2)*pow6(Mst1) + pow6(Mst2)))/pow2(Mst1))/(729.*
        pow4(Mst2))) + (288*Mst1*pow2(Mt)*(2*Dmglst1*(735 - 6*B4 - 3*DN - 78*
        lmMst1 + 6*pow2(lmMst1) + 6*lmMst2*(85 - 6*lmMst1 + pow2(lmMst1)) + 18*
        (7 + lmMst1)*pow2(lmMst2) + 24*lmMt*(2 + 3*lmMst2 - lmMst1*(3 + 2*
        lmMst2) + pow2(lmMst1) + pow2(lmMst2)) - 26*pow3(lmMst1) + 2*pow3(
        lmMst2)) + Mst1*(417 - 6*B4 - 3*DN - 162*lmMst1 + 6*lmMst2*(87 - 8*
        lmMst1 + pow2(lmMst1)) + 18*(8 + lmMst1)*pow2(lmMst2) + 24*lmMt*pow2(1
        - lmMst1 + lmMst2) - 26*pow3(lmMst1) + 2*pow3(lmMst2))) - 8*z4*(-3*
        pow2(Mst1)*(540*pow2(Mt) + 13*pow2(Mst2)*pow2(s2t)) - 54*Dmglst1*(Mt*
        s2t*(507*pow2(Mst1) - 241*pow2(Mst2)) + 60*Mst1*pow2(Mt) - 265*Mst1*(
        pow2(Mst1) + pow2(Mst2))*pow2(s2t)) - 135*xDmglst1*pow2(Dmglst1)*(492*
        Mst1*Mt*s2t + 12*pow2(Mt) - (5*pow2(Mst1) + 53*pow2(Mst2))*pow2(s2t)) +
        54*Mt*s2t*(241*Mst1*pow2(Mst2) + 313*pow3(Mst1)) + pow2(s2t)*(127*pow4(
        Mst1) - 3267*pow4(Mst2))) - (2*z3*(4*pow2(Mst1)*(108*Mst1*Mt*s2t*(23*(
        25 + 9*lmMst1 - 9*lmMst2)*pow2(Mst1) + (133 - 9*lmMst1 + 9*lmMst2)*
        pow2(Mst2)) - pow2(Mst1)*(216*(23 - 3*lmMst1 + 3*lmMst2)*pow2(Mt) + 3*(
        8783 + 1134*lmMst1 - 1134*lmMst2)*pow2(Mst2)*pow2(s2t)) + pow2(s2t)*(-
        14*(4829 + 243*lmMst1 - 243*lmMst2)*pow4(Mst1) - 4428*pow4(Mst2))) -
        27*xDmglst1*pow2(Dmglst1)*(216*Mst1*Mt*s2t*((51 - 52*lmMst1 + 52*
        lmMst2)*pow2(Mst1) + 25*pow2(Mst2)) - 8*pow2(Mst1)*(4*(181 + 3*lmMst1 -
        3*lmMst2)*pow2(Mt) - 3*(14 + 23*lmMst1 - 23*lmMst2)*pow2(Mst2)*pow2(
        s2t)) + pow2(s2t)*(12*(1145 + 14*lmMst1 - 14*lmMst2)*pow4(Mst1) - 1143*
        pow4(Mst2))) + 54*Dmglst1*Mst1*(12*pow2(Mst1)*((52 + 8*lmMst1 - 8*
        lmMst2)*pow2(Mt) + (-35 - 46*lmMst1 + 46*lmMst2)*pow2(Mst2)*pow2(s2t))
        + Mt*s2t*(8*(58 - 9*lmMst1 + 9*lmMst2)*Mst1*pow2(Mst2) + 24*(154 + 225*
        lmMst1 - 225*lmMst2)*pow3(Mst1)) + pow2(s2t)*(-8*(764 + 69*lmMst1 - 69*
        lmMst2)*pow4(Mst1) + 195*pow4(Mst2)))))/pow2(Mst1) + z2*(144*Mt*(Mst1*(
        3*Dmglst1*(8*(55 + 16*lmMst2)*Mt + (-687 + 540*lmMsq + 92*lmMst1 - 860*
        lmMst2)*Mst1*s2t) + Mst1*(12*(55 + 16*lmMst2)*Mt + (2407 + 180*lmMsq -
        408*lmMst1 + 408*lmMst2)*Mst1*s2t)) + (3*s2t*(Mst1*pow2(Mst2)*(-40*
        pow2(Msq)*pow2(Mst1) + (627 - 60*lmMsq - 64*lmMst1 + 192*lmMst2)*pow4(
        Msq) - 15*pow4(Mst1)) + Dmglst1*(pow2(Mst2)*((627 - 60*lmMsq - 64*
        lmMst1 + 192*lmMst2)*pow4(Msq) - 75*pow4(Mst1)) + 120*pow2(Msq)*(-(
        pow2(Mst1)*pow2(Mst2)) + 2*pow4(Mst1)))))/pow4(Msq)) + (pow2(s2t)*(36*
        Dmglst1*(150*pow4(Mst2)*pow5(Mst1) + pow4(Msq)*(12*(390 - 90*lmMsq -
        139*lmMst1 + 315*lmMst2)*pow2(Mst2)*pow3(Mst1) + 594*Mst1*pow4(Mst2) +
        (1613 - 1080*lmMsq - 3684*lmMst1 + 5796*lmMst2)*pow5(Mst1)) - 60*pow2(
        Msq)*(-3*pow3(Mst1)*pow4(Mst2) + 14*pow2(Mst2)*pow5(Mst1))) + 1350*
        pow4(Mst2)*pow6(Mst1) - 1080*pow2(Msq)*(-3*pow4(Mst1)*pow4(Mst2) + 4*
        pow2(Mst2)*pow6(Mst1)) - pow4(Msq)*(12*(311 + 1080*lmMsq - 3348*lmMst2
        - 108*(1 + 2*lmMst2)*shiftst3 + 36*lmMst1*(61 + 6*shiftst3))*pow2(Mst2)
        *pow4(Mst1) + 108*(339 - 30*lmMsq + 149*lmMst2 - 12*(1 + lmMst2)*
        shiftst3 + lmMst1*(-37 + 12*shiftst3))*pow2(Mst1)*pow4(Mst2) + (69349 +
        12960*lmMsq - 107568*lmMst2 - 1296*(1 + 3*lmMst2)*shiftst3 + 432*
        lmMst1*(217 + 9*shiftst3))*pow6(Mst1) - 648*(-1 + shiftst3)*pow6(Mst2))
        ))/(pow2(Mst1)*pow4(Msq))) + s2t*(12*pow2(z2)*(72*xDmglst1*pow2(
        Dmglst1)*(120*Mst1*Mt - s2t*(pow2(Mst1) + pow2(Mst2))) + 576*Mt*(2*
        Mst1*pow2(Mst2) + 11*pow3(Mst1)) - 144*Dmglst1*(-84*Mt*pow2(Mst1) - 8*
        Mt*pow2(Mst2) + Mst1*s2t*pow2(Mst2) + s2t*pow3(Mst1)) - s2t*(852*pow2(
        Mst1)*pow2(Mst2) + 1018*pow4(Mst1) + 315*pow4(Mst2))) + (2*Mt*(3*
        Dmglst1*(2*(pow2(Mst1)*(13969 - 11880*B4 + 108*DN - 4320*lmMsq + (-
        16728 - 3240*lmMsq*(1 + lmMsq))*lmMst1 + 2160*pow2(lmMsq) + 54*(-407 +
        300*lmMsq)*pow2(lmMst1) + 12*lmMst2*(4094 + 2205*lmMst1 - 90*lmMsq*(1 +
        24*lmMst1) + 270*pow2(lmMsq) + 804*pow2(lmMst1)) + 18*(583 + 540*lmMsq
        + 680*lmMst1)*pow2(lmMst2) - 16848*pow3(lmMst1) - 5040*pow3(lmMst2)) +
        3*pow2(Mst2)*(11337 - 408*B4 - 12*DN - 6480*lmMsq + 720*pow2(lmMsq) -
        30*lmMst1*(-73 - 36*lmMsq + 12*pow2(lmMsq)) + 6*lmMst2*(2331 - 420*
        lmMsq - 118*lmMst1 + 60*pow2(lmMsq) - 60*pow2(lmMst1)) + 18*(-13 + 20*
        lmMsq)*pow2(lmMst1) - 6*(-797 + 60*lmMsq + 4*lmMst1)*pow2(lmMst2) -
        200*pow3(lmMst1) + 584*pow3(lmMst2)))*pow4(Msq) - 720*pow2(Msq)*(2*(41
        - 10*lmMst1 + 6*lmMsq*(-2 + lmMst1 - lmMst2) + 22*lmMst2 - 3*pow2(
        lmMst1) + 3*pow2(lmMst2))*pow2(Mst1)*pow2(Mst2) + 2*(29 + 48*lmMst2 -
        24*pow2(lmMst1) - 6*(lmMst1*(6 - 5*lmMst2) + lmMsq*(2 - 3*lmMst1 + 3*
        lmMst2) + pow2(lmMst2)))*pow4(Mst1) - pow4(Mst2)) - 5*pow2(Mst2)*((-251
        + 12*lmMsq - 12*lmMst2)*pow2(Mst1)*pow2(Mst2) + 6*(818 - 117*lmMst1 +
        180*lmMsq*(-2 + lmMst1 - lmMst2) + 477*lmMst2 - 90*pow2(lmMst1) + 90*
        pow2(lmMst2))*pow4(Mst1) + (11 - 12*lmMsq + 12*lmMst2)*pow4(Mst2))) +
        Mst1*(2*(9*pow2(Mst2)*(13317 - 408*B4 - 12*DN - 4320*lmMsq + 720*pow2(
        lmMsq) - 6*lmMst1*(79 - 180*lmMsq + 60*pow2(lmMsq)) + 6*lmMst2*(2135 -
        420*lmMsq - 118*lmMst1 + 60*pow2(lmMsq) - 60*pow2(lmMst1)) + 6*(-167 +
        60*lmMsq)*pow2(lmMst1) - 6*(-797 + 60*lmMsq + 4*lmMst1)*pow2(lmMst2) -
        200*pow3(lmMst1) + 584*pow3(lmMst2)) + pow2(Mst1)*(107299 - 19224*B4 -
        108*DN - 32400*lmMsq + 6480*pow2(lmMsq) - 24*lmMst1*(6103 - 1755*lmMsq
        + 405*pow2(lmMsq)) + 18*(-4379 + 1260*lmMsq)*pow2(lmMst1) + 12*lmMst2*(
        20882 + 699*lmMst1 - 270*lmMsq*(17 + 8*lmMst1) + 810*pow2(lmMsq) + 54*
        pow2(lmMst1)) + 18*(5641 + 180*lmMsq + 540*lmMst1)*pow2(lmMst2) -
        18072*pow3(lmMst1) + 7704*pow3(lmMst2)))*pow4(Msq) - 720*pow2(Msq)*(2*(
        34 - 8*lmMst1 + 6*lmMsq*(-2 + lmMst1 - lmMst2) + 20*lmMst2 - 3*pow2(
        lmMst1) + 3*pow2(lmMst2))*pow2(Mst1)*pow2(Mst2) + 4*(14 - 21*lmMst1 +
        lmMsq*(-6 + 9*lmMst1 - 9*lmMst2) + 27*lmMst2 + 9*lmMst1*lmMst2 - 9*
        pow2(lmMst1))*pow4(Mst1) - 3*pow4(Mst2)) - 15*(6*(134 - 9*lmMst1 + 36*
        lmMsq*(-2 + lmMst1 - lmMst2) + 81*lmMst2 - 18*pow2(lmMst1) + 18*pow2(
        lmMst2))*pow2(Mst2)*pow4(Mst1) + (-107 + 12*lmMsq - 12*lmMst2)*pow2(
        Mst1)*pow4(Mst2) + (11 - 12*lmMsq + 12*lmMst2)*pow6(Mst2)))))/pow4(Msq)
        ))/pow4(Mst2))))/1944.;
}

/**
 * @return The diagonal (2, 2) matrix element of the Higgs mass matrix as a double for the hierarchy 'H5'
 */
double H5::getS2() const {
   return -(oneLoopFlag*((4*Mt*MuSUSY*s2t*(4*(lmMst1 - lmMst2)*pow2(Mt) + ((-2 -
        lmMst1 + lmMst2)*pow2(Mst1) + (2 - lmMst1 + lmMst2)*pow2(Mst2))*pow2(
        s2t)))/Tbeta + pow2(Mt)*pow2(s2t)*(8*(lmMst1 - lmMst2)*(pow2(Mst1) -
        pow2(Mst2)) + 4*pow2(MuSUSY)*(2 - lmMst1 + lmMst2 - (2*(lmMst1 -
        lmMst2)*pow2(Mst1)*(pow2(Mst1) + pow2(Mst2)))/pow4(Mst2)) - (4*pow2(
        MuSUSY)*(2 - lmMst1 + lmMst2 - (2*(lmMst1 - lmMst2)*pow2(Mst1)*(pow2(
        Mst1) + pow2(Mst2)))/pow4(Mst2)))/pow2(Sbeta)) + 16*(lmMst1 + lmMst2 -
        2*lmMt)*pow4(Mt) - (-4*pow2(Mst1)*pow2(Mst2) + (2 + lmMst1 - lmMst2)*
        pow4(Mst1) + (2 - lmMst1 + lmMst2)*pow4(Mst2))*pow4(s2t)))/32. - Al4p*(
        (twoLoopFlag*(-36*s2t*Tbeta*pow2(Mt)*(-8*pow2(Mst1)*pow2(Mst2)*(-(
        Dmglst1*(1 + 2*(lmMst1 + lmMst2))*s2t*pow2(Mst2)) + Mt*(8 - 2*lmMst1 +
        6*lmMst2 - pow2(lmMst1) + pow2(lmMst2))*pow2(MuSUSY))*pow2(Sbeta) +
        pow2(MuSUSY)*(8*Dmglst1*Mt*(Mst1*(12 - 2*lmMst1 + 6*lmMst2 - pow2(
        lmMst1) + pow2(lmMst2))*pow2(Mst2) + (10 + 12*lmMst2 + 4*lmMst1*(-2 +
        3*lmMst2) - 9*pow2(lmMst1) - 3*pow2(lmMst2))*pow3(Mst1)) + 8*Mst1*Mt*(
        Mst1*(8 - 2*lmMst1 + 6*lmMst2 - pow2(lmMst1) + pow2(lmMst2))*pow2(Mst2)
        + (6 + 4*lmMst1*(-3 + lmMst2) + 16*lmMst2 - 5*pow2(lmMst1) + pow2(
        lmMst2))*pow3(Mst1)) - 8*Dmglst1*s2t*((-(lmMst2*(1 + lmMst2)) + pow2(
        lmMst1))*pow2(Mst1)*pow2(Mst2) + (2 + 3*lmMst1 - 3*lmMst2 + pow2(
        lmMst1) - pow2(lmMst2))*pow4(Mst1) - (-1 + lmMst1)*pow4(Mst2)) + Mst1*
        s2t*(4*(2 + 11*lmMst2 - 2*lmMst1*(5 + 4*lmMst2) + pow2(lmMst1) + 7*
        pow2(lmMst2))*pow2(Mst1)*pow2(Mst2) + (1 + 46*lmMst2 - 2*lmMst1*(23 +
        32*lmMst2) + 20*pow2(lmMst1) + 44*pow2(lmMst2))*pow4(Mst1) + 4*(7 + 7*
        lmMst2 - 2*lmMst1*(1 + lmMst2) + 2*pow2(lmMst2))*pow4(Mst2)))) - 9*
        Tbeta*pow2(Sbeta)*pow4(Mst2)*pow4(s2t)*(4*(-2*lmMst1*(3 + 2*lmMst2) +
        pow2(lmMst1) + 3*(-4 - lmMst2 + pow2(lmMst2)))*pow2(Mst2)*pow3(Mst1) +
        4*Mst1*(7 + 7*lmMst2 - 2*lmMst1*(1 + lmMst2) + 2*pow2(lmMst2))*pow4(
        Mst2) + 8*Dmglst1*((2 - 2*lmMst1 + lmMst2 - pow2(lmMst1) + pow2(lmMst2)
        )*pow2(Mst1)*pow2(Mst2) + (-3 - 2*lmMst1 + lmMst2 + pow2(lmMst1) -
        pow2(lmMst2))*pow4(Mst1) + (-1 + lmMst1)*pow4(Mst2)) + (13 + lmMst1*(26
        - 8*lmMst2) - 14*lmMst2 + 12*pow2(lmMst1) - 4*pow2(lmMst2))*pow5(Mst1))
        + 4*Tbeta*pow2(Mt)*pow2(Sbeta)*(8*pow2(Mt)*(-36*Dmglst1*(1 + lmMt)*
        pow2(Mst1)*pow2(Mst2) - 18*(1 + lmMt)*pow2(Mst2)*pow3(Mst1) + 36*
        Dmglst1*(1 + lmMst2*(5 - 2*lmMt) - 2*lmMst1*(1 + lmMst2 - lmMt) - 3*
        lmMt + 2*pow2(lmMst2))*pow4(Mst1) - 2*Dmglst1*(1 + 12*lmMst1 + 6*lmMt)*
        pow4(Mst2) + 9*Mst1*(3 + lmMst1 - lmMt + 2*lmMst1*lmMt + 2*lmMst2*(1 +
        lmMt) + pow2(lmMst1) - 6*pow2(lmMt))*pow4(Mst2) + 9*Mst1*pow2(lmMst2)*(
        2*pow4(Mst1) + pow4(Mst2)) - 18*(lmMst1*(1 + lmMst2 - lmMt) + lmMst2*(-
        2 + lmMt) + lmMt)*pow5(Mst1)) + 72*Mst1*Mt*s2t*((2*(lmMst1 - lmMst2 +
        2*lmMst1*lmMst2 - 2*lmMst1*lmMt + 2*lmMst2*lmMt)*pow2(Mst2) + (6 + 4*
        lmMst1*(-3 + lmMst2) + 16*lmMst2 - 5*pow2(lmMst1))*pow2(MuSUSY) + pow2(
        lmMst2)*(-4*pow2(Mst2) + pow2(MuSUSY)))*pow3(Mst1) + Dmglst1*((12 - 2*
        lmMst1 + 6*lmMst2 - pow2(lmMst1) + pow2(lmMst2))*pow2(Mst2)*pow2(
        MuSUSY) + pow2(Mst1)*(2*(-2 + 3*lmMst1 - 3*lmMst2)*(1 + 2*lmMst2 - 2*
        lmMt)*pow2(Mst2) + (10 - 8*lmMst1 + 12*(1 + lmMst1)*lmMst2 - 9*pow2(
        lmMst1) - 3*pow2(lmMst2))*pow2(MuSUSY)) + (20 + 10*lmMst1*(5 + 2*lmMst2
        - 2*lmMt) + 8*lmMt + lmMst2*(-58 + 20*lmMt) - 20*pow2(lmMst2))*pow4(
        Mst1) + 2*(8 - 5*lmMst1 + lmMst2 + (4 - 2*lmMst1 + 2*lmMst2)*lmMt)*
        pow4(Mst2)) + 2*((3 + lmMst2 + 2*(1 + lmMst2)*lmMt - lmMst1*(3 + 2*
        lmMt))*Mst1*pow4(Mst2) + (4 + lmMst1*(5 + 2*lmMst2 - 2*lmMt) + lmMst2*(
        -5 + 2*lmMt) - 2*pow2(lmMst2))*pow5(Mst1))) - 9*pow2(s2t)*(4*(Mst1*(-4*
        lmMst2*pow2(Mst2) + pow2(lmMst1)*pow2(Mst2) - 3*pow2(lmMst2)*pow2(Mst2)
        - 7*pow2(MuSUSY) - 7*lmMst2*pow2(MuSUSY) - 2*pow2(lmMst2)*pow2(MuSUSY)
        + 2*lmMst1*(1 + lmMst2)*(pow2(Mst2) + pow2(MuSUSY)))*pow4(Mst2) + pow3(
        Mst1)*(-((2 + 11*lmMst2 - 2*lmMst1*(5 + 4*lmMst2) + pow2(lmMst1) + 7*
        pow2(lmMst2))*pow2(Mst2)*pow2(MuSUSY)) + (2 + 2*lmMst1*(-1 + lmMst2) +
        6*lmMst2 - 3*pow2(lmMst1) + pow2(lmMst2))*pow4(Mst2))) - (8*(1 +
        lmMst1)*pow2(Mst2) + (1 + 46*lmMst2 - 2*lmMst1*(23 + 32*lmMst2) + 20*
        pow2(lmMst1) + 44*pow2(lmMst2))*pow2(MuSUSY))*pow5(Mst1) + 8*Dmglst1*(
        lmMst2*(2*pow2(Mst2) - 3*pow2(MuSUSY))*pow4(Mst1) + pow2(MuSUSY)*((-(
        lmMst2*(1 + lmMst2)) + pow2(lmMst1))*pow2(Mst1)*pow2(Mst2) + (2 + pow2(
        lmMst1) - pow2(lmMst2))*pow4(Mst1) + pow4(Mst2)) - lmMst1*((4*pow2(
        Mst2) - 3*pow2(MuSUSY))*pow4(Mst1) + (2*pow2(Mst2) + pow2(MuSUSY))*
        pow4(Mst2)) + pow6(Mst2)))) + 36*Mt*pow2(Sbeta)*(-2*Mst1*Tbeta*pow2(
        Mst2)*pow3(s2t)*(-2*(10 + 8*lmMst1 - 4*(1 + lmMst1)*lmMst2 + 3*pow2(
        lmMst1) + pow2(lmMst2))*pow2(Mst2)*pow3(Mst1) + 2*Mst1*(8 - 2*lmMst1 +
        6*lmMst2 - pow2(lmMst1) + pow2(lmMst2))*pow4(Mst2) + Dmglst1*(-2*(14 +
        lmMst1*(4 - 12*lmMst2) + 7*pow2(lmMst1) + 5*pow2(lmMst2))*pow2(Mst1)*
        pow2(Mst2) + (19 + 18*lmMst1 - 18*lmMst2)*pow4(Mst1) + 2*(12 - 2*lmMst1
        + 6*lmMst2 - pow2(lmMst1) + pow2(lmMst2))*pow4(Mst2)) + (3 + 2*lmMst1 -
        2*lmMst2)*pow5(Mst1)) + MuSUSY*(-(Mst1*(4*Mst1*Mt*pow2(Mst2)*(8*(2 +
        lmMst2 + lmMt + lmMst2*lmMt - lmMst1*(2 + lmMt))*pow2(Mt) + 3*(-8 + 2*
        lmMst1 - 6*lmMst2 + pow2(lmMst1) - pow2(lmMst2))*pow2(Mst2)*pow2(s2t))
        + 8*Mt*(4*(2 + lmMst2 + lmMst1*(-2 + lmMst2 - 2*lmMt) + lmMt + 2*
        lmMst2*lmMt - pow2(lmMst2))*pow2(Mt) + 3*(1 + lmMst1*(5 - 2*lmMst2) -
        5*lmMst2 + 2*pow2(lmMst1))*pow2(Mst2)*pow2(s2t))*pow3(Mst1) + s2t*(
        lmMst1*(16 - 32*lmMst2)*pow2(Mt) + (7 + lmMst1*(6 + 32*lmMst2))*pow2(
        Mst2)*pow2(s2t) + 16*(pow2(lmMst1) + pow2(lmMst2))*(pow2(Mt) - pow2(
        Mst2)*pow2(s2t)) - 2*lmMst2*(8*pow2(Mt) + pow2(Mst2)*pow2(s2t)))*pow4(
        Mst1) - 4*s2t*((-8*lmMst2*pow2(Mt) + 2*pow2(lmMst1)*pow2(Mt) - 6*pow2(
        lmMst2)*pow2(Mt) + 7*pow2(Mst2)*pow2(s2t) + 7*lmMst2*pow2(Mst2)*pow2(
        s2t) + 2*pow2(lmMst2)*pow2(Mst2)*pow2(s2t) + 2*lmMst1*(1 + lmMst2)*(2*
        pow2(Mt) - pow2(Mst2)*pow2(s2t)))*pow4(Mst2) + pow2(Mst1)*(-4*(-1 -
        lmMst2 - 2*lmMst1*lmMst2 + pow2(lmMst1) + pow2(lmMst2))*pow2(Mst2)*
        pow2(Mt) + (-5 - 8*lmMst1 + 4*lmMst2 - 6*lmMst1*lmMst2 + pow2(lmMst1) +
        5*pow2(lmMst2))*pow2(s2t)*pow4(Mst2))) + 6*Mt*(1 + lmMst1*(18 - 8*
        lmMst2) - 18*lmMst2 + 8*pow2(lmMst1))*pow2(s2t)*pow5(Mst1))) + 2*
        Dmglst1*(4*s2t*pow2(Mst1)*pow2(Mst2)*(4*(1 + lmMst2)*pow2(Mt) + (1 -
        lmMst1 + lmMst2 - pow2(lmMst1) + pow2(lmMst2))*pow2(Mst2)*pow2(s2t)) -
        4*Mt*(-4*(-5 + lmMst2 - 4*lmMt - 4*lmMst2*lmMt + lmMst1*(3 - 3*lmMst2 +
        4*lmMt) + 3*pow2(lmMst2))*pow2(Mt) + 3*(1 + lmMst1*(3 - 6*lmMst2) - 3*
        lmMst2 + 4*pow2(lmMst1) + 2*pow2(lmMst2))*pow2(Mst2)*pow2(s2t))*pow3(
        Mst1) - 16*(5 + lmMst2 + (2 + lmMst2)*lmMt - lmMst1*(3 + lmMt))*Mst1*
        pow2(Mst2)*pow3(Mt) + 4*s2t*((4 - 8*lmMst1 + 8*lmMst2)*pow2(Mt) + (-2 -
        3*lmMst1 + 2*lmMst2)*pow2(Mst2)*pow2(s2t))*pow4(Mst1) + 8*(1 - 2*
        lmMst1)*s2t*pow2(Mt)*pow4(Mst2) + 6*Mst1*Mt*(12 - 2*lmMst1 + 6*lmMst2 -
        pow2(lmMst1) + pow2(lmMst2))*pow2(s2t)*pow4(Mst2) + 3*Mt*(15 - 6*lmMst2
        + 6*lmMst1*(1 + 4*lmMst2) - 16*pow2(lmMst1) - 8*pow2(lmMst2))*pow2(s2t)
        *pow5(Mst1) + 4*(-1 + lmMst1)*pow3(s2t)*pow6(Mst2))))))/(216.*Mst1*
        Tbeta*pow2(Sbeta)*pow4(Mst2)) + z2*((twoLoopFlag*(Mt*(8*(5*Dmglst1*
        pow2(Mst1) - (Dmglst1 + Mst1)*pow2(Mst2) + pow3(Mst1))*pow3(s2t) + (4*
        MuSUSY*((4*Dmglst1*Mst1 + pow2(Mst1) - pow2(Mst2))*pow3(s2t) + 6*Mt*
        pow2(s2t)*(Dmglst1 + Mst1 - (4*Dmglst1*pow2(Mst1)*(pow2(Mst1) + pow2(
        Mst2)))/pow4(Mst2)) + (16*(3*Dmglst1 + Mst1)*pow2(Mst1)*pow3(Mt))/pow4(
        Mst2)))/Tbeta) + (4*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*(4*Dmglst1*Mst1*(
        pow2(Mst1) + pow2(Mst2)) - pow4(Mst2)))/pow4(Mst2) + (pow2(Mt)*(Mt*(32*
        (2*Dmglst1*(2*Mt - 5*Mst1*s2t) + Mst1*(Mt - 2*Mst1*s2t))*pow3(Mst1) +
        16*s2t*((Dmglst1 + Mst1)*pow2(Mst2)*pow2(MuSUSY) - 3*Dmglst1*pow2(Mst1)
        *(4*pow2(Mst2) + pow2(MuSUSY)) + (-4*pow2(Mst2) + pow2(MuSUSY))*pow3(
        Mst1))) - (4*s2t*pow2(MuSUSY)*(4*Mst1*Mt*pow2(Mst2) + 4*Mt*pow3(Mst1) +
        4*Dmglst1*(-3*Mt*pow2(Mst1) + Mt*pow2(Mst2) + Mst1*s2t*pow2(Mst2) +
        s2t*pow3(Mst1)) - s2t*pow4(Mst2)))/pow2(Sbeta)))/pow4(Mst2) + (pow2(
        Mst1) - pow2(Mst2))*(4*Dmglst1*Mst1 + pow2(Mst1) - pow2(Mst2))*pow4(
        s2t)))/6. - xDmglst1*pow2(Dmglst1)*((twoLoopFlag*(8*s2t*(9*MuSUSY*s2t +
        40*Mt*Tbeta)*pow2(Mt)*pow2(Sbeta)*pow3(Mst1) + pow2(Mst2)*pow2(s2t)*(-
        4*Tbeta*pow2(Mt)*pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 4*Mt*MuSUSY*s2t*
        pow2(Mst2)*pow2(Sbeta) + Tbeta*pow2(s2t)*pow2(Sbeta)*pow4(Mst2)) - 24*
        Mst1*Mt*(-3*Mt*MuSUSY*pow2(Mst2)*pow2(s2t)*pow2(Sbeta) - 2*s2t*Tbeta*
        pow2(Mt)*(pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 2*pow2(Mst2)*pow2(Sbeta)) +
        4*MuSUSY*pow2(Sbeta)*pow3(Mt) + Tbeta*pow2(Sbeta)*pow3(s2t)*pow4(Mst2))
        - Tbeta*pow2(Mst1)*(4*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta)
        ) + 96*pow2(Sbeta)*pow4(Mt) + pow2(Sbeta)*pow4(Mst2)*pow4(s2t))))/(3.*
        Tbeta*pow2(Sbeta)*pow4(Mst2)) + Al4p*threeLoopFlag*(-(Mt*pow3(s2t)*((
        527.5555555555555 - 160*lmMsq - 76*lmMst1 + 332*lmMst2)*Mst1 - (240*
        pow3(Mst1))/pow2(Msq) + ((2.2962962962962963 - 312*lmMst1 + 312*lmMst2)
        *pow3(Mst1))/pow2(Mst2) + (20*pow2(Mst2)*(4*Mst1*pow2(Msq) + 5*pow3(
        Mst1)))/(3.*pow4(Msq)))) + pow2(Mt)*pow2(s2t)*(227.7037037037037 + (80*
        lmMst1)/9. + (304*lmMst2)/9. + (4*(-1861 + 906*lmMst1 - 1866*lmMst2)*
        pow2(Mst1))/(9.*pow2(Mst2)) - (22*pow2(Mst2))/pow2(Mst1) - pow2(MuSUSY)
        *(10/(3.*pow2(Msq)) + 11/pow2(Mst1) + (86.66666666666667 - 20*lmMsq - (
        278*lmMst1)/9. + 70*lmMst2 - (60*pow2(Mst1))/pow2(Msq))/pow2(Mst2) +
        pow2(Mst1)*(25/(3.*pow4(Msq)) - (113.5 + 20*lmMsq + (598*lmMst1)/9. - (
        950*lmMst2)/9.)/pow4(Mst2)))) - (Mt*MuSUSY*(pow3(s2t)*(
        75.66666666666667 - 20*lmMsq - (278*lmMst1)/9. + 70*lmMst2 - pow2(Mst1)
        *(190/(3.*pow2(Msq)) + (3603 + 640*lmMst1 - 640*lmMst2)/(18.*pow2(Mst2)
        )) + pow2(Mst2)*(10/(3.*pow2(Msq)) + 11/pow2(Mst1) + (25*pow2(Mst1))/(
        3.*pow4(Msq)))) - Mt*pow2(s2t)*((80*Mst1)/pow2(Msq) + (Mst1*(
        1582.6666666666667 - 480*lmMsq - 228*lmMst1 + 996*lmMst2 - (640*pow2(
        Mst1))/pow2(Msq)))/pow2(Mst2) + pow3(Mst1)*(100/pow4(Msq) - (2*(-7153 +
        2160*lmMsq + 5238*lmMst1 - 8694*lmMst2))/(9.*pow4(Mst2)))) + (4*pow2(
        Mt)*(s2t*((2777 + 120*lmMst1 + 456*lmMst2)*pow2(Mst1)*pow2(Mst2) + (-
        8389 + 5556*lmMst1 - 10740*lmMst2)*pow4(Mst1) - 297*pow4(Mst2)) + (12*
        Mt*pow3(Mst1)*(-600*pow2(Msq)*pow2(Mst1) + (1391 - 360*lmMsq - 66*
        lmMst1 + 813*lmMst2 + 21*lmMt)*pow4(Msq) + 30*pow4(Mst2)))/pow4(Msq)))/
        (27.*pow2(Mst1)*pow4(Mst2))))/Tbeta + (4*(360/pow2(Msq) + 297/pow2(
        Mst1) + (604 + 432*lmMst2 - 432*lmMt)/pow2(Mst2) + (90*(11*pow2(Mst1) +
        pow2(Mst2)))/pow4(Msq) + (12*(-1584 + 450*lmMsq + 211*lmMst1 - 736*
        lmMst2 - 3*lmMt)*pow2(Mst1))/pow4(Mst2) - (12*(55 + 16*lmMst2)*pow2(
        MuSUSY))/pow4(Mst2))*pow4(Mt))/27. + ((4*pow2(Mst2)*(291 - 90*lmMsq -
        139*lmMst1 + 315*lmMst2 + (15*pow2(Mst2))/pow2(Msq)) + (198*pow4(Mst2))
        /pow2(Mst1) + pow2(Mst1)*(-4965 + 360*lmMsq - 84*lmMst1 - 620*lmMst2 -
        (1200*pow2(Mst2))/pow2(Msq) + (150*pow4(Mst2))/pow4(Msq)))*pow4(s2t))/
        72. + (pow2(Mt)*((3*pow2(MuSUSY)*(150*s2t*pow2(Mst2)*(-8*Mst1*Mt + s2t*
        pow2(Mst2))*pow4(Mst1) + 60*s2t*pow2(Msq)*pow2(Mst1)*(-16*Mst1*Mt*pow2(
        Mst2) - 18*s2t*pow2(Mst1)*pow2(Mst2) + 112*Mt*pow3(Mst1) + s2t*pow4(
        Mst2)) + pow4(Msq)*(4*pow2(Mst1)*(8*(55 + 16*lmMst2)*pow2(Mt) + (390 -
        90*lmMsq - 139*lmMst1 + 315*lmMst2)*pow2(Mst2)*pow2(s2t)) + 16*(-1187 +
        360*lmMsq + 171*lmMst1 - 747*lmMst2)*Mt*s2t*pow3(Mst1) - (2043 + 360*
        lmMsq + 1196*lmMst1 - 1900*lmMst2)*pow2(s2t)*pow4(Mst1) + 198*pow2(s2t)
        *pow4(Mst2))))/pow2(Sbeta) + 16*Mt*s2t*pow3(Mst1)*(180*pow2(Msq)*(pow2(
        Mst2)*pow2(MuSUSY) - pow2(Mst1)*(20*pow2(Mst2) + 7*pow2(MuSUSY))) + ((
        42755 - 10800*lmMsq - 1500*lmMst1 + 20334*lmMst2 + 366*lmMt)*pow2(Mst1)
        + 3*(2761 - 720*lmMsq - 132*lmMst1 + 1626*lmMst2 + 42*lmMt)*pow2(Mst2)
        + 3*(1187 - 360*lmMsq - 171*lmMst1 + 747*lmMst2)*pow2(MuSUSY))*pow4(
        Msq) + 45*(pow2(Mst1)*(5*pow2(Mst2)*pow2(MuSUSY) - 4*pow4(Mst2)) + 4*
        pow6(Mst2)))))/(54.*pow2(Mst1)*pow4(Msq)*pow4(Mst2))))) + xDmglst1*
        pow2(Dmglst1)*((twoLoopFlag*((-72*Mst1*Mt*pow3(s2t)*(-((8 + 7*lmMst2 -
        lmMst1*(7 + 6*lmMst2) + 3*(pow2(lmMst1) + pow2(lmMst2)))*pow2(Mst1)*
        pow2(Mst2)) + 6*(1 + 2*lmMst1 - 2*lmMst2)*pow4(Mst1) + 2*pow4(Mst2)))/
        pow2(Mst2) + 9*(-((8 - 2*lmMst1 + lmMst2 - pow2(lmMst1) + pow2(lmMst2))
        *pow2(Mst1)*pow2(Mst2)) + (6 + 8*lmMst1 - 7*lmMst2 - pow2(lmMst1) +
        pow2(lmMst2))*pow4(Mst1) - (-2 + lmMst1)*pow4(Mst2))*pow4(s2t) + (4*Mt*
        (12*pow3(Mt)*(-3*(1 + lmMt)*pow2(Mst1)*pow2(Mst2) + 3*(7 + lmMst2*(19 -
        6*lmMt) - 6*lmMst1*(1 + lmMst2 - lmMt) - 13*lmMt + 6*pow2(lmMst2))*
        pow4(Mst1) + (4 - 3*lmMst1)*pow4(Mst2)) + Mt*s2t*(-2*Mst1*Mt*(-36*pow2(
        Mst2)*pow2(MuSUSY) + 18*pow2(Mst1)*((5 + lmMst2*(68 - 20*lmMt) - 10*
        lmMst1*(5 + 2*lmMst2 - 2*lmMt) - 18*lmMt + 20*pow2(lmMst2))*pow2(Mst1)
        + (5 - 3*lmMst1 + 3*lmMst2)*(1 + 2*lmMst2 - 2*lmMt)*pow2(Mst2) + (4 +
        7*lmMst2 - lmMst1*(7 + 6*lmMst2) + 3*(pow2(lmMst1) + pow2(lmMst2)))*
        pow2(MuSUSY)) + (-107 + 30*lmMst1 - 30*lmMt)*pow4(Mst2)) - (9*pow2(
        MuSUSY)*(8*Mst1*Mt*pow2(Mst2) + s2t*(4 + lmMst2 - pow2(lmMst1) + pow2(
        lmMst2))*pow2(Mst1)*pow2(Mst2) - 4*Mt*(4 + 7*lmMst2 - lmMst1*(7 + 6*
        lmMst2) + 3*pow2(lmMst1) + 3*pow2(lmMst2))*pow3(Mst1) + s2t*(4 - 9*
        lmMst1 + 9*lmMst2 - pow2(lmMst1) + pow2(lmMst2))*pow4(Mst1) + (-2 +
        lmMst1)*s2t*pow4(Mst2)))/pow2(Sbeta)) - 3*Mt*pow2(s2t)*(3*(2*(6 - 6*
        lmMst1 + 5*lmMst2)*pow2(Mst2) + (-4 + 9*lmMst1 - 9*lmMst2 + pow2(
        lmMst1) - pow2(lmMst2))*pow2(MuSUSY))*pow4(Mst1) - 3*(-2 + lmMst1)*
        pow2(MuSUSY)*pow4(Mst2) - pow2(Mst1)*(3*(4 + lmMst2 - pow2(lmMst1) +
        pow2(lmMst2))*pow2(Mst2)*pow2(MuSUSY) + (5 - 6*(lmMst1 + lmMst2))*pow4(
        Mst2)) + (11 - 6*lmMst1)*pow6(Mst2)) + (MuSUSY*(-2*pow3(Mst1)*(27*Mt*(6
        + 7*lmMst2 - lmMst1*(7 + 6*lmMst2) + 3*pow2(lmMst1) + 3*pow2(lmMst2))*
        pow2(Mst2)*pow2(s2t) + 2*(61 - 90*lmMst2 + 3*lmMst1*(-5 + 18*lmMst2 -
        18*lmMt) + 105*lmMt + 54*lmMst2*lmMt - 54*pow2(lmMst2))*pow3(Mt)) - 9*(
        4*(-7 + 6*lmMst1 - 6*lmMst2)*s2t*pow2(Mt) + (9*lmMst1 - 8*lmMst2)*pow2(
        Mst2)*pow3(s2t))*pow4(Mst1) + 6*(11 - 6*lmMst1)*s2t*pow2(Mt)*pow4(Mst2)
        + 4*Mst1*((-61 + 15*lmMst1 - 15*lmMt)*pow2(Mst2)*pow3(Mt) + 27*Mt*pow2(
        s2t)*pow4(Mst2)) + 9*pow2(Mst1)*(4*(1 + lmMst2)*s2t*pow2(Mst2)*pow2(Mt)
        + (6 - lmMst1 + lmMst2 - pow2(lmMst1) + pow2(lmMst2))*pow3(s2t)*pow4(
        Mst2)) - 54*Mt*(lmMst2*(19 + 3*lmMst2) - lmMst1*(19 + 6*lmMst2) + 3*
        pow2(lmMst1))*pow2(s2t)*pow5(Mst1) + 9*(-2 + lmMst1)*pow3(s2t)*pow6(
        Mst2)))/Tbeta))/pow4(Mst2)))/(54.*pow2(Mst1)) + (Al4p*threeLoopFlag*(
        800*pow2(Mt)*(-(pow2(s2t)*(972470 - 1728*B4 - 864*DN + 9720*lmMsq + 60*
        (1225 - 324*lmMsq)*lmMst1 + 19440*pow2(lmMsq) + 13320*pow2(lmMst1) +
        12*lmMst2*(18367 - 1620*lmMsq + 660*lmMst1 + 432*pow2(lmMst1)) + 72*(
        809 + 168*lmMst1)*pow2(lmMst2) + 3456*lmMt*(17*(1 + lmMst2) - lmMst1*(
        17 + 6*lmMst2) + 3*(pow2(lmMst1) + pow2(lmMst2))) - 36*((30*(5 + 60*
        lmMsq - 66*lmMst1 + 6*lmMst2))/pow2(Msq) + (4717 - 540*lmMsq*(-2 +
        lmMst1) - 252*lmMst2 + 2*lmMst1*(-661 + 69*lmMst2) + 270*pow2(lmMsq) +
        75*pow2(lmMst1) + 27*pow2(lmMst2))/pow2(Mst1))*pow2(Mst2) - 14400*pow3(
        lmMst1) - 2880*pow3(lmMst2) - 6*pow2(Mst1)*((540*(35 - 40*lmMsq + 22*
        lmMst1 + 18*lmMst2))/pow2(Msq) + (2897 - 5178*lmMst2 + 540*lmMsq*(85 -
        60*lmMst1 + 54*lmMst2) + 1620*pow2(lmMsq) - 432*(-283 + 12*lmMst2 + 48*
        lmMt)*pow2(lmMst1) + 6*lmMst1*(2345 + 9696*lmMt + 144*lmMst2*(-305 +
        48*lmMt) - 9504*pow2(lmMst2)) + 149616*pow2(lmMst2) - 576*lmMt*(54 +
        101*lmMst2 + 36*pow2(lmMst2)) + 24768*pow3(lmMst1) + 37440*pow3(lmMst2)
        )/pow2(Mst2)) - 9*pow2(MuSUSY)*((6415 + 60*lmMsq*(37 - 18*lmMst1) -
        1314*lmMst2 + 10*lmMst1*(-97 + 66*lmMst2) + 540*pow2(lmMsq) - 234*pow2(
        lmMst1) + 54*pow2(lmMst2))/pow2(Mst1) + (20*(160 - 213*lmMst1 + 18*
        lmMsq*(8 + lmMst1 - lmMst2) + 69*lmMst2 - 9*pow2(lmMst1) + 9*pow2(
        lmMst2)) + (4*(-90*(25 - 11*lmMst1 + 2*lmMsq*(-5 + 9*lmMst1 - 9*lmMst2)
        + 21*lmMst2 - 9*pow2(lmMst1) + 9*pow2(lmMst2))*pow2(Mst1) + pow2(Msq)*(
        265 + 36*B4 + 30*DN - 1935*lmMsq + lmMst1*(2178 + 270*lmMsq - 270*pow2(
        lmMsq)) + 135*pow2(lmMsq) + 3*lmMst2*(353 - 180*lmMsq + 86*lmMst1 + 90*
        pow2(lmMsq) - 79*pow2(lmMst1)) + 30*(-13 + 9*lmMsq)*pow2(lmMst1) - 3*(-
        372 + 90*lmMsq + 49*lmMst1)*pow2(lmMst2) - 65*pow3(lmMst1) + 449*pow3(
        lmMst2))))/pow2(Mst2))/pow2(Msq) - (45*pow2(Mst2))/pow4(Msq) + 108*
        pow2(Mst1)*((10585 - 8070*lmMst1 + 5370*lmMst2 - 900*pow2(lmMst1) +
        900*(lmMsq*(3 + 2*lmMst1 - 2*lmMst2) + pow2(lmMst2)))/(216.*pow4(Msq))
        + (129.29128086419752 - (52*B4)/9. - (22*DN)/9. - (380*lmMsq)/3. -
        lmMst1*(150.47685185185185 - 160*lmMsq + 10*pow2(lmMsq)) + (5*(-57 + 8*
        lmMsq)*pow2(lmMst1))/4. + lmMst2*(324.6990740740741 - 160*lmMsq - (
        3383*lmMst1)/18. + 10*pow2(lmMsq) + (833*pow2(lmMst1))/9.) - ((-8947 +
        360*lmMsq + 4868*lmMst1)*pow2(lmMst2))/36. - (721*pow3(lmMst1))/27. + (
        1873*pow3(lmMst2))/27.)/pow4(Mst2))) - (135*(3*(59 + 148*lmMsq - 184*
        lmMst1 + 36*lmMst2)*pow2(Mst1)*pow2(Mst2) + (13 - 12*lmMsq + 12*lmMst2)
        *pow4(Mst2)))/pow4(Msq))) - (3*pow2(MuSUSY)*(3*pow2(s2t)*((6415 - 60*
        lmMsq*(-37 + 18*lmMst1) - 1314*lmMst2 + 10*lmMst1*(-97 + 66*lmMst2) +
        540*pow2(lmMsq) - 234*pow2(lmMst1) + 54*pow2(lmMst2))/pow2(Mst1) + (20*
        (160 - 213*lmMst1 + 18*lmMsq*(8 + lmMst1 - lmMst2) + 69*lmMst2 - 9*
        pow2(lmMst1) + 9*pow2(lmMst2)) + (4*(-90*(25 - 11*lmMst1 + 2*lmMsq*(-5
        + 9*lmMst1 - 9*lmMst2) + 21*lmMst2 - 9*pow2(lmMst1) + 9*pow2(lmMst2))*
        pow2(Mst1) + pow2(Msq)*(265 + 36*B4 + 30*DN - 1935*lmMsq + lmMst1*(2178
        + 270*lmMsq - 270*pow2(lmMsq)) + 135*pow2(lmMsq) + 3*lmMst2*(353 - 180*
        lmMsq + 86*lmMst1 + 90*pow2(lmMsq) - 79*pow2(lmMst1)) + 30*(-13 + 9*
        lmMsq)*pow2(lmMst1) - 3*(-372 + 90*lmMsq + 49*lmMst1)*pow2(lmMst2) -
        65*pow3(lmMst1) + 449*pow3(lmMst2))))/pow2(Mst2))/pow2(Msq) + 108*pow2(
        Mst1)*((5*(2117 - 1614*lmMst1 + 1074*lmMst2 - 180*pow2(lmMst1) + 180*(
        lmMsq*(3 + 2*lmMst1 - 2*lmMst2) + pow2(lmMst2))))/(216.*pow4(Msq)) + (
        129.29128086419752 - (52*B4)/9. - (22*DN)/9. - (380*lmMsq)/3. - lmMst1*
        (150.47685185185185 - 160*lmMsq + 10*pow2(lmMsq)) + (5*(-57 + 8*lmMsq)*
        pow2(lmMst1))/4. + lmMst2*(324.6990740740741 - 160*lmMsq - (3383*
        lmMst1)/18. + 10*pow2(lmMsq) + (833*pow2(lmMst1))/9.) - ((-8947 + 360*
        lmMsq + 4868*lmMst1)*pow2(lmMst2))/36. - (721*pow3(lmMst1))/27. + (
        1873*pow3(lmMst2))/27.)/pow4(Mst2))) + (8*Mt*(4*Mt*(5483 - 18*B4 - 9*DN
        + 418*lmMst1 + 102*pow2(lmMst1) + 2*lmMst2*(739 - 42*lmMst1 + 9*pow2(
        lmMst1)) + 54*(5 + lmMst1)*pow2(lmMst2) + 24*lmMt*(17 + 17*lmMst2 -
        lmMst1*(17 + 6*lmMst2) + 3*pow2(lmMst1) + 3*pow2(lmMst2)) - 78*pow3(
        lmMst1) + 6*pow3(lmMst2)) - (3*s2t*(((5875 + 1020*lmMsq - 478*lmMst1 -
        610*lmMst2 - 576*pow2(lmMst1))*pow2(Mst2) + 2*pow2(Mst1)*(6129 + 684*B4
        - 18*DN - 3630*lmMsq + 6*(-251 + 600*lmMsq)*lmMst1 + lmMst2*(6914 -
        1122*lmMst1 + 720*lmMsq*(-5 + 3*lmMst1) - 1179*pow2(lmMst1)) - 3*(997 +
        360*lmMsq)*pow2(lmMst1) - 9*(-393 + 120*lmMsq + 125*lmMst1)*pow2(
        lmMst2) + 1353*pow3(lmMst1) + 951*pow3(lmMst2)))*pow4(Msq) + 5*(983 -
        159*lmMst1 + 180*lmMsq*(-2 + lmMst1 - lmMst2) + 519*lmMst2 - 90*pow2(
        lmMst1) + 90*pow2(lmMst2))*pow2(Mst2)*pow4(Mst1) + 20*pow2(Msq)*((311 -
        66*lmMst1 + 36*lmMsq*(-2 + lmMst1 - lmMst2) + 138*lmMst2 - 18*pow2(
        lmMst1) + 18*pow2(lmMst2))*pow2(Mst1)*pow2(Mst2) + (29 - 6*lmMst1 + 36*
        lmMsq*(-2 + 3*lmMst1 - 3*lmMst2) + 78*lmMst2 + 360*lmMst1*lmMst2 - 234*
        pow2(lmMst1) - 126*pow2(lmMst2))*pow4(Mst1)) - 90*pow2(Mst1)*pow4(Mst2)
        ))/(Mst1*pow4(Msq))))/pow4(Mst2)))/pow2(Sbeta)) + (-75*pow4(Mst2)*(60*(
        2153 - 1614*lmMst1 + 1074*lmMst2 - 180*pow2(lmMst1) + 180*(lmMsq*(3 +
        2*lmMst1 - 2*lmMst2) + pow2(lmMst2)))*pow4(Mst1)*pow4(Mst2) + pow4(Msq)
        *(-48*pow2(Mst1)*pow2(Mst2)*(5885 - 72*B4 - 60*DN + 6090*lmMsq + 270*
        pow2(lmMsq) + 2*lmMst1*(-2663 - 810*lmMsq + 270*pow2(lmMsq)) - 6*
        lmMst2*(572 - 180*lmMsq - 24*lmMst1 + 90*pow2(lmMsq) - 79*pow2(lmMst1))
        - 6*(-91 + 90*lmMsq)*pow2(lmMst1) + 6*(-363 + 90*lmMsq + 49*lmMst1)*
        pow2(lmMst2) + 130*pow3(lmMst1) - 898*pow3(lmMst2)) + (438203 - 21888*
        B4 - 12096*DN + 96480*lmMsq - 12960*pow2(lmMsq) + 108*lmMst1*(-7699 +
        3120*lmMsq + 240*pow2(lmMsq)) - 12*lmMst2*(-50563 + 25920*lmMsq +
        43404*lmMst1 + 2160*pow2(lmMsq) - 23784*pow2(lmMst1)) - 72*(1603 + 360*
        lmMsq)*pow2(lmMst1) + 72*(5989 + 360*lmMsq - 4476*lmMst1)*pow2(lmMst2)
        - 56736*pow3(lmMst1) + 93600*pow3(lmMst2))*pow4(Mst1) + 24*(6415 + 60*
        lmMsq*(37 - 18*lmMst1) - 970*lmMst1 + 6*(-219 + 110*lmMst1)*lmMst2 +
        540*pow2(lmMsq) - 234*pow2(lmMst1) + 54*pow2(lmMst2))*pow4(Mst2)) -
        480*pow2(Msq)*(2*(385 - 312*lmMst1 + 18*lmMsq*(3 + 10*lmMst1 - 10*
        lmMst2) + 258*lmMst2 - 90*pow2(lmMst1) + 90*pow2(lmMst2))*pow2(Mst2)*
        pow4(Mst1) - (160 - 213*lmMst1 + 18*lmMsq*(8 + lmMst1 - lmMst2) + 69*
        lmMst2 - 9*pow2(lmMst1) + 9*pow2(lmMst2))*pow2(Mst1)*pow4(Mst2)))*pow4(
        s2t) + 400*Mst1*Mt*pow2(Mst2)*pow3(s2t)*(pow4(Msq)*(144*pow2(Mst1)*
        pow2(Mst2)*(254 + 684*B4 - 18*DN - 4650*lmMsq + 4*(-257 + 900*lmMsq)*
        lmMst1 + 3*lmMst2*(2508 - 374*lmMst1 + 240*lmMsq*(-5 + 3*lmMst1) - 393*
        pow2(lmMst1)) - 15*(161 + 72*lmMsq)*pow2(lmMst1) - 9*(-393 + 120*lmMsq
        + 125*lmMst1)*pow2(lmMst2) + 1353*pow3(lmMst1) + 951*pow3(lmMst2)) + (
        3124381 + 25920*lmMsq*(4 + 39*lmMst1 - 39*lmMst2) - 1047828*lmMst2 +
        72*(-8423 + 2196*lmMst2)*pow2(lmMst1) + 762120*pow2(lmMst2) - 12*
        lmMst1*(-59671 + 9516*lmMst2 + 13176*pow2(lmMst2)) - 52704*pow3(lmMst1)
        + 52704*pow3(lmMst2))*pow4(Mst1) + 72*(5875 + 1020*lmMsq - 478*lmMst1 -
        610*lmMst2 - 576*pow2(lmMst1))*pow4(Mst2)) - 1440*pow2(Msq)*((593 - 36*
        lmMsq*(2 + lmMst1 - lmMst2) - 18*lmMst1*(7 + 20*lmMst2) + 198*(lmMst2 +
        pow2(lmMst1)) + 162*pow2(lmMst2))*pow2(Mst2)*pow4(Mst1) - (311 - 66*
        lmMst1 + 36*lmMsq*(-2 + lmMst1 - lmMst2) + 138*lmMst2 - 18*pow2(lmMst1)
        + 18*pow2(lmMst2))*pow2(Mst1)*pow4(Mst2)) + 360*((1019 - 159*lmMst1 +
        180*lmMsq*(-2 + lmMst1 - lmMst2) + 519*lmMst2 - 90*pow2(lmMst1) + 90*
        pow2(lmMst2))*pow4(Mst1)*pow4(Mst2) - 18*pow2(Mst1)*pow6(Mst2))) + 64*
        pow3(Mt)*(-100*Mst1*s2t*(pow2(Mst1)*(180*pow2(Msq)*((311 - 66*lmMst1 +
        36*lmMsq*(-2 + lmMst1 - lmMst2) + 138*lmMst2 - 18*pow2(lmMst1) + 18*
        pow2(lmMst2))*pow2(Mst2)*pow2(MuSUSY) + pow2(Mst1)*(6*(-37 + 30*lmMst1
        - 30*lmMst2)*(1 + 2*lmMst2 - 2*lmMt)*pow2(Mst2) + (29 - 6*lmMst1 + 36*
        lmMsq*(-2 + 3*lmMst1 - 3*lmMst2) + 78*lmMst2 + 360*lmMst1*lmMst2 - 234*
        pow2(lmMst1) - 126*pow2(lmMst2))*pow2(MuSUSY)) + 2*(203 + 12*lmMst1*(-
        11 + 3*lmMsq - 3*lmMt) + 84*lmMt + lmMst2*(48 - 36*lmMsq + 36*lmMt))*
        pow4(Mst2)) - 45*pow2(Mst2)*(18*pow2(Mst2)*pow2(MuSUSY) - pow2(Mst1)*(
        2*(742 + 105*lmMst2 + 3*(89 + 60*lmMst2)*lmMt - 72*lmMst1*(5 + 3*lmMt)
        + 12*lmMsq*(-1 + 18*lmMst1 - 15*lmMst2 + 3*lmMt) - 36*pow2(lmMsq))*
        pow2(Mst2) + (983 - 159*lmMst1 + 180*lmMsq*(-2 + lmMst1 - lmMst2) +
        519*lmMst2 - 90*pow2(lmMst1) + 90*pow2(lmMst2))*pow2(MuSUSY)) + (224 -
        24*lmMsq + lmMst2*(102 - 72*lmMt) - 78*lmMt + 72*lmMsq*(lmMst2 + lmMt)
        - 72*pow2(lmMsq))*pow4(Mst2))) + pow4(Msq)*(9*(5875 + 1020*lmMsq - 478*
        lmMst1 - 610*lmMst2 - 576*pow2(lmMst1))*pow2(Mst2)*pow2(MuSUSY) + 2*(
        348947 + 544227*lmMst2 - 9*(27947 + 7566*lmMst2 - 10518*lmMt)*pow2(
        lmMst1) + 392355*pow2(lmMst2) - 1620*lmMsq*(65 + lmMst2*(224 - 60*lmMt)
        - 30*lmMst1*(5 + 2*lmMst2 - 2*lmMt) - 74*lmMt + 60*pow2(lmMst2)) - 18*
        lmMt*(8710 + 7595*lmMst2 + 1461*pow2(lmMst2)) - 2592*(9 + 10*lmMst2)*
        pow2(lmMt) - 3*lmMst1*(70397 + 6*(373 + 3798*lmMst2)*lmMt - 5292*pow2(
        lmMst2) - 8640*(lmMst2 + pow2(lmMt))) - 8856*pow3(lmMst1) + 61074*pow3(
        lmMst2))*pow4(Mst1) + 6*(pow2(Mst1)*(3*pow2(MuSUSY)*(6129 + 684*B4 -
        18*DN - 3630*lmMsq + 6*(-251 + 600*lmMsq)*lmMst1 + lmMst2*(6914 - 1122*
        lmMst1 + 720*lmMsq*(-5 + 3*lmMst1) - 1179*pow2(lmMst1)) - 3*(997 + 360*
        lmMsq)*pow2(lmMst1) - 9*(-393 + 120*lmMsq + 125*lmMst1)*pow2(lmMst2) +
        1353*pow3(lmMst1) + 951*pow3(lmMst2)) + pow2(Mst2)*(10810 + 28427*
        lmMst2 + 540*lmMsq*(-13 + 6*lmMst1 - 6*lmMst2)*(1 + 2*lmMst2 - 2*lmMt)
        - 3*(239 + 1134*lmMst2 - 2502*lmMt)*pow2(lmMst1) + 19413*pow2(lmMst2) -
        6*lmMt*(2272 + 907*lmMst2 + 189*pow2(lmMst2)) - 864*(5 + 3*lmMst2)*
        pow2(lmMt) + lmMst1*(4789 + lmMst2*(5208 - 6372*lmMt) - 12414*lmMt -
        324*pow2(lmMst2) + 2592*pow2(lmMt)) - 1368*pow3(lmMst1) + 5094*pow3(
        lmMst2))) + (29825 - 777*lmMst2 + lmMst1*(-9275 + 54*lmMst2 - 1578*
        lmMt) - 6*(-262 + 9*lmMst2)*lmMt + 360*lmMsq*(14 - 3*lmMst1 + 3*lmMt) +
        570*pow2(lmMst1) - 720*pow2(lmMt))*pow4(Mst2)))) + Mt*(2*pow4(Msq)*(4*
        pow2(Mst1)*(pow2(Mst2)*(2829116 - 165375*lmMst2 + 30375*lmMsq*(3 + 2*
        lmMt) - 15*lmMst1*(8071 + 4395*lmMt + 45*lmMst2*(49 + 4*lmMt)) - 30375*
        pow2(lmMsq) + 1350*(-70 + lmMst2 - 31*lmMt)*pow2(lmMst1) - 79425*pow2(
        lmMst2) + 45*lmMt*(15327 + 4745*lmMst2 + 990*pow2(lmMst2)) - 32400*
        pow2(lmMt) + 13500*pow3(lmMst1) - 14850*pow3(lmMst2)) + 150*pow2(
        MuSUSY)*(5483 - 18*B4 - 9*DN + 418*lmMst1 + 102*pow2(lmMst1) + 2*
        lmMst2*(739 - 42*lmMst1 + 9*pow2(lmMst1)) + 54*(5 + lmMst1)*pow2(
        lmMst2) + 24*lmMt*(17 + 17*lmMst2 - lmMst1*(17 + 6*lmMst2) + 3*pow2(
        lmMst1) + 3*pow2(lmMst2)) - 78*pow3(lmMst1) + 6*pow3(lmMst2))) + (
        22616341 + 6969450*lmMst2 + 900*(383 - 1500*lmMst2 + 750*lmMt)*pow2(
        lmMst1) + 4066200*pow2(lmMst2) - 81000*lmMsq*(47 + lmMst2*(107 - 30*
        lmMt) - 30*lmMst1*(1 + lmMst2 - lmMt) - 77*lmMt + 30*pow2(lmMst2)) -
        180*lmMt*(-15793 - 21005*lmMst2 + 1530*pow2(lmMst2)) - 129600*(13 + 6*
        lmMst2)*pow2(lmMt) + 90*lmMst1*(3669 + lmMst2*(27470 - 4440*lmMt) -
        96890*lmMt + 1380*pow2(lmMst2) + 8640*pow2(lmMt)) + 225000*pow3(lmMst1)
        + 1000800*pow3(lmMst2))*pow4(Mst1) - 6*(176669 + 3168*lmMst2 - 30*
        lmMst1*(3644 + 144*lmMst2 - 975*lmMt) + 30*(-203 + 9*lmMst2)*lmMt +
        180*lmMsq*(376 - 240*lmMst1 + 15*lmMt) + 20250*pow2(lmMsq) + 6975*pow2(
        lmMst1) + 2025*pow2(lmMst2) - 16200*pow2(lmMt))*pow4(Mst2)) + 360*pow2(
        Msq)*(750*(11 - 6*lmMsq + 6*lmMt)*pow2(Mst2)*pow4(Mst1) + (2648 + 3480*
        lmMst1 + 75*lmMst2 + (345 - 450*(lmMst1 + lmMst2))*lmMt + 150*lmMsq*(-
        26 + 3*lmMst1 + 3*lmMst2 + 6*lmMt) - 900*pow2(lmMsq))*pow2(Mst1)*pow4(
        Mst2)) + 45*((36871 + 8250*lmMst2 + 360*lmMst1*(91 - 30*lmMt) + (7890 -
        9000*lmMst2)*lmMt + 300*lmMsq*(-163 + 36*lmMst1 + 30*lmMst2 + 66*lmMt)
        - 19800*pow2(lmMsq))*pow4(Mst1)*pow4(Mst2) + 50*(115 - 48*lmMsq +
        lmMst2*(69 - 36*lmMt) - 21*lmMt + 36*lmMsq*(lmMst2 + lmMt) - 36*pow2(
        lmMsq))*pow2(Mst1)*pow6(Mst2)))))/(pow2(Mst1)*pow4(Msq)*pow4(Mst2)) + (
        7200*Mt*MuSUSY*((45*Mt*MuSUSY*Tbeta*pow2(Mst2)*pow2(s2t))/(pow2(Sbeta)*
        pow4(Msq)) - (40*Mst1*(185 + 66*lmMst2 - 6*(13 + 12*lmMst2)*lmMt + 12*
        lmMsq*(1 + 6*(lmMst2 + lmMt)) - 72*pow2(lmMsq))*pow3(Mt))/pow4(Msq) -
        108*Mt*pow2(s2t)*((5875 + 1020*lmMsq - 478*lmMst1 - 610*lmMst2 - 576*
        pow2(lmMst1))/(9.*Mst1) + Mst1*((20*(311 - 66*lmMst1 + 36*lmMsq*(-2 +
        lmMst1 - lmMst2) + 138*lmMst2 - 18*pow2(lmMst1) + 18*pow2(lmMst2)))/(9.
        *pow2(Msq)) + (709.2222222222222 + 152*B4 - 4*DN - 920*lmMsq + (-
        281.55555555555554 + 800*lmMsq)*lmMst1 + (2*lmMst2*(7219 - 1122*lmMst1
        + 720*lmMsq*(-5 + 3*lmMst1) - 1179*pow2(lmMst1)))/9. - (2*(901 + 360*
        lmMsq)*pow2(lmMst1))/3. + (786 - 240*lmMsq - 250*lmMst1)*pow2(lmMst2) -
        (40*(47 + 2*(5 + 6*lmMsq)*lmMst2 - 2*lmMst1*(5 + 6*lmMsq + 30*lmMst2) +
        36*pow2(lmMst1) + 24*pow2(lmMst2))*pow2(Mst1))/(3.*pow2(Msq)) + (902*
        pow3(lmMst1))/3. + (634*pow3(lmMst2))/3.)/pow2(Mst2) - (10*pow2(Mst2))/
        pow4(Msq)) + pow3(Mst1)*((5*(1001 - 159*lmMst1 + 180*lmMsq*(-2 + lmMst1
        - lmMst2) + 519*lmMst2 - 90*pow2(lmMst1) + 90*pow2(lmMst2)))/(9.*pow4(
        Msq)) + (5530.797839506173 + 152*B4 - 4*DN - 760*lmMsq + (
        823.4629629629629 + 2360*lmMsq)*lmMst1 - (1536.5555555555557 + 240*
        lmMsq)*pow2(lmMst1) - lmMst2*(12.796296296296296 + 40*lmMsq*(59 - 12*
        lmMst1) + (3830*lmMst1)/9. + 18*pow2(lmMst1)) + (1962.111111111111 -
        240*lmMsq - 494*lmMst1)*pow2(lmMst2) + (658*pow3(lmMst1))/3. + (878*
        pow3(lmMst2))/3.)/pow4(Mst2))) + (3*pow2(Mst2)*pow3(s2t)*(pow4(Msq)*(-
        24*pow2(Mst1)*pow2(Mst2)*(5355 - 144*B4 - 120*DN + 9960*lmMsq + 2*
        lmMst1*(-4841 - 1080*lmMsq + 540*pow2(lmMsq)) - 6*lmMst2*(925 - 360*
        lmMsq + 62*lmMst1 + 180*pow2(lmMsq) - 158*pow2(lmMst1)) - 6*(-221 +
        180*lmMsq)*pow2(lmMst1) + 6*(-735 + 180*lmMsq + 98*lmMst1)*pow2(lmMst2)
        + 260*pow3(lmMst1) - 1796*pow3(lmMst2)) + (309683 - 18432*B4 - 9216*DN
        - 142560*lmMsq + 12*(-49927 + 32400*lmMsq)*lmMst1 - 12960*pow2(lmMsq) -
        12*lmMst2*(-61663 + 30240*lmMsq + 42660*lmMst1 - 21888*pow2(lmMst1)) -
        147240*pow2(lmMst1) - 72*(-7459 + 4672*lmMst1)*pow2(lmMst2) - 62976*
        pow3(lmMst1) + 136704*pow3(lmMst2))*pow4(Mst1) + 24*(6415 + 60*lmMsq*(
        37 - 18*lmMst1) - 970*lmMst1 + 6*(-219 + 110*lmMst1)*lmMst2 + 540*pow2(
        lmMsq) - 234*pow2(lmMst1) + 54*pow2(lmMst2))*pow4(Mst2)) - 480*pow2(
        Msq)*((610 - 411*lmMst1 + 18*lmMsq*(-2 + 19*lmMst1 - 19*lmMst2) + 447*
        lmMst2 - 171*pow2(lmMst1) + 171*pow2(lmMst2))*pow2(Mst2)*pow4(Mst1) - (
        160 - 213*lmMst1 + 18*lmMsq*(8 + lmMst1 - lmMst2) + 69*lmMst2 - 9*pow2(
        lmMst1) + 9*pow2(lmMst2))*pow2(Mst1)*pow4(Mst2)) + 60*((2135 - 1614*
        lmMst1 + 1074*lmMst2 - 180*pow2(lmMst1) + 180*(lmMsq*(3 + 2*lmMst1 - 2*
        lmMst2) + pow2(lmMst2)))*pow4(Mst1)*pow4(Mst2) - 18*pow2(Mst1)*pow6(
        Mst2))) + 16*pow2(Mt)*(24*Mst1*Mt*(((30251 - 804*lmMst2 + 2*lmMst1*(-
        5032 + 27*lmMst2 - 789*lmMt) + (852 - 54*lmMst2)*lmMt + 180*lmMsq*(31 -
        6*lmMst1 + 6*lmMt) + 570*pow2(lmMst1) - 720*pow2(lmMt))*pow2(Mst2) +
        pow2(Mst1)*(32085 + 23606*lmMst2 + (3606 - 3402*lmMst2 + 7506*lmMt)*
        pow2(lmMst1) + 180*lmMsq*(31 + 6*lmMst1*(-1 + 6*lmMst2 - 6*lmMt) + 84*
        lmMt + 6*lmMst2*(-13 + 6*lmMt) - 36*pow2(lmMst2)) + 18846*pow2(lmMst2)
        - 6*lmMt*(2850 + 1348*lmMst2 + 189*pow2(lmMst2)) - 2*lmMst1*(5093 +
        5700*lmMt + 6*lmMst2*(-173 + 531*lmMt) + 162*pow2(lmMst2) - 1296*pow2(
        lmMt)) - 144*(35 + 18*lmMst2)*pow2(lmMt) - 1368*pow3(lmMst1) + 5094*
        pow3(lmMst2)))*pow4(Msq) + 135*(87 + 4*lmMst1*(-13 + 6*lmMsq - 6*lmMt)
        + 34*lmMt + 6*lmMst2*(3 - 4*lmMsq + 4*lmMt))*pow2(Mst2)*pow4(Mst1) +
        60*pow2(Msq)*((245 + 6*lmMst1*(-25 + 6*lmMsq - 6*lmMt) + 84*lmMt +
        lmMst2*(66 - 36*lmMsq + 36*lmMt))*pow2(Mst1)*pow2(Mst2) + (245 + 6*
        lmMst1*(-25 + 6*lmMsq + 30*lmMst2 - 36*lmMt) - 12*lmMst2*(13 + 3*lmMsq
        - 18*lmMt) + 306*lmMt - 180*pow2(lmMst2))*pow4(Mst1))) + s2t*(2*pow4(
        Msq)*(pow2(Mst1)*pow2(Mst2)*(416017 - 864*B4 - 432*DN - 14580*lmMsq +
        45858*lmMst1 + 4860*pow2(lmMsq) - 6*lmMst2*(-21571 + 1620*lmMsq + 618*
        lmMst1 - 432*pow2(lmMst1)) + 7902*pow2(lmMst1) + 18*(1735 + 336*lmMst1)
        *pow2(lmMst2) + 1728*lmMt*(17 + 17*lmMst2 - lmMst1*(17 + 6*lmMst2) + 3*
        pow2(lmMst1) + 3*pow2(lmMst2)) - 7200*pow3(lmMst1) - 1440*pow3(lmMst2))
        + 2*(226991 - 432*B4 - 216*DN - 76140*lmMsq + 24*(-1742 + 2025*lmMsq)*
        lmMst1 - 163881*pow2(lmMst1) + 6*lmMst2*(19352 - 8100*lmMsq + 60387*
        lmMst1 + 1512*pow2(lmMst1)) + 9*(-21473 + 9840*lmMst1)*pow2(lmMst2) +
        864*lmMt*(71 + 118*lmMst2 - 2*lmMst1*(59 + 39*lmMst2) + 39*pow2(lmMst1)
        + 39*pow2(lmMst2)) - 40752*pow3(lmMst1) - 56880*pow3(lmMst2))*pow4(
        Mst1) - 18*(4717 - 540*lmMsq*(-2 + lmMst1) - 252*lmMst2 + 2*lmMst1*(-
        661 + 69*lmMst2) + 270*pow2(lmMsq) + 75*pow2(lmMst1) + 27*pow2(lmMst2))
        *pow4(Mst2)) - 1080*pow2(Msq)*(10*(11 - 6*lmMsq + 6*lmMst2)*pow2(Mst2)*
        pow4(Mst1) + (5 + 60*lmMsq - 66*lmMst1 + 6*lmMst2)*pow2(Mst1)*pow4(
        Mst2)) - 135*(2*(95 + 216*lmMsq - 276*lmMst1 + 60*lmMst2)*pow4(Mst1)*
        pow4(Mst2) + (13 - 12*lmMsq + 12*lmMst2)*pow2(Mst1)*pow6(Mst2)))))/(72.
        *pow2(Mst1)*pow4(Msq)*pow4(Mst2))))/Tbeta))/777600.)) + threeLoopFlag*
        pow2(Al4p)*(((10501 - 4920*lmMsq + 18*lmMst2*(149 + 50*lmMsq - 19*
        lmMst1*(-2 + lmMst2) + lmMst2 - 30*pow2(lmMsq)) + 1980*pow2(lmMsq) - 4*
        lmMst1*(196 + 270*lmMsq + 135*pow2(lmMsq)) + 1914*pow2(lmMst1) - 6*
        lmMt*(77 + 2*lmMst1*(221 - 57*lmMst2) + 594*lmMst2 + 90*lmMsq*(7 - 2*(
        lmMst1 + lmMst2)) - 129*pow2(lmMst1) + 87*pow2(lmMst2)) - 108*(-28 +
        10*lmMsq + 29*lmMst1 + 17*lmMst2)*pow2(lmMt) + (2*Dmglst1*(-2215 - 254*
        lmMst2 + 2322*lmMt - 264*lmMst2*lmMt - 60*lmMsq*(5 + 42*lmMst1 + 12*
        lmMt) - 4*lmMst1*(206 + 15*lmMst2 + 246*lmMt) + 1620*pow2(lmMsq) +
        3018*pow2(lmMst1) + 162*pow2(lmMst2) + 1872*pow2(lmMt)))/Mst1 + (360*
        Dmglst1*Mst1*(-27 - lmMst2 + 3*lmMt + 6*lmMst2*lmMt - 6*lmMsq*(-2 +
        lmMst1 + lmMst2 + 2*lmMt) + 2*lmMst1*(-7 + 3*lmMt) + 12*pow2(lmMsq)))/
        pow2(Msq) + (108*(1 - 2*lmMst2)*pow2(Mst2))/pow2(Mst1) + 360*pow3(
        lmMsq) - 750*pow3(lmMst1) - (4*Mst1*(Mst1*(2459588 - 813345*lmMst2 +
        151875*lmMsq*(3 + 2*lmMt) - 151875*pow2(lmMsq) + 225*(1339 + 30*lmMst2
        - 900*lmMt)*pow2(lmMst1) - 370350*pow2(lmMst2) + 3375*lmMt*(66 + 303*
        lmMst2 + 68*pow2(lmMst2)) + 15*lmMst1*(17548 - 71775*lmMt - 15*lmMst2*(
        773 + 120*lmMt) + 450*pow2(lmMst2)) - 162000*pow2(lmMt) + 65250*pow3(
        lmMst1) - 78750*pow3(lmMst2)) + 125*Dmglst1*(64093 - 12414*lmMst2 +
        2430*lmMsq*(3 + 2*lmMt) - 6*lmMst1*(988 + 2367*lmMt + lmMst2*(393 + 36*
        lmMt)) - 2430*pow2(lmMsq) + 36*(38 + 3*lmMst2 - 93*lmMt)*pow2(lmMst1) -
        6354*pow2(lmMst2) + 18*lmMt*(1190 + 933*lmMst2 + 198*pow2(lmMst2)) -
        2592*pow2(lmMt) + 1080*pow3(lmMst1) - 1188*pow3(lmMst2))))/(375.*pow2(
        Mst2)) - 204*pow3(lmMst2) + 4968*pow3(lmMt) + (15*Mst1*(2*Dmglst1 +
        Mst1)*(-115 + 21*lmMt - 12*lmMsq*(-4 + 3*lmMst2 + 3*lmMt) + lmMst2*(-69
        + 36*lmMt) + 36*pow2(lmMsq))*pow2(Mst2))/pow4(Msq) + (pow3(Mst1)*((15*(
        171500*Dmglst1*(-2095 - 330*lmMst2 - 114*lmMt + 360*lmMst2*lmMt + 48*
        lmMst1*(-25 + 12*lmMt) - 12*lmMsq*(-137 + 48*lmMst1 + 30*lmMst2 + 78*
        lmMt) + 936*pow2(lmMsq)) + Mst1*(-160842737 - 14148750*lmMst2 +
        10547250*lmMt + 15435000*lmMst2*lmMt + 840*lmMst1*(-106067 + 51450*
        lmMt) - 420*lmMsq*(-220709 + 82740*lmMst1 + 36750*lmMst2 + 139650*lmMt)
        + 54419400*pow2(lmMsq) - 4233600*pow2(lmMst1))))/pow4(Msq) - (4*(
        257250*Dmglst1*(99985 + 34432*lmMst2 + (3156 - 5832*lmMst2 + 1296*lmMt)
        *pow2(lmMst1) + 17388*pow2(lmMst2) + 24*lmMt*(-191 + 871*lmMst2 + 6*
        pow2(lmMst2)) - 1080*lmMsq*(9 + lmMst2*(29 - 10*lmMt) - 10*lmMst1*(1 +
        lmMst2 - lmMt) - 19*lmMt + 10*pow2(lmMst2)) - 1728*(3 + 2*lmMst2)*pow2(
        lmMt) + 8*lmMst1*(517 + lmMst2*(168 - 180*lmMt) - 4557*lmMt + 27*pow2(
        lmMst2) + 432*pow2(lmMt)) + 1512*pow3(lmMst1) + 4104*pow3(lmMst2)) +
        Mst1*(6336529126 + 2073605310*lmMst2 - 44100*(-8128 + 8295*lmMst2 +
        1680*lmMt)*pow2(lmMst1) + 976616550*pow2(lmMst2) - 138915000*lmMsq*(2 +
        lmMst2*(12 - 5*lmMt) - 5*lmMst1*(1 + lmMst2 - lmMt) - 7*lmMt + 5*pow2(
        lmMst2)) + 1157625*lmMt*(-899 + 1578*lmMst2 + 128*pow2(lmMst2)) -
        222264000*(1 + lmMst2)*pow2(lmMt) + 315*lmMst1*(1895351 - 7915950*lmMt
        - 70*lmMst2*(16867 + 3360*lmMt) - 14700*pow2(lmMst2) + 705600*pow2(
        lmMt)) + 146632500*pow3(lmMst1) + 223807500*pow3(lmMst2))))/pow4(Mst2))
        )/514500. - (24*Mst1*pow2(MuSUSY)*(2*Dmglst1*(735 - 6*B4 - 3*DN - 78*
        lmMst1 + 6*pow2(lmMst1) + 6*lmMst2*(85 - 6*lmMst1 + pow2(lmMst1)) + 18*
        (7 + lmMst1)*pow2(lmMst2) + 24*lmMt*(2 + 3*lmMst2 - lmMst1*(3 + 2*
        lmMst2) + pow2(lmMst1) + pow2(lmMst2)) - 26*pow3(lmMst1) + 2*pow3(
        lmMst2)) + Mst1*(417 - 6*B4 - 3*DN - 162*lmMst1 + 6*lmMst2*(87 - 8*
        lmMst1 + pow2(lmMst1)) + 18*(8 + lmMst1)*pow2(lmMst2) + 24*lmMt*pow2(1
        - lmMst1 + lmMst2) - 26*pow3(lmMst1) + 2*pow3(lmMst2))))/pow4(Mst2) - (
        3*(13597579 + 210*lmMst2*(43831 - 14700*lmMt) - 2315250*lmMt + 420*
        lmMsq*(-16403 + 630*lmMst2 + 7350*lmMt) - 1675800*pow2(lmMsq) +
        1411200*pow2(lmMst2))*pow4(Mst2))/(34300.*pow4(Msq)) + (4*(3*(-12339 -
        4910*lmMst1 - 375*lmMst2 + 4125*lmMt + 3750*lmMst1*lmMt + 2250*lmMst2*
        lmMt - 10*lmMsq*(-116 + 285*lmMst1 + 225*lmMst2 + 600*lmMt) + 5550*
        pow2(lmMsq) - 450*pow2(lmMst1))*pow2(Mst1)*pow2(Mst2) + 7500*Dmglst1*(-
        11 + 6*lmMsq - 6*lmMt)*pow3(Mst1) + 1875*(-11 + 6*lmMsq - 6*lmMt)*pow4(
        Mst1) + (-13517 - 7230*lmMst2 + 6000*lmMt + 4500*lmMst2*lmMt - 30*
        lmMsq*(-41 + 60*lmMst2 + 150*lmMt) + 3150*pow2(lmMsq) - 1350*pow2(
        lmMst2))*pow4(Mst2)))/(25.*pow2(Msq)*pow2(Mst2)))*pow4(Mt))/81. - pow4(
        s2t)*(-(Mst1*pow2(Mst2)*((-5876588*Mst1)/225. + 60*B4*Mst1 + 90*DN*Mst1
        + 3975*lmMsq*Mst1 - 1575*Mst1*pow2(lmMsq) - (3*lmMst1*Mst1*(34697 -
        8250*lmMsq + 2250*pow2(lmMsq)))/5. + 18*(-83 + 25*lmMsq)*Mst1*pow2(
        lmMst1) + (3*lmMst2*Mst1*(3322 + 3000*lmMsq*(-1 + lmMst1) - 25870*
        lmMst1 + 2250*pow2(lmMsq) + 3025*pow2(lmMst1)))/5. - 3*(-3372 + 750*
        lmMsq + 2365*lmMst1)*Mst1*pow2(lmMst2) + (2*(Mst1*(-2068 + 646*lmMst1 -
        991*lmMst2 + 225*lmMst1*lmMst2 + lmMsq*(345 - 615*lmMst1 + 615*lmMst2)
        + 195*pow2(lmMst1) - 420*pow2(lmMst2)) + 50*Dmglst1*(-137 + 92*lmMst1 -
        98*lmMst2 + lmMsq*(6 - 60*lmMst1 + 60*lmMst2) + 30*pow2(lmMst1) - 30*
        pow2(lmMst2)))*pow2(Mst1))/pow2(Msq) - 205*Mst1*pow3(lmMst1) + 5485*
        Mst1*pow3(lmMst2) + 5*Dmglst1*(1151 + 72*B4 + 60*DN - 1890*lmMsq - 270*
        pow2(lmMsq) - 18*lmMst1*(-49 - 90*lmMsq + 30*pow2(lmMsq)) + 6*lmMst2*(
        342 - 180*lmMsq - 24*lmMst1 + 90*pow2(lmMsq) - 79*pow2(lmMst1)) + 6*(-
        347 + 90*lmMsq)*pow2(lmMst1) - 6*(-363 + 90*lmMsq + 49*lmMst1)*pow2(
        lmMst2) - 130*pow3(lmMst1) + 898*pow3(lmMst2))))/540. + Dmglst1*(
        49.72923096707819 + (2*B4)/3. + (5*DN)/9. - (45*lmMsq)/2. + lmMst1*(
        79.81095679012346 - 15*lmMsq - 5*pow2(lmMsq)) + (5*pow2(lmMsq))/2. +
        lmMst2*(10*lmMsq + 5*pow2(lmMsq) + (-31507 + 25692*lmMst1 - 23544*pow2(
        lmMst1))/1296.) + (1.2546296296296295 + 5*lmMsq)*pow2(lmMst1) - (
        2.8564814814814814 + 5*lmMsq - (455*lmMst1)/18.)*pow2(lmMst2) - (73*
        pow3(lmMst1))/54. - (311*pow3(lmMst2))/54.)*pow3(Mst1) + (
        32.221840780308305 + (10*B4)/9. + DN/18. - (275*lmMsq)/72. + lmMst1*(
        3.2322576530612244 + (65*lmMsq)/18. - (35*pow2(lmMsq))/12.) + (5*pow2(
        lmMsq))/12. + lmMst2*((-5*lmMsq*(8 + 3*lmMst1))/9. + (35*pow2(lmMsq))/
        12. + (5211957 + 20946380*lmMst1 - 38602200*pow2(lmMst1))/2.1168e6) - (
        17.322652116402118 - (15*lmMsq)/4.)*pow2(lmMst1) + (8.48290343915344 -
        (25*lmMsq)/12. + (1793*lmMst1)/72.)*pow2(lmMst2) + (95*pow3(lmMst1))/
        216. - (1535*pow3(lmMst2))/216.)*pow4(Mst1) + (2*OepS2*(-150*pow2(Mst1)
        *pow2(Mst2) + 11*pow4(Mst1) - 27*pow4(Mst2)))/729. - (S2*(-18*(197 +
        50*lmMst1 - 50*lmMst2)*pow2(Mst1)*pow2(Mst2) + (281 + 66*lmMst1 - 66*
        lmMst2)*pow4(Mst1) + 81*(11 - 2*lmMst1 + 2*lmMst2)*pow4(Mst2)))/108. -
        (pow4(Mst2)*(41160*(2714 + 30*lmMsq*(-17 + 6*lmMst1 - 6*lmMst2) + 1167*
        lmMst2 + 9*lmMst1*(-73 + 50*lmMst2) - 315*pow2(lmMst1) - 135*pow2(
        lmMst2))*pow2(Msq)*pow3(Mst1) + 8575*Mst1*(53749 + 2592*B4 - 288*DN -
        13320*lmMsq + 5400*pow2(lmMsq) - 6*lmMst1*(3781 - 300*lmMsq + 180*pow2(
        lmMsq)) + 6*lmMst2*(14209 - 3498*lmMst1 + 60*lmMsq*(-35 + 12*lmMst1) +
        180*pow2(lmMsq) - 18*pow2(lmMst1)) - 18*(-229 + 60*lmMsq)*pow2(lmMst1)
        - 18*(-2193 + 180*lmMsq + 442*lmMst1)*pow2(lmMst2) + 396*pow3(lmMst1) +
        7668*pow3(lmMst2))*pow4(Msq) + 51450*Dmglst1*(40*(82 - 93*lmMst1 + 6*
        lmMsq*(4 + 3*lmMst1 - 3*lmMst2) + 69*lmMst2 - 9*pow2(lmMst1) + 9*pow2(
        lmMst2))*pow2(Msq)*pow2(Mst1) + 6*(7 + 60*lmMsq*(5 - 6*lmMst1) - 2*
        lmMst1 + (-226 + 220*lmMst1)*lmMst2 + 180*pow2(lmMsq) + 178*pow2(
        lmMst1) + 18*pow2(lmMst2))*pow4(Msq) + 5*(521 - 346*lmMst1 + 12*lmMsq*(
        1 + 10*lmMst1 - 10*lmMst2) + 334*lmMst2 - 60*pow2(lmMst1) + 60*pow2(
        lmMst2))*pow4(Mst1)) + 3*(12119532 + 420*lmMsq*(-9224 + 6545*lmMst1 -
        5705*lmMst2) + 7553665*lmMst2 + 35*lmMst1*(-105131 + 35280*lmMst2) -
        176400*pow2(lmMsq) - 1991850*pow2(lmMst1) + 580650*pow2(lmMst2))*pow5(
        Mst1)))/(2.22264e7*Mst1*pow4(Msq))) + (Mt*s2t*(-(pow2(Mst2)*pow2(s2t)*(
        180*pow2(Mst1)*(Dmglst1*(1807 - 234*lmMst1 + 12*lmMsq*(-61 + 30*lmMst1
        - 30*lmMst2) + 966*lmMst2 - 180*pow2(lmMst1) + 180*pow2(lmMst2))*pow2(
        Mst1) + Dmglst1*(-91 + 12*lmMsq - 12*lmMst2)*pow2(Mst2) + (-43 + 12*
        lmMsq - 12*lmMst2)*Mst1*pow2(Mst2) + (343 - 18*lmMst1 + 12*lmMsq*(-13 +
        6*lmMst1 - 6*lmMst2) + 174*lmMst2 - 36*pow2(lmMst1) + 36*pow2(lmMst2))*
        pow3(Mst1))*pow4(Mst2) + 2880*pow2(Msq)*pow2(Mst2)*(2*(37 - 8*lmMst1 +
        6*lmMsq*(-2 + lmMst1 - lmMst2) + 20*lmMst2 - 3*pow2(lmMst1) + 3*pow2(
        lmMst2))*pow2(Mst2)*pow3(Mst1) - 3*Mst1*pow4(Mst2) - 3*Dmglst1*(-2*(42
        - 10*lmMst1 + 6*lmMsq*(-2 + lmMst1 - lmMst2) + 22*lmMst2 - 3*pow2(
        lmMst1) + 3*pow2(lmMst2))*pow2(Mst1)*pow2(Mst2) + (107 + 32*lmMst1 -
        12*lmMsq*(2 + lmMst1 - lmMst2) - 8*lmMst2 - 60*lmMst1*lmMst2 + 36*pow2(
        lmMst1) + 24*pow2(lmMst2))*pow4(Mst1) + pow4(Mst2)) - (83 + 52*lmMst1 -
        12*lmMsq*(2 + lmMst1 - lmMst2) - 4*(7 + 9*lmMst1)*lmMst2 + 24*pow2(
        lmMst1) + 12*pow2(lmMst2))*pow5(Mst1)) + pow4(Msq)*(24*Dmglst1*pow2(
        Mst1)*pow2(Mst2)*(54053 + 9432*B4 - 180*DN - 34560*lmMsq + 2160*pow2(
        lmMsq) + 12*lmMst1*(2489 + 810*lmMsq + 90*pow2(lmMsq)) - 54*(-381 +
        260*lmMsq)*pow2(lmMst1) - 12*lmMst2*(-2899 + 2559*lmMst1 - 90*lmMsq*(-
        13 + 24*lmMst1) + 90*pow2(lmMsq) + 984*pow2(lmMst1)) - 18*(-1011 + 660*
        lmMsq + 688*lmMst1)*pow2(lmMst2) + 15648*pow3(lmMst1) + 8544*pow3(
        lmMst2)) + 8*pow2(Mst2)*(132407 + 11880*B4 - 108*DN - 45360*lmMsq +
        6480*pow2(lmMsq) + 60*lmMst1*(2299 - 378*lmMsq + 54*pow2(lmMsq)) - 18*(
        -3377 + 900*lmMsq)*pow2(lmMst1) - 12*lmMst2*(1667 + 1761*lmMst1 - 270*
        lmMsq*(3 + 8*lmMst1) + 270*pow2(lmMsq) + 594*pow2(lmMst1)) - 18*(859 +
        540*lmMsq + 564*lmMst1)*pow2(lmMst2) + 14472*pow3(lmMst1) + 2808*pow3(
        lmMst2))*pow3(Mst1) + Dmglst1*(3058187 + 12960*lmMsq*(43 + 66*lmMst1 -
        66*lmMst2) - 960732*lmMst2 + 72*(-9407 + 6948*lmMst2)*pow2(lmMst1) +
        36*lmMst1*(4415 + 9980*lmMst2 - 10824*pow2(lmMst2)) + 400968*pow2(
        lmMst2) - 203616*pow3(lmMst1) + 93024*pow3(lmMst2))*pow4(Mst1) - 72*
        Dmglst1*(11337 - 408*B4 - 12*DN - 6480*lmMsq + 720*pow2(lmMsq) - 30*
        lmMst1*(-73 - 36*lmMsq + 12*pow2(lmMsq)) + 6*lmMst2*(2331 - 420*lmMsq -
        118*lmMst1 + 60*pow2(lmMsq) - 60*pow2(lmMst1)) + 18*(-13 + 20*lmMsq)*
        pow2(lmMst1) - 6*(-797 + 60*lmMsq + 4*lmMst1)*pow2(lmMst2) - 200*pow3(
        lmMst1) + 584*pow3(lmMst2))*pow4(Mst2) - 72*Mst1*(13317 - 408*B4 - 12*
        DN - 4320*lmMsq + 720*pow2(lmMsq) - 6*lmMst1*(79 - 180*lmMsq + 60*pow2(
        lmMsq)) + 6*lmMst2*(2135 - 420*lmMsq - 118*lmMst1 + 60*pow2(lmMsq) -
        60*pow2(lmMst1)) + 6*(-167 + 60*lmMsq)*pow2(lmMst1) - 6*(-797 + 60*
        lmMsq + 4*lmMst1)*pow2(lmMst2) - 200*pow3(lmMst1) + 584*pow3(lmMst2))*
        pow4(Mst2) + (622151 + 12960*lmMsq*(11 + 10*lmMst1 - 10*lmMst2) -
        135756*lmMst2 + 72*(-4499 + 3060*lmMst2)*pow2(lmMst1) - 185688*pow2(
        lmMst2) - 36*lmMst1*(69 - 14924*lmMst2 + 3048*pow2(lmMst2)) - 110304*
        pow3(lmMst1) - 288*pow3(lmMst2))*pow5(Mst1)))) - 8*pow2(Mt)*(Mst1*(2*
        pow4(Msq)*(9*pow2(Mst2)*pow2(MuSUSY)*(13317 - 408*B4 - 12*DN - 4320*
        lmMsq + 720*pow2(lmMsq) - 6*lmMst1*(79 - 180*lmMsq + 60*pow2(lmMsq)) +
        6*lmMst2*(2135 - 420*lmMsq - 118*lmMst1 + 60*pow2(lmMsq) - 60*pow2(
        lmMst1)) + 6*(-167 + 60*lmMsq)*pow2(lmMst1) - 6*(-797 + 60*lmMsq + 4*
        lmMst1)*pow2(lmMst2) - 200*pow3(lmMst1) + 584*pow3(lmMst2)) + pow2(
        Mst1)*(-12*pow2(Mst2)*(83 + 540*lmMsq*(-1 + 2*lmMst1 - 2*lmMst2)*(1 +
        2*lmMst2 - 2*lmMt) - 18*(163 + 49*lmMst2 - 116*lmMt)*pow2(lmMst1) +
        6282*pow2(lmMst2) - 36*lmMt*(31 + 160*lmMst2 + 22*pow2(lmMst2)) +
        lmMst2*(6591 - 864*pow2(lmMt)) - 3*lmMst1*(1249 - 1536*lmMt + 108*
        lmMst2*(5 + 4*lmMt) + 66*pow2(lmMst2) - 288*pow2(lmMt)) - 594*pow3(
        lmMst1) + 1674*pow3(lmMst2)) + pow2(MuSUSY)*(107299 - 19224*B4 - 108*DN
        - 32400*lmMsq + 6480*pow2(lmMsq) - 24*lmMst1*(6103 - 1755*lmMsq + 405*
        pow2(lmMsq)) + 18*(-4379 + 1260*lmMsq)*pow2(lmMst1) + 12*lmMst2*(20882
        + 699*lmMst1 - 270*lmMsq*(17 + 8*lmMst1) + 810*pow2(lmMsq) + 54*pow2(
        lmMst1)) + 18*(5641 + 180*lmMsq + 540*lmMst1)*pow2(lmMst2) - 18072*
        pow3(lmMst1) + 7704*pow3(lmMst2))) + 4*(20293 + 9*(3199 + 504*lmMst2 -
        852*lmMt)*pow2(lmMst1) - 1620*lmMsq*(7 + 3*lmMst1*(5 + 2*lmMst2 - 2*
        lmMt) + 2*lmMt + lmMst2*(-17 + 6*lmMt) - 6*pow2(lmMst2)) - 41841*pow2(
        lmMst2) + 54*lmMt*(83 + 412*lmMst2 + 82*pow2(lmMst2)) + 18*lmMst1*(4094
        + 533*lmMst2 - 1044*lmMt + 180*lmMst2*lmMt + 42*pow2(lmMst2) - 144*
        pow2(lmMt)) + 18*lmMst2*(-2903 + 144*pow2(lmMt)) + 1620*pow3(lmMst1) -
        6912*pow3(lmMst2))*pow4(Mst1) - 36*(371 + 156*lmMst2 + 90*lmMsq*(11 +
        6*lmMt + 4*lmMst2*(1 + lmMt) - 2*lmMst1*(5 + 2*lmMt)) + 180*(lmMst1 -
        lmMst2)*pow2(lmMsq) + 12*(60 + 5*lmMst2 + 14*lmMt)*pow2(lmMst1) - 96*
        pow2(lmMst2) - 6*lmMt*(217 + 154*lmMst2 + 36*pow2(lmMst2)) - 288*(1 +
        lmMst2)*pow2(lmMt) + 6*lmMst1*(-145 + 26*lmMt + 8*lmMst2*(5 + lmMt) +
        2*pow2(lmMst2) + 48*pow2(lmMt)) - 44*pow3(lmMst1) - 28*pow3(lmMst2))*
        pow4(Mst2)) + 15*pow2(Mst2)*(pow2(Mst1)*pow2(Mst2)*(8*(112 - 12*lmMsq +
        51*lmMst2 - 3*(13 + 12*lmMst2)*lmMt + 36*lmMsq*(lmMst2 + lmMt) - 36*
        pow2(lmMsq))*pow2(Mst2) + (107 - 12*lmMsq + 12*lmMst2)*pow2(MuSUSY)) -
        6*((400 + 42*lmMst2 + 6*lmMt - 144*lmMst1*lmMt + 72*lmMst2*lmMt + 24*
        lmMsq*(-2 + 6*lmMst1 - 3*lmMst2 + 3*lmMt) - 72*pow2(lmMsq))*pow2(Mst2)
        + (134 - 9*lmMst1 + 36*lmMsq*(-2 + lmMst1 - lmMst2) + 81*lmMst2 - 18*
        pow2(lmMst1) + 18*pow2(lmMst2))*pow2(MuSUSY))*pow4(Mst1) + (-11 + 12*
        lmMsq - 12*lmMst2)*pow2(MuSUSY)*pow4(Mst2) + 4*(169 - 48*lmMsq +
        lmMst2*(87 - 36*lmMt) - 39*lmMt + 36*lmMsq*(lmMst2 + lmMt) - 36*pow2(
        lmMsq))*pow6(Mst2)) + 240*pow2(Msq)*(2*pow2(Mst1)*pow2(Mst2)*((27*
        lmMst1*(3 + 2*lmMt) - 3*lmMsq*(7 + 18*lmMst1 - 12*lmMst2 + 6*lmMt) - 2*
        (55 + 6*lmMt + 6*lmMst2*(4 + 3*lmMt)) + 18*pow2(lmMsq))*pow2(Mst2) - 3*
        (34 - 8*lmMst1 + 6*lmMsq*(-2 + lmMst1 - lmMst2) + 20*lmMst2 - 3*pow2(
        lmMst1) + 3*pow2(lmMst2))*pow2(MuSUSY)) + 6*(3*(1 - 3*lmMst1 + 3*
        lmMst2)*(1 + 2*lmMst2 - 2*lmMt)*pow2(Mst2) + 2*(-14 + lmMst1*(21 - 9*
        lmMst2) - 27*lmMst2 + lmMsq*(6 - 9*lmMst1 + 9*lmMst2) + 9*pow2(lmMst1))
        *pow2(MuSUSY))*pow4(Mst1) + 9*pow2(MuSUSY)*pow4(Mst2) + (67 + lmMst2*(
        24 - 36*lmMt) - 66*lmMt + 6*lmMsq*(7 + 6*lmMst2 + 6*lmMt) - 36*pow2(
        lmMsq))*pow6(Mst2))) + Dmglst1*(-2*pow4(Msq)*(-9*pow2(Mst2)*pow2(
        MuSUSY)*(11337 - 408*B4 - 12*DN - 6480*lmMsq + 720*pow2(lmMsq) - 30*
        lmMst1*(-73 - 36*lmMsq + 12*pow2(lmMsq)) + 6*lmMst2*(2331 - 420*lmMsq -
        118*lmMst1 + 60*pow2(lmMsq) - 60*pow2(lmMst1)) + 18*(-13 + 20*lmMsq)*
        pow2(lmMst1) - 6*(-797 + 60*lmMsq + 4*lmMst1)*pow2(lmMst2) - 200*pow3(
        lmMst1) + 584*pow3(lmMst2)) + pow2(Mst1)*(-3*pow2(MuSUSY)*(13969 -
        11880*B4 + 108*DN - 4320*lmMsq + (-16728 - 3240*lmMsq*(1 + lmMsq))*
        lmMst1 + 2160*pow2(lmMsq) + 54*(-407 + 300*lmMsq)*pow2(lmMst1) + 12*
        lmMst2*(4094 + 2205*lmMst1 - 90*lmMsq*(1 + 24*lmMst1) + 270*pow2(lmMsq)
        + 804*pow2(lmMst1)) + 18*(583 + 540*lmMsq + 680*lmMst1)*pow2(lmMst2) -
        16848*pow3(lmMst1) - 5040*pow3(lmMst2)) + 4*pow2(Mst2)*(24721 + 76707*
        lmMst2 + 1620*lmMsq*(-7 + 6*lmMst1 - 6*lmMst2)*(1 + 2*lmMst2 - 2*lmMt)
        - 18*(838 + 525*lmMst2 - 1182*lmMt)*pow2(lmMst1) + 58824*pow2(lmMst2) -
        36*lmMt*(939 + 1045*lmMst2 + 129*pow2(lmMst2)) - 9*lmMst1*(1119 - 1108*
        lmMt + 12*lmMst2*(53 + 154*lmMt) + 138*pow2(lmMst2) - 864*pow2(lmMt)) -
        2592*(2 + 3*lmMst2)*pow2(lmMt) - 4518*pow3(lmMst1) + 15210*pow3(lmMst2)
        )) + 4*(54037 + 273408*lmMst2 - 9*(14497 + 3276*lmMst2 - 4920*lmMt)*
        pow2(lmMst1) + 1620*lmMsq*(5 + 15*lmMst1*(5 + 2*lmMst2 - 2*lmMt) + 22*
        lmMt + lmMst2*(-97 + 30*lmMt) - 30*pow2(lmMst2)) + 207819*pow2(lmMst2)
        - 18*lmMt*(3201 + 5066*lmMst2 + 900*pow2(lmMst2)) - 2592*(2 + 5*lmMst2)
        *pow2(lmMt) + 6*lmMst1*(-41437 + 8286*lmMt - 15*lmMst2*(341 + 312*lmMt)
        + 342*pow2(lmMst2) + 2160*pow2(lmMt)) - 5508*pow3(lmMst1) + 32940*pow3(
        lmMst2))*pow4(Mst1) + 12*(7007 - 290*lmMst2 + 30*lmMsq*(251 + 102*lmMt
        + 36*lmMst2*(1 + lmMt) - 6*lmMst1*(23 + 6*lmMt)) + 540*(lmMst1 -
        lmMst2)*pow2(lmMsq) + 12*(331 + 15*lmMst2 + 42*lmMt)*pow2(lmMst1) -
        288*pow2(lmMst2) - 6*lmMt*(817 + 464*lmMst2 + 108*pow2(lmMst2)) - 864*(
        2 + lmMst2)*pow2(lmMt) + 4*lmMst1*(-2438 - 408*lmMt + 3*lmMst2*(61 +
        12*lmMt) + 9*pow2(lmMst2) + 216*pow2(lmMt)) - 132*pow3(lmMst1) - 84*
        pow3(lmMst2))*pow4(Mst2)) + 15*pow2(Mst2)*(pow2(Mst1)*pow2(Mst2)*(24*(
        112 - 12*lmMsq + 51*lmMst2 - 3*(13 + 12*lmMst2)*lmMt + 36*lmMsq*(lmMst2
        + lmMt) - 36*pow2(lmMsq))*pow2(Mst2) + (251 - 12*lmMsq + 12*lmMst2)*
        pow2(MuSUSY)) - 6*(2*(736 + 105*lmMst2 + 3*(53 + 60*lmMst2)*lmMt + 4*
        lmMsq*(-10 + 66*lmMst1 - 45*lmMst2 + 21*lmMt) - 8*lmMst1*(28 + 33*lmMt)
        - 84*pow2(lmMsq))*pow2(Mst2) + (818 - 117*lmMst1 + 180*lmMsq*(-2 +
        lmMst1 - lmMst2) + 477*lmMst2 - 90*pow2(lmMst1) + 90*pow2(lmMst2))*
        pow2(MuSUSY))*pow4(Mst1) + (-11 + 12*lmMsq - 12*lmMst2)*pow2(MuSUSY)*
        pow4(Mst2) + 4*(169 - 48*lmMsq + lmMst2*(87 - 36*lmMt) - 39*lmMt + 36*
        lmMsq*(lmMst2 + lmMt) - 36*pow2(lmMsq))*pow6(Mst2)) + 240*pow2(Msq)*(6*
        pow2(Mst1)*pow2(Mst2)*((-(lmMsq*(7 + 42*lmMst1 - 36*lmMst2 + 6*lmMt)) +
        lmMst1*(103 + 42*lmMt) - 4*(32 + 9*lmMst2*lmMt + 12*(lmMst2 + lmMt)) +
        6*pow2(lmMsq))*pow2(Mst2) - 3*(41 - 10*lmMst1 + 6*lmMsq*(-2 + lmMst1 -
        lmMst2) + 22*lmMst2 - 3*pow2(lmMst1) + 3*pow2(lmMst2))*pow2(MuSUSY)) +
        18*((11 - 15*lmMst1 + 15*lmMst2)*(1 + 2*lmMst2 - 2*lmMt)*pow2(Mst2) - (
        29 + 48*lmMst2 - 24*pow2(lmMst1) - 6*(lmMst1*(6 - 5*lmMst2) + lmMsq*(2
        - 3*lmMst1 + 3*lmMst2) + pow2(lmMst2)))*pow2(MuSUSY))*pow4(Mst1) + 9*
        pow2(MuSUSY)*pow4(Mst2) + (67 + lmMst2*(24 - 36*lmMt) - 66*lmMt + 6*
        lmMsq*(7 + 6*lmMst2 + 6*lmMt) - 36*pow2(lmMsq))*pow6(Mst2))))))/(3888.*
        pow4(Msq)*pow4(Mst2)) + pow2(Mt)*(pow2(s2t)*((Dmglst1*Mst1*(156781 -
        864*B4 - 432*DN - 10260*lmMsq - 30*(95 + 324*lmMsq)*lmMst1 + 9720*pow2(
        lmMsq) + 7524*pow2(lmMst1) + 6*lmMst2*(21991 - 1620*lmMsq + 372*lmMst1
        + 432*pow2(lmMst1)) + 36*(1025 + 168*lmMst1)*pow2(lmMst2) + 5184*lmMt*(
        2 + 3*lmMst2 - lmMst1*(3 + 2*lmMst2) + pow2(lmMst1) + pow2(lmMst2)) -
        7200*pow3(lmMst1) - 1440*pow3(lmMst2)))/243. + pow2(Mst1)*(
        268.4533251028807 - (8*DN)/9. - (100*lmMsq)/3. + 20*pow2(lmMsq) - (2*
        lmMst1*(4703 + 500*lmMsq + 375*pow2(lmMsq)))/225. + (2*lmMst2*(33703 +
        1745*lmMst1 - 250*lmMsq*(16 + 3*lmMst1) + 375*pow2(lmMsq) - 1025*pow2(
        lmMst1)))/225. + ((-273 + 100*lmMsq)*pow2(lmMst1))/15. + ((4601 + 1770*
        lmMst1)*pow2(lmMst2))/45. + (32*lmMt*pow2(1 - lmMst1 + lmMst2))/3. - (
        64*pow3(lmMst1))/3. - (80*pow3(lmMst2))/9.) + ((-2*(750*Dmglst1*(131 -
        120*lmMsq + 78*lmMst1 + 42*lmMst2) + Mst1*(22642 + 15*lmMsq*(-1282 +
        195*lmMst1 - 75*lmMst2) + 7500*lmMst2 + 15*lmMst1*(782 + 75*lmMst2) -
        900*pow2(lmMsq) - 2025*pow2(lmMst1))))/(2025.*pow2(Msq)) - (Mst1*(
        238.25271580142933 - (50*lmMsq)/3. + (543.1291609977325 - (100*lmMsq)/
        3.)*lmMst1 + 10*pow2(lmMsq) + (266299*pow2(lmMst1))/945. + (lmMst2*(-
        14193447 + 441000*lmMsq - 9820930*lmMst1 + 455700*pow2(lmMst1)))/33075.
         + (48.90899470899471 - (700*lmMst1)/9.)*pow2(lmMst2) - (64*lmMt*(-
        lmMst1 + lmMst2 - 2*lmMst1*lmMst2 + pow2(lmMst1) + pow2(lmMst2)))/3. +
        (836*pow3(lmMst1))/27. + (892*pow3(lmMst2))/27.) + (Dmglst1*(63295 -
        222366*lmMst2 + 1620*lmMsq*(9 - 20*lmMst1 + 14*lmMst2) + 4860*pow2(
        lmMsq) - 144*(-1355 + 39*lmMst2 + 144*lmMt)*pow2(lmMst1) + 6*lmMst1*(
        42083 + 6048*lmMt + 48*lmMst2*(-953 + 144*lmMt) - 9432*pow2(lmMst2)) +
        97488*pow2(lmMst2) - 5184*lmMt*(2 + 7*lmMst2 + 4*pow2(lmMst2)) + 29520*
        pow3(lmMst1) + 32688*pow3(lmMst2)))/243.)/pow2(Mst2))*pow3(Mst1) + ((2*
        (5108 + 15*lmMsq*(-218 + 75*lmMst1 - 195*lmMst2) + 6270*lmMst2 - 375*
        lmMst1*(8 + 3*lmMst2) + 900*pow2(lmMsq) + 2025*pow2(lmMst2)))/(2025.*
        pow2(Msq)) + (2*(-1 + 2*lmMst2))/(3.*pow2(Mst1)) - (Mst1*((5*Dmglst1*(
        13 - 12*lmMsq + 12*lmMst2))/18. + Mst1*(-(lmMsq*(2.247165532879819 +
        lmMst1 - (11*lmMst2)/21.)) + lmMst1*(1.4 + lmMst2) + (5*pow2(lmMsq))/
        21. + (2*(54751 + 19614*lmMst2 - 17640*pow2(lmMst2)))/46305.)))/pow4(
        Msq))*pow4(Mst2) + (2*((8*OepS2 - 27*(37 + 6*lmMst1 - 6*lmMst2)*S2)*
        pow2(Mst2) + (-1480*OepS2 + 27*(2069 + 1110*lmMst1 - 1110*lmMst2)*S2)*
        pow2(MuSUSY))*pow4(Mst1) - 27*((8*OepS2 - 81*(1 + 2*lmMst1 - 2*lmMst2)*
        S2)*pow2(Mst2) + (8*OepS2 + 81*(11 - 2*lmMst1 + 2*lmMst2)*S2)*pow2(
        MuSUSY))*pow4(Mst2) - 3*pow2(Mst1)*(4*(136*OepS2 - 81*(49 + 34*lmMst1 -
        34*lmMst2)*S2)*pow2(Mst2)*pow2(MuSUSY) + (40*OepS2 - 81*(61 + 10*lmMst1
        - 10*lmMst2)*S2)*pow4(Mst2)))/(729.*pow4(Mst2)) - pow2(MuSUSY)*(
        83.61265432098766 + 4*B4 - (4*DN)/9. - (185*lmMsq)/9. + (25*pow2(lmMsq)
        )/3. - (lmMst1*(3781 - 300*lmMsq + 180*pow2(lmMsq)))/108. + (lmMst2*(
        14065 - 3498*lmMst1 + 60*lmMsq*(-35 + 12*lmMst1) + 180*pow2(lmMsq) -
        18*pow2(lmMst1)))/108. + (6.361111111111111 - (5*lmMsq)/3.)*pow2(
        lmMst1) - ((-2193 + 180*lmMsq + 442*lmMst1)*pow2(lmMst2))/36. +
        Dmglst1*((7 - 2*lmMst1 - 60*lmMsq*(-5 + 6*lmMst1) - 226*lmMst2 + 220*
        lmMst1*lmMst2 + 180*pow2(lmMsq) + 178*pow2(lmMst1) + 18*pow2(lmMst2))/(
        18.*Mst1) + (10*Mst1*(82 - 93*lmMst1 + 6*lmMsq*(4 + 3*lmMst1 - 3*
        lmMst2) + 69*lmMst2 - 9*pow2(lmMst1) + 9*pow2(lmMst2)))/(27.*pow2(Msq))
        ) + (11*pow3(lmMst1))/18. + (71*pow3(lmMst2))/6. + (Mst1*(-(Mst1*(
        26.575942386831276 - (76*B4)/9. + (2*DN)/9. + (35*lmMsq)/3. - 5*pow2(
        lmMsq) + lmMst1*(224.2274074074074 - (380*lmMsq)/9. + (40*pow2(lmMsq))/
        3.) - (149*pow2(lmMst1))/90. - (lmMst2*(372457 - 242670*lmMst1 + 1500*
        lmMsq*(-47 + 24*lmMst1) + 18000*pow2(lmMsq) + 17700*pow2(lmMst1)))/
        1350. + ((-17709 + 2400*lmMsq + 6940*lmMst1)*pow2(lmMst2))/90. + (8*
        pow3(lmMst1))/27. - (1736*pow3(lmMst2))/27.)) + (2*Dmglst1*(586 + 36*B4
        + 30*DN - 495*lmMsq + 135*pow2(lmMsq) - 6*lmMst1*(-73 - 45*lmMsq + 45*
        pow2(lmMsq)) + 3*lmMst2*(229 - 180*lmMsq + 86*lmMst1 + 90*pow2(lmMsq) -
        79*pow2(lmMst1)) + 18*(-43 + 15*lmMsq)*pow2(lmMst1) - 3*(-372 + 90*
        lmMsq + 49*lmMst1)*pow2(lmMst2) - 65*pow3(lmMst1) + 449*pow3(lmMst2)))/
        27.))/pow2(Mst2) + pow2(Mst2)*((1 - 2*lmMst2)/(3.*pow2(Mst1)) - ((5*
        Dmglst1*Mst1)/6. + (1.7499706655148832 + lmMsq*(-0.990854119425548 + (
        8*lmMst1)/9. - (52*lmMst2)/63.) + (10511*lmMst2)/4410. - (4*lmMst1*(47
        + 30*lmMst2))/135. - (2*pow2(lmMsq))/63. + (6*pow2(lmMst2))/7.)*pow2(
        Mst1) + (0.6658120973257028 + lmMsq*(-0.36171579743008314 + lmMst1/9. -
        lmMst2/7.) + (15647*lmMst2)/26460. - (lmMst1*(31 + 15*lmMst2))/135. +
        pow2(lmMsq)/63. + (8*pow2(lmMst2))/63.)*pow2(Mst2))/pow4(Msq)) + pow3(
        Mst1)*((514500*Dmglst1*(485 - 346*lmMst1 + 12*lmMsq*(1 + 10*lmMst1 -
        10*lmMst2) + 334*lmMst2 - 60*pow2(lmMst1) + 60*pow2(lmMst2)) - Mst1*(-
        41220947 + 420*lmMsq*(12479 + 4830*lmMst1 - 5670*lmMst2) + 1081710*
        lmMst2 - 210*lmMst1*(30109 + 123480*lmMst2) + 176400*pow2(lmMsq) +
        11951100*pow2(lmMst1) + 14156100*pow2(lmMst2)))/(1.11132e7*pow4(Msq)) -
        (Mst1*(265.6519022158834 - (76*B4)/9. + (2*DN)/9. - (25*lmMsq)/2. +
        lmMst1*(426.37458616780043 - (605*lmMsq)/9. + (40*pow2(lmMsq))/3.) - (
        66.24060846560846 - (40*lmMsq)/3.)*pow2(lmMst1) - lmMst2*(
        411.7079195011338 - (1386139*lmMst1)/3780. + (5*lmMsq*(-121 + 96*
        lmMst1))/9. + (40*pow2(lmMsq))/3. + (298*pow2(lmMst1))/3.) - (
        298.6850529100529 - 40*lmMsq - (2174*lmMst1)/9.)*pow2(lmMst2) + (80*
        pow3(lmMst1))/27. - (3920*pow3(lmMst2))/27.) + Dmglst1*(
        112.49099794238683 - (8*B4)/3. - (20*DN)/9. + lmMst1*(
        254.24382716049382 - 120*lmMsq + 20*pow2(lmMsq)) + (129.57407407407408
         - 20*lmMsq)*pow2(lmMst1) - lmMst2*(
        211.57716049382717 - 120*lmMsq - (1439*lmMst1)/27. + 20*pow2(lmMsq) + (
        338*pow2(lmMst1))/9.) + ((-9491 + 1080*lmMsq + 6636*lmMst1)*pow2(
        lmMst2))/54. + (38*pow3(lmMst1))/9. - (806*pow3(lmMst2))/9.))/pow4(
        Mst2)) + ((836 - 353*lmMst2 + 30*lmMsq*(7 - 4*lmMst1 + 4*lmMst2) +
        lmMst1*(143 + 750*lmMst2) - 315*pow2(lmMst1) - 435*pow2(lmMst2))*pow2(
        Mst1) - (100*Dmglst1*(55 + lmMst1 + 6*lmMsq*(-5 + 7*lmMst1 - 7*lmMst2)
        + 29*lmMst2 - 21*pow2(lmMst1) + 21*pow2(lmMst2))*pow3(Mst1) + (1525 +
        30*lmMsq*(-25 + 44*lmMst1 - 44*lmMst2) + 1928*lmMst2 - 2*lmMst1*(589 +
        900*lmMst2) + 240*pow2(lmMst1) + 1560*pow2(lmMst2))*pow4(Mst1) + (939 +
        30*lmMsq*(-12 + 5*lmMst1 - 5*lmMst2) + 760*lmMst2 - 50*lmMst1*(8 + 3*
        lmMst2) + 150*pow2(lmMst2))*pow4(Mst2))/pow2(Mst2))/(135.*pow2(Msq))) -
        (pow2(Mst2)*(2744*(3091 + 10020*lmMst2 + 30*lmMst1*(-491 + 150*lmMst2)
        + 30*lmMsq*(157 - 60*(lmMst1 + lmMst2)) + 1800*pow2(lmMsq) - 1350*(
        pow2(lmMst1) + pow2(lmMst2)))*pow2(Msq)*pow3(Mst1) + 8575*Mst1*(4501 -
        288*DN - 10896*lmMst1 + 1080*(3 - lmMst1 + lmMst2)*pow2(lmMsq) + 2952*
        pow2(lmMst1) - 24*lmMst2*(-733 + 336*lmMst1 + 36*pow2(lmMst1)) - 72*(-
        191 + 44*lmMst1)*pow2(lmMst2) - 360*lmMsq*(3 + lmMst1*(4 - 6*lmMst2) +
        14*lmMst2 + 6*pow2(lmMst2)) - 696*pow3(lmMst1) + 4728*pow3(lmMst2))*
        pow4(Msq) - 51450*Dmglst1*(120*(7 - 20*lmMsq + 26*lmMst1 - 6*lmMst2)*
        pow2(Msq)*pow2(Mst1) + 4*(158 + 60*lmMsq*(-4 + 9*lmMst1) + lmMst1*(70 -
        138*lmMst2) + 96*lmMst2 - 270*pow2(lmMsq) - 459*pow2(lmMst1) - 27*pow2(
        lmMst2))*pow4(Msq) + 5*(17 - 324*lmMsq + 408*lmMst1 - 84*lmMst2)*pow4(
        Mst1)) - 60*(130598 + lmMst1*(423822 - 46305*lmMst2) - 105*lmMsq*(2929
        + 231*lmMst1 - 441*lmMst2) - 116277*lmMst2 - 11025*pow2(lmMsq) + 35280*
        pow2(lmMst1))*pow5(Mst1)))/(2.7783e6*Mst1*pow4(Msq))) + (pow2(MuSUSY)*(
        pow2(s2t)*(83.61265432098766 + 4*B4 - (4*DN)/9. - (185*lmMsq)/9. + (25*
        pow2(lmMsq))/3. - (lmMst1*(3781 - 300*lmMsq + 180*pow2(lmMsq)))/108. +
        (lmMst2*(14065 - 3498*lmMst1 + 60*lmMsq*(-35 + 12*lmMst1) + 180*pow2(
        lmMsq) - 18*pow2(lmMst1)))/108. + ((229 - 60*lmMsq)*pow2(lmMst1))/36. -
        ((-2193 + 180*lmMsq + 442*lmMst1)*pow2(lmMst2))/36. + Dmglst1*((7 - 2*
        lmMst1 - 60*lmMsq*(-5 + 6*lmMst1) - 226*lmMst2 + 220*lmMst1*lmMst2 +
        180*pow2(lmMsq) + 178*pow2(lmMst1) + 18*pow2(lmMst2))/(18.*Mst1) + (10*
        Mst1*(82 - 93*lmMst1 + 6*lmMsq*(4 + 3*lmMst1 - 3*lmMst2) + 69*lmMst2 -
        9*pow2(lmMst1) + 9*pow2(lmMst2)))/(27.*pow2(Msq))) + (11*pow3(lmMst1))/
        18. + (71*pow3(lmMst2))/6. + (Mst1*(-(Mst1*(26.575942386831276 - (76*
        B4)/9. + (2*DN)/9. + (35*lmMsq)/3. - 5*pow2(lmMsq) + lmMst1*(
        224.2274074074074 - (380*lmMsq)/9. + (40*pow2(lmMsq))/3.) - (149*pow2(
        lmMst1))/90. - (lmMst2*(372457 - 242670*lmMst1 + 1500*lmMsq*(-47 + 24*
        lmMst1) + 18000*pow2(lmMsq) + 17700*pow2(lmMst1)))/1350. + ((-17709 +
        2400*lmMsq + 6940*lmMst1)*pow2(lmMst2))/90. + (8*pow3(lmMst1))/27. - (
        1736*pow3(lmMst2))/27.)) + (2*Dmglst1*(586 + 36*B4 + 30*DN - 495*lmMsq
        + 135*pow2(lmMsq) - 6*lmMst1*(-73 - 45*lmMsq + 45*pow2(lmMsq)) + 3*
        lmMst2*(229 - 180*lmMsq + 86*lmMst1 + 90*pow2(lmMsq) - 79*pow2(lmMst1))
        + 18*(-43 + 15*lmMsq)*pow2(lmMst1) - 3*(-372 + 90*lmMsq + 49*lmMst1)*
        pow2(lmMst2) - 65*pow3(lmMst1) + 449*pow3(lmMst2)))/27.))/pow2(Mst2) +
        pow2(Mst2)*((1 - 2*lmMst2)/(3.*pow2(Mst1)) - ((5*Dmglst1*Mst1)/6. + (
        1.7499706655148832 + lmMsq*(-0.990854119425548 + (8*lmMst1)/9. - (52*
        lmMst2)/63.) + (10511*lmMst2)/4410. - (4*lmMst1*(47 + 30*lmMst2))/135.
         - (2*pow2(lmMsq))/
        63. + (6*pow2(lmMst2))/7.)*pow2(Mst1) + (0.6658120973257028 + lmMsq*(-
        0.36171579743008314 + lmMst1/9. - lmMst2/7.) + (15647*lmMst2)/26460. -
        (lmMst1*(31 + 15*lmMst2))/135. + pow2(lmMsq)/63. + (8*pow2(lmMst2))/63.
        )*pow2(Mst2))/pow4(Msq)) + pow3(Mst1)*((514500*Dmglst1*(485 - 346*
        lmMst1 + 12*lmMsq*(1 + 10*lmMst1 - 10*lmMst2) + 334*lmMst2 - 60*pow2(
        lmMst1) + 60*pow2(lmMst2)) - Mst1*(-41220947 + 420*lmMsq*(12479 + 4830*
        lmMst1 - 5670*lmMst2) + 1081710*lmMst2 - 210*lmMst1*(30109 + 123480*
        lmMst2) + 176400*pow2(lmMsq) + 11951100*pow2(lmMst1) + 14156100*pow2(
        lmMst2)))/(1.11132e7*pow4(Msq)) - (Mst1*(265.6519022158834 - (76*B4)/9.
         + (2*DN)/9. - (25*lmMsq)/2. + lmMst1*(426.37458616780043 - (605*lmMsq)
        /9. + (40*pow2(lmMsq))/3.) - (66.24060846560846 - (40*lmMsq)/3.)*pow2(
        lmMst1) - lmMst2*(411.7079195011338 - (1386139*lmMst1)/3780. + (5*
        lmMsq*(-121 + 96*lmMst1))/9. + (40*pow2(lmMsq))/3. + (298*pow2(lmMst1))
        /3.) - (298.6850529100529 - 40*lmMsq - (2174*lmMst1)/9.)*pow2(lmMst2) +
        (80*pow3(lmMst1))/27. - (3920*pow3(lmMst2))/27.) + Dmglst1*(
        112.49099794238683 - (8*B4)/3. - (20*DN)/9. + lmMst1*(
        254.24382716049382 - 120*lmMsq + 20*pow2(lmMsq)) + (129.57407407407408
         - 20*lmMsq)*pow2(lmMst1) - lmMst2*(
        211.57716049382717 - 120*lmMsq - (1439*lmMst1)/27. + 20*pow2(lmMsq) + (
        338*pow2(lmMst1))/9.) + ((-9491 + 1080*lmMsq + 6636*lmMst1)*pow2(
        lmMst2))/54. + (38*pow3(lmMst1))/9. - (806*pow3(lmMst2))/9.))/pow4(
        Mst2)) + (S2*(81*(11 - 2*lmMst1 + 2*lmMst2) - (36*(49 + 34*lmMst1 - 34*
        lmMst2)*pow2(Mst1))/pow2(Mst2) - (2*(2069 + 1110*lmMst1 - 1110*lmMst2)*
        pow4(Mst1))/pow4(Mst2)))/27. + (8*OepS2*(204*pow2(Mst1)*pow2(Mst2) +
        370*pow4(Mst1) + 27*pow4(Mst2)))/(729.*pow4(Mst2)) + ((836 - 353*lmMst2
        + 30*lmMsq*(7 - 4*lmMst1 + 4*lmMst2) + lmMst1*(143 + 750*lmMst2) - 315*
        pow2(lmMst1) - 435*pow2(lmMst2))*pow2(Mst1) - (100*Dmglst1*(55 + lmMst1
        + 6*lmMsq*(-5 + 7*lmMst1 - 7*lmMst2) + 29*lmMst2 - 21*pow2(lmMst1) +
        21*pow2(lmMst2))*pow3(Mst1) + (1525 + 30*lmMsq*(-25 + 44*lmMst1 - 44*
        lmMst2) + 1928*lmMst2 - 2*lmMst1*(589 + 900*lmMst2) + 240*pow2(lmMst1)
        + 1560*pow2(lmMst2))*pow4(Mst1) + (939 + 30*lmMsq*(-12 + 5*lmMst1 - 5*
        lmMst2) + 760*lmMst2 - 50*lmMst1*(8 + 3*lmMst2) + 150*pow2(lmMst2))*
        pow4(Mst2))/pow2(Mst2))/(135.*pow2(Msq))) + (Mt*(144*Mst1*Mt*(2*
        Dmglst1*(735 - 6*B4 - 3*DN - 78*lmMst1 + 6*pow2(lmMst1) + 6*lmMst2*(85
        - 6*lmMst1 + pow2(lmMst1)) + 18*(7 + lmMst1)*pow2(lmMst2) + 24*lmMt*(2
        + 3*lmMst2 - lmMst1*(3 + 2*lmMst2) + pow2(lmMst1) + pow2(lmMst2)) - 26*
        pow3(lmMst1) + 2*pow3(lmMst2)) + Mst1*(417 - 6*B4 - 3*DN - 162*lmMst1 +
        6*lmMst2*(87 - 8*lmMst1 + pow2(lmMst1)) + 18*(8 + lmMst1)*pow2(lmMst2)
        + 24*lmMt*pow2(1 - lmMst1 + lmMst2) - 26*pow3(lmMst1) + 2*pow3(lmMst2))
        ) + (s2t*(3*Dmglst1*(2*(pow2(Mst1)*(13969 - 11880*B4 + 108*DN - 4320*
        lmMsq + (-16728 - 3240*lmMsq*(1 + lmMsq))*lmMst1 + 2160*pow2(lmMsq) +
        54*(-407 + 300*lmMsq)*pow2(lmMst1) + 12*lmMst2*(4094 + 2205*lmMst1 -
        90*lmMsq*(1 + 24*lmMst1) + 270*pow2(lmMsq) + 804*pow2(lmMst1)) + 18*(
        583 + 540*lmMsq + 680*lmMst1)*pow2(lmMst2) - 16848*pow3(lmMst1) - 5040*
        pow3(lmMst2)) + 3*pow2(Mst2)*(11337 - 408*B4 - 12*DN - 6480*lmMsq +
        720*pow2(lmMsq) - 30*lmMst1*(-73 - 36*lmMsq + 12*pow2(lmMsq)) + 6*
        lmMst2*(2331 - 420*lmMsq - 118*lmMst1 + 60*pow2(lmMsq) - 60*pow2(
        lmMst1)) + 18*(-13 + 20*lmMsq)*pow2(lmMst1) - 6*(-797 + 60*lmMsq + 4*
        lmMst1)*pow2(lmMst2) - 200*pow3(lmMst1) + 584*pow3(lmMst2)))*pow4(Msq)
        - 720*pow2(Msq)*(2*(41 - 10*lmMst1 + 6*lmMsq*(-2 + lmMst1 - lmMst2) +
        22*lmMst2 - 3*pow2(lmMst1) + 3*pow2(lmMst2))*pow2(Mst1)*pow2(Mst2) + 2*
        (29 + 48*lmMst2 - 24*pow2(lmMst1) - 6*(lmMst1*(6 - 5*lmMst2) + lmMsq*(2
        - 3*lmMst1 + 3*lmMst2) + pow2(lmMst2)))*pow4(Mst1) - pow4(Mst2)) - 5*
        pow2(Mst2)*((-251 + 12*lmMsq - 12*lmMst2)*pow2(Mst1)*pow2(Mst2) + 6*(
        818 - 117*lmMst1 + 180*lmMsq*(-2 + lmMst1 - lmMst2) + 477*lmMst2 - 90*
        pow2(lmMst1) + 90*pow2(lmMst2))*pow4(Mst1) + (11 - 12*lmMsq + 12*
        lmMst2)*pow4(Mst2))) + Mst1*(2*(9*pow2(Mst2)*(13317 - 408*B4 - 12*DN -
        4320*lmMsq + 720*pow2(lmMsq) - 6*lmMst1*(79 - 180*lmMsq + 60*pow2(
        lmMsq)) + 6*lmMst2*(2135 - 420*lmMsq - 118*lmMst1 + 60*pow2(lmMsq) -
        60*pow2(lmMst1)) + 6*(-167 + 60*lmMsq)*pow2(lmMst1) - 6*(-797 + 60*
        lmMsq + 4*lmMst1)*pow2(lmMst2) - 200*pow3(lmMst1) + 584*pow3(lmMst2)) +
        pow2(Mst1)*(107299 - 19224*B4 - 108*DN - 32400*lmMsq + 6480*pow2(lmMsq)
        - 24*lmMst1*(6103 - 1755*lmMsq + 405*pow2(lmMsq)) + 18*(-4379 + 1260*
        lmMsq)*pow2(lmMst1) + 12*lmMst2*(20882 + 699*lmMst1 - 270*lmMsq*(17 +
        8*lmMst1) + 810*pow2(lmMsq) + 54*pow2(lmMst1)) + 18*(5641 + 180*lmMsq +
        540*lmMst1)*pow2(lmMst2) - 18072*pow3(lmMst1) + 7704*pow3(lmMst2)))*
        pow4(Msq) - 720*pow2(Msq)*(2*(34 - 8*lmMst1 + 6*lmMsq*(-2 + lmMst1 -
        lmMst2) + 20*lmMst2 - 3*pow2(lmMst1) + 3*pow2(lmMst2))*pow2(Mst1)*pow2(
        Mst2) + 4*(14 - 21*lmMst1 + lmMsq*(-6 + 9*lmMst1 - 9*lmMst2) + 27*
        lmMst2 + 9*lmMst1*lmMst2 - 9*pow2(lmMst1))*pow4(Mst1) - 3*pow4(Mst2)) -
        15*(6*(134 - 9*lmMst1 + 36*lmMsq*(-2 + lmMst1 - lmMst2) + 81*lmMst2 -
        18*pow2(lmMst1) + 18*pow2(lmMst2))*pow2(Mst2)*pow4(Mst1) + (-107 + 12*
        lmMsq - 12*lmMst2)*pow2(Mst1)*pow4(Mst2) + (11 - 12*lmMsq + 12*lmMst2)*
        pow6(Mst2)))))/pow4(Msq)))/(486.*pow4(Mst2))))/pow2(Sbeta)) + (z2*(-36*
        Mt*pow3(s2t)*(6*Dmglst1*(1941 - 660*lmMsq - 220*lmMst1 + 1244*lmMst2)*
        pow2(Mst1) + 2*(1355 - 540*lmMsq + 24*lmMst1 + 744*lmMst2)*pow3(Mst1) +
        ((-480*(6*Dmglst1 + Mst1))/pow2(Msq) - (Dmglst1*(1397 + 5124*lmMst1 -
        5124*lmMst2) + (737 + 852*lmMst1 - 852*lmMst2)*Mst1)/pow2(Mst2))*pow4(
        Mst1) - 6*pow2(Mst2)*((627 - 60*lmMsq - 64*lmMst1 + 192*lmMst2)*(
        Dmglst1 + Mst1) - (40*(3*Dmglst1 + Mst1)*pow2(Mst1))/pow2(Msq) - (15*(
        5*Dmglst1 + Mst1)*pow4(Mst1))/pow4(Msq))) - (pow2(Mt)*pow2(MuSUSY)*(
        270*s2t*pow2(Mst2)*(pow2(Mst1)*(24*Mst1*Mt - 5*s2t*pow2(Mst2)) + 4*
        pow2(Msq)*(16*Mst1*Mt + 4*s2t*pow2(Mst1) - 3*s2t*pow2(Mst2)))*pow4(
        Mst1) + 36*Dmglst1*(Mst1*pow4(Msq)*(12*Mst1*Mt*s2t*((687 - 540*lmMsq -
        92*lmMst1 + 860*lmMst2)*pow2(Mst1) + (-627 + 60*lmMsq + 64*lmMst1 -
        192*lmMst2)*pow2(Mst2)) - 12*pow2(Mst1)*(8*(55 + 16*lmMst2)*pow2(Mt) +
        (390 - 90*lmMsq - 139*lmMst1 + 315*lmMst2)*pow2(Mst2)*pow2(s2t)) +
        pow2(s2t)*((-1613 + 1080*lmMsq + 3684*lmMst1 - 5796*lmMst2)*pow4(Mst1)
        - 594*pow4(Mst2))) + s2t*(-60*pow2(Msq)*pow3(Mst1)*(-24*Mst1*Mt*pow2(
        Mst2) - 14*s2t*pow2(Mst1)*pow2(Mst2) + 48*Mt*pow3(Mst1) + 3*s2t*pow4(
        Mst2)) + 150*pow2(Mst2)*(6*Mst1*Mt - s2t*pow2(Mst2))*pow5(Mst1))) -
        pow4(Msq)*(144*Mt*s2t*((2407 + 180*lmMsq - 408*lmMst1 + 408*lmMst2)*
        pow2(Mst1) + 3*(627 - 60*lmMsq - 64*lmMst1 + 192*lmMst2)*pow2(Mst2))*
        pow3(Mst1) + 12*(144*(55 + 16*lmMst2)*pow2(Mt) - (311 + 1080*lmMsq +
        2196*lmMst1 - 3348*lmMst2)*pow2(Mst2)*pow2(s2t))*pow4(Mst1) + pow2(s2t)
        *(108*(-339 + 30*lmMsq + 37*lmMst1 - 149*lmMst2)*pow2(Mst1)*pow4(Mst2)
        - (69349 + 12960*lmMsq + 93744*lmMst1 - 107568*lmMst2)*pow6(Mst1) -
        648*pow6(Mst2)))) - 162*shiftst3*pow4(Msq)*(-4*pow2(Mst1)*pow2(s2t)*
        pow4(Mst2)*(2*(-1 + lmMst1 - lmMst2)*pow2(Mt)*(-(pow2(MuSUSY)*(-1 +
        pow2(Sbeta))) + pow2(Mst2)*pow2(Sbeta)) + pow2(s2t)*pow2(Sbeta)*pow4(
        Mst2)) + pow2(Mt)*(8*(-1 + 3*lmMst1 - 3*lmMst2)*pow2(MuSUSY)*pow2(s2t)*
        (-1 + pow2(Sbeta))*pow6(Mst1) + 4*(-(pow2(MuSUSY)*pow2(s2t)*(-1 + pow2(
        Sbeta))) + 2*(2*pow2(Mt) - pow2(Mst2)*pow2(s2t))*pow2(Sbeta))*pow6(
        Mst2)) + pow4(Mst1)*(8*(-1 + 2*lmMst1 - 2*lmMst2)*pow2(Mst2)*pow2(Mt)*
        pow2(MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta)) + (3 + 2*lmMst1 - 2*lmMst2)*
        pow2(Sbeta)*pow4(s2t)*pow6(Mst2))) + pow2(Mt)*pow2(Sbeta)*(pow2(s2t)*(
        36*Dmglst1*(Mst1*pow4(Msq)*(-4*pow2(Mst1)*pow2(Mst2)*(2*(1537 + 60*
        lmMst1 + 228*lmMst2)*pow2(Mst2) + 3*(-390 + 90*lmMsq + 139*lmMst1 -
        315*lmMst2)*pow2(MuSUSY)) + (8*(923 - 894*lmMst1 + 1470*lmMst2)*pow2(
        Mst2) + (1613 - 1080*lmMsq - 3684*lmMst1 + 5796*lmMst2)*pow2(MuSUSY))*
        pow4(Mst1) + 594*(2*pow2(Mst2) + pow2(MuSUSY))*pow4(Mst2)) + 30*pow2(
        MuSUSY)*(5*pow4(Mst2)*pow5(Mst1) + pow2(Msq)*(6*pow3(Mst1)*pow4(Mst2) -
        28*pow2(Mst2)*pow5(Mst1)))) + 270*pow2(MuSUSY)*(5*pow4(Mst2)*pow6(Mst1)
        - 4*pow2(Msq)*(-3*pow4(Mst1)*pow4(Mst2) + 4*pow2(Mst2)*pow6(Mst1))) -
        pow4(Msq)*(6*pow4(Mst1)*(2*(311 + 1080*lmMsq + 2196*lmMst1 - 3348*
        lmMst2)*pow2(Mst2)*pow2(MuSUSY) + (40285 + 4032*lmMst1 + 2880*lmMst2)*
        pow4(Mst2)) + (20*(59 + 3564*lmMst1 - 3564*lmMst2)*pow2(Mst2) + (69349
        + 12960*lmMsq + 93744*lmMst1 - 107568*lmMst2)*pow2(MuSUSY))*pow6(Mst1)
        + 648*(2*pow2(Mst2) + pow2(MuSUSY))*pow6(Mst2) + 54*pow2(Mst1)*(2*(339
        - 30*lmMsq - 37*lmMst1 + 149*lmMst2)*pow2(MuSUSY)*pow4(Mst2) + (-611 +
        24*lmMst1 - 24*lmMst2)*pow6(Mst2)))) - 144*Mt*(s2t*pow2(Mst1)*(Mst1*
        pow4(Msq)*(3*(-627 + 60*lmMsq + 64*lmMst1 - 192*lmMst2)*pow2(Mst2)*
        pow2(MuSUSY) - pow2(Mst1)*(12*(-463 + 120*lmMsq + 27*lmMst1 - 283*
        lmMst2)*pow2(Mst2) + (2407 + 180*lmMsq - 408*lmMst1 + 408*lmMst2)*pow2(
        MuSUSY)) + (9178 - 2160*lmMsq - 552*lmMst1 + 4392*lmMst2)*pow4(Mst1) +
        4*(13 + 96*lmMst1 - 54*lmMst2 - 42*lmMt)*pow4(Mst2)) + Dmglst1*pow4(
        Msq)*(3*(-627 + 60*lmMsq + 64*lmMst1 - 192*lmMst2)*pow2(Mst2)*pow2(
        MuSUSY) + pow2(Mst1)*((16600 - 4320*lmMsq - 852*lmMst1 + 9900*lmMst2 +
        168*lmMt)*pow2(Mst2) + 3*(687 - 540*lmMsq - 92*lmMst1 + 860*lmMst2)*
        pow2(MuSUSY)) + 6*(7373 - 1800*lmMsq - 372*lmMst1 + 3544*lmMst2 + 28*
        lmMt)*pow4(Mst1) + 4*(13 + 96*lmMst1 - 54*lmMst2 - 42*lmMt)*pow4(Mst2))
        + 120*Dmglst1*pow2(Msq)*(-6*(5*pow2(Mst2) + pow2(MuSUSY))*pow4(Mst1) +
        pow2(Mst1)*(3*pow2(Mst2)*pow2(MuSUSY) - 2*pow4(Mst2)) + 2*pow6(Mst2)) +
        15*Mst1*pow2(Mst2)*(3*(-4*pow2(Mst2) + pow2(MuSUSY))*pow4(Mst1) - 8*
        pow2(Msq)*(pow2(Mst1)*(2*pow2(Mst2) - pow2(MuSUSY)) + 6*pow4(Mst1) - 2*
        pow4(Mst2)) + 8*pow2(Mst1)*pow4(Mst2) + 4*pow6(Mst2)) + 15*Dmglst1*
        pow2(Mst2)*((-28*pow2(Mst2) + 15*pow2(MuSUSY))*pow4(Mst1) + 24*pow2(
        Mst1)*pow4(Mst2) + 4*pow6(Mst2))) + Mt*(15*pow2(Mst1)*pow4(Mst2)*(6*
        pow2(Mst1)*pow2(Mst2) + 8*pow2(Msq)*(4*pow2(Mst1) + pow2(Mst2)) + 19*
        pow4(Mst1) + 3*pow4(Mst2)) - pow4(Msq)*(2*((-289 + 6*lmMst1 - 228*
        lmMst2 + 222*lmMt)*pow2(Mst2) + 6*(55 + 16*lmMst2)*pow2(MuSUSY))*pow4(
        Mst1) + (-392 - 360*lmMsq + 111*lmMst1 + 249*lmMt)*pow2(Mst1)*pow4(
        Mst2) + (2383 - 900*lmMsq - 474*lmMst1 + 1230*lmMst2 + 300*lmMt)*pow6(
        Mst1) + 18*pow6(Mst2)) + 2*Dmglst1*(360*pow2(Msq)*pow3(Mst1)*pow4(Mst2)
        + 390*pow4(Mst2)*pow5(Mst1) + pow4(Msq)*(4*((151 + 108*lmMst2 - 108*
        lmMt)*pow2(Mst2) - 3*(55 + 16*lmMst2)*pow2(MuSUSY))*pow3(Mst1) + 297*
        Mst1*pow4(Mst2) + 6*(-933 + 300*lmMsq + 150*lmMst1 - 462*lmMst2 - 40*
        lmMt)*pow5(Mst1)) + 90*pow3(Mst1)*pow6(Mst2))))))/(pow2(Mst1)*pow2(
        Sbeta)*pow4(Msq)*pow4(Mst2)) - 972*(pow4(s2t)*(Dmglst1*(
        66.23148148148148 - 10*lmMsq + (29*lmMst1)/9. + (49*lmMst2)/3.)*pow3(
        Mst1) - pow2(Mst2)*(Dmglst1*(32.333333333333336 - 10*lmMsq - (139*
        lmMst1)/9. + 35*lmMst2)*Mst1 + (17.70679012345679 - 5*lmMsq - (53*
        lmMst1)/6. + (335*lmMst2)/18.)*pow2(Mst1) - (25*(4*Dmglst1 + Mst1)*
        pow3(Mst1))/(9.*pow2(Msq))) + (25.333590534979425 - (25*lmMsq)/6. + (
        343*lmMst1)/36. - (103*lmMst2)/36.)*pow4(Mst1) - (pow4(Mst2)*(60*pow2(
        Msq)*pow3(Mst1) + (-654 + 60*lmMsq + 74*lmMst1 - 298*lmMst2)*Mst1*pow4(
        Msq) + 4*Dmglst1*(30*pow2(Msq)*pow2(Mst1) + 99*pow4(Msq) + 25*pow4(
        Mst1)) + 25*pow5(Mst1)))/(72.*Mst1*pow4(Msq))) + (MuSUSY*((pow2(Mt)*
        pow2(s2t)*(6*(627 - 60*lmMsq - 64*lmMst1 + 192*lmMst2)*(Dmglst1 + Mst1)
        + (pow2(Mst1)*(-240*(3*Dmglst1 + Mst1) + (4*(3*Dmglst1*(-657 + 300*
        lmMsq + 78*lmMst1 - 526*lmMst2)*pow2(Msq) + (263 + 180*lmMsq - 108*
        lmMst1 - 84*lmMst2)*Mst1*pow2(Msq) + 540*Dmglst1*pow2(Mst1) + 60*pow3(
        Mst1)))/pow2(Mst2)))/pow2(Msq) + pow4(Mst1)*((-90*(5*Dmglst1 + Mst1))/
        pow4(Msq) + (Dmglst1*(-6487 + 3600*lmMsq + 6060*lmMst1 - 11436*lmMst2)
        + (1789 + 720*lmMsq + 420*lmMst1 - 1188*lmMst2)*Mst1)/pow4(Mst2))))/9.
         - (2*(4*s2t*shiftst3*pow2(Mst2)*pow3(Mt) + Mt*pow3(s2t)*pow4(Mst2)))/(
        3.*pow2(Mst1)) + Mt*pow3(s2t)*(Dmglst1*(151.33333333333334 - 40*lmMsq -
        (556*lmMst1)/9. + 140*lmMst2)*Mst1 + ((2740 - 1350*lmMsq - 2529*lmMst1
        + 4689*lmMst2)*pow2(Mst1))/81. + (2*shiftst3*pow2(Mst2)*((-1 - lmMst1 +
        lmMst2)*pow2(Mst1) + pow2(Mst2)))/(3.*pow2(Mst1)) + ((-10*(34*Dmglst1 +
        7*Mst1))/(9.*pow2(Msq)) - (36*Dmglst1*(3067 + 2016*lmMst1 - 2016*
        lmMst2) + (65617 + 67392*lmMst1 - 67392*lmMst2)*Mst1)/(972.*pow2(Mst2))
        )*pow3(Mst1) + (pow2(Mst2)*(60*pow2(Msq)*pow3(Mst1) + (-666 + 60*lmMsq
        + 74*lmMst1 - 298*lmMst2)*Mst1*pow4(Msq) + 4*Dmglst1*(30*pow2(Msq)*
        pow2(Mst1) + 99*pow4(Msq) + 25*pow4(Mst1)) + 25*pow5(Mst1)))/(18.*Mst1*
        pow4(Msq))) + (pow3(Mt)*((s2t*(6*(17501 + 2124*lmMst1 + 1332*lmMst2)*
        pow2(Mst2)*pow4(Mst1) + 27*(-587 + 24*lmMst1 - 24*lmMst2)*pow2(Mst1)*
        pow4(Mst2) + 72*Dmglst1*((2777 + 120*lmMst1 + 456*lmMst2)*pow2(Mst2)*
        pow3(Mst1) - 297*Mst1*pow4(Mst2) + (931 + 1908*lmMst1 - 2484*lmMst2)*
        pow5(Mst1)) + 4*(26399 + 12096*lmMst1 - 6912*lmMst2)*pow6(Mst1) + 648*
        pow6(Mst2)))/pow2(Mst1) + (144*Mt*(Mst1*(((1381 - 360*lmMsq + 15*lmMst1
        + 795*lmMst2 - 42*lmMt)*pow2(Mst1) - 2*(4 - 48*lmMst1 + 27*lmMst2 + 21*
        lmMt)*pow2(Mst2))*pow4(Msq) - 60*pow2(Msq)*(3*pow4(Mst1) - pow4(Mst2))
        + 15*(3*pow2(Mst1)*pow4(Mst2) + pow6(Mst2))) + Dmglst1*(-(((-4163 +
        1080*lmMsq + 117*lmMst1 - 2421*lmMst2)*pow2(Mst1) + 2*(4 - 48*lmMst1 +
        27*lmMst2 + 21*lmMt)*pow2(Mst2))*pow4(Msq)) + 60*pow2(Msq)*(-15*pow4(
        Mst1) + pow4(Mst2)) + 15*(7*pow2(Mst1)*pow4(Mst2) + pow6(Mst2)))))/
        pow4(Msq)))/(243.*pow4(Mst2))))/Tbeta)))/972. + (MuSUSY*(((3*Mst1*(4*(
        pow2(Mst1)*(1747 - 1929*lmMst2 - 540*(lmMst1 - lmMst2)*pow2(lmMsq) +
        18*(-29 + 39*lmMst2 - 144*lmMt)*pow2(lmMst1) - 540*lmMsq*(7 + lmMst2 +
        4*lmMst1*lmMst2 + 5*lmMt + 6*lmMst2*lmMt - 6*lmMst1*(1 + lmMt) - 4*
        pow2(lmMst2)) - 5274*pow2(lmMst2) + 18*lmMt*(327 + 570*lmMst2 + 80*
        pow2(lmMst2)) + 3*lmMst1*(985 - 2268*lmMt + 12*lmMst2*(41 + 32*lmMt) +
        54*pow2(lmMst2) - 576*pow2(lmMt)) + 864*(1 + 2*lmMst2)*pow2(lmMt) +
        726*pow3(lmMst1) - 1590*pow3(lmMst2)) + 6*pow2(Mst2)*(212 + 225*lmMst2
        - 90*lmMsq*(7 + 2*lmMst2*lmMt - 2*lmMst1*(3 + lmMt) + 3*(lmMst2 + lmMt)
        ) - 90*(lmMst1 - lmMst2)*pow2(lmMsq) - 6*(67 + 5*lmMst2 + 14*lmMt)*
        pow2(lmMst1) + 102*pow2(lmMst2) + 3*lmMt*(265 + 202*lmMst2 + 36*pow2(
        lmMst2)) + 144*(1 + lmMst2)*pow2(lmMt) - 6*lmMst1*(-54 + 37*lmMt +
        lmMst2*(22 + 4*lmMt) + pow2(lmMst2) + 24*pow2(lmMt)) + 22*pow3(lmMst1)
        + 14*pow3(lmMst2)))*pow4(Msq) + 5*pow2(Mst2)*(3*(223 - 12*lmMsq + 18*
        lmMst2*(5 - 4*lmMt) - 78*lmMt + 72*lmMsq*(lmMst2 + lmMt) - 72*pow2(
        lmMsq))*pow2(Mst1)*pow2(Mst2) + 36*(-15 + lmMst2 + 12*lmMsq*lmMst2 - 7*
        lmMt - 12*lmMst2*lmMt + lmMst1*(6 - 12*lmMsq + 12*lmMt))*pow4(Mst1) + (
        299 - 60*lmMsq + 6*lmMst2*(23 - 12*lmMt) - 78*lmMt + 72*lmMsq*(lmMst2 +
        lmMt) - 72*pow2(lmMsq))*pow4(Mst2)) - 80*pow2(Msq)*(9*(11 + 6*lmMst1*(-
        2 + lmMsq - lmMt) + 5*lmMt + lmMst2*(7 - 6*lmMsq + 6*lmMt))*pow2(Mst1)*
        pow2(Mst2) + 9*(11 + 6*lmMst1*(-2 + lmMsq + lmMst2 - 2*lmMt) + 7*lmMt +
        lmMst2*(5 - 6*lmMsq + 12*lmMt) - 6*pow2(lmMst2))*pow4(Mst1) + (-17 +
        33*lmMt - 6*lmMsq*(5 + 3*lmMst2 + 3*lmMt) + 3*lmMst2*(-1 + 6*lmMt) +
        18*pow2(lmMsq))*pow4(Mst2))) - Dmglst1*(4*(6*pow2(Mst2)*(1846 - 1057*
        lmMst2 + 30*lmMsq*(151 + 51*lmMt + 9*lmMst2*(3 + 2*lmMt) - 6*lmMst1*(13
        + 3*lmMt)) + 270*(lmMst1 - lmMst2)*pow2(lmMsq) + 6*(352 + 15*lmMst2 +
        42*lmMt)*pow2(lmMst1) - 306*pow2(lmMst2) - 3*lmMt*(1105 + 608*lmMst2 +
        108*pow2(lmMst2)) - 432*(2 + lmMst2)*pow2(lmMt) + 2*lmMst1*(-2534 -
        192*lmMt + 3*lmMst2*(67 + 12*lmMt) + 9*pow2(lmMst2) + 216*pow2(lmMt)) -
        66*pow3(lmMst1) - 42*pow3(lmMst2)) + pow2(Mst1)*(16303 + 47667*lmMst2 +
        1620*(lmMst1 - lmMst2)*pow2(lmMsq) + (8226 - 8910*lmMst2 + 22788*lmMt)*
        pow2(lmMst1) + 180*lmMsq*(151 - 99*lmMst2 + 6*lmMst1*(-13 + 18*lmMst2 -
        21*lmMt) + 3*(59 + 42*lmMst2)*lmMt - 108*pow2(lmMst2)) + 54666*pow2(
        lmMst2) - 18*lmMt*(3271 + 3130*lmMst2 + 366*pow2(lmMst2)) - 3*lmMst1*(
        10535 - 5148*lmMt + 12*lmMst2*(323 + 450*lmMt) + 378*pow2(lmMst2) -
        3456*pow2(lmMt)) - 10368*(1 + lmMst2)*pow2(lmMt) - 4914*pow3(lmMst1) +
        14958*pow3(lmMst2)))*pow4(Msq) + 240*pow2(Msq)*((439 + 6*lmMst1*(-62 +
        21*lmMsq - 21*lmMt) + 177*lmMt + 3*lmMst2*(65 - 42*lmMsq + 42*lmMt))*
        pow2(Mst1)*pow2(Mst2) + (439 - 3*lmMst2*(1 + 42*lmMsq - 132*lmMt) + 6*
        lmMst1*(-62 + 21*lmMsq + 45*lmMst2 - 66*lmMt) + 375*lmMt - 270*pow2(
        lmMst2))*pow4(Mst1) - (17 + 3*lmMst2*(1 - 6*lmMt) - 33*lmMt + 6*lmMsq*(
        5 + 3*(lmMst2 + lmMt)) - 18*pow2(lmMsq))*pow4(Mst2)) + 15*pow2(Mst2)*(-
        ((1409 + 6*lmMst2*(89 - 84*lmMt) - 546*lmMt + 12*lmMsq*(1 + 42*(lmMst2
        + lmMt)) - 504*pow2(lmMsq))*pow2(Mst1)*pow2(Mst2)) + 4*(871 + 6*lmMst1*
        (-89 + 66*lmMsq - 66*lmMt) + 375*lmMt + lmMst2*(159 - 396*lmMsq + 396*
        lmMt))*pow4(Mst1) + (-299 + 78*lmMt - 12*lmMsq*(-5 + 6*lmMst2 + 6*lmMt)
        + 6*lmMst2*(-23 + 12*lmMt) + 72*pow2(lmMsq))*pow4(Mst2))))*pow4(Mt))/(
        243.*pow4(Msq)*pow4(Mst2)) - Mt*pow3(s2t)*(-(pow2(Mst1)*(
        110.18859670781893 - (40*B4)/9. - (2*DN)/9. - (80*lmMsq)/9. + (10*pow2(
        lmMsq))/3. + lmMst1*(189.21814814814815 - (355*lmMsq)/9. + (35*pow2(
        lmMsq))/3.) + (4.705555555555556 - (5*lmMsq)/3.)*pow2(lmMst1) - (
        lmMst2*(393289 - 397890*lmMst1 + 1500*lmMsq*(-59 + 36*lmMst1) + 31500*
        pow2(lmMsq) + 35850*pow2(lmMst1)))/2700. + ((-8151 + 1300*lmMsq + 3890*
        lmMst1)*pow2(lmMst2))/60. + (49*pow3(lmMst1))/54. - (2833*pow3(lmMst2))
        /54.)) + (Dmglst1*Mst1*(2323 + 144*B4 + 120*DN - 2880*lmMsq + 54*
        lmMst1*(32.55555555555556 + 40*lmMsq - 20*pow2(lmMsq)) + 6*lmMst2*(571
        - 360*lmMsq + 62*lmMst1 + 180*pow2(lmMsq) - 158*pow2(lmMst1)) + 30*(-
        121 + 36*lmMsq)*pow2(lmMst1) - 6*(-735 + 180*lmMsq + 98*lmMst1)*pow2(
        lmMst2) - 260*pow3(lmMst1) + 1796*pow3(lmMst2)))/54. + (-(50*Dmglst1*(
        192 - 91*lmMst1 + 6*lmMsq*(-6 + 17*lmMst1 - 17*lmMst2) + 127*lmMst2 -
        51*pow2(lmMst1) + 51*pow2(lmMst2)) + 3*Mst1*(787 + 20*lmMsq*(-9 + 20*
        lmMst1 - 20*lmMst2) + 525*lmMst2 - 5*lmMst1*(69 + 70*lmMst2) - 25*pow2(
        lmMst1) + 375*pow2(lmMst2)))/(135.*pow2(Msq)) - ((Dmglst1*(606133 -
        12960*lmMsq*(11 + 30*lmMst1 - 24*lmMst2) - 624756*lmMst2 + 38880*pow2(
        lmMsq) + 72*(3901 - 2976*lmMst2)*pow2(lmMst1) - 361944*pow2(lmMst2) +
        12*lmMst1*(92887 + 23460*lmMst2 + 36288*pow2(lmMst2)) - 2304*pow3(
        lmMst1) - 218880*pow3(lmMst2)))/3888. + Mst1*(239.07595982905215 - (
        71872687*lmMst2)/529200. + 5*pow2(lmMsq) + (-64.58505291005291 - (776*
        lmMst2)/9.)*pow2(lmMst1) - (770503*pow2(lmMst2))/7560. + (5*lmMsq*(-29
        + 18*lmMst2 - 2*lmMst1*(15 + 16*lmMst2) + 16*pow2(lmMst1) + 16*pow2(
        lmMst2)))/6. + (37*lmMst1*(2891251 + 2673860*lmMst2 + 2352000*pow2(
        lmMst2)))/529200. + (8*pow3(lmMst1))/3. - (728*pow3(lmMst2))/9.))/pow2(
        Mst2))*pow3(Mst1) - ((939 + 30*lmMsq*(-12 + 5*lmMst1 - 5*lmMst2) + 760*
        lmMst2 - 50*lmMst1*(8 + 3*lmMst2) + 150*pow2(lmMst2))/(135.*pow2(Msq))
        + (-1 + 2*lmMst2)/(3.*pow2(Mst1)) + ((5*Dmglst1*Mst1)/6. + (
        1.0841585681891803 - (157*lmMst1)/135. + lmMsq*(-0.6291383219954648 + (
        7*lmMst1)/9. - (43*lmMst2)/63.) + (47419*lmMst2)/26460. - (7*lmMst1*
        lmMst2)/9. - pow2(lmMsq)/21. + (46*pow2(lmMst2))/63.)*pow2(Mst1))/pow4(
        Msq))*pow4(Mst2) + (177*(8*OepS2 - 81*(5 + 2*lmMst1 - 2*lmMst2)*S2)*
        pow2(Mst1)*pow2(Mst2) + 2*(664*OepS2 - 27*(1187 + 498*lmMst1 - 498*
        lmMst2)*S2)*pow4(Mst1) + 27*(8*OepS2 + 81*(11 - 2*lmMst1 + 2*lmMst2)*
        S2)*pow4(Mst2))/(729.*pow2(Mst2)) + (pow2(Mst2)*(82320*(1775 + 30*
        lmMsq*(-5 + lmMst1 - lmMst2) + 407*lmMst2 + lmMst1*(-257 + 600*lmMst2)
        - 315*pow2(lmMst1) - 285*pow2(lmMst2))*pow2(Msq)*pow3(Mst1) + 17150*
        Mst1*(53965 + 2592*B4 - 288*DN - 13320*lmMsq + 5400*pow2(lmMsq) - 6*
        lmMst1*(3781 - 300*lmMsq + 180*pow2(lmMsq)) + 6*lmMst2*(14137 - 3498*
        lmMst1 + 60*lmMsq*(-35 + 12*lmMst1) + 180*pow2(lmMsq) - 18*pow2(lmMst1)
        ) - 18*(-229 + 60*lmMsq)*pow2(lmMst1) - 18*(-2193 + 180*lmMsq + 442*
        lmMst1)*pow2(lmMst2) + 396*pow3(lmMst1) + 7668*pow3(lmMst2))*pow4(Msq)
        + 102900*Dmglst1*(40*(82 - 93*lmMst1 + 6*lmMsq*(4 + 3*lmMst1 - 3*
        lmMst2) + 69*lmMst2 - 9*pow2(lmMst1) + 9*pow2(lmMst2))*pow2(Msq)*pow2(
        Mst1) + 6*(7 + 60*lmMsq*(5 - 6*lmMst1) - 2*lmMst1 + (-226 + 220*lmMst1)
        *lmMst2 + 180*pow2(lmMsq) + 178*pow2(lmMst1) + 18*pow2(lmMst2))*pow4(
        Msq) + 5*(503 - 346*lmMst1 + 12*lmMsq*(1 + 10*lmMst1 - 10*lmMst2) +
        334*lmMst2 - 60*pow2(lmMst1) + 60*pow2(lmMst2))*pow4(Mst1)) + 9*(
        6740969 + 140*lmMsq*(-12899 + 6230*lmMst1 - 5390*lmMst2) + 2822890*
        lmMst2 + 70*lmMst1*(-14529 + 25480*lmMst2) - 58800*pow2(lmMsq) -
        1327900*pow2(lmMst1) - 514500*pow2(lmMst2))*pow5(Mst1)))/(1.11132e7*
        Mst1*pow4(Msq))) - pow2(Mt)*pow2(s2t)*((-20*(Mst1*(71 - 16*lmMst1 + 12*
        lmMsq*(-2 + lmMst1 - lmMst2) + 40*lmMst2 - 6*pow2(lmMst1) + 6*pow2(
        lmMst2)) + 3*Dmglst1*(83 - 20*lmMst1 + 12*lmMsq*(-2 + lmMst1 - lmMst2)
        + 44*lmMst2 - 6*pow2(lmMst1) + 6*pow2(lmMst2)))*pow2(Mst1))/(9.*pow2(
        Msq)) + Dmglst1*(629.8333333333334 - (68*B4)/3. - (2*DN)/3. - 360*lmMsq
        + lmMst1*(121.66666666666667 + 60*lmMsq - 20*pow2(lmMsq)) + 40*pow2(
        lmMsq) + lmMst2*(777 - 140*lmMsq - (118*lmMst1)/3. + 20*pow2(lmMsq) -
        20*pow2(lmMst1)) + (-13 + 20*lmMsq)*pow2(lmMst1) - ((-797 + 60*lmMsq +
        4*lmMst1)*pow2(lmMst2))/3. - (100*pow3(lmMst1))/9. + (292*pow3(lmMst2))
        /9.) + (Mst1*(13317 - 408*B4 - 12*DN - 4320*lmMsq + 720*pow2(lmMsq) -
        6*lmMst1*(79 - 180*lmMsq + 60*pow2(lmMsq)) + 6*(-167 + 60*lmMsq)*pow2(
        lmMst1) - 6*lmMst2*(-2135 - 797*lmMst2 + 60*lmMsq*(7 + lmMst2) + 2*
        lmMst1*(59 + 2*lmMst2) - 60*pow2(lmMsq) + 60*pow2(lmMst1)) - 200*pow3(
        lmMst1) + 584*pow3(lmMst2)))/18. - (pow2(Mst1)*((Mst1*(6277 + 7776*B4 +
        71103*lmMst1 + 3240*(lmMst1 - lmMst2)*pow2(lmMsq) + 34902*pow2(lmMst1)
        - 3*lmMst2*(22549 + 2460*lmMst1 + 648*pow2(lmMst1)) - 18*(1625 + 276*
        lmMst1)*pow2(lmMst2) - 3240*lmMsq*(1 + lmMst1*(5 - 4*lmMst2) - 5*lmMst2
        + 3*pow2(lmMst1) + pow2(lmMst2)) + 8136*pow3(lmMst1) - 1224*pow3(
        lmMst2)))/81. + Dmglst1*(371.14814814814815 + (592*B4)/3. - (8*DN)/3. -
        280*lmMsq + lmMst1*(431.44444444444446 + 120*lmMsq + 40*pow2(lmMsq)) +
        (394 - 280*lmMsq)*pow2(lmMst1) - lmMst2*(120*lmMsq*(1 - 4*lmMst1) + 40*
        pow2(lmMsq) + (1195 + 4764*lmMst1 + 1788*pow2(lmMst1))/9.) + (
        71.33333333333333 - 200*lmMsq - 228*lmMst1)*pow2(lmMst2) + (2708*pow3(
        lmMst1))/9. + (1132*pow3(lmMst2))/9.)))/pow2(Mst2) + (5*pow2(Mst2)*(2*
        Dmglst1*(131 - 12*lmMsq + 12*lmMst2)*pow2(Mst1) + (-11 + 12*lmMsq - 12*
        lmMst2)*(Dmglst1 + Mst1)*pow2(Mst2) + 2*(59 - 12*lmMsq + 12*lmMst2)*
        pow3(Mst1)))/(108.*pow4(Msq)) + pow4(Mst1)*((-5*(Mst1*(911 - 54*lmMst1
        + 12*lmMsq*(-37 + 18*lmMst1 - 18*lmMst2) + 498*lmMst2 - 108*pow2(
        lmMst1) + 108*pow2(lmMst2)) + Dmglst1*(5159 - 702*lmMst1 + 12*lmMsq*(-
        181 + 90*lmMst1 - 90*lmMst2) + 2874*lmMst2 - 540*pow2(lmMst1) + 540*
        pow2(lmMst2))))/(108.*pow4(Msq)) - (Mst1*(96*B4 + 40*(lmMst1 - lmMst2)*
        pow2(lmMsq) + 10*lmMsq*(7 + 10*lmMst2 + 2*lmMst1*(-5 + 8*lmMst2) - 12*
        pow2(lmMst1) - 4*pow2(lmMst2)) + (240861 - 406036*lmMst2 + 24*(3257 +
        2628*lmMst2)*pow2(lmMst1) + 4*lmMst1*(94597 + 34932*lmMst2 - 15768*
        pow2(lmMst2)) - 217896*pow2(lmMst2) + 6624*pow3(lmMst1) - 6624*pow3(
        lmMst2))/432.) + Dmglst1*(2730.8603395061727 + (592*B4)/3. - (8*DN)/3.
         + 150*lmMsq + lmMst1*(
        554.0833333333334 + 780*lmMsq + 40*pow2(lmMsq)) - lmMst2*(60*lmMsq*(13
        - 8*lmMst1) + 40*pow2(lmMsq) + (31467 + 9076*lmMst1 - 6744*pow2(lmMst1)
        )/36.) - (5*(463 + 1008*lmMsq)*pow2(lmMst1))/18. - ((-6853 + 3600*lmMsq
        + 9516*lmMst1)*pow2(lmMst2))/18. + (1294*pow3(lmMst1))/9. + (1778*pow3(
        lmMst2))/9.))/pow4(Mst2)) + (20*(3*Mst1*pow4(Mst2) + 3*Dmglst1*(2*(12 -
        26*lmMst2 + 12*lmMsq*lmMst2 - 2*lmMst1*(-13 + 6*lmMsq + 15*lmMst2) +
        21*pow2(lmMst1) + 9*pow2(lmMst2))*pow4(Mst1) + pow4(Mst2)) + 2*(6 - 34*
        lmMst2 + 12*lmMsq*lmMst2 - 2*lmMst1*(-17 + 6*lmMsq + 9*lmMst2) + 15*
        pow2(lmMst1) + 3*pow2(lmMst2))*pow5(Mst1)))/(9.*pow2(Msq)*pow2(Mst2)))
        + Mt*s2t*(-((-1 + 2*lmMst2)*shiftst3*pow2(Mst2)*(-4*pow2(Mt) + ((-1 -
        lmMst1 + lmMst2)*pow2(Mst1) + pow2(Mst2))*pow2(s2t)))/(3.*pow2(Mst1)) +
        pow2(Mt)*(29.117283950617285 - (16*DN)/9. - (20*lmMsq)/3. + 20*pow2(
        lmMsq) - (4*lmMst1*(454 + 60*lmMsq + 45*pow2(lmMsq)))/27. + (4*lmMst2*(
        715 - 336*lmMst1 + 30*lmMsq*(-7 + 3*lmMst1) + 45*pow2(lmMsq) - 36*pow2(
        lmMst1)))/27. + (164*pow2(lmMst1))/9. - (4*(-191 + 30*lmMsq + 44*
        lmMst1)*pow2(lmMst2))/9. - (4*(2017 + 15*lmMst1*(782 - 375*lmMst2) +
        15*lmMsq*(-532 + 195*lmMst1 - 75*lmMst2) - 3750*lmMst2 - 900*pow2(
        lmMsq) + 1350*pow2(lmMst1) + 3375*pow2(lmMst2))*pow2(Mst1))/(2025.*
        pow2(Msq)) - (4*Dmglst1*(158 + 70*lmMst1 + 60*lmMsq*(-4 + 9*lmMst1) + (
        96 - 138*lmMst1)*lmMst2 - 270*pow2(lmMsq) - 459*pow2(lmMst1) - 27*pow2(
        lmMst2) + (30*(7 - 20*lmMsq + 26*lmMst1 - 6*lmMst2)*pow2(Mst1))/pow2(
        Msq)))/(27.*Mst1) - (116*pow3(lmMst1))/27. - (2*Mst1*(Mst1*(7874051 +
        7815060*lmMst2 + 303750*pow2(lmMsq) - 675*(989 + 290*lmMst2 - 480*lmMt)
        *pow2(lmMst1) - 90*lmMst1*(6359 - 10035*lmMst2 + 7200*(1 + lmMst2)*lmMt
        - 16575*pow2(lmMst2)) + 1978425*pow2(lmMst2) + 101250*lmMsq*(-9 - 6*
        lmMst2 - 4*lmMst1*lmMst2 + 2*pow2(lmMst1) + 2*pow2(lmMst2)) + 324000*
        lmMt*pow2(1 + lmMst2) - 582750*pow3(lmMst1) - 713250*pow3(lmMst2)) +
        125*Dmglst1*(164809 - 864*B4 - 432*DN - 14580*lmMsq - 9366*lmMst1 +
        4860*pow2(lmMsq) - 6*lmMst2*(-23575 + 1620*lmMsq + 906*lmMst1 - 432*
        pow2(lmMst1)) + 1854*pow2(lmMst1) + 18*(2167 + 336*lmMst1)*pow2(lmMst2)
        + 5184*lmMt*(2 + 3*lmMst2 - lmMst1*(3 + 2*lmMst2) + pow2(lmMst1) +
        pow2(lmMst2)) - 7200*pow3(lmMst1) - 1440*pow3(lmMst2))))/(30375.*pow2(
        Mst2)) + (788*pow3(lmMst2))/27. + pow2(Mst2)*((4 - 8*lmMst2)/(3.*pow2(
        Mst1)) + (51450*Dmglst1*(13 - 12*lmMsq + 12*lmMst2)*Mst1 + 1715*(195 -
        60*lmMsq*(3 + 2*lmMst1 - 2*lmMst2) + 4*lmMst2 + 8*lmMst1*(22 + 15*
        lmMst2) - 120*pow2(lmMst2))*pow2(Mst1) - (103583 + 420*lmMsq*(-256 +
        49*lmMst1 - 259*lmMst2) + 150052*lmMst2 - 1372*lmMst1*(31 + 15*lmMst2)
        + 44100*pow2(lmMsq) + 64680*pow2(lmMst2))*pow2(Mst2))/(92610.*pow4(Msq)
        )) + ((34300*Dmglst1*(11 + 144*lmMsq - 204*lmMst1 + 60*lmMst2) - Mst1*(
        187967 + 28*lmMst1*(49766 - 13965*lmMst2) + 420*lmMsq*(-2194 + 259*
        lmMst1 - 49*lmMst2) - 471968*lmMst2 - 44100*pow2(lmMsq) + 141120*pow2(
        lmMst1) + 205800*pow2(lmMst2)))*pow3(Mst1))/(92610.*pow4(Msq)) + S2*(-
        6*(1 + 2*lmMst1 - 2*lmMst2) - (28*(5 + 2*lmMst1 - 2*lmMst2)*pow2(Mst1))
        /(3.*pow2(Mst2)) - (8*(139 + 60*lmMst1 - 60*lmMst2)*pow4(Mst1))/(27.*
        pow4(Mst2))) - ((Mst1*(1311202751 + 43575800310*lmMst2 - 132300*(138494
        + 9555*lmMst2 - 15120*lmMt)*pow2(lmMst1) - 1890*lmMst1*(18939979 +
        1411200*lmMt + 280*lmMst2*(-36067 + 7560*lmMt) - 4196850*pow2(lmMst2))
        + 1681003800*pow2(lmMst2) + 416745000*lmMsq*(-2 + lmMst1*(5 - 2*lmMst2)
        - 5*lmMst2 + pow2(lmMst1) + pow2(lmMst2)) + 666792000*lmMt*(1 + 4*
        lmMst2 + 3*pow2(lmMst2)) - 3134848500*pow3(lmMst1) - 3533071500*pow3(
        lmMst2)) + 1543500*Dmglst1*(17783 - 144*B4 - 72*DN - 4860*lmMsq + 12*(-
        3889 + 450*lmMsq)*lmMst1 - 30483*pow2(lmMst1) + 6*lmMst2*(10610 - 900*
        lmMsq + 6897*lmMst1 + 228*pow2(lmMst1)) + 9*(-891 + 1160*lmMst1)*pow2(
        lmMst2) + 864*lmMt*(4 + 10*lmMst2 - 10*lmMst1*(1 + lmMst2) + 5*pow2(
        lmMst1) + 5*pow2(lmMst2)) - 6120*pow3(lmMst1) - 5688*pow3(lmMst2)))*
        pow3(Mst1))/(3.1255875e7*pow4(Mst2)) + (16*OepS2*(42*pow2(Mst1)*pow2(
        Mst2) + 40*pow4(Mst1) + 27*pow4(Mst2)))/(729.*pow4(Mst2)) + (4*(7500*
        Dmglst1*(11 - 6*lmMsq + 6*lmMst2)*pow3(Mst1) + 375*(55 - 30*lmMsq + 30*
        lmMst2 + 18*lmMst1*lmMst2 - 9*pow2(lmMst1) - 9*pow2(lmMst2))*pow4(Mst1)
        - (5108 + 15*lmMsq*(-218 + 75*lmMst1 - 195*lmMst2) + 6270*lmMst2 - 375*
        lmMst1*(8 + 3*lmMst2) + 900*pow2(lmMsq) + 2025*pow2(lmMst2))*pow4(Mst2)
        ))/(2025.*pow2(Msq)*pow2(Mst2))) - (4*T1ep*(2*(40*pow2(Mt) - 83*pow2(
        Mst2)*pow2(s2t))*pow4(Mst1) + 54*pow2(Mt)*pow4(Mst2) + pow2(Mst1)*(84*
        pow2(Mst2)*pow2(Mt) - 177*pow2(s2t)*pow4(Mst2)) - 27*pow2(s2t)*pow6(
        Mst2)))/(243.*pow4(Mst2)))) - (3240*xMsq*z2*pow2(Msq)*(-(pow2(s2t)*
        pow4(Mst1)*(8*(-shiftst1 + lmMst1*(-1 + 2*shiftst1 - shiftst2) +
        shiftst2 + lmMst2*(1 - 2*shiftst1 + shiftst2))*Tbeta*pow2(Mt)*pow2(
        MuSUSY)*(-1 + pow2(Sbeta)) - 4*Mt*(-1 + shiftst2)*(-(MuSUSY*s2t) + 2*
        Mt*Tbeta)*pow2(Mst2)*pow2(Sbeta) + (1 + 2*(lmMst1 - lmMst2)*(-1 +
        shiftst1) + 3*shiftst1 - 4*shiftst2)*Tbeta*pow2(s2t)*pow2(Sbeta)*pow4(
        Mst2))) + (1 - shiftst1)*pow4(Mst2)*(-4*Tbeta*pow2(Mt)*pow2(s2t)*(pow2(
        MuSUSY)*(-1 + pow2(Sbeta)) + 2*pow2(Mst2)*pow2(Sbeta)) + pow2(Sbeta)*(
        16*MuSUSY*s2t*pow3(Mt) - 4*Mt*MuSUSY*pow2(Mst2)*pow3(s2t) + 16*Tbeta*
        pow4(Mt) + Tbeta*pow4(Mst2)*pow4(s2t))) + pow2(Mst1)*pow2(Mst2)*(4*
        Tbeta*pow2(Mt)*pow2(s2t)*((1 + 2*(lmMst1 - lmMst2)*(-1 + shiftst1) - 2*
        shiftst1 + shiftst2)*pow2(MuSUSY)*(1 - pow2(Sbeta)) + 2*(2 + (-1 +
        lmMst1 - lmMst2)*shiftst1 + (-1 - lmMst1 + lmMst2)*shiftst2)*pow2(Mst2)
        *pow2(Sbeta)) - (-1 + shiftst2)*Tbeta*pow2(Sbeta)*pow4(Mst1)*pow4(s2t)
        - pow2(Sbeta)*(-16*MuSUSY*s2t*(-1 + shiftst2)*pow3(Mt) + 4*Mt*MuSUSY*(
        shiftst1 - shiftst2 + (lmMst1 - lmMst2)*(-2 + shiftst1 + shiftst2))*
        pow2(Mst2)*pow3(s2t) + 16*(-1 + shiftst2)*Tbeta*pow4(Mt) + (1 - 4*
        shiftst1 - 2*(lmMst1 - lmMst2)*(-1 + shiftst2) + 3*shiftst2)*Tbeta*
        pow4(Mst2)*pow4(s2t)))) + z3*(Mst1*Tbeta*pow2(Sbeta)*(-864*s2t*pow2(
        Mst1)*pow2(Mst2)*pow2(Mt)*(Dmglst1*(-373 - 24*lmMst1 + 24*lmMst2)*s2t*
        pow2(Mst2) + 2*(133 - 9*lmMst1 + 9*lmMst2)*Mt*pow2(MuSUSY)) - 432*pow4(
        Mt)*(8*pow2(Mst1)*(386*Dmglst1*pow2(Mst2) + 72*Mst1*pow2(Mst2) + 3*
        Dmglst1*(13 + 2*lmMst1 - 2*lmMst2)*pow2(MuSUSY) + (-23 + 3*lmMst1 - 3*
        lmMst2)*Mst1*pow2(MuSUSY)) + 5264*Dmglst1*pow4(Mst1) + (-107*Dmglst1 -
        4*(163 - 6*(lmMst1 + lmMst2) + 12*lmMt)*Mst1)*pow4(Mst2) + 1268*pow5(
        Mst1)) - 6*pow4(Mst2)*pow4(s2t)*(2*(27*Dmglst1*(135 + 92*lmMst1 - 92*
        lmMst2) + 7*(833 + 162*lmMst1 - 162*lmMst2)*Mst1)*pow2(Mst1)*pow2(Mst2)
        + 9*Dmglst1*(5077 - 552*lmMst1 + 552*lmMst2)*pow4(Mst1) - 9*(195*
        Dmglst1 - 328*Mst1)*pow4(Mst2) + 324*(39.78600823045267 - 7*lmMst1 + 7*
        lmMst2)*pow5(Mst1)) + 432*Mst1*Mt*s2t*(-(pow2(Mst2)*pow2(s2t)*(-6*(103
        + 75*lmMst1 - 75*lmMst2)*pow2(Mst2)*pow3(Mst1) + 2991*Dmglst1*pow4(
        Mst1) + 2*(-133 + 9*lmMst1 - 9*lmMst2)*Mst1*pow4(Mst2) - 2*Dmglst1*((
        346 + 693*lmMst1 - 693*lmMst2)*pow2(Mst1)*pow2(Mst2) + (58 - 9*lmMst1 +
        9*lmMst2)*pow4(Mst2)) + 519*pow5(Mst1))) + 4*pow2(Mt)*((366*pow2(Mst2)
        - 23*(25 + 9*lmMst1 - 9*lmMst2)*pow2(MuSUSY))*pow3(Mst1) + 12*(-9 + 2*
        lmMst1 - 2*lmMst2)*Mst1*pow4(Mst2) + Dmglst1*((-58 + 9*lmMst1 - 9*
        lmMst2)*pow2(Mst2)*pow2(MuSUSY) + 3*pow2(Mst1)*(338*pow2(Mst2) + (-154
        - 225*lmMst1 + 225*lmMst2)*pow2(MuSUSY)) + 4392*pow4(Mst1) + (-74 + 24*
        lmMst1 - 24*lmMst2)*pow4(Mst2)) + 896*pow5(Mst1))) + 8*pow2(Mt)*pow2(
        s2t)*(6*pow2(Mst2)*(323*pow2(Mst2) + (8783 + 1134*lmMst1 - 1134*lmMst2)
        *pow2(MuSUSY))*pow3(Mst1) - 27*Dmglst1*(-12*(35 + 46*lmMst1 - 46*
        lmMst2)*pow2(Mst1)*pow2(Mst2)*pow2(MuSUSY) + 8*(199*pow2(Mst2) + (-764
        - 69*lmMst1 + 69*lmMst2)*pow2(MuSUSY))*pow4(Mst1) + 5*(44*pow2(Mst2) +
        39*pow2(MuSUSY))*pow4(Mst2)) + 2*(-662*pow2(Mst2) + 14*(4829 + 243*
        lmMst1 - 243*lmMst2)*pow2(MuSUSY))*pow5(Mst1) + 54*Mst1*(164*pow2(
        MuSUSY)*pow4(Mst2) + 43*pow6(Mst2)))) + 8*Mst1*Mt*MuSUSY*(Mt*MuSUSY*
        Tbeta*(2*Mst1*(108*Mst1*Mt*s2t*(23*(25 + 9*lmMst1 - 9*lmMst2)*pow2(
        Mst1) + (133 - 9*lmMst1 + 9*lmMst2)*pow2(Mst2)) - pow2(Mst1)*(216*(23 -
        3*lmMst1 + 3*lmMst2)*pow2(Mt) + 3*(8783 + 1134*lmMst1 - 1134*lmMst2)*
        pow2(Mst2)*pow2(s2t)) + pow2(s2t)*(-14*(4829 + 243*lmMst1 - 243*lmMst2)
        *pow4(Mst1) - 4428*pow4(Mst2))) + 27*Dmglst1*(12*pow2(Mst1)*((52 + 8*
        lmMst1 - 8*lmMst2)*pow2(Mt) + (-35 - 46*lmMst1 + 46*lmMst2)*pow2(Mst2)*
        pow2(s2t)) + Mt*s2t*(8*(58 - 9*lmMst1 + 9*lmMst2)*Mst1*pow2(Mst2) + 24*
        (154 + 225*lmMst1 - 225*lmMst2)*pow3(Mst1)) + pow2(s2t)*(-8*(764 + 69*
        lmMst1 - 69*lmMst2)*pow4(Mst1) + 195*pow4(Mst2)))) + pow2(Sbeta)*(-2*
        Mst1*((2936*s2t*pow2(Mt) - 41257*pow2(Mst2)*pow3(s2t))*pow4(Mst1) + 81*
        Mst1*Mt*(4*pow2(Mst1)*((86 + 8*lmMst1 - 8*lmMst2)*pow2(Mt) + (221 +
        108*lmMst1 - 108*lmMst2)*pow2(Mst2)*pow2(s2t)) + 2*pow2(Mst2)*(8*(-9 +
        2*lmMst1 - 2*lmMst2)*pow2(Mt) + (133 - 9*lmMst1 + 9*lmMst2)*pow2(Mst2)*
        pow2(s2t)) + (365 + 432*lmMst1 - 432*lmMst2)*pow2(s2t)*pow4(Mst1)) +
        54*s2t*(43*pow2(Mt) - 82*pow2(Mst2)*pow2(s2t))*pow4(Mst2) + pow2(Mst1)*
        (4260*s2t*pow2(Mst2)*pow2(Mt) - 3*(7307 + 1134*lmMst1 - 1134*lmMst2)*
        pow3(s2t)*pow4(Mst2))) + 27*Dmglst1*(-16*pow3(Mst1)*(3*(101 + 171*
        lmMst1 - 171*lmMst2)*Mt*pow2(Mst2)*pow2(s2t) + 2*(235 + 6*lmMst1 - 6*
        lmMst2)*pow3(Mt)) + (-64*(-10 + 3*lmMst1 - 3*lmMst2)*s2t*pow2(Mt) +
        5692*pow2(Mst2)*pow3(s2t))*pow4(Mst1) + 440*s2t*pow2(Mt)*pow4(Mst2) +
        4*Mst1*(4*(37 - 12*lmMst1 + 12*lmMst2)*pow2(Mst2)*pow3(Mt) + 3*(-58 +
        9*lmMst1 - 9*lmMst2)*Mt*pow2(s2t)*pow4(Mst2)) + 3*pow2(Mst1)*(-16*(53 +
        4*lmMst1 - 4*lmMst2)*s2t*pow2(Mst2)*pow2(Mt) + (205 + 184*lmMst1 - 184*
        lmMst2)*pow3(s2t)*pow4(Mst2)) + 6*(2183 - 1368*lmMst1 + 1368*lmMst2)*
        Mt*pow2(s2t)*pow5(Mst1) - 195*pow3(s2t)*pow6(Mst2)))) + 27*xDmglst1*
        pow2(Dmglst1)*(pow4(Mst2)*(-12*Tbeta*pow2(Mt)*pow2(s2t)*(381*pow2(
        MuSUSY)*(-1 + pow2(Sbeta)) + 1013*pow2(Mst2)*pow2(Sbeta)) + pow2(Sbeta)
        *(24312*MuSUSY*s2t*pow3(Mt) - 4572*Mt*MuSUSY*pow2(Mst2)*pow3(s2t) +
        14440*Tbeta*pow4(Mt) + 1143*Tbeta*pow4(Mst2)*pow4(s2t))) + pow4(Mst1)*(
        -48*Tbeta*pow2(Mt)*pow2(s2t)*(-((1145 + 14*lmMst1 - 14*lmMst2)*pow2(
        MuSUSY)*(-1 + pow2(Sbeta))) + 1194*pow2(Mst2)*pow2(Sbeta)) + pow2(
        Sbeta)*(256*(134 - 3*lmMst1 + 3*lmMst2)*MuSUSY*s2t*pow3(Mt) + 48*(1117
        - 32*lmMst1 + 32*lmMst2)*Mt*MuSUSY*pow2(Mst2)*pow3(s2t) - 129664*Tbeta*
        pow4(Mt) + 9*(-1325 + 104*lmMst1 - 104*lmMst2)*Tbeta*pow4(Mst2)*pow4(
        s2t))) + Mt*(-16*(Mst1*pow2(Mst2)*(-2*s2t*Tbeta*pow2(Mt)*(675*pow2(
        MuSUSY)*(-1 + pow2(Sbeta)) + 1843*pow2(Mst2)*pow2(Sbeta)) + MuSUSY*
        pow2(Sbeta)*(-2025*Mt*pow2(Mst2)*pow2(s2t) + 3686*pow3(Mt)) + 675*
        Tbeta*pow2(Sbeta)*pow3(s2t)*pow4(Mst2)) + pow3(Mst1)*(MuSUSY*pow2(
        Sbeta)*(2106*(-1 + 2*lmMst1 - 2*lmMst2)*Mt*pow2(Mst2)*pow2(s2t) + 7190*
        pow3(Mt)) + Tbeta*(-6*s2t*pow2(Mt)*(-9*(-51 + 52*lmMst1 - 52*lmMst2)*
        pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 584*pow2(Mst2)*pow2(Sbeta)) + 27*(1 -
        52*lmMst1 + 52*lmMst2)*pow2(Sbeta)*pow3(s2t)*pow4(Mst2)))) + 32*s2t*(3*
        (3635 - 702*lmMst1 + 702*lmMst2)*Mt*MuSUSY*s2t + 4*Tbeta*(4222*pow2(Mt)
        - 821*pow2(Mst2)*pow2(s2t)))*pow2(Sbeta)*pow5(Mst1)) - 2*pow2(Mst1)*(
        16*pow2(Mst2)*pow2(Mt)*(-3*(14 + 23*lmMst1 - 23*lmMst2)*Tbeta*pow2(
        MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta)) + 2*Mt*((1255 + 12*lmMst1 - 12*
        lmMst2)*MuSUSY*s2t + 1554*Mt*Tbeta)*pow2(Sbeta)) - 2*Mt*(3*(493 + 184*
        lmMst1 - 184*lmMst2)*MuSUSY*s2t + (13079 + 96*lmMst1 - 96*lmMst2)*Mt*
        Tbeta)*pow2(s2t)*pow2(Sbeta)*pow4(Mst2) + Tbeta*(64*(181 + 3*lmMst1 -
        3*lmMst2)*pow2(MuSUSY)*(-1 + pow2(Sbeta))*pow4(Mt) + 69*(19 + 4*lmMst1
        - 4*lmMst2)*pow2(Sbeta)*pow4(s2t)*pow6(Mst2))))) + pow2(Mst1)*(24*pow2(
        z2)*(Mt*MuSUSY*pow2(Sbeta)*(576*(Dmglst1 + Mst1)*(pow2(Mst1) + pow2(
        Mst2))*pow3(Mt) + 2*s2t*pow2(Mt)*(42*pow2(Mst1)*pow2(Mst2) + 40*pow4(
        Mst1) - 45*pow4(Mst2)) - pow2(Mst2)*pow3(s2t)*(3*Mst1*(48*Dmglst1 +
        179*Mst1)*pow2(Mst2) + 166*pow4(Mst1) + 315*pow4(Mst2)) + 864*Mt*pow2(
        s2t)*(19*Dmglst1*pow2(Mst1)*(pow2(Mst1) + pow2(Mst2)) + 9*(pow2(Mst1) +
        pow2(Mst2))*pow3(Mst1) + 2*Dmglst1*pow4(Mst2) + 2*Mst1*pow4(Mst2))) -
        18*s2t*xDmglst1*pow2(Dmglst1)*(-4*MuSUSY*pow2(Mt)*(MuSUSY*Tbeta*(120*
        Mst1*Mt - s2t*(pow2(Mst1) + pow2(Mst2)))*pow2(Sbeta) + s2t*(MuSUSY*
        Tbeta*pow2(Mst2) + 180*Mst1*(pow2(Mst1) + pow2(Mst2))*pow2(Sbeta))) +
        Tbeta*(240*Mst1*Mt + s2t*pow2(Mst1) - s2t*pow2(Mst2))*pow2(s2t)*pow2(
        Sbeta)*pow4(Mst2) + 4*Mt*MuSUSY*(Mst1*Mt*MuSUSY*(120*Mt - Mst1*s2t)*
        Tbeta + pow2(s2t)*pow2(Sbeta)*pow4(Mst2))) + Tbeta*(s2t*pow2(Mt)*(-9*
        Mt*pow2(Sbeta)*(-64*pow2(MuSUSY)*(21*Dmglst1*pow2(Mst1) + 2*(Dmglst1 +
        Mst1)*pow2(Mst2) + 11*pow3(Mst1)) + 64*(Dmglst1 + Mst1)*pow4(Mst2)) -
        pow2(MuSUSY)*(576*Mt*(2*Mst1*pow2(Mst2) + 11*pow3(Mst1)) - 144*Dmglst1*
        (-84*Mt*pow2(Mst1) - 8*Mt*pow2(Mst2) + Mst1*s2t*pow2(Mst2) + s2t*pow3(
        Mst1)) - s2t*(852*pow2(Mst1)*pow2(Mst2) + 1018*pow4(Mst1) + 315*pow4(
        Mst2)))) + pow2(Sbeta)*(pow4(Mst2)*(-288*Mt*(17*Dmglst1*pow2(Mst1) + 2*
        (Dmglst1 + Mst1)*pow2(Mst2) + 7*pow3(Mst1))*pow3(s2t) - ((-222*pow2(
        Mst1)*pow2(Mst2) + 144*Dmglst1*(-(Mst1*pow2(Mst2)) + pow3(Mst1)) + 371*
        pow4(Mst1) - 315*pow4(Mst2))*pow4(s2t))/4.) + pow2(Mt)*pow2(s2t)*(2*
        pow2(Mst2)*pow4(Mst1) - 87*pow2(Mst1)*pow4(Mst2) - pow2(MuSUSY)*(852*
        pow2(Mst1)*pow2(Mst2) + 144*Dmglst1*Mst1*(pow2(Mst1) + pow2(Mst2)) +
        1018*pow4(Mst1) + 315*pow4(Mst2)) + 45*pow6(Mst2))))) - 4*z4*(2*Mt*
        MuSUSY*pow2(Sbeta)*(864*(Dmglst1 + Mst1)*(pow2(Mst1) + pow2(Mst2))*
        pow3(Mt) + 2*pow2(Mst2)*pow3(s2t)*(6*Mst1*(2385*Dmglst1 + 538*Mst1)*
        pow2(Mst2) + 166*pow4(Mst1) - 3267*pow4(Mst2)) - 4*s2t*pow2(Mt)*(42*
        pow2(Mst1)*pow2(Mst2) + 3240*Dmglst1*Mst1*(pow2(Mst1) + pow2(Mst2)) +
        40*pow4(Mst1) - 243*pow4(Mst2)) - 162*Mt*pow2(s2t)*(4*(187*Dmglst1 -
        18*Mst1)*pow2(Mst1)*pow2(Mst2) + 748*Dmglst1*pow4(Mst1) - 241*(Dmglst1
        + Mst1)*pow4(Mst2) - 72*pow5(Mst1))) - 135*xDmglst1*pow2(Dmglst1)*(4*
        MuSUSY*s2t*pow2(Mst2)*pow2(Mt)*(-53*MuSUSY*s2t*Tbeta*(-1 + pow2(Sbeta))
        + 24*Mt*pow2(Sbeta)) + 984*Mst1*Mt*s2t*(2*Tbeta*pow2(Mt)*pow2(MuSUSY)*(
        -1 + pow2(Sbeta)) + s2t*pow2(Mst2)*(3*Mt*MuSUSY - s2t*Tbeta*pow2(Mst2))
        *pow2(Sbeta)) - 4*Mt*(53*MuSUSY*s2t + 12*Mt*Tbeta)*pow2(s2t)*pow2(
        Sbeta)*pow4(Mst2) + s2t*pow2(Mst1)*(-20*s2t*Tbeta*pow2(Mt)*pow2(MuSUSY)
        *(-1 + pow2(Sbeta)) + 96*MuSUSY*pow2(Sbeta)*(2*Mt*pow2(Mst2)*pow2(s2t)
        + pow3(Mt)) - 101*Tbeta*pow2(Sbeta)*pow3(s2t)*pow4(Mst2)) + 48*Tbeta*
        pow2(MuSUSY)*(-1 + pow2(Sbeta))*pow4(Mt) + pow2(Sbeta)*(2952*MuSUSY*
        pow2(Mt)*pow2(s2t)*pow3(Mst1) + 53*Tbeta*pow4(s2t)*pow6(Mst2))) +
        Tbeta*(4*pow2(Mt)*pow2(MuSUSY)*(3*pow2(Mst1)*(540*pow2(Mt) + 13*pow2(
        Mst2)*pow2(s2t)) + 54*Dmglst1*(Mt*s2t*(507*pow2(Mst1) - 241*pow2(Mst2))
        + 60*Mst1*pow2(Mt) - 265*Mst1*(pow2(Mst1) + pow2(Mst2))*pow2(s2t)) -
        1620*Mst1*(2*Dmglst1 + Mst1)*pow2(Mt)*pow2(Sbeta) - 54*Mt*s2t*(241*
        Mst1*pow2(Mst2) + 313*pow3(Mst1)) - pow2(s2t)*(127*pow4(Mst1) - 3267*
        pow4(Mst2))) + s2t*pow2(Sbeta)*(-216*pow3(Mt)*(pow2(MuSUSY)*(507*
        Dmglst1*pow2(Mst1) - 241*(Dmglst1 + Mst1)*pow2(Mst2) - 313*pow3(Mst1))
        + 8*(Dmglst1 + Mst1)*pow4(Mst2)) + pow2(s2t)*pow4(Mst2)*(108*Mt*(989*
        Dmglst1*pow2(Mst1) - 241*(Dmglst1 + Mst1)*pow2(Mst2) + 169*pow3(Mst1))
        + s2t*(-6495*pow2(Mst1)*pow2(Mst2) + 14310*Dmglst1*(-(Mst1*pow2(Mst2))
        + pow3(Mst1)) + 3062*pow4(Mst1) + 3267*pow4(Mst2))) - 4*s2t*pow2(Mt)*(
        2*pow2(Mst2)*pow4(Mst1) - pow2(MuSUSY)*(-39*pow2(Mst1)*pow2(Mst2) +
        14310*Dmglst1*Mst1*(pow2(Mst1) + pow2(Mst2)) + 127*pow4(Mst1) - 3267*
        pow4(Mst2)) - 15*Mst1*(216*Dmglst1 + 19*Mst1)*pow4(Mst2) + 243*pow6(
        Mst2)))))))/(1944.*pow2(Mst1)*pow2(Sbeta)*pow4(Mst2)))/Tbeta + (pow2(
        s2t)*((-2*pow2(Mt)*(pow2(MuSUSY)*(81*(-1 + 3*lmMst1 - 3*lmMst2)*(-1 +
        2*lmMst2)*shiftst3 - 740*T1ep*(-1 + pow2(Sbeta))) + 4*T1ep*pow2(Mst2)*
        pow2(Sbeta)) + 11*T1ep*pow2(s2t)*pow2(Sbeta)*pow4(Mst2))*pow6(Mst1) +
        81*(-1 + 2*lmMst2)*shiftst3*pow2(Mt)*pow2(MuSUSY)*pow6(Mst2) - 6*pow4(
        Mst1)*(pow2(Mst2)*pow2(Mt)*(pow2(MuSUSY)*(27*shiftst3*(1 - 2*lmMst1 +
        4*lmMst1*lmMst2 - 4*pow2(lmMst2)) - 136*T1ep*(-1 + pow2(Sbeta))) - 10*
        T1ep*pow2(Mst2)*pow2(Sbeta)) + 25*T1ep*pow2(s2t)*pow2(Sbeta)*pow6(Mst2)
        ) - 27*pow2(Mst1)*(2*pow2(Mt)*(pow2(MuSUSY)*(3*(-1 + lmMst1 - lmMst2)*(
        -1 + 2*lmMst2)*shiftst3 - 2*T1ep*(-1 + pow2(Sbeta))) - 2*T1ep*pow2(
        Mst2)*pow2(Sbeta))*pow4(Mst2) + T1ep*pow2(s2t)*pow2(Sbeta)*pow8(Mst2)))
        )/(243.*pow2(Mst1)*pow2(Sbeta)*pow4(Mst2)) + ((1 - 2*lmMst2)*shiftst3*
        Tbeta*pow2(Sbeta)*(-((8*(1 - lmMst1 + lmMst2)*pow2(Mst1)*pow2(Mt)*pow2(
        s2t) + 16*pow4(Mt) + pow2(Mst1)*((3 + 2*lmMst1 - 2*lmMst2)*pow2(Mst1) -
        4*pow2(Mst2))*pow4(s2t))*pow6(Mst2)) + 4*pow2(Mt)*pow2(s2t)*(pow2(
        MuSUSY)*pow6(Mst2) + 2*(pow2(MuSUSY)*((1 - 2*lmMst1 + 2*lmMst2)*pow2(
        Mst2)*pow4(Mst1) + (1 - lmMst1 + lmMst2)*pow2(Mst1)*pow4(Mst2) + (1 -
        3*lmMst1 + 3*lmMst2)*pow6(Mst1)) + pow8(Mst2)))) + 10*(1 - 2*lmMsq)*
        xMsq*pow2(Msq)*(4*Mt*MuSUSY*s2t*(Mt*MuSUSY*s2t*Tbeta*pow2(Mst1)*(pow2(
        Mst2) - 2*(lmMst1 - lmMst2)*(pow2(Mst1) + pow2(Mst2))) + pow2(Mst2)*
        pow2(Sbeta)*((-1 + shiftst1)*pow2(Mst2)*(-4*pow2(Mt) + pow2(Mst2)*pow2(
        s2t)) - pow2(Mst1)*(-4*(-1 + shiftst2)*pow2(Mt) + (shiftst1 - shiftst2
        + (lmMst1 - lmMst2)*(-2 + shiftst1 + shiftst2))*pow2(Mst2)*pow2(s2t)) -
        (-1 + shiftst2)*pow2(s2t)*pow4(Mst1))) + Tbeta*(-8*pow2(Mst1)*pow2(
        Mst2)*pow2(Mt)*((1 - lmMst1 + lmMst2)*shiftst1*pow2(MuSUSY)*pow2(s2t) +
        2*shiftst2*pow2(Mt)*pow2(Sbeta)) + pow2(Sbeta)*(-4*pow2(Mt)*pow2(s2t)*(
        2*(pow2(Mst2) + (-lmMst1 + lmMst2)*pow2(MuSUSY))*pow4(Mst1) + pow2(
        Mst1)*((1 - 2*lmMst1 + 2*lmMst2)*pow2(Mst2)*pow2(MuSUSY) - 4*pow4(Mst2)
        ) + (2*pow2(Mst2) + pow2(MuSUSY))*pow4(Mst2)) + pow2(Mst2)*(16*(pow2(
        Mst1) + pow2(Mst2))*pow4(Mt) + pow2(Mst1)*((-1 + 2*lmMst1 - 2*lmMst2)*
        pow2(Mst1)*pow2(Mst2) + pow4(Mst1) + (-1 - 2*lmMst1 + 2*lmMst2)*pow4(
        Mst2))*pow4(s2t))) - shiftst1*(4*pow2(Mst1)*pow2(Mst2)*pow2(s2t)*pow2(
        Sbeta)*(2*(1 - lmMst1 + lmMst2)*pow2(Mt)*(pow2(Mst2) - pow2(MuSUSY)) -
        pow2(s2t)*pow4(Mst2)) + pow4(Mst1)*(8*(-1 + 2*lmMst1 - 2*lmMst2)*pow2(
        Mt)*pow2(MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta)) + (3 + 2*lmMst1 - 2*
        lmMst2)*pow2(Sbeta)*pow4(Mst2)*pow4(s2t)) + pow4(Mst2)*(-4*pow2(Mt)*
        pow2(s2t)*(pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 2*pow2(Mst2)*pow2(Sbeta))
        + pow2(Sbeta)*(16*pow4(Mt) + pow4(Mst2)*pow4(s2t)))) + pow2(s2t)*(4*
        pow2(Mt)*pow2(MuSUSY)*pow4(Mst2) - shiftst2*pow2(Mst1)*(4*pow2(Mst2)*
        pow2(Mt)*(pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 2*(1 + lmMst1 - lmMst2)*
        pow2(Mst2)*pow2(Sbeta)) - 4*pow2(Mst1)*(2*pow2(Mt)*((-1 + lmMst1 -
        lmMst2)*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + pow2(Mst2)*pow2(Sbeta)) +
        pow2(s2t)*pow2(Sbeta)*pow4(Mst2)) + pow2(s2t)*pow2(Sbeta)*(pow2(Mst2)*
        pow4(Mst1) + (3 - 2*lmMst1 + 2*lmMst2)*pow6(Mst2))) + pow2(s2t)*pow2(
        Sbeta)*pow8(Mst2)))))/(12.*Tbeta*pow2(Mst1)*pow2(Sbeta)*pow4(Mst2)));
}

/**
 * @return The off-diagonal (1, 2) = (2, 1) matrix element of the Higgs mass matrix as a double for the hierarchy 'H5'
 */
double H5::getS12() const {
   return -(MuSUSY*(12*threeLoopFlag*pow2(Al4p)*(-324*Mt*pow2(s2t)*((Dmglst1*Mst1*
        s2t*(2323 + 144*B4 + 120*DN - 2880*lmMsq - 6*lmMst1*(-293 - 360*lmMsq +
        180*pow2(lmMsq)) + 6*lmMst2*(571 - 360*lmMsq + 62*lmMst1 + 180*pow2(
        lmMsq) - 158*pow2(lmMst1)) + 30*(-121 + 36*lmMsq)*pow2(lmMst1) - 6*(-
        735 + 180*lmMsq + 98*lmMst1)*pow2(lmMst2) - 260*pow3(lmMst1) + 1796*
        pow3(lmMst2)))/108. + Mt*((10*(Dmglst1 + Mst1)*pow2(Mst2))/(3.*pow2(
        Msq)) + Dmglst1*(314.9166666666667 - (34*B4)/3. - DN/3. - 180*lmMsq +
        lmMst1*(60.833333333333336 + 30*lmMsq - 10*pow2(lmMsq)) + 20*pow2(
        lmMsq) + lmMst2*(388.5 - 70*lmMsq - (59*lmMst1)/3. + 10*pow2(lmMsq) -
        10*pow2(lmMst1)) + ((-13 + 20*lmMsq)*pow2(lmMst1))/2. - ((-797 + 60*
        lmMsq + 4*lmMst1)*pow2(lmMst2))/6. - (50*pow3(lmMst1))/9. + (146*pow3(
        lmMst2))/9.) + (Mst1*(13317 - 408*B4 - 12*DN - 4320*lmMsq + 720*pow2(
        lmMsq) - 6*lmMst1*(79 - 180*lmMsq + 60*pow2(lmMsq)) + 6*lmMst2*(2135 -
        420*lmMsq - 118*lmMst1 + 60*pow2(lmMsq) - 60*pow2(lmMst1)) + 6*(-167 +
        60*lmMsq)*pow2(lmMst1) - 6*(-797 + 60*lmMsq + 4*lmMst1)*pow2(lmMst2) -
        200*pow3(lmMst1) + 584*pow3(lmMst2)))/36.)) - 72*T1ep*(2*s2t*pow3(Mt) -
        Mt*pow2(Mst2)*pow3(s2t)) + s2t*(-(Mt*pow2(Mst2)*(216*(1 - 2*lmMst2)*
        shiftst3*(4*pow2(Mt) + (1 + lmMst1 - lmMst2)*pow2(Mst1)*pow2(s2t)) + (
        Mst1*pow2(s2t)*(Dmglst1*(36*(7 - 2*lmMst1 - 60*lmMsq*(-5 + 6*lmMst1) -
        226*lmMst2 + 220*lmMst1*lmMst2 + 180*pow2(lmMsq) + 178*pow2(lmMst1) +
        18*pow2(lmMst2))*pow2(Msq) + 240*(82 - 93*lmMst1 + 6*lmMsq*(4 + 3*
        lmMst1 - 3*lmMst2) + 69*lmMst2 - 9*pow2(lmMst1) + 9*pow2(lmMst2))*pow2(
        Mst1)) + Mst1*pow2(Msq)*(53965 + 2592*B4 - 288*DN - 13320*lmMsq + 192*
        OepS2 - 1944*(-11 + 2*lmMst1 - 2*lmMst2)*S2 + 5400*pow2(lmMsq) - 6*
        lmMst1*(3781 - 300*lmMsq + 180*pow2(lmMsq)) + 6*lmMst2*(14137 - 3498*
        lmMst1 + 60*lmMsq*(-35 + 12*lmMst1) + 180*pow2(lmMsq) - 18*pow2(lmMst1)
        ) - 18*(-229 + 60*lmMsq)*pow2(lmMst1) - 18*(-2193 + 180*lmMsq + 442*
        lmMst1)*pow2(lmMst2) + 396*pow3(lmMst1) + 7668*pow3(lmMst2))))/pow2(
        Msq)))/(4.*pow2(Mst1)) + 324*pow3(Mt)*(14.558641975308642 - (8*DN)/9. -
        (10*lmMsq)/3. + (8*OepS2)/27. - 3*(1 + 2*lmMst1 - 2*lmMst2)*S2 + 10*
        pow2(lmMsq) - (2*lmMst1*(454 + 60*lmMsq + 45*pow2(lmMsq)))/27. + (2*
        lmMst2*(715 - 336*lmMst1 + 30*lmMsq*(-7 + 3*lmMst1) + 45*pow2(lmMsq) -
        36*pow2(lmMst1)))/27. + (82*pow2(lmMst1))/9. - (2*(-191 + 30*lmMsq +
        44*lmMst1)*pow2(lmMst2))/9. - (58*pow3(lmMst1))/27. - Dmglst1*((2*(158
        + 60*lmMsq*(-4 + 9*lmMst1) + lmMst1*(70 - 138*lmMst2) + 96*lmMst2 -
        270*pow2(lmMsq) - 459*pow2(lmMst1) - 27*pow2(lmMst2)))/(27.*Mst1) +
        Mst1*((20*(7 - 20*lmMsq + 26*lmMst1 - 6*lmMst2))/(9.*pow2(Msq)) + (
        164809 - 864*B4 - 432*DN - 14580*lmMsq - 9366*lmMst1 + 4860*pow2(lmMsq)
        + 1854*pow2(lmMst1) + 6*lmMst2*(23575 - 1620*lmMsq - 906*lmMst1 + 432*
        pow2(lmMst1)) + 18*(2167 + 336*lmMst1)*pow2(lmMst2) + 5184*lmMt*(2 + 3*
        lmMst2 - lmMst1*(3 + 2*lmMst2) + pow2(lmMst1) + pow2(lmMst2)) - 7200*
        pow3(lmMst1) - 1440*pow3(lmMst2))/(243.*pow2(Mst2)))) + (394*pow3(
        lmMst2))/27. - (pow2(Mst2)*((4*(5108 + 15*lmMsq*(-218 + 75*lmMst1 -
        195*lmMst2) + 6270*lmMst2 - 375*lmMst1*(8 + 3*lmMst2) + 900*pow2(lmMsq)
        + 2025*pow2(lmMst2)))/pow2(Msq) + (2700*(-1 + 2*lmMst2))/pow2(Mst1) + (
        1125*Dmglst1*(-13 + 12*lmMsq - 12*lmMst2)*Mst1)/pow4(Msq)))/4050.)) + (
        9*Mt*(60*s2t*xMsq*(1 - 2*(lmMsq + z2))*pow2(Msq)*pow2(Mst2)*((-1 +
        shiftst1)*pow2(Mst2)*(-2*Mt*MuSUSY*s2t - 4*Tbeta*pow2(Mt) + Tbeta*pow2(
        Mst2)*pow2(s2t)) + pow2(Mst1)*(2*Mt*MuSUSY*s2t*(1 + 2*(lmMst1 - lmMst2)
        *(-1 + shiftst1) - 2*shiftst1 + shiftst2) + Tbeta*(4*(-1 + shiftst2)*
        pow2(Mt) - (shiftst1 - shiftst2 + (lmMst1 - lmMst2)*(-2 + shiftst1 +
        shiftst2))*pow2(Mst2)*pow2(s2t)))) + 2*pow2(Mst1)*(-2*Mst1*Mt*pow2(
        Mst2)*(4*Mt*MuSUSY*s2t*(241*z4 - 32*pow2(z2)) + 16*Tbeta*pow2(Mt)*(-z4
        + 2*pow2(z2)) + 3*Tbeta*pow2(Mst2)*pow2(s2t)*(-241*z4 + 32*pow2(z2))) +
        s2t*(Mt*MuSUSY*s2t*(484*z4 - 70*pow2(z2)) + 2*Tbeta*pow2(Mt)*(18*z4 +
        5*pow2(z2)) + Tbeta*pow2(Mst2)*pow2(s2t)*(-242*z4 + 35*pow2(z2)))*pow4(
        Mst2) + 2*Dmglst1*(-4*Mt*pow2(Mst2)*(Mst1*MuSUSY*pow2(s2t)*(265*z4 + 4*
        pow2(z2)) + Tbeta*pow2(Mt)*(-4*z4 + 8*pow2(z2)) + Mt*s2t*(241*MuSUSY*z4
        + 60*Mst1*Tbeta*z4 - 32*MuSUSY*pow2(z2))) + 240*Mst1*MuSUSY*z4*pow3(Mt)
        + Tbeta*pow2(s2t)*(723*Mt*z4 + 530*Mst1*s2t*z4 - 96*Mt*pow2(z2) + 8*
        Mst1*s2t*pow2(z2))*pow4(Mst2)) + 2*xDmglst1*pow2(Dmglst1)*(30*Mst1*Mt*
        s2t*(4*Mt*MuSUSY - 3*s2t*Tbeta*pow2(Mst2))*(41*z4 + 8*pow2(z2)) - 2*Mt*
        s2t*pow2(Mst2)*(265*MuSUSY*s2t*z4 + 60*Mt*Tbeta*z4 + 4*MuSUSY*s2t*pow2(
        z2)) + 120*MuSUSY*z4*pow3(Mt) + Tbeta*(265*z4 + 4*pow2(z2))*pow3(s2t)*
        pow4(Mst2))) + z3*(8*pow2(Mst1)*pow2(Mst2)*(s2t*pow2(Mst2)*(164*Mt*
        MuSUSY*s2t + 43*Tbeta*pow2(Mt) - 82*Tbeta*pow2(Mst2)*pow2(s2t)) + Mst1*
        Mt*(4*(-133 + 9*lmMst1 - 9*lmMst2)*Mt*MuSUSY*s2t + 24*(-9 + 2*lmMst1 -
        2*lmMst2)*Tbeta*pow2(Mt) + 3*(133 - 9*lmMst1 + 9*lmMst2)*Tbeta*pow2(
        Mst2)*pow2(s2t))) - 2*Dmglst1*Mst1*(4*Mst1*Mt*pow2(Mst2)*(4*(58 - 9*
        lmMst1 + 9*lmMst2)*Mt*MuSUSY*s2t + 4*(37 - 12*lmMst1 + 12*lmMst2)*
        Tbeta*pow2(Mt) + 3*(-58 + 9*lmMst1 - 9*lmMst2)*Tbeta*pow2(Mst2)*pow2(
        s2t)) + 5*s2t*(78*Mt*MuSUSY*s2t + 88*Tbeta*pow2(Mt) - 39*Tbeta*pow2(
        Mst2)*pow2(s2t))*pow4(Mst2) + 3*pow2(Mst1)*(-16*(53 + 4*lmMst1 - 4*
        lmMst2)*s2t*Tbeta*pow2(Mst2)*pow2(Mt) - 8*(35 + 46*lmMst1 - 46*lmMst2)*
        Mt*MuSUSY*pow2(Mst2)*pow2(s2t) + 32*(13 + 2*lmMst1 - 2*lmMst2)*MuSUSY*
        pow3(Mt) + (205 + 184*lmMst1 - 184*lmMst2)*Tbeta*pow3(s2t)*pow4(Mst2)))
        - xDmglst1*pow2(Dmglst1)*(Mt*(4*Mst1*pow2(Mst2)*(-2700*Mt*MuSUSY*s2t -
        3686*Tbeta*pow2(Mt) + 2025*Tbeta*pow2(Mst2)*pow2(s2t)) - 8*(54*(51 -
        52*lmMst1 + 52*lmMst2)*Mt*MuSUSY*s2t + 3595*Tbeta*pow2(Mt) + 1053*(-1 +
        2*lmMst1 - 2*lmMst2)*Tbeta*pow2(Mst2)*pow2(s2t))*pow3(Mst1)) + 3*s2t*(
        762*Mt*MuSUSY*s2t + 2026*Tbeta*pow2(Mt) - 381*Tbeta*pow2(Mst2)*pow2(
        s2t))*pow4(Mst2) + pow2(Mst1)*(-16*(1255 + 12*lmMst1 - 12*lmMst2)*s2t*
        Tbeta*pow2(Mst2)*pow2(Mt) - 48*(14 + 23*lmMst1 - 23*lmMst2)*Mt*MuSUSY*
        pow2(Mst2)*pow2(s2t) + 64*(181 + 3*lmMst1 - 3*lmMst2)*MuSUSY*pow3(Mt) +
        3*(493 + 184*lmMst1 - 184*lmMst2)*Tbeta*pow3(s2t)*pow4(Mst2))))))/(
        Tbeta*pow2(Mst1)*pow4(Mst2)) + 2*((MuSUSY*pow2(Mt)*((pow2(s2t)*(270905
        + 12960*B4 - 66600*lmMsq + 960*OepS2 + 3240*(33 - 6*lmMst1 + 6*lmMst2)*
        S2 - 1440*(DN + T1ep) + 27000*pow2(lmMsq) - 30*lmMst1*(3781 - 300*lmMsq
        + 180*pow2(lmMsq)) + 30*lmMst2*(14065 - 3498*lmMst1 + 60*lmMsq*(-35 +
        12*lmMst1) + 180*pow2(lmMsq) - 18*pow2(lmMst1)) + 90*(229 - 60*lmMsq)*
        pow2(lmMst1) - 90*(-2193 + 180*lmMsq + 442*lmMst1)*pow2(lmMst2) + (
        1080*(-1 + 2*lmMst2)*shiftst3*((2 - 2*lmMst1 + 2*lmMst2)*pow2(Mst1) +
        pow2(Mst2)))/pow2(Mst1) + 1980*pow3(lmMst1) + 38340*pow3(lmMst2) + (60*
        Dmglst1*(3*(7 - 60*lmMsq*(-5 + 6*lmMst1) - 226*lmMst2 + lmMst1*(-2 +
        220*lmMst2) + 180*pow2(lmMsq) + 178*pow2(lmMst1) + 18*pow2(lmMst2)) +
        pow2(Mst1)*((20*(82 - 93*lmMst1 + 6*lmMsq*(4 + 3*lmMst1 - 3*lmMst2) +
        69*lmMst2 - 9*pow2(lmMst1) + 9*pow2(lmMst2)))/pow2(Msq) + (4*(586 + 36*
        B4 + 30*DN - 495*lmMsq + 135*pow2(lmMsq) - 6*lmMst1*(-73 - 45*lmMsq +
        45*pow2(lmMsq)) + 3*lmMst2*(229 - 180*lmMsq + 86*lmMst1 + 90*pow2(
        lmMsq) - 79*pow2(lmMst1)) + 18*(-43 + 15*lmMsq)*pow2(lmMst1) - 3*(-372
        + 90*lmMsq + 49*lmMst1)*pow2(lmMst2) - 65*pow3(lmMst1) + 449*pow3(
        lmMst2)))/pow2(Mst2))))/Mst1 - 12*pow2(Mst2)*((2*(939 + 30*lmMsq*(-12 +
        5*lmMst1 - 5*lmMst2) + 760*lmMst2 - 50*lmMst1*(8 + 3*lmMst2) + 150*
        pow2(lmMst2)))/pow2(Msq) + (90*(-1 + 2*lmMst2))/pow2(Mst1) + (225*
        Dmglst1*Mst1)/pow4(Msq))))/20. + Mt*((720*Dmglst1*s2t)/pow2(Msq) + (
        720*Mst1*s2t)/pow2(Msq) + (6*s2t*(Dmglst1*(11337 - 408*B4 - 12*DN -
        6480*lmMsq + 720*pow2(lmMsq) - 30*lmMst1*(-73 - 36*lmMsq + 12*pow2(
        lmMsq)) + 6*lmMst2*(2331 - 420*lmMsq - 118*lmMst1 + 60*pow2(lmMsq) -
        60*pow2(lmMst1)) + 18*(-13 + 20*lmMsq)*pow2(lmMst1) - 6*(-797 + 60*
        lmMsq + 4*lmMst1)*pow2(lmMst2) - 200*pow3(lmMst1) + 584*pow3(lmMst2)) +
        Mst1*(13317 - 408*B4 - 12*DN - 4320*lmMsq + 720*pow2(lmMsq) - 6*lmMst1*
        (79 - 180*lmMsq + 60*pow2(lmMsq)) + 6*lmMst2*(2135 - 420*lmMsq - 118*
        lmMst1 + 60*pow2(lmMsq) - 60*pow2(lmMst1)) + 6*(-167 + 60*lmMsq)*pow2(
        lmMst1) - 6*(-797 + 60*lmMsq + 4*lmMst1)*pow2(lmMst2) - 200*pow3(
        lmMst1) + 584*pow3(lmMst2))))/pow2(Mst2) + (5*(-11 + 12*lmMsq - 12*
        lmMst2)*(Dmglst1 + Mst1)*s2t*pow2(Mst2))/pow4(Msq) + (96*Dmglst1*Mst1*
        Mt*(735 - 6*B4 - 3*DN - 78*lmMst1 + 6*pow2(lmMst1) + 6*lmMst2*(85 - 6*
        lmMst1 + pow2(lmMst1)) + 18*(7 + lmMst1)*pow2(lmMst2) + 24*lmMt*(2 + 3*
        lmMst2 - lmMst1*(3 + 2*lmMst2) + pow2(lmMst1) + pow2(lmMst2)) - 26*
        pow3(lmMst1) + 2*pow3(lmMst2)))/pow4(Mst2))))/Tbeta + ((80*Dmglst1*(17
        + 3*lmMst2*(1 - 6*lmMt) - 33*lmMt + 6*lmMsq*(5 + 3*(lmMst2 + lmMt)) -
        18*pow2(lmMsq))*pow2(Msq)*pow2(Mst2) + 80*Mst1*(17 + 3*lmMst2*(1 - 6*
        lmMt) - 33*lmMt + 6*lmMsq*(5 + 3*(lmMst2 + lmMt)) - 18*pow2(lmMsq))*
        pow2(Msq)*pow2(Mst2) - 8*Dmglst1*(1846 - 1057*lmMst2 + 30*lmMsq*(151 +
        51*lmMt + 9*lmMst2*(3 + 2*lmMt) - 6*lmMst1*(13 + 3*lmMt)) + 270*(lmMst1
        - lmMst2)*pow2(lmMsq) + 6*(352 + 15*lmMst2 + 42*lmMt)*pow2(lmMst1) -
        306*pow2(lmMst2) - 3*lmMt*(1105 + 608*lmMst2 + 108*pow2(lmMst2)) - 432*
        (2 + lmMst2)*pow2(lmMt) + 2*lmMst1*(-2534 - 192*lmMt + 3*lmMst2*(67 +
        12*lmMt) + 9*pow2(lmMst2) + 216*pow2(lmMt)) - 66*pow3(lmMst1) - 42*
        pow3(lmMst2))*pow4(Msq) + 24*Mst1*(212 + 225*lmMst2 - 90*lmMsq*(7 + 3*
        lmMt + lmMst2*(3 - lmMsq + 2*lmMt) + lmMst1*(lmMsq - 2*(3 + lmMt))) -
        6*(67 + 5*lmMst2 + 14*lmMt)*pow2(lmMst1) + 102*pow2(lmMst2) + 3*lmMt*(
        265 + 202*lmMst2 + 36*pow2(lmMst2)) + 144*(1 + lmMst2)*pow2(lmMt) - 6*
        lmMst1*(-54 + 37*lmMt + lmMst2*(22 + 4*lmMt) + pow2(lmMst2) + 24*pow2(
        lmMt)) + 22*pow3(lmMst1) + 14*pow3(lmMst2))*pow4(Msq) + 5*Dmglst1*(299
        - 60*lmMsq + 6*lmMst2*(23 - 12*lmMt) - 78*lmMt + 72*lmMsq*(lmMst2 +
        lmMt) - 72*pow2(lmMsq))*pow4(Mst2) + 5*Mst1*(299 - 60*lmMsq + 6*lmMst2*
        (23 - 12*lmMt) - 78*lmMt + 72*lmMsq*(lmMst2 + lmMt) - 72*pow2(lmMsq))*
        pow4(Mst2))*pow4(Mt))/(pow2(Mst2)*pow4(Msq))) + 324*z2*(s2t*((Mt*pow2(
        Mst2)*((24*shiftst3*pow2(Mt))/pow2(Mst1) + (333 - 30*lmMsq - 37*lmMst1
        + 149*lmMst2 - (198*Dmglst1)/Mst1 + 6*shiftst3 + 6*lmMst1*shiftst3 - 6*
        lmMst2*shiftst3 - (60*Dmglst1*Mst1)/pow2(Msq))*pow2(s2t)))/18. + (
        32.611111111111114 + (4*lmMst2)/3. + (44*Dmglst1)/Mst1 - (4*Dmglst1*(
        2777 + 120*lmMst1 + 456*lmMst2)*Mst1)/(27.*pow2(Mst2)) - (4*(lmMst1 +
        pow2(Mst2)/pow2(Mst1)))/3.)*pow3(Mt)) + (Dmglst1*(-681 + 180*lmMsq +
        278*lmMst1 - 630*lmMst2)*Mst1*Mt*pow3(s2t))/9. + (MuSUSY*pow2(Mt)*(-((
        339 - 30*lmMsq - 37*lmMst1 + 149*lmMst2 - (198*Dmglst1)/Mst1 + 12*
        lmMst1*shiftst3 - 12*(1 + lmMst2)*shiftst3 + Dmglst1*Mst1*(-60/pow2(
        Msq) + (4*(-390 + 90*lmMsq + 139*lmMst1 - 315*lmMst2))/pow2(Mst2)) + (
        6*pow2(Mst2))/pow2(Mst1) - (6*shiftst3*pow2(Mst2))/pow2(Mst1))*pow2(
        s2t)) + (4*Mt*(8*Dmglst1*(55 + 16*lmMst2)*Mst1*Mt + Dmglst1*(627 - 60*
        lmMsq - 64*lmMst1 + 192*lmMst2)*s2t*pow2(Mst2) + (627 - 60*lmMsq - 64*
        lmMst1 + 192*lmMst2)*Mst1*s2t*pow2(Mst2)))/pow4(Mst2)))/(9.*Tbeta) + (
        Dmglst1 + Mst1)*(((-627 + 60*lmMsq + 64*lmMst1 - 192*lmMst2)*pow2(Mt)*
        pow2(s2t))/3. + (8*(-60*pow2(Msq)*pow2(Mst2) + (8 - 96*lmMst1 + 54*
        lmMst2 + 42*lmMt)*pow4(Msq) - 15*pow4(Mst2))*pow4(Mt))/(27.*pow2(Mst2)*
        pow4(Msq))))) + Mt*((243*oneLoopFlag*s2t*(2*(2 - lmMst1 + lmMst2)*Mt*
        MuSUSY*s2t + 4*(-lmMst1 + lmMst2)*Tbeta*pow2(Mt) + (-2 + lmMst1 -
        lmMst2)*Tbeta*pow2(Mst2)*pow2(s2t)))/Tbeta + Al4p*(1296*twoLoopFlag*(-
        2*s2t*(2*lmMst1*(1 + lmMst2) - lmMst2*(4 + 3*lmMst2) + (Dmglst1*(2 - 4*
        lmMst1))/Mst1 + pow2(lmMst1) + (4*Dmglst1*(1 + lmMst2)*Mst1)/pow2(Mst2)
        )*pow2(Mt) + (3*Dmglst1*Mt*(-12 + 2*lmMst1 - 6*lmMst2 + pow2(lmMst1) -
        pow2(lmMst2)) + 3*Mst1*Mt*(-8 + 2*lmMst1 - 6*lmMst2 + pow2(lmMst1) -
        pow2(lmMst2)) + 2*Dmglst1*Mst1*s2t*(-1 + lmMst1 - lmMst2 + pow2(lmMst1)
        - pow2(lmMst2)))*pow2(s2t) + ((2*Mt*MuSUSY*s2t*(2*Dmglst1*Mst1*Mt*(12 -
        2*lmMst1 + 6*lmMst2 - pow2(lmMst1) + pow2(lmMst2)) + 2*Mt*(8 - 2*lmMst1
        + 6*lmMst2 - pow2(lmMst1) + pow2(lmMst2))*pow2(Mst1) + Mst1*s2t*((7 -
        2*lmMst1)*(1 + lmMst2) + 2*pow2(lmMst2))*pow2(Mst2) - 2*Dmglst1*s2t*((-
        (lmMst2*(1 + lmMst2)) + pow2(lmMst1))*pow2(Mst1) - (-1 + lmMst1)*pow2(
        Mst2))))/(Mst1*Tbeta) + 8*(Dmglst1*(5 + lmMst2 + (2 + lmMst2)*lmMt -
        lmMst1*(3 + lmMt)) + (2 + lmMst2 + lmMt + lmMst2*lmMt - lmMst1*(2 +
        lmMt))*Mst1)*pow3(Mt))/pow2(Mst2) + ((-2*Dmglst1*(-1 + lmMst1) - Mst1*(
        (7 - 2*lmMst1)*(1 + lmMst2) + 2*pow2(lmMst2)))*pow2(Mst2)*pow3(s2t))/
        Mst1) - 72*z2*(18*s2t*twoLoopFlag*(6*(Dmglst1 + Mst1)*Mt*s2t + (2*Mt*
        MuSUSY*s2t)/Tbeta - (8*Mt*MuSUSY*(Mst1*Mt + Dmglst1*(Mt + Mst1*s2t)))/(
        Tbeta*pow2(Mst2)) + (4*Dmglst1*Mst1 - pow2(Mst2))*pow2(s2t)) +
        xDmglst1*pow2(Dmglst1)*((36*twoLoopFlag*(pow2(Mst2)*(-2*Mt*MuSUSY +
        s2t*Tbeta*pow2(Mst2))*pow2(s2t) + 6*Mst1*Mt*(4*Mt*MuSUSY*s2t + 4*Tbeta*
        pow2(Mt) - 3*Tbeta*pow2(Mst2)*pow2(s2t))))/(Tbeta*pow4(Mst2)) + Al4p*
        threeLoopFlag*((1440*Mst1*pow3(Mt))/pow4(Msq) + (6*Mt*pow2(Mst2)*(2*(-
        390 + 90*lmMsq + 139*lmMst1 - 315*lmMst2)*MuSUSY + 6*(-1187 + 360*lmMsq
        + 171*lmMst1 - 747*lmMst2)*Mst1*Tbeta - (30*(MuSUSY + 12*Mst1*Tbeta)*
        pow2(Mst2))/pow2(Msq) - (99*MuSUSY*pow2(Mst2))/pow2(Mst1))*pow2(s2t) -
        48*(2*(55 + 16*lmMst2)*MuSUSY + (-1391 + 360*lmMsq + 66*lmMst1 - 813*
        lmMst2 - 21*lmMt)*Mst1*Tbeta)*pow3(Mt) + 3*Tbeta*(681 - 180*lmMsq -
        278*lmMst1 + 630*lmMst2 + (30*pow2(Mst2))/pow2(Msq) + (99*pow2(Mst2))/
        pow2(Mst1))*pow3(s2t)*pow4(Mst2) + (4*s2t*pow2(Mt)*((2777 + 120*lmMst1
        + 456*lmMst2)*Tbeta*pow2(Mst1)*pow2(Mst2) + 12*MuSUSY*(1187 - 360*lmMsq
        - 171*lmMst1 + 747*lmMst2 + (60*pow2(Mst2))/pow2(Msq))*pow3(Mst1) -
        297*Tbeta*pow4(Mst2)))/pow2(Mst1))/(Tbeta*pow4(Mst2))))) + 2*xDmglst1*
        pow2(Dmglst1)*(Al4p*threeLoopFlag*(9*(5355 - 144*B4 - 120*DN + 9960*
        lmMsq + 2*lmMst1*(-4841 - 1080*lmMsq + 540*pow2(lmMsq)) - 6*lmMst2*(925
        - 360*lmMsq + 62*lmMst1 + 180*pow2(lmMsq) - 158*pow2(lmMst1)) - 6*(-221
        + 180*lmMsq)*pow2(lmMst1) + 6*(-735 + 180*lmMsq + 98*lmMst1)*pow2(
        lmMst2) + ((-20*(160 - 213*lmMst1 + 18*lmMsq*(8 + lmMst1 - lmMst2) +
        69*lmMst2 - 9*pow2(lmMst1) + 9*pow2(lmMst2)))/pow2(Msq) + (-6415 + 970*
        lmMst1 + 60*lmMsq*(-37 + 18*lmMst1) + 1314*lmMst2 - 660*lmMst1*lmMst2 -
        540*pow2(lmMsq) + 234*pow2(lmMst1) - 54*pow2(lmMst2))/pow2(Mst1))*pow2(
        Mst2) + 260*pow3(lmMst1) - 1796*pow3(lmMst2))*pow3(s2t) + pow2(Mt)*((-
        720*(4*(245 + 6*lmMst1*(-25 + 6*lmMsq - 6*lmMt) + 84*lmMt + lmMst2*(66
        - 36*lmMsq + 36*lmMt))*Mst1*Mt - 3*(5 + 60*lmMsq - 66*lmMst1 + 6*
        lmMst2)*s2t*pow2(Mst2)))/(pow2(Msq)*pow2(Mst2)) + 2*s2t*((36*(4717 -
        540*lmMsq*(-2 + lmMst1) - 252*lmMst2 + 2*lmMst1*(-661 + 69*lmMst2) +
        270*pow2(lmMsq) + 75*pow2(lmMst1) + 27*pow2(lmMst2)))/pow2(Mst1) - (2*(
        416017 - 864*B4 - 432*DN - 14580*lmMsq + 45858*lmMst1 + 4860*pow2(
        lmMsq) + 7902*pow2(lmMst1) + 6*lmMst2*(21571 - 1620*lmMsq - 618*lmMst1
        + 432*pow2(lmMst1)) + 18*(1735 + 336*lmMst1)*pow2(lmMst2) + 1728*lmMt*(
        17*(1 + lmMst2) - lmMst1*(17 + 6*lmMst2) + 3*(pow2(lmMst1) + pow2(
        lmMst2))) - 7200*pow3(lmMst1) - 1440*pow3(lmMst2)))/pow2(Mst2) + (135*(
        13 - 12*lmMsq + 12*lmMst2)*pow2(Mst2))/pow4(Msq))) + Mt*(1944*pow2(s2t)
        *((5875 + 1020*lmMsq - 478*lmMst1 - 610*lmMst2 - 576*pow2(lmMst1))/(18.
        *Mst1) + Mst1*((10*(311 - 66*lmMst1 + 36*lmMsq*(-2 + lmMst1 - lmMst2) +
        138*lmMst2 - 18*pow2(lmMst1) + 18*pow2(lmMst2)))/(9.*pow2(Msq)) + (
        354.6111111111111 + 76*B4 - 2*DN - 460*lmMsq + (-140.77777777777777 +
        400*lmMsq)*lmMst1 + lmMst2*(802.1111111111111 - (374*lmMst1)/3. + 80*
        lmMsq*(-5 + 3*lmMst1) - 131*pow2(lmMst1)) - ((901 + 360*lmMsq)*pow2(
        lmMst1))/3. + (393 - 120*lmMsq - 125*lmMst1)*pow2(lmMst2) + (451*pow3(
        lmMst1))/3. + (317*pow3(lmMst2))/3.)/pow2(Mst2) - (5*pow2(Mst2))/pow4(
        Msq))) + (6*MuSUSY*(3*pow2(s2t)*((6415 - 60*lmMsq*(-37 + 18*lmMst1) -
        1314*lmMst2 + 10*lmMst1*(-97 + 66*lmMst2) + 540*pow2(lmMsq) - 234*pow2(
        lmMst1) + 54*pow2(lmMst2))/pow2(Mst1) + (4*(265 + 36*B4 + 30*DN - 1935*
        lmMsq + 18*lmMst1*(121 + 15*lmMsq - 15*pow2(lmMsq)) + 135*pow2(lmMsq) +
        3*lmMst2*(353 - 180*lmMsq + 86*lmMst1 + 90*pow2(lmMsq) - 79*pow2(
        lmMst1)) + 30*(-13 + 9*lmMsq)*pow2(lmMst1) - 3*(-372 + 90*lmMsq + 49*
        lmMst1)*pow2(lmMst2) - 65*pow3(lmMst1) + 449*pow3(lmMst2)))/pow2(Mst2)
        - (45*pow2(Mst2))/pow4(Msq)) + 12*s2t*((5*s2t*(160 - 213*lmMst1 + 18*
        lmMsq*(8 + lmMst1 - lmMst2) + 69*lmMst2 - 9*pow2(lmMst1) + 9*pow2(
        lmMst2)))/pow2(Msq) - (2*Mt*(5875 + 1020*lmMsq - 478*lmMst1 - 610*
        lmMst2 - 576*pow2(lmMst1)))/(Mst1*pow2(Mst2)) - (40*Mst1*Mt*(311 - 66*
        lmMst1 + 36*lmMsq*(-2 + lmMst1 - lmMst2) + 138*lmMst2 - 18*pow2(lmMst1)
        + 18*pow2(lmMst2)))/(pow2(Msq)*pow2(Mst2)) + 2*Mst1*Mt*(90/pow4(Msq) -
        (2*(6129 + 684*B4 - 18*DN - 3630*lmMsq + 6*(-251 + 600*lmMsq)*lmMst1 +
        lmMst2*(6914 - 1122*lmMst1 + 720*lmMsq*(-5 + 3*lmMst1) - 1179*pow2(
        lmMst1)) - 3*(997 + 360*lmMsq)*pow2(lmMst1) - 9*(-393 + 120*lmMsq +
        125*lmMst1)*pow2(lmMst2) + 1353*pow3(lmMst1) + 951*pow3(lmMst2)))/pow4(
        Mst2))) + (32*pow2(Mt)*(5483 - 18*B4 - 9*DN + 418*lmMst1 + 102*pow2(
        lmMst1) + 2*lmMst2*(739 - 42*lmMst1 + 9*pow2(lmMst1)) + 54*(5 + lmMst1)
        *pow2(lmMst2) + 24*lmMt*(17 + 17*lmMst2 - lmMst1*(17 + 6*lmMst2) + 3*
        pow2(lmMst1) + 3*pow2(lmMst2)) - 78*pow3(lmMst1) + 6*pow3(lmMst2)))/
        pow4(Mst2)))/Tbeta) + (24*pow3(Mt)*((15*(185 + lmMst2*(66 - 72*lmMt) -
        78*lmMt + 12*lmMsq*(1 + 6*(lmMst2 + lmMt)) - 72*pow2(lmMsq))*pow2(Mst1)
        )/pow4(Msq) - (2*((30251 - 804*lmMst2 + 2*lmMst1*(-5032 + 27*lmMst2 -
        789*lmMt) + (852 - 54*lmMst2)*lmMt + 180*lmMsq*(31 - 6*lmMst1 + 6*lmMt)
        + 570*pow2(lmMst1) - 720*pow2(lmMt))*pow2(Mst2) + pow2(Mst1)*(32085 +
        23606*lmMst2 + (3606 - 3402*lmMst2 + 7506*lmMt)*pow2(lmMst1) + 180*
        lmMsq*(31 + 6*lmMst1*(-1 + 6*lmMst2 - 6*lmMt) + 84*lmMt + 6*lmMst2*(-13
        + 6*lmMt) - 36*pow2(lmMst2)) + 18846*pow2(lmMst2) - 6*lmMt*(2850 +
        1348*lmMst2 + 189*pow2(lmMst2)) - 2*lmMst1*(5093 + 5700*lmMt + 6*
        lmMst2*(-173 + 531*lmMt) + 162*pow2(lmMst2) - 1296*pow2(lmMt)) - 144*(
        35 + 18*lmMst2)*pow2(lmMt) - 1368*pow3(lmMst1) + 5094*pow3(lmMst2))))/
        pow4(Mst2)))/Mst1) + (72*twoLoopFlag*(2*pow3(Mst1)*(27*Mt*(6 + 7*lmMst2
        - lmMst1*(7 + 6*lmMst2) + 3*(pow2(lmMst1) + pow2(lmMst2)))*pow2(Mst2)*
        pow2(s2t) + 2*(61 + 3*lmMst1*(-5 + 18*lmMst2 - 18*lmMt) + 105*lmMt +
        18*lmMst2*(-5 + 3*lmMt) - 54*pow2(lmMst2))*pow3(Mt)) + 6*(-11 + 6*
        lmMst1)*s2t*pow2(Mt)*pow4(Mst2) - (18*Mt*MuSUSY*s2t*(-8*Mst1*Mt*pow2(
        Mst2) - s2t*(4 + lmMst2 - pow2(lmMst1) + pow2(lmMst2))*pow2(Mst1)*pow2(
        Mst2) + 4*Mt*(4 + 7*lmMst2 - lmMst1*(7 + 6*lmMst2) + 3*(pow2(lmMst1) +
        pow2(lmMst2)))*pow3(Mst1) - (-2 + lmMst1)*s2t*pow4(Mst2)))/Tbeta + 4*
        Mst1*((61 - 15*lmMst1 + 15*lmMt)*pow2(Mst2)*pow3(Mt) - 27*Mt*pow2(s2t)*
        pow4(Mst2)) - 9*pow2(Mst1)*(4*(1 + lmMst2)*s2t*pow2(Mst2)*pow2(Mt) + (6
        - lmMst1 + lmMst2 - pow2(lmMst1) + pow2(lmMst2))*pow3(s2t)*pow4(Mst2))
        - 9*(-2 + lmMst1)*pow3(s2t)*pow6(Mst2)))/(pow2(Mst1)*pow4(Mst2)))))))/
        3888.;
}

/**
 * @return returns the susy log^0 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double H5::calc_coef_at_as2_no_sm_logs_log0() const {

   const double result =
      1555.7037037037037 - (17720*Dmglst1)/(27.*Mst1) + (112112*Dmglst1*s2t)/(
        27.*Mt) + (5936*Mst1*s2t)/(9.*Mt) + (6272*z2)/9. + (1056*Dmglst1*z2)/
        Mst1 + (832*Dmglst1*s2t*z2)/(9.*Mt) + (832*Mst1*s2t*z2)/(9.*Mt) - (
        5216*z3)/3. - (856*Dmglst1*z3)/(3.*Mst1) + (2368*Dmglst1*s2t*z3)/(3.*
        Mt) + (1152*Mst1*s2t*z3)/Mt - (128*Dmglst1*s2t*z4)/(3.*Mt) - (128*Mst1*
        s2t*z4)/(3.*Mt) + (477200*s2t*pow2(Dmglst1))/(27.*Mst1*Mt) - (29488*
        s2t*z3*pow2(Dmglst1))/(3.*Mst1*Mt) - (1440*Dmglst1*Mst1)/pow2(Msq) + (
        1280*Dmglst1*Mst1*z2)/pow2(Msq) - (42368*pow2(Dmglst1))/(45.*pow2(Msq))
        + (64960*Mst1*s2t*pow2(Dmglst1))/(9.*Mt*pow2(Msq)) + (640*z2*pow2(
        Dmglst1))/pow2(Msq) + (1413352*pow2(Dmglst1))/(675.*pow2(Mst1)) + (528*
        z2*pow2(Dmglst1))/pow2(Mst1) - (7220*z3*pow2(Dmglst1))/(3.*pow2(Mst1))
        + (160*pow2(Msq))/pow2(Mst1) - (320*z2*pow2(Msq))/pow2(Mst1) - (21936*
        pow2(Mst1))/(25.*pow2(Msq)) + (40960*Dmglst1*s2t*pow2(Mst1))/(9.*Mt*
        pow2(Msq)) + (2560*z2*pow2(Mst1))/(3.*pow2(Msq)) - (1280*Dmglst1*s2t*
        z2*pow2(Mst1))/(3.*Mt*pow2(Msq)) - (1025488*Dmglst1*Mst1)/(81.*pow2(
        Mst2)) + (19328*Dmglst1*Mst1*z2)/(9.*pow2(Mst2)) + (24704*Dmglst1*Mst1*
        z3)/(3.*pow2(Mst2)) - (45265856*pow2(Dmglst1))/(2025.*pow2(Mst2)) + (
        172960*Mst1*s2t*pow2(Dmglst1))/(27.*Mt*pow2(Mst2)) + (9664*z2*pow2(
        Dmglst1))/(9.*pow2(Mst2)) + (88352*Mst1*s2t*z2*pow2(Dmglst1))/(3.*Mt*
        pow2(Mst2)) + (16576*z3*pow2(Dmglst1))/pow2(Mst2) - (9344*Mst1*s2t*z3*
        pow2(Dmglst1))/(Mt*pow2(Mst2)) + (160*pow2(Msq))/pow2(Mst2) - (320*z2*
        pow2(Msq))/pow2(Mst2) - (39353408*pow2(Mst1))/(10125.*pow2(Mst2)) + (
        395536*Dmglst1*s2t*pow2(Mst1))/(81.*Mt*pow2(Mst2)) + (9248*z2*pow2(
        Mst1))/(9.*pow2(Mst2)) + (265600*Dmglst1*s2t*z2*pow2(Mst1))/(9.*Mt*
        pow2(Mst2)) + (1536*z3*pow2(Mst1))/pow2(Mst2) - (10816*Dmglst1*s2t*z3*
        pow2(Mst1))/(Mt*pow2(Mst2)) - (8800*pow2(Dmglst1)*pow2(Mst1))/(3.*pow2(
        Msq)*pow2(Mst2)) - (216272*pow2(Mst2))/(675.*pow2(Msq)) - (10720*
        Dmglst1*s2t*pow2(Mst2))/(27.*Mt*pow2(Msq)) - (10720*Mst1*s2t*pow2(Mst2)
        )/(27.*Mt*pow2(Msq)) + (640*z2*pow2(Mst2))/(3.*pow2(Msq)) + (1280*
        Dmglst1*s2t*z2*pow2(Mst2))/(3.*Mt*pow2(Msq)) + (1280*Mst1*s2t*z2*pow2(
        Mst2))/(3.*Mt*pow2(Msq)) + (16*pow2(Mst2))/pow2(Mst1) - (32*z2*pow2(
        Mst2))/pow2(Mst1) + (627124*Dmglst1*Mst1*pow2(s2t))/(81.*pow2(Mt)) - (
        128*B4*Dmglst1*Mst1*pow2(s2t))/(3.*pow2(Mt)) - (64*Dmglst1*DN*Mst1*
        pow2(s2t))/(3.*pow2(Mt)) + (49184*Dmglst1*Mst1*z2*pow2(s2t))/(9.*pow2(
        Mt)) - (5968*Dmglst1*Mst1*z3*pow2(s2t))/(3.*pow2(Mt)) + (320*Dmglst1*
        Mst1*z4*pow2(s2t))/pow2(Mt) + (972470*pow2(Dmglst1)*pow2(s2t))/(81.*
        pow2(Mt)) - (64*B4*pow2(Dmglst1)*pow2(s2t))/(3.*pow2(Mt)) - (32*DN*
        pow2(Dmglst1)*pow2(s2t))/(3.*pow2(Mt)) + (24592*z2*pow2(Dmglst1)*pow2(
        s2t))/(9.*pow2(Mt)) - (26158*z3*pow2(Dmglst1)*pow2(s2t))/(3.*pow2(Mt))
        + (160*z4*pow2(Dmglst1)*pow2(s2t))/pow2(Mt) + (160*pow2(Msq)*pow2(s2t))
        /pow2(Mt) - (320*z2*pow2(Msq)*pow2(s2t))/pow2(Mt) + (32617079*pow2(
        Mst1)*pow2(s2t))/(10125.*pow2(Mt)) - (32*DN*pow2(Mst1)*pow2(s2t))/(3.*
        pow2(Mt)) - (160*OepS2*pow2(Mst1)*pow2(s2t))/(81.*pow2(Mt)) + (244*S2*
        pow2(Mst1)*pow2(s2t))/pow2(Mt) + (80*T1ep*pow2(Mst1)*pow2(s2t))/(27.*
        pow2(Mt)) + (80570*z2*pow2(Mst1)*pow2(s2t))/(27.*pow2(Mt)) - (2584*z3*
        pow2(Mst1)*pow2(s2t))/(27.*pow2(Mt)) + (760*z4*pow2(Mst1)*pow2(s2t))/(
        27.*pow2(Mt)) - (1400*pow2(Dmglst1)*pow2(Mst1)*pow2(s2t))/(pow2(Msq)*
        pow2(Mt)) - (5794*pow2(Dmglst1)*pow2(Mst1)*pow2(s2t))/(27.*pow2(Mst2)*
        pow2(Mt)) - (29776*z2*pow2(Dmglst1)*pow2(Mst1)*pow2(s2t))/(3.*pow2(
        Mst2)*pow2(Mt)) + (9552*z3*pow2(Dmglst1)*pow2(Mst1)*pow2(s2t))/(pow2(
        Mst2)*pow2(Mt)) - (80*pow2(Msq)*pow2(Mst1)*pow2(s2t))/(pow2(Mst2)*pow2(
        Mt)) + (160*z2*pow2(Msq)*pow2(Mst1)*pow2(s2t))/(pow2(Mst2)*pow2(Mt)) -
        (4501*pow2(Mst2)*pow2(s2t))/(27.*pow2(Mt)) + (32*DN*pow2(Mst2)*pow2(
        s2t))/(3.*pow2(Mt)) + (1264*Dmglst1*pow2(Mst2)*pow2(s2t))/(9.*Mst1*
        pow2(Mt)) - (32*OepS2*pow2(Mst2)*pow2(s2t))/(9.*pow2(Mt)) + (36*S2*
        pow2(Mst2)*pow2(s2t))/pow2(Mt) + (16*T1ep*pow2(Mst2)*pow2(s2t))/(3.*
        pow2(Mt)) - (1222*z2*pow2(Mst2)*pow2(s2t))/(3.*pow2(Mt)) - (528*
        Dmglst1*z2*pow2(Mst2)*pow2(s2t))/(Mst1*pow2(Mt)) - (344*z3*pow2(Mst2)*
        pow2(s2t))/(3.*pow2(Mt)) + (880*Dmglst1*z3*pow2(Mst2)*pow2(s2t))/(3.*
        Mst1*pow2(Mt)) - (24*z4*pow2(Mst2)*pow2(s2t))/pow2(Mt) + (560*Dmglst1*
        Mst1*pow2(Mst2)*pow2(s2t))/(3.*pow2(Msq)*pow2(Mt)) - (200*pow2(Dmglst1)
        *pow2(Mst2)*pow2(s2t))/(3.*pow2(Msq)*pow2(Mt)) - (18868*pow2(Dmglst1)*
        pow2(Mst2)*pow2(s2t))/(9.*pow2(Mst1)*pow2(Mt)) - (264*z2*pow2(Dmglst1)*
        pow2(Mst2)*pow2(s2t))/(pow2(Mst1)*pow2(Mt)) + (2026*z3*pow2(Dmglst1)*
        pow2(Mst2)*pow2(s2t))/(pow2(Mst1)*pow2(Mt)) - (80*pow2(Msq)*pow2(Mst2)*
        pow2(s2t))/(pow2(Mst1)*pow2(Mt)) + (160*z2*pow2(Msq)*pow2(Mst2)*pow2(
        s2t))/(pow2(Mst1)*pow2(Mt)) - (24728*pow2(Mst1)*pow2(Mst2)*pow2(s2t))/(
        675.*pow2(Msq)*pow2(Mt)) + (256*Dmglst1*s2t*pow2(z2))/(3.*Mt) + (256*
        Mst1*s2t*pow2(z2))/(3.*Mt) + (116*pow2(Mst1)*pow2(s2t)*pow2(z2))/(9.*
        pow2(Mt)) - (20*pow2(Mst2)*pow2(s2t)*pow2(z2))/(3.*pow2(Mt)) + (35200*
        s2t*pow3(Mst1))/(27.*Mt*pow2(Msq)) - (1280*s2t*z2*pow3(Mst1))/(3.*Mt*
        pow2(Msq)) + (1328*s2t*pow3(Mst1))/(27.*Mt*pow2(Mst2)) + (29632*s2t*z2*
        pow3(Mst1))/(3.*Mt*pow2(Mst2)) - (3904*s2t*z3*pow3(Mst1))/(Mt*pow2(
        Mst2)) - (17600*Dmglst1*pow3(Mst1))/(9.*pow2(Msq)*pow2(Mst2)) - (11840*
        s2t*pow2(Dmglst1)*pow3(Mst1))/(3.*Mt*pow2(Msq)*pow2(Mst2)) - (12800*
        s2t*z2*pow2(Dmglst1)*pow3(Mst1))/(Mt*pow2(Msq)*pow2(Mst2)) - (10480*
        Dmglst1*pow2(s2t)*pow3(Mst1))/(9.*pow2(Msq)*pow2(Mt)) - (253180*
        Dmglst1*pow2(s2t)*pow3(Mst1))/(81.*pow2(Mst2)*pow2(Mt)) - (29536*
        Dmglst1*z2*pow2(s2t)*pow3(Mst1))/(9.*pow2(Mst2)*pow2(Mt)) + (6368*
        Dmglst1*z3*pow2(s2t)*pow3(Mst1))/(3.*pow2(Mst2)*pow2(Mt)) - (2032*Mst1*
        pow2(Dmglst1)*pow3(s2t))/(9.*pow3(Mt)) - (608*B4*Mst1*pow2(Dmglst1)*
        pow3(s2t))/pow3(Mt) + (16*DN*Mst1*pow2(Dmglst1)*pow3(s2t))/pow3(Mt) - (
        18992*Mst1*z2*pow2(Dmglst1)*pow3(s2t))/(3.*pow3(Mt)) + (72*Mst1*z3*
        pow2(Dmglst1)*pow3(s2t))/pow3(Mt) + (3280*Mst1*z4*pow2(Dmglst1)*pow3(
        s2t))/pow3(Mt) - (108106*Dmglst1*pow2(Mst1)*pow3(s2t))/(27.*pow3(Mt)) -
        (2096*B4*Dmglst1*pow2(Mst1)*pow3(s2t))/(3.*pow3(Mt)) + (40*Dmglst1*DN*
        pow2(Mst1)*pow3(s2t))/(3.*pow3(Mt)) - (5176*Dmglst1*z2*pow2(Mst1)*pow3(
        s2t))/pow3(Mt) - (5536*Dmglst1*z3*pow2(Mst1)*pow3(s2t))/(3.*pow3(Mt)) +
        (7912*Dmglst1*z4*pow2(Mst1)*pow3(s2t))/(3.*pow3(Mt)) + (7558*Dmglst1*
        pow2(Mst2)*pow3(s2t))/(3.*pow3(Mt)) - (272*B4*Dmglst1*pow2(Mst2)*pow3(
        s2t))/(3.*pow3(Mt)) - (8*Dmglst1*DN*pow2(Mst2)*pow3(s2t))/(3.*pow3(Mt))
        + (8878*Mst1*pow2(Mst2)*pow3(s2t))/(3.*pow3(Mt)) - (272*B4*Mst1*pow2(
        Mst2)*pow3(s2t))/(3.*pow3(Mt)) - (8*DN*Mst1*pow2(Mst2)*pow3(s2t))/(3.*
        pow3(Mt)) + (1672*Dmglst1*z2*pow2(Mst2)*pow3(s2t))/pow3(Mt) + (1672*
        Mst1*z2*pow2(Mst2)*pow3(s2t))/pow3(Mt) - (928*Dmglst1*z3*pow2(Mst2)*
        pow3(s2t))/(3.*pow3(Mt)) - (2128*Mst1*z3*pow2(Mst2)*pow3(s2t))/(3.*
        pow3(Mt)) - (1928*Dmglst1*z4*pow2(Mst2)*pow3(s2t))/(3.*pow3(Mt)) - (
        1928*Mst1*z4*pow2(Mst2)*pow3(s2t))/(3.*pow3(Mt)) - (23500*pow2(Dmglst1)
        *pow2(Mst2)*pow3(s2t))/(9.*Mst1*pow3(Mt)) + (1800*z3*pow2(Dmglst1)*
        pow2(Mst2)*pow3(s2t))/(Mst1*pow3(Mt)) - (24880*Mst1*pow2(Dmglst1)*pow2(
        Mst2)*pow3(s2t))/(9.*pow2(Msq)*pow3(Mt)) - (320*Mst1*z2*pow2(Dmglst1)*
        pow2(Mst2)*pow3(s2t))/(pow2(Msq)*pow3(Mt)) - (2240*Dmglst1*pow2(Mst1)*
        pow2(Mst2)*pow3(s2t))/(pow2(Msq)*pow3(Mt)) - (320*Dmglst1*z2*pow2(Mst1)
        *pow2(Mst2)*pow3(s2t))/(pow2(Msq)*pow3(Mt)) + (640*Mst1*pow2(Dmglst1)*
        pow2(z2)*pow3(s2t))/pow3(Mt) + (2176*Dmglst1*pow2(Mst1)*pow2(z2)*pow3(
        s2t))/(3.*pow3(Mt)) + (256*Dmglst1*pow2(Mst2)*pow2(z2)*pow3(s2t))/(3.*
        pow3(Mt)) + (256*Mst1*pow2(Mst2)*pow2(z2)*pow3(s2t))/(3.*pow3(Mt)) - (
        264814*pow3(Mst1)*pow3(s2t))/(81.*pow3(Mt)) - (880*B4*pow3(Mst1)*pow3(
        s2t))/(3.*pow3(Mt)) + (8*DN*pow3(Mst1)*pow3(s2t))/(3.*pow3(Mt)) - (
        10840*z2*pow3(Mst1)*pow3(s2t))/(9.*pow3(Mt)) - (1648*z3*pow3(Mst1)*
        pow3(s2t))/pow3(Mt) + (1352*z4*pow3(Mst1)*pow3(s2t))/(3.*pow3(Mt)) + (
        47440*pow2(Dmglst1)*pow3(Mst1)*pow3(s2t))/(9.*pow2(Msq)*pow3(Mt)) + (
        2880*z2*pow2(Dmglst1)*pow3(Mst1)*pow3(s2t))/(pow2(Msq)*pow3(Mt)) - (
        3124381*pow2(Dmglst1)*pow3(Mst1)*pow3(s2t))/(162.*pow2(Mst2)*pow3(Mt))
        - (248*z2*pow2(Dmglst1)*pow3(Mst1)*pow3(s2t))/(9.*pow2(Mst2)*pow3(Mt))
        + (52544*z3*pow2(Dmglst1)*pow3(Mst1)*pow3(s2t))/(3.*pow2(Mst2)*pow3(Mt)
        ) - (5920*pow2(Mst2)*pow3(Mst1)*pow3(s2t))/(9.*pow2(Msq)*pow3(Mt)) - (
        320*z2*pow2(Mst2)*pow3(Mst1)*pow3(s2t))/(3.*pow2(Msq)*pow3(Mt)) + (896*
        pow2(z2)*pow3(Mst1)*pow3(s2t))/(3.*pow3(Mt)) + (160*pow3(log(pow2(Mst1)
        /pow2(Msq))))/3. - (73742*pow2(Dmglst1)*pow2(Mst1))/(45.*pow4(Msq)) + (
        1760*z2*pow2(Dmglst1)*pow2(Mst1))/pow4(Msq) - (4600*Dmglst1*Mst1*pow2(
        Mst2))/(9.*pow4(Msq)) + (320*Dmglst1*Mst1*z2*pow2(Mst2))/pow4(Msq) - (
        2300*pow2(Dmglst1)*pow2(Mst2))/(9.*pow4(Msq)) - (8960*Mst1*s2t*pow2(
        Dmglst1)*pow2(Mst2))/(9.*Mt*pow4(Msq)) + (160*z2*pow2(Dmglst1)*pow2(
        Mst2))/pow4(Msq) + (640*Mst1*s2t*z2*pow2(Dmglst1)*pow2(Mst2))/(Mt*pow4(
        Msq)) - (2300*pow2(Mst1)*pow2(Mst2))/(9.*pow4(Msq)) - (8960*Dmglst1*
        s2t*pow2(Mst1)*pow2(Mst2))/(9.*Mt*pow4(Msq)) + (160*z2*pow2(Mst1)*pow2(
        Mst2))/pow4(Msq) + (640*Dmglst1*s2t*z2*pow2(Mst1)*pow2(Mst2))/(Mt*pow4(
        Msq)) - (295*pow2(Dmglst1)*pow2(Mst1)*pow2(Mst2)*pow2(s2t))/(pow2(Mt)*
        pow4(Msq)) - (41900*Dmglst1*pow3(Mst1))/(27.*pow4(Msq)) + (4160*
        Dmglst1*z2*pow3(Mst1))/(3.*pow4(Msq)) + (59360*s2t*pow2(Dmglst1)*pow3(
        Mst1))/(9.*Mt*pow4(Msq)) - (640*s2t*z2*pow2(Dmglst1)*pow3(Mst1))/(Mt*
        pow4(Msq)) - (8960*s2t*pow2(Mst2)*pow3(Mst1))/(27.*Mt*pow4(Msq)) + (
        640*s2t*z2*pow2(Mst2)*pow3(Mst1))/(3.*Mt*pow4(Msq)) + (170*Dmglst1*
        pow2(Mst2)*pow2(s2t)*pow3(Mst1))/(9.*pow2(Mt)*pow4(Msq)) - (20380*pow2(
        Dmglst1)*pow2(Mst2)*pow3(Mst1)*pow3(s2t))/(9.*pow3(Mt)*pow4(Msq)) - (
        400*z2*pow2(Dmglst1)*pow2(Mst2)*pow3(Mst1)*pow3(s2t))/(pow3(Mt)*pow4(
        Msq)) - (4400*pow4(Mst1))/(9.*pow2(Msq)*pow2(Mst2)) - (3520*Dmglst1*
        s2t*pow4(Mst1))/(3.*Mt*pow2(Msq)*pow2(Mst2)) - (6400*Dmglst1*s2t*z2*
        pow4(Mst1))/(Mt*pow2(Msq)*pow2(Mst2)) - (181136*pow2(s2t)*pow4(Mst1))/(
        675.*pow2(Msq)*pow2(Mt)) - (29787188414*pow2(s2t)*pow4(Mst1))/(
        1.0418625e7*pow2(Mst2)*pow2(Mt)) + (64*OepS2*pow2(s2t)*pow4(Mst1))/(
        243.*pow2(Mst2)*pow2(Mt)) - (296*S2*pow2(s2t)*pow4(Mst1))/(9.*pow2(
        Mst2)*pow2(Mt)) - (32*T1ep*pow2(s2t)*pow4(Mst1))/(81.*pow2(Mst2)*pow2(
        Mt)) + (1180*z2*pow2(s2t)*pow4(Mst1))/(81.*pow2(Mst2)*pow2(Mt)) + (
        5296*z3*pow2(s2t)*pow4(Mst1))/(81.*pow2(Mst2)*pow2(Mt)) - (16*z4*pow2(
        s2t)*pow4(Mst1))/(81.*pow2(Mst2)*pow2(Mt)) - (8*pow2(s2t)*pow2(z2)*
        pow4(Mst1))/(27.*pow2(Mst2)*pow2(Mt)) + (8560*Dmglst1*pow3(s2t)*pow4(
        Mst1))/(3.*pow2(Msq)*pow3(Mt)) + (1280*Dmglst1*z2*pow3(s2t)*pow4(Mst1))
        /(pow2(Msq)*pow3(Mt)) - (3058187*Dmglst1*pow3(s2t)*pow4(Mst1))/(324.*
        pow2(Mst2)*pow3(Mt)) + (5588*Dmglst1*z2*pow3(s2t)*pow4(Mst1))/(9.*pow2(
        Mst2)*pow3(Mt)) + (7976*Dmglst1*z3*pow3(s2t)*pow4(Mst1))/(pow2(Mst2)*
        pow3(Mt)) + (1930*Cbeta*MuSUSY*pow2(Dmglst1)*pow3(s2t)*pow4(Mst1))/(3.*
        Sbeta*pow2(Msq)*pow2(Mst2)*pow3(Mt)) - (160842737*pow4(Mst1))/(231525.*
        pow4(Msq)) + (29440*Dmglst1*s2t*pow4(Mst1))/(9.*Mt*pow4(Msq)) + (1520*
        z2*pow4(Mst1))/(3.*pow4(Msq)) - (2240*Dmglst1*s2t*z2*pow4(Mst1))/(3.*
        Mt*pow4(Msq)) - (4550*Cbeta*MuSUSY*s2t*pow2(Dmglst1)*pow4(Mst1))/(3.*
        Mt*Sbeta*pow2(Mst2)*pow4(Msq)) + (522392*pow2(Mst2)*pow2(s2t)*pow4(
        Mst1))/(15435.*pow2(Mt)*pow4(Msq)) - (2735*pow2(Dmglst1)*pow2(MuSUSY)*
        pow2(s2t)*pow4(Mst1))/(3.*pow2(Mst2)*pow2(Mt)*pow4(Msq)) - (1000*z2*
        pow2(Dmglst1)*pow2(MuSUSY)*pow2(s2t)*pow4(Mst1))/(pow2(Mst2)*pow2(Mt)*
        pow4(Msq)) + (2735*pow2(Dmglst1)*pow2(MuSUSY)*pow2(s2t)*pow4(Mst1))/(3.
        *pow2(Mst2)*pow2(Mt)*pow2(Sbeta)*pow4(Msq)) + (1000*z2*pow2(Dmglst1)*
        pow2(MuSUSY)*pow2(s2t)*pow4(Mst1))/(pow2(Mst2)*pow2(Mt)*pow2(Sbeta)*
        pow4(Msq)) - (9395*Cbeta*MuSUSY*pow2(Dmglst1)*pow3(s2t)*pow4(Mst1))/(9.
        *Sbeta*pow3(Mt)*pow4(Msq)) - (600*Cbeta*MuSUSY*z2*pow2(Dmglst1)*pow3(
        s2t)*pow4(Mst1))/(Sbeta*pow3(Mt)*pow4(Msq)) - (9035*Dmglst1*pow2(Mst2)*
        pow3(s2t)*pow4(Mst1))/(9.*pow3(Mt)*pow4(Msq)) - (200*Dmglst1*z2*pow2(
        Mst2)*pow3(s2t)*pow4(Mst1))/(pow3(Mt)*pow4(Msq)) - (90465364*pow2(
        Dmglst1)*pow2(Mst1))/(2025.*pow4(Mst2)) - (33792*z2*pow2(Dmglst1)*pow2(
        Mst1))/pow4(Mst2) + (64832*z3*pow2(Dmglst1)*pow2(Mst1))/(3.*pow4(Mst2))
        - (799880*Dmglst1*pow3(Mst1))/(27.*pow4(Mst2)) - (19904*Dmglst1*z2*
        pow3(Mst1))/pow4(Mst2) + (42112*Dmglst1*z3*pow3(Mst1))/(3.*pow4(Mst2))
        + (5583152*s2t*pow2(Dmglst1)*pow3(Mst1))/(81.*Mt*pow4(Mst2)) + (
        1368160*s2t*z2*pow2(Dmglst1)*pow3(Mst1))/(9.*Mt*pow4(Mst2)) - (270208*
        s2t*z3*pow2(Dmglst1)*pow3(Mst1))/(3.*Mt*pow4(Mst2)) - (25346116504*
        pow4(Mst1))/(3.472875e6*pow4(Mst2)) + (864592*Dmglst1*s2t*pow4(Mst1))/(
        81.*Mt*pow4(Mst2)) - (38128*z2*pow4(Mst1))/(9.*pow4(Mst2)) + (235936*
        Dmglst1*s2t*z2*pow4(Mst1))/(3.*Mt*pow4(Mst2)) + (10144*z3*pow4(Mst1))/(
        3.*pow4(Mst2)) - (46848*Dmglst1*s2t*z3*pow4(Mst1))/(Mt*pow4(Mst2)) - (
        17120*Cbeta*MuSUSY*s2t*pow2(Dmglst1)*pow4(Mst1))/(3.*Mt*Sbeta*pow2(Msq)
        *pow4(Mst2)) + (194*pow2(Dmglst1)*pow2(s2t)*pow4(Mst1))/(3.*pow2(Mt)*
        pow4(Mst2)) - (16*B4*pow2(Dmglst1)*pow2(s2t)*pow4(Mst1))/(pow2(Mt)*
        pow4(Mst2)) - (8*DN*pow2(Dmglst1)*pow2(s2t)*pow4(Mst1))/(pow2(Mt)*pow4(
        Mst2)) - (1400*z2*pow2(Dmglst1)*pow2(s2t)*pow4(Mst1))/(3.*pow2(Mt)*
        pow4(Mst2)) + (1048*z3*pow2(Dmglst1)*pow2(s2t)*pow4(Mst1))/(pow2(Mt)*
        pow4(Mst2)) + (120*z4*pow2(Dmglst1)*pow2(s2t)*pow4(Mst1))/(pow2(Mt)*
        pow4(Mst2)) - (2140*pow2(Dmglst1)*pow2(MuSUSY)*pow2(s2t)*pow4(Mst1))/(
        3.*pow2(Msq)*pow2(Mt)*pow4(Mst2)) - (1440*z2*pow2(Dmglst1)*pow2(MuSUSY)
        *pow2(s2t)*pow4(Mst1))/(pow2(Msq)*pow2(Mt)*pow4(Mst2)) + (2140*pow2(
        Dmglst1)*pow2(MuSUSY)*pow2(s2t)*pow4(Mst1))/(3.*pow2(Msq)*pow2(Mt)*
        pow2(Sbeta)*pow4(Mst2)) + (1440*z2*pow2(Dmglst1)*pow2(MuSUSY)*pow2(s2t)
        *pow4(Mst1))/(pow2(Msq)*pow2(Mt)*pow2(Sbeta)*pow4(Mst2)) - (388*pow2(
        Dmglst1)*pow2(s2t)*pow2(Sbeta)*pow4(Mst1))/(3.*pow2(Mt)*pow4(Mst2)) + (
        32*B4*pow2(Dmglst1)*pow2(s2t)*pow2(Sbeta)*pow4(Mst1))/(pow2(Mt)*pow4(
        Mst2)) + (16*DN*pow2(Dmglst1)*pow2(s2t)*pow2(Sbeta)*pow4(Mst1))/(pow2(
        Mt)*pow4(Mst2)) + (2800*z2*pow2(Dmglst1)*pow2(s2t)*pow2(Sbeta)*pow4(
        Mst1))/(3.*pow2(Mt)*pow4(Mst2)) - (2096*z3*pow2(Dmglst1)*pow2(s2t)*
        pow2(Sbeta)*pow4(Mst1))/(pow2(Mt)*pow4(Mst2)) - (240*z4*pow2(Dmglst1)*
        pow2(s2t)*pow2(Sbeta)*pow4(Mst1))/(pow2(Mt)*pow4(Mst2)) + (291380477*
        Cbeta*MuSUSY*pow2(Dmglst1)*pow3(s2t)*pow4(Mst1))/(16200.*Sbeta*pow3(Mt)
        *pow4(Mst2)) - (256*B4*Cbeta*MuSUSY*pow2(Dmglst1)*pow3(s2t)*pow4(Mst1))
        /(3.*Sbeta*pow3(Mt)*pow4(Mst2)) - (128*Cbeta*DN*MuSUSY*pow2(Dmglst1)*
        pow3(s2t)*pow4(Mst1))/(3.*Sbeta*pow3(Mt)*pow4(Mst2)) + (119038*Cbeta*
        MuSUSY*z2*pow2(Dmglst1)*pow3(s2t)*pow4(Mst1))/(9.*Sbeta*pow3(Mt)*pow4(
        Mst2)) - (31808*Cbeta*MuSUSY*z3*pow2(Dmglst1)*pow3(s2t)*pow4(Mst1))/(3.
        *Sbeta*pow3(Mt)*pow4(Mst2)) + (640*Cbeta*MuSUSY*z4*pow2(Dmglst1)*pow3(
        s2t)*pow4(Mst1))/(Sbeta*pow3(Mt)*pow4(Mst2)) + (40864*pow2(s2t)*pow4(
        Mst2))/(675.*pow2(Msq)*pow2(Mt)) - (8*pow2(s2t)*pow4(Mst2))/(pow2(Mst1)
        *pow2(Mt)) + (16*z2*pow2(s2t)*pow4(Mst2))/(pow2(Mst1)*pow2(Mt)) + (80*
        Dmglst1*pow3(s2t)*pow4(Mst2))/(3.*pow2(Msq)*pow3(Mt)) + (80*Mst1*pow3(
        s2t)*pow4(Mst2))/(3.*pow2(Msq)*pow3(Mt)) - (13597579*pow4(Mst2))/(
        77175.*pow4(Msq)) - (6760*Dmglst1*s2t*pow4(Mst2))/(27.*Mt*pow4(Msq)) -
        (6760*Mst1*s2t*pow4(Mst2))/(27.*Mt*pow4(Msq)) + (80*z2*pow4(Mst2))/
        pow4(Msq) + (320*Dmglst1*s2t*z2*pow4(Mst2))/(3.*Mt*pow4(Msq)) + (320*
        Mst1*s2t*z2*pow4(Mst2))/(3.*Mt*pow4(Msq)) - (130*Dmglst1*Mst1*pow2(s2t)
        *pow4(Mst2))/(3.*pow2(Mt)*pow4(Msq)) - (65*pow2(Dmglst1)*pow2(s2t)*
        pow4(Mst2))/(3.*pow2(Mt)*pow4(Msq)) - (438008*pow2(Mst1)*pow2(s2t)*
        pow4(Mst2))/(15435.*pow2(Mt)*pow4(Msq)) + (40*Mst1*pow2(Dmglst1)*pow3(
        s2t)*pow4(Mst2))/(pow3(Mt)*pow4(Msq)) + (455*Dmglst1*pow2(Mst1)*pow3(
        s2t)*pow4(Mst2))/(9.*pow3(Mt)*pow4(Msq)) + (215*pow3(Mst1)*pow3(s2t)*
        pow4(Mst2))/(9.*pow3(Mt)*pow4(Msq)) + (438203*pow2(Dmglst1)*pow2(Mst1)*
        pow4(s2t))/(864.*pow4(Mt)) - (76*B4*pow2(Dmglst1)*pow2(Mst1)*pow4(s2t))
        /(3.*pow4(Mt)) - (14*DN*pow2(Dmglst1)*pow2(Mst1)*pow4(s2t))/pow4(Mt) -
        (1655*z2*pow2(Dmglst1)*pow2(Mst1)*pow4(s2t))/(2.*pow4(Mt)) + (3975*z3*
        pow2(Dmglst1)*pow2(Mst1)*pow4(s2t))/(2.*pow4(Mt)) + (1010*z4*pow2(
        Dmglst1)*pow2(Mst1)*pow4(s2t))/(3.*pow4(Mt)) - (10*pow2(Msq)*pow2(Mst1)
        *pow4(s2t))/pow4(Mt) + (20*z2*pow2(Msq)*pow2(Mst1)*pow4(s2t))/pow4(Mt)
        + (1151*Dmglst1*Mst1*pow2(Mst2)*pow4(s2t))/(9.*pow4(Mt)) + (8*B4*
        Dmglst1*Mst1*pow2(Mst2)*pow4(s2t))/pow4(Mt) + (20*Dmglst1*DN*Mst1*pow2(
        Mst2)*pow4(s2t))/(3.*pow4(Mt)) + (388*Dmglst1*Mst1*z2*pow2(Mst2)*pow4(
        s2t))/pow4(Mt) + (270*Dmglst1*Mst1*z3*pow2(Mst2)*pow4(s2t))/pow4(Mt) -
        (1060*Dmglst1*Mst1*z4*pow2(Mst2)*pow4(s2t))/(3.*pow4(Mt)) - (5885*pow2(
        Dmglst1)*pow2(Mst2)*pow4(s2t))/(18.*pow4(Mt)) + (4*B4*pow2(Dmglst1)*
        pow2(Mst2)*pow4(s2t))/pow4(Mt) + (10*DN*pow2(Dmglst1)*pow2(Mst2)*pow4(
        s2t))/(3.*pow4(Mt)) + (194*z2*pow2(Dmglst1)*pow2(Mst2)*pow4(s2t))/pow4(
        Mt) + (437*z3*pow2(Dmglst1)*pow2(Mst2)*pow4(s2t))/pow4(Mt) - (530*z4*
        pow2(Dmglst1)*pow2(Mst2)*pow4(s2t))/(3.*pow4(Mt)) - (10*pow2(Msq)*pow2(
        Mst2)*pow4(s2t))/pow4(Mt) + (20*z2*pow2(Msq)*pow2(Mst2)*pow4(s2t))/
        pow4(Mt) - (5876588*pow2(Mst1)*pow2(Mst2)*pow4(s2t))/(10125.*pow4(Mt))
        + (4*B4*pow2(Mst1)*pow2(Mst2)*pow4(s2t))/(3.*pow4(Mt)) + (2*DN*pow2(
        Mst1)*pow2(Mst2)*pow4(s2t))/pow4(Mt) + (400*OepS2*pow2(Mst1)*pow2(Mst2)
        *pow4(s2t))/(81.*pow4(Mt)) - (394*S2*pow2(Mst1)*pow2(Mst2)*pow4(s2t))/
        pow4(Mt) - (200*T1ep*pow2(Mst1)*pow2(Mst2)*pow4(s2t))/(27.*pow4(Mt)) +
        (5737*z2*pow2(Mst1)*pow2(Mst2)*pow4(s2t))/(27.*pow4(Mt)) + (11662*z3*
        pow2(Mst1)*pow2(Mst2)*pow4(s2t))/(27.*pow4(Mt)) - (4330*z4*pow2(Mst1)*
        pow2(Mst2)*pow4(s2t))/(27.*pow4(Mt)) - (3850*pow2(Dmglst1)*pow2(Mst1)*
        pow2(Mst2)*pow4(s2t))/(9.*pow2(Msq)*pow4(Mt)) - (200*z2*pow2(Dmglst1)*
        pow2(Mst1)*pow2(Mst2)*pow4(s2t))/(pow2(Msq)*pow4(Mt)) + (8*pow2(
        Dmglst1)*pow2(Mst1)*pow2(z2)*pow4(s2t))/(3.*pow4(Mt)) - (16*Dmglst1*
        Mst1*pow2(Mst2)*pow2(z2)*pow4(s2t))/(3.*pow4(Mt)) - (8*pow2(Dmglst1)*
        pow2(Mst2)*pow2(z2)*pow4(s2t))/(3.*pow4(Mt)) - (74*pow2(Mst1)*pow2(
        Mst2)*pow2(z2)*pow4(s2t))/(9.*pow4(Mt)) - (773389*Dmglst1*pow3(Mst1)*
        pow4(s2t))/(1296.*pow4(Mt)) - (8*B4*Dmglst1*pow3(Mst1)*pow4(s2t))/pow4(
        Mt) - (20*Dmglst1*DN*pow3(Mst1)*pow4(s2t))/(3.*pow4(Mt)) - (7153*
        Dmglst1*z2*pow3(Mst1)*pow4(s2t))/(9.*pow4(Mt)) + (5077*Dmglst1*z3*pow3(
        Mst1)*pow4(s2t))/(3.*pow4(Mt)) + (1060*Dmglst1*z4*pow3(Mst1)*pow4(s2t))
        /(3.*pow4(Mt)) - (2740*Dmglst1*pow2(Mst2)*pow3(Mst1)*pow4(s2t))/(9.*
        pow2(Msq)*pow4(Mt)) - (400*Dmglst1*z2*pow2(Mst2)*pow3(Mst1)*pow4(s2t))/
        (3.*pow2(Msq)*pow4(Mt)) + (16*Dmglst1*pow2(z2)*pow3(Mst1)*pow4(s2t))/(
        3.*pow4(Mt)) - (257823187891*pow4(Mst1)*pow4(s2t))/(6.66792e8*pow4(Mt))
        - (40*B4*pow4(Mst1)*pow4(s2t))/(3.*pow4(Mt)) - (2*DN*pow4(Mst1)*pow4(
        s2t))/(3.*pow4(Mt)) - (88*OepS2*pow4(Mst1)*pow4(s2t))/(243.*pow4(Mt)) +
        (281*S2*pow4(Mst1)*pow4(s2t))/(9.*pow4(Mt)) + (44*T1ep*pow4(Mst1)*pow4(
        s2t))/(81.*pow4(Mt)) - (98497*z2*pow4(Mst1)*pow4(s2t))/(324.*pow4(Mt))
        + (38672*z3*pow4(Mst1)*pow4(s2t))/(81.*pow4(Mt)) + (6124*z4*pow4(Mst1)*
        pow4(s2t))/(81.*pow4(Mt)) + (10*pow2(Msq)*pow4(Mst1)*pow4(s2t))/(pow2(
        Mst2)*pow4(Mt)) - (20*z2*pow2(Msq)*pow4(Mst1)*pow4(s2t))/(pow2(Mst2)*
        pow4(Mt)) - (4136*pow2(Mst2)*pow4(Mst1)*pow4(s2t))/(45.*pow2(Msq)*pow4(
        Mt)) - (100*z2*pow2(Mst2)*pow4(Mst1)*pow4(s2t))/(3.*pow2(Msq)*pow4(Mt))
        + (371*pow2(z2)*pow4(Mst1)*pow4(s2t))/(27.*pow4(Mt)) + (53749*pow4(
        Mst2)*pow4(s2t))/(216.*pow4(Mt)) + (12*B4*pow4(Mst2)*pow4(s2t))/pow4(
        Mt) - (4*DN*pow4(Mst2)*pow4(s2t))/(3.*pow4(Mt)) + (7*Dmglst1*pow4(Mst2)
        *pow4(s2t))/(6.*Mst1*pow4(Mt)) + (8*OepS2*pow4(Mst2)*pow4(s2t))/(9.*
        pow4(Mt)) + (99*S2*pow4(Mst2)*pow4(s2t))/pow4(Mt) - (4*T1ep*pow4(Mst2)*
        pow4(s2t))/(3.*pow4(Mt)) - (109*z2*pow4(Mst2)*pow4(s2t))/pow4(Mt) + (
        66*Dmglst1*z2*pow4(Mst2)*pow4(s2t))/(Mst1*pow4(Mt)) + (328*z3*pow4(
        Mst2)*pow4(s2t))/(3.*pow4(Mt)) - (65*Dmglst1*z3*pow4(Mst2)*pow4(s2t))/(
        Mst1*pow4(Mt)) + (242*z4*pow4(Mst2)*pow4(s2t))/(3.*pow4(Mt)) + (820*
        Dmglst1*Mst1*pow4(Mst2)*pow4(s2t))/(9.*pow2(Msq)*pow4(Mt)) + (20*
        Dmglst1*Mst1*z2*pow4(Mst2)*pow4(s2t))/(pow2(Msq)*pow4(Mt)) + (800*pow2(
        Dmglst1)*pow4(Mst2)*pow4(s2t))/(9.*pow2(Msq)*pow4(Mt)) + (10*z2*pow2(
        Dmglst1)*pow4(Mst2)*pow4(s2t))/(pow2(Msq)*pow4(Mt)) + (6415*pow2(
        Dmglst1)*pow4(Mst2)*pow4(s2t))/(36.*pow2(Mst1)*pow4(Mt)) + (33*z2*pow2(
        Dmglst1)*pow4(Mst2)*pow4(s2t))/(pow2(Mst1)*pow4(Mt)) - (381*z3*pow2(
        Dmglst1)*pow4(Mst2)*pow4(s2t))/(2.*pow2(Mst1)*pow4(Mt)) + (10*pow2(Msq)
        *pow4(Mst2)*pow4(s2t))/(pow2(Mst1)*pow4(Mt)) - (20*z2*pow2(Msq)*pow4(
        Mst2)*pow4(s2t))/(pow2(Mst1)*pow4(Mt)) + (2714*pow2(Mst1)*pow4(Mst2)*
        pow4(s2t))/(45.*pow2(Msq)*pow4(Mt)) + (10*z2*pow2(Mst1)*pow4(Mst2)*
        pow4(s2t))/(pow2(Msq)*pow4(Mt)) - (35*pow2(z2)*pow4(Mst2)*pow4(s2t))/(
        3.*pow4(Mt)) + (10765*pow2(Dmglst1)*pow2(Mst1)*pow4(Mst2)*pow4(s2t))/(
        72.*pow4(Msq)*pow4(Mt)) + (25*z2*pow2(Dmglst1)*pow2(Mst1)*pow4(Mst2)*
        pow4(s2t))/(pow4(Msq)*pow4(Mt)) + (2605*Dmglst1*pow3(Mst1)*pow4(Mst2)*
        pow4(s2t))/(36.*pow4(Msq)*pow4(Mt)) + (50*Dmglst1*z2*pow3(Mst1)*pow4(
        Mst2)*pow4(s2t))/(3.*pow4(Msq)*pow4(Mt)) + (1009961*pow4(Mst1)*pow4(
        Mst2)*pow4(s2t))/(51450.*pow4(Msq)*pow4(Mt)) + (25*z2*pow4(Mst1)*pow4(
        Mst2)*pow4(s2t))/(6.*pow4(Msq)*pow4(Mt)) + (194*pow2(Dmglst1)*pow2(s2t)
        *pow4(Mst1)*pow4(Sbeta))/(3.*pow2(Mt)*pow4(Mst2)) - (16*B4*pow2(
        Dmglst1)*pow2(s2t)*pow4(Mst1)*pow4(Sbeta))/(pow2(Mt)*pow4(Mst2)) - (8*
        DN*pow2(Dmglst1)*pow2(s2t)*pow4(Mst1)*pow4(Sbeta))/(pow2(Mt)*pow4(Mst2)
        ) - (1400*z2*pow2(Dmglst1)*pow2(s2t)*pow4(Mst1)*pow4(Sbeta))/(3.*pow2(
        Mt)*pow4(Mst2)) + (1048*z3*pow2(Dmglst1)*pow2(s2t)*pow4(Mst1)*pow4(
        Sbeta))/(pow2(Mt)*pow4(Mst2)) + (120*z4*pow2(Dmglst1)*pow2(s2t)*pow4(
        Mst1)*pow4(Sbeta))/(pow2(Mt)*pow4(Mst2)) - (320*s2t*pow5(Mst1))/(3.*Mt*
        pow2(Msq)*pow2(Mst2)) - (1280*s2t*z2*pow5(Mst1))/(Mt*pow2(Msq)*pow2(
        Mst2)) + (6640*pow3(s2t)*pow5(Mst1))/(9.*pow2(Msq)*pow3(Mt)) + (640*z2*
        pow3(s2t)*pow5(Mst1))/(3.*pow2(Msq)*pow3(Mt)) - (622151*pow3(s2t)*pow5(
        Mst1))/(324.*pow2(Mst2)*pow3(Mt)) + (2948*z2*pow3(s2t)*pow5(Mst1))/(9.*
        pow2(Mst2)*pow3(Mt)) + (1384*z3*pow3(s2t)*pow5(Mst1))/(pow2(Mst2)*pow3(
        Mt)) + (6140*Cbeta*Dmglst1*MuSUSY*pow3(s2t)*pow5(Mst1))/(9.*Sbeta*pow2(
        Msq)*pow2(Mst2)*pow3(Mt)) + (8000*s2t*pow5(Mst1))/(9.*Mt*pow4(Msq)) - (
        320*s2t*z2*pow5(Mst1))/(Mt*pow4(Msq)) - (1820*Cbeta*Dmglst1*MuSUSY*s2t*
        pow5(Mst1))/(3.*Mt*Sbeta*pow2(Mst2)*pow4(Msq)) + (3560*Cbeta*MuSUSY*
        pow2(Dmglst1)*pow2(s2t)*pow5(Mst1))/(Sbeta*pow2(Mst2)*pow2(Mt)*pow4(
        Msq)) + (4320*Cbeta*MuSUSY*z2*pow2(Dmglst1)*pow2(s2t)*pow5(Mst1))/(
        Sbeta*pow2(Mst2)*pow2(Mt)*pow4(Msq)) - (1030*Dmglst1*pow2(MuSUSY)*pow2(
        s2t)*pow5(Mst1))/(3.*pow2(Mst2)*pow2(Mt)*pow4(Msq)) - (880*Dmglst1*z2*
        pow2(MuSUSY)*pow2(s2t)*pow5(Mst1))/(3.*pow2(Mst2)*pow2(Mt)*pow4(Msq)) +
        (1030*Dmglst1*pow2(MuSUSY)*pow2(s2t)*pow5(Mst1))/(3.*pow2(Mst2)*pow2(
        Mt)*pow2(Sbeta)*pow4(Msq)) + (880*Dmglst1*z2*pow2(MuSUSY)*pow2(s2t)*
        pow5(Mst1))/(3.*pow2(Mst2)*pow2(Mt)*pow2(Sbeta)*pow4(Msq)) - (3970*
        Cbeta*Dmglst1*MuSUSY*pow3(s2t)*pow5(Mst1))/(9.*Sbeta*pow3(Mt)*pow4(Msq)
        ) - (640*Cbeta*Dmglst1*MuSUSY*z2*pow3(s2t)*pow5(Mst1))/(3.*Sbeta*pow3(
        Mt)*pow4(Msq)) - (1715*pow2(Mst2)*pow3(s2t)*pow5(Mst1))/(9.*pow3(Mt)*
        pow4(Msq)) - (40*z2*pow2(Mst2)*pow3(s2t)*pow5(Mst1))/(pow3(Mt)*pow4(
        Msq)) - (324688*s2t*pow5(Mst1))/(81.*Mt*pow4(Mst2)) + (146848*s2t*z2*
        pow5(Mst1))/(9.*Mt*pow4(Mst2)) - (28672*s2t*z3*pow5(Mst1))/(3.*Mt*pow4(
        Mst2)) - (4160*Cbeta*Dmglst1*MuSUSY*s2t*pow5(Mst1))/(3.*Mt*Sbeta*pow2(
        Msq)*pow4(Mst2)) + (11116*Dmglst1*pow2(s2t)*pow5(Mst1))/(81.*pow2(Mt)*
        pow4(Mst2)) - (32*B4*Dmglst1*pow2(s2t)*pow5(Mst1))/(pow2(Mt)*pow4(Mst2)
        ) - (16*Dmglst1*DN*pow2(s2t)*pow5(Mst1))/(pow2(Mt)*pow4(Mst2)) - (3536*
        Dmglst1*z2*pow2(s2t)*pow5(Mst1))/(9.*pow2(Mt)*pow4(Mst2)) + (416*
        Dmglst1*z3*pow2(s2t)*pow5(Mst1))/(pow2(Mt)*pow4(Mst2)) + (240*Dmglst1*
        z4*pow2(s2t)*pow5(Mst1))/(pow2(Mt)*pow4(Mst2)) + (400*Cbeta*MuSUSY*
        pow2(Dmglst1)*pow2(s2t)*pow5(Mst1))/(Sbeta*pow2(Msq)*pow2(Mt)*pow4(
        Mst2)) + (7680*Cbeta*MuSUSY*z2*pow2(Dmglst1)*pow2(s2t)*pow5(Mst1))/(
        Sbeta*pow2(Msq)*pow2(Mt)*pow4(Mst2)) + (1160*Dmglst1*pow2(MuSUSY)*pow2(
        s2t)*pow5(Mst1))/(3.*pow2(Msq)*pow2(Mt)*pow4(Mst2)) - (2240*Dmglst1*z2*
        pow2(MuSUSY)*pow2(s2t)*pow5(Mst1))/(3.*pow2(Msq)*pow2(Mt)*pow4(Mst2)) -
        (1160*Dmglst1*pow2(MuSUSY)*pow2(s2t)*pow5(Mst1))/(3.*pow2(Msq)*pow2(Mt)
        *pow2(Sbeta)*pow4(Mst2)) + (2240*Dmglst1*z2*pow2(MuSUSY)*pow2(s2t)*
        pow5(Mst1))/(3.*pow2(Msq)*pow2(Mt)*pow2(Sbeta)*pow4(Mst2)) - (22232*
        Dmglst1*pow2(s2t)*pow2(Sbeta)*pow5(Mst1))/(81.*pow2(Mt)*pow4(Mst2)) + (
        64*B4*Dmglst1*pow2(s2t)*pow2(Sbeta)*pow5(Mst1))/(pow2(Mt)*pow4(Mst2)) +
        (32*Dmglst1*DN*pow2(s2t)*pow2(Sbeta)*pow5(Mst1))/(pow2(Mt)*pow4(Mst2))
        + (7072*Dmglst1*z2*pow2(s2t)*pow2(Sbeta)*pow5(Mst1))/(9.*pow2(Mt)*pow4(
        Mst2)) - (832*Dmglst1*z3*pow2(s2t)*pow2(Sbeta)*pow5(Mst1))/(pow2(Mt)*
        pow4(Mst2)) - (480*Dmglst1*z4*pow2(s2t)*pow2(Sbeta)*pow5(Mst1))/(pow2(
        Mt)*pow4(Mst2)) + (302282291*Cbeta*Dmglst1*MuSUSY*pow3(s2t)*pow5(Mst1))
        /(40500.*Sbeta*pow3(Mt)*pow4(Mst2)) + (63364*Cbeta*Dmglst1*MuSUSY*z2*
        pow3(s2t)*pow5(Mst1))/(45.*Sbeta*pow3(Mt)*pow4(Mst2)) - (1408*Cbeta*
        Dmglst1*MuSUSY*z3*pow3(s2t)*pow5(Mst1))/(3.*Sbeta*pow3(Mt)*pow4(Mst2))
        + (6960*Cbeta*MuSUSY*pow2(Dmglst1)*pow5(Mst1))/(Sbeta*pow4(Msq)*pow4(
        Mst2)) - (8960*Cbeta*MuSUSY*z2*pow2(Dmglst1)*pow5(Mst1))/(Sbeta*pow4(
        Msq)*pow4(Mst2)) - (35920*s2t*pow2(Dmglst1)*pow2(MuSUSY)*pow5(Mst1))/(
        9.*Mt*pow4(Msq)*pow4(Mst2)) + (4160*s2t*z2*pow2(Dmglst1)*pow2(MuSUSY)*
        pow5(Mst1))/(Mt*pow4(Msq)*pow4(Mst2)) + (35920*s2t*pow2(Dmglst1)*pow2(
        MuSUSY)*pow5(Mst1))/(9.*Mt*pow2(Sbeta)*pow4(Msq)*pow4(Mst2)) - (4160*
        s2t*z2*pow2(Dmglst1)*pow2(MuSUSY)*pow5(Mst1))/(Mt*pow2(Sbeta)*pow4(Msq)
        *pow4(Mst2)) + (11116*Dmglst1*pow2(s2t)*pow4(Sbeta)*pow5(Mst1))/(81.*
        pow2(Mt)*pow4(Mst2)) - (32*B4*Dmglst1*pow2(s2t)*pow4(Sbeta)*pow5(Mst1))
        /(pow2(Mt)*pow4(Mst2)) - (16*Dmglst1*DN*pow2(s2t)*pow4(Sbeta)*pow5(
        Mst1))/(pow2(Mt)*pow4(Mst2)) - (3536*Dmglst1*z2*pow2(s2t)*pow4(Sbeta)*
        pow5(Mst1))/(9.*pow2(Mt)*pow4(Mst2)) + (416*Dmglst1*z3*pow2(s2t)*pow4(
        Sbeta)*pow5(Mst1))/(pow2(Mt)*pow4(Mst2)) + (240*Dmglst1*z4*pow2(s2t)*
        pow4(Sbeta)*pow5(Mst1))/(pow2(Mt)*pow4(Mst2)) - (7399303*Cbeta*MuSUSY*
        pow3(s2t)*pow6(Mst2))/(926100.*Sbeta*pow3(Mt)*pow4(Msq)) - (pow2(log(
        pow2(Mst1)/pow2(Msq)))*(105*Sbeta*log(pow2(Mst1)/pow2(Mst2))*pow2(Mst1)
        *pow2(Mst2)*pow4(Msq)*(8*(-pow2(Mst1) + pow2(Mst2))*pow2(Mt)*pow2(s2t)
        + 64*(Dmglst1 + Mst1)*s2t*pow3(Mt) - 16*(Dmglst1 + Mst1)*Mt*(pow2(Mst1)
        + pow2(Mst2))*pow3(s2t) + 16*pow4(Mt) + (pow2(Mst1) - pow2(Mst2))*(12*
        Dmglst1*Mst1 + 6*pow2(Dmglst1) + 7*pow2(Mst1) + pow2(Mst2))*pow4(s2t))
        + Sbeta*(105*pow2(Dmglst1)*(16*pow2(Mst1)*pow2(Mst2)*(-11*Mt*pow2(Mst1)
        - Mt*pow2(Mst2) - 4*Mst1*s2t*pow2(Mst2) + 4*s2t*pow3(Mst1))*pow3(Mt) -
        64*pow2(Msq)*pow2(Mst1)*pow2(Mst2)*pow4(Mt) + 3*pow4(Msq)*(-pow2(-4*
        Mst2*pow2(Mt) + pow2(s2t)*pow3(Mst2)) + pow4(Mst1)*(8*pow2(Mt)*pow2(
        s2t) + pow2(Mst2)*pow4(s2t)) + pow2(Mst1)*(-16*pow2(Mst2)*pow2(Mt)*
        pow2(s2t) - 16*pow4(Mt) + pow4(Mst2)*pow4(s2t)))) + 70*Dmglst1*Mst1*(
        64*Mst1*pow2(Msq)*pow2(Mst2)*(-3*Mst1*Mt + s2t*pow2(Mst1) - s2t*pow2(
        Mst2))*pow3(Mt) + 16*Mst1*pow2(Mst2)*pow3(Mt)*(-3*Mst1*Mt*pow2(Mst2) -
        6*s2t*pow2(Mst1)*pow2(Mst2) - 13*Mt*pow3(Mst1) + 7*s2t*pow4(Mst1) -
        s2t*pow4(Mst2)) + 3*pow4(Msq)*(-3*pow2(-4*Mst2*pow2(Mt) + pow2(s2t)*
        pow3(Mst2)) + 16*Mt*pow2(Mst2)*pow3(Mst1)*pow3(s2t) - 16*Mst1*Mt*pow3(
        s2t)*pow4(Mst2) + 3*pow4(Mst1)*(8*pow2(Mt)*pow2(s2t) + pow2(Mst2)*pow4(
        s2t)) + pow2(Mst1)*(-48*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 48*pow4(Mt) +
        3*pow4(Mst2)*pow4(s2t)))) + pow2(Mst1)*(-224*pow2(Msq)*pow2(Mst2)*pow2(
        Mt)*(20*Mst1*Mt*s2t*pow2(Mst2) + 7*pow2(Mst2)*pow2(Mt) + pow2(Mst1)*(
        37*pow2(Mt) - 2*pow2(Mst2)*pow2(s2t)) - 20*Mt*s2t*pow3(Mst1) + pow2(
        s2t)*pow4(Mst1) + pow2(s2t)*pow4(Mst2)) + 2*pow2(Mst2)*(-1120*s2t*pow2(
        Mst2)*pow3(Mst1)*pow3(Mt) - 560*Mst1*s2t*pow3(Mt)*pow4(Mst2) - 228*
        pow4(Mst2)*pow4(Mt) + 30*pow2(Mst1)*(pow2(Mt)*pow2(s2t)*pow4(Mst2) -
        28*pow2(Mst2)*pow4(Mt)) + pow4(Mst1)*(30*pow2(Mst2)*pow2(Mt)*pow2(s2t)
        - 2468*pow4(Mt) + 3*pow4(Mst2)*pow4(s2t)) + 1680*s2t*pow3(Mt)*pow5(
        Mst1)) + 35*pow4(Msq)*(96*Mt*pow2(Mst2)*pow3(Mst1)*pow3(s2t) + 72*pow2(
        Mt)*pow2(s2t)*pow4(Mst2) - 96*Mst1*Mt*pow3(s2t)*pow4(Mst2) - 176*pow2(
        Mst2)*pow4(Mt) + 3*pow4(Mst1)*(24*pow2(Mt)*pow2(s2t) + pow2(Mst2)*pow4(
        s2t)) - 3*pow2(Mst1)*(48*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 48*pow4(Mt) -
        7*pow4(Mst2)*pow4(s2t)) - 15*pow4(s2t)*pow6(Mst2)))) + 4*Cbeta*Mt*
        MuSUSY*pow2(Mst1)*pow3(s2t)*pow8(Mst2)))/(21.*Sbeta*pow2(Mst1)*pow2(
        Mst2)*pow4(Msq)*pow4(Mt)) + (pow3(log(pow2(Mst1)/pow2(Mst2)))*(128*
        Cbeta*Dmglst1*(849*Dmglst1 + 431*Mst1)*Mt*MuSUSY*pow3(s2t)*pow4(Mst1) +
        Sbeta*(-96*pow3(Mst1)*(-186*s2t*pow2(Mst2)*pow3(Mt) + 13*Mt*pow3(s2t)*
        pow4(Mst2)) - 544*pow4(Mst2)*pow4(Mt) + pow4(Mst1)*(-7136*pow2(Mst2)*
        pow2(Mt)*pow2(s2t) - 4640*pow4(Mt) + 1535*pow4(Mst2)*pow4(s2t)) + 16*(
        1536*s2t*pow3(Mt) + Mt*pow2(Mst2)*pow3(s2t))*pow5(Mst1) - 3152*pow2(Mt)
        *pow2(s2t)*pow6(Mst2) + 32*Mst1*(-28*s2t*pow3(Mt)*pow4(Mst2) + 73*Mt*
        pow3(s2t)*pow6(Mst2)) + 2*pow2(Dmglst1)*(48*pow3(Mst1)*(2262*s2t*pow3(
        Mt) - 61*Mt*pow2(Mst2)*pow3(s2t)) + 456*pow2(Mt)*pow2(s2t)*pow2(-1 +
        pow2(Sbeta))*pow4(Mst1) - 320*pow2(Mt)*pow2(s2t)*pow4(Mst2) + Mst1*(
        27168*s2t*pow2(Mst2)*pow3(Mt) - 7608*Mt*pow3(s2t)*pow4(Mst2)) + 1056*
        pow2(Mst2)*pow4(Mt) + pow2(Mst1)*(-24960*pow2(Mst2)*pow2(Mt)*pow2(s2t)
        - 17792*pow4(Mt) + 975*pow4(Mst2)*pow4(s2t)) + 449*pow4(s2t)*pow6(Mst2)
        ) + 2*pow2(Mst1)*(-960*pow2(Mt)*pow2(s2t)*pow4(Mst2) + 1120*pow2(Mst2)*
        pow4(Mt) + 1097*pow4(s2t)*pow6(Mst2)) + 4*Dmglst1*(4*(7320*s2t*pow3(Mt)
        - 323*Mt*pow2(Mst2)*pow3(s2t))*pow4(Mst1) + 8*Mt*s2t*(-28*pow2(Mt) +
        73*pow2(Mst2)*pow2(s2t))*pow4(Mst2) - 16*pow2(Mst1)*(-845*s2t*pow2(
        Mst2)*pow3(Mt) + 178*Mt*pow3(s2t)*pow4(Mst2)) + pow3(Mst1)*(-7264*pow2(
        Mst2)*pow2(Mt)*pow2(s2t) - 5472*pow4(Mt) + 311*pow4(Mst2)*pow4(s2t)) +
        104*pow2(Mt)*pow2(s2t)*pow2(-1 + pow2(Sbeta))*pow5(Mst1) + Mst1*(-320*
        pow2(Mt)*pow2(s2t)*pow4(Mst2) + 1056*pow2(Mst2)*pow4(Mt) + 449*pow4(
        s2t)*pow6(Mst2))) + 639*pow4(s2t)*pow8(Mst2))))/(18.*Sbeta*pow4(Mst2)*
        pow4(Mt)) - (pow2(log(pow2(Mst1)/pow2(Mst2)))*(140*Dmglst1*(30*pow2(
        Mst2)*pow2(s2t)*pow5(Mst1)*(8*Mt*MuSUSY*pow2(Mst1)*(8*Cbeta*s2t*Sbeta*
        pow2(Mst2) + 11*Mt*MuSUSY*(-1 + pow2(Sbeta))) + 60*Mst1*Mt*s2t*pow2(
        Sbeta)*pow4(Mst2) - 5*pow2(s2t)*pow2(Sbeta)*pow6(Mst2)) - Mst1*Sbeta*
        pow4(Msq)*(93180*Cbeta*Mt*MuSUSY*pow3(s2t)*pow6(Mst1) + Sbeta*(54*pow2(
        -4*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow4(Mst2) - 8*pow3(Mst1)*(-26144*
        s2t*pow2(Mst2)*pow3(Mt) + 3033*Mt*pow3(s2t)*pow4(Mst2)) + pow4(Mst1)*(-
        86656*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 92736*pow4(Mt) + 617*pow4(Mst2)*
        pow4(s2t)) + 4*(184728*s2t*pow3(Mt) - 5569*Mt*pow2(Mst2)*pow3(s2t))*
        pow5(Mst1) - 800*pow2(Mt)*pow2(s2t)*pow2(-1 + pow2(Sbeta))*pow6(Mst1) +
        24*Mst1*(-128*s2t*pow3(Mt)*pow4(Mst2) + 797*Mt*pow3(s2t)*pow6(Mst2)) +
        4*pow2(Mst1)*(8200*pow2(Mt)*pow2(s2t)*pow4(Mst2) + 5648*pow2(Mst2)*
        pow4(Mt) + 1089*pow4(s2t)*pow6(Mst2)))) + 60*s2t*pow2(Msq)*pow3(Mst1)*(
        -192*Mt*pow2(Mst2)*(-5*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow2(Sbeta)*
        pow3(Mst1) + 112*s2t*pow2(Mt)*pow2(MuSUSY)*(-1 + pow2(Sbeta))*pow4(
        Mst1) + 48*Mst1*Mt*pow2(s2t)*pow2(Sbeta)*pow6(Mst2) + 20*pow2(Mst1)*
        pow2(Sbeta)*pow3(s2t)*pow6(Mst2) - 3*pow2(Sbeta)*pow3(s2t)*pow8(Mst2)))
        + 70*pow2(Dmglst1)*(-54*pow2(Sbeta)*pow2(-4*pow2(Mt) + pow2(Mst2)*pow2(
        s2t))*pow4(Msq)*pow4(Mst2) + 288*Mt*s2t*pow2(Msq)*pow2(Mst2)*pow2(
        Sbeta)*pow3(Mst1)*(pow2(Msq)*(-1438*pow2(Mt) + 393*pow2(Mst2)*pow2(s2t)
        ) + 20*pow2(s2t)*pow4(Mst2)) - 4*Mt*pow2(s2t)*(-900*MuSUSY*pow2(Mst2)*(
        3*Cbeta*s2t*Sbeta*pow2(Mst2) + 5*Mt*MuSUSY*(-1 + pow2(Sbeta))) - 6480*
        Mt*pow2(Msq)*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + Sbeta*(159493*Cbeta*
        MuSUSY*s2t + 552*Mt*Sbeta*pow2(-1 + pow2(Sbeta)))*pow4(Msq))*pow6(Mst1)
        - 80*Mt*s2t*pow2(Sbeta)*pow5(Mst1)*((34876*pow2(Mt) - 2117*pow2(Mst2)*
        pow2(s2t))*pow4(Msq) + 72*pow2(Msq)*(-40*pow2(Mst2)*pow2(Mt) + 9*pow2(
        s2t)*pow4(Mst2)) - 90*pow2(s2t)*pow6(Mst2)) - 4*pow2(Msq)*pow2(Mst1)*
        pow2(Mst2)*pow2(Sbeta)*(pow2(Msq)*(6472*pow2(Mst2)*pow2(Mt)*pow2(s2t) +
        5648*pow4(Mt) + 1089*pow4(Mst2)*pow4(s2t)) + 45*pow4(s2t)*pow6(Mst2)) +
        2880*MuSUSY*pow2(Mt)*(56*Cbeta*Sbeta*pow2(Mt) - 3*Cbeta*Sbeta*(16*pow2(
        Msq) + 9*pow2(Mst2))*pow2(s2t) - 26*Mt*MuSUSY*s2t*(-1 + pow2(Sbeta)))*
        pow7(Mst1) + 3*pow2(Sbeta)*pow4(Mst1)*(pow4(Msq)*(132992*pow2(Mst2)*
        pow2(Mt)*pow2(s2t) + 96384*pow4(Mt) - 5989*pow4(Mst2)*pow4(s2t)) +
        1200*pow2(Msq)*pow4(s2t)*pow6(Mst2) - 150*pow4(s2t)*pow8(Mst2))) +
        Sbeta*pow2(Mst1)*(-(Sbeta*(840*pow2(Msq)*(-1920*s2t*pow2(Mst2)*pow3(Mt)
        *pow5(Mst1) + 320*Mt*pow3(s2t)*pow4(Mst2)*pow5(Mst1) - 8*(-6*pow2(Mst1)
        *pow2(Mt)*pow2(s2t) + 20*Mt*pow3(Mst1)*pow3(s2t) + 12*pow4(Mt) + 7*
        pow4(Mst1)*pow4(s2t))*pow6(Mst2) + 9*(8*pow2(Mt)*pow2(s2t) - pow2(Mst1)
        *pow4(s2t))*pow8(Mst2)) - 30*(1680*Mt*pow3(s2t)*pow5(Mst1)*pow6(Mst2) +
        (-768*pow2(Mst1)*pow2(Mt)*pow2(s2t) + 1536*pow4(Mt) - 79*pow4(Mst1)*
        pow4(s2t))*pow8(Mst2)) + pow4(Msq)*(1120*pow3(Mst1)*(8376*s2t*pow2(
        Mst2)*pow3(Mt) + 859*Mt*pow3(s2t)*pow4(Mst2)) - pow4(Mst1)*(1479008*
        pow2(Mst2)*pow2(Mt)*pow2(s2t) + 2834624*pow4(Mt) + 256523*pow4(Mst2)*
        pow4(s2t)) + 560*(37192*s2t*pow3(Mt) + 2579*Mt*pow2(Mst2)*pow3(s2t))*
        pow5(Mst1) + 3360*Mst1*(-128*s2t*pow3(Mt)*pow4(Mst2) + 797*Mt*pow3(s2t)
        *pow6(Mst2)) + 224*pow2(Mst1)*(13803*pow2(Mt)*pow2(s2t)*pow4(Mst2) +
        6584*pow2(Mst2)*pow4(Mt) + 2529*pow4(s2t)*pow6(Mst2)) + 210*(32*pow4(
        Mst2)*pow4(Mt) - 6112*pow2(Mt)*pow2(s2t)*pow6(Mst2) + 2193*pow4(s2t)*
        pow8(Mst2))))) + 3840*Cbeta*Mt*MuSUSY*pow3(s2t)*power10(Mst2))))/(2520.
        *pow2(Mst1)*pow2(Sbeta)*pow4(Msq)*pow4(Mst2)*pow4(Mt)) - (log(pow2(
        Mst1)/pow2(Mst2))*(-196*Dmglst1*Mst1*(-3000*Mst1*pow2(Msq)*(96*Mt*s2t*
        pow2(Mst2)*(37*pow2(Mt) + 2*pow2(Mst2)*pow2(s2t))*pow2(Sbeta)*pow4(
        Mst1) + 96*Mt*s2t*pow2(Mst1)*(-16*pow2(Mt) + 11*pow2(Mst2)*pow2(s2t))*
        pow2(Sbeta)*pow4(Mst2) + 28*pow2(s2t)*(12*pow2(Mt) + 7*pow2(Mst2)*pow2(
        s2t))*pow2(Sbeta)*pow3(Mst1)*pow4(Mst2) + 3*Mst1*pow2(Sbeta)*pow4(Mst2)
        *(48*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 16*pow4(Mt) - 23*pow4(Mst2)*pow4(
        s2t)) + 4*Mt*MuSUSY*s2t*(504*Cbeta*Sbeta*pow2(Mt) + 129*Cbeta*Sbeta*
        pow2(Mst2)*pow2(s2t) + 374*Mt*MuSUSY*s2t*(-1 + pow2(Sbeta)))*pow5(Mst1)
        + 128*s2t*pow2(Sbeta)*pow3(Mt)*pow6(Mst2)) - 750*Mst1*pow2(Mst2)*(84*
        Mt*s2t*pow2(Mst2)*(-20*pow2(Mt) + 23*pow2(Mst2)*pow2(s2t))*pow2(Sbeta)*
        pow4(Mst1) - 24*Mt*s2t*pow2(Mst1)*(-68*pow2(Mt) + pow2(Mst2)*pow2(s2t))
        *pow2(Sbeta)*pow4(Mst2) + 48*Mst1*pow2(Mt)*(23*pow2(Mt) + 3*pow2(Mst2)*
        pow2(s2t))*pow2(Sbeta)*pow4(Mst2) + pow2(Mst2)*pow2(Sbeta)*pow3(Mst1)*(
        336*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 880*pow4(Mt) - 167*pow4(Mst2)*pow4(
        s2t)) + 8*Mt*MuSUSY*s2t*(252*Cbeta*Sbeta*pow2(Mt) + 41*Cbeta*Sbeta*
        pow2(Mst2)*pow2(s2t) - 85*Mt*MuSUSY*s2t*(-1 + pow2(Sbeta)))*pow5(Mst1)
        + 464*s2t*pow2(Sbeta)*pow3(Mt)*pow6(Mst2)) + Sbeta*pow4(Msq)*(4*Cbeta*
        Mt*MuSUSY*(-4795361 + 10152000*z2)*pow3(s2t)*pow6(Mst1) + 25*Sbeta*(-
        96*pow3(Mst1)*(-18*s2t*(947 + 1100*z2)*pow2(Mst2)*pow3(Mt) + Mt*(2899 +
        3732*z2 - 4158*z3)*pow3(s2t)*pow4(Mst2)) - 8128*pow4(Mst2)*pow4(Mt) +
        pow4(Mst1)*(-32*(-37061 + 17640*z2)*pow2(Mst2)*pow2(Mt)*pow2(s2t) -
        512*(2152 + 2079*z2)*pow4(Mt) + (31507 - 21168*z2 + 19872*z3)*pow4(
        Mst2)*pow4(s2t)) + 4*(256*s2t*(5696 + 3987*z2)*pow3(Mt) + 3*Mt*(26687 -
        20496*z2)*pow2(Mst2)*pow3(s2t))*pow5(Mst1) + 128*(79 + 45*z2 + 81*z3)*
        pow2(Mt)*pow2(s2t)*pow2(-1 + pow2(Sbeta))*pow6(Mst1) + 9216*pow2(Mt)*
        pow2(s2t)*pow6(Mst2) + 16*Mst1*(-8*s2t*(145 + 324*z2 - 216*z3)*pow3(Mt)
        *pow4(Mst2) + 27*Mt*(777 + 128*z2 - 12*z3)*pow3(s2t)*pow6(Mst2)) + 16*
        pow2(Mst1)*(2*(21991 + 2736*z2 + 432*z3)*pow2(Mt)*pow2(s2t)*pow4(Mst2)
        + 8*(2069 + 1296*z2)*pow2(Mst2)*pow4(Mt) + 27*(57 + 105*z2 - 46*z3)*
        pow4(s2t)*pow6(Mst2)) - 4068*pow4(s2t)*pow8(Mst2)))) + 98*pow2(Dmglst1)
        *(-4800*Mst1*Mt*s2t*(-518*pow2(Mt) + 305*pow2(Mst2)*pow2(s2t))*pow2(
        Sbeta)*pow4(Msq)*pow4(Mst2) + 3200*Mt*s2t*pow2(Mst2)*pow2(Sbeta)*pow3(
        Mst1)*((-((28427 + 29268*z2)*pow2(Mt)) + 27*(418 + 249*z2 - 234*z3)*
        pow2(Mst2)*pow2(s2t))*pow4(Msq) + 765*pow2(Mt)*pow4(Mst2) + 90*pow2(
        Msq)*(-32*pow2(Mst2)*pow2(Mt) + 23*pow2(s2t)*pow4(Mst2))) + 36*pow2(
        Sbeta)*pow4(Msq)*pow4(Mst2)*(-16800*pow2(Mst2)*pow2(Mt)*pow2(s2t) -
        5632*pow4(Mt) + 5475*pow4(Mst2)*pow4(s2t)) - 400*Mt*s2t*pow2(Sbeta)*((
        8*(181409 + 122004*z2)*pow2(Mt) + (87319 - 50544*z2)*pow2(Mst2)*pow2(
        s2t))*pow4(Msq) + 90*(140*pow2(Mt) - 173*pow2(Mst2)*pow2(s2t))*pow4(
        Mst2) + 720*pow2(Msq)*(-208*pow2(Mst2)*pow2(Mt) + 33*pow2(s2t)*pow4(
        Mst2)))*pow5(Mst1) + 60*Mt*s2t*(300*MuSUSY*pow2(Mst2)*(420*Cbeta*Sbeta*
        pow2(Mt) + 90*Cbeta*Sbeta*pow2(Mst2)*pow2(s2t) + Mt*MuSUSY*s2t*(-1 +
        pow2(Sbeta))) + 1800*MuSUSY*pow2(Msq)*(280*Cbeta*Sbeta*pow2(Mt) + 85*
        Cbeta*Sbeta*pow2(Mst2)*pow2(s2t) + 254*Mt*MuSUSY*s2t*(-1 + pow2(Sbeta))
        ) - s2t*Sbeta*(Cbeta*MuSUSY*s2t*(-1127077 + 2913600*z2 + 23040*z3) +
        480*Mt*Sbeta*(55 + 33*z2 + 9*z3)*pow2(-1 + pow2(Sbeta)))*pow4(Msq))*
        pow6(Mst1) - 200*pow2(Mst1)*pow2(Mst2)*pow2(Sbeta)*(2*pow4(Msq)*(2*(
        18367 + 2736*z2 + 432*z3)*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 72*(245 +
        144*z2)*pow4(Mt) + 9*(286 + 315*z2 - 138*z3)*pow4(Mst2)*pow4(s2t)) -
        180*(23*pow4(Mst2)*pow4(Mt) + 3*pow2(Mt)*pow2(s2t)*pow6(Mst2)) + 45*
        pow2(Msq)*(-48*pow2(Mt)*pow2(s2t)*pow4(Mst2) - 16*pow2(Mst2)*pow4(Mt) +
        23*pow4(s2t)*pow6(Mst2))) + 72000*MuSUSY*pow2(Mt)*(184*Cbeta*Sbeta*
        pow2(Mt) + 3*Cbeta*Sbeta*(-424*pow2(Msq) + 57*pow2(Mst2))*pow2(s2t) +
        574*Mt*MuSUSY*s2t*(-1 + pow2(Sbeta)))*pow7(Mst1) + 25*pow2(Sbeta)*pow4(
        Mst1)*(pow4(Msq)*(96*(-863 + 22392*z2)*pow2(Mst2)*pow2(Mt)*pow2(s2t) +
        64*(46463 + 52992*z2)*pow4(Mt) + 3*(-50563 + 7440*z2 - 11232*z3)*pow4(
        Mst2)*pow4(s2t)) + 1440*pow2(Msq)*(108*pow2(Mt)*pow2(s2t)*pow4(Mst2) +
        43*pow4(s2t)*pow6(Mst2)) - 90*(-880*pow4(Mst2)*pow4(Mt) - 432*pow2(Mt)*
        pow2(s2t)*pow6(Mst2) + 179*pow4(s2t)*pow8(Mst2)))) + Sbeta*(Sbeta*(30*
        pow2(Mst1)*pow4(Mst2)*(2273600*Mst1*s2t*pow3(Mt)*pow4(Mst2) - 39200*
        pow3(Mst1)*(-68*s2t*pow2(Mst2)*pow3(Mt) + 3*Mt*pow3(s2t)*pow4(Mst2)) +
        2103888*pow4(Mst2)*pow4(Mt) + 96*pow2(Mst1)*(1868*pow2(Mt)*pow2(s2t)*
        pow4(Mst2) + 28175*pow2(Mst2)*pow4(Mt)) + pow4(Mst1)*(531552*pow2(Mst2)
        *pow2(Mt)*pow2(s2t) + 1078000*pow4(Mt) - 215819*pow4(Mst2)*pow4(s2t)) -
        58800*(28*s2t*pow3(Mt) - 29*Mt*pow2(Mst2)*pow3(s2t))*pow5(Mst1)) -
        10584000*(-1 + 2*z2)*pow2(Mst1)*(pow2(Mst1) - pow2(Mst2))*pow4(Mst2)*
        pow4(s2t)*pow6(Msq) + 11760*pow2(Msq)*pow2(Mst1)*pow2(Mst2)*(6400*Mst1*
        s2t*pow3(Mt)*pow4(Mst2) + 3200*pow3(Mst1)*(-8*s2t*pow2(Mst2)*pow3(Mt) +
        5*Mt*pow3(s2t)*pow4(Mst2)) + 7712*pow4(Mst2)*pow4(Mt) + 2*pow4(Mst1)*(
        2000*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 991*pow4(Mst2)*pow4(s2t)) + 1600*(
        15*s2t*pow3(Mt) + 7*Mt*pow2(Mst2)*pow3(s2t))*pow5(Mst1) - 3344*pow2(Mt)
        *pow2(s2t)*pow6(Mst2) + pow2(Mst1)*(5344*pow2(Mt)*pow2(s2t)*pow4(Mst2)
        + 1200*pow2(Mst2)*pow4(Mt) - 1167*pow4(s2t)*pow6(Mst2))) - pow4(Msq)*(-
        156800*(-6*s2t*(2197 + 3396*z2)*pow2(Mst2)*pow3(Mt) + Mt*(-1667 + 2232*
        z2 - 4050*z3)*pow3(s2t)*pow4(Mst2))*pow5(Mst1) - 3*(-192*(4731149 +
        4900*S2 - 808500*z2)*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 128*(3291437 +
        3013500*z2)*pow4(Mt) + 3*(1737319 + 431200*S2 - 2018800*z2 - 4939200*
        z3)*pow4(Mst2)*pow4(s2t))*pow6(Mst1) + 705600*pow3(Mst1)*(16*s2t*(13 -
        18*z2 + 12*z3)*pow3(Mt)*pow4(Mst2) + Mt*(2135 + 384*z2 - 36*z3)*pow3(
        s2t)*pow6(Mst2)) + 1568*pow4(Mst1)*(-36*(-33703 + 375*S2 - 2000*z2)*
        pow2(Mt)*pow2(s2t)*pow4(Mst2) + 8*(54223 + 34200*z2)*pow2(Mst2)*pow4(
        Mt) + 9*(1661 + 3750*S2 + 8375*z2 - 3150*z3)*pow4(s2t)*pow6(Mst2)) +
        58800*(32*s2t*(2903 + 2196*z2)*pow3(Mt) + 3*Mt*(1257 - 1136*z2)*pow2(
        Mst2)*pow3(s2t))*pow7(Mst1) + 8467200*(-2*pow4(Mt)*pow6(Mst2) + pow2(
        Mt)*pow2(s2t)*pow8(Mst2)) + 14700*pow2(Mst1)*(96*(149 + 24*z3)*pow4(
        Mst2)*pow4(Mt) - 32*(733 + 81*S2 + 18*z2)*pow2(Mt)*pow2(s2t)*pow6(Mst2)
        + (14209 + 648*S2 - 1788*z2)*pow4(s2t)*pow8(Mst2)))) + 3755280*Cbeta*
        Mt*MuSUSY*pow2(Mst1)*pow3(s2t)*power10(Mst2))))/(529200.*pow2(Mst1)*
        pow2(Sbeta)*pow4(Msq)*pow4(Mst2)*pow4(Mt)) - (log(pow2(Mst1)/pow2(Msq))
        *(22050*pow2(Mst1)*pow2(Sbeta)*pow2(log(pow2(Mst1)/pow2(Mst2)))*pow4(
        Msq)*(pow3(Mst1)*(256*s2t*pow2(Mst2)*pow3(Mt) - 48*Mt*pow3(s2t)*pow4(
        Mst2)) - 5*pow4(Mst1)*(32*pow4(Mt) + pow4(Mst2)*pow4(s2t)) + 384*s2t*
        pow3(Mt)*pow5(Mst1) - 16*pow2(Mt)*pow2(s2t)*pow6(Mst2) + 16*Mst1*Mt*
        pow3(s2t)*pow6(Mst2) + 10*pow2(Mst1)*pow4(s2t)*pow6(Mst2) + 6*pow2(
        Dmglst1)*(640*s2t*pow3(Mst1)*pow3(Mt) - 32*Mst1*(-4*s2t*pow2(Mst2)*
        pow3(Mt) + Mt*pow3(s2t)*pow4(Mst2)) - pow2(Mst1)*(160*pow4(Mt) + pow4(
        Mst2)*pow4(s2t)) + pow4(s2t)*pow6(Mst2)) + 4*Dmglst1*(480*s2t*pow3(Mt)*
        pow4(Mst1) + pow2(Mst1)*(192*s2t*pow2(Mst2)*pow3(Mt) - 44*Mt*pow3(s2t)*
        pow4(Mst2)) - pow3(Mst1)*(160*pow4(Mt) + 3*pow4(Mst2)*pow4(s2t)) + 4*
        Mt*pow3(s2t)*pow6(Mst2) + 3*Mst1*pow4(s2t)*pow6(Mst2)) + 3*pow4(s2t)*
        pow8(Mst2)) - 2450*Dmglst1*Mst1*(32*Mst1*Mt*s2t*pow2(Sbeta)*pow4(Mst2)*
        (-14*pow2(Msq)*pow2(Mst2)*pow2(Mt) + (251*pow2(Mt) - 9*(9 + z2)*pow2(
        Mst2)*pow2(s2t))*pow4(Msq) + 4*pow2(Mt)*pow4(Mst2)) + 2*pow2(Sbeta)*
        pow4(Msq)*pow4(Mst2)*(-192*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 80*pow4(Mt)
        + 45*pow4(Mst2)*pow4(s2t)) + 24*Mt*MuSUSY*s2t*(pow2(Mst2)*(42*Cbeta*
        Sbeta*pow2(Mt) - 20*Cbeta*Sbeta*pow2(Msq)*pow2(s2t) + 21*Mt*MuSUSY*s2t*
        (-1 + pow2(Sbeta))) + 9*Cbeta*Sbeta*pow2(s2t)*pow4(Msq) + 10*Cbeta*
        Sbeta*pow2(s2t)*pow4(Mst2))*pow6(Mst1) - 4*Mt*s2t*pow2(Sbeta)*pow5(
        Mst1)*(18*(40*(-1 + 12*z2)*pow2(Mt) + 43*pow2(Mst2)*pow2(s2t))*pow4(
        Msq) + 80*pow2(Mt)*pow4(Mst2) + 288*pow2(Msq)*pow2(s2t)*pow4(Mst2) -
        183*pow2(s2t)*pow6(Mst2)) + 4*Mt*s2t*pow2(Mst2)*pow2(Sbeta)*pow3(Mst1)*
        (-72*(2*(7 + 24*z2)*pow2(Mt) - (16 + 11*z2)*pow2(Mst2)*pow2(s2t))*pow4(
        Msq) + 48*pow2(Mt)*pow4(Mst2) + 16*pow2(Msq)*(7*pow2(Mst2)*pow2(Mt) +
        18*pow2(s2t)*pow4(Mst2)) - 3*pow2(s2t)*pow6(Mst2)) + 6*pow2(Mst1)*pow2(
        Mst2)*pow2(Sbeta)*(64*pow4(Mst2)*pow4(Mt) - pow4(Msq)*(152*pow2(Mst2)*
        pow2(Mt)*pow2(s2t) + 432*pow4(Mt) + 9*(7 + 4*z2)*pow4(Mst2)*pow4(s2t))
        + 12*pow2(Mt)*pow2(s2t)*pow6(Mst2) + 8*pow2(Msq)*(-20*pow2(Mt)*pow2(
        s2t)*pow4(Mst2) + 24*pow2(Mst2)*pow4(Mt) + pow4(s2t)*pow6(Mst2))) +
        pow2(Sbeta)*pow4(Mst1)*(2192*pow4(Mst2)*pow4(Mt) + 18*pow4(Msq)*(-72*
        pow2(Mst2)*pow2(Mt)*pow2(s2t) + 32*(9 + 20*z2)*pow4(Mt) + 3*(9 + 4*z2)*
        pow4(Mst2)*pow4(s2t)) - 648*pow2(Mt)*pow2(s2t)*pow6(Mst2) + 24*pow2(
        Msq)*(80*pow2(Mt)*pow2(s2t)*pow4(Mst2) + 80*pow2(Mst2)*pow4(Mt) + pow4(
        s2t)*pow6(Mst2)) + 3*pow4(s2t)*pow8(Mst2))) - 147*pow2(Dmglst1)*(-800*
        Mst1*Mt*s2t*(-112*pow2(Mt) + 17*pow2(Mst2)*pow2(s2t))*pow2(Sbeta)*pow4(
        Msq)*pow4(Mst2) + 800*Mt*s2t*pow2(Mst2)*pow2(Sbeta)*pow3(Mst1)*((-12*(
        13 + 24*z2)*pow2(Mt) + (155 + 72*z2)*pow2(Mst2)*pow2(s2t))*pow4(Msq) +
        4*pow2(Mt)*pow4(Mst2) + 24*pow2(Msq)*pow2(s2t)*pow4(Mst2)) + 2*pow2(
        Sbeta)*pow4(Msq)*pow4(Mst2)*(-7200*pow2(Mst2)*pow2(Mt)*pow2(s2t) +
        12032*pow4(Mt) + 925*pow4(Mst2)*pow4(s2t)) + 600*Mt*MuSUSY*s2t*(5*pow2(
        Mst2)*(14*Cbeta*Sbeta*pow2(Mt) - 4*Cbeta*Sbeta*pow2(Msq)*pow2(s2t) + 7*
        Mt*MuSUSY*s2t*(-1 + pow2(Sbeta))) - 29*Cbeta*Sbeta*pow2(s2t)*pow4(Msq)
        + 10*Cbeta*Sbeta*pow2(s2t)*pow4(Mst2))*pow6(Mst1) - 1600*Mt*s2t*pow2(
        Sbeta)*pow5(Mst1)*(6*(5*(13 + 24*z2)*pow2(Mt) + 2*pow2(Mst2)*pow2(s2t))
        *pow4(Msq) + 2*pow2(Mt)*pow4(Mst2) + 12*pow2(Msq)*pow2(s2t)*pow4(Mst2)
        - 15*pow2(s2t)*pow6(Mst2)) + 50*pow2(Mst1)*pow2(Mst2)*pow2(Sbeta)*(64*
        pow4(Mst2)*pow4(Mt) - pow4(Msq)*(-72*pow2(Mst2)*pow2(Mt)*pow2(s2t) +
        432*pow4(Mt) + (203 + 36*z2)*pow4(Mst2)*pow4(s2t)) + 12*pow2(Mt)*pow2(
        s2t)*pow6(Mst2) + 16*pow2(Msq)*(-30*pow2(Mt)*pow2(s2t)*pow4(Mst2) + 52*
        pow2(Mst2)*pow4(Mt) + 3*pow4(s2t)*pow6(Mst2))) + 96000*s2t*pow2(MuSUSY)
        *(-1 + pow2(Sbeta))*pow3(Mt)*pow7(Mst1) + 25*pow2(Sbeta)*pow4(Mst1)*(
        2608*pow4(Mst2)*pow4(Mt) + 2*pow4(Msq)*(-2040*pow2(Mst2)*pow2(Mt)*pow2(
        s2t) + 96*(47 + 60*z2)*pow4(Mt) + (67 + 36*z2)*pow4(Mst2)*pow4(s2t)) -
        888*pow2(Mt)*pow2(s2t)*pow6(Mst2) - 24*pow2(Msq)*(-80*pow2(Mt)*pow2(
        s2t)*pow4(Mst2) - 80*pow2(Mst2)*pow4(Mt) + 3*pow4(s2t)*pow6(Mst2)) +
        45*pow4(s2t)*pow8(Mst2))) - Sbeta*(Sbeta*(8*pow2(Mst1)*pow4(Mst2)*(
        39200*Mst1*s2t*pow3(Mt)*pow4(Mst2) - 1225*pow3(Mst1)*(-16*s2t*pow2(
        Mst2)*pow3(Mt) + 3*Mt*pow3(s2t)*pow4(Mst2)) + 49209*pow4(Mst2)*pow4(Mt)
        + 15*pow2(Mst1)*(991*pow2(Mt)*pow2(s2t)*pow4(Mst2) + 3920*pow2(Mst2)*
        pow4(Mt)) + pow4(Mst1)*(-43935*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 220709*
        pow4(Mt) - 3459*pow4(Mst2)*pow4(s2t)) - 3675*(16*s2t*pow3(Mt) - 13*Mt*
        pow2(Mst2)*pow3(s2t))*pow5(Mst1)) - 88200*pow2(Mst2)*pow6(Msq)*(pow2(-
        4*Mst2*pow2(Mt) + pow2(s2t)*pow3(Mst2)) - pow4(Mst1)*(8*pow2(Mt)*pow2(
        s2t) + pow2(Mst2)*pow4(s2t)) + pow2(Mst1)*(16*pow2(Mst2)*pow2(Mt)*pow2(
        s2t) + 16*pow4(Mt) - pow4(Mst2)*pow4(s2t)) + pow4(s2t)*pow6(Mst1)) -
        196*pow2(Msq)*pow2(Mst1)*pow2(Mst2)*(5600*Mst1*s2t*pow3(Mt)*pow4(Mst2)
        - 800*pow3(Mst1)*(7*s2t*pow2(Mst2)*pow3(Mt) + 6*Mt*pow3(s2t)*pow4(Mst2)
        ) - 656*pow4(Mst2)*pow4(Mt) - pow4(Mst1)*(5128*pow2(Mst2)*pow2(Mt)*
        pow2(s2t) + 6000*pow4(Mt) + 345*pow4(Mst2)*pow4(s2t)) + 4800*Mt*pow2(
        Mst2)*pow3(s2t)*pow5(Mst1) + 872*pow2(Mt)*pow2(s2t)*pow6(Mst2) + pow2(
        Mst1)*(1256*pow2(Mt)*pow2(s2t)*pow4(Mst2) - 1856*pow2(Mst2)*pow4(Mt) +
        255*pow4(s2t)*pow6(Mst2))) - 1225*pow2(Mst1)*pow4(Msq)*(-576*Mt*s2t*
        pow2(Mst2)*(-2*(1 + 8*z2)*pow2(Mt) + (7 + 3*z2)*pow2(Mst2)*pow2(s2t))*
        pow3(Mst1) + 576*Mst1*Mt*s2t*(-11*pow2(Mt) + (6 + z2)*pow2(Mst2)*pow2(
        s2t))*pow4(Mst2) - 2*pow4(Mst2)*(72*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 32*
        (-41 + 36*z2)*pow4(Mt) + 3*(-37 + 6*z2)*pow4(Mst2)*pow4(s2t)) - 3*pow4(
        Mst1)*(240*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 384*(1 + 5*z2)*pow4(Mt) + 5*
        (11 + 12*z2)*pow4(Mst2)*pow4(s2t)) + 144*Mt*s2t*(8*(-7 + 12*z2)*pow2(
        Mt) + 11*pow2(Mst2)*pow2(s2t))*pow5(Mst1) + 6*pow2(Mst1)*(240*pow2(Mt)*
        pow2(s2t)*pow4(Mst2) + 432*pow2(Mst2)*pow4(Mt) + (-53 + 36*z2)*pow4(
        s2t)*pow6(Mst2)))) + 19142*Cbeta*Mt*MuSUSY*pow2(Mst1)*pow3(s2t)*
        power10(Mst2)) + 105*log(pow2(Mst1)/pow2(Mst2))*pow2(Mst1)*(140*
        Dmglst1*(-60*Mt*s2t*(-4*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow2(Sbeta)*
        pow4(Mst1)*pow4(Mst2) + 5*pow2(Sbeta)*pow3(Mst1)*pow4(Mst2)*(16*pow4(
        Mt) + pow4(Mst2)*pow4(s2t)) - 8*Mt*MuSUSY*pow2(Mst2)*pow2(s2t)*(8*
        Cbeta*s2t*Sbeta*pow2(Mst2) + 11*Mt*MuSUSY*(-1 + pow2(Sbeta)))*pow5(
        Mst1) + 96*s2t*pow2(Mst1)*pow2(Sbeta)*pow3(Mt)*pow6(Mst2) + 48*Mst1*
        pow2(Sbeta)*pow4(Mt)*pow6(Mst2) + pow2(Msq)*(-96*Mt*s2t*pow2(Mst1)*(-4*
        pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow2(Sbeta)*pow4(Mst2) - 96*Mt*pow2(
        Sbeta)*pow3(s2t)*pow4(Mst1)*pow4(Mst2) + 6*Mst1*pow2(Sbeta)*pow4(Mst2)*
        (16*pow4(Mt) + pow4(Mst2)*pow4(s2t)) - 224*pow2(Mt)*pow2(MuSUSY)*pow2(
        s2t)*(-1 + pow2(Sbeta))*pow5(Mst1) + 64*s2t*pow2(Sbeta)*pow3(Mt)*pow6(
        Mst2) - 40*pow2(Sbeta)*pow3(Mst1)*pow4(s2t)*pow6(Mst2)) + 12*Sbeta*
        pow4(Msq)*(18*Cbeta*Mt*MuSUSY*pow3(s2t)*pow5(Mst1) + Sbeta*((776*s2t*
        pow3(Mt) - 66*Mt*pow2(Mst2)*pow3(s2t))*pow4(Mst1) + 2*Mt*s2t*(-8*pow2(
        Mt) + 7*pow2(Mst2)*pow2(s2t))*pow4(Mst2) + pow2(Mst1)*(160*s2t*pow2(
        Mst2)*pow3(Mt) - 26*Mt*pow3(s2t)*pow4(Mst2)) + pow3(Mst1)*(28*pow2(
        Mst2)*pow2(Mt)*pow2(s2t) - 232*pow4(Mt) + 3*pow4(Mst2)*pow4(s2t)) + 3*
        Mst1*(4*pow2(Mt)*pow2(s2t)*pow4(Mst2) + pow4(s2t)*pow6(Mst2)))) + 16*
        s2t*pow2(Sbeta)*pow3(Mt)*pow8(Mst2)) + 210*pow2(Dmglst1)*(8*Mt*MuSUSY*
        pow2(s2t)*(-5*pow2(Mst2)*(3*Cbeta*s2t*Sbeta*pow2(Mst2) + 5*Mt*MuSUSY*(-
        1 + pow2(Sbeta))) - 36*Mt*MuSUSY*pow2(Msq)*(-1 + pow2(Sbeta)) + 55*
        Cbeta*s2t*Sbeta*pow4(Msq))*pow4(Mst1) + 64*Mst1*Mt*s2t*pow2(Mst2)*pow2(
        Sbeta)*(2*(16*pow2(Mt) - 5*pow2(Mst2)*pow2(s2t))*pow4(Msq) + pow2(Mt)*
        pow4(Mst2) + pow2(Msq)*(4*pow2(Mst2)*pow2(Mt) - pow2(s2t)*pow4(Mst2)))
        + 2*pow2(Sbeta)*pow4(Mst2)*(8*pow2(Mst2)*pow4(Mt) + 6*pow4(Msq)*(4*
        pow2(Mt)*pow2(s2t) + pow2(Mst2)*pow4(s2t)) + pow2(Msq)*(16*pow4(Mt) +
        pow4(Mst2)*pow4(s2t))) + 96*MuSUSY*pow2(Mt)*(4*Cbeta*Sbeta*pow2(Mt) -
        Cbeta*Sbeta*(4*pow2(Msq) + 5*pow2(Mst2))*pow2(s2t) - 10*Mt*MuSUSY*s2t*(
        -1 + pow2(Sbeta)))*pow5(Mst1) + 16*Mt*s2t*pow2(Sbeta)*pow3(Mst1)*((896*
        pow2(Mt) - 78*pow2(Mst2)*pow2(s2t))*pow4(Msq) + 20*pow2(Mt)*pow4(Mst2)
        - 4*pow2(Msq)*pow2(s2t)*pow4(Mst2) - 5*pow2(s2t)*pow6(Mst2)) + pow2(
        Mst1)*pow2(Sbeta)*(80*pow4(Mst2)*pow4(Mt) - 8*pow4(Msq)*(-54*pow2(Mst2)
        *pow2(Mt)*pow2(s2t) + 428*pow4(Mt) - 9*pow4(Mst2)*pow4(s2t)) - 40*pow2(
        Msq)*pow4(s2t)*pow6(Mst2) + 5*pow4(s2t)*pow8(Mst2))) + Sbeta*(Sbeta*(
        28*pow2(Msq)*pow4(Mst2)*(320*Mst1*s2t*pow2(Mst2)*pow3(Mt) + 160*pow3(
        Mst1)*(4*s2t*pow3(Mt) - Mt*pow2(Mst2)*pow3(s2t)) + 52*pow2(Mt)*pow2(
        s2t)*pow4(Mst2) + 64*pow2(Mst2)*pow4(Mt) - pow4(Mst1)*(20*pow2(Mt)*
        pow2(s2t) + 41*pow2(Mst2)*pow4(s2t)) + pow2(Mst1)*(-32*pow2(Mst2)*pow2(
        Mt)*pow2(s2t) + 240*pow4(Mt) + 6*pow4(Mst2)*pow4(s2t)) - 160*Mt*pow3(
        s2t)*pow5(Mst1)) + pow4(Mst2)*(4480*s2t*pow2(Mst2)*pow3(Mst1)*pow3(Mt)
        + 2240*Mst1*s2t*pow3(Mt)*pow4(Mst2) + 144*pow4(Mst2)*pow4(Mt) + 24*
        pow2(Mst1)*(11*pow2(Mt)*pow2(s2t)*pow4(Mst2) + 140*pow2(Mst2)*pow4(Mt))
        + pow4(Mst1)*(-504*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 2800*pow4(Mt) + 163*
        pow4(Mst2)*pow4(s2t)) + 1680*(4*s2t*pow3(Mt) - Mt*pow2(Mst2)*pow3(s2t))
        *pow5(Mst1)) + 1680*(-pow2(Mst1) + pow2(Mst2))*pow4(Mst2)*pow4(s2t)*
        pow6(Msq) + 70*pow4(Msq)*(48*pow3(Mst1)*(16*s2t*pow2(Mst2)*pow3(Mt) +
        3*Mt*pow3(s2t)*pow4(Mst2)) - 80*pow4(Mst2)*pow4(Mt) - 32*pow4(Mst1)*(-
        3*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 36*pow4(Mt) + pow4(Mst2)*pow4(s2t)) +
        48*(68*s2t*pow3(Mt) - 5*Mt*pow2(Mst2)*pow3(s2t))*pow5(Mst1) - 112*pow2(
        Mt)*pow2(s2t)*pow6(Mst2) + 48*Mst1*(-8*s2t*pow3(Mt)*pow4(Mst2) + 7*Mt*
        pow3(s2t)*pow6(Mst2)) + 8*pow2(Mst1)*(32*pow2(Mt)*pow2(s2t)*pow4(Mst2)
        + 3*pow4(s2t)*pow6(Mst2)) + 35*pow4(s2t)*pow8(Mst2))) - 72*Cbeta*Mt*
        MuSUSY*pow3(s2t)*power10(Mst2)))))/(4410.*pow2(Mst1)*pow2(Sbeta)*pow4(
        Msq)*pow4(Mst2)*pow4(Mt));

   return result;
}

/**
 * @return returns the susy log^1 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double H5::calc_coef_at_as2_no_sm_logs_log1() const {

   const double result =
      -447.7037037037037 - (11024*Dmglst1)/(27.*Mst1) - (40192*Dmglst1*s2t)/(
        27.*Mt) + (1472*Mst1*s2t)/(3.*Mt) + (1328*z2)/3. + (896*Dmglst1*s2t*z2)
        /(3.*Mt) + (896*Mst1*s2t*z2)/(3.*Mt) + 128*z3 - (80192*s2t*pow2(
        Dmglst1))/(27.*Mst1*Mt) - (160*Dmglst1*Mst1)/pow2(Msq) + (368*pow2(
        Dmglst1))/(3.*pow2(Msq)) - (8960*Mst1*s2t*pow2(Dmglst1))/(3.*Mt*pow2(
        Msq)) - (102592*pow2(Dmglst1))/(225.*pow2(Mst1)) - (320*pow2(Msq))/
        pow2(Mst1) - (880*pow2(Mst1))/(3.*pow2(Msq)) - (5120*Dmglst1*s2t*pow2(
        Mst1))/(3.*Mt*pow2(Msq)) + (19648*Dmglst1*Mst1)/(9.*pow2(Mst2)) + (
        1536*Dmglst1*Mst1*z2)/pow2(Mst2) + (208336*pow2(Dmglst1))/(135.*pow2(
        Mst2)) + (139712*Mst1*s2t*pow2(Dmglst1))/(9.*Mt*pow2(Mst2)) + (768*z2*
        pow2(Dmglst1))/pow2(Mst2) + (8256*Mst1*s2t*z2*pow2(Dmglst1))/(Mt*pow2(
        Mst2)) - (320*pow2(Msq))/pow2(Mst2) + (448*pow2(Mst1))/(3.*pow2(Mst2))
        + (32768*Dmglst1*s2t*pow2(Mst1))/(3.*Mt*pow2(Mst2)) + (2368*z2*pow2(
        Mst1))/(3.*pow2(Mst2)) + (25216*Dmglst1*s2t*z2*pow2(Mst1))/(3.*Mt*pow2(
        Mst2)) + (1600*pow2(Dmglst1)*pow2(Mst1))/(pow2(Msq)*pow2(Mst2)) - (
        1280*pow2(Mst2))/(9.*pow2(Msq)) - (3520*Dmglst1*s2t*pow2(Mst2))/(9.*Mt*
        pow2(Msq)) - (3520*Mst1*s2t*pow2(Mst2))/(9.*Mt*pow2(Msq)) - (32*pow2(
        Mst2))/pow2(Mst1) + (52816*Dmglst1*Mst1*pow2(s2t))/(9.*pow2(Mt)) + (
        1024*Dmglst1*Mst1*z2*pow2(s2t))/pow2(Mt) + (33736*pow2(Dmglst1)*pow2(
        s2t))/(9.*pow2(Mt)) + (512*z2*pow2(Dmglst1)*pow2(s2t))/pow2(Mt) - (320*
        pow2(Msq)*pow2(s2t))/pow2(Mt) + (8080*pow2(Mst1)*pow2(s2t))/(3.*pow2(
        Mt)) + (512*z2*pow2(Mst1)*pow2(s2t))/pow2(Mt) - (12176*pow2(Dmglst1)*
        pow2(Mst1)*pow2(s2t))/(3.*pow2(Mst2)*pow2(Mt)) - (5120*z2*pow2(Dmglst1)
        *pow2(Mst1)*pow2(s2t))/(pow2(Mst2)*pow2(Mt)) + (160*pow2(Msq)*pow2(
        Mst1)*pow2(s2t))/(pow2(Mst2)*pow2(Mt)) - (208*pow2(Mst2)*pow2(s2t))/
        pow2(Mt) - (592*Dmglst1*pow2(Mst2)*pow2(s2t))/(9.*Mst1*pow2(Mt)) + (
        1976*pow2(Dmglst1)*pow2(Mst2)*pow2(s2t))/(9.*pow2(Mst1)*pow2(Mt)) + (
        160*pow2(Msq)*pow2(Mst2)*pow2(s2t))/(pow2(Mst1)*pow2(Mt)) - (1280*s2t*
        pow3(Mst1))/(9.*Mt*pow2(Msq)) + (4096*s2t*pow3(Mst1))/(3.*Mt*pow2(Mst2)
        ) + (8704*s2t*z2*pow3(Mst1))/(3.*Mt*pow2(Mst2)) + (3200*Dmglst1*pow3(
        Mst1))/(3.*pow2(Msq)*pow2(Mst2)) - (23680*s2t*pow2(Dmglst1)*pow3(Mst1))
        /(3.*Mt*pow2(Msq)*pow2(Mst2)) - (2208*Dmglst1*pow2(s2t)*pow3(Mst1))/(
        pow2(Mst2)*pow2(Mt)) - (2048*Dmglst1*z2*pow2(s2t)*pow3(Mst1))/(pow2(
        Mst2)*pow2(Mt)) - (14768*Mst1*pow2(Dmglst1)*pow3(s2t))/(9.*pow3(Mt)) -
        (1152*Mst1*z2*pow2(Dmglst1)*pow3(s2t))/pow3(Mt) - (6688*Dmglst1*pow2(
        Mst1)*pow3(s2t))/(3.*pow3(Mt)) - (2912*Dmglst1*z2*pow2(Mst1)*pow3(s2t))
        /(3.*pow3(Mt)) + (6464*Dmglst1*pow2(Mst2)*pow3(s2t))/(3.*pow3(Mt)) + (
        5344*Mst1*pow2(Mst2)*pow3(s2t))/(3.*pow3(Mt)) + (544*Dmglst1*z2*pow2(
        Mst2)*pow3(s2t))/(3.*pow3(Mt)) + (544*Mst1*z2*pow2(Mst2)*pow3(s2t))/(3.
        *pow3(Mt)) + (272*pow2(Dmglst1)*pow2(Mst2)*pow3(s2t))/(9.*Mst1*pow3(Mt)
        ) - (1792*pow3(Mst1)*pow3(s2t))/pow3(Mt) - (608*z2*pow3(Mst1)*pow3(s2t)
        )/(3.*pow3(Mt)) + (1408*pow2(Dmglst1)*pow3(Mst1)*pow3(s2t))/(pow2(Mst2)
        *pow3(Mt)) + (1052*pow2(Dmglst1)*pow2(Mst1))/(3.*pow4(Msq)) - (280*
        Dmglst1*Mst1*pow2(Mst2))/(3.*pow4(Msq)) - (140*pow2(Dmglst1)*pow2(Mst2)
        )/(3.*pow4(Msq)) - (1040*Mst1*s2t*pow2(Dmglst1)*pow2(Mst2))/(3.*Mt*
        pow4(Msq)) - (140*pow2(Mst1)*pow2(Mst2))/(3.*pow4(Msq)) - (1040*
        Dmglst1*s2t*pow2(Mst1)*pow2(Mst2))/(3.*Mt*pow4(Msq)) + (760*Dmglst1*
        pow3(Mst1))/(9.*pow4(Msq)) - (7120*s2t*pow2(Dmglst1)*pow3(Mst1))/(3.*
        Mt*pow4(Msq)) - (1040*s2t*pow2(Mst2)*pow3(Mst1))/(9.*Mt*pow4(Msq)) + (
        800*pow4(Mst1))/(3.*pow2(Msq)*pow2(Mst2)) - (7040*Dmglst1*s2t*pow4(
        Mst1))/(3.*Mt*pow2(Msq)*pow2(Mst2)) - (1168*pow2(s2t)*pow4(Mst1))/(
        pow2(Mst2)*pow2(Mt)) + (2264*Dmglst1*pow3(s2t)*pow4(Mst1))/(3.*pow2(
        Mst2)*pow3(Mt)) - (410*pow4(Mst1))/(9.*pow4(Msq)) - (2120*Dmglst1*s2t*
        pow4(Mst1))/(3.*Mt*pow4(Msq)) - (931376*pow2(Dmglst1)*pow2(Mst1))/(135.
        *pow4(Mst2)) - (1600*z2*pow2(Dmglst1)*pow2(Mst1))/pow4(Mst2) - (76928*
        Dmglst1*pow3(Mst1))/(9.*pow4(Mst2)) - (256*Dmglst1*z2*pow3(Mst1))/pow4(
        Mst2) + (404864*s2t*pow2(Dmglst1)*pow3(Mst1))/(9.*Mt*pow4(Mst2)) + (
        85696*s2t*z2*pow2(Dmglst1)*pow3(Mst1))/(3.*Mt*pow4(Mst2)) - (2756*pow4(
        Mst1))/pow4(Mst2) + (6496*Dmglst1*s2t*pow4(Mst1))/(Mt*pow4(Mst2)) + (
        256*z2*pow4(Mst1))/pow4(Mst2) + (43904*Dmglst1*s2t*z2*pow4(Mst1))/(3.*
        Mt*pow4(Mst2)) - (172*Cbeta*MuSUSY*pow2(Dmglst1)*pow3(s2t)*pow4(Mst1))/
        (3.*Sbeta*pow3(Mt)*pow4(Mst2)) + (16*pow2(s2t)*pow4(Mst2))/(pow2(Mst1)*
        pow2(Mt)) - (30*pow4(Mst2))/pow4(Msq) - (520*Dmglst1*s2t*pow4(Mst2))/(
        9.*Mt*pow4(Msq)) - (520*Mst1*s2t*pow4(Mst2))/(9.*Mt*pow4(Msq)) - (1336*
        pow2(Dmglst1)*pow2(Mst1)*pow4(s2t))/(9.*pow4(Mt)) - (172*z2*pow2(
        Dmglst1)*pow2(Mst1)*pow4(s2t))/(3.*pow4(Mt)) + (20*pow2(Msq)*pow2(Mst1)
        *pow4(s2t))/pow4(Mt) + (116*Dmglst1*Mst1*pow2(Mst2)*pow4(s2t))/pow4(Mt)
        + (344*Dmglst1*Mst1*z2*pow2(Mst2)*pow4(s2t))/(3.*pow4(Mt)) + (1334*
        pow2(Dmglst1)*pow2(Mst2)*pow4(s2t))/(9.*pow4(Mt)) + (172*z2*pow2(
        Dmglst1)*pow2(Mst2)*pow4(s2t))/(3.*pow4(Mt)) + (20*pow2(Msq)*pow2(Mst2)
        *pow4(s2t))/pow4(Mt) - (330*pow2(Mst1)*pow2(Mst2)*pow4(s2t))/pow4(Mt) +
        (172*z2*pow2(Mst1)*pow2(Mst2)*pow4(s2t))/(3.*pow4(Mt)) - (396*Dmglst1*
        pow3(Mst1)*pow4(s2t))/pow4(Mt) - (344*Dmglst1*z2*pow3(Mst1)*pow4(s2t))/
        (3.*pow4(Mt)) - (45*pow4(Mst1)*pow4(s2t))/(2.*pow4(Mt)) - (30*z2*pow4(
        Mst1)*pow4(s2t))/pow4(Mt) - (20*pow2(Msq)*pow4(Mst1)*pow4(s2t))/(pow2(
        Mst2)*pow4(Mt)) + (228*pow4(Mst2)*pow4(s2t))/pow4(Mt) + (12*Dmglst1*
        pow4(Mst2)*pow4(s2t))/(Mst1*pow4(Mt)) - (82*z2*pow4(Mst2)*pow4(s2t))/(
        3.*pow4(Mt)) - (16*pow2(Dmglst1)*pow4(Mst2)*pow4(s2t))/(9.*pow2(Mst1)*
        pow4(Mt)) - (20*pow2(Msq)*pow4(Mst2)*pow4(s2t))/(pow2(Mst1)*pow4(Mt)) -
        (640*s2t*pow5(Mst1))/(3.*Mt*pow2(Msq)*pow2(Mst2)) - (40*pow3(s2t)*pow5(
        Mst1))/(3.*pow2(Mst2)*pow3(Mt)) - (40*s2t*pow5(Mst1))/(3.*Mt*pow4(Msq))
        - (5984*s2t*pow5(Mst1))/(3.*Mt*pow4(Mst2)) + (8960*s2t*z2*pow5(Mst1))/(
        3.*Mt*pow4(Mst2)) - (88*Cbeta*Dmglst1*MuSUSY*pow3(s2t)*pow5(Mst1))/(3.*
        Sbeta*pow3(Mt)*pow4(Mst2)) - (19840*Cbeta*MuSUSY*pow2(Dmglst1)*pow5(
        Mst1))/(3.*Sbeta*pow4(Msq)*pow4(Mst2)) - (8*log(pow2(Mst1)/pow2(Msq))*(
        60*log(pow2(Mst1)/pow2(Mst2))*pow2(Mst1)*pow4(Msq)*(6*Mst1*pow2(
        Dmglst1)*(-5*Mst1*Mt + 20*s2t*pow2(Mst1) + 4*s2t*pow2(Mst2)) + 8*s2t*
        pow2(Mst2)*pow3(Mst1) - 5*Mt*pow4(Mst1) + Mt*pow4(Mst2) + 4*Mst1*s2t*
        pow4(Mst2) + 4*Dmglst1*(6*s2t*pow2(Mst1)*pow2(Mst2) - 5*Mt*pow3(Mst1) +
        15*s2t*pow4(Mst1) + s2t*pow4(Mst2)) + 12*s2t*pow5(Mst1)) + 5*pow2(Mst1)
        *(-8*pow2(Msq)*(4*Mt*pow2(Mst1) + Mt*pow2(Mst2) + 2*Mst1*s2t*pow2(Mst2)
        - 2*s2t*pow3(Mst1))*pow4(Mst2) - pow4(Mst2)*(6*Mt*pow2(Mst1)*pow2(Mst2)
        + 8*s2t*pow2(Mst2)*pow3(Mst1) + 19*Mt*pow4(Mst1) + 3*Mt*pow4(Mst2) + 4*
        Mst1*s2t*pow4(Mst2) - 12*s2t*pow5(Mst1)) + 6*pow4(Msq)*(-6*Mt*pow2(
        Mst1)*pow2(Mst2) + 8*s2t*pow2(Mst2)*pow3(Mst1) - 14*Mt*pow4(Mst1) - 7*
        Mt*pow4(Mst2) + 12*Mst1*s2t*pow4(Mst2) + 8*s2t*pow5(Mst1))) + 20*
        Dmglst1*Mst1*(4*Mst1*pow2(Msq)*(-3*Mst1*Mt + s2t*pow2(Mst1) - s2t*pow2(
        Mst2))*pow4(Mst2) + Mst1*pow4(Mst2)*(-3*Mst1*Mt*pow2(Mst2) - 6*s2t*
        pow2(Mst1)*pow2(Mst2) - 13*Mt*pow3(Mst1) + 7*s2t*pow4(Mst1) - s2t*pow4(
        Mst2)) + 2*pow4(Msq)*(-9*Mt*pow2(Mst1)*pow2(Mst2) + 42*s2t*pow2(Mst2)*
        pow3(Mst1) - 57*Mt*pow4(Mst1) - 2*Mt*pow4(Mst2) + 17*Mst1*s2t*pow4(
        Mst2) + 66*s2t*pow5(Mst1))) + 6*pow2(Dmglst1)*(-20*Mt*pow2(Msq)*pow2(
        Mst1)*pow4(Mst2) + 5*pow2(Mst1)*(-11*Mt*pow2(Mst1) - Mt*pow2(Mst2) - 4*
        Mst1*s2t*pow2(Mst2) + 4*s2t*pow3(Mst1))*pow4(Mst2) + 2*pow4(Msq)*(-15*
        Mt*pow2(Mst1)*pow2(Mst2) + 260*s2t*pow2(Mst2)*pow3(Mst1) - 385*Mt*pow4(
        Mst1) + Mt*pow4(Mst2) + 20*Mst1*s2t*pow4(Mst2) + 740*s2t*pow5(Mst1)))))
        /(3.*Mt*pow2(Mst1)*pow4(Msq)*pow4(Mst2)) + (2*pow2(log(pow2(Mst1)/pow2(
        Mst2)))*(512*Cbeta*Dmglst1*(Dmglst1 + 2*Mst1)*Mt*MuSUSY*pow3(s2t)*pow4(
        Mst1) + Sbeta*(8*pow3(Mst1)*(296*s2t*pow2(Mst2)*pow3(Mt) + 53*Mt*pow3(
        s2t)*pow4(Mst2)) - 212*pow4(Mst2)*pow4(Mt) + 3*pow4(Mst1)*(-128*pow2(
        Mst2)*pow2(Mt)*pow2(s2t) + 16*pow4(Mt) - 9*pow4(Mst2)*pow4(s2t)) + 32*(
        95*s2t*pow3(Mt) + 16*Mt*pow2(Mst2)*pow3(s2t))*pow5(Mst1) - 492*pow2(Mt)
        *pow2(s2t)*pow6(Mst2) + 24*Mst1*(-8*s2t*pow3(Mt)*pow4(Mst2) + 19*Mt*
        pow3(s2t)*pow6(Mst2)) + pow2(Dmglst1)*(30192*s2t*pow3(Mst1)*pow3(Mt) +
        64*pow2(Mt)*pow2(s2t)*pow4(Mst2) + Mst1*(7536*s2t*pow2(Mst2)*pow3(Mt) -
        864*Mt*pow3(s2t)*pow4(Mst2)) + 528*pow2(Mst2)*pow4(Mt) - 3*pow2(Mst1)*(
        2048*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 688*pow4(Mt) + 9*pow4(Mst2)*pow4(
        s2t)) + 155*pow4(s2t)*pow6(Mst2)) + pow2(Mst1)*(228*pow2(Mt)*pow2(s2t)*
        pow4(Mst2) + 544*pow2(Mst2)*pow4(Mt) + 237*pow4(s2t)*pow6(Mst2)) + 2*
        Dmglst1*(64*(121*s2t*pow3(Mt) + 4*Mt*pow2(Mst2)*pow3(s2t))*pow4(Mst1) -
        96*s2t*pow3(Mt)*pow4(Mst2) + pow2(Mst1)*(3696*s2t*pow2(Mst2)*pow3(Mt) -
        76*Mt*pow3(s2t)*pow4(Mst2)) - 3*pow3(Mst1)*(512*pow2(Mst2)*pow2(Mt)*
        pow2(s2t) + 128*pow4(Mt) + 9*pow4(Mst2)*pow4(s2t)) + 228*Mt*pow3(s2t)*
        pow6(Mst2) + Mst1*(64*pow2(Mt)*pow2(s2t)*pow4(Mst2) + 528*pow2(Mst2)*
        pow4(Mt) + 155*pow4(s2t)*pow6(Mst2))) + 82*pow4(s2t)*pow8(Mst2))))/(3.*
        Sbeta*pow4(Mst2)*pow4(Mt)) - (log(pow2(Mst1)/pow2(Mst2))*(120*Cbeta*
        Dmglst1*Mt*MuSUSY*(-98*Mst1*pow3(s2t)*pow4(Msq) + 5*Dmglst1*(480*Mst1*
        pow3(Mt) - 67*pow3(s2t)*pow4(Msq)))*pow6(Mst1) + Sbeta*(-15*pow2(Mst1)*
        (-160*pow2(Msq)*pow2(Mst2)*pow3(Mt)*(3*Mt*pow2(Mst1)*pow2(Mst2) + 8*
        s2t*pow2(Mst2)*pow3(Mst1) + 2*Mt*pow4(Mst2) + 4*Mst1*s2t*pow4(Mst2) +
        12*s2t*pow5(Mst1)) - 40*pow3(Mt)*pow4(Mst2)*(6*Mt*pow2(Mst1)*pow2(Mst2)
        + 8*s2t*pow2(Mst2)*pow3(Mst1) + 5*Mt*pow4(Mst1) + 3*Mt*pow4(Mst2) + 4*
        Mst1*s2t*pow4(Mst2) + 12*s2t*pow5(Mst1)) + 120*(pow2(Mst1) - pow2(Mst2)
        )*pow4(Mst2)*pow4(s2t)*pow6(Msq) + pow4(Msq)*(64*pow3(Mst1)*(244*s2t*
        pow2(Mst2)*pow3(Mt) + 49*Mt*pow3(s2t)*pow4(Mst2)) + pow4(Mst1)*(6688*
        pow2(Mst2)*pow2(Mt)*pow2(s2t) + 296*pow4(Mt) - 807*pow4(Mst2)*pow4(s2t)
        ) + 10*pow4(Mst2)*(-160*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 72*pow4(Mt) +
        63*pow4(Mst2)*pow4(s2t)) + 16*(1724*s2t*pow3(Mt) - 21*Mt*pow2(Mst2)*
        pow3(s2t))*pow5(Mst1) + 128*Mst1*(17*s2t*pow3(Mt)*pow4(Mst2) + 33*Mt*
        pow3(s2t)*pow6(Mst2)) + 2*pow2(Mst1)*(3320*pow2(Mt)*pow2(s2t)*pow4(
        Mst2) + 2168*pow2(Mst2)*pow4(Mt) + 97*pow4(s2t)*pow6(Mst2)))) - 20*
        Dmglst1*Mst1*(-120*Mst1*pow3(Mt)*pow4(Mst2)*(3*Mst1*Mt*pow2(Mst2) + 6*
        s2t*pow2(Mst1)*pow2(Mst2) + 5*Mt*pow3(Mst1) + 15*s2t*pow4(Mst1) + s2t*
        pow4(Mst2)) - 240*Mst1*pow2(Msq)*pow2(Mst2)*pow3(Mt)*(3*Mst1*Mt*pow2(
        Mst2) + 12*s2t*pow2(Mst1)*pow2(Mst2) + 30*s2t*pow4(Mst1) + 2*s2t*pow4(
        Mst2)) + pow4(Msq)*(16*pow3(Mst1)*(2209*s2t*pow2(Mst2)*pow3(Mt) + 87*
        Mt*pow3(s2t)*pow4(Mst2)) + pow4(Mst1)*(6312*pow2(Mst2)*pow2(Mt)*pow2(
        s2t) - 3200*pow4(Mt) - 651*pow4(Mst2)*pow4(s2t)) + 16*pow4(Mst2)*(-24*
        pow2(Mst2)*pow2(Mt)*pow2(s2t) + 11*pow4(Mt) + 6*pow4(Mst2)*pow4(s2t)) +
        4*(25312*s2t*pow3(Mt) - 531*Mt*pow2(Mst2)*pow3(s2t))*pow5(Mst1) + 16*
        Mst1*(103*s2t*pow3(Mt)*pow4(Mst2) + 198*Mt*pow3(s2t)*pow6(Mst2)) +
        pow2(Mst1)*(7368*pow2(Mt)*pow2(s2t)*pow4(Mst2) + 6696*pow2(Mst2)*pow4(
        Mt) + 783*pow4(s2t)*pow6(Mst2)))) - 2*pow2(Dmglst1)*(-3600*pow2(Msq)*
        pow2(Mst1)*pow2(Mst2)*((Mt + 8*Mst1*s2t)*pow2(Mst2) + 40*s2t*pow3(Mst1)
        )*pow3(Mt) - 1800*(Mt + 4*Mst1*s2t)*pow2(Mst1)*(5*pow2(Mst1) + pow2(
        Mst2))*pow3(Mt)*pow4(Mst2) + pow4(Msq)*(720*Mst1*s2t*pow3(Mt)*pow4(
        Mst2) - 3920*pow3(Mst1)*(-91*s2t*pow2(Mst2)*pow3(Mt) + 12*Mt*pow3(s2t)*
        pow4(Mst2)) - 5*pow4(Mst1)*(21624*pow2(Mst2)*pow2(Mt)*pow2(s2t) +
        17224*pow4(Mt) - 159*pow4(Mst2)*pow4(s2t)) + 80*(24875*s2t*pow3(Mt) -
        693*Mt*pow2(Mst2)*pow3(s2t))*pow5(Mst1) + 5*pow2(Mst1)*(5832*pow2(Mt)*
        pow2(s2t)*pow4(Mst2) + 6824*pow2(Mst2)*pow4(Mt) + 783*pow4(s2t)*pow6(
        Mst2)) + 24*(-3*pow4(Mst2)*pow4(Mt) - 80*pow2(Mt)*pow2(s2t)*pow6(Mst2)
        + 20*pow4(s2t)*pow8(Mst2)))))))/(45.*Sbeta*pow2(Mst1)*pow4(Msq)*pow4(
        Mst2)*pow4(Mt));

   return result;
}

/**
 * @return returns the susy log^2 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double H5::calc_coef_at_as2_no_sm_logs_log2() const {

   const double result =
      (263040*s2t*Sbeta*pow2(Dmglst1)*pow2(Mst2)*pow3(Mst1)*pow3(Mt) - 22440*
        Sbeta*pow2(Dmglst1)*pow2(Mst2)*pow2(Mt)*pow2(s2t)*pow4(Mst1) + 136320*
        Dmglst1*s2t*Sbeta*pow2(Mst2)*pow3(Mt)*pow4(Mst1) + 33360*Sbeta*pow2(
        Dmglst1)*pow2(Mst1)*pow2(Mt)*pow2(s2t)*pow4(Mst2) + 82080*Dmglst1*
        Sbeta*pow2(Mt)*pow2(s2t)*pow3(Mst1)*pow4(Mst2) - 12160*Mst1*s2t*Sbeta*
        pow2(Dmglst1)*pow3(Mt)*pow4(Mst2) + 36160*Dmglst1*s2t*Sbeta*pow2(Mst1)*
        pow3(Mt)*pow4(Mst2) + 25920*s2t*Sbeta*pow3(Mst1)*pow3(Mt)*pow4(Mst2) +
        42960*Sbeta*pow2(Mt)*pow2(s2t)*pow4(Mst1)*pow4(Mst2) - 19680*Dmglst1*
        Mt*Sbeta*pow3(s2t)*pow4(Mst1)*pow4(Mst2) + 84400*Sbeta*pow2(Dmglst1)*
        pow2(Mst1)*pow2(Mst2)*pow4(Mt) + 86880*Dmglst1*Sbeta*pow2(Mst2)*pow3(
        Mst1)*pow4(Mt) - 57440*Sbeta*pow2(Dmglst1)*pow4(Mst1)*pow4(Mt) + 28080*
        Sbeta*pow2(Mst2)*pow4(Mst1)*pow4(Mt) + 29600*Dmglst1*Mst1*Sbeta*pow4(
        Mst2)*pow4(Mt) - 9744*Sbeta*pow2(Dmglst1)*pow4(Mst2)*pow4(Mt) + 29440*
        Sbeta*pow2(Mst1)*pow4(Mst2)*pow4(Mt) - 7200*Sbeta*log(pow2(Mst1)/pow2(
        Msq))*pow2(Mst1)*pow4(Mst2)*pow4(Mt) - 10005*Sbeta*pow2(Dmglst1)*pow4(
        Mst1)*pow4(Mst2)*pow4(s2t) - 29520*Dmglst1*Sbeta*pow2(Mst2)*pow2(Mt)*
        pow2(s2t)*pow5(Mst1) + 416640*s2t*Sbeta*pow2(Dmglst1)*pow3(Mt)*pow5(
        Mst1) + 17280*s2t*Sbeta*pow2(Mst2)*pow3(Mt)*pow5(Mst1) - 11520*Mt*
        Sbeta*pow2(Dmglst1)*pow2(Mst2)*pow3(s2t)*pow5(Mst1) - 19680*Mt*Sbeta*
        pow3(s2t)*pow4(Mst2)*pow5(Mst1) - 18240*Dmglst1*Sbeta*pow4(Mt)*pow5(
        Mst1) - 8490*Dmglst1*Sbeta*pow4(Mst2)*pow4(s2t)*pow5(Mst1) - 12840*
        Sbeta*pow2(Mst2)*pow2(Mt)*pow2(s2t)*pow6(Mst1) + 97920*Dmglst1*s2t*
        Sbeta*pow3(Mt)*pow6(Mst1) + 5760*Cbeta*Mt*MuSUSY*pow2(Dmglst1)*pow3(
        s2t)*pow6(Mst1) - 11520*Dmglst1*Mt*Sbeta*pow2(Mst2)*pow3(s2t)*pow6(
        Mst1) + 480*Sbeta*pow4(Mt)*pow6(Mst1) - 345*Sbeta*pow4(Mst2)*pow4(s2t)*
        pow6(Mst1) - 14160*Dmglst1*Mst1*Sbeta*pow2(Mt)*pow2(s2t)*pow6(Mst2) +
        600*Sbeta*pow2(Dmglst1)*pow2(Mt)*pow2(s2t)*pow6(Mst2) - 9000*Sbeta*
        pow2(Mst1)*pow2(Mt)*pow2(s2t)*pow6(Mst2) + 11520*Mst1*Mt*Sbeta*pow2(
        Dmglst1)*pow3(s2t)*pow6(Mst2) + 31200*Dmglst1*Mt*Sbeta*pow2(Mst1)*pow3(
        s2t)*pow6(Mst2) + 23520*Mt*Sbeta*pow3(Mst1)*pow3(s2t)*pow6(Mst2) +
        4395*Sbeta*pow2(Dmglst1)*pow2(Mst1)*pow4(s2t)*pow6(Mst2) + 1110*
        Dmglst1*Sbeta*pow3(Mst1)*pow4(s2t)*pow6(Mst2) - 5325*Sbeta*pow4(Mst1)*
        pow4(s2t)*pow6(Mst2) + 1920*s2t*Sbeta*pow3(Mt)*pow7(Mst1) + 3840*Cbeta*
        Dmglst1*Mt*MuSUSY*pow3(s2t)*pow7(Mst1) - 3840*Mt*Sbeta*pow2(Mst2)*pow3(
        s2t)*pow7(Mst1) + 1770*Dmglst1*Mst1*Sbeta*pow4(s2t)*pow8(Mst2) - 75*
        Sbeta*pow2(Dmglst1)*pow4(s2t)*pow8(Mst2) + 3585*Sbeta*pow2(Mst1)*pow4(
        s2t)*pow8(Mst2) + 30*Sbeta*log(pow2(Mst1)/pow2(Mst2))*pow2(Mst1)*(8*
        pow3(Mst1)*(176*s2t*pow2(Mst2)*pow3(Mt) + 73*Mt*pow3(s2t)*pow4(Mst2)) -
        408*pow4(Mst2)*pow4(Mt) + 4*pow4(Mst1)*(192*pow2(Mst2)*pow2(Mt)*pow2(
        s2t) + 172*pow4(Mt) - 33*pow4(Mst2)*pow4(s2t)) + 1472*s2t*pow3(Mt)*
        pow5(Mst1) - 328*pow2(Mt)*pow2(s2t)*pow6(Mst2) + 8*Mst1*(-60*s2t*pow3(
        Mt)*pow4(Mst2) + 41*Mt*pow3(s2t)*pow6(Mst2)) + pow2(Mst1)*(712*pow2(Mt)
        *pow2(s2t)*pow4(Mst2) + 512*pow2(Mst2)*pow4(Mt) + 91*pow4(s2t)*pow6(
        Mst2)) + pow2(Dmglst1)*(14720*s2t*pow3(Mst1)*pow3(Mt) + 384*pow2(Mt)*
        pow2(s2t)*pow4(Mst2) + 384*Mst1*(11*s2t*pow2(Mst2)*pow3(Mt) + 2*Mt*
        pow3(s2t)*pow4(Mst2)) + 512*pow2(Mst2)*pow4(Mt) + pow2(Mst1)*(768*pow2(
        Mst2)*pow2(Mt)*pow2(s2t) + 1568*pow4(Mt) - 91*pow4(Mst2)*pow4(s2t)) +
        91*pow4(s2t)*pow6(Mst2)) + 2*Dmglst1*(3680*s2t*pow3(Mt)*pow4(Mst1) + 4*
        Mt*s2t*(-60*pow2(Mt) + 41*pow2(Mst2)*pow2(s2t))*pow4(Mst2) + 4*pow2(
        Mst1)*(528*s2t*pow2(Mst2)*pow3(Mt) + 137*Mt*pow3(s2t)*pow4(Mst2)) +
        pow3(Mst1)*(768*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 864*pow4(Mt) - 91*pow4(
        Mst2)*pow4(s2t)) + Mst1*(384*pow2(Mt)*pow2(s2t)*pow4(Mst2) + 512*pow2(
        Mst2)*pow4(Mt) + 91*pow4(s2t)*pow6(Mst2))) + 41*pow4(s2t)*pow8(Mst2)))/
        (45.*Sbeta*pow2(Mst1)*pow4(Mst2)*pow4(Mt));

   return result;
}

/**
 * @return returns the susy log^3 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double H5::calc_coef_at_as2_no_sm_logs_log3() const {

   const double result =
      -298.6666666666667;

   return result;
}

}        // hierarchies
}        // himalaya
