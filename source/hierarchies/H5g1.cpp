// ====================================================================
// This file is part of Himalaya.
//
// Himalaya is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================

#include "H5g1.hpp"
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
 * @param lmMt a double log((<renormalization scale> / Mt)^2)
 * @param lmMst1 a double log((<renormalization scale> / Mst1)^2)
 * @param lmMst2 a double log((<renormalization scale> / Mst2)^2)
 * @param lmMsq a double log((<renormalization scale> / Msq)^2)
 * @param Mgl a double gluino mass
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
H5g1::H5g1(const ExpansionFlags_t& expansionDepth, double Al4p, double beta, double Dmglst1,
                 double lmMt, double lmMst1, double lmMst2, double lmMsq,
                 double Mgl, double Mt, double Mst1, double Mst2, double Msq, double MuSUSY,
                 double s2t,
                 int mdrFlag, int oneLoopFlag, int twoLoopFlag, int threeLoopFlag){
   // abbrev for tan(beta) and sin(beta)
   Tbeta = tan(beta);
   Sbeta = sin(beta);
   Cbeta = cos(beta);
   this -> Dmglst1 = Dmglst1;
   this -> lmMt = lmMt;
   this -> lmMst1 = lmMst1;
   this -> lmMst2 = lmMst2;
   this -> lmMsq = lmMsq;
   this -> Mgl = Mgl;
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
   xDmglst1 = expansionDepth.at(ExpansionDepth::Dmglst1);
   xMsq = expansionDepth.at(ExpansionDepth::Msq);
}

/**
 * @return The diagonal (1, 1) matrix element of the Higgs mass matrix as a double for the hierarchy 'H5g1'
 */
double H5g1::getS1() const {
   return (pow2(Mt)*pow2(MuSUSY)*(972*oneLoopFlag*pow2(s2t)*(2 - lmMst1 + lmMst2 -
        (2*(lmMst1 - lmMst2)*pow2(Mst1)*(pow2(Mst1) + pow2(Mst2)))/pow4(Mst2))
        + 16*threeLoopFlag*pow2(Al4p)*(486*pow2(s2t)*(83.61265432098766 + 4*B4
        - (4*DN)/9. - (185*lmMsq)/9. + (25*pow2(lmMsq))/3. - (lmMst1*(3781 -
        300*lmMsq + 180*pow2(lmMsq)))/108. + (lmMst2*(14065 - 3498*lmMst1 + 60*
        lmMsq*(-35 + 12*lmMst1) + 180*pow2(lmMsq) - 18*pow2(lmMst1)))/108. + (
        6.361111111111111 - (5*lmMsq)/3.)*pow2(lmMst1) - ((-2193 + 180*lmMsq +
        442*lmMst1)*pow2(lmMst2))/36. + (Dmglst1*(7 - 60*lmMsq*(-5 + 6*lmMst1)
        - 226*lmMst2 + lmMst1*(-2 + 220*lmMst2) + 180*pow2(lmMsq) + 178*pow2(
        lmMst1) + 18*pow2(lmMst2)))/(18.*Mgl) + ((1 - 2*lmMst2)*pow2(Mst2))/(3.
        *pow2(Mst1)) + (11*pow3(lmMst1))/18. + (71*pow3(lmMst2))/6. + pow2(
        Mst1)*((Mgl*(836 + 143*lmMst1 + (-353 + 750*lmMst1)*lmMst2 + 30*lmMsq*(
        7 - 4*lmMst1 + 4*lmMst2) - 315*pow2(lmMst1) - 435*pow2(lmMst2)) + 50*
        Dmglst1*(82 - 93*lmMst1 + 6*lmMsq*(4 + 3*lmMst1 - 3*lmMst2) + 69*lmMst2
        - 9*pow2(lmMst1) + 9*pow2(lmMst2)))/(135.*Mgl*pow2(Msq)) - (
        26.575942386831276 - (76*B4)/9. + (2*DN)/9. + (35*lmMsq)/3. - 5*pow2(
        lmMsq) + lmMst1*(224.2274074074074 - (380*lmMsq)/9. + (40*pow2(lmMsq))/
        3.) - (149*pow2(lmMst1))/90. - (lmMst2*(372457 - 242670*lmMst1 + 1500*
        lmMsq*(-47 + 24*lmMst1) + 18000*pow2(lmMsq) + 17700*pow2(lmMst1)))/
        1350. + ((-17709 + 2400*lmMsq + 6940*lmMst1)*pow2(lmMst2))/90. + (8*
        pow3(lmMst1))/27. - (1736*pow3(lmMst2))/27. - (2*Dmglst1*(586 + 36*B4 +
        30*DN - 495*lmMsq + 135*pow2(lmMsq) - 6*lmMst1*(-73 - 45*lmMsq + 45*
        pow2(lmMsq)) + 3*lmMst2*(229 - 180*lmMsq + 86*lmMst1 + 90*pow2(lmMsq) -
        79*pow2(lmMst1)) + 18*(-43 + 15*lmMsq)*pow2(lmMst1) - 3*(-372 + 90*
        lmMsq + 49*lmMst1)*pow2(lmMst2) - 65*pow3(lmMst1) + 449*pow3(lmMst2)))/
        (27.*Mgl))/pow2(Mst2)) + pow4(Mst1)*((514500*Dmglst1*(485 - 346*lmMst1
        + 12*lmMsq*(1 + 10*lmMst1 - 10*lmMst2) + 334*lmMst2 - 60*pow2(lmMst1) +
        60*pow2(lmMst2)) - Mgl*(-41220947 + 420*lmMsq*(12479 + 4830*lmMst1 -
        5670*lmMst2) + 1081710*lmMst2 - 210*lmMst1*(30109 + 123480*lmMst2) +
        176400*pow2(lmMsq) + 11951100*pow2(lmMst1) + 14156100*pow2(lmMst2)))/(
        1.11132e7*Mgl*pow4(Msq)) - (265.6519022158834 - (76*B4)/9. + (2*DN)/9.
         - (25*lmMsq)/
        2. + lmMst1*(426.37458616780043 - (605*lmMsq)/9. + (40*pow2(lmMsq))/3.)
        - (66.24060846560846 - (40*lmMsq)/3.)*pow2(lmMst1) - lmMst2*(
        411.7079195011338 - (1386139*lmMst1)/3780. + (5*lmMsq*(-121 + 96*
        lmMst1))/9. + (40*pow2(lmMsq))/3. + (298*pow2(lmMst1))/3.) - (
        298.6850529100529 - 40*lmMsq - (2174*lmMst1)/9.)*pow2(lmMst2) + (80*
        pow3(lmMst1))/27. + (Dmglst1*(112.49099794238683 - (8*B4)/3. - (20*DN)/
        9. + lmMst1*(254.24382716049382 - 120*lmMsq + 20*pow2(lmMsq)) + (
        129.57407407407408 - 20*lmMsq)*pow2(lmMst1) - lmMst2*(
        211.57716049382717 - 120*lmMsq - (1439*lmMst1)/27. + 20*pow2(lmMsq) + (
        338*pow2(lmMst1))/9.) + ((-9491 + 1080*lmMsq + 6636*lmMst1)*pow2(
        lmMst2))/54. + (38*pow3(lmMst1))/9. - (806*pow3(lmMst2))/9.))/Mgl - (
        3920*pow3(lmMst2))/27.)/pow4(Mst2)) - ((1.7499706655148832 + lmMsq*(-
        0.990854119425548 + (8*lmMst1)/9. - (52*lmMst2)/63.) + (10511*lmMst2)/
        4410. - (4*lmMst1*(47 + 30*lmMst2))/135. + (5*Dmglst1)/(6.*Mgl) - (2*
        pow2(lmMsq))/63. + (6*pow2(lmMst2))/7.)*pow2(Mst1)*pow2(Mst2) + (
        0.6658120973257028 + lmMsq*(-0.36171579743008314 + lmMst1/9. - lmMst2/
        7.) + (15647*lmMst2)/26460. - (lmMst1*(31 + 15*lmMst2))/135. + pow2(
        lmMsq)/63. + (8*pow2(lmMst2))/63.)*pow4(Mst2))/pow4(Msq) - (100*
        Dmglst1*(55 + lmMst1 + 6*lmMsq*(-5 + 7*lmMst1 - 7*lmMst2) + 29*lmMst2 -
        21*pow2(lmMst1) + 21*pow2(lmMst2))*pow4(Mst1) + Mgl*((1525 + 30*lmMsq*(
        -25 + 44*lmMst1 - 44*lmMst2) + 1928*lmMst2 - 2*lmMst1*(589 + 900*
        lmMst2) + 240*pow2(lmMst1) + 1560*pow2(lmMst2))*pow4(Mst1) + (939 + 30*
        lmMsq*(-12 + 5*lmMst1 - 5*lmMst2) + 760*lmMst2 - 50*lmMst1*(8 + 3*
        lmMst2) + 150*pow2(lmMst2))*pow4(Mst2)))/(135.*Mgl*pow2(Msq)*pow2(Mst2)
        ) + (8*OepS2*(204*pow2(Mst1)*pow2(Mst2) + 370*pow4(Mst1) + 27*pow4(
        Mst2)) - 12*T1ep*(204*pow2(Mst1)*pow2(Mst2) + 370*pow4(Mst1) + 27*pow4(
        Mst2)) - 27*S2*(36*(49 + 34*lmMst1 - 34*lmMst2)*pow2(Mst1)*pow2(Mst2) +
        (4138 + 2220*lmMst1 - 2220*lmMst2)*pow4(Mst1) + 81*(-11 + 2*lmMst1 - 2*
        lmMst2)*pow4(Mst2)) + (2430*xMsq*(1 - 2*(lmMsq + z2))*pow2(Msq)*(pow2(
        Mst1)*pow2(Mst2) + 2*lmMst1*(-1 + shiftst1)*pow2(Mst1)*pow2(Mst2) - 2*
        lmMst2*(-1 + shiftst1)*pow2(Mst1)*pow2(Mst2) + shiftst2*pow2(Mst1)*(2*
        pow2(Mst1) + pow2(Mst2)) - 2*lmMst1*(1 - 2*shiftst1 + shiftst2)*pow4(
        Mst1) + 2*lmMst2*(1 - 2*shiftst1 + shiftst2)*pow4(Mst1) + pow4(Mst2) -
        shiftst1*(2*pow2(Mst1)*pow2(Mst2) + 2*pow4(Mst1) + pow4(Mst2))))/pow2(
        Mst1) + (243*(-1 + 2*lmMst2)*shiftst3*(2*(1 - 2*lmMst1 + 2*lmMst2)*
        pow2(Mst2)*pow4(Mst1) + 2*(1 - lmMst1 + lmMst2)*pow2(Mst1)*pow4(Mst2) +
        (2 - 6*lmMst1 + 6*lmMst2)*pow6(Mst1) + pow6(Mst2)))/pow2(Mst1))/(729.*
        pow4(Mst2))) + 18*z2*((Mt*(4*(12*(55 + 16*lmMst2)*(2*Dmglst1 + Mgl)*Mt
        + (3*Dmglst1*(-687 + 540*lmMsq + 92*lmMst1 - 860*lmMst2) + (2407 + 180*
        lmMsq - 408*lmMst1 + 408*lmMst2)*Mgl)*Mst1*s2t)*pow2(Mst1) + (3*s2t*(
        960*Dmglst1*pow2(Msq)*pow5(Mst1) + pow2(Mst2)*(-160*(3*Dmglst1 + Mgl)*
        pow2(Msq)*pow3(Mst1) + 4*(627 - 60*lmMsq - 64*lmMst1 + 192*lmMst2)*(
        Dmglst1 + Mgl)*Mst1*pow4(Msq) - 60*(5*Dmglst1 + Mgl)*pow5(Mst1))))/
        pow4(Msq)))/(Mgl*pow4(Mst2)) - 3*pow2(s2t)*(339 - 30*lmMsq - 37*lmMst1
        + 149*lmMst2 + (-198*Dmglst1 + (pow2(Mst1)*(-30*(2*Dmglst1 + Mgl) + (4*
        Dmglst1*((-390 + 90*lmMsq + 139*lmMst1 - 315*lmMst2)*pow2(Msq) + 70*
        pow2(Mst1)))/pow2(Mst2) + (Mgl*((311 + 1080*lmMsq + 2196*lmMst1 - 3348*
        lmMst2)*pow2(Msq) + 360*pow2(Mst1)))/(9.*pow2(Mst2))))/pow2(Msq))/Mgl -
        (25*(4*Dmglst1 + Mgl)*pow4(Mst1))/(2.*Mgl*pow4(Msq)) + (9*(
        71.34670781893004 + (40*lmMsq)/3. + (868*lmMst1)/9. - (332*lmMst2)/3. +
        (Dmglst1*(-1613 + 1080*lmMsq + 3684*lmMst1 - 5796*lmMst2))/(27.*Mgl))*
        pow4(Mst1))/pow4(Mst2) + (6*pow2(Mst2) - (6*shiftst3*(2*(1 - 2*lmMst1 +
        2*lmMst2)*pow2(Mst2)*pow4(Mst1) + 2*(1 - lmMst1 + lmMst2)*pow2(Mst1)*
        pow4(Mst2) + (2 - 6*lmMst1 + 6*lmMst2)*pow6(Mst1) + pow6(Mst2)))/pow4(
        Mst2))/pow2(Mst1))) + (432*Dmglst1*Mst1*s2t*(-2*Mgl*Mst1*s2t*(pow2(
        Mst1) + pow2(Mst2)) - 3*Dmglst1*Mst1*s2t*xDmglst1*(pow2(Mst1) + pow2(
        Mst2)) + 16*Dmglst1*Mt*xDmglst1*(18*pow2(Mst1) + pow2(Mst2)) + 8*Mgl*
        Mt*(21*pow2(Mst1) + 2*pow2(Mst2)))*pow2(z2) + 144*Mgl*pow2(Mst1)*pow2(
        Mt)*(2*Dmglst1*(735 - 6*B4 - 3*DN - 78*lmMst1 + 6*pow2(lmMst1) + 6*
        lmMst2*(85 - 6*lmMst1 + pow2(lmMst1)) + 18*(7 + lmMst1)*pow2(lmMst2) +
        24*lmMt*(2 + 3*lmMst2 - lmMst1*(3 + 2*lmMst2) + pow2(lmMst1) + pow2(
        lmMst2)) - 26*pow3(lmMst1) + 2*pow3(lmMst2)) + Mgl*(417 - 6*B4 - 3*DN -
        162*lmMst1 + 6*lmMst2*(87 - 8*lmMst1 + pow2(lmMst1)) + 18*(8 + lmMst1)*
        pow2(lmMst2) + 24*lmMt*pow2(1 - lmMst1 + lmMst2) - 26*pow3(lmMst1) + 2*
        pow3(lmMst2))) + 108*Dmglst1*Mst1*z4*(2*Mgl*(Mt*s2t*(507*pow2(Mst1) -
        241*pow2(Mst2)) + 60*Mst1*pow2(Mt) - 265*Mst1*(pow2(Mst1) + pow2(Mst2))
        *pow2(s2t)) + Dmglst1*xDmglst1*(Mt*s2t*(3474*pow2(Mst1) - 482*pow2(
        Mst2)) + 180*Mst1*pow2(Mt) - 15*pow2(s2t)*(53*Mst1*pow2(Mst2) + 37*
        pow3(Mst1)))) - 4*z4*pow2(Mgl)*(-3*pow2(Mst1)*(540*pow2(Mt) + 13*pow2(
        Mst2)*pow2(s2t)) + 54*Mt*s2t*(241*Mst1*pow2(Mst2) + 313*pow3(Mst1)) +
        pow2(s2t)*(127*pow4(Mst1) - 3267*pow4(Mst2))) + 27*xDmglst1*z3*pow2(
        Dmglst1)*(-8*pow2(Mst1)*(4*(220 + 9*lmMst1 - 9*lmMst2)*pow2(Mt) - 3*(49
        + 69*lmMst1 - 69*lmMst2)*pow2(Mst2)*pow2(s2t)) + Mt*s2t*(8*(559 + 18*
        lmMst1 - 18*lmMst2)*Mst1*pow2(Mst2) + 24*(151 - 918*lmMst1 + 918*
        lmMst2)*pow3(Mst1)) + pow2(s2t)*(4*(6491 + 318*lmMst1 - 318*lmMst2)*
        pow4(Mst1) - 1533*pow4(Mst2))) - 6*s2t*pow2(Mgl)*pow2(z2)*(-1152*Mst1*
        Mt*pow2(Mst2) + 852*s2t*pow2(Mst1)*pow2(Mst2) - 6336*Mt*pow3(Mst1) +
        1018*s2t*pow4(Mst1) + 315*s2t*pow4(Mst2)) + 2*Mgl*z3*(-216*Mst1*Mt*(
        Mst1*(6*Dmglst1*(13 + 2*lmMst1 - 2*lmMst2)*Mt + (-46 + 6*lmMst1 - 6*
        lmMst2)*Mgl*Mt + 3*Dmglst1*(154 + 225*lmMst1 - 225*lmMst2)*Mst1*s2t +
        23*(25 + 9*lmMst1 - 9*lmMst2)*Mgl*Mst1*s2t) + (Dmglst1*(58 - 9*lmMst1 +
        9*lmMst2) + (133 - 9*lmMst1 + 9*lmMst2)*Mgl)*s2t*pow2(Mst2)) + pow2(
        s2t)*(6*(54*Dmglst1*(35 + 46*lmMst1 - 46*lmMst2) + (8783 + 1134*lmMst1
        - 1134*lmMst2)*Mgl)*pow2(Mst1)*pow2(Mst2) + 4*(54*Dmglst1*(764 + 69*
        lmMst1 - 69*lmMst2) + 7*(4829 + 243*lmMst1 - 243*lmMst2)*Mgl)*pow4(
        Mst1) + (-5265*Dmglst1 + 8856*Mgl)*pow4(Mst2))) + (Mgl*Mst1*Mt*s2t*(3*
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
        Mgl*(2*(9*pow2(Mst2)*(13317 - 408*B4 - 12*DN - 4320*lmMsq + 720*pow2(
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
        )/(pow2(Mgl)*pow4(Mst2))) + (Al4p*(1296*Mgl*s2t*twoLoopFlag*(8*Dmglst1*
        Mt*(Mst1*(12 - 2*lmMst1 + 6*lmMst2 - pow2(lmMst1) + pow2(lmMst2))*pow2(
        Mst2) + (10 + 12*lmMst2 + 4*lmMst1*(-2 + 3*lmMst2) - 9*pow2(lmMst1) -
        3*pow2(lmMst2))*pow3(Mst1)) + 8*Mgl*Mt*(Mst1*(8 - 2*lmMst1 + 6*lmMst2 -
        pow2(lmMst1) + pow2(lmMst2))*pow2(Mst2) + (6 + 4*lmMst1*(-3 + lmMst2) +
        16*lmMst2 - 5*pow2(lmMst1) + pow2(lmMst2))*pow3(Mst1)) - 8*Dmglst1*s2t*
        ((-(lmMst2*(1 + lmMst2)) + pow2(lmMst1))*pow2(Mst1)*pow2(Mst2) + (2 +
        3*lmMst1 - 3*lmMst2 + pow2(lmMst1) - pow2(lmMst2))*pow4(Mst1) - (-1 +
        lmMst1)*pow4(Mst2)) + Mgl*s2t*(4*(2 + 11*lmMst2 - 2*lmMst1*(5 + 4*
        lmMst2) + pow2(lmMst1) + 7*pow2(lmMst2))*pow2(Mst1)*pow2(Mst2) + (1 +
        46*lmMst2 - 2*lmMst1*(23 + 32*lmMst2) + 20*pow2(lmMst1) + 44*pow2(
        lmMst2))*pow4(Mst1) + 4*(7 + 7*lmMst2 - 2*lmMst1*(1 + lmMst2) + 2*pow2(
        lmMst2))*pow4(Mst2))) - (144*z2*(36*Mgl*s2t*twoLoopFlag*pow4(Msq)*(-4*(
        Dmglst1 + Mgl)*Mst1*Mt*pow2(Mst2) - 4*Dmglst1*s2t*pow2(Mst1)*(pow2(
        Mst1) + pow2(Mst2)) + 4*(3*Dmglst1 - Mgl)*Mt*pow3(Mst1) + Mgl*s2t*pow4(
        Mst2)) + xDmglst1*pow2(Dmglst1)*(-72*Mst1*s2t*twoLoopFlag*(-18*Mt*pow2(
        Mst1) + 2*Mt*pow2(Mst2) + 3*Mst1*s2t*pow2(Mst2) + 3*s2t*pow3(Mst1))*
        pow4(Msq) + Al4p*threeLoopFlag*(pow4(Msq)*(24*Mst1*Mt*s2t*((3061 -
        1260*lmMsq - 434*lmMst1 + 2354*lmMst2)*pow2(Mst1) + (-627 + 60*lmMsq +
        64*lmMst1 - 192*lmMst2)*pow2(Mst2)) - 36*pow2(Mst1)*(8*(55 + 16*lmMst2)
        *pow2(Mt) + (390 - 90*lmMsq - 139*lmMst1 + 315*lmMst2)*pow2(Mst2)*pow2(
        s2t)) + pow2(s2t)*((2903 + 3240*lmMsq + 10956*lmMst1 - 17292*lmMst2)*
        pow4(Mst1) - 1782*pow4(Mst2))) + s2t*(150*pow2(Mst2)*(36*Mst1*Mt - 5*
        s2t*pow2(Mst2))*pow4(Mst1) - 60*pow2(Msq)*pow2(Mst1)*(-96*Mst1*Mt*pow2(
        Mst2) - 82*s2t*pow2(Mst1)*pow2(Mst2) + 432*Mt*pow3(Mst1) + 9*s2t*pow4(
        Mst2)))))))/pow4(Msq) + xDmglst1*pow2(Dmglst1)*(5184*s2t*twoLoopFlag*(
        2*Mt*(Mst1*(16 - 2*lmMst1 + 6*lmMst2 - pow2(lmMst1) + pow2(lmMst2))*
        pow2(Mst2) + (2 - 2*lmMst2 + 6*lmMst1*(1 + 4*lmMst2) - 15*pow2(lmMst1)
        - 9*pow2(lmMst2))*pow3(Mst1)) + s2t*((4 + 3*lmMst2*(1 + lmMst2) - 3*
        pow2(lmMst1))*pow2(Mst1)*pow2(Mst2) - 3*(5*lmMst1 - lmMst2*(5 + lmMst2)
        + pow2(lmMst1))*pow4(Mst1) + (-4 + 3*lmMst1)*pow4(Mst2))) + (Al4p*
        threeLoopFlag*(768*pow2(Mst1)*pow2(Mt)*(9893 - 54*B4 - 27*DN - 50*
        lmMst1 + 138*pow2(lmMst1) + lmMst2*(4538 - 300*lmMst1 + 54*pow2(lmMst1)
        ) + 54*(19 + 3*lmMst1)*pow2(lmMst2) + 24*lmMt*(29 + 35*lmMst2 - lmMst1*
        (35 + 18*lmMst2) + 9*pow2(lmMst1) + 9*pow2(lmMst2)) - 234*pow3(lmMst1)
        + 18*pow3(lmMst2))*pow4(Msq) - 48*Mst1*Mt*s2t*(2*(3*pow2(Mst2)*(413 +
        408*B4 + 12*DN + 8520*lmMsq - 720*pow2(lmMsq) + 2*lmMst1*(-1573 - 540*
        lmMsq + 180*pow2(lmMsq)) - 18*(51 + 20*lmMsq)*pow2(lmMst1) + lmMst2*(-
        15206 + 2520*lmMsq + 708*lmMst1 - 360*pow2(lmMsq) + 360*pow2(lmMst1)) +
        6*(-797 + 60*lmMsq + 4*lmMst1)*pow2(lmMst2) + 200*pow3(lmMst1) - 584*
        pow3(lmMst2)) + pow2(Mst1)*(59579 + 20088*B4 - 324*DN - 39240*lmMsq -
        2160*pow2(lmMsq) + 24*lmMst1*(-56 + 1935*lmMsq + 135*pow2(lmMsq)) - 18*
        (773 + 1620*lmMsq)*pow2(lmMst1) - 36*lmMst2*(-940 + 1109*lmMst1 - 90*
        lmMsq*(-13 + 16*lmMst1) + 90*pow2(lmMsq) + 661*pow2(lmMst1)) - 90*(-355
        + 252*lmMsq + 286*lmMst1)*pow2(lmMst2) + 33084*pow3(lmMst1) + 16452*
        pow3(lmMst2)))*pow4(Msq) + 90*(928 - 145*lmMst1 + 180*lmMsq*(-2 +
        lmMst1 - lmMst2) + 505*lmMst2 - 90*pow2(lmMst1) + 90*pow2(lmMst2))*
        pow2(Mst2)*pow4(Mst1) + 240*pow2(Msq)*((557 - 126*lmMst1 + 72*lmMsq*(-2
        + lmMst1 - lmMst2) + 270*lmMst2 - 36*pow2(lmMst1) + 36*pow2(lmMst2))*
        pow2(Mst1)*pow2(Mst2) + (203 + 72*lmMsq*(-2 + 3*lmMst1 - 3*lmMst2) +
        366*lmMst2 + 6*lmMst1*(-37 + 90*lmMst2) - 378*pow2(lmMst1) - 162*pow2(
        lmMst2))*pow4(Mst1) - 3*pow4(Mst2)) + 5*(-467 + 12*lmMsq - 12*lmMst2)*
        pow2(Mst1)*pow4(Mst2) + 5*(11 - 12*lmMsq + 12*lmMst2)*pow6(Mst2)) +
        pow2(s2t)*(pow4(Msq)*(864*pow2(Mst1)*pow2(Mst2)*(479 + 36*B4 + 30*DN -
        975*lmMsq + lmMst1*(1018 + 270*lmMsq - 270*pow2(lmMsq)) + 135*pow2(
        lmMsq) + lmMst2*(811 - 540*lmMsq + 258*lmMst1 + 270*pow2(lmMsq) - 237*
        pow2(lmMst1)) + (-646 + 270*lmMsq)*pow2(lmMst1) - 3*(-372 + 90*lmMsq +
        49*lmMst1)*pow2(lmMst2) - 65*pow3(lmMst1) + 449*pow3(lmMst2)) + (130639
        - 24192*B4 - 1728*DN - 984960*lmMsq - 12*lmMst1*(262259 - 181440*lmMsq
        + 19440*pow2(lmMsq)) + 72*(-21689 + 3240*lmMsq)*pow2(lmMst1) + 12*
        lmMst2*(347507 - 181440*lmMsq - 156324*lmMst1 + 19440*pow2(lmMsq) +
        84312*pow2(lmMst1)) - 72*(-45823 + 3240*lmMsq + 27876*lmMst1)*pow2(
        lmMst2) - 240480*pow3(lmMst1) + 1235808*pow3(lmMst2))*pow4(Mst1) + 72*(
        6457 + 60*lmMsq*(67 - 54*lmMst1) - 982*lmMst1 + 30*(-89 + 66*lmMst1)*
        lmMst2 + 1620*pow2(lmMsq) + 834*pow2(lmMst1) + 162*pow2(lmMst2))*pow4(
        Mst2)) - 1440*pow2(Msq)*(2*(335 - 97*lmMst1 + 6*lmMsq*(-25 + 41*lmMst1
        - 41*lmMst2) + 247*lmMst2 - 123*pow2(lmMst1) + 123*pow2(lmMst2))*pow2(
        Mst2)*pow4(Mst1) - 3*(108 - 133*lmMst1 + 2*lmMsq*(32 + 9*lmMst1 - 9*
        lmMst2) + 69*lmMst2 - 9*pow2(lmMst1) + 9*pow2(lmMst2))*pow2(Mst1)*pow4(
        Mst2)) + 180*((3087 - 2306*lmMst1 + 12*lmMsq*(47 + 50*lmMst1 - 50*
        lmMst2) + 1742*lmMst2 - 300*pow2(lmMst1) + 300*pow2(lmMst2))*pow4(Mst1)
        *pow4(Mst2) - 54*pow2(Mst1)*pow6(Mst2)))))/pow4(Msq))))/(pow2(Mgl)*
        pow4(Mst2))))/7776.;
}

/**
 * @return The diagonal (2, 2) matrix element of the Higgs mass matrix as a double for the hierarchy 'H5g1'
 */
double H5g1::getS2() const {
   return -(oneLoopFlag*((4*Mt*MuSUSY*s2t*(4*(lmMst1 - lmMst2)*pow2(Mt) + ((-2 -
        lmMst1 + lmMst2)*pow2(Mst1) + (2 - lmMst1 + lmMst2)*pow2(Mst2))*pow2(
        s2t)))/Tbeta + pow2(Mt)*pow2(s2t)*(8*(lmMst1 - lmMst2)*(pow2(Mst1) -
        pow2(Mst2)) + 4*pow2(MuSUSY)*(2 - lmMst1 + lmMst2 - (2*(lmMst1 -
        lmMst2)*pow2(Mst1)*(pow2(Mst1) + pow2(Mst2)))/pow4(Mst2)) - (4*pow2(
        MuSUSY)*(2 - lmMst1 + lmMst2 - (2*(lmMst1 - lmMst2)*pow2(Mst1)*(pow2(
        Mst1) + pow2(Mst2)))/pow4(Mst2)))/pow2(Sbeta)) + 16*(lmMst1 + lmMst2 -
        2*lmMt)*pow4(Mt) - (-4*pow2(Mst1)*pow2(Mst2) + (2 + lmMst1 - lmMst2)*
        pow4(Mst1) + (2 - lmMst1 + lmMst2)*pow4(Mst2))*pow4(s2t)))/32. +
        threeLoopFlag*pow2(Al4p)*(-((pow2(Mst1)*pow2(Mst2)*(48.366979423868315
         - B4/
        9. - DN/6. - (265*lmMsq)/36. + (35*pow2(lmMsq))/12. + lmMst1*(
        38.55222222222222 - (55*lmMsq)/6. + (5*pow2(lmMsq))/2.) + (
        2.7666666666666666 - (5*lmMsq)/6.)*pow2(lmMst1) - (lmMst2*(3322 + 3000*
        lmMsq*(-1 + lmMst1) - 25870*lmMst1 + 2250*pow2(lmMsq) + 3025*pow2(
        lmMst1)))/900. + ((-3372 + 750*lmMsq + 2365*lmMst1)*pow2(lmMst2))/180.
         + (41*pow3(lmMst1))/
        108. - (1097*pow3(lmMst2))/108. + ((2*(50*Dmglst1*(137 - 92*lmMst1 +
        lmMsq*(-6 + 60*lmMst1 - 60*lmMst2) + 98*lmMst2 - 30*pow2(lmMst1) + 30*
        pow2(lmMst2)) + Mgl*(2068 - 646*lmMst1 + 15*lmMsq*(-23 + 41*lmMst1 -
        41*lmMst2) + 991*lmMst2 - 225*lmMst1*lmMst2 - 195*pow2(lmMst1) + 420*
        pow2(lmMst2)))*pow2(Mst1))/pow2(Msq) - 5*Dmglst1*(1151 + 72*B4 + 60*DN
        - 1890*lmMsq - 270*pow2(lmMsq) - 18*lmMst1*(-49 - 90*lmMsq + 30*pow2(
        lmMsq)) + 6*lmMst2*(342 - 180*lmMsq - 24*lmMst1 + 90*pow2(lmMsq) - 79*
        pow2(lmMst1)) + 6*(-347 + 90*lmMsq)*pow2(lmMst1) - 6*(-363 + 90*lmMsq +
        49*lmMst1)*pow2(lmMst2) - 130*pow3(lmMst1) + 898*pow3(lmMst2)))/(540.*
        Mgl)) + (32.221840780308305 + (10*B4)/9. + DN/18. - (275*lmMsq)/72. +
        lmMst1*(3.2322576530612244 + (65*lmMsq)/18. - (35*pow2(lmMsq))/12.) + (
        5*pow2(lmMsq))/12. + lmMst2*((-5*lmMsq*(8 + 3*lmMst1))/9. + (35*pow2(
        lmMsq))/12. + (5211957 + 20946380*lmMst1 - 38602200*pow2(lmMst1))/
        2.1168e6) - (17.322652116402118 - (15*lmMsq)/4.)*pow2(lmMst1) + (
        8.48290343915344 - (25*lmMsq)/12. + (1793*lmMst1)/72.)*pow2(lmMst2) + (
        95*pow3(lmMst1))/216. + (Dmglst1*(49.72923096707819 + (2*B4)/3. + (5*
        DN)/9. - (45*lmMsq)/2. + lmMst1*(79.81095679012346 - 15*lmMsq - 5*pow2(
        lmMsq)) + (5*pow2(lmMsq))/2. + lmMst2*(10*lmMsq + 5*pow2(lmMsq) + (-
        31507 + 25692*lmMst1 - 23544*pow2(lmMst1))/1296.) + (1.2546296296296295
         + 5*lmMsq)*pow2(lmMst1) - (2.8564814814814814 + 5*lmMsq - (455*lmMst1)
        /18.)*pow2(lmMst2) - (73*pow3(lmMst1))/54. - (311*pow3(lmMst2))/54.))/
        Mgl - (1535*pow3(lmMst2))/216.)*pow4(Mst1) + (2*OepS2*(-150*pow2(Mst1)*
        pow2(Mst2) + 11*pow4(Mst1) - 27*pow4(Mst2)))/729. - ((51450*Dmglst1*(
        40*(82 - 93*lmMst1 + 6*lmMsq*(4 + 3*lmMst1 - 3*lmMst2) + 69*lmMst2 - 9*
        pow2(lmMst1) + 9*pow2(lmMst2))*pow2(Msq)*pow2(Mst1) + 6*(7 + 60*lmMsq*(
        5 - 6*lmMst1) - 2*lmMst1 + (-226 + 220*lmMst1)*lmMst2 + 180*pow2(lmMsq)
        + 178*pow2(lmMst1) + 18*pow2(lmMst2))*pow4(Msq) + 5*(521 - 346*lmMst1 +
        12*lmMsq*(1 + 10*lmMst1 - 10*lmMst2) + 334*lmMst2 - 60*pow2(lmMst1) +
        60*pow2(lmMst2))*pow4(Mst1)) + Mgl*(41160*(2714 + 30*lmMsq*(-17 + 6*
        lmMst1 - 6*lmMst2) + 1167*lmMst2 + 9*lmMst1*(-73 + 50*lmMst2) - 315*
        pow2(lmMst1) - 135*pow2(lmMst2))*pow2(Msq)*pow2(Mst1) + 8575*(53749 +
        2592*B4 - 288*DN - 13320*lmMsq + 5400*pow2(lmMsq) - 6*lmMst1*(3781 -
        300*lmMsq + 180*pow2(lmMsq)) + 6*lmMst2*(14209 - 3498*lmMst1 + 60*
        lmMsq*(-35 + 12*lmMst1) + 180*pow2(lmMsq) - 18*pow2(lmMst1)) - 18*(-229
        + 60*lmMsq)*pow2(lmMst1) - 18*(-2193 + 180*lmMsq + 442*lmMst1)*pow2(
        lmMst2) + 396*pow3(lmMst1) + 7668*pow3(lmMst2))*pow4(Msq) + 3*(12119532
        + 420*lmMsq*(-9224 + 6545*lmMst1 - 5705*lmMst2) + 7553665*lmMst2 + 35*
        lmMst1*(-105131 + 35280*lmMst2) - 176400*pow2(lmMsq) - 1991850*pow2(
        lmMst1) + 580650*pow2(lmMst2))*pow4(Mst1)))*pow4(Mst2))/(2.22264e7*Mgl*
        pow4(Msq)) - (S2*(-18*(197 + 50*lmMst1 - 50*lmMst2)*pow2(Mst1)*pow2(
        Mst2) + (281 + 66*lmMst1 - 66*lmMst2)*pow4(Mst1) + 81*(11 - 2*lmMst1 +
        2*lmMst2)*pow4(Mst2)))/108.)*pow4(s2t)) + (pow4(Mt)*(10501 - 4920*lmMsq
        + 18*lmMst2*(149 + 50*lmMsq - 19*lmMst1*(-2 + lmMst2) + lmMst2 - 30*
        pow2(lmMsq)) + 1980*pow2(lmMsq) - 4*lmMst1*(196 + 270*lmMsq + 135*pow2(
        lmMsq)) + 1914*pow2(lmMst1) - 6*lmMt*(77 + 2*lmMst1*(221 - 57*lmMst2) +
        594*lmMst2 + 90*lmMsq*(7 - 2*(lmMst1 + lmMst2)) - 129*pow2(lmMst1) +
        87*pow2(lmMst2)) + 108*(-((-28 + 10*lmMsq + 29*lmMst1 + 17*lmMst2)*
        pow2(lmMt)) + ((1 - 2*lmMst2)*pow2(Mst2))/pow2(Mst1)) + 360*pow3(lmMsq)
        - 750*pow3(lmMst1) - 204*pow3(lmMst2) + 4968*pow3(lmMt) - (1029000*
        Dmglst1*(2215 + 254*lmMst2 - 2322*lmMt + 264*lmMst2*lmMt + 60*lmMsq*(5
        + 42*lmMst1 + 12*lmMt) + lmMst1*(824 + 60*lmMst2 + 984*lmMt) - 1620*
        pow2(lmMsq) - 3018*pow2(lmMst1) - 162*pow2(lmMst2) - 1872*pow2(lmMt))*
        pow4(Msq)*pow4(Mst2) + 5488*pow2(Msq)*pow2(Mst1)*(pow2(Msq)*pow2(Mst2)*
        (Mgl*(2459588 - 813345*lmMst2 + 151875*lmMsq*(3 + 2*lmMt) - 151875*
        pow2(lmMsq) + 225*(1339 + 30*lmMst2 - 900*lmMt)*pow2(lmMst1) - 370350*
        pow2(lmMst2) + 3375*lmMt*(66 + 303*lmMst2 + 68*pow2(lmMst2)) + 15*
        lmMst1*(17548 - 71775*lmMt - 15*lmMst2*(773 + 120*lmMt) + 450*pow2(
        lmMst2)) - 162000*pow2(lmMt) + 65250*pow3(lmMst1) - 78750*pow3(lmMst2))
        + 125*Dmglst1*(64093 - 12414*lmMst2 + 2430*lmMsq*(3 + 2*lmMt) - 6*
        lmMst1*(988 + 2367*lmMt + lmMst2*(393 + 36*lmMt)) - 2430*pow2(lmMsq) +
        36*(38 + 3*lmMst2 - 93*lmMt)*pow2(lmMst1) - 6354*pow2(lmMst2) + 18*
        lmMt*(1190 + 933*lmMst2 + 198*pow2(lmMst2)) - 2592*pow2(lmMt) + 1080*
        pow3(lmMst1) - 1188*pow3(lmMst2))) + 2250*pow2(Msq)*pow2(MuSUSY)*(2*
        Dmglst1*(735 - 6*B4 - 3*DN - 78*lmMst1 + 6*pow2(lmMst1) + 6*lmMst2*(85
        - 6*lmMst1 + pow2(lmMst1)) + 18*(7 + lmMst1)*pow2(lmMst2) + 24*lmMt*(2
        + 3*lmMst2 - lmMst1*(3 + 2*lmMst2) + pow2(lmMst1) + pow2(lmMst2)) - 26*
        pow3(lmMst1) + 2*pow3(lmMst2)) + Mgl*(417 - 6*B4 - 3*DN - 162*lmMst1 +
        6*lmMst2*(87 - 8*lmMst1 + pow2(lmMst1)) + 18*(8 + lmMst1)*pow2(lmMst2)
        + 24*lmMt*pow2(1 - lmMst1 + lmMst2) - 26*pow3(lmMst1) + 2*pow3(lmMst2))
        ) + 45*(750*Dmglst1*(27 + 14*lmMst1 + lmMst2 + (-3 - 6*(lmMst1 +
        lmMst2))*lmMt + 6*lmMsq*(-2 + lmMst1 + lmMst2 + 2*lmMt) - 12*pow2(
        lmMsq)) + Mgl*(12339 + 4910*lmMst1 + 375*lmMst2 - 375*(11 + 10*lmMst1 +
        6*lmMst2)*lmMt + 10*lmMsq*(-116 + 285*lmMst1 + 225*lmMst2 + 600*lmMt) -
        5550*pow2(lmMsq) + 450*pow2(lmMst1)))*pow4(Mst2)) + 82320*pow2(Msq)*
        pow2(Mst2)*(-1875*(-11 + 6*lmMsq - 6*lmMt)*(4*Dmglst1 + Mgl)*pow4(Mst1)
        + Mgl*(13517 + lmMst2*(7230 - 4500*lmMt) - 6000*lmMt + 30*lmMsq*(-41 +
        60*lmMst2 + 150*lmMt) - 3150*pow2(lmMsq) + 1350*pow2(lmMst2))*pow4(
        Mst2)) + pow4(Mst1)*(514500*Dmglst1*(2*(99985 + 34432*lmMst2 + (3156 -
        5832*lmMst2 + 1296*lmMt)*pow2(lmMst1) + 17388*pow2(lmMst2) + 24*lmMt*(-
        191 + 871*lmMst2 + 6*pow2(lmMst2)) - 1080*lmMsq*(9 + lmMst2*(29 - 10*
        lmMt) - 10*lmMst1*(1 + lmMst2 - lmMt) - 19*lmMt + 10*pow2(lmMst2)) -
        1728*(3 + 2*lmMst2)*pow2(lmMt) + 8*lmMst1*(517 + lmMst2*(168 - 180*
        lmMt) - 4557*lmMt + 27*pow2(lmMst2) + 432*pow2(lmMt)) + 1512*pow3(
        lmMst1) + 4104*pow3(lmMst2))*pow4(Msq) + 5*(2095 + 330*lmMst2 + 48*
        lmMst1*(25 - 12*lmMt) - 6*(-19 + 60*lmMst2)*lmMt + 12*lmMsq*(-137 + 48*
        lmMst1 + 30*lmMst2 + 78*lmMt) - 936*pow2(lmMsq))*pow4(Mst2)) + Mgl*(4*(
        6336529126 + 2073605310*lmMst2 - 44100*(-8128 + 8295*lmMst2 + 1680*
        lmMt)*pow2(lmMst1) + 976616550*pow2(lmMst2) - 138915000*lmMsq*(2 +
        lmMst2*(12 - 5*lmMt) - 5*lmMst1*(1 + lmMst2 - lmMt) - 7*lmMt + 5*pow2(
        lmMst2)) + 1157625*lmMt*(-899 + 1578*lmMst2 + 128*pow2(lmMst2)) -
        222264000*(1 + lmMst2)*pow2(lmMt) + 315*lmMst1*(1895351 - 7915950*lmMt
        - 70*lmMst2*(16867 + 3360*lmMt) - 14700*pow2(lmMst2) + 705600*pow2(
        lmMt)) + 146632500*pow3(lmMst1) + 223807500*pow3(lmMst2))*pow4(Msq) +
        15*(160842737 + 14148750*lmMst2 + 840*lmMst1*(106067 - 51450*lmMt) -
        257250*(41 + 60*lmMst2)*lmMt + 420*lmMsq*(-220709 + 82740*lmMst1 +
        36750*lmMst2 + 139650*lmMt) - 54419400*pow2(lmMsq) + 4233600*pow2(
        lmMst1))*pow4(Mst2))) + 45*(343000*Dmglst1*(115 - 48*lmMsq + lmMst2*(69
        - 36*lmMt) - 21*lmMt + 36*lmMsq*(lmMst2 + lmMt) - 36*pow2(lmMsq))*pow2(
        Mst1) + Mgl*(171500*(115 - 48*lmMsq + lmMst2*(69 - 36*lmMt) - 21*lmMt +
        36*lmMsq*(lmMst2 + lmMt) - 36*pow2(lmMsq))*pow2(Mst1) + (13597579 +
        210*lmMst2*(43831 - 14700*lmMt) - 2315250*lmMt + 420*lmMsq*(-16403 +
        630*lmMst2 + 7350*lmMt) - 1675800*pow2(lmMsq) + 1411200*pow2(lmMst2))*
        pow2(Mst2)))*pow6(Mst2))/(514500.*Mgl*pow4(Msq)*pow4(Mst2))))/81. + (
        Mst1*Mt*s2t*(-(pow2(s2t)*(8*pow2(Mst1)*(132407 + 11880*B4 - 108*DN -
        45360*lmMsq + 6480*pow2(lmMsq) + 60*lmMst1*(2299 - 378*lmMsq + 54*pow2(
        lmMsq)) + 18*(3377 - 900*lmMsq)*pow2(lmMst1) - 12*lmMst2*(1667 + 1761*
        lmMst1 - 270*lmMsq*(3 + 8*lmMst1) + 270*pow2(lmMsq) + 594*pow2(lmMst1))
        - 18*(859 + 540*lmMsq + 564*lmMst1)*pow2(lmMst2) + 14472*pow3(lmMst1) +
        2808*pow3(lmMst2) + (3*Dmglst1*(54053 + 9432*B4 - 180*DN - 34560*lmMsq
        + 2160*pow2(lmMsq) + 12*(lmMst1*(2489 + 810*lmMsq + 90*pow2(lmMsq)) +
        lmMst2*(2899 - 90*lmMsq*(13 + lmMsq - 24*lmMst1) - 2559*lmMst1 - 984*
        pow2(lmMst1))) + 54*(381 - 260*lmMsq)*pow2(lmMst1) - 18*(-1011 + 660*
        lmMsq + 688*lmMst1)*pow2(lmMst2) + 15648*pow3(lmMst1) + 8544*pow3(
        lmMst2)))/Mgl) + (pow2(Msq)*(-2880*Mgl*(83 + lmMst1*(52 - 36*lmMst2) -
        12*lmMsq*(2 + lmMst1 - lmMst2) - 28*lmMst2 + 24*pow2(lmMst1) + 12*pow2(
        lmMst2)) - 8640*Dmglst1*(107 + lmMst1*(32 - 60*lmMst2) - 12*lmMsq*(2 +
        lmMst1 - lmMst2) - 8*lmMst2 + 36*pow2(lmMst1) + 24*pow2(lmMst2)) + (
        pow2(Msq)*(Mgl*(622151 + 12960*lmMsq*(11 + 10*lmMst1 - 10*lmMst2) -
        135756*lmMst2 + 72*(-4499 + 3060*lmMst2)*pow2(lmMst1) - 185688*pow2(
        lmMst2) - 36*lmMst1*(69 - 14924*lmMst2 + 3048*pow2(lmMst2)) - 110304*
        pow3(lmMst1) - 288*pow3(lmMst2)) + Dmglst1*(3058187 + 12960*lmMsq*(43 +
        66*lmMst1 - 66*lmMst2) - 960732*lmMst2 + 72*(-9407 + 6948*lmMst2)*pow2(
        lmMst1) + 36*lmMst1*(4415 + 9980*lmMst2 - 10824*pow2(lmMst2)) + 400968*
        pow2(lmMst2) - 203616*pow3(lmMst1) + 93024*pow3(lmMst2))))/pow2(Mst2))*
        pow4(Mst1) - 36*pow2(Mst2)*(Mgl*(-160*(37 - 8*lmMst1 + 6*lmMsq*(-2 +
        lmMst1 - lmMst2) + 20*lmMst2 - 3*pow2(lmMst1) + 3*pow2(lmMst2))*pow2(
        Msq)*pow2(Mst1) + 2*(13317 - 408*B4 - 12*DN - 4320*lmMsq + 720*pow2(
        lmMsq) - 6*lmMst1*(79 - 180*lmMsq + 60*pow2(lmMsq)) + 6*lmMst2*(2135 -
        420*lmMsq - 118*lmMst1 + 60*pow2(lmMsq) - 60*pow2(lmMst1)) + 6*(-167 +
        60*lmMsq)*pow2(lmMst1) - 6*(-797 + 60*lmMsq + 4*lmMst1)*pow2(lmMst2) -
        200*pow3(lmMst1) + 584*pow3(lmMst2))*pow4(Msq) - 5*(343 - 18*lmMst1 +
        12*lmMsq*(-13 + 6*lmMst1 - 6*lmMst2) + 174*lmMst2 - 36*pow2(lmMst1) +
        36*pow2(lmMst2))*pow4(Mst1)) + Dmglst1*(-480*(42 - 10*lmMst1 + 6*lmMsq*
        (-2 + lmMst1 - lmMst2) + 22*lmMst2 - 3*pow2(lmMst1) + 3*pow2(lmMst2))*
        pow2(Msq)*pow2(Mst1) + 2*(11337 - 408*B4 - 12*DN - 6480*lmMsq + 720*
        pow2(lmMsq) - 30*lmMst1*(-73 - 36*lmMsq + 12*pow2(lmMsq)) + 6*lmMst2*(
        2331 - 420*lmMsq - 118*lmMst1 + 60*pow2(lmMsq) - 60*pow2(lmMst1)) + 18*
        (-13 + 20*lmMsq)*pow2(lmMst1) - 6*(-797 + 60*lmMsq + 4*lmMst1)*pow2(
        lmMst2) - 200*pow3(lmMst1) + 584*pow3(lmMst2))*pow4(Msq) - 5*(1807 -
        234*lmMst1 + 12*lmMsq*(-61 + 30*lmMst1 - 30*lmMst2) + 966*lmMst2 - 180*
        pow2(lmMst1) + 180*pow2(lmMst2))*pow4(Mst1))) - 180*(48*(Dmglst1 + Mgl)
        *pow2(Msq) + (Dmglst1*(91 - 12*lmMsq + 12*lmMst2) + (43 - 12*lmMsq +
        12*lmMst2)*Mgl)*pow2(Mst1))*pow4(Mst2))/(Mgl*pow4(Msq)))) + (8*pow2(Mt)
        *(24*(Dmglst1*(7007 - 290*lmMst2 + 30*lmMsq*(251 + 102*lmMt + 36*
        lmMst2*(1 + lmMt) - 6*lmMst1*(23 + 6*lmMt)) + 540*(lmMst1 - lmMst2)*
        pow2(lmMsq) + 12*(331 + 15*lmMst2 + 42*lmMt)*pow2(lmMst1) - 288*pow2(
        lmMst2) - 6*lmMt*(817 + 464*lmMst2 + 108*pow2(lmMst2)) - 864*(2 +
        lmMst2)*pow2(lmMt) + 4*lmMst1*(-2438 - 408*lmMt + 3*lmMst2*(61 + 12*
        lmMt) + 9*pow2(lmMst2) + 216*pow2(lmMt)) - 132*pow3(lmMst1) - 84*pow3(
        lmMst2)) + 3*Mgl*(371 + 156*lmMst2 + 90*lmMsq*(11 + 6*lmMt + 4*lmMst2*(
        1 + lmMt) - 2*lmMst1*(5 + 2*lmMt)) + 180*(lmMst1 - lmMst2)*pow2(lmMsq)
        + 12*(60 + 5*lmMst2 + 14*lmMt)*pow2(lmMst1) - 96*pow2(lmMst2) - 6*lmMt*
        (217 + 154*lmMst2 + 36*pow2(lmMst2)) - 288*(1 + lmMst2)*pow2(lmMt) + 6*
        lmMst1*(-145 + 26*lmMt + 8*lmMst2*(5 + lmMt) + 2*pow2(lmMst2) + 48*
        pow2(lmMt)) - 44*pow3(lmMst1) - 28*pow3(lmMst2))) + pow2(Mst1)*((-480*(
        3*Dmglst1*(-(lmMsq*(7 + 42*lmMst1 - 36*lmMst2 + 6*lmMt)) + lmMst1*(103
        + 42*lmMt) - 4*(32 + 9*lmMst2*lmMt + 12*(lmMst2 + lmMt)) + 6*pow2(
        lmMsq)) + Mgl*(27*lmMst1*(3 + 2*lmMt) - 3*lmMsq*(7 + 18*lmMst1 - 12*
        lmMst2 + 6*lmMt) - 2*(55 + 6*lmMt + 6*lmMst2*(4 + 3*lmMt)) + 18*pow2(
        lmMsq))))/pow2(Msq) + (8*(3*Mgl*(83 + 540*lmMsq*(-1 + 2*lmMst1 - 2*
        lmMst2)*(1 + 2*lmMst2 - 2*lmMt) - 18*(163 + 49*lmMst2 - 116*lmMt)*pow2(
        lmMst1) + 6282*pow2(lmMst2) - 36*lmMt*(31 + 160*lmMst2 + 22*pow2(
        lmMst2)) + lmMst2*(6591 - 864*pow2(lmMt)) - 3*lmMst1*(1249 - 1536*lmMt
        + 108*lmMst2*(5 + 4*lmMt) + 66*pow2(lmMst2) - 288*pow2(lmMt)) - 594*
        pow3(lmMst1) + 1674*pow3(lmMst2)) + Dmglst1*(24721 + 76707*lmMst2 +
        1620*lmMsq*(-7 + 6*lmMst1 - 6*lmMst2)*(1 + 2*lmMst2 - 2*lmMt) - 18*(838
        + 525*lmMst2 - 1182*lmMt)*pow2(lmMst1) + 58824*pow2(lmMst2) - 36*lmMt*(
        939 + 1045*lmMst2 + 129*pow2(lmMst2)) - 9*lmMst1*(1119 - 1108*lmMt +
        12*lmMst2*(53 + 154*lmMt) + 138*pow2(lmMst2) - 864*pow2(lmMt)) - 2592*(
        2 + 3*lmMst2)*pow2(lmMt) - 4518*pow3(lmMst1) + 15210*pow3(lmMst2))))/
        pow2(Mst2)) - (240*(18*(1 + 2*lmMst2 - 2*lmMt)*(Dmglst1*(11 - 15*lmMst1
        + 15*lmMst2) + (1 - 3*lmMst1 + 3*lmMst2)*Mgl)*pow4(Mst1) + (Dmglst1 +
        Mgl)*(67 + 12*lmMst2*(2 - 3*lmMt) - 66*lmMt + 6*lmMsq*(7 + 6*(lmMst2 +
        lmMt)) - 36*pow2(lmMsq))*pow4(Mst2)))/(pow2(Msq)*pow2(Mst2)) + (4*pow4(
        Mst1)*((-2*Mgl*(20293 + 9*(3199 + 504*lmMst2 - 852*lmMt)*pow2(lmMst1) -
        1620*lmMsq*(7 + 3*lmMst1*(5 + 2*lmMst2 - 2*lmMt) + 2*lmMt + lmMst2*(-17
        + 6*lmMt) - 6*pow2(lmMst2)) - 41841*pow2(lmMst2) + 54*lmMt*(83 + 412*
        lmMst2 + 82*pow2(lmMst2)) + 18*lmMst1*(4094 + 533*lmMst2 - 1044*lmMt +
        180*lmMst2*lmMt + 42*pow2(lmMst2) - 144*pow2(lmMt)) + 18*lmMst2*(-2903
        + 144*pow2(lmMt)) + 1620*pow3(lmMst1) - 6912*pow3(lmMst2)) + 2*Dmglst1*
        (54037 + 273408*lmMst2 - 9*(14497 + 3276*lmMst2 - 4920*lmMt)*pow2(
        lmMst1) + 1620*lmMsq*(5 + 15*lmMst1*(5 + 2*lmMst2 - 2*lmMt) + 22*lmMt +
        lmMst2*(-97 + 30*lmMt) - 30*pow2(lmMst2)) + 207819*pow2(lmMst2) - 18*
        lmMt*(3201 + 5066*lmMst2 + 900*pow2(lmMst2)) - 2592*(2 + 5*lmMst2)*
        pow2(lmMt) + 6*lmMst1*(-41437 + 8286*lmMt - 15*lmMst2*(341 + 312*lmMt)
        + 342*pow2(lmMst2) + 2160*pow2(lmMt)) - 5508*pow3(lmMst1) + 32940*pow3(
        lmMst2)))*pow4(Msq) + 45*(Dmglst1*(736 + 105*lmMst2 + 3*(53 + 60*
        lmMst2)*lmMt + 4*lmMsq*(-10 + 66*lmMst1 - 45*lmMst2 + 21*lmMt) - 8*
        lmMst1*(28 + 33*lmMt) - 84*pow2(lmMsq)) + Mgl*(200 + 3*lmMt - 72*
        lmMst1*lmMt + 12*lmMsq*(-2 + 6*lmMst1 - 3*lmMst2 + 3*lmMt) + 3*lmMst2*(
        7 + 12*lmMt) - 36*pow2(lmMsq)))*pow4(Mst2)) - 60*(2*(3*Dmglst1 + Mgl)*(
        112 - 12*lmMsq + lmMst2*(51 - 36*lmMt) - 39*lmMt + 36*lmMsq*(lmMst2 +
        lmMt) - 36*pow2(lmMsq))*pow2(Mst1) + (Dmglst1 + Mgl)*(169 - 48*lmMsq +
        lmMst2*(87 - 36*lmMt) - 39*lmMt + 36*lmMsq*(lmMst2 + lmMt) - 36*pow2(
        lmMsq))*pow2(Mst2))*pow6(Mst2) - pow2(MuSUSY)*(3*Dmglst1*(2*(pow2(Mst1)
        *(13969 - 11880*B4 + 108*DN - 4320*lmMsq + (-16728 - 3240*lmMsq*(1 +
        lmMsq))*lmMst1 + 2160*pow2(lmMsq) + 54*(-407 + 300*lmMsq)*pow2(lmMst1)
        + 12*lmMst2*(4094 + 2205*lmMst1 - 90*lmMsq*(1 + 24*lmMst1) + 270*pow2(
        lmMsq) + 804*pow2(lmMst1)) + 18*(583 + 540*lmMsq + 680*lmMst1)*pow2(
        lmMst2) - 16848*pow3(lmMst1) - 5040*pow3(lmMst2)) + 3*pow2(Mst2)*(11337
        - 408*B4 - 12*DN - 6480*lmMsq + 720*pow2(lmMsq) - 30*lmMst1*(-73 - 36*
        lmMsq + 12*pow2(lmMsq)) + 6*lmMst2*(2331 - 420*lmMsq - 118*lmMst1 + 60*
        pow2(lmMsq) - 60*pow2(lmMst1)) + 18*(-13 + 20*lmMsq)*pow2(lmMst1) - 6*(
        -797 + 60*lmMsq + 4*lmMst1)*pow2(lmMst2) - 200*pow3(lmMst1) + 584*pow3(
        lmMst2)))*pow4(Msq) - 720*pow2(Msq)*(2*(41 - 10*lmMst1 + 6*lmMsq*(-2 +
        lmMst1 - lmMst2) + 22*lmMst2 - 3*pow2(lmMst1) + 3*pow2(lmMst2))*pow2(
        Mst1)*pow2(Mst2) + 2*(29 + 48*lmMst2 - 24*pow2(lmMst1) - 6*(lmMst1*(6 -
        5*lmMst2) + lmMsq*(2 - 3*lmMst1 + 3*lmMst2) + pow2(lmMst2)))*pow4(Mst1)
        - pow4(Mst2)) - 5*pow2(Mst2)*((-251 + 12*lmMsq - 12*lmMst2)*pow2(Mst1)*
        pow2(Mst2) + 6*(818 - 117*lmMst1 + 180*lmMsq*(-2 + lmMst1 - lmMst2) +
        477*lmMst2 - 90*pow2(lmMst1) + 90*pow2(lmMst2))*pow4(Mst1) + (11 - 12*
        lmMsq + 12*lmMst2)*pow4(Mst2))) + Mgl*(2*(9*pow2(Mst2)*(13317 - 408*B4
        - 12*DN - 4320*lmMsq + 720*pow2(lmMsq) - 6*lmMst1*(79 - 180*lmMsq + 60*
        pow2(lmMsq)) + 6*lmMst2*(2135 - 420*lmMsq - 118*lmMst1 + 60*pow2(lmMsq)
        - 60*pow2(lmMst1)) + 6*(-167 + 60*lmMsq)*pow2(lmMst1) - 6*(-797 + 60*
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
        pow6(Mst2)))))/(pow4(Msq)*pow4(Mst2))))/Mgl))/3888. + pow2(Mt)*(pow2(
        s2t)*(pow2(Mst1)*(268.4533251028807 - (8*DN)/9. - (100*lmMsq)/3. + 20*
        pow2(lmMsq) - (2*lmMst1*(4703 + 500*lmMsq + 375*pow2(lmMsq)))/225. + (
        2*lmMst2*(33703 + 1745*lmMst1 - 250*lmMsq*(16 + 3*lmMst1) + 375*pow2(
        lmMsq) - 1025*pow2(lmMst1)))/225. + ((-273 + 100*lmMsq)*pow2(lmMst1))/
        15. + ((4601 + 1770*lmMst1)*pow2(lmMst2))/45. + (32*lmMt*pow2(1 -
        lmMst1 + lmMst2))/3. - (64*pow3(lmMst1))/3. + (Dmglst1*(156781 - 864*B4
        - 432*DN - 10260*lmMsq - 30*(95 + 324*lmMsq)*lmMst1 + 9720*pow2(lmMsq)
        + 7524*pow2(lmMst1) + 6*lmMst2*(21991 - 1620*lmMsq + 372*lmMst1 + 432*
        pow2(lmMst1)) + 36*(1025 + 168*lmMst1)*pow2(lmMst2) + 5184*lmMt*(2 + 3*
        lmMst2 - lmMst1*(3 + 2*lmMst2) + pow2(lmMst1) + pow2(lmMst2)) - 7200*
        pow3(lmMst1) - 1440*pow3(lmMst2)))/(243.*Mgl) - (80*pow3(lmMst2))/9.) -
        ((2*(750*Dmglst1*(131 - 120*lmMsq + 78*lmMst1 + 42*lmMst2) + Mgl*(22642
        + 15*lmMsq*(-1282 + 195*lmMst1 - 75*lmMst2) + 7500*lmMst2 + 15*lmMst1*(
        782 + 75*lmMst2) - 900*pow2(lmMsq) - 2025*pow2(lmMst1))))/(2025.*Mgl*
        pow2(Msq)) + (238.25271580142933 - (50*lmMsq)/3. + (543.1291609977325 -
        (100*lmMsq)/3.)*lmMst1 + 10*pow2(lmMsq) + (266299*pow2(lmMst1))/945. +
        (lmMst2*(-14193447 + 441000*lmMsq - 9820930*lmMst1 + 455700*pow2(
        lmMst1)))/33075. + (48.90899470899471 - (700*lmMst1)/9.)*pow2(lmMst2) -
        (64*lmMt*(-lmMst1 + lmMst2 - 2*lmMst1*lmMst2 + pow2(lmMst1) + pow2(
        lmMst2)))/3. + (836*pow3(lmMst1))/27. + (892*pow3(lmMst2))/27. + (
        Dmglst1*(63295 - 222366*lmMst2 + 1620*lmMsq*(9 - 20*lmMst1 + 14*lmMst2)
        + 4860*pow2(lmMsq) - 144*(-1355 + 39*lmMst2 + 144*lmMt)*pow2(lmMst1) +
        6*lmMst1*(42083 + 6048*lmMt + 48*lmMst2*(-953 + 144*lmMt) - 9432*pow2(
        lmMst2)) + 97488*pow2(lmMst2) - 5184*lmMt*(2 + 7*lmMst2 + 4*pow2(
        lmMst2)) + 29520*pow3(lmMst1) + 32688*pow3(lmMst2)))/(243.*Mgl))/pow2(
        Mst2))*pow4(Mst1) + (pow2(Mst2)*(51450*Dmglst1*(120*(7 - 20*lmMsq + 26*
        lmMst1 - 6*lmMst2)*pow2(Msq)*pow2(Mst1) + 4*(158 + 60*lmMsq*(-4 + 9*
        lmMst1) + lmMst1*(70 - 138*lmMst2) + 96*lmMst2 - 270*pow2(lmMsq) - 459*
        pow2(lmMst1) - 27*pow2(lmMst2))*pow4(Msq) + 5*(17 - 324*lmMsq + 408*
        lmMst1 - 84*lmMst2)*pow4(Mst1)) - Mgl*(2744*(3091 + 10020*lmMst2 + 30*
        lmMst1*(-491 + 150*lmMst2) + 30*lmMsq*(157 - 60*(lmMst1 + lmMst2)) +
        1800*pow2(lmMsq) - 1350*pow2(lmMst1) - 1350*pow2(lmMst2))*pow2(Msq)*
        pow2(Mst1) + 8575*(4501 - 288*DN - 10896*lmMst1 + 1080*(3 - lmMst1 +
        lmMst2)*pow2(lmMsq) + 2952*pow2(lmMst1) - 24*lmMst2*(-733 + 336*lmMst1
        + 36*pow2(lmMst1)) - 72*(-191 + 44*lmMst1)*pow2(lmMst2) - 360*lmMsq*(3
        + lmMst1*(4 - 6*lmMst2) + 14*lmMst2 + 6*pow2(lmMst2)) - 696*pow3(
        lmMst1) + 4728*pow3(lmMst2))*pow4(Msq) - 60*(130598 + lmMst1*(423822 -
        46305*lmMst2) - 105*lmMsq*(2929 + 231*lmMst1 - 441*lmMst2) - 116277*
        lmMst2 - 11025*pow2(lmMsq) + 35280*pow2(lmMst1))*pow4(Mst1))))/(2.7783e6
        *Mgl*pow4(Msq)) + ((2*(5108 + 15*lmMsq*(-218 + 75*lmMst1 - 195*
        lmMst2) + 6270*lmMst2 - 375*lmMst1*(8 + 3*lmMst2) + 900*pow2(lmMsq) +
        2025*pow2(lmMst2)))/(2025.*pow2(Msq)) + (2*(-1 + 2*lmMst2))/(3.*pow2(
        Mst1)) - ((2.3647986178598424 - (991*lmMsq)/441. + (1.4 - lmMsq)*lmMst1
        + (0.8471655328798186 + (11*lmMsq)/21. + lmMst1)*lmMst2 + (5*Dmglst1*(
        13 - 12*lmMsq + 12*lmMst2))/(18.*Mgl) + (5*pow2(lmMsq))/21. - (16*pow2(
        lmMst2))/21.)*pow2(Mst1))/pow4(Msq))*pow4(Mst2) + (2*((8*OepS2 - 27*(37
        + 6*lmMst1 - 6*lmMst2)*S2)*pow2(Mst2) + (-1480*OepS2 + 27*(2069 + 1110*
        lmMst1 - 1110*lmMst2)*S2)*pow2(MuSUSY))*pow4(Mst1) - 27*((8*OepS2 - 81*
        (1 + 2*lmMst1 - 2*lmMst2)*S2)*pow2(Mst2) + (8*OepS2 + 81*(11 - 2*lmMst1
        + 2*lmMst2)*S2)*pow2(MuSUSY))*pow4(Mst2) - 3*pow2(Mst1)*(4*(136*OepS2 -
        81*(49 + 34*lmMst1 - 34*lmMst2)*S2)*pow2(Mst2)*pow2(MuSUSY) + (40*OepS2
        - 81*(61 + 10*lmMst1 - 10*lmMst2)*S2)*pow4(Mst2)))/(729.*pow4(Mst2)) -
        pow2(MuSUSY)*(83.61265432098766 + 4*B4 - (4*DN)/9. - (185*lmMsq)/9. + (
        25*pow2(lmMsq))/3. - (lmMst1*(3781 - 300*lmMsq + 180*pow2(lmMsq)))/108.
         + (lmMst2*(14065 - 3498*lmMst1 + 60*lmMsq*(-35 + 12*lmMst1) + 180*
        pow2(lmMsq) - 18*pow2(lmMst1)))/108. + (6.361111111111111 - (5*lmMsq)/
        3.)*pow2(lmMst1) - ((-2193 + 180*lmMsq + 442*lmMst1)*pow2(lmMst2))/36.
         + ((Dmglst1*(7 - 60*lmMsq*(-5 + 6*lmMst1) - 226*lmMst2 + lmMst1*(-2 +
        220*lmMst2) + 180*pow2(lmMsq) + 178*pow2(lmMst1) + 18*pow2(lmMst2)))/
        18. + ((Mgl*(836 + 143*lmMst1 + (-353 + 750*lmMst1)*lmMst2 + 30*lmMsq*(
        7 - 4*lmMst1 + 4*lmMst2) - 315*pow2(lmMst1) - 435*pow2(lmMst2)) + 50*
        Dmglst1*(82 - 93*lmMst1 + 6*lmMsq*(4 + 3*lmMst1 - 3*lmMst2) + 69*lmMst2
        - 9*pow2(lmMst1) + 9*pow2(lmMst2)))*pow2(Mst1))/(135.*pow2(Msq)))/Mgl +
        ((1 - 2*lmMst2)*pow2(Mst2))/(3.*pow2(Mst1)) + (11*pow3(lmMst1))/18. + (
        71*pow3(lmMst2))/6. - (pow2(Mst1)*(26.575942386831276 - (76*B4)/9. + (
        2*DN)/9. + (35*lmMsq)/3. - 5*pow2(lmMsq) + lmMst1*(224.2274074074074 -
        (380*lmMsq)/9. + (40*pow2(lmMsq))/3.) - (149*pow2(lmMst1))/90. - (
        lmMst2*(372457 - 242670*lmMst1 + 1500*lmMsq*(-47 + 24*lmMst1) + 18000*
        pow2(lmMsq) + 17700*pow2(lmMst1)))/1350. + ((-17709 + 2400*lmMsq +
        6940*lmMst1)*pow2(lmMst2))/90. + (8*pow3(lmMst1))/27. - (1736*pow3(
        lmMst2))/27. - (2*Dmglst1*(586 + 36*B4 + 30*DN - 495*lmMsq + 135*pow2(
        lmMsq) - 6*lmMst1*(-73 - 45*lmMsq + 45*pow2(lmMsq)) + 3*lmMst2*(229 -
        180*lmMsq + 86*lmMst1 + 90*pow2(lmMsq) - 79*pow2(lmMst1)) + 18*(-43 +
        15*lmMsq)*pow2(lmMst1) - 3*(-372 + 90*lmMsq + 49*lmMst1)*pow2(lmMst2) -
        65*pow3(lmMst1) + 449*pow3(lmMst2)))/(27.*Mgl)))/pow2(Mst2) + pow4(
        Mst1)*((514500*Dmglst1*(485 - 346*lmMst1 + 12*lmMsq*(1 + 10*lmMst1 -
        10*lmMst2) + 334*lmMst2 - 60*pow2(lmMst1) + 60*pow2(lmMst2)) - Mgl*(-
        41220947 + 420*lmMsq*(12479 + 4830*lmMst1 - 5670*lmMst2) + 1081710*
        lmMst2 - 210*lmMst1*(30109 + 123480*lmMst2) + 176400*pow2(lmMsq) +
        11951100*pow2(lmMst1) + 14156100*pow2(lmMst2)))/(1.11132e7*Mgl*pow4(
        Msq)) - (265.6519022158834 - (76*B4)/9. + (2*DN)/9. - (25*lmMsq)/2. +
        lmMst1*(426.37458616780043 - (605*lmMsq)/9. + (40*pow2(lmMsq))/3.) - (
        66.24060846560846 - (40*lmMsq)/3.)*pow2(lmMst1) - lmMst2*(
        411.7079195011338 - (1386139*lmMst1)/3780. + (5*lmMsq*(-121 + 96*
        lmMst1))/9. + (40*pow2(lmMsq))/3. + (298*pow2(lmMst1))/3.) - (
        298.6850529100529 - 40*lmMsq - (2174*lmMst1)/9.)*pow2(lmMst2) + (80*
        pow3(lmMst1))/27. + (Dmglst1*(112.49099794238683 - (8*B4)/3. - (20*DN)/
        9. + lmMst1*(254.24382716049382 - 120*lmMsq + 20*pow2(lmMsq)) + (
        129.57407407407408 - 20*lmMsq)*pow2(lmMst1) - lmMst2*(
        211.57716049382717 - 120*lmMsq - (1439*lmMst1)/27. + 20*pow2(lmMsq) + (
        338*pow2(lmMst1))/9.) + ((-9491 + 1080*lmMsq + 6636*lmMst1)*pow2(
        lmMst2))/54. + (38*pow3(lmMst1))/9. - (806*pow3(lmMst2))/9.))/Mgl - (
        3920*pow3(lmMst2))/27.)/pow4(Mst2)) - ((1.7499706655148832 + lmMsq*(-
        0.990854119425548 + (8*lmMst1)/9. - (52*lmMst2)/63.) + (10511*lmMst2)/
        4410. - (4*lmMst1*(47 + 30*lmMst2))/135. + (5*Dmglst1)/(6.*Mgl) - (2*
        pow2(lmMsq))/63. + (6*pow2(lmMst2))/7.)*pow2(Mst1)*pow2(Mst2) + (
        0.6658120973257028 + lmMsq*(-0.36171579743008314 + lmMst1/9. - lmMst2/
        7.) + (15647*lmMst2)/26460. - (lmMst1*(31 + 15*lmMst2))/135. + pow2(
        lmMsq)/63. + (8*pow2(lmMst2))/63.)*pow4(Mst2))/pow4(Msq) - (100*
        Dmglst1*(55 + lmMst1 + 6*lmMsq*(-5 + 7*lmMst1 - 7*lmMst2) + 29*lmMst2 -
        21*pow2(lmMst1) + 21*pow2(lmMst2))*pow4(Mst1) + Mgl*((1525 + 30*lmMsq*(
        -25 + 44*lmMst1 - 44*lmMst2) + 1928*lmMst2 - 2*lmMst1*(589 + 900*
        lmMst2) + 240*pow2(lmMst1) + 1560*pow2(lmMst2))*pow4(Mst1) + (939 + 30*
        lmMsq*(-12 + 5*lmMst1 - 5*lmMst2) + 760*lmMst2 - 50*lmMst1*(8 + 3*
        lmMst2) + 150*pow2(lmMst2))*pow4(Mst2)))/(135.*Mgl*pow2(Msq)*pow2(Mst2)
        ))) + (pow2(MuSUSY)*(pow2(s2t)*(83.61265432098766 + 4*B4 - (4*DN)/9. -
        (185*lmMsq)/9. + (25*pow2(lmMsq))/3. - (lmMst1*(3781 - 300*lmMsq + 180*
        pow2(lmMsq)))/108. + (lmMst2*(14065 - 3498*lmMst1 + 60*lmMsq*(-35 + 12*
        lmMst1) + 180*pow2(lmMsq) - 18*pow2(lmMst1)))/108. + ((229 - 60*lmMsq)*
        pow2(lmMst1))/36. - ((-2193 + 180*lmMsq + 442*lmMst1)*pow2(lmMst2))/36.
         + ((Dmglst1*(7 - 60*lmMsq*(-5 + 6*lmMst1) - 226*lmMst2 + lmMst1*(-2 +
        220*lmMst2) + 180*pow2(lmMsq) + 178*pow2(lmMst1) + 18*pow2(lmMst2)))/
        18. + ((Mgl*(836 + 143*lmMst1 + (-353 + 750*lmMst1)*lmMst2 + 30*lmMsq*(
        7 - 4*lmMst1 + 4*lmMst2) - 315*pow2(lmMst1) - 435*pow2(lmMst2)) + 50*
        Dmglst1*(82 - 93*lmMst1 + 6*lmMsq*(4 + 3*lmMst1 - 3*lmMst2) + 69*lmMst2
        - 9*pow2(lmMst1) + 9*pow2(lmMst2)))*pow2(Mst1))/(135.*pow2(Msq)))/Mgl +
        ((1 - 2*lmMst2)*pow2(Mst2))/(3.*pow2(Mst1)) + (11*pow3(lmMst1))/18. + (
        71*pow3(lmMst2))/6. - (pow2(Mst1)*(26.575942386831276 - (76*B4)/9. + (
        2*DN)/9. + (35*lmMsq)/3. - 5*pow2(lmMsq) + lmMst1*(224.2274074074074 -
        (380*lmMsq)/9. + (40*pow2(lmMsq))/3.) - (149*pow2(lmMst1))/90. - (
        lmMst2*(372457 - 242670*lmMst1 + 1500*lmMsq*(-47 + 24*lmMst1) + 18000*
        pow2(lmMsq) + 17700*pow2(lmMst1)))/1350. + ((-17709 + 2400*lmMsq +
        6940*lmMst1)*pow2(lmMst2))/90. + (8*pow3(lmMst1))/27. - (1736*pow3(
        lmMst2))/27. - (2*Dmglst1*(586 + 36*B4 + 30*DN - 495*lmMsq + 135*pow2(
        lmMsq) - 6*lmMst1*(-73 - 45*lmMsq + 45*pow2(lmMsq)) + 3*lmMst2*(229 -
        180*lmMsq + 86*lmMst1 + 90*pow2(lmMsq) - 79*pow2(lmMst1)) + 18*(-43 +
        15*lmMsq)*pow2(lmMst1) - 3*(-372 + 90*lmMsq + 49*lmMst1)*pow2(lmMst2) -
        65*pow3(lmMst1) + 449*pow3(lmMst2)))/(27.*Mgl)))/pow2(Mst2) + pow4(
        Mst1)*((514500*Dmglst1*(485 - 346*lmMst1 + 12*lmMsq*(1 + 10*lmMst1 -
        10*lmMst2) + 334*lmMst2 - 60*pow2(lmMst1) + 60*pow2(lmMst2)) - Mgl*(-
        41220947 + 420*lmMsq*(12479 + 4830*lmMst1 - 5670*lmMst2) + 1081710*
        lmMst2 - 210*lmMst1*(30109 + 123480*lmMst2) + 176400*pow2(lmMsq) +
        11951100*pow2(lmMst1) + 14156100*pow2(lmMst2)))/(1.11132e7*Mgl*pow4(
        Msq)) - (265.6519022158834 - (76*B4)/9. + (2*DN)/9. - (25*lmMsq)/2. +
        lmMst1*(426.37458616780043 - (605*lmMsq)/9. + (40*pow2(lmMsq))/3.) - (
        66.24060846560846 - (40*lmMsq)/3.)*pow2(lmMst1) - lmMst2*(
        411.7079195011338 - (1386139*lmMst1)/3780. + (5*lmMsq*(-121 + 96*
        lmMst1))/9. + (40*pow2(lmMsq))/3. + (298*pow2(lmMst1))/3.) - (
        298.6850529100529 - 40*lmMsq - (2174*lmMst1)/9.)*pow2(lmMst2) + (80*
        pow3(lmMst1))/27. + (Dmglst1*(112.49099794238683 - (8*B4)/3. - (20*DN)/
        9. + lmMst1*(254.24382716049382 - 120*lmMsq + 20*pow2(lmMsq)) + (
        129.57407407407408 - 20*lmMsq)*pow2(lmMst1) - lmMst2*(
        211.57716049382717 - 120*lmMsq - (1439*lmMst1)/27. + 20*pow2(lmMsq) + (
        338*pow2(lmMst1))/9.) + ((-9491 + 1080*lmMsq + 6636*lmMst1)*pow2(
        lmMst2))/54. + (38*pow3(lmMst1))/9. - (806*pow3(lmMst2))/9.))/Mgl - (
        3920*pow3(lmMst2))/27.)/pow4(Mst2)) + (S2*(81*(11 - 2*lmMst1 + 2*
        lmMst2) - (36*(49 + 34*lmMst1 - 34*lmMst2)*pow2(Mst1))/pow2(Mst2) - (2*
        (2069 + 1110*lmMst1 - 1110*lmMst2)*pow4(Mst1))/pow4(Mst2)))/27. + (8*
        OepS2*(204*pow2(Mst1)*pow2(Mst2) + 370*pow4(Mst1) + 27*pow4(Mst2)))/(
        729.*pow4(Mst2)) - ((1.7499706655148832 + lmMsq*(-0.990854119425548 + (
        8*lmMst1)/9. - (52*lmMst2)/63.) + (10511*lmMst2)/4410. - (4*lmMst1*(47
        + 30*lmMst2))/135. + (5*Dmglst1)/(6.*Mgl) - (2*pow2(lmMsq))/63. + (6*
        pow2(lmMst2))/7.)*pow2(Mst1)*pow2(Mst2) + (0.6658120973257028 + lmMsq*(
        -0.36171579743008314 + lmMst1/9. - lmMst2/7.) + (15647*lmMst2)/26460. -
        (lmMst1*(31 + 15*lmMst2))/135. + pow2(lmMsq)/63. + (8*pow2(lmMst2))/63.
        )*pow4(Mst2))/pow4(Msq) - (100*Dmglst1*(55 + lmMst1 + 6*lmMsq*(-5 + 7*
        lmMst1 - 7*lmMst2) + 29*lmMst2 - 21*pow2(lmMst1) + 21*pow2(lmMst2))*
        pow4(Mst1) + Mgl*((1525 + 30*lmMsq*(-25 + 44*lmMst1 - 44*lmMst2) +
        1928*lmMst2 - 2*lmMst1*(589 + 900*lmMst2) + 240*pow2(lmMst1) + 1560*
        pow2(lmMst2))*pow4(Mst1) + (939 + 30*lmMsq*(-12 + 5*lmMst1 - 5*lmMst2)
        + 760*lmMst2 - 50*lmMst1*(8 + 3*lmMst2) + 150*pow2(lmMst2))*pow4(Mst2))
        )/(135.*Mgl*pow2(Msq)*pow2(Mst2))) + (Mst1*Mt*(144*Mst1*Mt*(2*Dmglst1*(
        735 - 6*B4 - 3*DN - 78*lmMst1 + 6*pow2(lmMst1) + 6*lmMst2*(85 - 6*
        lmMst1 + pow2(lmMst1)) + 18*(7 + lmMst1)*pow2(lmMst2) + 24*lmMt*(2 + 3*
        lmMst2 - lmMst1*(3 + 2*lmMst2) + pow2(lmMst1) + pow2(lmMst2)) - 26*
        pow3(lmMst1) + 2*pow3(lmMst2)) + Mgl*(417 - 6*B4 - 3*DN - 162*lmMst1 +
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
        lmMst2)*pow4(Mst2))) + Mgl*(2*(9*pow2(Mst2)*(13317 - 408*B4 - 12*DN -
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
        pow6(Mst2)))))/pow4(Msq)))/(486.*Mgl*pow4(Mst2))))/pow2(Sbeta)) + (
        MuSUSY*(s2t*pow3(Mt)*(29.117283950617285 - (16*DN)/9. - (20*lmMsq)/3. +
        20*pow2(lmMsq) - (4*lmMst1*(454 + 60*lmMsq + 45*pow2(lmMsq)))/27. + (4*
        lmMst2*(715 - 336*lmMst1 + 30*lmMsq*(-7 + 3*lmMst1) + 45*pow2(lmMsq) -
        36*pow2(lmMst1)))/27. + (164*pow2(lmMst1))/9. - (4*(-191 + 30*lmMsq +
        44*lmMst1)*pow2(lmMst2))/9. - (116*pow3(lmMst1))/27. + (788*pow3(
        lmMst2))/27. + pow2(Mst2)*((4 - 8*lmMst2)/(3.*pow2(Mst1)) + (1715*(30*
        Dmglst1*(13 - 12*lmMsq + 12*lmMst2) + Mgl*(195 - 60*lmMsq*(3 + 2*lmMst1
        - 2*lmMst2) + 4*lmMst2 + 8*lmMst1*(22 + 15*lmMst2) - 120*pow2(lmMst2)))
        *pow2(Mst1) - Mgl*(103583 + 420*lmMsq*(-256 + 49*lmMst1 - 259*lmMst2) +
        150052*lmMst2 - 1372*lmMst1*(31 + 15*lmMst2) + 44100*pow2(lmMsq) +
        64680*pow2(lmMst2))*pow2(Mst2))/(92610.*Mgl*pow4(Msq))) + S2*(-6*(1 +
        2*lmMst1 - 2*lmMst2) - (28*(5 + 2*lmMst1 - 2*lmMst2)*pow2(Mst1))/(3.*
        pow2(Mst2)) - (8*(139 + 60*lmMst1 - 60*lmMst2)*pow4(Mst1))/(27.*pow4(
        Mst2))) + (16*OepS2*(42*pow2(Mst1)*pow2(Mst2) + 40*pow4(Mst1) + 27*
        pow4(Mst2)))/(729.*pow4(Mst2)) + ((4*Dmglst1*(-158 - 60*lmMsq*(-4 + 9*
        lmMst1) - 96*lmMst2 + 2*lmMst1*(-35 + 69*lmMst2) + 270*pow2(lmMsq) +
        459*pow2(lmMst1) + 27*pow2(lmMst2)))/27. + (2*pow2(Mst1)*((-30*(2250*
        Dmglst1*(7 - 20*lmMsq + 26*lmMst1 - 6*lmMst2) + Mgl*(2017 + 15*lmMst1*(
        782 - 375*lmMst2) + 15*lmMsq*(-532 + 195*lmMst1 - 75*lmMst2) - 3750*
        lmMst2 - 900*pow2(lmMsq) + 1350*pow2(lmMst1) + 3375*pow2(lmMst2))))/
        pow2(Msq) - (Mgl*(7874051 + 7815060*lmMst2 + 303750*pow2(lmMsq) - 675*(
        989 + 290*lmMst2 - 480*lmMt)*pow2(lmMst1) - 90*lmMst1*(6359 - 10035*
        lmMst2 + 7200*(1 + lmMst2)*lmMt - 16575*pow2(lmMst2)) + 1978425*pow2(
        lmMst2) - 101250*lmMsq*(9 + (6 + 4*lmMst1)*lmMst2 - 2*(pow2(lmMst1) +
        pow2(lmMst2))) + 324000*lmMt*pow2(1 + lmMst2) - 582750*pow3(lmMst1) -
        713250*pow3(lmMst2)) + 125*Dmglst1*(164809 - 864*B4 - 432*DN - 14580*
        lmMsq - 9366*lmMst1 + 4860*pow2(lmMsq) - 6*lmMst2*(-23575 + 1620*lmMsq
        + 906*lmMst1 - 432*pow2(lmMst1)) + 1854*pow2(lmMst1) + 18*(2167 + 336*
        lmMst1)*pow2(lmMst2) + 5184*lmMt*(2 + 3*lmMst2 - lmMst1*(3 + 2*lmMst2)
        + pow2(lmMst1) + pow2(lmMst2)) - 7200*pow3(lmMst1) - 1440*pow3(lmMst2))
        )/pow2(Mst2)))/30375. + (4*(375*(20*Dmglst1*(11 - 6*lmMsq + 6*lmMst2) +
        Mgl*(55 - 30*lmMsq + 6*(5 + 3*lmMst1)*lmMst2 - 9*(pow2(lmMst1) + pow2(
        lmMst2))))*pow4(Mst1) - Mgl*(5108 + 15*lmMsq*(-218 + 75*lmMst1 - 195*
        lmMst2) + 6270*lmMst2 - 375*lmMst1*(8 + 3*lmMst2) + 900*pow2(lmMsq) +
        2025*pow2(lmMst2))*pow4(Mst2)))/(2025.*pow2(Msq)*pow2(Mst2)) - (pow4(
        Mst1)*(1543500*Dmglst1*(2*(17783 - 144*B4 - 72*DN - 4860*lmMsq + 12*(-
        3889 + 450*lmMsq)*lmMst1 - 30483*pow2(lmMst1) + 6*lmMst2*(10610 - 900*
        lmMsq + 6897*lmMst1 + 228*pow2(lmMst1)) + 9*(-891 + 1160*lmMst1)*pow2(
        lmMst2) + 864*lmMt*(4 + 10*lmMst2 - 10*lmMst1*(1 + lmMst2) + 5*pow2(
        lmMst1) + 5*pow2(lmMst2)) - 6120*pow3(lmMst1) - 5688*pow3(lmMst2))*
        pow4(Msq) - 15*(11 + 144*lmMsq - 204*lmMst1 + 60*lmMst2)*pow4(Mst2)) +
        Mgl*(2*(1311202751 + 43575800310*lmMst2 - 132300*(138494 + 9555*lmMst2
        - 15120*lmMt)*pow2(lmMst1) - 1890*lmMst1*(18939979 + 1411200*lmMt +
        280*lmMst2*(-36067 + 7560*lmMt) - 4196850*pow2(lmMst2)) + 1681003800*
        pow2(lmMst2) + 416745000*lmMsq*(-2 + lmMst1*(5 - 2*lmMst2) - 5*lmMst2 +
        pow2(lmMst1) + pow2(lmMst2)) + 666792000*lmMt*(1 + 4*lmMst2 + 3*pow2(
        lmMst2)) - 3134848500*pow3(lmMst1) - 3533071500*pow3(lmMst2))*pow4(Msq)
        + 675*(187967 + 28*lmMst1*(49766 - 13965*lmMst2) + 420*lmMsq*(-2194 +
        259*lmMst1 - 49*lmMst2) - 471968*lmMst2 - 44100*pow2(lmMsq) + 141120*
        pow2(lmMst1) + 205800*pow2(lmMst2))*pow4(Mst2))))/(6.251175e7*pow4(Msq)
        *pow4(Mst2)))/Mgl) + (Mst1*(3*Mgl*(4*(pow2(Mst1)*(1747 - 1929*lmMst2 +
        18*(-29 + 39*lmMst2 - 144*lmMt)*pow2(lmMst1) - 540*lmMsq*(7 + lmMst2 -
        lmMsq*lmMst2 + lmMst1*(-6 + lmMsq + 4*lmMst2 - 6*lmMt) + 5*lmMt + 6*
        lmMst2*lmMt - 4*pow2(lmMst2)) - 5274*pow2(lmMst2) + 18*lmMt*(327 + 570*
        lmMst2 + 80*pow2(lmMst2)) + 3*lmMst1*(985 - 2268*lmMt + 12*lmMst2*(41 +
        32*lmMt) + 54*pow2(lmMst2) - 576*pow2(lmMt)) + 864*(1 + 2*lmMst2)*pow2(
        lmMt) + 726*pow3(lmMst1) - 1590*pow3(lmMst2)) + 6*pow2(Mst2)*(212 +
        225*lmMst2 - 90*lmMsq*(7 + 3*lmMt + lmMst2*(3 - lmMsq + 2*lmMt) +
        lmMst1*(lmMsq - 2*(3 + lmMt))) - 6*(67 + 5*lmMst2 + 14*lmMt)*pow2(
        lmMst1) + 102*pow2(lmMst2) + 3*lmMt*(265 + 202*lmMst2 + 36*pow2(lmMst2)
        ) + 144*(1 + lmMst2)*pow2(lmMt) - 6*lmMst1*(-54 + 37*lmMt + lmMst2*(22
        + 4*lmMt) + pow2(lmMst2) + 24*pow2(lmMt)) + 22*pow3(lmMst1) + 14*pow3(
        lmMst2)))*pow4(Msq) + 5*pow2(Mst2)*(3*(223 - 12*lmMsq + 18*lmMst2*(5 -
        4*lmMt) - 78*lmMt + 72*lmMsq*(lmMst2 + lmMt) - 72*pow2(lmMsq))*pow2(
        Mst1)*pow2(Mst2) + 36*(-15 + lmMst2 + 12*lmMsq*lmMst2 - 7*lmMt - 12*
        lmMst2*lmMt + lmMst1*(6 - 12*lmMsq + 12*lmMt))*pow4(Mst1) + (299 - 60*
        lmMsq + 6*lmMst2*(23 - 12*lmMt) - 78*lmMt + 72*lmMsq*(lmMst2 + lmMt) -
        72*pow2(lmMsq))*pow4(Mst2)) - 80*pow2(Msq)*(9*(11 + 6*lmMst1*(-2 +
        lmMsq - lmMt) + 5*lmMt + lmMst2*(7 - 6*lmMsq + 6*lmMt))*pow2(Mst1)*
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
        243.*Mgl*pow4(Msq)*pow4(Mst2)) - Mst1*pow2(Mt)*pow2(s2t)*(
        739.8333333333334 - (68*B4)/3. - (2*DN)/3. - 240*lmMsq + 40*pow2(lmMsq)
        - lmMst1*(26.333333333333332 - 60*lmMsq + 20*pow2(lmMsq)) + ((-167 +
        60*lmMsq)*pow2(lmMst1))/3. - (lmMst2*(-2135 - 797*lmMst2 + 60*lmMsq*(7
        + lmMst2) + 2*lmMst1*(59 + 2*lmMst2) - 60*pow2(lmMsq) + 60*pow2(lmMst1)
        ))/3. - (100*pow3(lmMst1))/9. + (292*pow3(lmMst2))/9. + (Dmglst1*(
        629.8333333333334 - (68*B4)/3. - (2*DN)/3. - 360*lmMsq + lmMst1*(
        121.66666666666667 + 60*lmMsq - 20*pow2(lmMsq)) + 40*pow2(lmMsq) +
        lmMst2*(777 - 140*lmMsq - (118*lmMst1)/3. + 20*pow2(lmMsq) - 20*pow2(
        lmMst1)) + (-13 + 20*lmMsq)*pow2(lmMst1) - ((-797 + 60*lmMsq + 4*
        lmMst1)*pow2(lmMst2))/3. - (100*pow3(lmMst1))/9. + (292*pow3(lmMst2))/
        9.))/Mgl + pow2(Mst1)*((-20*(Mgl*(71 - 16*lmMst1 + 12*lmMsq*(-2 +
        lmMst1 - lmMst2) + 40*lmMst2 - 6*pow2(lmMst1) + 6*pow2(lmMst2)) + 3*
        Dmglst1*(83 - 20*lmMst1 + 12*lmMsq*(-2 + lmMst1 - lmMst2) + 44*lmMst2 -
        6*pow2(lmMst1) + 6*pow2(lmMst2))))/(9.*Mgl*pow2(Msq)) - (
        77.49382716049382 + 96*B4 - 40*lmMsq + lmMst1*(877.8148148148148 - 200*
        lmMsq + 40*pow2(lmMsq)) + (430.8888888888889 - 120*lmMsq)*pow2(lmMst1)
        - lmMst2*(835.1481481481482 + (820*lmMst1)/9. - 40*lmMsq*(5 + 4*lmMst1)
        + 40*pow2(lmMsq) + 24*pow2(lmMst1)) - (2*(1625 + 180*lmMsq + 276*
        lmMst1)*pow2(lmMst2))/9. + (904*pow3(lmMst1))/9. - (136*pow3(lmMst2))/
        9. + (Dmglst1*(371.14814814814815 + (592*B4)/3. - (8*DN)/3. - 280*lmMsq
        + lmMst1*(431.44444444444446 + 120*lmMsq + 40*pow2(lmMsq)) + (394 -
        280*lmMsq)*pow2(lmMst1) - lmMst2*(120*lmMsq*(1 - 4*lmMst1) + 40*pow2(
        lmMsq) + (1195 + 4764*lmMst1 + 1788*pow2(lmMst1))/9.) + (
        71.33333333333333 - 200*lmMsq - 228*lmMst1)*pow2(lmMst2) + (2708*pow3(
        lmMst1))/9. + (1132*pow3(lmMst2))/9.))/Mgl)/pow2(Mst2)) - pow4(Mst1)*((
        5*(Mgl*(911 - 54*lmMst1 + 12*lmMsq*(-37 + 18*lmMst1 - 18*lmMst2) + 498*
        lmMst2 - 108*pow2(lmMst1) + 108*pow2(lmMst2)) + Dmglst1*(5159 - 702*
        lmMst1 + 12*lmMsq*(-181 + 90*lmMst1 - 90*lmMst2) + 2874*lmMst2 - 540*
        pow2(lmMst1) + 540*pow2(lmMst2))))/(108.*Mgl*pow4(Msq)) + (
        557.5486111111111 + 96*B4 + 70*lmMsq + lmMst1*(875.8981481481482 - 100*
        lmMsq + 40*pow2(lmMsq)) - lmMst2*(939.8981481481482 - (2911*lmMst1)/9.
         - 20*lmMsq*(5 + 8*lmMst1) + 40*pow2(lmMsq) - 146*pow2(lmMst1)) + (
        180.94444444444446 - 120*lmMsq)*pow2(lmMst1) - (504.3888888888889 + 40*
        lmMsq + 146*lmMst1)*pow2(lmMst2) + (46*pow3(lmMst1))/3. - (46*pow3(
        lmMst2))/3. + (Dmglst1*(2730.8603395061727 + (592*B4)/3. - (8*DN)/3. +
        150*lmMsq + lmMst1*(554.0833333333334 + 780*lmMsq + 40*pow2(lmMsq)) -
        lmMst2*(60*lmMsq*(13 - 8*lmMst1) + 40*pow2(lmMsq) + (31467 + 9076*
        lmMst1 - 6744*pow2(lmMst1))/36.) - (5*(463 + 1008*lmMsq)*pow2(lmMst1))/
        18. - ((-6853 + 3600*lmMsq + 9516*lmMst1)*pow2(lmMst2))/18. + (1294*
        pow3(lmMst1))/9. + (1778*pow3(lmMst2))/9.))/Mgl)/pow4(Mst2)) + (5*(2*(
        Dmglst1*(131 - 12*lmMsq + 12*lmMst2) + (59 - 12*lmMsq + 12*lmMst2)*Mgl)
        *pow2(Mst1)*pow4(Mst2) + 48*pow2(Msq)*((2*Mgl*(6 - 34*lmMst2 + 12*
        lmMsq*lmMst2 - 2*lmMst1*(-17 + 6*lmMsq + 9*lmMst2) + 15*pow2(lmMst1) +
        3*pow2(lmMst2)) + 6*Dmglst1*(12 - 26*lmMst2 + 12*lmMsq*lmMst2 - 2*
        lmMst1*(-13 + 6*lmMsq + 15*lmMst2) + 21*pow2(lmMst1) + 9*pow2(lmMst2)))
        *pow4(Mst1) + 3*(Dmglst1 + Mgl)*pow4(Mst2)) + (-11 + 12*lmMsq - 12*
        lmMst2)*(Dmglst1 + Mgl)*pow6(Mst2)))/(108.*Mgl*pow2(Mst2)*pow4(Msq))) +
        Mt*s2t*(-((-1 + 2*lmMst2)*shiftst3*pow2(Mst2)*(-4*pow2(Mt) + ((-1 -
        lmMst1 + lmMst2)*pow2(Mst1) + pow2(Mst2))*pow2(s2t)))/(3.*pow2(Mst1)) +
        pow2(s2t)*(pow2(Mst1)*(110.18859670781893 - (40*B4)/9. - (2*DN)/9. - (
        80*lmMsq)/9. + (10*pow2(lmMsq))/3. + lmMst1*(189.21814814814815 - (355*
        lmMsq)/9. + (35*pow2(lmMsq))/3.) + (4.705555555555556 - (5*lmMsq)/3.)*
        pow2(lmMst1) - (lmMst2*(393289 - 397890*lmMst1 + 1500*lmMsq*(-59 + 36*
        lmMst1) + 31500*pow2(lmMsq) + 35850*pow2(lmMst1)))/2700. + ((-8151 +
        1300*lmMsq + 3890*lmMst1)*pow2(lmMst2))/60. + (49*pow3(lmMst1))/54. - (
        2833*pow3(lmMst2))/54. - (Dmglst1*(2323 + 144*B4 + 120*DN - 2880*lmMsq
        + 54*lmMst1*(32.55555555555556 + 40*lmMsq - 20*pow2(lmMsq)) + 6*lmMst2*
        (571 - 360*lmMsq + 62*lmMst1 + 180*pow2(lmMsq) - 158*pow2(lmMst1)) +
        30*(-121 + 36*lmMsq)*pow2(lmMst1) - 6*(-735 + 180*lmMsq + 98*lmMst1)*
        pow2(lmMst2) - 260*pow3(lmMst1) + 1796*pow3(lmMst2)))/(54.*Mgl)) + ((
        50*Dmglst1*(192 - 91*lmMst1 + 6*lmMsq*(-6 + 17*lmMst1 - 17*lmMst2) +
        127*lmMst2 - 51*pow2(lmMst1) + 51*pow2(lmMst2)) + 3*Mgl*(787 + 20*
        lmMsq*(-9 + 20*lmMst1 - 20*lmMst2) + 525*lmMst2 - 5*lmMst1*(69 + 70*
        lmMst2) - 25*pow2(lmMst1) + 375*pow2(lmMst2)))/(135.*Mgl*pow2(Msq)) + (
        239.07595982905215 - (71872687*lmMst2)/529200. - ((488263 + 651840*
        lmMst2)*pow2(lmMst1))/7560. - (770503*pow2(lmMst2))/7560. + (37*lmMst1*
        (2891251 + 2673860*lmMst2 + 2352000*pow2(lmMst2)))/529200. + (8*pow3(
        lmMst1))/3. - (728*pow3(lmMst2))/9. + (5*(2*Dmglst1 + Mgl)*pow2(lmMsq)
        - (5*lmMsq*(4*Dmglst1*(11 + 30*lmMst1 - 24*lmMst2) + Mgl*(29 + 30*
        lmMst1 + 2*(-9 + 16*lmMst1)*lmMst2 - 16*(pow2(lmMst1) + pow2(lmMst2))))
        )/6. - (Dmglst1*(-606133 + 624756*lmMst2 + 72*(-3901 + 2976*lmMst2)*
        pow2(lmMst1) + 361944*pow2(lmMst2) - 12*lmMst1*(92887 + 23460*lmMst2 +
        36288*pow2(lmMst2)) + 2304*pow3(lmMst1) + 218880*pow3(lmMst2)))/3888.)/
        Mgl)/pow2(Mst2))*pow4(Mst1) - (pow2(Mst2)*(Mgl*(82320*(1775 + 30*lmMsq*
        (-5 + lmMst1 - lmMst2) + 407*lmMst2 + lmMst1*(-257 + 600*lmMst2) - 315*
        pow2(lmMst1) - 285*pow2(lmMst2))*pow2(Msq)*pow2(Mst1) + 17150*(53965 +
        2592*B4 - 288*DN - 13320*lmMsq + 5400*pow2(lmMsq) - 6*lmMst1*(3781 -
        300*lmMsq + 180*pow2(lmMsq)) + 6*lmMst2*(14137 - 3498*lmMst1 + 60*
        lmMsq*(-35 + 12*lmMst1) + 180*pow2(lmMsq) - 18*pow2(lmMst1)) - 18*(-229
        + 60*lmMsq)*pow2(lmMst1) - 18*(-2193 + 180*lmMsq + 442*lmMst1)*pow2(
        lmMst2) + 396*pow3(lmMst1) + 7668*pow3(lmMst2))*pow4(Msq) + 9*(6740969
        + 140*lmMsq*(-12899 + 6230*lmMst1 - 5390*lmMst2) + 2822890*lmMst2 + 70*
        lmMst1*(-14529 + 25480*lmMst2) - 58800*pow2(lmMsq) - 1327900*pow2(
        lmMst1) - 514500*pow2(lmMst2))*pow4(Mst1)) + 102900*Dmglst1*(40*(82 -
        93*lmMst1 + 6*lmMsq*(4 + 3*lmMst1 - 3*lmMst2) + 69*lmMst2 - 9*pow2(
        lmMst1) + 9*pow2(lmMst2))*pow2(Msq)*pow2(Mst1) + 6*(7 + 60*lmMsq*(5 -
        6*lmMst1) - 2*lmMst1 + (-226 + 220*lmMst1)*lmMst2 + 180*pow2(lmMsq) +
        178*pow2(lmMst1) + 18*pow2(lmMst2))*pow4(Msq) + 5*(503 - 346*lmMst1 +
        12*lmMsq*(1 + 10*lmMst1 - 10*lmMst2) + 334*lmMst2 - 60*pow2(lmMst1) +
        60*pow2(lmMst2))*pow4(Mst1))))/(1.11132e7*Mgl*pow4(Msq)) + ((939 + 30*
        lmMsq*(-12 + 5*lmMst1 - 5*lmMst2) + 760*lmMst2 - 50*lmMst1*(8 + 3*
        lmMst2) + 150*pow2(lmMst2))/(135.*pow2(Msq)) + (-1 + 2*lmMst2)/(3.*
        pow2(Mst1)) + ((1.0841585681891803 + lmMsq*(-0.6291383219954648 + (7*
        lmMst1)/9. - (43*lmMst2)/63.) + (47419*lmMst2)/26460. - (lmMst1*(157 +
        105*lmMst2))/135. + (5*Dmglst1)/(6.*Mgl) - pow2(lmMsq)/21. + (46*pow2(
        lmMst2))/63.)*pow2(Mst1))/pow4(Msq))*pow4(Mst2) + (-177*(8*OepS2 - 81*(
        5 + 2*lmMst1 - 2*lmMst2)*S2)*pow2(Mst1)*pow2(Mst2) + (-1328*OepS2 + 54*
        (1187 + 498*lmMst1 - 498*lmMst2)*S2)*pow4(Mst1) - 27*(8*OepS2 + 81*(11
        - 2*lmMst1 + 2*lmMst2)*S2)*pow4(Mst2))/(729.*pow2(Mst2))) - (4*T1ep*(2*
        (40*pow2(Mt) - 83*pow2(Mst2)*pow2(s2t))*pow4(Mst1) + 54*pow2(Mt)*pow4(
        Mst2) + pow2(Mst1)*(84*pow2(Mst2)*pow2(Mt) - 177*pow2(s2t)*pow4(Mst2))
        - 27*pow2(s2t)*pow6(Mst2)))/(243.*pow4(Mst2)))) - (324*xMsq*z2*pow2(
        Mgl)*pow2(Msq)*(20*Mt*MuSUSY*s2t*pow2(Mst2)*pow2(Sbeta)*(4*(-1 +
        shiftst2)*pow2(Mst1)*pow2(Mt) - (shiftst1 - shiftst2 + (lmMst1 -
        lmMst2)*(-2 + shiftst1 + shiftst2))*pow2(Mst1)*pow2(Mst2)*pow2(s2t) + (
        -1 + shiftst1)*pow2(Mst2)*(-4*pow2(Mt) + pow2(Mst2)*pow2(s2t)) - (-1 +
        shiftst2)*pow2(s2t)*pow4(Mst1)) + Tbeta*(pow2(Mt)*pow2(s2t)*(-20*(2*(
        lmMst1 - lmMst2)*pow2(Mst1) - pow2(Mst2))*(pow2(Mst1) + pow2(Mst2))*
        pow2(MuSUSY) - pow2(Sbeta)*(-20*(2*(lmMst1 - lmMst2)*pow2(Mst1) - pow2(
        Mst2))*(pow2(Mst1) + pow2(Mst2))*pow2(MuSUSY) + 40*pow2(Mst2)*(pow4(
        Mst1) + pow4(Mst2)))) + pow2(Mst2)*pow2(Sbeta)*(80*(pow2(Mst1)*pow2(
        Mst2)*pow2(Mt)*pow2(s2t) + (pow2(Mst1) + pow2(Mst2))*pow4(Mt)) + 5*(
        pow2(Mst1) - pow2(Mst2))*(2*(lmMst1 - lmMst2)*pow2(Mst1)*pow2(Mst2) +
        pow4(Mst1) - pow4(Mst2))*pow4(s2t)) - 5*shiftst1*(4*pow2(Mt)*pow2(
        MuSUSY)*pow2(s2t)*(2*(1 - lmMst1 + lmMst2)*pow2(Mst1)*pow2(Mst2) + (2 -
        4*lmMst1 + 4*lmMst2)*pow4(Mst1) + pow4(Mst2)) + pow2(Sbeta)*(-4*pow2(
        Mt)*pow2(s2t)*(2*(-1 + lmMst1 - lmMst2)*pow2(Mst1)*pow2(Mst2)*(pow2(
        Mst2) - pow2(MuSUSY)) + 2*(1 - 2*lmMst1 + 2*lmMst2)*pow2(MuSUSY)*pow4(
        Mst1) + (2*pow2(Mst2) + pow2(MuSUSY))*pow4(Mst2)) + pow4(Mst2)*(16*
        pow4(Mt) + (-4*pow2(Mst1)*pow2(Mst2) + (3 + 2*lmMst1 - 2*lmMst2)*pow4(
        Mst1) + pow4(Mst2))*pow4(s2t)))) + shiftst2*pow2(Mst1)*(20*(2*(1 -
        lmMst1 + lmMst2)*pow2(Mst1) + pow2(Mst2))*pow2(Mt)*pow2(MuSUSY)*pow2(
        s2t) + pow2(Sbeta)*(-(pow2(Mt)*(20*pow2(Mst2)*(2*(1 + lmMst1 - lmMst2)*
        pow2(Mst2) + pow2(MuSUSY)) - 40*pow2(Mst1)*(pow2(Mst2) + (-1 + lmMst1 -
        lmMst2)*pow2(MuSUSY)))*pow2(s2t)) + pow2(Mst2)*(-80*pow4(Mt) - 5*(-4*
        pow2(Mst1)*pow2(Mst2) + pow4(Mst1) + (3 - 2*lmMst1 + 2*lmMst2)*pow4(
        Mst2))*pow4(s2t)))))) + pow2(Mst1)*(-(z3*(54*xDmglst1*pow2(Dmglst1)*(
        Tbeta*(pow2(Mt)*pow2(MuSUSY)*(8*Mst1*Mt*(-4*(220 + 9*lmMst1 - 9*lmMst2)
        *Mst1*Mt + 3*(151 - 918*lmMst1 + 918*lmMst2)*s2t*pow2(Mst1) + (559 +
        18*lmMst1 - 18*lmMst2)*s2t*pow2(Mst2)) + pow2(s2t)*(24*(49 + 69*lmMst1
        - 69*lmMst2)*pow2(Mst1)*pow2(Mst2) + 4*(6491 + 318*lmMst1 - 318*lmMst2)
        *pow4(Mst1) - 1533*pow4(Mst2))) + pow2(Sbeta)*(pow2(Mt)*pow2(s2t)*(4*(
        4378*pow2(Mst2) + (-6491 - 318*lmMst1 + 318*lmMst2)*pow2(MuSUSY))*pow4(
        Mst1) + 7*(497*pow2(Mst2) + 219*pow2(MuSUSY))*pow4(Mst2) - pow2(Mst1)*(
        24*(49 + 69*lmMst1 - 69*lmMst2)*pow2(Mst2)*pow2(MuSUSY) + (16063 + 288*
        lmMst1 - 288*lmMst2)*pow4(Mst2))) + Mst1*(4*Mt*pow2(Mst2)*pow3(s2t)*(-
        5*(133 + 558*lmMst1 - 558*lmMst2)*pow2(Mst1)*pow2(Mst2) + 9559*pow4(
        Mst1) + (559 + 18*lmMst1 - 18*lmMst2)*pow4(Mst2)) - 8*s2t*pow3(Mt)*((
        559 + 18*lmMst1 - 18*lmMst2)*pow2(Mst2)*pow2(MuSUSY) + 3*pow2(Mst1)*(
        1260*pow2(Mst2) + (151 - 918*lmMst1 + 918*lmMst2)*pow2(MuSUSY)) +
        25672*pow4(Mst1) + 3*(565 + 16*lmMst1 - 16*lmMst2)*pow4(Mst2))) + 2*(
        16*pow2(Mst1)*(1163*pow2(Mst2) + (220 + 9*lmMst1 - 9*lmMst2)*pow2(
        MuSUSY)) + 26736*pow4(Mst1) - 2019*pow4(Mst2))*pow4(Mt) + ((6*(707 +
        276*lmMst1 - 276*lmMst2)*pow2(Mst1)*pow2(Mst2) + (22079 - 2040*lmMst1 +
        2040*lmMst2)*pow4(Mst1) - 1533*pow4(Mst2))*pow4(Mst2)*pow4(s2t))/4.)) +
        Mt*MuSUSY*pow2(Sbeta)*(24*((565 + 16*lmMst1 - 16*lmMst2)*Mst1*pow2(
        Mst2) + (1825 + 16*lmMst1 - 16*lmMst2)*pow3(Mst1))*pow3(Mt) - 3*pow2(
        Mst2)*pow3(s2t)*(3*(301 + 184*lmMst1 - 184*lmMst2)*pow2(Mst1)*pow2(
        Mst2) + (4*(6197 - 96*lmMst1 + 96*lmMst2)*pow4(Mst1))/3. - 511*pow4(
        Mst2)) - 2*s2t*pow2(Mt)*(-8*(1573 + 36*lmMst1 - 36*lmMst2)*pow2(Mst1)*
        pow2(Mst2) + 32*(154 - 9*lmMst1 + 9*lmMst2)*pow4(Mst1) + 3479*pow4(
        Mst2)) - 18*Mt*pow2(s2t)*((-4*(53 + 1386*lmMst1 - 1386*lmMst2)*pow2(
        Mst2)*pow3(Mst1))/3. + (2*(559 + 18*lmMst1 - 18*lmMst2)*Mst1*pow4(Mst2)
        )/3. + (6302 - 1848*lmMst1 + 1848*lmMst2)*pow5(Mst1)))) + Mgl*(4*Mt*
        MuSUSY*pow2(Sbeta)*(432*Mst1*((2*Dmglst1*(235 + 6*lmMst1 - 6*lmMst2) +
        3*(43 + 4*lmMst1 - 4*lmMst2)*Mgl)*pow2(Mst1) - (Dmglst1*(37 - 12*lmMst1
        + 12*lmMst2) + 6*(9 - 2*lmMst1 + 2*lmMst2)*Mgl)*pow2(Mst2))*pow3(Mt) -
        pow2(Mst2)*pow3(s2t)*(3*(27*Dmglst1*(205 + 184*lmMst1 - 184*lmMst2) +
        2*(7307 + 1134*lmMst1 - 1134*lmMst2)*Mgl)*pow2(Mst1)*pow2(Mst2) + 2*(
        76842*Dmglst1 + 41257*Mgl)*pow4(Mst1) - 27*(195*Dmglst1 - 328*Mgl)*
        pow4(Mst2)) - 4*s2t*pow2(Mt)*(-6*(54*Dmglst1*(53 + 4*lmMst1 - 4*lmMst2)
        + 355*Mgl)*pow2(Mst1)*pow2(Mst2) + 4*(108*Dmglst1*(10 - 3*lmMst1 + 3*
        lmMst2) - 367*Mgl)*pow4(Mst1) + 27*(110*Dmglst1 - 43*Mgl)*pow4(Mst2)) -
        81*Mt*pow2(s2t)*(-8*(Dmglst1*(202 + 342*lmMst1 - 342*lmMst2) + (221 +
        108*lmMst1 - 108*lmMst2)*Mgl)*pow2(Mst2)*pow3(Mst1) - 4*(Dmglst1*(58 -
        9*lmMst1 + 9*lmMst2) + (133 - 9*lmMst1 + 9*lmMst2)*Mgl)*Mst1*pow4(Mst2)
        + (Dmglst1*(4366 - 2736*lmMst1 + 2736*lmMst2) - 2*(365 + 432*lmMst1 -
        432*lmMst2)*Mgl)*pow5(Mst1))) + Tbeta*(-4*pow2(Mt)*pow2(MuSUSY)*(216*
        Mst1*Mt*(Mst1*(6*Dmglst1*(13 + 2*lmMst1 - 2*lmMst2)*Mt + (-46 + 6*
        lmMst1 - 6*lmMst2)*Mgl*Mt + 3*Dmglst1*(154 + 225*lmMst1 - 225*lmMst2)*
        Mst1*s2t + 23*(25 + 9*lmMst1 - 9*lmMst2)*Mgl*Mst1*s2t) + (Dmglst1*(58 -
        9*lmMst1 + 9*lmMst2) + (133 - 9*lmMst1 + 9*lmMst2)*Mgl)*s2t*pow2(Mst2))
        - pow2(s2t)*(6*(54*Dmglst1*(35 + 46*lmMst1 - 46*lmMst2) + (8783 + 1134*
        lmMst1 - 1134*lmMst2)*Mgl)*pow2(Mst1)*pow2(Mst2) + 4*(54*Dmglst1*(764 +
        69*lmMst1 - 69*lmMst2) + 7*(4829 + 243*lmMst1 - 243*lmMst2)*Mgl)*pow4(
        Mst1) + (-5265*Dmglst1 + 8856*Mgl)*pow4(Mst2))) + pow2(Sbeta)*(-108*
        Mst1*s2t*pow3(Mt)*(-8*(Dmglst1*(58 - 9*lmMst1 + 9*lmMst2) + (133 - 9*
        lmMst1 + 9*lmMst2)*Mgl)*pow2(Mst2)*pow2(MuSUSY) + pow2(Mst1)*(48*(169*
        Dmglst1 + 61*Mgl)*pow2(Mst2) - 8*(Dmglst1*(462 + 675*lmMst1 - 675*
        lmMst2) + 23*(25 + 9*lmMst1 - 9*lmMst2)*Mgl)*pow2(MuSUSY)) + 64*(549*
        Dmglst1 + 112*Mgl)*pow4(Mst1) - 16*(Dmglst1*(37 - 12*lmMst1 + 12*
        lmMst2) + 6*(9 - 2*lmMst1 + 2*lmMst2)*Mgl)*pow4(Mst2)) + pow4(Mst2)*(6*
        (27*Dmglst1*(135 + 92*lmMst1 - 92*lmMst2) + 7*(833 + 162*lmMst1 - 162*
        lmMst2)*Mgl)*pow2(Mst1)*pow2(Mst2) + 4*((27*Dmglst1*(5077 - 552*lmMst1
        + 552*lmMst2))/4. + (9668 - 1701*lmMst1 + 1701*lmMst2)*Mgl)*pow4(Mst1)
        - 27*(195*Dmglst1 - 328*Mgl)*pow4(Mst2))*pow4(s2t) + 108*((2*Dmglst1*(
        8*pow2(Mst1)*(386*pow2(Mst2) + 3*(13 + 2*lmMst1 - 2*lmMst2)*pow2(
        MuSUSY)) + 5264*pow4(Mst1) - 107*pow4(Mst2)) + 8*Mgl*(2*pow2(Mst1)*(72*
        pow2(Mst2) + (-23 + 3*lmMst1 - 3*lmMst2)*pow2(MuSUSY)) + 317*pow4(Mst1)
        + (-163 + 6*lmMst1 + 6*lmMst2 - 12*lmMt)*pow4(Mst2)))*pow4(Mt) + Mt*
        pow2(Mst2)*pow3(s2t)*(-4*((Dmglst1*(346 + 693*lmMst1 - 693*lmMst2) + 3*
        (103 + 75*lmMst1 - 75*lmMst2)*Mgl)*pow2(Mst2)*pow3(Mst1) + (Dmglst1*(58
        - 9*lmMst1 + 9*lmMst2) + (133 - 9*lmMst1 + 9*lmMst2)*Mgl)*Mst1*pow4(
        Mst2)) + 6*(997*Dmglst1 + 173*Mgl)*pow5(Mst1))) + 4*pow2(Mt)*pow2(s2t)*
        (4*(10746*Dmglst1 + 331*Mgl)*pow2(Mst2)*pow4(Mst1) - 6*(18*Dmglst1*(373
        + 24*lmMst1 - 24*lmMst2) + 323*Mgl)*pow2(Mst1)*pow4(Mst2) - pow2(
        MuSUSY)*(6*(54*Dmglst1*(35 + 46*lmMst1 - 46*lmMst2) + (8783 + 1134*
        lmMst1 - 1134*lmMst2)*Mgl)*pow2(Mst1)*pow2(Mst2) + 4*(54*Dmglst1*(764 +
        69*lmMst1 - 69*lmMst2) + 7*(4829 + 243*lmMst1 - 243*lmMst2)*Mgl)*pow4(
        Mst1) + (-5265*Dmglst1 + 8856*Mgl)*pow4(Mst2)) + 54*(110*Dmglst1 - 43*
        Mgl)*pow6(Mst2))))))) + 3*pow2(z2)*(36*xDmglst1*pow2(Dmglst1)*(8*Mt*
        MuSUSY*pow2(Sbeta)*(8*Mst1*(pow2(Mst1) + pow2(Mst2))*pow3(Mt) - 3*pow2(
        Mst1)*pow3(s2t)*pow4(Mst2) + 24*Mt*pow2(s2t)*(17*(pow2(Mst1) + pow2(
        Mst2))*pow3(Mst1) + Mst1*pow4(Mst2))) + 2*Mst1*s2t*Tbeta*(s2t*pow2(
        Sbeta)*(-12*Mst1*(pow2(Mst1) + pow2(Mst2))*pow2(Mt)*pow2(MuSUSY) + s2t*
        (-512*Mt*pow2(Mst1) - 32*Mt*pow2(Mst2) + 3*Mst1*s2t*pow2(Mst2) - 3*s2t*
        pow3(Mst1))*pow4(Mst2)) + 4*pow2(Mt)*(-((-3*Mst1*s2t*(pow2(Mst1) +
        pow2(Mst2)) + 16*Mt*(18*pow2(Mst1) + pow2(Mst2)))*pow2(MuSUSY)) - 8*Mt*
        pow2(Sbeta)*(-2*(18*pow2(Mst1) + pow2(Mst2))*pow2(MuSUSY) + pow4(Mst2))
        ))) + Mgl*(4*Mt*MuSUSY*pow2(Sbeta)*(576*(Dmglst1 + Mgl)*Mst1*(pow2(
        Mst1) + pow2(Mst2))*pow3(Mt) + 2*Mgl*s2t*pow2(Mt)*(42*pow2(Mst1)*pow2(
        Mst2) + 40*pow4(Mst1) - 45*pow4(Mst2)) + 864*Mt*pow2(s2t)*((19*Dmglst1
        + 9*Mgl)*(pow2(Mst1) + pow2(Mst2))*pow3(Mst1) + 2*(Dmglst1 + Mgl)*Mst1*
        pow4(Mst2)) - pow2(Mst2)*pow3(s2t)*(3*(48*Dmglst1 + 179*Mgl)*pow2(Mst1)
        *pow2(Mst2) + Mgl*(166*pow4(Mst1) + 315*pow4(Mst2)))) + s2t*Tbeta*(4*
        pow2(Mt)*(-576*Mt*pow2(Sbeta)*(-(pow2(MuSUSY)*(2*(Dmglst1 + Mgl)*Mst1*
        pow2(Mst2) + (21*Dmglst1 + 11*Mgl)*pow3(Mst1))) + (Dmglst1 + Mgl)*Mst1*
        pow4(Mst2)) + pow2(MuSUSY)*(-576*Mt*(2*(Dmglst1 + Mgl)*Mst1*pow2(Mst2)
        + (21*Dmglst1 + 11*Mgl)*pow3(Mst1)) + s2t*(12*(12*Dmglst1 + 71*Mgl)*
        pow2(Mst1)*pow2(Mst2) + 2*(72*Dmglst1 + 509*Mgl)*pow4(Mst1) + 315*Mgl*
        pow4(Mst2)))) + s2t*pow2(Sbeta)*(s2t*pow4(Mst2)*(-1152*Mt*(2*(Dmglst1 +
        Mgl)*Mst1*pow2(Mst2) + (17*Dmglst1 + 7*Mgl)*pow3(Mst1)) + s2t*(6*(24*
        Dmglst1 + 37*Mgl)*pow2(Mst1)*pow2(Mst2) - (144*Dmglst1 + 371*Mgl)*pow4(
        Mst1) + 315*Mgl*pow4(Mst2))) + 4*pow2(Mt)*(-(pow2(MuSUSY)*(12*(12*
        Dmglst1 + 71*Mgl)*pow2(Mst1)*pow2(Mst2) + 2*(72*Dmglst1 + 509*Mgl)*
        pow4(Mst1) + 315*Mgl*pow4(Mst2))) + Mgl*(2*pow2(Mst2)*pow4(Mst1) - 87*
        pow2(Mst1)*pow4(Mst2) + 45*pow6(Mst2))))))) - 2*z4*(54*xDmglst1*pow2(
        Dmglst1)*(2*Mt*MuSUSY*(Mst1*Mt*MuSUSY*Tbeta*(Mt*s2t*(3474*pow2(Mst1) -
        482*pow2(Mst2)) + 180*Mst1*pow2(Mt) - 15*pow2(s2t)*(53*Mst1*pow2(Mst2)
        + 37*pow3(Mst1))) + pow2(Sbeta)*(-8*Mst1*(-2*Mt + 45*Mst1*s2t)*(pow2(
        Mst1) + pow2(Mst2))*pow2(Mt) + 15*pow2(Mst2)*pow3(s2t)*(53*pow2(Mst1)*
        pow2(Mst2) - 16*pow4(Mst1)) + 3*Mt*pow2(s2t)*(-1978*(pow2(Mst1) + pow2(
        Mst2))*pow3(Mst1) + 241*Mst1*pow4(Mst2)))) + Tbeta*pow2(Sbeta)*((Mst1*(
        8876*Mt*pow2(Mst1) - 964*Mt*pow2(Mst2) - 795*Mst1*s2t*pow2(Mst2) +
        1035*s2t*pow3(Mst1))*pow3(s2t)*pow4(Mst2))/2. - 4*s2t*pow3(Mt)*(pow2(
        MuSUSY)*(-241*Mst1*pow2(Mst2) + 1737*pow3(Mst1)) + 8*Mst1*pow4(Mst2)) +
        pow2(Mst1)*(30*pow2(Mt)*pow2(s2t)*((37*pow2(Mst1) + 53*pow2(Mst2))*
        pow2(MuSUSY) + 12*pow4(Mst2)) - 360*pow2(MuSUSY)*pow4(Mt)))) + Mgl*(4*
        Mt*MuSUSY*pow2(Sbeta)*(432*(Dmglst1 + Mgl)*Mst1*(pow2(Mst1) + pow2(
        Mst2))*pow3(Mt) + pow2(Mst2)*pow3(s2t)*(6*(2385*Dmglst1 + 538*Mgl)*
        pow2(Mst1)*pow2(Mst2) + Mgl*(166*pow4(Mst1) - 3267*pow4(Mst2))) - 2*
        s2t*pow2(Mt)*(6*(540*Dmglst1 + 7*Mgl)*pow2(Mst1)*pow2(Mst2) + 40*(81*
        Dmglst1 + Mgl)*pow4(Mst1) - 243*Mgl*pow4(Mst2)) + 81*Mt*pow2(s2t)*(-4*(
        187*Dmglst1 - 18*Mgl)*(pow2(Mst1) + pow2(Mst2))*pow3(Mst1) + 241*(
        Dmglst1 + Mgl)*Mst1*pow4(Mst2))) + Tbeta*(pow2(MuSUSY)*(-4*pow2(Mt)*(
        54*Mst1*Mt*(-(Mst1*(30*(2*Dmglst1 + Mgl)*Mt + (507*Dmglst1 - 313*Mgl)*
        Mst1*s2t)) + 241*(Dmglst1 + Mgl)*s2t*pow2(Mst2)) + pow2(s2t)*(14310*
        Dmglst1*pow2(Mst1)*(pow2(Mst1) + pow2(Mst2)) + Mgl*(-39*pow2(Mst1)*
        pow2(Mst2) + 127*pow4(Mst1) - 3267*pow4(Mst2)))) - 6480*(2*Dmglst1 +
        Mgl)*pow2(Mst1)*pow2(Sbeta)*pow4(Mt)) + s2t*pow2(Sbeta)*(-216*pow3(Mt)*
        (pow2(MuSUSY)*(-241*(Dmglst1 + Mgl)*Mst1*pow2(Mst2) + (507*Dmglst1 -
        313*Mgl)*pow3(Mst1)) + 8*(Dmglst1 + Mgl)*Mst1*pow4(Mst2)) + pow2(s2t)*
        pow4(Mst2)*(108*Mt*(-241*(Dmglst1 + Mgl)*Mst1*pow2(Mst2) + (989*Dmglst1
        + 169*Mgl)*pow3(Mst1)) + s2t*(-15*(954*Dmglst1 + 433*Mgl)*pow2(Mst1)*
        pow2(Mst2) + 2*(7155*Dmglst1 + 1531*Mgl)*pow4(Mst1) + 3267*Mgl*pow4(
        Mst2))) - 4*s2t*pow2(Mt)*(-(pow2(MuSUSY)*(14310*Dmglst1*pow2(Mst1)*(
        pow2(Mst1) + pow2(Mst2)) + Mgl*(-39*pow2(Mst1)*pow2(Mst2) + 127*pow4(
        Mst1) - 3267*pow4(Mst2)))) - 15*(216*Dmglst1 + 19*Mgl)*pow2(Mst1)*pow4(
        Mst2) + Mgl*(2*pow2(Mst2)*pow4(Mst1) + 243*pow6(Mst2)))))))))/(972.*
        pow2(Mgl)*pow2(Mst1)*pow2(Sbeta)*pow4(Mst2)))/Tbeta + (pow2(s2t)*((-2*
        pow2(Mt)*(pow2(MuSUSY)*(81*(-1 + 3*lmMst1 - 3*lmMst2)*(-1 + 2*lmMst2)*
        shiftst3 - 740*T1ep*(-1 + pow2(Sbeta))) + 4*T1ep*pow2(Mst2)*pow2(Sbeta)
        ) + 11*T1ep*pow2(s2t)*pow2(Sbeta)*pow4(Mst2))*pow6(Mst1) + 81*(-1 + 2*
        lmMst2)*shiftst3*pow2(Mt)*pow2(MuSUSY)*pow6(Mst2) - 6*pow4(Mst1)*(pow2(
        Mst2)*pow2(Mt)*(pow2(MuSUSY)*(27*shiftst3*(1 - 2*lmMst1 + 4*lmMst1*
        lmMst2 - 4*pow2(lmMst2)) - 136*T1ep*(-1 + pow2(Sbeta))) - 10*T1ep*pow2(
        Mst2)*pow2(Sbeta)) + 25*T1ep*pow2(s2t)*pow2(Sbeta)*pow6(Mst2)) - 27*
        pow2(Mst1)*(2*pow2(Mt)*(pow2(MuSUSY)*(3*(-1 + lmMst1 - lmMst2)*(-1 + 2*
        lmMst2)*shiftst3 - 2*T1ep*(-1 + pow2(Sbeta))) - 2*T1ep*pow2(Mst2)*pow2(
        Sbeta))*pow4(Mst2) + T1ep*pow2(s2t)*pow2(Sbeta)*pow8(Mst2))))/(243.*
        pow2(Mst1)*pow2(Sbeta)*pow4(Mst2)) + ((1 - 2*lmMst2)*shiftst3*Tbeta*
        pow2(Sbeta)*(-((8*(1 - lmMst1 + lmMst2)*pow2(Mst1)*pow2(Mt)*pow2(s2t) +
        16*pow4(Mt) + pow2(Mst1)*((3 + 2*lmMst1 - 2*lmMst2)*pow2(Mst1) - 4*
        pow2(Mst2))*pow4(s2t))*pow6(Mst2)) + 4*pow2(Mt)*pow2(s2t)*(pow2(MuSUSY)
        *pow6(Mst2) + 2*(pow2(MuSUSY)*((1 - 2*lmMst1 + 2*lmMst2)*pow2(Mst2)*
        pow4(Mst1) + (1 - lmMst1 + lmMst2)*pow2(Mst1)*pow4(Mst2) + (1 - 3*
        lmMst1 + 3*lmMst2)*pow6(Mst1)) + pow8(Mst2)))) + 10*(1 - 2*lmMsq)*xMsq*
        pow2(Msq)*(4*Mt*MuSUSY*s2t*(Mt*MuSUSY*s2t*Tbeta*pow2(Mst1)*(pow2(Mst2)
        - 2*(lmMst1 - lmMst2)*(pow2(Mst1) + pow2(Mst2))) + pow2(Mst2)*pow2(
        Sbeta)*((-1 + shiftst1)*pow2(Mst2)*(-4*pow2(Mt) + pow2(Mst2)*pow2(s2t))
        - pow2(Mst1)*(-4*(-1 + shiftst2)*pow2(Mt) + (shiftst1 - shiftst2 + (
        lmMst1 - lmMst2)*(-2 + shiftst1 + shiftst2))*pow2(Mst2)*pow2(s2t)) - (-
        1 + shiftst2)*pow2(s2t)*pow4(Mst1))) + Tbeta*(-8*pow2(Mst1)*pow2(Mst2)*
        pow2(Mt)*((1 - lmMst1 + lmMst2)*shiftst1*pow2(MuSUSY)*pow2(s2t) + 2*
        shiftst2*pow2(Mt)*pow2(Sbeta)) + pow2(Sbeta)*(-4*pow2(Mt)*pow2(s2t)*(2*
        (pow2(Mst2) + (-lmMst1 + lmMst2)*pow2(MuSUSY))*pow4(Mst1) + pow2(Mst1)*
        ((1 - 2*lmMst1 + 2*lmMst2)*pow2(Mst2)*pow2(MuSUSY) - 4*pow4(Mst2)) + (
        2*pow2(Mst2) + pow2(MuSUSY))*pow4(Mst2)) + pow2(Mst2)*(16*(pow2(Mst1) +
        pow2(Mst2))*pow4(Mt) + pow2(Mst1)*((-1 + 2*lmMst1 - 2*lmMst2)*pow2(
        Mst1)*pow2(Mst2) + pow4(Mst1) + (-1 - 2*lmMst1 + 2*lmMst2)*pow4(Mst2))*
        pow4(s2t))) - shiftst1*(4*pow2(Mst1)*pow2(Mst2)*pow2(s2t)*pow2(Sbeta)*(
        2*(1 - lmMst1 + lmMst2)*pow2(Mt)*(pow2(Mst2) - pow2(MuSUSY)) - pow2(
        s2t)*pow4(Mst2)) + pow4(Mst1)*(8*(-1 + 2*lmMst1 - 2*lmMst2)*pow2(Mt)*
        pow2(MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta)) + (3 + 2*lmMst1 - 2*lmMst2)*
        pow2(Sbeta)*pow4(Mst2)*pow4(s2t)) + pow4(Mst2)*(-4*pow2(Mt)*pow2(s2t)*(
        pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 2*pow2(Mst2)*pow2(Sbeta)) + pow2(
        Sbeta)*(16*pow4(Mt) + pow4(Mst2)*pow4(s2t)))) + pow2(s2t)*(4*pow2(Mt)*
        pow2(MuSUSY)*pow4(Mst2) - shiftst2*pow2(Mst1)*(4*pow2(Mst2)*pow2(Mt)*(
        pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 2*(1 + lmMst1 - lmMst2)*pow2(Mst2)*
        pow2(Sbeta)) - 4*pow2(Mst1)*(2*pow2(Mt)*((-1 + lmMst1 - lmMst2)*pow2(
        MuSUSY)*(-1 + pow2(Sbeta)) + pow2(Mst2)*pow2(Sbeta)) + pow2(s2t)*pow2(
        Sbeta)*pow4(Mst2)) + pow2(s2t)*pow2(Sbeta)*(pow2(Mst2)*pow4(Mst1) + (3
        - 2*lmMst1 + 2*lmMst2)*pow6(Mst2))) + pow2(s2t)*pow2(Sbeta)*pow8(Mst2))
        )))/(12.*Tbeta*pow2(Mst1)*pow2(Sbeta)*pow4(Mst2)) + z2*(pow2(Mt)*pow2(
        s2t)*((pow2(Mst1)*(3*(48*Dmglst1*(1537 + 60*lmMst1 + 228*lmMst2) + (
        40285 + 4032*lmMst1 + 2880*lmMst2)*Mgl) + (2*(72*Dmglst1*(-923 + 894*
        lmMst1 - 1470*lmMst2) + 5*(59 + 3564*lmMst1 - 3564*lmMst2)*Mgl)*pow2(
        Mst1))/pow2(Mst2)))/(486.*Mgl) - (33.94444444444444 - (4*lmMst1)/3. + (
        4*lmMst2)/3. + (44*Dmglst1)/Mgl)*pow2(Mst2) + (pow2(MuSUSY)*(339 - 30*
        lmMsq - 37*lmMst1 + 149*lmMst2 + (-198*Dmglst1 + (pow2(Mst1)*(-30*(2*
        Dmglst1 + Mgl) + (4*Dmglst1*((-390 + 90*lmMsq + 139*lmMst1 - 315*
        lmMst2)*pow2(Msq) + 70*pow2(Mst1)))/pow2(Mst2) + (Mgl*((311 + 1080*
        lmMsq + 2196*lmMst1 - 3348*lmMst2)*pow2(Msq) + 360*pow2(Mst1)))/(9.*
        pow2(Mst2))))/pow2(Msq))/Mgl + (6*pow2(Mst2))/pow2(Mst1) - (25*(4*
        Dmglst1 + Mgl)*pow4(Mst1))/(2.*Mgl*pow4(Msq)) + (9*(71.34670781893004 +
        (40*lmMsq)/3. + (868*lmMst1)/9. - (332*lmMst2)/3. + (Dmglst1*(-1613 +
        1080*lmMsq + 3684*lmMst1 - 5796*lmMst2))/(27.*Mgl))*pow4(Mst1))/pow4(
        Mst2)))/9. + (4*pow4(Mst2))/(3.*pow2(Mst1))) + ((pow2(Mt)*pow2(MuSUSY)*
        (-3*pow2(s2t)*(339 - 30*lmMsq - 37*lmMst1 + 149*lmMst2 + (-198*Dmglst1
        + (pow2(Mst1)*(-30*(2*Dmglst1 + Mgl) + (4*Dmglst1*((-390 + 90*lmMsq +
        139*lmMst1 - 315*lmMst2)*pow2(Msq) + 70*pow2(Mst1)))/pow2(Mst2) + (Mgl*
        ((311 + 1080*lmMsq + 2196*lmMst1 - 3348*lmMst2)*pow2(Msq) + 360*pow2(
        Mst1)))/(9.*pow2(Mst2))))/pow2(Msq))/Mgl + (6*pow2(Mst2))/pow2(Mst1) -
        (25*(4*Dmglst1 + Mgl)*pow4(Mst1))/(2.*Mgl*pow4(Msq)) + (9*(
        71.34670781893004 + (40*lmMsq)/3. + (868*lmMst1)/9. - (332*lmMst2)/3. +
        (Dmglst1*(-1613 + 1080*lmMsq + 3684*lmMst1 - 5796*lmMst2))/(27.*Mgl))*
        pow4(Mst1))/pow4(Mst2)) + (Mt*(4*(12*(55 + 16*lmMst2)*(2*Dmglst1 + Mgl)
        *Mt + (3*Dmglst1*(-687 + 540*lmMsq + 92*lmMst1 - 860*lmMst2) + (2407 +
        180*lmMsq - 408*lmMst1 + 408*lmMst2)*Mgl)*Mst1*s2t)*pow2(Mst1) + (3*
        s2t*(960*Dmglst1*pow2(Msq)*pow5(Mst1) + pow2(Mst2)*(-160*(3*Dmglst1 +
        Mgl)*pow2(Msq)*pow3(Mst1) + 4*(627 - 60*lmMsq - 64*lmMst1 + 192*lmMst2)
        *(Dmglst1 + Mgl)*Mst1*pow4(Msq) - 60*(5*Dmglst1 + Mgl)*pow5(Mst1))))/
        pow4(Msq)))/(Mgl*pow4(Mst2))))/27. + (shiftst3*(-4*pow2(Mst1)*pow2(s2t)
        *pow4(Mst2)*(2*(-1 + lmMst1 - lmMst2)*pow2(Mt)*(-(pow2(MuSUSY)*(-1 +
        pow2(Sbeta))) + pow2(Mst2)*pow2(Sbeta)) + pow2(s2t)*pow2(Sbeta)*pow4(
        Mst2)) + pow2(Mt)*(8*(-1 + 3*lmMst1 - 3*lmMst2)*pow2(MuSUSY)*pow2(s2t)*
        (-1 + pow2(Sbeta))*pow6(Mst1) + 4*(-(pow2(MuSUSY)*pow2(s2t)*(-1 + pow2(
        Sbeta))) + 2*(2*pow2(Mt) - pow2(Mst2)*pow2(s2t))*pow2(Sbeta))*pow6(
        Mst2)) + pow4(Mst1)*(8*(-1 + 2*lmMst1 - 2*lmMst2)*pow2(Mst2)*pow2(Mt)*
        pow2(MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta)) + (3 + 2*lmMst1 - 2*lmMst2)*
        pow2(Sbeta)*pow4(s2t)*pow6(Mst2))))/(6.*pow2(Mst1)*pow4(Mst2)))/pow2(
        Sbeta) - (MuSUSY*(Mt*pow3(s2t)*((33.82716049382716 - (50*lmMsq)/3. - (
        281*lmMst1)/9. + (521*lmMst2)/9. + (Dmglst1*(151.33333333333334 - 40*
        lmMsq - (556*lmMst1)/9. + 140*lmMst2))/Mgl)*pow2(Mst1) + (2*shiftst3*
        pow2(Mst2)*((-1 - lmMst1 + lmMst2)*pow2(Mst1) + pow2(Mst2)))/(3.*pow2(
        Mst1)) - ((10*(34*Dmglst1 + 7*Mgl))/(9.*Mgl*pow2(Msq)) + (
        67.50720164609054 + (208*lmMst1)/3. - (208*lmMst2)/3. + (Dmglst1*(3067
        + 2016*lmMst1 - 2016*lmMst2))/(27.*Mgl))/pow2(Mst2))*pow4(Mst1) - (
        pow2(Mst2)*(666 - 60*lmMsq - 74*lmMst1 + 298*lmMst2 - (396*Dmglst1 + (
        60*(2*Dmglst1 + Mgl)*pow2(Mst1))/pow2(Msq) + (25*(4*Dmglst1 + Mgl)*
        pow4(Mst1))/pow4(Msq))/Mgl))/18.) - (2*(4*s2t*shiftst3*pow2(Mst2)*pow3(
        Mt) + Mt*pow3(s2t)*pow4(Mst2)))/(3.*pow2(Mst1)) + pow2(Mt)*pow2(s2t)*((
        2*Mst1*((627 - 60*lmMsq - 64*lmMst1 + 192*lmMst2)*(Dmglst1 + Mgl) - (
        40*(3*Dmglst1 + Mgl)*pow2(Mst1))/pow2(Msq)))/(3.*Mgl) + (4*(263 + 180*
        lmMsq - 108*lmMst1 - 84*lmMst2 + (3*Dmglst1*(-657 + 300*lmMsq + 78*
        lmMst1 - 526*lmMst2) + (60*(9*Dmglst1 + Mgl)*pow2(Mst1))/pow2(Msq))/
        Mgl)*pow3(Mst1))/(9.*pow2(Mst2)) - ((10*(5*Dmglst1 + Mgl))/(Mgl*pow4(
        Msq)) - (1789 + 720*lmMsq + 420*lmMst1 - 1188*lmMst2 + (Dmglst1*(-6487
        + 3600*lmMsq + 6060*lmMst1 - 11436*lmMst2))/Mgl)/(9.*pow4(Mst2)))*pow5(
        Mst1)) + (pow3(Mt)*(s2t*(72*Dmglst1*((2777 + 120*lmMst1 + 456*lmMst2)*
        pow2(Mst2)*pow4(Mst1) - 297*pow2(Mst1)*pow4(Mst2) + (931 + 1908*lmMst1
        - 2484*lmMst2)*pow6(Mst1)) + Mgl*(6*(17501 + 2124*lmMst1 + 1332*lmMst2)
        *pow2(Mst2)*pow4(Mst1) + 27*(-587 + 24*lmMst1 - 24*lmMst2)*pow2(Mst1)*
        pow4(Mst2) + 4*(26399 + 12096*lmMst1 - 6912*lmMst2)*pow6(Mst1) + 648*
        pow6(Mst2))) + (144*Mt*pow3(Mst1)*(Mgl*(((1381 - 360*lmMsq + 15*lmMst1
        + 795*lmMst2 - 42*lmMt)*pow2(Mst1) - 2*(4 - 48*lmMst1 + 27*lmMst2 + 21*
        lmMt)*pow2(Mst2))*pow4(Msq) - 60*pow2(Msq)*(3*pow4(Mst1) - pow4(Mst2))
        + 15*(3*pow2(Mst1)*pow4(Mst2) + pow6(Mst2))) + Dmglst1*(-(((-4163 +
        1080*lmMsq + 117*lmMst1 - 2421*lmMst2)*pow2(Mst1) + 2*(4 - 48*lmMst1 +
        27*lmMst2 + 21*lmMt)*pow2(Mst2))*pow4(Msq)) + 60*pow2(Msq)*(-15*pow4(
        Mst1) + pow4(Mst2)) + 15*(7*pow2(Mst1)*pow4(Mst2) + pow6(Mst2)))))/
        pow4(Msq)))/(243.*Mgl*pow2(Mst1)*pow4(Mst2))))/Tbeta + (-(((36*Dmglst1*
        (7153 - 1080*lmMsq + 348*lmMst1 + 1764*lmMst2) + (98497 - 16200*lmMsq +
        37044*lmMst1 - 11124*lmMst2)*Mgl)*pow4(Msq)*pow4(Mst1) - 12*pow2(Msq)*
        pow2(Mst2)*((36*Dmglst1*(291 - 90*lmMsq - 139*lmMst1 + 315*lmMst2) + (
        5737 - 1620*lmMsq - 2862*lmMst1 + 6030*lmMst2)*Mgl)*pow2(Msq)*pow2(
        Mst1) - 900*(4*Dmglst1 + Mgl)*pow4(Mst1)) - 54*(4*Dmglst1*(30*pow2(Msq)
        *pow2(Mst1) + 99*pow4(Msq) + 25*pow4(Mst1)) + Mgl*(60*pow2(Msq)*pow2(
        Mst1) + (-654 + 60*lmMsq + 74*lmMst1 - 298*lmMst2)*pow4(Msq) + 25*pow4(
        Mst1)))*pow4(Mst2))*pow4(s2t))/3888. - (Mt*pow3(s2t)*(2*(Dmglst1*(5823
        - 1980*lmMsq - 660*lmMst1 + 3732*lmMst2) + (1355 - 540*lmMsq + 24*
        lmMst1 + 744*lmMst2)*Mgl)*pow2(Mst2)*pow3(Mst1)*pow4(Msq) - pow2(Msq)*(
        (Dmglst1*(1397 + 5124*lmMst1 - 5124*lmMst2) + (737 + 852*lmMst1 - 852*
        lmMst2)*Mgl)*pow2(Msq) + 480*(6*Dmglst1 + Mgl)*pow2(Mst2))*pow5(Mst1) -
        3*pow4(Mst2)*(-80*(3*Dmglst1 + Mgl)*pow2(Msq)*pow3(Mst1) + 2*(627 - 60*
        lmMsq - 64*lmMst1 + 192*lmMst2)*(Dmglst1 + Mgl)*Mst1*pow4(Msq) - 30*(5*
        Dmglst1 + Mgl)*pow5(Mst1))))/(27.*pow2(Mst2)) + (4*((pow4(Mt)*(Mgl*(
        120*pow2(Msq)*pow2(Mst1)*(4*pow2(Mst1) + pow2(Mst2))*pow4(Mst2) + 15*
        pow2(Mst1)*pow4(Mst2)*(6*pow2(Mst1)*pow2(Mst2) + 19*pow4(Mst1) + 3*
        pow4(Mst2)) + pow4(Msq)*(-2*((-289 + 6*lmMst1 - 228*lmMst2 + 222*lmMt)*
        pow2(Mst2) + 6*(55 + 16*lmMst2)*pow2(MuSUSY))*pow4(Mst1) + (392 + 360*
        lmMsq - 111*lmMst1 - 249*lmMt)*pow2(Mst1)*pow4(Mst2) + (-2383 + 900*
        lmMsq + 474*lmMst1 - 1230*lmMst2 - 300*lmMt)*pow6(Mst1) - 18*pow6(Mst2)
        )) + 2*Dmglst1*(360*pow2(Msq)*pow4(Mst1)*pow4(Mst2) + 390*pow4(Mst2)*
        pow6(Mst1) + pow4(Msq)*(4*((151 + 108*lmMst2 - 108*lmMt)*pow2(Mst2) -
        3*(55 + 16*lmMst2)*pow2(MuSUSY))*pow4(Mst1) + 297*pow2(Mst1)*pow4(Mst2)
        + 6*(-933 + 300*lmMsq + 150*lmMst1 - 462*lmMst2 - 40*lmMt)*pow6(Mst1))
        + 90*pow4(Mst1)*pow6(Mst2))))/pow2(Mst1) + Mst1*s2t*pow3(Mt)*(120*Mgl*
        pow2(Msq)*pow2(Mst2)*(pow2(Mst1)*(-2*pow2(Mst2) + pow2(MuSUSY)) - 6*
        pow4(Mst1) + 2*pow4(Mst2)) + Mgl*pow4(Msq)*(3*(-627 + 60*lmMsq + 64*
        lmMst1 - 192*lmMst2)*pow2(Mst2)*pow2(MuSUSY) - pow2(Mst1)*(12*(-463 +
        120*lmMsq + 27*lmMst1 - 283*lmMst2)*pow2(Mst2) + (2407 + 180*lmMsq -
        408*lmMst1 + 408*lmMst2)*pow2(MuSUSY)) + (9178 - 2160*lmMsq - 552*
        lmMst1 + 4392*lmMst2)*pow4(Mst1) + 4*(13 + 96*lmMst1 - 54*lmMst2 - 42*
        lmMt)*pow4(Mst2)) + Dmglst1*pow4(Msq)*(3*(-627 + 60*lmMsq + 64*lmMst1 -
        192*lmMst2)*pow2(Mst2)*pow2(MuSUSY) + pow2(Mst1)*(4*(4150 - 1080*lmMsq
        - 213*lmMst1 + 2475*lmMst2 + 42*lmMt)*pow2(Mst2) + 3*(687 - 540*lmMsq -
        92*lmMst1 + 860*lmMst2)*pow2(MuSUSY)) + 6*(7373 - 1800*lmMsq - 372*
        lmMst1 + 3544*lmMst2 + 28*lmMt)*pow4(Mst1) + 4*(13 + 96*lmMst1 - 54*
        lmMst2 - 42*lmMt)*pow4(Mst2)) + 120*Dmglst1*pow2(Msq)*(-6*(5*pow2(Mst2)
        + pow2(MuSUSY))*pow4(Mst1) + pow2(Mst1)*(3*pow2(Mst2)*pow2(MuSUSY) - 2*
        pow4(Mst2)) + 2*pow6(Mst2)) + 15*Mgl*(3*pow2(Mst2)*(-4*pow2(Mst2) +
        pow2(MuSUSY))*pow4(Mst1) + 8*pow2(Mst1)*pow6(Mst2) + 4*pow8(Mst2)) +
        15*Dmglst1*(pow4(Mst1)*(15*pow2(Mst2)*pow2(MuSUSY) - 28*pow4(Mst2)) +
        24*pow2(Mst1)*pow6(Mst2) + 4*pow8(Mst2)))))/(27.*pow4(Mst2)))/(Mgl*
        pow4(Msq)))) - Al4p*((twoLoopFlag*((-72*Mst1*Mt*pow3(s2t)*(Mgl*(-2*(10
        + 8*lmMst1 - 4*(1 + lmMst1)*lmMst2 + 3*pow2(lmMst1) + pow2(lmMst2))*
        pow2(Mst1)*pow2(Mst2) + (3 + 2*lmMst1 - 2*lmMst2)*pow4(Mst1) + 2*(8 -
        2*lmMst1 + 6*lmMst2 - pow2(lmMst1) + pow2(lmMst2))*pow4(Mst2)) +
        Dmglst1*(-2*(14 + lmMst1*(4 - 12*lmMst2) + 7*pow2(lmMst1) + 5*pow2(
        lmMst2))*pow2(Mst1)*pow2(Mst2) + (19 + 18*lmMst1 - 18*lmMst2)*pow4(
        Mst1) + 2*(12 - 2*lmMst1 + 6*lmMst2 - pow2(lmMst1) + pow2(lmMst2))*
        pow4(Mst2))))/pow2(Mst2) + 9*(-8*Dmglst1*((2 - 2*lmMst1 + lmMst2 -
        pow2(lmMst1) + pow2(lmMst2))*pow2(Mst1)*pow2(Mst2) + (-3 - 2*lmMst1 +
        lmMst2 + pow2(lmMst1) - pow2(lmMst2))*pow4(Mst1) + (-1 + lmMst1)*pow4(
        Mst2)) + Mgl*(-4*(-2*lmMst1*(3 + 2*lmMst2) + pow2(lmMst1) + 3*(-4 -
        lmMst2 + pow2(lmMst2)))*pow2(Mst1)*pow2(Mst2) + (-13 + 14*lmMst2 +
        lmMst1*(-26 + 8*lmMst2) - 12*pow2(lmMst1) + 4*pow2(lmMst2))*pow4(Mst1)
        + 4*(-7 - 7*lmMst2 + 2*lmMst1*(1 + lmMst2) - 2*pow2(lmMst2))*pow4(Mst2)
        ))*pow4(s2t) + (4*Mt*(72*Mst1*s2t*pow2(Mt)*(Mgl*((8 - 2*lmMst1 + 6*
        lmMst2 - pow2(lmMst1) + pow2(lmMst2))*pow2(Mst2)*pow2(MuSUSY) + pow2(
        Mst1)*(2*(lmMst1 - lmMst2 + 2*lmMst1*lmMst2 - 2*lmMst1*lmMt + 2*lmMst2*
        lmMt)*pow2(Mst2) + (6 + 4*lmMst1*(-3 + lmMst2) + 16*lmMst2 - 5*pow2(
        lmMst1))*pow2(MuSUSY) + pow2(lmMst2)*(-4*pow2(Mst2) + pow2(MuSUSY))) +
        2*(4 + lmMst1*(5 + 2*lmMst2 - 2*lmMt) + lmMst2*(-5 + 2*lmMt) - 2*pow2(
        lmMst2))*pow4(Mst1) + (6 - 6*lmMst1 + 2*lmMst2 - 4*lmMst1*lmMt + 4*(1 +
        lmMst2)*lmMt)*pow4(Mst2)) + Dmglst1*((12 - 2*lmMst1 + 6*lmMst2 - pow2(
        lmMst1) + pow2(lmMst2))*pow2(Mst2)*pow2(MuSUSY) + pow2(Mst1)*(2*(-2 +
        3*lmMst1 - 3*lmMst2)*(1 + 2*lmMst2 - 2*lmMt)*pow2(Mst2) + (10 - 8*
        lmMst1 + 12*(1 + lmMst1)*lmMst2 - 9*pow2(lmMst1) - 3*pow2(lmMst2))*
        pow2(MuSUSY)) + (20 + 10*lmMst1*(5 + 2*lmMst2 - 2*lmMt) + 8*lmMt +
        lmMst2*(-58 + 20*lmMt) - 20*pow2(lmMst2))*pow4(Mst1) + 2*(8 - 5*lmMst1
        + lmMst2 + (4 - 2*lmMst1 + 2*lmMst2)*lmMt)*pow4(Mst2))) + 8*pow3(Mt)*(
        2*Dmglst1*(-18*(1 + lmMt)*pow2(Mst1)*pow2(Mst2) + 18*(1 + lmMst2*(5 -
        2*lmMt) - 2*lmMst1*(1 + lmMst2 - lmMt) - 3*lmMt + 2*pow2(lmMst2))*pow4(
        Mst1) - (1 + 12*lmMst1 + 6*lmMt)*pow4(Mst2)) - 9*Mgl*(2*(1 + lmMt)*
        pow2(Mst1)*pow2(Mst2) + 2*(lmMst1*(1 + lmMst2 - lmMt) + lmMst2*(-2 +
        lmMt) + lmMt)*pow4(Mst1) - (3 + lmMst1 - lmMt + 2*lmMst1*lmMt + 2*
        lmMst2*(1 + lmMt) + pow2(lmMst1) - 6*pow2(lmMt))*pow4(Mst2) - pow2(
        lmMst2)*(2*pow4(Mst1) + pow4(Mst2)))) - 9*MuSUSY*((Mt*MuSUSY*s2t*(8*
        Dmglst1*Mt*(Mst1*(12 - 2*lmMst1 + 6*lmMst2 - pow2(lmMst1) + pow2(
        lmMst2))*pow2(Mst2) + (10 + 12*lmMst2 + 4*lmMst1*(-2 + 3*lmMst2) - 9*
        pow2(lmMst1) - 3*pow2(lmMst2))*pow3(Mst1)) + 8*Mgl*Mt*(Mst1*(8 - 2*
        lmMst1 + 6*lmMst2 - pow2(lmMst1) + pow2(lmMst2))*pow2(Mst2) + (6 + 4*
        lmMst1*(-3 + lmMst2) + 16*lmMst2 - 5*pow2(lmMst1) + pow2(lmMst2))*pow3(
        Mst1)) - 8*Dmglst1*s2t*((-(lmMst2*(1 + lmMst2)) + pow2(lmMst1))*pow2(
        Mst1)*pow2(Mst2) + (2 + 3*lmMst1 - 3*lmMst2 + pow2(lmMst1) - pow2(
        lmMst2))*pow4(Mst1) - (-1 + lmMst1)*pow4(Mst2)) + Mgl*s2t*(4*(2 + 11*
        lmMst2 - 2*lmMst1*(5 + 4*lmMst2) + pow2(lmMst1) + 7*pow2(lmMst2))*pow2(
        Mst1)*pow2(Mst2) + (1 + 46*lmMst2 - 2*lmMst1*(23 + 32*lmMst2) + 20*
        pow2(lmMst1) + 44*pow2(lmMst2))*pow4(Mst1) + 4*(7 + 7*lmMst2 - 2*
        lmMst1*(1 + lmMst2) + 2*pow2(lmMst2))*pow4(Mst2))))/pow2(Sbeta) + (
        pow2(Mst2)*pow3(s2t)*(8*Dmglst1*((-1 + lmMst1 - lmMst2 + pow2(lmMst1) -
        pow2(lmMst2))*pow2(Mst1)*pow2(Mst2) + (2 + 3*lmMst1 - 2*lmMst2)*pow4(
        Mst1) - (-1 + lmMst1)*pow4(Mst2)) + Mgl*(-4*(-5 + 4*lmMst2 - 2*lmMst1*(
        4 + 3*lmMst2) + pow2(lmMst1) + 5*pow2(lmMst2))*pow2(Mst1)*pow2(Mst2) +
        (7 - 2*lmMst2 + lmMst1*(6 + 32*lmMst2) - 16*(pow2(lmMst1) + pow2(
        lmMst2)))*pow4(Mst1) - 4*((7 - 2*lmMst1)*(1 + lmMst2) + 2*pow2(lmMst2))
        *pow4(Mst2))) - 8*s2t*pow2(Mt)*(Dmglst1*(4*(1 + lmMst2)*pow2(Mst1)*
        pow2(Mst2) + (4 - 8*lmMst1 + 8*lmMst2)*pow4(Mst1) + 2*(1 - 2*lmMst1)*
        pow4(Mst2)) + Mgl*(2*pow2(Mst1)*pow2(Mst2) + 2*lmMst2*(pow2(Mst1)*pow2(
        Mst2) + pow4(Mst1) - 2*pow4(Mst2)) + pow2(lmMst1)*(-2*pow2(Mst1)*pow2(
        Mst2) - 2*pow4(Mst1) + pow4(Mst2)) - pow2(lmMst2)*(2*pow2(Mst1)*pow2(
        Mst2) + 2*pow4(Mst1) + 3*pow4(Mst2)) + 2*lmMst1*(2*lmMst2*pow2(Mst1)*
        pow2(Mst2) + (-1 + 2*lmMst2)*pow4(Mst1) + (1 + lmMst2)*pow4(Mst2)))) +
        Mst1*(32*((Dmglst1*(5 - lmMst2 + lmMst1*(-3 + 3*lmMst2 - 4*lmMt) + 4*
        lmMt + 4*lmMst2*lmMt - 3*pow2(lmMst2)) + Mgl*(2 + lmMst2 + lmMst1*(-2 +
        lmMst2 - 2*lmMt) + lmMt + 2*lmMst2*lmMt - pow2(lmMst2)))*pow2(Mst1) + (
        Dmglst1*(5 + lmMst2 + (2 + lmMst2)*lmMt - lmMst1*(3 + lmMt)) + (2 +
        lmMst2 + lmMt + lmMst2*lmMt - lmMst1*(2 + lmMt))*Mgl)*pow2(Mst2))*pow3(
        Mt) - 6*Mt*pow2(s2t)*(-(Mgl*(4*(1 + 5*lmMst1 - 5*lmMst2 - 2*lmMst1*
        lmMst2 + 2*pow2(lmMst1))*pow2(Mst1)*pow2(Mst2) + (1 + 18*lmMst1 - 18*
        lmMst2 - 8*lmMst1*lmMst2 + 8*pow2(lmMst1))*pow4(Mst1) + 2*(-8 + 2*
        lmMst1 - 6*lmMst2 + pow2(lmMst1) - pow2(lmMst2))*pow4(Mst2))) +
        Dmglst1*(-4*(1 + 3*lmMst1 - 3*lmMst2 - 6*lmMst1*lmMst2 + 4*pow2(lmMst1)
        + 2*pow2(lmMst2))*pow2(Mst1)*pow2(Mst2) + (15 + 6*lmMst1 - 6*lmMst2 +
        24*lmMst1*lmMst2 - 16*pow2(lmMst1) - 8*pow2(lmMst2))*pow4(Mst1) + 2*(12
        - 2*lmMst1 + 6*lmMst2 - pow2(lmMst1) + pow2(lmMst2))*pow4(Mst2)))))/
        Tbeta) + 9*Mt*pow2(s2t)*(Mgl*((8*(1 + lmMst1)*pow2(Mst2) + (1 + 46*
        lmMst2 - 2*lmMst1*(23 + 32*lmMst2) + 20*pow2(lmMst1) + 44*pow2(lmMst2))
        *pow2(MuSUSY))*pow4(Mst1) - 4*(((-(lmMst2*(4 + 3*lmMst2)) + pow2(
        lmMst1))*pow2(Mst2) - (7 + 7*lmMst2 + 2*pow2(lmMst2))*pow2(MuSUSY) + 2*
        lmMst1*(1 + lmMst2)*(pow2(Mst2) + pow2(MuSUSY)))*pow4(Mst2) + pow2(
        Mst1)*(-((2 + 11*lmMst2 - 2*lmMst1*(5 + 4*lmMst2) + pow2(lmMst1) + 7*
        pow2(lmMst2))*pow2(Mst2)*pow2(MuSUSY)) + (2 + 2*lmMst1*(-1 + lmMst2) +
        6*lmMst2 - 3*pow2(lmMst1) + pow2(lmMst2))*pow4(Mst2)))) - 8*Dmglst1*(
        pow2(MuSUSY)*((pow2(lmMst1) - pow2(lmMst2))*pow2(Mst1)*(pow2(Mst1) +
        pow2(Mst2)) + 2*pow4(Mst1) + pow4(Mst2)) - lmMst1*((4*pow2(Mst2) - 3*
        pow2(MuSUSY))*pow4(Mst1) - 2*pow2(Mst1)*pow4(Mst2) + (2*pow2(Mst2) +
        pow2(MuSUSY))*pow4(Mst2)) + pow2(Mst1)*(pow4(Mst2) + lmMst2*(pow2(Mst1)
        *(2*pow2(Mst2) - 3*pow2(MuSUSY)) - pow2(Mst2)*pow2(MuSUSY) + 2*pow4(
        Mst2))) + pow6(Mst2)))))/pow4(Mst2)))/(216.*Mgl) + ((z2*(6*Mgl*
        twoLoopFlag*pow4(Msq)*(Tbeta*pow2(Sbeta)*(pow3(Mt)*(16*s2t*((Dmglst1 +
        Mgl)*Mst1*pow2(Mst2)*pow2(MuSUSY) - (4*(3*Dmglst1 + Mgl)*pow2(Mst2) + (
        3*Dmglst1 - Mgl)*pow2(MuSUSY))*pow3(Mst1)) + 32*((4*Dmglst1 + Mgl)*Mt -
        2*(5*Dmglst1 + Mgl)*Mst1*s2t)*pow4(Mst1)) - 4*pow2(Mt)*pow2(MuSUSY)*
        pow2(s2t)*(-4*Dmglst1*pow2(Mst1)*(pow2(Mst1) + pow2(Mst2)) + Mgl*pow4(
        Mst2)) + pow3(s2t)*pow4(Mst2)*(-8*(Dmglst1 + Mgl)*Mst1*Mt*pow2(Mst2) -
        2*(2*Dmglst1 + Mgl)*s2t*pow2(Mst1)*pow2(Mst2) + 8*(5*Dmglst1 + Mgl)*Mt*
        pow3(Mst1) + (4*Dmglst1 + Mgl)*s2t*pow4(Mst1) + Mgl*s2t*pow4(Mst2))) +
        4*Mt*MuSUSY*(Mt*MuSUSY*s2t*Tbeta*(-4*(Dmglst1 + Mgl)*Mst1*Mt*pow2(Mst2)
        - 4*Dmglst1*s2t*pow2(Mst1)*(pow2(Mst1) + pow2(Mst2)) + 4*(3*Dmglst1 -
        Mgl)*Mt*pow3(Mst1) + Mgl*s2t*pow4(Mst2)) + pow2(Sbeta)*(16*(3*Dmglst1 +
        Mgl)*pow3(Mst1)*pow3(Mt) + ((4*Dmglst1 + Mgl)*pow2(Mst1) - Mgl*pow2(
        Mst2))*pow3(s2t)*pow4(Mst2) + 6*Mt*pow2(s2t)*(-4*Dmglst1*(pow2(Mst1) +
        pow2(Mst2))*pow3(Mst1) + (Dmglst1 + Mgl)*Mst1*pow4(Mst2))))) +
        xDmglst1*pow2(Dmglst1)*(12*twoLoopFlag*pow4(Msq)*(12*Mt*MuSUSY*pow2(
        Sbeta)*(16*pow3(Mst1)*pow3(Mt) + pow2(Mst1)*pow3(s2t)*pow4(Mst2) + Mt*
        pow2(s2t)*(-10*(pow2(Mst1) + pow2(Mst2))*pow3(Mst1) + Mst1*pow4(Mst2)))
        + Tbeta*(4*Mst1*s2t*pow2(Mt)*((2*Mt*(9*pow2(Mst1) - pow2(Mst2)) - 3*
        Mst1*s2t*(pow2(Mst1) + pow2(Mst2)))*pow2(MuSUSY) - 6*Mt*pow2(Mst1)*(8*
        pow2(Mst2) + 3*pow2(MuSUSY))*pow2(Sbeta)) + pow2(Sbeta)*(160*(Mt - 3*
        Mst1*s2t)*pow3(Mt)*pow4(Mst1) + 4*s2t*pow2(Mt)*pow2(MuSUSY)*(Mst1*(2*Mt
        + 3*Mst1*s2t)*pow2(Mst2) + 3*s2t*pow4(Mst1)) + Mst1*(44*Mt*pow2(Mst1) -
        4*Mt*pow2(Mst2) - 3*Mst1*s2t*pow2(Mst2) + 3*s2t*pow3(Mst1))*pow3(s2t)*
        pow4(Mst2)))) + (Al4p*threeLoopFlag*(Tbeta*pow2(Sbeta)*(4*pow2(Mt)*
        pow2(s2t)*(pow4(Msq)*(-12*pow2(Mst1)*pow2(Mst2)*(2*(1537 + 60*lmMst1 +
        228*lmMst2)*pow2(Mst2) + 3*(-390 + 90*lmMsq + 139*lmMst1 - 315*lmMst2)*
        pow2(MuSUSY)) + (8*(7429 - 4506*lmMst1 + 8538*lmMst2)*pow2(Mst2) - (
        2903 + 3240*lmMsq + 10956*lmMst1 - 17292*lmMst2)*pow2(MuSUSY))*pow4(
        Mst1) + 1782*(2*pow2(Mst2) + pow2(MuSUSY))*pow4(Mst2)) + pow2(MuSUSY)*(
        750*pow4(Mst1)*pow4(Mst2) - 60*pow2(Msq)*(82*pow2(Mst2)*pow4(Mst1) - 9*
        pow2(Mst1)*pow4(Mst2)))) - pow4(Mst2)*(750*pow4(Mst1)*pow4(Mst2) +
        pow4(Msq)*(-36*(-291 + 90*lmMsq + 139*lmMst1 - 315*lmMst2)*pow2(Mst1)*
        pow2(Mst2) + (-29201 + 3240*lmMsq - 948*lmMst1 - 5388*lmMst2)*pow4(
        Mst1) + 1782*pow4(Mst2)) - 60*pow2(Msq)*(100*pow2(Mst2)*pow4(Mst1) - 9*
        pow2(Mst1)*pow4(Mst2)))*pow4(s2t) + 72*Mt*pow2(Mst2)*pow3(s2t)*(150*
        pow4(Mst2)*pow5(Mst1) - (pow4(Msq)*(2*(-4315 + 1380*lmMsq + 562*lmMst1
        - 2738*lmMst2)*pow2(Mst2)*pow3(Mst1) + 2*(627 - 60*lmMsq - 64*lmMst1 +
        192*lmMst2)*Mst1*pow4(Mst2) + (445 + 4516*lmMst1 - 4516*lmMst2)*pow5(
        Mst1)))/3. - 80*pow2(Msq)*(-2*pow3(Mst1)*pow4(Mst2) + 13*pow2(Mst2)*
        pow5(Mst1))) - 32*Mst1*s2t*pow3(Mt)*(pow4(Msq)*(3*(-627 + 60*lmMsq +
        64*lmMst1 - 192*lmMst2)*pow2(Mst2)*pow2(MuSUSY) + pow2(Mst1)*((33166 -
        8640*lmMsq - 1644*lmMst1 + 19656*lmMst2 + 420*lmMt)*pow2(Mst2) + 3*(
        3061 - 1260*lmMsq - 434*lmMst1 + 2354*lmMst2)*pow2(MuSUSY)) + (129748 -
        32400*lmMsq - 5232*lmMst1 + 61932*lmMst2 + 900*lmMt)*pow4(Mst1) + 4*(13
        + 96*lmMst1 - 54*lmMst2 - 42*lmMt)*pow4(Mst2)) - 120*pow2(Msq)*(9*(10*
        pow2(Mst2) + 3*pow2(MuSUSY))*pow4(Mst1) + 2*pow2(Mst1)*(-3*pow2(Mst2)*
        pow2(MuSUSY) + pow4(Mst2)) - 2*pow6(Mst2)) + 15*pow2(Mst2)*((-52*pow2(
        Mst2) + 45*pow2(MuSUSY))*pow4(Mst1) + 48*pow2(Mst1)*pow4(Mst2) + 4*
        pow6(Mst2))) - 24*pow4(Mt)*(1440*pow2(Msq)*pow2(Mst1)*pow4(Mst2) + 4*
        pow4(Msq)*(4*pow2(Mst1)*((151 + 108*lmMst2 - 108*lmMt)*pow2(Mst2) - 3*(
        55 + 16*lmMst2)*pow2(MuSUSY)) + 4*(-2517 + 750*lmMsq + 361*lmMst1 -
        1198*lmMst2 - 43*lmMt)*pow4(Mst1) + 297*pow4(Mst2)) + 40*(59*pow4(Mst1)
        *pow4(Mst2) + 9*pow2(Mst1)*pow6(Mst2)))) + 4*Mt*MuSUSY*(Mt*MuSUSY*
        Tbeta*(pow4(Msq)*(24*Mst1*Mt*s2t*((3061 - 1260*lmMsq - 434*lmMst1 +
        2354*lmMst2)*pow2(Mst1) + (-627 + 60*lmMsq + 64*lmMst1 - 192*lmMst2)*
        pow2(Mst2)) - 36*pow2(Mst1)*(8*(55 + 16*lmMst2)*pow2(Mt) + (390 - 90*
        lmMsq - 139*lmMst1 + 315*lmMst2)*pow2(Mst2)*pow2(s2t)) + pow2(s2t)*((
        2903 + 3240*lmMsq + 10956*lmMst1 - 17292*lmMst2)*pow4(Mst1) - 1782*
        pow4(Mst2))) + s2t*(150*pow2(Mst2)*(36*Mst1*Mt - 5*s2t*pow2(Mst2))*
        pow4(Mst1) - 60*pow2(Msq)*pow2(Mst1)*(-96*Mst1*Mt*pow2(Mst2) - 82*s2t*
        pow2(Mst1)*pow2(Mst2) + 432*Mt*pow3(Mst1) + 9*s2t*pow4(Mst2)))) + pow2(
        Sbeta)*(-(pow4(Msq)*(pow2(Mst2)*(4*Mst1*Mt*(16*(4 - 48*lmMst1 + 27*
        lmMst2 + 21*lmMt)*pow2(Mt) + 9*(-627 + 60*lmMsq + 64*lmMst1 - 192*
        lmMst2)*pow2(Mst2)*pow2(s2t)) - 6*s2t*pow2(Mst1)*(4*(2777 + 120*lmMst1
        + 456*lmMst2)*pow2(Mt) + 3*(681 - 180*lmMsq - 278*lmMst1 + 630*lmMst2)*
        pow2(Mst2)*pow2(s2t))) + s2t*(8*(6527 - 9372*lmMst1 + 15708*lmMst2)*
        pow2(Mt) + (16943 + 5952*lmMst1 - 5952*lmMst2)*pow2(Mst2)*pow2(s2t))*
        pow4(Mst1) + Mt*(-8*(4*(8336 - 2160*lmMsq - 315*lmMst1 + 4860*lmMst2 +
        63*lmMt)*pow2(Mt) + 9*(-1844 + 660*lmMsq + 249*lmMst1 - 1273*lmMst2)*
        pow2(Mst2)*pow2(s2t))*pow3(Mst1) + 18*(6931 - 2640*lmMsq - 5512*lmMst1
        + 9608*lmMst2)*pow2(s2t)*pow5(Mst1)) - 1782*(-4*s2t*pow2(Mt)*pow4(Mst2)
        + pow3(s2t)*pow6(Mst2)))) + Mst1*(30*(16*(13*pow2(Mst1) + pow2(Mst2))*
        pow3(Mt) + 25*pow2(Mst2)*pow3(Mst1)*pow3(s2t) - 270*Mt*pow2(s2t)*pow4(
        Mst1))*pow4(Mst2) + 60*pow2(Msq)*(-72*(-11*Mt*pow2(Mst2)*pow2(s2t) +
        20*pow3(Mt))*pow4(Mst1) + (-144*Mt*pow2(Mst1)*pow2(s2t) + 32*pow3(Mt) -
        91*pow3(Mst1)*pow3(s2t))*pow4(Mst2) + 9*Mst1*pow3(s2t)*pow6(Mst2)))))))
        /6.)))/(36.*Tbeta*pow2(Sbeta)*pow4(Msq)*pow4(Mst2)) + xDmglst1*pow2(
        Dmglst1)*((twoLoopFlag*(9*(-3*(4 - 2*lmMst1 + lmMst2 - pow2(lmMst1) +
        pow2(lmMst2))*pow2(Mst1)*pow2(Mst2) + 3*(4 + 4*lmMst1 - 3*lmMst2 -
        pow2(lmMst1) + pow2(lmMst2))*pow4(Mst1) + (4 - 3*lmMst1)*pow4(Mst2))*
        pow4(s2t) - (18*Mt*pow3(s2t)*(-2*(30 + 14*lmMst2 - 2*lmMst1*(5 + 12*
        lmMst2) + 13*pow2(lmMst1) + 11*pow2(lmMst2))*pow2(Mst2)*pow3(Mst1) + 2*
        Mst1*(16 - 2*lmMst1 + 6*lmMst2 - pow2(lmMst1) + pow2(lmMst2))*pow4(
        Mst2) + (43 + 66*lmMst1 - 66*lmMst2)*pow5(Mst1)))/pow2(Mst2) + (2*Mt*(
        8*pow3(Mt)*(-27*(1 + lmMt)*pow2(Mst1)*pow2(Mst2) + 9*(9 + lmMst2*(29 -
        10*lmMt) - 10*lmMst1*(1 + lmMst2 - lmMt) - 19*lmMt + 10*pow2(lmMst2))*
        pow4(Mst1) + (11 - 21*lmMst1 - 6*lmMt)*pow4(Mst2)) + 2*Mt*s2t*(2*Mst1*
        Mt*(9*(16 - 2*lmMst1 + 6*lmMst2 - pow2(lmMst1) + pow2(lmMst2))*pow2(
        Mst2)*pow2(MuSUSY) - 9*pow2(Mst1)*(2*(7 - 6*lmMst1 + 6*lmMst2)*(1 + 2*
        lmMst2 - 2*lmMt)*pow2(Mst2) + (-2 + 2*lmMst2 - 6*lmMst1*(1 + 4*lmMst2)
        + 15*pow2(lmMst1) + 9*pow2(lmMst2))*pow2(MuSUSY)) + 18*(5 + 15*lmMst1*(
        5 + 2*lmMst2 - 2*lmMt) + 22*lmMt + lmMst2*(-97 + 30*lmMt) - 30*pow2(
        lmMst2))*pow4(Mst1) + (251 + 102*lmMt + 18*lmMst2*(1 + 2*lmMt) - 12*
        lmMst1*(10 + 3*lmMt))*pow4(Mst2)) + (9*pow2(MuSUSY)*(2*Mst1*Mt*(-16 +
        2*lmMst1 - 6*lmMst2 + pow2(lmMst1) - pow2(lmMst2))*pow2(Mst2) + s2t*(-4
        - 3*lmMst2 + 3*pow2(lmMst1) - 3*pow2(lmMst2))*pow2(Mst1)*pow2(Mst2) +
        2*Mt*(-2 + 2*lmMst2 - 6*lmMst1*(1 + 4*lmMst2) + 15*pow2(lmMst1) + 9*
        pow2(lmMst2))*pow3(Mst1) + 3*s2t*(5*lmMst1 - lmMst2*(5 + lmMst2) +
        pow2(lmMst1))*pow4(Mst1) + (4 - 3*lmMst1)*s2t*pow4(Mst2)))/pow2(Sbeta))
        - (MuSUSY*(8*(151 + 51*lmMt + 18*lmMst2*(1 + lmMt) - 3*lmMst1*(23 + 6*
        lmMt))*Mst1*pow2(Mst2)*pow3(Mt) + 4*pow3(Mst1)*(27*Mt*(7 + 4*lmMst2 -
        4*lmMst1*(1 + 3*lmMst2) + 7*pow2(lmMst1) + 5*pow2(lmMst2))*pow2(Mst2)*
        pow2(s2t) + 2*(151 + 3*lmMst1*(-23 + 36*lmMst2 - 42*lmMt) + 177*lmMt +
        18*lmMst2*(-6 + 7*lmMt) - 108*pow2(lmMst2))*pow3(Mt)) + 18*(4*(-9 + 10*
        lmMst1 - 10*lmMst2)*s2t*pow2(Mt) + (4 + 15*lmMst1 - 12*lmMst2)*pow2(
        Mst2)*pow3(s2t))*pow4(Mst1) + 54*Mst1*Mt*(-16 + 2*lmMst1 - 6*lmMst2 +
        pow2(lmMst1) - pow2(lmMst2))*pow2(s2t)*pow4(Mst2) + 6*s2t*((-34 + 36*
        lmMst1)*pow2(Mt) + 3*(4 - 3*lmMst1)*pow2(Mst2)*pow2(s2t))*pow4(Mst2) -
        18*pow2(Mst1)*(12*(1 + lmMst2)*s2t*pow2(Mst2)*pow2(Mt) + (8 - 3*lmMst1
        + 3*lmMst2 - 3*pow2(lmMst1) + 3*pow2(lmMst2))*pow3(s2t)*pow4(Mst2)) +
        27*Mt*(-15 + 82*lmMst2 - 2*lmMst1*(41 + 24*lmMst2) + 28*pow2(lmMst1) +
        20*pow2(lmMst2))*pow2(s2t)*pow5(Mst1)))/Tbeta - 6*Mt*pow2(s2t)*(((36 -
        60*lmMst1 + 42*lmMst2)*pow2(Mst2) + 9*(5*lmMst1 - lmMst2*(5 + lmMst2) +
        pow2(lmMst1))*pow2(MuSUSY))*pow4(Mst1) + 3*(4 - 3*lmMst1)*pow2(MuSUSY)*
        pow4(Mst2) + pow2(Mst1)*(-3*(4 + 3*lmMst2*(1 + lmMst2) - 3*pow2(lmMst1)
        )*pow2(Mst2)*pow2(MuSUSY) + (1 + 18*(lmMst1 + lmMst2))*pow4(Mst2)) + (
        17 - 18*lmMst1)*pow6(Mst2))))/pow4(Mst2)))/54. + (Al4p*threeLoopFlag*(
        100*Mt*MuSUSY*pow2(Sbeta)*(-(pow2(Mst2)*pow3(s2t)*(pow4(Msq)*(72*pow2(
        Mst1)*pow2(Mst2)*(709 - 432*B4 - 360*DN + 15720*lmMsq + 2*lmMst1*(-6599
        - 3240*lmMsq + 1620*pow2(lmMsq)) - 18*lmMst2*(689 - 360*lmMsq + 62*
        lmMst1 + 180*pow2(lmMsq) - 158*pow2(lmMst1)) - 162*(-53 + 20*lmMsq)*
        pow2(lmMst1) + 18*(-735 + 180*lmMsq + 98*lmMst1)*pow2(lmMst2) + 780*
        pow3(lmMst1) - 5388*pow3(lmMst2)) + (283217 + 55296*B4 + 27648*DN +
        142560*lmMsq - 60*(-67111 + 32400*lmMsq)*lmMst1 + 116640*pow2(lmMsq) +
        12*lmMst2*(-289115 + 142560*lmMsq + 174900*lmMst1 - 101376*pow2(lmMst1)
        ) + 1003464*pow2(lmMst1) + 72*(-32431 + 26112*lmMst1)*pow2(lmMst2) +
        184320*pow3(lmMst1) - 847872*pow3(lmMst2))*pow4(Mst1) - 72*(6457 + 60*
        lmMsq*(67 - 54*lmMst1) - 982*lmMst1 + 30*(-89 + 66*lmMst1)*lmMst2 +
        1620*pow2(lmMsq) + 834*pow2(lmMst1) + 162*pow2(lmMst2))*pow4(Mst2)) +
        1440*pow2(Msq)*((994 - 593*lmMst1 + 6*lmMsq*(-18 + 91*lmMst1 - 91*
        lmMst2) + 701*lmMst2 - 273*pow2(lmMst1) + 273*pow2(lmMst2))*pow2(Mst2)*
        pow4(Mst1) - 3*(108 - 133*lmMst1 + 2*lmMsq*(32 + 9*lmMst1 - 9*lmMst2) +
        69*lmMst2 - 9*pow2(lmMst1) + 9*pow2(lmMst2))*pow2(Mst1)*pow4(Mst2)) -
        180*((3141 - 2306*lmMst1 + 12*lmMsq*(47 + 50*lmMst1 - 50*lmMst2) +
        1742*lmMst2 - 300*pow2(lmMst1) + 300*pow2(lmMst2))*pow4(Mst1)*pow4(
        Mst2) - 54*pow2(Mst1)*pow6(Mst2)))) + 16*s2t*pow2(Mt)*(2*pow4(Msq)*(3*
        pow2(Mst1)*pow2(Mst2)*(248545 - 864*B4 - 432*DN - 14580*lmMsq + 9042*
        lmMst1 + 4860*pow2(lmMsq) - 6*lmMst2*(-22907 + 1620*lmMsq + 810*lmMst1
        - 432*pow2(lmMst1)) + 3870*pow2(lmMst1) + 126*(289 + 48*lmMst1)*pow2(
        lmMst2) + 576*lmMt*(29 + 35*lmMst2 - lmMst1*(35 + 18*lmMst2) + 9*pow2(
        lmMst1) + 9*pow2(lmMst2)) - 7200*pow3(lmMst1) - 1440*pow3(lmMst2)) + 2*
        (333689 - 1296*B4 - 648*DN - 105300*lmMsq + 24*(-13409 + 3375*lmMsq)*
        lmMst1 - 346779*pow2(lmMst1) + lmMst2*(498072 - 81000*lmMsq + 610614*
        lmMst1 + 17280*pow2(lmMst1)) + 9*(-26819 + 16800*lmMst1)*pow2(lmMst2) +
        864*lmMt*(95 + 178*lmMst2 - 2*lmMst1*(89 + 69*lmMst2) + 69*pow2(lmMst1)
        + 69*pow2(lmMst2)) - 77472*pow3(lmMst1) - 91008*pow3(lmMst2))*pow4(
        Mst1) - 18*(4401 + 60*lmMsq*(26 - 27*lmMst1) - 1462*lmMst1 + 6*(-74 +
        69*lmMst1)*lmMst2 + 810*pow2(lmMsq) + 993*pow2(lmMst1) + 81*pow2(
        lmMst2))*pow4(Mst2)) - 360*pow2(Msq)*(50*(11 - 6*lmMsq + 6*lmMst2)*
        pow2(Mst2)*pow4(Mst1) + 3*(-9 + 100*lmMsq - 118*lmMst1 + 18*lmMst2)*
        pow2(Mst1)*pow4(Mst2)) - 45*(2*(307 + 936*lmMsq - 1236*lmMst1 + 300*
        lmMst2)*pow4(Mst1)*pow4(Mst2) + 9*(13 - 12*lmMsq + 12*lmMst2)*pow2(
        Mst1)*pow6(Mst2))) + Mst1*(32*pow3(Mt)*(4*(3*pow2(Mst2)*(33943 - 2918*
        lmMst2 + 60*lmMsq*(244 + 69*lmMt + 9*lmMst2*(3 + 2*lmMt) - 6*lmMst1*(16
        + 3*lmMt)) + 540*(lmMst1 - lmMst2)*pow2(lmMsq) + 6*(799 + 30*lmMst2 +
        84*lmMt)*pow2(lmMst1) - 612*pow2(lmMst2) - 6*lmMt*(963 + 617*lmMst2 +
        108*pow2(lmMst2)) - 144*(17 + 6*lmMst2)*pow2(lmMt) + 2*lmMst1*(-10100 -
        1173*lmMt + lmMst2*(429 + 72*lmMt) + 18*pow2(lmMst2) + 432*pow2(lmMt))
        - 132*pow3(lmMst1) - 84*pow3(lmMst2)) + pow2(Mst1)*(112558 + 118485*
        lmMst2 + 1620*(lmMst1 - lmMst2)*pow2(lmMsq) + (19044 - 19116*lmMst2 +
        45306*lmMt)*pow2(lmMst1) + 180*lmMsq*(244 + 6*lmMst1*(-16 + 36*lmMst2 -
        39*lmMt) + 429*lmMt + 9*lmMst2*(-37 + 26*lmMt) - 216*pow2(lmMst2)) +
        111204*pow2(lmMst2) - 18*lmMt*(6121 + 4478*lmMst2 + 555*pow2(lmMst2)) -
        9*lmMst1*(6907 + 2084*lmMt + 12*lmMst2*(50 + 327*lmMt) + 234*pow2(
        lmMst2) - 2016*pow2(lmMt)) - 432*(59 + 42*lmMst2)*pow2(lmMt) - 9018*
        pow3(lmMst1) + 30240*pow3(lmMst2)))*pow4(Msq) + 240*pow2(Msq)*((1174 +
        6*lmMst1*(-137 + 39*lmMsq - 39*lmMt) + 429*lmMt + lmMst2*(393 - 234*
        lmMsq + 234*lmMt))*pow2(Mst1)*pow2(Mst2) + (1174 - 3*lmMst2*(157 + 78*
        lmMsq - 348*lmMt) + 6*lmMst1*(-137 + 39*lmMsq + 135*lmMst2 - 174*lmMt)
        + 1293*lmMt - 810*pow2(lmMst2))*pow4(Mst1) - (17 + 3*lmMst2*(1 - 6*
        lmMt) - 33*lmMt + 6*lmMsq*(5 + 3*(lmMst2 + lmMt)) - 18*pow2(lmMsq))*
        pow4(Mst2)) + 15*(4*(3220 + 6*lmMst1*(-323 + 174*lmMsq - 174*lmMt) +
        1293*lmMt + lmMst2*(645 - 1044*lmMsq + 1044*lmMt))*pow2(Mst2)*pow4(
        Mst1) - (2519 + lmMst2*(930 - 936*lmMt) - 1014*lmMt + 12*lmMsq*(7 + 78*
        (lmMst2 + lmMt)) - 936*pow2(lmMsq))*pow2(Mst1)*pow4(Mst2) + (-299 + 78*
        lmMt - 12*lmMsq*(-5 + 6*lmMst2 + 6*lmMt) + 6*lmMst2*(-23 + 12*lmMt) +
        72*pow2(lmMsq))*pow6(Mst2))) - 72*Mt*pow2(s2t)*(-720*pow2(Msq)*pow2(
        Mst2)*(2*(59 + 8*(-2 + 3*lmMsq)*lmMst2 - 2*lmMst1*(-8 + 12*lmMsq + 45*
        lmMst2) + 57*pow2(lmMst1) + 33*pow2(lmMst2))*pow4(Mst1) + pow4(Mst2)) +
        pow4(Mst1)*(108*(8261.658179012345 + (1048*B4)/3. - (20*DN)/3. - 610*
        lmMsq + (5*lmMst1*(29755 + 67824*lmMsq + 864*pow2(lmMsq)))/108. -
        lmMst2*(886.8796296296297 + lmMsq*(3140 - 960*lmMst1) + (2033*lmMst1)/
        3. + 40*pow2(lmMsq) - (508*pow2(lmMst1))/3.) - (1665.1666666666667 +
        520*lmMsq)*pow2(lmMst1) - ((-14057 + 2640*lmMsq + 6136*lmMst1)*pow2(
        lmMst2))/6. + (3268*pow3(lmMst1))/9. + (4412*pow3(lmMst2))/9.)*pow4(
        Msq) + 5*(17171 - 2610*lmMst1 + 12*lmMsq*(-541 + 270*lmMst1 - 270*
        lmMst2) + 9102*lmMst2 - 1620*pow2(lmMst1) + 1620*pow2(lmMst2))*pow4(
        Mst2)) + pow4(Msq)*(4*pow2(Mst1)*pow2(Mst2)*(29170 + 9432*B4 - 32400*
        lmMsq + 3*(lmMst1*(1349 + 8280*lmMsq + 360*pow2(lmMsq)) + lmMst2*(13243
        - 360*lmMsq*(23 + lmMsq - 24*lmMst1) - 7008*lmMst1 - 4146*pow2(lmMst1))
        ) - 180*(DN + (31 + 78*lmMsq)*pow2(lmMst1)) - 18*(-1286 + 660*lmMsq +
        717*lmMst1)*pow2(lmMst2) + 16242*pow3(lmMst1) + 9102*pow3(lmMst2)) + 6*
        (413 + 408*B4 + 12*DN + 8520*lmMsq - 720*pow2(lmMsq) - 2*(lmMst1*(1573
        + 540*lmMsq - 180*pow2(lmMsq)) + lmMst2*(7603 - 1260*lmMsq - 354*lmMst1
        + 180*pow2(lmMsq) - 180*pow2(lmMst1))) - 18*(51 + 20*lmMsq)*pow2(
        lmMst1) + 6*(-797 + 60*lmMsq + 4*lmMst1)*pow2(lmMst2) + 200*pow3(
        lmMst1) - 584*pow3(lmMst2))*pow4(Mst2)) + pow2(Mst1)*(480*(280 - 63*
        lmMst1 + 36*lmMsq*(-2 + lmMst1 - lmMst2) + 135*lmMst2 - 18*pow2(lmMst1)
        + 18*pow2(lmMst2))*pow2(Msq)*pow4(Mst2) - 5*(478 - 24*lmMsq + 24*
        lmMst2)*pow6(Mst2)) + 5*(11 - 12*lmMsq + 12*lmMst2)*pow8(Mst2)))) +
        Tbeta*(-100*pow2(Mt)*pow2(MuSUSY)*(768*pow2(Mst1)*pow2(Mt)*(9893 - 54*
        B4 - 27*DN - 50*lmMst1 + 138*pow2(lmMst1) + lmMst2*(4538 - 300*lmMst1 +
        54*pow2(lmMst1)) + 54*(19 + 3*lmMst1)*pow2(lmMst2) + 24*lmMt*(29 + 35*
        lmMst2 - lmMst1*(35 + 18*lmMst2) + 9*pow2(lmMst1) + 9*pow2(lmMst2)) -
        234*pow3(lmMst1) + 18*pow3(lmMst2))*pow4(Msq) - 48*Mst1*Mt*s2t*(2*(3*
        pow2(Mst2)*(413 + 408*B4 + 12*DN + 8520*lmMsq - 720*pow2(lmMsq) + 2*
        lmMst1*(-1573 - 540*lmMsq + 180*pow2(lmMsq)) - 18*(51 + 20*lmMsq)*pow2(
        lmMst1) + lmMst2*(-15206 + 2520*lmMsq + 708*lmMst1 - 360*pow2(lmMsq) +
        360*pow2(lmMst1)) + 6*(-797 + 60*lmMsq + 4*lmMst1)*pow2(lmMst2) + 200*
        pow3(lmMst1) - 584*pow3(lmMst2)) + pow2(Mst1)*(59579 + 20088*B4 - 324*
        DN - 39240*lmMsq - 2160*pow2(lmMsq) + 24*lmMst1*(-56 + 1935*lmMsq +
        135*pow2(lmMsq)) - 18*(773 + 1620*lmMsq)*pow2(lmMst1) - 36*lmMst2*(-940
        + 1109*lmMst1 - 90*lmMsq*(-13 + 16*lmMst1) + 90*pow2(lmMsq) + 661*pow2(
        lmMst1)) - 90*(-355 + 252*lmMsq + 286*lmMst1)*pow2(lmMst2) + 33084*
        pow3(lmMst1) + 16452*pow3(lmMst2)))*pow4(Msq) + 90*(928 - 145*lmMst1 +
        180*lmMsq*(-2 + lmMst1 - lmMst2) + 505*lmMst2 - 90*pow2(lmMst1) + 90*
        pow2(lmMst2))*pow2(Mst2)*pow4(Mst1) + 240*pow2(Msq)*((557 - 126*lmMst1
        + 72*lmMsq*(-2 + lmMst1 - lmMst2) + 270*lmMst2 - 36*pow2(lmMst1) + 36*
        pow2(lmMst2))*pow2(Mst1)*pow2(Mst2) + (203 + 72*lmMsq*(-2 + 3*lmMst1 -
        3*lmMst2) + 366*lmMst2 + 6*lmMst1*(-37 + 90*lmMst2) - 378*pow2(lmMst1)
        - 162*pow2(lmMst2))*pow4(Mst1) - 3*pow4(Mst2)) + 5*(-467 + 12*lmMsq -
        12*lmMst2)*pow2(Mst1)*pow4(Mst2) + 5*(11 - 12*lmMsq + 12*lmMst2)*pow6(
        Mst2)) + pow2(s2t)*(pow4(Msq)*(864*pow2(Mst1)*pow2(Mst2)*(479 + 36*B4 +
        30*DN - 975*lmMsq + lmMst1*(1018 + 270*lmMsq - 270*pow2(lmMsq)) + 135*
        pow2(lmMsq) + lmMst2*(811 - 540*lmMsq + 258*lmMst1 + 270*pow2(lmMsq) -
        237*pow2(lmMst1)) + (-646 + 270*lmMsq)*pow2(lmMst1) - 3*(-372 + 90*
        lmMsq + 49*lmMst1)*pow2(lmMst2) - 65*pow3(lmMst1) + 449*pow3(lmMst2)) +
        (130639 - 24192*B4 - 1728*DN - 984960*lmMsq - 12*lmMst1*(262259 -
        181440*lmMsq + 19440*pow2(lmMsq)) + 72*(-21689 + 3240*lmMsq)*pow2(
        lmMst1) + 12*lmMst2*(347507 - 181440*lmMsq - 156324*lmMst1 + 19440*
        pow2(lmMsq) + 84312*pow2(lmMst1)) - 72*(-45823 + 3240*lmMsq + 27876*
        lmMst1)*pow2(lmMst2) - 240480*pow3(lmMst1) + 1235808*pow3(lmMst2))*
        pow4(Mst1) + 72*(6457 + 60*lmMsq*(67 - 54*lmMst1) - 982*lmMst1 + 30*(-
        89 + 66*lmMst1)*lmMst2 + 1620*pow2(lmMsq) + 834*pow2(lmMst1) + 162*
        pow2(lmMst2))*pow4(Mst2)) - 1440*pow2(Msq)*(2*(335 - 97*lmMst1 + 6*
        lmMsq*(-25 + 41*lmMst1 - 41*lmMst2) + 247*lmMst2 - 123*pow2(lmMst1) +
        123*pow2(lmMst2))*pow2(Mst2)*pow4(Mst1) - 3*(108 - 133*lmMst1 + 2*
        lmMsq*(32 + 9*lmMst1 - 9*lmMst2) + 69*lmMst2 - 9*pow2(lmMst1) + 9*pow2(
        lmMst2))*pow2(Mst1)*pow4(Mst2)) + 180*((3087 - 2306*lmMst1 + 12*lmMsq*(
        47 + 50*lmMst1 - 50*lmMst2) + 1742*lmMst2 - 300*pow2(lmMst1) + 300*
        pow2(lmMst2))*pow4(Mst1)*pow4(Mst2) - 54*pow2(Mst1)*pow6(Mst2)))) +
        pow2(Sbeta)*(-25*pow4(Mst2)*(180*(3195 - 2306*lmMst1 + 12*lmMsq*(47 +
        50*lmMst1 - 50*lmMst2) + 1742*lmMst2 - 300*pow2(lmMst1) + 300*pow2(
        lmMst2))*pow4(Mst1)*pow4(Mst2) - pow4(Msq)*(144*pow2(Mst1)*pow2(Mst2)*(
        3583 - 216*B4 - 180*DN + 9870*lmMsq + 810*pow2(lmMsq) + 10*lmMst1*(-709
        - 486*lmMsq + 162*pow2(lmMsq)) - 6*lmMst2*(1256 - 540*lmMsq - 72*lmMst1
        + 270*pow2(lmMsq) - 237*pow2(lmMst1)) - 30*(-157 + 54*lmMsq)*pow2(
        lmMst1) + 18*(-363 + 90*lmMsq + 49*lmMst1)*pow2(lmMst2) + 390*pow3(
        lmMst1) - 2694*pow3(lmMst2)) + (232169 + 86400*B4 + 53568*DN - 989280*
        lmMsq + 116640*pow2(lmMsq) - 12*lmMst1*(-414743 + 123120*lmMsq + 19440*
        pow2(lmMsq)) + 12*lmMst2*(-214703 + 103680*lmMsq + 181596*lmMst1 +
        19440*pow2(lmMsq) - 118440*pow2(lmMst1)) + 72*(5351 + 3240*lmMsq)*pow2(
        lmMst1) - 72*(19201 + 3240*lmMsq - 24348*lmMst1)*pow2(lmMst2) + 128160*
        pow3(lmMst1) - 459936*pow3(lmMst2))*pow4(Mst1) - 72*(6457 + 60*lmMsq*(
        67 - 54*lmMst1) - 982*lmMst1 + 30*(-89 + 66*lmMst1)*lmMst2 + 1620*pow2(
        lmMsq) + 834*pow2(lmMst1) + 162*pow2(lmMst2))*pow4(Mst2)) - 1440*pow2(
        Msq)*(2*(659 - 496*lmMst1 + 6*lmMsq*(7 + 50*lmMst1 - 50*lmMst2) + 454*
        lmMst2 - 150*pow2(lmMst1) + 150*pow2(lmMst2))*pow2(Mst2)*pow4(Mst1) -
        3*(108 - 133*lmMst1 + 2*lmMsq*(32 + 9*lmMst1 - 9*lmMst2) + 69*lmMst2 -
        9*pow2(lmMst1) + 9*pow2(lmMst2))*pow2(Mst1)*pow4(Mst2)))*pow4(s2t) +
        64*pow4(Mt)*(2*pow4(Msq)*(12*pow2(Mst1)*(pow2(Mst2)*(1477147 - 158575*
        lmMst2 + 30375*lmMsq*(3 + 2*lmMt) - 5*lmMst1*(17951 + 28065*lmMt + 15*
        lmMst2*(409 + 36*lmMt)) - 30375*pow2(lmMsq) + 150*(-134 + 9*lmMst2 -
        279*lmMt)*pow2(lmMst1) - 79425*pow2(lmMst2) + 15*lmMt*(27227 + 14075*
        lmMst2 + 2970*pow2(lmMst2)) - 32400*pow2(lmMt) + 13500*pow3(lmMst1) -
        14850*pow3(lmMst2)) + 50*pow2(MuSUSY)*(9893 - 54*B4 - 27*DN - 50*lmMst1
        + 138*pow2(lmMst1) + lmMst2*(4538 - 300*lmMst1 + 54*pow2(lmMst1)) + 54*
        (19 + 3*lmMst1)*pow2(lmMst2) + 24*lmMt*(29 + 35*lmMst2 - lmMst1*(35 +
        18*lmMst2) + 9*pow2(lmMst1) + 9*pow2(lmMst2)) - 234*pow3(lmMst1) + 18*
        pow3(lmMst2))) + (37614091 + 12134250*lmMst2 + 2700*(303 - 824*lmMst2 +
        322*lmMt)*pow2(lmMst1) + 6674400*pow2(lmMst2) - 405000*lmMsq*(13 +
        lmMst2*(33 - 10*lmMt) - 10*lmMst1*(1 + lmMst2 - lmMt) - 23*lmMt + 10*
        pow2(lmMst2)) - 180*lmMt*(-11973 - 38425*lmMst2 + 1410*pow2(lmMst2)) -
        129600*(19 + 10*lmMst2)*pow2(lmMt) + 30*lmMst1*(31687 - 472950*lmMt -
        30*lmMst2*(-2971 + 684*lmMt) + 5220*pow2(lmMst2) + 43200*pow2(lmMt)) +
        451800*pow3(lmMst1) + 1616400*pow3(lmMst2))*pow4(Mst1) - 6*(121294 -
        3182*lmMst2 - 10*lmMst1*(12992 + 582*lmMst2 - 465*lmMt) + 51960*lmMt -
        6330*lmMst2*lmMt - 60*lmMsq*(-1003 + 1770*lmMst1 + 255*lmMt) + 60750*
        pow2(lmMsq) + 82425*pow2(lmMst1) + 6075*pow2(lmMst2) + 30600*pow2(lmMt)
        )*pow4(Mst2)) + 360*pow2(Msq)*(1250*(11 - 6*lmMsq + 6*lmMt)*pow2(Mst2)*
        pow4(Mst1) + (6698 + 225*lmMst2 + 90*lmMst1*(62 - 15*lmMt) - 15*(7 +
        90*lmMst2)*lmMt + 150*lmMsq*(-38 + 9*lmMst1 + 9*lmMst2 + 18*lmMt) -
        2700*pow2(lmMsq))*pow2(Mst1)*pow4(Mst2)) + 15*((215363 + 41250*lmMst2 +
        120*lmMst1*(1319 - 510*lmMt) - 30*(-979 + 1500*lmMst2)*lmMt + 300*
        lmMsq*(-763 + 204*lmMst1 + 150*lmMst2 + 354*lmMt) - 106200*pow2(lmMsq))
        *pow4(Mst1)*pow4(Mst2) + 450*(115 - 48*lmMsq + lmMst2*(69 - 36*lmMt) -
        21*lmMt + 36*lmMsq*(lmMst2 + lmMt) - 36*pow2(lmMsq))*pow2(Mst1)*pow6(
        Mst2))) - 100*pow2(Mt)*pow2(s2t)*(-16*pow2(Msq)*pow2(Mst2)*(180*(577 -
        600*lmMsq + 354*lmMst1 + 246*lmMst2)*pow2(Mst2) + pow2(Msq)*(135281 -
        460266*lmMst2 + 1620*lmMsq*(103 - 100*lmMst1 + 82*lmMst2) + 14580*pow2(
        lmMsq) - 144*(-5257 + 186*lmMst2 + 720*lmMt)*pow2(lmMst1) + 6*lmMst1*(
        91201 + 41184*lmMt + 48*lmMst2*(-4651 + 720*lmMt) - 47376*pow2(lmMst2))
        + 643824*pow2(lmMst2) - 1728*lmMt*(66 + 143*lmMst2 + 60*pow2(lmMst2)) +
        133344*pow3(lmMst1) + 177696*pow3(lmMst2)))*pow4(Mst1) - 72*(120*(-9 +
        100*lmMsq - 118*lmMst1 + 18*lmMst2)*pow2(Msq)*pow2(Mst1) + 4*(4401 +
        60*lmMsq*(26 - 27*lmMst1) - 1462*lmMst1 + 6*(-74 + 69*lmMst1)*lmMst2 +
        810*pow2(lmMsq) + 993*pow2(lmMst1) + 81*pow2(lmMst2))*pow4(Msq) + 5*(
        497 + 1980*lmMsq - 2472*lmMst1 + 492*lmMst2)*pow4(Mst1))*pow6(Mst2) -
        pow2(MuSUSY)*(pow4(Msq)*(864*pow2(Mst1)*pow2(Mst2)*(479 + 36*B4 + 30*DN
        - 975*lmMsq + lmMst1*(1018 + 270*lmMsq - 270*pow2(lmMsq)) + 135*pow2(
        lmMsq) + lmMst2*(811 - 540*lmMsq + 258*lmMst1 + 270*pow2(lmMsq) - 237*
        pow2(lmMst1)) + (-646 + 270*lmMsq)*pow2(lmMst1) - 3*(-372 + 90*lmMsq +
        49*lmMst1)*pow2(lmMst2) - 65*pow3(lmMst1) + 449*pow3(lmMst2)) + (130639
        - 24192*B4 - 1728*DN - 984960*lmMsq - 12*lmMst1*(262259 - 181440*lmMsq
        + 19440*pow2(lmMsq)) + 72*(-21689 + 3240*lmMsq)*pow2(lmMst1) + 12*
        lmMst2*(347507 - 181440*lmMsq - 156324*lmMst1 + 19440*pow2(lmMsq) +
        84312*pow2(lmMst1)) - 72*(-45823 + 3240*lmMsq + 27876*lmMst1)*pow2(
        lmMst2) - 240480*pow3(lmMst1) + 1235808*pow3(lmMst2))*pow4(Mst1) + 72*(
        6457 + 60*lmMsq*(67 - 54*lmMst1) - 982*lmMst1 + 30*(-89 + 66*lmMst1)*
        lmMst2 + 1620*pow2(lmMsq) + 834*pow2(lmMst1) + 162*pow2(lmMst2))*pow4(
        Mst2)) - 1440*pow2(Msq)*(2*(335 - 97*lmMst1 + 6*lmMsq*(-25 + 41*lmMst1
        - 41*lmMst2) + 247*lmMst2 - 123*pow2(lmMst1) + 123*pow2(lmMst2))*pow2(
        Mst2)*pow4(Mst1) - 3*(108 - 133*lmMst1 + 2*lmMsq*(32 + 9*lmMst1 - 9*
        lmMst2) + 69*lmMst2 - 9*pow2(lmMst1) + 9*pow2(lmMst2))*pow2(Mst1)*pow4(
        Mst2)) + 180*((3087 - 2306*lmMst1 + 12*lmMsq*(47 + 50*lmMst1 - 50*
        lmMst2) + 1742*lmMst2 - 300*pow2(lmMst1) + 300*pow2(lmMst2))*pow4(Mst1)
        *pow4(Mst2) - 54*pow2(Mst1)*pow6(Mst2))) + pow2(Mst1)*(48*(266599 -
        864*B4 - 432*DN - 5220*lmMsq + 90*(115 - 108*lmMsq)*lmMst1 + 9720*pow2(
        lmMsq) + 7236*pow2(lmMst1) + 6*lmMst2*(20783 - 1620*lmMsq + 468*lmMst1
        + 432*pow2(lmMst1)) + 36*(953 + 168*lmMst1)*pow2(lmMst2) + 576*lmMt*(29
        + 35*lmMst2 - lmMst1*(35 + 18*lmMst2) + 9*(pow2(lmMst1) + pow2(lmMst2))
        ) - 7200*pow3(lmMst1) - 1440*pow3(lmMst2))*pow4(Msq)*pow4(Mst2) + 3240*
        (-13 + 12*lmMsq - 12*lmMst2)*pow8(Mst2))) + 200*Mst1*Mt*s2t*(-(pow2(
        Mst2)*pow2(s2t)*(180*pow2(Mst1)*(-3*(1961 - 290*lmMst1 + 4*lmMsq*(-181
        + 90*lmMst1 - 90*lmMst2) + 1014*lmMst2 - 180*pow2(lmMst1) + 180*pow2(
        lmMst2))*pow2(Mst1) + (163 - 12*lmMsq + 12*lmMst2)*pow2(Mst2))*pow4(
        Mst2) + 2880*pow2(Msq)*pow2(Mst2)*(-((563 - 126*lmMst1 + 72*lmMsq*(-2 +
        lmMst1 - lmMst2) + 270*lmMst2 - 36*pow2(lmMst1) + 36*pow2(lmMst2))*
        pow2(Mst1)*pow2(Mst2)) + (914 - 72*lmMsq*(2 + lmMst1 - lmMst2) + 174*
        lmMst2 - 30*lmMst1*(1 + 18*lmMst2) + 306*pow2(lmMst1) + 234*pow2(
        lmMst2))*pow4(Mst1) + 3*pow4(Mst2)) - pow4(Msq)*(24*pow2(Mst1)*pow2(
        Mst2)*(57101 + 17640*B4 - 396*DN - 90360*lmMsq + 2160*pow2(lmMsq) + 36*
        lmMst1*(487 + 1470*lmMsq + 30*pow2(lmMsq)) - 18*(467 + 1500*lmMsq)*
        pow2(lmMst1) - 12*lmMst2*(-10423 + 3681*lmMst1 - 90*lmMsq*(-53 + 48*
        lmMst1) + 90*pow2(lmMsq) + 2163*pow2(lmMst1)) - 18*(-3369 + 1380*lmMsq
        + 1438*lmMst1)*pow2(lmMst2) + 31884*pow3(lmMst1) + 19956*pow3(lmMst2))
        + (9306949 + 12960*lmMsq*(59 + 222*lmMst1 - 222*lmMst2) - 3056388*
        lmMst2 + 648*(-2917 + 1260*lmMst2)*pow2(lmMst1) + lmMst1*(1591044 +
        130896*lmMst2 - 705888*pow2(lmMst2)) + 1925208*pow2(lmMst2) - 309024*
        pow3(lmMst1) + 198432*pow3(lmMst2))*pow4(Mst1) + 72*(413 + 408*B4 + 12*
        DN + 8520*lmMsq - 720*pow2(lmMsq) + 2*lmMst1*(-1573 - 540*lmMsq + 180*
        pow2(lmMsq)) - 18*(51 + 20*lmMsq)*pow2(lmMst1) + lmMst2*(-15206 + 2520*
        lmMsq + 708*lmMst1 - 360*pow2(lmMsq) + 360*pow2(lmMst1)) + 6*(-797 +
        60*lmMsq + 4*lmMst1)*pow2(lmMst2) + 200*pow3(lmMst1) - 584*pow3(lmMst2)
        )*pow4(Mst2)))) + 8*pow2(Mt)*(-2*pow4(Msq)*(9*pow2(Mst2)*pow2(MuSUSY)*(
        413 + 408*B4 + 12*DN + 8520*lmMsq - 720*pow2(lmMsq) + 2*lmMst1*(-1573 -
        540*lmMsq + 180*pow2(lmMsq)) - 18*(51 + 20*lmMsq)*pow2(lmMst1) +
        lmMst2*(-15206 + 2520*lmMsq + 708*lmMst1 - 360*pow2(lmMsq) + 360*pow2(
        lmMst1)) + 6*(-797 + 60*lmMsq + 4*lmMst1)*pow2(lmMst2) + 200*pow3(
        lmMst1) - 584*pow3(lmMst2)) + pow2(Mst1)*(3*pow2(MuSUSY)*(59579 +
        20088*B4 - 324*DN - 39240*lmMsq - 2160*pow2(lmMsq) + 24*lmMst1*(-56 +
        1935*lmMsq + 135*pow2(lmMsq)) - 18*(773 + 1620*lmMsq)*pow2(lmMst1) -
        36*lmMst2*(-940 + 1109*lmMst1 - 90*lmMsq*(-13 + 16*lmMst1) + 90*pow2(
        lmMsq) + 661*pow2(lmMst1)) - 90*(-355 + 252*lmMsq + 286*lmMst1)*pow2(
        lmMst2) + 33084*pow3(lmMst1) + 16452*pow3(lmMst2)) + 4*pow2(Mst2)*(
        57151 + 161988*lmMst2 + 6480*lmMsq*(-5 + 3*lmMst1 - 3*lmMst2)*(1 + 2*
        lmMst2 - 2*lmMt) - 9*(1915 + 2184*lmMst2 - 4866*lmMt)*pow2(lmMst1) +
        117063*pow2(lmMst2) - 18*lmMt*(4150 + 2997*lmMst2 + 447*pow2(lmMst2)) -
        2592*(7 + 6*lmMst2)*pow2(lmMt) + 6*lmMst1*(716 + 6*lmMst2*(275 - 993*
        lmMt) - 4545*lmMt - 369*pow2(lmMst2) + 2592*pow2(lmMt)) - 8622*pow3(
        lmMst1) + 30492*pow3(lmMst2))) + 12*((134328 + 272545*lmMst2 - 18*(7074
        + 1807*lmMst2 - 2573*lmMt)*pow2(lmMst1) + 200058*pow2(lmMst2) - 1620*
        lmMsq*(20 + lmMst2*(107 - 30*lmMt) - 15*lmMst1*(5 + 2*lmMst2 - 2*lmMt)
        - 32*lmMt + 30*pow2(lmMst2)) - 6*lmMt*(11911 + 12661*lmMst2 + 2361*
        pow2(lmMst2)) - 864*(11 + 15*lmMst2)*pow2(lmMt) + lmMst1*(-153271 +
        14334*lmMt - 6*lmMst2*(265 + 5358*lmMt) + 5976*pow2(lmMst2) + 12960*
        pow2(lmMt)) - 4788*pow3(lmMst1) + 31338*pow3(lmMst2))*pow4(Mst1) + (
        36832 - 1067*lmMst2 + 30*lmMsq*(419 + 138*lmMt + 36*lmMst2*(1 + lmMt) -
        6*lmMst1*(29 + 6*lmMt)) + 540*(lmMst1 - lmMst2)*pow2(lmMsq) + 6*(757 +
        30*lmMst2 + 84*lmMt)*pow2(lmMst1) - 288*pow2(lmMst2) - 6*lmMt*(555 +
        473*lmMst2 + 108*pow2(lmMst2)) - 144*(17 + 6*lmMst2)*pow2(lmMt) +
        lmMst1*(-19027 - 3210*lmMt + 6*lmMst2*(131 + 24*lmMt) + 36*pow2(lmMst2)
        + 864*pow2(lmMt)) - 132*pow3(lmMst1) - 84*pow3(lmMst2))*pow4(Mst2))) +
        240*pow2(Msq)*(3*(18*(16 - 15*lmMst1 + 15*lmMst2)*(1 + 2*lmMst2 - 2*
        lmMt)*pow2(Mst2) - (203 + 72*lmMsq*(-2 + 3*lmMst1 - 3*lmMst2) + 366*
        lmMst2 + 6*lmMst1*(-37 + 90*lmMst2) - 378*pow2(lmMst1) - 162*pow2(
        lmMst2))*pow2(MuSUSY))*pow4(Mst1) + 9*pow2(MuSUSY)*pow4(Mst2) - 3*pow2(
        Mst1)*((557 - 126*lmMst1 + 72*lmMsq*(-2 + lmMst1 - lmMst2) + 270*lmMst2
        - 36*pow2(lmMst1) + 36*pow2(lmMst2))*pow2(Mst2)*pow2(MuSUSY) + 2*(331 -
        235*lmMst1 + 96*lmMst2 + (132 - 78*lmMst1 + 72*lmMst2)*lmMt + lmMsq*(7
        + 78*lmMst1 - 72*lmMst2 + 6*lmMt) - 6*pow2(lmMsq))*pow4(Mst2)) + (67 +
        12*lmMst2*(2 - 3*lmMt) - 66*lmMt + 6*lmMsq*(7 + 6*(lmMst2 + lmMt)) -
        36*pow2(lmMsq))*pow6(Mst2)) + 15*(6*pow2(Mst2)*(2*(lmMsq*(64 - 696*
        lmMst1 + 540*lmMst2 - 156*lmMt) + 8*lmMst1*(118 + 87*lmMt) - 3*(740 +
        231*lmMt + 15*lmMst2*(7 + 12*lmMt)) + 156*pow2(lmMsq))*pow2(Mst2) - 3*(
        928 - 145*lmMst1 + 180*lmMsq*(-2 + lmMst1 - lmMst2) + 505*lmMst2 - 90*
        pow2(lmMst1) + 90*pow2(lmMst2))*pow2(MuSUSY))*pow4(Mst1) + pow2(Mst1)*(
        48*(112 - 12*lmMsq + 51*lmMst2 - 3*(13 + 12*lmMst2)*lmMt + 36*lmMsq*(
        lmMst2 + lmMt) - 36*pow2(lmMsq))*pow2(Mst2) + (467 - 12*lmMsq + 12*
        lmMst2)*pow2(MuSUSY))*pow4(Mst2) + (-11 + 12*lmMsq - 12*lmMst2)*pow2(
        MuSUSY)*pow6(Mst2) + 4*(169 - 48*lmMsq + lmMst2*(87 - 36*lmMt) - 39*
        lmMt + 36*lmMsq*(lmMst2 + lmMt) - 36*pow2(lmMsq))*pow8(Mst2))))))))/(
        777600.*Tbeta*pow2(Sbeta)*pow4(Msq)*pow4(Mst2))))/pow2(Mgl));
}

/**
 * @return The off-diagonal (1, 2) = (2, 1) matrix element of the Higgs mass matrix as a double for the hierarchy 'H5g1'
 */
double H5g1::getS12() const {
   return -(MuSUSY*((Mt*(9*oneLoopFlag*s2t*(2*(2 - lmMst1 + lmMst2)*Mt*MuSUSY*s2t +
        4*(-lmMst1 + lmMst2)*Tbeta*pow2(Mt) + (-2 + lmMst1 - lmMst2)*Tbeta*
        pow2(Mst2)*pow2(s2t)) - (48*Al4p*twoLoopFlag*(-(Mst1*Mt*(4*Mgl*Mt*
        MuSUSY*s2t*(8 - 2*lmMst1 + 6*lmMst2 - pow2(lmMst1) + pow2(lmMst2)) + 4*
        Dmglst1*Mt*MuSUSY*s2t*(12 - 2*lmMst1 + 6*lmMst2 - pow2(lmMst1) + pow2(
        lmMst2)) + Dmglst1*Tbeta*(8*(5 + lmMst2 + (2 + lmMst2)*lmMt - lmMst1*(3
        + lmMt))*pow2(Mt) + 3*(-12 + 2*lmMst1 - 6*lmMst2 + pow2(lmMst1) - pow2(
        lmMst2))*pow2(Mst2)*pow2(s2t)) + Mgl*Tbeta*(8*(2 + lmMst2 + lmMt +
        lmMst2*lmMt - lmMst1*(2 + lmMt))*pow2(Mt) + 3*(-8 + 2*lmMst1 - 6*lmMst2
        + pow2(lmMst1) - pow2(lmMst2))*pow2(Mst2)*pow2(s2t)))) + s2t*(-2*
        Dmglst1*pow2(Mst1)*(2*lmMst2*Mt*MuSUSY*s2t + s2t*(pow2(lmMst1) - pow2(
        lmMst2))*(-2*Mt*MuSUSY + s2t*Tbeta*pow2(Mst2)) - 4*Tbeta*pow2(Mt) + (-1
        + lmMst1)*Tbeta*pow2(Mst2)*pow2(s2t) - lmMst2*Tbeta*(4*pow2(Mt) + pow2(
        Mst2)*pow2(s2t))) + pow2(Mst2)*(Dmglst1*(-4*(-1 + lmMst1)*Mt*MuSUSY*s2t
        + (4 - 8*lmMst1)*Tbeta*pow2(Mt) + 2*(-1 + lmMst1)*Tbeta*pow2(Mst2)*
        pow2(s2t)) + Mgl*(2*Mt*MuSUSY*s2t*(-7 - 7*lmMst2 + 2*lmMst1*(1 +
        lmMst2) - 2*pow2(lmMst2)) + 2*Tbeta*(2*lmMst1*(1 + lmMst2) - lmMst2*(4
        + 3*lmMst2) + pow2(lmMst1))*pow2(Mt) + Tbeta*(7 + 7*lmMst2 - 2*lmMst1*(
        1 + lmMst2) + 2*pow2(lmMst2))*pow2(Mst2)*pow2(s2t))))))/(Mgl*pow2(Mst2)
        )))/Tbeta + 4*threeLoopFlag*pow2(Al4p)*((-6*(1 - 2*lmMst2)*Mt*s2t*
        shiftst3*pow2(Mst2)*(4*pow2(Mt) + (1 + lmMst1 - lmMst2)*pow2(Mst1)*
        pow2(s2t)))/pow2(Mst1) + 8*Mt*s2t*T1ep*(-2*pow2(Mt) + pow2(Mst2)*pow2(
        s2t)) - 36*Mst1*pow2(Mt)*pow2(s2t)*(369.9166666666667 - (34*B4)/3. -
        DN/3. - 120*lmMsq + 20*pow2(lmMsq) - lmMst1*(13.166666666666666 - 30*
        lmMsq + 10*pow2(lmMsq)) + lmMst2*(355.8333333333333 - 70*lmMsq - (59*
        lmMst1)/3. + 10*pow2(lmMsq) - 10*pow2(lmMst1)) + ((-167 + 60*lmMsq)*
        pow2(lmMst1))/6. - ((-797 + 60*lmMsq + 4*lmMst1)*pow2(lmMst2))/6. - (
        50*pow3(lmMst1))/9. + (146*pow3(lmMst2))/9. + ((10*(Dmglst1 + Mgl)*
        pow2(Mst2))/(3.*pow2(Msq)) + Dmglst1*(314.9166666666667 - (34*B4)/3. -
        DN/3. - 180*lmMsq + lmMst1*(60.833333333333336 + 30*lmMsq - 10*pow2(
        lmMsq)) + 20*pow2(lmMsq) + lmMst2*(388.5 - 70*lmMsq - (59*lmMst1)/3. +
        10*pow2(lmMsq) - 10*pow2(lmMst1)) + ((-13 + 20*lmMsq)*pow2(lmMst1))/2.
         - ((-797 + 60*lmMsq + 4*lmMst1)*pow2(lmMst2))/
        6. - (50*pow3(lmMst1))/9. + (146*pow3(lmMst2))/9.))/Mgl) - (Dmglst1*Mt*
        pow2(Mst1)*(2323 + 144*B4 + 120*DN - 2880*lmMsq - 6*lmMst1*(-293 - 360*
        lmMsq + 180*pow2(lmMsq)) + 6*lmMst2*(571 - 360*lmMsq + 62*lmMst1 + 180*
        pow2(lmMsq) - 158*pow2(lmMst1)) + 30*(-121 + 36*lmMsq)*pow2(lmMst1) -
        6*(-735 + 180*lmMsq + 98*lmMst1)*pow2(lmMst2) - 260*pow3(lmMst1) +
        1796*pow3(lmMst2))*pow3(s2t))/(3.*Mgl) - (Mt*pow2(Mst2)*(Dmglst1*(36*(7
        - 2*lmMst1 - 60*lmMsq*(-5 + 6*lmMst1) - 226*lmMst2 + 220*lmMst1*lmMst2
        + 180*pow2(lmMsq) + 178*pow2(lmMst1) + 18*pow2(lmMst2))*pow2(Msq) +
        240*(82 - 93*lmMst1 + 6*lmMsq*(4 + 3*lmMst1 - 3*lmMst2) + 69*lmMst2 -
        9*pow2(lmMst1) + 9*pow2(lmMst2))*pow2(Mst1)) + Mgl*pow2(Msq)*(53965 +
        2592*B4 - 288*DN - 13320*lmMsq + 192*OepS2 - 1944*(-11 + 2*lmMst1 - 2*
        lmMst2)*S2 + 5400*pow2(lmMsq) - 6*lmMst1*(3781 - 300*lmMsq + 180*pow2(
        lmMsq)) + 6*lmMst2*(14137 - 3498*lmMst1 + 60*lmMsq*(-35 + 12*lmMst1) +
        180*pow2(lmMsq) - 18*pow2(lmMst1)) - 18*(-229 + 60*lmMsq)*pow2(lmMst1)
        - 18*(-2193 + 180*lmMsq + 442*lmMst1)*pow2(lmMst2) + 396*pow3(lmMst1) +
        7668*pow3(lmMst2)))*pow3(s2t))/(36.*Mgl*pow2(Msq)) + (s2t*pow3(Mt)*(
        4717 - 288*DN - 1080*lmMsq + 96*OepS2 - 972*(1 + 2*lmMst1 - 2*lmMst2)*
        S2 + 3240*pow2(lmMsq) - 24*lmMst1*(454 + 60*lmMsq + 45*pow2(lmMsq)) +
        24*lmMst2*(715 - 336*lmMst1 + 30*lmMsq*(-7 + 3*lmMst1) + 45*pow2(lmMsq)
        - 36*pow2(lmMst1)) + 2952*pow2(lmMst1) - 72*(-191 + 30*lmMsq + 44*
        lmMst1)*pow2(lmMst2) - 696*pow3(lmMst1) + 4728*pow3(lmMst2) + (324*(-(
        Dmglst1*(11.703703703703704 - (160*lmMsq)/9. + (20*(7 + 54*lmMsq)*
        lmMst1)/27. + (4*(16 - 23*lmMst1)*lmMst2)/9. - 20*pow2(lmMsq) - 34*
        pow2(lmMst1) - 2*pow2(lmMst2) + pow2(Mst1)*((20*(7 - 20*lmMsq + 26*
        lmMst1 - 6*lmMst2))/(9.*pow2(Msq)) + (164809 - 864*B4 - 432*DN - 14580*
        lmMsq - 9366*lmMst1 + 4860*pow2(lmMsq) + 1854*pow2(lmMst1) + 6*lmMst2*(
        23575 - 1620*lmMsq - 906*lmMst1 + 432*pow2(lmMst1)) + 18*(2167 + 336*
        lmMst1)*pow2(lmMst2) + 5184*lmMt*(2 + 3*lmMst2 - lmMst1*(3 + 2*lmMst2)
        + pow2(lmMst1) + pow2(lmMst2)) - 7200*pow3(lmMst1) - 1440*pow3(lmMst2))
        /(243.*pow2(Mst2))))) + (pow2(Mst2)*((-4*Mgl*(5108 + 15*lmMsq*(-218 +
        75*lmMst1 - 195*lmMst2) + 6270*lmMst2 - 375*lmMst1*(8 + 3*lmMst2) +
        900*pow2(lmMsq) + 2025*pow2(lmMst2)))/pow2(Msq) + (2700*(1 - 2*lmMst2)*
        Mgl)/pow2(Mst1) + (1125*Dmglst1*(13 - 12*lmMsq + 12*lmMst2)*pow2(Mst1))
        /pow4(Msq)))/4050.))/Mgl))/9. + (MuSUSY*pow2(Mt)*(pow2(s2t)*(270905 +
        12960*B4 - 66600*lmMsq + 960*OepS2 + 3240*(33 - 6*lmMst1 + 6*lmMst2)*S2
        - 1440*(DN + T1ep) + 27000*pow2(lmMsq) - 30*lmMst1*(3781 - 300*lmMsq +
        180*pow2(lmMsq)) + 30*lmMst2*(14065 - 3498*lmMst1 + 60*lmMsq*(-35 + 12*
        lmMst1) + 180*pow2(lmMsq) - 18*pow2(lmMst1)) + 90*(229 - 60*lmMsq)*
        pow2(lmMst1) - 90*(-2193 + 180*lmMsq + 442*lmMst1)*pow2(lmMst2) + (
        1080*(-1 + 2*lmMst2)*shiftst3*((2 - 2*lmMst1 + 2*lmMst2)*pow2(Mst1) +
        pow2(Mst2)))/pow2(Mst1) + 1980*pow3(lmMst1) + 38340*pow3(lmMst2) + (60*
        Dmglst1*(20*(82 - 93*lmMst1 + 6*lmMsq*(4 + 3*lmMst1 - 3*lmMst2) + 69*
        lmMst2 - 9*pow2(lmMst1) + 9*pow2(lmMst2))*pow2(Mst1)*pow2(Mst2) + pow2(
        Msq)*(3*(7 + 60*lmMsq*(5 - 6*lmMst1) - 2*lmMst1 + (-226 + 220*lmMst1)*
        lmMst2 + 180*pow2(lmMsq) + 178*pow2(lmMst1) + 18*pow2(lmMst2))*pow2(
        Mst2) + 4*pow2(Mst1)*(586 + 36*B4 + 30*DN - 495*lmMsq + 135*pow2(lmMsq)
        - 6*lmMst1*(-73 - 45*lmMsq + 45*pow2(lmMsq)) + 3*lmMst2*(229 - 180*
        lmMsq + 86*lmMst1 + 90*pow2(lmMsq) - 79*pow2(lmMst1)) + 18*(-43 + 15*
        lmMsq)*pow2(lmMst1) - 3*(-372 + 90*lmMsq + 49*lmMst1)*pow2(lmMst2) -
        65*pow3(lmMst1) + 449*pow3(lmMst2)))))/(Mgl*pow2(Msq)*pow2(Mst2)) - 12*
        pow2(Mst2)*((2*(939 + 30*lmMsq*(-12 + 5*lmMst1 - 5*lmMst2) + 760*lmMst2
        - 50*lmMst1*(8 + 3*lmMst2) + 150*pow2(lmMst2)))/pow2(Msq) + (90*(-1 +
        2*lmMst2))/pow2(Mst1) + (225*Dmglst1*pow2(Mst1))/(Mgl*pow4(Msq)))) + (
        20*Mst1*Mt*(s2t*((6*(Dmglst1*(11337 - 408*B4 - 12*DN - 6480*lmMsq +
        720*pow2(lmMsq) - 30*lmMst1*(-73 - 36*lmMsq + 12*pow2(lmMsq)) + 6*
        lmMst2*(2331 - 420*lmMsq - 118*lmMst1 + 60*pow2(lmMsq) - 60*pow2(
        lmMst1)) + 18*(-13 + 20*lmMsq)*pow2(lmMst1) - 6*(-797 + 60*lmMsq + 4*
        lmMst1)*pow2(lmMst2) - 200*pow3(lmMst1) + 584*pow3(lmMst2)) + Mgl*(
        13317 - 408*B4 - 12*DN - 4320*lmMsq + 720*pow2(lmMsq) - 6*lmMst1*(79 -
        180*lmMsq + 60*pow2(lmMsq)) + 6*lmMst2*(2135 - 420*lmMsq - 118*lmMst1 +
        60*pow2(lmMsq) - 60*pow2(lmMst1)) + 6*(-167 + 60*lmMsq)*pow2(lmMst1) -
        6*(-797 + 60*lmMsq + 4*lmMst1)*pow2(lmMst2) - 200*pow3(lmMst1) + 584*
        pow3(lmMst2))))/pow2(Mst2) + (5*(Dmglst1 + Mgl)*(144*pow2(Msq) + (-11 +
        12*lmMsq - 12*lmMst2)*pow2(Mst2)))/pow4(Msq)) + (96*Dmglst1*Mst1*Mt*(
        735 - 6*B4 - 3*DN - 78*lmMst1 + 6*pow2(lmMst1) + 6*lmMst2*(85 - 6*
        lmMst1 + pow2(lmMst1)) + 18*(7 + lmMst1)*pow2(lmMst2) + 24*lmMt*(2 + 3*
        lmMst2 - lmMst1*(3 + 2*lmMst2) + pow2(lmMst1) + pow2(lmMst2)) - 26*
        pow3(lmMst1) + 2*pow3(lmMst2)))/pow4(Mst2)))/Mgl))/(90.*Tbeta) + (Mt*(
        60*s2t*xMsq*(1 - 2*(lmMsq + z2))*pow2(Mgl)*pow2(Msq)*pow2(Mst2)*((-1 +
        shiftst1)*pow2(Mst2)*(-2*Mt*MuSUSY*s2t - 4*Tbeta*pow2(Mt) + Tbeta*pow2(
        Mst2)*pow2(s2t)) + pow2(Mst1)*(2*Mt*MuSUSY*s2t*(1 + 2*(lmMst1 - lmMst2)
        *(-1 + shiftst1) - 2*shiftst1 + shiftst2) + Tbeta*(4*(-1 + shiftst2)*
        pow2(Mt) - (shiftst1 - shiftst2 + (lmMst1 - lmMst2)*(-2 + shiftst1 +
        shiftst2))*pow2(Mst2)*pow2(s2t)))) + 2*pow2(Mst1)*pow2(z2)*(Mgl*pow2(
        Mst2)*(Mgl*(-64*Mst1*Mt*(-4*Mt*MuSUSY*s2t + Tbeta*pow2(Mt) + 3*Tbeta*
        pow2(Mst2)*pow2(s2t)) + 5*s2t*pow2(Mst2)*(-14*Mt*MuSUSY*s2t + 2*Tbeta*
        pow2(Mt) + 7*Tbeta*pow2(Mst2)*pow2(s2t))) + 16*Dmglst1*Mst1*(16*MuSUSY*
        s2t*pow2(Mt) - 2*Mt*(Mst1*MuSUSY + 6*Tbeta*pow2(Mst2))*pow2(s2t) - 4*
        Tbeta*pow3(Mt) + Mst1*Tbeta*pow2(Mst2)*pow3(s2t))) + 8*Mst1*xDmglst1*
        pow2(Dmglst1)*(60*Mt*s2t*pow2(Mst1)*(4*Mt*MuSUSY - 3*s2t*Tbeta*pow2(
        Mst2)) + pow2(Mst2)*(32*MuSUSY*s2t*pow2(Mt) - 6*Mt*(Mst1*MuSUSY + 4*
        Tbeta*pow2(Mst2))*pow2(s2t) - 8*Tbeta*pow3(Mt) + 3*Mst1*Tbeta*pow2(
        Mst2)*pow3(s2t)))) - z3*pow2(Mst1)*(4*Mgl*Mt*MuSUSY*(8*Mgl*s2t*pow2(
        Mst2)*((133 - 9*lmMst1 + 9*lmMst2)*Mst1*Mt - 41*s2t*pow2(Mst2)) +
        Dmglst1*(8*(58 - 9*lmMst1 + 9*lmMst2)*Mst1*Mt*s2t*pow2(Mst2) + 12*pow2(
        Mst1)*((52 + 8*lmMst1 - 8*lmMst2)*pow2(Mt) + (-35 - 46*lmMst1 + 46*
        lmMst2)*pow2(Mst2)*pow2(s2t)) + 195*pow2(s2t)*pow4(Mst2))) - 2*Mgl*
        Tbeta*pow2(Mst2)*(Dmglst1*(5*s2t*pow2(Mst2)*(-88*pow2(Mt) + 39*pow2(
        Mst2)*pow2(s2t)) + 4*Mst1*Mt*(4*(-37 + 12*lmMst1 - 12*lmMst2)*pow2(Mt)
        + 3*(58 - 9*lmMst1 + 9*lmMst2)*pow2(Mst2)*pow2(s2t)) + 3*s2t*pow2(Mst1)
        *(16*(53 + 4*lmMst1 - 4*lmMst2)*pow2(Mt) + (-205 - 184*lmMst1 + 184*
        lmMst2)*pow2(Mst2)*pow2(s2t))) - 4*Mgl*(-43*s2t*pow2(Mst2)*pow2(Mt) -
        3*Mst1*Mt*(8*(-9 + 2*lmMst1 - 2*lmMst2)*pow2(Mt) + (133 - 9*lmMst1 + 9*
        lmMst2)*pow2(Mst2)*pow2(s2t)) + 82*pow3(s2t)*pow4(Mst2))) + xDmglst1*
        pow2(Dmglst1)*(Mt*(-4*Mst1*pow2(Mst2)*(4*(559 + 18*lmMst1 - 18*lmMst2)*
        Mt*MuSUSY*s2t + 6*(565 + 16*lmMst1 - 16*lmMst2)*Tbeta*pow2(Mt) - 3*(559
        + 18*lmMst1 - 18*lmMst2)*Tbeta*pow2(Mst2)*pow2(s2t)) - 8*(54*(51 - 52*
        lmMst1 + 52*lmMst2)*Mt*MuSUSY*s2t + 3595*Tbeta*pow2(Mt) + 1053*(-1 + 2*
        lmMst1 - 2*lmMst2)*Tbeta*pow2(Mst2)*pow2(s2t))*pow3(Mst1)) + 7*s2t*(
        438*Mt*MuSUSY*s2t + 994*Tbeta*pow2(Mt) - 219*Tbeta*pow2(Mst2)*pow2(s2t)
        )*pow4(Mst2) + pow2(Mst1)*(-16*(1573 + 36*lmMst1 - 36*lmMst2)*s2t*
        Tbeta*pow2(Mst2)*pow2(Mt) - 48*(49 + 69*lmMst1 - 69*lmMst2)*Mt*MuSUSY*
        pow2(Mst2)*pow2(s2t) + 64*(220 + 9*lmMst1 - 9*lmMst2)*MuSUSY*pow3(Mt) +
        9*(301 + 184*lmMst1 - 184*lmMst2)*Tbeta*pow3(s2t)*pow4(Mst2)))) + 4*z4*
        pow2(Mst1)*(pow2(Mgl)*pow2(Mst2)*(s2t*pow2(Mst2)*(242*Mt*MuSUSY*s2t +
        18*Tbeta*pow2(Mt) - 121*Tbeta*pow2(Mst2)*pow2(s2t)) + Mst1*Mt*(-964*Mt*
        MuSUSY*s2t + 16*Tbeta*pow2(Mt) + 723*Tbeta*pow2(Mst2)*pow2(s2t))) +
        Dmglst1*Mst1*(Mgl*(Mt*pow2(Mst2)*(-964*Mt*MuSUSY*s2t + 16*Tbeta*pow2(
        Mt) + 723*Tbeta*pow2(Mst2)*pow2(s2t)) + 10*Mst1*(-24*s2t*Tbeta*pow2(
        Mst2)*pow2(Mt) - 106*Mt*MuSUSY*pow2(Mst2)*pow2(s2t) + 24*MuSUSY*pow3(
        Mt) + 53*Tbeta*pow3(s2t)*pow4(Mst2))) + Dmglst1*xDmglst1*(1230*Mt*s2t*
        pow2(Mst1)*(4*Mt*MuSUSY - 3*s2t*Tbeta*pow2(Mst2)) + Mt*pow2(Mst2)*(-
        964*Mt*MuSUSY*s2t + 16*Tbeta*pow2(Mt) + 723*Tbeta*pow2(Mst2)*pow2(s2t))
        + 15*Mst1*(-24*s2t*Tbeta*pow2(Mst2)*pow2(Mt) - 106*Mt*MuSUSY*pow2(Mst2)
        *pow2(s2t) + 24*MuSUSY*pow3(Mt) + 53*Tbeta*pow3(s2t)*pow4(Mst2)))))))/(
        Tbeta*pow2(Mgl)*pow2(Mst1)*pow4(Mst2)) + (2*Mst1*(80*Dmglst1*(17 + 3*
        lmMst2*(1 - 6*lmMt) - 33*lmMt + 6*lmMsq*(5 + 3*(lmMst2 + lmMt)) - 18*
        pow2(lmMsq))*pow2(Msq)*pow2(Mst2) + 80*Mgl*(17 + 3*lmMst2*(1 - 6*lmMt)
        - 33*lmMt + 6*lmMsq*(5 + 3*(lmMst2 + lmMt)) - 18*pow2(lmMsq))*pow2(Msq)
        *pow2(Mst2) - 8*Dmglst1*(1846 - 1057*lmMst2 + 30*lmMsq*(151 + 51*lmMt +
        9*lmMst2*(3 + 2*lmMt) - 6*lmMst1*(13 + 3*lmMt)) + 270*(lmMst1 - lmMst2)
        *pow2(lmMsq) + 6*(352 + 15*lmMst2 + 42*lmMt)*pow2(lmMst1) - 306*pow2(
        lmMst2) - 3*lmMt*(1105 + 608*lmMst2 + 108*pow2(lmMst2)) - 432*(2 +
        lmMst2)*pow2(lmMt) + 2*lmMst1*(-2534 - 192*lmMt + 3*lmMst2*(67 + 12*
        lmMt) + 9*pow2(lmMst2) + 216*pow2(lmMt)) - 66*pow3(lmMst1) - 42*pow3(
        lmMst2))*pow4(Msq) + 24*Mgl*(212 + 225*lmMst2 - 90*lmMsq*(7 + 2*lmMst2*
        lmMt - 2*lmMst1*(3 + lmMt) + 3*(lmMst2 + lmMt)) - 90*(lmMst1 - lmMst2)*
        pow2(lmMsq) - 6*(67 + 5*lmMst2 + 14*lmMt)*pow2(lmMst1) + 102*pow2(
        lmMst2) + 3*lmMt*(265 + 202*lmMst2 + 36*pow2(lmMst2)) + 144*(1 +
        lmMst2)*pow2(lmMt) - 6*lmMst1*(-54 + 37*lmMt + lmMst2*(22 + 4*lmMt) +
        pow2(lmMst2) + 24*pow2(lmMt)) + 22*pow3(lmMst1) + 14*pow3(lmMst2))*
        pow4(Msq) + 5*Dmglst1*(299 - 60*lmMsq + 6*lmMst2*(23 - 12*lmMt) - 78*
        lmMt + 72*lmMsq*(lmMst2 + lmMt) - 72*pow2(lmMsq))*pow4(Mst2) + 5*Mgl*(
        299 - 60*lmMsq + 6*lmMst2*(23 - 12*lmMt) - 78*lmMt + 72*lmMsq*(lmMst2 +
        lmMt) - 72*pow2(lmMsq))*pow4(Mst2))*pow4(Mt))/(9.*Mgl*pow2(Mst2)*pow4(
        Msq)) + 2*z2*(Mt*s2t*pow2(Mst2)*((24*shiftst3*pow2(Mt))/pow2(Mst1) + (
        333 - 30*lmMsq - 37*lmMst1 + 149*lmMst2 + 6*(1 + lmMst1)*shiftst3 - 6*
        lmMst2*shiftst3 + (Dmglst1*(-198 - (60*pow2(Mst1))/pow2(Msq)))/Mgl)*
        pow2(s2t)) + 18*s2t*(32.611111111111114 + (4*lmMst2)/3. + (Dmglst1*(44
        - (4*(2777 + 120*lmMst1 + 456*lmMst2)*pow2(Mst1))/(27.*pow2(Mst2))))/
        Mgl - (4*(lmMst1 + pow2(Mst2)/pow2(Mst1)))/3.)*pow3(Mt) + (2*MuSUSY*
        pow2(Mt)*(-((339 - 30*lmMsq - 37*lmMst1 + 149*lmMst2 - (9*Dmglst1*(22 +
        (4*pow2(Mst1)*(15/pow2(Msq) + (390 - 90*lmMsq - 139*lmMst1 + 315*
        lmMst2)/pow2(Mst2)))/9.))/Mgl + (6*pow2(Mst2))/pow2(Mst1) + shiftst3*(
        12*(-1 + lmMst1 - lmMst2) - (6*pow2(Mst2))/pow2(Mst1)))*pow2(s2t)) + (
        4*Mst1*Mt*(8*Dmglst1*(55 + 16*lmMst2)*Mst1*Mt + (627 - 60*lmMsq - 64*
        lmMst1 + 192*lmMst2)*(Dmglst1 + Mgl)*s2t*pow2(Mst2)))/(Mgl*pow4(Mst2)))
        )/Tbeta + (18*Mst1*((Dmglst1*(-681 + 180*lmMsq + 278*lmMst1 - 630*
        lmMst2)*Mst1*Mt*pow3(s2t))/9. + (Dmglst1 + Mgl)*(((-627 + 60*lmMsq +
        64*lmMst1 - 192*lmMst2)*pow2(Mt)*pow2(s2t))/3. + (8*(-60*pow2(Msq)*
        pow2(Mst2) + (8 - 96*lmMst1 + 54*lmMst2 + 42*lmMt)*pow4(Msq) - 15*pow4(
        Mst2))*pow4(Mt))/(27.*pow2(Mst2)*pow4(Msq)))))/Mgl)) + (16*Al4p*(-9*
        xDmglst1*pow2(Dmglst1)*(Mst1*twoLoopFlag*(16 - 2*lmMst1 + 6*lmMst2 -
        pow2(lmMst1) + pow2(lmMst2) - (2*(6 + 7*lmMst2 - lmMst1*(7 + 6*lmMst2)
        + 3*pow2(lmMst1) + 3*pow2(lmMst2))*pow2(Mst1))/pow2(Mst2))*pow2(Mt)*
        pow2(s2t) + s2t*twoLoopFlag*(3.7777777777777777 - 4*lmMst1 + (4*(1 +
        lmMst2)*pow2(Mst1))/pow2(Mst2))*pow3(Mt) + Mt*twoLoopFlag*((
        2.6666666666666665 - lmMst1 + lmMst2 - pow2(lmMst1) + pow2(lmMst2))*
        pow2(Mst1) + ((-4 + 3*lmMst1)*pow2(Mst2))/3.)*pow3(s2t) + (2*MuSUSY*
        s2t*twoLoopFlag*pow2(Mt)*(2*Mst1*Mt*(-16 + 2*lmMst1 - 6*lmMst2 + pow2(
        lmMst1) - pow2(lmMst2))*pow2(Mst2) - s2t*(4 + 3*lmMst2*(1 + lmMst2) -
        3*pow2(lmMst1))*pow2(Mst1)*pow2(Mst2) + 4*Mt*(4 + 7*lmMst2 - lmMst1*(7
        + 6*lmMst2) + 3*(pow2(lmMst1) + pow2(lmMst2)))*pow3(Mst1) + (4 - 3*
        lmMst1)*s2t*pow4(Mst2)))/(3.*Tbeta*pow4(Mst2)) - (4*twoLoopFlag*((151 +
        51*lmMt + 18*lmMst2*(1 + lmMt) - 3*lmMst1*(23 + 6*lmMt))*Mst1*pow2(
        Mst2) + (61 + 3*lmMst1*(-5 + 18*lmMst2 - 18*lmMt) + 105*lmMt + 18*
        lmMst2*(-5 + 3*lmMt) - 54*pow2(lmMst2))*pow3(Mst1))*pow4(Mt))/(27.*
        pow4(Mst2)) - (Al4p*MuSUSY*threeLoopFlag*pow2(Mt)*(3*pow2(s2t)*(6457 +
        4020*lmMsq - 2*(491 + 1620*lmMsq)*lmMst1 + 30*(-89 + 66*lmMst1)*lmMst2
        + 1620*pow2(lmMsq) + 834*pow2(lmMst1) + 162*pow2(lmMst2) + (3*pow2(
        Mst1)*(20*(108 - 133*lmMst1 + 2*lmMsq*(32 + 9*lmMst1 - 9*lmMst2) + 69*
        lmMst2 - 9*pow2(lmMst1) + 9*pow2(lmMst2))*pow2(Msq)*pow2(Mst2) + 4*(479
        + 36*B4 + 30*DN - 975*lmMsq + lmMst1*(1018 + 270*lmMsq - 270*pow2(
        lmMsq)) + 135*pow2(lmMsq) + lmMst2*(811 - 540*lmMsq + 258*lmMst1 + 270*
        pow2(lmMsq) - 237*pow2(lmMst1)) + (-646 + 270*lmMsq)*pow2(lmMst1) - 3*(
        -372 + 90*lmMsq + 49*lmMst1)*pow2(lmMst2) - 65*pow3(lmMst1) + 449*pow3(
        lmMst2))*pow4(Msq) - 45*pow4(Mst2)))/(pow2(Mst2)*pow4(Msq))) + (2*Mst1*
        Mt*(16*Mst1*Mt*(9893 - 54*B4 - 27*DN - 50*lmMst1 + 138*pow2(lmMst1) +
        lmMst2*(4538 - 300*lmMst1 + 54*pow2(lmMst1)) + 54*(19 + 3*lmMst1)*pow2(
        lmMst2) + 24*lmMt*(29 + 35*lmMst2 - lmMst1*(35 + 18*lmMst2) + 9*pow2(
        lmMst1) + 9*pow2(lmMst2)) - 234*pow3(lmMst1) + 18*pow3(lmMst2)) + (s2t*
        (-6*(pow2(Mst2)*(413 + 408*B4 + 12*DN + 8520*lmMsq - 720*pow2(lmMsq) +
        2*lmMst1*(-1573 - 540*lmMsq + 180*pow2(lmMsq)) - 18*(51 + 20*lmMsq)*
        pow2(lmMst1) + lmMst2*(-15206 + 2520*lmMsq + 708*lmMst1 - 360*pow2(
        lmMsq) + 360*pow2(lmMst1)) + 6*(-797 + 60*lmMsq + 4*lmMst1)*pow2(
        lmMst2) + 200*pow3(lmMst1) - 584*pow3(lmMst2)) + 4*pow2(Mst1)*(6129 +
        684*B4 - 18*DN - 3630*lmMsq + 6*(-251 + 600*lmMsq)*lmMst1 + lmMst2*(
        6914 - 1122*lmMst1 + 720*lmMsq*(-5 + 3*lmMst1) - 1179*pow2(lmMst1)) -
        3*(997 + 360*lmMsq)*pow2(lmMst1) - 9*(-393 + 120*lmMsq + 125*lmMst1)*
        pow2(lmMst2) + 1353*pow3(lmMst1) + 951*pow3(lmMst2)))*pow4(Msq) - 240*
        pow2(Msq)*((311 - 66*lmMst1 + 36*lmMsq*(-2 + lmMst1 - lmMst2) + 138*
        lmMst2 - 18*pow2(lmMst1) + 18*pow2(lmMst2))*pow2(Mst1)*pow2(Mst2) - 3*
        pow4(Mst2)) + 1080*pow2(Mst1)*pow4(Mst2) + 5*(-11 + 12*lmMsq - 12*
        lmMst2)*pow6(Mst2)))/pow4(Msq)))/pow4(Mst2)))/(324.*Tbeta) + Al4p*Mt*
        threeLoopFlag*(-((-((((6457 + 60*lmMsq*(67 - 54*lmMst1) - 982*lmMst1 +
        30*(-89 + 66*lmMst1)*lmMst2 + 1620*pow2(lmMsq) + 834*pow2(lmMst1) +
        162*pow2(lmMst2))*pow2(Msq) + 60*(108 - 133*lmMst1 + 2*lmMsq*(32 + 9*
        lmMst1 - 9*lmMst2) + 69*lmMst2 - 9*pow2(lmMst1) + 9*pow2(lmMst2))*pow2(
        Mst1))*pow2(Mst2))/pow2(Msq)) + pow2(Mst1)*(709 - 432*B4 - 360*DN +
        15720*lmMsq + 2*lmMst1*(-6599 - 3240*lmMsq + 1620*pow2(lmMsq)) - 18*
        lmMst2*(689 - 360*lmMsq + 62*lmMst1 + 180*pow2(lmMsq) - 158*pow2(
        lmMst1)) + 162*(53 - 20*lmMsq)*pow2(lmMst1) + 18*(-735 + 180*lmMsq +
        98*lmMst1)*pow2(lmMst2) + 780*pow3(lmMst1) - 5388*pow3(lmMst2)))*pow3(
        s2t))/216. - s2t*pow2(Mt)*(163 + (520*lmMsq)/9. - (2*(731 + 810*lmMsq)*
        lmMst1)/27. + (2*(-74 + 69*lmMst1)*lmMst2)/9. + 30*pow2(lmMsq) + (331*
        pow2(lmMst1))/9. + 3*pow2(lmMst2) - pow2(Mst1)*((10*(9 - 100*lmMsq +
        118*lmMst1 - 18*lmMst2))/(9.*pow2(Msq)) + (1534.2283950617284 - (16*B4)
        /3. - (8*DN)/3. - 90*lmMsq + (1507*lmMst1)/27. + 30*pow2(lmMsq) + (215*
        pow2(lmMst1))/9. + lmMst2*(848.4074074074074 - 60*lmMsq - 30*lmMst1 +
        16*pow2(lmMst1)) + (7*(289 + 48*lmMst1)*pow2(lmMst2))/9. + (32*lmMt*(29
        + 35*lmMst2 - lmMst1*(35 + 18*lmMst2) + 9*pow2(lmMst1) + 9*pow2(lmMst2)
        ))/9. - (400*pow3(lmMst1))/9. - (80*pow3(lmMst2))/9.)/pow2(Mst2) - ((
        5.416666666666667 - 5*lmMsq + 5*lmMst2)*pow2(Mst2))/pow4(Msq))) - Mt*
        pow2(s2t)*((Mst1*(413 + 408*B4 + 12*DN + 8520*lmMsq - 720*pow2(lmMsq) +
        2*lmMst1*(-1573 - 540*lmMsq + 180*pow2(lmMsq)) - 2*lmMst2*(7603 - 1260*
        lmMsq - 354*lmMst1 + 180*pow2(lmMsq) - 180*pow2(lmMst1)) - 18*(51 + 20*
        lmMsq)*pow2(lmMst1) + 6*(-797 + 60*lmMsq + 4*lmMst1)*pow2(lmMst2) +
        200*pow3(lmMst1) - 584*pow3(lmMst2)))/36. + ((10*(311 - 66*lmMst1 + 36*
        lmMsq*(-2 + lmMst1 - lmMst2) + 138*lmMst2 - 18*pow2(lmMst1) + 18*pow2(
        lmMst2)))/(9.*pow2(Msq)) + (354.6111111111111 + 76*B4 - 2*DN - 460*
        lmMsq + (-140.77777777777777 + 400*lmMsq)*lmMst1 + lmMst2*(
        802.1111111111111 - (374*lmMst1)/3. + 80*lmMsq*(-5 + 3*lmMst1) - 131*
        pow2(lmMst1)) - ((901 + 360*lmMsq)*pow2(lmMst1))/3. + (393 - 120*lmMsq
        - 125*lmMst1)*pow2(lmMst2) + (451*pow3(lmMst1))/3. + (317*pow3(lmMst2))
        /3.)/pow2(Mst2))*pow3(Mst1) - (5*pow2(Mst2)*(2*Mst1*pow2(Msq) + 3*pow3(
        Mst1)))/(3.*pow4(Msq))) + (Mst1*pow3(Mt)*(4*(pow2(Mst2)*(33943 - 2918*
        lmMst2 + 60*lmMsq*(244 + 69*lmMt + 9*lmMst2*(3 + 2*lmMt) - 6*lmMst1*(16
        + 3*lmMt)) + 540*(lmMst1 - lmMst2)*pow2(lmMsq) + 6*(799 + 30*lmMst2 +
        84*lmMt)*pow2(lmMst1) - 612*pow2(lmMst2) - 6*lmMt*(963 + 617*lmMst2 +
        108*pow2(lmMst2)) - 144*(17 + 6*lmMst2)*pow2(lmMt) + 2*lmMst1*(-10100 -
        1173*lmMt + lmMst2*(429 + 72*lmMt) + 18*pow2(lmMst2) + 432*pow2(lmMt))
        - 132*pow3(lmMst1) - 84*pow3(lmMst2)) + pow2(Mst1)*(32085 + 23606*
        lmMst2 + (3606 - 3402*lmMst2 + 7506*lmMt)*pow2(lmMst1) + 180*lmMsq*(31
        + 6*lmMst1*(-1 + 6*lmMst2 - 6*lmMt) + 84*lmMt + 6*lmMst2*(-13 + 6*lmMt)
        - 36*pow2(lmMst2)) + 18846*pow2(lmMst2) - 6*lmMt*(2850 + 1348*lmMst2 +
        189*pow2(lmMst2)) - 2*lmMst1*(5093 + 5700*lmMt + 6*lmMst2*(-173 + 531*
        lmMt) + 162*pow2(lmMst2) - 1296*pow2(lmMt)) - 144*(35 + 18*lmMst2)*
        pow2(lmMt) - 1368*pow3(lmMst1) + 5094*pow3(lmMst2)))*pow4(Msq) + 80*
        pow2(Msq)*(3*(245 + 6*lmMst1*(-25 + 6*lmMsq - 6*lmMt) + 84*lmMt +
        lmMst2*(66 - 36*lmMsq + 36*lmMt))*pow2(Mst1)*pow2(Mst2) + (-17 + 33*
        lmMt - 6*lmMsq*(5 + 3*lmMst2 + 3*lmMt) + 3*lmMst2*(-1 + 6*lmMt) + 18*
        pow2(lmMsq))*pow4(Mst2)) - 5*(6*(185 + 66*lmMst2 - 6*(13 + 12*lmMst2)*
        lmMt + 12*lmMsq*(1 + 6*(lmMst2 + lmMt)) - 72*pow2(lmMsq))*pow2(Mst1)*
        pow4(Mst2) + (299 - 60*lmMsq + 138*lmMst2 - 6*(13 + 12*lmMst2)*lmMt +
        72*lmMsq*(lmMst2 + lmMt) - 72*pow2(lmMsq))*pow6(Mst2))))/(162.*pow4(
        Msq)*pow4(Mst2)))) + Mt*z2*((3*Mgl*s2t*twoLoopFlag*(Mgl*(2*Mst1*Mt*(4*
        Mt*MuSUSY - 3*s2t*Tbeta*pow2(Mst2)) + s2t*pow2(Mst2)*(-2*Mt*MuSUSY +
        s2t*Tbeta*pow2(Mst2))) + 2*Dmglst1*Mst1*(Mt*s2t*(4*Mst1*MuSUSY - 3*
        Tbeta*pow2(Mst2)) + 4*MuSUSY*pow2(Mt) - 2*Mst1*Tbeta*pow2(Mst2)*pow2(
        s2t))))/(Tbeta*pow2(Mst2)) - xDmglst1*pow2(Dmglst1)*((Mst1*twoLoopFlag*
        ((12*Mt*MuSUSY*s2t*(12*Mt*pow2(Mst1) - 2*Mt*pow2(Mst2) - 3*Mst1*s2t*
        pow2(Mst2)))/Tbeta + 18*(pow2(Mst1)*(-6*Mt*pow2(Mst2)*pow2(s2t) + 8*
        pow3(Mt)) + (Mt + Mst1*s2t)*pow2(s2t)*pow4(Mst2))))/pow4(Mst2) + Al4p*
        threeLoopFlag*((Mt*MuSUSY*(30*s2t*pow2(Mst1)*pow2(Mst2)*(16*Mst1*Mt -
        3*s2t*pow2(Mst2)) + pow2(Msq)*(-6*pow2(Mst1)*(8*(55 + 16*lmMst2)*pow2(
        Mt) + (390 - 90*lmMsq - 139*lmMst1 + 315*lmMst2)*pow2(Mst2)*pow2(s2t))
        + Mt*s2t*(4*(-627 + 60*lmMsq + 64*lmMst1 - 192*lmMst2)*Mst1*pow2(Mst2)
        + 8*(1187 - 360*lmMsq - 171*lmMst1 + 747*lmMst2)*pow3(Mst1)) - 297*
        pow2(s2t)*pow4(Mst2))))/(Tbeta*pow2(Msq)*pow4(Mst2)) + 9*(-(s2t*(66 - (
        2*(2777 + 120*lmMst1 + 456*lmMst2)*pow2(Mst1))/(9.*pow2(Mst2)))*pow2(
        Mt)) + ((33*pow2(Mst2))/2. + pow2(Mst1)*(113.5 - 30*lmMsq - (139*
        lmMst1)/3. + 105*lmMst2 + (5*pow2(Mst2))/pow2(Msq)))*pow3(s2t) + Mst1*(
        Mt*(209 - 20*lmMsq - (64*lmMst1)/3. + 64*lmMst2 - pow2(Mst1)*(40/pow2(
        Msq) + (791.3333333333334 - 240*lmMsq - 114*lmMst1 + 498*lmMst2)/pow2(
        Mst2)))*pow2(s2t) + (8*pow3(Mt)*((3*(1391 - 360*lmMsq - 66*lmMst1 +
        813*lmMst2 + 21*lmMt)*pow2(Mst1) - 2*(4 - 48*lmMst1 + 27*lmMst2 + 21*
        lmMt)*pow2(Mst2))*pow4(Msq) + 60*pow2(Msq)*pow4(Mst2) + 15*(6*pow2(
        Mst1)*pow4(Mst2) + pow6(Mst2))))/(27.*pow4(Msq)*pow4(Mst2)))))))))/
        pow2(Mgl)))/144.;
}

/**
 * @return returns the susy log^0 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double H5g1::calc_coef_at_as2_no_sm_logs_log0() const {

   const double result =
      1555.7037037037037 - (17720*Dmglst1)/(27.*Mgl) + (5936*Mst1*s2t)/(9.*Mt)
        + (112112*Dmglst1*Mst1*s2t)/(27.*Mgl*Mt) + (6272*z2)/9. + (1056*
        Dmglst1*z2)/Mgl + (832*Mst1*s2t*z2)/(9.*Mt) + (832*Dmglst1*Mst1*s2t*z2)
        /(9.*Mgl*Mt) - (5216*z3)/3. - (856*Dmglst1*z3)/(3.*Mgl) + (1152*Mst1*
        s2t*z3)/Mt + (2368*Dmglst1*Mst1*s2t*z3)/(3.*Mgl*Mt) - (128*Mst1*s2t*z4)
        /(3.*Mt) - (128*Dmglst1*Mst1*s2t*z4)/(3.*Mgl*Mt) + (970352*pow2(
        Dmglst1))/(675.*pow2(Mgl)) + (589312*Mst1*s2t*pow2(Dmglst1))/(27.*Mt*
        pow2(Mgl)) + (1584*z2*pow2(Dmglst1))/pow2(Mgl) + (832*Mst1*s2t*z2*pow2(
        Dmglst1))/(9.*Mt*pow2(Mgl)) - (2692*z3*pow2(Dmglst1))/pow2(Mgl) - (
        9040*Mst1*s2t*z3*pow2(Dmglst1))/(Mt*pow2(Mgl)) - (128*Mst1*s2t*z4*pow2(
        Dmglst1))/(3.*Mt*pow2(Mgl)) + (160*pow2(Msq))/pow2(Mst1) - (320*z2*
        pow2(Msq))/pow2(Mst1) - (21936*pow2(Mst1))/(25.*pow2(Msq)) - (1440*
        Dmglst1*pow2(Mst1))/(Mgl*pow2(Msq)) - (32272*Cbeta*MuSUSY*s2t*pow2(
        Mst1))/(675.*Mt*Sbeta*pow2(Msq)) + (2560*z2*pow2(Mst1))/(3.*pow2(Msq))
        + (1280*Dmglst1*z2*pow2(Mst1))/(Mgl*pow2(Msq)) - (107168*pow2(Dmglst1)*
        pow2(Mst1))/(45.*pow2(Mgl)*pow2(Msq)) + (1920*z2*pow2(Dmglst1)*pow2(
        Mst1))/(pow2(Mgl)*pow2(Msq)) + (160*pow2(Msq))/pow2(Mst2) - (320*z2*
        pow2(Msq))/pow2(Mst2) - (39353408*pow2(Mst1))/(10125.*pow2(Mst2)) - (
        1025488*Dmglst1*pow2(Mst1))/(81.*Mgl*pow2(Mst2)) - (62992408*Cbeta*
        MuSUSY*s2t*pow2(Mst1))/(10125.*Mt*Sbeta*pow2(Mst2)) + (896*Cbeta*
        MuSUSY*OepS2*s2t*pow2(Mst1))/(81.*Mt*Sbeta*pow2(Mst2)) - (560*Cbeta*
        MuSUSY*S2*s2t*pow2(Mst1))/(Mt*Sbeta*pow2(Mst2)) - (448*Cbeta*MuSUSY*
        s2t*T1ep*pow2(Mst1))/(27.*Mt*Sbeta*pow2(Mst2)) + (9248*z2*pow2(Mst1))/(
        9.*pow2(Mst2)) + (19328*Dmglst1*z2*pow2(Mst1))/(9.*Mgl*pow2(Mst2)) - (
        140008*Cbeta*MuSUSY*s2t*z2*pow2(Mst1))/(27.*Mt*Sbeta*pow2(Mst2)) + (
        1536*z3*pow2(Mst1))/pow2(Mst2) + (24704*Dmglst1*z3*pow2(Mst1))/(3.*Mgl*
        pow2(Mst2)) + (11360*Cbeta*MuSUSY*s2t*z3*pow2(Mst1))/(27.*Mt*Sbeta*
        pow2(Mst2)) - (224*Cbeta*MuSUSY*s2t*z4*pow2(Mst1))/(27.*Mt*Sbeta*pow2(
        Mst2)) - (23634352*pow2(Dmglst1)*pow2(Mst1))/(675.*pow2(Mgl)*pow2(Mst2)
        ) + (9664*z2*pow2(Dmglst1)*pow2(Mst1))/(3.*pow2(Mgl)*pow2(Mst2)) + (
        74432*z3*pow2(Dmglst1)*pow2(Mst1))/(3.*pow2(Mgl)*pow2(Mst2)) - (216272*
        pow2(Mst2))/(675.*pow2(Msq)) - (10720*Mst1*s2t*pow2(Mst2))/(27.*Mt*
        pow2(Msq)) - (10720*Dmglst1*Mst1*s2t*pow2(Mst2))/(27.*Mgl*Mt*pow2(Msq))
        + (640*z2*pow2(Mst2))/(3.*pow2(Msq)) + (1280*Mst1*s2t*z2*pow2(Mst2))/(
        3.*Mt*pow2(Msq)) + (1280*Dmglst1*Mst1*s2t*z2*pow2(Mst2))/(3.*Mgl*Mt*
        pow2(Msq)) - (10720*Mst1*s2t*pow2(Dmglst1)*pow2(Mst2))/(27.*Mt*pow2(
        Mgl)*pow2(Msq)) + (1280*Mst1*s2t*z2*pow2(Dmglst1)*pow2(Mst2))/(3.*Mt*
        pow2(Mgl)*pow2(Msq)) + (16*pow2(Mst2))/pow2(Mst1) - (32*z2*pow2(Mst2))/
        pow2(Mst1) + (160*pow2(Msq)*pow2(s2t))/pow2(Mt) - (320*z2*pow2(Msq)*
        pow2(s2t))/pow2(Mt) + (32617079*pow2(Mst1)*pow2(s2t))/(10125.*pow2(Mt))
        - (32*DN*pow2(Mst1)*pow2(s2t))/(3.*pow2(Mt)) + (627124*Dmglst1*pow2(
        Mst1)*pow2(s2t))/(81.*Mgl*pow2(Mt)) - (128*B4*Dmglst1*pow2(Mst1)*pow2(
        s2t))/(3.*Mgl*pow2(Mt)) - (64*Dmglst1*DN*pow2(Mst1)*pow2(s2t))/(3.*Mgl*
        pow2(Mt)) - (160*OepS2*pow2(Mst1)*pow2(s2t))/(81.*pow2(Mt)) + (244*S2*
        pow2(Mst1)*pow2(s2t))/pow2(Mt) + (80*T1ep*pow2(Mst1)*pow2(s2t))/(27.*
        pow2(Mt)) + (80570*z2*pow2(Mst1)*pow2(s2t))/(27.*pow2(Mt)) + (49184*
        Dmglst1*z2*pow2(Mst1)*pow2(s2t))/(9.*Mgl*pow2(Mt)) - (2584*z3*pow2(
        Mst1)*pow2(s2t))/(27.*pow2(Mt)) - (5968*Dmglst1*z3*pow2(Mst1)*pow2(s2t)
        )/(3.*Mgl*pow2(Mt)) + (760*z4*pow2(Mst1)*pow2(s2t))/(27.*pow2(Mt)) + (
        320*Dmglst1*z4*pow2(Mst1)*pow2(s2t))/(Mgl*pow2(Mt)) + (533198*pow2(
        Dmglst1)*pow2(Mst1)*pow2(s2t))/(27.*pow2(Mgl)*pow2(Mt)) - (64*B4*pow2(
        Dmglst1)*pow2(Mst1)*pow2(s2t))/(pow2(Mgl)*pow2(Mt)) - (32*DN*pow2(
        Dmglst1)*pow2(Mst1)*pow2(s2t))/(pow2(Mgl)*pow2(Mt)) + (24592*z2*pow2(
        Dmglst1)*pow2(Mst1)*pow2(s2t))/(3.*pow2(Mgl)*pow2(Mt)) - (32126*z3*
        pow2(Dmglst1)*pow2(Mst1)*pow2(s2t))/(3.*pow2(Mgl)*pow2(Mt)) + (480*z4*
        pow2(Dmglst1)*pow2(Mst1)*pow2(s2t))/(pow2(Mgl)*pow2(Mt)) - (80*pow2(
        Msq)*pow2(Mst1)*pow2(s2t))/(pow2(Mst2)*pow2(Mt)) + (160*z2*pow2(Msq)*
        pow2(Mst1)*pow2(s2t))/(pow2(Mst2)*pow2(Mt)) - (4501*pow2(Mst2)*pow2(
        s2t))/(27.*pow2(Mt)) + (32*DN*pow2(Mst2)*pow2(s2t))/(3.*pow2(Mt)) + (
        1264*Dmglst1*pow2(Mst2)*pow2(s2t))/(9.*Mgl*pow2(Mt)) - (32*OepS2*pow2(
        Mst2)*pow2(s2t))/(9.*pow2(Mt)) + (36*S2*pow2(Mst2)*pow2(s2t))/pow2(Mt)
        + (16*T1ep*pow2(Mst2)*pow2(s2t))/(3.*pow2(Mt)) - (1222*z2*pow2(Mst2)*
        pow2(s2t))/(3.*pow2(Mt)) - (528*Dmglst1*z2*pow2(Mst2)*pow2(s2t))/(Mgl*
        pow2(Mt)) - (344*z3*pow2(Mst2)*pow2(s2t))/(3.*pow2(Mt)) + (880*Dmglst1*
        z3*pow2(Mst2)*pow2(s2t))/(3.*Mgl*pow2(Mt)) - (24*z4*pow2(Mst2)*pow2(
        s2t))/pow2(Mt) - (1956*pow2(Dmglst1)*pow2(Mst2)*pow2(s2t))/(pow2(Mgl)*
        pow2(Mt)) - (792*z2*pow2(Dmglst1)*pow2(Mst2)*pow2(s2t))/(pow2(Mgl)*
        pow2(Mt)) + (6958*z3*pow2(Dmglst1)*pow2(Mst2)*pow2(s2t))/(3.*pow2(Mgl)*
        pow2(Mt)) - (80*pow2(Msq)*pow2(Mst2)*pow2(s2t))/(pow2(Mst1)*pow2(Mt)) +
        (160*z2*pow2(Msq)*pow2(Mst2)*pow2(s2t))/(pow2(Mst1)*pow2(Mt)) - (24728*
        pow2(Mst1)*pow2(Mst2)*pow2(s2t))/(675.*pow2(Msq)*pow2(Mt)) + (560*
        Dmglst1*pow2(Mst1)*pow2(Mst2)*pow2(s2t))/(3.*Mgl*pow2(Msq)*pow2(Mt)) +
        (120*pow2(Dmglst1)*pow2(Mst1)*pow2(Mst2)*pow2(s2t))/(pow2(Mgl)*pow2(
        Msq)*pow2(Mt)) - (6688*pow2(Mst1)*pow2(MuSUSY)*pow2(s2t))/(45.*pow2(
        Msq)*pow2(Mt)) - (80*z2*pow2(Mst1)*pow2(MuSUSY)*pow2(s2t))/(pow2(Msq)*
        pow2(Mt)) + (6457954*pow2(Mst1)*pow2(MuSUSY)*pow2(s2t))/(10125.*pow2(
        Mst2)*pow2(Mt)) - (608*B4*pow2(Mst1)*pow2(MuSUSY)*pow2(s2t))/(3.*pow2(
        Mst2)*pow2(Mt)) + (16*DN*pow2(Mst1)*pow2(MuSUSY)*pow2(s2t))/(3.*pow2(
        Mst2)*pow2(Mt)) - (4352*OepS2*pow2(Mst1)*pow2(MuSUSY)*pow2(s2t))/(81.*
        pow2(Mst2)*pow2(Mt)) + (1568*S2*pow2(Mst1)*pow2(MuSUSY)*pow2(s2t))/(
        pow2(Mst2)*pow2(Mt)) + (2176*T1ep*pow2(Mst1)*pow2(MuSUSY)*pow2(s2t))/(
        27.*pow2(Mst2)*pow2(Mt)) + (2488*z2*pow2(Mst1)*pow2(MuSUSY)*pow2(s2t))/
        (27.*pow2(Mst2)*pow2(Mt)) - (140528*z3*pow2(Mst1)*pow2(MuSUSY)*pow2(
        s2t))/(27.*pow2(Mst2)*pow2(Mt)) - (208*z4*pow2(Mst1)*pow2(MuSUSY)*pow2(
        s2t))/(27.*pow2(Mst2)*pow2(Mt)) + (6688*pow2(Mst1)*pow2(MuSUSY)*pow2(
        s2t))/(45.*pow2(Msq)*pow2(Mt)*pow2(Sbeta)) + (80*z2*pow2(Mst1)*pow2(
        MuSUSY)*pow2(s2t))/(pow2(Msq)*pow2(Mt)*pow2(Sbeta)) - (6457954*pow2(
        Mst1)*pow2(MuSUSY)*pow2(s2t))/(10125.*pow2(Mst2)*pow2(Mt)*pow2(Sbeta))
        + (608*B4*pow2(Mst1)*pow2(MuSUSY)*pow2(s2t))/(3.*pow2(Mst2)*pow2(Mt)*
        pow2(Sbeta)) - (16*DN*pow2(Mst1)*pow2(MuSUSY)*pow2(s2t))/(3.*pow2(Mst2)
        *pow2(Mt)*pow2(Sbeta)) + (4352*OepS2*pow2(Mst1)*pow2(MuSUSY)*pow2(s2t))
        /(81.*pow2(Mst2)*pow2(Mt)*pow2(Sbeta)) - (1568*S2*pow2(Mst1)*pow2(
        MuSUSY)*pow2(s2t))/(pow2(Mst2)*pow2(Mt)*pow2(Sbeta)) - (2176*T1ep*pow2(
        Mst1)*pow2(MuSUSY)*pow2(s2t))/(27.*pow2(Mst2)*pow2(Mt)*pow2(Sbeta)) - (
        2488*z2*pow2(Mst1)*pow2(MuSUSY)*pow2(s2t))/(27.*pow2(Mst2)*pow2(Mt)*
        pow2(Sbeta)) + (140528*z3*pow2(Mst1)*pow2(MuSUSY)*pow2(s2t))/(27.*pow2(
        Mst2)*pow2(Mt)*pow2(Sbeta)) + (208*z4*pow2(Mst1)*pow2(MuSUSY)*pow2(s2t)
        )/(27.*pow2(Mst2)*pow2(Mt)*pow2(Sbeta)) + (256*Mst1*s2t*pow2(z2))/(3.*
        Mt) + (256*Dmglst1*Mst1*s2t*pow2(z2))/(3.*Mgl*Mt) + (256*Mst1*s2t*pow2(
        Dmglst1)*pow2(z2))/(3.*Mt*pow2(Mgl)) - (112*Cbeta*MuSUSY*s2t*pow2(Mst1)
        *pow2(z2))/(9.*Mt*Sbeta*pow2(Mst2)) + (116*pow2(Mst1)*pow2(s2t)*pow2(
        z2))/(9.*pow2(Mt)) - (20*pow2(Mst2)*pow2(s2t)*pow2(z2))/(3.*pow2(Mt)) +
        (2272*pow2(Mst1)*pow2(MuSUSY)*pow2(s2t)*pow2(z2))/(9.*pow2(Mst2)*pow2(
        Mt)) - (2272*pow2(Mst1)*pow2(MuSUSY)*pow2(s2t)*pow2(z2))/(9.*pow2(Mst2)
        *pow2(Mt)*pow2(Sbeta)) + (35200*s2t*pow3(Mst1))/(27.*Mt*pow2(Msq)) + (
        40960*Dmglst1*s2t*pow3(Mst1))/(9.*Mgl*Mt*pow2(Msq)) - (1280*s2t*z2*
        pow3(Mst1))/(3.*Mt*pow2(Msq)) - (1280*Dmglst1*s2t*z2*pow3(Mst1))/(3.*
        Mgl*Mt*pow2(Msq)) + (105920*s2t*pow2(Dmglst1)*pow3(Mst1))/(9.*Mt*pow2(
        Mgl)*pow2(Msq)) - (1280*s2t*z2*pow2(Dmglst1)*pow3(Mst1))/(3.*Mt*pow2(
        Mgl)*pow2(Msq)) + (1328*s2t*pow3(Mst1))/(27.*Mt*pow2(Mst2)) + (395536*
        Dmglst1*s2t*pow3(Mst1))/(81.*Mgl*Mt*pow2(Mst2)) + (29632*s2t*z2*pow3(
        Mst1))/(3.*Mt*pow2(Mst2)) + (265600*Dmglst1*s2t*z2*pow3(Mst1))/(9.*Mgl*
        Mt*pow2(Mst2)) - (3904*s2t*z3*pow3(Mst1))/(Mt*pow2(Mst2)) - (10816*
        Dmglst1*s2t*z3*pow3(Mst1))/(Mgl*Mt*pow2(Mst2)) + (914416*s2t*pow2(
        Dmglst1)*pow3(Mst1))/(81.*Mt*pow2(Mgl)*pow2(Mst2)) + (530656*s2t*z2*
        pow2(Dmglst1)*pow3(Mst1))/(9.*Mt*pow2(Mgl)*pow2(Mst2)) - (20160*s2t*z3*
        pow2(Dmglst1)*pow3(Mst1))/(Mt*pow2(Mgl)*pow2(Mst2)) - (3520*Cbeta*
        MuSUSY*pow3(Mst1))/(3.*Sbeta*pow2(Msq)*pow2(Mst2)) - (140480*Cbeta*
        Dmglst1*MuSUSY*pow3(Mst1))/(27.*Mgl*Sbeta*pow2(Msq)*pow2(Mst2)) - (
        140480*Cbeta*MuSUSY*pow2(Dmglst1)*pow3(Mst1))/(27.*Sbeta*pow2(Mgl)*
        pow2(Msq)*pow2(Mst2)) + (21760*s2t*pow2(MuSUSY)*pow3(Mst1))/(9.*Mt*
        pow2(Msq)*pow2(Mst2)) + (26240*Dmglst1*s2t*pow2(MuSUSY)*pow3(Mst1))/(3.
        *Mgl*Mt*pow2(Msq)*pow2(Mst2)) + (1280*s2t*z2*pow2(MuSUSY)*pow3(Mst1))/(
        3.*Mt*pow2(Msq)*pow2(Mst2)) + (1280*Dmglst1*s2t*z2*pow2(MuSUSY)*pow3(
        Mst1))/(Mgl*Mt*pow2(Msq)*pow2(Mst2)) + (26240*s2t*pow2(Dmglst1)*pow2(
        MuSUSY)*pow3(Mst1))/(3.*Mt*pow2(Mgl)*pow2(Msq)*pow2(Mst2)) + (1280*s2t*
        z2*pow2(Dmglst1)*pow2(MuSUSY)*pow3(Mst1))/(Mt*pow2(Mgl)*pow2(Msq)*pow2(
        Mst2)) + (5680*Cbeta*MuSUSY*pow2(s2t)*pow3(Mst1))/(3.*Sbeta*pow2(Msq)*
        pow2(Mt)) + (6640*Cbeta*Dmglst1*MuSUSY*pow2(s2t)*pow3(Mst1))/(Mgl*
        Sbeta*pow2(Msq)*pow2(Mt)) + (320*Cbeta*MuSUSY*z2*pow2(s2t)*pow3(Mst1))/
        (Sbeta*pow2(Msq)*pow2(Mt)) + (960*Cbeta*Dmglst1*MuSUSY*z2*pow2(s2t)*
        pow3(Mst1))/(Mgl*Sbeta*pow2(Msq)*pow2(Mt)) + (6640*Cbeta*MuSUSY*pow2(
        Dmglst1)*pow2(s2t)*pow3(Mst1))/(Sbeta*pow2(Mgl)*pow2(Msq)*pow2(Mt)) + (
        960*Cbeta*MuSUSY*z2*pow2(Dmglst1)*pow2(s2t)*pow3(Mst1))/(Sbeta*pow2(
        Mgl)*pow2(Msq)*pow2(Mt)) + (25108*Cbeta*MuSUSY*pow2(s2t)*pow3(Mst1))/(
        27.*Sbeta*pow2(Mst2)*pow2(Mt)) + (1152*B4*Cbeta*MuSUSY*pow2(s2t)*pow3(
        Mst1))/(Sbeta*pow2(Mst2)*pow2(Mt)) + (40084*Cbeta*Dmglst1*MuSUSY*pow2(
        s2t)*pow3(Mst1))/(9.*Mgl*Sbeta*pow2(Mst2)*pow2(Mt)) + (2368*B4*Cbeta*
        Dmglst1*MuSUSY*pow2(s2t)*pow3(Mst1))/(Mgl*Sbeta*pow2(Mst2)*pow2(Mt)) -
        (32*Cbeta*Dmglst1*DN*MuSUSY*pow2(s2t)*pow3(Mst1))/(Mgl*Sbeta*pow2(Mst2)
        *pow2(Mt)) - (4208*Cbeta*MuSUSY*z2*pow2(s2t)*pow3(Mst1))/(3.*Sbeta*
        pow2(Mst2)*pow2(Mt)) + (10512*Cbeta*Dmglst1*MuSUSY*z2*pow2(s2t)*pow3(
        Mst1))/(Mgl*Sbeta*pow2(Mst2)*pow2(Mt)) + (7072*Cbeta*MuSUSY*z3*pow2(
        s2t)*pow3(Mst1))/(Sbeta*pow2(Mst2)*pow2(Mt)) + (6464*Cbeta*Dmglst1*
        MuSUSY*z3*pow2(s2t)*pow3(Mst1))/(Mgl*Sbeta*pow2(Mst2)*pow2(Mt)) + (576*
        Cbeta*MuSUSY*z4*pow2(s2t)*pow3(Mst1))/(Sbeta*pow2(Mst2)*pow2(Mt)) - (
        5984*Cbeta*Dmglst1*MuSUSY*z4*pow2(s2t)*pow3(Mst1))/(Mgl*Sbeta*pow2(
        Mst2)*pow2(Mt)) + (40084*Cbeta*MuSUSY*pow2(Dmglst1)*pow2(s2t)*pow3(
        Mst1))/(9.*Sbeta*pow2(Mgl)*pow2(Mst2)*pow2(Mt)) + (2368*B4*Cbeta*
        MuSUSY*pow2(Dmglst1)*pow2(s2t)*pow3(Mst1))/(Sbeta*pow2(Mgl)*pow2(Mst2)*
        pow2(Mt)) - (32*Cbeta*DN*MuSUSY*pow2(Dmglst1)*pow2(s2t)*pow3(Mst1))/(
        Sbeta*pow2(Mgl)*pow2(Mst2)*pow2(Mt)) + (10512*Cbeta*MuSUSY*z2*pow2(
        Dmglst1)*pow2(s2t)*pow3(Mst1))/(Sbeta*pow2(Mgl)*pow2(Mst2)*pow2(Mt)) +
        (6464*Cbeta*MuSUSY*z3*pow2(Dmglst1)*pow2(s2t)*pow3(Mst1))/(Sbeta*pow2(
        Mgl)*pow2(Mst2)*pow2(Mt)) - (5984*Cbeta*MuSUSY*z4*pow2(Dmglst1)*pow2(
        s2t)*pow3(Mst1))/(Sbeta*pow2(Mgl)*pow2(Mst2)*pow2(Mt)) - (21760*s2t*
        pow2(MuSUSY)*pow3(Mst1))/(9.*Mt*pow2(Msq)*pow2(Mst2)*pow2(Sbeta)) - (
        26240*Dmglst1*s2t*pow2(MuSUSY)*pow3(Mst1))/(3.*Mgl*Mt*pow2(Msq)*pow2(
        Mst2)*pow2(Sbeta)) - (1280*s2t*z2*pow2(MuSUSY)*pow3(Mst1))/(3.*Mt*pow2(
        Msq)*pow2(Mst2)*pow2(Sbeta)) - (1280*Dmglst1*s2t*z2*pow2(MuSUSY)*pow3(
        Mst1))/(Mgl*Mt*pow2(Msq)*pow2(Mst2)*pow2(Sbeta)) - (26240*s2t*pow2(
        Dmglst1)*pow2(MuSUSY)*pow3(Mst1))/(3.*Mt*pow2(Mgl)*pow2(Msq)*pow2(Mst2)
        *pow2(Sbeta)) - (1280*s2t*z2*pow2(Dmglst1)*pow2(MuSUSY)*pow3(Mst1))/(
        Mt*pow2(Mgl)*pow2(Msq)*pow2(Mst2)*pow2(Sbeta)) - (1152*Cbeta*MuSUSY*
        pow2(s2t)*pow2(z2)*pow3(Mst1))/(Sbeta*pow2(Mst2)*pow2(Mt)) - (2432*
        Cbeta*Dmglst1*MuSUSY*pow2(s2t)*pow2(z2)*pow3(Mst1))/(Mgl*Sbeta*pow2(
        Mst2)*pow2(Mt)) - (2432*Cbeta*MuSUSY*pow2(Dmglst1)*pow2(s2t)*pow2(z2)*
        pow3(Mst1))/(Sbeta*pow2(Mgl)*pow2(Mst2)*pow2(Mt)) + (26775829*Cbeta*
        MuSUSY*pow2(Mst1)*pow3(s2t))/(20250.*Sbeta*pow3(Mt)) - (160*B4*Cbeta*
        MuSUSY*pow2(Mst1)*pow3(s2t))/(3.*Sbeta*pow3(Mt)) - (8*Cbeta*DN*MuSUSY*
        pow2(Mst1)*pow3(s2t))/(3.*Sbeta*pow3(Mt)) - (1888*Cbeta*MuSUSY*OepS2*
        pow2(Mst1)*pow3(s2t))/(81.*Sbeta*pow3(Mt)) + (1180*Cbeta*MuSUSY*S2*
        pow2(Mst1)*pow3(s2t))/(Sbeta*pow3(Mt)) + (944*Cbeta*MuSUSY*T1ep*pow2(
        Mst1)*pow3(s2t))/(27.*Sbeta*pow3(Mt)) - (10960*Cbeta*MuSUSY*z2*pow2(
        Mst1)*pow3(s2t))/(27.*Sbeta*pow3(Mt)) - (58456*Cbeta*MuSUSY*z3*pow2(
        Mst1)*pow3(s2t))/(27.*Sbeta*pow3(Mt)) + (8608*Cbeta*MuSUSY*z4*pow2(
        Mst1)*pow3(s2t))/(27.*Sbeta*pow3(Mt)) + (40*Cbeta*MuSUSY*pow2(Msq)*
        pow2(Mst1)*pow3(s2t))/(Sbeta*pow2(Mst2)*pow3(Mt)) - (80*Cbeta*MuSUSY*
        z2*pow2(Msq)*pow2(Mst1)*pow3(s2t))/(Sbeta*pow2(Mst2)*pow3(Mt)) + (8878*
        Mst1*pow2(Mst2)*pow3(s2t))/(3.*pow3(Mt)) - (272*B4*Mst1*pow2(Mst2)*
        pow3(s2t))/(3.*pow3(Mt)) - (8*DN*Mst1*pow2(Mst2)*pow3(s2t))/(3.*pow3(
        Mt)) + (7558*Dmglst1*Mst1*pow2(Mst2)*pow3(s2t))/(3.*Mgl*pow3(Mt)) - (
        272*B4*Dmglst1*Mst1*pow2(Mst2)*pow3(s2t))/(3.*Mgl*pow3(Mt)) - (8*
        Dmglst1*DN*Mst1*pow2(Mst2)*pow3(s2t))/(3.*Mgl*pow3(Mt)) + (1672*Mst1*
        z2*pow2(Mst2)*pow3(s2t))/pow3(Mt) + (1672*Dmglst1*Mst1*z2*pow2(Mst2)*
        pow3(s2t))/(Mgl*pow3(Mt)) - (2128*Mst1*z3*pow2(Mst2)*pow3(s2t))/(3.*
        pow3(Mt)) - (928*Dmglst1*Mst1*z3*pow2(Mst2)*pow3(s2t))/(3.*Mgl*pow3(Mt)
        ) - (1928*Mst1*z4*pow2(Mst2)*pow3(s2t))/(3.*pow3(Mt)) - (1928*Dmglst1*
        Mst1*z4*pow2(Mst2)*pow3(s2t))/(3.*Mgl*pow3(Mt)) - (826*Mst1*pow2(
        Dmglst1)*pow2(Mst2)*pow3(s2t))/(9.*pow2(Mgl)*pow3(Mt)) - (272*B4*Mst1*
        pow2(Dmglst1)*pow2(Mst2)*pow3(s2t))/(3.*pow2(Mgl)*pow3(Mt)) - (8*DN*
        Mst1*pow2(Dmglst1)*pow2(Mst2)*pow3(s2t))/(3.*pow2(Mgl)*pow3(Mt)) + (
        1672*Mst1*z2*pow2(Dmglst1)*pow2(Mst2)*pow3(s2t))/(pow2(Mgl)*pow3(Mt)) +
        (4472*Mst1*z3*pow2(Dmglst1)*pow2(Mst2)*pow3(s2t))/(3.*pow2(Mgl)*pow3(
        Mt)) - (1928*Mst1*z4*pow2(Dmglst1)*pow2(Mst2)*pow3(s2t))/(3.*pow2(Mgl)*
        pow3(Mt)) - (1420*Cbeta*MuSUSY*pow2(Mst1)*pow2(Mst2)*pow3(s2t))/(9.*
        Sbeta*pow2(Msq)*pow3(Mt)) - (40*Cbeta*MuSUSY*z2*pow2(Mst1)*pow2(Mst2)*
        pow3(s2t))/(Sbeta*pow2(Msq)*pow3(Mt)) + (716*Cbeta*MuSUSY*pow2(Mst1)*
        pow2(z2)*pow3(s2t))/(9.*Sbeta*pow3(Mt)) + (256*Mst1*pow2(Mst2)*pow2(z2)
        *pow3(s2t))/(3.*pow3(Mt)) + (256*Dmglst1*Mst1*pow2(Mst2)*pow2(z2)*pow3(
        s2t))/(3.*Mgl*pow3(Mt)) + (256*Mst1*pow2(Dmglst1)*pow2(Mst2)*pow2(z2)*
        pow3(s2t))/(3.*pow2(Mgl)*pow3(Mt)) - (264814*pow3(Mst1)*pow3(s2t))/(81.
        *pow3(Mt)) - (880*B4*pow3(Mst1)*pow3(s2t))/(3.*pow3(Mt)) + (8*DN*pow3(
        Mst1)*pow3(s2t))/(3.*pow3(Mt)) - (108106*Dmglst1*pow3(Mst1)*pow3(s2t))/
        (27.*Mgl*pow3(Mt)) - (2096*B4*Dmglst1*pow3(Mst1)*pow3(s2t))/(3.*Mgl*
        pow3(Mt)) + (40*Dmglst1*DN*pow3(Mst1)*pow3(s2t))/(3.*Mgl*pow3(Mt)) - (
        10840*z2*pow3(Mst1)*pow3(s2t))/(9.*pow3(Mt)) - (5176*Dmglst1*z2*pow3(
        Mst1)*pow3(s2t))/(Mgl*pow3(Mt)) - (1648*z3*pow3(Mst1)*pow3(s2t))/pow3(
        Mt) - (5536*Dmglst1*z3*pow3(Mst1)*pow3(s2t))/(3.*Mgl*pow3(Mt)) + (1352*
        z4*pow3(Mst1)*pow3(s2t))/(3.*pow3(Mt)) + (7912*Dmglst1*z4*pow3(Mst1)*
        pow3(s2t))/(3.*Mgl*pow3(Mt)) - (114202*pow2(Dmglst1)*pow3(Mst1)*pow3(
        s2t))/(27.*pow2(Mgl)*pow3(Mt)) - (3920*B4*pow2(Dmglst1)*pow3(Mst1)*
        pow3(s2t))/(3.*pow2(Mgl)*pow3(Mt)) + (88*DN*pow2(Dmglst1)*pow3(Mst1)*
        pow3(s2t))/(3.*pow2(Mgl)*pow3(Mt)) - (34520*z2*pow2(Dmglst1)*pow3(Mst1)
        *pow3(s2t))/(3.*pow2(Mgl)*pow3(Mt)) - (5320*z3*pow2(Dmglst1)*pow3(Mst1)
        *pow3(s2t))/(3.*pow2(Mgl)*pow3(Mt)) + (17752*z4*pow2(Dmglst1)*pow3(
        Mst1)*pow3(s2t))/(3.*pow2(Mgl)*pow3(Mt)) - (5920*pow2(Mst2)*pow3(Mst1)*
        pow3(s2t))/(9.*pow2(Msq)*pow3(Mt)) - (2240*Dmglst1*pow2(Mst2)*pow3(
        Mst1)*pow3(s2t))/(Mgl*pow2(Msq)*pow3(Mt)) - (320*z2*pow2(Mst2)*pow3(
        Mst1)*pow3(s2t))/(3.*pow2(Msq)*pow3(Mt)) - (320*Dmglst1*z2*pow2(Mst2)*
        pow3(Mst1)*pow3(s2t))/(Mgl*pow2(Msq)*pow3(Mt)) - (45040*pow2(Dmglst1)*
        pow2(Mst2)*pow3(Mst1)*pow3(s2t))/(9.*pow2(Mgl)*pow2(Msq)*pow3(Mt)) - (
        640*z2*pow2(Dmglst1)*pow2(Mst2)*pow3(Mst1)*pow3(s2t))/(pow2(Mgl)*pow2(
        Msq)*pow3(Mt)) + (896*pow2(z2)*pow3(Mst1)*pow3(s2t))/(3.*pow3(Mt)) + (
        2176*Dmglst1*pow2(z2)*pow3(Mst1)*pow3(s2t))/(3.*Mgl*pow3(Mt)) + (4096*
        pow2(Dmglst1)*pow2(z2)*pow3(Mst1)*pow3(s2t))/(3.*pow2(Mgl)*pow3(Mt)) +
        (160*pow3(log(pow2(Mst1)/pow2(Msq))))/3. - (2300*pow2(Mst1)*pow2(Mst2))
        /(9.*pow4(Msq)) - (4600*Dmglst1*pow2(Mst1)*pow2(Mst2))/(9.*Mgl*pow4(
        Msq)) + (130*Cbeta*MuSUSY*s2t*pow2(Mst1)*pow2(Mst2))/(3.*Mt*Sbeta*pow4(
        Msq)) + (160*z2*pow2(Mst1)*pow2(Mst2))/pow4(Msq) + (320*Dmglst1*z2*
        pow2(Mst1)*pow2(Mst2))/(Mgl*pow4(Msq)) - (2300*pow2(Dmglst1)*pow2(Mst1)
        *pow2(Mst2))/(3.*pow2(Mgl)*pow4(Msq)) + (480*z2*pow2(Dmglst1)*pow2(
        Mst1)*pow2(Mst2))/(pow2(Mgl)*pow4(Msq)) + (9723887*pow2(Mst1)*pow2(
        Mst2)*pow2(MuSUSY)*pow2(s2t))/(231525.*pow2(Mt)*pow4(Msq)) - (9723887*
        pow2(Mst1)*pow2(Mst2)*pow2(MuSUSY)*pow2(s2t))/(231525.*pow2(Mt)*pow2(
        Sbeta)*pow4(Msq)) + (4460*Cbeta*MuSUSY*pow3(Mst1))/(9.*Sbeta*pow4(Msq))
        + (28180*Cbeta*Dmglst1*MuSUSY*pow3(Mst1))/(27.*Mgl*Sbeta*pow4(Msq)) - (
        320*Cbeta*MuSUSY*z2*pow3(Mst1))/(Sbeta*pow4(Msq)) - (2240*Cbeta*
        Dmglst1*MuSUSY*z2*pow3(Mst1))/(3.*Mgl*Sbeta*pow4(Msq)) + (28180*Cbeta*
        MuSUSY*pow2(Dmglst1)*pow3(Mst1))/(27.*Sbeta*pow2(Mgl)*pow4(Msq)) - (
        2240*Cbeta*MuSUSY*z2*pow2(Dmglst1)*pow3(Mst1))/(3.*Sbeta*pow2(Mgl)*
        pow4(Msq)) - (8960*s2t*pow2(Mst2)*pow3(Mst1))/(27.*Mt*pow4(Msq)) - (
        8960*Dmglst1*s2t*pow2(Mst2)*pow3(Mst1))/(9.*Mgl*Mt*pow4(Msq)) + (640*
        s2t*z2*pow2(Mst2)*pow3(Mst1))/(3.*Mt*pow4(Msq)) + (640*Dmglst1*s2t*z2*
        pow2(Mst2)*pow3(Mst1))/(Mgl*Mt*pow4(Msq)) - (17920*s2t*pow2(Dmglst1)*
        pow2(Mst2)*pow3(Mst1))/(9.*Mt*pow2(Mgl)*pow4(Msq)) + (1280*s2t*z2*pow2(
        Dmglst1)*pow2(Mst2)*pow3(Mst1))/(Mt*pow2(Mgl)*pow4(Msq)) - (2140*s2t*
        pow2(MuSUSY)*pow3(Mst1))/(27.*Mt*pow4(Msq)) - (5020*Dmglst1*s2t*pow2(
        MuSUSY)*pow3(Mst1))/(27.*Mgl*Mt*pow4(Msq)) - (5020*s2t*pow2(Dmglst1)*
        pow2(MuSUSY)*pow3(Mst1))/(27.*Mt*pow2(Mgl)*pow4(Msq)) - (590*Cbeta*
        MuSUSY*pow2(Mst2)*pow2(s2t)*pow3(Mst1))/(9.*Sbeta*pow2(Mt)*pow4(Msq)) -
        (1310*Cbeta*Dmglst1*MuSUSY*pow2(Mst2)*pow2(s2t)*pow3(Mst1))/(9.*Mgl*
        Sbeta*pow2(Mt)*pow4(Msq)) - (1310*Cbeta*MuSUSY*pow2(Dmglst1)*pow2(Mst2)
        *pow2(s2t)*pow3(Mst1))/(9.*Sbeta*pow2(Mgl)*pow2(Mt)*pow4(Msq)) + (2140*
        s2t*pow2(MuSUSY)*pow3(Mst1))/(27.*Mt*pow2(Sbeta)*pow4(Msq)) + (5020*
        Dmglst1*s2t*pow2(MuSUSY)*pow3(Mst1))/(27.*Mgl*Mt*pow2(Sbeta)*pow4(Msq))
        + (5020*s2t*pow2(Dmglst1)*pow2(MuSUSY)*pow3(Mst1))/(27.*Mt*pow2(Mgl)*
        pow2(Sbeta)*pow4(Msq)) - (4400*pow4(Mst1))/(9.*pow2(Msq)*pow2(Mst2)) -
        (17600*Dmglst1*pow4(Mst1))/(9.*Mgl*pow2(Msq)*pow2(Mst2)) + (4400*Cbeta*
        MuSUSY*s2t*pow4(Mst1))/(9.*Mt*Sbeta*pow2(Msq)*pow2(Mst2)) + (17600*
        Cbeta*Dmglst1*MuSUSY*s2t*pow4(Mst1))/(9.*Mgl*Mt*Sbeta*pow2(Msq)*pow2(
        Mst2)) - (44000*pow2(Dmglst1)*pow4(Mst1))/(9.*pow2(Mgl)*pow2(Msq)*pow2(
        Mst2)) + (44000*Cbeta*MuSUSY*s2t*pow2(Dmglst1)*pow4(Mst1))/(9.*Mt*
        Sbeta*pow2(Mgl)*pow2(Msq)*pow2(Mst2)) - (181136*pow2(s2t)*pow4(Mst1))/(
        675.*pow2(Msq)*pow2(Mt)) - (10480*Dmglst1*pow2(s2t)*pow4(Mst1))/(9.*
        Mgl*pow2(Msq)*pow2(Mt)) - (23080*pow2(Dmglst1)*pow2(s2t)*pow4(Mst1))/(
        9.*pow2(Mgl)*pow2(Msq)*pow2(Mt)) - (29787188414*pow2(s2t)*pow4(Mst1))/(
        1.0418625e7*pow2(Mst2)*pow2(Mt)) - (253180*Dmglst1*pow2(s2t)*pow4(Mst1)
        )/(81.*Mgl*pow2(Mst2)*pow2(Mt)) + (64*OepS2*pow2(s2t)*pow4(Mst1))/(243.
        *pow2(Mst2)*pow2(Mt)) - (296*S2*pow2(s2t)*pow4(Mst1))/(9.*pow2(Mst2)*
        pow2(Mt)) - (32*T1ep*pow2(s2t)*pow4(Mst1))/(81.*pow2(Mst2)*pow2(Mt)) +
        (1180*z2*pow2(s2t)*pow4(Mst1))/(81.*pow2(Mst2)*pow2(Mt)) - (29536*
        Dmglst1*z2*pow2(s2t)*pow4(Mst1))/(9.*Mgl*pow2(Mst2)*pow2(Mt)) + (5296*
        z3*pow2(s2t)*pow4(Mst1))/(81.*pow2(Mst2)*pow2(Mt)) + (6368*Dmglst1*z3*
        pow2(s2t)*pow4(Mst1))/(3.*Mgl*pow2(Mst2)*pow2(Mt)) - (16*z4*pow2(s2t)*
        pow4(Mst1))/(81.*pow2(Mst2)*pow2(Mt)) - (270562*pow2(Dmglst1)*pow2(s2t)
        *pow4(Mst1))/(81.*pow2(Mgl)*pow2(Mst2)*pow2(Mt)) - (118864*z2*pow2(
        Dmglst1)*pow2(s2t)*pow4(Mst1))/(9.*pow2(Mgl)*pow2(Mst2)*pow2(Mt)) + (
        35024*z3*pow2(Dmglst1)*pow2(s2t)*pow4(Mst1))/(3.*pow2(Mgl)*pow2(Mst2)*
        pow2(Mt)) + (2440*pow2(MuSUSY)*pow2(s2t)*pow4(Mst1))/(9.*pow2(Msq)*
        pow2(Mst2)*pow2(Mt)) + (8800*Dmglst1*pow2(MuSUSY)*pow2(s2t)*pow4(Mst1))
        /(9.*Mgl*pow2(Msq)*pow2(Mst2)*pow2(Mt)) + (320*z2*pow2(MuSUSY)*pow2(
        s2t)*pow4(Mst1))/(3.*pow2(Msq)*pow2(Mst2)*pow2(Mt)) + (2240*Dmglst1*z2*
        pow2(MuSUSY)*pow2(s2t)*pow4(Mst1))/(3.*Mgl*pow2(Msq)*pow2(Mst2)*pow2(
        Mt)) + (26800*pow2(Dmglst1)*pow2(MuSUSY)*pow2(s2t)*pow4(Mst1))/(9.*
        pow2(Mgl)*pow2(Msq)*pow2(Mst2)*pow2(Mt)) + (6560*z2*pow2(Dmglst1)*pow2(
        MuSUSY)*pow2(s2t)*pow4(Mst1))/(3.*pow2(Mgl)*pow2(Msq)*pow2(Mst2)*pow2(
        Mt)) - (2440*pow2(MuSUSY)*pow2(s2t)*pow4(Mst1))/(9.*pow2(Msq)*pow2(
        Mst2)*pow2(Mt)*pow2(Sbeta)) - (8800*Dmglst1*pow2(MuSUSY)*pow2(s2t)*
        pow4(Mst1))/(9.*Mgl*pow2(Msq)*pow2(Mst2)*pow2(Mt)*pow2(Sbeta)) - (320*
        z2*pow2(MuSUSY)*pow2(s2t)*pow4(Mst1))/(3.*pow2(Msq)*pow2(Mst2)*pow2(Mt)
        *pow2(Sbeta)) - (2240*Dmglst1*z2*pow2(MuSUSY)*pow2(s2t)*pow4(Mst1))/(3.
        *Mgl*pow2(Msq)*pow2(Mst2)*pow2(Mt)*pow2(Sbeta)) - (26800*pow2(Dmglst1)*
        pow2(MuSUSY)*pow2(s2t)*pow4(Mst1))/(9.*pow2(Mgl)*pow2(Msq)*pow2(Mst2)*
        pow2(Mt)*pow2(Sbeta)) - (6560*z2*pow2(Dmglst1)*pow2(MuSUSY)*pow2(s2t)*
        pow4(Mst1))/(3.*pow2(Mgl)*pow2(Msq)*pow2(Mst2)*pow2(Mt)*pow2(Sbeta)) -
        (8*pow2(s2t)*pow2(z2)*pow4(Mst1))/(27.*pow2(Mst2)*pow2(Mt)) + (3148*
        Cbeta*MuSUSY*pow3(s2t)*pow4(Mst1))/(15.*Sbeta*pow2(Msq)*pow3(Mt)) + (
        2560*Cbeta*Dmglst1*MuSUSY*pow3(s2t)*pow4(Mst1))/(3.*Mgl*Sbeta*pow2(Msq)
        *pow3(Mt)) + (280*Cbeta*MuSUSY*z2*pow3(s2t)*pow4(Mst1))/(3.*Sbeta*pow2(
        Msq)*pow3(Mt)) + (1360*Cbeta*Dmglst1*MuSUSY*z2*pow3(s2t)*pow4(Mst1))/(
        3.*Mgl*Sbeta*pow2(Msq)*pow3(Mt)) + (19880*Cbeta*MuSUSY*pow2(Dmglst1)*
        pow3(s2t)*pow4(Mst1))/(9.*Sbeta*pow2(Mgl)*pow2(Msq)*pow3(Mt)) + (3640*
        Cbeta*MuSUSY*z2*pow2(Dmglst1)*pow3(s2t)*pow4(Mst1))/(3.*Sbeta*pow2(Mgl)
        *pow2(Msq)*pow3(Mt)) + (478241812219*Cbeta*MuSUSY*pow3(s2t)*pow4(Mst1))
        /(1.66698e8*Sbeta*pow2(Mst2)*pow3(Mt)) + (606133*Cbeta*Dmglst1*MuSUSY*
        pow3(s2t)*pow4(Mst1))/(324.*Mgl*Sbeta*pow2(Mst2)*pow3(Mt)) - (5312*
        Cbeta*MuSUSY*OepS2*pow3(s2t)*pow4(Mst1))/(243.*Sbeta*pow2(Mst2)*pow3(
        Mt)) + (9496*Cbeta*MuSUSY*S2*pow3(s2t)*pow4(Mst1))/(9.*Sbeta*pow2(Mst2)
        *pow3(Mt)) + (2656*Cbeta*MuSUSY*T1ep*pow3(s2t)*pow4(Mst1))/(81.*Sbeta*
        pow2(Mst2)*pow3(Mt)) + (65617*Cbeta*MuSUSY*z2*pow3(s2t)*pow4(Mst1))/(
        81.*Sbeta*pow2(Mst2)*pow3(Mt)) + (12268*Cbeta*Dmglst1*MuSUSY*z2*pow3(
        s2t)*pow4(Mst1))/(9.*Mgl*Sbeta*pow2(Mst2)*pow3(Mt)) - (330056*Cbeta*
        MuSUSY*z3*pow3(s2t)*pow4(Mst1))/(81.*Sbeta*pow2(Mst2)*pow3(Mt)) - (
        22768*Cbeta*Dmglst1*MuSUSY*z3*pow3(s2t)*pow4(Mst1))/(3.*Mgl*Sbeta*pow2(
        Mst2)*pow3(Mt)) + (1328*Cbeta*MuSUSY*z4*pow3(s2t)*pow4(Mst1))/(81.*
        Sbeta*pow2(Mst2)*pow3(Mt)) + (283217*Cbeta*MuSUSY*pow2(Dmglst1)*pow3(
        s2t)*pow4(Mst1))/(648.*Sbeta*pow2(Mgl)*pow2(Mst2)*pow3(Mt)) + (256*B4*
        Cbeta*MuSUSY*pow2(Dmglst1)*pow3(s2t)*pow4(Mst1))/(3.*Sbeta*pow2(Mgl)*
        pow2(Mst2)*pow3(Mt)) + (128*Cbeta*DN*MuSUSY*pow2(Dmglst1)*pow3(s2t)*
        pow4(Mst1))/(3.*Sbeta*pow2(Mgl)*pow2(Mst2)*pow3(Mt)) + (33886*Cbeta*
        MuSUSY*z2*pow2(Dmglst1)*pow3(s2t)*pow4(Mst1))/(9.*Sbeta*pow2(Mgl)*pow2(
        Mst2)*pow3(Mt)) - (49576*Cbeta*MuSUSY*z3*pow2(Dmglst1)*pow3(s2t)*pow4(
        Mst1))/(3.*Sbeta*pow2(Mgl)*pow2(Mst2)*pow3(Mt)) - (640*Cbeta*MuSUSY*z4*
        pow2(Dmglst1)*pow3(s2t)*pow4(Mst1))/(Sbeta*pow2(Mgl)*pow2(Mst2)*pow3(
        Mt)) + (664*Cbeta*MuSUSY*pow2(z2)*pow3(s2t)*pow4(Mst1))/(27.*Sbeta*
        pow2(Mst2)*pow3(Mt)) - (160842737*pow4(Mst1))/(231525.*pow4(Msq)) - (
        41900*Dmglst1*pow4(Mst1))/(27.*Mgl*pow4(Msq)) - (375934*Cbeta*MuSUSY*
        s2t*pow4(Mst1))/(15435.*Mt*Sbeta*pow4(Msq)) + (440*Cbeta*Dmglst1*
        MuSUSY*s2t*pow4(Mst1))/(9.*Mgl*Mt*Sbeta*pow4(Msq)) + (1520*z2*pow4(
        Mst1))/(3.*pow4(Msq)) + (4160*Dmglst1*z2*pow4(Mst1))/(3.*Mgl*pow4(Msq))
        - (430726*pow2(Dmglst1)*pow4(Mst1))/(135.*pow2(Mgl)*pow4(Msq)) + (6140*
        Cbeta*MuSUSY*s2t*pow2(Dmglst1)*pow4(Mst1))/(9.*Mt*Sbeta*pow2(Mgl)*pow4(
        Msq)) + (9440*z2*pow2(Dmglst1)*pow4(Mst1))/(3.*pow2(Mgl)*pow4(Msq)) + (
        522392*pow2(Mst2)*pow2(s2t)*pow4(Mst1))/(15435.*pow2(Mt)*pow4(Msq)) + (
        170*Dmglst1*pow2(Mst2)*pow2(s2t)*pow4(Mst1))/(9.*Mgl*pow2(Mt)*pow4(Msq)
        ) - (2485*pow2(Dmglst1)*pow2(Mst2)*pow2(s2t)*pow4(Mst1))/(9.*pow2(Mgl)*
        pow2(Mt)*pow4(Msq)) - (41220947*pow2(MuSUSY)*pow2(s2t)*pow4(Mst1))/(
        463050.*pow2(Mt)*pow4(Msq)) - (4850*Dmglst1*pow2(MuSUSY)*pow2(s2t)*
        pow4(Mst1))/(9.*Mgl*pow2(Mt)*pow4(Msq)) - (100*z2*pow2(MuSUSY)*pow2(
        s2t)*pow4(Mst1))/(3.*pow2(Mt)*pow4(Msq)) - (400*Dmglst1*z2*pow2(MuSUSY)
        *pow2(s2t)*pow4(Mst1))/(3.*Mgl*pow2(Mt)*pow4(Msq)) - (1715*pow2(
        Dmglst1)*pow2(MuSUSY)*pow2(s2t)*pow4(Mst1))/(pow2(Mgl)*pow2(Mt)*pow4(
        Msq)) - (1000*z2*pow2(Dmglst1)*pow2(MuSUSY)*pow2(s2t)*pow4(Mst1))/(3.*
        pow2(Mgl)*pow2(Mt)*pow4(Msq)) + (41220947*pow2(MuSUSY)*pow2(s2t)*pow4(
        Mst1))/(463050.*pow2(Mt)*pow2(Sbeta)*pow4(Msq)) + (4850*Dmglst1*pow2(
        MuSUSY)*pow2(s2t)*pow4(Mst1))/(9.*Mgl*pow2(Mt)*pow2(Sbeta)*pow4(Msq)) +
        (100*z2*pow2(MuSUSY)*pow2(s2t)*pow4(Mst1))/(3.*pow2(Mt)*pow2(Sbeta)*
        pow4(Msq)) + (400*Dmglst1*z2*pow2(MuSUSY)*pow2(s2t)*pow4(Mst1))/(3.*
        Mgl*pow2(Mt)*pow2(Sbeta)*pow4(Msq)) + (1715*pow2(Dmglst1)*pow2(MuSUSY)*
        pow2(s2t)*pow4(Mst1))/(pow2(Mgl)*pow2(Mt)*pow2(Sbeta)*pow4(Msq)) + (
        1000*z2*pow2(Dmglst1)*pow2(MuSUSY)*pow2(s2t)*pow4(Mst1))/(3.*pow2(Mgl)*
        pow2(Mt)*pow2(Sbeta)*pow4(Msq)) - (6740969*Cbeta*MuSUSY*pow2(Mst2)*
        pow3(s2t)*pow4(Mst1))/(102900.*Sbeta*pow3(Mt)*pow4(Msq)) - (2515*Cbeta*
        Dmglst1*MuSUSY*pow2(Mst2)*pow3(s2t)*pow4(Mst1))/(9.*Mgl*Sbeta*pow3(Mt)*
        pow4(Msq)) - (50*Cbeta*MuSUSY*z2*pow2(Mst2)*pow3(s2t)*pow4(Mst1))/(3.*
        Sbeta*pow3(Mt)*pow4(Msq)) - (200*Cbeta*Dmglst1*MuSUSY*z2*pow2(Mst2)*
        pow3(s2t)*pow4(Mst1))/(3.*Mgl*Sbeta*pow3(Mt)*pow4(Msq)) - (1745*Cbeta*
        MuSUSY*pow2(Dmglst1)*pow2(Mst2)*pow3(s2t)*pow4(Mst1))/(2.*Sbeta*pow2(
        Mgl)*pow3(Mt)*pow4(Msq)) - (500*Cbeta*MuSUSY*z2*pow2(Dmglst1)*pow2(
        Mst2)*pow3(s2t)*pow4(Mst1))/(3.*Sbeta*pow2(Mgl)*pow3(Mt)*pow4(Msq)) - (
        8896*pow2(Mst1)*pow2(MuSUSY))/(3.*pow4(Mst2)) + (128*B4*pow2(Mst1)*
        pow2(MuSUSY))/(3.*pow4(Mst2)) + (64*DN*pow2(Mst1)*pow2(MuSUSY))/(3.*
        pow4(Mst2)) - (7040*z2*pow2(Mst1)*pow2(MuSUSY))/(3.*pow4(Mst2)) - (
        2944*z3*pow2(Mst1)*pow2(MuSUSY))/(3.*pow4(Mst2)) - (320*z4*pow2(Mst1)*
        pow2(MuSUSY))/pow4(Mst2) + (8896*pow2(Mst1)*pow2(MuSUSY))/(3.*pow2(
        Sbeta)*pow4(Mst2)) - (128*B4*pow2(Mst1)*pow2(MuSUSY))/(3.*pow2(Sbeta)*
        pow4(Mst2)) - (64*DN*pow2(Mst1)*pow2(MuSUSY))/(3.*pow2(Sbeta)*pow4(
        Mst2)) + (7040*z2*pow2(Mst1)*pow2(MuSUSY))/(3.*pow2(Sbeta)*pow4(Mst2))
        + (2944*z3*pow2(Mst1)*pow2(MuSUSY))/(3.*pow2(Sbeta)*pow4(Mst2)) + (320*
        z4*pow2(Mst1)*pow2(MuSUSY))/(pow2(Sbeta)*pow4(Mst2)) + (27952*Cbeta*
        MuSUSY*pow3(Mst1))/(27.*Sbeta*pow4(Mst2)) - (260848*Cbeta*Dmglst1*
        MuSUSY*pow3(Mst1))/(81.*Mgl*Sbeta*pow4(Mst2)) - (88384*Cbeta*MuSUSY*z2*
        pow3(Mst1))/(9.*Sbeta*pow4(Mst2)) - (266432*Cbeta*Dmglst1*MuSUSY*z2*
        pow3(Mst1))/(9.*Mgl*Sbeta*pow4(Mst2)) + (2752*Cbeta*MuSUSY*z3*pow3(
        Mst1))/(Sbeta*pow4(Mst2)) + (30080*Cbeta*Dmglst1*MuSUSY*z3*pow3(Mst1))/
        (3.*Mgl*Sbeta*pow4(Mst2)) + (128*Cbeta*MuSUSY*z4*pow3(Mst1))/(3.*Sbeta*
        pow4(Mst2)) + (128*Cbeta*Dmglst1*MuSUSY*z4*pow3(Mst1))/(3.*Mgl*Sbeta*
        pow4(Mst2)) - (260848*Cbeta*MuSUSY*pow2(Dmglst1)*pow3(Mst1))/(81.*
        Sbeta*pow2(Mgl)*pow4(Mst2)) - (266432*Cbeta*MuSUSY*z2*pow2(Dmglst1)*
        pow3(Mst1))/(9.*Sbeta*pow2(Mgl)*pow4(Mst2)) + (30080*Cbeta*MuSUSY*z3*
        pow2(Dmglst1)*pow3(Mst1))/(3.*Sbeta*pow2(Mgl)*pow4(Mst2)) + (128*Cbeta*
        MuSUSY*z4*pow2(Dmglst1)*pow3(Mst1))/(3.*Sbeta*pow2(Mgl)*pow4(Mst2)) - (
        858392*s2t*pow2(MuSUSY)*pow3(Mst1))/(81.*Mt*pow4(Mst2)) + (5696*B4*s2t*
        pow2(MuSUSY)*pow3(Mst1))/(3.*Mt*pow4(Mst2)) + (32*DN*s2t*pow2(MuSUSY)*
        pow3(Mst1))/(3.*Mt*pow4(Mst2)) - (111752*Dmglst1*s2t*pow2(MuSUSY)*pow3(
        Mst1))/(27.*Mgl*Mt*pow4(Mst2)) + (3520*B4*Dmglst1*s2t*pow2(MuSUSY)*
        pow3(Mst1))/(Mgl*Mt*pow4(Mst2)) - (32*Dmglst1*DN*s2t*pow2(MuSUSY)*pow3(
        Mst1))/(Mgl*Mt*pow4(Mst2)) - (77024*s2t*z2*pow2(MuSUSY)*pow3(Mst1))/(9.
        *Mt*pow4(Mst2)) + (7328*Dmglst1*s2t*z2*pow2(MuSUSY)*pow3(Mst1))/(Mgl*
        Mt*pow4(Mst2)) + (36800*s2t*z3*pow2(MuSUSY)*pow3(Mst1))/(3.*Mt*pow4(
        Mst2)) + (9856*Dmglst1*s2t*z3*pow2(MuSUSY)*pow3(Mst1))/(Mgl*Mt*pow4(
        Mst2)) + (10016*s2t*z4*pow2(MuSUSY)*pow3(Mst1))/(3.*Mt*pow4(Mst2)) - (
        5408*Dmglst1*s2t*z4*pow2(MuSUSY)*pow3(Mst1))/(Mgl*Mt*pow4(Mst2)) - (
        111752*s2t*pow2(Dmglst1)*pow2(MuSUSY)*pow3(Mst1))/(27.*Mt*pow2(Mgl)*
        pow4(Mst2)) + (3520*B4*s2t*pow2(Dmglst1)*pow2(MuSUSY)*pow3(Mst1))/(Mt*
        pow2(Mgl)*pow4(Mst2)) - (32*DN*s2t*pow2(Dmglst1)*pow2(MuSUSY)*pow3(
        Mst1))/(Mt*pow2(Mgl)*pow4(Mst2)) + (7328*s2t*z2*pow2(Dmglst1)*pow2(
        MuSUSY)*pow3(Mst1))/(Mt*pow2(Mgl)*pow4(Mst2)) + (9856*s2t*z3*pow2(
        Dmglst1)*pow2(MuSUSY)*pow3(Mst1))/(Mt*pow2(Mgl)*pow4(Mst2)) - (5408*
        s2t*z4*pow2(Dmglst1)*pow2(MuSUSY)*pow3(Mst1))/(Mt*pow2(Mgl)*pow4(Mst2))
        + (858392*s2t*pow2(MuSUSY)*pow3(Mst1))/(81.*Mt*pow2(Sbeta)*pow4(Mst2))
        - (5696*B4*s2t*pow2(MuSUSY)*pow3(Mst1))/(3.*Mt*pow2(Sbeta)*pow4(Mst2))
        - (32*DN*s2t*pow2(MuSUSY)*pow3(Mst1))/(3.*Mt*pow2(Sbeta)*pow4(Mst2)) +
        (111752*Dmglst1*s2t*pow2(MuSUSY)*pow3(Mst1))/(27.*Mgl*Mt*pow2(Sbeta)*
        pow4(Mst2)) - (3520*B4*Dmglst1*s2t*pow2(MuSUSY)*pow3(Mst1))/(Mgl*Mt*
        pow2(Sbeta)*pow4(Mst2)) + (32*Dmglst1*DN*s2t*pow2(MuSUSY)*pow3(Mst1))/(
        Mgl*Mt*pow2(Sbeta)*pow4(Mst2)) + (77024*s2t*z2*pow2(MuSUSY)*pow3(Mst1))
        /(9.*Mt*pow2(Sbeta)*pow4(Mst2)) - (7328*Dmglst1*s2t*z2*pow2(MuSUSY)*
        pow3(Mst1))/(Mgl*Mt*pow2(Sbeta)*pow4(Mst2)) - (36800*s2t*z3*pow2(
        MuSUSY)*pow3(Mst1))/(3.*Mt*pow2(Sbeta)*pow4(Mst2)) - (9856*Dmglst1*s2t*
        z3*pow2(MuSUSY)*pow3(Mst1))/(Mgl*Mt*pow2(Sbeta)*pow4(Mst2)) - (10016*
        s2t*z4*pow2(MuSUSY)*pow3(Mst1))/(3.*Mt*pow2(Sbeta)*pow4(Mst2)) + (5408*
        Dmglst1*s2t*z4*pow2(MuSUSY)*pow3(Mst1))/(Mgl*Mt*pow2(Sbeta)*pow4(Mst2))
        + (111752*s2t*pow2(Dmglst1)*pow2(MuSUSY)*pow3(Mst1))/(27.*Mt*pow2(Mgl)*
        pow2(Sbeta)*pow4(Mst2)) - (3520*B4*s2t*pow2(Dmglst1)*pow2(MuSUSY)*pow3(
        Mst1))/(Mt*pow2(Mgl)*pow2(Sbeta)*pow4(Mst2)) + (32*DN*s2t*pow2(Dmglst1)
        *pow2(MuSUSY)*pow3(Mst1))/(Mt*pow2(Mgl)*pow2(Sbeta)*pow4(Mst2)) - (
        7328*s2t*z2*pow2(Dmglst1)*pow2(MuSUSY)*pow3(Mst1))/(Mt*pow2(Mgl)*pow2(
        Sbeta)*pow4(Mst2)) - (9856*s2t*z3*pow2(Dmglst1)*pow2(MuSUSY)*pow3(Mst1)
        )/(Mt*pow2(Mgl)*pow2(Sbeta)*pow4(Mst2)) + (5408*s2t*z4*pow2(Dmglst1)*
        pow2(MuSUSY)*pow3(Mst1))/(Mt*pow2(Mgl)*pow2(Sbeta)*pow4(Mst2)) - (256*
        Cbeta*MuSUSY*pow2(z2)*pow3(Mst1))/(3.*Sbeta*pow4(Mst2)) - (256*Cbeta*
        Dmglst1*MuSUSY*pow2(z2)*pow3(Mst1))/(3.*Mgl*Sbeta*pow4(Mst2)) - (256*
        Cbeta*MuSUSY*pow2(Dmglst1)*pow2(z2)*pow3(Mst1))/(3.*Sbeta*pow2(Mgl)*
        pow4(Mst2)) - (5632*s2t*pow2(MuSUSY)*pow2(z2)*pow3(Mst1))/(3.*Mt*pow4(
        Mst2)) - (3584*Dmglst1*s2t*pow2(MuSUSY)*pow2(z2)*pow3(Mst1))/(Mgl*Mt*
        pow4(Mst2)) - (3584*s2t*pow2(Dmglst1)*pow2(MuSUSY)*pow2(z2)*pow3(Mst1))
        /(Mt*pow2(Mgl)*pow4(Mst2)) + (5632*s2t*pow2(MuSUSY)*pow2(z2)*pow3(Mst1)
        )/(3.*Mt*pow2(Sbeta)*pow4(Mst2)) + (3584*Dmglst1*s2t*pow2(MuSUSY)*pow2(
        z2)*pow3(Mst1))/(Mgl*Mt*pow2(Sbeta)*pow4(Mst2)) + (3584*s2t*pow2(
        Dmglst1)*pow2(MuSUSY)*pow2(z2)*pow3(Mst1))/(Mt*pow2(Mgl)*pow2(Sbeta)*
        pow4(Mst2)) - (25346116504*pow4(Mst1))/(3.472875e6*pow4(Mst2)) - (
        799880*Dmglst1*pow4(Mst1))/(27.*Mgl*pow4(Mst2)) - (5244811004*Cbeta*
        MuSUSY*s2t*pow4(Mst1))/(1.0418625e7*Mt*Sbeta*pow4(Mst2)) - (284528*
        Cbeta*Dmglst1*MuSUSY*s2t*pow4(Mst1))/(27.*Mgl*Mt*Sbeta*pow4(Mst2)) + (
        256*B4*Cbeta*Dmglst1*MuSUSY*s2t*pow4(Mst1))/(3.*Mgl*Mt*Sbeta*pow4(Mst2)
        ) + (128*Cbeta*Dmglst1*DN*MuSUSY*s2t*pow4(Mst1))/(3.*Mgl*Mt*Sbeta*pow4(
        Mst2)) + (2560*Cbeta*MuSUSY*OepS2*s2t*pow4(Mst1))/(243.*Mt*Sbeta*pow4(
        Mst2)) - (4448*Cbeta*MuSUSY*S2*s2t*pow4(Mst1))/(9.*Mt*Sbeta*pow4(Mst2))
        - (1280*Cbeta*MuSUSY*s2t*T1ep*pow4(Mst1))/(81.*Mt*Sbeta*pow4(Mst2)) - (
        38128*z2*pow4(Mst1))/(9.*pow4(Mst2)) - (19904*Dmglst1*z2*pow4(Mst1))/(
        Mgl*pow4(Mst2)) - (422384*Cbeta*MuSUSY*s2t*z2*pow4(Mst1))/(81.*Mt*
        Sbeta*pow4(Mst2)) - (29792*Cbeta*Dmglst1*MuSUSY*s2t*z2*pow4(Mst1))/(9.*
        Mgl*Mt*Sbeta*pow4(Mst2)) + (10144*z3*pow4(Mst1))/(3.*pow4(Mst2)) + (
        42112*Dmglst1*z3*pow4(Mst1))/(3.*Mgl*pow4(Mst2)) + (23488*Cbeta*MuSUSY*
        s2t*z3*pow4(Mst1))/(81.*Mt*Sbeta*pow4(Mst2)) - (2560*Cbeta*Dmglst1*
        MuSUSY*s2t*z3*pow4(Mst1))/(3.*Mgl*Mt*Sbeta*pow4(Mst2)) - (640*Cbeta*
        MuSUSY*s2t*z4*pow4(Mst1))/(81.*Mt*Sbeta*pow4(Mst2)) - (640*Cbeta*
        Dmglst1*MuSUSY*s2t*z4*pow4(Mst1))/(Mgl*Mt*Sbeta*pow4(Mst2)) - (
        150456364*pow2(Dmglst1)*pow4(Mst1))/(2025.*pow2(Mgl)*pow4(Mst2)) - (
        2669512*Cbeta*MuSUSY*s2t*pow2(Dmglst1)*pow4(Mst1))/(81.*Mt*Sbeta*pow2(
        Mgl)*pow4(Mst2)) + (128*B4*Cbeta*MuSUSY*s2t*pow2(Dmglst1)*pow4(Mst1))/(
        Mt*Sbeta*pow2(Mgl)*pow4(Mst2)) + (64*Cbeta*DN*MuSUSY*s2t*pow2(Dmglst1)*
        pow4(Mst1))/(Mt*Sbeta*pow2(Mgl)*pow4(Mst2)) - (53696*z2*pow2(Dmglst1)*
        pow4(Mst1))/(pow2(Mgl)*pow4(Mst2)) + (104432*Cbeta*MuSUSY*s2t*z2*pow2(
        Dmglst1)*pow4(Mst1))/(9.*Mt*Sbeta*pow2(Mgl)*pow4(Mst2)) + (35648*z3*
        pow2(Dmglst1)*pow4(Mst1))/(pow2(Mgl)*pow4(Mst2)) - (19712*Cbeta*MuSUSY*
        s2t*z3*pow2(Dmglst1)*pow4(Mst1))/(3.*Mt*Sbeta*pow2(Mgl)*pow4(Mst2)) - (
        960*Cbeta*MuSUSY*s2t*z4*pow2(Dmglst1)*pow4(Mst1))/(Mt*Sbeta*pow2(Mgl)*
        pow4(Mst2)) + (531403689547*pow2(MuSUSY)*pow2(s2t)*pow4(Mst1))/(8.3349e7
        *pow2(Mt)*pow4(Mst2)) - (608*B4*pow2(MuSUSY)*pow2(s2t)*pow4(Mst1))/(
        3.*pow2(Mt)*pow4(Mst2)) + (16*DN*pow2(MuSUSY)*pow2(s2t)*pow4(Mst1))/(3.
        *pow2(Mt)*pow4(Mst2)) + (437365*Dmglst1*pow2(MuSUSY)*pow2(s2t)*pow4(
        Mst1))/(162.*Mgl*pow2(Mt)*pow4(Mst2)) - (64*B4*Dmglst1*pow2(MuSUSY)*
        pow2(s2t)*pow4(Mst1))/(Mgl*pow2(Mt)*pow4(Mst2)) - (160*Dmglst1*DN*pow2(
        MuSUSY)*pow2(s2t)*pow4(Mst1))/(3.*Mgl*pow2(Mt)*pow4(Mst2)) - (23680*
        OepS2*pow2(MuSUSY)*pow2(s2t)*pow4(Mst1))/(243.*pow2(Mt)*pow4(Mst2)) + (
        33104*S2*pow2(MuSUSY)*pow2(s2t)*pow4(Mst1))/(9.*pow2(Mt)*pow4(Mst2)) +
        (11840*T1ep*pow2(MuSUSY)*pow2(s2t)*pow4(Mst1))/(81.*pow2(Mt)*pow4(Mst2)
        ) + (138698*z2*pow2(MuSUSY)*pow2(s2t)*pow4(Mst1))/(81.*pow2(Mt)*pow4(
        Mst2)) - (12904*Dmglst1*z2*pow2(MuSUSY)*pow2(s2t)*pow4(Mst1))/(9.*Mgl*
        pow2(Mt)*pow4(Mst2)) - (1081696*z3*pow2(MuSUSY)*pow2(s2t)*pow4(Mst1))/(
        81.*pow2(Mt)*pow4(Mst2)) - (48896*Dmglst1*z3*pow2(MuSUSY)*pow2(s2t)*
        pow4(Mst1))/(3.*Mgl*pow2(Mt)*pow4(Mst2)) + (2032*z4*pow2(MuSUSY)*pow2(
        s2t)*pow4(Mst1))/(81.*pow2(Mt)*pow4(Mst2)) + (8480*Dmglst1*z4*pow2(
        MuSUSY)*pow2(s2t)*pow4(Mst1))/(3.*Mgl*pow2(Mt)*pow4(Mst2)) - (130639*
        pow2(Dmglst1)*pow2(MuSUSY)*pow2(s2t)*pow4(Mst1))/(324.*pow2(Mgl)*pow2(
        Mt)*pow4(Mst2)) + (224*B4*pow2(Dmglst1)*pow2(MuSUSY)*pow2(s2t)*pow4(
        Mst1))/(3.*pow2(Mgl)*pow2(Mt)*pow4(Mst2)) + (16*DN*pow2(Dmglst1)*pow2(
        MuSUSY)*pow2(s2t)*pow4(Mst1))/(3.*pow2(Mgl)*pow2(Mt)*pow4(Mst2)) + (
        11612*z2*pow2(Dmglst1)*pow2(MuSUSY)*pow2(s2t)*pow4(Mst1))/(9.*pow2(Mgl)
        *pow2(Mt)*pow4(Mst2)) - (103856*z3*pow2(Dmglst1)*pow2(MuSUSY)*pow2(s2t)
        *pow4(Mst1))/(3.*pow2(Mgl)*pow2(Mt)*pow4(Mst2)) + (2960*z4*pow2(
        Dmglst1)*pow2(MuSUSY)*pow2(s2t)*pow4(Mst1))/(pow2(Mgl)*pow2(Mt)*pow4(
        Mst2)) - (531403689547*pow2(MuSUSY)*pow2(s2t)*pow4(Mst1))/(8.3349e7*
        pow2(Mt)*pow2(Sbeta)*pow4(Mst2)) + (608*B4*pow2(MuSUSY)*pow2(s2t)*pow4(
        Mst1))/(3.*pow2(Mt)*pow2(Sbeta)*pow4(Mst2)) - (16*DN*pow2(MuSUSY)*pow2(
        s2t)*pow4(Mst1))/(3.*pow2(Mt)*pow2(Sbeta)*pow4(Mst2)) - (437365*
        Dmglst1*pow2(MuSUSY)*pow2(s2t)*pow4(Mst1))/(162.*Mgl*pow2(Mt)*pow2(
        Sbeta)*pow4(Mst2)) + (64*B4*Dmglst1*pow2(MuSUSY)*pow2(s2t)*pow4(Mst1))/
        (Mgl*pow2(Mt)*pow2(Sbeta)*pow4(Mst2)) + (160*Dmglst1*DN*pow2(MuSUSY)*
        pow2(s2t)*pow4(Mst1))/(3.*Mgl*pow2(Mt)*pow2(Sbeta)*pow4(Mst2)) + (
        23680*OepS2*pow2(MuSUSY)*pow2(s2t)*pow4(Mst1))/(243.*pow2(Mt)*pow2(
        Sbeta)*pow4(Mst2)) - (33104*S2*pow2(MuSUSY)*pow2(s2t)*pow4(Mst1))/(9.*
        pow2(Mt)*pow2(Sbeta)*pow4(Mst2)) - (11840*T1ep*pow2(MuSUSY)*pow2(s2t)*
        pow4(Mst1))/(81.*pow2(Mt)*pow2(Sbeta)*pow4(Mst2)) - (138698*z2*pow2(
        MuSUSY)*pow2(s2t)*pow4(Mst1))/(81.*pow2(Mt)*pow2(Sbeta)*pow4(Mst2)) + (
        12904*Dmglst1*z2*pow2(MuSUSY)*pow2(s2t)*pow4(Mst1))/(9.*Mgl*pow2(Mt)*
        pow2(Sbeta)*pow4(Mst2)) + (1081696*z3*pow2(MuSUSY)*pow2(s2t)*pow4(Mst1)
        )/(81.*pow2(Mt)*pow2(Sbeta)*pow4(Mst2)) + (48896*Dmglst1*z3*pow2(
        MuSUSY)*pow2(s2t)*pow4(Mst1))/(3.*Mgl*pow2(Mt)*pow2(Sbeta)*pow4(Mst2))
        - (2032*z4*pow2(MuSUSY)*pow2(s2t)*pow4(Mst1))/(81.*pow2(Mt)*pow2(Sbeta)
        *pow4(Mst2)) - (8480*Dmglst1*z4*pow2(MuSUSY)*pow2(s2t)*pow4(Mst1))/(3.*
        Mgl*pow2(Mt)*pow2(Sbeta)*pow4(Mst2)) + (130639*pow2(Dmglst1)*pow2(
        MuSUSY)*pow2(s2t)*pow4(Mst1))/(324.*pow2(Mgl)*pow2(Mt)*pow2(Sbeta)*
        pow4(Mst2)) - (224*B4*pow2(Dmglst1)*pow2(MuSUSY)*pow2(s2t)*pow4(Mst1))/
        (3.*pow2(Mgl)*pow2(Mt)*pow2(Sbeta)*pow4(Mst2)) - (16*DN*pow2(Dmglst1)*
        pow2(MuSUSY)*pow2(s2t)*pow4(Mst1))/(3.*pow2(Mgl)*pow2(Mt)*pow2(Sbeta)*
        pow4(Mst2)) - (11612*z2*pow2(Dmglst1)*pow2(MuSUSY)*pow2(s2t)*pow4(Mst1)
        )/(9.*pow2(Mgl)*pow2(Mt)*pow2(Sbeta)*pow4(Mst2)) + (103856*z3*pow2(
        Dmglst1)*pow2(MuSUSY)*pow2(s2t)*pow4(Mst1))/(3.*pow2(Mgl)*pow2(Mt)*
        pow2(Sbeta)*pow4(Mst2)) - (2960*z4*pow2(Dmglst1)*pow2(MuSUSY)*pow2(s2t)
        *pow4(Mst1))/(pow2(Mgl)*pow2(Mt)*pow2(Sbeta)*pow4(Mst2)) - (320*Cbeta*
        MuSUSY*s2t*pow2(z2)*pow4(Mst1))/(27.*Mt*Sbeta*pow4(Mst2)) + (8144*pow2(
        MuSUSY)*pow2(s2t)*pow2(z2)*pow4(Mst1))/(27.*pow2(Mt)*pow4(Mst2)) + (
        128*Dmglst1*pow2(MuSUSY)*pow2(s2t)*pow2(z2)*pow4(Mst1))/(3.*Mgl*pow2(
        Mt)*pow4(Mst2)) + (64*pow2(Dmglst1)*pow2(MuSUSY)*pow2(s2t)*pow2(z2)*
        pow4(Mst1))/(pow2(Mgl)*pow2(Mt)*pow4(Mst2)) - (8144*pow2(MuSUSY)*pow2(
        s2t)*pow2(z2)*pow4(Mst1))/(27.*pow2(Mt)*pow2(Sbeta)*pow4(Mst2)) - (128*
        Dmglst1*pow2(MuSUSY)*pow2(s2t)*pow2(z2)*pow4(Mst1))/(3.*Mgl*pow2(Mt)*
        pow2(Sbeta)*pow4(Mst2)) - (64*pow2(Dmglst1)*pow2(MuSUSY)*pow2(s2t)*
        pow2(z2)*pow4(Mst1))/(pow2(Mgl)*pow2(Mt)*pow2(Sbeta)*pow4(Mst2)) + (
        40864*pow2(s2t)*pow4(Mst2))/(675.*pow2(Msq)*pow2(Mt)) - (8*pow2(s2t)*
        pow4(Mst2))/(pow2(Mst1)*pow2(Mt)) + (16*z2*pow2(s2t)*pow4(Mst2))/(pow2(
        Mst1)*pow2(Mt)) + (80*Mst1*pow3(s2t)*pow4(Mst2))/(3.*pow2(Msq)*pow3(Mt)
        ) + (80*Dmglst1*Mst1*pow3(s2t)*pow4(Mst2))/(3.*Mgl*pow2(Msq)*pow3(Mt))
        + (1252*Cbeta*MuSUSY*pow3(s2t)*pow4(Mst2))/(15.*Sbeta*pow2(Msq)*pow3(
        Mt)) + (80*Mst1*pow2(Dmglst1)*pow3(s2t)*pow4(Mst2))/(3.*pow2(Mgl)*pow2(
        Msq)*pow3(Mt)) - (4*Cbeta*MuSUSY*pow3(s2t)*pow4(Mst2))/(Sbeta*pow2(
        Mst1)*pow3(Mt)) + (8*Cbeta*MuSUSY*z2*pow3(s2t)*pow4(Mst2))/(Sbeta*pow2(
        Mst1)*pow3(Mt)) - (13597579*pow4(Mst2))/(77175.*pow4(Msq)) - (6760*
        Mst1*s2t*pow4(Mst2))/(27.*Mt*pow4(Msq)) - (6760*Dmglst1*Mst1*s2t*pow4(
        Mst2))/(27.*Mgl*Mt*pow4(Msq)) - (207166*Cbeta*MuSUSY*s2t*pow4(Mst2))/(
        15435.*Mt*Sbeta*pow4(Msq)) + (80*z2*pow4(Mst2))/pow4(Msq) + (320*Mst1*
        s2t*z2*pow4(Mst2))/(3.*Mt*pow4(Msq)) + (320*Dmglst1*Mst1*s2t*z2*pow4(
        Mst2))/(3.*Mgl*Mt*pow4(Msq)) - (6760*Mst1*s2t*pow2(Dmglst1)*pow4(Mst2))
        /(27.*Mt*pow2(Mgl)*pow4(Msq)) + (320*Mst1*s2t*z2*pow2(Dmglst1)*pow4(
        Mst2))/(3.*Mt*pow2(Mgl)*pow4(Msq)) + (55*Cbeta*Mst1*MuSUSY*pow2(s2t)*
        pow4(Mst2))/(9.*Sbeta*pow2(Mt)*pow4(Msq)) + (55*Cbeta*Dmglst1*Mst1*
        MuSUSY*pow2(s2t)*pow4(Mst2))/(9.*Mgl*Sbeta*pow2(Mt)*pow4(Msq)) + (55*
        Cbeta*Mst1*MuSUSY*pow2(Dmglst1)*pow2(s2t)*pow4(Mst2))/(9.*Sbeta*pow2(
        Mgl)*pow2(Mt)*pow4(Msq)) - (438008*pow2(Mst1)*pow2(s2t)*pow4(Mst2))/(
        15435.*pow2(Mt)*pow4(Msq)) - (130*Dmglst1*pow2(Mst1)*pow2(s2t)*pow4(
        Mst2))/(3.*Mgl*pow2(Mt)*pow4(Msq)) - (65*pow2(Dmglst1)*pow2(Mst1)*pow2(
        s2t)*pow4(Mst2))/(pow2(Mgl)*pow2(Mt)*pow4(Msq)) + (7399303*pow2(MuSUSY)
        *pow2(s2t)*pow4(Mst2))/(463050.*pow2(Mt)*pow4(Msq)) - (7399303*pow2(
        MuSUSY)*pow2(s2t)*pow4(Mst2))/(463050.*pow2(Mt)*pow2(Sbeta)*pow4(Msq))
        + (1338719*Cbeta*MuSUSY*pow2(Mst1)*pow3(s2t)*pow4(Mst2))/(102900.*
        Sbeta*pow3(Mt)*pow4(Msq)) + (10*Cbeta*Dmglst1*MuSUSY*pow2(Mst1)*pow3(
        s2t)*pow4(Mst2))/(Mgl*Sbeta*pow3(Mt)*pow4(Msq)) + (15*Cbeta*MuSUSY*
        pow2(Dmglst1)*pow2(Mst1)*pow3(s2t)*pow4(Mst2))/(Sbeta*pow2(Mgl)*pow3(
        Mt)*pow4(Msq)) + (215*pow3(Mst1)*pow3(s2t)*pow4(Mst2))/(9.*pow3(Mt)*
        pow4(Msq)) + (455*Dmglst1*pow3(Mst1)*pow3(s2t)*pow4(Mst2))/(9.*Mgl*
        pow3(Mt)*pow4(Msq)) + (815*pow2(Dmglst1)*pow3(Mst1)*pow3(s2t)*pow4(
        Mst2))/(9.*pow2(Mgl)*pow3(Mt)*pow4(Msq)) - (10*pow2(Msq)*pow2(Mst1)*
        pow4(s2t))/pow4(Mt) + (20*z2*pow2(Msq)*pow2(Mst1)*pow4(s2t))/pow4(Mt) -
        (10*pow2(Msq)*pow2(Mst2)*pow4(s2t))/pow4(Mt) + (20*z2*pow2(Msq)*pow2(
        Mst2)*pow4(s2t))/pow4(Mt) - (5876588*pow2(Mst1)*pow2(Mst2)*pow4(s2t))/(
        10125.*pow4(Mt)) + (4*B4*pow2(Mst1)*pow2(Mst2)*pow4(s2t))/(3.*pow4(Mt))
        + (2*DN*pow2(Mst1)*pow2(Mst2)*pow4(s2t))/pow4(Mt) + (1151*Dmglst1*pow2(
        Mst1)*pow2(Mst2)*pow4(s2t))/(9.*Mgl*pow4(Mt)) + (8*B4*Dmglst1*pow2(
        Mst1)*pow2(Mst2)*pow4(s2t))/(Mgl*pow4(Mt)) + (20*Dmglst1*DN*pow2(Mst1)*
        pow2(Mst2)*pow4(s2t))/(3.*Mgl*pow4(Mt)) + (400*OepS2*pow2(Mst1)*pow2(
        Mst2)*pow4(s2t))/(81.*pow4(Mt)) - (394*S2*pow2(Mst1)*pow2(Mst2)*pow4(
        s2t))/pow4(Mt) - (200*T1ep*pow2(Mst1)*pow2(Mst2)*pow4(s2t))/(27.*pow4(
        Mt)) + (5737*z2*pow2(Mst1)*pow2(Mst2)*pow4(s2t))/(27.*pow4(Mt)) + (388*
        Dmglst1*z2*pow2(Mst1)*pow2(Mst2)*pow4(s2t))/(Mgl*pow4(Mt)) + (11662*z3*
        pow2(Mst1)*pow2(Mst2)*pow4(s2t))/(27.*pow4(Mt)) + (270*Dmglst1*z3*pow2(
        Mst1)*pow2(Mst2)*pow4(s2t))/(Mgl*pow4(Mt)) - (4330*z4*pow2(Mst1)*pow2(
        Mst2)*pow4(s2t))/(27.*pow4(Mt)) - (1060*Dmglst1*z4*pow2(Mst1)*pow2(
        Mst2)*pow4(s2t))/(3.*Mgl*pow4(Mt)) - (3583*pow2(Dmglst1)*pow2(Mst1)*
        pow2(Mst2)*pow4(s2t))/(18.*pow2(Mgl)*pow4(Mt)) + (12*B4*pow2(Dmglst1)*
        pow2(Mst1)*pow2(Mst2)*pow4(s2t))/(pow2(Mgl)*pow4(Mt)) + (10*DN*pow2(
        Dmglst1)*pow2(Mst1)*pow2(Mst2)*pow4(s2t))/(pow2(Mgl)*pow4(Mt)) + (582*
        z2*pow2(Dmglst1)*pow2(Mst1)*pow2(Mst2)*pow4(s2t))/(pow2(Mgl)*pow4(Mt))
        + (707*z3*pow2(Dmglst1)*pow2(Mst1)*pow2(Mst2)*pow4(s2t))/(pow2(Mgl)*
        pow4(Mt)) - (530*z4*pow2(Dmglst1)*pow2(Mst1)*pow2(Mst2)*pow4(s2t))/(
        pow2(Mgl)*pow4(Mt)) - (74*pow2(Mst1)*pow2(Mst2)*pow2(z2)*pow4(s2t))/(9.
        *pow4(Mt)) - (16*Dmglst1*pow2(Mst1)*pow2(Mst2)*pow2(z2)*pow4(s2t))/(3.*
        Mgl*pow4(Mt)) - (8*pow2(Dmglst1)*pow2(Mst1)*pow2(Mst2)*pow2(z2)*pow4(
        s2t))/(pow2(Mgl)*pow4(Mt)) - (257823187891*pow4(Mst1)*pow4(s2t))/(
        6.66792e8*pow4(Mt)) - (40*B4*pow4(Mst1)*pow4(s2t))/(3.*pow4(Mt)) - (2*
        DN*pow4(Mst1)*pow4(s2t))/(3.*pow4(Mt)) - (773389*Dmglst1*pow4(Mst1)*
        pow4(s2t))/(1296.*Mgl*pow4(Mt)) - (8*B4*Dmglst1*pow4(Mst1)*pow4(s2t))/(
        Mgl*pow4(Mt)) - (20*Dmglst1*DN*pow4(Mst1)*pow4(s2t))/(3.*Mgl*pow4(Mt))
        - (88*OepS2*pow4(Mst1)*pow4(s2t))/(243.*pow4(Mt)) + (281*S2*pow4(Mst1)*
        pow4(s2t))/(9.*pow4(Mt)) + (44*T1ep*pow4(Mst1)*pow4(s2t))/(81.*pow4(Mt)
        ) - (98497*z2*pow4(Mst1)*pow4(s2t))/(324.*pow4(Mt)) - (7153*Dmglst1*z2*
        pow4(Mst1)*pow4(s2t))/(9.*Mgl*pow4(Mt)) + (38672*z3*pow4(Mst1)*pow4(
        s2t))/(81.*pow4(Mt)) + (5077*Dmglst1*z3*pow4(Mst1)*pow4(s2t))/(3.*Mgl*
        pow4(Mt)) + (6124*z4*pow4(Mst1)*pow4(s2t))/(81.*pow4(Mt)) + (1060*
        Dmglst1*z4*pow4(Mst1)*pow4(s2t))/(3.*Mgl*pow4(Mt)) - (232169*pow2(
        Dmglst1)*pow4(Mst1)*pow4(s2t))/(2592.*pow2(Mgl)*pow4(Mt)) - (100*B4*
        pow2(Dmglst1)*pow4(Mst1)*pow4(s2t))/(3.*pow2(Mgl)*pow4(Mt)) - (62*DN*
        pow2(Dmglst1)*pow4(Mst1)*pow4(s2t))/(3.*pow2(Mgl)*pow4(Mt)) - (29201*
        z2*pow2(Dmglst1)*pow4(Mst1)*pow4(s2t))/(18.*pow2(Mgl)*pow4(Mt)) + (
        22079*z3*pow2(Dmglst1)*pow4(Mst1)*pow4(s2t))/(6.*pow2(Mgl)*pow4(Mt)) +
        (690*z4*pow2(Dmglst1)*pow4(Mst1)*pow4(s2t))/(pow2(Mgl)*pow4(Mt)) + (10*
        pow2(Msq)*pow4(Mst1)*pow4(s2t))/(pow2(Mst2)*pow4(Mt)) - (20*z2*pow2(
        Msq)*pow4(Mst1)*pow4(s2t))/(pow2(Mst2)*pow4(Mt)) - (4136*pow2(Mst2)*
        pow4(Mst1)*pow4(s2t))/(45.*pow2(Msq)*pow4(Mt)) - (2740*Dmglst1*pow2(
        Mst2)*pow4(Mst1)*pow4(s2t))/(9.*Mgl*pow2(Msq)*pow4(Mt)) - (100*z2*pow2(
        Mst2)*pow4(Mst1)*pow4(s2t))/(3.*pow2(Msq)*pow4(Mt)) - (400*Dmglst1*z2*
        pow2(Mst2)*pow4(Mst1)*pow4(s2t))/(3.*Mgl*pow2(Msq)*pow4(Mt)) - (6590*
        pow2(Dmglst1)*pow2(Mst2)*pow4(Mst1)*pow4(s2t))/(9.*pow2(Mgl)*pow2(Msq)*
        pow4(Mt)) - (1000*z2*pow2(Dmglst1)*pow2(Mst2)*pow4(Mst1)*pow4(s2t))/(3.
        *pow2(Mgl)*pow2(Msq)*pow4(Mt)) + (371*pow2(z2)*pow4(Mst1)*pow4(s2t))/(
        27.*pow4(Mt)) + (16*Dmglst1*pow2(z2)*pow4(Mst1)*pow4(s2t))/(3.*Mgl*
        pow4(Mt)) + (8*pow2(Dmglst1)*pow2(z2)*pow4(Mst1)*pow4(s2t))/(pow2(Mgl)*
        pow4(Mt)) + (53749*pow4(Mst2)*pow4(s2t))/(216.*pow4(Mt)) + (12*B4*pow4(
        Mst2)*pow4(s2t))/pow4(Mt) - (4*DN*pow4(Mst2)*pow4(s2t))/(3.*pow4(Mt)) +
        (7*Dmglst1*pow4(Mst2)*pow4(s2t))/(6.*Mgl*pow4(Mt)) + (8*OepS2*pow4(
        Mst2)*pow4(s2t))/(9.*pow4(Mt)) + (99*S2*pow4(Mst2)*pow4(s2t))/pow4(Mt)
        - (4*T1ep*pow4(Mst2)*pow4(s2t))/(3.*pow4(Mt)) - (109*z2*pow4(Mst2)*
        pow4(s2t))/pow4(Mt) + (66*Dmglst1*z2*pow4(Mst2)*pow4(s2t))/(Mgl*pow4(
        Mt)) + (328*z3*pow4(Mst2)*pow4(s2t))/(3.*pow4(Mt)) - (65*Dmglst1*z3*
        pow4(Mst2)*pow4(s2t))/(Mgl*pow4(Mt)) + (242*z4*pow4(Mst2)*pow4(s2t))/(
        3.*pow4(Mt)) + (6457*pow2(Dmglst1)*pow4(Mst2)*pow4(s2t))/(36.*pow2(Mgl)
        *pow4(Mt)) + (99*z2*pow2(Dmglst1)*pow4(Mst2)*pow4(s2t))/(pow2(Mgl)*
        pow4(Mt)) - (511*z3*pow2(Dmglst1)*pow4(Mst2)*pow4(s2t))/(2.*pow2(Mgl)*
        pow4(Mt)) + (10*pow2(Msq)*pow4(Mst2)*pow4(s2t))/(pow2(Mst1)*pow4(Mt)) -
        (20*z2*pow2(Msq)*pow4(Mst2)*pow4(s2t))/(pow2(Mst1)*pow4(Mt)) + (2714*
        pow2(Mst1)*pow4(Mst2)*pow4(s2t))/(45.*pow2(Msq)*pow4(Mt)) + (820*
        Dmglst1*pow2(Mst1)*pow4(Mst2)*pow4(s2t))/(9.*Mgl*pow2(Msq)*pow4(Mt)) +
        (10*z2*pow2(Mst1)*pow4(Mst2)*pow4(s2t))/(pow2(Msq)*pow4(Mt)) + (20*
        Dmglst1*z2*pow2(Mst1)*pow4(Mst2)*pow4(s2t))/(Mgl*pow2(Msq)*pow4(Mt)) +
        (180*pow2(Dmglst1)*pow2(Mst1)*pow4(Mst2)*pow4(s2t))/(pow2(Mgl)*pow2(
        Msq)*pow4(Mt)) + (30*z2*pow2(Dmglst1)*pow2(Mst1)*pow4(Mst2)*pow4(s2t))/
        (pow2(Mgl)*pow2(Msq)*pow4(Mt)) - (35*pow2(z2)*pow4(Mst2)*pow4(s2t))/(3.
        *pow4(Mt)) + (1009961*pow4(Mst1)*pow4(Mst2)*pow4(s2t))/(51450.*pow4(
        Msq)*pow4(Mt)) + (2605*Dmglst1*pow4(Mst1)*pow4(Mst2)*pow4(s2t))/(36.*
        Mgl*pow4(Msq)*pow4(Mt)) + (25*z2*pow4(Mst1)*pow4(Mst2)*pow4(s2t))/(6.*
        pow4(Msq)*pow4(Mt)) + (50*Dmglst1*z2*pow4(Mst1)*pow4(Mst2)*pow4(s2t))/(
        3.*Mgl*pow4(Msq)*pow4(Mt)) + (1775*pow2(Dmglst1)*pow4(Mst1)*pow4(Mst2)*
        pow4(s2t))/(8.*pow2(Mgl)*pow4(Msq)*pow4(Mt)) + (125*z2*pow2(Dmglst1)*
        pow4(Mst1)*pow4(Mst2)*pow4(s2t))/(3.*pow2(Mgl)*pow4(Msq)*pow4(Mt)) - (
        320*s2t*pow5(Mst1))/(3.*Mt*pow2(Msq)*pow2(Mst2)) - (3520*Dmglst1*s2t*
        pow5(Mst1))/(3.*Mgl*Mt*pow2(Msq)*pow2(Mst2)) - (1280*s2t*z2*pow5(Mst1))
        /(Mt*pow2(Msq)*pow2(Mst2)) - (6400*Dmglst1*s2t*z2*pow5(Mst1))/(Mgl*Mt*
        pow2(Msq)*pow2(Mst2)) - (5120*s2t*pow2(Dmglst1)*pow5(Mst1))/(Mt*pow2(
        Mgl)*pow2(Msq)*pow2(Mst2)) - (19200*s2t*z2*pow2(Dmglst1)*pow5(Mst1))/(
        Mt*pow2(Mgl)*pow2(Msq)*pow2(Mst2)) - (320*Cbeta*MuSUSY*pow2(s2t)*pow5(
        Mst1))/(Sbeta*pow2(Msq)*pow2(Mst2)*pow2(Mt)) - (1920*Cbeta*Dmglst1*
        MuSUSY*pow2(s2t)*pow5(Mst1))/(Mgl*Sbeta*pow2(Msq)*pow2(Mst2)*pow2(Mt))
        - (320*Cbeta*MuSUSY*z2*pow2(s2t)*pow5(Mst1))/(Sbeta*pow2(Msq)*pow2(
        Mst2)*pow2(Mt)) - (2880*Cbeta*Dmglst1*MuSUSY*z2*pow2(s2t)*pow5(Mst1))/(
        Mgl*Sbeta*pow2(Msq)*pow2(Mst2)*pow2(Mt)) - (9440*Cbeta*MuSUSY*pow2(
        Dmglst1)*pow2(s2t)*pow5(Mst1))/(Sbeta*pow2(Mgl)*pow2(Msq)*pow2(Mst2)*
        pow2(Mt)) - (10560*Cbeta*MuSUSY*z2*pow2(Dmglst1)*pow2(s2t)*pow5(Mst1))/
        (Sbeta*pow2(Mgl)*pow2(Msq)*pow2(Mst2)*pow2(Mt)) + (6640*pow3(s2t)*pow5(
        Mst1))/(9.*pow2(Msq)*pow3(Mt)) + (8560*Dmglst1*pow3(s2t)*pow5(Mst1))/(
        3.*Mgl*pow2(Msq)*pow3(Mt)) + (640*z2*pow3(s2t)*pow5(Mst1))/(3.*pow2(
        Msq)*pow3(Mt)) + (1280*Dmglst1*z2*pow3(s2t)*pow5(Mst1))/(Mgl*pow2(Msq)*
        pow3(Mt)) + (73120*pow2(Dmglst1)*pow3(s2t)*pow5(Mst1))/(9.*pow2(Mgl)*
        pow2(Msq)*pow3(Mt)) + (4160*z2*pow2(Dmglst1)*pow3(s2t)*pow5(Mst1))/(
        pow2(Mgl)*pow2(Msq)*pow3(Mt)) - (622151*pow3(s2t)*pow5(Mst1))/(324.*
        pow2(Mst2)*pow3(Mt)) - (3058187*Dmglst1*pow3(s2t)*pow5(Mst1))/(324.*
        Mgl*pow2(Mst2)*pow3(Mt)) + (2948*z2*pow3(s2t)*pow5(Mst1))/(9.*pow2(
        Mst2)*pow3(Mt)) + (5588*Dmglst1*z2*pow3(s2t)*pow5(Mst1))/(9.*Mgl*pow2(
        Mst2)*pow3(Mt)) + (1384*z3*pow3(s2t)*pow5(Mst1))/(pow2(Mst2)*pow3(Mt))
        + (7976*Dmglst1*z3*pow3(s2t)*pow5(Mst1))/(Mgl*pow2(Mst2)*pow3(Mt)) - (
        9306949*pow2(Dmglst1)*pow3(s2t)*pow5(Mst1))/(324.*pow2(Mgl)*pow2(Mst2)*
        pow3(Mt)) + (1780*z2*pow2(Dmglst1)*pow3(s2t)*pow5(Mst1))/(3.*pow2(Mgl)*
        pow2(Mst2)*pow3(Mt)) + (76472*z3*pow2(Dmglst1)*pow3(s2t)*pow5(Mst1))/(
        3.*pow2(Mgl)*pow2(Mst2)*pow3(Mt)) + (8000*s2t*pow5(Mst1))/(9.*Mt*pow4(
        Msq)) + (29440*Dmglst1*s2t*pow5(Mst1))/(9.*Mgl*Mt*pow4(Msq)) - (320*
        s2t*z2*pow5(Mst1))/(Mt*pow4(Msq)) - (2240*Dmglst1*s2t*z2*pow5(Mst1))/(
        3.*Mgl*Mt*pow4(Msq)) + (29600*s2t*pow2(Dmglst1)*pow5(Mst1))/(3.*Mt*
        pow2(Mgl)*pow4(Msq)) - (4160*s2t*z2*pow2(Dmglst1)*pow5(Mst1))/(3.*Mt*
        pow2(Mgl)*pow4(Msq)) - (400*Cbeta*MuSUSY*pow5(Mst1))/(Sbeta*pow2(Mst2)*
        pow4(Msq)) - (69680*Cbeta*Dmglst1*MuSUSY*pow5(Mst1))/(27.*Mgl*Sbeta*
        pow2(Mst2)*pow4(Msq)) - (257600*Cbeta*MuSUSY*pow2(Dmglst1)*pow5(Mst1))/
        (27.*Sbeta*pow2(Mgl)*pow2(Mst2)*pow4(Msq)) + (5360*s2t*pow2(MuSUSY)*
        pow5(Mst1))/(9.*Mt*pow2(Mst2)*pow4(Msq)) + (32720*Dmglst1*s2t*pow2(
        MuSUSY)*pow5(Mst1))/(9.*Mgl*Mt*pow2(Mst2)*pow4(Msq)) + (160*s2t*z2*
        pow2(MuSUSY)*pow5(Mst1))/(Mt*pow2(Mst2)*pow4(Msq)) + (800*Dmglst1*s2t*
        z2*pow2(MuSUSY)*pow5(Mst1))/(Mgl*Mt*pow2(Mst2)*pow4(Msq)) + (37120*s2t*
        pow2(Dmglst1)*pow2(MuSUSY)*pow5(Mst1))/(3.*Mt*pow2(Mgl)*pow2(Mst2)*
        pow4(Msq)) + (2400*s2t*z2*pow2(Dmglst1)*pow2(MuSUSY)*pow5(Mst1))/(Mt*
        pow2(Mgl)*pow2(Mst2)*pow4(Msq)) + (4555*Cbeta*MuSUSY*pow2(s2t)*pow5(
        Mst1))/(9.*Sbeta*pow2(Mt)*pow4(Msq)) + (25795*Cbeta*Dmglst1*MuSUSY*
        pow2(s2t)*pow5(Mst1))/(9.*Mgl*Sbeta*pow2(Mt)*pow4(Msq)) + (120*Cbeta*
        MuSUSY*z2*pow2(s2t)*pow5(Mst1))/(Sbeta*pow2(Mt)*pow4(Msq)) + (600*
        Cbeta*Dmglst1*MuSUSY*z2*pow2(s2t)*pow5(Mst1))/(Mgl*Sbeta*pow2(Mt)*pow4(
        Msq)) + (85855*Cbeta*MuSUSY*pow2(Dmglst1)*pow2(s2t)*pow5(Mst1))/(9.*
        Sbeta*pow2(Mgl)*pow2(Mt)*pow4(Msq)) + (1800*Cbeta*MuSUSY*z2*pow2(
        Dmglst1)*pow2(s2t)*pow5(Mst1))/(Sbeta*pow2(Mgl)*pow2(Mt)*pow4(Msq)) - (
        5360*s2t*pow2(MuSUSY)*pow5(Mst1))/(9.*Mt*pow2(Mst2)*pow2(Sbeta)*pow4(
        Msq)) - (32720*Dmglst1*s2t*pow2(MuSUSY)*pow5(Mst1))/(9.*Mgl*Mt*pow2(
        Mst2)*pow2(Sbeta)*pow4(Msq)) - (160*s2t*z2*pow2(MuSUSY)*pow5(Mst1))/(
        Mt*pow2(Mst2)*pow2(Sbeta)*pow4(Msq)) - (800*Dmglst1*s2t*z2*pow2(MuSUSY)
        *pow5(Mst1))/(Mgl*Mt*pow2(Mst2)*pow2(Sbeta)*pow4(Msq)) - (37120*s2t*
        pow2(Dmglst1)*pow2(MuSUSY)*pow5(Mst1))/(3.*Mt*pow2(Mgl)*pow2(Mst2)*
        pow2(Sbeta)*pow4(Msq)) - (2400*s2t*z2*pow2(Dmglst1)*pow2(MuSUSY)*pow5(
        Mst1))/(Mt*pow2(Mgl)*pow2(Mst2)*pow2(Sbeta)*pow4(Msq)) - (1715*pow2(
        Mst2)*pow3(s2t)*pow5(Mst1))/(9.*pow3(Mt)*pow4(Msq)) - (9035*Dmglst1*
        pow2(Mst2)*pow3(s2t)*pow5(Mst1))/(9.*Mgl*pow3(Mt)*pow4(Msq)) - (40*z2*
        pow2(Mst2)*pow3(s2t)*pow5(Mst1))/(pow3(Mt)*pow4(Msq)) - (200*Dmglst1*
        z2*pow2(Mst2)*pow3(s2t)*pow5(Mst1))/(Mgl*pow3(Mt)*pow4(Msq)) - (9805*
        pow2(Dmglst1)*pow2(Mst2)*pow3(s2t)*pow5(Mst1))/(3.*pow2(Mgl)*pow3(Mt)*
        pow4(Msq)) - (600*z2*pow2(Dmglst1)*pow2(Mst2)*pow3(s2t)*pow5(Mst1))/(
        pow2(Mgl)*pow3(Mt)*pow4(Msq)) - (324688*s2t*pow5(Mst1))/(81.*Mt*pow4(
        Mst2)) + (864592*Dmglst1*s2t*pow5(Mst1))/(81.*Mgl*Mt*pow4(Mst2)) + (
        146848*s2t*z2*pow5(Mst1))/(9.*Mt*pow4(Mst2)) + (235936*Dmglst1*s2t*z2*
        pow5(Mst1))/(3.*Mgl*Mt*pow4(Mst2)) - (28672*s2t*z3*pow5(Mst1))/(3.*Mt*
        pow4(Mst2)) - (46848*Dmglst1*s2t*z3*pow5(Mst1))/(Mgl*Mt*pow4(Mst2)) + (
        716416*s2t*pow2(Dmglst1)*pow5(Mst1))/(9.*Mt*pow2(Mgl)*pow4(Mst2)) + (
        2075968*s2t*z2*pow2(Dmglst1)*pow5(Mst1))/(9.*Mt*pow2(Mgl)*pow4(Mst2)) -
        (410752*s2t*z3*pow2(Dmglst1)*pow5(Mst1))/(3.*Mt*pow2(Mgl)*pow4(Mst2)) -
        (3520*Cbeta*MuSUSY*pow5(Mst1))/(3.*Sbeta*pow2(Msq)*pow4(Mst2)) - (
        140480*Cbeta*Dmglst1*MuSUSY*pow5(Mst1))/(27.*Mgl*Sbeta*pow2(Msq)*pow4(
        Mst2)) + (1280*Cbeta*MuSUSY*z2*pow5(Mst1))/(Sbeta*pow2(Msq)*pow4(Mst2))
        + (6400*Cbeta*Dmglst1*MuSUSY*z2*pow5(Mst1))/(Mgl*Sbeta*pow2(Msq)*pow4(
        Mst2)) - (375680*Cbeta*MuSUSY*pow2(Dmglst1)*pow5(Mst1))/(27.*Sbeta*
        pow2(Mgl)*pow2(Msq)*pow4(Mst2)) + (19200*Cbeta*MuSUSY*z2*pow2(Dmglst1)*
        pow5(Mst1))/(Sbeta*pow2(Mgl)*pow2(Msq)*pow4(Mst2)) + (17920*s2t*pow2(
        MuSUSY)*pow5(Mst1))/(9.*Mt*pow2(Msq)*pow4(Mst2)) + (18560*Dmglst1*s2t*
        pow2(MuSUSY)*pow5(Mst1))/(3.*Mgl*Mt*pow2(Msq)*pow4(Mst2)) - (2560*
        Dmglst1*s2t*z2*pow2(MuSUSY)*pow5(Mst1))/(Mgl*Mt*pow2(Msq)*pow4(Mst2)) +
        (64960*s2t*pow2(Dmglst1)*pow2(MuSUSY)*pow5(Mst1))/(9.*Mt*pow2(Mgl)*
        pow2(Msq)*pow4(Mst2)) - (11520*s2t*z2*pow2(Dmglst1)*pow2(MuSUSY)*pow5(
        Mst1))/(Mt*pow2(Mgl)*pow2(Msq)*pow4(Mst2)) + (80287*Cbeta*MuSUSY*pow2(
        s2t)*pow5(Mst1))/(12.*Sbeta*pow2(Mt)*pow4(Mst2)) + (1152*B4*Cbeta*
        MuSUSY*pow2(s2t)*pow5(Mst1))/(Sbeta*pow2(Mt)*pow4(Mst2)) + (3539195*
        Cbeta*Dmglst1*MuSUSY*pow2(s2t)*pow5(Mst1))/(108.*Mgl*Sbeta*pow2(Mt)*
        pow4(Mst2)) + (2368*B4*Cbeta*Dmglst1*MuSUSY*pow2(s2t)*pow5(Mst1))/(Mgl*
        Sbeta*pow2(Mt)*pow4(Mst2)) - (32*Cbeta*Dmglst1*DN*MuSUSY*pow2(s2t)*
        pow5(Mst1))/(Mgl*Sbeta*pow2(Mt)*pow4(Mst2)) - (7156*Cbeta*MuSUSY*z2*
        pow2(s2t)*pow5(Mst1))/(3.*Sbeta*pow2(Mt)*pow4(Mst2)) + (25948*Cbeta*
        Dmglst1*MuSUSY*z2*pow2(s2t)*pow5(Mst1))/(3.*Mgl*Sbeta*pow2(Mt)*pow4(
        Mst2)) + (2920*Cbeta*MuSUSY*z3*pow2(s2t)*pow5(Mst1))/(Sbeta*pow2(Mt)*
        pow4(Mst2)) - (17464*Cbeta*Dmglst1*MuSUSY*z3*pow2(s2t)*pow5(Mst1))/(
        Mgl*Sbeta*pow2(Mt)*pow4(Mst2)) + (576*Cbeta*MuSUSY*z4*pow2(s2t)*pow5(
        Mst1))/(Sbeta*pow2(Mt)*pow4(Mst2)) - (5984*Cbeta*Dmglst1*MuSUSY*z4*
        pow2(s2t)*pow5(Mst1))/(Mgl*Sbeta*pow2(Mt)*pow4(Mst2)) + (10707109*
        Cbeta*MuSUSY*pow2(Dmglst1)*pow2(s2t)*pow5(Mst1))/(108.*Sbeta*pow2(Mgl)*
        pow2(Mt)*pow4(Mst2)) + (4192*B4*Cbeta*MuSUSY*pow2(Dmglst1)*pow2(s2t)*
        pow5(Mst1))/(Sbeta*pow2(Mgl)*pow2(Mt)*pow4(Mst2)) - (80*Cbeta*DN*
        MuSUSY*pow2(Dmglst1)*pow2(s2t)*pow5(Mst1))/(Sbeta*pow2(Mgl)*pow2(Mt)*
        pow4(Mst2)) + (27724*Cbeta*MuSUSY*z2*pow2(Dmglst1)*pow2(s2t)*pow5(Mst1)
        )/(Sbeta*pow2(Mgl)*pow2(Mt)*pow4(Mst2)) - (75624*Cbeta*MuSUSY*z3*pow2(
        Dmglst1)*pow2(s2t)*pow5(Mst1))/(Sbeta*pow2(Mgl)*pow2(Mt)*pow4(Mst2)) -
        (15824*Cbeta*MuSUSY*z4*pow2(Dmglst1)*pow2(s2t)*pow5(Mst1))/(Sbeta*pow2(
        Mgl)*pow2(Mt)*pow4(Mst2)) - (17920*s2t*pow2(MuSUSY)*pow5(Mst1))/(9.*Mt*
        pow2(Msq)*pow2(Sbeta)*pow4(Mst2)) - (18560*Dmglst1*s2t*pow2(MuSUSY)*
        pow5(Mst1))/(3.*Mgl*Mt*pow2(Msq)*pow2(Sbeta)*pow4(Mst2)) + (2560*
        Dmglst1*s2t*z2*pow2(MuSUSY)*pow5(Mst1))/(Mgl*Mt*pow2(Msq)*pow2(Sbeta)*
        pow4(Mst2)) - (64960*s2t*pow2(Dmglst1)*pow2(MuSUSY)*pow5(Mst1))/(9.*Mt*
        pow2(Mgl)*pow2(Msq)*pow2(Sbeta)*pow4(Mst2)) + (11520*s2t*z2*pow2(
        Dmglst1)*pow2(MuSUSY)*pow5(Mst1))/(Mt*pow2(Mgl)*pow2(Msq)*pow2(Sbeta)*
        pow4(Mst2)) - (1152*Cbeta*MuSUSY*pow2(s2t)*pow2(z2)*pow5(Mst1))/(Sbeta*
        pow2(Mt)*pow4(Mst2)) - (2432*Cbeta*Dmglst1*MuSUSY*pow2(s2t)*pow2(z2)*
        pow5(Mst1))/(Mgl*Sbeta*pow2(Mt)*pow4(Mst2)) - (4352*Cbeta*MuSUSY*pow2(
        Dmglst1)*pow2(s2t)*pow2(z2)*pow5(Mst1))/(Sbeta*pow2(Mgl)*pow2(Mt)*pow4(
        Mst2)) + (pow3(log(pow2(Mst1)/pow2(Mst2)))*(4*Dmglst1*Mgl*Mst1*(4*Mt*
        s2t*Sbeta*(2667*Cbeta*Mt*MuSUSY*s2t + 7320*Sbeta*pow2(Mt) - 323*Sbeta*
        pow2(Mst2)*pow2(s2t))*pow4(Mst1) + 8*Mt*s2t*(-28*pow2(Mt) + 73*pow2(
        Mst2)*pow2(s2t))*pow2(Sbeta)*pow4(Mst2) - 8*Mt*pow2(Mst1)*(-849*Cbeta*
        Mt*MuSUSY*Sbeta*pow2(Mst2)*pow2(s2t) - 10*s2t*pow2(Mt)*(84*pow2(MuSUSY)
        *(-1 + pow2(Sbeta)) + 169*pow2(Mst2)*pow2(Sbeta)) + 1662*Cbeta*MuSUSY*
        Sbeta*pow3(Mt) + 356*pow2(Sbeta)*pow3(s2t)*pow4(Mst2)) + Mst1*pow2(
        Mst2)*pow2(Sbeta)*(-320*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 1056*pow4(Mt) +
        449*pow4(Mst2)*pow4(s2t)) + pow3(Mst1)*(-8*pow2(Mt)*pow2(s2t)*(1209*
        pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 908*pow2(Mst2)*pow2(Sbeta)) + 15168*
        Cbeta*MuSUSY*s2t*Sbeta*pow3(Mt) - 3040*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)
        *pow3(s2t) - 5472*pow2(Sbeta)*pow4(Mt) + 311*pow2(Sbeta)*pow4(Mst2)*
        pow4(s2t))) + 2*Mst1*pow2(Dmglst1)*(8*Mt*s2t*Sbeta*(6618*Cbeta*Mt*
        MuSUSY*s2t + 20892*Sbeta*pow2(Mt) - 689*Sbeta*pow2(Mst2)*pow2(s2t))*
        pow4(Mst1) + 16*Mt*s2t*(-28*pow2(Mt) + 73*pow2(Mst2)*pow2(s2t))*pow2(
        Sbeta)*pow4(Mst2) - 8*Mt*pow2(Mst1)*(-1698*Cbeta*Mt*MuSUSY*Sbeta*pow2(
        Mst2)*pow2(s2t) - 56*s2t*pow2(Mt)*(30*pow2(MuSUSY)*(-1 + pow2(Sbeta)) +
        121*pow2(Mst2)*pow2(Sbeta)) + 3324*Cbeta*MuSUSY*Sbeta*pow3(Mt) + 1663*
        pow2(Sbeta)*pow3(s2t)*pow4(Mst2)) + 3*Mst1*pow2(Mst2)*pow2(Sbeta)*(-
        320*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 1056*pow4(Mt) + 449*pow4(Mst2)*
        pow4(s2t)) + pow3(Mst1)*(-8*pow2(Mt)*pow2(s2t)*(4291*pow2(MuSUSY)*(-1 +
        pow2(Sbeta)) + 4936*pow2(Mst2)*pow2(Sbeta)) + 80896*Cbeta*MuSUSY*s2t*
        Sbeta*pow3(Mt) - 11776*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*pow3(s2t) -
        28736*pow2(Sbeta)*pow4(Mt) + 1597*pow2(Sbeta)*pow4(Mst2)*pow4(s2t))) +
        pow2(Mgl)*(32*Mst1*Mt*s2t*(-28*pow2(Mt) + 73*pow2(Mst2)*pow2(s2t))*
        pow2(Sbeta)*pow4(Mst2) - 32*Mt*pow3(Mst1)*(102*Cbeta*Mt*MuSUSY*Sbeta*
        pow2(Mst2)*pow2(s2t) + 2*s2t*pow2(Mt)*(214*pow2(MuSUSY)*(-1 + pow2(
        Sbeta)) - 279*pow2(Mst2)*pow2(Sbeta)) + 530*Cbeta*MuSUSY*Sbeta*pow3(Mt)
        + 39*pow2(Sbeta)*pow3(s2t)*pow4(Mst2)) + pow2(Sbeta)*pow4(Mst2)*(-3152*
        pow2(Mst2)*pow2(Mt)*pow2(s2t) - 544*pow4(Mt) + 639*pow4(Mst2)*pow4(s2t)
        ) + pow4(Mst1)*(-32*pow2(Mt)*pow2(s2t)*(1960*pow2(MuSUSY)*(-1 + pow2(
        Sbeta)) + 223*pow2(Mst2)*pow2(Sbeta)) + 24416*Cbeta*MuSUSY*s2t*Sbeta*
        pow3(Mt) - 17472*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*pow3(s2t) - 4640*
        pow2(Sbeta)*pow4(Mt) + 1535*pow2(Sbeta)*pow4(Mst2)*pow4(s2t)) + 16*Mt*
        s2t*Sbeta*(-207*Cbeta*Mt*MuSUSY*s2t + 1536*Sbeta*pow2(Mt) + Sbeta*pow2(
        Mst2)*pow2(s2t))*pow5(Mst1) - 2*pow2(Mst1)*(64*pow2(Mst2)*pow2(Mt)*
        pow2(s2t)*(217*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 15*pow2(Mst2)*pow2(
        Sbeta)) - 5072*Cbeta*MuSUSY*s2t*Sbeta*pow2(Mst2)*pow3(Mt) + 5666*Cbeta*
        Mt*MuSUSY*Sbeta*pow3(s2t)*pow4(Mst2) + 32*(4*pow2(MuSUSY)*(-1 + pow2(
        Sbeta)) - 35*pow2(Mst2)*pow2(Sbeta))*pow4(Mt) - 1097*pow2(Sbeta)*pow4(
        s2t)*pow6(Mst2)))))/(18.*pow2(Mgl)*pow2(Sbeta)*pow4(Mst2)*pow4(Mt)) - (
        pow2(log(pow2(Mst1)/pow2(Msq)))*(pow2(Mgl)*(-224*Sbeta*pow2(Msq)*pow2(
        Mt)*pow4(Mst2)*(2*Cbeta*Mt*MuSUSY*s2t*pow2(Mst1) + Sbeta*(20*Mst1*Mt*
        s2t*pow2(Mst2) + 7*pow2(Mst2)*pow2(Mt) + pow2(Mst1)*(37*pow2(Mt) - 2*
        pow2(Mst2)*pow2(s2t)) - 20*Mt*s2t*pow3(Mst1) + pow2(s2t)*pow4(Mst1) +
        pow2(s2t)*pow4(Mst2))) - 2*pow4(Mst2)*(-560*Sbeta*(3*Cbeta*Mt*MuSUSY -
        2*s2t*Sbeta*pow2(Mst2))*pow3(Mst1)*pow3(Mt) - 2*Mt*pow2(Mst1)*pow2(
        Mst2)*(Mt*pow2(s2t)*(4*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 15*pow2(Mst2)*
        pow2(Sbeta)) - 420*pow2(Sbeta)*pow3(Mt) + 3*Cbeta*MuSUSY*Sbeta*pow2(
        Mst2)*pow3(s2t)) + 4*pow2(Mt)*(-15*Cbeta*Mt*MuSUSY*s2t*Sbeta + pow2(
        MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta)) + 57*pow2(Mt)*pow2(Sbeta))*pow4(
        Mst2) + 560*Mst1*s2t*pow2(Sbeta)*pow3(Mt)*pow4(Mst2) + pow4(Mst1)*(2*
        pow2(Mt)*pow2(s2t)*(2*pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 15*pow2(Mst2)*
        pow2(Sbeta)) + 60*Cbeta*MuSUSY*s2t*Sbeta*pow3(Mt) + 6*Cbeta*Mt*MuSUSY*
        Sbeta*pow2(Mst2)*pow3(s2t) + 2468*pow2(Sbeta)*pow4(Mt) - 3*pow2(Sbeta)*
        pow4(Mst2)*pow4(s2t)) - 1680*s2t*pow2(Sbeta)*pow3(Mt)*pow5(Mst1)) + 35*
        pow4(Msq)*(3*Sbeta*pow2(Mst2)*pow2(s2t)*(-12*Cbeta*Mt*MuSUSY*s2t + 24*
        Sbeta*pow2(Mt) + Sbeta*pow2(Mst2)*pow2(s2t))*pow4(Mst1) + 96*Mt*s2t*
        pow3(Mst1)*(4*pow2(Mt)*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + pow2(s2t)*
        pow2(Sbeta)*pow4(Mst2)) + pow2(Sbeta)*pow4(Mst2)*(72*pow2(Mst2)*pow2(
        Mt)*pow2(s2t) - 176*pow4(Mt) - 15*pow4(Mst2)*pow4(s2t)) + 3*pow2(Mst1)*
        pow2(Mst2)*(24*pow2(Mt)*pow2(s2t)*(pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 2*
        pow2(Mst2)*pow2(Sbeta)) + 48*Cbeta*MuSUSY*s2t*Sbeta*pow3(Mt) - 8*Cbeta*
        Mt*MuSUSY*Sbeta*pow2(Mst2)*pow3(s2t) - 48*pow2(Sbeta)*pow4(Mt) + 7*
        pow2(Sbeta)*pow4(Mst2)*pow4(s2t)) - 96*Mst1*Mt*pow2(Sbeta)*pow3(s2t)*
        pow6(Mst2))) + 35*pow2(Dmglst1)*(64*Mst1*pow2(Msq)*(-9*Mst1*Mt + 2*s2t*
        pow2(Mst1) - 2*s2t*pow2(Mst2))*pow2(Sbeta)*pow3(Mt)*pow4(Mst2) + 16*
        Mst1*Sbeta*pow3(Mt)*pow4(Mst2)*(14*Cbeta*Mt*MuSUSY*pow2(Mst1) + Sbeta*(
        -9*Mst1*Mt*pow2(Mst2) - 24*s2t*pow2(Mst1)*pow2(Mst2) - 59*Mt*pow3(Mst1)
        + 26*s2t*pow4(Mst1) - 2*s2t*pow4(Mst2))) + 3*pow4(Msq)*(9*Sbeta*pow2(
        Mst2)*pow2(s2t)*(-4*Cbeta*Mt*MuSUSY*s2t + 8*Sbeta*pow2(Mt) + Sbeta*
        pow2(Mst2)*pow2(s2t))*pow4(Mst1) - 9*pow2(Sbeta)*pow2(-4*pow2(Mt) +
        pow2(Mst2)*pow2(s2t))*pow4(Mst2) + 32*Mt*s2t*pow3(Mst1)*(4*pow2(Mt)*
        pow2(MuSUSY)*(-1 + pow2(Sbeta)) + pow2(s2t)*pow2(Sbeta)*pow4(Mst2)) +
        9*pow2(Mst1)*pow2(Mst2)*pow2(Sbeta)*(-16*pow2(Mst2)*pow2(Mt)*pow2(s2t)
        - 16*pow4(Mt) + pow4(Mst2)*pow4(s2t)) - 32*Mst1*Mt*pow2(Sbeta)*pow3(
        s2t)*pow6(Mst2))) + 70*Dmglst1*Mgl*(64*Mst1*pow2(Msq)*(-3*Mst1*Mt +
        s2t*pow2(Mst1) - s2t*pow2(Mst2))*pow2(Sbeta)*pow3(Mt)*pow4(Mst2) + 16*
        Mst1*Sbeta*pow3(Mt)*pow4(Mst2)*(7*Cbeta*Mt*MuSUSY*pow2(Mst1) + Sbeta*(-
        3*Mst1*Mt*pow2(Mst2) - 6*s2t*pow2(Mst1)*pow2(Mst2) - 13*Mt*pow3(Mst1) +
        7*s2t*pow4(Mst1) - s2t*pow4(Mst2))) + 3*pow4(Msq)*(3*Sbeta*pow2(Mst2)*
        pow2(s2t)*(-4*Cbeta*Mt*MuSUSY*s2t + 8*Sbeta*pow2(Mt) + Sbeta*pow2(Mst2)
        *pow2(s2t))*pow4(Mst1) - 3*pow2(Sbeta)*pow2(-4*pow2(Mt) + pow2(Mst2)*
        pow2(s2t))*pow4(Mst2) + 16*Mt*s2t*pow3(Mst1)*(4*pow2(Mt)*pow2(MuSUSY)*(
        -1 + pow2(Sbeta)) + pow2(s2t)*pow2(Sbeta)*pow4(Mst2)) + 3*pow2(Mst1)*
        pow2(Mst2)*pow2(Sbeta)*(-16*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 16*pow4(Mt)
        + pow4(Mst2)*pow4(s2t)) - 16*Mst1*Mt*pow2(Sbeta)*pow3(s2t)*pow6(Mst2)))
        + 105*log(pow2(Mst1)/pow2(Mst2))*pow4(Msq)*(pow2(Mgl)*(-16*Mst1*Mt*s2t*
        (-4*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow2(Sbeta)*pow4(Mst2) - 2*pow2(
        Mst1)*pow2(Mst2)*pow2(s2t)*(-14*Cbeta*Mt*MuSUSY*s2t*Sbeta*pow2(Mst2) +
        4*pow2(Mt)*(-8*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + pow2(Mst2)*pow2(Sbeta)
        ) + 3*pow2(s2t)*pow2(Sbeta)*pow4(Mst2)) + pow2(s2t)*pow4(Mst1)*(64*
        pow2(Mt)*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 7*pow2(s2t)*pow2(Sbeta)*
        pow4(Mst2)) - 16*Mt*pow3(Mst1)*(-6*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*
        pow2(s2t) - 12*s2t*pow2(Mt)*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 4*Cbeta*
        MuSUSY*Sbeta*pow3(Mt) + pow2(Sbeta)*pow3(s2t)*pow4(Mst2)) + pow2(Sbeta)
        *pow4(Mst2)*(8*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 16*pow4(Mt) - pow4(Mst2)
        *pow4(s2t)) + 96*Cbeta*MuSUSY*Sbeta*pow2(Mt)*pow2(s2t)*pow5(Mst1)) + 2*
        Mst1*pow2(Dmglst1)*(48*Cbeta*MuSUSY*Sbeta*pow2(Mt)*pow2(s2t)*pow4(Mst1)
        - 8*Mt*s2t*(-4*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow2(Sbeta)*pow4(Mst2)
        + 9*pow2(s2t)*pow3(Mst1)*(8*pow2(Mt)*pow2(MuSUSY)*(-1 + pow2(Sbeta)) +
        pow2(s2t)*pow2(Sbeta)*pow4(Mst2)) - 8*Mt*pow2(Mst1)*(-6*Cbeta*Mt*
        MuSUSY*Sbeta*pow2(Mst2)*pow2(s2t) - 12*s2t*pow2(Mt)*pow2(MuSUSY)*(-1 +
        pow2(Sbeta)) + 4*Cbeta*MuSUSY*Sbeta*pow3(Mt) + pow2(Sbeta)*pow3(s2t)*
        pow4(Mst2)) - 9*Mst1*pow2(Sbeta)*pow4(s2t)*pow6(Mst2)) + 4*Dmglst1*Mgl*
        Mst1*(24*Cbeta*MuSUSY*Sbeta*pow2(Mt)*pow2(s2t)*pow4(Mst1) - 4*Mt*s2t*(-
        4*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow2(Sbeta)*pow4(Mst2) + 3*pow2(s2t)
        *pow3(Mst1)*(8*pow2(Mt)*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + pow2(s2t)*
        pow2(Sbeta)*pow4(Mst2)) - 4*Mt*pow2(Mst1)*(-6*Cbeta*Mt*MuSUSY*Sbeta*
        pow2(Mst2)*pow2(s2t) - 12*s2t*pow2(Mt)*pow2(MuSUSY)*(-1 + pow2(Sbeta))
        + 4*Cbeta*MuSUSY*Sbeta*pow3(Mt) + pow2(Sbeta)*pow3(s2t)*pow4(Mst2)) -
        3*Mst1*pow2(Sbeta)*pow4(s2t)*pow6(Mst2)))))/(21.*pow2(Mgl)*pow2(Sbeta)*
        pow4(Msq)*pow4(Mst2)*pow4(Mt)) + (log(pow2(Mst1)/pow2(Mst2))*(98*pow2(
        Dmglst1)*pow2(Mst1)*(pow4(Msq)*(800*Mst1*Mt*s2t*(-4*(1067 + 648*z2 -
        432*z3)*pow2(Mt) + 3*(7603 + 1152*z2 - 108*z3)*pow2(Mst2)*pow2(s2t))*
        pow2(Sbeta)*pow4(Mst2) - 800*Mt*pow3(Mst1)*(-9*Cbeta*Mt*MuSUSY*Sbeta*(-
        1195 + 6312*z2 - 8208*z3)*pow2(Mst2)*pow2(s2t) - 16*s2t*pow2(Mt)*(3*(-
        2047 + 1290*z2 - 2025*z3)*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + (13499 +
        14742*z2)*pow2(Mst2)*pow2(Sbeta)) + 4*Cbeta*MuSUSY*Sbeta*(15889 +
        29052*z2 + 432*z3)*pow3(Mt) + 6*(10423 + 8214*z2 - 8370*z3)*pow2(Sbeta)
        *pow3(s2t)*pow4(Mst2)) - 4*pow2(Sbeta)*pow4(Mst2)*(-266400*pow2(Mst2)*
        pow2(Mt)*pow2(s2t) + 50912*pow4(Mt) + 100125*pow4(Mst2)*pow4(s2t)) +
        400*pow2(Mst1)*pow2(Mst2)*pow2(Sbeta)*(6*(20783 + 2736*z2 + 432*z3)*
        pow2(Mst2)*pow2(Mt)*pow2(s2t) + 8*(6343 + 3888*z2)*pow4(Mt) + 9*(628 +
        945*z2 - 414*z3)*pow4(Mst2)*pow4(s2t)) - 25*pow4(Mst1)*(8*pow2(Mt)*
        pow2(s2t)*((347507 + 207504*z2 - 45792*z3)*pow2(MuSUSY)*(-1 + pow2(
        Sbeta)) + 4*(-76711 + 102456*z2)*pow2(Mst2)*pow2(Sbeta)) - 512*Cbeta*
        MuSUSY*s2t*Sbeta*(-20753 + 11781*z2 - 162*z3)*pow3(Mt) + 4*Cbeta*Mt*
        MuSUSY*Sbeta*(289115 + 71424*z2 + 13824*z3)*pow2(Mst2)*pow3(s2t) + 192*
        (26965 + 28752*z2)*pow2(Sbeta)*pow4(Mt) + (-214703 + 64656*z2 - 73440*
        z3)*pow2(Sbeta)*pow4(Mst2)*pow4(s2t)) + 200*Mt*s2t*Sbeta*(3*Cbeta*Mt*
        MuSUSY*s2t*(-95783 + 345888*z2 - 199584*z3) + 208*Sbeta*(20965 + 14292*
        z2)*pow2(Mt) + Sbeta*(254699 - 162576*z2)*pow2(Mst2)*pow2(s2t))*pow5(
        Mst1)) - 750*Mst1*pow2(Mst2)*(144*Mst1*pow2(Mt)*(23*pow2(Mt) + 3*pow2(
        Mst2)*pow2(s2t))*pow2(Sbeta)*pow4(Mst2) - 16*Mt*pow2(Mst1)*pow2(Mst2)*(
        -6*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*pow2(s2t) - 4*s2t*pow2(Mt)*(pow2(
        MuSUSY)*(-1 + pow2(Sbeta)) + 102*pow2(Mst2)*pow2(Sbeta)) + 178*Cbeta*
        MuSUSY*Sbeta*pow3(Mt) + 3*pow2(Sbeta)*pow3(s2t)*pow4(Mst2)) + 8*Mt*
        pow4(Mst1)*(-4551*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*pow2(s2t) - 60*s2t*
        pow2(Mt)*(101*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 21*pow2(Mst2)*pow2(
        Sbeta)) + 1720*Cbeta*MuSUSY*Sbeta*pow3(Mt) + 1521*pow2(Sbeta)*pow3(s2t)
        *pow4(Mst2)) + pow2(Mst2)*pow3(Mst1)*(8*pow2(Mt)*pow2(s2t)*(871*pow2(
        MuSUSY)*(-1 + pow2(Sbeta)) + 246*pow2(Mst2)*pow2(Sbeta)) - 4800*Cbeta*
        MuSUSY*s2t*Sbeta*pow3(Mt) + 3484*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*pow3(
        s2t) + 4400*pow2(Sbeta)*pow4(Mt) - 871*pow2(Sbeta)*pow4(Mst2)*pow4(s2t)
        ) + 16*s2t*Sbeta*(-3*Cbeta*MuSUSY*s2t + 58*Mt*Sbeta)*pow2(Mt)*pow6(
        Mst2)) + 3000*Mst1*pow2(Msq)*(32*Mt*pow4(Mst1)*(144*Cbeta*Mt*MuSUSY*
        Sbeta*pow2(Mst2)*pow2(s2t) + 6*s2t*pow2(Mt)*(122*pow2(MuSUSY)*(-1 +
        pow2(Sbeta)) - 141*pow2(Mst2)*pow2(Sbeta)) + 314*Cbeta*MuSUSY*Sbeta*
        pow3(Mt) + 87*pow2(Sbeta)*pow3(s2t)*pow4(Mst2)) - 32*Mt*pow2(Mst1)*
        pow2(Mst2)*(-198*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*pow2(s2t) - 24*s2t*
        pow2(Mt)*(11*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 8*pow2(Mst2)*pow2(Sbeta)
        ) + 130*Cbeta*MuSUSY*Sbeta*pow3(Mt) + 135*pow2(Sbeta)*pow3(s2t)*pow4(
        Mst2)) - 4*s2t*pow2(Mst2)*pow3(Mst1)*(-701*Cbeta*Mt*MuSUSY*Sbeta*pow2(
        Mst2)*pow2(s2t) + pow2(Mt)*(-988*s2t*pow2(MuSUSY)*(-1 + pow2(Sbeta)) +
        492*s2t*pow2(Mst2)*pow2(Sbeta)) - 1200*Cbeta*MuSUSY*Sbeta*pow3(Mt) +
        227*pow2(Sbeta)*pow3(s2t)*pow4(Mst2)) + 9*Mst1*pow2(Sbeta)*pow4(Mst2)*(
        -48*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 16*pow4(Mt) + 23*pow4(Mst2)*pow4(
        s2t)) - 256*s2t*pow2(Sbeta)*pow3(Mt)*pow6(Mst2))) + 4900*Dmglst1*Mgl*
        pow2(Mst1)*(pow4(Msq)*(16*Mst1*Mt*s2t*(-8*(145 + 324*z2 - 216*z3)*pow2(
        Mt) + 27*(777 + 128*z2 - 12*z3)*pow2(Mst2)*pow2(s2t))*pow2(Sbeta)*pow4(
        Mst2) - 16*Mt*pow3(Mst1)*(-9*Cbeta*Mt*MuSUSY*Sbeta*(-1195 + 6312*z2 -
        8208*z3)*pow2(Mst2)*pow2(s2t) - 12*s2t*pow2(Mt)*(4*(-2047 + 1290*z2 -
        2025*z3)*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 9*(947 + 1100*z2)*pow2(Mst2)
        *pow2(Sbeta)) + 4*Cbeta*MuSUSY*Sbeta*(15889 + 29052*z2 + 432*z3)*pow3(
        Mt) + 6*(2899 + 3732*z2 - 4158*z3)*pow2(Sbeta)*pow3(s2t)*pow4(Mst2)) -
        4*pow2(Sbeta)*pow4(Mst2)*(-2304*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 2032*
        pow4(Mt) + 1017*pow4(Mst2)*pow4(s2t)) + 16*pow2(Mst1)*pow2(Mst2)*pow2(
        Sbeta)*(2*(21991 + 2736*z2 + 432*z3)*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 8*
        (2069 + 1296*z2)*pow4(Mt) + 27*(57 + 105*z2 - 46*z3)*pow4(Mst2)*pow4(
        s2t)) - pow4(Mst1)*(8*pow2(Mt)*pow2(s2t)*((68551 + 69552*z2 - 19872*z3)
        *pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 4*(-37061 + 17640*z2)*pow2(Mst2)*
        pow2(Sbeta)) - 768*Cbeta*MuSUSY*s2t*Sbeta*(-5305 + 1242*z2 - 36*z3)*
        pow3(Mt) + 4*Cbeta*Mt*MuSUSY*Sbeta*(52063 + 24192*z2)*pow2(Mst2)*pow3(
        s2t) + 512*(2152 + 2079*z2)*pow2(Sbeta)*pow4(Mt) + (-31507 + 21168*z2 -
        19872*z3)*pow2(Sbeta)*pow4(Mst2)*pow4(s2t)) + 4*Mt*s2t*Sbeta*(27*Cbeta*
        Mt*MuSUSY*s2t*(-10489 + 15248*z2 - 10944*z3) + 256*Sbeta*(5696 + 3987*
        z2)*pow2(Mt) + 3*Sbeta*(26687 - 20496*z2)*pow2(Mst2)*pow2(s2t))*pow5(
        Mst1)) - 30*Mst1*pow2(Mst2)*(48*Mst1*pow2(Mt)*(23*pow2(Mt) + 3*pow2(
        Mst2)*pow2(s2t))*pow2(Sbeta)*pow4(Mst2) - 8*Mt*pow2(Mst1)*pow2(Mst2)*(-
        6*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*pow2(s2t) - 4*s2t*pow2(Mt)*(pow2(
        MuSUSY)*(-1 + pow2(Sbeta)) + 51*pow2(Mst2)*pow2(Sbeta)) + 178*Cbeta*
        MuSUSY*Sbeta*pow3(Mt) + 3*pow2(Sbeta)*pow3(s2t)*pow4(Mst2)) + 4*Mt*
        pow4(Mst1)*(-1437*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*pow2(s2t) - 12*s2t*
        pow2(Mt)*(159*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 35*pow2(Mst2)*pow2(
        Sbeta)) + 424*Cbeta*MuSUSY*Sbeta*pow3(Mt) + 483*pow2(Sbeta)*pow3(s2t)*
        pow4(Mst2)) + pow2(Mst2)*pow3(Mst1)*(8*pow2(Mt)*pow2(s2t)*(167*pow2(
        MuSUSY)*(-1 + pow2(Sbeta)) + 42*pow2(Mst2)*pow2(Sbeta)) - 960*Cbeta*
        MuSUSY*s2t*Sbeta*pow3(Mt) + 668*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*pow3(
        s2t) + 880*pow2(Sbeta)*pow4(Mt) - 167*pow2(Sbeta)*pow4(Mst2)*pow4(s2t))
        + 8*s2t*Sbeta*(-3*Cbeta*MuSUSY*s2t + 58*Mt*Sbeta)*pow2(Mt)*pow6(Mst2))
        + 120*Mst1*pow2(Msq)*(32*Mt*pow4(Mst1)*(117*Cbeta*Mt*MuSUSY*Sbeta*pow2(
        Mst2)*pow2(s2t) + 3*s2t*pow2(Mt)*(96*pow2(MuSUSY)*(-1 + pow2(Sbeta)) -
        37*pow2(Mst2)*pow2(Sbeta)) + Cbeta*MuSUSY*Sbeta*pow3(Mt) - 6*pow2(
        Sbeta)*pow3(s2t)*pow4(Mst2)) - 32*Mt*pow2(Mst1)*pow2(Mst2)*(-99*Cbeta*
        Mt*MuSUSY*Sbeta*pow2(Mst2)*pow2(s2t) - 12*s2t*pow2(Mt)*(11*pow2(MuSUSY)
        *(-1 + pow2(Sbeta)) + 4*pow2(Mst2)*pow2(Sbeta)) + 65*Cbeta*MuSUSY*
        Sbeta*pow3(Mt) + 33*pow2(Sbeta)*pow3(s2t)*pow4(Mst2)) - 4*s2t*pow2(
        Mst2)*pow3(Mst1)*(-127*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*pow2(s2t) + 4*
        s2t*pow2(Mt)*(-29*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 21*pow2(Mst2)*pow2(
        Sbeta)) - 240*Cbeta*MuSUSY*Sbeta*pow3(Mt) + 49*pow2(Sbeta)*pow3(s2t)*
        pow4(Mst2)) + 3*Mst1*pow2(Sbeta)*pow4(Mst2)*(-48*pow2(Mst2)*pow2(Mt)*
        pow2(s2t) - 16*pow4(Mt) + 23*pow4(Mst2)*pow4(s2t)) - 128*s2t*pow2(
        Sbeta)*pow3(Mt)*pow6(Mst2))) + pow2(Mgl)*(10584000*(-1 + 2*z2)*pow2(
        Mst1)*pow2(s2t)*pow6(Msq)*(pow2(Mst1)*(8*pow2(Mt)*pow2(MuSUSY)*(-1 +
        pow2(Sbeta)) + pow2(s2t)*pow2(Sbeta)*pow4(Mst2)) - pow2(s2t)*pow2(
        Sbeta)*pow6(Mst2)) - 30*pow2(Mst1)*pow2(Mst2)*(8*Mt*pow2(Mst1)*(-1960*
        Cbeta*MuSUSY*s2t*Sbeta*pow2(Mt) + 12*Mt*pow2(s2t)*(-10511*pow2(MuSUSY)*
        (-1 + pow2(Sbeta)) + 1868*pow2(Mst2)*pow2(Sbeta)) + 338100*pow2(Sbeta)*
        pow3(Mt) - 47419*Cbeta*MuSUSY*Sbeta*pow2(Mst2)*pow3(s2t))*pow4(Mst2) -
        39200*Mt*pow2(Mst2)*pow3(Mst1)*(-6*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*
        pow2(s2t) - 4*s2t*pow2(Mt)*(pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 17*pow2(
        Mst2)*pow2(Sbeta)) + 90*Cbeta*MuSUSY*Sbeta*pow3(Mt) + 3*pow2(Sbeta)*
        pow3(s2t)*pow4(Mst2)) + pow2(Mst2)*pow4(Mst1)*(24*pow2(Mt)*pow2(s2t)*(-
        1717*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 22148*pow2(Mst2)*pow2(Sbeta)) -
        1078784*Cbeta*MuSUSY*s2t*Sbeta*pow3(Mt) + 483924*Cbeta*Mt*MuSUSY*Sbeta*
        pow2(Mst2)*pow3(s2t) + 1078000*pow2(Sbeta)*pow4(Mt) - 215819*pow2(
        Sbeta)*pow4(Mst2)*pow4(s2t)) - 58800*Mt*(83*Cbeta*Mt*MuSUSY*Sbeta*pow2(
        Mst2)*pow2(s2t) + 4*s2t*pow2(Mt)*(27*pow2(MuSUSY)*(-1 + pow2(Sbeta)) +
        7*pow2(Mst2)*pow2(Sbeta)) + 8*Cbeta*MuSUSY*Sbeta*pow3(Mt) - 29*pow2(
        Sbeta)*pow3(s2t)*pow4(Mst2))*pow5(Mst1) + 39200*Mst1*s2t*Sbeta*(-3*
        Cbeta*MuSUSY*s2t + 58*Mt*Sbeta)*pow2(Mt)*pow6(Mst2) + 16*pow2(Mt)*(
        21436*Cbeta*Mt*MuSUSY*s2t*Sbeta - 15647*pow2(MuSUSY)*pow2(s2t)*(-1 +
        pow2(Sbeta)) + 131493*pow2(Mt)*pow2(Sbeta))*pow6(Mst2)) - 11760*pow2(
        Msq)*pow2(Mst1)*(1600*Mt*pow2(Mst2)*pow3(Mst1)*(-30*Cbeta*Mt*MuSUSY*
        Sbeta*pow2(Mst2)*pow2(s2t) - 8*s2t*pow2(Mt)*(5*pow2(MuSUSY)*(-1 + pow2(
        Sbeta)) + 2*pow2(Mst2)*pow2(Sbeta)) + 21*Cbeta*MuSUSY*Sbeta*pow3(Mt) +
        10*pow2(Sbeta)*pow3(s2t)*pow4(Mst2)) + 2*s2t*pow2(Mst2)*pow4(Mst1)*(-
        3150*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*pow2(s2t) + 16*s2t*pow2(Mt)*(-
        482*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 125*pow2(Mst2)*pow2(Sbeta)) -
        6000*Cbeta*MuSUSY*Sbeta*pow3(Mt) + 991*pow2(Sbeta)*pow3(s2t)*pow4(Mst2)
        ) + pow2(Mst1)*pow4(Mst2)*(-8*pow2(Mt)*pow2(s2t)*(353*pow2(MuSUSY)*(-1
        + pow2(Sbeta)) - 668*pow2(Mst2)*pow2(Sbeta)) - 4000*Cbeta*MuSUSY*s2t*
        Sbeta*pow3(Mt) + 1628*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*pow3(s2t) +
        1200*pow2(Sbeta)*pow4(Mt) - 1167*pow2(Sbeta)*pow4(Mst2)*pow4(s2t)) +
        1600*Mt*(-51*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*pow2(s2t) - 3*s2t*pow2(
        Mt)*(36*pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 5*pow2(Mst2)*pow2(Sbeta)) +
        15*Cbeta*MuSUSY*Sbeta*pow3(Mt) + 7*pow2(Sbeta)*pow3(s2t)*pow4(Mst2))*
        pow5(Mst1) + 6400*Mst1*s2t*pow2(Sbeta)*pow3(Mt)*pow6(Mst2) - 16*Mt*
        Sbeta*(209*Mt*Sbeta*pow2(Mst2)*pow2(s2t) - 482*Sbeta*pow3(Mt) + 190*
        Cbeta*MuSUSY*pow2(Mst2)*pow3(s2t))*pow6(Mst2)) + pow4(Msq)*(705600*Mt*
        s2t*((208 - 288*z2 + 192*z3)*pow2(Mt) + (2135 + 384*z2 - 36*z3)*pow2(
        Mst2)*pow2(s2t))*pow2(Sbeta)*pow3(Mst1)*pow4(Mst2) + 14700*pow2(Mst1)*
        pow2(Sbeta)*pow4(Mst2)*(-32*(733 + 81*S2 + 18*z2)*pow2(Mst2)*pow2(Mt)*
        pow2(s2t) + 96*(149 + 24*z3)*pow4(Mt) + (14209 + 648*S2 - 1788*z2)*
        pow4(Mst2)*pow4(s2t)) - 78400*Mt*(3*Cbeta*Mt*MuSUSY*Sbeta*(22549 -
        1008*z2 + 7776*z3)*pow2(Mst2)*pow2(s2t) + 4*s2t*pow2(Mt)*(4*(10441 +
        612*z2 + 1863*z3)*pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 3*(2197 + 3396*z2)*
        pow2(Mst2)*pow2(Sbeta)) + 12*Cbeta*MuSUSY*Sbeta*(643 + 3180*z2 + 144*
        z3)*pow3(Mt) + 2*(-1667 + 2232*z2 - 4050*z3)*pow2(Sbeta)*pow3(s2t)*
        pow4(Mst2))*pow5(Mst1) - 3*(24*pow2(Mt)*pow2(s2t)*((72625277 +
        14504000*S2 + 19521600*z2 - 4939200*z3)*pow2(MuSUSY)*(-1 + pow2(Sbeta))
        - 8*(4731149 + 4900*S2 - 808500*z2)*pow2(Mst2)*pow2(Sbeta)) - 128*
        Cbeta*MuSUSY*s2t*Sbeta*(-23055979 + 294000*S2 + 1881600*z2)*pow3(Mt) +
        4*Cbeta*Mt*MuSUSY*Sbeta*(71872687 + 19521600*S2 + 36691200*z2)*pow2(
        Mst2)*pow3(s2t) + 128*(3291437 + 3013500*z2)*pow2(Sbeta)*pow4(Mt) + 3*(
        1737319 + 431200*S2 - 2018800*z2 - 4939200*z3)*pow2(Sbeta)*pow4(Mst2)*
        pow4(s2t))*pow6(Mst1) + 4233600*Mt*Sbeta*(2*Mt*Sbeta*pow2(Mst2)*pow2(
        s2t) - 4*Sbeta*pow3(Mt) + Cbeta*MuSUSY*pow2(Mst2)*pow3(s2t))*pow6(Mst2)
        - 784*pow4(Mst1)*(12*pow2(Mst2)*pow2(Mt)*pow2(s2t)*((372457 + 61200*S2
        + 55800*z2 - 37800*z3)*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 6*(-33703 +
        375*S2 - 2000*z2)*pow2(Mst2)*pow2(Sbeta)) - 96*Cbeta*MuSUSY*s2t*Sbeta*(
        -43417 + 1575*S2 - 2775*z2)*pow2(Mst2)*pow3(Mt) + 3*Cbeta*Mt*MuSUSY*
        Sbeta*(393289 + 106200*S2 + 156300*z2 - 75600*z3)*pow3(s2t)*pow4(Mst2)
        + 16*(1800*(87 + 16*z2 + 3*z3)*pow2(MuSUSY)*(-1 + pow2(Sbeta)) - (54223
        + 34200*z2)*pow2(Mst2)*pow2(Sbeta))*pow4(Mt) - 18*(1661 + 3750*S2 +
        8375*z2 - 3150*z3)*pow2(Sbeta)*pow4(s2t)*pow6(Mst2)) + 58800*Mt*s2t*
        Sbeta*(Cbeta*Mt*MuSUSY*s2t*(-101509 + 14256*z2 - 31104*z3) + 32*Sbeta*(
        2903 + 2196*z2)*pow2(Mt) + 3*Sbeta*(1257 - 1136*z2)*pow2(Mst2)*pow2(
        s2t))*pow7(Mst1)))))/(529200.*pow2(Mgl)*pow2(Mst1)*pow2(Sbeta)*pow4(
        Msq)*pow4(Mst2)*pow4(Mt)) - (log(pow2(Mst1)/pow2(Msq))*(2450*Dmglst1*
        Mgl*pow2(Mst1)*(2*pow4(Msq)*(16*Mst1*Mt*s2t*(-251*pow2(Mt) + 9*(9 + z2)
        *pow2(Mst2)*pow2(s2t))*pow2(Sbeta)*pow4(Mst2) - 16*Mt*pow3(Mst1)*(-27*
        Cbeta*Mt*MuSUSY*Sbeta*(7 + 10*z2)*pow2(Mst2)*pow2(s2t) - 18*s2t*pow2(
        Mt)*(2*(-2 + 9*z2)*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + (7 + 24*z2)*pow2(
        Mst2)*pow2(Sbeta)) + 2*Cbeta*MuSUSY*Sbeta*(-151 + 216*z2)*pow3(Mt) + 9*
        (16 + 11*z2)*pow2(Sbeta)*pow3(s2t)*pow4(Mst2)) + pow2(Sbeta)*pow4(Mst2)
        *(192*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 80*pow4(Mt) - 45*pow4(Mst2)*pow4(
        s2t)) + 3*pow2(Mst1)*pow2(Mst2)*pow2(Sbeta)*(152*pow2(Mst2)*pow2(Mt)*
        pow2(s2t) + 432*pow4(Mt) + 9*(7 + 4*z2)*pow4(Mst2)*pow4(s2t)) - 9*pow4(
        Mst1)*(4*Cbeta*Mt*MuSUSY*s2t*Sbeta*(72*pow2(Mt) - 11*pow2(Mst2)*pow2(
        s2t)) - 24*pow2(Mt)*pow2(s2t)*(-4*z2*pow2(MuSUSY)*(-1 + pow2(Sbeta)) +
        3*pow2(Mst2)*pow2(Sbeta)) + 32*(9 + 20*z2)*pow2(Sbeta)*pow4(Mt) + 3*(9
        + 4*z2)*pow2(Sbeta)*pow4(Mst2)*pow4(s2t)) + 36*Mt*s2t*Sbeta*(15*Cbeta*
        Mt*MuSUSY*s2t*(-3 + 8*z2) + 40*Sbeta*(-1 + 12*z2)*pow2(Mt) + 43*Sbeta*
        pow2(Mst2)*pow2(s2t))*pow5(Mst1)) - Mst1*pow2(Mst2)*(24*Mst1*pow2(Mt)*(
        16*pow2(Mt) + 3*pow2(Mst2)*pow2(s2t))*pow2(Sbeta)*pow4(Mst2) - 4*Mt*
        s2t*pow4(Mst1)*(543*Cbeta*Mt*MuSUSY*s2t*Sbeta*pow2(Mst2) + 80*pow2(Mt)*
        (9*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + pow2(Mst2)*pow2(Sbeta)) - 183*
        pow2(s2t)*pow2(Sbeta)*pow4(Mst2)) - 4*Mt*pow2(Mst1)*pow2(Mst2)*(-6*
        Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*pow2(s2t) - 4*s2t*pow2(Mt)*(pow2(
        MuSUSY)*(-1 + pow2(Sbeta)) + 12*pow2(Mst2)*pow2(Sbeta)) - 4*Cbeta*
        MuSUSY*Sbeta*pow3(Mt) + 3*pow2(Sbeta)*pow3(s2t)*pow4(Mst2)) + pow2(
        Mst2)*pow3(Mst1)*(-24*pow2(Mt)*pow2(s2t)*(pow2(MuSUSY)*(-1 + pow2(
        Sbeta)) + 27*pow2(Mst2)*pow2(Sbeta)) + 1152*Cbeta*MuSUSY*s2t*Sbeta*
        pow3(Mt) - 12*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*pow3(s2t) + 2192*pow2(
        Sbeta)*pow4(Mt) + 3*pow2(Sbeta)*pow4(Mst2)*pow4(s2t)) + 4*s2t*Sbeta*(-
        3*Cbeta*MuSUSY*s2t + 32*Mt*Sbeta)*pow2(Mt)*pow6(Mst2)) + 8*Mst1*pow2(
        Msq)*(-8*Mt*s2t*pow2(Mst1)*pow2(Mst2)*(-54*Cbeta*Mt*MuSUSY*s2t*Sbeta*
        pow2(Mst2) + pow2(Mt)*(-72*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 7*pow2(
        Mst2)*pow2(Sbeta)) + 18*pow2(s2t)*pow2(Sbeta)*pow4(Mst2)) + 144*pow4(
        Mst1)*(4*s2t*pow2(MuSUSY)*(-1 + pow2(Sbeta))*pow3(Mt) + Mt*pow2(Sbeta)*
        pow3(s2t)*pow4(Mst2)) - 6*Mst1*pow2(Sbeta)*pow4(Mst2)*(-20*pow2(Mst2)*
        pow2(Mt)*pow2(s2t) + 24*pow4(Mt) + pow4(Mst2)*pow4(s2t)) - 3*pow2(Mst2)
        *pow3(Mst1)*(40*pow2(Mt)*pow2(s2t)*(-(pow2(MuSUSY)*(-1 + pow2(Sbeta)))
        + 2*pow2(Mst2)*pow2(Sbeta)) - 80*Cbeta*MuSUSY*s2t*Sbeta*pow3(Mt) - 12*
        Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*pow3(s2t) + 80*pow2(Sbeta)*pow4(Mt) +
        pow2(Sbeta)*pow4(Mst2)*pow4(s2t)) + 56*s2t*pow2(Sbeta)*pow3(Mt)*pow6(
        Mst2))) + 49*pow2(Dmglst1)*pow2(Mst1)*(2*pow4(Msq)*(400*Mst1*Mt*s2t*(-
        838*pow2(Mt) + 3*(71 + 6*z2)*pow2(Mst2)*pow2(s2t))*pow2(Sbeta)*pow4(
        Mst2) - 400*Mt*pow3(Mst1)*(-54*Cbeta*Mt*MuSUSY*Sbeta*(7 + 10*z2)*pow2(
        Mst2)*pow2(s2t) - 72*s2t*pow2(Mt)*((-2 + 9*z2)*pow2(MuSUSY)*(-1 + pow2(
        Sbeta)) + 2*(5 + 12*z2)*pow2(Mst2)*pow2(Sbeta)) + 4*Cbeta*MuSUSY*Sbeta*
        (-151 + 216*z2)*pow3(Mt) + 3*(251 + 138*z2)*pow2(Sbeta)*pow3(s2t)*pow4(
        Mst2)) + pow2(Sbeta)*pow4(Mst2)*(31200*pow2(Mst2)*pow2(Mt)*pow2(s2t) -
        32096*pow4(Mt) - 5025*pow4(Mst2)*pow4(s2t)) + 75*pow2(Mst1)*pow2(Mst2)*
        pow2(Sbeta)*(232*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 1296*pow4(Mt) + (329 +
        108*z2)*pow4(Mst2)*pow4(s2t)) - 75*pow4(Mst1)*(24*pow2(Mt)*pow2(s2t)*(
        4*(19 + 9*z2)*pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 103*pow2(Mst2)*pow2(
        Sbeta)) + 6240*Cbeta*MuSUSY*s2t*Sbeta*pow3(Mt) + 132*Cbeta*Mt*MuSUSY*
        Sbeta*pow2(Mst2)*pow3(s2t) + 480*(13 + 20*z2)*pow2(Sbeta)*pow4(Mt) + (
        229 + 108*z2)*pow2(Sbeta)*pow4(Mst2)*pow4(s2t)) + 1800*Mt*s2t*Sbeta*(3*
        Cbeta*Mt*MuSUSY*s2t*(61 + 88*z2) + 480*Sbeta*(1 + 3*z2)*pow2(Mt) + 59*
        Sbeta*pow2(Mst2)*pow2(s2t))*pow5(Mst1)) - 25*Mst1*pow2(Mst2)*(72*Mst1*
        pow2(Mt)*(16*pow2(Mt) + 3*pow2(Mst2)*pow2(s2t))*pow2(Sbeta)*pow4(Mst2)
        - 8*Mt*s2t*pow4(Mst1)*(1623*Cbeta*Mt*MuSUSY*s2t*Sbeta*pow2(Mst2) + 16*
        pow2(Mt)*(135*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 8*pow2(Mst2)*pow2(
        Sbeta)) - 543*pow2(s2t)*pow2(Sbeta)*pow4(Mst2)) - 8*Mt*pow2(Mst1)*pow2(
        Mst2)*(-6*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*pow2(s2t) - 4*s2t*pow2(Mt)*(
        pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 24*pow2(Mst2)*pow2(Sbeta)) - 4*Cbeta*
        MuSUSY*Sbeta*pow3(Mt) + 3*pow2(Sbeta)*pow3(s2t)*pow4(Mst2)) + pow2(
        Mst2)*pow3(Mst1)*(-24*pow2(Mt)*pow2(s2t)*(47*pow2(MuSUSY)*(-1 + pow2(
        Sbeta)) + 165*pow2(Mst2)*pow2(Sbeta)) + 7488*Cbeta*MuSUSY*s2t*Sbeta*
        pow3(Mt) - 564*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*pow3(s2t) + 12208*pow2(
        Sbeta)*pow4(Mt) + 141*pow2(Sbeta)*pow4(Mst2)*pow4(s2t)) + 8*s2t*Sbeta*(
        -3*Cbeta*MuSUSY*s2t + 32*Mt*Sbeta)*pow2(Mt)*pow6(Mst2)) + 200*Mst1*
        pow2(Msq)*(576*Mt*s2t*pow4(Mst1)*(4*pow2(Mt)*pow2(MuSUSY)*(-1 + pow2(
        Sbeta)) + pow2(s2t)*pow2(Sbeta)*pow4(Mst2)) - 16*Mt*s2t*pow2(Mst1)*
        pow2(Mst2)*(-54*Cbeta*Mt*MuSUSY*s2t*Sbeta*pow2(Mst2) + pow2(Mt)*(-72*
        pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 7*pow2(Mst2)*pow2(Sbeta)) + 36*pow2(
        s2t)*pow2(Sbeta)*pow4(Mst2)) - 24*Mst1*pow2(Sbeta)*pow4(Mst2)*(-25*
        pow2(Mst2)*pow2(Mt)*pow2(s2t) + 38*pow4(Mt) + 2*pow4(Mst2)*pow4(s2t)) +
        3*pow2(Mst2)*pow3(Mst1)*(200*pow2(Mt)*pow2(s2t)*(pow2(MuSUSY)*(-1 +
        pow2(Sbeta)) - 2*pow2(Mst2)*pow2(Sbeta)) + 400*Cbeta*MuSUSY*s2t*Sbeta*
        pow3(Mt) + 36*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*pow3(s2t) - 400*pow2(
        Sbeta)*pow4(Mt) + 7*pow2(Sbeta)*pow4(Mst2)*pow4(s2t)) + 112*s2t*pow2(
        Sbeta)*pow3(Mt)*pow6(Mst2))) + pow2(Mgl)*(1225*pow2(Mst1)*pow4(Msq)*(
        576*Mst1*Mt*s2t*(-11*pow2(Mt) + (6 + z2)*pow2(Mst2)*pow2(s2t))*pow2(
        Sbeta)*pow4(Mst2) - 576*Mt*pow3(Mst1)*(-3*Cbeta*Mt*MuSUSY*Sbeta*(1 + 2*
        z2)*pow2(Mst2)*pow2(s2t) - 2*s2t*pow2(Mt)*(2*(-5 + z2)*pow2(MuSUSY)*(-1
        + pow2(Sbeta)) + (1 + 8*z2)*pow2(Mst2)*pow2(Sbeta)) + 2*Cbeta*MuSUSY*
        Sbeta*(-7 + 8*z2)*pow3(Mt) + (7 + 3*z2)*pow2(Sbeta)*pow3(s2t)*pow4(
        Mst2)) - 2*pow2(Sbeta)*pow4(Mst2)*(72*pow2(Mst2)*pow2(Mt)*pow2(s2t) +
        32*(-41 + 36*z2)*pow4(Mt) + 3*(-37 + 6*z2)*pow4(Mst2)*pow4(s2t)) - 3*
        pow4(Mst1)*(24*pow2(Mt)*pow2(s2t)*((-15 + 16*z2)*pow2(MuSUSY)*(-1 +
        pow2(Sbeta)) + 10*pow2(Mst2)*pow2(Sbeta)) + 384*Cbeta*MuSUSY*s2t*Sbeta*
        pow3(Mt) - 348*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*pow3(s2t) + 384*(1 + 5*
        z2)*pow2(Sbeta)*pow4(Mt) + 5*(11 + 12*z2)*pow2(Sbeta)*pow4(Mst2)*pow4(
        s2t)) + 6*pow2(Mst1)*pow2(Mst2)*(-24*pow2(Mt)*pow2(s2t)*((7 + 8*z2)*
        pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 10*pow2(Mst2)*pow2(Sbeta)) - 432*
        Cbeta*MuSUSY*s2t*Sbeta*pow3(Mt) - 8*Cbeta*Mt*MuSUSY*Sbeta*(-8 + 15*z2)*
        pow2(Mst2)*pow3(s2t) + 432*pow2(Sbeta)*pow4(Mt) + (-53 + 36*z2)*pow2(
        Sbeta)*pow4(Mst2)*pow4(s2t)) + 144*Mt*s2t*Sbeta*(3*Cbeta*Mt*MuSUSY*s2t*
        (-7 + 8*z2) + 8*Sbeta*(-7 + 12*z2)*pow2(Mt) + 11*Sbeta*pow2(Mst2)*pow2(
        s2t))*pow5(Mst1)) + 88200*Sbeta*pow2(Mst2)*pow6(Msq)*(4*Cbeta*Mt*
        MuSUSY*pow3(s2t)*pow4(Mst1) + Sbeta*(pow2(-4*Mst2*pow2(Mt) + pow2(s2t)*
        pow3(Mst2)) - pow4(Mst1)*(8*pow2(Mt)*pow2(s2t) + pow2(Mst2)*pow4(s2t))
        + pow2(Mst1)*(16*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 16*pow4(Mt) - pow4(
        Mst2)*pow4(s2t)) + pow4(s2t)*pow6(Mst1))) - 2*pow2(Mst1)*pow2(Mst2)*(
        Mt*pow2(Mst1)*(-88200*Cbeta*MuSUSY*s2t*Sbeta*pow2(Mt) + 4*Mt*pow2(s2t)*
        (-13109*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 14865*pow2(Mst2)*pow2(Sbeta))
        + 235200*pow2(Sbeta)*pow3(Mt) - 16647*Cbeta*MuSUSY*Sbeta*pow2(Mst2)*
        pow3(s2t))*pow4(Mst2) - 4900*Mt*pow2(Mst2)*pow3(Mst1)*(-6*Cbeta*Mt*
        MuSUSY*Sbeta*pow2(Mst2)*pow2(s2t) - 4*s2t*pow2(Mt)*(pow2(MuSUSY)*(-1 +
        pow2(Sbeta)) + 4*pow2(Mst2)*pow2(Sbeta)) + 12*Cbeta*MuSUSY*Sbeta*pow3(
        Mt) + 3*pow2(Sbeta)*pow3(s2t)*pow4(Mst2)) + pow2(Mst2)*pow4(Mst1)*(2*
        pow2(Mt)*pow2(s2t)*(12479*pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 87870*pow2(
        Mst2)*pow2(Sbeta)) + 263280*Cbeta*MuSUSY*s2t*Sbeta*pow3(Mt) + 38697*
        Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*pow3(s2t) + 882836*pow2(Sbeta)*pow4(
        Mt) - 13836*pow2(Sbeta)*pow4(Mst2)*pow4(s2t)) - 14700*Mt*s2t*(37*Cbeta*
        Mt*MuSUSY*s2t*Sbeta*pow2(Mst2) + 16*pow2(Mt)*(3*pow2(MuSUSY)*(-1 +
        pow2(Sbeta)) + pow2(Mst2)*pow2(Sbeta)) - 13*pow2(s2t)*pow2(Sbeta)*pow4(
        Mst2))*pow5(Mst1) + 4900*Mst1*s2t*Sbeta*(-3*Cbeta*MuSUSY*s2t + 32*Mt*
        Sbeta)*pow2(Mt)*pow6(Mst2) + 2*pow2(Mt)*(15360*Cbeta*Mt*MuSUSY*s2t*
        Sbeta - 9571*pow2(MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta)) + 98418*pow2(Mt)
        *pow2(Sbeta))*pow6(Mst2)) + 196*pow2(Msq)*pow2(Mst1)*(-800*Mt*s2t*pow2(
        Mst2)*pow3(Mst1)*(-18*Cbeta*Mt*MuSUSY*s2t*Sbeta*pow2(Mst2) + pow2(Mt)*(
        -24*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 7*pow2(Mst2)*pow2(Sbeta)) + 6*
        pow2(s2t)*pow2(Sbeta)*pow4(Mst2)) + pow2(Mst2)*pow4(Mst1)*(8*pow2(Mt)*
        pow2(s2t)*(375*pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 641*pow2(Mst2)*pow2(
        Sbeta)) + 6000*Cbeta*MuSUSY*s2t*Sbeta*pow3(Mt) + 1080*Cbeta*Mt*MuSUSY*
        Sbeta*pow2(Mst2)*pow3(s2t) - 6000*pow2(Sbeta)*pow4(Mt) - 345*pow2(
        Sbeta)*pow4(Mst2)*pow4(s2t)) + pow2(Mst1)*pow4(Mst2)*(8*pow2(Mt)*pow2(
        s2t)*(105*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 157*pow2(Mst2)*pow2(Sbeta))
        - 4256*Cbeta*MuSUSY*s2t*Sbeta*pow3(Mt) - 300*Cbeta*Mt*MuSUSY*Sbeta*
        pow2(Mst2)*pow3(s2t) - 1856*pow2(Sbeta)*pow4(Mt) + 255*pow2(Sbeta)*
        pow4(Mst2)*pow4(s2t)) + 4800*(4*s2t*pow2(MuSUSY)*(-1 + pow2(Sbeta))*
        pow3(Mt) + Mt*pow2(Sbeta)*pow3(s2t)*pow4(Mst2))*pow5(Mst1) + 5600*Mst1*
        s2t*pow2(Sbeta)*pow3(Mt)*pow6(Mst2) + 8*Mt*Sbeta*(109*Mt*Sbeta*pow2(
        Mst2)*pow2(s2t) - 82*Sbeta*pow3(Mt) + 90*Cbeta*MuSUSY*pow2(Mst2)*pow3(
        s2t))*pow6(Mst2))) + 22050*pow2(Mst1)*pow2(log(pow2(Mst1)/pow2(Mst2)))*
        pow4(Msq)*(pow2(Mgl)*(-16*Mt*pow3(Mst1)*(-6*Cbeta*Mt*MuSUSY*Sbeta*pow2(
        Mst2)*pow2(s2t) - 4*s2t*pow2(Mt)*(pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 4*
        pow2(Mst2)*pow2(Sbeta)) + 16*Cbeta*MuSUSY*Sbeta*pow3(Mt) + 3*pow2(
        Sbeta)*pow3(s2t)*pow4(Mst2)) + 2*s2t*pow2(Mst1)*pow2(Mst2)*(-26*Cbeta*
        Mt*MuSUSY*Sbeta*pow2(Mst2)*pow2(s2t) - 64*s2t*pow2(Mt)*pow2(MuSUSY)*(-1
        + pow2(Sbeta)) + 16*Cbeta*MuSUSY*Sbeta*pow3(Mt) + 5*pow2(Sbeta)*pow3(
        s2t)*pow4(Mst2)) - pow4(Mst1)*(192*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*(-1
        + pow2(Sbeta)) - 32*Cbeta*MuSUSY*s2t*Sbeta*pow3(Mt) + 32*Cbeta*Mt*
        MuSUSY*Sbeta*pow2(Mst2)*pow3(s2t) + 160*pow2(Sbeta)*pow4(Mt) + 5*pow2(
        Sbeta)*pow4(Mst2)*pow4(s2t)) + 96*s2t*Sbeta*(Cbeta*MuSUSY*s2t + 4*Mt*
        Sbeta)*pow2(Mt)*pow5(Mst1) + pow2(s2t)*(-16*pow2(Mt) + 3*pow2(Mst2)*
        pow2(s2t))*pow2(Sbeta)*pow6(Mst2) + 16*Mst1*Mt*pow2(Sbeta)*pow3(s2t)*
        pow6(Mst2)) + 4*Dmglst1*Mgl*Mst1*(120*s2t*Sbeta*(Cbeta*MuSUSY*s2t + 4*
        Mt*Sbeta)*pow2(Mt)*pow4(Mst1) - 4*Mt*pow2(Mst1)*(-30*Cbeta*Mt*MuSUSY*
        Sbeta*pow2(Mst2)*pow2(s2t) - 12*s2t*pow2(Mt)*(3*pow2(MuSUSY)*(-1 +
        pow2(Sbeta)) + 4*pow2(Mst2)*pow2(Sbeta)) + 48*Cbeta*MuSUSY*Sbeta*pow3(
        Mt) + 11*pow2(Sbeta)*pow3(s2t)*pow4(Mst2)) - pow3(Mst1)*(24*pow2(Mt)*
        pow2(MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta)) + 160*pow2(Sbeta)*pow4(Mt) +
        3*pow2(Sbeta)*pow4(Mst2)*pow4(s2t)) + 4*Mt*pow2(Sbeta)*pow3(s2t)*pow6(
        Mst2) + 3*Mst1*pow2(Sbeta)*pow4(s2t)*pow6(Mst2)) + 2*Mst1*pow2(Dmglst1)
        *(48*s2t*Sbeta*(11*Cbeta*MuSUSY*s2t + 60*Mt*Sbeta)*pow2(Mt)*pow4(Mst1)
        - 8*Mt*pow2(Mst1)*(-30*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*pow2(s2t) - 12*
        s2t*pow2(Mt)*(3*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 8*pow2(Mst2)*pow2(
        Sbeta)) + 48*Cbeta*MuSUSY*Sbeta*pow3(Mt) + 23*pow2(Sbeta)*pow3(s2t)*
        pow4(Mst2)) - pow3(Mst1)*(72*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*(-1 +
        pow2(Sbeta)) + 800*pow2(Sbeta)*pow4(Mt) + 9*pow2(Sbeta)*pow4(Mst2)*
        pow4(s2t)) + 8*Mt*pow2(Sbeta)*pow3(s2t)*pow6(Mst2) + 9*Mst1*pow2(Sbeta)
        *pow4(s2t)*pow6(Mst2))) + 105*log(pow2(Mst1)/pow2(Mst2))*pow2(Mst1)*(
        140*Dmglst1*Mgl*Mst1*(16*Sbeta*pow2(Mst1)*(-7*Cbeta*Mt*MuSUSY + 6*s2t*
        Sbeta*pow2(Mst2))*pow3(Mt)*pow4(Mst2) - 4*Mt*pow2(Mst2)*pow4(Mst1)*(-
        45*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*pow2(s2t) - 60*s2t*pow2(Mt)*(pow2(
        MuSUSY)*(-1 + pow2(Sbeta)) + pow2(Mst2)*pow2(Sbeta)) + 88*Cbeta*MuSUSY*
        Sbeta*pow3(Mt) + 15*pow2(Sbeta)*pow3(s2t)*pow4(Mst2)) + 5*pow3(Mst1)*
        pow4(Mst2)*(-8*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta)) - 4*
        Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*pow3(s2t) + 16*pow2(Sbeta)*pow4(Mt) +
        pow2(Sbeta)*pow4(Mst2)*pow4(s2t)) + 12*pow4(Msq)*(2*Mt*s2t*Sbeta*(117*
        Cbeta*Mt*MuSUSY*s2t + 388*Sbeta*pow2(Mt) - 33*Sbeta*pow2(Mst2)*pow2(
        s2t))*pow4(Mst1) + 3*Mst1*pow2(s2t)*(4*pow2(Mt) + pow2(Mst2)*pow2(s2t))
        *pow2(Sbeta)*pow4(Mst2) + 2*Mt*s2t*(-8*pow2(Mt) + 7*pow2(Mst2)*pow2(
        s2t))*pow2(Sbeta)*pow4(Mst2) + 2*Mt*pow2(Mst1)*(18*Cbeta*Mt*MuSUSY*
        Sbeta*pow2(Mst2)*pow2(s2t) + 4*s2t*pow2(Mt)*(pow2(MuSUSY) + 20*pow2(
        Mst2)*pow2(Sbeta) - pow2(MuSUSY)*pow2(Sbeta)) - 44*Cbeta*MuSUSY*Sbeta*
        pow3(Mt) - 13*pow2(Sbeta)*pow3(s2t)*pow4(Mst2)) + pow3(Mst1)*(4*pow2(
        Mt)*pow2(s2t)*(-18*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 7*pow2(Mst2)*pow2(
        Sbeta)) - 80*Cbeta*MuSUSY*s2t*Sbeta*pow3(Mt) - 24*Cbeta*Mt*MuSUSY*
        Sbeta*pow2(Mst2)*pow3(s2t) - 232*pow2(Sbeta)*pow4(Mt) + 3*pow2(Sbeta)*
        pow4(Mst2)*pow4(s2t))) + 48*Mst1*pow2(Sbeta)*pow4(Mt)*pow6(Mst2) +
        pow2(Msq)*(-8*pow2(Mst2)*pow2(s2t)*pow3(Mst1)*(-17*Cbeta*Mt*MuSUSY*s2t*
        Sbeta*pow2(Mst2) - 28*pow2(Mt)*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 5*
        pow2(s2t)*pow2(Sbeta)*pow4(Mst2)) - 32*Mt*pow4(Mst1)*(-18*Cbeta*Mt*
        MuSUSY*Sbeta*pow2(Mst2)*pow2(s2t) - 36*s2t*pow2(Mt)*pow2(MuSUSY)*(-1 +
        pow2(Sbeta)) + 14*Cbeta*MuSUSY*Sbeta*pow3(Mt) + 3*pow2(Sbeta)*pow3(s2t)
        *pow4(Mst2)) - 32*Mt*pow2(Mst1)*pow2(Mst2)*(-9*Cbeta*Mt*MuSUSY*Sbeta*
        pow2(Mst2)*pow2(s2t) - 12*s2t*pow2(Mt)*(pow2(MuSUSY)*(-1 + pow2(Sbeta))
        + pow2(Mst2)*pow2(Sbeta)) + 14*Cbeta*MuSUSY*Sbeta*pow3(Mt) + 3*pow2(
        Sbeta)*pow3(s2t)*pow4(Mst2)) + 6*Mst1*pow2(Sbeta)*pow4(Mst2)*(16*pow4(
        Mt) + pow4(Mst2)*pow4(s2t)) + 64*s2t*pow2(Sbeta)*pow3(Mt)*pow6(Mst2)) +
        16*s2t*pow2(Sbeta)*pow3(Mt)*pow8(Mst2)) + 70*Mst1*pow2(Dmglst1)*(32*
        Sbeta*pow2(Mst1)*(-7*Cbeta*Mt*MuSUSY + 12*s2t*Sbeta*pow2(Mst2))*pow3(
        Mt)*pow4(Mst2) - 8*Mt*pow2(Mst2)*pow4(Mst1)*(-135*Cbeta*Mt*MuSUSY*
        Sbeta*pow2(Mst2)*pow2(s2t) - 180*s2t*pow2(Mt)*(pow2(MuSUSY)*(-1 + pow2(
        Sbeta)) + pow2(Mst2)*pow2(Sbeta)) + 232*Cbeta*MuSUSY*Sbeta*pow3(Mt) +
        45*pow2(Sbeta)*pow3(s2t)*pow4(Mst2)) + 25*pow3(Mst1)*pow4(Mst2)*(-8*
        pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta)) - 4*Cbeta*Mt*MuSUSY*
        Sbeta*pow2(Mst2)*pow3(s2t) + 16*pow2(Sbeta)*pow4(Mt) + pow2(Sbeta)*
        pow4(Mst2)*pow4(s2t)) + 12*pow4(Msq)*(12*Mt*s2t*Sbeta*(157*Cbeta*Mt*
        MuSUSY*s2t + 428*Sbeta*pow2(Mt) - 37*Sbeta*pow2(Mst2)*pow2(s2t))*pow4(
        Mst1) + 9*Mst1*pow2(s2t)*(4*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow2(
        Sbeta)*pow4(Mst2) + 4*Mt*s2t*(-8*pow2(Mt) + 7*pow2(Mst2)*pow2(s2t))*
        pow2(Sbeta)*pow4(Mst2) - 4*Mt*pow2(Mst1)*(-18*Cbeta*Mt*MuSUSY*Sbeta*
        pow2(Mst2)*pow2(s2t) + 4*s2t*pow2(Mt)*(pow2(MuSUSY)*(-1 + pow2(Sbeta))
        - 52*pow2(Mst2)*pow2(Sbeta)) + 44*Cbeta*MuSUSY*Sbeta*pow3(Mt) + 53*
        pow2(Sbeta)*pow3(s2t)*pow4(Mst2)) - 4*pow3(Mst1)*(pow2(Mt)*pow2(s2t)*(
        84*pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 41*pow2(Mst2)*pow2(Sbeta)) + 100*
        Cbeta*MuSUSY*s2t*Sbeta*pow3(Mt) + 33*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*
        pow3(s2t) + 330*pow2(Sbeta)*pow4(Mt) - 6*pow2(Sbeta)*pow4(Mst2)*pow4(
        s2t))) + 144*Mst1*pow2(Sbeta)*pow4(Mt)*pow6(Mst2) - 2*pow2(Msq)*(4*
        pow2(Mst2)*pow2(s2t)*pow3(Mst1)*(-91*Cbeta*Mt*MuSUSY*s2t*Sbeta*pow2(
        Mst2) - 164*pow2(Mt)*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 25*pow2(s2t)*
        pow2(Sbeta)*pow4(Mst2)) + 64*Mt*pow4(Mst1)*(-18*Cbeta*Mt*MuSUSY*Sbeta*
        pow2(Mst2)*pow2(s2t) - 36*s2t*pow2(Mt)*pow2(MuSUSY)*(-1 + pow2(Sbeta))
        + 13*Cbeta*MuSUSY*Sbeta*pow3(Mt) + 3*pow2(Sbeta)*pow3(s2t)*pow4(Mst2))
        + 32*Mt*pow2(Mst1)*pow2(Mst2)*(-9*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*
        pow2(s2t) - 12*s2t*pow2(Mt)*(pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 2*pow2(
        Mst2)*pow2(Sbeta)) + 14*Cbeta*MuSUSY*Sbeta*pow3(Mt) + 6*pow2(Sbeta)*
        pow3(s2t)*pow4(Mst2)) - 9*Mst1*pow2(Sbeta)*pow4(Mst2)*(16*pow4(Mt) +
        pow4(Mst2)*pow4(s2t)) - 64*s2t*pow2(Sbeta)*pow3(Mt)*pow6(Mst2)) + 32*
        s2t*pow2(Sbeta)*pow3(Mt)*pow8(Mst2)) + pow2(Mgl)*(2240*Sbeta*(-3*Cbeta*
        Mt*MuSUSY + 2*s2t*Sbeta*pow2(Mst2))*pow3(Mst1)*pow3(Mt)*pow4(Mst2) +
        pow4(Mst1)*pow4(Mst2)*(72*pow2(Mt)*pow2(s2t)*(3*pow2(MuSUSY)*(-1 +
        pow2(Sbeta)) - 7*pow2(Mst2)*pow2(Sbeta)) - 112*Cbeta*MuSUSY*s2t*Sbeta*
        pow3(Mt) - 308*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*pow3(s2t) + 2800*pow2(
        Sbeta)*pow4(Mt) + 163*pow2(Sbeta)*pow4(Mst2)*pow4(s2t)) - 1680*Mt*pow2(
        Mst2)*(-3*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*pow2(s2t) - 4*s2t*pow2(Mt)*(
        pow2(MuSUSY)*(-1 + pow2(Sbeta)) + pow2(Mst2)*pow2(Sbeta)) + 8*Cbeta*
        MuSUSY*Sbeta*pow3(Mt) + pow2(Sbeta)*pow3(s2t)*pow4(Mst2))*pow5(Mst1) +
        70*pow4(Msq)*(48*Mst1*Mt*s2t*(-8*pow2(Mt) + 7*pow2(Mst2)*pow2(s2t))*
        pow2(Sbeta)*pow4(Mst2) + 48*Mt*pow3(Mst1)*(-30*Cbeta*Mt*MuSUSY*Sbeta*
        pow2(Mst2)*pow2(s2t) + 4*s2t*pow2(Mt)*(-17*pow2(MuSUSY)*(-1 + pow2(
        Sbeta)) + 4*pow2(Mst2)*pow2(Sbeta)) + 4*Cbeta*MuSUSY*Sbeta*pow3(Mt) +
        3*pow2(Sbeta)*pow3(s2t)*pow4(Mst2)) + 4*s2t*pow2(Mst1)*pow2(Mst2)*(-59*
        Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*pow2(s2t) + 4*s2t*pow2(Mt)*(-47*pow2(
        MuSUSY)*(-1 + pow2(Sbeta)) + 16*pow2(Mst2)*pow2(Sbeta)) - 72*Cbeta*
        MuSUSY*Sbeta*pow3(Mt) + 6*pow2(Sbeta)*pow3(s2t)*pow4(Mst2)) + pow2(
        Sbeta)*pow4(Mst2)*(-112*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 80*pow4(Mt) +
        35*pow4(Mst2)*pow4(s2t)) - 4*pow4(Mst1)*(2*pow2(Mt)*pow2(s2t)*(121*
        pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 12*pow2(Mst2)*pow2(Sbeta)) + 120*
        Cbeta*MuSUSY*s2t*Sbeta*pow3(Mt) + 27*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*
        pow3(s2t) + 288*pow2(Sbeta)*pow4(Mt) + 8*pow2(Sbeta)*pow4(Mst2)*pow4(
        s2t)) + 48*Mt*s2t*Sbeta*(-15*Cbeta*Mt*MuSUSY*s2t + 68*Sbeta*pow2(Mt) -
        5*Sbeta*pow2(Mst2)*pow2(s2t))*pow5(Mst1)) + 8*Mt*pow2(Mst1)*(-140*
        Cbeta*MuSUSY*s2t*Sbeta*pow2(Mt) + Mt*pow2(s2t)*(104*pow2(MuSUSY)*(-1 +
        pow2(Sbeta)) + 33*pow2(Mst2)*pow2(Sbeta)) + 420*pow2(Sbeta)*pow3(Mt) +
        43*Cbeta*MuSUSY*Sbeta*pow2(Mst2)*pow3(s2t))*pow6(Mst2) - 28*pow2(Msq)*(
        pow2(Mst2)*pow2(s2t)*pow4(Mst1)*(-160*Cbeta*Mt*MuSUSY*s2t*Sbeta*pow2(
        Mst2) - 4*pow2(Mt)*(88*pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 5*pow2(Mst2)*
        pow2(Sbeta)) + 41*pow2(s2t)*pow2(Sbeta)*pow4(Mst2)) + 160*Mt*pow2(Mst2)
        *pow3(Mst1)*(-3*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*pow2(s2t) - 4*s2t*
        pow2(Mt)*(pow2(MuSUSY)*(-1 + pow2(Sbeta)) + pow2(Mst2)*pow2(Sbeta)) +
        6*Cbeta*MuSUSY*Sbeta*pow3(Mt) + pow2(Sbeta)*pow3(s2t)*pow4(Mst2)) + 2*
        pow2(Mst1)*pow4(Mst2)*(16*pow2(Mt)*pow2(s2t)*(-(pow2(MuSUSY)*(-1 +
        pow2(Sbeta))) + pow2(Mst2)*pow2(Sbeta)) + 20*Cbeta*MuSUSY*s2t*Sbeta*
        pow3(Mt) + 2*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*pow3(s2t) - 120*pow2(
        Sbeta)*pow4(Mt) - 3*pow2(Sbeta)*pow4(Mst2)*pow4(s2t)) + 160*Mt*(-6*
        Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*pow2(s2t) - 12*s2t*pow2(Mt)*pow2(
        MuSUSY)*(-1 + pow2(Sbeta)) + 6*Cbeta*MuSUSY*Sbeta*pow3(Mt) + pow2(
        Sbeta)*pow3(s2t)*pow4(Mst2))*pow5(Mst1) - 320*Mst1*s2t*pow2(Sbeta)*
        pow3(Mt)*pow6(Mst2) - 4*Mt*Sbeta*(13*Mt*Sbeta*pow2(Mst2)*pow2(s2t) +
        16*Sbeta*pow3(Mt) + 5*Cbeta*MuSUSY*pow2(Mst2)*pow3(s2t))*pow6(Mst2)) -
        1680*pow6(Msq)*(pow2(Mst1)*pow2(s2t)*(8*pow2(Mt)*pow2(MuSUSY)*(-1 +
        pow2(Sbeta)) + pow2(s2t)*pow2(Sbeta)*pow4(Mst2)) - pow2(Sbeta)*pow4(
        s2t)*pow6(Mst2)) + 16*pow2(Mt)*(-37*Cbeta*Mt*MuSUSY*s2t*Sbeta + 9*pow2(
        MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta)) + 9*pow2(Mt)*pow2(Sbeta))*pow8(
        Mst2) + 2240*Mst1*s2t*pow2(Sbeta)*pow3(Mt)*pow8(Mst2)))))/(4410.*pow2(
        Mgl)*pow2(Mst1)*pow2(Sbeta)*pow4(Msq)*pow4(Mst2)*pow4(Mt)) + (pow2(log(
        pow2(Mst1)/pow2(Mst2)))*(pow2(Mgl)*(30*pow2(Mst2)*(32*Mt*s2t*pow2(Mst1)
        *(-70*Cbeta*MuSUSY*Sbeta*pow2(Mt) + 23*Cbeta*MuSUSY*Sbeta*pow2(Mst2)*
        pow2(s2t) + 6*Mt*s2t*(9*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 4*pow2(Mst2)*
        pow2(Sbeta)))*pow4(Mst2) + s2t*pow2(Mst2)*pow4(Mst1)*(420*Cbeta*Mt*
        MuSUSY*Sbeta*pow2(Mst2)*pow2(s2t) + 2568*s2t*pow2(Mt)*pow2(MuSUSY)*(-1
        + pow2(Sbeta)) - 2240*Cbeta*MuSUSY*Sbeta*pow3(Mt) + 79*pow2(Sbeta)*
        pow3(s2t)*pow4(Mst2)) + 1680*Mt*s2t*(3*Cbeta*Mt*MuSUSY*s2t*Sbeta*pow2(
        Mst2) + 4*pow2(Mt)*pow2(MuSUSY)*(-1 + pow2(Sbeta)) - pow2(s2t)*pow2(
        Sbeta)*pow4(Mst2))*pow5(Mst1) - 64*pow2(Mt)*(11*Cbeta*Mt*MuSUSY*s2t*
        Sbeta - 4*pow2(MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta)) + 24*pow2(Mt)*pow2(
        Sbeta))*pow6(Mst2)) + 840*pow2(Msq)*(-160*Mt*s2t*pow2(Mst2)*pow3(Mst1)*
        (-3*Cbeta*Mt*MuSUSY*s2t*Sbeta*pow2(Mst2) - 4*pow2(Mt)*pow2(MuSUSY)*(-1
        + pow2(Sbeta)) + pow2(s2t)*pow2(Sbeta)*pow4(Mst2)) + s2t*pow2(Mst1)*
        pow4(Mst2)*(76*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*pow2(s2t) + 8*s2t*pow2(
        Mt)*(29*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 6*pow2(Mst2)*pow2(Sbeta)) -
        240*Cbeta*MuSUSY*Sbeta*pow3(Mt) - 9*pow2(Sbeta)*pow3(s2t)*pow4(Mst2)) -
        4*s2t*pow2(Mst2)*pow4(Mst1)*(-75*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*pow2(
        s2t) - 208*s2t*pow2(Mt)*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 60*Cbeta*
        MuSUSY*Sbeta*pow3(Mt) + 14*pow2(Sbeta)*pow3(s2t)*pow4(Mst2)) + 160*Mt*
        Sbeta*(3*Cbeta*Mt*MuSUSY*(4*pow2(Mt) - pow2(Mst2)*pow2(s2t)) + 2*s2t*
        Sbeta*pow2(Mst2)*(-6*pow2(Mt) + pow2(Mst2)*pow2(s2t)))*pow5(Mst1) + 8*
        Mt*Sbeta*(9*Mt*Sbeta*pow2(Mst2)*pow2(s2t) - 12*Sbeta*pow3(Mt) + 5*
        Cbeta*MuSUSY*pow2(Mst2)*pow3(s2t))*pow6(Mst2)) + pow4(Msq)*(3360*Mst1*
        Mt*s2t*(-128*pow2(Mt) + 797*pow2(Mst2)*pow2(s2t))*pow2(Sbeta)*pow4(
        Mst2) - 1120*Mt*pow3(Mst1)*(9750*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*pow2(
        s2t) + 4*s2t*pow2(Mt)*(5641*pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 2094*
        pow2(Mst2)*pow2(Sbeta)) + 7032*Cbeta*MuSUSY*Sbeta*pow3(Mt) - 859*pow2(
        Sbeta)*pow3(s2t)*pow4(Mst2)) + 210*pow2(Sbeta)*pow4(Mst2)*(-6112*pow2(
        Mst2)*pow2(Mt)*pow2(s2t) + 32*pow4(Mt) + 2193*pow4(Mst2)*pow4(s2t)) -
        pow4(Mst1)*(8*pow2(Mt)*pow2(s2t)*(2258059*pow2(MuSUSY)*(-1 + pow2(
        Sbeta)) + 184876*pow2(Mst2)*pow2(Sbeta)) + 1626368*Cbeta*MuSUSY*s2t*
        Sbeta*pow3(Mt) + 3082012*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*pow3(s2t) +
        2834624*pow2(Sbeta)*pow4(Mt) + 256523*pow2(Sbeta)*pow4(Mst2)*pow4(s2t))
        + 560*Mt*s2t*Sbeta*(-27237*Cbeta*Mt*MuSUSY*s2t + 37192*Sbeta*pow2(Mt) +
        2579*Sbeta*pow2(Mst2)*pow2(s2t))*pow5(Mst1) - 56*pow2(Mst1)*(-12*pow2(
        Mst2)*pow2(Mt)*pow2(s2t)*(-17709*pow2(MuSUSY)*(-1 + pow2(Sbeta)) +
        4601*pow2(Mst2)*pow2(Sbeta)) + 70344*Cbeta*MuSUSY*s2t*Sbeta*pow2(Mst2)*
        pow3(Mt) + 73359*Cbeta*Mt*MuSUSY*Sbeta*pow3(s2t)*pow4(Mst2) + 32*(1440*
        pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 823*pow2(Mst2)*pow2(Sbeta))*pow4(Mt)
        - 10116*pow2(Sbeta)*pow4(s2t)*pow6(Mst2)))) + 140*Dmglst1*Mgl*(150*s2t*
        pow2(Mst2)*pow4(Mst1)*(12*Mst1*Mt*(3*Cbeta*Mt*MuSUSY*s2t*Sbeta*pow2(
        Mst2) + 4*pow2(Mt)*pow2(MuSUSY)*(-1 + pow2(Sbeta)) - pow2(s2t)*pow2(
        Sbeta)*pow4(Mst2)) + s2t*pow2(Mst2)*(-4*Cbeta*Mt*MuSUSY*s2t*Sbeta*pow2(
        Mst2) - 8*pow2(Mt)*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + pow2(s2t)*pow2(
        Sbeta)*pow4(Mst2))) + pow4(Msq)*(24*Mst1*Mt*s2t*(-128*pow2(Mt) + 797*
        pow2(Mst2)*pow2(s2t))*pow2(Sbeta)*pow4(Mst2) + 54*pow2(Sbeta)*pow2(-4*
        pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow4(Mst2) - 8*Mt*pow3(Mst1)*(-1926*
        Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*pow2(s2t) + 4*s2t*pow2(Mt)*(1749*pow2(
        MuSUSY)*(-1 + pow2(Sbeta)) - 6536*pow2(Mst2)*pow2(Sbeta)) + 24296*
        Cbeta*MuSUSY*Sbeta*pow3(Mt) + 3033*pow2(Sbeta)*pow3(s2t)*pow4(Mst2)) +
        4*pow2(Mst1)*pow2(Mst2)*pow2(Sbeta)*(8200*pow2(Mst2)*pow2(Mt)*pow2(s2t)
        + 5648*pow4(Mt) + 1089*pow4(Mst2)*pow4(s2t)) + pow4(Mst1)*(-8*pow2(Mt)*
        pow2(s2t)*(9491*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 10832*pow2(Mst2)*
        pow2(Sbeta)) + 85536*Cbeta*MuSUSY*s2t*Sbeta*pow3(Mt) - 20108*Cbeta*Mt*
        MuSUSY*Sbeta*pow2(Mst2)*pow3(s2t) - 92736*pow2(Sbeta)*pow4(Mt) + 617*
        pow2(Sbeta)*pow4(Mst2)*pow4(s2t)) + 4*Mt*s2t*Sbeta*(20559*Cbeta*Mt*
        MuSUSY*s2t + 184728*Sbeta*pow2(Mt) - 5569*Sbeta*pow2(Mst2)*pow2(s2t))*
        pow5(Mst1)) + 60*pow2(Msq)*pow2(Mst1)*(-48*Mst1*Mt*s2t*pow2(Mst2)*(-3*
        Cbeta*Mt*MuSUSY*s2t*Sbeta*pow2(Mst2) - 4*pow2(Mt)*pow2(MuSUSY)*(-1 +
        pow2(Sbeta)) + pow2(s2t)*pow2(Sbeta)*pow4(Mst2)) - 4*pow2(Mst1)*pow2(
        Mst2)*pow2(s2t)*(-17*Cbeta*Mt*MuSUSY*s2t*Sbeta*pow2(Mst2) - 28*pow2(Mt)
        *pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 5*pow2(s2t)*pow2(Sbeta)*pow4(Mst2))
        + 48*Mt*pow3(Mst1)*(-9*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*pow2(s2t) - 4*
        s2t*pow2(Mt)*(2*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 5*pow2(Mst2)*pow2(
        Sbeta)) + 20*Cbeta*MuSUSY*Sbeta*pow3(Mt) + 4*pow2(Sbeta)*pow3(s2t)*
        pow4(Mst2)) + 3*pow2(Sbeta)*pow4(s2t)*pow8(Mst2))) + 70*pow2(Dmglst1)*(
        150*s2t*pow2(Mst2)*pow4(Mst1)*(72*Mst1*Mt*(3*Cbeta*Mt*MuSUSY*s2t*Sbeta*
        pow2(Mst2) + 4*pow2(Mt)*pow2(MuSUSY)*(-1 + pow2(Sbeta)) - pow2(s2t)*
        pow2(Sbeta)*pow4(Mst2)) + 5*s2t*pow2(Mst2)*(-4*Cbeta*Mt*MuSUSY*s2t*
        Sbeta*pow2(Mst2) - 8*pow2(Mt)*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + pow2(
        s2t)*pow2(Sbeta)*pow4(Mst2))) + pow4(Msq)*(48*Mst1*Mt*s2t*(-128*pow2(
        Mt) + 797*pow2(Mst2)*pow2(s2t))*pow2(Sbeta)*pow4(Mst2) + 162*pow2(
        Sbeta)*pow2(-4*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow4(Mst2) - 16*Mt*
        pow3(Mst1)*(-1926*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*pow2(s2t) + 4*s2t*
        pow2(Mt)*(1749*pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 13007*pow2(Mst2)*pow2(
        Sbeta)) + 24296*Cbeta*MuSUSY*Sbeta*pow3(Mt) + 10107*pow2(Sbeta)*pow3(
        s2t)*pow4(Mst2)) + 12*pow2(Mst1)*pow2(Mst2)*pow2(Sbeta)*(7624*pow2(
        Mst2)*pow2(Mt)*pow2(s2t) + 5648*pow4(Mt) + 1089*pow4(Mst2)*pow4(s2t)) +
        pow4(Mst1)*(-8*pow2(Mt)*pow2(s2t)*(45823*pow2(MuSUSY)*(-1 + pow2(Sbeta)
        ) + 71536*pow2(Mst2)*pow2(Sbeta)) + 858208*Cbeta*MuSUSY*s2t*Sbeta*pow3(
        Mt) - 129724*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*pow3(s2t) - 474624*pow2(
        Sbeta)*pow4(Mt) + 19201*pow2(Sbeta)*pow4(Mst2)*pow4(s2t)) + 8*Mt*s2t*
        Sbeta*(126513*Cbeta*Mt*MuSUSY*s2t + 533488*Sbeta*pow2(Mt) - 26739*
        Sbeta*pow2(Mst2)*pow2(s2t))*pow5(Mst1)) + 60*pow2(Msq)*pow2(Mst1)*(-96*
        Mst1*Mt*s2t*pow2(Mst2)*(-3*Cbeta*Mt*MuSUSY*s2t*Sbeta*pow2(Mst2) - 4*
        pow2(Mt)*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 2*pow2(s2t)*pow2(Sbeta)*
        pow4(Mst2)) - 4*pow2(Mst1)*pow2(Mst2)*pow2(s2t)*(-91*Cbeta*Mt*MuSUSY*
        s2t*Sbeta*pow2(Mst2) - 164*pow2(Mt)*pow2(MuSUSY)*(-1 + pow2(Sbeta)) +
        25*pow2(s2t)*pow2(Sbeta)*pow4(Mst2)) + 96*Mt*pow3(Mst1)*(-33*Cbeta*Mt*
        MuSUSY*Sbeta*pow2(Mst2)*pow2(s2t) - 12*s2t*pow2(Mt)*(3*pow2(MuSUSY)*(-1
        + pow2(Sbeta)) + 5*pow2(Mst2)*pow2(Sbeta)) + 60*Cbeta*MuSUSY*Sbeta*
        pow3(Mt) + 13*pow2(Sbeta)*pow3(s2t)*pow4(Mst2)) + 9*pow2(Sbeta)*pow4(
        s2t)*pow8(Mst2)))))/(2520.*pow2(Mgl)*pow2(Sbeta)*pow4(Msq)*pow4(Mst2)*
        pow4(Mt));

   return result;
}

/**
 * @return returns the susy log^1 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double H5g1::calc_coef_at_as2_no_sm_logs_log1() const {

   const double result =
      -447.7037037037037 - (11024*Dmglst1)/(27.*Mgl) + (1472*Mst1*s2t)/(3.*Mt)
        - (40192*Dmglst1*Mst1*s2t)/(27.*Mgl*Mt) + (1328*z2)/3. + (896*Mst1*s2t*
        z2)/(3.*Mt) + (896*Dmglst1*Mst1*s2t*z2)/(3.*Mgl*Mt) + 128*z3 - (583376*
        pow2(Dmglst1))/(675.*pow2(Mgl)) - (13376*Mst1*s2t*pow2(Dmglst1))/(3.*
        Mt*pow2(Mgl)) + (896*Mst1*s2t*z2*pow2(Dmglst1))/(3.*Mt*pow2(Mgl)) - (
        320*pow2(Msq))/pow2(Mst1) - (880*pow2(Mst1))/(3.*pow2(Msq)) - (160*
        Dmglst1*pow2(Mst1))/(Mgl*pow2(Msq)) - (112*pow2(Dmglst1)*pow2(Mst1))/(
        3.*pow2(Mgl)*pow2(Msq)) - (320*pow2(Msq))/pow2(Mst2) + (448*pow2(Mst1))
        /(3.*pow2(Mst2)) + (19648*Dmglst1*pow2(Mst1))/(9.*Mgl*pow2(Mst2)) - (
        15008*Cbeta*MuSUSY*s2t*pow2(Mst1))/(3.*Mt*Sbeta*pow2(Mst2)) + (2368*z2*
        pow2(Mst1))/(3.*pow2(Mst2)) + (1536*Dmglst1*z2*pow2(Mst1))/(Mgl*pow2(
        Mst2)) - (1024*Cbeta*MuSUSY*s2t*z2*pow2(Mst1))/(Mt*Sbeta*pow2(Mst2)) +
        (503056*pow2(Dmglst1)*pow2(Mst1))/(135.*pow2(Mgl)*pow2(Mst2)) + (2304*
        z2*pow2(Dmglst1)*pow2(Mst1))/(pow2(Mgl)*pow2(Mst2)) - (1280*pow2(Mst2))
        /(9.*pow2(Msq)) - (3520*Mst1*s2t*pow2(Mst2))/(9.*Mt*pow2(Msq)) - (3520*
        Dmglst1*Mst1*s2t*pow2(Mst2))/(9.*Mgl*Mt*pow2(Msq)) - (3520*Mst1*s2t*
        pow2(Dmglst1)*pow2(Mst2))/(9.*Mt*pow2(Mgl)*pow2(Msq)) - (32*pow2(Mst2))
        /pow2(Mst1) - (320*pow2(Msq)*pow2(s2t))/pow2(Mt) + (8080*pow2(Mst1)*
        pow2(s2t))/(3.*pow2(Mt)) + (52816*Dmglst1*pow2(Mst1)*pow2(s2t))/(9.*
        Mgl*pow2(Mt)) + (512*z2*pow2(Mst1)*pow2(s2t))/pow2(Mt) + (1024*Dmglst1*
        z2*pow2(Mst1)*pow2(s2t))/(Mgl*pow2(Mt)) + (86552*pow2(Dmglst1)*pow2(
        Mst1)*pow2(s2t))/(9.*pow2(Mgl)*pow2(Mt)) + (1536*z2*pow2(Dmglst1)*pow2(
        Mst1)*pow2(s2t))/(pow2(Mgl)*pow2(Mt)) + (160*pow2(Msq)*pow2(Mst1)*pow2(
        s2t))/(pow2(Mst2)*pow2(Mt)) - (208*pow2(Mst2)*pow2(s2t))/pow2(Mt) - (
        592*Dmglst1*pow2(Mst2)*pow2(s2t))/(9.*Mgl*pow2(Mt)) + (1384*pow2(
        Dmglst1)*pow2(Mst2)*pow2(s2t))/(9.*pow2(Mgl)*pow2(Mt)) + (160*pow2(Msq)
        *pow2(Mst2)*pow2(s2t))/(pow2(Mst1)*pow2(Mt)) - (960*pow2(Mst1)*pow2(
        MuSUSY)*pow2(s2t))/(pow2(Mst2)*pow2(Mt)) - (64*z2*pow2(Mst1)*pow2(
        MuSUSY)*pow2(s2t))/(3.*pow2(Mst2)*pow2(Mt)) + (960*pow2(Mst1)*pow2(
        MuSUSY)*pow2(s2t))/(pow2(Mst2)*pow2(Mt)*pow2(Sbeta)) + (64*z2*pow2(
        Mst1)*pow2(MuSUSY)*pow2(s2t))/(3.*pow2(Mst2)*pow2(Mt)*pow2(Sbeta)) - (
        1280*s2t*pow3(Mst1))/(9.*Mt*pow2(Msq)) - (5120*Dmglst1*s2t*pow3(Mst1))/
        (3.*Mgl*Mt*pow2(Msq)) - (14080*s2t*pow2(Dmglst1)*pow3(Mst1))/(3.*Mt*
        pow2(Mgl)*pow2(Msq)) + (4096*s2t*pow3(Mst1))/(3.*Mt*pow2(Mst2)) + (
        32768*Dmglst1*s2t*pow3(Mst1))/(3.*Mgl*Mt*pow2(Mst2)) + (8704*s2t*z2*
        pow3(Mst1))/(3.*Mt*pow2(Mst2)) + (25216*Dmglst1*s2t*z2*pow3(Mst1))/(3.*
        Mgl*Mt*pow2(Mst2)) + (238016*s2t*pow2(Dmglst1)*pow3(Mst1))/(9.*Mt*pow2(
        Mgl)*pow2(Mst2)) + (49984*s2t*z2*pow2(Dmglst1)*pow3(Mst1))/(3.*Mt*pow2(
        Mgl)*pow2(Mst2)) + (1600*Cbeta*MuSUSY*pow3(Mst1))/(3.*Sbeta*pow2(Msq)*
        pow2(Mst2)) + (18880*Cbeta*Dmglst1*MuSUSY*pow3(Mst1))/(9.*Mgl*Sbeta*
        pow2(Msq)*pow2(Mst2)) + (18880*Cbeta*MuSUSY*pow2(Dmglst1)*pow3(Mst1))/(
        9.*Sbeta*pow2(Mgl)*pow2(Msq)*pow2(Mst2)) + (32*Cbeta*MuSUSY*pow2(s2t)*
        pow3(Mst1))/(Sbeta*pow2(Mst2)*pow2(Mt)) + (224*Cbeta*Dmglst1*MuSUSY*
        pow2(s2t)*pow3(Mst1))/(Mgl*Sbeta*pow2(Mst2)*pow2(Mt)) + (64*Cbeta*
        MuSUSY*z2*pow2(s2t)*pow3(Mst1))/(Sbeta*pow2(Mst2)*pow2(Mt)) + (2368*
        Cbeta*Dmglst1*MuSUSY*z2*pow2(s2t)*pow3(Mst1))/(Mgl*Sbeta*pow2(Mst2)*
        pow2(Mt)) + (224*Cbeta*MuSUSY*pow2(Dmglst1)*pow2(s2t)*pow3(Mst1))/(
        Sbeta*pow2(Mgl)*pow2(Mst2)*pow2(Mt)) + (2368*Cbeta*MuSUSY*z2*pow2(
        Dmglst1)*pow2(s2t)*pow3(Mst1))/(Sbeta*pow2(Mgl)*pow2(Mst2)*pow2(Mt)) +
        (416*Cbeta*MuSUSY*pow2(Mst1)*pow3(s2t))/(Sbeta*pow3(Mt)) - (120*Cbeta*
        MuSUSY*z2*pow2(Mst1)*pow3(s2t))/(Sbeta*pow3(Mt)) - (80*Cbeta*MuSUSY*
        pow2(Msq)*pow2(Mst1)*pow3(s2t))/(Sbeta*pow2(Mst2)*pow3(Mt)) + (5344*
        Mst1*pow2(Mst2)*pow3(s2t))/(3.*pow3(Mt)) + (6464*Dmglst1*Mst1*pow2(
        Mst2)*pow3(s2t))/(3.*Mgl*pow3(Mt)) + (544*Mst1*z2*pow2(Mst2)*pow3(s2t))
        /(3.*pow3(Mt)) + (544*Dmglst1*Mst1*z2*pow2(Mst2)*pow3(s2t))/(3.*Mgl*
        pow3(Mt)) + (19664*Mst1*pow2(Dmglst1)*pow2(Mst2)*pow3(s2t))/(9.*pow2(
        Mgl)*pow3(Mt)) + (544*Mst1*z2*pow2(Dmglst1)*pow2(Mst2)*pow3(s2t))/(3.*
        pow2(Mgl)*pow3(Mt)) - (1792*pow3(Mst1)*pow3(s2t))/pow3(Mt) - (6688*
        Dmglst1*pow3(Mst1)*pow3(s2t))/(3.*Mgl*pow3(Mt)) - (608*z2*pow3(Mst1)*
        pow3(s2t))/(3.*pow3(Mt)) - (2912*Dmglst1*z2*pow3(Mst1)*pow3(s2t))/(3.*
        Mgl*pow3(Mt)) - (34832*pow2(Dmglst1)*pow3(Mst1)*pow3(s2t))/(9.*pow2(
        Mgl)*pow3(Mt)) - (6368*z2*pow2(Dmglst1)*pow3(Mst1)*pow3(s2t))/(3.*pow2(
        Mgl)*pow3(Mt)) - (140*pow2(Mst1)*pow2(Mst2))/(3.*pow4(Msq)) - (280*
        Dmglst1*pow2(Mst1)*pow2(Mst2))/(3.*Mgl*pow4(Msq)) - (140*pow2(Dmglst1)*
        pow2(Mst1)*pow2(Mst2))/(pow2(Mgl)*pow4(Msq)) + (520*Cbeta*MuSUSY*pow3(
        Mst1))/(3.*Sbeta*pow4(Msq)) + (3640*Cbeta*Dmglst1*MuSUSY*pow3(Mst1))/(
        9.*Mgl*Sbeta*pow4(Msq)) + (3640*Cbeta*MuSUSY*pow2(Dmglst1)*pow3(Mst1))/
        (9.*Sbeta*pow2(Mgl)*pow4(Msq)) - (1040*s2t*pow2(Mst2)*pow3(Mst1))/(9.*
        Mt*pow4(Msq)) - (1040*Dmglst1*s2t*pow2(Mst2)*pow3(Mst1))/(3.*Mgl*Mt*
        pow4(Msq)) - (2080*s2t*pow2(Dmglst1)*pow2(Mst2)*pow3(Mst1))/(3.*Mt*
        pow2(Mgl)*pow4(Msq)) + (800*pow4(Mst1))/(3.*pow2(Msq)*pow2(Mst2)) + (
        3200*Dmglst1*pow4(Mst1))/(3.*Mgl*pow2(Msq)*pow2(Mst2)) + (8000*pow2(
        Dmglst1)*pow4(Mst1))/(3.*pow2(Mgl)*pow2(Msq)*pow2(Mst2)) - (1168*pow2(
        s2t)*pow4(Mst1))/(pow2(Mst2)*pow2(Mt)) - (2208*Dmglst1*pow2(s2t)*pow4(
        Mst1))/(Mgl*pow2(Mst2)*pow2(Mt)) - (2048*Dmglst1*z2*pow2(s2t)*pow4(
        Mst1))/(Mgl*pow2(Mst2)*pow2(Mt)) - (18800*pow2(Dmglst1)*pow2(s2t)*pow4(
        Mst1))/(3.*pow2(Mgl)*pow2(Mst2)*pow2(Mt)) - (7168*z2*pow2(Dmglst1)*
        pow2(s2t)*pow4(Mst1))/(pow2(Mgl)*pow2(Mst2)*pow2(Mt)) + (506*Cbeta*
        MuSUSY*pow3(s2t)*pow4(Mst1))/(Sbeta*pow2(Mst2)*pow3(Mt)) + (1072*Cbeta*
        Dmglst1*MuSUSY*pow3(s2t)*pow4(Mst1))/(Mgl*Sbeta*pow2(Mst2)*pow3(Mt)) +
        (1080*Cbeta*MuSUSY*pow2(Dmglst1)*pow3(s2t)*pow4(Mst1))/(Sbeta*pow2(Mgl)
        *pow2(Mst2)*pow3(Mt)) - (410*pow4(Mst1))/(9.*pow4(Msq)) + (760*Dmglst1*
        pow4(Mst1))/(9.*Mgl*pow4(Msq)) + (3916*pow2(Dmglst1)*pow4(Mst1))/(9.*
        pow2(Mgl)*pow4(Msq)) - (2560*pow2(Mst1)*pow2(MuSUSY))/pow4(Mst2) - (
        2048*z2*pow2(Mst1)*pow2(MuSUSY))/(3.*pow4(Mst2)) + (2560*pow2(Mst1)*
        pow2(MuSUSY))/(pow2(Sbeta)*pow4(Mst2)) + (2048*z2*pow2(Mst1)*pow2(
        MuSUSY))/(3.*pow2(Sbeta)*pow4(Mst2)) - (1632*Cbeta*MuSUSY*pow3(Mst1))/(
        Sbeta*pow4(Mst2)) - (230624*Cbeta*Dmglst1*MuSUSY*pow3(Mst1))/(27.*Mgl*
        Sbeta*pow4(Mst2)) - (3200*Cbeta*MuSUSY*z2*pow3(Mst1))/(Sbeta*pow4(Mst2)
        ) - (8704*Cbeta*Dmglst1*MuSUSY*z2*pow3(Mst1))/(Mgl*Sbeta*pow4(Mst2)) -
        (230624*Cbeta*MuSUSY*pow2(Dmglst1)*pow3(Mst1))/(27.*Sbeta*pow2(Mgl)*
        pow4(Mst2)) - (8704*Cbeta*MuSUSY*z2*pow2(Dmglst1)*pow3(Mst1))/(Sbeta*
        pow2(Mgl)*pow4(Mst2)) - (21248*s2t*pow2(MuSUSY)*pow3(Mst1))/(3.*Mt*
        pow4(Mst2)) - (8320*Dmglst1*s2t*pow2(MuSUSY)*pow3(Mst1))/(Mgl*Mt*pow4(
        Mst2)) - (640*s2t*z2*pow2(MuSUSY)*pow3(Mst1))/(Mt*pow4(Mst2)) + (2432*
        Dmglst1*s2t*z2*pow2(MuSUSY)*pow3(Mst1))/(Mgl*Mt*pow4(Mst2)) - (8320*
        s2t*pow2(Dmglst1)*pow2(MuSUSY)*pow3(Mst1))/(Mt*pow2(Mgl)*pow4(Mst2)) +
        (2432*s2t*z2*pow2(Dmglst1)*pow2(MuSUSY)*pow3(Mst1))/(Mt*pow2(Mgl)*pow4(
        Mst2)) + (21248*s2t*pow2(MuSUSY)*pow3(Mst1))/(3.*Mt*pow2(Sbeta)*pow4(
        Mst2)) + (8320*Dmglst1*s2t*pow2(MuSUSY)*pow3(Mst1))/(Mgl*Mt*pow2(Sbeta)
        *pow4(Mst2)) + (640*s2t*z2*pow2(MuSUSY)*pow3(Mst1))/(Mt*pow2(Sbeta)*
        pow4(Mst2)) - (2432*Dmglst1*s2t*z2*pow2(MuSUSY)*pow3(Mst1))/(Mgl*Mt*
        pow2(Sbeta)*pow4(Mst2)) + (8320*s2t*pow2(Dmglst1)*pow2(MuSUSY)*pow3(
        Mst1))/(Mt*pow2(Mgl)*pow2(Sbeta)*pow4(Mst2)) - (2432*s2t*z2*pow2(
        Dmglst1)*pow2(MuSUSY)*pow3(Mst1))/(Mt*pow2(Mgl)*pow2(Sbeta)*pow4(Mst2))
        - (2756*pow4(Mst1))/pow4(Mst2) - (76928*Dmglst1*pow4(Mst1))/(9.*Mgl*
        pow4(Mst2)) - (8000*Cbeta*MuSUSY*s2t*pow4(Mst1))/(3.*Mt*Sbeta*pow4(
        Mst2)) - (21568*Cbeta*Dmglst1*MuSUSY*s2t*pow4(Mst1))/(3.*Mgl*Mt*Sbeta*
        pow4(Mst2)) + (256*z2*pow4(Mst1))/pow4(Mst2) - (256*Dmglst1*z2*pow4(
        Mst1))/(Mgl*pow4(Mst2)) - (1024*Cbeta*MuSUSY*s2t*z2*pow4(Mst1))/(Mt*
        Sbeta*pow4(Mst2)) + (2048*Cbeta*Dmglst1*MuSUSY*s2t*z2*pow4(Mst1))/(Mgl*
        Mt*Sbeta*pow4(Mst2)) - (2085296*pow2(Dmglst1)*pow4(Mst1))/(135.*pow2(
        Mgl)*pow4(Mst2)) - (7008*Cbeta*MuSUSY*s2t*pow2(Dmglst1)*pow4(Mst1))/(
        Mt*Sbeta*pow2(Mgl)*pow4(Mst2)) - (1856*z2*pow2(Dmglst1)*pow4(Mst1))/(
        pow2(Mgl)*pow4(Mst2)) + (11264*Cbeta*MuSUSY*s2t*z2*pow2(Dmglst1)*pow4(
        Mst1))/(Mt*Sbeta*pow2(Mgl)*pow4(Mst2)) + (52*pow2(MuSUSY)*pow2(s2t)*
        pow4(Mst1))/(pow2(Mt)*pow4(Mst2)) + (1024*Dmglst1*pow2(MuSUSY)*pow2(
        s2t)*pow4(Mst1))/(Mgl*pow2(Mt)*pow4(Mst2)) - (64*z2*pow2(MuSUSY)*pow2(
        s2t)*pow4(Mst1))/(3.*pow2(Mt)*pow4(Mst2)) - (2752*Dmglst1*z2*pow2(
        MuSUSY)*pow2(s2t)*pow4(Mst1))/(3.*Mgl*pow2(Mt)*pow4(Mst2)) - (352*pow2(
        Dmglst1)*pow2(MuSUSY)*pow2(s2t)*pow4(Mst1))/(3.*pow2(Mgl)*pow2(Mt)*
        pow4(Mst2)) - (1376*z2*pow2(Dmglst1)*pow2(MuSUSY)*pow2(s2t)*pow4(Mst1))
        /(pow2(Mgl)*pow2(Mt)*pow4(Mst2)) - (52*pow2(MuSUSY)*pow2(s2t)*pow4(
        Mst1))/(pow2(Mt)*pow2(Sbeta)*pow4(Mst2)) - (1024*Dmglst1*pow2(MuSUSY)*
        pow2(s2t)*pow4(Mst1))/(Mgl*pow2(Mt)*pow2(Sbeta)*pow4(Mst2)) + (64*z2*
        pow2(MuSUSY)*pow2(s2t)*pow4(Mst1))/(3.*pow2(Mt)*pow2(Sbeta)*pow4(Mst2))
        + (2752*Dmglst1*z2*pow2(MuSUSY)*pow2(s2t)*pow4(Mst1))/(3.*Mgl*pow2(Mt)*
        pow2(Sbeta)*pow4(Mst2)) + (352*pow2(Dmglst1)*pow2(MuSUSY)*pow2(s2t)*
        pow4(Mst1))/(3.*pow2(Mgl)*pow2(Mt)*pow2(Sbeta)*pow4(Mst2)) + (1376*z2*
        pow2(Dmglst1)*pow2(MuSUSY)*pow2(s2t)*pow4(Mst1))/(pow2(Mgl)*pow2(Mt)*
        pow2(Sbeta)*pow4(Mst2)) + (16*pow2(s2t)*pow4(Mst2))/(pow2(Mst1)*pow2(
        Mt)) + (8*Cbeta*MuSUSY*pow3(s2t)*pow4(Mst2))/(Sbeta*pow2(Mst1)*pow3(Mt)
        ) - (30*pow4(Mst2))/pow4(Msq) - (520*Mst1*s2t*pow4(Mst2))/(9.*Mt*pow4(
        Msq)) - (520*Dmglst1*Mst1*s2t*pow4(Mst2))/(9.*Mgl*Mt*pow4(Msq)) - (520*
        Mst1*s2t*pow2(Dmglst1)*pow4(Mst2))/(9.*Mt*pow2(Mgl)*pow4(Msq)) + (20*
        pow2(Msq)*pow2(Mst1)*pow4(s2t))/pow4(Mt) + (20*pow2(Msq)*pow2(Mst2)*
        pow4(s2t))/pow4(Mt) - (330*pow2(Mst1)*pow2(Mst2)*pow4(s2t))/pow4(Mt) +
        (116*Dmglst1*pow2(Mst1)*pow2(Mst2)*pow4(s2t))/(Mgl*pow4(Mt)) + (172*z2*
        pow2(Mst1)*pow2(Mst2)*pow4(s2t))/(3.*pow4(Mt)) + (344*Dmglst1*z2*pow2(
        Mst1)*pow2(Mst2)*pow4(s2t))/(3.*Mgl*pow4(Mt)) + (2378*pow2(Dmglst1)*
        pow2(Mst1)*pow2(Mst2)*pow4(s2t))/(9.*pow2(Mgl)*pow4(Mt)) + (172*z2*
        pow2(Dmglst1)*pow2(Mst1)*pow2(Mst2)*pow4(s2t))/(pow2(Mgl)*pow4(Mt)) - (
        45*pow4(Mst1)*pow4(s2t))/(2.*pow4(Mt)) - (396*Dmglst1*pow4(Mst1)*pow4(
        s2t))/(Mgl*pow4(Mt)) - (30*z2*pow4(Mst1)*pow4(s2t))/pow4(Mt) - (344*
        Dmglst1*z2*pow4(Mst1)*pow4(s2t))/(3.*Mgl*pow4(Mt)) - (4900*pow2(
        Dmglst1)*pow4(Mst1)*pow4(s2t))/(9.*pow2(Mgl)*pow4(Mt)) - (172*z2*pow2(
        Dmglst1)*pow4(Mst1)*pow4(s2t))/(pow2(Mgl)*pow4(Mt)) - (20*pow2(Msq)*
        pow4(Mst1)*pow4(s2t))/(pow2(Mst2)*pow4(Mt)) + (228*pow4(Mst2)*pow4(s2t)
        )/pow4(Mt) + (12*Dmglst1*pow4(Mst2)*pow4(s2t))/(Mgl*pow4(Mt)) - (82*z2*
        pow4(Mst2)*pow4(s2t))/(3.*pow4(Mt)) + (92*pow2(Dmglst1)*pow4(Mst2)*
        pow4(s2t))/(9.*pow2(Mgl)*pow4(Mt)) - (20*pow2(Msq)*pow4(Mst2)*pow4(s2t)
        )/(pow2(Mst1)*pow4(Mt)) - (640*s2t*pow5(Mst1))/(3.*Mt*pow2(Msq)*pow2(
        Mst2)) - (7040*Dmglst1*s2t*pow5(Mst1))/(3.*Mgl*Mt*pow2(Msq)*pow2(Mst2))
        - (10240*s2t*pow2(Dmglst1)*pow5(Mst1))/(Mt*pow2(Mgl)*pow2(Msq)*pow2(
        Mst2)) - (40*pow3(s2t)*pow5(Mst1))/(3.*pow2(Mst2)*pow3(Mt)) + (2264*
        Dmglst1*pow3(s2t)*pow5(Mst1))/(3.*Mgl*pow2(Mst2)*pow3(Mt)) + (6488*
        pow2(Dmglst1)*pow3(s2t)*pow5(Mst1))/(3.*pow2(Mgl)*pow2(Mst2)*pow3(Mt))
        - (40*s2t*pow5(Mst1))/(3.*Mt*pow4(Msq)) - (2120*Dmglst1*s2t*pow5(Mst1))
        /(3.*Mgl*Mt*pow4(Msq)) - (3080*s2t*pow2(Dmglst1)*pow5(Mst1))/(Mt*pow2(
        Mgl)*pow4(Msq)) + (560*Cbeta*MuSUSY*pow5(Mst1))/(3.*Sbeta*pow2(Mst2)*
        pow4(Msq)) + (10000*Cbeta*Dmglst1*MuSUSY*pow5(Mst1))/(9.*Mgl*Sbeta*
        pow2(Mst2)*pow4(Msq)) + (34480*Cbeta*MuSUSY*pow2(Dmglst1)*pow5(Mst1))/(
        9.*Sbeta*pow2(Mgl)*pow2(Mst2)*pow4(Msq)) - (5984*s2t*pow5(Mst1))/(3.*
        Mt*pow4(Mst2)) + (6496*Dmglst1*s2t*pow5(Mst1))/(Mgl*Mt*pow4(Mst2)) + (
        8960*s2t*z2*pow5(Mst1))/(3.*Mt*pow4(Mst2)) + (43904*Dmglst1*s2t*z2*
        pow5(Mst1))/(3.*Mgl*Mt*pow4(Mst2)) + (463328*s2t*pow2(Dmglst1)*pow5(
        Mst1))/(9.*Mt*pow2(Mgl)*pow4(Mst2)) + (43200*s2t*z2*pow2(Dmglst1)*pow5(
        Mst1))/(Mt*pow2(Mgl)*pow4(Mst2)) + (2240*Cbeta*MuSUSY*pow5(Mst1))/(3.*
        Sbeta*pow2(Msq)*pow4(Mst2)) + (40000*Cbeta*Dmglst1*MuSUSY*pow5(Mst1))/(
        9.*Mgl*Sbeta*pow2(Msq)*pow4(Mst2)) + (137920*Cbeta*MuSUSY*pow2(Dmglst1)
        *pow5(Mst1))/(9.*Sbeta*pow2(Mgl)*pow2(Msq)*pow4(Mst2)) + (72*Cbeta*
        MuSUSY*pow2(s2t)*pow5(Mst1))/(Sbeta*pow2(Mt)*pow4(Mst2)) - (2040*Cbeta*
        Dmglst1*MuSUSY*pow2(s2t)*pow5(Mst1))/(Mgl*Sbeta*pow2(Mt)*pow4(Mst2)) +
        (64*Cbeta*MuSUSY*z2*pow2(s2t)*pow5(Mst1))/(Sbeta*pow2(Mt)*pow4(Mst2)) +
        (2368*Cbeta*Dmglst1*MuSUSY*z2*pow2(s2t)*pow5(Mst1))/(Mgl*Sbeta*pow2(Mt)
        *pow4(Mst2)) - (1432*Cbeta*MuSUSY*pow2(Dmglst1)*pow2(s2t)*pow5(Mst1))/(
        Sbeta*pow2(Mgl)*pow2(Mt)*pow4(Mst2)) + (5824*Cbeta*MuSUSY*z2*pow2(
        Dmglst1)*pow2(s2t)*pow5(Mst1))/(Sbeta*pow2(Mgl)*pow2(Mt)*pow4(Mst2)) +
        (8*log(pow2(Mst1)/pow2(Msq))*(20*Cbeta*Mt*MuSUSY*pow3(Mst1)*(Dmglst1*
        Mgl*(118*pow4(Msq) - 7*pow4(Mst2)) + pow2(Dmglst1)*(118*pow4(Msq) - 7*
        pow4(Mst2)) + 3*pow2(Mgl)*(10*pow4(Msq) - pow4(Mst2))) + Sbeta*(pow2(
        Dmglst1)*(-40*Mst1*pow2(Msq)*(-9*Mst1*Mt + 2*s2t*pow2(Mst1) - 2*s2t*
        pow2(Mst2))*pow4(Mst2) + 10*Mst1*pow4(Mst2)*(9*Mst1*Mt*pow2(Mst2) + 24*
        s2t*pow2(Mst1)*pow2(Mst2) + 59*Mt*pow3(Mst1) - 26*s2t*pow4(Mst1) + 2*
        s2t*pow4(Mst2)) + pow4(Msq)*(540*Mt*pow2(Mst1)*pow2(Mst2) - 4800*s2t*
        pow2(Mst2)*pow3(Mst1) + 6900*Mt*pow4(Mst1) + 68*Mt*pow4(Mst2) - 920*
        Mst1*s2t*pow4(Mst2) - 11520*s2t*pow5(Mst1))) - 5*pow2(Mgl)*(-8*pow2(
        Msq)*(4*Mt*pow2(Mst1) + Mt*pow2(Mst2) + 2*Mst1*s2t*pow2(Mst2) - 2*s2t*
        pow3(Mst1))*pow4(Mst2) - pow4(Mst2)*(6*Mt*pow2(Mst1)*pow2(Mst2) + 8*
        s2t*pow2(Mst2)*pow3(Mst1) + 19*Mt*pow4(Mst1) + 3*Mt*pow4(Mst2) + 4*
        Mst1*s2t*pow4(Mst2) - 12*s2t*pow5(Mst1)) + 6*pow4(Msq)*(-6*Mt*pow2(
        Mst1)*pow2(Mst2) + 8*s2t*pow2(Mst2)*pow3(Mst1) - 14*Mt*pow4(Mst1) - 7*
        Mt*pow4(Mst2) + 12*Mst1*s2t*pow4(Mst2) + 8*s2t*pow5(Mst1))) - 20*
        Dmglst1*Mgl*(4*Mst1*pow2(Msq)*(-3*Mst1*Mt + s2t*pow2(Mst1) - s2t*pow2(
        Mst2))*pow4(Mst2) + Mst1*pow4(Mst2)*(-3*Mst1*Mt*pow2(Mst2) - 6*s2t*
        pow2(Mst1)*pow2(Mst2) - 13*Mt*pow3(Mst1) + 7*s2t*pow4(Mst1) - s2t*pow4(
        Mst2)) + 2*pow4(Msq)*(-9*Mt*pow2(Mst1)*pow2(Mst2) + 42*s2t*pow2(Mst2)*
        pow3(Mst1) - 57*Mt*pow4(Mst1) - 2*Mt*pow4(Mst2) + 17*Mst1*s2t*pow4(
        Mst2) + 66*s2t*pow5(Mst1)))) + 60*log(pow2(Mst1)/pow2(Mst2))*pow4(Msq)*
        (4*Cbeta*Mt*MuSUSY*(7*Dmglst1*Mgl + 7*pow2(Dmglst1) + 3*pow2(Mgl))*
        pow3(Mst1) - Sbeta*(4*Dmglst1*Mgl*Mst1*(6*s2t*pow2(Mst1)*pow2(Mst2) -
        5*Mt*pow3(Mst1) + 15*s2t*pow4(Mst1) + s2t*pow4(Mst2)) + pow2(Mgl)*(8*
        s2t*pow2(Mst2)*pow3(Mst1) - 5*Mt*pow4(Mst1) + Mt*pow4(Mst2) + 4*Mst1*
        s2t*pow4(Mst2) + 12*s2t*pow5(Mst1)) + 2*pow2(Dmglst1)*(24*s2t*pow2(
        Mst2)*pow3(Mst1) - 25*Mt*pow4(Mst1) + 2*Mst1*s2t*pow4(Mst2) + 90*s2t*
        pow5(Mst1))))))/(3.*Mt*Sbeta*pow2(Mgl)*pow4(Msq)*pow4(Mst2)) + (2*pow2(
        log(pow2(Mst1)/pow2(Mst2)))*(2*Dmglst1*Mgl*Mst1*(8*Mt*s2t*Sbeta*(-153*
        Cbeta*Mt*MuSUSY*s2t + 968*Sbeta*pow2(Mt) + 32*Sbeta*pow2(Mst2)*pow2(
        s2t))*pow4(Mst1) + 12*Mt*s2t*(-8*pow2(Mt) + 19*pow2(Mst2)*pow2(s2t))*
        pow2(Sbeta)*pow4(Mst2) - 4*Mt*pow2(Mst1)*(114*Cbeta*Mt*MuSUSY*Sbeta*
        pow2(Mst2)*pow2(s2t) + 4*s2t*pow2(Mt)*(95*pow2(MuSUSY)*(-1 + pow2(
        Sbeta)) - 231*pow2(Mst2)*pow2(Sbeta)) + 900*Cbeta*MuSUSY*Sbeta*pow3(Mt)
        + 19*pow2(Sbeta)*pow3(s2t)*pow4(Mst2)) + Mst1*pow2(Mst2)*pow2(Sbeta)*(
        64*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 528*pow4(Mt) + 155*pow4(Mst2)*pow4(
        s2t)) - pow3(Mst1)*(8*pow2(Mt)*pow2(s2t)*(283*pow2(MuSUSY)*(-1 + pow2(
        Sbeta)) + 192*pow2(Mst2)*pow2(Sbeta)) - 2944*Cbeta*MuSUSY*s2t*Sbeta*
        pow3(Mt) + 512*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*pow3(s2t) + 384*pow2(
        Sbeta)*pow4(Mt) + 27*pow2(Sbeta)*pow4(Mst2)*pow4(s2t))) + Mst1*pow2(
        Dmglst1)*(16*Mt*s2t*Sbeta*(9*Cbeta*Mt*MuSUSY*s2t + 2855*Sbeta*pow2(Mt)
        + 32*Sbeta*pow2(Mst2)*pow2(s2t))*pow4(Mst1) + 24*Mt*s2t*(-8*pow2(Mt) +
        19*pow2(Mst2)*pow2(s2t))*pow2(Sbeta)*pow4(Mst2) - 8*Mt*pow2(Mst1)*(114*
        Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*pow2(s2t) + 2*s2t*pow2(Mt)*(190*pow2(
        MuSUSY)*(-1 + pow2(Sbeta)) - 933*pow2(Mst2)*pow2(Sbeta)) + 900*Cbeta*
        MuSUSY*Sbeta*pow3(Mt) + 127*pow2(Sbeta)*pow3(s2t)*pow4(Mst2)) + 3*Mst1*
        pow2(Mst2)*pow2(Sbeta)*(64*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 528*pow4(Mt)
        + 155*pow4(Mst2)*pow4(s2t)) - 3*pow3(Mst1)*(8*pow2(Mt)*pow2(s2t)*(283*
        pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 384*pow2(Mst2)*pow2(Sbeta)) - 6016*
        Cbeta*MuSUSY*s2t*Sbeta*pow3(Mt) + 512*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*
        pow3(s2t) + 944*pow2(Sbeta)*pow4(Mt) + 27*pow2(Sbeta)*pow4(Mst2)*pow4(
        s2t))) + pow2(Mgl)*(24*Mst1*Mt*s2t*(-8*pow2(Mt) + 19*pow2(Mst2)*pow2(
        s2t))*pow2(Sbeta)*pow4(Mst2) - 8*Mt*pow3(Mst1)*(330*Cbeta*Mt*MuSUSY*
        Sbeta*pow2(Mst2)*pow2(s2t) + 4*s2t*pow2(Mt)*(167*pow2(MuSUSY)*(-1 +
        pow2(Sbeta)) - 74*pow2(Mst2)*pow2(Sbeta)) + 272*Cbeta*MuSUSY*Sbeta*
        pow3(Mt) - 53*pow2(Sbeta)*pow3(s2t)*pow4(Mst2)) + 2*pow2(Sbeta)*pow4(
        Mst2)*(-246*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 106*pow4(Mt) + 41*pow4(
        Mst2)*pow4(s2t)) + pow4(Mst1)*(-24*pow2(Mt)*pow2(s2t)*(231*pow2(MuSUSY)
        *(-1 + pow2(Sbeta)) + 16*pow2(Mst2)*pow2(Sbeta)) + 1296*Cbeta*MuSUSY*
        s2t*Sbeta*pow3(Mt) - 1168*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*pow3(s2t) +
        48*pow2(Sbeta)*pow4(Mt) - 27*pow2(Sbeta)*pow4(Mst2)*pow4(s2t)) + 16*Mt*
        s2t*Sbeta*(-261*Cbeta*Mt*MuSUSY*s2t + 190*Sbeta*pow2(Mt) + 32*Sbeta*
        pow2(Mst2)*pow2(s2t))*pow5(Mst1) + pow2(Mst1)*(4*pow2(Mst2)*pow2(Mt)*
        pow2(s2t)*(-802*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 57*pow2(Mst2)*pow2(
        Sbeta)) + 528*Cbeta*MuSUSY*s2t*Sbeta*pow2(Mst2)*pow3(Mt) - 1276*Cbeta*
        Mt*MuSUSY*Sbeta*pow3(s2t)*pow4(Mst2) - 32*(8*pow2(MuSUSY)*(-1 + pow2(
        Sbeta)) - 17*pow2(Mst2)*pow2(Sbeta))*pow4(Mt) + 237*pow2(Sbeta)*pow4(
        s2t)*pow6(Mst2)))))/(3.*pow2(Mgl)*pow2(Sbeta)*pow4(Mst2)*pow4(Mt)) + (
        log(pow2(Mst1)/pow2(Mst2))*(2*pow2(Dmglst1)*(-600*Mst1*Sbeta*pow2(Mst2)
        *pow3(Mt)*(-2*Cbeta*Mt*MuSUSY*pow2(Mst1)*(58*pow2(Mst1) + 7*pow2(Mst2))
        + Sbeta*pow2(Mst2)*(9*Mst1*Mt*pow2(Mst2) + 24*s2t*pow2(Mst1)*pow2(Mst2)
        + 25*Mt*pow3(Mst1) + 90*s2t*pow4(Mst1) + 2*s2t*pow4(Mst2))) + 1200*
        Mst1*Sbeta*pow2(Msq)*pow3(Mt)*(4*Cbeta*Mt*MuSUSY*pow2(Mst1)*(58*pow2(
        Mst1) + 7*pow2(Mst2)) - Sbeta*pow2(Mst2)*(9*Mst1*Mt*pow2(Mst2) + 48*
        s2t*pow2(Mst1)*pow2(Mst2) + 180*s2t*pow4(Mst1) + 4*s2t*pow4(Mst2))) +
        pow4(Msq)*(80*Mst1*Mt*s2t*(215*pow2(Mt) + 396*pow2(Mst2)*pow2(s2t))*
        pow2(Sbeta)*pow4(Mst2) - 80*Mt*pow3(Mst1)*(1710*Cbeta*Mt*MuSUSY*Sbeta*
        pow2(Mst2)*pow2(s2t) + 3*s2t*pow2(Mt)*(1288*pow2(MuSUSY)*(-1 + pow2(
        Sbeta)) - 2959*pow2(Mst2)*pow2(Sbeta)) + 4438*Cbeta*MuSUSY*Sbeta*pow3(
        Mt) + 414*pow2(Sbeta)*pow3(s2t)*pow4(Mst2)) + 8*pow2(Sbeta)*pow4(Mst2)*
        (-720*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 211*pow4(Mt) + 180*pow4(Mst2)*
        pow4(s2t)) + 5*pow2(Mst1)*pow2(Mst2)*pow2(Sbeta)*(20568*pow2(Mst2)*
        pow2(Mt)*pow2(s2t) + 20216*pow4(Mt) + 2349*pow4(Mst2)*pow4(s2t)) - 5*
        pow4(Mst1)*(72*pow2(Mt)*pow2(s2t)*(491*pow2(MuSUSY)*(-1 + pow2(Sbeta))
        + 125*pow2(Mst2)*pow2(Sbeta)) + 20832*Cbeta*MuSUSY*s2t*Sbeta*pow3(Mt) +
        5976*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*pow3(s2t) + 23624*pow2(Sbeta)*
        pow4(Mt) + 1143*pow2(Sbeta)*pow4(Mst2)*pow4(s2t)) + 40*Mt*s2t*Sbeta*(
        5859*Cbeta*Mt*MuSUSY*s2t + 75062*Sbeta*pow2(Mt) - 1917*Sbeta*pow2(Mst2)
        *pow2(s2t))*pow5(Mst1))) + 20*Dmglst1*Mgl*(-120*Mst1*Sbeta*pow2(Mst2)*
        pow3(Mt)*(-(Cbeta*Mt*MuSUSY*pow2(Mst1)*(22*pow2(Mst1) + 7*pow2(Mst2)))
        + Sbeta*pow2(Mst2)*(3*Mst1*Mt*pow2(Mst2) + 6*s2t*pow2(Mst1)*pow2(Mst2)
        + 5*Mt*pow3(Mst1) + 15*s2t*pow4(Mst1) + s2t*pow4(Mst2))) + 240*Mst1*
        Sbeta*pow2(Msq)*pow3(Mt)*(2*Cbeta*Mt*MuSUSY*pow2(Mst1)*(22*pow2(Mst1) +
        7*pow2(Mst2)) - Sbeta*pow2(Mst2)*(3*Mst1*Mt*pow2(Mst2) + 12*s2t*pow2(
        Mst1)*pow2(Mst2) + 30*s2t*pow4(Mst1) + 2*s2t*pow4(Mst2))) + pow4(Msq)*(
        16*Mst1*Mt*s2t*(103*pow2(Mt) + 198*pow2(Mst2)*pow2(s2t))*pow2(Sbeta)*
        pow4(Mst2) - 16*Mt*pow3(Mst1)*(855*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*
        pow2(s2t) + s2t*pow2(Mt)*(1932*pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 2209*
        pow2(Mst2)*pow2(Sbeta)) + 2219*Cbeta*MuSUSY*Sbeta*pow3(Mt) - 87*pow2(
        Sbeta)*pow3(s2t)*pow4(Mst2)) + 16*pow2(Sbeta)*pow4(Mst2)*(-24*pow2(
        Mst2)*pow2(Mt)*pow2(s2t) + 11*pow4(Mt) + 6*pow4(Mst2)*pow4(s2t)) + 3*
        pow2(Mst1)*pow2(Mst2)*pow2(Sbeta)*(2456*pow2(Mst2)*pow2(Mt)*pow2(s2t) +
        2232*pow4(Mt) + 261*pow4(Mst2)*pow4(s2t)) - pow4(Mst1)*(24*pow2(Mt)*
        pow2(s2t)*(401*pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 263*pow2(Mst2)*pow2(
        Sbeta)) + 26592*Cbeta*MuSUSY*s2t*Sbeta*pow3(Mt) + 912*Cbeta*Mt*MuSUSY*
        Sbeta*pow2(Mst2)*pow3(s2t) + 3200*pow2(Sbeta)*pow4(Mt) + 651*pow2(
        Sbeta)*pow4(Mst2)*pow4(s2t)) + 4*Mt*s2t*Sbeta*(-1827*Cbeta*Mt*MuSUSY*
        s2t + 25312*Sbeta*pow2(Mt) - 531*Sbeta*pow2(Mst2)*pow2(s2t))*pow5(Mst1)
        )) + 15*pow2(Mgl)*(-160*Sbeta*pow2(Msq)*pow3(Mt)*(-12*Cbeta*Mt*MuSUSY*(
        2*pow2(Mst1) + pow2(Mst2))*pow3(Mst1) + Sbeta*pow2(Mst2)*(3*Mt*pow2(
        Mst1)*pow2(Mst2) + 8*s2t*pow2(Mst2)*pow3(Mst1) + 2*Mt*pow4(Mst2) + 4*
        Mst1*s2t*pow4(Mst2) + 12*s2t*pow5(Mst1))) - 40*Sbeta*pow2(Mst2)*pow3(
        Mt)*(-12*Cbeta*Mt*MuSUSY*(2*pow2(Mst1) + pow2(Mst2))*pow3(Mst1) +
        Sbeta*pow2(Mst2)*(6*Mt*pow2(Mst1)*pow2(Mst2) + 8*s2t*pow2(Mst2)*pow3(
        Mst1) + 5*Mt*pow4(Mst1) + 3*Mt*pow4(Mst2) + 4*Mst1*s2t*pow4(Mst2) + 12*
        s2t*pow5(Mst1))) + 120*pow2(s2t)*pow6(Msq)*(pow2(Mst1)*(8*pow2(Mt)*
        pow2(MuSUSY)*(-1 + pow2(Sbeta)) + pow2(s2t)*pow2(Sbeta)*pow4(Mst2)) -
        pow2(s2t)*pow2(Sbeta)*pow6(Mst2)) + pow4(Msq)*(128*Mst1*Mt*s2t*(17*
        pow2(Mt) + 33*pow2(Mst2)*pow2(s2t))*pow2(Sbeta)*pow4(Mst2) - 64*Mt*
        pow3(Mst1)*(345*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*pow2(s2t) + 4*s2t*
        pow2(Mt)*(181*pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 61*pow2(Mst2)*pow2(
        Sbeta)) + 267*Cbeta*MuSUSY*Sbeta*pow3(Mt) - 49*pow2(Sbeta)*pow3(s2t)*
        pow4(Mst2)) + 10*pow2(Sbeta)*pow4(Mst2)*(-160*pow2(Mst2)*pow2(Mt)*pow2(
        s2t) + 72*pow4(Mt) + 63*pow4(Mst2)*pow4(s2t)) + pow4(Mst1)*(-8*pow2(Mt)
        *pow2(s2t)*(1471*pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 836*pow2(Mst2)*pow2(
        Sbeta)) - 23456*Cbeta*MuSUSY*s2t*Sbeta*pow3(Mt) - 68*Cbeta*Mt*MuSUSY*
        Sbeta*pow2(Mst2)*pow3(s2t) + 296*pow2(Sbeta)*pow4(Mt) - 807*pow2(Sbeta)
        *pow4(Mst2)*pow4(s2t)) + 16*Mt*s2t*Sbeta*(-1317*Cbeta*Mt*MuSUSY*s2t +
        1724*Sbeta*pow2(Mt) - 21*Sbeta*pow2(Mst2)*pow2(s2t))*pow5(Mst1) - 2*
        pow2(Mst1)*(-8*pow2(Mst2)*pow2(Mt)*pow2(s2t)*(-727*pow2(MuSUSY)*(-1 +
        pow2(Sbeta)) + 415*pow2(Mst2)*pow2(Sbeta)) + 5040*Cbeta*MuSUSY*s2t*
        Sbeta*pow2(Mst2)*pow3(Mt) + 1648*Cbeta*Mt*MuSUSY*Sbeta*pow3(s2t)*pow4(
        Mst2) + 8*(320*pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 271*pow2(Mst2)*pow2(
        Sbeta))*pow4(Mt) - 97*pow2(Sbeta)*pow4(s2t)*pow6(Mst2))))))/(45.*pow2(
        Mgl)*pow2(Sbeta)*pow4(Msq)*pow4(Mst2)*pow4(Mt));

   return result;
}

/**
 * @return returns the susy log^2 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double H5g1::calc_coef_at_as2_no_sm_logs_log2() const {

   const double result =
      (14760*pow2(Mgl)*pow2(Mst1)*pow2(Mst2)*pow2(Mt)*pow2(MuSUSY)*pow2(s2t) -
        14760*pow2(Mgl)*pow2(Mst1)*pow2(Mst2)*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*
        pow2(Sbeta) - 34560*Cbeta*Dmglst1*Mgl*MuSUSY*Sbeta*pow2(Mst2)*pow2(Mt)*
        pow2(s2t)*pow3(Mst1) - 34560*Cbeta*MuSUSY*Sbeta*pow2(Dmglst1)*pow2(
        Mst2)*pow2(Mt)*pow2(s2t)*pow3(Mst1) - 11520*Cbeta*MuSUSY*Sbeta*pow2(
        Mgl)*pow2(Mst2)*pow2(Mt)*pow2(s2t)*pow3(Mst1) - 67920*Cbeta*MuSUSY*s2t*
        Sbeta*pow2(Mgl)*pow2(Mst1)*pow2(Mst2)*pow3(Mt) + 170880*Dmglst1*Mgl*
        s2t*pow2(MuSUSY)*pow3(Mst1)*pow3(Mt) + 170880*s2t*pow2(Dmglst1)*pow2(
        MuSUSY)*pow3(Mst1)*pow3(Mt) + 109440*s2t*pow2(Mgl)*pow2(MuSUSY)*pow3(
        Mst1)*pow3(Mt) + 136320*Dmglst1*Mgl*s2t*pow2(Mst2)*pow2(Sbeta)*pow3(
        Mst1)*pow3(Mt) + 399360*s2t*pow2(Dmglst1)*pow2(Mst2)*pow2(Sbeta)*pow3(
        Mst1)*pow3(Mt) + 17280*s2t*pow2(Mgl)*pow2(Mst2)*pow2(Sbeta)*pow3(Mst1)*
        pow3(Mt) - 170880*Dmglst1*Mgl*s2t*pow2(MuSUSY)*pow2(Sbeta)*pow3(Mst1)*
        pow3(Mt) - 170880*s2t*pow2(Dmglst1)*pow2(MuSUSY)*pow2(Sbeta)*pow3(Mst1)
        *pow3(Mt) - 109440*s2t*pow2(Mgl)*pow2(MuSUSY)*pow2(Sbeta)*pow3(Mst1)*
        pow3(Mt) - 7680*Dmglst1*Mgl*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*pow4(Mst1)
        - 19200*pow2(Dmglst1)*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*pow4(Mst1) -
        1920*pow2(Mgl)*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*pow4(Mst1) - 29520*
        Dmglst1*Mgl*pow2(Mst2)*pow2(Mt)*pow2(s2t)*pow2(Sbeta)*pow4(Mst1) -
        51960*pow2(Dmglst1)*pow2(Mst2)*pow2(Mt)*pow2(s2t)*pow2(Sbeta)*pow4(
        Mst1) - 12840*pow2(Mgl)*pow2(Mst2)*pow2(Mt)*pow2(s2t)*pow2(Sbeta)*pow4(
        Mst1) + 7680*Dmglst1*Mgl*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*pow2(Sbeta)*
        pow4(Mst1) + 19200*pow2(Dmglst1)*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*pow2(
        Sbeta)*pow4(Mst1) + 1920*pow2(Mgl)*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*
        pow2(Sbeta)*pow4(Mst1) - 76800*Cbeta*Dmglst1*Mgl*MuSUSY*s2t*Sbeta*pow3(
        Mt)*pow4(Mst1) - 99840*Cbeta*MuSUSY*s2t*Sbeta*pow2(Dmglst1)*pow3(Mt)*
        pow4(Mst1) - 42240*Cbeta*MuSUSY*s2t*Sbeta*pow2(Mgl)*pow3(Mt)*pow4(Mst1)
        + 22440*Cbeta*Dmglst1*Mgl*Mt*MuSUSY*Sbeta*pow2(Mst2)*pow3(s2t)*pow4(
        Mst1) + 45180*Cbeta*Mt*MuSUSY*Sbeta*pow2(Dmglst1)*pow2(Mst2)*pow3(s2t)*
        pow4(Mst1) + 8340*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mgl)*pow2(Mst2)*pow3(s2t)*
        pow4(Mst1) + 82080*Dmglst1*Mgl*pow2(Mst1)*pow2(Mt)*pow2(s2t)*pow2(
        Sbeta)*pow4(Mst2) + 115440*pow2(Dmglst1)*pow2(Mst1)*pow2(Mt)*pow2(s2t)*
        pow2(Sbeta)*pow4(Mst2) + 42960*pow2(Mgl)*pow2(Mst1)*pow2(Mt)*pow2(s2t)*
        pow2(Sbeta)*pow4(Mst2) + 36160*Dmglst1*Mgl*Mst1*s2t*pow2(Sbeta)*pow3(
        Mt)*pow4(Mst2) + 24000*Mst1*s2t*pow2(Dmglst1)*pow2(Sbeta)*pow3(Mt)*
        pow4(Mst2) + 25920*Mst1*s2t*pow2(Mgl)*pow2(Sbeta)*pow3(Mt)*pow4(Mst2) +
        6960*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mgl)*pow2(Mst1)*pow3(s2t)*pow4(Mst2) -
        19680*Dmglst1*Mgl*Mt*pow2(Sbeta)*pow3(Mst1)*pow3(s2t)*pow4(Mst2) -
        19680*Mt*pow2(Dmglst1)*pow2(Sbeta)*pow3(Mst1)*pow3(s2t)*pow4(Mst2) -
        19680*Mt*pow2(Mgl)*pow2(Sbeta)*pow3(Mst1)*pow3(s2t)*pow4(Mst2) + 30720*
        pow2(Mgl)*pow2(Mst1)*pow2(MuSUSY)*pow4(Mt) + 86880*Dmglst1*Mgl*pow2(
        Mst1)*pow2(Mst2)*pow2(Sbeta)*pow4(Mt) + 171280*pow2(Dmglst1)*pow2(Mst1)
        *pow2(Mst2)*pow2(Sbeta)*pow4(Mt) + 28080*pow2(Mgl)*pow2(Mst1)*pow2(
        Mst2)*pow2(Sbeta)*pow4(Mt) - 30720*pow2(Mgl)*pow2(Mst1)*pow2(MuSUSY)*
        pow2(Sbeta)*pow4(Mt) - 172480*Cbeta*Dmglst1*Mgl*MuSUSY*Sbeta*pow3(Mst1)
        *pow4(Mt) - 172480*Cbeta*MuSUSY*Sbeta*pow2(Dmglst1)*pow3(Mst1)*pow4(Mt)
        - 43200*Cbeta*MuSUSY*Sbeta*pow2(Mgl)*pow3(Mst1)*pow4(Mt) - 18240*
        Dmglst1*Mgl*pow2(Sbeta)*pow4(Mst1)*pow4(Mt) - 75680*pow2(Dmglst1)*pow2(
        Sbeta)*pow4(Mst1)*pow4(Mt) + 480*pow2(Mgl)*pow2(Sbeta)*pow4(Mst1)*pow4(
        Mt) + 29600*Dmglst1*Mgl*pow2(Sbeta)*pow4(Mst2)*pow4(Mt) + 19856*pow2(
        Dmglst1)*pow2(Sbeta)*pow4(Mst2)*pow4(Mt) + 29440*pow2(Mgl)*pow2(Sbeta)*
        pow4(Mst2)*pow4(Mt) - 7200*log(pow2(Mst1)/pow2(Msq))*pow2(Mgl)*pow2(
        Sbeta)*pow4(Mst2)*pow4(Mt) - 8490*Dmglst1*Mgl*pow2(Sbeta)*pow4(Mst1)*
        pow4(Mst2)*pow4(s2t) - 18495*pow2(Dmglst1)*pow2(Sbeta)*pow4(Mst1)*pow4(
        Mst2)*pow4(s2t) - 345*pow2(Mgl)*pow2(Sbeta)*pow4(Mst1)*pow4(Mst2)*pow4(
        s2t) + 97920*Dmglst1*Mgl*s2t*pow2(Sbeta)*pow3(Mt)*pow5(Mst1) + 514560*
        s2t*pow2(Dmglst1)*pow2(Sbeta)*pow3(Mt)*pow5(Mst1) + 1920*s2t*pow2(Mgl)*
        pow2(Sbeta)*pow3(Mt)*pow5(Mst1) - 11520*Dmglst1*Mgl*Mt*pow2(Mst2)*pow2(
        Sbeta)*pow3(s2t)*pow5(Mst1) - 23040*Mt*pow2(Dmglst1)*pow2(Mst2)*pow2(
        Sbeta)*pow3(s2t)*pow5(Mst1) - 3840*Mt*pow2(Mgl)*pow2(Mst2)*pow2(Sbeta)*
        pow3(s2t)*pow5(Mst1) - 14160*Dmglst1*Mgl*pow2(Mt)*pow2(s2t)*pow2(Sbeta)
        *pow6(Mst2) - 13560*pow2(Dmglst1)*pow2(Mt)*pow2(s2t)*pow2(Sbeta)*pow6(
        Mst2) - 9000*pow2(Mgl)*pow2(Mt)*pow2(s2t)*pow2(Sbeta)*pow6(Mst2) +
        31200*Dmglst1*Mgl*Mst1*Mt*pow2(Sbeta)*pow3(s2t)*pow6(Mst2) + 42720*
        Mst1*Mt*pow2(Dmglst1)*pow2(Sbeta)*pow3(s2t)*pow6(Mst2) + 23520*Mst1*Mt*
        pow2(Mgl)*pow2(Sbeta)*pow3(s2t)*pow6(Mst2) + 1110*Dmglst1*Mgl*pow2(
        Mst1)*pow2(Sbeta)*pow4(s2t)*pow6(Mst2) + 5505*pow2(Dmglst1)*pow2(Mst1)*
        pow2(Sbeta)*pow4(s2t)*pow6(Mst2) - 5325*pow2(Mgl)*pow2(Mst1)*pow2(
        Sbeta)*pow4(s2t)*pow6(Mst2) + 30*log(pow2(Mst1)/pow2(Mst2))*(Mst1*pow2(
        Dmglst1)*(48*s2t*Sbeta*(-137*Cbeta*MuSUSY*s2t + 460*Mt*Sbeta)*pow2(Mt)*
        pow4(Mst1) + 8*Mt*s2t*(-60*pow2(Mt) + 41*pow2(Mst2)*pow2(s2t))*pow2(
        Sbeta)*pow4(Mst2) - 8*Mt*pow2(Mst1)*(534*Cbeta*Mt*MuSUSY*Sbeta*pow2(
        Mst2)*pow2(s2t) + 12*s2t*pow2(Mt)*(73*pow2(MuSUSY)*(-1 + pow2(Sbeta)) -
        88*pow2(Mst2)*pow2(Sbeta)) + 468*Cbeta*MuSUSY*Sbeta*pow3(Mt) - 233*
        pow2(Sbeta)*pow3(s2t)*pow4(Mst2)) + 3*Mst1*pow2(Mst2)*pow2(Sbeta)*(384*
        pow2(Mst2)*pow2(Mt)*pow2(s2t) + 512*pow4(Mt) + 91*pow4(Mst2)*pow4(s2t))
        + pow3(Mst1)*(-24*pow2(Mt)*pow2(s2t)*(91*pow2(MuSUSY)*(-1 + pow2(Sbeta)
        ) - 96*pow2(Mst2)*pow2(Sbeta)) - 6912*Cbeta*MuSUSY*s2t*Sbeta*pow3(Mt) +
        3296*pow2(Sbeta)*pow4(Mt) - 273*pow2(Sbeta)*pow4(Mst2)*pow4(s2t))) + 2*
        Dmglst1*Mgl*Mst1*(8*s2t*Sbeta*(-267*Cbeta*MuSUSY*s2t + 460*Mt*Sbeta)*
        pow2(Mt)*pow4(Mst1) + 4*Mt*s2t*(-60*pow2(Mt) + 41*pow2(Mst2)*pow2(s2t))
        *pow2(Sbeta)*pow4(Mst2) - 4*Mt*pow2(Mst1)*(534*Cbeta*Mt*MuSUSY*Sbeta*
        pow2(Mst2)*pow2(s2t) + 12*s2t*pow2(Mt)*(73*pow2(MuSUSY)*(-1 + pow2(
        Sbeta)) - 44*pow2(Mst2)*pow2(Sbeta)) + 468*Cbeta*MuSUSY*Sbeta*pow3(Mt)
        - 137*pow2(Sbeta)*pow3(s2t)*pow4(Mst2)) + Mst1*pow2(Mst2)*pow2(Sbeta)*(
        384*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 512*pow4(Mt) + 91*pow4(Mst2)*pow4(
        s2t)) + pow3(Mst1)*(-8*pow2(Mt)*pow2(s2t)*(91*pow2(MuSUSY)*(-1 + pow2(
        Sbeta)) - 96*pow2(Mst2)*pow2(Sbeta)) - 2304*Cbeta*MuSUSY*s2t*Sbeta*
        pow3(Mt) + 864*pow2(Sbeta)*pow4(Mt) - 91*pow2(Sbeta)*pow4(Mst2)*pow4(
        s2t))) + pow2(Mgl)*(8*Mst1*Mt*s2t*(-60*pow2(Mt) + 41*pow2(Mst2)*pow2(
        s2t))*pow2(Sbeta)*pow4(Mst2) - 8*Mt*pow3(Mst1)*(342*Cbeta*Mt*MuSUSY*
        Sbeta*pow2(Mst2)*pow2(s2t) + 4*s2t*pow2(Mt)*(155*pow2(MuSUSY)*(-1 +
        pow2(Sbeta)) - 44*pow2(Mst2)*pow2(Sbeta)) + 116*Cbeta*MuSUSY*Sbeta*
        pow3(Mt) - 73*pow2(Sbeta)*pow3(s2t)*pow4(Mst2)) + pow2(Sbeta)*pow4(
        Mst2)*(-328*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 408*pow4(Mt) + 41*pow4(
        Mst2)*pow4(s2t)) + 4*pow4(Mst1)*(2*pow2(Mt)*pow2(s2t)*(-173*pow2(
        MuSUSY)*(-1 + pow2(Sbeta)) + 96*pow2(Mst2)*pow2(Sbeta)) - 576*Cbeta*
        MuSUSY*s2t*Sbeta*pow3(Mt) + 172*pow2(Sbeta)*pow4(Mt) - 33*pow2(Sbeta)*
        pow4(Mst2)*pow4(s2t)) + 16*s2t*Sbeta*(-171*Cbeta*MuSUSY*s2t + 92*Mt*
        Sbeta)*pow2(Mt)*pow5(Mst1) + pow2(Mst1)*(8*pow2(Mst2)*pow2(Mt)*pow2(
        s2t)*(-173*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 89*pow2(Mst2)*pow2(Sbeta))
        - 768*Cbeta*MuSUSY*s2t*Sbeta*pow2(Mst2)*pow3(Mt) - 528*Cbeta*Mt*MuSUSY*
        Sbeta*pow3(s2t)*pow4(Mst2) + 512*(-(pow2(MuSUSY)*(-1 + pow2(Sbeta))) +
        pow2(Mst2)*pow2(Sbeta))*pow4(Mt) + 91*pow2(Sbeta)*pow4(s2t)*pow6(Mst2))
        )) + 1770*Dmglst1*Mgl*pow2(Sbeta)*pow4(s2t)*pow8(Mst2) + 1695*pow2(
        Dmglst1)*pow2(Sbeta)*pow4(s2t)*pow8(Mst2) + 3585*pow2(Mgl)*pow2(Sbeta)*
        pow4(s2t)*pow8(Mst2))/(45.*pow2(Mgl)*pow2(Sbeta)*pow4(Mst2)*pow4(Mt));

   return result;
}

/**
 * @return returns the susy log^3 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double H5g1::calc_coef_at_as2_no_sm_logs_log3() const {

   const double result =
      -298.6666666666667;

   return result;
}

} // namespace hierarchies
} // namespace himalaya
