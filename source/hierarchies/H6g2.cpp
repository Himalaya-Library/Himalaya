// ====================================================================
// This file is part of Himalaya.
//
// Himalaya is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================

#include "H6g2.hpp"
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
H6g2::H6g2(const ExpansionFlags_t& expansionDepth, double Al4p, double beta, double Dmglst2,
                 double lmMt, double lmMst1, double lmMst2, double lmMsq,
                 double Mgl, double Mt, double Mst1, double Mst2, double Msq, double MuSUSY,
                 double s2t,
                 int mdrFlag, int oneLoopFlag, int twoLoopFlag, int threeLoopFlag){
   // abbrev for tan(beta) and sin(beta)
   Tbeta = tan(beta);
   Sbeta = sin(beta);
   this -> Dmglst2 = Dmglst2;
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
   xDR2DRMOD = mdrFlag;
   // expansion flags
   xDmglst2 = expansionDepth.at(ExpansionDepth::Dmglst2);
   xMsq = expansionDepth.at(ExpansionDepth::Msq);
   xMst = expansionDepth.at(ExpansionDepth::Mst);
}

/**
 * @return The diagonal (1, 1) matrix element of the Higgs mass matrix as a double for the hierarchy 'H6g2'
 */
double H6g2::getS1() const {
   return -(pow2(Mt)*((270*s2t*xMst*(27*(lmMst1 - lmMst2)*Mst2*oneLoopFlag*s2t*
        pow2(MuSUSY) - (2*Al4p*twoLoopFlag*(4*(18*(lmMst1 - lmMst2)*(-2 + 3*
        lmMst2)*Mst2*s2t*xDmglst2*xDR2DRMOD*pow2(Dmglst2) + Dmglst2*Mgl*(Mt*(
        785 + 6*lmMst1*(85 - 42*lmMst2) - 438*lmMst2 + 252*pow2(lmMst2)) -
        Mst2*s2t*(49 - 84*lmMst2 + lmMst1*(84 - 36*lmMst2*(-1 + xDR2DRMOD)) +
        36*(-1 + xDR2DRMOD)*pow2(lmMst2))) + (Mt*(193 + 474*lmMst2 - 6*lmMst1*(
        67 + 42*lmMst2) + 252*pow2(lmMst2)) - Mst2*s2t*(1 + 3*lmMst2*(-37 + 6*
        xDR2DRMOD) - 3*lmMst1*(-37 + 6*lmMst2*(-12 + xDR2DRMOD) + 6*xDR2DRMOD)
        - 81*pow2(lmMst1) + 9*(-15 + 2*xDR2DRMOD)*pow2(lmMst2)))*pow2(Mgl))*
        pow2(MuSUSY) - (xDmglst2*pow2(Dmglst2)*(4*Mt*Tbeta*(-2501 + 222*lmMst2
        + 42*lmMst1*(-7 + 6*lmMst2) - 252*pow2(lmMst2))*pow2(MuSUSY) + 6*Mst2*
        MuSUSY*(MuSUSY*s2t*Tbeta*(85 - 60*lmMst1 + 60*lmMst2 + 36*lmMst1*lmMst2
        - 36*pow2(lmMst2)) + 10*(-43 + 60*lmMst1 - 60*lmMst2)*Mt*pow2(Sbeta)) +
        15*(-43 + 60*lmMst1 - 60*lmMst2)*s2t*Tbeta*(-1 + pow2(Sbeta))*pow2(
        Sbeta)*pow3(Mst2)))/Tbeta))/pow2(Mgl))*pow6(Mst1))/pow7(Mst2) - pow2(
        MuSUSY)*(3645*oneLoopFlag*pow2(s2t)*(2 - lmMst1 + lmMst2 - (2*(lmMst1 -
        lmMst2)*pow2(Mst1)*(pow2(Mst1) + pow2(Mst2)))/pow4(Mst2)) + (19440*
        Al4p*pow2(s2t)*(((2 - 3*lmMst2)*twoLoopFlag*xDmglst2*xDR2DRMOD*pow2(
        Dmglst2)*((pow2(Mst1) + pow2(Mst2))*pow4(Mst2) - 2*(lmMst1 - lmMst2)*
        pow2(Mst1)*(pow2(Mst1)*pow2(Mst2) + pow4(Mst1) + pow4(Mst2))))/pow2(
        Mgl) + 5*Al4p*(1 - 2*lmMsq)*threeLoopFlag*xMsq*pow2(Msq)*((1 - 2*
        lmMst2*(-1 + shiftst1))*pow2(Mst1)*pow2(Mst2) + 2*lmMst1*(-1 +
        shiftst1)*pow2(Mst1)*pow2(Mst2) + shiftst2*pow2(Mst1)*(2*pow2(Mst1) +
        pow2(Mst2)) - 2*lmMst1*(1 - 2*shiftst1 + shiftst2)*pow4(Mst1) + 2*
        lmMst2*(1 - 2*shiftst1 + shiftst2)*pow4(Mst1) + pow4(Mst2) - shiftst1*(
        2*pow2(Mst1)*pow2(Mst2) + 2*pow4(Mst1) + pow4(Mst2)))))/(pow2(Mst1)*
        pow4(Mst2)) + 2*threeLoopFlag*pow2(Al4p)*((180*pow2(Mt)*(3*Mgl*(8*pow2(
        Mst2)*(436 - 6*B4 + 6*D3 - 3*DN - 48*lmMst1 + (408 - 96*lmMst1)*lmMst2
        + 24*lmMt - 972*S2 - 48*(-4 + lmMst1)*pow2(lmMst2) + 48*pow3(lmMst2)) +
        pow2(Mst1)*(10667 - 96*B4 + 96*D3 - 48*DN - 3072*lmMst1 - 384*(-25 + 9*
        lmMst1)*lmMst2 - 384*lmMst1*lmMt + 384*(1 + lmMst2)*lmMt - 224*OepS2 +
        324*(-43 + 14*lmMst1 - 14*lmMst2)*S2 - 384*(-13 + 4*lmMst1)*pow2(
        lmMst2) + 1536*pow3(lmMst2))) + 2*Dmglst2*(pow2(Mst1)*(28405 - 288*B4 +
        288*D3 - 144*DN + 5472*lmMst1 - 288*(-89 + 2*lmMst1)*lmMst2 + 576*(-2 +
        lmMst1 - lmMst2)*lmMt + 224*OepS2 - 324*(65 + 14*lmMst1 - 14*lmMst2)*S2
        - 576*(-9 + 8*lmMst1)*pow2(lmMst2) + 4608*pow3(lmMst2)) + 72*pow2(Mst2)
        *(180 - 2*B4 + 2*D3 - DN + 144*lmMst2 - 216*S2 - 16*(-2 + lmMst1)*pow2(
        lmMst2) + 16*(lmMst1 + pow3(lmMst2))))))/(Mgl*pow4(Mst2)) + 90*pow2(
        s2t)*(15829 - 720*B4 + 72*D3 - 36*DN - 3330*lmMsq + 1350*pow2(lmMsq) -
        3*lmMst1*(5279 - 1950*lmMsq + 630*pow2(lmMsq)) - 954*pow2(lmMst1) + 3*
        lmMst2*(10421 - 6558*lmMst1 + 30*lmMsq*(-95 + 42*lmMst1) + 630*pow2(
        lmMsq) + 234*pow2(lmMst1)) - 18*(-1460 + 210*lmMsq + 399*lmMst1)*pow2(
        lmMst2) - 288*pow3(lmMst1) + 6768*pow3(lmMst2) + (6*pow2(Mst1)*((Mgl*(
        111 + 2668*lmMst1 + 2*(-1514 + 675*lmMst1)*lmMst2 + 90*lmMsq*(4 - 13*
        lmMst1 + 13*lmMst2) - 90*pow2(lmMst1) - 1260*pow2(lmMst2)) + 100*
        Dmglst2*(69 + (-53 + 42*lmMsq)*lmMst2 + lmMst1*(53 - 42*lmMsq + 42*
        lmMst2) - 42*pow2(lmMst2)))/(Mgl*pow2(Msq)) - (135*(204.20053497942388
         + (76*B4)/
        9. - (2*DN)/9. - (50*lmMsq)/3. + (2*lmMst1*(57971 - 14625*lmMsq + 2700*
        pow2(lmMsq)))/405. + ((-1331 + 180*lmMsq)*pow2(lmMst1))/27. - (2*
        lmMst2*(52436 - 70455*lmMst1 + 225*lmMsq*(-65 + 36*lmMst1) + 2700*pow2(
        lmMsq) + 8280*pow2(lmMst1)))/405. + ((-8063 + 900*lmMsq + 3792*lmMst1)*
        pow2(lmMst2))/27. - (62*pow3(lmMst1))/27. - (2626*pow3(lmMst2))/27. - (
        Dmglst2*(109.11799176954733 - (8*B4)/3. + (32*D3)/9. - (20*DN)/9. + 80*
        lmMsq - lmMst1*(78.19061728395062 + 20*pow2(lmMsq)) - (2888*pow2(
        lmMst1))/135. + lmMst2*(40*lmMsq*lmMst1 + 20*pow2(lmMsq) + (4*(-21616 -
        64515*lmMst1 + 31275*pow2(lmMst1)))/2025.) - (4*(-5023 + 1350*lmMsq +
        6285*lmMst1)*pow2(lmMst2))/135. + (20*pow3(lmMst1))/27. + (3340*pow3(
        lmMst2))/27.))/Mgl))/pow2(Mst2)))/5. + 162*pow4(Mst1)*((
        1.0702990137854083 + (9571*lmMsq)/26460. + (4.249508692365835 - (169*
        lmMsq)/63.)*lmMst1 + (-4.6112244897959185 + (19*lmMsq)/7. + (31*lmMst1)
        /9.)*lmMst2 - pow2(lmMsq)/63. - (8*pow2(lmMst1))/21. + (5*Dmglst2*(216
        + (-95 + 132*lmMsq)*lmMst2 + lmMst1*(95 - 132*lmMsq + 132*lmMst2) -
        132*pow2(lmMst2)))/(54.*Mgl) - (194*pow2(lmMst2))/63.)/pow4(Msq) - (
        363.3804294212688 + (76*B4)/9. - (2*DN)/9. - (35*lmMsq)/2. + lmMst1*(
        211.3489518770471 - (695*lmMsq)/9. + (40*pow2(lmMsq))/3.) - (
        214.87936507936507 - 20*lmMsq)*pow2(lmMst1) - lmMst2*(
        190.46006298815823 - (71398*lmMst1)/105. + (5*lmMsq*(-139 + 120*lmMst1)
        )/9. + (40*pow2(lmMsq))/3. + (334*pow2(lmMst1))/3.) + ((-146507 +
        14700*lmMsq + 91070*lmMst1)*pow2(lmMst2))/315. - (44*pow3(lmMst1))/9. -
        (1556*pow3(lmMst2))/9. - (Dmglst2*(536.1152102791342 - (8*B4)/3. + (32*
        D3)/9. - (20*DN)/9. + 90*lmMsq - (123.11224321827497 + 20*lmMsq*(1 +
        lmMsq))*lmMst1 - lmMst2*(17.33220122616948 - 20*lmMsq*(1 + lmMsq) + (
        133.04550264550264 - 40*lmMsq)*lmMst1 - (1180*pow2(lmMst1))/9.) - (
        15886*pow2(lmMst1))/945. + (149.85608465608465 - 40*lmMsq - (2812*
        lmMst1)/9.)*pow2(lmMst2) - (92*pow3(lmMst1))/27. + (4988*pow3(lmMst2))/
        27.))/Mgl)/pow4(Mst2)) + (162*(Dmglst2*(19.734567901234566 - (8*B4)/3.
         + (32*D3)/
        9. - (20*DN)/9. + (170*lmMsq)/3. + lmMst1*(9.333333333333334 - 20*pow2(
        lmMsq)) + 10*pow2(lmMsq) - (26*pow2(lmMst1))/3. + (2*lmMst2*(-291 -
        464*lmMst1 + 90*lmMsq*(-1 + 2*lmMst1) + 90*pow2(lmMsq) + 32*pow2(
        lmMst1)))/9. - (2*(-607 + 180*lmMsq + 336*lmMst1)*pow2(lmMst2))/9. + (
        608*pow3(lmMst2))/9.) - (16*(1 + lmMst2)*(4*Dmglst2*lmMst2 + Mgl +
        lmMst2*Mgl)*pow4(Mst2))/(9.*pow4(Mst1)) + (pow2(Mst2)*((Mgl*(2*(4167613
        - 19932360*lmMst2 + 20580*lmMst1*(701 + 540*lmMst2) + 420*lmMsq*(13109
        - 26460*lmMst1 + 25620*lmMst2) + 176400*pow2(lmMsq) - 10936800*pow2(
        lmMst2))*pow2(Mst1) + (41220947 - 420*lmMsq*(12479 + 69090*lmMst1 -
        69930*lmMst2) - 21234990*lmMst2 + 10290*lmMst1*(2573 + 2820*lmMst2) -
        176400*pow2(lmMsq) - 29194200*pow2(lmMst2))*pow2(Mst2)))/1.11132e7 + (
        5*Dmglst2*(2*(219 + (-95 + 132*lmMsq)*lmMst2 + lmMst1*(95 - 132*lmMsq +
        132*lmMst2) - 132*pow2(lmMst2))*pow2(Mst1) + (557 - 224*lmMst2 + 4*
        lmMst1*(53 + 96*lmMst2) + lmMsq*(12 - 384*lmMst1 + 384*lmMst2) - 384*
        pow2(lmMst2))*pow2(Mst2)))/108. - (40*(8*Dmglst2*(14 - 15*lmMsq + 15*
        lmMst2) + (43 - 30*lmMsq + 30*lmMst2)*Mgl)*pow2(Msq)*pow2(Mst2) + (-12*
        Mgl*(341 + 642*lmMst2 + 64*lmMst1*(1 + lmMst2) - 30*lmMsq*(5 + 6*
        lmMst2) + 90*pow2(lmMsq) + 272*pow2(lmMst2)) - 24*Dmglst2*(-115 + (366
        + 64*lmMst1)*lmMst2 - 30*lmMsq*(-1 + 6*lmMst2) + 90*pow2(lmMsq) + 336*
        pow2(lmMst2)))*pow4(Msq) + (90*Dmglst2*(13 - 28*lmMsq + 28*lmMst2) + 5*
        (67 - 84*lmMsq + 84*lmMst2)*Mgl)*pow4(Mst2))/(216.*pow2(Mst1))))/pow4(
        Msq) - (2*Dmglst2*(54*(344*OepS2 + 9*(15643 - 774*lmMst1 + 774*lmMst2)*
        S2)*pow2(Mst1)*pow2(Mst2) + 4*(17308*OepS2 + 27*(93919 - 12981*lmMst1 +
        12981*lmMst2)*S2)*pow4(Mst1) + 27*(56*OepS2 - 81*(-1677 + 14*lmMst1 -
        14*lmMst2)*S2)*pow4(Mst2)) + 3*Mgl*(-3*(3896*OepS2 - 81*(9473 + 974*
        lmMst1 - 974*lmMst2)*S2)*pow2(Mst1)*pow2(Mst2) + (-29428*OepS2 + 27*(
        160997 + 22071*lmMst1 - 22071*lmMst2)*S2)*pow4(Mst1) + 432*(-4*OepS2 +
        81*(22 + lmMst1 - lmMst2)*S2)*pow4(Mst2)))/(2187.*pow4(Mst2)) + (Mgl*((
        1725 + (-7006 + 2640*lmMsq)*lmMst2 + lmMst1*(7006 - 2640*lmMsq + 4800*
        lmMst2) - 1080*pow2(lmMst1) - 3720*pow2(lmMst2))*pow4(Mst1) + 2*(836 -
        2235*lmMst2 + 75*lmMst1*(27 + 16*lmMst2) + 30*lmMsq*(7 - 40*lmMst1 +
        40*lmMst2) - 1200*pow2(lmMst2))*pow4(Mst2)) + 50*Dmglst2*((291 + 2*(-
        103 + 84*lmMsq)*lmMst2 + 2*lmMst1*(103 - 84*lmMsq + 84*lmMst2) - 168*
        pow2(lmMst2))*pow4(Mst1) + 2*(118 + 109*lmMst1 + (-133 + 102*lmMst1)*
        lmMst2 + 6*lmMsq*(4 - 17*lmMst1 + 17*lmMst2) - 102*pow2(lmMst2))*pow4(
        Mst2)))/(270.*pow2(Msq)*pow2(Mst2))))/Mgl + (54*(-1 + 2*lmMst2)*
        shiftst3*(2*(1 - 2*lmMst1 + 2*lmMst2)*pow2(Mst2)*pow4(Mst1) + 2*(1 -
        lmMst1 + lmMst2)*pow2(Mst1)*pow4(Mst2) + (2 - 6*lmMst1 + 6*lmMst2)*
        pow6(Mst1) + pow6(Mst2)))/(pow2(Mst1)*pow4(Mst2))) + ((10*(3*pow2(Mgl)*
        (-6*Mst2*pow2(Mst1)*(-36*pow2(Mt)*(218*z3 + 50*z4 + 21*pow2(z2)) +
        pow2(Mst2)*pow2(s2t)*((-56978 - 216*lmMst1 + 216*lmMst2)*z3 + 1136*z4 +
        813*pow2(z2)) + 4*Mst2*Mt*s2t*(-11396*z3 - 1135*z4 + 2388*pow2(z2))) +
        18*(8*(230*z3 + 27*z4)*pow2(Mt) + 5*Mst2*Mt*s2t*(2098*z3 + 102*z4 -
        333*pow2(z2)) + 3*pow2(Mst2)*pow2(s2t)*((1913 + 12*lmMst1 - 12*lmMst2)*
        z3 - 10*z4 + 48*pow2(z2)))*pow3(Mst2) + s2t*(Mt*(268928*z3 + 47056*z4 -
        82020*pow2(z2)) + Mst2*s2t*(8*(76813 + 162*lmMst1 - 162*lmMst2)*z3 -
        15686*z4 - 18183*pow2(z2)))*pow4(Mst1)) + 2*Dmglst2*Mgl*Mst2*(27*pow2(
        Mst2)*(16*(-83*z3 + 27*z4)*pow2(Mt) + Mst2*Mt*s2t*(37474*z3 + 542*z4 -
        5289*pow2(z2)) + pow2(Mst2)*pow2(s2t)*(16175*z3 + 424*z4 + 474*pow2(z2)
        )) - 36*pow2(Mst1)*(6*pow2(Mt)*(250*z3 - 94*z4 + 21*pow2(z2)) - 3*pow2(
        Mst2)*pow2(s2t)*(9207*z3 + 185*z4 + 237*pow2(z2)) + Mst2*Mt*s2t*(-
        59336*z3 - 709*z4 + 8535*pow2(z2))) + 4*pow2(s2t)*(338536*z3 + 11327*z4
        + 15897*pow2(z2))*pow4(Mst1)) + Mst2*xDmglst2*pow2(Dmglst2)*(9*pow2(
        Mst2)*(48*(-1990*z3 + 81*z4)*pow2(Mt) + 2*Mst2*Mt*s2t*(256474*z3 +
        2498*z4 - 31083*pow2(z2)) + pow2(Mst2)*pow2(s2t)*(226447*z3 + 8888*z4 +
        4098*pow2(z2))) - 12*pow2(Mst1)*(36*pow2(Mt)*(250*z3 - 94*z4 + 21*pow2(
        z2)) - pow2(Mst2)*pow2(s2t)*(409448*z3 + 10819*z4 + 3471*pow2(z2)) +
        Mst2*Mt*s2t*(-699292*z3 - 12344*z4 + 88647*pow2(z2))) + 8*pow2(s2t)*(
        338536*z3 + 11327*z4 + 15897*pow2(z2))*pow4(Mst1))))/pow2(Mgl) + (40*
        T1ep*(2*Dmglst2*Mst2*(-36*pow2(Mst1)*(253*Mst2*Mt*s2t + 42*pow2(Mt) -
        129*pow2(Mst2)*pow2(s2t)) + 189*s2t*(-5*Mt + 2*Mst2*s2t)*pow3(Mst2) +
        17308*pow2(s2t)*pow4(Mst1)) + 3*Mgl*(6*Mst2*pow2(Mst1)*(272*Mst2*Mt*s2t
        + 252*pow2(Mt) - 487*pow2(Mst2)*pow2(s2t)) + s2t*(3764*Mt - 7357*Mst2*
        s2t)*pow4(Mst1) + 54*s2t*(7*Mt - 8*Mst2*s2t)*pow4(Mst2))))/Mgl + Mt*
        s2t*(14580*pow2(Mst1)*pow2(Mst2)*(1035.3004115226338 + (1016*B4)/9. - (
        4*DN)/9. - 240*lmMsq + (80*pow2(lmMsq))/3. - (8*lmMst1*(422 - 45*lmMsq
        + 45*pow2(lmMsq)))/9. - (176*pow2(lmMst1))/9. + (8*lmMst2*(1096 + 15*
        lmMsq*(-7 + 6*lmMst1 - 6*lmMst2) + 717*lmMst2 - lmMst1*(551 + 248*
        lmMst2) + 45*pow2(lmMsq) + 8*pow2(lmMst1)))/9. + (640*pow3(lmMst2))/3.
         + ((-80*(Mgl*(55 + 6*lmMsq*(-2 + 5*lmMst1 - 5*lmMst2) + 52*lmMst2 -
        10*lmMst1*(4 + 3*lmMst2) + 30*pow2(lmMst2)) + Dmglst2*(321 + 18*lmMsq*(
        -2 + 5*lmMst1 - 5*lmMst2) + 96*lmMst2 - 30*lmMst1*(2 + 3*lmMst2) + 90*
        pow2(lmMst2)))*pow2(Mst1))/(27.*pow2(Msq)) + Dmglst2*(881.6139917695473
         + 248*B4 - 4*DN - (2560*lmMsq)/3. + lmMst1*(130.96296296296296 - 120*
        lmMsq - 40*pow2(lmMsq)) + (80*pow2(lmMsq))/3. + (176*pow2(lmMst1))/9. +
        (8*lmMst2*(4364 - 573*lmMst1 + 45*lmMsq*(5 + 6*lmMst1) + 135*pow2(
        lmMsq) + 24*pow2(lmMst1)))/27. - (8*(-377 + 90*lmMsq + 376*lmMst1)*
        pow2(lmMst2))/9. + (2944*pow3(lmMst2))/9.))/Mgl) + 10*(2552929 +
        257904*B4 - 648*DN - 456840*lmMsq + 38880*pow2(lmMsq) - 216*lmMst1*(
        4591 - 360*lmMsq + 450*pow2(lmMsq)) + 41904*pow2(lmMst1) + 216*lmMst2*(
        9211 - 6466*lmMst1 + 180*lmMsq*(-4 + 5*lmMst1) + 450*pow2(lmMsq) + 576*
        pow2(lmMst1)) - 864*(-1784 + 225*lmMsq + 840*lmMst1)*pow2(lmMst2) -
        3456*pow3(lmMst1) + 604800*pow3(lmMst2))*pow4(Mst1) + ((9*(15*Mgl*(-
        640*(23 + lmMsq*(-6 + 9*lmMst1 - 9*lmMst2) + 18*lmMst2 - 3*lmMst1*(4 +
        3*lmMst2) + 9*pow2(lmMst2))*pow2(Msq)*pow2(Mst1) + 3*(18201 + 1760*B4 -
        16*DN - 5760*lmMsq + 960*pow2(lmMsq) - 16*lmMst1*(199 - 30*lmMsq + 30*
        pow2(lmMsq)) + 16*lmMst2*(1291 - 322*lmMst1 + 30*lmMsq*(-5 + 2*lmMst1)
        + 30*pow2(lmMsq) - 16*pow2(lmMst1)) - 256*pow2(lmMst1) - 32*(-313 + 30*
        lmMsq + 72*lmMst1)*pow2(lmMst2) + 2560*pow3(lmMst2))*pow4(Msq) - 20*(
        233 + 36*lmMsq*(-2 + 5*lmMst1 - 5*lmMst2) + 207*lmMst2 - 45*lmMst1*(3 +
        4*lmMst2) + 180*pow2(lmMst2))*pow4(Mst1)) + Dmglst2*(-14400*(77 + 6*
        lmMsq*(-2 + 3*lmMst1 - 3*lmMst2) + 24*lmMst2 - 6*lmMst1*(2 + 3*lmMst2)
        + 18*pow2(lmMst2))*pow2(Msq)*pow2(Mst1) + (779917 + 188640*B4 - 3600*DN
        - 648000*lmMsq + 43200*pow2(lmMsq) - 720*lmMst1*(-173 + 90*lmMsq + 30*
        pow2(lmMsq)) + 720*lmMst2*(1623 - 130*lmMst1 + 30*lmMsq*(-1 + 2*lmMst1)
        + 30*pow2(lmMsq) - 16*pow2(lmMst1)) + 11520*pow2(lmMst1) - 1440*(-265 +
        30*lmMsq + 136*lmMst1)*pow2(lmMst2) + 207360*pow3(lmMst2))*pow4(Msq) -
        300*(1961 + 180*lmMsq*(-2 + 5*lmMst1 - 5*lmMst2) + 675*lmMst2 - 45*
        lmMst1*(7 + 20*lmMst2) + 900*pow2(lmMst2))*pow4(Mst1)))*pow4(Mst2))/
        pow4(Msq) - 160*OepS2*((-3036*Dmglst2 + 816*Mgl)*pow2(Mst1)*pow2(Mst2)
        + 1882*Mgl*pow4(Mst1) - 63*(5*Dmglst2 - 3*Mgl)*pow4(Mst2)) - 108*S2*(3*
        Dmglst2*pow2(Mst2)*(8*(2489 + 3795*lmMst1 - 3795*lmMst2)*pow2(Mst1) +
        9*(-453 + 350*lmMst1 - 350*lmMst2)*pow2(Mst2)) - 5*Mgl*(36*(169 + 136*
        lmMst1 - 136*lmMst2)*pow2(Mst1)*pow2(Mst2) + 2*(9185 + 5646*lmMst1 -
        5646*lmMst2)*pow4(Mst1) + 81*(-1 + 14*lmMst1 - 14*lmMst2)*pow4(Mst2)))
        + (90*(-48*(15*Dmglst2*(95 + 12*lmMsq*(-2 + lmMst1 - lmMst2) + 32*
        lmMst2 - 4*lmMst1*(2 + 3*lmMst2) + 12*pow2(lmMst2)) + 5*Mgl*(71 + 12*
        lmMsq*(-2 + lmMst1 - lmMst2) + 40*lmMst2 - 4*lmMst1*(4 + 3*lmMst2) +
        12*pow2(lmMst2)))*pow2(Msq) - 5*(Mgl*(1153 + 12*lmMsq*(-35 + 54*lmMst1
        - 54*lmMst2) + 906*lmMst2 - 162*lmMst1*(3 + 4*lmMst2) + 648*pow2(
        lmMst2)) + Dmglst2*(8713 + 12*lmMsq*(-179 + 270*lmMst1 - 270*lmMst2) +
        3282*lmMst2 - 162*lmMst1*(7 + 20*lmMst2) + 3240*pow2(lmMst2)))*pow2(
        Mst1) - 5*(Mgl*(911 + 12*lmMsq*(-37 + 18*lmMst1 - 18*lmMst2) + 606*
        lmMst2 - 54*lmMst1*(3 + 4*lmMst2) + 216*pow2(lmMst2)) + Dmglst2*(5591 +
        12*lmMsq*(-181 + 90*lmMst1 - 90*lmMst2) + 2550*lmMst2 - 54*lmMst1*(7 +
        20*lmMst2) + 1080*pow2(lmMst2)))*pow2(Mst2) + (2304*(1 + lmMst2)*(-
        Dmglst2 + 3*Dmglst2*lmMst2 + Mgl + lmMst2*Mgl)*pow4(Msq))/pow2(Mst1))*
        pow6(Mst2))/pow4(Msq))/Mgl))/pow5(Mst2)) + (2*Al4p*xDmglst2*pow2(
        Dmglst2)*(Al4p*threeLoopFlag*((72*pow2(Mt)*(5*pow2(Mst1)*(28405 - 288*
        B4 + 288*D3 - 144*DN + 5472*lmMst1 - 288*(-89 + 2*lmMst1)*lmMst2 + 576*
        (-2 + lmMst1 - lmMst2)*lmMt + 224*OepS2 - 324*(65 + 14*lmMst1 - 14*
        lmMst2)*S2 - 576*(-9 + 8*lmMst1)*pow2(lmMst2) + 4608*pow3(lmMst2)) + 4*
        pow2(Mst2)*(50134 - 270*B4 + 270*D3 - 135*DN + 120*(271 + 24*lmMst1)*
        lmMst2 - 120*lmMt - 29646*S2 - 720*(-2 + 3*lmMst1)*pow2(lmMst2) + 2160*
        (lmMst1 + pow3(lmMst2)))) + 80*T1ep*(6*pow2(Mst1)*(1555*Mst2*Mt*s2t -
        252*pow2(Mt) + 185*pow2(Mst2)*pow2(s2t)) + 63*s2t*(Mt + 5*Mst2*s2t)*
        pow3(Mst2) + 17308*pow2(s2t)*pow4(Mst1)))/pow4(Mst2) + 14580*pow2(s2t)*
        (112.1664609053498 - (100*B4)/9. + (112*D3)/9. - (62*DN)/9. + (1055*
        lmMsq)/9. + 15*pow2(lmMsq) - (2*lmMst1*(-544 + 210*lmMsq + 135*pow2(
        lmMsq)))/9. - (103*pow2(lmMst1))/9. + (lmMst2*(-7921 - 1476*lmMst1 +
        90*lmMsq*(5 + 18*lmMst1) + 810*pow2(lmMsq) + 288*pow2(lmMst1)))/27. + (
        93.66666666666667 - 60*lmMsq - 112*lmMst1)*pow2(lmMst2) + pow2(Mst1)*((
        10*(405 + (-187 + 246*lmMsq)*lmMst2 + lmMst1*(187 - 246*lmMsq + 246*
        lmMst2) - 246*pow2(lmMst2)))/(27.*pow2(Msq)) - (73.78472976680384 + (
        164*B4)/9. - (176*D3)/9. + (94*DN)/9. - (260*lmMsq)/3. - 20*lmMsq*(2 +
        3*lmMst1)*lmMst2 - 30*lmMst2*pow2(lmMsq) + lmMst1*(-96.55703703703703 +
        40*lmMsq + 30*pow2(lmMsq)) + (7556*pow2(lmMst1))/135. - (2*lmMst2*(-
        99788 + 2005*lmMst1 + 32775*pow2(lmMst1)))/675. + (2*(-3377 + 4050*
        lmMsq + 19155*lmMst1)*pow2(lmMst2))/135. + (10*pow3(lmMst1))/27. - (
        5050*pow3(lmMst2))/27.)/pow2(Mst2)) + (304*pow3(lmMst2))/3. + pow4(
        Mst1)*((5*(216 + (-95 + 132*lmMsq)*lmMst2 + lmMst1*(95 - 132*lmMsq +
        132*lmMst2) - 132*pow2(lmMst2)))/(54.*pow4(Msq)) + (536.1152102791342 -
        (8*B4)/3. + (32*D3)/9. - (20*DN)/9. + 90*lmMsq - (123.11224321827497 +
        20*lmMsq*(1 + lmMsq))*lmMst1 - lmMst2*(17.33220122616948 - 20*lmMsq*(1
        + lmMsq) + (133.04550264550264 - 40*lmMsq)*lmMst1 - (1180*pow2(lmMst1))
        /9.) - (15886*pow2(lmMst1))/945. + (149.85608465608465 - 40*lmMsq - (
        2812*lmMst1)/9.)*pow2(lmMst2) - (92*pow3(lmMst1))/27. + (4988*pow3(
        lmMst2))/27.)/pow4(Mst2)) - (32*(-2 + lmMst2 + 5*pow2(lmMst2))*pow4(
        Mst2))/(9.*pow4(Mst1)) - (6*(7400*OepS2 + 27*(1089707 - 5550*lmMst1 +
        5550*lmMst2)*S2)*pow2(Mst1)*pow2(Mst2) + 40*(17308*OepS2 + 27*(93919 -
        12981*lmMst1 + 12981*lmMst2)*S2)*pow4(Mst1) + 9*(1400*OepS2 + 81*(
        116129 - 350*lmMst1 + 350*lmMst2)*S2)*pow4(Mst2))/(10935.*pow4(Mst2)) +
        (5*((291 + 2*(-103 + 84*lmMsq)*lmMst2 + 2*lmMst1*(103 - 84*lmMsq + 84*
        lmMst2) - 168*pow2(lmMst2))*pow4(Mst1) + (684 + 347*lmMst1 + 7*(-77 +
        78*lmMst1)*lmMst2 + 6*lmMsq*(32 - 91*lmMst1 + 91*lmMst2) - 546*pow2(
        lmMst2))*pow4(Mst2)))/(27.*pow2(Msq)*pow2(Mst2)) + (pow2(Mst2)*(800*(-8
        + 15*lmMsq - 15*lmMst2)*pow2(Msq)*pow2(Mst2) + 12*(-1149 + 266*lmMst2 +
        64*lmMst1*(-2 + 3*lmMst2) - 90*lmMsq*(-5 + 6*lmMst2) + 270*pow2(lmMsq)
        + 1264*pow2(lmMst2))*pow4(Msq) - 5*((-4239 + 12*lmMsq*(-47 + 244*lmMst1
        - 244*lmMst2) + 1324*lmMst2 - 8*lmMst1*(95 + 366*lmMst2) + 2928*pow2(
        lmMst2))*pow2(Mst1)*pow2(Mst2) + (-3834 + lmMst1*(-890 + 2328*lmMsq -
        2328*lmMst2) + (890 - 2328*lmMsq)*lmMst2 + 2328*pow2(lmMst2))*pow4(
        Mst1) - 63*(-5 + 28*lmMsq - 28*lmMst2)*pow4(Mst2))))/(216.*pow2(Mst1)*
        pow4(Msq))) - (Mt*s2t*(pow4(Msq)*(-3*pow2(Mst1)*pow2(Mst2)*(433447 +
        1058400*B4 - 23760*DN - 4255200*lmMsq - 1120*OepS2 + 324*(-1387 + 70*
        lmMst1 - 70*lmMst2)*S2 + 129600*pow2(lmMsq) - 2160*lmMst1*(-359 + 150*
        lmMsq + 30*pow2(lmMsq)) + 720*lmMst2*(8503 + 666*lmMst1 + 90*lmMsq*(1 +
        2*lmMst1) + 90*pow2(lmMsq) - 48*pow2(lmMst1)) + 69120*pow2(lmMst1) -
        4320*(-177 + 30*lmMsq + 232*lmMst1)*pow2(lmMst2) + 1036800*pow3(lmMst2)
        ) - 2*(10274911 + 3285360*B4 - 68040*DN - 13381200*lmMsq - 248800*OepS2
        + 108*(-20803 + 46650*lmMst1 - 46650*lmMst2)*S2 + 194400*pow2(lmMsq) -
        3600*lmMst1*(167 + 405*lmMsq + 81*pow2(lmMsq)) - 103680*pow2(lmMst1) +
        720*lmMst2*(30469 + 2709*lmMst1 + 135*lmMsq*(11 + 6*lmMst1) + 405*pow2(
        lmMsq) + 72*pow2(lmMst1)) - 6480*(-19 + 90*lmMsq + 568*lmMst1)*pow2(
        lmMst2) + 3628800*pow3(lmMst2))*pow4(Mst1) + 414720*(2 + lmMst2 - 3*
        pow2(lmMst2))*pow4(Mst2)) + 21600*pow2(Msq)*((1445 + 72*lmMsq*(-2 + 3*
        lmMst1 - 3*lmMst2) + 180*lmMst2 - 36*lmMst1*(1 + 6*lmMst2) + 216*pow2(
        lmMst2))*pow2(Mst2)*pow4(Mst1) + 2*(397 + 36*lmMsq*(-2 + lmMst1 -
        lmMst2) + 78*lmMst2 - 6*lmMst1*(1 + 6*lmMst2) + 36*pow2(lmMst2))*pow2(
        Mst1)*pow4(Mst2) + 6*(107 + 6*lmMsq*(-2 + 5*lmMst1 - 5*lmMst2) + 32*
        lmMst2 - 10*lmMst1*(2 + 3*lmMst2) + 30*pow2(lmMst2))*pow6(Mst1)) + 450*
        ((37213 + 12*lmMsq*(-539 + 810*lmMst1 - 810*lmMst2) + 6630*lmMst2 -
        162*lmMst1*(1 + 60*lmMst2) + 9720*pow2(lmMst2))*pow4(Mst1)*pow4(Mst2) +
        6*(1961 + 180*lmMsq*(-2 + 5*lmMst1 - 5*lmMst2) + 675*lmMst2 - 45*
        lmMst1*(7 + 20*lmMst2) + 900*pow2(lmMst2))*pow2(Mst2)*pow6(Mst1) + (
        21275 + 12*lmMsq*(-541 + 270*lmMst1 - 270*lmMst2) + 6546*lmMst2 - 54*
        lmMst1*(1 + 60*lmMst2) + 3240*pow2(lmMst2))*pow2(Mst1)*pow6(Mst2))))/(
        pow2(Mst1)*pow3(Mst2)*pow4(Msq))) + (4860*s2t*twoLoopFlag*(-2*(4*Mt*(-
        31 + 3*lmMst1*(-2 + lmMst2) + 4*lmMst2 - 3*pow2(lmMst2)) + 3*Mst2*s2t*(
        4 + lmMst2 + lmMst1*(-1 + 2*lmMst2) - 2*pow2(lmMst2)))*pow2(Mst2)*pow4(
        Mst1) - 2*(Mst2*s2t*(10 + lmMst1*(-3 + 6*lmMst2) - 6*pow2(lmMst2)) + 4*
        Mt*(-14 + lmMst1*(-2 + lmMst2) - pow2(lmMst2)))*pow2(Mst1)*pow4(Mst2) +
        (-3*Mst2*s2t*(9 + 4*lmMst1*(-1 + lmMst2) + 4*lmMst2 - 4*pow2(lmMst2)) +
        Mt*(398 + lmMst1*(68 - 40*lmMst2) - 52*lmMst2 + 40*pow2(lmMst2)))*pow6(
        Mst1) + 2*(-2 + 3*lmMst2)*s2t*pow7(Mst2)))/(pow2(Mst1)*pow5(Mst2))))/
        pow2(Mgl) + Al4p*s2t*((4860*twoLoopFlag*(2*Mgl*(-(Mst2*s2t*(-14 - 20*
        lmMst2 + 2*lmMst1*(5 + 3*lmMst2) + pow2(lmMst1) - 7*pow2(lmMst2))) + 8*
        Mt*(4 + 3*lmMst2 - lmMst1*(1 + lmMst2) + pow2(lmMst2)))*pow2(Mst1)*
        pow4(Mst2) + 4*Dmglst2*(Mst2*s2t*(-5 + 8*lmMst2 - 4*lmMst1*(2 + lmMst2)
        + 4*pow2(lmMst2)) + Mt*(65 + lmMst1*(34 - 20*lmMst2) - 26*lmMst2 + 20*
        pow2(lmMst2)))*pow6(Mst1) + Mgl*(Mst2*s2t*(-1 + 50*lmMst2 - 2*lmMst1*(
        25 + 32*lmMst2) + 20*pow2(lmMst1) + 44*pow2(lmMst2)) + Mt*(84 + 152*
        lmMst2 - 40*lmMst1*(3 + 2*lmMst2) + 80*pow2(lmMst2)))*pow6(Mst1) + 8*
        Dmglst2*((Mst2*s2t*(-2 + 3*lmMst2 - lmMst1*(3 + 2*lmMst2) + 2*pow2(
        lmMst2)) + Mt*(22 + lmMst1*(8 - 6*lmMst2) - 4*lmMst2 + 6*pow2(lmMst2)))
        *pow2(Mst2)*pow4(Mst1) + (-(Mst2*s2t*(1 + lmMst1 - 2*lmMst2 + 2*lmMst1*
        lmMst2 - 2*pow2(lmMst2))) + 2*Mt*(6 + lmMst1 + lmMst2 - lmMst1*lmMst2 +
        pow2(lmMst2)))*pow2(Mst1)*pow4(Mst2) + lmMst2*s2t*pow7(Mst2)) + 4*Mgl*(
        (4*Mt*(5 + 6*lmMst2 - lmMst1*(4 + 3*lmMst2) + 3*pow2(lmMst2)) + Mst2*
        s2t*(-1 + 13*lmMst2 - lmMst1*(13 + 8*lmMst2) + pow2(lmMst1) + 7*pow2(
        lmMst2)))*pow2(Mst2)*pow4(Mst1) + (1 + lmMst2)*s2t*pow7(Mst2))))/(Mgl*
        pow2(Mst1)*pow5(Mst2)) - 135*xDR2DRMOD*((144*(2*Dmglst2*lmMst2 + Mgl +
        lmMst2*Mgl)*s2t*twoLoopFlag*((pow2(Mst1) + pow2(Mst2))*pow4(Mst2) - 2*(
        lmMst1 - lmMst2)*pow2(Mst1)*(pow2(Mst1)*pow2(Mst2) + pow4(Mst1) + pow4(
        Mst2))))/(Mgl*pow2(Mst1)*pow4(Mst2)) - Al4p*threeLoopFlag*(-(s2t*(3636
        - 1800*lmMsq + 1080*pow2(lmMsq) - 24*lmMst1*(431 - 150*lmMsq + 90*pow2(
        lmMsq)) + 768*pow2(lmMst1) + 48*lmMst2*(363 + 30*lmMsq*(-4 + 3*lmMst1 -
        3*lmMst2) + 451*lmMst2 - lmMst1*(407 + 168*lmMst2) + 45*pow2(lmMsq) +
        16*pow2(lmMst1)) + 7296*pow3(lmMst2) + (24*Dmglst2*(45 + 30*lmMsq*(1 +
        3*lmMsq - 6*lmMst2)*(1 - 2*lmMst1 + 2*lmMst2) + 8*lmMst2*(29 + 8*pow2(
        lmMst1)) + 1068*pow2(lmMst2) - 2*lmMst1*(-83 + 430*lmMst2 + 336*pow2(
        lmMst2)) + 608*pow3(lmMst2)) - (384*(1 + lmMst2)*(4*Dmglst2*lmMst2 +
        Mgl + lmMst2*Mgl)*pow4(Mst2))/pow4(Mst1) + (8*pow2(Mst1)*(80*Dmglst2*(
        lmMst1 - lmMst2)*(14 - 15*lmMsq + 15*lmMst2)*pow2(Mst2) + 10*(lmMst1 -
        lmMst2)*(43 - 30*lmMsq + 30*lmMst2)*Mgl*pow2(Mst2) + 6*Dmglst2*pow2(
        Msq)*(96 + lmMst2*(173 + 30*lmMsq + 90*pow2(lmMsq) + 288*pow2(lmMst1))
        + (526 - 180*lmMsq)*pow2(lmMst2) - lmMst1*(13 + 30*lmMsq*(1 - 6*lmMst2)
        + 526*lmMst2 + 90*pow2(lmMsq) + 848*pow2(lmMst2)) + 560*pow3(lmMst2)) +
        3*Mgl*pow2(Msq)*(64 + 15*lmMst2*(33 - 10*lmMsq + 6*pow2(lmMsq)) + 288*(
        1 + lmMst2)*pow2(lmMst1) + 6*(173 - 30*lmMsq)*pow2(lmMst2) - lmMst1*(
        431 + 1326*lmMst2 - 30*lmMsq*(5 + 6*lmMst2) + 90*pow2(lmMsq) + 848*
        pow2(lmMst2)) + 560*pow3(lmMst2))) - 40*(8*Dmglst2*(14 - 15*lmMsq + 15*
        lmMst2) + (43 - 30*lmMsq + 30*lmMst2)*Mgl)*(pow4(Mst2) - 2*(lmMst1 -
        lmMst2)*(pow4(Mst1) + pow4(Mst2))))/(pow2(Msq)*pow2(Mst2)) + ((2*pow4(
        Mst1)*((24*Dmglst2*(96 + lmMst2*(285 + 30*lmMsq + 90*pow2(lmMsq) + 544*
        pow2(lmMst1)) + (590 - 180*lmMsq)*pow2(lmMst2) - lmMst1*(109 + 30*
        lmMsq*(1 - 6*lmMst2) + 590*lmMst2 + 90*pow2(lmMsq) + 1360*pow2(lmMst2))
        + 816*pow3(lmMst2)) + 12*Mgl*(80 + lmMst2*(479 - 150*lmMsq + 90*pow2(
        lmMsq)) + 544*(1 + lmMst2)*pow2(lmMst1) + 2*(631 - 90*lmMsq)*pow2(
        lmMst2) - lmMst1*(399 + 1806*lmMst2 - 30*lmMsq*(5 + 6*lmMst2) + 90*
        pow2(lmMsq) + 1360*pow2(lmMst2)) + 816*pow3(lmMst2)))*pow4(Msq) + (
        lmMst1 - lmMst2)*(90*Dmglst2*(13 - 28*lmMsq + 28*lmMst2) + 5*(67 - 84*
        lmMsq + 84*lmMst2)*Mgl)*pow4(Mst2)))/pow4(Mst2) + pow2(Mst2)*(-5*(18*
        Dmglst2*(13 - 28*lmMsq + 28*lmMst2) + (67 - 84*lmMsq + 84*lmMst2)*Mgl)*
        (pow2(Mst2) - 2*(lmMst1 - lmMst2)*(pow2(Mst1) + pow2(Mst2))) - (40*(8*
        Dmglst2*(14 - 15*lmMsq + 15*lmMst2) + (43 - 30*lmMsq + 30*lmMst2)*Mgl)*
        pow2(Msq)*pow2(Mst2) + (-12*Mgl*(335 + 654*lmMst2 + 64*lmMst1*(1 +
        lmMst2) - 30*lmMsq*(5 + 6*lmMst2) + 90*pow2(lmMsq) + 272*pow2(lmMst2))
        - 24*Dmglst2*(-115 + (366 + 64*lmMst1)*lmMst2 - 30*lmMsq*(-1 + 6*
        lmMst2) + 90*pow2(lmMsq) + 336*pow2(lmMst2)))*pow4(Msq) + (90*Dmglst2*(
        13 - 28*lmMsq + 28*lmMst2) + 5*(67 - 84*lmMsq + 84*lmMst2)*Mgl)*pow4(
        Mst2))/pow2(Mst1)))/pow4(Msq))/Mgl)) + (1536*Mt*(2*Dmglst2*pow2(Mst2)*(
        pow2(Mst1)*pow2(Mst2)*(8 + 7*lmMst2 - 11*pow2(lmMst2) + lmMst1*(-3 + 5*
        lmMst2 + 6*pow2(lmMst2)) - 6*pow3(lmMst2)) + 2*(8 + 13*lmMst2 - 8*pow2(
        lmMst2) + lmMst1*(-5 + 5*lmMst2 + 6*pow2(lmMst2)) - 6*pow3(lmMst2))*
        pow4(Mst1) + (1 - 2*lmMst2 - 3*pow2(lmMst2))*pow4(Mst2)) + (1 + lmMst2)
        *Mgl*(2*(1 - 10*lmMst2 + 4*lmMst1*(2 + lmMst2) - 4*pow2(lmMst2))*pow2(
        Mst2)*pow4(Mst1) - 2*(1 + 5*lmMst2 - lmMst1*(3 + 2*lmMst2) + 2*pow2(
        lmMst2))*pow2(Mst1)*pow4(Mst2) + (7 - 32*lmMst2 + 4*lmMst1*(7 + 3*
        lmMst2) - 12*pow2(lmMst2))*pow6(Mst1) - 2*(1 + lmMst2)*pow6(Mst2))))/(
        Mgl*pow2(Mst1)*pow5(Mst2)) + (xDmglst2*pow2(Dmglst2)*(s2t*pow2(Mst1)*(
        160*pow2(Msq)*pow2(Mst2)*(5*(8 - 15*lmMsq)*(pow2(Mst1) + pow2(Mst2))*
        pow4(Mst2) + 30*pow2(lmMst2)*(5*pow2(Mst2)*pow4(Mst1) + 5*pow2(Mst1)*
        pow4(Mst2) + 2*pow6(Mst1)) - 2*lmMst1*(5*(8 - 15*lmMsq + 15*lmMst2)*
        pow2(Mst1)*pow2(Mst2)*(pow2(Mst1) + pow2(Mst2)) + (28 - 30*lmMsq + 30*
        lmMst2)*pow6(Mst1)) + lmMst2*(10*(8 - 15*lmMsq)*pow2(Mst2)*pow4(Mst1) +
        5*(31 - 30*lmMsq)*pow2(Mst1)*pow4(Mst2) + (56 - 60*lmMsq)*pow6(Mst1) +
        75*pow6(Mst2))) + 45*pow4(Mst2)*(7*(5 - 28*lmMsq)*(pow2(Mst1) + pow2(
        Mst2))*pow4(Mst2) + 56*pow2(lmMst2)*(7*pow2(Mst2)*pow4(Mst1) + 7*pow2(
        Mst1)*pow4(Mst2) + 2*pow6(Mst1)) - 2*lmMst1*(7*(5 - 28*lmMsq + 28*
        lmMst2)*pow2(Mst1)*pow2(Mst2)*(pow2(Mst1) + pow2(Mst2)) + (26 - 56*
        lmMsq + 56*lmMst2)*pow6(Mst1)) + lmMst2*(14*(5 - 28*lmMsq)*pow2(Mst2)*
        pow4(Mst1) + 14*(19 - 28*lmMsq)*pow2(Mst1)*pow4(Mst2) + (52 - 112*
        lmMsq)*pow6(Mst1) + 196*pow6(Mst2)))) - 4*pow4(Msq)*(-((Mst2*s2t*(3035
        + 5240*lmMst2 - 270*lmMsq*(5 + 3*lmMsq - 6*lmMst2)*(1 - 2*lmMst1 + 2*
        lmMst2) + 192*(2 - 3*lmMst2)*pow2(lmMst1) - 4620*pow2(lmMst2) + 6*
        lmMst1*(-1161 + 458*lmMst2 + 1008*pow2(lmMst2)) - 5472*pow3(lmMst2)) -
        128*Mt*(-53 - 191*lmMst2 + lmMst1*(60 + 18*lmMst2 - 72*pow2(lmMst2)) +
        54*pow2(lmMst2) + 72*pow3(lmMst2)))*pow3(Mst2)*pow4(Mst1)) - 3*(Mst2*
        s2t*(1065 + 64*lmMst1*(2 - 3*lmMst2) - 266*lmMst2 + 90*lmMsq*(-5 + 6*
        lmMst2) - 270*pow2(lmMsq) - 1264*pow2(lmMst2)) + 512*Mt*(2 + lmMst2 -
        3*pow2(lmMst2)))*pow2(Mst1)*pow5(Mst2) - 2*Mst2*(3*Mst2*s2t*(480 - 9*
        lmMst2*(-129 + 50*lmMsq + 30*pow2(lmMsq)) + 288*(2 - 3*lmMst2)*pow2(
        lmMst1) + 10*(-17 + 54*lmMsq)*pow2(lmMst2) + lmMst1*(-1385 - 406*lmMst2
        - 90*lmMsq*(-5 + 6*lmMst2) + 270*pow2(lmMsq) + 2544*pow2(lmMst2)) -
        1680*pow3(lmMst2)) + 64*Mt*(11 + 371*lmMst2 - 42*pow2(lmMst2) + 6*
        lmMst1*(-21 - 5*lmMst2 + 24*pow2(lmMst2)) - 144*pow3(lmMst2)))*pow6(
        Mst1) + s2t*(12*(96 + lmMst2*(285 + 30*lmMsq + 90*pow2(lmMsq) + 544*
        pow2(lmMst1)) + (590 - 180*lmMsq)*pow2(lmMst2) - lmMst1*(109 + 30*
        lmMsq*(1 - 6*lmMst2) + 590*lmMst2 + 90*pow2(lmMsq) + 1360*pow2(lmMst2))
        + 816*pow3(lmMst2))*pow8(Mst1) - 192*(-2 + lmMst2 + 5*pow2(lmMst2))*
        pow8(Mst2)))))/(pow2(Mgl)*pow4(Msq)*pow4(Mst1)*pow4(Mst2))))) + (2*
        Al4p*z2*(-9720*s2t*twoLoopFlag*(Mgl*pow2(Mst2)*(4*Mt*((9*Dmglst2 + Mgl)
        *pow2(Mst1)*pow2(Mst2) + (13*Dmglst2 + Mgl)*pow4(Mst1) + (5*Dmglst2 +
        Mgl)*pow4(Mst2)) - Mst2*s2t*(Mgl*pow4(Mst2) + 4*Dmglst2*(pow2(Mst1)*
        pow2(Mst2) + pow4(Mst1) + pow4(Mst2)))) + 2*xMst*(Dmglst2*Mgl*(34*Mt -
        2*Mst2*s2t) + (82*Mt - 3*Mst2*s2t)*xDmglst2*pow2(Dmglst2) + 2*Mt*pow2(
        Mgl))*pow6(Mst1)) + (pow2(Mst2)*(2*xDmglst2*pow2(Dmglst2)*(-9720*s2t*
        twoLoopFlag*pow2(Mst1)*pow4(Msq)*(-3*(-14*Mt + Mst2*s2t)*pow2(Mst1)*
        pow2(Mst2) + (62*Mt - 3*Mst2*s2t)*pow4(Mst1) + (22*Mt - 3*Mst2*s2t)*
        pow4(Mst2)) + Al4p*Mst2*threeLoopFlag*(4050*Mst2*s2t*pow2(Mst1)*((372*
        Mt - 22*Mst2*s2t)*pow2(Mst2)*pow4(Mst1) + (1252*Mt - 97*Mst2*s2t)*pow2(
        Mst1)*pow4(Mst2) + 2*pow2(Msq)*(2*(576*Mt - 41*Mst2*s2t)*pow2(Mst1)*
        pow2(Mst2) + 4*(120*Mt - 7*Mst2*s2t)*pow4(Mst1) + 13*(48*Mt - 7*Mst2*
        s2t)*pow4(Mst2)) + 2*(358*Mt - 61*Mst2*s2t)*pow6(Mst2)) + pow4(Msq)*(-
        9*pow2(Mst1)*pow2(Mst2)*((1340537 - 496800*lmMsq + 4680*lmMst1 +
        778680*lmMst2)*Mst2*Mt*s2t + 192*(617 + 2040*lmMst2)*pow2(Mt) + (-
        115931 + 48600*lmMsq + 53280*lmMst1 - 148320*lmMst2)*pow2(Mst2)*pow2(
        s2t)) - 6*((3255761 - 1458000*lmMsq - 351000*lmMst1 + 2631960*lmMst2)*
        Mst2*Mt*s2t + 180*(3643 - 120*lmMst1 + 3192*lmMst2)*pow2(Mt) + 5*(-
        47953 + 14580*lmMsq + 25650*lmMst1 - 54162*lmMst2)*pow2(Mst2)*pow2(s2t)
        )*pow4(Mst1) + pow2(s2t)*(4*(2286439 - 72900*lmMsq - 307800*lmMst1 +
        450360*lmMst2)*pow6(Mst1) + 233280*pow6(Mst2))))) - Al4p*threeLoopFlag*
        (12960*s2t*xDR2DRMOD*pow4(Msq)*(Mst2*xDmglst2*pow2(Dmglst2)*((32*(-1 +
        6*lmMst2)*Mt + (49 - 66*lmMst1 + 42*lmMst2)*Mst2*s2t)*pow2(Mst1)*pow3(
        Mst2) + 2*Mst2*(32*(-1 + 6*lmMst2)*Mt + (8 - 33*lmMst1 + 21*lmMst2)*
        Mst2*s2t)*pow4(Mst1) + (-46*lmMst1*s2t + 30*lmMst2*s2t)*pow6(Mst1) +
        33*s2t*pow6(Mst2)) + Mgl*(32*Mt*pow2(Mst1)*(Dmglst2*(1 + 3*lmMst2)*
        pow2(Mst2)*(2*pow2(Mst1) + pow2(Mst2)) + (1 + lmMst2)*Mgl*(2*pow2(Mst1)
        *pow2(Mst2) + 3*pow4(Mst1) + pow4(Mst2))) + Mst2*s2t*(-2*(Dmglst2*(23*
        lmMst1 - 15*lmMst2) + (4 + 13*lmMst1 - 9*lmMst2)*Mgl)*(pow2(Mst1) +
        pow2(Mst2))*pow4(Mst1) + (Dmglst2*(23 - 46*lmMst1 + 30*lmMst2) + (5 -
        26*lmMst1 + 18*lmMst2)*Mgl)*pow2(Mst1)*pow4(Mst2) + (23*Dmglst2 + 13*
        Mgl)*pow6(Mst2)))) + Mgl*(90*Mt*pow2(Mst1)*(6*Mst2*Mt*(2*(Dmglst2*(7286
        - 240*lmMst1 + 6384*lmMst2) + (-353 + 72*lmMst1 + 696*lmMst2)*Mgl)*
        pow2(Mst1) + 32*(24*Dmglst2*(5 + 6*lmMst2) + (53 + 24*lmMst2)*Mgl)*
        pow2(Mst2))*pow4(Msq) + s2t*((2*(412663 - 226800*lmMsq + 396360*lmMst1
        - 119880*lmMst2)*Mgl*pow4(Msq)*pow4(Mst1))/15. + 9*(Dmglst2*(-6720*
        pow2(Msq)*pow2(Mst1) + (18556.466666666667 - 5280*lmMsq + 248*lmMst1 +
        7944*lmMst2)*pow4(Msq) - 3720*pow4(Mst1)) + Mgl*(-960*pow2(Msq)*pow2(
        Mst1) + (3469 - 1440*lmMsq + 664*lmMst1 + 1384*lmMst2)*pow4(Msq) - 360*
        pow4(Mst1)))*pow4(Mst2) + 2*pow2(Msq)*pow2(Mst2)*((2*(Dmglst2*(534391 -
        113400*lmMsq - 23760*lmMst1 + 196560*lmMst2) + 10*(5827 - 2700*lmMsq +
        2988*lmMst1 + 468*lmMst2)*Mgl)*pow2(Msq)*pow2(Mst1))/5. - 2880*(15*
        Dmglst2 + 2*Mgl)*pow4(Mst1) - 2880*(6*Dmglst2 + Mgl)*pow4(Mst2)) - 360*
        ((67*Dmglst2 + 7*Mgl)*Mst2*pow2(Mst1) + (41*Dmglst2 + 5*Mgl)*pow3(Mst2)
        )*pow5(Mst2))) - Mst2*pow2(s2t)*(-3*Mgl*(32400*xMsq*((1 - 2*lmMst2*(-1
        + shiftst1))*pow2(Mst1)*pow2(Mst2) + 2*lmMst1*(-1 + shiftst1)*pow2(
        Mst1)*pow2(Mst2) + shiftst2*pow2(Mst1)*(2*pow2(Mst1) + pow2(Mst2)) - 2*
        lmMst1*(1 - 2*shiftst1 + shiftst2)*pow4(Mst1) + 2*lmMst2*(1 - 2*
        shiftst1 + shiftst2)*pow4(Mst1) + pow4(Mst2) - shiftst1*(2*pow2(Mst1)*
        pow2(Mst2) + 2*pow4(Mst1) + pow4(Mst2)))*pow6(Msq) - pow4(Msq)*(-6*(
        240379 + 10800*lmMsq - 21420*lmMst2 - 1080*(1 + 2*lmMst2)*shiftst3 +
        180*lmMst1*(55 + 12*shiftst3))*pow2(Mst2)*pow4(Mst1) + 180*(59 - 450*
        lmMsq + 1362*lmMst2 + 36*(1 + lmMst2)*shiftst3 - 6*lmMst1*(107 + 6*
        shiftst3))*pow2(Mst1)*pow4(Mst2) - (3647353 + 64800*lmMsq - 52380*
        lmMst2 - 6480*(1 + 3*lmMst2)*shiftst3 + 540*lmMst1*(-31 + 36*shiftst3))
        *pow6(Mst1) + 1080*(49 + 3*shiftst3)*pow6(Mst2)) + 5400*pow2(Msq)*(4*
        pow4(Mst1)*pow4(Mst2) + 4*pow2(Mst2)*pow6(Mst1) + 7*pow2(Mst1)*pow6(
        Mst2)) + 1350*(4*pow4(Mst2)*pow6(Mst1) + 4*pow4(Mst1)*pow6(Mst2) + 9*
        pow2(Mst1)*pow8(Mst2))) + 4*Dmglst2*(pow4(Msq)*(162*(9053 - 900*lmMsq -
        1810*lmMst1 + 3570*lmMst2)*pow2(Mst2)*pow4(Mst1) + 135*(4151 - 1080*
        lmMsq - 1368*lmMst1 + 3480*lmMst2)*pow2(Mst1)*pow4(Mst2) + 2*(2286439 -
        72900*lmMsq - 307800*lmMst1 + 450360*lmMst2)*pow6(Mst1) + 74520*pow6(
        Mst2)) - 8100*pow2(Msq)*(14*pow4(Mst1)*pow4(Mst2) + 14*pow2(Mst2)*pow6(
        Mst1) + 17*pow2(Mst1)*pow6(Mst2)) - 4050*(11*pow4(Mst2)*pow6(Mst1) +
        11*pow4(Mst1)*pow6(Mst2) + 16*pow2(Mst1)*pow8(Mst2))))))))/(pow2(Mst1)*
        pow4(Msq))))/(pow2(Mgl)*pow7(Mst2)))))/29160.;
}

/**
 * @return The diagonal (2, 2) matrix element of the Higgs mass matrix as a double for the hierarchy 'H6g2'
 */
double H6g2::getS2() const {
   return -(oneLoopFlag*((4*Mt*MuSUSY*s2t*(4*(lmMst1 - lmMst2)*pow2(Mt) + ((-2 -
        lmMst1 + lmMst2)*pow2(Mst1) + (2 - lmMst1 + lmMst2)*pow2(Mst2))*pow2(
        s2t)))/Tbeta + 4*pow2(Mt)*pow2(s2t)*(2*(lmMst1 - lmMst2)*(pow2(Mst1) -
        pow2(Mst2)) + pow2(MuSUSY)*(2 - lmMst1 + lmMst2 - (2*(lmMst1 - lmMst2)*
        pow2(Mst1)*(pow2(Mst1) + pow2(Mst2)))/pow4(Mst2))) - (4*pow2(Mt)*pow2(
        MuSUSY)*pow2(s2t)*(2 - lmMst1 + lmMst2 - (2*(lmMst1 - lmMst2)*pow2(
        Mst1)*(pow2(Mst1) + pow2(Mst2)))/pow4(Mst2)))/pow2(Sbeta) + 16*(lmMst1
        + lmMst2 - 2*lmMt)*pow4(Mt) - (-4*pow2(Mst1)*pow2(Mst2) + (2 + lmMst1 -
        lmMst2)*pow4(Mst1) + (2 - lmMst1 + lmMst2)*pow4(Mst2))*pow4(s2t)))/32.
         + (threeLoopFlag*pow2(Al4p)*(Dmglst2*Mst2*(-2*Mgl*(pow4(Mst1)*(4*pow2(
        Mt)*pow2(s2t)*(4*Tbeta*pow2(MuSUSY)*(-1 + pow2(Sbeta))*(338536*z3 +
        11327*z4 + 15897*pow2(z2)) + pow2(Sbeta)*(36*Tbeta*pow2(Mst2)*(-734*z3
        + 1538*z4 + 2307*pow2(z2)) - 3*Mst2*MuSUSY*(-844756*z3 + 10126*z4 +
        105585*pow2(z2)))) + pow2(Sbeta)*(-32*s2t*(9*MuSUSY*(2806*z3 + 1514*z4
        + 2757*pow2(z2)) + Mst2*Tbeta*(-614702*z3 + 25730*z4 + 38595*pow2(z2)))
        *pow3(Mt) + 2*Mt*pow2(Mst2)*(8*MuSUSY*(89947*z3 + 6332*z4 + 9498*pow2(
        z2)) + Mst2*Tbeta*(-565214*z3 + 31142*z4 + 46713*pow2(z2)))*pow3(s2t))
        + Tbeta*pow2(Sbeta)*(64*(-562226*z3 + 5582*z4 + 8373*pow2(z2))*pow4(Mt)
        + (197843*z3 - 16796*z4 - 25194*pow2(z2))*pow4(Mst2)*pow4(s2t))) + 18*
        pow2(Mst1)*(pow2(Mst2)*pow2(Mt)*(16*Tbeta*pow2(Mt)*pow2(Sbeta)*(-24526*
        z3 + 322*z4 + 483*pow2(z2)) + 3*MuSUSY*pow2(s2t)*(Mst2*pow2(Sbeta)*(
        124922*z3 + 1210*z4 - 18273*pow2(z2)) + 8*MuSUSY*Tbeta*(-1 + pow2(
        Sbeta))*(9207*z3 + 185*z4 + 237*pow2(z2))) - 12*Mt*s2t*pow2(Sbeta)*(8*
        MuSUSY*(590*z3 - 4*z4 + 75*pow2(z2)) + Mst2*Tbeta*(-21022*z3 + 922*z4 +
        1383*pow2(z2)))) + 8*Mst2*MuSUSY*(24*Mt*pow2(Sbeta)*(-955*z3 + 26*z4 +
        93*pow2(z2)) - MuSUSY*s2t*Tbeta*(-1 + pow2(Sbeta))*(-59336*z3 - 709*z4
        + 8535*pow2(z2)))*pow3(Mt) + 6*Mt*pow2(s2t)*pow2(Sbeta)*(4*Mt*Tbeta*(
        741*z3 + 100*z4 + 150*pow2(z2)) + MuSUSY*s2t*(20653*z3 + 316*z4 + 474*
        pow2(z2)))*pow4(Mst2) + Tbeta*(-48*pow2(MuSUSY)*(-1 + pow2(Sbeta))*(
        250*z3 - 94*z4 + 21*pow2(z2))*pow4(Mt) + pow2(Sbeta)*(3*Mst2*s2t*(-
        2239*z3 + 54*z4) + 2*Mt*(-6250*z3 + 208*z4 + 1203*pow2(z2)))*pow3(s2t)*
        pow5(Mst2))) + 27*pow2(Mst2)*(4*Mt*pow2(Mst2)*(pow2(Mst2)*pow2(s2t)*
        pow2(Sbeta)*(4*Mt*Tbeta*(439*z3 - 108*z4) + MuSUSY*s2t*(16175*z3 + 424*
        z4 + 474*pow2(z2))) + Mt*(Tbeta*pow2(MuSUSY)*pow2(s2t)*(-1 + pow2(
        Sbeta))*(16175*z3 + 424*z4 + 474*pow2(z2)) + 2*Mt*pow2(Sbeta)*(4*
        MuSUSY*s2t*(-439*z3 + 108*z4) + Mt*Tbeta*(11897*z3 + 56*z4 + 84*pow2(
        z2))))) - 2*s2t*pow2(Mt)*pow2(Sbeta)*(4*Mt*Tbeta*(5742*z3 - 506*z4 +
        105*pow2(z2)) + 3*MuSUSY*s2t*(-37474*z3 - 542*z4 + 5289*pow2(z2)))*
        pow3(Mst2) + 4*Mst2*MuSUSY*(2*Mt*pow2(Sbeta)*(5742*z3 - 506*z4 + 105*
        pow2(z2)) - MuSUSY*s2t*Tbeta*(-1 + pow2(Sbeta))*(-37474*z3 - 542*z4 +
        5289*pow2(z2)))*pow3(Mt) + Tbeta*(-64*(83*z3 - 27*z4)*pow2(MuSUSY)*(-1
        + pow2(Sbeta))*pow4(Mt) + pow2(Sbeta)*(-(Mst2*s2t*(16175*z3 + 424*z4 +
        474*pow2(z2))) + 2*Mt*(-37474*z3 - 542*z4 + 5289*pow2(z2)))*pow3(s2t)*
        pow5(Mst2)))) - Dmglst2*xDmglst2*(pow4(Mst1)*(-8*pow2(Mt)*pow2(s2t)*(-
        4*Tbeta*pow2(MuSUSY)*(-1 + pow2(Sbeta))*(338536*z3 + 11327*z4 + 15897*
        pow2(z2)) + 6*Tbeta*pow2(Mst2)*pow2(Sbeta)*(4404*z3 + 105742*xDR2DRMOD*
        z3 - 9228*z4 + 4838*xDR2DRMOD*z4 + 3*(-4614 + 2419*xDR2DRMOD)*pow2(z2))
        - 3*Mst2*MuSUSY*pow2(Sbeta)*(844756*z3 - 408118*xDR2DRMOD*z3 - 10126*z4
        + 83590*xDR2DRMOD*z4 + 3*(-35195 + 17009*xDR2DRMOD)*pow2(z2))) + 4*Mt*
        pow2(Mst2)*pow2(Sbeta)*(Mst2*Tbeta*(2*(-282607 + 789712*xDR2DRMOD)*z3 +
        2*(15571 - 63244*xDR2DRMOD)*z4 + (46713 - 189732*xDR2DRMOD)*pow2(z2)) +
        4*MuSUSY*((179894 + 614777*xDR2DRMOD)*z3 + 8*(1583 - 898*xDR2DRMOD)*z4
        + (18996 - 28272*xDR2DRMOD)*pow2(z2)))*pow3(s2t) + pow2(Sbeta)*(-64*
        s2t*(9*MuSUSY*(2806*z3 + 1514*z4 + 2757*pow2(z2)) + Mst2*Tbeta*(-
        614702*z3 + 1442476*xDR2DRMOD*z3 + 25730*z4 - 81328*xDR2DRMOD*z4 + (
        38595 - 121992*xDR2DRMOD)*pow2(z2)))*pow3(Mt) + 128*Tbeta*(-562226*z3 +
        5582*z4 + 8373*pow2(z2))*pow4(Mt) + Tbeta*(395686*z3 - 699017*
        xDR2DRMOD*z3 - 33592*z4 + 61508*xDR2DRMOD*z4 + 6*(-8398 + 15377*
        xDR2DRMOD)*pow2(z2))*pow4(Mst2)*pow4(s2t))) + 6*pow2(Mst1)*(-8*pow2(
        Mst2)*pow2(Mt)*(-(Tbeta*pow2(MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta))*(
        409448*z3 + 10819*z4 + 3471*pow2(z2))) + 2*Mt*pow2(Sbeta)*(36*MuSUSY*
        s2t*(3124*z3 + 29*z4 + 165*pow2(z2)) + Mt*Tbeta*(-316834*z3 + 9406*z4 +
        14109*pow2(z2)))) + 2*Mt*s2t*pow2(Sbeta)*(Mst2*s2t*(MuSUSY*s2t*(958451*
        z3 + 16612*z4 + 1590*pow2(z2)) + 9*Mt*Tbeta*(-9125*z3 + 1592*z4 + 2388*
        pow2(z2))) + Mt*(3*MuSUSY*s2t*(629162*z3 + 17194*z4 - 84045*pow2(z2)) +
        4*Mt*Tbeta*(-68146*z3 + 8614*z4 + 12921*pow2(z2))))*pow3(Mst2) + 8*
        Mst2*MuSUSY*(-2*Mt*pow2(Sbeta)*(-118001*z3 + 6920*z4 + 6492*pow2(z2)) -
        MuSUSY*s2t*Tbeta*(-1 + pow2(Sbeta))*(-699292*z3 - 12344*z4 + 88647*
        pow2(z2)))*pow3(Mt) + Tbeta*(-288*pow2(MuSUSY)*(-1 + pow2(Sbeta))*(250*
        z3 - 94*z4 + 21*pow2(z2))*pow4(Mt) + pow2(Sbeta)*(-8*Mt*(-35065*z3 +
        2425*z4 + 2301*pow2(z2)) + Mst2*s2t*(-139555*z3 + 5026*z4 + 5352*pow2(
        z2)))*pow3(s2t)*pow5(Mst2))) + 9*pow2(Mst2)*(-4*pow2(Mst2)*pow2(Mt)*(-(
        Tbeta*pow2(MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta))*(226447*z3 + 8888*z4 +
        4098*pow2(z2))) + 2*Mt*pow2(Sbeta)*(9*MuSUSY*s2t*(19703*z3 - 376*z4 +
        84*pow2(z2)) + Mt*Tbeta*(-225923*z3 + 728*z4 + 1092*pow2(z2)))) + 4*Mt*
        s2t*pow2(Sbeta)*(Mt*(MuSUSY*s2t*(769422*z3 + 7494*z4 - 93249*pow2(z2))
        + 4*Mt*Tbeta*(-55952*z3 + 1742*z4 + 21*pow2(z2))) + Mst2*s2t*(9*Mt*
        Tbeta*(19703*z3 - 376*z4 + 84*pow2(z2)) + MuSUSY*s2t*(226447*z3 + 8888*
        z4 + 4098*pow2(z2))))*pow3(Mst2) + 8*Mst2*MuSUSY*(-2*Mt*pow2(Sbeta)*(-
        55952*z3 + 1742*z4 + 21*pow2(z2)) - MuSUSY*s2t*Tbeta*(-1 + pow2(Sbeta))
        *(-256474*z3 - 2498*z4 + 31083*pow2(z2)))*pow3(Mt) + Tbeta*(-192*(1990*
        z3 - 81*z4)*pow2(MuSUSY)*(-1 + pow2(Sbeta))*pow4(Mt) + pow2(Sbeta)*(-(
        Mst2*s2t*(226447*z3 + 8888*z4 + 4098*pow2(z2))) + 4*Mt*(-256474*z3 -
        2498*z4 + 31083*pow2(z2)))*pow3(s2t)*pow5(Mst2))))) + 3*pow2(Mgl)*(6*
        Mst2*pow2(Mst1)*(pow2(Mt)*(-4*pow2(Mst2)*(2*Tbeta*pow2(Mt)*pow2(Sbeta)*
        (69050*z3 + 418*z4 + 627*pow2(z2)) - Tbeta*pow2(MuSUSY)*pow2(s2t)*(-1 +
        pow2(Sbeta))*((-56978 - 216*lmMst1 + 216*lmMst2)*z3 + 1136*z4 + 813*
        pow2(z2)) + Mt*MuSUSY*s2t*pow2(Sbeta)*(7078*z3 + 2990*z4 + 4485*pow2(
        z2))) - 2*s2t*pow2(Sbeta)*(3*MuSUSY*s2t*(14114*z3 + 3010*z4 - 4557*
        pow2(z2)) + 4*Mt*Tbeta*(-48454*z3 + 754*z4 + 1131*pow2(z2)))*pow3(Mst2)
        ) - 16*Mst2*MuSUSY*(-4*Mt*pow2(Sbeta)*(-8219*z3 + 326*z4 + 165*pow2(z2)
        ) - MuSUSY*s2t*Tbeta*(-1 + pow2(Sbeta))*(-11396*z3 - 1135*z4 + 2388*
        pow2(z2)))*pow3(Mt) + 4*Mt*pow2(s2t)*pow2(Sbeta)*(4*Mt*Tbeta*((-940 +
        54*lmMst1 - 54*lmMst2)*z3 + 214*z4 + 483*pow2(z2)) + MuSUSY*s2t*((-
        39761 - 108*lmMst1 + 108*lmMst2)*z3 + 1046*z4 + 1245*pow2(z2)))*pow4(
        Mst2) + Tbeta*(-144*pow2(MuSUSY)*(-1 + pow2(Sbeta))*(218*z3 + 50*z4 +
        21*pow2(z2))*pow4(Mt) + pow2(Sbeta)*(Mst2*s2t*(22544*z3 - 956*z4 -
        1677*pow2(z2)) + 4*Mt*(-8678*z3 + 740*z4 + 219*pow2(z2)))*pow3(s2t)*
        pow5(Mst2))) - pow4(Mst1)*(-4*Mst2*pow2(Mt)*pow2(s2t)*(Tbeta*pow2(
        MuSUSY)*(-1 + pow2(Sbeta))*(-8*(76813 + 162*lmMst1 - 162*lmMst2)*z3 +
        15686*z4 + 18183*pow2(z2)) + Mst2*pow2(Sbeta)*(6*MuSUSY*(1144*z3 -
        4954*z4 + 6177*pow2(z2)) + Mst2*Tbeta*(-66166*z3 + 13834*z4 + 20751*
        pow2(z2)))) + 16*s2t*(4*Tbeta*pow2(Mst2)*pow2(Sbeta)*(-57758*z3 + 1370*
        z4 + 2055*pow2(z2)) + Mst2*MuSUSY*pow2(Sbeta)*(-22466*z3 + 11402*z4 +
        17103*pow2(z2)) - Tbeta*pow2(MuSUSY)*(-1 + pow2(Sbeta))*(-67232*z3 -
        11764*z4 + 20505*pow2(z2)))*pow3(Mt) + pow2(Sbeta)*(-4*Mt*(Mst2*Tbeta*(
        -44630*z3 + 878*z4 + 1317*pow2(z2)) + MuSUSY*(-272636*z3 + 8870*z4 +
        13305*pow2(z2)))*pow3(Mst2)*pow3(s2t) + 64*(MuSUSY*(107072*z3 - 3326*z4
        - 3045*pow2(z2)) + 4*Mst2*Tbeta*(24241*z3 + 146*z4 + 219*pow2(z2)))*
        pow4(Mt) + Tbeta*((-34070 + 648*lmMst1 - 648*lmMst2)*z3 + 2594*z4 +
        5835*pow2(z2))*pow4(s2t)*pow5(Mst2))) + 18*pow3(Mst2)*(12*pow2(Mst2)*
        pow2(Mt)*(-(Mt*MuSUSY*s2t*pow2(Sbeta)*((1622 - 48*lmMst1 + 48*lmMst2)*
        z3 + 142*z4 + 69*pow2(z2))) + Tbeta*(-2*pow2(Mt)*pow2(Sbeta)*((3902 -
        24*(lmMst1 + lmMst2) + 48*lmMt)*z3 + 14*z4 + 21*pow2(z2)) + pow2(
        MuSUSY)*pow2(s2t)*(1 - pow2(Sbeta))*((1913 + 12*lmMst1 - 12*lmMst2)*z3
        - 10*z4 + 48*pow2(z2)))) + 6*s2t*pow2(Mt)*pow2(Sbeta)*(-4*Mt*Tbeta*(-
        1922*z3 + 206*z4 + 21*pow2(z2)) + 5*MuSUSY*s2t*(-2098*z3 - 102*z4 +
        333*pow2(z2)))*pow3(Mst2) - 4*Mst2*MuSUSY*(-6*Mt*pow2(Sbeta)*(-1922*z3
        + 206*z4 + 21*pow2(z2)) - 5*MuSUSY*s2t*Tbeta*(-1 + pow2(Sbeta))*(-2098*
        z3 - 102*z4 + 333*pow2(z2)))*pow3(Mt) - 6*Mt*pow2(s2t)*pow2(Sbeta)*(2*
        MuSUSY*s2t*((1913 + 12*lmMst1 - 12*lmMst2)*z3 - 10*z4 + 48*pow2(z2)) -
        Mt*Tbeta*((1622 - 48*lmMst1 + 48*lmMst2)*z3 + 142*z4 + 69*pow2(z2)))*
        pow4(Mst2) + Tbeta*(-32*(230*z3 + 27*z4)*pow2(MuSUSY)*(-1 + pow2(Sbeta)
        )*pow4(Mt) + 10*Mt*pow2(Sbeta)*(2098*z3 + 102*z4 - 333*pow2(z2))*pow3(
        s2t)*pow5(Mst2) + 3*pow2(Sbeta)*((1913 + 12*lmMst1 - 12*lmMst2)*z3 -
        10*z4 + 48*pow2(z2))*pow4(s2t)*pow6(Mst2))))))/(5832.*Tbeta*pow2(Mgl)*
        pow2(Sbeta)*pow5(Mst2)) + (4*Al4p*twoLoopFlag*xMst*z2*pow2(Mt)*(
        xDmglst2*pow2(Dmglst2)*(16*(7*MuSUSY + 15*Mst2*Tbeta)*pow2(Mt)*pow2(
        Sbeta) + 3*Mst2*MuSUSY*pow2(s2t)*(-(MuSUSY*Tbeta*(-1 + pow2(Sbeta))) +
        10*Mst2*pow2(Sbeta)) + 2*Mt*s2t*Tbeta*(41*pow2(MuSUSY)*(-1 + pow2(
        Sbeta)) - 40*pow2(Mst2)*pow2(Sbeta))) - 2*Dmglst2*Mgl*(16*(2*MuSUSY +
        3*Mst2*Tbeta)*pow2(Mt)*pow2(Sbeta) + Mst2*MuSUSY*pow2(s2t)*(MuSUSY*
        Tbeta*(-1 + pow2(Sbeta)) - 6*Mst2*pow2(Sbeta)) + Mt*s2t*Tbeta*(-17*
        pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 20*pow2(Mst2)*pow2(Sbeta))) + 2*Mt*
        pow2(Mgl)*(8*Mt*(2*MuSUSY + Mst2*Tbeta)*pow2(Sbeta) + s2t*Tbeta*(pow2(
        MuSUSY)*(-1 + pow2(Sbeta)) - 4*pow2(Mst2)*pow2(Sbeta))))*pow6(Mst1))/(
        3.*Tbeta*pow2(Mgl)*pow2(Sbeta)*pow7(Mst2)) - (xMst*(Mgl*Mst2*Tbeta*(-
        27*(lmMst1 - lmMst2)*Mgl*oneLoopFlag*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*(-
        1 + pow2(Sbeta)) - 32*Al4p*twoLoopFlag*(-(Mgl*(65 - 159*lmMst2 + 6*
        lmMst1*(25 + 6*lmMst2 - 6*lmMt) + 9*lmMt + 36*lmMst2*lmMt - 36*pow2(
        lmMst2))) + 18*Dmglst2*(4 - 61*lmMst2 + 6*lmMst1*(9 + 2*lmMst2 - 2*
        lmMt) + 7*lmMt + 12*lmMst2*lmMt - 12*pow2(lmMst2)))*pow2(Sbeta)*pow4(
        Mt)) - Al4p*twoLoopFlag*(144*(lmMst1 - lmMst2)*Mst2*Tbeta*xDR2DRMOD*(2*
        Dmglst2*lmMst2*Mgl + (-2 + 3*lmMst2)*xDmglst2*pow2(Dmglst2) + (1 +
        lmMst2)*pow2(Mgl))*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*(1 - pow2(Sbeta)) +
        Mgl*(8*s2t*Tbeta*(Dmglst2*Mst2*s2t*(-49 + 84*lmMst2 - 12*lmMst1*(7 + 3*
        lmMst2) + 36*pow2(lmMst2)) + Mgl*Mst2*s2t*(-1 + 111*lmMst2 - 3*lmMst1*(
        37 + 72*lmMst2) + 81*pow2(lmMst1) + 135*pow2(lmMst2)) + Dmglst2*Mt*(785
        - 438*lmMst2 - 6*lmMst1*(-85 + 42*lmMst2) + 252*pow2(lmMst2)) + Mgl*Mt*
        (193 + 474*lmMst2 - 6*lmMst1*(67 + 42*lmMst2) + 252*pow2(lmMst2)))*
        pow2(Mt)*pow2(MuSUSY) + 8*Mst2*Tbeta*pow2(Mt)*(Mgl*(1 - 111*lmMst2 + 3*
        lmMst1*(37 + 72*lmMst2) - 81*pow2(lmMst1) - 135*pow2(lmMst2))*pow2(
        MuSUSY) - Dmglst2*(36*(lmMst1 - lmMst2)*pow2(Mst2) + (-49 + 84*lmMst2 -
        12*lmMst1*(7 + 3*lmMst2) + 36*pow2(lmMst2))*pow2(MuSUSY)))*pow2(s2t)*
        pow2(Sbeta) - 8*s2t*Tbeta*((4*Dmglst2*(95 - 447*lmMst2 + lmMst1*(411 +
        90*lmMst2 - 90*lmMt) + 36*lmMt + 90*lmMst2*lmMt - 90*pow2(lmMst2)) - 4*
        Mgl*(49 - 75*lmMst2 + 3*lmMst1*(25 + 6*lmMst2 - 6*lmMt) + 18*lmMst2*
        lmMt - 18*pow2(lmMst2)))*pow2(Mst2) + (Dmglst2*(785 + 6*lmMst1*(85 -
        42*lmMst2) - 438*lmMst2 + 252*pow2(lmMst2)) + Mgl*(193 + 474*lmMst2 -
        6*lmMst1*(67 + 42*lmMst2) + 252*pow2(lmMst2)))*pow2(MuSUSY))*pow2(
        Sbeta)*pow3(Mt) + 2*Mt*MuSUSY*pow2(Sbeta)*(8*Dmglst2*(36*(-1 + 3*lmMst1
        - 3*lmMst2)*Mst2*s2t*pow2(Mt) - 3*Mt*(50 + lmMst1*(51 - 18*lmMst2) -
        51*lmMst2 + 18*pow2(lmMst2))*pow2(Mst2)*pow2(s2t) + 4*(7 - 381*lmMst2 +
        lmMst1*(327 + 72*lmMst2 - 81*lmMt) + 54*lmMt + 81*lmMst2*lmMt - 72*
        pow2(lmMst2))*pow3(Mt) + 2*(1 + 3*lmMst1 - 3*lmMst2)*pow3(Mst2)*pow3(
        s2t)) - Mgl*(-144*Mst2*s2t*(-lmMst1 + lmMst2 - 2*lmMst1*lmMst2 + pow2(
        lmMst1) + pow2(lmMst2))*pow2(Mt) + 24*Mt*(1 + 33*lmMst2 - 3*lmMst1*(11
        + 6*lmMst2) + 18*pow2(lmMst2))*pow2(Mst2)*pow2(s2t) + 32*(83 - 96*
        lmMst2 + 3*lmMst1*(29 + 12*lmMst2 - 9*lmMt) + 9*lmMt + 27*lmMst2*lmMt -
        36*pow2(lmMst2))*pow3(Mt) + (5 + lmMst1*(6 - 288*lmMst2) - 6*lmMst2 +
        144*pow2(lmMst1) + 144*pow2(lmMst2))*pow3(Mst2)*pow3(s2t))) + 4*(
        Dmglst2*(11 + 42*lmMst1 - 42*lmMst2) + (-5 - 6*lmMst1 + 6*lmMst2)*Mgl)*
        Mt*Tbeta*pow2(Sbeta)*pow3(s2t)*pow4(Mst2) + (2*Dmglst2*(5 + 6*lmMst1 -
        6*lmMst2) + (-11 - 6*lmMst1 + 6*lmMst2)*Mgl)*Tbeta*pow2(Sbeta)*pow4(
        s2t)*pow5(Mst2)) + xDmglst2*pow2(Dmglst2)*(-12*Mst2*pow2(Mt)*pow2(s2t)*
        (Tbeta*(-85 + lmMst1*(60 - 36*lmMst2) - 60*lmMst2 + 36*pow2(lmMst2))*
        pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 2*Mst2*(6*(4 - 9*lmMst1 + 9*lmMst2)*
        Mst2*Tbeta + MuSUSY*(355 + 6*lmMst2 - 6*lmMst1*(1 + 6*lmMst2) + 36*
        pow2(lmMst2)))*pow2(Sbeta)) + 8*s2t*(Tbeta*(2501 + 42*lmMst1*(7 - 6*
        lmMst2) - 222*lmMst2 + 252*pow2(lmMst2))*pow2(MuSUSY)*(1 - pow2(Sbeta))
        + 4*Mst2*((99 - 135*lmMst1 + 135*lmMst2)*MuSUSY + Mst2*Tbeta*(-257 +
        12*lmMst1*(76 + 15*lmMst2 - 15*lmMt) + 162*lmMt + 6*lmMst2*(-179 + 30*
        lmMt) - 180*pow2(lmMst2)))*pow2(Sbeta))*pow3(Mt) + pow2(Sbeta)*(4*Mt*(
        12*(1 - 6*lmMst1 + 6*lmMst2)*MuSUSY + (35 - 102*lmMst1 + 102*lmMst2)*
        Mst2*Tbeta)*pow3(Mst2)*pow3(s2t) + 32*(9*Mst2*Tbeta*(102 + lmMst2*(353
        - 60*lmMt) - 6*lmMst1*(49 + 10*lmMst2 - 10*lmMt) - 59*lmMt + 60*pow2(
        lmMst2)) + MuSUSY*(412 - 6*lmMst1*(176 + 42*lmMst2 - 39*lmMt) + 9*
        lmMst2*(147 - 26*lmMt) - 267*lmMt + 252*pow2(lmMst2)))*pow4(Mt) + 3*(5
        + 6*lmMst1 - 6*lmMst2)*Tbeta*pow4(s2t)*pow5(Mst2)))))*pow6(Mst1))/(108.
        *Tbeta*pow2(Mgl)*pow2(Sbeta)*pow7(Mst2)) - (Al4p*twoLoopFlag*((-72*Mt*
        pow3(s2t)*(-(Mgl*(4*(3 + lmMst1*(2 + lmMst2) - pow2(lmMst2))*pow2(Mst1)
        *pow2(Mst2) + (3 + 2*lmMst1 - 2*lmMst2)*pow4(Mst1) + 4*(-4 + lmMst1 -
        3*lmMst2 + lmMst1*lmMst2 - pow2(lmMst2))*pow4(Mst2))) + Dmglst2*(-4*(1
        + lmMst1*(-2 + lmMst2) + 4*lmMst2 - pow2(lmMst2))*pow2(Mst1)*pow2(Mst2)
        + (1 + 6*lmMst1 - 6*lmMst2)*pow4(Mst1) + 4*(6 + lmMst1 + lmMst2 -
        lmMst1*lmMst2 + pow2(lmMst2))*pow4(Mst2))))/Mst2 + (32*pow4(Mt)*(2*
        Dmglst2*((53 - 18*lmMst1 + 24*lmMst2 - 24*lmMt)*pow2(Mst1)*pow4(Mst2) +
        18*((5 + lmMst2*(11 - 2*lmMt) - 2*lmMst1*(4 + lmMst2 - lmMt) - 3*lmMt +
        2*pow2(lmMst2))*pow2(Mst2)*pow4(Mst1) + (3 - 6*lmMst2*(-5 + lmMt) - 5*
        lmMt + lmMst1*(-25 - 6*lmMst2 + 6*lmMt) + 6*pow2(lmMst2))*pow6(Mst1)) -
        18*lmMst2*pow6(Mst2)) + 9*Mgl*(2*(2 + 2*lmMst1*(3 + lmMst2 - lmMt) +
        lmMt + lmMst2*(-7 + 2*lmMt) - 2*pow2(lmMst2))*pow2(Mst2)*pow4(Mst1) + (
        3 + 2*lmMst1*(3 + lmMst2) + lmMt + lmMst2*(-5 + 4*lmMt) + pow2(lmMst1)
        - pow2(lmMst2) - 6*pow2(lmMt))*pow2(Mst1)*pow4(Mst2) + (9 + lmMst1*(22
        + 6*lmMst2 - 6*lmMt) + 6*lmMst2*(-4 + lmMt) + 2*lmMt - 6*pow2(lmMst2))*
        pow6(Mst1) - 2*(1 + lmMst2)*pow6(Mst2))))/(pow2(Mst1)*pow4(Mst2)) - (
        144*s2t*pow3(Mt)*(Dmglst2*(4*pow2(Mst1)*pow2(Mst2)*((6 + 13*lmMst2 - 4*
        lmMt - 2*lmMst2*lmMt + lmMst1*(-9 - 2*lmMst2 + 2*lmMt) + 2*pow2(lmMst2)
        )*pow2(Mst2) + (-11 + 2*lmMst2 + lmMst1*(-4 + 3*lmMst2) - 3*pow2(
        lmMst2))*pow2(MuSUSY)) - (4*(1 - 29*lmMst2 + lmMst1*(25 + 6*lmMst2 - 6*
        lmMt) + 4*lmMt + 6*lmMst2*lmMt - 6*pow2(lmMst2))*pow2(Mst2) + (65 +
        lmMst1*(34 - 20*lmMst2) - 26*lmMst2 + 20*pow2(lmMst2))*pow2(MuSUSY))*
        pow4(Mst1) - 4*(6 + lmMst1 + lmMst2 - lmMst1*lmMst2 + pow2(lmMst2))*
        pow2(MuSUSY)*pow4(Mst2) + 8*(5 + lmMst1*(-1 + lmMst2) + lmMst2 - pow2(
        lmMst2))*pow6(Mst2)) + Mgl*((4*(5 - 7*lmMst2 + lmMst1*(7 + 2*lmMst2 -
        2*lmMt) + 2*lmMst2*lmMt - 2*pow2(lmMst2))*pow2(Mst2) + (-21 - 38*lmMst2
        + 10*lmMst1*(3 + 2*lmMst2) - 20*pow2(lmMst2))*pow2(MuSUSY))*pow4(Mst1)
        + 4*((-4 + lmMst1 - 3*lmMst2 + lmMst1*lmMst2 - pow2(lmMst2))*pow2(
        MuSUSY)*pow4(Mst2) + pow2(Mst1)*((-5 - 6*lmMst2 + lmMst1*(4 + 3*lmMst2)
        - 3*pow2(lmMst2))*pow2(Mst2)*pow2(MuSUSY) + (4 - 5*lmMst2 + lmMst1*(5 +
        2*lmMst2 - 2*lmMt) + 2*lmMst2*lmMt - 2*pow2(lmMst2))*pow4(Mst2)) + (3 -
        4*lmMst2 + 2*(lmMst1 + lmMst1*lmMst2 + lmMt) - 2*pow2(lmMst2))*pow6(
        Mst2)))))/pow5(Mst2) + (36*Mt*MuSUSY*(16*pow3(Mt)*(Mgl*(13 + 2*lmMst1*(
        6 + 3*lmMst2 - 2*lmMt) + 2*lmMt + 2*lmMst2*(-7 + 2*lmMt) - 6*pow2(
        lmMst2))*pow4(Mst1) + Dmglst2*(11 - 6*lmMst1*lmMst2 - 8*lmMst2*(-5 +
        lmMt) + 8*lmMst1*(-4 + lmMt) - 8*lmMt + 6*pow2(lmMst2))*pow4(Mst1) + 2*
        Dmglst2*((7 + 7*lmMst2 + lmMst1*(-5 + lmMt) - 2*lmMt - lmMst2*lmMt)*
        pow2(Mst1)*pow2(Mst2) + (5 + lmMst1*(-1 + lmMst2) + lmMst2 - pow2(
        lmMst2))*pow4(Mst2)) + 2*Mgl*((4 + lmMst1*(3 + 2*lmMst2 - lmMt) +
        lmMst2*(-4 + lmMt) + lmMt - 2*pow2(lmMst2))*pow2(Mst1)*pow2(Mst2) + (2
        + lmMst1 + (-2 + lmMst1)*lmMst2 + lmMt - pow2(lmMst2))*pow4(Mst2))) +
        6*Mt*pow2(Mst2)*pow2(s2t)*(Dmglst2*(4*(5 + lmMst1*(3 - 2*lmMst2) - 3*
        lmMst2 + 2*pow2(lmMst2))*pow2(Mst1)*pow2(Mst2) + (21 + lmMst1*(18 - 8*
        lmMst2) - 18*lmMst2 + 8*pow2(lmMst2))*pow4(Mst1) + 4*(6 + lmMst1 +
        lmMst2 - lmMst1*lmMst2 + pow2(lmMst2))*pow4(Mst2)) + Mgl*(4*(1 + 3*
        lmMst2 - lmMst1*(3 + 2*lmMst2) + 2*pow2(lmMst2))*pow2(Mst1)*pow2(Mst2)
        + (1 + 14*lmMst2 - 2*lmMst1*(7 + 4*lmMst2) + 8*pow2(lmMst2))*pow4(Mst1)
        + 4*(4 + 3*lmMst2 - lmMst1*(1 + lmMst2) + pow2(lmMst2))*pow4(Mst2))) +
        (pow3(Mst2)*pow3(s2t)*(-4*Dmglst2*(2*(1 + 2*lmMst1 - lmMst2)*pow2(Mst2)
        *pow4(Mst1) + 2*(1 + lmMst1 - lmMst2 + 2*lmMst1*lmMst2 - 2*pow2(lmMst2)
        )*pow2(Mst1)*pow4(Mst2) + (1 + 2*lmMst1 - 2*lmMst2)*pow6(Mst1) - 2*
        lmMst2*pow6(Mst2)) + Mgl*(2*(-16 + 6*lmMst2 - 2*lmMst1*(8 + 5*lmMst2) +
        3*pow2(lmMst1) + 7*pow2(lmMst2))*pow2(Mst2)*pow4(Mst1) - 2*(-12 - 18*
        lmMst2 + 2*lmMst1*(5 + 3*lmMst2) + pow2(lmMst1) - 7*pow2(lmMst2))*pow2(
        Mst1)*pow4(Mst2) + (3 + lmMst1*(2 - 32*lmMst2) - 2*lmMst2 + 16*pow2(
        lmMst1) + 16*pow2(lmMst2))*pow6(Mst1) + 4*(1 + lmMst2)*pow6(Mst2))))/
        pow2(Mst1) + (8*Mst2*s2t*pow2(Mt)*(Dmglst2*(4*(1 - lmMst1 + lmMst2)*
        pow2(Mst2)*pow4(Mst1) + 2*(1 + 2*lmMst2)*pow2(Mst1)*pow4(Mst2) + (4 -
        8*lmMst1 + 8*lmMst2)*pow6(Mst1) - 4*lmMst2*pow6(Mst2)) - Mgl*(-2*
        lmMst1*((1 + 2*lmMst2)*(pow2(Mst1) + pow2(Mst2))*pow4(Mst1) + (3 +
        lmMst2)*pow2(Mst1)*pow4(Mst2)) + pow2(lmMst1)*(2*pow2(Mst2)*pow4(Mst1)
        - pow2(Mst1)*pow4(Mst2) + 2*pow6(Mst1)) + pow2(lmMst2)*(2*pow2(Mst2)*
        pow4(Mst1) + 3*pow2(Mst1)*pow4(Mst2) + 2*pow6(Mst1)) + 2*pow6(Mst2) +
        2*lmMst2*(pow2(Mst2)*pow4(Mst1) + 2*pow2(Mst1)*pow4(Mst2) + pow6(Mst1)
        + pow6(Mst2)))))/pow2(Mst1)))/(Tbeta*pow5(Mst2)) - 36*pow2(Mt)*pow2(
        s2t)*(4*(Dmglst2*(2 - 4*lmMst1) + Mgl*(2*lmMst1*(-2 + lmMst2) + lmMst2*
        (2 + lmMst2) - 3*pow2(lmMst1)))*pow2(Mst1) + 4*(Dmglst2*(2 + 8*lmMst2)
        + Mgl*(2 - 2*lmMst2 + 2*lmMst1*(3 + lmMst2) + pow2(lmMst1) - 3*pow2(
        lmMst2)))*pow2(Mst2) - (16*Dmglst2*(lmMst1 - lmMst2)*pow4(Mst1))/pow2(
        Mst2) + (-8*(2*Dmglst2*lmMst2 + Mgl + lmMst2*Mgl)*pow4(Mst2) + (pow2(
        MuSUSY)*(4*Dmglst2*(2*(2 - 3*lmMst2 + lmMst1*(3 + 2*lmMst2) - 2*pow2(
        lmMst2))*pow2(Mst2)*pow4(Mst1) + 2*(1 + lmMst1 - 2*lmMst2 + 2*lmMst1*
        lmMst2 - 2*pow2(lmMst2))*pow2(Mst1)*pow4(Mst2) + (5 - 8*lmMst2 + 4*
        lmMst1*(2 + lmMst2) - 4*pow2(lmMst2))*pow6(Mst1) - 2*lmMst2*pow6(Mst2))
        + Mgl*(-4*(-1 - 13*lmMst1 + 13*lmMst2 - 8*lmMst1*lmMst2 + pow2(lmMst1)
        + 7*pow2(lmMst2))*pow2(Mst2)*pow4(Mst1) + 2*(-14 + 10*lmMst1 - 20*
        lmMst2 + 6*lmMst1*lmMst2 + pow2(lmMst1) - 7*pow2(lmMst2))*pow2(Mst1)*
        pow4(Mst2) + (1 + 50*lmMst1 - 50*lmMst2 + 64*lmMst1*lmMst2 - 20*pow2(
        lmMst1) - 44*pow2(lmMst2))*pow6(Mst1) - 4*(1 + lmMst2)*pow6(Mst2))))/
        pow4(Mst2))/pow2(Mst1)) + (-9*pow4(s2t)*(4*Dmglst2*(-2*(lmMst1 - 2*
        lmMst1*lmMst2 + 2*pow2(lmMst2))*pow2(Mst2)*pow4(Mst1) - 2*(1 + lmMst1 +
        2*lmMst1*lmMst2 - 2*pow2(lmMst2))*pow2(Mst1)*pow4(Mst2) + (1 + 2*
        lmMst1)*pow6(Mst1) + 2*lmMst2*pow6(Mst2)) + Mgl*(-4*(14 + 6*lmMst2 +
        lmMst1*(3 + 2*lmMst2) - 2*pow2(lmMst1))*pow2(Mst2)*pow4(Mst1) - 2*(-10
        - 16*lmMst2 + 2*lmMst1*(5 + 3*lmMst2) + pow2(lmMst1) - 7*pow2(lmMst2))*
        pow2(Mst1)*pow4(Mst2) + (35 + lmMst1*(34 - 12*lmMst2) - 14*lmMst2 + 10*
        pow2(lmMst1) + 2*pow2(lmMst2))*pow6(Mst1) + 4*(1 + lmMst2)*pow6(Mst2)))
        - (36*s2t*pow2(Mt)*pow2(MuSUSY)*(Mgl*(4*(4*Mt*(5 + 6*lmMst2 - lmMst1*(4
        + 3*lmMst2) + 3*pow2(lmMst2)) + Mst2*s2t*(-1 + 13*lmMst2 - lmMst1*(13 +
        8*lmMst2) + pow2(lmMst1) + 7*pow2(lmMst2)))*pow2(Mst2)*pow4(Mst1) + 2*(
        -(Mst2*s2t*(-14 - 20*lmMst2 + 2*lmMst1*(5 + 3*lmMst2) + pow2(lmMst1) -
        7*pow2(lmMst2))) + 8*Mt*(4 + 3*lmMst2 - lmMst1*(1 + lmMst2) + pow2(
        lmMst2)))*pow2(Mst1)*pow4(Mst2) + (Mst2*s2t*(-1 + 50*lmMst2 - 2*lmMst1*
        (25 + 32*lmMst2) + 20*pow2(lmMst1) + 44*pow2(lmMst2)) + Mt*(84 + 152*
        lmMst2 - 40*lmMst1*(3 + 2*lmMst2) + 80*pow2(lmMst2)))*pow6(Mst1) + 4*(1
        + lmMst2)*s2t*pow7(Mst2)) + 4*Dmglst2*((Mst2*s2t*(-5 + 8*lmMst2 - 4*
        lmMst1*(2 + lmMst2) + 4*pow2(lmMst2)) + Mt*(65 + lmMst1*(34 - 20*
        lmMst2) - 26*lmMst2 + 20*pow2(lmMst2)))*pow6(Mst1) + 2*((Mst2*s2t*(-2 +
        3*lmMst2 - lmMst1*(3 + 2*lmMst2) + 2*pow2(lmMst2)) + Mt*(22 + lmMst1*(8
        - 6*lmMst2) - 4*lmMst2 + 6*pow2(lmMst2)))*pow2(Mst2)*pow4(Mst1) + (-(
        Mst2*s2t*(1 + lmMst1 - 2*lmMst2 + 2*lmMst1*lmMst2 - 2*pow2(lmMst2))) +
        2*Mt*(6 + lmMst1 + lmMst2 - lmMst1*lmMst2 + pow2(lmMst2)))*pow2(Mst1)*
        pow4(Mst2) + lmMst2*s2t*pow7(Mst2)))))/(pow2(Sbeta)*pow5(Mst2)))/pow2(
        Mst1)))/(216.*Mgl) + threeLoopFlag*z2*pow2(Al4p)*((-5*xMsq*pow2(Msq)*(-
        (pow2(s2t)*pow4(Mst1)*(8*(-shiftst1 + lmMst1*(-1 + 2*shiftst1 -
        shiftst2) + shiftst2 + lmMst2*(1 - 2*shiftst1 + shiftst2))*Tbeta*pow2(
        Mt)*pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 4*Mt*(-1 + shiftst2)*(-(MuSUSY*
        s2t) + 2*Mt*Tbeta)*pow2(Mst2)*pow2(Sbeta) + (1 + 2*(lmMst1 - lmMst2)*(-
        1 + shiftst1) + 3*shiftst1 - 4*shiftst2)*Tbeta*pow2(s2t)*pow2(Sbeta)*
        pow4(Mst2))) + (1 - shiftst1)*pow4(Mst2)*(-4*Tbeta*pow2(Mt)*pow2(s2t)*(
        pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 2*pow2(Mst2)*pow2(Sbeta)) + pow2(
        Sbeta)*(16*MuSUSY*s2t*pow3(Mt) - 4*Mt*MuSUSY*pow2(Mst2)*pow3(s2t) + 16*
        Tbeta*pow4(Mt) + Tbeta*pow4(Mst2)*pow4(s2t))) - pow2(Mst1)*pow2(Mst2)*(
        -4*Tbeta*pow2(Mt)*pow2(s2t)*((1 + 2*(lmMst1 - lmMst2)*(-1 + shiftst1) -
        2*shiftst1 + shiftst2)*pow2(MuSUSY)*(1 - pow2(Sbeta)) + 2*(2 + (-1 +
        lmMst1 - lmMst2)*shiftst1 + (-1 - lmMst1 + lmMst2)*shiftst2)*pow2(Mst2)
        *pow2(Sbeta)) + pow2(Sbeta)*(-16*MuSUSY*s2t*(-1 + shiftst2)*pow3(Mt) +
        4*Mt*MuSUSY*(shiftst1 - shiftst2 + lmMst1*(-2 + shiftst1 + shiftst2) -
        lmMst2*(-2 + shiftst1 + shiftst2))*pow2(Mst2)*pow3(s2t) + 16*(-1 +
        shiftst2)*Tbeta*pow4(Mt) + (1 - 4*shiftst1 - 2*(lmMst1 - lmMst2)*(-1 +
        shiftst2) + 3*shiftst2)*Tbeta*pow4(Mst2)*pow4(s2t))) - (-1 + shiftst2)*
        Tbeta*pow2(Mst2)*pow2(Sbeta)*pow4(s2t)*pow6(Mst1)))/(3.*Tbeta*pow2(
        Mst1)*pow2(Sbeta)*pow4(Mst2)) + (xDR2DRMOD*(8*pow2(Mt)*pow2(s2t)*((2*(
        23*Dmglst2 + 13*Mgl)*pow2(pow2(Mst1) - pow2(Mst2)))/(Mgl*pow2(Mst1)) +
        pow2(MuSUSY)*(5 - 26*lmMst1 + 18*lmMst2 + (Dmglst2*(23 - 46*lmMst1 +
        30*lmMst2) + ((23*Dmglst2 + 13*Mgl)*pow2(Mst2))/pow2(Mst1) - (2*(
        Dmglst2*(23*lmMst1 - 15*lmMst2) + (4 + 13*lmMst1 - 9*lmMst2)*Mgl)*pow2(
        Mst1)*(pow2(Mst1) + pow2(Mst2)))/pow4(Mst2))/Mgl)) + (32*pow3(Mt)*(8*
        s2t*((1 + lmMst2)*Mgl*(2*pow2(Mst2) + pow2(MuSUSY))*(2*pow2(Mst1)*pow2(
        Mst2) + 3*pow4(Mst1) + pow4(Mst2)) - Dmglst2*pow2(Mst2)*(-((1 + 3*
        lmMst2)*pow2(Mst2)*pow2(MuSUSY)) + pow2(Mst1)*(4*(3 + lmMst2)*pow2(
        Mst2) - 2*(1 + 3*lmMst2)*pow2(MuSUSY)) + 6*(5 + 3*lmMst2)*pow4(Mst1) -
        2*(-1 + lmMst2)*pow4(Mst2))) + (Mst2*Mt*(Dmglst2*(96*(2 + lmMst2)*pow2(
        Mst2)*pow4(Mst1) + 9*pow2(Mst1)*pow4(Mst2) + 192*(3 + 2*lmMst2)*pow6(
        Mst1) - 23*pow6(Mst2)) - Mgl*(48*(1 + lmMst2)*(2*pow2(Mst1) + pow2(
        Mst2))*pow4(Mst1) + (29 + 16*lmMst2)*pow2(Mst1)*pow4(Mst2) + 13*pow6(
        Mst2))))/pow2(Mst1)))/(Mgl*pow5(Mst2)) + (9*xDmglst2*pow2(Dmglst2)*(Mt*
        pow3(s2t)*((128*(1 - 6*lmMst2)*pow3(Mst2))/9. + ((2345.557475994513 -
        120*lmMst1 + 120*lmMst2)*pow4(Mst1))/Mst2) + (pow2(Mt)*pow2(s2t)*(528*
        pow2(Mst1) - 1056*pow2(Mst2) + 8*(49 - 66*lmMst1 + 42*lmMst2)*pow2(
        MuSUSY) + (16*(8 - 33*lmMst1 + 21*lmMst2)*pow2(Mst1)*pow2(MuSUSY))/
        pow2(Mst2) + (264*pow2(Mst2)*(2*pow2(Mst2) + pow2(MuSUSY)))/pow2(Mst1)
        + ((41.00423280423281 + 3956*lmMst1 - 7796*lmMst2 + (540*pow2(MuSUSY))/
        pow2(Msq))*pow4(Mst1))/pow2(Mst2) + pow2(MuSUSY)*pow4(Mst1)*(375/pow4(
        Msq) + (-368*lmMst1 + 240*lmMst2)/pow4(Mst2))))/9. - (Mt*MuSUSY*((-352*
        s2t*(pow2(Mst1) - pow2(Mst2))*pow2(Mt))/(3.*pow2(Mst1)) + (64*pow3(Mt)*
        (8*lmMst2*(pow2(Mst1) + pow2(Mst2))*pow4(Msq) + 30*pow2(Msq)*pow4(Mst1)
        + pow2(Mst2)*(-16*pow4(Msq) + 45*pow4(Mst1))))/(9.*pow3(Mst2)*pow4(Msq)
        ) + Mt*pow2(s2t)*((128*pow2(Mst1)*(1 - 6*lmMst2 + (15*pow2(Mst1))/pow2(
        Msq)))/(3.*Mst2) + ((7954.436625514403 + 480*lmMsq - (344*lmMst1)/3. -
        (2728*lmMst2)/3.)*pow4(Mst1))/pow3(Mst2) + Mst2*(42.666666666666664 -
        256*lmMst2 + (360*pow4(Mst1))/pow4(Msq))) + pow3(s2t)*((88*pow2(Mst1))/
        3. + (16*(-8 + 33*lmMst1 - 21*lmMst2)*pow2(Mst2))/9. - ((
        1275.4644718792867 - 150*lmMst1 + 150*lmMst2)*pow4(Mst1))/pow2(Mst2) -
        (88*pow4(Mst2))/(3.*pow2(Mst1)))))/Tbeta - (8*s2t*pow3(Mt)*(42525*(100*
        pow2(Msq) + 59*pow2(Mst2))*pow2(MuSUSY)*pow4(Mst1) + pow4(Msq)*(-
        181440*pow2(Mst1)*(2*pow2(Mst2) + (-1 + 6*lmMst2)*pow2(MuSUSY)) + (
        91451239 + 2041200*lmMsq - 5210730*lmMst1 + 6316380*lmMst2 - 561330*
        lmMt)*pow4(Mst1) - 90720*((-1 + 6*lmMst2)*pow2(Mst2)*pow2(MuSUSY) + 2*(
        -2 + lmMst2)*pow4(Mst2)))))/(25515.*pow3(Mst2)*pow4(Msq)) + (16*pow4(
        Mt)*(-96*(4 + lmMst2)*pow2(Mst2)*pow4(Mst1) - 34*pow2(Mst1)*pow4(Mst2)
        + (1152 + 768*lmMst2 + (45*pow2(Mst2))/pow2(Msq))*pow6(Mst1) - 66*pow6(
        Mst2)))/(9.*pow2(Mst1)*pow4(Mst2)) - (pow4(s2t)*(2*(66*lmMst1 - 7*(7 +
        6*lmMst2))*pow2(Mst1)*pow2(Mst2) - 2*(17 + 66*lmMst1 - 42*lmMst2)*pow4(
        Mst2) + pow4(Mst1)*(2498.3904320987654 - 360*lmMst1 + 360*lmMst2 + (15*
        pow2(Mst2))/(2.*pow2(Msq)) + (75*pow4(Mst2))/(4.*pow4(Msq))) + (66*
        pow6(Mst2))/pow2(Mst1)))/9. - (s2t*pow2(Mt)*pow2(MuSUSY)*(60*Mst2*(-
        200*Mt + 9*Mst2*s2t)*pow2(Msq)*pow6(Mst1) + 15*(-472*Mt + 25*Mst2*s2t)*
        pow3(Mst2)*pow6(Mst1) + 8*pow4(Msq)*((32*(-1 + 6*lmMst2)*Mt + (49 - 66*
        lmMst1 + 42*lmMst2)*Mst2*s2t)*pow2(Mst1)*pow3(Mst2) + 2*Mst2*(32*(-1 +
        6*lmMst2)*Mt + (8 - 33*lmMst1 + 21*lmMst2)*Mst2*s2t)*pow4(Mst1) + (-46*
        lmMst1*s2t + 30*lmMst2*s2t)*pow6(Mst1) + 33*s2t*pow6(Mst2))))/(9.*pow2(
        Mst1)*pow2(Sbeta)*pow4(Msq)*pow4(Mst2))))/pow2(Mgl) + (-128*(Dmglst2 +
        3*Dmglst2*lmMst2 + Mgl + lmMst2*Mgl)*Mt*pow3(Mst2)*pow3(s2t) + (-2*(
        pow2(Mst1) - pow2(Mst2))*(Dmglst2*(2*(23*lmMst1 - 15*lmMst2)*pow2(Mst1)
        *pow2(Mst2) + 23*pow4(Mst1) - 23*pow4(Mst2)) + Mgl*(2*(4 + 13*lmMst1 -
        9*lmMst2)*pow2(Mst1)*pow2(Mst2) + 13*pow4(Mst1) - 13*pow4(Mst2)))*pow4(
        s2t) + (Mt*MuSUSY*((384*(Dmglst2 + 3*Dmglst2*lmMst2 + Mgl + lmMst2*Mgl)
        *Mt*pow2(Mst1)*pow2(Mst2)*pow2(s2t)*(pow2(Mst1)*pow2(Mst2) + pow4(Mst1)
        + pow4(Mst2)) + 512*pow2(Mst1)*pow3(Mt)*(Dmglst2*(7 + lmMst2)*pow2(
        Mst1)*pow2(Mst2) - Dmglst2*(-1 + lmMst2)*pow4(Mst2) - (1 + lmMst2)*Mgl*
        (3*pow2(Mst1)*pow2(Mst2) + 6*pow4(Mst1) + pow4(Mst2))) + 32*(23*Dmglst2
        + 13*Mgl)*s2t*(pow2(Mst1) - pow2(Mst2))*pow2(Mt)*pow5(Mst2) - pow3(s2t)
        *(16*(Dmglst2*(23*lmMst1 - 15*lmMst2) + (4 + 13*lmMst1 - 9*lmMst2)*Mgl)
        *pow2(Mst1)*pow2(Mst2) + 8*(23*Dmglst2 + 13*Mgl)*pow4(Mst1) - 8*(23*
        Dmglst2 + 13*Mgl)*pow4(Mst2))*pow5(Mst2))/Tbeta - (8*Mt*MuSUSY*s2t*(
        Dmglst2*Mst2*((32*(1 + 3*lmMst2)*Mt + (23 - 46*lmMst1 + 30*lmMst2)*
        Mst2*s2t)*pow2(Mst1)*pow3(Mst2) + 2*Mst2*(32*(1 + 3*lmMst2)*Mt + (-23*
        lmMst1 + 15*lmMst2)*Mst2*s2t)*pow4(Mst1) + (-46*lmMst1*s2t + 30*lmMst2*
        s2t)*pow6(Mst1) + 23*s2t*pow6(Mst2)) + Mgl*(2*(32*(1 + lmMst2)*Mt + (-4
        - 13*lmMst1 + 9*lmMst2)*Mst2*s2t)*pow2(Mst2)*pow4(Mst1) + (32*(1 +
        lmMst2)*Mt + (5 - 26*lmMst1 + 18*lmMst2)*Mst2*s2t)*pow2(Mst1)*pow4(
        Mst2) + 2*(48*(1 + lmMst2)*Mt + (-4 - 13*lmMst1 + 9*lmMst2)*Mst2*s2t)*
        pow6(Mst1) + 13*s2t*pow7(Mst2))))/pow2(Sbeta)))/pow5(Mst2))/pow2(Mst1))
        /Mgl))/9.) + pow2(Al4p)*((5*(1 - 2*lmMsq)*threeLoopFlag*xMsq*pow2(Msq)*
        (4*Mt*MuSUSY*s2t*(Mt*MuSUSY*s2t*Tbeta*pow2(Mst1)*(pow2(Mst2) - 2*(
        lmMst1 - lmMst2)*(pow2(Mst1) + pow2(Mst2))) + pow2(Mst2)*pow2(Sbeta)*((
        -1 + shiftst1)*pow2(Mst2)*(-4*pow2(Mt) + pow2(Mst2)*pow2(s2t)) - pow2(
        Mst1)*(-4*(-1 + shiftst2)*pow2(Mt) + (shiftst1 - shiftst2 + lmMst1*(-2
        + shiftst1 + shiftst2) - lmMst2*(-2 + shiftst1 + shiftst2))*pow2(Mst2)*
        pow2(s2t)) - (-1 + shiftst2)*pow2(s2t)*pow4(Mst1))) + Tbeta*(-8*pow2(
        Mst1)*pow2(Mst2)*pow2(Mt)*((1 - lmMst1 + lmMst2)*shiftst1*pow2(MuSUSY)*
        pow2(s2t) + 2*shiftst2*pow2(Mt)*pow2(Sbeta)) + 4*pow2(Mt)*pow2(MuSUSY)*
        pow2(s2t)*pow4(Mst2) + pow2(Sbeta)*(-4*pow2(Mt)*pow2(s2t)*(2*(pow2(
        Mst2) + (-lmMst1 + lmMst2)*pow2(MuSUSY))*pow4(Mst1) + pow2(Mst1)*((1 -
        2*lmMst1 + 2*lmMst2)*pow2(Mst2)*pow2(MuSUSY) - 4*pow4(Mst2)) + (2*pow2(
        Mst2) + pow2(MuSUSY))*pow4(Mst2)) + 16*pow2(Mst2)*(pow2(Mst1) + pow2(
        Mst2))*pow4(Mt) + pow2(Mst1)*pow2(Mst2)*((-1 + 2*lmMst1 - 2*lmMst2)*
        pow2(Mst1)*pow2(Mst2) + pow4(Mst1) + (-1 - 2*lmMst1 + 2*lmMst2)*pow4(
        Mst2))*pow4(s2t)) - shiftst1*(4*pow2(Mst1)*pow2(Mst2)*pow2(s2t)*pow2(
        Sbeta)*(2*(1 - lmMst1 + lmMst2)*pow2(Mt)*(pow2(Mst2) - pow2(MuSUSY)) -
        pow2(s2t)*pow4(Mst2)) + pow4(Mst1)*(8*(-1 + 2*lmMst1 - 2*lmMst2)*pow2(
        Mt)*pow2(MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta)) + (3 + 2*lmMst1 - 2*
        lmMst2)*pow2(Sbeta)*pow4(Mst2)*pow4(s2t)) + pow4(Mst2)*(-4*pow2(Mt)*
        pow2(s2t)*(pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 2*pow2(Mst2)*pow2(Sbeta))
        + pow2(Sbeta)*(16*pow4(Mt) + pow4(Mst2)*pow4(s2t)))) - shiftst2*pow2(
        Mst1)*pow2(s2t)*(4*pow2(Mst2)*pow2(Mt)*(pow2(MuSUSY)*(-1 + pow2(Sbeta))
        + 2*(1 + lmMst1 - lmMst2)*pow2(Mst2)*pow2(Sbeta)) - 4*pow2(Mst1)*(2*
        pow2(Mt)*((-1 + lmMst1 - lmMst2)*pow2(MuSUSY)*(-1 + pow2(Sbeta)) +
        pow2(Mst2)*pow2(Sbeta)) + pow2(s2t)*pow2(Sbeta)*pow4(Mst2)) + pow2(s2t)
        *pow2(Sbeta)*(pow2(Mst2)*pow4(Mst1) + (3 - 2*lmMst1 + 2*lmMst2)*pow6(
        Mst2))) + pow2(Sbeta)*pow4(s2t)*pow8(Mst2))))/(6.*Tbeta*pow2(Mst1)*
        pow2(Sbeta)*pow4(Mst2)) + threeLoopFlag*((402.71234567901234 + (520*
        lmMsq)/27. + (220*pow2(lmMsq))/9. - (10*lmMst1*(85 - 10*lmMsq + 6*pow2(
        lmMsq)))/9. - (392*pow2(lmMst1))/9. + (4*lmMst2*(722 - 1449*lmMst1 +
        270*lmMsq*lmMst1 - 135*pow2(lmMsq) + 63*pow2(lmMst1)))/81. + (2*lmMt*(
        1417 - 688*lmMst2 + 90*lmMsq*(-9 + 4*lmMst2) + 12*lmMst1*(27 + 32*
        lmMst2) + 18*pow2(lmMst1) - 246*pow2(lmMst2)))/27. - (4*(-803 + 90*
        lmMsq + 189*lmMst1)*pow2(lmMst2))/27. - (4*(-12 + 10*lmMsq + lmMst1 +
        45*lmMst2)*pow2(lmMt))/3. + (40*pow3(lmMsq))/9. - (68*pow3(lmMst1))/9.
         + (4*(7*Mgl*(22483 + 2250*lmMst2 + 30*lmMst1*(-91 + 150*lmMst2) - 750*
        lmMt + 4500*lmMst2*lmMt - 30*lmMsq*(-41 + 60*lmMst1 + 150*lmMt) + 3150*
        pow2(lmMsq) - 1350*pow2(lmMst1) - 4500*pow2(lmMst2))*pow2(Mst1) + 5*
        Dmglst2*(-6300*(1 + 2*lmMst1 - 4*lmMst2 + 2*lmMt)*pow2(Mst1) + pow2(
        Msq)*(496858 - 397915*lmMst2 + 219975*lmMt + 109200*lmMst2*lmMt - 1050*
        lmMsq*(23 + 18*lmMst1 + 6*lmMst2 + 30*lmMt) + 28350*pow2(lmMsq) - 4410*
        pow2(lmMst1) - 221970*pow2(lmMst2) + 630*lmMst1*(269 + 254*lmMst2 - 64*
        lmMt + 64*pow2(lmMst2)) + 63000*pow2(lmMt) - 40320*pow3(lmMst2)))))/(
        14175.*Mgl*pow2(Msq)) + (148*pow3(lmMst2))/9. + (184*pow3(lmMt))/3. +
        pow4(Mst1)*((16.98404168016413 + (16403*lmMsq)/2205. - (
        4.383446712018141 + (2*lmMsq)/7.)*lmMst1 + (10*lmMst1*lmMst2)/3. - (5*(
        11 + 12*lmMsq - 12*lmMst2)*lmMt)/18. - (20*Dmglst2*(lmMst1 - 2*lmMst2 +
        lmMt))/(3.*Mgl) + (38*pow2(lmMsq))/21. - (32*pow2(lmMst1))/21. - (10*
        pow2(lmMst2))/3.)/pow4(Msq) + (1405.4626938375586 + 140*lmMsq + (
        514.9222171831696 - 520*lmMsq)*lmMst1 - (137992*pow2(lmMst1))/945. -
        lmMst2*(963.0333282942806 - (550874*lmMst1)/945. + 40*lmMsq*(-16 + 3*
        lmMst1) + (376*pow2(lmMst1))/9.) - (579.1343915343915 - 120*lmMsq - (
        656*lmMst1)/3.)*pow2(lmMst2) + (lmMt*(3425 + 1080*lmMsq*(-1 + lmMst1 -
        lmMst2) + 2618*lmMst2 - 2*lmMst1*(573 + 1112*lmMst2) + 312*pow2(lmMst1)
        + 1912*pow2(lmMst2)))/9. + (64*(-1 + 3*lmMst1 - 3*lmMst2)*pow2(lmMt))/
        3. + (64*pow3(lmMst1))/27. - (Dmglst2*(3879.568149934479 + 1840*lmMsq +
        (581.2930880994373 - 2320*lmMsq)*lmMst1 - (126272*pow2(lmMst1))/189. -
        (4*lmMst2*(57911521 - 22348620*lmMst1 + 2381400*lmMsq*(-19 + 3*lmMst1)
        + 2487240*pow2(lmMst1)))/59535. + (16*(-15893 + 5670*lmMsq + 4284*
        lmMst1)*pow2(lmMst2))/189. + (16*lmMt*(2105 + 405*lmMsq*(-3 + 2*lmMst1
        - 2*lmMst2) + 2247*lmMst2 - 3*lmMst1*(341 + 364*lmMst2) + 234*pow2(
        lmMst1) + 858*pow2(lmMst2)))/27. + (128*(-5 + 6*lmMst1 - 6*lmMst2)*
        pow2(lmMt))/3. + (256*pow3(lmMst1))/27. - (5536*pow3(lmMst2))/27.))/Mgl
        - (4840*pow3(lmMst2))/27.)/pow4(Mst2)) + ((2*(2*Dmglst2*(-243 + 2*(55 +
        32*lmMst1)*lmMst2 - 30*lmMsq*(-1 + 6*lmMst2) + 90*pow2(lmMsq) + 208*
        pow2(lmMst2)) + Mgl*(277 + 514*lmMst2 + 64*lmMst1*(1 + lmMst2) - 30*
        lmMsq*(5 + 6*lmMst2) + 90*pow2(lmMsq) + 208*pow2(lmMst2)))*pow2(Mst2))/
        (9.*pow2(Mst1)) + (pow2(Mst1)*(4*Dmglst2*(23048051 + 10812060*lmMst2 +
        12600*(1213 + 285*lmMst2 - 225*lmMt)*pow2(lmMst1) - 58646700*pow2(
        lmMst2) - 2835000*lmMsq*(27 + lmMst2*(19 - 2*lmMt) - 2*lmMst1*(6 +
        lmMst2 - lmMt) - 7*lmMt + 2*pow2(lmMst2)) + 31500*lmMt*(-622 - 389*
        lmMst2 + 6*pow2(lmMst2)) + 1260*lmMst1*(58569 - 1075*lmMt + 5*lmMst2*(
        6883 + 420*lmMt) + 9300*pow2(lmMst2) - 7200*pow2(lmMt)) + 4536000*(3 +
        2*lmMst2)*pow2(lmMt) - 252000*pow3(lmMst1) - 15057000*pow3(lmMst2)) +
        35*Mgl*(6906919 - 4511760*lmMst2 - 1800*(373 + 114*lmMst2 - 90*lmMt)*
        pow2(lmMst1) - 352800*pow2(lmMst2) + 162000*lmMsq*(12 + lmMst2*(15 - 2*
        lmMt) - 2*lmMst1*(5 + lmMst2 - lmMt) - 5*lmMt + 2*pow2(lmMst2)) +
        16200*lmMt*(194 + 71*lmMst2 + 42*pow2(lmMst2)) - 259200*(1 + 2*lmMst2)*
        pow2(lmMt) + 120*lmMst1*(3938 + 1215*lmMt - 15*lmMst2*(7 + 468*lmMt) +
        5940*pow2(lmMst2) + 4320*pow2(lmMt)) + 14400*pow3(lmMst1) - 522000*
        pow3(lmMst2))))/(425250.*pow2(Mst2)) + (pow2(Mst2)*(171500*Dmglst2*(6*(
        83 - 21*lmMst1 + 6*(-1 + 6*lmMst1)*lmMst2 + 3*(-7 + 12*lmMst2)*lmMt +
        12*lmMsq*(4 - 3*(lmMst1 + lmMt)) + 36*pow2(lmMsq) - 36*pow2(lmMst2))*
        pow2(Mst1) - (1267 + 1632*lmMst2 - 6*lmMst1*(11 + 60*lmMst2) + 78*lmMt
        - 936*lmMst2*lmMt + 12*lmMsq*(-137 + 30*lmMst1 + 48*lmMst2 + 78*lmMt) -
        936*pow2(lmMsq) + 360*pow2(lmMst2))*pow2(Mst2)) + Mgl*(514500*(5 - 78*
        lmMst2 + 3*lmMst1*(5 + 12*lmMst2) + 3*(5 + 12*lmMst2)*lmMt + 12*lmMsq*(
        4 - 3*(lmMst1 + lmMt)) + 36*pow2(lmMsq) - 36*pow2(lmMst2))*pow2(Mst1) -
        (142320737 + 123053280*lmMst2 - 257250*lmMst1*(41 + 60*lmMst2) -
        257250*(77 + 228*lmMst2)*lmMt + 420*lmMsq*(-220709 + 36750*lmMst1 +
        82740*lmMst2 + 139650*lmMt) - 54419400*pow2(lmMsq) + 19668600*pow2(
        lmMst2))*pow2(Mst2))))/(2.7783e6*pow4(Msq)) + ((-5*((64*Dmglst2*(14 -
        15*lmMsq + 15*lmMst2) + 8*(43 - 30*lmMsq + 30*lmMst2)*Mgl)*pow2(Msq) +
        (18*Dmglst2*(13 - 28*lmMsq + 28*lmMst2) + (67 - 84*lmMsq + 84*lmMst2)*
        Mgl)*pow2(Mst2)))/(54.*pow2(Mst1)*pow4(Msq)) - (64*(1 + lmMst2)*(4*
        Dmglst2*lmMst2 + Mgl + lmMst2*Mgl))/(9.*pow4(Mst1)))*pow4(Mst2) + (
        3000*Dmglst2*((-87 + lmMst1*(40 + 6*lmMst2 - 6*lmMt) + 6*lmMst2*(-10 +
        lmMt) + 20*lmMt - 6*pow2(lmMst2))*pow4(Mst1) + (-5 + lmMst1 - 12*lmMst2
        + 6*lmMst1*lmMst2 - lmMt + 12*lmMst2*lmMt - 6*lmMsq*(-2 + lmMst1 +
        lmMst2 + 2*lmMt) + 12*pow2(lmMsq) - 6*pow2(lmMst2))*pow4(Mst2)) - 4*
        Mgl*(375*(-39 + lmMst1*(34 + 6*lmMst2 - 6*lmMt) + 6*lmMst2*(-8 + lmMt)
        + 14*lmMt - 6*pow2(lmMst2))*pow4(Mst1) + (7839 + 8660*lmMst2 - 375*
        lmMst1*(7 + 6*lmMst2) - 4875*lmMt - 6000*lmMst2*lmMt + 10*lmMsq*(-116 +
        225*lmMst1 + 285*lmMst2 + 600*lmMt) - 5550*pow2(lmMsq) + 2700*pow2(
        lmMst2))*pow4(Mst2)))/(675.*pow2(Msq)*pow2(Mst2)) + (-945*pow2(MuSUSY)*
        (3*Mgl*(8*pow2(Mst2)*(436 - 6*B4 + 6*D3 - 3*DN - 48*lmMst1 + (408 - 96*
        lmMst1)*lmMst2 + 24*lmMt - 48*(-4 + lmMst1)*pow2(lmMst2) + 48*pow3(
        lmMst2)) + pow2(Mst1)*(10667 - 96*B4 + 96*D3 - 48*DN - 3072*lmMst1 -
        384*(-25 + 9*lmMst1)*lmMst2 - 384*(-1 + lmMst1 - lmMst2)*lmMt - 384*(-
        13 + 4*lmMst1)*pow2(lmMst2) + 1536*pow3(lmMst2))) + 2*Dmglst2*(72*pow2(
        Mst2)*(180 - 2*B4 + 2*D3 - DN + 16*lmMst1 + 144*lmMst2 - 16*(-2 +
        lmMst1)*pow2(lmMst2) + 16*pow3(lmMst2)) + pow2(Mst1)*(28405 - 288*B4 +
        288*D3 - 144*DN + 5472*lmMst1 - 288*(-89 + 2*lmMst1)*lmMst2 + 576*(-2 +
        lmMst1 - lmMst2)*lmMt - 576*(-9 + 8*lmMst1)*pow2(lmMst2) + 4608*pow3(
        lmMst2)))) + 560*OepS2*(3*Mgl*(pow2(Mst1)*(627*pow2(Mst2) + 378*pow2(
        MuSUSY)) + 1168*pow4(Mst1) + 189*pow4(Mst2)) + 4*Dmglst2*(63*pow2(Mst1)
        *(23*pow2(Mst2) - 3*pow2(MuSUSY)) + 5582*pow4(Mst1) + 189*pow4(Mst2)))
        - 27*S2*(42*Mgl*(9*pow2(Mst1)*((4163 + 2090*lmMst1 - 2090*lmMst2)*pow2(
        Mst2) + 90*(-43 + 14*lmMst1 - 14*lmMst2)*pow2(MuSUSY)) + 8*(8653 +
        4380*lmMst1 - 4380*lmMst2)*pow4(Mst1) + 81*(-240*pow2(Mst2)*pow2(
        MuSUSY) + (109 + 70*lmMst1 - 70*lmMst2)*pow4(Mst2))) + 4*Dmglst2*(45*
        pow2(Mst1)*(23*(115 + 588*lmMst1 - 588*lmMst2)*pow2(Mst2) - 126*(65 +
        14*lmMst1 - 14*lmMst2)*pow2(MuSUSY)) + (2001242 + 2344440*lmMst1 -
        2344440*lmMst2)*pow4(Mst1) - 81*(3360*pow2(Mst2)*pow2(MuSUSY) + (1593 -
        980*lmMst1 + 980*lmMst2)*pow4(Mst2)))))/(76545.*pow4(Mst2)))/Mgl)*pow4(
        Mt) + (-(pow2(Mst1)*pow2(Mst2)*(95.16896090534979 - B4/9. + (2*D3)/9. -
        DN/6. - (445*lmMsq)/36. + lmMst1*(21.80061728395062 - (5*pow2(lmMsq))/
        2.) + (35*pow2(lmMsq))/12. + (lmMst2*(18499 + 20550*lmMst1 + 675*lmMsq*
        (-7 + 2*lmMst1) + 2025*pow2(lmMsq) - 6525*pow2(lmMst1)))/810. + ((-1649
        + 180*lmMsq)*pow2(lmMst1))/108. + ((289 - 360*lmMsq + 1398*lmMst1)*
        pow2(lmMst2))/108. - (79*pow3(lmMst1))/54. + (((50*Dmglst2*(25 + 12*
        lmMsq*(-4 + 3*lmMst1 - 3*lmMst2) + (48 - 36*lmMst1)*lmMst2 + 36*pow2(
        lmMst2)) + Mgl*(-2953 + 60*lmMsq*(17 + 6*lmMst1 - 6*lmMst2) - 636*
        lmMst2 - 24*lmMst1*(16 + 75*lmMst2) + 720*pow2(lmMst1) + 1080*pow2(
        lmMst2)))*pow2(Mst1))/(1080.*pow2(Msq)) - Dmglst2*(14.217769547325103 +
        (2*B4)/3. - (8*D3)/9. + (5*DN)/9. - (15*lmMsq)/2. - (5*pow2(lmMsq))/2.
         + lmMst1*(-
        24.21432098765432 + 5*pow2(lmMsq)) - lmMst2*(5*lmMsq*(-1 + 2*lmMst1) +
        5*pow2(lmMsq) + (-128893 - 86970*lmMst1 - 48150*pow2(lmMst1))/4050.) -
        (137*pow2(lmMst1))/135. + (10*lmMsq - (83*(34 + 15*lmMst1))/135.)*pow2(
        lmMst2) + (5*pow3(lmMst1))/27. - (77*pow3(lmMst2))/27.))/Mgl - (185*
        pow3(lmMst2))/54.)) + (35.682629270197204 + B4 + D3/9. - DN/9. - (655*
        lmMsq)/72. + lmMst1*(65.86121882086168 - (70*lmMsq)/9. + (5*pow2(lmMsq)
        )/12.) + (25*pow2(lmMsq))/12. - lmMst2*(33.61121882086168 + (99767*
        lmMst1)/3780. - (5*lmMsq*(13 + 9*lmMst1))/18. + (5*pow2(lmMsq))/12. - (
        305*pow2(lmMst1))/36.) + (27.599470899470898 - (5*lmMsq)/3.)*pow2(
        lmMst1) + (7.516137566137566 - (5*lmMsq)/6. - (473*lmMst1)/36.)*pow2(
        lmMst2) - (10*pow3(lmMst1))/27. + (136*pow3(lmMst2))/27. + (Dmglst2*(
        84.40344866031853 - (5*lmMsq*(2 + 3*lmMst1))/3. + (1011403*lmMst2)/
        1.1907e6 + (5*pow2(lmMsq))/2. + (4.326984126984127 + (11*lmMst2)/3.)*
        pow2(lmMst1) + lmMst1*(10.650581170739901 + (1474*lmMst2)/315. - (11*
        pow2(lmMst2))/3.) - (113*pow2(lmMst2))/35. - (11*pow3(lmMst1))/9. + (
        11*pow3(lmMst2))/9.))/Mgl)*pow4(Mst1) + (((Mgl*(-82320*(1318 - 15*
        lmMsq*(23 + 41*lmMst1 - 41*lmMst2) - 346*lmMst2 + lmMst1*(691 + 525*
        lmMst2) + 45*pow2(lmMst1) - 570*pow2(lmMst2))*pow2(Msq)*pow2(Mst1) +
        34300*(9403 - 720*B4 + 72*D3 - 36*DN - 630*lmMsq - 270*pow2(lmMsq) - 3*
        lmMst1*(5663 - 1950*lmMsq + 630*pow2(lmMsq)) - 954*pow2(lmMst1) + 3*
        lmMst2*(6377 - 6942*lmMst1 + 30*lmMsq*(-59 + 42*lmMst1) + 630*pow2(
        lmMsq) + 234*pow2(lmMst1)) - 18*(-1172 + 210*lmMsq + 399*lmMst1)*pow2(
        lmMst2) - 288*pow3(lmMst1) + 6768*pow3(lmMst2))*pow4(Msq) + 3*(6074157
        - 420*lmMsq*(9224 + 5705*lmMst1 - 6545*lmMst2) + 1208165*lmMst2 + 35*
        lmMst1*(76169 + 108780*lmMst2) - 176400*pow2(lmMsq) - 705600*pow2(
        lmMst1) - 3278100*pow2(lmMst2))*pow4(Mst1)) + 17150*Dmglst2*(240*(-77 +
        lmMsq*(6 + 60*lmMst1 - 60*lmMst2) + 50*lmMst2 - 4*lmMst1*(14 + 15*
        lmMst2) + 60*pow2(lmMst2))*pow2(Msq)*pow2(Mst1) + 2*(7337 - 432*B4 +
        576*D3 - 360*DN + 8100*lmMsq - 1620*pow2(lmMsq) - 216*lmMst1*(-7 + 15*
        pow2(lmMsq)) - 1404*pow2(lmMst1) + 36*lmMst2*(-689 - 528*lmMst1 + 90*
        lmMsq*(1 + 2*lmMst1) + 90*pow2(lmMsq) + 32*pow2(lmMst1)) - 36*(-239 +
        180*lmMsq + 336*lmMst1)*pow2(lmMst2) + 10944*pow3(lmMst2))*pow4(Msq) +
        15*(113 - 34*lmMst2 + 12*lmMsq*(1 - 10*lmMst1 + 10*lmMst2) + 2*lmMst1*(
        11 + 60*lmMst2) - 120*pow2(lmMst2))*pow4(Mst1)) - (pow2(Mst2)*(Mgl*(-
        41160*(2986 - 30*lmMsq*(43 + 40*lmMst1 - 40*lmMst2) - 735*lmMst2 + 75*
        lmMst1*(27 + 16*lmMst2) - 1200*pow2(lmMst2))*pow2(Msq)*pow2(Mst1) -
        308700*(405 + 770*lmMst2 + 64*lmMst1*(1 + lmMst2) - 30*lmMsq*(5 + 6*
        lmMst2) + 90*pow2(lmMsq) + 336*pow2(lmMst2))*pow4(Msq) + (45671209 -
        420*lmMsq*(51313 + 42630*lmMst1 - 44310*lmMst2) + 9501870*lmMst2 +
        10290*lmMst1*(1171 + 1740*lmMst2) - 352800*pow2(lmMsq) - 18257400*pow2(
        lmMst2))*pow4(Mst1)) + 51450*Dmglst2*(40*(-230 + 6*lmMsq*(16 + 17*
        lmMst1 - 17*lmMst2) + 13*lmMst2 - lmMst1*(109 + 102*lmMst2) + 102*pow2(
        lmMst2))*pow2(Msq)*pow2(Mst1) - 12*(-115 + (494 + 64*lmMst1)*lmMst2 -
        30*lmMsq*(-1 + 6*lmMst2) + 90*pow2(lmMsq) + 464*pow2(lmMst2))*pow4(Msq)
        - 5*(-793 + 12*lmMsq*(19 + 42*lmMst1 - 42*lmMst2) + 6*lmMst2 - 18*
        lmMst1*(13 + 28*lmMst2) + 504*pow2(lmMst2))*pow4(Mst1))))/pow2(Mst1))*
        pow4(Mst2))/(2.22264e7*pow4(Msq)) - (OepS2*(4*Dmglst2*(1944*pow2(Mst1)*
        pow2(Mst2) + 4199*pow4(Mst1) + 189*pow4(Mst2)) - 3*Mgl*(2058*pow2(Mst1)
        *pow2(Mst2) + 1945*pow4(Mst1) + 432*pow4(Mst2))))/2187. + (S2*(2*
        Dmglst2*(36*(-275 + 324*lmMst1 - 324*lmMst2)*pow2(Mst1)*pow2(Mst2) + (
        51635 + 25194*lmMst1 - 25194*lmMst2)*pow4(Mst1) + 81*(-1677 + 14*lmMst1
        - 14*lmMst2)*pow4(Mst2)) - 3*Mgl*(9*(3137 + 686*lmMst1 - 686*lmMst2)*
        pow2(Mst1)*pow2(Mst2) + 5*(3799 + 1167*lmMst1 - 1167*lmMst2)*pow4(Mst1)
        + 1296*(22 + lmMst1 - lmMst2)*pow4(Mst2))))/324.)/Mgl)*pow4(s2t) - (Mt*
        pow3(s2t)*(-87480*pow2(Mst1)*pow2(Mst2)*(19.17798353909465 + (68*B4)/9.
         + (2*DN)/9. + 40*lmMsq - (40*pow2(lmMsq))/3. - (4*lmMst1*(223 - 15*
        lmMsq + 15*pow2(lmMsq)))/9. - (8*pow2(lmMst1))/3. + (4*lmMst2*(-163 -
        229*lmMst1 + 15*lmMsq*(3 + 2*lmMst1) + 15*pow2(lmMsq) + 24*pow2(lmMst1)
        ))/9. - (4*(-107 + 30*lmMsq + 104*lmMst1)*pow2(lmMst2))/9. - ((20*(
        Dmglst2 - Mgl)*pow2(Mst1))/(9.*pow2(Msq)) + Dmglst2*(47.73436213991769
         - (68*B4)/
        9. - (2*DN)/9. + (80*lmMsq)/3. + (40*pow2(lmMsq))/3. + (4*lmMst1*(77 +
        135*lmMsq + 45*pow2(lmMsq)))/27. + (4*lmMst2*(409 + 183*lmMst1 - 45*
        lmMsq*(7 + 2*lmMst1) - 45*pow2(lmMsq) - 72*pow2(lmMst1)))/27. - (8*
        pow2(lmMst1))/3. + (4*(105 + 30*lmMsq + 104*lmMst1)*pow2(lmMst2))/9. -
        (320*pow3(lmMst2))/9.))/Mgl + (320*pow3(lmMst2))/9.) + ((Dmglst2*(
        35542243 - 8755920*lmMst2 + 291600*lmMsq*(7 - 6*lmMst1 + 6*lmMst2) +
        38880*(37 - 40*lmMst2)*pow2(lmMst1) + 1438560*pow2(lmMst2) + 720*
        lmMst1*(8705 - 3996*lmMst2 + 5616*pow2(lmMst2)) - 311040*pow3(lmMst1) -
        2177280*pow3(lmMst2)) - 15*Mgl*(542267 - 29376*lmMst2 + 19440*lmMsq*(1
        - 2*lmMst1 + 2*lmMst2) + 864*(205 + 216*lmMst2)*pow2(lmMst1) + 177120*
        pow2(lmMst2) - 1728*lmMst1*(31 + 205*lmMst2 + 204*pow2(lmMst2)) - 6912*
        pow3(lmMst1) + 172800*pow3(lmMst2)))*pow4(Mst1))/Mgl + ((-27*(15*Mgl*(
        320*(25 - 6*lmMsq*(2 + lmMst1 - lmMst2) + 4*lmMst2 + lmMst1*(8 + 6*
        lmMst2) - 6*pow2(lmMst2))*pow2(Msq)*pow2(Mst1) + 3*(17177 + 1760*B4 -
        16*DN - 5760*lmMsq + 960*pow2(lmMsq) - 16*lmMst1*(199 - 30*lmMsq + 30*
        pow2(lmMsq)) + 16*lmMst2*(1163 - 322*lmMst1 + 30*lmMsq*(-5 + 2*lmMst1)
        + 30*pow2(lmMsq) - 16*pow2(lmMst1)) - 256*pow2(lmMst1) - 32*(-281 + 30*
        lmMsq + 72*lmMst1)*pow2(lmMst2) + 2560*pow3(lmMst2))*pow4(Msq) + 10*(-1
        + 12*lmMsq - 12*lmMst2)*pow4(Mst1)) + Dmglst2*(28800*(9 - 3*lmMsq*(2 +
        lmMst1 - lmMst2) + 4*lmMst2 + lmMst1*(2 + 3*lmMst2) - 3*pow2(lmMst2))*
        pow2(Msq)*pow2(Mst1) + (825997 + 188640*B4 - 3600*DN - 648000*lmMsq +
        43200*pow2(lmMsq) - 720*lmMst1*(-173 + 90*lmMsq + 30*pow2(lmMsq)) +
        720*lmMst2*(1495 - 130*lmMst1 + 30*lmMsq*(-1 + 2*lmMst1) + 30*pow2(
        lmMsq) - 16*pow2(lmMst1)) + 11520*pow2(lmMst1) - 1440*(-169 + 30*lmMsq
        + 136*lmMst1)*pow2(lmMst2) + 207360*pow3(lmMst2))*pow4(Msq) + 150*(23 +
        12*lmMsq - 12*lmMst2)*pow4(Mst1)))*pow4(Mst2))/pow4(Msq) + 160*OepS2*(
        3*Mgl*(438*pow2(Mst1)*pow2(Mst2) + 439*pow4(Mst1) + 189*pow4(Mst2)) -
        Dmglst2*(7218*pow2(Mst1)*pow2(Mst2) + 15571*pow4(Mst1) + 945*pow4(Mst2)
        )) + 108*S2*(Dmglst2*(18*(14033 + 12030*lmMst1 - 12030*lmMst2)*pow2(
        Mst1)*pow2(Mst2) + (836021 + 467130*lmMst1 - 467130*lmMst2)*pow4(Mst1)
        + 81*(-453 + 350*lmMst1 - 350*lmMst2)*pow4(Mst2)) - 15*Mgl*(18*(347 +
        146*lmMst1 - 146*lmMst2)*pow2(Mst1)*pow2(Mst2) + (6121 + 2634*lmMst1 -
        2634*lmMst2)*pow4(Mst1) + 81*(-1 + 14*lmMst1 - 14*lmMst2)*pow4(Mst2)))
        + (270*(48*(15*Dmglst2*(95 + 12*lmMsq*(-2 + lmMst1 - lmMst2) + 32*
        lmMst2 - 4*lmMst1*(2 + 3*lmMst2) + 12*pow2(lmMst2)) + 5*Mgl*(71 + 12*
        lmMsq*(-2 + lmMst1 - lmMst2) + 40*lmMst2 - 4*lmMst1*(4 + 3*lmMst2) +
        12*pow2(lmMst2)))*pow2(Msq) - 15*Mgl*(223 - 12*lmMsq*(13 + 6*lmMst1 -
        6*lmMst2) + 102*lmMst2 + 18*lmMst1*(3 + 4*lmMst2) - 72*pow2(lmMst2))*
        pow2(Mst1) + 5*Dmglst2*(-3*(823 - 12*lmMsq*(61 + 30*lmMst1 - 30*lmMst2)
        + 606*lmMst2 + 18*lmMst1*(7 + 20*lmMst2) - 360*pow2(lmMst2))*pow2(Mst1)
        + (5591 + 12*lmMsq*(-181 + 90*lmMst1 - 90*lmMst2) + 2550*lmMst2 - 54*
        lmMst1*(7 + 20*lmMst2) + 1080*pow2(lmMst2))*pow2(Mst2)) - (2304*(1 +
        lmMst2)*(-Dmglst2 + 3*Dmglst2*lmMst2 + Mgl + lmMst2*Mgl)*pow4(Msq))/
        pow2(Mst1))*pow6(Mst2))/pow4(Msq))/Mgl))/(87480.*Mst2) + pow2(Mt)*pow2(
        s2t)*((2*pow2(Mst1)*(5*Mgl*(517702 - 2700*D3 + 2700*DN + 10125*lmMsq -
        64800*(lmMst1 - lmMst2)*lmMt - 30375*pow2(lmMsq) - 15*lmMst1*(43538 -
        5850*lmMsq + 675*pow2(lmMsq)) + 15*lmMst2*(47993 - 10020*lmMst1 - 450*
        lmMsq*(4 + 3*lmMst1) + 675*pow2(lmMsq) - 1035*pow2(lmMst1)) + 225*(-887
        + 90*lmMsq)*pow2(lmMst1) - 225*(-1195 + 339*lmMst1)*pow2(lmMst2) -
        29250*pow3(lmMst1) + 121050*pow3(lmMst2)) + Dmglst2*(2416741 + 1789710*
        lmMst2 + 67500*lmMsq*(8 - 3*lmMst1 + 12*lmMst2) - 324000*(2 + lmMst2)*
        lmMt - 303750*pow2(lmMsq) + 450*(-341 + 30*lmMst2)*pow2(lmMst1) + 60*
        lmMst1*(34484 - 12885*lmMst2 + 5400*lmMt - 21825*pow2(lmMst2)) +
        224550*pow2(lmMst2) - 4500*pow3(lmMst1) + 1300500*pow3(lmMst2))))/(
        30375.*Mgl) + ((9000*Dmglst2*(4 - lmMst1 + lmMst2) + 2*Mgl*(5108 + 15*
        lmMst1*(718 - 75*lmMst2) - 15*lmMsq*(218 + 195*lmMst1 - 75*lmMst2) -
        7500*lmMst2 + 900*pow2(lmMsq) + 2025*pow2(lmMst1)))/(2025.*Mgl*pow2(
        Msq)) + (196.55862427463637 + (-293.24826908541195 + (40*lmMsq)/3.)*
        lmMst1 + lmMst2*(305.24826908541195 - (40*lmMsq)/3. - (32894*lmMst1)/
        945. - (286*pow2(lmMst1))/9.) - (50753*pow2(lmMst1))/945. + (32*lmMt*(-
        2*lmMst1*(1 + lmMst2) + lmMst2*(2 + lmMst2) + pow2(lmMst1)))/3. + (
        88.51534391534392 - (2*lmMst1)/9.)*pow2(lmMst2) + (190*pow3(lmMst1))/
        27. - (Dmglst2*(190.3527287429963 - (160*lmMsq)/3. + (-
        535.9278609221467 + (200*lmMsq)/3.)*lmMst1 + lmMst2*(315.4834164777022
         - (200*lmMsq)/
        3. + (84592*lmMst1)/315. - (208*pow2(lmMst1))/9.) - (35576*pow2(lmMst1)
        )/315. + (8*(-6127 + 5110*lmMst1)*pow2(lmMst2))/315. + (64*lmMt*(2 + 5*
        lmMst2 - lmMst1*(5 + 2*lmMst2) + pow2(lmMst1) + pow2(lmMst2)))/3. + (
        16*pow3(lmMst1))/27. - (2896*pow3(lmMst2))/27.))/Mgl + (674*pow3(
        lmMst2))/27.)/pow2(Mst2))*pow4(Mst1) + pow2(Mst2)*(252.5348765432099 +
        (8*D3)/9. - (8*DN)/9. - (100*lmMsq)/3. + (32*lmMt)/3. + 20*pow2(lmMsq)
        + (2*lmMst1*(763 - 300*lmMsq + 45*pow2(lmMsq)))/27. + (124*pow2(lmMst1)
        )/9. - (2*lmMst2*(30*lmMsq*(8 + 3*lmMst1) + 45*pow2(lmMsq) + 13*(-209 -
        30*lmMst1 + 9*pow2(lmMst1))))/27. + (2*(256 + 30*lmMsq + 31*lmMst1)*
        pow2(lmMst2))/9. + (32*pow3(lmMst1))/9. - (16*pow3(lmMst2))/9. + (2*(((
        11250*Dmglst2*(3 - 4*lmMsq + 4*lmMst2) - Mgl*(3091 - 10230*lmMst2 + 60*
        lmMst1*(92 + 75*lmMst2) + 30*lmMsq*(157 - 60*(lmMst1 + lmMst2)) + 1800*
        pow2(lmMsq) - 1350*pow2(lmMst1) - 1350*pow2(lmMst2)))*pow2(Mst1))/pow2(
        Msq) + 15*Dmglst2*(28373 - 240*B4 + 240*D3 - 120*DN + 750*lmMsq + 1380*
        lmMst1 + 2700*pow2(lmMsq) + 690*pow2(lmMst1) - 20*lmMst2*(-1676 + 270*
        lmMsq - 27*lmMst1 + 48*pow2(lmMst1)) - 30*(-391 + 32*lmMst1)*pow2(
        lmMst2) + 1920*pow3(lmMst2))))/(2025.*Mgl) - ((2.3647986178598424 - (
        991*lmMsq)/441. + ((-194 + 385*lmMsq)*lmMst1)/735. + (2.511111111111111
         - lmMsq + lmMst1)*lmMst2 - (5*Dmglst2*(1 + 4*lmMsq - 4*lmMst2))/(6.*
        Mgl) + (5*pow2(lmMsq))/21. - (16*pow2(lmMst1))/21.)*pow4(Mst1))/pow4(
        Msq)) + ((15*(8575*Dmglst2*((65 - 324*lmMsq + 324*lmMst2)*pow2(Mst1) +
        (-203 + 540*lmMsq - 540*lmMst2)*pow2(Mst2)) - 2*Mgl*((9261*lmMst1*(7 +
        5*lmMst2) + 105*lmMsq*(2929 - 441*lmMst1 + 231*lmMst2) + 11025*pow2(
        lmMsq) - 2*(65299 + 186186*lmMst2 + 17640*pow2(lmMst2)))*pow2(Mst1) + (
        242073 - 105*lmMsq*(3909 + 49*lmMst1 - 259*lmMst2) + 399812*lmMst2 +
        343*lmMst1*(31 + 15*lmMst2) - 11025*pow2(lmMsq) - 16170*pow2(lmMst2))*
        pow2(Mst2))) - (686*pow2(Msq)*(Mgl*(750*(-43 + 30*lmMsq - 30*lmMst2)*
        pow2(Mst1)*pow2(Mst2) + 225*pow2(Msq)*((309 + 578*lmMst2 + 64*lmMst1*(1
        + lmMst2) - 30*lmMsq*(5 + 6*lmMst2) + 90*pow2(lmMsq) + 240*pow2(lmMst2)
        )*pow2(Mst1) - 32*pow2(1 + lmMst2)*pow2(Mst2)) + 2*(22642 - 15*lmMsq*(
        1282 + 75*lmMst1 - 195*lmMst2) + 16230*lmMst2 + 375*lmMst1*(8 + 3*
        lmMst2) - 900*pow2(lmMsq) - 2025*pow2(lmMst2))*pow4(Mst1)) + 150*
        Dmglst2*(3*pow2(Msq)*((-179 + 30*lmMsq + 238*lmMst2 - 180*lmMsq*lmMst2
        + 64*lmMst1*lmMst2 + 90*pow2(lmMsq) + 272*pow2(lmMst2))*pow2(Mst1) -
        64*lmMst2*(1 + lmMst2)*pow2(Mst2)) - 10*(4*(14 - 15*lmMsq + 15*lmMst2)*
        pow2(Mst1)*pow2(Mst2) + (-107 + 120*lmMsq - 120*lmMst2)*pow4(Mst1)))))/
        pow4(Mst1))*pow4(Mst2))/(1.38915e6*Mgl*pow4(Msq)) - (pow2(MuSUSY)*(
        15829 - 720*B4 + 72*D3 - 36*DN - 3330*lmMsq + 1350*pow2(lmMsq) - 3*
        lmMst1*(5279 - 1950*lmMsq + 630*pow2(lmMsq)) - 954*pow2(lmMst1) + 3*
        lmMst2*(10421 - 6558*lmMst1 + 30*lmMsq*(-95 + 42*lmMst1) + 630*pow2(
        lmMsq) + 234*pow2(lmMst1)) - 18*(-1460 + 210*lmMsq + 399*lmMst1)*pow2(
        lmMst2) - 288*pow3(lmMst1) + 6768*pow3(lmMst2) + (162*(((Mgl*(111 -
        3028*lmMst2 + 90*lmMsq*(4 - 13*lmMst1 + 13*lmMst2) + 2*lmMst1*(1334 +
        675*lmMst2) - 90*pow2(lmMst1) - 1260*pow2(lmMst2)) + 100*Dmglst2*(69 +
        (-53 + 42*lmMsq)*lmMst2 + lmMst1*(53 - 42*lmMsq + 42*lmMst2) - 42*pow2(
        lmMst2)))*pow2(Mst1))/(135.*pow2(Msq)) + Dmglst2*(19.734567901234566 -
        (8*B4)/3. + (32*D3)/9. - (20*DN)/9. + (170*lmMsq)/3. + lmMst1*(
        9.333333333333334 - 20*pow2(lmMsq)) + 10*pow2(lmMsq) - (26*pow2(lmMst1)
        )/3. + (2*lmMst2*(-291 - 464*lmMst1 + 90*lmMsq*(-1 + 2*lmMst1) + 90*
        pow2(lmMsq) + 32*pow2(lmMst1)))/9. - (2*(-607 + 180*lmMsq + 336*lmMst1)
        *pow2(lmMst2))/9. + (608*pow3(lmMst2))/9.)))/Mgl - (162*pow2(Mst1)*(
        204.20053497942388 + (76*B4)/9. - (2*DN)/9. - (50*lmMsq)/3. + (2*
        lmMst1*(57971 - 14625*lmMsq + 2700*pow2(lmMsq)))/405. + ((-1331 + 180*
        lmMsq)*pow2(lmMst1))/27. - (2*lmMst2*(52436 - 70455*lmMst1 + 225*lmMsq*
        (-65 + 36*lmMst1) + 2700*pow2(lmMsq) + 8280*pow2(lmMst1)))/405. + ((-
        8063 + 900*lmMsq + 3792*lmMst1)*pow2(lmMst2))/27. - (62*pow3(lmMst1))/
        27. - (2626*pow3(lmMst2))/27. - (Dmglst2*(109.11799176954733 - (8*B4)/
        3. + (32*D3)/9. - (20*DN)/9. + 80*lmMsq - lmMst1*(78.19061728395062 +
        20*pow2(lmMsq)) - (2888*pow2(lmMst1))/135. + lmMst2*(40*lmMsq*lmMst1 +
        20*pow2(lmMsq) + (4*(-21616 - 64515*lmMst1 + 31275*pow2(lmMst1)))/2025.
        ) - (4*(-5023 + 1350*lmMsq + 6285*lmMst1)*pow2(lmMst2))/135. + (20*
        pow3(lmMst1))/27. + (3340*pow3(lmMst2))/27.))/Mgl))/pow2(Mst2) + 162*
        pow4(Mst1)*((1.0702990137854083 + (9571*lmMsq)/26460. + (
        4.249508692365835 - (169*lmMsq)/63.)*lmMst1 + (-4.6112244897959185 + (
        19*lmMsq)/7. + (31*lmMst1)/9.)*lmMst2 - pow2(lmMsq)/63. - (8*pow2(
        lmMst1))/21. + (5*Dmglst2*(216 + (-95 + 132*lmMsq)*lmMst2 + lmMst1*(95
        - 132*lmMsq + 132*lmMst2) - 132*pow2(lmMst2)))/(54.*Mgl) - (194*pow2(
        lmMst2))/63.)/pow4(Msq) - (363.3804294212688 + (76*B4)/9. - (2*DN)/9. -
        (35*lmMsq)/2. + lmMst1*(211.3489518770471 - (695*lmMsq)/9. + (40*pow2(
        lmMsq))/3.) - (214.87936507936507 - 20*lmMsq)*pow2(lmMst1) - lmMst2*(
        190.46006298815823 - (71398*lmMst1)/105. + (5*lmMsq*(-139 + 120*lmMst1)
        )/9. + (40*pow2(lmMsq))/3. + (334*pow2(lmMst1))/3.) + ((-146507 +
        14700*lmMsq + 91070*lmMst1)*pow2(lmMst2))/315. - (44*pow3(lmMst1))/9. -
        (1556*pow3(lmMst2))/9. - (Dmglst2*(536.1152102791342 - (8*B4)/3. + (32*
        D3)/9. - (20*DN)/9. + 90*lmMsq - (123.11224321827497 + 20*lmMsq*(1 +
        lmMsq))*lmMst1 - lmMst2*(17.33220122616948 - 20*lmMsq*(1 + lmMsq) + (
        133.04550264550264 - 40*lmMsq)*lmMst1 - (1180*pow2(lmMst1))/9.) - (
        15886*pow2(lmMst1))/945. + (149.85608465608465 - 40*lmMsq - (2812*
        lmMst1)/9.)*pow2(lmMst2) - (92*pow3(lmMst1))/27. + (4988*pow3(lmMst2))/
        27.))/Mgl)/pow4(Mst2)) + (162*(((Mgl*(341 + 642*lmMst2 + 64*lmMst1*(1 +
        lmMst2) - 30*lmMsq*(5 + 6*lmMst2) + 90*pow2(lmMsq) + 272*pow2(lmMst2))
        + 2*Dmglst2*(-115 + (366 + 64*lmMst1)*lmMst2 - 30*lmMsq*(-1 + 6*lmMst2)
        + 90*pow2(lmMsq) + 336*pow2(lmMst2)))*pow2(Mst2))/(18.*pow2(Mst1)) + (
        pow2(Mst2)*(Mgl*(2*(4167613 - 19932360*lmMst2 + 20580*lmMst1*(701 +
        540*lmMst2) + 420*lmMsq*(13109 - 26460*lmMst1 + 25620*lmMst2) + 176400*
        pow2(lmMsq) - 10936800*pow2(lmMst2))*pow2(Mst1) + (41220947 - 420*
        lmMsq*(12479 + 69090*lmMst1 - 69930*lmMst2) - 21234990*lmMst2 + 10290*
        lmMst1*(2573 + 2820*lmMst2) - 176400*pow2(lmMsq) - 29194200*pow2(
        lmMst2))*pow2(Mst2)) + 514500*Dmglst2*(2*(219 + (-95 + 132*lmMsq)*
        lmMst2 + lmMst1*(95 - 132*lmMsq + 132*lmMst2) - 132*pow2(lmMst2))*pow2(
        Mst1) + (557 - 224*lmMst2 + 12*lmMsq*(1 - 32*lmMst1 + 32*lmMst2) + 4*
        lmMst1*(53 + 96*lmMst2) - 384*pow2(lmMst2))*pow2(Mst2))))/(1.11132e7*
        pow4(Msq)) + ((-5*((64*Dmglst2*(14 - 15*lmMsq + 15*lmMst2) + 8*(43 -
        30*lmMsq + 30*lmMst2)*Mgl)*pow2(Msq) + (18*Dmglst2*(13 - 28*lmMsq + 28*
        lmMst2) + (67 - 84*lmMsq + 84*lmMst2)*Mgl)*pow2(Mst2)))/(216.*pow2(
        Mst1)*pow4(Msq)) - (16*(1 + lmMst2)*(4*Dmglst2*lmMst2 + Mgl + lmMst2*
        Mgl))/(9.*pow4(Mst1)))*pow4(Mst2) + (Mgl*((1725 + (-7006 + 2640*lmMsq)*
        lmMst2 + lmMst1*(7006 - 2640*lmMsq + 4800*lmMst2) - 1080*pow2(lmMst1) -
        3720*pow2(lmMst2))*pow4(Mst1) + 2*(836 - 2235*lmMst2 + 75*lmMst1*(27 +
        16*lmMst2) + 30*lmMsq*(7 - 40*lmMst1 + 40*lmMst2) - 1200*pow2(lmMst2))*
        pow4(Mst2)) + 50*Dmglst2*((291 + 2*(-103 + 84*lmMsq)*lmMst2 + 2*lmMst1*
        (103 - 84*lmMsq + 84*lmMst2) - 168*pow2(lmMst2))*pow4(Mst1) + 2*(118 +
        109*lmMst1 + (-133 + 102*lmMst1)*lmMst2 + 6*lmMsq*(4 - 17*lmMst1 + 17*
        lmMst2) - 102*pow2(lmMst2))*pow4(Mst2)))/(270.*pow2(Msq)*pow2(Mst2))))/
        Mgl))/162. + (-3*Mgl*(6*pow2(Mst1)*pow2(Mst2)*((25760*OepS2 - 162*(4783
        + 3220*lmMst1 - 3220*lmMst2)*S2)*pow2(Mst2) + 5*(3896*OepS2 - 81*(9473
        + 974*lmMst1 - 974*lmMst2)*S2)*pow2(MuSUSY)) + ((276680*OepS2 - 27*(
        406627 + 207510*lmMst1 - 207510*lmMst2)*S2)*pow2(Mst2) + 10*(29428*
        OepS2 - 27*(160997 + 22071*lmMst1 - 22071*lmMst2)*S2)*pow2(MuSUSY))*
        pow4(Mst1) + 27*((920*OepS2 - 81*(-891 + 230*lmMst1 - 230*lmMst2)*S2)*
        pow2(Mst2) + 160*(4*OepS2 - 81*(22 + lmMst1 - lmMst2)*S2)*pow2(MuSUSY))
        *pow4(Mst2)) + 4*Dmglst2*(54*pow2(Mst1)*pow2(Mst2)*((2000*OepS2 - 162*(
        31 + 250*lmMst1 - 250*lmMst2)*S2)*pow2(Mst2) + 5*(344*OepS2 + 9*(15643
        - 774*lmMst1 + 774*lmMst2)*S2)*pow2(MuSUSY)) + 2*(9*(30760*OepS2 - 27*(
        28283 + 23070*lmMst1 - 23070*lmMst2)*S2)*pow2(Mst2) + 10*(17308*OepS2 +
        27*(93919 - 12981*lmMst1 + 12981*lmMst2)*S2)*pow2(MuSUSY))*pow4(Mst1) +
        135*(56*OepS2 - 81*(-1677 + 14*lmMst1 - 14*lmMst2)*S2)*pow2(MuSUSY)*
        pow4(Mst2) - 2768742*S2*pow6(Mst2)))/(21870.*Mgl*pow4(Mst2))) + (pow2(
        Mt)*pow2(MuSUSY)*((pow2(Mt)*(3*Mgl*(8*pow2(Mst2)*(436 - 6*B4 + 6*D3 -
        3*DN - 48*lmMst1 + (408 - 96*lmMst1)*lmMst2 + 24*lmMt - 972*S2 - 48*(-4
        + lmMst1)*pow2(lmMst2) + 48*pow3(lmMst2)) + pow2(Mst1)*(10667 - 96*B4 +
        96*D3 - 48*DN - 3072*lmMst1 - 384*(-25 + 9*lmMst1)*lmMst2 - 384*(-1 +
        lmMst1 - lmMst2)*lmMt - 224*OepS2 + 324*(-43 + 14*lmMst1 - 14*lmMst2)*
        S2 - 384*(-13 + 4*lmMst1)*pow2(lmMst2) + 1536*pow3(lmMst2))) + 2*
        Dmglst2*(72*pow2(Mst2)*(180 - 2*B4 + 2*D3 - DN + 16*lmMst1 + 144*lmMst2
        - 216*S2 - 16*(-2 + lmMst1)*pow2(lmMst2) + 16*pow3(lmMst2)) + pow2(
        Mst1)*(28405 - 288*B4 + 288*D3 - 144*DN + 5472*lmMst1 - 288*(-89 + 2*
        lmMst1)*lmMst2 + 576*(-2 + lmMst1 - lmMst2)*lmMt + 224*OepS2 - 324*(65
        + 14*lmMst1 - 14*lmMst2)*S2 - 576*(-9 + 8*lmMst1)*pow2(lmMst2) + 4608*
        pow3(lmMst2)))))/(81.*Mgl*pow4(Mst2)) + pow2(s2t)*(97.70987654320987 -
        (40*B4)/9. + (4*D3)/9. - (2*DN)/9. - (185*lmMsq)/9. + (25*pow2(lmMsq))/
        3. - (lmMst1*(5279 - 1950*lmMsq + 630*pow2(lmMsq)))/54. - (53*pow2(
        lmMst1))/9. + (lmMst2*(10421 - 6558*lmMst1 + 30*lmMsq*(-95 + 42*lmMst1)
        + 630*pow2(lmMsq) + 234*pow2(lmMst1)))/54. - ((-1460 + 210*lmMsq + 399*
        lmMst1)*pow2(lmMst2))/9. - (16*pow3(lmMst1))/9. + (376*pow3(lmMst2))/9.
         - (pow2(Mst1)*(204.20053497942388 + (76*B4)/9. - (2*DN)/9. - (50*
        lmMsq)/3. + (2*lmMst1*(57971 - 14625*lmMsq + 2700*pow2(lmMsq)))/405. +
        ((-1331 + 180*lmMsq)*pow2(lmMst1))/27. - (2*lmMst2*(52436 - 70455*
        lmMst1 + 225*lmMsq*(-65 + 36*lmMst1) + 2700*pow2(lmMsq) + 8280*pow2(
        lmMst1)))/405. + ((-8063 + 900*lmMsq + 3792*lmMst1)*pow2(lmMst2))/27. -
        (62*pow3(lmMst1))/27. - (2626*pow3(lmMst2))/27. - (Dmglst2*(
        109.11799176954733 - (8*B4)/3. + (32*D3)/9. - (20*DN)/9. + 80*lmMsq -
        lmMst1*(78.19061728395062 + 20*pow2(lmMsq)) - (2888*pow2(lmMst1))/135.
         + lmMst2*(40*lmMsq*lmMst1 + 20*pow2(lmMsq) + (4*(-21616 - 64515*lmMst1
        + 31275*pow2(lmMst1)))/2025.) - (4*(-5023 + 1350*lmMsq + 6285*lmMst1)*
        pow2(lmMst2))/135. + (20*pow3(lmMst1))/27. + (3340*pow3(lmMst2))/27.))/
        Mgl))/pow2(Mst2) + ((6*(Mgl*(111 - 3028*lmMst2 + 90*lmMsq*(4 - 13*
        lmMst1 + 13*lmMst2) + 2*lmMst1*(1334 + 675*lmMst2) - 90*pow2(lmMst1) -
        1260*pow2(lmMst2)) + 100*Dmglst2*(69 + (-53 + 42*lmMsq)*lmMst2 +
        lmMst1*(53 - 42*lmMsq + 42*lmMst2) - 42*pow2(lmMst2)))*pow2(Mst1))/
        pow2(Msq) + 5*Dmglst2*(3197 - 432*B4 + 576*D3 - 360*DN + 9180*lmMsq +
        216*lmMst1*(7 - 15*pow2(lmMsq)) + 1620*pow2(lmMsq) - 1404*pow2(lmMst1)
        + 36*lmMst2*(-291 - 464*lmMst1 + 90*lmMsq*(-1 + 2*lmMst1) + 90*pow2(
        lmMsq) + 32*pow2(lmMst1)) - 36*(-607 + 180*lmMsq + 336*lmMst1)*pow2(
        lmMst2) + 10944*pow3(lmMst2)))/(810.*Mgl) - S2*(1056 + 48*lmMst1 - 48*
        lmMst2 + ((9473 + 974*lmMst1 - 974*lmMst2 + (3*Dmglst2*(
        6952.444444444444 - 344*lmMst1 + 344*lmMst2))/Mgl)*pow2(Mst1))/(3.*
        pow2(Mst2)) + (Dmglst2*(3354 - 28*lmMst1 + 28*lmMst2) + ((8*Dmglst2*(
        93919 - 12981*lmMst1 + 12981*lmMst2) + 3*(160997 + 22071*lmMst1 -
        22071*lmMst2)*Mgl)*pow4(Mst1))/(81.*pow4(Mst2)))/Mgl) + pow4(Mst1)*((
        1.0702990137854083 + (9571*lmMsq)/26460. + ((56221 - 35490*lmMsq)*
        lmMst1)/13230. + (-4.6112244897959185 + (19*lmMsq)/7. + (31*lmMst1)/9.)
        *lmMst2 - pow2(lmMsq)/63. - (8*pow2(lmMst1))/21. + (5*Dmglst2*(216 + (-
        95 + 132*lmMsq)*lmMst2 + lmMst1*(95 - 132*lmMsq + 132*lmMst2) - 132*
        pow2(lmMst2)))/(54.*Mgl) - (194*pow2(lmMst2))/63.)/pow4(Msq) - (
        363.3804294212688 + (76*B4)/9. - (2*DN)/9. - (35*lmMsq)/2. + lmMst1*(
        211.3489518770471 - (695*lmMsq)/9. + (40*pow2(lmMsq))/3.) - (
        214.87936507936507 - 20*lmMsq)*pow2(lmMst1) - lmMst2*(
        190.46006298815823 - (71398*lmMst1)/105. + (5*lmMsq*(-139 + 120*lmMst1)
        )/9. + (40*pow2(lmMsq))/3. + (334*pow2(lmMst1))/3.) + ((-146507 +
        14700*lmMsq + 91070*lmMst1)*pow2(lmMst2))/315. - (44*pow3(lmMst1))/9. -
        (1556*pow3(lmMst2))/9. - (Dmglst2*(536.1152102791342 - (8*B4)/3. + (32*
        D3)/9. - (20*DN)/9. + 90*lmMsq - (123.11224321827497 + 20*lmMsq*(1 +
        lmMsq))*lmMst1 - lmMst2*(17.33220122616948 - 20*lmMsq*(1 + lmMsq) + (
        133.04550264550264 - 40*lmMsq)*lmMst1 - (1180*pow2(lmMst1))/9.) - (
        15886*pow2(lmMst1))/945. + (149.85608465608465 - 40*lmMsq - (2812*
        lmMst1)/9.)*pow2(lmMst2) - (92*pow3(lmMst1))/27. + (4988*pow3(lmMst2))/
        27.))/Mgl)/pow4(Mst2)) + (((Mgl*(341 + 642*lmMst2 + 64*lmMst1*(1 +
        lmMst2) - 30*lmMsq*(5 + 6*lmMst2) + 90*pow2(lmMsq) + 272*pow2(lmMst2))
        + 2*Dmglst2*(-115 + (366 + 64*lmMst1)*lmMst2 - 30*lmMsq*(-1 + 6*lmMst2)
        + 90*pow2(lmMsq) + 336*pow2(lmMst2)))*pow2(Mst2))/(18.*pow2(Mst1)) + (
        pow2(Mst2)*(Mgl*(2*(4167613 - 19932360*lmMst2 + 20580*lmMst1*(701 +
        540*lmMst2) + 420*lmMsq*(13109 - 26460*lmMst1 + 25620*lmMst2) + 176400*
        pow2(lmMsq) - 10936800*pow2(lmMst2))*pow2(Mst1) + (41220947 - 420*
        lmMsq*(12479 + 69090*lmMst1 - 69930*lmMst2) - 21234990*lmMst2 + 10290*
        lmMst1*(2573 + 2820*lmMst2) - 176400*pow2(lmMsq) - 29194200*pow2(
        lmMst2))*pow2(Mst2)) + 514500*Dmglst2*(2*(219 + (-95 + 132*lmMsq)*
        lmMst2 + lmMst1*(95 - 132*lmMsq + 132*lmMst2) - 132*pow2(lmMst2))*pow2(
        Mst1) + (557 - 224*lmMst2 + 12*lmMsq*(1 - 32*lmMst1 + 32*lmMst2) + 4*
        lmMst1*(53 + 96*lmMst2) - 384*pow2(lmMst2))*pow2(Mst2))))/(1.11132e7*
        pow4(Msq)) + ((-5*((64*Dmglst2*(14 - 15*lmMsq + 15*lmMst2) + 8*(43 -
        30*lmMsq + 30*lmMst2)*Mgl)*pow2(Msq) + (18*Dmglst2*(13 - 28*lmMsq + 28*
        lmMst2) + (67 - 84*lmMsq + 84*lmMst2)*Mgl)*pow2(Mst2)))/(216.*pow2(
        Mst1)*pow4(Msq)) - (16*(1 + lmMst2)*(4*Dmglst2*lmMst2 + Mgl + lmMst2*
        Mgl))/(9.*pow4(Mst1)))*pow4(Mst2) - (4*OepS2*(4*Dmglst2*(2322*pow2(
        Mst1)*pow2(Mst2) + 8654*pow4(Mst1) + 189*pow4(Mst2)) - 3*Mgl*(2922*
        pow2(Mst1)*pow2(Mst2) + 7357*pow4(Mst1) + 432*pow4(Mst2))))/(2187.*
        pow4(Mst2)) + (Mgl*((1725 + (-7006 + 2640*lmMsq)*lmMst2 + lmMst1*(7006
        - 2640*lmMsq + 4800*lmMst2) - 1080*pow2(lmMst1) - 3720*pow2(lmMst2))*
        pow4(Mst1) + 2*(836 - 2235*lmMst2 + 75*lmMst1*(27 + 16*lmMst2) + 30*
        lmMsq*(7 - 40*lmMst1 + 40*lmMst2) - 1200*pow2(lmMst2))*pow4(Mst2)) +
        50*Dmglst2*((291 + 2*(-103 + 84*lmMsq)*lmMst2 + 2*lmMst1*(103 - 84*
        lmMsq + 84*lmMst2) - 168*pow2(lmMst2))*pow4(Mst1) + 2*(118 + 109*lmMst1
        + (-133 + 102*lmMst1)*lmMst2 + 6*lmMsq*(4 - 17*lmMst1 + 17*lmMst2) -
        102*pow2(lmMst2))*pow4(Mst2)))/(270.*pow2(Msq)*pow2(Mst2)))/Mgl) + (Mt*
        s2t*(14580*pow2(Mst1)*pow2(Mst2)*(1035.3004115226338 + (1016*B4)/9. - (
        4*DN)/9. - 240*lmMsq + (80*pow2(lmMsq))/3. - (8*lmMst1*(422 - 45*lmMsq
        + 45*pow2(lmMsq)))/9. - (176*pow2(lmMst1))/9. + (8*lmMst2*(1096 + 15*
        lmMsq*(-7 + 6*lmMst1 - 6*lmMst2) + 717*lmMst2 - lmMst1*(551 + 248*
        lmMst2) + 45*pow2(lmMsq) + 8*pow2(lmMst1)))/9. + (640*pow3(lmMst2))/3.
         + ((-80*(Mgl*(55 + 6*lmMsq*(-2 + 5*lmMst1 - 5*lmMst2) + 52*lmMst2 -
        10*lmMst1*(4 + 3*lmMst2) + 30*pow2(lmMst2)) + Dmglst2*(321 + 18*lmMsq*(
        -2 + 5*lmMst1 - 5*lmMst2) + 96*lmMst2 - 30*lmMst1*(2 + 3*lmMst2) + 90*
        pow2(lmMst2)))*pow2(Mst1))/(27.*pow2(Msq)) + Dmglst2*(881.6139917695473
         + 248*B4 - 4*DN - (2560*lmMsq)/3. + lmMst1*(130.96296296296296 - 120*
        lmMsq - 40*pow2(lmMsq)) + (80*pow2(lmMsq))/3. + (176*pow2(lmMst1))/9. +
        (8*lmMst2*(4364 - 573*lmMst1 + 45*lmMsq*(5 + 6*lmMst1) + 135*pow2(
        lmMsq) + 24*pow2(lmMst1)))/27. - (8*(-377 + 90*lmMsq + 376*lmMst1)*
        pow2(lmMst2))/9. + (2944*pow3(lmMst2))/9.))/Mgl) + 10*(2552929 +
        257904*B4 - 648*DN - 456840*lmMsq + 38880*pow2(lmMsq) - 216*lmMst1*(
        4591 - 360*lmMsq + 450*pow2(lmMsq)) + 41904*pow2(lmMst1) + 216*lmMst2*(
        9211 - 6466*lmMst1 + 180*lmMsq*(-4 + 5*lmMst1) + 450*pow2(lmMsq) + 576*
        pow2(lmMst1)) - 864*(-1784 + 225*lmMsq + 840*lmMst1)*pow2(lmMst2) -
        3456*pow3(lmMst1) + 604800*pow3(lmMst2))*pow4(Mst1) + ((9*(4800*(-2*
        Mgl*(23 + lmMsq*(-6 + 9*lmMst1 - 9*lmMst2) + 18*lmMst2 - 3*lmMst1*(4 +
        3*lmMst2) + 9*pow2(lmMst2)) - 3*Dmglst2*(77 + 6*lmMsq*(-2 + 3*lmMst1 -
        3*lmMst2) + 24*lmMst2 - 6*lmMst1*(2 + 3*lmMst2) + 18*pow2(lmMst2)))*
        pow2(Msq)*pow2(Mst1) + (45*Mgl*(18201 + 1760*B4 - 16*DN - 5760*lmMsq +
        960*pow2(lmMsq) - 16*lmMst1*(199 - 30*lmMsq + 30*pow2(lmMsq)) + 16*
        lmMst2*(1291 - 322*lmMst1 + 30*lmMsq*(-5 + 2*lmMst1) + 30*pow2(lmMsq) -
        16*pow2(lmMst1)) - 256*pow2(lmMst1) - 32*(-313 + 30*lmMsq + 72*lmMst1)*
        pow2(lmMst2) + 2560*pow3(lmMst2)) + Dmglst2*(779917 + 188640*B4 - 3600*
        DN - 648000*lmMsq + 43200*pow2(lmMsq) - 720*lmMst1*(-173 + 90*lmMsq +
        30*pow2(lmMsq)) + 720*lmMst2*(1623 - 130*lmMst1 + 30*lmMsq*(-1 + 2*
        lmMst1) + 30*pow2(lmMsq) - 16*pow2(lmMst1)) + 11520*pow2(lmMst1) -
        1440*(-265 + 30*lmMsq + 136*lmMst1)*pow2(lmMst2) + 207360*pow3(lmMst2))
        )*pow4(Msq) - 300*(Mgl*(233 + 36*lmMsq*(-2 + 5*lmMst1 - 5*lmMst2) +
        207*lmMst2 - 45*lmMst1*(3 + 4*lmMst2) + 180*pow2(lmMst2)) + Dmglst2*(
        1961 + 180*lmMsq*(-2 + 5*lmMst1 - 5*lmMst2) + 675*lmMst2 - 45*lmMst1*(7
        + 20*lmMst2) + 900*pow2(lmMst2)))*pow4(Mst1))*pow4(Mst2))/pow4(Msq) -
        160*OepS2*((-3036*Dmglst2 + 816*Mgl)*pow2(Mst1)*pow2(Mst2) + 1882*Mgl*
        pow4(Mst1) - 63*(5*Dmglst2 - 3*Mgl)*pow4(Mst2)) - 108*S2*(3*Dmglst2*
        pow2(Mst2)*(8*(2489 + 3795*lmMst1 - 3795*lmMst2)*pow2(Mst1) + 9*(-453 +
        350*lmMst1 - 350*lmMst2)*pow2(Mst2)) - 5*Mgl*(36*(169 + 136*lmMst1 -
        136*lmMst2)*pow2(Mst1)*pow2(Mst2) + 2*(9185 + 5646*lmMst1 - 5646*
        lmMst2)*pow4(Mst1) + 81*(-1 + 14*lmMst1 - 14*lmMst2)*pow4(Mst2))) + (
        90*(-48*(15*Dmglst2*(95 + 12*lmMsq*(-2 + lmMst1 - lmMst2) + 32*lmMst2 -
        4*lmMst1*(2 + 3*lmMst2) + 12*pow2(lmMst2)) + 5*Mgl*(71 + 12*lmMsq*(-2 +
        lmMst1 - lmMst2) + 40*lmMst2 - 4*lmMst1*(4 + 3*lmMst2) + 12*pow2(
        lmMst2)))*pow2(Msq) - 5*(Mgl*(1153 + 12*lmMsq*(-35 + 54*lmMst1 - 54*
        lmMst2) + 906*lmMst2 - 162*lmMst1*(3 + 4*lmMst2) + 648*pow2(lmMst2)) +
        Dmglst2*(8713 + 12*lmMsq*(-179 + 270*lmMst1 - 270*lmMst2) + 3282*lmMst2
        - 162*lmMst1*(7 + 20*lmMst2) + 3240*pow2(lmMst2)))*pow2(Mst1) - 5*(Mgl*
        (911 + 12*lmMsq*(-37 + 18*lmMst1 - 18*lmMst2) + 606*lmMst2 - 54*lmMst1*
        (3 + 4*lmMst2) + 216*pow2(lmMst2)) + Dmglst2*(5591 + 12*lmMsq*(-181 +
        90*lmMst1 - 90*lmMst2) + 2550*lmMst2 - 54*lmMst1*(7 + 20*lmMst2) +
        1080*pow2(lmMst2)))*pow2(Mst2) + (2304*(1 + lmMst2)*(-Dmglst2 + 3*
        Dmglst2*lmMst2 + Mgl + lmMst2*Mgl)*pow4(Msq))/pow2(Mst1))*pow6(Mst2))/
        pow4(Msq))/Mgl))/(14580.*pow5(Mst2))))/pow2(Sbeta) + (T1ep*(4*Dmglst2*
        Mst2*(pow4(Mst1)*(-8*pow2(Mt)*pow2(s2t)*(4327*pow2(MuSUSY)*(-1 + pow2(
        Sbeta)) + 6921*pow2(Mst2)*pow2(Sbeta)) + 205840*Mst2*s2t*pow2(Sbeta)*
        pow3(Mt) - 15571*Mt*pow2(Sbeta)*pow3(Mst2)*pow3(s2t) - 89312*pow2(
        Sbeta)*pow4(Mt) + 4199*pow2(Sbeta)*pow4(Mst2)*pow4(s2t)) + 189*pow3(
        Mst2)*(-4*Mst2*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta)) + 10*
        s2t*(pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 2*pow2(Mst2)*pow2(Sbeta))*pow3(
        Mt) - 5*Mt*pow2(Sbeta)*pow3(s2t)*pow4(Mst2) - 16*Mst2*pow2(Sbeta)*pow4(
        Mt) + pow2(Sbeta)*pow4(s2t)*pow5(Mst2)) + 18*pow2(Mst1)*(-12*pow2(Mst2)
        *pow2(Mt)*pow2(s2t)*(43*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 50*pow2(Mst2)
        *pow2(Sbeta)) + 2*Mst2*s2t*(506*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 1383*
        pow2(Mst2)*pow2(Sbeta))*pow3(Mt) + 56*(3*pow2(MuSUSY)*(-1 + pow2(Sbeta)
        ) - 23*pow2(Mst2)*pow2(Sbeta))*pow4(Mt) - 401*Mt*pow2(Sbeta)*pow3(s2t)*
        pow5(Mst2) + 108*pow2(Sbeta)*pow4(s2t)*pow6(Mst2))) - 3*Mgl*(108*pow4(
        Mst2)*(Mst2*pow2(Mt)*pow2(s2t)*(-16*pow2(MuSUSY)*(-1 + pow2(Sbeta)) -
        23*pow2(Mst2)*pow2(Sbeta)) + 14*s2t*(pow2(MuSUSY)*(-1 + pow2(Sbeta)) +
        2*pow2(Mst2)*pow2(Sbeta))*pow3(Mt) - 7*Mt*pow2(Sbeta)*pow3(s2t)*pow4(
        Mst2) + 28*Mst2*pow2(Sbeta)*pow4(Mt) + 4*pow2(Sbeta)*pow4(s2t)*pow5(
        Mst2)) + pow4(Mst1)*(-4*Mst2*pow2(Mt)*pow2(s2t)*(7357*pow2(MuSUSY)*(-1
        + pow2(Sbeta)) + 6917*pow2(Mst2)*pow2(Sbeta)) + 16*s2t*(941*pow2(
        MuSUSY)*(-1 + pow2(Sbeta)) + 2740*pow2(Mst2)*pow2(Sbeta))*pow3(Mt) -
        1756*Mt*pow2(Sbeta)*pow3(s2t)*pow4(Mst2) + 18688*Mst2*pow2(Sbeta)*pow4(
        Mt) + 1945*pow2(Sbeta)*pow4(s2t)*pow5(Mst2)) + 6*Mst2*pow2(Mst1)*(-4*
        pow2(Mst2)*pow2(Mt)*pow2(s2t)*(487*pow2(MuSUSY)*(-1 + pow2(Sbeta)) +
        644*pow2(Mst2)*pow2(Sbeta)) + 8*Mst2*s2t*(136*pow2(MuSUSY)*(-1 + pow2(
        Sbeta)) + 377*pow2(Mst2)*pow2(Sbeta))*pow3(Mt) + 8*(126*pow2(MuSUSY)*(-
        1 + pow2(Sbeta)) + 209*pow2(Mst2)*pow2(Sbeta))*pow4(Mt) - 292*Mt*pow2(
        Sbeta)*pow3(s2t)*pow5(Mst2) + 343*pow2(Sbeta)*pow4(s2t)*pow6(Mst2)))))/
        (1458.*Mgl*pow2(Sbeta)*pow5(Mst2)) + (MuSUSY*(-(pow2(Mt)*pow2(s2t)*((
        pow2(Mst1)*(794.5756172839506 + 96*B4 - 120*lmMsq - 10*lmMst1*(43 - 4*
        lmMsq + 4*pow2(lmMsq)) - (56*pow2(lmMst1))/3. + (2*lmMst2*(901 - 780*
        lmMst1 + 60*lmMsq*(-1 + 2*lmMst1) + 60*pow2(lmMsq) + 32*pow2(lmMst1)))/
        3. - (16*(-101 + 15*lmMsq + 44*lmMst1)*pow2(lmMst2))/3. + (640*pow3(
        lmMst2))/3. + ((-40*(6*Dmglst2*(15 + lmMst1*(-4 + 6*lmMsq - 6*lmMst2) +
        4*lmMst2 - 6*lmMsq*lmMst2 + 6*pow2(lmMst2)) + Mgl*(9 + 4*lmMst1*(-4 +
        3*lmMsq - 3*lmMst2) + 16*lmMst2 - 12*lmMsq*lmMst2 + 12*pow2(lmMst2)))*
        pow2(Mst1))/(9.*pow2(Msq)) + Dmglst2*(600.2756172839506 + (592*B4)/3. -
        (8*DN)/3. - 680*lmMsq + 40*lmMsq*(3 + lmMsq + 2*lmMst1)*lmMst2 +
        lmMst1*(81.11111111111111 - 120*lmMsq - 40*pow2(lmMsq)) + (56*pow2(
        lmMst1))/3. + (2*lmMst2*(3859 - 756*lmMst1 + 96*pow2(lmMst1)))/9. - (
        16*(-28 + 15*lmMsq + 60*lmMst1)*pow2(lmMst2))/3. + (896*pow3(lmMst2))/
        3.))/Mgl))/Mst2 + ((1073.5195473251028 + 96*B4 - 110*lmMsq - lmMst1*(
        457.55555555555554 - 20*lmMsq + 40*pow2(lmMsq)) + (652*pow2(lmMst1))/9.
         + (2*lmMst2*(2635 - 3160*lmMst1 + 90*lmMsq*(-1 + 4*lmMst1) + 180*pow2(
        lmMsq) + 528*pow2(lmMst1)))/9. + (629.7777777777778 - 80*lmMsq - 416*
        lmMst1)*pow2(lmMst2) - (32*pow3(lmMst1))/9. - (Dmglst2*(
        618.5941700960219 - (592*B4)/3. + (8*DN)/3. + 750*lmMsq + (20*lmMst1*(
        542 + 243*lmMsq + 162*pow2(lmMsq)))/81. + (92*pow2(lmMst1))/3. -
        lmMst2*(20*lmMsq*(3 + 4*lmMst1) + 40*pow2(lmMsq) + (8*(11723 - 702*
        lmMst1 + 756*pow2(lmMst1)))/81.) + (4*(-75 + 60*lmMsq + 344*lmMst1)*
        pow2(lmMst2))/3. - (32*pow3(lmMst1))/3. - (1120*pow3(lmMst2))/3.))/Mgl
        + (2720*pow3(lmMst2))/9.)*pow4(Mst1))/pow3(Mst2) - S2*((4.5 - 63*lmMst1
        + 63*lmMst2 + (3*Dmglst2*(-453 + 350*lmMst1 - 350*lmMst2))/(10.*Mgl))*
        Mst2 - ((342.5 + 209*lmMst1 - 209*lmMst2 - (Dmglst2*(799.6333333333333
         + 907*lmMst1 - 907*lmMst2))/Mgl)*pow2(Mst1))/Mst2 + ((Dmglst2*(525961
        + 356010*lmMst1 - 356010*lmMst2) - 15*(6143 + 3198*lmMst1 - 3198*
        lmMst2)*Mgl)*pow4(Mst1))/(135.*Mgl*pow3(Mst2))) + Mst2*(
        737.0416666666666 + (220*B4)/3. - (2*DN)/3. - 240*lmMsq + 40*pow2(
        lmMsq) - (2*lmMst1*(199 - 30*lmMsq + 30*pow2(lmMsq)))/3. + (2*lmMst2*(
        1227 - 322*lmMst1 + 30*lmMsq*(-5 + 2*lmMst1) + 30*pow2(lmMsq) - 16*
        pow2(lmMst1)))/3. - (32*pow2(lmMst1))/3. - 4*(-99 + 10*lmMsq + 24*
        lmMst1)*pow2(lmMst2) + (320*pow3(lmMst2))/3. + ((-20*(3*Dmglst2*(59 +
        8*lmMst1*(-2 + 3*lmMsq - 3*lmMst2) + 16*lmMst2 - 24*lmMsq*lmMst2 + 24*
        pow2(lmMst2)) + Mgl*(21 + 8*lmMst1*(-4 + 3*lmMsq - 3*lmMst2) + 32*
        lmMst2 - 24*lmMsq*lmMst2 + 24*pow2(lmMst2)))*pow2(Mst1))/(9.*pow2(Msq))
        + Dmglst2*(743.4787037037037 + (524*B4)/3. - (10*DN)/3. - 600*lmMsq +
        lmMst1*(115.33333333333333 - 60*lmMsq - 20*pow2(lmMsq)) + 40*pow2(
        lmMsq) + (2*lmMst2*(1559 - 130*lmMst1 + 30*lmMsq*(-1 + 2*lmMst1) + 30*
        pow2(lmMsq) - 16*pow2(lmMst1)))/3. + (32*pow2(lmMst1))/3. - (4*(-217 +
        30*lmMsq + 136*lmMst1)*pow2(lmMst2))/3. + 192*pow3(lmMst2)))/Mgl - ((
        11.342592592592593 - (5*lmMsq)/9. + 5*(-3 + 4*lmMsq)*lmMst1 + (20*(7 -
        9*(lmMsq + lmMst1))*lmMst2)/9. + 20*pow2(lmMst2) + (Dmglst2*(
        141.34259259259258 - (5*lmMsq)/9. + 5*(-7 + 20*lmMsq)*lmMst1 + (
        35.55555555555556 - 100*(lmMsq + lmMst1))*lmMst2 + 100*pow2(lmMst2)))/
        Mgl)*pow4(Mst1))/pow4(Msq)) + (pow3(Mst2)*((4*((-5*(3*Dmglst2*(95 + 12*
        lmMsq*(-2 + lmMst1 - lmMst2) + 32*lmMst2 - 4*lmMst1*(2 + 3*lmMst2) +
        12*pow2(lmMst2)) + Mgl*(71 + 12*lmMsq*(-2 + lmMst1 - lmMst2) + 40*
        lmMst2 - 4*lmMst1*(4 + 3*lmMst2) + 12*pow2(lmMst2))))/pow2(Msq) + (48*(
        1 + lmMst2)*(-Dmglst2 + 3*Dmglst2*lmMst2 + Mgl + lmMst2*Mgl))/pow2(
        Mst1)))/9. - (5*((Mgl*(242 + 24*lmMsq*(1 + 18*lmMst1 - 18*lmMst2) +
        300*lmMst2 - 108*lmMst1*(3 + 4*lmMst2) + 432*pow2(lmMst2)) + Dmglst2*(
        3122 + 24*lmMsq*(1 + 90*lmMst1 - 90*lmMst2) + 732*lmMst2 - 108*lmMst1*(
        7 + 20*lmMst2) + 2160*pow2(lmMst2)))*pow2(Mst1) + (Mgl*(911 + 12*lmMsq*
        (-37 + 18*lmMst1 - 18*lmMst2) + 606*lmMst2 - 54*lmMst1*(3 + 4*lmMst2) +
        216*pow2(lmMst2)) + Dmglst2*(5591 + 12*lmMsq*(-181 + 90*lmMst1 - 90*
        lmMst2) + 2550*lmMst2 - 54*lmMst1*(7 + 20*lmMst2) + 1080*pow2(lmMst2)))
        *pow2(Mst2)))/(108.*pow4(Msq))) + (4*OepS2*(-3*Mgl*(627*pow2(Mst1)*
        pow2(Mst2) + 1066*pow4(Mst1) + 189*pow4(Mst2)) + Dmglst2*(8163*pow2(
        Mst1)*pow2(Mst2) + 23734*pow4(Mst1) + 945*pow4(Mst2))))/(729.*pow3(
        Mst2)))/Mgl)) - s2t*pow3(Mt)*(454.18086419753087 + (16*D3)/9. - (16*DN)
        /9. - (100*lmMsq)/3. + (64*lmMt)/3. + 20*pow2(lmMsq) + (4*lmMst1*(667 -
        300*lmMsq + 45*pow2(lmMsq)))/27. + (4*lmMst2*(1946 + lmMsq*(30 - 90*
        lmMst1) + 294*lmMst1 - 45*pow2(lmMsq) - 117*pow2(lmMst1)))/27. + (248*
        pow2(lmMst1))/9. + (4*(152 + 30*lmMsq + 31*lmMst1)*pow2(lmMst2))/9. + (
        64*pow3(lmMst1))/9. - (32*pow3(lmMst2))/9. + (4*(-(((4500*Dmglst2 +
        Mgl*(9608 + 8520*lmMst1 - 15*lmMsq*(218 + 195*lmMst1 - 75*lmMst2) +
        375*(-14 + 15*lmMst1)*lmMst2 + 900*pow2(lmMsq) - 1350*pow2(lmMst1) -
        3375*pow2(lmMst2)))*pow2(Mst1))/pow2(Msq)) + 30*Dmglst2*(15529 - 120*B4
        + 120*D3 - 60*DN + 150*lmMsq + 690*lmMst1 + 675*pow2(lmMsq) + 345*pow2(
        lmMst1) - 5*lmMst2*(-3091 + 270*lmMsq + 42*lmMst1 + 96*pow2(lmMst1)) +
        (4305 - 480*lmMst1)*pow2(lmMst2) + 960*pow3(lmMst2))))/(2025.*Mgl) -
        pow4(Mst1)*((3.340708346830796 - (2*lmMsq*(256 + 259*lmMst1 - 49*
        lmMst2))/441. - (212*lmMst2)/135. + lmMst1*(2.731368102796674 + (38*
        lmMst2)/9.) + (20*Dmglst2)/(9.*Mgl) + (10*pow2(lmMsq))/21. - (32*pow2(
        lmMst1))/21. - (20*pow2(lmMst2))/9.)/pow4(Msq) - (1188.1718576027706 -
        (80*lmMsq)/3. + (-960.3533282942807 + 40*lmMsq)*lmMst1 + lmMst2*(
        1415.4644394053919 - (138308*lmMst1)/945. - (40*lmMsq*(3 + 2*lmMst1))/
        3. - (820*pow2(lmMst1))/9.) + (2*(-94783 + 6300*lmMsq)*pow2(lmMst1))/
        945. + (432.2899470899471 + (40*lmMsq)/3. - (332*lmMst1)/9.)*pow2(
        lmMst2) + (64*lmMt*(1 + 4*lmMst2 - 2*lmMst1*(2 + lmMst2) + pow2(lmMst1)
        + pow2(lmMst2)))/3. + (52*pow3(lmMst1))/27. + (3404*pow3(lmMst2))/27. +
        (Dmglst2*(772.4522050654477 - (64*B4)/9. + (64*D3)/9. - (32*DN)/9. + (
        560*lmMsq)/3. + (1513.2107835726883 - 160*lmMsq)*lmMst1 + (38704*pow2(
        lmMst1))/189. + lmMst2*(160*lmMsq + (4*(1947623 - 3019800*lmMst1 +
        97020*pow2(lmMst1)))/19845.) + (16*(6787 - 5439*lmMst1)*pow2(lmMst2))/
        189. - (128*lmMt*(4 + 6*lmMst2 - 2*lmMst1*(3 + lmMst2) + pow2(lmMst1) +
        pow2(lmMst2)))/3. - (16*pow3(lmMst1))/9. + (1328*pow3(lmMst2))/3.))/
        Mgl)/pow4(Mst2)) + ((-2*(2*Dmglst2*(-179 + 2*(87 + 32*lmMst1)*lmMst2 -
        30*lmMsq*(-1 + 6*lmMst2) + 90*pow2(lmMsq) + 208*pow2(lmMst2)) + Mgl*(
        277 + 514*lmMst2 + 64*lmMst1*(1 + lmMst2) - 30*lmMsq*(5 + 6*lmMst2) +
        90*pow2(lmMsq) + 208*pow2(lmMst2)))*pow2(Mst2))/(9.*pow2(Mst1)) + (
        pow2(Mst1)*(16*Dmglst2*(9080791 - 54000*B4 + 54000*D3 - 27000*DN +
        607500*lmMsq - 60*(-42359 + 3375*lmMsq)*lmMst1 + 324000*(-2 + lmMst1 -
        lmMst2)*lmMt + 60*lmMst2*(143041 + 3375*lmMsq - 14460*lmMst1 - 3375*
        pow2(lmMst1)) + 1800*pow2(lmMst1) - 900*(-2402 + 1695*lmMst1)*pow2(
        lmMst2) - 4500*pow3(lmMst1) + 1732500*pow3(lmMst2)) + 5*Mgl*(19319827 +
        19042320*lmMst2 + 518400*(1 + 2*lmMst2)*lmMt - 3600*(701 + 186*lmMst2)*
        pow2(lmMst1) + 5943600*pow2(lmMst2) + 324000*lmMsq*(-2 + lmMst1 -
        lmMst2 - 2*lmMst1*lmMst2 + pow2(lmMst1) + pow2(lmMst2)) - 240*lmMst1*(
        35693 + 5610*lmMst2 + 4320*lmMt + 3690*pow2(lmMst2)) - 295200*pow3(
        lmMst1) + 1850400*pow3(lmMst2))))/(121500.*pow2(Mst2)) + (pow2(Mst2)*(-
        17150*Dmglst2*(3*(7 + 12*lmMsq - 12*lmMst2)*pow2(Mst1) + 2*(43 - 144*
        lmMsq + 144*lmMst2)*pow2(Mst2)) + Mgl*(1715*(75 + 60*lmMsq*(-3 + 2*
        lmMst1 - 2*lmMst2) + 356*lmMst2 - 8*lmMst1*(22 + 15*lmMst2) + 120*pow2(
        lmMst2))*pow2(Mst1) + (-393767 + 420*lmMsq*(2194 + 49*lmMst1 - 259*
        lmMst2) - 878948*lmMst2 - 1372*lmMst1*(31 + 15*lmMst2) + 44100*pow2(
        lmMsq) + 64680*pow2(lmMst2))*pow2(Mst2))))/(92610.*pow4(Msq)) + ((5*((
        64*Dmglst2*(14 - 15*lmMsq + 15*lmMst2) + 8*(43 - 30*lmMsq + 30*lmMst2)*
        Mgl)*pow2(Msq) + (18*Dmglst2*(13 - 28*lmMsq + 28*lmMst2) + (67 - 84*
        lmMsq + 84*lmMst2)*Mgl)*pow2(Mst2)))/(54.*pow2(Mst1)*pow4(Msq)) + (64*(
        1 + lmMst2)*(4*Dmglst2*lmMst2 + Mgl + lmMst2*Mgl))/(9.*pow4(Mst1)))*
        pow4(Mst2) + (-9000*Dmglst2*(2*(-3 + lmMst1 - lmMst2)*pow4(Mst1) + (17
        - 20*lmMsq + 20*lmMst2)*pow4(Mst2)) + 4*Mgl*(1125*(-4 + lmMst1*(2 - 6*
        lmMst2) - 2*lmMst2 + 3*pow2(lmMst1) + 3*pow2(lmMst2))*pow4(Mst1) + (-
        6517 + 15*lmMsq*(532 + 75*lmMst1 - 195*lmMst2) - 4980*lmMst2 - 375*
        lmMst1*(8 + 3*lmMst2) + 900*pow2(lmMsq) + 2025*pow2(lmMst2))*pow4(Mst2)
        ))/(2025.*pow2(Msq)*pow2(Mst2)) + ((8*OepS2*(24*Dmglst2*(150*pow2(Mst1)
        *pow2(Mst2) + 919*pow4(Mst1)) - Mgl*(4485*pow2(Mst1)*pow2(Mst2) +
        11402*pow4(Mst1) + 621*pow4(Mst2))))/729. - (S2*(24*Dmglst2*(45*(139 +
        100*lmMst1 - 100*lmMst2)*pow2(Mst1)*pow2(Mst2) + (34538 + 27570*lmMst1
        - 27570*lmMst2)*pow4(Mst1) + 5697*pow4(Mst2)) - Mgl*(9*(11113 + 14950*
        lmMst1 - 14950*lmMst2)*pow2(Mst1)*pow2(Mst2) + 4*(126661 + 85515*lmMst1
        - 85515*lmMst2)*pow4(Mst1) + 81*(-891 + 230*lmMst1 - 230*lmMst2)*pow4(
        Mst2))))/135.)/pow4(Mst2))/Mgl) + (pow4(Mt)*(20*(567143 - 125172*lmMst2
        + 648*lmMst2*(-38 + 5*lmMst2)*lmMt - 9720*(lmMsq*(lmMst1*(9 - 2*lmMt) +
        2*lmMst2*(-6 + lmMt) + 3*(-4 + lmMt)) + lmMt) - 9720*(lmMst1 - lmMst2)*
        pow2(lmMsq) + 240300*pow2(lmMst2) - 108*lmMst1*(2111 + lmMst1*(457 -
        12*lmMst2 - 126*lmMt) - 420*lmMt + 4*lmMst2*(454 + 39*lmMt) + 750*pow2(
        lmMst2) - 288*pow2(lmMt)) - 15552*(1 + 2*lmMst2)*pow2(lmMt) + 216*pow3(
        lmMst1) + 79488*pow3(lmMst2))*pow4(Mst1) + ((-24*pow2(Mst1)*pow2(Mst2)*
        (2700*Mgl*(-25 + lmMst1*(-5 + 6*lmMsq - 6*lmMst2) + (4 - 6*lmMsq)*
        lmMst2 + lmMt + 6*pow2(lmMst2))*pow2(Mst1) + 300*Dmglst2*(-313 + 3*
        lmMst1*(-23 + 42*lmMsq - 42*lmMst2) - 42*(-2 + 3*lmMsq)*lmMst2 - 15*
        lmMt + 126*pow2(lmMst2))*pow2(Mst1) + 2*Mgl*pow2(Msq)*(-71744 + 28755*
        lmMst2 + 15255*lmMt + 17820*lmMst2*lmMt + 4050*(lmMst1 - lmMst2)*pow2(
        lmMsq) - 540*(-15 + 4*lmMt)*pow2(lmMst1) - 89910*pow2(lmMst2) + 6480*
        lmMt*pow2(lmMst2) + 4050*lmMsq*(-11 + lmMst1 - 2*lmMst2 - 2*lmMst1*
        lmMst2 + lmMt + 2*pow2(lmMst2)) + 270*lmMst1*(221 + 279*lmMst2 - 66*
        lmMt - 16*lmMst2*lmMt + 112*pow2(lmMst2) - 24*pow2(lmMt)) + 6480*pow2(
        lmMt) + 6480*lmMst2*pow2(lmMt) - 30240*pow3(lmMst2)) + 3*Dmglst2*pow2(
        Msq)*(93487 + 196530*lmMst2 + 9030*lmMt - 7380*lmMst2*lmMt + 2700*(
        lmMst1 - lmMst2)*pow2(lmMsq) - 360*(23 + 8*lmMst2 - 4*lmMt)*pow2(
        lmMst1) + 68040*pow2(lmMst2) + 1440*lmMt*pow2(lmMst2) + 300*lmMsq*(29 +
        lmMst1*(9 - 18*lmMst2) + 6*lmMst2 - 15*lmMt + 18*pow2(lmMst2)) - 8640*
        pow2(lmMt) - 4320*lmMst2*pow2(lmMt) + 20*lmMst1*(-5846 + 321*lmMt - 12*
        lmMst2*(209 + 12*lmMt) + 576*pow2(lmMst2) + 216*pow2(lmMt)) - 8640*
        pow3(lmMst2))))/pow2(Msq) + (9*(9*Mgl*(800*(19 + (-8 + 6*lmMsq)*lmMst2
        + lmMst1*(7 - 6*lmMsq + 6*lmMst2) + lmMt - 6*pow2(lmMst2))*pow2(Msq)*
        pow2(Mst1) + (813 - 3360*lmMst2 - 2400*(lmMst1 - lmMst2)*pow2(lmMsq) -
        1280*(1 + lmMst2)*pow2(lmMst1) + 2400*lmMsq*(7 + lmMst1 - 2*lmMst2 + 2*
        lmMst1*lmMst2 + lmMt - 2*pow2(lmMst2)) + 30880*pow2(lmMst2) - 80*lmMt*(
        175 + 96*lmMst2 + 32*pow2(lmMst2)) - 80*lmMst1*(167 + 258*lmMst2 - 32*(
        1 + lmMst2)*lmMt + 112*pow2(lmMst2)) - 3840*pow2(lmMt) + 10240*pow3(
        lmMst2))*pow4(Msq) + 200*(38 + 4*(-2 + 3*lmMsq)*lmMst2 + lmMst1*(7 -
        12*lmMsq + 12*lmMst2) + lmMt - 12*pow2(lmMst2))*pow4(Mst1)) - 5*
        Dmglst2*(-160*(439 + 6*(-8 + 21*lmMsq)*lmMst2 + 3*lmMst1*(17 - 42*lmMsq
        + 42*lmMst2) - 3*lmMt - 126*pow2(lmMst2))*pow2(Msq)*pow2(Mst1) + (7335
        + 169856*lmMst2 + 4320*(lmMst1 - lmMst2)*pow2(lmMsq) + 2304*(-1 +
        lmMst2)*pow2(lmMst1) + 21984*pow2(lmMst2) + 480*lmMsq*(-79 - 30*lmMst2
        - 9*lmMst1*(-3 + 2*lmMst2) + 3*lmMt + 18*pow2(lmMst2)) + 48*lmMt*(31 -
        128*lmMst2 + 96*pow2(lmMst2)) + 144*lmMst1*(-269 + 32*lmMt - 2*lmMst2*(
        63 + 16*lmMt) + 112*pow2(lmMst2)) - 18432*pow3(lmMst2))*pow4(Msq) + 40*
        (-1330 + lmMst1*(-51 + 396*lmMsq - 396*lmMst2) + (48 - 396*lmMsq)*
        lmMst2 + 3*lmMt + 396*pow2(lmMst2))*pow4(Mst1)))*pow4(Mst2))/pow4(Msq)
        + 108*S2*((-48*Dmglst2*(601 + 2790*lmMst1 - 2790*lmMst2) + 72*(457 +
        550*lmMst1 - 550*lmMst2)*Mgl)*pow2(Mst1)*pow2(Mst2) + 40*(4082 + 3045*
        lmMst1 - 3045*lmMst2)*Mgl*pow4(Mst1) + (675*Dmglst2*(27 - 14*lmMst1 +
        14*lmMst2) + 81*(-79 + 70*lmMst1 - 70*lmMst2)*Mgl)*pow4(Mst2)) + 160*
        OepS2*(9*Dmglst2*pow2(Mst2)*(496*pow2(Mst1) + 35*pow2(Mst2)) - Mgl*(
        1320*pow2(Mst1)*pow2(Mst2) + 4060*pow4(Mst1) + 189*pow4(Mst2))) + (90*(
        -32*(-5*Mgl*(58 - 36*lmMst2 + 6*lmMst1*(4 + 3*lmMst2) - 3*lmMt - 9*
        lmMst2*lmMt + 3*lmMsq*(5 - 6*lmMst1 + 9*lmMst2 + 3*lmMt) - 9*pow2(
        lmMsq) - 18*pow2(lmMst2)) - 15*Dmglst2*(64 - 20*lmMst2 + 6*lmMst1*(2 +
        3*lmMst2) + 3*lmMt - 3*lmMst2*lmMt + lmMsq*(5 - 18*lmMst1 + 21*lmMst2 +
        3*lmMt) - 3*pow2(lmMsq) - 18*pow2(lmMst2)))*pow2(Msq) + 5*(3*Dmglst2*(
        1693 - 136*lmMst2 + 6*lmMst1*(19 + 84*lmMst2) - 4*lmMsq*(5 + 126*lmMst1
        - 132*lmMst2 - 6*lmMt) + 42*lmMt - 24*lmMst2*lmMt - 24*pow2(lmMsq) -
        504*pow2(lmMst2)) + Mgl*(1271 - 216*lmMst2 + 6*lmMst1*(41 + 60*lmMst2)
        - 12*lmMsq*(5 + 30*lmMst1 - 36*lmMst2 - 6*lmMt) + 30*lmMt - 72*lmMst2*
        lmMt - 72*pow2(lmMsq) - 360*pow2(lmMst2)))*pow2(Mst1) + 15*(Dmglst2*(
        1535 - 184*lmMst2 + 18*lmMst1*(7 + 20*lmMst2) - 6*(-9 + 28*lmMst2)*lmMt
        + 4*lmMsq*(1 - 90*lmMst1 + 132*lmMst2 + 42*lmMt) - 168*pow2(lmMsq) -
        360*pow2(lmMst2)) + Mgl*(403 - 24*lmMst2 + 18*lmMst1*(3 + 4*lmMst2) -
        12*lmMsq*(1 + 6*lmMst1 - 12*lmMst2 - 6*lmMt) - 18*(1 + 4*lmMst2)*lmMt -
        72*pow2(lmMsq) - 72*pow2(lmMst2)))*pow2(Mst2) + (2304*(1 + lmMst2)*(-
        Dmglst2 + 3*Dmglst2*lmMst2 + Mgl + lmMst2*Mgl)*pow4(Msq))/pow2(Mst1))*
        pow6(Mst2))/pow4(Msq))/Mgl))/(7290.*pow5(Mst2)) + Mt*(pow3(s2t)*(pow2(
        Mst1)*(301.91041152263375 + 4*B4 + (4*D3)/9. - (4*DN)/9. - (335*lmMsq)/
        9. + (25*pow2(lmMsq))/3. + lmMst1*(188.5172839506173 - (325*lmMsq)/9. +
        (5*pow2(lmMsq))/3.) + (10*(-149 + 18*lmMsq)*pow2(lmMst1))/27. - (
        lmMst2*(53429 - 183450*lmMst1 + 2250*lmMsq*(-7 + 6*lmMst1) + 1350*pow2(
        lmMsq) + 29610*pow2(lmMst1)))/810. + ((-3683 + 270*lmMsq + 2595*lmMst1)
        *pow2(lmMst2))/27. - (110*pow3(lmMst1))/27. - (1498*pow3(lmMst2))/27. -
        (Dmglst2*(5430043 + 1334580*lmMst2 + 202500*lmMsq*(7 + 6*lmMst2) -
        607500*pow2(lmMsq) + 900*(-859 + 3690*lmMst2)*pow2(lmMst1) + 846900*
        pow2(lmMst2) - 120*lmMst1*(44309 + 12315*lmMst2 + 56475*pow2(lmMst2)) +
        45000*pow3(lmMst1) + 3411000*pow3(lmMst2)))/(60750.*Mgl)) - (((150*
        Dmglst2*(5 - 2*lmMst1 + 2*lmMst2) + Mgl*(1503 + 1670*lmMst1 - 60*lmMsq*
        (12 + 5*lmMst1 - 5*lmMst2) + 50*(-19 + 42*lmMst1)*lmMst2 - 900*pow2(
        lmMst1) - 1200*pow2(lmMst2)))/(270.*pow2(Msq)) + ((5*lmMsq*(1 - 2*
        lmMst1 + 2*lmMst2)*(12*Dmglst2 + Mgl + 8*lmMst1*Mgl - 8*lmMst2*Mgl))/6.
         + (-3*Mgl*(79604910131 + 34248010860*lmMst2 - 1058400*(78238 + 33285*
        lmMst2)*pow2(lmMst1) - 83251627200*pow2(lmMst2) + 1260*lmMst1*(-
        29738761 + 131792640*lmMst2 + 59005800*pow2(lmMst2)) - 1296540000*pow3(
        lmMst1) - 37821924000*pow3(lmMst2)) + 20*Dmglst2*(32030812049 +
        1902812940*lmMst2 + 793800*(433 + 6552*lmMst2)*pow2(lmMst1) + 76998600*
        pow2(lmMst2) - 1260*lmMst1*(2674409 + 333900*lmMst2 + 7514640*pow2(
        lmMst2)) - 311169600*pow3(lmMst1) + 4578638400*pow3(lmMst2)))/1.500282e9
        )/pow2(Mst2))*pow4(Mst1))/Mgl - (pow2(Mst2)*(12760 - 720*B4 + 72*D3 -
        36*DN - 1980*lmMsq + 540*pow2(lmMsq) - 3*lmMst1*(5471 - 1950*lmMsq +
        630*pow2(lmMsq)) - 954*pow2(lmMst1) + 3*lmMst2*(8495 - 6750*lmMst1 +
        210*lmMsq*(-11 + 6*lmMst1) + 630*pow2(lmMsq) + 234*pow2(lmMst1)) - 18*(
        -1324 + 210*lmMsq + 399*lmMst1)*pow2(lmMst2) - 288*pow3(lmMst1) + 6768*
        pow3(lmMst2) + (162*(((Mgl*(-725 + 30*lmMsq*(5 + lmMst1 - lmMst2) -
        793*lmMst2 + lmMst1*(643 + 150*lmMst2) - 90*pow2(lmMst1) - 60*pow2(
        lmMst2)) + 50*Dmglst2*(20 + 6*lmMsq*(-4 + 3*lmMst1 - 3*lmMst2) + 27*
        lmMst2 - 3*lmMst1*(1 + 6*lmMst2) + 18*pow2(lmMst2)))*pow2(Mst1))/(135.*
        pow2(Msq)) + Dmglst2*(32.51234567901235 - (8*B4)/3. + (32*D3)/9. - (20*
        DN)/9. + (160*lmMsq)/3. + lmMst1*(9.333333333333334 - 20*pow2(lmMsq)) -
        (26*pow2(lmMst1))/3. + (4*lmMst2*(-237 - 248*lmMst1 + 90*lmMsq*lmMst1 +
        45*pow2(lmMsq) + 16*pow2(lmMst1)))/9. - (2*(-439 + 180*lmMsq + 336*
        lmMst1)*pow2(lmMst2))/9. + (608*pow3(lmMst2))/9.)))/Mgl + (162*(
        0.32026967930029154 + lmMsq*(-0.6291383219954648 - (43*lmMst1)/63. + (
        7*lmMst2)/9.) - (553*lmMst2)/540. + lmMst1*(1.6532123960695388 + (13*
        lmMst2)/9.) - (5*Dmglst2)/(18.*Mgl) - pow2(lmMsq)/21. - (8*pow2(lmMst1)
        )/21. - (10*pow2(lmMst2))/9.)*pow4(Mst1))/pow4(Msq)))/162. + ((1 - 2*
        lmMst2)*shiftst3*pow4(Mst2))/(3.*pow2(Mst1)) + (((Mgl*(9*(3653969 -
        140*lmMsq*(12899 + 5390*lmMst1 - 6230*lmMst2) + 2069970*lmMst2 + 37730*
        lmMst1*(-7 + 20*lmMst2) - 58800*pow2(lmMsq) - 813400*pow2(lmMst2))*
        pow2(Mst1) - (58456697 - 420*lmMsq*(63929 + 69090*lmMst1 - 69930*
        lmMst2) + 374010*lmMst2 + 10290*lmMst1*(2573 + 2820*lmMst2) - 176400*
        pow2(lmMsq) - 29194200*pow2(lmMst2))*pow2(Mst2)) + 514500*Dmglst2*((119
        - 34*lmMst2 + 12*lmMsq*(1 - 10*lmMst1 + 10*lmMst2) + 2*lmMst1*(11 + 60*
        lmMst2) - 120*pow2(lmMst2))*pow2(Mst1) - 2*(337 - 24*lmMsq*(5 + 8*
        lmMst1 - 8*lmMst2) + 14*lmMst2 + 2*lmMst1*(53 + 96*lmMst2) - 192*pow2(
        lmMst2))*pow2(Mst2)))*pow4(Mst2))/(1.11132e7*pow4(Msq)) + ((-30*
        Dmglst2*pow2(Msq)*((-115 + 430*lmMst2 + 64*lmMst1*lmMst2 - 30*lmMsq*(-1
        + 6*lmMst2) + 90*pow2(lmMsq) + 400*pow2(lmMst2))*pow2(Mst1) - 64*
        lmMst2*(1 + lmMst2)*pow2(Mst2)) + 100*Dmglst2*pow2(Mst1)*((-174 - 109*
        lmMst1 + 6*lmMsq*(6 + 17*lmMst1 - 17*lmMst2) + 73*lmMst2 - 102*lmMst1*
        lmMst2 + 102*pow2(lmMst2))*pow2(Mst1) + 4*(14 - 15*lmMsq + 15*lmMst2)*
        pow2(Mst2)) + Mgl*(50*(43 - 30*lmMsq + 30*lmMst2)*pow2(Mst1)*pow2(Mst2)
        - 15*pow2(Msq)*((373 + 706*lmMst2 + 64*lmMst1*(1 + lmMst2) - 30*lmMsq*(
        5 + 6*lmMst2) + 90*pow2(lmMsq) + 304*pow2(lmMst2))*pow2(Mst1) - 32*
        pow2(1 + lmMst2)*pow2(Mst2)) + 6*(-637 + 20*lmMsq*(9 + 20*lmMst1 - 20*
        lmMst2) + 495*lmMst2 - 25*lmMst1*(27 + 16*lmMst2) + 400*pow2(lmMst2))*
        pow4(Mst1)))*pow4(Mst2))/(270.*pow2(Msq)*pow4(Mst1)) + (-45*Mgl*(664*
        OepS2 - 81*(1261 + 166*lmMst1 - 166*lmMst2)*S2)*pow2(Mst1)*pow2(Mst2) +
        54*Dmglst2*(632*OepS2 + 9*(16193 - 1422*lmMst1 + 1422*lmMst2)*S2)*pow2(
        Mst1)*pow2(Mst2) - 15*Mgl*(3548*OepS2 - 27*(15148 + 2661*lmMst1 - 2661*
        lmMst2)*S2)*pow4(Mst1) + 4*Dmglst2*(25328*OepS2 + 27*(47051 - 18996*
        lmMst1 + 18996*lmMst2)*S2)*pow4(Mst1) + 3024*Dmglst2*OepS2*pow4(Mst2) -
        5184*Mgl*OepS2*pow4(Mst2))/(2187.*pow2(Mst2)))/Mgl) + ((s2t*pow2(Mst2)*
        ((4*(-1 + 2*lmMst2)*Mgl*shiftst3*pow2(Mt))/pow2(Mst1) + 6*Dmglst2*(1677
        - 14*lmMst1 + 14*lmMst2)*S2*pow2(s2t) + Mgl*(144*(22 + lmMst1 - lmMst2)
        *S2 + (1 + lmMst1 - lmMst2)*(-1 + 2*lmMst2)*shiftst3)*pow2(s2t)))/3. +
        (2*T1ep*(3*Mgl*(3*pow2(Mst1)*pow2(Mst2)*(-2990*Mst2*s2t*pow2(Mt) - 627*
        Mt*pow2(Mst2)*pow2(s2t) + 1760*pow3(Mt) + 830*pow3(Mst2)*pow3(s2t)) + (
        -22804*Mst2*s2t*pow2(Mt) - 3198*Mt*pow2(Mst2)*pow2(s2t) + 16240*pow3(
        Mt) + 4435*pow3(Mst2)*pow3(s2t))*pow4(Mst1) + 27*(-46*Mst2*s2t*pow2(Mt)
        - 21*Mt*pow2(Mst2)*pow2(s2t) + 28*pow3(Mt) + 16*pow3(Mst2)*pow3(s2t))*
        pow4(Mst2)) + Dmglst2*Mst2*(-27*Mst2*pow2(Mst1)*(-800*Mst2*s2t*pow2(Mt)
        - 907*Mt*pow2(Mst2)*pow2(s2t) + 1984*pow3(Mt) + 316*pow3(Mst2)*pow3(
        s2t)) + 2*s2t*(35601*Mst2*Mt*s2t + 66168*pow2(Mt) - 12664*pow2(Mst2)*
        pow2(s2t))*pow4(Mst1) - 189*(20*pow3(Mst2)*pow3(Mt) - 15*Mt*pow2(s2t)*
        pow5(Mst2) + 4*pow3(s2t)*pow6(Mst2)))))/(729.*pow5(Mst2)))/Mgl)))/Tbeta
        + (s2t*pow3(Mt)*(18*pow2(Msq)*(-21600*(Dmglst2*(-15 + 2*lmMst1 - 4*
        lmMst2 + 2*lmMt) + (7 - 2*lmMst1 + 4*lmMst2 - 2*lmMt)*Mgl)*pow2(Mst1) +
        pow2(Msq)*(3*Dmglst2*(678821 + 737120*lmMst2 + 21600*lmMsq*(13 - 2*
        lmMst1 + 4*lmMst2 - 2*lmMt) - 2880*(21 + 12*lmMst2 - 4*lmMt)*pow2(
        lmMst1) + 440160*pow2(lmMst2) - 480*lmMt*(-279 - 13*lmMst2 + 24*pow2(
        lmMst2)) - 34560*(2 + lmMst2)*pow2(lmMt) + 160*lmMst1*(-4724 - 1941*
        lmMst2 - 39*lmMt + 72*pow2(lmMst2) + 216*pow2(lmMt)) + 23040*pow3(
        lmMst2)) - Mgl*(1058993 - 330480*lmMst2 + 64800*lmMsq*(5 - 2*lmMst1 +
        4*lmMst2 - 2*lmMt) + 8640*(-13 + 4*lmMst2 + 4*lmMt)*pow2(lmMst1) +
        622080*pow2(lmMst2) - 4320*lmMt*(-31 - 6*lmMst2 + 8*pow2(lmMst2)) -
        2160*lmMst1*(325 + 300*lmMst2 - 52*lmMt + 112*pow2(lmMst2) - 48*pow2(
        lmMt)) - 103680*lmMst2*pow2(lmMt) + 207360*pow3(lmMst2))))*pow4(Mst1)*
        pow4(Mst2) + 8*pow2(Mst2)*(Dmglst2*(24208447 + 1263240*lmMst2 - 9720*(
        266 + 74*lmMst2 - 55*lmMt)*pow2(lmMst1) - 1140480*pow2(lmMst2) +
        145800*lmMsq*(28 + lmMst2*(41 - 6*lmMt) - 10*lmMt + lmMst1*(-31 - 6*
        lmMst2 + 6*lmMt) + 6*pow2(lmMst2)) + 3240*lmMt*(846 + 1042*lmMst2 +
        309*pow2(lmMst2)) - 233280*(2 + 3*lmMst2)*pow2(lmMt) + 360*lmMst1*(-
        21491 + lmMst2*(8190 - 4266*lmMt) - 5922*lmMt + 3267*pow2(lmMst2) +
        1944*pow2(lmMt)) + 9720*pow3(lmMst1) - 466560*pow3(lmMst2)) - 3*Mgl*(
        1907557 - 326160*lmMst2 + 1080*(-158 + 6*lmMst2 + 39*lmMt)*pow2(lmMst1)
        + 75600*pow2(lmMst2) + 48600*lmMsq*(2 + 11*lmMst2 - lmMst1*(9 + 2*
        lmMst2 - 2*lmMt) - 2*(1 + lmMst2)*lmMt + 2*pow2(lmMst2)) + 1620*lmMt*(
        83 + 104*lmMst2 + 58*pow2(lmMst2)) - 540*lmMst1*(797 + 120*lmMt + 4*
        lmMst2*(4 + 63*lmMt) + 78*pow2(lmMst2) - 144*pow2(lmMt)) - 77760*
        lmMst2*pow2(lmMt) + 1080*pow3(lmMst1) + 34560*pow3(lmMst2)))*pow4(Msq)*
        pow6(Mst1) + 54*pow2(Mst1)*(5*Dmglst2*(-80*(131 + 126*lmMst2 + 6*
        lmMst1*(-7 + 6*lmMst2) + 6*(-7 + 6*lmMst2)*lmMt - 6*lmMsq*(7 + 6*(
        lmMst1 + lmMt)) + 36*pow2(lmMsq) - 36*pow2(lmMst2))*pow2(Msq)*pow2(
        Mst1) + (1983 + 182144*lmMst2 + 4320*(lmMst1 - lmMst2)*pow2(lmMsq) +
        2304*(-1 + lmMst2)*pow2(lmMst1) + 33504*pow2(lmMst2) + 240*lmMsq*(-161
        + lmMst1*(54 - 36*lmMst2) - 60*lmMst2 + 6*lmMt + 36*pow2(lmMst2)) + 48*
        lmMt*(31 - 128*lmMst2 + 96*pow2(lmMst2)) + 144*lmMst1*(-285 + 32*lmMt -
        2*lmMst2*(55 + 16*lmMt) + 112*pow2(lmMst2)) - 18432*pow3(lmMst2))*pow4(
        Msq) - 20*(155 + 90*lmMst2 + 3*lmMst1*(-23 + 12*lmMst2) + (-69 + 36*
        lmMst2)*lmMt + 12*lmMsq*(4 - 3*(lmMst1 + lmMt)) + 36*pow2(lmMsq) - 36*
        pow2(lmMst2))*pow4(Mst1)) - Mgl*(400*(95 - 18*lmMst2 + 6*lmMst1*(5 + 6*
        lmMst2) + 6*(5 + 6*lmMst2)*lmMt - 6*lmMsq*(7 + 6*(lmMst1 + lmMt)) + 36*
        pow2(lmMsq) - 36*pow2(lmMst2))*pow2(Msq)*pow2(Mst1) + 9*(5253 - 4640*
        lmMst2 + 1200*lmMsq*(13 - 4*lmMst2*(1 + lmMst2) + lmMst1*(2 + 4*lmMst2)
        + 2*lmMt) - 2400*(lmMst1 - lmMst2)*pow2(lmMsq) - 1280*(1 + lmMst2)*
        pow2(lmMst1) + 29600*pow2(lmMst2) - 80*lmMt*(127 + 96*lmMst2 + 32*pow2(
        lmMst2)) - 80*lmMst1*(183 + 274*lmMst2 - 32*(1 + lmMst2)*lmMt + 112*
        pow2(lmMst2)) - 3840*pow2(lmMt) + 10240*pow3(lmMst2))*pow4(Msq) + 100*(
        47 - 54*lmMst2 + (3 + 36*lmMst2)*(lmMst1 + lmMt) + 12*lmMsq*(4 - 3*(
        lmMst1 + lmMt)) + 36*pow2(lmMsq) - 36*pow2(lmMst2))*pow4(Mst1)))*pow6(
        Mst2) - 3*pow2(MuSUSY)*(60*pow2(Msq)*pow2(Mst2)*(-720*(Mgl*(55 + 6*
        lmMsq*(-2 + 5*lmMst1 - 5*lmMst2) + 52*lmMst2 - 10*lmMst1*(4 + 3*lmMst2)
        + 30*pow2(lmMst2)) + Dmglst2*(321 + 18*lmMsq*(-2 + 5*lmMst1 - 5*lmMst2)
        + 96*lmMst2 - 30*lmMst1*(2 + 3*lmMst2) + 90*pow2(lmMst2)))*pow2(Mst1) +
        pow2(Msq)*(243*Dmglst2*(881.6139917695473 + 248*B4 - 4*DN - (2560*
        lmMsq)/3. + lmMst1*(130.96296296296296 - 120*lmMsq - 40*pow2(lmMsq)) +
        (80*pow2(lmMsq))/3. + (176*pow2(lmMst1))/9. + (8*lmMst2*(4364 - 573*
        lmMst1 + 45*lmMsq*(5 + 6*lmMst1) + 135*pow2(lmMsq) + 24*pow2(lmMst1)))/
        27. - (8*(-377 + 90*lmMsq + 376*lmMst1)*pow2(lmMst2))/9. + (2944*pow3(
        lmMst2))/9.) + Mgl*(251578 + 27432*B4 - 108*DN - 58320*lmMsq + 6480*
        pow2(lmMsq) - 216*lmMst1*(422 - 45*lmMsq + 45*pow2(lmMsq)) - 4752*pow2(
        lmMst1) + 216*lmMst2*(1096 + 15*lmMsq*(-7 + 6*lmMst1 - 6*lmMst2) + 717*
        lmMst2 - lmMst1*(551 + 248*lmMst2) + 45*pow2(lmMsq) + 8*pow2(lmMst1)) +
        51840*pow3(lmMst2))))*pow4(Mst1) + 9*pow2(Mst1)*(15*Mgl*(-640*(23 +
        lmMsq*(-6 + 9*lmMst1 - 9*lmMst2) + 18*lmMst2 - 3*lmMst1*(4 + 3*lmMst2)
        + 9*pow2(lmMst2))*pow2(Msq)*pow2(Mst1) + 3*(18201 + 1760*B4 - 16*DN -
        5760*lmMsq + 960*pow2(lmMsq) - 16*lmMst1*(199 - 30*lmMsq + 30*pow2(
        lmMsq)) + 16*lmMst2*(1291 - 322*lmMst1 + 30*lmMsq*(-5 + 2*lmMst1) + 30*
        pow2(lmMsq) - 16*pow2(lmMst1)) - 256*pow2(lmMst1) - 32*(-313 + 30*lmMsq
        + 72*lmMst1)*pow2(lmMst2) + 2560*pow3(lmMst2))*pow4(Msq) - 20*(233 +
        36*lmMsq*(-2 + 5*lmMst1 - 5*lmMst2) + 207*lmMst2 - 45*lmMst1*(3 + 4*
        lmMst2) + 180*pow2(lmMst2))*pow4(Mst1)) + Dmglst2*(-14400*(77 + 6*
        lmMsq*(-2 + 3*lmMst1 - 3*lmMst2) + 24*lmMst2 - 6*lmMst1*(2 + 3*lmMst2)
        + 18*pow2(lmMst2))*pow2(Msq)*pow2(Mst1) + (779917 + 188640*B4 - 3600*DN
        - 648000*lmMsq + 43200*pow2(lmMsq) - 720*lmMst1*(-173 + 90*lmMsq + 30*
        pow2(lmMsq)) + 720*lmMst2*(1623 - 130*lmMst1 + 30*lmMsq*(-1 + 2*lmMst1)
        + 30*pow2(lmMsq) - 16*pow2(lmMst1)) + 11520*pow2(lmMst1) - 1440*(-265 +
        30*lmMsq + 136*lmMst1)*pow2(lmMst2) + 207360*pow3(lmMst2))*pow4(Msq) -
        300*(1961 + 180*lmMsq*(-2 + 5*lmMst1 - 5*lmMst2) + 675*lmMst2 - 45*
        lmMst1*(7 + 20*lmMst2) + 900*pow2(lmMst2))*pow4(Mst1)))*pow4(Mst2) +
        10*Mgl*(2552929 + 257904*B4 - 648*DN - 456840*lmMsq + 38880*pow2(lmMsq)
        - 216*lmMst1*(4591 - 360*lmMsq + 450*pow2(lmMsq)) + 41904*pow2(lmMst1)
        + 216*lmMst2*(9211 - 6466*lmMst1 + 180*lmMsq*(-4 + 5*lmMst1) + 450*
        pow2(lmMsq) + 576*pow2(lmMst1)) - 864*(-1784 + 225*lmMsq + 840*lmMst1)*
        pow2(lmMst2) - 3456*pow3(lmMst1) + 604800*pow3(lmMst2))*pow4(Msq)*pow6(
        Mst1) + 90*(-48*(15*Dmglst2*(95 + 12*lmMsq*(-2 + lmMst1 - lmMst2) + 32*
        lmMst2 - 4*lmMst1*(2 + 3*lmMst2) + 12*pow2(lmMst2)) + 5*Mgl*(71 + 12*
        lmMsq*(-2 + lmMst1 - lmMst2) + 40*lmMst2 - 4*lmMst1*(4 + 3*lmMst2) +
        12*pow2(lmMst2)))*pow2(Msq)*pow2(Mst1) - 5*(Mgl*(911 + 12*lmMsq*(-37 +
        18*lmMst1 - 18*lmMst2) + 606*lmMst2 - 54*lmMst1*(3 + 4*lmMst2) + 216*
        pow2(lmMst2)) + Dmglst2*(5591 + 12*lmMsq*(-181 + 90*lmMst1 - 90*lmMst2)
        + 2550*lmMst2 - 54*lmMst1*(7 + 20*lmMst2) + 1080*pow2(lmMst2)))*pow2(
        Mst1)*pow2(Mst2) + 2304*(1 + lmMst2)*(-Dmglst2 + 3*Dmglst2*lmMst2 + Mgl
        + lmMst2*Mgl)*pow4(Msq) - 5*(Mgl*(1153 + 12*lmMsq*(-35 + 54*lmMst1 -
        54*lmMst2) + 906*lmMst2 - 162*lmMst1*(3 + 4*lmMst2) + 648*pow2(lmMst2))
        + Dmglst2*(8713 + 12*lmMsq*(-179 + 270*lmMst1 - 270*lmMst2) + 3282*
        lmMst2 - 162*lmMst1*(7 + 20*lmMst2) + 3240*pow2(lmMst2)))*pow4(Mst1))*
        pow6(Mst2)) + pow2(Mst1)*pow4(Msq)*(-4*Dmglst2*pow2(Mst2)*(18*pow2(
        Mst1)*(3*(18440*OepS2 - 9*(15691 + 41490*lmMst1 - 41490*lmMst2)*S2)*
        pow2(Mst2) + 4*(5060*OepS2 - 27*(2489 + 3795*lmMst1 - 3795*lmMst2)*S2)*
        pow2(MuSUSY)) + 27*pow2(Mst2)*(50*(56*OepS2 - 81*(-27 + 14*lmMst1 - 14*
        lmMst2)*S2)*pow2(Mst2) + (1400*OepS2 - 81*(-453 + 350*lmMst1 - 350*
        lmMst2)*S2)*pow2(MuSUSY)) + 16*(257300*OepS2 - 27*(220117 + 192975*
        lmMst1 - 192975*lmMst2)*S2)*pow4(Mst1)) + 12*Mgl*(2*(8*(13700*OepS2 -
        27*(16297 + 10275*lmMst1 - 10275*lmMst2)*S2)*pow2(Mst2) + 5*(7528*OepS2
        - 27*(9185 + 5646*lmMst1 - 5646*lmMst2)*S2)*pow2(MuSUSY))*pow4(Mst1) +
        135*(56*OepS2 + 81*(1 - 14*lmMst1 + 14*lmMst2)*S2)*pow2(MuSUSY)*pow4(
        Mst2) + 6*pow2(Mst1)*(10*(544*OepS2 - 81*(169 + 136*lmMst1 - 136*
        lmMst2)*S2)*pow2(Mst2)*pow2(MuSUSY) + (15080*OepS2 - 81*(4367 + 3770*
        lmMst1 - 3770*lmMst2)*S2)*pow4(Mst2)) + 54*(280*OepS2 - 81*(-79 + 70*
        lmMst1 - 70*lmMst2)*S2)*pow6(Mst2))) + 1080*(8*(-5*Mgl*(119 - 63*lmMst2
        + 12*lmMst1*(4 + 3*lmMst2) - 6*(1 + 3*lmMst2)*lmMt + 3*lmMsq*(7 - 12*
        lmMst1 + 18*lmMst2 + 6*lmMt) - 18*pow2(lmMsq) - 36*pow2(lmMst2)) - 15*
        Dmglst2*(125 - 37*lmMst2 + 12*lmMst1*(2 + 3*lmMst2) + 6*lmMt - 6*
        lmMst2*lmMt + lmMsq*(7 - 36*lmMst1 + 42*lmMst2 + 6*lmMt) - 6*pow2(
        lmMsq) - 36*pow2(lmMst2)))*pow2(Msq)*pow2(Mst1) - 15*(Dmglst2*(754 -
        50*lmMst2 + 9*lmMst1*(7 + 20*lmMst2) - 4*lmMsq*(10 + 45*lmMst1 - 66*
        lmMst2 - 21*lmMt) + 27*lmMt - 84*lmMst2*lmMt - 84*pow2(lmMsq) - 180*
        pow2(lmMst2)) + Mgl*(206 + 6*lmMst2 + 9*lmMst1*(3 + 4*lmMst2) - 12*
        lmMsq*(2 + 3*lmMst1 - 6*lmMst2 - 3*lmMt) - 9*(1 + 4*lmMst2)*lmMt - 36*
        pow2(lmMsq) - 36*pow2(lmMst2)))*pow2(Mst1)*pow2(Mst2) - 1152*(1 +
        lmMst2)*(-Dmglst2 + 3*Dmglst2*lmMst2 + Mgl + lmMst2*Mgl)*pow4(Msq) - 5*
        (2*Mgl*(5 - 54*lmMst2 + 3*lmMst1*(7 + 12*lmMst2) + 3*(7 + 12*lmMst2)*
        lmMt + 12*lmMsq*(1 - 3*(lmMst1 + lmMt)) + 36*pow2(lmMsq) - 36*pow2(
        lmMst2)) + 6*Dmglst2*(41 - 6*lmMst2 + 3*lmMst1*(-1 + 12*lmMst2) + (-3 +
        36*lmMst2)*lmMt + 12*lmMsq*(1 - 3*(lmMst1 + lmMt)) + 36*pow2(lmMsq) -
        36*pow2(lmMst2)))*pow4(Mst1))*pow8(Mst2)))/(43740.*Mgl*pow2(Mst1)*pow4(
        Msq)*pow5(Mst2)) + ((1 - 2*lmMst2)*shiftst3*((8*(-1 + lmMst1 - lmMst2)*
        pow2(Mst1)*pow2(Mt)*pow2(s2t) - 16*pow4(Mt) - (-4*pow2(Mst1)*pow2(Mst2)
        + (3 + 2*lmMst1 - 2*lmMst2)*pow4(Mst1) + pow4(Mst2))*pow4(s2t))*pow6(
        Mst2) - (4*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*(2*(1 - 2*lmMst1 + 2*lmMst2)
        *pow2(Mst2)*pow4(Mst1) + 2*(1 - lmMst1 + lmMst2)*pow2(Mst1)*pow4(Mst2)
        + (2 - 6*lmMst1 + 6*lmMst2)*pow6(Mst1) + pow6(Mst2)))/pow2(Sbeta) + 4*
        pow2(Mt)*pow2(s2t)*(pow2(MuSUSY)*pow6(Mst2) + 2*(pow2(MuSUSY)*((1 - 2*
        lmMst1 + 2*lmMst2)*pow2(Mst2)*pow4(Mst1) + (1 - lmMst1 + lmMst2)*pow2(
        Mst1)*pow4(Mst2) + (1 - 3*lmMst1 + 3*lmMst2)*pow6(Mst1)) + pow8(Mst2)))
        ))/(12.*pow2(Mst1)*pow4(Mst2)))) + (Al4p*xDmglst2*pow2(Dmglst2)*((
        twoLoopFlag*((36*Mt*pow3(s2t)*(4*(3 - lmMst1*(-2 + lmMst2) - 4*lmMst2 +
        pow2(lmMst2))*pow2(Mst1)*pow2(Mst2) + (7 - 6*lmMst1 + 6*lmMst2)*pow4(
        Mst1) + 4*(14 - lmMst1*(-2 + lmMst2) + pow2(lmMst2))*pow4(Mst2)))/Mst2
        + (32*pow4(Mt)*(-9*(-27 - 19*lmMst2 + 2*lmMst1*(6 + lmMst2 - lmMt) + 7*
        lmMt + 2*lmMst2*lmMt - 2*pow2(lmMst2))*pow2(Mst2)*pow4(Mst1) + (-2 + 9*
        lmMst1 + 3*lmMst2 + 15*lmMt)*pow2(Mst1)*pow4(Mst2) + 27*(23 + 38*lmMst2
        - 9*lmMt - 6*lmMst2*lmMt + lmMst1*(-29 - 6*lmMst2 + 6*lmMt) + 6*pow2(
        lmMst2))*pow6(Mst1) + 9*(-2 + 3*lmMst2)*pow6(Mst2)))/(pow2(Mst1)*pow4(
        Mst2)) + (8*s2t*pow3(Mt)*(9*(4*(-28 - 41*lmMst2 + lmMst1*(31 + 6*lmMst2
        - 6*lmMt) + 10*lmMt + 6*lmMst2*lmMt - 6*pow2(lmMst2))*pow2(Mst2) + (-
        199 + 26*lmMst2 + lmMst1*(-34 + 20*lmMst2) - 20*pow2(lmMst2))*pow2(
        MuSUSY))*pow4(Mst1) + 36*(-14 + lmMst1*(-2 + lmMst2) - pow2(lmMst2))*
        pow2(MuSUSY)*pow4(Mst2) + 36*pow2(Mst1)*((-31 + 3*lmMst1*(-2 + lmMst2)
        + 4*lmMst2 - 3*pow2(lmMst2))*pow2(Mst2)*pow2(MuSUSY) + (-13 + 2*lmMst1
        - 4*lmMst2 + 2*lmMt)*pow4(Mst2)) + 2*(161 + 36*lmMst1*(-2 + lmMst2) +
        78*lmMst2 - 6*lmMt - 36*pow2(lmMst2))*pow6(Mst2)))/pow5(Mst2) + (12*
        pow2(Mt)*pow2(s2t)*(2*pow4(Mst1)*(9*(4 - lmMst1 + lmMst2 + 2*lmMst1*
        lmMst2 - 2*pow2(lmMst2))*pow2(Mst2)*pow2(MuSUSY) + (-7 + 6*lmMst1 - 24*
        lmMst2)*pow4(Mst2)) + (12*(-4 + 5*lmMst1 - 5*lmMst2)*pow2(Mst2) + 9*(9
        + 4*lmMst1*(-1 + lmMst2) + 4*lmMst2 - 4*pow2(lmMst2))*pow2(MuSUSY))*
        pow6(Mst1) - 6*(-2 + 3*lmMst2)*(2*pow2(Mst2) + pow2(MuSUSY))*pow6(Mst2)
        + 2*pow2(Mst1)*(3*(10 - 3*lmMst1 + 6*lmMst1*lmMst2 - 6*pow2(lmMst2))*
        pow2(MuSUSY)*pow4(Mst2) + (-23 + 36*lmMst2)*pow6(Mst2))))/(pow2(Mst1)*
        pow4(Mst2)) + (9*pow4(s2t)*(-6*(-2 + lmMst1 - 2*lmMst1*lmMst2 + 2*pow2(
        lmMst2))*pow2(Mst2)*pow4(Mst1) + 6*(-2 + lmMst1 - 2*lmMst2 - 2*lmMst1*
        lmMst2 + 2*pow2(lmMst2))*pow2(Mst1)*pow4(Mst2) + (1 + 6*lmMst1)*pow6(
        Mst1) + 2*(-2 + 3*lmMst2)*pow6(Mst2)) + (4*Mt*MuSUSY*((9*Mt*MuSUSY*s2t*
        (-2*(4*Mt*(-31 + 3*lmMst1*(-2 + lmMst2) + 4*lmMst2 - 3*pow2(lmMst2)) +
        3*Mst2*s2t*(4 + lmMst2 + lmMst1*(-1 + 2*lmMst2) - 2*pow2(lmMst2)))*
        pow2(Mst2)*pow4(Mst1) - 2*(Mst2*s2t*(10 + lmMst1*(-3 + 6*lmMst2) - 6*
        pow2(lmMst2)) + 4*Mt*(-14 + lmMst1*(-2 + lmMst2) - pow2(lmMst2)))*pow2(
        Mst1)*pow4(Mst2) + (-3*Mst2*s2t*(9 + 4*lmMst1*(-1 + lmMst2) + 4*lmMst2
        - 4*pow2(lmMst2)) + Mt*(398 + lmMst1*(68 - 40*lmMst2) - 52*lmMst2 + 40*
        pow2(lmMst2)))*pow6(Mst1) + 2*(-2 + 3*lmMst2)*s2t*pow7(Mst2)))/pow2(
        Sbeta) - (-2*pow2(Mst2)*(36*(3 - lmMst1 + lmMst2)*Mst2*s2t*pow2(Mt) +
        54*Mt*(-17 + 2*lmMst1*(-2 + lmMst2) + 4*lmMst2 - 2*pow2(lmMst2))*pow2(
        Mst2)*pow2(s2t) + 4*(29 - 18*lmMst1*(-1 + lmMst2) - 3*lmMst2 - 15*lmMt
        + 18*pow2(lmMst2))*pow3(Mt) + 9*(2 + 3*lmMst2)*pow3(Mst2)*pow3(s2t))*
        pow4(Mst1) + 2*pow2(Mst1)*(6*(-11 + 18*lmMst2)*Mst2*s2t*pow2(Mt) + 54*
        Mt*(14 - lmMst1*(-2 + lmMst2) + pow2(lmMst2))*pow2(Mst2)*pow2(s2t) + 4*
        (79 + 18*lmMst1*(-2 + lmMst2) + 39*lmMst2 - 3*lmMt - 18*pow2(lmMst2))*
        pow3(Mt) - 9*(8 + 3*lmMst2 + lmMst1*(-3 + 6*lmMst2) - 6*pow2(lmMst2))*
        pow3(Mst2)*pow3(s2t))*pow4(Mst2) - (72*(7 - 6*lmMst1 + 6*lmMst2)*Mst2*
        s2t*pow2(Mt) - 27*Mt*(75 + lmMst1*(10 - 8*lmMst2) - 10*lmMst2 + 8*pow2(
        lmMst2))*pow2(Mst2)*pow2(s2t) + 8*(236 + 339*lmMst2 - 18*lmMst1*(13 +
        4*lmMst2 - 3*lmMt) - 3*(35 + 18*lmMst2)*lmMt + 72*pow2(lmMst2))*pow3(
        Mt) + 27*(1 - 2*lmMst1 + 2*lmMst2)*pow3(Mst2)*pow3(s2t))*pow6(Mst1) +
        18*(-2 + 3*lmMst2)*s2t*(-4*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow7(Mst2))
        /Tbeta))/pow5(Mst2))/pow2(Mst1)))/108. + Al4p*threeLoopFlag*((
        3676.7919341563784 - (976*lmMsq)/135. - (20*(103 + 6*lmMsq)*lmMst1)/9.
         - (4*(-350909 + 39600*lmMsq + 450*lmMst1)*lmMst2)/
        2025. + (4*(6293 - 960*(lmMsq + lmMst1) + 3704*lmMst2)*lmMt)/135. + 60*
        pow2(lmMsq) - (28*pow2(lmMst1))/9. + (4*(-629 + 2880*lmMst1)*pow2(
        lmMst2))/135. + (464*pow2(lmMt))/9. - (40*(9 + 2*lmMst1 - 4*lmMst2 + 2*
        lmMt)*pow2(Mst1))/(9.*pow2(Msq)) + (2*(-1021 - 502*lmMst2 + 64*lmMst1*(
        -2 + 3*lmMst2) - 90*lmMsq*(-5 + 6*lmMst2) + 270*pow2(lmMsq) + 624*pow2(
        lmMst2))*pow2(Mst2))/(9.*pow2(Mst1)) + (pow2(Mst1)*(3582403037 +
        2058687540*lmMst2 - 37800*(1783 + 285*lmMst2 - 225*lmMt)*pow2(lmMst1) +
        592080300*pow2(lmMst2) + 8505000*lmMsq*(65 + lmMst2*(27 - 2*lmMt) - 2*
        lmMst1*(8 + lmMst2 - lmMt) - 11*lmMt + 2*pow2(lmMst2)) - 56700*lmMt*(-
        5825 + 1311*lmMst2 + 10*pow2(lmMst2)) - 13608000*(7 + 2*lmMst2)*pow2(
        lmMt) + 3780*lmMst1*(-371868 + 19265*lmMt - 15*lmMst2*(7547 + 140*lmMt)
        + 9900*pow2(lmMst2) + 7200*pow2(lmMt)) + 756000*pow3(lmMst1) -
        27405000*pow3(lmMst2)))/(637875.*pow2(Mst2)) - (256*pow3(lmMst2))/3. -
        pow4(Mst1)*((20*(lmMst1 - 2*lmMst2 + lmMt))/(3.*pow4(Msq)) + (
        3879.568149934479 + 1840*lmMsq + (581.2930880994373 - 2320*lmMsq)*
        lmMst1 - (126272*pow2(lmMst1))/189. - (4*lmMst2*(57911521 - 22348620*
        lmMst1 + 2381400*lmMsq*(-19 + 3*lmMst1) + 2487240*pow2(lmMst1)))/59535.
         + (16*(-15893 + 5670*lmMsq + 4284*lmMst1)*pow2(lmMst2))/189. + (16*
        lmMt*(2105 + 405*lmMsq*(-3 + 2*lmMst1 - 2*lmMst2) + 2247*lmMst2 - 3*
        lmMst1*(341 + 364*lmMst2) + 234*pow2(lmMst1) + 858*pow2(lmMst2)))/27. +
        (128*(-5 + 6*lmMst1 - 6*lmMst2)*pow2(lmMt))/3. + (256*pow3(lmMst1))/27.
         - (5536*pow3(lmMst2))/27.)/pow4(Mst2)) - (128*(-2 + lmMst2 + 5*pow2(
        lmMst2))*pow4(Mst2))/(9.*pow4(Mst1)) - (4*(150*(87 - 6*lmMst2*(-10 +
        lmMt) - 20*lmMt + lmMst1*(-40 - 6*lmMst2 + 6*lmMt) + 6*pow2(lmMst2))*
        pow4(Mst1) + (2048 + 3330*lmMst2 - 675*lmMst1*(-1 + 2*lmMst2) + 1695*
        lmMt - 2700*lmMst2*lmMt + 150*lmMsq*(-38 + 9*lmMst1 + 9*lmMst2 + 18*
        lmMt) - 2700*pow2(lmMsq) + 1350*pow2(lmMst2))*pow4(Mst2)))/(135.*pow2(
        Msq)*pow2(Mst2)) + (pow2(Mst2)*((-90263 - 176280*lmMst2 + 750*lmMst1*(-
        13 + 60*lmMst2) - 42870*lmMt + 106200*lmMst2*lmMt - 300*lmMsq*(-763 +
        150*lmMst1 + 204*lmMst2 + 354*lmMt) + 106200*pow2(lmMsq) - 45000*pow2(
        lmMst2))*pow2(Mst1)*pow2(Mst2) + 1350*(29 + 14*lmMst2 + 3*lmMst1*(-5 +
        4*lmMst2) - 15*lmMt + 12*lmMst2*lmMt - 4*lmMsq*(-4 + 3*lmMst1 + 3*lmMt)
        + 12*pow2(lmMsq) - 12*pow2(lmMst2))*pow4(Mst1) + 75*(160*(-8 + 15*lmMsq
        - 15*lmMst2)*pow2(Msq)*pow2(Mst2) + 63*(-5 + 28*lmMsq - 28*lmMst2)*
        pow4(Mst2))))/(810.*pow2(Mst1)*pow4(Msq)) + (2*(-189*pow2(MuSUSY)*(4*
        pow2(Mst2)*(50134 - 270*B4 + 270*D3 - 135*DN + 2160*lmMst1 + 120*(271 +
        24*lmMst1)*lmMst2 - 120*lmMt - 720*(-2 + 3*lmMst1)*pow2(lmMst2) + 2160*
        pow3(lmMst2)) + 5*pow2(Mst1)*(28405 - 288*B4 + 288*D3 - 144*DN + 5472*
        lmMst1 - 288*(-89 + 2*lmMst1)*lmMst2 + 576*(-2 + lmMst1 - lmMst2)*lmMt
        - 576*(-9 + 8*lmMst1)*pow2(lmMst2) + 4608*pow3(lmMst2))) + 560*OepS2*(-
        3*pow2(Mst1)*(4703*pow2(Mst2) + 126*pow2(MuSUSY)) + 11164*pow4(Mst1) -
        819*pow4(Mst2)) - 27*S2*(3*pow2(Mst1)*((128797 - 1975260*lmMst1 +
        1975260*lmMst2)*pow2(Mst2) - 3780*(65 + 14*lmMst1 - 14*lmMst2)*pow2(
        MuSUSY)) + (4002484 + 4688880*lmMst1 - 4688880*lmMst2)*pow4(Mst1) +
        189*(-4392*pow2(Mst2)*pow2(MuSUSY) + (3869 - 1820*lmMst1 + 1820*lmMst2)
        *pow4(Mst2)))))/(76545.*pow4(Mst2)))*pow4(Mt) + (s2t*pow3(Mt)*((-3*
        pow2(Mst1)*(90713501 + 5443200*lmMsq + 61803000*lmMst2 - 15120*(-652 +
        7*lmMst2)*lmMt + 725760*(-5 + 2*lmMst2)*pow2(lmMst1) + 4460400*pow2(
        lmMst2) - 7560*lmMst1*(4679 - 370*lmMst2 + 178*lmMt + 384*pow2(lmMst2))
        - 2177280*pow2(lmMt) + (453600*(-15 + 2*lmMst1 - 4*lmMst2 + 2*lmMt)*
        pow2(Mst1))/pow2(Msq) + 1451520*pow3(lmMst2)))/Mst2 - 90720*((5*(259 -
        31*lmMst2 + 12*lmMst1*(1 + 6*lmMst2) - 6*(-2 + lmMst2)*lmMt + lmMsq*(7
        - 72*lmMst1 + 78*lmMst2 + 6*lmMt) - 6*pow2(lmMsq) - 72*pow2(lmMst2)))/
        pow2(Msq) - (96*(2 + lmMst2 - 3*pow2(lmMst2)))/pow2(Mst1))*pow3(Mst2) +
        (28*(24208447 + 1263240*lmMst2 - 9720*(266 + 74*lmMst2 - 55*lmMt)*pow2(
        lmMst1) - 1140480*pow2(lmMst2) + 145800*lmMsq*(28 + lmMst2*(41 - 6*
        lmMt) - 10*lmMt + lmMst1*(-31 - 6*lmMst2 + 6*lmMt) + 6*pow2(lmMst2)) +
        3240*lmMt*(846 + 1042*lmMst2 + 309*pow2(lmMst2)) - 233280*(2 + 3*
        lmMst2)*pow2(lmMt) + 360*lmMst1*(-21491 + lmMst2*(8190 - 4266*lmMt) -
        5922*lmMt + 3267*pow2(lmMst2) + 1944*pow2(lmMt)) + 9720*pow3(lmMst1) -
        466560*pow3(lmMst2))*pow4(Mst1))/pow3(Mst2) - (9*Mst2*(8400*(5 + 198*
        lmMst2 + 6*lmMst1*(-13 + 6*lmMst2) + 6*(-13 + 6*lmMst2)*lmMt - 6*lmMsq*
        (7 + 6*(lmMst1 + lmMt)) + 36*pow2(lmMsq) - 36*pow2(lmMst2))*pow2(Msq)*
        pow2(Mst1) + (36729361 - 13527360*lmMst2 - 453600*(lmMst1 - lmMst2)*
        pow2(lmMsq) - 241920*(-2 + lmMst2)*pow2(lmMst1) + 15120*lmMst1*(135 -
        64*lmMt + lmMst2*(302 + 32*lmMt) - 112*pow2(lmMst2)) + 25200*lmMsq*(131
        + 96*lmMst2 + 18*lmMst1*(-5 + 2*lmMst2) - 6*lmMt - 36*pow2(lmMst2)) -
        8316000*pow2(lmMst2) - 5040*lmMt*(-7 - 240*lmMst2 + 96*pow2(lmMst2)) +
        120960*pow2(lmMt) + 1935360*pow3(lmMst2))*pow4(Msq) + 2100*(155 + 90*
        lmMst2 + 3*lmMst1*(-23 + 12*lmMst2) + (-69 + 36*lmMst2)*lmMt + 12*
        lmMsq*(4 - 3*(lmMst1 + lmMt)) + 36*pow2(lmMsq) - 36*pow2(lmMst2))*pow4(
        Mst1)))/pow4(Msq) - 153090*pow2(MuSUSY)*((-16*Mst2*((5*(397 + 36*lmMsq*
        (-2 + lmMst1 - lmMst2) + 78*lmMst2 - 6*lmMst1*(1 + 6*lmMst2) + 36*pow2(
        lmMst2)))/pow2(Msq) + (48*(2 + lmMst2 - 3*pow2(lmMst2)))/pow2(Mst1)))/
        27. + (pow2(Mst1)*(1409.4528120713305 + (1352*B4)/3. - (28*DN)/3. - (
        16520*lmMsq)/9. + (80*pow2(lmMsq))/3. - (40*lmMst1*(167 + 405*lmMsq +
        81*pow2(lmMsq)))/81. - (128*pow2(lmMst1))/9. + (8*lmMst2*(30469 + 2709*
        lmMst1 + 135*lmMsq*(11 + 6*lmMst1) + 405*pow2(lmMsq) + 72*pow2(lmMst1))
        )/81. - (8*(-19 + 90*lmMsq + 568*lmMst1)*pow2(lmMst2))/9. - (80*(107 +
        6*lmMsq*(-2 + 5*lmMst1 - 5*lmMst2) + 32*lmMst2 - 10*lmMst1*(2 + 3*
        lmMst2) + 30*pow2(lmMst2))*pow2(Mst1))/(9.*pow2(Msq)) + (4480*pow3(
        lmMst2))/9.))/pow3(Mst2) - (5*(Mst2*(37213 + 12*lmMsq*(-539 + 810*
        lmMst1 - 810*lmMst2) + 6630*lmMst2 - 162*lmMst1*(1 + 60*lmMst2) + 9720*
        pow2(lmMst2))*pow2(Mst1) + (21275 + 12*lmMsq*(-541 + 270*lmMst1 - 270*
        lmMst2) + 6546*lmMst2 - 54*lmMst1*(1 + 60*lmMst2) + 3240*pow2(lmMst2))*
        pow3(Mst2)))/(162.*pow4(Msq)) + (89.18662551440329 + (1960*B4)/9. - (
        44*DN)/9. - (7880*lmMsq)/9. + (80*pow2(lmMsq))/3. - (4*lmMst1*(-359 +
        150*lmMsq + 30*pow2(lmMsq)))/9. + (4*lmMst2*(8503 + 666*lmMst1 + 90*
        lmMsq*(1 + 2*lmMst1) + 90*pow2(lmMsq) - 48*pow2(lmMst1)))/27. + (128*
        pow2(lmMst1))/9. - (8*(-177 + 30*lmMsq + 232*lmMst1)*pow2(lmMst2))/9. -
        (40*(1445 + 72*lmMsq*(-2 + 3*lmMst1 - 3*lmMst2) + 180*lmMst2 - 36*
        lmMst1*(1 + 6*lmMst2) + 216*pow2(lmMst2))*pow2(Mst1))/(27.*pow2(Msq)) +
        (640*pow3(lmMst2))/3. - (5*(1961 + 180*lmMsq*(-2 + 5*lmMst1 - 5*lmMst2)
        + 675*lmMst2 - 45*lmMst1*(7 + 20*lmMst2) + 900*pow2(lmMst2))*pow4(Mst1)
        )/(27.*pow4(Msq)))/Mst2) + (-560*OepS2*(-63*pow2(Mst2)*(2*pow2(Mst2) +
        pow2(MuSUSY)) - 6*pow2(Mst1)*(4307*pow2(Mst2) + 1555*pow2(MuSUSY)) +
        102920*pow4(Mst1)) + 54*S2*(6*pow2(Mst1)*((760703 - 904470*lmMst1 +
        904470*lmMst2)*pow2(Mst2) + 7*(20803 - 46650*lmMst1 + 46650*lmMst2)*
        pow2(MuSUSY)) + 112*(220117 + 192975*lmMst1 - 192975*lmMst2)*pow4(Mst1)
        + 27*(7*(1387 - 70*lmMst1 + 70*lmMst2)*pow2(Mst2)*pow2(MuSUSY) + (18722
        - 980*lmMst1 + 980*lmMst2)*pow4(Mst2))))/pow3(Mst2) - (56700*(4*(44 -
        21*lmMst1 + 6*(5 + 6*lmMst1)*lmMst2 + 3*(-7 + 12*lmMst2)*lmMt + 12*
        lmMsq*(1 - 3*(lmMst1 + lmMt)) + 36*pow2(lmMsq) - 36*pow2(lmMst2))*pow2(
        Mst1)*pow3(Mst2) + (2100 - 62*lmMst2 + 9*lmMst1*(1 + 60*lmMst2) + 117*
        lmMt - 156*lmMst2*lmMt + lmMsq*(-64 - 540*lmMst1 + 696*lmMst2 + 156*
        lmMt) - 156*pow2(lmMsq) - 540*pow2(lmMst2))*pow5(Mst2)))/pow4(Msq)))/
        153090. + ((T1ep*(-2*pow4(Mst1)*(8*pow2(Mt)*pow2(s2t)*(4327*pow2(
        MuSUSY)*(-1 + pow2(Sbeta)) + 6921*pow2(Mst2)*pow2(Sbeta)) - 205840*
        Mst2*s2t*pow2(Sbeta)*pow3(Mt) + 15571*Mt*pow2(Sbeta)*pow3(Mst2)*pow3(
        s2t) + 89312*pow2(Sbeta)*pow4(Mt) - 4199*pow2(Sbeta)*pow4(Mst2)*pow4(
        s2t)) + 63*pow3(Mst2)*(-4*Mst2*pow2(Mt)*pow2(s2t)*(5*pow2(MuSUSY)*(-1 +
        pow2(Sbeta)) + 18*pow2(Mst2)*pow2(Sbeta)) - 4*s2t*(pow2(MuSUSY)*(-1 +
        pow2(Sbeta)) + 2*pow2(Mst2)*pow2(Sbeta))*pow3(Mt) + 2*Mt*pow2(Sbeta)*
        pow3(s2t)*pow4(Mst2) + 208*Mst2*pow2(Sbeta)*pow4(Mt) + 5*pow2(Sbeta)*
        pow4(s2t)*pow5(Mst2)) + 24*pow2(Mst1)*(pow2(Mst2)*pow2(Mt)*pow2(s2t)*(-
        185*pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 1791*pow2(Mst2)*pow2(Sbeta)) +
        Mst2*s2t*(-1555*pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 4307*pow2(Mst2)*pow2(
        Sbeta))*pow3(Mt) + (252*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 9406*pow2(
        Mst2)*pow2(Sbeta))*pow4(Mt) + 767*Mt*pow2(Sbeta)*pow3(s2t)*pow5(Mst2) +
        20*pow2(Sbeta)*pow4(s2t)*pow6(Mst2))))/(729.*pow4(Mst2)) + pow2(Mt)*
        pow2(MuSUSY)*((2*pow2(Mt)*(4*pow2(Mst2)*(50134 - 270*B4 + 270*D3 - 135*
        DN + 2160*lmMst1 + 120*(271 + 24*lmMst1)*lmMst2 - 120*lmMt - 29646*S2 -
        720*(-2 + 3*lmMst1)*pow2(lmMst2) + 2160*pow3(lmMst2)) + 5*pow2(Mst1)*(
        28405 - 288*B4 + 288*D3 - 144*DN + 5472*lmMst1 - 288*(-89 + 2*lmMst1)*
        lmMst2 + 576*(-2 + lmMst1 - lmMst2)*lmMt + 224*OepS2 - 324*(65 + 14*
        lmMst1 - 14*lmMst2)*S2 - 576*(-9 + 8*lmMst1)*pow2(lmMst2) + 4608*pow3(
        lmMst2))))/(405.*pow4(Mst2)) + pow2(s2t)*(112.1664609053498 - (100*B4)/
        9. + (112*D3)/9. - (62*DN)/9. + (1055*lmMsq)/9. + 15*pow2(lmMsq) - (2*
        lmMst1*(-544 + 210*lmMsq + 135*pow2(lmMsq)))/9. - (103*pow2(lmMst1))/9.
         + (lmMst2*(-7921 - 1476*lmMst1 + 90*lmMsq*(5 + 18*lmMst1) + 810*pow2(
        lmMsq) + 288*pow2(lmMst1)))/27. + (93.66666666666667 - 60*lmMsq - 112*
        lmMst1)*pow2(lmMst2) + (10*(405 + (-187 + 246*lmMsq)*lmMst2 + lmMst1*(
        187 - 246*lmMsq + 246*lmMst2) - 246*pow2(lmMst2))*pow2(Mst1))/(27.*
        pow2(Msq)) + ((-1149 + 266*lmMst2 + 64*lmMst1*(-2 + 3*lmMst2) - 90*
        lmMsq*(-5 + 6*lmMst2) + 270*pow2(lmMsq) + 1264*pow2(lmMst2))*pow2(Mst2)
        )/(18.*pow2(Mst1)) - (pow2(Mst1)*(73.78472976680384 + (164*B4)/9. - (
        176*D3)/9. + (94*DN)/9. - (260*lmMsq)/3. + lmMst1*(-96.55703703703703 +
        40*lmMsq + 30*pow2(lmMsq)) + (7556*pow2(lmMst1))/135. - (2*lmMst2*(-
        99788 + 2005*lmMst1 + 6750*lmMsq*(2 + 3*lmMst1) + 10125*pow2(lmMsq) +
        32775*pow2(lmMst1)))/675. + (2*(-3377 + 4050*lmMsq + 19155*lmMst1)*
        pow2(lmMst2))/135. + (10*pow3(lmMst1))/27. - (5050*pow3(lmMst2))/27.))/
        pow2(Mst2) + (304*pow3(lmMst2))/3. + pow4(Mst1)*((5*(216 + (-95 + 132*
        lmMsq)*lmMst2 + lmMst1*(95 - 132*lmMsq + 132*lmMst2) - 132*pow2(lmMst2)
        ))/(54.*pow4(Msq)) + (536.1152102791342 - (8*B4)/3. + (32*D3)/9. - (20*
        DN)/9. + 90*lmMsq - (123.11224321827497 + 20*lmMsq*(1 + lmMsq))*lmMst1
        - lmMst2*(17.33220122616948 - 20*lmMsq*(1 + lmMsq) + (
        133.04550264550264 - 40*lmMsq)*lmMst1 - (1180*pow2(lmMst1))/9.) - (
        15886*pow2(lmMst1))/945. + (149.85608465608465 - 40*lmMsq - (2812*
        lmMst1)/9.)*pow2(lmMst2) - (92*pow3(lmMst1))/27. + (4988*pow3(lmMst2))/
        27.)/pow4(Mst2)) - (32*(-2 + lmMst2 + 5*pow2(lmMst2))*pow4(Mst2))/(9.*
        pow4(Mst1)) - (6*(7400*OepS2 + 27*(1089707 - 5550*lmMst1 + 5550*lmMst2)
        *S2)*pow2(Mst1)*pow2(Mst2) + 40*(17308*OepS2 + 27*(93919 - 12981*lmMst1
        + 12981*lmMst2)*S2)*pow4(Mst1) + 9*(1400*OepS2 + 81*(116129 - 350*
        lmMst1 + 350*lmMst2)*S2)*pow4(Mst2))/(10935.*pow4(Mst2)) + (5*((291 +
        2*(-103 + 84*lmMsq)*lmMst2 + 2*lmMst1*(103 - 84*lmMsq + 84*lmMst2) -
        168*pow2(lmMst2))*pow4(Mst1) + (684 + 347*lmMst1 + 7*(-77 + 78*lmMst1)*
        lmMst2 + 6*lmMsq*(32 - 91*lmMst1 + 91*lmMst2) - 546*pow2(lmMst2))*pow4(
        Mst2)))/(27.*pow2(Msq)*pow2(Mst2)) + (5*((3834 + (-890 + 2328*lmMsq)*
        lmMst2 + lmMst1*(890 - 2328*lmMsq + 2328*lmMst2) - 2328*pow2(lmMst2))*
        pow2(Mst2)*pow4(Mst1) + 160*(-8 + 15*lmMsq - 15*lmMst2)*pow2(Msq)*pow4(
        Mst2) + (4239 - 1324*lmMst2 + 8*lmMst1*(95 + 366*lmMst2) + lmMsq*(564 -
        2928*lmMst1 + 2928*lmMst2) - 2928*pow2(lmMst2))*pow2(Mst1)*pow4(Mst2) +
        63*(-5 + 28*lmMsq - 28*lmMst2)*pow6(Mst2)))/(216.*pow2(Mst1)*pow4(Msq))
        ) - (Mt*s2t*(pow4(Msq)*(-3*pow2(Mst1)*pow2(Mst2)*(433447 + 1058400*B4 -
        23760*DN - 4255200*lmMsq - 1120*OepS2 + 324*(-1387 + 70*lmMst1 - 70*
        lmMst2)*S2 + 129600*pow2(lmMsq) - 2160*lmMst1*(-359 + 150*lmMsq + 30*
        pow2(lmMsq)) + 720*lmMst2*(8503 + 666*lmMst1 + 90*lmMsq*(1 + 2*lmMst1)
        + 90*pow2(lmMsq) - 48*pow2(lmMst1)) + 69120*pow2(lmMst1) - 4320*(-177 +
        30*lmMsq + 232*lmMst1)*pow2(lmMst2) + 1036800*pow3(lmMst2)) - 2*(
        10274911 + 3285360*B4 - 68040*DN - 13381200*lmMsq - 248800*OepS2 + 108*
        (-20803 + 46650*lmMst1 - 46650*lmMst2)*S2 + 194400*pow2(lmMsq) - 3600*
        lmMst1*(167 + 405*lmMsq + 81*pow2(lmMsq)) - 103680*pow2(lmMst1) + 720*
        lmMst2*(30469 + 2709*lmMst1 + 135*lmMsq*(11 + 6*lmMst1) + 405*pow2(
        lmMsq) + 72*pow2(lmMst1)) - 6480*(-19 + 90*lmMsq + 568*lmMst1)*pow2(
        lmMst2) + 3628800*pow3(lmMst2))*pow4(Mst1) + 414720*(2 + lmMst2 - 3*
        pow2(lmMst2))*pow4(Mst2)) + 21600*pow2(Msq)*((1445 + 72*lmMsq*(-2 + 3*
        lmMst1 - 3*lmMst2) + 180*lmMst2 - 36*lmMst1*(1 + 6*lmMst2) + 216*pow2(
        lmMst2))*pow2(Mst2)*pow4(Mst1) + 2*(397 + 36*lmMsq*(-2 + lmMst1 -
        lmMst2) + 78*lmMst2 - 6*lmMst1*(1 + 6*lmMst2) + 36*pow2(lmMst2))*pow2(
        Mst1)*pow4(Mst2) + 6*(107 + 6*lmMsq*(-2 + 5*lmMst1 - 5*lmMst2) + 32*
        lmMst2 - 10*lmMst1*(2 + 3*lmMst2) + 30*pow2(lmMst2))*pow6(Mst1)) + 450*
        ((37213 + 12*lmMsq*(-539 + 810*lmMst1 - 810*lmMst2) + 6630*lmMst2 -
        162*lmMst1*(1 + 60*lmMst2) + 9720*pow2(lmMst2))*pow4(Mst1)*pow4(Mst2) +
        6*(1961 + 180*lmMsq*(-2 + 5*lmMst1 - 5*lmMst2) + 675*lmMst2 - 45*
        lmMst1*(7 + 20*lmMst2) + 900*pow2(lmMst2))*pow2(Mst2)*pow6(Mst1) + (
        21275 + 12*lmMsq*(-541 + 270*lmMst1 - 270*lmMst2) + 6546*lmMst2 - 54*
        lmMst1*(1 + 60*lmMst2) + 3240*pow2(lmMst2))*pow2(Mst1)*pow6(Mst2))))/(
        14580.*pow2(Mst1)*pow3(Mst2)*pow4(Msq))))/pow2(Sbeta) - Mt*pow3(s2t)*(-
        (Mst2*pow2(Mst1)*(587.0953360768176 + (68*B4)/9. + (2*DN)/9. - (380*
        lmMsq)/9. - (40*pow2(lmMsq))/3. - (4*lmMst1*(4066 + 675*lmMsq + 135*
        pow2(lmMsq)))/81. - (64*pow2(lmMst1))/3. + (4*lmMst2*(4672 + 711*lmMst1
        + 135*lmMsq*(9 + 2*lmMst1) + 135*pow2(lmMsq) + 216*pow2(lmMst1)))/81. -
        (4*(239 + 30*lmMsq + 104*lmMst1)*pow2(lmMst2))/9. + (320*pow3(lmMst2))/
        9.)) + (20*Mst2*pow4(Mst1))/(9.*pow2(Msq)) + ((406.2899291266575 - (
        24322*lmMst2)/243. + lmMsq*(23.333333333333332 - 20*lmMst1 + 20*lmMst2)
        + (4*(37 - 40*lmMst2)*pow2(lmMst1))/9. + (148*pow2(lmMst2))/9. + (2*
        lmMst1*(8705 - 3996*lmMst2 + 5616*pow2(lmMst2)))/243. - (32*pow3(
        lmMst1))/9. - (224*pow3(lmMst2))/9.)*pow4(Mst1))/Mst2 - pow3(Mst2)*(
        101.48220164609053 + (980*B4)/9. - (22*DN)/9. - (3940*lmMsq)/9. + (40*
        pow2(lmMsq))/3. - (2*lmMst1*(-359 + 150*lmMsq + 30*pow2(lmMsq)))/9. + (
        2*lmMst2*(8887 + 666*lmMst1 + 90*lmMsq*(1 + 2*lmMst1) + 90*pow2(lmMsq)
        - 48*pow2(lmMst1)))/27. + (64*pow2(lmMst1))/9. - (4*(15 + 30*lmMsq +
        232*lmMst1)*pow2(lmMst2))/9. + (20*(143 - 72*lmMsq*(2 + lmMst1 -
        lmMst2) + 132*lmMst2 + 12*lmMst1*(1 + 6*lmMst2) - 72*pow2(lmMst2))*
        pow2(Mst1))/(27.*pow2(Msq)) + (320*pow3(lmMst2))/3. + (5*(23 + 12*lmMsq
        - 12*lmMst2)*pow4(Mst1))/(108.*pow4(Msq))) + (1560*(236*OepS2 + 27*(32
        - 177*lmMst1 + 177*lmMst2)*S2)*pow2(Mst1)*pow2(Mst2) + (-622840*OepS2 +
        27*(836021 + 467130*lmMst1 - 467130*lmMst2)*S2)*pow4(Mst1) + 9*(280*
        OepS2 - 81*(-1387 + 70*lmMst1 - 70*lmMst2)*S2)*pow4(Mst2))/(21870.*
        Mst2) + (8*((5*(397 + 36*lmMsq*(-2 + lmMst1 - lmMst2) + 78*lmMst2 - 6*
        lmMst1*(1 + 6*lmMst2) + 36*pow2(lmMst2)))/pow2(Msq) + (48*(2 + lmMst2 -
        3*pow2(lmMst2)))/pow2(Mst1))*pow5(Mst2))/27. + (5*(9*(-593 + 4*lmMsq*(
        181 + 90*lmMst1 - 90*lmMst2) - 718*lmMst2 - 6*lmMst1*(1 + 60*lmMst2) +
        360*pow2(lmMst2))*pow2(Mst1)*pow5(Mst2) + (21275 + 12*lmMsq*(-541 +
        270*lmMst1 - 270*lmMst2) + 6546*lmMst2 - 54*lmMst1*(1 + 60*lmMst2) +
        3240*pow2(lmMst2))*pow7(Mst2)))/(324.*pow4(Msq))) + pow2(Mt)*pow2(s2t)*
        (pow2(Mst1)*(330.8944644326867 - (760*lmMsq)/9. + (4*(211261 + 3375*
        lmMsq)*lmMst1)/2025. + (32*(6 + lmMst1 - lmMst2)*lmMt)/9. - 30*pow2(
        lmMsq) + (2*lmMst2*(548953 + 54000*lmMsq + 174270*lmMst1 - 450*pow2(
        lmMst1)))/2025. - (578*pow2(lmMst1))/135. - (2*(15049 + 8610*lmMst1)*
        pow2(lmMst2))/135. + (4*pow3(lmMst1))/27. + (3452*pow3(lmMst2))/27.) +
        ((40*(4 - lmMst1 + lmMst2))/(9.*pow2(Msq)) - (190.3527287429963 - (160*
        lmMsq)/3. + (-535.9278609221467 + (200*lmMsq)/3.)*lmMst1 + lmMst2*(
        315.4834164777022 - (200*lmMsq)/3. + (84592*lmMst1)/315. - (208*pow2(
        lmMst1))/9.) - (35576*pow2(lmMst1))/315. + (8*(-6127 + 5110*lmMst1)*
        pow2(lmMst2))/315. + (64*lmMt*(2 + 5*lmMst2 - lmMst1*(5 + 2*lmMst2) +
        pow2(lmMst1) + pow2(lmMst2)))/3. + (16*pow3(lmMst1))/27. - (2896*pow3(
        lmMst2))/27.)/pow2(Mst2))*pow4(Mst1) + pow2(Mst2)*(1427.379365079365 -
        (16*B4)/3. + (16*D3)/3. - (8*DN)/3. + (910*lmMsq)/9. + (608*lmMst1)/9.
         - (32*lmMt)/
        9. + 60*pow2(lmMsq) + (4*lmMst2*(5987 - 810*lmMsq + 603*lmMst1 - 144*
        pow2(lmMst1)))/27. + (58*pow2(lmMst1))/3. + (2*(743 - 96*lmMst1)*pow2(
        lmMst2))/9. + (10*(49 - 100*lmMsq + 100*lmMst2)*pow2(Mst1))/(9.*pow2(
        Msq)) + (128*pow3(lmMst2))/3. + (5*(1 + 4*lmMsq - 4*lmMst2)*pow4(Mst1))
        /(6.*pow4(Msq))) - (((10*(313 - 600*lmMsq + 600*lmMst2))/pow2(Msq) + (
        3*(-1085 - 118*lmMst2 + 64*lmMst1*(-2 + 3*lmMst2) - 90*lmMsq*(-5 + 6*
        lmMst2) + 270*pow2(lmMsq) + 944*pow2(lmMst2)))/pow2(Mst1))*pow4(Mst2))/
        27. - S2*((807.2952380952381 + 796*lmMst1 - 796*lmMst2)*pow2(Mst1) + (
        596.0571428571428 + 84*lmMst1 - 84*lmMst2)*pow2(Mst2) + (4*(28283 +
        23070*lmMst1 - 23070*lmMst2)*pow4(Mst1))/(45.*pow2(Mst2)) - (pow2(
        MuSUSY)*(6*(1089707 - 5550*lmMst1 + 5550*lmMst2)*pow2(Mst1)*pow2(Mst2)
        + 40*(93919 - 12981*lmMst1 + 12981*lmMst2)*pow4(Mst1) + 27*(116129 -
        350*lmMst1 + 350*lmMst2)*pow4(Mst2)))/(405.*pow4(Mst2))) + (64*(-2 +
        lmMst2 + 5*pow2(lmMst2))*pow6(Mst2))/(9.*pow4(Mst1)) + (8*OepS2*(4*(
        6921*pow2(Mst2) + 4327*pow2(MuSUSY))*pow4(Mst1) + 6*pow2(Mst1)*(185*
        pow2(Mst2)*pow2(MuSUSY) + 1791*pow4(Mst2)) + 63*(5*pow2(MuSUSY)*pow4(
        Mst2) + 18*pow6(Mst2))))/(2187.*pow4(Mst2)) - pow2(MuSUSY)*(
        112.1664609053498 - (100*B4)/9. + (112*D3)/9. - (62*DN)/9. + (1055*
        lmMsq)/9. + 15*pow2(lmMsq) - (2*lmMst1*(-544 + 210*lmMsq + 135*pow2(
        lmMsq)))/9. - (103*pow2(lmMst1))/9. + (lmMst2*(-7921 - 1476*lmMst1 +
        90*lmMsq*(5 + 18*lmMst1) + 810*pow2(lmMsq) + 288*pow2(lmMst1)))/27. + (
        93.66666666666667 - 60*lmMsq - 112*lmMst1)*pow2(lmMst2) + (10*(405 + (-
        187 + 246*lmMsq)*lmMst2 + lmMst1*(187 - 246*lmMsq + 246*lmMst2) - 246*
        pow2(lmMst2))*pow2(Mst1))/(27.*pow2(Msq)) + ((-1149 + 266*lmMst2 + 64*
        lmMst1*(-2 + 3*lmMst2) - 90*lmMsq*(-5 + 6*lmMst2) + 270*pow2(lmMsq) +
        1264*pow2(lmMst2))*pow2(Mst2))/(18.*pow2(Mst1)) - (pow2(Mst1)*(
        73.78472976680384 + (164*B4)/9. - (176*D3)/9. + (94*DN)/9. - (260*
        lmMsq)/3. - 20*lmMsq*(2 + 3*lmMst1)*lmMst2 - 30*lmMst2*pow2(lmMsq) +
        lmMst1*(-96.55703703703703 + 40*lmMsq + 30*pow2(lmMsq)) + (7556*pow2(
        lmMst1))/135. - (2*lmMst2*(-99788 + 2005*lmMst1 + 32775*pow2(lmMst1)))/
        675. + (2*(-3377 + 4050*lmMsq + 19155*lmMst1)*pow2(lmMst2))/135. + (10*
        pow3(lmMst1))/27. - (5050*pow3(lmMst2))/27.))/pow2(Mst2) + (304*pow3(
        lmMst2))/3. + pow4(Mst1)*((5*(216 + (-95 + 132*lmMsq)*lmMst2 + lmMst1*(
        95 - 132*lmMsq + 132*lmMst2) - 132*pow2(lmMst2)))/(54.*pow4(Msq)) + (
        536.1152102791342 - (8*B4)/3. + (32*D3)/9. - (20*DN)/9. + 90*lmMsq - (
        123.11224321827497 + 20*lmMsq*(1 + lmMsq))*lmMst1 - lmMst2*(
        17.33220122616948 - 20*lmMsq*(1 + lmMsq) + (133.04550264550264 - 40*
        lmMsq)*lmMst1 - (1180*pow2(lmMst1))/9.) - (15886*pow2(lmMst1))/945. + (
        149.85608465608465 - 40*lmMsq - (2812*lmMst1)/9.)*pow2(lmMst2) - (92*
        pow3(lmMst1))/27. + (4988*pow3(lmMst2))/27.)/pow4(Mst2)) - (32*(-2 +
        lmMst2 + 5*pow2(lmMst2))*pow4(Mst2))/(9.*pow4(Mst1)) + (5*((291 + 2*(-
        103 + 84*lmMsq)*lmMst2 + 2*lmMst1*(103 - 84*lmMsq + 84*lmMst2) - 168*
        pow2(lmMst2))*pow4(Mst1) + (684 + 347*lmMst1 + 7*(-77 + 78*lmMst1)*
        lmMst2 + 6*lmMsq*(32 - 91*lmMst1 + 91*lmMst2) - 546*pow2(lmMst2))*pow4(
        Mst2)))/(27.*pow2(Msq)*pow2(Mst2)) + (5*((3834 + (-890 + 2328*lmMsq)*
        lmMst2 + lmMst1*(890 - 2328*lmMsq + 2328*lmMst2) - 2328*pow2(lmMst2))*
        pow2(Mst2)*pow4(Mst1) + 160*(-8 + 15*lmMsq - 15*lmMst2)*pow2(Msq)*pow4(
        Mst2) + (4239 - 1324*lmMst2 + 8*lmMst1*(95 + 366*lmMst2) + lmMsq*(564 -
        2928*lmMst1 + 2928*lmMst2) - 2928*pow2(lmMst2))*pow2(Mst1)*pow4(Mst2) +
        63*(-5 + 28*lmMsq - 28*lmMst2)*pow6(Mst2)))/(216.*pow2(Mst1)*pow4(Msq))
        ) - (5*((-31 + 1980*lmMsq - 1980*lmMst2)*pow4(Mst1)*pow4(Mst2) + 160*(-
        8 + 15*lmMsq - 15*lmMst2)*pow2(Msq)*pow6(Mst2) + (481 - 3636*lmMsq +
        3636*lmMst2)*pow2(Mst1)*pow6(Mst2) + 9*(-9 + 140*lmMsq - 140*lmMst2)*
        pow8(Mst2)))/(108.*pow2(Mst1)*pow4(Msq))) + pow4(s2t)*(-(pow2(Mst1)*
        pow2(Mst2)*(90.4877462277092 - B4 + (4*D3)/3. - (5*DN)/6. + (1105*
        lmMsq)/36. + lmMst1*(38.08296296296297 - (40*lmMsq)/3. - (15*pow2(
        lmMsq))/2.) + (15*pow2(lmMsq))/4. + (2233*pow2(lmMst1))/270. - (lmMst2*
        (206449 + 85010*lmMst1 - 2250*lmMsq*(7 + 18*lmMst1) - 20250*pow2(lmMsq)
        + 51150*pow2(lmMst1)))/2700. + (16.77037037037037 - 15*lmMsq + (269*
        lmMst1)/18.)*pow2(lmMst2) + (5*(25 + 12*lmMsq*(-4 + 3*lmMst1 - 3*
        lmMst2) + (48 - 36*lmMst1)*lmMst2 + 36*pow2(lmMst2))*pow2(Mst1))/(108.*
        pow2(Msq)) + (5*pow3(lmMst1))/54. + (211*pow3(lmMst2))/54.)) + (
        84.40344866031853 - (5*lmMsq*(2 + 3*lmMst1))/3. + (1011403*lmMst2)/
        1.1907e6 + (5*pow2(lmMsq))/2. + (4.326984126984127 + (11*lmMst2)/3.)*
        pow2(lmMst1) + lmMst1*(10.650581170739901 + (1474*lmMst2)/315. - (11*
        pow2(lmMst2))/3.) - (113*pow2(lmMst2))/35. - (11*pow3(lmMst1))/9. + (
        11*pow3(lmMst2))/9.)*pow4(Mst1) + (61.73605967078189 - (25*B4)/9. + (
        28*D3)/9. - (31*DN)/18. + (605*lmMsq)/36. - (15*pow2(lmMsq))/4. - (
        lmMst1*(-608 + 210*lmMsq + 135*pow2(lmMsq)))/18. - (103*pow2(lmMst1))/
        36. + (lmMst2*(-8815 - 2052*lmMst1 + 90*lmMsq*(23 + 18*lmMst1) + 810*
        pow2(lmMsq) + 288*pow2(lmMst1)))/108. - (16.13888888888889 + 15*lmMsq +
        28*lmMst1)*pow2(lmMst2) - (5*(359 - 202*lmMst2 + 20*lmMst1*(8 + 15*
        lmMst2) + 6*lmMsq*(7 - 50*lmMst1 + 50*lmMst2) - 300*pow2(lmMst2))*pow2(
        Mst1))/(54.*pow2(Msq)) + (76*pow3(lmMst2))/3. + (5*(113 - 34*lmMst2 +
        12*lmMsq*(1 - 10*lmMst1 + 10*lmMst2) + 2*lmMst1*(11 + 60*lmMst2) - 120*
        pow2(lmMst2))*pow4(Mst1))/(432.*pow4(Msq)))*pow4(Mst2) - (2*OepS2*(480*
        pow2(Mst1)*pow2(Mst2) + 8398*pow4(Mst1) + 315*pow4(Mst2)))/2187. + (S2*
        (12*(-22273 + 1200*lmMst1 - 1200*lmMst2)*pow2(Mst1)*pow2(Mst2) + 10*(
        51635 + 25194*lmMst1 - 25194*lmMst2)*pow4(Mst1) + 27*(-116129 + 350*
        lmMst1 - 350*lmMst2)*pow4(Mst2)))/1620. + (((10*(1004 - 6*lmMsq*(68 +
        91*lmMst1 - 91*lmMst2) + 61*lmMst2 + lmMst1*(347 + 546*lmMst2) - 546*
        pow2(lmMst2)))/pow2(Msq) + (3*(-1405 + 394*lmMst2 + 64*lmMst1*(-2 + 3*
        lmMst2) - 90*lmMsq*(-5 + 6*lmMst2) + 270*pow2(lmMsq) + 1904*pow2(
        lmMst2)))/pow2(Mst1))*pow6(Mst2))/216. + ((8*(2 + lmMst2 - 3*pow2(
        lmMst2))*pow2(Msq) + 5*(-4 + 15*lmMsq - 15*lmMst2)*pow2(Mst1))*pow8(
        Mst2))/(9.*pow2(Msq)*pow4(Mst1)) + (5*(3*(-1653 + 4*lmMsq*(53 + 294*
        lmMst1 - 294*lmMst2) - 2*lmMst2 - 42*lmMst1*(5 + 28*lmMst2) + 1176*
        pow2(lmMst2))*pow2(Mst1)*pow6(Mst2) + (3287 - 180*lmMsq*(11 + 12*lmMst1
        - 12*lmMst2) + 1644*lmMst2 + 48*lmMst1*(7 + 45*lmMst2) - 2160*pow2(
        lmMst2))*pow8(Mst2)))/(864.*pow4(Msq))) + (MuSUSY*(pow2(Mt)*pow2(s2t)*(
        -((pow2(Mst1)*(1980.3992798353909 + (1048*B4)/3. - (20*DN)/3. - 1440*
        lmMsq - lmMst1*(363.037037037037 + 200*lmMsq + 40*pow2(lmMsq)) - (128*
        pow2(lmMst1))/3. + lmMst2*(40*lmMsq*(5 + lmMsq + 2*lmMst1) + (2*(35429
        + 3420*lmMst1 + 288*pow2(lmMst1)))/27.) - (8*(79 + 30*lmMsq + 168*
        lmMst1)*pow2(lmMst2))/3. - (80*(15 + lmMst1*(-4 + 6*lmMsq - 6*lmMst2) +
        (4 - 6*lmMsq)*lmMst2 + 6*pow2(lmMst2))*pow2(Mst1))/(3.*pow2(Msq)) + (
        1280*pow3(lmMst2))/3.))/Mst2) + (8*((5*(397 + 36*lmMsq*(-2 + lmMst1 -
        lmMst2) + 78*lmMst2 - 6*lmMst1*(1 + 6*lmMst2) + 36*pow2(lmMst2)))/pow2(
        Msq) + (48*(2 + lmMst2 - 3*pow2(lmMst2)))/pow2(Mst1))*pow3(Mst2))/9. +
        ((618.5941700960219 - (592*B4)/3. + (8*DN)/3. + 750*lmMsq + (20*lmMst1*
        (542 + 243*lmMsq + 162*pow2(lmMsq)))/81. + (92*pow2(lmMst1))/3. -
        lmMst2*(20*lmMsq*(3 + 4*lmMst1) + 40*pow2(lmMsq) + (8*(11723 - 702*
        lmMst1 + 756*pow2(lmMst1)))/81.) + (4*(-75 + 60*lmMsq + 344*lmMst1)*
        pow2(lmMst2))/3. - (32*pow3(lmMst1))/3. - (1120*pow3(lmMst2))/3.)*pow4(
        Mst1))/pow3(Mst2) - Mst2*(219.11327160493826 + (980*B4)/3. - (22*DN)/3.
         - (3940*lmMsq)/3. + lmMst1*(239.33333333333334 - 100*lmMsq - 20*pow2(
        lmMsq)) + 40*pow2(lmMsq) + lmMst2*(1932.2222222222222 + 20*lmMsq*(1 +
        lmMsq) + 4*(37 + 10*lmMsq)*lmMst1 - (32*pow2(lmMst1))/3.) + (64*pow2(
        lmMst1))/3. - (4*(-81 + 30*lmMsq + 232*lmMst1)*pow2(lmMst2))/3. - (20*(
        217 + 8*lmMst1*(-1 + 6*lmMsq - 6*lmMst2) + (8 - 48*lmMsq)*lmMst2 + 48*
        pow2(lmMst2))*pow2(Mst1))/(3.*pow2(Msq)) + 320*pow3(lmMst2) - ((
        141.34259259259258 - (5*lmMsq)/9. + 5*(-7 + 20*lmMsq)*lmMst1 + (
        35.55555555555556 - 100*(lmMsq + lmMst1))*lmMst2 + 100*pow2(lmMst2))*
        pow4(Mst1))/pow4(Msq)) + ((-4*OepS2*(-9267*pow2(Mst1)*pow2(Mst2) +
        23734*pow4(Mst1) - 63*pow4(Mst2)))/729. + (S2*(3*(29123 - 92670*lmMst1
        + 92670*lmMst2)*pow2(Mst1)*pow2(Mst2) + 2*(525961 + 356010*lmMst1 -
        356010*lmMst2)*pow4(Mst1) + 27*(1387 - 70*lmMst1 + 70*lmMst2)*pow4(
        Mst2)))/270.)/pow3(Mst2) + (5*(2*(7969 + 12*lmMsq*(1 + 270*lmMst1 -
        270*lmMst2) + 42*lmMst2 - 54*lmMst1*(1 + 60*lmMst2) + 3240*pow2(lmMst2)
        )*pow2(Mst1)*pow3(Mst2) + (21275 + 12*lmMsq*(-541 + 270*lmMst1 - 270*
        lmMst2) + 6546*lmMst2 - 54*lmMst1*(1 + 60*lmMst2) + 3240*pow2(lmMst2))*
        pow5(Mst2)))/(108.*pow4(Msq))) - s2t*pow3(Mt)*(3063.8698412698413 - (
        32*B4)/3. + (32*D3)/3. - (16*DN)/3. + (920*lmMsq)/9. + (1472*lmMst1)/9.
         - (64*lmMt)/9. + 60*pow2(lmMsq) + (4*lmMst2*(12247 - 810*lmMsq + 918*
        lmMst1 - 288*pow2(lmMst1)))/27. + (116*pow2(lmMst1))/3. + (4*(431 - 96*
        lmMst1)*pow2(lmMst2))/9. - (40*pow2(Mst1))/(9.*pow2(Msq)) + (2*(957 +
        64*lmMst1*(2 - 3*lmMst2) + 182*lmMst2 + 90*lmMsq*(-5 + 6*lmMst2) - 270*
        pow2(lmMsq) - 624*pow2(lmMst2))*pow2(Mst2))/(9.*pow2(Mst1)) + (256*
        pow3(lmMst2))/3. + (pow2(Mst1)*(3746.992103468548 - (32*B4)/3. + (32*
        D3)/3. - (16*DN)/3. - (200*lmMsq)/3. + (8*(253561 + 3375*lmMsq)*lmMst1)
        /2025. + (64*(5 + lmMst1 - lmMst2)*lmMt)/9. + (8*lmMst2*(732839 - 3375*
        lmMsq + 121560*lmMst1 - 11025*pow2(lmMst1)))/2025. + (4064*pow2(lmMst1)
        )/135. - (8*(4292 + 5025*lmMst1)*pow2(lmMst2))/135. + (8*pow3(lmMst1))/
        27. + (9208*pow3(lmMst2))/27.))/pow2(Mst2) - (5*pow2(Mst2)*(27*(5 + 4*
        lmMsq - 4*lmMst2)*pow2(Mst1) + 2*(83 - 936*lmMsq + 936*lmMst2)*pow2(
        Mst2)))/(54.*pow4(Msq)) - (20*pow4(Mst1))/(9.*pow4(Msq)) + ((
        772.4522050654477 - (64*B4)/9. + (64*D3)/9. - (32*DN)/9. + (560*lmMsq)/
        3. + (1513.2107835726883 - 160*lmMsq)*lmMst1 + (38704*pow2(lmMst1))/
        189. + lmMst2*(160*lmMsq + (4*(1947623 - 3019800*lmMst1 + 97020*pow2(
        lmMst1)))/19845.) + (16*(6787 - 5439*lmMst1)*pow2(lmMst2))/189. - (128*
        lmMt*(4 + 6*lmMst2 - 2*lmMst1*(3 + lmMst2) + pow2(lmMst1) + pow2(
        lmMst2)))/3. - (16*pow3(lmMst1))/9. + (1328*pow3(lmMst2))/3.)*pow4(
        Mst1))/pow4(Mst2) + (128*(-2 + lmMst2 + 5*pow2(lmMst2))*pow4(Mst2))/(9.
        *pow4(Mst1)) + (20*(4*(3 - lmMst1 + lmMst2)*pow4(Mst1) + (-51 + 100*
        lmMsq - 100*lmMst2)*pow4(Mst2)))/(9.*pow2(Msq)*pow2(Mst2)) + ((32*
        OepS2*(660*pow2(Mst1)*pow2(Mst2) + 1838*pow4(Mst1) + 63*pow4(Mst2)))/
        243. - (4*S2*(12*(18419 + 11550*lmMst1 - 11550*lmMst2)*pow2(Mst1)*pow2(
        Mst2) + 28*(17269 + 13785*lmMst1 - 13785*lmMst2)*pow4(Mst1) + 27*(3477
        + 490*lmMst1 - 490*lmMst2)*pow4(Mst2)))/315.)/pow4(Mst2) + (800*(8 -
        15*lmMsq + 15*lmMst2)*pow2(Msq)*pow4(Mst2) + 315*(5 - 28*lmMsq + 28*
        lmMst2)*pow6(Mst2))/(54.*pow2(Mst1)*pow4(Msq))) + (pow4(Mt)*(-(pow4(
        Msq)*(-3*pow2(Mst1)*pow2(Mst2)*(34872121 - 13890240*lmMst2 - 7840*OepS2
        - 324*(9361 + 490*lmMst2)*S2 - 453600*(lmMst1 - lmMst2)*pow2(lmMsq) -
        241920*(-2 + lmMst2)*pow2(lmMst1) + 7560*lmMst1*(206 - 128*lmMt +
        lmMst2*(636 + 64*lmMt) + 21*S2 - 224*pow2(lmMst2)) + 50400*lmMsq*(64 +
        48*lmMst2 + 9*lmMst1*(-5 + 2*lmMst2) - 3*lmMt - 18*pow2(lmMst2)) -
        5654880*pow2(lmMst2) - 5040*lmMt*(-31 - 240*lmMst2 + 96*pow2(lmMst2)) +
        120960*pow2(lmMt) + 1935360*pow3(lmMst2)) - 8*(24896293 + 2509920*
        lmMst2 - 605920*OepS2 - 108*(105619 + 113610*lmMst2)*S2 - 170100*(
        lmMst1 - lmMst2)*pow2(lmMsq) + 90720*(-3 + lmMst2)*pow2(lmMst1) +
        18900*lmMsq*(100 + 48*lmMst2 + 9*lmMst1*(-5 + 2*lmMst2) - 3*lmMt - 18*
        pow2(lmMst2)) - 1563030*pow2(lmMst2) - 1890*lmMt*(-539 - 233*lmMst2 +
        96*pow2(lmMst2)) - 1890*lmMst1*(2075 + 281*lmMt - lmMst2*(1139 + 96*
        lmMt) - 6492*S2 + 528*pow2(lmMst2)) - 226800*pow2(lmMt) + 907200*pow3(
        lmMst2))*pow4(Mst1) + 2903040*(2 + lmMst2 - 3*pow2(lmMst2))*pow4(Mst2))
        ) + 50400*pow2(Msq)*((778 - 3*lmMst1*(1 + 78*lmMsq - 78*lmMst2) + 6*(1
        + 39*lmMsq)*lmMst2 - 3*lmMt - 234*pow2(lmMst2))*pow2(Mst2)*pow4(Mst1) +
        3*(265 - 34*lmMst2 + 12*lmMst1*(1 + 6*lmMst2) - 6*(-2 + lmMst2)*lmMt +
        2*lmMsq*(5 - 36*lmMst1 + 39*lmMst2 + 3*lmMt) - 6*pow2(lmMsq) - 72*pow2(
        lmMst2))*pow2(Mst1)*pow4(Mst2) + (313 + 42*(-2 + 3*lmMsq)*lmMst2 + 3*
        lmMst1*(23 - 42*lmMsq + 42*lmMst2) + 15*lmMt - 126*pow2(lmMst2))*pow6(
        Mst1)) + 3150*(4*(1330 + 12*(-4 + 33*lmMsq)*lmMst2 + lmMst1*(51 - 396*
        lmMsq + 396*lmMst2) - 3*lmMt - 396*pow2(lmMst2))*pow2(Mst2)*pow6(Mst1)
        + 3*((4585 + 104*lmMst2 + 6*lmMst1*(-25 + 228*lmMst2) - 4*lmMsq*(5 +
        342*lmMst1 - 348*lmMst2 - 6*lmMt) + 66*lmMt - 24*lmMst2*lmMt - 24*pow2(
        lmMsq) - 1368*pow2(lmMst2))*pow4(Mst1)*pow4(Mst2) + (4317 - 280*lmMst2
        + 18*lmMst1*(1 + 60*lmMst2) - 78*(-3 + 4*lmMst2)*lmMt + 4*lmMsq*(7 -
        270*lmMst1 + 348*lmMst2 + 78*lmMt) - 312*pow2(lmMsq) - 1080*pow2(
        lmMst2))*pow2(Mst1)*pow6(Mst2)))))/(51030.*pow2(Mst1)*pow3(Mst2)*pow4(
        Msq)) + Mt*((2*T1ep*(3*Mst2*pow2(Mst1)*(15840*Mst2*s2t*pow2(Mt) - 9267*
        Mt*pow2(Mst2)*pow2(s2t) + 17312*pow3(Mt) - 530*pow3(Mst2)*pow3(s2t)) +
        63*pow3(Mst2)*(72*Mst2*s2t*pow2(Mt) - 3*Mt*pow2(Mst2)*pow2(s2t) + 4*
        pow3(Mt) - 10*pow3(Mst2)*pow3(s2t)) + 2*s2t*(35601*Mst2*Mt*s2t + 66168*
        pow2(Mt) - 12664*pow2(Mst2)*pow2(s2t))*pow4(Mst1)))/(729.*pow4(Mst2)) +
        pow3(s2t)*(pow2(Mst1)*(185.95119067215364 + (64*B4)/9. - (64*D3)/9. + (
        32*DN)/9. + (275*lmMsq)/9. + (4*(4106 - 1125*lmMsq)*lmMst1)/675. + 15*
        pow2(lmMsq) + (6011*pow2(lmMst1))/135. - (lmMst2*(-1551 + 15750*lmMsq +
        40910*lmMst1 + 58350*pow2(lmMst1)))/675. + ((5891 + 23190*lmMst1)*pow2(
        lmMst2))/135. + (10*pow3(lmMst1))/27. - (2314*pow3(lmMst2))/27.) - ((5*
        (5 - 2*lmMst1 + 2*lmMst2))/(9.*pow2(Msq)) + (426.99721850958684 + (
        1510169*lmMst2)/59535. + lmMsq*(10 - 20*lmMst1 + 20*lmMst2) + (2*(433 +
        6552*lmMst2)*pow2(lmMst1))/189. + (194*pow2(lmMst2))/189. - (lmMst1*(
        2674409 + 333900*lmMst2 + 7514640*pow2(lmMst2)))/59535. - (112*pow3(
        lmMst1))/27. + (1648*pow3(lmMst2))/27.)/pow2(Mst2))*pow4(Mst1) - pow2(
        Mst2)*(175.99979423868314 - (100*B4)/9. + (112*D3)/9. - (62*DN)/9. + (
        830*lmMsq)/9. - (2*lmMst1*(-192 + 70*lmMsq + 45*pow2(lmMsq)))/3. - (
        103*pow2(lmMst1))/9. + (2*lmMst2*(-4160 - 882*lmMst1 + 90*lmMsq*(7 + 9*
        lmMst1) + 405*pow2(lmMsq) + 144*pow2(lmMst1)))/27. + (
        23.444444444444443 - 60*lmMsq - 112*lmMst1)*pow2(lmMst2) + (5*(42 +
        lmMst1*(9 - 18*lmMst2) + 2*lmMsq*(-32 + 9*lmMst1 - 9*lmMst2) + 55*
        lmMst2 + 18*pow2(lmMst2))*pow2(Mst1))/(9.*pow2(Msq)) + (304*pow3(
        lmMst2))/3. - (5*pow4(Mst1))/(18.*pow4(Msq))) - (((10*(844 - 6*lmMsq*(
        18 + 91*lmMst1 - 91*lmMst2) - 239*lmMst2 + lmMst1*(347 + 546*lmMst2) -
        546*pow2(lmMst2)))/pow2(Msq) + (3*(-1277 + 330*lmMst2 + 64*lmMst1*(-2 +
        3*lmMst2) - 90*lmMsq*(-5 + 6*lmMst2) + 270*pow2(lmMsq) + 1584*pow2(
        lmMst2)))/pow2(Mst1))*pow4(Mst2))/54. + (159*(200*OepS2 + 27*(21401 -
        150*lmMst1 + 150*lmMst2)*S2)*pow2(Mst1)*pow2(Mst2) + 20*(25328*OepS2 +
        27*(47051 - 18996*lmMst1 + 18996*lmMst2)*S2)*pow4(Mst1) + 9*(1400*OepS2
        + 81*(116129 - 350*lmMst1 + 350*lmMst2)*S2)*pow4(Mst2))/(10935.*pow2(
        Mst2)) + (32*(-2 + lmMst2 + 5*pow2(lmMst2))*pow6(Mst2))/(9.*pow4(Mst1))
        + (5*((405 - 434*lmMst2 + 10*lmMst1*(-13 + 60*lmMst2) + lmMsq*(564 -
        600*lmMst1 + 600*lmMst2) - 600*pow2(lmMst2))*pow4(Mst1)*pow4(Mst2) -
        160*(-8 + 15*lmMsq - 15*lmMst2)*pow2(Msq)*pow6(Mst2) + 2*(-2277 + 24*
        lmMsq*(25 + 61*lmMst1 - 61*lmMst2) - 220*lmMst2 - 4*lmMst1*(95 + 366*
        lmMst2) + 1464*pow2(lmMst2))*pow2(Mst1)*pow6(Mst2) + 9*(9 - 140*lmMsq +
        140*lmMst2)*pow8(Mst2)))/(216.*pow2(Mst1)*pow4(Msq))))))/Tbeta)))/pow2(
        Mgl) + (Al4p*xDmglst2*z2*pow2(Dmglst2)*((twoLoopFlag*(4*pow2(Mt)*(8*(4*
        MuSUSY + 9*Mst2*Tbeta)*pow2(Mt)*pow2(Sbeta) + 3*Mst2*MuSUSY*pow2(s2t)*(
        -(MuSUSY*Tbeta*(-1 + pow2(Sbeta))) + 10*Mst2*pow2(Sbeta)) + 2*Mt*s2t*
        Tbeta*(31*pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 12*pow2(Mst2)*pow2(Sbeta)))
        *pow4(Mst1) + pow2(Mst1)*pow2(Mst2)*(12*Mst2*MuSUSY*pow2(Mt)*pow2(s2t)*
        (MuSUSY*Tbeta + 10*Mst2*pow2(Sbeta) - MuSUSY*Tbeta*pow2(Sbeta)) + 168*
        s2t*Tbeta*pow2(MuSUSY)*(-1 + pow2(Sbeta))*pow3(Mt) + 4*Mt*Tbeta*pow2(
        Sbeta)*pow3(s2t)*pow4(Mst2) + 32*(MuSUSY + Mst2*Tbeta)*pow2(Sbeta)*
        pow4(Mt) - 3*Tbeta*pow2(Sbeta)*pow4(s2t)*pow5(Mst2)) + pow4(Mst2)*(12*
        Mst2*MuSUSY*pow2(Mt)*pow2(s2t)*(-(MuSUSY*Tbeta*(-1 + pow2(Sbeta))) +
        11*Mst2*pow2(Sbeta)) + 8*s2t*Tbeta*(11*pow2(MuSUSY)*(-1 + pow2(Sbeta))
        - 4*pow2(Mst2)*pow2(Sbeta))*pow3(Mt) - 4*Mt*(3*MuSUSY + 11*Mst2*Tbeta)*
        pow2(Sbeta)*pow3(Mst2)*pow3(s2t) + 32*MuSUSY*pow2(Sbeta)*pow4(Mt) + 3*
        Tbeta*pow2(Sbeta)*pow4(s2t)*pow5(Mst2))))/(3.*Tbeta*pow2(Sbeta)*pow5(
        Mst2)) + (Al4p*threeLoopFlag*((7*pow4(s2t)*(pow4(Msq)*(12*(-34574 +
        36450*lmMsq + 15795*lmMst1 - 87075*lmMst2)*pow2(Mst2)*pow4(Mst1) + 9*(
        64091 - 48600*lmMsq - 53280*lmMst1 + 148320*lmMst2)*pow2(Mst1)*pow4(
        Mst2) + (4400182 - 427680*lmMst1 + 427680*lmMst2)*pow6(Mst1) + 233280*
        pow6(Mst2)) - 8100*pow2(Msq)*(-100*pow4(Mst1)*pow4(Mst2) + 6*pow2(Mst2)
        *pow6(Mst1) + 91*pow2(Mst1)*pow6(Mst2)) - 4050*(10*pow4(Mst2)*pow6(
        Mst1) - 147*pow4(Mst1)*pow6(Mst2) + 90*pow2(Mst1)*pow8(Mst2))))/(pow2(
        Mst1)*pow4(Msq)) + (2*Mt*(4*pow3(Mt)*(MuSUSY*(-48*(-1161746 + 56700*
        lmMsq + 6615*lmMst1 - 222075*lmMst2 + 34020*lmMt)*Mst2*pow2(Mst1) - 9*(
        718271 + 302400*lmMsq + 32760*lmMst1 - 955080*lmMst2 - 45360*lmMt)*
        pow3(Mst2)) + Tbeta*(30*(6261919 + 45360*lmMsq - 240408*lmMst1 +
        406728*lmMst2 - 287280*lmMt)*pow2(Mst1)*pow2(Mst2) - (756*(5*(-3643 +
        120*lmMst1 - 3192*lmMst2)*pow2(Mst1) - 8*(617 + 2040*lmMst2)*pow2(Mst2)
        )*pow2(MuSUSY)*(-1 + pow2(Sbeta)))/pow2(Sbeta) - 16*(17026487 +
        1530900*lmMsq - 3163860*lmMst1 + 4320540*lmMst2 - 340200*lmMt)*pow4(
        Mst1) - 630*(-19211 + 648*lmMst1 + 5328*lmMst2 + 7848*lmMt)*pow4(Mst2)
        + (3265920*pow6(Mst2))/pow2(Mst1)) - (226800*(22*MuSUSY*pow3(Mst2)*
        pow4(Mst1) + 58*MuSUSY*pow2(Mst1)*pow5(Mst2) + pow2(Msq)*(52*MuSUSY*
        pow2(Mst1)*pow3(Mst2) + Mst2*(28*MuSUSY + 6*Mst2*Tbeta)*pow4(Mst1) + (
        52*MuSUSY - 9*Mst2*Tbeta)*pow5(Mst2)) + (58*MuSUSY - 17*Mst2*Tbeta)*
        pow7(Mst2)))/pow4(Msq)) + (s2t*(7*pow2(Mst2)*pow2(s2t)*(-(pow4(Msq)*(6*
        ((131737 - 96660*lmMst1 + 96660*lmMst2)*MuSUSY + 10*(-76585 + 3240*
        lmMsq - 36504*lmMst1 + 29592*lmMst2)*Mst2*Tbeta)*pow2(Mst2)*pow4(Mst1)
        + 9*(2*(90011 - 48600*lmMsq - 53280*lmMst1 + 148320*lmMst2)*MuSUSY + (
        1340537 - 496800*lmMsq + 4680*lmMst1 + 778680*lmMst2)*Mst2*Tbeta)*pow2(
        Mst1)*pow4(Mst2) + (16*(776573 - 80595*lmMst1 + 80595*lmMst2)*MuSUSY +
        (8194483 - 2329560*lmMst1 + 2329560*lmMst2)*Mst2*Tbeta)*pow6(Mst1) +
        466560*MuSUSY*pow6(Mst2))) + 16200*pow2(Msq)*(-3*(3*MuSUSY + 16*Mst2*
        Tbeta)*pow4(Mst1)*pow4(Mst2) + 13*(7*MuSUSY + 24*Mst2*Tbeta)*pow2(Mst1)
        *pow6(Mst2)) - 8100*(5*(5*MuSUSY + 18*Mst2*Tbeta)*pow4(Mst1)*pow6(Mst2)
        - 2*(61*MuSUSY + 179*Mst2*Tbeta)*pow2(Mst1)*pow8(Mst2))) + (Mt*(-56700*
        MuSUSY*s2t*pow2(Mst1)*pow2(Mst2)*(2*pow2(Mst2)*(-11*MuSUSY*Tbeta*(-1 +
        pow2(Sbeta)) + 78*Mst2*pow2(Sbeta))*pow4(Mst1) + pow2(Mst1)*(-97*
        MuSUSY*Tbeta*(-1 + pow2(Sbeta)) + 804*Mst2*pow2(Sbeta))*pow4(Mst2) + 2*
        pow2(Msq)*(2*pow2(Mst1)*pow2(Mst2)*(-41*MuSUSY*Tbeta*(-1 + pow2(Sbeta))
        + 396*Mst2*pow2(Sbeta)) + 4*(-7*MuSUSY*Tbeta*(-1 + pow2(Sbeta)) + 54*
        Mst2*pow2(Sbeta))*pow4(Mst1) + 13*(-7*MuSUSY*Tbeta*(-1 + pow2(Sbeta)) +
        72*Mst2*pow2(Sbeta))*pow4(Mst2)) + 2*(-61*MuSUSY*Tbeta*(-1 + pow2(
        Sbeta)) + 537*Mst2*pow2(Sbeta))*pow6(Mst2)) - s2t*pow4(Msq)*(-3*pow2(
        Mst2)*(140*(-47953 + 14580*lmMsq + 25650*lmMst1 - 54162*lmMst2)*Tbeta*
        pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 3*Mst2*(7*(2489911 - 1425600*lmMsq -
        716040*lmMst1 + 2927880*lmMst2)*MuSUSY + 60*(5208*lmMst1 - 19*(13297 +
        6216*lmMst2))*Mst2*Tbeta)*pow2(Sbeta))*pow4(Mst1) - 9*pow2(Mst1)*(14*(-
        115931 + 48600*lmMsq + 53280*lmMst1 - 148320*lmMst2)*Tbeta*pow2(MuSUSY)
        *(-1 + pow2(Sbeta)) + 3*Mst2*(7*(1340537 - 496800*lmMsq + 4680*lmMst1 +
        778680*lmMst2)*MuSUSY + 12*(-53927 + 1820*lmMst1 - 230300*lmMst2)*Mst2*
        Tbeta)*pow2(Sbeta))*pow4(Mst2) - 14*(4*(-2286439 + 72900*lmMsq +
        307800*lmMst1 - 450360*lmMst2)*Tbeta*pow2(MuSUSY)*(-1 + pow2(Sbeta)) +
        15*Mst2*((1991719 - 194400*lmMsq - 328536*lmMst1 + 618840*lmMst2)*
        MuSUSY - 12*(78533 + 6732*lmMst1 + 180*lmMst2)*Mst2*Tbeta)*pow2(Sbeta))
        *pow6(Mst1) + 3265920*Tbeta*(pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 2*pow2(
        Mst2)*pow2(Sbeta))*pow6(Mst2)) + 2*Mt*pow4(Msq)*(9*pow2(Mst1)*(-7*(-
        1340537 + 496800*lmMsq - 4680*lmMst1 - 778680*lmMst2)*Tbeta*pow2(
        MuSUSY)*(-1 + pow2(Sbeta)) + 2*Mst2*(18*(74087 - 1820*lmMst1 + 230300*
        lmMst2)*MuSUSY + (740951 + 302400*lmMsq + 32760*lmMst1 - 955080*lmMst2
        - 45360*lmMt)*Mst2*Tbeta)*pow2(Sbeta))*pow3(Mst2) + 6*Mst2*(-7*(-
        3255761 + 1458000*lmMsq + 351000*lmMst1 - 2631960*lmMst2)*Tbeta*pow2(
        MuSUSY)*(-1 + pow2(Sbeta)) + Mst2*(72*(371369 - 7875*lmMst1 + 320355*
        lmMst2)*MuSUSY + (-21082949 + 7560*lmMst1 - 687960*lmMst2 + 680400*
        lmMt)*Mst2*Tbeta)*pow2(Sbeta))*pow4(Mst1) + pow2(Sbeta)*(56*(9*(520939
        + 44280*lmMst1 + 76680*lmMst2)*MuSUSY + (9227767 + 291600*lmMsq -
        754920*lmMst1 + 1088640*lmMst2 - 22680*lmMt)*Mst2*Tbeta)*pow6(Mst1) +
        6531840*MuSUSY*pow6(Mst2))) + 226800*Mst2*Mt*Tbeta*pow2(Mst1)*(-(pow2(
        Mst2)*pow2(MuSUSY)*(-1 + pow2(Sbeta))*(313*pow2(Mst1)*pow2(Mst2) + 93*
        pow4(Mst1) + 179*pow4(Mst2))) + 8*pow2(Msq)*(-3*pow2(MuSUSY)*(-1 +
        pow2(Sbeta))*(24*pow2(Mst1)*pow2(Mst2) + 10*pow4(Mst1) + 13*pow4(Mst2))
        + 26*pow2(Sbeta)*pow6(Mst2)) + 232*pow2(Sbeta)*pow8(Mst2))))/pow2(
        Sbeta)))/(pow2(Mst1)*pow4(Msq))))/(Tbeta*pow4(Mst2))))/204120.))/pow2(
        Mgl) + xDR2DRMOD*((threeLoopFlag*pow2(Al4p)*(-12*pow2(Mt)*pow2(s2t)*((
        8*(2*Dmglst2*(889 + 694*lmMst2 - 576*lmMst1*(1 + lmMst2) + 90*lmMsq*(-1
        + 6*lmMst2) - 270*pow2(lmMsq) - 48*pow2(lmMst2)) - 3*Mgl*(207 + 494*
        lmMst2 - 30*lmMsq*(5 + 6*lmMst2) + 90*pow2(lmMsq) + 240*pow2(lmMst2)))*
        pow2(Mst1) + (1536*((lmMst1 - lmMst2)*(1 + lmMst2)*Mgl + 2*Dmglst2*(2 +
        7*lmMst2 - lmMst1*(5 + 4*lmMst2) + 4*pow2(lmMst2)))*pow4(Mst1))/pow2(
        Mst2) + (10*(18*Dmglst2*(13 - 28*lmMsq + 28*lmMst2) + (67 - 84*lmMsq +
        84*lmMst2)*Mgl)*(pow2(Mst1) - 2*pow2(Mst2))*pow4(Mst2))/pow4(Msq) + (8*
        pow2(Mst2)*((160*Dmglst2*(-14 + 15*lmMsq - 15*lmMst2) + 20*(-43 + 30*
        lmMsq - 30*lmMst2)*Mgl - (3*(Mgl*(303 + 590*lmMst2 + 64*lmMst1*(1 +
        lmMst2) - 30*lmMsq*(5 + 6*lmMst2) + 90*pow2(lmMsq) + 240*pow2(lmMst2))
        + 2*Dmglst2*(-179 + (238 + 64*lmMst1)*lmMst2 - 30*lmMsq*(-1 + 6*lmMst2)
        + 90*pow2(lmMsq) + 272*pow2(lmMst2)))*pow2(Msq))/pow2(Mst1))*pow2(Mst2)
        - 2*((40*Dmglst2*(-14 + 15*lmMsq - 15*lmMst2) + 5*(-43 + 30*lmMsq - 30*
        lmMst2)*Mgl)*pow2(Mst1) + pow2(Msq)*(2*Dmglst2*(329 + 90*lmMsq*(-1 + 6*
        lmMst2) + 96*lmMst1*(1 - 2*lmMst2*(1 + lmMst2)) - 270*pow2(lmMsq) +
        lmMst2*(-586 + 96*pow2(lmMst1)) - 528*pow2(lmMst2) + 96*pow3(lmMst2)) -
        3*Mgl*(207 + 414*lmMst2 - 30*lmMsq*(5 + 6*lmMst2) + 90*pow2(lmMsq) +
        128*pow2(lmMst2) + 32*lmMst1*(3 + 5*lmMst2 + 2*pow2(lmMst2)) - 32*((1 +
        lmMst2)*pow2(lmMst1) + pow3(lmMst2))))) + (2*(Mgl*(48*pow2(1 + lmMst2)*
        pow2(Msq) + 5*(43 - 30*lmMsq + 30*lmMst2)*pow2(Mst1)) + 8*Dmglst2*(24*
        pow2(lmMst2)*pow2(Msq) + 5*(14 - 15*lmMsq)*pow2(Mst1) + 3*lmMst2*(8*
        pow2(Msq) + 25*pow2(Mst1))))*pow4(Mst2))/pow4(Mst1)))/pow2(Msq))/Mgl -
        pow2(MuSUSY)*(3636 - 1800*lmMsq + 1080*pow2(lmMsq) - 24*lmMst1*(431 -
        150*lmMsq + 90*pow2(lmMsq)) + 768*pow2(lmMst1) + 48*lmMst2*(363 + 30*
        lmMsq*(-4 + 3*lmMst1 - 3*lmMst2) + 451*lmMst2 - lmMst1*(407 + 168*
        lmMst2) + 45*pow2(lmMsq) + 16*pow2(lmMst1)) + 7296*pow3(lmMst2) + ((12*
        (Mgl*(335 + 654*lmMst2 + 64*lmMst1*(1 + lmMst2) - 30*lmMsq*(5 + 6*
        lmMst2) + 90*pow2(lmMsq) + 272*pow2(lmMst2)) + 2*Dmglst2*(-115 + (366 +
        64*lmMst1)*lmMst2 - 30*lmMsq*(-1 + 6*lmMst2) + 90*pow2(lmMsq) + 336*
        pow2(lmMst2)))*pow2(Mst2))/pow2(Mst1) + 24*Dmglst2*(45 + 30*lmMsq*(1 -
        6*lmMst2)*(1 - 2*lmMst1 + 2*lmMst2) + 90*(1 - 2*lmMst1 + 2*lmMst2)*
        pow2(lmMsq) + 8*lmMst2*(29 + 8*pow2(lmMst1)) + 1068*pow2(lmMst2) - 2*
        lmMst1*(-83 + 430*lmMst2 + 336*pow2(lmMst2)) + 608*pow3(lmMst2)) +
        pow2(Mst1)*((80*(lmMst1 - lmMst2)*(8*Dmglst2*(14 - 15*lmMsq + 15*
        lmMst2) + (43 - 30*lmMsq + 30*lmMst2)*Mgl))/pow2(Msq) + (24*(2*Dmglst2*
        (96 + lmMst2*(173 + 30*lmMsq + 90*pow2(lmMsq) + 288*pow2(lmMst1)) + (
        526 - 180*lmMsq)*pow2(lmMst2) - lmMst1*(13 + 30*lmMsq*(1 - 6*lmMst2) +
        526*lmMst2 + 90*pow2(lmMsq) + 848*pow2(lmMst2)) + 560*pow3(lmMst2)) +
        Mgl*(64 + 15*lmMst2*(33 - 10*lmMsq + 6*pow2(lmMsq)) + 288*(1 + lmMst2)*
        pow2(lmMst1) + 6*(173 - 30*lmMsq)*pow2(lmMst2) - lmMst1*(431 + 1326*
        lmMst2 - 30*lmMsq*(5 + 6*lmMst2) + 90*pow2(lmMsq) + 848*pow2(lmMst2)) +
        560*pow3(lmMst2))))/pow2(Mst2)) + ((-5*((64*Dmglst2*(14 - 15*lmMsq +
        15*lmMst2) + 8*(43 - 30*lmMsq + 30*lmMst2)*Mgl)*pow2(Msq) + (18*
        Dmglst2*(13 - 28*lmMsq + 28*lmMst2) + (67 - 84*lmMsq + 84*lmMst2)*Mgl)*
        pow2(Mst2)))/(pow2(Mst1)*pow4(Msq)) - (384*(1 + lmMst2)*(4*Dmglst2*
        lmMst2 + Mgl + lmMst2*Mgl))/pow4(Mst1))*pow4(Mst2) - (40*(8*Dmglst2*(14
        - 15*lmMsq + 15*lmMst2) + (43 - 30*lmMsq + 30*lmMst2)*Mgl)*(pow4(Mst2)
        - 2*(lmMst1 - lmMst2)*(pow4(Mst1) + pow4(Mst2))))/(pow2(Msq)*pow2(Mst2)
        ) + (-5*(18*Dmglst2*(13 - 28*lmMsq + 28*lmMst2) + (67 - 84*lmMsq + 84*
        lmMst2)*Mgl)*pow2(Mst2)*(pow2(Mst2) - 2*(lmMst1 - lmMst2)*(pow2(Mst1) +
        pow2(Mst2))) + (2*pow4(Mst1)*((24*Dmglst2*(96 + lmMst2*(285 + 30*lmMsq
        + 90*pow2(lmMsq) + 544*pow2(lmMst1)) + (590 - 180*lmMsq)*pow2(lmMst2) -
        lmMst1*(109 + 30*lmMsq*(1 - 6*lmMst2) + 590*lmMst2 + 90*pow2(lmMsq) +
        1360*pow2(lmMst2)) + 816*pow3(lmMst2)) + 12*Mgl*(80 + lmMst2*(479 -
        150*lmMsq + 90*pow2(lmMsq)) + 544*(1 + lmMst2)*pow2(lmMst1) + 2*(631 -
        90*lmMsq)*pow2(lmMst2) - lmMst1*(399 + 1806*lmMst2 - 30*lmMsq*(5 + 6*
        lmMst2) + 90*pow2(lmMsq) + 1360*pow2(lmMst2)) + 816*pow3(lmMst2)))*
        pow4(Msq) + (lmMst1 - lmMst2)*(90*Dmglst2*(13 - 28*lmMsq + 28*lmMst2) +
        5*(67 - 84*lmMsq + 84*lmMst2)*Mgl)*pow4(Mst2)))/pow4(Mst2))/pow4(Msq))/
        Mgl)) + (2592*s2t*pow2(Mt)*pow2(MuSUSY)*(-(s2t*(3636 - 1800*lmMsq +
        1080*pow2(lmMsq) - 24*lmMst1*(431 - 150*lmMsq + 90*pow2(lmMsq)) + 768*
        pow2(lmMst1) + 48*lmMst2*(363 + 30*lmMsq*(-4 + 3*lmMst1 - 3*lmMst2) +
        451*lmMst2 - lmMst1*(407 + 168*lmMst2) + 45*pow2(lmMsq) + 16*pow2(
        lmMst1)) + 7296*pow3(lmMst2) + ((12*(Mgl*(335 + 654*lmMst2 + 64*lmMst1*
        (1 + lmMst2) - 30*lmMsq*(5 + 6*lmMst2) + 90*pow2(lmMsq) + 272*pow2(
        lmMst2)) + 2*Dmglst2*(-115 + (366 + 64*lmMst1)*lmMst2 - 30*lmMsq*(-1 +
        6*lmMst2) + 90*pow2(lmMsq) + 336*pow2(lmMst2)))*pow2(Mst2))/pow2(Mst1)
        + 24*Dmglst2*(45 + 30*lmMsq*(1 - 6*lmMst2)*(1 - 2*lmMst1 + 2*lmMst2) +
        90*(1 - 2*lmMst1 + 2*lmMst2)*pow2(lmMsq) + 8*lmMst2*(29 + 8*pow2(
        lmMst1)) + 1068*pow2(lmMst2) - 2*lmMst1*(-83 + 430*lmMst2 + 336*pow2(
        lmMst2)) + 608*pow3(lmMst2)) + pow2(Mst1)*((80*(lmMst1 - lmMst2)*(8*
        Dmglst2*(14 - 15*lmMsq + 15*lmMst2) + (43 - 30*lmMsq + 30*lmMst2)*Mgl))
        /pow2(Msq) + (24*(2*Dmglst2*(96 + lmMst2*(173 + 30*lmMsq + 90*pow2(
        lmMsq) + 288*pow2(lmMst1)) + (526 - 180*lmMsq)*pow2(lmMst2) - lmMst1*(
        13 + 30*lmMsq*(1 - 6*lmMst2) + 526*lmMst2 + 90*pow2(lmMsq) + 848*pow2(
        lmMst2)) + 560*pow3(lmMst2)) + Mgl*(64 + 15*lmMst2*(33 - 10*lmMsq + 6*
        pow2(lmMsq)) + 288*(1 + lmMst2)*pow2(lmMst1) + 6*(173 - 30*lmMsq)*pow2(
        lmMst2) - lmMst1*(431 + 1326*lmMst2 - 30*lmMsq*(5 + 6*lmMst2) + 90*
        pow2(lmMsq) + 848*pow2(lmMst2)) + 560*pow3(lmMst2))))/pow2(Mst2)) + ((-
        5*((64*Dmglst2*(14 - 15*lmMsq + 15*lmMst2) + 8*(43 - 30*lmMsq + 30*
        lmMst2)*Mgl)*pow2(Msq) + (18*Dmglst2*(13 - 28*lmMsq + 28*lmMst2) + (67
        - 84*lmMsq + 84*lmMst2)*Mgl)*pow2(Mst2)))/(pow2(Mst1)*pow4(Msq)) - (
        384*(1 + lmMst2)*(4*Dmglst2*lmMst2 + Mgl + lmMst2*Mgl))/pow4(Mst1))*
        pow4(Mst2) - (40*(8*Dmglst2*(14 - 15*lmMsq + 15*lmMst2) + (43 - 30*
        lmMsq + 30*lmMst2)*Mgl)*(pow4(Mst2) - 2*(lmMst1 - lmMst2)*(pow4(Mst1) +
        pow4(Mst2))))/(pow2(Msq)*pow2(Mst2)) + (-5*(18*Dmglst2*(13 - 28*lmMsq +
        28*lmMst2) + (67 - 84*lmMsq + 84*lmMst2)*Mgl)*pow2(Mst2)*(pow2(Mst2) -
        2*(lmMst1 - lmMst2)*(pow2(Mst1) + pow2(Mst2))) + (2*pow4(Mst1)*((24*
        Dmglst2*(96 + lmMst2*(285 + 30*lmMsq + 90*pow2(lmMsq) + 544*pow2(
        lmMst1)) + (590 - 180*lmMsq)*pow2(lmMst2) - lmMst1*(109 + 30*lmMsq*(1 -
        6*lmMst2) + 590*lmMst2 + 90*pow2(lmMsq) + 1360*pow2(lmMst2)) + 816*
        pow3(lmMst2)) + 12*Mgl*(80 + lmMst2*(479 - 150*lmMsq + 90*pow2(lmMsq))
        + 544*(1 + lmMst2)*pow2(lmMst1) + 2*(631 - 90*lmMsq)*pow2(lmMst2) -
        lmMst1*(399 + 1806*lmMst2 - 30*lmMsq*(5 + 6*lmMst2) + 90*pow2(lmMsq) +
        1360*pow2(lmMst2)) + 816*pow3(lmMst2)))*pow4(Msq) + (lmMst1 - lmMst2)*(
        90*Dmglst2*(13 - 28*lmMsq + 28*lmMst2) + 5*(67 - 84*lmMsq + 84*lmMst2)*
        Mgl)*pow4(Mst2)))/pow4(Mst2))/pow4(Msq))/Mgl))/216. + (64*Mt*(2*
        Dmglst2*pow2(Mst2)*(pow2(Mst1)*pow2(Mst2)*(8 + 7*lmMst2 - 11*pow2(
        lmMst2) + lmMst1*(-3 + 5*lmMst2 + 6*pow2(lmMst2)) - 6*pow3(lmMst2)) +
        2*(8 + 13*lmMst2 - 8*pow2(lmMst2) + lmMst1*(-5 + 5*lmMst2 + 6*pow2(
        lmMst2)) - 6*pow3(lmMst2))*pow4(Mst1) + (1 - 2*lmMst2 - 3*pow2(lmMst2))
        *pow4(Mst2)) + (1 + lmMst2)*Mgl*(2*(1 - 10*lmMst2 + 4*lmMst1*(2 +
        lmMst2) - 4*pow2(lmMst2))*pow2(Mst2)*pow4(Mst1) - 2*(1 + 5*lmMst2 -
        lmMst1*(3 + 2*lmMst2) + 2*pow2(lmMst2))*pow2(Mst1)*pow4(Mst2) + (7 -
        32*lmMst2 + 4*lmMst1*(7 + 3*lmMst2) - 12*pow2(lmMst2))*pow6(Mst1) - 2*(
        1 + lmMst2)*pow6(Mst2))))/(9.*Mgl*pow2(Mst1)*pow5(Mst2))))/pow2(Sbeta)
        + (3*pow3(s2t)*((3072*Mt*((1 + lmMst2)*Mgl*(2*(2 + 2*lmMst1 - lmMst2)*
        pow2(Mst2)*pow4(Mst1) + 2*(1 - 3*lmMst2 + lmMst1*(3 + 2*lmMst2) - 2*
        pow2(lmMst2))*pow2(Mst1)*pow4(Mst2) + (1 + 2*lmMst1 - 2*lmMst2)*pow6(
        Mst1) - 2*(1 + lmMst2)*pow6(Mst2)) + Dmglst2*(2*(1 - 4*lmMst1 + 10*
        lmMst2 + 3*pow2(lmMst2))*pow2(Mst2)*pow4(Mst1) + 2*pow2(Mst1)*(6 + 11*
        lmMst2 - 5*pow2(lmMst2) + lmMst1*(-3 + 5*lmMst2 + 6*pow2(lmMst2)) - 6*
        pow3(lmMst2))*pow4(Mst2) + (1 + 13*lmMst2 - 2*lmMst1*(5 + 3*lmMst2) +
        6*pow2(lmMst2))*pow6(Mst1) - 2*(-1 + 2*lmMst2 + 3*pow2(lmMst2))*pow6(
        Mst2))))/(Mst2*pow2(Mst1)) + s2t*((4*pow2(Mst1)*pow2(Mst2)*(80*Dmglst2*
        (14 - 15*lmMsq + 15*lmMst2)*pow2(Mst1) + 10*(43 - 30*lmMsq + 30*lmMst2)
        *Mgl*pow2(Mst1) + 3*Mgl*pow2(Msq)*(143 + 1260*lmMst2*(1 + lmMst2) + 30*
        lmMsq*(-1 + 2*lmMst1 - 2*lmMst2)*(5 + 6*lmMst2) + 90*(1 - 2*lmMst1 + 2*
        lmMst2)*pow2(lmMsq) - 448*(1 + lmMst2)*pow2(lmMst1) + 2*lmMst1*(-463 -
        334*lmMst2 + 176*pow2(lmMst2)) + 96*pow3(lmMst2)) + 6*Dmglst2*pow2(Msq)
        *(13 + 30*lmMsq*(1 - 6*lmMst2)*(1 - 2*lmMst1 + 2*lmMst2) + 90*(1 - 2*
        lmMst1 + 2*lmMst2)*pow2(lmMsq) - 8*lmMst2*(31 + 56*pow2(lmMst1)) + 748*
        pow2(lmMst2) + lmMst1*(358 - 732*lmMst2 + 352*pow2(lmMst2)) + 96*pow3(
        lmMst2))))/pow2(Msq) + 12*(-2*Dmglst2*(-147 + 2*(55 + 32*lmMst1)*lmMst2
        - 30*lmMsq*(-1 + 6*lmMst2) + 90*pow2(lmMsq) + 144*pow2(lmMst2)) - Mgl*(
        207 + 430*lmMst2 + 64*lmMst1*(1 + lmMst2) - 30*lmMsq*(5 + 6*lmMst2) +
        90*pow2(lmMsq) + 176*pow2(lmMst2)))*pow4(Mst1) + (-((40*(1 - 2*lmMst1 +
        2*lmMst2)*(8*Dmglst2*(14 - 15*lmMsq + 15*lmMst2) + (43 - 30*lmMsq + 30*
        lmMst2)*Mgl)*pow2(Msq)*pow2(Mst1) + (-12*Mgl*(399 - 80*lmMst2 - 30*
        lmMsq*(1 + 2*lmMst1 - 2*lmMst2)*(5 + 6*lmMst2) + 90*(1 + 2*lmMst1 - 2*
        lmMst2)*pow2(lmMsq) - 64*(1 + lmMst2)*pow2(lmMst1) - 1228*pow2(lmMst2)
        + 2*lmMst1*(495 + 878*lmMst2 + 336*pow2(lmMst2)) - 608*pow3(lmMst2)) +
        24*Dmglst2*(275 - 564*lmMst2 + 30*lmMsq*(1 + 2*lmMst1 - 2*lmMst2)*(-1 +
        6*lmMst2) - 90*(1 + 2*lmMst1 - 2*lmMst2)*pow2(lmMsq) + 64*lmMst2*pow2(
        lmMst1) + 332*pow2(lmMst2) - 2*lmMst1*(-83 + 494*lmMst2 + 336*pow2(
        lmMst2)) + 608*pow3(lmMst2)))*pow4(Msq) + (90*Dmglst2*(-13 + 28*lmMsq -
        28*lmMst2) + 5*(-67 + 84*lmMsq - 84*lmMst2)*Mgl)*pow4(Mst1))*pow4(Mst2)
        ) + ((40*(1 + 2*lmMst1 - 2*lmMst2)*(8*Dmglst2*(-14 + 15*lmMsq - 15*
        lmMst2) + (-43 + 30*lmMsq - 30*lmMst2)*Mgl)*pow2(Msq)*pow2(Mst1) + (-
        12*Mgl*(399 + 782*lmMst2 + 64*lmMst1*(1 + lmMst2) - 30*lmMsq*(5 + 6*
        lmMst2) + 90*pow2(lmMsq) + 336*pow2(lmMst2)) - 24*Dmglst2*(-115 + (494
        + 64*lmMst1)*lmMst2 - 30*lmMsq*(-1 + 6*lmMst2) + 90*pow2(lmMsq) + 464*
        pow2(lmMst2)))*pow4(Msq) + (1 - 2*lmMst1 + 2*lmMst2)*(90*Dmglst2*(-13 +
        28*lmMsq - 28*lmMst2) + 5*(-67 + 84*lmMsq - 84*lmMst2)*Mgl)*pow4(Mst1))
        *pow6(Mst2))/pow2(Mst1))/pow4(Msq))))/Mgl + (16*pow3(Mt)*(128*s2t*pow2(
        Mst1)*(-(Dmglst2*pow2(Mst2)*(pow2(Mst1)*pow2(Mst2)*(pow2(Mst2)*(253 +
        5*lmMst2*(107 - 6*lmMt) - 102*lmMt - 18*lmMst1*(9 + lmMst2 - 2*lmMt +
        2*lmMst2*lmMt - 2*pow2(lmMst2)) + 12*(10 + 3*lmMt)*pow2(lmMst2) - 36*
        pow3(lmMst2)) + 18*pow2(MuSUSY)*(8 + 7*lmMst2 - 11*pow2(lmMst2) +
        lmMst1*(-3 + 5*lmMst2 + 6*pow2(lmMst2)) - 6*pow3(lmMst2))) + 18*((2*
        pow2(MuSUSY)*(8 + 13*lmMst2 - 8*pow2(lmMst2) + lmMst1*(-5 + 5*lmMst2 +
        6*pow2(lmMst2)) - 6*pow3(lmMst2)) + pow2(Mst2)*(23 + lmMst2*(93 - 22*
        lmMt) - 14*lmMt - 4*(-11 + lmMt)*pow2(lmMst2) - 2*lmMst1*(25 + lmMst2*(
        17 - 2*lmMt) - 6*lmMt + 2*pow2(lmMst2)) + 4*pow3(lmMst2)))*pow4(Mst1) +
        (1 - 2*lmMst2 - 3*pow2(lmMst2))*(2*pow2(Mst2) + pow2(MuSUSY))*pow4(
        Mst2) + (14 + lmMst2*(189 - 48*lmMt) - 22*lmMt + lmMst1*(-137 + 30*lmMt
        + 9*lmMst2*(-13 + 2*lmMt) - 18*pow2(lmMst2)) - 9*(-15 + 2*lmMt)*pow2(
        lmMst2) + 18*pow3(lmMst2))*pow6(Mst1)))) + 9*(1 + lmMst2)*Mgl*(-2*pow2(
        Mst2)*((3 + 2*lmMst1*(7 + 2*lmMst2 - 2*lmMt) + 4*lmMst2*(-4 + lmMt) +
        2*lmMt - 4*pow2(lmMst2))*pow2(Mst2) + (1 - 10*lmMst2 + 4*lmMst1*(2 +
        lmMst2) - 4*pow2(lmMst2))*pow2(MuSUSY))*pow4(Mst1) + pow2(Mst1)*((1 -
        2*lmMst1*(5 + 2*lmMst2 - 2*lmMt) - 4*lmMst2*(-3 + lmMt) - 6*lmMt + 4*
        pow2(lmMst2))*pow2(Mst2) + 2*(1 + 5*lmMst2 - lmMst1*(3 + 2*lmMst2) + 2*
        pow2(lmMst2))*pow2(MuSUSY))*pow4(Mst2) - (2*(8 - 27*lmMst2 + lmMst1*(25
        + 6*lmMst2 - 6*lmMt) + 2*lmMt + 6*lmMst2*lmMt - 6*pow2(lmMst2))*pow2(
        Mst2) + (7 - 32*lmMst2 + 4*lmMst1*(7 + 3*lmMst2) - 12*pow2(lmMst2))*
        pow2(MuSUSY))*pow6(Mst1) + 2*(1 + lmMst2)*(2*pow2(Mst2) + pow2(MuSUSY))
        *pow6(Mst2))) + Mst2*Mt*(8*Dmglst2*((9115 + 12498*lmMst2 + 270*lmMsq*(-
        1 + 6*lmMst2) - 1152*lmMst1*(4 + 2*lmMst2 - lmMt) - 96*(29 + 27*lmMst2)
        *lmMt - 810*pow2(lmMsq) + 1872*pow2(lmMst2))*pow4(Mst1)*pow4(Mst2) - 9*
        (-243 + 2*(55 + 32*lmMst1)*lmMst2 - 30*lmMsq*(-1 + 6*lmMst2) + 90*pow2(
        lmMsq) + 208*pow2(lmMst2))*pow2(Mst1)*pow6(Mst2) + 576*(pow2(Mst2)*(31
        + lmMst2*(95 - 23*lmMt) - 16*lmMt + (51 - 6*lmMt)*pow2(lmMst2) - 2*
        lmMst1*(25 + lmMst2*(20 - 3*lmMt) - 6*lmMt + 3*pow2(lmMst2)) + 6*pow3(
        lmMst2))*pow6(Mst1) + (42 + lmMst2*(242 - 62*lmMt) - 33*lmMt + 6*(29 -
        4*lmMt)*pow2(lmMst2) - 2*lmMst1*(81 + lmMst2*(74 - 12*lmMt) - 18*lmMt +
        12*pow2(lmMst2)) + 24*pow3(lmMst2))*pow8(Mst1)) + 576*lmMst2*(1 +
        lmMst2)*pow8(Mst2)) + (pow2(Mst1)*(pow2(Mst1) + pow2(Mst2))*(120*(8*
        Dmglst2*(14 - 15*lmMsq + 15*lmMst2) + (43 - 30*lmMsq + 30*lmMst2)*Mgl)*
        pow2(Msq)*pow6(Mst2) + (270*Dmglst2*(13 - 28*lmMsq + 28*lmMst2) + 15*(
        67 - 84*lmMsq + 84*lmMst2)*Mgl)*pow8(Mst2)))/pow4(Msq) + 4*Mgl*(-((3559
        + 10942*lmMst2 - 270*lmMsq*(5 + 6*lmMst2) - 1152*lmMst1*(1 + lmMst2)*(3
        + lmMst2 - lmMt) + 810*pow2(lmMsq) + 8112*pow2(lmMst2) - 192*lmMt*(7 +
        13*lmMst2 + 6*pow2(lmMst2)) + 1152*pow3(lmMst2))*pow4(Mst1)*pow4(Mst2))
        + 576*(1 + lmMst2)*((5 - 57*lmMst2 + 2*lmMst1*(25 + 6*lmMst2 - 6*lmMt)
        + 7*lmMt + 12*lmMst2*lmMt - 12*pow2(lmMst2))*pow2(Mst1) + (-2 - 27*
        lmMst2 + lmMst1*(22 + 6*lmMst2 - 6*lmMt) + 5*lmMt + 6*lmMst2*lmMt - 6*
        pow2(lmMst2))*pow2(Mst2))*pow6(Mst1) - 9*(271 + 526*lmMst2 + 64*lmMst1*
        (1 + lmMst2) - 30*lmMsq*(5 + 6*lmMst2) + 90*pow2(lmMsq) + 208*pow2(
        lmMst2))*pow2(Mst1)*pow6(Mst2) + 288*pow2(1 + lmMst2)*pow8(Mst2)))))/(
        Mgl*pow4(Mst1)*pow5(Mst2)) + (4*Mt*MuSUSY*(1024*pow2(Mst1)*pow3(Mt)*(
        Dmglst2*pow2(Mst2)*(pow2(Mst1)*pow2(Mst2)*(119 + lmMst2*(224 - 15*lmMt)
        - 51*lmMt + 3*(5 + 6*lmMt)*pow2(lmMst2) + 18*lmMst1*(-4 + lmMt -
        lmMst2*(1 + lmMt) + pow2(lmMst2)) - 18*pow3(lmMst2)) + (263 + lmMst2*(
        962 - 213*lmMt) - 177*lmMt + (393 - 18*lmMt)*pow2(lmMst2) - 18*lmMst1*(
        26 - lmMst2*(-17 + lmMt) - 7*lmMt + pow2(lmMst2)) + 18*pow3(lmMst2))*
        pow4(Mst1) - 18*(-1 + 2*lmMst2 + 3*pow2(lmMst2))*pow4(Mst2)) + 9*(1 +
        lmMst2)*Mgl*((3 - 21*lmMst2 + 2*lmMst1*(8 + 3*lmMst2 - 3*lmMt) + 5*lmMt
        + 6*lmMst2*lmMt - 6*pow2(lmMst2))*pow2(Mst2)*pow4(Mst1) + (-1 - 7*
        lmMst2 + 2*lmMst1*(2 + lmMst2 - lmMt) + 3*lmMt + 2*lmMst2*lmMt - 2*
        pow2(lmMst2))*pow2(Mst1)*pow4(Mst2) + (12 - 45*lmMst2 + 2*lmMst1*(19 +
        6*lmMst2 - 6*lmMt) + 7*lmMt + 12*lmMst2*lmMt - 12*pow2(lmMst2))*pow6(
        Mst1) - 2*(1 + lmMst2)*pow6(Mst2))) - 6912*Mt*pow2(Mst1)*pow2(Mst2)*
        pow2(s2t)*((1 + lmMst2)*Mgl*(2*(2 - 5*lmMst2 + lmMst1*(5 + 2*lmMst2) -
        2*pow2(lmMst2))*pow2(Mst2)*pow4(Mst1) + 2*(-2*lmMst2*(2 + lmMst2) +
        lmMst1*(3 + 2*lmMst2))*pow2(Mst1)*pow4(Mst2) + (5 - 12*lmMst2 + 4*
        lmMst1*(3 + lmMst2) - 4*pow2(lmMst2))*pow6(Mst1) - 2*(1 + lmMst2)*pow6(
        Mst2)) + Dmglst2*(2*pow2(Mst2)*(8 + 19*lmMst2 - 5*pow2(lmMst2) +
        lmMst1*(-7 + 5*lmMst2 + 6*pow2(lmMst2)) - 6*pow3(lmMst2))*pow4(Mst1) +
        2*pow2(Mst1)*(7 + 9*lmMst2 - 8*pow2(lmMst2) + lmMst1*(-3 + 5*lmMst2 +
        6*pow2(lmMst2)) - 6*pow3(lmMst2))*pow4(Mst2) + (17 + 51*lmMst2 - 4*
        pow2(lmMst2) + 4*lmMst1*(-6 + lmMst2 + 3*pow2(lmMst2)) - 12*pow3(
        lmMst2))*pow6(Mst1) - 2*(-1 + 2*lmMst2 + 3*pow2(lmMst2))*pow6(Mst2))) -
        (12*Mst2*s2t*pow2(Mt)*(pow2(Mst1)*(pow2(Mst1) - pow2(Mst2))*(40*(8*
        Dmglst2*(14 - 15*lmMsq + 15*lmMst2) + (43 - 30*lmMsq + 30*lmMst2)*Mgl)*
        pow2(Msq)*pow6(Mst2) + (90*Dmglst2*(13 - 28*lmMsq + 28*lmMst2) + 5*(67
        - 84*lmMsq + 84*lmMst2)*Mgl)*pow8(Mst2)) + pow4(Msq)*(-8*Dmglst2*(-((
        121 - 650*lmMst2 + 90*lmMsq*(-1 + 6*lmMst2) - 270*pow2(lmMsq) + 192*
        lmMst2*pow2(lmMst1) - 432*pow2(lmMst2) - 192*lmMst1*(-1 + lmMst2 + 2*
        pow2(lmMst2)) + 192*pow3(lmMst2))*pow4(Mst1)*pow4(Mst2)) - 192*(lmMst2*
        pow2(lmMst1) + lmMst1*(4 + 2*lmMst2 - 2*pow2(lmMst2)) + (-4 + lmMst2)*
        pow2(1 + lmMst2))*pow2(Mst2)*pow6(Mst1) - 3*(-179 + 2*(87 + 32*lmMst1)*
        lmMst2 - 30*lmMsq*(-1 + 6*lmMst2) + 90*pow2(lmMsq) + 208*pow2(lmMst2))*
        pow2(Mst1)*pow6(Mst2) - 192*(-6 - 14*lmMst2 + lmMst2*pow2(lmMst1) +
        lmMst1*(9 + 6*lmMst2 - 2*pow2(lmMst2)) - 6*pow2(lmMst2) + pow3(lmMst2))
        *pow8(Mst1) + 192*lmMst2*(1 + lmMst2)*pow8(Mst2)) + 12*Mgl*(-((143 +
        302*lmMst2 - 30*lmMsq*(5 + 6*lmMst2) + 90*pow2(lmMsq) + 48*pow2(lmMst2)
        + 128*lmMst1*pow2(1 + lmMst2) - 64*((1 + lmMst2)*pow2(lmMst1) + pow3(
        lmMst2)))*pow4(Mst1)*pow4(Mst2)) + (271 + 526*lmMst2 + 64*lmMst1*(1 +
        lmMst2) - 30*lmMsq*(5 + 6*lmMst2) + 90*pow2(lmMsq) + 208*pow2(lmMst2))*
        pow2(Mst1)*pow6(Mst2) + 64*(1 + lmMst2)*(pow2(1 - lmMst1 + lmMst2)*
        pow2(Mst2)*pow6(Mst1) + (1 + 3*lmMst2 - lmMst1*(3 + 2*lmMst2) + pow2(
        lmMst1) + pow2(lmMst2))*pow8(Mst1)) - 32*pow2(1 + lmMst2)*pow8(Mst2))))
        )/pow4(Msq) + (3*pow3(Mst2)*pow3(s2t)*(Mgl*(40*(43 - 30*lmMsq + 30*
        lmMst2)*pow2(Msq)*pow2(Mst1)*(2*(lmMst1 - lmMst2)*pow2(Mst1)*pow2(Mst2)
        + pow4(Mst1) - pow4(Mst2))*pow4(Mst2) + 5*(67 - 84*lmMsq + 84*lmMst2)*(
        pow2(Mst1) + 2*(lmMst1 - lmMst2)*pow2(Mst2))*pow4(Mst1)*pow6(Mst2) +
        12*pow4(Msq)*(-2*(16 - 3*lmMst2*(133 - 50*lmMsq + 30*pow2(lmMsq)) - 32*
        (1 + lmMst2)*pow2(lmMst1) + 2*(-383 + 90*lmMsq)*pow2(lmMst2) + lmMst1*(
        463 + 846*lmMst2 - 30*lmMsq*(5 + 6*lmMst2) + 90*pow2(lmMsq) + 336*pow2(
        lmMst2)) - 304*pow3(lmMst2))*pow4(Mst1)*pow4(Mst2) - pow2(Mst2)*(175 +
        462*lmMst2 + 1024*lmMst1*lmMst2*(1 + lmMst2) - 30*lmMsq*(5 + 6*lmMst2)
        + 90*pow2(lmMsq) - 272*pow2(lmMst2) - 512*((1 + lmMst2)*pow2(lmMst1) +
        pow3(lmMst2)))*pow6(Mst1) + (367 + 718*lmMst2 + 64*lmMst1*(1 + lmMst2)
        - 30*lmMsq*(5 + 6*lmMst2) + 90*pow2(lmMsq) + 304*pow2(lmMst2))*pow2(
        Mst1)*pow6(Mst2) + 32*(1 + lmMst2)*(1 + lmMst1*(2 - 32*lmMst2) - 2*
        lmMst2 + 16*pow2(lmMst1) + 16*pow2(lmMst2))*pow8(Mst1) - 32*pow2(1 +
        lmMst2)*pow8(Mst2))) + 2*Dmglst2*(384*lmMst2*pow4(Msq)*pow4(Mst1)*(2*
        pow2(lmMst1)*(8*pow2(Mst1)*pow2(Mst2) + 8*pow4(Mst1) + pow4(Mst2)) +
        pow2(lmMst2)*(16*pow2(Mst1)*(pow2(Mst1) + pow2(Mst2)) + 19*pow4(Mst2)))
        + pow2(Mst1)*pow2(Mst2)*(45*(13 - 28*lmMsq)*pow4(Mst1)*pow4(Mst2) +
        160*(-14 + 15*lmMsq)*pow2(Msq)*pow2(Mst2)*(-pow4(Mst1) + pow4(Mst2)) +
        12*pow4(Msq)*(160*pow2(Mst1)*pow2(Mst2) - 3*(-49 + 10*lmMsq + 30*pow2(
        lmMsq))*pow4(Mst1) + 5*(-23 + 6*lmMsq + 18*pow2(lmMsq))*pow4(Mst2))) -
        24*pow2(lmMst2)*(200*pow2(Msq)*pow4(Mst1)*pow6(Mst2) + 105*pow4(Mst1)*
        pow8(Mst2) + pow4(Msq)*(6*(-61 + 30*lmMsq)*pow4(Mst1)*pow4(Mst2) + 8*
        pow2(Mst2)*pow6(Mst1) - 200*pow2(Mst1)*pow6(Mst2) - 64*pow8(Mst1) + 32*
        pow8(Mst2))) + 2*lmMst1*(-12*pow4(Msq)*((-83 + 462*lmMst2 - 30*lmMsq*(-
        1 + 6*lmMst2) + 90*pow2(lmMsq) + 336*pow2(lmMst2))*pow4(Mst1)*pow4(
        Mst2) + 32*(3 + 3*lmMst2 + 16*pow2(lmMst2))*pow2(Mst2)*pow6(Mst1) - 32*
        lmMst2*pow2(Mst1)*pow6(Mst2) + 32*(3 + 2*lmMst2 + 16*pow2(lmMst2))*
        pow8(Mst1)) + pow4(Mst1)*(160*(14 - 15*lmMsq + 15*lmMst2)*pow2(Msq)*
        pow6(Mst2) + 45*(13 - 28*lmMsq + 28*lmMst2)*pow8(Mst2))) + 2*lmMst2*(
        630*pow6(Mst1)*pow6(Mst2) + 12*pow4(Msq)*((-67 + 30*lmMsq + 90*pow2(
        lmMsq))*pow4(Mst1)*pow4(Mst2) + 3*(19 + 30*lmMsq)*pow2(Mst2)*pow6(Mst1)
        + 5*(43 - 18*lmMsq)*pow2(Mst1)*pow6(Mst2) + 112*pow8(Mst1) - 32*pow8(
        Mst2)) + 45*(-13 + 28*lmMsq)*pow4(Mst1)*pow8(Mst2) + 80*pow2(Msq)*(15*
        pow4(Mst2)*pow6(Mst1) + 2*(-14 + 15*lmMsq)*pow4(Mst1)*pow6(Mst2) - 15*
        pow2(Mst1)*pow8(Mst2))))))/pow4(Msq)))/(Mgl*Tbeta*pow4(Mst1)*pow5(Mst2)
        )))/2592. + Al4p*(-((2*Dmglst2*lmMst2 + Mgl + lmMst2*Mgl)*twoLoopFlag*(
        -(pow2(Mst2)*pow2(s2t)*pow4(Mst1)*(-8*(lmMst1 - lmMst2)*Tbeta*pow2(Mt)*
        pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 4*Mt*(-(MuSUSY*s2t) + 2*Mt*Tbeta)*
        pow2(Mst2)*pow2(Sbeta) + (1 - 2*lmMst1 + 2*lmMst2)*Tbeta*pow2(s2t)*
        pow2(Sbeta)*pow4(Mst2))) - pow2(Mst1)*pow4(Mst2)*(-4*Tbeta*pow2(Mt)*
        pow2(s2t)*((-1 + 2*lmMst1 - 2*lmMst2)*pow2(MuSUSY)*(-1 + pow2(Sbeta)) +
        4*pow2(Mst2)*pow2(Sbeta)) + pow2(Sbeta)*(16*MuSUSY*s2t*pow3(Mt) + 8*(-
        lmMst1 + lmMst2)*Mt*MuSUSY*pow2(Mst2)*pow3(s2t) - 16*Tbeta*pow4(Mt) + (
        1 + 2*lmMst1 - 2*lmMst2)*Tbeta*pow4(Mst2)*pow4(s2t))) + Tbeta*pow2(s2t)
        *(8*(lmMst1 - lmMst2)*pow2(Mt)*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + pow2(
        s2t)*pow2(Sbeta)*pow4(Mst2))*pow6(Mst1) + (-4*Tbeta*pow2(Mt)*pow2(s2t)*
        (pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 2*pow2(Mst2)*pow2(Sbeta)) + pow2(
        Sbeta)*(16*MuSUSY*s2t*pow3(Mt) - 4*Mt*MuSUSY*pow2(Mst2)*pow3(s2t) + 16*
        Tbeta*pow4(Mt) + Tbeta*pow4(Mst2)*pow4(s2t)))*pow6(Mst2)))/(6.*Mgl*
        Tbeta*pow2(Mst1)*pow2(Sbeta)*pow4(Mst2)) + (xDmglst2*pow2(Dmglst2)*(((2
        - 3*lmMst2)*twoLoopFlag*(4*Mt*MuSUSY*s2t*((-4*pow2(Mst1)*pow2(Mt) + 4*
        pow2(Mst2)*pow2(Mt) + 2*(lmMst1 - lmMst2)*pow2(Mst1)*pow2(Mst2)*pow2(
        s2t) + pow2(s2t)*(pow4(Mst1) - pow4(Mst2)))/Tbeta + (Mt*MuSUSY*s2t*((
        pow2(Mst1) + pow2(Mst2))*pow4(Mst2) - 2*(lmMst1 - lmMst2)*pow2(Mst1)*(
        pow2(Mst1)*pow2(Mst2) + pow4(Mst1) + pow4(Mst2))))/(pow2(Sbeta)*pow4(
        Mst2))) + 16*(pow2(Mst1) + pow2(Mst2))*pow4(Mt) + (pow2(Mst1) - pow2(
        Mst2))*(2*(lmMst1 - lmMst2)*pow2(Mst1)*pow2(Mst2) + pow4(Mst1) - pow4(
        Mst2))*pow4(s2t) - (4*pow2(Mt)*pow2(s2t)*(2*pow4(Mst1)*((-lmMst1 +
        lmMst2)*pow2(Mst2)*pow2(MuSUSY) + pow4(Mst2)) + 2*(-lmMst1 + lmMst2)*
        pow2(MuSUSY)*pow6(Mst1) + pow2(Mst1)*((1 - 2*lmMst1 + 2*lmMst2)*pow2(
        MuSUSY)*pow4(Mst2) - 4*pow6(Mst2)) + (2*pow2(Mst2) + pow2(MuSUSY))*
        pow6(Mst2)))/pow4(Mst2)))/(6.*pow2(Mst1)) + Al4p*threeLoopFlag*((Mt*
        pow3(s2t)*(25920*(-71 + 6*lmMst1*(-1 + lmMst2) + lmMst2 + 30*pow2(
        lmMst2))*pow2(Mst2)*pow4(Mst1) + 25920*pow2(Mst1)*(5 + 167*lmMst2 + 18*
        pow2(lmMst2) + 6*lmMst1*(-10 - 3*lmMst2 + 12*pow2(lmMst2)) - 72*pow3(
        lmMst2))*pow4(Mst2) + (48600*pow2(Mst2)*pow6(Mst1))/pow2(Msq) - (-
        26558551 + 291600*lmMsq*(-5 + 3*lmMst1 - 3*lmMst2) - 5328630*lmMst2 +
        2529760*OepS2 - 66700746*S2 + 51227640*lmMst2*S2 - 19440*(-61 + 24*
        lmMst2)*pow2(lmMst1) + 1652400*pow2(lmMst2) + 90*lmMst1*(79943 - 31536*
        lmMst2 - 569196*S2 + 5184*pow2(lmMst2)) + 155520*pow3(lmMst1) - 155520*
        pow3(lmMst2))*pow6(Mst1) + (12150*pow4(Mst2)*pow6(Mst1))/pow4(Msq) -
        311040*(-2 - lmMst2 + 3*pow2(lmMst2))*pow6(Mst2)))/(21870.*Mst2*pow2(
        Mst1)) - pow4(s2t)*(-((615080*OepS2 - 27*(669989 + 461310*lmMst1 -
        461310*lmMst2)*S2)*pow4(Mst1))/43740. + (146.7968136853605 - (80299571*
        lmMst2)/2.3814e6 - (5*lmMsq*(78*lmMst1 - 7*(5 + 6*lmMst2)))/36. + (5*
        pow2(lmMsq))/2. + (1.4137566137566138 + (143*lmMst2)/18.)*pow2(lmMst1)
        + lmMst1*(27.663925002099607 + (13708*lmMst2)/945. - (143*pow2(lmMst2))
        /18.) - (9584*pow2(lmMst2))/945. - (143*pow3(lmMst1))/54. + (143*pow3(
        lmMst2))/54.)*pow4(Mst1) - (pow2(Mst2)*(pow2(Msq)*pow2(Mst1)*(5 - 4312*
        lmMst2 + 270*lmMsq*(5 - 6*lmMst2)*(1 - 2*lmMst1 + 2*lmMst2) + 810*(1 -
        2*lmMst1 + 2*lmMst2)*pow2(lmMsq) + 1344*(2 - 3*lmMst2)*pow2(lmMst1) +
        4428*pow2(lmMst2) + 6*lmMst1*(1001 - 1418*lmMst2 + 528*pow2(lmMst2)) +
        864*pow3(lmMst2)) + 5*(91 - 12*lmMsq*(26 + 3*lmMst1 - 3*lmMst2) + 360*
        lmMst2 + 12*lmMst1*(-4 + 3*lmMst2) - 36*pow2(lmMst2))*pow4(Mst1)))/(
        216.*pow2(Msq)) - ((400*(-8 + 15*lmMsq - 15*lmMst2)*(1 - 2*lmMst1 + 2*
        lmMst2)*pow2(Msq)*pow2(Mst1) + 2*(-3739 + 7028*lmMst2 - 270*lmMsq*(1 +
        2*lmMst1 - 2*lmMst2)*(-5 + 6*lmMst2) + 810*(1 + 2*lmMst1 - 2*lmMst2)*
        pow2(lmMsq) - 192*(-2 + 3*lmMst2)*pow2(lmMst1) + 3924*pow2(lmMst2) + 6*
        lmMst1*(-1289 + 650*lmMst2 + 1008*pow2(lmMst2)) - 5472*pow3(lmMst2))*
        pow4(Msq) - 5*(-244 + lmMst1*(87 - 180*lmMst2) + 36*lmMsq*(17 + 5*
        lmMst1 - 5*lmMst2) - 699*lmMst2 + 180*pow2(lmMst2))*pow4(Mst1))*pow4(
        Mst2))/(432.*pow4(Msq)) + (((200*(1 + 2*lmMst1 - 2*lmMst2)*(8 - 15*
        lmMsq + 15*lmMst2))/pow2(Msq) + (3*(-1321 + 394*lmMst2 + 64*lmMst1*(-2
        + 3*lmMst2) - 90*lmMsq*(-5 + 6*lmMst2) + 270*pow2(lmMsq) + 1904*pow2(
        lmMst2)))/pow2(Mst1))*pow6(Mst2))/216. + ((8*(2 + lmMst2 - 3*pow2(
        lmMst2))*pow2(Msq) + 5*(-4 + 15*lmMsq - 15*lmMst2)*pow2(Mst1))*pow8(
        Mst2))/(9.*pow2(Msq)*pow4(Mst1)) + (5*(7*(1 - 2*lmMst1 + 2*lmMst2)*(5 -
        28*lmMsq + 28*lmMst2)*pow2(Mst1)*pow6(Mst2) + (1 + 2*lmMst1 - 2*lmMst2)
        *(9 - 140*lmMsq + 140*lmMst2)*pow8(Mst2)))/(96.*pow4(Msq))) + (pow3(Mt)
        *(-20*Mst2*s2t*pow2(Mst1)*(14175*(-5*(688 + 27*lmMst1 + 36*lmMsq*(-2 +
        5*lmMst1 - 5*lmMst2) + 45*lmMst2 - 180*lmMst1*lmMst2 + 180*pow2(lmMst2)
        )*pow2(Mst2)*pow2(MuSUSY) - 8*pow2(Msq)*(3*(-19 + 2*lmMst1 - 4*lmMst2 +
        2*lmMt)*pow2(Mst2) + (727 + 30*lmMst1 + 18*lmMsq*(-2 + 5*lmMst1 - 5*
        lmMst2) + 6*lmMst2 - 90*lmMst1*lmMst2 + 90*pow2(lmMst2))*pow2(MuSUSY))
        - 6*(5 + 2*lmMst1 - 4*lmMst2 + 2*lmMt)*pow4(Mst2))*pow6(Mst1) + pow4(
        Msq)*(-30240*pow2(Mst1)*pow2(Mst2)*(pow2(Mst2)*(300 + 72*lmMt + lmMst2*
        (-341 + 78*lmMt) + 18*lmMst1*(8 - 4*lmMt + lmMst2*(5 + 2*lmMt) - 2*
        pow2(lmMst2)) - 12*(26 + 3*lmMt)*pow2(lmMst2) + 36*pow3(lmMst2)) + 3*
        pow2(MuSUSY)*(-53 - 191*lmMst2 + lmMst1*(60 + 18*lmMst2 - 72*pow2(
        lmMst2)) + 54*pow2(lmMst2) + 72*pow3(lmMst2))) - 90720*(6*(93 + 88*
        lmMst2 - 22*lmMt - 8*lmMst2*lmMt + lmMst1*(-34 - 8*lmMst2 + 4*lmMt) +
        12*pow2(lmMst2))*pow2(Mst2) + pow2(MuSUSY)*(-11 - 371*lmMst2 + 42*pow2(
        lmMst2) - 6*lmMst1*(-21 - 5*lmMst2 + 24*pow2(lmMst2)) + 144*pow3(
        lmMst2)))*pow4(Mst1) + 1088640*(2 + lmMst2 - 3*pow2(lmMst2))*(2*pow2(
        Mst2) + pow2(MuSUSY))*pow4(Mst2) + (377667593 - 198080820*lmMst2 -
        45543680*OepS2 - 4968*(-135731 + 185640*lmMst2)*S2 - 136080*(339 + 62*
        lmMst2 - 59*lmMt)*pow2(lmMst1) - 136987200*pow2(lmMst2) + 1020600*
        lmMsq*(97 + 94*lmMst2 - 4*lmMst1*(17 + 3*lmMst2 - 3*lmMt) - 2*(13 + 6*
        lmMst2)*lmMt + 12*pow2(lmMst2)) + 5670*lmMt*(19949 + 19864*lmMst2 +
        4296*pow2(lmMst2)) - 1632960*(7 + 6*lmMst2)*pow2(lmMt) + 630*lmMst1*(-
        19783 - 114840*lmMt - 72*lmMst2*(-3401 + 714*lmMt) + 1463904*S2 +
        36504*pow2(lmMst2) + 15552*pow2(lmMt)) + 136080*pow3(lmMst1) -
        14696640*pow3(lmMst2))*pow6(Mst1))) + 189*Mt*(3000*pow2(Msq)*(20*(8 -
        15*lmMsq + 15*lmMst2)*pow2(Mst1)*(pow2(Mst1) + pow2(Mst2))*pow6(Mst2) +
        9*(127 + lmMst2*(68 - 6*lmMt) - 24*lmMt + lmMst1*(-44 - 6*lmMst2 + 6*
        lmMt) + 6*pow2(lmMst2))*pow2(Mst2)*pow8(Mst1)) + 4*pow4(Msq)*(-((316871
        + 20250*lmMsq*(5 - 6*lmMst2) - 480*(137 + 60*lmMst1)*lmMt + 6*lmMst2*(-
        18049 + 9600*lmMst1 + 12240*lmMt) + 60750*pow2(lmMsq) + 9360*pow2(
        lmMst2))*pow4(Mst1)*pow4(Mst2)) - 14400*pow2(Mst2)*(221 + lmMst2*(299 -
        47*lmMt) - 62*lmMt + (87 - 6*lmMt)*pow2(lmMst2) - 2*lmMst1*(65 +
        lmMst2*(32 - 3*lmMt) - 12*lmMt + 3*pow2(lmMst2)) + 6*pow3(lmMst2))*
        pow6(Mst1) + 225*(937 + 64*lmMst1*(2 - 3*lmMst2) + 502*lmMst2 + 90*
        lmMsq*(-5 + 6*lmMst2) - 270*pow2(lmMsq) - 624*pow2(lmMst2))*pow2(Mst1)*
        pow6(Mst2) + 28800*(42 + lmMst2*(242 - 62*lmMt) - 33*lmMt + 6*(29 - 4*
        lmMt)*pow2(lmMst2) - 2*lmMst1*(81 + lmMst2*(74 - 12*lmMt) - 18*lmMt +
        12*pow2(lmMst2)) + 24*pow3(lmMst2))*pow8(Mst1) + 14400*(-2 + lmMst2 +
        5*pow2(lmMst2))*pow8(Mst2)) - 3375*(-4*(-4 + lmMst1 - 2*lmMst2 + lmMt)*
        pow4(Mst2)*pow8(Mst1) + 7*(-5 + 28*lmMsq - 28*lmMst2)*pow2(Mst1)*(pow2(
        Mst1) + pow2(Mst2))*pow8(Mst2)))))/(765450.*pow4(Msq)*pow4(Mst1)*pow4(
        Mst2)) - (s2t*(8*Mst2*T1ep*(-58056*Mst2*s2t*pow2(Mt) - 126488*Mt*pow2(
        Mst2)*pow2(s2t) + 1301248*pow3(Mt) + 15377*pow3(Mst2)*pow3(s2t))*pow8(
        Mst1) + (27*pow2(Mt)*pow2(MuSUSY)*(5*pow2(Mst1)*pow3(Mst2)*(63*(-5 +
        28*lmMsq - 28*lmMst2)*s2t*(pow2(Mst1) + pow2(Mst2))*(-2*lmMst1*pow2(
        Mst1) + 2*lmMst2*pow2(Mst1) + pow2(Mst2))*pow3(Mst2) + 4*(-3*Mst2*s2t*(
        247 + 2*lmMst1*(-5 + 72*lmMsq - 72*lmMst2) + 10*lmMst2 - 144*lmMsq*
        lmMst2 + 144*pow2(lmMst2)) + 20*Mt*(688 + 27*lmMst1 + 36*lmMsq*(-2 + 5*
        lmMst1 - 5*lmMst2) + 45*lmMst2 - 180*lmMst1*lmMst2 + 180*pow2(lmMst2)))
        *pow6(Mst1)) + 20*pow2(Msq)*(40*(-8 + 15*lmMsq - 15*lmMst2)*s2t*pow2(
        Mst1)*(pow2(Mst1) + pow2(Mst2))*(-2*lmMst1*pow2(Mst1) + 2*lmMst2*pow2(
        Mst1) + pow2(Mst2))*pow4(Mst2) + Mst2*(Mst2*s2t*(-1011 - 298*lmMst2 +
        552*lmMsq*lmMst2 + lmMst1*(298 - 552*lmMsq + 552*lmMst2) - 552*pow2(
        lmMst2)) + 32*Mt*(727 + 30*lmMst1 + 18*lmMsq*(-2 + 5*lmMst1 - 5*lmMst2)
        + 6*lmMst2 - 90*lmMst1*lmMst2 + 90*pow2(lmMst2)))*pow8(Mst1)) + 4*pow4(
        Msq)*(-((Mst2*s2t*(3035 + 5240*lmMst2 + 270*lmMsq*(1 - 2*lmMst1 + 2*
        lmMst2)*(-5 + 6*lmMst2) + 810*(-1 + 2*lmMst1 - 2*lmMst2)*pow2(lmMsq) +
        192*(2 - 3*lmMst2)*pow2(lmMst1) - 4620*pow2(lmMst2) + 6*lmMst1*(-1161 +
        458*lmMst2 + 1008*pow2(lmMst2)) - 5472*pow3(lmMst2)) - 128*Mt*(-53 -
        191*lmMst2 + lmMst1*(60 + 18*lmMst2 - 72*pow2(lmMst2)) + 54*pow2(
        lmMst2) + 72*pow3(lmMst2)))*pow3(Mst2)*pow4(Mst1)) - 3*(Mst2*s2t*(1065
        + 64*lmMst1*(2 - 3*lmMst2) - 266*lmMst2 + 90*lmMsq*(-5 + 6*lmMst2) -
        270*pow2(lmMsq) - 1264*pow2(lmMst2)) + 512*Mt*(2 + lmMst2 - 3*pow2(
        lmMst2)))*pow2(Mst1)*pow5(Mst2) - 2*Mst2*(3*Mst2*s2t*(480 - 9*lmMst2*(-
        129 + 50*lmMsq + 30*pow2(lmMsq)) + 288*(2 - 3*lmMst2)*pow2(lmMst1) +
        10*(-17 + 54*lmMsq)*pow2(lmMst2) + lmMst1*(-1385 - 406*lmMst2 - 90*
        lmMsq*(-5 + 6*lmMst2) + 270*pow2(lmMsq) + 2544*pow2(lmMst2)) - 1680*
        pow3(lmMst2)) + 64*Mt*(11 + 371*lmMst2 - 42*pow2(lmMst2) + 6*lmMst1*(-
        21 - 5*lmMst2 + 24*pow2(lmMst2)) - 144*pow3(lmMst2)))*pow6(Mst1) + s2t*
        (12*(96 + lmMst2*(285 + 30*lmMsq + 90*pow2(lmMsq) + 544*pow2(lmMst1)) +
        (590 - 180*lmMsq)*pow2(lmMst2) - lmMst1*(109 + 30*lmMsq*(1 - 6*lmMst2)
        + 590*lmMst2 + 90*pow2(lmMsq) + 1360*pow2(lmMst2)) + 816*pow3(lmMst2))*
        pow8(Mst1) - 192*(-2 + lmMst2 + 5*pow2(lmMst2))*pow8(Mst2)))))/(pow2(
        Sbeta)*pow4(Msq))))/(5832.*pow4(Mst1)*pow4(Mst2)) + (pow2(Mt)*pow2(s2t)
        *(8*(1061 + 270*lmMsq*(5 - 6*lmMst2) + 192*lmMst1*(-11 + lmMst2) +
        2398*lmMst2 + 810*pow2(lmMsq) + 1296*pow2(lmMst2))*pow2(Mst1) - (32*(
        338660*OepS2 + 27*(182471 - 253995*lmMst1 + 253995*lmMst2)*S2)*pow4(
        Mst1))/(945.*pow2(Mst2)) - 216*((20*(4 - lmMst1 + lmMst2))/(3.*pow2(
        Msq)) - (781.8382197186929 - 240*lmMsq + (-23.457324263038547 + 180*
        lmMsq)*lmMst1 + lmMst2*(930.1239909297052 - 180*lmMsq + (520424*lmMst1)
        /945. - (304*pow2(lmMst1))/9.) - (162772*pow2(lmMst1))/945. + (4*(-
        89413 + 5460*lmMst1)*pow2(lmMst2))/945. + (32*lmMt*(54 + 65*lmMst2 -
        lmMst1*(65 + 18*lmMst2) + 9*pow2(lmMst1) + 9*pow2(lmMst2)))/9. + (16*
        pow3(lmMst1))/27. + (272*pow3(lmMst2))/27.)/pow2(Mst2))*pow4(Mst1) + (
        2*pow2(Mst2)*(800*(-8 + 15*lmMsq - 15*lmMst2)*pow2(Msq)*pow2(Mst1) + 8*
        (2507 + 514*lmMst2 + 270*lmMsq*(-5 + 6*lmMst2) - 810*pow2(lmMsq) + 96*(
        -2 + 3*lmMst2)*pow2(lmMst1) - 1584*pow2(lmMst2) - 96*lmMst1*(-5 + 2*
        lmMst2 + 6*pow2(lmMst2)) + 288*pow3(lmMst2))*pow4(Msq) + 15*(35 + 12*
        lmMsq - 12*lmMst2)*pow4(Mst1)))/pow4(Msq) + 8*((400*(8 - 15*lmMsq + 15*
        lmMst2))/pow2(Msq) + (3*(-1001 - 118*lmMst2 + 64*lmMst1*(-2 + 3*lmMst2)
        - 90*lmMsq*(-5 + 6*lmMst2) + 270*pow2(lmMsq) + 944*pow2(lmMst2)))/pow2(
        Mst1))*pow4(Mst2) - (1536*(-2 + lmMst2 + 5*pow2(lmMst2))*pow6(Mst2))/
        pow4(Mst1) - (630*(5 - 28*lmMsq + 28*lmMst2)*(pow2(Mst1) - 2*pow2(Mst2)
        )*pow4(Mst2) + (2*(800*(8 - 15*lmMsq + 15*lmMst2)*pow2(Msq)*pow6(Mst2)
        + 45*(9 - 140*lmMsq + 140*lmMst2)*pow8(Mst2)))/pow2(Mst1) + (pow2(
        MuSUSY)*(15*pow2(Mst1)*pow4(Mst2)*(21*(5 - 28*lmMsq + 28*lmMst2)*pow2(
        Mst2)*(pow2(Mst1) + pow2(Mst2))*(-2*lmMst1*pow2(Mst1) + 2*lmMst2*pow2(
        Mst1) + pow2(Mst2)) + (988 + 8*lmMst1*(-5 + 72*lmMsq - 72*lmMst2) + 8*(
        5 - 72*lmMsq)*lmMst2 + 576*pow2(lmMst2))*pow6(Mst1)) + 20*pow2(Msq)*(-
        40*(-8 + 15*lmMsq - 15*lmMst2)*pow2(Mst1)*(pow2(Mst1) + pow2(Mst2))*(-
        2*lmMst1*pow2(Mst1) + 2*lmMst2*pow2(Mst1) + pow2(Mst2))*pow4(Mst2) + (
        1011 + lmMst1*(-298 + 552*lmMsq - 552*lmMst2) + (298 - 552*lmMsq)*
        lmMst2 + 552*pow2(lmMst2))*pow2(Mst2)*pow8(Mst1)) - 4*pow4(Msq)*(-((
        3035 + 5240*lmMst2 + 270*lmMsq*(1 - 2*lmMst1 + 2*lmMst2)*(-5 + 6*
        lmMst2) + 810*(-1 + 2*lmMst1 - 2*lmMst2)*pow2(lmMsq) + 192*(2 - 3*
        lmMst2)*pow2(lmMst1) - 4620*pow2(lmMst2) + 6*lmMst1*(-1161 + 458*lmMst2
        + 1008*pow2(lmMst2)) - 5472*pow3(lmMst2))*pow4(Mst1)*pow4(Mst2)) - 6*
        pow2(Mst2)*(480 - 9*lmMst2*(-129 + 50*lmMsq + 30*pow2(lmMsq)) + 288*(2
        - 3*lmMst2)*pow2(lmMst1) + 10*(-17 + 54*lmMsq)*pow2(lmMst2) + lmMst1*(-
        1385 - 406*lmMst2 - 90*lmMsq*(-5 + 6*lmMst2) + 270*pow2(lmMsq) + 2544*
        pow2(lmMst2)) - 1680*pow3(lmMst2))*pow6(Mst1) + 3*(-1065 + 266*lmMst2 +
        64*lmMst1*(-2 + 3*lmMst2) - 90*lmMsq*(-5 + 6*lmMst2) + 270*pow2(lmMsq)
        + 1264*pow2(lmMst2))*pow2(Mst1)*pow6(Mst2) + 12*(96 + lmMst2*(285 + 30*
        lmMsq + 90*pow2(lmMsq) + 544*pow2(lmMst1)) + (590 - 180*lmMsq)*pow2(
        lmMst2) - lmMst1*(109 + 30*lmMsq*(1 - 6*lmMst2) + 590*lmMst2 + 90*pow2(
        lmMsq) + 1360*pow2(lmMst2)) + 816*pow3(lmMst2))*pow8(Mst1) - 192*(-2 +
        lmMst2 + 5*pow2(lmMst2))*pow8(Mst2))))/(pow4(Mst1)*pow4(Mst2)))/pow4(
        Msq)))/216. + (MuSUSY*(-((pow2(Mt)*pow2(s2t)*(-((322696*OepS2 - 81*(
        86631 + 80674*lmMst1 - 80674*lmMst2)*S2)*pow4(Mst1))/729. + (
        4873.938614540466 + 152*B4 - 4*DN - 560*lmMsq - (1665.7654320987654 +
        200*lmMsq)*lmMst1 - 224*pow2(lmMst1) + lmMst2*(3137.7654320987654 +
        200*lmMsq + 768*lmMst1 + 64*pow2(lmMst1)) + 32*(-17 + 2*lmMst1)*pow2(
        lmMst2) - (64*pow3(lmMst1))/3. - (320*pow3(lmMst2))/3.)*pow4(Mst1) - (
        4*pow2(Mst2)*(16*pow2(Msq)*pow2(Mst1)*(7 - 30*lmMst2 + lmMst1*(11 + 2*
        lmMst2 - 12*pow2(lmMst2)) - 2*pow2(lmMst2) + 12*pow3(lmMst2)) + 5*(157
        + 8*lmMst1*(1 + 3*lmMsq - 3*lmMst2) - 8*(1 + 3*lmMsq)*lmMst2 + 24*pow2(
        lmMst2))*pow4(Mst1)))/(3.*pow2(Msq)) + ((32*(29 + 179*lmMst2 - 18*pow2(
        lmMst2) + 6*lmMst1*(-10 - 3*lmMst2 + 12*pow2(lmMst2)) - 72*pow3(lmMst2)
        ) - (15*(355 + 6*lmMst1*(3 + 20*lmMsq - 20*lmMst2) - 6*(3 + 20*lmMsq)*
        lmMst2 + 120*pow2(lmMst2))*pow4(Mst1))/pow4(Msq))*pow4(Mst2))/9. + (
        128*(2 + lmMst2 - 3*pow2(lmMst2))*pow6(Mst2))/(3.*pow2(Mst1))))/pow3(
        Mst2)) + (pow3(Mt)*(8*Mst2*Mt*pow2(Mst1)*(-32*pow4(Msq)*(pow2(Mst1)*
        pow2(Mst2)*(96 + 36*lmMt + lmMst2*(-187 + 39*lmMt) - 3*(19 + 6*lmMt)*
        pow2(lmMst2) - 18*lmMst1*(-3 + 2*lmMt - lmMst2*(3 + lmMt) + pow2(
        lmMst2)) + 18*pow3(lmMst2)) + (834 + lmMst2*(569 - 33*lmMt) - 162*lmMt
        + (51 - 18*lmMt)*pow2(lmMst2) - 18*lmMst1*(13 + lmMst2 - lmMst2*lmMt +
        pow2(lmMst2)) + 18*pow3(lmMst2))*pow4(Mst1) - 36*(2 + lmMst2 - 3*pow2(
        lmMst2))*pow4(Mst2)) + (60*(167 - 12*lmMst1*(2 + 3*lmMsq - 3*lmMst2) +
        6*(5 + 6*lmMsq)*lmMst2 - 6*lmMt - 36*pow2(lmMst2))*pow2(Msq) + 45*(235
        - 24*lmMst1*(1 + 3*lmMsq - 3*lmMst2) + 24*(1 + 3*lmMsq)*lmMst2 - 72*
        pow2(lmMst2))*pow2(Mst2))*pow6(Mst1)) - 3*s2t*(60*pow4(Mst2)*pow8(Mst1)
        - 80*pow2(Msq)*(10*(-8 + 15*lmMsq - 15*lmMst2)*pow2(Mst1)*(pow2(Mst1) -
        pow2(Mst2))*pow6(Mst2) + (33 - 9*lmMst1 + 9*lmMst2)*pow2(Mst2)*pow8(
        Mst1)) + 315*(5 - 28*lmMsq + 28*lmMst2)*pow2(Mst1)*(pow2(Mst1) - pow2(
        Mst2))*pow8(Mst2) - 4*pow4(Msq)*(-((2395 + 482*lmMst2 + 270*lmMsq*(-5 +
        6*lmMst2) - 810*pow2(lmMsq) + 192*(-2 + 3*lmMst2)*pow2(lmMst1) + 192*
        lmMst1*(3 + lmMst2 - 6*pow2(lmMst2)) - 1296*pow2(lmMst2) + 576*pow3(
        lmMst2))*pow4(Mst1)*pow4(Mst2)) + 192*pow2(Mst2)*((2 - 3*lmMst2)*pow2(
        lmMst1) + lmMst1*(8 - 2*lmMst2 + 6*pow2(lmMst2)) - 3*(6 + 5*lmMst2 +
        pow3(lmMst2)))*pow6(Mst1) + 3*(873 + 64*lmMst1*(2 - 3*lmMst2) + 182*
        lmMst2 + 90*lmMsq*(-5 + 6*lmMst2) - 270*pow2(lmMsq) - 624*pow2(lmMst2))
        *pow2(Mst1)*pow6(Mst2) - 384*(-6 - 14*lmMst2 + lmMst2*pow2(lmMst1) +
        lmMst1*(9 + 6*lmMst2 - 2*pow2(lmMst2)) - 6*pow2(lmMst2) + pow3(lmMst2))
        *pow8(Mst1) + 192*(-2 + lmMst2 + 5*pow2(lmMst2))*pow8(Mst2)))))/(162.*
        pow4(Msq)*pow4(Mst1)*pow4(Mst2)) + (Mt*pow2(s2t)*((-16*(Mst2*s2t*(
        376960*OepS2 - 27*(1307689 + 282720*lmMst1 - 282720*lmMst2)*S2 -
        565440*T1ep) + 3630330*Mt*T1ep)*pow4(Mst1))/pow3(Mst2) + 405*s2t*(4*
        pow2(Mst1)*(155 - 1726*lmMst2 + 270*lmMsq*(-5 + 6*lmMst2) - 810*pow2(
        lmMsq) + 1536*(-2 + 3*lmMst2)*pow2(lmMst1) + 192*lmMst1*(7 + 27*lmMst2
        - 48*pow2(lmMst2)) - 3600*pow2(lmMst2) + 4608*pow3(lmMst2)) + (2*pow2(
        Mst2)*(-400*(-8 + 15*lmMsq - 15*lmMst2)*pow2(Msq)*pow2(Mst1) + 4*(80 +
        lmMst2*(-3019 + 1350*lmMsq + 810*pow2(lmMsq)) + 96*(-2 + 3*lmMst2)*
        pow2(lmMst1) - 18*(-23 + 90*lmMsq)*pow2(lmMst2) - 3*lmMst1*(-1225 +
        554*lmMst2 - 90*lmMsq*(-5 + 6*lmMst2) + 270*pow2(lmMsq) + 1008*pow2(
        lmMst2)) + 2736*pow3(lmMst2))*pow4(Msq) - 15*pow4(Mst1)))/pow4(Msq) +
        216*((5*(19 - 6*lmMst1 + 6*lmMst2))/(18.*pow2(Msq)) + (
        865.3922396522787 + (64*B4)/9. - (64*D3)/9. + (32*DN)/9. + (145*lmMsq)/
        3. + (72.3524901318552 - 50*lmMsq)*lmMst1 - (3673*pow2(lmMst1))/189. +
        lmMst2*(-142.57471235407743 + 50*lmMsq + (13058*lmMst1)/189. + (256*
        pow2(lmMst1))/3.) - (5*(1877 + 5376*lmMst1)*pow2(lmMst2))/189. - (256*
        pow3(lmMst1))/27. + (1792*pow3(lmMst2))/27.)/pow2(Mst2))*pow4(Mst1) + (
        (1600*(lmMst1 - lmMst2)*(8 - 15*lmMsq + 15*lmMst2))/pow2(Msq) + (12*(-
        1193 + 330*lmMst2 + 64*lmMst1*(-2 + 3*lmMst2) - 90*lmMsq*(-5 + 6*
        lmMst2) + 270*pow2(lmMsq) + 1584*pow2(lmMst2)))/pow2(Mst1) + (315*(5 -
        28*lmMsq + 28*lmMst2)*pow2(Mst1))/pow4(Msq))*pow4(Mst2) + ((2*(400*(-8
        + 15*lmMsq - 15*lmMst2)*pow2(Msq)*pow2(Mst1) - 384*(-2 + lmMst2 + 5*
        pow2(lmMst2))*pow4(Msq) + 315*(lmMst1 - lmMst2)*(5 - 28*lmMsq + 28*
        lmMst2)*pow4(Mst1))*pow6(Mst2))/pow4(Mst1) + (45*(-9 + 140*lmMsq - 140*
        lmMst2)*pow8(Mst2))/pow2(Mst1))/pow4(Msq))))/87480.))/Tbeta)))/pow2(
        Mgl))) - (z2*(-(Al4p*twoLoopFlag*(-4*Dmglst2*(4*pow2(Mt)*(12*(MuSUSY +
        2*Mst2*Tbeta)*pow2(Mt)*pow2(Sbeta) + Mst2*MuSUSY*pow2(s2t)*(MuSUSY*
        Tbeta*(-1 + pow2(Sbeta)) - 6*Mst2*pow2(Sbeta)) + Mt*s2t*Tbeta*(-13*
        pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 12*pow2(Mst2)*pow2(Sbeta)))*pow4(
        Mst1) + pow4(Mst2)*(-2*Mst2*MuSUSY*pow2(Mt)*pow2(s2t)*(-2*MuSUSY*Tbeta*
        (-1 + pow2(Sbeta)) + 15*Mst2*pow2(Sbeta)) - 4*s2t*Tbeta*(5*pow2(MuSUSY)
        *(-1 + pow2(Sbeta)) - 4*pow2(Mst2)*pow2(Sbeta))*pow3(Mt) + 2*Mt*(2*
        MuSUSY + 5*Mst2*Tbeta)*pow2(Sbeta)*pow3(Mst2)*pow3(s2t) - 16*MuSUSY*
        pow2(Sbeta)*pow4(Mt) - Tbeta*pow2(Sbeta)*pow4(s2t)*pow5(Mst2)) + pow2(
        Mst1)*pow2(Mst2)*(-8*s2t*(3*MuSUSY*s2t + 2*Mt*Tbeta)*pow2(Mst2)*pow2(
        Mt)*pow2(Sbeta) + 4*Mst2*Tbeta*pow2(Mt)*(pow2(MuSUSY)*pow2(s2t)*(-1 +
        pow2(Sbeta)) + 8*pow2(Mt)*pow2(Sbeta)) - 36*s2t*Tbeta*pow2(MuSUSY)*(-1
        + pow2(Sbeta))*pow3(Mt) - 2*Mt*Tbeta*pow2(Sbeta)*pow3(s2t)*pow4(Mst2) +
        Tbeta*pow2(Sbeta)*pow4(s2t)*pow5(Mst2))) + Mgl*(-2*pow2(Mst1)*pow2(
        Mst2)*(8*s2t*Tbeta*(pow2(MuSUSY) + 4*pow2(Mst2)*pow2(Sbeta) - pow2(
        MuSUSY)*pow2(Sbeta))*pow3(Mt) - 2*Mt*(MuSUSY + 2*Mst2*Tbeta)*pow2(
        Sbeta)*pow3(Mst2)*pow3(s2t) - 32*(2*MuSUSY + Mst2*Tbeta)*pow2(Sbeta)*
        pow4(Mt) + Tbeta*pow2(Sbeta)*pow4(s2t)*pow5(Mst2)) + pow4(Mst2)*(4*
        Mst2*MuSUSY*pow2(Mt)*pow2(s2t)*(-(MuSUSY*Tbeta*(-1 + pow2(Sbeta))) + 6*
        Mst2*pow2(Sbeta)) + 16*s2t*Tbeta*(pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 4*
        pow2(Mst2)*pow2(Sbeta))*pow3(Mt) - 4*Mt*(MuSUSY + 2*Mst2*Tbeta)*pow2(
        Sbeta)*pow3(Mst2)*pow3(s2t) + 32*(2*MuSUSY + Mst2*Tbeta)*pow2(Sbeta)*
        pow4(Mt) + Tbeta*pow2(Sbeta)*pow4(s2t)*pow5(Mst2)) + pow4(Mst1)*(16*
        s2t*Tbeta*(pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 4*pow2(Mst2)*pow2(Sbeta))*
        pow3(Mt) + 96*(2*MuSUSY + Mst2*Tbeta)*pow2(Sbeta)*pow4(Mt) + Tbeta*
        pow2(Sbeta)*pow4(s2t)*pow5(Mst2)))))/6. + (threeLoopFlag*pow2(Al4p)*(
        84*Mt*MuSUSY*pow2(Sbeta)*(4*pow2(Mst1)*pow3(Mt)*(9*Dmglst2*pow2(Mst2)*(
        (8*(57943 + 1800*lmMsq - 2850*lmMst1 + 2190*lmMst2 - 1260*lmMt)*pow2(
        Mst1) + 15*(1805 + 960*lmMsq + 248*lmMst1 - 2872*lmMst2 - 448*lmMt)*
        pow2(Mst2))*pow4(Msq) + 2400*(14*pow2(Msq) + 11*pow2(Mst2))*(pow2(Mst1)
        *pow2(Mst2) + pow4(Mst1) + pow4(Mst2))) + Mgl*(64800*(2*pow2(Msq)*pow2(
        Mst2) + pow4(Mst2))*(pow2(Mst1)*pow2(Mst2) + pow4(Mst1) + pow4(Mst2)) -
        pow4(Msq)*(48*(45986 - 2700*lmMsq - 7245*lmMst1 + 21015*lmMst2 + 630*
        lmMt)*pow2(Mst1)*pow2(Mst2) + 4*(1126591 - 184140*lmMst1 + 383940*
        lmMst2 + 7560*lmMt)*pow4(Mst1) + 9*(88373 - 14400*lmMsq - 9960*lmMst1 +
        52680*lmMst2 + 3360*lmMt)*pow4(Mst2)))) + pow3(Mst2)*pow3(s2t)*(180*
        pow2(Msq)*((Dmglst2*(3599 - 1080*lmMsq - 1368*lmMst1 + 3480*lmMst2) - (
        235 + 450*lmMsq + 642*lmMst1 - 1362*lmMst2)*Mgl)*pow2(Msq)*pow2(Mst1) +
        6*((92*Dmglst2 + 49*Mgl)*pow2(Msq) - 5*(34*Dmglst2 + 7*Mgl)*pow2(Mst1))
        *pow2(Mst2) + 90*(2*Dmglst2 + Mgl)*pow4(Mst1))*pow4(Mst2) + pow4(Msq)*(
        6*(6*Dmglst2*(33563 - 4020*lmMst1 + 4020*lmMst2) + (-242149 + 2700*
        lmMsq + 9360*lmMst1 - 19440*lmMst2)*Mgl)*pow2(Mst2)*pow4(Mst1) + ((16*
        Dmglst2*(776573 - 80595*lmMst1 + 80595*lmMst2) + 3*(-2205079 + 76140*
        lmMst1 - 76140*lmMst2)*Mgl)*pow6(Mst1))/3.) + 1350*pow2(Mst1)*(5*(4*
        Dmglst2 + Mgl)*pow2(Mst1) - (64*Dmglst2 + 9*Mgl)*pow2(Mst2))*pow6(Mst2)
        ) - 3*Mt*pow2(Mst1)*pow2(Mst2)*pow2(s2t)*(15*pow2(Msq)*pow2(Mst1)*((
        Dmglst2*(260504.6 - 43200*lmMsq - 21240*lmMst1 + 85752*lmMst2) + (15395
        - 8640*lmMsq + 17928*lmMst1 - 8712*lmMst2)*Mgl)*pow2(Msq) - 2880*(9*
        Dmglst2 + Mgl)*pow2(Mst1))*pow2(Mst2) + 2*(Dmglst2*(3.3195316666666665e6
        - 324000*lmMsq - 547560*lmMst1 + 1031400*lmMst2) + (63043 - 64800*
        lmMsq + 217080*lmMst1 - 147960*lmMst2)*Mgl)*pow4(Msq)*pow4(Mst1) + 135*
        (-320*(9*Dmglst2 + Mgl)*pow2(Msq)*pow2(Mst1) + (Dmglst2*(
        18556.466666666667 - 5280*lmMsq + 248*lmMst1 + 7944*lmMst2) + (3469 -
        1440*lmMsq + 664*lmMst1 + 1384*lmMst2)*Mgl)*pow4(Msq) - 80*(13*Dmglst2
        + Mgl)*pow4(Mst1))*pow4(Mst2) - 5400*(Mgl*(16*pow2(Msq) + 2*pow2(Mst1)
        + 5*pow2(Mst2)) + Dmglst2*(96*pow2(Msq) + 26*pow2(Mst1) + 41*pow2(Mst2)
        ))*pow6(Mst2)) - 2*Mst2*s2t*pow2(Mt)*pow4(Msq)*(-(Mgl*(3*(502399 +
        65160*lmMst1 - 134280*lmMst2)*pow2(Mst2)*pow4(Mst1) + 9*(-29683 +
        13560*lmMst1 - 36600*lmMst2)*pow2(Mst1)*pow4(Mst2) + (4828526 + 406080*
        lmMst1 - 613440*lmMst2)*pow6(Mst1) - 105840*pow6(Mst2))) + 24*Dmglst2*(
        6*(21379 + 1770*lmMst1 + 12630*lmMst2)*pow2(Mst2)*pow4(Mst1) + 6*(4879
        - 270*lmMst1 + 8910*lmMst2)*pow2(Mst1)*pow4(Mst2) + (520939 + 44280*
        lmMst1 + 76680*lmMst2)*pow6(Mst1) + 8280*pow6(Mst2))) - 3240*Mgl*s2t*
        shiftst3*(4*pow2(Mt) + ((1 + lmMst1 - lmMst2)*pow2(Mst1) - pow2(Mst2))*
        pow2(s2t))*pow4(Msq)*pow7(Mst2)) + Tbeta*(-28*pow2(Mt)*pow2(MuSUSY)*(2*
        Dmglst2*Mst2*(8100*Mst2*s2t*pow2(Mst1)*((186*Mt - 11*Mst2*s2t)*pow2(
        Mst2)*pow4(Mst1) + (134*Mt - 11*Mst2*s2t)*pow2(Mst1)*pow4(Mst2) + 2*
        pow2(Msq)*(-14*(-12*Mt + Mst2*s2t)*pow2(Mst1)*pow2(Mst2) + 2*(120*Mt -
        7*Mst2*s2t)*pow4(Mst1) + (96*Mt - 17*Mst2*s2t)*pow4(Mst2)) + 2*(41*Mt -
        8*Mst2*s2t)*pow6(Mst2)) + pow4(Msq)*(-27*pow2(Mst1)*pow2(Mst2)*((278347
        - 79200*lmMsq + 3720*lmMst1 + 119160*lmMst2)*Mst2*Mt*s2t + 7680*(5 + 6*
        lmMst2)*pow2(Mt) + 10*(-4151 + 1080*lmMsq + 1368*lmMst1 - 3480*lmMst2)*
        pow2(Mst2)*pow2(s2t)) - 36*Mst2*s2t*((534391 - 113400*lmMsq - 23760*
        lmMst1 + 196560*lmMst2)*Mt + 9*(-9053 + 900*lmMsq + 1810*lmMst1 - 3570*
        lmMst2)*Mst2*s2t)*pow4(Mst1) + pow2(s2t)*(4*(2286439 - 72900*lmMsq -
        307800*lmMst1 + 450360*lmMst2)*pow6(Mst1) + 149040*pow6(Mst2)))) + 3*
        Mgl*(1350*s2t*pow2(Mst1)*pow2(Mst2)*((72*Mt - 4*Mst2*s2t)*pow2(Mst2)*
        pow4(Mst1) - 4*(-14*Mt + Mst2*s2t)*pow2(Mst1)*pow4(Mst2) + 4*pow2(Msq)*
        (-4*(-12*Mt + Mst2*s2t)*pow2(Mst1)*pow2(Mst2) + (64*Mt - 4*Mst2*s2t)*
        pow4(Mst1) + (32*Mt - 7*Mst2*s2t)*pow4(Mst2)) + (40*Mt - 9*Mst2*s2t)*
        pow6(Mst2)) + pow4(Msq)*(-90*pow2(Mst1)*(3*(3469 - 1440*lmMsq + 664*
        lmMst1 + 1384*lmMst2)*Mst2*Mt*s2t + 64*(53 + 24*lmMst2)*pow2(Mt) + 2*(-
        59 + 450*lmMsq + 642*lmMst1 - 1362*lmMst2)*pow2(Mst2)*pow2(s2t))*pow3(
        Mst2) + s2t*(-6*(40*(5827 - 2700*lmMsq + 2988*lmMst1 + 468*lmMst2)*Mt +
        (240379 + 10800*lmMsq + 9900*lmMst1 - 21420*lmMst2)*Mst2*s2t)*pow2(
        Mst2)*pow4(Mst1) - (4*(412663 - 226800*lmMsq + 396360*lmMst1 - 119880*
        lmMst2)*Mt + (3647353 + 64800*lmMsq - 16740*lmMst1 - 52380*lmMst2)*
        Mst2*s2t)*pow6(Mst1)) + 52920*pow2(s2t)*pow7(Mst2)))) + 3780*Mst2*(-18*
        Mgl*shiftst3*pow4(Msq)*(-4*pow2(Mst1)*pow2(s2t)*pow4(Mst2)*(2*(-1 +
        lmMst1 - lmMst2)*pow2(Mt)*(-(pow2(MuSUSY)*(-1 + pow2(Sbeta))) + pow2(
        Mst2)*pow2(Sbeta)) + pow2(s2t)*pow2(Sbeta)*pow4(Mst2)) + 8*(-1 + 3*
        lmMst1 - 3*lmMst2)*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta))*
        pow6(Mst1) + (-4*pow2(Mt)*pow2(s2t)*(pow2(MuSUSY)*(-1 + pow2(Sbeta)) +
        2*pow2(Mst2)*pow2(Sbeta)) + pow2(Sbeta)*(16*pow4(Mt) + pow4(Mst2)*pow4(
        s2t)))*pow6(Mst2) + pow4(Mst1)*(8*(-1 + 2*lmMst1 - 2*lmMst2)*pow2(Mst2)
        *pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta)) + (3 + 2*lmMst1 -
        2*lmMst2)*pow2(Sbeta)*pow4(s2t)*pow6(Mst2))) + pow4(Mst1)*(8*(Dmglst2*(
        7286 - 240*lmMst1 + 6384*lmMst2) + (-353 + 72*lmMst1 + 696*lmMst2)*Mgl)
        *pow2(MuSUSY)*pow4(Msq)*pow4(Mt) - 105*(6*Dmglst2 + Mgl)*pow2(Sbeta)*
        pow4(s2t)*power10(Mst2))) + pow2(Sbeta)*(-21*pow4(s2t)*pow5(Mst2)*(6*
        pow2(Msq)*((24*Dmglst2*(3892 + 1350*lmMsq + 705*lmMst1 - 3345*lmMst2) +
        (-235099 + 16200*lmMsq + 28620*lmMst1 - 60300*lmMst2)*Mgl)*pow2(Msq) -
        2700*(2*Dmglst2 + Mgl)*pow2(Mst1))*pow2(Mst2)*pow4(Mst1) - 90*pow2(
        Mst1)*(Mgl*(-600*pow2(Msq)*pow2(Mst1) + 2*(529 + 450*lmMsq + 642*lmMst1
        - 1362*lmMst2)*pow4(Msq) + 75*pow4(Mst1)) + Dmglst2*(-2400*pow2(Msq)*
        pow2(Mst1) + (-6094 + 2160*lmMsq + 2736*lmMst1 - 6960*lmMst2)*pow4(Msq)
        + 300*pow4(Mst1)))*pow4(Mst2) + ((Dmglst2*(8800364 - 855360*lmMst1 +
        855360*lmMst2) - 15*(150437 + 3240*lmMsq - 3996*lmMst1 - 8100*lmMst2)*
        Mgl)*pow4(Msq)*pow6(Mst1))/3. + 1080*pow2(Msq)*((92*Dmglst2 + 49*Mgl)*
        pow2(Msq) - 5*(34*Dmglst2 + 7*Mgl)*pow2(Mst1))*pow6(Mst2)) + 84*Mt*
        pow2(Mst1)*pow3(s2t)*pow4(Mst2)*(pow4(Msq)*(30*(Dmglst2*(46748.2 +
        2160*lmMsq - 11736*lmMst1 + 7128*lmMst2) + (-7913 + 2160*lmMsq + 5976*
        lmMst1 - 10584*lmMst2)*Mgl)*pow2(Mst1)*pow2(Mst2) + ((Dmglst2*(8194483
        - 2329560*lmMst1 + 2329560*lmMst2) + 51*(-6167 + 9720*lmMst1 - 9720*
        lmMst2)*Mgl)*pow4(Mst1))/3.) + 135*pow2(Msq)*((Dmglst2*(
        18556.466666666667 - 5280*lmMsq + 248*lmMst1 + 7944*lmMst2) + (3469 -
        1440*lmMsq + 664*lmMst1 + 1384*lmMst2)*Mgl)*pow2(Msq) + 320*(3*Dmglst2
        + Mgl)*pow2(Mst1) - 640*(6*Dmglst2 + Mgl)*pow2(Mst2))*pow4(Mst2) +
        5400*(3*(5*Dmglst2 + Mgl)*pow2(Mst1) - 41*Dmglst2*pow2(Mst2))*pow6(
        Mst2)) - 28*Mst2*pow2(Mt)*pow2(s2t)*(-3*Mgl*pow4(Msq)*(-6*pow2(Mst2)*(
        4*(73931 + 3060*lmMst1 - 3060*lmMst2)*pow2(Mst2) + (240379 + 10800*
        lmMsq + 9900*lmMst1 - 21420*lmMst2)*pow2(MuSUSY))*pow4(Mst1) + 9*pow2(
        Mst1)*((17923 - 13560*lmMst1 + 36600*lmMst2)*pow2(Mst2) + 20*(59 - 450*
        lmMsq - 642*lmMst1 + 1362*lmMst2)*pow2(MuSUSY))*pow4(Mst2) - ((3321329
        + 210600*lmMst1 - 210600*lmMst2)*pow2(Mst2) + (3647353 + 64800*lmMsq -
        16740*lmMst1 - 52380*lmMst2)*pow2(MuSUSY))*pow6(Mst1) + 52920*(2*pow2(
        Mst2) + pow2(MuSUSY))*pow6(Mst2)) - 4*Dmglst2*pow4(Msq)*(162*pow2(Mst2)
        *(40*(275 + 34*lmMst1 + 62*lmMst2)*pow2(Mst2) + (9053 - 900*lmMsq -
        1810*lmMst1 + 3570*lmMst2)*pow2(MuSUSY))*pow4(Mst1) + 27*pow2(Mst1)*(4*
        (3499 - 270*lmMst1 + 8910*lmMst2)*pow2(Mst2) + 5*(4151 - 1080*lmMsq -
        1368*lmMst1 + 3480*lmMst2)*pow2(MuSUSY))*pow4(Mst2) + 2*(45*(78533 +
        6732*lmMst1 + 180*lmMst2)*pow2(Mst2) + (2286439 - 72900*lmMsq - 307800*
        lmMst1 + 450360*lmMst2)*pow2(MuSUSY))*pow6(Mst1) + 74520*(2*pow2(Mst2)
        + pow2(MuSUSY))*pow6(Mst2)) + 4050*Mgl*pow2(Mst1)*pow2(MuSUSY)*(4*pow4(
        Mst1)*pow4(Mst2) + 4*pow2(Mst1)*pow6(Mst2) + 4*pow2(Msq)*(4*pow2(Mst2)*
        pow4(Mst1) + 4*pow2(Mst1)*pow4(Mst2) + 7*pow6(Mst2)) + 9*pow8(Mst2)) +
        16200*Dmglst2*pow2(Mst1)*pow2(MuSUSY)*(11*pow4(Mst1)*pow4(Mst2) + 11*
        pow2(Mst1)*pow6(Mst2) + pow2(Msq)*(28*pow2(Mst2)*pow4(Mst1) + 28*pow2(
        Mst1)*pow4(Mst2) + 34*pow6(Mst2)) + 16*pow8(Mst2))) - 16*Mst2*pow4(Mt)*
        (-4*Dmglst2*pow4(Msq)*(27*pow2(Mst1)*pow2(Mst2)*((114841 - 1680*lmMst1
        + 31920*lmMst2 + 23520*lmMt)*pow2(Mst2) - 6720*(5 + 6*lmMst2)*pow2(
        MuSUSY)) + 9*((2655197 + 75600*lmMsq - 325080*lmMst1 + 778680*lmMst2 -
        10080*lmMt)*pow2(Mst2) + 105*(-3643 + 120*lmMst1 - 3192*lmMst2)*pow2(
        MuSUSY))*pow4(Mst1) + 4*(17026487 + 1530900*lmMsq - 3163860*lmMst1 +
        4320540*lmMst2 - 340200*lmMt)*pow6(Mst1) - 521640*pow6(Mst2)) - 453600*
        Dmglst2*pow2(Mst1)*(3*pow2(Msq)*pow2(Mst2)*(pow4(Mst1) - pow4(Mst2)) -
        4*pow8(Mst2)) + 21*Mgl*(pow4(Msq)*(9*pow2(Mst1)*pow2(Mst2)*((54383 +
        3600*lmMsq - 12480*lmMst1 + 25020*lmMst2 - 13980*lmMt)*pow2(Mst2) +
        160*(53 + 24*lmMst2)*pow2(MuSUSY)) + 3*((310039 + 21600*lmMsq - 91080*
        lmMst1 + 123480*lmMst2 - 43920*lmMt)*pow2(Mst2) + 30*(-353 + 72*lmMst1
        + 696*lmMst2)*pow2(MuSUSY))*pow4(Mst1) + 2*(610669 + 145800*lmMsq -
        300240*lmMst1 + 249480*lmMst2 - 78840*lmMt)*pow6(Mst1) + 52920*pow6(
        Mst2)) + 10800*pow2(Msq)*(3*pow2(Mst2)*pow6(Mst1) + 5*pow2(Mst1)*pow6(
        Mst2)) + 37800*pow2(Mst1)*pow8(Mst2))) - 56*s2t*pow2(Mst1)*pow3(Mt)*(
        Dmglst2*pow2(Mst2)*(pow4(Msq)*(27*pow2(Mst2)*(30*(2029 + 960*lmMsq +
        248*lmMst1 - 2872*lmMst2 - 448*lmMt)*pow2(Mst2) + (278347 - 79200*lmMsq
        + 3720*lmMst1 + 119160*lmMst2)*pow2(MuSUSY)) + 18*pow2(Mst1)*(3*(438149
        - 26520*lmMst1 + 60600*lmMst2 - 3360*lmMt)*pow2(Mst2) + 2*(534391 -
        113400*lmMsq - 23760*lmMst1 + 196560*lmMst2)*pow2(MuSUSY)) + 8*(9227767
        + 291600*lmMsq - 754920*lmMst1 + 1088640*lmMst2 - 22680*lmMt)*pow4(
        Mst1)) + 129600*pow2(Msq)*(-3*pow2(MuSUSY)*(7*pow2(Mst1)*pow2(Mst2) +
        10*pow4(Mst1) + 4*pow4(Mst2)) + 14*pow6(Mst2)) + 16200*(-(pow2(MuSUSY)*
        (93*pow2(Mst2)*pow4(Mst1) + 67*pow2(Mst1)*pow4(Mst2) + 41*pow6(Mst2)))
        + 88*pow8(Mst2))) + 3*Mgl*(-(pow4(Msq)*(6*pow2(Mst1)*pow2(Mst2)*((
        470657 - 86040*lmMst1 + 178200*lmMst2)*pow2(Mst2) - 20*(5827 - 2700*
        lmMsq + 2988*lmMst1 + 468*lmMst2)*pow2(MuSUSY)) + 2*(4*(574759 + 32400*
        lmMsq - 97200*lmMst1 + 131760*lmMst2)*pow2(Mst2) + (-412663 + 226800*
        lmMsq - 396360*lmMst1 + 119880*lmMst2)*pow2(MuSUSY))*pow4(Mst1) + 9*(2*
        (86693 - 14400*lmMsq - 9960*lmMst1 + 52680*lmMst2 + 3360*lmMt)*pow2(
        Mst2) - 15*(3469 - 1440*lmMsq + 664*lmMst1 + 1384*lmMst2)*pow2(MuSUSY))
        *pow4(Mst2))) + 43200*pow2(Msq)*(-(pow2(MuSUSY)*(4*pow2(Mst2)*pow4(
        Mst1) + 3*pow2(Mst1)*pow4(Mst2) + 2*pow6(Mst2))) + 6*pow8(Mst2)) +
        5400*(-(pow2(MuSUSY)*(9*pow4(Mst1)*pow4(Mst2) + 7*pow2(Mst1)*pow6(Mst2)
        + 5*pow8(Mst2))) + 24*power10(Mst2))))))))/(408240.*pow2(Mst1)*pow4(
        Msq))))/(Mgl*Tbeta*pow2(Sbeta)*pow5(Mst2));
}

/**
 * @return The off-diagonal (1, 2) = (2, 1) matrix element of the Higgs mass matrix as a double for the hierarchy 'H6g2'
 */
double H6g2::getS12() const {
   return (-(Al4p*MuSUSY*z2*(-68040*Mt*twoLoopFlag*pow2(Mst1)*pow4(Msq)*(Mgl*pow2(
        Mst2)*(2*Mt*MuSUSY*s2t*(4*Mt*((9*Dmglst2 + Mgl)*pow2(Mst1)*pow2(Mst2) +
        (13*Dmglst2 + Mgl)*pow4(Mst1) + (5*Dmglst2 + Mgl)*pow4(Mst2)) - Mst2*
        s2t*(Mgl*pow4(Mst2) + 4*Dmglst2*(pow2(Mst1)*pow2(Mst2) + pow4(Mst1) +
        pow4(Mst2)))) + Tbeta*(16*pow3(Mt)*(-2*Mgl*pow2(Mst1)*pow2(Mst2) + 3*(
        Dmglst2 - Mgl)*pow4(Mst1) - (Dmglst2 + Mgl)*pow4(Mst2)) - 6*Mt*pow2(
        Mst2)*pow2(s2t)*(4*Dmglst2*pow2(Mst1)*(pow2(Mst1) + pow2(Mst2)) + (5*
        Dmglst2 + Mgl)*pow4(Mst2)) + (4*Dmglst2*pow2(Mst2) + Mgl*(-pow2(Mst1) +
        pow2(Mst2)))*pow3(s2t)*pow5(Mst2))) + 4*Mt*xMst*(xDmglst2*pow2(Dmglst2)
        *(82*Mt*MuSUSY*s2t - 56*Tbeta*pow2(Mt) - 3*Mst2*(MuSUSY + 5*Mst2*Tbeta)
        *pow2(s2t)) + 2*Mgl*(Mgl*Mt*(MuSUSY*s2t - 8*Mt*Tbeta) + Dmglst2*(17*Mt*
        MuSUSY*s2t + 16*Tbeta*pow2(Mt) - Mst2*(MuSUSY + 3*Mst2*Tbeta)*pow2(s2t)
        )))*pow6(Mst1)) + pow2(Mst2)*(9*Mt*xDmglst2*pow2(Dmglst2)*(-15120*
        twoLoopFlag*pow2(Mst1)*pow4(Msq)*(-6*Mst2*Mt*MuSUSY*pow2(s2t)*(pow2(
        Mst1)*pow2(Mst2) + pow4(Mst1) + pow4(Mst2)) - 8*Tbeta*pow3(Mt)*(pow2(
        Mst1)*pow2(Mst2) + 4*pow4(Mst1) + pow4(Mst2)) - 3*Mt*Tbeta*pow2(Mst2)*
        pow2(s2t)*(10*pow2(Mst1)*(pow2(Mst1) + pow2(Mst2)) + 11*pow4(Mst2)) +
        4*MuSUSY*s2t*pow2(Mt)*(21*pow2(Mst1)*pow2(Mst2) + 31*pow4(Mst1) + 11*
        pow4(Mst2)) + 3*Tbeta*pow3(s2t)*pow7(Mst2)) + (Al4p*Mst2*threeLoopFlag*
        (113400*pow2(Msq)*(-(pow3(Mst2)*(-2304*MuSUSY*s2t*pow2(Mt) + 4*Mst2*Mt*
        (41*MuSUSY + 198*Mst2*Tbeta)*pow2(s2t) + 416*Tbeta*pow3(Mt) + 9*Tbeta*
        pow3(Mst2)*pow3(s2t))*pow4(Mst1)) + 13*pow2(Mst1)*(96*MuSUSY*s2t*pow2(
        Mt) - 2*Mst2*Mt*(7*MuSUSY + 36*Mst2*Tbeta)*pow2(s2t) - 32*Tbeta*pow3(
        Mt) + 7*Tbeta*pow3(Mst2)*pow3(s2t))*pow5(Mst2) - 8*Mst2*Mt*(-120*Mt*
        MuSUSY*s2t + 28*Tbeta*pow2(Mt) + Mst2*(7*MuSUSY + 27*Mst2*Tbeta)*pow2(
        s2t))*pow6(Mst1)) + pow4(Msq)*(-3*pow4(Mst1)*(8*Mst2*s2t*(7*(3255761 -
        1458000*lmMsq - 351000*lmMst1 + 2631960*lmMst2)*MuSUSY + 36*(-371369 +
        7875*lmMst1 - 320355*lmMst2)*Mst2*Tbeta)*pow2(Mt) - 7*Mt*(40*(47953 -
        14580*lmMsq - 25650*lmMst1 + 54162*lmMst2)*MuSUSY + 3*(2489911 -
        1425600*lmMsq - 716040*lmMst1 + 2927880*lmMst2)*Mst2*Tbeta)*pow2(Mst2)*
        pow2(s2t) + 32*(315*(3643 - 120*lmMst1 + 3192*lmMst2)*MuSUSY + 2*(-
        1161746 + 56700*lmMsq + 6615*lmMst1 - 222075*lmMst2 + 34020*lmMt)*Mst2*
        Tbeta)*pow3(Mt) + 14*(131737 - 96660*lmMst1 + 96660*lmMst2)*Tbeta*pow3(
        s2t)*pow4(Mst2)) - 9*pow2(Mst1)*pow2(Mst2)*(4*(Mst2*s2t*(7*(1340537 -
        496800*lmMsq + 4680*lmMst1 + 778680*lmMst2)*MuSUSY + 18*(-74087 + 1820*
        lmMst1 - 230300*lmMst2)*Mst2*Tbeta) + Mt*(1344*(617 + 2040*lmMst2)*
        MuSUSY + (718271 + 302400*lmMsq + 32760*lmMst1 - 955080*lmMst2 - 45360*
        lmMt)*Mst2*Tbeta))*pow2(Mt) - 7*Mt*(4*(115931 - 48600*lmMsq - 53280*
        lmMst1 + 148320*lmMst2)*MuSUSY + 3*(1340537 - 496800*lmMsq + 4680*
        lmMst1 + 778680*lmMst2)*Mst2*Tbeta)*pow2(Mst2)*pow2(s2t) + 14*(90011 -
        48600*lmMsq - 53280*lmMst1 + 148320*lmMst2)*Tbeta*pow3(s2t)*pow4(Mst2))
        + s2t*(14*(8*(2286439 - 72900*lmMsq - 307800*lmMst1 + 450360*lmMst2)*
        Mt*MuSUSY*s2t + 15*(1991719 - 194400*lmMsq - 328536*lmMst1 + 618840*
        lmMst2)*Mst2*Mt*s2t*Tbeta + 72*(520939 + 44280*lmMst1 + 76680*lmMst2)*
        Tbeta*pow2(Mt) + 8*(-776573 + 80595*lmMst1 - 80595*lmMst2)*Tbeta*pow2(
        Mst2)*pow2(s2t))*pow6(Mst1) - 3265920*(-2*Mt*MuSUSY*s2t - 4*Tbeta*pow2(
        Mt) + Tbeta*pow2(Mst2)*pow2(s2t))*pow6(Mst2))) + 56700*(-((-2504*
        MuSUSY*s2t*pow2(Mt) + 2*Mst2*Mt*(97*MuSUSY + 402*Mst2*Tbeta)*pow2(s2t)
        + 928*Tbeta*pow3(Mt) + 25*Tbeta*pow3(Mst2)*pow3(s2t))*pow4(Mst1)*pow5(
        Mst2)) - 4*Mt*(-186*Mt*MuSUSY*s2t + 88*Tbeta*pow2(Mt) + Mst2*(11*MuSUSY
        + 39*Mst2*Tbeta)*pow2(s2t))*pow3(Mst2)*pow6(Mst1) + 2*pow2(Mst1)*(716*
        MuSUSY*s2t*pow2(Mt) - Mst2*Mt*(122*MuSUSY + 537*Mst2*Tbeta)*pow2(s2t) -
        464*Tbeta*pow3(Mt) + 61*Tbeta*pow3(Mst2)*pow3(s2t))*pow7(Mst2))))/9.) +
        7*Al4p*threeLoopFlag*(-97200*Mst2*Mt*s2t*xMsq*pow2(Mgl)*(Tbeta*pow2(
        Mst2)*(4*(-1 + shiftst2)*pow2(Mst1)*pow2(Mt) - (shiftst1 - shiftst2 + (
        lmMst1 - lmMst2)*(-2 + shiftst1 + shiftst2))*pow2(Mst1)*pow2(Mst2)*
        pow2(s2t) + (-1 + shiftst1)*pow2(Mst2)*(-4*pow2(Mt) + pow2(Mst2)*pow2(
        s2t)) - (-1 + shiftst2)*pow2(s2t)*pow4(Mst1)) + 2*Mt*MuSUSY*s2t*((1 -
        2*lmMst2*(-1 + shiftst1))*pow2(Mst1)*pow2(Mst2) + 2*lmMst1*(-1 +
        shiftst1)*pow2(Mst1)*pow2(Mst2) + shiftst2*pow2(Mst1)*(2*pow2(Mst1) +
        pow2(Mst2)) - 2*lmMst1*(1 - 2*shiftst1 + shiftst2)*pow4(Mst1) + 2*
        lmMst2*(1 - 2*shiftst1 + shiftst2)*pow4(Mst1) + pow4(Mst2) - shiftst1*(
        2*pow2(Mst1)*pow2(Mst2) + 2*pow4(Mst1) + pow4(Mst2))))*pow6(Msq) -
        3240*Mt*xDR2DRMOD*pow4(Msq)*(Mst2*xDmglst2*pow2(Dmglst2)*(-(Mst2*Tbeta*
        (64*Mt*pow2(Mst1)*(4*(2*pow2(Mst2) - lmMst2*(pow2(Mst1) + pow2(Mst2)))*
        pow2(Mt) + 3*pow2(s2t)*((-1 + 6*lmMst2)*pow2(Mst2)*(pow2(Mst1) + pow2(
        Mst2)) + (1 + 3*lmMst2)*pow4(Mst1))) + pow3(Mst2)*(528*s2t*(pow2(Mst1)
        - pow2(Mst2))*pow2(Mt) - 3*pow3(s2t)*((8*(-8 + 33*lmMst1 - 21*lmMst2)*
        pow2(Mst1)*pow2(Mst2))/3. + 44*pow4(Mst1) - 44*pow4(Mst2))))) + 8*Mt*
        MuSUSY*s2t*((32*(-1 + 6*lmMst2)*Mt + (49 - 66*lmMst1 + 42*lmMst2)*Mst2*
        s2t)*pow2(Mst1)*pow3(Mst2) + 2*Mst2*(32*(-1 + 6*lmMst2)*Mt + (8 - 33*
        lmMst1 + 21*lmMst2)*Mst2*s2t)*pow4(Mst1) + (-46*lmMst1*s2t + 30*lmMst2*
        s2t)*pow6(Mst1) + 33*s2t*pow6(Mst2))) + Mgl*(-(Tbeta*(64*Mt*pow2(Mst1)*
        (3*(Dmglst2 + 3*Dmglst2*lmMst2 + Mgl + lmMst2*Mgl)*pow2(Mst2)*pow2(s2t)
        *(pow2(Mst1)*pow2(Mst2) + pow4(Mst1) + pow4(Mst2)) + 4*pow2(Mt)*(
        Dmglst2*(7 + lmMst2)*pow2(Mst1)*pow2(Mst2) - Dmglst2*(-1 + lmMst2)*
        pow4(Mst2) - (1 + lmMst2)*Mgl*(3*pow2(Mst1)*pow2(Mst2) + 6*pow4(Mst1) +
        pow4(Mst2)))) + (16*(23*Dmglst2 + 13*Mgl)*s2t*(pow2(Mst1) - pow2(Mst2))
        *pow2(Mt) - pow3(s2t)*(8*(Dmglst2*(23*lmMst1 - 15*lmMst2) + (4 + 13*
        lmMst1 - 9*lmMst2)*Mgl)*pow2(Mst1)*pow2(Mst2) + 4*(23*Dmglst2 + 13*Mgl)
        *(pow4(Mst1) - pow4(Mst2))))*pow5(Mst2))) + Mt*MuSUSY*s2t*(256*Mt*pow2(
        Mst1)*(Dmglst2*(1 + 3*lmMst2)*pow2(Mst2)*(2*pow2(Mst1) + pow2(Mst2)) +
        (1 + lmMst2)*Mgl*(2*pow2(Mst1)*pow2(Mst2) + 3*pow4(Mst1) + pow4(Mst2)))
        + 8*Mst2*s2t*(-2*(Dmglst2*(23*lmMst1 - 15*lmMst2) + (4 + 13*lmMst1 - 9*
        lmMst2)*Mgl)*(pow2(Mst1) + pow2(Mst2))*pow4(Mst1) + (Dmglst2*(23 - 46*
        lmMst1 + 30*lmMst2) + (5 - 26*lmMst1 + 18*lmMst2)*Mgl)*pow2(Mst1)*pow4(
        Mst2) + (23*Dmglst2 + 13*Mgl)*pow6(Mst2))))) + Mgl*(-3*Tbeta*(pow2(
        Mst1)*(4*(9*Dmglst2*pow2(Mst2)*((8*(57943 + 1800*lmMsq - 2850*lmMst1 +
        2190*lmMst2 - 1260*lmMt)*pow2(Mst1) + 15*(1805 + 960*lmMsq + 248*lmMst1
        - 2872*lmMst2 - 448*lmMt)*pow2(Mst2))*pow4(Msq) + 2400*(14*pow2(Msq) +
        11*pow2(Mst2))*(pow2(Mst1)*pow2(Mst2) + pow4(Mst1) + pow4(Mst2))) +
        Mgl*(64800*(2*pow2(Msq)*pow2(Mst2) + pow4(Mst2))*(pow2(Mst1)*pow2(Mst2)
        + pow4(Mst1) + pow4(Mst2)) - pow4(Msq)*(48*(45986 - 2700*lmMsq - 7245*
        lmMst1 + 21015*lmMst2 + 630*lmMt)*pow2(Mst1)*pow2(Mst2) + 4*(1126591 -
        184140*lmMst1 + 383940*lmMst2 + 7560*lmMt)*pow4(Mst1) + 9*(88373 -
        14400*lmMsq - 9960*lmMst1 + 52680*lmMst2 + 3360*lmMt)*pow4(Mst2))))*
        pow4(Mt) - 3*pow2(Mst2)*pow2(Mt)*pow2(s2t)*(2*(Dmglst2*(
        3.3195316666666665e6 - 324000*lmMsq - 547560*lmMst1 + 1031400*lmMst2) +
        (63043 - 64800*lmMsq + 217080*lmMst1 - 147960*lmMst2)*Mgl)*pow4(Msq)*
        pow4(Mst1) + 135*(-320*(9*Dmglst2 + Mgl)*pow2(Msq)*pow2(Mst1) + (
        Dmglst2*(18556.466666666667 - 5280*lmMsq + 248*lmMst1 + 7944*lmMst2) +
        (3469 - 1440*lmMsq + 664*lmMst1 + 1384*lmMst2)*Mgl)*pow4(Msq) - 80*(13*
        Dmglst2 + Mgl)*pow4(Mst1))*pow4(Mst2) + 15*pow2(Msq)*pow2(Mst2)*((
        Dmglst2*(260504.6 - 43200*lmMsq - 21240*lmMst1 + 85752*lmMst2) + (15395
        - 8640*lmMsq + 17928*lmMst1 - 8712*lmMst2)*Mgl)*pow2(Msq)*pow2(Mst1) -
        2880*(9*Dmglst2 + Mgl)*pow4(Mst1) - 5760*(6*Dmglst2 + Mgl)*pow4(Mst2))
        - 5400*pow3(Mst2)*(2*(13*Dmglst2 + Mgl)*pow2(Mst1)*pow3(Mst2) + (41*
        Dmglst2 + 5*Mgl)*pow5(Mst2)))) + Mt*pow3(Mst2)*pow3(s2t)*(180*pow2(Msq)
        *((Dmglst2*(3599 - 1080*lmMsq - 1368*lmMst1 + 3480*lmMst2) - (235 +
        450*lmMsq + 642*lmMst1 - 1362*lmMst2)*Mgl)*pow2(Msq)*pow2(Mst1) + 6*((
        92*Dmglst2 + 49*Mgl)*pow2(Msq) - 5*(34*Dmglst2 + 7*Mgl)*pow2(Mst1))*
        pow2(Mst2) + 90*(2*Dmglst2 + Mgl)*pow4(Mst1))*pow4(Mst2) + pow4(Msq)*(
        6*(6*Dmglst2*(33563 - 4020*lmMst1 + 4020*lmMst2) + (-242149 + 2700*
        lmMsq + 9360*lmMst1 - 19440*lmMst2)*Mgl)*pow2(Mst2)*pow4(Mst1) + ((16*
        Dmglst2*(776573 - 80595*lmMst1 + 80595*lmMst2) + 3*(-2205079 + 76140*
        lmMst1 - 76140*lmMst2)*Mgl)*pow6(Mst1))/3.) + 1350*pow2(Mst1)*(5*(4*
        Dmglst2 + Mgl)*pow2(Mst1) - (64*Dmglst2 + 9*Mgl)*pow2(Mst2))*pow6(Mst2)
        ) + pow4(Msq)*(3240*Mgl*shiftst3*(-(pow2(Mst2)*(4*s2t*pow3(Mt) + (1 +
        lmMst1 - lmMst2)*Mt*pow2(Mst1)*pow3(s2t))) + Mt*pow3(s2t)*pow4(Mst2))*
        pow5(Mst2) - 2*Mst2*s2t*pow3(Mt)*(-(Mgl*(3*(502399 + 65160*lmMst1 -
        134280*lmMst2)*pow2(Mst2)*pow4(Mst1) + 9*(-29683 + 13560*lmMst1 -
        36600*lmMst2)*pow2(Mst1)*pow4(Mst2) + (4828526 + 406080*lmMst1 -
        613440*lmMst2)*pow6(Mst1) - 105840*pow6(Mst2))) + 24*Dmglst2*(6*(21379
        + 1770*lmMst1 + 12630*lmMst2)*pow2(Mst2)*pow4(Mst1) + 6*(4879 - 270*
        lmMst1 + 8910*lmMst2)*pow2(Mst1)*pow4(Mst2) + (520939 + 44280*lmMst1 +
        76680*lmMst2)*pow6(Mst1) + 8280*pow6(Mst2))))) + 2*MuSUSY*pow2(Mt)*(2*
        Dmglst2*Mst2*(8100*Mst2*s2t*pow2(Mst1)*((186*Mt - 11*Mst2*s2t)*pow2(
        Mst2)*pow4(Mst1) + (134*Mt - 11*Mst2*s2t)*pow2(Mst1)*pow4(Mst2) + 2*
        pow2(Msq)*(-14*(-12*Mt + Mst2*s2t)*pow2(Mst1)*pow2(Mst2) + 2*(120*Mt -
        7*Mst2*s2t)*pow4(Mst1) + (96*Mt - 17*Mst2*s2t)*pow4(Mst2)) + 2*(41*Mt -
        8*Mst2*s2t)*pow6(Mst2)) + pow4(Msq)*(-27*pow2(Mst1)*pow2(Mst2)*((278347
        - 79200*lmMsq + 3720*lmMst1 + 119160*lmMst2)*Mst2*Mt*s2t + 7680*(5 + 6*
        lmMst2)*pow2(Mt) + 10*(-4151 + 1080*lmMsq + 1368*lmMst1 - 3480*lmMst2)*
        pow2(Mst2)*pow2(s2t)) - 36*((534391 - 113400*lmMsq - 23760*lmMst1 +
        196560*lmMst2)*Mst2*Mt*s2t + 30*(3643 - 120*lmMst1 + 3192*lmMst2)*pow2(
        Mt) + 9*(-9053 + 900*lmMsq + 1810*lmMst1 - 3570*lmMst2)*pow2(Mst2)*
        pow2(s2t))*pow4(Mst1) + pow2(s2t)*(4*(2286439 - 72900*lmMsq - 307800*
        lmMst1 + 450360*lmMst2)*pow6(Mst1) + 149040*pow6(Mst2)))) + 3*Mgl*(
        1350*s2t*pow2(Mst1)*pow2(Mst2)*((72*Mt - 4*Mst2*s2t)*pow2(Mst2)*pow4(
        Mst1) - 4*(-14*Mt + Mst2*s2t)*pow2(Mst1)*pow4(Mst2) + 4*pow2(Msq)*(-4*(
        -12*Mt + Mst2*s2t)*pow2(Mst1)*pow2(Mst2) + (64*Mt - 4*Mst2*s2t)*pow4(
        Mst1) + (32*Mt - 7*Mst2*s2t)*pow4(Mst2)) + (40*Mt - 9*Mst2*s2t)*pow6(
        Mst2)) + pow4(Msq)*(-90*pow2(Mst1)*(3*(3469 - 1440*lmMsq + 664*lmMst1 +
        1384*lmMst2)*Mst2*Mt*s2t + 64*(53 + 24*lmMst2)*pow2(Mt) - 2*(59 - 450*
        lmMsq + 1362*lmMst2 + 36*(1 + lmMst2)*shiftst3 - 6*lmMst1*(107 + 6*
        shiftst3))*pow2(Mst2)*pow2(s2t))*pow3(Mst2) + 6*Mst2*(-40*(5827 - 2700*
        lmMsq + 2988*lmMst1 + 468*lmMst2)*Mst2*Mt*s2t - 60*(-353 + 72*lmMst1 +
        696*lmMst2)*pow2(Mt) - (240379 + 10800*lmMsq - 21420*lmMst2 - 1080*(1 +
        2*lmMst2)*shiftst3 + 180*lmMst1*(55 + 12*shiftst3))*pow2(Mst2)*pow2(
        s2t))*pow4(Mst1) - s2t*(4*(412663 - 226800*lmMsq + 396360*lmMst1 -
        119880*lmMst2)*Mt + Mst2*s2t*(3647353 + 64800*lmMsq - 52380*lmMst2 -
        6480*(1 + 3*lmMst2)*shiftst3 + 540*lmMst1*(-31 + 36*shiftst3)))*pow6(
        Mst1) + 1080*(49 + 3*shiftst3)*pow2(s2t)*pow7(Mst2))))))))) + 1890*Mt*
        xMst*(27*(lmMst1 - lmMst2)*Mst2*Mt*oneLoopFlag*pow2(Mgl)*pow2(MuSUSY)*
        pow2(s2t) - Al4p*twoLoopFlag*(8*Mt*s2t*(18*(lmMst1 - lmMst2)*(-2 + 3*
        lmMst2)*Mst2*s2t*xDmglst2*xDR2DRMOD*pow2(Dmglst2) + Dmglst2*Mgl*(Mt*(
        785 + 6*lmMst1*(85 - 42*lmMst2) - 438*lmMst2 + 252*pow2(lmMst2)) -
        Mst2*s2t*(49 - 84*lmMst2 + lmMst1*(84 - 36*lmMst2*(-1 + xDR2DRMOD)) +
        36*(-1 + xDR2DRMOD)*pow2(lmMst2))) + (Mt*(193 + 474*lmMst2 - 6*lmMst1*(
        67 + 42*lmMst2) + 252*pow2(lmMst2)) - Mst2*s2t*(1 + 3*lmMst2*(-37 + 6*
        xDR2DRMOD) - 3*lmMst1*(-37 + 6*lmMst2*(-12 + xDR2DRMOD) + 6*xDR2DRMOD)
        - 81*pow2(lmMst1) + 9*(-15 + 2*xDR2DRMOD)*pow2(lmMst2)))*pow2(Mgl))*
        pow2(MuSUSY) + Mgl*MuSUSY*Tbeta*(8*Dmglst2*(36*(-1 + 3*lmMst1 - 3*
        lmMst2)*Mst2*s2t*pow2(Mt) - 3*Mt*(50 + lmMst1*(51 - 18*lmMst2) - 51*
        lmMst2 + 18*pow2(lmMst2))*pow2(Mst2)*pow2(s2t) + 4*(7 - 381*lmMst2 +
        lmMst1*(327 + 72*lmMst2 - 81*lmMt) + 54*lmMt + 81*lmMst2*lmMt - 72*
        pow2(lmMst2))*pow3(Mt) + 2*(1 + 3*lmMst1 - 3*lmMst2)*pow3(Mst2)*pow3(
        s2t)) - Mgl*(-144*Mst2*s2t*(-lmMst1 + lmMst2 - 2*lmMst1*lmMst2 + pow2(
        lmMst1) + pow2(lmMst2))*pow2(Mt) + 24*Mt*(1 + 33*lmMst2 - 3*lmMst1*(11
        + 6*lmMst2) + 18*pow2(lmMst2))*pow2(Mst2)*pow2(s2t) + 32*(83 + 3*
        lmMst1*(29 + 12*lmMst2 - 9*lmMt) + 9*lmMt + 3*lmMst2*(-32 + 9*lmMt) -
        36*pow2(lmMst2))*pow3(Mt) + (5 + lmMst1*(6 - 288*lmMst2) - 6*lmMst2 +
        144*(pow2(lmMst1) + pow2(lmMst2)))*pow3(Mst2)*pow3(s2t))) + 2*xDmglst2*
        pow2(Dmglst2)*(-3*Mst2*Mt*pow2(s2t)*(2*Mst2*MuSUSY*Tbeta*(355 + 6*
        lmMst2 - 6*lmMst1*(1 + 6*lmMst2) + 36*pow2(lmMst2)) + 2*(85 + 60*lmMst2
        + 12*lmMst1*(-5 + 3*lmMst2) - 36*pow2(lmMst2))*pow2(MuSUSY) + 5*(-43 +
        60*lmMst1 - 60*lmMst2)*pow2(Mst2)*(-1 + pow2(Sbeta))*pow2(Sbeta)) +
        MuSUSY*(2*s2t*pow2(Mt)*(MuSUSY*(5002 + 84*lmMst1*(7 - 6*lmMst2) - 444*
        lmMst2 + 504*pow2(lmMst2)) + 3*Mst2*Tbeta*(347 - 430*pow2(Sbeta) + 120*
        (lmMst1 - lmMst2)*(-4 + 5*pow2(Sbeta)))) + Tbeta*(8*(412 - 6*lmMst1*(
        176 + 42*lmMst2 - 39*lmMt) + 9*lmMst2*(147 - 26*lmMt) - 267*lmMt + 252*
        pow2(lmMst2))*pow3(Mt) + 12*(1 - 6*lmMst1 + 6*lmMst2)*pow3(Mst2)*pow3(
        s2t))))))*pow4(Msq)*pow8(Mst1))/(204120.*Tbeta*pow2(Mgl)*pow2(Mst1)*
        pow4(Msq)*pow7(Mst2)) + MuSUSY*((Mt*oneLoopFlag*s2t*(4*(lmMst1 -
        lmMst2)*pow2(Mt) - ((2 + lmMst1 - lmMst2)*pow2(Mst1) + (-2 + lmMst1 -
        lmMst2)*pow2(Mst2))*pow2(s2t) - (2*Mt*MuSUSY*s2t*(2 - lmMst1 + lmMst2 -
        (2*(lmMst1 - lmMst2)*pow2(Mst1)*(pow2(Mst1) + pow2(Mst2)))/pow4(Mst2)))
        /Tbeta))/16. - (Mt*threeLoopFlag*pow2(Al4p)*((24300*(1 - 2*lmMsq)*s2t*
        xMsq*pow2(Msq)*(Tbeta*pow2(Mst2)*((-1 + shiftst1)*pow2(Mst2)*(-4*pow2(
        Mt) + pow2(Mst2)*pow2(s2t)) + pow2(Mst1)*(4*(-1 + shiftst2)*pow2(Mt) +
        (-shiftst1 + shiftst2)*pow2(Mst2)*pow2(s2t)) - (-1 + shiftst2)*pow2(
        s2t)*pow4(Mst1)) + s2t*(Mt*MuSUSY*(2*(1 - 2*shiftst1 + shiftst2)*pow2(
        Mst1)*pow2(Mst2) - 4*(shiftst1 - shiftst2)*pow4(Mst1) - 2*(-1 +
        shiftst1)*pow4(Mst2)) + (-lmMst1 + lmMst2)*pow2(Mst1)*(4*Mt*MuSUSY*(1 -
        2*shiftst1 + shiftst2)*pow2(Mst1) - 4*Mt*MuSUSY*(-1 + shiftst1)*pow2(
        Mst2) + s2t*(-2 + shiftst1 + shiftst2)*Tbeta*pow4(Mst2)))))/(pow2(Mst1)
        *pow4(Mst2)) - (5*(Dmglst2*Mgl*Mst2*(-2*s2t*(72*Tbeta*pow2(Mt)*(2806*z3
        + 1514*z4 + 2757*pow2(z2)) - 4*Tbeta*pow2(Mst2)*pow2(s2t)*(89947*z3 +
        6332*z4 + 9498*pow2(z2)) + Mt*s2t*(8*MuSUSY*(338536*z3 + 11327*z4 +
        15897*pow2(z2)) + 3*Mst2*Tbeta*(-844756*z3 + 10126*z4 + 105585*pow2(z2)
        )))*pow4(Mst1) + 9*pow2(Mst1)*(-16*Mst2*s2t*pow2(Mt)*(MuSUSY*(59336*z3
        + 709*z4 - 8535*pow2(z2)) + 6*Mst2*Tbeta*(590*z3 - 4*z4 + 75*pow2(z2)))
        - 3*Mt*pow2(Mst2)*pow2(s2t)*(16*MuSUSY*(9207*z3 + 185*z4 + 237*pow2(z2)
        ) + Mst2*Tbeta*(-124922*z3 - 1210*z4 + 18273*pow2(z2))) + 96*(MuSUSY*(
        250*z3 - 94*z4 + 21*pow2(z2)) + 2*Mst2*Tbeta*(-955*z3 + 26*z4 + 93*
        pow2(z2)))*pow3(Mt) + 6*Tbeta*(20653*z3 + 316*z4 + 474*pow2(z2))*pow3(
        s2t)*pow4(Mst2)) + 27*pow2(Mst2)*(-4*Mt*s2t*pow2(Mst2)*(4*Mt*Tbeta*(
        439*z3 - 108*z4) + MuSUSY*s2t*(16175*z3 + 424*z4 + 474*pow2(z2))) + 4*
        Mst2*pow2(Mt)*(Mt*Tbeta*(5742*z3 - 506*z4 + 105*pow2(z2)) + MuSUSY*s2t*
        (-37474*z3 - 542*z4 + 5289*pow2(z2))) - 3*Mt*Tbeta*pow2(s2t)*(-37474*z3
        - 542*z4 + 5289*pow2(z2))*pow3(Mst2) + 64*MuSUSY*(83*z3 - 27*z4)*pow3(
        Mt) + 2*Tbeta*(16175*z3 + 424*z4 + 474*pow2(z2))*pow3(s2t)*pow4(Mst2)))
        - 3*pow2(Mgl)*((4*s2t*pow2(Mt)*(MuSUSY*(134464*z3 + 23528*z4 - 41010*
        pow2(z2)) + Mst2*Tbeta*(22466*z3 - 11402*z4 - 17103*pow2(z2))) + 2*
        Mst2*Mt*pow2(s2t)*(MuSUSY*(8*(76813 + 162*lmMst1 - 162*lmMst2)*z3 -
        15686*z4 - 18183*pow2(z2)) + 3*Mst2*Tbeta*(1144*z3 - 4954*z4 + 6177*
        pow2(z2))) + 16*Tbeta*(-107072*z3 + 3326*z4 + 3045*pow2(z2))*pow3(Mt) +
        Tbeta*(-272636*z3 + 8870*z4 + 13305*pow2(z2))*pow3(Mst2)*pow3(s2t))*
        pow4(Mst1) - 9*pow3(Mst2)*(-4*Mst2*pow2(Mt)*(5*MuSUSY*s2t*(2098*z3 +
        102*z4 - 333*pow2(z2)) + 3*Mt*Tbeta*(-1922*z3 + 206*z4 + 21*pow2(z2)))
        + 6*Mt*s2t*pow2(Mst2)*(-2*MuSUSY*s2t*((1913 + 12*lmMst1 - 12*lmMst2)*z3
        - 10*z4 + 48*pow2(z2)) + Mt*Tbeta*((1622 - 48*lmMst1 + 48*lmMst2)*z3 +
        142*z4 + 69*pow2(z2))) + 15*Mt*Tbeta*pow2(s2t)*(2098*z3 + 102*z4 - 333*
        pow2(z2))*pow3(Mst2) - 32*MuSUSY*(230*z3 + 27*z4)*pow3(Mt) + 6*Tbeta*((
        1913 + 12*lmMst1 - 12*lmMst2)*z3 - 10*z4 + 48*pow2(z2))*pow3(s2t)*pow4(
        Mst2)) + 3*Mst2*pow2(Mst1)*(-2*Mst2*s2t*pow2(Mt)*(8*MuSUSY*(-11396*z3 -
        1135*z4 + 2388*pow2(z2)) + Mst2*Tbeta*(7078*z3 + 2990*z4 + 4485*pow2(
        z2))) + Mt*pow2(Mst2)*pow2(s2t)*(-4*MuSUSY*((-56978 - 216*lmMst1 + 216*
        lmMst2)*z3 + 1136*z4 + 813*pow2(z2)) + 3*Mst2*Tbeta*(-14114*z3 - 3010*
        z4 + 4557*pow2(z2))) + 16*(9*MuSUSY*(218*z3 + 50*z4 + 21*pow2(z2)) + 2*
        Mst2*Tbeta*(-8219*z3 + 326*z4 + 165*pow2(z2)))*pow3(Mt) + 2*Tbeta*((-
        39761 - 108*lmMst1 + 108*lmMst2)*z3 + 1046*z4 + 1245*pow2(z2))*pow3(
        s2t)*pow4(Mst2))) + Mst2*xDmglst2*pow2(Dmglst2)*(-2*s2t*(72*Tbeta*pow2(
        Mt)*(2806*z3 + 1514*z4 + 2757*pow2(z2)) - 4*Tbeta*pow2(Mst2)*pow2(s2t)*
        (89947*z3 + 6332*z4 + 9498*pow2(z2)) + Mt*s2t*(8*MuSUSY*(338536*z3 +
        11327*z4 + 15897*pow2(z2)) + 3*Mst2*Tbeta*(-844756*z3 + 10126*z4 +
        105585*pow2(z2))))*pow4(Mst1) + 3*pow2(Mst1)*(-8*Mst2*s2t*pow2(Mt)*(
        MuSUSY*(699292*z3 + 12344*z4 - 88647*pow2(z2)) + 36*Mst2*Tbeta*(3124*z3
        + 29*z4 + 165*pow2(z2))) - Mt*pow2(Mst2)*pow2(s2t)*(8*MuSUSY*(409448*z3
        + 10819*z4 + 3471*pow2(z2)) + 3*Mst2*Tbeta*(-629162*z3 - 17194*z4 +
        84045*pow2(z2))) + 8*(Mst2*Tbeta*(118001*z3 - 6920*z4 - 6492*pow2(z2))
        + 36*MuSUSY*(250*z3 - 94*z4 + 21*pow2(z2)))*pow3(Mt) + Tbeta*(958451*z3
        + 16612*z4 + 1590*pow2(z2))*pow3(s2t)*pow4(Mst2)) + 9*pow2(Mst2)*(-4*
        Mst2*pow2(Mt)*(MuSUSY*s2t*(256474*z3 + 2498*z4 - 31083*pow2(z2)) + Mt*
        Tbeta*(-55952*z3 + 1742*z4 + 21*pow2(z2))) - 2*Mt*s2t*pow2(Mst2)*(9*Mt*
        Tbeta*(19703*z3 - 376*z4 + 84*pow2(z2)) + MuSUSY*s2t*(226447*z3 + 8888*
        z4 + 4098*pow2(z2))) + 3*Mt*Tbeta*pow2(s2t)*(256474*z3 + 2498*z4 -
        31083*pow2(z2))*pow3(Mst2) + 96*MuSUSY*(1990*z3 - 81*z4)*pow3(Mt) +
        Tbeta*(226447*z3 + 8888*z4 + 4098*pow2(z2))*pow3(s2t)*pow4(Mst2)))))/(
        pow2(Mgl)*pow5(Mst2)) + Mt*MuSUSY*((20*(9*Mst2*pow2(Mt)*(3*Mgl*(8*pow2(
        Mst2)*(436 - 6*B4 + 6*D3 - 3*DN - 48*lmMst1 + (408 - 96*lmMst1)*lmMst2
        + 24*lmMt - 972*S2 - 48*(-4 + lmMst1)*pow2(lmMst2) + 48*pow3(lmMst2)) +
        pow2(Mst1)*(10667 - 96*B4 + 96*D3 - 48*DN - 3072*lmMst1 - 384*(-25 + 9*
        lmMst1)*lmMst2 - 384*lmMst1*lmMt + 384*(1 + lmMst2)*lmMt - 224*OepS2 +
        324*(-43 + 14*lmMst1 - 14*lmMst2)*S2 - 384*(-13 + 4*lmMst1)*pow2(
        lmMst2) + 1536*pow3(lmMst2))) + 2*Dmglst2*(pow2(Mst1)*(28405 - 288*B4 +
        288*D3 - 144*DN + 5472*lmMst1 - 288*(-89 + 2*lmMst1)*lmMst2 + 576*(-2 +
        lmMst1 - lmMst2)*lmMt + 224*OepS2 - 324*(65 + 14*lmMst1 - 14*lmMst2)*S2
        - 576*(-9 + 8*lmMst1)*pow2(lmMst2) + 4608*pow3(lmMst2)) + 72*pow2(Mst2)
        *(180 - 2*B4 + 2*D3 - DN + 144*lmMst2 - 216*S2 - 16*(-2 + lmMst1)*pow2(
        lmMst2) + 16*(lmMst1 + pow3(lmMst2))))) + 4*Dmglst2*Mst2*T1ep*(-36*
        pow2(Mst1)*(253*Mst2*Mt*s2t + 42*pow2(Mt) - 129*pow2(Mst2)*pow2(s2t)) +
        189*s2t*(-5*Mt + 2*Mst2*s2t)*pow3(Mst2) + 17308*pow2(s2t)*pow4(Mst1)) +
        6*Mgl*T1ep*(6*Mst2*pow2(Mst1)*(272*Mst2*Mt*s2t + 252*pow2(Mt) - 487*
        pow2(Mst2)*pow2(s2t)) + s2t*(3764*Mt - 7357*Mst2*s2t)*pow4(Mst1) + 54*
        s2t*(7*Mt - 8*Mst2*s2t)*pow4(Mst2))))/(Mgl*pow5(Mst2)) + 90*pow2(s2t)*(
        15829 - 720*B4 + 72*D3 - 36*DN - 3330*lmMsq + 1350*pow2(lmMsq) - 3*
        lmMst1*(5279 - 1950*lmMsq + 630*pow2(lmMsq)) - 954*pow2(lmMst1) + 3*
        lmMst2*(10421 - 6558*lmMst1 + 30*lmMsq*(-95 + 42*lmMst1) + 630*pow2(
        lmMsq) + 234*pow2(lmMst1)) - 18*(-1460 + 210*lmMsq + 399*lmMst1)*pow2(
        lmMst2) - 288*pow3(lmMst1) + 6768*pow3(lmMst2) + (6*pow2(Mst1)*((Mgl*(
        111 - 3028*lmMst2 + 90*lmMsq*(4 - 13*lmMst1 + 13*lmMst2) + 2*lmMst1*(
        1334 + 675*lmMst2) - 90*pow2(lmMst1) - 1260*pow2(lmMst2)) + 100*
        Dmglst2*(69 + (-53 + 42*lmMsq)*lmMst2 + lmMst1*(53 - 42*lmMsq + 42*
        lmMst2) - 42*pow2(lmMst2)))/(Mgl*pow2(Msq)) - (135*(204.20053497942388
         + (76*B4)/
        9. - (2*DN)/9. - (50*lmMsq)/3. + (2*lmMst1*(57971 - 14625*lmMsq + 2700*
        pow2(lmMsq)))/405. - (2*lmMst2*(52436 - 70455*lmMst1 + 225*lmMsq*(-65 +
        36*lmMst1) + 2700*pow2(lmMsq) + 8280*pow2(lmMst1)))/405. + ((-1331 +
        180*lmMsq)*pow2(lmMst1) + (-8063 + 900*lmMsq)*pow2(lmMst2) + 3792*
        lmMst1*pow2(lmMst2))/27. - (62*pow3(lmMst1))/27. - (2626*pow3(lmMst2))/
        27. - (Dmglst2*(109.11799176954733 - (8*B4)/3. + (32*D3)/9. - (20*DN)/
        9. + 80*lmMsq - lmMst1*(78.19061728395062 + 20*pow2(lmMsq)) - (2888*
        pow2(lmMst1))/135. + lmMst2*(40*lmMsq*lmMst1 + 20*pow2(lmMsq) + (4*(-
        21616 - 64515*lmMst1 + 31275*pow2(lmMst1)))/2025.) - (4*(-5023 + 1350*
        lmMsq + 6285*lmMst1)*pow2(lmMst2))/135. + (20*pow3(lmMst1))/27. + (
        3340*pow3(lmMst2))/27.))/Mgl))/pow2(Mst2)))/5. - 162*S2*(1056 + 48*
        lmMst1 - 48*lmMst2 + ((9473 + 974*lmMst1 - 974*lmMst2 + (3*Dmglst2*(
        6952.444444444444 - 344*lmMst1 + 344*lmMst2))/Mgl)*pow2(Mst1))/(3.*
        pow2(Mst2)) + (Dmglst2*(3354 - 28*lmMst1 + 28*lmMst2) + ((8*Dmglst2*(
        93919 - 12981*lmMst1 + 12981*lmMst2) + 3*(160997 + 22071*lmMst1 -
        22071*lmMst2)*Mgl)*pow4(Mst1))/(81.*pow4(Mst2)))/Mgl) + 162*pow4(Mst1)*
        ((1.0702990137854083 + (9571*lmMsq)/26460. + ((56221 - 35490*lmMsq)*
        lmMst1)/13230. + (-4.6112244897959185 + (19*lmMsq)/7. + (31*lmMst1)/9.)
        *lmMst2 - pow2(lmMsq)/63. - (8*pow2(lmMst1))/21. + (5*Dmglst2*(216 + (-
        95 + 132*lmMsq)*lmMst2 + lmMst1*(95 - 132*lmMsq + 132*lmMst2) - 132*
        pow2(lmMst2)))/(54.*Mgl) - (194*pow2(lmMst2))/63.)/pow4(Msq) - (
        363.3804294212688 + (76*B4)/9. - (2*DN)/9. - (35*lmMsq)/2. + lmMst1*(
        211.3489518770471 - (695*lmMsq)/9. + (40*pow2(lmMsq))/3.) - (
        214.87936507936507 - 20*lmMsq)*pow2(lmMst1) - lmMst2*(
        190.46006298815823 - (71398*lmMst1)/105. + (5*lmMsq*(-139 + 120*lmMst1)
        )/9. + (40*pow2(lmMsq))/3. + (334*pow2(lmMst1))/3.) + ((-146507 +
        14700*lmMsq + 91070*lmMst1)*pow2(lmMst2))/315. - (44*pow3(lmMst1))/9. -
        (1556*pow3(lmMst2))/9. - (Dmglst2*(536.1152102791342 - (8*B4)/3. + (32*
        D3)/9. - (20*DN)/9. + 90*lmMsq - (123.11224321827497 + 20*lmMsq*(1 +
        lmMsq))*lmMst1 - lmMst2*(17.33220122616948 - 20*lmMsq*(1 + lmMsq) + (
        133.04550264550264 - 40*lmMsq)*lmMst1 - (1180*pow2(lmMst1))/9.) - (
        15886*pow2(lmMst1))/945. + (149.85608465608465 - 40*lmMsq - (2812*
        lmMst1)/9.)*pow2(lmMst2) - (92*pow3(lmMst1))/27. + (4988*pow3(lmMst2))/
        27.))/Mgl)/pow4(Mst2)) + (162*((Dmglst2*(3197 - 432*B4 + 576*D3 - 360*
        DN + 9180*lmMsq + 216*lmMst1*(7 - 15*pow2(lmMsq)) + 1620*pow2(lmMsq) -
        1404*pow2(lmMst1) + 36*lmMst2*(-291 - 464*lmMst1 + 90*lmMsq*(-1 + 2*
        lmMst1) + 90*pow2(lmMsq) + 32*pow2(lmMst1)) - 36*(-607 + 180*lmMsq +
        336*lmMst1)*pow2(lmMst2) + 10944*pow3(lmMst2)))/162. - (16*(1 + lmMst2)
        *(4*Dmglst2*lmMst2 + Mgl + lmMst2*Mgl)*pow4(Mst2))/(9.*pow4(Mst1)) - (
        4*OepS2*(4*Dmglst2*(2322*pow2(Mst1)*pow2(Mst2) + 8654*pow4(Mst1) + 189*
        pow4(Mst2)) - 3*Mgl*(2922*pow2(Mst1)*pow2(Mst2) + 7357*pow4(Mst1) +
        432*pow4(Mst2))))/(2187.*pow4(Mst2)) + (pow2(Mst2)*((Mgl*(2*(4167613 -
        19932360*lmMst2 + 20580*lmMst1*(701 + 540*lmMst2) + 420*lmMsq*(13109 -
        26460*lmMst1 + 25620*lmMst2) + 176400*pow2(lmMsq) - 10936800*pow2(
        lmMst2))*pow2(Mst1) + (41220947 - 420*lmMsq*(12479 + 69090*lmMst1 -
        69930*lmMst2) - 21234990*lmMst2 + 10290*lmMst1*(2573 + 2820*lmMst2) -
        176400*pow2(lmMsq) - 29194200*pow2(lmMst2))*pow2(Mst2)))/1.11132e7 + (
        5*Dmglst2*(2*(219 + (-95 + 132*lmMsq)*lmMst2 + lmMst1*(95 - 132*lmMsq +
        132*lmMst2) - 132*pow2(lmMst2))*pow2(Mst1) + (557 - 224*lmMst2 + 4*
        lmMst1*(53 + 96*lmMst2) + lmMsq*(12 - 384*lmMst1 + 384*lmMst2) - 384*
        pow2(lmMst2))*pow2(Mst2)))/108. - (40*(8*Dmglst2*(14 - 15*lmMsq + 15*
        lmMst2) + (43 - 30*lmMsq + 30*lmMst2)*Mgl)*pow2(Msq)*pow2(Mst2) + (-12*
        Mgl*(341 + 642*lmMst2 + 64*lmMst1*(1 + lmMst2) - 30*lmMsq*(5 + 6*
        lmMst2) + 90*pow2(lmMsq) + 272*pow2(lmMst2)) - 24*Dmglst2*(-115 + (366
        + 64*lmMst1)*lmMst2 - 30*lmMsq*(-1 + 6*lmMst2) + 90*pow2(lmMsq) + 336*
        pow2(lmMst2)))*pow4(Msq) + (90*Dmglst2*(13 - 28*lmMsq + 28*lmMst2) + 5*
        (67 - 84*lmMsq + 84*lmMst2)*Mgl)*pow4(Mst2))/(216.*pow2(Mst1))))/pow4(
        Msq) + (Mgl*((1725 + (-7006 + 2640*lmMsq)*lmMst2 + lmMst1*(7006 - 2640*
        lmMsq + 4800*lmMst2) - 1080*pow2(lmMst1) - 3720*pow2(lmMst2))*pow4(
        Mst1) + 2*(836 - 2235*lmMst2 + 75*lmMst1*(27 + 16*lmMst2) + 30*lmMsq*(7
        - 40*lmMst1 + 40*lmMst2) - 1200*pow2(lmMst2))*pow4(Mst2)) + 50*Dmglst2*
        ((291 + 2*(-103 + 84*lmMsq)*lmMst2 + 2*lmMst1*(103 - 84*lmMsq + 84*
        lmMst2) - 168*pow2(lmMst2))*pow4(Mst1) + 2*(118 + 109*lmMst1 + (-133 +
        102*lmMst1)*lmMst2 + 6*lmMsq*(4 - 17*lmMst1 + 17*lmMst2) - 102*pow2(
        lmMst2))*pow4(Mst2)))/(270.*pow2(Msq)*pow2(Mst2))))/Mgl + (54*(-1 + 2*
        lmMst2)*shiftst3*(2*(1 - 2*lmMst1 + 2*lmMst2)*pow2(Mst2)*pow4(Mst1) +
        2*(1 - lmMst1 + lmMst2)*pow2(Mst1)*pow4(Mst2) + (2 - 6*lmMst1 + 6*
        lmMst2)*pow6(Mst1) + pow6(Mst2)))/(pow2(Mst1)*pow4(Mst2))) + (Mt*s2t*(
        14580*pow2(Mst1)*pow2(Mst2)*(1035.3004115226338 + (1016*B4)/9. - (4*DN)
        /9. - 240*lmMsq + (80*pow2(lmMsq))/3. - (8*lmMst1*(422 - 45*lmMsq + 45*
        pow2(lmMsq)))/9. - (176*pow2(lmMst1))/9. + (8*lmMst2*(1096 + 15*lmMsq*(
        -7 + 6*lmMst1 - 6*lmMst2) + 717*lmMst2 - lmMst1*(551 + 248*lmMst2) +
        45*pow2(lmMsq) + 8*pow2(lmMst1)))/9. + (640*pow3(lmMst2))/3. + ((-80*(
        Mgl*(55 + 6*lmMsq*(-2 + 5*lmMst1 - 5*lmMst2) + 52*lmMst2 - 10*lmMst1*(4
        + 3*lmMst2) + 30*pow2(lmMst2)) + Dmglst2*(321 + 18*lmMsq*(-2 + 5*lmMst1
        - 5*lmMst2) + 96*lmMst2 - 30*lmMst1*(2 + 3*lmMst2) + 90*pow2(lmMst2)))*
        pow2(Mst1))/(27.*pow2(Msq)) + Dmglst2*(881.6139917695473 + 248*B4 - 4*
        DN - (2560*lmMsq)/3. + lmMst1*(130.96296296296296 - 120*lmMsq - 40*
        pow2(lmMsq)) + (80*pow2(lmMsq))/3. + (176*pow2(lmMst1))/9. + (8*lmMst2*
        (4364 - 573*lmMst1 + 45*lmMsq*(5 + 6*lmMst1) + 135*pow2(lmMsq) + 24*
        pow2(lmMst1)))/27. - (8*(-377 + 90*lmMsq + 376*lmMst1)*pow2(lmMst2))/9.
         + (2944*pow3(lmMst2))/9.))/Mgl) + 10*(2552929 + 257904*B4 - 648*DN -
        456840*lmMsq + 38880*pow2(lmMsq) - 216*lmMst1*(4591 - 360*lmMsq + 450*
        pow2(lmMsq)) + 41904*pow2(lmMst1) + 216*lmMst2*(9211 - 6466*lmMst1 +
        180*lmMsq*(-4 + 5*lmMst1) + 450*pow2(lmMsq) + 576*pow2(lmMst1)) - 864*(
        -1784 + 225*lmMsq + 840*lmMst1)*pow2(lmMst2) - 3456*pow3(lmMst1) +
        604800*pow3(lmMst2))*pow4(Mst1) + ((9*(4800*(-2*Mgl*(23 + lmMsq*(-6 +
        9*lmMst1 - 9*lmMst2) + 18*lmMst2 - 3*lmMst1*(4 + 3*lmMst2) + 9*pow2(
        lmMst2)) - 3*Dmglst2*(77 + 6*lmMsq*(-2 + 3*lmMst1 - 3*lmMst2) + 24*
        lmMst2 - 6*lmMst1*(2 + 3*lmMst2) + 18*pow2(lmMst2)))*pow2(Msq)*pow2(
        Mst1) + (45*Mgl*(18201 + 1760*B4 - 16*DN - 5760*lmMsq + 960*pow2(lmMsq)
        - 16*lmMst1*(199 - 30*lmMsq + 30*pow2(lmMsq)) + 16*lmMst2*(1291 - 322*
        lmMst1 + 30*lmMsq*(-5 + 2*lmMst1) + 30*pow2(lmMsq) - 16*pow2(lmMst1)) -
        256*pow2(lmMst1) - 32*(-313 + 30*lmMsq + 72*lmMst1)*pow2(lmMst2) +
        2560*pow3(lmMst2)) + Dmglst2*(779917 + 188640*B4 - 3600*DN - 648000*
        lmMsq + 43200*pow2(lmMsq) - 720*lmMst1*(-173 + 90*lmMsq + 30*pow2(
        lmMsq)) + 720*lmMst2*(1623 - 130*lmMst1 + 30*lmMsq*(-1 + 2*lmMst1) +
        30*pow2(lmMsq) - 16*pow2(lmMst1)) + 11520*pow2(lmMst1) - 1440*(-265 +
        30*lmMsq + 136*lmMst1)*pow2(lmMst2) + 207360*pow3(lmMst2)))*pow4(Msq) -
        300*(Mgl*(233 + 36*lmMsq*(-2 + 5*lmMst1 - 5*lmMst2) + 207*lmMst2 - 45*
        lmMst1*(3 + 4*lmMst2) + 180*pow2(lmMst2)) + Dmglst2*(1961 + 180*lmMsq*(
        -2 + 5*lmMst1 - 5*lmMst2) + 675*lmMst2 - 45*lmMst1*(7 + 20*lmMst2) +
        900*pow2(lmMst2)))*pow4(Mst1))*pow4(Mst2))/pow4(Msq) - 160*OepS2*((-
        3036*Dmglst2 + 816*Mgl)*pow2(Mst1)*pow2(Mst2) + 1882*Mgl*pow4(Mst1) -
        63*(5*Dmglst2 - 3*Mgl)*pow4(Mst2)) - 108*S2*(3*Dmglst2*pow2(Mst2)*(8*(
        2489 + 3795*lmMst1 - 3795*lmMst2)*pow2(Mst1) + 9*(-453 + 350*lmMst1 -
        350*lmMst2)*pow2(Mst2)) - 5*Mgl*(36*(169 + 136*lmMst1 - 136*lmMst2)*
        pow2(Mst1)*pow2(Mst2) + 2*(9185 + 5646*lmMst1 - 5646*lmMst2)*pow4(Mst1)
        + 81*(-1 + 14*lmMst1 - 14*lmMst2)*pow4(Mst2))) + (90*(-48*(15*Dmglst2*(
        95 + 12*lmMsq*(-2 + lmMst1 - lmMst2) + 32*lmMst2 - 4*lmMst1*(2 + 3*
        lmMst2) + 12*pow2(lmMst2)) + 5*Mgl*(71 + 12*lmMsq*(-2 + lmMst1 -
        lmMst2) + 40*lmMst2 - 4*lmMst1*(4 + 3*lmMst2) + 12*pow2(lmMst2)))*pow2(
        Msq) - 5*(Mgl*(1153 + 12*lmMsq*(-35 + 54*lmMst1 - 54*lmMst2) + 906*
        lmMst2 - 162*lmMst1*(3 + 4*lmMst2) + 648*pow2(lmMst2)) + Dmglst2*(8713
        + 12*lmMsq*(-179 + 270*lmMst1 - 270*lmMst2) + 3282*lmMst2 - 162*lmMst1*
        (7 + 20*lmMst2) + 3240*pow2(lmMst2)))*pow2(Mst1) - 5*(Mgl*(911 + 12*
        lmMsq*(-37 + 18*lmMst1 - 18*lmMst2) + 606*lmMst2 - 54*lmMst1*(3 + 4*
        lmMst2) + 216*pow2(lmMst2)) + Dmglst2*(5591 + 12*lmMsq*(-181 + 90*
        lmMst1 - 90*lmMst2) + 2550*lmMst2 - 54*lmMst1*(7 + 20*lmMst2) + 1080*
        pow2(lmMst2)))*pow2(Mst2) + (2304*(1 + lmMst2)*(-Dmglst2 + 3*Dmglst2*
        lmMst2 + Mgl + lmMst2*Mgl)*pow4(Msq))/pow2(Mst1))*pow6(Mst2))/pow4(Msq)
        )/Mgl))/pow5(Mst2))))/(14580.*Tbeta) - threeLoopFlag*pow2(Al4p)*(-(s2t*
        pow3(Mt)*(227.09043209876543 + (8*D3)/9. - (8*DN)/9. - (50*lmMsq)/3. +
        (32*lmMt)/3. + 10*pow2(lmMsq) + (2*lmMst1*(667 - 300*lmMsq + 45*pow2(
        lmMsq)))/27. + (2*lmMst2*(1946 + lmMsq*(30 - 90*lmMst1) + 294*lmMst1 -
        45*pow2(lmMsq) - 117*pow2(lmMst1)))/27. + (124*pow2(lmMst1))/9. + (2*(
        152 + 30*lmMsq + 31*lmMst1)*pow2(lmMst2))/9. + (32*pow3(lmMst1))/9. - (
        16*pow3(lmMst2))/9. + (4*Dmglst2*(15529 - 120*B4 + 120*D3 - 60*DN +
        150*lmMsq + 690*lmMst1 + 675*pow2(lmMsq) + 345*pow2(lmMst1) - 5*lmMst2*
        (-3091 + 270*lmMsq + 42*lmMst1 + 96*pow2(lmMst1)) + (4305 - 480*lmMst1)
        *pow2(lmMst2) + 960*pow3(lmMst2)))/(135.*Mgl) - S2*(267.3 - 69*lmMst1 +
        69*lmMst2 - ((370.43333333333334 + (1495*lmMst1)/3. - (1495*lmMst2)/3.
         - (4*Dmglst2*(139 + 100*lmMst1 - 100*lmMst2))/Mgl)*pow2(Mst1))/pow2(
        Mst2) + ((2532*Dmglst2)/5. + (2*(12*Dmglst2*(17269 + 13785*lmMst1 -
        13785*lmMst2) + (-126661 - 85515*lmMst1 + 85515*lmMst2)*Mgl)*pow4(Mst1)
        )/(135.*pow4(Mst2)))/Mgl) - pow4(Mst1)*((1.670354173415398 - (lmMsq*(
        256 + 259*lmMst1 - 49*lmMst2))/441. - (106*lmMst2)/135. + lmMst1*(
        1.365684051398337 + (19*lmMst2)/9.) + (10*Dmglst2)/(9.*Mgl) + (5*pow2(
        lmMsq))/21. - (16*pow2(lmMst1))/21. - (10*pow2(lmMst2))/9.)/pow4(Msq) -
        (594.0859288013853 - (40*lmMsq)/3. + (-480.17666414714034 + 20*lmMsq)*
        lmMst1 + lmMst2*(707.7322197026959 - (69154*lmMst1)/945. - (20*lmMsq*(3
        + 2*lmMst1))/3. - (410*pow2(lmMst1))/9.) + ((-94783 + 6300*lmMsq)*pow2(
        lmMst1))/945. + (216.14497354497354 + (20*lmMsq)/3. - (166*lmMst1)/9.)*
        pow2(lmMst2) + (32*lmMt*(1 + 4*lmMst2 - 2*lmMst1*(2 + lmMst2) + pow2(
        lmMst1) + pow2(lmMst2)))/3. + (26*pow3(lmMst1))/27. + (1702*pow3(
        lmMst2))/27. + (Dmglst2*(386.22610253272387 - (32*B4)/9. + (32*D3)/9. -
        (16*DN)/9. + (280*lmMsq)/3. + (756.6053917863442 - 80*lmMsq)*lmMst1 + (
        19352*pow2(lmMst1))/189. + lmMst2*(196.28349710254471 + 80*lmMsq - (
        57520*lmMst1)/189. + (88*pow2(lmMst1))/9.) + (8*(6787 - 5439*lmMst1)*
        pow2(lmMst2))/189. - (64*lmMt*(4 + 6*lmMst2 - 2*lmMst1*(3 + lmMst2) +
        pow2(lmMst1) + pow2(lmMst2)))/3. - (8*pow3(lmMst1))/9. + (664*pow3(
        lmMst2))/3.))/Mgl)/pow4(Mst2)) + ((32*(1 + lmMst2)*(4*Dmglst2*lmMst2 +
        Mgl + lmMst2*Mgl)*pow4(Mst2))/(9.*pow4(Mst1)) + (4*OepS2*(24*Dmglst2*(
        150*pow2(Mst1)*pow2(Mst2) + 919*pow4(Mst1)) - Mgl*(4485*pow2(Mst1)*
        pow2(Mst2) + 11402*pow4(Mst1) + 621*pow4(Mst2))))/(729.*pow4(Mst2)) + (
        pow2(Mst2)*((-5*Dmglst2*(3*(7 + 12*lmMsq - 12*lmMst2)*pow2(Mst1) + 2*(
        43 - 144*lmMsq + 144*lmMst2)*pow2(Mst2)))/54. + (Mgl*(1715*(75 + 60*
        lmMsq*(-3 + 2*lmMst1 - 2*lmMst2) + 356*lmMst2 - 8*lmMst1*(22 + 15*
        lmMst2) + 120*pow2(lmMst2))*pow2(Mst1) + (-393767 + 420*lmMsq*(2194 +
        49*lmMst1 - 259*lmMst2) - 878948*lmMst2 - 1372*lmMst1*(31 + 15*lmMst2)
        + 44100*pow2(lmMsq) + 64680*pow2(lmMst2))*pow2(Mst2)))/185220. + (40*(
        8*Dmglst2*(14 - 15*lmMsq + 15*lmMst2) + (43 - 30*lmMsq + 30*lmMst2)*
        Mgl)*pow2(Msq)*pow2(Mst2) + (-24*Dmglst2*(-179 + 2*(87 + 32*lmMst1)*
        lmMst2 - 30*lmMsq*(-1 + 6*lmMst2) + 90*pow2(lmMsq) + 208*pow2(lmMst2))
        - 12*Mgl*(277 + 514*lmMst2 + 64*lmMst1*(1 + lmMst2) - 30*lmMsq*(5 + 6*
        lmMst2) + 90*pow2(lmMsq) + 208*pow2(lmMst2)))*pow4(Msq) + (90*Dmglst2*(
        13 - 28*lmMsq + 28*lmMst2) + 5*(67 - 84*lmMsq + 84*lmMst2)*Mgl)*pow4(
        Mst2))/(108.*pow2(Mst1))))/pow4(Msq) - (16*Dmglst2*(pow2(Msq)*pow2(
        Mst1)*(-9080791 + 54000*B4 - 54000*D3 + 27000*DN - 607500*lmMsq -
        2541540*lmMst1 + 202500*lmMsq*lmMst1 - 8582460*lmMst2 - 202500*lmMsq*
        lmMst2 + 867600*lmMst1*lmMst2 + 648000*lmMt - 324000*lmMst1*lmMt +
        324000*lmMst2*lmMt - 1800*pow2(lmMst1) + 202500*lmMst2*pow2(lmMst1) -
        2161800*pow2(lmMst2) + 1525500*lmMst1*pow2(lmMst2) + 4500*pow3(lmMst1)
        - 1732500*pow3(lmMst2)) + 33750*(2*pow2(Mst1)*pow2(Mst2) + 2*(-3 +
        lmMst1 - lmMst2)*pow4(Mst1) + (17 - 20*lmMsq + 20*lmMst2)*pow4(Mst2)))
        + 5*Mgl*(pow2(Msq)*pow2(Mst1)*(-19319827 - 19042320*lmMst2 - 518400*
        lmMt - 1036800*lmMst2*lmMt + 3600*(701 + 186*lmMst2)*pow2(lmMst1) -
        5943600*pow2(lmMst2) - 324000*lmMsq*(-2 + lmMst1 - lmMst2 - 2*lmMst1*
        lmMst2 + pow2(lmMst1) + pow2(lmMst2)) + 240*lmMst1*(35693 + 5610*lmMst2
        + 4320*lmMt + 3690*pow2(lmMst2)) + 295200*pow3(lmMst1) - 1850400*pow3(
        lmMst2)) + 48*((9608 - 15*lmMsq*(218 + 195*lmMst1 - 75*lmMst2) - 5250*
        lmMst2 + 15*lmMst1*(568 + 375*lmMst2) + 900*pow2(lmMsq) - 1350*pow2(
        lmMst1) - 3375*pow2(lmMst2))*pow2(Mst1)*pow2(Mst2) - 1125*(-4 + lmMst1*
        (2 - 6*lmMst2) - 2*lmMst2 + 3*pow2(lmMst1) + 3*pow2(lmMst2))*pow4(Mst1)
        + (6517 - 15*lmMsq*(532 + 75*lmMst1 - 195*lmMst2) + 4980*lmMst2 + 375*
        lmMst1*(8 + 3*lmMst2) - 900*pow2(lmMsq) - 2025*pow2(lmMst2))*pow4(Mst2)
        )))/(243000.*pow2(Msq)*pow2(Mst2)))/Mgl)) - (pow2(Mt)*pow2(s2t)*(pow2(
        Mst1)*pow2(Mst2)*(397.2878086419753 + 48*B4 - 60*lmMsq - 5*lmMst1*(43 -
        4*lmMsq + 4*pow2(lmMsq)) - (28*pow2(lmMst1))/3. + (lmMst2*(901 - 780*
        lmMst1 + 60*lmMsq*(-1 + 2*lmMst1) + 60*pow2(lmMsq) + 32*pow2(lmMst1)))/
        3. - (8*(-101 + 15*lmMsq + 44*lmMst1)*pow2(lmMst2))/3. + (320*pow3(
        lmMst2))/3. + ((-20*(6*Dmglst2*(15 + lmMst1*(-4 + 6*lmMsq - 6*lmMst2) +
        4*lmMst2 - 6*lmMsq*lmMst2 + 6*pow2(lmMst2)) + Mgl*(9 + 4*lmMst1*(-4 +
        3*lmMsq - 3*lmMst2) + 16*lmMst2 - 12*lmMsq*lmMst2 + 12*pow2(lmMst2)))*
        pow2(Mst1))/(9.*pow2(Msq)) + Dmglst2*(300.1378086419753 + (296*B4)/3. -
        (4*DN)/3. - 340*lmMsq + lmMst1*(40.55555555555556 - 60*lmMsq - 20*pow2(
        lmMsq)) + (28*pow2(lmMst1))/3. + lmMst2*(428.77777777777777 - 84*lmMst1
        + 20*lmMsq*(3 + 2*lmMst1) + 20*pow2(lmMsq) + (32*pow2(lmMst1))/3.) - (
        8*(-28 + 15*lmMsq + 60*lmMst1)*pow2(lmMst2))/3. + (448*pow3(lmMst2))/3.
        ))/Mgl) + (536.7597736625514 + 48*B4 - 55*lmMsq - lmMst1*(
        228.77777777777777 - 10*lmMsq + 20*pow2(lmMsq)) + (326*pow2(lmMst1))/9.
         + (lmMst2*(2635 - 3160*lmMst1 + 90*lmMsq*(-1 + 4*lmMst1) + 180*pow2(
        lmMsq) + 528*pow2(lmMst1)))/9. + (314.8888888888889 - 40*lmMsq - 208*
        lmMst1)*pow2(lmMst2) - (16*pow3(lmMst1))/9. - (Dmglst2*(
        309.29708504801096 - (296*B4)/3. + (4*DN)/3. + 375*lmMsq + lmMst1*(
        66.91358024691358 + 30*lmMsq + 20*pow2(lmMsq)) + (46*pow2(lmMst1))/3. -
        lmMst2*(10*lmMsq*(3 + 4*lmMst1) + 20*pow2(lmMsq) + (4*(11723 - 702*
        lmMst1 + 756*pow2(lmMst1)))/81.) + (2*(-75 + 60*lmMsq + 344*lmMst1)*
        pow2(lmMst2))/3. - (16*pow3(lmMst1))/3. - (560*pow3(lmMst2))/3.))/Mgl +
        (1360*pow3(lmMst2))/9.)*pow4(Mst1) + (368.5208333333333 + (110*B4)/3. -
        DN/3. - 120*lmMsq + 20*pow2(lmMsq) - (lmMst1*(199 - 30*lmMsq + 30*pow2(
        lmMsq)))/3. + (lmMst2*(1227 - 322*lmMst1 + 30*lmMsq*(-5 + 2*lmMst1) +
        30*pow2(lmMsq) - 16*pow2(lmMst1)))/3. - (16*pow2(lmMst1))/3. - 2*(-99 +
        10*lmMsq + 24*lmMst1)*pow2(lmMst2) + (160*pow3(lmMst2))/3. + ((-10*(3*
        Dmglst2*(59 + 8*lmMst1*(-2 + 3*lmMsq - 3*lmMst2) + 16*lmMst2 - 24*
        lmMsq*lmMst2 + 24*pow2(lmMst2)) + Mgl*(21 + 8*lmMst1*(-4 + 3*lmMsq - 3*
        lmMst2) + 32*lmMst2 - 24*lmMsq*lmMst2 + 24*pow2(lmMst2)))*pow2(Mst1))/(
        9.*pow2(Msq)) + Dmglst2*(371.73935185185184 + (262*B4)/3. - (5*DN)/3. -
        300*lmMsq + lmMst1*(57.666666666666664 - 30*lmMsq - 10*pow2(lmMsq)) +
        20*pow2(lmMsq) + (lmMst2*(1559 - 130*lmMst1 + 30*lmMsq*(-1 + 2*lmMst1)
        + 30*pow2(lmMsq) - 16*pow2(lmMst1)))/3. + (16*pow2(lmMst1))/3. - (2*(-
        217 + 30*lmMsq + 136*lmMst1)*pow2(lmMst2))/3. + 96*pow3(lmMst2)))/Mgl -
        (5*(Mgl*(245 + 12*lmMsq*(-1 + 36*lmMst1 - 36*lmMst2) + 336*lmMst2 -
        108*lmMst1*(3 + 4*lmMst2) + 432*pow2(lmMst2)) + Dmglst2*(3053 + 12*
        lmMsq*(-1 + 180*lmMst1 - 180*lmMst2) + 768*lmMst2 - 108*lmMst1*(7 + 20*
        lmMst2) + 2160*pow2(lmMst2)))*pow4(Mst1))/(216.*Mgl*pow4(Msq)))*pow4(
        Mst2) + ((2*OepS2*(-3*Mgl*(627*pow2(Mst1)*pow2(Mst2) + 1066*pow4(Mst1)
        + 189*pow4(Mst2)) + Dmglst2*(8163*pow2(Mst1)*pow2(Mst2) + 23734*pow4(
        Mst1) + 945*pow4(Mst2))))/729. - (S2*(Dmglst2*(9*(23989 + 27210*lmMst1
        - 27210*lmMst2)*pow2(Mst1)*pow2(Mst2) + 2*(525961 + 356010*lmMst1 -
        356010*lmMst2)*pow4(Mst1) + 81*(-453 + 350*lmMst1 - 350*lmMst2)*pow4(
        Mst2)) - 15*Mgl*(9*(685 + 418*lmMst1 - 418*lmMst2)*pow2(Mst1)*pow2(
        Mst2) + 2*(6143 + 3198*lmMst1 - 3198*lmMst2)*pow4(Mst1) + 81*(-1 + 14*
        lmMst1 - 14*lmMst2)*pow4(Mst2))))/540. + ((2*((-5*(3*Dmglst2*(95 + 12*
        lmMsq*(-2 + lmMst1 - lmMst2) + 32*lmMst2 - 4*lmMst1*(2 + 3*lmMst2) +
        12*pow2(lmMst2)) + Mgl*(71 + 12*lmMsq*(-2 + lmMst1 - lmMst2) + 40*
        lmMst2 - 4*lmMst1*(4 + 3*lmMst2) + 12*pow2(lmMst2))))/pow2(Msq) + (48*(
        1 + lmMst2)*(-Dmglst2 + 3*Dmglst2*lmMst2 + Mgl + lmMst2*Mgl))/pow2(
        Mst1)))/9. - (5*((Mgl*(242 + 24*lmMsq*(1 + 18*lmMst1 - 18*lmMst2) +
        300*lmMst2 - 108*lmMst1*(3 + 4*lmMst2) + 432*pow2(lmMst2)) + Dmglst2*(
        3122 + 24*lmMsq*(1 + 90*lmMst1 - 90*lmMst2) + 732*lmMst2 - 108*lmMst1*(
        7 + 20*lmMst2) + 2160*pow2(lmMst2)))*pow2(Mst1) + (Mgl*(911 + 12*lmMsq*
        (-37 + 18*lmMst1 - 18*lmMst2) + 606*lmMst2 - 54*lmMst1*(3 + 4*lmMst2) +
        216*pow2(lmMst2)) + Dmglst2*(5591 + 12*lmMsq*(-181 + 90*lmMst1 - 90*
        lmMst2) + 2550*lmMst2 - 54*lmMst1*(7 + 20*lmMst2) + 1080*pow2(lmMst2)))
        *pow2(Mst2)))/(216.*pow4(Msq)))*pow6(Mst2))/Mgl))/pow3(Mst2) + (pow4(
        Mt)*(20*(567143 - 125172*lmMst2 + 648*lmMst2*(-38 + 5*lmMst2)*lmMt -
        9720*(lmMsq*(lmMst1*(9 - 2*lmMt) + 2*lmMst2*(-6 + lmMt) + 3*(-4 + lmMt)
        ) + lmMt) - 9720*(lmMst1 - lmMst2)*pow2(lmMsq) + 240300*pow2(lmMst2) -
        108*lmMst1*(2111 + lmMst1*(457 - 12*lmMst2 - 126*lmMt) - 420*lmMt + 4*
        lmMst2*(454 + 39*lmMt) + 750*pow2(lmMst2) - 288*pow2(lmMt)) - 15552*(1
        + 2*lmMst2)*pow2(lmMt) + 216*pow3(lmMst1) + 79488*pow3(lmMst2))*pow4(
        Mst1) + (108*S2*(-48*Dmglst2*(601 + 2790*lmMst1 - 2790*lmMst2)*pow2(
        Mst1)*pow2(Mst2) + 72*(457 + 550*lmMst1 - 550*lmMst2)*Mgl*pow2(Mst1)*
        pow2(Mst2) + 40*(4082 + 3045*lmMst1 - 3045*lmMst2)*Mgl*pow4(Mst1) -
        675*Dmglst2*(-27 + 14*lmMst1 - 14*lmMst2)*pow4(Mst2) + 81*(-79 + 70*
        lmMst1 - 70*lmMst2)*Mgl*pow4(Mst2)))/Mgl + ((-24*pow2(Mst1)*pow2(Mst2)*
        (2700*Mgl*(-25 + lmMst1*(-5 + 6*lmMsq - 6*lmMst2) + (4 - 6*lmMsq)*
        lmMst2 + lmMt + 6*pow2(lmMst2))*pow2(Mst1) + 300*Dmglst2*(-313 + 3*
        lmMst1*(-23 + 42*lmMsq - 42*lmMst2) - 42*(-2 + 3*lmMsq)*lmMst2 - 15*
        lmMt + 126*pow2(lmMst2))*pow2(Mst1) + 2*Mgl*pow2(Msq)*(-71744 + 28755*
        lmMst2 + 15255*lmMt + 17820*lmMst2*lmMt + 4050*(lmMst1 - lmMst2)*pow2(
        lmMsq) - 540*(-15 + 4*lmMt)*pow2(lmMst1) - 89910*pow2(lmMst2) + 6480*
        lmMt*pow2(lmMst2) + 4050*lmMsq*(-11 + lmMst1 - 2*lmMst2 - 2*lmMst1*
        lmMst2 + lmMt + 2*pow2(lmMst2)) + 270*lmMst1*(221 + 279*lmMst2 - 66*
        lmMt - 16*lmMst2*lmMt + 112*pow2(lmMst2) - 24*pow2(lmMt)) + 6480*pow2(
        lmMt) + 6480*lmMst2*pow2(lmMt) - 30240*pow3(lmMst2)) + 3*Dmglst2*pow2(
        Msq)*(93487 + 196530*lmMst2 + 9030*lmMt - 7380*lmMst2*lmMt + 2700*(
        lmMst1 - lmMst2)*pow2(lmMsq) - 360*(23 + 8*lmMst2 - 4*lmMt)*pow2(
        lmMst1) + 68040*pow2(lmMst2) + 1440*lmMt*pow2(lmMst2) + 300*lmMsq*(29 +
        lmMst1*(9 - 18*lmMst2) + 6*lmMst2 - 15*lmMt + 18*pow2(lmMst2)) - 8640*
        pow2(lmMt) - 4320*lmMst2*pow2(lmMt) + 20*lmMst1*(-5846 + 321*lmMt - 12*
        lmMst2*(209 + 12*lmMt) + 576*pow2(lmMst2) + 216*pow2(lmMt)) - 8640*
        pow3(lmMst2))))/pow2(Msq) + (9*(9*Mgl*(800*(19 + (-8 + 6*lmMsq)*lmMst2
        + lmMst1*(7 - 6*lmMsq + 6*lmMst2) + lmMt - 6*pow2(lmMst2))*pow2(Msq)*
        pow2(Mst1) + (813 - 3360*lmMst2 - 2400*(lmMst1 - lmMst2)*pow2(lmMsq) -
        1280*(1 + lmMst2)*pow2(lmMst1) + 2400*lmMsq*(7 + lmMst1 - 2*lmMst2 + 2*
        lmMst1*lmMst2 + lmMt - 2*pow2(lmMst2)) + 30880*pow2(lmMst2) - 80*lmMt*(
        175 + 96*lmMst2 + 32*pow2(lmMst2)) - 80*lmMst1*(167 + 258*lmMst2 - 32*(
        1 + lmMst2)*lmMt + 112*pow2(lmMst2)) - 3840*pow2(lmMt) + 10240*pow3(
        lmMst2))*pow4(Msq) + 200*(38 + 4*(-2 + 3*lmMsq)*lmMst2 + lmMst1*(7 -
        12*lmMsq + 12*lmMst2) + lmMt - 12*pow2(lmMst2))*pow4(Mst1)) - 5*
        Dmglst2*(-160*(439 + 6*(-8 + 21*lmMsq)*lmMst2 + 3*lmMst1*(17 - 42*lmMsq
        + 42*lmMst2) - 3*lmMt - 126*pow2(lmMst2))*pow2(Msq)*pow2(Mst1) + (7335
        + 169856*lmMst2 + 4320*(lmMst1 - lmMst2)*pow2(lmMsq) + 2304*(-1 +
        lmMst2)*pow2(lmMst1) + 21984*pow2(lmMst2) + 480*lmMsq*(-79 - 30*lmMst2
        - 9*lmMst1*(-3 + 2*lmMst2) + 3*lmMt + 18*pow2(lmMst2)) + 48*lmMt*(31 -
        128*lmMst2 + 96*pow2(lmMst2)) + 144*lmMst1*(-269 + 32*lmMt - 2*lmMst2*(
        63 + 16*lmMt) + 112*pow2(lmMst2)) - 18432*pow3(lmMst2))*pow4(Msq) + 40*
        (-1330 + lmMst1*(-51 + 396*lmMsq - 396*lmMst2) + (48 - 396*lmMsq)*
        lmMst2 + 3*lmMt + 396*pow2(lmMst2))*pow4(Mst1)))*pow4(Mst2))/pow4(Msq)
        + 160*OepS2*(9*Dmglst2*pow2(Mst2)*(496*pow2(Mst1) + 35*pow2(Mst2)) -
        Mgl*(1320*pow2(Mst1)*pow2(Mst2) + 4060*pow4(Mst1) + 189*pow4(Mst2))) +
        (90*(-32*(-5*Mgl*(58 - 36*lmMst2 + 6*lmMst1*(4 + 3*lmMst2) - 3*lmMt -
        9*lmMst2*lmMt + 3*lmMsq*(5 - 6*lmMst1 + 9*lmMst2 + 3*lmMt) - 9*pow2(
        lmMsq) - 18*pow2(lmMst2)) - 15*Dmglst2*(64 - 20*lmMst2 + 6*lmMst1*(2 +
        3*lmMst2) + 3*lmMt - 3*lmMst2*lmMt + lmMsq*(5 - 18*lmMst1 + 21*lmMst2 +
        3*lmMt) - 3*pow2(lmMsq) - 18*pow2(lmMst2)))*pow2(Msq) + 5*(3*Dmglst2*(
        1693 - 136*lmMst2 + 6*lmMst1*(19 + 84*lmMst2) - 4*lmMsq*(5 + 126*lmMst1
        - 132*lmMst2 - 6*lmMt) + 42*lmMt - 24*lmMst2*lmMt - 24*pow2(lmMsq) -
        504*pow2(lmMst2)) + Mgl*(1271 - 216*lmMst2 + 6*lmMst1*(41 + 60*lmMst2)
        - 12*lmMsq*(5 + 30*lmMst1 - 36*lmMst2 - 6*lmMt) + 30*lmMt - 72*lmMst2*
        lmMt - 72*pow2(lmMsq) - 360*pow2(lmMst2)))*pow2(Mst1) + 15*(Dmglst2*(
        1535 - 184*lmMst2 + 18*lmMst1*(7 + 20*lmMst2) - 6*(-9 + 28*lmMst2)*lmMt
        + 4*lmMsq*(1 - 90*lmMst1 + 132*lmMst2 + 42*lmMt) - 168*pow2(lmMsq) -
        360*pow2(lmMst2)) + Mgl*(403 - 24*lmMst2 + 18*lmMst1*(3 + 4*lmMst2) -
        12*lmMsq*(1 + 6*lmMst1 - 12*lmMst2 - 6*lmMt) - 18*(1 + 4*lmMst2)*lmMt -
        72*pow2(lmMsq) - 72*pow2(lmMst2)))*pow2(Mst2) + (2304*(1 + lmMst2)*(-
        Dmglst2 + 3*Dmglst2*lmMst2 + Mgl + lmMst2*Mgl)*pow4(Msq))/pow2(Mst1))*
        pow6(Mst2))/pow4(Msq))/Mgl))/(14580.*pow5(Mst2)) + Mt*((s2t*pow2(Mst2)*
        ((6*(Dmglst2*(1677 - 14*lmMst1 + 14*lmMst2) + 24*(22 + lmMst1 - lmMst2)
        *Mgl)*S2*pow2(s2t))/Mgl + ((-1 + 2*lmMst2)*shiftst3*(4*pow2(Mt) + (1 +
        lmMst1 - lmMst2)*pow2(Mst1)*pow2(s2t)))/pow2(Mst1)))/6. + pow3(s2t)*(
        pow2(Mst1)*(150.95520576131688 + 2*B4 + (2*D3)/9. - (2*DN)/9. - (335*
        lmMsq)/18. + lmMst1*(94.25864197530865 - (325*lmMsq)/18. + (5*pow2(
        lmMsq))/6.) + (25*pow2(lmMsq))/6. + (5*(-149 + 18*lmMsq)*pow2(lmMst1))/
        27. - (lmMst2*(53429 - 183450*lmMst1 + 2250*lmMsq*(-7 + 6*lmMst1) +
        1350*pow2(lmMsq) + 29610*pow2(lmMst1)))/1620. + ((-3683 + 270*lmMsq +
        2595*lmMst1)*pow2(lmMst2))/54. - (55*pow3(lmMst1))/27. + (Dmglst2*((-5*
        lmMsq*(7 + 6*lmMst2))/3. + 5*pow2(lmMsq) + (-5430043 - 1334580*lmMst2 +
        900*(859 - 3690*lmMst2)*pow2(lmMst1) - 846900*pow2(lmMst2) + 120*
        lmMst1*(44309 + 12315*lmMst2 + 56475*pow2(lmMst2)) - 45000*pow3(lmMst1)
        - 3411000*pow3(lmMst2))/121500.))/Mgl - (749*pow3(lmMst2))/27.) - (((
        5556600*(150*Dmglst2*(5 - 2*lmMst1 + 2*lmMst2) + Mgl*(1503 - 60*lmMsq*(
        12 + 5*lmMst1 - 5*lmMst2) - 950*lmMst2 + 10*lmMst1*(167 + 210*lmMst2) -
        900*pow2(lmMst1) - 1200*pow2(lmMst2))))/pow2(Msq) + (1250235000*lmMsq*(
        1 - 2*lmMst1 + 2*lmMst2)*(12*Dmglst2 + Mgl + 8*lmMst1*Mgl - 8*lmMst2*
        Mgl) - 3*Mgl*(79604910131 + 34248010860*lmMst2 - 1058400*(78238 +
        33285*lmMst2)*pow2(lmMst1) - 83251627200*pow2(lmMst2) + 1260*lmMst1*(-
        29738761 + 131792640*lmMst2 + 59005800*pow2(lmMst2)) - 1296540000*pow3(
        lmMst1) - 37821924000*pow3(lmMst2)) + 20*Dmglst2*(32030812049 +
        1902812940*lmMst2 + 793800*(433 + 6552*lmMst2)*pow2(lmMst1) + 76998600*
        pow2(lmMst2) - 1260*lmMst1*(2674409 + 333900*lmMst2 + 7514640*pow2(
        lmMst2)) - 311169600*pow3(lmMst1) + 4578638400*pow3(lmMst2)))/pow2(
        Mst2))*pow4(Mst1))/(3.000564e9*Mgl) - (pow2(Mst2)*(12760 - 720*B4 + 72*
        D3 - 36*DN - 1980*lmMsq + 540*pow2(lmMsq) - 3*lmMst1*(5471 - 1950*lmMsq
        + 630*pow2(lmMsq)) - 954*pow2(lmMst1) + 3*lmMst2*(8495 - 6750*lmMst1 +
        210*lmMsq*(-11 + 6*lmMst1) + 630*pow2(lmMsq) + 234*pow2(lmMst1)) - 18*(
        -1324 + 210*lmMsq + 399*lmMst1)*pow2(lmMst2) - 288*pow3(lmMst1) + 6768*
        pow3(lmMst2) + (324*(((Mgl*(-725 + 30*lmMsq*(5 + lmMst1 - lmMst2) -
        793*lmMst2 + lmMst1*(643 + 150*lmMst2) - 90*pow2(lmMst1) - 60*pow2(
        lmMst2)) + 50*Dmglst2*(20 + 6*lmMsq*(-4 + 3*lmMst1 - 3*lmMst2) + 27*
        lmMst2 - 3*lmMst1*(1 + 6*lmMst2) + 18*pow2(lmMst2)))*pow2(Mst1))/(270.*
        pow2(Msq)) + Dmglst2*(16.256172839506174 - (4*B4)/3. + (16*D3)/9. - (
        10*DN)/9. + (80*lmMsq)/3. + lmMst1*(4.666666666666667 - 10*pow2(lmMsq))
        - (13*pow2(lmMst1))/3. + (2*lmMst2*(-237 - 248*lmMst1 + 90*lmMsq*lmMst1
        + 45*pow2(lmMsq) + 16*pow2(lmMst1)))/9. - ((-439 + 180*lmMsq + 336*
        lmMst1)*pow2(lmMst2))/9. + (304*pow3(lmMst2))/9.)))/Mgl + (324*(
        0.16013483965014577 - (5549*lmMsq)/17640. + (0.8266061980347694 - (43*
        lmMsq)/126.)*lmMst1 + ((-553 + 420*lmMsq + 780*lmMst1)*lmMst2)/1080. -
        (5*Dmglst2)/(36.*Mgl) - pow2(lmMsq)/42. - (4*pow2(lmMst1))/21. - (5*
        pow2(lmMst2))/9.)*pow4(Mst1))/pow4(Msq)))/324. + (-((15*Mgl*(373 + 706*
        lmMst2 + 64*lmMst1*(1 + lmMst2) - 30*lmMsq*(5 + 6*lmMst2) - 6*shiftst3
        + 12*lmMst2*shiftst3 + 90*pow2(lmMsq) + 304*pow2(lmMst2))*pow2(Msq) +
        30*Dmglst2*(-115 + (430 + 64*lmMst1)*lmMst2 - 30*lmMsq*(-1 + 6*lmMst2)
        + 90*pow2(lmMsq) + 400*pow2(lmMst2))*pow2(Msq) + 6*Mgl*(637 - 20*lmMsq*
        (9 + 20*lmMst1 - 20*lmMst2) - 495*lmMst2 + 25*lmMst1*(27 + 16*lmMst2) -
        400*pow2(lmMst2))*pow2(Mst1) - 100*Dmglst2*(-174 + 6*lmMsq*(6 + 17*
        lmMst1 - 17*lmMst2) + 73*lmMst2 - lmMst1*(109 + 102*lmMst2) + 102*pow2(
        lmMst2))*pow2(Mst1))*pow4(Mst2))/(540.*pow2(Msq)*pow2(Mst1)) + ((Mgl*(
        9*(3653969 - 140*lmMsq*(12899 + 5390*lmMst1 - 6230*lmMst2) + 2069970*
        lmMst2 + 37730*lmMst1*(-7 + 20*lmMst2) - 58800*pow2(lmMsq) - 813400*
        pow2(lmMst2))*pow2(Mst1) - (58456697 - 420*lmMsq*(63929 + 69090*lmMst1
        - 69930*lmMst2) + 374010*lmMst2 + 10290*lmMst1*(2573 + 2820*lmMst2) -
        176400*pow2(lmMsq) - 29194200*pow2(lmMst2))*pow2(Mst2)) + 514500*
        Dmglst2*((119 - 34*lmMst2 + 12*lmMsq*(1 - 10*lmMst1 + 10*lmMst2) + 2*
        lmMst1*(11 + 60*lmMst2) - 120*pow2(lmMst2))*pow2(Mst1) - 2*(337 - 24*
        lmMsq*(5 + 8*lmMst1 - 8*lmMst2) + 14*lmMst2 + 2*lmMst1*(53 + 96*lmMst2)
        - 192*pow2(lmMst2))*pow2(Mst2)))*pow4(Mst2))/(2.22264e7*pow4(Msq)) + (-
        45*Mgl*(664*OepS2 - 81*(1261 + 166*lmMst1 - 166*lmMst2)*S2)*pow2(Mst1)*
        pow2(Mst2) + 54*Dmglst2*(632*OepS2 + 9*(16193 - 1422*lmMst1 + 1422*
        lmMst2)*S2)*pow2(Mst1)*pow2(Mst2) - 15*Mgl*(3548*OepS2 - 27*(15148 +
        2661*lmMst1 - 2661*lmMst2)*S2)*pow4(Mst1) + 4*Dmglst2*(25328*OepS2 +
        27*(47051 - 18996*lmMst1 + 18996*lmMst2)*S2)*pow4(Mst1) + 3024*Dmglst2*
        OepS2*pow4(Mst2) - 5184*Mgl*OepS2*pow4(Mst2))/(4374.*pow2(Mst2)) + ((
        Mgl*(48*pow2(1 + lmMst2)*pow2(Msq) + 5*(43 - 30*lmMsq + 30*lmMst2)*
        pow2(Mst1)) + 8*Dmglst2*(24*pow2(lmMst2)*pow2(Msq) + 5*(14 - 15*lmMsq)*
        pow2(Mst1) + 3*lmMst2*(8*pow2(Msq) + 25*pow2(Mst1))))*pow6(Mst2))/(54.*
        pow2(Msq)*pow4(Mst1)))/Mgl) + (T1ep*(3*Mgl*(3*pow2(Mst1)*pow2(Mst2)*(-
        2990*Mst2*s2t*pow2(Mt) - 627*Mt*pow2(Mst2)*pow2(s2t) + 1760*pow3(Mt) +
        830*pow3(Mst2)*pow3(s2t)) + (-22804*Mst2*s2t*pow2(Mt) - 3198*Mt*pow2(
        Mst2)*pow2(s2t) + 16240*pow3(Mt) + 4435*pow3(Mst2)*pow3(s2t))*pow4(
        Mst1) + 27*(-46*Mst2*s2t*pow2(Mt) - 21*Mt*pow2(Mst2)*pow2(s2t) + 28*
        pow3(Mt) + 16*pow3(Mst2)*pow3(s2t))*pow4(Mst2)) + Dmglst2*Mst2*(-27*
        Mst2*pow2(Mst1)*(-800*Mst2*s2t*pow2(Mt) - 907*Mt*pow2(Mst2)*pow2(s2t) +
        1984*pow3(Mt) + 316*pow3(Mst2)*pow3(s2t)) + 2*s2t*(35601*Mst2*Mt*s2t +
        66168*pow2(Mt) - 12664*pow2(Mst2)*pow2(s2t))*pow4(Mst1) - 189*(20*pow3(
        Mst2)*pow3(Mt) - 15*Mt*pow2(s2t)*pow5(Mst2) + 4*pow3(s2t)*pow6(Mst2))))
        )/(729.*Mgl*pow5(Mst2)))) + (Al4p*Mt*twoLoopFlag*(12*Mt*pow2(s2t)*(2*
        Mst2*(4 + 3*lmMst2 - lmMst1*(1 + lmMst2) + pow2(lmMst2) + (Dmglst2*(6 +
        lmMst1 + lmMst2 - lmMst1*lmMst2 + pow2(lmMst2)))/Mgl) + ((2 - 6*lmMst1
        + (6 - 4*lmMst1)*lmMst2 + 4*pow2(lmMst2) + (Dmglst2*(10 + lmMst1*(6 -
        4*lmMst2) - 6*lmMst2 + 4*pow2(lmMst2)))/Mgl)*pow2(Mst1))/Mst2 + ((0.5 -
        7*lmMst1 + (7 - 4*lmMst1)*lmMst2 + 4*pow2(lmMst2) + (Dmglst2*(10.5 +
        lmMst1*(9 - 4*lmMst2) - 9*lmMst2 + 4*pow2(lmMst2)))/Mgl)*pow4(Mst1))/
        pow3(Mst2)) + 8*s2t*pow2(Mt)*(6*lmMst1 + 2*(-2 + lmMst1)*lmMst2 + pow2(
        lmMst1) - 3*pow2(lmMst2) + (2*(Dmglst2 + 2*Dmglst2*lmMst2 + ((1 -
        lmMst1 + lmMst2)*(2*Dmglst2 + (lmMst1 - lmMst2)*Mgl)*pow2(Mst1))/pow2(
        Mst2) - ((2*Dmglst2*lmMst2 + Mgl + lmMst2*Mgl)*pow2(Mst2))/pow2(Mst1) -
        ((Dmglst2*(-2 + 4*lmMst1 - 4*lmMst2) + Mgl*(-lmMst1 + lmMst2 - 2*
        lmMst1*lmMst2 + pow2(lmMst1) + pow2(lmMst2)))*pow4(Mst1))/pow4(Mst2)))/
        Mgl) + (16*Tbeta*pow2(Mst1)*pow3(Mt)*(Mgl*(13 + 2*lmMst1*(6 + 3*lmMst2
        - 2*lmMt) + 2*lmMt + 2*lmMst2*(-7 + 2*lmMt) - 6*pow2(lmMst2))*pow4(
        Mst1) + Dmglst2*(11 - 6*lmMst1*lmMst2 - 8*lmMst2*(-5 + lmMt) + 8*
        lmMst1*(-4 + lmMt) - 8*lmMt + 6*pow2(lmMst2))*pow4(Mst1) + 2*Dmglst2*((
        7 + 7*lmMst2 + lmMst1*(-5 + lmMt) - 2*lmMt - lmMst2*lmMt)*pow2(Mst1)*
        pow2(Mst2) + (5 + lmMst1*(-1 + lmMst2) + lmMst2 - pow2(lmMst2))*pow4(
        Mst2)) + 2*Mgl*((4 + lmMst1*(3 + 2*lmMst2 - lmMt) + lmMst2*(-4 + lmMt)
        + lmMt - 2*pow2(lmMst2))*pow2(Mst1)*pow2(Mst2) + (2 + lmMst1 + (-2 +
        lmMst1)*lmMst2 + lmMt - pow2(lmMst2))*pow4(Mst2))) + Tbeta*pow3(Mst2)*
        pow3(s2t)*(-4*Dmglst2*(2*(1 + 2*lmMst1 - lmMst2)*pow2(Mst2)*pow4(Mst1)
        + 2*(1 + lmMst1 - lmMst2 + 2*lmMst1*lmMst2 - 2*pow2(lmMst2))*pow2(Mst1)
        *pow4(Mst2) + (1 + 2*lmMst1 - 2*lmMst2)*pow6(Mst1) - 2*lmMst2*pow6(
        Mst2)) + Mgl*(2*(-16 + 6*lmMst2 - 2*lmMst1*(8 + 5*lmMst2) + 3*pow2(
        lmMst1) + 7*pow2(lmMst2))*pow2(Mst2)*pow4(Mst1) - 2*(-12 - 18*lmMst2 +
        2*lmMst1*(5 + 3*lmMst2) + pow2(lmMst1) - 7*pow2(lmMst2))*pow2(Mst1)*
        pow4(Mst2) + (3 + lmMst1*(2 - 32*lmMst2) - 2*lmMst2 + 16*(pow2(lmMst1)
        + pow2(lmMst2)))*pow6(Mst1) + 4*(1 + lmMst2)*pow6(Mst2))) - 2*Mt*
        MuSUSY*s2t*(2*Mgl*(-(Mst2*s2t*(-14 - 20*lmMst2 + 2*lmMst1*(5 + 3*
        lmMst2) + pow2(lmMst1) - 7*pow2(lmMst2))) + 8*Mt*(4 + 3*lmMst2 -
        lmMst1*(1 + lmMst2) + pow2(lmMst2)))*pow2(Mst1)*pow4(Mst2) + 4*Dmglst2*
        (Mst2*s2t*(-5 + 8*lmMst2 - 4*lmMst1*(2 + lmMst2) + 4*pow2(lmMst2)) +
        Mt*(65 + lmMst1*(34 - 20*lmMst2) - 26*lmMst2 + 20*pow2(lmMst2)))*pow6(
        Mst1) + Mgl*(Mst2*s2t*(-1 + 50*lmMst2 - 2*lmMst1*(25 + 32*lmMst2) + 20*
        pow2(lmMst1) + 44*pow2(lmMst2)) + Mt*(84 + 152*lmMst2 - 40*lmMst1*(3 +
        2*lmMst2) + 80*pow2(lmMst2)))*pow6(Mst1) + 8*Dmglst2*((Mst2*s2t*(-2 +
        3*lmMst2 - lmMst1*(3 + 2*lmMst2) + 2*pow2(lmMst2)) + Mt*(22 + lmMst1*(8
        - 6*lmMst2) - 4*lmMst2 + 6*pow2(lmMst2)))*pow2(Mst2)*pow4(Mst1) + (-(
        Mst2*s2t*(1 + lmMst1 - 2*lmMst2 + 2*lmMst1*lmMst2 - 2*pow2(lmMst2))) +
        2*Mt*(6 + lmMst1 + lmMst2 - lmMst1*lmMst2 + pow2(lmMst2)))*pow2(Mst1)*
        pow4(Mst2) + lmMst2*s2t*pow7(Mst2)) + 4*Mgl*((4*Mt*(5 + 6*lmMst2 -
        lmMst1*(4 + 3*lmMst2) + 3*pow2(lmMst2)) + Mst2*s2t*(-1 + 13*lmMst2 -
        lmMst1*(13 + 8*lmMst2) + pow2(lmMst1) + 7*pow2(lmMst2)))*pow2(Mst2)*
        pow4(Mst1) + (1 + lmMst2)*s2t*pow7(Mst2))))/(Mgl*Tbeta*pow2(Mst1)*pow5(
        Mst2))))/12. - (Al4p*xDmglst2*pow2(Dmglst2)*(twoLoopFlag*(Mt*(-(Mt*
        pow2(s2t)*(Mst2*(28 + 4*lmMst1 - 2*lmMst1*lmMst2 + 2*pow2(lmMst2)) + ((
        34 + 8*lmMst1 - 4*(2 + lmMst1)*lmMst2 + 4*pow2(lmMst2))*pow2(Mst1))/
        Mst2 + ((37.5 + lmMst1*(5 - 4*lmMst2) - 5*lmMst2 + 4*pow2(lmMst2))*
        pow4(Mst1))/pow3(Mst2))) + (2*s2t*pow2(Mt)*(11 - 18*lmMst2 + (6*(3 -
        lmMst1 + lmMst2)*pow2(Mst1))/pow2(Mst2) + (6*(-2 + 3*lmMst2)*pow2(Mst2)
        )/pow2(Mst1) + (6*(7 - 6*lmMst1 + 6*lmMst2)*pow4(Mst1))/pow4(Mst2)))/9.
         + pow3(s2t)*(((2 + 3*lmMst2)*pow2(Mst1))/3. + (2.6666666666666665 +
        lmMst2 + lmMst1*(-1 + 2*lmMst2) - 2*pow2(lmMst2))*pow2(Mst2) + ((0.5 -
        lmMst1 + lmMst2)*pow4(Mst1))/pow2(Mst2) + ((0.6666666666666666 -
        lmMst2)*pow4(Mst2))/pow2(Mst1)) + (4*pow3(Mt)*((29 - 18*lmMst1*(-1 +
        lmMst2) - 3*lmMst2 - 15*lmMt + 18*pow2(lmMst2))*pow2(Mst1)*pow2(Mst2) +
        (236 + 339*lmMst2 - 18*lmMst1*(13 + 4*lmMst2 - 3*lmMt) - 105*lmMt - 54*
        lmMst2*lmMt + 72*pow2(lmMst2))*pow4(Mst1) + (-79 - 18*lmMst1*(-2 +
        lmMst2) - 39*lmMst2 + 3*lmMt + 18*pow2(lmMst2))*pow4(Mst2)))/(27.*pow5(
        Mst2))) + (MuSUSY*s2t*pow2(Mt)*(-2*(4*Mt*(-31 + 3*lmMst1*(-2 + lmMst2)
        + 4*lmMst2 - 3*pow2(lmMst2)) + 3*Mst2*s2t*(4 + lmMst2 + lmMst1*(-1 + 2*
        lmMst2) - 2*pow2(lmMst2)))*pow2(Mst2)*pow4(Mst1) - 2*(Mst2*s2t*(10 +
        lmMst1*(-3 + 6*lmMst2) - 6*pow2(lmMst2)) + 4*Mt*(-14 + lmMst1*(-2 +
        lmMst2) - pow2(lmMst2)))*pow2(Mst1)*pow4(Mst2) + (-3*Mst2*s2t*(9 + 4*
        lmMst1*(-1 + lmMst2) + 4*lmMst2 - 4*pow2(lmMst2)) + Mt*(398 + lmMst1*(
        68 - 40*lmMst2) - 52*lmMst2 + 40*pow2(lmMst2)))*pow6(Mst1) + 2*(-2 + 3*
        lmMst2)*s2t*pow7(Mst2)))/(3.*Tbeta*pow2(Mst1)*pow5(Mst2))) + Al4p*
        threeLoopFlag*(-(s2t*pow3(Mt)*(1531.9349206349207 - (16*B4)/3. + (16*
        D3)/3. - (8*DN)/3. + (460*lmMsq)/9. + (736*lmMst1)/9. - (32*lmMt)/9. +
        30*pow2(lmMsq) + lmMst2*(907.1851851851852 - 60*lmMsq + 68*lmMst1 - (
        64*pow2(lmMst1))/3.) + (58*pow2(lmMst1))/3. + (2*(431 - 96*lmMst1)*
        pow2(lmMst2))/9. - (20*pow2(Mst1))/(9.*pow2(Msq)) + (128*pow3(lmMst2))/
        3. + (pow2(Mst1)*(1873.496051734274 - (16*B4)/3. + (16*D3)/3. - (8*DN)/
        3. - (100*lmMsq)/3. + (4*(253561 + 3375*lmMsq)*lmMst1)/2025. + (32*(5 +
        lmMst1 - lmMst2)*lmMt)/9. + (4*lmMst2*(732839 - 3375*lmMsq + 121560*
        lmMst1 - 11025*pow2(lmMst1)))/2025. + (2032*pow2(lmMst1))/135. - (4*(
        4292 + 5025*lmMst1)*pow2(lmMst2))/135. + (4*pow3(lmMst1))/27. + (4604*
        pow3(lmMst2))/27.))/pow2(Mst2) - (5*pow2(Mst2)*(27*(5 + 4*lmMsq - 4*
        lmMst2)*pow2(Mst1) + 2*(83 - 936*lmMsq + 936*lmMst2)*pow2(Mst2)))/(108.
        *pow4(Msq)) - (10*pow4(Mst1))/(9.*pow4(Msq)) - S2*(596.0571428571428 +
        84*lmMst1 - 84*lmMst2 + ((1403.352380952381 + 880*lmMst1 - 880*lmMst2)*
        pow2(Mst1))/pow2(Mst2) + (8*(17269 + 13785*lmMst1 - 13785*lmMst2)*pow4(
        Mst1))/(45.*pow4(Mst2))) + ((386.22610253272387 - (32*B4)/9. + (32*D3)/
        9. - (16*DN)/9. + (280*lmMsq)/3. + (756.6053917863442 - 80*lmMsq)*
        lmMst1 + (19352*pow2(lmMst1))/189. + lmMst2*(196.28349710254471 + 80*
        lmMsq - (57520*lmMst1)/189. + (88*pow2(lmMst1))/9.) + (8*(6787 - 5439*
        lmMst1)*pow2(lmMst2))/189. - (64*lmMt*(4 + 6*lmMst2 - 2*lmMst1*(3 +
        lmMst2) + pow2(lmMst1) + pow2(lmMst2)))/3. - (8*pow3(lmMst1))/9. + (
        664*pow3(lmMst2))/3.)*pow4(Mst1))/pow4(Mst2) + (64*(-2 + lmMst2 + 5*
        pow2(lmMst2))*pow4(Mst2))/(9.*pow4(Mst1)) + (16*OepS2*(660*pow2(Mst1)*
        pow2(Mst2) + 1838*pow4(Mst1) + 63*pow4(Mst2)))/(243.*pow4(Mst2)) + (10*
        (4*(3 - lmMst1 + lmMst2)*pow4(Mst1) + (-51 + 100*lmMsq - 100*lmMst2)*
        pow4(Mst2)))/(9.*pow2(Msq)*pow2(Mst2)) + (12*(957 + 64*lmMst1*(2 - 3*
        lmMst2) + 182*lmMst2 + 90*lmMsq*(-5 + 6*lmMst2) - 270*pow2(lmMsq) -
        624*pow2(lmMst2))*pow2(Mst2)*pow4(Msq) + 800*(8 - 15*lmMsq + 15*lmMst2)
        *pow2(Msq)*pow4(Mst2) + 315*(5 - 28*lmMsq + 28*lmMst2)*pow6(Mst2))/(
        108.*pow2(Mst1)*pow4(Msq)))) + (pow4(Mt)*(-(pow4(Msq)*(-3*pow2(Mst1)*
        pow2(Mst2)*(34872121 - 13890240*lmMst2 - 7840*OepS2 - 324*(9361 + 490*
        lmMst2)*S2 - 453600*(lmMst1 - lmMst2)*pow2(lmMsq) - 241920*(-2 +
        lmMst2)*pow2(lmMst1) + 7560*lmMst1*(206 - 128*lmMt + lmMst2*(636 + 64*
        lmMt) + 21*S2 - 224*pow2(lmMst2)) + 50400*lmMsq*(64 + 48*lmMst2 + 9*
        lmMst1*(-5 + 2*lmMst2) - 3*lmMt - 18*pow2(lmMst2)) - 5654880*pow2(
        lmMst2) - 5040*lmMt*(-31 - 240*lmMst2 + 96*pow2(lmMst2)) + 120960*pow2(
        lmMt) + 1935360*pow3(lmMst2)) - 8*(24896293 + 2509920*lmMst2 - 605920*
        OepS2 - 108*(105619 + 113610*lmMst2)*S2 - 170100*(lmMst1 - lmMst2)*
        pow2(lmMsq) + 90720*(-3 + lmMst2)*pow2(lmMst1) + 18900*lmMsq*(100 + 48*
        lmMst2 + 9*lmMst1*(-5 + 2*lmMst2) - 3*lmMt - 18*pow2(lmMst2)) -
        1563030*pow2(lmMst2) - 1890*lmMt*(-539 - 233*lmMst2 + 96*pow2(lmMst2))
        - 1890*lmMst1*(2075 + 281*lmMt - lmMst2*(1139 + 96*lmMt) - 6492*S2 +
        528*pow2(lmMst2)) - 226800*pow2(lmMt) + 907200*pow3(lmMst2))*pow4(Mst1)
        + 2903040*(2 + lmMst2 - 3*pow2(lmMst2))*pow4(Mst2))) + 50400*pow2(Msq)*
        ((778 + 6*(1 + 39*lmMsq)*lmMst2 - 3*(lmMst1 + 78*lmMsq*lmMst1 - 78*
        lmMst1*lmMst2 + lmMt) - 234*pow2(lmMst2))*pow2(Mst2)*pow4(Mst1) + 3*(
        265 - 34*lmMst2 + 12*lmMst1*(1 + 6*lmMst2) - 6*(-2 + lmMst2)*lmMt + 2*
        lmMsq*(5 - 36*lmMst1 + 39*lmMst2 + 3*lmMt) - 6*pow2(lmMsq) - 72*pow2(
        lmMst2))*pow2(Mst1)*pow4(Mst2) + (313 + 42*(-2 + 3*lmMsq)*lmMst2 + 3*
        lmMst1*(23 - 42*lmMsq + 42*lmMst2) + 15*lmMt - 126*pow2(lmMst2))*pow6(
        Mst1)) + 3150*(4*(1330 + 12*(-4 + 33*lmMsq)*lmMst2 + lmMst1*(51 - 396*
        lmMsq + 396*lmMst2) - 3*lmMt - 396*pow2(lmMst2))*pow2(Mst2)*pow6(Mst1)
        + 3*((4585 + 104*lmMst2 + 6*lmMst1*(-25 + 228*lmMst2) - 4*lmMsq*(5 +
        342*lmMst1 - 348*lmMst2 - 6*lmMt) + 66*lmMt - 24*lmMst2*lmMt - 24*pow2(
        lmMsq) - 1368*pow2(lmMst2))*pow4(Mst1)*pow4(Mst2) + (4317 - 280*lmMst2
        + 18*lmMst1*(1 + 60*lmMst2) - 78*(-3 + 4*lmMst2)*lmMt + 4*lmMsq*(7 -
        270*lmMst1 + 348*lmMst2 + 78*lmMt) - 312*pow2(lmMsq) - 1080*pow2(
        lmMst2))*pow2(Mst1)*pow6(Mst2)))))/(102060.*pow2(Mst1)*pow3(Mst2)*pow4(
        Msq)) + pow2(Mt)*(pow2(s2t)*(-((pow2(Mst1)*(990.1996399176954 + (524*
        B4)/3. - (10*DN)/3. - 720*lmMsq - lmMst1*(181.5185185185185 + 100*lmMsq
        + 20*pow2(lmMsq)) - (64*pow2(lmMst1))/3. + lmMst2*(20*lmMsq*(5 + lmMsq
        + 2*lmMst1) + (35429 + 3420*lmMst1 + 288*pow2(lmMst1))/27.) - (4*(79 +
        30*lmMsq + 168*lmMst1)*pow2(lmMst2))/3. - (40*(15 + lmMst1*(-4 + 6*
        lmMsq - 6*lmMst2) + (4 - 6*lmMsq)*lmMst2 + 6*pow2(lmMst2))*pow2(Mst1))/
        (3.*pow2(Msq)) + (640*pow3(lmMst2))/3.))/Mst2) + (4*((5*(397 + 36*
        lmMsq*(-2 + lmMst1 - lmMst2) + 78*lmMst2 - 6*lmMst1*(1 + 6*lmMst2) +
        36*pow2(lmMst2)))/pow2(Msq) + (48*(2 + lmMst2 - 3*pow2(lmMst2)))/pow2(
        Mst1))*pow3(Mst2))/9. - Mst2*(109.55663580246913 + (490*B4)/3. - (11*
        DN)/3. - (1970*lmMsq)/3. + lmMst1*(119.66666666666667 - 50*lmMsq - 10*
        pow2(lmMsq)) + 20*pow2(lmMsq) + lmMst2*(966.1111111111111 + 10*lmMsq*(1
        + lmMsq) + (74 + 20*lmMsq)*lmMst1 - (16*pow2(lmMst1))/3.) + (32*pow2(
        lmMst1))/3. + (54 - 20*lmMsq - (464*lmMst1)/3.)*pow2(lmMst2) - (10*(217
        + 8*lmMst1*(-1 + 6*lmMsq - 6*lmMst2) + (8 - 48*lmMsq)*lmMst2 + 48*pow2(
        lmMst2))*pow2(Mst1))/(3.*pow2(Msq)) + 160*pow3(lmMst2) - ((
        70.67129629629629 - (5*lmMsq)/18. + (5*(-7 + 20*lmMsq)*lmMst1)/2. + (
        17.77777777777778 - 50*(lmMsq + lmMst1))*lmMst2 + 50*pow2(lmMst2))*
        pow4(Mst1))/pow4(Msq)) + ((309.29708504801096 - (296*B4)/3. + (4*DN)/3.
         + 375*lmMsq + lmMst1*(66.91358024691358 + 30*lmMsq + 20*pow2(lmMsq)) +
        (46*pow2(lmMst1))/3. - lmMst2*(10*lmMsq*(3 + 4*lmMst1) + 20*pow2(lmMsq)
        + (4*(11723 - 702*lmMst1 + 756*pow2(lmMst1)))/81.) + (2*(-75 + 60*lmMsq
        + 344*lmMst1)*pow2(lmMst2))/3. - (16*pow3(lmMst1))/3. - (560*pow3(
        lmMst2))/3.)*pow4(Mst1) - (2*OepS2*(-9267*pow2(Mst1)*pow2(Mst2) +
        23734*pow4(Mst1) - 63*pow4(Mst2)))/729. + (S2*(3*(29123 - 92670*lmMst1
        + 92670*lmMst2)*pow2(Mst1)*pow2(Mst2) + 2*(525961 + 356010*lmMst1 -
        356010*lmMst2)*pow4(Mst1) + 27*(1387 - 70*lmMst1 + 70*lmMst2)*pow4(
        Mst2)))/540.)/pow3(Mst2) + (5*(2*(7969 + 12*lmMsq*(1 + 270*lmMst1 -
        270*lmMst2) + 42*lmMst2 - 54*lmMst1*(1 + 60*lmMst2) + 3240*pow2(lmMst2)
        )*pow2(Mst1)*pow3(Mst2) + (21275 + 12*lmMsq*(-541 + 270*lmMst1 - 270*
        lmMst2) + 6546*lmMst2 - 54*lmMst1*(1 + 60*lmMst2) + 3240*pow2(lmMst2))*
        pow5(Mst2)))/(216.*pow4(Msq))) + (MuSUSY*((72*pow2(Mt)*(5*pow2(Mst1)*(
        28405 - 288*B4 + 288*D3 - 144*DN + 5472*lmMst1 - 288*(-89 + 2*lmMst1)*
        lmMst2 + 576*(-2 + lmMst1 - lmMst2)*lmMt + 224*OepS2 - 324*(65 + 14*
        lmMst1 - 14*lmMst2)*S2 - 576*(-9 + 8*lmMst1)*pow2(lmMst2) + 4608*pow3(
        lmMst2)) + 4*pow2(Mst2)*(50134 - 270*B4 + 270*D3 - 135*DN + 120*(271 +
        24*lmMst1)*lmMst2 - 120*lmMt - 29646*S2 - 720*(-2 + 3*lmMst1)*pow2(
        lmMst2) + 2160*(lmMst1 + pow3(lmMst2)))) + 80*T1ep*(6*pow2(Mst1)*(1555*
        Mst2*Mt*s2t - 252*pow2(Mt) + 185*pow2(Mst2)*pow2(s2t)) + 63*s2t*(Mt +
        5*Mst2*s2t)*pow3(Mst2) + 17308*pow2(s2t)*pow4(Mst1)))/pow4(Mst2) +
        14580*pow2(s2t)*(112.1664609053498 - (100*B4)/9. + (112*D3)/9. - (62*
        DN)/9. + (1055*lmMsq)/9. + 15*pow2(lmMsq) - (2*lmMst1*(-544 + 210*lmMsq
        + 135*pow2(lmMsq)))/9. - (103*pow2(lmMst1))/9. + (lmMst2*(-7921 - 1476*
        lmMst1 + 90*lmMsq*(5 + 18*lmMst1) + 810*pow2(lmMsq) + 288*pow2(lmMst1))
        )/27. + (93.66666666666667 - 60*lmMsq - 112*lmMst1)*pow2(lmMst2) +
        pow2(Mst1)*((10*(405 + (-187 + 246*lmMsq)*lmMst2 + lmMst1*(187 - 246*
        lmMsq + 246*lmMst2) - 246*pow2(lmMst2)))/(27.*pow2(Msq)) - (
        73.78472976680384 + (164*B4)/9. - (176*D3)/9. + (94*DN)/9. - (260*
        lmMsq)/3. + lmMst1*(-96.55703703703703 + 40*lmMsq + 30*pow2(lmMsq)) + (
        7556*pow2(lmMst1))/135. - (2*lmMst2*(-99788 + 2005*lmMst1 + 6750*lmMsq*
        (2 + 3*lmMst1) + 10125*pow2(lmMsq) + 32775*pow2(lmMst1)))/675. + (2*(-
        3377 + 4050*lmMsq + 19155*lmMst1)*pow2(lmMst2))/135. + (10*pow3(lmMst1)
        )/27. - (5050*pow3(lmMst2))/27.)/pow2(Mst2)) + (304*pow3(lmMst2))/3. +
        pow4(Mst1)*((5*(216 + (-95 + 132*lmMsq)*lmMst2 + lmMst1*(95 - 132*lmMsq
        + 132*lmMst2) - 132*pow2(lmMst2)))/(54.*pow4(Msq)) + (536.1152102791342
         - (8*B4)/3. + (32*D3)/9. - (20*DN)/9. + 90*lmMsq - (123.11224321827497
         + 20*lmMsq*(1 + lmMsq))*lmMst1 - lmMst2*(17.33220122616948 - 20*lmMsq*
        (1 + lmMsq) + (133.04550264550264 - 40*lmMsq)*lmMst1 - (1180*pow2(
        lmMst1))/9.) - (15886*pow2(lmMst1))/945. + (149.85608465608465 - 40*
        lmMsq - (2812*lmMst1)/9.)*pow2(lmMst2) - (92*pow3(lmMst1))/27. + (4988*
        pow3(lmMst2))/27.)/pow4(Mst2)) - (32*(-2 + lmMst2 + 5*pow2(lmMst2))*
        pow4(Mst2))/(9.*pow4(Mst1)) - (6*(7400*OepS2 + 27*(1089707 - 5550*
        lmMst1 + 5550*lmMst2)*S2)*pow2(Mst1)*pow2(Mst2) + 40*(17308*OepS2 + 27*
        (93919 - 12981*lmMst1 + 12981*lmMst2)*S2)*pow4(Mst1) + 9*(1400*OepS2 +
        81*(116129 - 350*lmMst1 + 350*lmMst2)*S2)*pow4(Mst2))/(10935.*pow4(
        Mst2)) + (5*((291 + 2*(-103 + 84*lmMsq)*lmMst2 + 2*lmMst1*(103 - 84*
        lmMsq + 84*lmMst2) - 168*pow2(lmMst2))*pow4(Mst1) + (684 + 347*lmMst1 +
        7*(-77 + 78*lmMst1)*lmMst2 + 6*lmMsq*(32 - 91*lmMst1 + 91*lmMst2) -
        546*pow2(lmMst2))*pow4(Mst2)))/(27.*pow2(Msq)*pow2(Mst2)) + (pow2(Mst2)
        *(800*(-8 + 15*lmMsq - 15*lmMst2)*pow2(Msq)*pow2(Mst2) + 12*(-1149 +
        266*lmMst2 + 64*lmMst1*(-2 + 3*lmMst2) - 90*lmMsq*(-5 + 6*lmMst2) +
        270*pow2(lmMsq) + 1264*pow2(lmMst2))*pow4(Msq) - 5*((-4239 + 12*lmMsq*(
        -47 + 244*lmMst1 - 244*lmMst2) + 1324*lmMst2 - 8*lmMst1*(95 + 366*
        lmMst2) + 2928*pow2(lmMst2))*pow2(Mst1)*pow2(Mst2) + (-3834 + lmMst1*(-
        890 + 2328*lmMsq - 2328*lmMst2) + (890 - 2328*lmMsq)*lmMst2 + 2328*
        pow2(lmMst2))*pow4(Mst1) - 63*(-5 + 28*lmMsq - 28*lmMst2)*pow4(Mst2))))
        /(216.*pow2(Mst1)*pow4(Msq))) - (Mt*s2t*(pow4(Msq)*(-3*pow2(Mst1)*pow2(
        Mst2)*(433447 + 1058400*B4 - 23760*DN - 4255200*lmMsq - 1120*OepS2 +
        324*(-1387 + 70*lmMst1 - 70*lmMst2)*S2 + 129600*pow2(lmMsq) - 2160*
        lmMst1*(-359 + 150*lmMsq + 30*pow2(lmMsq)) + 720*lmMst2*(8503 + 666*
        lmMst1 + 90*lmMsq*(1 + 2*lmMst1) + 90*pow2(lmMsq) - 48*pow2(lmMst1)) +
        69120*pow2(lmMst1) - 4320*(-177 + 30*lmMsq + 232*lmMst1)*pow2(lmMst2) +
        1036800*pow3(lmMst2)) - 2*(10274911 + 3285360*B4 - 68040*DN - 13381200*
        lmMsq - 248800*OepS2 + 108*(-20803 + 46650*lmMst1 - 46650*lmMst2)*S2 +
        194400*pow2(lmMsq) - 3600*lmMst1*(167 + 405*lmMsq + 81*pow2(lmMsq)) -
        103680*pow2(lmMst1) + 720*lmMst2*(30469 + 2709*lmMst1 + 135*lmMsq*(11 +
        6*lmMst1) + 405*pow2(lmMsq) + 72*pow2(lmMst1)) - 6480*(-19 + 90*lmMsq +
        568*lmMst1)*pow2(lmMst2) + 3628800*pow3(lmMst2))*pow4(Mst1) + 414720*(2
        + lmMst2 - 3*pow2(lmMst2))*pow4(Mst2)) + 21600*pow2(Msq)*((1445 + 72*
        lmMsq*(-2 + 3*lmMst1 - 3*lmMst2) + 180*lmMst2 - 36*lmMst1*(1 + 6*
        lmMst2) + 216*pow2(lmMst2))*pow2(Mst2)*pow4(Mst1) + 2*(397 + 36*lmMsq*(
        -2 + lmMst1 - lmMst2) + 78*lmMst2 - 6*lmMst1*(1 + 6*lmMst2) + 36*pow2(
        lmMst2))*pow2(Mst1)*pow4(Mst2) + 6*(107 + 6*lmMsq*(-2 + 5*lmMst1 - 5*
        lmMst2) + 32*lmMst2 - 10*lmMst1*(2 + 3*lmMst2) + 30*pow2(lmMst2))*pow6(
        Mst1)) + 450*((37213 + 12*lmMsq*(-539 + 810*lmMst1 - 810*lmMst2) +
        6630*lmMst2 - 162*lmMst1*(1 + 60*lmMst2) + 9720*pow2(lmMst2))*pow4(
        Mst1)*pow4(Mst2) + 6*(1961 + 180*lmMsq*(-2 + 5*lmMst1 - 5*lmMst2) +
        675*lmMst2 - 45*lmMst1*(7 + 20*lmMst2) + 900*pow2(lmMst2))*pow2(Mst2)*
        pow6(Mst1) + (21275 + 12*lmMsq*(-541 + 270*lmMst1 - 270*lmMst2) + 6546*
        lmMst2 - 54*lmMst1*(1 + 60*lmMst2) + 3240*pow2(lmMst2))*pow2(Mst1)*
        pow6(Mst2))))/(pow2(Mst1)*pow3(Mst2)*pow4(Msq))))/(14580.*Tbeta)) + Mt*
        ((T1ep*(3*Mst2*pow2(Mst1)*(15840*Mst2*s2t*pow2(Mt) - 9267*Mt*pow2(Mst2)
        *pow2(s2t) + 17312*pow3(Mt) - 530*pow3(Mst2)*pow3(s2t)) + 63*pow3(Mst2)
        *(72*Mst2*s2t*pow2(Mt) - 3*Mt*pow2(Mst2)*pow2(s2t) + 4*pow3(Mt) - 10*
        pow3(Mst2)*pow3(s2t)) + 2*s2t*(35601*Mst2*Mt*s2t + 66168*pow2(Mt) -
        12664*pow2(Mst2)*pow2(s2t))*pow4(Mst1)))/(729.*pow4(Mst2)) + pow3(s2t)*
        (pow2(Mst1)*(92.97559533607682 + (32*B4)/9. - (32*D3)/9. + (16*DN)/9. +
        (275*lmMsq)/18. + (12.165925925925926 - (10*lmMsq)/3.)*lmMst1 + (15*
        pow2(lmMsq))/2. + (6011*pow2(lmMst1))/270. - (lmMst2*(-1551 + 15750*
        lmMsq + 40910*lmMst1 + 58350*pow2(lmMst1)))/1350. + ((5891 + 23190*
        lmMst1)*pow2(lmMst2))/270. + (5*pow3(lmMst1))/27. - (1157*pow3(lmMst2))
        /27.) - ((5*(5 - 2*lmMst1 + 2*lmMst2))/(18.*pow2(Msq)) + (
        213.49860925479342 + (1510169*lmMst2)/119070. + lmMsq*(5 - 10*lmMst1 +
        10*lmMst2) + ((433 + 6552*lmMst2)*pow2(lmMst1))/189. + (97*pow2(lmMst2)
        )/189. - (lmMst1*(2674409 + 333900*lmMst2 + 7514640*pow2(lmMst2)))/
        119070. - (56*pow3(lmMst1))/27. + (824*pow3(lmMst2))/27.)/pow2(Mst2))*
        pow4(Mst1) - pow2(Mst2)*(87.99989711934157 - (50*B4)/9. + (56*D3)/9. -
        (31*DN)/9. + (415*lmMsq)/9. + lmMst1*(64 - (70*lmMsq)/3. - 15*pow2(
        lmMsq)) - (103*pow2(lmMst1))/18. + (lmMst2*(-4160 - 882*lmMst1 + 90*
        lmMsq*(7 + 9*lmMst1) + 405*pow2(lmMsq) + 144*pow2(lmMst1)))/27. + (
        11.722222222222221 - 30*lmMsq - 56*lmMst1)*pow2(lmMst2) + (5*(42 +
        lmMst1*(9 - 18*lmMst2) + 2*lmMsq*(-32 + 9*lmMst1 - 9*lmMst2) + 55*
        lmMst2 + 18*pow2(lmMst2))*pow2(Mst1))/(18.*pow2(Msq)) + (152*pow3(
        lmMst2))/3. - (5*pow4(Mst1))/(36.*pow4(Msq))) - (((10*(844 - 6*lmMsq*(
        18 + 91*lmMst1 - 91*lmMst2) - 239*lmMst2 + lmMst1*(347 + 546*lmMst2) -
        546*pow2(lmMst2)))/pow2(Msq) + (3*(-1277 + 330*lmMst2 + 64*lmMst1*(-2 +
        3*lmMst2) - 90*lmMsq*(-5 + 6*lmMst2) + 270*pow2(lmMsq) + 1584*pow2(
        lmMst2)))/pow2(Mst1))*pow4(Mst2))/108. + (159*(200*OepS2 + 27*(21401 -
        150*lmMst1 + 150*lmMst2)*S2)*pow2(Mst1)*pow2(Mst2) + 20*(25328*OepS2 +
        27*(47051 - 18996*lmMst1 + 18996*lmMst2)*S2)*pow4(Mst1) + 9*(1400*OepS2
        + 81*(116129 - 350*lmMst1 + 350*lmMst2)*S2)*pow4(Mst2))/(21870.*pow2(
        Mst2)) + (16*(-2 + lmMst2 + 5*pow2(lmMst2))*pow6(Mst2))/(9.*pow4(Mst1))
        + (5*((405 - 434*lmMst2 + 10*lmMst1*(-13 + 60*lmMst2) + lmMsq*(564 -
        600*lmMst1 + 600*lmMst2) - 600*pow2(lmMst2))*pow4(Mst1)*pow4(Mst2) -
        160*(-8 + 15*lmMsq - 15*lmMst2)*pow2(Msq)*pow6(Mst2) + 2*(-2277 + 24*
        lmMsq*(25 + 61*lmMst1 - 61*lmMst2) - 220*lmMst2 - 4*lmMst1*(95 + 366*
        lmMst2) + 1464*pow2(lmMst2))*pow2(Mst1)*pow6(Mst2) + 9*(9 - 140*lmMsq +
        140*lmMst2)*pow8(Mst2)))/(432.*pow2(Mst1)*pow4(Msq)))))))/pow2(Mgl) - (
        Al4p*Mt*xDR2DRMOD*(Al4p*threeLoopFlag*((1296*Mt*MuSUSY*s2t*(-(s2t*(3636
        - 1800*lmMsq + 1080*pow2(lmMsq) - 24*lmMst1*(431 - 150*lmMsq + 90*pow2(
        lmMsq)) + 768*pow2(lmMst1) + 48*lmMst2*(363 + 30*lmMsq*(-4 + 3*lmMst1 -
        3*lmMst2) + 451*lmMst2 - lmMst1*(407 + 168*lmMst2) + 45*pow2(lmMsq) +
        16*pow2(lmMst1)) + 7296*pow3(lmMst2) + (24*Dmglst2*(45 + 30*lmMsq*(1 +
        3*lmMsq - 6*lmMst2)*(1 - 2*lmMst1 + 2*lmMst2) + 8*lmMst2*(29 + 8*pow2(
        lmMst1)) + 1068*pow2(lmMst2) - 2*lmMst1*(-83 + 430*lmMst2 + 336*pow2(
        lmMst2)) + 608*pow3(lmMst2)) - (384*(1 + lmMst2)*(4*Dmglst2*lmMst2 +
        Mgl + lmMst2*Mgl)*pow4(Mst2))/pow4(Mst1) + (8*pow2(Mst1)*(80*Dmglst2*(
        lmMst1 - lmMst2)*(14 - 15*lmMsq + 15*lmMst2)*pow2(Mst2) + 10*(lmMst1 -
        lmMst2)*(43 - 30*lmMsq + 30*lmMst2)*Mgl*pow2(Mst2) + 6*Dmglst2*pow2(
        Msq)*(96 + lmMst2*(173 + 30*lmMsq + 90*pow2(lmMsq) + 288*pow2(lmMst1))
        + (526 - 180*lmMsq)*pow2(lmMst2) - lmMst1*(13 + 30*lmMsq*(1 - 6*lmMst2)
        + 526*lmMst2 + 90*pow2(lmMsq) + 848*pow2(lmMst2)) + 560*pow3(lmMst2)) +
        3*Mgl*pow2(Msq)*(64 + 15*lmMst2*(33 - 10*lmMsq + 6*pow2(lmMsq)) + 288*(
        1 + lmMst2)*pow2(lmMst1) + 6*(173 - 30*lmMsq)*pow2(lmMst2) - lmMst1*(
        431 + 1326*lmMst2 - 30*lmMsq*(5 + 6*lmMst2) + 90*pow2(lmMsq) + 848*
        pow2(lmMst2)) + 560*pow3(lmMst2))) - 40*(8*Dmglst2*(14 - 15*lmMsq + 15*
        lmMst2) + (43 - 30*lmMsq + 30*lmMst2)*Mgl)*(pow4(Mst2) - 2*(lmMst1 -
        lmMst2)*(pow4(Mst1) + pow4(Mst2))))/(pow2(Msq)*pow2(Mst2)) + ((2*pow4(
        Mst1)*((24*Dmglst2*(96 + lmMst2*(285 + 30*lmMsq + 90*pow2(lmMsq) + 544*
        pow2(lmMst1)) + (590 - 180*lmMsq)*pow2(lmMst2) - lmMst1*(109 + 30*
        lmMsq*(1 - 6*lmMst2) + 590*lmMst2 + 90*pow2(lmMsq) + 1360*pow2(lmMst2))
        + 816*pow3(lmMst2)) + 12*Mgl*(80 + lmMst2*(479 - 150*lmMsq + 90*pow2(
        lmMsq)) + 544*(1 + lmMst2)*pow2(lmMst1) + 2*(631 - 90*lmMsq)*pow2(
        lmMst2) - lmMst1*(399 + 1806*lmMst2 - 30*lmMsq*(5 + 6*lmMst2) + 90*
        pow2(lmMsq) + 1360*pow2(lmMst2)) + 816*pow3(lmMst2)))*pow4(Msq) + (
        lmMst1 - lmMst2)*(90*Dmglst2*(13 - 28*lmMsq + 28*lmMst2) + 5*(67 - 84*
        lmMsq + 84*lmMst2)*Mgl)*pow4(Mst2)))/pow4(Mst2) + pow2(Mst2)*(-5*(18*
        Dmglst2*(13 - 28*lmMsq + 28*lmMst2) + (67 - 84*lmMsq + 84*lmMst2)*Mgl)*
        (pow2(Mst2) - 2*(lmMst1 - lmMst2)*(pow2(Mst1) + pow2(Mst2))) - (40*(8*
        Dmglst2*(14 - 15*lmMsq + 15*lmMst2) + (43 - 30*lmMsq + 30*lmMst2)*Mgl)*
        pow2(Msq)*pow2(Mst2) + (-12*Mgl*(335 + 654*lmMst2 + 64*lmMst1*(1 +
        lmMst2) - 30*lmMsq*(5 + 6*lmMst2) + 90*pow2(lmMsq) + 272*pow2(lmMst2))
        - 24*Dmglst2*(-115 + (366 + 64*lmMst1)*lmMst2 - 30*lmMsq*(-1 + 6*
        lmMst2) + 90*pow2(lmMsq) + 336*pow2(lmMst2)))*pow4(Msq) + (90*Dmglst2*(
        13 - 28*lmMsq + 28*lmMst2) + 5*(67 - 84*lmMsq + 84*lmMst2)*Mgl)*pow4(
        Mst2))/pow2(Mst1)))/pow4(Msq))/Mgl))/216. + (64*Mt*(2*Dmglst2*pow2(
        Mst2)*(pow2(Mst1)*pow2(Mst2)*(8 + 7*lmMst2 - 11*pow2(lmMst2) + lmMst1*(
        -3 + 5*lmMst2 + 6*pow2(lmMst2)) - 6*pow3(lmMst2)) + 2*(8 + 13*lmMst2 -
        8*pow2(lmMst2) + lmMst1*(-5 + 5*lmMst2 + 6*pow2(lmMst2)) - 6*pow3(
        lmMst2))*pow4(Mst1) + (1 - 2*lmMst2 - 3*pow2(lmMst2))*pow4(Mst2)) + (1
        + lmMst2)*Mgl*(2*(1 - 10*lmMst2 + 4*lmMst1*(2 + lmMst2) - 4*pow2(
        lmMst2))*pow2(Mst2)*pow4(Mst1) - 2*(1 + 5*lmMst2 - lmMst1*(3 + 2*
        lmMst2) + 2*pow2(lmMst2))*pow2(Mst1)*pow4(Mst2) + (7 - 32*lmMst2 + 4*
        lmMst1*(7 + 3*lmMst2) - 12*pow2(lmMst2))*pow6(Mst1) - 2*(1 + lmMst2)*
        pow6(Mst2))))/(9.*Mgl*pow2(Mst1)*pow5(Mst2))))/Tbeta + (3*pow3(s2t)*(
        12*pow2(Mst1)*(2*Dmglst2*(147 + 114*lmMst2 + 30*lmMsq*(-1 + 6*lmMst2) -
        90*pow2(lmMsq) + 512*lmMst2*pow2(lmMst1) - 16*pow2(lmMst2) - 64*lmMst1*
        (3 + 3*lmMst2 + 16*pow2(lmMst2)) + 512*pow3(lmMst2)) - Mgl*(175 - 30*
        lmMsq*(5 + 6*lmMst2) + lmMst2*(462 + 1024*lmMst1*(1 + lmMst2)) + 90*
        pow2(lmMsq) - 272*pow2(lmMst2) - 512*((1 + lmMst2)*pow2(lmMst1) + pow3(
        lmMst2)))) + (384*((1 + lmMst2)*Mgl*(1 + lmMst1*(2 - 32*lmMst2) - 2*
        lmMst2 + 16*(pow2(lmMst1) + pow2(lmMst2))) + 2*Dmglst2*(-2*lmMst1*(3 +
        2*lmMst2 + 16*pow2(lmMst2)) + lmMst2*(7 + 4*lmMst2 + 16*(pow2(lmMst1) +
        pow2(lmMst2)))))*pow4(Mst1))/pow2(Mst2) + (5*(18*Dmglst2*(13 - 28*lmMsq
        + 28*lmMst2) + (67 - 84*lmMsq + 84*lmMst2)*Mgl)*(pow2(Mst1) + 2*(lmMst1
        - lmMst2)*pow2(Mst2))*pow4(Mst2))/pow4(Msq) + (8*pow2(Mst2)*(40*
        Dmglst2*(14 - 15*lmMsq + 15*lmMst2)*pow2(Mst1) + 5*(43 - 30*lmMsq + 30*
        lmMst2)*Mgl*pow2(Mst1) - 3*Mgl*pow2(Msq)*(16 - 3*lmMst2*(133 - 50*lmMsq
        + 30*pow2(lmMsq)) - 32*(1 + lmMst2)*pow2(lmMst1) + 2*(-383 + 90*lmMsq)*
        pow2(lmMst2) + lmMst1*(463 + 846*lmMst2 - 30*lmMsq*(5 + 6*lmMst2) + 90*
        pow2(lmMsq) + 336*pow2(lmMst2)) - 304*pow3(lmMst2)) + 6*Dmglst2*pow2(
        Msq)*(80 + lmMst2*(-67 + 30*lmMsq + 90*pow2(lmMsq) + 32*pow2(lmMst1)) +
        6*(61 - 30*lmMsq)*pow2(lmMst2) - lmMst1*(-83 + 462*lmMst2 - 30*lmMsq*(-
        1 + 6*lmMst2) + 90*pow2(lmMsq) + 336*pow2(lmMst2)) + 304*pow3(lmMst2)))
        - (4*((-3*Mgl*(367 + 718*lmMst2 + 64*lmMst1*(1 + lmMst2) - 30*lmMsq*(5
        + 6*lmMst2) + 90*pow2(lmMsq) + 304*pow2(lmMst2)) - 6*Dmglst2*(-115 + (
        430 + 64*lmMst1)*lmMst2 - 30*lmMsq*(-1 + 6*lmMst2) + 90*pow2(lmMsq) +
        400*pow2(lmMst2)))*pow2(Msq) + 20*(lmMst1 - lmMst2)*(8*Dmglst2*(-14 +
        15*lmMsq - 15*lmMst2) + (-43 + 30*lmMsq - 30*lmMst2)*Mgl)*pow2(Mst1))*
        pow4(Mst2))/pow2(Mst1) - (8*(Mgl*(48*pow2(1 + lmMst2)*pow2(Msq) + 5*(43
        - 30*lmMsq + 30*lmMst2)*pow2(Mst1)) + 8*Dmglst2*(24*pow2(lmMst2)*pow2(
        Msq) + 5*(14 - 15*lmMsq)*pow2(Mst1) + 3*lmMst2*(8*pow2(Msq) + 25*pow2(
        Mst1))))*pow6(Mst2))/pow4(Mst1))/pow2(Msq)) + (256*Mt*(4*pow2(Mt)*(
        Dmglst2*pow2(Mst2)*(pow2(Mst1)*pow2(Mst2)*(119 + lmMst2*(224 - 15*lmMt)
        - 51*lmMt + 3*(5 + 6*lmMt)*pow2(lmMst2) + 18*lmMst1*(-4 + lmMt -
        lmMst2*(1 + lmMt) + pow2(lmMst2)) - 18*pow3(lmMst2)) + (263 + lmMst2*(
        962 - 213*lmMt) - 177*lmMt + (393 - 18*lmMt)*pow2(lmMst2) - 18*lmMst1*(
        26 - lmMst2*(-17 + lmMt) - 7*lmMt + pow2(lmMst2)) + 18*pow3(lmMst2))*
        pow4(Mst1) - 18*(-1 + 2*lmMst2 + 3*pow2(lmMst2))*pow4(Mst2)) + 9*(1 +
        lmMst2)*Mgl*((3 - 21*lmMst2 + 2*lmMst1*(8 + 3*lmMst2 - 3*lmMt) + 5*lmMt
        + 6*lmMst2*lmMt - 6*pow2(lmMst2))*pow2(Mst2)*pow4(Mst1) + (-1 - 7*
        lmMst2 + 2*lmMst1*(2 + lmMst2 - lmMt) + 3*lmMt + 2*lmMst2*lmMt - 2*
        pow2(lmMst2))*pow2(Mst1)*pow4(Mst2) + (12 - 45*lmMst2 + 2*lmMst1*(19 +
        6*lmMst2 - 6*lmMt) + 7*lmMt + 12*lmMst2*lmMt - 12*pow2(lmMst2))*pow6(
        Mst1) - 2*(1 + lmMst2)*pow6(Mst2))) - 27*pow2(Mst2)*pow2(s2t)*((1 +
        lmMst2)*Mgl*(2*(2 - 5*lmMst2 + lmMst1*(5 + 2*lmMst2) - 2*pow2(lmMst2))*
        pow2(Mst2)*pow4(Mst1) + 2*(-2*lmMst2*(2 + lmMst2) + lmMst1*(3 + 2*
        lmMst2))*pow2(Mst1)*pow4(Mst2) + (5 - 12*lmMst2 + 4*lmMst1*(3 + lmMst2)
        - 4*pow2(lmMst2))*pow6(Mst1) - 2*(1 + lmMst2)*pow6(Mst2)) + Dmglst2*(2*
        pow2(Mst2)*(8 + 19*lmMst2 - 5*pow2(lmMst2) + lmMst1*(-7 + 5*lmMst2 + 6*
        pow2(lmMst2)) - 6*pow3(lmMst2))*pow4(Mst1) + 2*pow2(Mst1)*(7 + 9*lmMst2
        - 8*pow2(lmMst2) + lmMst1*(-3 + 5*lmMst2 + 6*pow2(lmMst2)) - 6*pow3(
        lmMst2))*pow4(Mst2) + (17 + 51*lmMst2 - 4*pow2(lmMst2) + 4*lmMst1*(-6 +
        lmMst2 + 3*pow2(lmMst2)) - 12*pow3(lmMst2))*pow6(Mst1) - 2*(-1 + 2*
        lmMst2 + 3*pow2(lmMst2))*pow6(Mst2)))))/(pow2(Mst1)*pow5(Mst2)) - (12*
        s2t*pow2(Mt)*(pow2(Mst1)*(pow2(Mst1) - pow2(Mst2))*(40*(8*Dmglst2*(14 -
        15*lmMsq + 15*lmMst2) + (43 - 30*lmMsq + 30*lmMst2)*Mgl)*pow2(Msq)*
        pow6(Mst2) + (90*Dmglst2*(13 - 28*lmMsq + 28*lmMst2) + 5*(67 - 84*lmMsq
        + 84*lmMst2)*Mgl)*pow8(Mst2)) + pow4(Msq)*(12*Mgl*(-((143 + 302*lmMst2
        - 30*lmMsq*(5 + 6*lmMst2) + 90*pow2(lmMsq) + 48*pow2(lmMst2) + 128*
        lmMst1*pow2(1 + lmMst2) - 64*((1 + lmMst2)*pow2(lmMst1) + pow3(lmMst2))
        )*pow4(Mst1)*pow4(Mst2)) + (271 + 526*lmMst2 + 64*lmMst1*(1 + lmMst2) -
        30*lmMsq*(5 + 6*lmMst2) + 90*pow2(lmMsq) + 208*pow2(lmMst2))*pow2(Mst1)
        *pow6(Mst2) + 64*(1 + lmMst2)*(pow2(1 - lmMst1 + lmMst2)*pow2(Mst2)*
        pow6(Mst1) + (1 + 3*lmMst2 - lmMst1*(3 + 2*lmMst2) + pow2(lmMst1) +
        pow2(lmMst2))*pow8(Mst1)) - 32*pow2(1 + lmMst2)*pow8(Mst2)) - 8*
        Dmglst2*(-((121 + 90*lmMsq*(-1 + 6*lmMst2) - 270*pow2(lmMsq) + 2*
        lmMst2*(-325 + 96*pow2(lmMst1)) - 432*pow2(lmMst2) + 192*(-(lmMst1*(-1
        + lmMst2 + 2*pow2(lmMst2))) + pow3(lmMst2)))*pow4(Mst1)*pow4(Mst2)) -
        192*(lmMst2*pow2(lmMst1) + lmMst1*(4 + 2*lmMst2 - 2*pow2(lmMst2)) + (-4
        + lmMst2)*pow2(1 + lmMst2))*pow2(Mst2)*pow6(Mst1) - 3*(-179 + 2*(87 +
        32*lmMst1)*lmMst2 - 30*lmMsq*(-1 + 6*lmMst2) + 90*pow2(lmMsq) + 208*
        pow2(lmMst2))*pow2(Mst1)*pow6(Mst2) + 192*(-((-6 - 14*lmMst2 + lmMst2*
        pow2(lmMst1) + lmMst1*(9 + 6*lmMst2 - 2*pow2(lmMst2)) - 6*pow2(lmMst2)
        + pow3(lmMst2))*pow8(Mst1)) + lmMst2*(1 + lmMst2)*pow8(Mst2))))))/(
        pow4(Msq)*pow4(Mst1)*pow4(Mst2)))/Mgl) + (432*Mgl*(2*Dmglst2*lmMst2 +
        Mgl + lmMst2*Mgl)*s2t*twoLoopFlag*pow2(Mst1)*(-((2*Mt*(MuSUSY*s2t - 2*
        Mt*Tbeta)*pow2(Mst1) + pow2(Mst2)*(2*Mt*MuSUSY*s2t + 4*Tbeta*pow2(Mt) -
        Tbeta*pow2(Mst2)*pow2(s2t)) + Tbeta*pow2(s2t)*pow4(Mst1))*pow4(Mst2)) +
        2*(lmMst1 - lmMst2)*s2t*pow2(Mst1)*(2*Mt*MuSUSY*(pow2(Mst1)*pow2(Mst2)
        + pow4(Mst1) + pow4(Mst2)) - s2t*Tbeta*pow6(Mst2))) + (xDmglst2*pow2(
        Dmglst2)*(432*(2 - 3*lmMst2)*s2t*twoLoopFlag*pow2(Mst1)*pow4(Msq)*(
        Tbeta*(-4*pow2(Mst1)*pow2(Mt) + 4*pow2(Mst2)*pow2(Mt) + 2*(lmMst1 -
        lmMst2)*pow2(Mst1)*pow2(Mst2)*pow2(s2t) + pow2(s2t)*(pow4(Mst1) - pow4(
        Mst2)))*pow4(Mst2) + 2*Mt*MuSUSY*s2t*((pow2(Mst1) + pow2(Mst2))*pow4(
        Mst2) - 2*(lmMst1 - lmMst2)*pow2(Mst1)*(pow2(Mst1)*pow2(Mst2) + pow4(
        Mst1) + pow4(Mst2)))) + Al4p*threeLoopFlag*(6*Mt*MuSUSY*s2t*(s2t*pow2(
        Mst1)*(160*pow2(Msq)*pow2(Mst2)*(5*(8 - 15*lmMsq)*(pow2(Mst1) + pow2(
        Mst2))*pow4(Mst2) + 30*pow2(lmMst2)*(5*pow2(Mst2)*pow4(Mst1) + 5*pow2(
        Mst1)*pow4(Mst2) + 2*pow6(Mst1)) - 2*lmMst1*(5*(8 - 15*lmMsq + 15*
        lmMst2)*pow2(Mst1)*pow2(Mst2)*(pow2(Mst1) + pow2(Mst2)) + (28 - 30*
        lmMsq + 30*lmMst2)*pow6(Mst1)) + lmMst2*(10*(8 - 15*lmMsq)*pow2(Mst2)*
        pow4(Mst1) + 5*(31 - 30*lmMsq)*pow2(Mst1)*pow4(Mst2) + (56 - 60*lmMsq)*
        pow6(Mst1) + 75*pow6(Mst2))) + 45*pow4(Mst2)*(7*(5 - 28*lmMsq)*(pow2(
        Mst1) + pow2(Mst2))*pow4(Mst2) + 56*pow2(lmMst2)*(7*pow2(Mst2)*pow4(
        Mst1) + 7*pow2(Mst1)*pow4(Mst2) + 2*pow6(Mst1)) - 2*lmMst1*(7*(5 - 28*
        lmMsq + 28*lmMst2)*pow2(Mst1)*pow2(Mst2)*(pow2(Mst1) + pow2(Mst2)) + (
        26 - 56*lmMsq + 56*lmMst2)*pow6(Mst1)) + lmMst2*(14*(5 - 28*lmMsq)*
        pow2(Mst2)*pow4(Mst1) + 14*(19 - 28*lmMsq)*pow2(Mst1)*pow4(Mst2) + (52
        - 112*lmMsq)*pow6(Mst1) + 196*pow6(Mst2)))) - 4*pow4(Msq)*(-((Mst2*s2t*
        (3035 + 5240*lmMst2 - 270*lmMsq*(5 + 3*lmMsq - 6*lmMst2)*(1 - 2*lmMst1
        + 2*lmMst2) + 192*(2 - 3*lmMst2)*pow2(lmMst1) - 4620*pow2(lmMst2) + 6*
        lmMst1*(-1161 + 458*lmMst2 + 1008*pow2(lmMst2)) - 5472*pow3(lmMst2)) -
        128*Mt*(-53 - 191*lmMst2 + lmMst1*(60 + 18*lmMst2 - 72*pow2(lmMst2)) +
        54*pow2(lmMst2) + 72*pow3(lmMst2)))*pow3(Mst2)*pow4(Mst1)) - 3*(Mst2*
        s2t*(1065 + 64*lmMst1*(2 - 3*lmMst2) - 266*lmMst2 + 90*lmMsq*(-5 + 6*
        lmMst2) - 270*pow2(lmMsq) - 1264*pow2(lmMst2)) + 512*Mt*(2 + lmMst2 -
        3*pow2(lmMst2)))*pow2(Mst1)*pow5(Mst2) - 2*Mst2*(3*Mst2*s2t*(480 - 9*
        lmMst2*(-129 + 50*lmMsq + 30*pow2(lmMsq)) + 288*(2 - 3*lmMst2)*pow2(
        lmMst1) + 10*(-17 + 54*lmMsq)*pow2(lmMst2) + lmMst1*(-1385 - 406*lmMst2
        - 90*lmMsq*(-5 + 6*lmMst2) + 270*pow2(lmMsq) + 2544*pow2(lmMst2)) -
        1680*pow3(lmMst2)) + 64*Mt*(11 + 371*lmMst2 - 42*pow2(lmMst2) + 6*
        lmMst1*(-21 - 5*lmMst2 + 24*pow2(lmMst2)) - 144*pow3(lmMst2)))*pow6(
        Mst1) + s2t*(12*(96 + lmMst2*(285 + 30*lmMsq + 90*pow2(lmMsq) + 544*
        pow2(lmMst1)) + (590 - 180*lmMsq)*pow2(lmMst2) - lmMst1*(109 + 30*
        lmMsq*(1 - 6*lmMst2) + 590*lmMst2 + 90*pow2(lmMsq) + 1360*pow2(lmMst2))
        + 816*pow3(lmMst2))*pow8(Mst1) - 192*(-2 + lmMst2 + 5*pow2(lmMst2))*
        pow8(Mst2)))) - Tbeta*(256*Mst2*Mt*pow2(Mst1)*pow4(Msq)*(4*pow2(Mt)*(
        pow2(Mst1)*pow2(Mst2)*(96 + 36*lmMt + lmMst2*(-187 + 39*lmMt) - 3*(19 +
        6*lmMt)*pow2(lmMst2) + 18*(lmMst1*(3 - 2*lmMt + lmMst2*(3 + lmMt) -
        pow2(lmMst2)) + pow3(lmMst2))) + (834 + lmMst2*(569 - 33*lmMt) - 162*
        lmMt + (51 - 18*lmMt)*pow2(lmMst2) - 18*lmMst1*(13 + lmMst2 - lmMst2*
        lmMt + pow2(lmMst2)) + 18*pow3(lmMst2))*pow4(Mst1) - 36*(2 + lmMst2 -
        3*pow2(lmMst2))*pow4(Mst2)) + 9*pow2(s2t)*(6*pow2(Mst2)*(-7 + 30*lmMst2
        + 2*pow2(lmMst2) + lmMst1*(-11 - 2*lmMst2 + 12*pow2(lmMst2)) - 12*pow3(
        lmMst2))*pow4(Mst1) + pow2(Mst1)*(29 + 179*lmMst2 - 18*pow2(lmMst2) +
        6*lmMst1*(-10 - 3*lmMst2 + 12*pow2(lmMst2)) - 72*pow3(lmMst2))*pow4(
        Mst2) + 3*(17 + 51*lmMst2 - 4*pow2(lmMst2) + 4*lmMst1*(-6 + lmMst2 + 3*
        pow2(lmMst2)) - 12*pow3(lmMst2))*pow6(Mst1) + 12*(2 + lmMst2 - 3*pow2(
        lmMst2))*pow6(Mst2))) - 3*pow2(Mst2)*pow3(s2t)*(8*pow2(Msq)*(100*(8 -
        15*lmMsq + 15*lmMst2)*pow2(Mst1) + pow2(Msq)*(80 + lmMst2*(-3019 +
        1350*lmMsq + 810*pow2(lmMsq)) + 96*(-2 + 3*lmMst2)*pow2(lmMst1) + 18*(
        23 - 90*lmMsq)*pow2(lmMst2) - 3*lmMst1*(-1225 + 554*lmMst2 - 90*lmMsq*(
        -5 + 6*lmMst2) + 270*pow2(lmMsq) + 1008*pow2(lmMst2)) + 2736*pow3(
        lmMst2)))*pow4(Mst1)*pow4(Mst2) + pow2(Mst1)*(1600*(lmMst1 - lmMst2)*(8
        - 15*lmMsq + 15*lmMst2)*pow2(Msq)*pow2(Mst1) + 12*(-1193 + 330*lmMst2 +
        64*lmMst1*(-2 + 3*lmMst2) - 90*lmMsq*(-5 + 6*lmMst2) + 270*pow2(lmMsq)
        + 1584*pow2(lmMst2))*pow4(Msq) + 315*(5 - 28*lmMsq + 28*lmMst2)*pow4(
        Mst1) + 45*(-9 + 140*lmMsq - 140*lmMst2)*pow4(Mst2))*pow6(Mst2) + pow4(
        Msq)*(4*pow2(Mst2)*(155 - 1726*lmMst2 + 270*lmMsq*(-5 + 6*lmMst2) -
        810*pow2(lmMsq) + 1536*(-2 + 3*lmMst2)*pow2(lmMst1) + 192*lmMst1*(7 +
        27*lmMst2 - 48*pow2(lmMst2)) - 3600*pow2(lmMst2) + 4608*pow3(lmMst2))*
        pow6(Mst1) + 768*(-2*lmMst1*(3 + 2*lmMst2 + 16*pow2(lmMst2)) + lmMst2*(
        7 + 4*lmMst2 + 16*(pow2(lmMst1) + pow2(lmMst2))))*pow8(Mst1)) + 2*(400*
        (-8 + 15*lmMsq - 15*lmMst2)*pow2(Msq)*pow2(Mst1) - 384*(-2 + lmMst2 +
        5*pow2(lmMst2))*pow4(Msq) + 315*(lmMst1 - lmMst2)*(5 - 28*lmMsq + 28*
        lmMst2)*pow4(Mst1))*pow8(Mst2)) + 12*s2t*pow2(Mt)*(pow2(Mst1)*(pow2(
        Mst1) - pow2(Mst2))*(800*(8 - 15*lmMsq + 15*lmMst2)*pow2(Msq)*pow6(
        Mst2) + 315*(5 - 28*lmMsq + 28*lmMst2)*pow8(Mst2)) - 4*pow4(Msq)*(-((
        2395 + 482*lmMst2 + 270*lmMsq*(-5 + 6*lmMst2) - 810*pow2(lmMsq) + 192*(
        -2 + 3*lmMst2)*pow2(lmMst1) + 192*lmMst1*(3 + lmMst2 - 6*pow2(lmMst2))
        - 1296*pow2(lmMst2) + 576*pow3(lmMst2))*pow4(Mst1)*pow4(Mst2)) + 192*
        pow2(Mst2)*((2 - 3*lmMst2)*pow2(lmMst1) + lmMst1*(8 - 2*lmMst2 + 6*
        pow2(lmMst2)) - 3*(6 + 5*lmMst2 + pow3(lmMst2)))*pow6(Mst1) + 3*(873 +
        64*lmMst1*(2 - 3*lmMst2) + 182*lmMst2 + 90*lmMsq*(-5 + 6*lmMst2) - 270*
        pow2(lmMsq) - 624*pow2(lmMst2))*pow2(Mst1)*pow6(Mst2) - 384*(-6 - 14*
        lmMst2 + lmMst2*pow2(lmMst1) + lmMst1*(9 + 6*lmMst2 - 2*pow2(lmMst2)) -
        6*pow2(lmMst2) + pow3(lmMst2))*pow8(Mst1) + 192*(-2 + lmMst2 + 5*pow2(
        lmMst2))*pow8(Mst2)))))))/pow4(Msq))/(Tbeta*pow2(Mgl)*pow4(Mst1)*pow4(
        Mst2))))/1296.);
}

/**
 * @return returns the susy log^0 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double H6g2::calc_coef_at_as2_no_sm_logs_log0() const {

   const double result =
      4832.548148148148 + (7949728*Dmglst2)/(945.*Mgl) + (448*OepS2)/9. + (
        1792*Dmglst2*OepS2)/(27.*Mgl) - (7848*S2)/5. + (76464*Dmglst2*S2)/(35.*
        Mgl) - (3502*Mst2*s2t)/(5.*Mt) + (1322*Dmglst2*Mst2*s2t)/(9.*Mgl*Mt) +
        (448*Mst2*OepS2*s2t)/(9.*Mt) - (2240*Dmglst2*Mst2*OepS2*s2t)/(27.*Mgl*
        Mt) + (5688*Mst2*S2*s2t)/(5.*Mt) - (3240*Dmglst2*Mst2*S2*s2t)/(Mgl*Mt)
        - (224*T1ep)/3. - (896*Dmglst2*T1ep)/(9.*Mgl) - (224*Mst2*s2t*T1ep)/(3.
        *Mt) + (1120*Dmglst2*Mst2*s2t*T1ep)/(9.*Mgl*Mt) + (217532*z2)/45. - (
        1837456*Dmglst2*z2)/(315.*Mgl) - (346772*Mst2*s2t*z2)/(45.*Mt) + (8116*
        Dmglst2*Mst2*s2t*z2)/(3.*Mgl*Mt) - (31216*z3)/3. - (95176*Dmglst2*z3)/(
        9.*Mgl) + (15376*Mst2*s2t*z3)/(3.*Mt) + (5104*Dmglst2*Mst2*s2t*z3)/(
        Mgl*Mt) - (112*z4)/3. - (448*Dmglst2*z4)/(9.*Mgl) - (1648*Mst2*s2t*z4)/
        (3.*Mt) - (4048*Dmglst2*Mst2*s2t*z4)/(9.*Mgl*Mt) + (89346044*pow2(
        Dmglst2))/(2025.*pow2(Mgl)) - (11648*OepS2*pow2(Dmglst2))/(81.*pow2(
        Mgl)) - (30952*S2*pow2(Dmglst2))/(5.*pow2(Mgl)) - (73458722*Mst2*s2t*
        pow2(Dmglst2))/(2835.*Mt*pow2(Mgl)) + (448*Mst2*OepS2*s2t*pow2(Dmglst2)
        )/(81.*Mt*pow2(Mgl)) + (74888*Mst2*S2*s2t*pow2(Dmglst2))/(35.*Mt*pow2(
        Mgl)) + (5824*T1ep*pow2(Dmglst2))/(27.*pow2(Mgl)) - (224*Mst2*s2t*T1ep*
        pow2(Dmglst2))/(27.*Mt*pow2(Mgl)) + (153688*z2*pow2(Dmglst2))/(27.*
        pow2(Mgl)) + (2963804*Mst2*s2t*z2*pow2(Dmglst2))/(945.*Mt*pow2(Mgl)) -
        (903692*z3*pow2(Dmglst2))/(27.*pow2(Mgl)) + (447616*Mst2*s2t*z3*pow2(
        Dmglst2))/(27.*Mt*pow2(Mgl)) + (2912*z4*pow2(Dmglst2))/(27.*pow2(Mgl))
        - (13936*Mst2*s2t*z4*pow2(Dmglst2))/(27.*Mt*pow2(Mgl)) + (160*pow2(Msq)
        )/pow2(Mst1) - (320*z2*pow2(Msq))/pow2(Mst1) - (2117986*s2t*pow2(Mst1))
        /(405.*Mst2*Mt) + (1357642*Dmglst2*s2t*pow2(Mst1))/(135.*Mgl*Mst2*Mt) +
        (24128*OepS2*s2t*pow2(Mst1))/(81.*Mst2*Mt) - (29504*Dmglst2*OepS2*s2t*
        pow2(Mst1))/(27.*Mgl*Mst2*Mt) - (34936*S2*s2t*pow2(Mst1))/(5.*Mst2*Mt)
        + (125528*Dmglst2*S2*s2t*pow2(Mst1))/(15.*Mgl*Mst2*Mt) - (12064*s2t*
        T1ep*pow2(Mst1))/(27.*Mst2*Mt) + (14752*Dmglst2*s2t*T1ep*pow2(Mst1))/(
        9.*Mgl*Mst2*Mt) - (1882628*s2t*z2*pow2(Mst1))/(135.*Mst2*Mt) + (
        1752596*Dmglst2*s2t*z2*pow2(Mst1))/(45.*Mgl*Mst2*Mt) + (387632*s2t*z3*
        pow2(Mst1))/(27.*Mst2*Mt) - (168176*Dmglst2*s2t*z3*pow2(Mst1))/(9.*Mgl*
        Mst2*Mt) - (6032*s2t*z4*pow2(Mst1))/(27.*Mst2*Mt) + (7376*Dmglst2*s2t*
        z4*pow2(Mst1))/(9.*Mgl*Mst2*Mt) - (181427002*s2t*pow2(Dmglst2)*pow2(
        Mst1))/(8505.*Mst2*Mt*pow2(Mgl)) + (275648*OepS2*s2t*pow2(Dmglst2)*
        pow2(Mst1))/(243.*Mst2*Mt*pow2(Mgl)) + (6085624*S2*s2t*pow2(Dmglst2)*
        pow2(Mst1))/(315.*Mst2*Mt*pow2(Mgl)) - (137824*s2t*T1ep*pow2(Dmglst2)*
        pow2(Mst1))/(81.*Mst2*Mt*pow2(Mgl)) - (84331796*s2t*z2*pow2(Dmglst2)*
        pow2(Mst1))/(2835.*Mst2*Mt*pow2(Mgl)) + (545168*s2t*z3*pow2(Dmglst2)*
        pow2(Mst1))/(81.*Mst2*Mt*pow2(Mgl)) - (68912*s2t*z4*pow2(Dmglst2)*pow2(
        Mst1))/(81.*Mst2*Mt*pow2(Mgl)) + (359728*pow2(Mst1))/(675.*pow2(Msq)) -
        (320*Dmglst2*pow2(Mst1))/(3.*Mgl*pow2(Msq)) - (15200*Mst2*s2t*pow2(
        Mst1))/(27.*Mt*pow2(Msq)) - (20960*Dmglst2*Mst2*s2t*pow2(Mst1))/(27.*
        Mgl*Mt*pow2(Msq)) - (480*pow2(Dmglst2)*pow2(Mst1))/(pow2(Mgl)*pow2(Msq)
        ) - (800*Mst2*s2t*pow2(Dmglst2)*pow2(Mst1))/(27.*Mt*pow2(Mgl)*pow2(Msq)
        ) + (160*pow2(Msq))/pow2(Mst2) - (320*z2*pow2(Msq))/pow2(Mst2) + (
        13813838*pow2(Mst1))/(2025.*pow2(Mst2)) + (184384408*Dmglst2*pow2(Mst1)
        )/(70875.*Mgl*pow2(Mst2)) + (13376*OepS2*pow2(Mst1))/(81.*pow2(Mst2)) +
        (41216*Dmglst2*OepS2*pow2(Mst1))/(81.*Mgl*pow2(Mst2)) - (33304*S2*pow2(
        Mst1))/(5.*pow2(Mst2)) - (42320*Dmglst2*S2*pow2(Mst1))/(21.*Mgl*pow2(
        Mst2)) - (6688*T1ep*pow2(Mst1))/(27.*pow2(Mst2)) - (20608*Dmglst2*T1ep*
        pow2(Mst1))/(27.*Mgl*pow2(Mst2)) + (1240156*z2*pow2(Mst1))/(135.*pow2(
        Mst2)) - (42483152*Dmglst2*z2*pow2(Mst1))/(945.*Mgl*pow2(Mst2)) - (
        552400*z3*pow2(Mst1))/(27.*pow2(Mst2)) + (784832*Dmglst2*z3*pow2(Mst1))
        /(27.*Mgl*pow2(Mst2)) - (3344*z4*pow2(Mst1))/(27.*pow2(Mst2)) - (10304*
        Dmglst2*z4*pow2(Mst1))/(27.*Mgl*pow2(Mst2)) + (14329612148*pow2(
        Dmglst2)*pow2(Mst1))/(212625.*pow2(Mgl)*pow2(Mst2)) - (601984*OepS2*
        pow2(Dmglst2)*pow2(Mst1))/(243.*pow2(Mgl)*pow2(Mst2)) - (1030376*S2*
        pow2(Dmglst2)*pow2(Mst1))/(315.*pow2(Mgl)*pow2(Mst2)) + (300992*T1ep*
        pow2(Dmglst2)*pow2(Mst1))/(81.*pow2(Mgl)*pow2(Mst2)) + (50095352*z2*
        pow2(Dmglst2)*pow2(Mst1))/(567.*pow2(Mgl)*pow2(Mst2)) - (5069344*z3*
        pow2(Dmglst2)*pow2(Mst1))/(81.*pow2(Mgl)*pow2(Mst2)) + (150496*z4*pow2(
        Dmglst2)*pow2(Mst1))/(81.*pow2(Mgl)*pow2(Mst2)) - (13936*pow2(Mst2))/(
        25.*pow2(Msq)) - (800*Dmglst2*pow2(Mst2))/(3.*Mgl*pow2(Msq)) + (1600*
        z2*pow2(Mst2))/(3.*pow2(Msq)) + (640*Dmglst2*z2*pow2(Mst2))/(Mgl*pow2(
        Msq)) - (32768*pow2(Dmglst2)*pow2(Mst2))/(45.*pow2(Mgl)*pow2(Msq)) + (
        960*z2*pow2(Dmglst2)*pow2(Mst2))/(pow2(Mgl)*pow2(Msq)) + (2216*pow2(
        Mst2))/(3.*pow2(Mst1)) - (1296*Dmglst2*pow2(Mst2))/(Mgl*pow2(Mst1)) + (
        1568*z2*pow2(Mst2))/(3.*pow2(Mst1)) + (2944*Dmglst2*z2*pow2(Mst2))/(3.*
        Mgl*pow2(Mst1)) - (8168*pow2(Dmglst2)*pow2(Mst2))/(3.*pow2(Mgl)*pow2(
        Mst1)) + (1536*z2*pow2(Dmglst2)*pow2(Mst2))/(pow2(Mgl)*pow2(Mst1)) + (
        160*pow2(Msq)*pow2(s2t))/pow2(Mt) - (320*z2*pow2(Msq)*pow2(s2t))/pow2(
        Mt) + (4141616*pow2(Mst1)*pow2(s2t))/(2025.*pow2(Mt)) - (32*D3*pow2(
        Mst1)*pow2(s2t))/(3.*pow2(Mt)) + (32*DN*pow2(Mst1)*pow2(s2t))/(3.*pow2(
        Mt)) + (19333928*Dmglst2*pow2(Mst1)*pow2(s2t))/(10125.*Mgl*pow2(Mt)) -
        (20608*OepS2*pow2(Mst1)*pow2(s2t))/(81.*pow2(Mt)) + (6400*Dmglst2*
        OepS2*pow2(Mst1)*pow2(s2t))/(27.*Mgl*pow2(Mt)) + (38264*S2*pow2(Mst1)*
        pow2(s2t))/(5.*pow2(Mt)) - (2976*Dmglst2*S2*pow2(Mst1)*pow2(s2t))/(5.*
        Mgl*pow2(Mt)) + (10304*T1ep*pow2(Mst1)*pow2(s2t))/(27.*pow2(Mt)) - (
        3200*Dmglst2*T1ep*pow2(Mst1)*pow2(s2t))/(9.*Mgl*pow2(Mt)) + (591448*z2*
        pow2(Mst1)*pow2(s2t))/(135.*pow2(Mt)) - (17600*Dmglst2*z2*pow2(Mst1)*
        pow2(s2t))/(3.*Mgl*pow2(Mt)) - (15040*z3*pow2(Mst1)*pow2(s2t))/(27.*
        pow2(Mt)) - (3952*Dmglst2*z3*pow2(Mst1)*pow2(s2t))/(3.*Mgl*pow2(Mt)) +
        (3424*z4*pow2(Mst1)*pow2(s2t))/(27.*pow2(Mt)) - (1600*Dmglst2*z4*pow2(
        Mst1)*pow2(s2t))/(9.*Mgl*pow2(Mt)) + (281425742*pow2(Dmglst2)*pow2(
        Mst1)*pow2(s2t))/(70875.*pow2(Mgl)*pow2(Mt)) + (12736*OepS2*pow2(
        Dmglst2)*pow2(Mst1)*pow2(s2t))/(27.*pow2(Mgl)*pow2(Mt)) - (339064*S2*
        pow2(Dmglst2)*pow2(Mst1)*pow2(s2t))/(35.*pow2(Mgl)*pow2(Mt)) - (6368*
        T1ep*pow2(Dmglst2)*pow2(Mst1)*pow2(s2t))/(9.*pow2(Mgl)*pow2(Mt)) - (
        1010572*z2*pow2(Dmglst2)*pow2(Mst1)*pow2(s2t))/(63.*pow2(Mgl)*pow2(Mt))
        + (18250*z3*pow2(Dmglst2)*pow2(Mst1)*pow2(s2t))/(9.*pow2(Mgl)*pow2(Mt))
        - (3184*z4*pow2(Dmglst2)*pow2(Mst1)*pow2(s2t))/(9.*pow2(Mgl)*pow2(Mt))
        - (80*pow2(Msq)*pow2(Mst1)*pow2(s2t))/(pow2(Mst2)*pow2(Mt)) + (160*z2*
        pow2(Msq)*pow2(Mst1)*pow2(s2t))/(pow2(Mst2)*pow2(Mt)) + (818213*pow2(
        Mst2)*pow2(s2t))/(270.*pow2(Mt)) + (32*D3*pow2(Mst2)*pow2(s2t))/(3.*
        pow2(Mt)) - (32*DN*pow2(Mst2)*pow2(s2t))/(3.*pow2(Mt)) + (226984*
        Dmglst2*pow2(Mst2)*pow2(s2t))/(45.*Mgl*pow2(Mt)) - (128*B4*Dmglst2*
        pow2(Mst2)*pow2(s2t))/(3.*Mgl*pow2(Mt)) + (128*D3*Dmglst2*pow2(Mst2)*
        pow2(s2t))/(3.*Mgl*pow2(Mt)) - (64*Dmglst2*DN*pow2(Mst2)*pow2(s2t))/(3.
        *Mgl*pow2(Mt)) - (368*OepS2*pow2(Mst2)*pow2(s2t))/(9.*pow2(Mt)) - (
        16038*S2*pow2(Mst2)*pow2(s2t))/(5.*pow2(Mt)) - (30384*Dmglst2*S2*pow2(
        Mst2)*pow2(s2t))/(5.*Mgl*pow2(Mt)) + (184*T1ep*pow2(Mst2)*pow2(s2t))/(
        3.*pow2(Mt)) - (17923*z2*pow2(Mst2)*pow2(s2t))/(45.*pow2(Mt)) - (55984*
        Dmglst2*z2*pow2(Mst2)*pow2(s2t))/(45.*Mgl*pow2(Mt)) + (3244*z3*pow2(
        Mst2)*pow2(s2t))/(3.*pow2(Mt)) - (7024*Dmglst2*z3*pow2(Mst2)*pow2(s2t))
        /(9.*Mgl*pow2(Mt)) + (284*z4*pow2(Mst2)*pow2(s2t))/(3.*pow2(Mt)) + (
        192*Dmglst2*z4*pow2(Mst2)*pow2(s2t))/(Mgl*pow2(Mt)) + (1798498*pow2(
        Dmglst2)*pow2(Mst2)*pow2(s2t))/(105.*pow2(Mgl)*pow2(Mt)) - (64*B4*pow2(
        Dmglst2)*pow2(Mst2)*pow2(s2t))/(pow2(Mgl)*pow2(Mt)) + (64*D3*pow2(
        Dmglst2)*pow2(Mst2)*pow2(s2t))/(pow2(Mgl)*pow2(Mt)) - (32*DN*pow2(
        Dmglst2)*pow2(Mst2)*pow2(s2t))/(pow2(Mgl)*pow2(Mt)) + (448*OepS2*pow2(
        Dmglst2)*pow2(Mst2)*pow2(s2t))/(9.*pow2(Mgl)*pow2(Mt)) - (250344*S2*
        pow2(Dmglst2)*pow2(Mst2)*pow2(s2t))/(35.*pow2(Mgl)*pow2(Mt)) - (224*
        T1ep*pow2(Dmglst2)*pow2(Mst2)*pow2(s2t))/(3.*pow2(Mgl)*pow2(Mt)) - (
        215708*z2*pow2(Dmglst2)*pow2(Mst2)*pow2(s2t))/(105.*pow2(Mgl)*pow2(Mt))
        - (39406*z3*pow2(Dmglst2)*pow2(Mst2)*pow2(s2t))/(3.*pow2(Mgl)*pow2(Mt))
        + (752*z4*pow2(Dmglst2)*pow2(Mst2)*pow2(s2t))/(3.*pow2(Mgl)*pow2(Mt)) -
        (80*pow2(Msq)*pow2(Mst2)*pow2(s2t))/(pow2(Mst1)*pow2(Mt)) + (160*z2*
        pow2(Msq)*pow2(Mst2)*pow2(s2t))/(pow2(Mst1)*pow2(Mt)) - (24728*pow2(
        Mst1)*pow2(Mst2)*pow2(s2t))/(675.*pow2(Msq)*pow2(Mt)) + (400*Dmglst2*
        pow2(Mst1)*pow2(Mst2)*pow2(s2t))/(Mgl*pow2(Msq)*pow2(Mt)) + (1960*pow2(
        Dmglst2)*pow2(Mst1)*pow2(Mst2)*pow2(s2t))/(3.*pow2(Mgl)*pow2(Msq)*pow2(
        Mt)) - 56*pow2(z2) - (224*Dmglst2*pow2(z2))/(3.*Mgl) - (56*Mst2*s2t*
        pow2(z2))/Mt + (280*Dmglst2*Mst2*s2t*pow2(z2))/(3.*Mgl*Mt) + (1456*
        pow2(Dmglst2)*pow2(z2))/(9.*pow2(Mgl)) - (56*Mst2*s2t*pow2(Dmglst2)*
        pow2(z2))/(9.*Mt*pow2(Mgl)) - (3016*s2t*pow2(Mst1)*pow2(z2))/(9.*Mst2*
        Mt) + (3688*Dmglst2*s2t*pow2(Mst1)*pow2(z2))/(3.*Mgl*Mst2*Mt) - (34456*
        s2t*pow2(Dmglst2)*pow2(Mst1)*pow2(z2))/(27.*Mst2*Mt*pow2(Mgl)) - (1672*
        pow2(Mst1)*pow2(z2))/(9.*pow2(Mst2)) - (5152*Dmglst2*pow2(Mst1)*pow2(
        z2))/(9.*Mgl*pow2(Mst2)) + (75248*pow2(Dmglst2)*pow2(Mst1)*pow2(z2))/(
        27.*pow2(Mgl)*pow2(Mst2)) + (2576*pow2(Mst1)*pow2(s2t)*pow2(z2))/(9.*
        pow2(Mt)) - (800*Dmglst2*pow2(Mst1)*pow2(s2t)*pow2(z2))/(3.*Mgl*pow2(
        Mt)) - (1592*pow2(Dmglst2)*pow2(Mst1)*pow2(s2t)*pow2(z2))/(3.*pow2(Mgl)
        *pow2(Mt)) + (46*pow2(Mst2)*pow2(s2t)*pow2(z2))/pow2(Mt) - (56*pow2(
        Dmglst2)*pow2(Mst2)*pow2(s2t)*pow2(z2))/(pow2(Mgl)*pow2(Mt)) - (38080*
        s2t*pow3(Mst2))/(27.*Mt*pow2(Msq)) - (40000*Dmglst2*s2t*pow3(Mst2))/(9.
        *Mgl*Mt*pow2(Msq)) + (1280*s2t*z2*pow3(Mst2))/(Mt*pow2(Msq)) + (8960*
        Dmglst2*s2t*z2*pow3(Mst2))/(3.*Mgl*Mt*pow2(Msq)) - (82880*s2t*pow2(
        Dmglst2)*pow3(Mst2))/(9.*Mt*pow2(Mgl)*pow2(Msq)) + (16640*s2t*z2*pow2(
        Dmglst2)*pow3(Mst2))/(3.*Mt*pow2(Mgl)*pow2(Msq)) - (1024*s2t*pow3(Mst2)
        )/(3.*Mt*pow2(Mst1)) + (1024*Dmglst2*s2t*pow3(Mst2))/(3.*Mgl*Mt*pow2(
        Mst1)) + (4096*s2t*pow2(Dmglst2)*pow3(Mst2))/(3.*Mt*pow2(Mgl)*pow2(
        Mst1)) + (18641*Mst2*pow2(Mst1)*pow3(s2t))/(81.*pow3(Mt)) + (272*B4*
        Mst2*pow2(Mst1)*pow3(s2t))/(3.*pow3(Mt)) + (8*DN*Mst2*pow2(Mst1)*pow3(
        s2t))/(3.*pow3(Mt)) - (231989*Dmglst2*Mst2*pow2(Mst1)*pow3(s2t))/(405.*
        Mgl*pow3(Mt)) + (272*B4*Dmglst2*Mst2*pow2(Mst1)*pow3(s2t))/(3.*Mgl*
        pow3(Mt)) + (8*Dmglst2*DN*Mst2*pow2(Mst1)*pow3(s2t))/(3.*Mgl*pow3(Mt))
        - (2336*Mst2*OepS2*pow2(Mst1)*pow3(s2t))/(81.*pow3(Mt)) + (12832*
        Dmglst2*Mst2*OepS2*pow2(Mst1)*pow3(s2t))/(81.*Mgl*pow3(Mt)) + (1388*
        Mst2*S2*pow2(Mst1)*pow3(s2t))/pow3(Mt) - (56132*Dmglst2*Mst2*S2*pow2(
        Mst1)*pow3(s2t))/(15.*Mgl*pow3(Mt)) + (1168*Mst2*T1ep*pow2(Mst1)*pow3(
        s2t))/(27.*pow3(Mt)) - (6416*Dmglst2*Mst2*T1ep*pow2(Mst1)*pow3(s2t))/(
        27.*Mgl*pow3(Mt)) + (15826*Mst2*z2*pow2(Mst1)*pow3(s2t))/(27.*pow3(Mt))
        - (467482*Dmglst2*Mst2*z2*pow2(Mst1)*pow3(s2t))/(135.*Mgl*pow3(Mt)) - (
        34712*Mst2*z3*pow2(Mst1)*pow3(s2t))/(27.*pow3(Mt)) + (25000*Dmglst2*
        Mst2*z3*pow2(Mst1)*pow3(s2t))/(27.*Mgl*pow3(Mt)) + (2960*Mst2*z4*pow2(
        Mst1)*pow3(s2t))/(27.*pow3(Mt)) - (832*Dmglst2*Mst2*z4*pow2(Mst1)*pow3(
        s2t))/(27.*Mgl*pow3(Mt)) + (1711970*Mst2*pow2(Dmglst2)*pow2(Mst1)*pow3(
        s2t))/(243.*pow2(Mgl)*pow3(Mt)) + (272*B4*Mst2*pow2(Dmglst2)*pow2(Mst1)
        *pow3(s2t))/(3.*pow2(Mgl)*pow3(Mt)) + (8*DN*Mst2*pow2(Dmglst2)*pow2(
        Mst1)*pow3(s2t))/(3.*pow2(Mgl)*pow3(Mt)) - (49088*Mst2*OepS2*pow2(
        Dmglst2)*pow2(Mst1)*pow3(s2t))/(243.*pow2(Mgl)*pow3(Mt)) - (6656*Mst2*
        S2*pow2(Dmglst2)*pow2(Mst1)*pow3(s2t))/(9.*pow2(Mgl)*pow3(Mt)) + (
        24544*Mst2*T1ep*pow2(Dmglst2)*pow2(Mst1)*pow3(s2t))/(81.*pow2(Mgl)*
        pow3(Mt)) + (306340*Mst2*z2*pow2(Dmglst2)*pow2(Mst1)*pow3(s2t))/(81.*
        pow2(Mgl)*pow3(Mt)) - (280520*Mst2*z3*pow2(Dmglst2)*pow2(Mst1)*pow3(
        s2t))/(81.*pow2(Mgl)*pow3(Mt)) + (19400*Mst2*z4*pow2(Dmglst2)*pow2(
        Mst1)*pow3(s2t))/(81.*pow2(Mgl)*pow3(Mt)) + (292*Mst2*pow2(Mst1)*pow2(
        z2)*pow3(s2t))/(9.*pow3(Mt)) - (1604*Dmglst2*Mst2*pow2(Mst1)*pow2(z2)*
        pow3(s2t))/(9.*Mgl*pow3(Mt)) + (6136*Mst2*pow2(Dmglst2)*pow2(Mst1)*
        pow2(z2)*pow3(s2t))/(27.*pow2(Mgl)*pow3(Mt)) + (17177*pow3(Mst2)*pow3(
        s2t))/(6.*pow3(Mt)) + (880*B4*pow3(Mst2)*pow3(s2t))/(3.*pow3(Mt)) - (8*
        DN*pow3(Mst2)*pow3(s2t))/(3.*pow3(Mt)) + (825997*Dmglst2*pow3(Mst2)*
        pow3(s2t))/(270.*Mgl*pow3(Mt)) + (2096*B4*Dmglst2*pow3(Mst2)*pow3(s2t))
        /(3.*Mgl*pow3(Mt)) - (40*Dmglst2*DN*pow3(Mst2)*pow3(s2t))/(3.*Mgl*pow3(
        Mt)) - (112*OepS2*pow3(Mst2)*pow3(s2t))/(9.*pow3(Mt)) + (560*Dmglst2*
        OepS2*pow3(Mst2)*pow3(s2t))/(27.*Mgl*pow3(Mt)) - (18*S2*pow3(Mst2)*
        pow3(s2t))/pow3(Mt) + (2718*Dmglst2*S2*pow3(Mst2)*pow3(s2t))/(5.*Mgl*
        pow3(Mt)) + (56*T1ep*pow3(Mst2)*pow3(s2t))/(3.*pow3(Mt)) - (280*
        Dmglst2*T1ep*pow3(Mst2)*pow3(s2t))/(9.*Mgl*pow3(Mt)) - (3469*z2*pow3(
        Mst2)*pow3(s2t))/(3.*pow3(Mt)) - (278347*Dmglst2*z2*pow3(Mst2)*pow3(
        s2t))/(45.*Mgl*pow3(Mt)) + (20980*z3*pow3(Mst2)*pow3(s2t))/(9.*pow3(Mt)
        ) + (74948*Dmglst2*z3*pow3(Mst2)*pow3(s2t))/(9.*Mgl*pow3(Mt)) + (340*
        z4*pow3(Mst2)*pow3(s2t))/(3.*pow3(Mt)) + (1084*Dmglst2*z4*pow3(Mst2)*
        pow3(s2t))/(9.*Mgl*pow3(Mt)) + (986407*pow2(Dmglst2)*pow3(Mst2)*pow3(
        s2t))/(810.*pow2(Mgl)*pow3(Mt)) + (3920*B4*pow2(Dmglst2)*pow3(Mst2)*
        pow3(s2t))/(3.*pow2(Mgl)*pow3(Mt)) - (88*DN*pow2(Dmglst2)*pow3(Mst2)*
        pow3(s2t))/(3.*pow2(Mgl)*pow3(Mt)) - (112*OepS2*pow2(Dmglst2)*pow3(
        Mst2)*pow3(s2t))/(81.*pow2(Mgl)*pow3(Mt)) - (2774*S2*pow2(Dmglst2)*
        pow3(Mst2)*pow3(s2t))/(5.*pow2(Mgl)*pow3(Mt)) + (56*T1ep*pow2(Dmglst2)*
        pow3(Mst2)*pow3(s2t))/(27.*pow2(Mgl)*pow3(Mt)) - (1340537*z2*pow2(
        Dmglst2)*pow3(Mst2)*pow3(s2t))/(135.*pow2(Mgl)*pow3(Mt)) + (512948*z3*
        pow2(Dmglst2)*pow3(Mst2)*pow3(s2t))/(27.*pow2(Mgl)*pow3(Mt)) + (4996*
        z4*pow2(Dmglst2)*pow3(Mst2)*pow3(s2t))/(27.*pow2(Mgl)*pow3(Mt)) + (
        4000*pow2(Mst1)*pow3(Mst2)*pow3(s2t))/(9.*pow2(Msq)*pow3(Mt)) + (960*
        Dmglst2*pow2(Mst1)*pow3(Mst2)*pow3(s2t))/(Mgl*pow2(Msq)*pow3(Mt)) - (
        320*z2*pow2(Mst1)*pow3(Mst2)*pow3(s2t))/(3.*pow2(Msq)*pow3(Mt)) - (320*
        Dmglst2*z2*pow2(Mst1)*pow3(Mst2)*pow3(s2t))/(Mgl*pow2(Msq)*pow3(Mt)) +
        (11440*pow2(Dmglst2)*pow2(Mst1)*pow3(Mst2)*pow3(s2t))/(9.*pow2(Mgl)*
        pow2(Msq)*pow3(Mt)) - (640*z2*pow2(Dmglst2)*pow2(Mst1)*pow3(Mst2)*pow3(
        s2t))/(pow2(Mgl)*pow2(Msq)*pow3(Mt)) - (370*pow2(z2)*pow3(Mst2)*pow3(
        s2t))/pow3(Mt) - (3526*Dmglst2*pow2(z2)*pow3(Mst2)*pow3(s2t))/(3.*Mgl*
        pow3(Mt)) - (20722*pow2(Dmglst2)*pow2(z2)*pow3(Mst2)*pow3(s2t))/(9.*
        pow2(Mgl)*pow3(Mt)) + (160*pow3(log(pow2(Mst1)/pow2(Msq))))/3. + (100*
        pow2(Mst1)*pow2(Mst2))/(9.*pow4(Msq)) + (3320*Dmglst2*pow2(Mst1)*pow2(
        Mst2))/(9.*Mgl*pow4(Msq)) + (580*pow2(Dmglst2)*pow2(Mst1)*pow2(Mst2))/(
        pow2(Mgl)*pow4(Msq)) - (400*s2t*pow2(Mst1)*pow3(Mst2))/(27.*Mt*pow4(
        Msq)) - (3280*Dmglst2*s2t*pow2(Mst1)*pow3(Mst2))/(9.*Mgl*Mt*pow4(Msq))
        - (7040*s2t*pow2(Dmglst2)*pow2(Mst1)*pow3(Mst2))/(9.*Mt*pow2(Mgl)*pow4(
        Msq)) - (2240*s2t*pow4(Mst1))/(3.*Mst2*Mt*pow2(Msq)) + (1600*Dmglst2*
        s2t*pow4(Mst1))/(Mgl*Mst2*Mt*pow2(Msq)) + (1600*s2t*pow2(Dmglst2)*pow4(
        Mst1))/(Mst2*Mt*pow2(Mgl)*pow2(Msq)) + (1040*pow4(Mst1))/(pow2(Msq)*
        pow2(Mst2)) - (4640*Dmglst2*pow4(Mst1))/(Mgl*pow2(Msq)*pow2(Mst2)) + (
        320*z2*pow4(Mst1))/(pow2(Msq)*pow2(Mst2)) - (640*Dmglst2*z2*pow4(Mst1))
        /(Mgl*pow2(Msq)*pow2(Mst2)) - (4640*pow2(Dmglst2)*pow4(Mst1))/(pow2(
        Mgl)*pow2(Msq)*pow2(Mst2)) - (640*z2*pow2(Dmglst2)*pow4(Mst1))/(pow2(
        Mgl)*pow2(Msq)*pow2(Mst2)) + (40864*pow2(s2t)*pow4(Mst1))/(675.*pow2(
        Msq)*pow2(Mt)) + (640*Dmglst2*pow2(s2t)*pow4(Mst1))/(3.*Mgl*pow2(Msq)*
        pow2(Mt)) + (640*pow2(Dmglst2)*pow2(s2t)*pow4(Mst1))/(3.*pow2(Mgl)*
        pow2(Msq)*pow2(Mt)) + (24574447162*pow2(s2t)*pow4(Mst1))/(1.0418625e7*
        pow2(Mst2)*pow2(Mt)) - (7932854794*Dmglst2*pow2(s2t)*pow4(Mst1))/(
        3.472875e6*Mgl*pow2(Mst2)*pow2(Mt)) - (110672*OepS2*pow2(s2t)*pow4(
        Mst1))/(243.*pow2(Mst2)*pow2(Mt)) + (98432*Dmglst2*OepS2*pow2(s2t)*
        pow4(Mst1))/(81.*Mgl*pow2(Mst2)*pow2(Mt)) + (813254*S2*pow2(s2t)*pow4(
        Mst1))/(45.*pow2(Mst2)*pow2(Mt)) - (452528*Dmglst2*S2*pow2(s2t)*pow4(
        Mst1))/(15.*Mgl*pow2(Mst2)*pow2(Mt)) + (55336*T1ep*pow2(s2t)*pow4(Mst1)
        )/(81.*pow2(Mst2)*pow2(Mt)) - (49216*Dmglst2*T1ep*pow2(s2t)*pow4(Mst1))
        /(27.*Mgl*pow2(Mst2)*pow2(Mt)) + (3321329*z2*pow2(s2t)*pow4(Mst1))/(
        405.*pow2(Mst2)*pow2(Mt)) - (628264*Dmglst2*z2*pow2(s2t)*pow4(Mst1))/(
        27.*Mgl*pow2(Mst2)*pow2(Mt)) - (132332*z3*pow2(s2t)*pow4(Mst1))/(81.*
        pow2(Mst2)*pow2(Mt)) + (11744*Dmglst2*z3*pow2(s2t)*pow4(Mst1))/(27.*
        Mgl*pow2(Mst2)*pow2(Mt)) + (27668*z4*pow2(s2t)*pow4(Mst1))/(81.*pow2(
        Mst2)*pow2(Mt)) - (24608*Dmglst2*z4*pow2(s2t)*pow4(Mst1))/(27.*Mgl*
        pow2(Mst2)*pow2(Mt)) - (7932854794*pow2(Dmglst2)*pow2(s2t)*pow4(Mst1))/
        (3.472875e6*pow2(Mgl)*pow2(Mst2)*pow2(Mt)) + (98432*OepS2*pow2(Dmglst2)
        *pow2(s2t)*pow4(Mst1))/(81.*pow2(Mgl)*pow2(Mst2)*pow2(Mt)) - (452528*
        S2*pow2(Dmglst2)*pow2(s2t)*pow4(Mst1))/(15.*pow2(Mgl)*pow2(Mst2)*pow2(
        Mt)) - (49216*T1ep*pow2(Dmglst2)*pow2(s2t)*pow4(Mst1))/(27.*pow2(Mgl)*
        pow2(Mst2)*pow2(Mt)) - (628264*z2*pow2(Dmglst2)*pow2(s2t)*pow4(Mst1))/(
        27.*pow2(Mgl)*pow2(Mst2)*pow2(Mt)) + (11744*z3*pow2(Dmglst2)*pow2(s2t)*
        pow4(Mst1))/(27.*pow2(Mgl)*pow2(Mst2)*pow2(Mt)) - (24608*z4*pow2(
        Dmglst2)*pow2(s2t)*pow4(Mst1))/(27.*pow2(Mgl)*pow2(Mst2)*pow2(Mt)) + (
        13834*pow2(s2t)*pow2(z2)*pow4(Mst1))/(27.*pow2(Mst2)*pow2(Mt)) - (
        12304*Dmglst2*pow2(s2t)*pow2(z2)*pow4(Mst1))/(9.*Mgl*pow2(Mst2)*pow2(
        Mt)) - (12304*pow2(Dmglst2)*pow2(s2t)*pow2(z2)*pow4(Mst1))/(9.*pow2(
        Mgl)*pow2(Mst2)*pow2(Mt)) - (15260456*s2t*pow4(Mst1))/(1215.*Mt*pow3(
        Mst2)) + (193667576*Dmglst2*s2t*pow4(Mst1))/(3645.*Mgl*Mt*pow3(Mst2)) +
        (175360*OepS2*s2t*pow4(Mst1))/(243.*Mt*pow3(Mst2)) - (3293440*Dmglst2*
        OepS2*s2t*pow4(Mst1))/(729.*Mgl*Mt*pow3(Mst2)) - (1043008*S2*s2t*pow4(
        Mst1))/(45.*Mt*pow3(Mst2)) + (14087488*Dmglst2*S2*s2t*pow4(Mst1))/(135.
        *Mgl*Mt*pow3(Mst2)) - (87680*s2t*T1ep*pow4(Mst1))/(81.*Mt*pow3(Mst2)) +
        (1646720*Dmglst2*s2t*T1ep*pow4(Mst1))/(243.*Mgl*Mt*pow3(Mst2)) - (
        9196144*s2t*z2*pow4(Mst1))/(405.*Mt*pow3(Mst2)) + (147644272*Dmglst2*
        s2t*z2*pow4(Mst1))/(1215.*Mgl*Mt*pow3(Mst2)) + (1848256*s2t*z3*pow4(
        Mst1))/(81.*Mt*pow3(Mst2)) - (19670464*Dmglst2*s2t*z3*pow4(Mst1))/(243.
        *Mgl*Mt*pow3(Mst2)) - (43840*s2t*z4*pow4(Mst1))/(81.*Mt*pow3(Mst2)) + (
        823360*Dmglst2*s2t*z4*pow4(Mst1))/(243.*Mgl*Mt*pow3(Mst2)) + (
        193667576*s2t*pow2(Dmglst2)*pow4(Mst1))/(3645.*Mt*pow2(Mgl)*pow3(Mst2))
        - (3293440*OepS2*s2t*pow2(Dmglst2)*pow4(Mst1))/(729.*Mt*pow2(Mgl)*pow3(
        Mst2)) + (14087488*S2*s2t*pow2(Dmglst2)*pow4(Mst1))/(135.*Mt*pow2(Mgl)*
        pow3(Mst2)) + (1646720*s2t*T1ep*pow2(Dmglst2)*pow4(Mst1))/(243.*Mt*
        pow2(Mgl)*pow3(Mst2)) + (147644272*s2t*z2*pow2(Dmglst2)*pow4(Mst1))/(
        1215.*Mt*pow2(Mgl)*pow3(Mst2)) - (19670464*s2t*z3*pow2(Dmglst2)*pow4(
        Mst1))/(243.*Mt*pow2(Mgl)*pow3(Mst2)) + (823360*s2t*z4*pow2(Dmglst2)*
        pow4(Mst1))/(243.*Mt*pow2(Mgl)*pow3(Mst2)) - (21920*s2t*pow2(z2)*pow4(
        Mst1))/(27.*Mt*pow3(Mst2)) + (411680*Dmglst2*s2t*pow2(z2)*pow4(Mst1))/(
        81.*Mgl*Mt*pow3(Mst2)) + (411680*s2t*pow2(Dmglst2)*pow2(z2)*pow4(Mst1))
        /(81.*Mt*pow2(Mgl)*pow3(Mst2)) + (542267*pow3(s2t)*pow4(Mst1))/(486.*
        Mst2*pow3(Mt)) - (35542243*Dmglst2*pow3(s2t)*pow4(Mst1))/(7290.*Mgl*
        Mst2*pow3(Mt)) - (7024*OepS2*pow3(s2t)*pow4(Mst1))/(243.*Mst2*pow3(Mt))
        + (249136*Dmglst2*OepS2*pow3(s2t)*pow4(Mst1))/(729.*Mgl*Mst2*pow3(Mt))
        + (12242*S2*pow3(s2t)*pow4(Mst1))/(9.*Mst2*pow3(Mt)) - (1672042*
        Dmglst2*S2*pow3(s2t)*pow4(Mst1))/(135.*Mgl*Mst2*pow3(Mt)) + (3512*T1ep*
        pow3(s2t)*pow4(Mst1))/(81.*Mst2*pow3(Mt)) - (124568*Dmglst2*T1ep*pow3(
        s2t)*pow4(Mst1))/(243.*Mgl*Mst2*pow3(Mt)) + (104839*z2*pow3(s2t)*pow4(
        Mst1))/(405.*Mst2*pow3(Mt)) - (8194483*Dmglst2*z2*pow3(s2t)*pow4(Mst1))
        /(1215.*Mgl*Mst2*pow3(Mt)) - (89260*z3*pow3(s2t)*pow4(Mst1))/(81.*Mst2*
        pow3(Mt)) + (1130428*Dmglst2*z3*pow3(s2t)*pow4(Mst1))/(243.*Mgl*Mst2*
        pow3(Mt)) + (1756*z4*pow3(s2t)*pow4(Mst1))/(81.*Mst2*pow3(Mt)) - (
        62284*Dmglst2*z4*pow3(s2t)*pow4(Mst1))/(243.*Mgl*Mst2*pow3(Mt)) - (
        35542243*pow2(Dmglst2)*pow3(s2t)*pow4(Mst1))/(7290.*Mst2*pow2(Mgl)*
        pow3(Mt)) + (249136*OepS2*pow2(Dmglst2)*pow3(s2t)*pow4(Mst1))/(729.*
        Mst2*pow2(Mgl)*pow3(Mt)) - (1672042*S2*pow2(Dmglst2)*pow3(s2t)*pow4(
        Mst1))/(135.*Mst2*pow2(Mgl)*pow3(Mt)) - (124568*T1ep*pow2(Dmglst2)*
        pow3(s2t)*pow4(Mst1))/(243.*Mst2*pow2(Mgl)*pow3(Mt)) - (8194483*z2*
        pow2(Dmglst2)*pow3(s2t)*pow4(Mst1))/(1215.*Mst2*pow2(Mgl)*pow3(Mt)) + (
        1130428*z3*pow2(Dmglst2)*pow3(s2t)*pow4(Mst1))/(243.*Mst2*pow2(Mgl)*
        pow3(Mt)) - (62284*z4*pow2(Dmglst2)*pow3(s2t)*pow4(Mst1))/(243.*Mst2*
        pow2(Mgl)*pow3(Mt)) + (80*Mst2*pow3(s2t)*pow4(Mst1))/(3.*pow2(Msq)*
        pow3(Mt)) - (80*Dmglst2*Mst2*pow3(s2t)*pow4(Mst1))/(3.*Mgl*pow2(Msq)*
        pow3(Mt)) - (80*Mst2*pow2(Dmglst2)*pow3(s2t)*pow4(Mst1))/(3.*pow2(Mgl)*
        pow2(Msq)*pow3(Mt)) + (878*pow2(z2)*pow3(s2t)*pow4(Mst1))/(27.*Mst2*
        pow3(Mt)) - (31142*Dmglst2*pow2(z2)*pow3(s2t)*pow4(Mst1))/(81.*Mgl*
        Mst2*pow3(Mt)) - (31142*pow2(Dmglst2)*pow2(z2)*pow3(s2t)*pow4(Mst1))/(
        81.*Mst2*pow2(Mgl)*pow3(Mt)) + (15728921*pow4(Mst1))/(77175.*pow4(Msq))
        - (1880*Mst2*s2t*pow4(Mst1))/(27.*Mt*pow4(Msq)) - (6200*Dmglst2*Mst2*
        s2t*pow4(Mst1))/(27.*Mgl*Mt*pow4(Msq)) - (6200*Mst2*s2t*pow2(Dmglst2)*
        pow4(Mst1))/(27.*Mt*pow2(Mgl)*pow4(Msq)) - (438008*pow2(Mst2)*pow2(s2t)
        *pow4(Mst1))/(15435.*pow2(Mt)*pow4(Msq)) + (10*Dmglst2*pow2(Mst2)*pow2(
        s2t)*pow4(Mst1))/(Mgl*pow2(Mt)*pow4(Msq)) + (10*pow2(Dmglst2)*pow2(
        Mst2)*pow2(s2t)*pow4(Mst1))/(pow2(Mgl)*pow2(Mt)*pow4(Msq)) - (5*pow3(
        Mst2)*pow3(s2t)*pow4(Mst1))/(9.*pow3(Mt)*pow4(Msq)) + (115*Dmglst2*
        pow3(Mst2)*pow3(s2t)*pow4(Mst1))/(9.*Mgl*pow3(Mt)*pow4(Msq)) + (115*
        pow2(Dmglst2)*pow3(Mst2)*pow3(s2t)*pow4(Mst1))/(9.*pow2(Mgl)*pow3(Mt)*
        pow4(Msq)) + (175715865103*pow4(Mst1))/(1.0418625e7*pow4(Mst2)) - (
        291022313156*Dmglst2*pow4(Mst1))/(6.251175e6*Mgl*pow4(Mst2)) + (74752*
        OepS2*pow4(Mst1))/(243.*pow4(Mst2)) + (1428992*Dmglst2*OepS2*pow4(Mst1)
        )/(729.*Mgl*pow4(Mst2)) - (553792*S2*pow4(Mst1))/(45.*pow4(Mst2)) - (
        32019872*Dmglst2*S2*pow4(Mst1))/(945.*Mgl*pow4(Mst2)) - (37376*T1ep*
        pow4(Mst1))/(81.*pow4(Mst2)) - (714496*Dmglst2*T1ep*pow4(Mst1))/(243.*
        Mgl*pow4(Mst2)) + (4885352*z2*pow4(Mst1))/(405.*pow4(Mst2)) - (
        1089695168*Dmglst2*z2*pow4(Mst1))/(8505.*Mgl*pow4(Mst2)) - (3102848*z3*
        pow4(Mst1))/(81.*pow4(Mst2)) + (35982464*Dmglst2*z3*pow4(Mst1))/(243.*
        Mgl*pow4(Mst2)) - (18688*z4*pow4(Mst1))/(81.*pow4(Mst2)) - (357248*
        Dmglst2*z4*pow4(Mst1))/(243.*Mgl*pow4(Mst2)) - (291022313156*pow2(
        Dmglst2)*pow4(Mst1))/(6.251175e6*pow2(Mgl)*pow4(Mst2)) + (1428992*
        OepS2*pow2(Dmglst2)*pow4(Mst1))/(729.*pow2(Mgl)*pow4(Mst2)) - (
        32019872*S2*pow2(Dmglst2)*pow4(Mst1))/(945.*pow2(Mgl)*pow4(Mst2)) - (
        714496*T1ep*pow2(Dmglst2)*pow4(Mst1))/(243.*pow2(Mgl)*pow4(Mst2)) - (
        1089695168*z2*pow2(Dmglst2)*pow4(Mst1))/(8505.*pow2(Mgl)*pow4(Mst2)) +
        (35982464*z3*pow2(Dmglst2)*pow4(Mst1))/(243.*pow2(Mgl)*pow4(Mst2)) - (
        357248*z4*pow2(Dmglst2)*pow4(Mst1))/(243.*pow2(Mgl)*pow4(Mst2)) - (
        9344*pow2(z2)*pow4(Mst1))/(27.*pow4(Mst2)) - (178624*Dmglst2*pow2(z2)*
        pow4(Mst1))/(81.*Mgl*pow4(Mst2)) - (178624*pow2(Dmglst2)*pow2(z2)*pow4(
        Mst1))/(81.*pow2(Mgl)*pow4(Mst2)) - (3440*pow4(Mst2))/(9.*pow2(Msq)*
        pow2(Mst1)) - (8960*Dmglst2*pow4(Mst2))/(9.*Mgl*pow2(Msq)*pow2(Mst1)) -
        (12800*pow2(Dmglst2)*pow4(Mst2))/(9.*pow2(Mgl)*pow2(Msq)*pow2(Mst1)) -
        (181136*pow2(s2t)*pow4(Mst2))/(675.*pow2(Msq)*pow2(Mt)) - (8560*
        Dmglst2*pow2(s2t)*pow4(Mst2))/(9.*Mgl*pow2(Msq)*pow2(Mt)) - (12520*
        pow2(Dmglst2)*pow2(s2t)*pow4(Mst2))/(9.*pow2(Mgl)*pow2(Msq)*pow2(Mt)) -
        (412*pow2(s2t)*pow4(Mst2))/(pow2(Mst1)*pow2(Mt)) + (1432*Dmglst2*pow2(
        s2t)*pow4(Mst2))/(3.*Mgl*pow2(Mst1)*pow2(Mt)) - (784*z2*pow2(s2t)*pow4(
        Mst2))/(3.*pow2(Mst1)*pow2(Mt)) - (1472*Dmglst2*z2*pow2(s2t)*pow4(Mst2)
        )/(3.*Mgl*pow2(Mst1)*pow2(Mt)) + (4340*pow2(Dmglst2)*pow2(s2t)*pow4(
        Mst2))/(3.*pow2(Mgl)*pow2(Mst1)*pow2(Mt)) - (768*z2*pow2(Dmglst2)*pow2(
        s2t)*pow4(Mst2))/(pow2(Mgl)*pow2(Mst1)*pow2(Mt)) - (142320737*pow4(
        Mst2))/(231525.*pow4(Msq)) - (25340*Dmglst2*pow4(Mst2))/(27.*Mgl*pow4(
        Msq)) + (1120*z2*pow4(Mst2))/(3.*pow4(Msq)) + (2560*Dmglst2*z2*pow4(
        Mst2))/(3.*Mgl*pow4(Msq)) - (180526*pow2(Dmglst2)*pow4(Mst2))/(135.*
        pow2(Mgl)*pow4(Msq)) + (5440*z2*pow2(Dmglst2)*pow4(Mst2))/(3.*pow2(Mgl)
        *pow4(Msq)) + (522392*pow2(Mst1)*pow2(s2t)*pow4(Mst2))/(15435.*pow2(Mt)
        *pow4(Msq)) + (650*Dmglst2*pow2(Mst1)*pow2(s2t)*pow4(Mst2))/(9.*Mgl*
        pow2(Mt)*pow4(Msq)) + (155*pow2(Dmglst2)*pow2(Mst1)*pow2(s2t)*pow4(
        Mst2))/(9.*pow2(Mgl)*pow2(Mt)*pow4(Msq)) - (256*pow4(Mst2))/(3.*pow4(
        Mst1)) + (1024*pow2(Dmglst2)*pow4(Mst2))/(3.*pow2(Mgl)*pow4(Mst1)) - (
        10*pow2(Msq)*pow2(Mst1)*pow4(s2t))/pow4(Mt) + (20*z2*pow2(Msq)*pow2(
        Mst1)*pow4(s2t))/pow4(Mt) - (10*pow2(Msq)*pow2(Mst2)*pow4(s2t))/pow4(
        Mt) + (20*z2*pow2(Msq)*pow2(Mst2)*pow4(s2t))/pow4(Mt) - (9250423*pow2(
        Mst1)*pow2(Mst2)*pow4(s2t))/(8100.*pow4(Mt)) + (4*B4*pow2(Mst1)*pow2(
        Mst2)*pow4(s2t))/(3.*pow4(Mt)) - (8*D3*pow2(Mst1)*pow2(Mst2)*pow4(s2t))
        /(3.*pow4(Mt)) + (2*DN*pow2(Mst1)*pow2(Mst2)*pow4(s2t))/pow4(Mt) + (
        1727459*Dmglst2*pow2(Mst1)*pow2(Mst2)*pow4(s2t))/(10125.*Mgl*pow4(Mt))
        + (8*B4*Dmglst2*pow2(Mst1)*pow2(Mst2)*pow4(s2t))/(Mgl*pow4(Mt)) - (32*
        D3*Dmglst2*pow2(Mst1)*pow2(Mst2)*pow4(s2t))/(3.*Mgl*pow4(Mt)) + (20*
        Dmglst2*DN*pow2(Mst1)*pow2(Mst2)*pow4(s2t))/(3.*Mgl*pow4(Mt)) + (2744*
        OepS2*pow2(Mst1)*pow2(Mst2)*pow4(s2t))/(81.*pow4(Mt)) - (128*Dmglst2*
        OepS2*pow2(Mst1)*pow2(Mst2)*pow4(s2t))/(3.*Mgl*pow4(Mt)) - (3137*S2*
        pow2(Mst1)*pow2(Mst2)*pow4(s2t))/pow4(Mt) - (2200*Dmglst2*S2*pow2(Mst1)
        *pow2(Mst2)*pow4(s2t))/(3.*Mgl*pow4(Mt)) - (1372*T1ep*pow2(Mst1)*pow2(
        Mst2)*pow4(s2t))/(27.*pow4(Mt)) + (64*Dmglst2*T1ep*pow2(Mst1)*pow2(
        Mst2)*pow4(s2t))/(Mgl*pow4(Mt)) - (235099*z2*pow2(Mst1)*pow2(Mst2)*
        pow4(s2t))/(270.*pow4(Mt)) + (15568*Dmglst2*z2*pow2(Mst1)*pow2(Mst2)*
        pow4(s2t))/(45.*Mgl*pow4(Mt)) + (22544*z3*pow2(Mst1)*pow2(Mst2)*pow4(
        s2t))/(27.*pow4(Mt)) + (4478*Dmglst2*z3*pow2(Mst1)*pow2(Mst2)*pow4(s2t)
        )/(9.*Mgl*pow4(Mt)) - (956*z4*pow2(Mst1)*pow2(Mst2)*pow4(s2t))/(27.*
        pow4(Mt)) - (12*Dmglst2*z4*pow2(Mst1)*pow2(Mst2)*pow4(s2t))/(Mgl*pow4(
        Mt)) - (65965567*pow2(Dmglst2)*pow2(Mst1)*pow2(Mst2)*pow4(s2t))/(60750.
        *pow2(Mgl)*pow4(Mt)) + (12*B4*pow2(Dmglst2)*pow2(Mst1)*pow2(Mst2)*pow4(
        s2t))/(pow2(Mgl)*pow4(Mt)) - (16*D3*pow2(Dmglst2)*pow2(Mst1)*pow2(Mst2)
        *pow4(s2t))/(pow2(Mgl)*pow4(Mt)) + (10*DN*pow2(Dmglst2)*pow2(Mst1)*
        pow2(Mst2)*pow4(s2t))/(pow2(Mgl)*pow4(Mt)) - (1280*OepS2*pow2(Dmglst2)*
        pow2(Mst1)*pow2(Mst2)*pow4(s2t))/(243.*pow2(Mgl)*pow4(Mt)) - (89092*S2*
        pow2(Dmglst2)*pow2(Mst1)*pow2(Mst2)*pow4(s2t))/(45.*pow2(Mgl)*pow4(Mt))
        + (640*T1ep*pow2(Dmglst2)*pow2(Mst1)*pow2(Mst2)*pow4(s2t))/(81.*pow2(
        Mgl)*pow4(Mt)) - (69148*z2*pow2(Dmglst2)*pow2(Mst1)*pow2(Mst2)*pow4(
        s2t))/(405.*pow2(Mgl)*pow4(Mt)) + (139555*z3*pow2(Dmglst2)*pow2(Mst1)*
        pow2(Mst2)*pow4(s2t))/(81.*pow2(Mgl)*pow4(Mt)) - (5026*z4*pow2(Dmglst2)
        *pow2(Mst1)*pow2(Mst2)*pow4(s2t))/(81.*pow2(Mgl)*pow4(Mt)) - (559*pow2(
        Mst1)*pow2(Mst2)*pow2(z2)*pow4(s2t))/(9.*pow4(Mt)) - (1784*pow2(
        Dmglst2)*pow2(Mst1)*pow2(Mst2)*pow2(z2)*pow4(s2t))/(27.*pow2(Mgl)*pow4(
        Mt)) + (71378675209*pow4(Mst1)*pow4(s2t))/(1.66698e8*pow4(Mt)) + (12*
        B4*pow4(Mst1)*pow4(s2t))/pow4(Mt) + (4*D3*pow4(Mst1)*pow4(s2t))/(3.*
        pow4(Mt)) - (4*DN*pow4(Mst1)*pow4(s2t))/(3.*pow4(Mt)) + (126628974763*
        Dmglst2*pow4(Mst1)*pow4(s2t))/(1.250235e8*Mgl*pow4(Mt)) + (7780*OepS2*
        pow4(Mst1)*pow4(s2t))/(243.*pow4(Mt)) - (67184*Dmglst2*OepS2*pow4(Mst1)
        *pow4(s2t))/(729.*Mgl*pow4(Mt)) - (18995*S2*pow4(Mst1)*pow4(s2t))/(9.*
        pow4(Mt)) + (103270*Dmglst2*S2*pow4(Mst1)*pow4(s2t))/(27.*Mgl*pow4(Mt))
        - (3890*T1ep*pow4(Mst1)*pow4(s2t))/(81.*pow4(Mt)) + (33592*Dmglst2*
        T1ep*pow4(Mst1)*pow4(s2t))/(243.*Mgl*pow4(Mt)) - (150437*z2*pow4(Mst1)*
        pow4(s2t))/(324.*pow4(Mt)) + (2200091*Dmglst2*z2*pow4(Mst1)*pow4(s2t))/
        (1215.*Mgl*pow4(Mt)) + (17035*z3*pow4(Mst1)*pow4(s2t))/(81.*pow4(Mt)) -
        (197843*Dmglst2*z3*pow4(Mst1)*pow4(s2t))/(243.*Mgl*pow4(Mt)) - (1297*
        z4*pow4(Mst1)*pow4(s2t))/(81.*pow4(Mt)) + (16796*Dmglst2*z4*pow4(Mst1)*
        pow4(s2t))/(243.*Mgl*pow4(Mt)) + (126628974763*pow2(Dmglst2)*pow4(Mst1)
        *pow4(s2t))/(1.250235e8*pow2(Mgl)*pow4(Mt)) - (67184*OepS2*pow2(
        Dmglst2)*pow4(Mst1)*pow4(s2t))/(729.*pow2(Mgl)*pow4(Mt)) + (103270*S2*
        pow2(Dmglst2)*pow4(Mst1)*pow4(s2t))/(27.*pow2(Mgl)*pow4(Mt)) + (33592*
        T1ep*pow2(Dmglst2)*pow4(Mst1)*pow4(s2t))/(243.*pow2(Mgl)*pow4(Mt)) + (
        2200091*z2*pow2(Dmglst2)*pow4(Mst1)*pow4(s2t))/(1215.*pow2(Mgl)*pow4(
        Mt)) - (197843*z3*pow2(Dmglst2)*pow4(Mst1)*pow4(s2t))/(243.*pow2(Mgl)*
        pow4(Mt)) + (16796*z4*pow2(Dmglst2)*pow4(Mst1)*pow4(s2t))/(243.*pow2(
        Mgl)*pow4(Mt)) + (10*pow2(Msq)*pow4(Mst1)*pow4(s2t))/(pow2(Mst2)*pow4(
        Mt)) - (20*z2*pow2(Msq)*pow4(Mst1)*pow4(s2t))/(pow2(Mst2)*pow4(Mt)) + (
        2953*pow2(Mst2)*pow4(Mst1)*pow4(s2t))/(90.*pow2(Msq)*pow4(Mt)) - (125*
        Dmglst2*pow2(Mst2)*pow4(Mst1)*pow4(s2t))/(9.*Mgl*pow2(Msq)*pow4(Mt)) -
        (10*z2*pow2(Mst2)*pow4(Mst1)*pow4(s2t))/(pow2(Msq)*pow4(Mt)) - (20*
        Dmglst2*z2*pow2(Mst2)*pow4(Mst1)*pow4(s2t))/(Mgl*pow2(Msq)*pow4(Mt)) -
        (125*pow2(Dmglst2)*pow2(Mst2)*pow4(Mst1)*pow4(s2t))/(9.*pow2(Mgl)*pow2(
        Msq)*pow4(Mt)) - (20*z2*pow2(Dmglst2)*pow2(Mst2)*pow4(Mst1)*pow4(s2t))/
        (pow2(Mgl)*pow2(Msq)*pow4(Mt)) - (1945*pow2(z2)*pow4(Mst1)*pow4(s2t))/(
        54.*pow4(Mt)) + (8398*Dmglst2*pow2(z2)*pow4(Mst1)*pow4(s2t))/(81.*Mgl*
        pow4(Mt)) + (8398*pow2(Dmglst2)*pow2(z2)*pow4(Mst1)*pow4(s2t))/(81.*
        pow2(Mgl)*pow4(Mt)) + (9403*pow4(Mst2)*pow4(s2t))/(54.*pow4(Mt)) - (40*
        B4*pow4(Mst2)*pow4(s2t))/(3.*pow4(Mt)) + (4*D3*pow4(Mst2)*pow4(s2t))/(
        3.*pow4(Mt)) - (2*DN*pow4(Mst2)*pow4(s2t))/(3.*pow4(Mt)) + (7337*
        Dmglst2*pow4(Mst2)*pow4(s2t))/(54.*Mgl*pow4(Mt)) - (8*B4*Dmglst2*pow4(
        Mst2)*pow4(s2t))/(Mgl*pow4(Mt)) + (32*D3*Dmglst2*pow4(Mst2)*pow4(s2t))/
        (3.*Mgl*pow4(Mt)) - (20*Dmglst2*DN*pow4(Mst2)*pow4(s2t))/(3.*Mgl*pow4(
        Mt)) + (64*OepS2*pow4(Mst2)*pow4(s2t))/(9.*pow4(Mt)) - (112*Dmglst2*
        OepS2*pow4(Mst2)*pow4(s2t))/(27.*Mgl*pow4(Mt)) - (3168*S2*pow4(Mst2)*
        pow4(s2t))/pow4(Mt) - (10062*Dmglst2*S2*pow4(Mst2)*pow4(s2t))/(Mgl*
        pow4(Mt)) - (32*T1ep*pow4(Mst2)*pow4(s2t))/(3.*pow4(Mt)) + (56*Dmglst2*
        T1ep*pow4(Mst2)*pow4(s2t))/(9.*Mgl*pow4(Mt)) - (529*z2*pow4(Mst2)*pow4(
        s2t))/(9.*pow4(Mt)) + (3047*Dmglst2*z2*pow4(Mst2)*pow4(s2t))/(9.*Mgl*
        pow4(Mt)) + (1913*z3*pow4(Mst2)*pow4(s2t))/(3.*pow4(Mt)) + (16175*
        Dmglst2*z3*pow4(Mst2)*pow4(s2t))/(9.*Mgl*pow4(Mt)) - (10*z4*pow4(Mst2)*
        pow4(s2t))/(3.*pow4(Mt)) + (424*Dmglst2*z4*pow4(Mst2)*pow4(s2t))/(9.*
        Mgl*pow4(Mt)) + (1200149*pow2(Dmglst2)*pow4(Mst2)*pow4(s2t))/(1620.*
        pow2(Mgl)*pow4(Mt)) - (100*B4*pow2(Dmglst2)*pow4(Mst2)*pow4(s2t))/(3.*
        pow2(Mgl)*pow4(Mt)) + (112*D3*pow2(Dmglst2)*pow4(Mst2)*pow4(s2t))/(3.*
        pow2(Mgl)*pow4(Mt)) - (62*DN*pow2(Dmglst2)*pow4(Mst2)*pow4(s2t))/(3.*
        pow2(Mgl)*pow4(Mt)) - (280*OepS2*pow2(Dmglst2)*pow4(Mst2)*pow4(s2t))/(
        81.*pow2(Mgl)*pow4(Mt)) - (116129*S2*pow2(Dmglst2)*pow4(Mst2)*pow4(s2t)
        )/(5.*pow2(Mgl)*pow4(Mt)) + (140*T1ep*pow2(Dmglst2)*pow4(Mst2)*pow4(
        s2t))/(27.*pow2(Mgl)*pow4(Mt)) + (64091*z2*pow2(Dmglst2)*pow4(Mst2)*
        pow4(s2t))/(270.*pow2(Mgl)*pow4(Mt)) + (226447*z3*pow2(Dmglst2)*pow4(
        Mst2)*pow4(s2t))/(54.*pow2(Mgl)*pow4(Mt)) + (4444*z4*pow2(Dmglst2)*
        pow4(Mst2)*pow4(s2t))/(27.*pow2(Mgl)*pow4(Mt)) + (10*pow2(Msq)*pow4(
        Mst2)*pow4(s2t))/(pow2(Mst1)*pow4(Mt)) - (20*z2*pow2(Msq)*pow4(Mst2)*
        pow4(s2t))/(pow2(Mst1)*pow4(Mt)) - (2636*pow2(Mst1)*pow4(Mst2)*pow4(
        s2t))/(45.*pow2(Msq)*pow4(Mt)) - (1540*Dmglst2*pow2(Mst1)*pow4(Mst2)*
        pow4(s2t))/(9.*Mgl*pow2(Msq)*pow4(Mt)) + (100*z2*pow2(Mst1)*pow4(Mst2)*
        pow4(s2t))/(3.*pow2(Msq)*pow4(Mt)) + (400*Dmglst2*z2*pow2(Mst1)*pow4(
        Mst2)*pow4(s2t))/(3.*Mgl*pow2(Msq)*pow4(Mt)) - (3590*pow2(Dmglst2)*
        pow2(Mst1)*pow4(Mst2)*pow4(s2t))/(9.*pow2(Mgl)*pow2(Msq)*pow4(Mt)) + (
        1000*z2*pow2(Dmglst2)*pow2(Mst1)*pow4(Mst2)*pow4(s2t))/(3.*pow2(Mgl)*
        pow2(Msq)*pow4(Mt)) + (16*pow2(z2)*pow4(Mst2)*pow4(s2t))/pow4(Mt) + (
        158*Dmglst2*pow2(z2)*pow4(Mst2)*pow4(s2t))/(3.*Mgl*pow4(Mt)) + (683*
        pow2(Dmglst2)*pow2(z2)*pow4(Mst2)*pow4(s2t))/(9.*pow2(Mgl)*pow4(Mt)) +
        (2024719*pow4(Mst1)*pow4(Mst2)*pow4(s2t))/(205800.*pow4(Msq)*pow4(Mt))
        + (565*Dmglst2*pow4(Mst1)*pow4(Mst2)*pow4(s2t))/(36.*Mgl*pow4(Msq)*
        pow4(Mt)) - (25*z2*pow4(Mst1)*pow4(Mst2)*pow4(s2t))/(6.*pow4(Msq)*pow4(
        Mt)) - (50*Dmglst2*z2*pow4(Mst1)*pow4(Mst2)*pow4(s2t))/(3.*Mgl*pow4(
        Msq)*pow4(Mt)) + (565*pow2(Dmglst2)*pow4(Mst1)*pow4(Mst2)*pow4(s2t))/(
        36.*pow2(Mgl)*pow4(Msq)*pow4(Mt)) - (50*z2*pow2(Dmglst2)*pow4(Mst1)*
        pow4(Mst2)*pow4(s2t))/(3.*pow2(Mgl)*pow4(Msq)*pow4(Mt)) + (2*pow3(log(
        pow2(Mst1)/pow2(Mst2)))*(2*Dmglst2*Mgl*(pow2(Mst1)*pow2(Mst2)*(2312*
        pow2(Mst2)*pow2(Mt)*pow2(s2t) + 768*Mst2*s2t*pow3(Mt) + 960*Mt*pow3(
        Mst2)*pow3(s2t) - 3824*pow4(Mt) - 77*pow4(Mst2)*pow4(s2t)) + 24*pow4(
        Mst2)*(32*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 128*Mst2*s2t*pow3(Mt) + 72*
        Mt*pow3(Mst2)*pow3(s2t) - 64*pow4(Mt) + 19*pow4(Mst2)*pow4(s2t)) +
        pow4(Mst1)*(2896*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 2304*Mst2*s2t*pow3(Mt)
        + 672*Mt*pow3(Mst2)*pow3(s2t) + 5536*pow4(Mt) + 33*pow4(Mst2)*pow4(s2t)
        )) + pow2(Dmglst2)*(pow2(Mst1)*pow2(Mst2)*(6904*pow2(Mst2)*pow2(Mt)*
        pow2(s2t) - 1536*Mst2*s2t*pow3(Mt) + 1920*Mt*pow3(Mst2)*pow3(s2t) -
        2320*pow4(Mt) - 211*pow4(Mst2)*pow4(s2t)) + 2*pow4(Mst1)*(2896*pow2(
        Mst2)*pow2(Mt)*pow2(s2t) - 2304*Mst2*s2t*pow3(Mt) + 672*Mt*pow3(Mst2)*
        pow3(s2t) + 5536*pow4(Mt) + 33*pow4(Mst2)*pow4(s2t)) + 24*pow4(Mst2)*(
        96*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 256*Mst2*s2t*pow3(Mt) + 240*Mt*pow3(
        Mst2)*pow3(s2t) - 192*pow4(Mt) + 57*pow4(Mst2)*pow4(s2t))) + pow2(Mgl)*
        (-4*pow4(Mst1)*(-337*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 256*Mst2*s2t*pow3(
        Mt) - 400*Mt*pow3(Mst2)*pow3(s2t) + 2420*pow4(Mt) - 68*pow4(Mst2)*pow4(
        s2t)) + 12*pow4(Mst2)*(-8*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 512*Mst2*s2t*
        pow3(Mt) + 160*Mt*pow3(Mst2)*pow3(s2t) + 74*pow4(Mt) + 47*pow4(Mst2)*
        pow4(s2t)) + pow2(Mst1)*pow2(Mst2)*(2152*pow2(Mst2)*pow2(Mt)*pow2(s2t)
        - 4608*Mst2*s2t*pow3(Mt) + 1920*Mt*pow3(Mst2)*pow3(s2t) - 2320*pow4(Mt)
        + 185*pow4(Mst2)*pow4(s2t)))))/(9.*pow2(Mgl)*pow4(Mst2)*pow4(Mt)) - (
        5680*pow3(s2t)*pow5(Mst2))/(9.*pow2(Msq)*pow3(Mt)) - (7600*Dmglst2*
        pow3(s2t)*pow5(Mst2))/(3.*Mgl*pow2(Msq)*pow3(Mt)) + (640*z2*pow3(s2t)*
        pow5(Mst2))/(3.*pow2(Msq)*pow3(Mt)) + (1280*Dmglst2*z2*pow3(s2t)*pow5(
        Mst2))/(Mgl*pow2(Msq)*pow3(Mt)) - (63520*pow2(Dmglst2)*pow3(s2t)*pow5(
        Mst2))/(9.*pow2(Mgl)*pow2(Msq)*pow3(Mt)) + (4160*z2*pow2(Dmglst2)*pow3(
        s2t)*pow5(Mst2))/(pow2(Mgl)*pow2(Msq)*pow3(Mt)) + (256*pow3(s2t)*pow5(
        Mst2))/(3.*pow2(Mst1)*pow3(Mt)) - (256*Dmglst2*pow3(s2t)*pow5(Mst2))/(
        3.*Mgl*pow2(Mst1)*pow3(Mt)) - (1024*pow2(Dmglst2)*pow3(s2t)*pow5(Mst2))
        /(3.*pow2(Mgl)*pow2(Mst1)*pow3(Mt)) - (8240*s2t*pow5(Mst2))/(9.*Mt*
        pow4(Msq)) - (30160*Dmglst2*s2t*pow5(Mst2))/(9.*Mgl*Mt*pow4(Msq)) + (
        640*s2t*z2*pow5(Mst2))/(Mt*pow4(Msq)) + (7040*Dmglst2*s2t*z2*pow5(Mst2)
        )/(3.*Mgl*Mt*pow4(Msq)) - (28000*s2t*pow2(Dmglst2)*pow5(Mst2))/(3.*Mt*
        pow2(Mgl)*pow4(Msq)) + (18560*s2t*z2*pow2(Dmglst2)*pow5(Mst2))/(3.*Mt*
        pow2(Mgl)*pow4(Msq)) + (1115*pow2(Mst1)*pow3(s2t)*pow5(Mst2))/(9.*pow3(
        Mt)*pow4(Msq)) + (4115*Dmglst2*pow2(Mst1)*pow3(s2t)*pow5(Mst2))/(9.*
        Mgl*pow3(Mt)*pow4(Msq)) - (40*z2*pow2(Mst1)*pow3(s2t)*pow5(Mst2))/(
        pow3(Mt)*pow4(Msq)) - (200*Dmglst2*z2*pow2(Mst1)*pow3(s2t)*pow5(Mst2))/
        (Mgl*pow3(Mt)*pow4(Msq)) + (2965*pow2(Dmglst2)*pow2(Mst1)*pow3(s2t)*
        pow5(Mst2))/(3.*pow2(Mgl)*pow3(Mt)*pow4(Msq)) - (600*z2*pow2(Dmglst2)*
        pow2(Mst1)*pow3(s2t)*pow5(Mst2))/(pow2(Mgl)*pow3(Mt)*pow4(Msq)) + (
        pow2(log(pow2(Mst1)/pow2(Msq)))*(-105*log(pow2(Mst1)/pow2(Mst2))*pow2(
        Mst1)*pow4(Msq)*(2*Mst2*s2t*pow2(Dmglst2)*(-8*Mt*(pow2(Mst1) + pow2(
        Mst2))*pow2(s2t) + 32*pow3(Mt) + 9*Mst2*(pow2(Mst1) - pow2(Mst2))*pow3(
        s2t)) - 4*Dmglst2*Mgl*Mst2*s2t*(4*Mt*(pow2(Mst1) + pow2(Mst2))*pow2(
        s2t) - 16*pow3(Mt) + 3*Mst2*(-pow2(Mst1) + pow2(Mst2))*pow3(s2t)) +
        pow2(Mgl)*(8*(-pow2(Mst1) + pow2(Mst2))*pow2(Mt)*pow2(s2t) + 64*Mst2*
        s2t*pow3(Mt) - 16*Mst2*Mt*(pow2(Mst1) + pow2(Mst2))*pow3(s2t) + 16*
        pow4(Mt) + (6*pow2(Mst1)*pow2(Mst2) + pow4(Mst1) - 7*pow4(Mst2))*pow4(
        s2t))) + 70*Dmglst2*Mgl*(-64*Mst2*pow2(Msq)*pow2(Mst1)*(-3*Mst2*Mt +
        s2t*pow2(Mst1) - s2t*pow2(Mst2))*pow3(Mt) - 16*Mst2*pow2(Mst1)*pow3(Mt)
        *(3*Mst2*(-Mt + 2*Mst2*s2t)*pow2(Mst1) - (13*Mt + 7*Mst2*s2t)*pow3(
        Mst2) + s2t*pow4(Mst1)) + 3*pow4(Msq)*(3*pow2(-4*Mst2*pow2(Mt) + pow2(
        s2t)*pow3(Mst2)) - pow2(s2t)*(16*Mst2*Mt*s2t + 24*pow2(Mt) + 3*pow2(
        Mst2)*pow2(s2t))*pow4(Mst1) + pow2(Mst1)*(48*pow2(Mst2)*pow2(Mt)*pow2(
        s2t) + 16*Mt*pow3(Mst2)*pow3(s2t) + 48*pow4(Mt) - 3*pow4(Mst2)*pow4(
        s2t)) + 3*pow4(s2t)*pow6(Mst1))) + 35*pow2(Dmglst2)*(-64*Mst2*pow2(Msq)
        *pow2(Mst1)*(-9*Mst2*Mt + 2*s2t*pow2(Mst1) - 2*s2t*pow2(Mst2))*pow3(Mt)
        - 16*Mst2*pow2(Mst1)*pow3(Mt)*(3*Mst2*(-3*Mt + 8*Mst2*s2t)*pow2(Mst1) -
        (59*Mt + 26*Mst2*s2t)*pow3(Mst2) + 2*s2t*pow4(Mst1)) + 3*pow4(Msq)*(9*
        pow2(-4*Mst2*pow2(Mt) + pow2(s2t)*pow3(Mst2)) - pow2(s2t)*(32*Mst2*Mt*
        s2t + 72*pow2(Mt) + 9*pow2(Mst2)*pow2(s2t))*pow4(Mst1) + pow2(Mst1)*(
        144*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 32*Mt*pow3(Mst2)*pow3(s2t) + 144*
        pow4(Mt) - 9*pow4(Mst2)*pow4(s2t)) + 6*pow4(s2t)*pow6(Mst1))) + pow2(
        Mgl)*(224*pow2(Msq)*pow2(Mst1)*pow2(Mt)*(pow2(Mst1)*(-20*Mst2*Mt*s2t +
        7*pow2(Mt) - 2*pow2(Mst2)*pow2(s2t)) + pow2(Mst2)*(20*Mst2*Mt*s2t + 37*
        pow2(Mt) + pow2(Mst2)*pow2(s2t)) + pow2(s2t)*pow4(Mst1)) + 4*pow2(Mst1)
        *pow2(Mt)*(840*Mst2*Mt*s2t + 1234*pow2(Mt) + 15*pow2(Mst2)*pow2(s2t))*
        pow4(Mst2) + 4*pow2(Mst2)*pow4(Mst1)*(-15*pow2(Mst2)*pow2(Mt)*pow2(s2t)
        - 560*Mst2*s2t*pow3(Mt) + 420*pow4(Mt) + pow4(Mst2)*pow4(s2t)) + 2*(-
        30*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 560*Mst2*s2t*pow3(Mt) + 228*pow4(Mt)
        - 3*pow4(Mst2)*pow4(s2t))*pow6(Mst1) + 35*pow4(Msq)*(9*pow2(-4*Mst2*
        pow2(Mt) + pow2(s2t)*pow3(Mst2)) - 3*pow2(s2t)*(32*Mst2*Mt*s2t + 24*
        pow2(Mt) + 7*pow2(Mst2)*pow2(s2t))*pow4(Mst1) + pow2(Mst1)*(144*pow2(
        Mst2)*pow2(Mt)*pow2(s2t) + 96*Mt*pow3(Mst2)*pow3(s2t) + 176*pow4(Mt) -
        3*pow4(Mst2)*pow4(s2t)) + 15*pow4(s2t)*pow6(Mst1)))))/(21.*pow2(Mgl)*
        pow2(Mst1)*pow4(Msq)*pow4(Mt)) + (1720*pow2(s2t)*pow6(Mst2))/(9.*pow2(
        Msq)*pow2(Mst1)*pow2(Mt)) + (4480*Dmglst2*pow2(s2t)*pow6(Mst2))/(9.*
        Mgl*pow2(Msq)*pow2(Mst1)*pow2(Mt)) + (6400*pow2(Dmglst2)*pow2(s2t)*
        pow6(Mst2))/(9.*pow2(Mgl)*pow2(Msq)*pow2(Mst1)*pow2(Mt)) - (670*pow6(
        Mst2))/(9.*pow2(Mst1)*pow4(Msq)) - (260*Dmglst2*pow6(Mst2))/(Mgl*pow2(
        Mst1)*pow4(Msq)) - (350*pow2(Dmglst2)*pow6(Mst2))/(pow2(Mgl)*pow2(Mst1)
        *pow4(Msq)) - (107588*pow2(s2t)*pow6(Mst2))/(1715.*pow2(Mt)*pow4(Msq))
        - (2030*Dmglst2*pow2(s2t)*pow6(Mst2))/(9.*Mgl*pow2(Mt)*pow4(Msq)) - (
        2405*pow2(Dmglst2)*pow2(s2t)*pow6(Mst2))/(9.*pow2(Mgl)*pow2(Mt)*pow4(
        Msq)) + (128*pow2(s2t)*pow6(Mst2))/(3.*pow2(Mt)*pow4(Mst1)) - (512*
        pow2(Dmglst2)*pow2(s2t)*pow6(Mst2))/(3.*pow2(Mgl)*pow2(Mt)*pow4(Mst1))
        + (2986*pow4(s2t)*pow6(Mst2))/(45.*pow2(Msq)*pow4(Mt)) + (2300*Dmglst2*
        pow4(s2t)*pow6(Mst2))/(9.*Mgl*pow2(Msq)*pow4(Mt)) - (70*z2*pow4(s2t)*
        pow6(Mst2))/(3.*pow2(Msq)*pow4(Mt)) - (340*Dmglst2*z2*pow4(s2t)*pow6(
        Mst2))/(3.*Mgl*pow2(Msq)*pow4(Mt)) + (5020*pow2(Dmglst2)*pow4(s2t)*
        pow6(Mst2))/(9.*pow2(Mgl)*pow2(Msq)*pow4(Mt)) - (910*z2*pow2(Dmglst2)*
        pow4(s2t)*pow6(Mst2))/(3.*pow2(Mgl)*pow2(Msq)*pow4(Mt)) + (135*pow4(
        s2t)*pow6(Mst2))/(2.*pow2(Mst1)*pow4(Mt)) - (115*Dmglst2*pow4(s2t)*
        pow6(Mst2))/(3.*Mgl*pow2(Mst1)*pow4(Mt)) + (98*z2*pow4(s2t)*pow6(Mst2))
        /(3.*pow2(Mst1)*pow4(Mt)) + (184*Dmglst2*z2*pow4(s2t)*pow6(Mst2))/(3.*
        Mgl*pow2(Mst1)*pow4(Mt)) - (1405*pow2(Dmglst2)*pow4(s2t)*pow6(Mst2))/(
        6.*pow2(Mgl)*pow2(Mst1)*pow4(Mt)) + (96*z2*pow2(Dmglst2)*pow4(s2t)*
        pow6(Mst2))/(pow2(Mgl)*pow2(Mst1)*pow4(Mt)) - (45671209*pow2(Mst1)*
        pow4(s2t)*pow6(Mst2))/(1.8522e6*pow4(Msq)*pow4(Mt)) - (3965*Dmglst2*
        pow2(Mst1)*pow4(s2t)*pow6(Mst2))/(36.*Mgl*pow4(Msq)*pow4(Mt)) + (35*z2*
        pow2(Mst1)*pow4(s2t)*pow6(Mst2))/(3.*pow4(Msq)*pow4(Mt)) + (70*Dmglst2*
        z2*pow2(Mst1)*pow4(s2t)*pow6(Mst2))/(Mgl*pow4(Msq)*pow4(Mt)) - (2755*
        pow2(Dmglst2)*pow2(Mst1)*pow4(s2t)*pow6(Mst2))/(8.*pow2(Mgl)*pow4(Msq)*
        pow4(Mt)) + (245*z2*pow2(Dmglst2)*pow2(Mst1)*pow4(s2t)*pow6(Mst2))/(
        pow2(Mgl)*pow4(Msq)*pow4(Mt)) - (27955*Dmglst2*pow3(s2t)*pow7(Mst2))/(
        27.*Mgl*pow3(Mt)*pow4(Msq)) + (1640*Dmglst2*z2*pow3(s2t)*pow7(Mst2))/(
        3.*Mgl*pow3(Mt)*pow4(Msq)) - (106375*pow2(Dmglst2)*pow3(s2t)*pow7(Mst2)
        )/(27.*pow2(Mgl)*pow3(Mt)*pow4(Msq)) + (7160*z2*pow2(Dmglst2)*pow3(s2t)
        *pow7(Mst2))/(3.*pow2(Mgl)*pow3(Mt)*pow4(Msq)) + (45*pow2(Dmglst2)*
        pow2(s2t)*pow8(Mst2))/(pow2(Mgl)*pow2(Mst1)*pow2(Mt)*pow4(Msq)) - (80*
        pow2(Dmglst2)*pow4(s2t)*pow8(Mst2))/(3.*pow2(Mgl)*pow2(Msq)*pow2(Mst1)*
        pow4(Mt)) + (16435*pow2(Dmglst2)*pow4(s2t)*pow8(Mst2))/(72.*pow2(Mgl)*
        pow4(Msq)*pow4(Mt)) - (150*z2*pow2(Dmglst2)*pow4(s2t)*pow8(Mst2))/(
        pow2(Mgl)*pow4(Msq)*pow4(Mt)) + (64*pow2(Dmglst2)*pow4(s2t)*pow8(Mst2))
        /(3.*pow2(Mgl)*pow4(Mst1)*pow4(Mt)) + (pow2(log(pow2(Mst1)/pow2(Mst2)))
        *(4*Dmglst2*Mgl*(-525*pow4(Mst1)*(20*Mt*pow3(Mst2)*(-12*Mst2*s2t*pow2(
        Mt) + 4*pow3(Mt) + 3*pow3(Mst2)*pow3(s2t)) + pow4(Mst1)*(-16*s2t*pow3(
        Mt) + 5*pow3(Mst2)*pow4(s2t)) + 3*Mst2*pow2(Mst1)*(-32*Mst2*s2t*pow3(
        Mt) + 20*Mt*pow3(Mst2)*pow3(s2t) + 16*pow4(Mt) - 7*pow4(Mst2)*pow4(s2t)
        ))*pow5(Mst2) - 1050*pow2(Msq)*pow2(Mst2)*pow4(Mst1)*(3*pow4(Mst1)*(16*
        pow4(Mt) + pow4(Mst2)*pow4(s2t)) + pow4(Mst2)*(-192*Mst2*s2t*pow3(Mt) +
        48*Mt*pow3(Mst2)*pow3(s2t) + 48*pow4(Mt) + 17*pow4(Mst2)*pow4(s2t)) -
        4*pow2(Mst1)*(8*s2t*pow3(Mst2)*pow3(Mt) - 12*Mt*pow3(s2t)*pow5(Mst2) +
        5*pow4(s2t)*pow6(Mst2))) + pow4(Msq)*(35*pow4(Mst1)*pow4(Mst2)*(9384*
        pow2(Mst2)*pow2(Mt)*pow2(s2t) + 11168*Mst2*s2t*pow3(Mt) + 4056*Mt*pow3(
        Mst2)*pow3(s2t) - 16912*pow4(Mt) + 717*pow4(Mst2)*pow4(s2t)) - 28*pow2(
        Mst2)*(-998*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 36680*Mst2*s2t*pow3(Mt) +
        3150*Mt*pow3(Mst2)*pow3(s2t) + 37236*pow4(Mt) + 1411*pow4(Mst2)*pow4(
        s2t))*pow6(Mst1) + 840*pow2(Mst1)*(-136*pow2(Mst2)*pow2(Mt)*pow2(s2t) -
        192*Mst2*s2t*pow3(Mt) + 48*Mt*pow3(Mst2)*pow3(s2t) + 208*pow4(Mt) + 29*
        pow4(Mst2)*pow4(s2t))*pow6(Mst2) + (294096*pow2(Mst2)*pow2(Mt)*pow2(
        s2t) - 394240*Mst2*s2t*pow3(Mt) - 31080*Mt*pow3(Mst2)*pow3(s2t) +
        2542880*pow4(Mt) - 6102*pow4(Mst2)*pow4(s2t))*pow8(Mst1) + 26880*pow2(
        Mt)*(-2*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow8(Mst2))) + pow2(Mgl)*(-
        840*pow2(Msq)*pow2(Mst2)*pow4(Mst1)*(pow2(Mst1)*pow2(Mst2)*(-12*pow2(
        Mst2)*pow2(Mt)*pow2(s2t) - 160*Mst2*s2t*pow3(Mt) + 80*Mt*pow3(Mst2)*
        pow3(s2t) + 80*pow4(Mt) - 19*pow4(Mst2)*pow4(s2t)) + pow4(Mst1)*(-120*
        pow4(Mt) + 9*pow4(Mst2)*pow4(s2t)) + 2*pow4(Mst2)*(-9*pow2(Mst2)*pow2(
        Mt)*pow2(s2t) - 160*Mst2*s2t*pow3(Mt) + 40*Mt*pow3(Mst2)*pow3(s2t) +
        72*pow4(Mt) + 10*pow4(Mst2)*pow4(s2t))) + 30*(192*pow2(Mst2)*pow2(Mt)*
        pow2(s2t) + 2240*Mst2*s2t*pow3(Mt) - 840*Mt*pow3(Mst2)*pow3(s2t) -
        1680*pow4(Mt) + 207*pow4(Mst2)*pow4(s2t))*pow6(Mst1)*pow6(Mst2) - 15*
        pow4(Mst2)*(-2240*Mst2*s2t*pow3(Mt) + 1680*pow4(Mt) + 223*pow4(Mst2)*
        pow4(s2t))*pow8(Mst1) + 240*pow2(Mt)*(420*Mst2*Mt*s2t - 223*pow2(Mt) +
        11*pow2(Mst2)*pow2(s2t))*pow4(Mst1)*pow8(Mst2) + pow4(Msq)*(280*pow4(
        Mst1)*pow4(Mst2)*(1536*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 8880*Mst2*s2t*
        pow3(Mt) + 3372*Mt*pow3(Mst2)*pow3(s2t) + 3212*pow4(Mt) + 879*pow4(
        Mst2)*pow4(s2t)) - 70*pow2(Mst2)*(-9560*pow2(Mst2)*pow2(Mt)*pow2(s2t) +
        27648*Mst2*s2t*pow3(Mt) - 5136*Mt*pow3(Mst2)*pow3(s2t) + 3136*pow4(Mt)
        + 289*pow4(Mst2)*pow4(s2t))*pow6(Mst1) + 1680*pow2(Mst1)*(-120*pow2(
        Mst2)*pow2(Mt)*pow2(s2t) - 128*Mst2*s2t*pow3(Mt) + 32*Mt*pow3(Mst2)*
        pow3(s2t) + 208*pow4(Mt) + 21*pow4(Mst2)*pow4(s2t))*pow6(Mst2) + (
        669176*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 313600*Mst2*s2t*pow3(Mt) +
        229600*Mt*pow3(Mst2)*pow3(s2t) - 4378256*pow4(Mt) + 56822*pow4(Mst2)*
        pow4(s2t))*pow8(Mst1) + 26880*pow2(Mt)*(-2*pow2(Mt) + pow2(Mst2)*pow2(
        s2t))*pow8(Mst2))) + 2*pow2(Dmglst2)*(-525*pow4(Mst1)*(pow4(Mst1)*(-32*
        s2t*pow3(Mt) + 10*pow3(Mst2)*pow4(s2t)) + 10*pow3(Mst2)*(-144*Mst2*s2t*
        pow3(Mt) + 36*Mt*pow3(Mst2)*pow3(s2t) + 40*pow4(Mt) + 9*pow4(Mst2)*
        pow4(s2t)) - 3*Mst2*pow2(Mst1)*(128*Mst2*s2t*pow3(Mt) - 120*Mt*pow3(
        Mst2)*pow3(s2t) - 48*pow4(Mt) + 49*pow4(Mst2)*pow4(s2t)))*pow5(Mst2) -
        1050*pow2(Msq)*pow2(Mst2)*pow4(Mst1)*(6*pow4(Mst1)*(16*pow4(Mt) + pow4(
        Mst2)*pow4(s2t)) + pow4(Mst2)*(-768*Mst2*s2t*pow3(Mt) + 192*Mt*pow3(
        Mst2)*pow3(s2t) + 144*pow4(Mt) + 91*pow4(Mst2)*pow4(s2t)) - 4*pow2(
        Mst1)*(16*s2t*pow3(Mst2)*pow3(Mt) - 48*Mt*pow3(s2t)*pow5(Mst2) + 25*
        pow4(s2t)*pow6(Mst2))) + pow4(Msq)*(-7*pow4(Mst1)*pow4(Mst2)*(-89160*
        pow2(Mst2)*pow2(Mt)*pow2(s2t) - 264000*Mst2*s2t*pow3(Mt) + 3600*Mt*
        pow3(Mst2)*pow3(s2t) + 10064*pow4(Mt) + 8715*pow4(Mst2)*pow4(s2t)) -
        56*pow2(Mst2)*(15049*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 5900*Mst2*s2t*
        pow3(Mt) + 7170*Mt*pow3(Mst2)*pow3(s2t) - 62654*pow4(Mt) + 1132*pow4(
        Mst2)*pow4(s2t))*pow6(Mst1) + 840*pow2(Mst1)*(-472*pow2(Mst2)*pow2(Mt)*
        pow2(s2t) - 768*Mst2*s2t*pow3(Mt) + 192*Mt*pow3(Mst2)*pow3(s2t) + 624*
        pow4(Mt) + 119*pow4(Mst2)*pow4(s2t))*pow6(Mst2) + 4*(147048*pow2(Mst2)*
        pow2(Mt)*pow2(s2t) - 197120*Mst2*s2t*pow3(Mt) - 15540*Mt*pow3(Mst2)*
        pow3(s2t) + 1271440*pow4(Mt) - 3051*pow4(Mst2)*pow4(s2t))*pow8(Mst1) -
        3360*(-40*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 80*pow4(Mt) + 3*pow4(Mst2)*
        pow4(s2t))*pow8(Mst2)))))/(630.*pow2(Mgl)*pow4(Msq)*pow4(Mst1)*pow4(
        Mst2)*pow4(Mt)) - (log(pow2(Mst1)/pow2(Msq))*(2450*Dmglst2*Mgl*(-(pow5(
        Mst2)*(4*Mt*pow2(Mst1)*pow3(Mst2)*(80*Mst2*s2t*pow2(Mt) + 270*Mt*pow2(
        Mst2)*pow2(s2t) + 548*pow3(Mt) + 181*pow3(Mst2)*pow3(s2t)) + 3*Mst2*
        pow4(Mst1)*(-216*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 64*Mst2*s2t*pow3(Mt) -
        244*Mt*pow3(Mst2)*pow3(s2t) + 128*pow4(Mt) + 19*pow4(Mst2)*pow4(s2t)) +
        1008*pow4(Mt)*pow5(Mst2) + s2t*(72*Mst2*s2t*pow2(Mt) + 12*Mt*pow2(Mst2)
        *pow2(s2t) - 128*pow3(Mt) + 3*pow3(Mst2)*pow3(s2t))*pow6(Mst1))) - 8*
        pow2(Msq)*pow5(Mst2)*(-120*pow2(Mt)*(-2*pow2(Mt) + pow2(Mst2)*pow2(s2t)
        )*pow3(Mst2) + s2t*(-120*Mst2*s2t*pow2(Mt) - 144*Mt*pow2(Mst2)*pow2(
        s2t) + 56*pow3(Mt) + 3*pow3(Mst2)*pow3(s2t))*pow4(Mst1) - 8*Mst2*pow2(
        Mst1)*(-30*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 7*Mst2*s2t*pow3(Mt) - 18*Mt*
        pow3(Mst2)*pow3(s2t) - 18*pow4(Mt) + 3*pow4(Mst2)*pow4(s2t)) + 6*Mst2*
        pow4(s2t)*pow6(Mst1)) + 2*pow4(Msq)*(pow2(Mst1)*pow4(Mst2)*(-120*pow2(
        Mst2)*pow2(Mt)*pow2(s2t) - 16*Mst2*s2t*(-161 + 72*z2)*pow3(Mt) - 144*
        Mt*(-15 + 11*z2)*pow3(Mst2)*pow3(s2t) + 368*pow4(Mt) + 27*(-5 + 4*z2)*
        pow4(Mst2)*pow4(s2t)) - 3*pow2(Mst2)*pow4(Mst1)*(128*pow2(Mst2)*pow2(
        Mt)*pow2(s2t) + 1248*Mst2*s2t*pow3(Mt) - 48*Mt*(2 + z2)*pow3(Mst2)*
        pow3(s2t) - 96*(27 + 4*z2)*pow4(Mt) + 9*(-3 + 4*z2)*pow4(Mst2)*pow4(
        s2t)) + 36*(-16*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 32*Mst2*s2t*(7 + 3*z2)*
        pow3(Mt) + 7*Mt*pow3(Mst2)*pow3(s2t) + 24*(23 + 12*z2)*pow4(Mt) + pow4(
        Mst2)*pow4(s2t))*pow6(Mst1) - 9*pow2(-4*pow2(Mt) + pow2(Mst2)*pow2(s2t)
        )*pow6(Mst2))) + pow2(Mgl)*(88200*pow2(Mst2)*pow6(Msq)*(pow2(-4*Mst2*
        pow2(Mt) + pow2(s2t)*pow3(Mst2)) - pow4(Mst1)*(8*pow2(Mt)*pow2(s2t) +
        pow2(Mst2)*pow4(s2t)) + pow2(Mst1)*(16*pow2(Mst2)*pow2(Mt)*pow2(s2t) +
        16*pow4(Mt) - pow4(Mst2)*pow4(s2t)) + pow4(s2t)*pow6(Mst1)) - 1225*
        pow4(Msq)*(-2*pow2(Mst1)*pow4(Mst2)*(720*pow2(Mst2)*pow2(Mt)*pow2(s2t)
        - 288*Mst2*s2t*(-13 + 8*z2)*pow3(Mt) - 864*Mt*(-2 + z2)*pow3(Mst2)*
        pow3(s2t) - 32*(13 + 18*z2)*pow4(Mt) + 3*(7 + 30*z2)*pow4(Mst2)*pow4(
        s2t)) + 6*pow2(Mst2)*pow4(Mst1)*(24*pow2(Mst2)*pow2(Mt)*pow2(s2t) -
        960*Mst2*s2t*pow3(Mt) - 96*Mt*(-3 + z2)*pow3(Mst2)*pow3(s2t) + 384*(3 +
        z2)*pow4(Mt) + (89 + 36*z2)*pow4(Mst2)*pow4(s2t)) + 3*(-768*Mst2*s2t*(1
        + 2*z2)*pow3(Mt) + 48*Mt*pow3(Mst2)*pow3(s2t) + 288*(7 + 12*z2)*pow4(
        Mt) - (131 + 12*z2)*pow4(Mst2)*pow4(s2t))*pow6(Mst1) - 90*pow2(-4*pow2(
        Mt) + pow2(Mst2)*pow2(s2t))*pow6(Mst2)) - pow4(Mst2)*(8*pow2(Mst1)*
        pow2(Mt)*(58800*Mst2*Mt*s2t + 220709*pow2(Mt) + 58635*pow2(Mst2)*pow2(
        s2t))*pow4(Mst2) + pow2(Mst2)*pow4(Mst1)*(-351480*pow2(Mst2)*pow2(Mt)*
        pow2(s2t) - 156800*Mst2*s2t*pow3(Mt) - 382200*Mt*pow3(Mst2)*pow3(s2t) +
        470400*pow4(Mt) + 51313*pow4(Mst2)*pow4(s2t)) + 8*(14865*pow2(Mst2)*
        pow2(Mt)*pow2(s2t) - 39200*Mst2*s2t*pow3(Mt) + 3675*Mt*pow3(Mst2)*pow3(
        s2t) + 49209*pow4(Mt) - 3459*pow4(Mst2)*pow4(s2t))*pow6(Mst1) + 411600*
        pow4(Mt)*pow6(Mst2)) + 196*pow2(Msq)*pow4(Mst2)*(pow4(Mst1)*(1256*pow2(
        Mst2)*pow2(Mt)*pow2(s2t) - 5600*Mst2*s2t*pow3(Mt) + 4800*Mt*pow3(Mst2)*
        pow3(s2t) - 656*pow4(Mt) - 345*pow4(Mst2)*pow4(s2t)) + pow2(Mst1)*pow2(
        Mst2)*(-5128*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 5600*Mst2*s2t*pow3(Mt) -
        4800*Mt*pow3(Mst2)*pow3(s2t) - 1856*pow4(Mt) + 645*pow4(Mst2)*pow4(s2t)
        ) + (872*pow2(Mt)*pow2(s2t) + 255*pow2(Mst2)*pow4(s2t))*pow6(Mst1) +
        3000*(-2*pow4(Mst2)*pow4(Mt) + pow2(Mt)*pow2(s2t)*pow6(Mst2)))) + 49*
        pow2(Dmglst2)*(-25*pow5(Mst2)*(pow2(Mst1)*pow3(Mst2)*(7272*pow2(Mst2)*
        pow2(Mt)*pow2(s2t) + 1024*Mst2*s2t*pow3(Mt) + 4328*Mt*pow3(Mst2)*pow3(
        s2t) + 12208*pow4(Mt) - 495*pow4(Mst2)*pow4(s2t)) + 3*Mst2*pow4(Mst1)*(
        -1320*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 256*Mst2*s2t*pow3(Mt) - 1448*Mt*
        pow3(Mst2)*pow3(s2t) + 384*pow4(Mt) + 53*pow4(Mst2)*pow4(s2t)) + 504*
        pow2(Mt)*(14*pow2(Mt) - 5*pow2(Mst2)*pow2(s2t))*pow5(Mst2) + 2*s2t*(72*
        Mst2*s2t*pow2(Mt) + 12*Mt*pow2(Mst2)*pow2(s2t) - 128*pow3(Mt) + 3*pow3(
        Mst2)*pow3(s2t))*pow6(Mst1)) + 2*pow4(Msq)*(-75*pow2(Mst2)*pow4(Mst1)*(
        -608*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 768*Mst2*s2t*pow3(Mt) - 16*Mt*(19
        + 6*z2)*pow3(Mst2)*pow3(s2t) + 96*(65 + 4*z2)*pow4(Mt) + (-221 + 108*
        z2)*pow4(Mst2)*pow4(s2t)) + pow2(Mst1)*pow4(Mst2)*(-54600*pow2(Mst2)*
        pow2(Mt)*pow2(s2t) - 800*Mst2*s2t*(-131 + 72*z2)*pow3(Mt) - 1200*Mt*(-
        197 + 138*z2)*pow3(Mst2)*pow3(s2t) + 3904*pow4(Mt) + 75*(-121 + 108*z2)
        *pow4(Mst2)*pow4(s2t)) + 1800*(-16*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 32*
        Mst2*s2t*(7 + 3*z2)*pow3(Mt) + 7*Mt*pow3(Mst2)*pow3(s2t) + 24*(23 + 12*
        z2)*pow4(Mt) + pow4(Mst2)*pow4(s2t))*pow6(Mst1) - 3375*pow2(-4*pow2(Mt)
        + pow2(Mst2)*pow2(s2t))*pow6(Mst2)) - 200*pow2(Msq)*pow5(Mst2)*(-(s2t*(
        600*Mst2*s2t*pow2(Mt) + 576*Mt*pow2(Mst2)*pow2(s2t) - 112*pow3(Mt) +
        21*pow3(Mst2)*pow3(s2t))*pow4(Mst1)) - 2*Mst2*pow2(Mst1)*(-600*pow2(
        Mst2)*pow2(Mt)*pow2(s2t) + 56*Mst2*s2t*pow3(Mt) - 288*Mt*pow3(Mst2)*
        pow3(s2t) - 456*pow4(Mt) + 51*pow4(Mst2)*pow4(s2t)) + 12*Mst2*pow4(s2t)
        *pow6(Mst1) + 15*(80*pow3(Mst2)*pow4(Mt) - 40*pow2(Mt)*pow2(s2t)*pow5(
        Mst2) + 3*pow4(s2t)*pow7(Mst2)))) + 44100*pow2(Mst1)*pow2(log(pow2(
        Mst1)/pow2(Mst2)))*pow4(Msq)*(pow2(Mgl)*(pow4(Mst1)*(64*Mst2*s2t*pow3(
        Mt) - 144*pow4(Mt) + pow4(Mst2)*pow4(s2t)) + pow4(Mst2)*(-8*pow2(Mst2)*
        pow2(Mt)*pow2(s2t) - 64*Mst2*s2t*pow3(Mt) + 16*Mt*pow3(Mst2)*pow3(s2t)
        + 16*pow4(Mt) + 7*pow4(Mst2)*pow4(s2t)) - 4*pow2(Mst1)*(8*pow2(Mst2)*
        pow4(Mt) - 4*Mt*pow3(s2t)*pow5(Mst2) + pow4(s2t)*pow6(Mst2))) + 4*
        Dmglst2*Mgl*(48*(3*Mt - Mst2*s2t)*pow3(Mt)*pow4(Mst1) - 16*s2t*pow3(Mt)
        *pow5(Mst2) + pow2(Mst1)*(16*pow2(Mst2)*pow4(Mt) + 4*Mt*pow3(s2t)*pow5(
        Mst2) - 3*pow4(s2t)*pow6(Mst2)) + 4*Mt*pow3(s2t)*pow7(Mst2) + 3*pow4(
        s2t)*pow8(Mst2)) + 2*pow2(Dmglst2)*(96*(3*Mt - Mst2*s2t)*pow3(Mt)*pow4(
        Mst1) - 32*s2t*pow3(Mt)*pow5(Mst2) + pow2(Mst1)*(-16*pow2(Mst2)*pow4(
        Mt) + 8*Mt*pow3(s2t)*pow5(Mst2) - 9*pow4(s2t)*pow6(Mst2)) + 8*Mt*pow3(
        s2t)*pow7(Mst2) + 9*pow4(s2t)*pow8(Mst2))) + 105*log(pow2(Mst1)/pow2(
        Mst2))*(pow2(Mgl)*(-1680*pow2(Mst1)*(pow2(Mst1) - pow2(Mst2))*pow4(
        Mst2)*pow4(s2t)*pow6(Msq) + pow2(Mst1)*(1680*pow3(Mst2)*(8*s2t*pow3(Mt)
        - Mt*pow2(Mst1)*pow3(s2t)) - 504*pow2(Mt)*pow2(s2t)*pow4(Mst1) + pow4(
        Mst2)*(296*pow2(Mt)*pow2(s2t) + 422*pow2(Mst1)*pow4(s2t)) + pow2(Mst2)*
        (264*pow2(Mst1)*pow2(Mt)*pow2(s2t) + 6304*pow4(Mt) - 187*pow4(Mst1)*
        pow4(s2t)))*pow6(Mst2) - 70*pow4(Msq)*(2*pow2(Mst2)*pow4(Mst1)*(-32*
        pow2(Mst2)*pow2(Mt)*pow2(s2t) - 384*Mst2*s2t*pow3(Mt) + 72*Mt*pow3(
        Mst2)*pow3(s2t) + 720*pow4(Mt) + 21*pow4(Mst2)*pow4(s2t)) - s2t*pow2(
        Mst1)*(128*Mst2*s2t*pow2(Mt) + 240*Mt*pow2(Mst2)*pow2(s2t) - 384*pow3(
        Mt) + 59*pow3(Mst2)*pow3(s2t))*pow5(Mst2) + (-96*pow2(Mst2)*pow2(Mt)*
        pow2(s2t) - 2112*Mst2*s2t*pow3(Mt) + 48*Mt*pow3(Mst2)*pow3(s2t) + 4608*
        pow4(Mt) + 26*pow4(Mst2)*pow4(s2t))*pow6(Mst1) - 18*pow2(-4*pow2(Mt) +
        pow2(Mst2)*pow2(s2t))*pow6(Mst2)) - 28*pow2(Msq)*pow2(Mst1)*pow4(Mst2)*
        (-160*pow3(Mst2)*(6*s2t*pow3(Mt) - Mt*pow2(Mst1)*pow3(s2t)) + 20*pow2(
        Mt)*pow2(s2t)*pow4(Mst1) - pow4(Mst2)*(52*pow2(Mt)*pow2(s2t) + 41*pow2(
        Mst1)*pow4(s2t)) + pow2(Mst2)*(32*pow2(Mst1)*pow2(Mt)*pow2(s2t) - 304*
        pow4(Mt) + 6*pow4(Mst1)*pow4(s2t)) + 160*Mt*pow3(s2t)*pow5(Mst2) + 40*
        pow4(s2t)*pow6(Mst2))) + 140*Dmglst2*Mgl*(-2*pow2(Msq)*pow2(Mst1)*(-
        224*Mst2*s2t*pow3(Mt) + 48*Mst2*Mt*(pow2(Mst1) + pow2(Mst2))*pow3(s2t)
        - 48*pow4(Mt) + (-20*pow2(Mst1)*pow2(Mst2) + 3*pow4(Mst1) + 17*pow4(
        Mst2))*pow4(s2t))*pow6(Mst2) + 2*pow4(Msq)*(pow2(Mst1)*pow4(Mst2)*(144*
        pow2(Mst2)*pow2(Mt)*pow2(s2t) + 160*Mst2*s2t*pow3(Mt) + 12*Mt*pow3(
        Mst2)*pow3(s2t) + 16*pow4(Mt) - 9*pow4(Mst2)*pow4(s2t)) - 3*pow2(Mst2)*
        pow4(Mst1)*(32*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 64*Mst2*s2t*pow3(Mt) +
        28*Mt*pow3(Mst2)*pow3(s2t) - 304*pow4(Mt) + 3*pow4(Mst2)*pow4(s2t)) +
        12*Mt*(-164*Mst2*s2t*pow2(Mt) - 10*Mt*pow2(Mst2)*pow2(s2t) + 456*pow3(
        Mt) + 3*pow3(Mst2)*pow3(s2t))*pow6(Mst1) + 9*pow2(-4*pow2(Mt) + pow2(
        Mst2)*pow2(s2t))*pow6(Mst2)) + pow2(Mst1)*(352*Mst2*s2t*pow3(Mt) - 60*
        Mst2*Mt*(pow2(Mst1) + pow2(Mst2))*pow3(s2t) + 128*pow4(Mt) + pow2(Mst1)
        *(-5*pow2(Mst1) + 21*pow2(Mst2))*pow4(s2t))*pow8(Mst2)) + 14*pow2(
        Dmglst2)*(-10*pow2(Msq)*pow2(Mst1)*(-832*Mst2*s2t*pow3(Mt) + 192*Mst2*
        Mt*(pow2(Mst1) + pow2(Mst2))*pow3(s2t) - 144*pow4(Mt) + (-100*pow2(
        Mst1)*pow2(Mst2) + 6*pow4(Mst1) + 91*pow4(Mst2))*pow4(s2t))*pow6(Mst2)
        + 2*pow4(Msq)*(pow2(Mst1)*pow4(Mst2)*(2160*pow2(Mst2)*pow2(Mt)*pow2(
        s2t) + 2560*Mst2*s2t*pow3(Mt) - 120*Mt*pow3(Mst2)*pow3(s2t) + 1408*
        pow4(Mt) - 345*pow4(Mst2)*pow4(s2t)) + 120*Mt*(-164*Mst2*s2t*pow2(Mt) -
        10*Mt*pow2(Mst2)*pow2(s2t) + 456*pow3(Mt) + 3*pow3(Mst2)*pow3(s2t))*
        pow6(Mst1) + 135*pow2(-4*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow6(Mst2) +
        15*pow4(Mst1)*(-64*pow2(Mt)*pow2(s2t)*pow4(Mst2) - 432*pow2(Mst2)*pow4(
        Mt) - 72*Mt*pow3(s2t)*pow5(Mst2) + 7*pow4(s2t)*pow6(Mst2))) - 5*pow2(
        Mst1)*(-1856*Mst2*s2t*pow3(Mt) + 360*Mst2*Mt*(pow2(Mst1) + pow2(Mst2))*
        pow3(s2t) - 544*pow4(Mt) + (-147*pow2(Mst1)*pow2(Mst2) + 10*pow4(Mst1)
        + 90*pow4(Mst2))*pow4(s2t))*pow8(Mst2)))))/(4410.*pow2(Mgl)*pow2(Mst1)*
        pow4(Msq)*pow4(Mst2)*pow4(Mt)) + (log(pow2(Mst1)/pow2(Mst2))*(4*
        Dmglst2*Mgl*(-55125*pow2(Mst1)*pow4(Mst2)*(4*Mt*pow2(Mst1)*(-200*Mst2*
        s2t*pow2(Mt) + 540*Mt*pow2(Mst2)*pow2(s2t) + 1088*pow3(Mt) + 425*pow3(
        Mst2)*pow3(s2t))*pow4(Mst2) - 3*pow2(Mst2)*pow4(Mst1)*(432*pow2(Mst2)*
        pow2(Mt)*pow2(s2t) + 64*Mst2*s2t*pow3(Mt) + 404*Mt*pow3(Mst2)*pow3(s2t)
        - 32*pow4(Mt) + pow4(Mst2)*pow4(s2t)) + (144*pow2(Mst2)*pow2(Mt)*pow2(
        s2t) + 480*Mst2*s2t*pow3(Mt) + 24*Mt*pow3(Mst2)*pow3(s2t) - 576*pow4(
        Mt) + 17*pow4(Mst2)*pow4(s2t))*pow6(Mst1) + 2016*pow4(Mt)*pow6(Mst2)) +
        2*pow4(Msq)*(-7350*pow4(Mst1)*pow4(Mst2)*(48*(-1676 + 891*z2)*pow2(
        Mst2)*pow2(Mt)*pow2(s2t) + 8*Mst2*s2t*(-22768 + 2835*S2 + 6462*z2)*
        pow3(Mt) - 18*Mt*(2990 + 315*S2 - 1986*z2)*pow3(Mst2)*pow3(s2t) + (
        90952 - 18144*S2 + 21888*z2)*pow4(Mt) + 9*(689 + 126*S2 - 580*z2)*pow4(
        Mst2)*pow4(s2t)) - 294*pow2(Mst2)*(-8*(59657 + 202500*S2 - 55800*z2)*
        pow2(Mst2)*pow2(Mt)*pow2(s2t) + 200*Mst2*s2t*(-18428 + 37341*S2 - 9090*
        z2)*pow3(Mt) - 300*Mt*(-818 + 3609*S2 - 594*z2)*pow3(Mst2)*pow3(s2t) -
        48*(8581 + 72450*S2 - 92700*z2)*pow4(Mt) + (-128893 + 291600*S2 +
        100350*z2)*pow4(Mst2)*pow4(s2t))*pow6(Mst1) + 66150*pow2(Mst1)*(-952*
        pow2(Mst2)*pow2(Mt)*pow2(s2t) - 1024*Mst2*s2t*pow3(Mt) + 256*Mt*pow3(
        Mst2)*pow3(s2t) + 880*pow4(Mt) + 247*pow4(Mst2)*pow4(s2t))*pow6(Mst2) +
        (72*(-5217307 + 33912900*S2 - 73500*z2)*pow2(Mst2)*pow2(Mt)*pow2(s2t) -
        78400*Mst2*s2t*(-3509 + 115785*S2 - 18144*z2)*pow3(Mt) + 4900*Mt*(24322
        + 140139*S2 - 38826*z2)*pow3(Mst2)*pow3(s2t) + 80*(57911521 + 49233240*
        S2 - 40325040*z2)*pow4(Mt) + (1011403 - 185175900*S2 + 17463600*z2)*
        pow4(Mst2)*pow4(s2t))*pow8(Mst1) + 16934400*pow2(Mt)*(-2*pow2(Mt) +
        pow2(Mst2)*pow2(s2t))*pow8(Mst2)) - 220500*pow2(Msq)*(-4*pow4(Mst2)*(
        120*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 168*Mst2*s2t*pow3(Mt) + 96*Mt*pow3(
        Mst2)*pow3(s2t) + 96*pow4(Mt) + 25*pow4(Mst2)*pow4(s2t))*pow6(Mst1) +
        pow4(Mst1)*(960*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 1184*Mst2*s2t*pow3(Mt)
        + 768*Mt*pow3(Mst2)*pow3(s2t) + 576*pow4(Mt) + 13*pow4(Mst2)*pow4(s2t))
        *pow6(Mst2) + 24*pow2(Mst2)*(-2*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 16*
        Mst2*s2t*pow3(Mt) + 120*pow4(Mt) + pow4(Mst2)*pow4(s2t))*pow8(Mst1) -
        480*pow2(Mst1)*pow2(Mt)*(-2*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow8(Mst2)
        )) + 3*pow2(Mgl)*(5292000*(-1 + 2*z2)*(pow2(Mst1) - pow2(Mst2))*pow4(
        Mst1)*pow4(Mst2)*pow4(s2t)*pow6(Msq) + 15*pow2(Mst1)*pow5(Mst2)*(-128*
        pow2(Mst1)*pow2(Mt)*(3675*Mst2*Mt*s2t + 73246*pow2(Mt) + 14279*pow2(
        Mst2)*pow2(s2t))*pow3(Mst2) + 2*Mst2*pow4(Mst1)*(851136*pow2(Mst2)*
        pow2(Mt)*pow2(s2t) + 1411200*Mst2*s2t*pow3(Mt) + 499800*Mt*pow3(Mst2)*
        pow3(s2t) - 1528800*pow4(Mt) - 45247*pow4(Mst2)*pow4(s2t)) - 1646400*
        pow4(Mt)*pow5(Mst2) + s2t*(-531552*Mst2*s2t*pow2(Mt) - 117600*Mt*pow2(
        Mst2)*pow2(s2t) + 1411200*pow3(Mt) + 34519*pow3(Mst2)*pow3(s2t))*pow6(
        Mst1)) + 2*pow4(Msq)*(2450*pow4(Mst1)*pow4(Mst2)*(-24*(-5434 + 1863*S2
        + 1830*z2 - 72*z3)*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 288*Mst2*s2t*(116 +
        189*S2 - 878*z2)*pow3(Mt) - 72*Mt*(-2326 + 189*S2 + 346*z2)*pow3(Mst2)*
        pow3(s2t) + 32*(722 + 1701*S2 + 3753*z2 + 108*z3)*pow4(Mt) + 3*(6377 +
        2592*S2 + 2724*z2 - 72*z3)*pow4(Mst2)*pow4(s2t)) + 1960*pow2(Mst2)*(-4*
        (-47993 + 86940*S2 + 3060*z2 + 540*z3)*pow2(Mst2)*pow2(Mt)*pow2(s2t) +
        1080*Mst2*s2t*(102 + 377*S2 - 330*z2)*pow3(Mt) - 180*Mt*(326 + 219*S2 -
        294*z2)*pow3(Mst2)*pow3(s2t) + 8*(-37598 + 28215*S2 + 30870*z2)*pow4(
        Mt) + (-18499 + 46305*S2 - 15075*z2)*pow4(Mst2)*pow4(s2t))*pow6(Mst1) +
        44100*pow2(Mst1)*(-2312*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 2048*Mst2*s2t*
        pow3(Mt) + 512*Mt*pow3(Mst2)*pow3(s2t) + 4112*pow4(Mt) + 385*pow4(Mst2)
        *pow4(s2t))*pow6(Mst2) + (-24*(-20192173 + 50839950*S2 + 2866500*z2)*
        pow2(Mst2)*pow2(Mt)*pow2(s2t) + 940800*Mst2*s2t*(302 + 2055*S2 - 732*
        z2)*pow3(Mt) - 58800*Mt*(136 + 1317*S2 - 918*z2)*pow3(Mst2)*pow3(s2t) +
        32*(-47778491 + 25754400*S2 + 20374200*z2)*pow4(Mt) + 9*(-5929019 +
        9530500*S2 + 367500*z2 + 58800*z3)*pow4(Mst2)*pow4(s2t))*pow8(Mst1) +
        11289600*pow2(Mt)*(-2*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow8(Mst2)) +
        5880*pow2(Msq)*pow2(Mst1)*pow2(Mst2)*(4*pow2(Mst2)*pow4(Mst1)*(1364*
        pow2(Mst2)*pow2(Mt)*pow2(s2t) + 1200*Mst2*s2t*pow3(Mt) + 800*Mt*pow3(
        Mst2)*pow3(s2t) + 600*pow4(Mt) + 173*pow4(Mst2)*pow4(s2t)) - pow2(Mst1)
        *pow4(Mst2)*(8656*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 33600*Mst2*s2t*pow3(
        Mt) + 16000*Mt*pow3(Mst2)*pow3(s2t) + 27712*pow4(Mt) + 735*pow4(Mst2)*
        pow4(s2t)) + 2*(-2000*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 9600*Mst2*s2t*
        pow3(Mt) + 28800*pow4(Mt) + 159*pow4(Mst2)*pow4(s2t))*pow6(Mst1) +
        6000*(-2*pow4(Mt)*pow6(Mst2) + pow2(Mt)*pow2(s2t)*pow8(Mst2)))) + 2*
        pow2(Dmglst2)*(4*pow4(Msq)*(-147*pow4(Mst1)*pow4(Mst2)*(-1200*(5987 +
        567*S2 - 4935*z2)*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 1200*Mst2*s2t*(5368 +
        63*S2 - 2274*z2)*pow3(Mt) + 300*Mt*(-17774 + 63*S2 + 12978*z2)*pow3(
        Mst2)*pow3(s2t) + 16*(-350909 + 122850*S2 + 66600*z2)*pow4(Mt) + 75*(
        8815 + 630*S2 - 4944*z2)*pow4(Mst2)*pow4(s2t)) - 147*pow2(Mst2)*(-8*(
        548953 + 805950*S2 - 632700*z2)*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 3600*
        Mst2*s2t*(-2725 + 4307*S2 - 182*z2)*pow3(Mt) + 400*Mt*(-4672 + 6903*S2
        + 2466*z2)*pow3(Mst2)*pow3(s2t) + 16*(-1633879 + 2116350*S2 - 242100*
        z2)*pow4(Mt) + 3*(-206449 + 24000*S2 + 96750*z2)*pow4(Mst2)*pow4(s2t))*
        pow6(Mst1) + 33075*pow2(Mst1)*(472*pow2(Mst2)*pow2(Mt)*pow2(s2t) +
        2048*Mst2*s2t*pow3(Mt) - 512*Mt*pow3(Mst2)*pow3(s2t) - 4016*pow4(Mt) +
        197*pow4(Mst2)*pow4(s2t))*pow6(Mst2) + (72*(-5217307 + 33912900*S2 -
        73500*z2)*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 78400*Mst2*s2t*(-3509 +
        115785*S2 - 18144*z2)*pow3(Mt) + 4900*Mt*(24322 + 140139*S2 - 38826*z2)
        *pow3(Mst2)*pow3(s2t) + 80*(57911521 + 49233240*S2 - 40325040*z2)*pow4(
        Mt) + (1011403 - 185175900*S2 + 17463600*z2)*pow4(Mst2)*pow4(s2t))*
        pow8(Mst1) + 1058400*(8*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 16*pow4(Mt) +
        pow4(Mst2)*pow4(s2t))*pow8(Mst2)) - 44100*pow2(Msq)*(-20*pow4(Mst2)*(
        600*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 528*Mst2*s2t*pow3(Mt) + 528*Mt*
        pow3(Mst2)*pow3(s2t) + 96*pow4(Mt) + 101*pow4(Mst2)*pow4(s2t))*pow6(
        Mst1) + pow4(Mst1)*(24000*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 9920*Mst2*
        s2t*pow3(Mt) + 12480*Mt*pow3(Mst2)*pow3(s2t) + 10656*pow4(Mt) - 305*
        pow4(Mst2)*pow4(s2t))*pow6(Mst2) + 240*pow2(Mst2)*(-2*pow2(Mst2)*pow2(
        Mt)*pow2(s2t) - 16*Mst2*s2t*pow3(Mt) + 120*pow4(Mt) + pow4(Mst2)*pow4(
        s2t))*pow8(Mst1) + 300*pow2(Mst1)*(-40*pow2(Mst2)*pow2(Mt)*pow2(s2t) +
        80*pow4(Mt) + 3*pow4(Mst2)*pow4(s2t))*pow8(Mst2)) - 11025*pow2(Mst1)*
        pow4(Mst2)*(15*pow2(Mst2)*pow4(Mst1)*(-2640*pow2(Mst2)*pow2(Mt)*pow2(
        s2t) + 1280*Mst2*s2t*pow3(Mt) - 2872*Mt*pow3(Mst2)*pow3(s2t) - 672*
        pow4(Mt) + pow4(Mst2)*pow4(s2t)) + 10*(144*pow2(Mst2)*pow2(Mt)*pow2(
        s2t) + 480*Mst2*s2t*pow3(Mt) + 24*Mt*pow3(Mst2)*pow3(s2t) - 576*pow4(
        Mt) + 17*pow4(Mst2)*pow4(s2t))*pow6(Mst1) - 5040*(-14*pow4(Mt)*pow6(
        Mst2) + 5*pow2(Mt)*pow2(s2t)*pow8(Mst2)) + pow2(Mst1)*(94016*pow4(Mst2)
        *pow4(Mt) - 9920*s2t*pow3(Mt)*pow5(Mst2) + 72720*pow2(Mt)*pow2(s2t)*
        pow6(Mst2) + 43640*Mt*pow3(s2t)*pow7(Mst2) - 4110*pow4(s2t)*pow8(Mst2))
        ))))/(793800.*pow2(Mgl)*pow4(Msq)*pow4(Mst1)*pow4(Mst2)*pow4(Mt));

   return result;
}

/**
 * @return returns the susy log^1 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double H6g2::calc_coef_at_as2_no_sm_logs_log1() const {

   const double result =
      -474.3703703703704 - (115472*Dmglst2)/(27.*Mgl) + (1472*Mst2*s2t)/(3.*Mt)
        + (204928*Dmglst2*Mst2*s2t)/(27.*Mgl*Mt) + (4304*z2)/3. - (1536*
        Dmglst2*z2)/Mgl - (7552*Mst2*s2t*z2)/(3.*Mt) - (6656*Dmglst2*Mst2*s2t*
        z2)/(3.*Mgl*Mt) + 128*z3 + (3701984*pow2(Dmglst2))/(675.*pow2(Mgl)) + (
        51968*Mst2*s2t*pow2(Dmglst2))/(9.*Mt*pow2(Mgl)) - (5312*z2*pow2(
        Dmglst2))/(3.*pow2(Mgl)) - (2624*Mst2*s2t*z2*pow2(Dmglst2))/(Mt*pow2(
        Mgl)) - (320*pow2(Msq))/pow2(Mst1) + (10496*s2t*pow2(Mst1))/(3.*Mst2*
        Mt) + (11648*Dmglst2*s2t*pow2(Mst1))/(3.*Mgl*Mst2*Mt) - (8192*s2t*z2*
        pow2(Mst1))/(3.*Mst2*Mt) + (9088*Dmglst2*s2t*z2*pow2(Mst1))/(3.*Mgl*
        Mst2*Mt) - (67456*s2t*pow2(Dmglst2)*pow2(Mst1))/(9.*Mst2*Mt*pow2(Mgl))
        - (960*s2t*z2*pow2(Dmglst2)*pow2(Mst1))/(Mst2*Mt*pow2(Mgl)) + (160*
        pow2(Mst1))/(9.*pow2(Msq)) + (640*Dmglst2*pow2(Mst1))/(3.*Mgl*pow2(Msq)
        ) + (1600*Mst2*s2t*pow2(Mst1))/(9.*Mt*pow2(Msq)) - (2240*Dmglst2*Mst2*
        s2t*pow2(Mst1))/(9.*Mgl*Mt*pow2(Msq)) + (320*pow2(Dmglst2)*pow2(Mst1))/
        (3.*pow2(Mgl)*pow2(Msq)) - (4160*Mst2*s2t*pow2(Dmglst2)*pow2(Mst1))/(9.
        *Mt*pow2(Mgl)*pow2(Msq)) - (320*pow2(Msq))/pow2(Mst2) - (6208*pow2(
        Mst1))/(3.*pow2(Mst2)) + (8192*Dmglst2*pow2(Mst1))/(9.*Mgl*pow2(Mst2))
        + (1600*z2*pow2(Mst1))/pow2(Mst2) - (8960*Dmglst2*z2*pow2(Mst1))/(Mgl*
        pow2(Mst2)) + (612496*pow2(Dmglst2)*pow2(Mst1))/(27.*pow2(Mgl)*pow2(
        Mst2)) + (8960*z2*pow2(Dmglst2)*pow2(Mst1))/(3.*pow2(Mgl)*pow2(Mst2)) -
        (1040*pow2(Mst2))/(3.*pow2(Msq)) + (160*Dmglst2*pow2(Mst2))/(3.*Mgl*
        pow2(Msq)) + (1808*pow2(Dmglst2)*pow2(Mst2))/(3.*pow2(Mgl)*pow2(Msq)) +
        (3424*pow2(Mst2))/(3.*pow2(Mst1)) + (2240*Dmglst2*pow2(Mst2))/(3.*Mgl*
        pow2(Mst1)) - (480*pow2(Dmglst2)*pow2(Mst2))/(pow2(Mgl)*pow2(Mst1)) - (
        320*pow2(Msq)*pow2(s2t))/pow2(Mt) + (304*pow2(Mst1)*pow2(s2t))/pow2(Mt)
        + (31280*Dmglst2*pow2(Mst1)*pow2(s2t))/(9.*Mgl*pow2(Mt)) - (2048*
        Dmglst2*z2*pow2(Mst1)*pow2(s2t))/(Mgl*pow2(Mt)) + (94504*pow2(Dmglst2)*
        pow2(Mst1)*pow2(s2t))/(9.*pow2(Mgl)*pow2(Mt)) - (7168*z2*pow2(Dmglst2)*
        pow2(Mst1)*pow2(s2t))/(pow2(Mgl)*pow2(Mt)) + (160*pow2(Msq)*pow2(Mst1)*
        pow2(s2t))/(pow2(Mst2)*pow2(Mt)) + (8080*pow2(Mst2)*pow2(s2t))/(3.*
        pow2(Mt)) + (57040*Dmglst2*pow2(Mst2)*pow2(s2t))/(9.*Mgl*pow2(Mt)) - (
        512*z2*pow2(Mst2)*pow2(s2t))/pow2(Mt) - (3072*Dmglst2*z2*pow2(Mst2)*
        pow2(s2t))/(Mgl*pow2(Mt)) + (114008*pow2(Dmglst2)*pow2(Mst2)*pow2(s2t))
        /(9.*pow2(Mgl)*pow2(Mt)) - (8704*z2*pow2(Dmglst2)*pow2(Mst2)*pow2(s2t))
        /(pow2(Mgl)*pow2(Mt)) + (160*pow2(Msq)*pow2(Mst2)*pow2(s2t))/(pow2(
        Mst1)*pow2(Mt)) - (640*s2t*pow3(Mst2))/(9.*Mt*pow2(Msq)) + (640*
        Dmglst2*s2t*pow3(Mst2))/(3.*Mgl*Mt*pow2(Msq)) + (1280*s2t*pow2(Dmglst2)
        *pow3(Mst2))/(3.*Mt*pow2(Mgl)*pow2(Msq)) - (2048*s2t*pow3(Mst2))/(3.*
        Mt*pow2(Mst1)) - (2048*Dmglst2*s2t*pow3(Mst2))/(3.*Mgl*Mt*pow2(Mst1)) +
        (2048*s2t*pow2(Dmglst2)*pow3(Mst2))/(3.*Mt*pow2(Mgl)*pow2(Mst1)) - (
        4736*Mst2*pow2(Mst1)*pow3(s2t))/(3.*pow3(Mt)) - (1184*Dmglst2*Mst2*
        pow2(Mst1)*pow3(s2t))/(Mgl*pow3(Mt)) + (544*Mst2*z2*pow2(Mst1)*pow3(
        s2t))/(3.*pow3(Mt)) + (544*Dmglst2*Mst2*z2*pow2(Mst1)*pow3(s2t))/(3.*
        Mgl*pow3(Mt)) - (1328*Mst2*pow2(Dmglst2)*pow2(Mst1)*pow3(s2t))/(9.*
        pow2(Mgl)*pow3(Mt)) + (544*Mst2*z2*pow2(Dmglst2)*pow2(Mst1)*pow3(s2t))/
        (3.*pow2(Mgl)*pow3(Mt)) + (4832*pow3(Mst2)*pow3(s2t))/(3.*pow3(Mt)) + (
        2048*Dmglst2*pow3(Mst2)*pow3(s2t))/(Mgl*pow3(Mt)) - (608*z2*pow3(Mst2)*
        pow3(s2t))/(3.*pow3(Mt)) - (2912*Dmglst2*z2*pow3(Mst2)*pow3(s2t))/(3.*
        Mgl*pow3(Mt)) + (32432*pow2(Dmglst2)*pow3(Mst2)*pow3(s2t))/(9.*pow2(
        Mgl)*pow3(Mt)) - (6368*z2*pow2(Dmglst2)*pow3(Mst2)*pow3(s2t))/(3.*pow2(
        Mgl)*pow3(Mt)) - (100*pow2(Mst1)*pow2(Mst2))/(3.*pow4(Msq)) + (280*
        Dmglst2*pow2(Mst1)*pow2(Mst2))/(3.*Mgl*pow4(Msq)) + (300*pow2(Dmglst2)*
        pow2(Mst1)*pow2(Mst2))/(pow2(Mgl)*pow4(Msq)) + (560*s2t*pow2(Mst1)*
        pow3(Mst2))/(9.*Mt*pow4(Msq)) - (80*Dmglst2*s2t*pow2(Mst1)*pow3(Mst2))/
        (3.*Mgl*Mt*pow4(Msq)) - (1120*s2t*pow2(Dmglst2)*pow2(Mst1)*pow3(Mst2))/
        (3.*Mt*pow2(Mgl)*pow4(Msq)) - (640*s2t*pow4(Mst1))/(3.*Mst2*Mt*pow2(
        Msq)) + (640*Dmglst2*s2t*pow4(Mst1))/(3.*Mgl*Mst2*Mt*pow2(Msq)) + (640*
        s2t*pow2(Dmglst2)*pow4(Mst1))/(3.*Mst2*Mt*pow2(Mgl)*pow2(Msq)) + (1120*
        pow4(Mst1))/(3.*pow2(Msq)*pow2(Mst2)) - (3200*Dmglst2*pow4(Mst1))/(3.*
        Mgl*pow2(Msq)*pow2(Mst2)) - (3200*pow2(Dmglst2)*pow4(Mst1))/(3.*pow2(
        Mgl)*pow2(Msq)*pow2(Mst2)) + (144*pow2(s2t)*pow4(Mst1))/(pow2(Mst2)*
        pow2(Mt)) + (9856*Dmglst2*pow2(s2t)*pow4(Mst1))/(3.*Mgl*pow2(Mst2)*
        pow2(Mt)) - (2048*Dmglst2*z2*pow2(s2t)*pow4(Mst1))/(Mgl*pow2(Mst2)*
        pow2(Mt)) + (9856*pow2(Dmglst2)*pow2(s2t)*pow4(Mst1))/(3.*pow2(Mgl)*
        pow2(Mst2)*pow2(Mt)) - (2048*z2*pow2(Dmglst2)*pow2(s2t)*pow4(Mst1))/(
        pow2(Mgl)*pow2(Mst2)*pow2(Mt)) + (13024*s2t*pow4(Mst1))/(3.*Mt*pow3(
        Mst2)) - (5248*Dmglst2*s2t*pow4(Mst1))/(Mgl*Mt*pow3(Mst2)) - (7936*s2t*
        z2*pow4(Mst1))/(3.*Mt*pow3(Mst2)) + (24704*Dmglst2*s2t*z2*pow4(Mst1))/(
        3.*Mgl*Mt*pow3(Mst2)) - (5248*s2t*pow2(Dmglst2)*pow4(Mst1))/(Mt*pow2(
        Mgl)*pow3(Mst2)) + (24704*s2t*z2*pow2(Dmglst2)*pow4(Mst1))/(3.*Mt*pow2(
        Mgl)*pow3(Mst2)) - (392*pow3(s2t)*pow4(Mst1))/(3.*Mst2*pow3(Mt)) + (
        184*Dmglst2*pow3(s2t)*pow4(Mst1))/(3.*Mgl*Mst2*pow3(Mt)) + (184*pow2(
        Dmglst2)*pow3(s2t)*pow4(Mst1))/(3.*Mst2*pow2(Mgl)*pow3(Mt)) + (110*
        pow4(Mst1))/(3.*pow4(Msq)) + (80*Dmglst2*pow4(Mst1))/(Mgl*pow4(Msq)) +
        (40*Mst2*s2t*pow4(Mst1))/(9.*Mt*pow4(Msq)) - (920*Dmglst2*Mst2*s2t*
        pow4(Mst1))/(9.*Mgl*Mt*pow4(Msq)) + (80*pow2(Dmglst2)*pow4(Mst1))/(
        pow2(Mgl)*pow4(Msq)) - (920*Mst2*s2t*pow2(Dmglst2)*pow4(Mst1))/(9.*Mt*
        pow2(Mgl)*pow4(Msq)) - (11092*pow4(Mst1))/(3.*pow4(Mst2)) + (158720*
        Dmglst2*pow4(Mst1))/(9.*Mgl*pow4(Mst2)) + (5632*z2*pow4(Mst1))/(3.*
        pow4(Mst2)) - (20224*Dmglst2*z2*pow4(Mst1))/(Mgl*pow4(Mst2)) + (158720*
        pow2(Dmglst2)*pow4(Mst1))/(9.*pow2(Mgl)*pow4(Mst2)) - (20224*z2*pow2(
        Dmglst2)*pow4(Mst1))/(pow2(Mgl)*pow4(Mst2)) - (656*pow2(s2t)*pow4(Mst2)
        )/(pow2(Mst1)*pow2(Mt)) - (2144*Dmglst2*pow2(s2t)*pow4(Mst2))/(3.*Mgl*
        pow2(Mst1)*pow2(Mt)) - (272*pow2(Dmglst2)*pow2(s2t)*pow4(Mst2))/(pow2(
        Mgl)*pow2(Mst1)*pow2(Mt)) - (770*pow4(Mst2))/(9.*pow4(Msq)) + (520*
        Dmglst2*pow4(Mst2))/(9.*Mgl*pow4(Msq)) + (5716*pow2(Dmglst2)*pow4(Mst2)
        )/(9.*pow2(Mgl)*pow4(Msq)) - (512*pow4(Mst2))/(3.*pow4(Mst1)) - (1024*
        Dmglst2*pow4(Mst2))/(3.*Mgl*pow4(Mst1)) - (512*pow2(Dmglst2)*pow4(Mst2)
        )/(3.*pow2(Mgl)*pow4(Mst1)) + (20*pow2(Msq)*pow2(Mst1)*pow4(s2t))/pow4(
        Mt) + (20*pow2(Msq)*pow2(Mst2)*pow4(s2t))/pow4(Mt) - (1162*pow2(Mst1)*
        pow2(Mst2)*pow4(s2t))/(3.*pow4(Mt)) + (4*Dmglst2*pow2(Mst1)*pow2(Mst2)*
        pow4(s2t))/(3.*Mgl*pow4(Mt)) - (172*z2*pow2(Mst1)*pow2(Mst2)*pow4(s2t))
        /(3.*pow4(Mt)) - (344*Dmglst2*z2*pow2(Mst1)*pow2(Mst2)*pow4(s2t))/(3.*
        Mgl*pow4(Mt)) + (830*pow2(Dmglst2)*pow2(Mst1)*pow2(Mst2)*pow4(s2t))/(9.
        *pow2(Mgl)*pow4(Mt)) - (172*z2*pow2(Dmglst2)*pow2(Mst1)*pow2(Mst2)*
        pow4(s2t))/(pow2(Mgl)*pow4(Mt)) + (1667*pow4(Mst1)*pow4(s2t))/(6.*pow4(
        Mt)) + (98*Dmglst2*pow4(Mst1)*pow4(s2t))/(Mgl*pow4(Mt)) + (82*z2*pow4(
        Mst1)*pow4(s2t))/(3.*pow4(Mt)) + (98*pow2(Dmglst2)*pow4(Mst1)*pow4(s2t)
        )/(pow2(Mgl)*pow4(Mt)) - (20*pow2(Msq)*pow4(Mst1)*pow4(s2t))/(pow2(
        Mst2)*pow4(Mt)) + (28*pow4(Mst2)*pow4(s2t))/pow4(Mt) - (844*Dmglst2*
        pow4(Mst2)*pow4(s2t))/(3.*Mgl*pow4(Mt)) + (30*z2*pow4(Mst2)*pow4(s2t))/
        pow4(Mt) + (344*Dmglst2*z2*pow4(Mst2)*pow4(s2t))/(3.*Mgl*pow4(Mt)) - (
        3352*pow2(Dmglst2)*pow4(Mst2)*pow4(s2t))/(9.*pow2(Mgl)*pow4(Mt)) + (
        172*z2*pow2(Dmglst2)*pow4(Mst2)*pow4(s2t))/(pow2(Mgl)*pow4(Mt)) - (20*
        pow2(Msq)*pow4(Mst2)*pow4(s2t))/(pow2(Mst1)*pow4(Mt)) + (pow2(log(pow2(
        Mst1)/pow2(Mst2)))*(-8*Dmglst2*Mgl*(8*Mt*(-67*Mst2*s2t*pow2(Mt) - 108*
        Mt*pow2(Mst2)*pow2(s2t) + 128*pow3(Mt) - 16*pow3(Mst2)*pow3(s2t))*pow4(
        Mst1) + pow4(Mst2)*(-320*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 848*Mst2*s2t*
        pow3(Mt) - 532*Mt*pow3(Mst2)*pow3(s2t) + 512*pow4(Mt) - 99*pow4(Mst2)*
        pow4(s2t)) + pow2(Mst1)*pow2(Mst2)*(-768*pow2(Mst2)*pow2(Mt)*pow2(s2t)
        - 448*Mst2*s2t*pow3(Mt) - 212*Mt*pow3(Mst2)*pow3(s2t) + 1656*pow4(Mt) +
        35*pow4(Mst2)*pow4(s2t))) - 4*pow2(Dmglst2)*(16*Mt*(-67*Mst2*s2t*pow2(
        Mt) - 108*Mt*pow2(Mst2)*pow2(s2t) + 128*pow3(Mt) - 16*pow3(Mst2)*pow3(
        s2t))*pow4(Mst1) + pow4(Mst2)*(-960*pow2(Mst2)*pow2(Mt)*pow2(s2t) +
        1696*Mst2*s2t*pow3(Mt) - 1832*Mt*pow3(Mst2)*pow3(s2t) + 1536*pow4(Mt) -
        297*pow4(Mst2)*pow4(s2t)) + pow2(Mst1)*pow2(Mst2)*(-2304*pow2(Mst2)*
        pow2(Mt)*pow2(s2t) + 256*Mst2*s2t*pow3(Mt) - 424*Mt*pow3(Mst2)*pow3(
        s2t) + 392*pow4(Mt) + 105*pow4(Mst2)*pow4(s2t))) + pow2(Mgl)*(8*pow2(
        Mst1)*pow2(Mst2)*(425*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 704*Mst2*s2t*
        pow3(Mt) + 212*Mt*pow3(Mst2)*pow3(s2t) - 196*pow4(Mt) + 3*pow4(Mst2)*
        pow4(s2t)) + pow4(Mst1)*(2688*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 3136*
        Mst2*s2t*pow3(Mt) + 1024*Mt*pow3(Mst2)*pow3(s2t) - 7168*pow4(Mt) + 41*
        pow4(Mst2)*pow4(s2t)) + pow4(Mst2)*(296*pow2(Mst2)*pow2(Mt)*pow2(s2t) -
        6784*Mst2*s2t*pow3(Mt) + 2208*Mt*pow3(Mst2)*pow3(s2t) + 288*pow4(Mt) +
        519*pow4(Mst2)*pow4(s2t)))))/(3.*pow2(Mgl)*pow4(Mst2)*pow4(Mt)) + (512*
        pow3(s2t)*pow5(Mst2))/(3.*pow2(Mst1)*pow3(Mt)) + (512*Dmglst2*pow3(s2t)
        *pow5(Mst2))/(3.*Mgl*pow2(Mst1)*pow3(Mt)) - (512*pow2(Dmglst2)*pow3(
        s2t)*pow5(Mst2))/(3.*pow2(Mgl)*pow2(Mst1)*pow3(Mt)) - (40*s2t*pow5(
        Mst2))/(Mt*pow4(Msq)) + (120*Dmglst2*s2t*pow5(Mst2))/(Mgl*Mt*pow4(Msq))
        + (520*s2t*pow2(Dmglst2)*pow5(Mst2))/(Mt*pow2(Mgl)*pow4(Msq)) + (8*log(
        pow2(Mst1)/pow2(Msq))*(-60*log(pow2(Mst1)/pow2(Mst2))*pow4(Msq)*(4*
        Dmglst2*Mgl*pow2(Mst1)*((9*Mt - 3*Mst2*s2t)*pow2(Mst1) + Mt*pow2(Mst2))
        - 2*pow2(Dmglst2)*(Mt*pow2(Mst1)*pow2(Mst2) + 6*(-3*Mt + Mst2*s2t)*
        pow4(Mst1)) + pow2(Mgl)*(-2*Mt*pow2(Mst1)*pow2(Mst2) + (-9*Mt + 4*Mst2*
        s2t)*pow4(Mst1) + 2*Mt*pow4(Mst2))) + 5*pow2(Mgl)*(8*pow2(Msq)*((Mt -
        2*Mst2*s2t)*pow2(Mst1) + 2*(2*Mt + Mst2*s2t)*pow2(Mst2))*pow4(Mst2) +
        6*pow4(Msq)*(2*(5*Mt - 4*Mst2*s2t)*pow2(Mst1)*pow2(Mst2) + 2*(9*Mt - 4*
        Mst2*s2t)*pow4(Mst1) + (9*Mt + 4*Mst2*s2t)*pow4(Mst2)) + pow4(Mst2)*(
        pow2(Mst1)*(6*Mt*pow2(Mst2) - 8*s2t*pow3(Mst2)) + (3*Mt - 4*Mst2*s2t)*
        pow4(Mst1) + (19*Mt + 12*Mst2*s2t)*pow4(Mst2))) + 2*pow2(Dmglst2)*(
        pow4(Msq)*(330*Mt*pow2(Mst1)*pow2(Mst2) - 60*(27*Mt - 10*Mst2*s2t)*
        pow4(Mst1) + 4*(16*Mt - 5*Mst2*s2t)*pow4(Mst2)) + 20*pow2(Msq)*(9*Mst2*
        Mt - 2*s2t*pow2(Mst1) + 2*s2t*pow2(Mst2))*pow5(Mst2) + 5*(3*Mst2*(3*Mt
        - 8*Mst2*s2t)*pow2(Mst1) + (59*Mt + 26*Mst2*s2t)*pow3(Mst2) - 2*s2t*
        pow4(Mst1))*pow5(Mst2)) - 20*Dmglst2*Mgl*(2*pow4(Msq)*(3*(7*Mt - 2*
        Mst2*s2t)*pow2(Mst1)*pow2(Mst2) + (81*Mt - 30*Mst2*s2t)*pow4(Mst1) + (-
        5*Mt + Mst2*s2t)*pow4(Mst2)) - 4*pow2(Msq)*(3*Mst2*Mt - s2t*pow2(Mst1)
        + s2t*pow2(Mst2))*pow5(Mst2) + (3*Mst2*(-Mt + 2*Mst2*s2t)*pow2(Mst1) -
        (13*Mt + 7*Mst2*s2t)*pow3(Mst2) + s2t*pow4(Mst1))*pow5(Mst2))))/(3.*Mt*
        pow2(Mgl)*pow4(Msq)*pow4(Mst2)) + (256*pow2(s2t)*pow6(Mst2))/(3.*pow2(
        Mt)*pow4(Mst1)) + (512*Dmglst2*pow2(s2t)*pow6(Mst2))/(3.*Mgl*pow2(Mt)*
        pow4(Mst1)) + (256*pow2(Dmglst2)*pow2(s2t)*pow6(Mst2))/(3.*pow2(Mgl)*
        pow2(Mt)*pow4(Mst1)) + (114*pow4(s2t)*pow6(Mst2))/(pow2(Mst1)*pow4(Mt))
        + (524*Dmglst2*pow4(s2t)*pow6(Mst2))/(3.*Mgl*pow2(Mst1)*pow4(Mt)) + (
        358*pow2(Dmglst2)*pow4(s2t)*pow6(Mst2))/(3.*pow2(Mgl)*pow2(Mst1)*pow4(
        Mt)) + (32*pow2(Dmglst2)*pow4(s2t)*pow8(Mst2))/(3.*pow2(Mgl)*pow4(Mst1)
        *pow4(Mt)) - (log(pow2(Mst1)/pow2(Mst2))*(20*Dmglst2*Mgl*(240*pow2(Msq)
        *pow2(Mst2)*pow3(Mt)*pow4(Mst1)*(-2*s2t*pow2(Mst1)*pow3(Mst2) + 3*Mt*
        pow4(Mst1) + 2*(3*Mt + Mst2*s2t)*pow4(Mst2)) - 120*pow3(Mt)*pow4(Mst1)*
        (3*Mst2*(-Mt + 2*Mst2*s2t)*pow2(Mst1) - (13*Mt + 7*Mst2*s2t)*pow3(Mst2)
        + s2t*pow4(Mst1))*pow5(Mst2) + pow4(Msq)*(-4*pow4(Mst1)*pow4(Mst2)*(
        1860*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 1532*Mst2*s2t*pow3(Mt) + 774*Mt*
        pow3(Mst2)*pow3(s2t) - 2764*pow4(Mt) + 15*pow4(Mst2)*pow4(s2t)) + 2*
        pow2(Mst2)*(-432*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 10936*Mst2*s2t*pow3(
        Mt) + 996*Mt*pow3(Mst2)*pow3(s2t) + 16228*pow4(Mt) + 207*pow4(Mst2)*
        pow4(s2t))*pow6(Mst1) - 3*pow2(Mst1)*(-856*pow2(Mst2)*pow2(Mt)*pow2(
        s2t) - 1536*Mst2*s2t*pow3(Mt) + 384*Mt*pow3(Mst2)*pow3(s2t) + 1200*
        pow4(Mt) + 203*pow4(Mst2)*pow4(s2t))*pow6(Mst2) + 4*(-738*pow2(Mst2)*
        pow2(Mt)*pow2(s2t) - 8204*Mst2*s2t*pow3(Mt) + 135*Mt*pow3(Mst2)*pow3(
        s2t) + 12492*pow4(Mt) + 12*pow4(Mst2)*pow4(s2t))*pow8(Mst1) - 768*pow2(
        Mt)*(-2*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow8(Mst2))) - 5*pow2(Mgl)*(
        480*pow2(Msq)*pow2(Mst2)*pow3(Mt)*pow4(Mst1)*(2*(-Mt + 2*Mst2*s2t)*
        pow2(Mst1)*pow2(Mst2) + 3*Mt*pow4(Mst1) - 4*(2*Mt + Mst2*s2t)*pow4(
        Mst2)) + 120*pow3(Mt)*pow4(Mst1)*pow4(Mst2)*(2*(-3*Mt + 4*Mst2*s2t)*
        pow2(Mst1)*pow2(Mst2) + (-3*Mt + 4*Mst2*s2t)*pow4(Mst1) - (19*Mt + 12*
        Mst2*s2t)*pow4(Mst2)) + 360*(pow2(Mst1) - pow2(Mst2))*pow4(Mst1)*pow4(
        Mst2)*pow4(s2t)*pow6(Msq) + pow4(Msq)*(4*pow4(Mst1)*pow4(Mst2)*(3372*
        pow2(Mst2)*pow2(Mt)*pow2(s2t) - 9744*Mst2*s2t*pow3(Mt) + 3912*Mt*pow3(
        Mst2)*pow3(s2t) + 4492*pow4(Mt) + 669*pow4(Mst2)*pow4(s2t)) - 48*pow2(
        Mst2)*(-267*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 792*Mst2*s2t*pow3(Mt) - 30*
        Mt*pow3(Mst2)*pow3(s2t) - 317*pow4(Mt) + 56*pow4(Mst2)*pow4(s2t))*pow6(
        Mst1) + 6*pow2(Mst1)*(-728*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 1024*Mst2*
        s2t*pow3(Mt) + 256*Mt*pow3(Mst2)*pow3(s2t) + 1200*pow4(Mt) + 139*pow4(
        Mst2)*pow4(s2t))*pow6(Mst2) + 3*(4640*pow2(Mst2)*pow2(Mt)*pow2(s2t) -
        13376*Mst2*s2t*pow3(Mt) + 240*Mt*pow3(Mst2)*pow3(s2t) + 2328*pow4(Mt) -
        279*pow4(Mst2)*pow4(s2t))*pow8(Mst1) + 768*pow2(Mt)*(-2*pow2(Mt) +
        pow2(Mst2)*pow2(s2t))*pow8(Mst2))) + 2*pow2(Dmglst2)*(2400*pow2(Msq)*
        pow2(Mst2)*pow3(Mt)*pow4(Mst1)*(-2*s2t*pow2(Mst1)*pow3(Mst2) + 3*Mt*
        pow4(Mst1) + (9*Mt + 2*Mst2*s2t)*pow4(Mst2)) - 600*pow3(Mt)*pow4(Mst1)*
        (3*Mst2*(-3*Mt + 8*Mst2*s2t)*pow2(Mst1) - (59*Mt + 26*Mst2*s2t)*pow3(
        Mst2) + 2*s2t*pow4(Mst1))*pow5(Mst2) + pow4(Msq)*(2*pow4(Mst1)*pow4(
        Mst2)*(-40440*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 76560*Mst2*s2t*pow3(Mt) -
        5760*Mt*pow3(Mst2)*pow3(s2t) + 15712*pow4(Mt) + 4335*pow4(Mst2)*pow4(
        s2t)) + 2*pow2(Mst2)*(29760*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 31000*Mst2*
        s2t*pow3(Mt) + 15840*Mt*pow3(Mst2)*pow3(s2t) - 208652*pow4(Mt) + 1065*
        pow4(Mst2)*pow4(s2t))*pow6(Mst1) - 15*pow2(Mst1)*(-3080*pow2(Mst2)*
        pow2(Mt)*pow2(s2t) - 6144*Mst2*s2t*pow3(Mt) + 1536*Mt*pow3(Mst2)*pow3(
        s2t) + 3600*pow4(Mt) + 865*pow4(Mst2)*pow4(s2t))*pow6(Mst2) + 40*(-738*
        pow2(Mst2)*pow2(Mt)*pow2(s2t) - 8204*Mst2*s2t*pow3(Mt) + 135*Mt*pow3(
        Mst2)*pow3(s2t) + 12492*pow4(Mt) + 12*pow4(Mst2)*pow4(s2t))*pow8(Mst1)
        + 480*(-40*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 80*pow4(Mt) + 3*pow4(Mst2)*
        pow4(s2t))*pow8(Mst2)))))/(45.*pow2(Mgl)*pow4(Msq)*pow4(Mst1)*pow4(
        Mst2)*pow4(Mt));

   return result;
}

/**
 * @return returns the susy log^2 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double H6g2::calc_coef_at_as2_no_sm_logs_log2() const {

   const double result =
      (-48160*Dmglst2*Mgl*pow4(Mst1)*pow4(Mst2)*pow4(Mt) - 29264*pow2(Dmglst2)*
        pow4(Mst1)*pow4(Mst2)*pow4(Mt) + 21280*pow2(Mgl)*pow4(Mst1)*pow4(Mst2)*
        pow4(Mt) - 7200*log(pow2(Mst1)/pow2(Msq))*pow2(Mgl)*pow4(Mst1)*pow4(
        Mst2)*pow4(Mt) - 1440*pow12(Mst2)*pow2(Dmglst2)*pow4(s2t) + 46400*
        Dmglst2*Mgl*s2t*pow3(Mt)*pow4(Mst1)*pow5(Mst2) + 98880*s2t*pow2(
        Dmglst2)*pow3(Mt)*pow4(Mst1)*pow5(Mst2) - 24000*s2t*pow2(Mgl)*pow3(Mt)*
        pow4(Mst1)*pow5(Mst2) + 74880*Dmglst2*Mgl*s2t*pow3(Mst2)*pow3(Mt)*pow6(
        Mst1) - 38400*s2t*pow2(Dmglst2)*pow3(Mst2)*pow3(Mt)*pow6(Mst1) + 1920*
        s2t*pow2(Mgl)*pow3(Mst2)*pow3(Mt)*pow6(Mst1) - 14160*Dmglst2*Mgl*pow2(
        Mt)*pow2(s2t)*pow4(Mst2)*pow6(Mst1) - 13560*pow2(Dmglst2)*pow2(Mt)*
        pow2(s2t)*pow4(Mst2)*pow6(Mst1) - 9000*pow2(Mgl)*pow2(Mt)*pow2(s2t)*
        pow4(Mst2)*pow6(Mst1) - 100800*Dmglst2*Mgl*pow2(Mst2)*pow4(Mt)*pow6(
        Mst1) + 161120*pow2(Dmglst2)*pow2(Mst2)*pow4(Mt)*pow6(Mst1) - 10080*
        pow2(Mgl)*pow2(Mst2)*pow4(Mt)*pow6(Mst1) - 31200*Dmglst2*Mgl*Mt*pow3(
        s2t)*pow5(Mst2)*pow6(Mst1) - 42720*Mt*pow2(Dmglst2)*pow3(s2t)*pow5(
        Mst2)*pow6(Mst1) - 23520*Mt*pow2(Mgl)*pow3(s2t)*pow5(Mst2)*pow6(Mst1) +
        82080*Dmglst2*Mgl*pow2(Mt)*pow2(s2t)*pow4(Mst1)*pow6(Mst2) + 115440*
        pow2(Dmglst2)*pow2(Mt)*pow2(s2t)*pow4(Mst1)*pow6(Mst2) + 42960*pow2(
        Mgl)*pow2(Mt)*pow2(s2t)*pow4(Mst1)*pow6(Mst2) + 43680*Dmglst2*Mgl*pow2(
        Mst1)*pow4(Mt)*pow6(Mst2) + 65520*pow2(Dmglst2)*pow2(Mst1)*pow4(Mt)*
        pow6(Mst2) + 21840*pow2(Mgl)*pow2(Mst1)*pow4(Mt)*pow6(Mst2) + 1110*
        Dmglst2*Mgl*pow4(s2t)*pow6(Mst1)*pow6(Mst2) + 5505*pow2(Dmglst2)*pow4(
        s2t)*pow6(Mst1)*pow6(Mst2) - 5325*pow2(Mgl)*pow4(s2t)*pow6(Mst1)*pow6(
        Mst2) - 30*log(pow2(Mst1)/pow2(Mst2))*pow4(Mst1)*(pow2(Mgl)*(-4*pow4(
        Mst2)*(14*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 456*Mst2*s2t*pow3(Mt) + 146*
        Mt*pow3(Mst2)*pow3(s2t) - 106*pow4(Mt) + 33*pow4(Mst2)*pow4(s2t)) +
        pow4(Mst1)*(-768*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 1216*Mst2*s2t*pow3(Mt)
        + 400*pow4(Mt) + 41*pow4(Mst2)*pow4(s2t)) + pow2(Mst1)*pow2(Mst2)*(-
        1096*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 1280*Mst2*s2t*pow3(Mt) - 328*Mt*
        pow3(Mst2)*pow3(s2t) + 32*pow4(Mt) + 91*pow4(Mst2)*pow4(s2t))) + 2*
        Dmglst2*Mgl*(32*pow2(Mt)*(-57*Mst2*Mt*s2t + 119*pow2(Mt) - 24*pow2(
        Mst2)*pow2(s2t))*pow4(Mst1) + pow4(Mst2)*(-384*pow2(Mst2)*pow2(Mt)*
        pow2(s2t) + 912*Mst2*s2t*pow3(Mt) - 548*Mt*pow3(Mst2)*pow3(s2t) + 512*
        pow4(Mt) - 91*pow4(Mst2)*pow4(s2t)) + pow2(Mst1)*pow2(Mst2)*(-768*pow2(
        Mst2)*pow2(Mt)*pow2(s2t) - 640*Mst2*s2t*pow3(Mt) - 164*Mt*pow3(Mst2)*
        pow3(s2t) + 2016*pow4(Mt) + 91*pow4(Mst2)*pow4(s2t))) + pow2(Dmglst2)*(
        64*pow2(Mt)*(-57*Mst2*Mt*s2t + 119*pow2(Mt) - 24*pow2(Mst2)*pow2(s2t))*
        pow4(Mst1) + pow4(Mst2)*(-1152*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 1824*
        Mst2*s2t*pow3(Mt) - 1864*Mt*pow3(Mst2)*pow3(s2t) + 1536*pow4(Mt) - 273*
        pow4(Mst2)*pow4(s2t)) + pow2(Mst1)*(-2304*pow2(Mt)*pow2(s2t)*pow4(Mst2)
        + 32*pow2(Mst2)*pow4(Mt) - 328*Mt*pow3(s2t)*pow5(Mst2) + 273*pow4(s2t)*
        pow6(Mst2)))) - 46080*Dmglst2*Mgl*s2t*pow2(Mst1)*pow3(Mt)*pow7(Mst2) -
        92160*s2t*pow2(Dmglst2)*pow2(Mst1)*pow3(Mt)*pow7(Mst2) - 15360*s2t*
        pow2(Mgl)*pow2(Mst1)*pow3(Mt)*pow7(Mst2) + 19680*Dmglst2*Mgl*Mt*pow3(
        s2t)*pow4(Mst1)*pow7(Mst2) + 19680*Mt*pow2(Dmglst2)*pow3(s2t)*pow4(
        Mst1)*pow7(Mst2) + 19680*Mt*pow2(Mgl)*pow3(s2t)*pow4(Mst1)*pow7(Mst2) +
        67200*Dmglst2*Mgl*Mst2*s2t*pow3(Mt)*pow8(Mst1) + 67200*Mst2*s2t*pow2(
        Dmglst2)*pow3(Mt)*pow8(Mst1) + 1920*Mst2*s2t*pow2(Mgl)*pow3(Mt)*pow8(
        Mst1) - 112320*Dmglst2*Mgl*pow4(Mt)*pow8(Mst1) - 112320*pow2(Dmglst2)*
        pow4(Mt)*pow8(Mst1) - 12000*pow2(Mgl)*pow4(Mt)*pow8(Mst1) + 1770*
        Dmglst2*Mgl*pow4(Mst2)*pow4(s2t)*pow8(Mst1) + 1770*pow2(Dmglst2)*pow4(
        Mst2)*pow4(s2t)*pow8(Mst1) + 3585*pow2(Mgl)*pow4(Mst2)*pow4(s2t)*pow8(
        Mst1) - 29520*Dmglst2*Mgl*pow2(Mst1)*pow2(Mt)*pow2(s2t)*pow8(Mst2) -
        51960*pow2(Dmglst2)*pow2(Mst1)*pow2(Mt)*pow2(s2t)*pow8(Mst2) - 12840*
        pow2(Mgl)*pow2(Mst1)*pow2(Mt)*pow2(s2t)*pow8(Mst2) - 15360*Dmglst2*Mgl*
        pow4(Mt)*pow8(Mst2) - 38400*pow2(Dmglst2)*pow4(Mt)*pow8(Mst2) - 3840*
        pow2(Mgl)*pow4(Mt)*pow8(Mst2) - 8490*Dmglst2*Mgl*pow4(Mst1)*pow4(s2t)*
        pow8(Mst2) - 18495*pow2(Dmglst2)*pow4(Mst1)*pow4(s2t)*pow8(Mst2) - 345*
        pow2(Mgl)*pow4(Mst1)*pow4(s2t)*pow8(Mst2) + 11520*Dmglst2*Mgl*Mt*pow2(
        Mst1)*pow3(s2t)*pow9(Mst2) + 23040*Mt*pow2(Dmglst2)*pow2(Mst1)*pow3(
        s2t)*pow9(Mst2) + 3840*Mt*pow2(Mgl)*pow2(Mst1)*pow3(s2t)*pow9(Mst2) +
        7680*Dmglst2*Mgl*pow2(Mt)*pow2(s2t)*power10(Mst2) + 19200*pow2(Dmglst2)
        *pow2(Mt)*pow2(s2t)*power10(Mst2) + 1920*pow2(Mgl)*pow2(Mt)*pow2(s2t)*
        power10(Mst2) + 6570*Dmglst2*Mgl*pow2(Mst1)*pow4(s2t)*power10(Mst2) +
        13695*pow2(Dmglst2)*pow2(Mst1)*pow4(s2t)*power10(Mst2) + 2325*pow2(Mgl)
        *pow2(Mst1)*pow4(s2t)*power10(Mst2))/(45.*pow2(Mgl)*pow4(Mst1)*pow4(
        Mst2)*pow4(Mt));

   return result;
}

/**
 * @return returns the susy log^3 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double H6g2::calc_coef_at_as2_no_sm_logs_log3() const {

   const double result =
      -298.6666666666667;

   return result;
}

}        // hierarchies
}        // himalaya
