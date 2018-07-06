// ====================================================================
// This file is part of Himalaya.
//
// Himalaya is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================

#include "H9q2.hpp"
#include "Hierarchies.hpp"
#include "Constants.hpp"
#include "Utils.hpp"
#include <cmath>
#include <type_traits>

namespace himalaya{
namespace hierarchies{

/**
 * 	Constructor
 * 	@param flagMap the flagMap for the truncation of expansion variables
 * 	@param Al4p a double alpha_s/4/Pi
 * 	@param beta a double which is the mixing angle beta
 * 	@param Dmst12 a double Mst1^2 - Mst2^2
 * 	@param Dmsqst1  a double Msq1^2 - Mst1^2
 * 	@param lmMt a double log((<renormalization scale> / Mt)^2)
 * 	@param lmMgl a double log((<renormalization scale> / Mgl)^2)
 * 	@param lmMst1 a double log((<renormalization scale> / Mst1)^2)
 * 	@param Mgl a double gluino mass
 * 	@param Mt a double top/bottom quark mass
 * 	@param Mst1 a double stop 1 mass
 * 	@param Mst2 a double stop 2 mass
 * 	@param Msq a double average squark mass w/o the stop quark
 * 	@param MuSUSY a double mu parameter
 * 	@param s2t a double 2 times the sine of the stop/sbottom quark mixing angle
 * 	@param mdrFlag an int 0 for DR and 1 for MDR scheme
 * 	@param oneLoopFlag an int flag to consider the one-loop expansion terms
 * 	@param twoLoopFlag an int flag to consider the two-loop expansion terms
 * 	@param threeLoopFlag an int flag to consider the three-loop expansion terms
 */
H9q2::H9q2(const std::map<unsigned int, unsigned int>& flagMap, double Al4p, double beta, double Dmst12, double Dmsqst1,
		 double lmMt, double lmMgl, double lmMst1,
		 double Mgl, double Mt, double Mst1, double Mst2, double Msq, double MuSUSY,
		 double s2t,
		 int mdrFlag, int oneLoopFlag, int twoLoopFlag, int threeLoopFlag){
   // abbrev for tan(beta) and sin(beta)
   Tbeta = tan(beta);
   Sbeta = sin(beta);
   this -> Dmst12 = Dmst12;
   this -> Dmsqst1 = Dmsqst1;
   this -> lmMt = lmMt;
   this -> lmMgl = lmMgl;
   this -> lmMst1 = lmMst1;
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
   x = flagMap.at(ExpansionDepth::xx);
   xDmst12 = flagMap.at(ExpansionDepth::xxDmglst1);
   xDmsqst1 = flagMap.at(ExpansionDepth::xxDmsqst1);
   xMgl = flagMap.at(ExpansionDepth::xxMgl);
}

/**
 * 	@return The diagonal (1, 1) matrix element of the Higgs mass matrix as a double for the hierarchy 'H9q2'
 */
double H9q2::getS1() const {
   return -(pow2(MuSUSY)*((Al4p*(10*xMgl*pow4(Mgl)*(2*xDmst12*pow2(Mt)*(216*(1 - 3*
        lmMgl + 3*lmMst1)*twoLoopFlag*pow2(Mst1)*pow2(s2t) + Al4p*
        threeLoopFlag*(288*(262 + 335*lmMst1 - 96*lmMt - 72*lmMst1*lmMt +
        lmMgl*(-199 - 120*lmMst1 + 72*lmMt) + 120*pow2(lmMst1))*pow2(Mt) + (
        41653 + 11784*lmMst1 + 48*lmMgl*(-313 + 180*lmMst1) - 12636*pow2(lmMgl)
        + 540*pow2(lmMst1))*pow2(Mst1)*pow2(s2t)))*pow3(Dmst12) + pow2(Mst2)*(
        216*(-1 + 3*lmMgl - 3*lmMst1)*twoLoopFlag*pow2(Dmst12)*pow2(Mst1)*pow2(
        Mt)*pow2(s2t) + Al4p*threeLoopFlag*(-(pow2(Dmst12)*pow2(Mt)*(576*(262 +
        335*lmMst1 - 96*lmMt - 72*lmMst1*lmMt + lmMgl*(-199 - 120*lmMst1 + 72*
        lmMt) + 120*pow2(lmMst1))*pow2(Mt) + (41653 + 11784*lmMst1 + 48*lmMgl*(
        -313 + 180*lmMst1) - 12636*pow2(lmMgl) + 540*pow2(lmMst1))*pow2(Mst1)*
        pow2(s2t))) + 576*pow2(Mst2)*(Dmst12*(262 + 335*lmMst1 - 96*lmMt - 72*
        lmMst1*lmMt + lmMgl*(-199 - 120*lmMst1 + 72*lmMt) + 120*pow2(lmMst1)) +
        2*(59 + 85*lmMst1 - 24*lmMt - 24*lmMst1*lmMt + lmMgl*(-53 - 40*lmMst1 +
        24*lmMt) + 40*pow2(lmMst1))*pow2(Mst2))*pow4(Mt)))) + (9*Al4p*
        threeLoopFlag*(2*pow2(Mst1)*pow2(Mst2)*pow2(Mt)*(5*Dmst12*s2t*(-4*s2t*
        pow2(Msq)*pow2(Mst1)*(15*Dmsqst1*(48*Dmst12*pow2(Mgl) + pow2(Mst1)*(
        Dmst12*(3 - 6*shiftst1 - 6*shiftst2 + 4*lmMst1*(-2 + shiftst1 +
        shiftst2)) + 4*(3 - 2*lmMst1)*(shiftst1 - shiftst2)*pow2(Mst2))) +
        pow2(Msq)*(2*Dmst12*(200 + 18*lmMgl + 265*lmMst1 + 91*pow2(lmMst1))*
        pow2(Mgl) + pow2(Mst1)*(-(Dmst12*(32 + 30*(shiftst1 + shiftst2) + 2*
        lmMst1*(115 - 30*(shiftst1 + shiftst2) - 9*shiftst3) + 21*shiftst3 +
        82*pow2(lmMst1))) + 6*(1 - 2*lmMst1)*(10*shiftst1 - 10*shiftst2 +
        shiftst3)*pow2(Mst2)))) - 15*xDmsqst1*pow2(Dmsqst1)*(32*pow2(Mst2)*(20*
        Mgl*Mt*pow2(Mst1) + 8*(3 - lmMgl + lmMst1)*Mt*pow3(Mgl) - (-2 + lmMst1)
        *s2t*(shiftst1 - shiftst2)*pow4(Mst1)) + Dmst12*(192*Mgl*Mt*pow2(Mst1)
        + 192*s2t*pow2(Mgl)*pow2(Mst1) + 256*(-3 + lmMgl - lmMst1)*Mt*pow3(Mgl)
        + s2t*(33 + 16*lmMst1*(-2 + shiftst1 + shiftst2) - 32*(shiftst1 +
        shiftst2))*pow4(Mst1)))) + 32*Mgl*Mt*pow2(Msq)*(5*Dmst12*s2t*(-30*
        Dmsqst1*(Dmst12*(4*(-3 + lmMgl - lmMst1)*pow2(Mgl) + 7*pow2(Mst1)) + 2*
        ((6 - 2*lmMgl + 2*lmMst1)*pow2(Mgl) + 3*pow2(Mst1))*pow2(Mst2)) + pow2(
        Msq)*(Dmst12*(2*(-249 + 164*lmMgl - 176*lmMst1 - 20*lmMgl*lmMst1 + 27*
        pow2(lmMgl) + pow2(lmMst1))*pow2(Mgl) - (149 + 41*lmMst1)*pow2(Mst1)) +
        (4*(-11 + 154*lmMgl - 148*lmMst1 - 20*lmMgl*lmMst1 + 27*pow2(lmMgl) +
        pow2(lmMst1))*pow2(Mgl) - (215 + 118*lmMst1 + 41*pow2(lmMst1))*pow2(
        Mst1))*pow2(Mst2))) + Mgl*Mt*pow2(Msq)*(pow2(Dmst12)*(163 + 15*lmMgl -
        27*lmMst1 + 20*lmMt - 8*pow2(lmMst1)) + 20*Dmst12*(16 - 9*lmMgl + 41*
        lmMst1 - 12*lmMt + 4*pow2(lmMst1))*pow2(Mst2) - 10*(5 + 18*lmMgl - 74*
        lmMst1 + 24*lmMt - 8*pow2(lmMst1))*pow4(Mst2)))) - 16*z2*pow2(Mt)*(
        xDmst12*pow2(Mst1)*pow3(Dmst12)*(pow2(Msq)*(4*Mgl*Mt*(498*Mgl*Mt*pow2(
        Msq) - 40*s2t*pow2(Mgl)*(15*Dmsqst1 + (-94 + 45*lmMgl - 45*lmMst1)*
        pow2(Msq)) + 3*s2t*(500*Dmsqst1 - 79*pow2(Msq))*pow2(Mst1)) + pow2(
        Mst1)*(30*Dmsqst1*(30*pow2(Mgl) - (3*shiftst1 + 2*shiftst2)*pow2(Mst1))
        - 2*pow2(Msq)*(455*pow2(Mgl) + (10 + 45*shiftst1 + 30*shiftst2 + 27*
        shiftst3)*pow2(Mst1)))*pow2(s2t)) - 30*s2t*xDmsqst1*pow2(Dmsqst1)*(-
        140*Mgl*Mt*pow2(Mst1) - 30*s2t*pow2(Mgl)*pow2(Mst1) + 80*Mt*pow3(Mgl) +
        s2t*(-15 + 3*shiftst1 + 2*shiftst2)*pow4(Mst1))) + pow2(Mst2)*(4*Mt*
        pow2(Msq)*pow2(Mst1)*(Mt*pow2(Mgl)*pow2(Msq)*(11*pow2(Dmst12) - 520*
        pow2(Mst2)*(Dmst12 + pow2(Mst2))) + 5*Dmst12*s2t*(Dmst12*(-30*Mgl*(4*
        Dmsqst1 + pow2(Msq))*pow2(Mst1) + (120*Dmsqst1 + (-263 + 90*lmMgl - 90*
        lmMst1)*pow2(Msq))*pow3(Mgl)) - 2*pow2(Mst2)*(Mgl*(30*Dmsqst1 - 17*
        pow2(Msq))*pow2(Mst1) + (60*Dmsqst1 + (113 - 90*lmMgl + 90*lmMst1)*
        pow2(Msq))*pow3(Mgl)))) + Dmst12*(-5*pow2(Msq)*(30*Dmsqst1*(3*Dmst12*
        pow2(Mgl) + pow2(Mst1)*(-(Dmst12*(-3 + shiftst1 + shiftst2)) + 2*(
        shiftst1 - shiftst2)*pow2(Mst2))) - pow2(Msq)*(274*Dmst12*pow2(Mgl) +
        pow2(Mst1)*(-41*Dmst12 + 30*Dmst12*(shiftst1 + shiftst2) + 9*Dmst12*
        shiftst3 - 6*(10*shiftst1 - 10*shiftst2 + shiftst3)*pow2(Mst2))))*pow2(
        s2t)*pow4(Mst1) - 75*s2t*xDmsqst1*pow2(Dmsqst1)*pow2(Mst1)*(Dmst12*(8*
        Mgl*Mt*pow2(Mst1) + 6*s2t*pow2(Mgl)*pow2(Mst1) - 32*Mt*pow3(Mgl) + s2t*
        (9 - 2*shiftst1 - 2*shiftst2)*pow4(Mst1)) + 4*pow2(Mst2)*(10*Mgl*Mt*
        pow2(Mst1) + 8*Mt*pow3(Mgl) + s2t*(shiftst1 - shiftst2)*pow4(Mst1)))))
        + 20*xMgl*pow4(Mgl)*pow4(Msq)*(-6*xDmst12*(12*(17 - 4*lmMgl + 4*lmMst1)
        *pow2(Mt) + (49 - 17*lmMgl + 17*lmMst1)*pow2(Mst1)*pow2(s2t))*pow3(
        Dmst12) + pow2(Mst2)*(3*(49 - 17*lmMgl + 17*lmMst1)*pow2(Dmst12)*pow2(
        Mst1)*pow2(s2t) + 8*pow2(Mt)*(9*Dmst12*(17 - 4*lmMgl + 4*lmMst1)*(
        Dmst12 - pow2(Mst2)) + 2*(-47 + 12*lmMgl - 12*lmMst1)*pow4(Mst2))))) +
        3*z3*(3*pow2(Dmst12)*pow2(Mst1)*pow2(Mt)*(-32*Dmst12*Mgl*Mt*s2t*
        xDmst12*pow2(Mst1)*pow4(Msq) + 40*pow2(Mgl)*(-14*Dmst12*xDmst12*pow2(
        Mt) + 15*pow2(Mst1)*(-(Dmst12*xDmst12) + 2*pow2(Mst2))*pow2(s2t))*pow4(
        Msq) + 34560*Dmst12*Mt*s2t*xDmst12*pow3(Mgl)*pow4(Msq) - 5*pow2(s2t)*(
        35*Dmsqst1*pow2(Msq)*(Dmst12*xDmst12 + 4*pow2(Mst2)) + 35*xDmsqst1*
        pow2(Dmsqst1)*(-5*Dmst12*xDmst12 + 7*pow2(Mst2)) + 2*(11*Dmst12*xDmst12
        + 24*pow2(Mst2))*pow4(Msq))*pow4(Mst1)) + 40*Mgl*pow4(Msq)*(pow2(Mst1)*
        pow2(Mst2)*(56*Dmst12*s2t*pow2(Mst1)*pow2(Mst2) - Mgl*(648*Dmst12*Mgl*
        s2t*(Dmst12 + 2*pow2(Mst2)) + 7*Mt*(pow2(Dmst12) - 8*pow2(Mst2)*(Dmst12
        + pow2(Mst2)))))*pow3(Mt) + xMgl*pow3(Mgl)*(-6*xDmst12*pow3(Dmst12)*(
        259*pow2(Mst1)*pow2(Mt)*pow2(s2t) + 484*pow4(Mt)) + pow2(Mst2)*(777*
        pow2(Dmst12)*pow2(Mst1)*pow2(Mt)*pow2(s2t) + 968*(3*pow2(Dmst12) - 3*
        Dmst12*pow2(Mst2) - 2*pow4(Mst2))*pow4(Mt)))))))/pow4(Msq)))/pow6(Mst1)
        + pow2(Mt)*(-3*xDmst12*pow3(Dmst12)*(135*oneLoopFlag*pow2(s2t) + (2*
        threeLoopFlag*pow2(Al4p)*(16*Mgl*Mt*(12*Mgl*(20*Mgl*s2t*(-130 + lmMgl*(
        159 - 20*lmMst1) - 162*lmMst1 + 27*pow2(lmMgl) + pow2(lmMst1)) + Mt*(
        323 - 75*lmMgl + 383*lmMst1 - 100*lmMt + 32*pow2(lmMst1))) + (3600*
        Dmsqst1*(3 - lmMgl + lmMst1)*s2t*pow2(Mgl))/pow2(Msq) + s2t*(365 -
        1049*lmMst1 - 123*pow2(lmMst1) - (18000*Dmsqst1)/pow2(Msq))*pow2(Mst1))
        - (3*pow2(Mst1)*(15*Dmsqst1*(1920*pow2(Mgl) + (5 - 72*shiftst1 - 48*
        shiftst2 + 16*lmMst1*(-5 + 3*shiftst1 + 2*shiftst2))*pow2(Mst1)) + 2*
        pow2(Msq)*(10*(-425 + 18*lmMgl + 83*lmMst1 + 91*pow2(lmMst1))*pow2(Mgl)
        + (605 - 180*shiftst1 - 120*shiftst2 - 348*shiftst3 + 8*lmMst1*(-110 +
        45*shiftst1 + 30*shiftst2 + 27*shiftst3) - 820*pow2(lmMst1))*pow2(Mst1)
        ))*pow2(s2t))/pow2(Msq) - (45*s2t*xDmsqst1*pow2(Dmsqst1)*(5120*Mgl*Mt*
        pow2(Mst1) + 1920*s2t*pow2(Mgl)*pow2(Mst1) + 1280*(-3 + lmMgl - lmMst1)
        *Mt*pow3(Mgl) + s2t*(215 - 72*shiftst1 - 208*shiftst2 + 16*lmMst1*(-5 +
        3*shiftst1 + 2*shiftst2))*pow4(Mst1)))/pow4(Msq)))/pow4(Mst1) - (72*
        Al4p*s2t*twoLoopFlag*(160*(1 - lmMgl + lmMst1)*Mt*pow2(Mst1)*pow3(Mgl)
        + 2*(9 + 2*lmMst1)*Mgl*Mt*pow4(Mst1) + 5*(1 + 2*lmMst1)*s2t*pow2(Mgl)*
        pow4(Mst1) + 40*(-1 + 3*lmMgl - 3*lmMst1)*Mt*x*pow5(Mgl) - 20*lmMst1*
        s2t*pow6(Mst1)))/pow6(Mst1)) + 135*Dmst12*s2t*pow2(Mst2)*(3*Dmst12*
        oneLoopFlag*s2t - (16*Al4p*twoLoopFlag*(4*Mgl*Mt*pow2(Mst2)*((1 - 2*
        lmMgl + 2*lmMst1)*pow2(Mgl)*pow2(Mst1) + (1 - 3*lmMgl + 3*lmMst1)*x*
        pow4(Mgl) + (2 + lmMst1)*pow4(Mst1)) + Dmst12*(2*(3 - 2*lmMgl + 2*
        lmMst1)*Mt*pow2(Mst1)*pow3(Mgl) + 2*Mgl*Mt*pow4(Mst1) + (3 + 2*lmMst1)*
        s2t*pow2(Mgl)*pow4(Mst1) + 4*(-1 + 3*lmMgl - 3*lmMst1)*Mt*x*pow5(Mgl) -
        (1 + 2*lmMst1)*s2t*pow6(Mst1))))/pow6(Mst1)))))/(19440.*pow6(Mst2));
}

/**
 * 	@return The diagonal (2, 2) matrix element of the Higgs mass matrix as a double for the hierarchy 'H9q2'
 */
double H9q2::getS2() const {
   return (576*Al4p*z2*(1440*twoLoopFlag*pow2(Sbeta)*pow3(Mgl)*pow3(Mt)*pow4(Msq)*(
        xDmst12*pow2(Mst1)*(3*x*pow2(Mgl)*(Mt*(MuSUSY - Mgl*Tbeta) + s2t*Tbeta*
        pow2(Mgl)) - (Mt*MuSUSY + 4*s2t*Tbeta*x*pow2(Mgl))*pow2(Mst1))*pow3(
        Dmst12) + (Dmst12*s2t*Tbeta*pow2(Mst1) + Mt*MuSUSY*(Dmst12 + pow2(Mst2)
        ))*pow4(Mst1)*pow4(Mst2) + x*pow2(Mgl)*pow2(Mst2)*(Dmst12*pow2(Mst1)*(
        Dmst12*(-3*Mt*(MuSUSY - Mgl*Tbeta) + s2t*Tbeta*(-3*pow2(Mgl) + pow2(
        Mst1))) + (3*Mt*(MuSUSY - Mgl*Tbeta) + s2t*Tbeta*(3*pow2(Mgl) + 2*pow2(
        Mst1)))*pow2(Mst2)) + Mt*(MuSUSY - Mgl*Tbeta)*(3*pow2(Mgl) + 2*pow2(
        Mst1))*pow4(Mst2))) - 5*pow2(Mst1)*(Mt*xMgl*pow4(Mgl)*(pow2(Mst1)*pow2(
        Sbeta)*(-9*Mt*xDmst12*pow3(Dmst12)*(400*Al4p*threeLoopFlag*xDmsqst1*
        pow2(Dmsqst1)*(Mt*MuSUSY*s2t + 6*(4 + lmMst1 - lmMt)*Tbeta*pow2(Mt) -
        Tbeta*pow2(Mst1)*pow2(s2t)) + 64*Tbeta*twoLoopFlag*pow2(Mt)*pow4(Msq))
        + 144*Tbeta*twoLoopFlag*pow2(Mst2)*pow3(Mt)*pow4(Msq)*(pow2(Dmst12) +
        2*Dmst12*pow2(Mst2) + 2*pow4(Mst2))) - Al4p*threeLoopFlag*(xDmst12*
        pow2(Msq)*pow3(Dmst12)*(3*Mt*MuSUSY*s2t*pow2(Mst1)*(8*(-49 + 17*lmMgl -
        17*lmMst1)*MuSUSY*s2t*Tbeta*pow2(Msq) + 16*Mt*(75*Dmsqst1 + 2*(170 -
        27*lmMgl + 27*lmMst1)*pow2(Msq))*pow2(Sbeta)) + Tbeta*pow2(Sbeta)*(-3*
        Mt*pow2(Mst1)*(4*(300*Dmsqst1 + (-242 + 27*lmMgl - 27*lmMst1)*pow2(Msq)
        )*pow2(Mst1) + 8*(-49 + 17*lmMgl - 17*lmMst1)*pow2(Msq)*pow2(MuSUSY))*
        pow2(s2t) + (21600*Dmsqst1*(4 + lmMst1 - lmMt)*pow2(Mst1) - 16*pow2(
        Msq)*((2363 + 324*lmMgl + 564*lmMst1 - 2004*lmMt)*pow2(Mst1) + 18*(-17
        + 4*lmMgl - 4*lmMst1)*pow2(MuSUSY)))*pow3(Mt)) + 12*MuSUSY*pow2(Msq)*(
        24*(-17 + 4*lmMgl - 4*lmMst1)*MuSUSY*Tbeta*pow3(Mt) + (49 - 17*lmMgl +
        17*lmMst1)*pow2(Sbeta)*pow3(s2t)*pow4(Mst1))) + 1800*Mt*xDmsqst1*pow2(
        Dmsqst1)*pow2(Mst1)*pow2(Mst2)*pow2(Sbeta)*(2*Dmst12*Mt*(MuSUSY*s2t +
        6*(4 + lmMst1 - lmMt)*Mt*Tbeta)*pow2(Mst2) + pow2(Dmst12)*(-2*Mt*
        MuSUSY*s2t - 12*(4 + lmMst1 - lmMt)*Tbeta*pow2(Mt) + Tbeta*pow2(Mst1)*
        pow2(s2t)) + 2*(19 + 5*lmMst1 - 5*lmMt)*Tbeta*pow2(Mt)*pow4(Mst2)))) +
        Al4p*threeLoopFlag*(-(xMgl*pow2(Msq)*pow2(Mst2)*pow2(Mt)*pow4(Mgl)*(3*
        Dmst12*s2t*pow2(Mst1)*(-600*Dmsqst1*(2*Dmst12*Mt*MuSUSY - Dmst12*s2t*
        Tbeta*pow2(Mst1) - 2*Mt*MuSUSY*pow2(Mst2))*pow2(Sbeta) + 4*pow2(Msq)*(
        2*Dmst12*(-146 + 27*lmMgl - 27*lmMst1)*Mt*MuSUSY*pow2(Sbeta) + 4*(-194
        + 27*lmMgl - 27*lmMst1)*Mt*MuSUSY*pow2(Mst2)*pow2(Sbeta) + Dmst12*s2t*
        Tbeta*((-49 + 17*lmMgl - 17*lmMst1)*pow2(MuSUSY)*(-1 + pow2(Sbeta)) +
        2*(-194 + 27*lmMgl - 27*lmMst1)*pow2(Mst1)*pow2(Sbeta)))) + Tbeta*pow2(
        Mt)*(32*pow2(Msq)*pow2(MuSUSY)*(9*Dmst12*(17 - 4*lmMgl + 4*lmMst1)*(
        Dmst12 - pow2(Mst2)) + 2*(-47 + 12*lmMgl - 12*lmMst1)*pow4(Mst2)) -
        pow2(Sbeta)*(32*pow2(Msq)*pow2(MuSUSY)*(9*Dmst12*(17 - 4*lmMgl + 4*
        lmMst1)*(Dmst12 - pow2(Mst2)) + 2*(-47 + 12*lmMgl - 12*lmMst1)*pow4(
        Mst2)) + pow2(Mst1)*(7200*Dmsqst1*(4 + lmMst1 - lmMt)*(3*pow2(Dmst12) -
        3*Dmst12*pow2(Mst2) - pow4(Mst2)) - 4*pow2(Msq)*(3*(2997 + 108*lmMgl +
        798*lmMst1 - 1214*lmMt)*pow2(Dmst12) + 2*Dmst12*(-4265 + 324*lmMgl -
        1266*lmMst1 - 366*lmMt)*pow2(Mst2) + 2*(-5207 + 324*lmMgl - 1716*lmMst1
        + 84*lmMt)*pow4(Mst2))))))) + 3*xDmsqst1*xDmst12*pow2(Dmsqst1)*pow2(
        Mst1)*pow3(Dmst12)*(Mt*(10*Mgl*pow2(Mst1)*(4*s2t*pow2(Mt)*(7*Tbeta*
        pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 63*Mgl*MuSUSY*pow2(Sbeta) - 12*(13 +
        4*lmMst1 - 4*lmMt)*Tbeta*pow2(Mst1)*pow2(Sbeta)) - 6*Mt*pow2(s2t)*(-(
        Mgl*Tbeta*pow2(MuSUSY)*(-1 + pow2(Sbeta))) + 4*(MuSUSY - 3*Mgl*Tbeta)*
        pow2(Mst1)*pow2(Sbeta)) - 12*((78 + 20*lmMst1 - 20*lmMt)*MuSUSY - 3*(81
        + 22*lmMst1 - 22*lmMt)*Mgl*Tbeta)*pow2(Sbeta)*pow3(Mt) + pow2(Mst1)*(3*
        Mgl*MuSUSY + 10*Tbeta*pow2(Mst1))*pow2(Sbeta)*pow3(s2t)) + 80*pow3(Mgl)
        *(MuSUSY*pow2(Sbeta)*(-6*Mt*pow2(Mst1)*pow2(s2t) + 72*(4 + lmMst1 -
        lmMt)*pow3(Mt)) + Tbeta*(-2*s2t*pow2(Mt)*(pow2(MuSUSY)*(-1 + pow2(
        Sbeta)) + 3*(65 + 14*lmMst1 - 14*lmMt)*pow2(Mst1)*pow2(Sbeta)) + pow2(
        Sbeta)*pow3(s2t)*pow4(Mst1)))) + pow4(Mst1)*(2*Tbeta*pow2(Mt)*pow2(s2t)
        *(-((-15 + 3*shiftst1 + 2*shiftst2)*pow2(MuSUSY)*(-1 + pow2(Sbeta))) +
        10*(9 - 8*shiftst1 + 2*shiftst2)*pow2(Mst1)*pow2(Sbeta)) + 5*pow2(
        Sbeta)*(96*MuSUSY*s2t*pow3(Mt) + Mt*MuSUSY*(9 - 8*shiftst1 + 4*
        shiftst2)*pow2(Mst1)*pow3(s2t) + 16*(11 + 6*lmMst1 - 6*lmMt)*Tbeta*
        pow4(Mt) + 2*(shiftst1 - shiftst2)*Tbeta*pow4(Mst1)*pow4(s2t)))))) -
        Al4p*threeLoopFlag*pow4(Mst1)*(75*xDmsqst1*pow2(Dmsqst1)*pow2(Mst2)*(
        Tbeta*pow2(Mt)*(Dmst12*s2t*(8*Mgl*Mt*(4*pow2(Mgl) + 5*pow2(Mst1))*pow2(
        Mst2)*pow2(MuSUSY) - Dmst12*(-8*Mgl*Mt*pow2(Mst1)*pow2(MuSUSY) + 6*s2t*
        pow2(Mgl)*pow2(Mst1)*(pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 30*pow2(Mst1)*
        pow2(Sbeta)) + 32*Mt*pow2(MuSUSY)*pow3(Mgl) + 9*s2t*(pow2(MuSUSY)*(-1 +
        pow2(Sbeta)) + 8*pow2(Mst1)*pow2(Sbeta))*pow4(Mst1))) - 2*pow4(Mst1)*(
        pow2(Dmst12)*pow2(s2t)*(-((shiftst1 + shiftst2)*pow2(MuSUSY)*(-1 +
        pow2(Sbeta))) + 6*(-3*shiftst1 + shiftst2)*pow2(Mst1)*pow2(Sbeta)) + 2*
        Dmst12*pow2(Mst2)*((shiftst1 - shiftst2)*pow2(MuSUSY)*pow2(s2t)*(-1 +
        pow2(Sbeta)) + 12*(shiftst2*pow2(Mt) + (shiftst1 - shiftst2)*pow2(Mst1)
        *pow2(s2t))*pow2(Sbeta)) + 24*(shiftst1 + shiftst2)*pow2(Mt)*pow2(
        Sbeta)*pow4(Mst2)) + 4*Mt*pow2(Sbeta)*(2*Dmst12*s2t*(Mgl*(Dmst12*(4*
        pow2(Mgl) - pow2(Mst1)) - (4*pow2(Mgl) + 5*pow2(Mst1))*pow2(Mst2))*
        pow2(MuSUSY) + 12*pow2(Mst1)*(Mgl*((43 + 10*lmMst1 - 10*lmMt)*pow2(Mgl)
        + (15 + 4*lmMst1 - 4*lmMt)*pow2(Mst1))*pow2(Mst2) + Dmst12*((11 + 3*
        lmMst1 - 3*lmMt)*Mgl*pow2(Mst1) + (11 + 2*lmMst1 - 2*lmMt)*pow3(Mgl))))
        - 9*Mt*pow2(Mst1)*(pow2(Dmst12)*(4*(9 + 2*lmMst1 - 2*lmMt)*pow2(Mgl) +
        (5 + 2*lmMst1 - 2*lmMt)*pow2(Mst1)) + 2*Dmst12*((45 + 14*lmMst1 - 14*
        lmMt)*pow2(Mgl) + (3 + 2*lmMst1 - 2*lmMt)*pow2(Mst1))*pow2(Mst2) + 4*((
        13 + 4*lmMst1 - 4*lmMt)*pow2(Mgl) + (1 + lmMst1 - lmMt)*pow2(Mst1))*
        pow4(Mst2)))) + 6*MuSUSY*pow2(Sbeta)*(pow3(Mt)*(4*Dmst12*pow2(Mst2)*(2*
        (45 + 14*lmMst1 - 14*lmMt)*Mgl*Mt*pow2(Mst1) - 15*s2t*pow2(Mgl)*pow2(
        Mst1) + 48*(4 + lmMst1 - lmMt)*Mt*pow3(Mgl) + 2*s2t*(-3 + shiftst2)*
        pow4(Mst1)) + 8*(2*(13 + 4*lmMst1 - 4*lmMt)*Mgl*Mt*pow2(Mst1) + 4*(19 +
        5*lmMst1 - 5*lmMt)*Mt*pow3(Mgl) + s2t*(-shiftst1 + shiftst2)*pow4(Mst1)
        )*pow4(Mst2)) - Mt*pow2(Dmst12)*(4*Mgl*((6 + 4*lmMst1 - 4*lmMt)*Mt + 3*
        Mgl*s2t)*pow2(Mst1)*pow2(Mt) + 8*pow3(Mgl)*(-(Mt*pow2(Mst1)*pow2(s2t))
        + 24*(4 + lmMst1 - lmMt)*pow3(Mt)) + 2*Mt*s2t*(6*Mt - 5*Mgl*s2t)*pow4(
        Mst1) - (shiftst1 - shiftst2)*pow3(s2t)*pow6(Mst1)))) - xDmst12*pow2(
        Msq)*pow3(Dmst12)*(Tbeta*(4*Mgl*(498*Mgl*Mt*pow2(Msq) - 40*s2t*pow2(
        Mgl)*(15*Dmsqst1 + (-94 + 45*lmMgl - 45*lmMst1)*pow2(Msq)) + 3*s2t*(
        500*Dmsqst1 - 79*pow2(Msq))*pow2(Mst1))*pow2(MuSUSY)*pow3(Mt) + pow2(
        s2t)*(5*pow2(Mst1)*pow2(Mt)*(2*(90*Dmsqst1*pow2(Mgl) - pow2(Msq)*(91*
        pow2(Mgl) + 2*pow2(Mst1)))*pow2(MuSUSY) + 6*shiftst1*(Dmsqst1 + pow2(
        Msq))*pow2(Mst1)*(80*pow2(Mst1) + 3*pow2(MuSUSY))*pow2(Sbeta)) - 30*(
        Dmsqst1 + pow2(Msq))*pow4(Mst1)*(pow2(Mt)*(20*shiftst2*pow2(Mst1)*pow2(
        Sbeta) + pow2(MuSUSY)*(3*shiftst1 + 2*shiftst2 - 2*shiftst2*pow2(Sbeta)
        )) + 5*(shiftst1 - shiftst2)*pow2(s2t)*pow2(Sbeta)*pow4(Mst1))) + Mt*
        pow2(Sbeta)*(-20*(Mgl*(30*Dmsqst1 - 17*pow2(Msq))*pow2(Mst1) + (60*
        Dmsqst1 + (113 - 90*lmMgl + 90*lmMst1)*pow2(Msq))*pow3(Mgl))*pow3(s2t)*
        pow4(Mst1) + 10*Mt*pow2(Mst1)*pow2(s2t)*(pow2(Msq)*(2*pow2(Mst1)*pow2(
        MuSUSY) + pow2(Mgl)*(1332*pow2(Mst1) + 91*pow2(MuSUSY)) - 90*pow4(Mst1)
        ) - 90*Dmsqst1*(pow2(Mgl)*(-6*pow2(Mst1) + pow2(MuSUSY)) + 2*pow4(Mst1)
        )) - 4*Mgl*s2t*pow2(Mt)*(300*Dmsqst1*(5*pow2(Mst1)*pow2(MuSUSY) - 2*
        pow2(Mgl)*(30*(9 + 2*lmMst1 - 2*lmMt)*pow2(Mst1) + pow2(MuSUSY)) + 6*
        pow4(Mst1)) + pow2(Msq)*(-237*pow2(Mst1)*pow2(MuSUSY) + 20*pow2(Mgl)*(
        21*(56 + 11*lmMst1 - 11*lmMt)*pow2(Mst1) + 2*(94 - 45*lmMgl + 45*
        lmMst1)*pow2(MuSUSY)) + 60*(32 + 5*lmMst1 - 5*lmMt)*pow4(Mst1))) -
        pow3(Mt)*(1200*Dmsqst1*(63*(7 + 2*lmMst1 - 2*lmMt)*pow2(Mgl)*pow2(Mst1)
        - pow4(Mst1)) - 8*pow2(Msq)*(pow2(Mgl)*(5*(929 + 210*lmMst1 - 210*lmMt)
        *pow2(Mst1) - 249*pow2(MuSUSY)) + 5*(17 - 24*lmMst1 + 24*lmMt)*pow4(
        Mst1)))) - 3*shiftst3*pow2(Msq)*pow4(Mst1)*(-2*pow2(Mt)*pow2(s2t)*(9*
        pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 145*pow2(Mst1)*pow2(Sbeta)) + 5*pow2(
        Sbeta)*(24*pow4(Mt) + pow4(Mst1)*pow4(s2t)))) + Mt*MuSUSY*pow2(Sbeta)*(
        pow2(Msq)*(-80*Mgl*Mt*pow2(Mst1)*((16 - 9*lmMst1 + 9*lmMt)*pow2(Mt) +
        24*pow2(Mst1)*pow2(s2t)) + 10*Mt*(32*(380 + 147*lmMst1 - 249*lmMt)*
        pow2(Mt) - 3*(37 + 90*lmMgl - 90*lmMst1)*pow2(Mst1)*pow2(s2t))*pow3(
        Mgl) + 5*s2t*(-8*(-16 + 9*shiftst3)*pow2(Mt) + (-41 + 120*shiftst1 -
        60*shiftst2 + 21*shiftst3)*pow2(Mst1)*pow2(s2t))*pow4(Mst1) + 2*pow2(
        Mgl)*(4556*s2t*pow2(Mst1)*pow2(Mt) + 685*pow3(s2t)*pow4(Mst1))) - 150*
        Dmsqst1*(12*Mgl*Mt*pow2(Mst1)*(-16*(7 + 2*lmMst1 - 2*lmMt)*pow2(Mt) +
        pow2(Mst1)*pow2(s2t)) + 48*pow3(Mgl)*(-(Mt*pow2(Mst1)*pow2(s2t)) + 12*(
        4 + lmMst1 - lmMt)*pow3(Mt)) + 3*pow2(Mgl)*(120*s2t*pow2(Mst1)*pow2(Mt)
        + pow3(s2t)*pow4(Mst1)) + (3 - 4*shiftst1 + 2*shiftst2)*pow3(s2t)*pow6(
        Mst1)))))) + threeLoopFlag*pow2(Al4p)*(-576*pow4(Mst1)*(900*xDmsqst1*
        xDmst12*xMgl*pow2(Dmsqst1)*pow2(Mt)*(2*(29 - 10*lmMgl + 10*lmMst1)*Mt*
        MuSUSY*s2t + Tbeta*(307 + 4*lmMst1*(67 - 5*lmMt) - 10*lmMgl*(3 + 2*
        lmMst1 - 2*lmMt) - 238*lmMt + 20*pow2(lmMst1))*pow2(Mt) + 2*(-29 + 10*
        lmMgl - 10*lmMst1)*Tbeta*pow2(Mst1)*pow2(s2t))*pow2(Sbeta)*pow3(Dmst12)
        *pow4(Mgl) - z2*pow2(Msq)*pow2(Mst2)*(5*MuSUSY*pow2(Sbeta)*(12*Dmst12*
        s2t*pow2(Mst1)*pow2(Mt)*(6*Mt*(10*Dmsqst1*(Dmst12*(6*pow2(Mgl) + pow2(
        Mst1)) + (3*pow2(Mgl) + 2*pow2(Mst1))*pow2(Mst2)) + pow2(Msq)*(5*
        Dmst12*(3*pow2(Mgl) + pow2(Mst1)) + (-59*pow2(Mgl) + 10*pow2(Mst1))*
        pow2(Mst2))) - Dmst12*s2t*(Mgl*(30*Dmsqst1 - 17*pow2(Msq))*pow2(Mst1) +
        (60*Dmsqst1 + (113 - 90*lmMgl + 90*lmMst1)*pow2(Msq))*pow3(Mgl))) +
        pow4(Mst1)*(90*(Dmsqst1 + pow2(Msq))*(-8*Dmst12*s2t*shiftst2*pow2(Mst2)
        *pow3(Mt) + Mt*(-shiftst1 + shiftst2)*pow2(Dmst12)*pow2(Mst1)*pow3(s2t)
        + 8*s2t*(shiftst1 - shiftst2)*pow3(Mt)*pow4(Mst2)) + 9*shiftst3*pow2(
        Msq)*(-(Mt*pow2(Dmst12)*pow2(Mst1)*pow3(s2t)) + 8*s2t*pow3(Mt)*(pow2(
        Dmst12) - Dmst12*pow2(Mst2) + pow4(Mst2)))) + (-16*Dmst12*Mgl*(180*
        Dmsqst1*(6*(4 + lmMst1 - lmMt)*pow2(Mgl) + (7 + 2*lmMst1 - 2*lmMt)*
        pow2(Mst1)) - pow2(Msq)*((1627 + 234*lmMst1 + 174*lmMt)*pow2(Mgl) + (-
        457 - 87*lmMst1 + 87*lmMt)*pow2(Mst1)))*pow2(Mst2) + pow2(Dmst12)*(-8*
        Mgl*(180*Dmsqst1*(7 + 2*lmMst1 - 2*lmMt) + (691 + 216*lmMst1 - 216*
        lmMt)*pow2(Msq))*pow2(Mst1) + 24*(720*Dmsqst1*(4 + lmMst1 - lmMt) + (-
        1049 - 274*lmMst1 + 274*lmMt)*pow2(Msq))*pow3(Mgl)) - 32*(90*Dmsqst1*(4
        + lmMst1 - lmMt)*Mgl*(2*pow2(Mgl) + pow2(Mst1)) + Mgl*pow2(Msq)*((-1070
        - 207*lmMst1 + 3*lmMt)*pow2(Mgl) + (100 - 3*lmMst1 + 3*lmMt)*pow2(Mst1)
        ))*pow4(Mst2))*pow4(Mt)) + Tbeta*(pow2(Mt)*pow2(MuSUSY)*(4*pow2(Mgl)*
        pow2(Msq)*(11*pow2(Dmst12) - 520*pow2(Mst2)*(Dmst12 + pow2(Mst2)))*
        pow2(Mt) - 5*pow2(Dmst12)*pow2(Mst1)*(90*Dmsqst1*(pow2(Mgl) + pow2(
        Mst1)) + pow2(Msq)*(-274*pow2(Mgl) + 41*pow2(Mst1)))*pow2(s2t) + 20*
        Dmst12*Mt*s2t*(Dmst12*(-30*Mgl*(4*Dmsqst1 + pow2(Msq))*pow2(Mst1) + (
        120*Dmsqst1 + (-263 + 90*lmMgl - 90*lmMst1)*pow2(Msq))*pow3(Mgl)) - 2*
        pow2(Mst2)*(Mgl*(30*Dmsqst1 - 17*pow2(Msq))*pow2(Mst1) + (60*Dmsqst1 +
        (113 - 90*lmMgl + 90*lmMst1)*pow2(Msq))*pow3(Mgl)))) + pow2(Mt)*pow2(
        Sbeta)*(-5*Dmst12*Mt*s2t*(4*pow2(MuSUSY)*(Dmst12*(-30*Mgl*(4*Dmsqst1 +
        pow2(Msq))*pow2(Mst1) + (120*Dmsqst1 + (-263 + 90*lmMgl - 90*lmMst1)*
        pow2(Msq))*pow3(Mgl)) - 2*pow2(Mst2)*(Mgl*(30*Dmsqst1 - 17*pow2(Msq))*
        pow2(Mst1) + (60*Dmsqst1 + (113 - 90*lmMgl + 90*lmMst1)*pow2(Msq))*
        pow3(Mgl))) + pow2(Mst1)*(16*Mgl*(90*Dmsqst1*(9 + 2*lmMst1 - 2*lmMt)*(
        2*pow2(Mgl) + pow2(Mst1)) + pow2(Msq)*((-2143 - 414*lmMst1 + 6*lmMt)*
        pow2(Mgl) + (197 - 6*lmMst1 + 6*lmMt)*pow2(Mst1)))*pow2(Mst2) + Dmst12*
        (8*Mgl*(180*Dmsqst1*(7 + 2*lmMst1 - 2*lmMt) + (607 + 186*lmMst1 - 186*
        lmMt)*pow2(Msq))*pow2(Mst1) + 144*(40*Dmsqst1*(9 + 2*lmMst1 - 2*lmMt) +
        (67 + 20*lmMst1 - 20*lmMt)*pow2(Msq))*pow3(Mgl)))) + 5*pow2(Dmst12)*
        pow2(Mst1)*pow2(s2t)*(36*pow2(Mgl)*(30*Dmsqst1 - 59*pow2(Msq))*pow2(
        Mst1) + (90*Dmsqst1*(pow2(Mgl) + pow2(Mst1)) + pow2(Msq)*(-274*pow2(
        Mgl) + 41*pow2(Mst1)))*pow2(MuSUSY) + 360*(2*Dmsqst1 + pow2(Msq))*pow4(
        Mst1)) + pow2(Mt)*(-4*pow2(Mgl)*pow2(Msq)*(11*pow2(Dmst12) - 520*pow2(
        Mst2)*(Dmst12 + pow2(Mst2)))*pow2(MuSUSY) + (7200*Dmsqst1*(1 + lmMst1 -
        lmMt) + 40*(223 + 180*lmMst1 - 180*lmMt)*pow2(Msq))*pow4(Mst1)*pow4(
        Mst2) + 40*pow2(Mst1)*(pow2(Dmst12)*(pow2(Msq)*(4*(613 + 165*lmMst1 -
        165*lmMt)*pow2(Mgl) + 23*(5 + 3*lmMst1 - 3*lmMt)*pow2(Mst1)) + 45*
        Dmsqst1*(15*(7 + 2*lmMst1 - 2*lmMt)*pow2(Mgl) + (5 + 2*lmMst1 - 2*lmMt)
        *pow2(Mst1))) + 2*Dmst12*(45*Dmsqst1*(6*(7 + 2*lmMst1 - 2*lmMt)*pow2(
        Mgl) + (3 + 2*lmMst1 - 2*lmMt)*pow2(Mst1)) + pow2(Msq)*((679 + 150*
        lmMst1 - 150*lmMt)*pow2(Mgl) + 45*(2 + lmMst1 - lmMt)*pow2(Mst1)))*
        pow2(Mst2) + 4*pow2(Mgl)*(135*Dmsqst1*(4 + lmMst1 - lmMt) + (124 + 15*
        lmMst1 - 15*lmMt)*pow2(Msq))*pow4(Mst2)))) + pow4(Mst1)*(150*(Dmsqst1 +
        pow2(Msq))*pow2(Mt)*(-(pow2(Dmst12)*pow2(s2t)*((shiftst1 + shiftst2)*
        pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 6*(3*shiftst1 - shiftst2)*pow2(Mst1)*
        pow2(Sbeta))) + 2*Dmst12*pow2(Mst2)*((shiftst1 - shiftst2)*pow2(MuSUSY)
        *pow2(s2t)*(-1 + pow2(Sbeta)) + 12*(shiftst2*pow2(Mt) + (shiftst1 -
        shiftst2)*pow2(Mst1)*pow2(s2t))*pow2(Sbeta)) + 24*(shiftst1 + shiftst2)
        *pow2(Mt)*pow2(Sbeta)*pow4(Mst2)) + 15*shiftst3*pow2(Msq)*(-(Dmst12*
        pow2(Mt)*pow2(s2t)*((3*Dmst12 - 2*pow2(Mst2))*pow2(MuSUSY)*(-1 + pow2(
        Sbeta)) + 6*pow2(Mst1)*(7*Dmst12 - 4*pow2(Mst2))*pow2(Sbeta))) + 24*
        pow2(Sbeta)*(pow2(Dmst12) - Dmst12*pow2(Mst2) + pow4(Mst2))*pow4(Mt))))
        )) + 27*Mt*z3*pow2(Mst1)*(75*Dmsqst1*pow2(Mst1)*(Dmsqst1*xDmsqst1*(-4*
        Mt*pow2(Dmst12)*pow2(Mst2)*((3*Mt*(70*MuSUSY*s2t + 449*Mt*Tbeta)*pow2(
        Sbeta) + 7*Tbeta*pow2(s2t)*(7*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 30*
        pow2(Mst1)*pow2(Sbeta)))*pow4(Mst1) - 2048*Mgl*Mt*pow2(Sbeta)*(-(Mt*
        MuSUSY*pow2(Mst1)) - 3*Mgl*Mt*Tbeta*pow2(Mst1) + 2*pow2(Mgl)*(-6*Mt*
        MuSUSY + s2t*Tbeta*pow2(Mst1)) + 15*Mt*Tbeta*xMgl*pow3(Mgl) + 3*s2t*
        Tbeta*pow4(Mst1))) + xDmst12*pow3(Dmst12)*((140*Mt*Tbeta*pow2(s2t)*(
        pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 3*pow2(Mst1)*pow2(Sbeta)) + 7184*
        Tbeta*pow2(Sbeta)*pow3(Mt) + 7*MuSUSY*pow2(Sbeta)*(159*s2t*pow2(Mt) +
        28*pow2(Mst1)*pow3(s2t)))*pow4(Mst1) - 4096*Mgl*pow2(Mt)*pow2(Sbeta)*(
        10*Mt*MuSUSY*pow2(Mst1) - 33*Mgl*Mt*Tbeta*pow2(Mst1) - 4*pow2(Mgl)*(6*
        Mt*MuSUSY - 7*s2t*Tbeta*pow2(Mst1)) + 30*Mt*Tbeta*xMgl*pow3(Mgl) + 8*
        s2t*Tbeta*pow4(Mst1))) + 16*pow2(Mt)*pow2(Sbeta)*(1024*Dmst12*(6*Mt*
        MuSUSY + 5*s2t*Tbeta*pow2(Mst1))*pow3(Mgl) - 768*Dmst12*Mt*Tbeta*(7*
        pow2(Mgl)*pow2(Mst1) + 10*xMgl*pow4(Mgl)) - 3*Dmst12*(35*MuSUSY*s2t +
        214*Mt*Tbeta)*pow4(Mst1) - 4*Mt*pow2(Mst2)*(-512*Mgl*MuSUSY*pow2(Mst1)
        + 768*Tbeta*pow2(Mgl)*pow2(Mst1) - 1280*MuSUSY*pow3(Mgl) + 1600*Tbeta*
        xMgl*pow4(Mgl) + 171*Tbeta*pow4(Mst1)) + 512*Dmst12*Mgl*(7*Mt*MuSUSY*
        pow2(Mst1) + 4*s2t*Tbeta*pow4(Mst1)))*pow4(Mst2)) + pow2(Msq)*(-16*Mt*
        pow2(Dmst12)*pow2(Mst2)*(-256*Mgl*Mt*pow2(Sbeta)*(-15*Mgl*Mt*Tbeta*
        pow2(Mst1) - 8*pow2(Mgl)*(3*Mt*MuSUSY - 2*s2t*Tbeta*pow2(Mst1)) + 4*
        pow2(Mst1)*(Mt*MuSUSY + s2t*Tbeta*pow2(Mst1)) + 30*Mt*Tbeta*xMgl*pow3(
        Mgl)) + (Mt*(49*MuSUSY*s2t + 214*Mt*Tbeta)*pow2(Sbeta) + 7*Tbeta*pow2(
        s2t)*(pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 4*pow2(Mst1)*pow2(Sbeta)))*
        pow4(Mst1)) - xDmst12*pow3(Dmst12)*(8192*Mgl*pow2(Mt)*pow2(Sbeta)*(8*
        Mt*MuSUSY*pow2(Mst1) - 21*Mgl*Mt*Tbeta*pow2(Mst1) - 4*pow2(Mgl)*(3*Mt*
        MuSUSY - 5*s2t*Tbeta*pow2(Mst1)) + 15*Mt*Tbeta*xMgl*pow3(Mgl)) - 7*(4*
        Mt*Tbeta*pow2(s2t)*(pow2(MuSUSY) + (2*pow2(Mst1) - pow2(MuSUSY))*pow2(
        Sbeta)) - 24*Tbeta*pow2(Sbeta)*pow3(Mt) + MuSUSY*pow2(Sbeta)*(31*s2t*
        pow2(Mt) + 16*pow2(Mst1)*pow3(s2t)))*pow4(Mst1)) + 64*pow2(Mt)*pow2(
        Sbeta)*(-128*Mt*pow2(Mst2)*(-2*Mgl*MuSUSY*pow2(Mst1) + 3*Tbeta*pow2(
        Mgl)*pow2(Mst1) - 4*MuSUSY*pow3(Mgl) + 5*Tbeta*xMgl*pow4(Mgl) + Tbeta*
        pow4(Mst1)) - Dmst12*(-512*(3*Mt*MuSUSY + s2t*Tbeta*pow2(Mst1))*pow3(
        Mgl) + 384*Mt*Tbeta*(2*pow2(Mgl)*pow2(Mst1) + 5*xMgl*pow4(Mgl)) + (14*
        MuSUSY*s2t + 107*Mt*Tbeta)*pow4(Mst1) - 256*Mgl*(2*Mt*MuSUSY*pow2(Mst1)
        + s2t*Tbeta*pow4(Mst1))))*pow4(Mst2))) - 8*pow4(Msq)*(-80*Dmst12*pow2(
        Mt)*pow4(Mst2)*(2*pow2(Mgl)*pow2(Mst1)*(7*Mt*Tbeta*pow2(MuSUSY)*(-1 +
        pow2(Sbeta)) + 200*Mgl*Mt*MuSUSY*pow2(Sbeta) + 4*(3*MuSUSY*s2t - 199*
        Mt*Tbeta)*pow2(Mst1)*pow2(Sbeta) - 2*Mgl*s2t*Tbeta*(81*pow2(MuSUSY)*(-1
        + pow2(Sbeta)) + 380*pow2(Mst1)*pow2(Sbeta))) + 2*(xMgl*(-363*Mt*Tbeta*
        pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 2*(279*MuSUSY*s2t + 58*Mt*Tbeta)*
        pow2(Mst1)*pow2(Sbeta))*pow4(Mgl) + Mgl*(496*Mt*MuSUSY*pow2(Sbeta) +
        s2t*Tbeta*(7*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 32*pow2(Mst1)*pow2(
        Sbeta)))*pow4(Mst1)) + 3*(31*MuSUSY*s2t - 152*Mt*Tbeta)*pow2(Sbeta)*
        pow6(Mst1)) + 20*Mt*pow2(Dmst12)*pow2(Mst2)*(pow2(Mgl)*pow2(Mst1)*(-8*
        Mgl*(MuSUSY*(2192*pow2(Mt) + 243*pow2(Mst1)*pow2(s2t))*pow2(Sbeta) + 3*
        Mt*s2t*Tbeta*(-27*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 320*pow2(Mst1)*
        pow2(Sbeta))) + Tbeta*(-6*pow2(Mst1)*pow2(s2t)*(15*pow2(MuSUSY)*(-1 +
        pow2(Sbeta)) + 8*pow2(Mst1)*pow2(Sbeta)) + pow2(Mt)*(7*pow2(MuSUSY)*(-1
        + pow2(Sbeta)) + 14080*pow2(Mst1)*pow2(Sbeta)))) - 4*Mgl*(928*Mt*s2t*
        Tbeta*pow2(Mst1) + 1136*MuSUSY*pow2(Mt) - 21*MuSUSY*pow2(Mst1)*pow2(
        s2t))*pow2(Sbeta)*pow4(Mst1) - xMgl*pow4(Mgl)*(pow2(Mst1)*(777*Tbeta*
        pow2(MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta)) - 8*Mt*(-279*MuSUSY*s2t +
        4442*Mt*Tbeta)*pow2(Sbeta)) + Tbeta*(2904*pow2(Mt)*pow2(MuSUSY)*(-1 +
        pow2(Sbeta)) + 2232*pow2(s2t)*pow2(Sbeta)*pow4(Mst1))) + (Mt*(396*
        MuSUSY*s2t + 1141*Mt*Tbeta)*pow2(Sbeta) - 6*Tbeta*pow2(s2t)*(-3*pow2(
        MuSUSY)*(-1 + pow2(Sbeta)) + 31*pow2(Mst1)*pow2(Sbeta)))*pow6(Mst1)) +
        xDmst12*pow3(Dmst12)*(-20*xMgl*pow4(Mgl)*(-2*Mt*pow2(Mst1)*(777*Tbeta*
        pow2(MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta)) + 16*Mt*(279*MuSUSY*s2t -
        2192*Mt*Tbeta)*pow2(Sbeta)) - 2904*Tbeta*pow2(MuSUSY)*(-1 + pow2(Sbeta)
        )*pow3(Mt) - 3*(259*MuSUSY*s2t + 372*Mt*Tbeta)*pow2(s2t)*pow2(Sbeta)*
        pow4(Mst1)) + 16*Mgl*pow4(Mst1)*(-(s2t*Tbeta*pow2(Mt)*(-3*pow2(MuSUSY)*
        (-1 + pow2(Sbeta)) + 1040*pow2(Mst1)*pow2(Sbeta))) + 5*MuSUSY*pow2(
        Sbeta)*(-21*Mt*pow2(Mst1)*pow2(s2t) + 88*pow3(Mt)) + 35*Tbeta*pow2(
        Sbeta)*pow3(s2t)*pow4(Mst1)) + 20*pow2(Mgl)*pow2(Mst1)*(-(Mt*pow2(Mst1)
        *(-45*Tbeta*pow2(MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta)) + 14*Mt*(MuSUSY*
        s2t - 320*Mt*Tbeta)*pow2(Sbeta))) + 42*Tbeta*pow2(MuSUSY)*(-1 + pow2(
        Sbeta))*pow3(Mt) + 6*(15*MuSUSY*s2t + 8*Mt*Tbeta)*pow2(s2t)*pow2(Sbeta)
        *pow4(Mst1) + 4*Mgl*(-8*s2t*Tbeta*pow2(Mt)*(81*pow2(MuSUSY)*(-1 + pow2(
        Sbeta)) + 308*pow2(Mst1)*pow2(Sbeta)) + 3*MuSUSY*pow2(Sbeta)*(81*Mt*
        pow2(Mst1)*pow2(s2t) + 3056*pow3(Mt)) - 162*Tbeta*pow2(Sbeta)*pow3(s2t)
        *pow4(Mst1))) - 15*(-(Mt*Tbeta*pow2(s2t)*(11*pow2(MuSUSY)*(-1 + pow2(
        Sbeta)) + 512*pow2(Mst1)*pow2(Sbeta))) + 500*Tbeta*pow2(Sbeta)*pow3(Mt)
        + MuSUSY*pow2(Sbeta)*(85*s2t*pow2(Mt) + 24*pow2(Mst1)*pow3(s2t)))*pow6(
        Mst1)) + 160*pow3(Mt)*(8*MuSUSY*pow2(Sbeta)*(95*pow2(Mst1)*pow3(Mgl) -
        4*Mgl*pow4(Mst1)) + Tbeta*(2*xMgl*(121*pow2(MuSUSY)*(-1 + pow2(Sbeta))
        - 658*pow2(Mst1)*pow2(Sbeta))*pow4(Mgl) + pow2(Mgl)*(-7*pow2(Mst1)*
        pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 152*pow2(Sbeta)*pow4(Mst1)) + 2*(395
        - 12*lmMst1 + 12*lmMt)*pow2(Sbeta)*pow6(Mst1)))*pow6(Mst2)))) + pow2(
        Mst2)*(72*threeLoopFlag*pow2(Al4p)*pow4(Mst1)*(-4*pow2(Msq)*(Tbeta*
        pow2(Mt)*(-(pow2(MuSUSY)*(-40*Dmst12*Mgl*Mt*(-4*Mgl*Mt*(-16 + 9*lmMgl -
        41*lmMst1 + 12*lmMt - 4*pow2(lmMst1))*pow2(Msq) + 4*s2t*pow2(Mgl)*(30*
        Dmsqst1*(-3 + lmMgl - lmMst1) + (-11 + 154*lmMgl - 148*lmMst1 - 20*
        lmMgl*lmMst1 + 27*pow2(lmMgl) + pow2(lmMst1))*pow2(Msq)) - s2t*(180*
        Dmsqst1 + (215 + 118*lmMst1 + 41*pow2(lmMst1))*pow2(Msq))*pow2(Mst1))*
        pow2(Mst2) - pow2(Dmst12)*(-2*pow2(Mgl)*(1800*Dmsqst1*pow2(Mst1)*pow2(
        s2t) + pow2(Msq)*(-4*(163 + 15*lmMgl - 27*lmMst1 + 20*lmMt - 8*pow2(
        lmMst1))*pow2(Mt) + 5*(200 + 18*lmMgl + 265*lmMst1 + 91*pow2(lmMst1))*
        pow2(Mst1)*pow2(s2t))) + Mt*s2t*(-40*Mgl*(210*Dmsqst1 + (149 + 41*
        lmMst1)*pow2(Msq))*pow2(Mst1) + 80*(60*Dmsqst1*(3 - lmMgl + lmMst1) + (
        -249 - 176*lmMst1 - 4*lmMgl*(-41 + 5*lmMst1) + 27*pow2(lmMgl) + pow2(
        lmMst1))*pow2(Msq))*pow3(Mgl)) - 5*(15*Dmsqst1*(3 - 8*lmMst1) - 2*(16 +
        115*lmMst1 + 41*pow2(lmMst1))*pow2(Msq))*pow2(s2t)*pow4(Mst1)) + 80*(5
        + 18*lmMgl - 74*lmMst1 + 24*lmMt - 8*pow2(lmMst1))*pow2(Mgl)*pow2(Msq)*
        pow2(Mt)*pow4(Mst2))) + pow4(Mst1)*(15*shiftst3*pow2(Msq)*(-2*Dmst12*
        pow2(Mst2)*((-1 + 2*lmMst1)*pow2(MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta)) +
        12*((3 - 2*lmMst1)*pow2(Mt) + (-1 + 2*lmMst1)*pow2(Mst1)*pow2(s2t))*
        pow2(Sbeta)) + pow2(Dmst12)*((-7 + 6*lmMst1)*pow2(MuSUSY)*pow2(s2t)*(-1
        + pow2(Sbeta)) + 6*(-8*(-2 + lmMst1)*pow2(Mt) + (-15 + 14*lmMst1)*pow2(
        Mst1)*pow2(s2t))*pow2(Sbeta)) + 24*(1 - 2*lmMst1)*pow2(Mt)*pow2(Sbeta)*
        pow4(Mst2)) + 150*(Dmsqst1*(3 - 2*lmMst1) + (1 - 2*lmMst1)*pow2(Msq))*(
        pow2(Dmst12)*pow2(s2t)*(-((shiftst1 + shiftst2)*pow2(MuSUSY)*(-1 +
        pow2(Sbeta))) + 6*(-3*shiftst1 + shiftst2)*pow2(Mst1)*pow2(Sbeta)) + 2*
        Dmst12*pow2(Mst2)*((shiftst1 - shiftst2)*pow2(MuSUSY)*pow2(s2t)*(-1 +
        pow2(Sbeta)) + 12*(shiftst2*pow2(Mt) + (shiftst1 - shiftst2)*pow2(Mst1)
        *pow2(s2t))*pow2(Sbeta)) + 24*(shiftst1 + shiftst2)*pow2(Mt)*pow2(
        Sbeta)*pow4(Mst2))) + pow2(Sbeta)*(5*pow2(Dmst12)*pow2(Mst1)*pow2(s2t)*
        (15*Dmsqst1*((3 - 8*lmMst1)*pow2(Mst1)*pow2(MuSUSY) + 48*pow2(Mgl)*(9*
        pow2(Mst1) + pow2(MuSUSY)) + (60 - 48*lmMst1)*pow4(Mst1)) + 2*pow2(Msq)
        *(6*(525 - 18*lmMgl + 504*lmMst1 - 48*lmMt + 107*pow2(lmMst1))*pow2(
        Mgl)*pow2(Mst1) + (200 + 18*lmMgl + 265*lmMst1 + 91*pow2(lmMst1))*pow2(
        Mgl)*pow2(MuSUSY) - (16 + 115*lmMst1 + 41*pow2(lmMst1))*pow2(Mst1)*
        pow2(MuSUSY) - 12*(-6 + 79*lmMst1 + 41*pow2(lmMst1))*pow4(Mst1))) - 40*
        Dmst12*Mgl*Mt*s2t*(30*Dmsqst1*(Dmst12*(-7*pow2(Mst1)*pow2(MuSUSY) + 4*
        pow2(Mgl)*(-3*(-24 + lmMgl - 15*lmMst1 + 14*lmMt)*pow2(Mst1) + (3 -
        lmMgl + lmMst1)*pow2(MuSUSY)) + 6*(5 + 6*lmMst1 - 6*lmMt)*pow4(Mst1)) +
        2*pow2(Mst2)*(-3*pow2(Mst1)*pow2(MuSUSY) - 2*pow2(Mgl)*((-57 - 39*
        lmMst1 + lmMgl*(9 + 6*lmMst1 - 6*lmMt) + 30*lmMt + 6*lmMst1*lmMt - 6*
        pow2(lmMst1))*pow2(Mst1) + (3 - lmMgl + lmMst1)*pow2(MuSUSY)) + 24*(1 +
        lmMst1 - lmMt)*pow4(Mst1))) + pow2(Msq)*(pow2(Mst2)*(-((215 + 118*
        lmMst1 + 41*pow2(lmMst1))*pow2(Mst1)*pow2(MuSUSY)) + 2*pow2(Mgl)*(2*(-
        11 + lmMgl*(154 - 20*lmMst1) - 148*lmMst1 + 27*pow2(lmMgl) + pow2(
        lmMst1))*pow2(MuSUSY) + pow2(Mst1)*(2086 - 996*lmMt - 81*(9 + 4*lmMst1
        - 6*lmMt)*pow2(lmMgl) + (663 + 366*lmMt)*pow2(lmMst1) - 144*pow2(lmMt)
        - 6*lmMst1*(-521 + 74*lmMt + 24*pow2(lmMt)) + 6*lmMgl*(-238 + lmMst1*(
        73 - 142*lmMt) + 52*lmMt + 91*pow2(lmMst1) + 24*pow2(lmMt)) - 54*pow3(
        lmMgl) - 168*pow3(lmMst1))) - 4*(-407 + 276*lmMt + 3*lmMst1*(-65 + 2*
        lmMt) + 45*pow2(lmMst1) + 72*pow2(lmMt))*pow4(Mst1)) + Dmst12*(-((149 +
        41*lmMst1)*pow2(Mst1)*pow2(MuSUSY)) + 2*pow2(Mgl)*(3*(1091 + 611*lmMst1
        + lmMgl*(-17 + 122*lmMst1 - 82*lmMt) + (-374 + 62*lmMst1)*lmMt - 54*
        pow2(lmMgl) - 24*(pow2(lmMst1) + pow2(lmMt)))*pow2(Mst1) + (-249 - 176*
        lmMst1 - 4*lmMgl*(-41 + 5*lmMst1) + 27*pow2(lmMgl) + pow2(lmMst1))*
        pow2(MuSUSY)) + 2*(148 - 90*lmMt + lmMst1*(75 + 6*lmMt) + 45*pow2(
        lmMst1) + 72*pow2(lmMt))*pow4(Mst1)))) + 2*pow2(Mt)*(5*Dmst12*pow2(
        Mst1)*(Dmst12*(45*Dmsqst1*(12*(19 + 33*lmMst1 - 33*lmMt)*pow2(Mgl) + (-
        3 + 28*lmMst1 - 28*lmMt)*pow2(Mst1)) + pow2(Msq)*((5308 - 7972*lmMt +
        4*lmMst1*(1977 + 64*lmMt) - 128*pow2(lmMst1))*pow2(Mgl) - (131 + 1500*
        lmMt - 60*lmMst1*(29 + 2*lmMt) + 306*(pow2(lmMst1) + pow2(lmMt)))*pow2(
        Mst1))) + (90*Dmsqst1*(24*(3 + 7*lmMst1 - 7*lmMt)*pow2(Mgl) + (-19 +
        32*lmMst1 - 24*lmMt)*pow2(Mst1)) + 4*pow2(Msq)*((3049 - 1908*lmMt - 54*
        lmMgl*(1 + 12*lmMt) + 6*lmMst1*(215 + 49*lmMt) + 324*pow2(lmMgl) - 195*
        pow2(lmMst1) - 144*pow2(lmMt))*pow2(Mgl) + (-196 + lmMst1*(288 - 60*
        lmMt) + 54*lmMt + 153*pow2(lmMst1) + 153*pow2(lmMt))*pow2(Mst1)))*pow2(
        Mst2)) + (40*pow2(Mgl)*(540*Dmsqst1*(1 + 2*lmMst1 - 2*lmMt) + (1909 +
        600*lmMst1 + 6*(-187 + 49*lmMst1)*lmMt - 54*lmMgl*(1 + 12*lmMt) + 324*
        pow2(lmMgl) - 195*pow2(lmMst1) - 144*pow2(lmMt))*pow2(Msq))*pow2(Mst1)
        + (-200*Dmsqst1*(92 + 66*lmMt - 6*lmMst1*(17 + 3*lmMt) + 9*(pow2(
        lmMst1) + pow2(lmMt))) + 30*pow2(Msq)*(-1375 + 64*lmMt - 24*(2 + 3*
        lmMt)*pow2(lmMgl) - 4*(37 + 25*lmMt)*pow2(lmMst1) - 120*pow2(lmMt) + 8*
        lmMgl*(26 + 12*lmMt + 9*pow2(lmMt)) + 8*lmMst1*(44 + 64*lmMt + 33*pow2(
        lmMt)) + 24*pow3(lmMgl) + 88*pow3(lmMst1) - 276*pow3(lmMt)))*pow4(Mst1)
        )*pow4(Mst2) - 4*pow2(Mgl)*pow2(Msq)*pow2(MuSUSY)*(pow2(Dmst12)*(163 +
        15*lmMgl - 27*lmMst1 + 20*lmMt - 8*pow2(lmMst1)) + 20*Dmst12*(16 - 9*
        lmMgl + 41*lmMst1 - 12*lmMt + 4*pow2(lmMst1))*pow2(Mst2) - 10*(5 + 18*
        lmMgl - 74*lmMst1 + 24*lmMt - 8*pow2(lmMst1))*pow4(Mst2))))) + 5*
        MuSUSY*pow2(Sbeta)*(3*Dmst12*s2t*pow2(Mst1)*(-4*Dmst12*Mgl*s2t*(60*
        Dmsqst1*(2*(3 - lmMgl + lmMst1)*pow2(Mgl) + 3*pow2(Mst1)) + pow2(Msq)*(
        -4*(-11 + lmMgl*(154 - 20*lmMst1) - 148*lmMst1 + 27*pow2(lmMgl) + pow2(
        lmMst1))*pow2(Mgl) + (215 + 118*lmMst1 + 41*pow2(lmMst1))*pow2(Mst1)))
        + Mt*(15*Dmsqst1*(Dmst12*(336*pow2(Mgl) + 31*pow2(Mst1)) + 8*(36*pow2(
        Mgl) + (5 - 4*lmMst1)*pow2(Mst1))*pow2(Mst2)) + 4*pow2(Msq)*(Dmst12*(2*
        (522 + 107*lmMst1)*pow2(Mgl) + (-95 - 66*lmMst1 + 82*pow2(lmMst1))*
        pow2(Mst1)) + 2*((501 - 18*lmMgl + 504*lmMst1 - 48*lmMt + 107*pow2(
        lmMst1))*pow2(Mgl) - 2*(-6 + 79*lmMst1 + 41*pow2(lmMst1))*pow2(Mst1))*
        pow2(Mst2))))*pow2(Mt) + pow4(Mst1)*(-9*Mt*s2t*shiftst3*pow2(Msq)*(-8*
        Dmst12*(-3 + 2*lmMst1)*pow2(Mst2)*pow2(Mt) + pow2(Dmst12)*(16*(-2 +
        lmMst1)*pow2(Mt) + (1 - 2*lmMst1)*pow2(Mst1)*pow2(s2t)) + 8*(-1 + 2*
        lmMst1)*pow2(Mt)*pow4(Mst2)) + 90*Mt*s2t*(Dmsqst1*(3 - 2*lmMst1) + (1 -
        2*lmMst1)*pow2(Msq))*(-8*Dmst12*shiftst2*pow2(Mst2)*pow2(Mt) + (-
        shiftst1 + shiftst2)*pow2(Dmst12)*pow2(Mst1)*pow2(s2t) + 8*(shiftst1 -
        shiftst2)*pow2(Mt)*pow4(Mst2))) + 16*Mgl*(30*Dmsqst1*(pow2(Dmst12)*(6*(
        31 + 27*lmMst1 - 24*lmMt - 2*lmMst1*lmMt + lmMgl*(-3 - 2*lmMst1 + 2*
        lmMt) + 2*pow2(lmMst1))*pow2(Mgl) + (-13 - 19*lmMst1 + 19*lmMt)*pow2(
        Mst1)) - 6*Dmst12*((31 + 27*lmMst1 - 24*lmMt - 2*lmMst1*lmMt + lmMgl*(-
        3 - 2*lmMst1 + 2*lmMt) + 2*pow2(lmMst1))*pow2(Mgl) + (3 + 7*lmMst1 - 7*
        lmMt)*pow2(Mst1))*pow2(Mst2) - 12*((7 + 6*lmMst1 - 5*lmMt - lmMst1*lmMt
        + lmMgl*(-1 - lmMst1 + lmMt) + pow2(lmMst1))*pow2(Mgl) + (1 + 2*lmMst1
        - 2*lmMt)*pow2(Mst1))*pow4(Mst2)) - pow2(Msq)*(pow2(Dmst12)*(-((-3749 -
        2374*lmMst1 + 1770*lmMt - 117*lmMst1*lmMt + lmMgl*(96 - 237*lmMst1 +
        157*lmMt) + 108*pow2(lmMgl) + 41*pow2(lmMst1) + 48*pow2(lmMt))*pow2(
        Mgl)) + (lmMst1*(297 + 2*lmMt) + 15*pow2(lmMst1) + 3*(47 - 79*lmMt + 8*
        pow2(lmMt)))*pow2(Mst1)) + Dmst12*pow2(Mst2)*(-((-605 + 6*lmMst1*(-77 +
        lmMt) + 714*lmMt + 45*pow2(lmMst1) + 72*pow2(lmMt))*pow2(Mst1)) + 2*
        pow2(Mgl)*(2096 - 1167*lmMt - 81*(4 + 2*lmMst1 - 3*lmMt)*pow2(lmMgl) +
        3*(129 + 61*lmMt)*pow2(lmMst1) - 108*pow2(lmMt) - 3*lmMst1*(-793 + 67*
        lmMt + 24*pow2(lmMt)) + 3*lmMgl*(-229 + lmMst1*(63 - 142*lmMt) + 35*
        lmMt + 91*pow2(lmMst1) + 24*pow2(lmMt)) - 27*pow3(lmMgl) - 84*pow3(
        lmMst1))) + 2*(-((-233 + 6*lmMst1*(-32 + lmMt) + 348*lmMt + 45*pow2(
        lmMst1) + 72*pow2(lmMt))*pow2(Mst1)) + pow2(Mgl)*(758 - 570*lmMt - 81*(
        3 + 2*lmMst1 - 3*lmMt)*pow2(lmMgl) + 3*(141 + 61*lmMt)*pow2(lmMst1) -
        72*pow2(lmMt) - 6*lmMst1*(-236 + 49*lmMt + 12*pow2(lmMt)) + 3*lmMgl*(-
        200 + 2*lmMst1 + 76*lmMt - 142*lmMst1*lmMt + 91*pow2(lmMst1) + 24*pow2(
        lmMt)) - 27*pow3(lmMgl) - 84*pow3(lmMst1)))*pow4(Mst2)))*pow4(Mt))) -
        25*xDmsqst1*pow2(Dmsqst1)*(Tbeta*pow2(Mt)*(3*Dmst12*s2t*(128*Mgl*Mt*(2*
        (-3 + lmMgl - lmMst1)*pow2(Mgl) - 5*pow2(Mst1))*pow2(Mst2)*pow2(MuSUSY)
        + Dmst12*(-192*Mgl*Mt*pow2(Mst1)*pow2(MuSUSY) + 192*s2t*pow2(Mgl)*pow2(
        Mst1)*(pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 15*pow2(Mst1)*pow2(Sbeta)) -
        256*(-3 + lmMgl - lmMst1)*Mt*pow2(MuSUSY)*pow3(Mgl) + s2t*(-((-33 + 32*
        lmMst1)*pow2(MuSUSY)*(-1 + pow2(Sbeta))) + 6*(87 - 32*lmMst1)*pow2(
        Mst1)*pow2(Sbeta))*pow4(Mst1))) + 48*(2 - lmMst1)*pow4(Mst1)*(pow2(
        Dmst12)*pow2(s2t)*(-((shiftst1 + shiftst2)*pow2(MuSUSY)*(-1 + pow2(
        Sbeta))) + 6*(-3*shiftst1 + shiftst2)*pow2(Mst1)*pow2(Sbeta)) + 2*
        Dmst12*pow2(Mst2)*((shiftst1 - shiftst2)*pow2(MuSUSY)*pow2(s2t)*(-1 +
        pow2(Sbeta)) + 12*(shiftst2*pow2(Mt) + (shiftst1 - shiftst2)*pow2(Mst1)
        *pow2(s2t))*pow2(Sbeta)) + 24*(shiftst1 + shiftst2)*pow2(Mt)*pow2(
        Sbeta)*pow4(Mst2)) + Mt*pow2(Sbeta)*(-96*Dmst12*Mgl*s2t*(Dmst12*(-6*
        pow2(Mst1)*pow2(MuSUSY) + 4*pow2(Mgl)*((96 + 51*lmMst1 - 3*lmMgl*(5 +
        2*lmMst1 - 2*lmMt) - 36*lmMt - 6*lmMst1*lmMt + 6*pow2(lmMst1))*pow2(
        Mst1) + 2*(3 - lmMgl + lmMst1)*pow2(MuSUSY)) + 3*(33 + 38*lmMst1 - 38*
        lmMt)*pow4(Mst1)) + 4*pow2(Mst2)*(-5*pow2(Mst1)*pow2(MuSUSY) + pow2(
        Mgl)*(3*(18 - lmMgl + lmMst1)*(3 + 2*lmMst1 - 2*lmMt)*pow2(Mst1) + 2*(-
        3 + lmMgl - lmMst1)*pow2(MuSUSY)) + 3*(9 + 14*lmMst1 - 14*lmMt)*pow4(
        Mst1))) + Mt*pow2(Mst1)*(pow2(Dmst12)*(1728*(10 + 7*lmMst1 - 7*lmMt)*
        pow2(Mgl) + (-425 + 2928*lmMst1 - 2928*lmMt)*pow2(Mst1)) + 72*Dmst12*(
        12*(11 + 47*lmMst1 - 47*lmMt)*pow2(Mgl) + (-73 + 92*lmMst1 - 76*lmMt)*
        pow2(Mst1))*pow2(Mst2) + 16*(216*(1 + 7*lmMst1 - 7*lmMt)*pow2(Mgl) - (
        559 + 384*lmMt - 12*lmMst1*(44 + 3*lmMt) + 18*(pow2(lmMst1) + pow2(
        lmMt)))*pow2(Mst1))*pow4(Mst2)))) + 6*MuSUSY*pow2(Sbeta)*(pow3(Mt)*(-6*
        Dmst12*pow2(Mst2)*(16*(11 + 47*lmMst1 - 47*lmMt)*Mgl*Mt*pow2(Mst1) -
        480*s2t*pow2(Mgl)*pow2(Mst1) + 64*Mt*(31 + 27*lmMst1 - 24*lmMt - 2*
        lmMst1*lmMt + lmMgl*(-3 - 2*lmMst1 + 2*lmMt) + 2*pow2(lmMst1))*pow3(
        Mgl) + s2t*(-87 - 32*lmMst1*(-1 + shiftst2) + 64*shiftst2)*pow4(Mst1))
        - 192*(2*Mt*((1 + 7*lmMst1 - 7*lmMt)*Mgl*pow2(Mst1) + (18 - lmMgl +
        lmMst1)*(1 + lmMst1 - lmMt)*pow3(Mgl)) + (-2 + lmMst1)*s2t*(shiftst1 -
        shiftst2)*pow4(Mst1))*pow4(Mst2)) + Mt*pow2(Dmst12)*(864*s2t*pow2(Mgl)*
        pow2(Mst1)*pow2(Mt) - 192*pow3(Mgl)*((3 - lmMgl + lmMst1)*Mt*pow2(Mst1)
        *pow2(s2t) + (-62 - 54*lmMst1 + lmMgl*(6 + 4*lmMst1 - 4*lmMt) + 48*lmMt
        + 4*lmMst1*lmMt - 4*pow2(lmMst1))*pow3(Mt)) + 321*s2t*pow2(Mt)*pow4(
        Mst1) - 32*Mgl*((29 - 19*lmMst1 + 19*lmMt)*pow2(Mst1)*pow3(Mt) + 15*Mt*
        pow2(s2t)*pow4(Mst1)) + 24*(-2 + lmMst1)*(shiftst1 - shiftst2)*pow3(
        s2t)*pow6(Mst1))))) + 540*pow2(Mt)*pow4(Msq)*(16*Al4p*twoLoopFlag*(
        MuSUSY*pow4(Mst1)*(Dmst12*MuSUSY*s2t*Tbeta*(4*Mt*pow2(Mst2)*((2 +
        lmMst1)*Mgl*pow2(Mst1) + (1 - 2*lmMgl + 2*lmMst1)*pow3(Mgl)) + Dmst12*
        s2t*((3 + 2*lmMst1)*pow2(Mgl)*pow2(Mst1) - (1 + 2*lmMst1)*pow4(Mst1)))
        + 2*pow2(Sbeta)*(12*Dmst12*Mt*pow2(Mst2)*(-((1 + lmMt)*Mgl*Mt*pow2(
        Mst1)) - (2 + lmMst1)*s2t*pow2(Mgl)*pow2(Mst1) + Mt*(2 + 6*lmMst1 - 3*
        lmMt - 2*lmMst1*lmMt + lmMgl*(-3 - 2*lmMst1 + 2*lmMt) + 2*pow2(lmMst1))
        *pow3(Mgl) + (1 + 2*lmMst1)*s2t*pow4(Mst1)) - pow2(Dmst12)*(6*Mt*s2t*
        pow2(Mgl)*pow2(Mst1) + (8*(-3 + lmMgl - 2*lmMst1 + lmMt)*pow2(Mt) + 3*(
        -1 + 2*lmMgl - 2*lmMst1)*pow2(Mst1)*pow2(s2t))*pow3(Mgl) + 6*(-1 + 2*
        lmMst1)*Mt*s2t*pow4(Mst1) - Mgl*(4*(1 + lmMt)*pow2(Mst1)*pow2(Mt) + 3*(
        2 + lmMst1)*pow2(s2t)*pow4(Mst1))) + 24*Mgl*((-(lmMst1*(-2 + lmMt)) -
        lmMt + lmMgl*(-1 - lmMst1 + lmMt) + pow2(lmMst1))*pow2(Mgl) - (1 +
        lmMt)*pow2(Mst1))*pow2(Mt)*pow4(Mst2))) + Tbeta*pow4(Mst1)*(s2t*pow2(
        Dmst12)*(2*Mt*pow2(MuSUSY)*(Mgl*pow2(Mst1) + (3 - 2*lmMgl + 2*lmMst1)*
        pow3(Mgl)) + s2t*pow2(Sbeta)*((1 + 2*lmMst1)*(12*pow2(Mst1) + pow2(
        MuSUSY))*pow4(Mst1) - pow2(Mgl)*((3 + 2*lmMst1)*pow2(Mst1)*pow2(MuSUSY)
        + 12*(2 + lmMst1)*pow4(Mst1)))) + 2*Mt*pow2(Sbeta)*(6*(lmMst1 + lmMt)*
        Mt*pow2(Dmst12)*pow4(Mst1) + Mgl*s2t*pow2(Dmst12)*(-(pow2(Mst1)*pow2(
        MuSUSY)) + pow2(Mgl)*(6*(5 + 4*lmMst1 - 2*(lmMgl + lmMt))*pow2(Mst1) +
        (-3 + 2*lmMgl - 2*lmMst1)*pow2(MuSUSY)) + 6*(1 + 2*lmMt)*pow4(Mst1)) -
        2*Dmst12*Mgl*s2t*pow2(Mst2)*((2 + lmMst1)*pow2(Mst1)*pow2(MuSUSY) +
        pow2(Mgl)*(6*(-1 - 5*lmMst1 + lmMgl*(3 + 2*lmMst1 - 2*lmMt) + 2*lmMt +
        2*lmMst1*lmMt - 2*pow2(lmMst1))*pow2(Mst1) + (1 - 2*lmMgl + 2*lmMst1)*
        pow2(MuSUSY)) + 6*(1 + 2*lmMt)*pow4(Mst1)) + 12*Mt*pow2(Mst1)*(Dmst12*(
        (1 + lmMt)*pow2(Mgl) - (1 + lmMst1 + lmMt)*pow2(Mst1))*pow2(Mst2) + (2*
        (1 + lmMt)*pow2(Mgl) - (2 + 2*lmMst1*(1 + lmMt) + pow2(lmMst1) - 3*
        pow2(lmMt))*pow2(Mst1))*pow4(Mst2)))) - x*(3*MuSUSY*pow2(Sbeta)*(2*
        Dmst12*(-1 + 3*lmMgl - 3*lmMst1)*s2t*pow2(Mst1)*(4*Dmst12*Mgl*Mt +
        Dmst12*s2t*pow2(Mst1) - 4*Mgl*Mt*pow2(Mst2)) + 8*pow2(Mt)*(Dmst12*(
        lmMst1*(22 - 6*lmMt) - 5*lmMt + lmMgl*(-17 - 6*lmMst1 + 6*lmMt) + 6*
        pow2(lmMst1))*pow2(Mst1)*(Dmst12 - pow2(Mst2)) + ((9 + lmMgl*(20 + 6*
        lmMst1 - 6*lmMt) + 2*lmMt + lmMst1*(-22 + 6*lmMt) - 6*pow2(lmMst1))*
        pow2(Mgl) + 2*(2 + lmMgl*(5 + 2*lmMst1 - 2*lmMt) + 2*lmMst1*(-3 + lmMt)
        + lmMt - 2*pow2(lmMst1))*pow2(Mst1))*pow4(Mst2))) + Tbeta*(4*Dmst12*(1
        - 3*lmMgl + 3*lmMst1)*s2t*pow2(Mst1)*(Dmst12*Mt*pow2(MuSUSY) - Mt*pow2(
        Mst2)*pow2(MuSUSY) + 3*Dmst12*Mgl*s2t*pow2(Mst1)*pow2(Sbeta)) + 4*Mt*
        pow2(Sbeta)*(Dmst12*s2t*pow2(Mst1)*(pow2(Mst2)*(6*(8 - 25*lmMst1 +
        lmMgl*(23 + 6*lmMst1 - 6*lmMt) + 2*lmMt + 6*lmMst1*lmMt - 6*pow2(
        lmMst1))*pow2(Mgl) + 6*(3 - 14*lmMst1 + 4*lmMgl*(3 + lmMst1 - lmMt) +
        2*lmMt + 4*lmMst1*lmMt - 4*pow2(lmMst1))*pow2(Mst1) + (1 - 3*lmMgl + 3*
        lmMst1)*pow2(MuSUSY)) - Dmst12*(6*(8 - 25*lmMst1 + lmMgl*(23 + 6*lmMst1
        - 6*lmMt) + 2*lmMt + 6*lmMst1*lmMt - 6*pow2(lmMst1))*pow2(Mgl) + 3*(11
        + lmMst1*(22 - 4*lmMt) - 4*lmMgl*(4 + lmMst1 - lmMt) - 6*lmMt + 4*pow2(
        lmMst1))*pow2(Mst1) + (1 - 3*lmMgl + 3*lmMst1)*pow2(MuSUSY))) - 6*Mgl*
        Mt*(Dmst12*(lmMst1*(22 - 6*lmMt) - 5*lmMt + lmMgl*(-17 - 6*lmMst1 + 6*
        lmMt) + 6*pow2(lmMst1))*pow2(Mst1)*(Dmst12 - pow2(Mst2)) + ((9 + lmMgl*
        (20 + 6*lmMst1 - 6*lmMt) + 2*lmMt + lmMst1*(-22 + 6*lmMt) - 6*pow2(
        lmMst1))*pow2(Mgl) + 2*(2 + lmMgl*(5 + 2*lmMst1 - 2*lmMt) + 2*lmMst1*(-
        3 + lmMt) + lmMt - 2*pow2(lmMst1))*pow2(Mst1))*pow4(Mst2)))))*pow5(Mgl)
        ) + 3*oneLoopFlag*(s2t*pow2(Dmst12)*(-12*Mt*MuSUSY*pow2(Sbeta) + s2t*
        Tbeta*(pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 12*pow2(Mst1)*pow2(Sbeta))) +
        12*Mt*pow2(Sbeta)*(2*Dmst12*MuSUSY*s2t*pow2(Mst2) + Mt*Tbeta*(pow2(
        Dmst12) - 2*Dmst12*pow2(Mst2) + 4*(-lmMst1 + lmMt)*pow4(Mst2))))*pow8(
        Mst1))) + pow2(Mst1)*(40*Al4p*Mt*xMgl*pow4(Mgl)*(216*twoLoopFlag*
        xDmst12*pow2(Mst1)*pow3(Dmst12)*(Tbeta*(-2*Mt*pow2(s2t)*((-1 + 3*lmMgl
        - 3*lmMst1)*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 3*(1 + 2*lmMgl - 2*
        lmMst1)*pow2(Mst1)*pow2(Sbeta)) + 48*(5 - 2*lmMgl*(2 + lmMst1 - lmMt) +
        2*lmMst1*(4 + lmMst1 - lmMt) - 4*lmMt)*pow2(Sbeta)*pow3(Mt)) + MuSUSY*
        pow2(Sbeta)*(96*(1 - lmMgl + lmMst1)*s2t*pow2(Mt) + (1 - 3*lmMgl + 3*
        lmMst1)*pow2(Mst1)*pow3(s2t)))*pow4(Msq) + Al4p*threeLoopFlag*(xDmst12*
        pow2(Msq)*pow3(Dmst12)*(-2*Mt*MuSUSY*s2t*pow2(Mst1)*(MuSUSY*s2t*Tbeta*(
        41653 + 11784*lmMst1 + 48*lmMgl*(-313 + 180*lmMst1) - 12636*pow2(lmMgl)
        + 540*pow2(lmMst1))*pow2(Msq) + 12960*Dmsqst1*(29 - 10*lmMgl + 10*
        lmMst1)*Mt*pow2(Sbeta) + 864*Mt*(-1797 - 1711*lmMst1 + lmMgl*(1315 -
        198*lmMst1 - 96*lmMt) + 144*lmMt + 96*lmMst1*lmMt + 333*pow2(lmMgl) -
        103*pow2(lmMst1))*pow2(Msq)*pow2(Sbeta)) + Tbeta*pow2(Sbeta)*(2*Mt*
        pow2(Mst1)*(12960*Dmsqst1*(29 - 10*lmMgl + 10*lmMst1)*pow2(Mst1) - 108*
        (1695 + 8*lmMst1*(-125 + 12*lmMt) - 2*lmMgl*(-605 + 99*lmMst1 + 48*
        lmMt) + 333*pow2(lmMgl) - 103*pow2(lmMst1))*pow2(Msq)*pow2(Mst1) + (
        41653 + 11784*lmMst1 + 48*lmMgl*(-313 + 180*lmMst1) - 12636*pow2(lmMgl)
        + 540*pow2(lmMst1))*pow2(Msq)*pow2(MuSUSY))*pow2(s2t) - 288*(45*
        Dmsqst1*(307 + 4*lmMst1*(67 - 5*lmMt) - 10*lmMgl*(3 + 2*lmMst1 - 2*
        lmMt) - 238*lmMt + 20*pow2(lmMst1))*pow2(Mst1) + pow2(Msq)*(-2*(262 +
        lmMst1*(335 - 72*lmMt) - 96*lmMt + lmMgl*(-199 - 120*lmMst1 + 72*lmMt)
        + 120*pow2(lmMst1))*pow2(MuSUSY) + pow2(Mst1)*(72*(49 + 18*lmMst1 - 39*
        lmMt)*pow2(lmMgl) - 16*(181 + 108*lmMt)*pow2(lmMst1) + 8*lmMgl*(868 +
        287*lmMt + 9*lmMst1*(-59 + 63*lmMt) - 306*pow2(lmMst1) - 72*pow2(lmMt))
        + 8*lmMst1*(-4817 - 3*lmMt + 72*pow2(lmMt)) + 3*(-17921 + 7796*lmMt +
        384*pow2(lmMt)) + 504*pow3(lmMgl) + 648*pow3(lmMst1))))*pow3(Mt)) +
        MuSUSY*pow2(Msq)*(576*MuSUSY*Tbeta*(-262 - 335*lmMst1 + lmMgl*(199 +
        120*lmMst1 - 72*lmMt) + 96*lmMt + 72*lmMst1*lmMt - 120*pow2(lmMst1))*
        pow3(Mt) + (41653 + 11784*lmMst1 + 48*lmMgl*(-313 + 180*lmMst1) -
        12636*pow2(lmMgl) + 540*pow2(lmMst1))*pow2(Sbeta)*pow3(s2t)*pow4(Mst1))
        ) - 12960*Mt*xDmsqst1*pow2(Dmsqst1)*pow2(Mst1)*pow2(Mst2)*pow2(Sbeta)*(
        Dmst12*Mt*(2*(29 - 10*lmMgl + 10*lmMst1)*MuSUSY*s2t + Mt*Tbeta*(307 +
        268*lmMst1 - 10*lmMgl*(3 + 2*lmMst1 - 2*lmMt) - 238*lmMt - 20*lmMst1*
        lmMt + 20*pow2(lmMst1)))*pow2(Mst2) - pow2(Dmst12)*(2*(29 - 10*lmMgl +
        10*lmMst1)*Mt*MuSUSY*s2t + Tbeta*(307 + 4*lmMst1*(67 - 5*lmMt) - 10*
        lmMgl*(3 + 2*lmMst1 - 2*lmMt) - 238*lmMt + 20*pow2(lmMst1))*pow2(Mt) +
        (-29 + 10*lmMgl - 10*lmMst1)*Tbeta*pow2(Mst1)*pow2(s2t)) + (179 - 10*
        lmMgl + 10*lmMst1)*(1 + lmMst1 - lmMt)*Tbeta*pow2(Mt)*pow4(Mst2))) +
        Mt*pow2(Mst2)*(-216*twoLoopFlag*pow2(Mst1)*(24*Mt*pow2(Mst2)*(Dmst12*((
        1 - 2*lmMgl + 2*lmMst1)*MuSUSY*s2t + Mt*Tbeta*(2 + 6*lmMst1 - 3*lmMt -
        2*lmMst1*lmMt + lmMgl*(-3 - 2*lmMst1 + 2*lmMt) + 2*pow2(lmMst1))) - 2*
        Mt*Tbeta*(lmMgl*(1 + lmMst1 - lmMt) + lmMst1*(-2 + lmMt) + lmMt - pow2(
        lmMst1))*pow2(Mst2))*pow2(Sbeta) + pow2(Dmst12)*(12*Mt*((3 - 2*lmMgl +
        2*lmMst1)*MuSUSY*s2t + Mt*Tbeta*(8 - 2*lmMst1*(-5 + lmMt) - 5*lmMt +
        lmMgl*(-5 - 2*lmMst1 + 2*lmMt) + 2*pow2(lmMst1)))*pow2(Sbeta) + Tbeta*
        pow2(s2t)*(-((-1 + 3*lmMgl - 3*lmMst1)*pow2(MuSUSY)*(-1 + pow2(Sbeta)))
        + 12*(1 - 2*lmMgl + 2*lmMst1)*pow2(Mst1)*pow2(Sbeta))))*pow4(Msq) +
        Al4p*threeLoopFlag*pow2(Msq)*(Dmst12*s2t*pow2(Mst1)*(-(Dmst12*s2t*
        Tbeta*(12960*Dmsqst1*(29 - 10*lmMgl + 10*lmMst1)*pow2(Mst1) + pow2(Msq)
        *(432*(681 + 1506*lmMst1 - 96*(1 + lmMst1)*lmMt + 2*lmMgl*(-656 + 99*
        lmMst1 + 48*lmMt) - 333*pow2(lmMgl) + 103*pow2(lmMst1))*pow2(Mst1) + (
        41653 + 11784*lmMst1 + 48*lmMgl*(-313 + 180*lmMst1) - 12636*pow2(lmMgl)
        + 540*pow2(lmMst1))*pow2(MuSUSY)))*pow2(Sbeta)) + MuSUSY*(Dmst12*
        MuSUSY*s2t*Tbeta*(41653 + 11784*lmMst1 + 48*lmMgl*(-313 + 180*lmMst1) -
        12636*pow2(lmMgl) + 540*pow2(lmMst1))*pow2(Msq) + 432*Mt*(-60*Dmsqst1*(
        -29 + 10*lmMgl - 10*lmMst1)*(Dmst12 - pow2(Mst2)) + pow2(Msq)*(Dmst12*(
        -2961 - 1964*lmMst1 + 192*lmMt + 96*lmMst1*lmMt - 2*lmMgl*(-683 + 99*
        lmMst1 + 48*lmMt) + 333*pow2(lmMgl) - 103*pow2(lmMst1)) + 2*(-633 -
        1458*lmMst1 + 96*lmMt + 96*lmMst1*lmMt - 2*lmMgl*(-632 + 99*lmMst1 +
        48*lmMt) + 333*pow2(lmMgl) - 103*pow2(lmMst1))*pow2(Mst2)))*pow2(Sbeta)
        )) + 144*Tbeta*pow2(Mt)*(4*pow2(Msq)*pow2(MuSUSY)*(Dmst12*(262 +
        lmMst1*(335 - 72*lmMt) - 96*lmMt + lmMgl*(-199 - 120*lmMst1 + 72*lmMt)
        + 120*pow2(lmMst1))*(Dmst12 - pow2(Mst2)) - 2*(59 + 85*lmMst1 - lmMgl*(
        53 + 40*lmMst1 - 24*lmMt) - 24*(1 + lmMst1)*lmMt + 40*pow2(lmMst1))*
        pow4(Mst2)) + pow2(Sbeta)*(90*Dmsqst1*pow2(Mst1)*(Dmst12*(307 + 4*
        lmMst1*(67 - 5*lmMt) - 10*lmMgl*(3 + 2*lmMst1 - 2*lmMt) - 238*lmMt +
        20*pow2(lmMst1))*(Dmst12 - pow2(Mst2)) - 2*(69 + lmMst1*(59 - 10*lmMt)
        - 10*lmMgl*(1 + lmMst1 - lmMt) - 49*lmMt + 10*pow2(lmMst1))*pow4(Mst2))
        - pow2(Msq)*(pow2(Dmst12)*(4*(262 + lmMst1*(335 - 72*lmMt) - 96*lmMt +
        lmMgl*(-199 - 120*lmMst1 + 72*lmMt) + 120*pow2(lmMst1))*pow2(MuSUSY) +
        pow2(Mst1)*(40279 - 17847*lmMt - 6*(303 + 108*lmMst1 - 203*lmMt)*pow2(
        lmMgl) + (952 + 870*lmMt)*pow2(lmMst1) + lmMst1*(26599 + 804*lmMt -
        288*pow2(lmMt)) - 720*pow2(lmMt) + 2*lmMgl*(-1474 - 977*lmMt - 12*
        lmMst1*(-118 + 87*lmMt) + 615*pow2(lmMst1) + 144*pow2(lmMt)) - 190*
        pow3(lmMgl) - 392*pow3(lmMst1))) + 2*(Dmst12*pow2(Mst2)*(-2*(262 +
        lmMst1*(335 - 72*lmMt) - 96*lmMt + lmMgl*(-199 - 120*lmMst1 + 72*lmMt)
        + 120*pow2(lmMst1))*pow2(MuSUSY) + pow2(Mst1)*(13484 - 5541*lmMt + 6*(
        324 + 143*lmMt)*pow2(lmMst1) - 6*lmMgl*(666 + lmMgl*(285 + 108*lmMst1 -
        265*lmMt) + 57*lmMt + 4*lmMst1*(-59 + 102*lmMt) - 203*pow2(lmMst1) -
        48*pow2(lmMt)) - 432*pow2(lmMt) - 3*lmMst1*(-3979 + 260*lmMt + 96*pow2(
        lmMt)) - 314*pow3(lmMgl) - 256*pow3(lmMst1))) + (-4*(59 + lmMst1*(85 -
        24*lmMt) - 24*lmMt + lmMgl*(-53 - 40*lmMst1 + 24*lmMt) + 40*pow2(
        lmMst1))*pow2(MuSUSY) + pow2(Mst1)*(6677 - 2742*lmMt + 6*(313 + 143*
        lmMt)*pow2(lmMst1) - 6*lmMgl*(645 + lmMgl*(231 + 108*lmMst1 - 265*lmMt)
        - 72*lmMt + 12*lmMst1*(-9 + 34*lmMt) - 203*pow2(lmMst1) - 48*pow2(lmMt)
        ) - 288*pow2(lmMt) - 12*lmMst1*(-617 + 99*lmMt + 24*pow2(lmMt)) - 314*
        pow3(lmMgl) - 256*pow3(lmMst1)))*pow4(Mst2)))))))) + 6*xDmst12*pow3(
        Dmst12)*(Mt*pow4(Msq)*(270*oneLoopFlag*(-(Mt*Tbeta*pow2(s2t)*(pow2(
        MuSUSY)*(-1 + pow2(Sbeta)) + 18*pow2(Mst1)*pow2(Sbeta))) - pow2(Sbeta)*
        (-8*MuSUSY*s2t*pow2(Mt) + 8*Tbeta*pow3(Mt) + MuSUSY*pow2(Mst1)*pow3(
        s2t)))*pow6(Mst1) - 144*Al4p*twoLoopFlag*(-20*x*(Tbeta*((-12*Mgl*(-((
        lmMgl*(17 + 6*lmMst1 - 6*lmMt) + 5*lmMt + lmMst1*(-22 + 6*lmMt))*Mt) +
        (8 - 25*lmMst1 + lmMgl*(23 + 6*lmMst1 - 6*lmMt) + 2*lmMt + 6*lmMst1*
        lmMt)*Mgl*s2t + 6*(Mt - Mgl*s2t)*pow2(lmMst1))*pow2(Mt) + (1 - 3*lmMgl
        + 3*lmMst1)*pow2(Mst1)*(12*Mgl*Mt + s2t*pow2(Mst1))*pow2(s2t))*pow2(
        Sbeta) + s2t*pow2(Mt)*(2*(-1 + 3*lmMgl - 3*lmMst1)*pow2(MuSUSY)*(-1 +
        pow2(Sbeta)) + 48*(-2 - 9*lmMst1 + lmMgl*(7 + 2*lmMst1 - 2*lmMt) + 2*
        lmMt + 2*lmMst1*lmMt - 2*pow2(lmMst1))*pow2(Mst1)*pow2(Sbeta))) +
        MuSUSY*pow2(Sbeta)*(6*(-1 + 3*lmMgl - 3*lmMst1)*Mt*s2t*(2*Mgl*Mt + s2t*
        pow2(Mst1)) + 12*(lmMst1*(22 - 6*lmMt) - 5*lmMt + lmMgl*(-17 - 6*lmMst1
        + 6*lmMt) + 6*pow2(lmMst1))*pow3(Mt)))*pow5(Mgl) + Tbeta*pow2(Sbeta)*(
        2*Mgl*s2t*pow2(Mst1)*pow2(Mt)*(-((9 + 2*lmMst1)*pow2(Mst1)*pow2(MuSUSY)
        ) + 20*pow2(Mgl)*((1 + 4*lmMst1 - 2*(lmMgl + lmMt))*pow2(Mst1) + 4*(-1
        + lmMgl - lmMst1)*pow2(MuSUSY)) + 40*(1 + 2*lmMt)*pow4(Mst1)) - 5*Mt*
        pow2(s2t)*pow4(Mst1)*(pow2(Mgl)*(12*(3 + 2*lmMst1)*pow2(Mst1) + (1 + 2*
        lmMst1)*pow2(MuSUSY)) - 4*(lmMst1*pow2(Mst1)*pow2(MuSUSY) + 3*(1 + 6*
        lmMst1)*pow4(Mst1))) + (40*(-1 + 2*lmMst1 + 2*lmMt)*pow3(Mt) - 20*((2 +
        lmMst1)*Mgl*pow2(Mst1) + (1 - 2*lmMgl + 2*lmMst1)*pow3(Mgl))*pow3(s2t))
        *pow6(Mst1)) + MuSUSY*pow2(Mst1)*(Mt*MuSUSY*s2t*Tbeta*(2*(9 + 2*lmMst1)
        *Mgl*Mt*pow2(Mst1) + 5*(1 + 2*lmMst1)*s2t*pow2(Mgl)*pow2(Mst1) + 160*(1
        - lmMgl + lmMst1)*Mt*pow3(Mgl) - 20*lmMst1*s2t*pow4(Mst1)) + 10*pow2(
        Sbeta)*(pow3(Mgl)*(-3*(1 + 2*lmMgl - 2*lmMst1)*Mt*pow2(Mst1)*pow2(s2t)
        + 8*(18 + 26*lmMst1 - 13*lmMt - 6*lmMst1*lmMt + lmMgl*(-13 - 6*lmMst1 +
        6*lmMt) + 6*pow2(lmMst1))*pow3(Mt)) + 2*(9 - 8*lmMst1)*s2t*pow2(Mt)*
        pow4(Mst1) + Mgl*(4*(1 + lmMt)*pow2(Mst1)*pow3(Mt) + 3*(3 + 2*lmMst1)*
        Mt*pow2(s2t)*pow4(Mst1)) - pow2(Mgl)*(4*s2t*pow2(Mst1)*pow2(Mt) + (3 +
        2*lmMst1)*pow3(s2t)*pow4(Mst1)) + (1 + 2*lmMst1)*pow3(s2t)*pow6(Mst1)))
        )) + threeLoopFlag*pow2(Al4p)*pow2(Mst1)*(5*xDmsqst1*pow2(Dmsqst1)*(
        5760*Mgl*Mt*(pow2(Mst1)*(4*s2t*pow2(Mt)*(8*Tbeta*pow2(MuSUSY)*(-1 +
        pow2(Sbeta)) + 72*Mgl*MuSUSY*pow2(Sbeta) + (-15 - 38*lmMst1 + 38*lmMt)*
        Tbeta*pow2(Mst1)*pow2(Sbeta)) - 3*Mt*pow2(s2t)*(-4*Mgl*Tbeta*pow2(
        MuSUSY)*(-1 + pow2(Sbeta)) + 7*(MuSUSY - 3*Mgl*Tbeta)*pow2(Mst1)*pow2(
        Sbeta)) + (-2*(91 + 103*lmMst1 - 103*lmMt)*MuSUSY + 27*(17 + 25*lmMst1
        - 25*lmMt)*Mgl*Tbeta)*pow2(Sbeta)*pow3(Mt) + 2*pow2(Mst1)*(3*Mgl*MuSUSY
        + 5*Tbeta*pow2(Mst1))*pow2(Sbeta)*pow3(s2t)) + 4*pow2(Mgl)*(MuSUSY*
        pow2(Sbeta)*(6*(-3 + lmMgl - lmMst1)*Mt*pow2(Mst1)*pow2(s2t) + 6*(31 +
        27*lmMst1 - 24*lmMt - 2*lmMst1*lmMt + lmMgl*(-3 - 2*lmMst1 + 2*lmMt) +
        2*pow2(lmMst1))*pow3(Mt)) + Tbeta*(-(s2t*pow2(Mt)*(-2*(-3 + lmMgl -
        lmMst1)*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 3*(118 + lmMst1*(73 - 6*lmMt)
        - 60*lmMt + lmMgl*(-13 - 6*lmMst1 + 6*lmMt) + 6*pow2(lmMst1))*pow2(
        Mst1)*pow2(Sbeta))) + (3 - lmMgl + lmMst1)*pow2(Sbeta)*pow3(s2t)*pow4(
        Mst1)))) + pow4(Mst1)*(36*Tbeta*pow2(Mt)*pow2(s2t)*((215 - 72*shiftst1
        - 208*shiftst2 + 16*lmMst1*(-5 + 3*shiftst1 + 2*shiftst2))*pow2(MuSUSY)
        *(-1 + pow2(Sbeta)) + 5*(201 - 512*shiftst1 + 128*shiftst2 - 64*lmMst1*
        (3 - 4*shiftst1 + shiftst2))*pow2(Mst1)*pow2(Sbeta)) + 5*pow2(Sbeta)*(
        27*MuSUSY*s2t*(1389 - 256*shiftst2)*pow3(Mt) + 36*Mt*MuSUSY*(33 - 128*
        shiftst1 + 64*shiftst2 - 32*lmMst1*(1 - 2*shiftst1 + shiftst2))*pow2(
        Mst1)*pow3(s2t) + 16*(-2657 + 2928*lmMst1 - 2928*lmMt + 432*shiftst2)*
        Tbeta*pow4(Mt) - 576*(-2 + lmMst1)*(shiftst1 - shiftst2)*Tbeta*pow4(
        Mst1)*pow4(s2t)))) + pow2(Msq)*(3*MuSUSY*pow2(Sbeta)*(s2t*pow2(Mst1)*(
        Mt*(225*Dmsqst1*(15360*pow2(Mgl) - 83*pow2(Mst1)) - 8*pow2(Msq)*(8*(
        4113 + 90*lmMgl - 1792*lmMst1 + 120*lmMt - 48*pow2(lmMst1))*pow2(Mgl) +
        5*(2615 + 2600*lmMst1 - 1312*pow2(lmMst1))*pow2(Mst1))) - 960*Mgl*s2t*(
        -30*Dmsqst1*(8*(-3 + lmMgl - lmMst1)*pow2(Mgl) + pow2(Mst1)) + pow2(
        Msq)*(-2*(227 + 144*lmMgl - 120*lmMst1 - 20*lmMgl*lmMst1 + 27*pow2(
        lmMgl) + pow2(lmMst1))*pow2(Mgl) + (66 + 77*lmMst1 + 41*pow2(lmMst1))*
        pow2(Mst1))))*pow2(Mt) + 80*(Mt*(15*Dmsqst1*(48*pow2(Mgl) + (3 - 24*
        shiftst1 + 12*shiftst2 - 8*lmMst1*(1 - 2*shiftst1 + shiftst2))*pow2(
        Mst1)) + 2*pow2(Msq)*((200 + 18*lmMgl + 265*lmMst1 + 91*pow2(lmMst1))*
        pow2(Mgl) - (16 + 60*shiftst1 - 30*shiftst2 + lmMst1*(115 - 120*
        shiftst1 + 60*shiftst2) + 41*pow2(lmMst1))*pow2(Mst1)))*pow3(s2t) + 3*
        shiftst3*pow2(Msq)*(16*(7 - 3*lmMst1)*s2t*pow3(Mt) + (-13 + 14*lmMst1)*
        Mt*pow2(Mst1)*pow3(s2t)))*pow4(Mst1) + 640*Mgl*(360*Dmsqst1*(31 +
        lmMst1*(27 - 2*lmMt) - 24*lmMt + lmMgl*(-3 - 2*lmMst1 + 2*lmMt) + 2*
        pow2(lmMst1))*pow2(Mgl) - 240*Dmsqst1*(11 + 20*lmMst1 - 20*lmMt)*pow2(
        Mst1) - (-233 + 47*lmMt + lmMst1*(-3 + 2*lmMt) + 15*pow2(lmMst1) + 24*
        pow2(lmMt))*pow2(Msq)*pow2(Mst1) - 4*pow2(Mgl)*pow2(Msq)*(5845 - 2937*
        lmMt - 27*(16 + lmMgl + 6*lmMst1 - 9*lmMt)*pow2(lmMgl) + (346 + 183*
        lmMt)*pow2(lmMst1) + lmMst1*(4753 - 84*lmMt - 72*pow2(lmMt)) - 156*
        pow2(lmMt) + lmMgl*(-783 - 426*lmMst1*(-1 + lmMt) - 52*lmMt + 273*pow2(
        lmMst1) + 72*pow2(lmMt)) - 84*pow3(lmMst1)))*pow4(Mt)) + Tbeta*(-60*
        pow2(Mst1)*pow2(Mt)*pow2(s2t)*((15*Dmsqst1*(384*pow2(Mgl) + (1 - 16*
        lmMst1)*pow2(Mst1)) + 2*pow2(Msq)*(2*(-425 + 18*lmMgl + 83*lmMst1 + 91*
        pow2(lmMst1))*pow2(Mgl) + (121 - 176*lmMst1 - 164*pow2(lmMst1))*pow2(
        Mst1)))*pow2(MuSUSY) + 24*shiftst1*(Dmsqst1*(3 - 2*lmMst1) + (1 - 2*
        lmMst1)*pow2(Msq))*pow2(Mst1)*(80*pow2(Mst1) + 3*pow2(MuSUSY))*pow2(
        Sbeta)) + 64*Mgl*(12*Mgl*Mt*(323 - 75*lmMgl + 383*lmMst1 - 100*lmMt +
        32*pow2(lmMst1))*pow2(Msq) + 240*s2t*pow2(Mgl)*(15*Dmsqst1*(3 - lmMgl +
        lmMst1) + (-130 + lmMgl*(159 - 20*lmMst1) - 162*lmMst1 + 27*pow2(lmMgl)
        + pow2(lmMst1))*pow2(Msq)) - s2t*(18000*Dmsqst1 + (-365 + 1049*lmMst1 +
        123*pow2(lmMst1))*pow2(Msq))*pow2(Mst1))*pow2(MuSUSY)*pow3(Mt) + 144*
        pow4(Mst1)*(shiftst3*pow2(Msq)*(2*pow2(Mt)*pow2(s2t)*((-29 + 18*lmMst1)
        *pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 5*(-83 + 58*lmMst1)*pow2(Mst1)*pow2(
        Sbeta)) - 80*(-7 + 3*lmMst1)*pow2(Sbeta)*pow4(Mt) + 5*(1 - 2*lmMst1)*
        pow2(Sbeta)*pow4(Mst1)*pow4(s2t)) + 10*(Dmsqst1*(3 - 2*lmMst1) + (1 -
        2*lmMst1)*pow2(Msq))*(pow2(Mt)*pow2(s2t)*(20*shiftst2*pow2(Mst1)*pow2(
        Sbeta) + pow2(MuSUSY)*(3*shiftst1 + 2*shiftst2 - 2*shiftst2*pow2(Sbeta)
        )) - 5*shiftst2*pow2(Sbeta)*pow4(Mst1)*pow4(s2t))) + pow2(Sbeta)*(60*
        pow2(Mst1)*pow2(Mt)*pow2(s2t)*(15*Dmsqst1*(-96*pow2(Mgl)*(3*pow2(Mst1)
        - 4*pow2(MuSUSY)) + (1 - 16*lmMst1)*pow2(Mst1)*pow2(MuSUSY) + 6*(9 -
        32*lmMst1)*pow4(Mst1)) + 2*pow2(Msq)*((121 - 176*lmMst1 - 164*pow2(
        lmMst1))*pow2(Mst1)*pow2(MuSUSY) + 2*pow2(Mgl)*(12*(3 - 18*lmMgl + 397*
        lmMst1 - 48*lmMt + 107*pow2(lmMst1))*pow2(Mst1) + (-425 + 18*lmMgl +
        83*lmMst1 + 91*pow2(lmMst1))*pow2(MuSUSY)) - 12*(-119 + 250*lmMst1 +
        246*pow2(lmMst1))*pow4(Mst1))) + Mgl*(960*Mt*(60*Dmsqst1*(2*(3 - lmMgl
        + lmMst1)*pow2(Mgl) + 3*pow2(Mst1)) + pow2(Msq)*(-4*(-11 + lmMgl*(154 -
        20*lmMst1) - 148*lmMst1 + 27*pow2(lmMgl) + pow2(lmMst1))*pow2(Mgl) + (
        215 + 118*lmMst1 + 41*pow2(lmMst1))*pow2(Mst1)))*pow3(s2t)*pow4(Mst1) +
        64*s2t*pow3(Mt)*(900*Dmsqst1*(20*pow2(Mst1)*pow2(MuSUSY) + 4*pow2(Mgl)*
        (3*(-67 - 43*lmMst1 + lmMgl*(5 + 2*lmMst1 - 2*lmMt) + 38*lmMt + 2*
        lmMst1*lmMt - 2*pow2(lmMst1))*pow2(Mst1) + (-3 + lmMgl - lmMst1)*pow2(
        MuSUSY)) + (15 + 2*lmMst1 - 2*lmMt)*pow4(Mst1)) + pow2(Msq)*((-365 +
        1049*lmMst1 + 123*pow2(lmMst1))*pow2(Mst1)*pow2(MuSUSY) + 30*pow2(Mgl)*
        ((2722 + 965*lmMst1 - 6*(200 + 23*lmMst1)*lmMt + lmMgl*(67 - 258*lmMst1
        + 178*lmMt) + 108*pow2(lmMgl) + 62*pow2(lmMst1) + 48*pow2(lmMt))*pow2(
        Mst1) - 8*(-130 + lmMgl*(159 - 20*lmMst1) - 162*lmMst1 + 27*pow2(lmMgl)
        + pow2(lmMst1))*pow2(MuSUSY)) - 30*(-223 + 342*lmMt + lmMst1*(-442 + 8*
        lmMt) + 60*pow2(lmMst1) + 96*pow2(lmMt))*pow4(Mst1)))) + 8*(25*Dmsqst1*
        (2592*(25 + 47*lmMst1 - 47*lmMt)*pow2(Mgl)*pow2(Mst1) + (-907 + 96*
        lmMst1 - 96*lmMt)*pow4(Mst1)) - 4*pow2(Msq)*(12*pow2(Mgl)*(5*(1140 +
        lmMst1*(747 - 32*lmMt) - 747*lmMt + 16*pow2(lmMst1))*pow2(Mst1) + 2*(
        323 - 75*lmMgl + 383*lmMst1 - 100*lmMt + 32*pow2(lmMst1))*pow2(MuSUSY))
        + 5*(2243 + 2052*lmMt - 24*lmMst1*(129 + 10*lmMt) + 612*(pow2(lmMst1) +
        pow2(lmMt)))*pow4(Mst1)))*pow4(Mt) + 7200*shiftst1*(Dmsqst1*(3 - 2*
        lmMst1) + (1 - 2*lmMst1)*pow2(Msq))*pow4(s2t)*pow8(Mst1))))))))/(77760.
        *Tbeta*pow2(Sbeta)*pow4(Msq)*pow6(Mst2)*pow8(Mst1));
}

/**
 * 	@return The off-diagonal (1, 2) = (2, 1) matrix element of the Higgs mass matrix as a double for the hierarchy 'H9q2'
 */
double H9q2::getS12() const {
   return (MuSUSY*(-6*xDmst12*pow3(Dmst12)*(threeLoopFlag*pow2(Al4p)*(-3*s2t*(
        104600 + 104000*lmMst1 - 52480*pow2(lmMst1) + (18675*Dmsqst1)/pow2(Msq)
        + (64*pow2(Mgl)*(4113 + 90*lmMgl - 1792*lmMst1 + 120*lmMt - 48*pow2(
        lmMst1) - (54000*Dmsqst1)/pow2(Msq)))/pow2(Mst1))*pow3(Mt) + 720*
        shiftst3*(16*(7 - 3*lmMst1)*s2t*pow3(Mt) + (-13 + 14*lmMst1)*Mt*pow2(
        Mst1)*pow3(s2t)) + (240*Mt*((15*Dmsqst1*(48*pow2(Mgl) + (3 - 24*
        shiftst1 + 12*shiftst2 - 8*lmMst1*(1 - 2*shiftst1 + shiftst2))*pow2(
        Mst1)) + 2*pow2(Msq)*((200 + 18*lmMgl + 265*lmMst1 + 91*pow2(lmMst1))*
        pow2(Mgl) - (16 + 60*shiftst1 - 30*shiftst2 + lmMst1*(115 - 120*
        shiftst1 + 60*shiftst2) + 41*pow2(lmMst1))*pow2(Mst1)))*pow3(s2t) + (4*
        Mgl*Mt*(-3*pow2(Mst1)*(-30*Dmsqst1*(8*(-3 + lmMgl - lmMst1)*pow2(Mgl) +
        pow2(Mst1)) + pow2(Msq)*(-2*(227 + 144*lmMgl - 120*lmMst1 - 20*lmMgl*
        lmMst1 + 27*pow2(lmMgl) + pow2(lmMst1))*pow2(Mgl) + (66 + 77*lmMst1 +
        41*pow2(lmMst1))*pow2(Mst1)))*pow2(s2t) + 2*pow2(Mt)*(360*Dmsqst1*(31 +
        lmMst1*(27 - 2*lmMt) - 24*lmMt + lmMgl*(-3 - 2*lmMst1 + 2*lmMt) + 2*
        pow2(lmMst1))*pow2(Mgl) - 240*Dmsqst1*(11 + 20*lmMst1 - 20*lmMt)*pow2(
        Mst1) - (-233 + 47*lmMt + lmMst1*(-3 + 2*lmMt) + 15*pow2(lmMst1) + 24*
        pow2(lmMt))*pow2(Msq)*pow2(Mst1) - 4*pow2(Mgl)*pow2(Msq)*(5845 - 2937*
        lmMt - 27*(16 + lmMgl + 6*lmMst1 - 9*lmMt)*pow2(lmMgl) + (346 + 183*
        lmMt)*pow2(lmMst1) + lmMst1*(4753 - 84*lmMt - 72*pow2(lmMt)) - 156*
        pow2(lmMt) + lmMgl*(-783 - 426*lmMst1*(-1 + lmMt) - 52*lmMt + 273*pow2(
        lmMst1) + 72*pow2(lmMt)) - 84*pow3(lmMst1)))))/pow4(Mst1)))/pow2(Msq) +
        (Mt*(8*Mt*MuSUSY*pow2(Msq)*(-12*pow2(Mgl)*(7200*Dmsqst1*pow2(Mst1)*
        pow2(s2t) + pow2(Msq)*(16*(-323 + 75*lmMgl - 383*lmMst1 + 100*lmMt -
        32*pow2(lmMst1))*pow2(Mt) + 5*(-425 + 18*lmMgl + 83*lmMst1 + 91*pow2(
        lmMst1))*pow2(Mst1)*pow2(s2t))) + Mt*s2t*(-16*Mgl*(18000*Dmsqst1 + (-
        365 + 1049*lmMst1 + 123*pow2(lmMst1))*pow2(Msq))*pow2(Mst1) + 3840*(15*
        Dmsqst1*(3 - lmMgl + lmMst1) + (-130 + lmMgl*(159 - 20*lmMst1) - 162*
        lmMst1 + 27*pow2(lmMgl) + pow2(lmMst1))*pow2(Msq))*pow3(Mgl)) - 3*(15*
        Dmsqst1*(5 - 72*shiftst1 - 48*shiftst2 + 16*lmMst1*(-5 + 3*shiftst1 +
        2*shiftst2)) - 2*(-605 + 180*shiftst1 + 120*shiftst2 + 348*shiftst3 -
        8*lmMst1*(-110 + 45*shiftst1 + 30*shiftst2 + 27*shiftst3) + 820*pow2(
        lmMst1))*pow2(Msq))*pow2(s2t)*pow4(Mst1)) - 45*xDmsqst1*pow2(Dmsqst1)*(
        -5120*Mt*(2*(3 - lmMgl + lmMst1)*Mt*MuSUSY*s2t + 3*Tbeta*(31 + 27*
        lmMst1 - 24*lmMt - 2*lmMst1*lmMt + lmMgl*(-3 - 2*lmMst1 + 2*lmMt) + 2*
        pow2(lmMst1))*pow2(Mt) + 3*(-3 + lmMgl - lmMst1)*Tbeta*pow2(Mst1)*pow2(
        s2t))*pow3(Mgl) + 640*Mgl*pow2(Mst1)*(32*s2t*(2*MuSUSY - 9*Mgl*Tbeta)*
        pow2(Mt) + 3*Mt*(8*Mgl*MuSUSY + 7*Tbeta*pow2(Mst1))*pow2(s2t) + 2*(91 +
        103*lmMst1 - 103*lmMt)*Tbeta*pow3(Mt) - 6*Mgl*Tbeta*pow2(Mst1)*pow3(
        s2t)) + s2t*(8*Mt*MuSUSY*s2t*(215 - 72*shiftst1 - 208*shiftst2 + 16*
        lmMst1*(-5 + 3*shiftst1 + 2*shiftst2)) + Tbeta*(15*(-1389 + 256*
        shiftst2)*pow2(Mt) - 20*(33 - 128*shiftst1 + 64*shiftst2 - 32*lmMst1*(1
        - 2*shiftst1 + shiftst2))*pow2(Mst1)*pow2(s2t)))*pow4(Mst1))))/(Tbeta*
        pow4(Msq)*pow4(Mst1))) + (Mt*(-270*oneLoopFlag*s2t*(-2*Mt*MuSUSY*s2t -
        8*Tbeta*pow2(Mt) + Tbeta*pow2(Mst1)*pow2(s2t)) - (288*Al4p*twoLoopFlag*
        (5*Mt*pow2(Mst1)*(32*(1 - lmMgl + lmMst1)*Mt*MuSUSY*s2t + 8*Tbeta*(18 +
        lmMst1*(26 - 6*lmMt) - 13*lmMt + lmMgl*(-13 - 6*lmMst1 + 6*lmMt) + 6*
        pow2(lmMst1))*pow2(Mt) - 3*(1 + 2*lmMgl - 2*lmMst1)*Tbeta*pow2(Mst1)*
        pow2(s2t))*pow3(Mgl) + 5*s2t*pow2(Mgl)*((1 + 2*lmMst1)*Mt*MuSUSY*s2t -
        4*Tbeta*pow2(Mt) - (3 + 2*lmMst1)*Tbeta*pow2(Mst1)*pow2(s2t))*pow4(
        Mst1) + Mgl*Mt*(2*(9 + 2*lmMst1)*Mt*MuSUSY*s2t + 20*(1 + lmMt)*Tbeta*
        pow2(Mt) + 15*(3 + 2*lmMst1)*Tbeta*pow2(Mst1)*pow2(s2t))*pow4(Mst1) -
        20*Mt*x*(2*(1 - 3*lmMgl + 3*lmMst1)*Mt*MuSUSY*s2t - 6*Tbeta*(-22*lmMst1
        + lmMgl*(17 + 6*lmMst1 - 6*lmMt) + 5*lmMt + 6*lmMst1*lmMt - 6*pow2(
        lmMst1))*pow2(Mt) + 3*(-1 + 3*lmMgl - 3*lmMst1)*Tbeta*pow2(Mst1)*pow2(
        s2t))*pow5(Mgl) + 120*(1 - 3*lmMgl + 3*lmMst1)*s2t*Tbeta*x*pow2(Mt)*
        pow6(Mgl) - 5*s2t*(4*lmMst1*Mt*MuSUSY*s2t + 2*(-9 + 8*lmMst1)*Tbeta*
        pow2(Mt) - (1 + 2*lmMst1)*Tbeta*pow2(Mst1)*pow2(s2t))*pow6(Mst1)))/
        pow6(Mst1)))/Tbeta) + ((-9*threeLoopFlag*pow2(Al4p)*(Mt*(1200*xDmsqst1*
        pow2(Dmsqst1)*pow2(Mst1)*pow2(Mst2)*(-2*Dmst12*Mt*pow2(Mst2)*(Mt*(16*
        Mgl*(-20*MuSUSY*s2t - 3*(11 + 47*lmMst1 - 47*lmMt)*Mt*Tbeta + 90*Mgl*
        s2t*Tbeta)*pow2(Mst1) - 64*(2*(3 - lmMgl + lmMst1)*MuSUSY*s2t + 3*Mt*
        Tbeta*(31 + 27*lmMst1 - 24*lmMt - 2*lmMst1*lmMt + lmMgl*(-3 - 2*lmMst1
        + 2*lmMt) + 2*pow2(lmMst1)))*pow3(Mgl)) - s2t*(-16*(-2 + lmMst1)*
        MuSUSY*s2t*(shiftst1 - shiftst2) + 3*Mt*(-87 - 32*lmMst1*(-1 +
        shiftst2) + 64*shiftst2)*Tbeta)*pow4(Mst1)) + pow2(Dmst12)*(Mt*(32*Mgl*
        pow2(Mst1)*(3*Mt*s2t*(2*MuSUSY - 9*Mgl*Tbeta) + (29 - 19*lmMst1 + 19*
        lmMt)*Tbeta*pow2(Mt) + 3*(2*Mgl*MuSUSY + 5*Tbeta*pow2(Mst1))*pow2(s2t))
        - 64*(4*(3 - lmMgl + lmMst1)*Mt*MuSUSY*s2t + 6*Tbeta*(31 + 27*lmMst1 -
        24*lmMt - 2*lmMst1*lmMt + lmMgl*(-3 - 2*lmMst1 + 2*lmMt) + 2*pow2(
        lmMst1))*pow2(Mt) + 3*(-3 + lmMgl - lmMst1)*Tbeta*pow2(Mst1)*pow2(s2t))
        *pow3(Mgl)) + s2t*(Mt*MuSUSY*s2t*(33 + 16*lmMst1*(-2 + shiftst1 +
        shiftst2) - 32*(shiftst1 + shiftst2)) + Tbeta*(-321*pow2(Mt) - 24*(-2 +
        lmMst1)*(shiftst1 - shiftst2)*pow2(Mst1)*pow2(s2t)))*pow4(Mst1)) + 192*
        Tbeta*pow2(Mt)*(2*Mt*((1 + 7*lmMst1 - 7*lmMt)*Mgl*pow2(Mst1) + (18 -
        lmMgl + lmMst1)*(1 + lmMst1 - lmMt)*pow3(Mgl)) + (-2 + lmMst1)*s2t*(
        shiftst1 - shiftst2)*pow4(Mst1))*pow4(Mst2)) - 3*z3*(80*Dmst12*pow2(Mt)
        *pow4(Mst2)*(-96*xMgl*(121*Mt*MuSUSY + 93*s2t*Tbeta*pow2(Mst1))*pow4(
        Mgl)*pow4(Msq) + 32*pow2(Mgl)*pow2(Mst1)*((7*Mt*MuSUSY - 6*s2t*Tbeta*
        pow2(Mst1))*pow4(Msq) - 2*Mgl*(1440*Dmsqst1*Mt*Tbeta*(Dmsqst1*xDmsqst1
        + pow2(Msq)) + (81*MuSUSY*s2t + 50*Mt*Tbeta)*pow4(Msq))) - 32*Mgl*(
        1680*Mt*Tbeta*xDmsqst1*pow2(Dmsqst1) + 960*Dmsqst1*Mt*Tbeta*pow2(Msq) +
        (-7*MuSUSY*s2t + 248*Mt*Tbeta)*pow4(Msq))*pow4(Mst1) + 3*s2t*Tbeta*(
        525*xDmsqst1*pow2(Dmsqst1) + 280*Dmsqst1*pow2(Msq) - 248*pow4(Msq))*
        pow6(Mst1)) - 40*Mt*pow2(Dmst12)*pow2(Mst2)*(-24*xMgl*(-372*Mt*s2t*
        Tbeta*pow2(Mst1) + 968*MuSUSY*pow2(Mt) + 259*MuSUSY*pow2(Mst1)*pow2(
        s2t))*pow4(Mgl)*pow4(Msq) + 8*pow2(Mgl)*pow2(Mst1)*(-23040*Dmsqst1*Mgl*
        Tbeta*(Dmsqst1*xDmsqst1 + pow2(Msq))*pow2(Mt) + MuSUSY*(7*pow2(Mt) -
        90*pow2(Mst1)*pow2(s2t))*pow4(Msq) + 4*Mgl*(162*Mt*MuSUSY*s2t + 2192*
        Tbeta*pow2(Mt) + 243*Tbeta*pow2(Mst1)*pow2(s2t))*pow4(Msq)) + 16*Mgl*
        Tbeta*(-960*Dmsqst1*(Dmsqst1*xDmsqst1 - 2*pow2(Msq))*pow2(Mt) + (1136*
        pow2(Mt) - 21*pow2(Mst1)*pow2(s2t))*pow4(Msq))*pow4(Mst1) + 3*s2t*(35*(
        7*MuSUSY*s2t - 15*Mt*Tbeta)*xDmsqst1*pow2(Dmsqst1) + 70*Dmsqst1*(2*
        MuSUSY*s2t - 7*Mt*Tbeta)*pow2(Msq) + 48*(MuSUSY*s2t - 11*Mt*Tbeta)*
        pow4(Msq))*pow6(Mst1)) - xDmst12*pow3(Dmst12)*(Mt*(1920*pow2(Mst1)*
        pow3(Mgl)*(3840*Dmsqst1*Tbeta*(Dmsqst1*xDmsqst1 + pow2(Msq))*pow2(Mt) -
        (432*Mt*MuSUSY*s2t + 3056*Tbeta*pow2(Mt) + 81*Tbeta*pow2(Mst1)*pow2(
        s2t))*pow4(Msq)) - 128*Mgl*(4800*Dmsqst1*Tbeta*(5*Dmsqst1*xDmsqst1 + 8*
        pow2(Msq))*pow2(Mt) + (-6*Mt*MuSUSY*s2t + 440*Tbeta*pow2(Mt) - 105*
        Tbeta*pow2(Mst1)*pow2(s2t))*pow4(Msq))*pow4(Mst1)) + 15*s2t*(35*
        Dmsqst1*pow2(Msq)*(8*Mt*MuSUSY*s2t + 31*Tbeta*pow2(Mt) + 16*Tbeta*pow2(
        Mst1)*pow2(s2t)) + 35*xDmsqst1*pow2(Dmsqst1)*(-40*Mt*MuSUSY*s2t + 159*
        Tbeta*pow2(Mt) + 28*Tbeta*pow2(Mst1)*pow2(s2t)) + 8*(22*Mt*MuSUSY*s2t +
        85*Tbeta*pow2(Mt) + 24*Tbeta*pow2(Mst1)*pow2(s2t))*pow4(Msq))*pow6(
        Mst1) + pow4(Msq)*(480*xMgl*pow4(Mgl)*(-2976*s2t*Tbeta*pow2(Mst1)*pow2(
        Mt) + 1036*Mt*MuSUSY*pow2(Mst1)*pow2(s2t) + 1936*MuSUSY*pow3(Mt) - 259*
        Tbeta*pow3(s2t)*pow4(Mst1)) - 320*pow2(Mgl)*(-42*MuSUSY*pow2(Mst1)*
        pow3(Mt) - Mt*s2t*(45*MuSUSY*s2t + 7*Mt*Tbeta)*pow4(Mst1) + 45*Tbeta*
        pow3(s2t)*pow6(Mst1)))) + 2560*Mgl*pow3(Mt)*(-20*Tbeta*pow2(Mgl)*pow2(
        Mst1)*(120*xDmsqst1*pow2(Dmsqst1) + 48*Dmsqst1*pow2(Msq) - 19*pow4(Msq)
        ) + MuSUSY*(7*Mgl*pow2(Mst1) - 242*xMgl*pow3(Mgl))*pow4(Msq) - 16*
        Tbeta*(60*xDmsqst1*pow2(Dmsqst1) + 30*Dmsqst1*pow2(Msq) + pow4(Msq))*
        pow4(Mst1))*pow6(Mst2)) - 320*z2*(3*xDmsqst1*xDmst12*pow2(Dmsqst1)*
        pow2(Mst1)*pow3(Dmst12)*(30*s2t*pow2(Mgl)*pow2(Mst1)*(-4*Mt*MuSUSY*s2t
        + 84*Tbeta*pow2(Mt) + Tbeta*pow2(Mst1)*pow2(s2t)) - 80*Mgl*Mt*pow2(
        Mst1)*(7*Mt*MuSUSY*s2t + 3*(39 + 10*lmMst1 - 10*lmMt)*Tbeta*pow2(Mt) +
        3*Tbeta*pow2(Mst1)*pow2(s2t)) + 160*Mt*(2*Mt*MuSUSY*s2t + 36*(4 +
        lmMst1 - lmMt)*Tbeta*pow2(Mt) - 3*Tbeta*pow2(Mst1)*pow2(s2t))*pow3(Mgl)
        + s2t*(4*Mt*MuSUSY*s2t*(-15 + 3*shiftst1 + 2*shiftst2) + 480*Tbeta*
        pow2(Mt) + 5*(9 - 8*shiftst1 + 4*shiftst2)*Tbeta*pow2(Mst1)*pow2(s2t))*
        pow4(Mst1)) + 4*xMgl*pow4(Mgl)*(6*Mt*pow2(Dmst12)*pow2(Mst2)*(150*
        Dmsqst1*Mt*s2t*Tbeta*(Dmsqst1*xDmsqst1 + pow2(Msq))*pow2(Mst1) + (146 -
        27*lmMgl + 27*lmMst1)*Mt*s2t*Tbeta*pow2(Mst1)*pow4(Msq) - MuSUSY*(24*(
        17 - 4*lmMgl + 4*lmMst1)*pow2(Mt) + (49 - 17*lmMgl + 17*lmMst1)*pow2(
        Mst1)*pow2(s2t))*pow4(Msq)) - 3*xDmst12*pow3(Dmst12)*(300*Dmsqst1*s2t*
        Tbeta*(Dmsqst1*xDmsqst1 + pow2(Msq))*pow2(Mst1)*pow2(Mt) - pow4(Msq)*(
        8*(-170 + 27*lmMgl - 27*lmMst1)*s2t*Tbeta*pow2(Mst1)*pow2(Mt) + 4*(49 -
        17*lmMgl + 17*lmMst1)*Mt*MuSUSY*pow2(Mst1)*pow2(s2t) + 48*(17 - 4*lmMgl
        + 4*lmMst1)*MuSUSY*pow3(Mt) + (-49 + 17*lmMgl - 17*lmMst1)*Tbeta*pow3(
        s2t)*pow4(Mst1))) - 12*Dmst12*pow2(Mt)*(75*Dmsqst1*s2t*Tbeta*(Dmsqst1*
        xDmsqst1 + pow2(Msq))*pow2(Mst1) - (12*(17 - 4*lmMgl + 4*lmMst1)*Mt*
        MuSUSY + (194 - 27*lmMgl + 27*lmMst1)*s2t*Tbeta*pow2(Mst1))*pow4(Msq))*
        pow4(Mst2) + 32*(47 - 12*lmMgl + 12*lmMst1)*MuSUSY*pow3(Mt)*pow4(Msq)*
        pow6(Mst2)))) + pow2(Mst1)*(64*MuSUSY*pow2(Msq)*pow2(Mst2)*pow2(Mt)*(-
        10*Dmst12*pow2(Mst2)*(16*(16 - 9*lmMgl + 41*lmMst1 - 12*lmMt + 26*z2 +
        4*pow2(lmMst1))*pow2(Mgl)*pow2(Msq)*pow2(Mt) + Mt*s2t*(4*Mgl*pow2(Msq)*
        (2*(-22 - 296*lmMst1 + lmMgl*(308 - 40*lmMst1 - 90*z2) + 113*z2 + 90*
        lmMst1*z2 + 54*pow2(lmMgl) + 2*pow2(lmMst1))*pow2(Mgl) - (215 + 118*
        lmMst1 + 34*z2 + 41*pow2(lmMst1))*pow2(Mst1)) + 240*Dmsqst1*(Mgl*(-3 +
        z2)*pow2(Mst1) + 2*(-3 + lmMgl - lmMst1 + z2)*pow3(Mgl))) - 3*(10*
        Dmsqst1*(shiftst1 - shiftst2)*(3 - 2*(lmMst1 + z2)) + (10*shiftst1 -
        10*shiftst2 + shiftst3)*(1 - 2*(lmMst1 + z2))*pow2(Msq))*pow2(s2t)*
        pow4(Mst1)) + pow2(Dmst12)*(2*pow2(Mgl)*(-450*Dmsqst1*(-4 + z2)*pow2(
        Mst1)*pow2(s2t) + pow2(Msq)*(-4*(163 + 15*lmMgl - 27*lmMst1 + 20*lmMt -
        11*z2 - 8*pow2(lmMst1))*pow2(Mt) + 5*(200 + 18*lmMgl + 265*lmMst1 +
        274*z2 + 91*pow2(lmMst1))*pow2(Mst1)*pow2(s2t))) + Mt*s2t*(40*Mgl*(30*
        Dmsqst1*(7 - 4*z2) + (149 + 41*lmMst1 - 30*z2)*pow2(Msq))*pow2(Mst1) -
        40*(-120*Dmsqst1*(-3 + lmMgl - lmMst1 + z2) + (-498 - 352*lmMst1 +
        lmMgl*(328 - 40*lmMst1 - 90*z2) + 263*z2 + 90*lmMst1*z2 + 54*pow2(
        lmMgl) + 2*pow2(lmMst1))*pow2(Msq))*pow3(Mgl)) + 5*(15*Dmsqst1*(3 + 4*
        lmMst1*(-2 + shiftst1 + shiftst2) - 6*(shiftst1 + shiftst2) + 4*(-3 +
        shiftst1 + shiftst2)*z2) - (32 + 30*(shiftst1 + shiftst2) + 2*lmMst1*(
        115 - 30*(shiftst1 + shiftst2) - 9*shiftst3) + 21*shiftst3 + (82 - 60*(
        shiftst1 + shiftst2) - 18*shiftst3)*z2 + 82*pow2(lmMst1))*pow2(Msq))*
        pow2(s2t)*pow4(Mst1)) + 80*(5 + 18*lmMgl - 74*lmMst1 + 24*lmMt - 52*z2
        - 8*pow2(lmMst1))*pow2(Mgl)*pow2(Msq)*pow2(Mt)*pow4(Mst2)) - 160*Tbeta*
        (Mt*pow2(Msq)*pow2(Mst2)*(3*Dmst12*Mt*s2t*pow2(Mst1)*(-4*Dmst12*Mgl*
        s2t*(60*Dmsqst1*(2*(3 - lmMgl + lmMst1)*pow2(Mgl) + 3*pow2(Mst1)) +
        pow2(Msq)*(-4*(-11 + lmMgl*(154 - 20*lmMst1) - 148*lmMst1 + 27*pow2(
        lmMgl) + pow2(lmMst1))*pow2(Mgl) + (215 + 118*lmMst1 + 41*pow2(lmMst1))
        *pow2(Mst1))) + Mt*(15*Dmsqst1*(Dmst12*(336*pow2(Mgl) + 31*pow2(Mst1))
        + 8*(36*pow2(Mgl) + (5 - 4*lmMst1)*pow2(Mst1))*pow2(Mst2)) + 4*pow2(
        Msq)*(Dmst12*(2*(522 + 107*lmMst1)*pow2(Mgl) + (-95 - 66*lmMst1 + 82*
        pow2(lmMst1))*pow2(Mst1)) + 2*((501 - 18*lmMgl + 504*lmMst1 - 48*lmMt +
        107*pow2(lmMst1))*pow2(Mgl) - 2*(-6 + 79*lmMst1 + 41*pow2(lmMst1))*
        pow2(Mst1))*pow2(Mst2)))) - 9*s2t*shiftst3*pow2(Msq)*pow4(Mst1)*(-8*
        Dmst12*(-3 + 2*lmMst1)*pow2(Mst2)*pow2(Mt) + pow2(Dmst12)*(16*(-2 +
        lmMst1)*pow2(Mt) + (1 - 2*lmMst1)*pow2(Mst1)*pow2(s2t)) + 8*(-1 + 2*
        lmMst1)*pow2(Mt)*pow4(Mst2)) + 90*s2t*(Dmsqst1*(3 - 2*lmMst1) + (1 - 2*
        lmMst1)*pow2(Msq))*pow4(Mst1)*(-8*Dmst12*shiftst2*pow2(Mst2)*pow2(Mt) +
        (-shiftst1 + shiftst2)*pow2(Dmst12)*pow2(Mst1)*pow2(s2t) + 8*(shiftst1
        - shiftst2)*pow2(Mt)*pow4(Mst2)) + 16*Mgl*pow3(Mt)*(30*Dmsqst1*(pow2(
        Dmst12)*(6*(31 + 27*lmMst1 - 24*lmMt - 2*lmMst1*lmMt + lmMgl*(-3 - 2*
        lmMst1 + 2*lmMt) + 2*pow2(lmMst1))*pow2(Mgl) + (-13 - 19*lmMst1 + 19*
        lmMt)*pow2(Mst1)) - 6*Dmst12*((31 + 27*lmMst1 - 24*lmMt - 2*lmMst1*lmMt
        + lmMgl*(-3 - 2*lmMst1 + 2*lmMt) + 2*pow2(lmMst1))*pow2(Mgl) + (3 + 7*
        lmMst1 - 7*lmMt)*pow2(Mst1))*pow2(Mst2) - 12*((7 + 6*lmMst1 - 5*lmMt -
        lmMst1*lmMt + lmMgl*(-1 - lmMst1 + lmMt) + pow2(lmMst1))*pow2(Mgl) + (1
        + 2*lmMst1 - 2*lmMt)*pow2(Mst1))*pow4(Mst2)) - pow2(Msq)*(pow2(Dmst12)*
        (-((-3749 - 2374*lmMst1 + 1770*lmMt - 117*lmMst1*lmMt + lmMgl*(96 -
        237*lmMst1 + 157*lmMt) + 108*pow2(lmMgl) + 41*pow2(lmMst1) + 48*pow2(
        lmMt))*pow2(Mgl)) + (lmMst1*(297 + 2*lmMt) + 15*pow2(lmMst1) + 3*(47 -
        79*lmMt + 8*pow2(lmMt)))*pow2(Mst1)) + Dmst12*pow2(Mst2)*(-((-605 + 6*
        lmMst1*(-77 + lmMt) + 714*lmMt + 45*pow2(lmMst1) + 72*pow2(lmMt))*pow2(
        Mst1)) + 2*pow2(Mgl)*(2096 - 1167*lmMt - 81*(4 + 2*lmMst1 - 3*lmMt)*
        pow2(lmMgl) + 3*(129 + 61*lmMt)*pow2(lmMst1) - 108*pow2(lmMt) - 3*
        lmMst1*(-793 + 67*lmMt + 24*pow2(lmMt)) + 3*lmMgl*(-229 + lmMst1*(63 -
        142*lmMt) + 35*lmMt + 91*pow2(lmMst1) + 24*pow2(lmMt)) - 27*pow3(lmMgl)
        - 84*pow3(lmMst1))) + 2*(-((-233 + 6*lmMst1*(-32 + lmMt) + 348*lmMt +
        45*pow2(lmMst1) + 72*pow2(lmMt))*pow2(Mst1)) + pow2(Mgl)*(758 - 570*
        lmMt - 81*(3 + 2*lmMst1 - 3*lmMt)*pow2(lmMgl) + 3*(141 + 61*lmMt)*pow2(
        lmMst1) - 72*pow2(lmMt) - 6*lmMst1*(-236 + 49*lmMt + 12*pow2(lmMt)) +
        3*lmMgl*(-200 + 2*lmMst1 + 76*lmMt - 142*lmMst1*lmMt + 91*pow2(lmMst1)
        + 24*pow2(lmMt)) - 27*pow3(lmMgl) - 84*pow3(lmMst1)))*pow4(Mst2))) + 2*
        z2*(90*Dmsqst1*(8*pow2(Mt)*(Dmst12*pow2(Mst2)*(4*(7 + 2*lmMst1 - 2*
        lmMt)*Mgl*Mt*pow2(Mst1) - 3*s2t*pow2(Mgl)*pow2(Mst1) + 24*(4 + lmMst1 -
        lmMt)*Mt*pow3(Mgl) + s2t*(-2 + shiftst2)*pow4(Mst1)) + pow2(Mst1)*(4*(4
        + lmMst1 - lmMt)*Mgl*Mt + s2t*(-shiftst1 + shiftst2)*pow2(Mst1))*pow4(
        Mst2)) - pow2(Dmst12)*(8*pow3(Mgl)*(-(Mt*pow2(Mst1)*pow2(s2t)) + 24*(4
        + lmMst1 - lmMt)*pow3(Mt)) + 8*s2t*pow2(Mt)*(6*pow2(Mgl)*pow2(Mst1) +
        pow4(Mst1)) - 4*Mgl*(4*(7 + 2*lmMst1 - 2*lmMt)*pow2(Mst1)*pow3(Mt) +
        Mt*pow2(s2t)*pow4(Mst1)) - (shiftst1 - shiftst2)*pow3(s2t)*pow6(Mst1)))
        + pow2(Msq)*(pow2(Mt)*(-8*Dmst12*pow2(Mst2)*(-2*(457 + 87*lmMst1 - 87*
        lmMt)*Mgl*Mt*pow2(Mst1) - 531*s2t*pow2(Mgl)*pow2(Mst1) + 2*(1627 + 234*
        lmMst1 + 174*lmMt)*Mt*pow3(Mgl) - 9*s2t*(-10 + 10*shiftst2 + shiftst3)*
        pow4(Mst1)) + 8*pow2(Mst1)*(4*(100 - 3*lmMst1 + 3*lmMt)*Mgl*Mt - 9*s2t*
        (10*shiftst1 - 10*shiftst2 + shiftst3)*pow2(Mst1))*pow4(Mst2)) - pow2(
        Dmst12)*(4*Mgl*Mt*pow2(Mst1)*(270*Mgl*Mt*s2t - 2*(691 + 216*lmMst1 -
        216*lmMt)*pow2(Mt) + 51*pow2(Mst1)*pow2(s2t)) - 12*Mt*((2098 + 548*
        lmMst1 - 548*lmMt)*pow2(Mt) + (113 - 90*lmMgl + 90*lmMst1)*pow2(Mst1)*
        pow2(s2t))*pow3(Mgl) + 72*s2t*(5 + shiftst3)*pow2(Mt)*pow4(Mst1) - 9*(
        10*shiftst1 - 10*shiftst2 + shiftst3)*pow3(s2t)*pow6(Mst1))))) + 16*
        pow3(Mgl)*pow3(Mt)*(45*(29 - 10*lmMgl + 10*lmMst1)*Mgl*s2t*xDmsqst1*
        xDmst12*xMgl*pow2(Dmsqst1)*pow3(Dmst12) + Mt*z2*(720*Dmsqst1*(4 +
        lmMst1 - lmMt)*pow2(Msq) - 4*(1070 + 207*lmMst1 - 3*lmMt)*pow4(Msq))*
        pow6(Mst2))))))/(pow4(Msq)*pow6(Mst1)) + 1080*pow2(Mst2)*pow2(Mt)*(3*
        Dmst12*oneLoopFlag*s2t*(Dmst12*MuSUSY*s2t + 6*Dmst12*Mt*Tbeta - 12*Mt*
        Tbeta*pow2(Mst2)) - (16*Al4p*twoLoopFlag*(pow4(Mst1)*(Dmst12*s2t*(4*
        Mgl*Mt*MuSUSY*((1 - 2*lmMgl + 2*lmMst1)*pow2(Mgl) + (2 + lmMst1)*pow2(
        Mst1))*pow2(Mst2) + Dmst12*((3 + 2*lmMst1)*MuSUSY*s2t*pow2(Mgl)*pow2(
        Mst1) + (2*(3 - 2*lmMgl + 2*lmMst1)*Mt*MuSUSY + 3*(1 - 2*lmMgl + 2*
        lmMst1)*s2t*Tbeta*pow2(Mst1))*pow3(Mgl) - (1 + 2*lmMst1)*MuSUSY*s2t*
        pow4(Mst1) + Mgl*(2*Mt*MuSUSY*pow2(Mst1) + 3*(2 + lmMst1)*s2t*Tbeta*
        pow4(Mst1)))) + Mt*Tbeta*(-3*Dmst12*s2t*pow2(Mst1)*(2*Dmst12*(pow2(Mgl)
        + (-1 + 2*lmMst1)*pow2(Mst1)) + 4*((2 + lmMst1)*pow2(Mgl) - (1 + 2*
        lmMst1)*pow2(Mst1))*pow2(Mst2)) + Mt*(pow2(Dmst12)*(4*(1 + lmMt)*Mgl*
        pow2(Mst1) - 8*(-3 + lmMgl - 2*lmMst1 + lmMt)*pow3(Mgl)) + 3*Dmst12*
        pow2(Mst2)*(-4*(1 + lmMt)*Mgl*pow2(Mst1) + 4*(2 + 6*lmMst1 - 3*lmMt -
        2*lmMst1*lmMt + lmMgl*(-3 - 2*lmMst1 + 2*lmMt) + 2*pow2(lmMst1))*pow3(
        Mgl)) + (-24*(1 + lmMt)*Mgl*pow2(Mst1) - 24*(lmMgl*(1 + lmMst1 - lmMt)
        + lmMst1*(-2 + lmMt) + lmMt - pow2(lmMst1))*pow3(Mgl))*pow4(Mst2)))) -
        x*(Dmst12*(-1 + 3*lmMgl - 3*lmMst1)*s2t*pow2(Mst1)*(-4*Dmst12*Mt*MuSUSY
        + 3*Dmst12*s2t*Tbeta*pow2(Mst1) + 4*Mt*MuSUSY*pow2(Mst2)) + 3*Mt*Tbeta*
        (4*Mt*(9 + lmMgl*(20 + 6*lmMst1 - 6*lmMt) + 2*lmMt + lmMst1*(-22 + 6*
        lmMt) - 6*pow2(lmMst1))*pow2(Mgl)*pow4(Mst2) + pow2(Mst1)*(4*Dmst12*(-1
        + 3*lmMgl - 3*lmMst1)*Mgl*s2t*(Dmst12 - pow2(Mst2)) + 4*Dmst12*Mt*(
        lmMst1*(22 - 6*lmMt) - 5*lmMt + lmMgl*(-17 - 6*lmMst1 + 6*lmMt) + 6*
        pow2(lmMst1))*(Dmst12 - pow2(Mst2)) + 8*Mt*(2 + lmMgl*(5 + 2*lmMst1 -
        2*lmMt) + 2*lmMst1*(-3 + lmMt) + lmMt - 2*pow2(lmMst1))*pow4(Mst2))))*
        pow5(Mgl)))/pow8(Mst1)))/Tbeta + (8*Al4p*((5*xMgl*pow2(Mst1)*pow4(Mgl)*
        (216*Dmst12*Mt*s2t*twoLoopFlag*pow2(Mst1)*(-2*Dmst12*Mt*((1 - 3*lmMgl +
        3*lmMst1)*MuSUSY*s2t + 6*(-3 + 2*lmMgl - 2*lmMst1)*Mt*Tbeta)*pow2(Mst2)
        + xDmst12*pow2(Dmst12)*(4*(1 - 3*lmMgl + 3*lmMst1)*Mt*MuSUSY*s2t + 96*(
        -1 + lmMgl - lmMst1)*Tbeta*pow2(Mt) + (-1 + 3*lmMgl - 3*lmMst1)*Tbeta*
        pow2(Mst1)*pow2(s2t)) + 24*(1 - 2*lmMgl + 2*lmMst1)*Tbeta*pow2(Mt)*
        pow4(Mst2)) + (Al4p*threeLoopFlag*(xDmst12*pow2(Msq)*pow3(Dmst12)*(
        25920*Dmsqst1*(29 - 10*lmMgl + 10*lmMst1)*s2t*Tbeta*pow2(Mst1)*pow3(Mt)
        + Mt*pow2(Msq)*(1728*s2t*Tbeta*(-1797 - 1711*lmMst1 + lmMgl*(1315 -
        198*lmMst1 - 96*lmMt) + 144*lmMt + 96*lmMst1*lmMt + 333*pow2(lmMgl) -
        103*pow2(lmMst1))*pow2(Mst1)*pow2(Mt) + 4*Mt*MuSUSY*(41653 + 11784*
        lmMst1 + 48*lmMgl*(-313 + 180*lmMst1) - 12636*pow2(lmMgl) + 540*pow2(
        lmMst1))*pow2(Mst1)*pow2(s2t) + 1152*MuSUSY*(262 + 335*lmMst1 - 96*lmMt
        - 72*lmMst1*lmMt + lmMgl*(-199 - 120*lmMst1 + 72*lmMt) + 120*pow2(
        lmMst1))*pow3(Mt) + Tbeta*(-41653 - 11784*lmMst1 - 48*lmMgl*(-313 +
        180*lmMst1) + 12636*pow2(lmMgl) - 540*pow2(lmMst1))*pow3(s2t)*pow4(
        Mst1))) + pow2(Mst2)*(-432*Dmst12*s2t*Tbeta*pow2(Mst1)*pow3(Mt)*(60*
        Dmsqst1*(29 - 10*lmMgl + 10*lmMst1)*(Dmsqst1*xDmsqst1 + pow2(Msq))*(
        Dmst12 - pow2(Mst2)) - (Dmst12*(2961 + lmMst1*(1964 - 96*lmMt) - 192*
        lmMt + 2*lmMgl*(-683 + 99*lmMst1 + 48*lmMt) - 333*pow2(lmMgl) + 103*
        pow2(lmMst1)) + 2*(633 + 1458*lmMst1 - 96*(1 + lmMst1)*lmMt + 2*lmMgl*(
        -632 + 99*lmMst1 + 48*lmMt) - 333*pow2(lmMgl) + 103*pow2(lmMst1))*pow2(
        Mst2))*pow4(Msq)) + 2*MuSUSY*pow4(Msq)*(-(pow2(Dmst12)*pow2(Mt)*(576*(
        262 + 335*lmMst1 - 96*lmMt - 72*lmMst1*lmMt + lmMgl*(-199 - 120*lmMst1
        + 72*lmMt) + 120*pow2(lmMst1))*pow2(Mt) + (41653 + 11784*lmMst1 + 48*
        lmMgl*(-313 + 180*lmMst1) - 12636*pow2(lmMgl) + 540*pow2(lmMst1))*pow2(
        Mst1)*pow2(s2t))) + 576*pow2(Mst2)*(Dmst12*(262 + 335*lmMst1 - 96*lmMt
        - 72*lmMst1*lmMt + lmMgl*(-199 - 120*lmMst1 + 72*lmMt) + 120*pow2(
        lmMst1)) + 2*(59 + 85*lmMst1 - 24*lmMt - 24*lmMst1*lmMt + lmMgl*(-53 -
        40*lmMst1 + 24*lmMt) + 40*pow2(lmMst1))*pow2(Mst2))*pow4(Mt)))))/pow4(
        Msq)))/Tbeta + 72*z2*((Al4p*threeLoopFlag*pow4(Mst1)*(-(Mt*xDmst12*
        pow2(Msq)*pow3(Dmst12)*(Mt*(8*Mgl*pow2(Mst1)*(75*Dmsqst1*(20*Mt*MuSUSY*
        s2t + 48*(7 + 2*lmMst1 - 2*lmMt)*Tbeta*pow2(Mt) - 3*Tbeta*pow2(Mst1)*
        pow2(s2t)) - pow2(Msq)*(237*Mt*MuSUSY*s2t + 10*(16 - 9*lmMst1 + 9*lmMt)
        *Tbeta*pow2(Mt) + 240*Tbeta*pow2(Mst1)*pow2(s2t))) - 10*(240*Dmsqst1*(
        2*Mt*MuSUSY*s2t + 36*(4 + lmMst1 - lmMt)*Tbeta*pow2(Mt) - 3*Tbeta*pow2(
        Mst1)*pow2(s2t)) - pow2(Msq)*(32*Mt*((94 - 45*lmMgl + 45*lmMst1)*
        MuSUSY*s2t + (380 + 147*lmMst1 - 249*lmMt)*Mt*Tbeta) - 3*(37 + 90*lmMgl
        - 90*lmMst1)*Tbeta*pow2(Mst1)*pow2(s2t)))*pow3(Mgl)) - s2t*(30*Dmsqst1*
        s2t*(Mt*MuSUSY*(6*shiftst1 + 4*shiftst2) + 5*s2t*(3 - 4*shiftst1 + 2*
        shiftst2)*Tbeta*pow2(Mst1)) + pow2(Msq)*(4*Mt*MuSUSY*s2t*(10 + 45*
        shiftst1 + 30*shiftst2 + 27*shiftst3) + 40*(-16 + 9*shiftst3)*Tbeta*
        pow2(Mt) + 5*(41 - 120*shiftst1 + 60*shiftst2 - 21*shiftst3)*Tbeta*
        pow2(Mst1)*pow2(s2t)))*pow4(Mst1) + 2*pow2(Mgl)*(-225*Dmsqst1*s2t*pow2(
        Mst1)*(-4*Mt*MuSUSY*s2t + 120*Tbeta*pow2(Mt) + Tbeta*pow2(Mst1)*pow2(
        s2t)) + pow2(Msq)*(4556*s2t*Tbeta*pow2(Mst1)*pow2(Mt) - 910*Mt*MuSUSY*
        pow2(Mst1)*pow2(s2t) + 1992*MuSUSY*pow3(Mt) + 685*Tbeta*pow3(s2t)*pow4(
        Mst1))))) + 30*xDmsqst1*pow2(Dmsqst1)*pow2(Mst2)*(5*Dmst12*MuSUSY*s2t*
        pow2(Mt)*(Dmst12*(8*Mgl*Mt*pow2(Mst1) + 6*s2t*pow2(Mgl)*pow2(Mst1) -
        32*Mt*pow3(Mgl) + s2t*(9 - 2*shiftst1 - 2*shiftst2)*pow4(Mst1)) + 4*
        pow2(Mst2)*(10*Mgl*Mt*pow2(Mst1) + 8*Mt*pow3(Mgl) + s2t*(shiftst1 -
        shiftst2)*pow4(Mst1))) + 15*Tbeta*(pow3(Mt)*(4*Dmst12*pow2(Mst2)*(2*(45
        + 14*lmMst1 - 14*lmMt)*Mgl*Mt*pow2(Mst1) - 15*s2t*pow2(Mgl)*pow2(Mst1)
        + 48*(4 + lmMst1 - lmMt)*Mt*pow3(Mgl) + 2*s2t*(-3 + shiftst2)*pow4(
        Mst1)) + 8*(2*(13 + 4*lmMst1 - 4*lmMt)*Mgl*Mt*pow2(Mst1) + 4*(19 + 5*
        lmMst1 - 5*lmMt)*Mt*pow3(Mgl) + s2t*(-shiftst1 + shiftst2)*pow4(Mst1))*
        pow4(Mst2)) - Mt*pow2(Dmst12)*(4*Mgl*((6 + 4*lmMst1 - 4*lmMt)*Mt + 3*
        Mgl*s2t)*pow2(Mst1)*pow2(Mt) + 8*pow3(Mgl)*(-(Mt*pow2(Mst1)*pow2(s2t))
        + 24*(4 + lmMst1 - lmMt)*pow3(Mt)) + 2*Mt*s2t*(6*Mt - 5*Mgl*s2t)*pow4(
        Mst1) - (shiftst1 - shiftst2)*pow3(s2t)*pow6(Mst1))))))/(Tbeta*pow4(
        Msq)) - 1440*twoLoopFlag*pow4(Mt)*((Dmst12 + pow2(Mst2))*pow3(Mgl)*
        pow4(Mst1)*pow4(Mst2) + x*pow2(Mst1)*pow2(Mst2)*(-3*pow2(Dmst12) + 3*
        Dmst12*pow2(Mst2) + 2*pow4(Mst2))*pow5(Mgl) - xDmst12*pow2(Mst1)*pow3(
        Dmst12)*(pow2(Mst1)*pow3(Mgl) - 3*x*pow5(Mgl)) + 3*x*pow6(Mst2)*pow7(
        Mgl)))))/pow8(Mst1)))/(155520.*pow6(Mst2));
}

/**
 * 	@return returns the susy log^0 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double H9q2::calc_coef_at_as2_no_sm_logs_log0() const {

   const double result =
      (-432*pow2(log(pow2(Mst1)/pow2(Mgl)))*pow4(Msq)*(-2*Mt*pow2(Dmst12)*pow2(
        Mst2)*(-72*Mt*s2t*pow2(Mst1) + 202*Mgl*pow2(Mt) + 111*Mgl*pow2(Mst1)*
        pow2(s2t))*pow3(Mgl) + pow3(Dmst12)*pow3(Mgl)*(-48*s2t*pow2(Mst1)*pow2(
        Mt) + 111*Mgl*Mt*pow2(Mst1)*pow2(s2t) + 1568*Mgl*pow3(Mt) + 24*pow3(
        s2t)*pow4(Mst1)) - 8*Dmst12*pow2(Mgl)*(95*Mt*pow2(Mgl) - 18*Mt*pow2(
        Mst1) - 81*Mgl*s2t*pow2(Mst1))*pow2(Mt)*pow4(Mst2) - 8*pow3(Mt)*(-36*
        pow2(Mgl)*pow2(Mst1) + 77*pow4(Mgl) + 4*pow4(Mst1))*pow6(Mst2)) - 192*
        pow2(Mt)*pow3(log(pow2(Mst1)/pow2(Mgl)))*pow4(Msq)*(-95*Mt*pow2(Dmst12)
        *pow2(Mst2)*pow4(Mgl) + 504*Mt*pow3(Dmst12)*pow4(Mgl) + 2*Dmst12*(-157*
        Mgl*Mt + 54*s2t*pow2(Mst1))*pow3(Mgl)*pow4(Mst2) + 2*Mt*(-157*pow4(Mgl)
        + 18*pow4(Mst1))*pow6(Mst2)) + 5*pow2(Dmsqst1)*(216*Dmst12*pow2(Mt)*(
        24*Mt*(-11 + 90*z2 - 112*z3)*pow2(Mgl)*pow2(Mst1) - 64*s2t*(-27 + 43*z2
        - 40*z3)*pow2(Mst1)*pow3(Mgl) + 8*Mt*(-307 + 480*z2 - 480*z3)*pow4(Mgl)
        + Mt*(146 + 144*z2 - 321*z3)*pow4(Mst1) - 32*Mgl*s2t*(-9 + 30*z2 - 32*
        z3)*pow4(Mst1))*pow4(Mst2) - 3*Mt*pow2(Dmst12)*pow2(Mst2)*(4608*Mt*s2t*
        (11*z2 - 8*(2 + z3))*pow2(Mst1)*pow3(Mgl) + 576*((-307 + 480*z2 - 480*
        z3)*pow2(Mt) + (29 - 10*z2)*pow2(Mst1)*pow2(s2t))*pow4(Mgl) + 576*Mgl*
        Mt*s2t*(-33 + 88*z2 - 96*z3)*pow4(Mst1) + (-850 - 8640*z2 + 12123*z3)*
        pow2(Mt)*pow4(Mst1) - 1728*pow2(Mgl)*(4*(-5 + 9*z2 - 8*z3)*pow2(Mst1)*
        pow2(Mt) + 5*(-2 + z2)*pow2(s2t)*pow4(Mst1)) - 54*(-58 + 64*z2 - 35*z3)
        *pow2(s2t)*pow6(Mst1)) + pow3(Dmst12)*(1728*(2*Mt*(29 - 10*z2)*pow2(
        Mst1)*pow2(s2t) + (-307 + 480*z2 - 480*z3)*pow3(Mt))*pow4(Mgl) - 4*(
        5314 + 6336*z2 - 12123*z3)*pow3(Mt)*pow4(Mst1) - 5184*pow2(Mgl)*((-51 +
        162*z2 - 176*z3)*pow2(Mst1)*pow3(Mt) + Mt*(-7 + 4*z2)*pow2(s2t)*pow4(
        Mst1)) - 2304*pow3(Mgl)*(-6*s2t*(-59 + 65*z2 - 56*z3)*pow2(Mst1)*pow2(
        Mt) + (-3 + z2)*pow3(s2t)*pow4(Mst1)) + 27*Mt*(134 - 192*z2 + 105*z3)*
        pow2(s2t)*pow6(Mst1) - 576*Mgl*(12*s2t*(5 - 26*z2 + 32*z3)*pow2(Mt)*
        pow4(Mst1) + 5*(-2 + z2)*pow3(s2t)*pow6(Mst1))) + 48*pow3(Mt)*(432*(-1
        + 13*z2 - 16*z3)*pow2(Mgl)*pow2(Mst1) + 36*(-179 + 380*z2 - 400*z3)*
        pow4(Mgl) + (1118 + 432*z2 - 1539*z3)*pow4(Mst1))*pow6(Mst2)) + 10*
        Dmsqst1*pow2(Msq)*(216*Dmst12*pow2(Mt)*(48*Mt*(-3 + 14*z2 - 16*z3)*
        pow2(Mgl)*pow2(Mst1) + 32*s2t*(19 - 18*z2 + 16*z3)*pow2(Mst1)*pow3(Mgl)
        + 4*Mt*(-307 + 480*z2 - 480*z3)*pow4(Mgl) + Mt*(38 + 48*z2 - 107*z3)*
        pow4(Mst1) - 32*Mgl*s2t*(-4 + 9*z2 - 8*z3)*pow4(Mst1))*pow4(Mst2) -
        108*Mt*pow2(Dmst12)*pow2(Mst2)*(256*Mt*s2t*(-6 + 9*z2 - 8*z3)*pow2(
        Mst1)*pow3(Mgl) + 8*((-307 + 480*z2 - 480*z3)*pow2(Mt) + (29 - 10*z2)*
        pow2(Mst1)*pow2(s2t))*pow4(Mgl) + 32*Mgl*Mt*s2t*(-5 + 14*z2 - 16*z3)*
        pow4(Mst1) + (-6 - 80*z2 + 107*z3)*pow2(Mt)*pow4(Mst1) - 24*pow2(Mgl)*(
        (-19 + 70*z2 - 80*z3)*pow2(Mst1)*pow2(Mt) + 2*(-3 + z2)*pow2(s2t)*pow4(
        Mst1)) + 2*(10 - 16*z2 + 7*z3)*pow2(s2t)*pow6(Mst1)) + pow3(Dmst12)*(
        864*(2*Mt*(29 - 10*z2)*pow2(Mst1)*pow2(s2t) + (-307 + 480*z2 - 480*z3)*
        pow3(Mt))*pow4(Mgl) + (-1814 + 1152*z2 - 567*z3)*pow3(Mt)*pow4(Mst1) +
        2592*pow2(Mgl)*(2*(25 - 98*z2 + 112*z3)*pow2(Mst1)*pow3(Mt) + Mt*(-1 +
        2*z2)*pow2(s2t)*pow4(Mst1)) - 1152*pow3(Mgl)*(-6*s2t*(-67 + 90*z2 - 80*
        z3)*pow2(Mst1)*pow2(Mt) + (-3 + z2)*pow3(s2t)*pow4(Mst1)) + 27*Mt*(18 -
        64*z2 + 7*z3)*pow2(s2t)*pow6(Mst1) - 576*Mgl*(3*s2t*(-5 + 4*z2)*pow2(
        Mt)*pow4(Mst1) + (-3 + z2)*pow3(s2t)*pow6(Mst1))) + 192*pow3(Mt)*(108*(
        -1 + 4*z2 - 4*z3)*pow2(Mgl)*pow2(Mst1) + 9*(-69 + 80*z2 - 80*z3)*pow4(
        Mgl) + 4*(23 + 9*z2 - 36*z3)*pow4(Mst1))*pow6(Mst2)) + 8*pow4(Msq)*(-
        24*Dmst12*pow2(Mt)*(Mt*(3049 - 2716*z2 + 2388*z3)*pow2(Mgl)*pow2(Mst1)
        - 4*s2t*(1043 + 2143*z2 - 570*z3)*pow2(Mst1)*pow3(Mgl) + 2*Mt*(6742 +
        4265*z2 - 174*z3)*pow4(Mgl) - 4*Mt*(49 + 90*z2 - 171*z3)*pow4(Mst1) +
        4*Mgl*s2t*(-407 + 197*z2 - 24*z3)*pow4(Mst1))*pow4(Mst2) + 3*Mt*pow2(
        Dmst12)*pow2(Mst2)*(-48*Mt*s2t*(-1091 + 402*z2 - 480*z3)*pow2(Mst1)*
        pow3(Mgl) + 4*((-40279 + 17982*z2 - 26652*z3)*pow2(Mt) - 3*(681 + 776*
        z2 - 558*z3)*pow2(Mst1)*pow2(s2t))*pow4(Mgl) - 16*Mgl*Mt*s2t*(-148 +
        607*z2 - 696*z3)*pow4(Mst1) + (262 + 1840*z2 - 3423*z3)*pow2(Mt)*pow4(
        Mst1) - 4*pow2(Mgl)*(2*(1327 - 4904*z2 + 5280*z3)*pow2(Mst1)*pow2(Mt) +
        9*(175 + 118*z2 - 4*z3)*pow2(s2t)*pow4(Mst1)) + 18*(-8 + 40*z2 + 31*z3)
        *pow2(s2t)*pow6(Mst1)) + pow3(Dmst12)*(-6*(3*Mt*(1695 - 968*z2 + 558*
        z3)*pow2(Mst1)*pow2(s2t) + 4*(-53763 + 9452*z2 - 26304*z3)*pow3(Mt))*
        pow4(Mgl) + (-4486 + 816*z2 + 3375*z3)*pow3(Mt)*pow4(Mst1) + 12*pow2(
        Mgl)*(4*(-570 + 929*z2 - 840*z3)*pow2(Mst1)*pow3(Mt) + 9*Mt*(1 + 148*z2
        - 4*z3)*pow2(s2t)*pow4(Mst1)) - 24*pow3(Mgl)*(s2t*(-2722 + 4704*z2 -
        3696*z3)*pow2(Mst1)*pow2(Mt) + (-22 + 113*z2 - 243*z3)*pow3(s2t)*pow4(
        Mst1)) - 18*Mt*(-119 + 60*z2 + 192*z3)*pow2(s2t)*pow6(Mst1) + 12*Mgl*(
        s2t*(446 - 768*z2 + 624*z3)*pow2(Mt)*pow4(Mst1) + (215 + 34*z2 - 21*z3)
        *pow3(s2t)*pow6(Mst1))) - 12*pow3(Mt)*(4*(1909 - 496*z2 + 228*z3)*pow2(
        Mgl)*pow2(Mst1) + 2*(6677 + 10414*z2 - 3948*z3)*pow4(Mgl) + (-4125 -
        892*z2 + 4740*z3)*pow4(Mst1))*pow6(Mst2)) + 96*log(pow2(Mst1)/pow2(Mgl)
        )*(60*pow2(Dmsqst1)*pow3(Mgl)*(-15*Mt*pow2(Dmst12)*pow2(Mst2)*(2*Mt*
        s2t*pow2(Mst1) + 3*Mgl*pow2(Mt) - Mgl*pow2(Mst1)*pow2(s2t)) + pow3(
        Dmst12)*(78*s2t*pow2(Mst1)*pow2(Mt) - 30*Mgl*Mt*pow2(Mst1)*pow2(s2t) +
        45*Mgl*pow3(Mt) - 2*pow3(s2t)*pow4(Mst1)) + 9*Dmst12*(5*Mgl*Mt - 2*s2t*
        pow2(Mst1))*pow2(Mt)*pow4(Mst2) + 15*Mgl*pow3(Mt)*pow6(Mst2)) + 60*
        Dmsqst1*pow2(Msq)*pow3(Mgl)*(-3*Mt*pow2(Dmst12)*pow2(Mst2)*(4*Mt*s2t*
        pow2(Mst1) + 15*Mgl*pow2(Mt) - 5*Mgl*pow2(Mst1)*pow2(s2t)) + pow3(
        Dmst12)*(60*s2t*pow2(Mst1)*pow2(Mt) - 30*Mgl*Mt*pow2(Mst1)*pow2(s2t) +
        45*Mgl*pow3(Mt) - 2*pow3(s2t)*pow4(Mst1)) + 9*Dmst12*(5*Mgl*Mt - 4*s2t*
        pow2(Mst1))*pow2(Mt)*pow4(Mst2) + 30*Mgl*pow3(Mt)*pow6(Mst2)) + pow4(
        Msq)*(2*Mt*pow2(Dmst12)*pow2(Mgl)*pow2(Mst2)*(-102*Mgl*Mt*s2t*pow2(
        Mst1) + 2*pow2(Mgl)*((737 + 162*z2)*pow2(Mt) + 3*(328 + 27*z2)*pow2(
        Mst1)*pow2(s2t)) + 27*pow2(s2t)*pow4(Mst1)) - pow2(Mgl)*pow3(Dmst12)*(
        Mt*pow2(Mgl)*(32*(434 + 81*z2)*pow2(Mt) + 3*(605 + 54*z2)*pow2(Mst1)*
        pow2(s2t)) + 54*Mt*pow2(s2t)*pow4(Mst1) - 2*Mgl*(67*s2t*pow2(Mst1)*
        pow2(Mt) + 2*(-154 + 45*z2)*pow3(s2t)*pow4(Mst1))) + 12*Dmst12*pow2(
        Mgl)*(18*Mt*(37 + 6*z2)*pow2(Mgl) + 9*Mt*pow2(Mst1) - 476*Mgl*s2t*pow2(
        Mst1))*pow2(Mt)*pow4(Mst2) + 12*pow3(Mt)*(18*pow2(Mgl)*pow2(Mst1) + 3*(
        215 + 36*z2)*pow4(Mgl) - 52*pow4(Mst1))*pow6(Mst2))))/(108.*pow3(Mt)*
        pow4(Msq)*pow4(Mst1)*pow6(Mst2));

   return result;
}

/**
 * 	@return returns the susy log^1 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double H9q2::calc_coef_at_as2_no_sm_logs_log1() const {

   const double result =
      (8*(-6*pow2(Mt)*pow2(log(pow2(Mst1)/pow2(Mgl)))*pow4(Msq)*(-203*Mt*pow2(
        Dmst12)*pow2(Mst2)*pow4(Mgl) + 936*Mt*pow3(Dmst12)*pow4(Mgl) + 2*
        Dmst12*(-265*Mgl*Mt + 162*s2t*pow2(Mst1))*pow3(Mgl)*pow4(Mst2) + 2*Mt*(
        -265*pow4(Mgl) + 18*pow4(Mst1))*pow6(Mst2)) + 10*Dmsqst1*Mt*pow2(Msq)*(
        18*Dmst12*Mt*(6*Mt*(-7 + 4*z2)*pow2(Mgl)*pow2(Mst1) - 8*s2t*(-5 + 2*z2)
        *pow2(Mst1)*pow3(Mgl) + Mt*(-119 + 60*z2)*pow4(Mgl) + 4*Mt*(-2 + z2)*
        pow4(Mst1) - 8*Mgl*s2t*(-2 + z2)*pow4(Mst1))*pow4(Mst2) + 2*pow3(
        Dmst12)*(-27*(-47 + 28*z2)*pow2(Mgl)*pow2(Mst1)*pow2(Mt) + 72*Mt*s2t*(-
        19 + 10*z2)*pow2(Mst1)*pow3(Mgl) + 9*(-119 + 60*z2)*pow2(Mt)*pow4(Mgl)
        + 6*Mgl*Mt*s2t*pow4(Mst1) + pow2(Mt)*pow4(Mst1) - 9*pow2(s2t)*pow6(
        Mst1)) - 9*pow2(Dmst12)*pow2(Mst2)*(3*(33 - 20*z2)*pow2(Mgl)*pow2(Mst1)
        *pow2(Mt) + 16*Mt*s2t*(-7 + 4*z2)*pow2(Mst1)*pow3(Mgl) + 2*(-119 + 60*
        z2)*pow2(Mt)*pow4(Mgl) + 8*Mgl*Mt*s2t*(-3 + 2*z2)*pow4(Mst1) + (7 - 4*
        z2)*pow2(Mt)*pow4(Mst1) - 2*pow2(s2t)*pow6(Mst1)) + 6*pow2(Mt)*(36*(-2
        + z2)*pow2(Mgl)*pow2(Mst1) + 3*(-49 + 20*z2)*pow4(Mgl) + 2*(-17 + 6*z2)
        *pow4(Mst1))*pow6(Mst2)) + 5*Mt*pow2(Dmsqst1)*(2*pow3(Dmst12)*(-27*(-75
        + 44*z2)*pow2(Mgl)*pow2(Mst1)*pow2(Mt) + 144*Mt*s2t*(-15 + 7*z2)*pow2(
        Mst1)*pow3(Mgl) + 18*(-119 + 60*z2)*pow2(Mt)*pow4(Mgl) + 24*Mgl*Mt*s2t*
        (-19 + 12*z2)*pow4(Mst1) - 2*((-61 + 36*z2)*pow2(Mt) + 9*pow2(Mst1)*
        pow2(s2t))*pow4(Mst1)) + 18*Dmst12*Mt*(3*Mt*(-47 + 28*z2)*pow2(Mgl)*
        pow2(Mst1) - 16*s2t*(-9 + 5*z2)*pow2(Mst1)*pow3(Mgl) + 2*Mt*(-119 + 60*
        z2)*pow4(Mgl) - 8*Mgl*s2t*(-7 + 4*z2)*pow4(Mst1) + Mt*(-23 + 12*z2)*
        pow4(Mst1))*pow4(Mst2) - 3*pow2(Dmst12)*pow2(Mst2)*(-36*(-7 + 4*z2)*
        pow2(Mgl)*pow2(Mst1)*pow2(Mt) + 96*Mt*s2t*(-3 + z2)*pow2(Mst1)*pow3(
        Mgl) + 12*(-119 + 60*z2)*pow2(Mt)*pow4(Mgl) + 12*Mgl*Mt*s2t*(-19 + 12*
        z2)*pow4(Mst1) + (61 - 36*z2)*pow2(Mt)*pow4(Mst1) - 12*pow2(s2t)*pow6(
        Mst1)) + 6*pow2(Mt)*(36*(-7 + 4*z2)*pow2(Mgl)*pow2(Mst1) + (-537 + 300*
        z2)*pow4(Mgl) + 4*(-22 + 9*z2)*pow4(Mst1))*pow6(Mst2)) + pow4(Msq)*(Mt*
        pow2(Dmst12)*pow2(Mst2)*(-6*pow2(Mgl)*pow2(Mst1)*((659 - 440*z2)*pow2(
        Mt) + 243*pow2(Mst1)*pow2(s2t)) - 72*Mt*s2t*(-99 + 20*z2)*pow2(Mst1)*
        pow3(Mgl) + ((-23651 + 5436*z2)*pow2(Mt) - 582*pow2(Mst1)*pow2(s2t))*
        pow4(Mgl) + 12*Mgl*Mt*s2t*(25 - 62*z2)*pow4(Mst1) + 6*((-145 + 46*z2)*
        pow2(Mt) + 79*pow2(Mst1)*pow2(s2t))*pow4(Mst1)) + pow3(Dmst12)*(3*Mt*
        pow2(Mgl)*pow2(Mst1)*((-498 + 280*z2)*pow2(Mt) + 379*pow2(Mst1)*pow2(
        s2t)) - 24*s2t*pow2(Mst1)*((-86 + 77*z2)*pow2(Mt) + pow2(Mst1)*pow2(
        s2t))*pow3(Mgl) + (-315*Mt*pow2(Mst1)*pow2(s2t) + (63184 - 7104*z2)*
        pow3(Mt))*pow4(Mgl) + 2*Mgl*s2t*((442 - 60*z2)*pow2(Mt) + 59*pow2(Mst1)
        *pow2(s2t))*pow4(Mst1) - 3*Mt*(4*(-43 + 8*z2)*pow2(Mt) + 125*pow2(Mst1)
        *pow2(s2t))*pow4(Mst1)) - 6*Dmst12*pow2(Mt)*(4*Mt*(103 - 50*z2)*pow2(
        Mgl)*pow2(Mst1) - 4*s2t*(283 + 138*z2)*pow2(Mst1)*pow3(Mgl) + Mt*(2647
        + 628*z2)*pow4(Mgl) + 12*Mt*(8 - 5*z2)*pow4(Mst1) - 4*Mgl*s2t*(65 + 2*
        z2)*pow4(Mst1))*pow4(Mst2) - 12*pow3(Mt)*(2*(91 - 10*z2)*pow2(Mgl)*
        pow2(Mst1) + (589 + 464*z2)*pow4(Mgl) - 4*(-35 + 15*z2 + 3*z3)*pow4(
        Mst1))*pow6(Mst2)) + 2*log(pow2(Mst1)/pow2(Mgl))*(180*Dmsqst1*pow2(Msq)
        *pow2(Mt)*pow3(Mgl)*(-5*Mgl*Mt*pow2(Dmst12)*pow2(Mst2) + (5*Mgl*Mt + 4*
        s2t*pow2(Mst1))*pow3(Dmst12) + Dmst12*(5*Mgl*Mt - 4*s2t*pow2(Mst1))*
        pow4(Mst2) + 5*Mgl*Mt*pow6(Mst2)) + 90*pow2(Dmsqst1)*pow2(Mt)*pow3(Mgl)
        *(-2*pow2(Dmst12)*(5*Mgl*Mt + 2*s2t*pow2(Mst1))*pow2(Mst2) + 2*(5*Mgl*
        Mt + 6*s2t*pow2(Mst1))*pow3(Dmst12) + 2*Dmst12*(5*Mgl*Mt - 2*s2t*pow2(
        Mst1))*pow4(Mst2) + 5*Mgl*Mt*pow6(Mst2)) + pow4(Msq)*(6*Mt*pow2(Dmst12)
        *pow2(Mst2)*(14*Mt*s2t*pow2(Mst1) + 67*Mgl*pow2(Mt) + 117*Mgl*pow2(
        Mst1)*pow2(s2t))*pow3(Mgl) - pow3(Dmst12)*pow3(Mgl)*(42*s2t*pow2(Mst1)*
        pow2(Mt) + 351*Mgl*(Mt*pow2(Mst1)*pow2(s2t) + 8*pow3(Mt)) + 68*pow3(
        s2t)*pow4(Mst1)) + 12*Dmst12*pow2(Mgl)*(167*Mt*pow2(Mgl) - 54*Mt*pow2(
        Mst1) - 170*Mgl*s2t*pow2(Mst1))*pow2(Mt)*pow4(Mst2) + 36*pow3(Mt)*(-36*
        pow2(Mgl)*pow2(Mst1) + 59*pow4(Mgl) + 4*pow4(Mst1))*pow6(Mst2)))))/(9.*
        pow3(Mt)*pow4(Msq)*pow4(Mst1)*pow6(Mst2));

   return result;
}

/**
 * 	@return returns the susy log^2 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double H9q2::calc_coef_at_as2_no_sm_logs_log2() const {

   const double result =
      (8*(Mt*pow2(Dmst12)*pow2(Mst2)*pow4(Msq)*(528*Mt*s2t*pow2(Mst1)*pow3(Mgl)
        + (-1966*pow2(Mt) + 96*pow2(Mst1)*pow2(s2t))*pow4(Mgl) + 180*Mgl*Mt*
        s2t*pow4(Mst1) + 3*(51*pow2(Mt) + 82*pow2(Mst1)*pow2(s2t))*pow4(Mst1) +
        pow2(Mgl)*(64*pow2(Mst1)*pow2(Mt) - 321*pow2(s2t)*pow4(Mst1))) - 6*
        Dmst12*pow2(Mt)*pow4(Msq)*(43*Mt*pow2(Mgl)*pow2(Mst1) - 248*s2t*pow2(
        Mst1)*pow3(Mgl) + 550*Mt*pow4(Mgl) + 51*Mt*pow4(Mst1) + 60*Mgl*s2t*
        pow4(Mst1))*pow4(Mst2) - pow3(Dmst12)*pow4(Msq)*((48*Mt*pow2(Mst1)*
        pow2(s2t) - 7232*pow3(Mt))*pow4(Mgl) + 102*pow3(Mt)*pow4(Mst1) + pow2(
        Mgl)*(32*pow2(Mst1)*pow3(Mt) - 321*Mt*pow2(s2t)*pow4(Mst1)) + 16*pow3(
        Mgl)*(11*s2t*pow2(Mst1)*pow2(Mt) + 2*pow3(s2t)*pow4(Mst1)) + 369*Mt*
        pow2(s2t)*pow6(Mst1) + Mgl*(120*s2t*pow2(Mt)*pow4(Mst1) - 41*pow3(s2t)*
        pow6(Mst1))) + 6*pow3(Mt)*(-86*pow2(Mgl)*pow2(Mst1)*pow4(Msq) - 380*
        pow4(Mgl)*pow4(Msq) + (15*pow2(Dmsqst1) + 30*Dmsqst1*pow2(Msq) + 98*
        pow4(Msq))*pow4(Mst1))*pow6(Mst2) - 12*log(pow2(Mst1)/pow2(Mgl))*pow2(
        Mt)*pow4(Msq)*(-53*Mt*pow2(Dmst12)*pow2(Mst2)*pow4(Mgl) + 276*Mt*pow3(
        Dmst12)*pow4(Mgl) + 2*Dmst12*(-85*Mgl*Mt + 44*s2t*pow2(Mst1))*pow3(Mgl)
        *pow4(Mst2) + 2*Mt*(-85*pow4(Mgl) + 9*pow4(Mst1))*pow6(Mst2))))/(9.*
        pow3(Mt)*pow4(Msq)*pow4(Mst1)*pow6(Mst2));

   return result;
}

/**
 * 	@return returns the susy log^3 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double H9q2::calc_coef_at_as2_no_sm_logs_log3() const {

   const double result =
      -298.6666666666667;

   return result;
}

}	// hierarchies
}	// himalaya
