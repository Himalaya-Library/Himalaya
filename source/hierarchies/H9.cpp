// ====================================================================
// This file is part of Himalaya.
//
// Himalaya is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================

#include "H9.hpp"
#include "enums.hpp"
#include "Constants.hpp"
#include "power.hpp"
#include <cmath>

namespace himalaya{
namespace hierarchies{

/**
 * Constructor
 * @param flagMap the flagMap for the truncation of expansion variables
 * @param Al4p a double alpha_s/4/Pi
 * @param beta a double which is the mixing angle beta
 * @param Dmst12 a double Mst1^2 - Mst2^2
 * @param Dmsqst1  a double Msq1^2 - Mst1^2
 * @param lmMt a double log((<renormalization scale> / Mt)^2)
 * @param lmMgl a double log((<renormalization scale> / Mgl)^2)
 * @param lmMst1 a double log((<renormalization scale> / Mst1)^2)
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
H9::H9(const ExpansionFlag_t& flagMap, double Al4p, double beta, double Dmst12, double Dmsqst1,
                 double lmMt, double lmMgl, double lmMst1,
                 double Mgl, double Mt,  double Mst1, double Mst2, double MuSUSY,
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
 * @return The diagonal (1, 1) matrix element of the Higgs mass matrix as a double for the hierarchy 'H9'
 */
double H9::getS1() const {
   return -(pow2(MuSUSY)*((10*Al4p*xMgl*pow4(Mgl)*(2*xDmst12*pow2(Mt)*(216*(1 - 3*
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
        24*lmMt) + 40*pow2(lmMst1))*pow2(Mst2))*pow4(Mt)))))/9. +
        threeLoopFlag*pow2(Al4p)*(-16*z2*pow2(Mt)*(-2*xDmst12*pow3(Dmst12)*(Mt*
        s2t*(80*(15*Dmsqst1 + (-94 + 45*lmMgl - 45*lmMst1)*pow2(Mst1))*pow3(
        Mgl) + 6*Mgl*(150*xDmsqst1*pow2(Dmsqst1) - 500*Dmsqst1*pow2(Mst1) + 79*
        pow4(Mst1))) + pow2(Mst1)*(pow2(Mgl)*(-996*pow2(Mt) + 5*(-90*Dmsqst1 +
        91*pow2(Mst1))*pow2(s2t)) + pow2(s2t)*(-225*xDmsqst1*pow2(Dmsqst1) +
        15*Dmsqst1*(3*shiftst1 + 2*shiftst2)*pow2(Mst1) + (10 + 45*shiftst1 +
        30*shiftst2 + 27*shiftst3)*pow4(Mst1)))) + pow2(Mst2)*(4*pow2(Mgl)*
        pow2(Mst1)*(11*pow2(Dmst12) - 520*pow2(Mst2)*(Dmst12 + pow2(Mst2)))*
        pow2(Mt) + Dmst12*s2t*(225*xDmsqst1*pow2(Dmsqst1)*(8*Dmst12*Mgl*Mt -
        Dmst12*s2t*pow2(Mst1) - 8*Mgl*Mt*pow2(Mst2)) + 5*s2t*pow2(Mst1)*(274*
        Dmst12*pow2(Mgl)*pow2(Mst1) - 30*Dmsqst1*(3*Dmst12*pow2(Mgl) - Dmst12*(
        -3 + shiftst1 + shiftst2)*pow2(Mst1) + 2*(shiftst1 - shiftst2)*pow2(
        Mst1)*pow2(Mst2)) + Dmst12*(-41 + 30*shiftst1 + 30*shiftst2 + 9*
        shiftst3)*pow4(Mst1) - 6*(10*shiftst1 - 10*shiftst2 + shiftst3)*pow2(
        Mst2)*pow4(Mst1)) - 5*Mt*(8*Mgl*pow2(Mst2)*((113 - 90*lmMgl + 90*
        lmMst1)*pow2(Mgl)*pow2(Mst1) + 30*Dmsqst1*(2*pow2(Mgl) + pow2(Mst1)) -
        17*pow4(Mst1)) + 4*Dmst12*((263 - 90*lmMgl + 90*lmMst1)*pow2(Mst1)*
        pow3(Mgl) - 120*Dmsqst1*(-(Mgl*pow2(Mst1)) + pow3(Mgl)) + 30*Mgl*pow4(
        Mst1))))) + 20*xMgl*pow4(Mgl)*(-6*xDmst12*(12*(17 - 4*lmMgl + 4*lmMst1)
        *pow2(Mt) + (49 - 17*lmMgl + 17*lmMst1)*pow2(Mst1)*pow2(s2t))*pow3(
        Dmst12) + pow2(Mst2)*(3*(49 - 17*lmMgl + 17*lmMst1)*pow2(Dmst12)*pow2(
        Mst1)*pow2(s2t) + 8*pow2(Mt)*(9*Dmst12*(17 - 4*lmMgl + 4*lmMst1)*(
        Dmst12 - pow2(Mst2)) + 2*(-47 + 12*lmMgl - 12*lmMst1)*pow4(Mst2))))) +
        pow2(Mst2)*(10*Dmst12*s2t*pow2(Mt)*(-15*xDmsqst1*pow2(Dmsqst1)*(-256*
        Dmst12*Mgl*Mt + Dmst12*s2t*(21 - 8*(shiftst1 + shiftst2))*pow2(Mst1) +
        256*Mgl*Mt*pow2(Mst2) + 16*s2t*(shiftst1 - shiftst2)*pow2(Mst1)*pow2(
        Mst2)) + 32*Mgl*Mt*(-30*Dmsqst1*(Dmst12*(4*(-3 + lmMgl - lmMst1)*pow2(
        Mgl) + 7*pow2(Mst1)) + 2*((6 - 2*lmMgl + 2*lmMst1)*pow2(Mgl) + 3*pow2(
        Mst1))*pow2(Mst2)) - pow2(Mst1)*(Dmst12*(-2*(-249 + 164*lmMgl - 176*
        lmMst1 - 20*lmMgl*lmMst1 + 27*pow2(lmMgl) + pow2(lmMst1))*pow2(Mgl) + (
        149 + 41*lmMst1)*pow2(Mst1)) + (-4*(-11 + 154*lmMgl - 148*lmMst1 - 20*
        lmMgl*lmMst1 + 27*pow2(lmMgl) + pow2(lmMst1))*pow2(Mgl) + (215 + 118*
        lmMst1 + 41*pow2(lmMst1))*pow2(Mst1))*pow2(Mst2))) - 4*s2t*pow2(Mst1)*(
        2*Dmst12*(200 + 18*lmMgl + 265*lmMst1 + 91*pow2(lmMst1))*pow2(Mgl)*
        pow2(Mst1) + 15*Dmsqst1*(48*Dmst12*pow2(Mgl) + pow2(Mst1)*(Dmst12*(3 -
        6*shiftst1 - 6*shiftst2 + 4*lmMst1*(-2 + shiftst1 + shiftst2)) + 4*(3 -
        2*lmMst1)*(shiftst1 - shiftst2)*pow2(Mst2))) + (-(Dmst12*(32 + 30*(
        shiftst1 + shiftst2) + 2*lmMst1*(115 - 30*(shiftst1 + shiftst2) - 9*
        shiftst3) + 21*shiftst3 + 82*pow2(lmMst1))) + 6*(1 - 2*lmMst1)*(10*
        shiftst1 - 10*shiftst2 + shiftst3)*pow2(Mst2))*pow4(Mst1))) + 64*pow2(
        Mgl)*pow2(Mst1)*(pow2(Dmst12)*(163 + 15*lmMgl - 27*lmMst1 + 20*lmMt -
        8*pow2(lmMst1)) + 20*Dmst12*(16 - 9*lmMgl + 41*lmMst1 - 12*lmMt + 4*
        pow2(lmMst1))*pow2(Mst2) - 10*(5 + 18*lmMgl - 74*lmMst1 + 24*lmMt - 8*
        pow2(lmMst1))*pow4(Mst2))*pow4(Mt)) + 3*z3*(40*Mgl*pow2(Mst1)*pow2(
        Mst2)*(56*Dmst12*s2t*pow2(Mst1)*pow2(Mst2) - Mgl*(648*Dmst12*Mgl*s2t*(
        Dmst12 + 2*pow2(Mst2)) + 7*Mt*(pow2(Dmst12) - 8*pow2(Mst2)*(Dmst12 +
        pow2(Mst2)))))*pow3(Mt) + 3*pow2(Dmst12)*pow2(Mst1)*pow2(Mt)*(-5*(105*
        xDmsqst1*pow2(Dmsqst1) + 140*Dmsqst1*pow2(Mst1) + 48*pow2(Mst1)*(-5*
        pow2(Mgl) + pow2(Mst1)))*pow2(Mst2)*pow2(s2t) - 2*Dmst12*xDmst12*(8*
        Mgl*Mt*(35*Mgl*Mt - 2160*s2t*pow2(Mgl) + 2*s2t*pow2(Mst1)) + 5*pow2(
        s2t)*(-105*xDmsqst1*pow2(Dmsqst1) + (5*(7*Dmsqst1 + 24*pow2(Mgl))*pow2(
        Mst1))/2. + 11*pow4(Mst1)))) + 40*xMgl*pow4(Mgl)*(-6*xDmst12*pow3(
        Dmst12)*(259*pow2(Mst1)*pow2(Mt)*pow2(s2t) + 484*pow4(Mt)) + pow2(Mst2)
        *(777*pow2(Dmst12)*pow2(Mst1)*pow2(Mt)*pow2(s2t) + 968*(3*pow2(Dmst12)
        - 3*Dmst12*pow2(Mst2) - 2*pow4(Mst2))*pow4(Mt))))) + pow2(Mt)*(-(
        xDmst12*pow3(Dmst12)*(Al4p*(-12*pow2(Mst1)*(Al4p*threeLoopFlag*(640*
        Mgl*Mt*s2t*(75*Dmsqst1 - (-130 + lmMgl*(159 - 20*lmMst1) - 162*lmMst1 +
        27*pow2(lmMgl) + pow2(lmMst1))*pow2(Mgl)) + 75*(21 - 16*shiftst2)*
        xDmsqst1*pow2(Dmsqst1)*pow2(s2t) - 32*pow2(Mgl)*((323 - 75*lmMgl + 383*
        lmMst1 - 100*lmMt + 32*pow2(lmMst1))*pow2(Mt) - 450*Dmsqst1*pow2(s2t)))
        + 960*(1 - lmMgl + lmMst1)*Mt*s2t*twoLoopFlag*pow3(Mgl)) + s2t*(2880*
        Mgl*Mt*(40*Al4p*Dmsqst1*threeLoopFlag*(Dmsqst1*xDmsqst1 + (3 - lmMgl +
        lmMst1)*pow2(Mgl)) + (1 - 3*lmMgl + 3*lmMst1)*twoLoopFlag*x*pow4(Mgl))
        + 2*(-36*Mgl*(2*(9 + 2*lmMst1)*Mt + 5*(1 + 2*lmMst1)*Mgl*s2t)*
        twoLoopFlag - Al4p*threeLoopFlag*(45*Dmsqst1*s2t*(5 - 72*shiftst1 - 48*
        shiftst2 + 16*lmMst1*(-5 + 3*shiftst1 + 2*shiftst2)) + 16*Mgl*Mt*(-365
        + 1049*lmMst1 + 123*pow2(lmMst1)) + 60*s2t*(-425 + 18*lmMgl + 83*lmMst1
        + 91*pow2(lmMst1))*pow2(Mgl)))*pow4(Mst1))) + 3*(45*oneLoopFlag + 4*
        Al4p*(120*lmMst1*twoLoopFlag + Al4p*threeLoopFlag*(-605 + 180*shiftst1
        + 120*shiftst2 + 348*shiftst3 - 8*lmMst1*(-110 + 45*shiftst1 + 30*
        shiftst2 + 27*shiftst3) + 820*pow2(lmMst1))))*pow2(s2t)*pow6(Mst1)))/3.
         + 15*Dmst12*s2t*pow2(Mst2)*(3*Dmst12*oneLoopFlag*s2t*pow6(Mst1) - 16*
        Al4p*twoLoopFlag*(4*Mgl*Mt*pow2(Mst2)*((1 - 2*lmMgl + 2*lmMst1)*pow2(
        Mgl)*pow2(Mst1) + (1 - 3*lmMgl + 3*lmMst1)*x*pow4(Mgl) + (2 + lmMst1)*
        pow4(Mst1)) + Dmst12*(2*(3 - 2*lmMgl + 2*lmMst1)*Mt*pow2(Mst1)*pow3(
        Mgl) + 2*Mgl*Mt*pow4(Mst1) + (3 + 2*lmMst1)*s2t*pow2(Mgl)*pow4(Mst1) +
        4*(-1 + 3*lmMgl - 3*lmMst1)*Mt*x*pow5(Mgl) - (1 + 2*lmMst1)*s2t*pow6(
        Mst1)))))))/(2160.*pow6(Mst1)*pow6(Mst2));
}

/**
 * @return The diagonal (2, 2) matrix element of the Higgs mass matrix as a double for the hierarchy 'H9'
 */
double H9::getS2() const {
   return -(-5*xDmst12*pow2(Mst1)*pow3(Dmst12)*(Mt*(54*oneLoopFlag*(-(Mt*Tbeta*
        pow2(s2t)*(pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 18*pow2(Mst1)*pow2(Sbeta))
        ) - 8*Tbeta*pow2(Sbeta)*pow3(Mt) + MuSUSY*pow2(Sbeta)*(8*s2t*pow2(Mt) -
        pow2(Mst1)*pow3(s2t)))*pow6(Mst1) - (144*Al4p*twoLoopFlag*(MuSUSY*pow2(
        Mst1)*(Mt*MuSUSY*s2t*Tbeta*(2*(9 + 2*lmMst1)*Mgl*Mt*pow2(Mst1) + 5*(1 +
        2*lmMst1)*s2t*pow2(Mgl)*pow2(Mst1) + 160*(1 - lmMgl + lmMst1)*Mt*pow3(
        Mgl) - 20*lmMst1*s2t*pow4(Mst1)) + 10*pow2(Sbeta)*(Mt*s2t*pow2(Mst1)*(-
        4*Mt*pow2(Mgl) + 2*(9 - 8*lmMst1)*Mt*pow2(Mst1) + 3*(3 + 2*lmMst1)*Mgl*
        s2t*pow2(Mst1) - 3*(1 + 2*lmMgl - 2*lmMst1)*s2t*pow3(Mgl)) + (4*(1 +
        lmMt)*Mgl*pow2(Mst1) + 8*(18 + 26*lmMst1 - 13*lmMt - 6*lmMst1*lmMt +
        lmMgl*(-13 - 6*lmMst1 + 6*lmMt) + 6*pow2(lmMst1))*pow3(Mgl))*pow3(Mt) +
        (-((3 + 2*lmMst1)*pow2(Mgl)) + (1 + 2*lmMst1)*pow2(Mst1))*pow3(s2t)*
        pow4(Mst1))) - 20*x*(Tbeta*((-12*Mgl*(-((lmMgl*(17 + 6*lmMst1 - 6*lmMt)
        + 5*lmMt + lmMst1*(-22 + 6*lmMt))*Mt) + (8 - 25*lmMst1 + lmMgl*(23 + 6*
        lmMst1 - 6*lmMt) + 2*lmMt + 6*lmMst1*lmMt)*Mgl*s2t + 6*(Mt - Mgl*s2t)*
        pow2(lmMst1))*pow2(Mt) + (1 - 3*lmMgl + 3*lmMst1)*pow2(Mst1)*(12*Mgl*Mt
        + s2t*pow2(Mst1))*pow2(s2t))*pow2(Sbeta) + s2t*pow2(Mt)*(2*(-1 + 3*
        lmMgl - 3*lmMst1)*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 48*(-2 - 9*lmMst1 +
        lmMgl*(7 + 2*lmMst1 - 2*lmMt) + 2*lmMt + 2*lmMst1*lmMt - 2*pow2(lmMst1)
        )*pow2(Mst1)*pow2(Sbeta))) + MuSUSY*pow2(Sbeta)*(6*(-1 + 3*lmMgl - 3*
        lmMst1)*Mt*s2t*(2*Mgl*Mt + s2t*pow2(Mst1)) + 12*(lmMst1*(22 - 6*lmMt) -
        5*lmMt + lmMgl*(-17 - 6*lmMst1 + 6*lmMt) + 6*pow2(lmMst1))*pow3(Mt)))*
        pow5(Mgl) + Tbeta*pow2(Sbeta)*(5*Mt*(4*pow2(Mst1)*(3*(1 + 6*lmMst1)*
        pow2(Mst1) + lmMst1*pow2(MuSUSY)) - pow2(Mgl)*(12*(3 + 2*lmMst1)*pow2(
        Mst1) + (1 + 2*lmMst1)*pow2(MuSUSY)))*pow2(s2t)*pow4(Mst1) + 2*Mgl*s2t*
        pow2(Mst1)*pow2(Mt)*(-((9 + 2*lmMst1)*pow2(Mst1)*pow2(MuSUSY)) + 20*
        pow2(Mgl)*((1 + 4*lmMst1 - 2*(lmMgl + lmMt))*pow2(Mst1) + 4*(-1 + lmMgl
        - lmMst1)*pow2(MuSUSY)) + 40*(1 + 2*lmMt)*pow4(Mst1)) + (40*(-1 + 2*
        lmMst1 + 2*lmMt)*pow3(Mt) - 20*((2 + lmMst1)*Mgl*pow2(Mst1) + (1 - 2*
        lmMgl + 2*lmMst1)*pow3(Mgl))*pow3(s2t))*pow6(Mst1))))/5.) +
        threeLoopFlag*pow2(Al4p)*(36*xDmsqst1*pow2(Dmsqst1)*(480*pow2(Mgl)*
        pow2(Mt)*(4*Mt*s2t*(-6*MuSUSY + Mgl*Tbeta*(16 + 13*lmMst1 + lmMgl*(3 +
        2*lmMst1 - 2*lmMt) - 16*lmMt + 2*lmMst1*lmMt - 2*pow2(lmMst1))) + 3*(1
        - 19*lmMst1 + 19*lmMt)*Tbeta*pow2(Mt) + 24*Tbeta*pow2(Mst1)*pow2(s2t))*
        pow2(Sbeta) - 320*Mgl*Mt*(MuSUSY*pow2(Sbeta)*(12*Mt*pow2(Mst1)*pow2(
        s2t) + (3 - 57*lmMst1 + 57*lmMt)*pow3(Mt)) + Tbeta*(s2t*pow2(Mt)*(4*
        pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 3*(15 + 26*lmMst1 - 26*lmMt)*pow2(
        Mst1)*pow2(Sbeta)) - 2*pow2(Sbeta)*pow3(s2t)*pow4(Mst1))) + (5*pow2(
        Mst1)*(-3*Tbeta*pow2(Mt)*pow2(s2t)*(2*(-21 + 16*shiftst2)*pow2(MuSUSY)*
        (-1 + pow2(Sbeta)) + (-147 + 128*shiftst1 - 32*shiftst2)*pow2(Mst1)*
        pow2(Sbeta)) + pow2(Sbeta)*(144*MuSUSY*s2t*(23 - 4*shiftst2)*pow3(Mt) +
        3*Mt*MuSUSY*(21 - 32*shiftst1 + 16*shiftst2)*pow2(Mst1)*pow3(s2t) + 2*(
        -1469 + 1920*lmMst1 - 1920*lmMt + 288*shiftst2)*Tbeta*pow4(Mt) + 24*(
        shiftst1 - shiftst2)*Tbeta*pow4(Mst1)*pow4(s2t))))/3.) + (-3*Mt*MuSUSY*
        pow2(Sbeta)*(640*Mt*pow2(Mst1)*(-3*(227 + 4*lmMgl*(36 - 5*lmMst1) -
        120*lmMst1 + 27*pow2(lmMgl) + pow2(lmMst1))*pow2(Mst1)*pow2(s2t) + 4*
        pow2(Mt)*(5845 - 2937*lmMt - 27*(16 + lmMgl + 6*lmMst1 - 9*lmMt)*pow2(
        lmMgl) + (346 + 183*lmMt)*pow2(lmMst1) + lmMst1*(4753 - 84*lmMt - 72*
        pow2(lmMt)) - 156*pow2(lmMt) + lmMgl*(-783 - 426*lmMst1*(-1 + lmMt) -
        52*lmMt + 273*pow2(lmMst1) + 72*pow2(lmMt)) - 84*pow3(lmMst1)))*pow3(
        Mgl) - 8*(-40*Mgl*(3*Mt*(66 + 77*lmMst1 + 41*pow2(lmMst1))*pow2(Mst1)*
        pow2(s2t) + (-466 - 6*lmMst1 + 94*lmMt + 4*lmMst1*lmMt + 30*pow2(
        lmMst1) + 48*pow2(lmMt))*pow3(Mt)) - 4*pow2(Mgl)*(2*s2t*(4113 + 90*
        lmMgl - 1792*lmMst1 + 120*lmMt - 48*pow2(lmMst1))*pow2(Mt) - 5*(200 +
        18*lmMgl + 265*lmMst1 + 91*pow2(lmMst1))*pow2(Mst1)*pow3(s2t)))*pow4(
        Mst1) + 40*s2t*((2615 - 672*shiftst3 + 8*lmMst1*(325 + 36*shiftst3) -
        1312*pow2(lmMst1))*pow2(Mt) + 2*(32 + 120*shiftst1 - 60*shiftst2 +
        lmMst1*(230 - 240*shiftst1 + 120*shiftst2 - 42*shiftst3) + 39*shiftst3
        + 82*pow2(lmMst1))*pow2(Mst1)*pow2(s2t))*pow6(Mst1) + 75*Dmsqst1*(-
        3072*pow3(Mgl)*((-3 + lmMgl - lmMst1)*Mt*pow2(Mst1)*pow2(s2t) + (31 +
        27*lmMst1 - 24*lmMt - 2*lmMst1*lmMt + lmMgl*(-3 - 2*lmMst1 + 2*lmMt) +
        2*pow2(lmMst1))*pow3(Mt)) + 249*s2t*pow2(Mt)*pow4(Mst1) + 128*Mgl*(16*(
        11 + 20*lmMst1 - 20*lmMt)*pow2(Mst1)*pow3(Mt) - 3*Mt*pow2(s2t)*pow4(
        Mst1)) - 768*pow2(Mgl)*(60*s2t*pow2(Mst1)*pow2(Mt) + pow3(s2t)*pow4(
        Mst1)) - 16*(3 - 24*shiftst1 + 12*shiftst2 - 8*lmMst1*(1 - 2*shiftst1 +
        shiftst2))*pow3(s2t)*pow6(Mst1))) - Tbeta*(-64*Mgl*pow2(MuSUSY)*pow3(
        Mt)*(3600*Dmsqst1*(3 - lmMgl + lmMst1)*s2t*pow2(Mgl) - 18000*Dmsqst1*
        s2t*pow2(Mst1) + 12*Mgl*Mt*(323 - 75*lmMgl + 383*lmMst1 - 100*lmMt +
        32*pow2(lmMst1))*pow2(Mst1) + 240*s2t*(-130 + lmMgl*(159 - 20*lmMst1) -
        162*lmMst1 + 27*pow2(lmMgl) + pow2(lmMst1))*pow2(Mgl)*pow2(Mst1) + s2t*
        (365 - 1049*lmMst1 - 123*pow2(lmMst1))*pow4(Mst1)) + 60*pow2(Mst1)*
        pow2(s2t)*(-24*pow2(Mst1)*(Dmsqst1*(3 - 2*lmMst1) + (1 - 2*lmMst1)*
        pow2(Mst1))*(pow2(Mt)*(20*shiftst2*pow2(Mst1)*pow2(Sbeta) + pow2(
        MuSUSY)*(3*shiftst1 + 2*shiftst2 - 2*shiftst2*pow2(Sbeta))) - 5*
        shiftst2*pow2(s2t)*pow2(Sbeta)*pow4(Mst1)) + pow2(Mt)*(24*shiftst1*
        pow2(Mst1)*(Dmsqst1*(3 - 2*lmMst1) + (1 - 2*lmMst1)*pow2(Mst1))*(80*
        pow2(Mst1) + 3*pow2(MuSUSY))*pow2(Sbeta) - pow2(MuSUSY)*(-4*(-425 + 18*
        lmMgl + 83*lmMst1 + 91*pow2(lmMst1))*pow2(Mgl)*pow2(Mst1) - 15*Dmsqst1*
        (384*pow2(Mgl) + (1 - 16*lmMst1)*pow2(Mst1)) + 2*(-121 + 176*lmMst1 +
        164*pow2(lmMst1))*pow4(Mst1)))) - 144*shiftst3*(2*pow2(Mt)*pow2(s2t)*((
        -29 + 18*lmMst1)*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 5*(-83 + 58*lmMst1)*
        pow2(Mst1)*pow2(Sbeta)) - 80*(-7 + 3*lmMst1)*pow2(Sbeta)*pow4(Mt) + 5*(
        1 - 2*lmMst1)*pow2(Sbeta)*pow4(Mst1)*pow4(s2t))*pow6(Mst1) + pow2(
        Sbeta)*(-480*Mgl*Mt*pow3(s2t)*pow4(Mst1)*(360*Dmsqst1*pow2(Mst1) + 8*
        pow2(Mgl)*(30*Dmsqst1*(3 - lmMgl + lmMst1) - (-11 + lmMgl*(154 - 20*
        lmMst1) - 148*lmMst1 + 27*pow2(lmMgl) + pow2(lmMst1))*pow2(Mst1)) + 2*(
        215 + 118*lmMst1 + 41*pow2(lmMst1))*pow4(Mst1)) - 7200*shiftst1*(
        Dmsqst1*(3 - 2*lmMst1)*pow2(Mst1) + (1 - 2*lmMst1)*pow4(Mst1))*pow4(
        s2t)*pow6(Mst1) + pow2(Mst1)*(8*(25*Dmsqst1*(-2592*(25 + 47*lmMst1 -
        47*lmMt)*pow2(Mgl) + (907 - 96*lmMst1 + 96*lmMt)*pow2(Mst1)) + 48*pow2(
        Mgl)*(5*(1140 + lmMst1*(747 - 32*lmMt) - 747*lmMt + 16*pow2(lmMst1))*
        pow2(Mst1) + 2*(323 - 75*lmMgl + 383*lmMst1 - 100*lmMt + 32*pow2(
        lmMst1))*pow2(MuSUSY)) + 20*(2243 + 2052*lmMt - 24*lmMst1*(129 + 10*
        lmMt) + 612*(pow2(lmMst1) + pow2(lmMt)))*pow4(Mst1))*pow4(Mt) - 120*
        pow2(Mt)*pow2(s2t)*(-24*pow2(Mgl)*pow2(Mst1)*(90*Dmsqst1 + (-3 + 18*
        lmMgl - 397*lmMst1 + 48*lmMt - 107*pow2(lmMst1))*pow2(Mst1)) + 45*
        Dmsqst1*(9 - 32*lmMst1)*pow4(Mst1) + (pow2(MuSUSY)*(4*(-425 + 18*lmMgl
        + 83*lmMst1 + 91*pow2(lmMst1))*pow2(Mgl)*pow2(Mst1) + 15*Dmsqst1*(384*
        pow2(Mgl) + (1 - 16*lmMst1)*pow2(Mst1)) - 2*(-121 + 176*lmMst1 + 164*
        pow2(lmMst1))*pow4(Mst1)))/2. - 12*(-119 + 250*lmMst1 + 246*pow2(
        lmMst1))*pow6(Mst1))) - 64*Mgl*s2t*pow3(Mt)*((-365 + 1049*lmMst1 + 123*
        pow2(lmMst1))*pow2(MuSUSY)*pow4(Mst1) + 900*Dmsqst1*(20*pow2(Mst1)*
        pow2(MuSUSY) + 4*pow2(Mgl)*(3*(-67 - 43*lmMst1 + lmMgl*(5 + 2*lmMst1 -
        2*lmMt) + 38*lmMt + 2*lmMst1*lmMt - 2*pow2(lmMst1))*pow2(Mst1) + (-3 +
        lmMgl - lmMst1)*pow2(MuSUSY)) + (15 + 2*lmMst1 - 2*lmMt)*pow4(Mst1)) +
        30*(pow2(Mgl)*pow2(Mst1)*((2722 + 965*lmMst1 - 6*(200 + 23*lmMst1)*lmMt
        + lmMgl*(67 - 258*lmMst1 + 178*lmMt) + 108*pow2(lmMgl) + 62*pow2(
        lmMst1) + 48*pow2(lmMt))*pow2(Mst1) - 8*(-130 + lmMgl*(159 - 20*lmMst1)
        - 162*lmMst1 + 27*pow2(lmMgl) + pow2(lmMst1))*pow2(MuSUSY)) - (-223 +
        342*lmMt + lmMst1*(-442 + 8*lmMt) + 60*pow2(lmMst1) + 96*pow2(lmMt))*
        pow6(Mst1))))))/5.)) - (9*Mt*threeLoopFlag*z3*pow2(Al4p)*(75*Dmsqst1*
        pow2(Mst1)*(16*Mt*pow2(Dmst12)*pow2(Mst2)*(1024*Mgl*Mt*pow2(Mst1)*(Mt*
        MuSUSY + s2t*Tbeta*pow2(Mst1))*pow2(Sbeta) - 2048*Mt*(3*Mt*MuSUSY - 2*
        s2t*Tbeta*pow2(Mst1))*pow2(Sbeta)*pow3(Mgl) + 7680*Tbeta*xMgl*pow2(Mt)*
        pow2(Sbeta)*pow4(Mgl) - 7*s2t*(7*Mt*MuSUSY*pow2(Sbeta) + s2t*Tbeta*(
        pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 4*pow2(Mst1)*pow2(Sbeta)))*pow4(Mst1)
        ) + xDmst12*pow3(Dmst12)*(32768*(3*Mt*MuSUSY - 5*s2t*Tbeta*pow2(Mst1))*
        pow2(Mt)*pow2(Sbeta)*pow3(Mgl) - 65536*Mgl*MuSUSY*pow2(Mst1)*pow2(
        Sbeta)*pow3(Mt) + 172032*Tbeta*pow2(Mgl)*pow2(Mst1)*pow2(Sbeta)*pow3(
        Mt) - 122880*Tbeta*xMgl*pow2(Sbeta)*pow3(Mt)*pow4(Mgl) + 7*(31*MuSUSY*
        s2t*pow2(Mt)*pow2(Sbeta) + 4*Mt*Tbeta*pow2(s2t)*(pow2(MuSUSY) + 2*pow2(
        Mst1)*pow2(Sbeta) - pow2(MuSUSY)*pow2(Sbeta)) - 24*Tbeta*pow2(Sbeta)*
        pow3(Mt) + 16*MuSUSY*pow2(Mst1)*pow2(Sbeta)*pow3(s2t))*pow4(Mst1)) -
        128*Dmst12*pow2(Mt)*pow2(Sbeta)*(-256*(3*Mt*MuSUSY + s2t*Tbeta*pow2(
        Mst1))*pow3(Mgl) + 960*Mt*Tbeta*xMgl*pow4(Mgl) + 7*MuSUSY*s2t*pow4(
        Mst1) - 128*Mgl*(2*Mt*MuSUSY*pow2(Mst1) + s2t*Tbeta*pow4(Mst1)))*pow4(
        Mst2) - 8192*Mgl*pow2(Sbeta)*(-4*MuSUSY*pow2(Mgl) - 2*MuSUSY*pow2(Mst1)
        + 5*Tbeta*xMgl*pow3(Mgl))*pow3(Mt)*pow6(Mst2)) + 300*xDmsqst1*pow2(
        Dmsqst1)*(-(Mt*pow2(Dmst12)*pow2(Mst1)*pow2(Mst2)*(2048*Mgl*Mt*(3*Mt*
        MuSUSY - s2t*Tbeta*pow2(Mst1))*pow2(Sbeta) - 9216*Tbeta*pow2(Mgl)*pow2(
        Mt)*pow2(Sbeta) + pow2(Mst1)*(14*Mt*MuSUSY*s2t*pow2(Sbeta) + 491*Tbeta*
        pow2(Mt)*pow2(Sbeta) + 7*Tbeta*pow2(s2t)*(3*pow2(MuSUSY)*(-1 + pow2(
        Sbeta)) + 14*pow2(Mst1)*pow2(Sbeta))) + 12288*Mt*s2t*Tbeta*pow2(Sbeta)*
        pow3(Mgl))) + xDmst12*pow2(Mst1)*pow3(Dmst12)*(2048*Mgl*(3*Mt*MuSUSY -
        4*s2t*Tbeta*pow2(Mst1))*pow2(Mt)*pow2(Sbeta) + 12288*s2t*Tbeta*pow2(Mt)
        *pow2(Sbeta)*pow3(Mgl) - 9216*Tbeta*pow2(Mgl)*pow2(Sbeta)*pow3(Mt) +
        pow2(Mst1)*(224*MuSUSY*s2t*pow2(Mt)*pow2(Sbeta) + 7*Mt*Tbeta*pow2(s2t)*
        (6*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 13*pow2(Mst1)*pow2(Sbeta)) + 1838*
        Tbeta*pow2(Sbeta)*pow3(Mt) + 21*MuSUSY*pow2(Mst1)*pow2(Sbeta)*pow3(s2t)
        )) + 4*Dmst12*pow2(Mst1)*pow2(Mt)*pow2(Sbeta)*(-2304*Mt*Tbeta*pow2(Mgl)
        - (49*MuSUSY*s2t + 214*Mt*Tbeta)*pow2(Mst1) + 512*Mgl*(3*Mt*MuSUSY + 2*
        s2t*Tbeta*pow2(Mst1)) + 3072*s2t*Tbeta*pow3(Mgl))*pow4(Mst2) - 16*pow2(
        Sbeta)*pow3(Mt)*(-256*Mgl*MuSUSY*pow2(Mst1) + 384*Tbeta*pow2(Mgl)*pow2(
        Mst1) - 768*MuSUSY*pow3(Mgl) + 960*Tbeta*xMgl*pow4(Mgl) + 43*Tbeta*
        pow4(Mst1))*pow6(Mst2)) - 8*pow2(Mst1)*(-80*Dmst12*pow2(Mt)*pow4(Mst2)*
        (2*MuSUSY*pow2(Mgl)*pow2(Mst1)*(-7*Mt*MuSUSY*Tbeta + 12*s2t*pow2(Mst1)*
        pow2(Sbeta)) - 4*pow2(Mst1)*(-100*Mt*MuSUSY*pow2(Sbeta) + s2t*Tbeta*(
        81*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 380*pow2(Mst1)*pow2(Sbeta)))*pow3(
        Mgl) + 2*xMgl*(-363*Mt*Tbeta*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 2*(279*
        MuSUSY*s2t + 58*Mt*Tbeta)*pow2(Mst1)*pow2(Sbeta))*pow4(Mgl) + 2*Mgl*(
        496*Mt*MuSUSY*pow2(Sbeta) + s2t*Tbeta*(7*pow2(MuSUSY)*(-1 + pow2(Sbeta)
        ) + 32*pow2(Mst1)*pow2(Sbeta)))*pow4(Mst1) + 93*MuSUSY*s2t*pow2(Sbeta)*
        pow6(Mst1)) - 20*Mt*pow2(Dmst12)*pow2(Mst2)*(Tbeta*pow2(Mgl)*pow2(Mst1)
        *(7*pow2(Mt)*pow2(MuSUSY) + 6*pow2(Mst1)*pow2(s2t)*(15*pow2(MuSUSY)*(-1
        + pow2(Sbeta)) + 8*pow2(Mst1)*pow2(Sbeta))) + 8*pow2(Mst1)*(2192*
        MuSUSY*pow2(Mt)*pow2(Sbeta) + 243*MuSUSY*pow2(Mst1)*pow2(s2t)*pow2(
        Sbeta) + 3*Mt*s2t*Tbeta*(-27*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 320*
        pow2(Mst1)*pow2(Sbeta)))*pow3(Mgl) + 4*Mgl*(928*Mt*s2t*Tbeta*pow2(Mst1)
        + 1136*MuSUSY*pow2(Mt) - 21*MuSUSY*pow2(Mst1)*pow2(s2t))*pow2(Sbeta)*
        pow4(Mst1) + xMgl*pow4(Mgl)*(2904*Tbeta*pow2(Mt)*pow2(MuSUSY)*(-1 +
        pow2(Sbeta)) + pow2(Mst1)*(777*Tbeta*pow2(MuSUSY)*pow2(s2t)*(-1 + pow2(
        Sbeta)) + 2232*Mt*MuSUSY*s2t*pow2(Sbeta) - 35536*Tbeta*pow2(Mt)*pow2(
        Sbeta)) + 2232*Tbeta*pow2(s2t)*pow2(Sbeta)*pow4(Mst1)) + 6*s2t*(-66*Mt*
        MuSUSY*pow2(Sbeta) + s2t*Tbeta*(-3*pow2(MuSUSY)*(-1 + pow2(Sbeta)) +
        31*pow2(Mst1)*pow2(Sbeta)))*pow6(Mst1)) + xDmst12*pow3(Dmst12)*(20*
        pow2(Mgl)*pow2(Mst1)*(Mt*pow2(Mst1)*(45*Tbeta*pow2(MuSUSY)*pow2(s2t)*(-
        1 + pow2(Sbeta)) - 14*Mt*MuSUSY*s2t*pow2(Sbeta) + 4480*Tbeta*pow2(Mt)*
        pow2(Sbeta)) + 42*Tbeta*pow2(MuSUSY)*(-1 + pow2(Sbeta))*pow3(Mt) + 6*(
        15*MuSUSY*s2t + 8*Mt*Tbeta)*pow2(s2t)*pow2(Sbeta)*pow4(Mst1)) + 20*
        xMgl*pow4(Mgl)*(-2*Mt*pow2(Mst1)*(-777*Tbeta*pow2(MuSUSY)*pow2(s2t)*(-1
        + pow2(Sbeta)) - 4464*Mt*MuSUSY*s2t*pow2(Sbeta) + 35072*Tbeta*pow2(Mt)*
        pow2(Sbeta)) + 2904*Tbeta*pow2(MuSUSY)*(-1 + pow2(Sbeta))*pow3(Mt) + 3*
        (259*MuSUSY*s2t + 372*Mt*Tbeta)*pow2(s2t)*pow2(Sbeta)*pow4(Mst1)) + 16*
        Mgl*pow4(Mst1)*(-105*Mt*MuSUSY*pow2(Mst1)*pow2(s2t)*pow2(Sbeta) + s2t*
        Tbeta*pow2(Mt)*(3*pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 1040*pow2(Mst1)*
        pow2(Sbeta)) + 440*MuSUSY*pow2(Sbeta)*pow3(Mt) + 35*Tbeta*pow2(Sbeta)*
        pow3(s2t)*pow4(Mst1)) - 80*pow2(Mst1)*pow3(Mgl)*(-243*Mt*MuSUSY*pow2(
        Mst1)*pow2(s2t)*pow2(Sbeta) + 8*s2t*Tbeta*pow2(Mt)*(81*pow2(MuSUSY)*(-1
        + pow2(Sbeta)) + 308*pow2(Mst1)*pow2(Sbeta)) - 9168*MuSUSY*pow2(Sbeta)*
        pow3(Mt) + 162*Tbeta*pow2(Sbeta)*pow3(s2t)*pow4(Mst1)) + 15*(-85*
        MuSUSY*s2t*pow2(Mt)*pow2(Sbeta) + Mt*Tbeta*pow2(s2t)*(11*pow2(MuSUSY)*(
        -1 + pow2(Sbeta)) + 512*pow2(Mst1)*pow2(Sbeta)) - 500*Tbeta*pow2(Sbeta)
        *pow3(Mt) - 24*MuSUSY*pow2(Mst1)*pow2(Sbeta)*pow3(s2t))*pow6(Mst1)) -
        160*Mgl*pow3(Mt)*(-7*Mgl*Tbeta*pow2(Mst1)*pow2(MuSUSY) - 760*MuSUSY*
        pow2(Mgl)*pow2(Mst1)*pow2(Sbeta) + 2*Tbeta*xMgl*(-121*pow2(MuSUSY)*(-1
        + pow2(Sbeta)) + 658*pow2(Mst1)*pow2(Sbeta))*pow3(Mgl) + 32*MuSUSY*
        pow2(Sbeta)*pow4(Mst1))*pow6(Mst2))))/2. + Al4p*(-720*Tbeta*pow2(Sbeta)
        *pow4(Mst1)*(384*twoLoopFlag*xDmst12*xMgl*z2*pow3(Dmst12)*pow4(Mgl) -
        Al4p*threeLoopFlag*z3*pow2(Mst2)*(8*Dmst12*pow2(Mst2)*(pow2(Mgl)*(796*
        pow2(Mst1) - 7*pow2(MuSUSY)) + 228*pow4(Mst1)) + pow2(Dmst12)*(pow2(
        Mgl)*(14080*pow2(Mst1) + 7*pow2(MuSUSY)) + 1141*pow4(Mst1)) + 8*(pow2(
        Mgl)*(152*pow2(Mst1) - 7*pow2(MuSUSY)) + 2*(395 - 12*lmMst1 + 12*lmMt)*
        pow4(Mst1))*pow4(Mst2) + 15*Dmsqst1*(pow2(Dmst12)*(1920*pow2(Mgl) +
        107*pow2(Mst1)) + 2*Dmst12*(768*pow2(Mgl) + 107*pow2(Mst1))*pow2(Mst2)
        + 256*(3*pow2(Mgl) + pow2(Mst1))*pow4(Mst2))))*pow4(Mt) + 480*xMgl*
        pow4(Mgl)*(-3*Mt*twoLoopFlag*xDmst12*pow3(Dmst12)*(-2*(-1 + 3*lmMgl -
        3*lmMst1)*Mt*Tbeta*pow2(MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta)) + pow2(
        Sbeta)*(96*(1 - lmMgl + lmMst1)*MuSUSY*s2t*pow2(Mt) - 6*(1 + 2*lmMgl -
        2*lmMst1)*Mt*Tbeta*pow2(Mst1)*pow2(s2t) + 48*(5 - 2*lmMgl*(2 + lmMst1 -
        lmMt) + 2*lmMst1*(4 + lmMst1 - lmMt) - 4*lmMt)*Tbeta*pow3(Mt) + (1 - 3*
        lmMgl + 3*lmMst1)*MuSUSY*pow2(Mst1)*pow3(s2t)))*pow4(Mst1) + pow2(Mst2)
        *pow2(Mt)*(3*twoLoopFlag*pow4(Mst1)*(s2t*pow2(Dmst12)*(12*(3 - 2*lmMgl
        + 2*lmMst1)*Mt*MuSUSY*pow2(Sbeta) + s2t*Tbeta*(-((-1 + 3*lmMgl - 3*
        lmMst1)*pow2(MuSUSY)*(-1 + pow2(Sbeta))) + 12*(1 - 2*lmMgl + 2*lmMst1)*
        pow2(Mst1)*pow2(Sbeta))) + 12*Mt*pow2(Sbeta)*(2*Dmst12*(1 - 2*lmMgl +
        2*lmMst1)*MuSUSY*s2t*pow2(Mst2) + Mt*Tbeta*(pow2(Dmst12)*(8 - 2*lmMst1*
        (-5 + lmMt) - 5*lmMt + lmMgl*(-5 - 2*lmMst1 + 2*lmMt) + 2*pow2(lmMst1))
        + 2*Dmst12*(2 - 2*lmMst1*(-3 + lmMt) - 3*lmMt + lmMgl*(-3 - 2*lmMst1 +
        2*lmMt) + 2*pow2(lmMst1))*pow2(Mst2) - 4*(lmMgl*(1 + lmMst1 - lmMt) +
        lmMst1*(-2 + lmMt) + lmMt - pow2(lmMst1))*pow4(Mst2)))) + (Al4p*
        threeLoopFlag*pow2(Mst1)*(864*Dmst12*Mt*MuSUSY*s2t*(30*Dmsqst1*(29 -
        10*lmMgl + 10*lmMst1) + (633 + 1458*lmMst1 - 96*(1 + lmMst1)*lmMt + 2*
        lmMgl*(-632 + 99*lmMst1 + 48*lmMt) - 333*pow2(lmMgl) + 103*pow2(lmMst1)
        )*pow2(Mst1))*pow2(Mst2)*pow2(Sbeta) + s2t*pow2(Dmst12)*(25920*Dmsqst1*
        (-29 + 10*lmMgl - 10*lmMst1)*Mt*MuSUSY*pow2(Sbeta) + s2t*Tbeta*pow2(
        Mst1)*(12960*Dmsqst1*(29 - 10*lmMgl + 10*lmMst1) + 432*(681 + 1506*
        lmMst1 - 96*(1 + lmMst1)*lmMt + 2*lmMgl*(-656 + 99*lmMst1 + 48*lmMt) -
        333*pow2(lmMgl) + 103*pow2(lmMst1))*pow2(Mst1) + (41653 + 11784*lmMst1
        + 48*lmMgl*(-313 + 180*lmMst1) - 12636*pow2(lmMgl) + 540*pow2(lmMst1))*
        pow2(MuSUSY))*pow2(Sbeta) + MuSUSY*pow2(Mst1)*(MuSUSY*s2t*Tbeta*(-41653
        - 11784*lmMst1 - 48*lmMgl*(-313 + 180*lmMst1) + 12636*pow2(lmMgl) -
        540*pow2(lmMst1)) + 432*Mt*(2961 + lmMst1*(1964 - 96*lmMt) - 192*lmMt +
        2*lmMgl*(-683 + 99*lmMst1 + 48*lmMt) - 333*pow2(lmMgl) + 103*pow2(
        lmMst1))*pow2(Sbeta))) + Tbeta*pow2(Mt)*(-576*pow2(MuSUSY)*(Dmst12*(262
        + lmMst1*(335 - 72*lmMt) - 96*lmMt + lmMgl*(-199 - 120*lmMst1 + 72*
        lmMt) + 120*pow2(lmMst1))*(Dmst12 - pow2(Mst2)) - 2*(59 + 85*lmMst1 -
        lmMgl*(53 + 40*lmMst1 - 24*lmMt) - 24*(1 + lmMst1)*lmMt + 40*pow2(
        lmMst1))*pow4(Mst2)) - 72*pow2(Sbeta)*(2*pow2(Dmst12)*(90*Dmsqst1*(307
        + 4*lmMst1*(67 - 5*lmMt) - 10*lmMgl*(3 + 2*lmMst1 - 2*lmMt) - 238*lmMt
        + 20*pow2(lmMst1)) - pow2(Mst1)*(40279 - 17847*lmMt - 6*(303 + 108*
        lmMst1 - 203*lmMt)*pow2(lmMgl) + (952 + 870*lmMt)*pow2(lmMst1) +
        lmMst1*(26599 + 804*lmMt - 288*pow2(lmMt)) - 720*pow2(lmMt) + 2*lmMgl*(
        -1474 - 977*lmMt - 12*lmMst1*(-118 + 87*lmMt) + 615*pow2(lmMst1) + 144*
        pow2(lmMt)) - 190*pow3(lmMgl) - 392*pow3(lmMst1))) - 4*Dmst12*pow2(
        Mst2)*(45*Dmsqst1*(307 + 4*lmMst1*(67 - 5*lmMt) - 10*lmMgl*(3 + 2*
        lmMst1 - 2*lmMt) - 238*lmMt + 20*pow2(lmMst1)) + pow2(Mst1)*(13484 -
        5541*lmMt + 6*(324 + 143*lmMt)*pow2(lmMst1) - 6*lmMgl*(666 + lmMgl*(285
        + 108*lmMst1 - 265*lmMt) + 57*lmMt + 4*lmMst1*(-59 + 102*lmMt) - 203*
        pow2(lmMst1) - 48*pow2(lmMt)) - 432*pow2(lmMt) - 3*lmMst1*(-3979 + 260*
        lmMt + 96*pow2(lmMt)) - 314*pow3(lmMgl) - 256*pow3(lmMst1))) + (-360*
        Dmsqst1*(69 + lmMst1*(59 - 10*lmMt) - 10*lmMgl*(1 + lmMst1 - lmMt) -
        49*lmMt + 10*pow2(lmMst1)) - 4*pow2(Mst1)*(6677 - 2742*lmMt + 6*(313 +
        143*lmMt)*pow2(lmMst1) - 6*lmMgl*(645 + lmMgl*(231 + 108*lmMst1 - 265*
        lmMt) - 72*lmMt + 12*lmMst1*(-9 + 34*lmMt) - 203*pow2(lmMst1) - 48*
        pow2(lmMt)) - 288*pow2(lmMt) - 12*lmMst1*(-617 + 99*lmMt + 24*pow2(
        lmMt)) - 314*pow3(lmMgl) - 256*pow3(lmMst1)))*pow4(Mst2) - 8*pow2(
        MuSUSY)*(Dmst12*(262 + lmMst1*(335 - 72*lmMt) - 96*lmMt + lmMgl*(-199 -
        120*lmMst1 + 72*lmMt) + 120*pow2(lmMst1))*(Dmst12 - pow2(Mst2)) - 2*(59
        + 85*lmMst1 - lmMgl*(53 + 40*lmMst1 - 24*lmMt) - 24*(1 + lmMst1)*lmMt +
        40*pow2(lmMst1))*pow4(Mst2))))))/72.) + Al4p*threeLoopFlag*((xDmst12*
        pow2(Mst1)*pow3(Dmst12)*(Mt*MuSUSY*s2t*pow2(Sbeta)*(1728*(15*Dmsqst1*(
        29 - 10*lmMgl + 10*lmMst1) + (-1797 - 1711*lmMst1 + lmMgl*(1315 - 198*
        lmMst1 - 96*lmMt) + 144*lmMt + 96*lmMst1*lmMt + 333*pow2(lmMgl) - 103*
        pow2(lmMst1))*pow2(Mst1))*pow2(Mt) + (-41653 - 11784*lmMst1 - 48*lmMgl*
        (-313 + 180*lmMst1) + 12636*pow2(lmMgl) - 540*pow2(lmMst1))*pow2(s2t)*
        pow4(Mst1)) + Tbeta*(2*pow2(MuSUSY)*((41653 + 11784*lmMst1 + 48*lmMgl*(
        -313 + 180*lmMst1) - 12636*pow2(lmMgl) + 540*pow2(lmMst1))*pow2(Mst1)*
        pow2(Mt)*pow2(s2t) + 288*(262 + 335*lmMst1 - 96*lmMt - 72*lmMst1*lmMt +
        lmMgl*(-199 - 120*lmMst1 + 72*lmMt) + 120*pow2(lmMst1))*pow4(Mt)) +
        pow2(Sbeta)*(-2*pow2(Mst1)*pow2(Mt)*(12960*Dmsqst1*(29 - 10*lmMgl + 10*
        lmMst1) - 108*(1695 + 8*lmMst1*(-125 + 12*lmMt) - 2*lmMgl*(-605 + 99*
        lmMst1 + 48*lmMt) + 333*pow2(lmMgl) - 103*pow2(lmMst1))*pow2(Mst1) + (
        41653 + 11784*lmMst1 + 48*lmMgl*(-313 + 180*lmMst1) - 12636*pow2(lmMgl)
        + 540*pow2(lmMst1))*pow2(MuSUSY))*pow2(s2t) + 288*(45*Dmsqst1*(307 + 4*
        lmMst1*(67 - 5*lmMt) - 10*lmMgl*(3 + 2*lmMst1 - 2*lmMt) - 238*lmMt +
        20*pow2(lmMst1)) - 2*(262 + lmMst1*(335 - 72*lmMt) - 96*lmMt + lmMgl*(-
        199 - 120*lmMst1 + 72*lmMt) + 120*pow2(lmMst1))*pow2(MuSUSY) + pow2(
        Mst1)*(72*(49 + 18*lmMst1 - 39*lmMt)*pow2(lmMgl) - 16*(181 + 108*lmMt)*
        pow2(lmMst1) + 8*lmMgl*(868 + 287*lmMt + 9*lmMst1*(-59 + 63*lmMt) -
        306*pow2(lmMst1) - 72*pow2(lmMt)) + 8*lmMst1*(-4817 - 3*lmMt + 72*pow2(
        lmMt)) + 3*(-17921 + 7796*lmMt + 384*pow2(lmMt)) + 504*pow3(lmMgl) +
        648*pow3(lmMst1)))*pow4(Mt)))))/72. + 180*Tbeta*xDmsqst1*pow2(Dmsqst1)*
        (41 + 10*lmMgl*(1 + lmMst1 - lmMt) - 81*lmMt + lmMst1*(71 + 10*lmMt) -
        10*pow2(lmMst1))*pow2(Sbeta)*pow4(Mt)*pow6(Mst2))) - 96*z2*(pow2(Mst1)*
        (1440*twoLoopFlag*xDmst12*(3*x*pow2(Mgl)*(Mt*(MuSUSY - Mgl*Tbeta) +
        s2t*Tbeta*pow2(Mgl)) - (Mt*MuSUSY + 4*s2t*Tbeta*x*pow2(Mgl))*pow2(Mst1)
        )*pow2(Sbeta)*pow3(Dmst12)*pow3(Mgl)*pow3(Mt) - 5*Al4p*threeLoopFlag*(
        45*xDmsqst1*xDmst12*pow2(Dmsqst1)*pow3(Dmst12)*(Mt*MuSUSY*pow2(Sbeta)*(
        -8*s2t*(9*pow2(Mgl) - 4*pow2(Mst1))*pow2(Mt) - 24*Mgl*Mt*pow2(Mst1)*
        pow2(s2t) + 16*(17 + 6*lmMst1 - 6*lmMt)*Mgl*pow3(Mt) + pow3(s2t)*pow4(
        Mst1)) + Tbeta*(-2*pow2(Mst1)*pow2(Mt)*pow2(MuSUSY)*pow2(s2t) + 2*pow2(
        Mt)*(-4*(3*(17 + 6*lmMst1 - 6*lmMt)*pow2(Mgl) - 4*(2 + lmMst1 - lmMt)*
        pow2(Mst1))*pow2(Mt) + 4*Mgl*Mt*s2t*(4*(25 + 6*lmMst1 - 6*lmMt)*pow2(
        Mgl) - 8*(7 + 2*lmMst1 - 2*lmMt)*pow2(Mst1) - pow2(MuSUSY)) + pow2(
        Mst1)*(36*pow2(Mgl) + 2*pow2(Mst1) + pow2(MuSUSY))*pow2(s2t))*pow2(
        Sbeta) + 4*Mgl*(2*s2t*pow2(MuSUSY)*pow3(Mt) + Mt*pow2(Sbeta)*pow3(s2t)*
        pow4(Mst1)))) - xMgl*pow2(Mst2)*pow2(Mt)*pow4(Mgl)*(9*Dmst12*Mt*MuSUSY*
        s2t*(400*Dmsqst1 + (16*(-194 + 27*lmMgl - 27*lmMst1)*pow2(Mst1))/3.)*
        pow2(Mst2)*pow2(Sbeta) + 12*s2t*pow2(Dmst12)*(-150*Dmsqst1*(2*Mt*MuSUSY
        - s2t*Tbeta*pow2(Mst1))*pow2(Sbeta) + pow2(Mst1)*(2*(-146 + 27*lmMgl -
        27*lmMst1)*Mt*MuSUSY*pow2(Sbeta) + s2t*Tbeta*((-49 + 17*lmMgl - 17*
        lmMst1)*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 2*(-194 + 27*lmMgl - 27*
        lmMst1)*pow2(Mst1)*pow2(Sbeta)))) + Tbeta*pow2(Mt)*(32*pow2(MuSUSY)*(9*
        Dmst12*(17 - 4*lmMgl + 4*lmMst1)*(Dmst12 - pow2(Mst2)) + 2*(-47 + 12*
        lmMgl - 12*lmMst1)*pow4(Mst2)) - pow2(Sbeta)*(-12*pow2(Dmst12)*((2997 +
        108*lmMgl + 798*lmMst1 - 1214*lmMt)*pow2(Mst1) + 24*(-17 + 4*lmMgl - 4*
        lmMst1)*pow2(MuSUSY)) + 7200*Dmsqst1*(4 + lmMst1 - lmMt)*(3*pow2(
        Dmst12) - 3*Dmst12*pow2(Mst2) - pow4(Mst2)) + 8*(Dmst12*pow2(Mst2)*((
        4265 - 324*lmMgl + 1266*lmMst1 + 366*lmMt)*pow2(Mst1) + 36*(-17 + 4*
        lmMgl - 4*lmMst1)*pow2(MuSUSY)) + ((5207 - 324*lmMgl + 1716*lmMst1 -
        84*lmMt)*pow2(Mst1) + 8*(-47 + 12*lmMgl - 12*lmMst1)*pow2(MuSUSY))*
        pow4(Mst2))))))) + pow2(Mst2)*(1440*twoLoopFlag*pow2(Sbeta)*pow3(Mgl)*
        pow3(Mt)*(pow2(Mst2)*(Dmst12*s2t*Tbeta*pow2(Mst1) + Mt*MuSUSY*(Dmst12 +
        pow2(Mst2)))*pow4(Mst1) + x*pow2(Mgl)*(Dmst12*pow2(Mst1)*(Dmst12*(-3*
        Mt*(MuSUSY - Mgl*Tbeta) + s2t*Tbeta*(-3*pow2(Mgl) + pow2(Mst1))) + (3*
        Mt*(MuSUSY - Mgl*Tbeta) + s2t*Tbeta*(3*pow2(Mgl) + 2*pow2(Mst1)))*pow2(
        Mst2)) + Mt*(MuSUSY - Mgl*Tbeta)*(3*pow2(Mgl) + 2*pow2(Mst1))*pow4(
        Mst2))) + Al4p*threeLoopFlag*pow2(Mst1)*(15*MuSUSY*pow2(Sbeta)*(4*
        Dmst12*s2t*pow2(Mst1)*pow2(Mt)*(-(Dmst12*Mgl*s2t*((113 - 90*lmMgl + 90*
        lmMst1)*pow2(Mgl)*pow2(Mst1) + 30*Dmsqst1*(2*pow2(Mgl) + pow2(Mst1)) -
        17*pow4(Mst1))) + 6*Mt*(5*Dmst12*(3*pow2(Mgl)*pow2(Mst1) + 2*Dmsqst1*(
        6*pow2(Mgl) + pow2(Mst1)) + pow4(Mst1)) + pow2(Mst2)*(30*Dmsqst1*pow2(
        Mgl) + 20*Dmsqst1*pow2(Mst1) - 59*pow2(Mgl)*pow2(Mst1) + 10*pow4(Mst1))
        )) + 30*Mt*s2t*(Dmsqst1 + pow2(Mst1))*pow4(Mst1)*(-8*Dmst12*shiftst2*
        pow2(Mst2)*pow2(Mt) + (-shiftst1 + shiftst2)*pow2(Dmst12)*pow2(Mst1)*
        pow2(s2t) + 8*(shiftst1 - shiftst2)*pow2(Mt)*pow4(Mst2)) - (8*Mgl*(-
        180*Dmsqst1*(pow2(Dmst12)*(12*(4 + lmMst1 - lmMt)*pow2(Mgl) + (-7 - 2*
        lmMst1 + 2*lmMt)*pow2(Mst1)) - 2*Dmst12*(6*(4 + lmMst1 - lmMt)*pow2(
        Mgl) + (7 + 2*lmMst1 - 2*lmMt)*pow2(Mst1))*pow2(Mst2) - 2*(4 + lmMst1 -
        lmMt)*(2*pow2(Mgl) + pow2(Mst1))*pow4(Mst2)) + pow2(Mst1)*(pow2(Dmst12)
        *((3147 + 822*lmMst1 - 822*lmMt)*pow2(Mgl) + (691 + 216*lmMst1 - 216*
        lmMt)*pow2(Mst1)) - 2*Dmst12*((1627 + 234*lmMst1 + 174*lmMt)*pow2(Mgl)
        + (-457 - 87*lmMst1 + 87*lmMt)*pow2(Mst1))*pow2(Mst2) + 4*((-1070 -
        207*lmMst1 + 3*lmMt)*pow2(Mgl) + (100 - 3*lmMst1 + 3*lmMt)*pow2(Mst1))*
        pow4(Mst2)))*pow4(Mt))/3. + 3*Mt*s2t*shiftst3*(-(pow2(Dmst12)*pow2(
        Mst1)*pow2(s2t)) + 8*pow2(Mt)*(pow2(Dmst12) - Dmst12*pow2(Mst2) + pow4(
        Mst2)))*pow6(Mst1)) + Tbeta*(pow2(Mt)*pow2(MuSUSY)*(pow2(Mst1)*(4*pow2(
        Mgl)*(11*pow2(Dmst12) - 520*pow2(Mst2)*(Dmst12 + pow2(Mst2)))*pow2(Mt)
        + 5*pow2(Dmst12)*pow2(s2t)*(274*pow2(Mgl)*pow2(Mst1) - 90*Dmsqst1*(
        pow2(Mgl) + pow2(Mst1)) - 41*pow4(Mst1))) - 5*Dmst12*Mt*s2t*(8*Mgl*
        pow2(Mst2)*((113 - 90*lmMgl + 90*lmMst1)*pow2(Mgl)*pow2(Mst1) + 30*
        Dmsqst1*(2*pow2(Mgl) + pow2(Mst1)) - 17*pow4(Mst1)) + 4*Dmst12*((263 -
        90*lmMgl + 90*lmMst1)*pow2(Mst1)*pow3(Mgl) - 120*Dmsqst1*(-(Mgl*pow2(
        Mst1)) + pow3(Mgl)) + 30*Mgl*pow4(Mst1)))) + 150*(Dmsqst1 + pow2(Mst1))
        *pow2(Mt)*pow4(Mst1)*(pow2(Dmst12)*pow2(s2t)*(-((shiftst1 + shiftst2)*
        pow2(MuSUSY)*(-1 + pow2(Sbeta))) + 6*(-3*shiftst1 + shiftst2)*pow2(
        Mst1)*pow2(Sbeta)) + 2*Dmst12*pow2(Mst2)*((shiftst1 - shiftst2)*pow2(
        MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta)) + 12*(shiftst2*pow2(Mt) + (
        shiftst1 - shiftst2)*pow2(Mst1)*pow2(s2t))*pow2(Sbeta)) + 24*(shiftst1
        + shiftst2)*pow2(Mt)*pow2(Sbeta)*pow4(Mst2)) + 15*shiftst3*(-(Dmst12*
        pow2(Mt)*pow2(s2t)*((3*Dmst12 - 2*pow2(Mst2))*pow2(MuSUSY)*(-1 + pow2(
        Sbeta)) + 6*pow2(Mst1)*(7*Dmst12 - 4*pow2(Mst2))*pow2(Sbeta))) + 24*
        pow2(Sbeta)*(pow2(Dmst12) - Dmst12*pow2(Mst2) + pow4(Mst2))*pow4(Mt))*
        pow6(Mst1) + pow2(Sbeta)*(-5*Dmst12*s2t*pow3(Mt)*(Mgl*pow2(Mst1)*(8*
        Dmst12*(18*(67 + 20*lmMst1 - 20*lmMt)*pow2(Mgl)*pow2(Mst1) + 180*
        Dmsqst1*(4*(9 + 2*lmMst1 - 2*lmMt)*pow2(Mgl) + (7 + 2*lmMst1 - 2*lmMt)*
        pow2(Mst1)) + (607 + 186*lmMst1 - 186*lmMt)*pow4(Mst1)) - 16*pow2(Mst2)
        *((2143 + 414*lmMst1 - 6*lmMt)*pow2(Mgl)*pow2(Mst1) - 90*Dmsqst1*(9 +
        2*lmMst1 - 2*lmMt)*(2*pow2(Mgl) + pow2(Mst1)) + (-197 + 6*lmMst1 - 6*
        lmMt)*pow4(Mst1))) - pow2(MuSUSY)*(8*Mgl*pow2(Mst2)*((113 - 90*lmMgl +
        90*lmMst1)*pow2(Mgl)*pow2(Mst1) + 30*Dmsqst1*(2*pow2(Mgl) + pow2(Mst1))
        - 17*pow4(Mst1)) + 4*Dmst12*((263 - 90*lmMgl + 90*lmMst1)*pow2(Mst1)*
        pow3(Mgl) - 120*Dmsqst1*(-(Mgl*pow2(Mst1)) + pow3(Mgl)) + 30*Mgl*pow4(
        Mst1)))) + pow2(Mst1)*((4*pow2(Dmst12)*(pow2(Mgl)*(40*(613 + 165*lmMst1
        - 165*lmMt)*pow2(Mst1) - 11*pow2(MuSUSY)) + 230*(5 + 3*lmMst1 - 3*lmMt)
        *pow4(Mst1)) + 80*Dmst12*pow2(Mst2)*(pow2(Mgl)*((679 + 150*lmMst1 -
        150*lmMt)*pow2(Mst1) + 26*pow2(MuSUSY)) + 45*(2 + lmMst1 - lmMt)*pow4(
        Mst1)) + 40*(pow2(Mgl)*((496 + 60*lmMst1 - 60*lmMt)*pow2(Mst1) + 52*
        pow2(MuSUSY)) + (223 + 180*lmMst1 - 180*lmMt)*pow4(Mst1))*pow4(Mst2) +
        1800*Dmsqst1*(pow2(Dmst12)*(15*(7 + 2*lmMst1 - 2*lmMt)*pow2(Mgl) + (5 +
        2*lmMst1 - 2*lmMt)*pow2(Mst1)) + 2*Dmst12*(6*(7 + 2*lmMst1 - 2*lmMt)*
        pow2(Mgl) + (3 + 2*lmMst1 - 2*lmMt)*pow2(Mst1))*pow2(Mst2) + 4*(3*(4 +
        lmMst1 - lmMt)*pow2(Mgl) + (1 + lmMst1 - lmMt)*pow2(Mst1))*pow4(Mst2)))
        *pow4(Mt) + 5*pow2(Dmst12)*pow2(Mt)*pow2(s2t)*(36*pow2(Mgl)*(30*Dmsqst1
        - 59*pow2(Mst1))*pow2(Mst1) - pow2(MuSUSY)*(274*pow2(Mgl)*pow2(Mst1) -
        90*Dmsqst1*(pow2(Mgl) + pow2(Mst1)) - 41*pow4(Mst1)) + 720*Dmsqst1*
        pow4(Mst1) + 360*pow6(Mst1))))))) + Al4p*threeLoopFlag*(-45*xDmsqst1*
        pow2(Dmsqst1)*pow2(Mst2)*pow2(Mt)*(Tbeta*pow2(Mst1)*(5*Dmst12*s2t*((
        Dmst12*s2t*pow2(Mst1) + 8*Mgl*Mt*pow2(Mst2))*pow2(MuSUSY) - Dmst12*(8*
        Mgl*Mt*pow2(MuSUSY) + s2t*pow2(Mst1)*(36*pow2(Mgl) + 8*pow2(Mst1) +
        pow2(MuSUSY))*pow2(Sbeta))) + 20*Mt*pow2(Sbeta)*(pow2(Dmst12)*(6*(17 +
        6*lmMst1 - 6*lmMt)*Mt*pow2(Mgl) + (-5 - 2*lmMst1 + 2*lmMt)*Mt*pow2(
        Mst1) + 2*Mgl*s2t*(4*(4 + lmMst1 - lmMt)*pow2(Mst1) + pow2(MuSUSY)) -
        8*(25 + 6*lmMst1 - 6*lmMt)*s2t*pow3(Mgl)) + 2*Dmst12*pow2(Mst2)*(-3*(17
        + 6*lmMst1 - 6*lmMt)*Mt*pow2(Mgl) + (-3 - 2*lmMst1 + 2*lmMt)*Mt*pow2(
        Mst1) + Mgl*s2t*(8*(3 + lmMst1 - lmMt)*pow2(Mst1) - pow2(MuSUSY)) + 4*(
        25 + 6*lmMst1 - 6*lmMt)*s2t*pow3(Mgl)) - 4*Mt*(3*(5 + 2*lmMst1 - 2*
        lmMt)*pow2(Mgl) + (1 + lmMst1 - lmMt)*pow2(Mst1))*pow4(Mst2))) +
        MuSUSY*pow2(Sbeta)*(-40*Dmst12*Mt*s2t*pow2(Mst1)*(Dmst12*(-9*pow2(Mgl)
        + pow2(Mst1)) + (9*pow2(Mgl) + 2*pow2(Mst1))*pow2(Mst2)) + 60*Mgl*pow2(
        Dmst12)*pow2(s2t)*pow4(Mst1) + pow2(Mt)*(320*(11 + 3*lmMst1 - 3*lmMt)*
        pow3(Mgl)*pow4(Mst2) - 80*Mgl*pow2(Mst1)*(Dmst12*(17 + 6*lmMst1 - 6*
        lmMt)*(Dmst12 - pow2(Mst2)) - 2*(5 + 2*lmMst1 - 2*lmMt)*pow4(Mst2)))))
        + xDmst12*pow2(Mst1)*pow3(Dmst12)*(MuSUSY*pow2(Sbeta)*(2*s2t*pow2(Mst1)
        *pow2(Mt)*(900*Dmsqst1*Mgl*(-30*Mgl*Mt + 4*s2t*pow2(Mgl) - s2t*pow2(
        Mst1)) + pow2(Mst1)*(4556*Mt*pow2(Mgl) + 320*Mt*pow2(Mst1) - 960*Mgl*
        s2t*pow2(Mst1) - 15*(37 + 90*lmMgl - 90*lmMst1)*s2t*pow3(Mgl))) + 5*Mt*
        pow3(s2t)*pow4(Mst1)*(274*pow2(Mgl)*pow2(Mst1) - 30*Dmsqst1*(3*pow2(
        Mgl) + (3 - 4*shiftst1 + 2*shiftst2)*pow2(Mst1)) + (-41 + 120*shiftst1
        - 60*shiftst2)*pow4(Mst1)) + 80*Mgl*(4*(380 + 147*lmMst1 - 249*lmMt)*
        pow2(Mgl)*pow2(Mst1) - 360*Dmsqst1*(3*(4 + lmMst1 - lmMt)*pow2(Mgl) + (
        -7 - 2*lmMst1 + 2*lmMt)*pow2(Mst1)) + (-16 + 9*lmMst1 - 9*lmMt)*pow4(
        Mst1))*pow4(Mt) + 15*shiftst3*(-24*s2t*pow3(Mt) + 7*Mt*pow2(Mst1)*pow3(
        s2t))*pow6(Mst1)) + Tbeta*(-30*(Dmsqst1 + pow2(Mst1))*pow2(s2t)*pow4(
        Mst1)*(pow2(Mt)*(20*shiftst2*pow2(Mst1)*pow2(Sbeta) + pow2(MuSUSY)*(3*
        shiftst1 + 2*shiftst2 - 2*shiftst2*pow2(Sbeta))) - 5*shiftst2*pow2(s2t)
        *pow2(Sbeta)*pow4(Mst1)) - 3*shiftst3*(-2*pow2(Mt)*pow2(s2t)*(9*pow2(
        MuSUSY)*(-1 + pow2(Sbeta)) + 145*pow2(Mst1)*pow2(Sbeta)) + 5*pow2(
        Sbeta)*(24*pow4(Mt) + pow4(Mst1)*pow4(s2t)))*pow6(Mst1) + 2*pow2(Mt)*
        pow2(MuSUSY)*(6*Mgl*Mt*s2t*(500*Dmsqst1 - 79*pow2(Mst1))*pow2(Mst1) +
        pow2(Mgl)*pow2(Mst1)*(996*pow2(Mt) + 5*(90*Dmsqst1 - 91*pow2(Mst1))*
        pow2(s2t)) - 80*Mt*s2t*(15*Dmsqst1 + (-94 + 45*lmMgl - 45*lmMst1)*pow2(
        Mst1))*pow3(Mgl) - 10*pow2(s2t)*pow6(Mst1)) + pow2(Sbeta)*(10*pow2(s2t)
        *pow4(Mst1)*(-2*Mgl*Mt*s2t*((113 - 90*lmMgl + 90*lmMst1)*pow2(Mgl)*
        pow2(Mst1) + 30*Dmsqst1*(2*pow2(Mgl) + pow2(Mst1)) - 17*pow4(Mst1)) +
        3*shiftst1*(Dmsqst1 + pow2(Mst1))*(pow2(Mt)*(80*pow2(Mst1) + 3*pow2(
        MuSUSY)) - 5*pow2(s2t)*pow4(Mst1))) - 3*s2t*pow3(Mt)*((80*pow2(Mst1)*(
        21*(56 + 11*lmMst1 - 11*lmMt)*pow2(Mst1) + 2*(94 - 45*lmMgl + 45*
        lmMst1)*pow2(MuSUSY))*pow3(Mgl))/3. + Mgl*(-316*pow2(MuSUSY)*pow4(Mst1)
        + 400*Dmsqst1*(5*pow2(Mst1)*pow2(MuSUSY) - 2*pow2(Mgl)*(30*(9 + 2*
        lmMst1 - 2*lmMt)*pow2(Mst1) + pow2(MuSUSY)) + 6*pow4(Mst1)) + 80*(32 +
        5*lmMst1 - 5*lmMt)*pow6(Mst1))) + pow2(Mst1)*(8*(-150*Dmsqst1*(63*(7 +
        2*lmMst1 - 2*lmMt)*pow2(Mgl) - pow2(Mst1)) + pow2(Mgl)*(5*(929 + 210*
        lmMst1 - 210*lmMt)*pow2(Mst1) - 249*pow2(MuSUSY)) + 5*(17 - 24*lmMst1 +
        24*lmMt)*pow4(Mst1))*pow4(Mt) - 10*pow2(Mt)*pow2(s2t)*(-2*pow2(MuSUSY)*
        pow4(Mst1) - pow2(Mgl)*(91*pow2(Mst1)*pow2(MuSUSY) + 1332*pow4(Mst1)) +
        90*(Dmsqst1*pow2(Mgl)*(-6*pow2(Mst1) + pow2(MuSUSY)) + 2*Dmsqst1*pow4(
        Mst1) + pow6(Mst1)))))))) - 45*xMgl*pow4(Mgl)*(16*Tbeta*twoLoopFlag*
        pow2(Mst2)*pow2(Sbeta)*pow4(Mst1)*(pow2(Dmst12) + 2*Dmst12*pow2(Mst2) +
        2*pow4(Mst2))*pow4(Mt) + (4*Al4p*Mt*threeLoopFlag*(xDmst12*pow2(Mst1)*
        pow3(Dmst12)*(-2*Mt*pow2(Mst1)*(-3*(-49 + 17*lmMgl - 17*lmMst1)*Tbeta*
        pow2(MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta)) + 2*Mt*(6*(170 - 27*lmMgl +
        27*lmMst1)*MuSUSY*s2t - (2363 + 324*lmMgl + 564*lmMst1 - 2004*lmMt)*Mt*
        Tbeta)*pow2(Sbeta)) + 72*(-17 + 4*lmMgl - 4*lmMst1)*Tbeta*pow2(MuSUSY)*
        (-1 + pow2(Sbeta))*pow3(Mt) + pow2(Sbeta)*(-900*Dmsqst1*Mt*(Mt*MuSUSY*
        s2t + 6*(4 + lmMst1 - lmMt)*Tbeta*pow2(Mt) - Tbeta*pow2(Mst1)*pow2(s2t)
        ) - 3*((49 - 17*lmMgl + 17*lmMst1)*MuSUSY*s2t + (242 - 27*lmMgl + 27*
        lmMst1)*Mt*Tbeta)*pow2(s2t)*pow4(Mst1))) - 900*(11 + 3*lmMst1 - 3*lmMt)
        *Tbeta*xDmsqst1*pow2(Dmsqst1)*pow2(Sbeta)*pow3(Mt)*pow6(Mst2)))/9.))) +
        pow2(Mst2)*(90*pow2(Mt)*(-16*Al4p*twoLoopFlag*(Tbeta*pow4(Mst1)*(s2t*
        pow2(Dmst12)*(s2t*pow2(Mst1)*((1 + 2*lmMst1)*pow2(Mst1)*(12*pow2(Mst1)
        + pow2(MuSUSY)) - pow2(Mgl)*(12*(2 + lmMst1)*pow2(Mst1) + (3 + 2*
        lmMst1)*pow2(MuSUSY)))*pow2(Sbeta) + Mt*pow2(MuSUSY)*(2*Mgl*pow2(Mst1)
        + (6 - 4*lmMgl + 4*lmMst1)*pow3(Mgl))) + Mt*pow2(Sbeta)*(Dmst12*s2t*(-(
        pow2(MuSUSY)*(4*Mgl*((1 - 2*lmMgl + 2*lmMst1)*pow2(Mgl) + (2 + lmMst1)*
        pow2(Mst1))*pow2(Mst2) + Dmst12*(2*Mgl*pow2(Mst1) + (6 - 4*lmMgl + 4*
        lmMst1)*pow3(Mgl)))) + pow2(Mst1)*(12*Dmst12*((1 + 2*lmMt)*Mgl*pow2(
        Mst1) + (5 + 4*lmMst1 - 2*(lmMgl + lmMt))*pow3(Mgl)) + 24*pow2(Mst2)*(-
        ((1 + 2*lmMt)*Mgl*pow2(Mst1)) + (1 + 5*lmMst1 - lmMgl*(3 + 2*lmMst1 -
        2*lmMt) - 2*(1 + lmMst1)*lmMt + 2*pow2(lmMst1))*pow3(Mgl)))) + 12*Mt*
        pow2(Mst1)*((lmMst1 + lmMt)*pow2(Dmst12)*pow2(Mst1) + 2*Dmst12*((1 +
        lmMt)*pow2(Mgl) - (1 + lmMst1 + lmMt)*pow2(Mst1))*pow2(Mst2) + 2*(2*(1
        + lmMt)*pow2(Mgl) - (2 + 2*lmMst1*(1 + lmMt) + pow2(lmMst1) - 3*pow2(
        lmMt))*pow2(Mst1))*pow4(Mst2)))) + MuSUSY*pow4(Mst1)*(Dmst12*MuSUSY*
        s2t*Tbeta*(Dmst12*s2t*pow2(Mst1)*((3 + 2*lmMst1)*pow2(Mgl) - (1 + 2*
        lmMst1)*pow2(Mst1)) + 4*Mgl*Mt*((1 - 2*lmMgl + 2*lmMst1)*pow2(Mgl) + (2
        + lmMst1)*pow2(Mst1))*pow2(Mst2)) + pow2(Sbeta)*(3*Dmst12*s2t*pow2(
        Mst1)*(8*Mt*(-((2 + lmMst1)*pow2(Mgl)) + (1 + 2*lmMst1)*pow2(Mst1))*
        pow2(Mst2) + 2*Dmst12*(-2*Mt*pow2(Mgl) + 2*(1 - 2*lmMst1)*Mt*pow2(Mst1)
        + (2 + lmMst1)*Mgl*s2t*pow2(Mst1) + (1 - 2*lmMgl + 2*lmMst1)*s2t*pow3(
        Mgl))) + pow2(Mt)*(pow2(Dmst12)*(8*(1 + lmMt)*Mgl*pow2(Mst1) - 16*(-3 +
        lmMgl - 2*lmMst1 + lmMt)*pow3(Mgl)) + 3*Dmst12*pow2(Mst2)*(-8*(1 +
        lmMt)*Mgl*pow2(Mst1) + 8*(2 + 6*lmMst1 - 3*lmMt - 2*lmMst1*lmMt +
        lmMgl*(-3 - 2*lmMst1 + 2*lmMt) + 2*pow2(lmMst1))*pow3(Mgl)) + (-48*(1 +
        lmMt)*Mgl*pow2(Mst1) - 48*(lmMgl*(1 + lmMst1 - lmMt) + lmMst1*(-2 +
        lmMt) + lmMt - pow2(lmMst1))*pow3(Mgl))*pow4(Mst2)))) - x*(3*MuSUSY*
        pow2(Sbeta)*(2*Dmst12*(-1 + 3*lmMgl - 3*lmMst1)*s2t*pow2(Mst1)*(4*
        Dmst12*Mgl*Mt + Dmst12*s2t*pow2(Mst1) - 4*Mgl*Mt*pow2(Mst2))*pow5(Mgl)
        + 8*pow2(Mt)*(pow2(Mst1)*(Dmst12*(lmMst1*(22 - 6*lmMt) - 5*lmMt +
        lmMgl*(-17 - 6*lmMst1 + 6*lmMt) + 6*pow2(lmMst1))*(Dmst12 - pow2(Mst2))
        + 2*(2 + lmMgl*(5 + 2*lmMst1 - 2*lmMt) + 2*lmMst1*(-3 + lmMt) + lmMt -
        2*pow2(lmMst1))*pow4(Mst2))*pow5(Mgl) + (9 + lmMgl*(20 + 6*lmMst1 - 6*
        lmMt) + 2*lmMt + lmMst1*(-22 + 6*lmMt) - 6*pow2(lmMst1))*pow4(Mst2)*
        pow7(Mgl))) + Tbeta*(4*Dmst12*(1 - 3*lmMgl + 3*lmMst1)*s2t*pow2(Mst1)*(
        Dmst12*Mt*pow2(MuSUSY) - Mt*pow2(Mst2)*pow2(MuSUSY) + 3*Dmst12*Mgl*s2t*
        pow2(Mst1)*pow2(Sbeta))*pow5(Mgl) + 4*Mt*pow2(Sbeta)*(Dmst12*s2t*pow2(
        Mst1)*(pow2(Mst2)*(6*(8 - 25*lmMst1 + lmMgl*(23 + 6*lmMst1 - 6*lmMt) +
        2*lmMt + 6*lmMst1*lmMt - 6*pow2(lmMst1))*pow2(Mgl) + 6*(3 - 14*lmMst1 +
        4*lmMgl*(3 + lmMst1 - lmMt) + 2*lmMt + 4*lmMst1*lmMt - 4*pow2(lmMst1))*
        pow2(Mst1) + (1 - 3*lmMgl + 3*lmMst1)*pow2(MuSUSY)) - Dmst12*(6*(8 -
        25*lmMst1 + lmMgl*(23 + 6*lmMst1 - 6*lmMt) + 2*lmMt + 6*lmMst1*lmMt -
        6*pow2(lmMst1))*pow2(Mgl) + 3*(11 + lmMst1*(22 - 4*lmMt) - 4*lmMgl*(4 +
        lmMst1 - lmMt) - 6*lmMt + 4*pow2(lmMst1))*pow2(Mst1) + (1 - 3*lmMgl +
        3*lmMst1)*pow2(MuSUSY)))*pow5(Mgl) - 6*Mt*(pow2(Mst1)*(Dmst12*(lmMst1*(
        22 - 6*lmMt) - 5*lmMt + lmMgl*(-17 - 6*lmMst1 + 6*lmMt) + 6*pow2(
        lmMst1))*(Dmst12 - pow2(Mst2)) + 2*(2 + lmMgl*(5 + 2*lmMst1 - 2*lmMt) +
        2*lmMst1*(-3 + lmMt) + lmMt - 2*pow2(lmMst1))*pow4(Mst2))*pow6(Mgl) + (
        9 + lmMgl*(20 + 6*lmMst1 - 6*lmMt) + 2*lmMt + lmMst1*(-22 + 6*lmMt) -
        6*pow2(lmMst1))*pow4(Mst2)*pow8(Mgl)))))) + 3*oneLoopFlag*(12*Dmst12*
        Mt*MuSUSY*s2t*(Dmst12 - 2*pow2(Mst2))*pow2(Sbeta) - Tbeta*pow2(Dmst12)*
        pow2(s2t)*(pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 12*pow2(Mst1)*pow2(Sbeta))
        - 12*Tbeta*pow2(Mt)*pow2(Sbeta)*(pow2(Dmst12) - 2*Dmst12*pow2(Mst2) +
        4*(-lmMst1 + lmMt)*pow4(Mst2)))*pow8(Mst1)) + 60*threeLoopFlag*pow2(
        Al4p)*(5*Mt*xDmsqst1*pow2(Dmsqst1)*(18*MuSUSY*pow2(Sbeta)*(Dmst12*pow2(
        Mst1)*(2*(16*(1 - 19*lmMst1 + 19*lmMt)*Mgl*Mt + 192*s2t*pow2(Mgl) +
        s2t*(47 - 16*shiftst2)*pow2(Mst1))*pow2(Mst2)*pow2(Mt) - Dmst12*(s2t*(
        384*pow2(Mgl) - 45*pow2(Mst1))*pow2(Mt) + 32*Mgl*Mt*((1 - 19*lmMst1 +
        19*lmMt)*pow2(Mt) + 2*pow2(Mst1)*pow2(s2t)) + 4*(shiftst1 - shiftst2)*
        pow3(s2t)*pow4(Mst1))) - 32*pow2(Mt)*(4*(-1 + 3*lmMst1 - 3*lmMt)*Mgl*
        Mt*pow2(Mst1) + 4*Mt*(4 + lmMgl*(1 + lmMst1 - lmMt) - 8*lmMt + lmMst1*(
        7 + lmMt) - pow2(lmMst1))*pow3(Mgl) + s2t*(-shiftst1 + shiftst2)*pow4(
        Mst1))*pow4(Mst2)) + Mt*Tbeta*pow2(Mst1)*(3*Dmst12*s2t*(256*Dmst12*Mgl*
        Mt*pow2(MuSUSY) - (21*Dmst12*s2t*pow2(Mst1) + 256*Mgl*Mt*pow2(Mst2))*
        pow2(MuSUSY) + 3*Dmst12*s2t*pow2(Mst1)*(384*pow2(Mgl) + 94*pow2(Mst1) +
        7*pow2(MuSUSY))*pow2(Sbeta)) + 24*pow2(Mst1)*(pow2(Dmst12)*pow2(s2t)*(-
        ((shiftst1 + shiftst2)*pow2(MuSUSY)*(-1 + pow2(Sbeta))) + 6*(-3*
        shiftst1 + shiftst2)*pow2(Mst1)*pow2(Sbeta)) + 2*Dmst12*pow2(Mst2)*((
        shiftst1 - shiftst2)*pow2(MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta)) + 12*(
        shiftst2*pow2(Mt) + (shiftst1 - shiftst2)*pow2(Mst1)*pow2(s2t))*pow2(
        Sbeta)) + 24*(shiftst1 + shiftst2)*pow2(Mt)*pow2(Sbeta)*pow4(Mst2)) +
        Mt*pow2(Sbeta)*(96*Dmst12*s2t*(-8*Mgl*(Dmst12 - pow2(Mst2))*pow2(
        MuSUSY) - 12*pow2(Mst2)*((1 + 6*lmMst1 - 6*lmMt)*Mgl*pow2(Mst1) + (16 +
        13*lmMst1 + lmMgl*(3 + 2*lmMst1 - 2*lmMt) + 2*(-8 + lmMst1)*lmMt - 2*
        pow2(lmMst1))*pow3(Mgl)) + 3*Dmst12*(-((13 + 14*lmMst1 - 14*lmMt)*Mgl*
        pow2(Mst1)) + 4*(16 + 13*lmMst1 + lmMgl*(3 + 2*lmMst1 - 2*lmMt) + 2*(-8
        + lmMst1)*lmMt - 2*pow2(lmMst1))*pow3(Mgl))) + Mt*(pow2(Dmst12)*(864*(1
        - 19*lmMst1 + 19*lmMt)*pow2(Mgl) + 19*(-11 + 48*lmMst1 - 48*lmMt)*pow2(
        Mst1)) - 72*Dmst12*(12*(1 - 19*lmMst1 + 19*lmMt)*pow2(Mgl) + 7*(5 - 4*
        lmMst1 + 4*lmMt)*pow2(Mst1))*pow2(Mst2) - 16*(216*(1 - 3*lmMst1 + 3*
        lmMt)*pow2(Mgl) + (191 - 120*lmMst1 + 12*(10 + 3*lmMst1)*lmMt - 18*(
        pow2(lmMst1) + pow2(lmMt)))*pow2(Mst1))*pow4(Mst2))))) + 4*pow2(Mst1)*(
        Tbeta*(pow2(Mt)*(-(pow2(MuSUSY)*(-40*Dmst12*Mgl*Mt*(-4*Mgl*Mt*(-16 + 9*
        lmMgl - 41*lmMst1 + 12*lmMt - 4*pow2(lmMst1))*pow2(Mst1) + 4*s2t*pow2(
        Mgl)*(30*Dmsqst1*(-3 + lmMgl - lmMst1) + (-11 + 154*lmMgl - 148*lmMst1
        - 20*lmMgl*lmMst1 + 27*pow2(lmMgl) + pow2(lmMst1))*pow2(Mst1)) - s2t*
        pow2(Mst1)*(180*Dmsqst1 + (215 + 118*lmMst1 + 41*pow2(lmMst1))*pow2(
        Mst1)))*pow2(Mst2) + pow2(Dmst12)*(-2*pow2(Mgl)*pow2(Mst1)*(4*(163 +
        15*lmMgl - 27*lmMst1 + 20*lmMt - 8*pow2(lmMst1))*pow2(Mt) - 5*(360*
        Dmsqst1 + (200 + 18*lmMgl + 265*lmMst1 + 91*pow2(lmMst1))*pow2(Mst1))*
        pow2(s2t)) + Mt*s2t*(40*Mgl*pow2(Mst1)*(210*Dmsqst1 + (149 + 41*lmMst1)
        *pow2(Mst1)) - 80*(60*Dmsqst1*(3 - lmMgl + lmMst1) + (-249 - 176*lmMst1
        - 4*lmMgl*(-41 + 5*lmMst1) + 27*pow2(lmMgl) + pow2(lmMst1))*pow2(Mst1))
        *pow3(Mgl)) + 5*(15*Dmsqst1*(3 - 8*lmMst1) - 2*(16 + 115*lmMst1 + 41*
        pow2(lmMst1))*pow2(Mst1))*pow2(s2t)*pow4(Mst1)) + 80*(5 + 18*lmMgl -
        74*lmMst1 + 24*lmMt - 8*pow2(lmMst1))*pow2(Mgl)*pow2(Mst1)*pow2(Mt)*
        pow4(Mst2)))/5. + 30*(Dmsqst1*(3 - 2*lmMst1) + (1 - 2*lmMst1)*pow2(
        Mst1))*pow4(Mst1)*(pow2(Dmst12)*pow2(s2t)*(-((shiftst1 + shiftst2)*
        pow2(MuSUSY)*(-1 + pow2(Sbeta))) + 6*(-3*shiftst1 + shiftst2)*pow2(
        Mst1)*pow2(Sbeta)) + 2*Dmst12*pow2(Mst2)*((shiftst1 - shiftst2)*pow2(
        MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta)) + 12*(shiftst2*pow2(Mt) + (
        shiftst1 - shiftst2)*pow2(Mst1)*pow2(s2t))*pow2(Sbeta)) + 24*(shiftst1
        + shiftst2)*pow2(Mt)*pow2(Sbeta)*pow4(Mst2))) + 18*shiftst3*((Dmst12*
        pow2(Mt)*pow2(s2t)*(-2*(-1 + 2*lmMst1)*pow2(Mst2)*(pow2(MuSUSY)*(-1 +
        pow2(Sbeta)) + 12*pow2(Mst1)*pow2(Sbeta)) + Dmst12*((-7 + 6*lmMst1)*
        pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 6*(-15 + 14*lmMst1)*pow2(Mst1)*pow2(
        Sbeta))))/6. + 4*pow2(Sbeta)*(-2*(-2 + lmMst1)*pow2(Dmst12) + Dmst12*(-
        3 + 2*lmMst1)*pow2(Mst2) + (1 - 2*lmMst1)*pow4(Mst2))*pow4(Mt))*pow6(
        Mst1) + pow2(Sbeta)*(-8*Dmst12*Mgl*s2t*pow3(Mt)*((-30*Dmsqst1*(Dmst12*(
        4*(-3 + lmMgl - lmMst1)*pow2(Mgl) + 7*pow2(Mst1)) + 2*((6 - 2*lmMgl +
        2*lmMst1)*pow2(Mgl) + 3*pow2(Mst1))*pow2(Mst2)) - pow2(Mst1)*(Dmst12*(-
        2*(-249 + 164*lmMgl - 176*lmMst1 - 20*lmMgl*lmMst1 + 27*pow2(lmMgl) +
        pow2(lmMst1))*pow2(Mgl) + (149 + 41*lmMst1)*pow2(Mst1)) + (-4*(-11 +
        154*lmMgl - 148*lmMst1 - 20*lmMgl*lmMst1 + 27*pow2(lmMgl) + pow2(
        lmMst1))*pow2(Mgl) + (215 + 118*lmMst1 + 41*pow2(lmMst1))*pow2(Mst1))*
        pow2(Mst2)))*pow2(MuSUSY) + 2*pow2(Mst1)*(pow2(Mst2)*(180*Dmsqst1*((19
        + 13*lmMst1 - 10*lmMt - 2*lmMst1*lmMt + lmMgl*(-3 - 2*lmMst1 + 2*lmMt)
        + 2*pow2(lmMst1))*pow2(Mgl) + 4*(1 + lmMst1 - lmMt)*pow2(Mst1)) + pow2(
        Mgl)*pow2(Mst1)*(2086 - 996*lmMt - 81*(9 + 4*lmMst1 - 6*lmMt)*pow2(
        lmMgl) + (663 + 366*lmMt)*pow2(lmMst1) - 144*pow2(lmMt) - 6*lmMst1*(-
        521 + 74*lmMt + 24*pow2(lmMt)) + 6*lmMgl*(-238 + lmMst1*(73 - 142*lmMt)
        + 52*lmMt + 91*pow2(lmMst1) + 24*pow2(lmMt)) - 54*pow3(lmMgl) - 168*
        pow3(lmMst1)) - 2*(-407 + 276*lmMt + 3*lmMst1*(-65 + 2*lmMt) + 45*pow2(
        lmMst1) + 72*pow2(lmMt))*pow4(Mst1)) + Dmst12*(3*(1091 + 611*lmMst1 +
        lmMgl*(-17 + 122*lmMst1 - 82*lmMt) + (-374 + 62*lmMst1)*lmMt - 54*pow2(
        lmMgl) - 24*(pow2(lmMst1) + pow2(lmMt)))*pow2(Mgl)*pow2(Mst1) + 90*
        Dmsqst1*(-2*(-24 + lmMgl - 15*lmMst1 + 14*lmMt)*pow2(Mgl) + (5 + 6*
        lmMst1 - 6*lmMt)*pow2(Mst1)) + (148 - 90*lmMt + lmMst1*(75 + 6*lmMt) +
        45*pow2(lmMst1) + 72*pow2(lmMt))*pow4(Mst1)))) + pow2(Mst1)*((2*(-20*
        Dmst12*pow2(Mst2)*(-(pow2(Mgl)*((3049 - 1908*lmMt - 54*lmMgl*(1 + 12*
        lmMt) + 6*lmMst1*(215 + 49*lmMt) + 324*pow2(lmMgl) - 195*pow2(lmMst1) -
        144*pow2(lmMt))*pow2(Mst1) + 4*(-16 + 9*lmMgl - 41*lmMst1 + 12*lmMt -
        4*pow2(lmMst1))*pow2(MuSUSY))) + (196 - 54*lmMt + 12*lmMst1*(-24 + 5*
        lmMt) - 153*(pow2(lmMst1) + pow2(lmMt)))*pow4(Mst1)) - pow2(Dmst12)*(-
        4*pow2(Mgl)*((6635 - 9965*lmMt + 5*lmMst1*(1977 + 64*lmMt) - 160*pow2(
        lmMst1))*pow2(Mst1) + (-163 - 15*lmMgl + 27*lmMst1 - 20*lmMt + 8*pow2(
        lmMst1))*pow2(MuSUSY)) + 5*(131 + 1500*lmMt - 60*lmMst1*(29 + 2*lmMt) +
        306*(pow2(lmMst1) + pow2(lmMt)))*pow4(Mst1)) - 10*(-4*pow2(Mgl)*((1909
        + 600*lmMst1 + 6*(-187 + 49*lmMst1)*lmMt - 54*lmMgl*(1 + 12*lmMt) +
        324*pow2(lmMgl) - 195*pow2(lmMst1) - 144*pow2(lmMt))*pow2(Mst1) + (5 +
        18*lmMgl - 74*lmMst1 + 24*lmMt - 8*pow2(lmMst1))*pow2(MuSUSY)) + 3*(
        1375 - 64*lmMt + 24*(2 + 3*lmMt)*pow2(lmMgl) + 4*(37 + 25*lmMt)*pow2(
        lmMst1) + 120*pow2(lmMt) - 8*lmMgl*(26 + 12*lmMt + 9*pow2(lmMt)) - 8*
        lmMst1*(44 + 64*lmMt + 33*pow2(lmMt)) - 24*pow3(lmMgl) - 88*pow3(
        lmMst1) + 276*pow3(lmMt))*pow4(Mst1))*pow4(Mst2) + 25*Dmsqst1*(9*pow2(
        Dmst12)*(12*(19 + 33*lmMst1 - 33*lmMt)*pow2(Mgl) + (-3 + 28*lmMst1 -
        28*lmMt)*pow2(Mst1)) + 18*Dmst12*(24*(3 + 7*lmMst1 - 7*lmMt)*pow2(Mgl)
        + (-19 + 32*lmMst1 - 24*lmMt)*pow2(Mst1))*pow2(Mst2) + 8*(108*(1 + 2*
        lmMst1 - 2*lmMt)*pow2(Mgl) - (92 + 66*lmMt - 6*lmMst1*(17 + 3*lmMt) +
        9*(pow2(lmMst1) + pow2(lmMt)))*pow2(Mst1))*pow4(Mst2)))*pow4(Mt))/5. +
        2*pow2(Dmst12)*pow2(Mt)*pow2(s2t)*(6*pow2(Mgl)*pow2(Mst1)*(540*Dmsqst1
        + (525 - 18*lmMgl + 504*lmMst1 - 48*lmMt + 107*pow2(lmMst1))*pow2(Mst1)
        ) + 90*Dmsqst1*(5 - 4*lmMst1)*pow4(Mst1) + (pow2(MuSUSY)*(2*(200 + 18*
        lmMgl + 265*lmMst1 + 91*pow2(lmMst1))*pow2(Mgl)*pow2(Mst1) + 15*
        Dmsqst1*(48*pow2(Mgl) + (3 - 8*lmMst1)*pow2(Mst1)) - 2*(16 + 115*lmMst1
        + 41*pow2(lmMst1))*pow4(Mst1)))/2. - 12*(-6 + 79*lmMst1 + 41*pow2(
        lmMst1))*pow6(Mst1))))) + Mt*MuSUSY*pow2(Sbeta)*(3*Dmst12*Mt*s2t*pow2(
        Mst1)*(Mt*(15*Dmsqst1*(Dmst12*(336*pow2(Mgl) + 31*pow2(Mst1)) + 8*(36*
        pow2(Mgl) + (5 - 4*lmMst1)*pow2(Mst1))*pow2(Mst2)) + 4*pow2(Mst1)*(
        Dmst12*(2*(522 + 107*lmMst1)*pow2(Mgl) + (-95 - 66*lmMst1 + 82*pow2(
        lmMst1))*pow2(Mst1)) + 2*((501 - 18*lmMgl + 504*lmMst1 - 48*lmMt + 107*
        pow2(lmMst1))*pow2(Mgl) - 2*(-6 + 79*lmMst1 + 41*pow2(lmMst1))*pow2(
        Mst1))*pow2(Mst2))) - 4*Dmst12*Mgl*s2t*(-4*(-11 + lmMgl*(154 - 20*
        lmMst1) - 148*lmMst1 + 27*pow2(lmMgl) + pow2(lmMst1))*pow2(Mgl)*pow2(
        Mst1) + 60*Dmsqst1*(2*(3 - lmMgl + lmMst1)*pow2(Mgl) + 3*pow2(Mst1)) +
        (215 + 118*lmMst1 + 41*pow2(lmMst1))*pow4(Mst1))) + 90*s2t*(Dmsqst1*(-3
        + 2*lmMst1) + (-1 + 2*lmMst1)*pow2(Mst1))*pow4(Mst1)*(8*Dmst12*
        shiftst2*pow2(Mst2)*pow2(Mt) + (shiftst1 - shiftst2)*pow2(Dmst12)*pow2(
        Mst1)*pow2(s2t) + 8*(-shiftst1 + shiftst2)*pow2(Mt)*pow4(Mst2)) + 16*
        Mgl*pow3(Mt)*(30*Dmsqst1*(pow2(Dmst12)*(6*(31 + 27*lmMst1 - 24*lmMt -
        2*lmMst1*lmMt + lmMgl*(-3 - 2*lmMst1 + 2*lmMt) + 2*pow2(lmMst1))*pow2(
        Mgl) + (-13 - 19*lmMst1 + 19*lmMt)*pow2(Mst1)) - 6*Dmst12*((31 + 27*
        lmMst1 - 24*lmMt - 2*lmMst1*lmMt + lmMgl*(-3 - 2*lmMst1 + 2*lmMt) + 2*
        pow2(lmMst1))*pow2(Mgl) + (3 + 7*lmMst1 - 7*lmMt)*pow2(Mst1))*pow2(
        Mst2) - 12*((7 + 6*lmMst1 - 5*lmMt - lmMst1*lmMt + lmMgl*(-1 - lmMst1 +
        lmMt) + pow2(lmMst1))*pow2(Mgl) + (1 + 2*lmMst1 - 2*lmMt)*pow2(Mst1))*
        pow4(Mst2)) - pow2(Mst1)*(pow2(Dmst12)*(-((-3749 - 2374*lmMst1 + 1770*
        lmMt - 117*lmMst1*lmMt + lmMgl*(96 - 237*lmMst1 + 157*lmMt) + 108*pow2(
        lmMgl) + 41*pow2(lmMst1) + 48*pow2(lmMt))*pow2(Mgl)) + (lmMst1*(297 +
        2*lmMt) + 15*pow2(lmMst1) + 3*(47 - 79*lmMt + 8*pow2(lmMt)))*pow2(Mst1)
        ) + Dmst12*pow2(Mst2)*(-((-605 + 6*lmMst1*(-77 + lmMt) + 714*lmMt + 45*
        pow2(lmMst1) + 72*pow2(lmMt))*pow2(Mst1)) + 2*pow2(Mgl)*(2096 - 1167*
        lmMt - 81*(4 + 2*lmMst1 - 3*lmMt)*pow2(lmMgl) + 3*(129 + 61*lmMt)*pow2(
        lmMst1) - 108*pow2(lmMt) - 3*lmMst1*(-793 + 67*lmMt + 24*pow2(lmMt)) +
        3*lmMgl*(-229 + lmMst1*(63 - 142*lmMt) + 35*lmMt + 91*pow2(lmMst1) +
        24*pow2(lmMt)) - 27*pow3(lmMgl) - 84*pow3(lmMst1))) + 2*(-((-233 + 6*
        lmMst1*(-32 + lmMt) + 348*lmMt + 45*pow2(lmMst1) + 72*pow2(lmMt))*pow2(
        Mst1)) + pow2(Mgl)*(758 - 570*lmMt - 81*(3 + 2*lmMst1 - 3*lmMt)*pow2(
        lmMgl) + 3*(141 + 61*lmMt)*pow2(lmMst1) - 72*pow2(lmMt) - 6*lmMst1*(-
        236 + 49*lmMt + 12*pow2(lmMt)) + 3*lmMgl*(-200 + 2*lmMst1 + 76*lmMt -
        142*lmMst1*lmMt + 91*pow2(lmMst1) + 24*pow2(lmMt)) - 27*pow3(lmMgl) -
        84*pow3(lmMst1)))*pow4(Mst2))) + shiftst3*(72*s2t*pow2(Mt)*(-2*(-2 +
        lmMst1)*pow2(Dmst12) + Dmst12*(-3 + 2*lmMst1)*pow2(Mst2) + (1 - 2*
        lmMst1)*pow4(Mst2))*pow6(Mst1) + 9*(-1 + 2*lmMst1)*pow2(Dmst12)*pow3(
        s2t)*pow8(Mst1)))))))/(12960.*Tbeta*pow2(Sbeta)*pow6(Mst2)*pow8(Mst1));
}

/**
 * @return The off-diagonal (1, 2) = (2, 1) matrix element of the Higgs mass matrix as a double for the hierarchy 'H9'
 */
double H9::getS12() const {
   return (MuSUSY*(Mt*pow2(Mst1)*((20*Al4p*xMgl*pow4(Mgl)*(xDmst12*pow3(Dmst12)*(
        216*s2t*twoLoopFlag*pow2(Mst1)*(4*(1 - 3*lmMgl + 3*lmMst1)*Mt*MuSUSY*
        s2t + 96*(-1 + lmMgl - lmMst1)*Tbeta*pow2(Mt) + (-1 + 3*lmMgl - 3*
        lmMst1)*Tbeta*pow2(Mst1)*pow2(s2t)) + Al4p*threeLoopFlag*(1728*s2t*
        Tbeta*(15*Dmsqst1*(29 - 10*lmMgl + 10*lmMst1) + (-1797 - 1711*lmMst1 +
        lmMgl*(1315 - 198*lmMst1 - 96*lmMt) + 144*lmMt + 96*lmMst1*lmMt + 333*
        pow2(lmMgl) - 103*pow2(lmMst1))*pow2(Mst1))*pow2(Mt) + 4*Mt*MuSUSY*(
        41653 + 11784*lmMst1 + 48*lmMgl*(-313 + 180*lmMst1) - 12636*pow2(lmMgl)
        + 540*pow2(lmMst1))*pow2(Mst1)*pow2(s2t) + 1152*MuSUSY*(262 + lmMst1*(
        335 - 72*lmMt) - 96*lmMt + lmMgl*(-199 - 120*lmMst1 + 72*lmMt) + 120*
        pow2(lmMst1))*pow3(Mt) + Tbeta*(-41653 - 11784*lmMst1 - 48*lmMgl*(-313
        + 180*lmMst1) + 12636*pow2(lmMgl) - 540*pow2(lmMst1))*pow3(s2t)*pow4(
        Mst1))) + Mt*pow2(Mst2)*(-432*Dmst12*s2t*twoLoopFlag*pow2(Mst1)*(
        Dmst12*(1 - 3*lmMgl + 3*lmMst1)*MuSUSY*s2t + 6*Dmst12*(-3 + 2*lmMgl -
        2*lmMst1)*Mt*Tbeta + 12*(-1 + 2*lmMgl - 2*lmMst1)*Mt*Tbeta*pow2(Mst2))
        + 2*Al4p*threeLoopFlag*(144*Dmst12*Mt*(4*Mt*MuSUSY*(262 + lmMst1*(335 -
        72*lmMt) - 96*lmMt + lmMgl*(-199 - 120*lmMst1 + 72*lmMt) + 120*pow2(
        lmMst1)) + 3*s2t*Tbeta*(30*Dmsqst1*(29 - 10*lmMgl + 10*lmMst1) + (633 +
        1458*lmMst1 - 96*(1 + lmMst1)*lmMt + 2*lmMgl*(-632 + 99*lmMst1 + 48*
        lmMt) - 333*pow2(lmMgl) + 103*pow2(lmMst1))*pow2(Mst1)))*pow2(Mst2) -
        pow2(Dmst12)*(-216*Mt*s2t*Tbeta*(60*Dmsqst1*(-29 + 10*lmMgl - 10*
        lmMst1) + (2961 + 1964*lmMst1 - 192*lmMt - 96*lmMst1*lmMt + 2*lmMgl*(-
        683 + 99*lmMst1 + 48*lmMt) - 333*pow2(lmMgl) + 103*pow2(lmMst1))*pow2(
        Mst1)) + 576*MuSUSY*(262 + 335*lmMst1 - 96*lmMt - 72*lmMst1*lmMt +
        lmMgl*(-199 - 120*lmMst1 + 72*lmMt) + 120*pow2(lmMst1))*pow2(Mt) +
        MuSUSY*(41653 + 11784*lmMst1 + 48*lmMgl*(-313 + 180*lmMst1) - 12636*
        pow2(lmMgl) + 540*pow2(lmMst1))*pow2(Mst1)*pow2(s2t)) + 1152*MuSUSY*(59
        + lmMst1*(85 - 24*lmMt) - 24*lmMt + lmMgl*(-53 - 40*lmMst1 + 24*lmMt) +
        40*pow2(lmMst1))*pow2(Mt)*pow4(Mst2)))))/9. - (xDmst12*pow3(Dmst12)*(-
        270*oneLoopFlag*s2t*(-2*Mt*MuSUSY*s2t - 8*Tbeta*pow2(Mt) + Tbeta*pow2(
        Mst1)*pow2(s2t))*pow6(Mst1) - 288*Al4p*twoLoopFlag*(5*Mt*pow2(Mst1)*(
        32*(1 - lmMgl + lmMst1)*Mt*MuSUSY*s2t + 8*Tbeta*(18 + lmMst1*(26 - 6*
        lmMt) - 13*lmMt + lmMgl*(-13 - 6*lmMst1 + 6*lmMt) + 6*pow2(lmMst1))*
        pow2(Mt) - 3*(1 + 2*lmMgl - 2*lmMst1)*Tbeta*pow2(Mst1)*pow2(s2t))*pow3(
        Mgl) + 5*s2t*pow2(Mgl)*((1 + 2*lmMst1)*Mt*MuSUSY*s2t - 4*Tbeta*pow2(Mt)
        - (3 + 2*lmMst1)*Tbeta*pow2(Mst1)*pow2(s2t))*pow4(Mst1) + Mgl*Mt*(2*(9
        + 2*lmMst1)*Mt*MuSUSY*s2t + 20*(1 + lmMt)*Tbeta*pow2(Mt) + 15*(3 + 2*
        lmMst1)*Tbeta*pow2(Mst1)*pow2(s2t))*pow4(Mst1) - 20*Mt*x*(2*(1 - 3*
        lmMgl + 3*lmMst1)*Mt*MuSUSY*s2t - 6*Tbeta*(-22*lmMst1 + lmMgl*(17 + 6*
        lmMst1 - 6*lmMt) + 5*lmMt + 6*lmMst1*lmMt - 6*pow2(lmMst1))*pow2(Mt) +
        3*(-1 + 3*lmMgl - 3*lmMst1)*Tbeta*pow2(Mst1)*pow2(s2t))*pow5(Mgl) +
        120*(1 - 3*lmMgl + 3*lmMst1)*s2t*Tbeta*x*pow2(Mt)*pow6(Mgl) - 5*s2t*(4*
        lmMst1*Mt*MuSUSY*s2t + 2*(-9 + 8*lmMst1)*Tbeta*pow2(Mt) - (1 + 2*
        lmMst1)*Tbeta*pow2(Mst1)*pow2(s2t))*pow6(Mst1)) + threeLoopFlag*pow2(
        Al4p)*(64*Mgl*Mt*(900*xDmsqst1*pow2(Dmsqst1)*(8*Mt*MuSUSY*s2t + 3*(-1 +
        19*lmMst1 - 19*lmMt)*Tbeta*pow2(Mt) - 12*Tbeta*pow2(Mst1)*pow2(s2t)) -
        450*Dmsqst1*pow2(Mst1)*(80*Mt*MuSUSY*s2t + 16*(11 + 20*lmMst1 - 20*
        lmMt)*Tbeta*pow2(Mt) - 3*Tbeta*pow2(Mst1)*pow2(s2t)) + 30*pow2(Mgl)*(
        120*Dmsqst1*(2*(3 - lmMgl + lmMst1)*Mt*MuSUSY*s2t + 3*Tbeta*(31 + 27*
        lmMst1 - 24*lmMt - 2*lmMst1*lmMt + lmMgl*(-3 - 2*lmMst1 + 2*lmMt) + 2*
        pow2(lmMst1))*pow2(Mt) + 3*(-3 + lmMgl - lmMst1)*Tbeta*pow2(Mst1)*pow2(
        s2t)) - pow2(Mst1)*(-16*Mt*MuSUSY*s2t*(-130 + lmMgl*(159 - 20*lmMst1) -
        162*lmMst1 + 27*pow2(lmMgl) + pow2(lmMst1)) + Tbeta*(-3*(227 + 4*lmMgl*
        (36 - 5*lmMst1) - 120*lmMst1 + 27*pow2(lmMgl) + pow2(lmMst1))*pow2(
        Mst1)*pow2(s2t) + 4*pow2(Mt)*(5845 - 2937*lmMt - 27*(16 + lmMgl + 6*
        lmMst1 - 9*lmMt)*pow2(lmMgl) + (346 + 183*lmMt)*pow2(lmMst1) + lmMst1*(
        4753 - 84*lmMt - 72*pow2(lmMt)) - 156*pow2(lmMt) + lmMgl*(-783 - 426*
        lmMst1*(-1 + lmMt) - 52*lmMt + 273*pow2(lmMst1) + 72*pow2(lmMt)) - 84*
        pow3(lmMst1))))) + (-2*Mt*MuSUSY*s2t*(-365 + 1049*lmMst1 + 123*pow2(
        lmMst1)) - 30*Tbeta*(-233 + 47*lmMt + lmMst1*(-3 + 2*lmMt) + 15*pow2(
        lmMst1) + 24*pow2(lmMt))*pow2(Mt) - 45*Tbeta*(66 + 77*lmMst1 + 41*pow2(
        lmMst1))*pow2(Mst1)*pow2(s2t))*pow4(Mst1)) - 3*s2t*pow2(Mst1)*(300*
        xDmsqst1*pow2(Dmsqst1)*(4*Mt*MuSUSY*s2t*(21 - 16*shiftst2) + 48*(-23 +
        4*shiftst2)*Tbeta*pow2(Mt) + (-21 + 32*shiftst1 - 16*shiftst2)*Tbeta*
        pow2(Mst1)*pow2(s2t)) + 15*Dmsqst1*pow2(Mst1)*(8*Mt*MuSUSY*s2t*(5 - 72*
        shiftst1 - 48*shiftst2 + 16*lmMst1*(-5 + 3*shiftst1 + 2*shiftst2)) +
        Tbeta*(1245*pow2(Mt) - 80*(3 - 24*shiftst1 + 12*shiftst2 - 8*lmMst1*(1
        - 2*shiftst1 + shiftst2))*pow2(Mst1)*pow2(s2t))) + 8*(-2*Mt*MuSUSY*s2t*
        (-605 + 180*shiftst1 + 120*shiftst2 + 348*shiftst3 - 8*lmMst1*(-110 +
        45*shiftst1 + 30*shiftst2 + 27*shiftst3) + 820*pow2(lmMst1)) - 5*Tbeta*
        (-2615 + 672*shiftst3 - 8*lmMst1*(325 + 36*shiftst3) + 1312*pow2(
        lmMst1))*pow2(Mt) + 10*Tbeta*(32 + 120*shiftst1 - 60*shiftst2 + lmMst1*
        (230 - 240*shiftst1 + 120*shiftst2 - 42*shiftst3) + 39*shiftst3 + 82*
        pow2(lmMst1))*pow2(Mst1)*pow2(s2t))*pow4(Mst1)) + 96*pow2(Mgl)*(16*Mt*
        pow2(Mst1)*(6750*Dmsqst1*Mt*s2t*Tbeta + MuSUSY*(323 - 75*lmMgl + 383*
        lmMst1 - 100*lmMt + 32*pow2(lmMst1))*pow2(Mt) - 450*Dmsqst1*MuSUSY*
        pow2(s2t)) + s2t*(-21600*Tbeta*xDmsqst1*pow2(Dmsqst1)*pow2(Mt) + (-5*
        Mt*MuSUSY*s2t*(-425 + 18*lmMgl + 83*lmMst1 + 91*pow2(lmMst1)) - 2*
        Tbeta*(4113 + 90*lmMgl - 1792*lmMst1 + 120*lmMt - 48*pow2(lmMst1))*
        pow2(Mt) + 1800*Dmsqst1*Tbeta*pow2(s2t))*pow4(Mst1)) + 5*Tbeta*(200 +
        18*lmMgl + 265*lmMst1 + 91*pow2(lmMst1))*pow3(s2t)*pow6(Mst1)))))/3.) -
        32*Al4p*z2*(Al4p*threeLoopFlag*(xDmst12*pow2(Mst1)*pow3(Dmst12)*(4*
        MuSUSY*pow2(Mt)*(6*Mgl*Mt*s2t*(500*Dmsqst1 - 79*pow2(Mst1))*pow2(Mst1)
        + pow2(Mgl)*pow2(Mst1)*(996*pow2(Mt) + 5*(90*Dmsqst1 - 91*pow2(Mst1))*
        pow2(s2t)) - 80*Mt*s2t*(15*Dmsqst1 + (-94 + 45*lmMgl - 45*lmMst1)*pow2(
        Mst1))*pow3(Mgl) - (15*Dmsqst1*(3*shiftst1 + 2*shiftst2) + (10 + 45*
        shiftst1 + 30*shiftst2 + 27*shiftst3)*pow2(Mst1))*pow2(s2t)*pow4(Mst1))
        + Tbeta*(2*s2t*pow2(Mst1)*pow2(Mt)*(900*Dmsqst1*Mgl*(-30*Mgl*Mt + 4*
        s2t*pow2(Mgl) - s2t*pow2(Mst1)) + pow2(Mst1)*(4556*Mt*pow2(Mgl) + 320*
        Mt*pow2(Mst1) - 960*Mgl*s2t*pow2(Mst1) - 15*(37 + 90*lmMgl - 90*lmMst1)
        *s2t*pow3(Mgl))) + 5*Mt*pow3(s2t)*pow4(Mst1)*(274*pow2(Mgl)*pow2(Mst1)
        - 30*Dmsqst1*(3*pow2(Mgl) + (3 - 4*shiftst1 + 2*shiftst2)*pow2(Mst1)) +
        (-41 + 120*shiftst1 - 60*shiftst2)*pow4(Mst1)) + 80*Mgl*(4*(380 + 147*
        lmMst1 - 249*lmMt)*pow2(Mgl)*pow2(Mst1) - 360*Dmsqst1*(3*(4 + lmMst1 -
        lmMt)*pow2(Mgl) + (-7 - 2*lmMst1 + 2*lmMt)*pow2(Mst1)) + (-16 + 9*
        lmMst1 - 9*lmMt)*pow4(Mst1))*pow4(Mt) + 15*shiftst3*(-24*s2t*pow3(Mt) +
        7*Mt*pow2(Mst1)*pow3(s2t))*pow6(Mst1))) + Mt*pow2(Mst2)*(-90*Mt*
        xDmsqst1*pow2(Dmsqst1)*(5*Dmst12*pow2(Mst1)*(Dmst12*(-8*Mgl*Mt*(MuSUSY*
        s2t + (17 + 6*lmMst1 - 6*lmMt)*Mt*Tbeta) + s2t*(MuSUSY*s2t - 4*Mt*
        Tbeta)*pow2(Mst1) + 6*Mgl*s2t*Tbeta*(6*Mgl*Mt + s2t*pow2(Mst1))) - 4*
        Mt*(-2*Mgl*(MuSUSY*s2t + (17 + 6*lmMst1 - 6*lmMt)*Mt*Tbeta) + 9*s2t*
        Tbeta*pow2(Mgl) + 2*s2t*Tbeta*pow2(Mst1))*pow2(Mst2)) + 80*Mgl*Tbeta*((
        22 + 6*lmMst1 - 6*lmMt)*pow2(Mgl) + (5 + 2*lmMst1 - 2*lmMt)*pow2(Mst1))
        *pow2(Mt)*pow4(Mst2)) + pow2(Mst1)*(5*Tbeta*(6*Dmst12*Mt*s2t*pow2(Mst1)
        *(-(Dmst12*Mgl*s2t*(2*(113 - 90*lmMgl + 90*lmMst1)*pow2(Mgl)*pow2(Mst1)
        + 60*Dmsqst1*(2*pow2(Mgl) + pow2(Mst1)) - 34*pow4(Mst1))) + 12*Mt*(5*
        Dmst12*(3*pow2(Mgl)*pow2(Mst1) + 2*Dmsqst1*(6*pow2(Mgl) + pow2(Mst1)) +
        pow4(Mst1)) + pow2(Mst2)*(30*Dmsqst1*pow2(Mgl) + 20*Dmsqst1*pow2(Mst1)
        - 59*pow2(Mgl)*pow2(Mst1) + 10*pow4(Mst1)))) + 90*s2t*(Dmsqst1 + pow2(
        Mst1))*pow4(Mst1)*(-8*Dmst12*shiftst2*pow2(Mst2)*pow2(Mt) + (-shiftst1
        + shiftst2)*pow2(Dmst12)*pow2(Mst1)*pow2(s2t) + 8*(shiftst1 - shiftst2)
        *pow2(Mt)*pow4(Mst2)) - 8*Mgl*pow3(Mt)*(-180*Dmsqst1*(pow2(Dmst12)*(12*
        (4 + lmMst1 - lmMt)*pow2(Mgl) + (-7 - 2*lmMst1 + 2*lmMt)*pow2(Mst1)) -
        2*Dmst12*(6*(4 + lmMst1 - lmMt)*pow2(Mgl) + (7 + 2*lmMst1 - 2*lmMt)*
        pow2(Mst1))*pow2(Mst2) - 2*(4 + lmMst1 - lmMt)*(2*pow2(Mgl) + pow2(
        Mst1))*pow4(Mst2)) + pow2(Mst1)*(pow2(Dmst12)*((3147 + 822*lmMst1 -
        822*lmMt)*pow2(Mgl) + (691 + 216*lmMst1 - 216*lmMt)*pow2(Mst1)) - 2*
        Dmst12*((1627 + 234*lmMst1 + 174*lmMt)*pow2(Mgl) + (-457 - 87*lmMst1 +
        87*lmMt)*pow2(Mst1))*pow2(Mst2) + 4*((-1070 - 207*lmMst1 + 3*lmMt)*
        pow2(Mgl) + (100 - 3*lmMst1 + 3*lmMt)*pow2(Mst1))*pow4(Mst2))) + 9*s2t*
        shiftst3*(-(pow2(Dmst12)*pow2(Mst1)*pow2(s2t)) + 8*pow2(Mt)*(pow2(
        Dmst12) - Dmst12*pow2(Mst2) + pow4(Mst2)))*pow6(Mst1)) + 2*Mt*MuSUSY*(
        4*pow2(Mgl)*pow2(Mst1)*(11*pow2(Dmst12) - 520*pow2(Mst2)*(Dmst12 +
        pow2(Mst2)))*pow2(Mt) + 5*Dmst12*s2t*(-8*Mgl*Mt*pow2(Mst2)*((113 - 90*
        lmMgl + 90*lmMst1)*pow2(Mgl)*pow2(Mst1) + 30*Dmsqst1*(2*pow2(Mgl) +
        pow2(Mst1)) - 17*pow4(Mst1)) - 30*s2t*(Dmsqst1 + pow2(Mst1))*(-(Dmst12*
        shiftst2) + 2*(shiftst1 - shiftst2)*pow2(Mst2))*pow4(Mst1) - 4*Dmst12*
        Mt*((263 - 90*lmMgl + 90*lmMst1)*pow2(Mst1)*pow3(Mgl) - 120*Dmsqst1*(-(
        Mgl*pow2(Mst1)) + pow3(Mgl)) + 30*Mgl*pow4(Mst1)) + Dmst12*s2t*pow2(
        Mst1)*(274*pow2(Mgl)*pow2(Mst1) + 30*Dmsqst1*(-3*pow2(Mgl) + (-3 +
        shiftst1)*pow2(Mst1)) + (-41 + 30*shiftst1)*pow4(Mst1)) + 3*s2t*
        shiftst3*(3*Dmst12 - 2*pow2(Mst2))*pow6(Mst1)))))) + 1440*Tbeta*
        twoLoopFlag*pow4(Mt)*((Dmst12 + pow2(Mst2))*pow3(Mgl)*pow4(Mst1)*pow4(
        Mst2) + x*pow2(Mst1)*pow2(Mst2)*(-3*pow2(Dmst12) + 3*Dmst12*pow2(Mst2)
        + 2*pow4(Mst2))*pow5(Mgl) - xDmst12*pow2(Mst1)*pow3(Dmst12)*(pow2(Mst1)
        *pow3(Mgl) - 3*x*pow5(Mgl)) + 3*x*pow6(Mst2)*pow7(Mgl))) + pow2(Mst2)*
        pow2(Mt)*(-960*Al4p*twoLoopFlag*(pow4(Mst1)*(Dmst12*s2t*(4*Mgl*Mt*
        MuSUSY*((1 - 2*lmMgl + 2*lmMst1)*pow2(Mgl) + (2 + lmMst1)*pow2(Mst1))*
        pow2(Mst2) + Dmst12*((3 + 2*lmMst1)*MuSUSY*s2t*pow2(Mgl)*pow2(Mst1) + (
        2*(3 - 2*lmMgl + 2*lmMst1)*Mt*MuSUSY + 3*(1 - 2*lmMgl + 2*lmMst1)*s2t*
        Tbeta*pow2(Mst1))*pow3(Mgl) - (1 + 2*lmMst1)*MuSUSY*s2t*pow4(Mst1) +
        Mgl*(2*Mt*MuSUSY*pow2(Mst1) + 3*(2 + lmMst1)*s2t*Tbeta*pow4(Mst1)))) +
        Mt*Tbeta*(-3*Dmst12*s2t*pow2(Mst1)*(2*Dmst12*(pow2(Mgl) + (-1 + 2*
        lmMst1)*pow2(Mst1)) + 4*((2 + lmMst1)*pow2(Mgl) - (1 + 2*lmMst1)*pow2(
        Mst1))*pow2(Mst2)) + Mt*(pow2(Dmst12)*(4*(1 + lmMt)*Mgl*pow2(Mst1) - 8*
        (-3 + lmMgl - 2*lmMst1 + lmMt)*pow3(Mgl)) + 3*Dmst12*pow2(Mst2)*(-4*(1
        + lmMt)*Mgl*pow2(Mst1) + 4*(2 + 6*lmMst1 - 3*lmMt - 2*lmMst1*lmMt +
        lmMgl*(-3 - 2*lmMst1 + 2*lmMt) + 2*pow2(lmMst1))*pow3(Mgl)) + (-24*(1 +
        lmMt)*Mgl*pow2(Mst1) - 24*(lmMgl*(1 + lmMst1 - lmMt) + lmMst1*(-2 +
        lmMt) + lmMt - pow2(lmMst1))*pow3(Mgl))*pow4(Mst2)))) - x*(Dmst12*(-1 +
        3*lmMgl - 3*lmMst1)*s2t*pow2(Mst1)*(-4*Dmst12*Mt*MuSUSY + 3*Dmst12*s2t*
        Tbeta*pow2(Mst1) + 4*Mt*MuSUSY*pow2(Mst2)) + 3*Mt*Tbeta*(4*Mt*(9 +
        lmMgl*(20 + 6*lmMst1 - 6*lmMt) + 2*lmMt + lmMst1*(-22 + 6*lmMt) - 6*
        pow2(lmMst1))*pow2(Mgl)*pow4(Mst2) + pow2(Mst1)*(4*Dmst12*(-1 + 3*lmMgl
        - 3*lmMst1)*Mgl*s2t*(Dmst12 - pow2(Mst2)) + 4*Dmst12*Mt*(lmMst1*(22 -
        6*lmMt) - 5*lmMt + lmMgl*(-17 - 6*lmMst1 + 6*lmMt) + 6*pow2(lmMst1))*(
        Dmst12 - pow2(Mst2)) + 8*Mt*(2 + lmMgl*(5 + 2*lmMst1 - 2*lmMt) + 2*
        lmMst1*(-3 + lmMt) + lmMt - 2*pow2(lmMst1))*pow4(Mst2))))*pow5(Mgl)) +
        180*Dmst12*oneLoopFlag*s2t*(Dmst12*MuSUSY*s2t + 6*Dmst12*Mt*Tbeta - 12*
        Mt*Tbeta*pow2(Mst2))*pow8(Mst1)) - 2*threeLoopFlag*pow2(Al4p)*(80*z2*
        pow2(Mst1)*(45*Mt*xDmsqst1*xDmst12*pow2(Dmsqst1)*(72*s2t*Tbeta*pow2(
        Mgl)*pow2(Mt) - 8*Mgl*Mt*(2*Mt*MuSUSY*s2t + 2*(17 + 6*lmMst1 - 6*lmMt)*
        Tbeta*pow2(Mt) - 3*Tbeta*pow2(Mst1)*pow2(s2t)) - s2t*pow2(Mst1)*(-4*Mt*
        MuSUSY*s2t + 32*Tbeta*pow2(Mt) + Tbeta*pow2(Mst1)*pow2(s2t)))*pow3(
        Dmst12) + 2*xMgl*pow4(Mgl)*(-6*Mt*xDmst12*pow3(Dmst12)*(-4*s2t*Tbeta*(
        75*Dmsqst1 + 2*(170 - 27*lmMgl + 27*lmMst1)*pow2(Mst1))*pow2(Mt) + 4*(
        49 - 17*lmMgl + 17*lmMst1)*Mt*MuSUSY*pow2(Mst1)*pow2(s2t) + 48*(17 - 4*
        lmMgl + 4*lmMst1)*MuSUSY*pow3(Mt) + (-49 + 17*lmMgl - 17*lmMst1)*Tbeta*
        pow3(s2t)*pow4(Mst1)) + pow2(Mst2)*(-12*Dmst12*s2t*Tbeta*(150*Dmsqst1*(
        Dmst12 - pow2(Mst2)) + pow2(Mst1)*(Dmst12*(146 - 27*lmMgl + 27*lmMst1)
        + 2*(194 - 27*lmMgl + 27*lmMst1)*pow2(Mst2)))*pow3(Mt) + MuSUSY*(12*(49
        - 17*lmMgl + 17*lmMst1)*pow2(Dmst12)*pow2(Mst1)*pow2(Mt)*pow2(s2t) +
        32*(9*Dmst12*(17 - 4*lmMgl + 4*lmMst1)*(Dmst12 - pow2(Mst2)) + 2*(-47 +
        12*lmMgl - 12*lmMst1)*pow4(Mst2))*pow4(Mt))))) + (3*Mt*z3*(80*Dmst12*
        pow2(Mst1)*pow2(Mt)*(64*(1440*Dmsqst1*Mt*Tbeta + (81*MuSUSY*s2t + 50*
        Mt*Tbeta)*pow2(Mst1))*pow3(Mgl) + 96*xMgl*(121*Mt*MuSUSY + 93*s2t*
        Tbeta*pow2(Mst1))*pow4(Mgl) + 3*s2t*Tbeta*pow2(Mst1)*(-245*xDmsqst1*
        pow2(Dmsqst1) - 280*Dmsqst1*pow2(Mst1) + 248*pow4(Mst1)) + 32*pow2(Mgl)
        *(-7*Mt*MuSUSY*pow2(Mst1) + 6*s2t*Tbeta*pow4(Mst1)) + 32*Mgl*(720*Mt*
        Tbeta*xDmsqst1*pow2(Dmsqst1) + 960*Dmsqst1*Mt*Tbeta*pow2(Mst1) + (-7*
        MuSUSY*s2t + 248*Mt*Tbeta)*pow4(Mst1)))*pow4(Mst2) + 40*Mt*pow2(Dmst12)
        *pow2(Mst1)*pow2(Mst2)*(-8*MuSUSY*pow2(Mgl)*pow2(Mst1)*(-7*pow2(Mt) +
        90*pow2(Mst1)*pow2(s2t)) - 24*xMgl*(-372*Mt*s2t*Tbeta*pow2(Mst1) + 968*
        MuSUSY*pow2(Mt) + 259*MuSUSY*pow2(Mst1)*pow2(s2t))*pow4(Mgl) + 3*s2t*
        pow2(Mst1)*(35*(3*MuSUSY*s2t - Mt*Tbeta)*xDmsqst1*pow2(Dmsqst1) + 70*
        Dmsqst1*(2*MuSUSY*s2t - 7*Mt*Tbeta)*pow2(Mst1) + 48*(MuSUSY*s2t - 11*
        Mt*Tbeta)*pow4(Mst1)) + 32*pow3(Mgl)*(2*Mt*(81*MuSUSY*s2t + 1096*Mt*
        Tbeta)*pow2(Mst1) - 5760*Dmsqst1*Tbeta*pow2(Mt) + 243*Tbeta*pow2(s2t)*
        pow4(Mst1)) - 16*Mgl*Tbeta*(2880*xDmsqst1*pow2(Dmsqst1)*pow2(Mt) -
        1920*Dmsqst1*pow2(Mst1)*pow2(Mt) - 1136*pow2(Mt)*pow4(Mst1) + 21*pow2(
        s2t)*pow6(Mst1))) + xDmst12*pow2(Mst1)*pow3(Dmst12)*(-1920*Mt*pow3(Mgl)
        *(16*Mt*(27*MuSUSY*s2t + 191*Mt*Tbeta)*pow2(Mst1) - 3840*Dmsqst1*Tbeta*
        pow2(Mt) + 81*Tbeta*pow2(s2t)*pow4(Mst1)) + 15*s2t*pow2(Mst1)*(140*
        xDmsqst1*pow2(Dmsqst1)*(-12*Mt*MuSUSY*s2t + 32*Tbeta*pow2(Mt) + 3*
        Tbeta*pow2(Mst1)*pow2(s2t)) + 35*Dmsqst1*pow2(Mst1)*(8*Mt*MuSUSY*s2t +
        31*Tbeta*pow2(Mt) + 16*Tbeta*pow2(Mst1)*pow2(s2t)) + 8*(22*Mt*MuSUSY*
        s2t + 85*Tbeta*pow2(Mt) + 24*Tbeta*pow2(Mst1)*pow2(s2t))*pow4(Mst1)) +
        480*xMgl*pow4(Mgl)*(-2976*s2t*Tbeta*pow2(Mst1)*pow2(Mt) + 1036*Mt*
        MuSUSY*pow2(Mst1)*pow2(s2t) + 1936*MuSUSY*pow3(Mt) - 259*Tbeta*pow3(
        s2t)*pow4(Mst1)) + 128*Mgl*Mt*(14400*Tbeta*xDmsqst1*pow2(Dmsqst1)*pow2(
        Mt) - 38400*Dmsqst1*Tbeta*pow2(Mst1)*pow2(Mt) + (6*Mt*MuSUSY*s2t - 440*
        Tbeta*pow2(Mt))*pow4(Mst1) + 105*Tbeta*pow2(s2t)*pow6(Mst1)) - 320*
        pow2(Mgl)*(-42*MuSUSY*pow2(Mst1)*pow3(Mt) - Mt*s2t*(45*MuSUSY*s2t + 7*
        Mt*Tbeta)*pow4(Mst1) + 45*Tbeta*pow3(s2t)*pow6(Mst1))) + 2560*Mgl*pow3(
        Mt)*(242*MuSUSY*xMgl*pow2(Mst1)*pow3(Mgl) + 20*Tbeta*pow2(Mgl)*(72*
        xDmsqst1*pow2(Dmsqst1) + 48*Dmsqst1*pow2(Mst1) - 19*pow4(Mst1)) - 7*
        Mgl*MuSUSY*pow4(Mst1) + 16*Tbeta*pow2(Mst1)*(30*xDmsqst1*pow2(Dmsqst1)
        + 30*Dmsqst1*pow2(Mst1) + pow4(Mst1)))*pow6(Mst2)))/4. + 20*Mt*pow2(
        Mst2)*(15*xDmsqst1*pow2(Dmsqst1)*(Dmst12*Mt*MuSUSY*s2t*pow2(Mst1)*(-
        256*Dmst12*Mgl*Mt + Dmst12*s2t*(21 - 8*(shiftst1 + shiftst2))*pow2(
        Mst1) + 256*Mgl*Mt*pow2(Mst2) + 16*s2t*(shiftst1 - shiftst2)*pow2(Mst1)
        *pow2(Mst2)) - 3*Dmst12*Tbeta*pow2(Mst1)*(2*(16*(1 - 19*lmMst1 + 19*
        lmMt)*Mgl*Mt + 192*s2t*pow2(Mgl) + s2t*(47 - 16*shiftst2)*pow2(Mst1))*
        pow2(Mst2)*pow2(Mt) - Dmst12*(s2t*(384*pow2(Mgl) - 45*pow2(Mst1))*pow2(
        Mt) + 32*Mgl*Mt*((1 - 19*lmMst1 + 19*lmMt)*pow2(Mt) + 2*pow2(Mst1)*
        pow2(s2t)) + 4*(shiftst1 - shiftst2)*pow3(s2t)*pow4(Mst1))) + 96*Tbeta*
        pow2(Mt)*(4*(-1 + 3*lmMst1 - 3*lmMt)*Mgl*Mt*pow2(Mst1) + 4*Mt*(4 +
        lmMgl*(1 + lmMst1 - lmMt) - 8*lmMt + lmMst1*(7 + lmMt) - pow2(lmMst1))*
        pow3(Mgl) + s2t*(-shiftst1 + shiftst2)*pow4(Mst1))*pow4(Mst2)) + 2*
        pow2(Mst1)*(4*Mt*MuSUSY*((pow2(Dmst12)*(-2*pow2(Mgl)*pow2(Mst1)*(4*(163
        + 15*lmMgl - 27*lmMst1 + 20*lmMt - 8*pow2(lmMst1))*pow2(Mt) - 5*(360*
        Dmsqst1 + (200 + 18*lmMgl + 265*lmMst1 + 91*pow2(lmMst1))*pow2(Mst1))*
        pow2(s2t)) + Mt*s2t*(40*Mgl*pow2(Mst1)*(210*Dmsqst1 + (149 + 41*lmMst1)
        *pow2(Mst1)) - 80*(60*Dmsqst1*(3 - lmMgl + lmMst1) + (-249 - 176*lmMst1
        - 4*lmMgl*(-41 + 5*lmMst1) + 27*pow2(lmMgl) + pow2(lmMst1))*pow2(Mst1))
        *pow3(Mgl)) + 5*(15*Dmsqst1*(3 - 6*shiftst1 - 6*shiftst2 + 4*lmMst1*(-2
        + shiftst1 + shiftst2)) - (32 + 30*(shiftst1 + shiftst2) + 2*lmMst1*(
        115 - 30*(shiftst1 + shiftst2) - 9*shiftst3) + 21*shiftst3 + 82*pow2(
        lmMst1))*pow2(Mst1))*pow2(s2t)*pow4(Mst1)))/10. - Dmst12*pow2(Mst2)*(
        16*(16 - 9*lmMgl + 41*lmMst1 - 12*lmMt + 4*pow2(lmMst1))*pow2(Mgl)*
        pow2(Mst1)*pow2(Mt) - 3*(10*Dmsqst1*(3 - 2*lmMst1)*(shiftst1 -
        shiftst2) + (1 - 2*lmMst1)*(10*shiftst1 - 10*shiftst2 + shiftst3)*pow2(
        Mst1))*pow2(s2t)*pow4(Mst1) + Mt*s2t*(16*(-11 + 154*lmMgl - 148*lmMst1
        - 20*lmMgl*lmMst1 + 27*pow2(lmMgl) + pow2(lmMst1))*pow2(Mst1)*pow3(Mgl)
        + 240*Dmsqst1*(-3*Mgl*pow2(Mst1) + 2*(-3 + lmMgl - lmMst1)*pow3(Mgl)) -
        4*Mgl*(215 + 118*lmMst1 + 41*pow2(lmMst1))*pow4(Mst1))) + 8*(5 + 18*
        lmMgl - 74*lmMst1 + 24*lmMt - 8*pow2(lmMst1))*pow2(Mgl)*pow2(Mst1)*
        pow2(Mt)*pow4(Mst2)) - Tbeta*(3*Dmst12*Mt*s2t*pow2(Mst1)*(-4*Dmst12*
        Mgl*s2t*(-4*(-11 + lmMgl*(154 - 20*lmMst1) - 148*lmMst1 + 27*pow2(
        lmMgl) + pow2(lmMst1))*pow2(Mgl)*pow2(Mst1) + 60*Dmsqst1*((6 - 2*lmMgl
        + 2*lmMst1)*pow2(Mgl) + 3*pow2(Mst1)) + (215 + 118*lmMst1 + 41*pow2(
        lmMst1))*pow4(Mst1)) + Mt*(4*pow2(Mst2)*(2*(501 - 18*lmMgl + 504*lmMst1
        - 48*lmMt + 107*pow2(lmMst1))*pow2(Mgl)*pow2(Mst1) + 30*Dmsqst1*(36*
        pow2(Mgl) + (5 - 4*lmMst1)*pow2(Mst1)) - 4*(-6 + 79*lmMst1 + 41*pow2(
        lmMst1))*pow4(Mst1)) + Dmst12*(8*(522 + 107*lmMst1)*pow2(Mgl)*pow2(
        Mst1) + 15*Dmsqst1*(336*pow2(Mgl) + 31*pow2(Mst1)) + 4*(-95 - 66*lmMst1
        + 82*pow2(lmMst1))*pow4(Mst1)))) + 90*s2t*(Dmsqst1*(-3 + 2*lmMst1) + (-
        1 + 2*lmMst1)*pow2(Mst1))*pow4(Mst1)*(8*Dmst12*shiftst2*pow2(Mst2)*
        pow2(Mt) + (shiftst1 - shiftst2)*pow2(Dmst12)*pow2(Mst1)*pow2(s2t) + 8*
        (-shiftst1 + shiftst2)*pow2(Mt)*pow4(Mst2)) + 16*Mgl*pow3(Mt)*(30*
        Dmsqst1*(pow2(Dmst12)*(6*(31 + 27*lmMst1 - 24*lmMt - 2*lmMst1*lmMt +
        lmMgl*(-3 - 2*lmMst1 + 2*lmMt) + 2*pow2(lmMst1))*pow2(Mgl) + (-13 - 19*
        lmMst1 + 19*lmMt)*pow2(Mst1)) - 6*Dmst12*((31 + 27*lmMst1 - 24*lmMt -
        2*lmMst1*lmMt + lmMgl*(-3 - 2*lmMst1 + 2*lmMt) + 2*pow2(lmMst1))*pow2(
        Mgl) + (3 + 7*lmMst1 - 7*lmMt)*pow2(Mst1))*pow2(Mst2) - 12*((7 + 6*
        lmMst1 - 5*lmMt - lmMst1*lmMt + lmMgl*(-1 - lmMst1 + lmMt) + pow2(
        lmMst1))*pow2(Mgl) + (1 + 2*lmMst1 - 2*lmMt)*pow2(Mst1))*pow4(Mst2)) -
        pow2(Mst1)*(pow2(Dmst12)*(-((-3749 - 2374*lmMst1 + 1770*lmMt - 117*
        lmMst1*lmMt + lmMgl*(96 - 237*lmMst1 + 157*lmMt) + 108*pow2(lmMgl) +
        41*pow2(lmMst1) + 48*pow2(lmMt))*pow2(Mgl)) + (lmMst1*(297 + 2*lmMt) +
        15*pow2(lmMst1) + 3*(47 - 79*lmMt + 8*pow2(lmMt)))*pow2(Mst1)) +
        Dmst12*pow2(Mst2)*(-((-605 + 6*lmMst1*(-77 + lmMt) + 714*lmMt + 45*
        pow2(lmMst1) + 72*pow2(lmMt))*pow2(Mst1)) + 2*pow2(Mgl)*(2096 - 1167*
        lmMt - 81*(4 + 2*lmMst1 - 3*lmMt)*pow2(lmMgl) + 3*(129 + 61*lmMt)*pow2(
        lmMst1) - 108*pow2(lmMt) - 3*lmMst1*(-793 + 67*lmMt + 24*pow2(lmMt)) +
        3*lmMgl*(-229 + lmMst1*(63 - 142*lmMt) + 35*lmMt + 91*pow2(lmMst1) +
        24*pow2(lmMt)) - 27*pow3(lmMgl) - 84*pow3(lmMst1))) + 2*(-((-233 + 6*
        lmMst1*(-32 + lmMt) + 348*lmMt + 45*pow2(lmMst1) + 72*pow2(lmMt))*pow2(
        Mst1)) + pow2(Mgl)*(758 - 570*lmMt - 81*(3 + 2*lmMst1 - 3*lmMt)*pow2(
        lmMgl) + 3*(141 + 61*lmMt)*pow2(lmMst1) - 72*pow2(lmMt) - 6*lmMst1*(-
        236 + 49*lmMt + 12*pow2(lmMt)) + 3*lmMgl*(-200 + 2*lmMst1 + 76*lmMt -
        142*lmMst1*lmMt + 91*pow2(lmMst1) + 24*pow2(lmMt)) - 27*pow3(lmMgl) -
        84*pow3(lmMst1)))*pow4(Mst2))) + shiftst3*(72*s2t*pow2(Mt)*(-2*(-2 +
        lmMst1)*pow2(Dmst12) + Dmst12*(-3 + 2*lmMst1)*pow2(Mst2) + (1 - 2*
        lmMst1)*pow4(Mst2))*pow6(Mst1) + 9*(-1 + 2*lmMst1)*pow2(Dmst12)*pow3(
        s2t)*pow8(Mst1))))))))/(8640.*Tbeta*pow6(Mst2)*pow8(Mst1));
}

/**
 * @return returns the susy log^0 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double H9::calc_coef_at_as2_no_sm_logs_log0() const {

   const double result =
      -(432*pow2(log(pow2(Mst1)/pow2(Mgl)))*pow4(Mst1)*(-2*Mt*pow2(Dmst12)*
        pow2(Mst2)*(-72*Mt*s2t*pow2(Mst1) + 202*Mgl*pow2(Mt) + 111*Mgl*pow2(
        Mst1)*pow2(s2t))*pow3(Mgl) + pow3(Dmst12)*pow3(Mgl)*(-48*s2t*pow2(Mst1)
        *pow2(Mt) + 111*Mgl*Mt*pow2(Mst1)*pow2(s2t) + 1568*Mgl*pow3(Mt) + 24*
        pow3(s2t)*pow4(Mst1)) - 8*Dmst12*pow2(Mgl)*(95*Mt*pow2(Mgl) - 18*Mt*
        pow2(Mst1) - 81*Mgl*s2t*pow2(Mst1))*pow2(Mt)*pow4(Mst2) - 8*pow3(Mt)*(-
        36*pow2(Mgl)*pow2(Mst1) + 77*pow4(Mgl) + 4*pow4(Mst1))*pow6(Mst2)) +
        192*pow2(Mt)*pow3(log(pow2(Mst1)/pow2(Mgl)))*pow4(Mst1)*(-95*Mt*pow2(
        Dmst12)*pow2(Mst2)*pow4(Mgl) + 504*Mt*pow3(Dmst12)*pow4(Mgl) + 2*
        Dmst12*(-157*Mgl*Mt + 54*s2t*pow2(Mst1))*pow3(Mgl)*pow4(Mst2) + 2*Mt*(-
        157*pow4(Mgl) + 18*pow4(Mst1))*pow6(Mst2)) + 15*pow2(Dmsqst1)*(Mt*pow2(
        Dmst12)*pow2(Mst1)*pow2(Mst2)*(576*Mgl*Mt*s2t*(-13 + 32*z2 - 32*z3)*
        pow2(Mst1) + (-418 - 2880*z2 + 4419*z3)*pow2(Mst1)*pow2(Mt) + 1728*
        pow2(Mgl)*((1 + 34*z2 - 48*z3)*pow2(Mt) + (4 - 3*z2)*pow2(Mst1)*pow2(
        s2t)) - 4608*Mt*s2t*(-8 + 25*z2 - 24*z3)*pow3(Mgl) + 18*(94 - 64*z2 +
        49*z3)*pow2(s2t)*pow4(Mst1)) + pow2(Mst1)*pow3(Dmst12)*(4608*s2t*(-8 +
        25*z2 - 24*z3)*pow2(Mt)*pow3(Mgl) + 2*(2938 + 4608*z2 - 8271*z3)*pow2(
        Mst1)*pow3(Mt) - 1728*pow2(Mgl)*(2*Mt*(4 - 3*z2)*pow2(Mst1)*pow2(s2t) +
        (1 + 34*z2 - 48*z3)*pow3(Mt)) + 9*Mt*(-98 + 64*z2 - 91*z3)*pow2(s2t)*
        pow4(Mst1) + 192*Mgl*(6*s2t*(15 - 56*z2 + 64*z3)*pow2(Mst1)*pow2(Mt) +
        (-4 + 3*z2)*pow3(s2t)*pow4(Mst1))) - 72*Dmst12*pow2(Mst1)*pow2(Mt)*(24*
        Mt*(1 + 34*z2 - 48*z3)*pow2(Mgl) + Mt*(70 + 48*z2 - 107*z3)*pow2(Mst1)
        + 32*Mgl*s2t*(1 - 12*z2 + 16*z3)*pow2(Mst1) - 64*s2t*(-8 + 25*z2 - 24*
        z3)*pow3(Mgl))*pow4(Mst2) - 16*pow3(Mt)*(432*(1 + 5*z2 - 8*z3)*pow2(
        Mgl)*pow2(Mst1) + 36*(-41 + 220*z2 - 240*z3)*pow4(Mgl) + (382 + 144*z2
        - 387*z3)*pow4(Mst1))*pow6(Mst2)) + 10*Dmsqst1*pow2(Mst1)*(-216*Dmst12*
        pow2(Mt)*(48*Mt*(-3 + 14*z2 - 16*z3)*pow2(Mgl)*pow2(Mst1) + 32*s2t*(19
        - 18*z2 + 16*z3)*pow2(Mst1)*pow3(Mgl) + 4*Mt*(-307 + 480*z2 - 480*z3)*
        pow4(Mgl) + Mt*(38 + 48*z2 - 107*z3)*pow4(Mst1) - 32*Mgl*s2t*(-4 + 9*z2
        - 8*z3)*pow4(Mst1))*pow4(Mst2) + 108*Mt*pow2(Dmst12)*pow2(Mst2)*(256*
        Mt*s2t*(-6 + 9*z2 - 8*z3)*pow2(Mst1)*pow3(Mgl) + 8*((-307 + 480*z2 -
        480*z3)*pow2(Mt) + (29 - 10*z2)*pow2(Mst1)*pow2(s2t))*pow4(Mgl) + 32*
        Mgl*Mt*s2t*(-5 + 14*z2 - 16*z3)*pow4(Mst1) + (-6 - 80*z2 + 107*z3)*
        pow2(Mt)*pow4(Mst1) - 24*pow2(Mgl)*((-19 + 70*z2 - 80*z3)*pow2(Mst1)*
        pow2(Mt) + 2*(-3 + z2)*pow2(s2t)*pow4(Mst1)) + 2*(10 - 16*z2 + 7*z3)*
        pow2(s2t)*pow6(Mst1)) + pow3(Dmst12)*(-864*(2*Mt*(29 - 10*z2)*pow2(
        Mst1)*pow2(s2t) + (-307 + 480*z2 - 480*z3)*pow3(Mt))*pow4(Mgl) + (1814
        - 1152*z2 + 567*z3)*pow3(Mt)*pow4(Mst1) - 2592*pow2(Mgl)*(2*(25 - 98*z2
        + 112*z3)*pow2(Mst1)*pow3(Mt) + Mt*(-1 + 2*z2)*pow2(s2t)*pow4(Mst1)) +
        1152*pow3(Mgl)*(-6*s2t*(-67 + 90*z2 - 80*z3)*pow2(Mst1)*pow2(Mt) + (-3
        + z2)*pow3(s2t)*pow4(Mst1)) + 27*Mt*(-18 + 64*z2 - 7*z3)*pow2(s2t)*
        pow6(Mst1) + 576*Mgl*(3*s2t*(-5 + 4*z2)*pow2(Mt)*pow4(Mst1) + (-3 + z2)
        *pow3(s2t)*pow6(Mst1))) - 192*pow3(Mt)*(108*(-1 + 4*z2 - 4*z3)*pow2(
        Mgl)*pow2(Mst1) + 9*(-69 + 80*z2 - 80*z3)*pow4(Mgl) + 4*(23 + 9*z2 -
        36*z3)*pow4(Mst1))*pow6(Mst2)) - 8*pow4(Mst1)*(-24*Dmst12*pow2(Mt)*(Mt*
        (3049 - 2716*z2 + 2388*z3)*pow2(Mgl)*pow2(Mst1) - 4*s2t*(1043 + 2143*z2
        - 570*z3)*pow2(Mst1)*pow3(Mgl) + 2*Mt*(6742 + 4265*z2 - 174*z3)*pow4(
        Mgl) - 4*Mt*(49 + 90*z2 - 171*z3)*pow4(Mst1) + 4*Mgl*s2t*(-407 + 197*z2
        - 24*z3)*pow4(Mst1))*pow4(Mst2) + 3*Mt*pow2(Dmst12)*pow2(Mst2)*(-48*Mt*
        s2t*(-1091 + 402*z2 - 480*z3)*pow2(Mst1)*pow3(Mgl) + 4*((-40279 +
        17982*z2 - 26652*z3)*pow2(Mt) - 3*(681 + 776*z2 - 558*z3)*pow2(Mst1)*
        pow2(s2t))*pow4(Mgl) - 16*Mgl*Mt*s2t*(-148 + 607*z2 - 696*z3)*pow4(
        Mst1) + (262 + 1840*z2 - 3423*z3)*pow2(Mt)*pow4(Mst1) - 4*pow2(Mgl)*(2*
        (1327 - 4904*z2 + 5280*z3)*pow2(Mst1)*pow2(Mt) + 9*(175 + 118*z2 - 4*
        z3)*pow2(s2t)*pow4(Mst1)) + 18*(-8 + 40*z2 + 31*z3)*pow2(s2t)*pow6(
        Mst1)) + pow3(Dmst12)*(-6*(3*Mt*(1695 - 968*z2 + 558*z3)*pow2(Mst1)*
        pow2(s2t) + 4*(-53763 + 9452*z2 - 26304*z3)*pow3(Mt))*pow4(Mgl) + (-
        4486 + 816*z2 + 3375*z3)*pow3(Mt)*pow4(Mst1) + 12*pow2(Mgl)*(4*(-570 +
        929*z2 - 840*z3)*pow2(Mst1)*pow3(Mt) + 9*Mt*(1 + 148*z2 - 4*z3)*pow2(
        s2t)*pow4(Mst1)) - 24*pow3(Mgl)*(s2t*(-2722 + 4704*z2 - 3696*z3)*pow2(
        Mst1)*pow2(Mt) + (-22 + 113*z2 - 243*z3)*pow3(s2t)*pow4(Mst1)) - 18*Mt*
        (-119 + 60*z2 + 192*z3)*pow2(s2t)*pow6(Mst1) + 12*Mgl*(s2t*(446 - 768*
        z2 + 624*z3)*pow2(Mt)*pow4(Mst1) + (215 + 34*z2 - 21*z3)*pow3(s2t)*
        pow6(Mst1))) - 12*pow3(Mt)*(4*(1909 - 496*z2 + 228*z3)*pow2(Mgl)*pow2(
        Mst1) + 2*(6677 + 10414*z2 - 3948*z3)*pow4(Mgl) + (-4125 - 892*z2 +
        4740*z3)*pow4(Mst1))*pow6(Mst2)) - 96*log(pow2(Mst1)/pow2(Mgl))*(-180*
        pow2(Dmsqst1)*pow2(Mt)*pow3(Mgl)*(-6*Dmst12*s2t*pow2(Mst1)*(pow2(
        Dmst12) - Dmst12*pow2(Mst2) + pow4(Mst2)) + 5*Mgl*Mt*pow6(Mst2)) - 60*
        Dmsqst1*pow2(Mst1)*pow3(Mgl)*(3*Mt*pow2(Dmst12)*pow2(Mst2)*(4*Mt*s2t*
        pow2(Mst1) + 15*Mgl*pow2(Mt) - 5*Mgl*pow2(Mst1)*pow2(s2t)) + pow3(
        Dmst12)*(-60*s2t*pow2(Mst1)*pow2(Mt) + 30*Mgl*Mt*pow2(Mst1)*pow2(s2t) -
        45*Mgl*pow3(Mt) + 2*pow3(s2t)*pow4(Mst1)) - 9*Dmst12*(5*Mgl*Mt - 4*s2t*
        pow2(Mst1))*pow2(Mt)*pow4(Mst2) - 30*Mgl*pow3(Mt)*pow6(Mst2)) + pow4(
        Mst1)*(2*Mt*pow2(Dmst12)*pow2(Mgl)*pow2(Mst2)*(-102*Mgl*Mt*s2t*pow2(
        Mst1) + 2*pow2(Mgl)*((737 + 162*z2)*pow2(Mt) + 3*(328 + 27*z2)*pow2(
        Mst1)*pow2(s2t)) + 27*pow2(s2t)*pow4(Mst1)) - pow2(Mgl)*pow3(Dmst12)*(
        Mt*pow2(Mgl)*(32*(434 + 81*z2)*pow2(Mt) + 3*(605 + 54*z2)*pow2(Mst1)*
        pow2(s2t)) + 54*Mt*pow2(s2t)*pow4(Mst1) - 2*Mgl*(67*s2t*pow2(Mst1)*
        pow2(Mt) + 2*(-154 + 45*z2)*pow3(s2t)*pow4(Mst1))) + 12*Dmst12*pow2(
        Mgl)*(18*Mt*(37 + 6*z2)*pow2(Mgl) + 9*Mt*pow2(Mst1) - 476*Mgl*s2t*pow2(
        Mst1))*pow2(Mt)*pow4(Mst2) + 12*pow3(Mt)*(18*pow2(Mgl)*pow2(Mst1) + 3*(
        215 + 36*z2)*pow4(Mgl) - 52*pow4(Mst1))*pow6(Mst2))))/(108.*pow3(Mt)*
        pow6(Mst2)*pow8(Mst1));

   return result;
}

/**
 * @return returns the susy log^1 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double H9::calc_coef_at_as2_no_sm_logs_log1() const {

   const double result =
      (-8*(6*pow2(Mt)*pow2(log(pow2(Mst1)/pow2(Mgl)))*pow4(Mst1)*(-203*Mt*pow2(
        Dmst12)*pow2(Mst2)*pow4(Mgl) + 936*Mt*pow3(Dmst12)*pow4(Mgl) + 2*
        Dmst12*(-265*Mgl*Mt + 162*s2t*pow2(Mst1))*pow3(Mgl)*pow4(Mst2) + 2*Mt*(
        -265*pow4(Mgl) + 18*pow4(Mst1))*pow6(Mst2)) + 15*pow2(Dmsqst1)*pow2(Mt)
        *(pow2(Dmst12)*pow2(Mst1)*pow2(Mst2)*(18*Mt*(-19 + 12*z2)*pow2(Mgl) +
        Mt*(19 - 12*z2)*pow2(Mst1) + 12*Mgl*s2t*(-7 + 4*z2)*pow2(Mst1) - 96*
        s2t*(-4 + 3*z2)*pow3(Mgl)) + 2*pow2(Mst1)*pow3(Dmst12)*(-9*Mt*(-19 +
        12*z2)*pow2(Mgl) + 12*Mgl*s2t*(13 - 8*z2)*pow2(Mst1) + 8*Mt*(-5 + 3*z2)
        *pow2(Mst1) + 48*s2t*(-4 + 3*z2)*pow3(Mgl)) + 6*Dmst12*pow2(Mst1)*(Mt*(
        57 - 36*z2)*pow2(Mgl) + Mt*(7 - 4*z2)*pow2(Mst1) + 8*Mgl*s2t*(-3 + 2*
        z2)*pow2(Mst1) + 16*s2t*(-4 + 3*z2)*pow3(Mgl))*pow4(Mst2) + 2*Mt*(-36*(
        -3 + 2*z2)*pow2(Mgl)*pow2(Mst1) - 9*(-27 + 20*z2)*pow4(Mgl) + 4*(5 - 3*
        z2)*pow4(Mst1))*pow6(Mst2)) - 10*Dmsqst1*Mt*pow2(Mst1)*(18*Dmst12*Mt*(
        6*Mt*(-7 + 4*z2)*pow2(Mgl)*pow2(Mst1) - 8*s2t*(-5 + 2*z2)*pow2(Mst1)*
        pow3(Mgl) + Mt*(-119 + 60*z2)*pow4(Mgl) + 4*Mt*(-2 + z2)*pow4(Mst1) -
        8*Mgl*s2t*(-2 + z2)*pow4(Mst1))*pow4(Mst2) + 2*pow3(Dmst12)*(-27*(-47 +
        28*z2)*pow2(Mgl)*pow2(Mst1)*pow2(Mt) + 72*Mt*s2t*(-19 + 10*z2)*pow2(
        Mst1)*pow3(Mgl) + 9*(-119 + 60*z2)*pow2(Mt)*pow4(Mgl) + 6*Mgl*Mt*s2t*
        pow4(Mst1) + pow2(Mt)*pow4(Mst1) - 9*pow2(s2t)*pow6(Mst1)) - 9*pow2(
        Dmst12)*pow2(Mst2)*(3*(33 - 20*z2)*pow2(Mgl)*pow2(Mst1)*pow2(Mt) + 16*
        Mt*s2t*(-7 + 4*z2)*pow2(Mst1)*pow3(Mgl) + 2*(-119 + 60*z2)*pow2(Mt)*
        pow4(Mgl) + 8*Mgl*Mt*s2t*(-3 + 2*z2)*pow4(Mst1) + (7 - 4*z2)*pow2(Mt)*
        pow4(Mst1) - 2*pow2(s2t)*pow6(Mst1)) + 6*pow2(Mt)*(36*(-2 + z2)*pow2(
        Mgl)*pow2(Mst1) + 3*(-49 + 20*z2)*pow4(Mgl) + 2*(-17 + 6*z2)*pow4(Mst1)
        )*pow6(Mst2)) + pow4(Mst1)*(Mt*pow2(Dmst12)*pow2(Mst2)*(6*pow2(Mgl)*
        pow2(Mst1)*((659 - 440*z2)*pow2(Mt) + 243*pow2(Mst1)*pow2(s2t)) + 72*
        Mt*s2t*(-99 + 20*z2)*pow2(Mst1)*pow3(Mgl) + ((23651 - 5436*z2)*pow2(Mt)
        + 582*pow2(Mst1)*pow2(s2t))*pow4(Mgl) + 12*Mgl*Mt*s2t*(-25 + 62*z2)*
        pow4(Mst1) - 6*((-145 + 46*z2)*pow2(Mt) + 79*pow2(Mst1)*pow2(s2t))*
        pow4(Mst1)) + 6*Dmst12*pow2(Mt)*(4*Mt*(103 - 50*z2)*pow2(Mgl)*pow2(
        Mst1) - 4*s2t*(283 + 138*z2)*pow2(Mst1)*pow3(Mgl) + Mt*(2647 + 628*z2)*
        pow4(Mgl) + 12*Mt*(8 - 5*z2)*pow4(Mst1) - 4*Mgl*s2t*(65 + 2*z2)*pow4(
        Mst1))*pow4(Mst2) + pow3(Dmst12)*(-3*Mt*pow2(Mgl)*pow2(Mst1)*((-498 +
        280*z2)*pow2(Mt) + 379*pow2(Mst1)*pow2(s2t)) + 24*s2t*pow2(Mst1)*((-86
        + 77*z2)*pow2(Mt) + pow2(Mst1)*pow2(s2t))*pow3(Mgl) + (315*Mt*pow2(
        Mst1)*pow2(s2t) + 16*(-3949 + 444*z2)*pow3(Mt))*pow4(Mgl) - 2*Mgl*s2t*(
        (442 - 60*z2)*pow2(Mt) + 59*pow2(Mst1)*pow2(s2t))*pow4(Mst1) + 12*(-43
        + 8*z2)*pow3(Mt)*pow4(Mst1) + 375*Mt*pow2(s2t)*pow6(Mst1)) + 12*pow3(
        Mt)*(2*(91 - 10*z2)*pow2(Mgl)*pow2(Mst1) + (589 + 464*z2)*pow4(Mgl) -
        4*(-35 + 15*z2 + 3*z3)*pow4(Mst1))*pow6(Mst2)) - 2*log(pow2(Mst1)/pow2(
        Mgl))*(180*Dmsqst1*pow2(Mst1)*pow2(Mt)*pow3(Mgl)*(-5*Mgl*Mt*pow2(
        Dmst12)*pow2(Mst2) + (5*Mgl*Mt + 4*s2t*pow2(Mst1))*pow3(Dmst12) +
        Dmst12*(5*Mgl*Mt - 4*s2t*pow2(Mst1))*pow4(Mst2) + 5*Mgl*Mt*pow6(Mst2))
        - 90*pow2(Dmsqst1)*pow2(Mt)*pow3(Mgl)*(-4*Dmst12*s2t*pow2(Mst1)*(pow2(
        Dmst12) - Dmst12*pow2(Mst2) + pow4(Mst2)) + 5*Mgl*Mt*pow6(Mst2)) +
        pow4(Mst1)*(6*Mt*pow2(Dmst12)*pow2(Mst2)*(14*Mt*s2t*pow2(Mst1) + 67*
        Mgl*pow2(Mt) + 117*Mgl*pow2(Mst1)*pow2(s2t))*pow3(Mgl) - pow3(Dmst12)*
        pow3(Mgl)*(42*s2t*pow2(Mst1)*pow2(Mt) + 351*Mgl*(Mt*pow2(Mst1)*pow2(
        s2t) + 8*pow3(Mt)) + 68*pow3(s2t)*pow4(Mst1)) + 12*Dmst12*pow2(Mgl)*(
        167*Mt*pow2(Mgl) - 54*Mt*pow2(Mst1) - 170*Mgl*s2t*pow2(Mst1))*pow2(Mt)*
        pow4(Mst2) + 36*pow3(Mt)*(-36*pow2(Mgl)*pow2(Mst1) + 59*pow4(Mgl) + 4*
        pow4(Mst1))*pow6(Mst2)))))/(9.*pow3(Mt)*pow6(Mst2)*pow8(Mst1));

   return result;
}

/**
 * @return returns the susy log^2 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double H9::calc_coef_at_as2_no_sm_logs_log2() const {

   const double result =
      (-8*(6*Dmst12*pow2(Mt)*(43*Mt*pow2(Mgl)*pow2(Mst1) - 248*s2t*pow2(Mst1)*
        pow3(Mgl) + 550*Mt*pow4(Mgl) + 51*Mt*pow4(Mst1) + 60*Mgl*s2t*pow4(Mst1)
        )*pow4(Mst2) + Mt*pow2(Dmst12)*pow2(Mst2)*(-528*Mt*s2t*pow2(Mst1)*pow3(
        Mgl) + 2*(983*pow2(Mt) - 48*pow2(Mst1)*pow2(s2t))*pow4(Mgl) - 180*Mgl*
        Mt*s2t*pow4(Mst1) + pow2(Mgl)*(-64*pow2(Mst1)*pow2(Mt) + 321*pow2(s2t)*
        pow4(Mst1)) - 3*(51*pow2(Mt)*pow4(Mst1) + 82*pow2(s2t)*pow6(Mst1))) +
        pow3(Dmst12)*((48*Mt*pow2(Mst1)*pow2(s2t) - 7232*pow3(Mt))*pow4(Mgl) +
        102*pow3(Mt)*pow4(Mst1) + pow2(Mgl)*(32*pow2(Mst1)*pow3(Mt) - 321*Mt*
        pow2(s2t)*pow4(Mst1)) + 16*pow3(Mgl)*(11*s2t*pow2(Mst1)*pow2(Mt) + 2*
        pow3(s2t)*pow4(Mst1)) + 369*Mt*pow2(s2t)*pow6(Mst1) + Mgl*(120*s2t*
        pow2(Mt)*pow4(Mst1) - 41*pow3(s2t)*pow6(Mst1))) + 6*pow3(Mt)*(15*pow2(
        Dmsqst1) - 30*Dmsqst1*pow2(Mst1) + 86*pow2(Mgl)*pow2(Mst1) + 380*pow4(
        Mgl) - 98*pow4(Mst1))*pow6(Mst2) + 12*log(pow2(Mst1)/pow2(Mgl))*pow2(
        Mt)*(-53*Mt*pow2(Dmst12)*pow2(Mst2)*pow4(Mgl) + 276*Mt*pow3(Dmst12)*
        pow4(Mgl) + 2*Dmst12*(-85*Mgl*Mt + 44*s2t*pow2(Mst1))*pow3(Mgl)*pow4(
        Mst2) + 2*Mt*(-85*pow4(Mgl) + 9*pow4(Mst1))*pow6(Mst2))))/(9.*pow3(Mt)*
        pow4(Mst1)*pow6(Mst2));

   return result;
}

/**
 * @return returns the susy log^3 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double H9::calc_coef_at_as2_no_sm_logs_log3() const {

   const double result =
      -298.6666666666667;

   return result;
}

}        // hierarchies
}        // himalaya
