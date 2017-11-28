#define Pi M_PI

#include "H9q2.hpp"
#include "HierarchyCalculator.hpp"
#include "Utils.hpp"
#include <type_traits>
#include <cmath>

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
himalaya::H9q2::H9q2(std::map<unsigned int, unsigned int> flagMap, double Al4p, double beta, double Dmst12, double Dmsqst1,
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
   // zeta functions
   z2 = pow2(Pi) / 6.;
   z3 = 1.202056903159594;
   // mdr flags, indicates if one wants to shift the dr stop mass to the mdr stop mass
   shiftst1 = mdrFlag;
   shiftst2 = mdrFlag;
   shiftst3 = mdrFlag;
   // expansion flags
   x = flagMap.at(HierarchyCalculator::xx);
   xDmst12 = flagMap.at(HierarchyCalculator::xxDmglst1);
   xDmsqst1 = flagMap.at(HierarchyCalculator::xxDmsqst1);
   xMgl = flagMap.at(HierarchyCalculator::xxMgl);
   
   s1 = 
   #include "../hierarchies/h9q2/sigS1Full.inc"
   ;
   s2 = 
   #include "../hierarchies/h9q2/sigS2Full.inc"
   ;
   s12 = 
   #include "../hierarchies/h9q2/sigS12Full.inc"
   ;
}

/**
 * 	@return The diagonal (1, 1) matrix element of the Higgs mass matrix as a double for the hierarchy 'H9q2'
 */
double himalaya::H9q2::getS1(){
   return s1;
}

/**
 * 	@return The diagonal (2, 2) matrix element of the Higgs mass matrix as a double for the hierarchy 'H9q2'
 */
double himalaya::H9q2::getS2(){
   return s2;
}

/**
 * 	@return The off-diagonal (1, 2) = (2, 1) matrix element of the Higgs mass matrix as a double for the hierarchy 'H9q2'
 */
double himalaya::H9q2::getS12(){
   return s12;
}

/**
 * 	@return returns the constant term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double himalaya::H9q2::calc_at_as2_no_logs(){

   const double result =
      ((Mt*pow2(Sbeta)*(48*pow3(Mt)*(2750 + (3680*Dmsqst1)/(3.*pow2(Msq)) + (
        5590*pow2(Dmsqst1))/(3.*pow4(Msq)) + (pow2(Dmst12)*(1048 + (1080*
        Dmsqst1)/pow2(Msq) + (2125*pow2(Dmsqst1))/pow4(Msq) - (32*pow2(Mgl)*(
        1327 + (2565*Dmsqst1)/pow2(Msq) + (2700*pow2(Dmsqst1))/pow4(Msq) + (
        2835*pow3(Dmsqst1))/pow6(Msq)))/pow2(Mst1) - (8*pow4(Mgl)*(80558 - (
        55260*Dmsqst1*(pow2(Dmsqst1) + Dmsqst1*pow2(Msq) + pow4(Msq)))/pow6(
        Msq)))/pow4(Mst1) + (3170*pow3(Dmsqst1))/pow6(Msq)))/(24.*pow4(Mst2)) +
        (2500*pow3(Dmsqst1) - (8*pow2(Mgl)*(270*pow2(Dmsqst1)*pow2(Msq) + 540*
        Dmsqst1*pow4(Msq) + 1909*pow6(Msq)))/(3.*pow2(Mst1)) - (4*pow4(Mgl)*(
        8055*pow2(Dmsqst1)*pow2(Msq) + 9900*pow3(Dmsqst1) + 6210*Dmsqst1*pow4(
        Msq) + 6677*pow6(Msq)))/(3.*pow4(Mst1)))/pow6(Msq) + (Dmst12*(784 + (
        1710*Dmsqst1)/pow2(Msq) + (3285*pow2(Dmsqst1))/pow4(Msq) + (4860*pow3(
        Dmsqst1))/pow6(Msq) - 4*((pow4(Mgl)*(13484 + (13815*Dmsqst1*(pow2(
        Dmsqst1) + Dmsqst1*pow2(Msq) + pow4(Msq)))/pow6(Msq)))/pow4(Mst1) + (
        pow2(Mgl)*(1485*pow2(Dmsqst1)*pow2(Msq) + 1350*pow3(Dmsqst1) + 1620*
        Dmsqst1*pow4(Msq) + 3049*pow6(Msq)))/(pow2(Mst1)*pow6(Msq)))))/(3.*
        pow2(Mst2)) - (pow3(Dmst12)*(8972 + (26570*pow2(Dmsqst1))/pow4(Msq) -
        48*((pow4(Mgl)*(53763 - (13815*Dmsqst1*(pow2(Dmsqst1) + Dmsqst1*pow2(
        Msq) + pow4(Msq)))/pow6(Msq)))/pow4(Mst1) + (15*pow2(Mgl)*(459*pow2(
        Dmsqst1)*pow2(Msq) + 468*pow3(Dmsqst1) + 450*Dmsqst1*pow4(Msq) - 76*
        pow6(Msq)))/(pow2(Mst1)*pow6(Msq))) + (5*(9721*pow3(Dmsqst1) + 907*
        Dmsqst1*pow4(Msq)))/pow6(Msq)))/(36.*pow6(Mst2))) - (6*Mt*(320*pow3(
        Dmsqst1)*pow4(Mst1) - (pow2(Dmst12)*pow2(s2t)*(30*Dmsqst1*pow4(Msq)*(-
        8*pow2(Mst2)*(36*pow2(Mgl)*pow2(Mst1) + 58*pow4(Mgl) + 5*pow4(Mst1)) +
        Dmst12*(-48*pow2(Mgl)*pow2(Mst1) + 928*pow4(Mgl) + 9*pow4(Mst1))) + 60*
        pow3(Dmsqst1)*(Dmst12*(360*pow2(Mgl)*pow2(Mst1) + 464*pow4(Mgl) + 29*
        pow4(Mst1)) - pow2(Mst2)*(336*pow2(Mgl)*pow2(Mst1) + 232*pow4(Mgl) +
        67*pow4(Mst1))) + 15*pow2(Dmsqst1)*pow2(Msq)*(Dmst12*(672*pow2(Mgl)*
        pow2(Mst1) + 1856*pow4(Mgl) + 67*pow4(Mst1)) - 2*pow2(Mst2)*(480*pow2(
        Mgl)*pow2(Mst1) + 464*pow4(Mgl) + 87*pow4(Mst1))) - 8*(Dmst12*(-6*pow2(
        Mgl)*pow2(Mst1) + 1695*pow4(Mgl) - 119*pow4(Mst1)) + 6*pow2(Mst2)*(175*
        pow2(Mgl)*pow2(Mst1) + 227*pow4(Mgl) + 4*pow4(Mst1)))*pow6(Msq)))/pow6(
        Mst2)))/(pow2(Mst1)*pow6(Msq)) + ((Mgl*(32*pow3(Dmst12)*pow3(s2t)*(60*
        Dmsqst1*(Dmsqst1*pow2(Msq)*(6*pow2(Mgl) + 5*pow2(Mst1)) + pow2(Dmsqst1)
        *(6*pow2(Mgl) + 7*pow2(Mst1)) + 3*(2*pow2(Mgl) + pow2(Mst1))*pow4(Msq))
        + (44*pow2(Mgl) + 215*pow2(Mst1))*pow6(Msq)) - (64*Dmst12*s2t*pow2(Mt)*
        (45*pow2(Dmsqst1)*pow2(Msq)*(4*pow2(Dmst12)*(118*pow2(Mgl) + 5*pow2(
        Mst1)) - Dmst12*(128*pow2(Mgl) + 33*pow2(Mst1))*pow2(Mst2) - 36*(6*
        pow2(Mgl) + pow2(Mst1))*pow4(Mst2)) + 90*pow3(Dmsqst1)*(pow2(Dmst12)*(
        204*pow2(Mgl) + 25*pow2(Mst1)) - Dmst12*(32*pow2(Mgl) + 23*pow2(Mst1))*
        pow2(Mst2) - 20*(7*pow2(Mgl) + pow2(Mst1))*pow4(Mst2)) + 90*Dmsqst1*
        pow4(Msq)*(pow2(Dmst12)*(268*pow2(Mgl) - 5*pow2(Mst1)) - 2*Dmst12*(48*
        pow2(Mgl) + 5*pow2(Mst1))*pow2(Mst2) - 4*(19*pow2(Mgl) + 4*pow2(Mst1))*
        pow4(Mst2)) - (pow2(Dmst12)*(2722*pow2(Mgl) + 223*pow2(Mst1)) + 2*
        Dmst12*(3273*pow2(Mgl) + 148*pow2(Mst1))*pow2(Mst2) + 4*(1043*pow2(Mgl)
        + 407*pow2(Mst1))*pow4(Mst2))*pow6(Msq)))/pow2(Mst1)) + 64*z2*(-(Mgl*
        pow3(Dmst12)*pow3(s2t)*(15*pow2(Dmsqst1)*pow2(Msq)*(4*pow2(Mgl) + 5*
        pow2(Mst1)) + 60*(pow2(Mgl) + 2*pow2(Mst1))*pow3(Dmsqst1) + 30*Dmsqst1*
        (2*pow2(Mgl) + pow2(Mst1))*pow4(Msq) + (113*pow2(Mgl) - 17*pow2(Mst1))*
        pow6(Msq))) - (Dmst12*Mt*s2t*(90*Dmsqst1*pow4(Msq)*(Dmst12*pow2(Mst2)*(
        28*Mgl*Mt*pow2(Mst1) - 3*s2t*pow2(Mgl)*pow2(Mst1) + 144*Mt*pow3(Mgl) -
        5*s2t*pow4(Mgl) - 2*s2t*pow4(Mst1)) + pow2(Dmst12)*(4*Mgl*Mt*pow2(Mst1)
        - 3*s2t*pow2(Mgl)*pow2(Mst1) - 360*Mt*pow3(Mgl) + 10*s2t*pow4(Mgl) +
        s2t*pow4(Mst1)) + 36*Mgl*Mt*(2*pow2(Mgl) + pow2(Mst1))*pow4(Mst2)) +
        45*pow2(Dmsqst1)*pow2(Msq)*(Dmst12*pow2(Mst2)*(88*Mgl*Mt*pow2(Mst1) -
        15*s2t*pow2(Mgl)*pow2(Mst1) + 88*Mt*pow3(Mgl) - 10*s2t*pow4(Mgl) - 6*
        s2t*pow4(Mst1)) + pow2(Dmst12)*(-104*Mgl*Mt*pow2(Mst1) + 12*s2t*pow2(
        Mgl)*pow2(Mst1) - 520*Mt*pow3(Mgl) + 20*s2t*pow4(Mgl) + 3*s2t*pow4(
        Mst1)) + 8*Mgl*Mt*(43*pow2(Mgl) + 15*pow2(Mst1))*pow4(Mst2)) + 90*pow3(
        Dmsqst1)*(pow2(Dmst12)*(-108*Mgl*Mt*pow2(Mst1) + 15*s2t*pow2(Mgl)*pow2(
        Mst1) - 160*Mt*pow3(Mgl) + 10*s2t*pow4(Mgl) + 2*s2t*pow4(Mst1)) -
        Dmst12*pow2(Mst2)*(-60*Mgl*Mt*pow2(Mst1) + 12*s2t*pow2(Mgl)*pow2(Mst1)
        + 56*Mt*pow3(Mgl) + 5*s2t*pow4(Mgl) + 4*s2t*pow4(Mst1)) + 4*Mgl*Mt*(68*
        pow2(Mgl) + 21*pow2(Mst1))*pow4(Mst2)) + (Dmst12*pow2(Mst2)*(1214*Mgl*
        Mt*pow2(Mst1) + 531*s2t*pow2(Mgl)*pow2(Mst1) + 2412*Mt*pow3(Mgl) +
        1164*s2t*pow4(Mgl) - 90*s2t*pow4(Mst1)) + pow2(Dmst12)*(384*Mgl*Mt*
        pow2(Mst1) - 666*s2t*pow2(Mgl)*pow2(Mst1) + 4704*Mt*pow3(Mgl) - 726*
        s2t*pow4(Mgl) + 45*s2t*pow4(Mst1)) + 4*Mgl*Mt*(-2143*pow2(Mgl) + 197*
        pow2(Mst1))*pow4(Mst2))*pow6(Msq)))/pow2(Mst1) + (pow3(Mt)*(30*Dmsqst1*
        pow4(Msq)*(15*pow2(Dmst12)*pow2(Mst2)*(21*pow2(Mgl)*pow2(Mst1) - 48*
        pow4(Mgl) + pow4(Mst1)) + 2*pow3(Dmst12)*(-441*pow2(Mgl)*pow2(Mst1) +
        360*pow4(Mgl) + pow4(Mst1)) + 18*Dmst12*(14*pow2(Mgl)*pow2(Mst1) + 40*
        pow4(Mgl) + pow4(Mst1))*pow4(Mst2) + 12*(12*pow2(Mgl)*pow2(Mst1) + 20*
        pow4(Mgl) + pow4(Mst1))*pow6(Mst2)) + 30*pow3(Dmsqst1)*(pow3(Dmst12)*(-
        576*pow2(Mgl)*pow2(Mst1) + 720*pow4(Mgl) - 46*pow4(Mst1)) + 3*pow2(
        Dmst12)*pow2(Mst2)*(3*pow2(Mgl)*pow2(Mst1) - 240*pow4(Mgl) + 10*pow4(
        Mst1)) + 18*Dmst12*(31*pow2(Mgl)*pow2(Mst1) + 40*pow4(Mgl) + 2*pow4(
        Mst1))*pow4(Mst2) + 12*(27*pow2(Mgl)*pow2(Mst1) + 75*pow4(Mgl) + 2*
        pow4(Mst1))*pow6(Mst2)) + 15*pow2(Dmsqst1)*pow2(Msq)*(2*pow3(Dmst12)*(-
        729*pow2(Mgl)*pow2(Mst1) + 720*pow4(Mgl) - 22*pow4(Mst1)) - 9*pow2(
        Dmst12)*pow2(Mst2)*(-36*pow2(Mgl)*pow2(Mst1) + 160*pow4(Mgl) - 5*pow4(
        Mst1)) + 18*Dmst12*(45*pow2(Mgl)*pow2(Mst1) + 80*pow4(Mgl) + 3*pow4(
        Mst1))*pow4(Mst2) + 12*(39*pow2(Mgl)*pow2(Mst1) + 95*pow4(Mgl) + 3*
        pow4(Mst1))*pow6(Mst2)) + pow6(Msq)*(pow3(Dmst12)*(1858*pow2(Mgl)*pow2(
        Mst1) - 9452*pow4(Mgl) + 34*pow4(Mst1)) + pow2(Dmst12)*pow2(Mst2)*(
        4904*pow2(Mgl)*pow2(Mst1) + 8991*pow4(Mgl) + 230*pow4(Mst1)) + 2*
        Dmst12*(1358*pow2(Mgl)*pow2(Mst1) - 4265*pow4(Mgl) + 180*pow4(Mst1))*
        pow4(Mst2) + 2*(496*pow2(Mgl)*pow2(Mst1) - 5207*pow4(Mgl) + 223*pow4(
        Mst1))*pow6(Mst2))))/pow4(Mst1)))/pow6(Msq) - 3*z3*(-32*pow3(Dmst12)*(-
        7*Mgl*pow2(Mst1) + 162*pow3(Mgl))*pow3(s2t) + 144*Mt*pow2(Dmst12)*pow2(
        s2t)*(Dmst12*((8*pow2(Mgl))/3. + (62*pow4(Mgl))/pow2(Mst1) - (pow2(
        Mst1)*(525*pow2(Dmsqst1)*pow2(Msq) + 980*pow3(Dmsqst1) + 70*Dmsqst1*
        pow4(Msq) - 1024*pow6(Msq)))/(48.*pow6(Msq))) - pow2(Mst2)*((8*pow2(
        Mgl))/3. + (124*pow4(Mgl))/pow2(Mst1) - (pow2(Mst1)*(525*pow2(Dmsqst1)*
        pow2(Msq) + 770*pow3(Dmsqst1) + 280*Dmsqst1*pow4(Msq) - 248*pow6(Msq)))
        /(24.*pow6(Msq)))) + (pow2(Mt)*(512*Dmst12*Mgl*s2t*pow2(Mst1)*(60*pow2(
        Dmsqst1)*pow2(Msq)*(2*pow2(Dmst12)*(7*pow2(Mgl) + 2*pow2(Mst1)) -
        Dmst12*(2*pow2(Mgl) + 3*pow2(Mst1))*pow2(Mst2) - 2*(5*pow2(Mgl) + 2*
        pow2(Mst1))*pow4(Mst2)) + 120*(Dmsqst1*pow4(Msq)*(10*pow2(Dmst12)*pow2(
        Mgl) - Dmst12*(4*pow2(Mgl) + pow2(Mst1))*pow2(Mst2) - (2*pow2(Mgl) +
        pow2(Mst1))*pow4(Mst2)) + pow3(Dmsqst1)*(4*pow2(Dmst12)*(pow2(Mgl) +
        pow2(Mst1)) + 2*Dmst12*(pow2(Mgl) - pow2(Mst1))*pow2(Mst2) - (8*pow2(
        Mgl) + 3*pow2(Mst1))*pow4(Mst2))) - (pow2(Dmst12)*(154*pow2(Mgl) + 13*
        pow2(Mst1)) + 2*Dmst12*(60*pow2(Mgl) + 29*pow2(Mst1))*pow2(Mst2) + (-
        95*pow2(Mgl) + 4*pow2(Mst1))*pow4(Mst2))*pow6(Msq)) + Mt*(-8*pow6(Msq)*
        (pow3(Dmst12)*(-4480*pow2(Mgl)*pow2(Mst1) + 70144*pow4(Mgl) + 375*pow4(
        Mst1)) - pow2(Dmst12)*pow2(Mst2)*(14080*pow2(Mgl)*pow2(Mst1) + 35536*
        pow4(Mgl) + 1141*pow4(Mst1)) + 32*Dmst12*(-199*pow2(Mgl)*pow2(Mst1) +
        29*pow4(Mgl) - 57*pow4(Mst1))*pow4(Mst2) + 16*(-76*pow2(Mgl)*pow2(Mst1)
        + 658*pow4(Mgl) - 395*pow4(Mst1))*pow6(Mst2)) + 15*pow2(Dmsqst1)*pow2(
        Msq)*(4*pow3(Dmst12)*(-8448*pow2(Mgl)*pow2(Mst1) + 7680*pow4(Mgl) -
        449*pow4(Mst1)) + 3*pow2(Dmst12)*pow2(Mst2)*(2048*pow2(Mgl)*pow2(Mst1)
        - 10240*pow4(Mgl) + 449*pow4(Mst1)) + 24*Dmst12*(896*pow2(Mgl)*pow2(
        Mst1) + 1280*pow4(Mgl) + 107*pow4(Mst1))*pow4(Mst2) + 16*(768*pow2(Mgl)
        *pow2(Mst1) + 1600*pow4(Mgl) + 171*pow4(Mst1))*pow6(Mst2)) + 30*(
        Dmsqst1*pow4(Msq)*(3*pow3(Dmst12)*(-7168*pow2(Mgl)*pow2(Mst1) + 5120*
        pow4(Mgl) + 7*pow4(Mst1)) + 4*pow2(Dmst12)*pow2(Mst2)*(1920*pow2(Mgl)*
        pow2(Mst1) - 3840*pow4(Mgl) + 107*pow4(Mst1)) + 8*Dmst12*(768*pow2(Mgl)
        *pow2(Mst1) + 1920*pow4(Mgl) + 107*pow4(Mst1))*pow4(Mst2) + 1024*(3*
        pow2(Mgl)*pow2(Mst1) + 5*pow4(Mgl) + pow4(Mst1))*pow6(Mst2)) + pow3(
        Dmsqst1)*(pow3(Dmst12)*(-12288*pow2(Mgl)*pow2(Mst1) + 15360*pow4(Mgl) -
        1817*pow4(Mst1)) - pow2(Dmst12)*pow2(Mst2)*(1536*pow2(Mgl)*pow2(Mst1) +
        15360*pow4(Mgl) - 919*pow4(Mst1)) + 16*(Dmst12*(960*pow2(Mgl)*(pow2(
        Mgl) + pow2(Mst1)) + 107*pow4(Mst1))*pow4(Mst2) + (576*pow2(Mgl)*pow2(
        Mst1) + 1280*pow4(Mgl) + 107*pow4(Mst1))*pow6(Mst2)))))))/(pow4(Mst1)*
        pow6(Msq))))/pow6(Mst2)))/432.)/pow4(Mt)/pow2(Sbeta)*12.;
 
   return result;
}

/**
 * 	@return returns the susy log^0 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double himalaya::H9q2::calc_coef_at_as2_no_sm_logs_log0(){
   
   const double result =
      ((Mt*pow2(Sbeta)*(-432*pow2(log(pow2(Mgl)/pow2(Mst1)))*pow6(Msq)*(-2*Mt*
        pow2(Dmst12)*pow2(Mst2)*(-72*Mt*s2t*pow2(Mst1) + 202*Mgl*pow2(Mt) +
        111*Mgl*pow2(Mst1)*pow2(s2t))*pow3(Mgl) + pow3(Dmst12)*pow3(Mgl)*(-48*
        s2t*pow2(Mst1)*pow2(Mt) + 111*Mgl*Mt*pow2(Mst1)*pow2(s2t) + 1568*Mgl*
        pow3(Mt) + 24*pow3(s2t)*pow4(Mst1)) - 8*Dmst12*pow2(Mgl)*(95*Mt*pow2(
        Mgl) - 18*Mt*pow2(Mst1) - 81*Mgl*s2t*pow2(Mst1))*pow2(Mt)*pow4(Mst2) -
        8*pow3(Mt)*(-36*pow2(Mgl)*pow2(Mst1) + 77*pow4(Mgl) + 4*pow4(Mst1))*
        pow6(Mst2)) + 192*pow2(Mt)*pow3(log(pow2(Mgl)/pow2(Mst1)))*pow6(Msq)*(-
        95*Mt*pow2(Dmst12)*pow2(Mst2)*pow4(Mgl) + 504*Mt*pow3(Dmst12)*pow4(Mgl)
        + 2*Dmst12*(-157*Mgl*Mt + 54*s2t*pow2(Mst1))*pow3(Mgl)*pow4(Mst2) + 2*
        Mt*(-157*pow4(Mgl) + 18*pow4(Mst1))*pow6(Mst2)) + 5*pow2(Dmsqst1)*pow2(
        Msq)*(-3*Mt*pow2(Dmst12)*pow2(Mst2)*(-1728*pow2(Mgl)*pow2(Mst1)*(4*(-5
        + 9*z2 - 8*z3)*pow2(Mt) + 5*(-2 + z2)*pow2(Mst1)*pow2(s2t)) + 4608*Mt*
        s2t*(11*z2 - 8*(2 + z3))*pow2(Mst1)*pow3(Mgl) + 576*((-307 + 480*z2 -
        480*z3)*pow2(Mt) + (29 - 10*z2)*pow2(Mst1)*pow2(s2t))*pow4(Mgl) + 576*
        Mgl*Mt*s2t*(-33 + 88*z2 - 96*z3)*pow4(Mst1) + ((-850 - 8640*z2 + 12123*
        z3)*pow2(Mt) - 54*(-58 + 64*z2 - 35*z3)*pow2(Mst1)*pow2(s2t))*pow4(
        Mst1)) + pow3(Dmst12)*(-5184*Mt*pow2(Mgl)*pow2(Mst1)*((-51 + 162*z2 -
        176*z3)*pow2(Mt) + (-7 + 4*z2)*pow2(Mst1)*pow2(s2t)) - 2304*s2t*pow2(
        Mst1)*(-6*(-59 + 65*z2 - 56*z3)*pow2(Mt) + (-3 + z2)*pow2(Mst1)*pow2(
        s2t))*pow3(Mgl) + 1728*Mt*((-307 + 480*z2 - 480*z3)*pow2(Mt) + 2*(29 -
        10*z2)*pow2(Mst1)*pow2(s2t))*pow4(Mgl) - 576*Mgl*s2t*((60 - 312*z2 +
        384*z3)*pow2(Mt) + 5*(-2 + z2)*pow2(Mst1)*pow2(s2t))*pow4(Mst1) + Mt*(-
        4*(5314 + 6336*z2 - 12123*z3)*pow2(Mt) + 27*(134 - 192*z2 + 105*z3)*
        pow2(Mst1)*pow2(s2t))*pow4(Mst1)) + 216*Dmst12*pow2(Mt)*(24*Mt*(-11 +
        90*z2 - 112*z3)*pow2(Mgl)*pow2(Mst1) - 64*s2t*(-27 + 43*z2 - 40*z3)*
        pow2(Mst1)*pow3(Mgl) + 8*Mt*(-307 + 480*z2 - 480*z3)*pow4(Mgl) + Mt*(
        146 + 144*z2 - 321*z3)*pow4(Mst1) - 32*Mgl*s2t*(-9 + 30*z2 - 32*z3)*
        pow4(Mst1))*pow4(Mst2) + 48*pow3(Mt)*(432*(-1 + 13*z2 - 16*z3)*pow2(
        Mgl)*pow2(Mst1) + 36*(-179 + 380*z2 - 400*z3)*pow4(Mgl) + (1118 + 432*
        z2 - 1539*z3)*pow4(Mst1))*pow6(Mst2)) + 10*Dmsqst1*pow4(Msq)*(pow3(
        Dmst12)*(2592*Mt*pow2(Mgl)*pow2(Mst1)*((50 - 196*z2 + 224*z3)*pow2(Mt)
        + (-1 + 2*z2)*pow2(Mst1)*pow2(s2t)) - 1152*s2t*pow2(Mst1)*(-6*(-67 +
        90*z2 - 80*z3)*pow2(Mt) + (-3 + z2)*pow2(Mst1)*pow2(s2t))*pow3(Mgl) +
        864*Mt*((-307 + 480*z2 - 480*z3)*pow2(Mt) + 2*(29 - 10*z2)*pow2(Mst1)*
        pow2(s2t))*pow4(Mgl) - 576*Mgl*s2t*(3*(-5 + 4*z2)*pow2(Mt) + (-3 + z2)*
        pow2(Mst1)*pow2(s2t))*pow4(Mst1) + Mt*((-1814 + 1152*z2 - 567*z3)*pow2(
        Mt) + 27*(18 - 64*z2 + 7*z3)*pow2(Mst1)*pow2(s2t))*pow4(Mst1)) - 108*
        Mt*pow2(Dmst12)*pow2(Mst2)*(-24*pow2(Mgl)*pow2(Mst1)*((-19 + 70*z2 -
        80*z3)*pow2(Mt) + 2*(-3 + z2)*pow2(Mst1)*pow2(s2t)) + 256*Mt*s2t*(-6 +
        9*z2 - 8*z3)*pow2(Mst1)*pow3(Mgl) + 8*((-307 + 480*z2 - 480*z3)*pow2(
        Mt) + (29 - 10*z2)*pow2(Mst1)*pow2(s2t))*pow4(Mgl) + 32*Mgl*Mt*s2t*(-5
        + 14*z2 - 16*z3)*pow4(Mst1) + ((-6 - 80*z2 + 107*z3)*pow2(Mt) + 2*(10 -
        16*z2 + 7*z3)*pow2(Mst1)*pow2(s2t))*pow4(Mst1)) + 216*Dmst12*pow2(Mt)*(
        48*Mt*(-3 + 14*z2 - 16*z3)*pow2(Mgl)*pow2(Mst1) + 32*s2t*(19 - 18*z2 +
        16*z3)*pow2(Mst1)*pow3(Mgl) + 4*Mt*(-307 + 480*z2 - 480*z3)*pow4(Mgl) +
        Mt*(38 + 48*z2 - 107*z3)*pow4(Mst1) - 32*Mgl*s2t*(-4 + 9*z2 - 8*z3)*
        pow4(Mst1))*pow4(Mst2) + 192*pow3(Mt)*(108*(-1 + 4*z2 - 4*z3)*pow2(Mgl)
        *pow2(Mst1) + 9*(-69 + 80*z2 - 80*z3)*pow4(Mgl) + 4*(23 + 9*z2 - 36*z3)
        *pow4(Mst1))*pow6(Mst2)) + 8*pow6(Msq)*(3*Mt*pow2(Dmst12)*pow2(Mst2)*(-
        4*pow2(Mgl)*pow2(Mst1)*((2654 - 9808*z2 + 10560*z3)*pow2(Mt) + 9*(175 +
        118*z2 - 4*z3)*pow2(Mst1)*pow2(s2t)) - 48*Mt*s2t*(-1091 + 402*z2 - 480*
        z3)*pow2(Mst1)*pow3(Mgl) + 4*((-40279 + 17982*z2 - 26652*z3)*pow2(Mt) -
        3*(681 + 776*z2 - 558*z3)*pow2(Mst1)*pow2(s2t))*pow4(Mgl) - 16*Mgl*Mt*
        s2t*(-148 + 607*z2 - 696*z3)*pow4(Mst1) + ((262 + 1840*z2 - 3423*z3)*
        pow2(Mt) + 18*(-8 + 40*z2 + 31*z3)*pow2(Mst1)*pow2(s2t))*pow4(Mst1)) +
        pow3(Dmst12)*(12*Mt*pow2(Mgl)*pow2(Mst1)*(4*(-570 + 929*z2 - 840*z3)*
        pow2(Mt) + 9*(1 + 148*z2 - 4*z3)*pow2(Mst1)*pow2(s2t)) - 24*s2t*pow2(
        Mst1)*((-2722 + 4704*z2 - 3696*z3)*pow2(Mt) + (-22 + 113*z2 - 243*z3)*
        pow2(Mst1)*pow2(s2t))*pow3(Mgl) - 6*Mt*(4*(-53763 + 9452*z2 - 26304*z3)
        *pow2(Mt) + 3*(1695 - 968*z2 + 558*z3)*pow2(Mst1)*pow2(s2t))*pow4(Mgl)
        + 12*Mgl*s2t*((446 - 768*z2 + 624*z3)*pow2(Mt) + (215 + 34*z2 - 21*z3)*
        pow2(Mst1)*pow2(s2t))*pow4(Mst1) + Mt*((-4486 + 816*z2 + 3375*z3)*pow2(
        Mt) - 18*(-119 + 60*z2 + 192*z3)*pow2(Mst1)*pow2(s2t))*pow4(Mst1)) -
        24*Dmst12*pow2(Mt)*(Mt*(3049 - 2716*z2 + 2388*z3)*pow2(Mgl)*pow2(Mst1)
        - 4*s2t*(1043 + 2143*z2 - 570*z3)*pow2(Mst1)*pow3(Mgl) + 2*Mt*(6742 +
        4265*z2 - 174*z3)*pow4(Mgl) - 4*Mt*(49 + 90*z2 - 171*z3)*pow4(Mst1) +
        4*Mgl*s2t*(-407 + 197*z2 - 24*z3)*pow4(Mst1))*pow4(Mst2) - 12*pow3(Mt)*
        (4*(1909 - 496*z2 + 228*z3)*pow2(Mgl)*pow2(Mst1) + 2*(6677 + 10414*z2 -
        3948*z3)*pow4(Mgl) + (-4125 - 892*z2 + 4740*z3)*pow4(Mst1))*pow6(Mst2))
        - 10*pow3(Dmsqst1)*(18*pow6(Mst1)*(3*Mt*(134 - 128*z2 + 77*z3)*pow2(
        Dmst12)*pow2(Mst2)*pow2(s2t) + (32*Mgl*s2t*(-7 + 4*z2) + 3*Mt*(-58 +
        64*z2 - 49*z3))*pow2(s2t)*pow3(Dmst12) + 32*Mt*pow6(Mst2)) - 864*pow3(
        Mt)*pow4(Mgl)*((307 - 480*z2 + 480*z3)*pow2(Dmst12)*pow2(Mst2) + (-307
        + 480*z2 - 480*z3)*pow3(Dmst12) + Dmst12*(-307 + 480*z2 - 480*z3)*pow4(
        Mst2) + 20*(-11 + 30*z2 - 32*z3)*pow6(Mst2)) + 864*Mt*pow2(Mgl)*pow2(
        Mst1)*(-(pow2(Dmst12)*pow2(Mst2)*(16*Mgl*Mt*s2t*(4 + 7*z2 - 8*z3) + (-
        63 + 6*z2 + 48*z3)*pow2(Mt) + (-29 + 10*z2)*pow2(Mgl)*pow2(s2t))) + 2*(
        -4*Mgl*Mt*s2t*(-51 + 40*z2 - 32*z3) + 6*(-13 + 32*z2 - 32*z3)*pow2(Mt)
        + (-29 + 10*z2)*pow2(Mgl)*pow2(s2t))*pow3(Dmst12) + 2*Dmst12*Mt*(-3*Mt*
        (-5 + 62*z2 - 80*z3) + 4*Mgl*s2t*(-35 + 68*z2 - 64*z3))*pow4(Mst2) +
        72*(-3*z2 + 4*z3)*pow2(Mt)*pow6(Mst2)) + pow4(Mst1)*(-3*Mt*pow2(Dmst12)
        *pow2(Mst2)*(-576*Mgl*Mt*s2t*(-23 + 60*z2 - 64*z3) + (634 + 5760*z2 -
        8271*z3)*pow2(Mt) + 1728*(-7 + 4*z2)*pow2(Mgl)*pow2(s2t)) + pow3(
        Dmst12)*(-1728*Mgl*s2t*(-25 + 108*z2 - 128*z3)*pow2(Mt) + 12960*Mt*(-3
        + 2*z2)*pow2(Mgl)*pow2(s2t) + (19442 + 26496*z2 - 49059*z3)*pow3(Mt) +
        1152*(-3 + z2)*pow3(Mgl)*pow3(s2t)) - 432*Dmst12*(Mt*(54 + 48*z2 - 107*
        z3) + 16*Mgl*s2t*(5 - 21*z2 + 24*z3))*pow2(Mt)*pow4(Mst2) - 144*(250 +
        96*z2 - 321*z3)*pow3(Mt)*pow6(Mst2))) - 96*log(pow2(Mgl)/pow2(Mst1))*(
        60*Dmst12*pow3(Dmsqst1)*pow3(Mgl)*(-3*Dmst12*Mt*pow2(Mst2)*(16*Mt*s2t*
        pow2(Mst1) + 15*Mgl*pow2(Mt) - 5*Mgl*pow2(Mst1)*pow2(s2t)) + pow2(
        Dmst12)*(96*s2t*pow2(Mst1)*pow2(Mt) - 30*Mgl*Mt*pow2(Mst1)*pow2(s2t) +
        45*Mgl*pow3(Mt) - 2*pow3(s2t)*pow4(Mst1)) + 45*Mgl*pow3(Mt)*pow4(Mst2))
        + 60*pow2(Dmsqst1)*pow2(Msq)*pow3(Mgl)*(-15*Mt*pow2(Dmst12)*pow2(Mst2)*
        (2*Mt*s2t*pow2(Mst1) + 3*Mgl*pow2(Mt) - Mgl*pow2(Mst1)*pow2(s2t)) +
        pow3(Dmst12)*(78*s2t*pow2(Mst1)*pow2(Mt) - 30*Mgl*Mt*pow2(Mst1)*pow2(
        s2t) + 45*Mgl*pow3(Mt) - 2*pow3(s2t)*pow4(Mst1)) + 9*Dmst12*(5*Mgl*Mt -
        2*s2t*pow2(Mst1))*pow2(Mt)*pow4(Mst2) + 15*Mgl*pow3(Mt)*pow6(Mst2)) +
        60*Dmsqst1*pow3(Mgl)*pow4(Msq)*(-3*Mt*pow2(Dmst12)*pow2(Mst2)*(4*Mt*
        s2t*pow2(Mst1) + 15*Mgl*pow2(Mt) - 5*Mgl*pow2(Mst1)*pow2(s2t)) + pow3(
        Dmst12)*(60*s2t*pow2(Mst1)*pow2(Mt) - 30*Mgl*Mt*pow2(Mst1)*pow2(s2t) +
        45*Mgl*pow3(Mt) - 2*pow3(s2t)*pow4(Mst1)) + 9*Dmst12*(5*Mgl*Mt - 4*s2t*
        pow2(Mst1))*pow2(Mt)*pow4(Mst2) + 30*Mgl*pow3(Mt)*pow6(Mst2)) + pow6(
        Msq)*(2*Mt*pow2(Dmst12)*pow2(Mgl)*pow2(Mst2)*(-102*Mgl*Mt*s2t*pow2(
        Mst1) + 2*pow2(Mgl)*((737 + 162*z2)*pow2(Mt) + 3*(328 + 27*z2)*pow2(
        Mst1)*pow2(s2t)) + 27*pow2(s2t)*pow4(Mst1)) - pow2(Mgl)*pow3(Dmst12)*(
        Mt*pow2(Mgl)*(32*(434 + 81*z2)*pow2(Mt) + 3*(605 + 54*z2)*pow2(Mst1)*
        pow2(s2t)) + 54*Mt*pow2(s2t)*pow4(Mst1) - 2*Mgl*(67*s2t*pow2(Mst1)*
        pow2(Mt) + 2*(-154 + 45*z2)*pow3(s2t)*pow4(Mst1))) + 12*Dmst12*pow2(
        Mgl)*(18*Mt*(37 + 6*z2)*pow2(Mgl) + 9*Mt*pow2(Mst1) - 476*Mgl*s2t*pow2(
        Mst1))*pow2(Mt)*pow4(Mst2) + 12*pow3(Mt)*(18*pow2(Mgl)*pow2(Mst1) + 3*(
        215 + 36*z2)*pow4(Mgl) - 52*pow4(Mst1))*pow6(Mst2)))))/(1296.*pow4(
        Mst1)*pow6(Msq)*pow6(Mst2)))/pow4(Mt)/pow2(Sbeta)*12.;

   return result;
}

/**
 * 	@return returns the susy log^1 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double himalaya::H9q2::calc_coef_at_as2_no_sm_logs_log1(){

   const double result =
      ((2*Mt*pow2(Sbeta)*(-6*pow2(Mt)*pow2(log(pow2(Mgl)/pow2(Mst1)))*pow6(Msq)*
        (-203*Mt*pow2(Dmst12)*pow2(Mst2)*pow4(Mgl) + 936*Mt*pow3(Dmst12)*pow4(
        Mgl) + 2*Dmst12*(-265*Mgl*Mt + 162*s2t*pow2(Mst1))*pow3(Mgl)*pow4(Mst2)
        + 2*Mt*(-265*pow4(Mgl) + 18*pow4(Mst1))*pow6(Mst2)) + 10*Mt*pow3(
        Dmsqst1)*(2*pow3(Dmst12)*(-108*(-7 + 4*z2)*pow2(Mgl)*pow2(Mst1)*pow2(
        Mt) + 72*Mt*s2t*(-11 + 4*z2)*pow2(Mst1)*pow3(Mgl) + 9*(-119 + 60*z2)*
        pow2(Mt)*pow4(Mgl) + 6*Mgl*Mt*s2t*(-77 + 48*z2)*pow4(Mst1) + ((121 -
        72*z2)*pow2(Mt) - 9*pow2(Mst1)*pow2(s2t))*pow4(Mst1)) - 3*pow2(Dmst12)*
        pow2(Mst2)*(9*(-5 + 4*z2)*pow2(Mgl)*pow2(Mst1)*pow2(Mt) + 48*Mt*s2t*(1
        - 2*z2)*pow2(Mst1)*pow3(Mgl) + 6*(-119 + 60*z2)*pow2(Mt)*pow4(Mgl) +
        12*Mgl*Mt*s2t*(-13 + 8*z2)*pow4(Mst1) - 2*(4*(-5 + 3*z2)*pow2(Mt) + 3*
        pow2(Mst1)*pow2(s2t))*pow4(Mst1)) + 18*Dmst12*Mt*(3*Mt*(-33 + 20*z2)*
        pow2(Mgl)*pow2(Mst1) - 8*s2t*(-13 + 8*z2)*pow2(Mst1)*pow3(Mgl) + Mt*(-
        119 + 60*z2)*pow4(Mgl) - 8*Mgl*s2t*(-5 + 3*z2)*pow4(Mst1) + Mt*(-15 +
        8*z2)*pow4(Mst1))*pow4(Mst2) + 36*pow2(Mt)*(6*(-5 + 3*z2)*pow2(Mgl)*
        pow2(Mst1) + 5*(-13 + 8*z2)*pow4(Mgl) + (-9 + 4*z2)*pow4(Mst1))*pow6(
        Mst2)) + 10*Dmsqst1*Mt*pow4(Msq)*(-9*pow2(Dmst12)*pow2(Mst2)*(3*(33 -
        20*z2)*pow2(Mgl)*pow2(Mst1)*pow2(Mt) + 16*Mt*s2t*(-7 + 4*z2)*pow2(Mst1)
        *pow3(Mgl) + 2*(-119 + 60*z2)*pow2(Mt)*pow4(Mgl) + 8*Mgl*Mt*s2t*(-3 +
        2*z2)*pow4(Mst1) + ((7 - 4*z2)*pow2(Mt) - 2*pow2(Mst1)*pow2(s2t))*pow4(
        Mst1)) + 18*Dmst12*Mt*(6*Mt*(-7 + 4*z2)*pow2(Mgl)*pow2(Mst1) - 8*s2t*(-
        5 + 2*z2)*pow2(Mst1)*pow3(Mgl) + Mt*(-119 + 60*z2)*pow4(Mgl) + 4*Mt*(-2
        + z2)*pow4(Mst1) - 8*Mgl*s2t*(-2 + z2)*pow4(Mst1))*pow4(Mst2) + 2*pow3(
        Dmst12)*(-27*(-47 + 28*z2)*pow2(Mgl)*pow2(Mst1)*pow2(Mt) + 72*Mt*s2t*(-
        19 + 10*z2)*pow2(Mst1)*pow3(Mgl) + 9*(-119 + 60*z2)*pow2(Mt)*pow4(Mgl)
        + 6*Mgl*Mt*s2t*pow4(Mst1) + pow2(Mt)*pow4(Mst1) - 9*pow2(s2t)*pow6(
        Mst1)) + 6*pow2(Mt)*(36*(-2 + z2)*pow2(Mgl)*pow2(Mst1) + 3*(-49 + 20*
        z2)*pow4(Mgl) + 2*(-17 + 6*z2)*pow4(Mst1))*pow6(Mst2)) + 5*Mt*pow2(
        Dmsqst1)*pow2(Msq)*(-3*pow2(Dmst12)*pow2(Mst2)*(-36*(-7 + 4*z2)*pow2(
        Mgl)*pow2(Mst1)*pow2(Mt) + 96*Mt*s2t*(-3 + z2)*pow2(Mst1)*pow3(Mgl) +
        12*(-119 + 60*z2)*pow2(Mt)*pow4(Mgl) + 12*Mgl*Mt*s2t*(-19 + 12*z2)*
        pow4(Mst1) + ((61 - 36*z2)*pow2(Mt) - 12*pow2(Mst1)*pow2(s2t))*pow4(
        Mst1)) + 2*pow3(Dmst12)*(-27*(-75 + 44*z2)*pow2(Mgl)*pow2(Mst1)*pow2(
        Mt) + 144*Mt*s2t*(-15 + 7*z2)*pow2(Mst1)*pow3(Mgl) + 18*(-119 + 60*z2)*
        pow2(Mt)*pow4(Mgl) + 24*Mgl*Mt*s2t*(-19 + 12*z2)*pow4(Mst1) - 2*((-61 +
        36*z2)*pow2(Mt) + 9*pow2(Mst1)*pow2(s2t))*pow4(Mst1)) + 18*Dmst12*Mt*(
        3*Mt*(-47 + 28*z2)*pow2(Mgl)*pow2(Mst1) - 16*s2t*(-9 + 5*z2)*pow2(Mst1)
        *pow3(Mgl) + 2*Mt*(-119 + 60*z2)*pow4(Mgl) - 8*Mgl*s2t*(-7 + 4*z2)*
        pow4(Mst1) + Mt*(-23 + 12*z2)*pow4(Mst1))*pow4(Mst2) + 6*pow2(Mt)*(36*(
        -7 + 4*z2)*pow2(Mgl)*pow2(Mst1) + (-537 + 300*z2)*pow4(Mgl) + 4*(-22 +
        9*z2)*pow4(Mst1))*pow6(Mst2)) + pow6(Msq)*(Mt*pow2(Dmst12)*pow2(Mst2)*(
        -6*pow2(Mgl)*pow2(Mst1)*((659 - 440*z2)*pow2(Mt) + 243*pow2(Mst1)*pow2(
        s2t)) - 72*Mt*s2t*(-99 + 20*z2)*pow2(Mst1)*pow3(Mgl) + ((-23651 + 5436*
        z2)*pow2(Mt) - 582*pow2(Mst1)*pow2(s2t))*pow4(Mgl) + 12*Mgl*Mt*s2t*(25
        - 62*z2)*pow4(Mst1) + 6*((-145 + 46*z2)*pow2(Mt) + 79*pow2(Mst1)*pow2(
        s2t))*pow4(Mst1)) + pow3(Dmst12)*(3*Mt*pow2(Mgl)*pow2(Mst1)*((-498 +
        280*z2)*pow2(Mt) + 379*pow2(Mst1)*pow2(s2t)) - 24*s2t*pow2(Mst1)*((-86
        + 77*z2)*pow2(Mt) + pow2(Mst1)*pow2(s2t))*pow3(Mgl) + (-315*Mt*pow2(
        Mst1)*pow2(s2t) + (63184 - 7104*z2)*pow3(Mt))*pow4(Mgl) + 2*Mgl*s2t*((
        442 - 60*z2)*pow2(Mt) + 59*pow2(Mst1)*pow2(s2t))*pow4(Mst1) - 3*Mt*(4*(
        -43 + 8*z2)*pow2(Mt) + 125*pow2(Mst1)*pow2(s2t))*pow4(Mst1)) - 6*
        Dmst12*pow2(Mt)*(4*Mt*(103 - 50*z2)*pow2(Mgl)*pow2(Mst1) - 4*s2t*(283 +
        138*z2)*pow2(Mst1)*pow3(Mgl) + Mt*(2647 + 628*z2)*pow4(Mgl) + 12*Mt*(8
        - 5*z2)*pow4(Mst1) - 4*Mgl*s2t*(65 + 2*z2)*pow4(Mst1))*pow4(Mst2) - 12*
        pow3(Mt)*(2*(91 - 10*z2)*pow2(Mgl)*pow2(Mst1) + (589 + 464*z2)*pow4(
        Mgl) - 4*(-35 + 15*z2 + 3*z3)*pow4(Mst1))*pow6(Mst2)) - 2*log(pow2(Mgl)
        /pow2(Mst1))*(180*Dmst12*pow2(Mt)*pow3(Dmsqst1)*pow3(Mgl)*(pow2(Dmst12)
        *(5*Mgl*Mt + 8*s2t*pow2(Mst1)) - Dmst12*(5*Mgl*Mt + 4*s2t*pow2(Mst1))*
        pow2(Mst2) + 5*Mgl*Mt*pow4(Mst2)) + 180*Dmsqst1*pow2(Mt)*pow3(Mgl)*
        pow4(Msq)*(-5*Mgl*Mt*pow2(Dmst12)*pow2(Mst2) + (5*Mgl*Mt + 4*s2t*pow2(
        Mst1))*pow3(Dmst12) + Dmst12*(5*Mgl*Mt - 4*s2t*pow2(Mst1))*pow4(Mst2) +
        5*Mgl*Mt*pow6(Mst2)) + 90*pow2(Dmsqst1)*pow2(Msq)*pow2(Mt)*pow3(Mgl)*(-
        2*pow2(Dmst12)*(5*Mgl*Mt + 2*s2t*pow2(Mst1))*pow2(Mst2) + 2*(5*Mgl*Mt +
        6*s2t*pow2(Mst1))*pow3(Dmst12) + 2*Dmst12*(5*Mgl*Mt - 2*s2t*pow2(Mst1))
        *pow4(Mst2) + 5*Mgl*Mt*pow6(Mst2)) + pow6(Msq)*(6*Mt*pow2(Dmst12)*pow2(
        Mst2)*(14*Mt*s2t*pow2(Mst1) + 67*Mgl*pow2(Mt) + 117*Mgl*pow2(Mst1)*
        pow2(s2t))*pow3(Mgl) - pow3(Dmst12)*pow3(Mgl)*(42*s2t*pow2(Mst1)*pow2(
        Mt) + 351*Mgl*(Mt*pow2(Mst1)*pow2(s2t) + 8*pow3(Mt)) + 68*pow3(s2t)*
        pow4(Mst1)) + 12*Dmst12*pow2(Mgl)*(167*Mt*pow2(Mgl) - 54*Mt*pow2(Mst1)
        - 170*Mgl*s2t*pow2(Mst1))*pow2(Mt)*pow4(Mst2) + 36*pow3(Mt)*(-36*pow2(
        Mgl)*pow2(Mst1) + 59*pow4(Mgl) + 4*pow4(Mst1))*pow6(Mst2)))))/(27.*
        pow4(Mst1)*pow6(Msq)*pow6(Mst2)))/pow4(Mt)/pow2(Sbeta)*12.; 

   return result;
}

/**
 * 	@return returns the susy log^2 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double himalaya::H9q2::calc_coef_at_as2_no_sm_logs_log2(){

   const double result =
      ((2*Mt*pow2(Sbeta)*(Mt*pow2(Dmst12)*pow2(Mst2)*pow4(Msq)*(528*Mt*s2t*pow2(
        Mst1)*pow3(Mgl) + (-1966*pow2(Mt) + 96*pow2(Mst1)*pow2(s2t))*pow4(Mgl)
        + 180*Mgl*Mt*s2t*pow4(Mst1) + 3*(51*pow2(Mt) + 82*pow2(Mst1)*pow2(s2t))
        *pow4(Mst1) + pow2(Mgl)*(64*pow2(Mst1)*pow2(Mt) - 321*pow2(s2t)*pow4(
        Mst1))) - 6*Dmst12*pow2(Mt)*pow4(Msq)*(43*Mt*pow2(Mgl)*pow2(Mst1) -
        248*s2t*pow2(Mst1)*pow3(Mgl) + 550*Mt*pow4(Mgl) + 51*Mt*pow4(Mst1) +
        60*Mgl*s2t*pow4(Mst1))*pow4(Mst2) - pow3(Dmst12)*pow4(Msq)*((48*Mt*
        pow2(Mst1)*pow2(s2t) - 7232*pow3(Mt))*pow4(Mgl) + 102*pow3(Mt)*pow4(
        Mst1) + pow2(Mgl)*(32*pow2(Mst1)*pow3(Mt) - 321*Mt*pow2(s2t)*pow4(Mst1)
        ) + 16*pow3(Mgl)*(11*s2t*pow2(Mst1)*pow2(Mt) + 2*pow3(s2t)*pow4(Mst1))
        + 369*Mt*pow2(s2t)*pow6(Mst1) + Mgl*(120*s2t*pow2(Mt)*pow4(Mst1) - 41*
        pow3(s2t)*pow6(Mst1))) + 6*pow3(Mt)*(-86*pow2(Mgl)*pow2(Mst1)*pow4(Msq)
        - 380*pow4(Mgl)*pow4(Msq) + (15*pow2(Dmsqst1) + 30*Dmsqst1*pow2(Msq) +
        98*pow4(Msq))*pow4(Mst1))*pow6(Mst2) + 12*log(pow2(Mgl)/pow2(Mst1))*
        pow2(Mt)*pow4(Msq)*(-53*Mt*pow2(Dmst12)*pow2(Mst2)*pow4(Mgl) + 276*Mt*
        pow3(Dmst12)*pow4(Mgl) + 2*Dmst12*(-85*Mgl*Mt + 44*s2t*pow2(Mst1))*
        pow3(Mgl)*pow4(Mst2) + 2*Mt*(-85*pow4(Mgl) + 9*pow4(Mst1))*pow6(Mst2)))
        )/(27.*pow4(Msq)*pow4(Mst1)*pow6(Mst2)))/pow4(Mt)/pow2(Sbeta)*12.; 
   
   return result;
}

/**
 * 	@return returns the susy log^3 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double himalaya::H9q2::calc_coef_at_as2_no_sm_logs_log3(){

   const double result =
      ((-224*pow2(Sbeta)*pow4(Mt))/9.)/pow4(Mt)/pow2(Sbeta)*12.; 

   return result;
}



