#include "H6b2qg2.hpp"
#include "HierarchyCalculator.hpp"
#include "Hierarchies.hpp"
#include "Constants.hpp"
#include "Utils.hpp"
#include <cmath>
#include <type_traits>

/**
 * 	Constructor
 * 	@param flagMap the flagMap for the truncation of expansion variables
 * 	@param Al4p a double alpha_s/4/Pi
 * 	@param beta a double which is the mixing angle beta
 * 	@param Dmglst2 a double Mgl - Mst2
 * 	@param Dmsqst2 a double Msq - Mst2
 * 	@param lmMt a double log((<renormalization scale> / Mt)^2)
 * 	@param lmMst1 a double log((<renormalization scale> / Mst1)^2)
 * 	@param lmMst2 a double log((<renormalization scale> / Mst2)^2)
 * 	@param Mgl a double gluino mass
 * 	@param Mt a double top/bottom quark mass
 * 	@param Mst1 a double stop 1 mass
 * 	@param Mst2 a double stop 2 mass
 * 	@param MuSUSY a double mu parameter
 * 	@param s2t a double 2 times the sine of the stop/sbottom quark mixing angle
 * 	@param mdrFlag an int 0 for DR and 1 for MDR scheme
 * 	@param oneLoopFlag an int flag to consider the one-loop expansion terms
 * 	@param twoLoopFlag an int flag to consider the two-loop expansion terms
 * 	@param threeLoopFlag an int flag to consider the three-loop expansion terms
 */
himalaya::H6b2qg2::H6b2qg2(const std::map<unsigned int, unsigned int>& flagMap, double Al4p, double beta, double Dmglst2,
		 double Dmsqst2, double lmMt, double lmMst1, double lmMst2,
		 double Mgl, double Mt, double Mst1, double Mst2, double MuSUSY,
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
   this -> Mgl = Mgl;
   this -> Mt = Mt;
   this -> Mst1 = Mst1;
   this -> Mst2 = Mst2;
   this -> MuSUSY = MuSUSY;
   this -> s2t = s2t;

   // mdr flags, indicates if one wants to shift the dr stop mass to the mdr stop mass
   shiftst1 = mdrFlag;
   shiftst2 = mdrFlag;
   shiftst3 = mdrFlag;
   xDR2DRMOD = mdrFlag;
   // expansion flags
   xDmglst2 = flagMap.at(ExpansionDepth::xxDmglst2);
   xDmsqst2 = flagMap.at(ExpansionDepth::xxDmsqst2);
   xMst = flagMap.at(ExpansionDepth::xxMst);
   
   s1 = 
   #include "../hierarchies/h6b2qg2/sigS1Full.inc"
   ;
   s2 = 
   #include "../hierarchies/h6b2qg2/sigS2Full.inc"
   ;
   s12 = 
   #include "../hierarchies/h6b2qg2/sigS12Full.inc"
   ;
}

/**
 * 	@return The diagonal (1, 1) matrix element of the Higgs mass matrix as a double for the hierarchy 'H6b2qg2'
 */
double himalaya::H6b2qg2::getS1() const {
   return s1;
}

/**
 * 	@return The diagonal (2, 2) matrix element of the Higgs mass matrix as a double for the hierarchy 'H6b2qg2'
 */
double himalaya::H6b2qg2::getS2() const {
   return s2;
}

/**
 * 	@return The off-diagonal (1, 2) = (2, 1) matrix element of the Higgs mass matrix as a double for the hierarchy 'H6b2qg2'
 */
double himalaya::H6b2qg2::getS12() const {
   return s12;
}

/**
 * 	@return returns the susy log^0 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double himalaya::H6b2qg2::calc_coef_at_as2_no_sm_logs_log0() const {

   const double result =
     (((-257250*z3*(23328*pow4(Mt)*(913.7777777777778 + (18034*
        Dmglst2)/(27.*Mgl) + (16*log(pow2(Mst2)/pow2(Mst1)))/3. + (60*pow2(Mgl)
        *(432*Dmsqst2 + 6905*pow2(Mst1)) - 48*Dmglst2*Mgl*(2700*Dmsqst2 +
        14423*pow2(Mst1)) + pow2(Dmglst2)*(-129600*Dmsqst2 + 1785736*pow2(Mst1)
        + 781449*pow2(Mst2)))/(243.*pow2(Mgl)*pow2(Mst2)) + ((125*pow2(Dmsqst2)
        )/2. + (640*Dmsqst2*(-9*Dmglst2*Mgl - 9*pow2(Dmglst2) + pow2(Mgl))*
        pow2(Mst1))/(3.*pow2(Mgl)) + (32*(72723 - (295693*Dmglst2*(Dmglst2 +
        Mgl))/pow2(Mgl))*pow4(Mst1))/729.)/pow4(Mst2) + (160*Dmsqst2*pow2(Mst1)
        *(5*Dmsqst2 + 6*pow2(Mst1)))/(3.*pow6(Mst2))) + (s2t*(-24*s2t*pow2(Mt)*
        (2*(3*pow2(Dmglst2)*(27375*pow2(Mst1)*pow2(Mst2) + 5872*pow4(Mst1) -
        177327*pow4(Mst2)) - 2*pow2(Mgl)*(660*pow2(Mst1)*pow2(Mst2) + 648*log(
        pow2(Mst2)/pow2(Mst1))*pow2(Mst2)*(-pow2(Mst1) + pow2(Mst2)) + 32663*
        pow4(Mst1) - 6237*pow4(Mst2)) + 24*Dmglst2*Mgl*(-2223*pow2(Mst1)*pow2(
        Mst2) + 734*pow4(Mst1) - 1317*pow4(Mst2)))*pow4(Mst2) + 15*pow2(
        Dmsqst2)*pow2(Mgl)*(-7*pow2(Mst1)*pow2(Mst2) + 884*pow4(Mst1) + 903*
        pow4(Mst2)) - 80*Dmsqst2*pow2(Mst2)*(405*Dmglst2*(Dmglst2 + Mgl)*pow4(
        Mst1) + pow2(Mgl)*(126*pow2(Mst1)*pow2(Mst2) + 56*pow4(Mst1) + 783*
        pow4(Mst2)))) - pow2(Mst2)*pow3(s2t)*(240*Dmsqst2*pow2(Mst2)*(pow2(Mgl)
        *(1557*pow2(Mst1)*pow2(Mst2) + 520*pow4(Mst1) - 621*pow4(Mst2)) - 324*
        Dmglst2*(Dmglst2 + Mgl)*(-2*pow2(Mst1)*pow2(Mst2) + 5*pow4(Mst1) - 3*
        pow4(Mst2))) + 15*pow2(Dmsqst2)*pow2(Mgl)*(-11370*pow2(Mst1)*pow2(Mst2)
        + 2453*pow4(Mst1) - 39627*pow4(Mst2)) - 4*pow4(Mst2)*(pow2(Dmglst2)*(-
        740130*pow2(Mst1)*pow2(Mst2) + 434566*pow4(Mst1) - 2057463*pow4(Mst2))
        + Dmglst2*Mgl*(-202932*pow2(Mst1)*pow2(Mst2) + 434566*pow4(Mst1) -
        951210*pow4(Mst2)) - 12*pow2(Mgl)*(35931*pow2(Mst1)*pow2(Mst2) + 5225*
        pow4(Mst1) - 162*log(pow2(Mst2)/pow2(Mst1))*(pow4(Mst1) - pow4(Mst2)) +
        35613*pow4(Mst2)))) + 32*Mst2*Mt*(-4*pow2(Mt)*(19440*pow2(Dmsqst2)*
        pow2(Mgl)*(2*pow2(Mst1) + pow2(Mst2)) + 38880*Dmsqst2*(pow2(Mgl)*(pow2(
        Mst1)*pow2(Mst2) + pow4(Mst1) + pow4(Mst2)) - Dmglst2*(Dmglst2 + Mgl)*(
        7*pow2(Mst1)*pow2(Mst2) + 11*pow4(Mst1) + 3*pow4(Mst2))) - pow2(Mst2)*(
        pow2(Dmglst2)*(-374379*pow2(Mst1)*pow2(Mst2) + 1307164*pow4(Mst1) -
        290664*pow4(Mst2)) + Dmglst2*Mgl*(361557*pow2(Mst1)*pow2(Mst2) +
        1307164*pow4(Mst1) + 243*pow4(Mst2)) - 3*pow2(Mgl)*(72681*pow2(Mst1)*
        pow2(Mst2) + 115516*pow4(Mst1) + 25947*pow4(Mst2)))) + pow2(s2t)*(2430*
        pow2(Dmsqst2)*pow2(Mgl)*(-4*pow2(Mst2)*(pow2(Mst1) + pow2(Mst2)) + 15*
        pow4(Mst1)) + pow4(Mst2)*(19440*Dmsqst2*(Dmglst2*(Dmglst2 + Mgl)*(3*
        pow2(Mst1) + pow2(Mst2)) - pow2(Mgl)*(pow2(Mst1) + 3*pow2(Mst2))) + 3*
        pow2(Mgl)*(19554*pow2(Mst1)*pow2(Mst2) + 22315*pow4(Mst1) - 53685*pow4(
        Mst2)) - Dmglst2*Mgl*(36810*pow2(Mst1)*pow2(Mst2) + 282607*pow4(Mst1) +
        408699*pow4(Mst2)) - pow2(Dmglst2)*(-171510*pow2(Mst1)*pow2(Mst2) +
        282607*pow4(Mst1) + 998613*pow4(Mst2)))))))/(pow2(Mgl)*pow6(Mst2))) -
        2058000*z4*(-((540*pow2(Dmsqst2) + (60*Dmsqst2*pow2(Mst1)*(24*Dmsqst2 -
        7*(2*pow2(Mst1) + 9*pow2(Mst2))))/pow2(Mst2) - pow4(Mst1)*(3441 - (
        16796*Dmglst2*(Dmglst2 + Mgl))/pow2(Mgl) - (840*pow2(Dmsqst2))/pow4(
        Mst2)) + (-6*(270*Dmsqst2*pow2(Mgl) + (486*Dmglst2*Mgl + 2513*pow2(
        Dmglst2) + 2334*pow2(Mgl))*pow2(Mst1))*pow2(Mst2) + 36*(318*Dmglst2*Mgl
        + 1111*pow2(Dmglst2) - 90*pow2(Mgl))*pow4(Mst2))/pow2(Mgl))*pow4(s2t))
        + (4*Mt*(4*pow2(Mst2)*pow2(Mt)*(pow2(Dmglst2)*(3*(-9406*Mt + 4307*Mst2*
        s2t)*pow2(Mst1)*pow2(Mst2) + 4*(5582*Mt - 12865*Mst2*s2t)*pow4(Mst1) +
        117*(-14*Mt + 67*Mst2*s2t)*pow4(Mst2)) + 3*pow2(Mgl)*(3*(209*Mt + 377*
        Mst2*s2t)*pow2(Mst1)*pow2(Mst2) + 4*(292*Mt + 685*Mst2*s2t)*pow4(Mst1)
        + 27*(7*Mt + 103*Mst2*s2t)*pow4(Mst2)) + Dmglst2*Mgl*(-9*(-644*Mt +
        1383*Mst2*s2t)*pow2(Mst1)*pow2(Mst2) + 4*(5582*Mt - 12865*Mst2*s2t)*
        pow4(Mst1) + 27*(28*Mt + 253*Mst2*s2t)*pow4(Mst2))) + 3*Mt*pow2(s2t)*(
        20*Dmsqst2*pow2(Mgl)*(3*Dmsqst2*(5*pow2(Mst1)*pow2(Mst2) + 5*pow4(Mst1)
        + 3*pow4(Mst2)) - pow2(Mst2)*(18*pow2(Mst1)*pow2(Mst2) + 8*pow4(Mst1) +
        27*pow4(Mst2))) + pow4(Mst2)*(12*pow2(Dmglst2)*(597*pow2(Mst1)*pow2(
        Mst2) + 1538*pow4(Mst1) - 423*pow4(Mst2)) + 24*Dmglst2*Mgl*(150*pow2(
        Mst1)*pow2(Mst2) + 769*pow4(Mst1) - 162*pow4(Mst2)) - pow2(Mgl)*(2478*
        pow2(Mst1)*pow2(Mst2) + 6857*pow4(Mst1) + 2727*pow4(Mst2)))) + pow3(
        s2t)*(pow2(Dmglst2)*(-14550*pow2(Mst1)*pow2(Mst2) + 15571*pow4(Mst1) -
        11241*pow4(Mst2)) + Dmglst2*Mgl*(1872*pow2(Mst1)*pow2(Mst2) + 15571*
        pow4(Mst1) - 7317*pow4(Mst2)) - 3*pow2(Mgl)*(2220*pow2(Mst1)*pow2(Mst2)
        + 439*pow4(Mst1) + 2295*pow4(Mst2)))*pow5(Mst2)))/(pow2(Mgl)*pow6(Mst2)
        )) + (1029000*(4*Dmglst2*Mgl*pow2(Mst2)*(-((4*T1ep + 3*pow2(z2))*pow4(
        Mst1)*(55368*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 205840*Mst2*s2t*pow3(Mt) +
        15571*Mt*pow3(Mst2)*pow3(s2t) + 89312*pow4(Mt) - 4199*pow4(Mst2)*pow4(
        s2t))) + 18*pow2(Mst1)*pow2(Mst2)*(-600*pow2(Mst2)*pow2(Mt)*pow2(s2t)*(
        4*T1ep + 3*pow2(z2)) + 2766*Mst2*s2t*(4*T1ep + 3*pow2(z2))*pow3(Mt) -
        401*Mt*(4*T1ep + 3*pow2(z2))*pow3(Mst2)*pow3(s2t) - 1288*(4*T1ep + 3*
        pow2(z2))*pow4(Mt) + 432*T1ep*pow4(Mst2)*pow4(s2t)) + 27*pow4(Mst2)*(
        140*Mst2*s2t*(4*T1ep + 3*pow2(z2))*pow3(Mt) - Mt*(140*T1ep + 5289*pow2(
        z2))*pow3(Mst2)*pow3(s2t) - 112*(4*T1ep + 3*pow2(z2))*pow4(Mt) + (28*
        T1ep + 237*pow2(z2))*pow4(Mst2)*pow4(s2t))) + 2*pow2(Dmglst2)*pow2(
        Mst2)*(-2*(4*T1ep + 3*pow2(z2))*pow4(Mst1)*(55368*pow2(Mst2)*pow2(Mt)*
        pow2(s2t) - 205840*Mst2*s2t*pow3(Mt) + 15571*Mt*pow3(Mst2)*pow3(s2t) +
        89312*pow4(Mt) - 4199*pow4(Mst2)*pow4(s2t)) + 24*pow2(Mst1)*pow2(Mst2)*
        (-1791*pow2(Mst2)*pow2(Mt)*pow2(s2t)*(4*T1ep + 3*pow2(z2)) - 4307*Mst2*
        s2t*(4*T1ep + 3*pow2(z2))*pow3(Mt) + 767*Mt*(4*T1ep + 3*pow2(z2))*pow3(
        Mst2)*pow3(s2t) + 9406*(4*T1ep + 3*pow2(z2))*pow4(Mt) + (80*T1ep - 669*
        pow2(z2))*pow4(Mst2)*pow4(s2t)) + 9*pow4(Mst2)*(-504*pow2(Mst2)*pow2(
        Mt)*pow2(s2t)*(4*T1ep + 3*pow2(z2)) - 56*Mst2*s2t*(4*T1ep + 3*pow2(z2))
        *pow3(Mt) + 2*Mt*(28*T1ep - 31083*pow2(z2))*pow3(Mst2)*pow3(s2t) +
        1456*(4*T1ep + 3*pow2(z2))*pow4(Mt) + (140*T1ep + 2049*pow2(z2))*pow4(
        Mst2)*pow4(s2t))) - 3*pow2(Mgl)*((4*T1ep + 3*pow2(z2))*pow4(Mst1)*(
        1200*pow2(Dmsqst2)*pow2(Mt)*pow2(s2t) + 43840*s2t*pow3(Mst2)*pow3(Mt) -
        4*pow4(Mst2)*(6857*pow2(Mt)*pow2(s2t) - 70*Dmsqst2*pow4(s2t)) + 8*pow2(
        Mst2)*(-80*Dmsqst2*pow2(Mt)*pow2(s2t) + 2336*pow4(Mt) - 35*pow2(
        Dmsqst2)*pow4(s2t)) - 1756*Mt*pow3(s2t)*pow5(Mst2) + 1795*pow4(s2t)*
        pow6(Mst2)) + 18*pow4(Mst2)*(40*pow2(Dmsqst2)*pow2(Mt)*pow2(s2t)*(4*
        T1ep + 3*pow2(z2)) + 168*s2t*(4*T1ep + 3*pow2(z2))*pow3(Mst2)*pow3(Mt)
        + 6*pow2(s2t)*(-53*pow2(Mt) + 5*Dmsqst2*pow2(s2t))*(4*T1ep + 3*pow2(z2)
        )*pow4(Mst2) + 2*pow2(Mst2)*(4*T1ep + 3*pow2(z2))*(-60*Dmsqst2*pow2(Mt)
        *pow2(s2t) + 84*pow4(Mt) - 5*pow2(Dmsqst2)*pow4(s2t)) - 6*Mt*(28*T1ep -
        555*pow2(z2))*pow3(s2t)*pow5(Mst2) + 3*(92*T1ep - 3*pow2(z2))*pow4(s2t)
        *pow6(Mst2)) + 6*pow2(Mst1)*pow2(Mst2)*(200*pow2(Dmsqst2)*pow2(Mt)*
        pow2(s2t)*(4*T1ep + 3*pow2(z2)) + 3016*s2t*(4*T1ep + 3*pow2(z2))*pow3(
        Mst2)*pow3(Mt) + 2*pow2(s2t)*(-1258*pow2(Mt) + 105*Dmsqst2*pow2(s2t))*(
        4*T1ep + 3*pow2(z2))*pow4(Mst2) + 8*pow2(Mst2)*(4*T1ep + 3*pow2(z2))*(-
        30*Dmsqst2*pow2(Mt)*pow2(s2t) + 209*pow4(Mt) - 10*pow2(Dmsqst2)*pow4(
        s2t)) - 292*Mt*(4*T1ep + 3*pow2(z2))*pow3(s2t)*pow5(Mst2) + (2572*T1ep
        + 2577*pow2(z2))*pow4(s2t)*pow6(Mst2)))))/(pow2(Mgl)*pow6(Mst2)) + (
        68600*Mt*pow3(s2t)*(72900*pow2(Dmsqst2)*pow2(Mgl)*pow2(Mst1)*(16*(-2 +
        z2)*pow2(Mst1)*pow2(Mst2) + (-69 + 90*z2)*pow4(Mst1) - 48*(-1 + z2)*
        pow4(Mst2)) - 12960*pow2(Mst1)*pow3(log(pow2(Mst2)/pow2(Mst1)))*pow4(
        Mst2)*(5*pow2(Mgl)*(39*pow2(Mst1)*pow2(Mst2) + 40*pow4(Mst1) + 39*pow4(
        Mst2)) + 3*Dmglst2*Mgl*(65*pow2(Mst1)*pow2(Mst2) + 56*pow4(Mst1) + 129*
        pow4(Mst2)) + 3*pow2(Dmglst2)*(65*pow2(Mst1)*pow2(Mst2) + 56*pow4(Mst1)
        + 225*pow4(Mst2))) - 1166400*Dmsqst2*pow2(Mst1)*pow2(Mst2)*(Dmglst2*
        Mgl*(2*(-12 + 5*z2)*pow2(Mst1)*pow2(Mst2) + (-20 + 11*z2)*pow4(Mst1) -
        4*(-1 + z2)*pow4(Mst2)) + pow2(Dmglst2)*(2*(-12 + 5*z2)*pow2(Mst1)*
        pow2(Mst2) + (-20 + 11*z2)*pow4(Mst1) - 4*(-1 + z2)*pow4(Mst2)) + pow2(
        Mgl)*(2*(1 + z2)*pow2(Mst1)*pow2(Mst2) - (-2 + z2)*pow4(Mst1) + 2*(-1 +
        2*z2)*pow4(Mst2))) + 25920*pow2(log(pow2(Mst2)/pow2(Mst1)))*pow4(Mst2)*
        (-3*pow2(Dmglst2)*(67*pow2(Mst2)*pow4(Mst1) - 15*pow2(Mst1)*pow4(Mst2)
        + 41*pow6(Mst1) - 48*pow6(Mst2)) - 3*Dmglst2*Mgl*(15*pow2(Mst2)*pow4(
        Mst1) - 92*pow2(Mst1)*pow4(Mst2) + 41*pow6(Mst1) - 24*pow6(Mst2)) +
        pow2(Mgl)*(183*pow2(Mst2)*pow4(Mst1) + 354*pow2(Mst1)*pow4(Mst2) + 125*
        pow6(Mst1) + 24*pow6(Mst2))) + pow4(Mst2)*(pow2(Dmglst2)*(60*(280885 +
        11016*B4 + 324*DN - 24544*OepS2 - 89856*S2 + 858030*z2)*pow2(Mst2)*
        pow4(Mst1) + 9*(-2048393 + 1058400*B4 - 23760*DN - 1120*OepS2 - 449388*
        S2 - 5386422*z2)*pow2(Mst1)*pow4(Mst2) + (-31897243 + 2491360*OepS2 -
        90290268*S2 - 51499698*z2)*pow6(Mst1) - 2488320*pow6(Mst2)) + Dmglst2*
        Mgl*(18*(27211 + 36720*B4 + 1080*DN + 64160*OepS2 - 1515564*S2 -
        1240446*z2)*pow2(Mst2)*pow4(Mst1) + 27*(69997 + 188640*B4 - 3600*DN +
        5600*OepS2 + 146772*S2 - 1043682*z2)*pow2(Mst1)*pow4(Mst2) + (-31897243
        + 2491360*OepS2 - 90290268*S2 - 51499698*z2)*pow6(Mst1) - 622080*pow6(
        Mst2)) + 3*pow2(Mgl)*(30*(51041 + 7344*B4 + 216*DN - 2336*OepS2 +
        112428*S2 + 28038*z2)*pow2(Mst2)*pow4(Mst1) + 135*(25611 + 5280*B4 -
        48*DN - 224*OepS2 - 324*S2 - 13614*z2)*pow2(Mst1)*pow4(Mst2) + (2759935
        - 70240*OepS2 + 3305340*S2 + 629034*z2)*pow6(Mst1) + 207360*pow6(Mst2))
        ) + 180*log(pow2(Mst2)/pow2(Mst1))*(405*pow2(Dmsqst2)*pow2(Mgl)*pow2(
        Mst1)*(-16*(-2 + z2)*pow2(Mst1)*pow2(Mst2) + (-121 + 60*z2)*pow4(Mst1)
        - 16*(-2 + z2)*pow4(Mst2)) - 3240*Dmsqst2*pow2(Mst1)*pow2(Mst2)*(
        Dmglst2*Mgl*(-4*(-5 + 3*z2)*pow2(Mst1)*pow2(Mst2) + 3*pow4(Mst1) + 4*(5
        - 3*z2)*pow4(Mst2)) + pow2(Dmglst2)*(-4*(-5 + 3*z2)*pow2(Mst1)*pow2(
        Mst2) + 3*pow4(Mst1) + 4*(5 - 3*z2)*pow4(Mst2)) - pow2(Mgl)*(-4*(-1 +
        z2)*pow2(Mst1)*pow2(Mst2) + pow4(Mst1) - 4*(-1 + z2)*pow4(Mst2))) + 2*
        pow4(Mst2)*(9*pow2(Mgl)*(6*(221 + 219*S2 - 234*z2)*pow2(Mst2)*pow4(
        Mst1) + 21*(-208 + 27*S2 - 2*z2)*pow2(Mst1)*pow4(Mst2) + (46 + 1317*S2
        - 918*z2)*pow6(Mst1) - 384*pow6(Mst2)) - Dmglst2*Mgl*(18*(-503 + 3609*
        S2 - 1494*z2)*pow2(Mst2)*pow4(Mst1) + 27*(1640 + 315*S2 - 1146*z2)*
        pow2(Mst1)*pow4(Mst2) + (18652 + 140139*S2 - 38826*z2)*pow6(Mst1) +
        3456*pow6(Mst2)) + pow2(Dmglst2)*(6*(-3449 + 13806*S2 + 4392*z2)*pow2(
        Mst2)*pow4(Mst1) + 9*(-3164 + 63*S2 + 3978*z2)*pow2(Mst1)*pow4(Mst2) +
        (-18652 - 140139*S2 + 38826*z2)*pow6(Mst1) + 3456*pow6(Mst2))))))/(
        pow2(Mgl)*pow2(Mst1)*pow5(Mst2)) - (39200*s2t*pow3(Mt)*(-1020600*pow2(
        Dmsqst2)*pow2(Mgl)*pow2(Mst1)*(44*pow2(Mst1) + 15*pow2(Mst2)) + 90720*
        pow2(Mst1)*pow2(Mst2)*pow3(log(pow2(Mst2)/pow2(Mst1)))*(3*pow2(Dmglst2)
        *(-16*pow2(Mst1)*pow2(Mst2) + 42*pow4(Mst1) - 49*pow4(Mst2)) + 3*
        Dmglst2*Mgl*(16*pow2(Mst1)*pow2(Mst2) + 42*pow4(Mst1) - 49*pow4(Mst2))
        - pow2(Mgl)*(144*pow2(Mst1)*pow2(Mst2) + 122*pow4(Mst1) + 147*pow4(
        Mst2))) + 226800*Dmsqst2*pow2(Mst1)*(-9*pow2(Mgl)*(24*pow2(Mst1)*pow2(
        Mst2) + 35*pow4(Mst1) - 14*pow4(Mst2)) + Dmglst2*Mgl*(2016*pow2(Mst1)*
        pow2(Mst2) + 4185*pow4(Mst1) + 578*pow4(Mst2)) + pow2(Dmglst2)*(2016*
        pow2(Mst1)*pow2(Mst2) + 4185*pow4(Mst1) + 578*pow4(Mst2))) - pow2(Mst2)
        *(-7*Dmglst2*Mgl*(27*(-322421 + 73760*OepS2 - 564876*S2)*pow2(Mst2)*
        pow4(Mst1) + 135*(-30183 + 1120*OepS2 + 43740*S2)*pow2(Mst1)*pow4(Mst2)
        + 4*(-20964397 + 2058400*OepS2 - 47545272*S2)*pow6(Mst1) - 622080*pow6(
        Mst2)) + 21*pow2(Mgl)*(3*(-1156193 + 60320*OepS2 - 1414908*S2)*pow2(
        Mst2)*pow4(Mst1) + 27*(-60759 + 1120*OepS2 + 25596*S2)*pow2(Mst1)*pow4(
        Mst2) + 4*(-2016907 + 109600*OepS2 - 3520152*S2)*pow6(Mst1) - 207360*
        pow6(Mst2)) + pow2(Dmglst2)*(3*(79386499 + 4823840*OepS2 + 82155924*S2)
        *pow2(Mst2)*pow4(Mst1) + 9*(-6779161 + 7840*OepS2 + 3032964*S2)*pow2(
        Mst1)*pow4(Mst2) - 28*(-20964397 + 2058400*OepS2 - 47545272*S2)*pow6(
        Mst1) + 17418240*pow6(Mst2))) + 2520*log(pow2(Mst2)/pow2(Mst1))*(2430*
        pow2(Dmsqst2)*pow2(Mgl)*pow2(Mst1)*(8*pow2(Mst1) + 5*pow2(Mst2)) - 540*
        Dmsqst2*pow2(Mst1)*(-3*pow2(Mgl)*(16*pow2(Mst1)*pow2(Mst2) + 25*pow4(
        Mst1) + 10*pow4(Mst2)) + Dmglst2*Mgl*(288*pow2(Mst1)*pow2(Mst2) + 561*
        pow4(Mst1) + 118*pow4(Mst2)) + pow2(Dmglst2)*(288*pow2(Mst1)*pow2(Mst2)
        + 561*pow4(Mst1) + 118*pow4(Mst2))) + pow2(Mst2)*(Dmglst2*Mgl*(3*(11408
        - 37341*S2)*pow2(Mst2)*pow4(Mst1) + 3*(8398 - 2835*S2)*pow2(Mst1)*pow4(
        Mst2) + (18896 - 463140*S2)*pow6(Mst1) - 3456*pow6(Mst2)) + 9*pow2(Mgl)
        *(9*(42 + 377*S2)*pow2(Mst2)*pow4(Mst1) + 3*(26 + 189*S2)*pow2(Mst1)*
        pow4(Mst2) + (668 + 8220*S2)*pow6(Mst1) - 384*pow6(Mst2)) + pow2(
        Dmglst2)*(1971*(35 + 59*S2)*pow2(Mst2)*pow4(Mst1) + 9*(8618 + 63*S2)*
        pow2(Mst1)*pow4(Mst2) + (18896 - 463140*S2)*pow6(Mst1) + 3456*pow6(
        Mst2)))) + 45360*pow2(log(pow2(Mst2)/pow2(Mst1)))*(2*pow2(Mgl)*(612*
        pow4(Mst1)*pow4(Mst2) - 90*Dmsqst2*pow6(Mst1) + 565*pow2(Mst2)*pow6(
        Mst1) + 465*pow2(Mst1)*pow6(Mst2) + 48*pow8(Mst2)) + 2*Dmglst2*Mgl*(-
        1097*pow4(Mst1)*pow4(Mst2) + 270*Dmsqst2*pow6(Mst1) - 1493*pow2(Mst2)*
        pow6(Mst1) - 199*pow2(Mst1)*pow6(Mst2) + 144*pow8(Mst2)) + pow2(
        Dmglst2)*(295*pow4(Mst1)*pow4(Mst2) + 540*Dmsqst2*pow6(Mst1) - 2986*
        pow2(Mst2)*pow6(Mst1) - 1170*pow2(Mst1)*pow6(Mst2) + 576*pow8(Mst2)))))
        /(pow2(Mgl)*pow2(Mst1)*pow5(Mst2)) + (pow4(s2t)*(-55566000*pow3(log(
        pow2(Mst2)/pow2(Mst1)))*(4*Dmglst2*Mgl*(58*pow2(Mst1)*pow2(Mst2) + 33*
        pow4(Mst1) + 321*pow4(Mst2)) + pow2(Mgl)*(460*pow2(Mst1)*pow2(Mst2) +
        409*pow4(Mst1) + 813*pow4(Mst2)) + 2*pow2(Dmglst2)*(194*pow2(Mst1)*
        pow2(Mst2) + 66*pow4(Mst1) + 963*pow4(Mst2))) + (-60*pow2(Dmsqst2)*
        pow2(Mgl)*pow2(Mst1)*(16464*(-120796 + 4000*OepS2 - 205875*S2 - 45375*
        z2)*pow2(Mst2)*pow4(Mst1) + 385875*(1061 + 64*OepS2 - 27864*S2 - 2076*
        z2)*pow2(Mst1)*pow4(Mst2) + 5*(38390869 + 7683200*OepS2 - 418597200*S2
        - 103826100*z2)*pow6(Mst1) - 83349000*(-2 + z2)*pow6(Mst2)) + 120*
        Dmsqst2*pow2(Mst1)*pow2(Mst2)*(4630500*Dmglst2*Mgl*(9*(-75 + 32*z2)*
        pow2(Mst2)*pow4(Mst1) - 72*(-1 + 2*z2)*pow2(Mst1)*pow4(Mst2) + 2*(38 +
        63*z2)*pow6(Mst1) - 36*(-3 + z2)*pow6(Mst2)) + 4630500*pow2(Dmglst2)*(
        9*(-75 + 32*z2)*pow2(Mst2)*pow4(Mst1) - 72*(-1 + 2*z2)*pow2(Mst1)*pow4(
        Mst2) + 2*(38 + 63*z2)*pow6(Mst1) - 36*(-3 + z2)*pow6(Mst2)) + pow2(
        Mgl)*(3087*(52199 + 28000*OepS2 - 3415500*S2 - 557250*z2)*pow2(Mst2)*
        pow4(Mst1) + 1157625*(-1213 + 32*OepS2 + 4860*S2 - 462*z2)*pow2(Mst1)*
        pow4(Mst2) + 2*(378483467 + 9604000*OepS2 - 754771500*S2 - 306899250*
        z2)*pow6(Mst1) + 41674500*(1 + 2*z2)*pow6(Mst2))) + pow4(Mst2)*(-4*
        Dmglst2*Mgl*pow2(Mst1)*(-12348*(-95041 + 81000*B4 - 108000*D3 + 67500*
        DN - 432000*OepS2 - 7425000*S2 + 5325300*z2)*pow2(Mst2)*pow4(Mst1) +
        2315250*(-15707 + 432*B4 - 576*D3 + 360*DN + 224*OepS2 + 543348*S2 -
        11802*z2)*pow2(Mst1)*pow4(Mst2) + (-119405394763 + 11522056000*OepS2 -
        478191735000*S2 - 233890773900*z2)*pow6(Mst1) + 333396000*(5 + 7*z2)*
        pow6(Mst2)) + 3*pow2(Mgl)*pow2(Mst1)*(20580*(-7680127 + 10800*B4 -
        21600*D3 + 16200*DN + 514400*OepS2 - 46307700*S2 - 10868970*z2)*pow2(
        Mst2)*pow4(Mst1) - 1543500*(25289 + 1440*B4 - 144*D3 + 72*DN - 2208*
        OepS2 + 298404*S2 + 36318*z2)*pow2(Mst1)*pow4(Mst2) + (93508520089 +
        2000376000*B4 + 222264000*D3 - 222264000*DN + 4925480000*OepS2 -
        312095700000*S2 - 77106571500*z2)*pow6(Mst1) + 55566000*(135 + 38*z2)*
        pow6(Mst2)) - 4*pow2(Dmglst2)*(77175*(-1395899 + 54000*B4 - 60480*D3 +
        33480*DN + 5600*OepS2 + 37625796*S2 - 481746*z2)*pow4(Mst1)*pow4(Mst2)
        - 2058*(-22023067 + 729000*B4 - 972000*D3 + 607500*DN - 320000*OepS2 -
        120274200*S2 - 23129700*z2)*pow2(Mst2)*pow6(Mst1) - 166698000*(-245 +
        12*z2)*pow2(Mst1)*pow6(Mst2) + (-119405394763 + 11522056000*OepS2 -
        478191735000*S2 - 233890773900*z2)*pow8(Mst1) - 2667168000*pow8(Mst2)))
        )/(pow4(Mst1)*pow4(Mst2)) + (1260*log(pow2(Mst2)/pow2(Mst1))*(30*pow2(
        Dmsqst2)*pow2(Mgl)*pow4(Mst1)*(294*(-1837 + 7200*S2 + 900*z2)*pow2(
        Mst1)*pow2(Mst2) + 5*(8753 + 246960*S2 - 79380*z2)*pow4(Mst1) + 44100*(
        -1 + 18*S2 + 3*z2)*pow4(Mst2)) + 120*Dmsqst2*pow2(Mst1)*pow2(Mst2)*(
        11025*Dmglst2*Mgl*pow2(Mst1)*(6*(13 - 4*z2)*pow2(Mst1)*pow2(Mst2) + (-
        83 + 60*z2)*pow4(Mst1) + 12*(2 - 3*z2)*pow4(Mst2)) + 11025*pow2(
        Dmglst2)*pow2(Mst1)*(6*(13 - 4*z2)*pow2(Mst1)*pow2(Mst2) + (-83 + 60*
        z2)*pow4(Mst1) + 12*(2 - 3*z2)*pow4(Mst2)) + pow2(Mgl)*(-441*(-317 +
        3150*S2 - 300*z2)*pow2(Mst2)*pow4(Mst1) - 22050*(17 + 27*S2)*pow2(Mst1)
        *pow4(Mst2) + (621671 - 308700*S2 - 132300*z2)*pow6(Mst1) + 66150*pow6(
        Mst2))) - pow4(Mst2)*(-4*Dmglst2*Mgl*pow2(Mst1)*(294*(-81643 + 291600*
        S2 + 5850*z2)*pow2(Mst2)*pow4(Mst1) + 132300*(277 + 63*S2 - 140*z2)*
        pow2(Mst1)*pow4(Mst2) + (-8287903 + 185175900*S2 - 9525600*z2)*pow6(
        Mst1) - 17331300*pow6(Mst2)) + 3*pow2(Mgl)*pow2(Mst1)*(1960*(-16051 +
        86805*S2 - 8325*z2)*pow2(Mst2)*pow4(Mst1) + 29400*(1148 + 1863*S2 +
        276*z2)*pow2(Mst1)*pow4(Mst2) + 3*(-23468297 + 26386500*S2 + 661500*z2)
        *pow6(Mst1) + 11025000*pow6(Mst2)) - 4*pow2(Dmglst2)*(22050*(3050 +
        315*S2 - 222*z2)*pow4(Mst1)*pow4(Mst2) + 441*(-165199 + 24000*S2 +
        20250*z2)*pow2(Mst2)*pow6(Mst1) - 13957650*pow2(Mst1)*pow6(Mst2) + (-
        8287903 + 185175900*S2 - 9525600*z2)*pow8(Mst1) - 1058400*pow8(Mst2))))
        )/(pow4(Mst1)*pow4(Mst2)) + (1587600*pow2(log(pow2(Mst2)/pow2(Mst1)))*(
        90*Dmsqst2*pow2(Mgl)*pow2(Mst2)*pow4(Mst1)*(217*pow2(Mst1)*pow2(Mst2) +
        106*pow4(Mst1) - 140*pow4(Mst2)) - 270*pow2(Dmsqst2)*pow2(Mgl)*(25*
        pow2(Mst1) + 14*pow2(Mst2))*pow6(Mst1) + pow4(Mst2)*(pow2(Mgl)*pow2(
        Mst1)*(12880*pow2(Mst2)*pow4(Mst1) + 77910*pow2(Mst1)*pow4(Mst2) +
        47056*pow6(Mst1) + 12915*pow6(Mst2)) + 2*Dmglst2*Mgl*pow2(Mst1)*(-
        34783*pow2(Mst2)*pow4(Mst1) + 29820*pow2(Mst1)*pow4(Mst2) - 1377*pow6(
        Mst1) + 19635*pow6(Mst2)) - pow2(Dmglst2)*(2730*pow4(Mst1)*pow4(Mst2) +
        99617*pow2(Mst2)*pow6(Mst1) - 85785*pow2(Mst1)*pow6(Mst2) + 2754*pow8(
        Mst1) + 10080*pow8(Mst2)))))/(pow4(Mst1)*pow4(Mst2))))/pow2(Mgl) + (16*
        pow4(Mt)*(-185220*pow2(Dmsqst2)*pow2(Mgl)*pow2(Mst1)*(386875*pow2(Mst1)
        *pow2(Mst2) + 1058838*pow4(Mst1) + 54000*pow4(Mst2)) + 111132000*pow2(
        Mst2)*pow3(log(pow2(Mst2)/pow2(Mst1)))*pow4(Mst1)*(pow2(Mgl)*(55*pow2(
        Mst1)*pow2(Mst2) + 200*pow4(Mst1) - 3*pow4(Mst2)) + 2*Dmglst2*Mgl*(329*
        pow2(Mst1)*pow2(Mst2) + 464*pow4(Mst1) + 96*pow4(Mst2)) + pow2(Dmglst2)
        *(55*pow2(Mst1)*pow2(Mst2) + 928*pow4(Mst1) + 288*pow4(Mst2))) - 2160*
        Dmsqst2*pow2(Mst1)*(-85750*Dmglst2*Mgl*pow2(Mst2)*(2954*pow2(Mst1)*
        pow2(Mst2) + 10827*pow4(Mst1) + 324*pow4(Mst2)) - 85750*pow2(Dmglst2)*
        pow2(Mst2)*(2954*pow2(Mst1)*pow2(Mst2) + 10827*pow4(Mst1) + 324*pow4(
        Mst2)) + 3*pow2(Mgl)*(31323103*pow2(Mst2)*pow4(Mst1) + 8832250*pow2(
        Mst1)*pow4(Mst2) + 52662272*pow6(Mst1) - 771750*pow6(Mst2))) + pow2(
        Mst2)*(4*Dmglst2*Mgl*pow2(Mst1)*(882*(98884301 + 4508000*OepS2 -
        17853750*S2)*pow2(Mst2)*pow4(Mst1) + 132300*(773533 + 3920*OepS2 +
        129033*S2)*pow2(Mst1)*pow4(Mst2) + 5*(-54946675289 + 3063401600*OepS2 -
        52952863320*S2)*pow6(Mst1) - 7001316000*pow6(Mst2)) + 3*pow2(Mgl)*(
        154350*(389597 + 3360*OepS2 - 105948*S2)*pow4(Mst1)*pow4(Mst2) + 10290*
        (7319011 + 167200*OepS2 - 6744060*S2)*pow2(Mst2)*pow6(Mst1) +
        3945186000*pow2(Mst1)*pow6(Mst2) + (187545955063 + 3204992000*OepS2 -
        128216692800*S2)*pow8(Mst1) - 889056000*pow8(Mst2)) + 4*pow2(Dmglst2)*(
        -15435*(-13292051 + 72800*OepS2 + 3133890*S2)*pow4(Mst1)*pow4(Mst2) -
        147*(55610713 + 131684000*OepS2 + 173875950*S2)*pow2(Mst2)*pow6(Mst1) -
        32839506000*pow2(Mst1)*pow6(Mst2) + 5*(-54946675289 + 3063401600*OepS2
        - 52952863320*S2)*pow8(Mst1) + 2667168000*pow8(Mst2))) + 793800*pow2(
        log(pow2(Mst2)/pow2(Mst1)))*(630*pow2(Dmsqst2)*pow2(Mgl)*(8*pow2(Mst1)
        - 5*pow2(Mst2))*pow4(Mst1) - 180*Dmsqst2*(-140*Dmglst2*Mgl*pow2(Mst2)*
        pow6(Mst1) - 140*pow2(Dmglst2)*pow2(Mst2)*pow6(Mst1) + pow2(Mgl)*(-35*
        pow4(Mst1)*pow4(Mst2) + 63*pow2(Mst2)*pow6(Mst1) + 342*pow8(Mst1))) +
        pow2(Mst2)*(pow2(Mgl)*(67760*pow4(Mst1)*pow4(Mst2) + 79520*pow2(Mst2)*
        pow6(Mst1) + 12390*pow2(Mst1)*pow6(Mst2) + 34339*pow8(Mst1) - 3360*
        pow8(Mst2)) - 4*Dmglst2*Mgl*(33320*pow4(Mst1)*pow4(Mst2) + 125013*pow2(
        Mst2)*pow6(Mst1) - 6195*pow2(Mst1)*pow6(Mst2) + 203320*pow8(Mst1) +
        3360*pow8(Mst2)) - 2*pow2(Dmglst2)*(8708*pow4(Mst1)*pow4(Mst2) -
        304339*pow2(Mst2)*pow6(Mst1) - 18585*pow2(Mst1)*pow6(Mst2) + 406640*
        pow8(Mst1) + 16800*pow8(Mst2)))) + 2520*log(pow2(Mst2)/pow2(Mst1))*(
        2205*pow2(Dmsqst2)*pow2(Mgl)*(23226*pow2(Mst1) + 9125*pow2(Mst2))*pow4(
        Mst1) + 180*Dmsqst2*(-1225*Dmglst2*Mgl*pow2(Mst2)*(2466*pow2(Mst1) +
        571*pow2(Mst2))*pow4(Mst1) - 1225*pow2(Dmglst2)*pow2(Mst2)*(2466*pow2(
        Mst1) + 571*pow2(Mst2))*pow4(Mst1) + 3*pow2(Mgl)*(56350*pow4(Mst1)*
        pow4(Mst2) + 150283*pow2(Mst2)*pow6(Mst1) + 7350*pow2(Mst1)*pow6(Mst2)
        + 350242*pow8(Mst1))) + pow2(Mst2)*(-2*Dmglst2*Mgl*(3675*(-6569 + 2268*
        S2)*pow4(Mst1)*pow4(Mst2) + 882*(8581 + 72450*S2)*pow2(Mst2)*pow6(Mst1)
        + 4630500*pow2(Mst1)*pow6(Mst2) + 5*(42300121 + 49233240*S2)*pow8(Mst1)
        - 2116800*pow8(Mst2)) + 3*pow2(Mgl)*(-17150*(-79 + 243*S2)*pow4(Mst1)*
        pow4(Mst2) - 5390*(-2032 + 2565*S2)*pow2(Mst2)*pow6(Mst1) - 2690100*
        pow2(Mst1)*pow6(Mst2) + (35585111 - 25754400*S2)*pow8(Mst1) + 705600*
        pow8(Mst2)) + pow2(Dmglst2)*(294*(50551 + 122850*S2)*pow4(Mst1)*pow4(
        Mst2) + 294*(-297379 + 2116350*S2)*pow2(Mst2)*pow6(Mst1) + 1719900*
        pow2(Mst1)*pow6(Mst2) - 10*(42300121 + 49233240*S2)*pow8(Mst1) +
        2116800*pow8(Mst2))))))/(pow2(Mgl)*pow4(Mst1)*pow6(Mst2)) + (24*pow2(
        Mt)*(-3333960000*pow2(Mst2)*pow3(Dmsqst2) + (pow2(s2t)*(-18522000*pow3(
        log(pow2(Mst2)/pow2(Mst1)))*pow4(Mst1)*pow4(Mst2)*(pow2(Mgl)*(583*pow2(
        Mst1)*pow2(Mst2) + 337*pow4(Mst1) + 21*pow4(Mst2)) + 4*Dmglst2*Mgl*(
        289*pow2(Mst1)*pow2(Mst2) + 362*pow4(Mst1) + 96*pow4(Mst2)) + 2*pow2(
        Dmglst2)*(863*pow2(Mst1)*pow2(Mst2) + 724*pow4(Mst1) + 288*pow4(Mst2)))
        + 20*Dmsqst2*pow2(Mst1)*pow2(Mst2)*(pow2(Mgl)*(-6174*(36449 + 4000*
        OepS2 - 310500*S2)*pow2(Mst2)*pow4(Mst1) - 1157625*(-1173 + 32*OepS2 -
        2916*S2)*pow2(Mst1)*pow4(Mst2) + (44860502 - 10976000*OepS2 +
        648270000*S2)*pow6(Mst1) - 83349000*pow6(Mst2)) + 2315250*Dmglst2*Mgl*(
        -1356*pow2(Mst2)*pow4(Mst1) - 552*pow2(Mst1)*pow4(Mst2) + 1931*pow6(
        Mst1) - 432*pow6(Mst2)) + 2315250*pow2(Dmglst2)*(-1356*pow2(Mst2)*pow4(
        Mst1) - 552*pow2(Mst1)*pow4(Mst2) + 1931*pow6(Mst1) - 432*pow6(Mst2)))
        + 30*pow2(Dmsqst2)*pow2(Mgl)*pow2(Mst1)*(1372*(121757 + 10000*OepS2 -
        479250*S2)*pow2(Mst2)*pow4(Mst1) + 1029000*(-20 + 8*OepS2 - 81*S2)*
        pow2(Mst1)*pow4(Mst2) + 125*(109760*OepS2 - 9*(427143 + 573496*S2))*
        pow6(Mst1) + 111132000*pow6(Mst2)) + pow4(Mst2)*(12*Dmglst2*Mgl*pow2(
        Mst1)*(5488*(608404 + 75000*OepS2 - 188325*S2)*pow2(Mst2)*pow4(Mst1) -
        617400*(-11674 + 120*B4 - 120*D3 + 60*DN + 17091*S2)*pow2(Mst1)*pow4(
        Mst2) + (-3650009897 + 2110136000*OepS2 - 52385772600*S2)*pow6(Mst1) +
        481572000*pow6(Mst2)) + 12*pow2(Dmglst2)*(-11025*(-2818497 + 10080*B4 -
        10080*D3 + 5040*DN - 7840*OepS2 + 1126548*S2)*pow4(Mst1)*pow4(Mst2) +
        49*(273721621 + 16716000*OepS2 - 343302300*S2)*pow2(Mst2)*pow6(Mst1) +
        3797010000*pow2(Mst1)*pow6(Mst2) + (-3650009897 + 2110136000*OepS2 -
        52385772600*S2)*pow8(Mst1) - 296352000*pow8(Mst2)) + pow2(Mgl)*(77175*(
        1013263 + 2880*D3 - 2880*DN - 25440*OepS2 - 137052*S2)*pow4(Mst1)*pow4(
        Mst2) - 10290*(-2367149 + 21600*D3 - 21600*DN + 503200*OepS2 -
        14889420*S2)*pow2(Mst2)*pow6(Mst1) - 4834242000*pow2(Mst1)*pow6(Mst2) +
        (50004748544 - 9407804000*OepS2 + 369631514700*S2)*pow8(Mst1) +
        889056000*pow8(Mst2))) + 264600*pow2(log(pow2(Mst2)/pow2(Mst1)))*(90*
        pow2(Dmsqst2)*pow2(Mgl)*(45*pow2(Mst1) + 49*pow2(Mst2))*pow6(Mst1) -
        90*Dmsqst2*pow2(Mgl)*pow2(Mst2)*(37*pow2(Mst1) + 84*pow2(Mst2))*pow6(
        Mst1) + pow4(Mst2)*(4*Dmglst2*Mgl*(31605*pow4(Mst1)*pow4(Mst2) + 11368*
        pow2(Mst2)*pow6(Mst1) - 9555*pow2(Mst1)*pow6(Mst2) + 52512*pow8(Mst1) +
        3360*pow8(Mst2)) + pow2(Mgl)*(55860*pow4(Mst1)*pow4(Mst2) + 67060*pow2(
        Mst2)*pow6(Mst1) - 15750*pow2(Mst1)*pow6(Mst2) + 73657*pow8(Mst1) +
        3360*pow8(Mst2)) + 2*pow2(Dmglst2)*(49665*pow4(Mst1)*pow4(Mst2) -
        94318*pow2(Mst2)*pow6(Mst1) - 35385*pow2(Mst1)*pow6(Mst2) + 105024*
        pow8(Mst1) + 16800*pow8(Mst2)))) - 630*log(pow2(Mst2)/pow2(Mst1))*(15*
        pow2(Dmsqst2)*pow2(Mgl)*pow4(Mst1)*(392*(-887 + 2250*S2)*pow2(Mst1)*
        pow2(Mst2) + 25*(-20869 + 35280*S2)*pow4(Mst1) + 58800*(2 + 9*S2)*pow4(
        Mst2)) + 20*Dmsqst2*pow2(Mst1)*pow2(Mst2)*(11025*Dmglst2*Mgl*(37*pow2(
        Mst1) - 24*pow2(Mst2))*pow4(Mst1) + 11025*pow2(Dmglst2)*(-24*pow2(Mst2)
        *pow4(Mst1) + 37*pow6(Mst1)) - 4*pow2(Mgl)*(441*(-517 + 450*S2)*pow2(
        Mst2)*pow4(Mst1) + 11025*(16 + 27*S2)*pow2(Mst1)*pow4(Mst2) + 2*(-58853
        + 44100*S2)*pow6(Mst1) - 66150*pow6(Mst2))) + 2*pow4(Mst2)*(4*pow2(
        Dmglst2)*(3675*(13339 + 1134*S2)*pow4(Mst1)*pow4(Mst2) + 49*(463453 +
        805950*S2)*pow2(Mst2)*pow6(Mst1) - 1830150*pow2(Mst1)*pow6(Mst2) + 3*(-
        4261807 + 33912900*S2)*pow8(Mst1) + 352800*pow8(Mst2)) + 4*Dmglst2*Mgl*
        (25188450*pow4(Mst1)*pow4(Mst2) + 98*(77657 + 202500*S2)*pow2(Mst2)*
        pow6(Mst1) - 2954700*pow2(Mst1)*pow6(Mst2) + 3*(-4261807 + 33912900*S2)
        *pow8(Mst1) + 705600*pow8(Mst2)) + pow2(Mgl)*(-7350*(-4534 + 4293*S2)*
        pow4(Mst1)*pow4(Mst2) - 980*(-52322 + 84915*S2)*pow2(Mst2)*pow6(Mst1) -
        6791400*pow2(Mst1)*pow6(Mst2) + (58223589 - 151196850*S2)*pow8(Mst1) +
        1411200*pow8(Mst2))))))/(pow2(Mgl)*pow4(Mst1)) - (2450*z2*(11340*log(
        pow2(Mst2)/pow2(Mst1))*pow2(Mst1)*(15*pow2(Dmsqst2)*pow2(Mgl)*(32*Mt*(
        5*Mt - 4*Mst2*s2t)*pow2(Mst1) + 16*Mt*(3*Mt - 4*Mst2*s2t)*pow2(Mst2) +
        3*pow2(s2t)*pow4(Mst1)) + 60*Dmsqst2*(Dmglst2*Mgl*Mst2*(32*Mst2*Mt*(-9*
        Mt + 7*Mst2*s2t)*pow2(Mst1) + 16*Mt*(-5*Mt + 6*Mst2*s2t)*pow3(Mst2) +
        s2t*(352*Mt - 5*Mst2*s2t)*pow4(Mst1)) + Mst2*pow2(Dmglst2)*(32*Mst2*Mt*
        (-9*Mt + 7*Mst2*s2t)*pow2(Mst1) + 16*Mt*(-5*Mt + 6*Mst2*s2t)*pow3(Mst2)
        + s2t*(352*Mt - 5*Mst2*s2t)*pow4(Mst1)) + 16*Mt*pow2(Mgl)*(2*(Mt -
        Mst2*s2t)*pow2(Mst1)*pow2(Mst2) + (3*Mt - 2*Mst2*s2t)*pow4(Mst1) + (Mt
        - 2*Mst2*s2t)*pow4(Mst2))) - 2*pow2(Mst2)*(4*Dmglst2*Mgl*(pow2(Mst1)*
        pow2(Mst2)*(-985*Mst2*Mt*s2t + 1836*pow2(Mt) + 124*pow2(Mst2)*pow2(s2t)
        ) + (-2184*Mst2*Mt*s2t + 4848*pow2(Mt) + 5*pow2(Mst2)*pow2(s2t))*pow4(
        Mst1) + (-181*Mst2*Mt*s2t + 392*pow2(Mt) + 297*pow2(Mst2)*pow2(s2t))*
        pow4(Mst2)) + pow2(Mgl)*(4*pow2(Mst1)*pow2(Mst2)*(495*Mst2*Mt*s2t -
        403*pow2(Mt) + 17*pow2(Mst2)*pow2(s2t)) + (2432*Mst2*Mt*s2t - 2928*
        pow2(Mt) + 195*pow2(Mst2)*pow2(s2t))*pow4(Mst1) + (1516*Mst2*Mt*s2t -
        894*pow2(Mt) + 305*pow2(Mst2)*pow2(s2t))*pow4(Mst2)) + 2*pow2(Dmglst2)*
        (2*pow2(Mst1)*pow2(Mst2)*(1771*Mst2*Mt*s2t - 2998*pow2(Mt) + 703*pow2(
        Mst2)*pow2(s2t)) + 2*(-2184*Mst2*Mt*s2t + 4848*pow2(Mt) + 5*pow2(Mst2)*
        pow2(s2t))*pow4(Mst1) + (1118*Mst2*Mt*s2t - 664*pow2(Mt) + 1645*pow2(
        Mst2)*pow2(s2t))*pow4(Mst2)))) + 3150*pow2(Dmsqst2)*pow2(Mgl)*(6*pow2(
        Mst1)*pow2(Mst2)*(720*Mst2*Mt*s2t - 936*pow2(Mt) + 173*pow2(Mst2)*pow2(
        s2t)) - 2*(-6048*Mst2*Mt*s2t + 9720*pow2(Mt) - 721*pow2(Mst2)*pow2(s2t)
        )*pow4(Mst1) + 65*pow2(s2t)*pow6(Mst1) + 216*(-2*pow2(Mt)*pow4(Mst2) +
        pow2(s2t)*pow6(Mst2))) - 2100*Dmsqst2*(162*Dmglst2*Mgl*Mst2*(16*pow2(
        Mst1)*(19*Mst2*Mt*s2t - 23*pow2(Mt) + 2*pow2(Mst2)*pow2(s2t))*pow3(
        Mst2) + 16*Mst2*(51*Mst2*Mt*s2t - 77*pow2(Mt) + 2*pow2(Mst2)*pow2(s2t))
        *pow4(Mst1) + 8*(-2*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow5(Mst2) + s2t*(
        1392*Mt - 61*Mst2*s2t)*pow6(Mst1)) + 162*Mst2*pow2(Dmglst2)*(16*pow2(
        Mst1)*(19*Mst2*Mt*s2t - 23*pow2(Mt) + 2*pow2(Mst2)*pow2(s2t))*pow3(
        Mst2) + 16*Mst2*(51*Mst2*Mt*s2t - 77*pow2(Mt) + 2*pow2(Mst2)*pow2(s2t))
        *pow4(Mst1) + 8*(-2*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow5(Mst2) + s2t*(
        1392*Mt - 61*Mst2*s2t)*pow6(Mst1)) + pow2(Mgl)*(18*pow2(Mst2)*(-720*
        Mst2*Mt*s2t + 864*pow2(Mt) + 101*pow2(Mst2)*pow2(s2t))*pow4(Mst1) + 27*
        pow2(Mst1)*(-96*Mst2*Mt*s2t + 288*pow2(Mt) + 77*pow2(Mst2)*pow2(s2t))*
        pow4(Mst2) + 8*(-1620*Mst2*Mt*s2t + 2268*pow2(Mt) + 101*pow2(Mst2)*
        pow2(s2t))*pow6(Mst1) - 648*(-2*pow2(Mt)*pow6(Mst2) + pow2(s2t)*pow8(
        Mst2)))) + pow2(Mst2)*(4*Dmglst2*Mgl*(9*pow2(Mst2)*(-10713129*Mst2*Mt*
        s2t + 12737588*pow2(Mt) + 1272600*pow2(Mst2)*pow2(s2t))*pow4(Mst1) +
        27*pow2(Mst1)*(-439845*Mst2*Mt*s2t + 837364*pow2(Mt) + 135772*pow2(
        Mst2)*pow2(s2t))*pow4(Mst2) + (-277428676*Mst2*Mt*s2t + 306443792*pow2(
        Mt) + 49475790*pow2(Mst2)*pow2(s2t))*pow6(Mst1) - 317520*(-2*pow2(Mt)*
        pow6(Mst2) + pow2(s2t)*pow8(Mst2))) - 21*pow2(Mgl)*(6*pow2(Mst2)*(-
        941314*Mst2*Mt*s2t + 620078*pow2(Mt) + 309749*pow2(Mst2)*pow2(s2t))*
        pow4(Mst1) + 9*pow2(Mst1)*(-274772*Mst2*Mt*s2t + 195932*pow2(Mt) +
        5927*pow2(Mst2)*pow2(s2t))*pow4(Mst2) + (-9196144*Mst2*Mt*s2t +
        5014952*pow2(Mt) + 3291029*pow2(Mst2)*pow2(s2t))*pow6(Mst1) - 41040*(-
        2*pow2(Mt)*pow6(Mst2) + pow2(s2t)*pow8(Mst2))) + 4*pow2(Dmglst2)*(3*
        pow2(Mst2)*(47845349*Mst2*Mt*s2t - 105257590*pow2(Mt) + 12219435*pow2(
        Mst2)*pow2(s2t))*pow4(Mst1) + 9*pow2(Mst1)*(2661049*Mst2*Mt*s2t -
        4633370*pow2(Mt) + 1109043*pow2(Mst2)*pow2(s2t))*pow4(Mst2) + (-
        277428676*Mst2*Mt*s2t + 306443792*pow2(Mt) + 49475790*pow2(Mst2)*pow2(
        s2t))*pow6(Mst1) + 272160*(-2*pow2(Mt)*pow6(Mst2) + pow2(s2t)*pow8(
        Mst2))))))/(pow2(Mgl)*pow2(Mst1))))/pow6(Mst2)))/6.001128e9)/pow4(Mt)
        *12.;

   return result;
}

/**
 * 	@return returns the susy log^1 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double himalaya::H6b2qg2::calc_coef_at_as2_no_sm_logs_log1() const {

   const double result =
     (((6000*pow2(Dmsqst2)*pow2(Mgl)*(72*(2*Mst2*s2t*(3 - 2*z2) +
        Mt*(-7 + 5*z2))*pow2(Mst1) + (36*Mst2*s2t*(7 - 4*z2) + Mt*(-161 + 108*
        z2))*pow2(Mst2))*pow3(Mt)*pow4(Mst1) + 450*pow2(Mst2)*pow2(log(pow2(
        Mst2)/pow2(Mst1)))*pow4(Mst1)*(-8*Dmglst2*Mgl*(8*Mt*(-157*Mst2*s2t*
        pow2(Mt) - 108*Mt*pow2(Mst2)*pow2(s2t) + 398*pow3(Mt) - 16*pow3(Mst2)*
        pow3(s2t))*pow4(Mst1) + pow4(Mst2)*(-320*pow2(Mst2)*pow2(Mt)*pow2(s2t)
        + 848*Mst2*s2t*pow3(Mt) - 532*Mt*pow3(Mst2)*pow3(s2t) + 512*pow4(Mt) -
        99*pow4(Mst2)*pow4(s2t)) + pow2(Mst1)*pow2(Mst2)*(-768*pow2(Mst2)*pow2(
        Mt)*pow2(s2t) - 448*Mst2*s2t*pow3(Mt) - 212*Mt*pow3(Mst2)*pow3(s2t) +
        1896*pow4(Mt) + 35*pow4(Mst2)*pow4(s2t))) - 4*pow2(Dmglst2)*(16*Mt*(-
        157*Mst2*s2t*pow2(Mt) - 108*Mt*pow2(Mst2)*pow2(s2t) + 398*pow3(Mt) -
        16*pow3(Mst2)*pow3(s2t))*pow4(Mst1) + pow4(Mst2)*(-960*pow2(Mst2)*pow2(
        Mt)*pow2(s2t) + 1696*Mst2*s2t*pow3(Mt) - 1832*Mt*pow3(Mst2)*pow3(s2t) +
        1536*pow4(Mt) - 297*pow4(Mst2)*pow4(s2t)) + pow2(Mst1)*pow2(Mst2)*(-
        2304*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 256*Mst2*s2t*pow3(Mt) - 424*Mt*
        pow3(Mst2)*pow3(s2t) + 152*pow4(Mt) + 105*pow4(Mst2)*pow4(s2t))) +
        pow2(Mgl)*(8*pow2(Mst1)*pow2(Mst2)*(425*pow2(Mst2)*pow2(Mt)*pow2(s2t) -
        704*Mst2*s2t*pow3(Mt) + 212*Mt*pow3(Mst2)*pow3(s2t) - 76*pow4(Mt) + 3*
        pow4(Mst2)*pow4(s2t)) + pow4(Mst1)*(2688*pow2(Mst2)*pow2(Mt)*pow2(s2t)
        - 5056*Mst2*s2t*pow3(Mt) + 1024*Mt*pow3(Mst2)*pow3(s2t) - 2848*pow4(Mt)
        + 41*pow4(Mst2)*pow4(s2t)) + pow4(Mst2)*(296*pow2(Mst2)*pow2(Mt)*pow2(
        s2t) - 6784*Mst2*s2t*pow3(Mt) + 2208*Mt*pow3(Mst2)*pow3(s2t) - 672*
        pow4(Mt) + 519*pow4(Mst2)*pow4(s2t)))) + 3000*Dmsqst2*pow2(Mst1)*(8*
        Dmglst2*Mgl*Mst2*pow2(Mst1)*pow3(Mt)*(36*Mst2*(Mt*(33 - 18*z2) + 2*
        Mst2*s2t*(-12 + 7*z2))*pow2(Mst1) + (Mt*(283 - 180*z2) + 12*Mst2*s2t*(-
        29 + 18*z2))*pow3(Mst2) + 18*s2t*(-83 + 44*z2)*pow4(Mst1)) + 8*Mst2*
        pow2(Dmglst2)*pow2(Mst1)*pow3(Mt)*(36*Mst2*(Mt*(33 - 18*z2) + 2*Mst2*
        s2t*(-12 + 7*z2))*pow2(Mst1) + (Mt*(283 - 180*z2) + 12*Mst2*s2t*(-29 +
        18*z2))*pow3(Mst2) + 18*s2t*(-83 + 44*z2)*pow4(Mst1)) + 3*pow2(Mgl)*(3*
        pow2(Mst2)*pow4(Mst1)*(8*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 64*Mst2*s2t*(-
        2 + z2)*pow3(Mt) + 16*(-9 + 4*z2)*pow4(Mt) + pow4(Mst2)*pow4(s2t)) +
        pow2(Mst1)*pow4(Mst2)*(-48*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 96*Mst2*s2t*
        (-3 + 2*z2)*pow3(Mt) + 32*(-8 + 3*z2)*pow4(Mt) + 3*pow4(Mst2)*pow4(s2t)
        ) + 3*(-16*Mst2*s2t*(-9 + 4*z2)*pow3(Mt) + 8*(-31 + 12*z2)*pow4(Mt) -
        pow4(Mst2)*pow4(s2t))*pow6(Mst1) - 3*pow2(-4*pow2(Mt) + pow2(Mst2)*
        pow2(s2t))*pow6(Mst2))) + pow2(Mst2)*(-100*Dmglst2*Mgl*(2*pow4(Mst1)*
        pow4(Mst2)*(12*(-3565 + 1728*z2)*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 8*
        Mst2*s2t*(-3659 + 252*z2)*pow3(Mt) + 72*Mt*(-192 + 91*z2)*pow3(Mst2)*
        pow3(s2t) + 44*(427 + 432*z2)*pow4(Mt) + 9*(211 - 86*z2)*pow4(Mst2)*
        pow4(s2t)) + 6*pow2(Mst2)*(4*(-1955 + 1152*z2)*pow2(Mst2)*pow2(Mt)*
        pow2(s2t) - 48*Mst2*s2t*(-13 + 262*z2)*pow3(Mt) + 24*Mt*(111 - 17*z2)*
        pow3(Mst2)*pow3(s2t) + 8*(-1471 + 3240*z2)*pow4(Mt) + 3*(-1 + 86*z2)*
        pow4(Mst2)*pow4(s2t))*pow6(Mst1) - 18*pow2(Mst1)*(-536*pow2(Mst2)*pow2(
        Mt)*pow2(s2t) - 512*Mst2*s2t*pow3(Mt) + 128*Mt*pow3(Mst2)*pow3(s2t) +
        560*pow4(Mt) + 131*pow4(Mst2)*pow4(s2t))*pow6(Mst2) + 3*(192*(-77 + 48*
        z2)*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 48*Mst2*s2t*(-867 + 1012*z2)*pow3(
        Mt) - 276*Mt*pow3(Mst2)*pow3(s2t) + 16*(-6445 + 6768*z2)*pow4(Mt) -
        441*pow4(Mst2)*pow4(s2t))*pow8(Mst1) - 2304*pow2(Mt)*(-2*pow2(Mt) +
        pow2(Mst2)*pow2(s2t))*pow8(Mst2)) + 25*pow2(Mgl)*(4*pow4(Mst1)*pow4(
        Mst2)*(-72*(-445 + 96*z2)*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 144*Mst2*s2t*
        (-121 + 296*z2)*pow3(Mt) - 144*Mt*(-151 + 19*z2)*pow3(Mst2)*pow3(s2t) +
        8*(-1993 + 2151*z2 + 216*z3)*pow4(Mt) + 81*(8 + 5*z2)*pow4(Mst2)*pow4(
        s2t)) - 36*pow2(Mst2)*(-696*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 32*Mst2*
        s2t*(-179 + 128*z2)*pow3(Mt) + 16*Mt*(148 - 17*z2)*pow3(Mst2)*pow3(s2t)
        - 8*(-403 + 300*z2)*pow4(Mt) + (551 + 86*z2)*pow4(Mst2)*pow4(s2t))*
        pow6(Mst1) + 36*pow2(Mst1)*(-744*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 1024*
        Mst2*s2t*pow3(Mt) + 256*Mt*pow3(Mst2)*pow3(s2t) + 1232*pow4(Mt) + 141*
        pow4(Mst2)*pow4(s2t))*pow6(Mst2) + 9*(864*pow2(Mst2)*pow2(Mt)*pow2(s2t)
        - 512*Mst2*s2t*(-49 + 31*z2)*pow3(Mt) - 784*Mt*pow3(Mst2)*pow3(s2t) +
        8*(-2503 + 1408*z2)*pow4(Mt) + (1547 + 164*z2)*pow4(Mst2)*pow4(s2t))*
        pow8(Mst1) + 4608*pow2(Mt)*(-2*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow8(
        Mst2)) - 4*pow2(Dmglst2)*(-2*pow4(Mst1)*pow4(Mst2)*(-150*(-14251 +
        9792*z2)*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 600*Mst2*s2t*(-2509 + 1278*z2)
        *pow3(Mt) - 300*Mt*(-2027 + 1194*z2)*pow3(Mst2)*pow3(s2t) + 4*(25949 +
        33300*z2)*pow4(Mt) + 75*(-838 + 387*z2)*pow4(Mst2)*pow4(s2t)) + 25*
        pow2(Mst2)*(12*(-11813 + 8064*z2)*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 96*
        Mst2*s2t*(-1016 + 1395*z2)*pow3(Mt) - 24*Mt*(-83 + 102*z2)*pow3(Mst2)*
        pow3(s2t) - 64*(347 + 3330*z2)*pow4(Mt) + 3*(-415 + 774*z2)*pow4(Mst2)*
        pow4(s2t))*pow6(Mst1) - 225*pow2(Mst1)*(-408*pow2(Mst2)*pow2(Mt)*pow2(
        s2t) + 1024*Mst2*s2t*pow3(Mt) - 256*Mt*pow3(Mst2)*pow3(s2t) - 720*pow4(
        Mt) + 179*pow4(Mst2)*pow4(s2t))*pow6(Mst2) + 75*(192*(-77 + 48*z2)*
        pow2(Mst2)*pow2(Mt)*pow2(s2t) - 48*Mst2*s2t*(-867 + 1012*z2)*pow3(Mt) -
        276*Mt*pow3(Mst2)*pow3(s2t) + 16*(-6445 + 6768*z2)*pow4(Mt) - 441*pow4(
        Mst2)*pow4(s2t))*pow8(Mst1) - 3600*(8*pow2(Mst2)*pow2(Mt)*pow2(s2t) -
        16*pow4(Mt) + pow4(Mst2)*pow4(s2t))*pow8(Mst2))) - 30*log(pow2(Mst2)/
        pow2(Mst1))*(7200*pow2(Dmsqst2)*pow2(Mgl)*(pow2(Mst1) - pow2(Mst2))*
        pow4(Mst1)*pow4(Mt) - 1800*Dmsqst2*pow4(Mst1)*(16*Dmglst2*Mgl*Mst2*
        pow2(Mst1)*(-(Mst2*Mt) + 3*s2t*pow2(Mst1))*pow3(Mt) + 16*Mst2*pow2(
        Dmglst2)*pow2(Mst1)*(-(Mst2*Mt) + 3*s2t*pow2(Mst1))*pow3(Mt) + pow2(
        Mgl)*(8*(5*Mt - 2*Mst2*s2t)*pow3(Mt)*pow4(Mst1) - 8*pow4(Mst2)*pow4(Mt)
        + pow2(Mst1)*(8*pow2(Mst2)*pow4(Mt) - pow4(s2t)*pow6(Mst2)) + pow4(s2t)
        *pow8(Mst2))) + pow2(Mst2)*(-20*Dmglst2*Mgl*(-4*pow4(Mst1)*pow4(Mst2)*(
        1860*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 1472*Mst2*s2t*pow3(Mt) + 774*Mt*
        pow3(Mst2)*pow3(s2t) - 2464*pow4(Mt) + 15*pow4(Mst2)*pow4(s2t)) + 2*
        pow2(Mst2)*(-432*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 11656*Mst2*s2t*pow3(
        Mt) + 996*Mt*pow3(Mst2)*pow3(s2t) + 18748*pow4(Mt) + 207*pow4(Mst2)*
        pow4(s2t))*pow6(Mst1) - 3*pow2(Mst1)*(-856*pow2(Mst2)*pow2(Mt)*pow2(
        s2t) - 1536*Mst2*s2t*pow3(Mt) + 384*Mt*pow3(Mst2)*pow3(s2t) + 1200*
        pow4(Mt) + 203*pow4(Mst2)*pow4(s2t))*pow6(Mst2) + 4*(-738*pow2(Mst2)*
        pow2(Mt)*pow2(s2t) - 10004*Mst2*s2t*pow3(Mt) + 135*Mt*pow3(Mst2)*pow3(
        s2t) + 17532*pow4(Mt) + 12*pow4(Mst2)*pow4(s2t))*pow8(Mst1) - 768*pow2(
        Mt)*(-2*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow8(Mst2)) + 5*pow2(Mgl)*(4*
        pow4(Mst1)*pow4(Mst2)*(3372*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 9024*Mst2*
        s2t*pow3(Mt) + 3912*Mt*pow3(Mst2)*pow3(s2t) + 6112*pow4(Mt) + 579*pow4(
        Mst2)*pow4(s2t)) - 24*pow2(Mst2)*(-534*pow2(Mst2)*pow2(Mt)*pow2(s2t) +
        1824*Mst2*s2t*pow3(Mt) - 60*Mt*pow3(Mst2)*pow3(s2t) - 934*pow4(Mt) +
        97*pow4(Mst2)*pow4(s2t))*pow6(Mst1) + 6*pow2(Mst1)*(-728*pow2(Mst2)*
        pow2(Mt)*pow2(s2t) - 1024*Mst2*s2t*pow3(Mt) + 256*Mt*pow3(Mst2)*pow3(
        s2t) + 1200*pow4(Mt) + 139*pow4(Mst2)*pow4(s2t))*pow6(Mst2) + 3*(4640*
        pow2(Mst2)*pow2(Mt)*pow2(s2t) - 15296*Mst2*s2t*pow3(Mt) + 240*Mt*pow3(
        Mst2)*pow3(s2t) + 7128*pow4(Mt) - 279*pow4(Mst2)*pow4(s2t))*pow8(Mst1)
        + 768*pow2(Mt)*(-2*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow8(Mst2)) - 2*
        pow2(Dmglst2)*(2*pow4(Mst1)*pow4(Mst2)*(-40440*pow2(Mst2)*pow2(Mt)*
        pow2(s2t) - 75360*Mst2*s2t*pow3(Mt) - 5760*Mt*pow3(Mst2)*pow3(s2t) +
        11872*pow4(Mt) + 4335*pow4(Mst2)*pow4(s2t)) + 2*pow2(Mst2)*(29760*pow2(
        Mst2)*pow2(Mt)*pow2(s2t) + 31000*Mst2*s2t*pow3(Mt) + 15840*Mt*pow3(
        Mst2)*pow3(s2t) - 228452*pow4(Mt) + 1065*pow4(Mst2)*pow4(s2t))*pow6(
        Mst1) - 15*pow2(Mst1)*(-3080*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 6144*Mst2*
        s2t*pow3(Mt) + 1536*Mt*pow3(Mst2)*pow3(s2t) + 3600*pow4(Mt) + 865*pow4(
        Mst2)*pow4(s2t))*pow6(Mst2) + 40*(-738*pow2(Mst2)*pow2(Mt)*pow2(s2t) -
        10004*Mst2*s2t*pow3(Mt) + 135*Mt*pow3(Mst2)*pow3(s2t) + 17532*pow4(Mt)
        + 12*pow4(Mst2)*pow4(s2t))*pow8(Mst1) + 480*(-40*pow2(Mst2)*pow2(Mt)*
        pow2(s2t) + 80*pow4(Mt) + 3*pow4(Mst2)*pow4(s2t))*pow8(Mst2))))))/(
        16200.*pow2(Mgl)*pow4(Mst1)*pow6(Mst2)))/pow4(Mt)*12.; 

   return result;
}

/**
 * 	@return returns the susy log^2 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double himalaya::H6b2qg2::calc_coef_at_as2_no_sm_logs_log2() const {

   const double result =
      (((-3600*pow2(Dmsqst2)*pow2(Mgl)*pow4(Mst1)*pow4(Mt) + 7200*
        Dmsqst2*pow2(Mgl)*pow2(Mst2)*pow4(Mst1)*pow4(Mt) - 48160*Dmglst2*Mgl*
        pow4(Mst1)*pow4(Mst2)*pow4(Mt) - 29264*pow2(Dmglst2)*pow4(Mst1)*pow4(
        Mst2)*pow4(Mt) + 21280*pow2(Mgl)*pow4(Mst1)*pow4(Mst2)*pow4(Mt) - 1440*
        pow12(Mst2)*pow2(Dmglst2)*pow4(s2t) + 46400*Dmglst2*Mgl*s2t*pow3(Mt)*
        pow4(Mst1)*pow5(Mst2) + 98880*s2t*pow2(Dmglst2)*pow3(Mt)*pow4(Mst1)*
        pow5(Mst2) - 24000*s2t*pow2(Mgl)*pow3(Mt)*pow4(Mst1)*pow5(Mst2) +
        74880*Dmglst2*Mgl*s2t*pow3(Mst2)*pow3(Mt)*pow6(Mst1) - 38400*s2t*pow2(
        Dmglst2)*pow3(Mst2)*pow3(Mt)*pow6(Mst1) + 1920*s2t*pow2(Mgl)*pow3(Mst2)
        *pow3(Mt)*pow6(Mst1) - 14160*Dmglst2*Mgl*pow2(Mt)*pow2(s2t)*pow4(Mst2)*
        pow6(Mst1) - 13560*pow2(Dmglst2)*pow2(Mt)*pow2(s2t)*pow4(Mst2)*pow6(
        Mst1) - 9000*pow2(Mgl)*pow2(Mt)*pow2(s2t)*pow4(Mst2)*pow6(Mst1) -
        100800*Dmglst2*Mgl*pow2(Mst2)*pow4(Mt)*pow6(Mst1) + 161120*pow2(
        Dmglst2)*pow2(Mst2)*pow4(Mt)*pow6(Mst1) - 10080*pow2(Mgl)*pow2(Mst2)*
        pow4(Mt)*pow6(Mst1) - 31200*Dmglst2*Mgl*Mt*pow3(s2t)*pow5(Mst2)*pow6(
        Mst1) - 42720*Mt*pow2(Dmglst2)*pow3(s2t)*pow5(Mst2)*pow6(Mst1) - 23520*
        Mt*pow2(Mgl)*pow3(s2t)*pow5(Mst2)*pow6(Mst1) + 82080*Dmglst2*Mgl*pow2(
        Mt)*pow2(s2t)*pow4(Mst1)*pow6(Mst2) + 115440*pow2(Dmglst2)*pow2(Mt)*
        pow2(s2t)*pow4(Mst1)*pow6(Mst2) + 42960*pow2(Mgl)*pow2(Mt)*pow2(s2t)*
        pow4(Mst1)*pow6(Mst2) + 43680*Dmglst2*Mgl*pow2(Mst1)*pow4(Mt)*pow6(
        Mst2) + 65520*pow2(Dmglst2)*pow2(Mst1)*pow4(Mt)*pow6(Mst2) + 21840*
        pow2(Mgl)*pow2(Mst1)*pow4(Mt)*pow6(Mst2) + 1110*Dmglst2*Mgl*pow4(s2t)*
        pow6(Mst1)*pow6(Mst2) + 5505*pow2(Dmglst2)*pow4(s2t)*pow6(Mst1)*pow6(
        Mst2) - 5325*pow2(Mgl)*pow4(s2t)*pow6(Mst1)*pow6(Mst2) + 30*log(pow2(
        Mst2)/pow2(Mst1))*pow4(Mst1)*(pow2(Mgl)*(-4*pow4(Mst2)*(14*pow2(Mst2)*
        pow2(Mt)*pow2(s2t) - 456*Mst2*s2t*pow3(Mt) + 146*Mt*pow3(Mst2)*pow3(
        s2t) - 166*pow4(Mt) + 33*pow4(Mst2)*pow4(s2t)) + pow4(Mst1)*(-768*pow2(
        Mst2)*pow2(Mt)*pow2(s2t) + 1216*Mst2*s2t*pow3(Mt) + 400*pow4(Mt) + 41*
        pow4(Mst2)*pow4(s2t)) + pow2(Mst1)*pow2(Mst2)*(-1096*pow2(Mst2)*pow2(
        Mt)*pow2(s2t) + 1280*Mst2*s2t*pow3(Mt) - 328*Mt*pow3(Mst2)*pow3(s2t) +
        32*pow4(Mt) + 91*pow4(Mst2)*pow4(s2t))) + 2*Dmglst2*Mgl*(32*pow2(Mt)*(-
        57*Mst2*Mt*s2t + 119*pow2(Mt) - 24*pow2(Mst2)*pow2(s2t))*pow4(Mst1) +
        pow4(Mst2)*(-384*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 912*Mst2*s2t*pow3(Mt)
        - 548*Mt*pow3(Mst2)*pow3(s2t) + 512*pow4(Mt) - 91*pow4(Mst2)*pow4(s2t))
        + pow2(Mst1)*pow2(Mst2)*(-768*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 640*Mst2*
        s2t*pow3(Mt) - 164*Mt*pow3(Mst2)*pow3(s2t) + 2016*pow4(Mt) + 91*pow4(
        Mst2)*pow4(s2t))) + pow2(Dmglst2)*(64*pow2(Mt)*(-57*Mst2*Mt*s2t + 119*
        pow2(Mt) - 24*pow2(Mst2)*pow2(s2t))*pow4(Mst1) + pow4(Mst2)*(-1152*
        pow2(Mst2)*pow2(Mt)*pow2(s2t) + 1824*Mst2*s2t*pow3(Mt) - 1864*Mt*pow3(
        Mst2)*pow3(s2t) + 1536*pow4(Mt) - 273*pow4(Mst2)*pow4(s2t)) + pow2(
        Mst1)*(-2304*pow2(Mt)*pow2(s2t)*pow4(Mst2) + 32*pow2(Mst2)*pow4(Mt) -
        328*Mt*pow3(s2t)*pow5(Mst2) + 273*pow4(s2t)*pow6(Mst2)))) - 46080*
        Dmglst2*Mgl*s2t*pow2(Mst1)*pow3(Mt)*pow7(Mst2) - 92160*s2t*pow2(
        Dmglst2)*pow2(Mst1)*pow3(Mt)*pow7(Mst2) - 15360*s2t*pow2(Mgl)*pow2(
        Mst1)*pow3(Mt)*pow7(Mst2) + 19680*Dmglst2*Mgl*Mt*pow3(s2t)*pow4(Mst1)*
        pow7(Mst2) + 19680*Mt*pow2(Dmglst2)*pow3(s2t)*pow4(Mst1)*pow7(Mst2) +
        19680*Mt*pow2(Mgl)*pow3(s2t)*pow4(Mst1)*pow7(Mst2) + 67200*Dmglst2*Mgl*
        Mst2*s2t*pow3(Mt)*pow8(Mst1) + 67200*Mst2*s2t*pow2(Dmglst2)*pow3(Mt)*
        pow8(Mst1) + 1920*Mst2*s2t*pow2(Mgl)*pow3(Mt)*pow8(Mst1) - 112320*
        Dmglst2*Mgl*pow4(Mt)*pow8(Mst1) - 112320*pow2(Dmglst2)*pow4(Mt)*pow8(
        Mst1) - 12000*pow2(Mgl)*pow4(Mt)*pow8(Mst1) + 1770*Dmglst2*Mgl*pow4(
        Mst2)*pow4(s2t)*pow8(Mst1) + 1770*pow2(Dmglst2)*pow4(Mst2)*pow4(s2t)*
        pow8(Mst1) + 3585*pow2(Mgl)*pow4(Mst2)*pow4(s2t)*pow8(Mst1) - 29520*
        Dmglst2*Mgl*pow2(Mst1)*pow2(Mt)*pow2(s2t)*pow8(Mst2) - 51960*pow2(
        Dmglst2)*pow2(Mst1)*pow2(Mt)*pow2(s2t)*pow8(Mst2) - 12840*pow2(Mgl)*
        pow2(Mst1)*pow2(Mt)*pow2(s2t)*pow8(Mst2) - 15360*Dmglst2*Mgl*pow4(Mt)*
        pow8(Mst2) - 38400*pow2(Dmglst2)*pow4(Mt)*pow8(Mst2) - 3840*pow2(Mgl)*
        pow4(Mt)*pow8(Mst2) - 8490*Dmglst2*Mgl*pow4(Mst1)*pow4(s2t)*pow8(Mst2)
        - 18495*pow2(Dmglst2)*pow4(Mst1)*pow4(s2t)*pow8(Mst2) - 345*pow2(Mgl)*
        pow4(Mst1)*pow4(s2t)*pow8(Mst2) + 11520*Dmglst2*Mgl*Mt*pow2(Mst1)*pow3(
        s2t)*pow9(Mst2) + 23040*Mt*pow2(Dmglst2)*pow2(Mst1)*pow3(s2t)*pow9(
        Mst2) + 3840*Mt*pow2(Mgl)*pow2(Mst1)*pow3(s2t)*pow9(Mst2) + 7680*
        Dmglst2*Mgl*pow2(Mt)*pow2(s2t)*power10(Mst2) + 19200*pow2(Dmglst2)*
        pow2(Mt)*pow2(s2t)*power10(Mst2) + 1920*pow2(Mgl)*pow2(Mt)*pow2(s2t)*
        power10(Mst2) + 6570*Dmglst2*Mgl*pow2(Mst1)*pow4(s2t)*power10(Mst2) +
        13695*pow2(Dmglst2)*pow2(Mst1)*pow4(s2t)*power10(Mst2) + 2325*pow2(Mgl)
        *pow2(Mst1)*pow4(s2t)*power10(Mst2)))/(540.*pow2(Mgl)*pow4(Mst1)*pow4(
        Mst2)))/pow4(Mt)*12.; 

   return result;
}

/**
 * 	@return returns the susy log^3 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double himalaya::H6b2qg2::calc_coef_at_as2_no_sm_logs_log3() const {

   const double result =
      ((-224*pow4(Mt))/9.)/pow4(Mt)*12.; 

   return result;
}


