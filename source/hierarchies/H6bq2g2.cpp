#include "H6bq2g2.hpp"
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
 * 	@param Msq a double average squark mass w/o the stop quark
 * 	@param MuSUSY a double mu parameter
 * 	@param s2t a double 2 times the sine of the stop/sbottom quark mixing angle
 * 	@param mdrFlag an int 0 for DR and 1 for MDR scheme
 * 	@param oneLoopFlag an int flag to consider the one-loop expansion terms
 * 	@param twoLoopFlag an int flag to consider the two-loop expansion terms
 * 	@param threeLoopFlag an int flag to consider the three-loop expansion terms
 */
himalaya::H6bq2g2::H6bq2g2(const std::map<unsigned int, unsigned int>& flagMap, double Al4p, double beta, double Dmglst2,
		 double Dmsqst2, double lmMt, double lmMst1, double lmMst2,
		 double Mgl, double Mt, double Mst1, double Mst2, double Msq, double MuSUSY,
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
   this -> Msq = Msq;
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
   #include "../hierarchies/h6bq2g2/sigS1Full.inc"
   ;
   s2 = 
   #include "../hierarchies/h6bq2g2/sigS2Full.inc"
   ;
   s12 = 
   #include "../hierarchies/h6bq2g2/sigS12Full.inc"
   ;
}

/**
 * 	@return The diagonal (1, 1) matrix element of the Higgs mass matrix as a double for the hierarchy 'H6bq2g2'
 */
double himalaya::H6bq2g2::getS1() const {
   return s1;
}

/**
 * 	@return The diagonal (2, 2) matrix element of the Higgs mass matrix as a double for the hierarchy 'H6bq2g2'
 */
double himalaya::H6bq2g2::getS2() const {
   return s2;
}

/**
 * 	@return The off-diagonal (1, 2) = (2, 1) matrix element of the Higgs mass matrix as a double for the hierarchy 'H6bq2g2'
 */
double himalaya::H6bq2g2::getS12() const {
   return s12;
}

/**
 * 	@return returns the susy log^0 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double himalaya::H6bq2g2::calc_coef_at_as2_no_sm_logs_log0() const {

   const double result =
      ((-(T1ep*(2*Dmglst2*(2*Mgl*(pow4(Mst1)*(55368*pow2(Mst2)*pow2(
        Mt)*pow2(s2t) - 205840*Mst2*s2t*pow3(Mt) + 15571*Mt*pow3(Mst2)*pow3(
        s2t) + 89312*pow4(Mt) - 4199*pow4(Mst2)*pow4(s2t)) + 18*pow2(Mst1)*
        pow2(Mst2)*(600*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 2766*Mst2*s2t*pow3(Mt)
        + 401*Mt*pow3(Mst2)*pow3(s2t) + 1288*pow4(Mt) - 108*pow4(Mst2)*pow4(
        s2t)) - 189*pow4(Mst2)*(20*Mst2*s2t*pow3(Mt) - 5*Mt*pow3(Mst2)*pow3(
        s2t) - 16*pow4(Mt) + pow4(Mst2)*pow4(s2t))) + Dmglst2*(2*pow4(Mst1)*(
        55368*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 205840*Mst2*s2t*pow3(Mt) + 15571*
        Mt*pow3(Mst2)*pow3(s2t) + 89312*pow4(Mt) - 4199*pow4(Mst2)*pow4(s2t)) -
        63*pow4(Mst2)*(-72*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 8*Mst2*s2t*pow3(Mt)
        + 2*Mt*pow3(Mst2)*pow3(s2t) + 208*pow4(Mt) + 5*pow4(Mst2)*pow4(s2t)) -
        24*pow2(Mst1)*pow2(Mst2)*(-1791*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 4307*
        Mst2*s2t*pow3(Mt) + 767*Mt*pow3(Mst2)*pow3(s2t) + 9406*pow4(Mt) + 20*
        pow4(Mst2)*pow4(s2t))))*pow8(Msq) + 3*pow2(Mgl)*(20*Dmsqst2*pow2(Mst2)*
        pow2(s2t)*(pow2(Mst1)*pow3(Dmsqst2)*(108*pow2(Mst2)*pow2(Mt) + 4*pow2(
        Mst1)*(37*pow2(Mt) - 7*pow2(Mst2)*pow2(s2t)) - 9*pow2(s2t)*pow4(Mst2))
        + pow2(Dmsqst2)*pow2(Msq)*(2*(44*pow2(Mt) - 7*pow2(Mst2)*pow2(s2t))*
        pow4(Mst1) + 9*(-4*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow4(Mst2) + 3*
        pow2(Mst1)*(16*pow2(Mst2)*pow2(Mt) + 5*pow2(s2t)*pow4(Mst2))) +
        Dmsqst2*pow4(Msq)*(28*pow2(Mt)*pow4(Mst1) + 18*(-4*pow2(Mt) + pow2(
        Mst2)*pow2(s2t))*pow4(Mst2) + 3*pow2(Mst1)*(-4*pow2(Mst2)*pow2(Mt) +
        13*pow2(s2t)*pow4(Mst2))) + pow6(Msq)*((-32*pow2(Mt) + 14*pow2(Mst2)*
        pow2(s2t))*pow4(Mst1) + 9*pow2(Mst1)*(-8*pow2(Mst2)*pow2(Mt) + 7*pow2(
        s2t)*pow4(Mst2)) + 27*(-4*pow2(Mt)*pow4(Mst2) + pow2(s2t)*pow6(Mst2))))
        + (54*pow4(Mst2)*(-106*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 56*(Mt + Mst2*
        s2t)*pow3(Mt) - 14*Mt*pow3(Mst2)*pow3(s2t) + 23*pow4(Mst2)*pow4(s2t)) +
        6*pow2(Mst1)*pow2(Mst2)*(-2516*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 3016*
        Mst2*s2t*pow3(Mt) - 292*Mt*pow3(Mst2)*pow3(s2t) + 1672*pow4(Mt) + 643*
        pow4(Mst2)*pow4(s2t)) + pow4(Mst1)*(-27428*pow2(Mst2)*pow2(Mt)*pow2(
        s2t) + 43840*Mst2*s2t*pow3(Mt) - 1756*Mt*pow3(Mst2)*pow3(s2t) + 18688*
        pow4(Mt) + 1795*pow4(Mst2)*pow4(s2t)))*pow8(Msq))))/(1458.*pow2(Mgl)*
        pow4(Mst2)*pow8(Msq)) - (pow2(z2)*(pow4(s2t)*(-18*pow4(Mst2)*(9 + (474*
        Dmglst2)/Mgl + (683*pow2(Dmglst2))/pow2(Mgl) - (90*Dmsqst2)/pow2(Msq) -
        (60*pow2(Dmsqst2))/pow4(Msq) - (30*pow3(Dmsqst2))/pow6(Msq)) + pow4(
        Mst1)*(5385 - (16796*Dmglst2*(Dmglst2 + Mgl))/pow2(Mgl) + (840*Dmsqst2)
        /pow2(Msq) - (840*pow3(Dmsqst2))/pow6(Msq) - (1680*pow4(Dmsqst2))/pow8(
        Msq)) + 6*pow2(Mst1)*pow2(Mst2)*(2577 + (1784*pow2(Dmglst2))/pow2(Mgl)
        + (630*Dmsqst2)/pow2(Msq) + (390*pow2(Dmsqst2))/pow4(Msq) + (150*pow3(
        Dmsqst2))/pow6(Msq) - (90*pow4(Dmsqst2))/pow8(Msq))) + (4*Mt*(4*pow3(
        Mt)*(2*pow2(Dmglst2)*(-14109*pow2(Mst1)*pow2(Mst2) + 11164*pow4(Mst1) -
        819*pow4(Mst2)) + 3*pow2(Mgl)*(627*pow2(Mst1)*pow2(Mst2) + 1168*pow4(
        Mst1) + 189*pow4(Mst2)) + 4*Dmglst2*Mgl*(1449*pow2(Mst1)*pow2(Mst2) +
        5582*pow4(Mst1) + 189*pow4(Mst2))) - 4*Mst2*s2t*pow2(Mt)*(pow2(Dmglst2)
        *(-12921*pow2(Mst1)*pow2(Mst2) + 51460*pow4(Mst1) - 63*pow4(Mst2)) - 3*
        pow2(Mgl)*(1131*pow2(Mst1)*pow2(Mst2) + 2740*pow4(Mst1) + 189*pow4(
        Mst2)) + Dmglst2*Mgl*(12447*pow2(Mst1)*pow2(Mst2) + 51460*pow4(Mst1) +
        945*pow4(Mst2))) - pow3(Mst2)*pow3(s2t)*(pow2(Dmglst2)*(9204*pow2(Mst1)
        *pow2(Mst2) - 15571*pow4(Mst1) - 93249*pow4(Mst2)) + 3*pow2(Mgl)*(438*
        pow2(Mst1)*pow2(Mst2) + 439*pow4(Mst1) - 4995*pow4(Mst2)) - Dmglst2*
        Mgl*(7218*pow2(Mst1)*pow2(Mst2) + 15571*pow4(Mst1) + 47601*pow4(Mst2)))
        + (3*Mt*pow2(Mst2)*pow2(s2t)*(20*Dmsqst2*pow2(Mgl)*(pow3(Dmsqst2)*(27*
        pow2(Mst1)*pow2(Mst2) + 37*pow4(Mst1)) + Dmsqst2*pow4(Msq)*(-3*pow2(
        Mst1)*pow2(Mst2) + 7*pow4(Mst1) - 18*pow4(Mst2)) + pow2(Dmsqst2)*pow2(
        Msq)*(12*pow2(Mst1)*pow2(Mst2) + 22*pow4(Mst1) - 9*pow4(Mst2)) - (18*
        pow2(Mst1)*pow2(Mst2) + 8*pow4(Mst1) + 27*pow4(Mst2))*pow6(Msq)) + (24*
        Dmglst2*Mgl*pow2(Mst1)*(769*pow2(Mst1) + 150*pow2(Mst2)) + 12*pow2(
        Dmglst2)*(597*pow2(Mst1)*pow2(Mst2) + 1538*pow4(Mst1) + 63*pow4(Mst2))
        - pow2(Mgl)*(3774*pow2(Mst1)*pow2(Mst2) + 6857*pow4(Mst1) + 1431*pow4(
        Mst2)))*pow8(Msq)))/pow8(Msq)))/(pow2(Mgl)*pow4(Mst2))))/1944. - (z4*(
        pow4(s2t)*(-36*pow4(Mst2)*((318*Dmglst2)/Mgl + (1111*pow2(Dmglst2))/
        pow2(Mgl) - (15*(2*pow2(Dmsqst2)*pow2(Msq) + pow3(Dmsqst2) + 3*Dmsqst2*
        pow4(Msq) + 6*pow6(Msq)))/pow6(Msq)) + 6*pow2(Mst1)*pow2(Mst2)*((486*
        Dmglst2)/Mgl + (2513*pow2(Dmglst2))/pow2(Mgl) + 6*(389 + (105*Dmsqst2)/
        pow2(Msq) + (65*pow2(Dmsqst2))/pow4(Msq) + (25*pow3(Dmsqst2))/pow6(Msq)
        - (15*pow4(Dmsqst2))/pow8(Msq))) + pow4(Mst1)*(3441 - (16796*Dmglst2*(
        Dmglst2 + Mgl))/pow2(Mgl) + (840*Dmsqst2)/pow2(Msq) - (840*pow3(
        Dmsqst2))/pow6(Msq) - (1680*pow4(Dmsqst2))/pow8(Msq))) + (4*Mt*(pow3(
        Mst2)*pow3(s2t)*(pow2(Dmglst2)*(-14550*pow2(Mst1)*pow2(Mst2) + 15571*
        pow4(Mst1) - 11241*pow4(Mst2)) + Dmglst2*Mgl*(1872*pow2(Mst1)*pow2(
        Mst2) + 15571*pow4(Mst1) - 7317*pow4(Mst2)) - 3*pow2(Mgl)*(2220*pow2(
        Mst1)*pow2(Mst2) + 439*pow4(Mst1) + 2295*pow4(Mst2))) + 4*pow2(Mt)*(
        pow2(Dmglst2)*(3*(-9406*Mt + 4307*Mst2*s2t)*pow2(Mst1)*pow2(Mst2) + 4*(
        5582*Mt - 12865*Mst2*s2t)*pow4(Mst1) + 117*(-14*Mt + 67*Mst2*s2t)*pow4(
        Mst2)) + 3*pow2(Mgl)*(3*(209*Mt + 377*Mst2*s2t)*pow2(Mst1)*pow2(Mst2) +
        4*(292*Mt + 685*Mst2*s2t)*pow4(Mst1) + 27*(7*Mt + 103*Mst2*s2t)*pow4(
        Mst2)) + Dmglst2*Mgl*(-9*(-644*Mt + 1383*Mst2*s2t)*pow2(Mst1)*pow2(
        Mst2) + 4*(5582*Mt - 12865*Mst2*s2t)*pow4(Mst1) + 27*(28*Mt + 253*Mst2*
        s2t)*pow4(Mst2))) + (3*Mt*pow2(Mst2)*pow2(s2t)*(20*Dmsqst2*pow2(Mgl)*(
        pow3(Dmsqst2)*(27*pow2(Mst1)*pow2(Mst2) + 37*pow4(Mst1)) + Dmsqst2*
        pow4(Msq)*(-3*pow2(Mst1)*pow2(Mst2) + 7*pow4(Mst1) - 18*pow4(Mst2)) +
        pow2(Dmsqst2)*pow2(Msq)*(12*pow2(Mst1)*pow2(Mst2) + 22*pow4(Mst1) - 9*
        pow4(Mst2)) - (18*pow2(Mst1)*pow2(Mst2) + 8*pow4(Mst1) + 27*pow4(Mst2))
        *pow6(Msq)) + (12*pow2(Dmglst2)*(597*pow2(Mst1)*pow2(Mst2) + 1538*pow4(
        Mst1) - 423*pow4(Mst2)) + 24*Dmglst2*Mgl*(150*pow2(Mst1)*pow2(Mst2) +
        769*pow4(Mst1) - 162*pow4(Mst2)) - pow2(Mgl)*(2478*pow2(Mst1)*pow2(
        Mst2) + 6857*pow4(Mst1) + 2727*pow4(Mst2)))*pow8(Msq)))/pow8(Msq)))/(
        pow2(Mgl)*pow4(Mst2))))/2916. + pow4(Mt)*(480.98395061728394 + (2212*
        log(pow2(Mst2)/pow2(Mst1)))/81. + (3872*pow2(log(pow2(Mst2)/pow2(Mst1))
        ))/27. + (40*Dmsqst2*(Dmglst2*(Dmglst2 + Mgl)*(2954 - 1713*log(pow2(
        Mst2)/pow2(Mst1))) + 3*pow2(Mgl)*(-103 + 138*log(pow2(Mst2)/pow2(Mst1))
        + 9*pow2(log(pow2(Mst2)/pow2(Mst1))))))/(81.*pow2(Mgl)*pow2(Msq)) - (8*
        pow3(log(pow2(Mst2)/pow2(Mst1))))/9. + (4*Dmglst2*(773533 + 229915*log(
        pow2(Mst2)/pow2(Mst1)) - 199920*pow2(log(pow2(Mst2)/pow2(Mst1))) +
        40320*pow3(log(pow2(Mst2)/pow2(Mst1)))))/(2835.*Mgl) + (pow2(Dmglst2)*(
        13292051 + 606612*log(pow2(Mst2)/pow2(Mst1)) - 223920*pow2(log(pow2(
        Mst2)/pow2(Mst1))) + 518400*pow3(log(pow2(Mst2)/pow2(Mst1)))))/(6075.*
        pow2(Mgl)) - (5*pow2(Dmsqst2)*(5567 - 5502*log(pow2(Mst2)/pow2(Mst1)) +
        (8*Dmglst2*(Dmglst2 + Mgl)*(-2954 + 1713*log(pow2(Mst2)/pow2(Mst1))))/
        pow2(Mgl) - 108*pow2(log(pow2(Mst2)/pow2(Mst1)))))/(81.*pow4(Msq)) - (
        2*S2*(2*Dmglst2*Mgl*(119025*pow2(Mst1)*pow2(Mst2) + 2001242*pow4(Mst1)
        - 129033*pow4(Mst2)) + 21*pow2(Mgl)*(37467*pow2(Mst1)*pow2(Mst2) +
        69224*pow4(Mst1) + 8829*pow4(Mst2)) + pow2(Dmglst2)*(386391*pow2(Mst1)*
        pow2(Mst2) + 4002484*pow4(Mst1) + 731241*pow4(Mst2)) + 210*log(pow2(
        Mst2)/pow2(Mst1))*(2*pow2(Dmglst2)*(-14109*pow2(Mst1)*pow2(Mst2) +
        11164*pow4(Mst1) - 819*pow4(Mst2)) + 3*pow2(Mgl)*(627*pow2(Mst1)*pow2(
        Mst2) + 1168*pow4(Mst1) + 189*pow4(Mst2)) + 4*Dmglst2*Mgl*(1449*pow2(
        Mst1)*pow2(Mst2) + 5582*pow4(Mst1) + 189*pow4(Mst2)))))/(2835.*pow2(
        Mgl)*pow4(Mst2)) - (10*(4331 - 3846*log(pow2(Mst2)/pow2(Mst1)) + (4*
        Dmglst2*(Dmglst2 + Mgl)*(-2954 + 1713*log(pow2(Mst2)/pow2(Mst1))))/
        pow2(Mgl))*pow3(Dmsqst2))/(81.*pow6(Msq)) + (pow2(Mst1)*(
        602.3877366255144 + (89408*log(pow2(Mst2)/pow2(Mst1)))/405. + (4544*
        pow2(log(pow2(Mst2)/pow2(Mst1))))/27. - (4*Dmsqst2*(91321 - 92010*log(
        pow2(Mst2)/pow2(Mst1)) + 4050*pow2(log(pow2(Mst2)/pow2(Mst1))) - (2250*
        Dmglst2*(Dmglst2 + Mgl)*(401 - 274*log(pow2(Mst2)/pow2(Mst1)) + 4*pow2(
        log(pow2(Mst2)/pow2(Mst1)))))/pow2(Mgl)))/(675.*pow2(Msq)) + (440*pow3(
        log(pow2(Mst2)/pow2(Mst1))))/27. + (pow2(Dmglst2)*(-55610713 -
        374697540*log(pow2(Mst2)/pow2(Mst1)) + 821715300*pow2(log(pow2(Mst2)/
        pow2(Mst1))) + 10395000*pow3(log(pow2(Mst2)/pow2(Mst1)))))/(637875.*
        pow2(Mgl)) + (2*Dmglst2*(98884301 - 10812060*log(pow2(Mst2)/pow2(Mst1))
        - 112511700*pow2(log(pow2(Mst2)/pow2(Mst1))) + 20727000*pow3(log(pow2(
        Mst2)/pow2(Mst1)))))/(212625.*Mgl) - (2*pow2(Dmsqst2)*(23941 - 20010*
        log(pow2(Mst2)/pow2(Mst1)) + 300*pow2(log(pow2(Mst2)/pow2(Mst1))) - (
        300*Dmglst2*(Dmglst2 + Mgl)*(401 - 274*log(pow2(Mst2)/pow2(Mst1)) + 4*
        pow2(log(pow2(Mst2)/pow2(Mst1)))))/pow2(Mgl)))/(45.*pow4(Msq)) - (8*(
        133897 - 104070*log(pow2(Mst2)/pow2(Mst1)) + 225*pow2(log(pow2(Mst2)/
        pow2(Mst1))) - (1125*Dmglst2*(Dmglst2 + Mgl)*(401 - 274*log(pow2(Mst2)/
        pow2(Mst1)) + 4*pow2(log(pow2(Mst2)/pow2(Mst1)))))/pow2(Mgl))*pow3(
        Dmsqst2))/(675.*pow6(Msq)) - ((2109.8103703703705 - (70988*log(pow2(
        Mst2)/pow2(Mst1)))/45. - 8*pow2(log(pow2(Mst2)/pow2(Mst1))) - (40*
        Dmglst2*(Dmglst2 + Mgl)*(401 - 274*log(pow2(Mst2)/pow2(Mst1)) + 4*pow2(
        log(pow2(Mst2)/pow2(Mst1)))))/(3.*pow2(Mgl)))*pow4(Dmsqst2))/pow8(Msq))
        )/pow2(Mst2) + (pow4(Mst1)*(562637865189 + 269023439160*log(pow2(Mst2)/
        pow2(Mst1)) + 27258298200*pow2(log(pow2(Mst2)/pow2(Mst1))) - (20*
        Dmglst2*(Dmglst2 + Mgl)*(54946675289 + 53298152460*log(pow2(Mst2)/pow2(
        Mst1)) + 32279083200*pow2(log(pow2(Mst2)/pow2(Mst1))) - 5156524800*
        pow3(log(pow2(Mst2)/pow2(Mst1)))))/pow2(Mgl) + 22226400000*pow3(log(
        pow2(Mst2)/pow2(Mst1))) - (25920*Dmsqst2*(13165568 - 18387705*log(pow2(
        Mst2)/pow2(Mst1)) + 1885275*pow2(log(pow2(Mst2)/pow2(Mst1))))*(pow2(
        Dmsqst2)*pow2(Msq) + pow3(Dmsqst2) + Dmsqst2*pow4(Msq) + pow6(Msq)))/
        pow8(Msq)))/(3.750705e8*pow4(Mst2)) - (5*(11757 - 9882*log(pow2(Mst2)/
        pow2(Mst1)) + (8*Dmglst2*(Dmglst2 + Mgl)*(-2954 + 1713*log(pow2(Mst2)/
        pow2(Mst1))))/pow2(Mgl) + 108*pow2(log(pow2(Mst2)/pow2(Mst1))))*pow4(
        Dmsqst2))/(81.*pow8(Msq)) + ((-64*(-4*pow2(Dmglst2) + pow2(Mgl) - 2*
        log(pow2(Mst2)/pow2(Mst1))*pow2(Dmglst2 + Mgl) + (4*Dmglst2*Mgl + 10*
        pow2(Dmglst2) + pow2(Mgl))*pow2(log(pow2(Mst2)/pow2(Mst1))))*pow4(Mst2)
        )/(9.*pow4(Mst1)) + (16*OepS2*(2*pow2(Dmglst2)*(-14109*pow2(Mst1)*pow2(
        Mst2) + 11164*pow4(Mst1) - 819*pow4(Mst2)) + 3*pow2(Mgl)*(627*pow2(
        Mst1)*pow2(Mst2) + 1168*pow4(Mst1) + 189*pow4(Mst2)) + 4*Dmglst2*Mgl*(
        1449*pow2(Mst1)*pow2(Mst2) + 5582*pow4(Mst1) + 189*pow4(Mst2))))/(2187.
        *pow4(Mst2)) + (4*pow2(Mst2)*(4*pow2(Dmglst2)*(90*pow2(Msq)*pow3(
        Dmsqst2) + 90*pow4(Dmsqst2) + 90*pow2(Dmsqst2)*pow4(Msq) + 90*Dmsqst2*
        pow6(Msq) - 197*pow8(Msq)) + 24*Dmglst2*Mgl*(15*pow2(Msq)*pow3(Dmsqst2)
        + 15*pow4(Dmsqst2) + 15*pow2(Dmsqst2)*pow4(Msq) + 15*Dmsqst2*pow6(Msq)
        - 7*pow8(Msq)) + 59*(2*Dmglst2*Mgl + 3*pow2(Dmglst2) + pow2(Mgl))*pow2(
        log(pow2(Mst2)/pow2(Mst1)))*pow8(Msq) + pow2(Mgl)*(-90*pow2(Msq)*pow3(
        Dmsqst2) - 150*pow4(Dmsqst2) - 30*pow2(Dmsqst2)*pow4(Msq) + 30*Dmsqst2*
        pow6(Msq) + 71*pow8(Msq)) + 2*log(pow2(Mst2)/pow2(Mst1))*(30*pow2(Mgl)*
        pow2(Msq)*pow3(Dmsqst2) + 30*pow2(Mgl)*pow4(Dmsqst2) + 30*pow2(Dmsqst2)
        *pow2(Mgl)*pow4(Msq) + 30*Dmsqst2*pow2(Mgl)*pow6(Msq) + (-70*Dmglst2*
        Mgl + 13*pow2(Dmglst2) - 61*pow2(Mgl))*pow8(Msq))))/(9.*pow2(Mst1)*
        pow8(Msq)))/pow2(Mgl)) + Mt*pow3(s2t)*(-(S2*(Dmglst2*Mgl*(252594*pow2(
        Mst1)*pow2(Mst2) + 836021*pow4(Mst1) - 36693*pow4(Mst2)) - 15*pow2(Mgl)
        *(6246*pow2(Mst1)*pow2(Mst2) + 6121*pow4(Mst1) - 81*pow4(Mst2)) + pow2(
        Dmglst2)*(49920*pow2(Mst1)*pow2(Mst2) + 836021*pow4(Mst1) + 37449*pow4(
        Mst2)) + 30*log(pow2(Mst2)/pow2(Mst1))*(pow2(Dmglst2)*(-9204*pow2(Mst1)
        *pow2(Mst2) + 15571*pow4(Mst1) - 63*pow4(Mst2)) - 3*pow2(Mgl)*(438*
        pow2(Mst1)*pow2(Mst2) + 439*pow4(Mst1) + 189*pow4(Mst2)) + Dmglst2*Mgl*
        (7218*pow2(Mst1)*pow2(Mst2) + 15571*pow4(Mst1) + 945*pow4(Mst2)))))/(
        810.*Mst2*pow2(Mgl)) + (pow4(Mst1)*(94.64797668038409 + (46*log(pow2(
        Mst2)/pow2(Mst1)))/27. + (20*Dmsqst2*(Dmglst2*(Dmglst2 + Mgl)*(40 - 3*
        log(pow2(Mst2)/pow2(Mst1))) + (-4 + log(pow2(Mst2)/pow2(Mst1)))*pow2(
        Mgl)))/(3.*pow2(Mgl)*pow2(Msq)) + (1000*pow2(log(pow2(Mst2)/pow2(Mst1))
        ))/27. - (800*pow3(log(pow2(Mst2)/pow2(Mst1))))/27. - (Dmglst2*(Dmglst2
        + Mgl)*(31897243 + 6714720*log(pow2(Mst2)/pow2(Mst1)) + 3188160*pow2(
        log(pow2(Mst2)/pow2(Mst1))) + 2177280*pow3(log(pow2(Mst2)/pow2(Mst1))))
        )/(87480.*pow2(Mgl)) - (pow2(Dmsqst2)*(505 + 565*log(pow2(Mst2)/pow2(
        Mst1)) + (40*Dmglst2*(Dmglst2 + Mgl)*(-40 + 3*log(pow2(Mst2)/pow2(Mst1)
        )))/pow2(Mgl)))/(6.*pow4(Msq)) - (5*(85 + 117*log(pow2(Mst2)/pow2(Mst1)
        ) + (4*Dmglst2*(Dmglst2 + Mgl)*(-40 + 3*log(pow2(Mst2)/pow2(Mst1))))/
        pow2(Mgl))*pow3(Dmsqst2))/(3.*pow6(Msq)) - (5*(239 + 355*log(pow2(Mst2)
        /pow2(Mst1)) + (8*Dmglst2*(Dmglst2 + Mgl)*(-40 + 3*log(pow2(Mst2)/pow2(
        Mst1))))/pow2(Mgl))*pow4(Dmsqst2))/(6.*pow8(Msq))))/Mst2 + Mst2*pow2(
        Mst1)*(52.51131687242798 + (68*B4)/9. + (2*DN)/9. + (442*log(pow2(Mst2)
        /pow2(Mst1)))/9. + (488*pow2(log(pow2(Mst2)/pow2(Mst1))))/9. + (pow2(
        Dmglst2)*(192.65089163237312 + (68*B4)/9. + (2*DN)/9. - (6898*log(pow2(
        Mst2)/pow2(Mst1)))/81. - (536*pow2(log(pow2(Mst2)/pow2(Mst1))))/9. - (
        260*pow3(log(pow2(Mst2)/pow2(Mst1))))/9.))/pow2(Mgl) + (Dmglst2*(
        5.598971193415638 + (68*B4)/9. + (2*DN)/9. + (1006*log(pow2(Mst2)/pow2(
        Mst1)))/27. - (40*pow2(log(pow2(Mst2)/pow2(Mst1))))/3. - (260*pow3(log(
        pow2(Mst2)/pow2(Mst1))))/9.))/Mgl - (260*pow3(log(pow2(Mst2)/pow2(Mst1)
        )))/9. + (80*Dmsqst2*(-(Dmglst2*Mgl*(-12 + 5*log(pow2(Mst2)/pow2(Mst1))
        )*(Dmsqst2 + pow2(Msq))) - (-12 + 5*log(pow2(Mst2)/pow2(Mst1)))*pow2(
        Dmglst2)*(Dmsqst2 + pow2(Msq)) + (-1 + log(pow2(Mst2)/pow2(Mst1)))*
        pow2(Mgl)*(2*Dmsqst2 + pow2(Msq))))/(3.*pow2(Mgl)*pow4(Msq)) - (80*(3 -
        3*log(pow2(Mst2)/pow2(Mst1)) + (Dmglst2*(Dmglst2 + Mgl)*(-12 + 5*log(
        pow2(Mst2)/pow2(Mst1))))/pow2(Mgl))*pow3(Dmsqst2))/(3.*pow6(Msq)) + (
        80*(Dmglst2*(Dmglst2 + Mgl)*(12 - 5*log(pow2(Mst2)/pow2(Mst1))) + 4*(-1
        + log(pow2(Mst2)/pow2(Mst1)))*pow2(Mgl))*pow4(Dmsqst2))/(3.*pow2(Mgl)*
        pow8(Msq))) + ((4*OepS2*(pow2(Dmglst2)*(-9204*pow2(Mst1)*pow2(Mst2) +
        15571*pow4(Mst1) - 63*pow4(Mst2)) - 3*pow2(Mgl)*(438*pow2(Mst1)*pow2(
        Mst2) + 439*pow4(Mst1) + 189*pow4(Mst2)) + Dmglst2*Mgl*(7218*pow2(Mst1)
        *pow2(Mst2) + 15571*pow4(Mst1) + 945*pow4(Mst2))))/(2187.*Mst2) + (64*(
        -(Dmglst2*Mgl) - 4*pow2(Dmglst2) + 2*log(pow2(Mst2)/pow2(Mst1))*(-(
        Dmglst2*Mgl) + pow2(Dmglst2) - pow2(Mgl)) + pow2(Mgl) + (3*Dmglst2*Mgl
        + 6*pow2(Dmglst2) + pow2(Mgl))*pow2(log(pow2(Mst2)/pow2(Mst1))))*pow5(
        Mst2))/(9.*pow2(Mst1)) - (pow3(Mst2)*(3*Dmglst2*Mgl*(86400*Dmsqst2*(2 +
        5*log(pow2(Mst2)/pow2(Mst1)))*(pow2(Dmsqst2)*pow2(Msq) + pow3(Dmsqst2)
        + Dmsqst2*pow4(Msq) + pow6(Msq)) - (69997 + 188640*B4 - 3600*DN -
        590400*log(pow2(Mst2)/pow2(Mst1)) + 264960*pow2(log(pow2(Mst2)/pow2(
        Mst1))) - 185760*pow3(log(pow2(Mst2)/pow2(Mst1))))*pow8(Msq)) - 135*
        pow2(Mgl)*(7680*pow2(Msq)*pow3(Dmsqst2) + 10560*pow4(Dmsqst2) + 4800*
        pow2(Dmsqst2)*pow4(Msq) + 1920*Dmsqst2*pow6(Msq) + 128*log(pow2(Mst2)/
        pow2(Mst1))*(45*pow2(Msq)*pow3(Dmsqst2) + 60*pow4(Dmsqst2) + 30*pow2(
        Dmsqst2)*pow4(Msq) + 15*Dmsqst2*pow6(Msq) - 91*pow8(Msq)) + (8537 +
        1760*B4 - 16*DN)*pow8(Msq) + 7552*pow2(log(pow2(Mst2)/pow2(Mst1)))*
        pow8(Msq) - 2080*pow3(log(pow2(Mst2)/pow2(Mst1)))*pow8(Msq)) + pow2(
        Dmglst2)*(259200*Dmsqst2*(2 + 5*log(pow2(Mst2)/pow2(Mst1)))*(pow2(
        Dmsqst2)*pow2(Msq) + pow3(Dmsqst2) + Dmsqst2*pow4(Msq) + pow6(Msq)) + (
        2048393 - 1058400*B4 + 23760*DN + 1139040*log(pow2(Mst2)/pow2(Mst1)) -
        129600*pow2(log(pow2(Mst2)/pow2(Mst1))) + 972000*pow3(log(pow2(Mst2)/
        pow2(Mst1))))*pow8(Msq))))/(9720.*pow8(Msq)))/pow2(Mgl)) + z2*(-(s2t*
        pow3(Mt)*((4*pow4(Mst1)*(1724277 - 492480*log(pow2(Mst2)/pow2(Mst1)) +
        (Dmglst2*(Dmglst2 + Mgl)*(-9908167 + 1769040*log(pow2(Mst2)/pow2(Mst1))
        ) - (48600*Dmsqst2*(Dmglst2*(Dmglst2 + Mgl)*(87 - 44*log(pow2(Mst2)/
        pow2(Mst1))) + (-5 + 4*log(pow2(Mst2)/pow2(Mst1)))*pow2(Mgl))*(pow2(
        Dmsqst2)*pow2(Msq) + pow3(Dmsqst2) + Dmsqst2*pow4(Msq) + pow6(Msq)))/
        pow8(Msq))/pow2(Mgl)))/(3645.*pow3(Mst2)) + Mst2*(508.837037037037 - (
        3032*log(pow2(Mst2)/pow2(Mst1)))/9. + (Dmglst2*(-4189 + 1448*log(pow2(
        Mst2)/pow2(Mst1))))/(9.*Mgl) + ((938.6416225749559 - (4472*log(pow2(
        Mst2)/pow2(Mst1)))/9.)*pow2(Dmglst2))/pow2(Mgl) - (160*Dmsqst2*(
        Dmglst2*(Dmglst2 + Mgl)*(19 - 12*log(pow2(Mst2)/pow2(Mst1))) + (-1 + 4*
        log(pow2(Mst2)/pow2(Mst1)))*pow2(Mgl)))/(3.*pow2(Mgl)*pow2(Msq)) + (80*
        pow2(Dmsqst2)*(7 - 12*log(pow2(Mst2)/pow2(Mst1)) + (2*Dmglst2*(Dmglst2
        + Mgl)*(-19 + 12*log(pow2(Mst2)/pow2(Mst1))))/pow2(Mgl)))/(3.*pow4(Msq)
        ) + (160*(6 - 8*log(pow2(Mst2)/pow2(Mst1)) + (Dmglst2*(Dmglst2 + Mgl)*(
        -19 + 12*log(pow2(Mst2)/pow2(Mst1))))/pow2(Mgl))*pow3(Dmsqst2))/(3.*
        pow6(Msq)) + (80*(17 - 20*log(pow2(Mst2)/pow2(Mst1)) + (2*Dmglst2*(
        Dmglst2 + Mgl)*(-19 + 12*log(pow2(Mst2)/pow2(Mst1))))/pow2(Mgl))*pow4(
        Dmsqst2))/(3.*pow8(Msq))) + (pow2(Mst1)*(1162.116049382716 - 440*log(
        pow2(Mst2)/pow2(Mst1)) + (Dmglst2*(-510149 + 118200*log(pow2(Mst2)/
        pow2(Mst1))))/(135.*Mgl) + ((5625.55543797766 - (14168*log(pow2(Mst2)/
        pow2(Mst1)))/9.)*pow2(Dmglst2))/pow2(Mgl) + (160*Dmsqst2*(-51*Dmglst2*
        Mgl*(Dmsqst2 + pow2(Msq)) - 51*pow2(Dmglst2)*(Dmsqst2 + pow2(Msq)) +
        pow2(Mgl)*(12*Dmsqst2 + 5*pow2(Msq)) + 4*log(pow2(Mst2)/pow2(Mst1))*(7*
        Dmglst2*Mgl*(Dmsqst2 + pow2(Msq)) + 7*pow2(Dmglst2)*(Dmsqst2 + pow2(
        Msq)) - pow2(Mgl)*(2*Dmsqst2 + pow2(Msq)))))/(3.*pow2(Mgl)*pow4(Msq)) +
        (160*(19 - 12*log(pow2(Mst2)/pow2(Mst1)) + (Dmglst2*(Dmglst2 + Mgl)*(-
        51 + 28*log(pow2(Mst2)/pow2(Mst1))))/pow2(Mgl))*pow3(Dmsqst2))/(3.*
        pow6(Msq)) - (160*(Dmglst2*(Dmglst2 + Mgl)*(51 - 28*log(pow2(Mst2)/
        pow2(Mst1))) + 2*(-13 + 8*log(pow2(Mst2)/pow2(Mst1)))*pow2(Mgl))*pow4(
        Dmsqst2))/(3.*pow2(Mgl)*pow8(Msq))))/Mst2)) + pow4(Mt)*(
        362.837037037037 - (596*log(pow2(Mst2)/pow2(Mst1)))/3. + (4*Dmglst2*(-
        209341 + 82320*log(pow2(Mst2)/pow2(Mst1))))/(945.*Mgl) + (2*(66191 -
        11952*log(pow2(Mst2)/pow2(Mst1)))*pow2(Dmglst2))/(81.*pow2(Mgl)) + (
        pow2(Mst1)*(765.5283950617284 - (Dmglst2*(4492.976366843033 - 1632*log(
        pow2(Mst2)/pow2(Mst1))))/Mgl - (3224*log(pow2(Mst2)/pow2(Mst1)))/9. + (
        (12375.965902410348 - (23984*log(pow2(Mst2)/pow2(Mst1)))/9.)*pow2(
        Dmglst2))/pow2(Mgl) + (160*Dmsqst2*(6 - 4*log(pow2(Mst2)/pow2(Mst1)) +
        (Dmglst2*(Dmglst2 + Mgl)*(-77 + 36*log(pow2(Mst2)/pow2(Mst1))))/pow2(
        Mgl)))/(3.*pow2(Msq)) + (40*pow2(Dmsqst2)*(69 - 36*log(pow2(Mst2)/pow2(
        Mst1)) + (4*Dmglst2*(Dmglst2 + Mgl)*(-77 + 36*log(pow2(Mst2)/pow2(Mst1)
        )))/pow2(Mgl)))/(3.*pow4(Msq)) + (80*(57 - 28*log(pow2(Mst2)/pow2(Mst1)
        ) + (2*Dmglst2*(Dmglst2 + Mgl)*(-77 + 36*log(pow2(Mst2)/pow2(Mst1))))/
        pow2(Mgl))*pow3(Dmsqst2))/(3.*pow6(Msq)) + (40*(159 - 76*log(pow2(Mst2)
        /pow2(Mst1)) + (4*Dmglst2*(Dmglst2 + Mgl)*(-77 + 36*log(pow2(Mst2)/
        pow2(Mst1))))/pow2(Mgl))*pow4(Dmsqst2))/(3.*pow8(Msq))))/pow2(Mst2) + (
        2*pow4(Mst1)*(13164249 - 8300880*log(pow2(Mst2)/pow2(Mst1)) + (8*
        Dmglst2*(Dmglst2 + Mgl)*(-19152737 + 6872040*log(pow2(Mst2)/pow2(Mst1))
        ))/pow2(Mgl) + (680400*Dmsqst2*(7 - 6*log(pow2(Mst2)/pow2(Mst1)))*(
        pow2(Dmsqst2)*pow2(Msq) + pow3(Dmsqst2) + Dmsqst2*pow4(Msq) + pow6(Msq)
        ))/pow8(Msq)))/(25515.*pow4(Mst2)) - (40*Dmsqst2*(-4*Dmglst2*Mgl*(-23 +
        10*log(pow2(Mst2)/pow2(Mst1)))*(pow2(Dmsqst2)*pow2(Msq) + pow3(Dmsqst2)
        + Dmsqst2*pow4(Msq) + pow6(Msq)) - 4*(-23 + 10*log(pow2(Mst2)/pow2(
        Mst1)))*pow2(Dmglst2)*(pow2(Dmsqst2)*pow2(Msq) + pow3(Dmsqst2) +
        Dmsqst2*pow4(Msq) + pow6(Msq)) + pow2(Mgl)*(-38*pow2(Dmsqst2)*pow2(Msq)
        - 51*pow3(Dmsqst2) - 25*Dmsqst2*pow4(Msq) - 12*pow6(Msq) + 2*log(pow2(
        Mst2)/pow2(Mst1))*(10*pow2(Dmsqst2)*pow2(Msq) + 13*pow3(Dmsqst2) + 7*
        Dmsqst2*pow4(Msq) + 4*pow6(Msq)))))/(3.*pow2(Mgl)*pow8(Msq)) - (8*pow2(
        Mst2)*(12*pow2(Dmglst2)*(5*Dmsqst2*(pow2(Dmsqst2)*pow2(Msq) + pow3(
        Dmsqst2) + Dmsqst2*pow4(Msq) + pow6(Msq)) - 2*pow8(Msq)) + 4*Dmglst2*
        Mgl*(15*Dmsqst2*(pow2(Dmsqst2)*pow2(Msq) + pow3(Dmsqst2) + Dmsqst2*
        pow4(Msq) + pow6(Msq)) + 7*pow8(Msq)) - pow2(Mgl)*(60*pow2(Msq)*pow3(
        Dmsqst2) + 75*pow4(Dmsqst2) + 45*pow2(Dmsqst2)*pow4(Msq) + 30*Dmsqst2*
        pow6(Msq) + 19*pow8(Msq))))/(9.*pow2(Mgl)*pow2(Mst1)*pow8(Msq))) -
        pow4(s2t)*((pow2(Mst1)*pow2(Mst2)*(362299 - 33300*log(pow2(Mst2)/pow2(
        Mst1)) - (24*Dmglst2*(5917 + 195*log(pow2(Mst2)/pow2(Mst1))))/Mgl - (4*
        (-77099 + 18225*log(pow2(Mst2)/pow2(Mst1)))*pow2(Dmglst2))/(3.*pow2(
        Mgl)) - (150*Dmsqst2*(-144*Dmglst2*(Dmglst2 + Mgl)*(-4 + log(pow2(Mst2)
        /pow2(Mst1))) + (-743 + 72*log(pow2(Mst2)/pow2(Mst1)))*pow2(Mgl)))/(
        pow2(Mgl)*pow2(Msq)) + (3240*pow2(Dmsqst2)*(26.929012345679013 - 5*log(
        pow2(Mst2)/pow2(Mst1)) + (20*Dmglst2*(Dmglst2 + Mgl)*(-4 + log(pow2(
        Mst2)/pow2(Mst1))))/(3.*pow2(Mgl))))/pow4(Msq) + (50*(1261 - 432*log(
        pow2(Mst2)/pow2(Mst1)) + (432*Dmglst2*(Dmglst2 + Mgl)*(-4 + log(pow2(
        Mst2)/pow2(Mst1))))/pow2(Mgl))*pow3(Dmsqst2))/pow6(Msq) + (150*(259 -
        180*log(pow2(Mst2)/pow2(Mst1)) + (144*Dmglst2*(Dmglst2 + Mgl)*(-4 +
        log(pow2(Mst2)/pow2(Mst1))))/pow2(Mgl))*pow4(Dmsqst2))/pow8(Msq)))/
        3240. + (pow4(Mst2)*(6053 + 1104*log(pow2(Mst2)/pow2(Mst1)) + (14*
        Dmglst2*(-281 + 240*log(pow2(Mst2)/pow2(Mst1))))/Mgl + ((-80291 +
        13320*log(pow2(Mst2)/pow2(Mst1)))*pow2(Dmglst2))/(15.*pow2(Mgl)) + (30*
        Dmsqst2*(77 + (24*Dmglst2*(Dmglst2 + Mgl)*(4 + 3*log(pow2(Mst2)/pow2(
        Mst1))))/pow2(Mgl)))/pow2(Msq) + (20*pow2(Dmsqst2)*(29 - 9*log(pow2(
        Mst2)/pow2(Mst1)) + (36*Dmglst2*(Dmglst2 + Mgl)*(4 + 3*log(pow2(Mst2)/
        pow2(Mst1))))/pow2(Mgl)))/pow4(Msq) - (10*(115 + 36*log(pow2(Mst2)/
        pow2(Mst1)) - (72*Dmglst2*(Dmglst2 + Mgl)*(4 + 3*log(pow2(Mst2)/pow2(
        Mst1))))/pow2(Mgl))*pow3(Dmsqst2))/pow6(Msq) - (180*(16 + 3*log(pow2(
        Mst2)/pow2(Mst1)) - (4*Dmglst2*(Dmglst2 + Mgl)*(4 + 3*log(pow2(Mst2)/
        pow2(Mst1))))/pow2(Mgl))*pow4(Dmsqst2))/pow8(Msq)))/216. + pow4(Mst1)*(
        38.546039094650205 + (5*log(pow2(Mst2)/pow2(Mst1)))/4. - (Dmglst2*(
        Dmglst2 + Mgl)*(155.89787379972566 - 8*log(pow2(Mst2)/pow2(Mst1))))/
        pow2(Mgl) + (Dmsqst2*(12.27366255144033 + (10*log(pow2(Mst2)/pow2(Mst1)
        ))/3. - (5*Dmglst2*(Dmglst2 + Mgl)*(7 + 10*log(pow2(Mst2)/pow2(Mst1))))
        /(3.*pow2(Mgl))))/pow2(Msq) + (5*pow2(Dmsqst2)*(17 + 14*log(pow2(Mst2)/
        pow2(Mst1)) - (4*Dmglst2*(Dmglst2 + Mgl)*(7 + 10*log(pow2(Mst2)/pow2(
        Mst1))))/pow2(Mgl)))/(12.*pow4(Msq)) + (5*(92 + 405*log(pow2(Mst2)/
        pow2(Mst1)) - (81*Dmglst2*(Dmglst2 + Mgl)*(7 + 10*log(pow2(Mst2)/pow2(
        Mst1))))/pow2(Mgl))*pow3(Dmsqst2))/(243.*pow6(Msq)) - ((
        3.2973251028806585 - (65*log(pow2(Mst2)/pow2(Mst1)))/6. + (5*Dmglst2*(
        Dmglst2 + Mgl)*(7 + 10*log(pow2(Mst2)/pow2(Mst1))))/(3.*pow2(Mgl)))*
        pow4(Dmsqst2))/pow8(Msq)) + (pow6(Mst2)*(12*pow2(Dmglst2)*(5*Dmsqst2*(
        pow2(Dmsqst2)*pow2(Msq) + pow3(Dmsqst2) + Dmsqst2*pow4(Msq) + pow6(Msq)
        ) - 2*pow8(Msq)) + 4*Dmglst2*Mgl*(15*Dmsqst2*(pow2(Dmsqst2)*pow2(Msq) +
        pow3(Dmsqst2) + Dmsqst2*pow4(Msq) + pow6(Msq)) + 7*pow8(Msq)) - pow2(
        Mgl)*(60*pow2(Msq)*pow3(Dmsqst2) + 75*pow4(Dmsqst2) + 45*pow2(Dmsqst2)*
        pow4(Msq) + 30*Dmsqst2*pow6(Msq) + 19*pow8(Msq))))/(18.*pow2(Mgl)*pow2(
        Mst1)*pow8(Msq))) + pow2(Mt)*pow2(s2t)*((pow4(Mst1)*(677.1664609053498
         + (130*log(pow2(Mst2)/pow2(Mst1)))/
        3. + (2*Dmglst2*(Dmglst2 + Mgl)*(-78533 + 180*log(pow2(Mst2)/pow2(Mst1)
        )))/(81.*pow2(Mgl)) + (Dmsqst2*(16.625514403292183 + (10*Dmglst2*(
        Dmglst2 + Mgl)*(-61 + 10*log(pow2(Mst2)/pow2(Mst1))))/(3.*pow2(Mgl))))/
        pow2(Msq) + (pow2(Dmsqst2)*(14.619341563786008 - 5*log(pow2(Mst2)/pow2(
        Mst1)) + (10*Dmglst2*(Dmglst2 + Mgl)*(-61 + 10*log(pow2(Mst2)/pow2(
        Mst1))))/(3.*pow2(Mgl))))/pow4(Msq) + ((12.613168724279836 - 10*log(
        pow2(Mst2)/pow2(Mst1)) + (10*Dmglst2*(Dmglst2 + Mgl)*(-61 + 10*log(
        pow2(Mst2)/pow2(Mst1))))/(3.*pow2(Mgl)))*pow3(Dmsqst2))/pow6(Msq) + ((
        10.606995884773662 - 15*log(pow2(Mst2)/pow2(Mst1)) + (10*Dmglst2*(
        Dmglst2 + Mgl)*(-61 + 10*log(pow2(Mst2)/pow2(Mst1))))/(3.*pow2(Mgl)))*
        pow4(Dmsqst2))/pow8(Msq)))/pow2(Mst2) + pow2(Mst2)*(10.975925925925926
         + (610*log(pow2(Mst2)/pow2(Mst1)))/
        9. + (4*Dmglst2*(-4849 + 8910*log(pow2(Mst2)/pow2(Mst1))))/(135.*Mgl) +
        ((-123227 + 230300*log(pow2(Mst2)/pow2(Mst1)))*pow2(Dmglst2))/(315.*
        pow2(Mgl)) + (5*Dmsqst2*(576*Dmglst2*Mgl*(Dmsqst2 + pow2(Msq)) + 576*
        pow2(Dmglst2)*(Dmsqst2 + pow2(Msq)) + pow2(Mgl)*(58*Dmsqst2 + 231*pow2(
        Msq))))/(27.*pow2(Mgl)*pow4(Msq)) - ((21.296296296296298 - (320*
        Dmglst2*(Dmglst2 + Mgl))/(3.*pow2(Mgl)))*pow3(Dmsqst2))/pow6(Msq) - (
        160*(-2*Dmglst2*Mgl - 2*pow2(Dmglst2) + pow2(Mgl))*pow4(Dmsqst2))/(3.*
        pow2(Mgl)*pow8(Msq))) + (4*pow4(Mst2)*(12*pow2(Dmglst2)*(5*Dmsqst2*(
        pow2(Dmsqst2)*pow2(Msq) + pow3(Dmsqst2) + Dmsqst2*pow4(Msq) + pow6(Msq)
        ) - 2*pow8(Msq)) + 4*Dmglst2*Mgl*(15*Dmsqst2*(pow2(Dmsqst2)*pow2(Msq) +
        pow3(Dmsqst2) + Dmsqst2*pow4(Msq) + pow6(Msq)) + 7*pow8(Msq)) - pow2(
        Mgl)*(60*pow2(Msq)*pow3(Dmsqst2) + 75*pow4(Dmsqst2) + 45*pow2(Dmsqst2)*
        pow4(Msq) + 30*Dmsqst2*pow6(Msq) + 19*pow8(Msq))))/(9.*pow2(Mgl)*pow2(
        Mst1)*pow8(Msq)) + (pow2(Mst1)*(5040*Dmglst2*Mgl*(120*Dmsqst2*(pow2(
        Dmsqst2)*pow2(Msq) + pow3(Dmsqst2) + Dmsqst2*pow4(Msq) + pow6(Msq)) + (
        -505 + 124*log(pow2(Mst2)/pow2(Mst1)))*pow8(Msq)) + 30*pow2(Dmglst2)*(
        20160*Dmsqst2*(pow2(Dmsqst2)*pow2(Msq) + pow3(Dmsqst2) + Dmsqst2*pow4(
        Msq) + pow6(Msq)) + 37*(-7339 + 3192*log(pow2(Mst2)/pow2(Mst1)))*pow8(
        Msq)) - 7*pow2(Mgl)*(41800*pow2(Msq)*pow3(Dmsqst2) + 77850*pow4(
        Dmsqst2) + 5750*pow2(Dmsqst2)*pow4(Msq) - 30300*Dmsqst2*pow6(Msq) - (
        309749 + 12240*log(pow2(Mst2)/pow2(Mst1)))*pow8(Msq))))/(5670.*pow2(
        Mgl)*pow8(Msq))) + Mt*pow3(s2t)*((pow4(Mst1)*(21.57181069958848 - 34*
        log(pow2(Mst2)/pow2(Mst1)) + ((Dmglst2*(Dmglst2 + Mgl)*(-8583283 +
        2329560*log(pow2(Mst2)/pow2(Mst1))))/14580. + (40*Dmsqst2*(-11*Dmglst2*
        Mgl - 11*pow2(Dmglst2) + pow2(Mgl)))/(3.*pow2(Msq)) - (5*pow2(Dmsqst2)*
        (88*Dmglst2*(Dmglst2 + Mgl) - (53 + 30*log(pow2(Mst2)/pow2(Mst1)))*
        pow2(Mgl)))/(3.*pow4(Msq)) - (10*(44*Dmglst2*(Dmglst2 + Mgl) - (49 +
        30*log(pow2(Mst2)/pow2(Mst1)))*pow2(Mgl))*pow3(Dmsqst2))/(3.*pow6(Msq))
        - (5*(88*Dmglst2*(Dmglst2 + Mgl) - (143 + 90*log(pow2(Mst2)/pow2(Mst1))
        )*pow2(Mgl))*pow4(Dmsqst2))/(3.*pow8(Msq)))/pow2(Mgl)))/Mst2 + Mst2*
        pow2(Mst1)*(28.84567901234568 - 52*log(pow2(Mst2)/pow2(Mst1)) + (
        Dmglst2*(-206741 + 89640*log(pow2(Mst2)/pow2(Mst1))))/(810.*Mgl) + ((
        143005 + 26352*log(pow2(Mst2)/pow2(Mst1)))*pow2(Dmglst2))/(243.*pow2(
        Mgl)) - (80*Dmsqst2*(Dmglst2*(Dmglst2 + Mgl)*(5 - 3*log(pow2(Mst2)/
        pow2(Mst1))) + (1 + log(pow2(Mst2)/pow2(Mst1)))*pow2(Mgl)))/(3.*pow2(
        Mgl)*pow2(Msq)) - (40*pow2(Dmsqst2)*(1 + 3*log(pow2(Mst2)/pow2(Mst1)) +
        (2*Dmglst2*(Dmglst2 + Mgl)*(5 - 3*log(pow2(Mst2)/pow2(Mst1))))/pow2(
        Mgl)))/(3.*pow4(Msq)) - (80*(2*log(pow2(Mst2)/pow2(Mst1)) + (Dmglst2*(
        Dmglst2 + Mgl)*(5 - 3*log(pow2(Mst2)/pow2(Mst1))))/pow2(Mgl))*pow3(
        Dmsqst2))/(3.*pow6(Msq)) + (40*(1 - 5*log(pow2(Mst2)/pow2(Mst1)) + (2*
        Dmglst2*(Dmglst2 + Mgl)*(-5 + 3*log(pow2(Mst2)/pow2(Mst1))))/pow2(Mgl))
        *pow4(Dmsqst2))/(3.*pow8(Msq))) + (pow3(Mst2)*(3*Dmglst2*Mgl*(14400*
        Dmsqst2*(2 + 3*log(pow2(Mst2)/pow2(Mst1)))*(pow2(Dmsqst2)*pow2(Msq) +
        pow3(Dmsqst2) + Dmsqst2*pow4(Msq) + pow6(Msq)) + (-173947 + 68760*log(
        pow2(Mst2)/pow2(Mst1)))*pow8(Msq)) + pow2(Dmglst2)*(43200*Dmsqst2*(2 +
        3*log(pow2(Mst2)/pow2(Mst1)))*(pow2(Dmsqst2)*pow2(Msq) + pow3(Dmsqst2)
        + Dmsqst2*pow4(Msq) + pow6(Msq)) + (-897737 + 238680*log(pow2(Mst2)/
        pow2(Mst1)))*pow8(Msq)) - 45*pow2(Mgl)*(4800*pow2(Msq)*pow3(Dmsqst2) +
        6240*pow4(Dmsqst2) + 3360*pow2(Dmsqst2)*pow4(Msq) + 1920*Dmsqst2*pow6(
        Msq) + 2269*pow8(Msq) + 8*log(pow2(Mst2)/pow2(Mst1))*(240*pow2(Msq)*
        pow3(Dmsqst2) + 300*pow4(Dmsqst2) + 180*pow2(Dmsqst2)*pow4(Msq) + 120*
        Dmsqst2*pow6(Msq) + 7*pow8(Msq)))))/(1620.*pow2(Mgl)*pow8(Msq)))) -
        s2t*pow3(Mt)*((2*S2*((3*(-329511*Dmglst2*Mgl - 760703*pow2(Dmglst2) +
        275121*pow2(Mgl) + 210*log(pow2(Mst2)/pow2(Mst1))*(-4149*Dmglst2*Mgl +
        4307*pow2(Dmglst2) + 1131*pow2(Mgl)))*pow2(Mst1)*pow2(Mst2))/pow2(Mgl)
        + 56*(48891 + 30825*log(pow2(Mst2)/pow2(Mst1)) - (Dmglst2*(Dmglst2 +
        Mgl)*(220117 + 192975*log(pow2(Mst2)/pow2(Mst1))))/pow2(Mgl))*pow4(
        Mst1) + (27*(14175*Dmglst2*Mgl - 9361*pow2(Dmglst2) - 4977*pow2(Mgl) +
        490*log(pow2(Mst2)/pow2(Mst1))*(-15*Dmglst2*Mgl + pow2(Dmglst2) + 9*
        pow2(Mgl)))*pow4(Mst2))/pow2(Mgl)))/(2835.*pow3(Mst2)) + (2*pow4(Mst1)*
        (6050721 + 541080*log(pow2(Mst2)/pow2(Mst1)) + 1830600*pow2(log(pow2(
        Mst2)/pow2(Mst1))) - 395280*pow3(log(pow2(Mst2)/pow2(Mst1))) + (
        Dmglst2*(Dmglst2 + Mgl)*(-20964397 + 1700640*log(pow2(Mst2)/pow2(Mst1))
        - 4837320*pow2(log(pow2(Mst2)/pow2(Mst1))) + 408240*pow3(log(pow2(Mst2)
        /pow2(Mst1)))) + (72900*Dmsqst2*(-(pow2(Mgl)*(35 - 50*log(pow2(Mst2)/
        pow2(Mst1)) + 4*pow2(log(pow2(Mst2)/pow2(Mst1))))) + Dmglst2*(Dmglst2 +
        Mgl)*(465 - 374*log(pow2(Mst2)/pow2(Mst1)) + 12*pow2(log(pow2(Mst2)/
        pow2(Mst1)))))*(pow2(Dmsqst2)*pow2(Msq) + pow3(Dmsqst2) + Dmsqst2*pow4(
        Msq) + pow6(Msq)))/pow8(Msq))/pow2(Mgl)))/(10935.*pow3(Mst2)) + (pow2(
        Mst1)*(475.7995884773662 + 56*log(pow2(Mst2)/pow2(Mst1)) + (1088*pow2(
        log(pow2(Mst2)/pow2(Mst1))))/3. - (256*pow3(log(pow2(Mst2)/pow2(Mst1)))
        )/3. + (Dmglst2*(-322421 + 456320*log(pow2(Mst2)/pow2(Mst1)) - 526560*
        pow2(log(pow2(Mst2)/pow2(Mst1))) + 23040*pow3(log(pow2(Mst2)/pow2(Mst1)
        ))))/(810.*Mgl) - (pow2(Dmglst2)*(79386499 - 57947400*log(pow2(Mst2)/
        pow2(Mst1)) - 4460400*pow2(log(pow2(Mst2)/pow2(Mst1))) + 1451520*pow3(
        log(pow2(Mst2)/pow2(Mst1)))))/(51030.*pow2(Mgl)) + (80*Dmsqst2*(112*
        Dmglst2*Mgl*(pow2(Dmsqst2)*pow2(Msq) + pow3(Dmsqst2) + Dmsqst2*pow4(
        Msq) + pow6(Msq)) + 112*pow2(Dmglst2)*(pow2(Dmsqst2)*pow2(Msq) + pow3(
        Dmsqst2) + Dmsqst2*pow4(Msq) + pow6(Msq)) - pow2(Mgl)*(34*pow2(Dmsqst2)
        *pow2(Msq) + 45*pow3(Dmsqst2) + 23*Dmsqst2*pow4(Msq) + 12*pow6(Msq)) -
        4*log(pow2(Mst2)/pow2(Mst1))*(24*Dmglst2*Mgl*(pow2(Dmsqst2)*pow2(Msq) +
        pow3(Dmsqst2) + Dmsqst2*pow4(Msq) + pow6(Msq)) + 24*pow2(Dmglst2)*(
        pow2(Dmsqst2)*pow2(Msq) + pow3(Dmsqst2) + Dmsqst2*pow4(Msq) + pow6(Msq)
        ) - pow2(Mgl)*(10*pow2(Dmsqst2)*pow2(Msq) + 13*pow3(Dmsqst2) + 7*
        Dmsqst2*pow4(Msq) + 4*pow6(Msq)))))/(3.*pow2(Mgl)*pow8(Msq))))/Mst2 + (
        (256*(-(Dmglst2*Mgl) - 4*pow2(Dmglst2) + 2*log(pow2(Mst2)/pow2(Mst1))*(
        -(Dmglst2*Mgl) + pow2(Dmglst2) - pow2(Mgl)) + pow2(Mgl) + (3*Dmglst2*
        Mgl + 6*pow2(Dmglst2) + pow2(Mgl))*pow2(log(pow2(Mst2)/pow2(Mst1))))*
        pow3(Mst2))/(9.*pow2(Mst1)) + (16*OepS2*(pow2(Dmglst2)*(-12921*pow2(
        Mst1)*pow2(Mst2) + 51460*pow4(Mst1) - 63*pow4(Mst2)) - 3*pow2(Mgl)*(
        1131*pow2(Mst1)*pow2(Mst2) + 2740*pow4(Mst1) + 189*pow4(Mst2)) +
        Dmglst2*Mgl*(12447*pow2(Mst1)*pow2(Mst2) + 51460*pow4(Mst1) + 945*pow4(
        Mst2))))/(2187.*pow3(Mst2)) + (Mst2*(pow2(Dmglst2)*(50400*Dmsqst2*(289
        - 354*log(pow2(Mst2)/pow2(Mst1)))*(pow2(Dmsqst2)*pow2(Msq) + pow3(
        Dmsqst2) + Dmsqst2*pow4(Msq) + pow6(Msq)) + (6779161 + 21717360*log(
        pow2(Mst2)/pow2(Mst1)) - 5896800*pow2(log(pow2(Mst2)/pow2(Mst1))) -
        1481760*pow3(log(pow2(Mst2)/pow2(Mst1))))*pow8(Msq)) + 105*Dmglst2*Mgl*
        (480*Dmsqst2*(289 - 354*log(pow2(Mst2)/pow2(Mst1)))*(pow2(Dmsqst2)*
        pow2(Msq) + pow3(Dmsqst2) + Dmsqst2*pow4(Msq) + pow6(Msq)) - (30183 -
        67184*log(pow2(Mst2)/pow2(Mst1)) + 19104*pow2(log(pow2(Mst2)/pow2(Mst1)
        )) + 14112*pow3(log(pow2(Mst2)/pow2(Mst1))))*pow8(Msq)) - 189*pow2(Mgl)
        *(3*(400*pow2(Msq)*pow3(Dmsqst2) + 3400*pow4(Dmsqst2) - 2600*pow2(
        Dmsqst2)*pow4(Msq) - 5600*Dmsqst2*pow6(Msq) - 6751*pow8(Msq)) - 24800*
        pow2(log(pow2(Mst2)/pow2(Mst1)))*pow8(Msq) + 7840*pow3(log(pow2(Mst2)/
        pow2(Mst1)))*pow8(Msq) - 80*log(pow2(Mst2)/pow2(Mst1))*(750*pow2(Msq)*
        pow3(Dmsqst2) + 975*pow4(Dmsqst2) + 525*pow2(Dmsqst2)*pow4(Msq) + 300*
        Dmsqst2*pow6(Msq) + 13*pow8(Msq)))))/(17010.*pow8(Msq)))/pow2(Mgl)) - (
        z3*(16*pow4(Mt)*(1332288 + (973836*Dmglst2)/Mgl + 7776*log(pow2(Mst2)/
        pow2(Mst1)) + (4688694*pow2(Dmglst2))/pow2(Mgl) - (1215*Dmsqst2*(640*
        Dmglst2*Mgl*(Dmsqst2 + pow2(Msq)) + 640*pow2(Dmglst2)*(Dmsqst2 + pow2(
        Msq)) - pow2(Mgl)*(203*Dmsqst2 + 128*pow2(Msq))))/(pow2(Mgl)*pow4(Msq))
        + (2430*(139 - (320*Dmglst2*(Dmglst2 + Mgl))/pow2(Mgl))*pow3(Dmsqst2))/
        pow6(Msq) + (6*pow2(Mst1)*(414300 - (692304*Dmglst2)/Mgl + (1785736*
        pow2(Dmglst2))/pow2(Mgl) - (12960*Dmsqst2*(36*Dmglst2*Mgl*(pow2(
        Dmsqst2) + Dmsqst2*pow2(Msq) + pow4(Msq)) + 36*pow2(Dmglst2)*(pow2(
        Dmsqst2) + Dmsqst2*pow2(Msq) + pow4(Msq)) - pow2(Mgl)*(14*pow2(Dmsqst2)
        + 9*Dmsqst2*pow2(Msq) + 4*pow4(Msq))))/(pow2(Mgl)*pow6(Msq)) + (243*(
        1013.3333333333334 - (1920*Dmglst2*(Dmglst2 + Mgl))/pow2(Mgl))*pow4(
        Dmsqst2))/pow8(Msq)))/pow2(Mst2) + (64*pow4(Mst1)*(72723 - (295693*
        Dmglst2*(Dmglst2 + Mgl))/pow2(Mgl) + (7290*Dmsqst2*(pow2(Dmsqst2)*pow2(
        Msq) + pow3(Dmsqst2) + Dmsqst2*pow4(Msq) + pow6(Msq)))/pow8(Msq)))/
        pow4(Mst2) - (1215*(640*Dmglst2*(Dmglst2 + Mgl) - 353*pow2(Mgl))*pow4(
        Dmsqst2))/(pow2(Mgl)*pow8(Msq))) + (pow2(s2t)*((32*Mt*s2t*(2430*(pow4(
        Dmsqst2)*(8*Dmglst2*(Dmglst2 + Mgl)*pow2(Mst2)*(3*pow2(Mst1) + pow2(
        Mst2)) + pow2(Mgl)*(-20*pow2(Mst1)*pow2(Mst2) + 45*pow4(Mst1) - 36*
        pow4(Mst2))) + pow2(Dmsqst2)*pow4(Msq)*(8*Dmglst2*(Dmglst2 + Mgl)*pow2(
        Mst2)*(3*pow2(Mst1) + pow2(Mst2)) + pow2(Mgl)*(-12*pow2(Mst1)*pow2(
        Mst2) + 15*pow4(Mst1) - 28*pow4(Mst2)))) + 4860*pow2(Msq)*pow3(Dmsqst2)
        *(4*Dmglst2*(Dmglst2 + Mgl)*pow2(Mst2)*(3*pow2(Mst1) + pow2(Mst2)) +
        pow2(Mgl)*(-8*pow2(Mst1)*pow2(Mst2) + 15*pow4(Mst1) - 16*pow4(Mst2))) +
        19440*Dmsqst2*pow2(Mst2)*(Dmglst2*(Dmglst2 + Mgl)*(3*pow2(Mst1) + pow2(
        Mst2)) - pow2(Mgl)*(pow2(Mst1) + 3*pow2(Mst2)))*pow6(Msq) + (pow2(
        Dmglst2)*(171510*pow2(Mst1)*pow2(Mst2) - 282607*pow4(Mst1) - 998613*
        pow4(Mst2)) + 3*pow2(Mgl)*(19554*pow2(Mst1)*pow2(Mst2) + 22315*pow4(
        Mst1) - 53685*pow4(Mst2)) - Dmglst2*Mgl*(36810*pow2(Mst1)*pow2(Mst2) +
        282607*pow4(Mst1) + 408699*pow4(Mst2)))*pow8(Msq)))/Mst2 + (24*pow2(Mt)
        *(6*pow2(Dmglst2)*(5400*Dmsqst2*pow4(Mst1)*(pow2(Dmsqst2)*pow2(Msq) +
        pow3(Dmsqst2) + Dmsqst2*pow4(Msq) + pow6(Msq)) - (27375*pow2(Mst1)*
        pow2(Mst2) + 5872*pow4(Mst1) - 177327*pow4(Mst2))*pow8(Msq)) - pow2(
        Mgl)*(5*pow2(Dmsqst2)*pow4(Msq)*(-2037*pow2(Mst1)*pow2(Mst2) + 1756*
        pow4(Mst1) - 9819*pow4(Mst2)) + 5*pow4(Dmsqst2)*(-2079*pow2(Mst1)*pow2(
        Mst2) + 7060*pow4(Mst1) - 4401*pow4(Mst2)) + 10*pow2(Msq)*pow3(Dmsqst2)
        *(-1029*pow2(Mst1)*pow2(Mst2) + 2204*pow4(Mst1) - 3555*pow4(Mst2)) -
        80*Dmsqst2*(126*pow2(Mst1)*pow2(Mst2) + 56*pow4(Mst1) + 783*pow4(Mst2))
        *pow6(Msq) - 4*(660*pow2(Mst1)*pow2(Mst2) + 648*log(pow2(Mst2)/pow2(
        Mst1))*pow2(Mst2)*(-pow2(Mst1) + pow2(Mst2)) + 32663*pow4(Mst1) - 6237*
        pow4(Mst2))*pow8(Msq)) + 48*Dmglst2*Mgl*(675*Dmsqst2*pow4(Mst1)*(pow2(
        Dmsqst2)*pow2(Msq) + pow3(Dmsqst2) + Dmsqst2*pow4(Msq) + pow6(Msq)) - (
        -2223*pow2(Mst1)*pow2(Mst2) + 734*pow4(Mst1) - 1317*pow4(Mst2))*pow8(
        Msq))))/pow2(Mst2) + pow2(s2t)*(4*pow2(Dmglst2)*(19440*Dmsqst2*(-2*
        pow2(Mst1)*pow2(Mst2) + 5*pow4(Mst1) - 3*pow4(Mst2))*(pow2(Dmsqst2)*
        pow2(Msq) + pow3(Dmsqst2) + Dmsqst2*pow4(Msq) + pow6(Msq)) + (-740130*
        pow2(Mst1)*pow2(Mst2) + 434566*pow4(Mst1) - 2057463*pow4(Mst2))*pow8(
        Msq)) + 8*Dmglst2*Mgl*(9720*Dmsqst2*(-2*pow2(Mst1)*pow2(Mst2) + 5*pow4(
        Mst1) - 3*pow4(Mst2))*(pow2(Dmsqst2)*pow2(Msq) + pow3(Dmsqst2) +
        Dmsqst2*pow4(Msq) + pow6(Msq)) + (-101466*pow2(Mst1)*pow2(Mst2) +
        217283*pow4(Mst1) - 475605*pow4(Mst2))*pow8(Msq)) - 3*pow2(Mgl)*(5*
        pow4(Dmsqst2)*(-9198*pow2(Mst1)*pow2(Mst2) + 15679*pow4(Mst1) - 128817*
        pow4(Mst2)) + 10*pow2(Msq)*pow3(Dmsqst2)*(1086*pow2(Mst1)*pow2(Mst2) +
        6613*pow4(Mst1) - 44595*pow4(Mst2)) + 15*pow2(Dmsqst2)*pow4(Msq)*(4514*
        pow2(Mst1)*pow2(Mst2) + 3591*pow4(Mst1) - 16521*pow4(Mst2)) + 80*
        Dmsqst2*(1557*pow2(Mst1)*pow2(Mst2) + 520*pow4(Mst1) - 621*pow4(Mst2))*
        pow6(Msq) + 16*(35931*pow2(Mst1)*pow2(Mst2) + 5225*pow4(Mst1) - 162*
        log(pow2(Mst2)/pow2(Mst1))*(pow4(Mst1) - pow4(Mst2)) + 35613*pow4(Mst2)
        )*pow8(Msq)))))/(pow2(Mgl)*pow8(Msq)) - (128*s2t*pow3(Mt)*(4*pow4(Mst1)
        *(86637 + (-326791*Dmglst2*(Dmglst2 + Mgl) + (9720*Dmsqst2*(-11*
        Dmglst2*Mgl - 11*pow2(Dmglst2) + pow2(Mgl))*(pow2(Dmsqst2)*pow2(Msq) +
        pow3(Dmsqst2) + Dmsqst2*pow4(Msq) + pow6(Msq)))/pow8(Msq))/pow2(Mgl)) +
        (-9*pow4(Mst2)*(8*pow2(Dmglst2)*(1620*Dmsqst2*(pow2(Dmsqst2)*pow2(Msq)
        + pow3(Dmsqst2) + Dmsqst2*pow4(Msq) + pow6(Msq)) - 4037*pow8(Msq)) +
        27*Dmglst2*Mgl*(480*Dmsqst2*(pow2(Dmsqst2)*pow2(Msq) + pow3(Dmsqst2) +
        Dmsqst2*pow4(Msq) + pow6(Msq)) + pow8(Msq)) - 9*pow2(Mgl)*(960*pow2(
        Msq)*pow3(Dmsqst2) + 1200*pow4(Dmsqst2) + 720*pow2(Dmsqst2)*pow4(Msq) +
        480*Dmsqst2*pow6(Msq) + 961*pow8(Msq))) - 3*pow2(Mst1)*pow2(Mst2)*(
        pow2(Dmglst2)*(90720*Dmsqst2*(pow2(Dmsqst2)*pow2(Msq) + pow3(Dmsqst2) +
        Dmsqst2*pow4(Msq) + pow6(Msq)) - 124793*pow8(Msq)) + 63*Dmglst2*Mgl*(
        1440*Dmsqst2*(pow2(Dmsqst2)*pow2(Msq) + pow3(Dmsqst2) + Dmsqst2*pow4(
        Msq) + pow6(Msq)) + 1913*pow8(Msq)) - 3*pow2(Mgl)*(12960*pow2(Msq)*
        pow3(Dmsqst2) + 17280*pow4(Dmsqst2) + 8640*pow2(Dmsqst2)*pow4(Msq) +
        4320*Dmsqst2*pow6(Msq) + 24227*pow8(Msq))))/(pow2(Mgl)*pow8(Msq))))/
        pow3(Mst2)))/23328. - pow2(Mt)*((40*(3*Dmsqst2 + pow2(Msq))*pow2(Mst2)*
        pow3(Dmsqst2))/(3.*pow8(Msq)) - pow2(s2t)*(S2*((pow2(Mst1)*(27573 - (
        2232*Dmglst2)/Mgl - (254298*pow2(Dmglst2))/(7.*pow2(Mgl)) + (6900*
        Dmsqst2)/pow2(Msq) + (3350*pow2(Dmsqst2))/pow4(Msq) - (200*pow3(
        Dmsqst2))/pow6(Msq) + 30*log(pow2(Mst2)/pow2(Mst1))*(629 - (600*
        Dmglst2)/Mgl - (1194*pow2(Dmglst2))/pow2(Mgl) + (60*Dmsqst2)/pow2(Msq)
        + (10*pow2(Dmsqst2))/pow4(Msq) - (40*pow3(Dmsqst2))/pow6(Msq) - (90*
        pow4(Dmsqst2))/pow8(Msq)) - (3750*pow4(Dmsqst2))/pow8(Msq)))/45. -
        pow2(Mst2)*(42.3 + (2532*Dmglst2)/(5.*Mgl) + (20862*pow2(Dmglst2))/(35.
        *pow2(Mgl)) - (270*Dmsqst2)/pow2(Msq) - (260*pow2(Dmsqst2))/pow4(Msq) +
        log(pow2(Mst2)/pow2(Mst1))*(-159 + (84*pow2(Dmglst2))/pow2(Mgl) - (60*
        Dmsqst2)/pow2(Msq) - (40*pow2(Dmsqst2))/pow4(Msq) - (20*pow3(Dmsqst2))/
        pow6(Msq)) - (250*pow3(Dmsqst2))/pow6(Msq) - (240*pow4(Dmsqst2))/pow8(
        Msq)) + (pow4(Mst1)*(399127 + 205710*log(pow2(Mst2)/pow2(Mst1)) - (24*
        Dmglst2*(Dmglst2 + Mgl)*(28283 + 23070*log(pow2(Mst2)/pow2(Mst1))))/
        pow2(Mgl) + (400*Dmsqst2*(35 + 12*log(pow2(Mst2)/pow2(Mst1))))/pow2(
        Msq) - (300*(23 + 14*log(pow2(Mst2)/pow2(Mst1)))*pow2(Dmsqst2))/pow4(
        Msq) - (200*(139 + 66*log(pow2(Mst2)/pow2(Mst1)))*pow3(Dmsqst2))/pow6(
        Msq) - (100*(487 + 222*log(pow2(Mst2)/pow2(Mst1)))*pow4(Dmsqst2))/pow8(
        Msq)))/(270.*pow2(Mst2))) + (pow4(Mst1)*(199.98139767323744 - (
        19407863*log(pow2(Mst2)/pow2(Mst1)))/66150. + (Dmsqst2*(
        3.5881655848700444 - (470824*log(pow2(Mst2)/pow2(Mst1)))/19845. + (5*
        Dmglst2*(Dmglst2 + Mgl)*(1931 - 111*log(pow2(Mst2)/pow2(Mst1))))/(27.*
        pow2(Mgl)) - (74*pow2(log(pow2(Mst2)/pow2(Mst1))))/21.))/pow2(Msq) + (
        73657*pow2(log(pow2(Mst2)/pow2(Mst1))))/945. - (674*pow3(log(pow2(Mst2)
        /pow2(Mst1))))/27. - (Dmglst2*(Dmglst2 + Mgl)*(3650009897 - 5369876820*
        log(pow2(Mst2)/pow2(Mst1)) - 4631558400*pow2(log(pow2(Mst2)/pow2(Mst1))
        ) + 2234988000*pow3(log(pow2(Mst2)/pow2(Mst1)))))/(2.083725e7*pow2(Mgl)
        ) - (pow2(Dmsqst2)*(54.06530056349406 + (318121*log(pow2(Mst2)/pow2(
        Mst1)))/79380. + (5*Dmglst2*(Dmglst2 + Mgl)*(-1931 + 111*log(pow2(Mst2)
        /pow2(Mst1))))/(27.*pow2(Mgl)) - (16*pow2(log(pow2(Mst2)/pow2(Mst1))))/
        21.))/pow4(Msq) - ((111.71876671185817 - (623527*log(pow2(Mst2)/pow2(
        Mst1)))/39690. + (5*Dmglst2*(Dmglst2 + Mgl)*(-1931 + 111*log(pow2(Mst2)
        /pow2(Mst1))))/(27.*pow2(Mgl)) - (106*pow2(log(pow2(Mst2)/pow2(Mst1))))
        /21.)*pow3(Dmsqst2))/pow6(Msq) - ((169.37223286022228 - (401747*log(
        pow2(Mst2)/pow2(Mst1)))/11340. + (5*Dmglst2*(Dmglst2 + Mgl)*(-1931 +
        111*log(pow2(Mst2)/pow2(Mst1))))/(27.*pow2(Mgl)) - (28*pow2(log(pow2(
        Mst2)/pow2(Mst1))))/3.)*pow4(Dmsqst2))/pow8(Msq)))/pow2(Mst2) + pow2(
        Mst1)*(97.4135390946502 - (8*D3)/9. + (8*DN)/9. - (104644*log(pow2(
        Mst2)/pow2(Mst1)))/405. + (1916*pow2(log(pow2(Mst2)/pow2(Mst1))))/27. -
        (Dmsqst2*(17.999506172839506 + (2068*log(pow2(Mst2)/pow2(Mst1)))/45. +
        (20*Dmglst2*(Dmglst2 + Mgl)*(113 - 6*log(pow2(Mst2)/pow2(Mst1))))/(9.*
        pow2(Mgl)) + 8*pow2(log(pow2(Mst2)/pow2(Mst1)))))/pow2(Msq) + (4*
        Dmglst2*(1216808 - 1164855*log(pow2(Mst2)/pow2(Mst1)) + 365400*pow2(
        log(pow2(Mst2)/pow2(Mst1))) - 650250*pow3(log(pow2(Mst2)/pow2(Mst1)))))
        /(30375.*Mgl) - (1166*pow3(log(pow2(Mst2)/pow2(Mst1))))/27. - (pow2(
        Dmglst2)*(-273721621 + 194650260*log(pow2(Mst2)/pow2(Mst1)) + 84886200*
        pow2(log(pow2(Mst2)/pow2(Mst1))) + 54369000*pow3(log(pow2(Mst2)/pow2(
        Mst1)))))/(425250.*pow2(Mgl)) + (2*pow2(Dmsqst2)*(1241 - 19935*log(
        pow2(Mst2)/pow2(Mst1)) + (1350*Dmglst2*(Dmglst2 + Mgl)*(-113 + 6*log(
        pow2(Mst2)/pow2(Mst1))))/pow2(Mgl) - 2025*pow2(log(pow2(Mst2)/pow2(
        Mst1)))))/(1215.*pow4(Msq)) + ((22.085102880658436 - (2656*log(pow2(
        Mst2)/pow2(Mst1)))/135. + (20*Dmglst2*(Dmglst2 + Mgl)*(-113 + 6*log(
        pow2(Mst2)/pow2(Mst1))))/(9.*pow2(Mgl)) + (4*pow2(log(pow2(Mst2)/pow2(
        Mst1))))/3.)*pow3(Dmsqst2))/pow6(Msq) + ((42.12740740740741 - (98*log(
        pow2(Mst2)/pow2(Mst1)))/15. + (20*Dmglst2*(Dmglst2 + Mgl)*(-113 + 6*
        log(pow2(Mst2)/pow2(Mst1))))/(9.*pow2(Mgl)) + 6*pow2(log(pow2(Mst2)/
        pow2(Mst1))))*pow4(Dmsqst2))/pow8(Msq)) - (pow2(Mst2)*(-5040*(602*
        Dmglst2*Mgl + 473*pow2(Dmglst2) + 266*pow2(Mgl))*pow2(log(pow2(Mst2)/
        pow2(Mst1)))*pow8(Msq) + 5040*(128*Dmglst2*Mgl + 192*pow2(Dmglst2) + 7*
        pow2(Mgl))*pow3(log(pow2(Mst2)/pow2(Mst1)))*pow8(Msq) - 7*pow2(Mgl)*(
        335900*pow2(Msq)*pow3(Dmsqst2) + 327900*pow4(Dmsqst2) + 343900*pow2(
        Dmsqst2)*pow4(Msq) + 351900*Dmsqst2*pow6(Msq) + (1013263 + 2880*D3 -
        2880*DN)*pow8(Msq)) + 1344*Dmglst2*Mgl*(1725*pow2(Msq)*pow3(Dmsqst2) +
        1725*pow4(Dmsqst2) + 1725*pow2(Dmsqst2)*pow4(Msq) + 1725*Dmsqst2*pow6(
        Msq) + (-5837 + 60*B4 - 60*D3 + 30*DN)*pow8(Msq)) + 36*pow2(Dmglst2)*(
        64400*pow2(Msq)*pow3(Dmsqst2) + 64400*pow4(Dmsqst2) + 64400*pow2(
        Dmsqst2)*pow4(Msq) + 64400*Dmsqst2*pow6(Msq) + (-939499 + 3360*B4 -
        3360*D3 + 1680*DN)*pow8(Msq)) - 1680*log(pow2(Mst2)/pow2(Mst1))*(360*
        pow2(Mgl)*pow2(Msq)*pow3(Dmsqst2) + 300*pow2(Mgl)*pow4(Dmsqst2) + 420*
        pow2(Dmsqst2)*pow2(Mgl)*pow4(Msq) + 480*Dmsqst2*pow2(Mgl)*pow6(Msq) - (
        6854*Dmglst2*Mgl + 13339*pow2(Dmglst2) + 2267*pow2(Mgl))*pow8(Msq))))/(
        22680.*pow2(Mgl)*pow8(Msq)) + ((32*(-4*pow2(Dmglst2) + pow2(Mgl) - 2*
        log(pow2(Mst2)/pow2(Mst1))*pow2(Dmglst2 + Mgl) + (4*Dmglst2*Mgl + 10*
        pow2(Dmglst2) + pow2(Mgl))*pow2(log(pow2(Mst2)/pow2(Mst1))))*pow6(Mst2)
        )/(9.*pow4(Mst1)) + (4*OepS2*(20*Dmsqst2*pow2(Mgl)*(pow3(Dmsqst2)*(27*
        pow2(Mst1)*pow2(Mst2) + 37*pow4(Mst1)) + Dmsqst2*pow4(Msq)*(-3*pow2(
        Mst1)*pow2(Mst2) + 7*pow4(Mst1) - 18*pow4(Mst2)) + pow2(Dmsqst2)*pow2(
        Msq)*(12*pow2(Mst1)*pow2(Mst2) + 22*pow4(Mst1) - 9*pow4(Mst2)) - (18*
        pow2(Mst1)*pow2(Mst2) + 8*pow4(Mst1) + 27*pow4(Mst2))*pow6(Msq)) + (24*
        Dmglst2*Mgl*pow2(Mst1)*(769*pow2(Mst1) + 150*pow2(Mst2)) + 12*pow2(
        Dmglst2)*(597*pow2(Mst1)*pow2(Mst2) + 1538*pow4(Mst1) + 63*pow4(Mst2))
        - pow2(Mgl)*(3774*pow2(Mst1)*pow2(Mst2) + 6857*pow4(Mst1) + 1431*pow4(
        Mst2)))*pow8(Msq)))/(729.*pow2(Mst2)*pow8(Msq)) - (2*pow4(Mst2)*(20*
        pow2(Dmglst2)*(18*pow2(Msq)*pow3(Dmsqst2) + 18*pow4(Dmsqst2) + 18*pow2(
        Dmsqst2)*pow4(Msq) + 18*Dmsqst2*pow6(Msq) - 41*pow8(Msq)) - 3*pow2(Mgl)
        *(30*pow2(Msq)*pow3(Dmsqst2) + 50*pow4(Dmsqst2) + 10*pow2(Dmsqst2)*
        pow4(Msq) - 10*Dmsqst2*pow6(Msq) - 29*pow8(Msq)) + 8*Dmglst2*Mgl*(45*
        pow2(Msq)*pow3(Dmsqst2) + 45*pow4(Dmsqst2) + 45*pow2(Dmsqst2)*pow4(Msq)
        + 45*Dmsqst2*pow6(Msq) - 13*pow8(Msq)) + (182*Dmglst2*Mgl + 337*pow2(
        Dmglst2) + 75*pow2(Mgl))*pow2(log(pow2(Mst2)/pow2(Mst1)))*pow8(Msq) +
        2*log(pow2(Mst2)/pow2(Mst1))*(30*pow2(Mgl)*pow2(Msq)*pow3(Dmsqst2) +
        30*pow2(Mgl)*pow4(Dmsqst2) + 30*pow2(Dmsqst2)*pow2(Mgl)*pow4(Msq) + 30*
        Dmsqst2*pow2(Mgl)*pow6(Msq) - (134*Dmglst2*Mgl + 83*pow2(Dmglst2) + 77*
        pow2(Mgl))*pow8(Msq))))/(9.*pow2(Mst1)*pow8(Msq)))/pow2(Mgl))) + pow4(
        s2t)*(pow4(Mst1)*(46.745471895783595 + B4 + D3/9. - DN/9. + (23468297*
        log(pow2(Mst2)/pow2(Mst1)))/529200. + (11764*pow2(log(pow2(Mst2)/pow2(
        Mst1))))/945. + (Dmsqst2*(15.13649301931237 + (621671*log(pow2(Mst2)/
        pow2(Mst1)))/39690. + (5*Dmglst2*(Dmglst2 + Mgl)*(76 - 249*log(pow2(
        Mst2)/pow2(Mst1))))/(54.*pow2(Mgl)) + (53*pow2(log(pow2(Mst2)/pow2(
        Mst1))))/21.))/pow2(Msq) + (Dmglst2*(Dmglst2 + Mgl)*(79.58863384550371
         - (8287903*log(pow2(Mst2)/pow2(Mst1)))/
        1.1907e6 - (51*pow2(log(pow2(Mst2)/pow2(Mst1))))/70. - (11*pow3(log(
        pow2(Mst2)/pow2(Mst1))))/9.))/pow2(Mgl) - (409*pow3(log(pow2(Mst2)/
        pow2(Mst1))))/108. + (pow2(Dmsqst2)*(13.217310375649378 + (281161*log(
        pow2(Mst2)/pow2(Mst1)))/17640. + (5*Dmglst2*(Dmglst2 + Mgl)*(76 - 249*
        log(pow2(Mst2)/pow2(Mst1))))/(54.*pow2(Mgl)) + (31*pow2(log(pow2(Mst2)/
        pow2(Mst1))))/42.))/pow4(Msq) + ((11.298127731986387 + (1287107*log(
        pow2(Mst2)/pow2(Mst1)))/79380. + (5*Dmglst2*(Dmglst2 + Mgl)*(76 - 249*
        log(pow2(Mst2)/pow2(Mst1))))/(54.*pow2(Mgl)) - (22*pow2(log(pow2(Mst2)/
        pow2(Mst1))))/21.)*pow3(Dmsqst2))/pow6(Msq) + ((9.378945088323395 + (
        373997*log(pow2(Mst2)/pow2(Mst1)))/22680. + (5*Dmglst2*(Dmglst2 + Mgl)*
        (76 - 249*log(pow2(Mst2)/pow2(Mst1))))/(54.*pow2(Mgl)) - (17*pow2(log(
        pow2(Mst2)/pow2(Mst1))))/6.)*pow4(Dmsqst2))/pow8(Msq)) - pow2(Mst1)*
        pow2(Mst2)*(79.01365226337448 - B4/9. + (2*D3)/9. - DN/6. - (16051*log(
        pow2(Mst2)/pow2(Mst1)))/810. - (92*pow2(log(pow2(Mst2)/pow2(Mst1))))/
        27. - (Dmsqst2*(3.2221604938271606 + (317*log(pow2(Mst2)/pow2(Mst1)))/
        90. + (5*Dmglst2*(Dmglst2 + Mgl)*(-75 + 26*log(pow2(Mst2)/pow2(Mst1))))
        /(6.*pow2(Mgl)) + (31*pow2(log(pow2(Mst2)/pow2(Mst1))))/6.))/pow2(Msq)
        + (115*pow3(log(pow2(Mst2)/pow2(Mst1))))/27. + (Dmglst2*(
        0.7822304526748971 - (2*B4)/3. + (8*D3)/9. - (5*DN)/9. + (81643*log(
        pow2(Mst2)/pow2(Mst1)))/4050. + (4969*pow2(log(pow2(Mst2)/pow2(Mst1))))
        /270. + (58*pow3(log(pow2(Mst2)/pow2(Mst1))))/27.))/Mgl + (pow2(
        Dmglst2)*(30.209968449931413 - B4 + (4*D3)/3. - (5*DN)/6. + (165199*
        log(pow2(Mst2)/pow2(Mst1)))/2700. + (14231*pow2(log(pow2(Mst2)/pow2(
        Mst1))))/540. + (97*pow3(log(pow2(Mst2)/pow2(Mst1))))/27.))/pow2(Mgl) -
        (pow2(Dmsqst2)*(23.10627572016461 + (13*log(pow2(Mst2)/pow2(Mst1)))/
        108. + (5*Dmglst2*(Dmglst2 + Mgl)*(-75 + 26*log(pow2(Mst2)/pow2(Mst1)))
        )/(6.*pow2(Mgl)) + (25*pow2(log(pow2(Mst2)/pow2(Mst1))))/6.))/pow4(Msq)
        - ((42.99039094650206 - (443*log(pow2(Mst2)/pow2(Mst1)))/135. + (5*
        Dmglst2*(Dmglst2 + Mgl)*(-75 + 26*log(pow2(Mst2)/pow2(Mst1))))/(6.*
        pow2(Mgl)) + (19*pow2(log(pow2(Mst2)/pow2(Mst1))))/6.)*pow3(Dmsqst2))/
        pow6(Msq) - ((62.8745061728395 - (401*log(pow2(Mst2)/pow2(Mst1)))/60. +
        (5*Dmglst2*(Dmglst2 + Mgl)*(-75 + 26*log(pow2(Mst2)/pow2(Mst1))))/(6.*
        pow2(Mgl)) + (13*pow2(log(pow2(Mst2)/pow2(Mst1))))/6.)*pow4(Dmsqst2))/
        pow8(Msq)) - S2*((pow4(Mst1)*(50550 + 16155*log(pow2(Mst2)/pow2(Mst1))
        - (2*Dmglst2*(Dmglst2 + Mgl)*(51635 + 25194*log(pow2(Mst2)/pow2(Mst1)))
        )/pow2(Mgl) + (60*Dmsqst2*(163 + 42*log(pow2(Mst2)/pow2(Mst1))))/pow2(
        Msq) + (3000*pow2(Dmsqst2))/pow4(Msq) - (1260*(3 + 2*log(pow2(Mst2)/
        pow2(Mst1)))*pow3(Dmsqst2))/pow6(Msq) - (240*(44 + 21*log(pow2(Mst2)/
        pow2(Mst1)))*pow4(Dmsqst2))/pow8(Msq)))/324. + pow4(Mst2)*((1677*
        Dmglst2)/(2.*Mgl) + (116129*pow2(Dmglst2))/(60.*pow2(Mgl)) + log(pow2(
        Mst2)/pow2(Mst1))*(34.5 - (7*Dmglst2)/Mgl - (35*pow2(Dmglst2))/(6.*
        pow2(Mgl)) + (15*Dmsqst2)/pow2(Msq) + (10*pow2(Dmsqst2))/pow4(Msq) + (
        5*pow3(Dmsqst2))/pow6(Msq)) - (1310*pow2(Msq)*pow3(Dmsqst2) + 1740*
        pow4(Dmsqst2) + 880*pow2(Dmsqst2)*pow4(Msq) + 450*Dmsqst2*pow6(Msq) -
        921*pow8(Msq))/(4.*pow8(Msq))) + pow2(Mst1)*pow2(Mst2)*(log(pow2(Mst2)/
        pow2(Mst1))*(107.16666666666667 - (72*Dmglst2)/Mgl - (80*pow2(Dmglst2))
        /(9.*pow2(Mgl)) + (35*Dmsqst2)/pow2(Msq) + (65*pow2(Dmsqst2))/(3.*pow4(
        Msq)) + (25*pow3(Dmsqst2))/(3.*pow6(Msq)) - (5*pow4(Dmsqst2))/pow8(Msq)
        ) + ((33000*Dmglst2)/Mgl + (89092*pow2(Dmglst2))/pow2(Mgl) + (15*(5150*
        pow2(Msq)*pow3(Dmsqst2) + 3930*pow4(Dmsqst2) + 6370*pow2(Dmsqst2)*pow4(
        Msq) + 7590*Dmsqst2*pow6(Msq) + 17151*pow8(Msq)))/pow8(Msq))/540.)) + (
        (OepS2*(60*Dmsqst2*pow2(Mgl)*(3*Dmsqst2*pow2(Mst2)*(13*pow2(Mst1) + 6*
        pow2(Mst2))*pow4(Msq) - pow3(Dmsqst2)*(9*pow2(Mst1)*pow2(Mst2) + 28*
        pow4(Mst1)) + pow2(Dmsqst2)*pow2(Msq)*(15*pow2(Mst1)*pow2(Mst2) - 14*
        pow4(Mst1) + 9*pow4(Mst2)) + (63*pow2(Mst1)*pow2(Mst2) + 14*pow4(Mst1)
        + 27*pow4(Mst2))*pow6(Msq)) - (4*Dmglst2*Mgl*(1944*pow2(Mst1)*pow2(
        Mst2) + 4199*pow4(Mst1) + 189*pow4(Mst2)) + 2*pow2(Dmglst2)*(480*pow2(
        Mst1)*pow2(Mst2) + 8398*pow4(Mst1) + 315*pow4(Mst2)) - 3*pow2(Mgl)*(
        3858*pow2(Mst1)*pow2(Mst2) + 1795*pow4(Mst1) + 1242*pow4(Mst2)))*pow8(
        Msq)))/(2187.*pow8(Msq)) + (pow4(Mst2)*(pow2(Dmglst2)*(129600*Dmsqst2*(
        1 + log(pow2(Mst2)/pow2(Mst1)))*(pow2(Dmsqst2)*pow2(Msq) + pow3(
        Dmsqst2) + Dmsqst2*pow4(Msq) + pow6(Msq)) + (1395899 - 54000*B4 +
        60480*D3 - 33480*DN + 1098000*log(pow2(Mst2)/pow2(Mst1)) - 14040*pow2(
        log(pow2(Mst2)/pow2(Mst1))) - 346680*pow3(log(pow2(Mst2)/pow2(Mst1))))*
        pow8(Msq)) + 30*Dmglst2*Mgl*(4320*Dmsqst2*(1 + log(pow2(Mst2)/pow2(
        Mst1)))*(pow2(Dmsqst2)*pow2(Msq) + pow3(Dmsqst2) + Dmsqst2*pow4(Msq) +
        pow6(Msq)) + (15707 - 432*B4 + 576*D3 - 360*DN + 19944*log(pow2(Mst2)/
        pow2(Mst1)) + 10224*pow2(log(pow2(Mst2)/pow2(Mst1))) - 7704*pow3(log(
        pow2(Mst2)/pow2(Mst1))))*pow8(Msq)) - 15*pow2(Mgl)*(47000*pow2(Msq)*
        pow3(Dmsqst2) + 52305*pow4(Dmsqst2) + 41695*pow2(Dmsqst2)*pow4(Msq) +
        36390*Dmsqst2*pow6(Msq) + 72*pow2(log(pow2(Mst2)/pow2(Mst1)))*(60*pow2(
        Msq)*pow3(Dmsqst2) + 60*pow4(Dmsqst2) + 60*pow2(Dmsqst2)*pow4(Msq) +
        60*Dmsqst2*pow6(Msq) - 371*pow8(Msq)) + (25289 + 1440*B4 - 144*D3 + 72*
        DN)*pow8(Msq) + 9756*pow3(log(pow2(Mst2)/pow2(Mst1)))*pow8(Msq) + 24*
        log(pow2(Mst2)/pow2(Mst1))*(540*pow2(Msq)*pow3(Dmsqst2) + 555*pow4(
        Dmsqst2) + 525*pow2(Dmsqst2)*pow4(Msq) + 510*Dmsqst2*pow6(Msq) + 1148*
        pow8(Msq))) + (540*pow2(Mst2)*(5*(4*pow2(Dmglst2)*(18*pow2(Msq)*pow3(
        Dmsqst2) + 18*pow4(Dmsqst2) + 18*pow2(Dmsqst2)*pow4(Msq) + 18*Dmsqst2*
        pow6(Msq) - 49*pow8(Msq)) - 3*pow2(Mgl)*(6*pow2(Msq)*pow3(Dmsqst2) +
        10*pow4(Dmsqst2) + 2*pow2(Dmsqst2)*pow4(Msq) - 2*Dmsqst2*pow6(Msq) - 9*
        pow8(Msq)) + 8*Dmglst2*Mgl*(9*pow2(Msq)*pow3(Dmsqst2) + 9*pow4(Dmsqst2)
        + 9*pow2(Dmsqst2)*pow4(Msq) + 9*Dmsqst2*pow6(Msq) - pow8(Msq))) + (374*
        Dmglst2*Mgl + 817*pow2(Dmglst2) + 123*pow2(Mgl))*pow2(log(pow2(Mst2)/
        pow2(Mst1)))*pow8(Msq) + log(pow2(Mst2)/pow2(Mst1))*(60*pow2(Mgl)*pow2(
        Msq)*pow3(Dmsqst2) + 60*pow2(Mgl)*pow4(Dmsqst2) + 60*pow2(Dmsqst2)*
        pow2(Mgl)*pow4(Msq) + 60*Dmsqst2*pow2(Mgl)*pow6(Msq) - 2*(262*Dmglst2*
        Mgl + 211*pow2(Dmglst2) + 125*pow2(Mgl))*pow8(Msq))))/pow2(Mst1)))/(
        19440.*pow8(Msq)) - (8*pow2(Dmglst2)*(-2 + log(pow2(Mst2)/pow2(Mst1)) +
        3*pow2(log(pow2(Mst2)/pow2(Mst1))))*pow8(Mst2))/(9.*pow4(Mst1)))/pow2(
        Mgl))))/pow4(Mt)*12.;

   return result;
}

/**
 * 	@return returns the susy log^1 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double himalaya::H6bq2g2::calc_coef_at_as2_no_sm_logs_log1() const {

   const double result =
      (((3000*pow2(Msq)*pow2(Mst1)*pow3(Dmsqst2)*(8*Dmglst2*Mgl*
        Mst2*pow2(Mst1)*pow3(Mt)*(36*Mst2*(Mt*(33 - 18*z2) + 2*Mst2*s2t*(-12 +
        7*z2))*pow2(Mst1) + (Mt*(283 - 180*z2) + 12*Mst2*s2t*(-29 + 18*z2))*
        pow3(Mst2) + 18*s2t*(-83 + 44*z2)*pow4(Mst1)) + 8*Mst2*pow2(Dmglst2)*
        pow2(Mst1)*pow3(Mt)*(36*Mst2*(Mt*(33 - 18*z2) + 2*Mst2*s2t*(-12 + 7*z2)
        )*pow2(Mst1) + (Mt*(283 - 180*z2) + 12*Mst2*s2t*(-29 + 18*z2))*pow3(
        Mst2) + 18*s2t*(-83 + 44*z2)*pow4(Mst1)) + pow2(Mgl)*(9*pow2(Mst2)*
        pow4(Mst1)*(8*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 64*Mst2*s2t*(-5 + 3*z2)*
        pow3(Mt) + 16*(-23 + 14*z2)*pow4(Mt) + pow4(Mst2)*pow4(s2t)) + pow2(
        Mst1)*pow4(Mst2)*(-144*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 144*Mst2*s2t*(-
        13 + 8*z2)*pow3(Mt) + 4*(-353 + 180*z2)*pow4(Mt) + 9*pow4(Mst2)*pow4(
        s2t)) + 9*(-16*Mst2*s2t*(-9 + 4*z2)*pow3(Mt) + 8*(-31 + 12*z2)*pow4(Mt)
        - pow4(Mst2)*pow4(s2t))*pow6(Mst1) - 9*pow2(-4*pow2(Mt) + pow2(Mst2)*
        pow2(s2t))*pow6(Mst2))) + 3000*pow2(Dmsqst2)*pow2(Mst1)*pow4(Msq)*(8*
        Dmglst2*Mgl*Mst2*pow2(Mst1)*pow3(Mt)*(36*Mst2*(Mt*(33 - 18*z2) + 2*
        Mst2*s2t*(-12 + 7*z2))*pow2(Mst1) + (Mt*(283 - 180*z2) + 12*Mst2*s2t*(-
        29 + 18*z2))*pow3(Mst2) + 18*s2t*(-83 + 44*z2)*pow4(Mst1)) + 8*Mst2*
        pow2(Dmglst2)*pow2(Mst1)*pow3(Mt)*(36*Mst2*(Mt*(33 - 18*z2) + 2*Mst2*
        s2t*(-12 + 7*z2))*pow2(Mst1) + (Mt*(283 - 180*z2) + 12*Mst2*s2t*(-29 +
        18*z2))*pow3(Mst2) + 18*s2t*(-83 + 44*z2)*pow4(Mst1)) + pow2(Mgl)*(9*
        pow2(Mst2)*pow4(Mst1)*(8*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 32*Mst2*s2t*(-
        7 + 4*z2)*pow3(Mt) + 16*(-16 + 9*z2)*pow4(Mt) + pow4(Mst2)*pow4(s2t)) +
        pow2(Mst1)*pow4(Mst2)*(-144*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 72*Mst2*
        s2t*(-19 + 12*z2)*pow3(Mt) + 2*(-545 + 252*z2)*pow4(Mt) + 9*pow4(Mst2)*
        pow4(s2t)) + 9*(-16*Mst2*s2t*(-9 + 4*z2)*pow3(Mt) + 8*(-31 + 12*z2)*
        pow4(Mt) - pow4(Mst2)*pow4(s2t))*pow6(Mst1) - 9*pow2(-4*pow2(Mt) +
        pow2(Mst2)*pow2(s2t))*pow6(Mst2))) + 3000*Dmsqst2*pow2(Mst1)*pow6(Msq)*
        (8*Dmglst2*Mgl*Mst2*pow2(Mst1)*pow3(Mt)*(36*Mst2*(Mt*(33 - 18*z2) + 2*
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
        pow2(s2t))*pow6(Mst2))) + 3000*pow2(Mst1)*pow4(Dmsqst2)*(8*Dmglst2*Mgl*
        Mst2*pow2(Mst1)*pow3(Mt)*(36*Mst2*(Mt*(33 - 18*z2) + 2*Mst2*s2t*(-12 +
        7*z2))*pow2(Mst1) + (Mt*(283 - 180*z2) + 12*Mst2*s2t*(-29 + 18*z2))*
        pow3(Mst2) + 18*s2t*(-83 + 44*z2)*pow4(Mst1)) + 8*Mst2*pow2(Dmglst2)*
        pow2(Mst1)*pow3(Mt)*(36*Mst2*(Mt*(33 - 18*z2) + 2*Mst2*s2t*(-12 + 7*z2)
        )*pow2(Mst1) + (Mt*(283 - 180*z2) + 12*Mst2*s2t*(-29 + 18*z2))*pow3(
        Mst2) + 18*s2t*(-83 + 44*z2)*pow4(Mst1)) + 3*pow2(Mgl)*(3*pow2(Mst2)*
        pow4(Mst1)*(8*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 32*Mst2*s2t*(-13 + 8*z2)*
        pow3(Mt) + 16*(-30 + 19*z2)*pow4(Mt) + pow4(Mst2)*pow4(s2t)) + pow2(
        Mst1)*pow4(Mst2)*(-48*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 24*Mst2*s2t*(33 -
        20*z2)*pow3(Mt) + (-578 + 312*z2)*pow4(Mt) + 3*pow4(Mst2)*pow4(s2t)) +
        3*(-16*Mst2*s2t*(-9 + 4*z2)*pow3(Mt) + 8*(-31 + 12*z2)*pow4(Mt) - pow4(
        Mst2)*pow4(s2t))*pow6(Mst1) - 3*pow2(-4*pow2(Mt) + pow2(Mst2)*pow2(s2t)
        )*pow6(Mst2))) - 450*pow2(log(pow2(Mst2)/pow2(Mst1)))*pow4(Mst1)*(pow2(
        Mgl)*(pow4(Mst2)*(-296*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 6784*Mst2*s2t*
        pow3(Mt) - 2208*Mt*pow3(Mst2)*pow3(s2t) + 672*pow4(Mt) - 519*pow4(Mst2)
        *pow4(s2t)) + pow4(Mst1)*(-2688*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 5056*
        Mst2*s2t*pow3(Mt) - 1024*Mt*pow3(Mst2)*pow3(s2t) + 2848*pow4(Mt) - 41*
        pow4(Mst2)*pow4(s2t)) - 8*pow2(Mst1)*pow2(Mst2)*(425*pow2(Mst2)*pow2(
        Mt)*pow2(s2t) - 704*Mst2*s2t*pow3(Mt) + 212*Mt*pow3(Mst2)*pow3(s2t) -
        76*pow4(Mt) + 3*pow4(Mst2)*pow4(s2t))) + 8*Dmglst2*Mgl*(8*Mt*(-157*
        Mst2*s2t*pow2(Mt) - 108*Mt*pow2(Mst2)*pow2(s2t) + 398*pow3(Mt) - 16*
        pow3(Mst2)*pow3(s2t))*pow4(Mst1) + pow4(Mst2)*(-320*pow2(Mst2)*pow2(Mt)
        *pow2(s2t) + 848*Mst2*s2t*pow3(Mt) - 532*Mt*pow3(Mst2)*pow3(s2t) + 512*
        pow4(Mt) - 99*pow4(Mst2)*pow4(s2t)) + pow2(Mst1)*pow2(Mst2)*(-768*pow2(
        Mst2)*pow2(Mt)*pow2(s2t) - 448*Mst2*s2t*pow3(Mt) - 212*Mt*pow3(Mst2)*
        pow3(s2t) + 1896*pow4(Mt) + 35*pow4(Mst2)*pow4(s2t))) + 4*pow2(Dmglst2)
        *(16*Mt*(-157*Mst2*s2t*pow2(Mt) - 108*Mt*pow2(Mst2)*pow2(s2t) + 398*
        pow3(Mt) - 16*pow3(Mst2)*pow3(s2t))*pow4(Mst1) + pow4(Mst2)*(-960*pow2(
        Mst2)*pow2(Mt)*pow2(s2t) + 1696*Mst2*s2t*pow3(Mt) - 1832*Mt*pow3(Mst2)*
        pow3(s2t) + 1536*pow4(Mt) - 297*pow4(Mst2)*pow4(s2t)) + pow2(Mst1)*
        pow2(Mst2)*(-2304*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 256*Mst2*s2t*pow3(Mt)
        - 424*Mt*pow3(Mst2)*pow3(s2t) + 152*pow4(Mt) + 105*pow4(Mst2)*pow4(s2t)
        )))*pow8(Msq) + pow8(Msq)*(-100*Dmglst2*Mgl*(2*pow4(Mst1)*pow4(Mst2)*(
        12*(-3565 + 1728*z2)*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 8*Mst2*s2t*(-3659
        + 252*z2)*pow3(Mt) + 72*Mt*(-192 + 91*z2)*pow3(Mst2)*pow3(s2t) + 44*(
        427 + 432*z2)*pow4(Mt) + 9*(211 - 86*z2)*pow4(Mst2)*pow4(s2t)) + 6*
        pow2(Mst2)*(4*(-1955 + 1152*z2)*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 48*
        Mst2*s2t*(-13 + 262*z2)*pow3(Mt) + 24*Mt*(111 - 17*z2)*pow3(Mst2)*pow3(
        s2t) + 8*(-1471 + 3240*z2)*pow4(Mt) + 3*(-1 + 86*z2)*pow4(Mst2)*pow4(
        s2t))*pow6(Mst1) - 18*pow2(Mst1)*(-536*pow2(Mst2)*pow2(Mt)*pow2(s2t) -
        512*Mst2*s2t*pow3(Mt) + 128*Mt*pow3(Mst2)*pow3(s2t) + 560*pow4(Mt) +
        131*pow4(Mst2)*pow4(s2t))*pow6(Mst2) + 3*(192*(-77 + 48*z2)*pow2(Mst2)*
        pow2(Mt)*pow2(s2t) - 48*Mst2*s2t*(-867 + 1012*z2)*pow3(Mt) - 276*Mt*
        pow3(Mst2)*pow3(s2t) + 16*(-6445 + 6768*z2)*pow4(Mt) - 441*pow4(Mst2)*
        pow4(s2t))*pow8(Mst1) - 2304*pow2(Mt)*(-2*pow2(Mt) + pow2(Mst2)*pow2(
        s2t))*pow8(Mst2)) + 25*pow2(Mgl)*(4*pow4(Mst1)*pow4(Mst2)*(-72*(-445 +
        96*z2)*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 144*Mst2*s2t*(-121 + 296*z2)*
        pow3(Mt) - 144*Mt*(-151 + 19*z2)*pow3(Mst2)*pow3(s2t) + 8*(-1993 +
        2151*z2 + 216*z3)*pow4(Mt) + 81*(8 + 5*z2)*pow4(Mst2)*pow4(s2t)) - 36*
        pow2(Mst2)*(-696*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 32*Mst2*s2t*(-179 +
        128*z2)*pow3(Mt) + 16*Mt*(148 - 17*z2)*pow3(Mst2)*pow3(s2t) - 8*(-403 +
        300*z2)*pow4(Mt) + (551 + 86*z2)*pow4(Mst2)*pow4(s2t))*pow6(Mst1) + 36*
        pow2(Mst1)*(-744*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 1024*Mst2*s2t*pow3(Mt)
        + 256*Mt*pow3(Mst2)*pow3(s2t) + 1232*pow4(Mt) + 141*pow4(Mst2)*pow4(
        s2t))*pow6(Mst2) + 9*(864*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 512*Mst2*s2t*
        (-49 + 31*z2)*pow3(Mt) - 784*Mt*pow3(Mst2)*pow3(s2t) + 8*(-2503 + 1408*
        z2)*pow4(Mt) + (1547 + 164*z2)*pow4(Mst2)*pow4(s2t))*pow8(Mst1) + 4608*
        pow2(Mt)*(-2*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow8(Mst2)) - 4*pow2(
        Dmglst2)*(-2*pow4(Mst1)*pow4(Mst2)*(-150*(-14251 + 9792*z2)*pow2(Mst2)*
        pow2(Mt)*pow2(s2t) - 600*Mst2*s2t*(-2509 + 1278*z2)*pow3(Mt) - 300*Mt*(
        -2027 + 1194*z2)*pow3(Mst2)*pow3(s2t) + 4*(25949 + 33300*z2)*pow4(Mt) +
        75*(-838 + 387*z2)*pow4(Mst2)*pow4(s2t)) + 25*pow2(Mst2)*(12*(-11813 +
        8064*z2)*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 96*Mst2*s2t*(-1016 + 1395*z2)*
        pow3(Mt) - 24*Mt*(-83 + 102*z2)*pow3(Mst2)*pow3(s2t) - 64*(347 + 3330*
        z2)*pow4(Mt) + 3*(-415 + 774*z2)*pow4(Mst2)*pow4(s2t))*pow6(Mst1) -
        225*pow2(Mst1)*(-408*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 1024*Mst2*s2t*
        pow3(Mt) - 256*Mt*pow3(Mst2)*pow3(s2t) - 720*pow4(Mt) + 179*pow4(Mst2)*
        pow4(s2t))*pow6(Mst2) + 75*(192*(-77 + 48*z2)*pow2(Mst2)*pow2(Mt)*pow2(
        s2t) - 48*Mst2*s2t*(-867 + 1012*z2)*pow3(Mt) - 276*Mt*pow3(Mst2)*pow3(
        s2t) + 16*(-6445 + 6768*z2)*pow4(Mt) - 441*pow4(Mst2)*pow4(s2t))*pow8(
        Mst1) - 3600*(8*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 16*pow4(Mt) + pow4(
        Mst2)*pow4(s2t))*pow8(Mst2))) + 30*log(pow2(Mst2)/pow2(Mst1))*(1800*
        pow2(Msq)*pow3(Dmsqst2)*pow4(Mst1)*(16*Dmglst2*Mgl*Mst2*pow2(Mst1)*(-(
        Mst2*Mt) + 3*s2t*pow2(Mst1))*pow3(Mt) + 16*Mst2*pow2(Dmglst2)*pow2(
        Mst1)*(-(Mst2*Mt) + 3*s2t*pow2(Mst1))*pow3(Mt) + pow2(Mgl)*(8*(5*Mt -
        2*Mst2*s2t)*pow3(Mt)*pow4(Mst1) - pow2(Mst1)*pow4(s2t)*pow6(Mst2) +
        pow4(s2t)*pow8(Mst2))) + 1800*pow2(Dmsqst2)*pow4(Msq)*pow4(Mst1)*(16*
        Dmglst2*Mgl*Mst2*pow2(Mst1)*(-(Mst2*Mt) + 3*s2t*pow2(Mst1))*pow3(Mt) +
        16*Mst2*pow2(Dmglst2)*pow2(Mst1)*(-(Mst2*Mt) + 3*s2t*pow2(Mst1))*pow3(
        Mt) + pow2(Mgl)*(8*(5*Mt - 2*Mst2*s2t)*pow3(Mt)*pow4(Mst1) - 4*pow4(
        Mst2)*pow4(Mt) + pow2(Mst1)*(4*pow2(Mst2)*pow4(Mt) - pow4(s2t)*pow6(
        Mst2)) + pow4(s2t)*pow8(Mst2))) + 1800*Dmsqst2*pow4(Mst1)*pow6(Msq)*(
        16*Dmglst2*Mgl*Mst2*pow2(Mst1)*(-(Mst2*Mt) + 3*s2t*pow2(Mst1))*pow3(Mt)
        + 16*Mst2*pow2(Dmglst2)*pow2(Mst1)*(-(Mst2*Mt) + 3*s2t*pow2(Mst1))*
        pow3(Mt) + pow2(Mgl)*(8*(5*Mt - 2*Mst2*s2t)*pow3(Mt)*pow4(Mst1) - 8*
        pow4(Mst2)*pow4(Mt) + pow2(Mst1)*(8*pow2(Mst2)*pow4(Mt) - pow4(s2t)*
        pow6(Mst2)) + pow4(s2t)*pow8(Mst2))) + 1800*pow4(Dmsqst2)*pow4(Mst1)*(
        16*Dmglst2*Mgl*Mst2*pow2(Mst1)*(-(Mst2*Mt) + 3*s2t*pow2(Mst1))*pow3(Mt)
        + 16*Mst2*pow2(Dmglst2)*pow2(Mst1)*(-(Mst2*Mt) + 3*s2t*pow2(Mst1))*
        pow3(Mt) + pow2(Mgl)*(8*(5*Mt - 2*Mst2*s2t)*pow3(Mt)*pow4(Mst1) + 4*
        pow4(Mst2)*pow4(Mt) - pow2(Mst1)*(4*pow2(Mst2)*pow4(Mt) + pow4(s2t)*
        pow6(Mst2)) + pow4(s2t)*pow8(Mst2))) + pow8(Msq)*(20*Dmglst2*Mgl*(-4*
        pow4(Mst1)*pow4(Mst2)*(1860*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 1472*Mst2*
        s2t*pow3(Mt) + 774*Mt*pow3(Mst2)*pow3(s2t) - 2464*pow4(Mt) + 15*pow4(
        Mst2)*pow4(s2t)) + 2*pow2(Mst2)*(-432*pow2(Mst2)*pow2(Mt)*pow2(s2t) -
        11656*Mst2*s2t*pow3(Mt) + 996*Mt*pow3(Mst2)*pow3(s2t) + 18748*pow4(Mt)
        + 207*pow4(Mst2)*pow4(s2t))*pow6(Mst1) - 3*pow2(Mst1)*(-856*pow2(Mst2)*
        pow2(Mt)*pow2(s2t) - 1536*Mst2*s2t*pow3(Mt) + 384*Mt*pow3(Mst2)*pow3(
        s2t) + 1200*pow4(Mt) + 203*pow4(Mst2)*pow4(s2t))*pow6(Mst2) + 4*(-738*
        pow2(Mst2)*pow2(Mt)*pow2(s2t) - 10004*Mst2*s2t*pow3(Mt) + 135*Mt*pow3(
        Mst2)*pow3(s2t) + 17532*pow4(Mt) + 12*pow4(Mst2)*pow4(s2t))*pow8(Mst1)
        - 768*pow2(Mt)*(-2*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow8(Mst2)) - 5*
        pow2(Mgl)*(4*pow4(Mst1)*pow4(Mst2)*(3372*pow2(Mst2)*pow2(Mt)*pow2(s2t)
        - 9024*Mst2*s2t*pow3(Mt) + 3912*Mt*pow3(Mst2)*pow3(s2t) + 6112*pow4(Mt)
        + 579*pow4(Mst2)*pow4(s2t)) - 24*pow2(Mst2)*(-534*pow2(Mst2)*pow2(Mt)*
        pow2(s2t) + 1824*Mst2*s2t*pow3(Mt) - 60*Mt*pow3(Mst2)*pow3(s2t) - 934*
        pow4(Mt) + 97*pow4(Mst2)*pow4(s2t))*pow6(Mst1) + 6*pow2(Mst1)*(-728*
        pow2(Mst2)*pow2(Mt)*pow2(s2t) - 1024*Mst2*s2t*pow3(Mt) + 256*Mt*pow3(
        Mst2)*pow3(s2t) + 1200*pow4(Mt) + 139*pow4(Mst2)*pow4(s2t))*pow6(Mst2)
        + 3*(4640*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 15296*Mst2*s2t*pow3(Mt) +
        240*Mt*pow3(Mst2)*pow3(s2t) + 7128*pow4(Mt) - 279*pow4(Mst2)*pow4(s2t))
        *pow8(Mst1) + 768*pow2(Mt)*(-2*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow8(
        Mst2)) + 2*pow2(Dmglst2)*(2*pow4(Mst1)*pow4(Mst2)*(-40440*pow2(Mst2)*
        pow2(Mt)*pow2(s2t) - 75360*Mst2*s2t*pow3(Mt) - 5760*Mt*pow3(Mst2)*pow3(
        s2t) + 11872*pow4(Mt) + 4335*pow4(Mst2)*pow4(s2t)) + 2*pow2(Mst2)*(
        29760*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 31000*Mst2*s2t*pow3(Mt) + 15840*
        Mt*pow3(Mst2)*pow3(s2t) - 228452*pow4(Mt) + 1065*pow4(Mst2)*pow4(s2t))*
        pow6(Mst1) - 15*pow2(Mst1)*(-3080*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 6144*
        Mst2*s2t*pow3(Mt) + 1536*Mt*pow3(Mst2)*pow3(s2t) + 3600*pow4(Mt) + 865*
        pow4(Mst2)*pow4(s2t))*pow6(Mst2) + 40*(-738*pow2(Mst2)*pow2(Mt)*pow2(
        s2t) - 10004*Mst2*s2t*pow3(Mt) + 135*Mt*pow3(Mst2)*pow3(s2t) + 17532*
        pow4(Mt) + 12*pow4(Mst2)*pow4(s2t))*pow8(Mst1) + 480*(-40*pow2(Mst2)*
        pow2(Mt)*pow2(s2t) + 80*pow4(Mt) + 3*pow4(Mst2)*pow4(s2t))*pow8(Mst2)))
        )))/(16200.*pow2(Mgl)*pow4(Mst1)*pow4(Mst2)*pow8(Msq)))/pow4(Mt)
        *12.; 

   return result;
}

/**
 * 	@return returns the susy log^2 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double himalaya::H6bq2g2::calc_coef_at_as2_no_sm_logs_log2() const {

   const double result =
      (((30*log(pow2(Mst2)/pow2(Mst1))*pow4(Mst1)*(pow2(Mgl)*(-4*
        pow4(Mst2)*(14*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 456*Mst2*s2t*pow3(Mt) +
        146*Mt*pow3(Mst2)*pow3(s2t) - 166*pow4(Mt) + 33*pow4(Mst2)*pow4(s2t)) +
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
        pow6(Mst2))))*pow8(Msq) - 10*Dmglst2*Mgl*pow8(Msq)*(pow4(Mst1)*pow4(
        Mst2)*(-8208*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 4640*Mst2*s2t*pow3(Mt) -
        1968*Mt*pow3(Mst2)*pow3(s2t) + 4816*pow4(Mt) + 849*pow4(Mst2)*pow4(s2t)
        ) - 3*pow2(Mst2)*(-472*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 2496*Mst2*s2t*
        pow3(Mt) - 1040*Mt*pow3(Mst2)*pow3(s2t) - 3360*pow4(Mt) + 37*pow4(Mst2)
        *pow4(s2t))*pow6(Mst1) - 3*pow2(Mst1)*(-984*pow2(Mst2)*pow2(Mt)*pow2(
        s2t) - 1536*Mst2*s2t*pow3(Mt) + 384*Mt*pow3(Mst2)*pow3(s2t) + 1456*
        pow4(Mt) + 219*pow4(Mst2)*pow4(s2t))*pow6(Mst2) + 3*(-2240*Mst2*s2t*
        pow3(Mt) + 3744*pow4(Mt) - 59*pow4(Mst2)*pow4(s2t))*pow8(Mst1) - 768*
        pow2(Mt)*(-2*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow8(Mst2)) + pow2(
        Dmglst2)*pow8(Msq)*(pow4(Mst1)*pow4(Mst2)*(115440*pow2(Mst2)*pow2(Mt)*
        pow2(s2t) + 98880*Mst2*s2t*pow3(Mt) + 19680*Mt*pow3(Mst2)*pow3(s2t) -
        29264*pow4(Mt) - 18495*pow4(Mst2)*pow4(s2t)) + 5*pow2(Mst2)*(-2712*
        pow2(Mst2)*pow2(Mt)*pow2(s2t) - 7680*Mst2*s2t*pow3(Mt) - 8544*Mt*pow3(
        Mst2)*pow3(s2t) + 32224*pow4(Mt) + 1101*pow4(Mst2)*pow4(s2t))*pow6(
        Mst1) + 15*pow2(Mst1)*(-3464*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 6144*Mst2*
        s2t*pow3(Mt) + 1536*Mt*pow3(Mst2)*pow3(s2t) + 4368*pow4(Mt) + 913*pow4(
        Mst2)*pow4(s2t))*pow6(Mst2) - 30*(-2240*Mst2*s2t*pow3(Mt) + 3744*pow4(
        Mt) - 59*pow4(Mst2)*pow4(s2t))*pow8(Mst1) - 480*(-40*pow2(Mst2)*pow2(
        Mt)*pow2(s2t) + 80*pow4(Mt) + 3*pow4(Mst2)*pow4(s2t))*pow8(Mst2)) - 5*
        pow2(Mgl)*(720*pow4(Dmsqst2)*pow4(Mst1)*pow4(Mst2)*pow4(Mt) - 720*pow2(
        Dmsqst2)*pow4(Msq)*pow4(Mst1)*pow4(Mst2)*pow4(Mt) - 1440*Dmsqst2*pow4(
        Mst1)*pow4(Mst2)*pow4(Mt)*pow6(Msq) + pow8(Msq)*(pow4(Mst1)*pow4(Mst2)*
        (-8592*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 4800*Mst2*s2t*pow3(Mt) - 3936*
        Mt*pow3(Mst2)*pow3(s2t) - 4256*pow4(Mt) + 69*pow4(Mst2)*pow4(s2t)) + 3*
        pow2(Mst2)*(600*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 128*Mst2*s2t*pow3(Mt) +
        1568*Mt*pow3(Mst2)*pow3(s2t) + 672*pow4(Mt) + 355*pow4(Mst2)*pow4(s2t))
        *pow6(Mst1) - 3*pow2(Mst1)*(-856*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 1024*
        Mst2*s2t*pow3(Mt) + 256*Mt*pow3(Mst2)*pow3(s2t) + 1456*pow4(Mt) + 155*
        pow4(Mst2)*pow4(s2t))*pow6(Mst2) + 3*(-128*Mst2*s2t*pow3(Mt) + 800*
        pow4(Mt) - 239*pow4(Mst2)*pow4(s2t))*pow8(Mst1) - 384*pow2(Mt)*(-2*
        pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow8(Mst2)))))/(540.*pow2(Mgl)*pow4(
        Mst1)*pow4(Mst2)*pow8(Msq)))/pow4(Mt)*12.; 

   return result;
}

/**
 * 	@return returns the susy log^3 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double himalaya::H6bq2g2::calc_coef_at_as2_no_sm_logs_log3() const {

   const double result =
      ((-224*pow4(Mt))/9.)/pow4(Mt)*12.; 

   return result;
}

