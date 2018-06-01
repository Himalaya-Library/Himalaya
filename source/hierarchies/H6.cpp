#include "H6.hpp"
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
 * 	@param lmMt a double log((<renormalization scale> / Mt)^2)
 * 	@param lmMst1 a double log((<renormalization scale> / Mst1)^2)
 * 	@param lmMst2 a double log((<renormalization scale> / Mst2)^2)
 * 	@param lmMsq a double log((<renormalization scale> / Msq)^2)
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
himalaya::H6::H6(const std::map<unsigned int, unsigned int>& flagMap, double Al4p, double beta, double Dmglst2,
		 double lmMt, double lmMst1, double lmMst2, double lmMsq,
		 double Mt, double Mst1, double Mst2, double Msq, double MuSUSY,
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
   xMsq = flagMap.at(ExpansionDepth::xxMsq);
   xMst = flagMap.at(ExpansionDepth::xxMst);
   
   s1 = 
   #include "h6/sigS1Full.inc"
   ;
   s2 = 
   #include "h6/sigS2Full.inc"
   ;
   s12 = 
   #include "h6/sigS12Full.inc"
   ;
}

/**
 * 	@return The diagonal (1, 1) matrix element of the Higgs mass matrix as a double for the hierarchy 'H6'
 */
double himalaya::H6::getS1() const {
   return s1;
}

/**
 * 	@return The diagonal (2, 2) matrix element of the Higgs mass matrix as a double for the hierarchy 'H6'
 */
double himalaya::H6::getS2() const {
   return s2;
}

/**
 * 	@return The off-diagonal (1, 2) = (2, 1) matrix element of the Higgs mass matrix as a double for the hierarchy 'H6'
 */
double himalaya::H6::getS12() const {
   return s12;
}

double himalaya::H6::calc_coef_at_as2_no_sm_logs_log0() const {

   const double result =
      ((-(3*Mst2*pow2(Dmglst2)*(3*pow2(Mst2)*(12*pow2(Mst2)*pow2(Mt)
        *pow2(s2t)*(336*T1ep + 55597*z3 - 264*z4 + 252*pow2(z2)) + 32*Mst2*s2t*
        (224*T1ep - 19363*z3 + 112*z4 + 168*pow2(z2))*pow3(Mt) - 16*Mt*(112*
        T1ep + 36013*z3 + 218*z4 - 3804*pow2(z2))*pow3(Mst2)*pow3(s2t) - 8*(
        2128*T1ep - 154541*z3 + 1064*z4 + 1596*pow2(z2))*pow4(Mt) + (56*T1ep -
        129397*z3 - 6344*z4 - 1254*pow2(z2))*pow4(Mst2)*pow4(s2t)) + 2*pow2(
        Mst1)*(18*pow2(Mst2)*pow2(Mt)*pow2(s2t)*(1584*T1ep - 15053*z3 + 792*z4
        + 1188*pow2(z2)) + 64*Mst2*s2t*(4228*T1ep - 32168*z3 + 2114*z4 + 3171*
        pow2(z2))*pow3(Mt) - 4*Mt*(10948*T1ep - 88880*z3 + 5474*z4 + 8211*pow2(
        z2))*pow3(Mst2)*pow3(s2t) - 16*(22676*T1ep - 463990*z3 + 11338*z4 +
        17007*pow2(z2))*pow4(Mt) + (4544*T1ep - 99253*z3 + 4054*z4 + 5352*pow2(
        z2))*pow4(Mst2)*pow4(s2t))) + 3*Mst2*(18*pow4(Mst2)*(-6*pow2(Mst2)*
        pow2(Mt)*pow2(s2t)*(92*T1ep + 1622*z3 + 142*z4 + 69*pow2(z2)) + 24*
        Mst2*s2t*(28*T1ep - 1922*z3 + 206*z4 + 21*pow2(z2))*pow3(Mt) - 2*Mt*(
        84*T1ep + 10490*z3 + 510*z4 - 1665*pow2(z2))*pow3(Mst2)*pow3(s2t) + 24*
        (28*T1ep + 3902*z3 + 14*z4 + 21*pow2(z2))*pow4(Mt) + 3*(32*T1ep - 1913*
        z3 + 10*z4 - 48*pow2(z2))*pow4(Mst2)*pow4(s2t)) + 6*pow2(Mst1)*pow2(
        Mst2)*(-16*pow2(Mst2)*pow2(Mt)*pow2(s2t)*(644*T1ep - 940*z3 + 214*z4 +
        483*pow2(z2)) + 8*Mst2*s2t*(1508*T1ep - 48454*z3 + 754*z4 + 1131*pow2(
        z2))*pow3(Mt) - 4*Mt*(292*T1ep - 8678*z3 + 740*z4 + 219*pow2(z2))*pow3(
        Mst2)*pow3(s2t) + 8*(836*T1ep + 69050*z3 + 418*z4 + 627*pow2(z2))*pow4(
        Mt) + (1372*T1ep - 22544*z3 + 956*z4 + 1677*pow2(z2))*pow4(Mst2)*pow4(
        s2t)) + pow4(Mst1)*(-4*pow2(Mst2)*pow2(Mt)*pow2(s2t)*(27668*T1ep -
        66166*z3 + 13834*z4 + 20751*pow2(z2)) + 64*Mst2*s2t*(2740*T1ep - 57758*
        z3 + 1370*z4 + 2055*pow2(z2))*pow3(Mt) - 4*Mt*(1756*T1ep - 44630*z3 +
        878*z4 + 1317*pow2(z2))*pow3(Mst2)*pow3(s2t) + 256*(292*T1ep + 24241*z3
        + 146*z4 + 219*pow2(z2))*pow4(Mt) + (7780*T1ep - 34070*z3 + 2594*z4 +
        5835*pow2(z2))*pow4(Mst2)*pow4(s2t))) + 2*Dmglst2*(-18*pow2(Mst1)*pow2(
        Mst2)*(-24*pow2(Mst2)*pow2(Mt)*pow2(s2t)*(200*T1ep + 741*z3 + 100*z4 +
        150*pow2(z2)) + 12*Mst2*s2t*(1844*T1ep - 21022*z3 + 922*z4 + 1383*pow2(
        z2))*pow3(Mt) - 2*Mt*(1604*T1ep - 6250*z3 + 208*z4 + 1203*pow2(z2))*
        pow3(Mst2)*pow3(s2t) - 16*(644*T1ep - 24526*z3 + 322*z4 + 483*pow2(z2))
        *pow4(Mt) + 3*(288*T1ep + 2239*z3 - 54*z4)*pow4(Mst2)*pow4(s2t)) - 27*
        pow4(Mst2)*(16*(-439*z3 + 108*z4)*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 8*
        Mst2*s2t*(140*T1ep + 5742*z3 - 506*z4 + 105*pow2(z2))*pow3(Mt) - 2*Mt*(
        140*T1ep - 37474*z3 - 542*z4 + 5289*pow2(z2))*pow3(Mst2)*pow3(s2t) - 8*
        (112*T1ep + 11897*z3 + 56*z4 + 84*pow2(z2))*pow4(Mt) + (56*T1ep +
        16175*z3 + 424*z4 + 474*pow2(z2))*pow4(Mst2)*pow4(s2t)) + pow4(Mst1)*(
        144*pow2(Mst2)*pow2(Mt)*pow2(s2t)*(3076*T1ep - 734*z3 + 1538*z4 + 2307*
        pow2(z2)) - 32*Mst2*s2t*(51460*T1ep - 614702*z3 + 25730*z4 + 38595*
        pow2(z2))*pow3(Mt) + 2*Mt*(62284*T1ep - 565214*z3 + 31142*z4 + 46713*
        pow2(z2))*pow3(Mst2)*pow3(s2t) + 64*(11164*T1ep - 562226*z3 + 5582*z4 +
        8373*pow2(z2))*pow4(Mt) - (33592*T1ep - 197843*z3 + 16796*z4 + 25194*
        pow2(z2))*pow4(Mst2)*pow4(s2t))) - 1944*z3*log(pow2(Mst2)/pow2(Mst1))*(
        8*(pow2(Mst1) - pow2(Mst2))*pow2(Mt)*pow2(s2t) - 16*pow4(Mt) + (-pow4(
        Mst1) + pow4(Mst2))*pow4(s2t))*pow5(Mst2))/(5832.*pow5(Mst2)) - pow4(
        s2t)*(pow2(Mst1)*((5*(1 + 2*log(pow2(Msq)/pow2(Mst1)))*(1 - 2*log(pow2(
        Mst2)/pow2(Mst1)))*pow2(Msq))/6. + pow2(Dmglst2)*(104.7055157750343 -
        B4/3. + (4*D3)/9. - (5*DN)/18. + (361561*log(pow2(Mst2)/pow2(Mst1)))/
        8100. - (5*(-1 + 2*log(pow2(Mst2)/pow2(Mst1)))*pow2(log(pow2(Msq)/pow2(
        Mst1))))/4. - (62*pow2(log(pow2(Mst2)/pow2(Mst1))))/15. + (5*log(pow2(
        Msq)/pow2(Mst1))*(-167 + 78*log(pow2(Mst2)/pow2(Mst1)) + 36*pow2(log(
        pow2(Mst2)/pow2(Mst1)))))/36. - (19*pow3(log(pow2(Mst2)/pow2(Mst1))))/
        18.)) - ((5*(6*Dmglst2 + 5*Mst2 + Mst2*log(pow2(Mst2)/pow2(Mst1)))*
        pow2(log(pow2(Msq)/pow2(Mst1))))/(12.*Mst2) + (5*log(pow2(Msq)/pow2(
        Mst1))*(131 + (48*Dmglst2)/Mst2 + 52*log(pow2(Mst2)/pow2(Mst1)) + (24*
        pow2(Msq))/pow2(Mst2) + 12*pow2(log(pow2(Mst2)/pow2(Mst1)))))/72. + (
        Mst2*(506515899052*Dmglst2 + 3*(71378675209 + 2000376000*B4 +
        222264000*D3 - 222264000*DN)*Mst2) - 1260*(4045612*Dmglst2 - 160083513*
        Mst2)*Mst2*log(pow2(Mst2)/pow2(Mst1)) + 5000940000*pow2(Msq) - 1587600*
        (12204*Dmglst2 - 28411*Mst2)*Mst2*pow2(log(pow2(Mst2)/pow2(Mst1))) -
        222264000*Mst2*(33*Dmglst2 + 136*Mst2)*pow3(log(pow2(Mst2)/pow2(Mst1)))
        )/(6.001128e9*pow2(Mst2)))*pow4(Mst1) + (pow2(Mst2)*(-4900195*pow2(
        Dmglst2)*pow2(Msq) + 205200*B4*pow2(Dmglst2)*pow2(Msq) - 216000*D3*
        pow2(Dmglst2)*pow2(Msq) + 113400*DN*pow2(Dmglst2)*pow2(Msq) + 1845000*
        pow2(Dmglst2)*pow2(Mst1) + 9250423*pow2(Msq)*pow2(Mst1) - 10800*B4*
        pow2(Msq)*pow2(Mst1) + 21600*D3*pow2(Msq)*pow2(Mst1) - 16200*DN*pow2(
        Msq)*pow2(Mst1) + 40500*pow2(Msq)*(3*pow2(Dmglst2) + 6*log(pow2(Mst2)/
        pow2(Mst1))*(pow2(Dmglst2) - pow2(Mst1)) + 7*pow2(Mst1))*pow2(log(pow2(
        Msq)/pow2(Mst1))) + 1800*pow2(Msq)*(456*pow2(Dmglst2) + 185*pow2(Mst1))
        *pow3(log(pow2(Mst2)/pow2(Mst1))) + 81000*pow4(Msq) - 265770*pow4(Mst1)
        + 900*pow2(log(pow2(Mst2)/pow2(Mst1)))*(9*pow2(Dmglst2)*(353*pow2(Msq)
        - 200*pow2(Mst1)) + 289*pow2(Msq)*pow2(Mst1) + 108*pow4(Mst1)) + 60*
        log(pow2(Mst2)/pow2(Mst1))*(-36998*pow2(Msq)*pow2(Mst1) + pow2(Dmglst2)
        *(-70215*pow2(Msq) + 15300*pow2(Mst1)) + 2700*pow4(Msq) + 954*pow4(
        Mst1)) - 2700*log(pow2(Msq)/pow2(Mst1))*(-5*pow2(Dmglst2)*(31*pow2(Msq)
        - 36*pow2(Mst1)) - 445*pow2(Msq)*pow2(Mst1) + 60*pow2(Msq)*(3*pow2(
        Dmglst2) - 2*pow2(Mst1))*pow2(log(pow2(Mst2)/pow2(Mst1))) - 60*pow4(
        Msq) + 34*pow4(Mst1) + 6*log(pow2(Mst2)/pow2(Mst1))*(5*pow2(Dmglst2)*(
        17*pow2(Msq) - 20*pow2(Mst1)) + 35*pow2(Msq)*pow2(Mst1) - 20*pow4(Msq)
        + 2*pow4(Mst1)))))/(97200.*pow2(Msq)) + (-3*Mst2*pow2(Dmglst2)*(4*(
        11360*OepS2 - 378621*S2)*pow2(Mst1) + 3*(280*OepS2 - 5331339*S2)*pow2(
        Mst2)) - 15*Mst2*(3*(2744*OepS2 - 254097*S2)*pow2(Mst1)*pow2(Mst2) + (
        7780*OepS2 - 512865*S2)*pow4(Mst1) + 864*(2*OepS2 - 891*S2)*pow4(Mst2))
        + 10*Dmglst2*(972*(16*OepS2 + 275*S2)*pow2(Mst1)*pow2(Mst2) + (33592*
        OepS2 - 1394145*S2)*pow4(Mst1) + 27*(56*OepS2 + 135837*S2)*pow4(Mst2))
        + 405*S2*log(pow2(Mst2)/pow2(Mst1))*(6174*pow2(Mst1)*pow3(Mst2) + 6*
        pow2(Dmglst2)*(1136*Mst2*pow2(Mst1) + 21*pow3(Mst2)) + 5835*Mst2*pow4(
        Mst1) - 4*Dmglst2*(1944*pow2(Mst1)*pow2(Mst2) + 4199*pow4(Mst1) + 189*
        pow4(Mst2)) + 1296*pow5(Mst2)))/(43740.*Mst2) + (Dmglst2*Mst2*(-
        5502750*pow2(Mst1)*pow2(Mst2)*pow4(Msq) + 324000*B4*pow2(Mst1)*pow2(
        Mst2)*pow4(Msq) - 432000*D3*pow2(Mst1)*pow2(Mst2)*pow4(Msq) + 270000*
        DN*pow2(Mst1)*pow2(Mst2)*pow4(Msq) - 18000*pow2(Mst1)*(77*pow2(Mst1) -
        456*pow2(Mst2))*pow3(log(pow2(Mst2)/pow2(Mst1)))*pow4(Msq) + 6930000*
        pow2(Msq)*pow2(Mst2)*pow4(Mst1) - 6909836*pow4(Msq)*pow4(Mst1) -
        324000*B4*pow4(Msq)*pow4(Mst1) + 432000*D3*pow4(Msq)*pow4(Mst1) -
        270000*DN*pow4(Msq)*pow4(Mst1) - 10350000*pow2(Msq)*pow2(Mst1)*pow4(
        Mst2) + 1552500*pow4(Msq)*pow4(Mst2) + 4460625*pow4(Mst1)*pow4(Mst2) -
        1215000*pow2(log(pow2(Msq)/pow2(Mst1)))*pow4(Msq)*(-(pow2(Mst1)*pow2(
        Mst2)) - pow4(Mst1) + 2*log(pow2(Mst2)/pow2(Mst1))*(-(pow2(Mst1)*pow2(
        Mst2)) + pow4(Mst1)) + pow4(Mst2)) + 562500*pow2(Msq)*pow6(Mst1) -
        635625*pow2(Mst2)*pow6(Mst1) + 30*log(pow2(Mst2)/pow2(Mst1))*(375*pow2(
        Mst2)*(-17*pow2(Mst1) + 3*pow2(Mst2))*pow4(Mst1) + 4*pow4(Msq)*(-
        155025*pow2(Mst1)*pow2(Mst2) + 128893*pow4(Mst1) + 55575*pow4(Mst2)) -
        1500*pow2(Msq)*(-100*pow2(Mst2)*pow4(Mst1) + 13*pow2(Mst1)*pow4(Mst2) +
        24*pow6(Mst1))) + 1800*pow2(log(pow2(Mst2)/pow2(Mst1)))*(pow4(Msq)*(-
        3585*pow2(Mst1)*pow2(Mst2) + 5644*pow4(Mst1) - 3480*pow4(Mst2)) + 150*
        pow2(Msq)*(-20*pow2(Mst2)*pow4(Mst1) + 17*pow2(Mst1)*pow4(Mst2) + 3*
        pow6(Mst1)) + 75*(-21*pow4(Mst1)*pow4(Mst2) + 5*pow2(Mst2)*pow6(Mst1)))
        + 67500*log(pow2(Msq)/pow2(Mst1))*(72*pow2(Mst1)*(pow2(Mst1) - pow2(
        Mst2))*pow2(log(pow2(Mst2)/pow2(Mst1)))*pow4(Msq) + pow2(Mst2)*(pow2(
        Mst1) + 19*pow2(Mst2))*pow4(Mst1) + 6*pow4(Msq)*(15*pow2(Mst1)*pow2(
        Mst2) - 9*pow4(Mst1) + pow4(Mst2)) + 8*pow2(Msq)*(pow2(Mst2)*pow4(Mst1)
        - 8*pow2(Mst1)*pow4(Mst2) + 2*pow6(Mst1)) - 2*log(pow2(Mst2)/pow2(Mst1)
        )*(18*pow4(Msq)*(pow2(Mst1)*pow2(Mst2) + pow4(Mst1) - pow4(Mst2)) - 21*
        pow4(Mst1)*pow4(Mst2) + 5*pow2(Mst2)*pow6(Mst1) + pow2(Msq)*(-40*pow2(
        Mst2)*pow4(Mst1) + 34*pow2(Mst1)*pow4(Mst2) + 6*pow6(Mst1))))))/(
        486000.*pow2(Mst1)*pow4(Msq)) + (pow4(Mst2)*(49392000*pow2(Dmglst2)*
        pow2(Msq)*pow2(Mst1)*pow2(Mst2) + 362722500*pow2(Dmglst2)*pow2(Mst1)*
        pow4(Msq) - 39513600*pow2(Dmglst2)*pow2(Mst2)*pow4(Msq) - 125023500*
        pow2(Mst1)*pow2(Mst2)*pow4(Msq) - 559776000*pow2(Dmglst2)*pow2(Msq)*
        pow4(Mst1) - 422790375*pow2(Dmglst2)*pow2(Mst2)*pow4(Mst1) - 122903760*
        pow2(Msq)*pow2(Mst2)*pow4(Mst1) - 322522900*pow4(Msq)*pow4(Mst1) +
        24696000*B4*pow4(Msq)*pow4(Mst1) - 2469600*D3*pow4(Msq)*pow4(Mst1) +
        1234800*DN*pow4(Msq)*pow4(Mst1) + 232142400*pow3(log(pow2(Mst2)/pow2(
        Mst1)))*pow4(Msq)*pow4(Mst1) - 18522000*pow2(Mst1)*pow6(Msq) + 420*log(
        pow2(Msq)/pow2(Mst1))*pow2(Mst1)*(-126420*pow2(Msq)*pow2(Mst1)*pow2(
        Mst2) - 51450*pow2(Mst1)*pow4(Msq) - 110250*pow2(Mst2)*pow4(Msq) -
        308700*pow2(Mst1)*pow2(log(pow2(Mst2)/pow2(Mst1)))*pow4(Msq) + 11025*
        pow2(Dmglst2)*(5*pow2(Mst1)*(pow2(Mst1) - 11*pow2(Mst2)) - 8*pow2(Msq)*
        (6*pow2(Mst1) - 5*pow2(Mst2)) + 26*pow4(Msq)) + 67620*pow2(Msq)*pow4(
        Mst1) + 51313*pow2(Mst2)*pow4(Mst1) - 88200*pow6(Msq) + 105*log(pow2(
        Mst2)/pow2(Mst1))*(70*(59*pow2(Mst1) + 18*pow2(Mst2))*pow4(Msq) + 422*
        pow2(Mst2)*pow4(Mst1) + 210*pow2(Dmglst2)*(-38*pow2(Msq)*pow2(Mst1) -
        30*pow2(Mst1)*pow2(Mst2) + 6*pow4(Msq) + 35*pow4(Mst1)) + 28*pow2(Msq)*
        (-40*pow2(Mst1)*pow2(Mst2) + 41*pow4(Mst1)) - 187*pow6(Mst1)) - 27672*
        pow6(Mst1)) + 433852125*pow2(Dmglst2)*pow6(Mst1) + 108497760*pow2(Msq)*
        pow6(Mst1) + 45671209*pow2(Mst2)*pow6(Mst1) + 88200*pow2(Mst1)*pow2(
        log(pow2(Msq)/pow2(Mst1)))*(-315*pow2(Dmglst2)*pow4(Msq) + 105*pow2(
        Mst1)*pow4(Msq) + 735*log(pow2(Mst2)/pow2(Mst1))*pow2(Mst1)*pow4(Msq) -
        315*pow2(Mst2)*pow4(Msq) - 4*pow2(Mst2)*pow4(Mst1) + 6*pow6(Mst1)) -
        105*log(pow2(Mst2)/pow2(Mst1))*(-6860*pow4(Msq)*(330*pow2(Mst1)*pow2(
        Mst2) + 911*pow4(Mst1)) + 1470*pow2(Dmglst2)*(4*(297*pow2(Mst1) - 32*
        pow2(Mst2))*pow4(Msq) + 5*(3*pow2(Mst1) - 274*pow2(Mst2))*pow4(Mst1) -
        20*pow2(Msq)*(-60*pow2(Mst1)*pow2(Mst2) + 29*pow4(Mst1))) + 90494*pow2(
        Mst2)*pow6(Mst1) - 392*pow2(Msq)*(-735*pow2(Mst2)*pow4(Mst1) + 692*
        pow6(Mst1)) - 34519*pow8(Mst1)) - 44100*pow2(log(pow2(Mst2)/pow2(Mst1))
        )*(56*pow4(Msq)*(42*pow2(Mst1)*pow2(Mst2) + 293*pow4(Mst1)) + 414*pow2(
        Mst2)*pow6(Mst1) + 56*pow2(Msq)*(-20*pow2(Mst2)*pow4(Mst1) + 19*pow6(
        Mst1)) + 14*pow2(Dmglst2)*(8*(61*pow2(Mst1) - 12*pow2(Mst2))*pow4(Msq)
        - 570*pow2(Msq)*pow4(Mst1) + 75*(-6*pow2(Mst2)*pow4(Mst1) + 7*pow6(
        Mst1))) - 223*pow8(Mst1)) - 18222471*pow8(Mst1)))/(2.22264e7*pow4(Msq)*
        pow4(Mst1))) + pow2(Mt)*pow2(s2t)*((40*(1 + 2*log(pow2(Msq)/pow2(Mst1))
        )*pow2(Msq))/3. + pow2(Dmglst2)*(1007.0386243386243 - (16*B4)/9. + (16*
        D3)/9. - (8*DN)/9. - (10540*log(pow2(Mst2)/pow2(Mst1)))/27. - 10*log(
        pow2(Msq)/pow2(Mst1))*(9 + 4*log(pow2(Mst2)/pow2(Mst1))) + 20*pow2(log(
        pow2(Msq)/pow2(Mst1))) - (26*pow2(log(pow2(Mst2)/pow2(Mst1))))/3. - (
        128*pow3(log(pow2(Mst2)/pow2(Mst1))))/9.) + pow2(Mst1)*((2*(517702 -
        2700*D3 + 2700*DN + log(pow2(Mst2)/pow2(Mst1))*(-719895 - (202500*pow2(
        Dmglst2))/pow2(Msq)) + (64125*pow2(Dmglst2))/pow2(Msq) - (3375*log(
        pow2(Msq)/pow2(Mst1))*(-60*pow2(Dmglst2) + 3*pow2(Msq) + 8*log(pow2(
        Mst2)/pow2(Mst1))*pow2(Msq)))/pow2(Msq) - 10125*(3 + log(pow2(Mst2)/
        pow2(Mst1)))*pow2(log(pow2(Msq)/pow2(Mst1))) + 268875*pow2(log(pow2(
        Mst2)/pow2(Mst1))) - 121050*pow3(log(pow2(Mst2)/pow2(Mst1)))))/6075. +
        ((24348041*pow2(Dmglst2))/141750. - (20*pow2(Msq))/3. + (12287223581*
        pow2(Mst1))/6.251175e7 - (log(pow2(Mst2)/pow2(Mst1))*(28069748*pow2(
        Dmglst2) + 20192173*pow2(Mst1)))/66150. + (40*log(pow2(Msq)/pow2(Mst1))
        *(9*pow2(Dmglst2) - pow2(Msq) - log(pow2(Mst2)/pow2(Mst1))*pow2(Mst1)))
        /3. - 10*pow2(Dmglst2)*pow2(log(pow2(Msq)/pow2(Mst1))) + ((-3566*pow2(
        Dmglst2))/15. + (83647*pow2(Mst1))/945.)*pow2(log(pow2(Mst2)/pow2(Mst1)
        )) - (2*(570*pow2(Dmglst2) + 337*pow2(Mst1))*pow3(log(pow2(Mst2)/pow2(
        Mst1))))/27.)/pow2(Mst2)) + pow2(Mst2)*(252.5348765432099 + (8*D3)/9. -
        (8*DN)/9. + (100*log(pow2(Msq)/pow2(Mst1)))/3. + 20*pow2(log(pow2(Msq)/
        pow2(Mst1))) + (2*log(pow2(Mst2)/pow2(Mst1))*(-2717 - 240*log(pow2(Msq)
        /pow2(Mst1)) + 45*pow2(log(pow2(Msq)/pow2(Mst1)))))/27. + (2*(256 - 30*
        log(pow2(Msq)/pow2(Mst1)))*pow2(log(pow2(Mst2)/pow2(Mst1))))/9. + (16*
        pow3(log(pow2(Mst2)/pow2(Mst1))))/9. - (pow2(Mst1)*(12733875*pow2(
        Dmglst2) + 8481704*pow2(Msq) + 20580*log(pow2(Mst2)/pow2(Mst1))*(8325*
        pow2(Dmglst2) + 1364*pow2(Msq) - 339*pow2(Mst1)) + 6570120*pow2(Mst1) -
        420*log(pow2(Msq)/pow2(Mst1))*(407925*pow2(Dmglst2) + 30772*pow2(Msq) -
        14865*pow2(Mst1) + 735*log(pow2(Mst2)/pow2(Mst1))*(16*pow2(Msq) + 9*
        pow2(Mst1))) + 44100*(112*pow2(Msq) + 15*pow2(Mst1))*pow2(log(pow2(Msq)
        /pow2(Mst1))) - 3704400*pow2(Msq)*pow2(log(pow2(Mst2)/pow2(Mst1)))))/(
        2.7783e6*pow4(Msq)) + (pow2(Dmglst2)*(-30*(11 + 40*log(pow2(Msq)/pow2(
        Mst1)) - 40*log(pow2(Mst2)/pow2(Mst1)))*pow2(Mst1) + pow2(Msq)*(727 -
        594*log(pow2(Mst2)/pow2(Mst1)) + 30*log(pow2(Msq)/pow2(Mst1))*(13 + 6*
        log(pow2(Mst2)/pow2(Mst1))) - 90*pow2(log(pow2(Msq)/pow2(Mst1))) - 400*
        pow2(log(pow2(Mst2)/pow2(Mst1))))) - 60*(1 + 2*log(pow2(Msq)/pow2(Mst1)
        ))*pow4(Msq))/(9.*pow2(Msq)*pow2(Mst1))) + (-((-20580*(5108 + 7500*log(
        pow2(Mst2)/pow2(Mst1)) + 15*log(pow2(Msq)/pow2(Mst1))*(218 + 75*log(
        pow2(Mst2)/pow2(Mst1))) + 900*pow2(log(pow2(Msq)/pow2(Mst1))))*pow3(
        Mst2) + Dmglst2*pow2(Msq)*(3966427397 - 6573806820*log(pow2(Mst2)/pow2(
        Mst1)) - 277830000*log(pow2(Msq)/pow2(Mst1))*(-4 + 5*log(pow2(Mst2)/
        pow2(Mst1))) - 3242408400*pow2(log(pow2(Mst2)/pow2(Mst1))) +
        2234988000*pow3(log(pow2(Mst2)/pow2(Mst1)))))*pow4(Mst1))/(2.083725e7*
        pow2(Msq)) + (4*OepS2*(108*pow2(Dmglst2)*(33*Mst2*pow2(Mst1) + 7*pow3(
        Mst2)) + 24*Dmglst2*(150*pow2(Mst1)*pow2(Mst2) + 769*pow4(Mst1)) -
        Mst2*(3864*pow2(Mst1)*pow2(Mst2) + 6917*pow4(Mst1) + 621*pow4(Mst2))))/
        729. - (S2*(36*pow2(Dmglst2)*(39779*Mst2*pow2(Mst1) + 4707*pow3(Mst2))
        + 168*Dmglst2*(558*pow2(Mst1)*pow2(Mst2) + 28283*pow4(Mst1) + 5697*
        pow4(Mst2)) + 210*log(pow2(Mst2)/pow2(Mst1))*(108*pow2(Dmglst2)*(33*
        Mst2*pow2(Mst1) + 7*pow3(Mst2)) + 24*Dmglst2*(150*pow2(Mst1)*pow2(Mst2)
        + 769*pow4(Mst1)) - Mst2*(3864*pow2(Mst1)*pow2(Mst2) + 6917*pow4(Mst1)
        + 621*pow4(Mst2))) - 7*(172188*pow2(Mst1)*pow3(Mst2) + 406627*Mst2*
        pow4(Mst1) - 72171*pow5(Mst2))))/1890.)/pow3(Mst2) + (17517490200*
        Dmglst2*pow2(Mst2)*pow4(Msq)*pow4(Mst1) - 148176000*B4*Dmglst2*pow2(
        Mst2)*pow4(Msq)*pow4(Mst1) + 148176000*D3*Dmglst2*pow2(Mst2)*pow4(Msq)*
        pow4(Mst1) - 74088000*Dmglst2*DN*pow2(Mst2)*pow4(Msq)*pow4(Mst1) -
        12348000*Dmglst2*(289*pow2(Mst1) + 96*pow2(Mst2))*pow3(log(pow2(Mst2)/
        pow2(Mst1)))*pow4(Msq)*pow4(Mst1) + 1657719000*Dmglst2*pow2(Mst1)*pow4(
        Msq)*pow4(Mst2) - 3303090000*Dmglst2*pow2(Msq)*pow4(Mst1)*pow4(Mst2) +
        740880000*pow2(Dmglst2)*pow2(Msq)*pow2(Mst1)*pow5(Mst2) - 592704000*
        pow2(Dmglst2)*pow4(Msq)*pow5(Mst2) - 1430824500*pow2(Mst1)*pow4(Msq)*
        pow5(Mst2) - 144703125*pow2(Dmglst2)*pow4(Mst1)*pow5(Mst2) - 931944720*
        pow2(Msq)*pow4(Mst1)*pow5(Mst2) - 661500*pow2(Mst1)*pow2(log(pow2(Msq)/
        pow2(Mst1)))*(1260*Dmglst2*pow2(pow2(Mst1) - pow2(Mst2))*pow4(Msq) + (-
        56*pow2(Msq)*pow2(Mst1) + 630*pow4(Msq) + 15*(-(pow2(Mst1)*pow2(Mst2))
        + pow4(Mst1)))*pow5(Mst2)) - 6300*log(pow2(Msq)/pow2(Mst1))*pow2(Mst1)*
        (7350*Dmglst2*(-40*pow2(Msq)*pow2(-(Mst2*pow2(Mst1)) + pow3(Mst2)) + 2*
        pow4(Msq)*(5*pow2(Mst1)*pow2(Mst2) + 16*pow4(Mst1) - 3*pow4(Mst2)) + 3*
        pow2(Mst1)*pow2(Mst2)*(-9*pow2(Mst1)*pow2(Mst2) + pow4(Mst1) + 15*pow4(
        Mst2))) - 11025*pow2(Dmglst2)*(40*pow2(Msq) - 71*pow2(Mst1) + 35*pow2(
        Mst2))*pow5(Mst2) + (196*pow2(Msq)*(641*pow2(Mst1) - 375*pow2(Mst2)) +
        58635*pow2(Mst1)*pow2(Mst2) + 110250*pow4(Msq) - 43935*pow4(Mst1))*
        pow5(Mst2) - 105*log(pow2(Mst2)/pow2(Mst1))*(840*Dmglst2*pow4(Msq)*(-6*
        pow2(Mst1)*pow2(Mst2) + 4*pow4(Mst1) + 3*pow4(Mst2)) - (182*pow2(Msq)*
        pow2(Mst1) + 37*pow2(Mst1)*pow2(Mst2) - 1260*pow4(Msq) + 33*pow4(Mst1))
        *pow5(Mst2))) + 1389150000*Dmglst2*pow2(Msq)*pow2(Mst2)*pow6(Mst1) +
        6631537304*Dmglst2*pow4(Msq)*pow6(Mst1) + 250818750*Dmglst2*pow4(Mst2)*
        pow6(Mst1) + 117538200*pow5(Mst2)*pow6(Mst1) + 1728720000*Dmglst2*pow2(
        Msq)*pow2(Mst1)*pow6(Mst2) - 783326250*Dmglst2*pow4(Mst1)*pow6(Mst2) +
        88200*pow2(log(pow2(Mst2)/pow2(Mst1)))*(10080*pow2(Dmglst2)*pow4(Msq)*
        pow5(Mst2) + 15*pow5(Mst2)*(-56*(15*pow2(Mst1) - 2*pow2(Mst2))*pow4(
        Msq) + 63*pow2(Msq)*pow4(Mst1) + 11*pow2(Mst2)*pow4(Mst1) + 24*pow6(
        Mst1)) + 14*Dmglst2*pow4(Msq)*(5865*pow2(Mst2)*pow4(Mst1) - 2040*pow2(
        Mst1)*pow4(Mst2) + 499*pow6(Mst1) + 480*pow6(Mst2))) + 156279375*pow2(
        Dmglst2)*pow2(Mst1)*pow7(Mst2) + 663705000*pow2(Msq)*pow2(Mst1)*pow7(
        Mst2) + 148176000*pow4(Msq)*pow7(Mst2) - 217865700*pow4(Mst1)*pow7(
        Mst2) + 740880000*Dmglst2*pow2(Msq)*pow8(Mst1) + 34728750*Dmglst2*pow2(
        Mst2)*pow8(Mst1) + 420*log(pow2(Mst2)/pow2(Mst1))*(11025*pow2(Dmglst2)*
        (-600*pow2(Msq)*pow2(Mst1) - 525*pow2(Mst1)*pow2(Mst2) + 64*pow4(Msq) +
        1065*pow4(Mst1))*pow5(Mst2) + 30*pow5(Mst2)*(735*(289*pow2(Mst1) - 32*
        pow2(Mst2))*pow4(Msq) + 28558*pow2(Mst2)*pow4(Mst1) + 98*pow2(Msq)*(-
        375*pow2(Mst1)*pow2(Mst2) + 541*pow4(Mst1)) - 26598*pow6(Mst1)) - 98*
        Dmglst2*(-3375*pow2(Mst2)*pow4(Mst1)*(-9*pow2(Mst1)*pow2(Mst2) + pow4(
        Mst1) + 15*pow4(Mst2)) + 2*pow4(Msq)*(251400*pow2(Mst2)*pow4(Mst1) -
        26775*pow2(Mst1)*pow4(Mst2) + 59657*pow6(Mst1) + 7200*pow6(Mst2)) +
        4500*pow2(Msq)*(-20*pow4(Mst1)*pow4(Mst2) + 10*pow2(Mst2)*pow6(Mst1) +
        10*pow2(Mst1)*pow6(Mst2) + pow8(Mst1)))))/(4.16745e7*Mst2*pow4(Msq)*
        pow4(Mst1))) - (Mt*pow3(s2t)*(583200*(Dmglst2 + Mst2)*pow2(Mst1)*pow2(
        Mst2)*(2*(pow2(Mst1) - pow2(Mst2)) + log(pow2(Mst2)/pow2(Mst1))*(pow2(
        Mst1) + pow2(Mst2)))*pow2(log(pow2(Msq)/pow2(Mst1)))*pow4(Msq) +
        13424256*pow2(Dmglst2)*pow2(Mst1)*pow3(Mst2)*pow4(Msq) - 4432320*B4*
        pow2(Dmglst2)*pow2(Mst1)*pow3(Mst2)*pow4(Msq) + 116640*DN*pow2(Dmglst2)
        *pow2(Mst1)*pow3(Mst2)*pow4(Msq) + 161280*OepS2*pow2(Dmglst2)*pow2(
        Mst1)*pow3(Mst2)*pow4(Msq) + 8007336*S2*pow2(Dmglst2)*pow2(Mst1)*pow3(
        Mst2)*pow4(Msq) - 2268000*pow2(Dmglst2)*pow2(Msq)*pow3(Mst2)*pow4(Mst1)
        - 55534902*Mst2*pow2(Dmglst2)*pow4(Msq)*pow4(Mst1) + 2627520*Mst2*
        OepS2*pow2(Dmglst2)*pow4(Msq)*pow4(Mst1) - 21888792*Mst2*S2*pow2(
        Dmglst2)*pow4(Msq)*pow4(Mst1) + 4175802*Dmglst2*pow2(Mst2)*pow4(Msq)*
        pow4(Mst1) - 660960*B4*Dmglst2*pow2(Mst2)*pow4(Msq)*pow4(Mst1) - 19440*
        Dmglst2*DN*pow2(Mst2)*pow4(Msq)*pow4(Mst1) - 1154880*Dmglst2*OepS2*
        pow2(Mst2)*pow4(Msq)*pow4(Mst1) + 27280152*Dmglst2*S2*pow2(Mst2)*pow4(
        Msq)*pow4(Mst1) - 1677690*pow3(Mst2)*pow4(Msq)*pow4(Mst1) - 660960*B4*
        pow3(Mst2)*pow4(Msq)*pow4(Mst1) - 19440*DN*pow3(Mst2)*pow4(Msq)*pow4(
        Mst1) + 210240*OepS2*pow3(Mst2)*pow4(Msq)*pow4(Mst1) - 10118520*S2*
        pow3(Mst2)*pow4(Msq)*pow4(Mst1) - 22301919*Dmglst2*pow2(Mst1)*pow4(Msq)
        *pow4(Mst2) - 5093280*B4*Dmglst2*pow2(Mst1)*pow4(Msq)*pow4(Mst2) +
        97200*Dmglst2*DN*pow2(Mst1)*pow4(Msq)*pow4(Mst2) - 151200*Dmglst2*
        OepS2*pow2(Mst1)*pow4(Msq)*pow4(Mst2) - 3962844*Dmglst2*S2*pow2(Mst1)*
        pow4(Msq)*pow4(Mst2) - 6998400*Dmglst2*pow2(Msq)*pow4(Mst1)*pow4(Mst2)
        + 32983200*pow2(Dmglst2)*pow2(Msq)*pow2(Mst1)*pow5(Mst2) + 1866240*
        pow2(Dmglst2)*pow4(Msq)*pow5(Mst2) - 20870055*pow2(Mst1)*pow4(Msq)*
        pow5(Mst2) - 2138400*B4*pow2(Mst1)*pow4(Msq)*pow5(Mst2) + 19440*DN*
        pow2(Mst1)*pow4(Msq)*pow5(Mst2) + 90720*OepS2*pow2(Mst1)*pow4(Msq)*
        pow5(Mst2) + 131220*S2*pow2(Mst1)*pow4(Msq)*pow5(Mst2) - 3871800*pow2(
        Dmglst2)*pow4(Mst1)*pow5(Mst2) - 3240000*pow2(Msq)*pow4(Mst1)*pow5(
        Mst2) + 103680*pow2(Mst1)*pow3(log(pow2(Mst2)/pow2(Mst1)))*pow4(Msq)*(
        36*pow2(Dmglst2)*pow3(Mst2) + 3*Dmglst2*(10*pow2(Mst1)*pow2(Mst2) + 7*
        pow4(Mst1) + 18*pow4(Mst2)) + 5*(6*pow2(Mst1)*pow3(Mst2) + 5*Mst2*pow4(
        Mst1) + 6*pow5(Mst2))) + 194400*Dmglst2*pow2(Msq)*pow2(Mst2)*pow6(Mst1)
        - 194400*pow2(Msq)*pow3(Mst2)*pow6(Mst1) + 35542243*Dmglst2*pow4(Msq)*
        pow6(Mst1) - 8134005*Mst2*pow4(Msq)*pow6(Mst1) - 2491360*Dmglst2*OepS2*
        pow4(Msq)*pow6(Mst1) + 210720*Mst2*OepS2*pow4(Msq)*pow6(Mst1) +
        90290268*Dmglst2*S2*pow4(Msq)*pow6(Mst1) - 9916020*Mst2*S2*pow4(Msq)*
        pow6(Mst1) - 93150*Dmglst2*pow4(Mst2)*pow6(Mst1) + 4050*pow5(Mst2)*
        pow6(Mst1) + 18468000*Dmglst2*pow2(Msq)*pow2(Mst1)*pow6(Mst2) + 622080*
        Dmglst2*pow4(Msq)*pow6(Mst2) - 3333150*Dmglst2*pow4(Mst1)*pow6(Mst2) +
        21173400*pow2(Dmglst2)*pow2(Mst1)*pow7(Mst2) + 4600800*pow2(Msq)*pow2(
        Mst1)*pow7(Mst2) - 622080*pow4(Msq)*pow7(Mst2) - 903150*pow4(Mst1)*
        pow7(Mst2) + 6480*pow2(log(pow2(Mst2)/pow2(Mst1)))*(120*pow2(Msq)*pow2(
        Mst1)*(pow2(Mst1) + pow2(Mst2))*pow5(Mst2) + 6*pow2(Dmglst2)*(60*pow2(
        Msq)*pow2(Mst1)*(pow2(Mst1) + pow2(Mst2))*pow3(Mst2) + 2*pow4(Msq)*(92*
        pow2(Mst1)*pow3(Mst2) + 67*Mst2*pow4(Mst1) - 24*pow5(Mst2)) + 75*pow2(
        Mst1)*(pow2(Mst1) + pow2(Mst2))*pow5(Mst2)) + 3*Dmglst2*(120*pow2(Msq)*
        pow2(Mst1)*(pow2(Mst1) + pow2(Mst2))*pow4(Mst2) + pow4(Msq)*(210*pow2(
        Mst2)*pow4(Mst1) - 338*pow2(Mst1)*pow4(Mst2) + 74*pow6(Mst1) - 96*pow6(
        Mst2)) + 75*pow2(Mst1)*(pow2(Mst1) + pow2(Mst2))*pow6(Mst2)) + 45*pow4(
        Mst1)*pow7(Mst2) - 2*pow4(Msq)*(321*pow3(Mst2)*pow4(Mst1) + 843*pow2(
        Mst1)*pow5(Mst2) + 205*Mst2*pow6(Mst1) + 48*pow7(Mst2))) - 16200*log(
        pow2(Msq)/pow2(Mst1))*pow2(Mst1)*(72*(Dmglst2 + Mst2)*pow2(Mst2)*(pow2(
        Mst1) + pow2(Mst2))*pow2(log(pow2(Mst2)/pow2(Mst1)))*pow4(Msq) + 12*
        pow2(Dmglst2)*(24*pow2(Msq)*(pow2(Mst1) - pow2(Mst2))*pow3(Mst2) + (7*
        Mst2*pow2(Mst1) + 107*pow3(Mst2))*pow4(Msq) + 30*(pow2(Mst1) - pow2(
        Mst2))*pow5(Mst2)) - 3*(6*Mst2*pow4(Msq)*(12*pow2(Mst1)*pow2(Mst2) +
        pow4(Mst1) - 24*pow4(Mst2)) + pow2(Mst1)*(pow2(Mst1) - 13*pow2(Mst2))*
        pow5(Mst2) + 32*pow2(Msq)*(-pow2(Mst1) + pow2(Mst2))*pow5(Mst2)) + 6*
        log(pow2(Mst2)/pow2(Mst1))*(6*Mst2*pow4(Msq)*(3*pow2(Mst1)*pow2(Mst2) +
        pow4(Mst1) - 5*pow4(Mst2)) + 6*Mst2*pow2(Dmglst2)*(pow2(Mst1) + pow2(
        Mst2))*(4*pow2(Msq)*pow2(Mst2) + 2*pow4(Msq) + 5*pow4(Mst2)) + 8*pow2(
        Msq)*(pow2(Mst1) + pow2(Mst2))*pow5(Mst2) - 3*Dmglst2*(-8*pow2(Msq)*(
        pow2(Mst1) + pow2(Mst2))*pow4(Mst2) + 2*pow4(Msq)*(-7*pow2(Mst1)*pow2(
        Mst2) + 3*pow4(Mst1) + pow4(Mst2)) - 5*(pow2(Mst1) + pow2(Mst2))*pow6(
        Mst2)) + 3*pow2(Mst1)*pow7(Mst2)) + Dmglst2*(288*pow2(Msq)*(pow2(Mst1)
        - pow2(Mst2))*pow4(Mst2) - 3*pow4(Mst1)*pow4(Mst2) + 18*pow4(Msq)*(8*
        pow2(Mst1)*pow2(Mst2) + 7*pow4(Mst1) + 60*pow4(Mst2)) + 183*pow2(Mst1)*
        pow6(Mst2) - 181*pow8(Mst2))) + 7547850*Dmglst2*pow2(Mst1)*pow8(Mst2) -
        180*log(pow2(Mst2)/pow2(Mst1))*(6*Mst2*pow2(Dmglst2)*(-5805*pow4(Mst1)*
        pow4(Mst2) + pow4(Msq)*(12*(-2201 + 252*S2)*pow2(Mst1)*pow2(Mst2) + 34*
        (-694 + 1449*S2)*pow4(Mst1) + 2304*pow4(Mst2)) - 720*pow2(Msq)*(7*pow2(
        Mst2)*pow4(Mst1) - 5*pow2(Mst1)*pow4(Mst2)) + 4995*pow2(Mst1)*pow6(
        Mst2)) + 9*Mst2*(-320*pow2(Msq)*pow2(Mst1)*(pow2(Mst1) - 5*pow2(Mst2))*
        pow4(Mst2) + 30*pow4(Mst2)*pow6(Mst1) + 2*pow4(Msq)*(6*(326 + 219*S2)*
        pow2(Mst2)*pow4(Mst1) + 3*(-2326 + 189*S2)*pow2(Mst1)*pow4(Mst2) + (136
        + 1317*S2)*pow6(Mst1) - 384*pow6(Mst2)) - 255*pow4(Mst1)*pow6(Mst2)) -
        Dmglst2*(17280*pow2(Msq)*pow2(Mst1)*(pow2(Mst1) - 2*pow2(Mst2))*pow4(
        Mst2) + 2*pow4(Msq)*(18*(-818 + 3609*S2)*pow2(Mst2)*pow4(Mst1) + 135*(
        598 + 63*S2)*pow2(Mst1)*pow4(Mst2) + (24322 + 140139*S2)*pow6(Mst1) +
        3456*pow6(Mst2)) - 45*(6*pow4(Mst2)*pow6(Mst1) - 303*pow4(Mst1)*pow6(
        Mst2) + 425*pow2(Mst1)*pow8(Mst2))))))/(87480.*pow2(Mst1)*pow2(Mst2)*
        pow4(Msq)) + z2*(pow4(Mt)*(402.837037037037 - (80*log(pow2(Msq)/pow2(
        Mst1)))/3. - (556*log(pow2(Mst2)/pow2(Mst1)))/3. + pow2(Dmglst2)*(80/(
        3.*pow2(Msq)) + 416/(9.*pow2(Mst1))) - (80*pow2(Msq))/(3.*pow2(Mst1)) +
        (8*pow2(Mst2)*(50/pow2(Msq) + 49/pow2(Mst1) + (90*pow2(Dmglst2))/pow4(
        Msq)))/9. + (40*(16*Dmglst2 + 7*Mst2)*pow3(Mst2))/(9.*pow4(Msq)) + ((
        2722862*pow2(Dmglst2))/2835. - (80*pow2(Msq))/3. + (310039*pow2(Mst1))/
        405. - (160*log(pow2(Msq)/pow2(Mst1))*pow2(Mst1))/3. - (8*log(pow2(
        Mst2)/pow2(Mst1))*(4*pow2(Dmglst2) + 343*pow2(Mst1)))/9. + (80*pow4(
        Mst1))/(3.*pow2(Msq)))/pow2(Mst2) - (2*pow2(Mst1)*(-47240777*pow2(
        Dmglst2) - 4274683*pow2(Mst1) + 340200*log(pow2(Msq)/pow2(Mst1))*(2*
        pow2(Dmglst2) + 3*pow2(Mst1)) + 7560*log(pow2(Mst2)/pow2(Mst1))*(887*
        pow2(Dmglst2) + 231*pow2(Mst1))))/(8505.*pow4(Mst2)) - Dmglst2*((4*(
        114841 - 31920*log(pow2(Mst2)/pow2(Mst1))))/(945.*Mst2) - (32*Mst2*(15/
        pow2(Msq) + 23/pow2(Mst1)))/9. + (4*pow2(Mst1)*(2655197 - 75600*log(
        pow2(Msq)/pow2(Mst1)) - 778680*log(pow2(Mst2)/pow2(Mst1)) + (37800*
        pow2(Mst1))/pow2(Msq)))/(2835.*pow3(Mst2)) + ((10677.005369390554 -
        960*log(pow2(Msq)/pow2(Mst1)) - (8128*log(pow2(Mst2)/pow2(Mst1)))/3.)*
        pow4(Mst1))/pow5(Mst2))) + Mt*pow3(s2t)*((pow2(Mst1)*(2934146*pow2(
        Dmglst2) + 1080*log(pow2(Mst2)/pow2(Mst1))*(350*pow2(Dmglst2) - 153*
        pow2(Mst1)) + 104839*pow2(Mst1)))/(4860.*Mst2) + Mst2*((-126374*pow2(
        Dmglst2))/405. + (4*log(pow2(Mst2)/pow2(Mst1))*(195*pow2(Dmglst2) - 49*
        pow2(Mst1)))/3. + (7913*pow2(Mst1))/162. - (80*pow2(Dmglst2)*pow2(Mst1)
        )/(3.*pow2(Msq)) + (40*log(pow2(Msq)/pow2(Mst1))*(-12*pow2(Dmglst2) +
        pow2(Mst1)))/3.) - (pow3(Mst2)*(3469 + 1440*log(pow2(Msq)/pow2(Mst1)) -
        1384*log(pow2(Mst2)/pow2(Mst1)) + (320*pow2(Mst1))/pow2(Msq) - (240*
        pow2(Dmglst2)*(36*pow2(Msq) - 5*pow2(Mst1)))/pow4(Msq)))/36. - Dmglst2*
        ((288.56913580246913 - (40*log(pow2(Msq)/pow2(Mst1)))/3. - 44*log(pow2(
        Mst2)/pow2(Mst1)))*pow2(Mst1) + (515.4574074074075 + (440*log(pow2(Msq)
        /pow2(Mst1)))/3. - (662*log(pow2(Mst2)/pow2(Mst1)))/3. + (80*pow2(Mst1)
        )/(3.*pow2(Msq)))*pow2(Mst2) + ((562.0358710562414 - (1438*log(pow2(
        Mst2)/pow2(Mst1)))/9.)*pow4(Mst1))/pow2(Mst2) - (10*(96*pow2(Msq) - 15*
        pow2(Mst1) + 41*pow2(Mst2))*pow4(Mst2))/(9.*pow4(Msq))) + (10*(138*
        pow2(Dmglst2) + 16*pow2(Msq) - 3*pow2(Mst1))*pow5(Mst2))/(9.*pow4(Msq))
        ) + (pow2(Mt)*pow2(s2t)*(7560*log(pow2(Mst2)/pow2(Mst1))*pow2(Mst1)*(
        68*pow2(Mst1)*pow3(Mst2) + pow2(Dmglst2)*(2316*Mst2*pow2(Mst1) + 2102*
        pow3(Mst2)) + 195*Mst2*pow4(Mst1) + 4*Dmglst2*(124*pow2(Mst1)*pow2(
        Mst2) + 5*pow4(Mst1) + 297*pow4(Mst2)) + 305*pow5(Mst2)) - 36*pow2(
        Dmglst2)*(63809*pow2(Mst1)*pow3(Mst2) + 801215*Mst2*pow4(Mst1) + 21840*
        pow5(Mst2)) + 7*Mst2*(64800*pow2(Msq)*pow2(pow2(Mst1) - pow2(Mst2)) +
        1774344*pow2(Mst2)*pow4(Mst1) - 161307*pow2(Mst1)*pow4(Mst2) + 3321329*
        pow6(Mst1) - 105840*pow6(Mst2)) - 168*Dmglst2*(99000*pow2(Mst2)*pow4(
        Mst1) + 20994*pow2(Mst1)*pow4(Mst2) + 392665*pow6(Mst1) + 8280*pow6(
        Mst2))))/(34020.*pow2(Mst1)*pow3(Mst2)) + (s2t*((-7*pow2(Mst2)*pow3(
        s2t)*(4050*pow2(Mst1)*pow3(Mst2)*(pow2(Mst2)*(-30*pow2(Dmglst2)*(7*
        Mst2*pow2(Mst1) - 6*pow3(Mst2)) - 14*pow2(Mst1)*pow3(Mst2) + 5*Mst2*
        pow4(Mst1) + 4*Dmglst2*(-21*pow2(Mst1)*pow2(Mst2) + 5*pow4(Mst1))) + 4*
        pow2(Msq)*(-10*pow2(Mst1)*pow3(Mst2) + pow2(Dmglst2)*(-60*Mst2*pow2(
        Mst1) + 57*pow3(Mst2)) + 3*Mst2*pow4(Mst1) + Dmglst2*(-40*pow2(Mst1)*
        pow2(Mst2) + 6*pow4(Mst1) + 34*pow4(Mst2)) + 7*pow5(Mst2))) + 97200*(
        pow2(Mst1) - pow2(Mst2))*(2*log(pow2(Mst2)/pow2(Mst1))*pow2(Mst1)*pow2(
        Mst2) + pow4(Mst1) - pow4(Mst2))*pow6(Msq) + Mst2*pow4(Msq)*(-48600*
        Mst2*log(pow2(Msq)/pow2(Mst1))*pow2(Mst1)*(-12*Dmglst2*Mst2 - 6*pow2(
        Dmglst2) + pow2(Mst1) - 5*pow2(Mst2))*(pow2(Mst1) - pow2(Mst2)) +
        491742*pow2(Dmglst2)*pow2(Mst1)*pow3(Mst2) + 2511120*Mst2*pow2(Dmglst2)
        *pow4(Mst1) - 1681344*Dmglst2*pow2(Mst2)*pow4(Mst1) + 4231782*pow3(
        Mst2)*pow4(Mst1) - 1645380*Dmglst2*pow2(Mst1)*pow4(Mst2) - 168480*pow2(
        Dmglst2)*pow5(Mst2) + 285660*pow2(Mst1)*pow5(Mst2) + 1620*log(pow2(
        Mst2)/pow2(Mst1))*pow2(Mst1)*(-670*pow2(Mst1)*pow3(Mst2) + pow2(
        Dmglst2)*(-398*Mst2*pow2(Mst1) + 488*pow3(Mst2)) + 75*Mst2*pow4(Mst1) +
        4*Dmglst2*(-223*pow2(Mst1)*pow2(Mst2) + 132*pow4(Mst1) + 290*pow4(Mst2)
        ) + 454*pow5(Mst2)) - 8800364*Dmglst2*pow6(Mst1) + 2256555*Mst2*pow6(
        Mst1) - 298080*Dmglst2*pow6(Mst2) - 158760*pow7(Mst2))))/pow2(Mst1) +
        16*pow3(Mt)*(7*Dmglst2*(pow4(Msq)*(11830023*pow2(Mst1)*pow2(Mst2) +
        36911068*pow4(Mst1) - 3240*log(pow2(Mst2)/pow2(Mst1))*(505*pow2(Mst1)*
        pow2(Mst2) + 1344*pow4(Mst1) - 359*pow4(Mst2)) + 821745*pow4(Mst2) -
        388800*log(pow2(Msq)/pow2(Mst1))*(3*pow4(Mst1) + pow4(Mst2))) + 907200*
        pow2(Msq)*pow6(Mst2) + 712800*pow8(Mst2)) + Mst2*(24*pow2(Dmglst2)*((-
        6085792*pow2(Mst1) + 38181*pow2(Mst2) + 3780*log(pow2(Mst2)/pow2(Mst1))
        *(149*pow2(Mst1) + 5*pow2(Mst2)))*pow4(Msq) + 226800*pow2(Msq)*pow4(
        Mst2) + 340200*pow6(Mst2)) + 21*(-(pow4(Msq)*(1411971*pow2(Mst1)*pow2(
        Mst2) + 2299036*pow4(Mst1) - 129600*log(pow2(Msq)/pow2(Mst1))*(pow4(
        Mst1) - pow4(Mst2)) + 780237*pow4(Mst2) - 1080*log(pow2(Mst2)/pow2(
        Mst1))*(495*pow2(Mst1)*pow2(Mst2) + 488*pow4(Mst1) + 439*pow4(Mst2))))
        + 129600*pow2(Msq)*pow6(Mst2) + 64800*pow8(Mst2))))))/(408240.*pow4(
        Msq)*pow4(Mst2))) + (s2t*pow3(Mt)*(-332438184*pow2(Dmglst2)*pow2(Mst1)*
        pow3(Mst2)*pow4(Msq) + 1128960*OepS2*pow2(Dmglst2)*pow2(Mst1)*pow3(
        Mst2)*pow4(Msq) + 68630976*S2*pow2(Dmglst2)*pow2(Mst1)*pow3(Mst2)*pow4(
        Msq) + 9525600*pow2(Dmglst2)*pow2(Msq)*pow3(Mst2)*pow4(Mst1) -
        400437672*Mst2*pow2(Dmglst2)*pow4(Msq)*pow4(Mst1) + 28412160*Mst2*
        OepS2*pow2(Dmglst2)*pow4(Msq)*pow4(Mst1) + 139706208*Mst2*S2*pow2(
        Dmglst2)*pow4(Msq)*pow4(Mst1) + 128297169*Dmglst2*pow2(Mst2)*pow4(Msq)*
        pow4(Mst1) - 13940640*Dmglst2*OepS2*pow2(Mst2)*pow4(Msq)*pow4(Mst1) +
        106761564*Dmglst2*S2*pow2(Mst2)*pow4(Msq)*pow4(Mst1) - 66716559*pow3(
        Mst2)*pow4(Msq)*pow4(Mst1) + 3800160*OepS2*pow3(Mst2)*pow4(Msq)*pow4(
        Mst1) - 89139204*S2*pow3(Mst2)*pow4(Msq)*pow4(Mst1) + 1873935*Dmglst2*
        pow2(Mst1)*pow4(Msq)*pow4(Mst2) - 1058400*Dmglst2*OepS2*pow2(Mst1)*
        pow4(Msq)*pow4(Mst2) - 41334300*Dmglst2*S2*pow2(Mst1)*pow4(Msq)*pow4(
        Mst2) - 680400*pow2(Mst1)*pow2(log(pow2(Msq)/pow2(Mst1)))*((pow2(Mst1)
        - pow2(Mst2))*(6*Mst2*pow2(Dmglst2) + Mst2*(4*pow2(Msq) + pow2(Mst1) +
        3*pow2(Mst2)) + Dmglst2*(4*pow2(Msq) + pow2(Mst1) + 7*pow2(Mst2))) - 6*
        (Dmglst2 + Mst2)*log(pow2(Mst2)/pow2(Mst1))*pow4(Msq))*pow4(Mst2) -
        9903600*Dmglst2*pow2(Msq)*pow4(Mst1)*pow4(Mst2) - 60782400*pow2(
        Dmglst2)*pow2(Msq)*pow2(Mst1)*pow5(Mst2) + 13063680*pow2(Dmglst2)*pow4(
        Msq)*pow5(Mst2) - 8935353*pow2(Mst1)*pow4(Msq)*pow5(Mst2) + 635040*
        OepS2*pow2(Mst1)*pow4(Msq)*pow5(Mst2) + 14512932*S2*pow2(Mst1)*pow4(
        Msq)*pow5(Mst2) - 5329800*pow2(Dmglst2)*pow4(Mst1)*pow5(Mst2) -
        7182000*pow2(Msq)*pow4(Mst1)*pow5(Mst2) + 1451520*pow2(Mst1)*pow3(log(
        pow2(Mst2)/pow2(Mst1)))*pow4(Msq)*(6*Mst2*pow2(Dmglst2)*pow2(Mst1) + 9*
        pow2(Mst1)*pow3(Mst2) + 2*Mst2*pow4(Mst1) + 3*Dmglst2*(-(pow2(Mst1)*
        pow2(Mst2)) + 3*pow4(Mst1) + 4*pow4(Mst2)) + 12*pow5(Mst2)) + 20412000*
        Dmglst2*pow2(Msq)*pow2(Mst2)*pow6(Mst1) - 9525600*pow2(Msq)*pow3(Mst2)*
        pow6(Mst1) + 677836516*Dmglst2*pow4(Msq)*pow6(Mst1) - 160234788*Mst2*
        pow4(Msq)*pow6(Mst1) - 57635200*Dmglst2*OepS2*pow4(Msq)*pow6(Mst1) +
        9206400*Mst2*OepS2*pow4(Msq)*pow6(Mst1) + 1331267616*Dmglst2*S2*pow4(
        Msq)*pow6(Mst1) - 295692768*Mst2*S2*pow4(Msq)*pow6(Mst1) - 2929500*
        Dmglst2*pow4(Mst2)*pow6(Mst1) - 888300*pow5(Mst2)*pow6(Mst1) -
        56700000*Dmglst2*pow2(Msq)*pow2(Mst1)*pow6(Mst2) + 4354560*Dmglst2*
        pow4(Msq)*pow6(Mst2) - 4649400*Dmglst2*pow4(Mst1)*pow6(Mst2) - 2520*
        log(pow2(Mst2)/pow2(Mst1))*(3*Mst2*pow2(Dmglst2)*(-360*pow2(Msq)*pow2(
        Mst1)*pow2(Mst2)*(2*pow2(Mst1) + pow2(Mst2)) + 90*pow2(Mst1)*(-11*pow2(
        Mst1) + pow2(Mst2))*pow4(Mst2) + pow4(Msq)*(56*(-119 + 54*S2)*pow2(
        Mst1)*pow2(Mst2) + (-42953 + 76104*S2)*pow4(Mst1) + 2304*pow4(Mst2))) +
        9*Mst2*(15*pow2(Mst1)*(6*pow2(Mst1)*pow2(Mst2) + 3*pow4(Mst1) - pow4(
        Mst2))*pow4(Mst2) + pow4(Msq)*(9*(102 + 377*S2)*pow2(Mst2)*pow4(Mst1) +
        3*(116 + 189*S2)*pow2(Mst1)*pow4(Mst2) + 4*(302 + 2055*S2)*pow6(Mst1) -
        384*pow6(Mst2)) - 60*pow2(Msq)*(-(pow4(Mst1)*pow4(Mst2)) + 4*pow2(Mst2)
        *pow6(Mst1) - 7*pow2(Mst1)*pow6(Mst2))) + Dmglst2*(45*pow2(Mst1)*pow4(
        Mst2)*(6*pow2(Mst1)*pow2(Mst2) - 15*pow4(Mst1) + 25*pow4(Mst2)) + pow4(
        Msq)*(3*(18428 - 37341*S2)*pow2(Mst2)*pow4(Mst1) + 3*(22768 - 2835*S2)*
        pow2(Mst1)*pow4(Mst2) + (14036 - 463140*S2)*pow6(Mst1) - 3456*pow6(
        Mst2)) + 180*pow2(Msq)*(-21*pow4(Mst1)*pow4(Mst2) + 12*pow2(Mst2)*pow6(
        Mst1) + 37*pow2(Mst1)*pow6(Mst2)))) - 76318200*pow2(Dmglst2)*pow2(Mst1)
        *pow7(Mst2) - 17992800*pow2(Msq)*pow2(Mst1)*pow7(Mst2) - 4354560*pow4(
        Msq)*pow7(Mst2) - 189000*pow4(Mst1)*pow7(Mst2) - 45360*pow2(log(pow2(
        Mst2)/pow2(Mst1)))*(-60*pow2(Msq)*pow2(Mst1)*(pow2(Mst1) + 2*pow2(Mst2)
        )*pow5(Mst2) - 15*pow2(Mst1)*(2*pow2(Mst1)*pow2(Mst2) + pow4(Mst1) + 3*
        pow4(Mst2))*pow5(Mst2) + pow2(Dmglst2)*(-360*pow2(Msq)*pow2(Mst1)*pow5(
        Mst2) - 90*pow2(Mst1)*(pow2(Mst1) + 5*pow2(Mst2))*pow5(Mst2) + pow4(
        Msq)*(-952*pow2(Mst1)*pow3(Mst2) + 2129*Mst2*pow4(Mst1) + 288*pow5(
        Mst2))) + Dmglst2*(-60*pow2(Msq)*pow2(Mst1)*(pow2(Mst1) + 6*pow2(Mst2))
        *pow4(Mst2) - 15*pow2(Mst1)*pow4(Mst2)*(6*pow2(Mst1)*pow2(Mst2) + pow4(
        Mst1) + 15*pow4(Mst2)) + 2*pow4(Msq)*(-917*pow2(Mst2)*pow4(Mst1) - 349*
        pow2(Mst1)*pow4(Mst2) + 352*pow6(Mst1) + 144*pow6(Mst2))) + 2*pow4(Msq)
        *(432*pow3(Mst2)*pow4(Mst1) + 555*pow2(Mst1)*pow5(Mst2) + 70*Mst2*pow6(
        Mst1) + 48*pow7(Mst2))) - 42751800*Dmglst2*pow2(Mst1)*pow8(Mst2) -
        11680200*pow2(Mst1)*pow9(Mst2) - 226800*log(pow2(Msq)/pow2(Mst1))*pow2(
        Mst1)*(36*pow2(log(pow2(Mst2)/pow2(Mst1)))*pow4(Msq)*(-(Mst2*pow4(Mst1)
        ) + Dmglst2*(3*pow4(Mst1) + pow4(Mst2)) + pow5(Mst2)) + 6*pow2(Dmglst2)
        *((-51*Mst2*pow2(Mst1) + 5*pow3(Mst2))*pow4(Msq) - pow2(Mst1)*pow5(
        Mst2) + pow7(Mst2)) - Mst2*(14*pow2(Msq)*(-pow2(Mst1) + pow2(Mst2))*
        pow4(Mst2) + 4*pow4(Mst1)*pow4(Mst2) + 9*pow4(Msq)*(10*pow2(Mst1)*pow2(
        Mst2) + 4*pow4(Mst1) + 13*pow4(Mst2)) + 2*pow2(Mst1)*pow6(Mst2) - 6*
        pow8(Mst2)) + Dmglst2*(pow4(Msq)*(234*pow2(Mst1)*pow2(Mst2) + 504*pow4(
        Mst1) - 161*pow4(Mst2)) + 14*pow2(Msq)*(pow2(Mst1) - pow2(Mst2))*pow4(
        Mst2) - 4*pow4(Mst1)*pow4(Mst2) - 6*pow2(Mst1)*pow6(Mst2) + 10*pow8(
        Mst2)) + 6*log(pow2(Mst2)/pow2(Mst1))*(6*pow2(Dmglst2)*((2*Mst2*pow2(
        Mst1) + pow3(Mst2))*pow4(Msq) + 2*pow2(Msq)*pow5(Mst2) + 3*pow7(Mst2))
        + Dmglst2*(pow4(Msq)*(-12*pow2(Mst1)*pow2(Mst2) - 123*pow4(Mst1) + 10*
        pow4(Mst2)) + 14*pow2(Msq)*pow6(Mst2) + 11*pow8(Mst2)) + 3*(pow4(Msq)*(
        4*pow2(Mst1)*pow3(Mst2) + 11*Mst2*pow4(Mst1) - 2*pow5(Mst2)) + 2*pow2(
        Msq)*pow7(Mst2) + pow9(Mst2))))))/(153090.*pow2(Mst1)*pow4(Msq)*pow4(
        Mst2)) - (pow4(Mt)*(2326826250*pow11(Mst2)*pow2(Mst1) - 1116118935540*
        pow2(Dmglst2)*pow3(Mst2)*pow4(Msq)*pow4(Mst1) + 6569136000*OepS2*pow2(
        Dmglst2)*pow3(Mst2)*pow4(Msq)*pow4(Mst1) + 261770632200*S2*pow2(
        Dmglst2)*pow3(Mst2)*pow4(Msq)*pow4(Mst1) - 262937253600*Dmglst2*pow4(
        Msq)*pow4(Mst1)*pow4(Mst2) - 2074464000*Dmglst2*OepS2*pow4(Msq)*pow4(
        Mst1)*pow4(Mst2) - 68284263600*Dmglst2*S2*pow4(Msq)*pow4(Mst1)*pow4(
        Mst2) + 44591715000*pow2(Dmglst2)*pow2(Mst1)*pow4(Msq)*pow5(Mst2) +
        14424933600*pow2(Dmglst2)*pow2(Msq)*pow4(Mst1)*pow5(Mst2) -
        151045520850*pow4(Msq)*pow4(Mst1)*pow5(Mst2) - 1555848000*OepS2*pow4(
        Msq)*pow4(Mst1)*pow5(Mst2) + 49059221400*S2*pow4(Msq)*pow4(Mst1)*pow5(
        Mst2) + 1666980000*pow3(log(pow2(Msq)/pow2(Mst1)))*pow4(Msq)*pow4(Mst1)
        *pow5(Mst2) + 55566000*pow3(log(pow2(Mst2)/pow2(Mst1)))*pow4(Msq)*pow4(
        Mst1)*(6*pow2(Dmglst2)*(111*Mst2*pow2(Mst1) - 32*pow3(Mst2)) - 290*
        pow2(Mst1)*pow3(Mst2) - 1210*Mst2*pow4(Mst1) + 4*Dmglst2*(-239*pow2(
        Mst1)*pow2(Mst2) + 346*pow4(Mst1) - 96*pow4(Mst2)) + 111*pow5(Mst2)) -
        5000940000*pow3(Mst2)*pow4(Mst1)*pow6(Msq) - 5000940000*pow2(Mst1)*
        pow5(Mst2)*pow6(Msq) + 11668860000*pow2(Dmglst2)*pow2(Msq)*pow3(Mst2)*
        pow6(Mst1) - 2025139461828*Mst2*pow2(Dmglst2)*pow4(Msq)*pow6(Mst1) +
        93334416000*Mst2*OepS2*pow2(Dmglst2)*pow4(Msq)*pow6(Mst1) +
        39251028600*Mst2*S2*pow2(Dmglst2)*pow4(Msq)*pow6(Mst1) - 81313523928*
        Dmglst2*pow2(Mst2)*pow4(Msq)*pow6(Mst1) - 15904224000*Dmglst2*OepS2*
        pow2(Mst2)*pow4(Msq)*pow6(Mst1) + 62988030000*Dmglst2*S2*pow2(Mst2)*
        pow4(Msq)*pow6(Mst1) - 213216589530*pow3(Mst2)*pow4(Msq)*pow6(Mst1) -
        5161464000*OepS2*pow3(Mst2)*pow4(Msq)*pow6(Mst1) + 208189132200*S2*
        pow3(Mst2)*pow4(Msq)*pow6(Mst1) + 3333960000*Dmglst2*pow2(Msq)*pow4(
        Mst2)*pow6(Mst1) - 6598462500*pow2(Dmglst2)*pow5(Mst2)*pow6(Mst1) -
        16657205040*pow2(Msq)*pow5(Mst2)*pow6(Mst1) - 11907000*pow2(Mst1)*pow2(
        log(pow2(Msq)/pow2(Mst1)))*pow3(Mst2)*(210*log(pow2(Mst2)/pow2(Mst1))*
        pow2(Mst1)*pow2(Mst2)*pow4(Msq) + 210*pow2(Dmglst2)*(4*pow2(Msq)*pow2(
        Mst1)*pow2(Mst2) + pow2(Mst1)*pow2(Mst2)*(pow2(Mst1) + 11*pow2(Mst2)) +
        3*(pow2(Mst1) + pow2(Mst2))*pow4(Msq)) + 140*Dmglst2*(12*pow2(Msq)*
        pow2(Mst1)*pow3(Mst2) + 9*Mst2*(pow2(Mst1) + pow2(Mst2))*pow4(Msq) + 3*
        pow3(Mst2)*pow4(Mst1) + 13*pow2(Mst1)*pow5(Mst2)) + pow2(Mst2)*(70*(11*
        pow2(Mst1) + 9*pow2(Mst2))*pow4(Msq) + 210*pow2(Mst2)*pow4(Mst1) + 28*
        pow2(Msq)*(37*pow2(Mst1)*pow2(Mst2) + 7*pow4(Mst1)) + 617*pow2(Mst1)*
        pow4(Mst2) + 57*pow6(Mst1))) + 40507614000*Dmglst2*pow2(Mst1)*pow4(Msq)
        *pow6(Mst2) + 8334900000*Dmglst2*pow2(Msq)*pow4(Mst1)*pow6(Mst2) -
        11529945000*Dmglst2*pow6(Mst1)*pow6(Mst2) + 13335840000*pow2(Dmglst2)*
        pow2(Msq)*pow2(Mst1)*pow7(Mst2) - 10668672000*pow2(Dmglst2)*pow4(Msq)*
        pow7(Mst2) - 23087673000*pow2(Mst1)*pow4(Msq)*pow7(Mst2) + 12462064650*
        pow2(Dmglst2)*pow4(Mst1)*pow7(Mst2) + 17423274960*pow2(Msq)*pow4(Mst1)*
        pow7(Mst2) - 347287500*pow6(Mst1)*pow7(Mst2) + 145027260000*Dmglst2*
        pow2(Msq)*pow2(Mst2)*pow8(Mst1) - 32506110000*pow2(Msq)*pow3(Mst2)*
        pow8(Mst1) + 1455111565780*Dmglst2*pow4(Msq)*pow8(Mst1) - 527147595309*
        Mst2*pow4(Msq)*pow8(Mst1) - 61268032000*Dmglst2*OepS2*pow4(Msq)*pow8(
        Mst1) - 9614976000*Mst2*OepS2*pow4(Msq)*pow8(Mst1) + 1059057266400*
        Dmglst2*S2*pow4(Msq)*pow8(Mst1) + 384650078400*Mst2*S2*pow4(Msq)*pow8(
        Mst1) - 6370213005*pow5(Mst2)*pow8(Mst1) + 31116960000*Dmglst2*pow2(
        Msq)*pow2(Mst1)*pow8(Mst2) + 29334217500*Dmglst2*pow4(Mst1)*pow8(Mst2)
        - 2520*log(pow2(Mst2)/pow2(Mst1))*(147*Mst2*pow2(Dmglst2)*(2*pow4(Msq)*
        (2*(-317567 + 89775*S2)*pow2(Mst2)*pow4(Mst1) + 81225*pow2(Mst1)*pow4(
        Mst2) + (-1582393 + 2551050*S2)*pow6(Mst1) - 7200*pow6(Mst2)) + 900*
        pow2(Msq)*(51*pow4(Mst1)*pow4(Mst2) + 20*pow2(Mst2)*pow6(Mst1) + 150*
        pow2(Mst1)*pow6(Mst2)) - 225*(115*pow4(Mst2)*pow6(Mst1) - 526*pow4(
        Mst1)*pow6(Mst2) - 525*pow2(Mst1)*pow8(Mst2))) - 3*Mst2*(-15*pow2(Mst1)
        *(146492*pow2(Mst1)*pow2(Mst2) + 47775*pow4(Mst1) + 25725*pow4(Mst2))*
        pow6(Mst2) + pow4(Msq)*(2450*(722 + 1701*S2)*pow4(Mst1)*pow4(Mst2) +
        5390*(-3418 + 2565*S2)*pow2(Mst2)*pow6(Mst1) + 5666850*pow2(Mst1)*pow6(
        Mst2) + (-47778491 + 25754400*S2)*pow8(Mst1) - 705600*pow8(Mst2)) +
        2940*pow2(Msq)*(75*pow4(Mst2)*pow6(Mst1) - 866*pow4(Mst1)*pow6(Mst2) +
        1800*pow2(Mst2)*pow8(Mst1) - 375*pow2(Mst1)*pow8(Mst2))) - 2*Dmglst2*(
        55125*pow2(Mst1)*pow4(Mst2)*(-3*pow2(Mst2)*pow4(Mst1) - 136*pow2(Mst1)*
        pow4(Mst2) + 18*pow6(Mst1) - 63*pow6(Mst2)) + pow4(Msq)*(3675*(-11369 +
        2268*S2)*pow4(Mst1)*pow4(Mst2) + 882*(8581 + 72450*S2)*pow2(Mst2)*pow6(
        Mst1) + 3638250*pow2(Mst1)*pow6(Mst2) + 5*(57911521 + 49233240*S2)*
        pow8(Mst1) - 2116800*pow8(Mst2)) - 1323000*pow2(Msq)*(-2*pow4(Mst2)*
        pow6(Mst1) + 3*pow4(Mst1)*pow6(Mst2) + 15*pow2(Mst2)*pow8(Mst1) + 5*
        pow2(Mst1)*pow8(Mst2)))) - 793800*pow2(log(pow2(Mst2)/pow2(Mst1)))*(14*
        pow2(Dmglst2)*(-450*pow2(Msq)*pow4(Mst1)*pow5(Mst2) - 225*(pow2(Mst1) +
        5*pow2(Mst2))*pow4(Mst1)*pow5(Mst2) + pow4(Msq)*(9941*pow3(Mst2)*pow4(
        Mst1) + 1560*pow2(Mst1)*pow5(Mst2) + 49945*Mst2*pow6(Mst1) - 1440*pow7(
        Mst2))) - Mst2*(15*pow4(Mst1)*pow4(Mst2)*(210*pow2(Mst1)*pow2(Mst2) +
        105*pow4(Mst1) + 223*pow4(Mst2)) - 420*pow2(Msq)*(-10*pow4(Mst2)*pow6(
        Mst1) - 18*pow4(Mst1)*pow6(Mst2) + 15*pow2(Mst2)*pow8(Mst1)) + pow4(
        Msq)*(-56210*pow4(Mst1)*pow4(Mst2) + 13720*pow2(Mst2)*pow6(Mst1) -
        21840*pow2(Mst1)*pow6(Mst2) + 273641*pow8(Mst1) + 3360*pow8(Mst2))) +
        4*Dmglst2*(-3150*pow2(Msq)*pow2(Mst2)*pow4(Mst1)*(pow4(Mst1) + pow4(
        Mst2)) + pow4(Msq)*(-36995*pow4(Mst1)*pow4(Mst2) - 65163*pow2(Mst2)*
        pow6(Mst1) + 10920*pow2(Mst1)*pow6(Mst2) + 158930*pow8(Mst1) - 3360*
        pow8(Mst2)) - 525*(3*pow6(Mst1)*pow6(Mst2) + 5*pow4(Mst1)*pow8(Mst2))))
        + 2813028750*pow2(Dmglst2)*pow2(Mst1)*pow9(Mst2) + 11946690000*pow2(
        Msq)*pow2(Mst1)*pow9(Mst2) + 2667168000*pow4(Msq)*pow9(Mst2) +
        19213299495*pow4(Mst1)*pow9(Mst2) + 56700*log(pow2(Msq)/pow2(Mst1))*
        pow2(Mst1)*(51450*pow11(Mst2) + 1058400*pow3(Mst2)*pow4(Msq)*pow4(Mst1)
        + 127400*pow2(Mst1)*pow4(Msq)*pow5(Mst2) + 16072*pow2(Msq)*pow4(Mst1)*
        pow5(Mst2) - 88200*pow2(Mst1)*pow2(log(pow2(Mst2)/pow2(Mst1)))*pow4(
        Msq)*(-6*Mst2*pow2(Dmglst2)*pow2(Mst1) + 4*Dmglst2*pow2(Mst1)*(9*pow2(
        Mst1) + pow2(Mst2)) - 2*pow2(Mst1)*pow3(Mst2) - 9*Mst2*pow4(Mst1) +
        pow5(Mst2)) - 176400*pow2(Mst1)*pow3(Mst2)*pow6(Msq) - 176400*pow5(
        Mst2)*pow6(Msq) + 926100*Mst2*pow4(Msq)*pow6(Mst1) + 49209*pow5(Mst2)*
        pow6(Mst1) + 45472*pow2(Msq)*pow2(Mst1)*pow7(Mst2) - 220500*pow4(Msq)*
        pow7(Mst2) + 58800*pow4(Mst1)*pow7(Mst2) - 4900*Dmglst2*(2*pow4(Msq)*(
        486*pow2(Mst2)*pow4(Mst1) + 23*pow2(Mst1)*pow4(Mst2) + 1242*pow6(Mst1)
        - 9*pow6(Mst2)) - (137*pow2(Mst1)*pow2(Mst2) + 24*pow4(Mst1) + 63*pow4(
        Mst2))*pow6(Mst2) - 24*pow2(Msq)*(3*pow2(Mst1)*pow6(Mst2) + 5*pow8(
        Mst2))) + 147000*pow2(Msq)*pow9(Mst2) + 220709*pow2(Mst1)*pow9(Mst2) -
        420*log(pow2(Mst2)/pow2(Mst1))*(630*pow4(Msq)*(5*pow3(Mst2)*pow4(Mst1)
        + 16*Mst2*pow6(Mst1) - pow7(Mst2)) - 266*pow2(Msq)*pow2(Mst1)*pow7(
        Mst2) + 42*pow2(Dmglst2)*(pow4(Msq)*(-26*pow2(Mst1)*pow3(Mst2) + 325*
        Mst2*pow4(Mst1) - 15*pow5(Mst2)) - 5*pow2(Msq)*pow2(Mst1)*pow5(Mst2) -
        15*pow2(Mst1)*pow7(Mst2)) - 140*Dmglst2*(3*pow2(Msq)*pow2(Mst1)*pow6(
        Mst2) + pow4(Msq)*(57*pow2(Mst2)*pow4(Mst1) + pow2(Mst1)*pow4(Mst2) +
        342*pow6(Mst1) + 9*pow6(Mst2)) + 4*pow2(Mst1)*pow8(Mst2)) - 197*pow2(
        Mst1)*pow9(Mst2)) + 294*pow2(Dmglst2)*(2*pow4(Msq)*(302*pow2(Mst1)*
        pow3(Mst2) + 17850*Mst2*pow4(Mst1) + 975*pow5(Mst2)) + 200*pow2(Msq)*(
        13*pow2(Mst1)*pow5(Mst2) + 15*pow7(Mst2)) + 25*(8*pow4(Mst1)*pow5(Mst2)
        + 163*pow2(Mst1)*pow7(Mst2) + 105*pow9(Mst2)))) + 8126527500*Dmglst2*
        pow2(Mst1)*power10(Mst2)))/(3.750705e8*pow4(Msq)*pow4(Mst1)*pow5(Mst2))
        ))/pow4(Mt)*12.;

   return result;
}

double himalaya::H6::calc_coef_at_as2_no_sm_logs_log1() const {

   const double result =
     (((115200*pow11(Mst2)*pow2(Mt)*pow2(s2t)*pow4(Msq) - 2451200*
        s2t*pow2(Dmglst2)*pow3(Mt)*pow4(Msq)*pow4(Mst1)*pow4(Mst2) - 547200*
        s2t*z2*pow2(Dmglst2)*pow3(Mt)*pow4(Msq)*pow4(Mst1)*pow4(Mst2) +
        13177568*pow2(Dmglst2)*pow3(Mst2)*pow4(Msq)*pow4(Mst1)*pow4(Mt) -
        316800*z2*pow2(Dmglst2)*pow3(Mst2)*pow4(Msq)*pow4(Mst1)*pow4(Mt) -
        5773600*Dmglst2*pow4(Msq)*pow4(Mst1)*pow4(Mst2)*pow4(Mt) - 2073600*
        Dmglst2*z2*pow4(Msq)*pow4(Mst1)*pow4(Mst2)*pow4(Mt) + 14400*pow11(Mst2)
        *pow2(Dmglst2)*pow4(Msq)*pow4(s2t) + 153900*pow11(Mst2)*pow2(Mst1)*
        pow4(Msq)*pow4(s2t) + 8545200*pow2(Dmglst2)*pow2(Mt)*pow2(s2t)*pow4(
        Msq)*pow4(Mst1)*pow5(Mst2) - 7603200*z2*pow2(Dmglst2)*pow2(Mt)*pow2(
        s2t)*pow4(Msq)*pow4(Mst1)*pow5(Mst2) + 10246400*Dmglst2*s2t*pow3(Mt)*
        pow4(Msq)*pow4(Mst1)*pow5(Mst2) - 2995200*Dmglst2*s2t*z2*pow3(Mt)*pow4(
        Msq)*pow4(Mst1)*pow5(Mst2) - 1656000*pow2(Dmglst2)*pow2(Mst1)*pow4(Msq)
        *pow4(Mt)*pow5(Mst2) + 741600*pow2(Dmglst2)*pow2(Msq)*pow4(Mst1)*pow4(
        Mt)*pow5(Mst2) - 640400*pow4(Msq)*pow4(Mst1)*pow4(Mt)*pow5(Mst2) +
        1936800*z2*pow4(Msq)*pow4(Mst1)*pow4(Mt)*pow5(Mst2) + 172800*z3*pow4(
        Msq)*pow4(Mst1)*pow4(Mt)*pow5(Mst2) - 3600*log(pow2(Msq)/pow2(Mst1))*
        pow3(Mt)*pow4(Mst1)*(6*Mst2*pow2(Dmglst2)*(2*(5*(25*Mt - 4*Mst2*s2t)*
        pow2(Mst1) - 6*Mt*pow2(Mst2))*pow4(Msq) + 20*Mt*pow2(Msq)*pow4(Mst2) +
        5*((Mt - 4*Mst2*s2t)*pow2(Mst1) + (11*Mt + 4*Mst2*s2t)*pow2(Mst2))*
        pow4(Mst2)) - 60*log(pow2(Mst2)/pow2(Mst1))*pow4(Msq)*(6*Mst2*Mt*pow2(
        Dmglst2)*pow2(Mst1) - 4*Dmglst2*(Mt*pow2(Mst1)*pow2(Mst2) + (9*Mt - 3*
        Mst2*s2t)*pow4(Mst1)) + Mst2*(2*Mt*pow2(Mst1)*pow2(Mst2) + (9*Mt - 4*
        Mst2*s2t)*pow4(Mst1) - 2*Mt*pow4(Mst2))) + 5*Mst2*(8*pow2(Msq)*((Mt -
        2*Mst2*s2t)*pow2(Mst1) + 2*(2*Mt + Mst2*s2t)*pow2(Mst2))*pow4(Mst2) +
        6*pow4(Msq)*(2*(5*Mt - 4*Mst2*s2t)*pow2(Mst1)*pow2(Mst2) + 2*(9*Mt - 4*
        Mst2*s2t)*pow4(Mst1) + (9*Mt + 4*Mst2*s2t)*pow4(Mst2)) + pow4(Mst2)*(
        pow2(Mst1)*(6*Mt*pow2(Mst2) - 8*s2t*pow3(Mst2)) + (3*Mt - 4*Mst2*s2t)*
        pow4(Mst1) + (19*Mt + 12*Mst2*s2t)*pow4(Mst2))) - 20*Dmglst2*(2*pow4(
        Msq)*(3*(7*Mt - 2*Mst2*s2t)*pow2(Mst1)*pow2(Mst2) + (81*Mt - 30*Mst2*
        s2t)*pow4(Mst1) + (-5*Mt + Mst2*s2t)*pow4(Mst2)) - 4*pow2(Msq)*(3*Mst2*
        Mt - s2t*pow2(Mst1) + s2t*pow2(Mst2))*pow5(Mst2) + (3*Mst2*(-Mt + 2*
        Mst2*s2t)*pow2(Mst1) - (13*Mt + 7*Mst2*s2t)*pow3(Mst2) + s2t*pow4(Mst1)
        )*pow5(Mst2))) - 432000*pow3(Mst2)*pow4(Mst1)*pow4(Mt)*pow6(Msq) -
        432000*pow2(Mt)*pow2(s2t)*pow4(Mst1)*pow5(Mst2)*pow6(Msq) - 432000*
        pow2(Mst1)*pow4(Mt)*pow5(Mst2)*pow6(Msq) + 9483600*pow2(Dmglst2)*pow2(
        Mt)*pow2(s2t)*pow3(Mst2)*pow4(Msq)*pow6(Mst1) - 6912000*z2*pow2(
        Dmglst2)*pow2(Mt)*pow2(s2t)*pow3(Mst2)*pow4(Msq)*pow6(Mst1) - 15360000*
        s2t*pow2(Dmglst2)*pow2(Mst2)*pow3(Mt)*pow4(Msq)*pow6(Mst1) - 5385600*
        s2t*z2*pow2(Dmglst2)*pow2(Mst2)*pow3(Mt)*pow4(Msq)*pow6(Mst1) +
        5241600*Dmglst2*s2t*pow3(Mst2)*pow3(Mt)*pow4(Msq)*pow6(Mst1) + 4089600*
        Dmglst2*s2t*z2*pow3(Mst2)*pow3(Mt)*pow4(Msq)*pow6(Mst1) - 288000*s2t*
        pow2(Dmglst2)*pow2(Msq)*pow3(Mt)*pow4(Mst2)*pow6(Mst1) + 4692000*
        Dmglst2*pow2(Mt)*pow2(s2t)*pow4(Msq)*pow4(Mst2)*pow6(Mst1) - 2764800*
        Dmglst2*z2*pow2(Mt)*pow2(s2t)*pow4(Msq)*pow4(Mst2)*pow6(Mst1) +
        4723200*s2t*pow3(Mt)*pow4(Msq)*pow4(Mst2)*pow6(Mst1) - 3686400*s2t*z2*
        pow3(Mt)*pow4(Msq)*pow4(Mst2)*pow6(Mst1) + 1399200*Mt*pow2(Dmglst2)*
        pow3(s2t)*pow4(Msq)*pow4(Mst2)*pow6(Mst1) - 144000*pow2(Dmglst2)*pow2(
        Msq)*pow3(Mst2)*pow4(Mt)*pow6(Mst1) + 29396000*Mst2*pow2(Dmglst2)*pow4(
        Msq)*pow4(Mt)*pow6(Mst1) + 16128000*Mst2*z2*pow2(Dmglst2)*pow4(Msq)*
        pow4(Mt)*pow6(Mst1) + 1228800*Dmglst2*pow2(Mst2)*pow4(Msq)*pow4(Mt)*
        pow6(Mst1) - 12096000*Dmglst2*z2*pow2(Mst2)*pow4(Msq)*pow4(Mt)*pow6(
        Mst1) - 2793600*pow3(Mst2)*pow4(Msq)*pow4(Mt)*pow6(Mst1) + 2160000*z2*
        pow3(Mst2)*pow4(Msq)*pow4(Mt)*pow6(Mst1) + 288000*Dmglst2*pow2(Msq)*
        pow4(Mst2)*pow4(Mt)*pow6(Mst1) - 336000*Dmglst2*s2t*pow2(Msq)*pow3(Mt)*
        pow5(Mst2)*pow6(Mst1) + 410400*pow2(Mt)*pow2(s2t)*pow4(Msq)*pow5(Mst2)*
        pow6(Mst1) - 1598400*Dmglst2*Mt*pow3(s2t)*pow4(Msq)*pow5(Mst2)*pow6(
        Mst1) + 244800*Dmglst2*Mt*z2*pow3(s2t)*pow4(Msq)*pow5(Mst2)*pow6(Mst1)
        + 279000*pow2(Dmglst2)*pow4(Mt)*pow5(Mst2)*pow6(Mst1) + 24000*pow2(Msq)
        *pow4(Mt)*pow5(Mst2)*pow6(Mst1) + 122700*pow2(Dmglst2)*pow4(Msq)*pow4(
        s2t)*pow5(Mst2)*pow6(Mst1) - 77400*z2*pow2(Dmglst2)*pow4(Msq)*pow4(s2t)
        *pow5(Mst2)*pow6(Mst1) + 216000*pow2(Mt)*pow2(s2t)*pow3(Mst2)*pow6(Msq)
        *pow6(Mst1) + 27000*pow4(s2t)*pow5(Mst2)*pow6(Msq)*pow6(Mst1) +
        1843200*s2t*pow2(Dmglst2)*pow2(Mst1)*pow3(Mt)*pow4(Msq)*pow6(Mst2) +
        288000*s2t*pow2(Dmglst2)*pow2(Msq)*pow3(Mt)*pow4(Mst1)*pow6(Mst2) +
        8556000*Dmglst2*pow2(Mt)*pow2(s2t)*pow4(Msq)*pow4(Mst1)*pow6(Mst2) -
        4147200*Dmglst2*z2*pow2(Mt)*pow2(s2t)*pow4(Msq)*pow4(Mst1)*pow6(Mst2) +
        662400*s2t*pow3(Mt)*pow4(Msq)*pow4(Mst1)*pow6(Mst2) - 3398400*s2t*z2*
        pow3(Mt)*pow4(Msq)*pow4(Mst1)*pow6(Mst2) + 2100000*Mt*pow2(Dmglst2)*
        pow3(s2t)*pow4(Msq)*pow4(Mst1)*pow6(Mst2) - 1555200*Mt*z2*pow2(Dmglst2)
        *pow3(s2t)*pow4(Msq)*pow4(Mst1)*pow6(Mst2) + 1008000*Dmglst2*pow2(Mst1)
        *pow4(Msq)*pow4(Mt)*pow6(Mst2) + 72000*Dmglst2*pow2(Msq)*pow4(Mst1)*
        pow4(Mt)*pow6(Mst2) - 468000*s2t*pow2(Dmglst2)*pow3(Mt)*pow6(Mst1)*
        pow6(Mst2) + 240000*s2t*pow2(Msq)*pow3(Mt)*pow6(Mst1)*pow6(Mst2) -
        2131200*Mt*pow3(s2t)*pow4(Msq)*pow6(Mst1)*pow6(Mst2) + 244800*Mt*z2*
        pow3(s2t)*pow4(Msq)*pow6(Mst1)*pow6(Mst2) + 126000*Dmglst2*pow4(Mt)*
        pow6(Mst1)*pow6(Mst2) + 1800*Dmglst2*pow4(Msq)*pow4(s2t)*pow6(Mst1)*
        pow6(Mst2) - 154800*Dmglst2*z2*pow4(Msq)*pow4(s2t)*pow6(Mst1)*pow6(
        Mst2) + 450*pow2(log(pow2(Mst2)/pow2(Mst1)))*pow4(Msq)*pow4(Mst1)*(-8*
        Dmglst2*(8*Mt*(-67*Mst2*s2t*pow2(Mt) - 108*Mt*pow2(Mst2)*pow2(s2t) +
        128*pow3(Mt) - 16*pow3(Mst2)*pow3(s2t))*pow4(Mst1) + pow4(Mst2)*(-320*
        pow2(Mst2)*pow2(Mt)*pow2(s2t) + 848*Mst2*s2t*pow3(Mt) - 532*Mt*pow3(
        Mst2)*pow3(s2t) + 512*pow4(Mt) - 99*pow4(Mst2)*pow4(s2t)) + pow2(Mst1)*
        pow2(Mst2)*(-768*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 448*Mst2*s2t*pow3(Mt)
        - 212*Mt*pow3(Mst2)*pow3(s2t) + 1656*pow4(Mt) + 35*pow4(Mst2)*pow4(s2t)
        )) + Mst2*(8*pow2(Mst1)*pow2(Mst2)*(425*pow2(Mst2)*pow2(Mt)*pow2(s2t) -
        704*Mst2*s2t*pow3(Mt) + 212*Mt*pow3(Mst2)*pow3(s2t) - 196*pow4(Mt) + 3*
        pow4(Mst2)*pow4(s2t)) + pow4(Mst1)*(2688*pow2(Mst2)*pow2(Mt)*pow2(s2t)
        - 3136*Mst2*s2t*pow3(Mt) + 1024*Mt*pow3(Mst2)*pow3(s2t) - 7168*pow4(Mt)
        + 41*pow4(Mst2)*pow4(s2t)) + pow4(Mst2)*(296*pow2(Mst2)*pow2(Mt)*pow2(
        s2t) - 6784*Mst2*s2t*pow3(Mt) + 2208*Mt*pow3(Mst2)*pow3(s2t) + 288*
        pow4(Mt) + 519*pow4(Mst2)*pow4(s2t))) + 4*Mst2*pow2(Dmglst2)*(320*pow2(
        Mt)*pow2(s2t)*pow4(Mst2) - 512*pow2(Mst2)*pow4(Mt) + pow2(Mst1)*(768*
        pow2(Mst2)*pow2(Mt)*pow2(s2t) - 1152*Mst2*s2t*pow3(Mt) + 2920*pow4(Mt)
        - 35*pow4(Mst2)*pow4(s2t)) + 768*Mt*pow3(s2t)*pow5(Mst2) + 99*pow4(s2t)
        *pow6(Mst2))) + 597600*pow2(Dmglst2)*pow2(Mst1)*pow2(Mt)*pow2(s2t)*
        pow4(Msq)*pow7(Mst2) - 921600*Dmglst2*s2t*pow2(Mst1)*pow3(Mt)*pow4(Msq)
        *pow7(Mst2) + 288000*Dmglst2*s2t*pow2(Msq)*pow3(Mt)*pow4(Mst1)*pow7(
        Mst2) + 3636000*pow2(Mt)*pow2(s2t)*pow4(Msq)*pow4(Mst1)*pow7(Mst2) -
        691200*z2*pow2(Mt)*pow2(s2t)*pow4(Msq)*pow4(Mst1)*pow7(Mst2) + 2764800*
        Dmglst2*Mt*pow3(s2t)*pow4(Msq)*pow4(Mst1)*pow7(Mst2) - 1310400*Dmglst2*
        Mt*z2*pow3(s2t)*pow4(Msq)*pow4(Mst1)*pow7(Mst2) + 230400*pow2(Dmglst2)*
        pow4(Msq)*pow4(Mt)*pow7(Mst2) + 1540800*pow2(Mst1)*pow4(Msq)*pow4(Mt)*
        pow7(Mst2) + 779400*pow2(Dmglst2)*pow4(Mst1)*pow4(Mt)*pow7(Mst2) -
        468000*pow2(Msq)*pow4(Mst1)*pow4(Mt)*pow7(Mst2) - 123000*pow2(Dmglst2)*
        pow4(Msq)*pow4(Mst1)*pow4(s2t)*pow7(Mst2) + 77400*z2*pow2(Dmglst2)*
        pow4(Msq)*pow4(Mst1)*pow4(s2t)*pow7(Mst2) + 216000*pow2(Mst1)*pow2(Mt)*
        pow2(s2t)*pow6(Msq)*pow7(Mst2) + 27000*pow4(Mst1)*pow4(s2t)*pow6(Msq)*
        pow7(Mst2) - 36000*Dmglst2*s2t*pow3(Mt)*pow6(Mst1)*pow7(Mst2) - 45000*
        pow4(Mt)*pow6(Mst1)*pow7(Mst2) - 522900*pow4(Msq)*pow4(s2t)*pow6(Mst1)*
        pow7(Mst2) - 77400*z2*pow4(Msq)*pow4(s2t)*pow6(Mst1)*pow7(Mst2) +
        288000*Dmglst2*s2t*pow2(Msq)*pow3(Mst2)*pow3(Mt)*pow8(Mst1) + 4435200*
        Dmglst2*pow2(Mst2)*pow2(Mt)*pow2(s2t)*pow4(Msq)*pow8(Mst1) - 2764800*
        Dmglst2*z2*pow2(Mst2)*pow2(Mt)*pow2(s2t)*pow4(Msq)*pow8(Mst1) + 194400*
        pow2(Mt)*pow2(s2t)*pow3(Mst2)*pow4(Msq)*pow8(Mst1) - 7084800*Dmglst2*
        Mst2*s2t*pow3(Mt)*pow4(Msq)*pow8(Mst1) + 11116800*Dmglst2*Mst2*s2t*z2*
        pow3(Mt)*pow4(Msq)*pow8(Mst1) + 5860800*s2t*pow2(Mst2)*pow3(Mt)*pow4(
        Msq)*pow8(Mst1) - 3571200*s2t*z2*pow2(Mst2)*pow3(Mt)*pow4(Msq)*pow8(
        Mst1) + 82800*Dmglst2*Mt*pow3(Mst2)*pow3(s2t)*pow4(Msq)*pow8(Mst1) -
        288000*s2t*pow2(Msq)*pow3(Mt)*pow4(Mst2)*pow8(Mst1) - 176400*Mt*pow3(
        s2t)*pow4(Msq)*pow4(Mst2)*pow8(Mst1) - 1440000*Dmglst2*pow2(Msq)*pow2(
        Mst2)*pow4(Mt)*pow8(Mst1) + 504000*pow2(Msq)*pow3(Mst2)*pow4(Mt)*pow8(
        Mst1) + 23808000*Dmglst2*pow4(Msq)*pow4(Mt)*pow8(Mst1) - 4991400*Mst2*
        pow4(Msq)*pow4(Mt)*pow8(Mst1) - 27302400*Dmglst2*z2*pow4(Msq)*pow4(Mt)*
        pow8(Mst1) + 2534400*Mst2*z2*pow4(Msq)*pow4(Mt)*pow8(Mst1) + 108000*
        Dmglst2*pow4(Mst2)*pow4(Mt)*pow8(Mst1) + 132300*Dmglst2*pow4(Msq)*pow4(
        Mst2)*pow4(s2t)*pow8(Mst1) - 138000*Dmglst2*s2t*pow3(Mt)*pow5(Mst2)*
        pow8(Mst1) + 49500*pow4(Mt)*pow5(Mst2)*pow8(Mst1) + 375075*pow4(Msq)*
        pow4(s2t)*pow5(Mst2)*pow8(Mst1) + 36900*z2*pow4(Msq)*pow4(s2t)*pow5(
        Mst2)*pow8(Mst1) - 27000*pow3(Mst2)*pow4(s2t)*pow6(Msq)*pow8(Mst1) +
        6000*s2t*pow3(Mt)*pow6(Mst2)*pow8(Mst1) - 964800*Dmglst2*pow2(Mst1)*
        pow2(Mt)*pow2(s2t)*pow4(Msq)*pow8(Mst2) - 921600*s2t*pow2(Mst1)*pow3(
        Mt)*pow4(Msq)*pow8(Mst2) - 460800*Mt*pow2(Dmglst2)*pow2(Mst1)*pow3(s2t)
        *pow4(Msq)*pow8(Mst2) + 540000*s2t*pow2(Dmglst2)*pow3(Mt)*pow4(Mst1)*
        pow8(Mst2) - 96000*s2t*pow2(Msq)*pow3(Mt)*pow4(Mst1)*pow8(Mst2) +
        2174400*Mt*pow3(s2t)*pow4(Msq)*pow4(Mst1)*pow8(Mst2) - 273600*Mt*z2*
        pow3(s2t)*pow4(Msq)*pow4(Mst1)*pow8(Mst2) - 460800*Dmglst2*pow4(Msq)*
        pow4(Mt)*pow8(Mst2) + 78000*Dmglst2*pow4(Mst1)*pow4(Mt)*pow8(Mst2) -
        379800*Dmglst2*pow4(Msq)*pow4(Mst1)*pow4(s2t)*pow8(Mst2) + 154800*
        Dmglst2*z2*pow4(Msq)*pow4(Mst1)*pow4(s2t)*pow8(Mst2) + 84000*s2t*pow3(
        Mt)*pow6(Mst1)*pow8(Mst2) - 30*log(pow2(Mst2)/pow2(Mst1))*(-2*Mst2*
        pow2(Dmglst2)*(-1800*(-((Mt - 4*Mst2*s2t)*pow2(Mst1)) - (11*Mt + 4*
        Mst2*s2t)*pow2(Mst2))*pow3(Mt)*pow4(Mst1)*pow4(Mst2) + 7200*pow2(Msq)*
        pow4(Mst1)*pow4(Mst2)*pow4(Mt) + pow4(Msq)*(-45*pow2(Mst1)*pow4(Mst2)*(
        -456*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 1024*Mst2*s2t*pow3(Mt) + 256*Mt*
        pow3(Mst2)*pow3(s2t) + 400*pow4(Mt) + 153*pow4(Mst2)*pow4(s2t)) + 2*
        pow2(Mst2)*pow4(Mst1)*(-3240*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 45920*
        Mst2*s2t*pow3(Mt) + 9720*Mt*pow3(Mst2)*pow3(s2t) - 39568*pow4(Mt) +
        4635*pow4(Mst2)*pow4(s2t)) + (68160*pow2(Mst2)*pow2(Mt)*pow2(s2t) +
        280720*Mst2*s2t*pow3(Mt) + 11760*Mt*pow3(Mst2)*pow3(s2t) - 741864*pow4(
        Mt) - 2010*pow4(Mst2)*pow4(s2t))*pow6(Mst1) + 1440*pow2(-4*pow2(Mt) +
        pow2(Mst2)*pow2(s2t))*pow6(Mst2))) - 20*Dmglst2*(240*pow2(Msq)*pow2(
        Mst2)*pow3(Mt)*pow4(Mst1)*(-2*s2t*pow2(Mst1)*pow3(Mst2) + 3*Mt*pow4(
        Mst1) + 2*(3*Mt + Mst2*s2t)*pow4(Mst2)) - 120*pow3(Mt)*pow4(Mst1)*(3*
        Mst2*(-Mt + 2*Mst2*s2t)*pow2(Mst1) - (13*Mt + 7*Mst2*s2t)*pow3(Mst2) +
        s2t*pow4(Mst1))*pow5(Mst2) + pow4(Msq)*(-4*pow4(Mst1)*pow4(Mst2)*(1860*
        pow2(Mst2)*pow2(Mt)*pow2(s2t) + 1532*Mst2*s2t*pow3(Mt) + 774*Mt*pow3(
        Mst2)*pow3(s2t) - 2764*pow4(Mt) + 15*pow4(Mst2)*pow4(s2t)) + 2*pow2(
        Mst2)*(-432*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 10936*Mst2*s2t*pow3(Mt) +
        996*Mt*pow3(Mst2)*pow3(s2t) + 16228*pow4(Mt) + 207*pow4(Mst2)*pow4(s2t)
        )*pow6(Mst1) - 3*pow2(Mst1)*(-856*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 1536*
        Mst2*s2t*pow3(Mt) + 384*Mt*pow3(Mst2)*pow3(s2t) + 1200*pow4(Mt) + 203*
        pow4(Mst2)*pow4(s2t))*pow6(Mst2) + 4*(-738*pow2(Mst2)*pow2(Mt)*pow2(
        s2t) - 8204*Mst2*s2t*pow3(Mt) + 135*Mt*pow3(Mst2)*pow3(s2t) + 12492*
        pow4(Mt) + 12*pow4(Mst2)*pow4(s2t))*pow8(Mst1) - 768*pow2(Mt)*(-2*pow2(
        Mt) + pow2(Mst2)*pow2(s2t))*pow8(Mst2))) + 5*Mst2*(480*pow2(Msq)*pow2(
        Mst2)*pow3(Mt)*pow4(Mst1)*(2*(-Mt + 2*Mst2*s2t)*pow2(Mst1)*pow2(Mst2) +
        3*Mt*pow4(Mst1) - 4*(2*Mt + Mst2*s2t)*pow4(Mst2)) + 120*pow3(Mt)*pow4(
        Mst1)*pow4(Mst2)*(2*(-3*Mt + 4*Mst2*s2t)*pow2(Mst1)*pow2(Mst2) + (-3*Mt
        + 4*Mst2*s2t)*pow4(Mst1) - (19*Mt + 12*Mst2*s2t)*pow4(Mst2)) + 360*(
        pow2(Mst1) - pow2(Mst2))*pow4(Mst1)*pow4(Mst2)*pow4(s2t)*pow6(Msq) +
        pow4(Msq)*(4*pow4(Mst1)*pow4(Mst2)*(3372*pow2(Mst2)*pow2(Mt)*pow2(s2t)
        - 9744*Mst2*s2t*pow3(Mt) + 3912*Mt*pow3(Mst2)*pow3(s2t) + 4492*pow4(Mt)
        + 669*pow4(Mst2)*pow4(s2t)) - 48*pow2(Mst2)*(-267*pow2(Mst2)*pow2(Mt)*
        pow2(s2t) + 792*Mst2*s2t*pow3(Mt) - 30*Mt*pow3(Mst2)*pow3(s2t) - 317*
        pow4(Mt) + 56*pow4(Mst2)*pow4(s2t))*pow6(Mst1) + 6*pow2(Mst1)*(-728*
        pow2(Mst2)*pow2(Mt)*pow2(s2t) - 1024*Mst2*s2t*pow3(Mt) + 256*Mt*pow3(
        Mst2)*pow3(s2t) + 1200*pow4(Mt) + 139*pow4(Mst2)*pow4(s2t))*pow6(Mst2)
        + 3*(4640*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 13376*Mst2*s2t*pow3(Mt) +
        240*Mt*pow3(Mst2)*pow3(s2t) + 2328*pow4(Mt) - 279*pow4(Mst2)*pow4(s2t))
        *pow8(Mst1) + 768*pow2(Mt)*(-2*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow8(
        Mst2)))) - 115200*pow2(Dmglst2)*pow2(Mt)*pow2(s2t)*pow4(Msq)*pow9(Mst2)
        - 885600*pow2(Mst1)*pow2(Mt)*pow2(s2t)*pow4(Msq)*pow9(Mst2) + 230400*
        Dmglst2*Mt*pow2(Mst1)*pow3(s2t)*pow4(Msq)*pow9(Mst2) + 162000*Dmglst2*
        s2t*pow3(Mt)*pow4(Mst1)*pow9(Mst2) - 230400*pow4(Msq)*pow4(Mt)*pow9(
        Mst2) - 115500*pow4(Mst1)*pow4(Mt)*pow9(Mst2) - 74700*pow2(Dmglst2)*
        pow2(Mst1)*pow4(Msq)*pow4(s2t)*pow9(Mst2) + 37800*pow4(Msq)*pow4(Mst1)*
        pow4(s2t)*pow9(Mst2) + 40500*z2*pow4(Msq)*pow4(Mst1)*pow4(s2t)*pow9(
        Mst2) - 27000*pow2(Mst1)*pow4(s2t)*pow6(Msq)*pow9(Mst2) + 230400*
        Dmglst2*pow2(Mt)*pow2(s2t)*pow4(Msq)*power10(Mst2) + 230400*Mt*pow2(
        Mst1)*pow3(s2t)*pow4(Msq)*power10(Mst2) - 54000*s2t*pow3(Mt)*pow4(Mst1)
        *power10(Mst2) + 235800*Dmglst2*pow2(Mst1)*pow4(Msq)*pow4(s2t)*power10(
        Mst2)))/(16200.*pow4(Msq)*pow4(Mst1)*pow5(Mst2)))/pow4(Mt)*
        12.; 

   return result;
}

double himalaya::H6::calc_coef_at_as2_no_sm_logs_log2() const {

   const double result =
      (((1920*pow11(Mst2)*pow2(Mt)*pow2(s2t) + 52480*s2t*pow2(
        Dmglst2)*pow3(Mt)*pow4(Mst1)*pow4(Mst2) + 18896*pow2(Dmglst2)*pow3(
        Mst2)*pow4(Mst1)*pow4(Mt) - 48160*Dmglst2*pow4(Mst1)*pow4(Mst2)*pow4(
        Mt) - 1440*pow11(Mst2)*pow2(Dmglst2)*pow4(s2t) + 2325*pow11(Mst2)*pow2(
        Mst1)*pow4(s2t) + 33360*pow2(Dmglst2)*pow2(Mt)*pow2(s2t)*pow4(Mst1)*
        pow5(Mst2) + 46400*Dmglst2*s2t*pow3(Mt)*pow4(Mst1)*pow5(Mst2) + 21840*
        pow2(Dmglst2)*pow2(Mst1)*pow4(Mt)*pow5(Mst2) + 21280*pow4(Mst1)*pow4(
        Mt)*pow5(Mst2) + 7200*log(pow2(Msq)/pow2(Mst1))*pow4(Mst1)*pow4(Mt)*
        pow5(Mst2) + 600*pow2(Dmglst2)*pow2(Mt)*pow2(s2t)*pow3(Mst2)*pow6(Mst1)
        - 113280*s2t*pow2(Dmglst2)*pow2(Mst2)*pow3(Mt)*pow6(Mst1) + 74880*
        Dmglst2*s2t*pow3(Mst2)*pow3(Mt)*pow6(Mst1) - 14160*Dmglst2*pow2(Mt)*
        pow2(s2t)*pow4(Mst2)*pow6(Mst1) + 1920*s2t*pow3(Mt)*pow4(Mst2)*pow6(
        Mst1) - 11520*Mt*pow2(Dmglst2)*pow3(s2t)*pow4(Mst2)*pow6(Mst1) +
        261920*Mst2*pow2(Dmglst2)*pow4(Mt)*pow6(Mst1) - 100800*Dmglst2*pow2(
        Mst2)*pow4(Mt)*pow6(Mst1) - 10080*pow3(Mst2)*pow4(Mt)*pow6(Mst1) -
        9000*pow2(Mt)*pow2(s2t)*pow5(Mst2)*pow6(Mst1) - 31200*Dmglst2*Mt*pow3(
        s2t)*pow5(Mst2)*pow6(Mst1) + 4395*pow2(Dmglst2)*pow4(s2t)*pow5(Mst2)*
        pow6(Mst1) - 46080*s2t*pow2(Dmglst2)*pow2(Mst1)*pow3(Mt)*pow6(Mst2) +
        82080*Dmglst2*pow2(Mt)*pow2(s2t)*pow4(Mst1)*pow6(Mst2) - 24000*s2t*
        pow3(Mt)*pow4(Mst1)*pow6(Mst2) + 43680*Dmglst2*pow2(Mst1)*pow4(Mt)*
        pow6(Mst2) - 23520*Mt*pow3(s2t)*pow6(Mst1)*pow6(Mst2) + 1110*Dmglst2*
        pow4(s2t)*pow6(Mst1)*pow6(Mst2) - 30*log(pow2(Mst2)/pow2(Mst1))*pow4(
        Mst1)*(Mst2*(pow2(Mst1)*pow2(Mst2)*(1096*pow2(Mst2)*pow2(Mt)*pow2(s2t)
        - 1280*Mst2*s2t*pow3(Mt) + 328*Mt*pow3(Mst2)*pow3(s2t) - 32*pow4(Mt) -
        91*pow4(Mst2)*pow4(s2t)) + 4*pow4(Mst2)*(14*pow2(Mst2)*pow2(Mt)*pow2(
        s2t) - 456*Mst2*s2t*pow3(Mt) + 146*Mt*pow3(Mst2)*pow3(s2t) - 106*pow4(
        Mt) + 33*pow4(Mst2)*pow4(s2t)) - pow4(Mst1)*(-768*pow2(Mst2)*pow2(Mt)*
        pow2(s2t) + 1216*Mst2*s2t*pow3(Mt) + 400*pow4(Mt) + 41*pow4(Mst2)*pow4(
        s2t))) - 2*Dmglst2*(32*pow2(Mt)*(-57*Mst2*Mt*s2t + 119*pow2(Mt) - 24*
        pow2(Mst2)*pow2(s2t))*pow4(Mst1) + pow4(Mst2)*(-384*pow2(Mst2)*pow2(Mt)
        *pow2(s2t) + 912*Mst2*s2t*pow3(Mt) - 548*Mt*pow3(Mst2)*pow3(s2t) + 512*
        pow4(Mt) - 91*pow4(Mst2)*pow4(s2t)) + pow2(Mst1)*pow2(Mst2)*(-768*pow2(
        Mst2)*pow2(Mt)*pow2(s2t) - 640*Mst2*s2t*pow3(Mt) - 164*Mt*pow3(Mst2)*
        pow3(s2t) + 2016*pow4(Mt) + 91*pow4(Mst2)*pow4(s2t))) + Mst2*pow2(
        Dmglst2)*(384*pow2(Mt)*pow2(s2t)*pow4(Mst2) - 512*pow2(Mst2)*pow4(Mt) +
        pow2(Mst1)*(768*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 1280*Mst2*s2t*pow3(Mt)
        + 4000*pow4(Mt) - 91*pow4(Mst2)*pow4(s2t)) + 768*Mt*pow3(s2t)*pow5(
        Mst2) + 91*pow4(s2t)*pow6(Mst2))) - 22440*pow2(Dmglst2)*pow2(Mst1)*
        pow2(Mt)*pow2(s2t)*pow7(Mst2) - 46080*Dmglst2*s2t*pow2(Mst1)*pow3(Mt)*
        pow7(Mst2) + 42960*pow2(Mt)*pow2(s2t)*pow4(Mst1)*pow7(Mst2) + 19680*
        Dmglst2*Mt*pow3(s2t)*pow4(Mst1)*pow7(Mst2) - 23040*pow2(Dmglst2)*pow4(
        Mt)*pow7(Mst2) + 21840*pow2(Mst1)*pow4(Mt)*pow7(Mst2) - 10005*pow2(
        Dmglst2)*pow4(Mst1)*pow4(s2t)*pow7(Mst2) - 5325*pow4(s2t)*pow6(Mst1)*
        pow7(Mst2) + 67200*Dmglst2*Mst2*s2t*pow3(Mt)*pow8(Mst1) + 1920*s2t*
        pow2(Mst2)*pow3(Mt)*pow8(Mst1) - 112320*Dmglst2*pow4(Mt)*pow8(Mst1) -
        12000*Mst2*pow4(Mt)*pow8(Mst1) + 1770*Dmglst2*pow4(Mst2)*pow4(s2t)*
        pow8(Mst1) + 3585*pow4(s2t)*pow5(Mst2)*pow8(Mst1) - 29520*Dmglst2*pow2(
        Mst1)*pow2(Mt)*pow2(s2t)*pow8(Mst2) - 15360*s2t*pow2(Mst1)*pow3(Mt)*
        pow8(Mst2) + 11520*Mt*pow2(Dmglst2)*pow2(Mst1)*pow3(s2t)*pow8(Mst2) +
        19680*Mt*pow3(s2t)*pow4(Mst1)*pow8(Mst2) - 15360*Dmglst2*pow4(Mt)*pow8(
        Mst2) - 8490*Dmglst2*pow4(Mst1)*pow4(s2t)*pow8(Mst2) + 11520*pow2(
        Dmglst2)*pow2(Mt)*pow2(s2t)*pow9(Mst2) - 12840*pow2(Mst1)*pow2(Mt)*
        pow2(s2t)*pow9(Mst2) + 11520*Dmglst2*Mt*pow2(Mst1)*pow3(s2t)*pow9(Mst2)
        - 3840*pow4(Mt)*pow9(Mst2) + 7125*pow2(Dmglst2)*pow2(Mst1)*pow4(s2t)*
        pow9(Mst2) - 345*pow4(Mst1)*pow4(s2t)*pow9(Mst2) + 7680*Dmglst2*pow2(
        Mt)*pow2(s2t)*power10(Mst2) + 3840*Mt*pow2(Mst1)*pow3(s2t)*power10(
        Mst2) + 6570*Dmglst2*pow2(Mst1)*pow4(s2t)*power10(Mst2)))/(540.*pow4(
        Mst1)*pow5(Mst2)))/pow4(Mt)*12.; 

   return result;
}

double himalaya::H6::calc_coef_at_as2_no_sm_logs_log3() const {

   const double result =
      ((-224*pow4(Mt))/9.)/pow4(Mt)*12.; 

   return result;
}




