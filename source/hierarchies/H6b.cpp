#include "H6b.hpp"
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
himalaya::H6b::H6b(const std::map<unsigned int, unsigned int>& flagMap, double Al4p, double beta, double Dmglst2,
		 double Dmsqst2, double lmMt, double lmMst1, double lmMst2,
		 double Mt, double Mst1, double Mst2, double MuSUSY,
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
   #include "hierarchies/h6b/sigS1Full.inc"
   ;
   s2 = 
   #include "hierarchies/h6b/sigS2Full.inc"
   ;
   s12 = 
   #include "hierarchies/h6b/sigS12Full.inc"
   ;
}

/**
 * 	@return The diagonal (1, 1) matrix element of the Higgs mass matrix as a double for the hierarchy 'H6b'
 */
double himalaya::H6b::getS1() const {
   return s1;
}

/**
 * 	@return The diagonal (2, 2) matrix element of the Higgs mass matrix as a double for the hierarchy 'H6b'
 */
double himalaya::H6b::getS2() const {
   return s2;
}

/**
 * 	@return The off-diagonal (1, 2) = (2, 1) matrix element of the Higgs mass matrix as a double for the hierarchy 'H6b'
 */
double himalaya::H6b::getS12() const {
   return s12;
}

/**
 * 	@return returns the susy log^0 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double himalaya::H6b::calc_coef_at_as2_no_sm_logs_log0() const {

   const double result =
      (-(pow2(Sbeta)*(80015040000*pow2(Mst2)*pow2(Mt)*pow3(Dmsqst2) + (30*
        Dmsqst2*(Dmsqst2*(-166698000*(-2 + z2)*pow2(-4*pow2(Mt) + pow2(Mst2)*
        pow2(s2t))*pow4(Mst2) - 385875*pow2(Mst1)*pow2(Mst2)*(-8*pow2(Mst2)*
        pow2(Mt)*pow2(s2t)*(160 - 64*OepS2 + 648*S2 + 96*T1ep + 2076*z2 - 903*
        z3 + 48*z4 + 144*(2 + 9*S2)*log(pow2(Mst2)/pow2(Mst1)) + 72*pow2(z2)) +
        3456*Mst2*s2t*(15 - 20*z2 + 16*z3 + 2*(-15 + 8*z2)*log(pow2(Mst2)/pow2(
        Mst1)))*pow3(Mt) + 6912*Mt*(3 - 3*z2 + z3 - (-2 + z2)*log(pow2(Mst2)/
        pow2(Mst1)))*pow3(Mst2)*pow3(s2t) - 16*(6190 - 5616*z2 + 2025*z3 + 12*(
        -365 + 216*z2)*log(pow2(Mst2)/pow2(Mst1)) + 216*pow2(log(pow2(Mst2)/
        pow2(Mst1))))*pow4(Mt) - (2122 + 128*OepS2 - 55728*S2 - 192*T1ep -
        4152*z2 + 13209*z3 - 96*z4 - 144*(-1 + 18*S2 + 3*z2)*log(pow2(Mst2)/
        pow2(Mst1)) - 144*pow2(z2))*pow4(Mst2)*pow4(s2t)) - 2058*pow4(Mst1)*(4*
        pow2(Mst2)*pow2(Mt)*pow2(s2t)*(487028 + 40000*OepS2 - 1917000*S2 -
        60000*T1ep - 1081500*z2 - 2625*z3 - 30000*z4 - 45000*pow2(z2)) -
        2592000*Mst2*s2t*(-11 + 14*z2 - 8*z3)*pow3(Mt) + 1296000*Mt*(-2 + z2 +
        z3)*pow3(Mst2)*pow3(s2t) + 288*(-176473 + 202500*z2 - 90000*z3)*pow4(
        Mt) + (1932736 - 64000*OepS2 + 3294000*S2 + 96000*T1ep + 726000*z2 -
        710625*z3 + 48000*z4 + 72000*pow2(z2))*pow4(Mst2)*pow4(s2t) + 32400*
        pow2(log(pow2(Mst2)/pow2(Mst1)))*(14*pow2(Mst2)*pow2(Mt)*pow2(s2t) +
        32*pow4(Mt) - 3*pow4(Mst2)*pow4(s2t)) - 180*log(pow2(Mst2)/pow2(Mst1))*
        (8*(-887 + 2250*S2)*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 57600*Mst2*s2t*(-3
        + 2*z2)*pow3(Mt) + 7200*Mt*(-2 + z2)*pow3(Mst2)*pow3(s2t) + 48*(-3871 +
        3000*z2)*pow4(Mt) + (1837 - 7200*S2 - 900*z2)*pow4(Mst2)*pow4(s2t))) +
        5*pow2(s2t)*(33339600*Mst2*Mt*s2t*(69 - 90*z2 + 60*z3 + (121 - 60*z2)*
        log(pow2(Mst2)/pow2(Mst1))) + 600*pow2(Mt)*(-109760*OepS2 + 63*(-20869
        + 35280*S2 + 5292*z2)*log(pow2(Mst2)/pow2(Mst1)) + 3*(1281429 +
        1720488*S2 + 54880*T1ep + 44590*z2 - 303212*z3 + 27440*z4 + 41160*pow2(
        z2)) - 285768*pow2(log(pow2(Mst2)/pow2(Mst1)))) + pow2(Mst2)*pow2(s2t)*
        (76781738 + 15366400*OepS2 - 837194400*S2 - 23049600*T1ep - 207652200*
        z2 - 63103425*z3 - 11524800*z4 - 1260*(8753 + 246960*S2 - 79380*z2)*
        log(pow2(Mst2)/pow2(Mst1)) - 17287200*pow2(z2) + 71442000*pow2(log(
        pow2(Mst2)/pow2(Mst1)))))*pow6(Mst1)) - 4*(-385875*pow2(Mst1)*pow3(
        Mst2)*(16*Dmglst2*(36*(23 - 24*z2)*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 24*
        Mst2*s2t*(289 - 342*z2 + 216*z3 + 6*(-59 + 36*z2)*log(pow2(Mst2)/pow2(
        Mst1)))*pow3(Mt) + 216*Mt*(2 - 2*z2 + z3 + (5 - 3*z2)*log(pow2(Mst2)/
        pow2(Mst1)))*pow3(Mst2)*pow3(s2t) + 4*(-2954 + 2484*z2 - 1080*z3 - 3*(-
        571 + 360*z2)*log(pow2(Mst2)/pow2(Mst1)))*pow4(Mt) - 27*(2 - 4*z2 + 3*
        z3 + (2 - 3*z2)*log(pow2(Mst2)/pow2(Mst1)))*pow4(Mst2)*pow4(s2t)) + 3*
        Mst2*(-4*pow2(Mst2)*pow2(Mt)*pow2(s2t)*(-32*OepS2 + 24*(16 + 27*S2)*
        log(pow2(Mst2)/pow2(Mst1)) + 3*(391 + 972*S2 + 16*T1ep + 154*z2 - 232*
        z3 + 8*z4 + 12*pow2(z2))) + 1152*Mst2*s2t*(7 + 2*z2 - 8*z3 + (10 - 8*
        z2)*log(pow2(Mst2)/pow2(Mst1)))*pow3(Mt) + 1152*Mt*(-1 + 2*z2 - 3*z3 +
        (-1 + z2)*log(pow2(Mst2)/pow2(Mst1)))*pow3(Mst2)*pow3(s2t) + 64*(103 -
        108*z2 + 72*z3 + 6*(-23 + 12*z2)*log(pow2(Mst2)/pow2(Mst1)) - 9*pow2(
        log(pow2(Mst2)/pow2(Mst1))))*pow4(Mt) + (1213 - 32*OepS2 - 4860*S2 +
        48*T1ep + 462*z2 + 276*z3 + 24*z4 + 24*(17 + 27*S2)*log(pow2(Mst2)/
        pow2(Mst1)) + 36*pow2(z2) + 144*pow2(log(pow2(Mst2)/pow2(Mst1))))*pow4(
        Mst2)*pow4(s2t))) - 3087*Mst2*pow4(Mst1)*(4500*Dmglst2*(-8*(-113 + 48*
        z2 + 6*log(pow2(Mst2)/pow2(Mst1)))*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 192*
        Mst2*s2t*(-51*z2 + 28*(2 + z3) + 4*(-12 + 7*z2)*log(pow2(Mst2)/pow2(
        Mst1)))*pow3(Mt) + 96*Mt*(5*z2 + 3*(-4 + z3) + (5 - 3*z2)*log(pow2(
        Mst2)/pow2(Mst1)))*pow3(Mst2)*pow3(s2t) - 48*(401 - 308*z2 + 144*z3 +
        2*(-137 + 72*z2)*log(pow2(Mst2)/pow2(Mst1)) + 4*pow2(log(pow2(Mst2)/
        pow2(Mst1))))*pow4(Mt) + 3*(75 - 32*z2 - 8*z3 + (-26 + 8*z2)*log(pow2(
        Mst2)/pow2(Mst1)))*pow4(Mst2)*pow4(s2t)) + Mst2*(8*pow2(Mst2)*pow2(Mt)*
        pow2(s2t)*(36449 + 4000*OepS2 - 310500*S2 - 6000*T1ep - 75750*z2 +
        21000*z3 - 3000*z4 - 4500*pow2(z2)) + 864000*Mst2*s2t*(-6 + 5*z2 - 4*
        z3)*pow3(Mt) + 432000*Mt*(1 + z2 - z3)*pow3(Mst2)*pow3(s2t) - 96*(-
        91321 + 54000*z2 - 36000*z3)*pow4(Mt) + (-52199 - 28000*OepS2 +
        3415500*S2 + 42000*T1ep + 557250*z2 - 259500*z3 + 21000*z4 + 31500*
        pow2(z2))*pow4(Mst2)*pow4(s2t) + 2700*pow2(log(pow2(Mst2)/pow2(Mst1)))*
        (48*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 144*pow4(Mt) - 31*pow4(Mst2)*pow4(
        s2t)) + 180*log(pow2(Mst2)/pow2(Mst1))*(-8*(-517 + 450*S2)*pow2(Mst2)*
        pow2(Mt)*pow2(s2t) - 19200*Mst2*s2t*(-2 + z2)*pow3(Mt) + 2400*Mt*(-1 +
        z2)*pow3(Mst2)*pow3(s2t) + 16*(-3067 + 1200*z2)*pow4(Mt) + (-317 +
        3150*S2 - 300*z2)*pow4(Mst2)*pow4(s2t)))) + 41674500*(Mst2 - 4*Dmglst2*
        (-3 + z2) + 2*Mst2*z2 + 2*Mst2*log(pow2(Mst2)/pow2(Mst1)))*pow2(-4*
        pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow5(Mst2) + 2*(-4*Mst2*pow2(Mt)*pow2(
        s2t)*(1157625*Dmglst2*(-1931 + 1098*z2 + 180*z3) + Mst2*(-22430251 +
        5488000*OepS2 - 324135000*S2 - 8232000*T1ep - 103929000*z2 + 28812000*
        z3 - 4116000*z4 - 6174000*pow2(z2))) + 333396000*s2t*(Dmglst2*(-465 +
        348*z2 - 176*z3) + Mst2*(35 - 20*z2 + 16*z3))*pow3(Mt) + 333396000*Mt*(
        Dmglst2*(20 - 11*z2) + Mst2*(-2 + z2))*pow2(Mst2)*pow3(s2t) + 3456*(-
        6582784 + 2701125*z2 - 2315250*z3)*pow4(Mt) + (4630500*Dmglst2*(38 +
        63*z2 - 90*z3) + Mst2*(378483467 + 9604000*OepS2 - 754771500*S2 -
        14406000*T1ep - 306899250*z2 + 133770000*z3 - 7203000*z4 - 10804500*
        pow2(z2)))*pow3(Mst2)*pow4(s2t) - 630*log(pow2(Mst2)/pow2(Mst1))*(2*
        Mst2*(8*Mst2*(58853 - 44100*S2) - 11025*Dmglst2*(-37 + 60*z2))*pow2(Mt)
        *pow2(s2t) + 1058400*s2t*(Mst2*(25 - 8*z2) + 11*Dmglst2*(-17 + 8*z2))*
        pow3(Mt) + 264600*(3*Dmglst2 - Mst2)*Mt*pow2(Mst2)*pow3(s2t) + 288*(-
        175121 + 44100*z2)*pow4(Mt) + (-11025*Dmglst2*(-83 + 60*z2) + Mst2*(-
        621671 + 308700*S2 + 132300*z2))*pow3(Mst2)*pow4(s2t)) - 1190700*pow2(
        log(pow2(Mst2)/pow2(Mst1)))*(74*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 1120*(
        3*Dmglst2 - Mst2)*s2t*pow3(Mt) + 2736*pow4(Mt) - 53*pow4(Mst2)*pow4(
        s2t)))*pow6(Mst1))))/pow2(Mst1) - (Mst2*(-4*Dmglst2*(-66150*pow4(Mst1)*
        pow4(Mst2)*(-672*(-11674 + 120*B4 - 120*D3 + 60*DN + 17091*S2 + 4849*z2
        + 2195*z3 - 540*z4)*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 140*Mst2*s2t*(-
        30183 + 1120*OepS2 + 43740*S2 - 1680*T1ep - 75402*z2 + 216*z3 + 6072*z4
        - 1260*pow2(z2))*pow3(Mt) + 7*Mt*(69997 + 188640*B4 - 3600*DN + 5600*
        OepS2 + 146772*S2 - 8400*T1ep - 1043682*z2 + 1816440*z3 + 32520*z4 -
        317340*pow2(z2))*pow3(Mst2)*pow3(s2t) + 16*(1547066 + 7840*OepS2 +
        258066*S2 - 11760*T1ep - 1256046*z2 - 946785*z3 - 5880*z4 - 8820*pow2(
        z2))*pow4(Mt) - 35*(-15707 + 432*B4 - 576*D3 + 360*DN + 224*OepS2 +
        543348*S2 - 336*T1ep - 11802*z2 - 105690*z3 - 2544*z4 - 2844*pow2(z2))*
        pow4(Mst2)*pow4(s2t) - 1680*pow2(log(pow2(Mst2)/pow2(Mst1)))*(-1806*
        pow2(Mst2)*pow2(Mt)*pow2(s2t) - 1592*Mst2*s2t*pow3(Mt) - 1104*Mt*pow3(
        Mst2)*pow3(s2t) + 3808*pow4(Mt) - 213*pow4(Mst2)*pow4(s2t)) + 2520*
        pow3(log(pow2(Mst2)/pow2(Mst1)))*(-256*pow2(Mst2)*pow2(Mt)*pow2(s2t) +
        784*Mst2*s2t*pow3(Mt) - 516*Mt*pow3(Mst2)*pow3(s2t) + 512*pow4(Mt) -
        107*pow4(Mst2)*pow4(s2t)) - 280*log(pow2(Mst2)/pow2(Mst1))*(12*(3427 -
        1782*z2)*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 4*Mst2*s2t*(-8398 + 2835*S2 -
        3258*z2)*pow3(Mt) + 9*Mt*(1640 + 315*S2 - 1146*z2)*pow3(Mst2)*pow3(s2t)
        + 4*(-6569 + 2268*S2 - 7056*z2)*pow4(Mt) - 9*(277 + 63*S2 - 140*z2)*
        pow4(Mst2)*pow4(s2t))) + 1764*pow2(Mst2)*(-112*pow2(Mst2)*pow2(Mt)*
        pow2(s2t)*(1216808 + 150000*OepS2 - 376650*S2 - 225000*T1ep - 3408750*
        z2 - 833625*z3 - 112500*z4 - 168750*pow2(z2)) + 1050*Mst2*s2t*(-322421
        + 73760*OepS2 - 564876*S2 - 110640*T1ep - 3060894*z2 + 1606920*z3 -
        55320*z4 - 82980*pow2(z2))*pow3(Mt) - 175*Mt*(27211 + 36720*B4 + 1080*
        DN + 64160*OepS2 - 1515564*S2 - 96240*T1ep - 1240446*z2 + 245400*z3 -
        12480*z4 - 72180*pow2(z2))*pow3(Mst2)*pow3(s2t) - 8*(98884301 +
        4508000*OepS2 - 17853750*S2 - 6762000*T1ep - 477659550*z2 + 302883000*
        z3 - 3381000*z4 - 5071500*pow2(z2))*pow4(Mt) - 7*(-95041 + 81000*B4 -
        108000*D3 + 67500*DN - 432000*OepS2 - 7425000*S2 + 648000*T1ep +
        5325300*z2 + 4227750*z3 - 121500*z4)*pow4(Mst2)*pow4(s2t) - 63000*pow3(
        log(pow2(Mst2)/pow2(Mst1)))*(-1156*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 384*
        Mst2*s2t*pow3(Mt) - 390*Mt*pow3(Mst2)*pow3(s2t) + 2632*pow4(Mt) - 29*
        pow4(Mst2)*pow4(s2t)) + 3150*pow2(log(pow2(Mst2)/pow2(Mst1)))*(-12992*
        pow2(Mst2)*pow2(Mt)*pow2(s2t) - 175520*Mst2*s2t*pow3(Mt) + 3600*Mt*
        pow3(Mst2)*pow3(s2t) + 285744*pow4(Mt) + 4969*pow4(Mst2)*pow4(s2t)) +
        210*log(pow2(Mst2)/pow2(Mst1))*(8*(77657 + 202500*S2 - 55800*z2)*pow2(
        Mst2)*pow2(Mt)*pow2(s2t) - 200*Mst2*s2t*(-11408 + 37341*S2 - 17730*z2)*
        pow3(Mt) + 300*Mt*(-503 + 3609*S2 - 1494*z2)*pow3(Mst2)*pow3(s2t) + 48*
        (8581 + 72450*S2 - 137700*z2)*pow4(Mt) + (81643 - 291600*S2 - 5850*z2)*
        pow4(Mst2)*pow4(s2t)))*pow6(Mst1) - 83349000*pow2(Mst1)*(-4*(2*Mt +
        Mst2*s2t)*(14*Mt*(3 + z2) + Mst2*s2t*(5 + 7*z2))*pow2(-2*Mt + Mst2*s2t)
        - 2*log(pow2(Mst2)/pow2(Mst1))*(-536*pow2(Mst2)*pow2(Mt)*pow2(s2t) -
        512*Mst2*s2t*pow3(Mt) + 128*Mt*pow3(Mst2)*pow3(s2t) + 560*pow4(Mt) +
        131*pow4(Mst2)*pow4(s2t)) + pow2(log(pow2(Mst2)/pow2(Mst1)))*(-728*
        pow2(Mst2)*pow2(Mt)*pow2(s2t) - 1536*Mst2*s2t*pow3(Mt) + 384*Mt*pow3(
        Mst2)*pow3(s2t) + 944*pow4(Mt) + 187*pow4(Mst2)*pow4(s2t)))*pow6(Mst2)
        - (72*pow2(Mst2)*pow2(Mt)*pow2(s2t)*(-3650009897 + 2110136000*OepS2 -
        52385772600*S2 - 3165204000*T1ep - 40405228500*z2 + 755286000*z3 -
        1582602000*z4 - 2373903000*pow2(z2)) - 274400*Mst2*s2t*(-20964397 +
        2058400*OepS2 - 47545272*S2 - 3087600*T1ep - 59449002*z2 + 39214920*z3
        - 1543800*z4 - 2315700*pow2(z2))*pow3(Mt) + 17150*Mt*(-31897243 +
        2491360*OepS2 - 90290268*S2 - 3737040*T1ep - 51499698*z2 + 33912840*z3
        - 1868520*z4 - 2802780*pow2(z2))*pow3(Mst2)*pow3(s2t) + 80*(-
        54946675289 + 3063401600*OepS2 - 52952863320*S2 - 4595102400*T1ep -
        225236187120*z2 + 243414477600*z3 - 2297551200*z4 - 3446326800*pow2(z2)
        )*pow4(Mt) + (119405394763 - 11522056000*OepS2 + 478191735000*S2 +
        17283084000*T1ep + 233890773900*z2 - 111792103500*z3 + 8641542000*z4 +
        12962313000*pow2(z2))*pow4(Mst2)*pow4(s2t) + 55566000*pow3(log(pow2(
        Mst2)/pow2(Mst1)))*(-2896*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 2016*Mst2*
        s2t*pow3(Mt) - 672*Mt*pow3(Mst2)*pow3(s2t) + 7424*pow4(Mt) - 33*pow4(
        Mst2)*pow4(s2t)) - 793800*pow2(log(pow2(Mst2)/pow2(Mst1)))*(-420096*
        pow2(Mst2)*pow2(Mt)*pow2(s2t) - 1672160*Mst2*s2t*pow3(Mt) + 68880*Mt*
        pow3(Mst2)*pow3(s2t) + 3253120*pow4(Mt) + 1377*pow4(Mst2)*pow4(s2t)) -
        1260*log(pow2(Mst2)/pow2(Mst1))*(72*(-4261807 + 33912900*S2 - 73500*z2)
        *pow2(Mst2)*pow2(Mt)*pow2(s2t) - 78400*Mst2*s2t*(-4724 + 115785*S2 -
        29484*z2)*pow3(Mt) + 4900*Mt*(18652 + 140139*S2 - 38826*z2)*pow3(Mst2)*
        pow3(s2t) + 80*(42300121 + 49233240*S2 - 64139040*z2)*pow4(Mt) + (
        8287903 - 185175900*S2 + 9525600*z2)*pow4(Mst2)*pow4(s2t)))*pow8(Mst1)
        - 21337344000*(-1 + log(pow2(Mst2)/pow2(Mst1)))*log(pow2(Mst2)/pow2(
        Mst1))*pow2(Mt)*(-2*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow8(Mst2)) + 3*
        Mst2*(-308700*pow4(Mst1)*pow4(Mst2)*(-2*pow2(Mst2)*pow2(Mt)*pow2(s2t)*(
        1013263 + 2880*D3 - 2880*DN - 25440*OepS2 - 137052*S2 + 38160*T1ep +
        35562*z2 + 83160*z3 + 36360*z4 + 28620*pow2(z2)) - 24*Mst2*s2t*(-60759
        + 1120*OepS2 + 25596*S2 - 1680*T1ep - 137386*z2 + 115320*z3 - 12360*z4
        - 1260*pow2(z2))*pow3(Mt) - 30*Mt*(25611 + 5280*B4 - 48*DN - 224*OepS2
        - 324*S2 + 336*T1ep - 13614*z2 + 47720*z3 + 2040*z4 - 6660*pow2(z2))*
        pow3(Mst2)*pow3(s2t) - 8*(389597 + 3360*OepS2 - 105948*S2 - 5040*T1ep +
        293898*z2 - 740160*z3 - 2520*z4 - 3780*pow2(z2))*pow4(Mt) + 5*(25289 +
        1440*B4 - 144*D3 + 72*DN - 2208*OepS2 + 298404*S2 + 3312*T1ep + 36318*
        z2 - 94968*z3 + 1440*z4 - 108*pow2(z2))*pow4(Mst2)*pow4(s2t) + 180*
        pow3(log(pow2(Mst2)/pow2(Mst1)))*(56*pow2(Mst2)*pow2(Mt)*pow2(s2t) -
        3136*Mst2*s2t*pow3(Mt) + 1040*Mt*pow3(Mst2)*pow3(s2t) + 32*pow4(Mt) +
        271*pow4(Mst2)*pow4(s2t)) - 120*pow2(log(pow2(Mst2)/pow2(Mst1)))*(3192*
        pow2(Mst2)*pow2(Mt)*pow2(s2t) - 14880*Mst2*s2t*pow3(Mt) + 5664*Mt*pow3(
        Mst2)*pow3(s2t) + 7744*pow4(Mt) + 1113*pow4(Mst2)*pow4(s2t)) + 40*log(
        pow2(Mst2)/pow2(Mst1))*(-6*(-4534 + 4293*S2 + 1830*z2 - 72*z3)*pow2(
        Mst2)*pow2(Mt)*pow2(s2t) + 72*Mst2*s2t*(26 + 189*S2 - 758*z2)*pow3(Mt)
        - 126*Mt*(-208 + 27*S2 - 2*z2)*pow3(Mst2)*pow3(s2t) + 8*(-553 + 1701*S2
        + 4023*z2 + 108*z3)*pow4(Mt) + 3*(1148 + 1863*S2 + 276*z2 - 18*z3)*
        pow4(Mst2)*pow4(s2t))) - 20580*pow2(Mst2)*(4*pow2(Mst2)*pow2(Mt)*pow2(
        s2t)*(-2367149 + 21600*D3 - 21600*DN + 503200*OepS2 - 14889420*S2 -
        754800*T1ep - 9292470*z2 + 66000*z3 - 247800*z4 - 566100*pow2(z2)) -
        40*Mst2*s2t*(-1156193 + 60320*OepS2 - 1414908*S2 - 90480*T1ep -
        2823942*z2 + 2907240*z3 - 45240*z4 - 67860*pow2(z2))*pow3(Mt) - 100*Mt*
        (51041 + 7344*B4 + 216*DN - 2336*OepS2 + 112428*S2 + 3504*T1ep + 28038*
        z2 - 78216*z3 + 8880*z4 + 2628*pow2(z2))*pow3(Mst2)*pow3(s2t) - 8*(
        7319011 + 167200*OepS2 - 6744060*S2 - 250800*T1ep + 9301170*z2 -
        20715000*z3 - 125400*z4 - 188100*pow2(z2))*pow4(Mt) + (7680127 - 10800*
        B4 + 21600*D3 - 16200*DN - 514400*OepS2 + 46307700*S2 + 771600*T1ep +
        10868970*z2 - 7186200*z3 + 466800*z4 + 773100*pow2(z2))*pow4(Mst2)*
        pow4(s2t) - 3600*pow3(log(pow2(Mst2)/pow2(Mst1)))*(-1166*pow2(Mst2)*
        pow2(Mt)*pow2(s2t) + 2304*Mst2*s2t*pow3(Mt) - 780*Mt*pow3(Mst2)*pow3(
        s2t) + 440*pow4(Mt) - 115*pow4(Mst2)*pow4(s2t)) - 14400*pow2(log(pow2(
        Mst2)/pow2(Mst1)))*(479*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 2448*Mst2*s2t*
        pow3(Mt) + 366*Mt*pow3(Mst2)*pow3(s2t) + 1136*pow4(Mt) + 23*pow4(Mst2)*
        pow4(s2t)) + 120*log(pow2(Mst2)/pow2(Mst1))*(-4*(-52322 + 84915*S2 +
        3060*z2 + 540*z3)*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 1080*Mst2*s2t*(42 +
        377*S2 - 330*z2)*pow3(Mt) - 180*Mt*(221 + 219*S2 - 234*z2)*pow3(Mst2)*
        pow3(s2t) + 8*(-22352 + 28215*S2 + 36270*z2)*pow4(Mt) + (-16051 +
        86805*S2 - 8325*z2)*pow4(Mst2)*pow4(s2t)))*pow6(Mst1) - 55566000*pow2(
        Mst1)*(8*pow2(Mst2)*pow2(Mt)*pow2(s2t)*(87 + 38*z2 - 154*log(pow2(Mst2)
        /pow2(Mst1)) + 75*pow2(log(pow2(Mst2)/pow2(Mst1)))) - 256*pow2(-1 +
        log(pow2(Mst2)/pow2(Mst1)))*(-4*Mst2*s2t*pow3(Mt) + Mt*pow3(Mst2)*pow3(
        s2t)) - 16*(71 + 38*z2 - 122*log(pow2(Mst2)/pow2(Mst1)) + 59*pow2(log(
        pow2(Mst2)/pow2(Mst1))))*pow4(Mt) - (135 + 38*z2 - 250*log(pow2(Mst2)/
        pow2(Mst1)) + 123*pow2(log(pow2(Mst2)/pow2(Mst1))))*pow4(Mst2)*pow4(
        s2t))*pow6(Mst2) + 196*pow2(Dmglst2)*(283500*pow2(Mst1)*pow4(Mst2)*(-8*
        pow2(Mst2)*pow2(Mt)*pow2(s2t)*(-716 + 104*z2 + 102*log(pow2(Mst2)/pow2(
        Mst1)) + 155*pow2(log(pow2(Mst2)/pow2(Mst1)))) + 256*(-3 + 4*log(pow2(
        Mst2)/pow2(Mst1)) + 3*pow2(log(pow2(Mst2)/pow2(Mst1))))*(-4*Mst2*s2t*
        pow3(Mt) + Mt*pow3(Mst2)*pow3(s2t)) + 16*(-620 + 104*z2 + 166*log(pow2(
        Mst2)/pow2(Mst1)) + 59*pow2(log(pow2(Mst2)/pow2(Mst1))))*pow4(Mt) + (-
        940 + 104*z2 + 102*log(pow2(Mst2)/pow2(Mst1)) + 443*pow2(log(pow2(Mst2)
        /pow2(Mst1))))*pow4(Mst2)*pow4(s2t)) + 15*pow2(Mst2)*pow4(Mst1)*(-360*
        pow2(Mst2)*pow2(Mt)*pow2(s2t)*(-2164753 + 3360*B4 - 3360*D3 + 1680*DN -
        7840*OepS2 + 169452*S2 + 11760*T1ep + 467818*z2 + 1945895*z3 - 9240*z4
        + 8820*pow2(z2)) + 320*Mst2*s2t*(-1243547 + 15680*OepS2 + 953208*S2 -
        23520*T1ep - 2985438*z2 + 3393915*z3 - 11760*z4 - 17640*pow2(z2))*pow3(
        Mt) + 560*Mt*(-282298 + 61560*B4 - 1620*DN - 2240*OepS2 - 111213*S2 +
        3360*T1ep - 281922*z2 + 983190*z3 + 6540*z4 - 114120*pow2(z2))*pow3(
        Mst2)*pow3(s2t) - 16*(-46632377 + 744800*OepS2 + 29679210*S2 - 1117200*
        T1ep - 107181930*z2 + 108350025*z3 - 558600*z4 - 837900*pow2(z2))*pow4(
        Mt) - 35*(-924689 + 41040*B4 - 43200*D3 + 22680*DN - 1120*OepS2 +
        21325356*S2 + 1680*T1ep - 127686*z2 - 3687510*z3 - 190320*z4 - 37620*
        pow2(z2))*pow4(Mst2)*pow4(s2t) + 2520*pow2(log(pow2(Mst2)/pow2(Mst1)))*
        (-7740*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 61760*Mst2*s2t*pow3(Mt) - 18480*
        Mt*pow3(Mst2)*pow3(s2t) + 66208*pow4(Mt) - 4455*pow4(Mst2)*pow4(s2t)) +
        37800*pow3(log(pow2(Mst2)/pow2(Mst1)))*(-256*pow2(Mst2)*pow2(Mt)*pow2(
        s2t) - 768*Mt*pow3(Mst2)*pow3(s2t) + 512*pow4(Mt) - 107*pow4(Mst2)*
        pow4(s2t)) + 168*log(pow2(Mst2)/pow2(Mst1))*(-300*(6485 + 1134*S2 -
        6306*z2)*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 1600*Mst2*s2t*(2182 + 378*S2 -
        1665*z2)*pow3(Mt) + 600*Mt*(439 + 252*S2 + 135*z2)*pow3(Mst2)*pow3(s2t)
        + 16*(-56837 + 89775*S2 - 162900*z2)*pow4(Mt) - 75*(-1388 + 63*S2 -
        618*z2)*pow4(Mst2)*pow4(s2t))) - 2*(-36*pow2(Mst2)*pow2(Mt)*pow2(s2t)*(
        68526791 + 2772000*OepS2 - 107403300*S2 - 4158000*T1ep - 140027250*z2 +
        39514125*z3 - 2079000*z4 - 3118500*pow2(z2)) - 800*Mst2*s2t*(7384247 +
        1183840*OepS2 + 5821092*S2 - 1775760*T1ep - 59988552*z2 + 25757760*z3 -
        887880*z4 - 1331820*pow2(z2))*pow3(Mt) + 350*Mt*(-2727217 + 437920*
        OepS2 - 3648132*S2 - 656880*T1ep - 12301638*z2 + 4166400*z3 - 328440*z4
        - 492660*pow2(z2))*pow3(Mst2)*pow3(s2t) + 8*(648916519 + 158732000*
        OepS2 + 66753450*S2 - 238098000*T1ep - 10760276550*z2 + 6504855000*z3 -
        119049000*z4 - 178573500*pow2(z2))*pow4(Mt) - 7*(-21452821 + 243000*B4
        - 324000*D3 + 202500*DN + 2272000*OepS2 - 75724200*S2 - 3408000*T1ep -
        55081500*z2 + 67149750*z3 - 3040500*z4 - 4014000*pow2(z2))*pow4(Mst2)*
        pow4(s2t) - 9450*pow2(log(pow2(Mst2)/pow2(Mst1)))*(-133776*pow2(Mst2)*
        pow2(Mt)*pow2(s2t) - 398240*Mst2*s2t*pow3(Mt) - 24960*Mt*pow3(Mst2)*
        pow3(s2t) + 1267120*pow4(Mt) - 4293*pow4(Mst2)*pow4(s2t)) + 567000*
        pow3(log(pow2(Mst2)/pow2(Mst1)))*(380*pow2(Mst2)*pow2(Mt)*pow2(s2t) -
        512*Mst2*s2t*pow3(Mt) + 1608*pow4(Mt) + 13*pow4(Mst2)*pow4(s2t)) - 630*
        log(pow2(Mst2)/pow2(Mst1))*(-24*(102713 + 133650*S2 - 173700*z2)*pow2(
        Mst2)*pow2(Mt)*pow2(s2t) - 400*Mst2*s2t*(11587 + 76104*S2 - 49608*z2)*
        pow3(Mt) + 200*Mt*(-4958 + 24633*S2 - 90*z2)*pow3(Mst2)*pow3(s2t) + 16*
        (-245893 + 2551050*S2 - 2175300*z2)*pow4(Mt) + (-332311 - 511200*S2 +
        49050*z2)*pow4(Mst2)*pow4(s2t)))*pow6(Mst1) - 9072000*(-2 + log(pow2(
        Mst2)/pow2(Mst1)) + 3*pow2(log(pow2(Mst2)/pow2(Mst1))))*pow2(-4*pow2(
        Mt) + pow2(Mst2)*pow2(s2t))*pow6(Mst2)) + (-16*pow2(Mst2)*pow2(Mt)*
        pow2(s2t)*(-25002374272 + 4703902000*OepS2 - 184815757350*S2 -
        7055853000*T1ep - 84661721025*z2 + 16805113500*z3 - 3527926500*z4 -
        5291889750*pow2(z2)) + 1097600*Mst2*s2t*(-2016907 + 109600*OepS2 -
        3520152*S2 - 164400*T1ep - 3448554*z2 + 3465480*z3 - 82200*z4 - 123300*
        pow2(z2))*pow3(Mt) - 68600*Mt*(-2759935 + 70240*OepS2 - 3305340*S2 -
        105360*T1ep - 629034*z2 + 2677800*z3 - 52680*z4 - 79020*pow2(z2))*pow3(
        Mst2)*pow3(s2t) + 16*(187545955063 + 3204992000*OepS2 - 128216692800*S2
        - 4807488000*T1ep + 129009640200*z2 - 399103824000*z3 - 2403744000*z4 -
        3605616000*pow2(z2))*pow4(Mt) + (93508520089 + 2000376000*B4 +
        222264000*D3 - 222264000*DN + 4925480000*OepS2 - 312095700000*S2 -
        7388220000*T1ep - 77106571500*z2 + 21506100000*z3 - 2360526000*z4 -
        5541165000*pow2(z2))*pow4(Mst2)*pow4(s2t) + 18522000*pow3(log(pow2(
        Mst2)/pow2(Mst1)))*(-2696*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 7808*Mst2*
        s2t*pow3(Mt) - 3200*Mt*pow3(Mst2)*pow3(s2t) + 6400*pow4(Mt) - 409*pow4(
        Mst2)*pow4(s2t)) + 2116800*pow2(log(pow2(Mst2)/pow2(Mst1)))*(73657*
        pow2(Mst2)*pow2(Mt)*pow2(s2t) - 316400*Mst2*s2t*pow3(Mt) + 35000*Mt*
        pow3(Mst2)*pow3(s2t) + 68678*pow4(Mt) + 11764*pow4(Mst2)*pow4(s2t)) -
        1260*log(pow2(Mst2)/pow2(Mst1))*(-24*(-19407863 + 50398950*S2 +
        2866500*z2)*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 940800*Mst2*s2t*(167 +
        2055*S2 - 912*z2)*pow3(Mt) - 58800*Mt*(46 + 1317*S2 - 918*z2)*pow3(
        Mst2)*pow3(s2t) + 32*(-35585111 + 25754400*S2 + 32281200*z2)*pow4(Mt) +
        3*(-23468297 + 26386500*S2 + 661500*z2 + 176400*z3)*pow4(Mst2)*pow4(
        s2t)))*pow8(Mst1) + 7112448000*pow2(Mt)*(-2*pow2(Mt) + pow2(Mst2)*pow2(
        s2t))*pow2(-1 + log(pow2(Mst2)/pow2(Mst1)))*pow8(Mst2))))/pow4(Mst1)))/
        (6.001128e9*pow6(Mst2)))/pow4(Mt)/pow2(Sbeta)*12.;

   return result;
}

/**
 * 	@return returns the susy log^1 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double himalaya::H6b::calc_coef_at_as2_no_sm_logs_log1() const {

   const double result =
      ((pow2(Sbeta)*(6000*pow2(Dmsqst2)*(72*(2*Mst2*s2t*(3 - 2*z2) + Mt*(-7 + 5*
        z2))*pow2(Mst1) + (36*Mst2*s2t*(7 - 4*z2) + Mt*(-161 + 108*z2))*pow2(
        Mst2))*pow3(Mt)*pow4(Mst1) + 3000*Dmsqst2*pow2(Mst1)*(9*Mst2*pow4(Mst1)
        *(-32*Dmglst2*(2*Mst2*s2t*(12 - 7*z2) + 3*Mt*(-11 + 6*z2))*pow3(Mt) +
        Mst2*(8*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 64*Mst2*s2t*(-2 + z2)*pow3(Mt)
        + 16*(-9 + 4*z2)*pow4(Mt) + pow4(Mst2)*pow4(s2t))) + pow2(Mst1)*pow3(
        Mst2)*(-8*Dmglst2*(12*Mst2*s2t*(29 - 18*z2) + Mt*(-283 + 180*z2))*pow3(
        Mt) + 3*Mst2*(-48*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 96*Mst2*s2t*(-3 + 2*
        z2)*pow3(Mt) + 32*(-8 + 3*z2)*pow4(Mt) + 3*pow4(Mst2)*pow4(s2t))) + 9*(
        16*s2t*(-83*Dmglst2 + 9*Mst2 + 44*Dmglst2*z2 - 4*Mst2*z2)*pow3(Mt) + 8*
        (-31 + 12*z2)*pow4(Mt) - pow4(Mst2)*pow4(s2t))*pow6(Mst1) - 9*pow2(-4*
        pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow6(Mst2)) + 450*Mst2*pow2(log(pow2(
        Mst2)/pow2(Mst1)))*pow4(Mst1)*(-8*Dmglst2*(8*Mt*(-157*Mst2*s2t*pow2(Mt)
        - 108*Mt*pow2(Mst2)*pow2(s2t) + 398*pow3(Mt) - 16*pow3(Mst2)*pow3(s2t))
        *pow4(Mst1) + pow4(Mst2)*(-320*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 848*
        Mst2*s2t*pow3(Mt) - 532*Mt*pow3(Mst2)*pow3(s2t) + 512*pow4(Mt) - 99*
        pow4(Mst2)*pow4(s2t)) + pow2(Mst1)*pow2(Mst2)*(-768*pow2(Mst2)*pow2(Mt)
        *pow2(s2t) - 448*Mst2*s2t*pow3(Mt) - 212*Mt*pow3(Mst2)*pow3(s2t) +
        1896*pow4(Mt) + 35*pow4(Mst2)*pow4(s2t))) + Mst2*(8*pow2(Mst1)*pow2(
        Mst2)*(425*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 704*Mst2*s2t*pow3(Mt) + 212*
        Mt*pow3(Mst2)*pow3(s2t) - 76*pow4(Mt) + 3*pow4(Mst2)*pow4(s2t)) + pow4(
        Mst1)*(2688*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 5056*Mst2*s2t*pow3(Mt) +
        1024*Mt*pow3(Mst2)*pow3(s2t) - 2848*pow4(Mt) + 41*pow4(Mst2)*pow4(s2t))
        + pow4(Mst2)*(296*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 6784*Mst2*s2t*pow3(
        Mt) + 2208*Mt*pow3(Mst2)*pow3(s2t) - 672*pow4(Mt) + 519*pow4(Mst2)*
        pow4(s2t))) + 4*Mst2*pow2(Dmglst2)*(320*pow2(Mt)*pow2(s2t)*pow4(Mst2) -
        512*pow2(Mst2)*pow4(Mt) + pow2(Mst1)*(768*pow2(Mst2)*pow2(Mt)*pow2(s2t)
        - 1152*Mst2*s2t*pow3(Mt) + 3640*pow4(Mt) - 35*pow4(Mst2)*pow4(s2t)) +
        768*Mt*pow3(s2t)*pow5(Mst2) + 99*pow4(s2t)*pow6(Mst2))) + Mst2*(4*Mst2*
        pow2(Dmglst2)*(-225*pow2(Mst1)*pow4(Mst2)*(-664*pow2(Mst2)*pow2(Mt)*
        pow2(s2t) - 2048*Mst2*s2t*pow3(Mt) + 512*Mt*pow3(Mst2)*pow3(s2t) +
        1840*pow4(Mt) + 83*pow4(Mst2)*pow4(s2t)) + 2*pow2(Mst2)*pow4(Mst1)*(-
        150*(-7121 + 6336*z2)*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 400*Mst2*s2t*(-
        1934 + 1791*z2)*pow3(Mt) - 300*Mt*(-875 + 648*z2)*pow3(Mst2)*pow3(s2t)
        + 8*(71687 + 76050*z2)*pow4(Mt) + 75*(-205 + 129*z2)*pow4(Mst2)*pow4(
        s2t)) + 25*(12*(7903 - 5760*z2)*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 96*
        Mst2*s2t*(-1055 + 2181*z2)*pow3(Mt) + 13992*Mt*pow3(Mst2)*pow3(s2t) +
        80*(-605 + 4608*z2)*pow4(Mt) + 3*(409 - 258*z2)*pow4(Mst2)*pow4(s2t))*
        pow6(Mst1) + 3600*pow2(-4*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow6(Mst2))
        - 100*Dmglst2*(2*pow4(Mst1)*pow4(Mst2)*(12*(-3565 + 1728*z2)*pow2(Mst2)
        *pow2(Mt)*pow2(s2t) + 8*Mst2*s2t*(-3659 + 252*z2)*pow3(Mt) + 72*Mt*(-
        192 + 91*z2)*pow3(Mst2)*pow3(s2t) + 44*(427 + 432*z2)*pow4(Mt) + 9*(211
        - 86*z2)*pow4(Mst2)*pow4(s2t)) + 6*pow2(Mst2)*(4*(-1955 + 1152*z2)*
        pow2(Mst2)*pow2(Mt)*pow2(s2t) - 48*Mst2*s2t*(-13 + 262*z2)*pow3(Mt) +
        24*Mt*(111 - 17*z2)*pow3(Mst2)*pow3(s2t) + 8*(-1471 + 3240*z2)*pow4(Mt)
        + 3*(-1 + 86*z2)*pow4(Mst2)*pow4(s2t))*pow6(Mst1) - 18*pow2(Mst1)*(-
        536*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 512*Mst2*s2t*pow3(Mt) + 128*Mt*
        pow3(Mst2)*pow3(s2t) + 560*pow4(Mt) + 131*pow4(Mst2)*pow4(s2t))*pow6(
        Mst2) + 3*(192*(-77 + 48*z2)*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 48*Mst2*
        s2t*(-867 + 1012*z2)*pow3(Mt) - 276*Mt*pow3(Mst2)*pow3(s2t) + 16*(-6445
        + 6768*z2)*pow4(Mt) - 441*pow4(Mst2)*pow4(s2t))*pow8(Mst1) - 2304*pow2(
        Mt)*(-2*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow8(Mst2)) + 25*Mst2*(4*pow4(
        Mst1)*pow4(Mst2)*(-72*(-445 + 96*z2)*pow2(Mst2)*pow2(Mt)*pow2(s2t) -
        144*Mst2*s2t*(-121 + 296*z2)*pow3(Mt) - 144*Mt*(-151 + 19*z2)*pow3(
        Mst2)*pow3(s2t) + 8*(-1993 + 2151*z2 + 216*z3)*pow4(Mt) + 81*(8 + 5*z2)
        *pow4(Mst2)*pow4(s2t)) - 36*pow2(Mst2)*(-696*pow2(Mst2)*pow2(Mt)*pow2(
        s2t) + 32*Mst2*s2t*(-179 + 128*z2)*pow3(Mt) + 16*Mt*(148 - 17*z2)*pow3(
        Mst2)*pow3(s2t) - 8*(-403 + 300*z2)*pow4(Mt) + (551 + 86*z2)*pow4(Mst2)
        *pow4(s2t))*pow6(Mst1) + 36*pow2(Mst1)*(-744*pow2(Mst2)*pow2(Mt)*pow2(
        s2t) - 1024*Mst2*s2t*pow3(Mt) + 256*Mt*pow3(Mst2)*pow3(s2t) + 1232*
        pow4(Mt) + 141*pow4(Mst2)*pow4(s2t))*pow6(Mst2) + 9*(864*pow2(Mst2)*
        pow2(Mt)*pow2(s2t) - 512*Mst2*s2t*(-49 + 31*z2)*pow3(Mt) - 784*Mt*pow3(
        Mst2)*pow3(s2t) + 8*(-2503 + 1408*z2)*pow4(Mt) + (1547 + 164*z2)*pow4(
        Mst2)*pow4(s2t))*pow8(Mst1) + 4608*pow2(Mt)*(-2*pow2(Mt) + pow2(Mst2)*
        pow2(s2t))*pow8(Mst2))) - 30*log(pow2(Mst2)/pow2(Mst1))*(7200*pow2(
        Dmsqst2)*(pow2(Mst1) - pow2(Mst2))*pow4(Mst1)*pow4(Mt) - 1800*Dmsqst2*
        pow4(Mst1)*(8*(5*Mt + 6*Dmglst2*s2t - 2*Mst2*s2t)*pow3(Mt)*pow4(Mst1) -
        8*pow4(Mst2)*pow4(Mt) - Mst2*pow2(Mst1)*(16*Dmglst2*pow4(Mt) - 8*Mst2*
        pow4(Mt) + pow4(s2t)*pow5(Mst2)) + pow4(s2t)*pow8(Mst2)) + Mst2*(2*
        Mst2*pow2(Dmglst2)*(2*pow2(Mst2)*pow4(Mst1)*(3240*pow2(Mst2)*pow2(Mt)*
        pow2(s2t) + 45920*Mst2*s2t*pow3(Mt) - 9720*Mt*pow3(Mst2)*pow3(s2t) +
        37408*pow4(Mt) - 4635*pow4(Mst2)*pow4(s2t)) + 45*pow2(Mst1)*pow4(Mst2)*
        (-456*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 1024*Mst2*s2t*pow3(Mt) + 256*Mt*
        pow3(Mst2)*pow3(s2t) + 400*pow4(Mt) + 153*pow4(Mst2)*pow4(s2t)) + 2*(-
        34080*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 147560*Mst2*s2t*pow3(Mt) - 5880*
        Mt*pow3(Mst2)*pow3(s2t) + 415932*pow4(Mt) + 1005*pow4(Mst2)*pow4(s2t))*
        pow6(Mst1) - 1440*pow2(-4*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow6(Mst2))
        - 20*Dmglst2*(-4*pow4(Mst1)*pow4(Mst2)*(1860*pow2(Mst2)*pow2(Mt)*pow2(
        s2t) + 1472*Mst2*s2t*pow3(Mt) + 774*Mt*pow3(Mst2)*pow3(s2t) - 2464*
        pow4(Mt) + 15*pow4(Mst2)*pow4(s2t)) + 2*pow2(Mst2)*(-432*pow2(Mst2)*
        pow2(Mt)*pow2(s2t) - 11656*Mst2*s2t*pow3(Mt) + 996*Mt*pow3(Mst2)*pow3(
        s2t) + 18748*pow4(Mt) + 207*pow4(Mst2)*pow4(s2t))*pow6(Mst1) - 3*pow2(
        Mst1)*(-856*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 1536*Mst2*s2t*pow3(Mt) +
        384*Mt*pow3(Mst2)*pow3(s2t) + 1200*pow4(Mt) + 203*pow4(Mst2)*pow4(s2t))
        *pow6(Mst2) + 4*(-738*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 10004*Mst2*s2t*
        pow3(Mt) + 135*Mt*pow3(Mst2)*pow3(s2t) + 17532*pow4(Mt) + 12*pow4(Mst2)
        *pow4(s2t))*pow8(Mst1) - 768*pow2(Mt)*(-2*pow2(Mt) + pow2(Mst2)*pow2(
        s2t))*pow8(Mst2)) + 5*Mst2*(4*pow4(Mst1)*pow4(Mst2)*(3372*pow2(Mst2)*
        pow2(Mt)*pow2(s2t) - 9024*Mst2*s2t*pow3(Mt) + 3912*Mt*pow3(Mst2)*pow3(
        s2t) + 6112*pow4(Mt) + 579*pow4(Mst2)*pow4(s2t)) - 24*pow2(Mst2)*(-534*
        pow2(Mst2)*pow2(Mt)*pow2(s2t) + 1824*Mst2*s2t*pow3(Mt) - 60*Mt*pow3(
        Mst2)*pow3(s2t) - 934*pow4(Mt) + 97*pow4(Mst2)*pow4(s2t))*pow6(Mst1) +
        6*pow2(Mst1)*(-728*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 1024*Mst2*s2t*pow3(
        Mt) + 256*Mt*pow3(Mst2)*pow3(s2t) + 1200*pow4(Mt) + 139*pow4(Mst2)*
        pow4(s2t))*pow6(Mst2) + 3*(4640*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 15296*
        Mst2*s2t*pow3(Mt) + 240*Mt*pow3(Mst2)*pow3(s2t) + 7128*pow4(Mt) - 279*
        pow4(Mst2)*pow4(s2t))*pow8(Mst1) + 768*pow2(Mt)*(-2*pow2(Mt) + pow2(
        Mst2)*pow2(s2t))*pow8(Mst2))))))/(16200.*pow4(Mst1)*pow6(Mst2)))/
        pow4(Mt)/pow2(Sbeta)*12.; 

   return result;
}

/**
 * 	@return returns the susy log^2 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double himalaya::H6b::calc_coef_at_as2_no_sm_logs_log2() const {

   const double result =
      ((pow2(Sbeta)*(Mst2*pow2(Dmglst2)*(pow2(Mst2)*pow4(Mst1)*(33360*pow2(Mst2)
        *pow2(Mt)*pow2(s2t) + 52480*Mst2*s2t*pow3(Mt) + 18896*pow4(Mt) - 10005*
        pow4(Mst2)*pow4(s2t)) + 15*pow2(Mst1)*pow4(Mst2)*(-1496*pow2(Mst2)*
        pow2(Mt)*pow2(s2t) - 3072*Mst2*s2t*pow3(Mt) + 768*Mt*pow3(Mst2)*pow3(
        s2t) + 1456*pow4(Mt) + 475*pow4(Mst2)*pow4(s2t)) + 5*(120*pow2(Mst2)*
        pow2(Mt)*pow2(s2t) - 22656*Mst2*s2t*pow3(Mt) - 2304*Mt*pow3(Mst2)*pow3(
        s2t) + 52384*pow4(Mt) + 879*pow4(Mst2)*pow4(s2t))*pow6(Mst1) - 1440*
        pow2(-4*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow6(Mst2)) - 30*log(pow2(
        Mst2)/pow2(Mst1))*pow4(Mst1)*(Mst2*(pow2(Mst1)*pow2(Mst2)*(1096*pow2(
        Mst2)*pow2(Mt)*pow2(s2t) - 1280*Mst2*s2t*pow3(Mt) + 328*Mt*pow3(Mst2)*
        pow3(s2t) - 32*pow4(Mt) - 91*pow4(Mst2)*pow4(s2t)) + 4*pow4(Mst2)*(14*
        pow2(Mst2)*pow2(Mt)*pow2(s2t) - 456*Mst2*s2t*pow3(Mt) + 146*Mt*pow3(
        Mst2)*pow3(s2t) - 166*pow4(Mt) + 33*pow4(Mst2)*pow4(s2t)) - pow4(Mst1)*
        (-768*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 1216*Mst2*s2t*pow3(Mt) + 400*
        pow4(Mt) + 41*pow4(Mst2)*pow4(s2t))) - 2*Dmglst2*(32*pow2(Mt)*(-57*
        Mst2*Mt*s2t + 119*pow2(Mt) - 24*pow2(Mst2)*pow2(s2t))*pow4(Mst1) +
        pow4(Mst2)*(-384*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 912*Mst2*s2t*pow3(Mt)
        - 548*Mt*pow3(Mst2)*pow3(s2t) + 512*pow4(Mt) - 91*pow4(Mst2)*pow4(s2t))
        + pow2(Mst1)*pow2(Mst2)*(-768*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 640*Mst2*
        s2t*pow3(Mt) - 164*Mt*pow3(Mst2)*pow3(s2t) + 2016*pow4(Mt) + 91*pow4(
        Mst2)*pow4(s2t))) + Mst2*pow2(Dmglst2)*(384*pow2(Mt)*pow2(s2t)*pow4(
        Mst2) - 512*pow2(Mst2)*pow4(Mt) + pow2(Mst1)*(768*pow2(Mst2)*pow2(Mt)*
        pow2(s2t) - 1280*Mst2*s2t*pow3(Mt) + 4000*pow4(Mt) - 91*pow4(Mst2)*
        pow4(s2t)) + 768*Mt*pow3(s2t)*pow5(Mst2) + 91*pow4(s2t)*pow6(Mst2))) -
        10*Dmglst2*(pow4(Mst1)*pow4(Mst2)*(-8208*pow2(Mst2)*pow2(Mt)*pow2(s2t)
        - 4640*Mst2*s2t*pow3(Mt) - 1968*Mt*pow3(Mst2)*pow3(s2t) + 4816*pow4(Mt)
        + 849*pow4(Mst2)*pow4(s2t)) - 3*pow2(Mst2)*(-472*pow2(Mst2)*pow2(Mt)*
        pow2(s2t) + 2496*Mst2*s2t*pow3(Mt) - 1040*Mt*pow3(Mst2)*pow3(s2t) -
        3360*pow4(Mt) + 37*pow4(Mst2)*pow4(s2t))*pow6(Mst1) - 3*pow2(Mst1)*(-
        984*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 1536*Mst2*s2t*pow3(Mt) + 384*Mt*
        pow3(Mst2)*pow3(s2t) + 1456*pow4(Mt) + 219*pow4(Mst2)*pow4(s2t))*pow6(
        Mst2) + 3*(-2240*Mst2*s2t*pow3(Mt) + 3744*pow4(Mt) - 59*pow4(Mst2)*
        pow4(s2t))*pow8(Mst1) - 768*pow2(Mt)*(-2*pow2(Mt) + pow2(Mst2)*pow2(
        s2t))*pow8(Mst2)) + 5*Mst2*(-720*pow2(Dmsqst2)*pow4(Mst1)*pow4(Mt) +
        1440*Dmsqst2*pow2(Mst2)*pow4(Mst1)*pow4(Mt) + pow4(Mst1)*pow4(Mst2)*(
        8592*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 4800*Mst2*s2t*pow3(Mt) + 3936*Mt*
        pow3(Mst2)*pow3(s2t) + 4256*pow4(Mt) - 69*pow4(Mst2)*pow4(s2t)) - 3*
        pow2(Mst2)*(600*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 128*Mst2*s2t*pow3(Mt) +
        1568*Mt*pow3(Mst2)*pow3(s2t) + 672*pow4(Mt) + 355*pow4(Mst2)*pow4(s2t))
        *pow6(Mst1) + 3*pow2(Mst1)*(-856*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 1024*
        Mst2*s2t*pow3(Mt) + 256*Mt*pow3(Mst2)*pow3(s2t) + 1456*pow4(Mt) + 155*
        pow4(Mst2)*pow4(s2t))*pow6(Mst2) + (384*Mst2*s2t*pow3(Mt) - 2400*pow4(
        Mt) + 717*pow4(Mst2)*pow4(s2t))*pow8(Mst1) + 384*pow2(Mt)*(-2*pow2(Mt)
        + pow2(Mst2)*pow2(s2t))*pow8(Mst2))))/(540.*pow4(Mst1)*pow5(Mst2)))/
        pow4(Mt)/pow2(Sbeta)*12.; 

   return result;
}

/**
 * 	@return returns the susy log^3 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double himalaya::H6b::calc_coef_at_as2_no_sm_logs_log3() const {

   const double result =
      -896/3.; 

   return result;
}

