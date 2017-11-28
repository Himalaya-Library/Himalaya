#define Pi M_PI

#include "H6bq22g.hpp"
#include "HierarchyCalculator.hpp"
#include "Utils.hpp"
#include <type_traits>
#include <cmath>

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
 * 	@param Msq a double average squark mass w/o the stop quark
 * 	@param MuSUSY a double mu parameter
 * 	@param s2t a double 2 times the sine of the stop/sbottom quark mixing angle
 * 	@param mdrFlag an int 0 for DR and 1 for MDR scheme
 * 	@param oneLoopFlag an int flag to consider the one-loop expansion terms
 * 	@param twoLoopFlag an int flag to consider the two-loop expansion terms
 * 	@param threeLoopFlag an int flag to consider the three-loop expansion terms
 */
himalaya::H6bq22g::H6bq22g(std::map<unsigned int, unsigned int> flagMap, double Al4p, double beta, double Dmglst2,
		 double Dmsqst2, double lmMt, double lmMst1, double lmMst2,
		 double Mt, double Mst1, double Mst2, double Msq, double MuSUSY,
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
   this -> Msq = Msq;
   this -> MuSUSY = MuSUSY;
   this -> s2t = s2t;
   // zeta functions
   z2 = pow2(Pi)/6.;
   z3 = 1.202056903159594;
   z4 = pow4(Pi)/90.;
   // poly logs
   double pl412 = 0.51747906167389934317668576113647; // PolyLog[4,1/2]
   std::complex<double> pl2expPi3(0.27415567780803773941206919444, 1.014941606409653625021202554275); // PolyLog[2, Exp[I Pi / 3]]
   std::complex<double> pl3expPi6sqrt3(0.51928806536375962552715984277228, - 0.33358157526196370641686908633664); // PolyLog[3, Exp[- I Pi / 6] / Sqrt[3]]

   // polylog functions, checked
   B4 = (-4 * z2 * pow2(log(2)) + 2 / 3.* pow4(log(2)) - 13 / 2. * z4 + 16. * pl412);
   D3 = 6 * z3 - 15 / 4. * z4 - 6. * pow2(std::imag(pl2expPi3));
   DN = 6 * z3 - 4 * z2 * pow2(log(2)) + 2 / 3. * pow4(log(2)) - 21 / 2. * z4 + 16. * pl412;
   OepS2 = - 763 / 32. - (9 * Pi * sqrt(3) * pow2(log(3))) / 16. - (35 * pow3(Pi) * sqrt(3)) / 48.
      + 195 / 16. * z2 - 15 / 4. * z3 + 57 / 16. * z4 + 45 * sqrt(3) / 2. * std::imag(pl2expPi3)
      - 27 * sqrt(3) * std::imag(pl3expPi6sqrt3);
   S2 = 4 * std::imag(pl2expPi3) / (9. * sqrt(3));
   T1ep = - 45 / 2. - (Pi * sqrt(3) * pow2(log(3))) / 8. - (35 * pow3(Pi) * sqrt(3)) / 216. - 9 / 2. * z2 + z3 
      + 6. * sqrt(3) * std::imag(pl2expPi3) - 6. * sqrt(3) * std::imag(pl3expPi6sqrt3);
   // mdr flags, indicates if one wants to shift the dr stop mass to the mdr stop mass
   shiftst1 = mdrFlag;
   shiftst2 = mdrFlag;
   shiftst3 = mdrFlag;
   xDR2DRMOD = mdrFlag;
   // expansion flags
   xDmglst2 = flagMap.at(HierarchyCalculator::xxDmglst2);
   xDmsqst2 = flagMap.at(HierarchyCalculator::xxDmsqst2);
   xMst = flagMap.at(HierarchyCalculator::xxMst);
   
   s1 = 
   #include "../hierarchies/h6bq22g/sigS1Full.inc"
   ;
   s2 = 
   #include "../hierarchies/h6bq22g/sigS2Full.inc"
   ;
   s12 = 
   #include "../hierarchies/h6bq22g/sigS12Full.inc"
   ;
}

/**
 * 	@return The diagonal (1, 1) matrix element of the Higgs mass matrix as a double for the hierarchy 'H6bq22g'
 */
double himalaya::H6bq22g::getS1(){
   return s1;
}

/**
 * 	@return The diagonal (2, 2) matrix element of the Higgs mass matrix as a double for the hierarchy 'H6bq22g'
 */
double himalaya::H6bq22g::getS2(){
   return s2;
}

/**
 * 	@return The off-diagonal (1, 2) = (2, 1) matrix element of the Higgs mass matrix as a double for the hierarchy 'H6bq22g'
 */
double himalaya::H6bq22g::getS12(){
   return s12;
}

/**
 * 	@return returns the constant term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double himalaya::H6bq22g::calc_at_as2_no_logs(){

   const double result =
      ((pow2(Sbeta)*((30*Mst2*pow2(Dmsqst2)*pow2(Mst1)*(10290*Mst2*pow4(Mst1)*(-
        5400*Dmglst2*(8*(113 - 48*z2)*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 192*Mst2*
        s2t*(51*z2 - 28*(2 + z3))*pow3(Mt) + 96*Mt*(5*z2 + 3*(-4 + z3))*pow3(
        Mst2)*pow3(s2t) + 48*(-401 + 308*z2 - 144*z3)*pow4(Mt) - 3*(-75 + 32*z2
        + 8*z3)*pow4(Mst2)*pow4(s2t)) + Mst2*(-4*pow2(Mst2)*pow2(Mt)*pow2(s2t)*
        (-9928 + 1600*OepS2 - 361800*S2 - 2400*T1ep + 34500*z2 + 50925*z3 -
        1200*z4 - 1800*pow2(z2)) - 518400*Mst2*s2t*(-23 + 24*z2 - 16*z3)*pow3(
        Mt) - 259200*Mt*(4 + z2 - 3*z3)*pow3(Mst2)*pow3(s2t) + 864*(-23941 +
        20700*z2 - 10800*z3)*pow4(Mt) + (449186 + 20800*OepS2 - 3439800*S2 -
        31200*T1ep - 523500*z2 + 169275*z3 - 15600*z4 - 23400*pow2(z2))*pow4(
        Mst2)*pow4(s2t))) + 385875*pow2(Mst1)*pow3(Mst2)*(-64*Dmglst2*(36*(23 -
        24*z2)*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 24*Mst2*s2t*(289 - 342*z2 + 216*
        z3)*pow3(Mt) - 216*Mt*(-2 + 2*z2 - z3)*pow3(Mst2)*pow3(s2t) + 8*(-1477
        + 1242*z2 - 540*z3)*pow4(Mt) + 27*(-2 + 4*z2 - 3*z3)*pow4(Mst2)*pow4(
        s2t)) + Mst2*(-8*pow2(Mst2)*pow2(Mt)*pow2(s2t)*(-6878 + 128*OepS2 -
        16848*S2 - 192*T1ep - 696*z2 + 3273*z3 - 96*z4 - 144*pow2(z2)) - 3456*
        Mst2*s2t*(13 + 28*z2 - 48*z3)*pow3(Mt) - 6912*Mt*(-5 + 7*z2 - 7*z3)*
        pow3(Mst2)*pow3(s2t) + 16*(-11134 + 10800*z2 - 5481*z3)*pow4(Mt) + (-
        16678 + 256*OepS2 + 114048*S2 - 384*T1ep - 1392*z2 - 16521*z3 - 192*z4
        - 288*pow2(z2))*pow4(Mst2)*pow4(s2t))) - 166698000*(Mst2 + 4*Dmglst2*(-
        3 + z2) - 3*Mst2*z2)*pow2(-4*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow5(
        Mst2) + (-8*Mst2*pow2(Mt)*pow2(s2t)*(4630500*Dmglst2*(-1931 + 1098*z2 +
        180*z3) + Mst2*(1351886621 - 19208000*OepS2 + 639009000*S2 + 28812000*
        T1ep - 365552250*z2 - 225865500*z3 + 14406000*z4 + 21609000*pow2(z2)))
        + 2667168000*s2t*(Dmglst2*(-465 + 348*z2 - 176*z3) + Mst2*(35 - 20*z2 +
        16*z3))*pow3(Mt) - 166698000*Mt*(16*Dmglst2*(-20 + 11*z2) + Mst2*(101 -
        106*z2 + 60*z3))*pow2(Mst2)*pow3(s2t) + 27648*(-6582784 + 2701125*z2 -
        2315250*z3)*pow4(Mt) + 3*(12348000*Dmglst2*(38 + 63*z2 - 90*z3) + Mst2*
        (881319682 - 617400000*S2 - 472311000*z2 + 461892375*z3))*pow3(Mst2)*
        pow4(s2t))*pow6(Mst1)))/pow4(Msq) - (60*Mst2*pow2(Mst1)*pow3(Dmsqst2)*(
        -2058*Mst2*pow4(Mst1)*(-13500*Dmglst2*(8*(113 - 48*z2)*pow2(Mst2)*pow2(
        Mt)*pow2(s2t) - 192*Mst2*s2t*(51*z2 - 28*(2 + z3))*pow3(Mt) + 96*Mt*(5*
        z2 + 3*(-4 + z3))*pow3(Mst2)*pow3(s2t) + 48*(-401 + 308*z2 - 144*z3)*
        pow4(Mt) - 3*(-75 + 32*z2 + 8*z3)*pow4(Mst2)*pow4(s2t)) + Mst2*(4*pow2(
        Mst2)*pow2(Mt)*pow2(s2t)*(268334 + 16000*OepS2 - 54000*S2 - 24000*T1ep
        - 627000*z2 - 128625*z3 - 12000*z4 - 18000*pow2(z2)) - 2592000*Mst2*
        s2t*(-17 + 19*z2 - 12*z3)*pow3(Mt) + 1296000*Mt*(-3 + 2*z3)*pow3(Mst2)*
        pow3(s2t) + 576*(-133897 + 128250*z2 - 63000*z3)*pow4(Mt) + (2089333 +
        20000*OepS2 - 6952500*S2 - 30000*T1ep - 945750*z2 + 67875*z3 - 15000*z4
        - 22500*pow2(z2))*pow4(Mst2)*pow4(s2t))) + 83349000*(Mst2*(3 - 4*z2) +
        4*Dmglst2*(-3 + z2))*pow2(-4*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow5(
        Mst2) - 385875*pow2(Mst1)*pow3(Mst2)*(-1152*s2t*(Dmglst2*s2t*(23 - 24*
        z2) + Mt*(-3 + 72*z2 - 96*z3))*pow2(Mst2)*pow2(Mt) + 16*Mst2*(Mt*(-8662
        + 8208*z2 - 3753*z3) + 48*Dmglst2*s2t*(-289 + 342*z2 - 216*z3))*pow3(
        Mt) - 8*Mt*pow3(Mst2)*(Mt*(432 + pow2(s2t)*(-3359 + 32*OepS2 - 8100*S2
        - 48*T1ep + 690*z2 + 1185*z3 - 24*z4 - 36*pow2(z2))) - 864*Dmglst2*(-2
        + 2*z2 - z3)*pow3(s2t)) - 864*(8*Mt*(-4 + 5*z2 - 4*z3) + Dmglst2*s2t*(-
        2 + 4*z2 - 3*z3))*pow3(s2t)*pow4(Mst2) - 256*Dmglst2*(-1477 + 1242*z2 -
        540*z3)*pow4(Mt) + (-9400 + 64*OepS2 + 84888*S2 - 96*T1ep + 1380*z2 -
        14865*z3 - 48*z4 - 72*pow2(z2))*pow4(s2t)*pow5(Mst2)) + (8*Mst2*pow2(
        Mt)*pow2(s2t)*(2315250*Dmglst2*(-1931 + 1098*z2 + 180*z3) + Mst2*(
        1396747123 - 30184000*OepS2 + 1287279000*S2 + 45276000*T1ep -
        157694250*z2 - 283489500*z3 + 22638000*z4 + 33957000*pow2(z2))) -
        1333584000*s2t*(Dmglst2*(-465 + 348*z2 - 176*z3) + Mst2*(35 - 20*z2 +
        16*z3))*pow3(Mt) + 166698000*Mt*(8*Dmglst2*(-20 + 11*z2) + Mst2*(85 -
        98*z2 + 60*z3))*pow2(Mst2)*pow3(s2t) - 13824*(-6582784 + 2701125*z2 -
        2315250*z3)*pow4(Mt) + (-18522000*Dmglst2*(38 + 63*z2 - 90*z3) + Mst2*(
        -1130025178 + 38416000*OepS2 - 1166886000*S2 - 57624000*T1ep +
        189336000*z2 - 850597125*z3 - 28812000*z4 - 43218000*pow2(z2)))*pow3(
        Mst2)*pow4(s2t))*pow6(Mst1)))/pow6(Msq) + (120*Dmsqst2*Mst2*pow2(Mst1)*
        (3087*Mst2*pow4(Mst1)*(-4500*Dmglst2*(8*(113 - 48*z2)*pow2(Mst2)*pow2(
        Mt)*pow2(s2t) - 192*Mst2*s2t*(51*z2 - 28*(2 + z3))*pow3(Mt) + 96*Mt*(5*
        z2 + 3*(-4 + z3))*pow3(Mst2)*pow3(s2t) + 48*(-401 + 308*z2 - 144*z3)*
        pow4(Mt) - 3*(-75 + 32*z2 + 8*z3)*pow4(Mst2)*pow4(s2t)) + Mst2*(-8*
        pow2(Mst2)*pow2(Mt)*pow2(s2t)*(36449 + 4000*OepS2 - 310500*S2 - 6000*
        T1ep - 75750*z2 + 21000*z3 - 3000*z4 - 4500*pow2(z2)) - 864000*Mst2*
        s2t*(-6 + 5*z2 - 4*z3)*pow3(Mt) - 432000*Mt*(1 + z2 - z3)*pow3(Mst2)*
        pow3(s2t) + 96*(-91321 + 54000*z2 - 36000*z3)*pow4(Mt) + (52199 +
        28000*OepS2 - 3415500*S2 - 42000*T1ep - 557250*z2 + 259500*z3 - 21000*
        z4 - 31500*pow2(z2))*pow4(Mst2)*pow4(s2t))) + 385875*pow2(Mst1)*pow3(
        Mst2)*(-16*Dmglst2*(36*(23 - 24*z2)*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 24*
        Mst2*s2t*(289 - 342*z2 + 216*z3)*pow3(Mt) - 216*Mt*(-2 + 2*z2 - z3)*
        pow3(Mst2)*pow3(s2t) + 8*(-1477 + 1242*z2 - 540*z3)*pow4(Mt) + 27*(-2 +
        4*z2 - 3*z3)*pow4(Mst2)*pow4(s2t)) + 3*Mst2*(-4*pow2(Mst2)*pow2(Mt)*
        pow2(s2t)*(32*OepS2 - 3*(391 + 972*S2 + 16*T1ep + 154*z2 - 232*z3 + 8*
        z4 + 12*pow2(z2))) - 1152*Mst2*s2t*(7 + 2*z2 - 8*z3)*pow3(Mt) - 1152*
        Mt*(-1 + 2*z2 - 3*z3)*pow3(Mst2)*pow3(s2t) + 64*(-103 + 108*z2 - 72*z3)
        *pow4(Mt) + (-1213 + 32*OepS2 + 4860*S2 - 48*T1ep - 462*z2 - 276*z3 -
        24*z4 - 36*pow2(z2))*pow4(Mst2)*pow4(s2t))) + 41674500*(Mst2 - 4*
        Dmglst2*(-3 + z2) + 2*Mst2*z2)*pow2(-4*pow2(Mt) + pow2(Mst2)*pow2(s2t))
        *pow5(Mst2) + 2*(-4*Mst2*pow2(Mt)*pow2(s2t)*(1157625*Dmglst2*(-1931 +
        1098*z2 + 180*z3) + Mst2*(-22430251 + 5488000*OepS2 - 324135000*S2 -
        8232000*T1ep - 103929000*z2 + 28812000*z3 - 4116000*z4 - 6174000*pow2(
        z2))) + 333396000*s2t*(Dmglst2*(-465 + 348*z2 - 176*z3) + Mst2*(35 -
        20*z2 + 16*z3))*pow3(Mt) + 333396000*Mt*(Dmglst2*(20 - 11*z2) + Mst2*(-
        2 + z2))*pow2(Mst2)*pow3(s2t) + 3456*(-6582784 + 2701125*z2 - 2315250*
        z3)*pow4(Mt) + (4630500*Dmglst2*(38 + 63*z2 - 90*z3) + Mst2*(378483467
        + 9604000*OepS2 - 754771500*S2 - 14406000*T1ep - 306899250*z2 +
        133770000*z3 - 7203000*z4 - 10804500*pow2(z2)))*pow3(Mst2)*pow4(s2t))*
        pow6(Mst1)))/pow2(Msq) + 4*Dmglst2*pow2(Mst1)*(-1764*pow2(Mst2)*pow4(
        Mst1)*(-112*pow2(Mst2)*pow2(Mt)*pow2(s2t)*(1216808 + 150000*OepS2 -
        376650*S2 - 225000*T1ep - 3408750*z2 - 833625*z3 - 112500*z4 - 168750*
        pow2(z2)) + 1050*Mst2*s2t*(-322421 + 73760*OepS2 - 564876*S2 - 110640*
        T1ep - 3060894*z2 + 1606920*z3 - 55320*z4 - 82980*pow2(z2))*pow3(Mt) -
        175*Mt*(27211 + 36720*B4 + 1080*DN + 64160*OepS2 - 1515564*S2 - 96240*
        T1ep - 1240446*z2 + 245400*z3 - 12480*z4 - 72180*pow2(z2))*pow3(Mst2)*
        pow3(s2t) - 8*(98884301 + 4508000*OepS2 - 17853750*S2 - 6762000*T1ep -
        477659550*z2 + 302883000*z3 - 3381000*z4 - 5071500*pow2(z2))*pow4(Mt) -
        7*(-95041 + 81000*B4 - 108000*D3 + 67500*DN - 432000*OepS2 - 7425000*S2
        + 648000*T1ep + 5325300*z2 + 4227750*z3 - 121500*z4)*pow4(Mst2)*pow4(
        s2t)) - 66150*pow2(Mst1)*pow4(Mst2)*(672*(-11674 + 120*B4 - 120*D3 +
        60*DN + 17091*S2 + 4849*z2 + 2195*z3 - 540*z4)*pow2(Mst2)*pow2(Mt)*
        pow2(s2t) + 140*Mst2*s2t*(1120*OepS2 + 3*(-10061 + 14580*S2 - 560*T1ep
        - 25134*z2 + 72*z3 + 2024*z4 - 420*pow2(z2)))*pow3(Mt) - 7*Mt*(69997 +
        188640*B4 - 3600*DN + 5600*OepS2 + 146772*S2 - 8400*T1ep - 1043682*z2 +
        1816440*z3 + 32520*z4 - 317340*pow2(z2))*pow3(Mst2)*pow3(s2t) - 16*(
        1547066 + 7840*OepS2 + 258066*S2 - 11760*T1ep - 1256046*z2 - 946785*z3
        - 5880*z4 - 8820*pow2(z2))*pow4(Mt) + 35*(-15707 + 432*B4 - 576*D3 +
        360*DN + 224*OepS2 + 543348*S2 - 336*T1ep - 11802*z2 - 105690*z3 -
        2544*z4 - 2844*pow2(z2))*pow4(Mst2)*pow4(s2t)) + (72*pow2(Mst2)*pow2(
        Mt)*pow2(s2t)*(-3650009897 + 2110136000*OepS2 - 52385772600*S2 -
        3165204000*T1ep - 40405228500*z2 + 755286000*z3 - 1582602000*z4 -
        2373903000*pow2(z2)) - 274400*Mst2*s2t*(-20964397 + 2058400*OepS2 -
        47545272*S2 - 3087600*T1ep - 59449002*z2 + 39214920*z3 - 1543800*z4 -
        2315700*pow2(z2))*pow3(Mt) + 17150*Mt*(-31897243 + 2491360*OepS2 -
        90290268*S2 - 3737040*T1ep - 51499698*z2 + 33912840*z3 - 1868520*z4 -
        2802780*pow2(z2))*pow3(Mst2)*pow3(s2t) + 80*(-54946675289 + 3063401600*
        OepS2 - 52952863320*S2 - 4595102400*T1ep - 225236187120*z2 +
        243414477600*z3 - 2297551200*z4 - 3446326800*pow2(z2))*pow4(Mt) + (
        119405394763 - 11522056000*OepS2 + 478191735000*S2 + 17283084000*T1ep +
        233890773900*z2 - 111792103500*z3 + 8641542000*z4 + 12962313000*pow2(
        z2))*pow4(Mst2)*pow4(s2t))*pow6(Mst1) - 333396000*(2*Mt + Mst2*s2t)*(
        14*Mt*(3 + z2) + Mst2*s2t*(5 + 7*z2))*pow2(-2*Mt + Mst2*s2t)*pow6(Mst2)
        ) + 588*Mst2*pow2(Dmglst2)*(1134000*pow2(Mst1)*pow4(Mst2)*(-8*(-179 +
        26*z2)*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 768*Mst2*s2t*pow3(Mt) - 192*Mt*
        pow3(Mst2)*pow3(s2t) + 16*(-155 + 26*z2)*pow4(Mt) + (-235 + 26*z2)*
        pow4(Mst2)*pow4(s2t)) + 15*pow2(Mst2)*pow4(Mst1)*(-360*pow2(Mst2)*pow2(
        Mt)*pow2(s2t)*(-2164753 + 3360*B4 - 3360*D3 + 1680*DN - 7840*OepS2 +
        169452*S2 + 11760*T1ep + 467818*z2 + 1945895*z3 - 9240*z4 + 8820*pow2(
        z2)) + 320*Mst2*s2t*(-1243547 + 15680*OepS2 + 953208*S2 - 23520*T1ep -
        2985438*z2 + 3393915*z3 - 11760*z4 - 17640*pow2(z2))*pow3(Mt) + 560*Mt*
        (-282298 + 61560*B4 - 1620*DN - 2240*OepS2 - 111213*S2 + 3360*T1ep -
        281922*z2 + 983190*z3 + 6540*z4 - 114120*pow2(z2))*pow3(Mst2)*pow3(s2t)
        - 16*(-46632377 + 744800*OepS2 + 29679210*S2 - 1117200*T1ep -
        107181930*z2 + 108350025*z3 - 558600*z4 - 837900*pow2(z2))*pow4(Mt) -
        35*(-924689 + 41040*B4 - 43200*D3 + 22680*DN - 1120*OepS2 + 21325356*S2
        + 1680*T1ep - 127686*z2 - 3687510*z3 - 190320*z4 - 37620*pow2(z2))*
        pow4(Mst2)*pow4(s2t)) - 2*(-36*pow2(Mst2)*pow2(Mt)*pow2(s2t)*(68526791
        + 2772000*OepS2 - 107403300*S2 - 4158000*T1ep - 140027250*z2 +
        39514125*z3 - 2079000*z4 - 3118500*pow2(z2)) - 800*Mst2*s2t*(7384247 +
        1183840*OepS2 + 5821092*S2 - 1775760*T1ep - 59988552*z2 + 25757760*z3 -
        887880*z4 - 1331820*pow2(z2))*pow3(Mt) + 350*Mt*(-2727217 + 437920*
        OepS2 - 3648132*S2 - 656880*T1ep - 12301638*z2 + 4166400*z3 - 328440*z4
        - 492660*pow2(z2))*pow3(Mst2)*pow3(s2t) + 8*(648916519 + 158732000*
        OepS2 + 66753450*S2 - 238098000*T1ep - 10760276550*z2 + 6504855000*z3 -
        119049000*z4 - 178573500*pow2(z2))*pow4(Mt) - 7*(-21452821 + 243000*B4
        - 324000*D3 + 202500*DN + 2272000*OepS2 - 75724200*S2 - 3408000*T1ep -
        55081500*z2 + 67149750*z3 - 3040500*z4 - 4014000*pow2(z2))*pow4(Mst2)*
        pow4(s2t))*pow6(Mst1) + 18144000*pow2(-4*pow2(Mt) + pow2(Mst2)*pow2(
        s2t))*pow6(Mst2)) + (30*Mst2*pow2(Mst1)*pow4(Dmsqst2)*(-6174*Mst2*pow4(
        Mst1)*(9000*Dmglst2*(8*(113 - 48*z2)*pow2(Mst2)*pow2(Mt)*pow2(s2t) -
        192*Mst2*s2t*(51*z2 - 28*(2 + z3))*pow3(Mt) + 96*Mt*(5*z2 + 3*(-4 + z3)
        )*pow3(Mst2)*pow3(s2t) + 48*(-401 + 308*z2 - 144*z3)*pow4(Mt) - 3*(-75
        + 32*z2 + 8*z3)*pow4(Mst2)*pow4(s2t)) + Mst2*(-12*pow2(Mst2)*pow2(Mt)*
        pow2(s2t)*(113744 + 8000*OepS2 - 225000*S2 - 12000*T1ep - 259500*z2 -
        28875*z3 - 6000*z4 - 9000*pow2(z2)) + 864000*Mst2*s2t*(-45 + 52*z2 -
        32*z3)*pow3(Mt) - 432000*Mt*(-8 + z2 + 5*z3)*pow3(Mst2)*pow3(s2t) - 96*
        (-712061 + 715500*z2 - 342000*z3)*pow4(Mt) + (-2037134 + 8000*OepS2 +
        3537000*S2 - 12000*T1ep + 388500*z2 + 191625*z3 - 6000*z4 - 9000*pow2(
        z2))*pow4(Mst2)*pow4(s2t))) + 166698000*(-4*Dmglst2*(-3 + z2) + 5*Mst2*
        (-1 + z2))*pow2(-4*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow5(Mst2) +
        385875*pow2(Mst1)*pow3(Mst2)*(-1152*s2t*(2*Dmglst2*s2t*(23 - 24*z2) +
        3*Mt*(-17 + 68*z2 - 80*z3))*pow2(Mst2)*pow2(Mt) + 48*Mst2*(Mt*(-7838 +
        7344*z2 - 3177*z3) + 32*Dmglst2*s2t*(-289 + 342*z2 - 216*z3))*pow3(Mt)
        + 24*Mt*pow3(Mst2)*(Mt*(-864 + (2186 + 5184*S2 - 1152*z2 - 489*z3)*
        pow2(s2t)) + 576*Dmglst2*(-2 + 2*z2 - z3)*pow3(s2t)) - 1728*(4*Mt*(-11
        + 13*z2 - 9*z3) + Dmglst2*s2t*(-2 + 4*z2 - 3*z3))*pow3(s2t)*pow4(Mst2)
        - 512*Dmglst2*(-1477 + 1242*z2 - 540*z3)*pow4(Mt) + 3*(-6974 + 75168*S2
        + 2304*z2 - 14313*z3)*pow4(s2t)*pow5(Mst2)) + (-56*Mst2*pow2(Mt)*pow2(
        s2t)*(661500*Dmglst2*(-1931 + 1098*z2 + 180*z3) + Mst2*(605014553 -
        14504000*OepS2 + 644301000*S2 + 21756000*T1ep - 37889250*z2 -
        129727500*z3 + 10878000*z4 + 16317000*pow2(z2))) + 2667168000*s2t*(
        Dmglst2*(-465 + 348*z2 - 176*z3) + Mst2*(35 - 20*z2 + 16*z3))*pow3(Mt)
        - 166698000*Mt*(16*Dmglst2*(-20 + 11*z2) + Mst2*(239 - 286*z2 + 180*z3)
        )*pow2(Mst2)*pow3(s2t) + 27648*(-6582784 + 2701125*z2 - 2315250*z3)*
        pow4(Mt) + 7*(5292000*Dmglst2*(38 + 63*z2 - 90*z3) + Mst2*(268020238 -
        21952000*OepS2 + 931392000*S2 + 32928000*T1ep + 94227000*z2 +
        288101625*z3 + 16464000*z4 + 24696000*pow2(z2)))*pow3(Mst2)*pow4(s2t))*
        pow6(Mst1)))/pow8(Msq) + 3*Mst2*(308700*pow4(Mst1)*pow4(Mst2)*(2*pow2(
        Mst2)*pow2(Mt)*pow2(s2t)*(1013263 + 2880*D3 - 2880*DN - 25440*OepS2 -
        137052*S2 + 38160*T1ep + 35562*z2 + 83160*z3 + 36360*z4 + 28620*pow2(
        z2)) + 24*Mst2*s2t*(-60759 + 1120*OepS2 + 25596*S2 - 1680*T1ep -
        137386*z2 + 115320*z3 - 12360*z4 - 1260*pow2(z2))*pow3(Mt) + 30*Mt*(
        25611 + 5280*B4 - 48*DN - 224*OepS2 - 324*S2 + 336*T1ep - 13614*z2 +
        47720*z3 + 2040*z4 - 6660*pow2(z2))*pow3(Mst2)*pow3(s2t) + 8*(389597 +
        3360*OepS2 - 105948*S2 - 5040*T1ep + 293898*z2 - 740160*z3 - 2520*z4 -
        3780*pow2(z2))*pow4(Mt) - 5*(25289 + 1440*B4 - 144*D3 + 72*DN - 2208*
        OepS2 + 298404*S2 + 3312*T1ep + 36318*z2 - 94968*z3 + 1440*z4 - 108*
        pow2(z2))*pow4(Mst2)*pow4(s2t)) + 20580*pow2(Mst2)*(4*pow2(Mst2)*pow2(
        Mt)*pow2(s2t)*(2367149 - 21600*D3 + 21600*DN - 503200*OepS2 + 14889420*
        S2 + 754800*T1ep + 9292470*z2 - 66000*z3 + 247800*z4 + 566100*pow2(z2))
        + 40*Mst2*s2t*(-1156193 + 60320*OepS2 - 1414908*S2 - 90480*T1ep -
        2823942*z2 + 2907240*z3 - 45240*z4 - 67860*pow2(z2))*pow3(Mt) + 100*Mt*
        (51041 + 7344*B4 + 216*DN - 2336*OepS2 + 112428*S2 + 3504*T1ep + 28038*
        z2 - 78216*z3 + 8880*z4 + 2628*pow2(z2))*pow3(Mst2)*pow3(s2t) + 8*(
        7319011 + 167200*OepS2 - 6744060*S2 - 250800*T1ep + 9301170*z2 -
        20715000*z3 - 125400*z4 - 188100*pow2(z2))*pow4(Mt) + (-7680127 +
        10800*B4 - 21600*D3 + 16200*DN + 514400*OepS2 - 46307700*S2 - 771600*
        T1ep - 10868970*z2 + 7186200*z3 - 466800*z4 - 773100*pow2(z2))*pow4(
        Mst2)*pow4(s2t))*pow6(Mst1) + 55566000*pow2(Mst1)*(-8*(87 + 38*z2)*
        pow2(Mst2)*pow2(Mt)*pow2(s2t) - 1024*Mst2*s2t*pow3(Mt) + 256*Mt*pow3(
        Mst2)*pow3(s2t) + 16*(71 + 38*z2)*pow4(Mt) + (135 + 38*z2)*pow4(Mst2)*
        pow4(s2t))*pow6(Mst2) + (-16*pow2(Mst2)*pow2(Mt)*pow2(s2t)*(-
        25002374272 + 4703902000*OepS2 - 184815757350*S2 - 7055853000*T1ep -
        84661721025*z2 + 16805113500*z3 - 3527926500*z4 - 5291889750*pow2(z2))
        + 1097600*Mst2*s2t*(-2016907 + 109600*OepS2 - 3520152*S2 - 164400*T1ep
        - 3448554*z2 + 3465480*z3 - 82200*z4 - 123300*pow2(z2))*pow3(Mt) -
        68600*Mt*(-2759935 + 70240*OepS2 - 3305340*S2 - 105360*T1ep - 629034*z2
        + 2677800*z3 - 52680*z4 - 79020*pow2(z2))*pow3(Mst2)*pow3(s2t) + 16*(
        187545955063 + 3204992000*OepS2 - 128216692800*S2 - 4807488000*T1ep +
        129009640200*z2 - 399103824000*z3 - 2403744000*z4 - 3605616000*pow2(z2)
        )*pow4(Mt) + (93508520089 + 2000376000*B4 + 222264000*D3 - 222264000*DN
        + 4925480000*OepS2 - 312095700000*S2 - 7388220000*T1ep - 77106571500*z2
        + 21506100000*z3 - 2360526000*z4 - 5541165000*pow2(z2))*pow4(Mst2)*
        pow4(s2t))*pow8(Mst1) + 7112448000*pow2(Mt)*(-2*pow2(Mt) + pow2(Mst2)*
        pow2(s2t))*pow8(Mst2))))/(6.001128e9*pow4(Mst1)*pow5(Mst2)))/pow4(Mt)/
        pow2(Sbeta)*12.;
 
   return result;
}

/**
 * 	@return returns the susy log^0 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double himalaya::H6bq22g::calc_coef_at_as2_no_sm_logs_log0(){

   const double result =
      (pow2(Sbeta)*(-(z4*((4*Dmglst2*(pow4(Mst1)*(55368*pow2(Mst2)*pow2(Mt)*
        pow2(s2t) - 205840*Mst2*s2t*pow3(Mt) + 15571*Mt*pow3(Mst2)*pow3(s2t) +
        89312*pow4(Mt) - 4199*pow4(Mst2)*pow4(s2t)) + 9*pow2(Mst1)*pow2(Mst2)*(
        1200*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 5532*Mst2*s2t*pow3(Mt) + 208*Mt*
        pow3(Mst2)*pow3(s2t) + 2576*pow4(Mt) + 81*pow4(Mst2)*pow4(s2t)) - 27*
        pow4(Mst2)*(432*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 1012*Mst2*s2t*pow3(Mt)
        + 271*Mt*pow3(Mst2)*pow3(s2t) - 112*pow4(Mt) + 106*pow4(Mst2)*pow4(s2t)
        )) - 6*Mst2*pow2(Dmglst2)*(pow2(Mst1)*(-7128*pow2(Mst2)*pow2(Mt)*pow2(
        s2t) - 67648*Mst2*s2t*pow3(Mt) + 10948*Mt*pow3(Mst2)*pow3(s2t) + 90704*
        pow4(Mt) - 2027*pow4(Mst2)*pow4(s2t)) + 6*pow2(Mst2)*(396*pow2(Mst2)*
        pow2(Mt)*pow2(s2t) - 448*Mst2*s2t*pow3(Mt) + 436*Mt*pow3(Mst2)*pow3(
        s2t) + 1064*pow4(Mt) + 793*pow4(Mst2)*pow4(s2t))))*pow8(Msq) + 3*Mst2*(
        20*Dmsqst2*pow2(Mst2)*pow2(s2t)*(pow2(Mst1)*pow3(Dmsqst2)*(108*pow2(
        Mst2)*pow2(Mt) + 4*pow2(Mst1)*(37*pow2(Mt) - 7*pow2(Mst2)*pow2(s2t)) -
        9*pow2(s2t)*pow4(Mst2)) + pow2(Dmsqst2)*pow2(Msq)*(2*(44*pow2(Mt) - 7*
        pow2(Mst2)*pow2(s2t))*pow4(Mst1) + 9*(-4*pow2(Mt) + pow2(Mst2)*pow2(
        s2t))*pow4(Mst2) + 3*pow2(Mst1)*(16*pow2(Mst2)*pow2(Mt) + 5*pow2(s2t)*
        pow4(Mst2))) + Dmsqst2*pow4(Msq)*(28*pow2(Mt)*pow4(Mst1) + 18*(-4*pow2(
        Mt) + pow2(Mst2)*pow2(s2t))*pow4(Mst2) + 3*pow2(Mst1)*(-4*pow2(Mst2)*
        pow2(Mt) + 13*pow2(s2t)*pow4(Mst2))) + pow6(Msq)*((-32*pow2(Mt) + 14*
        pow2(Mst2)*pow2(s2t))*pow4(Mst1) + 9*pow2(Mst1)*(-8*pow2(Mst2)*pow2(Mt)
        + 7*pow2(s2t)*pow4(Mst2)) + 27*(-4*pow2(Mt)*pow4(Mst2) + pow2(s2t)*
        pow6(Mst2)))) + (108*pow4(Mst2)*(-101*pow2(Mst2)*pow2(Mt)*pow2(s2t) +
        412*Mst2*s2t*pow3(Mt) - 85*Mt*pow3(Mst2)*pow3(s2t) + 28*pow4(Mt) + 10*
        pow4(Mst2)*pow4(s2t)) + 12*pow2(Mst1)*pow2(Mst2)*(-826*pow2(Mst2)*pow2(
        Mt)*pow2(s2t) + 1508*Mst2*s2t*pow3(Mt) - 740*Mt*pow3(Mst2)*pow3(s2t) +
        836*pow4(Mt) + 389*pow4(Mst2)*pow4(s2t)) + pow4(Mst1)*(-27428*pow2(
        Mst2)*pow2(Mt)*pow2(s2t) + 43840*Mst2*s2t*pow3(Mt) - 1756*Mt*pow3(Mst2)
        *pow3(s2t) + 18688*pow4(Mt) + 1147*pow4(Mst2)*pow4(s2t)))*pow8(Msq))))/
        (2916.*pow5(Mst2)*pow8(Msq)) - (T1ep*((4*Dmglst2*(pow4(Mst1)*(55368*
        pow2(Mst2)*pow2(Mt)*pow2(s2t) - 205840*Mst2*s2t*pow3(Mt) + 15571*Mt*
        pow3(Mst2)*pow3(s2t) + 89312*pow4(Mt) - 4199*pow4(Mst2)*pow4(s2t)) -
        189*pow4(Mst2)*(20*Mst2*s2t*pow3(Mt) - 5*Mt*pow3(Mst2)*pow3(s2t) - 16*
        pow4(Mt) + pow4(Mst2)*pow4(s2t)) - 18*pow2(Mst1)*pow2(Mst2)*(-600*pow2(
        Mst2)*pow2(Mt)*pow2(s2t) + 2766*Mst2*s2t*pow3(Mt) - 401*Mt*pow3(Mst2)*
        pow3(s2t) - 1288*pow4(Mt) + 108*pow4(Mst2)*pow4(s2t))) + 6*Mst2*pow2(
        Dmglst2)*(21*pow2(Mst2)*(72*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 128*Mst2*
        s2t*pow3(Mt) - 32*Mt*pow3(Mst2)*pow3(s2t) - 304*pow4(Mt) + pow4(Mst2)*
        pow4(s2t)) + 4*pow2(Mst1)*(1782*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 16912*
        Mst2*s2t*pow3(Mt) - 2737*Mt*pow3(Mst2)*pow3(s2t) - 22676*pow4(Mt) +
        284*pow4(Mst2)*pow4(s2t))))*pow8(Msq) + 3*Mst2*(20*Dmsqst2*pow2(Mst2)*
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
        pow4(Mt) + 1795*pow4(Mst2)*pow4(s2t)))*pow8(Msq))))/(1458.*pow5(Mst2)*
        pow8(Msq)) - (pow2(z2)*(16*pow3(Mt)*(6*Mst2*pow2(Dmglst2)*((-5669*Mt +
        4228*Mst2*s2t)*pow2(Mst1) + 21*(-19*Mt + 8*Mst2*s2t)*pow2(Mst2)) +
        Dmglst2*(-9*(-644*Mt + 1383*Mst2*s2t)*pow2(Mst1)*pow2(Mst2) + 4*(5582*
        Mt - 12865*Mst2*s2t)*pow4(Mst1) + 189*(4*Mt - 5*Mst2*s2t)*pow4(Mst2)) +
        3*Mst2*(3*(209*Mt + 377*Mst2*s2t)*pow2(Mst1)*pow2(Mst2) + 4*(292*Mt +
        685*Mst2*s2t)*pow4(Mst1) + 189*(Mt + Mst2*s2t)*pow4(Mst2))) + 4*Mt*
        pow3(Mst2)*pow3(s2t)*(-6*pow2(Dmglst2)*(2737*Mst2*pow2(Mst1) - 7608*
        pow3(Mst2)) + Dmglst2*(7218*pow2(Mst1)*pow2(Mst2) + 15571*pow4(Mst1) +
        47601*pow4(Mst2)) - 3*(438*pow2(Mst1)*pow3(Mst2) + 439*Mst2*pow4(Mst1)
        - 4995*pow5(Mst2))) + (pow4(Mst2)*pow4(s2t)*(-2*Dmglst2*(-5352*Dmglst2*
        Mst2*pow2(Mst1) + 9*(209*Dmglst2 + 474*Mst2)*pow3(Mst2) + 8398*pow4(
        Mst1))*pow8(Msq) + 3*Mst2*(60*pow2(Dmsqst2)*pow2(Mst2)*(13*pow2(Mst1) +
        6*pow2(Mst2))*pow4(Msq) - 20*pow4(Dmsqst2)*(9*pow2(Mst1)*pow2(Mst2) +
        28*pow4(Mst1)) - 20*pow2(Msq)*pow3(Dmsqst2)*(-15*pow2(Mst1)*pow2(Mst2)
        + 14*pow4(Mst1) - 9*pow4(Mst2)) + 20*Dmsqst2*(63*pow2(Mst1)*pow2(Mst2)
        + 14*pow4(Mst1) + 27*pow4(Mst2))*pow6(Msq) + (5154*pow2(Mst1)*pow2(
        Mst2) + 1795*pow4(Mst1) - 54*pow4(Mst2))*pow8(Msq))))/pow8(Msq) + (12*
        pow2(Mst2)*pow2(Mt)*pow2(s2t)*(12*Dmglst2*(3*Mst2*(99*Dmglst2 + 100*
        Mst2)*pow2(Mst1) + 63*Dmglst2*pow3(Mst2) + 1538*pow4(Mst1))*pow8(Msq) +
        Mst2*(20*pow2(Dmsqst2)*(pow2(Dmsqst2)*(27*pow2(Mst1)*pow2(Mst2) + 37*
        pow4(Mst1)) + pow4(Msq)*(-3*pow2(Mst1)*pow2(Mst2) + 7*pow4(Mst1) - 18*
        pow4(Mst2)) + Dmsqst2*pow2(Msq)*(12*pow2(Mst1)*pow2(Mst2) + 22*pow4(
        Mst1) - 9*pow4(Mst2))) - 20*Dmsqst2*(18*pow2(Mst1)*pow2(Mst2) + 8*pow4(
        Mst1) + 27*pow4(Mst2))*pow6(Msq) - (3774*pow2(Mst1)*pow2(Mst2) + 6857*
        pow4(Mst1) + 1431*pow4(Mst2))*pow8(Msq))))/pow8(Msq)))/(1944.*pow5(
        Mst2)) + (s2t*pow3(Mt)*((-4354560*(-1 + log(pow2(Mst2)/pow2(Mst1)))*(
        Dmglst2 - Mst2 + (3*Dmglst2 + Mst2)*log(pow2(Mst2)/pow2(Mst1)))*pow2(
        Mst2))/pow2(Mst1) + (4*(24*Mst2*pow2(Dmglst2)*((295960*OepS2 + 1455273*
        S2)*pow2(Mst1) + 6*(1960*OepS2 + 119151*S2)*pow2(Mst2)) - 7*Dmglst2*(
        27*(18440*OepS2 - 141219*S2)*pow2(Mst1)*pow2(Mst2) + 8*(257300*OepS2 -
        5943159*S2)*pow4(Mst1) + 675*(56*OepS2 + 2187*S2)*pow4(Mst2)) + 21*
        Mst2*(3*(15080*OepS2 - 353727*S2)*pow2(Mst1)*pow2(Mst2) + 8*(13700*
        OepS2 - 440019*S2)*pow4(Mst1) + 27*(280*OepS2 + 6399*S2)*pow4(Mst2)) -
        5670*S2*log(pow2(Mst2)/pow2(Mst1))*(3393*pow2(Mst1)*pow3(Mst2) + 168*
        pow2(Dmglst2)*(151*Mst2*pow2(Mst1) + 6*pow3(Mst2)) + 8220*Mst2*pow4(
        Mst1) - Dmglst2*(12447*pow2(Mst1)*pow2(Mst2) + 51460*pow4(Mst1) + 945*
        pow4(Mst2)) + 567*pow5(Mst2))))/pow4(Mst2) + (12*pow2(Mst1)*(2*pow2(
        Dmglst2)*(7384247 - 3649905*log(pow2(Mst2)/pow2(Mst1)) - 4704210*pow2(
        log(pow2(Mst2)/pow2(Mst1))) + 362880*pow3(log(pow2(Mst2)/pow2(Mst1))))
        - 7*pow2(Mst1)*(2016907 + 180360*log(pow2(Mst2)/pow2(Mst1)) + 610200*
        pow2(log(pow2(Mst2)/pow2(Mst1))) - 131760*pow3(log(pow2(Mst2)/pow2(
        Mst1))) - (24300*Dmsqst2*(35 - 50*log(pow2(Mst2)/pow2(Mst1)) + 4*pow2(
        log(pow2(Mst2)/pow2(Mst1))))*(pow2(Dmsqst2)*pow2(Msq) + pow3(Dmsqst2) +
        Dmsqst2*pow4(Msq) + pow6(Msq)))/pow8(Msq))))/pow3(Mst2) + 7*Dmglst2*(
        4074705 - 9069840*log(pow2(Mst2)/pow2(Mst1)) + 2579040*pow2(log(pow2(
        Mst2)/pow2(Mst1))) + 1905120*pow3(log(pow2(Mst2)/pow2(Mst1))) + (27*
        pow2(Mst1)*(322421 - 456320*log(pow2(Mst2)/pow2(Mst1)) + 526560*pow2(
        log(pow2(Mst2)/pow2(Mst1))) - 23040*pow3(log(pow2(Mst2)/pow2(Mst1))) +
        (345600*Dmsqst2*(-7 + 6*log(pow2(Mst2)/pow2(Mst1)))*(pow2(Dmsqst2)*
        pow2(Msq) + pow3(Dmsqst2) + Dmsqst2*pow4(Msq) + pow6(Msq)))/pow8(Msq)))
        /pow2(Mst2) + (4*pow4(Mst1)*(20964397 - 1700640*log(pow2(Mst2)/pow2(
        Mst1)) + 4837320*pow2(log(pow2(Mst2)/pow2(Mst1))) - 408240*pow3(log(
        pow2(Mst2)/pow2(Mst1))) - (72900*Dmsqst2*(465 - 374*log(pow2(Mst2)/
        pow2(Mst1)) + 12*pow2(log(pow2(Mst2)/pow2(Mst1))))*(pow2(Dmsqst2)*pow2(
        Msq) + pow3(Dmsqst2) + Dmsqst2*pow4(Msq) + pow6(Msq)))/pow8(Msq)))/
        pow4(Mst2) + (64800*Dmsqst2*(-289 + 354*log(pow2(Mst2)/pow2(Mst1)))*(
        pow2(Dmsqst2)*pow2(Msq) + pow3(Dmsqst2) + Dmsqst2*pow4(Msq) + pow6(Msq)
        ))/pow8(Msq)) + (9*(8*pow2(Dmglst2)*(-1243547*pow2(Mst1) + 181440*pow2(
        Mst2))*pow8(Msq) + 30240*pow2(Mst1)*(48*pow2(Mst1) + 49*pow2(Mst2))*
        pow3(log(pow2(Mst2)/pow2(Mst1)))*pow8(Msq) - 10080*pow2(log(pow2(Mst2)/
        pow2(Mst1)))*(465*pow2(Mst1)*pow2(Mst2) + pow2(Dmglst2)*(-386*pow2(
        Mst1) + 144*pow2(Mst2)) + 612*pow4(Mst1))*pow8(Msq) + 7*pow2(Mst1)*(
        32400*pow2(Msq)*(68*pow2(Mst1) + pow2(Mst2))*pow3(Dmsqst2) + 16200*(
        180*pow2(Mst1) + 17*pow2(Mst2))*pow4(Dmsqst2) + 16200*pow2(Dmsqst2)*(
        92*pow2(Mst1) - 13*pow2(Mst2))*pow4(Msq) + 64800*Dmsqst2*(12*pow2(Mst1)
        - 7*pow2(Mst2))*pow6(Msq) - (1156193*pow2(Mst1) + 546831*pow2(Mst2))*
        pow8(Msq)) - 1680*log(pow2(Mst2)/pow2(Mst1))*(8*pow2(Dmglst2)*(1091*
        pow2(Mst1) + 144*pow2(Mst2))*pow8(Msq) + 9*pow2(Mst1)*(150*pow2(Msq)*(
        8*pow2(Mst1) + 5*pow2(Mst2))*pow3(Dmsqst2) + 195*(8*pow2(Mst1) + 5*
        pow2(Mst2))*pow4(Dmsqst2) + 105*pow2(Dmsqst2)*(8*pow2(Mst1) + 5*pow2(
        Mst2))*pow4(Msq) + 60*Dmsqst2*(8*pow2(Mst1) + 5*pow2(Mst2))*pow6(Msq) +
        (63*pow2(Mst1) + 13*pow2(Mst2))*pow8(Msq)))))/(Mst2*pow2(Mst1)*pow8(
        Msq))))/153090. + Mt*pow3(s2t)*((64*(-1 + log(pow2(Mst2)/pow2(Mst1)))*(
        Dmglst2 - Mst2 + (3*Dmglst2 + Mst2)*log(pow2(Mst2)/pow2(Mst1)))*pow4(
        Mst2))/(9.*pow2(Mst1)) + (-6*Mst2*pow2(Dmglst2)*(17*(6440*OepS2 -
        53649*S2)*pow2(Mst1) + 3*(2240*OepS2 + 111213*S2)*pow2(Mst2)) - 15*
        Mst2*(6*(584*OepS2 - 28107*S2)*pow2(Mst1)*pow2(Mst2) + (3512*OepS2 -
        165267*S2)*pow4(Mst1) + 27*(56*OepS2 + 81*S2)*pow4(Mst2)) + Dmglst2*(
        18*(16040*OepS2 - 378891*S2)*pow2(Mst1)*pow2(Mst2) + (622840*OepS2 -
        22572567*S2)*pow4(Mst1) + 27*(1400*OepS2 + 36693*S2)*pow4(Mst2)) + 810*
        S2*log(pow2(Mst2)/pow2(Mst1))*(42*pow2(Dmglst2)*(391*Mst2*pow2(Mst1) +
        24*pow3(Mst2)) + 3*Mst2*(438*pow2(Mst1)*pow2(Mst2) + 439*pow4(Mst1) +
        189*pow4(Mst2)) - Dmglst2*(7218*pow2(Mst1)*pow2(Mst2) + 15571*pow4(
        Mst1) + 945*pow4(Mst2))))/(21870.*pow2(Mst2)) + Dmglst2*(pow2(Mst1)*(
        5.598971193415638 + (68*B4)/9. + (2*DN)/9. + (1006*log(pow2(Mst2)/pow2(
        Mst1)))/27. - (40*pow2(log(pow2(Mst2)/pow2(Mst1))))/3. - (260*pow3(log(
        pow2(Mst2)/pow2(Mst1))))/9. + (80*Dmsqst2*(12 - 5*log(pow2(Mst2)/pow2(
        Mst1)))*(pow2(Dmsqst2)*pow2(Msq) + pow3(Dmsqst2) + Dmsqst2*pow4(Msq) +
        pow6(Msq)))/(3.*pow8(Msq))) - (pow4(Mst1)*(364.6232624599909 + (18652*
        log(pow2(Mst2)/pow2(Mst1)))/243. + (328*pow2(log(pow2(Mst2)/pow2(Mst1))
        ))/9. + (224*pow3(log(pow2(Mst2)/pow2(Mst1))))/9. + (20*Dmsqst2*(-40 +
        3*log(pow2(Mst2)/pow2(Mst1)))*(pow2(Dmsqst2)*pow2(Msq) + pow3(Dmsqst2)
        + Dmsqst2*pow4(Msq) + pow6(Msq)))/(3.*pow8(Msq))))/pow2(Mst2) + pow2(
        Mst2)*(21.604012345679013 + (524*B4)/9. - (10*DN)/9. - (1640*log(pow2(
        Mst2)/pow2(Mst1)))/9. + (736*pow2(log(pow2(Mst2)/pow2(Mst1))))/9. - (
        172*pow3(log(pow2(Mst2)/pow2(Mst1))))/3. - (80*Dmsqst2*(2 + 5*log(pow2(
        Mst2)/pow2(Mst1)))*(pow2(Dmsqst2)*pow2(Msq) + pow3(Dmsqst2) + Dmsqst2*
        pow4(Msq) + pow6(Msq)))/(3.*pow8(Msq)))) + (2*pow2(Dmglst2)*(24*(-
        141149 + 30780*B4 - 810*DN)*pow2(Mst1)*pow2(Mst2) + 2727217*pow4(Mst1)
        - 311040*pow4(Mst2))*pow8(Msq) - 4320*pow2(Mst1)*pow3(log(pow2(Mst2)/
        pow2(Mst1)))*(288*pow2(Dmglst2)*pow2(Mst2) + 195*pow2(Mst1)*pow2(Mst2)
        + 200*pow4(Mst1) + 195*pow4(Mst2))*pow8(Msq) - 8640*pow2(log(pow2(Mst2)
        /pow2(Mst1)))*(3*pow2(Dmglst2)*(77*pow2(Mst1)*pow2(Mst2) + 52*pow4(
        Mst1) - 24*pow4(Mst2)) - pow2(Mst1)*(183*pow2(Mst1)*pow2(Mst2) + 125*
        pow4(Mst1) + 354*pow4(Mst2)))*pow8(Msq) + 5*pow2(Mst1)*(-4860*pow4(
        Dmsqst2)*(128*pow2(Mst1)*pow2(Mst2) + 239*pow4(Mst1) - 176*pow4(Mst2))
        - 4860*pow2(Dmsqst2)*pow4(Msq)*(64*pow2(Mst1)*pow2(Mst2) + 101*pow4(
        Mst1) - 80*pow4(Mst2)) - 9720*pow2(Msq)*pow3(Dmsqst2)*(48*pow2(Mst1)*
        pow2(Mst2) + 85*pow4(Mst1) - 64*pow4(Mst2)) - 155520*Dmsqst2*(pow2(
        Mst1)*pow2(Mst2) + pow4(Mst1) - pow4(Mst2))*pow6(Msq) + (6*(51041 +
        7344*B4 + 216*DN)*pow2(Mst1)*pow2(Mst2) + 551987*pow4(Mst1) + 81*(8537
        + 1760*B4 - 16*DN)*pow4(Mst2))*pow8(Msq)) - 180*log(pow2(Mst2)/pow2(
        Mst1))*(8*pow2(Dmglst2)*(-1317*pow2(Mst1)*pow2(Mst2) + 2479*pow4(Mst1)
        - 576*pow4(Mst2))*pow8(Msq) + 3*pow2(Mst1)*(45*pow4(Dmsqst2)*(-128*
        pow2(Mst1)*pow2(Mst2) + 355*pow4(Mst1) - 128*pow4(Mst2)) + 45*pow2(
        Dmsqst2)*pow4(Msq)*(-64*pow2(Mst1)*pow2(Mst2) + 113*pow4(Mst1) - 64*
        pow4(Mst2)) + 270*pow2(Msq)*pow3(Dmsqst2)*(-16*pow2(Mst1)*pow2(Mst2) +
        39*pow4(Mst1) - 16*pow4(Mst2)) - 360*Dmsqst2*pow2(pow2(Mst1) + 2*pow2(
        Mst2))*pow6(Msq) - 4*(663*pow2(Mst1)*pow2(Mst2) + 23*pow4(Mst1) - 2184*
        pow4(Mst2))*pow8(Msq))))/(29160.*Mst2*pow2(Mst1)*pow8(Msq))) + pow4(Mt)
        *(480.98395061728394 + (2212*log(pow2(Mst2)/pow2(Mst1)))/81. + (3872*
        pow2(log(pow2(Mst2)/pow2(Mst1))))/27. + (40*Dmsqst2*(-103 + 138*log(
        pow2(Mst2)/pow2(Mst1)) + 9*pow2(log(pow2(Mst2)/pow2(Mst1)))))/(27.*
        pow2(Msq)) + (4*pow2(Dmglst2)*(-620 + 166*log(pow2(Mst2)/pow2(Mst1)) +
        59*pow2(log(pow2(Mst2)/pow2(Mst1)))))/(9.*pow2(Mst1)) - (8*pow3(log(
        pow2(Mst2)/pow2(Mst1))))/9. + (5*pow2(Dmsqst2)*(-5567 + 5502*log(pow2(
        Mst2)/pow2(Mst1)) + 108*pow2(log(pow2(Mst2)/pow2(Mst1)))))/(81.*pow4(
        Msq)) - (64*(-1 + log(pow2(Mst2)/pow2(Mst1)))*(-Mst2 + (4*Dmglst2 +
        Mst2)*log(pow2(Mst2)/pow2(Mst1)))*pow3(Mst2))/(9.*pow4(Mst1)) - (2*(3*
        Mst2*pow2(Dmglst2)*((3174640*OepS2 + 1335069*S2)*pow2(Mst1) + 3*(74480*
        OepS2 + 2967921*S2)*pow2(Mst2)) - 21*Mst2*(3*(8360*OepS2 - 337203*S2)*
        pow2(Mst1)*pow2(Mst2) + 8*(5840*OepS2 - 233631*S2)*pow4(Mst1) + 27*(
        280*OepS2 - 8829*S2)*pow4(Mst2)) - 2*Dmglst2*(1035*(784*OepS2 - 3105*
        S2)*pow2(Mst1)*pow2(Mst2) + (3125920*OepS2 - 54033534*S2)*pow4(Mst1) +
        27*(3920*OepS2 + 129033*S2)*pow4(Mst2)) - 5670*S2*log(pow2(Mst2)/pow2(
        Mst1))*(6*pow2(Dmglst2)*(5669*Mst2*pow2(Mst1) + 399*pow3(Mst2)) - 4*
        Dmglst2*(1449*pow2(Mst1)*pow2(Mst2) + 5582*pow4(Mst1) + 189*pow4(Mst2))
        - 3*(627*pow2(Mst1)*pow3(Mst2) + 1168*Mst2*pow4(Mst1) + 189*pow5(Mst2))
        )))/(76545.*pow5(Mst2)) + (10*(-4331 + 3846*log(pow2(Mst2)/pow2(Mst1)))
        *pow3(Dmsqst2))/(81.*pow6(Msq)) - (pow2(Dmglst2)*pow2(Mst1)*(
        1017.3098475406623 + (983572*log(pow2(Mst2)/pow2(Mst1)))/2025. - (
        63356*pow2(log(pow2(Mst2)/pow2(Mst1))))/27. + (536*pow3(log(pow2(Mst2)/
        pow2(Mst1))))/3.) - (pow4(Mst1)*(187545955063 + 89674479720*log(pow2(
        Mst2)/pow2(Mst1)) + 9086099400*pow2(log(pow2(Mst2)/pow2(Mst1))) +
        7408800000*pow3(log(pow2(Mst2)/pow2(Mst1))) - (8640*Dmsqst2*(13165568 -
        18387705*log(pow2(Mst2)/pow2(Mst1)) + 1885275*pow2(log(pow2(Mst2)/pow2(
        Mst1))))*(pow2(Dmsqst2)*pow2(Msq) + pow3(Dmsqst2) + Dmsqst2*pow4(Msq) +
        pow6(Msq)))/pow8(Msq)))/1.250235e8)/pow4(Mst2) - (5*(3919 - 3294*log(
        pow2(Mst2)/pow2(Mst1)) + 36*pow2(log(pow2(Mst2)/pow2(Mst1))))*pow4(
        Dmsqst2))/(27.*pow8(Msq)) - Dmglst2*(((54946675289 + 53298152460*log(
        pow2(Mst2)/pow2(Mst1)) + 32279083200*pow2(log(pow2(Mst2)/pow2(Mst1))) -
        5156524800*pow3(log(pow2(Mst2)/pow2(Mst1))))*pow4(Mst1))/(1.8753525e7*
        pow5(Mst2)) - (4*(773533 + 229915*log(pow2(Mst2)/pow2(Mst1)) - 199920*
        pow2(log(pow2(Mst2)/pow2(Mst1))) + 40320*pow3(log(pow2(Mst2)/pow2(Mst1)
        )) + (350*Dmsqst2*(2954 - 1713*log(pow2(Mst2)/pow2(Mst1)))*(pow2(
        Dmsqst2)*pow2(Msq) + pow3(Dmsqst2) + Dmsqst2*pow4(Msq) + pow6(Msq)))/
        pow8(Msq)))/(2835.*Mst2) - (pow2(Mst1)*(930.1286396237507 - (68648*log(
        pow2(Mst2)/pow2(Mst1)))/675. - (47624*pow2(log(pow2(Mst2)/pow2(Mst1))))
        /45. + (5264*pow3(log(pow2(Mst2)/pow2(Mst1))))/27. + (40*Dmsqst2*(401 -
        274*log(pow2(Mst2)/pow2(Mst1)) + 4*pow2(log(pow2(Mst2)/pow2(Mst1))))*(
        pow2(Dmsqst2)*pow2(Msq) + pow3(Dmsqst2) + Dmsqst2*pow4(Msq) + pow6(Msq)
        ))/(3.*pow8(Msq))))/pow3(Mst2) - (8*Mst2*(180*Dmsqst2*(pow2(Dmsqst2)*
        pow2(Msq) + pow3(Dmsqst2) + Dmsqst2*pow4(Msq) + pow6(Msq)) + (-84 - 70*
        log(pow2(Mst2)/pow2(Mst1)) + 59*pow2(log(pow2(Mst2)/pow2(Mst1))))*pow8(
        Msq)))/(9.*pow2(Mst1)*pow8(Msq))) + (25200*(96*pow2(Dmglst2) + 55*pow2(
        Mst1))*pow3(log(pow2(Mst2)/pow2(Mst1)))*pow4(Mst1)*pow8(Msq) + 2*pow2(
        Dmglst2)*(46632377*pow4(Mst1) + 1209600*pow4(Mst2))*pow8(Msq) - 7*pow2(
        Mst1)*(540*pow2(Dmsqst2)*pow4(Msq)*(23941*pow4(Mst1) + 300*pow4(Mst2))
        + 144*pow2(Msq)*pow3(Dmsqst2)*(133897*pow4(Mst1) + 3375*pow4(Mst2)) +
        36*pow4(Dmsqst2)*(712061*pow4(Mst1) + 22500*pow4(Mst2)) + 72*Dmsqst2*(
        91321*pow4(Mst1) - 2250*pow4(Mst2))*pow6(Msq) - (7319011*pow4(Mst1) +
        383400*pow4(Mst2))*pow8(Msq)) - 168*log(pow2(Mst2)/pow2(Mst1))*(2*pow2(
        Dmglst2)*(56837*pow4(Mst1) + 3600*pow4(Mst2))*pow8(Msq) - 5*pow2(Mst1)*
        (135*pow2(Dmsqst2)*pow4(Msq)*(667*pow4(Mst1) + 20*pow4(Mst2)) + 36*
        pow2(Msq)*pow3(Dmsqst2)*(3469*pow4(Mst1) + 75*pow4(Mst2)) + 9*pow4(
        Dmsqst2)*(17747*pow4(Mst1) + 300*pow4(Mst2)) + 18*Dmsqst2*(3067*pow4(
        Mst1) + 150*pow4(Mst2))*pow6(Msq) + 2*(11176*pow4(Mst1) - 2745*pow4(
        Mst2))*pow8(Msq))) + 2520*pow2(log(pow2(Mst2)/pow2(Mst1)))*(4*pow2(
        Dmglst2)*(2069*pow4(Mst1) - 360*pow4(Mst2))*pow8(Msq) + 5*pow2(Mst1)*(-
        18*pow2(Msq)*pow3(Dmsqst2)*pow4(Mst1) + 54*pow4(Dmsqst2)*pow4(Mst1) -
        90*pow2(Dmsqst2)*pow4(Msq)*pow4(Mst1) - 162*Dmsqst2*pow4(Mst1)*pow6(
        Msq) + (1136*pow4(Mst1) + 177*pow4(Mst2))*pow8(Msq))))/(85050.*pow2(
        Mst2)*pow4(Mst1)*pow8(Msq))) + pow4(s2t)*(pow4(Mst1)*(
        46.745471895783595 + B4 + D3/9. - DN/9. + (23468297*log(pow2(Mst2)/
        pow2(Mst1)))/529200. + (11764*pow2(log(pow2(Mst2)/pow2(Mst1))))/945. +
        (Dmsqst2*(378483467 + 391652730*log(pow2(Mst2)/pow2(Mst1)) + 63107100*
        pow2(log(pow2(Mst2)/pow2(Mst1)))))/(2.50047e7*pow2(Msq)) - (409*pow3(
        log(pow2(Mst2)/pow2(Mst1))))/108. + (pow2(Dmsqst2)*(440659841 +
        531394290*log(pow2(Mst2)/pow2(Mst1)) + 24607800*pow2(log(pow2(Mst2)/
        pow2(Mst1)))))/(3.33396e7*pow4(Msq)) + ((11.298127731986387 + (1287107*
        log(pow2(Mst2)/pow2(Mst1)))/79380. - (22*pow2(log(pow2(Mst2)/pow2(Mst1)
        )))/21.)*pow3(Dmsqst2))/pow6(Msq) + (Dmglst2*(79.58863384550371 - (
        8287903*log(pow2(Mst2)/pow2(Mst1)))/1.1907e6 - (51*pow2(log(pow2(Mst2)/
        pow2(Mst1))))/70. - (11*pow3(log(pow2(Mst2)/pow2(Mst1))))/9. + (5*
        Dmsqst2*(76 - 249*log(pow2(Mst2)/pow2(Mst1)))*(pow2(Dmsqst2)*pow2(Msq)
        + pow3(Dmsqst2) + Dmsqst2*pow4(Msq) + pow6(Msq)))/(54.*pow8(Msq))))/
        Mst2 + ((9.378945088323395 + (373997*log(pow2(Mst2)/pow2(Mst1)))/22680.
         - (17*pow2(log(pow2(Mst2)/pow2(Mst1))))/6.)*pow4(Dmsqst2))/pow8(Msq))
        + Dmglst2*(-(pow2(Mst1)*(Dmglst2*(29.427737997256514 - B4/3. + (4*D3)/
        9. - (5*DN)/18. + (332311*log(pow2(Mst2)/pow2(Mst1)))/8100. + (159*
        pow2(log(pow2(Mst2)/pow2(Mst1))))/20. + (13*pow3(log(pow2(Mst2)/pow2(
        Mst1))))/9.) + Mst2*(0.7822304526748971 - (2*B4)/3. + (8*D3)/9. - (5*
        DN)/9. + (81643*log(pow2(Mst2)/pow2(Mst1)))/4050. + (4969*pow2(log(
        pow2(Mst2)/pow2(Mst1))))/270. + (58*pow3(log(pow2(Mst2)/pow2(Mst1))))/
        27. + (5*Dmsqst2*(75 - 26*log(pow2(Mst2)/pow2(Mst1)))*(pow2(Dmsqst2)*
        pow2(Msq) + pow3(Dmsqst2) + Dmsqst2*pow4(Msq) + pow6(Msq)))/(6.*pow8(
        Msq))))) + (pow3(Mst2)*(2160*Dmsqst2*(2*(1 + log(pow2(Mst2)/pow2(Mst1))
        )*pow2(Mst1) + 3*pow2(Mst2))*(pow2(Dmsqst2)*pow2(Msq) + pow3(Dmsqst2) +
        Dmsqst2*pow4(Msq) + pow6(Msq)) + ((15707 - 432*B4 + 576*D3 - 360*DN)*
        pow2(Mst1) + 72*log(pow2(Mst2)/pow2(Mst1))*(277*pow2(Mst1) - 131*pow2(
        Mst2)) - 720*pow2(Mst2) + 36*(284*pow2(Mst1) + 187*pow2(Mst2))*pow2(
        log(pow2(Mst2)/pow2(Mst1))) - 7704*pow2(Mst1)*pow3(log(pow2(Mst2)/pow2(
        Mst1))))*pow8(Msq)))/(648.*pow2(Mst1)*pow8(Msq))) - (pow2(Mst2)*(5*(-
        924689 + 41040*B4 - 43200*D3 + 22680*DN)*pow2(Dmglst2)*pow8(Msq) +
        1800*(321*pow2(Dmglst2) + 230*pow2(Mst1))*pow3(log(pow2(Mst2)/pow2(
        Mst1)))*pow8(Msq) - pow2(Mst1)*(4178666*pow2(Msq)*pow3(Dmsqst2) +
        6111402*pow4(Dmsqst2) + 2245930*pow2(Dmsqst2)*pow4(Msq) + 313194*
        Dmsqst2*pow6(Msq) + (-7680127 + 10800*B4 - 21600*D3 + 16200*DN)*pow8(
        Msq)) + 1800*pow2(log(pow2(Mst2)/pow2(Mst1)))*(891*pow2(Dmglst2)*pow8(
        Msq) - pow2(Mst1)*(171*pow2(Msq)*pow3(Dmsqst2) + 117*pow4(Dmsqst2) +
        225*pow2(Dmsqst2)*pow4(Msq) + 279*Dmsqst2*pow6(Msq) + 184*pow8(Msq))) -
        60*log(pow2(Mst2)/pow2(Mst1))*(41640*pow2(Dmglst2)*pow8(Msq) + pow2(
        Mst1)*(-5316*pow2(Msq)*pow3(Dmsqst2) - 10827*pow4(Dmsqst2) + 195*pow2(
        Dmsqst2)*pow4(Msq) + 5706*Dmsqst2*pow6(Msq) + 32102*pow8(Msq)))))/(
        97200.*pow8(Msq)) - (pow4(Mst2)*(144*pow2(Dmglst2)*(235*pow2(Mst1) -
        16*pow2(Mst2))*pow8(Msq) + 9756*pow3(log(pow2(Mst2)/pow2(Mst1)))*pow4(
        Mst1)*pow8(Msq) + pow2(Mst1)*(40*pow2(Msq)*(1175*pow2(Mst1) + 81*pow2(
        Mst2))*pow3(Dmsqst2) + 15*(3487*pow2(Mst1) + 360*pow2(Mst2))*pow4(
        Dmsqst2) + 5*pow2(Dmsqst2)*(8339*pow2(Mst1) + 216*pow2(Mst2))*pow4(Msq)
        + 30*Dmsqst2*(1213*pow2(Mst1) - 36*pow2(Mst2))*pow6(Msq) + ((25289 +
        1440*B4 - 144*D3 + 72*DN)*pow2(Mst1) - 4860*pow2(Mst2))*pow8(Msq)) -
        36*pow2(log(pow2(Mst2)/pow2(Mst1)))*(-120*pow2(Msq)*pow3(Dmsqst2)*pow4(
        Mst1) - 120*pow4(Dmsqst2)*pow4(Mst1) - 120*pow2(Dmsqst2)*pow4(Msq)*
        pow4(Mst1) - 120*Dmsqst2*pow4(Mst1)*pow6(Msq) + pow2(Dmglst2)*(443*
        pow2(Mst1) - 96*pow2(Mst2))*pow8(Msq) + pow2(Mst1)*(742*pow2(Mst1) +
        123*pow2(Mst2))*pow8(Msq)) - 24*log(pow2(Mst2)/pow2(Mst1))*(3*pow2(
        Dmglst2)*(51*pow2(Mst1) - 16*pow2(Mst2))*pow8(Msq) - pow2(Mst1)*(90*
        pow2(Msq)*(6*pow2(Mst1) - pow2(Mst2))*pow3(Dmsqst2) + 15*(37*pow2(Mst1)
        - 6*pow2(Mst2))*pow4(Dmsqst2) + 15*pow2(Dmsqst2)*(35*pow2(Mst1) - 6*
        pow2(Mst2))*pow4(Msq) + 30*Dmsqst2*(17*pow2(Mst1) - 3*pow2(Mst2))*pow6(
        Msq) + (1148*pow2(Mst1) + 375*pow2(Mst2))*pow8(Msq)))))/(1296.*pow4(
        Mst1)*pow8(Msq)) + (3*Mst2*pow2(Dmglst2)*(4*(11360*OepS2 - 378621*S2)*
        pow2(Mst1) + 3*(280*OepS2 - 5331339*S2)*pow2(Mst2))*pow8(Msq) - 10*
        Dmglst2*(972*(16*OepS2 + 275*S2)*pow2(Mst1)*pow2(Mst2) + (33592*OepS2 -
        1394145*S2)*pow4(Mst1) + 27*(56*OepS2 + 135837*S2)*pow4(Mst2))*pow8(
        Msq) + 15*Mst2*(-10*pow4(Dmsqst2)*(9*(8*OepS2 + 3537*S2)*pow2(Mst1)*
        pow2(Mst2) + 32*(7*OepS2 - 297*S2)*pow4(Mst1) - 126846*S2*pow4(Mst2)) -
        30*pow2(Dmsqst2)*pow4(Msq)*(-13*(8*OepS2 - 1323*S2)*pow2(Mst1)*pow2(
        Mst2) + 900*S2*pow4(Mst1) - 24*(2*OepS2 + 891*S2)*pow4(Mst2)) - 10*
        pow2(Msq)*pow3(Dmsqst2)*(-15*(8*OepS2 - 2781*S2)*pow2(Mst1)*pow2(Mst2)
        + 14*(8*OepS2 - 243*S2)*pow4(Mst1) - 9*(8*OepS2 + 10611*S2)*pow4(Mst2))
        + 10*Dmsqst2*(9*(56*OepS2 - 6831*S2)*pow2(Mst1)*pow2(Mst2) + 2*(56*
        OepS2 - 4401*S2)*pow4(Mst1) + 27*(8*OepS2 + 1215*S2)*pow4(Mst2))*pow6(
        Msq) + (3*(5144*OepS2 - 463077*S2)*pow2(Mst1)*pow2(Mst2) + 10*(718*
        OepS2 - 45495*S2)*pow4(Mst1) + 27*(184*OepS2 - 24867*S2)*pow4(Mst2))*
        pow8(Msq)) - 405*S2*log(pow2(Mst2)/pow2(Mst1))*(6*Mst2*pow2(Dmglst2)*(
        1136*pow2(Mst1) + 21*pow2(Mst2))*pow8(Msq) - 4*Dmglst2*(1944*pow2(Mst1)
        *pow2(Mst2) + 4199*pow4(Mst1) + 189*pow4(Mst2))*pow8(Msq) + 3*Mst2*(60*
        pow2(Dmsqst2)*pow2(Mst2)*(13*pow2(Mst1) + 6*pow2(Mst2))*pow4(Msq) - 20*
        pow4(Dmsqst2)*(9*pow2(Mst1)*pow2(Mst2) + 28*pow4(Mst1)) - 20*pow2(Msq)*
        pow3(Dmsqst2)*(-15*pow2(Mst1)*pow2(Mst2) + 14*pow4(Mst1) - 9*pow4(Mst2)
        ) + 20*Dmsqst2*(63*pow2(Mst1)*pow2(Mst2) + 14*pow4(Mst1) + 27*pow4(
        Mst2))*pow6(Msq) + (3858*pow2(Mst1)*pow2(Mst2) + 1795*pow4(Mst1) +
        1242*pow4(Mst2))*pow8(Msq))))/(43740.*Mst2*pow8(Msq))) - z3*(Mt*pow3(
        s2t)*(((69440*pow2(Dmglst2)*pow2(Mst1))/243. + pow4(Mst1)*(
        91.83127572016461 + (50*pow2(Dmsqst2))/pow4(Msq) + (100*pow3(Dmsqst2))/
        pow6(Msq) + (150*pow4(Dmsqst2))/pow8(Msq)))/Mst2 - Dmglst2*((282607*
        pow4(Mst1))/(729.*pow2(Mst2)) + pow2(Mst1)*(50.49382716049383 - (80*
        Dmsqst2*(pow2(Dmsqst2)*pow2(Msq) + pow3(Dmsqst2) + Dmsqst2*pow4(Msq) +
        pow6(Msq)))/pow8(Msq)) + pow2(Mst2)*(560.6296296296297 - (80*Dmsqst2*(
        pow2(Dmsqst2)*pow2(Msq) + pow3(Dmsqst2) + Dmsqst2*pow4(Msq) + pow6(Msq)
        ))/(3.*pow8(Msq)))) - (Mst2*(4320*pow2(Msq)*(pow2(Mst1) + 2*pow2(Mst2))
        *pow3(Dmsqst2) + 1080*(5*pow2(Mst1) + 9*pow2(Mst2))*pow4(Dmsqst2) +
        1080*pow2(Dmsqst2)*(3*pow2(Mst1) + 7*pow2(Mst2))*pow4(Msq) + 2160*
        Dmsqst2*(pow2(Mst1) + 3*pow2(Mst2))*pow6(Msq) + 65546*pow2(Dmglst2)*
        pow8(Msq) + (-6518*pow2(Mst1) + 17895*pow2(Mst2))*pow8(Msq)))/(81.*
        pow8(Msq))) - (pow4(s2t)*(12*Mst2*pow2(Dmglst2)*(179066*pow2(Mst1) +
        368751*pow2(Mst2))*pow8(Msq) - 8*Dmglst2*(9720*Dmsqst2*(-2*pow2(Mst1)*
        pow2(Mst2) + 5*pow4(Mst1) - 3*pow4(Mst2))*(pow2(Dmsqst2)*pow2(Msq) +
        pow3(Dmsqst2) + Dmsqst2*pow4(Msq) + pow6(Msq)) + (-101466*pow2(Mst1)*
        pow2(Mst2) + 217283*pow4(Mst1) - 475605*pow4(Mst2))*pow8(Msq)) + 3*
        Mst2*(5*pow4(Dmsqst2)*(-9198*pow2(Mst1)*pow2(Mst2) + 15679*pow4(Mst1) -
        128817*pow4(Mst2)) + 10*pow2(Msq)*pow3(Dmsqst2)*(1086*pow2(Mst1)*pow2(
        Mst2) + 6613*pow4(Mst1) - 44595*pow4(Mst2)) + 15*pow2(Dmsqst2)*pow4(
        Msq)*(4514*pow2(Mst1)*pow2(Mst2) + 3591*pow4(Mst1) - 16521*pow4(Mst2))
        + 80*Dmsqst2*(1557*pow2(Mst1)*pow2(Mst2) + 520*pow4(Mst1) - 621*pow4(
        Mst2))*pow6(Msq) + 16*(35931*pow2(Mst1)*pow2(Mst2) + 5225*pow4(Mst1) -
        162*log(pow2(Mst2)/pow2(Mst1))*(pow4(Mst1) - pow4(Mst2)) + 35613*pow4(
        Mst2))*pow8(Msq))))/(23328.*Mst2*pow8(Msq)) + pow4(Mt)*(
        913.7777777777778 + (16*log(pow2(Mst2)/pow2(Mst1)))/3. + (320*Dmsqst2)/
        (3.*pow2(Msq)) + (1015*pow2(Dmsqst2))/(6.*pow4(Msq)) + (695*pow3(
        Dmsqst2))/(3.*pow6(Msq)) + ((2478040*pow2(Dmglst2)*pow2(Mst1))/243. +
        pow4(Mst1)*(3192.230452674897 + (320*Dmsqst2*(pow2(Dmsqst2)*pow2(Msq) +
        pow3(Dmsqst2) + Dmsqst2*pow4(Msq) + pow6(Msq)))/pow8(Msq)))/pow4(Mst2)
        - Dmglst2*((9462176*pow4(Mst1))/(729.*pow5(Mst2)) - (2*(9017 - (7200*
        Dmsqst2*(pow2(Dmsqst2)*pow2(Msq) + pow3(Dmsqst2) + Dmsqst2*pow4(Msq) +
        pow6(Msq)))/pow8(Msq)))/(27.*Mst2) + (pow2(Mst1)*(2848.9876543209875 +
        (1920*Dmsqst2*(pow2(Dmsqst2)*pow2(Msq) + pow3(Dmsqst2) + Dmsqst2*pow4(
        Msq) + pow6(Msq)))/pow8(Msq)))/pow3(Mst2)) + (1765*pow4(Dmsqst2))/(6.*
        pow8(Msq)) + (206381*pow2(Dmglst2)*pow8(Msq) + 20*pow2(Mst1)*(3024*
        pow2(Msq)*pow3(Dmsqst2) + 4104*pow4(Dmsqst2) + 1944*pow2(Dmsqst2)*pow4(
        Msq) + 864*Dmsqst2*pow6(Msq) + 6905*pow8(Msq)))/(81.*pow2(Mst2)*pow8(
        Msq))) + (s2t*pow2(Mt)*(-3*Mst2*s2t*(18*Mst2*pow2(Dmglst2)*(15053*pow2(
        Mst1) - 55597*pow2(Mst2))*pow8(Msq) + Mst2*(5*pow2(Dmsqst2)*pow4(Msq)*(
        -2037*pow2(Mst1)*pow2(Mst2) + 1756*pow4(Mst1) - 9819*pow4(Mst2)) + 5*
        pow4(Dmsqst2)*(-2079*pow2(Mst1)*pow2(Mst2) + 7060*pow4(Mst1) - 4401*
        pow4(Mst2)) + 10*pow2(Msq)*pow3(Dmsqst2)*(-1029*pow2(Mst1)*pow2(Mst2) +
        2204*pow4(Mst1) - 3555*pow4(Mst2)) - 80*Dmsqst2*(126*pow2(Mst1)*pow2(
        Mst2) + 56*pow4(Mst1) + 783*pow4(Mst2))*pow6(Msq) - 4*(660*pow2(Mst1)*
        pow2(Mst2) + 648*log(pow2(Mst2)/pow2(Mst1))*pow2(Mst2)*(-pow2(Mst1) +
        pow2(Mst2)) + 32663*pow4(Mst1) - 6237*pow4(Mst2))*pow8(Msq)) - 48*
        Dmglst2*(675*Dmsqst2*pow4(Mst1)*(pow2(Dmsqst2)*pow2(Msq) + pow3(
        Dmsqst2) + Dmsqst2*pow4(Msq) + pow6(Msq)) - (-2223*pow2(Mst1)*pow2(
        Mst2) + 734*pow4(Mst1) - 1317*pow4(Mst2))*pow8(Msq))) - 16*Mt*(3*Mst2*
        pow2(Dmglst2)*(245312*pow2(Mst1) + 96969*pow2(Mst2))*pow8(Msq) -
        Dmglst2*(38880*Dmsqst2*(7*pow2(Mst1)*pow2(Mst2) + 11*pow4(Mst1) + 3*
        pow4(Mst2))*(pow2(Dmsqst2)*pow2(Msq) + pow3(Dmsqst2) + Dmsqst2*pow4(
        Msq) + pow6(Msq)) + (361557*pow2(Mst1)*pow2(Mst2) + 1307164*pow4(Mst1)
        + 243*pow4(Mst2))*pow8(Msq)) + 3*Mst2*(6480*(pow2(Dmsqst2)*pow4(Msq)*(
        4*pow2(Mst1)*pow2(Mst2) + 2*pow4(Mst1) + 3*pow4(Mst2)) + pow4(Dmsqst2)*
        (8*pow2(Mst1)*pow2(Mst2) + 2*pow4(Mst1) + 5*pow4(Mst2))) + 12960*(pow2(
        Msq)*pow3(Dmsqst2)*(3*pow2(Mst1)*pow2(Mst2) + pow4(Mst1) + 2*pow4(Mst2)
        ) + Dmsqst2*(pow2(Mst1)*pow2(Mst2) + pow4(Mst1) + pow4(Mst2))*pow6(Msq)
        ) + (72681*pow2(Mst1)*pow2(Mst2) + 115516*pow4(Mst1) + 25947*pow4(Mst2)
        )*pow8(Msq)))))/(2916.*pow4(Mst2)*pow8(Msq))) + z2*(-(pow2(Mt)*pow2(
        s2t)*((247.52275132275133 - (4204*log(pow2(Mst2)/pow2(Mst1)))/9.)*pow2(
        Dmglst2) + Dmglst2*((8*pow2(Mst1)*(505 - 124*log(pow2(Mst2)/pow2(Mst1))
        - (120*Dmsqst2*(pow2(Dmsqst2)*pow2(Msq) + pow3(Dmsqst2) + Dmsqst2*pow4(
        Msq) + pow6(Msq)))/pow8(Msq)))/(9.*Mst2) + Mst2*(143.67407407407407 -
        264*log(pow2(Mst2)/pow2(Mst1)) - (320*Dmsqst2*(pow2(Dmsqst2)*pow2(Msq)
        + pow3(Dmsqst2) + Dmsqst2*pow4(Msq) + pow6(Msq)))/(3.*pow8(Msq))) + (2*
        pow4(Mst1)*(78533 - 180*log(pow2(Mst2)/pow2(Mst1)) + (135*Dmsqst2*(61 -
        10*log(pow2(Mst2)/pow2(Mst1)))*(pow2(Dmsqst2)*pow2(Msq) + pow3(Dmsqst2)
        + Dmsqst2*pow4(Msq) + pow6(Msq)))/pow8(Msq)))/(81.*pow3(Mst2))) - (
        pow2(Mst2)*(5927 + 36600*log(pow2(Mst2)/pow2(Mst1)) + (23100*Dmsqst2)/
        pow2(Msq) - (12480*pow2(Dmglst2))/pow2(Mst1) + (5800*pow2(Dmsqst2))/
        pow4(Msq) - (11500*pow3(Dmsqst2))/pow6(Msq) - (28800*pow4(Dmsqst2))/
        pow8(Msq)))/540. + (pow2(Mst1)*(41800*pow2(Msq)*pow3(Dmsqst2) + 77850*
        pow4(Dmsqst2) + 5750*pow2(Dmsqst2)*pow4(Msq) - 30300*Dmsqst2*pow6(Msq)
        - (309749 + 12240*log(pow2(Mst2)/pow2(Mst1)))*pow8(Msq)))/(810.*pow8(
        Msq)) - (4*pow3(Mst2)*(4*Dmglst2*(15*Dmsqst2*(pow2(Dmsqst2)*pow2(Msq) +
        pow3(Dmsqst2) + Dmsqst2*pow4(Msq) + pow6(Msq)) + 7*pow8(Msq)) - Mst2*(
        60*pow2(Msq)*pow3(Dmsqst2) + 75*pow4(Dmsqst2) + 45*pow2(Dmsqst2)*pow4(
        Msq) + 30*Dmsqst2*pow6(Msq) + 19*pow8(Msq))))/(9.*pow2(Mst1)*pow8(Msq))
        - (pow2(Mst1)*(-33606540*pow2(Dmglst2)*pow8(Msq) + 7*pow2(Mst1)*(61300*
        pow2(Msq)*pow3(Dmsqst2) + 51550*pow4(Dmsqst2) + 71050*pow2(Dmsqst2)*
        pow4(Msq) + 80800*Dmsqst2*pow6(Msq) + 3291029*pow8(Msq)) + 11340*log(
        pow2(Mst2)/pow2(Mst1))*(-5*pow2(Mst1)*(6*pow2(Msq)*pow3(Dmsqst2) + 9*
        pow4(Dmsqst2) + 3*pow2(Dmsqst2)*pow4(Msq) - 26*pow8(Msq)) + 1544*pow2(
        Dmglst2)*pow8(Msq))))/(34020.*pow2(Mst2)*pow8(Msq)))) + (Mt*pow3(s2t)*(
        -4860*Dmglst2*((pow4(Mst1)*(588.7025377229081 - (1438*log(pow2(Mst2)/
        pow2(Mst1)))/9. + (440*Dmsqst2*(pow2(Dmsqst2)*pow2(Msq) + pow3(Dmsqst2)
        + Dmsqst2*pow4(Msq) + pow6(Msq)))/(3.*pow8(Msq))))/pow2(Mst2) + pow2(
        Mst1)*(255.2358024691358 - (332*log(pow2(Mst2)/pow2(Mst1)))/3. + (80*
        Dmsqst2*(5 - 3*log(pow2(Mst2)/pow2(Mst1)))*(pow2(Dmsqst2)*pow2(Msq) +
        pow3(Dmsqst2) + Dmsqst2*pow4(Msq) + pow6(Msq)))/(3.*pow8(Msq))) + pow2(
        Mst2)*(322.1240740740741 - (382*log(pow2(Mst2)/pow2(Mst1)))/3. - (80*
        Dmsqst2*(2 + 3*log(pow2(Mst2)/pow2(Mst1)))*(pow2(Dmsqst2)*pow2(Msq) +
        pow3(Dmsqst2) + Dmsqst2*pow4(Msq) + pow6(Msq)))/(3.*pow8(Msq)))) + (
        pow2(Mst1)*(4100546*pow2(Dmglst2)*pow8(Msq) + pow2(Mst1)*(793800*pow2(
        Msq)*pow3(Dmsqst2) + 1158300*pow4(Dmsqst2) + 429300*pow2(Dmsqst2)*pow4(
        Msq) + 64800*Dmsqst2*pow6(Msq) + 104839*pow8(Msq)) - 1080*log(pow2(
        Mst2)/pow2(Mst1))*(-9*pow2(Mst1)*(50*pow2(Msq)*pow3(Dmsqst2) + 75*pow4(
        Dmsqst2) + 25*pow2(Dmsqst2)*pow4(Msq) - 17*pow8(Msq)) + 10*pow2(
        Dmglst2)*pow8(Msq))))/(Mst2*pow8(Msq)) - (3*Mst2*(375896*pow2(Dmglst2)*
        pow8(Msq) - 360*log(pow2(Mst2)/pow2(Mst1))*(-240*pow2(Msq)*(pow2(Mst1)
        + pow2(Mst2))*pow3(Dmsqst2) - 300*(pow2(Mst1) + pow2(Mst2))*pow4(
        Dmsqst2) - 180*pow2(Dmsqst2)*(pow2(Mst1) + pow2(Mst2))*pow4(Msq) - 120*
        Dmsqst2*(pow2(Mst1) + pow2(Mst2))*pow6(Msq) + 90*pow2(Dmglst2)*pow8(
        Msq) - (234*pow2(Mst1) + 7*pow2(Mst2))*pow8(Msq)) + 5*(43200*pow2(Msq)*
        pow2(Mst2)*pow3(Dmsqst2) - 4320*(pow2(Mst1) - 13*pow2(Mst2))*pow4(
        Dmsqst2) + 4320*pow2(Dmsqst2)*(pow2(Mst1) + 7*pow2(Mst2))*pow4(Msq) +
        8640*Dmsqst2*(pow2(Mst1) + 2*pow2(Mst2))*pow6(Msq) + (-9346*pow2(Mst1)
        + 20421*pow2(Mst2))*pow8(Msq))))/pow8(Msq)))/4860. + s2t*pow3(Mt)*(-(((
        -416*(-192271 + 50085*log(pow2(Mst2)/pow2(Mst1)))*pow2(Dmglst2)*pow2(
        Mst1))/8505. + pow4(Mst1)*(1892.2106995884774 - (4864*log(pow2(Mst2)/
        pow2(Mst1)))/9. + (160*Dmsqst2*(5 - 4*log(pow2(Mst2)/pow2(Mst1)))*(
        pow2(Dmsqst2)*pow2(Msq) + pow3(Dmsqst2) + Dmsqst2*pow4(Msq) + pow6(Msq)
        ))/(3.*pow8(Msq))))/pow3(Mst2)) + (Dmglst2*(4189 - 1448*log(pow2(Mst2)/
        pow2(Mst1)) + 9*((pow4(Mst1)*(10873.159945130315 - (5824*log(pow2(Mst2)
        /pow2(Mst1)))/3. + (160*Dmsqst2*(87 - 44*log(pow2(Mst2)/pow2(Mst1)))*(
        pow2(Dmsqst2)*pow2(Msq) + pow3(Dmsqst2) + Dmsqst2*pow4(Msq) + pow6(Msq)
        ))/(3.*pow8(Msq))))/pow4(Mst2) + (pow2(Mst1)*(3778.8814814814814 - (
        7880*log(pow2(Mst2)/pow2(Mst1)))/9. + (160*Dmsqst2*(51 - 28*log(pow2(
        Mst2)/pow2(Mst1)))*(pow2(Dmsqst2)*pow2(Msq) + pow3(Dmsqst2) + Dmsqst2*
        pow4(Msq) + pow6(Msq)))/(3.*pow8(Msq))))/pow2(Mst2)) + (480*Dmsqst2*(19
        - 12*log(pow2(Mst2)/pow2(Mst1)))*(pow2(Dmsqst2)*pow2(Msq) + pow3(
        Dmsqst2) + Dmsqst2*pow4(Msq) + pow6(Msq)))/pow8(Msq)))/9. + (-3980584*
        pow2(Dmglst2)*pow8(Msq) + 2520*log(pow2(Mst2)/pow2(Mst1))*(240*pow2(
        Msq)*(3*pow2(Mst1) + 2*pow2(Mst2))*pow3(Dmsqst2) + 120*(8*pow2(Mst1) +
        5*pow2(Mst2))*pow4(Dmsqst2) + 120*pow2(Dmsqst2)*(4*pow2(Mst1) + 3*pow2(
        Mst2))*pow4(Msq) + 240*Dmsqst2*(pow2(Mst1) + pow2(Mst2))*pow6(Msq) +
        740*pow2(Dmglst2)*pow8(Msq) + (495*pow2(Mst1) + 379*pow2(Mst2))*pow8(
        Msq)) - 7*(21600*pow2(Msq)*(19*pow2(Mst1) + 6*pow2(Mst2))*pow3(Dmsqst2)
        + 10800*(52*pow2(Mst1) + 17*pow2(Mst2))*pow4(Dmsqst2) + 10800*pow2(
        Dmsqst2)*(24*pow2(Mst1) + 7*pow2(Mst2))*pow4(Msq) + 21600*Dmsqst2*(5*
        pow2(Mst1) + pow2(Mst2))*pow6(Msq) + (470657*pow2(Mst1) + 206079*pow2(
        Mst2))*pow8(Msq)))/(2835.*Mst2*pow8(Msq))) - pow4(s2t)*(Dmglst2*pow2(
        Mst1)*((Dmglst2*(36721 - 2943*log(pow2(Mst2)/pow2(Mst1))))/486. - Mst2*
        (43.82962962962963 + (13*log(pow2(Mst2)/pow2(Mst1)))/9. - (20*Dmsqst2*(
        -4 + log(pow2(Mst2)/pow2(Mst1)))*(pow2(Dmsqst2)*pow2(Msq) + pow3(
        Dmsqst2) + Dmsqst2*pow4(Msq) + pow6(Msq)))/(3.*pow8(Msq)))) + (pow4(
        Mst1)*(149867 + 4860*log(pow2(Mst2)/pow2(Mst1)) + (40*Dmsqst2*(1193 +
        324*log(pow2(Mst2)/pow2(Mst1))))/pow2(Msq) + (1620*(17 + 14*log(pow2(
        Mst2)/pow2(Mst1)))*pow2(Dmsqst2))/pow4(Msq) + (80*(92 + 405*log(pow2(
        Mst2)/pow2(Mst1)))*pow3(Dmsqst2))/pow6(Msq) - (3888*Dmglst2*(
        155.89787379972566 - 8*log(pow2(Mst2)/pow2(Mst1)) + (5*Dmsqst2*(7 + 10*
        log(pow2(Mst2)/pow2(Mst1)))*(pow2(Dmsqst2)*pow2(Msq) + pow3(Dmsqst2) +
        Dmsqst2*pow4(Msq) + pow6(Msq)))/(3.*pow8(Msq))))/Mst2 + (20*(-641 +
        2106*log(pow2(Mst2)/pow2(Mst1)))*pow4(Dmsqst2))/pow8(Msq)))/3888. -
        Dmglst2*pow3(Mst2)*(18.212962962962962 - (140*log(pow2(Mst2)/pow2(Mst1)
        ))/9. - (10*Dmsqst2*(4 + 3*log(pow2(Mst2)/pow2(Mst1)))*(pow2(Dmsqst2)*
        pow2(Msq) + pow3(Dmsqst2) + Dmsqst2*pow4(Msq) + pow6(Msq)))/(3.*pow8(
        Msq))) + (pow4(Mst2)*(-(pow2(Mst1)*(1150*pow2(Msq)*pow3(Dmsqst2) +
        2880*pow4(Dmsqst2) - 580*pow2(Dmsqst2)*pow4(Msq) - 2310*Dmsqst2*pow6(
        Msq) + 12*log(pow2(Mst2)/pow2(Mst1))*(30*pow2(Msq)*pow3(Dmsqst2) + 45*
        pow4(Dmsqst2) + 15*pow2(Dmsqst2)*pow4(Msq) - 92*pow8(Msq)) - 6053*pow8(
        Msq))) - 624*pow2(Dmglst2)*pow8(Msq) + 12*Mst2*(4*Dmglst2*(15*Dmsqst2*(
        pow2(Dmsqst2)*pow2(Msq) + pow3(Dmsqst2) + Dmsqst2*pow4(Msq) + pow6(Msq)
        ) + 7*pow8(Msq)) - Mst2*(60*pow2(Msq)*pow3(Dmsqst2) + 75*pow4(Dmsqst2)
        + 45*pow2(Dmsqst2)*pow4(Msq) + 30*Dmsqst2*pow6(Msq) + 19*pow8(Msq)))))/
        (216.*pow2(Mst1)*pow8(Msq)) - (pow2(Mst2)*(21281*pow2(Dmglst2)*pow8(
        Msq) - pow2(Mst1)*(63050*pow2(Msq)*pow3(Dmsqst2) + 38850*pow4(Dmsqst2)
        + 87250*pow2(Dmsqst2)*pow4(Msq) + 111450*Dmsqst2*pow6(Msq) + 362299*
        pow8(Msq)) + 180*log(pow2(Mst2)/pow2(Mst1))*(206*pow2(Dmglst2)*pow8(
        Msq) + 5*pow2(Mst1)*(24*pow2(Msq)*pow3(Dmsqst2) + 30*pow4(Dmsqst2) +
        18*pow2(Dmsqst2)*pow4(Msq) + 12*Dmsqst2*pow6(Msq) + 37*pow8(Msq)))))/(
        3240.*pow8(Msq))) + pow4(Mt)*(362.837037037037 - (596*log(pow2(Mst2)/
        pow2(Mst1)))/3. - (160*Dmsqst2*(-3 + 2*log(pow2(Mst2)/pow2(Mst1))))/(3.
        *pow2(Msq)) + (40*(25 - 14*log(pow2(Mst2)/pow2(Mst1)))*pow2(Dmsqst2))/(
        3.*pow4(Msq)) + (80*(19 - 10*log(pow2(Mst2)/pow2(Mst1)))*pow3(Dmsqst2))
        /(3.*pow6(Msq)) - Dmglst2*(((12010.338702723888 - (12928*log(pow2(Mst2)
        /pow2(Mst1)))/3.)*pow4(Mst1))/pow5(Mst2) + (32*Mst2*(7 + (15*Dmsqst2*(
        pow2(Dmsqst2)*pow2(Msq) + pow3(Dmsqst2) + Dmsqst2*pow4(Msq) + pow6(Msq)
        ))/pow8(Msq)))/(9.*pow2(Mst1)) + (pow2(Mst1)*(4492.976366843033 - 1632*
        log(pow2(Mst2)/pow2(Mst1)) + (160*Dmsqst2*(77 - 36*log(pow2(Mst2)/pow2(
        Mst1)))*(pow2(Dmsqst2)*pow2(Msq) + pow3(Dmsqst2) + Dmsqst2*pow4(Msq) +
        pow6(Msq)))/(3.*pow8(Msq))))/pow3(Mst2) + (886.0994708994709 - (3136*
        log(pow2(Mst2)/pow2(Mst1)))/9. + (160*Dmsqst2*(23 - 10*log(pow2(Mst2)/
        pow2(Mst1)))*(pow2(Dmsqst2)*pow2(Msq) + pow3(Dmsqst2) + Dmsqst2*pow4(
        Msq) + pow6(Msq)))/(3.*pow8(Msq)))/Mst2) + ((16868.94226925338 - (
        38672*log(pow2(Mst2)/pow2(Mst1)))/9.)*pow2(Dmglst2)*pow2(Mst1) + pow4(
        Mst1)*(1031.8831275720165 - (1952*log(pow2(Mst2)/pow2(Mst1)))/3. + (
        160*Dmsqst2*(7 - 6*log(pow2(Mst2)/pow2(Mst1)))*(pow2(Dmsqst2)*pow2(Msq)
        + pow3(Dmsqst2) + Dmsqst2*pow4(Msq) + pow6(Msq)))/(3.*pow8(Msq))))/
        pow4(Mst2) + (40*(51 - 26*log(pow2(Mst2)/pow2(Mst1)))*pow4(Dmsqst2))/(
        3.*pow8(Msq)) + (8*(52*pow2(Dmglst2)*pow8(Msq) + pow2(Mst2)*(60*pow2(
        Msq)*pow3(Dmsqst2) + 75*pow4(Dmsqst2) + 45*pow2(Dmsqst2)*pow4(Msq) +
        30*Dmsqst2*pow6(Msq) + 19*pow8(Msq))))/(9.*pow2(Mst1)*pow8(Msq)) + (
        7145462*pow2(Dmglst2)*pow8(Msq) + 7*pow2(Mst1)*(615600*pow2(Msq)*pow3(
        Dmsqst2) + 858600*pow4(Dmsqst2) + 372600*pow2(Dmsqst2)*pow4(Msq) +
        129600*Dmsqst2*pow6(Msq) + 310039*pow8(Msq)) - 2520*log(pow2(Mst2)/
        pow2(Mst1))*(724*pow2(Dmglst2)*pow8(Msq) + pow2(Mst1)*(840*pow2(Msq)*
        pow3(Dmsqst2) + 1140*pow4(Dmsqst2) + 540*pow2(Dmsqst2)*pow4(Msq) + 240*
        Dmsqst2*pow6(Msq) + 403*pow8(Msq))))/(2835.*pow2(Mst2)*pow8(Msq)))) -
        pow2(Mt)*((40*(3*Dmsqst2 + pow2(Msq))*pow2(Mst2)*pow3(Dmsqst2))/(3.*
        pow8(Msq)) + (pow2(s2t)*(18522000*pow3(log(pow2(Mst2)/pow2(Mst1)))*
        pow4(Mst1)*(583*pow2(Mst1)*pow3(Mst2) + 6*pow2(Dmglst2)*(95*Mst2*pow2(
        Mst1) + 32*pow3(Mst2)) + 337*Mst2*pow4(Mst1) + 4*Dmglst2*(289*pow2(
        Mst1)*pow2(Mst2) + 362*pow4(Mst1) + 96*pow4(Mst2)) + 21*pow5(Mst2))*
        pow8(Msq) + 1764*Mst2*pow2(Dmglst2)*(75*(-2164753 + 3360*B4 - 3360*D3 +
        1680*DN - 7840*OepS2 + 169452*S2)*pow2(Mst2)*pow4(Mst1) - 22554000*
        pow2(Mst1)*pow4(Mst2) + (-68526791 - 2772000*OepS2 + 107403300*S2)*
        pow6(Mst1) + 2016000*pow6(Mst2))*pow8(Msq) - 12*Dmglst2*pow2(Mst1)*(
        3858750*pow2(Msq)*pow3(Dmsqst2)*(-1356*pow2(Mst2)*pow4(Mst1) - 552*
        pow2(Mst1)*pow4(Mst2) + 1931*pow6(Mst1) - 432*pow6(Mst2)) + 3858750*
        pow4(Dmsqst2)*(-1356*pow2(Mst2)*pow4(Mst1) - 552*pow2(Mst1)*pow4(Mst2)
        + 1931*pow6(Mst1) - 432*pow6(Mst2)) + 3858750*pow2(Dmsqst2)*pow4(Msq)*(
        -1356*pow2(Mst2)*pow4(Mst1) - 552*pow2(Mst1)*pow4(Mst2) + 1931*pow6(
        Mst1) - 432*pow6(Mst2)) + 3858750*Dmsqst2*pow6(Msq)*(-1356*pow2(Mst2)*
        pow4(Mst1) - 552*pow2(Mst1)*pow4(Mst2) + 1931*pow6(Mst1) - 432*pow6(
        Mst2)) + (5488*(608404 + 75000*OepS2 - 188325*S2)*pow2(Mst2)*pow4(Mst1)
        - 617400*(-11674 + 120*B4 - 120*D3 + 60*DN + 17091*S2)*pow2(Mst1)*pow4(
        Mst2) + (-3650009897 + 2110136000*OepS2 - 52385772600*S2)*pow6(Mst1) +
        481572000*pow6(Mst2))*pow8(Msq)) - Mst2*(-20*Dmsqst2*pow2(Mst1)*pow6(
        Msq)*(6174*(36449 + 4000*OepS2 - 310500*S2)*pow2(Mst2)*pow4(Mst1) +
        1157625*(-1173 + 32*OepS2 - 2916*S2)*pow2(Mst1)*pow4(Mst2) + 2*(-
        22430251 + 5488000*OepS2 - 324135000*S2)*pow6(Mst1) + 83349000*pow6(
        Mst2)) + 10*pow2(Dmsqst2)*pow2(Mst1)*pow4(Msq)*(-41160*(-1241 + 200*
        OepS2 - 45225*S2)*pow2(Mst2)*pow4(Mst1) - 771750*(-3439 + 64*OepS2 -
        8424*S2)*pow2(Mst1)*pow4(Mst2) + (-1351886621 + 19208000*OepS2 -
        639009000*S2)*pow6(Mst1) + 166698000*pow6(Mst2)) + 20*pow2(Msq)*pow2(
        Mst1)*pow3(Dmsqst2)*(2058*(134167 + 8000*OepS2 - 27000*S2)*pow2(Mst2)*
        pow4(Mst1) - 385875*(-3359 + 32*OepS2 - 8100*S2)*pow2(Mst1)*pow4(Mst2)
        + (-1396747123 + 30184000*OepS2 - 1287279000*S2)*pow6(Mst1) +
        250047000*pow6(Mst2)) + 70*pow4(Dmsqst2)*(330750*(1093 + 2592*S2)*pow4(
        Mst1)*pow4(Mst2) + 10584*(14218 + 1000*OepS2 - 28125*S2)*pow2(Mst2)*
        pow6(Mst1) + 119070000*pow2(Mst1)*pow6(Mst2) + (-605014553 + 14504000*
        OepS2 - 644301000*S2)*pow8(Mst1)) + pow8(Msq)*(77175*(1013263 + 2880*D3
        - 2880*DN - 25440*OepS2 - 137052*S2)*pow4(Mst1)*pow4(Mst2) - 10290*(-
        2367149 + 21600*D3 - 21600*DN + 503200*OepS2 - 14889420*S2)*pow2(Mst2)*
        pow6(Mst1) - 4834242000*pow2(Mst1)*pow6(Mst2) + (50004748544 -
        9407804000*OepS2 + 369631514700*S2)*pow8(Mst1) + 889056000*pow8(Mst2)))
        + 264600*pow2(log(pow2(Mst2)/pow2(Mst1)))*(42*Mst2*pow2(Dmglst2)*(645*
        pow2(Mst2)*pow4(Mst1) + 775*pow2(Mst1)*pow4(Mst2) + 5574*pow6(Mst1) -
        480*pow6(Mst2))*pow8(Msq) - 4*Dmglst2*pow8(Msq)*(31605*pow4(Mst1)*pow4(
        Mst2) + 11368*pow2(Mst2)*pow6(Mst1) - 9555*pow2(Mst1)*pow6(Mst2) +
        52512*pow8(Mst1) + 3360*pow8(Mst2)) - Mst2*(90*pow2(Msq)*(53*pow2(Mst1)
        + 14*pow2(Mst2))*pow3(Dmsqst2)*pow6(Mst1) + 90*pow2(Dmsqst2)*(8*pow2(
        Mst1) - 35*pow2(Mst2))*pow4(Msq)*pow6(Mst1) - 90*Dmsqst2*(37*pow2(Mst1)
        + 84*pow2(Mst2))*pow6(Msq)*pow6(Mst1) + 630*pow4(Dmsqst2)*(9*pow2(Mst2)
        *pow6(Mst1) + 14*pow8(Mst1)) + pow8(Msq)*(55860*pow4(Mst1)*pow4(Mst2) +
        67060*pow2(Mst2)*pow6(Mst1) - 15750*pow2(Mst1)*pow6(Mst2) + 73657*pow8(
        Mst1) + 3360*pow8(Mst2)))) + 630*log(pow2(Mst2)/pow2(Mst1))*(1176*Mst2*
        pow2(Dmglst2)*(25*(6485 + 1134*S2)*pow2(Mst2)*pow4(Mst1) + 7650*pow2(
        Mst1)*pow4(Mst2) + (102713 + 133650*S2)*pow6(Mst1) - 2400*pow6(Mst2))*
        pow8(Msq) + Mst2*(-80*Dmsqst2*pow2(Mst1)*pow6(Msq)*(441*(-517 + 450*S2)
        *pow2(Mst2)*pow4(Mst1) + 11025*(16 + 27*S2)*pow2(Mst1)*pow4(Mst2) + 2*(
        -58853 + 44100*S2)*pow6(Mst1) - 66150*pow6(Mst2)) + 10*pow2(Msq)*pow2(
        Mst1)*pow3(Dmsqst2)*(4704*(166 + 225*S2)*pow2(Mst2)*pow4(Mst1) -
        264600*(4 + 3*S2)*pow2(Mst1)*pow4(Mst2) + (-623527 + 1940400*S2)*pow6(
        Mst1) + 529200*pow6(Mst2)) + 5*pow2(Dmsqst2)*pow2(Mst1)*pow4(Msq)*(-
        5880*(-443 + 90*S2)*pow2(Mst2)*pow4(Mst1) - 352800*(7 + 9*S2)*pow2(
        Mst1)*pow4(Mst2) + (318121 + 1234800*S2)*pow6(Mst1) + 1058400*pow6(
        Mst2)) + 35*pow4(Dmsqst2)*(-252000*pow4(Mst1)*pow4(Mst2) + 1512*(49 +
        450*S2)*pow2(Mst2)*pow6(Mst1) + 151200*pow2(Mst1)*pow6(Mst2) + (-401747
        + 932400*S2)*pow8(Mst1)) - 2*pow8(Msq)*(7350*(-4534 + 4293*S2)*pow4(
        Mst1)*pow4(Mst2) + 980*(-52322 + 84915*S2)*pow2(Mst2)*pow6(Mst1) +
        6791400*pow2(Mst1)*pow6(Mst2) + 3*(-19407863 + 50398950*S2)*pow8(Mst1)
        - 1411200*pow8(Mst2))) + 4*Dmglst2*(55125*pow2(Msq)*(37*pow2(Mst1) -
        24*pow2(Mst2))*pow3(Dmsqst2)*pow6(Mst1) + 55125*pow2(Dmsqst2)*(37*pow2(
        Mst1) - 24*pow2(Mst2))*pow4(Msq)*pow6(Mst1) + 55125*Dmsqst2*(37*pow2(
        Mst1) - 24*pow2(Mst2))*pow6(Msq)*pow6(Mst1) + 55125*pow4(Dmsqst2)*(-24*
        pow2(Mst2)*pow6(Mst1) + 37*pow8(Mst1)) + 2*pow8(Msq)*(25188450*pow4(
        Mst1)*pow4(Mst2) + 98*(77657 + 202500*S2)*pow2(Mst2)*pow6(Mst1) -
        2954700*pow2(Mst1)*pow6(Mst2) + 3*(-4261807 + 33912900*S2)*pow8(Mst1) +
        705600*pow8(Mst2))))))/(2.50047e8*pow3(Mst2)*pow4(Mst1)*pow8(Msq)))))/
        pow4(Mt)/pow2(Sbeta)*12.;

   return result;
}

/**
 * 	@return returns the susy log^1 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double himalaya::H6bq22g::calc_coef_at_as2_no_sm_logs_log1(){

   const double result =
      ((pow2(Sbeta)*(4*Mst2*pow2(Dmglst2)*(-225*pow2(Mst1)*pow4(Mst2)*(-664*
        pow2(Mst2)*pow2(Mt)*pow2(s2t) - 2048*Mst2*s2t*pow3(Mt) + 512*Mt*pow3(
        Mst2)*pow3(s2t) + 1840*pow4(Mt) + 83*pow4(Mst2)*pow4(s2t)) + 2*pow2(
        Mst2)*pow4(Mst1)*(-150*(-7121 + 6336*z2)*pow2(Mst2)*pow2(Mt)*pow2(s2t)
        - 400*Mst2*s2t*(-1934 + 1791*z2)*pow3(Mt) - 300*Mt*(-875 + 648*z2)*
        pow3(Mst2)*pow3(s2t) + 8*(71687 + 76050*z2)*pow4(Mt) + 75*(-205 + 129*
        z2)*pow4(Mst2)*pow4(s2t)) + 25*(12*(7903 - 5760*z2)*pow2(Mst2)*pow2(Mt)
        *pow2(s2t) - 96*Mst2*s2t*(-1055 + 2181*z2)*pow3(Mt) + 13992*Mt*pow3(
        Mst2)*pow3(s2t) + 80*(-605 + 4608*z2)*pow4(Mt) + 3*(409 - 258*z2)*pow4(
        Mst2)*pow4(s2t))*pow6(Mst1) + 3600*pow2(-4*pow2(Mt) + pow2(Mst2)*pow2(
        s2t))*pow6(Mst2))*pow8(Msq) + 450*pow2(log(pow2(Mst2)/pow2(Mst1)))*
        pow4(Mst1)*(-8*Dmglst2*(8*Mt*(-157*Mst2*s2t*pow2(Mt) - 108*Mt*pow2(
        Mst2)*pow2(s2t) + 398*pow3(Mt) - 16*pow3(Mst2)*pow3(s2t))*pow4(Mst1) +
        pow4(Mst2)*(-320*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 848*Mst2*s2t*pow3(Mt)
        - 532*Mt*pow3(Mst2)*pow3(s2t) + 512*pow4(Mt) - 99*pow4(Mst2)*pow4(s2t))
        + pow2(Mst1)*pow2(Mst2)*(-768*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 448*Mst2*
        s2t*pow3(Mt) - 212*Mt*pow3(Mst2)*pow3(s2t) + 1896*pow4(Mt) + 35*pow4(
        Mst2)*pow4(s2t))) + Mst2*(8*pow2(Mst1)*pow2(Mst2)*(425*pow2(Mst2)*pow2(
        Mt)*pow2(s2t) - 704*Mst2*s2t*pow3(Mt) + 212*Mt*pow3(Mst2)*pow3(s2t) -
        76*pow4(Mt) + 3*pow4(Mst2)*pow4(s2t)) + pow4(Mst1)*(2688*pow2(Mst2)*
        pow2(Mt)*pow2(s2t) - 5056*Mst2*s2t*pow3(Mt) + 1024*Mt*pow3(Mst2)*pow3(
        s2t) - 2848*pow4(Mt) + 41*pow4(Mst2)*pow4(s2t)) + pow4(Mst2)*(296*pow2(
        Mst2)*pow2(Mt)*pow2(s2t) - 6784*Mst2*s2t*pow3(Mt) + 2208*Mt*pow3(Mst2)*
        pow3(s2t) - 672*pow4(Mt) + 519*pow4(Mst2)*pow4(s2t))) + 4*Mst2*pow2(
        Dmglst2)*(320*pow2(Mt)*pow2(s2t)*pow4(Mst2) - 512*pow2(Mst2)*pow4(Mt) +
        pow2(Mst1)*(768*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 1152*Mst2*s2t*pow3(Mt)
        + 3640*pow4(Mt) - 35*pow4(Mst2)*pow4(s2t)) + 768*Mt*pow3(s2t)*pow5(
        Mst2) + 99*pow4(s2t)*pow6(Mst2)))*pow8(Msq) - 100*Dmglst2*(-240*Mst2*
        pow2(Msq)*pow3(Dmsqst2)*pow3(Mt)*pow4(Mst1)*(36*Mst2*(Mt*(33 - 18*z2) +
        2*Mst2*s2t*(-12 + 7*z2))*pow2(Mst1) + (Mt*(283 - 180*z2) + 12*Mst2*s2t*
        (-29 + 18*z2))*pow3(Mst2) + 18*s2t*(-83 + 44*z2)*pow4(Mst1)) - 240*
        Mst2*pow3(Mt)*pow4(Dmsqst2)*pow4(Mst1)*(36*Mst2*(Mt*(33 - 18*z2) + 2*
        Mst2*s2t*(-12 + 7*z2))*pow2(Mst1) + (Mt*(283 - 180*z2) + 12*Mst2*s2t*(-
        29 + 18*z2))*pow3(Mst2) + 18*s2t*(-83 + 44*z2)*pow4(Mst1)) - 240*Mst2*
        pow2(Dmsqst2)*pow3(Mt)*pow4(Msq)*pow4(Mst1)*(36*Mst2*(Mt*(33 - 18*z2) +
        2*Mst2*s2t*(-12 + 7*z2))*pow2(Mst1) + (Mt*(283 - 180*z2) + 12*Mst2*s2t*
        (-29 + 18*z2))*pow3(Mst2) + 18*s2t*(-83 + 44*z2)*pow4(Mst1)) - 240*
        Dmsqst2*Mst2*pow3(Mt)*pow4(Mst1)*(36*Mst2*(Mt*(33 - 18*z2) + 2*Mst2*
        s2t*(-12 + 7*z2))*pow2(Mst1) + (Mt*(283 - 180*z2) + 12*Mst2*s2t*(-29 +
        18*z2))*pow3(Mst2) + 18*s2t*(-83 + 44*z2)*pow4(Mst1))*pow6(Msq) + pow8(
        Msq)*(2*pow4(Mst1)*pow4(Mst2)*(12*(-3565 + 1728*z2)*pow2(Mst2)*pow2(Mt)
        *pow2(s2t) + 8*Mst2*s2t*(-3659 + 252*z2)*pow3(Mt) + 72*Mt*(-192 + 91*
        z2)*pow3(Mst2)*pow3(s2t) + 44*(427 + 432*z2)*pow4(Mt) + 9*(211 - 86*z2)
        *pow4(Mst2)*pow4(s2t)) + 6*pow2(Mst2)*(4*(-1955 + 1152*z2)*pow2(Mst2)*
        pow2(Mt)*pow2(s2t) - 48*Mst2*s2t*(-13 + 262*z2)*pow3(Mt) + 24*Mt*(111 -
        17*z2)*pow3(Mst2)*pow3(s2t) + 8*(-1471 + 3240*z2)*pow4(Mt) + 3*(-1 +
        86*z2)*pow4(Mst2)*pow4(s2t))*pow6(Mst1) - 18*pow2(Mst1)*(-536*pow2(
        Mst2)*pow2(Mt)*pow2(s2t) - 512*Mst2*s2t*pow3(Mt) + 128*Mt*pow3(Mst2)*
        pow3(s2t) + 560*pow4(Mt) + 131*pow4(Mst2)*pow4(s2t))*pow6(Mst2) + 3*(
        192*(-77 + 48*z2)*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 48*Mst2*s2t*(-867 +
        1012*z2)*pow3(Mt) - 276*Mt*pow3(Mst2)*pow3(s2t) + 16*(-6445 + 6768*z2)*
        pow4(Mt) - 441*pow4(Mst2)*pow4(s2t))*pow8(Mst1) - 2304*pow2(Mt)*(-2*
        pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow8(Mst2))) - 25*Mst2*(-120*pow2(Msq)
        *pow2(Mst1)*pow3(Dmsqst2)*(9*pow2(Mst2)*pow4(Mst1)*(8*pow2(Mst2)*pow2(
        Mt)*pow2(s2t) - 64*Mst2*s2t*(-5 + 3*z2)*pow3(Mt) + 16*(-23 + 14*z2)*
        pow4(Mt) + pow4(Mst2)*pow4(s2t)) + pow2(Mst1)*pow4(Mst2)*(-144*pow2(
        Mst2)*pow2(Mt)*pow2(s2t) - 144*Mst2*s2t*(-13 + 8*z2)*pow3(Mt) + 4*(-353
        + 180*z2)*pow4(Mt) + 9*pow4(Mst2)*pow4(s2t)) + 9*(-16*Mst2*s2t*(-9 + 4*
        z2)*pow3(Mt) + 8*(-31 + 12*z2)*pow4(Mt) - pow4(Mst2)*pow4(s2t))*pow6(
        Mst1) - 9*pow2(-4*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow6(Mst2)) - 120*
        pow2(Dmsqst2)*pow2(Mst1)*pow4(Msq)*(9*pow2(Mst2)*pow4(Mst1)*(8*pow2(
        Mst2)*pow2(Mt)*pow2(s2t) - 32*Mst2*s2t*(-7 + 4*z2)*pow3(Mt) + 16*(-16 +
        9*z2)*pow4(Mt) + pow4(Mst2)*pow4(s2t)) + pow2(Mst1)*pow4(Mst2)*(-144*
        pow2(Mst2)*pow2(Mt)*pow2(s2t) - 72*Mst2*s2t*(-19 + 12*z2)*pow3(Mt) + 2*
        (-545 + 252*z2)*pow4(Mt) + 9*pow4(Mst2)*pow4(s2t)) + 9*(-16*Mst2*s2t*(-
        9 + 4*z2)*pow3(Mt) + 8*(-31 + 12*z2)*pow4(Mt) - pow4(Mst2)*pow4(s2t))*
        pow6(Mst1) - 9*pow2(-4*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow6(Mst2)) -
        360*Dmsqst2*pow2(Mst1)*pow6(Msq)*(3*pow2(Mst2)*pow4(Mst1)*(8*pow2(Mst2)
        *pow2(Mt)*pow2(s2t) - 64*Mst2*s2t*(-2 + z2)*pow3(Mt) + 16*(-9 + 4*z2)*
        pow4(Mt) + pow4(Mst2)*pow4(s2t)) + pow2(Mst1)*pow4(Mst2)*(-48*pow2(
        Mst2)*pow2(Mt)*pow2(s2t) - 96*Mst2*s2t*(-3 + 2*z2)*pow3(Mt) + 32*(-8 +
        3*z2)*pow4(Mt) + 3*pow4(Mst2)*pow4(s2t)) + 3*(-16*Mst2*s2t*(-9 + 4*z2)*
        pow3(Mt) + 8*(-31 + 12*z2)*pow4(Mt) - pow4(Mst2)*pow4(s2t))*pow6(Mst1)
        - 3*pow2(-4*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow6(Mst2)) - 360*pow2(
        Mst1)*pow4(Dmsqst2)*(3*pow2(Mst2)*pow4(Mst1)*(8*pow2(Mst2)*pow2(Mt)*
        pow2(s2t) - 32*Mst2*s2t*(-13 + 8*z2)*pow3(Mt) + 16*(-30 + 19*z2)*pow4(
        Mt) + pow4(Mst2)*pow4(s2t)) + pow2(Mst1)*pow4(Mst2)*(-48*pow2(Mst2)*
        pow2(Mt)*pow2(s2t) + 24*Mst2*s2t*(33 - 20*z2)*pow3(Mt) + (-578 + 312*
        z2)*pow4(Mt) + 3*pow4(Mst2)*pow4(s2t)) + 3*(-16*Mst2*s2t*(-9 + 4*z2)*
        pow3(Mt) + 8*(-31 + 12*z2)*pow4(Mt) - pow4(Mst2)*pow4(s2t))*pow6(Mst1)
        - 3*pow2(-4*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow6(Mst2)) - pow8(Msq)*(
        4*pow4(Mst1)*pow4(Mst2)*(-72*(-445 + 96*z2)*pow2(Mst2)*pow2(Mt)*pow2(
        s2t) - 144*Mst2*s2t*(-121 + 296*z2)*pow3(Mt) - 144*Mt*(-151 + 19*z2)*
        pow3(Mst2)*pow3(s2t) + 8*(-1993 + 2151*z2 + 216*z3)*pow4(Mt) + 81*(8 +
        5*z2)*pow4(Mst2)*pow4(s2t)) - 36*pow2(Mst2)*(-696*pow2(Mst2)*pow2(Mt)*
        pow2(s2t) + 32*Mst2*s2t*(-179 + 128*z2)*pow3(Mt) + 16*Mt*(148 - 17*z2)*
        pow3(Mst2)*pow3(s2t) - 8*(-403 + 300*z2)*pow4(Mt) + (551 + 86*z2)*pow4(
        Mst2)*pow4(s2t))*pow6(Mst1) + 36*pow2(Mst1)*(-744*pow2(Mst2)*pow2(Mt)*
        pow2(s2t) - 1024*Mst2*s2t*pow3(Mt) + 256*Mt*pow3(Mst2)*pow3(s2t) +
        1232*pow4(Mt) + 141*pow4(Mst2)*pow4(s2t))*pow6(Mst2) + 9*(864*pow2(
        Mst2)*pow2(Mt)*pow2(s2t) - 512*Mst2*s2t*(-49 + 31*z2)*pow3(Mt) - 784*
        Mt*pow3(Mst2)*pow3(s2t) + 8*(-2503 + 1408*z2)*pow4(Mt) + (1547 + 164*
        z2)*pow4(Mst2)*pow4(s2t))*pow8(Mst1) + 4608*pow2(Mt)*(-2*pow2(Mt) +
        pow2(Mst2)*pow2(s2t))*pow8(Mst2))) - 30*log(pow2(Mst2)/pow2(Mst1))*(2*
        Mst2*pow2(Dmglst2)*(2*pow2(Mst2)*pow4(Mst1)*(3240*pow2(Mst2)*pow2(Mt)*
        pow2(s2t) + 45920*Mst2*s2t*pow3(Mt) - 9720*Mt*pow3(Mst2)*pow3(s2t) +
        37408*pow4(Mt) - 4635*pow4(Mst2)*pow4(s2t)) + 45*pow2(Mst1)*pow4(Mst2)*
        (-456*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 1024*Mst2*s2t*pow3(Mt) + 256*Mt*
        pow3(Mst2)*pow3(s2t) + 400*pow4(Mt) + 153*pow4(Mst2)*pow4(s2t)) + 2*(-
        34080*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 147560*Mst2*s2t*pow3(Mt) - 5880*
        Mt*pow3(Mst2)*pow3(s2t) + 415932*pow4(Mt) + 1005*pow4(Mst2)*pow4(s2t))*
        pow6(Mst1) - 1440*pow2(-4*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow6(Mst2))*
        pow8(Msq) - 20*Dmglst2*(1440*Mst2*pow2(Msq)*(-(Mst2*Mt) + 3*s2t*pow2(
        Mst1))*pow3(Dmsqst2)*pow3(Mt)*pow6(Mst1) + 1440*Mst2*(-(Mst2*Mt) + 3*
        s2t*pow2(Mst1))*pow3(Mt)*pow4(Dmsqst2)*pow6(Mst1) + 1440*Mst2*pow2(
        Dmsqst2)*(-(Mst2*Mt) + 3*s2t*pow2(Mst1))*pow3(Mt)*pow4(Msq)*pow6(Mst1)
        + 1440*Dmsqst2*Mst2*(-(Mst2*Mt) + 3*s2t*pow2(Mst1))*pow3(Mt)*pow6(Msq)*
        pow6(Mst1) + pow8(Msq)*(-4*pow4(Mst1)*pow4(Mst2)*(1860*pow2(Mst2)*pow2(
        Mt)*pow2(s2t) + 1472*Mst2*s2t*pow3(Mt) + 774*Mt*pow3(Mst2)*pow3(s2t) -
        2464*pow4(Mt) + 15*pow4(Mst2)*pow4(s2t)) + 2*pow2(Mst2)*(-432*pow2(
        Mst2)*pow2(Mt)*pow2(s2t) - 11656*Mst2*s2t*pow3(Mt) + 996*Mt*pow3(Mst2)*
        pow3(s2t) + 18748*pow4(Mt) + 207*pow4(Mst2)*pow4(s2t))*pow6(Mst1) - 3*
        pow2(Mst1)*(-856*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 1536*Mst2*s2t*pow3(Mt)
        + 384*Mt*pow3(Mst2)*pow3(s2t) + 1200*pow4(Mt) + 203*pow4(Mst2)*pow4(
        s2t))*pow6(Mst2) + 4*(-738*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 10004*Mst2*
        s2t*pow3(Mt) + 135*Mt*pow3(Mst2)*pow3(s2t) + 17532*pow4(Mt) + 12*pow4(
        Mst2)*pow4(s2t))*pow8(Mst1) - 768*pow2(Mt)*(-2*pow2(Mt) + pow2(Mst2)*
        pow2(s2t))*pow8(Mst2))) + 5*Mst2*(pow8(Msq)*(4*pow4(Mst1)*pow4(Mst2)*(
        3372*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 9024*Mst2*s2t*pow3(Mt) + 3912*Mt*
        pow3(Mst2)*pow3(s2t) + 6112*pow4(Mt) + 579*pow4(Mst2)*pow4(s2t)) - 24*
        pow2(Mst2)*(-534*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 1824*Mst2*s2t*pow3(Mt)
        - 60*Mt*pow3(Mst2)*pow3(s2t) - 934*pow4(Mt) + 97*pow4(Mst2)*pow4(s2t))*
        pow6(Mst1) + 6*pow2(Mst1)*(-728*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 1024*
        Mst2*s2t*pow3(Mt) + 256*Mt*pow3(Mst2)*pow3(s2t) + 1200*pow4(Mt) + 139*
        pow4(Mst2)*pow4(s2t))*pow6(Mst2) + 3*(4640*pow2(Mst2)*pow2(Mt)*pow2(
        s2t) - 15296*Mst2*s2t*pow3(Mt) + 240*Mt*pow3(Mst2)*pow3(s2t) + 7128*
        pow4(Mt) - 279*pow4(Mst2)*pow4(s2t))*pow8(Mst1) + 768*pow2(Mt)*(-2*
        pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow8(Mst2)) - 360*pow2(Msq)*pow3(
        Dmsqst2)*pow4(Mst1)*(8*(5*Mt - 2*Mst2*s2t)*pow3(Mt)*pow4(Mst1) - pow2(
        Mst1)*pow4(s2t)*pow6(Mst2) + pow4(s2t)*pow8(Mst2)) - 360*pow2(Dmsqst2)*
        pow4(Msq)*pow4(Mst1)*(8*(5*Mt - 2*Mst2*s2t)*pow3(Mt)*pow4(Mst1) - 4*
        pow4(Mst2)*pow4(Mt) + pow2(Mst1)*(4*pow2(Mst2)*pow4(Mt) - pow4(s2t)*
        pow6(Mst2)) + pow4(s2t)*pow8(Mst2)) - 360*Dmsqst2*pow4(Mst1)*pow6(Msq)*
        (8*(5*Mt - 2*Mst2*s2t)*pow3(Mt)*pow4(Mst1) - 8*pow4(Mst2)*pow4(Mt) +
        pow2(Mst1)*(8*pow2(Mst2)*pow4(Mt) - pow4(s2t)*pow6(Mst2)) + pow4(s2t)*
        pow8(Mst2)) - 360*pow4(Dmsqst2)*pow4(Mst1)*(8*(5*Mt - 2*Mst2*s2t)*pow3(
        Mt)*pow4(Mst1) + 4*pow4(Mst2)*pow4(Mt) - pow2(Mst1)*(4*pow2(Mst2)*pow4(
        Mt) + pow4(s2t)*pow6(Mst2)) + pow4(s2t)*pow8(Mst2))))))/(16200.*pow4(
        Mst1)*pow5(Mst2)*pow8(Msq)))/pow4(Mt)/pow2(Sbeta)*12.; 

   return result;
}

/**
 * 	@return returns the susy log^2 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double himalaya::H6bq22g::calc_coef_at_as2_no_sm_logs_log2(){

   const double result =
      ((pow2(Sbeta)*(Mst2*pow2(Dmglst2)*(pow2(Mst2)*pow4(Mst1)*(33360*pow2(Mst2)
        *pow2(Mt)*pow2(s2t) + 52480*Mst2*s2t*pow3(Mt) + 18896*pow4(Mt) - 10005*
        pow4(Mst2)*pow4(s2t)) + 15*pow2(Mst1)*pow4(Mst2)*(-1496*pow2(Mst2)*
        pow2(Mt)*pow2(s2t) - 3072*Mst2*s2t*pow3(Mt) + 768*Mt*pow3(Mst2)*pow3(
        s2t) + 1456*pow4(Mt) + 475*pow4(Mst2)*pow4(s2t)) + 5*(120*pow2(Mst2)*
        pow2(Mt)*pow2(s2t) - 22656*Mst2*s2t*pow3(Mt) - 2304*Mt*pow3(Mst2)*pow3(
        s2t) + 52384*pow4(Mt) + 879*pow4(Mst2)*pow4(s2t))*pow6(Mst1) - 1440*
        pow2(-4*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow6(Mst2))*pow8(Msq) - 30*
        log(pow2(Mst2)/pow2(Mst1))*pow4(Mst1)*(Mst2*(pow2(Mst1)*pow2(Mst2)*(
        1096*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 1280*Mst2*s2t*pow3(Mt) + 328*Mt*
        pow3(Mst2)*pow3(s2t) - 32*pow4(Mt) - 91*pow4(Mst2)*pow4(s2t)) + 4*pow4(
        Mst2)*(14*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 456*Mst2*s2t*pow3(Mt) + 146*
        Mt*pow3(Mst2)*pow3(s2t) - 166*pow4(Mt) + 33*pow4(Mst2)*pow4(s2t)) -
        pow4(Mst1)*(-768*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 1216*Mst2*s2t*pow3(Mt)
        + 400*pow4(Mt) + 41*pow4(Mst2)*pow4(s2t))) - 2*Dmglst2*(32*pow2(Mt)*(-
        57*Mst2*Mt*s2t + 119*pow2(Mt) - 24*pow2(Mst2)*pow2(s2t))*pow4(Mst1) +
        pow4(Mst2)*(-384*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 912*Mst2*s2t*pow3(Mt)
        - 548*Mt*pow3(Mst2)*pow3(s2t) + 512*pow4(Mt) - 91*pow4(Mst2)*pow4(s2t))
        + pow2(Mst1)*pow2(Mst2)*(-768*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 640*Mst2*
        s2t*pow3(Mt) - 164*Mt*pow3(Mst2)*pow3(s2t) + 2016*pow4(Mt) + 91*pow4(
        Mst2)*pow4(s2t))) + Mst2*pow2(Dmglst2)*(384*pow2(Mt)*pow2(s2t)*pow4(
        Mst2) - 512*pow2(Mst2)*pow4(Mt) + pow2(Mst1)*(768*pow2(Mst2)*pow2(Mt)*
        pow2(s2t) - 1280*Mst2*s2t*pow3(Mt) + 4000*pow4(Mt) - 91*pow4(Mst2)*
        pow4(s2t)) + 768*Mt*pow3(s2t)*pow5(Mst2) + 91*pow4(s2t)*pow6(Mst2)))*
        pow8(Msq) - 10*Dmglst2*pow8(Msq)*(pow4(Mst1)*pow4(Mst2)*(-8208*pow2(
        Mst2)*pow2(Mt)*pow2(s2t) - 4640*Mst2*s2t*pow3(Mt) - 1968*Mt*pow3(Mst2)*
        pow3(s2t) + 4816*pow4(Mt) + 849*pow4(Mst2)*pow4(s2t)) - 3*pow2(Mst2)*(-
        472*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 2496*Mst2*s2t*pow3(Mt) - 1040*Mt*
        pow3(Mst2)*pow3(s2t) - 3360*pow4(Mt) + 37*pow4(Mst2)*pow4(s2t))*pow6(
        Mst1) - 3*pow2(Mst1)*(-984*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 1536*Mst2*
        s2t*pow3(Mt) + 384*Mt*pow3(Mst2)*pow3(s2t) + 1456*pow4(Mt) + 219*pow4(
        Mst2)*pow4(s2t))*pow6(Mst2) + 3*(-2240*Mst2*s2t*pow3(Mt) + 3744*pow4(
        Mt) - 59*pow4(Mst2)*pow4(s2t))*pow8(Mst1) - 768*pow2(Mt)*(-2*pow2(Mt) +
        pow2(Mst2)*pow2(s2t))*pow8(Mst2)) + 5*Mst2*(-720*pow4(Dmsqst2)*pow4(
        Mst1)*pow4(Mst2)*pow4(Mt) + 720*pow2(Dmsqst2)*pow4(Msq)*pow4(Mst1)*
        pow4(Mst2)*pow4(Mt) + 1440*Dmsqst2*pow4(Mst1)*pow4(Mst2)*pow4(Mt)*pow6(
        Msq) + pow8(Msq)*(pow4(Mst1)*pow4(Mst2)*(8592*pow2(Mst2)*pow2(Mt)*pow2(
        s2t) - 4800*Mst2*s2t*pow3(Mt) + 3936*Mt*pow3(Mst2)*pow3(s2t) + 4256*
        pow4(Mt) - 69*pow4(Mst2)*pow4(s2t)) - 3*pow2(Mst2)*(600*pow2(Mst2)*
        pow2(Mt)*pow2(s2t) - 128*Mst2*s2t*pow3(Mt) + 1568*Mt*pow3(Mst2)*pow3(
        s2t) + 672*pow4(Mt) + 355*pow4(Mst2)*pow4(s2t))*pow6(Mst1) + 3*pow2(
        Mst1)*(-856*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 1024*Mst2*s2t*pow3(Mt) +
        256*Mt*pow3(Mst2)*pow3(s2t) + 1456*pow4(Mt) + 155*pow4(Mst2)*pow4(s2t))
        *pow6(Mst2) + (384*Mst2*s2t*pow3(Mt) - 2400*pow4(Mt) + 717*pow4(Mst2)*
        pow4(s2t))*pow8(Mst1) + 384*pow2(Mt)*(-2*pow2(Mt) + pow2(Mst2)*pow2(
        s2t))*pow8(Mst2)))))/(540.*pow4(Mst1)*pow5(Mst2)*pow8(Msq)))/pow4(Mt)/
        pow2(Sbeta)*12.; 

   return result;
}

/**
 * 	@return returns the susy log^3 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double himalaya::H6bq22g::calc_coef_at_as2_no_sm_logs_log3(){

   const double result =
      ((-224*pow2(Sbeta)*pow4(Mt))/9.)/pow4(Mt)/pow2(Sbeta)*12.; 

   return result;
}
