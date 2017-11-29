#include "H6g2.hpp"
#include "HierarchyCalculator.hpp"
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
himalaya::H6g2::H6g2(std::map<unsigned int, unsigned int> flagMap, double Al4p, double beta, double Dmglst2,
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

   // mdr flags, indicates if one wants to shift the dr stop mass to the mdr stop mass
   shiftst1 = mdrFlag;
   shiftst2 = mdrFlag;
   shiftst3 = mdrFlag;
   xDR2DRMOD = mdrFlag;
   // expansion flags
   xDmglst2 = flagMap.at(HierarchyCalculator::xxDmglst2);
   xMsq = flagMap.at(HierarchyCalculator::xxMsq);
   xMst = flagMap.at(HierarchyCalculator::xxMst);
   
   s1 = 
   #include "../hierarchies/h6g2/sigS1Full.inc"
   ;
   s2 = 
   #include "../hierarchies/h6g2/sigS2Full.inc"
   ;
   s12 = 
   #include "../hierarchies/h6g2/sigS12Full.inc"
   ;
}

/**
 * 	@return The diagonal (1, 1) matrix element of the Higgs mass matrix as a double for the hierarchy 'H6g2'
 */
double himalaya::H6g2::getS1(){
   return s1;
}

/**
 * 	@return The diagonal (2, 2) matrix element of the Higgs mass matrix as a double for the hierarchy 'H6g2'
 */
double himalaya::H6g2::getS2(){
   return s2;
}

/**
 * 	@return The off-diagonal (1, 2) = (2, 1) matrix element of the Higgs mass matrix as a double for the hierarchy 'H6g2'
 */
double himalaya::H6g2::getS12(){
   return s12;
}

/**
 * 	@return returns the constant term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double himalaya::H6g2::calc_at_as2_no_logs(){

   const double result =
      ((pow2(Sbeta)*(-39200*Mst2*s2t*pow2(Mst1)*pow3(Mt)*(pow2(Dmglst2)*(18900*
        pow2(Mst1)*pow4(Mst2)*(528*pow2(Mst1)*pow2(Mst2) + 155*pow4(Mst1) + 36*
        (175 - 116*z2)*pow4(Mst2)) - 75600*pow2(Msq)*pow2(Mst1)*pow2(Mst2)*(-5*
        pow2(Mst1)*pow2(Mst2) + 270*pow4(Mst1) + 6*(-259 + 156*z2)*pow4(Mst2))
        + pow4(Msq)*(-3*pow2(Mst2)*(-90713501 + 4823840*OepS2 + 82155924*S2 -
        7235760*T1ep - 126497694*z2 + 28621320*z3 - 3617880*z4 - 5426820*pow2(
        z2))*pow4(Mst1) - 9*pow2(Mst1)*(-36729361 + 7840*OepS2 + 3032964*S2 -
        11760*T1ep + 4445706*z2 + 23499840*z3 - 731640*z4 - 8820*pow2(z2))*
        pow4(Mst2) + 28*(-24208447 + 2058400*OepS2 - 47545272*S2 - 3087600*T1ep
        - 55366602*z2 + 36882120*z3 - 1543800*z4 - 2315700*pow2(z2))*pow6(Mst1)
        - 17418240*pow6(Mst2))) + 7*Dmglst2*Mgl*(2700*pow2(Mst1)*pow4(Mst2)*(
        246*pow2(Mst1)*pow2(Mst2) + 155*pow4(Mst1) + 6*(377 - 264*z2)*pow4(
        Mst2)) - 10800*pow2(Msq)*pow2(Mst1)*pow2(Mst2)*(-131*pow2(Mst1)*pow2(
        Mst2) + 270*pow4(Mst1) + 6*(-125 + 84*z2)*pow4(Mst2)) + pow4(Msq)*(27*
        pow2(Mst2)*(-678821 + 73760*OepS2 - 564876*S2 - 110640*T1ep - 2628894*
        z2 + 1261320*z3 - 55320*z4 - 82980*pow2(z2))*pow4(Mst1) + 135*pow2(
        Mst1)*(-1983 + 1120*OepS2 + 43740*S2 - 1680*T1ep - 36522*z2 - 68904*z3
        + 6072*z4 - 1260*pow2(z2))*pow4(Mst2) + 4*(-24208447 + 2058400*OepS2 -
        47545272*S2 - 3087600*T1ep - 55366602*z2 + 36882120*z3 - 1543800*z4 -
        2315700*pow2(z2))*pow6(Mst1) - 622080*pow6(Mst2))) - 21*pow2(Mgl)*(-
        3600*pow2(Msq)*pow2(Mst1)*pow2(Mst2)*(95*pow2(Mst1)*pow2(Mst2) + 126*
        pow4(Mst1) + 2*(119 - 108*z2)*pow4(Mst2)) - 900*pow2(Mst1)*pow4(Mst2)*(
        10*pow2(Mst1)*pow2(Mst2) + 47*pow4(Mst1) + 6*(103 - 72*z2)*pow4(Mst2))
        + pow4(Msq)*(3*pow2(Mst2)*(-1058993 + 60320*OepS2 - 1414908*S2 - 90480*
        T1ep - 2823942*z2 + 2907240*z3 - 45240*z4 - 67860*pow2(z2))*pow4(Mst1)
        + 27*pow2(Mst1)*(-15759 + 1120*OepS2 + 25596*S2 - 1680*T1ep - 173386*z2
        + 115320*z3 - 12360*z4 - 1260*pow2(z2))*pow4(Mst2) + 4*(-1907557 +
        109600*OepS2 - 3520152*S2 - 164400*T1ep - 3448554*z2 + 3465480*z3 -
        82200*z4 - 123300*pow2(z2))*pow6(Mst1) - 207360*pow6(Mst2)))) + 68600*
        Mt*pow2(Mst1)*pow3(Mst2)*pow3(s2t)*(pow2(Dmglst2)*(-64800*pow2(Msq)*
        pow2(Mst1)*pow2(Mst2)*((-143 + 72*z2)*pow2(Mst1)*pow2(Mst2) + 3*pow4(
        Mst1) + 2*(397 - 234*z2)*pow4(Mst2)) + 1350*pow2(Mst1)*pow4(Mst2)*(-9*(
        -593 + 360*z2)*pow2(Mst1)*pow2(Mst2) + 69*pow4(Mst1) + (-21275 + 12888*
        z2)*pow4(Mst2)) + pow4(Msq)*(60*pow2(Mst2)*(855985 + 11016*B4 + 324*DN
        - 24544*OepS2 - 89856*S2 + 36816*T1ep + 459510*z2 - 420780*z3 + 29100*
        z4 + 27612*pow2(z2))*pow4(Mst1) + 9*pow2(Mst1)*(986407 + 1058400*B4 -
        23760*DN - 1120*OepS2 - 449388*S2 + 1680*T1ep - 8043222*z2 + 15388440*
        z3 + 149880*z4 - 1864980*pow2(z2))*pow4(Mst2) + (-35542243 + 2491360*
        OepS2 - 90290268*S2 - 3737040*T1ep - 49166898*z2 + 33912840*z3 -
        1868520*z4 - 2802780*pow2(z2))*pow6(Mst1) - 2488320*pow6(Mst2))) +
        Dmglst2*Mgl*(-194400*pow2(Msq)*pow2(Mst1)*pow2(Mst2)*(12*(-3 + z2)*
        pow2(Mst1)*pow2(Mst2) + pow4(Mst1) + (95 - 48*z2)*pow4(Mst2)) + 1350*
        pow2(Mst1)*pow4(Mst2)*(-3*(-823 + 360*z2)*pow2(Mst1)*pow2(Mst2) + 69*
        pow4(Mst1) + (-5591 + 2952*z2)*pow4(Mst2)) + pow4(Msq)*(18*pow2(Mst2)*(
        -231989 + 36720*B4 + 1080*DN + 64160*OepS2 - 1515564*S2 - 96240*T1ep -
        1402446*z2 + 375000*z3 - 12480*z4 - 72180*pow2(z2))*pow4(Mst1) + 27*
        pow2(Mst1)*(825997 + 188640*B4 - 3600*DN + 5600*OepS2 + 146772*S2 -
        8400*T1ep - 1670082*z2 + 2248440*z3 + 32520*z4 - 317340*pow2(z2))*pow4(
        Mst2) + (-35542243 + 2491360*OepS2 - 90290268*S2 - 3737040*T1ep -
        49166898*z2 + 33912840*z3 - 1868520*z4 - 2802780*pow2(z2))*pow6(Mst1) -
        622080*pow6(Mst2))) + 3*pow2(Mgl)*(-1350*(pow2(Mst1) + (-223 + 72*z2)*
        pow2(Mst2))*pow4(Mst1)*pow4(Mst2) + 21600*pow2(Msq)*pow2(Mst1)*pow2(
        Mst2)*(2*(25 - 6*z2)*pow2(Mst1)*pow2(Mst2) + 3*pow4(Mst1) + (-71 + 24*
        z2)*pow4(Mst2)) + pow4(Msq)*(30*pow2(Mst2)*(18641 + 7344*B4 + 216*DN -
        2336*OepS2 + 112428*S2 + 3504*T1ep + 47478*z2 - 104136*z3 + 8880*z4 +
        2628*pow2(z2))*pow4(Mst1) + 135*pow2(Mst1)*(51531 + 5280*B4 - 48*DN -
        224*OepS2 - 324*S2 + 336*T1ep - 20814*z2 + 41960*z3 + 2040*z4 - 6660*
        pow2(z2))*pow4(Mst2) + (2711335 - 70240*OepS2 + 3305340*S2 + 105360*
        T1ep + 629034*z2 - 2677800*z3 + 52680*z4 + 79020*pow2(z2))*pow6(Mst1) +
        207360*pow6(Mst2)))) + pow2(Mst2)*pow4(s2t)*(-4*Dmglst2*Mgl*pow2(Mst1)*
        pow2(Mst2)*(17364375*((-113 + 120*z2)*pow2(Mst1) + (793 - 504*z2)*pow2(
        Mst2))*pow4(Mst1)*pow4(Mst2) + 69457500*pow2(Msq)*pow2(Mst1)*pow2(Mst2)
        *(-4*(-77 + 60*z2)*pow2(Mst1)*pow2(Mst2) + (25 + 36*z2)*pow4(Mst1) + 4*
        (-115 + 51*z2)*pow4(Mst2)) + pow4(Msq)*(-12348*(1727459 + 81000*B4 -
        108000*D3 + 67500*DN - 432000*OepS2 - 7425000*S2 + 648000*T1ep +
        3502800*z2 + 5037750*z3 - 121500*z4)*pow2(Mst2)*pow4(Mst1) + 2315250*
        pow2(Mst1)*(-7337 + 432*B4 - 576*D3 + 360*DN + 224*OepS2 + 543348*S2 -
        336*T1ep - 18282*z2 - 97050*z3 - 2544*z4 - 2844*pow2(z2))*pow4(Mst2) +
        (-126628974763 + 11522056000*OepS2 - 478191735000*S2 - 17283084000*T1ep
        - 226389363900*z2 + 101790223500*z3 - 8641542000*z4 - 12962313000*pow2(
        z2))*pow6(Mst1) - 958513500*(-5 + 8*z2)*pow6(Mst2))) - 3*pow2(Mgl)*
        pow2(Mst1)*(1852200*pow2(Msq)*pow2(Mst1)*pow4(Mst2)*(-8*(-659 + 375*z2)
        *pow2(Mst1)*pow2(Mst2) + (-2953 + 900*z2)*pow4(Mst1) + 4*(-1493 + 525*
        z2)*pow4(Mst2)) + 1666980000*(-1 + 2*z2)*(pow2(Mst1) + pow2(Mst2))*
        pow2(pow2(Mst1) - pow2(Mst2))*pow6(Msq) + 90*(9*(-2024719 + 857500*z2)*
        pow2(Mst1) + (45671209 - 21609000*z2)*pow2(Mst2))*pow4(Mst1)*pow6(Mst2)
        + pow2(Mst2)*pow4(Msq)*(-20580*pow2(Mst2)*(-9250423 + 10800*B4 - 21600*
        D3 + 16200*DN + 274400*OepS2 - 25409700*S2 - 411600*T1ep - 7052970*z2 +
        6763200*z3 - 286800*z4 - 503100*pow2(z2))*pow4(Mst1) + 3087000*pow2(
        Mst1)*(-9403 + 720*B4 - 72*D3 + 36*DN - 384*OepS2 + 171072*S2 + 576*
        T1ep + 3174*z2 - 34434*z3 + 180*z4 - 864*pow2(z2))*pow4(Mst2) + (-
        71378675209 - 2000376000*B4 - 222264000*D3 + 222264000*DN - 5337080000*
        OepS2 + 351825390000*S2 + 8005620000*T1ep + 77399836500*z2 -
        35058030000*z3 + 2669226000*z4 + 6004215000*pow2(z2))*pow6(Mst1) -
        27783000*(405 + 196*z2)*pow6(Mst2))) + 2*pow2(Dmglst2)*pow2(Mst2)*(-
        17364375*pow4(Mst1)*pow4(Mst2)*(-9*(-551 + 392*z2)*pow2(Mst1)*pow2(
        Mst2) + (-226 + 240*z2)*pow4(Mst1) + (-3287 + 2160*z2)*pow4(Mst2)) -
        138915000*pow2(Msq)*pow2(Mst1)*pow2(Mst2)*(2*(359 - 300*z2)*pow2(Mst2)*
        pow4(Mst1) + 2*(-502 + 273*z2)*pow2(Mst1)*pow4(Mst2) + (25 + 36*z2)*
        pow6(Mst1) + 48*pow6(Mst2)) - 2*pow4(Msq)*(77175*(-1200149 + 54000*B4 -
        60480*D3 + 33480*DN + 5600*OepS2 + 37625796*S2 - 8400*T1ep - 384546*z2
        - 6793410*z3 - 266640*z4 - 122940*pow2(z2))*pow4(Mst1)*pow4(Mst2) -
        2058*pow2(Mst2)*(-65965567 + 729000*B4 - 972000*D3 + 607500*DN -
        320000*OepS2 - 120274200*S2 + 480000*T1ep - 10372200*z2 + 104666250*z3
        - 3769500*z4 - 4014000*pow2(z2))*pow6(Mst1) - 20837250*(-1405 + 576*z2)
        *pow2(Mst1)*pow6(Mst2) + (-126628974763 + 11522056000*OepS2 -
        478191735000*S2 - 17283084000*T1ep - 226389363900*z2 + 101790223500*z3
        - 8641542000*z4 - 12962313000*pow2(z2))*pow8(Mst1) - 2667168000*pow8(
        Mst2)))) + 16*pow4(Mt)*(4*Dmglst2*Mgl*pow2(Mst1)*(5788125*((-1267 +
        1152*z2)*pow2(Mst1)*pow2(Mst2) + 498*pow4(Mst1) - 351*pow4(Mst2))*pow6(
        Mst2) - 138915000*pow2(Msq)*pow2(Mst2)*(6*pow2(Mst2)*pow4(Mst1) + 3*(5
        - 12*z2)*pow2(Mst1)*pow4(Mst2) + 9*(29 + 4*z2)*pow6(Mst1) + 56*pow6(
        Mst2)) + pow4(Msq)*(882*pow2(Mst2)*(23048051 + 4508000*OepS2 -
        17853750*S2 - 6762000*T1ep - 398279550*z2 + 257523000*z3 - 3381000*z4 -
        5071500*pow2(z2))*pow4(Mst1) + 66150*pow2(Mst1)*(993716 + 7840*OepS2 +
        258066*S2 - 11760*T1ep - 689046*z2 - 1249185*z3 - 5880*z4 - 8820*pow2(
        z2))*pow4(Mst2) + 5*(-72755578289 + 3063401600*OepS2 - 52952863320*S2 -
        4595102400*T1ep - 200231487120*z2 + 231412221600*z3 - 2297551200*z4 -
        3446326800*pow2(z2))*pow6(Mst1) + 41674500*(-243 + 184*z2)*pow6(Mst2)))
        - 3*pow2(Mgl)*(1666980000*(-1 + 2*z2)*pow2(Mst1)*pow2(Mst2)*(pow2(Mst1)
        + pow2(Mst2))*pow6(Msq) - 45*pow2(Mst1)*pow4(Mst2)*(2572500*pow2(Mst2)*
        pow4(Mst1) + (-142320737 + 86436000*z2)*pow2(Mst1)*pow4(Mst2) +
        47186763*pow6(Mst1) - 17235750*pow6(Mst2)) - 246960*pow2(Msq)*pow2(
        Mst1)*pow2(Mst2)*(22483*pow2(Mst2)*pow4(Mst1) + 9*(-2613 + 2500*z2)*
        pow2(Mst1)*pow4(Mst2) + 3375*(13 + 4*z2)*pow6(Mst1) - 16125*pow6(Mst2))
        + pow4(Msq)*(-154350*(326197 + 3360*OepS2 - 105948*S2 - 5040*T1ep +
        326298*z2 - 702360*z3 - 2520*z4 - 3780*pow2(z2))*pow4(Mst1)*pow4(Mst2)
        - 10290*pow2(Mst2)*(6906919 + 167200*OepS2 - 6744060*S2 - 250800*T1ep +
        9301170*z2 - 20715000*z3 - 125400*z4 - 188100*pow2(z2))*pow6(Mst1) -
        27783000*(277 + 196*z2)*pow2(Mst1)*pow6(Mst2) + (-175715865103 -
        3204992000*OepS2 + 128216692800*S2 + 4807488000*T1ep - 125675680200*z2
        + 399103824000*z3 + 2403744000*z4 + 3605616000*pow2(z2))*pow8(Mst1) +
        889056000*pow8(Mst2))) + 2*pow2(Dmglst2)*(231525*pow2(Mst1)*((-90263 +
        122400*z2)*pow2(Mst1)*pow2(Mst2) + 39150*pow4(Mst1) - 23625*pow4(Mst2))
        *pow6(Mst2) - 11113200*pow2(Msq)*pow2(Mst1)*pow2(Mst2)*(675*pow2(Mst2)*
        pow4(Mst1) - 2*(-512 + 675*z2)*pow2(Mst1)*pow4(Mst2) + 225*(29 + 4*z2)*
        pow6(Mst1) + 2000*pow6(Mst2)) + 2*pow4(Msq)*(-15435*(-22336511 + 72800*
        OepS2 + 3133890*S2 - 109200*T1ep - 2881650*z2 + 16944225*z3 - 54600*z4
        - 81900*pow2(z2))*pow4(Mst1)*pow4(Mst2) - 147*pow2(Mst2)*(-3582403037 +
        131684000*OepS2 + 173875950*S2 - 197526000*T1ep - 4696439250*z2 +
        3326757000*z3 - 98763000*z4 - 148144500*pow2(z2))*pow6(Mst1) +
        20837250*(-1021 + 576*z2)*pow2(Mst1)*pow6(Mst2) + 5*(-72755578289 +
        3063401600*OepS2 - 52952863320*S2 - 4595102400*T1ep - 200231487120*z2 +
        231412221600*z3 - 2297551200*z4 - 3446326800*pow2(z2))*pow8(Mst1) +
        2667168000*pow8(Mst2)))) + 24*pow2(Mst2)*pow2(Mt)*pow2(s2t)*(12*
        Dmglst2*Mgl*pow2(Mst1)*(pow4(Msq)*(1372*pow2(Mst2)*(2416741 + 300000*
        OepS2 - 753300*S2 - 450000*T1ep - 7425000*z2 - 1667250*z3 - 225000*z4 -
        337500*pow2(z2))*pow4(Mst1) - 308700*(-28373 + 240*B4 - 240*D3 + 120*DN
        + 34182*S2 + 6998*z2 + 4390*z3 - 1080*z4)*pow2(Mst1)*pow4(Mst2) + (-
        3966427397 + 2110136000*OepS2 - 52385772600*S2 - 3165204000*T1ep -
        40405228500*z2 + 755286000*z3 - 1582602000*z4 - 2373903000*pow2(z2))*
        pow6(Mst1) - 4630500*(-179 + 184*z2)*pow6(Mst2)) + 15435000*pow2(Msq)*(
        45*pow4(Mst1)*pow4(Mst2) + 24*pow2(Mst2)*pow6(Mst1) - 107*pow2(Mst1)*
        pow6(Mst2) + 56*pow8(Mst2)) + 1929375*(9*pow4(Mst2)*pow6(Mst1) + 65*
        pow4(Mst1)*pow6(Mst2) - 203*pow2(Mst1)*pow8(Mst2))) - 6*pow2(Dmglst2)*(
        -1929375*pow2(Mst1)*pow4(Mst2)*(31*pow2(Mst2)*pow4(Mst1) - 481*pow2(
        Mst1)*pow4(Mst2) + 18*pow6(Mst1) + 81*pow6(Mst2)) + 2*pow4(Msq)*(11025*
        (-2697747 + 10080*B4 - 10080*D3 + 5040*DN - 7840*OepS2 + 1126548*S2 +
        11760*T1ep + 323562*z2 + 2068815*z3 - 39480*z4 + 8820*pow2(z2))*pow4(
        Mst1)*pow4(Mst2) - 49*pow2(Mst2)*(140712871 + 16716000*OepS2 -
        343302300*S2 - 25074000*T1ep - 568446750*z2 + 71859375*z3 - 12537000*z4
        - 18805500*pow2(z2))*pow6(Mst1) + 2315250*(-1085 + 576*z2)*pow2(Mst1)*
        pow6(Mst2) + (3966427397 - 2110136000*OepS2 + 52385772600*S2 +
        3165204000*T1ep + 40405228500*z2 - 755286000*z3 + 1582602000*z4 +
        2373903000*pow2(z2))*pow8(Mst1) + 296352000*pow8(Mst2)) - 15435000*
        pow2(Msq)*(147*pow4(Mst2)*pow6(Mst1) - 313*pow4(Mst1)*pow6(Mst2) + 48*
        pow2(Mst2)*pow8(Mst1) + 160*pow2(Mst1)*pow8(Mst2))) + pow2(Mgl)*(
        1666980000*(-1 + 2*z2)*pow2(-(Mst1*pow2(Mst2)) + pow3(Mst1))*pow6(Msq)
        + pow4(Msq)*(77175*(818213 + 2880*D3 - 2880*DN - 11040*OepS2 - 866052*
        S2 + 16560*T1ep - 107538*z2 + 291960*z3 + 25560*z4 + 12420*pow2(z2))*
        pow4(Mst1)*pow4(Mst2) - 82320*pow2(Mst2)*(-517702 + 2700*D3 - 2700*DN +
        64400*OepS2 - 1937115*S2 - 96600*T1ep - 1108965*z2 + 141000*z3 - 32100*
        z4 - 72450*pow2(z2))*pow6(Mst1) - 27783000*(309 + 196*z2)*pow2(Mst1)*
        pow6(Mst2) + 2*(24574447162 - 4745062000*OepS2 + 188288632350*S2 +
        7117593000*T1ep + 85441188525*z2 - 17021203500*z3 + 3558796500*z4 +
        5338194750*pow2(z2))*pow8(Mst1) + 889056000*pow8(Mst2)) + 246960*pow2(
        Msq)*(-3091*pow4(Mst2)*pow6(Mst1) - 22642*pow4(Mst1)*pow6(Mst2) + 5108*
        pow2(Mst2)*pow8(Mst1) + 16125*pow2(Mst1)*pow8(Mst2)) - 5400*(-130598*
        pow6(Mst1)*pow6(Mst2) + 109502*pow4(Mst2)*pow8(Mst1) + 242073*pow4(
        Mst1)*pow8(Mst2))))))/(6.001128e9*pow2(Mgl)*pow4(Msq)*pow4(Mst1)*pow4(
        Mst2)))/pow4(Mt)/pow2(Sbeta)*12.;
 
   return result;
}

/**
 * 	@return returns the susy log^0 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double himalaya::H6g2::calc_coef_at_as2_no_sm_logs_log0(){

   const double result =
      (pow2(Sbeta)*((pow2(Dmglst2)*(6*pow2(Mst1)*pow2(Mst2)*(-18*pow2(Mst2)*
        pow2(Mt)*pow2(s2t)*(3184*T1ep - 9125*z3 + 1592*z4 + 2388*pow2(z2)) - 8*
        Mst2*s2t*(17228*T1ep - 68146*z3 + 8614*z4 + 12921*pow2(z2))*pow3(Mt) +
        8*Mt*(3068*T1ep - 35065*z3 + 2425*z4 + 2301*pow2(z2))*pow3(Mst2)*pow3(
        s2t) + 16*(18812*T1ep - 316834*z3 + 9406*z4 + 14109*pow2(z2))*pow4(Mt)
        + (640*T1ep + 139555*z3 - 5026*z4 - 5352*pow2(z2))*pow4(Mst2)*pow4(s2t)
        ) + 9*pow4(Mst2)*(-36*pow2(Mst2)*pow2(Mt)*pow2(s2t)*(112*T1ep + 19703*
        z3 - 376*z4 + 84*pow2(z2)) - 16*Mst2*s2t*(28*T1ep - 55952*z3 + 1742*z4
        + 21*pow2(z2))*pow3(Mt) + 4*Mt*(28*T1ep + 256474*z3 + 2498*z4 - 31083*
        pow2(z2))*pow3(Mst2)*pow3(s2t) + 8*(1456*T1ep - 225923*z3 + 728*z4 +
        1092*pow2(z2))*pow4(Mt) + (280*T1ep + 226447*z3 + 8888*z4 + 4098*pow2(
        z2))*pow4(Mst2)*pow4(s2t)) - 2*pow4(Mst1)*(144*pow2(Mst2)*pow2(Mt)*
        pow2(s2t)*(3076*T1ep - 734*z3 + 1538*z4 + 2307*pow2(z2)) - 32*Mst2*s2t*
        (51460*T1ep - 614702*z3 + 25730*z4 + 38595*pow2(z2))*pow3(Mt) + 2*Mt*(
        62284*T1ep - 565214*z3 + 31142*z4 + 46713*pow2(z2))*pow3(Mst2)*pow3(
        s2t) + 64*(11164*T1ep - 562226*z3 + 5582*z4 + 8373*pow2(z2))*pow4(Mt) -
        (33592*T1ep - 197843*z3 + 16796*z4 + 25194*pow2(z2))*pow4(Mst2)*pow4(
        s2t))) - 2*Dmglst2*Mgl*(-18*pow2(Mst1)*pow2(Mst2)*(-24*pow2(Mst2)*pow2(
        Mt)*pow2(s2t)*(200*T1ep + 741*z3 + 100*z4 + 150*pow2(z2)) + 12*Mst2*
        s2t*(1844*T1ep - 21022*z3 + 922*z4 + 1383*pow2(z2))*pow3(Mt) - 2*Mt*(
        1604*T1ep - 6250*z3 + 208*z4 + 1203*pow2(z2))*pow3(Mst2)*pow3(s2t) -
        16*(644*T1ep - 24526*z3 + 322*z4 + 483*pow2(z2))*pow4(Mt) + 3*(288*T1ep
        + 2239*z3 - 54*z4)*pow4(Mst2)*pow4(s2t)) - 27*pow4(Mst2)*(16*(-439*z3 +
        108*z4)*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 8*Mst2*s2t*(140*T1ep + 5742*z3
        - 506*z4 + 105*pow2(z2))*pow3(Mt) - 2*Mt*(140*T1ep - 37474*z3 - 542*z4
        + 5289*pow2(z2))*pow3(Mst2)*pow3(s2t) - 8*(112*T1ep + 11897*z3 + 56*z4
        + 84*pow2(z2))*pow4(Mt) + (56*T1ep + 16175*z3 + 424*z4 + 474*pow2(z2))*
        pow4(Mst2)*pow4(s2t)) + pow4(Mst1)*(144*pow2(Mst2)*pow2(Mt)*pow2(s2t)*(
        3076*T1ep - 734*z3 + 1538*z4 + 2307*pow2(z2)) - 32*Mst2*s2t*(51460*T1ep
        - 614702*z3 + 25730*z4 + 38595*pow2(z2))*pow3(Mt) + 2*Mt*(62284*T1ep -
        565214*z3 + 31142*z4 + 46713*pow2(z2))*pow3(Mst2)*pow3(s2t) + 64*(
        11164*T1ep - 562226*z3 + 5582*z4 + 8373*pow2(z2))*pow4(Mt) - (33592*
        T1ep - 197843*z3 + 16796*z4 + 25194*pow2(z2))*pow4(Mst2)*pow4(s2t))) -
        3*pow2(Mgl)*(18*pow4(Mst2)*(-6*pow2(Mst2)*pow2(Mt)*pow2(s2t)*(92*T1ep +
        1622*z3 + 142*z4 + 69*pow2(z2)) + 24*Mst2*s2t*(28*T1ep - 1922*z3 + 206*
        z4 + 21*pow2(z2))*pow3(Mt) - 2*Mt*(84*T1ep + 10490*z3 + 510*z4 - 1665*
        pow2(z2))*pow3(Mst2)*pow3(s2t) + 24*(28*T1ep + 3902*z3 + 14*z4 + 21*
        pow2(z2))*pow4(Mt) + 3*(32*T1ep - 1913*z3 + 10*z4 - 48*pow2(z2))*pow4(
        Mst2)*pow4(s2t)) + 6*pow2(Mst1)*pow2(Mst2)*(-16*pow2(Mst2)*pow2(Mt)*
        pow2(s2t)*(644*T1ep - 940*z3 + 214*z4 + 483*pow2(z2)) + 8*Mst2*s2t*(
        1508*T1ep - 48454*z3 + 754*z4 + 1131*pow2(z2))*pow3(Mt) - 4*Mt*(292*
        T1ep - 8678*z3 + 740*z4 + 219*pow2(z2))*pow3(Mst2)*pow3(s2t) + 8*(836*
        T1ep + 69050*z3 + 418*z4 + 627*pow2(z2))*pow4(Mt) + (1372*T1ep - 22544*
        z3 + 956*z4 + 1677*pow2(z2))*pow4(Mst2)*pow4(s2t)) + pow4(Mst1)*(-4*
        pow2(Mst2)*pow2(Mt)*pow2(s2t)*(27668*T1ep - 66166*z3 + 13834*z4 +
        20751*pow2(z2)) + 64*Mst2*s2t*(2740*T1ep - 57758*z3 + 1370*z4 + 2055*
        pow2(z2))*pow3(Mt) - 4*Mt*(1756*T1ep - 44630*z3 + 878*z4 + 1317*pow2(
        z2))*pow3(Mst2)*pow3(s2t) + 256*(292*T1ep + 24241*z3 + 146*z4 + 219*
        pow2(z2))*pow4(Mt) + (7780*T1ep - 34070*z3 + 2594*z4 + 5835*pow2(z2))*
        pow4(Mst2)*pow4(s2t)) - 648*z3*log(pow2(Mst2)/pow2(Mst1))*pow4(Mst2)*(
        8*(pow2(Mst1) - pow2(Mst2))*pow2(Mt)*pow2(s2t) - 16*pow4(Mt) + (-pow4(
        Mst1) + pow4(Mst2))*pow4(s2t))))/(5832.*pow2(Mgl)*pow4(Mst2)) + Mt*
        pow3(s2t)*((Mst2*pow2(Mst1)*(-695967*Dmglst2*Mgl*pow2(Msq) + 110160*B4*
        Dmglst2*Mgl*pow2(Msq) + 3240*Dmglst2*DN*Mgl*pow2(Msq) + 8559850*pow2(
        Dmglst2)*pow2(Msq) + 110160*B4*pow2(Dmglst2)*pow2(Msq) + 3240*DN*pow2(
        Dmglst2)*pow2(Msq) - 720*log(pow2(Mst2)/pow2(Mst1))*(-1227*Dmglst2*Mgl
        + 4672*pow2(Dmglst2) - 1467*pow2(Mgl))*pow2(Msq) + 279615*pow2(Mgl)*
        pow2(Msq) + 110160*B4*pow2(Mgl)*pow2(Msq) + 3240*DN*pow2(Mgl)*pow2(Msq)
        - 32400*Dmglst2*Mgl*pow2(Mst1) - 32400*pow2(Dmglst2)*pow2(Mst1) +
        32400*pow2(Mgl)*pow2(Mst1) - 97200*(2 + log(pow2(Mst2)/pow2(Mst1)))*(
        Dmglst2*Mgl + pow2(Dmglst2) + pow2(Mgl))*pow2(Msq)*pow2(log(pow2(Msq)/
        pow2(Mst1))) - 6480*(105*Dmglst2*Mgl + 239*pow2(Dmglst2) - 107*pow2(
        Mgl))*pow2(Msq)*pow2(log(pow2(Mst2)/pow2(Mst1))) + 32400*log(pow2(Msq)/
        pow2(Mst1))*pow2(Msq)*(12*Dmglst2*Mgl + 19*pow2(Dmglst2) - 18*pow2(Mgl)
        + 3*log(pow2(Mst2)/pow2(Mst1))*(7*Dmglst2*Mgl + 9*pow2(Dmglst2) + 3*
        pow2(Mgl)) + 6*(Dmglst2*Mgl + pow2(Dmglst2) + pow2(Mgl))*pow2(log(pow2(
        Mst2)/pow2(Mst1)))) - 518400*(Dmglst2*Mgl + pow2(Dmglst2) + pow2(Mgl))*
        pow2(Msq)*pow3(log(pow2(Mst2)/pow2(Mst1)))))/(14580.*pow2(Mgl)*pow2(
        Msq)) + ((92.98131001371742 - (10*log(pow2(Msq)/pow2(Mst1)))/3. + (4*(
        34 + 45*log(pow2(Msq)/pow2(Mst1)))*log(pow2(Mst2)/pow2(Mst1)))/27. + (
        820*pow2(log(pow2(Mst2)/pow2(Mst1))))/27. - (800*pow3(log(pow2(Mst2)/
        pow2(Mst1))))/27. - (Dmglst2*(Dmglst2 + Mgl)*(35542243 + 8755920*log(
        pow2(Mst2)/pow2(Mst1)) + 291600*log(pow2(Msq)/pow2(Mst1))*(-7 + 6*log(
        pow2(Mst2)/pow2(Mst1))) + 1438560*pow2(log(pow2(Mst2)/pow2(Mst1))) +
        2177280*pow3(log(pow2(Mst2)/pow2(Mst1)))))/(87480.*pow2(Mgl)))*pow4(
        Mst1))/Mst2 + (pow3(Mst2)*(777600*Dmglst2*Mgl*pow2(Msq)*pow2(Mst1) +
        1029600*pow2(Dmglst2)*pow2(Msq)*pow2(Mst1) + 360000*pow2(Mgl)*pow2(Msq)
        *pow2(Mst1) - 4320*pow2(Msq)*(15*pow2(Dmglst2)*(pow2(Msq) + 8*pow2(
        Mst1)) + pow2(Mgl)*(-281*pow2(Msq) + 20*pow2(Mst1)) + Dmglst2*Mgl*(-
        169*pow2(Msq) + 60*pow2(Mst1)))*pow2(log(pow2(Mst2)/pow2(Mst1))) +
        2477991*Dmglst2*Mgl*pow4(Msq) + 565920*B4*Dmglst2*Mgl*pow4(Msq) -
        10800*Dmglst2*DN*Mgl*pow4(Msq) + 986407*pow2(Dmglst2)*pow4(Msq) +
        1058400*B4*pow2(Dmglst2)*pow4(Msq) - 23760*DN*pow2(Dmglst2)*pow4(Msq) +
        2318895*pow2(Mgl)*pow4(Msq) + 237600*B4*pow2(Mgl)*pow4(Msq) - 2160*DN*
        pow2(Mgl)*pow4(Msq) - 64800*(-2 + log(pow2(Mst2)/pow2(Mst1)))*(Dmglst2*
        Mgl + pow2(Dmglst2) + pow2(Mgl))*pow2(log(pow2(Msq)/pow2(Mst1)))*pow4(
        Msq) - 69120*(9*Dmglst2*Mgl + 15*pow2(Dmglst2) + 5*pow2(Mgl))*pow3(log(
        pow2(Mst2)/pow2(Mst1)))*pow4(Msq) - 360*log(pow2(Mst2)/pow2(Mst1))*(
        pow2(Mgl)*(160*pow2(Msq)*pow2(Mst1) + 6978*pow4(Msq) - 15*pow4(Mst1)) +
        pow2(Dmglst2)*(2640*pow2(Msq)*pow2(Mst1) + 17774*pow4(Msq) - 15*pow4(
        Mst1)) + 15*Dmglst2*Mgl*(64*pow2(Msq)*pow2(Mst1) + 598*pow4(Msq) -
        pow4(Mst1))) + 5400*log(pow2(Msq)/pow2(Mst1))*(4*log(pow2(Mst2)/pow2(
        Mst1))*pow2(Msq)*(-3*Dmglst2*Mgl*(pow2(Msq) - 4*pow2(Mst1)) + pow2(Mgl)
        *(-15*pow2(Msq) + 4*pow2(Mst1)) + 3*pow2(Dmglst2)*(pow2(Msq) + 8*pow2(
        Mst1))) + 24*(Dmglst2*Mgl + pow2(Dmglst2) + pow2(Mgl))*pow2(log(pow2(
        Mst2)/pow2(Mst1)))*pow4(Msq) + pow2(Mgl)*(32*pow2(Msq)*pow2(Mst1) +
        144*pow4(Msq) - pow4(Mst1)) + Dmglst2*Mgl*(96*pow2(Msq)*pow2(Mst1) +
        360*pow4(Msq) - pow4(Mst1)) + pow2(Dmglst2)*(192*pow2(Msq)*pow2(Mst1) +
        788*pow4(Msq) - pow4(Mst1))) + 10350*Dmglst2*Mgl*pow4(Mst1) + 10350*
        pow2(Dmglst2)*pow4(Mst1) - 450*pow2(Mgl)*pow4(Mst1)))/(9720.*pow2(Mgl)*
        pow4(Msq)) - (S2*(Dmglst2*Mgl*(252594*pow2(Mst1)*pow2(Mst2) + 836021*
        pow4(Mst1) - 36693*pow4(Mst2)) - 15*pow2(Mgl)*(6246*pow2(Mst1)*pow2(
        Mst2) + 6121*pow4(Mst1) - 81*pow4(Mst2)) + pow2(Dmglst2)*(49920*pow2(
        Mst1)*pow2(Mst2) + 836021*pow4(Mst1) + 37449*pow4(Mst2)) + 30*log(pow2(
        Mst2)/pow2(Mst1))*(pow2(Dmglst2)*(-9204*pow2(Mst1)*pow2(Mst2) + 15571*
        pow4(Mst1) - 63*pow4(Mst2)) - 3*pow2(Mgl)*(438*pow2(Mst1)*pow2(Mst2) +
        439*pow4(Mst1) + 189*pow4(Mst2)) + Dmglst2*Mgl*(7218*pow2(Mst1)*pow2(
        Mst2) + 15571*pow4(Mst1) + 945*pow4(Mst2)))))/(810.*Mst2*pow2(Mgl)) + (
        (4*OepS2*(pow2(Dmglst2)*(-9204*pow2(Mst1)*pow2(Mst2) + 15571*pow4(Mst1)
        - 63*pow4(Mst2)) - 3*pow2(Mgl)*(438*pow2(Mst1)*pow2(Mst2) + 439*pow4(
        Mst1) + 189*pow4(Mst2)) + Dmglst2*Mgl*(7218*pow2(Mst1)*pow2(Mst2) +
        15571*pow4(Mst1) + 945*pow4(Mst2))))/(2187.*Mst2) + ((-68400*Dmglst2*
        Mgl*pow2(Msq)*pow2(Mst1) - 190560*pow2(Dmglst2)*pow2(Msq)*pow2(Mst1) -
        17040*pow2(Mgl)*pow2(Msq)*pow2(Mst1) - 27955*Dmglst2*Mgl*pow2(Mst1)*
        pow2(Mst2) - 106375*pow2(Dmglst2)*pow2(Mst1)*pow2(Mst2) + 60*log(pow2(
        Msq)/pow2(Mst1))*pow2(Mst1)*(3*pow2(Mgl)*(-32*pow2(Msq) + 13*pow2(Mst1)
        ) + pow2(Dmglst2)*(-576*pow2(Msq) + 543*pow2(Mst1) - 541*pow2(Mst2)) +
        Dmglst2*Mgl*(-288*pow2(Msq) + 183*pow2(Mst1) - 181*pow2(Mst2)) + 6*log(
        pow2(Mst2)/pow2(Mst1))*(pow2(Mgl)*(8*pow2(Msq) + 3*pow2(Mst1)) + 3*
        Dmglst2*Mgl*(8*pow2(Msq) + 5*(pow2(Mst1) + pow2(Mst2))) + pow2(Dmglst2)
        *(48*pow2(Msq) + 45*(pow2(Mst1) + pow2(Mst2))))) - 2304*Dmglst2*Mgl*
        pow4(Msq) - 9216*pow2(Dmglst2)*pow4(Msq) + 2304*pow2(Mgl)*pow4(Msq) +
        6*log(pow2(Mst2)/pow2(Mst1))*(pow2(Dmglst2)*(6240*pow2(Msq)*pow2(Mst1)
        + 5455*pow2(Mst1)*pow2(Mst2) + 768*pow4(Msq) - 5385*pow4(Mst1)) +
        Dmglst2*Mgl*(3840*pow2(Msq)*pow2(Mst1) + 2125*pow2(Mst1)*pow2(Mst2) -
        768*pow4(Msq) - 1515*pow4(Mst1)) + pow2(Mgl)*(1600*pow2(Msq)*pow2(Mst1)
        - 768*pow4(Msq) - 255*pow4(Mst1))) + 72*pow2(log(pow2(Mst2)/pow2(Mst1))
        )*(3*Dmglst2*Mgl*(-40*pow2(Msq)*pow2(Mst1) - 25*pow2(Mst1)*(pow2(Mst1)
        + pow2(Mst2)) + 32*pow4(Msq)) + 3*pow2(Dmglst2)*(-80*pow2(Msq)*pow2(
        Mst1) - 75*pow2(Mst1)*(pow2(Mst1) + pow2(Mst2)) + 64*pow4(Msq)) + pow2(
        Mgl)*(-40*pow2(Msq)*pow2(Mst1) + 32*pow4(Msq) - 15*pow4(Mst1))) +
        12345*Dmglst2*Mgl*pow4(Mst1) + 26685*pow2(Dmglst2)*pow4(Mst1) + 3345*
        pow2(Mgl)*pow4(Mst1))*pow5(Mst2))/(324.*pow2(Mst1)*pow4(Msq)))/pow2(
        Mgl)) + pow2(Mt)*pow2(s2t)*((40*(1 + 2*log(pow2(Msq)/pow2(Mst1)))*pow2(
        Msq))/3. + S2*((-2*(2604*Dmglst2*Mgl + 42383*pow2(Dmglst2) + 70*log(
        pow2(Mst2)/pow2(Mst1))*(300*Dmglst2*Mgl + 597*pow2(Dmglst2) - 322*pow2(
        Mgl)) - 33481*pow2(Mgl))*pow2(Mst1))/(105.*pow2(Mgl)) - (3*(11816*
        Dmglst2*Mgl + 13908*pow2(Dmglst2) + 70*log(pow2(Mst2)/pow2(Mst1))*(28*
        pow2(Dmglst2) - 23*pow2(Mgl)) + 6237*pow2(Mgl))*pow2(Mst2))/(70.*pow2(
        Mgl)) + ((1506.025925925926 + (6917*log(pow2(Mst2)/pow2(Mst1)))/9. - (
        4*Dmglst2*(Dmglst2 + Mgl)*(28283 + 23070*log(pow2(Mst2)/pow2(Mst1))))/(
        45.*pow2(Mgl)))*pow4(Mst1))/pow2(Mst2)) - ((20*(1 + 2*log(pow2(Msq)/
        pow2(Mst1)))*pow2(Msq)*pow2(Mst1))/3. - (196.55862427463637 - (
        305.24826908541195 + (40*log(pow2(Msq)/pow2(Mst1)))/3.)*log(pow2(Mst2)/
        pow2(Mst1)) + (83647*pow2(log(pow2(Mst2)/pow2(Mst1))))/945. + (Dmglst2*
        (Dmglst2 + Mgl)*(-3966427397 + 6573806820*log(pow2(Mst2)/pow2(Mst1)) +
        277830000*log(pow2(Msq)/pow2(Mst1))*(-4 + 5*log(pow2(Mst2)/pow2(Mst1)))
        + 3242408400*pow2(log(pow2(Mst2)/pow2(Mst1))) - 2234988000*pow3(log(
        pow2(Mst2)/pow2(Mst1)))))/(2.083725e7*pow2(Mgl)) - (674*pow3(log(pow2(
        Mst2)/pow2(Mst1))))/27.)*pow4(Mst1))/pow2(Mst2) + pow2(Mst2)*(
        252.5348765432099 + (8*D3)/9. - (8*DN)/9. + (100*log(pow2(Msq)/pow2(
        Mst1)))/3. - (20*(1 + 2*log(pow2(Msq)/pow2(Mst1)))*pow2(Msq))/(3.*pow2(
        Mst1)) + 20*pow2(log(pow2(Msq)/pow2(Mst1))) + (2*log(pow2(Mst2)/pow2(
        Mst1))*(-2717 - 240*log(pow2(Msq)/pow2(Mst1)) + 45*pow2(log(pow2(Msq)/
        pow2(Mst1)))))/27. + (2*(256 - 30*log(pow2(Msq)/pow2(Mst1)))*pow2(log(
        pow2(Mst2)/pow2(Mst1))))/9. + ((2*pow2(Mst1)*(33750*Dmglst2*Mgl +
        55125*pow2(Dmglst2) - 3091*pow2(Mgl) - 30*log(pow2(Mst2)/pow2(Mst1))*(
        1500*Dmglst2*Mgl + 3750*pow2(Dmglst2) + 341*pow2(Mgl)) + 30*log(pow2(
        Msq)/pow2(Mst1))*(1500*Dmglst2*Mgl + 3750*pow2(Dmglst2) + 157*pow2(Mgl)
        + 60*log(pow2(Mst2)/pow2(Mst1))*pow2(Mgl)) - 1800*pow2(Mgl)*pow2(log(
        pow2(Msq)/pow2(Mst1))) + 1350*pow2(Mgl)*pow2(log(pow2(Mst2)/pow2(Mst1))
        )))/(2025.*pow2(Msq)) + pow2(Dmglst2)*(1427.379365079365 - (16*B4)/3. +
        (16*D3)/3. - (8*DN)/3. - (910*log(pow2(Msq)/pow2(Mst1)))/9. - (4*(5987
        + 810*log(pow2(Msq)/pow2(Mst1)))*log(pow2(Mst2)/pow2(Mst1)))/27. + 60*
        pow2(log(pow2(Msq)/pow2(Mst1))) + (1486*pow2(log(pow2(Mst2)/pow2(Mst1))
        ))/9. - (128*pow3(log(pow2(Mst2)/pow2(Mst1))))/3.))/pow2(Mgl) + (16*
        pow3(log(pow2(Mst2)/pow2(Mst1))))/9. - (2*Dmglst2*(-28373 + 240*B4 -
        240*D3 + 120*DN + 33520*log(pow2(Mst2)/pow2(Mst1)) + 150*log(pow2(Msq)/
        pow2(Mst1))*(5 + 36*log(pow2(Mst2)/pow2(Mst1))) - 2700*pow2(log(pow2(
        Msq)/pow2(Mst1))) - 11730*pow2(log(pow2(Mst2)/pow2(Mst1))) + 1920*pow3(
        log(pow2(Mst2)/pow2(Mst1)))))/(135.*Mgl) - ((2.3647986178598424 + (991*
        log(pow2(Msq)/pow2(Mst1)))/441. - (2.511111111111111 + log(pow2(Msq)/
        pow2(Mst1)))*log(pow2(Mst2)/pow2(Mst1)) + (5*Dmglst2*(Dmglst2 + Mgl)*(-
        1 + 4*log(pow2(Msq)/pow2(Mst1)) - 4*log(pow2(Mst2)/pow2(Mst1))))/(6.*
        pow2(Mgl)) + (5*pow2(log(pow2(Msq)/pow2(Mst1))))/21.)*pow4(Mst1))/pow4(
        Msq)) + (-(pow2(Mst1)*(-67668748*Dmglst2*Mgl - 140712871*pow2(Dmglst2)
        - 472500*log(pow2(Msq)/pow2(Mst1))*(-32*Dmglst2*Mgl + 76*pow2(Dmglst2)
        + 8*log(pow2(Mst2)/pow2(Mst1))*(6*Dmglst2*Mgl + 6*pow2(Dmglst2) - pow2(
        Mgl)) - 3*pow2(Mgl)) - 72478280*pow2(Mgl) + 378000*D3*pow2(Mgl) -
        378000*DN*pow2(Mgl) + 420*log(pow2(Mst2)/pow2(Mst1))*(119314*Dmglst2*
        Mgl + 548953*pow2(Dmglst2) + 239965*pow2(Mgl)) + 1417500*(6*Dmglst2*Mgl
        + 9*pow2(Dmglst2) + 3*pow2(Mgl) + log(pow2(Mst2)/pow2(Mst1))*pow2(Mgl))
        *pow2(log(pow2(Msq)/pow2(Mst1))) + 6300*(-998*Dmglst2*Mgl + 15049*pow2(
        Dmglst2) - 5975*pow2(Mgl))*pow2(log(pow2(Mst2)/pow2(Mst1))) + 63000*(
        578*Dmglst2*Mgl + 863*pow2(Dmglst2) + 269*pow2(Mgl))*pow3(log(pow2(
        Mst2)/pow2(Mst1)))))/425250. - ((-1114750*Dmglst2*Mgl*pow2(Mst1) -
        265825*pow2(Dmglst2)*pow2(Mst1) - 522392*pow2(Mgl)*pow2(Mst1) + 28*log(
        pow2(Mst2)/pow2(Mst1))*(4*pow2(Mgl)*(13299*pow2(Mst1) - 14279*pow2(
        Mst2)) + 11025*pow2(Dmglst2)*(55*pow2(Mst1) - 101*pow2(Mst2)) + 66150*
        Dmglst2*Mgl*(3*pow2(Mst1) - 5*pow2(Mst2))) + 3481450*Dmglst2*Mgl*pow2(
        Mst2) + 4124575*pow2(Dmglst2)*pow2(Mst2) + 968292*pow2(Mgl)*pow2(Mst2)
        - 420*log(pow2(Msq)/pow2(Mst1))*(pow2(Mgl)*(2929*pow2(Mst1) - 3909*
        pow2(Mst2)) + 735*pow2(Dmglst2)*(55*pow2(Mst1) - 101*pow2(Mst2)) +
        4410*Dmglst2*Mgl*(3*pow2(Mst1) - 5*pow2(Mst2)) - 7*log(pow2(Mst2)/pow2(
        Mst1))*pow2(Mgl)*(33*pow2(Mst1) + 37*pow2(Mst2))) + 44100*pow2(Mgl)*(
        pow2(Mst1) - pow2(Mst2))*pow2(log(pow2(Msq)/pow2(Mst1))) - 5880*pow2(
        Mgl)*(24*pow2(Mst1) + 11*pow2(Mst2))*pow2(log(pow2(Mst2)/pow2(Mst1))))*
        pow4(Mst2))/(185220.*pow4(Msq)) + (4*OepS2*(24*Dmglst2*Mgl*pow2(Mst1)*(
        769*pow2(Mst1) + 150*pow2(Mst2)) + 12*pow2(Dmglst2)*(597*pow2(Mst1)*
        pow2(Mst2) + 1538*pow4(Mst1) + 63*pow4(Mst2)) - pow2(Mgl)*(3864*pow2(
        Mst1)*pow2(Mst2) + 6917*pow4(Mst1) + 621*pow4(Mst2))))/(729.*pow2(Mst2)
        ) + ((80550*Dmglst2*Mgl*pow2(Msq) + 244125*pow2(Dmglst2)*pow2(Msq) -
        69525*pow2(Mgl)*pow2(Msq) - 160500*Dmglst2*Mgl*pow2(Mst1) - 234750*
        pow2(Dmglst2)*pow2(Mst1) - 45284*pow2(Mgl)*pow2(Mst1) - 30*log(pow2(
        Mst2)/pow2(Mst1))*(15*pow2(Dmglst2)*(59*pow2(Msq) - 1000*pow2(Mst1)) -
        30*Dmglst2*Mgl*(119*pow2(Msq) + 200*pow2(Mst1)) - pow2(Mgl)*(4335*pow2(
        Msq) + 1082*pow2(Mst1))) + 30*log(pow2(Msq)/pow2(Mst1))*(15*log(pow2(
        Mst2)/pow2(Mst1))*(180*Dmglst2*Mgl*pow2(Msq) + 270*pow2(Dmglst2)*pow2(
        Msq) + pow2(Mgl)*(90*pow2(Msq) - 13*pow2(Mst1))) + 150*Dmglst2*Mgl*(3*
        pow2(Msq) - 40*pow2(Mst1)) + 375*pow2(Dmglst2)*(9*pow2(Msq) - 40*pow2(
        Mst1)) - pow2(Mgl)*(1125*pow2(Msq) + 1282*pow2(Mst1))) - 450*(90*
        Dmglst2*Mgl*pow2(Msq) + 135*pow2(Dmglst2)*pow2(Msq) + pow2(Mgl)*(45*
        pow2(Msq) - 4*pow2(Mst1)))*pow2(log(pow2(Msq)/pow2(Mst1))) - 450*(272*
        Dmglst2*Mgl*pow2(Msq) + 472*pow2(Dmglst2)*pow2(Msq) + 3*pow2(Mgl)*(40*
        pow2(Msq) - 3*pow2(Mst1)))*pow2(log(pow2(Mst2)/pow2(Mst1))))*pow4(Mst2)
        + 2*(-4500*Dmglst2*(Dmglst2 + Mgl)*(-4 + log(pow2(Mst2)/pow2(Mst1))) +
        pow2(Mgl)*(5108 + 7500*log(pow2(Mst2)/pow2(Mst1)) + 15*log(pow2(Msq)/
        pow2(Mst1))*(218 + 75*log(pow2(Mst2)/pow2(Mst1))) + 900*pow2(log(pow2(
        Msq)/pow2(Mst1)))))*pow6(Mst1))/(2025.*pow2(Msq)*pow2(Mst1)) + ((5*(
        896*Dmglst2*Mgl*pow2(Msq) + 1280*pow2(Dmglst2)*pow2(Msq) + 344*pow2(
        Mgl)*pow2(Msq) + 81*pow2(Dmglst2)*pow2(Mst2) + 60*log(pow2(Msq)/pow2(
        Mst1))*(16*Dmglst2*Mgl*pow2(Msq) + 4*pow2(Mgl)*pow2(Msq) + pow2(
        Dmglst2)*(40*pow2(Msq) + 21*pow2(Mst2))) - 60*log(pow2(Mst2)/pow2(Mst1)
        )*(16*Dmglst2*Mgl*pow2(Msq) + 4*pow2(Mgl)*pow2(Msq) + pow2(Dmglst2)*(
        40*pow2(Msq) + 21*pow2(Mst2)))))/(108.*pow2(Mst1)*pow4(Msq)) + (32*(-4*
        pow2(Dmglst2) + pow2(Mgl) - 2*log(pow2(Mst2)/pow2(Mst1))*pow2(Dmglst2 +
        Mgl) + (4*Dmglst2*Mgl + 10*pow2(Dmglst2) + pow2(Mgl))*pow2(log(pow2(
        Mst2)/pow2(Mst1)))))/(9.*pow4(Mst1)))*pow6(Mst2))/pow2(Mgl)) - pow4(
        s2t)*((5*(1 + 2*log(pow2(Msq)/pow2(Mst1)))*(1 - 2*log(pow2(Mst2)/pow2(
        Mst1)))*pow2(Msq)*pow2(Mst1))/6. - (35.682629270197204 + B4 + D3/9. -
        DN/9. + (655*log(pow2(Msq)/pow2(Mst1)))/72. + (5*(1 + 2*log(pow2(Msq)/
        pow2(Mst1)))*pow2(Msq))/(6.*pow2(Mst2)) + log(pow2(Mst2)/pow2(Mst1))*(
        33.61121882086168 + (65*log(pow2(Msq)/pow2(Mst1)))/18. + (5*pow2(log(
        pow2(Msq)/pow2(Mst1))))/12.) + (25*pow2(log(pow2(Msq)/pow2(Mst1))))/12.
         + (7.516137566137566 + (5*log(pow2(Msq)/pow2(Mst1)))/6.)*pow2(log(
        pow2(Mst2)/pow2(Mst1))) + (Dmglst2*(Dmglst2 + Mgl)*(84.40344866031853 +
        (10*log(pow2(Msq)/pow2(Mst1)))/3. - (1011403*log(pow2(Mst2)/pow2(Mst1))
        )/1.1907e6 + (5*pow2(log(pow2(Msq)/pow2(Mst1))))/2. - (113*pow2(log(
        pow2(Mst2)/pow2(Mst1))))/35. - (11*pow3(log(pow2(Mst2)/pow2(Mst1))))/9.
        ))/pow2(Mgl) - (136*pow3(log(pow2(Mst2)/pow2(Mst1))))/27.)*pow4(Mst1) +
        pow2(Mst2)*((5*(1 + 2*log(pow2(Msq)/pow2(Mst1)))*(1 + 2*log(pow2(Mst2)/
        pow2(Mst1)))*pow2(Msq))/6. - (pow2(Mst1)*(20729508*Dmglst2*Mgl +
        972000*B4*Dmglst2*Mgl - 1296000*D3*Dmglst2*Mgl + 810000*Dmglst2*DN*Mgl
        - 131931134*pow2(Dmglst2) + 1458000*B4*pow2(Dmglst2) - 1944000*D3*pow2(
        Dmglst2) + 1215000*DN*pow2(Dmglst2) - 180*log(pow2(Mst2)/pow2(Mst1))*(
        257786*Dmglst2*Mgl + 619347*pow2(Dmglst2) - 184990*pow2(Mgl)) -
        138756345*pow2(Mgl) + 162000*B4*pow2(Mgl) - 324000*D3*pow2(Mgl) +
        243000*DN*pow2(Mgl) + 607500*(-6*Dmglst2*Mgl - 9*pow2(Dmglst2) - 7*
        pow2(Mgl) + 6*log(pow2(Mst2)/pow2(Mst1))*(2*Dmglst2*Mgl + 3*pow2(
        Dmglst2) + pow2(Mgl)))*pow2(log(pow2(Msq)/pow2(Mst1))) - 2700*(11288*
        Dmglst2*Mgl + 9056*pow2(Dmglst2) + 1445*pow2(Mgl))*pow2(log(pow2(Mst2)/
        pow2(Mst1))) - 202500*log(pow2(Msq)/pow2(Mst1))*(-54*Dmglst2*Mgl - 221*
        pow2(Dmglst2) + 6*log(pow2(Mst2)/pow2(Mst1))*(-6*Dmglst2*Mgl + 7*pow2(
        Dmglst2) - 7*pow2(Mgl)) + 89*pow2(Mgl) + 12*(6*Dmglst2*Mgl + 9*pow2(
        Dmglst2) + 2*pow2(Mgl))*pow2(log(pow2(Mst2)/pow2(Mst1)))) + 27000*(154*
        Dmglst2*Mgl + 211*pow2(Dmglst2) - 185*pow2(Mgl))*pow3(log(pow2(Mst2)/
        pow2(Mst1)))))/(1.458e6*pow2(Mgl)) - ((2953 + 1020*log(pow2(Msq)/pow2(
        Mst1)) + 12*(-53 + 30*log(pow2(Msq)/pow2(Mst1)))*log(pow2(Mst2)/pow2(
        Mst1)) + (50*Dmglst2*(Dmglst2 + Mgl)*(-25 + 48*log(pow2(Mst2)/pow2(
        Mst1)) + 12*log(pow2(Msq)/pow2(Mst1))*(-4 + 3*log(pow2(Mst2)/pow2(Mst1)
        )) - 36*pow2(log(pow2(Mst2)/pow2(Mst1)))))/pow2(Mgl) - 1080*pow2(log(
        pow2(Mst2)/pow2(Mst1))))*pow4(Mst1))/(1080.*pow2(Msq))) - (
        14.510802469135802 - (10*B4)/9. + D3/9. - DN/18. + (35*log(pow2(Msq)/
        pow2(Mst1)))/36. + (5*(1 + 2*log(pow2(Msq)/pow2(Mst1)))*pow2(Msq))/(6.*
        pow2(Mst1)) - (5*pow2(log(pow2(Msq)/pow2(Mst1))))/12. - (log(pow2(Mst2)
        /pow2(Mst1))*(6377 + 1770*log(pow2(Msq)/pow2(Mst1)) + 630*pow2(log(
        pow2(Msq)/pow2(Mst1)))))/216. + ((586 + 105*log(pow2(Msq)/pow2(Mst1)))*
        pow2(log(pow2(Mst2)/pow2(Mst1))))/18. - (pow2(Mst1)*(3850*Dmglst2*Mgl +
        8975*pow2(Dmglst2) + 1318*pow2(Mgl) + log(pow2(Mst2)/pow2(Mst1))*(2500*
        Dmglst2*Mgl + 5050*pow2(Dmglst2) + 346*pow2(Mgl)) + 15*log(pow2(Msq)/
        pow2(Mst1))*(20*Dmglst2*Mgl - 70*pow2(Dmglst2) + 23*pow2(Mgl) + log(
        pow2(Mst2)/pow2(Mst1))*(200*Dmglst2*Mgl + 500*pow2(Dmglst2) + 41*pow2(
        Mgl))) - 30*(100*Dmglst2*Mgl + 250*pow2(Dmglst2) + 19*pow2(Mgl))*pow2(
        log(pow2(Mst2)/pow2(Mst1)))))/(270.*pow2(Mgl)*pow2(Msq)) - (94*pow3(
        log(pow2(Mst2)/pow2(Mst1))))/9. - (Dmglst2*(-7337 + 432*B4 - 576*D3 +
        360*DN - 24804*log(pow2(Mst2)/pow2(Mst1)) + 1620*(1 + 2*log(pow2(Mst2)/
        pow2(Mst1)))*pow2(log(pow2(Msq)/pow2(Mst1))) - 8604*pow2(log(pow2(Mst2)
        /pow2(Mst1))) - 1620*log(pow2(Msq)/pow2(Mst1))*(-5 + 2*log(pow2(Mst2)/
        pow2(Mst1)) + 4*pow2(log(pow2(Mst2)/pow2(Mst1)))) + 10944*pow3(log(
        pow2(Mst2)/pow2(Mst1)))))/(648.*Mgl) - (pow2(Dmglst2)*(-1200149 +
        54000*B4 - 60480*D3 + 33480*DN - 1586700*log(pow2(Mst2)/pow2(Mst1)) +
        72900*(1 + 2*log(pow2(Mst2)/pow2(Mst1)))*pow2(log(pow2(Msq)/pow2(Mst1))
        ) + 313740*pow2(log(pow2(Mst2)/pow2(Mst1))) - 2700*log(pow2(Msq)/pow2(
        Mst1))*(-121 + 138*log(pow2(Mst2)/pow2(Mst1)) + 108*pow2(log(pow2(Mst2)
        /pow2(Mst1)))) + 492480*pow3(log(pow2(Mst2)/pow2(Mst1)))))/(19440.*
        pow2(Mgl)) + ((6074157 + 3874080*log(pow2(Msq)/pow2(Mst1)) + 35*(-34519
        + 78540*log(pow2(Msq)/pow2(Mst1)))*log(pow2(Mst2)/pow2(Mst1)) - 176400*
        pow2(log(pow2(Msq)/pow2(Mst1))) + (85750*Dmglst2*(Dmglst2 + Mgl)*(113 +
        34*log(pow2(Mst2)/pow2(Mst1)) + 12*log(pow2(Msq)/pow2(Mst1))*(-1 + 10*
        log(pow2(Mst2)/pow2(Mst1))) - 120*pow2(log(pow2(Mst2)/pow2(Mst1)))))/
        pow2(Mgl) - 3278100*pow2(log(pow2(Mst2)/pow2(Mst1))))*pow4(Mst1))/(
        7.4088e6*pow4(Msq)))*pow4(Mst2) + (S2*(((33000*Dmglst2*Mgl + 89092*
        pow2(Dmglst2) - 30*log(pow2(Mst2)/pow2(Mst1))*(1296*Dmglst2*Mgl + 160*
        pow2(Dmglst2) - 1029*pow2(Mgl)) + 141165*pow2(Mgl))*pow2(Mst1)*pow2(
        Mst2))/pow2(Mgl) + 5*(18995 + 5835*log(pow2(Mst2)/pow2(Mst1)) - (2*
        Dmglst2*(Dmglst2 + Mgl)*(51635 + 25194*log(pow2(Mst2)/pow2(Mst1))))/(3.
        *pow2(Mgl)))*pow4(Mst1) + (9*(50310*Dmglst2*Mgl + 116129*pow2(Dmglst2)
        - 10*log(pow2(Mst2)/pow2(Mst1))*(42*Dmglst2*Mgl + 35*pow2(Dmglst2) -
        72*pow2(Mgl)) + 15840*pow2(Mgl))*pow4(Mst2))/pow2(Mgl)))/540. + ((
        OepS2*(4*Dmglst2*Mgl*(1944*pow2(Mst1)*pow2(Mst2) + 4199*pow4(Mst1) +
        189*pow4(Mst2)) + 2*pow2(Dmglst2)*(480*pow2(Mst1)*pow2(Mst2) + 8398*
        pow4(Mst1) + 315*pow4(Mst2)) - 3*pow2(Mgl)*(2058*pow2(Mst1)*pow2(Mst2)
        + 1945*pow4(Mst1) + 432*pow4(Mst2))))/2187. + ((203999250*Dmglst2*Mgl*
        pow2(Mst1) + 637851375*pow2(Dmglst2)*pow2(Mst1) + 45671209*pow2(Mgl)*
        pow2(Mst1) - 210*log(pow2(Mst2)/pow2(Mst1))*(-7350*Dmglst2*Mgl*pow2(
        Mst1) + 45247*pow2(Mgl)*pow2(Mst1) + 3675*pow2(Dmglst2)*(pow2(Mst1) -
        274*pow2(Mst2))) + 420*log(pow2(Msq)/pow2(Mst1))*(139650*Dmglst2*Mgl*
        pow2(Mst1) + 51313*pow2(Mgl)*pow2(Mst1) + 210*log(pow2(Mst2)/pow2(Mst1)
        )*(1470*Dmglst2*Mgl*pow2(Mst1) + 211*pow2(Mgl)*pow2(Mst1) + 105*pow2(
        Dmglst2)*(49*pow2(Mst1) - 30*pow2(Mst2))) + 3675*pow2(Dmglst2)*(53*
        pow2(Mst1) - 165*pow2(Mst2))) - 422790375*pow2(Dmglst2)*pow2(Mst2) -
        352800*pow2(Mgl)*pow2(Mst1)*pow2(log(pow2(Msq)/pow2(Mst1))) - 264600*(
        490*Dmglst2*Mgl*pow2(Mst1) + 69*pow2(Mgl)*pow2(Mst1) + 35*pow2(Dmglst2)
        *(49*pow2(Mst1) - 30*pow2(Mst2)))*pow2(log(pow2(Mst2)/pow2(Mst1))))*
        pow6(Mst2))/(2.22264e7*pow4(Msq)) - ((-3450*Dmglst2*Mgl*pow2(Msq)*pow2(
        Mst1) - 21075*pow2(Dmglst2)*pow2(Msq)*pow2(Mst1) + 6075*pow2(Mgl)*pow2(
        Msq)*pow2(Mst1) + 1920*pow2(Dmglst2)*pow2(Msq)*pow2(Mst2) - 2400*pow2(
        Dmglst2)*pow2(Mst1)*pow2(Mst2) - 30*log(pow2(Msq)/pow2(Mst1))*pow2(
        Mst1)*(10*log(pow2(Mst2)/pow2(Mst1))*(pow2(Dmglst2)*(27*pow2(Msq) - 91*
        pow2(Mst1)) + 2*Dmglst2*Mgl*(9*pow2(Msq) - 17*pow2(Mst1)) + pow2(Mgl)*(
        9*pow2(Msq) - 8*pow2(Mst1))) + 10*Dmglst2*Mgl*(3*pow2(Msq) - 32*pow2(
        Mst1)) - pow2(Mgl)*(75*pow2(Msq) + 86*pow2(Mst1)) + 5*pow2(Dmglst2)*(
        45*pow2(Msq) - 136*pow2(Mst1) + 60*pow2(Mst2))) + 1350*(2*Dmglst2*Mgl +
        3*pow2(Dmglst2) + pow2(Mgl))*pow2(Msq)*pow2(Mst1)*pow2(log(pow2(Msq)/
        pow2(Mst1))) + 60*pow2(log(pow2(Mst2)/pow2(Mst1)))*(2*Dmglst2*Mgl*(116*
        pow2(Msq) - 85*pow2(Mst1))*pow2(Mst1) + 4*pow2(Mgl)*(21*pow2(Msq) - 10*
        pow2(Mst1))*pow2(Mst1) + pow2(Dmglst2)*(pow2(Msq)*(476*pow2(Mst1) - 48*
        pow2(Mst2)) - 455*pow4(Mst1))) + 23000*Dmglst2*Mgl*pow4(Mst1) + 50200*
        pow2(Dmglst2)*pow4(Mst1) + 5972*pow2(Mgl)*pow4(Mst1) - 10*log(pow2(
        Mst2)/pow2(Mst1))*(21*pow2(Mgl)*(55*pow2(Msq) - 7*pow2(Mst1))*pow2(
        Mst1) + 26*Dmglst2*Mgl*(57*pow2(Msq) - 5*pow2(Mst1))*pow2(Mst1) + pow2(
        Dmglst2)*(-900*pow2(Mst1)*pow2(Mst2) + pow2(Msq)*(591*pow2(Mst1) + 96*
        pow2(Mst2)) + 305*pow4(Mst1))))*pow6(Mst2))/(1080.*pow2(Msq)*pow4(Mst1)
        ))/pow2(Mgl)) + z2*(-(s2t*pow3(Mt)*((-3*(27603387*Dmglst2*Mgl -
        21082949*pow2(Dmglst2) - 9883797*pow2(Mgl) + 7560*log(pow2(Mst2)/pow2(
        Mst1))*(-505*Dmglst2*Mgl + 91*pow2(Dmglst2) + 495*pow2(Mgl)))*pow2(
        Mst1)*pow2(Mst2))/pow2(Mgl) + 28*(1724277 - 97200*log(pow2(Msq)/pow2(
        Mst1)) - 395280*log(pow2(Mst2)/pow2(Mst1)) + (Dmglst2*(Dmglst2 + Mgl)*(
        -9227767 + 291600*log(pow2(Msq)/pow2(Mst1)) + 1088640*log(pow2(Mst2)/
        pow2(Mst1))))/pow2(Mgl))*pow4(Mst1) - (9*pow4(Mst2)*(639135*Dmglst2*Mgl
        + 740951*pow2(Dmglst2) - 1820553*pow2(Mgl) - 302400*log(pow2(Msq)/pow2(
        Mst1))*(Dmglst2*Mgl + pow2(Dmglst2) + pow2(Mgl)) + 2520*log(pow2(Mst2)/
        pow2(Mst1))*(359*Dmglst2*Mgl + 379*pow2(Dmglst2) + 439*pow2(Mgl)) + (
        100800*(7*Dmglst2*Mgl + 13*pow2(Dmglst2) + 3*pow2(Mgl))*pow2(Mst2))/
        pow2(Msq) + (50400*(11*Dmglst2*Mgl + 29*pow2(Dmglst2) + 3*pow2(Mgl))*
        pow4(Mst2))/pow4(Msq)))/pow2(Mgl)))/(25515.*pow3(Mst2)) + (
        402.837037037037 - (80*log(pow2(Msq)/pow2(Mst1)))/3. - (556*log(pow2(
        Mst2)/pow2(Mst1)))/3. + (4*Dmglst2*(-114841 + 31920*log(pow2(Mst2)/
        pow2(Mst1))))/(945.*Mgl) - (80*pow2(Msq))/(3.*pow2(Mst1)) - ((80*pow2(
        Msq))/3. + (pow2(Mst1)*(31862364*Dmglst2*Mgl*pow2(Msq) - 62619190*pow2(
        Dmglst2)*pow2(Msq) + 453600*log(pow2(Msq)/pow2(Mst1))*pow2(Dmglst2 -
        Mgl)*pow2(Msq) - 6510819*pow2(Mgl)*pow2(Msq) + 7560*log(pow2(Mst2)/
        pow2(Mst1))*(-1236*Dmglst2*Mgl + 538*pow2(Dmglst2) + 343*pow2(Mgl))*
        pow2(Msq) + 453600*Dmglst2*Mgl*pow2(Mst1) + 453600*pow2(Dmglst2)*pow2(
        Mst1) - 226800*pow2(Mgl)*pow2(Mst1)))/(8505.*pow2(Mgl)*pow2(Msq)))/
        pow2(Mst2) + ((1005.2164609053498 - 240*log(pow2(Msq)/pow2(Mst1)) - (
        1232*log(pow2(Mst2)/pow2(Mst1)))/3. - (Dmglst2*(Dmglst2 + Mgl)*(
        10677.005369390554 - 960*log(pow2(Msq)/pow2(Mst1)) - (8128*log(pow2(
        Mst2)/pow2(Mst1)))/3.))/pow2(Mgl))*pow4(Mst1))/pow4(Mst2) + (2*((19211
        + 5328*log(pow2(Mst2)/pow2(Mst1)))*pow2(Dmglst2) + 36*((10*(6*Dmglst2*
        Mgl + 9*pow2(Dmglst2) + 5*pow2(Mgl)))/pow2(Msq) + (92*Dmglst2*Mgl +
        144*pow2(Dmglst2) + 49*pow2(Mgl))/pow2(Mst1))*pow2(Mst2) + (180*(16*
        Dmglst2*Mgl + 34*pow2(Dmglst2) + 7*pow2(Mgl))*pow4(Mst2))/pow4(Msq)))/(
        81.*pow2(Mgl)))*pow4(Mt) + Mt*pow3(s2t)*((Mst2*(-701223*Dmglst2*Mgl +
        765850*pow2(Dmglst2) + 1080*log(pow2(Mst2)/pow2(Mst1))*(99*Dmglst2*Mgl
        + 274*pow2(Dmglst2) - 147*pow2(Mgl)) + 118695*pow2(Mgl) + 32400*log(
        pow2(Msq)/pow2(Mst1))*(Dmglst2*Mgl + pow2(Dmglst2) + pow2(Mgl)))*pow2(
        Mst1))/(2430.*pow2(Mgl)) - ((835041*Dmglst2*Mgl*pow2(Msq) + 1340537*
        pow2(Dmglst2)*pow2(Msq) + 156105*pow2(Mgl)*pow2(Msq) + 21600*log(pow2(
        Msq)/pow2(Mst1))*(11*Dmglst2*Mgl + 23*pow2(Dmglst2) + 3*pow2(Mgl))*
        pow2(Msq) - 360*log(pow2(Mst2)/pow2(Mst1))*(993*Dmglst2*Mgl + 2163*
        pow2(Dmglst2) + 173*pow2(Mgl))*pow2(Msq) + 43200*Dmglst2*Mgl*pow2(Mst1)
        + 86400*pow2(Dmglst2)*pow2(Mst1) + 14400*pow2(Mgl)*pow2(Mst1))*pow3(
        Mst2))/(1620.*pow2(Mgl)*pow2(Msq)) + ((21.57181069958848 - 34*log(pow2(
        Mst2)/pow2(Mst1)) + (Dmglst2*(Dmglst2 + Mgl)*(-8194483 + 2329560*log(
        pow2(Mst2)/pow2(Mst1))))/(14580.*pow2(Mgl)))*pow4(Mst1))/Mst2 + (10*(
        pow2(Mgl)*(16*pow2(Msq) - 3*pow2(Mst1)) + Dmglst2*Mgl*(96*pow2(Msq) -
        15*pow2(Mst1) + 41*pow2(Mst2)) + pow2(Dmglst2)*(312*pow2(Msq) - 45*
        pow2(Mst1) + 179*pow2(Mst2)))*pow5(Mst2))/(9.*pow2(Mgl)*pow4(Msq))) +
        pow4(s2t)*((5*(1 - 2*log(pow2(Mst2)/pow2(Mst1)))*pow2(Msq)*pow2(Mst1))/
        3. - (38.69264403292181 - (5*log(pow2(Msq)/pow2(Mst1)))/6. + (25*log(
        pow2(Mst2)/pow2(Mst1)))/12. - (Dmglst2*(Dmglst2 + Mgl)*(
        150.89787379972566 - (44*log(pow2(Mst2)/pow2(Mst1)))/3.))/pow2(Mgl) + (
        5*pow2(Msq))/(3.*pow2(Mst2)))*pow4(Mst1) + pow2(Mst2)*((5*(1 + 2*log(
        pow2(Mst2)/pow2(Mst1)))*pow2(Msq))/3. - ((-280224*Dmglst2*Mgl + 138296*
        pow2(Dmglst2) + 705297*pow2(Mgl) + 48600*log(pow2(Msq)/pow2(Mst1))*(2*
        Dmglst2*Mgl + 3*pow2(Dmglst2) + pow2(Mgl)) - 540*log(pow2(Mst2)/pow2(
        Mst1))*(446*Dmglst2*Mgl + 645*pow2(Dmglst2) + 335*pow2(Mgl)))*pow2(
        Mst1))/(9720.*pow2(Mgl)) - (5*(2*Dmglst2*Mgl + 2*pow2(Dmglst2) + pow2(
        Mgl))*pow4(Mst1))/(6.*pow2(Mgl)*pow2(Msq))) - ((1058 - 900*log(pow2(
        Msq)/pow2(Mst1)) - (2*Dmglst2*(3047 + 1080*log(pow2(Msq)/pow2(Mst1)) -
        3480*log(pow2(Mst2)/pow2(Mst1))))/Mgl + 2724*log(pow2(Mst2)/pow2(Mst1))
        - ((64091 + 48600*log(pow2(Msq)/pow2(Mst1)) - 148320*log(pow2(Mst2)/
        pow2(Mst1)))*pow2(Dmglst2))/(15.*pow2(Mgl)) + (360*pow2(Msq))/pow2(
        Mst1) + (75*pow2(Mst1)*(-8*(4*Dmglst2*Mgl + 10*pow2(Dmglst2) + pow2(
        Mgl))*pow2(Msq) + pow2(2*Dmglst2 + Mgl)*pow2(Mst1)))/(pow2(Mgl)*pow4(
        Msq)))*pow4(Mst2))/216. - ((5*(34*Dmglst2*Mgl + 91*pow2(Dmglst2) + 7*
        pow2(Mgl)))/(18.*pow2(Mgl)*pow2(Msq)) - (2.7222222222222223 + (46*
        Dmglst2)/(9.*Mgl) + (8*pow2(Dmglst2))/pow2(Mgl))/pow2(Mst1))*pow6(Mst2)
        + (5*(7*Mgl*(6*Dmglst2 + Mgl)*pow2(Mst1) + 3*pow2(Dmglst2)*(49*pow2(
        Mst1) - 30*pow2(Mst2)))*pow6(Mst2))/(36.*pow2(Mgl)*pow4(Msq))) + (pow2(
        Mt)*pow2(s2t)*(7560*log(pow2(Mst2)/pow2(Mst1))*pow2(Mst1)*(4*Dmglst2*
        Mgl*(124*pow2(Mst1)*pow2(Mst2) + 5*pow4(Mst1) + 297*pow4(Mst2)) + pow2(
        Mgl)*(68*pow2(Mst1)*pow2(Mst2) + 195*pow4(Mst1) + 305*pow4(Mst2)) +
        pow2(Dmglst2)*(2812*pow2(Mst1)*pow2(Mst2) + 20*pow4(Mst1) + 3290*pow4(
        Mst2))) + 7*pow2(Mgl)*(64800*pow2(Msq)*pow2(pow2(Mst1) - pow2(Mst2)) +
        1774344*pow2(Mst2)*pow4(Mst1) - 161307*pow2(Mst1)*pow4(Mst2) + 3321329*
        pow6(Mst1) - 105840*pow6(Mst2)) - 168*Dmglst2*Mgl*(99000*pow2(Mst2)*
        pow4(Mst1) + 20994*pow2(Mst1)*pow4(Mst2) + 392665*pow6(Mst1) + 8280*
        pow6(Mst2)) - 12*pow2(Dmglst2)*(3789645*pow2(Mst2)*pow4(Mst1) + 485343*
        pow2(Mst1)*pow4(Mst2) + 5497310*pow6(Mst1) + 181440*pow6(Mst2))))/(
        34020.*pow2(Mgl)*pow2(Mst1)*pow2(Mst2))) + pow4(Mt)*(402.71234567901234
         - (520*log(pow2(Msq)/pow2(Mst1)))/27. + (220*pow2(log(pow2(Msq)/pow2(
        Mst1))))/9. + (4*log(pow2(Mst2)/pow2(Mst1))*(-722 + 135*pow2(log(pow2(
        Msq)/pow2(Mst1)))))/81. + (4*(803 + 90*log(pow2(Msq)/pow2(Mst1)))*pow2(
        log(pow2(Mst2)/pow2(Mst1))))/27. - (40*pow3(log(pow2(Msq)/pow2(Mst1))))
        /9. - (148*pow3(log(pow2(Mst2)/pow2(Mst1))))/9. + (4*Dmglst2*(496858 +
        1050*log(pow2(Msq)/pow2(Mst1))*(23 - 6*log(pow2(Mst2)/pow2(Mst1))) +
        397915*log(pow2(Mst2)/pow2(Mst1)) + 28350*pow2(log(pow2(Msq)/pow2(Mst1)
        )) - 221970*pow2(log(pow2(Mst2)/pow2(Mst1))) + 40320*pow3(log(pow2(
        Mst2)/pow2(Mst1)))))/(2835.*Mgl) + pow4(Mst1)*((16.98404168016413 - (
        16403*log(pow2(Msq)/pow2(Mst1)))/2205. - (40*Dmglst2*(Dmglst2 + Mgl)*
        log(pow2(Mst2)/pow2(Mst1)))/(3.*pow2(Mgl)) + (38*pow2(log(pow2(Msq)/
        pow2(Mst1))))/21. - (10*pow2(log(pow2(Mst2)/pow2(Mst1))))/3.)/pow4(Msq)
        + (1405.4626938375586 - 140*log(pow2(Msq)/pow2(Mst1)) + (
        963.0333282942806 + 640*log(pow2(Msq)/pow2(Mst1)))*log(pow2(Mst2)/pow2(
        Mst1)) - (579.1343915343915 + 120*log(pow2(Msq)/pow2(Mst1)))*pow2(log(
        pow2(Mst2)/pow2(Mst1))) + (Dmglst2*(Dmglst2 + Mgl)*(-72755578289 -
        72968516460*log(pow2(Mst2)/pow2(Mst1)) + 25231726800*pow2(log(pow2(
        Mst2)/pow2(Mst1))) + 1500282000*log(pow2(Msq)/pow2(Mst1))*(23 - 38*log(
        pow2(Mst2)/pow2(Mst1)) + 6*pow2(log(pow2(Mst2)/pow2(Mst1)))) -
        3845167200*pow3(log(pow2(Mst2)/pow2(Mst1)))))/(1.8753525e7*pow2(Mgl)) +
        (4840*pow3(log(pow2(Mst2)/pow2(Mst1))))/27.)/pow4(Mst2)) - S2*(130.8 +
        84*log(pow2(Mst2)/pow2(Mst1)) + (4*Dmglst2*(-1593 + 980*log(pow2(Mst2)/
        pow2(Mst1))))/(35.*Mgl) + (158700*Dmglst2*Mgl*pow2(Mst1) + 524538*pow2(
        Mgl)*pow2(Mst1) + pow2(Dmglst2)*(257594*pow2(Mst1) + 487494*pow2(Mst2))
        - 420*log(pow2(Mst2)/pow2(Mst1))*(-1932*Dmglst2*Mgl*pow2(Mst1) - 627*
        pow2(Mgl)*pow2(Mst1) + pow2(Dmglst2)*(9406*pow2(Mst1) + 546*pow2(Mst2))
        ))/(945.*pow2(Mgl)*pow2(Mst2)) + (8*(363426 + 183960*log(pow2(Mst2)/
        pow2(Mst1)) + (Dmglst2*(Dmglst2 + Mgl)*(1000621 + 1172220*log(pow2(
        Mst2)/pow2(Mst1))))/pow2(Mgl))*pow4(Mst1))/(2835.*pow4(Mst2))) + ((
        pow2(Dmglst2)*(22336511 - 4210908*log(pow2(Mst2)/pow2(Mst1)) - 720*log(
        pow2(Msq)/pow2(Mst1))*(-61 + 660*log(pow2(Mst2)/pow2(Mst1))) + 364500*
        pow2(log(pow2(Msq)/pow2(Mst1))) - 113220*pow2(log(pow2(Mst2)/pow2(Mst1)
        )) + 518400*pow3(log(pow2(Mst2)/pow2(Mst1)))))/6075. - (64*(-4*pow2(
        Dmglst2) + pow2(Mgl) - 2*log(pow2(Mst2)/pow2(Mst1))*pow2(Dmglst2 + Mgl)
        + (4*Dmglst2*Mgl + 10*pow2(Dmglst2) + pow2(Mgl))*pow2(log(pow2(Mst2)/
        pow2(Mst1))))*pow4(Mst2))/(9.*pow4(Mst1)) + (16*OepS2*(2*pow2(Dmglst2)*
        (-14109*pow2(Mst1)*pow2(Mst2) + 11164*pow4(Mst1) - 819*pow4(Mst2)) + 3*
        pow2(Mgl)*(627*pow2(Mst1)*pow2(Mst2) + 1168*pow4(Mst1) + 189*pow4(Mst2)
        ) + 4*Dmglst2*Mgl*(1449*pow2(Mst1)*pow2(Mst2) + 5582*pow4(Mst1) + 189*
        pow4(Mst2))))/(2187.*pow4(Mst2)) + (276576612*Dmglst2*Mgl*pow2(Msq)*
        pow2(Mst1) + 7164806074*pow2(Dmglst2)*pow2(Msq)*pow2(Mst1) + 725226495*
        pow2(Mgl)*pow2(Msq)*pow2(Mst1) - 11340000*Dmglst2*Mgl*pow2(Mst1)*pow2(
        Mst2) - 51030000*pow2(Dmglst2)*pow2(Mst1)*pow2(Mst2) + 56657160*pow2(
        Mgl)*pow2(Mst1)*pow2(Mst2) + 378000*(478*Dmglst2*Mgl + 145*pow2(
        Dmglst2) + 145*pow2(Mgl))*pow2(Msq)*pow2(Mst1)*pow3(log(pow2(Mst2)/
        pow2(Mst1))) + 17010000*pow2(Mgl)*pow4(Msq) - 493290000*Dmglst2*Mgl*
        pow4(Mst1) - 493290000*pow2(Dmglst2)*pow4(Mst1) + 110565000*pow2(Mgl)*
        pow4(Mst1) - 28350000*Dmglst2*Mgl*pow4(Mst2) - 77414400*pow2(Dmglst2)*
        pow4(Mst2) - 59262840*pow2(Mgl)*pow4(Mst2) + 1134000*pow2(log(pow2(Msq)
        /pow2(Mst1)))*(60*Dmglst2*Mgl*pow4(Mst2) + 90*pow2(Dmglst2)*pow4(Mst2)
        + pow2(Mgl)*(7*pow2(Mst1)*pow2(Mst2) + 37*pow4(Mst2))) + 75600*log(
        pow2(Msq)/pow2(Mst1))*(-450*pow2(Dmglst2 - Mgl)*pow2(Msq)*pow2(Mst1)*
        pow2(log(pow2(Mst2)/pow2(Mst1))) + pow2(Mgl)*(-2700*pow2(Msq)*pow2(
        Mst1) - 41*pow2(Mst1)*pow2(Mst2) + 450*pow4(Msq) - 116*pow4(Mst2)) +
        450*Dmglst2*Mgl*(27*pow2(Msq)*pow2(Mst1) - 2*pow4(Mst2)) - 75*pow2(
        Dmglst2)*(195*pow2(Msq)*pow2(Mst1) + 38*pow4(Mst2)) + 15*log(pow2(Mst2)
        /pow2(Mst1))*(pow2(Mgl)*(225*pow2(Msq)*pow2(Mst1) - 19*pow4(Mst2)) +
        45*pow2(Dmglst2)*(9*pow2(Msq)*pow2(Mst1) - pow4(Mst2)) - 30*Dmglst2*
        Mgl*(19*pow2(Msq)*pow2(Mst1) + pow4(Mst2)))) + 37800*pow2(log(pow2(
        Mst2)/pow2(Mst1)))*(-10*pow2(Mgl)*(98*pow2(Msq)*pow2(Mst1) + 30*pow2(
        Mst1)*pow2(Mst2) - 45*pow4(Mst1) + 54*pow4(Mst2)) - 6*Dmglst2*Mgl*(
        3103*pow2(Msq)*pow2(Mst1) + 150*(pow4(Mst1) + pow4(Mst2))) + pow2(
        Dmglst2)*(31327*pow2(Msq)*pow2(Mst1) - 450*(2*pow4(Mst1) + 3*pow4(Mst2)
        ))) - 2520*log(pow2(Mst2)/pow2(Mst1))*(10*pow2(Mgl)*(-18799*pow2(Msq)*
        pow2(Mst1) + 225*pow2(Mst1)*pow2(Mst2) + 5400*pow4(Mst1) - 2598*pow4(
        Mst2)) + 6*Dmglst2*Mgl*(8581*pow2(Msq)*pow2(Mst1) - 1500*(-2*pow2(Mst1)
        *pow2(Mst2) + 15*pow4(Mst1) + 3*pow4(Mst2))) + pow2(Dmglst2)*(1633879*
        pow2(Msq)*pow2(Mst1) - 450*(-20*pow2(Mst1)*pow2(Mst2) + 300*pow4(Mst1)
        + 111*pow4(Mst2)))))/(1.27575e6*pow2(Msq)*pow2(Mst2)) - (300056400*
        Dmglst2*Mgl*pow2(Mst2)*pow4(Msq) + 630365400*pow2(Dmglst2)*pow2(Mst2)*
        pow4(Msq) - 171019800*pow2(Mgl)*pow2(Mst2)*pow4(Msq) - 85407000*
        Dmglst2*Mgl*pow2(Mst2)*pow4(Mst1) - 134284500*pow2(Dmglst2)*pow2(Mst2)*
        pow4(Mst1) - 2572500*pow2(Mgl)*pow2(Mst2)*pow4(Mst1) - 88200*pow2(Mst2)
        *pow2(log(pow2(Msq)/pow2(Mst1)))*(140*Dmglst2*Mgl*(13*pow2(Mst1)*pow2(
        Mst2) + 9*pow4(Msq) + 3*pow4(Mst1)) + 70*pow2(Dmglst2)*(59*pow2(Mst1)*
        pow2(Mst2) + 27*pow4(Msq) + 9*pow4(Mst1)) + pow2(Mgl)*(617*pow2(Mst1)*
        pow2(Mst2) + 630*pow4(Msq) + 210*pow4(Mst1))) + 88200*pow2(Mst2)*pow2(
        log(pow2(Mst2)/pow2(Mst1)))*(28*Dmglst2*Mgl*(25*pow2(Mst1)*pow2(Mst2) -
        104*pow4(Msq) + 15*pow4(Mst1)) + pow2(Mgl)*(223*pow2(Mst1)*pow2(Mst2) -
        1456*pow4(Msq) + 210*pow4(Mst1)) + pow2(Dmglst2)*(-4368*pow4(Msq) + 70*
        (25*pow2(Mst1)*pow2(Mst2) + 9*pow4(Mst1)))) + 230496000*Dmglst2*Mgl*
        pow2(Msq)*pow4(Mst2) + 329280000*pow2(Dmglst2)*pow2(Msq)*pow4(Mst2) +
        88494000*pow2(Mgl)*pow2(Msq)*pow4(Mst2) + 217290500*Dmglst2*Mgl*pow2(
        Mst1)*pow4(Mst2) + 309602090*pow2(Dmglst2)*pow2(Mst1)*pow4(Mst2) +
        142320737*pow2(Mgl)*pow2(Mst1)*pow4(Mst2) - 840*log(pow2(Mst2)/pow2(
        Mst1))*pow2(Mst2)*(-2450*Dmglst2*Mgl*(-120*pow2(Msq)*pow2(Mst2) - 136*
        pow2(Mst1)*pow2(Mst2) + 66*pow4(Msq) - 3*pow4(Mst1) - 63*pow4(Mst2)) +
        245*pow2(Dmglst2)*(3000*pow2(Msq)*pow2(Mst2) + 2938*pow2(Mst1)*pow2(
        Mst2) + 1506*pow4(Msq) - 315*pow4(Mst1) + 2205*pow4(Mst2)) + pow2(Mgl)*
        (73500*pow2(Msq)*pow2(Mst2) + 146492*pow2(Mst1)*pow2(Mst2) - 377790*
        pow4(Msq) + 47775*pow4(Mst1) + 25725*pow4(Mst2))) - 37044000*pow2(Mgl)*
        pow6(Msq) - 420*log(pow2(Msq)/pow2(Mst1))*(-420*log(pow2(Mst2)/pow2(
        Mst1))*pow2(Mst2)*(140*Dmglst2*Mgl*(4*pow2(Mst1)*pow2(Mst2) + 9*pow4(
        Msq)) + 70*pow2(Dmglst2)*(17*pow2(Mst1)*pow2(Mst2) + 27*pow4(Msq)) +
        pow2(Mgl)*(197*pow2(Mst1)*pow2(Mst2) + 630*pow4(Msq))) - 4900*Dmglst2*
        Mgl*pow2(Mst2)*(120*pow2(Msq)*pow2(Mst2) + 137*pow2(Mst1)*pow2(Mst2) +
        18*pow4(Msq) + 24*pow4(Mst1) + 63*pow4(Mst2)) - 2450*pow2(Dmglst2)*
        pow2(Mst2)*(600*pow2(Msq)*pow2(Mst2) + 763*pow2(Mst1)*pow2(Mst2) + 270*
        pow4(Msq) + 72*pow4(Mst1) + 441*pow4(Mst2)) + pow2(Mgl)*(220500*pow2(
        Mst2)*pow4(Msq) - 58800*pow2(Mst2)*pow4(Mst1) - 147000*pow2(Msq)*pow4(
        Mst2) - 220709*pow2(Mst1)*pow4(Mst2) + 176400*pow6(Msq) - 51450*pow6(
        Mst2))) + 60196500*Dmglst2*Mgl*pow6(Mst2) + 81033750*pow2(Dmglst2)*
        pow6(Mst2) + 17235750*pow2(Mgl)*pow6(Mst2))/(2.7783e6*pow2(Mst1)*pow4(
        Msq)))/pow2(Mgl)) - (s2t*pow3(Mt)*(-128297169*Dmglst2*Mgl*pow2(Mst2)*
        pow4(Msq)*pow4(Mst1) + 13940640*Dmglst2*Mgl*OepS2*pow2(Mst2)*pow4(Msq)*
        pow4(Mst1) - 106761564*Dmglst2*Mgl*S2*pow2(Mst2)*pow4(Msq)*pow4(Mst1) +
        272140503*pow2(Dmglst2)*pow2(Mst2)*pow4(Msq)*pow4(Mst1) - 14471520*
        OepS2*pow2(Dmglst2)*pow2(Mst2)*pow4(Msq)*pow4(Mst1) - 246467772*S2*
        pow2(Dmglst2)*pow2(Mst2)*pow4(Msq)*pow4(Mst1) + 66716559*pow2(Mgl)*
        pow2(Mst2)*pow4(Msq)*pow4(Mst1) - 3800160*OepS2*pow2(Mgl)*pow2(Mst2)*
        pow4(Msq)*pow4(Mst1) + 89139204*S2*pow2(Mgl)*pow2(Mst2)*pow4(Msq)*pow4(
        Mst1) - 1873935*Dmglst2*Mgl*pow2(Mst1)*pow4(Msq)*pow4(Mst2) + 1058400*
        Dmglst2*Mgl*OepS2*pow2(Mst1)*pow4(Msq)*pow4(Mst2) + 41334300*Dmglst2*
        Mgl*S2*pow2(Mst1)*pow4(Msq)*pow4(Mst2) + 330564249*pow2(Dmglst2)*pow2(
        Mst1)*pow4(Msq)*pow4(Mst2) - 70560*OepS2*pow2(Dmglst2)*pow2(Mst1)*pow4(
        Msq)*pow4(Mst2) - 27296676*S2*pow2(Dmglst2)*pow2(Mst1)*pow4(Msq)*pow4(
        Mst2) + 8935353*pow2(Mgl)*pow2(Mst1)*pow4(Msq)*pow4(Mst2) - 635040*
        OepS2*pow2(Mgl)*pow2(Mst1)*pow4(Msq)*pow4(Mst2) - 14512932*S2*pow2(Mgl)
        *pow2(Mst1)*pow4(Msq)*pow4(Mst2) + 680400*pow2(Mst1)*pow2(log(pow2(Msq)
        /pow2(Mst1)))*((pow2(Mst1) - pow2(Mst2))*(pow2(Mgl)*(4*pow2(Msq) +
        pow2(Mst1) + 3*pow2(Mst2)) + Dmglst2*Mgl*(4*pow2(Msq) + pow2(Mst1) + 7*
        pow2(Mst2)) + pow2(Dmglst2)*(4*pow2(Msq) + pow2(Mst1) + 13*pow2(Mst2)))
        - 6*log(pow2(Mst2)/pow2(Mst1))*(Dmglst2*Mgl + pow2(Dmglst2) + pow2(Mgl)
        )*pow4(Msq))*pow4(Mst2) + 9903600*Dmglst2*Mgl*pow2(Msq)*pow4(Mst1)*
        pow4(Mst2) + 378000*pow2(Dmglst2)*pow2(Msq)*pow4(Mst1)*pow4(Mst2) +
        7182000*pow2(Mgl)*pow2(Msq)*pow4(Mst1)*pow4(Mst2) - 1451520*pow2(Mst1)*
        pow3(log(pow2(Mst2)/pow2(Mst1)))*pow4(Msq)*(3*Dmglst2*Mgl*(-(pow2(Mst1)
        *pow2(Mst2)) + 3*pow4(Mst1) + 4*pow4(Mst2)) + 3*pow2(Dmglst2)*(pow2(
        Mst1)*pow2(Mst2) + 3*pow4(Mst1) + 4*pow4(Mst2)) + pow2(Mgl)*(9*pow2(
        Mst1)*pow2(Mst2) + 2*pow4(Mst1) + 12*pow4(Mst2))) - 20412000*Dmglst2*
        Mgl*pow2(Msq)*pow2(Mst2)*pow6(Mst1) - 20412000*pow2(Dmglst2)*pow2(Msq)*
        pow2(Mst2)*pow6(Mst1) + 9525600*pow2(Mgl)*pow2(Msq)*pow2(Mst2)*pow6(
        Mst1) - 677836516*Dmglst2*Mgl*pow4(Msq)*pow6(Mst1) + 57635200*Dmglst2*
        Mgl*OepS2*pow4(Msq)*pow6(Mst1) - 1331267616*Dmglst2*Mgl*S2*pow4(Msq)*
        pow6(Mst1) - 677836516*pow2(Dmglst2)*pow4(Msq)*pow6(Mst1) + 57635200*
        OepS2*pow2(Dmglst2)*pow4(Msq)*pow6(Mst1) - 1331267616*S2*pow2(Dmglst2)*
        pow4(Msq)*pow6(Mst1) + 160234788*pow2(Mgl)*pow4(Msq)*pow6(Mst1) -
        9206400*OepS2*pow2(Mgl)*pow4(Msq)*pow6(Mst1) + 295692768*S2*pow2(Mgl)*
        pow4(Msq)*pow6(Mst1) + 2929500*Dmglst2*Mgl*pow4(Mst2)*pow6(Mst1) +
        2929500*pow2(Dmglst2)*pow4(Mst2)*pow6(Mst1) + 888300*pow2(Mgl)*pow4(
        Mst2)*pow6(Mst1) + 56700000*Dmglst2*Mgl*pow2(Msq)*pow2(Mst1)*pow6(Mst2)
        + 117482400*pow2(Dmglst2)*pow2(Msq)*pow2(Mst1)*pow6(Mst2) + 17992800*
        pow2(Mgl)*pow2(Msq)*pow2(Mst1)*pow6(Mst2) - 4354560*Dmglst2*Mgl*pow4(
        Msq)*pow6(Mst2) - 17418240*pow2(Dmglst2)*pow4(Msq)*pow6(Mst2) +
        4354560*pow2(Mgl)*pow4(Msq)*pow6(Mst2) + 4649400*Dmglst2*Mgl*pow4(Mst1)
        *pow6(Mst2) + 9979200*pow2(Dmglst2)*pow4(Mst1)*pow6(Mst2) + 189000*
        pow2(Mgl)*pow4(Mst1)*pow6(Mst2) + 45360*pow2(log(pow2(Mst2)/pow2(Mst1))
        )*(pow2(Mgl)*(-60*pow2(Msq)*pow2(Mst1)*(pow2(Mst1) + 2*pow2(Mst2))*
        pow4(Mst2) - 15*pow2(Mst1)*pow4(Mst2)*(2*pow2(Mst1)*pow2(Mst2) + pow4(
        Mst1) + 3*pow4(Mst2)) + 2*pow4(Msq)*(432*pow2(Mst2)*pow4(Mst1) + 555*
        pow2(Mst1)*pow4(Mst2) + 70*pow6(Mst1) + 48*pow6(Mst2))) + Dmglst2*Mgl*(
        -60*pow2(Msq)*pow2(Mst1)*(pow2(Mst1) + 6*pow2(Mst2))*pow4(Mst2) - 15*
        pow2(Mst1)*pow4(Mst2)*(6*pow2(Mst1)*pow2(Mst2) + pow4(Mst1) + 15*pow4(
        Mst2)) + 2*pow4(Msq)*(-917*pow2(Mst2)*pow4(Mst1) - 349*pow2(Mst1)*pow4(
        Mst2) + 352*pow6(Mst1) + 144*pow6(Mst2))) + pow2(Dmglst2)*(-60*pow2(
        Msq)*pow2(Mst1)*(pow2(Mst1) + 12*pow2(Mst2))*pow4(Mst2) - 15*pow2(Mst1)
        *pow4(Mst2)*(12*pow2(Mst1)*pow2(Mst2) + pow4(Mst1) + 45*pow4(Mst2)) +
        pow4(Msq)*(295*pow2(Mst2)*pow4(Mst1) - 1650*pow2(Mst1)*pow4(Mst2) +
        704*pow6(Mst1) + 576*pow6(Mst2)))) + 42751800*Dmglst2*Mgl*pow2(Mst1)*
        pow8(Mst2) + 119070000*pow2(Dmglst2)*pow2(Mst1)*pow8(Mst2) + 11680200*
        pow2(Mgl)*pow2(Mst1)*pow8(Mst2) + 226800*log(pow2(Msq)/pow2(Mst1))*
        pow2(Mst1)*(pow2(Dmglst2)*(pow4(Msq)*(-72*pow2(Mst1)*pow2(Mst2) + 504*
        pow4(Mst1) - 131*pow4(Mst2)) + 14*pow2(Msq)*(pow2(Mst1) - pow2(Mst2))*
        pow4(Mst2) - 4*(3*pow2(Mst1)*pow2(Mst2) + pow4(Mst1) - 4*pow4(Mst2))*
        pow4(Mst2)) + 36*pow2(log(pow2(Mst2)/pow2(Mst1)))*pow4(Msq)*(pow2(Mgl)*
        (-pow4(Mst1) + pow4(Mst2)) + Dmglst2*Mgl*(3*pow4(Mst1) + pow4(Mst2)) +
        pow2(Dmglst2)*(3*pow4(Mst1) + pow4(Mst2))) - pow2(Mgl)*(14*pow2(Msq)*(-
        pow2(Mst1) + pow2(Mst2))*pow4(Mst2) + 4*pow4(Mst1)*pow4(Mst2) + 9*pow4(
        Msq)*(10*pow2(Mst1)*pow2(Mst2) + 4*pow4(Mst1) + 13*pow4(Mst2)) + 2*
        pow2(Mst1)*pow6(Mst2) - 6*pow8(Mst2)) + Dmglst2*Mgl*(pow4(Msq)*(234*
        pow2(Mst1)*pow2(Mst2) + 504*pow4(Mst1) - 161*pow4(Mst2)) + 14*pow2(Msq)
        *(pow2(Mst1) - pow2(Mst2))*pow4(Mst2) - 4*pow4(Mst1)*pow4(Mst2) - 6*
        pow2(Mst1)*pow6(Mst2) + 10*pow8(Mst2)) + 6*log(pow2(Mst2)/pow2(Mst1))*(
        3*pow2(Mgl)*(pow4(Msq)*(4*pow2(Mst1)*pow2(Mst2) + 11*pow4(Mst1) - 2*
        pow4(Mst2)) + 2*pow2(Msq)*pow6(Mst2) + pow8(Mst2)) + Dmglst2*Mgl*(pow4(
        Msq)*(-12*pow2(Mst1)*pow2(Mst2) - 123*pow4(Mst1) + 10*pow4(Mst2)) + 14*
        pow2(Msq)*pow6(Mst2) + 11*pow8(Mst2)) + pow2(Dmglst2)*(pow4(Msq)*(-123*
        pow4(Mst1) + 16*pow4(Mst2)) + 26*pow2(Msq)*pow6(Mst2) + 29*pow8(Mst2)))
        ) - 2520*log(pow2(Mst2)/pow2(Mst1))*(-9*pow2(Mgl)*(15*pow2(Mst1)*(6*
        pow2(Mst1)*pow2(Mst2) + 3*pow4(Mst1) - pow4(Mst2))*pow4(Mst2) + pow4(
        Msq)*(9*(102 + 377*S2)*pow2(Mst2)*pow4(Mst1) + 3*(116 + 189*S2)*pow2(
        Mst1)*pow4(Mst2) + 4*(302 + 2055*S2)*pow6(Mst1) - 384*pow6(Mst2)) - 60*
        pow2(Msq)*(-(pow4(Mst1)*pow4(Mst2)) + 4*pow2(Mst2)*pow6(Mst1) - 7*pow2(
        Mst1)*pow6(Mst2))) + pow2(Dmglst2)*(pow4(Msq)*(-27*(-2725 + 4307*S2)*
        pow2(Mst2)*pow4(Mst1) - 9*(5368 + 63*S2)*pow2(Mst1)*pow4(Mst2) + 4*(-
        3509 + 115785*S2)*pow6(Mst1) - 3456*pow6(Mst2)) - 180*pow2(Msq)*(-33*
        pow4(Mst1)*pow4(Mst2) + 12*pow2(Mst2)*pow6(Mst1) + 31*pow2(Mst1)*pow6(
        Mst2)) + 45*(15*pow4(Mst2)*pow6(Mst1) + 60*pow4(Mst1)*pow6(Mst2) - 31*
        pow2(Mst1)*pow8(Mst2))) + Dmglst2*Mgl*(pow4(Msq)*(3*(-18428 + 37341*S2)
        *pow2(Mst2)*pow4(Mst1) + 3*(-22768 + 2835*S2)*pow2(Mst1)*pow4(Mst2) +
        4*(-3509 + 115785*S2)*pow6(Mst1) + 3456*pow6(Mst2)) - 180*pow2(Msq)*(-
        21*pow4(Mst1)*pow4(Mst2) + 12*pow2(Mst2)*pow6(Mst1) + 37*pow2(Mst1)*
        pow6(Mst2)) + 45*(15*pow4(Mst2)*pow6(Mst1) - 6*pow4(Mst1)*pow6(Mst2) -
        25*pow2(Mst1)*pow8(Mst2))))))/(153090.*pow2(Mgl)*pow2(Mst1)*pow3(Mst2)*
        pow4(Msq))))/pow4(Mt)/pow2(Sbeta)*12.;

   return result;
}

/**
 * 	@return returns the susy log^1 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double himalaya::H6g2::calc_coef_at_as2_no_sm_logs_log1(){

   const double result =
      (-(pow2(Sbeta)*(5773600*Dmglst2*Mgl*pow4(Msq)*pow4(Mst1)*pow4(Mst2)*pow4(
        Mt) + 2073600*Dmglst2*Mgl*z2*pow4(Msq)*pow4(Mst1)*pow4(Mst2)*pow4(Mt) -
        7403968*pow2(Dmglst2)*pow4(Msq)*pow4(Mst1)*pow4(Mst2)*pow4(Mt) +
        2390400*z2*pow2(Dmglst2)*pow4(Msq)*pow4(Mst1)*pow4(Mst2)*pow4(Mt) +
        640400*pow2(Mgl)*pow4(Msq)*pow4(Mst1)*pow4(Mst2)*pow4(Mt) - 1936800*z2*
        pow2(Mgl)*pow4(Msq)*pow4(Mst1)*pow4(Mst2)*pow4(Mt) - 172800*z3*pow2(
        Mgl)*pow4(Msq)*pow4(Mst1)*pow4(Mst2)*pow4(Mt) - 14400*pow12(Mst2)*pow2(
        Dmglst2)*pow4(Msq)*pow4(s2t) + 450*pow2(log(pow2(Mst2)/pow2(Mst1)))*
        pow4(Msq)*pow4(Mst1)*(8*Dmglst2*Mgl*(8*Mt*(-67*Mst2*s2t*pow2(Mt) - 108*
        Mt*pow2(Mst2)*pow2(s2t) + 128*pow3(Mt) - 16*pow3(Mst2)*pow3(s2t))*pow4(
        Mst1) + pow4(Mst2)*(-320*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 848*Mst2*s2t*
        pow3(Mt) - 532*Mt*pow3(Mst2)*pow3(s2t) + 512*pow4(Mt) - 99*pow4(Mst2)*
        pow4(s2t)) + pow2(Mst1)*pow2(Mst2)*(-768*pow2(Mst2)*pow2(Mt)*pow2(s2t)
        - 448*Mst2*s2t*pow3(Mt) - 212*Mt*pow3(Mst2)*pow3(s2t) + 1656*pow4(Mt) +
        35*pow4(Mst2)*pow4(s2t))) + 4*pow2(Dmglst2)*(16*Mt*(-67*Mst2*s2t*pow2(
        Mt) - 108*Mt*pow2(Mst2)*pow2(s2t) + 128*pow3(Mt) - 16*pow3(Mst2)*pow3(
        s2t))*pow4(Mst1) + pow4(Mst2)*(-960*pow2(Mst2)*pow2(Mt)*pow2(s2t) +
        1696*Mst2*s2t*pow3(Mt) - 1832*Mt*pow3(Mst2)*pow3(s2t) + 1536*pow4(Mt) -
        297*pow4(Mst2)*pow4(s2t)) + pow2(Mst1)*pow2(Mst2)*(-2304*pow2(Mst2)*
        pow2(Mt)*pow2(s2t) + 256*Mst2*s2t*pow3(Mt) - 424*Mt*pow3(Mst2)*pow3(
        s2t) + 392*pow4(Mt) + 105*pow4(Mst2)*pow4(s2t))) + pow2(Mgl)*(pow4(
        Mst1)*(-2688*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 3136*Mst2*s2t*pow3(Mt) -
        1024*Mt*pow3(Mst2)*pow3(s2t) + 7168*pow4(Mt) - 41*pow4(Mst2)*pow4(s2t))
        - 8*pow2(Mst1)*pow2(Mst2)*(425*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 704*
        Mst2*s2t*pow3(Mt) + 212*Mt*pow3(Mst2)*pow3(s2t) - 196*pow4(Mt) + 3*
        pow4(Mst2)*pow4(s2t)) - pow4(Mst2)*(296*pow2(Mst2)*pow2(Mt)*pow2(s2t) -
        6784*Mst2*s2t*pow3(Mt) + 2208*Mt*pow3(Mst2)*pow3(s2t) + 288*pow4(Mt) +
        519*pow4(Mst2)*pow4(s2t)))) - 10246400*Dmglst2*Mgl*s2t*pow3(Mt)*pow4(
        Msq)*pow4(Mst1)*pow5(Mst2) + 2995200*Dmglst2*Mgl*s2t*z2*pow3(Mt)*pow4(
        Msq)*pow4(Mst1)*pow5(Mst2) - 7795200*s2t*pow2(Dmglst2)*pow3(Mt)*pow4(
        Msq)*pow4(Mst1)*pow5(Mst2) + 3542400*s2t*z2*pow2(Dmglst2)*pow3(Mt)*
        pow4(Msq)*pow4(Mst1)*pow5(Mst2) - 662400*s2t*pow2(Mgl)*pow3(Mt)*pow4(
        Msq)*pow4(Mst1)*pow5(Mst2) + 3398400*s2t*z2*pow2(Mgl)*pow3(Mt)*pow4(
        Msq)*pow4(Mst1)*pow5(Mst2) + 3600*log(pow2(Msq)/pow2(Mst1))*pow3(Mt)*
        pow4(Mst1)*(60*log(pow2(Mst2)/pow2(Mst1))*pow4(Msq)*(4*Dmglst2*Mgl*
        pow2(Mst1)*((9*Mt - 3*Mst2*s2t)*pow2(Mst1) + Mt*pow2(Mst2)) - 2*pow2(
        Dmglst2)*(Mt*pow2(Mst1)*pow2(Mst2) + 6*(-3*Mt + Mst2*s2t)*pow4(Mst1)) +
        pow2(Mgl)*(-2*Mt*pow2(Mst1)*pow2(Mst2) + (-9*Mt + 4*Mst2*s2t)*pow4(
        Mst1) + 2*Mt*pow4(Mst2))) + 5*pow2(Mgl)*(8*pow2(Msq)*((Mt - 2*Mst2*s2t)
        *pow2(Mst1) + 2*(2*Mt + Mst2*s2t)*pow2(Mst2))*pow4(Mst2) + 6*pow4(Msq)*
        (2*(5*Mt - 4*Mst2*s2t)*pow2(Mst1)*pow2(Mst2) + 2*(9*Mt - 4*Mst2*s2t)*
        pow4(Mst1) + (9*Mt + 4*Mst2*s2t)*pow4(Mst2)) + pow4(Mst2)*(pow2(Mst1)*(
        6*Mt*pow2(Mst2) - 8*s2t*pow3(Mst2)) + (3*Mt - 4*Mst2*s2t)*pow4(Mst1) +
        (19*Mt + 12*Mst2*s2t)*pow4(Mst2))) + 2*pow2(Dmglst2)*(pow4(Msq)*(330*
        Mt*pow2(Mst1)*pow2(Mst2) - 60*(27*Mt - 10*Mst2*s2t)*pow4(Mst1) + 4*(16*
        Mt - 5*Mst2*s2t)*pow4(Mst2)) + 20*pow2(Msq)*(9*Mst2*Mt - 2*s2t*pow2(
        Mst1) + 2*s2t*pow2(Mst2))*pow5(Mst2) + 5*(3*Mst2*(3*Mt - 8*Mst2*s2t)*
        pow2(Mst1) + (59*Mt + 26*Mst2*s2t)*pow3(Mst2) - 2*s2t*pow4(Mst1))*pow5(
        Mst2)) - 20*Dmglst2*Mgl*(2*pow4(Msq)*(3*(7*Mt - 2*Mst2*s2t)*pow2(Mst1)*
        pow2(Mst2) + (81*Mt - 30*Mst2*s2t)*pow4(Mst1) + (-5*Mt + Mst2*s2t)*
        pow4(Mst2)) - 4*pow2(Msq)*(3*Mst2*Mt - s2t*pow2(Mst1) + s2t*pow2(Mst2))
        *pow5(Mst2) + (3*Mst2*(-Mt + 2*Mst2*s2t)*pow2(Mst1) - (13*Mt + 7*Mst2*
        s2t)*pow3(Mst2) + s2t*pow4(Mst1))*pow5(Mst2))) + 432000*pow2(Mgl)*pow2(
        Mt)*pow2(s2t)*pow4(Mst1)*pow4(Mst2)*pow6(Msq) + 432000*pow2(Mgl)*pow2(
        Mst2)*pow4(Mst1)*pow4(Mt)*pow6(Msq) + 432000*pow2(Mgl)*pow2(Mst1)*pow4(
        Mst2)*pow4(Mt)*pow6(Msq) - 5241600*Dmglst2*Mgl*s2t*pow3(Mst2)*pow3(Mt)*
        pow4(Msq)*pow6(Mst1) - 4089600*Dmglst2*Mgl*s2t*z2*pow3(Mst2)*pow3(Mt)*
        pow4(Msq)*pow6(Mst1) + 10118400*s2t*pow2(Dmglst2)*pow3(Mst2)*pow3(Mt)*
        pow4(Msq)*pow6(Mst1) + 1296000*s2t*z2*pow2(Dmglst2)*pow3(Mst2)*pow3(Mt)
        *pow4(Msq)*pow6(Mst1) - 4723200*s2t*pow2(Mgl)*pow3(Mst2)*pow3(Mt)*pow4(
        Msq)*pow6(Mst1) + 3686400*s2t*z2*pow2(Mgl)*pow3(Mst2)*pow3(Mt)*pow4(
        Msq)*pow6(Mst1) - 4692000*Dmglst2*Mgl*pow2(Mt)*pow2(s2t)*pow4(Msq)*
        pow4(Mst2)*pow6(Mst1) + 2764800*Dmglst2*Mgl*z2*pow2(Mt)*pow2(s2t)*pow4(
        Msq)*pow4(Mst2)*pow6(Mst1) - 14175600*pow2(Dmglst2)*pow2(Mt)*pow2(s2t)*
        pow4(Msq)*pow4(Mst2)*pow6(Mst1) + 9676800*z2*pow2(Dmglst2)*pow2(Mt)*
        pow2(s2t)*pow4(Msq)*pow4(Mst2)*pow6(Mst1) - 410400*pow2(Mgl)*pow2(Mt)*
        pow2(s2t)*pow4(Msq)*pow4(Mst2)*pow6(Mst1) - 1228800*Dmglst2*Mgl*pow2(
        Mst2)*pow4(Msq)*pow4(Mt)*pow6(Mst1) + 12096000*Dmglst2*Mgl*z2*pow2(
        Mst2)*pow4(Msq)*pow4(Mt)*pow6(Mst1) - 30624800*pow2(Dmglst2)*pow2(Mst2)
        *pow4(Msq)*pow4(Mt)*pow6(Mst1) - 4032000*z2*pow2(Dmglst2)*pow2(Mst2)*
        pow4(Msq)*pow4(Mt)*pow6(Mst1) + 2793600*pow2(Mgl)*pow2(Mst2)*pow4(Msq)*
        pow4(Mt)*pow6(Mst1) - 2160000*z2*pow2(Mgl)*pow2(Mst2)*pow4(Msq)*pow4(
        Mt)*pow6(Mst1) - 288000*Dmglst2*Mgl*pow2(Msq)*pow4(Mst2)*pow4(Mt)*pow6(
        Mst1) - 144000*pow2(Dmglst2)*pow2(Msq)*pow4(Mst2)*pow4(Mt)*pow6(Mst1) -
        24000*pow2(Mgl)*pow2(Msq)*pow4(Mst2)*pow4(Mt)*pow6(Mst1) + 336000*
        Dmglst2*Mgl*s2t*pow2(Msq)*pow3(Mt)*pow5(Mst2)*pow6(Mst1) + 624000*s2t*
        pow2(Dmglst2)*pow2(Msq)*pow3(Mt)*pow5(Mst2)*pow6(Mst1) - 240000*s2t*
        pow2(Mgl)*pow2(Msq)*pow3(Mt)*pow5(Mst2)*pow6(Mst1) + 1598400*Dmglst2*
        Mgl*Mt*pow3(s2t)*pow4(Msq)*pow5(Mst2)*pow6(Mst1) - 244800*Dmglst2*Mgl*
        Mt*z2*pow3(s2t)*pow4(Msq)*pow5(Mst2)*pow6(Mst1) + 199200*Mt*pow2(
        Dmglst2)*pow3(s2t)*pow4(Msq)*pow5(Mst2)*pow6(Mst1) - 244800*Mt*z2*pow2(
        Dmglst2)*pow3(s2t)*pow4(Msq)*pow5(Mst2)*pow6(Mst1) + 2131200*Mt*pow2(
        Mgl)*pow3(s2t)*pow4(Msq)*pow5(Mst2)*pow6(Mst1) - 244800*Mt*z2*pow2(Mgl)
        *pow3(s2t)*pow4(Msq)*pow5(Mst2)*pow6(Mst1) - 216000*pow2(Mgl)*pow2(
        Mst2)*pow2(Mt)*pow2(s2t)*pow6(Msq)*pow6(Mst1) - 27000*pow2(Mgl)*pow4(
        Mst2)*pow4(s2t)*pow6(Msq)*pow6(Mst1) - 8556000*Dmglst2*Mgl*pow2(Mt)*
        pow2(s2t)*pow4(Msq)*pow4(Mst1)*pow6(Mst2) + 4147200*Dmglst2*Mgl*z2*
        pow2(Mt)*pow2(s2t)*pow4(Msq)*pow4(Mst1)*pow6(Mst2) - 17101200*pow2(
        Dmglst2)*pow2(Mt)*pow2(s2t)*pow4(Msq)*pow4(Mst1)*pow6(Mst2) + 11750400*
        z2*pow2(Dmglst2)*pow2(Mt)*pow2(s2t)*pow4(Msq)*pow4(Mst1)*pow6(Mst2) -
        3636000*pow2(Mgl)*pow2(Mt)*pow2(s2t)*pow4(Msq)*pow4(Mst1)*pow6(Mst2) +
        691200*z2*pow2(Mgl)*pow2(Mt)*pow2(s2t)*pow4(Msq)*pow4(Mst1)*pow6(Mst2)
        - 1008000*Dmglst2*Mgl*pow2(Mst1)*pow4(Msq)*pow4(Mt)*pow6(Mst2) +
        648000*pow2(Dmglst2)*pow2(Mst1)*pow4(Msq)*pow4(Mt)*pow6(Mst2) -
        1540800*pow2(Mgl)*pow2(Mst1)*pow4(Msq)*pow4(Mt)*pow6(Mst2) - 72000*
        Dmglst2*Mgl*pow2(Msq)*pow4(Mst1)*pow4(Mt)*pow6(Mst2) - 813600*pow2(
        Dmglst2)*pow2(Msq)*pow4(Mst1)*pow4(Mt)*pow6(Mst2) + 468000*pow2(Mgl)*
        pow2(Msq)*pow4(Mst1)*pow4(Mt)*pow6(Mst2) - 216000*pow2(Mgl)*pow2(Mst1)*
        pow2(Mt)*pow2(s2t)*pow6(Msq)*pow6(Mst2) - 27000*pow2(Mgl)*pow4(Mst1)*
        pow4(s2t)*pow6(Msq)*pow6(Mst2) - 126000*Dmglst2*Mgl*pow4(Mt)*pow6(Mst1)
        *pow6(Mst2) - 405000*pow2(Dmglst2)*pow4(Mt)*pow6(Mst1)*pow6(Mst2) +
        45000*pow2(Mgl)*pow4(Mt)*pow6(Mst1)*pow6(Mst2) - 1800*Dmglst2*Mgl*pow4(
        Msq)*pow4(s2t)*pow6(Mst1)*pow6(Mst2) + 154800*Dmglst2*Mgl*z2*pow4(Msq)*
        pow4(s2t)*pow6(Mst1)*pow6(Mst2) - 124500*pow2(Dmglst2)*pow4(Msq)*pow4(
        s2t)*pow6(Mst1)*pow6(Mst2) + 232200*z2*pow2(Dmglst2)*pow4(Msq)*pow4(
        s2t)*pow6(Mst1)*pow6(Mst2) + 522900*pow2(Mgl)*pow4(Msq)*pow4(s2t)*pow6(
        Mst1)*pow6(Mst2) + 77400*z2*pow2(Mgl)*pow4(Msq)*pow4(s2t)*pow6(Mst1)*
        pow6(Mst2) + 921600*Dmglst2*Mgl*s2t*pow2(Mst1)*pow3(Mt)*pow4(Msq)*pow7(
        Mst2) - 921600*s2t*pow2(Dmglst2)*pow2(Mst1)*pow3(Mt)*pow4(Msq)*pow7(
        Mst2) + 921600*s2t*pow2(Mgl)*pow2(Mst1)*pow3(Mt)*pow4(Msq)*pow7(Mst2) -
        288000*Dmglst2*Mgl*s2t*pow2(Msq)*pow3(Mt)*pow4(Mst1)*pow7(Mst2) -
        576000*s2t*pow2(Dmglst2)*pow2(Msq)*pow3(Mt)*pow4(Mst1)*pow7(Mst2) +
        96000*s2t*pow2(Mgl)*pow2(Msq)*pow3(Mt)*pow4(Mst1)*pow7(Mst2) - 2764800*
        Dmglst2*Mgl*Mt*pow3(s2t)*pow4(Msq)*pow4(Mst1)*pow7(Mst2) + 1310400*
        Dmglst2*Mgl*Mt*z2*pow3(s2t)*pow4(Msq)*pow4(Mst1)*pow7(Mst2) - 4864800*
        Mt*pow2(Dmglst2)*pow3(s2t)*pow4(Msq)*pow4(Mst1)*pow7(Mst2) + 2865600*
        Mt*z2*pow2(Dmglst2)*pow3(s2t)*pow4(Msq)*pow4(Mst1)*pow7(Mst2) -
        2174400*Mt*pow2(Mgl)*pow3(s2t)*pow4(Msq)*pow4(Mst1)*pow7(Mst2) +
        273600*Mt*z2*pow2(Mgl)*pow3(s2t)*pow4(Msq)*pow4(Mst1)*pow7(Mst2) +
        36000*Dmglst2*Mgl*s2t*pow3(Mt)*pow6(Mst1)*pow7(Mst2) + 504000*s2t*pow2(
        Dmglst2)*pow3(Mt)*pow6(Mst1)*pow7(Mst2) - 84000*s2t*pow2(Mgl)*pow3(Mt)*
        pow6(Mst1)*pow7(Mst2) - 288000*Dmglst2*Mgl*s2t*pow2(Msq)*pow3(Mst2)*
        pow3(Mt)*pow8(Mst1) - 288000*s2t*pow2(Dmglst2)*pow2(Msq)*pow3(Mst2)*
        pow3(Mt)*pow8(Mst1) + 288000*s2t*pow2(Mgl)*pow2(Msq)*pow3(Mst2)*pow3(
        Mt)*pow8(Mst1) - 4435200*Dmglst2*Mgl*pow2(Mst2)*pow2(Mt)*pow2(s2t)*
        pow4(Msq)*pow8(Mst1) + 2764800*Dmglst2*Mgl*z2*pow2(Mst2)*pow2(Mt)*pow2(
        s2t)*pow4(Msq)*pow8(Mst1) - 4435200*pow2(Dmglst2)*pow2(Mst2)*pow2(Mt)*
        pow2(s2t)*pow4(Msq)*pow8(Mst1) + 2764800*z2*pow2(Dmglst2)*pow2(Mst2)*
        pow2(Mt)*pow2(s2t)*pow4(Msq)*pow8(Mst1) - 194400*pow2(Mgl)*pow2(Mst2)*
        pow2(Mt)*pow2(s2t)*pow4(Msq)*pow8(Mst1) + 7084800*Dmglst2*Mgl*Mst2*s2t*
        pow3(Mt)*pow4(Msq)*pow8(Mst1) - 11116800*Dmglst2*Mgl*Mst2*s2t*z2*pow3(
        Mt)*pow4(Msq)*pow8(Mst1) + 7084800*Mst2*s2t*pow2(Dmglst2)*pow3(Mt)*
        pow4(Msq)*pow8(Mst1) - 11116800*Mst2*s2t*z2*pow2(Dmglst2)*pow3(Mt)*
        pow4(Msq)*pow8(Mst1) - 5860800*Mst2*s2t*pow2(Mgl)*pow3(Mt)*pow4(Msq)*
        pow8(Mst1) + 3571200*Mst2*s2t*z2*pow2(Mgl)*pow3(Mt)*pow4(Msq)*pow8(
        Mst1) - 82800*Dmglst2*Mgl*Mt*pow3(Mst2)*pow3(s2t)*pow4(Msq)*pow8(Mst1)
        - 82800*Mt*pow2(Dmglst2)*pow3(Mst2)*pow3(s2t)*pow4(Msq)*pow8(Mst1) +
        176400*Mt*pow2(Mgl)*pow3(Mst2)*pow3(s2t)*pow4(Msq)*pow8(Mst1) +
        1440000*Dmglst2*Mgl*pow2(Msq)*pow2(Mst2)*pow4(Mt)*pow8(Mst1) + 1440000*
        pow2(Dmglst2)*pow2(Msq)*pow2(Mst2)*pow4(Mt)*pow8(Mst1) - 504000*pow2(
        Mgl)*pow2(Msq)*pow2(Mst2)*pow4(Mt)*pow8(Mst1) - 23808000*Dmglst2*Mgl*
        pow4(Msq)*pow4(Mt)*pow8(Mst1) + 27302400*Dmglst2*Mgl*z2*pow4(Msq)*pow4(
        Mt)*pow8(Mst1) - 23808000*pow2(Dmglst2)*pow4(Msq)*pow4(Mt)*pow8(Mst1) +
        27302400*z2*pow2(Dmglst2)*pow4(Msq)*pow4(Mt)*pow8(Mst1) + 4991400*pow2(
        Mgl)*pow4(Msq)*pow4(Mt)*pow8(Mst1) - 2534400*z2*pow2(Mgl)*pow4(Msq)*
        pow4(Mt)*pow8(Mst1) - 108000*Dmglst2*Mgl*pow4(Mst2)*pow4(Mt)*pow8(Mst1)
        - 108000*pow2(Dmglst2)*pow4(Mst2)*pow4(Mt)*pow8(Mst1) - 49500*pow2(Mgl)
        *pow4(Mst2)*pow4(Mt)*pow8(Mst1) - 132300*Dmglst2*Mgl*pow4(Msq)*pow4(
        Mst2)*pow4(s2t)*pow8(Mst1) - 132300*pow2(Dmglst2)*pow4(Msq)*pow4(Mst2)*
        pow4(s2t)*pow8(Mst1) - 375075*pow2(Mgl)*pow4(Msq)*pow4(Mst2)*pow4(s2t)*
        pow8(Mst1) - 36900*z2*pow2(Mgl)*pow4(Msq)*pow4(Mst2)*pow4(s2t)*pow8(
        Mst1) + 138000*Dmglst2*Mgl*s2t*pow3(Mt)*pow5(Mst2)*pow8(Mst1) + 138000*
        s2t*pow2(Dmglst2)*pow3(Mt)*pow5(Mst2)*pow8(Mst1) - 6000*s2t*pow2(Mgl)*
        pow3(Mt)*pow5(Mst2)*pow8(Mst1) + 27000*pow2(Mgl)*pow2(Mst2)*pow4(s2t)*
        pow6(Msq)*pow8(Mst1) + 964800*Dmglst2*Mgl*pow2(Mst1)*pow2(Mt)*pow2(s2t)
        *pow4(Msq)*pow8(Mst2) + 367200*pow2(Dmglst2)*pow2(Mst1)*pow2(Mt)*pow2(
        s2t)*pow4(Msq)*pow8(Mst2) + 885600*pow2(Mgl)*pow2(Mst1)*pow2(Mt)*pow2(
        s2t)*pow4(Msq)*pow8(Mst2) + 460800*Dmglst2*Mgl*pow4(Msq)*pow4(Mt)*pow8(
        Mst2) + 230400*pow2(Dmglst2)*pow4(Msq)*pow4(Mt)*pow8(Mst2) + 230400*
        pow2(Mgl)*pow4(Msq)*pow4(Mt)*pow8(Mst2) - 78000*Dmglst2*Mgl*pow4(Mst1)*
        pow4(Mt)*pow8(Mst2) - 857400*pow2(Dmglst2)*pow4(Mst1)*pow4(Mt)*pow8(
        Mst2) + 115500*pow2(Mgl)*pow4(Mst1)*pow4(Mt)*pow8(Mst2) + 379800*
        Dmglst2*Mgl*pow4(Msq)*pow4(Mst1)*pow4(s2t)*pow8(Mst2) - 154800*Dmglst2*
        Mgl*z2*pow4(Msq)*pow4(Mst1)*pow4(s2t)*pow8(Mst2) + 502800*pow2(Dmglst2)
        *pow4(Msq)*pow4(Mst1)*pow4(s2t)*pow8(Mst2) - 232200*z2*pow2(Dmglst2)*
        pow4(Msq)*pow4(Mst1)*pow4(s2t)*pow8(Mst2) - 37800*pow2(Mgl)*pow4(Msq)*
        pow4(Mst1)*pow4(s2t)*pow8(Mst2) - 40500*z2*pow2(Mgl)*pow4(Msq)*pow4(
        Mst1)*pow4(s2t)*pow8(Mst2) + 27000*pow2(Mgl)*pow2(Mst1)*pow4(s2t)*pow6(
        Msq)*pow8(Mst2) - 30*log(pow2(Mst2)/pow2(Mst1))*(20*Dmglst2*Mgl*(240*
        pow2(Msq)*pow2(Mst2)*pow3(Mt)*pow4(Mst1)*(-2*s2t*pow2(Mst1)*pow3(Mst2)
        + 3*Mt*pow4(Mst1) + 2*(3*Mt + Mst2*s2t)*pow4(Mst2)) - 120*pow3(Mt)*
        pow4(Mst1)*(3*Mst2*(-Mt + 2*Mst2*s2t)*pow2(Mst1) - (13*Mt + 7*Mst2*s2t)
        *pow3(Mst2) + s2t*pow4(Mst1))*pow5(Mst2) + pow4(Msq)*(-4*pow4(Mst1)*
        pow4(Mst2)*(1860*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 1532*Mst2*s2t*pow3(Mt)
        + 774*Mt*pow3(Mst2)*pow3(s2t) - 2764*pow4(Mt) + 15*pow4(Mst2)*pow4(s2t)
        ) + 2*pow2(Mst2)*(-432*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 10936*Mst2*s2t*
        pow3(Mt) + 996*Mt*pow3(Mst2)*pow3(s2t) + 16228*pow4(Mt) + 207*pow4(
        Mst2)*pow4(s2t))*pow6(Mst1) - 3*pow2(Mst1)*(-856*pow2(Mst2)*pow2(Mt)*
        pow2(s2t) - 1536*Mst2*s2t*pow3(Mt) + 384*Mt*pow3(Mst2)*pow3(s2t) +
        1200*pow4(Mt) + 203*pow4(Mst2)*pow4(s2t))*pow6(Mst2) + 4*(-738*pow2(
        Mst2)*pow2(Mt)*pow2(s2t) - 8204*Mst2*s2t*pow3(Mt) + 135*Mt*pow3(Mst2)*
        pow3(s2t) + 12492*pow4(Mt) + 12*pow4(Mst2)*pow4(s2t))*pow8(Mst1) - 768*
        pow2(Mt)*(-2*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow8(Mst2))) - 5*pow2(
        Mgl)*(480*pow2(Msq)*pow2(Mst2)*pow3(Mt)*pow4(Mst1)*(2*(-Mt + 2*Mst2*
        s2t)*pow2(Mst1)*pow2(Mst2) + 3*Mt*pow4(Mst1) - 4*(2*Mt + Mst2*s2t)*
        pow4(Mst2)) + 120*pow3(Mt)*pow4(Mst1)*pow4(Mst2)*(2*(-3*Mt + 4*Mst2*
        s2t)*pow2(Mst1)*pow2(Mst2) + (-3*Mt + 4*Mst2*s2t)*pow4(Mst1) - (19*Mt +
        12*Mst2*s2t)*pow4(Mst2)) + 360*(pow2(Mst1) - pow2(Mst2))*pow4(Mst1)*
        pow4(Mst2)*pow4(s2t)*pow6(Msq) + pow4(Msq)*(4*pow4(Mst1)*pow4(Mst2)*(
        3372*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 9744*Mst2*s2t*pow3(Mt) + 3912*Mt*
        pow3(Mst2)*pow3(s2t) + 4492*pow4(Mt) + 669*pow4(Mst2)*pow4(s2t)) - 48*
        pow2(Mst2)*(-267*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 792*Mst2*s2t*pow3(Mt)
        - 30*Mt*pow3(Mst2)*pow3(s2t) - 317*pow4(Mt) + 56*pow4(Mst2)*pow4(s2t))*
        pow6(Mst1) + 6*pow2(Mst1)*(-728*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 1024*
        Mst2*s2t*pow3(Mt) + 256*Mt*pow3(Mst2)*pow3(s2t) + 1200*pow4(Mt) + 139*
        pow4(Mst2)*pow4(s2t))*pow6(Mst2) + 3*(4640*pow2(Mst2)*pow2(Mt)*pow2(
        s2t) - 13376*Mst2*s2t*pow3(Mt) + 240*Mt*pow3(Mst2)*pow3(s2t) + 2328*
        pow4(Mt) - 279*pow4(Mst2)*pow4(s2t))*pow8(Mst1) + 768*pow2(Mt)*(-2*
        pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow8(Mst2))) + 2*pow2(Dmglst2)*(2400*
        pow2(Msq)*pow2(Mst2)*pow3(Mt)*pow4(Mst1)*(-2*s2t*pow2(Mst1)*pow3(Mst2)
        + 3*Mt*pow4(Mst1) + (9*Mt + 2*Mst2*s2t)*pow4(Mst2)) - 600*pow3(Mt)*
        pow4(Mst1)*(3*Mst2*(-3*Mt + 8*Mst2*s2t)*pow2(Mst1) - (59*Mt + 26*Mst2*
        s2t)*pow3(Mst2) + 2*s2t*pow4(Mst1))*pow5(Mst2) + pow4(Msq)*(2*pow4(
        Mst1)*pow4(Mst2)*(-40440*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 76560*Mst2*
        s2t*pow3(Mt) - 5760*Mt*pow3(Mst2)*pow3(s2t) + 15712*pow4(Mt) + 4335*
        pow4(Mst2)*pow4(s2t)) + 2*pow2(Mst2)*(29760*pow2(Mst2)*pow2(Mt)*pow2(
        s2t) + 31000*Mst2*s2t*pow3(Mt) + 15840*Mt*pow3(Mst2)*pow3(s2t) -
        208652*pow4(Mt) + 1065*pow4(Mst2)*pow4(s2t))*pow6(Mst1) - 15*pow2(Mst1)
        *(-3080*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 6144*Mst2*s2t*pow3(Mt) + 1536*
        Mt*pow3(Mst2)*pow3(s2t) + 3600*pow4(Mt) + 865*pow4(Mst2)*pow4(s2t))*
        pow6(Mst2) + 40*(-738*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 8204*Mst2*s2t*
        pow3(Mt) + 135*Mt*pow3(Mst2)*pow3(s2t) + 12492*pow4(Mt) + 12*pow4(Mst2)
        *pow4(s2t))*pow8(Mst1) + 480*(-40*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 80*
        pow4(Mt) + 3*pow4(Mst2)*pow4(s2t))*pow8(Mst2)))) - 230400*Dmglst2*Mgl*
        Mt*pow2(Mst1)*pow3(s2t)*pow4(Msq)*pow9(Mst2) + 230400*Mt*pow2(Dmglst2)*
        pow2(Mst1)*pow3(s2t)*pow4(Msq)*pow9(Mst2) - 230400*Mt*pow2(Mgl)*pow2(
        Mst1)*pow3(s2t)*pow4(Msq)*pow9(Mst2) - 162000*Dmglst2*Mgl*s2t*pow3(Mt)*
        pow4(Mst1)*pow9(Mst2) - 702000*s2t*pow2(Dmglst2)*pow3(Mt)*pow4(Mst1)*
        pow9(Mst2) + 54000*s2t*pow2(Mgl)*pow3(Mt)*pow4(Mst1)*pow9(Mst2) -
        230400*Dmglst2*Mgl*pow2(Mt)*pow2(s2t)*pow4(Msq)*power10(Mst2) - 115200*
        pow2(Dmglst2)*pow2(Mt)*pow2(s2t)*pow4(Msq)*power10(Mst2) - 115200*pow2(
        Mgl)*pow2(Mt)*pow2(s2t)*pow4(Msq)*power10(Mst2) - 235800*Dmglst2*Mgl*
        pow2(Mst1)*pow4(Msq)*pow4(s2t)*power10(Mst2) - 161100*pow2(Dmglst2)*
        pow2(Mst1)*pow4(Msq)*pow4(s2t)*power10(Mst2) - 153900*pow2(Mgl)*pow2(
        Mst1)*pow4(Msq)*pow4(s2t)*power10(Mst2)))/(16200.*pow2(Mgl)*pow4(Msq)*
        pow4(Mst1)*pow4(Mst2)))/pow4(Mt)/pow2(Sbeta)*12.; 

   return result;
}

/**
 * 	@return returns the susy log^2 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double himalaya::H6g2::calc_coef_at_as2_no_sm_logs_log2(){

   const double result =
      ((pow2(Sbeta)*(-48160*Dmglst2*Mgl*pow4(Mst1)*pow4(Mst2)*pow4(Mt) - 29264*
        pow2(Dmglst2)*pow4(Mst1)*pow4(Mst2)*pow4(Mt) + 21280*pow2(Mgl)*pow4(
        Mst1)*pow4(Mst2)*pow4(Mt) + 7200*log(pow2(Msq)/pow2(Mst1))*pow2(Mgl)*
        pow4(Mst1)*pow4(Mst2)*pow4(Mt) - 1440*pow12(Mst2)*pow2(Dmglst2)*pow4(
        s2t) + 46400*Dmglst2*Mgl*s2t*pow3(Mt)*pow4(Mst1)*pow5(Mst2) + 98880*
        s2t*pow2(Dmglst2)*pow3(Mt)*pow4(Mst1)*pow5(Mst2) - 24000*s2t*pow2(Mgl)*
        pow3(Mt)*pow4(Mst1)*pow5(Mst2) + 74880*Dmglst2*Mgl*s2t*pow3(Mst2)*pow3(
        Mt)*pow6(Mst1) - 38400*s2t*pow2(Dmglst2)*pow3(Mst2)*pow3(Mt)*pow6(Mst1)
        + 1920*s2t*pow2(Mgl)*pow3(Mst2)*pow3(Mt)*pow6(Mst1) - 14160*Dmglst2*
        Mgl*pow2(Mt)*pow2(s2t)*pow4(Mst2)*pow6(Mst1) - 13560*pow2(Dmglst2)*
        pow2(Mt)*pow2(s2t)*pow4(Mst2)*pow6(Mst1) - 9000*pow2(Mgl)*pow2(Mt)*
        pow2(s2t)*pow4(Mst2)*pow6(Mst1) - 100800*Dmglst2*Mgl*pow2(Mst2)*pow4(
        Mt)*pow6(Mst1) + 161120*pow2(Dmglst2)*pow2(Mst2)*pow4(Mt)*pow6(Mst1) -
        10080*pow2(Mgl)*pow2(Mst2)*pow4(Mt)*pow6(Mst1) - 31200*Dmglst2*Mgl*Mt*
        pow3(s2t)*pow5(Mst2)*pow6(Mst1) - 42720*Mt*pow2(Dmglst2)*pow3(s2t)*
        pow5(Mst2)*pow6(Mst1) - 23520*Mt*pow2(Mgl)*pow3(s2t)*pow5(Mst2)*pow6(
        Mst1) + 82080*Dmglst2*Mgl*pow2(Mt)*pow2(s2t)*pow4(Mst1)*pow6(Mst2) +
        115440*pow2(Dmglst2)*pow2(Mt)*pow2(s2t)*pow4(Mst1)*pow6(Mst2) + 42960*
        pow2(Mgl)*pow2(Mt)*pow2(s2t)*pow4(Mst1)*pow6(Mst2) + 43680*Dmglst2*Mgl*
        pow2(Mst1)*pow4(Mt)*pow6(Mst2) + 65520*pow2(Dmglst2)*pow2(Mst1)*pow4(
        Mt)*pow6(Mst2) + 21840*pow2(Mgl)*pow2(Mst1)*pow4(Mt)*pow6(Mst2) + 1110*
        Dmglst2*Mgl*pow4(s2t)*pow6(Mst1)*pow6(Mst2) + 5505*pow2(Dmglst2)*pow4(
        s2t)*pow6(Mst1)*pow6(Mst2) - 5325*pow2(Mgl)*pow4(s2t)*pow6(Mst1)*pow6(
        Mst2) + 30*log(pow2(Mst2)/pow2(Mst1))*pow4(Mst1)*(pow2(Mgl)*(-4*pow4(
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
        *pow2(Mst1)*pow4(s2t)*power10(Mst2)))/(540.*pow2(Mgl)*pow4(Mst1)*pow4(
        Mst2)))/pow4(Mt)/pow2(Sbeta)*12.; 

   return result;
}

/**
 * 	@return returns the susy log^3 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double himalaya::H6g2::calc_coef_at_as2_no_sm_logs_log3(){

   const double result =
      ((-224*pow2(Sbeta)*pow4(Mt))/9.)/pow4(Mt)/pow2(Sbeta)*12.; 

   return result;
}

