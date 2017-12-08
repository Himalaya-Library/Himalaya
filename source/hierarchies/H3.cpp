#include "H3.hpp"
#include "HierarchyCalculator.hpp"
#include "Constants.hpp"
#include "Utils.hpp"
#include <iomanip>      // std::setprecision
#include <cmath>
#include <type_traits>

/**
 * 	Constuctor
 * 	@param flagMap the flagMap for the truncation of expansion variables
 * 	@param Al4p a double alpha_s/4/Pi
 * 	@param beta a double which is the mixing angle beta
 * 	@param Dmglst1 a double Mgl - Mst1
 * 	@param Dmst12 a double Mst1^2 - Mst2^2
 * 	@param Dmsqst1 a double Msq^2 - Mst1^2
 * 	@param lmMt a double log((renormalization scale / Mt)^2)
 * 	@param lmMst1 a double log((renormalization scale / Mst1)^2)
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
himalaya::H3::H3 (const std::map<unsigned int, unsigned int>& flagMap, double Al4p, double beta,
                   double Dmglst1, double Dmst12, double Dmsqst1, double lmMt, double lmMst1,
                   double Mgl, double Mt, double Mst1, double Mst2, double Msq, double MuSUSY,
                   double s2t,
                   int mdrFlag, int oneLoopFlag, int twoLoopFlag, int threeLoopFlag) {
   // abbrev for tan(beta) and sin(beta)
   Tbeta = tan(beta);
   Sbeta = sin(beta);
   this -> Dmglst1 = Dmglst1;
   this -> Dmst12 = Dmst12;
   this -> Dmsqst1 = Dmsqst1;
   this -> lmMst1 = lmMst1;
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
   // expansion flags
   xDmglst1 = flagMap.at (HierarchyCalculator::xxDmglst1);
   xDmst12 = flagMap.at (HierarchyCalculator::xxDmglst1);
   xDmsqst1 = flagMap.at (HierarchyCalculator::xxDmsqst1);
   
   s1 = oneLoopFlag * calc_coeff_as_0_log_0_s1() 
      + twoLoopFlag * Al4p * calc_coeff_as_1_log_0_s1() 
      + threeLoopFlag * pow2(Al4p) * (calc_coeff_as_2_log_0_s1() 
	 + lmMt * calc_coeff_as_2_log_1_s1());
   s2 = oneLoopFlag * (calc_coeff_as_0_log_0_s2() + lmMt * calc_coeff_as_0_log_1_s2())
      + twoLoopFlag * Al4p * (calc_coeff_as_1_log_0_s2() + lmMt * calc_coeff_as_1_log_1_s2()
	+ pow2(lmMt) * calc_coeff_as_1_log_2_s2())
      + threeLoopFlag * pow2(Al4p) * (calc_coeff_as_2_log_0_s2() + lmMt * calc_coeff_as_2_log_1_s2()
	 + pow2(lmMt) * calc_coeff_as_2_log_2_s2() + pow3(lmMt) * calc_coeff_as_2_log_3_s2());
   s12 = oneLoopFlag * calc_coeff_as_0_log_0_s12()
      + twoLoopFlag * Al4p * (calc_coeff_as_1_log_0_s12() + lmMt * calc_coeff_as_1_log_1_s12())
      + threeLoopFlag * pow2(Al4p) * (calc_coeff_as_2_log_0_s12() + lmMt * calc_coeff_as_2_log_1_s12()
	 + pow2(lmMt) * calc_coeff_as_2_log_2_s12());
}

/**
 * 	@return The diagonal (1, 1) matrix element of the Higgs mass matrix as a double for the hierarchy 'H3'
 */
double himalaya::H3::getS1() {
   return s1;
}

/**
 * 	@return The diagonal (2, 2) matrix element of the Higgs mass matrix as a double for the hierarchy 'H3'
 */
double himalaya::H3::getS2() {
   return s2;
}

/**
 * 	@return The off-diagonal (1, 2) = (2, 1) matrix element of the Higgs mass matrix as a double for the hierarchy 'H3'
 */
double himalaya::H3::getS12() {
   return s12;
}

/**
 * 	@return returns the susy log^0 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double himalaya::H3::calc_coef_at_as2_no_sm_logs_log0(){

   const double result =
      ((Mt*pow2(Sbeta)*(1470*Dmsqst1*pow4(Msq)*(80*Dmglst1*pow2(Mgl)*(-2*Mt*
        pow2(Dmst12)*pow2(Mst2)*(-1516*Mst1*Mt*s2t + 7614*pow2(Mt) + 3525*pow2(
        Mst1)*pow2(s2t)) + pow3(Dmst12)*(-564*Mst1*s2t*pow2(Mt) + 14100*Mt*
        pow2(Mst1)*pow2(s2t) + 15228*pow3(Mt) + 125*pow3(Mst1)*pow3(s2t)) + 4*
        Dmst12*(3807*Mt - 1375*Mst1*s2t)*pow2(Mt)*pow4(Mst2) + 29200*pow3(Mt)*
        pow6(Mst2)) + 16*pow3(Dmglst1)*(-8*Mt*pow2(Dmst12)*pow2(Mst2)*(2729*
        Mst1*Mt*s2t + 10203*pow2(Mt) + 4575*pow2(Mst1)*pow2(s2t)) + pow3(
        Dmst12)*(34172*Mst1*s2t*pow2(Mt) + 73200*Mt*pow2(Mst1)*pow2(s2t) +
        81624*pow3(Mt) + 625*pow3(Mst1)*pow3(s2t)) + 12*Dmst12*(6802*Mt + 791*
        Mst1*s2t)*pow2(Mt)*pow4(Mst2) + 153928*pow3(Mt)*pow6(Mst2)) + 16*Mgl*
        pow2(Dmglst1)*(-(Mt*pow2(Dmst12)*pow2(Mst2)*(8660*Mst1*Mt*s2t + 78882*
        pow2(Mt) + 35925*pow2(Mst1)*pow2(s2t))) + pow3(Dmst12)*(21000*Mst1*s2t*
        pow2(Mt) + 71850*Mt*pow2(Mst1)*pow2(s2t) + 78882*pow3(Mt) + 625*pow3(
        Mst1)*pow3(s2t)) + 2*Dmst12*(39441*Mt - 1840*Mst1*s2t)*pow2(Mt)*pow4(
        Mst2) + 160280*pow3(Mt)*pow6(Mst2)) + 5*pow3(Mgl)*(4*Mt*pow2(Dmst12)*
        pow2(Mst2)*(19600*Mst1*Mt*s2t + (-9122 + 14175*z3)*pow2(Mt) + 450*(2 -
        21*z3)*pow2(Mst1)*pow2(s2t)) + pow3(Dmst12)*(20800*Mst1*s2t*pow2(Mt) +
        75*Mt*(274 + 63*z3)*pow2(Mst1)*pow2(s2t) - 16*(-6536 + 14175*z3)*pow3(
        Mt) + 8400*pow3(Mst1)*pow3(s2t)) + 200*Dmst12*(-888*Mst1*s2t + Mt*(-158
        + 567*z3))*pow2(Mt)*pow4(Mst2) + 201600*pow3(Mt)*pow6(Mst2))) + 1470*
        pow2(Dmsqst1)*pow2(Msq)*(80*Dmglst1*pow2(Mgl)*(-6*Mt*pow2(Dmst12)*pow2(
        Mst2)*(-347*Mst1*Mt*s2t + 2538*pow2(Mt) + 1175*pow2(Mst1)*pow2(s2t)) +
        pow3(Dmst12)*(386*Mst1*s2t*pow2(Mt) + 14100*Mt*pow2(Mst1)*pow2(s2t) +
        15228*pow3(Mt) + 125*pow3(Mst1)*pow3(s2t)) + 2*Dmst12*(7614*Mt - 2275*
        Mst1*s2t)*pow2(Mt)*pow4(Mst2) + 22600*pow3(Mt)*pow6(Mst2)) + 16*pow3(
        Dmglst1)*(-2*Mt*pow2(Dmst12)*pow2(Mst2)*(13291*Mst1*Mt*s2t + 40812*
        pow2(Mt) + 18300*pow2(Mst1)*pow2(s2t)) + pow3(Dmst12)*(38922*Mst1*s2t*
        pow2(Mt) + 73200*Mt*pow2(Mst1)*pow2(s2t) + 81624*pow3(Mt) + 625*pow3(
        Mst1)*pow3(s2t)) + 2*Dmst12*(40812*Mt + 7121*Mst1*s2t)*pow2(Mt)*pow4(
        Mst2) + 122500*pow3(Mt)*pow6(Mst2)) + 16*Mgl*pow2(Dmglst1)*(-3*Mt*pow2(
        Dmst12)*pow2(Mst2)*(4470*Mst1*Mt*s2t + 26294*pow2(Mt) + 11975*pow2(
        Mst1)*pow2(s2t)) + pow3(Dmst12)*(25750*Mst1*s2t*pow2(Mt) + 71850*Mt*
        pow2(Mst1)*pow2(s2t) + 78882*pow3(Mt) + 625*pow3(Mst1)*pow3(s2t)) + 2*
        Dmst12*(39441*Mt + 535*Mst1*s2t)*pow2(Mt)*pow4(Mst2) + 128066*pow3(Mt)*
        pow6(Mst2)) + 5*pow3(Mgl)*(Mt*pow2(Dmst12)*pow2(Mst2)*(36800*Mst1*Mt*
        s2t + 20112*pow2(Mt) + 75*(458 - 945*z3)*pow2(Mst1)*pow2(s2t)) + pow3(
        Dmst12)*(62400*Mst1*s2t*pow2(Mt) + 1575*Mt*(-26 + 45*z3)*pow2(Mst1)*
        pow2(s2t) + (47976 - 170100*z3)*pow3(Mt) + 8400*pow3(Mst1)*pow3(s2t)) +
        100*Dmst12*(-1360*Mst1*s2t + 63*Mt*(-14 + 27*z3))*pow2(Mt)*pow4(Mst2) +
        200*(410 + 567*z3)*pow3(Mt)*pow6(Mst2))) + pow6(Msq)*(-49*pow3(Mgl)*(-
        20*Mt*pow2(Dmst12)*pow2(Mst2)*(300*Mst1*Mt*s2t*(-17512 + 14805*z3) + (-
        375892 + 621675*z3)*pow2(Mt) - 900*(-1226 + 495*z3)*pow2(Mst1)*pow2(
        s2t)) + pow3(Dmst12)*(-30*Mst1*s2t*(-1117238 + 877575*z3)*pow2(Mt) -
        2250*Mt*(-5570 + 333*z3)*pow2(Mst1)*pow2(s2t) + (-41715182 + 38174625*
        z3)*pow3(Mt) + 3000*(-2722 + 2259*z3)*pow3(Mst1)*pow3(s2t)) - 6000*
        Dmst12*(81*Mt*(-406 + 285*z3) + 8*Mst1*s2t*(-694 + 477*z3))*pow2(Mt)*
        pow4(Mst2) + 48000*(623 + 963*z3)*pow3(Mt)*pow6(Mst2)) - 196*Dmglst1*
        pow2(Mgl)*(Mt*pow2(Dmst12)*pow2(Mst2)*(-40*Mst1*Mt*s2t*(-4511549 +
        3729375*z3) + (-45149198 + 35285625*z3)*pow2(Mt) + 47250*(-430 + 207*
        z3)*pow2(Mst1)*pow2(s2t)) + pow3(Dmst12)*(2*Mst1*s2t*(28188929 -
        23099625*z3)*pow2(Mt) - 225*Mt*(-160136 + 100785*z3)*pow2(Mst1)*pow2(
        s2t) + 6*(-19196074 + 14320125*z3)*pow3(Mt) + 1125*(-22174 + 18791*z3)*
        pow3(Mst1)*pow3(s2t)) - 40*Dmst12*(150*Mst1*s2t*(-19856 + 17667*z3) +
        Mt*(-5136871 + 3912300*z3))*pow2(Mt)*pow4(Mst2) - 6000*(-31142 + 22653*
        z3)*pow3(Mt)*pow6(Mst2)) + 2*Mgl*pow2(Dmglst1)*(Mt*pow2(Dmst12)*pow2(
        Mst2)*(196*Mst1*Mt*s2t*(-353948822 + 292953375*z3) + (31554389946 -
        25980577875*z3)*pow2(Mt) + 22050*(26498 + 49425*z3)*pow2(Mst1)*pow2(
        s2t)) + 5*pow3(Dmst12)*(-196*Mst1*s2t*(-148185343 + 123161625*z3)*pow2(
        Mt) + 4410*Mt*(-612347 + 438045*z3)*pow2(Mst1)*pow2(s2t) + (-2800003036
        + 2827317150*z3)*pow3(Mt) - 735*(-2573582 + 2144805*z3)*pow3(Mst1)*
        pow3(s2t)) + 392*Dmst12*(260*Mst1*s2t*(-579323 + 490725*z3) + Mt*(-
        125277461 + 96491250*z3))*pow2(Mt)*pow4(Mst2) + 3920*(-12894992 +
        9523575*z3)*pow3(Mt)*pow6(Mst2)) + 4*pow3(Dmglst1)*(2*Mt*pow2(Dmst12)*
        pow2(Mst2)*(Mst1*Mt*s2t*(-51549748862 + 42804507375*z3) + (15344105782
        - 12740986125*z3)*pow2(Mt) + 36750*(-109771 + 107379*z3)*pow2(Mst1)*
        pow2(s2t)) + pow3(Dmst12)*(2*Mst1*s2t*(137797425107 - 114549584625*z3)*
        pow2(Mt) - 3675*Mt*(-973342 + 1115325*z3)*pow2(Mst1)*pow2(s2t) + 12*(-
        1749149438 + 1631424375*z3)*pow3(Mt) - 6370*(-2386547 + 1986525*z3)*
        pow3(Mst1)*pow3(s2t)) + 8*Dmst12*(49*Mst1*s2t*(-244084964 + 203974875*
        z3) + Mt*(-5048328734 + 3923356500*z3))*pow2(Mt)*pow4(Mst2) + 2352*(-
        21126629 + 16012125*z3)*pow3(Mt)*pow6(Mst2))) + 1470*pow3(Dmsqst1)*(80*
        Dmglst1*pow2(Mgl)*(-2*Mt*pow2(Dmst12)*pow2(Mst2)*(-566*Mst1*Mt*s2t +
        7614*pow2(Mt) + 3525*pow2(Mst1)*pow2(s2t)) + pow3(Dmst12)*(1336*Mst1*
        s2t*pow2(Mt) + 14100*Mt*pow2(Mst1)*pow2(s2t) + 15228*pow3(Mt) + 125*
        pow3(Mst1)*pow3(s2t)) + 36*Dmst12*(423*Mt - 100*Mst1*s2t)*pow2(Mt)*
        pow4(Mst2) + 16000*pow3(Mt)*pow6(Mst2)) + 16*pow3(Dmglst1)*(-12*Mt*
        pow2(Dmst12)*pow2(Mst2)*(2611*Mst1*Mt*s2t + 6802*pow2(Mt) + 3050*pow2(
        Mst1)*pow2(s2t)) + pow3(Dmst12)*(43672*Mst1*s2t*pow2(Mt) + 73200*Mt*
        pow2(Mst1)*pow2(s2t) + 81624*pow3(Mt) + 625*pow3(Mst1)*pow3(s2t)) + 8*
        Dmst12*(10203*Mt + 2374*Mst1*s2t)*pow2(Mt)*pow4(Mst2) + 91072*pow3(Mt)*
        pow6(Mst2)) + 16*Mgl*pow2(Dmglst1)*(-(Mt*pow2(Dmst12)*pow2(Mst2)*(
        18160*Mst1*Mt*s2t + 78882*pow2(Mt) + 35925*pow2(Mst1)*pow2(s2t))) +
        pow3(Dmst12)*(30500*Mst1*s2t*pow2(Mt) + 71850*Mt*pow2(Mst1)*pow2(s2t) +
        78882*pow3(Mt) + 625*pow3(Mst1)*pow3(s2t)) + 6*Dmst12*(13147*Mt + 970*
        Mst1*s2t)*pow2(Mt)*pow4(Mst2) + 95852*pow3(Mt)*pow6(Mst2)) + 5*pow3(
        Mgl)*(8400*pow3(Dmst12)*pow3(Mst1)*pow3(s2t) + 1600*Dmst12*Mst1*s2t*
        pow2(Mt)*(65*pow2(Dmst12) - 3*Dmst12*pow2(Mst2) - 59*pow4(Mst2)) - 75*
        Mt*pow2(Mst1)*(14*(-62 + 99*z3)*pow2(Dmst12)*pow2(Mst2)*pow2(s2t) + (
        1366 - 1827*z3)*pow2(s2t)*pow3(Dmst12) + 192*pow6(Mst2)) - 4*pow3(Mt)*(
        (-19178 + 14175*z3)*pow2(Dmst12)*pow2(Mst2) + 14*(154 + 2025*z3)*pow3(
        Dmst12) - 100*Dmst12*(-362 + 567*z3)*pow4(Mst2) + 100*(94 - 567*z3)*
        pow6(Mst2))))))/(2.3814e7*pow3(Mgl)*pow6(Msq)*pow6(Mst2)))/pow4(Mt)/
        pow2(Sbeta)*12.;

   return result;
}

/**
 * 	@return returns the susy log^1 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double himalaya::H3::calc_coef_at_as2_no_sm_logs_log1(){

   const double result =
      ((Mt*pow2(Sbeta)*(2*pow3(Dmglst1)*(2940*pow2(Dmsqst1)*pow2(Msq)*pow2(Mt)*(
        -((42*Mt + 11*Mst1*s2t)*pow2(Dmst12)*pow2(Mst2)) + (42*Mt + 51*Mst1*
        s2t)*pow3(Dmst12) + Dmst12*(42*Mt - 29*Mst1*s2t)*pow4(Mst2) - 65*Mt*
        pow6(Mst2)) + 5880*pow2(Mt)*pow3(Dmsqst1)*(-3*(7*Mt + 6*Mst1*s2t)*pow2(
        Dmst12)*pow2(Mst2) + (21*Mt + 38*Mst1*s2t)*pow3(Dmst12) + Dmst12*(21*Mt
        - 2*Mst1*s2t)*pow4(Mst2) - 37*Mt*pow6(Mst2)) + 5880*Dmsqst1*pow2(Mt)*
        pow4(Msq)*(7*(-3*Mt + Mst1*s2t)*pow2(Dmst12)*pow2(Mst2) + (21*Mt + 13*
        Mst1*s2t)*pow3(Dmst12) + 3*Dmst12*(7*Mt - 9*Mst1*s2t)*pow4(Mst2) - 28*
        Mt*pow6(Mst2)) + pow6(Msq)*(-10*Mt*pow2(Dmst12)*pow2(Mst2)*(623998*
        Mst1*Mt*s2t + 46947*pow2(Mt) + 134505*pow2(Mst1)*pow2(s2t)) + pow3(
        Dmst12)*(17975555*Mst1*s2t*pow2(Mt) + 1956570*Mt*pow2(Mst1)*pow2(s2t) +
        2032034*pow3(Mt) + 1187760*pow3(Mst1)*pow3(s2t)) - 2*Dmst12*(546547*Mt
        + 2011058*Mst1*s2t)*pow2(Mt)*pow4(Mst2) - 11351536*pow3(Mt)*pow6(Mst2))
        ) - 98*Dmglst1*pow2(Mgl)*(1200*pow2(Mt)*pow3(Dmsqst1)*((-3*Mt + 2*Mst1*
        s2t)*pow2(Dmst12)*pow2(Mst2) + (3*Mt - 4*Mst1*s2t)*pow3(Dmst12) + 3*
        Dmst12*Mt*pow4(Mst2) + 5*Mt*pow6(Mst2)) + 600*Dmsqst1*pow2(Mt)*pow4(
        Msq)*(-((6*Mt + Mst1*s2t)*pow2(Dmst12)*pow2(Mst2)) + (6*Mt - 3*Mst1*
        s2t)*pow3(Dmst12) + Dmst12*(6*Mt + 5*Mst1*s2t)*pow4(Mst2) + 25*Mt*pow6(
        Mst2)) + 300*pow2(Dmsqst1)*pow2(Msq)*pow2(Mt)*(3*(-4*Mt + Mst1*s2t)*
        pow2(Dmst12)*pow2(Mst2) + (12*Mt - 11*Mst1*s2t)*pow3(Dmst12) + Dmst12*(
        12*Mt + 5*Mst1*s2t)*pow4(Mst2) + 35*Mt*pow6(Mst2)) + pow6(Msq)*(Mt*
        pow2(Dmst12)*pow2(Mst2)*(29758*Mst1*Mt*s2t + 6677*pow2(Mt) + 22350*
        pow2(Mst1)*pow2(s2t)) + pow3(Dmst12)*(-34587*Mst1*s2t*pow2(Mt) - 18045*
        Mt*pow2(Mst1)*pow2(s2t) + 22414*pow3(Mt) + 3325*pow3(Mst1)*pow3(s2t)) -
        8*Dmst12*(4471*Mt + 6875*Mst1*s2t)*pow2(Mt)*pow4(Mst2) + 50800*pow3(Mt)
        *pow6(Mst2))) - Mgl*pow2(Dmglst1)*(5880*pow2(Mt)*pow3(Dmsqst1)*((-9*Mt
        + 10*Mst1*s2t)*pow2(Dmst12)*pow2(Mst2) + (9*Mt - 50*Mst1*s2t)*pow3(
        Dmst12) + 3*Dmst12*(3*Mt + 10*Mst1*s2t)*pow4(Mst2) + 79*Mt*pow6(Mst2))
        + 5880*pow2(Dmsqst1)*pow2(Msq)*pow2(Mt)*(-3*(3*Mt + 5*Mst1*s2t)*pow2(
        Dmst12)*pow2(Mst2) + (9*Mt - 25*Mst1*s2t)*pow3(Dmst12) + Dmst12*(9*Mt +
        55*Mst1*s2t)*pow4(Mst2) + 112*Mt*pow6(Mst2)) + 5880*Dmsqst1*pow2(Mt)*
        pow4(Msq)*(-((9*Mt + 40*Mst1*s2t)*pow2(Dmst12)*pow2(Mst2)) + 9*Mt*pow3(
        Dmst12) + Dmst12*(9*Mt + 80*Mst1*s2t)*pow4(Mst2) + 145*Mt*pow6(Mst2)) +
        pow6(Msq)*(3*Mt*pow2(Dmst12)*pow2(Mst2)*(2334948*Mst1*Mt*s2t - 104761*
        pow2(Mt) + 112210*pow2(Mst1)*pow2(s2t)) + pow3(Dmst12)*(-10892014*Mst1*
        s2t*pow2(Mt) + 1366365*Mt*pow2(Mst1)*pow2(s2t) + 177178*pow3(Mt) -
        363580*pow3(Mst1)*pow3(s2t)) + 196*Dmst12*(2303*Mt - 30942*Mst1*s2t)*
        pow2(Mt)*pow4(Mst2) + 13177472*pow3(Mt)*pow6(Mst2))) - 49*pow3(Mgl)*(-
        150*Dmsqst1*Mt*pow4(Msq)*(pow2(Dmst12)*pow2(Mst2)*(-80*Mst1*Mt*s2t -
        17*pow2(Mt) + 180*pow2(Mst1)*pow2(s2t)) + 2*(20*Mst1*Mt*s2t + 187*pow2(
        Mt) - 90*pow2(Mst1)*pow2(s2t))*pow3(Dmst12) - 20*Dmst12*Mt*(17*Mt - 6*
        Mst1*s2t)*pow4(Mst2) - 1080*pow2(Mt)*pow6(Mst2)) - 150*Mt*pow2(Dmsqst1)
        *pow2(Msq)*(-4*pow2(Dmst12)*pow2(Mst2)*(10*Mst1*Mt*s2t + 3*pow2(Mt) -
        45*pow2(Mst1)*pow2(s2t)) + 9*(41*pow2(Mt) - 20*pow2(Mst1)*pow2(s2t))*
        pow3(Dmst12) + 5*Dmst12*Mt*(-69*Mt + 16*Mst1*s2t)*pow4(Mst2) - 1010*
        pow2(Mt)*pow6(Mst2)) - 150*Mt*pow3(Dmsqst1)*(pow2(Dmst12)*pow2(Mst2)*(-
        7*pow2(Mt) + 180*pow2(Mst1)*pow2(s2t)) + 4*(-10*Mst1*Mt*s2t + 91*pow2(
        Mt) - 45*pow2(Mst1)*pow2(s2t))*pow3(Dmst12) + 10*Dmst12*Mt*(-35*Mt + 4*
        Mst1*s2t)*pow4(Mst2) - 940*pow2(Mt)*pow6(Mst2)) - pow6(Msq)*(-2*Mt*
        pow2(Dmst12)*pow2(Mst2)*(25850*Mst1*Mt*s2t + 21033*pow2(Mt) + 75*pow2(
        Mst1)*pow2(s2t)) + pow3(Dmst12)*(68264*Mst1*s2t*pow2(Mt) + 5700*Mt*
        pow2(Mst1)*pow2(s2t) + 36107*pow3(Mt) + 250*pow3(Mst1)*pow3(s2t)) +
        200*Dmst12*(131*Mt + 100*Mst1*s2t)*pow2(Mt)*pow4(Mst2) + 400*(-533 +
        54*z3)*pow3(Mt)*pow6(Mst2)))))/(99225.*pow3(Mgl)*pow6(Msq)*pow6(Mst2))
	)/pow4(Mt)/pow2(Sbeta)*12.; 

   return result;
}

/**
 * 	@return returns the susy log^2 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double himalaya::H3::calc_coef_at_as2_no_sm_logs_log2(){

   const double result =
      ((2*Mt*pow2(Sbeta)*(7*Dmglst1*pow2(Mgl)*pow4(Msq)*(-2*Mt*pow2(Dmst12)*
        pow2(Mst2)*(-1304*Mst1*Mt*s2t + 888*pow2(Mt) + 645*pow2(Mst1)*pow2(s2t)
        ) + pow3(Dmst12)*(-1544*Mst1*s2t*pow2(Mt) + 2250*Mt*pow2(Mst1)*pow2(
        s2t) + 1640*pow3(Mt) - 275*pow3(Mst1)*pow3(s2t)) + 8*Dmst12*(239*Mt -
        35*Mst1*s2t)*pow2(Mt)*pow4(Mst2) + 6520*pow3(Mt)*pow6(Mst2)) + pow3(
        Dmglst1)*pow4(Msq)*(-4*Mt*pow2(Dmst12)*pow2(Mst2)*(-6158*Mst1*Mt*s2t +
        7818*pow2(Mt) - 5565*pow2(Mst1)*pow2(s2t)) + pow3(Dmst12)*(-109632*
        Mst1*s2t*pow2(Mt) - 8820*Mt*pow2(Mst1)*pow2(s2t) + 98096*pow3(Mt) -
        9765*pow3(Mst1)*pow3(s2t)) + 16*Dmst12*(-2222*Mt + 5257*Mst1*s2t)*pow2(
        Mt)*pow4(Mst2) + 36232*pow3(Mt)*pow6(Mst2)) + Mgl*pow2(Dmglst1)*pow4(
        Msq)*(-3*Mt*pow2(Dmst12)*pow2(Mst2)*(-7448*Mst1*Mt*s2t + 11772*pow2(Mt)
        + 35*pow2(Mst1)*pow2(s2t)) + pow3(Dmst12)*(-53536*Mst1*s2t*pow2(Mt) +
        16905*Mt*pow2(Mst1)*pow2(s2t) + 68252*pow3(Mt) - 5285*pow3(Mst1)*pow3(
        s2t)) + 28*Dmst12*(85*Mt + 1164*Mst1*s2t)*pow2(Mt)*pow4(Mst2) + 49588*
        pow3(Mt)*pow6(Mst2)) + 7*pow3(Mgl)*(Mt*pow2(Dmst12)*pow2(Mst2)*(1760*
        Mst1*Mt*s2t + 538*pow2(Mt) + 105*pow2(Mst1)*pow2(s2t))*pow4(Msq) -
        pow3(Dmst12)*(1032*Mst1*s2t*pow2(Mt) + 480*Mt*pow2(Mst1)*pow2(s2t) +
        644*pow3(Mt) - 45*pow3(Mst1)*pow3(s2t))*pow4(Msq) - 40*Dmst12*(11*Mt +
        61*Mst1*s2t)*pow2(Mt)*pow4(Msq)*pow4(Mst2) + 10*pow3(Mt)*(45*pow2(
        Dmsqst1) + 90*Dmsqst1*pow2(Msq) + 442*pow4(Msq))*pow6(Mst2))))/(945.*
        pow3(Mgl)*pow4(Msq)*pow6(Mst2)))/pow4(Mt)/pow2(Sbeta)*12.; 

   return result;
}

/**
 * 	@return returns the susy log^3 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double himalaya::H3::calc_coef_at_as2_no_sm_logs_log3(){

   const double result =
      ((-224*pow2(Sbeta)*pow4(Mt))/9.)/pow4(Mt)/pow2(Sbeta)*12.; 

   return result;
}

/// calc coefficient

/// calc coefficient O(as^0,log(mu^2/mt^2)^0) s1
double himalaya::H3::calc_coeff_as_0_log_0_s1 () {
   return -(pow2(Dmst12)*(-(Dmst12*xDmst12)+pow2(Mst2))*pow2(Mt)*pow2(MuSUSY)*
   pow2(s2t))/(48.*pow6(Mst2));
}

/// calc coefficient O(as^1,log(mu^2/mt^2)^0) s1
double himalaya::H3::calc_coeff_as_1_log_0_s1 () {
   return -(Dmst12*s2t*pow2(Mt)*pow2(MuSUSY)*(4*xDmglst1*pow3(Dmglst1)*((3*(9+10
   *lmMst1)*Mt-10*(13-12*lmMst1)*Mst1*s2t)*xDmst12*pow2(Dmst12)-Dmst12*(52*Mt-5*
   (13-12*lmMst1)*Mst1*s2t)*pow2(Mst2)+(77-30*lmMst1)*Mt*pow4(Mst2))+20*Dmglst1*
   pow2(Mgl)*(-(((5-6*lmMst1)*Mt+4*(1-3*lmMst1)*Mst1*s2t)*xDmst12*pow2(Dmst12))+
   2*Dmst12*(1-3*lmMst1)*Mst1*s2t*pow2(Mst2)+(5-6*lmMst1)*Mt*pow4(Mst2))-5*pow3(
   Mgl)*(-((4*(5+6*lmMst1)*Mt-(1+6*lmMst1)*Mst1*s2t)*xDmst12*pow2(Dmst12))+
   Dmst12*(4*Mt-3*Mst1*s2t)*pow2(Mst2)+12*(1+2*lmMst1)*Mt*pow4(Mst2))+10*Mgl*
   pow2(Dmglst1)*(-3*Dmst12*(4*Mt-(5-6*lmMst1)*Mst1*s2t)*pow2(Mst2)+2*((Mt+6*
   lmMst1*Mt-3*(5-6*lmMst1)*Mst1*s2t)*xDmst12*pow2(Dmst12)+(11-6*lmMst1)*Mt*
   pow4(Mst2)))))/(270.*Mst1*pow3(Mgl)*pow6(Mst2));
}

/// calc coefficient O(as^2,log(mu^2/mt^2)^0) s1
double himalaya::H3::calc_coeff_as_2_log_0_s1() {
   return -(pow2(Mt)*pow2(MuSUSY)*(-4*xDmst12*pow3(Dmst12)*(((72*pow2(Dmglst1)*(3891491+27200*lmMst1-
19200*pow2(lmMst1))+
200*Dmglst1*Mgl*(403559+384*lmMst1-4608*pow2(lmMst1))-
15*(1763661-47104*lmMst1-24576*pow2(lmMst1))*pow2(Mgl))*pow2(Mt))/(pow2(Mgl)*pow2(Mst1))+
((15*pow2(s2t)*(2400*pow2(Mgl)*((7+6*shiftst1+24*lmMst1*(1-shiftst2)+
30*shiftst2)*pow2(Dmsqst1)+
Dmsqst1*(7+24*lmMst1*(1-shiftst2)+36*shiftst2)*pow2(Msq))-
(8*pow2(Dmglst1)*(1732531+16896*lmMst1-24840*pow2(lmMst1))+
40*Dmglst1*Mgl*(84209+1264*lmMst1-240*pow2(lmMst1))-
(350605+4320*shiftst1+2880*shiftst2+
96*lmMst1*(115-90*shiftst1-60*shiftst2-
54*shiftst3)+8352*shiftst3-
2160*pow2(lmMst1))*pow2(Mgl))*pow4(Msq)))/pow2(Mgl)+(8*xDmglst1*pow3(Dmglst1)*(30000*Mst1*Mt*s2t*(pow2(Dmsqst1)+Dmsqst1*pow2(Msq))+
(2*Mst1*Mt*s2t*(37824007+770520*lmMst1-131400*pow2(lmMst1))
+(59957863+480000*lmMst1-230400*pow2(lmMst1))*pow2(Mt)-
15*(3044017+27472*lmMst1-48480*pow2(lmMst1))*pow2(Mst1)*pow2(s2t))*pow4(Msq)))/(pow2(Mst1)*pow3(Mgl)))/pow4(Msq)+
(s2t*(-240*Mt*(184315-400*lmMst1+
(Dmglst1*(1282471-7264*lmMst1-18120*pow2(lmMst1)))/Mgl-2760*pow2(lmMst1)+
(10*pow2(Dmglst1)*(32829-1852*lmMst1-660*pow2(lmMst1))-
(200*Dmsqst1*(5*(Dmglst1*Mgl+pow2(Dmglst1))+
21*pow2(Mgl))*(Dmsqst1+pow2(Msq)))/pow4(Msq))/pow2(Mgl))+
(12000*xDmsqst1*pow3(Dmsqst1)*(20*Mt*(Mgl*pow2(Dmglst1)+Dmglst1*pow2(Mgl)+
xDmglst1*pow3(Dmglst1))+
3*(28*Mt+
Mst1*s2t*(7+12*shiftst1+
24*(lmMst1*(1-shiftst2)+shiftst2)))*pow3(Mgl)))/(pow3(Mgl)*pow6(Msq))))/Mst1)+
pow2(Mst2)*(((480*Dmst12*Mt*s2t*(pow2(Dmglst1)*(2000*(pow2(Dmsqst1)+Dmsqst1*pow2(Msq))*(Dmst12-pow2(Mst2))+
(Dmst12*(958501+24456*lmMst1-11520*pow2(lmMst1))
-2*(1286791+5936*lmMst1-18120*pow2(lmMst1))*pow2(Mst2))*pow4(Msq))
+Dmglst1*Mgl*(2000*(pow2(Dmsqst1)+Dmsqst1*pow2(Msq))*(Dmst12-pow2(Mst2))-
(Dmst12*(949861-1944*lmMst1-11520*pow2(lmMst1))+
20*(33261-532*lmMst1-660*pow2(lmMst1))*pow2(Mst2))*pow4(Msq))+
5*pow2(Mgl)*(1680*(pow2(Dmsqst1)+Dmsqst1*pow2(Msq))*(Dmst12-pow2(Mst2))-
(3*Dmst12*(10473-40*lmMst1-256*pow2(lmMst1))+
8*(1361+10*lmMst1+54*pow2(lmMst1))*pow2(Mst2))*pow4(Msq))))/(Mst1*pow4(Msq))+
(2*pow2(Mt)*(16*pow2(Dmglst1)*((3891491+27200*lmMst1-19200*pow2(lmMst1))*(9*pow2(Dmst12)-9*Dmst12*pow2(Mst2))-
50*(345581+4896*lmMst1-3456*pow2(lmMst1))*pow4(Mst2))+
400*Dmglst1*Mgl*((403559+384*lmMst1-4608*pow2(lmMst1))*(pow2(Dmst12)-Dmst12*pow2(Mst2))-
24*(9631+16*lmMst1-192*pow2(lmMst1))*pow4(Mst2))-15*pow2(Mgl)*(pow2(Dmst12)*(852541+9216*lmMst1+6144*pow2(lmMst1))+
160*Dmst12*(11389-704*lmMst1-384*pow2(lmMst1))*pow2(Mst2)+
1920*(349-56*lmMst1-32*pow2(lmMst1))*pow4(Mst2))))/pow2(Mst1))/pow2(Mgl)+((-16*xDmglst1*pow3(Dmglst1)*(2*Dmst12*Mt*pow2(Mst2)*(30000*Mst1*s2t*(pow2(Dmsqst1)+Dmsqst1*pow2(Msq))+
(4*Mst1*s2t*(31025111+290880*lmMst1-
251100*pow2(lmMst1))+
Mt*(59957863+480000*lmMst1-
230400*pow2(lmMst1)))*pow4(Msq))-
pow2(Dmst12)*(60000*Mst1*Mt*s2t*(pow2(Dmsqst1)+Dmsqst1*pow2(Msq))+
(2*(Mst1*Mt*s2t*(99874229+1352280*lmMst1-
633600*pow2(lmMst1))+
(59957863+480000*lmMst1-
230400*pow2(lmMst1))*pow2(Mt))-
15*(3044017+27472*lmMst1-
48480*pow2(lmMst1))*pow2(Mst1)*pow2(s2t))*pow4(Msq))+
24*(3877891+46400*lmMst1-19200*pow2(lmMst1))*pow2(Mt)*pow4(Msq)*pow4(Mst2)))/(pow2(Mst1)*pow4(Msq))
+(24000*Dmst12*s2t*xDmsqst1*pow3(Dmsqst1)*(40*Mt*(Dmst12-pow2(Mst2))*(Mgl*pow2(Dmglst1)+Dmglst1*pow2(Mgl)+
xDmglst1*pow3(Dmglst1))-
3*(-(Dmst12*(56*Mt+
Mst1*s2t*(7+42*shiftst1+
12*lmMst1*(2-shiftst1-shiftst2)-
6*shiftst2)))+
4*(14*Mt+
3*(5-2*lmMst1)*Mst1*s2t*(shiftst1-shiftst2))*pow2(Mst2))*pow3(Mgl)))/(Mst1*pow6(Msq)))/pow3(Mgl))-
240*Dmst12*pow2(s2t)*(pow2(Mst2)*(-180*shiftst3*(Dmst12*(7-6*lmMst1)-
2*(1-2*lmMst1)*pow2(Mst2))-
(3600*z2*(-(Dmst12*(shiftst1+shiftst2))+
2*(shiftst1-shiftst2)*pow2(Mst2))*(pow2(Dmsqst1)+Dmsqst1*pow2(Msq)+pow4(Msq)))/pow4(Msq))
+(pow2(Mst2)*(-((Dmst12*(300*pow2(Mgl)*((7+12*lmMst1*(2-shiftst1)+
30*shiftst1)*pow2(Dmsqst1)+
Dmsqst1*(7+12*lmMst1*(2-shiftst1)+18*shiftst1)*pow2(Msq))-
(pow2(Dmglst1)*(1732531+16896*lmMst1-24840*pow2(lmMst1))
+5*Dmglst1*Mgl*(84209+1264*lmMst1-240*pow2(lmMst1))+
10*(1429-454*lmMst1+
(-180+360*lmMst1)*shiftst1+
24*pow2(lmMst1))*pow2(Mgl))*pow4(Msq)))/pow2(Mgl))-
1800*shiftst2*(pow2(Dmsqst1)*(Dmst12*(1-2*lmMst1)+4*(2-lmMst1)*pow2(Mst2))
+(Dmst12+2*pow2(Mst2))*(Dmsqst1*(3-2*lmMst1)*pow2(Msq)+
(1-2*lmMst1)*pow4(Msq))))+
3600*shiftst1*(2*(2-lmMst1)*pow2(Dmsqst1)+
Dmsqst1*(3-2*lmMst1)*pow2(Msq)+
(1-2*lmMst1)*pow4(Msq))*pow4(Mst2))/pow4(Msq)-
(72*z2*(10*pow4(Mst2)*(10*(shiftst1-shiftst2)*xDmsqst1*pow3(Dmsqst1)+
shiftst3*pow6(Msq))-
5*Dmst12*pow2(Mst2)*(10*(shiftst1+shiftst2)*xDmsqst1*pow3(Dmsqst1)+
3*shiftst3*pow6(Msq))+
2*xDmst12*pow2(Dmst12)*(100*shiftst2*(pow2(Dmsqst1)*pow2(Msq)+xDmsqst1*pow3(Dmsqst1)+
Dmsqst1*pow4(Msq))+
(15*shiftst1+10*shiftst2+9*shiftst3)*pow6(Msq))))/pow6(Msq))+(225*z3*(pow3(Mgl)*(2*xDmst12*pow3(Dmst12)*(10080*pow2(Mst1)*pow2(s2t)*(pow2(Dmsqst1)*pow2(Msq)+xDmsqst1*pow3(Dmsqst1)+
Dmsqst1*pow4(Msq))-
(324752*Mst1*Mt*s2t+197889*pow2(Mt)-
37669*pow2(Mst1)*pow2(s2t))*pow6(Msq))-
pow2(Dmst12)*pow2(Mst2)*(10080*pow2(Mst1)*pow2(s2t)*(pow2(Dmsqst1)*pow2(Msq)+xDmsqst1*pow3(Dmsqst1)+
Dmsqst1*pow4(Msq))-
(276560*Mst1*Mt*s2t+95889*pow2(Mt)+
20592*pow2(Mst1)*pow2(s2t))*pow6(Msq))+
pow6(Msq)*(96*Dmst12*Mt*(2125*Mt+1004*Mst1*s2t)*pow4(Mst2)+72192*pow2(Mt)*pow6(Mst2)))+
8*pow6(Msq)*(Dmglst1*pow2(Mgl)*(-(pow2(Dmst12)*pow2(Mst2)*(-210810*Mst1*Mt*s2t+147834*pow2(Mt)-
50485*pow2(Mst1)*pow2(s2t)))+
2*xDmst12*(-285974*Mst1*Mt*s2t+73917*pow2(Mt)-
50485*pow2(Mst1)*pow2(s2t))*pow3(Dmst12)+
86*Dmst12*Mt*(1719*Mt+1748*Mst1*s2t)*pow4(Mst2)+
83952*pow2(Mt)*pow6(Mst2))+
Mgl*pow2(Dmglst1)*(-15*pow2(Dmst12)*pow2(Mst2)*(14054*Mst1*Mt*s2t+34398*pow2(Mt)-
13341*pow2(Mst1)*pow2(s2t))+
2*xDmst12*(-75164*Mst1*Mt*s2t+257985*pow2(Mt)-
200115*pow2(Mst1)*pow2(s2t))*pow3(Dmst12)+
26*Dmst12*Mt*(19845*Mt+21998*Mst1*s2t)*pow4(Mst2)+
253692*pow2(Mt)*pow6(Mst2))+
xDmglst1*pow3(Dmglst1)*(-(pow2(Dmst12)*pow2(Mst2)*(1475294*Mst1*Mt*s2t+884106*pow2(Mt)-
349745*pow2(Mst1)*pow2(s2t)))+
2*xDmst12*(557078*Mst1*Mt*s2t+442053*pow2(Mt)-
349745*pow2(Mst1)*pow2(s2t))*pow3(Dmst12)+
18*Dmst12*Mt*(49117*Mt+102024*Mst1*s2t)*pow4(Mst2)+
687960*pow2(Mt)*pow6(Mst2)))))/(pow2(Mst1)*pow3(Mgl)*pow6(Msq))))/(777600.*pow6(Mst2));
}

/// calc coefficient O(as^2,log(mu^2/mt^2)^1) s1
double himalaya::H3::calc_coeff_as_2_log_1_s1(){
   return (16*pow2(MuSUSY)*pow4(Mt)*(-5*pow3(Mgl)*(2*(pow2(Dmst12)*pow2(Mst2)+xDmst12*pow3(Dmst12))-
6*Dmst12*pow4(Mst2)-9*pow6(Mst2))-
10*Dmglst1*pow2(Mgl)*(pow2(Dmst12)*pow2(Mst2)-
xDmst12*pow3(Dmst12)-Dmst12*pow4(Mst2)-3*pow6(Mst2))+
4*xDmglst1*pow3(Dmglst1)*(7*pow2(Dmst12)*pow2(Mst2)-
7*(xDmst12*pow3(Dmst12)+Dmst12*pow4(Mst2))-3*pow6(Mst2))+
Mgl*pow2(Dmglst1)*(9*pow2(Dmst12)*pow2(Mst2)-
9*(xDmst12*pow3(Dmst12)+Dmst12*pow4(Mst2))+5*pow6(Mst2))))/(405.*pow2(Mst1)*pow3(Mgl)*pow6(Mst2));
}

/// calc coefficient O(as^0,log(mu^2/mt^2)^0) s2
double himalaya::H3::calc_coeff_as_0_log_0_s2(){
   return -(Mt*(Mt*pow2(Dmst12)*pow2(Mst2)*((12*Mt*MuSUSY*s2t-12*Tbeta*pow2(Mt))*pow2(Sbeta)+
Tbeta*pow2(s2t)*(pow2(MuSUSY)+
(-12*pow2(Mst1)-pow2(MuSUSY))*pow2(Sbeta)))-
xDmst12*pow3(Dmst12)*(Tbeta*(-(Mt*pow2(s2t)*(-(pow2(MuSUSY)*(1-pow2(Sbeta)))+
18*pow2(Mst1)*pow2(Sbeta)))-8*pow2(Sbeta)*pow3(Mt))+
MuSUSY*pow2(Sbeta)*(8*s2t*pow2(Mt)-pow2(Mst1)*pow3(s2t)))+
pow2(Sbeta)*(-24*Dmst12*(MuSUSY*s2t-Mt*Tbeta)*pow2(Mt)*pow4(Mst2)+
48*lmMst1*Tbeta*pow3(Mt)*pow6(Mst2))))/(48.*Tbeta*pow2(Sbeta)*pow6(Mst2));
}

/// calc coefficient O(as^0,log(mu^2/mt^2)^1) s2
double himalaya::H3::calc_coeff_as_0_log_1_s2(){
   return pow4(Mt);
}

/// calc coefficient O(as^1,log(mu^2/mt^2)^0) s2
double himalaya::H3::calc_coeff_as_1_log_0_s2(){
   return -(Mt*(pow3(Mt)*(3600*(1+lmMst1+pow2(lmMst1))-
(2*(25*Dmst12*pow2(Mgl)*(Dmst12*(5+42*lmMst1)-8*(4+3*lmMst1)*pow2(Mst2))-
2*pow2(Dmglst1)*((83-120*lmMst1)*pow2(Dmst12)+
Dmst12*(1381-2490*lmMst1)*pow2(Mst2)+
200*(11-21*lmMst1)*pow4(Mst2))+
4*Dmglst1*Mgl*((11+60*lmMst1)*pow2(Dmst12)+
25*Dmst12*(1+30*lmMst1)*pow2(Mst2)+
100*(1+12*lmMst1)*pow4(Mst2))))/(pow2(Mgl)*pow4(Mst2)))
+((Dmst12*s2t*(-4*pow2(Mt)*(25*pow2(Mgl)*(Dmst12*(4*(1+3*lmMst1)*pow2(Mst1)-
pow2(MuSUSY))-
3*(1+2*lmMst1)*pow2(Mst2)*(6*pow2(Mst1)+pow2(MuSUSY)))-
25*Dmglst1*Mgl*(8*Dmst12*(2-3*lmMst1)*pow2(Mst1)-
pow2(Mst2)*((34-60*lmMst1)*pow2(Mst1)+
(5-6*lmMst1)*pow2(MuSUSY)))-
pow2(Dmglst1)*(6*Dmst12*((169-135*lmMst1)*pow2(Mst1)+
25*pow2(MuSUSY))-
25*pow2(Mst2)*(24*(4-3*lmMst1)*pow2(Mst1)+
(11-6*lmMst1)*pow2(MuSUSY))))-
(25*Mt*(-4*Dmglst1*Mgl*(2*Dmst12*(1-3*lmMst1)*Mst1*s2t+
(5-6*lmMst1)*Mt*pow2(Mst2))-
pow2(Dmglst1)*(6*Dmst12*(5-6*lmMst1)*Mst1*s2t+
4*(11-6*lmMst1)*Mt*pow2(Mst2))+
3*pow2(Mgl)*(-(Dmst12*Mst1*s2t)+
4*(Mt+2*lmMst1*Mt)*pow2(Mst2)))*pow2(MuSUSY))/pow2(Sbeta))+
Mt*((-25*s2t*pow2(Dmst12)*(8*Dmglst1*Mgl*Mst1*s2t*(3*(1-6*lmMst1)*pow2(Mst1)+
(1-3*lmMst1)*pow2(MuSUSY))*pow2(Sbeta)+
pow2(Mgl)*(4*Mt*pow2(MuSUSY)+
3*Mst1*s2t*(12*(1+2*lmMst1)*pow2(Mst1)+pow2(MuSUSY))*pow2(Sbeta))+
6*pow2(Dmglst1)*(4*Mt*pow2(MuSUSY)+
Mst1*s2t*((26-36*lmMst1)*pow2(Mst1)+
(5-6*lmMst1)*pow2(MuSUSY))*pow2(Sbeta))))/pow2(Sbeta)+
(2*MuSUSY*(-(Dmglst1*Mgl*(300*Dmst12*Mt*((3-6*lmMst1)*Mt+
2*(1-6*lmMst1)*Mst1*s2t)*pow2(Mst2)-
3*pow2(Dmst12)*(100*Mst1*Mt*s2t+
2*(41-90*lmMst1)*pow2(Mt)+
25*(5-6*lmMst1)*pow2(Mst1)*pow2(s2t))+
200*(7-15*lmMst1)*pow2(Mt)*pow4(Mst2)))-
3*pow2(Dmglst1)*(2*Dmst12*Mt*((477-330*lmMst1)*Mt+
50*(13-18*lmMst1)*Mst1*s2t)*pow2(Mst2)-
pow2(Dmst12)*(70*Mst1*Mt*s2t+
16*(46-15*lmMst1)*pow2(Mt)+
25*(11-6*lmMst1)*pow2(Mst1)*pow2(s2t))+
400*(4-3*lmMst1)*pow2(Mt)*pow4(Mst2))+
25*pow2(Mgl)*(4*Dmst12*Mt*(2*(5+6*lmMst1)*Mt-
9*(1+2*lmMst1)*Mst1*s2t)*pow2(Mst2)-
pow2(Dmst12)*(24*(1-3*lmMst1)*Mst1*Mt*s2t+
2*(5+6*lmMst1)*pow2(Mt)+
9*(1+2*lmMst1)*pow2(Mst1)*pow2(s2t))+
72*(1+lmMst1)*pow2(Mt)*pow4(Mst2))))/Tbeta))/(pow2(Mgl)*pow4(Mst2))+
((-2*Mt*xDmglst1*pow3(Dmglst1)*(2*Dmst12*Mt*pow2(Mst2)*(Mt*((2231-990*lmMst1)*MuSUSY-
12*(301-290*lmMst1)*Mst1*Tbeta)*pow2(Sbeta)+
s2t*(-5*(77-30*lmMst1)*Tbeta*pow2(MuSUSY)*(1-pow2(Sbeta))+
(60*(71-60*lmMst1)*Mst1*MuSUSY+
8*(481-240*lmMst1)*Tbeta*pow2(Mst1))*pow2(Sbeta)))-
pow2(Dmst12)*(8*((476-90*lmMst1)*MuSUSY-
3*(51+10*lmMst1)*Mst1*Tbeta)*pow2(Mt)*pow2(Sbeta)+
4*Mt*s2t*(-130*Tbeta*pow2(MuSUSY)*(1-pow2(Sbeta))+
(360*Mst1*MuSUSY+
3*(277-155*lmMst1)*Tbeta*pow2(Mst1))*pow2(Sbeta))+
5*Mst1*pow2(s2t)*(10*(13-12*lmMst1)*Tbeta*pow2(MuSUSY)*(1-pow2(Sbeta))+
(3*(77-30*lmMst1)*Mst1*MuSUSY-
12*(71-60*lmMst1)*Tbeta*pow2(Mst1))*pow2(Sbeta)))+
8*((977-480*lmMst1)*MuSUSY-
3*(509-510*lmMst1)*Mst1*Tbeta)*pow2(Mt)*pow2(Sbeta)*pow4(Mst2)))/pow4(Mst2)-
(xDmst12*pow3(Dmst12)*(-2*Dmglst1*pow2(Mgl)*(-25*Mst1*Mt*pow2(s2t)*(8*(1-3*lmMst1)*Tbeta*pow2(MuSUSY)*(1-pow2(Sbeta))+
(3*(5-6*lmMst1)*Mst1*MuSUSY-
18*(1-4*lmMst1)*Tbeta*pow2(Mst1))*pow2(Sbeta))
-2*s2t*pow2(Mt)*(Tbeta*(125*pow2(MuSUSY)*(1-pow2(Sbeta))+
262*pow2(Mst1)*pow2(Sbeta))+
30*lmMst1*(-5*Tbeta*pow2(MuSUSY)*(1-pow2(Sbeta))+
(60*Mst1*MuSUSY-11*Tbeta*pow2(Mst1))*pow2(Sbeta)))+
pow2(Sbeta)*(4*(6*(17-30*lmMst1)*MuSUSY+
(47+870*lmMst1)*Mst1*Tbeta)*pow3(Mt)+
25*(4*(1-3*lmMst1)*MuSUSY+
(5-6*lmMst1)*Mst1*Tbeta)*pow3(Mst1)*pow3(s2t)))+pow3(Mgl)*(-100*s2t*pow2(Mt)*((5+6*lmMst1)*Tbeta*pow2(MuSUSY)*(1-pow2(Sbeta))+
(24*(1-lmMst1)*Mst1*MuSUSY+
2*(1+3*lmMst1)*Tbeta*pow2(Mst1))*pow2(Sbeta))
-25*Mst1*Mt*pow2(s2t)*(-((1+6*lmMst1)*Tbeta*pow2(MuSUSY)*(1-pow2(Sbeta)))+
12*((1+3*lmMst1)*Mst1*MuSUSY+
(1+12*lmMst1)*Tbeta*pow2(Mst1))*pow2(Sbeta))+
pow2(Sbeta)*(4*(50*(5+6*lmMst1)*MuSUSY+
(137-330*lmMst1)*Mst1*Tbeta)*pow3(Mt)-
75*(MuSUSY-2*(1+2*lmMst1)*Mst1*Tbeta)*pow3(Mst1)*pow3(s2t)))+
2*(xDmglst1*pow3(Dmglst1)*(5*Mst1*Mt*pow2(s2t)*(20*(13-12*lmMst1)*Tbeta*pow2(MuSUSY)*(1-pow2(Sbeta))+
(9*(43-10*lmMst1)*Mst1*MuSUSY-
12*(83-60*lmMst1)*Tbeta*pow2(Mst1))*pow2(Sbeta))-
2*s2t*pow2(Mt)*(15*(9+10*lmMst1)*Tbeta*pow2(MuSUSY)*(1-pow2(Sbeta))+
(60*(47-60*lmMst1)*Mst1*MuSUSY+
2*(106+285*lmMst1)*Tbeta*pow2(Mst1))*pow2(Sbeta))+
pow2(Sbeta)*(2*((1577+270*lmMst1)*MuSUSY+
12*(199-310*lmMst1)*Mst1*Tbeta)*pow3(Mt)-
5*(10*(13-12*lmMst1)*MuSUSY+
(77-30*lmMst1)*Mst1*Tbeta)*pow3(Mst1)*pow3(s2t)))+
Mgl*pow2(Dmglst1)*(15*Mst1*Mt*pow2(s2t)*(10*(5-6*lmMst1)*Tbeta*pow2(MuSUSY)*(1-pow2(Sbeta))+
(5*(17-6*lmMst1)*Mst1*MuSUSY-
(137-180*lmMst1)*Tbeta*pow2(Mst1))*pow2(Sbeta))-10*s2t*pow2(Mt)*(5*(1+6*lmMst1)*Tbeta*pow2(MuSUSY)*(1-pow2(Sbeta))+
(12*(29-45*lmMst1)*Mst1*MuSUSY+
6*(2+15*lmMst1)*Tbeta*pow2(Mst1))*pow2(Sbeta))
+pow2(Sbeta)*(2*(3*(259+90*lmMst1)*MuSUSY+
91*(17-30*lmMst1)*Mst1*Tbeta)*pow3(Mt)-
25*(3*(5-6*lmMst1)*MuSUSY+
(11-6*lmMst1)*Mst1*Tbeta)*pow3(Mst1)*pow3(s2t))))))/pow6(Mst2))/(Tbeta*pow2(Sbeta)*pow3(Mgl)))/Mst1))/1350.;
}

/// calc coefficient O(as^1,log(mu^2/mt^2)^1) s2
double himalaya::H3::calc_coeff_as_1_log_1_s2(){
   return (-2*pow3(Mt)*(pow3(Mgl)*(-5*pow2(Dmst12)*(Mt*(2*MuSUSY+5*Mst1*Tbeta)+8*s2t*Tbeta*pow2(Mst1))*pow2(Mst2)+6*Mst1*(3*Mt+5*Mst1*s2t)*Tbeta*xDmst12*pow3(Dmst12)+20*Dmst12*(Mt*(MuSUSY+2*Mst1*Tbeta)+3*s2t*Tbeta*pow2(Mst1))*pow4(Mst2)
+60*Mt*(MuSUSY-(1-2*lmMst1)*Mst1*Tbeta)*pow6(Mst2))-
2*(Dmglst1*pow2(Mgl)*(pow2(Dmst12)*(Mt*(MuSUSY-4*Mst1*Tbeta)+10*s2t*Tbeta*pow2(Mst1))*pow2(Mst2)-xDmst12*(2*Mt*(MuSUSY+Mst1*Tbeta)+9*s2t*Tbeta*pow2(Mst1))*pow3(Dmst12)+10*Dmst12*Mst1*(Mt-Mst1*s2t)*Tbeta*pow4(Mst2)-10*Mt*(MuSUSY-4*Mst1*Tbeta)*pow6(Mst2))+
xDmglst1*pow3(Dmglst1)*(-(pow2(Dmst12)*(Mt*(2*MuSUSY+4*Mst1*Tbeta)+
s2t*Tbeta*pow2(Mst1))*pow2(Mst2))+
xDmst12*(Mt*(MuSUSY+4*Mst1*Tbeta)-s2t*Tbeta*pow2(Mst1))*pow3(Dmst12)+Dmst12*(Mt*(3*MuSUSY+4*Mst1*Tbeta)+4*s2t*Tbeta*pow2(Mst1))*pow4(Mst2)+4*Mt*(MuSUSY+9*Mst1*Tbeta)*pow6(Mst2))+
Mgl*pow2(Dmglst1)*(pow2(Dmst12)*(-2*Mt*(MuSUSY+2*Mst1*Tbeta)+3*s2t*Tbeta*pow2(Mst1))*pow2(Mst2)-xDmst12*(-(Mt*(MuSUSY+Mst1*Tbeta))+5*s2t*Tbeta*pow2(Mst1))*pow3(Dmst12)+Mt*(Dmst12*(3*MuSUSY+7*Mst1*Tbeta)*pow4(Mst2)+
40*Mst1*Tbeta*pow6(Mst2))))))/(45.*Mst1*Tbeta*pow3(Mgl)*pow6(Mst2));
}

/// calc coefficient O(as^1,log(mu^2/mt^2)^2) s2
double himalaya::H3::calc_coeff_as_1_log_2_s2(){
   return 8*pow4(Mt);
}

/// calc coefficient O(as^2,log(mu^2/mt^2)^0) s2
double himalaya::H3::calc_coeff_as_2_log_0_s2(){
   return -(((pow2(Dmst12)*pow2(Mt)*pow2(s2t)*(pow2(Dmglst1)*(114960*(pow2(Dmsqst1)+Dmsqst1*pow2(Msq))*pow2(Mst1)-(12*(13249-916*lmMst1-60*pow2(lmMst1))*pow2(Mst1)-
(1732531+16896*lmMst1-24840*pow2(lmMst1))*pow2(MuSUSY))*pow4(Msq))+
5*Dmglst1*Mgl*(22560*(pow2(Dmsqst1)+Dmsqst1*pow2(Msq))*pow2(Mst1)-(24*(4515-596*lmMst1-516*pow2(lmMst1))*pow2(Mst1)-
(84209+1264*lmMst1-240*pow2(lmMst1))*pow2(MuSUSY))*pow4(Msq))-10*pow2(Mgl)*(30*Dmsqst1*pow2(Msq)*(12*(1+12*lmMst1)*pow2(Mst1)+
(7+24*lmMst1)*pow2(MuSUSY))+
15*pow2(Dmsqst1)*((229+288*lmMst1)*pow2(Mst1)+
2*(7+24*lmMst1)*pow2(MuSUSY))+
(24*(613-lmMst1+21*pow2(lmMst1))*pow2(Mst1)-
(1429-454*lmMst1+24*pow2(lmMst1))*pow2(MuSUSY))*pow4(Msq))))/3240.-
(Dmst12*s2t*pow3(Mt)*(pow2(Dmglst1)*(1200*pow2(Dmsqst1)*(pow2(Mst2)*((107-330*lmMst1)*pow2(Mst1)-
125*pow2(MuSUSY))-
Dmst12*(9*(149-10*lmMst1)*pow2(Mst1)-
125*pow2(MuSUSY)))-
1200*Dmsqst1*pow2(Msq)*(Dmst12*((866-240*lmMst1)*pow2(Mst1)-
125*pow2(MuSUSY))+
pow2(Mst2)*(16*(23+30*lmMst1)*pow2(Mst1)+
125*pow2(MuSUSY)))-
(10*pow2(Mst2)*(8*(7531199-92826*lmMst1-104760*pow2(lmMst1))*pow2(Mst1)+
15*(1286791+5936*lmMst1-18120*pow2(lmMst1))*pow2(MuSUSY))+
Dmst12*((707897644+8577360*lmMst1-
5745600*pow2(lmMst1))*pow2(Mst1)-
75*(958501+24456*lmMst1-11520*pow2(lmMst1))*pow2(MuSUSY)))*pow4(Msq))+
5*Dmglst1*Mgl*(240*pow2(Dmsqst1)*(-25*pow2(Mst2)*((91+6*lmMst1)*pow2(Mst1)+5*pow2(MuSUSY))+
Dmst12*(3*(347-30*lmMst1)*pow2(Mst1)+
125*pow2(MuSUSY)))-
240*Dmsqst1*pow2(Msq)*(25*pow2(Mst2)*(2*(55+6*lmMst1)*pow2(Mst1)+5*pow2(MuSUSY))-
Dmst12*(4*(379+15*lmMst1)*pow2(Mst1)+
125*pow2(MuSUSY)))-
(Dmst12*((36092392+714192*lmMst1-
938880*pow2(lmMst1))*pow2(Mst1)+
15*(949861-1944*lmMst1-11520*pow2(lmMst1))*pow2(MuSUSY))+
300*pow2(Mst2)*(16*(4964-275*lmMst1+21*pow2(lmMst1))*pow2(Mst1)+
(33261-532*lmMst1-660*pow2(lmMst1))*pow2(MuSUSY)))*pow4(Msq))+375*pow2(Mgl)*(80*(Dmsqst1*pow2(Msq)*(-3*pow2(Mst2)*((74-12*lmMst1)*pow2(Mst1)+7*pow2(MuSUSY))+
Dmst12*((98-24*lmMst1)*pow2(Mst1)+
21*pow2(MuSUSY)))+
pow2(Dmsqst1)*(-(pow2(Mst2)*((170-24*lmMst1)*pow2(Mst1)+21*pow2(MuSUSY)))
+Dmst12*(2*(23-6*lmMst1)*pow2(Mst1)+21*pow2(MuSUSY))))-
(Dmst12*(176*(398+47*lmMst1-48*pow2(lmMst1))*pow2(Mst1)+
3*(10473-40*lmMst1-256*pow2(lmMst1))*pow2(MuSUSY))+
8*pow2(Mst2)*(8*(347-50*lmMst1+183*pow2(lmMst1))*pow2(Mst1)+
(1361+10*lmMst1+54*pow2(lmMst1))*pow2(MuSUSY)))*pow4(Msq))))/(121500.*Mst1))/(pow2(Mgl)*pow4(Msq)*pow4(Mst2)))-pow2(Mt)*((shiftst3*(pow2(Dmst12)*((7-6*lmMst1)*pow2(MuSUSY)*pow2(s2t)*(1-pow2(Sbeta))+
6*(8*(2-lmMst1)*pow2(Mt)-
(15-14*lmMst1)*pow2(Mst1)*pow2(s2t))*pow2(Sbeta))-
2*Dmst12*pow2(Mst2)*((1-2*lmMst1)*pow2(MuSUSY)*pow2(s2t)*(1-pow2(Sbeta))+
12*((3-2*lmMst1)*pow2(Mt)-
(1-2*lmMst1)*pow2(Mst1)*pow2(s2t))*pow2(Sbeta))+
24*(1-2*lmMst1)*pow2(Mt)*pow2(Sbeta)*pow4(Mst2)))/(18.*pow2(Sbeta)*pow4(Mst2))+
(pow2(MuSUSY)*(-120*pow2(Dmst12)*pow2(s2t)*(14290-4540*lmMst1+
(5*Dmglst1*(84209+1264*lmMst1-240*pow2(lmMst1)))/Mgl+
240*pow2(lmMst1)+
(pow2(Dmglst1)*(1732531+16896*lmMst1-
24840*pow2(lmMst1)))/pow2(Mgl)-
(300*Dmsqst1*(7+24*lmMst1)*(Dmsqst1+pow2(Msq)))/pow4(Msq))+((240*Dmst12*Mt*s2t*(pow2(Dmglst1)*(2000*(pow2(Dmsqst1)+Dmsqst1*pow2(Msq))*(Dmst12-pow2(Mst2))+
(Dmst12*(958501+24456*lmMst1-11520*pow2(lmMst1))
-2*(1286791+5936*lmMst1-18120*pow2(lmMst1))*pow2(Mst2))*pow4(Msq))
+Dmglst1*Mgl*(2000*(pow2(Dmsqst1)+Dmsqst1*pow2(Msq))*(Dmst12-pow2(Mst2))-
(Dmst12*(949861-1944*lmMst1-11520*pow2(lmMst1))+
20*(33261-532*lmMst1-660*pow2(lmMst1))*pow2(Mst2))*pow4(Msq))+
5*pow2(Mgl)*(1680*(pow2(Dmsqst1)+Dmsqst1*pow2(Msq))*(Dmst12-pow2(Mst2))-
(3*Dmst12*(10473-40*lmMst1-256*pow2(lmMst1))+
8*(1361+10*lmMst1+54*pow2(lmMst1))*pow2(Mst2))*pow4(Msq))))/(Mst1*pow4(Msq))+
(pow2(Mt)*(16*pow2(Dmglst1)*((3891491+27200*lmMst1-19200*pow2(lmMst1))*(9*pow2(Dmst12)-9*Dmst12*pow2(Mst2))-
50*(345581+4896*lmMst1-3456*pow2(lmMst1))*pow4(Mst2))+
400*Dmglst1*Mgl*((403559+384*lmMst1-4608*pow2(lmMst1))*(pow2(Dmst12)-Dmst12*pow2(Mst2))-
24*(9631+16*lmMst1-192*pow2(lmMst1))*pow4(Mst2))-15*pow2(Mgl)*(pow2(Dmst12)*(852541+9216*lmMst1+6144*pow2(lmMst1))+
160*Dmst12*(11389-704*lmMst1-384*pow2(lmMst1))*pow2(Mst2)+
1920*(349-56*lmMst1-32*pow2(lmMst1))*pow4(Mst2))))/pow2(Mst1))/pow2(Mgl)))/(388800.*pow2(Sbeta)*pow4(Mst2))+
(5*(shiftst1*(24*pow2(Mt)*(2*(2-lmMst1)*pow2(Dmsqst1)+
Dmsqst1*(3-2*lmMst1)*pow2(Msq)+
(1-2*lmMst1)*pow4(Msq))+
(Dmst12*pow2(s2t)*(pow2(Dmsqst1)*(4*(2-lmMst1)*pow2(Mst2)*(-(pow2(MuSUSY)*(1-pow2(Sbeta)))+
12*pow2(Mst1)*pow2(Sbeta))-
Dmst12*(-((5-2*lmMst1)*pow2(MuSUSY)*(1-pow2(Sbeta)))+
36*(2-lmMst1)*pow2(Mst1)*pow2(Sbeta)))+
(2*pow2(Mst2)*(-(pow2(MuSUSY)*(1-pow2(Sbeta)))+
12*pow2(Mst1)*pow2(Sbeta))-
Dmst12*(-(pow2(MuSUSY)*(1-pow2(Sbeta)))+
18*pow2(Mst1)*pow2(Sbeta)))*(Dmsqst1*(3-2*lmMst1)*pow2(Msq)+
(1-2*lmMst1)*pow4(Msq))))/(pow2(Sbeta)*pow4(Mst2)))
+(shiftst2*((Dmst12*pow2(s2t)*(-(pow2(Dmsqst1)*(4*(2-lmMst1)*pow2(Mst2)*(-(pow2(MuSUSY)*(1-pow2(Sbeta)))+
12*pow2(Mst1)*pow2(Sbeta))-
Dmst12*((1-2*lmMst1)*pow2(MuSUSY)*(1-pow2(Sbeta))+
12*(2-lmMst1)*pow2(Mst1)*pow2(Sbeta))))+
(Dmst12*(pow2(MuSUSY)*(1-pow2(Sbeta))+
6*pow2(Mst1)*pow2(Sbeta))-
2*pow2(Mst2)*(-(pow2(MuSUSY)*(1-pow2(Sbeta)))+
12*pow2(Mst1)*pow2(Sbeta)))*(Dmsqst1*(3-2*lmMst1)*pow2(Msq)+
(1-2*lmMst1)*pow4(Msq))))/pow2(Sbeta)+
24*pow2(Mt)*(pow2(Mst2)*(Dmst12+pow2(Mst2))*(Dmsqst1*(3-2*lmMst1)*pow2(Msq)+
(1-2*lmMst1)*pow4(Msq))-
pow2(Dmsqst1)*(pow2(Dmst12)-
2*(2-lmMst1)*(Dmst12*pow2(Mst2)+pow4(Mst2))))))/pow4(Mst2)))/(9.*pow4(Msq)))-
(4984/81.+(8528*lmMst1)/81.+
(8*Dmglst1*(15571+508*lmMst1-978*pow2(lmMst1)))/(81.*Mgl)-
(1768*pow2(lmMst1))/27.+
(88*pow2(Dmglst1)*(293068+9168*lmMst1-7245*pow2(lmMst1)))/(6075.*pow2(Mgl))+(224*pow3(lmMst1))/9.+
((-8*Dmsqst1*(50*Dmglst1*(146-15*lmMst1)*Mgl+
(8014-435*lmMst1)*pow2(Dmglst1)+
225*(14-18*lmMst1+3*pow2(lmMst1))*pow2(Mgl)))/(405.*pow2(Msq))-(2*pow2(Dmsqst1)*(500*Dmglst1*(226-21*lmMst1)*Mgl+
2*(64033-3360*lmMst1)*pow2(Dmglst1)+
125*(205-606*lmMst1+54*pow2(lmMst1))*pow2(Mgl)))/(2025.*pow4(Msq))+(-(Dmst12*(pow2(Dmglst1)*(180*(13147-90*lmMst1)*(pow2(Dmsqst1)+Dmsqst1*pow2(Msq))-
(125277461+138180*lmMst1-153000*pow2(lmMst1))*pow4(Msq))+
10*Dmglst1*Mgl*(540*(423-20*lmMst1)*(pow2(Dmsqst1)+Dmsqst1*pow2(Msq))-
(5136871-107304*lmMst1-86040*pow2(lmMst1))*pow4(Msq))-
375*pow2(Mgl)*(45*(49+46*lmMst1)*pow2(Dmsqst1)+
10*Dmsqst1*(79+204*lmMst1)*pow2(Msq)+
2*(16443-524*lmMst1+264*pow2(lmMst1))*pow4(Msq))))/(30375.*pow2(Mst2))+(pow2(Dmst12)*(3*pow2(Dmglst1)*(11760*(13147-90*lmMst1)*(pow2(Dmsqst1)+Dmsqst1*pow2(Msq))-
(5259064991+6285660*lmMst1-
148327200*pow2(lmMst1))*pow4(Msq))+
98*Dmglst1*Mgl*(10800*(423-20*lmMst1)*(pow2(Dmsqst1)+Dmsqst1*pow2(Msq))-
(22574599-400620*lmMst1-1598400*pow2(lmMst1))*pow4(Msq))-
980*pow2(Mgl)*(90*(419-60*lmMst1)*pow2(Dmsqst1)-
15*Dmsqst1*(4561+510*lmMst1)*pow2(Msq)-
(93973+126198*lmMst1-48420*pow2(lmMst1))*pow4(Msq))))/(5.9535e6*pow4(Mst2)))/pow4(Msq)-
(pow2(MuSUSY)*(16*pow2(Dmglst1)*((3891491+27200*lmMst1-19200*pow2(lmMst1))*(9*pow2(Dmst12)-9*Dmst12*pow2(Mst2))-
50*(345581+4896*lmMst1-3456*pow2(lmMst1))*pow4(Mst2))
+400*Dmglst1*Mgl*((403559+384*lmMst1-4608*pow2(lmMst1))*(pow2(Dmst12)-Dmst12*pow2(Mst2))-
24*(9631+16*lmMst1-192*pow2(lmMst1))*pow4(Mst2))-
15*pow2(Mgl)*(pow2(Dmst12)*(852541+9216*lmMst1+6144*pow2(lmMst1))+
160*Dmst12*(11389-704*lmMst1-384*pow2(lmMst1))*pow2(Mst2)+
1920*(349-56*lmMst1-32*pow2(lmMst1))*pow4(Mst2))))/(388800.*pow2(Mst1)*pow4(Mst2)))/pow2(Mgl))*pow4(Mt)-
(z2*((pow2(Mt)*(-(shiftst3*(3*pow2(Dmst12)*(pow2(MuSUSY)*pow2(s2t)*(1-pow2(Sbeta))+
2*(4*pow2(Mt)-7*pow2(Mst1)*pow2(s2t))*pow2(Sbeta))-
2*Dmst12*pow2(Mst2)*(pow2(MuSUSY)*pow2(s2t)*(1-pow2(Sbeta))+
12*(pow2(Mt)-pow2(Mst1)*pow2(s2t))*pow2(Sbeta))
+24*pow2(Mt)*pow2(Sbeta)*pow4(Mst2)))/9.-
(10*(pow2(Dmsqst1)+Dmsqst1*pow2(Msq)+pow4(Msq))*(shiftst2*(pow2(Dmst12)*pow2(s2t)*(pow2(MuSUSY)*(1-pow2(Sbeta))+
6*pow2(Mst1)*pow2(Sbeta))+
2*Dmst12*pow2(Mst2)*(pow2(MuSUSY)*pow2(s2t)*(1-pow2(Sbeta))+
12*(pow2(Mt)-pow2(Mst1)*pow2(s2t))*pow2(Sbeta))+24*pow2(Mt)*pow2(Sbeta)*pow4(Mst2))+
shiftst1*(pow2(s2t)*(2*Dmst12*pow2(Mst2)*(-(pow2(MuSUSY)*(1-pow2(Sbeta)))+
12*pow2(Mst1)*pow2(Sbeta))-
pow2(Dmst12)*(-(pow2(MuSUSY)*(1-pow2(Sbeta)))+
18*pow2(Mst1)*pow2(Sbeta)))+
24*pow2(Mt)*pow2(Sbeta)*pow4(Mst2))))/(9.*pow4(Msq))))/pow2(Sbeta)+(Mt*(-(MuSUSY*s2t*(10*(pow2(Dmsqst1)+Dmsqst1*pow2(Msq))*(-((shiftst1-shiftst2)*pow2(Dmst12)*pow2(Mst1)*pow2(s2t))+
pow2(Mt)*(-8*Dmst12*shiftst2*pow2(Mst2)+
8*(shiftst1-shiftst2)*pow4(Mst2)))+
pow4(Msq)*(pow2(Dmst12)*(8*shiftst3*pow2(Mt)-
(10*shiftst1-10*shiftst2+shiftst3)*pow2(Mst1)*pow2(s2t))+
pow2(Mt)*(-8*Dmst12*(10*shiftst2+shiftst3)*pow2(Mst2)+
8*(10*shiftst1-10*shiftst2+shiftst3)*pow4(Mst2)))))/(3.*pow4(Msq))+
(10*xDmsqst1*pow3(Dmsqst1)*(-(pow2(Dmst12)*pow2(s2t)*(Mt*(shiftst1+shiftst2)*Tbeta*pow2(MuSUSY)*(1-pow2(Sbeta))-
3*(MuSUSY*s2t*(shiftst1-shiftst2)+
2*Mt*(3*shiftst1-shiftst2)*Tbeta)*pow2(Mst1)*pow2(Sbeta)))+
2*Dmst12*Mt*pow2(Mst2)*(12*Mt*MuSUSY*s2t*shiftst2*pow2(Sbeta)+
Tbeta*(-12*shiftst2*pow2(Mt)*pow2(Sbeta)-
(shiftst1-shiftst2)*pow2(s2t)*(-(pow2(MuSUSY)*(1-pow2(Sbeta)))+
12*pow2(Mst1)*pow2(Sbeta))))-
24*(MuSUSY*s2t*(shiftst1-shiftst2)+
Mt*(shiftst1+shiftst2)*Tbeta)*pow2(Mt)*pow2(Sbeta)*pow4(Mst2)))/(9.*pow2(Sbeta)*pow6(Msq))))/Tbeta))/pow4(Mst2)-(-(MuSUSY*(225*Mst1*pow2(Dmst12)*pow2(Mt)*pow2(s2t)*(27220+200*lmMst1+
(10*Dmglst1*(33261-532*lmMst1-660*pow2(lmMst1)))/Mgl+
1080*pow2(lmMst1)+
(pow2(Dmglst1)*(1286791+5936*lmMst1-18120*pow2(lmMst1)))/pow2(Mgl)
+(200*Dmsqst1*(5*(Dmglst1*Mgl+pow2(Dmglst1))+21*pow2(Mgl))*(Dmsqst1+pow2(Msq)))/(pow2(Mgl)*pow4(Msq)))+
20250*Mt*s2t*shiftst3*(-(pow2(Dmst12)*(16*(2-lmMst1)*pow2(Mt)-
(1-2*lmMst1)*pow2(Mst1)*pow2(s2t)))+
pow2(Mt)*(8*Dmst12*(3-2*lmMst1)*pow2(Mst2)-
8*(1-2*lmMst1)*pow4(Mst2)))+
(-202500*shiftst1*(2*(2-lmMst1)*pow2(Dmsqst1)+
Dmsqst1*(3-2*lmMst1)*pow2(Msq)+
(1-2*lmMst1)*pow4(Msq))*(-(Mt*pow2(Dmst12)*pow2(Mst1)*pow3(s2t))+
8*s2t*pow3(Mt)*pow4(Mst2))+
202500*Mt*s2t*shiftst2*((Dmsqst1*(3-2*lmMst1)*pow2(Msq)+
(1-2*lmMst1)*pow4(Msq))*(-(pow2(Dmst12)*pow2(Mst1)*pow2(s2t))+
8*pow2(Mt)*(Dmst12*pow2(Mst2)+pow4(Mst2)))+
2*pow2(Dmsqst1)*(-(pow2(Dmst12)*(4*pow2(Mt)+(2-lmMst1)*pow2(Mst1)*pow2(s2t)))+8*(2-lmMst1)*pow2(Mt)*(Dmst12*pow2(Mst2)+pow4(Mst2)))))/pow4(Msq)
+(450*Dmst12*s2t*pow3(Mt)*(4*Dmglst1*Mgl*(4700*(pow2(Dmsqst1)+Dmsqst1*pow2(Msq))*(Dmst12-pow2(Mst2))-
(Dmst12*(17539+574*lmMst1-1920*pow2(lmMst1))-
5*(4539-596*lmMst1-516*pow2(lmMst1))*pow2(Mst2))*pow4(Msq))+
pow2(Dmglst1)*(19160*(pow2(Dmsqst1)+Dmsqst1*pow2(Msq))*(Dmst12-pow2(Mst2))-
(Dmst12*(586073+9268*lmMst1-19200*pow2(lmMst1))-
2*(13289-916*lmMst1-60*pow2(lmMst1))*pow2(Mst2))*pow4(Msq))+
5*pow2(Mgl)*(5*Dmsqst1*pow2(Msq)*(161*Dmst12+24*(1+12*lmMst1)*pow2(Mst2))-
5*pow2(Dmsqst1)*(44*Dmst12-(229+288*lmMst1)*pow2(Mst2))+
(Dmst12*(2071+296*lmMst1-600*pow2(lmMst1))+
8*(631-lmMst1+21*pow2(lmMst1))*pow2(Mst2))*pow4(Msq)))+
((15*pow2(Mgl)*(-(pow4(Msq)*(pow2(Dmst12)*(2284899-57376*lmMst1-87360*pow2(lmMst1))
+200*Dmst12*(11697+448*lmMst1+408*pow2(lmMst1))*pow2(Mst2)+
1600*(434-83*lmMst1+183*pow2(lmMst1))*pow4(Mst2)))+
4000*(pow2(Dmsqst1)*((65-6*lmMst1)*(pow2(Dmst12)-Dmst12*pow2(Mst2))-
(91-12*lmMst1)*pow4(Mst2))+
Dmsqst1*pow2(Msq)*((65-6*lmMst1)*(pow2(Dmst12)-Dmst12*pow2(Mst2))-
6*(20-3*lmMst1)*pow4(Mst2))))+
2*Dmglst1*Mgl*(-(pow4(Msq)*(21*pow2(Dmst12)*(5639279-7540*lmMst1-45600*pow2(lmMst1))
+40*Dmst12*(3746977-48798*lmMst1-52380*pow2(lmMst1))*pow2(Mst2)+
6000*(9961-670*lmMst1+42*pow2(lmMst1))*pow4(Mst2)))+
1200*(Dmsqst1*pow2(Msq)*((557+120*lmMst1)*(pow2(Dmst12)-Dmst12*pow2(Mst2))-
50*(26+3*lmMst1)*pow4(Mst2))+
pow2(Dmsqst1)*((557+120*lmMst1)*(pow2(Dmst12)-Dmst12*pow2(Mst2))-
25*(44+3*lmMst1)*pow4(Mst2))))+
2*pow2(Dmglst1)*(pow4(Msq)*(pow2(Dmst12)*(387579543+2205300*lmMst1-
4010400*pow2(lmMst1))-
2*Dmst12*(327941741+47520*lmMst1-
3531600*pow2(lmMst1))*pow2(Mst2)-
80*(3777727-57348*lmMst1-
52380*pow2(lmMst1))*pow4(Mst2))+
1200*(pow2(Dmsqst1)*((557+120*lmMst1)*(pow2(Dmst12)-Dmst12*pow2(Mst2))+
(136-165*lmMst1)*pow4(Mst2))+
Dmsqst1*pow2(Msq)*((557+120*lmMst1)*(pow2(Dmst12)-Dmst12*pow2(Mst2))-
16*(4+15*lmMst1)*pow4(Mst2)))))*pow4(Mt))/Mst1)/(pow2(Mgl)*pow4(Msq))))/(121500.*Tbeta)+
((xDmglst1*pow2(Mt)*pow3(Dmglst1)*(11907000*pow2(Dmst12)*pow2(s2t)*(((3044017+27472*lmMst1-48480*pow2(lmMst1))*pow2(MuSUSY))/3240.+
(pow2(Mst1)*(109771+2196*lmMst1-
3816*pow2(lmMst1)+
(2928*Dmsqst1*(Dmsqst1+pow2(Msq)))/pow4(Msq)))/81.)+((-2*Dmst12*Mt*s2t*(11760*Dmsqst1*pow2(Msq)*(pow2(Mst2)*(6*(791-270*lmMst1)*pow2(Mst1)-
625*pow2(MuSUSY))-
Dmst12*(4*(2729-105*lmMst1)*pow2(Mst1)-
625*pow2(MuSUSY)))-
11760*pow2(Dmsqst1)*(-(pow2(Mst2)*((7121-870*lmMst1)*pow2(Mst1)-
625*pow2(MuSUSY)))+
Dmst12*((13291+330*lmMst1)*pow2(Mst1)-
625*pow2(MuSUSY)))-
(Dmst12*((103099497724+748797600*lmMst1-
310363200*pow2(lmMst1))*pow2(Mst1)-
245*(99874229+1352280*lmMst1-
633600*pow2(lmMst1))*pow2(MuSUSY))+
196*pow2(Mst2)*(8*(61021241+307815*lmMst1-
675900*pow2(lmMst1))*pow2(Mst1)+
5*(31025111+290880*lmMst1-
251100*pow2(lmMst1))*pow2(MuSUSY)))*pow4(Msq)))/Mst1+(2*pow2(Mt)*(pow4(Msq)*(-(pow2(Dmst12)*(4*(7672052891-14084100*lmMst1-
98506800*pow2(lmMst1))*pow2(Mst1)+
245*(59957863+480000*lmMst1-
230400*pow2(lmMst1))*pow2(MuSUSY)))+
Dmst12*pow2(Mst2)*(16*(2524164367+8198205*lmMst1+
27997200*pow2(lmMst1))*pow2(Mst1)+
245*(59957863+480000*lmMst1-
230400*pow2(lmMst1))*pow2(MuSUSY))+
588*(4*(21126629+579160*lmMst1-
194100*pow2(lmMst1))*pow2(Mst1)+
5*(3877891+46400*lmMst1-
19200*pow2(lmMst1))*pow2(MuSUSY))*pow4(Mst2))
+pow2(Mst1)*(47040*Dmsqst1*pow2(Msq)*((3401+105*lmMst1)*(3*pow2(Dmst12)-3*Dmst12*pow2(Mst2))-
(19241-420*lmMst1)*pow4(Mst2))+
23520*pow2(Dmsqst1)*((3401+105*lmMst1)*(6*pow2(Dmst12)-6*Dmst12*pow2(Mst2))-
25*(1225-39*lmMst1)*pow4(Mst2)))))/pow2(Mst1))/pow4(Msq)+((-245*pow2(MuSUSY)*(2*Dmst12*Mt*pow2(Mst2)*(30000*Mst1*s2t*(pow2(Dmsqst1)+Dmsqst1*pow2(Msq))+
(4*Mst1*s2t*(31025111+290880*lmMst1-
251100*pow2(lmMst1))+
Mt*(59957863+480000*lmMst1-
230400*pow2(lmMst1)))*pow4(Msq))-
pow2(Dmst12)*(60000*Mst1*Mt*s2t*(pow2(Dmsqst1)+Dmsqst1*pow2(Msq))+
(2*(Mst1*Mt*s2t*(99874229+1352280*lmMst1-
633600*pow2(lmMst1))+
(59957863+480000*lmMst1-
230400*pow2(lmMst1))*pow2(Mt))-
15*(3044017+27472*lmMst1-
48480*pow2(lmMst1))*pow2(Mst1)*pow2(s2t))*pow4(Msq))+
24*(3877891+46400*lmMst1-19200*pow2(lmMst1))*pow2(Mt)*pow4(Msq)*pow4(Mst2)))/(pow2(Mst1)*pow2(Sbeta))-
(4*MuSUSY*(-(pow4(Msq)*(2*Dmst12*Mt*(Mt*(49723877243+296163780*lmMst1-
342543600*pow2(lmMst1))+
7350*Mst1*s2t*(548999+10980*lmMst1-19080*pow2(lmMst1)))*pow2(Mst2)-
15*pow2(Dmst12)*(-490*Mst1*Mt*s2t*(611423-9984*lmMst1-23040*pow2(lmMst1))+
5*(1150678153+9276404*lmMst1-
7140000*pow2(lmMst1))*pow2(Mt)+
49*(31025111+290880*lmMst1-
251100*pow2(lmMst1))*pow2(Mst1)*pow2(s2t))+
392*(122282257+499080*lmMst1-
1351800*pow2(lmMst1))*pow2(Mt)*pow4(Mst2)))+
2940*(pow2(Dmsqst1)*(-20*Dmst12*Mt*((557+120*lmMst1)*Mt+3660*Mst1*s2t)*pow2(Mst2)+
5*pow2(Dmst12)*(14640*Mst1*Mt*s2t+
4*(557+120*lmMst1)*pow2(Mt)+
375*pow2(Mst1)*pow2(s2t))+
4*(3778-435*lmMst1)*pow2(Mt)*pow4(Mst2))+
Dmsqst1*pow2(Msq)*(-20*Dmst12*Mt*((557+120*lmMst1)*Mt+3660*Mst1*s2t)*pow2(Mst2)+
5*pow2(Dmst12)*(14640*Mst1*Mt*s2t+
4*(557+120*lmMst1)*pow2(Mt)+
375*pow2(Mst1)*pow2(s2t))+
24*(463-135*lmMst1)*pow2(Mt)*pow4(Mst2)))))/(Mst1*Tbeta))/pow4(Msq)))/1.1907e7+
(Mt*xDmsqst1*pow3(Dmsqst1)*(-4*Mt*xDmglst1*pow3(Dmglst1)*(-2*Dmst12*Mt*pow2(Mst2)*(2*Mt*(5*(557+120*lmMst1)*MuSUSY-
6*(3401+105*lmMst1)*Mst1*Tbeta)*pow2(Sbeta)+
s2t*(-625*Tbeta*pow2(MuSUSY)*(1-pow2(Sbeta))+
(36600*Mst1*MuSUSY-
8*(1187-15*lmMst1)*Tbeta*pow2(Mst1))*pow2(Sbeta)))+
pow2(Dmst12)*((4*(5*(557+120*lmMst1)*MuSUSY-
6*(3401+105*lmMst1)*Mst1*Tbeta)*pow2(Mt)+
75*(25*MuSUSY-488*Mst1*Tbeta)*pow2(Mst1)*pow2(s2t))*pow2(Sbeta)+
2*Mt*s2t*(-625*Tbeta*pow2(MuSUSY)*(1-pow2(Sbeta))+
(36600*Mst1*MuSUSY-
6*(2611+180*lmMst1)*Tbeta*pow2(Mst1))*pow2(Sbeta)))+
8*((2389-30*lmMst1)*MuSUSY+
(11384-555*lmMst1)*Mst1*Tbeta)*pow2(Mt)*pow2(Sbeta)*pow4(Mst2))+
Mt*(-4*Mgl*pow2(Dmglst1)*(-2*Dmst12*Mt*pow2(Mst2)*(Mt*(10*(557+120*lmMst1)*MuSUSY-
3*(13147-90*lmMst1)*Mst1*Tbeta)*pow2(Sbeta)
+5*s2t*(-125*Tbeta*pow2(MuSUSY)*(1-pow2(Sbeta))+
(7185*Mst1*MuSUSY-
6*(97-30*lmMst1)*Tbeta*pow2(Mst1))*pow2(Sbeta)))+
pow2(Dmst12)*((2*(10*(557+120*lmMst1)*MuSUSY-
3*(13147-90*lmMst1)*Mst1*Tbeta)*pow2(Mt)+
75*(25*MuSUSY-479*Mst1*Tbeta)*pow2(Mst1)*pow2(s2t))*pow2(Sbeta)+
10*Mt*s2t*(-125*Tbeta*pow2(MuSUSY)*(1-pow2(Sbeta))+
(7185*Mst1*MuSUSY-
4*(454+15*lmMst1)*Tbeta*pow2(Mst1))*pow2(Sbeta)))+
4*(30*(56-15*lmMst1)*MuSUSY+
(23963-1185*lmMst1)*Mst1*Tbeta)*pow2(Mt)*pow2(Sbeta)*pow4(Mst2))+
20*Dmglst1*pow2(Mgl)*(2*Dmst12*Mt*pow2(Mst2)*(2*Mt*((557+120*lmMst1)*MuSUSY-
9*(423-20*lmMst1)*Mst1*Tbeta)*pow2(Sbeta)+
25*s2t*(-5*Tbeta*pow2(MuSUSY)*(1-pow2(Sbeta))+
(282*Mst1*MuSUSY+72*Tbeta*pow2(Mst1))*pow2(Sbeta)))-
pow2(Dmst12)*((4*((557+120*lmMst1)*MuSUSY-
9*(423-20*lmMst1)*Mst1*Tbeta)*pow2(Mt)+
75*(5*MuSUSY-94*Mst1*Tbeta)*pow2(Mst1)*pow2(s2t))*pow2(Sbeta)+
2*Mt*s2t*(-125*Tbeta*pow2(MuSUSY)*(1-pow2(Sbeta))+
(7050*Mst1*MuSUSY+
2*(283-120*lmMst1)*Tbeta*pow2(Mst1))*pow2(Sbeta)))+
400*(9*MuSUSY-(40-3*lmMst1)*Mst1*Tbeta)*pow2(Mt)*pow2(Sbeta)*pow4(Mst2)))-
5*pow3(Mgl)*(pow2(Dmst12)*(-75*s2t*pow2(Mt)*(56*Tbeta*pow2(MuSUSY)*(1-pow2(Sbeta))+
(3*Mst1*MuSUSY*(83+96*shiftst2)+
16*Tbeta*pow2(Mst1))*pow2(Sbeta))+
75*Mst1*Mt*pow2(s2t)*(-((7+42*shiftst1+
12*lmMst1*(2-shiftst1-shiftst2)-
6*shiftst2)*Tbeta*pow2(MuSUSY)*(1-pow2(Sbeta)))+
(84*Mst1*MuSUSY+
(217+492*shiftst1-132*shiftst2+
72*lmMst1*(2-3*shiftst1+shiftst2))*Tbeta*pow2(Mst1))*pow2(Sbeta))+
2*(200*(65-6*lmMst1)*MuSUSY+
Mst1*(9589-210*lmMst1+10800*shiftst2)*Tbeta)*pow2(Sbeta)*pow3(Mt)+
1350*(5-2*lmMst1)*MuSUSY*(shiftst1-shiftst2)*pow2(Sbeta)*pow3(Mst1)*pow3(s2t))+
Mt*(-50*Dmst12*pow2(Mst2)*(4*(2*(65-6*lmMst1)*MuSUSY+
Mst1*(181+3*lmMst1*(35-36*shiftst2)+
270*shiftst2)*Tbeta)*pow2(Mt)*pow2(Sbeta)+
18*Mst1*(shiftst1-shiftst2)*Tbeta*pow2(s2t)*(-((5-2*lmMst1)*pow2(MuSUSY)*(1-pow2(Sbeta)))+
8*(7-3*lmMst1)*pow2(Mst1)*pow2(Sbeta))-
Mt*s2t*(84*Tbeta*pow2(MuSUSY)*(1-pow2(Sbeta))+
(3*Mst1*MuSUSY*(217+144*lmMst1*(1-shiftst2)+
360*shiftst2)-
8*(59-6*lmMst1)*Tbeta*pow2(Mst1))*pow2(Sbeta)))-200*pow2(Sbeta)*(36*(7-3*lmMst1)*Mst1*Mt*MuSUSY*s2t*(shiftst1-shiftst2)+
(4*(31-3*lmMst1)*MuSUSY+
Mst1*(47+252*(shiftst1+shiftst2)+
6*lmMst1*(47-18*(shiftst1+shiftst2)))*Tbeta)*pow2(Mt)+18*Tbeta*pow3(Mst1))*pow4(Mst2)))))/(4050.*Mst1*Tbeta*pow2(Sbeta)*pow6(Msq)))/pow3(Mgl))/pow4(Mst2)+
(xDmst12*pow3(Dmst12)*((245*pow2(MuSUSY)*pow3(Mt)*((Mt*(72*pow2(Dmglst1)*(3891491+27200*lmMst1-19200*pow2(lmMst1))+
200*Dmglst1*Mgl*(403559+384*lmMst1-4608*pow2(lmMst1))-
15*(1763661-47104*lmMst1-24576*pow2(lmMst1))*pow2(Mgl)))/pow2(Mgl)-
240*Mst1*s2t*(184315-400*lmMst1+
(Dmglst1*(1282471-7264*lmMst1-18120*pow2(lmMst1)))/Mgl-2760*pow2(lmMst1)+
(10*pow2(Dmglst1)*(32829-1852*lmMst1-660*pow2(lmMst1)))/pow2(Mgl)-
(200*Dmsqst1*(5*(Dmglst1*Mgl+pow2(Dmglst1))+
21*pow2(Mgl))*(Dmsqst1+pow2(Msq)))/(pow2(Mgl)*pow4(Msq)))))/(pow2(Mst1)*pow2(Sbeta))+
29400*Mt*pow3(Mst1)*pow3(s2t)*(27220+200*lmMst1+(10*Dmglst1*(33261-532*lmMst1-660*pow2(lmMst1)))/Mgl+
1080*pow2(lmMst1)+(pow2(Dmglst1)*(1286791+5936*lmMst1-18120*pow2(lmMst1)))/pow2(Mgl)+
(200*Dmsqst1*(5*(Dmglst1*Mgl+pow2(Dmglst1))+21*pow2(Mgl))*(Dmsqst1+pow2(Msq)))/(pow2(Mgl)*pow4(Msq)))+
((392*s2t*pow3(Mt)*(-2*Dmglst1*Mgl*(-600*pow2(Dmsqst1)*((193+330*lmMst1)*pow2(Mst1)-125*pow2(MuSUSY))
+600*Dmsqst1*pow2(Msq)*(6*(47-30*lmMst1)*pow2(Mst1)+
125*pow2(MuSUSY))+
((28188929-2075220*lmMst1+
1389600*pow2(lmMst1))*pow2(Mst1)-
75*(1282471-7264*lmMst1-18120*pow2(lmMst1))*pow2(MuSUSY))*pow4(Msq))+
10*pow2(Dmglst1)*(3000*(Dmsqst1*pow2(Msq)*(84*pow2(Mst1)-5*pow2(MuSUSY))+
pow2(Dmsqst1)*((103+6*lmMst1)*pow2(Mst1)-5*pow2(MuSUSY)))
+((148185343+1333716*lmMst1-1376640*pow2(lmMst1))*pow2(Mst1)+
150*(32829-1852*lmMst1-660*pow2(lmMst1))*pow2(MuSUSY))*pow4(Msq))+
15*pow2(Mgl)*(2000*Dmsqst1*pow2(Msq)*(2*(13+6*lmMst1)*pow2(Mst1)-21*pow2(MuSUSY))+
6000*pow2(Dmsqst1)*(26*pow2(Mst1)-7*pow2(MuSUSY))-
((558619-273056*lmMst1+123840*pow2(lmMst1))*pow2(Mst1)-
50*(36863-80*lmMst1-552*pow2(lmMst1))*pow2(MuSUSY))*pow4(Msq))))/Mst1+
3675*pow2(Mt)*pow2(s2t)*(8*pow2(Dmglst1)*(114960*(pow2(Dmsqst1)+Dmsqst1*pow2(Msq))*pow2(Mst1)-
(3*(612347+7436*lmMst1-19320*pow2(lmMst1))*pow2(Mst1)-
(1732531+16896*lmMst1-24840*pow2(lmMst1))*pow2(MuSUSY))*pow4(Msq))-
5*pow2(Mgl)*(-120*Dmsqst1*pow2(Msq)*((137-288*lmMst1)*pow2(Mst1)-
4*(7+24*lmMst1)*pow2(MuSUSY))+
120*pow2(Dmsqst1)*(3*(91+96*lmMst1)*pow2(Mst1)+
4*(7+24*lmMst1)*pow2(MuSUSY))+
(24*(2785-304*lmMst1+768*pow2(lmMst1))*pow2(Mst1)+
(70121+2208*lmMst1-432*pow2(lmMst1))*pow2(MuSUSY))*pow4(Msq))+8*Dmglst1*Mgl*(112800*(pow2(Dmsqst1)+Dmsqst1*pow2(Msq))*pow2(Mst1)-
(24*(20017-1203*lmMst1-2250*pow2(lmMst1))*pow2(Mst1)-
5*(84209+1264*lmMst1-240*pow2(lmMst1))*pow2(MuSUSY))*pow4(Msq))))/(pow2(Mgl)*pow4(Msq))+
pow2(s2t)*(47628000*pow2(Mt)*((pow2(MuSUSY)*(350605+11040*lmMst1-
(40*Dmglst1*(84209+1264*lmMst1-240*pow2(lmMst1)))/Mgl-
2160*pow2(lmMst1)-
(8*pow2(Dmglst1)*(1732531+16896*lmMst1-24840*pow2(lmMst1)))/pow2(Mgl)+
(2400*Dmsqst1*(7+24*lmMst1)*(Dmsqst1+pow2(Msq)))/pow4(Msq)))/(12960.*pow2(Sbeta))-
(shiftst1*(80*Dmsqst1*(3-2*lmMst1)*pow2(Msq)*pow2(Mst1)+
10*pow2(Dmsqst1)*(2*(15-8*lmMst1)*pow2(Mst1)+pow2(MuSUSY))+
(1-2*lmMst1)*(80*pow2(Mst1)+3*pow2(MuSUSY))*pow4(Msq)))/(9.*pow4(Msq)))+
5292000*shiftst1*((pow2(Mt)*pow2(MuSUSY)*(3-6*lmMst1+(10*pow2(Dmsqst1))/pow4(Msq)))/pow2(Sbeta)
+(5*pow2(s2t)*(2*(2-lmMst1)*pow2(Dmsqst1)+
Dmsqst1*(3-2*lmMst1)*pow2(Msq)+
(1-2*lmMst1)*pow4(Msq))*pow4(Mst1))/pow4(Msq)))+
(4088087836+849236640*lmMst1-454406400*pow2(lmMst1)+
(4704*Dmglst1*(9598037-224140*lmMst1+246000*pow2(lmMst1)))/Mgl-(80*pow2(Dmglst1)*(700000759+1063068*lmMst1-85997520*pow2(lmMst1)))/pow2(Mgl)+((94080*Dmsqst1*(90*Dmglst1*(423-20*lmMst1)*Mgl+
3*(13147-90*lmMst1)*pow2(Dmglst1)+
5*(3268+2805*lmMst1)*pow2(Mgl)))/pow2(Msq)-
(245*(72*pow2(Dmglst1)*(3891491+27200*lmMst1-19200*pow2(lmMst1))+
200*Dmglst1*Mgl*(403559+384*lmMst1-4608*pow2(lmMst1))-
15*(1763661-47104*lmMst1-24576*pow2(lmMst1))*pow2(Mgl))*pow2(MuSUSY))/pow2(Mst1)+
(70560*pow2(Dmsqst1)*(120*Dmglst1*(423-20*lmMst1)*Mgl+
4*(13147-90*lmMst1)*pow2(Dmglst1)+
5*(1999+3690*lmMst1)*pow2(Mgl)))/pow4(Msq))/pow2(Mgl))*pow4(Mt)+(529200*shiftst3*(-2*pow2(Mt)*pow2(s2t)*(-((29-18*lmMst1)*pow2(MuSUSY)*(1-pow2(Sbeta)))+
5*(83-58*lmMst1)*pow2(Mst1)*pow2(Sbeta))+
pow2(Sbeta)*(80*(7-3*lmMst1)*pow4(Mt)+
5*(1-2*lmMst1)*pow4(Mst1)*pow4(s2t)))+
(5292000*shiftst2*(-((1-2*lmMst1)*pow2(s2t)*pow4(Msq)*(-2*pow2(Mt)*(pow2(MuSUSY)*(1-pow2(Sbeta))+
10*pow2(Mst1)*pow2(Sbeta))+
5*pow2(s2t)*pow2(Sbeta)*pow4(Mst1)))+
5*Dmsqst1*(3-2*lmMst1)*pow2(Msq)*(4*pow2(Mt)*pow2(s2t)*(pow2(MuSUSY)*(1-pow2(Sbeta))+
pow2(Mst1)*pow2(Sbeta))+
pow2(Sbeta)*(24*pow4(Mt)-pow4(Mst1)*pow4(s2t)))+
10*pow2(Dmsqst1)*(-(pow2(Mt)*pow2(s2t)*(-((5-4*lmMst1)*pow2(MuSUSY)*(1-pow2(Sbeta)))+
4*lmMst1*pow2(Mst1)*pow2(Sbeta)))+
pow2(Sbeta)*(24*(1-lmMst1)*pow4(Mt)-
(2-lmMst1)*pow4(Mst1)*pow4(s2t)))))/pow4(Msq))/pow2(Sbeta)+(Mt*(196*MuSUSY*(13500*s2t*shiftst3*(16*(7-3*lmMst1)*pow2(Mt)-
(13-14*lmMst1)*pow2(Mst1)*pow2(s2t))+
225*Mst1*Mt*pow2(s2t)*(102655-1000*lmMst1-6000*pow2(lmMst1)+
(Dmglst1*(284641+8696*lmMst1+1680*pow2(lmMst1)))/Mgl-
(pow2(Dmglst1)*(3532083+36328*lmMst1-47760*pow2(lmMst1)))/pow2(Mgl)-
(800*Dmsqst1*(5*(Dmglst1*Mgl+pow2(Dmglst1))+21*pow2(Mgl))*(Dmsqst1+pow2(Msq)))/(pow2(Mgl)*pow4(Msq)))+
(4*pow3(Mt)*(51818985-188640*lmMst1+
(2*Dmglst1*(193364399-1134300*lmMst1-
2005200*pow2(lmMst1)))/Mgl-
698400*pow2(lmMst1)-
(4*pow2(Dmglst1)*(29818901+1078890*lmMst1-
239400*pow2(lmMst1)))/pow2(Mgl)-
(1200*Dmsqst1*((557+120*lmMst1)*(Dmglst1*Mgl+pow2(Dmglst1))+
25*(65-6*lmMst1)*pow2(Mgl))*(Dmsqst1+pow2(Msq)))/(pow2(Mgl)*pow4(Msq))))/Mst1-(225*s2t*pow2(Mt)*(pow2(Mgl)*(300*(47+96*lmMst1)*pow2(Dmsqst1)+
200*Dmsqst1*(173+144*lmMst1)*pow2(Msq)+
(102747+10592*lmMst1-13888*pow2(lmMst1))*pow4(Msq))+
32*pow2(Dmglst1)*(2395*(pow2(Dmsqst1)+Dmsqst1*pow2(Msq))-
(143196+2546*lmMst1-4785*pow2(lmMst1))*pow4(Msq))+
16*Dmglst1*Mgl*(4700*(pow2(Dmsqst1)+Dmsqst1*pow2(Msq))-
(12383+4128*lmMst1-1260*pow2(lmMst1))*pow4(Msq))))/(pow2(Mgl)*pow4(Msq))+
((-75*pow2(Mst1)*pow3(s2t)*(300*pow2(Mgl)*((7+24*lmMst1*(1-2*shiftst1)+
108*shiftst1)*pow2(Dmsqst1)+
Dmsqst1*(7+24*lmMst1*(1-2*shiftst1)+
72*shiftst1)*pow2(Msq))-
(pow2(Dmglst1)*(1732531+16896*lmMst1-24840*pow2(lmMst1))
+5*Dmglst1*Mgl*(84209+1264*lmMst1-240*pow2(lmMst1))+
10*(1429-454*lmMst1+
(-720+1440*lmMst1)*shiftst1+
24*pow2(lmMst1))*pow2(Mgl))*pow4(Msq)))/pow2(Mgl)+
270000*s2t*shiftst2*(-(Dmsqst1*(3-2*lmMst1)*pow2(Msq)*(12*pow2(Mt)-pow2(Mst1)*pow2(s2t)))-
pow2(Dmsqst1)*(24*(1-lmMst1)*pow2(Mt)-
2*(3-lmMst1)*pow2(Mst1)*pow2(s2t))+
(1-2*lmMst1)*pow2(Mst1)*pow2(s2t)*pow4(Msq)))/pow4(Msq))+
(4*xDmglst1*pow3(Dmglst1)*(11760*Mst1*(Dmsqst1*pow2(Msq)*(-2*s2t*pow2(Mt)*(-625*Tbeta*pow2(MuSUSY)*(1-pow2(Sbeta))+
(36600*Mst1*MuSUSY-
2*(8543+390*lmMst1)*Tbeta*pow2(Mst1))*pow2(Sbeta))+
pow2(Sbeta)*(-150*Mt*(25*MuSUSY-488*Mst1*Tbeta)*pow2(Mst1)*pow2(s2t)-
4*(5*(557+120*lmMst1)*MuSUSY-
6*(3401+105*lmMst1)*Mst1*Tbeta)*pow3(Mt)+
625*Tbeta*pow3(s2t)*pow4(Mst1)))+
pow2(Dmsqst1)*(-2*s2t*pow2(Mt)*(-625*Tbeta*pow2(MuSUSY)*(1-pow2(Sbeta))+
(36600*Mst1*MuSUSY-
3*(6487+510*lmMst1)*Tbeta*pow2(Mst1))*pow2(Sbeta))+
pow2(Sbeta)*(-150*Mt*(25*MuSUSY-488*Mst1*Tbeta)*pow2(Mst1)*pow2(s2t)-
4*(5*(557+120*lmMst1)*MuSUSY-
6*(3401+105*lmMst1)*Mst1*Tbeta)*pow3(Mt)+
625*Tbeta*pow3(s2t)*pow4(Mst1))))+
pow4(Msq)*(-4*Mst1*MuSUSY*pow2(Mt)*(-245*MuSUSY*s2t*Tbeta*(37824007+770520*lmMst1-
131400*pow2(lmMst1))*(1-pow2(Sbeta))+
16*Mt*(4572123029+49945815*lmMst1-
24119550*pow2(lmMst1))*pow2(Sbeta))+
6*Mt*pow2(Mst1)*(-1225*Tbeta*(3044017+27472*lmMst1-48480*pow2(lmMst1))*pow2(MuSUSY)*pow2(s2t)*(1-pow2(Sbeta))+
(19600*Mt*MuSUSY*s2t*(580211+498*lmMst1-21060*pow2(lmMst1))-
8*Tbeta*(874574719-10160170*lmMst1-
51500400*pow2(lmMst1))*pow2(Mt))*pow2(Sbeta))+
pow2(Sbeta)*(-(Mt*s2t*(-4*Mt*Tbeta*(137797425107+1078533300*lmMst1-
690681600*pow2(lmMst1))+
735*MuSUSY*s2t*(223974673+2515800*lmMst1-
1638000*pow2(lmMst1)))*pow3(Mst1))+
3675*(MuSUSY*s2t*(3044017+27472*lmMst1-48480*pow2(lmMst1))
+4*Mt*Tbeta*(486671+31944*lmMst1-15120*pow2(lmMst1)))*pow2(s2t)*pow4(Mst1))+
Tbeta*(490*(59957863+480000*lmMst1-
230400*pow2(lmMst1))*pow2(MuSUSY)*(1-pow2(Sbeta))*pow3(Mt)+
980*(31025111+290880*lmMst1-
251100*pow2(lmMst1))*pow2(Sbeta)*pow3(s2t)*pow5(Mst1)))))/(pow2(Mst1)*pow2(Sbeta)*pow3(Mgl)*pow4(Msq)))+
(5880*xDmsqst1*pow3(Dmsqst1)*(8*Mt*xDmglst1*pow3(Dmglst1)*(-2*s2t*pow2(Mt)*(-625*Tbeta*pow2(MuSUSY)*(1-pow2(Sbeta))+
(36600*Mst1*MuSUSY-
4*(5459+570*lmMst1)*Tbeta*pow2(Mst1))*pow2(Sbeta))+
pow2(Sbeta)*(-150*Mt*(25*MuSUSY-488*Mst1*Tbeta)*pow2(Mst1)*pow2(s2t)-
4*(5*(557+120*lmMst1)*MuSUSY-
6*(3401+105*lmMst1)*Mst1*Tbeta)*pow3(Mt)+
625*Tbeta*pow3(s2t)*pow4(Mst1)))+
Mt*(40*Dmglst1*pow2(Mgl)*(-2*s2t*pow2(Mt)*(-125*Tbeta*pow2(MuSUSY)*(1-pow2(Sbeta))+
(7050*Mst1*MuSUSY-
4*(167+120*lmMst1)*Tbeta*pow2(Mst1))*pow2(Sbeta))+
pow2(Sbeta)*(-150*Mt*(5*MuSUSY-94*Mst1*Tbeta)*pow2(Mst1)*pow2(s2t)-
4*((557+120*lmMst1)*MuSUSY-
9*(423-20*lmMst1)*Mst1*Tbeta)*pow3(Mt)+
125*Tbeta*pow3(s2t)*pow4(Mst1)))+
8*Mgl*pow2(Dmglst1)*(-50*s2t*pow2(Mt)*(-25*Tbeta*pow2(MuSUSY)*(1-pow2(Sbeta))+
(1437*Mst1*MuSUSY-
10*(61+6*lmMst1)*Tbeta*pow2(Mst1))*pow2(Sbeta))+pow2(Sbeta)*(-150*Mt*(25*MuSUSY-479*Mst1*Tbeta)*pow2(Mst1)*pow2(s2t)-
2*(10*(557+120*lmMst1)*MuSUSY-
3*(13147-90*lmMst1)*Mst1*Tbeta)*pow3(Mt)+
625*Tbeta*pow3(s2t)*pow4(Mst1))))+
5*pow3(Mgl)*(-75*Mst1*pow2(Mt)*pow2(s2t)*(-4*(7+12*shiftst1+
24*(lmMst1*(1-shiftst2)+shiftst2))*Tbeta*pow2(MuSUSY)*(1-pow2(Sbeta))+
(336*Mst1*MuSUSY+
(683+768*shiftst1+240*shiftst2+
96*lmMst1*(3-4*shiftst1+shiftst2))*Tbeta*pow2(Mst1))*pow2(Sbeta))+
400*s2t*(21*Tbeta*pow2(MuSUSY)*(1-pow2(Sbeta))+
(6*Mst1*MuSUSY*(4-18*lmMst1*(1-shiftst2)-9*shiftst2)+
2*(65-6*lmMst1)*Tbeta*pow2(Mst1))*pow2(Sbeta))*pow3(Mt)-
8*(100*(65-6*lmMst1)*MuSUSY+
Mst1*(539-60*lmMst1*(91-90*shiftst2)-
2700*shiftst2)*Tbeta)*pow2(Sbeta)*pow4(Mt)+
pow2(Sbeta)*(-150*Mt*(MuSUSY*(7+144*shiftst1-108*shiftst2+
24*lmMst1*(1-2*shiftst1+shiftst2))-
28*Mst1*Tbeta)*pow3(Mst1)*pow3(s2t)+
900*(5-2*lmMst1)*(shiftst1-shiftst2)*Tbeta*pow4(s2t)*pow5(Mst1)))))/(Mst1*pow2(Sbeta)*pow3(Mgl)*pow6(Msq)))/Tbeta))/(4.7628e7*pow6(Mst2))-((xDmst12*z2*pow3(Dmst12)*(50*xDmsqst1*pow3(Dmsqst1)*(-2*Mt*(MuSUSY*s2t*(2*shiftst1-shiftst2)+
2*Mt*(4*shiftst1-shiftst2)*Tbeta)*pow2(Mst1)*pow2(s2t)*pow2(Sbeta)-4*shiftst2*pow2(Mt)*(-(Tbeta*pow2(MuSUSY)*pow2(s2t)*(1-pow2(Sbeta)))+
(6*Mt*MuSUSY*s2t-6*Tbeta*pow2(Mt))*pow2(Sbeta))+
(shiftst1-shiftst2)*Tbeta*pow2(Sbeta)*pow4(Mst1)*pow4(s2t))
+50*(pow2(Dmsqst1)*pow2(Msq)+Dmsqst1*pow4(Msq))*(-2*Mt*(MuSUSY*s2t*(2*shiftst1-shiftst2)+
2*Mt*(4*shiftst1-shiftst2)*Tbeta)*pow2(Mst1)*pow2(s2t)*pow2(Sbeta)-4*shiftst2*pow2(Mt)*(-(Tbeta*pow2(MuSUSY)*pow2(s2t)*(1-pow2(Sbeta)))+
(6*Mt*MuSUSY*s2t-6*Tbeta*pow2(Mt))*pow2(Sbeta))+
(shiftst1-shiftst2)*Tbeta*pow2(Sbeta)*pow4(Mst1)*pow4(s2t))
+(-5*Mt*(MuSUSY*s2t*(40*shiftst1-20*shiftst2+7*shiftst3)+
2*Mt*(80*shiftst1-20*shiftst2+29*shiftst3)*Tbeta)*pow2(Mst1)*pow2(s2t)*pow2(Sbeta)+
2*pow2(Mt)*((15*shiftst1+10*shiftst2+9*shiftst3)*Tbeta*pow2(MuSUSY)*pow2(s2t)*(1-pow2(Sbeta))+
60*shiftst3*(Mt*MuSUSY*s2t+Tbeta*pow2(Mt))*pow2(Sbeta))
+5*(10*shiftst1-10*shiftst2+shiftst3)*Tbeta*pow2(Sbeta)*pow4(Mst1)*pow4(s2t))*pow6(Msq)))/45.-
(Mt*z3*(pow3(Mgl)*(-96*Dmst12*pow2(Mt)*pow4(Mst2)*(210*(11*MuSUSY*s2t-12*Mt*Tbeta)*xDmsqst1*pow2(Mst1)*pow2(Sbeta)*pow3(Dmsqst1)+
pow2(Mst1)*pow2(Sbeta)*(315*(5*MuSUSY*s2t-6*Mt*Tbeta)*pow2(Dmsqst1)*pow2(Msq)+
420*Dmsqst1*(2*MuSUSY*s2t-3*Mt*Tbeta)*pow4(Msq))-
(-60*(22*MuSUSY*s2t-171*Mt*Tbeta)*pow2(Mst1)*pow2(Sbeta)+
4*Mst1*MuSUSY*(-251*MuSUSY*s2t*Tbeta*(1-pow2(Sbeta))+
2069*Mt*pow2(Sbeta))+
Tbeta*(-2125*Mt*pow2(MuSUSY)*(1-pow2(Sbeta))+
1696*s2t*pow2(Sbeta)*pow3(Mst1)))*pow6(Msq))-
Mt*pow2(Dmst12)*pow2(Mst2)*(-10080*xDmsqst1*pow2(Mst1)*((7*Mt*MuSUSY*s2t-6*Tbeta*pow2(Mt))*pow2(Sbeta)-
Tbeta*pow2(s2t)*(-(pow2(MuSUSY)*(1-pow2(Sbeta)))+
11*pow2(Mst1)*pow2(Sbeta)))*pow3(Dmsqst1)+
pow2(Mst1)*(5040*Tbeta*pow2(Dmsqst1)*pow2(Msq)*pow2(s2t)*(-2*pow2(MuSUSY)*(1-pow2(Sbeta))+
15*pow2(Mst1)*pow2(Sbeta))+
10080*Dmsqst1*((7*Mt*MuSUSY*s2t-6*Tbeta*pow2(Mt))*pow2(Sbeta)+
Tbeta*pow2(s2t)*(-(pow2(MuSUSY)*(1-pow2(Sbeta)))+
4*pow2(Mst1)*pow2(Sbeta)))*pow4(Msq))+
(-16*Mst1*Mt*MuSUSY*(-17285*MuSUSY*s2t*Tbeta*(1-pow2(Sbeta))+
51181*Mt*pow2(Sbeta))+
144*pow2(Mst1)*(143*Tbeta*pow2(MuSUSY)*pow2(s2t)*(1-pow2(Sbeta))+
(806*Mt*MuSUSY*s2t-614*Tbeta*pow2(Mt))*pow2(Sbeta))+
192*s2t*(753*MuSUSY*s2t-3290*Mt*Tbeta)*pow2(Sbeta)*pow3(Mst1)+
Tbeta*(95889*pow2(Mt)*pow2(MuSUSY)*(1-pow2(Sbeta))+
63360*pow2(s2t)*pow2(Sbeta)*pow4(Mst1)))*pow6(Msq))
+2*xDmst12*pow3(Dmst12)*(2520*xDmsqst1*pow2(Mst1)*pow3(Dmsqst1)*(Tbeta*(Mt*pow2(s2t)*(-4*pow2(MuSUSY)*(1-pow2(Sbeta))+
29*pow2(Mst1)*pow2(Sbeta))-
24*pow2(Sbeta)*pow3(Mt))+
MuSUSY*pow2(Sbeta)*(16*s2t*pow2(Mt)+2*pow2(Mst1)*pow3(s2t)))+
2520*pow2(Mst1)*(pow2(Dmsqst1)*pow2(Msq)*(Tbeta*(Mt*pow2(s2t)*(-4*pow2(MuSUSY)*(1-pow2(Sbeta))+
15*pow2(Mst1)*pow2(Sbeta))-
36*pow2(Sbeta)*pow3(Mt))+
MuSUSY*pow2(Sbeta)*(30*s2t*pow2(Mt)+2*pow2(Mst1)*pow3(s2t)))+
Dmsqst1*(Tbeta*(Mt*pow2(s2t)*(-4*pow2(MuSUSY)*(1-pow2(Sbeta))+
pow2(Mst1)*pow2(Sbeta))-
48*pow2(Sbeta)*pow3(Mt))+
MuSUSY*pow2(Sbeta)*(44*s2t*pow2(Mt)+2*pow2(Mst1)*pow3(s2t)))*pow4(Msq))-
(16*Mst1*MuSUSY*pow2(Mt)*(-20297*MuSUSY*s2t*Tbeta*(1-pow2(Sbeta))+
76009*Mt*pow2(Sbeta))-
Mt*pow2(Mst1)*(-37669*Tbeta*pow2(MuSUSY)*pow2(s2t)*(1-pow2(Sbeta))+
(131132*Mt*MuSUSY*s2t-135732*Tbeta*pow2(Mt))*pow2(Sbeta))+
pow2(Sbeta)*(4*Mt*s2t*(33783*MuSUSY*s2t-23402*Mt*Tbeta)*pow3(Mst1)+
72*(143*MuSUSY*s2t-37*Mt*Tbeta)*pow2(s2t)*pow4(Mst1))+
Tbeta*(-197889*pow2(MuSUSY)*(1-pow2(Sbeta))*pow3(Mt)+24096*pow2(Sbeta)*pow3(s2t)*pow5(Mst1)))*pow6(Msq))+384*pow3(Mt)*(315*Tbeta*pow2(Dmsqst1)*pow2(Msq)*pow2(Mst1)*pow2(Sbeta)+
630*Tbeta*xDmsqst1*pow2(Mst1)*pow2(Sbeta)*pow3(Dmsqst1)+
4*(-47*Tbeta*pow2(MuSUSY)*(1-pow2(Sbeta))+
(106*Mst1*MuSUSY-
2*(107-12*lmMst1)*Tbeta*pow2(Mst1))*pow2(Sbeta))*pow6(Msq))*pow6(Mst2))+
8*pow6(Msq)*(Dmglst1*pow2(Mgl)*(-(Mt*pow2(Dmst12)*pow2(Mst2)*(-2*Mst1*Mt*MuSUSY*(-105405*MuSUSY*s2t*Tbeta*(1-pow2(Sbeta))+
347332*Mt*pow2(Sbeta))-
pow2(Mst1)*(-50485*Tbeta*pow2(MuSUSY)*pow2(s2t)*(1-pow2(Sbeta))+
(91704*Mt*MuSUSY*s2t-125460*Tbeta*pow2(Mt))*pow2(Sbeta))+
12*s2t*(18791*MuSUSY*s2t-44200*Mt*Tbeta)*pow2(Sbeta)*pow3(Mst1)+
Tbeta*(-147834*pow2(Mt)*pow2(MuSUSY)*(1-pow2(Sbeta))+
34776*pow2(s2t)*pow2(Sbeta)*pow4(Mst1))))+
2*Dmst12*pow2(Mt)*(-34776*(MuSUSY*s2t-8*Mt*Tbeta)*pow2(Mst1)*pow2(Sbeta)+
4*Mst1*MuSUSY*(-18791*MuSUSY*s2t*Tbeta*(1-pow2(Sbeta))+
113412*Mt*pow2(Sbeta))+
Tbeta*(-73917*Mt*pow2(MuSUSY)*(1-pow2(Sbeta))+
188448*s2t*pow2(Sbeta)*pow3(Mst1)))*pow4(Mst2)-
xDmst12*pow3(Dmst12)*(4*Mst1*MuSUSY*pow2(Mt)*(-142987*MuSUSY*s2t*Tbeta*(1-pow2(Sbeta))+
574156*Mt*pow2(Sbeta))+
2*Mt*pow2(Mst1)*(-50485*Tbeta*pow2(MuSUSY)*pow2(s2t)*(1-pow2(Sbeta))+
(56928*Mt*MuSUSY*s2t+152748*Tbeta*pow2(Mt))*pow2(Sbeta))+
pow2(Sbeta)*(Mt*s2t*(90723*MuSUSY*s2t-164264*Mt*Tbeta)*pow3(Mst1)+
(50485*MuSUSY*s2t-80628*Mt*Tbeta)*pow2(s2t)*pow4(Mst1))+
Tbeta*(147834*pow2(MuSUSY)*(1-pow2(Sbeta))*pow3(Mt)+75164*pow2(Sbeta)*pow3(s2t)*pow5(Mst1)))
+48*(-1749*Tbeta*pow2(MuSUSY)*(1-pow2(Sbeta))+
(7852*Mst1*MuSUSY+10068*Tbeta*pow2(Mst1))*pow2(Sbeta))*pow3(Mt)*pow6(Mst2))+Mgl*pow2(Dmglst1)*(Mt*pow2(Dmst12)*pow2(Mst2)*(-2*Mst1*Mt*MuSUSY*(-105405*MuSUSY*s2t*Tbeta*(1-pow2(Sbeta))+
1147928*Mt*pow2(Sbeta))+
3*pow2(Mst1)*(-66705*Tbeta*pow2(MuSUSY)*pow2(s2t)*(1-pow2(Sbeta))+
(259984*Mt*MuSUSY*s2t-314202*Tbeta*pow2(Mt))*pow2(Sbeta))-
78*s2t*(10999*MuSUSY*s2t-26708*Mt*Tbeta)*pow2(Sbeta)*pow3(Mst1)+
Tbeta*(515970*pow2(Mt)*pow2(MuSUSY)*(1-pow2(Sbeta))+
39540*pow2(s2t)*pow2(Sbeta)*pow4(Mst1)))+
2*Dmst12*pow2(Mt)*(60*(659*MuSUSY*s2t+11436*Mt*Tbeta)*pow2(Mst1)*pow2(Sbeta)+
26*Mst1*MuSUSY*(-10999*MuSUSY*s2t*Tbeta*(1-pow2(Sbeta))+
74958*Mt*pow2(Sbeta))+
Tbeta*(-257985*Mt*pow2(MuSUSY)*(1-pow2(Sbeta))+
907296*s2t*pow2(Sbeta)*pow3(Mst1)))*pow4(Mst2)-
xDmst12*pow3(Dmst12)*(-8*Mst1*MuSUSY*pow2(Mt)*(18791*MuSUSY*s2t*Tbeta*(1-pow2(Sbeta))+
86737*Mt*pow2(Sbeta))+
6*Mt*pow2(Mst1)*(-66705*Tbeta*pow2(MuSUSY)*pow2(s2t)*(1-pow2(Sbeta))+
(273164*Mt*MuSUSY*s2t-85482*Tbeta*pow2(Mt))*pow2(Sbeta))+
pow2(Sbeta)*(-(Mt*s2t*(1174137*MuSUSY*s2t-4379080*Mt*Tbeta)*pow3(Mst1))+
3*(66705*MuSUSY*s2t-116812*Mt*Tbeta)*pow2(s2t)*pow4(Mst1))+
Tbeta*(515970*pow2(MuSUSY)*(1-pow2(Sbeta))*pow3(Mt)+285974*pow2(Sbeta)*pow3(s2t)*pow5(Mst1)))+12*(-21141*Tbeta*pow2(MuSUSY)*(1-pow2(Sbeta))+
(151216*Mst1*MuSUSY+112872*Tbeta*pow2(Mst1))*pow2(Sbeta))*pow3(Mt)*pow6(Mst2))+
xDmglst1*pow3(Dmglst1)*(Mt*pow2(Dmst12)*pow2(Mst2)*(-2*Mst1*Mt*MuSUSY*(-737647*MuSUSY*s2t*Tbeta*(1-pow2(Sbeta))+
5205992*Mt*pow2(Sbeta))+
pow2(Mst1)*(-349745*Tbeta*pow2(MuSUSY)*pow2(s2t)*(1-pow2(Sbeta))+
(550536*Mt*MuSUSY*s2t-1849032*Tbeta*pow2(Mt))*pow2(Sbeta))-
24*s2t*(114777*MuSUSY*s2t-258833*Mt*Tbeta)*pow2(Sbeta)*pow3(Mst1)+
Tbeta*(884106*pow2(Mt)*pow2(MuSUSY)*(1-pow2(Sbeta))+
572688*pow2(s2t)*pow2(Sbeta)*pow4(Mst1)))+
6*Dmst12*pow2(Mt)*(48*(3977*MuSUSY*s2t+7908*Mt*Tbeta)*pow2(Mst1)*pow2(Sbeta)+
4*Mst1*MuSUSY*(-76518*MuSUSY*s2t*Tbeta*(1-pow2(Sbeta))+
500581*Mt*pow2(Sbeta))+
Tbeta*(-147351*Mt*pow2(MuSUSY)*(1-pow2(Sbeta))+
966992*s2t*pow2(Sbeta)*pow3(Mst1)))*pow4(Mst2)-
xDmst12*pow3(Dmst12)*(-4*Mst1*MuSUSY*pow2(Mt)*(-278539*MuSUSY*s2t*Tbeta*(1-pow2(Sbeta))+
2202506*Mt*pow2(Sbeta))+
2*Mt*pow2(Mst1)*(-349745*Tbeta*pow2(MuSUSY)*pow2(s2t)*(1-pow2(Sbeta))+
(1123224*Mt*MuSUSY*s2t-710280*Tbeta*pow2(Mt))*pow2(Sbeta))+
pow2(Sbeta)*(-(Mt*s2t*(4967589*MuSUSY*s2t-16623976*Mt*Tbeta)*pow3(Mst1))+
5*(69949*MuSUSY*s2t+59484*Mt*Tbeta)*pow2(s2t)*pow4(Mst1))+
Tbeta*(884106*pow2(MuSUSY)*(1-pow2(Sbeta))*pow3(Mt)+
918216*pow2(Sbeta)*pow3(s2t)*pow5(Mst1)))+
24*(-28665*Tbeta*pow2(MuSUSY)*(1-pow2(Sbeta))+
(241748*Mst1*MuSUSY+113864*Tbeta*pow2(Mst1))*pow2(Sbeta))*pow3(Mt)*pow6(Mst2)))))/(3456.*pow2(Mst1)*pow3(Mgl)))/(Tbeta*pow2(Sbeta)*pow6(Msq)*pow6(Mst2));
}

/// calc coefficient O(as^2,log(mu^2/mt^2)^1) s2
double himalaya::H3::calc_coeff_as_2_log_1_s2(){
   return (pow2(Mt)*(pow2(Mst2)*(Mgl*(4200*pow2(Dmst12)*pow2(Dmglst1+3*Mgl)*pow2(Mst1)*pow2(s2t)+
(420*Dmst12*Mt*s2t*(-(pow2(Dmglst1)*(-2*pow2(Mst2)*(Mst1*Tbeta*(55*pow2(Dmsqst1)+80*Dmsqst1*pow2(Msq))+
(10*MuSUSY-(539-486*lmMst1)*Mst1*Tbeta)*pow4(Msq))+
Dmst12*(Mst1*Tbeta*(30*pow2(Dmsqst1)+80*Dmsqst1*pow2(Msq))+
(56*MuSUSY-3*(81-146*lmMst1)*Mst1*Tbeta)*pow4(Msq))))+
5*pow2(Mgl)*(Dmst12*(Mst1*Tbeta*(20*pow2(Dmsqst1)+40*Dmsqst1*pow2(Msq))-
(12*MuSUSY-(81-26*lmMst1)*Mst1*Tbeta)*pow4(Msq))-
4*pow2(Mst2)*(Mst1*Tbeta*(10*pow2(Dmsqst1)+15*Dmsqst1*pow2(Msq))-
(9*MuSUSY-22*(1-lmMst1)*Mst1*Tbeta)*pow4(Msq)))+Dmglst1*Mgl*(Dmst12*(Mst1*Tbeta*(30*pow2(Dmsqst1)-20*Dmsqst1*pow2(Msq))-
(80*MuSUSY+(109+314*lmMst1)*Mst1*Tbeta)*pow4(Msq))+
10*pow2(Mst2)*(Mst1*Tbeta*(5*pow2(Dmsqst1)+10*Dmsqst1*pow2(Msq))+
2*(6*MuSUSY+Mst1*(Tbeta+40*lmMst1*Tbeta))*pow4(Msq)))))/(Tbeta*pow4(Msq)))+
((2*xDmglst1*pow3(Dmglst1)*(-(pow4(Msq)*(2*Dmst12*Mt*pow2(Mst2)*(3920*Mt*Tbeta*pow2(MuSUSY)*(1-pow2(Sbeta))+
pow2(Sbeta)*(6*(22759-9285*lmMst1)*Mst1*Mt*MuSUSY+
(5040*MuSUSY*s2t+
(183521-23460*lmMst1)*Mt*Tbeta)*pow2(Mst1)+
70*(3341-1554*lmMst1)*s2t*Tbeta*pow3(Mst1)))
-pow2(Dmst12)*(7840*Tbeta*pow2(Mt)*pow2(MuSUSY)*(1-pow2(Sbeta))+
pow2(Sbeta)*(-2*Mt*(6720*MuSUSY*s2t-
(24499-24180*lmMst1)*Mt*Tbeta)*pow2(Mst1)+
3*(84421-13340*lmMst1)*Mst1*MuSUSY*pow2(Mt)+
Tbeta*(4*(44878-26535*lmMst1)*Mt*s2t*pow3(Mst1)-
5040*pow2(s2t)*pow4(Mst1))))+
28*pow2(Mt)*(120*Tbeta*pow2(MuSUSY)*(1-pow2(Sbeta))+
(5*(3197-1554*lmMst1)*Mst1*MuSUSY+
(21851-10110*lmMst1)*Tbeta*pow2(Mst1))*pow2(Sbeta))*pow4(Mst2)))+
Mst1*Mt*pow2(Sbeta)*(840*Dmsqst1*pow2(Msq)*(-(pow2(Dmst12)*(Mt*(20*MuSUSY-21*Mst1*Tbeta)+
7*s2t*Tbeta*pow2(Mst1)))+
Dmst12*(Mt*(20*MuSUSY-21*Mst1*Tbeta)+
27*s2t*Tbeta*pow2(Mst1))*pow2(Mst2)+
Mt*(27*MuSUSY+28*Mst1*Tbeta)*pow4(Mst2))+
420*pow2(Dmsqst1)*(-(pow2(Dmst12)*(Mt*(40*MuSUSY-42*Mst1*Tbeta)-
11*s2t*Tbeta*pow2(Mst1)))+
Dmst12*(Mt*(40*MuSUSY-42*Mst1*Tbeta)+
29*s2t*Tbeta*pow2(Mst1))*pow2(Mst2)+
Mt*(29*MuSUSY+65*Mst1*Tbeta)*pow4(Mst2)))))/(pow2(Mst1)*pow2(Sbeta)*pow4(Msq))+
(210*Mt*xDmsqst1*pow3(Dmsqst1)*(8*xDmglst1*pow3(Dmglst1)*(-(pow2(Dmst12)*(Mt*(20*MuSUSY-21*Mst1*Tbeta)-
18*s2t*Tbeta*pow2(Mst1)))+
Dmst12*(Mt*(20*MuSUSY-21*Mst1*Tbeta)+
2*s2t*Tbeta*pow2(Mst1))*pow2(Mst2)+
Mt*(2*MuSUSY+37*Mst1*Tbeta)*pow4(Mst2))+
4*Mgl*pow2(Dmglst1)*(-(pow2(Dmst12)*(Mt*(40*MuSUSY+9*Mst1*Tbeta)-
10*s2t*Tbeta*pow2(Mst1)))+
Dmst12*(Mt*(40*MuSUSY+9*Mst1*Tbeta)+
30*s2t*Tbeta*pow2(Mst1))*pow2(Mst2)+
Mt*(30*MuSUSY+79*Mst1*Tbeta)*pow4(Mst2))+
80*Dmglst1*pow2(Mgl)*(-(pow2(Dmst12)*(Mt*(2*MuSUSY+3*Mst1*Tbeta)-
2*s2t*Tbeta*pow2(Mst1)))+
Mt*(Dmst12*(2*MuSUSY+3*Mst1*Tbeta)*pow2(Mst2)+
5*Mst1*Tbeta*pow4(Mst2)))+
5*pow3(Mgl)*(-10*Dmst12*(Mt*(4*MuSUSY+Mst1*Tbeta)+
4*s2t*Tbeta*pow2(Mst1))*pow2(Mst2)+
Mt*((40*MuSUSY+7*Mst1*Tbeta)*pow2(Dmst12)-
20*(2*MuSUSY-11*Mst1*Tbeta)*pow4(Mst2)))))/(Mst1*pow6(Msq)))/Tbeta)+
((Mgl*pow2(Mst2)*pow2(Mt)*(14*Dmglst1*Mgl*(-(pow4(Msq)*(pow2(Dmst12)*(400*Tbeta*pow2(MuSUSY)*(1-pow2(Sbeta))+
(15*(11+292*lmMst1)*Mst1*MuSUSY-
3*(117-2960*lmMst1)*Tbeta*pow2(Mst1))*pow2(Sbeta))+
10*Dmst12*pow2(Mst2)*(-40*Tbeta*pow2(MuSUSY)*(1-pow2(Sbeta))+
3*((89-486*lmMst1)*Mst1*MuSUSY+
(331-574*lmMst1)*Tbeta*pow2(Mst1))*pow2(Sbeta))+
600*(-2*Tbeta*pow2(MuSUSY)*(1-pow2(Sbeta))+
((11-40*lmMst1)*Mst1*MuSUSY-
2*(71+4*lmMst1)*Tbeta*pow2(Mst1))*pow2(Sbeta))*pow4(Mst2)))+Mst1*pow2(Sbeta)*(600*Dmsqst1*pow2(Msq)*(-((4*MuSUSY+6*Mst1*Tbeta)*pow2(Dmst12))+
2*Dmst12*(2*MuSUSY+3*Mst1*Tbeta)*pow2(Mst2)+
5*(MuSUSY+5*Mst1*Tbeta)*pow4(Mst2))-
300*pow2(Dmsqst1)*((2*MuSUSY+3*Mst1*Tbeta)*(4*pow2(Dmst12)-4*Dmst12*pow2(Mst2))-
5*(MuSUSY+7*Mst1*Tbeta)*pow4(Mst2))))-
pow2(Dmglst1)*(pow4(Msq)*(-5*pow2(Dmst12)*(1008*Tbeta*pow2(MuSUSY)*(1-pow2(Sbeta))+
(42*(1337-388*lmMst1)*Mst1*MuSUSY+
(46889-21960*lmMst1)*Tbeta*pow2(Mst1))*pow2(Sbeta))+
42*Dmst12*pow2(Mst2)*(120*Tbeta*pow2(MuSUSY)*(1-pow2(Sbeta))+
(10*(763-534*lmMst1)*Mst1*MuSUSY+
(14751-3970*lmMst1)*Tbeta*pow2(Mst1))*pow2(Sbeta))+
280*(-10*Tbeta*pow2(MuSUSY)*(1-pow2(Sbeta))+
(3*(539-486*lmMst1)*Mst1*MuSUSY-
3*(49+614*lmMst1)*Tbeta*pow2(Mst1))*pow2(Sbeta))*pow4(Mst2))+
840*Mst1*pow2(Sbeta)*(Dmsqst1*pow2(Msq)*((40*MuSUSY+9*Mst1*Tbeta)*(pow2(Dmst12)-Dmst12*pow2(Mst2))-
5*(16*MuSUSY+29*Mst1*Tbeta)*pow4(Mst2))+
pow2(Dmsqst1)*((40*MuSUSY+9*Mst1*Tbeta)*(pow2(Dmst12)-Dmst12*pow2(Mst2))-
(55*MuSUSY+112*Mst1*Tbeta)*pow4(Mst2))))-
35*pow2(Mgl)*(pow4(Msq)*(-(pow2(Dmst12)*(-160*Tbeta*pow2(MuSUSY)*(1-pow2(Sbeta))+
(4*(623-222*lmMst1)*Mst1*MuSUSY+
(4087-4518*lmMst1)*Tbeta*pow2(Mst1))*pow2(Sbeta)))+
60*Dmst12*pow2(Mst2)*(-8*Tbeta*pow2(MuSUSY)*(1-pow2(Sbeta))+
((55-62*lmMst1)*Mst1*MuSUSY+
(49-176*lmMst1)*Tbeta*pow2(Mst1))*pow2(Sbeta))
+120*(4*(29-11*lmMst1)*Mst1*MuSUSY*pow2(Sbeta)+
Tbeta*(-6*pow2(MuSUSY)*(1-pow2(Sbeta))-
(185-281*lmMst1-36*z3+129*pow2(lmMst1))*pow2(Mst1)*pow2(Sbeta)))*pow4(Mst2))+
30*Mst1*pow2(Sbeta)*(pow2(Dmsqst1)*(-4*(10*MuSUSY+3*Mst1*Tbeta)*pow2(Dmst12)+
5*Dmst12*(8*MuSUSY+3*Mst1*Tbeta)*pow2(Mst2)+
10*(8*MuSUSY-(29-18*lmMst1)*Mst1*Tbeta)*pow4(Mst2))+
Dmsqst1*pow2(Msq)*(-((40*MuSUSY+17*Mst1*Tbeta)*pow2(Dmst12))+
20*Dmst12*(2*MuSUSY+Mst1*Tbeta)*pow2(Mst2)+
120*(MuSUSY-3*(1-lmMst1)*Mst1*Tbeta)*pow4(Mst2))))))/pow4(Msq)+(xDmst12*pow3(Dmst12)*(14*Dmglst1*pow2(Mgl)*(Mst1*Mt*pow2(Sbeta)*(-300*pow2(Dmsqst1)*pow2(Msq)*(-4*Mt*(2*MuSUSY+3*Mst1*Tbeta)+
11*s2t*Tbeta*pow2(Mst1))+
1200*xDmsqst1*(Mt*(2*MuSUSY+3*Mst1*Tbeta)-
4*s2t*Tbeta*pow2(Mst1))*pow3(Dmsqst1)+
600*Dmsqst1*(Mt*(4*MuSUSY+6*Mst1*Tbeta)-
3*s2t*Tbeta*pow2(Mst1))*pow4(Msq))+
(400*Tbeta*pow2(Mt)*pow2(MuSUSY)*(1-pow2(Sbeta))+
pow2(Sbeta)*(12*Mt*(100*MuSUSY*s2t+
(769+45*lmMst1)*Mt*Tbeta)*pow2(Mst1)+
60*(50-97*lmMst1)*Mst1*MuSUSY*pow2(Mt)+
Tbeta*(45*(53+112*lmMst1)*Mt*s2t*pow3(Mst1)-
3000*pow2(s2t)*pow4(Mst1))))*pow6(Msq))-
7*pow3(Mgl)*(Mst1*pow2(Sbeta)*(150*(40*MuSUSY+9*Mst1*Tbeta)*pow2(Dmsqst1)*pow2(Msq)*pow2(Mt)+
Mt*(-600*xDmsqst1*(-(Mt*(10*MuSUSY+Mst1*Tbeta))+
10*s2t*Tbeta*pow2(Mst1))*pow3(Dmsqst1)+
300*Dmsqst1*(Mt*(20*MuSUSY+7*Mst1*Tbeta)+
20*s2t*Tbeta*pow2(Mst1))*pow4(Msq)))+
(800*Tbeta*pow2(Mt)*pow2(MuSUSY)*(1-pow2(Sbeta))+
pow2(Sbeta)*(-(Mt*(1200*MuSUSY*s2t-
(17933-13380*lmMst1)*Mt*Tbeta)*pow2(Mst1))+
20*(421+486*lmMst1)*Mst1*MuSUSY*pow2(Mt)+
Tbeta*(1120*(17-3*lmMst1)*Mt*s2t*pow3(Mst1)+
7200*pow2(s2t)*pow4(Mst1))))*pow6(Msq))+
2*(xDmglst1*pow3(Dmglst1)*(Mst1*Mt*pow2(Sbeta)*(420*pow2(Dmsqst1)*pow2(Msq)*(Mt*(40*MuSUSY-42*Mst1*Tbeta)-
51*s2t*Tbeta*pow2(Mst1))+
840*(xDmsqst1*(Mt*(20*MuSUSY-21*Mst1*Tbeta)-
38*s2t*Tbeta*pow2(Mst1))*pow3(Dmsqst1)+
Dmsqst1*(Mt*(20*MuSUSY-21*Mst1*Tbeta)-
13*s2t*Tbeta*pow2(Mst1))*pow4(Msq)))-
(7840*Tbeta*pow2(Mt)*pow2(MuSUSY)*(1-pow2(Sbeta))+
pow2(Sbeta)*(-6*Mt*(6160*MuSUSY*s2t+
(44841+8300*lmMst1)*Mt*Tbeta)*pow2(Mst1)+
6*(38903+5230*lmMst1)*Mst1*MuSUSY*pow2(Mt)+
Tbeta*(-((83831+66120*lmMst1)*Mt*s2t*pow3(Mst1))+
1680*pow2(s2t)*pow4(Mst1))))*pow6(Msq))+
Mgl*pow2(Dmglst1)*(420*Mst1*pow2(Sbeta)*(Mt*(pow2(Dmsqst1)*pow2(Msq)*(Mt*(40*MuSUSY+9*Mst1*Tbeta)-
25*s2t*Tbeta*pow2(Mst1))+
xDmsqst1*(Mt*(40*MuSUSY+9*Mst1*Tbeta)-
50*s2t*Tbeta*pow2(Mst1))*pow3(Dmsqst1))+
Dmsqst1*(40*MuSUSY+9*Mst1*Tbeta)*pow2(Mt)*pow4(Msq))-
(2520*Tbeta*pow2(Mt)*pow2(MuSUSY)*(1-pow2(Sbeta))+
pow2(Sbeta)*(-2*Mt*(9660*MuSUSY*s2t+
(37663+13215*lmMst1)*Mt*Tbeta)*pow2(Mst1)+
420*(287+73*lmMst1)*Mst1*MuSUSY*pow2(Mt)+
Tbeta*(-105*(947+488*lmMst1)*Mt*s2t*pow3(Mst1)+
7980*pow2(s2t)*pow4(Mst1))))*pow6(Msq)))))/pow6(Msq))/(Tbeta*pow2(Mst1)*pow2(Sbeta))))/(14175.*pow3(Mgl)*pow6(Mst2));
}

/// calc coefficient O(as^2,log(mu^2/mt^2)^2) s2
double himalaya::H3::calc_coeff_as_2_log_2_s2(){
   return (-2*pow3(Mt)*(pow3(Mgl)*(pow4(Msq)*(-5*pow2(Dmst12)*(Mt*(16*MuSUSY+43*Mst1*Tbeta)+
64*s2t*Tbeta*pow2(Mst1))*pow2(Mst2)+
2*Mst1*(77*Mt+120*Mst1*s2t)*Tbeta*xDmst12*pow3(Dmst12)+
10*Dmst12*(Mt*(16*MuSUSY+35*Mst1*Tbeta)+
48*s2t*Tbeta*pow2(Mst1))*pow4(Mst2))-
30*Mt*(Mst1*Tbeta*(5*pow2(Dmsqst1)+10*Dmsqst1*pow2(Msq))-
4*(4*MuSUSY-(9-14*lmMst1)*Mst1*Tbeta)*pow4(Msq))*pow6(Mst2))
+pow4(Msq)*(-8*Mgl*pow2(Dmglst1)*(2*pow2(Dmst12)*(-2*Mt*(MuSUSY+2*Mst1*Tbeta)+
3*s2t*Tbeta*pow2(Mst1))*pow2(Mst2)-
2*xDmst12*(-(Mt*(MuSUSY+Mst1*Tbeta))+
5*s2t*Tbeta*pow2(Mst1))*pow3(Dmst12)+
Mt*(2*Dmst12*(3*MuSUSY+7*Mst1*Tbeta)*pow4(Mst2)+
125*Mst1*Tbeta*pow6(Mst2)))-
16*(Dmglst1*pow2(Mgl)*(pow2(Dmst12)*(Mt*(MuSUSY-4*Mst1*Tbeta)+
10*s2t*Tbeta*pow2(Mst1))*pow2(Mst2)-
xDmst12*(2*Mt*(MuSUSY+Mst1*Tbeta)+
9*s2t*Tbeta*pow2(Mst1))*pow3(Dmst12)+
10*Dmst12*Mst1*(Mt-Mst1*s2t)*Tbeta*pow4(Mst2)-
5*Mt*(2*MuSUSY-17*Mst1*Tbeta)*pow6(Mst2))+
xDmglst1*pow3(Dmglst1)*(-(pow2(Dmst12)*(Mt*(2*MuSUSY+4*Mst1*Tbeta)+
s2t*Tbeta*pow2(Mst1))*pow2(Mst2))+
xDmst12*(Mt*(MuSUSY+4*Mst1*Tbeta)-s2t*Tbeta*pow2(Mst1))*pow3(Dmst12)+
Dmst12*(Mt*(3*MuSUSY+4*Mst1*Tbeta)+
4*s2t*Tbeta*pow2(Mst1))*pow4(Mst2)+
Mt*(4*MuSUSY+51*Mst1*Tbeta)*pow6(Mst2))))))/(45.*Mst1*Tbeta*pow3(Mgl)*pow4(Msq)*pow6(Mst2));
}

/// calc coefficient O(as^2,log(mu^2/mt^2)^3) s2
double himalaya::H3::calc_coeff_as_2_log_3_s2(){
   return (184*pow4(Mt))/3.;
}

/// calc coefficient (as^0,log(mu^2/mt^2)^0) s12
double himalaya::H3::calc_coeff_as_0_log_0_s12(){
   return -(Dmst12*Mt*MuSUSY*s2t*(-2*Dmst12*Mt*(MuSUSY*s2t+6*Mt*Tbeta)*pow2(Mst2)+
pow2(Dmst12)*(2*Mt*MuSUSY*s2t+
Tbeta*(8*pow2(Mt)-pow2(Mst1)*pow2(s2t)))+
24*Tbeta*pow2(Mt)*pow4(Mst2)))/(96.*Tbeta*pow6(Mst2));
}

/// calc coeffcieint O(as^1,log(mu^2/mt^2)^0) s12
double himalaya::H3::calc_coeff_as_1_log_0_s12(){
   return (Mt*MuSUSY*(25*pow3(Mgl)*(-2*Mt*pow2(Dmst12)*pow2(Mst2)*(4*Mt*s2t*(MuSUSY+6*(1-3*lmMst1)*Mst1*Tbeta)+
2*(5+6*lmMst1)*Tbeta*pow2(Mt)-
3*Mst1*(MuSUSY-3*(1+2*lmMst1)*Mst1*Tbeta)*pow2(s2t))-
pow3(Dmst12)*(-8*s2t*((5+6*lmMst1)*MuSUSY+12*(1-lmMst1)*Mst1*Tbeta)*pow2(Mt)
+2*Mst1*Mt*((1+6*lmMst1)*MuSUSY-6*(1+3*lmMst1)*Mst1*Tbeta)*pow2(s2t)+Tbeta*(8*(5+6*lmMst1)*pow3(Mt)-3*pow3(Mst1)*pow3(s2t)))-
8*Dmst12*(3*MuSUSY*(s2t+2*lmMst1*s2t)-
(2*(5+6*lmMst1)*Mt-9*(1+2*lmMst1)*Mst1*s2t)*Tbeta)*pow2(Mt)*pow4(Mst2)+144*(1+lmMst1)*Tbeta*pow3(Mt)*pow6(Mst2))
-2*(xDmglst1*pow3(Dmglst1)*(-(Mt*pow2(Dmst12)*pow2(Mst2)*(-80*Mt*s2t*(13*MuSUSY-18*Mst1*Tbeta)+
16*(238-45*lmMst1)*Tbeta*pow2(Mt)+
5*Mst1*(20*(13-12*lmMst1)*MuSUSY+
3*(77-30*lmMst1)*Mst1*Tbeta)*pow2(s2t)))+
pow3(Dmst12)*(-60*s2t*((9+10*lmMst1)*MuSUSY+2*(47-60*lmMst1)*Mst1*Tbeta)*pow2(Mt)+5*Mst1*Mt*(40*(13-12*lmMst1)*MuSUSY+
9*(43-10*lmMst1)*Mst1*Tbeta)*pow2(s2t)+
Tbeta*(2*(1577+270*lmMst1)*pow3(Mt)-
50*(13-12*lmMst1)*pow3(Mst1)*pow3(s2t)))-
2*Dmst12*(10*(77-30*lmMst1)*MuSUSY*s2t-
((2231-990*lmMst1)*Mt+60*(71-60*lmMst1)*Mst1*s2t)*Tbeta)*pow2(Mt)*pow4(Mst2)+
8*(977-480*lmMst1)*Tbeta*pow3(Mt)*pow6(Mst2))+
Dmglst1*pow2(Mgl)*(-(Mt*pow2(Dmst12)*pow2(Mst2)*(100*Mst1*s2t*((2-6*lmMst1)*MuSUSY*s2t+3*Mt*Tbeta)+
Tbeta*(6*(41-90*lmMst1)*pow2(Mt)+
75*(5-6*lmMst1)*pow2(Mst1)*pow2(s2t))))-
pow3(Dmst12)*(-100*s2t*((5-6*lmMst1)*MuSUSY+36*lmMst1*Mst1*Tbeta)*pow2(Mt)-
25*Mst1*Mt*(16*(1-3*lmMst1)*MuSUSY+
3*(5-6*lmMst1)*Mst1*Tbeta)*pow2(s2t)+
Tbeta*(24*(17-30*lmMst1)*pow3(Mt)+
100*(1-3*lmMst1)*pow3(Mst1)*pow3(s2t)))-
100*Dmst12*((5-6*lmMst1)*MuSUSY*s2t-
3*((3-6*lmMst1)*Mt+2*(1-6*lmMst1)*Mst1*s2t)*Tbeta)*pow2(Mt)*pow4(Mst2)+
200*(7-15*lmMst1)*Tbeta*pow3(Mt)*pow6(Mst2))+
Mgl*pow2(Dmglst1)*(-3*Mt*pow2(Dmst12)*pow2(Mst2)*(-10*Mt*s2t*(20*MuSUSY-7*Mst1*Tbeta)+
16*(46-15*lmMst1)*Tbeta*pow2(Mt)+
25*Mst1*(2*(5-6*lmMst1)*MuSUSY+
(11-6*lmMst1)*Mst1*Tbeta)*pow2(s2t))+
pow3(Dmst12)*(-20*s2t*(5*(1+6*lmMst1)*MuSUSY+6*(29-45*lmMst1)*Mst1*Tbeta)*pow2(Mt)+75*Mst1*Mt*(4*(5-6*lmMst1)*MuSUSY+(17-6*lmMst1)*Mst1*Tbeta)*pow2(s2t)+Tbeta*(6*(259+90*lmMst1)*pow3(Mt)-
75*(5-6*lmMst1)*pow3(Mst1)*pow3(s2t)))-
2*Dmst12*(50*(11-6*lmMst1)*MuSUSY*s2t-
3*((477-330*lmMst1)*Mt+50*(13-18*lmMst1)*Mst1*s2t)*Tbeta)*pow2(Mt)*pow4(Mst2)+
1200*(4-3*lmMst1)*Tbeta*pow3(Mt)*pow6(Mst2)))))/(2700.*Mst1*Tbeta*pow3(Mgl)*pow6(Mst2));
}

/// calc coefficient O(as^1,log(mu^2/mt^2)^1) s12
double himalaya::H3::calc_coeff_as_1_log_1_s12(){
   return (-2 * MuSUSY * pow4(Mt) * (5 * pow2(Mst2) * pow3(Mgl) * (pow2(Dmst12) - 2 *
         Dmst12 * pow2(Mst2) - 6 * pow4(Mst2)) + Dmst12 * Mgl * pow2(Dmglst1) *
      (pow2(Dmst12) - 2 * Dmst12 * pow2(Mst2) + 3 * pow4(Mst2)) + xDmglst1 *
      pow3(Dmglst1) * (-2 * pow2(Dmst12) * pow2(Mst2) + pow3(Dmst12) + 3 *
         Dmst12 * pow4(Mst2) + 4 * pow6(Mst2)) - Dmglst1 * pow2(Mgl) * (-(
         pow2(Dmst12) * pow2(Mst2)) + 2 * pow3(Dmst12) + 10 * pow6(Mst2)))) / (
      45. * Mst1 * pow3(Mgl) * pow6(Mst2));
}

/// calc coefficient O(as^2,log(mu^2/mt^2)^0) s12
double himalaya::H3::calc_coeff_as_2_log_0_s12(){
   return -(Mt*MuSUSY*(((-490*Mt*MuSUSY*((Dmst12*(-30*pow2(s2t)*(-(pow4(Msq)*(Dmst12*(4*pow2(Dmglst1)*(1732531+16896*lmMst1-
24840*pow2(lmMst1))+
20*Dmglst1*Mgl*(84209+1264*lmMst1-240*pow2(lmMst1)))*(2*Dmst12-pow2(Mst2))-
pow2(Mgl)*(pow2(Dmst12)*(350605+2880*shiftst2+
96*lmMst1*(115-90*shiftst1-60*shiftst2-
54*shiftst3)+8352*shiftst3+
4320*shiftst1*(1-2*z2)+
(-5760*shiftst2-5184*shiftst3)*z2-
2160*pow2(lmMst1))+
40*Dmst12*(1429-
2*lmMst1*(227-180*(shiftst1+shiftst2)-
54*shiftst3)-126*shiftst3-
180*(shiftst2+shiftst1*(1-2*z2))+
(360*shiftst2+108*shiftst3)*z2+
24*pow2(lmMst1))*pow2(Mst2)+
1440*(10*shiftst1-10*shiftst2+shiftst3)*(1-2*(lmMst1+z2))*pow4(Mst2))))+
1200*pow2(Mgl)*(pow2(Dmsqst1)*(2*(7+6*shiftst1+
24*lmMst1*(1-shiftst2)+
shiftst2*(30-24*z2))*pow2(Dmst12)-
Dmst12*(7+30*shiftst1+
12*lmMst1*(2-shiftst1-shiftst2)+
6*shiftst2-12*(shiftst1+shiftst2)*z2)*pow2(Mst2)+
24*(shiftst1-shiftst2)*(2-lmMst1-z2)*pow4(Mst2))+
Dmsqst1*pow2(Msq)*(2*(7+24*lmMst1*(1-shiftst2)+
shiftst2*(36-24*z2))*pow2(Dmst12)-
Dmst12*(7+12*lmMst1*(2-shiftst1-shiftst2)+
(shiftst1+shiftst2)*(18-12*z2))*pow2(Mst2)
+12*(shiftst1-shiftst2)*(3-2*(lmMst1+z2))*pow4(Mst2))))-
(240*Mt*s2t*(pow2(Dmglst1)*(2000*(pow2(Dmsqst1)+Dmsqst1*pow2(Msq))*(pow2(Dmst12)-Dmst12*pow2(Mst2)+
pow4(Mst2))-
pow4(Msq)*(20*pow2(Dmst12)*(32829-1852*lmMst1-660*pow2(lmMst1))+
Dmst12*(958501+24456*lmMst1-11520*pow2(lmMst1))*pow2(Mst2)-
2*(1286791+5936*lmMst1-
18120*pow2(lmMst1))*pow4(Mst2)))+
Dmglst1*Mgl*(2000*(pow2(Dmsqst1)+Dmsqst1*pow2(Msq))*(pow2(Dmst12)-Dmst12*pow2(Mst2)+
pow4(Mst2))-
pow4(Msq)*(2*pow2(Dmst12)*(1282471-7264*lmMst1-18120*pow2(lmMst1))
-Dmst12*(949861-1944*lmMst1-11520*pow2(lmMst1))*pow2(Mst2)-
20*(33261-532*lmMst1-660*pow2(lmMst1))*pow4(Mst2)))+
5*pow2(Mgl)*(1680*(pow2(Dmsqst1)+Dmsqst1*pow2(Msq))*(pow2(Dmst12)-Dmst12*pow2(Mst2)+
pow4(Mst2))-
pow4(Msq)*(2*pow2(Dmst12)*(36863-80*lmMst1-552*pow2(lmMst1))-
3*Dmst12*(10473-40*lmMst1-256*pow2(lmMst1))*pow2(Mst2)-
8*(1361+10*lmMst1+54*pow2(lmMst1))*pow4(Mst2)))))/Mst1))/pow4(Msq)+
(pow2(Mt)*(-16*pow2(Dmglst1)*((3891491+27200*lmMst1-19200*pow2(lmMst1))*(-9*pow2(Dmst12)*pow2(Mst2)+
9*(pow3(Dmst12)+Dmst12*pow4(Mst2)))+
50*(345581+4896*lmMst1-3456*pow2(lmMst1))*pow6(Mst2))-
400*Dmglst1*Mgl*((403559+384*lmMst1-4608*pow2(lmMst1))*(-(pow2(Dmst12)*pow2(Mst2))+pow3(Dmst12)+
Dmst12*pow4(Mst2))+
24*(9631+16*lmMst1-192*pow2(lmMst1))*pow6(Mst2))+
15*pow2(Mgl)*(-(pow2(Dmst12)*(852541+9216*lmMst1+6144*pow2(lmMst1))*pow2(Mst2))+
(3527322-94208*lmMst1-49152*pow2(lmMst1))*pow3(Dmst12)-
160*Dmst12*(11389-704*lmMst1-384*pow2(lmMst1))*pow4(Mst2)-
1920*(349-56*lmMst1-32*pow2(lmMst1))*pow6(Mst2))))/pow2(Mst1)))/pow2(Mgl)+
((-23520*xDmsqst1*pow3(Dmsqst1)*(125*pow3(Mgl)*(-(pow2(Dmst12)*pow2(Mst2)*(-3*s2t*(112*MuSUSY+
3*Mst1*(83+96*shiftst2)*Tbeta)*pow2(Mt)+
6*Mst1*Mt*(42*Mst1*Tbeta-
MuSUSY*(7+12*lmMst1*(2-shiftst1-shiftst2)+
shiftst2*(-6-12*z2)+6*shiftst1*(7-2*z2)))*pow2(s2t)+Tbeta*(16*(65-6*lmMst1)*pow3(Mt)+
54*(shiftst1-shiftst2)*(5-2*(lmMst1+z2))*pow3(Mst1)*pow3(s2t))))
+pow3(Dmst12)*(-48*s2t*(7*MuSUSY+
Mst1*Tbeta*(4-18*lmMst1*(1-shiftst2)+
shiftst2*(-9+18*z2)))*pow2(Mt)+
12*Mst1*Mt*(42*Mst1*Tbeta-
MuSUSY*(7+12*shiftst1+
24*(lmMst1*(1-shiftst2)+shiftst2)-
24*shiftst2*z2))*pow2(s2t)+
Tbeta*(16*(65-6*lmMst1)*pow3(Mt)+
3*(7+
24*lmMst1*(1-2*shiftst1+shiftst2)+
48*shiftst1*(3-z2)+
shiftst2*(-108+24*z2))*pow3(Mst1)*pow3(s2t)))+2*Dmst12*Mt*(-3*Mt*s2t*(56*MuSUSY+
Mst1*Tbeta*(217+144*lmMst1*(1-shiftst2)+
72*shiftst2*(5-2*z2)))+
8*(65-6*lmMst1)*Tbeta*pow2(Mt)-
36*Mst1*MuSUSY*(shiftst1-shiftst2)*(5-2*(lmMst1+z2))*pow2(s2t))*pow4(Mst2)+
32*Tbeta*((31-3*lmMst1)*Mt+
9*Mst1*s2t*(shiftst1-shiftst2)*(7-3*(lmMst1+z2)))*pow2(Mt)*pow6(Mst2))+
Mt*(20*Dmglst1*pow2(Mgl)*(-(pow2(Dmst12)*pow2(Mst2)*(-100*Mt*s2t*(5*MuSUSY-141*Mst1*Tbeta)+
Tbeta*(4*(557+120*lmMst1)*pow2(Mt)+
375*pow2(Mst1)*pow2(s2t))))+
2*(-50*Mt*s2t*(5*MuSUSY-141*Mst1*Tbeta)+
Tbeta*(2*(557+120*lmMst1)*pow2(Mt)+
375*pow2(Mst1)*pow2(s2t)))*pow3(Dmst12)-
4*Dmst12*Mt*(125*MuSUSY*s2t-
((557+120*lmMst1)*Mt+3525*Mst1*s2t)*Tbeta)*pow4(Mst2)+3600*Tbeta*pow2(Mt)*pow6(Mst2))-
4*xDmglst1*pow3(Dmglst1)*(5*pow2(Dmst12)*pow2(Mst2)*(-20*Mt*s2t*(25*MuSUSY-732*Mst1*Tbeta)+
Tbeta*(4*(557+120*lmMst1)*pow2(Mt)+
375*pow2(Mst1)*pow2(s2t)))-
10*(-(Mt*s2t*(250*MuSUSY-7320*Mst1*Tbeta))+
Tbeta*(2*(557+120*lmMst1)*pow2(Mt)+
375*pow2(Mst1)*pow2(s2t)))*pow3(Dmst12)+
20*Dmst12*Mt*(125*MuSUSY*s2t-
((557+120*lmMst1)*Mt+3660*Mst1*s2t)*Tbeta)*pow4(Mst2)+8*(2389-30*lmMst1)*Tbeta*pow2(Mt)*pow6(Mst2))-
20*Mgl*pow2(Dmglst1)*(pow2(Dmst12)*pow2(Mst2)*(-(Mt*s2t*(500*MuSUSY-14370*Mst1*Tbeta))+
Tbeta*(4*(557+120*lmMst1)*pow2(Mt)+
375*pow2(Mst1)*pow2(s2t)))-
2*(-(Mt*s2t*(250*MuSUSY-7185*Mst1*Tbeta))+
Tbeta*(2*(557+120*lmMst1)*pow2(Mt)+
375*pow2(Mst1)*pow2(s2t)))*pow3(Dmst12)+
2*Dmst12*Mt*(250*MuSUSY*s2t-
((1114+240*lmMst1)*Mt+7185*Mst1*s2t)*Tbeta)*pow4(Mst2)+24*(56-15*lmMst1)*Tbeta*pow2(Mt)*pow6(Mst2)))))/(Mst1*pow6(Msq))+
((8*xDmglst1*pow3(Dmglst1)*(10*Mt*pow2(Dmst12)*pow2(Mst2)*(5880*Mst1*(pow2(Dmsqst1)+Dmsqst1*pow2(Msq))*(-20*Mt*s2t*(25*MuSUSY-732*Mst1*Tbeta)+
Tbeta*(4*(557+120*lmMst1)*pow2(Mt)+
375*pow2(Mst1)*pow2(s2t)))-
(98*Mst1*Mt*s2t*(MuSUSY*(99874229+1352280*lmMst1-
633600*pow2(lmMst1))+
30*Mst1*Tbeta*(611423-9984*lmMst1-23040*pow2(lmMst1)))
+(-30*Mst1*Tbeta*(1150678153+9276404*lmMst1-
7140000*pow2(lmMst1))+
98*MuSUSY*(59957863+480000*lmMst1-
230400*pow2(lmMst1)))*pow2(Mt)-
147*(2*Mst1*Tbeta*(31025111+290880*lmMst1-
251100*pow2(lmMst1))+
5*MuSUSY*(3044017+27472*lmMst1-48480*pow2(lmMst1)))*pow2(Mst1)*pow2(s2t))*pow4(Msq))+
pow3(Dmst12)*(-117600*Mst1*Mt*(pow2(Dmsqst1)+Dmsqst1*pow2(Msq))*(-(Mt*s2t*(250*MuSUSY-7320*Mst1*Tbeta))+
Tbeta*(2*(557+120*lmMst1)*pow2(Mt)+
375*pow2(Mst1)*pow2(s2t)))+
pow4(Msq)*(1960*Mst1*s2t*(MuSUSY*(37824007+770520*lmMst1-
131400*pow2(lmMst1))+
60*Mst1*Tbeta*(580211+498*lmMst1-21060*pow2(lmMst1)))*pow2(Mt)-
735*Mt*(Mst1*Tbeta*(223974673+2515800*lmMst1-
1638000*pow2(lmMst1))+
20*MuSUSY*(3044017+27472*lmMst1-48480*pow2(lmMst1)))*pow2(Mst1)*pow2(s2t)+4*(-16*Mst1*Tbeta*(4572123029+49945815*lmMst1-
24119550*pow2(lmMst1))+
245*MuSUSY*(59957863+480000*lmMst1-
230400*pow2(lmMst1)))*pow3(Mt)+
3675*Tbeta*(3044017+27472*lmMst1-48480*pow2(lmMst1))*pow3(s2t)*pow4(Mst1)))+
4*Dmst12*pow2(Mt)*(58800*Mst1*(125*MuSUSY*s2t-
((557+120*lmMst1)*Mt+3660*Mst1*s2t)*Tbeta)*(pow2(Dmsqst1)+Dmsqst1*pow2(Msq))+
(Mt*(-2*Mst1*Tbeta*(49723877243+296163780*lmMst1-
342543600*pow2(lmMst1))+
245*MuSUSY*(59957863+480000*lmMst1-
230400*pow2(lmMst1)))+
980*Mst1*s2t*(MuSUSY*(31025111+290880*lmMst1-
251100*pow2(lmMst1))-
15*Mst1*Tbeta*(548999+10980*lmMst1-19080*pow2(lmMst1))))*pow4(Msq))*pow4(Mst2)+784*pow3(Mt)*(Mst1*Tbeta*(60*(3778-435*lmMst1)*pow2(Dmsqst1)+
360*Dmsqst1*(463-135*lmMst1)*pow2(Msq))+
(-2*Mst1*Tbeta*(122282257+499080*lmMst1-
1351800*pow2(lmMst1))+
15*MuSUSY*(3877891+46400*lmMst1-19200*pow2(lmMst1)))*pow4(Msq))*pow6(Mst2)))/pow4(Msq)-
(55125*z3*(pow3(Mgl)*(48*Dmst12*pow2(Mt)*pow4(Mst2)*(s2t*Tbeta*pow2(Mst1)*(1575*pow2(Dmsqst1)*pow2(Msq)+
2310*xDmsqst1*pow3(Dmsqst1)+
840*Dmsqst1*pow4(Msq))+
2*(MuSUSY*(2125*Mt+1004*Mst1*s2t)+
Tbeta*(-4138*Mst1*Mt+660*s2t*pow2(Mst1)))*pow6(Msq))-
Mt*pow2(Dmst12)*pow2(Mst2)*(pow2(Mst1)*(10080*MuSUSY*pow2(Dmsqst1)*pow2(Msq)*pow2(s2t)+
5040*s2t*((2*MuSUSY*s2t+7*Mt*Tbeta)*xDmsqst1*pow3(Dmsqst1)+
Dmsqst1*(2*MuSUSY*s2t-7*Mt*Tbeta)*pow4(Msq)))-(16*Mst1*Mt*s2t*(17285*MuSUSY+3627*Mst1*Tbeta)+
(95889*MuSUSY-409448*Mst1*Tbeta)*pow2(Mt)+
144*(143*MuSUSY+502*Mst1*Tbeta)*pow2(Mst1)*pow2(s2t))*pow6(Msq))+
2*pow3(Dmst12)*(2520*s2t*pow2(Mst1)*(pow2(Dmsqst1)*pow2(Msq)*(4*Mt*MuSUSY*s2t+
Tbeta*(-15*pow2(Mt)-pow2(Mst1)*pow2(s2t)))
+xDmsqst1*(4*Mt*MuSUSY*s2t+
Tbeta*(-8*pow2(Mt)-pow2(Mst1)*pow2(s2t)))*pow3(Dmsqst1)+
Dmsqst1*(4*Mt*MuSUSY*s2t+
Tbeta*(-22*pow2(Mt)-pow2(Mst1)*pow2(s2t)))*pow4(Msq))+
(-2*Mst1*s2t*(162376*MuSUSY+32783*Mst1*Tbeta)*pow2(Mt)+
Mt*(37669*MuSUSY+67566*Mst1*Tbeta)*pow2(Mst1)*pow2(s2t)-
(197889*MuSUSY-608072*Mst1*Tbeta)*pow3(Mt)+5148*Tbeta*pow3(s2t)*pow4(Mst1))*pow6(Msq))+
1536*(47*MuSUSY-53*Mst1*Tbeta)*pow3(Mt)*pow6(Msq)*pow6(Mst2))+
4*pow6(Msq)*(Mgl*pow2(Dmglst1)*(-2*Mt*pow2(Dmst12)*pow2(Mst2)*(6*Mst1*Mt*s2t*(35135*MuSUSY+64996*Mst1*Tbeta)+
(515970*MuSUSY-1147928*Mst1*Tbeta)*pow2(Mt)-
3*(66705*MuSUSY+142987*Mst1*Tbeta)*pow2(Mst1)*pow2(s2t))+
pow3(Dmst12)*(-8*Mst1*s2t*(37582*MuSUSY-204873*Mst1*Tbeta)*pow2(Mt)-
3*Mt*(266820*MuSUSY+391379*Mst1*Tbeta)*pow2(Mst1)*pow2(s2t)+
28*(36855*MuSUSY-24782*Mst1*Tbeta)*pow3(Mt)+200115*Tbeta*pow3(s2t)*pow4(Mst1))+
4*Dmst12*(39*Mt*(6615*MuSUSY-24986*Mst1*Tbeta)+
2*Mst1*s2t*(142987*MuSUSY-9885*Mst1*Tbeta))*pow2(Mt)*pow4(Mst2)+
24*(21141*MuSUSY-75608*Mst1*Tbeta)*pow3(Mt)*pow6(Mst2))+
xDmglst1*pow3(Dmglst1)*(-2*Mt*pow2(Dmst12)*pow2(Mst2)*(2*Mst1*Mt*s2t*(737647*MuSUSY+137634*Mst1*Tbeta)+
(884106*MuSUSY-5205992*Mst1*Tbeta)*pow2(Mt)-
(349745*MuSUSY+1377324*Mst1*Tbeta)*pow2(Mst1)*pow2(s2t))+
pow3(Dmst12)*(8*Mst1*s2t*(278539*MuSUSY+280806*Mst1*Tbeta)*pow2(Mt)
-11*Mt*(127180*MuSUSY+451599*Mst1*Tbeta)*pow2(Mst1)*pow2(s2t)+
4*(442053*MuSUSY-2202506*Mst1*Tbeta)*pow3(Mt)+349745*Tbeta*pow3(s2t)*pow4(Mst1))
+12*Dmst12*(Mt*(147351*MuSUSY-1001162*Mst1*Tbeta)+
24*Mst1*s2t*(12753*MuSUSY-3977*Mst1*Tbeta))*pow2(Mt)*pow4(Mst2)+
624*(2205*MuSUSY-9298*Mst1*Tbeta)*pow3(Mt)*pow6(Mst2))+
Dmglst1*pow2(Mgl)*(2*Mt*pow2(Dmst12)*pow2(Mst2)*(6*Mst1*Mt*s2t*(35135*MuSUSY-7642*Mst1*Tbeta)-
2*(73917*MuSUSY+173666*Mst1*Tbeta)*pow2(Mt)+
23*(2195*MuSUSY+4902*Mst1*Tbeta)*pow2(Mst1)*pow2(s2t))+
pow3(Dmst12)*(-8*Mst1*s2t*(142987*MuSUSY-14232*Mst1*Tbeta)*pow2(Mt)-
Mt*(201940*MuSUSY-90723*Mst1*Tbeta)*pow2(Mst1)*pow2(s2t)+
4*(73917*MuSUSY+574156*Mst1*Tbeta)*pow3(Mt)+50485*Tbeta*pow3(s2t)*pow4(Mst1))+
4*Dmst12*(3*Mt*(24639*MuSUSY-75608*Mst1*Tbeta)+
92*Mst1*s2t*(817*MuSUSY+189*Mst1*Tbeta))*pow2(Mt)*pow4(Mst2)+
96*(1749*MuSUSY-3926*Mst1*Tbeta)*pow3(Mt)*pow6(Mst2)))))/pow6(Msq))/pow2(Mst1))/pow3(Mgl))/Tbeta-392*(-225*Mst1*Mt*pow2(Dmst12)*pow2(s2t)*(Dmst12*(102655-1000*lmMst1-6000*pow2(lmMst1)+
(Dmglst1*(284641+8696*lmMst1+1680*pow2(lmMst1)))/Mgl+
(-(pow2(Dmglst1)*(3532083+36328*lmMst1-47760*pow2(lmMst1)))-
(800*Dmsqst1*(5*(Dmglst1*Mgl+pow2(Dmglst1))+21*pow2(Mgl))*(Dmsqst1+pow2(Msq)))/pow4(Msq))/pow2(Mgl))+
2*pow2(Mst2)*(27220+
(10*Dmglst1*(33261-532*lmMst1-660*pow2(lmMst1)))/Mgl+
1080*pow2(lmMst1)+
(pow2(Dmglst1)*(1286791+5936*lmMst1-18120*pow2(lmMst1)))/pow2(Mgl)
+200*(lmMst1+(Dmsqst1*(5*(Dmglst1*Mgl+pow2(Dmglst1))+21*pow2(Mgl))*(Dmsqst1+pow2(Msq)))/(pow2(Mgl)*pow4(Msq)))))-
75*pow2(Mst1)*pow3(Dmst12)*pow3(s2t)*(14290-4540*lmMst1+
(5*Dmglst1*(84209+1264*lmMst1-240*pow2(lmMst1)))/Mgl+
240*pow2(lmMst1)+
(pow2(Dmglst1)*(1732531+16896*lmMst1-24840*pow2(lmMst1)))/pow2(Mgl)-(300*Dmsqst1*(7+24*lmMst1)*(Dmsqst1+pow2(Msq)))/pow4(Msq))+
s2t*(-13500*shiftst3*(3*pow2(Dmst12)*pow2(Mst2)*(-16*(2-lmMst1-z2)*pow2(Mt)+
(1-2*(lmMst1+z2))*pow2(Mst1)*pow2(s2t))-
(-16*(7-3*(lmMst1+z2))*pow2(Mt)+
(13-14*(lmMst1+z2))*pow2(Mst1)*pow2(s2t))*pow3(Dmst12)+
pow2(Mt)*(24*Dmst12*(3-2*(lmMst1+z2))*pow4(Mst2)-
24*(1-2*(lmMst1+z2))*pow6(Mst2)))+
(-135000*shiftst2*(Dmsqst1*(3-2*(lmMst1+z2))*pow2(Msq)*(-3*pow2(Dmst12)*pow2(Mst1)*pow2(Mst2)*pow2(s2t)-
(24*pow2(Mt)-2*pow2(Mst1)*pow2(s2t))*pow3(Dmst12)+
24*pow2(Mt)*(Dmst12*pow4(Mst2)+pow6(Mst2)))+
(1-2*(lmMst1+z2))*pow4(Msq)*(pow2(Mst1)*pow2(s2t)*(-3*pow2(Dmst12)*pow2(Mst2)+2*pow3(Dmst12))+
24*pow2(Mt)*(Dmst12*pow4(Mst2)+pow6(Mst2)))+
2*pow2(Dmsqst1)*(-3*pow2(Dmst12)*pow2(Mst2)*(4*pow2(Mt)+
(2-lmMst1-z2)*pow2(Mst1)*pow2(s2t))+
2*(-12*(1-lmMst1-z2)*pow2(Mt)+
(3-lmMst1-z2)*pow2(Mst1)*pow2(s2t))*pow3(Dmst12)+
24*(2-lmMst1-z2)*pow2(Mt)*(Dmst12*pow4(Mst2)+pow6(Mst2))))+
135000*shiftst1*((Dmsqst1*(3-2*(lmMst1+z2))*pow2(Msq)+
(1-2*(lmMst1+z2))*pow4(Msq))*(pow2(Mst1)*pow2(s2t)*(-3*pow2(Dmst12)*pow2(Mst2)+4*pow3(Dmst12))+
24*pow2(Mt)*pow6(Mst2))+
2*pow2(Dmsqst1)*((9-4*(lmMst1+z2))*pow2(Mst1)*pow2(s2t)*pow3(Dmst12)+
(2-lmMst1-z2)*(-3*pow2(Dmst12)*pow2(Mst1)*pow2(Mst2)*pow2(s2t)+
24*pow2(Mt)*pow6(Mst2)))))/pow4(Msq))+
(225*Dmst12*s2t*pow2(Mt)*(16*Dmglst1*Mgl*(4700*(pow2(Dmsqst1)+Dmsqst1*pow2(Msq))*(pow2(Dmst12)-Dmst12*pow2(Mst2)+pow4(Mst2))-
pow4(Msq)*(pow2(Dmst12)*(12383+4128*lmMst1-1260*pow2(lmMst1))-
Dmst12*(17539+574*lmMst1-1920*pow2(lmMst1))*pow2(Mst2)+
5*(4539-596*lmMst1-516*pow2(lmMst1))*pow4(Mst2)))+4*pow2(Dmglst1)*(19160*(pow2(Dmsqst1)+Dmsqst1*pow2(Msq))*(pow2(Dmst12)-Dmst12*pow2(Mst2)+pow4(Mst2))-
pow4(Msq)*(8*pow2(Dmst12)*(143196+2546*lmMst1-4785*pow2(lmMst1))-
Dmst12*(586073+9268*lmMst1-19200*pow2(lmMst1))*pow2(Mst2)+
2*(13289-916*lmMst1-60*pow2(lmMst1))*pow4(Mst2)))+pow2(Mgl)*(pow4(Msq)*(pow2(Dmst12)*(102747+10592*lmMst1-13888*pow2(lmMst1))-
20*Dmst12*(2071+296*lmMst1-600*pow2(lmMst1))*pow2(Mst2)-
160*(631-lmMst1+21*pow2(lmMst1))*pow4(Mst2))+
100*(Dmsqst1*pow2(Msq)*((346+288*lmMst1)*pow2(Dmst12)-
161*Dmst12*pow2(Mst2)-
24*(1+12*lmMst1)*pow4(Mst2))+
pow2(Dmsqst1)*(3*(47+96*lmMst1)*pow2(Dmst12)+
44*Dmst12*pow2(Mst2)-
(229+288*lmMst1)*pow4(Mst2)))))+
(2*pow3(Mt)*(15*pow2(Mgl)*(-(pow4(Msq)*(-(pow2(Dmst12)*(2284899-57376*lmMst1-87360*pow2(lmMst1))*pow2(Mst2))+
6*(1151533-4192*lmMst1-15520*pow2(lmMst1))*pow3(Dmst12)-
200*Dmst12*(11697+448*lmMst1+408*pow2(lmMst1))*pow4(Mst2)-
1600*(434-83*lmMst1+183*pow2(lmMst1))*pow6(Mst2)))+
4000*(pow2(Dmsqst1)*((65-6*lmMst1)*(-(pow2(Dmst12)*pow2(Mst2))+pow3(Dmst12)+
Dmst12*pow4(Mst2))+(91-12*lmMst1)*pow6(Mst2))+Dmsqst1*pow2(Msq)*((65-6*lmMst1)*(-(pow2(Dmst12)*pow2(Mst2))+pow3(Dmst12)+
Dmst12*pow4(Mst2))+6*(20-3*lmMst1)*pow6(Mst2))))+2*(Dmglst1*Mgl*(-(pow4(Msq)*(-21*pow2(Dmst12)*(5639279-7540*lmMst1-45600*pow2(lmMst1))*pow2(Mst2)+
(386728798-2268600*lmMst1-
4010400*pow2(lmMst1))*pow3(Dmst12)-
40*Dmst12*(3746977-48798*lmMst1-52380*pow2(lmMst1))*pow4(Mst2)-
6000*(9961-670*lmMst1+42*pow2(lmMst1))*pow6(Mst2)))+
1200*(Dmsqst1*pow2(Msq)*((557+120*lmMst1)*(-(pow2(Dmst12)*pow2(Mst2))+pow3(Dmst12)+
Dmst12*pow4(Mst2))+
50*(26+3*lmMst1)*pow6(Mst2))+
pow2(Dmsqst1)*((557+120*lmMst1)*(-(pow2(Dmst12)*pow2(Mst2))+pow3(Dmst12)+
Dmst12*pow4(Mst2))+
25*(44+3*lmMst1)*pow6(Mst2))))+
pow2(Dmglst1)*(pow4(Msq)*(-3*pow2(Dmst12)*(129193181+735100*lmMst1-
1336800*pow2(lmMst1))*pow2(Mst2)+
(119275604+4315560*lmMst1-
957600*pow2(lmMst1))*pow3(Dmst12)+
2*Dmst12*(327941741+47520*lmMst1-
3531600*pow2(lmMst1))*pow4(Mst2)+
80*(3777727-57348*lmMst1-52380*pow2(lmMst1))*pow6(Mst2))+
1200*(pow2(Dmsqst1)*((557+120*lmMst1)*(-(pow2(Dmst12)*pow2(Mst2))+pow3(Dmst12)+
Dmst12*pow4(Mst2))-
(136-165*lmMst1)*pow6(Mst2))+
Dmsqst1*pow2(Msq)*((557+120*lmMst1)*(-(pow2(Dmst12)*pow2(Mst2))+pow3(Dmst12)+
Dmst12*pow4(Mst2))+16*(4+15*lmMst1)*pow6(Mst2)))))))/Mst1)/(pow2(Mgl)*pow4(Msq)))))/(1.90512e8*pow6(Mst2));
}

/// calc coefficient O(as^2,log(mu^2/mt^2)^1) s12
double himalaya::H3::calc_coeff_as_2_log_1_s12(){
   return (MuSUSY*pow3(Mt)*((-((xDmglst1*pow3(Dmglst1)*(-(pow2(Dmst12)*pow2(Mst2)*(16800*Mst1*Mt*Tbeta*(pow2(Dmsqst1)+Dmsqst1*pow2(Msq))+
(-(Mt*(15680*MuSUSY+
3*(84421-13340*lmMst1)*Mst1*Tbeta))+
13440*s2t*Tbeta*pow2(Mst1))*pow4(Msq)))+
2*pow3(Dmst12)*(8400*Mst1*Mt*Tbeta*(pow2(Dmsqst1)+Dmsqst1*pow2(Msq))+
(-(Mt*(7840*MuSUSY+
3*(38903+5230*lmMst1)*Mst1*Tbeta))+
18480*s2t*Tbeta*pow2(Mst1))*pow4(Msq))+
4*Dmst12*(4200*Mst1*Mt*Tbeta*(pow2(Dmsqst1)+Dmsqst1*pow2(Msq))-
(Mt*(3920*MuSUSY+
3*(22759-9285*lmMst1)*Mst1*Tbeta)+
2520*s2t*Tbeta*pow2(Mst1))*pow4(Msq))*pow4(Mst2)+
140*Mt*(Mst1*Tbeta*(87*pow2(Dmsqst1)+162*Dmsqst1*pow2(Msq))-
(48*MuSUSY+(3197-1554*lmMst1)*Mst1*Tbeta)*pow4(Msq))*pow6(Mst2)))/pow4(Msq))+
560*Mgl*Mt*MuSUSY*(5*pow2(Mgl)*(2*(pow2(Dmst12)*pow2(Mst2)+pow3(Dmst12))-
6*Dmst12*pow4(Mst2)-9*pow6(Mst2))+
pow2(Dmglst1)*(-9*pow2(Dmst12)*pow2(Mst2)+
9*(pow3(Dmst12)+Dmst12*pow4(Mst2))-5*pow6(Mst2))-
10*Dmglst1*Mgl*(-(pow2(Dmst12)*pow2(Mst2))+pow3(Dmst12)+
Dmst12*pow4(Mst2)+3*pow6(Mst2))))/Tbeta-
(35*Mst1*(48*Mt*xDmglst1*xDmsqst1*pow3(Dmglst1)*pow3(Dmsqst1)*(-10*pow2(Dmst12)*pow2(Mst2)+
10*(pow3(Dmst12)+Dmst12*pow4(Mst2))+pow6(Mst2))-
2*pow3(Mgl)*(pow6(Msq)*(-(((623-222*lmMst1)*Mt-180*Mst1*s2t)*pow2(Dmst12)*pow2(Mst2))+
((421+486*lmMst1)*Mt-60*Mst1*s2t)*pow3(Dmst12)+
15*Dmst12*((55-62*lmMst1)*Mt-36*Mst1*s2t)*pow4(Mst2)+120*(29-11*lmMst1)*Mt*pow6(Mst2))+
300*Mt*(xDmsqst1*pow3(Dmsqst1)*(-(pow2(Dmst12)*pow2(Mst2))+pow3(Dmst12)+
Dmst12*pow4(Mst2)+pow6(Mst2))+
pow2(Dmsqst1)*pow2(Msq)*(-(pow2(Dmst12)*pow2(Mst2))+pow3(Dmst12)+
Dmst12*pow4(Mst2)+2*pow6(Mst2))+
Dmsqst1*pow4(Msq)*(-(pow2(Dmst12)*pow2(Mst2))+pow3(Dmst12)+
Dmst12*pow4(Mst2)+3*pow6(Mst2))))+
3*(Dmglst1*pow2(Mgl)*(-(pow6(Msq)*(((11+292*lmMst1)*Mt+160*Mst1*s2t)*pow2(Dmst12)*pow2(Mst2)-
((200-388*lmMst1)*Mt+80*Mst1*s2t)*pow3(Dmst12)+
2*Dmst12*((89-486*lmMst1)*Mt-120*Mst1*s2t)*pow4(Mst2)+40*(11-40*lmMst1)*Mt*pow6(Mst2)))+
Mt*(160*Dmst12*xDmsqst1*pow3(Dmsqst1)*(pow2(Dmst12)-Dmst12*pow2(Mst2)+pow4(Mst2))+
40*Dmsqst1*pow4(Msq)*(-4*pow2(Dmst12)*pow2(Mst2)+
4*(pow3(Dmst12)+Dmst12*pow4(Mst2))+5*pow6(Mst2))+20*pow2(Dmsqst1)*pow2(Msq)*(-8*pow2(Dmst12)*pow2(Mst2)+
8*(pow3(Dmst12)+Dmst12*pow4(Mst2))+5*pow6(Mst2))))+Mgl*pow2(Dmglst1)*(Mt*(160*Dmsqst1*pow4(Msq)*(-(pow2(Dmst12)*pow2(Mst2))+pow3(Dmst12)+
Dmst12*pow4(Mst2)+2*pow6(Mst2))+
40*xDmsqst1*pow3(Dmsqst1)*(-4*pow2(Dmst12)*pow2(Mst2)+
4*(pow3(Dmst12)+Dmst12*pow4(Mst2))+3*pow6(Mst2))+20*pow2(Dmsqst1)*pow2(Msq)*(-8*pow2(Dmst12)*pow2(Mst2)+
8*(pow3(Dmst12)+Dmst12*pow4(Mst2))+11*pow6(Mst2)))-pow6(Msq)*(-(((1337-388*lmMst1)*Mt-112*Mst1*s2t)*pow2(Dmst12)*pow2(Mst2))+
2*Dmst12*((763-534*lmMst1)*Mt-20*Mst1*s2t)*pow4(Mst2)+
4*(((287+73*lmMst1)*Mt-46*Mst1*s2t)*pow3(Dmst12)+
(539-486*lmMst1)*Mt*pow6(Mst2)))))))/pow6(Msq)))/(14175.*pow2(Mst1)*pow3(Mgl)*pow6(Mst2));
}

/// calc coefficient O(as^2,log(mu^2/mt^2)^2) s12
double himalaya::H3::calc_coeff_as_2_log_2_s12(){
   return (-16*MuSUSY*pow4(Mt)*(5*pow2(Mst2)*pow3(Mgl)*(pow2(Dmst12)-2*Dmst12*pow2(Mst2)-6*pow4(Mst2))+
Dmst12*Mgl*pow2(Dmglst1)*(pow2(Dmst12)-2*Dmst12*pow2(Mst2)+3*pow4(Mst2))+
xDmglst1*pow3(Dmglst1)*(-2*pow2(Dmst12)*pow2(Mst2)+pow3(Dmst12)+
3*Dmst12*pow4(Mst2)+4*pow6(Mst2))-
Dmglst1*pow2(Mgl)*(-(pow2(Dmst12)*pow2(Mst2))+2*pow3(Dmst12)+
10*pow6(Mst2))))/(45.*Mst1*pow3(Mgl)*pow6(Mst2));
}
