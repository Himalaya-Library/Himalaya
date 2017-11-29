#include "H5.hpp"
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
 * 	@param Dmglst1 a double Mgl - Mst1
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
himalaya::H5::H5(std::map<unsigned int, unsigned int> flagMap, double Al4p, double beta, double Dmglst1,
		 double lmMt, double lmMst1, double lmMst2, double lmMsq, double Mt, double Mst1,
		 double Mst2, double Msq, double MuSUSY,
		 double s2t,
		 int mdrFlag, int oneLoopFlag, int twoLoopFlag, int threeLoopFlag){
   // abbrev for tan(beta) and sin(beta)
   Tbeta = tan(beta);
   Cbeta = cos(beta);
   Sbeta = sin(beta);
   this -> Dmglst1 = Dmglst1;
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
   // expansion flags
   xDmglst1 = flagMap.at(HierarchyCalculator::xxDmglst1);
   xMsq = flagMap.at(HierarchyCalculator::xxMsq);
   
   s1 = 
   #include "../hierarchies/h5/sigS1Full.inc"
   ;
   s2 = 
   #include "../hierarchies/h5/sigS2Full.inc"
   ;
   s12 = 
   #include "../hierarchies/h5/sigS12Full.inc"
   ;
}

/**
 * 	@return The diagonal (1, 1) matrix element of the Higgs mass matrix as a double for the hierarchy 'H5'
 */
double himalaya::H5::getS1(){
   return s1;
}

/**
 * 	@return The diagonal (2, 2) matrix element of the Higgs mass matrix as a double for the hierarchy 'H5'
 */
double himalaya::H5::getS2(){
   return s2;
}

/**
 * 	@return The off-diagonal (1, 2) = (2, 1) matrix element of the Higgs mass matrix as a double for the hierarchy 'H5'
 */
double himalaya::H5::getS12(){
   return s12;
}

/**
 * 	@return returns the constant term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double himalaya::H5::calc_at_as2_no_logs(){

   const double result =
      ((-3888*Cbeta*MuSUSY*Sbeta*((s2t*pow3(Mt)*((-675*(334425*pow2(Mst1)*pow2(
        Mst2) + 377300*Dmglst1*pow3(Mst1) - 187967*pow4(Mst1) - 103583*pow4(
        Mst2)))/pow4(Msq) + (343000*pow2(Mst1)*((-160*OepS2 + 7506*S2)*pow2(
        Mst1) + 21*(-8*OepS2 + 405*S2)*pow2(Mst2)))/pow4(Mst2) + (pow2(Mst1)*(-
        3087000*Dmglst1*(-17783 + 144*B4 + 72*DN)*Mst1 - 514500*(-226991 + 432*
        B4 + 216*DN)*pow2(Dmglst1) + 2622405502*pow2(Mst1) + 32409593916*pow2(
        Mst2) - (123480*(82500*Dmglst1*Mst1 + 123750*pow2(Dmglst1) + 20625*
        pow2(Mst1) - 2017*pow2(Mst2))*pow2(Mst2))/pow2(Msq)))/pow4(Mst2)))/
        6.251175e7 - Mt*pow3(s2t)*(-((((-606133*Dmglst1*Mst1)/3888. + (
        119.47646604938272 - (64*B4)/9. - (32*DN)/9.)*pow2(Dmglst1) - (10*pow2(
        Msq))/3. - (478241812219*pow2(Mst1))/2.000376e9)*pow2(Mst1))/pow2(Mst2)
        ) + (pow2(Mst1)*((-1328*OepS2 + 64098*S2)*pow2(Mst1) + 177*(-8*OepS2 +
        405*S2)*pow2(Mst2)))/(729.*pow2(Mst2)) + ((3200*Dmglst1 + 787*Mst1)*
        pow3(Mst1))/(45.*pow2(Msq)) + pow2(Mst1)*(110.18859670781893 - (40*B4)/
        9. - (2*DN)/9. + (3050*pow2(Dmglst1))/(27.*pow2(Msq)) - (5*(2135*pow2(
        Dmglst1) + 568*pow2(Msq))*pow2(Mst2))/(216.*pow4(Msq))) + (pow2(Mst2)*(
        1543500*pow2(Dmglst1)*pow2(Mst2) - 86264500*Dmglst1*pow3(Mst1) -
        20222907*pow4(Mst1)))/(3.7044e6*pow4(Msq)) + (313/(45.*pow2(Msq)) - 1/(
        3.*pow2(Mst1)) + ((5*Dmglst1*Mst1)/6. + (1338719*pow2(Mst1))/1.2348e6)/
        pow4(Msq))*pow4(Mst2)) + (pow2(Mst1)*pow3(Mt)*(Dmglst1*Mt*(52260*pow2(
        Mst1)*pow2(Mst2) + 105360*pow2(Msq)*(pow2(Mst1) + pow2(Mst2)) + 65212*
        pow4(Msq) - 21135*pow4(Mst2)) + 3*Mst1*Mt*(2700*pow2(Mst1)*pow2(Mst2) +
        7920*pow2(Msq)*(pow2(Mst1) + pow2(Mst2)) - 6988*pow4(Msq) - 3345*pow4(
        Mst2)) + 45*pow2(Dmglst1)*(3920*Mst1*Mt*pow2(Msq) + 3132*Mst1*Mt*pow2(
        Mst2) - 285*s2t*pow4(Mst2))))/(243.*pow4(Msq)*pow4(Mst2)) - pow2(Mt)*
        pow2(s2t)*((1660*Dmglst1*pow2(Mst1))/(3.*pow2(Msq)) + ((-50760*Mst1*
        pow2(Dmglst1) + 3*Dmglst1*((10021 + 5328*B4 - 72*DN)*pow2(Msq) - 4320*
        pow2(Mst1)) + Mst1*((6277 + 7776*B4)*pow2(Msq) - 2160*pow2(Mst1)))*
        pow2(Mst1))/(81.*pow2(Msq)*pow2(Mst2)) + pow3(Mst1)*((5*(1001*pow2(
        Dmglst1) + 284*pow2(Msq)))/(9.*pow4(Msq)) + (Dmglst1*(
        2730.8603395061727 + (592*B4)/3. - (8*DN)/3.)*Mst1 + (5530.797839506173
         + 152*B4 - 4*DN)*pow2(Dmglst1) + (557.5486111111111 + 96*B4)*pow2(
        Mst1))/pow4(Mst2)) + (5*(-118*pow2(Mst2)*pow3(Mst1) + 11*Mst1*pow4(
        Mst2) + Dmglst1*(-262*pow2(Mst1)*pow2(Mst2) + 5159*pow4(Mst1) + 11*
        pow4(Mst2)) + 911*pow5(Mst1)))/(108.*pow4(Msq)))) - (pow2(Mt)*pow2(
        MuSUSY)*(pow2(Mst1)*(8232000*(Dmglst1*(-41907 + 35640*B4 - 324*DN) + (-
        107299 + 19224*B4 + 108*DN)*Mst1)*Mt*s2t + 1778112000*(-139 + 2*B4 +
        DN)*pow2(Mt) + (-514500*Dmglst1*(-437365 + 10368*B4 + 8640*DN)*Mst1 +
        771750*(-335123 + 14976*B4 + 6336*DN)*pow2(Dmglst1) + (531403689547 -
        16892064000*B4 + 444528000*DN - 8122240000*OepS2 + 306576144000*S2)*
        pow2(Mst1) - 16464*(-3228977 + 1026000*B4 - 27000*DN + 272000*OepS2 -
        7938000*S2)*pow2(Mst2))*pow2(s2t))*pow4(Msq) + 14817600*s2t*pow2(Msq)*
        pow2(Mst1)*(13600*Mst1*Mt*pow2(Mst2) + 1525*s2t*pow2(Mst1)*pow2(Mst2) +
        50*pow2(Dmglst1)*(116*Mst1*Mt + 225*s2t*pow2(Mst2)) + 100*Dmglst1*(348*
        Mt*pow2(Mst1) + 492*Mt*pow2(Mst2) + 55*Mst1*s2t*pow2(Mst2)) + 11200*Mt*
        pow3(Mst1) - 836*s2t*pow4(Mst2)) + 180*s2t*pow2(Mst2)*(257250*pow2(
        Dmglst1)*pow2(Mst1)*(15728*Mst1*Mt - 2117*s2t*pow2(Mst2)) + 171500*
        Dmglst1*pow2(Mst1)*(9816*Mt*pow2(Mst1) - 502*Mt*pow2(Mst2) - 1455*Mst1*
        s2t*pow2(Mst2)) - 36701000*Mt*pow2(Mst2)*pow3(Mst1) - 41220947*s2t*
        pow2(Mst2)*pow4(Mst1) + 19447774*s2t*pow2(Mst1)*pow4(Mst2) + 275772000*
        Mt*pow5(Mst1) + 7399303*s2t*pow6(Mst2))))/(257250.*pow4(Msq)*pow4(Mst2)
        ) - (4*(-108*Mst1*s2t*pow2(Sbeta)*(240*Dmglst1*s2t*z4*pow2(Mt) + 2*Mt*
        pow2(s2t)*(30*pow2(Dmglst1)*(41*z4 + 8*pow2(z2)) + pow2(Mst2)*(-241*z4
        + 32*pow2(z2))) + 32*(-z4 + 2*pow2(z2))*pow3(Mt) - Dmglst1*pow2(Mst2)*(
        265*z4 + 4*pow2(z2))*pow3(s2t))*pow4(Mst2) + 27*s2t*pow2(Sbeta)*pow4(
        Mst2)*(-32*Dmglst1*pow2(Mt)*(-4*Mt*z4 + 15*Dmglst1*s2t*z4 + 8*Mt*pow2(
        z2)) + 2*s2t*pow2(Mst2)*(4*Dmglst1*Mt*s2t*(241*z4 - 32*pow2(z2)) +
        pow2(Dmglst1)*pow2(s2t)*(265*z4 + 4*pow2(z2)) + pow2(Mt)*(-8*T1ep + 36*
        z4 + 10*pow2(z2))) + (4*T1ep - 242*z4 + 35*pow2(z2))*pow3(s2t)*pow4(
        Mst2)) + s2t*pow4(Mst1)*(-664*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*pow2(
        s2t)*(4*T1ep + 2*z4 + 3*pow2(z2)) - 8*s2t*pow2(Mt)*(-(pow2(Mst2)*pow2(
        Sbeta)*(4*T1ep + 2*z4 + 3*pow2(z2))) - 324*Cbeta*Dmglst1*MuSUSY*Sbeta*(
        187*z4 + 76*pow2(z2)) + 2*pow2(MuSUSY)*(-1 + pow2(Sbeta))*(740*T1ep +
        127*z4 + 1527*pow2(z2))) + 320*Cbeta*MuSUSY*Sbeta*(4*T1ep + 2*z4 + 3*
        pow2(z2))*pow3(Mt) - pow2(Sbeta)*(44*T1ep + 6124*z4 + 1113*pow2(z2))*
        pow3(s2t)*pow4(Mst2)) + 108*pow3(Mst1)*(4*MuSUSY*pow2(Mt)*pow2(s2t)*(
        108*Cbeta*Sbeta*pow2(Mst2)*(-z4 + 2*pow2(z2)) - 2*Dmglst1*MuSUSY*(-1 +
        pow2(Sbeta))*(265*z4 + 4*pow2(z2)) + 45*Cbeta*Sbeta*pow2(Dmglst1)*(41*
        z4 + 8*pow2(z2))) + 8*MuSUSY*s2t*(60*Cbeta*Dmglst1*Sbeta*z4 + MuSUSY*(-
        1 + pow2(Sbeta))*(-313*z4 + 176*pow2(z2)))*pow3(Mt) - 2*Mt*pow2(Sbeta)*
        (169*z4 + 112*pow2(z2))*pow3(s2t)*pow4(Mst2) + 32*Cbeta*MuSUSY*Sbeta*(-
        z4 + 2*pow2(z2))*pow4(Mt) - Dmglst1*pow2(Sbeta)*(265*z4 + 4*pow2(z2))*
        pow4(Mst2)*pow4(s2t)) + 46656*Cbeta*MuSUSY*Sbeta*pow2(Mt)*pow2(s2t)*(-
        z4 + 2*pow2(z2))*pow5(Mst1) + 6*pow2(Mst1)*(72*MuSUSY*pow2(Mt)*(60*
        MuSUSY*z4*pow2(Mt)*(-1 + pow2(Sbeta)) + s2t*pow2(Dmglst1)*(60*Cbeta*Mt*
        Sbeta*z4 - MuSUSY*s2t*(-1 + pow2(Sbeta))*(25*z4 + 4*pow2(z2))) + 2*
        Dmglst1*Mt*(4*Cbeta*Mt*Sbeta*(-z4 + 2*pow2(z2)) + 3*MuSUSY*s2t*(-1 +
        pow2(Sbeta))*(169*z4 + 112*pow2(z2)))) + 8*Mt*MuSUSY*s2t*pow2(Mst2)*(
        1080*Cbeta*Sbeta*z4*pow2(Dmglst1)*pow2(s2t) + 7*Cbeta*Sbeta*pow2(Mt)*(
        4*T1ep + 2*z4 + 3*pow2(z2)) + Mt*s2t*(54*Cbeta*Dmglst1*Sbeta*(187*z4 +
        76*pow2(z2)) - MuSUSY*(-1 + pow2(Sbeta))*(136*T1ep - 13*z4 + 426*pow2(
        z2)))) - Sbeta*pow2(s2t)*(2*Cbeta*Mt*MuSUSY*s2t*(236*T1ep + 2152*z4 +
        537*pow2(z2)) + Sbeta*(9*pow2(Dmglst1)*pow2(s2t)*(505*z4 + 4*pow2(z2))
        + 2*pow2(Mt)*(20*T1ep + 190*z4 + 87*pow2(z2)) + 36*Dmglst1*Mt*s2t*(989*
        z4 + 272*pow2(z2))))*pow4(Mst2) + pow2(Sbeta)*(100*T1ep + 2165*z4 +
        111*pow2(z2))*pow4(s2t)*pow6(Mst2))))/pow4(Mst2) + 16*z3*((2*Mt*MuSUSY*
        pow2(Mst1)*(-4*s2t*(31050*Mst1*MuSUSY + 108*Dmglst1*(231*MuSUSY + 20*
        Cbeta*Mst1*Sbeta) + 14472*Cbeta*Sbeta*pow2(Dmglst1) - 734*Cbeta*Sbeta*
        pow2(Mst1) - 1065*Cbeta*Sbeta*pow2(Mst2))*pow2(Mt) + Mt*pow2(s2t)*(810*
        (229*MuSUSY - 727*Cbeta*Mst1*Sbeta)*pow2(Dmglst1) + 135212*MuSUSY*pow2(
        Mst1) + 52698*MuSUSY*pow2(Mst2) + 71604*Cbeta*Mst1*Sbeta*pow2(Mst2) -
        27*Dmglst1*(-6112*Mst1*MuSUSY + 6549*Cbeta*Sbeta*pow2(Mst1) - 2424*
        Cbeta*Sbeta*pow2(Mst2)) + 29565*Cbeta*Sbeta*pow3(Mst1)) + 216*(46*
        MuSUSY + Cbeta*(470*Dmglst1 + 129*Mst1)*Sbeta)*pow3(Mt) - Cbeta*Sbeta*
        pow2(Mst2)*(76842*Dmglst1*Mst1 + 90477*pow2(Dmglst1) + 41257*pow2(Mst1)
        + 21921*pow2(Mst2))*pow3(s2t)))/pow4(Mst2) + 243*pow2(Sbeta)*(Mt*pow3(
        s2t)*(6*Mst1*pow2(Dmglst1) - (4*pow2(Mst1)*(346*Dmglst1 + 309*Mst1 - (
        3284*Mst1*pow2(Dmglst1))/pow2(Mst2)))/9. - ((232*Dmglst1)/9. + (532*
        Mst1)/9. - (150*pow2(Dmglst1))/Mst1)*pow2(Mst2) + (2*(997*Dmglst1 +
        173*Mst1)*pow4(Mst1))/(3.*pow2(Mst2))) + ((1325*pow2(Dmglst1)*pow2(
        Mst1))/8. + ((7290*Dmglst1*Mst1 + 11799*pow2(Dmglst1) + 11662*pow2(
        Mst1))*pow2(Mst2))/324. + (5077*Dmglst1*pow3(Mst1))/36. + (9668*pow4(
        Mst1))/243. + (9.11111111111111 - (65*Dmglst1)/(12.*Mst1) - (127*pow2(
        Dmglst1))/(8.*pow2(Mst1)))*pow4(Mst2))*pow4(s2t) + ((-4*s2t*pow3(Mt)*(
        pow2(Dmglst1)*(1752*pow2(Mst1)*pow2(Mst2) + 16888*pow4(Mst1) + 1843*
        pow4(Mst2)) + 4*Dmglst1*((507*pow2(Mst2) - 462*pow2(MuSUSY))*pow3(Mst1)
        - 37*Mst1*pow4(Mst2) + 2196*pow5(Mst1)) + 4*((183*pow2(Mst2) - 575*
        pow2(MuSUSY))*pow4(Mst1) - 54*pow2(Mst1)*pow4(Mst2) + 448*pow6(Mst1))))
        /(9.*Mst1) + ((pow4(Mt)*(pow2(Dmglst1)*(12432*pow2(Mst1)*pow2(Mst2) +
        16208*pow4(Mst1) - 1805*pow4(Mst2)) + 2*Dmglst1*(3088*pow2(Mst2)*pow3(
        Mst1) - 107*Mst1*pow4(Mst2) + 5264*pow5(Mst1)) + 8*(4*(36*pow2(Mst2) -
        23*pow2(MuSUSY))*pow4(Mst1) - 163*pow2(Mst1)*pow4(Mst2) + 317*pow6(
        Mst1))))/9. + (pow2(Mt)*pow2(s2t)*(27*pow2(Dmglst1)*(24*(597*pow2(Mst2)
        - 1145*pow2(MuSUSY))*pow4(Mst1) - 13079*pow2(Mst1)*pow4(Mst2) + 3039*
        pow6(Mst2)) + 216*Dmglst1*(-373*pow3(Mst1)*pow4(Mst2) + (398*pow2(Mst2)
        - 3056*pow2(MuSUSY))*pow5(Mst1) + 55*Mst1*pow6(Mst2)) - 4*(pow4(Mst1)*(
        52698*pow2(Mst2)*pow2(MuSUSY) + 969*pow4(Mst2)) + (-662*pow2(Mst2) +
        135212*pow2(MuSUSY))*pow6(Mst1) + 1161*pow2(Mst1)*pow6(Mst2))))/486.)/
        pow2(Mst1))/pow4(Mst2))) + (pow2(Sbeta)*(-2058000*Mt*pow3(s2t)*(576*(
        127 + 342*B4 - 9*DN)*Mst1*pow2(Dmglst1) + 24*Dmglst1*(54053 + 9432*B4 -
        180*DN)*pow2(Mst1) + (8*(-360*Mst1*(321*Dmglst1 + 83*Mst1) - 213480*
        pow2(Dmglst1) + (132407 + 11880*B4 - 108*DN)*pow2(Msq))*pow3(Mst1))/
        pow2(Msq) + ((3058187*Dmglst1*Mst1 + 6248762*pow2(Dmglst1) + 622151*
        pow2(Mst1))*pow3(Mst1))/pow2(Mst2) + (36*pow2(Mst2)*(pow2(Mst1)*(5*
        pow2(Mst1)*(343*pow2(Mst1) - 43*pow2(Mst2)) + 80*pow2(Msq)*(74*pow2(
        Mst1) - 3*pow2(Mst2)) + 6*(-4439 + 136*B4 + 4*DN)*pow4(Msq)) + Dmglst1*
        Mst1*(65*pow2(Mst1)*(139*pow2(Mst1) - 7*pow2(Mst2)) + 240*pow2(Msq)*(
        84*pow2(Mst1) - pow2(Mst2)) + 6*(-3779 + 136*B4 + 4*DN)*pow4(Msq)) +
        20*pow2(Dmglst1)*(1244*pow2(Msq)*pow2(Mst1) - 18*pow2(Mst1)*pow2(Mst2)
        + 1175*pow4(Msq) + 1019*pow4(Mst1))))/(Mst1*pow4(Msq))) + (pow4(s2t)*(-
        29635200*pow2(Msq)*pow2(Mst1)*(125*pow2(Dmglst1)*(77*pow2(Mst1) - 16*
        pow2(Mst2)) - 1357*pow2(Mst1)*pow2(Mst2) + 50*Dmglst1*(-41*Mst1*pow2(
        Mst2) + 137*pow3(Mst1)) + 2068*pow4(Mst1))*pow4(Mst2) + pow2(Mst2)*
        pow4(Msq)*(-771750*pow2(Dmglst1)*(-48*(-5885 + 72*B4 + 60*DN)*pow2(
        Mst1)*pow2(Mst2) + (-438203 + 21888*B4 + 12096*DN)*pow4(Mst1) - 153960*
        pow4(Mst2)) + pow2(Mst1)*(131712*(-2938294 + 6750*B4 + 10125*DN +
        25000*OepS2 - 1994625*S2)*pow2(Mst1)*pow2(Mst2) - (257823187891 +
        8890560000*B4 + 444528000*DN + 241472000*OepS2 - 20818728000*S2)*pow4(
        Mst1) + 3087000*(53749 + 2592*B4 - 288*DN + 192*OepS2 + 21384*S2)*pow4(
        Mst2)) - 514500*Dmglst1*(-144*(1151 + 72*B4 + 60*DN)*pow2(Mst2)*pow3(
        Mst1) - 1512*Mst1*pow4(Mst2) + (773389 + 10368*B4 + 8640*DN)*pow5(Mst1)
        )) + 6667920000*(pow2(Mst1) + pow2(Mst2))*pow2(pow2(Mst1) - pow2(Mst2))
        *pow6(Msq) + 1080*(44675750*Dmglst1*Mst1 + 92309875*pow2(Dmglst1) +
        12119532*pow2(Mst1))*pow4(Mst1)*pow6(Mst2)))/(pow2(Mst1)*pow2(Mst2)*
        pow4(Msq)) - (192*pow4(Mt)*(15*pow2(Mst1)*pow4(Mst2)*(59167500*pow2(
        Mst1)*pow2(Mst2) + 160842737*pow4(Mst1) + 40792737*pow4(Mst2)) -
        555660000*pow2(Mst2)*(pow2(Mst1) + pow2(Mst2))*pow6(Msq) + 4*pow4(Msq)*
        (5488*(614897*pow2(Mst2) - 3375*(-139 + 2*B4 + DN)*pow2(MuSUSY))*pow4(
        Mst1) - 1350691125*pow2(Mst1)*pow4(Mst2) + 6336529126*pow6(Mst1) -
        13891500*pow6(Mst2)) + 82320*pow2(Msq)*(37017*pow4(Mst1)*pow4(Mst2) +
        20625*pow2(Mst2)*pow6(Mst1) + 13517*pow2(Mst1)*pow6(Mst2)) + 3430*pow2(
        Dmglst1)*(pow4(Msq)*(22632928*pow2(Mst1)*pow2(Mst2) + 45232682*pow4(
        Mst1) - 2120028*pow4(Mst2)) + 720*pow2(Msq)*(4125*pow2(Mst2)*pow4(Mst1)
        + 1324*pow2(Mst1)*pow4(Mst2)) + 45*(36871*pow4(Mst1)*pow4(Mst2) + 5750*
        pow2(Mst1)*pow6(Mst2))) + 171500*Dmglst1*(2*pow4(Msq)*(128186*pow2(
        Mst2)*pow3(Mst1) + 6645*Mst1*pow4(Mst2) + 299955*pow5(Mst1)) + 360*
        pow2(Msq)*(81*pow3(Mst1)*pow4(Mst2) + 110*pow2(Mst2)*pow5(Mst1)) + 75*(
        419*pow4(Mst2)*pow5(Mst1) + 138*pow3(Mst1)*pow6(Mst2)))))/(pow2(Mst1)*
        pow4(Msq)*pow4(Mst2)) + (32928000*s2t*pow3(Mt)*(4*pow2(Dmglst1)*(pow4(
        Msq)*(32430*pow2(Mst1)*pow2(Mst2) + 348947*pow4(Mst1) + 89475*pow4(
        Mst2)) - 180*pow2(Msq)*((111*pow2(Mst2) - 29*pow2(MuSUSY))*pow4(Mst1) -
        203*pow2(Mst1)*pow4(Mst2)) + 45*(pow4(Mst1)*(983*pow2(Mst2)*pow2(
        MuSUSY) + 742*pow4(Mst2)) - 112*pow2(Mst1)*pow6(Mst2))) + pow2(Mst1)*(
        pow4(Msq)*(pow2(Mst1)*(996*pow2(Mst2) + 2*(-107299 + 19224*B4 + 108*DN)
        *pow2(MuSUSY)) - 81172*pow4(Mst1) + 13356*pow4(Mst2)) - 120*pow2(Msq)*(
        6*(3*pow2(Mst2) - 56*pow2(MuSUSY))*pow4(Mst1) - 4*pow2(Mst1)*(102*pow2(
        Mst2)*pow2(MuSUSY) + 55*pow4(Mst2)) + 67*pow6(Mst2)) + 15*(12*pow4(
        Mst1)*(67*pow2(Mst2)*pow2(MuSUSY) + 100*pow4(Mst2)) - pow2(Mst1)*(107*
        pow2(MuSUSY)*pow4(Mst2) + 448*pow6(Mst2)) - 338*pow8(Mst2))) + Dmglst1*
        Mst1*(pow4(Msq)*(pow2(Mst1)*(98884*pow2(Mst2) + 6*(-13969 + 11880*B4 -
        108*DN)*pow2(MuSUSY)) + 216148*pow4(Mst1) + 84084*pow4(Mst2)) - 120*
        pow2(Msq)*(18*(11*pow2(Mst2) - 58*pow2(MuSUSY))*pow4(Mst1) - 12*pow2(
        Mst1)*(123*pow2(Mst2)*pow2(MuSUSY) + 64*pow4(Mst2)) + 67*pow6(Mst2)) +
        15*(12*pow4(Mst1)*(409*pow2(Mst2)*pow2(MuSUSY) + 368*pow4(Mst2)) -
        pow2(Mst1)*(251*pow2(MuSUSY)*pow4(Mst2) + 1344*pow6(Mst2)) - 338*pow8(
        Mst2)))))/(Mst1*pow4(Msq)*pow4(Mst2)) - (8*pow2(Mt)*pow2(s2t)*(
        6667920000*pow2(-(Mst2*pow2(Mst1)) + pow3(Mst2))*pow6(Msq) + 514500*
        Dmglst1*(180*pow3(Mst1)*pow4(Mst2)*(pow2(Mst1)*(-17*pow2(Mst2) + 485*
        pow2(MuSUSY)) + 39*pow4(Mst2)) + pow4(Msq)*(8*(-156781 + 864*B4 + 432*
        DN)*pow3(Mst1)*pow4(Mst2) + (506360*pow2(Mst2) + (-437365 + 10368*B4 +
        8640*DN)*pow2(MuSUSY))*pow5(Mst1) - 22752*Mst1*pow6(Mst2)) + 1440*pow2(
        Msq)*((-110*pow2(Mst2)*pow2(MuSUSY) + 131*pow4(Mst2))*pow5(Mst1) - 21*
        pow3(Mst1)*pow6(Mst2))) + pow4(Msq)*(8232*pow2(Mst2)*((-32617079 +
        108000*DN + 20000*OepS2 - 2470500*S2)*pow2(Mst2) + 2*(-3228977 +
        1026000*B4 - 27000*DN + 272000*OepS2 - 7938000*S2)*pow2(MuSUSY))*pow4(
        Mst1) + ((238297507312 - 21952000*OepS2 + 2741256000*S2)*pow2(Mst2) + (
        -531403689547 + 16892064000*B4 - 444528000*DN + 8122240000*OepS2 -
        306576144000*S2)*pow2(MuSUSY))*pow6(Mst1) - 3087000*(-4501 + 288*DN -
        96*OepS2 + 972*S2)*pow2(Mst1)*pow6(Mst2) + 666792000*pow8(Mst2)) +
        987840*pow2(Msq)*((-22875*pow2(Mst2)*pow2(MuSUSY) + 22642*pow4(Mst2))*
        pow6(Mst1) + 11*pow4(Mst1)*(1140*pow2(MuSUSY)*pow4(Mst2) + 281*pow6(
        Mst2)) - 5108*pow2(Mst1)*pow8(Mst2)) - 180*(pow6(Mst1)*(-41220947*pow2(
        MuSUSY)*pow4(Mst2) + 15671760*pow6(Mst2)) + 7399303*pow2(Mst1)*pow2(
        MuSUSY)*pow8(Mst2) - 2*pow4(Mst1)*(-9723887*pow2(MuSUSY)*pow6(Mst2) +
        6570120*pow8(Mst2))) + 257250*pow2(Dmglst1)*(21600*pow2(Msq)*pow2(Mst1)
        *pow2(Mst2)*(3*pow2(Mst1)*(7*pow2(Mst2) - 10*pow2(MuSUSY)) + pow4(Mst2)
        ) + pow4(Msq)*(3*(23176*pow2(Mst2) + (335123 - 14976*B4 - 6336*DN)*
        pow2(MuSUSY))*pow4(Mst1) + 8*(-486235 + 864*B4 + 432*DN)*pow2(Mst1)*
        pow4(Mst2) + 679248*pow6(Mst2)) + 180*(pow4(Mst1)*(2117*pow2(MuSUSY)*
        pow4(Mst2) + 531*pow6(Mst2)) + 39*pow2(Mst1)*pow8(Mst2)))))/(pow2(Mst1)
        *pow4(Msq)*pow4(Mst2))))/2.058e6 + z2*(-(pow2(Sbeta)*(144*Mt*pow3(s2t)*
        (14244*Mst1*pow2(Dmglst1) - (10*(288*Dmglst1*Mst1 + 648*pow2(Dmglst1) -
        271*pow2(Msq) + 48*pow2(Mst1))*pow3(Mst1))/pow2(Msq) - ((1397*Dmglst1*
        Mst1 - 62*pow2(Dmglst1) + 737*pow2(Mst1))*pow3(Mst1))/pow2(Mst2) + 6*(
        1941*Dmglst1*pow2(Mst1) + pow2(Mst2)*((Dmglst1 + Mst1)*(-627 + (120*
        Dmglst1*Mst1)/pow2(Msq)) + (5*(15*Dmglst1*Mst1 + 30*pow2(Dmglst1) + 8*
        pow2(Msq) + 3*pow2(Mst1))*pow3(Mst1))/pow4(Msq)))) + ((pow2(Mst2)*pow4(
        s2t)*(1080*pow2(Msq)*pow2(Mst1)*(pow2(Dmglst1)*(60*pow2(Mst1) - 3*pow2(
        Mst2)) - 3*pow2(Mst1)*pow2(Mst2) + Dmglst1*(-6*Mst1*pow2(Mst2) + 40*
        pow3(Mst1)) + 10*pow4(Mst1))*pow4(Mst2) + 6480*(pow2(Mst1) + pow2(Mst2)
        )*pow2(pow2(Mst1) - pow2(Mst2))*pow6(Msq) + pow2(Mst2)*pow4(Msq)*(-
        68844*pow2(Mst2)*pow4(Mst1) + 162*pow2(Dmglst1)*(-388*pow2(Mst1)*pow2(
        Mst2) + 1655*pow4(Mst1) - 66*pow4(Mst2)) + 35316*pow2(Mst1)*pow4(Mst2)
        + 36*Dmglst1*(-3492*pow2(Mst2)*pow3(Mst1) - 594*Mst1*pow4(Mst2) + 7153*
        pow5(Mst1)) + 98497*pow6(Mst1)) - 1350*(4*Dmglst1*Mst1 + 6*pow2(
        Dmglst1) + pow2(Mst1))*pow4(Mst1)*pow6(Mst2)))/pow2(Mst1) + 1152*s2t*
        pow3(Mt)*(-(Mst1*(120*pow2(Msq)*pow2(Mst2)*(pow2(Mst1)*(-pow2(Mst2) +
        pow2(MuSUSY)) - 3*pow4(Mst1) + pow4(Mst2)) + pow4(Msq)*(pow2(Mst1)*(
        2778*pow2(Mst2) - 2407*pow2(MuSUSY)) + 4589*pow4(Mst1) + 26*pow4(Mst2))
        + 15*pow2(Mst2)*((-6*pow2(Mst2) + 3*pow2(MuSUSY))*pow4(Mst1) + 4*pow2(
        Mst1)*pow4(Mst2) + 2*pow6(Mst2)))) - Dmglst1*(pow4(Msq)*(pow2(Mst1)*(
        8300*pow2(Mst2) + 2061*pow2(MuSUSY)) + 22119*pow4(Mst1) + 26*pow4(Mst2)
        ) - 120*pow2(Msq)*(3*(5*pow2(Mst2) + 2*pow2(MuSUSY))*pow4(Mst1) + pow2(
        Mst1)*(-3*pow2(Mst2)*pow2(MuSUSY) + pow4(Mst2)) - pow6(Mst2)) + 15*
        pow2(Mst2)*((-14*pow2(Mst2) + 15*pow2(MuSUSY))*pow4(Mst1) + 12*pow2(
        Mst1)*pow4(Mst2) + 2*pow6(Mst2))) + pow2(Dmglst1)*(360*pow2(Msq)*(10*
        pow2(Mst2) + 7*pow2(MuSUSY))*pow3(Mst1) - (8283*Mst1*pow2(Mst2) +
        42755*pow3(Mst1))*pow4(Msq) + 90*(pow3(Mst1)*(-5*pow2(Mst2)*pow2(
        MuSUSY) + 2*pow4(Mst2)) - 2*Mst1*pow6(Mst2)))) + (576*pow4(Mt)*(-15*
        pow2(Mst1)*pow4(Mst2)*(6*pow2(Mst1)*pow2(Mst2) + 8*pow2(Msq)*(4*pow2(
        Mst1) + pow2(Mst2)) + 19*pow4(Mst1) + 3*pow4(Mst2)) - pow2(Dmglst1)*(-(
        pow4(Msq)*(-604*pow2(Mst1)*pow2(Mst2) + 19008*pow4(Mst1) - 297*pow4(
        Mst2))) + 90*pow2(Mst1)*(4*pow2(Msq) + 11*pow2(Mst1) + pow2(Mst2))*
        pow4(Mst2)) + 180*pow2(Mst2)*(pow2(Mst1) + pow2(Mst2))*pow6(Msq) +
        pow4(Msq)*((-578*pow2(Mst2) + 1320*pow2(MuSUSY))*pow4(Mst1) - 392*pow2(
        Mst1)*pow4(Mst2) + 2383*pow6(Mst1) + 18*pow6(Mst2)) - 2*Dmglst1*(360*
        pow2(Msq)*pow3(Mst1)*pow4(Mst2) + pow4(Msq)*(604*pow2(Mst2)*pow3(Mst1)
        + 297*Mst1*pow4(Mst2) - 5598*pow5(Mst1)) + 390*pow4(Mst2)*pow5(Mst1) +
        90*pow3(Mst1)*pow6(Mst2))) + 8*pow2(Mt)*pow2(s2t)*(-6480*pow2(-(Mst2*
        pow2(Mst1)) + pow3(Mst2))*pow6(Msq) + 270*pow2(MuSUSY)*(5*pow4(Mst2)*
        pow6(Mst1) - 4*pow2(Msq)*(-3*pow4(Mst1)*pow4(Mst2) + 4*pow2(Mst2)*pow6(
        Mst1))) - 18*pow2(Dmglst1)*(90*pow2(MuSUSY)*pow4(Mst1)*(36*pow2(Msq)*
        pow2(Mst2) - 5*pow4(Mst2)) + pow4(Msq)*((-22332*pow2(Mst2) + 6129*pow2(
        MuSUSY))*pow4(Mst1) + 6148*pow2(Mst1)*pow4(Mst2) - 594*pow6(Mst2))) -
        36*Dmglst1*(30*pow2(MuSUSY)*(28*pow2(Msq)*pow2(Mst2) - 5*pow4(Mst2))*
        pow5(Mst1) + pow4(Msq)*(6148*pow3(Mst1)*pow4(Mst2) - (3692*pow2(Mst2) +
        1613*pow2(MuSUSY))*pow5(Mst1) - 594*Mst1*pow6(Mst2))) + pow4(Msq)*(-3*
        pow4(Mst1)*(1244*pow2(Mst2)*pow2(MuSUSY) + 40285*pow4(Mst2)) - (590*
        pow2(Mst2) + 69349*pow2(MuSUSY))*pow6(Mst1) + 16497*pow2(Mst1)*pow6(
        Mst2) - 648*pow8(Mst2))))/pow2(Mst1))/(pow4(Msq)*pow4(Mst2)))) + (4*Mt*
        MuSUSY*(2*Mt*MuSUSY*(270*s2t*(-(pow2(Mst2)*(pow2(Mst1)*(24*Mst1*Mt - 5*
        s2t*pow2(Mst2)) + 20*Dmglst1*Mst1*(6*Mst1*Mt - s2t*pow2(Mst2)) + 30*
        pow2(Dmglst1)*(8*Mst1*Mt - s2t*pow2(Mst2)))) + 4*pow2(Msq)*(6*pow2(
        Dmglst1)*(56*Mst1*Mt - 9*s2t*pow2(Mst2)) + pow2(Mst2)*(-16*Mst1*Mt - 4*
        s2t*pow2(Mst1) + 3*s2t*pow2(Mst2)) + 4*Dmglst1*(24*Mt*pow2(Mst1) - 12*
        Mt*pow2(Mst2) - 7*Mst1*s2t*pow2(Mst2)))) + (-144*(2061*Dmglst1 - 2407*
        Mst1)*Mt*s2t + 95040*pow2(Mt) - (-58068*Dmglst1*Mst1 + 110322*pow2(
        Dmglst1) + 69349*pow2(Mst1) + 3732*pow2(Mst2))*pow2(s2t))*pow4(Msq))*
        pow4(Mst1) - Cbeta*Sbeta*(270*(-36*Mt*pow2(Mst1)*pow2(s2t) + 96*pow3(
        Mt) + 5*Mst1*pow2(Mst2)*pow3(s2t))*pow4(Mst2)*pow5(Mst1) + 18*pow4(
        Mst1)*(2*Dmglst1*((7448*Mst1*s2t*pow2(Mt) - 3*Mt*(6487*pow2(Mst1) +
        7884*pow2(Mst2))*pow2(s2t) + 66608*pow3(Mt) - 3067*Mst1*pow2(Mst2)*
        pow3(s2t))*pow4(Msq) + 30*(-45*Mt*pow2(Mst1)*pow2(s2t) + 56*pow3(Mt) +
        5*Mst1*pow2(Mst2)*pow3(s2t))*pow4(Mst2) - 60*pow2(Msq)*(12*pow2(Mst1)*(
        -9*Mt*pow2(Mst2)*pow2(s2t) + 20*pow3(Mt)) + (36*Mt + 17*Mst1*s2t)*pow2(
        s2t)*pow4(Mst2))) - pow2(Dmglst1)*(s2t*(85836*Mst1*Mt*s2t + 67112*pow2(
        Mt) + 10809*pow2(Mst2)*pow2(s2t))*pow4(Msq) - 450*(-12*Mst1*Mt + s2t*
        pow2(Mst2))*pow2(s2t)*pow4(Mst2) + 180*pow2(Msq)*(64*Mst1*(-3*Mt*pow2(
        Mst2)*pow2(s2t) + 5*pow3(Mt)) + 19*pow3(s2t)*pow4(Mst2))) + 360*pow2(
        Mst2)*pow3(s2t)*pow6(Msq) - 60*pow2(Msq)*(24*pow3(Mst1)*(-(Mt*pow2(
        Mst2)*pow2(s2t)) + 4*pow3(Mt)) + Mst1*(24*Mt + 7*Mst1*s2t)*pow2(s2t)*
        pow4(Mst2) - 3*pow3(s2t)*pow6(Mst2))) + pow4(Msq)*(s2t*pow4(Mst1)*(
        420024*pow2(Mst2)*pow2(Mt) + pow2(Mst1)*(422384*pow2(Mt) - 65617*pow2(
        Mst2)*pow2(s2t)) + 32880*pow2(s2t)*pow4(Mst2)) + 36*Mt*(22096*pow2(Mt)
        + 3*(1789*pow2(Mst1) + 1052*pow2(Mst2))*pow2(s2t))*pow5(Mst1) - 648*
        pow3(s2t)*pow8(Mst2)))))/(pow2(Mst1)*pow4(Msq)*pow4(Mst2))))/3888.)/
        pow4(Mt)/pow2(Sbeta)*12.;
 
   return result;
}

/**
 * 	@return returns the susy log^0 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double himalaya::H5::calc_coef_at_as2_no_sm_logs_log0(){

   const double result =
      (((-4*(-108*Mst1*s2t*pow2(Sbeta)*(240*Dmglst1*s2t*z4*pow2(Mt) + 2*Mt*pow2(
        s2t)*(30*pow2(Dmglst1)*(41*z4 + 8*pow2(z2)) + pow2(Mst2)*(-241*z4 + 32*
        pow2(z2))) + 32*(-z4 + 2*pow2(z2))*pow3(Mt) - Dmglst1*pow2(Mst2)*(265*
        z4 + 4*pow2(z2))*pow3(s2t))*pow4(Mst2) + 27*s2t*pow2(Sbeta)*pow4(Mst2)*
        (-32*Dmglst1*pow2(Mt)*(-4*Mt*z4 + 15*Dmglst1*s2t*z4 + 8*Mt*pow2(z2)) +
        2*s2t*pow2(Mst2)*(4*Dmglst1*Mt*s2t*(241*z4 - 32*pow2(z2)) + pow2(
        Dmglst1)*pow2(s2t)*(265*z4 + 4*pow2(z2)) + pow2(Mt)*(-8*T1ep + 36*z4 +
        10*pow2(z2))) + (4*T1ep - 242*z4 + 35*pow2(z2))*pow3(s2t)*pow4(Mst2)) +
        s2t*pow4(Mst1)*(-664*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*pow2(s2t)*(4*T1ep
        + 2*z4 + 3*pow2(z2)) - 8*s2t*pow2(Mt)*(-(pow2(Mst2)*pow2(Sbeta)*(4*T1ep
        + 2*z4 + 3*pow2(z2))) - 324*Cbeta*Dmglst1*MuSUSY*Sbeta*(187*z4 + 76*
        pow2(z2)) + 2*pow2(MuSUSY)*(-1 + pow2(Sbeta))*(740*T1ep + 127*z4 +
        1527*pow2(z2))) + 320*Cbeta*MuSUSY*Sbeta*(4*T1ep + 2*z4 + 3*pow2(z2))*
        pow3(Mt) - pow2(Sbeta)*(44*T1ep + 6124*z4 + 1113*pow2(z2))*pow3(s2t)*
        pow4(Mst2)) + 108*pow3(Mst1)*(4*MuSUSY*pow2(Mt)*pow2(s2t)*(108*Cbeta*
        Sbeta*pow2(Mst2)*(-z4 + 2*pow2(z2)) - 2*Dmglst1*MuSUSY*(-1 + pow2(
        Sbeta))*(265*z4 + 4*pow2(z2)) + 45*Cbeta*Sbeta*pow2(Dmglst1)*(41*z4 +
        8*pow2(z2))) + 8*MuSUSY*s2t*(60*Cbeta*Dmglst1*Sbeta*z4 + MuSUSY*(-1 +
        pow2(Sbeta))*(-313*z4 + 176*pow2(z2)))*pow3(Mt) - 2*Mt*pow2(Sbeta)*(
        169*z4 + 112*pow2(z2))*pow3(s2t)*pow4(Mst2) + 32*Cbeta*MuSUSY*Sbeta*(-
        z4 + 2*pow2(z2))*pow4(Mt) - Dmglst1*pow2(Sbeta)*(265*z4 + 4*pow2(z2))*
        pow4(Mst2)*pow4(s2t)) + 46656*Cbeta*MuSUSY*Sbeta*pow2(Mt)*pow2(s2t)*(-
        z4 + 2*pow2(z2))*pow5(Mst1) + 6*pow2(Mst1)*(72*MuSUSY*pow2(Mt)*(60*
        MuSUSY*z4*pow2(Mt)*(-1 + pow2(Sbeta)) + s2t*pow2(Dmglst1)*(60*Cbeta*Mt*
        Sbeta*z4 - MuSUSY*s2t*(-1 + pow2(Sbeta))*(25*z4 + 4*pow2(z2))) + 2*
        Dmglst1*Mt*(4*Cbeta*Mt*Sbeta*(-z4 + 2*pow2(z2)) + 3*MuSUSY*s2t*(-1 +
        pow2(Sbeta))*(169*z4 + 112*pow2(z2)))) + 8*Mt*MuSUSY*s2t*pow2(Mst2)*(
        1080*Cbeta*Sbeta*z4*pow2(Dmglst1)*pow2(s2t) + 7*Cbeta*Sbeta*pow2(Mt)*(
        4*T1ep + 2*z4 + 3*pow2(z2)) + Mt*s2t*(54*Cbeta*Dmglst1*Sbeta*(187*z4 +
        76*pow2(z2)) - MuSUSY*(-1 + pow2(Sbeta))*(136*T1ep - 13*z4 + 426*pow2(
        z2)))) - Sbeta*pow2(s2t)*(2*Cbeta*Mt*MuSUSY*s2t*(236*T1ep + 2152*z4 +
        537*pow2(z2)) + Sbeta*(9*pow2(Dmglst1)*pow2(s2t)*(505*z4 + 4*pow2(z2))
        + 2*pow2(Mt)*(20*T1ep + 190*z4 + 87*pow2(z2)) + 36*Dmglst1*Mt*s2t*(989*
        z4 + 272*pow2(z2))))*pow4(Mst2) + pow2(Sbeta)*(100*T1ep + 2165*z4 +
        111*pow2(z2))*pow4(s2t)*pow6(Mst2))))/pow4(Mst2) + (2*z3*(54*Dmglst1*
        Mst1*pow2(Sbeta)*pow4(Mst2)*(880*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 29488*
        Dmglst1*s2t*pow3(Mt) + 5400*Dmglst1*Mt*pow2(Mst2)*pow3(s2t) - 856*pow4(
        Mt) - 195*pow4(Mst2)*pow4(s2t)) - 27*pow2(Dmglst1)*pow2(Sbeta)*pow4(
        Mst2)*(-12156*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 14440*pow4(Mt) + 1143*
        pow4(Mst2)*pow4(s2t)) - 3*pow4(Mst1)*(16*pow2(Mt)*pow2(s2t)*(-21816*
        Cbeta*Dmglst1*MuSUSY*Sbeta*pow2(Mst2) + 54*pow2(Dmglst1)*(1145*pow2(
        MuSUSY)*(-1 + pow2(Sbeta)) - 597*pow2(Mst2)*pow2(Sbeta)) + pow2(Mst2)*(
        17566*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 323*pow2(Mst2)*pow2(Sbeta))) +
        64*s2t*(4824*Cbeta*MuSUSY*Sbeta*pow2(Dmglst1) - 355*Cbeta*MuSUSY*Sbeta*
        pow2(Mst2) + 54*Dmglst1*(-154*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 169*
        pow2(Mst2)*pow2(Sbeta)))*pow3(Mt) + 16*Mt*Sbeta*pow2(Mst2)*(6228*
        Dmglst1*Sbeta*pow2(Mst2) + Cbeta*MuSUSY*(30159*pow2(Dmglst1) + 7307*
        pow2(Mst2)))*pow3(s2t) + 1152*(-470*Cbeta*Dmglst1*MuSUSY*Sbeta + 46*
        pow2(MuSUSY)*(-1 + pow2(Sbeta)) - (1013*pow2(Dmglst1) + 72*pow2(Mst2))*
        pow2(Sbeta))*pow4(Mt) - (107325*pow2(Dmglst1) + 23324*pow2(Mst2))*pow2(
        Sbeta)*pow4(Mst2)*pow4(s2t)) + 108*pow2(Mst2)*pow2(Sbeta)*pow3(Mst1)*(
        1728*s2t*pow2(Mst2)*pow3(Mt) - 12*pow2(Dmglst1)*(1168*s2t*pow3(Mt) - 9*
        Mt*pow2(Mst2)*pow3(s2t)) - 1064*Mt*pow3(s2t)*pow4(Mst2) + Dmglst1*(-
        2984*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 12352*pow4(Mt) + 405*pow4(Mst2)*
        pow4(s2t))) + 54*(-16*pow2(Mt)*pow2(s2t)*(10905*Cbeta*MuSUSY*Sbeta*
        pow2(Dmglst1) - 1326*Cbeta*MuSUSY*Sbeta*pow2(Mst2) + Dmglst1*(3056*
        pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 398*pow2(Mst2)*pow2(Sbeta))) + 64*
        s2t*(-40*Cbeta*Dmglst1*MuSUSY*Sbeta + 575*pow2(MuSUSY)*(-1 + pow2(
        Sbeta)) - (4222*pow2(Dmglst1) + 183*pow2(Mst2))*pow2(Sbeta))*pow3(Mt) -
        16*Mt*Sbeta*pow2(Mst2)*(1423*Cbeta*Dmglst1*MuSUSY - 3284*Sbeta*pow2(
        Dmglst1) + 309*Sbeta*pow2(Mst2))*pow3(s2t) + 64*Sbeta*(129*Cbeta*MuSUSY
        + 658*Dmglst1*Sbeta)*pow4(Mt) + 5077*Dmglst1*pow2(Sbeta)*pow4(Mst2)*
        pow4(s2t))*pow5(Mst1) + 648*log(pow2(Mst2)/pow2(Mst1))*pow2(Mst1)*(2*
        Mst1*s2t*pow2(Sbeta)*(-16*Dmglst1*s2t*pow2(Mt) + 6*Mt*(-78*pow2(
        Dmglst1) + pow2(Mst2))*pow2(s2t) - 32*pow3(Mt) + 23*Dmglst1*pow2(Mst2)*
        pow3(s2t))*pow4(Mst2) - 3*pow2(s2t)*pow4(Mst1)*(8*MuSUSY*pow2(Mt)*(-
        114*Cbeta*Dmglst1*Sbeta + 7*MuSUSY*(-1 + pow2(Sbeta))) + 7*pow2(s2t)*
        pow2(Sbeta)*pow4(Mst2)) + pow2(Sbeta)*pow4(Mst2)*(-16*pow2(Dmglst1)*
        pow2(Mt)*pow2(s2t) - 64*Dmglst1*s2t*pow3(Mt) + 12*Dmglst1*Mt*pow2(Mst2)
        *pow3(s2t) - 16*pow4(Mt) + 23*pow2(Dmglst1)*pow2(Mst2)*pow4(s2t)) + 2*
        pow3(Mst1)*(4*MuSUSY*pow2(Mt)*pow2(s2t)*(351*Cbeta*Sbeta*pow2(Dmglst1)
        + 108*Cbeta*Sbeta*pow2(Mst2) - 46*Dmglst1*MuSUSY*(-1 + pow2(Sbeta))) +
        8*MuSUSY*s2t*(4*Cbeta*Dmglst1*Sbeta + 69*MuSUSY*(-1 + pow2(Sbeta)))*
        pow3(Mt) - 150*Mt*pow2(Sbeta)*pow3(s2t)*pow4(Mst2) + 32*Cbeta*MuSUSY*
        Sbeta*pow4(Mt) - 23*Dmglst1*pow2(Sbeta)*pow4(Mst2)*pow4(s2t)) + pow2(
        Mst1)*(-8*MuSUSY*pow2(Mt)*pow2(s2t)*(-342*Cbeta*Dmglst1*Sbeta*pow2(
        Mst2) + 7*MuSUSY*pow2(Dmglst1)*(-1 + pow2(Sbeta)) + 21*MuSUSY*pow2(
        Mst2)*(-1 + pow2(Sbeta))) + 16*Dmglst1*MuSUSY*s2t*(2*Cbeta*Dmglst1*
        Sbeta + 225*MuSUSY*(-1 + pow2(Sbeta)))*pow3(Mt) - 4*Mt*Sbeta*pow2(Mst2)
        *(-16*Cbeta*MuSUSY*pow2(Dmglst1) + 21*Cbeta*MuSUSY*pow2(Mst2) + 231*
        Dmglst1*Sbeta*pow2(Mst2))*pow3(s2t) + 32*MuSUSY*(2*Cbeta*Dmglst1*Sbeta
        + MuSUSY*(-1 + pow2(Sbeta)))*pow4(Mt) + 3*(-13*pow2(Dmglst1) + 7*pow2(
        Mst2))*pow2(Sbeta)*pow4(Mst2)*pow4(s2t)) + 864*Cbeta*MuSUSY*Sbeta*pow2(
        Mt)*pow2(s2t)*pow5(Mst1)) + 16*(pow2(Mt)*pow2(s2t)*(-176823*Cbeta*
        Dmglst1*MuSUSY*Sbeta - 135212*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 662*
        pow2(Mst2)*pow2(Sbeta)) + 8*s2t*Sbeta*(367*Cbeta*MuSUSY - 59292*
        Dmglst1*Sbeta)*pow3(Mt) + Mt*Sbeta*(-41257*Cbeta*MuSUSY + 80757*
        Dmglst1*Sbeta)*pow2(Mst2)*pow3(s2t) + 34236*pow2(Sbeta)*pow4(Mt) +
        4834*pow2(Sbeta)*pow4(Mst2)*pow4(s2t))*pow6(Mst1) + 54*pow2(Mst1)*pow2(
        Mst2)*pow2(Sbeta)*(-32*Dmglst1*(-74*s2t*pow2(Mst2)*pow3(Mt) + 29*Mt*
        pow3(s2t)*pow4(Mst2)) + pow2(Dmglst1)*(-26158*pow2(Mst2)*pow2(Mt)*pow2(
        s2t) + 49728*pow4(Mt) + 1311*pow4(Mst2)*pow4(s2t)) + 8*(-43*pow2(Mt)*
        pow2(s2t)*pow4(Mst2) - 652*pow2(Mst2)*pow4(Mt) + 41*pow4(s2t)*pow6(
        Mst2))) + 432*Mt*s2t*Sbeta*(1095*Cbeta*Mt*MuSUSY*s2t - 3584*Sbeta*pow2(
        Mt) + 519*Sbeta*pow2(Mst2)*pow2(s2t))*pow7(Mst1)))/(pow2(Mst1)*pow4(
        Mst2)) - 3888*Cbeta*MuSUSY*Sbeta*(-(Mt*pow3(s2t)*((pow2(Mst1)*((-1328*
        OepS2 + 54*S2*(1187 + 498*log(pow2(Mst2)/pow2(Mst1))))*pow2(Mst1) +
        177*(-8*OepS2 + 81*S2*(5 + 2*log(pow2(Mst2)/pow2(Mst1))))*pow2(Mst2)))/
        (729.*pow2(Mst2)) + ((9600*Dmglst1 + 2361*Mst1 - 25*(254*Dmglst1 + 63*
        Mst1)*log(pow2(Mst2)/pow2(Mst1)) - 60*log(pow2(Msq)/pow2(Mst1))*(-30*
        Dmglst1 - 9*Mst1 + 5*(17*Dmglst1 + 4*Mst1)*log(pow2(Mst2)/pow2(Mst1)))
        + 75*(34*Dmglst1 + 15*Mst1)*pow2(log(pow2(Mst2)/pow2(Mst1))))*pow3(
        Mst1))/(135.*pow2(Msq)) + (pow2(Mst1)*(311855428500*Dmglst1*Mst1 -
        238997855250*pow2(Dmglst1) + 14224896000*B4*pow2(Dmglst1) + 7112448000*
        DN*pow2(Dmglst1) + 6667920000*pow2(Msq) + 478241812219*pow2(Mst1) +
        1260*log(pow2(Mst2)/pow2(Mst1))*(255108700*Dmglst1*Mst1 + 453223050*
        pow2(Dmglst1) + 215618061*pow2(Mst1)) + 10001880000*pow2(Dmglst1 +
        Mst1)*pow2(log(pow2(Msq)/pow2(Mst1))) - 264600*(703780*Dmglst1*Mst1 +
        1566390*pow2(Dmglst1) + 770503*pow2(Mst1))*pow2(log(pow2(Mst2)/pow2(
        Mst1))) + 1666980000*log(pow2(Msq)/pow2(Mst1))*(44*Dmglst1*Mst1 - 66*
        pow2(Dmglst1) + 8*pow2(Msq) + 29*pow2(Mst1) + 6*log(pow2(Mst2)/pow2(
        Mst1))*(16*Dmglst1*Mst1 + 28*pow2(Dmglst1) + 3*pow2(Mst1)) - 16*pow2(
        Mst1)*pow2(log(pow2(Mst2)/pow2(Mst1)))) + 592704000*(190*Dmglst1*Mst1 +
        178*pow2(Dmglst1) + 273*pow2(Mst1))*pow3(log(pow2(Mst2)/pow2(Mst1)))))/
        (2.000376e9*pow2(Mst2)) - (pow2(Mst1)*(-27450000*pow2(Dmglst1)*pow2(
        Msq) + 12009375*pow2(Dmglst1)*pow2(Mst2) + 3195000*pow2(Msq)*pow2(Mst2)
        + 90*log(pow2(Mst2)/pow2(Mst1))*(375*pow2(Dmglst1)*(596*pow2(Msq) -
        179*pow2(Mst2)) - 8140*pow2(Msq)*pow2(Mst2) - 393289*pow4(Msq)) - 1350*
        pow2(log(pow2(Mst2)/pow2(Mst1)))*(150*pow2(Dmglst1)*(38*pow2(Msq) - 5*
        pow2(Mst2)) + 380*pow2(Msq)*pow2(Mst2) - 24453*pow4(Msq)) - 26775829*
        pow4(Msq) + 1080000*B4*pow4(Msq) + 54000*DN*pow4(Msq) - 405000*(2 + 7*
        log(pow2(Mst2)/pow2(Mst1)))*pow2(log(pow2(Msq)/pow2(Mst1)))*pow4(Msq) -
        12748500*pow3(log(pow2(Mst2)/pow2(Mst1)))*pow4(Msq) + 13500*log(pow2(
        Msq)/pow2(Mst1))*(390*pow2(log(pow2(Mst2)/pow2(Mst1)))*pow4(Msq) - 5*(-
        4*pow2(Msq)*pow2(Mst2) + 3*pow2(Dmglst1)*(8*pow2(Msq) + 15*pow2(Mst2))
        + 32*pow4(Msq)) - 2*log(pow2(Mst2)/pow2(Mst1))*(2*pow2(Msq)*pow2(Mst2)
        + pow2(Dmglst1)*(-570*pow2(Msq) + 75*pow2(Mst2)) + 295*pow4(Msq)))))/(
        243000.*pow4(Msq)) + (pow2(Mst2)*(1543500*pow2(Dmglst1)*pow2(Mst2) -
        86264500*Dmglst1*pow3(Mst1) + 3430*(16700*Dmglst1 + 2469*Mst1)*log(
        pow2(Mst2)/pow2(Mst1))*pow3(Mst1) + 420*log(pow2(Msq)/pow2(Mst1))*(
        4900*Dmglst1 - 12899*Mst1 + 490*(100*Dmglst1 + 11*Mst1)*log(pow2(Mst2)/
        pow2(Mst1)))*pow3(Mst1) - 514500*(20*Dmglst1 - 3*Mst1)*pow2(log(pow2(
        Mst2)/pow2(Mst1)))*pow3(Mst1) - 20222907*pow4(Mst1) + 176400*pow2(log(
        pow2(Msq)/pow2(Mst1)))*pow4(Mst1)))/(3.7044e6*pow4(Msq)) + (-(1 + 2*
        log(pow2(Mst2)/pow2(Mst1)))/(3.*pow2(Mst1)) + (939 - 760*log(pow2(Mst2)
        /pow2(Mst1)) - 30*log(pow2(Msq)/pow2(Mst1))*(-12 + 5*log(pow2(Mst2)/
        pow2(Mst1))) + 150*pow2(log(pow2(Mst2)/pow2(Mst1))))/(135.*pow2(Msq)) +
        ((5*Dmglst1*Mst1)/6. - (pow2(Mst1)*(-4016157 + 6638660*log(pow2(Mst2)/
        pow2(Mst1)) + 420*log(pow2(Msq)/pow2(Mst1))*(-5549 + 6020*log(pow2(
        Mst2)/pow2(Mst1))) + 176400*pow2(log(pow2(Msq)/pow2(Mst1))) - 2704800*
        pow2(log(pow2(Mst2)/pow2(Mst1)))))/3.7044e6)/pow4(Msq))*pow4(Mst2))) -
        (pow2(Mst1)*pow3(Mt)*(-176400*Mst1*Mt*pow2(Dmglst1)*pow2(Msq) - 105360*
        Dmglst1*Mt*pow2(Msq)*pow2(Mst1) - 140940*Mst1*Mt*pow2(Dmglst1)*pow2(
        Mst2) - 105360*Dmglst1*Mt*pow2(Msq)*pow2(Mst2) - 23760*Mst1*Mt*pow2(
        Msq)*pow2(Mst2) - 52260*Dmglst1*Mt*pow2(Mst1)*pow2(Mst2) - 23760*Mt*
        pow2(Msq)*pow3(Mst1) - 8100*Mt*pow2(Mst2)*pow3(Mst1) + 72*Mt*pow2(Msq)*
        pow2(log(pow2(Mst2)/pow2(Mst1)))*(1800*Mst1*pow2(Dmglst1) - 3037*
        Dmglst1*pow2(Msq) - 879*Mst1*pow2(Msq) + 900*Dmglst1*pow2(Mst1) + 180*
        pow3(Mst1)) - 65212*Dmglst1*Mt*pow4(Msq) + 20964*Mst1*Mt*pow4(Msq) +
        72*(831*Dmglst1 + 265*Mst1)*Mt*pow3(log(pow2(Mst2)/pow2(Mst1)))*pow4(
        Msq) + 21135*Dmglst1*Mt*pow4(Mst2) + 10035*Mst1*Mt*pow4(Mst2) + 12825*
        s2t*pow2(Dmglst1)*pow4(Mst2) - 1080*Mt*pow2(log(pow2(Msq)/pow2(Mst1)))*
        (6*(Dmglst1 + Mst1)*log(pow2(Mst2)/pow2(Mst1))*pow4(Msq) + (7*Dmglst1 +
        3*Mst1)*pow4(Mst2)) - 6*log(pow2(Mst2)/pow2(Mst1))*(-3*Mst1*Mt*(120*
        pow2(Msq)*(5*pow2(Mst1) + 7*pow2(Mst2)) - 15*pow2(Mst2)*(2*pow2(Mst1) +
        15*pow2(Mst2)) + 1286*pow4(Msq)) + Dmglst1*Mt*(120*pow2(Msq)*(pow2(
        Mst1) - 65*pow2(Mst2)) - 1590*pow2(Mst1)*pow2(Mst2) - 31778*pow4(Msq) +
        1335*pow4(Mst2)) + 90*pow2(Dmglst1)*(208*Mst1*Mt*pow2(Msq) - 54*Mst1*
        Mt*pow2(Mst2) + 15*s2t*pow4(Mst2))) + 180*log(pow2(Msq)/pow2(Mst1))*(-
        144*(3*Dmglst1 + Mst1)*Mt*pow2(log(pow2(Mst2)/pow2(Mst1)))*pow4(Msq) +
        Dmglst1*Mt*(604*pow4(Msq) - pow4(Mst2)) - 162*s2t*pow2(Dmglst1)*pow4(
        Mst2) + 3*Mst1*Mt*(84*pow4(Msq) + pow4(Mst2)) + 6*Mt*log(pow2(Mst2)/
        pow2(Mst1))*(12*Mst1*pow2(Dmglst1)*(2*pow2(Msq) + 3*pow2(Mst2)) + 3*
        Mst1*(2*pow2(Mst1)*pow2(Mst2) + 4*pow2(Msq)*(pow2(Mst1) + pow2(Mst2)) -
        2*pow4(Msq) + pow4(Mst2)) + Dmglst1*(22*pow2(Mst1)*pow2(Mst2) + 28*
        pow2(Msq)*(pow2(Mst1) + pow2(Mst2)) + 66*pow4(Msq) + 7*pow4(Mst2))))))/
        (243.*pow4(Msq)*pow4(Mst2)) + (s2t*pow3(Mt)*((343000*pow2(Mst1)*((-160*
        OepS2 + 54*S2*(139 + 60*log(pow2(Mst2)/pow2(Mst1))))*pow2(Mst1) + 21*(-
        8*OepS2 + 81*S2*(5 + 2*log(pow2(Mst2)/pow2(Mst1))))*pow2(Mst2)))/pow4(
        Mst2) + (675*(-334425*pow2(Mst1)*pow2(Mst2) - 377300*Dmglst1*pow3(Mst1)
        + 187967*pow4(Mst1) + 28*log(pow2(Mst2)/pow2(Mst1))*(245*pow2(Mst1)*
        pow2(Mst2) + 73500*Dmglst1*pow3(Mst1) + 16856*pow4(Mst1) - 5359*pow4(
        Mst2)) - 44100*pow2(log(pow2(Msq)/pow2(Mst1)))*(pow4(Mst1) - pow4(Mst2)
        ) + 103583*pow4(Mst2) + 5880*pow2(log(pow2(Mst2)/pow2(Mst1)))*(35*pow2(
        Mst1)*pow2(Mst2) + 35*pow4(Mst1) + 11*pow4(Mst2)) + 420*log(pow2(Msq)/
        pow2(Mst1))*(-735*pow2(Mst1)*pow2(Mst2) + 11760*Dmglst1*pow3(Mst1) +
        2194*pow4(Mst1) + 256*pow4(Mst2) - 7*log(pow2(Mst2)/pow2(Mst1))*(70*
        pow2(Mst1)*pow2(Mst2) + 7*pow4(Mst1) + 37*pow4(Mst2)))))/pow4(Msq) - (
        2*pow2(Mst1)*(-27448060500*Dmglst1*Mst1*pow2(Msq) + 222264000*B4*
        Dmglst1*Mst1*pow2(Msq) + 111132000*Dmglst1*DN*Mst1*pow2(Msq) -
        58393434750*pow2(Dmglst1)*pow2(Msq) + 111132000*B4*pow2(Dmglst1)*pow2(
        Msq) + 55566000*DN*pow2(Dmglst1)*pow2(Msq) - 1311202751*pow2(Msq)*pow2(
        Mst1) + 5093550000*Dmglst1*Mst1*pow2(Mst2) + 7640325000*pow2(Dmglst1)*
        pow2(Mst2) - 16204796958*pow2(Msq)*pow2(Mst2) + 1273387500*pow2(Mst1)*
        pow2(Mst2) + 66150*(187110*Dmglst1*Mst1*pow2(Msq) + 751555*pow2(
        Dmglst1)*pow2(Msq) - 3150*pow2(Mst2)*(pow2(Mst1) + pow2(Mst2)) - pow2(
        Msq)*(25412*pow2(Mst1) + 61551*pow2(Mst2)))*pow2(log(pow2(Mst2)/pow2(
        Mst1))) - 4630500*pow2(Msq)*(1896*Dmglst1*Mst1 + 3160*pow2(Dmglst1) +
        763*pow2(Mst1) + 317*pow2(Mst2))*pow3(log(pow2(Mst2)/pow2(Mst1))) -
        13891500*pow2(log(pow2(Msq)/pow2(Mst1)))*(45*pow2(Msq)*pow2(Mst2) - 4*
        pow4(Mst2)) - 124529580*pow4(Mst2) + 926100*log(pow2(Msq)/pow2(Mst1))*(
        -900*pow2(Msq)*pow2(Mst1) - 300*Dmglst1*Mst1*(27*pow2(Msq) - 10*pow2(
        Mst2)) - 450*pow2(Dmglst1)*(47*pow2(Msq) - 10*pow2(Mst2)) - 2025*pow2(
        Msq)*pow2(Mst2) + 750*pow2(Mst1)*pow2(Mst2) + 450*pow2(Msq)*(pow2(Mst1)
        + pow2(Mst2))*pow2(log(pow2(Mst2)/pow2(Mst1))) - 532*pow4(Mst2) + 75*
        log(pow2(Mst2)/pow2(Mst1))*(120*Dmglst1*Mst1*pow2(Msq) + 180*pow2(
        Dmglst1)*pow2(Msq) + 6*pow2(Msq)*(5*pow2(Mst1) + 3*pow2(Mst2)) + pow4(
        Mst2))) + 630*log(pow2(Mst2)/pow2(Mst1))*(9800*pow2(Dmglst1)*(4838*
        pow2(Msq) - 675*pow2(Mst2)) + 147000*Dmglst1*Mst1*(1061*pow2(Msq) - 30*
        pow2(Mst2)) + 3*(pow2(Msq)*(23055979*pow2(Mst1) + 8509732*pow2(Mst2)) -
        122500*(3*pow2(Mst1)*pow2(Mst2) + pow4(Mst2))))))/(pow2(Msq)*pow4(Mst2)
        )))/6.251175e7 - (pow2(Mt)*pow2(s2t)*(-812160*pow2(Dmglst1)*pow2(Msq)*
        pow2(Mst2)*pow3(Mst1) + 481008*Dmglst1*pow2(Mst1)*pow2(Mst2)*pow4(Msq)
        + 255744*B4*Dmglst1*pow2(Mst1)*pow2(Mst2)*pow4(Msq) - 3456*Dmglst1*DN*
        pow2(Mst1)*pow2(Mst2)*pow4(Msq) + 51840*(Dmglst1 + Mst1)*log(pow2(Mst2)
        /pow2(Mst1))*pow2(Mst1)*(pow2(Mst1) + pow2(Mst2))*pow2(log(pow2(Msq)/
        pow2(Mst1)))*pow4(Msq) + 7167914*pow2(Dmglst1)*pow3(Mst1)*pow4(Msq) +
        196992*B4*pow2(Dmglst1)*pow3(Mst1)*pow4(Msq) - 5184*DN*pow2(Dmglst1)*
        pow3(Mst1)*pow4(Msq) + 100432*pow2(Mst2)*pow3(Mst1)*pow4(Msq) + 124416*
        B4*pow2(Mst2)*pow3(Mst1)*pow4(Msq) + 288*pow2(Mst1)*(-1317*Mst1*pow2(
        Dmglst1) - 889*Dmglst1*pow2(Mst1) - 566*Dmglst1*pow2(Mst2) + 68*Mst1*
        pow2(Mst2) + 69*pow3(Mst1))*pow3(log(pow2(Mst2)/pow2(Mst1)))*pow4(Msq)
        - 207360*Dmglst1*pow2(Msq)*pow2(Mst2)*pow4(Mst1) + 3539195*Dmglst1*
        pow4(Msq)*pow4(Mst1) + 255744*B4*Dmglst1*pow4(Msq)*pow4(Mst1) - 3456*
        Dmglst1*DN*pow4(Msq)*pow4(Mst1) + 717120*Dmglst1*pow2(Msq)*pow2(Mst1)*
        pow4(Mst2) + 720720*pow2(Dmglst1)*pow3(Mst1)*pow4(Mst2) + 204480*pow2(
        Msq)*pow3(Mst1)*pow4(Mst2) + 309540*Dmglst1*pow4(Mst1)*pow4(Mst2) + 72*
        pow2(Mst1)*pow2(log(pow2(Mst2)/pow2(Mst1)))*(-240*Mst1*pow2(Msq)*(pow2(
        Mst1) - pow2(Mst2))*pow2(Mst2) - (6500*Mst1*pow2(Mst2) + 9079*pow3(
        Mst1))*pow4(Msq) + 90*pow3(Mst1)*pow4(Mst2) + 2*Mst1*pow2(Dmglst1)*(-
        2880*pow2(Msq)*pow2(Mst2) + 17659*pow4(Msq) + 450*pow4(Mst2)) +
        Dmglst1*((6853*pow2(Mst1) + 1284*pow2(Mst2))*pow4(Msq) - 720*pow2(Msq)*
        (3*pow2(Mst1)*pow2(Mst2) - pow4(Mst2)) + 450*pow2(Mst1)*pow4(Mst2))) -
        34560*pow2(Msq)*pow2(Mst2)*pow5(Mst1) + 722583*pow4(Msq)*pow5(Mst1) +
        124416*B4*pow4(Msq)*pow5(Mst1) + 54660*pow4(Mst2)*pow5(Mst1) - 15720*
        Dmglst1*pow2(Mst1)*pow6(Mst2) - 7080*pow3(Mst1)*pow6(Mst2) + 660*
        Dmglst1*pow8(Mst2) + 660*Mst1*pow8(Mst2) + 720*log(pow2(Msq)/pow2(Mst1)
        )*(72*pow2(Mst1)*(6*Mst1*pow2(Dmglst1) + 5*Dmglst1*(pow2(Mst1) + pow2(
        Mst2)) + Mst1*(pow2(Mst1) + pow2(Mst2)))*pow2(log(pow2(Mst2)/pow2(Mst1)
        ))*pow4(Msq) + 72*pow2(Mst2)*pow3(Mst1)*pow4(Msq) + 96*pow2(Msq)*pow3(
        Mst1)*pow4(Mst2) + 72*pow2(Dmglst1)*pow3(Mst1)*(19*pow4(Msq) + 5*pow4(
        Mst2)) - 6*log(pow2(Mst2)/pow2(Mst1))*pow2(Mst1)*(8*Mst1*pow2(Msq)*
        pow2(Mst2)*(2*pow2(Mst1) + pow2(Mst2)) - 30*(2*Mst1*pow2(Mst2) + pow3(
        Mst1))*pow4(Msq) + 3*pow3(Mst1)*pow4(Mst2) + 6*Mst1*pow2(Dmglst1)*(8*
        pow2(Msq)*pow2(Mst2) + 118*pow4(Msq) + 5*pow4(Mst2)) + 3*Dmglst1*(6*(
        13*pow2(Mst1) + 2*pow2(Mst2))*pow4(Msq) + 5*pow2(Mst1)*pow4(Mst2) + 8*
        pow2(Msq)*(2*pow2(Mst1)*pow2(Mst2) + pow4(Mst2)))) - 126*pow4(Msq)*
        pow5(Mst1) + 37*pow4(Mst2)*pow5(Mst1) - 2*pow3(Mst1)*pow6(Mst2) + Mst1*
        pow8(Mst2) + Dmglst1*(pow4(Msq)*(504*pow2(Mst1)*pow2(Mst2) - 270*pow4(
        Mst1)) + 288*pow2(Msq)*pow2(Mst1)*pow4(Mst2) + 181*pow4(Mst1)*pow4(
        Mst2) - 2*pow2(Mst1)*pow6(Mst2) + pow8(Mst2))) + 12*log(pow2(Mst2)/
        pow2(Mst1))*(2*pow2(Dmglst1)*pow3(Mst1)*(7200*pow2(Msq)*pow2(Mst2) +
        691*pow4(Msq) - 15570*pow4(Mst2)) + 3*Dmglst1*(pow4(Msq)*(4780*pow2(
        Mst1)*pow2(Mst2) + 31467*pow4(Mst1)) - 10*pow4(Mst2)*(-4*pow2(Mst1)*
        pow2(Mst2) + 479*pow4(Mst1) + 2*pow4(Mst2)) - 960*pow2(Msq)*(13*pow2(
        Mst2)*pow4(Mst1) + 11*pow2(Mst1)*pow4(Mst2))) + pow4(Msq)*(90196*pow2(
        Mst2)*pow3(Mst1) + 101509*pow5(Mst1)) - 960*pow2(Msq)*(10*pow3(Mst1)*
        pow4(Mst2) + 17*pow2(Mst2)*pow5(Mst1)) - 30*(83*pow4(Mst2)*pow5(Mst1) -
        4*pow3(Mst1)*pow6(Mst2) + 2*Mst1*pow8(Mst2)))))/(1296.*pow4(Msq)*pow4(
        Mst2))) + z2*(-(pow2(Sbeta)*(144*Mt*pow3(s2t)*(12*Mst1*(1187 + 360*log(
        pow2(Msq)/pow2(Mst1)) - 747*log(pow2(Mst2)/pow2(Mst1)))*pow2(Dmglst1) -
        (2*(1440*Dmglst1*Mst1 + 3240*pow2(Dmglst1) + (-1355 - 540*log(pow2(Msq)
        /pow2(Mst1)) + 744*log(pow2(Mst2)/pow2(Mst1)))*pow2(Msq) + 240*pow2(
        Mst1))*pow3(Mst1))/pow2(Msq) - ((1397*Dmglst1*Mst1 - 62*pow2(Dmglst1) +
        737*pow2(Mst1) + 12*log(pow2(Mst2)/pow2(Mst1))*(427*Dmglst1*Mst1 + 702*
        pow2(Dmglst1) + 71*pow2(Mst1)))*pow3(Mst1))/pow2(Mst2) + 6*(Dmglst1*(
        1941 + 660*log(pow2(Msq)/pow2(Mst1)) - 1244*log(pow2(Mst2)/pow2(Mst1)))
        *pow2(Mst1) + pow2(Mst2)*((Dmglst1 + Mst1)*(-627 - 60*log(pow2(Msq)/
        pow2(Mst1)) + 192*log(pow2(Mst2)/pow2(Mst1)) + (120*Dmglst1*Mst1)/pow2(
        Msq)) + (5*(15*Dmglst1*Mst1 + 30*pow2(Dmglst1) + 8*pow2(Msq) + 3*pow2(
        Mst1))*pow3(Mst1))/pow4(Msq)))) + ((pow2(Mst2)*pow4(s2t)*(1080*pow2(
        Msq)*pow2(Mst1)*(pow2(Dmglst1)*(60*pow2(Mst1) - 3*pow2(Mst2)) - 3*pow2(
        Mst1)*pow2(Mst2) + Dmglst1*(-6*Mst1*pow2(Mst2) + 40*pow3(Mst1)) + 10*
        pow4(Mst1))*pow4(Mst2) + 6480*(pow2(Mst1) - pow2(Mst2))*(2*log(pow2(
        Mst2)/pow2(Mst1))*pow2(Mst1)*pow2(Mst2) + pow4(Mst1) - pow4(Mst2))*
        pow6(Msq) + pow2(Mst2)*pow4(Msq)*(3240*log(pow2(Msq)/pow2(Mst1))*pow2(
        Mst1)*(pow2(Mst1) - pow2(Mst2))*(12*Dmglst1*Mst1 + 6*pow2(Dmglst1) + 5*
        pow2(Mst1) - pow2(Mst2)) - 62856*pow2(Dmglst1)*pow2(Mst1)*pow2(Mst2) -
        125712*Dmglst1*pow2(Mst2)*pow3(Mst1) + 268110*pow2(Dmglst1)*pow4(Mst1)
        - 68844*pow2(Mst2)*pow4(Mst1) + 108*log(pow2(Mst2)/pow2(Mst1))*pow2(
        Mst1)*(670*pow2(Mst1)*pow2(Mst2) + pow2(Dmglst1)*(-310*pow2(Mst1) +
        630*pow2(Mst2)) - 84*Dmglst1*(-15*Mst1*pow2(Mst2) + 7*pow3(Mst1)) +
        103*pow4(Mst1) - 149*pow4(Mst2)) - 21384*Dmglst1*Mst1*pow4(Mst2) -
        10692*pow2(Dmglst1)*pow4(Mst2) + 35316*pow2(Mst1)*pow4(Mst2) + 257508*
        Dmglst1*pow5(Mst1) + 98497*pow6(Mst1)) - 1350*(4*Dmglst1*Mst1 + 6*pow2(
        Dmglst1) + pow2(Mst1))*pow4(Mst1)*pow6(Mst2)))/pow2(Mst1) + (576*pow4(
        Mt)*(-15*pow2(Mst1)*pow4(Mst2)*(6*pow2(Mst1)*pow2(Mst2) + 8*pow2(Msq)*(
        4*pow2(Mst1) + pow2(Mst2)) + 19*pow4(Mst1) + 3*pow4(Mst2)) - pow2(
        Dmglst1)*(-(pow4(Msq)*(4*(-151 + 108*log(pow2(Mst2)/pow2(Mst1)))*pow2(
        Mst1)*pow2(Mst2) + 24*(792 + 225*log(pow2(Msq)/pow2(Mst1)) - 368*log(
        pow2(Mst2)/pow2(Mst1)))*pow4(Mst1) - 297*pow4(Mst2))) + 90*pow2(Mst1)*(
        4*pow2(Msq) + 11*pow2(Mst1) + pow2(Mst2))*pow4(Mst2)) + 180*pow2(Mst2)*
        (pow2(Mst1) + pow2(Mst2))*pow6(Msq) + pow4(Msq)*(-578*pow2(Mst2)*pow4(
        Mst1) + 1320*pow2(MuSUSY)*pow4(Mst1) - 6*log(pow2(Mst2)/pow2(Mst1))*(
        205*pow2(Mst1) - 76*pow2(Mst2) + 64*pow2(MuSUSY))*pow4(Mst1) - 392*
        pow2(Mst1)*pow4(Mst2) + 2383*pow6(Mst1) + 180*log(pow2(Msq)/pow2(Mst1))
        *(2*pow2(Mst1)*pow4(Mst2) + 5*pow6(Mst1)) + 18*pow6(Mst2)) - 2*Dmglst1*
        (360*pow2(Msq)*pow3(Mst1)*pow4(Mst2) + 390*pow4(Mst2)*pow5(Mst1) +
        pow4(Msq)*(4*(151 - 108*log(pow2(Mst2)/pow2(Mst1)))*pow2(Mst2)*pow3(
        Mst1) + 297*Mst1*pow4(Mst2) - 18*(311 + 100*log(pow2(Msq)/pow2(Mst1)) -
        154*log(pow2(Mst2)/pow2(Mst1)))*pow5(Mst1)) + 90*pow3(Mst1)*pow6(Mst2))
        ) + 8*pow2(Mt)*pow2(s2t)*(-6480*(pow2(-(Mst2*pow2(Mst1)) + pow3(Mst2))
        - 2*log(pow2(Mst2)/pow2(Mst1))*pow2(MuSUSY)*pow4(Mst1))*pow6(Msq) +
        270*pow2(MuSUSY)*(5*pow4(Mst2)*pow6(Mst1) - 4*pow2(Msq)*(-3*pow4(Mst1)*
        pow4(Mst2) + 4*pow2(Mst2)*pow6(Mst1))) - 18*pow2(Dmglst1)*(90*pow2(
        MuSUSY)*pow4(Mst1)*(36*pow2(Msq)*pow2(Mst2) - 5*pow4(Mst2)) + pow4(Msq)
        *(3*(-7444*pow2(Mst2) + 2043*pow2(MuSUSY) - 360*log(pow2(Msq)/pow2(
        Mst1))*pow2(MuSUSY) + 4*log(pow2(Mst2)/pow2(Mst1))*(1866*pow2(Mst2) +
        475*pow2(MuSUSY)))*pow4(Mst1) + 4*(1537 - 228*log(pow2(Mst2)/pow2(Mst1)
        ))*pow2(Mst1)*pow4(Mst2) - 594*pow6(Mst2))) - 36*Dmglst1*(30*pow2(
        MuSUSY)*(28*pow2(Msq)*pow2(Mst2) - 5*pow4(Mst2))*pow5(Mst1) + pow4(Msq)
        *(4*(1537 - 228*log(pow2(Mst2)/pow2(Mst1)))*pow3(Mst1)*pow4(Mst2) + (-
        3692*pow2(Mst2) - 1613*pow2(MuSUSY) - 1080*log(pow2(Msq)/pow2(Mst1))*
        pow2(MuSUSY) + 84*log(pow2(Mst2)/pow2(Mst1))*(70*pow2(Mst2) + 69*pow2(
        MuSUSY)))*pow5(Mst1) - 594*Mst1*pow6(Mst2))) + pow4(Msq)*(-3732*pow2(
        Mst2)*pow2(MuSUSY)*pow4(Mst1) + 12960*log(pow2(Msq)/pow2(Mst1))*(pow2(
        Mst1) + pow2(Mst2))*pow2(MuSUSY)*pow4(Mst1) - 120855*pow4(Mst1)*pow4(
        Mst2) - 590*pow2(Mst2)*pow6(Mst1) - 69349*pow2(MuSUSY)*pow6(Mst1) +
        16497*pow2(Mst1)*pow6(Mst2) - 216*log(pow2(Mst2)/pow2(Mst1))*(pow4(
        Mst1)*(186*pow2(Mst2)*pow2(MuSUSY) - 40*pow4(Mst2)) + 3*(55*pow2(Mst2)
        + 166*pow2(MuSUSY))*pow6(Mst1) + 3*pow2(Mst1)*pow6(Mst2)) - 648*pow8(
        Mst2))))/pow2(Mst1) + 1152*s2t*pow3(Mt)*(Mst1*pow2(Dmglst1)*(360*pow2(
        Msq)*pow2(Mst1)*(10*pow2(Mst2) + 7*pow2(MuSUSY)) + (-42755*pow2(Mst1) -
        8283*pow2(Mst2) - 2160*log(pow2(Msq)/pow2(Mst1))*(5*pow2(Mst1) + pow2(
        Mst2)) + 6*log(pow2(Mst2)/pow2(Mst1))*(3389*pow2(Mst1) + 813*pow2(Mst2)
        ))*pow4(Msq) + 90*(pow2(Mst1)*(-5*pow2(Mst2)*pow2(MuSUSY) + 2*pow4(
        Mst2)) - 2*pow6(Mst2))) + Dmglst1*(pow4(Msq)*(-8300*pow2(Mst1)*pow2(
        Mst2) - 2061*pow2(Mst1)*pow2(MuSUSY) - 540*log(pow2(Msq)/pow2(Mst1))*
        pow2(Mst1)*(10*pow2(Mst1) + 4*pow2(Mst2) + 3*pow2(MuSUSY)) - 22119*
        pow4(Mst1) + 6*log(pow2(Mst2)/pow2(Mst1))*(5*pow2(Mst1)*(165*pow2(Mst2)
        + 86*pow2(MuSUSY)) + 1772*pow4(Mst1) - 18*pow4(Mst2)) - 26*pow4(Mst2))
        + 120*pow2(Msq)*(3*(5*pow2(Mst2) + 2*pow2(MuSUSY))*pow4(Mst1) + pow2(
        Mst1)*(-3*pow2(Mst2)*pow2(MuSUSY) + pow4(Mst2)) - pow6(Mst2)) + 15*(
        pow4(Mst1)*(-15*pow2(Mst2)*pow2(MuSUSY) + 14*pow4(Mst2)) - 12*pow2(
        Mst1)*pow6(Mst2) - 2*pow8(Mst2))) + Mst1*(pow4(Msq)*(-2778*pow2(Mst1)*
        pow2(Mst2) + 2407*pow2(Mst1)*pow2(MuSUSY) - 180*log(pow2(Msq)/pow2(
        Mst1))*pow2(Mst1)*(6*pow2(Mst1) + 4*pow2(Mst2) + pow2(MuSUSY)) - 4589*
        pow4(Mst1) + 6*log(pow2(Mst2)/pow2(Mst1))*(pow2(Mst1)*(283*pow2(Mst2) -
        68*pow2(MuSUSY)) + 366*pow4(Mst1) - 18*pow4(Mst2)) - 26*pow4(Mst2)) -
        120*pow2(Msq)*pow2(Mst2)*(pow2(Mst1)*(-pow2(Mst2) + pow2(MuSUSY)) - 3*
        pow4(Mst1) + pow4(Mst2)) + 15*(pow4(Mst1)*(-3*pow2(Mst2)*pow2(MuSUSY) +
        6*pow4(Mst2)) - 4*pow2(Mst1)*pow6(Mst2) - 2*pow8(Mst2)))))/(pow4(Msq)*
        pow4(Mst2)))) + (4*Mt*MuSUSY*(2*Mt*MuSUSY*pow4(Mst1)*(270*s2t*(-(pow2(
        Mst2)*(pow2(Mst1)*(24*Mst1*Mt - 5*s2t*pow2(Mst2)) + 20*Dmglst1*Mst1*(6*
        Mst1*Mt - s2t*pow2(Mst2)) + 30*pow2(Dmglst1)*(8*Mst1*Mt - s2t*pow2(
        Mst2)))) + 4*pow2(Msq)*(6*pow2(Dmglst1)*(56*Mst1*Mt - 9*s2t*pow2(Mst2))
        + pow2(Mst2)*(-16*Mst1*Mt - 4*s2t*pow2(Mst1) + 3*s2t*pow2(Mst2)) + 4*
        Dmglst1*(24*Mt*pow2(Mst1) - 12*Mt*pow2(Mst2) - 7*Mst1*s2t*pow2(Mst2))))
        + (-296784*Dmglst1*Mt*s2t + 346608*Mst1*Mt*s2t + 6480*s2t*log(pow2(Msq)
        /pow2(Mst1))*(-36*Dmglst1*Mt - 4*Mst1*Mt + 6*Dmglst1*Mst1*s2t + 3*s2t*
        pow2(Dmglst1) + 2*s2t*pow2(Mst1) + 2*s2t*pow2(Mst2)) + 95040*pow2(Mt) +
        58068*Dmglst1*Mst1*pow2(s2t) - 110322*pow2(Dmglst1)*pow2(s2t) - 69349*
        pow2(Mst1)*pow2(s2t) - 3732*pow2(Mst2)*pow2(s2t) - 216*log(pow2(Mst2)/
        pow2(Mst1))*(-8*(215*Dmglst1 - 34*Mst1)*Mt*s2t + 128*pow2(Mt) + (966*
        Dmglst1*Mst1 + 475*pow2(Dmglst1) + 498*pow2(Mst1) + 186*pow2(Mst2))*
        pow2(s2t)))*pow4(Msq) + 12960*log(pow2(Mst2)/pow2(Mst1))*pow2(s2t)*
        pow6(Msq)) - Cbeta*Sbeta*(270*(-36*Mt*pow2(Mst1)*pow2(s2t) + 96*pow3(
        Mt) + 5*Mst1*pow2(Mst2)*pow3(s2t))*pow4(Mst2)*pow5(Mst1) + 18*pow4(
        Mst1)*(2*Dmglst1*((7448*Mst1*s2t*pow2(Mt) - 19461*Mt*pow2(Mst1)*pow2(
        s2t) - 23652*Mt*pow2(Mst2)*pow2(s2t) + 2160*Mt*log(pow2(Msq)/pow2(Mst1)
        )*(8*pow2(Mt) - 5*(pow2(Mst1) + pow2(Mst2))*pow2(s2t)) + 66608*pow3(Mt)
        - 3067*Mst1*pow2(Mst2)*pow3(s2t) - 36*log(pow2(Mst2)/pow2(Mst1))*(-552*
        Mst1*s2t*pow2(Mt) - Mt*(953*pow2(Mst1) + 526*pow2(Mst2))*pow2(s2t) +
        1076*pow3(Mt) + 56*Mst1*pow2(Mst2)*pow3(s2t)))*pow4(Msq) + 30*(-45*Mt*
        pow2(Mst1)*pow2(s2t) + 56*pow3(Mt) + 5*Mst1*pow2(Mst2)*pow3(s2t))*pow4(
        Mst2) - 60*pow2(Msq)*(12*pow2(Mst1)*(-9*Mt*pow2(Mst2)*pow2(s2t) + 20*
        pow3(Mt)) + (36*Mt + 17*Mst1*s2t)*pow2(s2t)*pow4(Mst2))) - pow2(
        Dmglst1)*(s2t*(85836*Mst1*Mt*s2t + 25920*Mst1*Mt*s2t*log(pow2(Msq)/
        pow2(Mst1)) + 67112*pow2(Mt) + 10809*pow2(Mst2)*pow2(s2t) - 24*log(
        pow2(Mst2)/pow2(Mst1))*(4347*Mst1*Mt*s2t + 3580*pow2(Mt) - 80*pow2(
        Mst2)*pow2(s2t)))*pow4(Msq) - 450*(-12*Mst1*Mt + s2t*pow2(Mst2))*pow2(
        s2t)*pow4(Mst2) + 180*pow2(Msq)*(64*Mst1*(-3*Mt*pow2(Mst2)*pow2(s2t) +
        5*pow3(Mt)) + 19*pow3(s2t)*pow4(Mst2))) + 360*pow2(Mst2)*pow3(s2t)*
        pow6(Msq) - 60*pow2(Msq)*(24*pow3(Mst1)*(-(Mt*pow2(Mst2)*pow2(s2t)) +
        4*pow3(Mt)) + Mst1*(24*Mt + 7*Mst1*s2t)*pow2(s2t)*pow4(Mst2) - 3*pow3(
        s2t)*pow6(Mst2))) + pow4(Msq)*(36*Mt*(22096*pow2(Mt) + 5367*pow2(Mst1)*
        pow2(s2t) + 3156*pow2(Mst2)*pow2(s2t) + 720*log(pow2(Msq)/pow2(Mst1))*(
        8*pow2(Mt) - 3*(pow2(Mst1) + pow2(Mst2))*pow2(s2t)) - 12*log(pow2(Mst2)
        /pow2(Mst1))*(1060*pow2(Mt) - 3*(99*pow2(Mst1) + 28*pow2(Mst2))*pow2(
        s2t)))*pow5(Mst1) + s2t*(12*pow2(Mst2)*(74*(473 - 36*log(pow2(Mst2)/
        pow2(Mst1)))*pow2(Mt) + (2740 + 1350*log(pow2(Msq)/pow2(Mst1)) - 4689*
        log(pow2(Mst2)/pow2(Mst1)))*pow2(Mst2)*pow2(s2t))*pow4(Mst1) + (422384*
        pow2(Mt) - 65617*pow2(Mst2)*pow2(s2t) + 1728*log(pow2(Mst2)/pow2(Mst1))
        *(64*pow2(Mt) - 39*pow2(Mst2)*pow2(s2t)))*pow6(Mst1)) - 648*pow3(s2t)*
        pow8(Mst2)))))/(pow2(Mst1)*pow4(Msq)*pow4(Mst2))) - 4*pow2(Mt)*pow2(
        MuSUSY)*((pow2(Mst1)*(83520*Mst1*Mt*s2t*pow2(Dmglst1)*pow2(Msq) +
        501120*Dmglst1*Mt*s2t*pow2(Msq)*pow2(Mst1) + 707760*Mst1*Mt*s2t*pow2(
        Dmglst1)*pow2(Mst2) + 708480*Dmglst1*Mt*s2t*pow2(Msq)*pow2(Mst2) +
        195840*Mst1*Mt*s2t*pow2(Msq)*pow2(Mst2) + 294480*Dmglst1*Mt*s2t*pow2(
        Mst1)*pow2(Mst2) + 161280*Mt*s2t*pow2(Msq)*pow3(Mst1) + 48240*Mt*s2t*
        pow2(Mst2)*pow3(Mst1) - 335256*Dmglst1*Mt*s2t*pow4(Msq) + 285120*B4*
        Dmglst1*Mt*s2t*pow4(Msq) - 2592*Dmglst1*DN*Mt*s2t*pow4(Msq) - 858392*
        Mst1*Mt*s2t*pow4(Msq) + 153792*B4*Mst1*Mt*s2t*pow4(Msq) + 864*DN*Mst1*
        Mt*s2t*pow4(Msq) - 240192*pow2(Mt)*pow4(Msq) + 3456*B4*pow2(Mt)*pow4(
        Msq) + 1728*DN*pow2(Mt)*pow4(Msq) + 25920*(Dmglst1 + Mst1)*Mt*s2t*(-2 +
        3*log(pow2(Mst2)/pow2(Mst1)))*pow2(log(pow2(Msq)/pow2(Mst1)))*pow4(Msq)
        + 576*Mt*(2*Mt - 210*Dmglst1*s2t + 107*Mst1*s2t)*pow3(log(pow2(Mst2)/
        pow2(Mst1)))*pow4(Msq) - 36*pow2(log(pow2(Mst2)/pow2(Mst1)))*(480*Mt*
        s2t*pow2(Msq)*(21*Mst1*pow2(Dmglst1) + 6*Dmglst1*pow2(Mst1) - 3*
        Dmglst1*pow2(Mst2) - Mst1*pow2(Mst2)) - 45*s2t*pow2(Mst2)*(20*Dmglst1*
        Mt*pow2(Mst1) + 5*pow2(Dmglst1)*(8*Mst1*Mt - s2t*pow2(Mst2)) + 4*Mt*
        pow3(Mst1)) + 4*Mt*(576*Mt + 1749*Dmglst1*s2t + 5641*Mst1*s2t)*pow4(
        Msq)) + 6*log(pow2(Mst2)/pow2(Mst1))*(-960*Mt*s2t*pow2(Msq)*(39*Mst1*
        pow2(Dmglst1) + 144*Dmglst1*pow2(Mst1) + 66*Dmglst1*pow2(Mst2) + 20*
        Mst1*pow2(Mst2) + 54*pow3(Mst1)) + 15*s2t*pow2(Mst2)*(8*Mst1*Mt*pow2(
        Mst2) + Dmglst1*(-1908*Mt*pow2(Mst1) + 8*Mt*pow2(Mst2)) + pow2(Dmglst1)
        *(-4152*Mst1*Mt + 537*s2t*pow2(Mst2)) - 324*Mt*pow3(Mst1)) + 32*Mt*(
        1566*Mt + 6141*Dmglst1*s2t + 10441*Mst1*s2t)*pow4(Msq)) - 15060*
        Dmglst1*Mt*s2t*pow4(Mst2) - 6420*Mst1*Mt*s2t*pow4(Mst2) - 95265*pow2(
        Dmglst1)*pow2(s2t)*pow4(Mst2) + 180*s2t*log(pow2(Msq)/pow2(Mst1))*(144*
        (9*Dmglst1 + Mst1)*Mt*pow2(log(pow2(Mst2)/pow2(Mst1)))*pow4(Msq) + 4*
        Mst1*Mt*(36*pow2(Mst1)*pow2(Mst2) + 96*pow2(Msq)*(pow2(Mst1) + pow2(
        Mst2)) - 360*pow4(Msq) - pow4(Mst2)) - 4*Dmglst1*Mt*(-180*pow2(Mst1)*
        pow2(Mst2) - 288*pow2(Msq)*(pow2(Mst1) + pow2(Mst2)) + 144*pow4(Msq) +
        pow4(Mst2)) + 9*pow2(Dmglst1)*(128*Mst1*Mt*pow2(Msq) + 160*Mst1*Mt*
        pow2(Mst2) + 15*s2t*pow4(Mst2)) - 6*log(pow2(Mst2)/pow2(Mst1))*(4*Mst1*
        Mt*(3*pow2(Mst1)*pow2(Mst2) + 8*pow2(Msq)*(3*pow2(Mst1) + pow2(Mst2)) -
        102*pow4(Msq)) + 12*Dmglst1*Mt*(5*pow2(Mst1)*pow2(Mst2) + 8*pow2(Msq)*(
        3*pow2(Mst1) + pow2(Mst2)) - 2*pow4(Msq)) + 3*pow2(Dmglst1)*(96*Mst1*
        Mt*pow2(Msq) + 40*Mst1*Mt*pow2(Mst2) - 5*s2t*pow4(Mst2))))))/(pow4(Msq)
        *pow4(Mst2)) + (pow2(s2t)*(166698000000*pow2(Dmglst1)*pow2(Msq)*pow2(
        Mst1)*pow2(Mst2) + 81496800000*Dmglst1*pow2(Msq)*pow2(Mst2)*pow3(Mst1)
        - 258631175250*pow2(Dmglst1)*pow2(Mst1)*pow4(Msq) + 11557728000*B4*
        pow2(Dmglst1)*pow2(Mst1)*pow4(Msq) + 4889808000*DN*pow2(Dmglst1)*pow2(
        Mst1)*pow4(Msq) + 53161877328*pow2(Mst1)*pow2(Mst2)*pow4(Msq) -
        16892064000*B4*pow2(Mst1)*pow2(Mst2)*pow4(Msq) + 444528000*DN*pow2(
        Mst1)*pow2(Mst2)*pow4(Msq) - 4478208000*OepS2*pow2(Mst1)*pow2(Mst2)*
        pow4(Msq) + 130691232000*S2*pow2(Mst1)*pow2(Mst2)*pow4(Msq) +
        225024292500*Dmglst1*pow3(Mst1)*pow4(Msq) - 5334336000*B4*Dmglst1*pow3(
        Mst1)*pow4(Msq) - 4445280000*Dmglst1*DN*pow3(Mst1)*pow4(Msq) +
        74088000*pow2(Mst1)*(2418*Dmglst1*Mst1 + 1873*pow2(Dmglst1) + 56*(70*
        pow2(Mst1) + 31*pow2(Mst2)))*pow3(log(pow2(Mst2)/pow2(Mst1)))*pow4(Msq)
        + 22596840000*pow2(Msq)*pow2(Mst2)*pow4(Mst1) + 531403689547*pow4(Msq)*
        pow4(Mst1) - 16892064000*B4*pow4(Msq)*pow4(Mst1) + 444528000*DN*pow4(
        Msq)*pow4(Mst1) - 8122240000*OepS2*pow4(Msq)*pow4(Mst1) + 306576144000*
        S2*pow4(Msq)*pow4(Mst1) - 12387513600*pow2(Msq)*pow2(Mst1)*pow4(Mst2) -
        44915850000*Dmglst1*pow3(Mst1)*pow4(Mst2) - 7419770460*pow4(Mst1)*pow4(
        Mst2) + 31752000*pow2(log(pow2(Msq)/pow2(Mst1)))*(-315*pow2(Mst1)*pow2(
        Mst2)*pow4(Msq) + 210*log(pow2(Mst2)/pow2(Mst1))*pow2(Mst1)*(6*Dmglst1*
        Mst1 + 3*pow2(Dmglst1) + 4*(pow2(Mst1) + pow2(Mst2)))*pow4(Msq) + pow2(
        pow2(Mst1) - pow2(Mst2))*pow4(Mst2)) + 3500599320*pow2(Mst1)*pow6(Mst2)
        - 75600*log(pow2(Msq)/pow2(Mst1))*(176400*pow2(Dmglst1)*pow2(Msq)*pow2(
        Mst1)*(19*pow2(Msq) - 5*pow2(Mst2)) - 14700*Dmglst1*pow2(Mst2)*(40*
        pow2(Msq) + pow2(Mst2))*pow3(Mst1) + 308700*pow2(Mst1)*pow2(Mst2)*pow4(
        Msq) + 88200*pow2(Mst1)*(6*Dmglst1*Mst1 + 3*pow2(Dmglst1) + 12*pow2(
        Mst1) + 8*pow2(Mst2))*pow2(log(pow2(Mst2)/pow2(Mst1)))*pow4(Msq) -
        147000*pow2(Msq)*pow2(Mst2)*pow4(Mst1) - 330750*pow4(Msq)*pow4(Mst1) -
        41160*pow2(Msq)*pow2(Mst1)*pow4(Mst2) + 12479*pow4(Mst1)*pow4(Mst2) -
        210*log(pow2(Mst2)/pow2(Mst1))*(2520*pow2(Dmglst1)*pow2(Msq)*pow2(Mst1)
        *(8*pow2(Msq) - 3*pow2(Mst2)) - 112*pow2(Msq)*pow2(Mst1)*pow2(Mst2)*(
        11*pow2(Mst1) + pow2(Mst2)) + 70*pow4(Msq)*(94*pow2(Mst1)*pow2(Mst2) +
        121*pow4(Mst1)) + 140*Dmglst1*pow3(Mst1)*(-28*pow2(Msq)*pow2(Mst2) +
        108*pow4(Msq) + 5*pow4(Mst2)) - pow4(Mst2)*(104*pow2(Mst1)*pow2(Mst2) +
        27*pow4(Mst1) + 18*pow4(Mst2)) + 1680*pow2(Mst1)*pow6(Msq)) - 26218*
        pow2(Mst1)*pow6(Mst2) - 9571*pow8(Mst2)) + 1331874540*pow8(Mst2) -
        264600*pow2(log(pow2(Mst2)/pow2(Mst1)))*(210*pow2(Dmglst1)*pow2(Msq)*
        pow2(Mst1)*(8947*pow2(Msq) - 1080*pow2(Mst2)) + pow4(Msq)*(1487556*
        pow2(Mst1)*pow2(Mst2) + 2258059*pow4(Mst1)) + 140*Dmglst1*pow3(Mst1)*(-
        840*pow2(Msq)*pow2(Mst2) + 9491*pow4(Msq) + 150*pow4(Mst2)) - 840*pow2(
        Msq)*(104*pow2(Mst2)*pow4(Mst1) + 29*pow2(Mst1)*pow4(Mst2)) - 30*(321*
        pow4(Mst1)*pow4(Mst2) + 216*pow2(Mst1)*pow6(Mst2) + 32*pow8(Mst2))) +
        1260*log(pow2(Mst2)/pow2(Mst1))*(36750*pow2(Dmglst1)*pow2(Msq)*pow2(
        Mst1)*(14027*pow2(Msq) - 3024*pow2(Mst2)) + 4900*Dmglst1*pow3(Mst1)*(-
        6960*pow2(Msq)*pow2(Mst2) + 68551*pow4(Msq) + 5010*pow4(Mst2)) + 3*(
        pow2(Mst1)*(3*(72625277 + 14504000*S2)*pow2(Mst1) + 392*(372457 +
        61200*S2)*pow2(Mst2))*pow4(Msq) - 3920*pow2(Msq)*(1928*pow2(Mst2)*pow4(
        Mst1) + 353*pow2(Mst1)*pow4(Mst2)) + 3528000*pow2(Mst1)*pow6(Msq) - 10*
        (5151*pow4(Mst1)*pow4(Mst2) + 126132*pow2(Mst1)*pow6(Mst2) + 31294*
        pow8(Mst2))))))/(1.029e6*pow4(Msq)*pow4(Mst2))) + 3888*pow2(Sbeta)*(
        pow4(s2t)*(-(Dmglst1*pow3(Mst1)*(773389 + 10368*B4 + 8640*DN + 378084*
        log(pow2(Mst2)/pow2(Mst1)) - 38880*(-1 + 2*log(pow2(Mst2)/pow2(Mst1)))*
        pow2(log(pow2(Msq)/pow2(Mst1))) - 44424*pow2(log(pow2(Mst2)/pow2(Mst1))
        ) + 38880*log(pow2(Msq)/pow2(Mst1))*(9 + 4*log(pow2(Mst2)/pow2(Mst1)) +
        2*pow2(log(pow2(Mst2)/pow2(Mst1)))) + 89568*pow3(log(pow2(Mst2)/pow2(
        Mst1)))))/15552. + pow2(Mst1)*((5*(1 + 2*log(pow2(Msq)/pow2(Mst1)))*(-1
        + 2*log(pow2(Mst2)/pow2(Mst1)))*pow2(Msq))/6. - (pow2(Dmglst1)*(-438203
        + 21888*B4 + 12096*DN + 606756*log(pow2(Mst2)/pow2(Mst1)) - 12960*(-1 +
        2*log(pow2(Mst2)/pow2(Mst1)))*pow2(log(pow2(Msq)/pow2(Mst1))) - 431208*
        pow2(log(pow2(Mst2)/pow2(Mst1))) + 1440*log(pow2(Msq)/pow2(Mst1))*(67 +
        216*log(pow2(Mst2)/pow2(Mst1)) + 18*pow2(log(pow2(Mst2)/pow2(Mst1)))) +
        93600*pow3(log(pow2(Mst2)/pow2(Mst1)))))/10368.) - (32.221840780308305
         + (10*B4)/
        9. + DN/18. - (1737319*log(pow2(Mst2)/pow2(Mst1)))/705600. - (5*pow2(
        Msq))/(6.*pow2(Mst2)) - (5*(-1 + 7*log(pow2(Mst2)/pow2(Mst1)))*pow2(
        log(pow2(Msq)/pow2(Mst1))))/12. + (256523*pow2(log(pow2(Mst2)/pow2(
        Mst1))))/30240. + (5*log(pow2(Msq)/pow2(Mst1))*(55 - 64*log(pow2(Mst2)/
        pow2(Mst1)) - (24*pow2(Msq))/pow2(Mst2) + 30*pow2(log(pow2(Mst2)/pow2(
        Mst1)))))/72. + (1535*pow3(log(pow2(Mst2)/pow2(Mst1))))/216.)*pow4(
        Mst1) - (pow2(Mst2)*(-2589750*Dmglst1*Mst1*pow2(Msq) - 162000*B4*
        Dmglst1*Mst1*pow2(Msq) - 135000*Dmglst1*DN*Mst1*pow2(Msq) + 6620625*
        pow2(Dmglst1)*pow2(Msq) - 81000*B4*pow2(Dmglst1)*pow2(Msq) - 67500*DN*
        pow2(Dmglst1)*pow2(Msq) + 8662500*pow2(Dmglst1)*pow2(Mst1) + 11753176*
        pow2(Msq)*pow2(Mst1) - 27000*B4*pow2(Msq)*pow2(Mst1) - 40500*DN*pow2(
        Msq)*pow2(Mst1) + 101250*pow2(Msq)*(6*Dmglst1*Mst1 + 3*pow2(Dmglst1) +
        7*pow2(Mst1) + 6*log(pow2(Mst2)/pow2(Mst1))*pow2(Dmglst1 + Mst1))*pow2(
        log(pow2(Msq)/pow2(Mst1))) + 6165000*Dmglst1*pow3(Mst1) + 2250*pow2(
        Msq)*(898*Dmglst1*Mst1 + 449*pow2(Dmglst1) + 1097*pow2(Mst1))*pow3(log(
        pow2(Mst2)/pow2(Mst1))) + 202500*pow4(Msq) + 180*log(pow2(Mst2)/pow2(
        Mst1))*(150*pow2(Dmglst1)*(143*pow2(Msq) - 215*pow2(Mst1)) + 4983*pow2(
        Msq)*pow2(Mst1) + 50*Dmglst1*(513*Mst1*pow2(Msq) - 490*pow3(Mst1)) +
        2250*pow4(Msq) - 4955*pow4(Mst1)) - 1350*pow2(log(pow2(Mst2)/pow2(Mst1)
        ))*(3630*Dmglst1*Mst1*pow2(Msq) + 1815*pow2(Dmglst1)*pow2(Msq) - 1500*
        pow2(Dmglst1)*pow2(Mst1) + 3372*pow2(Msq)*pow2(Mst1) - 1000*Dmglst1*
        pow3(Mst1) - 280*pow4(Mst1)) - 6750*log(pow2(Msq)/pow2(Mst1))*(-265*
        pow2(Msq)*pow2(Mst1) + 5*pow2(Dmglst1)*(203*pow2(Msq) + 36*pow2(Mst1))
        + 30*pow2(Msq)*(6*Dmglst1*Mst1 + 3*pow2(Dmglst1) + 5*pow2(Mst1))*pow2(
        log(pow2(Mst2)/pow2(Mst1))) + 10*Dmglst1*(63*Mst1*pow2(Msq) - 4*pow3(
        Mst1)) - 60*pow4(Msq) - 2*log(pow2(Mst2)/pow2(Mst1))*(30*pow2(Dmglst1)*
        (3*pow2(Msq) - 10*pow2(Mst1)) + 60*pow2(Msq)*pow2(Mst1) + 20*Dmglst1*(
        9*Mst1*pow2(Msq) - 10*pow3(Mst1)) + 60*pow4(Msq) - 41*pow4(Mst1)) - 46*
        pow4(Mst1)) + 1861200*pow4(Mst1)))/(243000.*pow2(Msq)) - (2*OepS2*(-
        150*pow2(Mst1)*pow2(Mst2) + 11*pow4(Mst1) - 27*pow4(Mst2)))/729. + (S2*
        (-3546*pow2(Mst1)*pow2(Mst2) + 281*pow4(Mst1) + 6*log(pow2(Mst2)/pow2(
        Mst1))*(-150*pow2(Mst1)*pow2(Mst2) + 11*pow4(Mst1) - 27*pow4(Mst2)) +
        891*pow4(Mst2)))/108. + (pow4(Mst2)*(164640000*pow2(Dmglst1)*pow2(Msq)*
        pow2(Mst1) + 168756000*Dmglst1*pow2(Msq)*pow3(Mst1) + 2160900*Dmglst1*
        Mst1*pow4(Msq) + 330051750*pow2(Dmglst1)*pow4(Msq) + 460897675*pow2(
        Mst1)*pow4(Msq) + 22226400*B4*pow2(Mst1)*pow4(Msq) - 2469600*DN*pow2(
        Mst1)*pow4(Msq) - 65753100*pow2(Mst1)*pow3(log(pow2(Mst2)/pow2(Mst1)))*
        pow4(Msq) + 276929625*pow2(Dmglst1)*pow4(Mst1) + 111708240*pow2(Msq)*
        pow4(Mst1) + 105*log(pow2(Mst2)/pow2(Mst1))*(1470*pow2(Dmglst1)*(-460*
        pow2(Msq)*pow2(Mst1) + 438*pow4(Msq) - 895*pow4(Mst1)) - pow2(Mst1)*(
        457464*pow2(Msq)*pow2(Mst1) + 6962410*pow4(Msq) + 215819*pow4(Mst1)) +
        980*Dmglst1*(-1380*pow2(Msq)*pow3(Mst1) + 678*Mst1*pow4(Msq) - 835*
        pow5(Mst1))) + 134027250*Dmglst1*pow5(Mst1) + 18522000*pow6(Msq) +
        264600*pow2(log(pow2(Msq)/pow2(Mst1)))*(210*Dmglst1*Mst1*pow4(Msq) +
        105*pow2(Dmglst1)*pow4(Msq) + 175*pow2(Mst1)*pow4(Msq) - 35*log(pow2(
        Mst2)/pow2(Mst1))*pow2(Mst1)*pow4(Msq) - 2*pow6(Mst1)) + 36358596*pow6(
        Mst1) + 22050*pow2(log(pow2(Mst2)/pow2(Mst1)))*(28*Dmglst1*Mst1*pow2(3*
        pow2(Msq) + 5*pow2(Mst1)) + 15351*pow2(Mst1)*pow4(Msq) - 252*pow2(Msq)*
        pow4(Mst1) + 42*pow2(Dmglst1)*(10*pow2(Msq)*pow2(Mst1) + 3*pow4(Msq) +
        25*pow4(Mst1)) + 79*pow6(Mst1)) - 1260*log(pow2(Msq)/pow2(Mst1))*(-
        22050*pow2(Mst1)*pow2(log(pow2(Mst2)/pow2(Mst1)))*pow4(Msq) + 1225*
        pow2(Dmglst1)*(96*pow2(Msq)*pow2(Mst1) + 74*pow4(Msq) + 45*pow4(Mst1))
        + 35*log(pow2(Mst2)/pow2(Mst1))*pow2(Mst1)*(168*pow2(Msq)*pow2(Mst1) +
        210*pow2(Dmglst1)*(2*pow2(Msq) + 5*pow2(Mst1)) + 140*Dmglst1*(6*Mst1*
        pow2(Msq) + 5*pow3(Mst1)) + 2450*pow4(Msq) + 163*pow4(Mst1)) + 2450*
        Dmglst1*(16*pow2(Msq)*pow3(Mst1) + 30*Mst1*pow4(Msq) + pow5(Mst1)) - 2*
        (45325*pow2(Mst1)*pow4(Msq) + 8330*pow2(Msq)*pow4(Mst1) + 14700*pow6(
        Msq) + 4612*pow6(Mst1)))))/(2.22264e7*pow2(Mst1)*pow4(Msq))) + (Mt*
        pow3(s2t)*(-73152*pow2(Dmglst1)*pow2(Mst1)*pow2(Mst2)*pow4(Msq) -
        196992*B4*pow2(Dmglst1)*pow2(Mst1)*pow2(Mst2)*pow4(Msq) + 5184*DN*pow2(
        Dmglst1)*pow2(Mst1)*pow2(Mst2)*pow4(Msq) - 25920*Mst1*(Dmglst1 + Mst1)*
        pow2(Mst2)*(2*(pow2(Mst1) - pow2(Mst2)) + log(pow2(Mst2)/pow2(Mst1))*(
        pow2(Mst1) + pow2(Mst2)))*pow2(log(pow2(Msq)/pow2(Mst1)))*pow4(Msq) -
        1297272*Dmglst1*pow2(Mst2)*pow3(Mst1)*pow4(Msq) - 226368*B4*Dmglst1*
        pow2(Mst2)*pow3(Mst1)*pow4(Msq) + 4320*Dmglst1*DN*pow2(Mst2)*pow3(Mst1)
        *pow4(Msq) + 1707840*pow2(Dmglst1)*pow2(Msq)*pow2(Mst2)*pow4(Mst1) -
        6248762*pow2(Dmglst1)*pow4(Msq)*pow4(Mst1) - 1059256*pow2(Mst2)*pow4(
        Msq)*pow4(Mst1) - 95040*B4*pow2(Mst2)*pow4(Msq)*pow4(Mst1) + 864*DN*
        pow2(Mst2)*pow4(Msq)*pow4(Mst1) - 895680*pow2(Dmglst1)*pow2(Msq)*pow2(
        Mst1)*pow4(Mst2) - 725760*Dmglst1*pow2(Msq)*pow3(Mst1)*pow4(Mst2) +
        816264*Dmglst1*Mst1*pow4(Msq)*pow4(Mst2) - 29376*B4*Dmglst1*Mst1*pow4(
        Msq)*pow4(Mst2) - 864*Dmglst1*DN*Mst1*pow4(Msq)*pow4(Mst2) - 846000*
        pow2(Dmglst1)*pow4(Msq)*pow4(Mst2) + 958824*pow2(Mst1)*pow4(Msq)*pow4(
        Mst2) - 29376*B4*pow2(Mst1)*pow4(Msq)*pow4(Mst2) - 864*DN*pow2(Mst1)*
        pow4(Msq)*pow4(Mst2) - 733680*pow2(Dmglst1)*pow4(Mst1)*pow4(Mst2) -
        213120*pow2(Msq)*pow4(Mst1)*pow4(Mst2) + 288*Mst1*pow3(log(pow2(Mst2)/
        pow2(Mst1)))*pow4(Msq)*(pow2(Dmglst1)*(951*Mst1*pow2(Mst2) + 366*pow3(
        Mst1)) + Dmglst1*(712*pow2(Mst1)*pow2(Mst2) + 323*pow4(Mst1) - 146*
        pow4(Mst2)) - Mst1*(-78*pow2(Mst1)*pow2(Mst2) + pow4(Mst1) + 146*pow4(
        Mst2))) + 924480*Dmglst1*pow2(Msq)*pow2(Mst2)*pow5(Mst1) - 3058187*
        Dmglst1*pow4(Msq)*pow5(Mst1) - 325260*Dmglst1*pow4(Mst2)*pow5(Mst1) -
        72*Mst1*pow2(log(pow2(Mst2)/pow2(Mst1)))*(2*Mst1*pow2(Dmglst1)*((10585*
        pow2(Mst1) + 7074*pow2(Mst2))*pow4(Msq) - 360*pow2(Msq)*(9*pow2(Mst1)*
        pow2(Mst2) - pow4(Mst2)) + 450*pow2(Mst1)*pow4(Mst2)) + Dmglst1*(pow4(
        Msq)*(6066*pow2(Mst1)*pow2(Mst2) + 5569*pow4(Mst1) - 4782*pow4(Mst2)) +
        450*pow4(Mst1)*pow4(Mst2) - 720*pow2(Msq)*(4*pow2(Mst2)*pow4(Mst1) -
        pow2(Mst1)*pow4(Mst2))) + 90*pow4(Mst2)*pow5(Mst1) - pow4(Msq)*(1718*
        pow2(Mst2)*pow3(Mst1) + 4782*Mst1*pow4(Mst2) + 2579*pow5(Mst1)) - 240*
        pow2(Msq)*(-(pow3(Mst1)*pow4(Mst2)) + 2*pow2(Mst2)*pow5(Mst1))) +
        239040*pow2(Msq)*pow2(Mst2)*pow6(Mst1) - 622151*pow4(Msq)*pow6(Mst1) -
        61740*pow4(Mst2)*pow6(Mst1) + 8640*Dmglst1*Mst1*pow2(Msq)*pow6(Mst2) +
        12960*pow2(Dmglst1)*pow2(Mst1)*pow6(Mst2) + 8640*pow2(Msq)*pow2(Mst1)*
        pow6(Mst2) + 16380*Dmglst1*pow3(Mst1)*pow6(Mst2) + 7740*pow4(Mst1)*
        pow6(Mst2) + 2160*log(pow2(Msq)/pow2(Mst1))*(-12*Mst1*pow2(Mst2)*pow2(
        log(pow2(Mst2)/pow2(Mst1)))*(12*Mst1*pow2(Dmglst1) + 11*Dmglst1*pow2(
        Mst1) - Dmglst1*pow2(Mst2) - Mst1*pow2(Mst2) + 3*pow3(Mst1))*pow4(Msq)
        - 168*pow2(Mst2)*pow4(Msq)*pow4(Mst1) + 144*pow2(Mst1)*pow4(Msq)*pow4(
        Mst2) - 32*pow2(Msq)*pow4(Mst1)*pow4(Mst2) + 4*pow2(Dmglst1)*(24*pow2(
        Msq)*pow2(Mst1)*(pow2(Mst1) - pow2(Mst2))*pow2(Mst2) - 30*pow4(Mst1)*
        pow4(Mst2) + pow4(Msq)*(-155*pow2(Mst1)*pow2(Mst2) + 24*pow4(Mst1) +
        17*pow4(Mst2))) + 2*Mst1*log(pow2(Mst2)/pow2(Mst1))*(8*pow2(Msq)*pow2(
        Mst2)*(pow2(Mst1) + pow2(Mst2))*pow3(Mst1) + 6*Mst1*pow2(Dmglst1)*(4*
        pow2(Msq)*pow2(Mst2)*(pow2(Mst1) + pow2(Mst2)) + (78*pow2(Mst1) + 40*
        pow2(Mst2))*pow4(Msq) + 5*pow2(Mst1)*pow4(Mst2)) + 3*Dmglst1*(8*pow2(
        Msq)*pow2(Mst1)*pow2(Mst2)*(pow2(Mst1) + pow2(Mst2)) + 2*pow4(Msq)*(13*
        pow2(Mst1)*pow2(Mst2) + 33*pow4(Mst1) - 7*pow4(Mst2)) + 5*pow4(Mst1)*
        pow4(Mst2)) + 3*pow4(Mst2)*pow5(Mst1) + 6*pow4(Msq)*(-3*pow2(Mst2)*
        pow3(Mst1) - 7*Mst1*pow4(Mst2) + 5*pow5(Mst1))) + 32*pow2(Msq)*pow2(
        Mst2)*pow6(Mst1) + 66*pow4(Msq)*pow6(Mst1) - 13*pow4(Mst2)*pow6(Mst1) +
        pow4(Mst1)*pow6(Mst2) + Dmglst1*(96*pow2(Msq)*(pow2(Mst1) - pow2(Mst2))
        *pow2(Mst2)*pow3(Mst1) - 61*pow4(Mst2)*pow5(Mst1) + 6*pow4(Msq)*(-64*
        pow2(Mst2)*pow3(Mst1) + 36*Mst1*pow4(Mst2) + 43*pow5(Mst1)) + pow3(
        Mst1)*pow6(Mst2))) - 12*log(pow2(Mst2)/pow2(Mst1))*(90*(-29*pow2(Mst1)
        + 2*pow2(Mst2))*pow4(Mst1)*pow4(Mst2) + 2*pow2(Dmglst1)*(-15570*pow4(
        Mst1)*pow4(Mst2) + pow4(Msq)*(-90288*pow2(Mst1)*pow2(Mst2) + 87319*
        pow4(Mst1) + 3660*pow4(Mst2)) + 720*pow2(Msq)*(33*pow2(Mst2)*pow4(Mst1)
        - 23*pow2(Mst1)*pow4(Mst2))) + pow4(Msq)*(13336*pow2(Mst2)*pow4(Mst1) +
        76860*pow2(Mst1)*pow4(Mst2) + 11313*pow6(Mst1)) - 960*pow2(Msq)*(10*
        pow4(Mst1)*pow4(Mst2) + 7*pow2(Mst2)*pow6(Mst1)) + 3*Dmglst1*(-4830*
        pow4(Mst2)*pow5(Mst1) + pow4(Msq)*(-23192*pow2(Mst2)*pow3(Mst1) +
        27972*Mst1*pow4(Mst2) + 26687*pow5(Mst1)) - 960*pow2(Msq)*(11*pow3(
        Mst1)*pow4(Mst2) + 2*pow2(Mst2)*pow5(Mst1)) + 60*pow3(Mst1)*pow6(Mst2))
        )))/(3888.*Mst1*pow2(Mst2)*pow4(Msq)) - (pow4(Mt)*(77630943040*pow2(
        Dmglst1)*pow2(Mst1)*pow2(Mst2)*pow4(Msq) + 43967798000*Dmglst1*pow2(
        Mst2)*pow3(Mst1)*pow4(Msq) + 10187100000*pow2(Dmglst1)*pow2(Msq)*pow2(
        Mst2)*pow4(Mst1) + 155148099260*pow2(Dmglst1)*pow4(Msq)*pow4(Mst1) +
        13498218944*pow2(Mst2)*pow4(Msq)*pow4(Mst1) + 10298232000*pow2(MuSUSY)*
        pow4(Msq)*pow4(Mst1) - 148176000*B4*pow2(MuSUSY)*pow4(Msq)*pow4(Mst1) -
        74088000*DN*pow2(MuSUSY)*pow4(Msq)*pow4(Mst1) + 3269750400*pow2(
        Dmglst1)*pow2(Msq)*pow2(Mst1)*pow4(Mst2) + 5000940000*Dmglst1*pow2(Msq)
        *pow3(Mst1)*pow4(Mst2) + 2279235000*Dmglst1*Mst1*pow4(Msq)*pow4(Mst2) -
        7271696040*pow2(Dmglst1)*pow4(Msq)*pow4(Mst2) - 5402764500*pow2(Mst1)*
        pow4(Msq)*pow4(Mst2) + 185220000*pow2(Mst1)*pow3(log(pow2(Msq)/pow2(
        Mst1)))*pow4(Msq)*pow4(Mst2) + 5691038850*pow2(Dmglst1)*pow4(Mst1)*
        pow4(Mst2) + 3047239440*pow2(Msq)*pow4(Mst1)*pow4(Mst2) - 6174000*pow2(
        Mst1)*pow3(log(pow2(Mst2)/pow2(Mst1)))*pow4(Msq)*(2*pow2(Dmglst1)*(556*
        pow2(Mst1) - 33*pow2(Mst2)) + pow2(Mst1)*(-70*pow2(Mst2) + 8*pow2(
        MuSUSY)) + 12*Dmglst1*(-11*Mst1*pow2(Mst2) + 57*pow3(Mst1)) + 145*pow4(
        Mst1) + 17*pow4(Mst2)) + 6791400000*Dmglst1*pow2(Msq)*pow2(Mst2)*pow5(
        Mst1) + 102884565000*Dmglst1*pow4(Msq)*pow5(Mst1) + 5389387500*Dmglst1*
        pow4(Mst2)*pow5(Mst1) - 555660000*pow2(Mst1)*pow2(Mst2)*pow6(Msq) -
        555660000*pow4(Mst2)*pow6(Msq) + 1697850000*pow2(Msq)*pow2(Mst2)*pow6(
        Mst1) + 25346116504*pow4(Msq)*pow6(Mst1) + 2412641055*pow4(Mst2)*pow6(
        Mst1) + 887512500*pow2(Dmglst1)*pow2(Mst1)*pow6(Mst2) + 1112719440*
        pow2(Msq)*pow2(Mst1)*pow6(Mst2) + 1775025000*Dmglst1*pow3(Mst1)*pow6(
        Mst2) - 55566000*pow4(Msq)*pow6(Mst2) + 887512500*pow4(Mst1)*pow6(Mst2)
        - 1323000*pow2(Mst2)*pow2(log(pow2(Msq)/pow2(Mst1)))*(210*log(pow2(
        Mst2)/pow2(Mst1))*pow2(Mst1)*pow2(Mst2)*pow4(Msq) + 210*pow2(Dmglst1)*(
        4*pow2(Msq)*pow2(Mst1)*pow2(Mst2) + pow2(Mst1)*pow2(Mst2)*(11*pow2(
        Mst1) + pow2(Mst2)) + 3*(pow2(Mst1) + pow2(Mst2))*pow4(Msq)) + 140*
        Dmglst1*(12*pow2(Msq)*pow2(Mst2)*pow3(Mst1) + 9*Mst1*(pow2(Mst1) +
        pow2(Mst2))*pow4(Msq) + 3*pow3(Mst1)*pow4(Mst2) + 13*pow2(Mst2)*pow5(
        Mst1)) + pow2(Mst1)*(70*(9*pow2(Mst1) + 11*pow2(Mst2))*pow4(Msq) + 617*
        pow2(Mst2)*pow4(Mst1) + 210*pow2(Mst1)*pow4(Mst2) + 28*pow2(Msq)*(37*
        pow2(Mst1)*pow2(Mst2) + 7*pow4(Mst2)) + 57*pow6(Mst2))) + 611891055*
        pow2(Mst1)*pow8(Mst2) + 6300*log(pow2(Msq)/pow2(Mst1))*(-396900*pow2(
        Mst2)*pow4(Msq)*pow4(Mst1) + 441000*(4*Dmglst1*Mst1 + 6*pow2(Dmglst1) +
        pow2(Mst1))*pow2(log(pow2(Mst2)/pow2(Mst1)))*pow4(Msq)*pow4(Mst1) -
        401800*pow2(Mst1)*pow4(Msq)*pow4(Mst2) + 45472*pow2(Msq)*pow4(Mst1)*
        pow4(Mst2) + 294*pow2(Dmglst1)*(25*pow2(Mst1)*(163*pow2(Mst1) + 8*pow2(
        Mst2))*pow4(Mst2) + 2*pow4(Msq)*(-675*pow2(Mst1)*pow2(Mst2) + 7050*
        pow4(Mst1) + 752*pow4(Mst2)) + 200*pow2(Msq)*(15*pow2(Mst2)*pow4(Mst1)
        + 13*pow2(Mst1)*pow4(Mst2))) - 176400*pow2(Mst1)*pow2(Mst2)*pow6(Msq) -
        176400*pow4(Mst2)*pow6(Msq) + 147000*pow2(Msq)*pow2(Mst2)*pow6(Mst1) +
        176400*pow4(Msq)*pow6(Mst1) + 220709*pow4(Mst2)*pow6(Mst1) + 16072*
        pow2(Msq)*pow2(Mst1)*pow6(Mst2) + 58800*pow4(Mst1)*pow6(Mst2) + 4900*
        Dmglst1*(137*pow4(Mst2)*pow5(Mst1) + 2*pow4(Msq)*(-81*pow2(Mst2)*pow3(
        Mst1) - 5*Mst1*pow4(Mst2) + 162*pow5(Mst1)) + 24*pow2(Msq)*(3*pow3(
        Mst1)*pow4(Mst2) + 5*pow2(Mst2)*pow5(Mst1)) + 24*pow3(Mst1)*pow6(Mst2))
        - 210*log(pow2(Mst2)/pow2(Mst1))*pow2(Mst1)*(70*pow4(Msq)*(72*pow4(
        Mst1) + 5*pow4(Mst2)) - pow4(Mst2)*(210*pow2(Mst1)*pow2(Mst2) + 175*
        pow4(Mst1) + 9*pow4(Mst2)) + 210*pow2(Dmglst1)*(214*pow2(Mst1)*pow4(
        Msq) - 2*pow2(Msq)*pow4(Mst2) - (5*pow2(Mst1) + pow2(Mst2))*pow4(Mst2))
        + 140*Dmglst1*Mst1*(174*pow2(Mst1)*pow4(Msq) - 6*pow2(Msq)*pow4(Mst2) -
        5*pow2(Mst1)*pow4(Mst2) - 3*pow6(Mst2)) - 28*pow2(Msq)*(15*pow2(Mst1)*
        pow4(Mst2) + 4*pow6(Mst2))) + 49209*pow2(Mst1)*pow8(Mst2)) + 88200*
        pow2(log(pow2(Mst2)/pow2(Mst1)))*(70*Dmglst1*Mst1*pow4(Msq)*(-706*pow2(
        Mst1)*pow2(Mst2) + 2898*pow4(Mst1) - 27*pow4(Mst2)) + 35*pow2(Dmglst1)*
        pow4(Msq)*(-706*pow2(Mst1)*pow2(Mst2) + 9036*pow4(Mst1) - 27*pow4(Mst2)
        ) + pow2(Mst1)*(pow4(Msq)*(-28*pow2(Mst1)*(823*pow2(Mst2) - 1440*pow2(
        MuSUSY)) + 44291*pow4(Mst1) - 105*pow4(Mst2)) + 1260*pow2(Msq)*pow6(
        Mst2) + 720*pow8(Mst2))) - 210*log(pow2(Mst2)/pow2(Mst1))*(4900*
        Dmglst1*(180*pow2(Msq)*pow3(Mst1)*pow4(Mst2) + 15*(55*pow2(Mst1) + 69*
        pow2(Mst2))*pow3(Mst1)*pow4(Mst2) + pow4(Msq)*(-8276*pow2(Mst2)*pow3(
        Mst1) + 254*Mst1*pow4(Mst2) + 34432*pow5(Mst1))) + 4*pow4(Msq)*(-98*(
        54223*pow2(Mst2) - 156600*pow2(MuSUSY))*pow4(Mst1) - 1642725*pow2(Mst1)
        *pow4(Mst2) + 9874311*pow6(Mst1) + 132300*pow6(Mst2)) + 5880*pow2(Msq)*
        (75*pow4(Mst1)*pow4(Mst2) + 482*pow2(Mst1)*pow6(Mst2)) + 98*pow2(
        Dmglst1)*(pow4(Msq)*(-220500*pow2(Mst1)*pow2(Mst2) + 2323150*pow4(Mst1)
        - 6336*pow4(Mst2)) + 4500*pow2(Msq)*pow2(Mst1)*pow4(Mst2) + 1125*(55*
        pow4(Mst1)*pow4(Mst2) + 23*pow2(Mst1)*pow6(Mst2))) + 15*(67375*pow4(
        Mst2)*pow6(Mst1) + 169050*pow4(Mst1)*pow6(Mst2) + 131493*pow2(Mst1)*
        pow8(Mst2)))))/(4.16745e7*pow2(Mst1)*pow4(Msq)*pow4(Mst2)) + (s2t*pow3(
        Mt)*(177120*Dmglst1*pow2(Msq)*pow2(Mst2)*pow2(MuSUSY)*pow3(Mst1) +
        129720*pow2(Dmglst1)*pow2(Mst1)*pow2(Mst2)*pow4(Msq) + 98884*Dmglst1*
        pow2(Mst2)*pow3(Mst1)*pow4(Msq) - 83814*Dmglst1*pow2(MuSUSY)*pow3(Mst1)
        *pow4(Msq) + 71280*B4*Dmglst1*pow2(MuSUSY)*pow3(Mst1)*pow4(Msq) - 648*
        Dmglst1*DN*pow2(MuSUSY)*pow3(Mst1)*pow4(Msq) - 79920*pow2(Dmglst1)*
        pow2(Msq)*pow2(Mst2)*pow4(Mst1) + 20880*pow2(Dmglst1)*pow2(Msq)*pow2(
        MuSUSY)*pow4(Mst1) + 176940*pow2(Dmglst1)*pow2(Mst2)*pow2(MuSUSY)*pow4(
        Mst1) + 48960*pow2(Msq)*pow2(Mst2)*pow2(MuSUSY)*pow4(Mst1) + 1395788*
        pow2(Dmglst1)*pow4(Msq)*pow4(Mst1) + 996*pow2(Mst2)*pow4(Msq)*pow4(
        Mst1) - 214598*pow2(MuSUSY)*pow4(Msq)*pow4(Mst1) + 38448*B4*pow2(
        MuSUSY)*pow4(Msq)*pow4(Mst1) + 216*DN*pow2(MuSUSY)*pow4(Msq)*pow4(Mst1)
        + 146160*pow2(Dmglst1)*pow2(Msq)*pow2(Mst1)*pow4(Mst2) + 92160*Dmglst1*
        pow2(Msq)*pow3(Mst1)*pow4(Mst2) - 3765*Dmglst1*pow2(MuSUSY)*pow3(Mst1)*
        pow4(Mst2) + 84084*Dmglst1*Mst1*pow4(Msq)*pow4(Mst2) + 357900*pow2(
        Dmglst1)*pow4(Msq)*pow4(Mst2) + 13356*pow2(Mst1)*pow4(Msq)*pow4(Mst2) +
        133560*pow2(Dmglst1)*pow4(Mst1)*pow4(Mst2) + 26400*pow2(Msq)*pow4(Mst1)
        *pow4(Mst2) - 1605*pow2(MuSUSY)*pow4(Mst1)*pow4(Mst2) - 23760*Dmglst1*
        pow2(Msq)*pow2(Mst2)*pow5(Mst1) + 125280*Dmglst1*pow2(Msq)*pow2(MuSUSY)
        *pow5(Mst1) + 73620*Dmglst1*pow2(Mst2)*pow2(MuSUSY)*pow5(Mst1) +
        216148*Dmglst1*pow4(Msq)*pow5(Mst1) + 66240*Dmglst1*pow4(Mst2)*pow5(
        Mst1) - 72*Mst1*pow3(log(pow2(Mst2)/pow2(Mst1)))*pow4(Msq)*((279*pow2(
        Mst2) - 214*pow2(MuSUSY))*pow3(Mst1) + pow2(Dmglst1)*(849*Mst1*pow2(
        Mst2) + 3393*pow3(Mst1)) + Dmglst1*(5*pow2(Mst1)*(169*pow2(Mst2) + 84*
        pow2(MuSUSY)) + 1830*pow4(Mst1) - 14*pow4(Mst2)) - 14*Mst1*pow4(Mst2) +
        384*pow5(Mst1)) + 36*Mst1*pow2(log(pow2(Mst2)/pow2(Mst1)))*(-120*pow2(
        Msq)*pow2(Mst2)*(3*pow2(Mst1) - pow2(MuSUSY))*pow3(Mst1) + pow2(
        Dmglst1)*(450*pow2(Mst2)*pow2(MuSUSY)*pow3(Mst1) - 360*pow2(Msq)*(10*
        pow2(Mst2) + 7*pow2(MuSUSY))*pow3(Mst1) + (6471*Mst1*pow2(Mst2) +
        43595*pow3(Mst1))*pow4(Msq)) + Dmglst1*(225*pow2(Mst2)*pow2(MuSUSY)*
        pow4(Mst1) - 360*pow2(Msq)*(-(pow2(Mst1)*pow2(Mst2)*pow2(MuSUSY)) + (5*
        pow2(Mst2) + 2*pow2(MuSUSY))*pow4(Mst1)) + pow4(Msq)*(pow2(Mst1)*(6536*
        pow2(Mst2) - 1749*pow2(MuSUSY)) + 23091*pow4(Mst1) - 96*pow4(Mst2))) +
        45*pow2(Mst2)*pow2(MuSUSY)*pow5(Mst1) + pow4(Msq)*((2094*pow2(Mst2) -
        5641*pow2(MuSUSY))*pow3(Mst1) - 96*Mst1*pow4(Mst2) + 4649*pow5(Mst1)))
        - 2160*pow2(Msq)*pow2(Mst2)*pow6(Mst1) + 40320*pow2(Msq)*pow2(MuSUSY)*
        pow6(Mst1) + 12060*pow2(Mst2)*pow2(MuSUSY)*pow6(Mst1) - 81172*pow4(Msq)
        *pow6(Mst1) + 18000*pow4(Mst2)*pow6(Mst1) - 8040*Dmglst1*Mst1*pow2(Msq)
        *pow6(Mst2) - 20160*pow2(Dmglst1)*pow2(Mst1)*pow6(Mst2) - 8040*pow2(
        Msq)*pow2(Mst1)*pow6(Mst2) - 20160*Dmglst1*pow3(Mst1)*pow6(Mst2) -
        6720*pow4(Mst1)*pow6(Mst2) - 6*log(pow2(Mst2)/pow2(Mst1))*(2*pow2(
        Dmglst1)*(pow4(Msq)*(28427*pow2(Mst1)*pow2(Mst2) + 181409*pow4(Mst1) -
        777*pow4(Mst2)) - 360*pow2(Msq)*(13*(4*pow2(Mst2) - pow2(MuSUSY))*pow4(
        Mst1) - 8*pow2(Mst1)*pow4(Mst2)) + 45*(pow4(Mst1)*(173*pow2(Mst2)*pow2(
        MuSUSY) + 35*pow4(Mst2)) - 17*pow2(Mst1)*pow6(Mst2))) + pow2(Mst1)*(2*
        pow4(Msq)*(pow2(Mst1)*(6591*pow2(Mst2) - 41764*pow2(MuSUSY)) + 17418*
        pow4(Mst1) + 468*pow4(Mst2)) - 120*pow2(Msq)*(3*(5*pow2(Mst2) - 36*
        pow2(MuSUSY))*pow4(Mst1) - 8*pow2(Mst1)*(5*pow2(Mst2)*pow2(MuSUSY) + 2*
        pow4(Mst2)) + 4*pow6(Mst2)) + 15*(-2*pow2(Mst1)*(17*pow2(Mst2) + pow2(
        MuSUSY))*pow4(Mst2) + 3*pow4(Mst1)*(27*pow2(Mst2)*pow2(MuSUSY) + 7*
        pow4(Mst2)) - 29*pow8(Mst2))) + Dmglst1*Mst1*(2*pow4(Msq)*(3*pow2(Mst1)
        *(8523*pow2(Mst2) - 8188*pow2(MuSUSY)) + 91136*pow4(Mst1) - 290*pow4(
        Mst2)) - 120*pow2(Msq)*(3*(37*pow2(Mst2) - 96*pow2(MuSUSY))*pow4(Mst1)
        - 12*pow2(Mst1)*(11*pow2(Mst2)*pow2(MuSUSY) + 4*pow4(Mst2)) + 4*pow6(
        Mst2)) + 15*(-2*pow2(Mst1)*(51*pow2(Mst2) + pow2(MuSUSY))*pow4(Mst2) +
        3*pow4(Mst1)*(159*pow2(Mst2)*pow2(MuSUSY) + 35*pow4(Mst2)) - 29*pow8(
        Mst2)))) - 1080*Mst1*pow2(log(pow2(Msq)/pow2(Mst1)))*(6*Mst1*pow2(
        Dmglst1)*(pow2(Mst1) - pow2(Mst2))*pow4(Mst2) - 6*(Dmglst1 + Mst1)*log(
        pow2(Mst2)/pow2(Mst1))*pow4(Msq)*(3*pow2(Mst1)*pow2(MuSUSY) + pow4(
        Mst2)) + Dmglst1*(12*pow2(Mst1)*pow2(MuSUSY)*pow4(Msq) + 4*pow2(Msq)*(
        pow2(Mst1) - pow2(Mst2))*pow4(Mst2) + 7*pow4(Mst1)*pow4(Mst2) - 6*pow2(
        Mst1)*pow6(Mst2) - pow8(Mst2)) + Mst1*(12*pow2(Mst1)*pow2(MuSUSY)*pow4(
        Msq) + 4*pow2(Msq)*(pow2(Mst1) - pow2(Mst2))*pow4(Mst2) + 3*pow4(Mst1)*
        pow4(Mst2) - 2*pow2(Mst1)*pow6(Mst2) - pow8(Mst2))) - 5070*Dmglst1*
        Mst1*pow8(Mst2) - 5070*pow2(Mst1)*pow8(Mst2) + 180*log(pow2(Msq)/pow2(
        Mst1))*(36*pow2(Mst1)*(12*pow2(Dmglst1)*(5*pow2(Mst1) + pow2(Mst2)) +
        pow2(Mst1)*(6*pow2(Mst1) + 4*pow2(Mst2) + pow2(MuSUSY)) + 3*Dmglst1*
        Mst1*(10*pow2(Mst1) + 4*pow2(Mst2) + 3*pow2(MuSUSY)))*pow2(log(pow2(
        Mst2)/pow2(Mst1)))*pow4(Msq) + 12*pow2(Dmglst1)*(24*pow2(Msq)*pow2(
        MuSUSY)*pow4(Mst1) + pow4(Msq)*(39*pow2(Mst1)*pow2(Mst2) + 195*pow4(
        Mst1) - 28*pow4(Mst2)) + pow4(Mst1)*(30*pow2(Mst2)*pow2(MuSUSY) + pow4(
        Mst2)) - pow2(Mst1)*pow6(Mst2)) - pow2(Mst1)*(pow2(Mst1)*(4*pow2(Mst2)
        + pow2(MuSUSY))*pow4(Mst2) - 12*pow4(Mst1)*(3*pow2(Mst2)*pow2(MuSUSY) +
        pow4(Mst2)) + 18*pow4(Msq)*(-2*pow2(Mst1)*(pow2(Mst2) - 10*pow2(MuSUSY)
        ) + 14*pow4(Mst1) + 11*pow4(Mst2)) - 4*pow2(Msq)*(24*pow2(MuSUSY)*pow4(
        Mst1) + pow2(Mst1)*(24*pow2(Mst2)*pow2(MuSUSY) - 7*pow4(Mst2)) + 7*
        pow6(Mst2)) + 8*pow8(Mst2)) - Dmglst1*Mst1*(pow2(Mst1)*(12*pow2(Mst2) +
        pow2(MuSUSY))*pow4(Mst2) - 20*pow4(Mst1)*(9*pow2(Mst2)*pow2(MuSUSY) +
        pow4(Mst2)) + 2*pow4(Msq)*(-18*pow2(Mst1)*(7*pow2(Mst2) - 4*pow2(
        MuSUSY)) + 90*pow4(Mst1) + 251*pow4(Mst2)) - 4*pow2(Msq)*(72*pow2(
        MuSUSY)*pow4(Mst1) + pow2(Mst1)*(72*pow2(Mst2)*pow2(MuSUSY) - 7*pow4(
        Mst2)) + 7*pow6(Mst2)) + 8*pow8(Mst2)) - 6*Mst1*log(pow2(Mst2)/pow2(
        Mst1))*(6*Mst1*pow2(Dmglst1)*(5*pow2(Mst1)*pow2(Mst2)*(pow2(Mst2) +
        pow2(MuSUSY)) + 32*(7*pow2(Mst1) + pow2(Mst2))*pow4(Msq) + 4*pow2(Msq)*
        (3*pow2(Mst1)*pow2(MuSUSY) + pow4(Mst2)) + pow6(Mst2)) + Mst1*(3*pow2(
        Mst2)*(pow2(Mst2) + pow2(MuSUSY))*pow4(Mst1) + 6*pow4(Msq)*(pow2(Mst1)*
        (4*pow2(Mst2) - 17*pow2(MuSUSY)) + 17*pow4(Mst1) - 2*pow4(Mst2)) + 2*
        pow2(Mst1)*pow6(Mst2) + 4*pow2(Msq)*(2*pow2(Mst1)*pow2(Mst2)*(pow2(
        Mst2) + pow2(MuSUSY)) + 6*pow2(MuSUSY)*pow4(Mst1) + pow6(Mst2)) + pow8(
        Mst2)) + Dmglst1*(15*pow2(Mst2)*(pow2(Mst2) + pow2(MuSUSY))*pow4(Mst1)
        + 6*pow4(Msq)*(pow2(Mst1)*(20*pow2(Mst2) - pow2(MuSUSY)) + 97*pow4(
        Mst1) - 2*pow4(Mst2)) + 6*pow2(Mst1)*pow6(Mst2) + 4*pow2(Msq)*(6*pow2(
        Mst1)*pow2(Mst2)*(pow2(Mst2) + pow2(MuSUSY)) + 18*pow2(MuSUSY)*pow4(
        Mst1) + pow6(Mst2)) + pow8(Mst2))))))/(243.*Mst1*pow4(Msq)*pow4(Mst2))
        + pow2(Mt)*pow2(s2t)*((40*(1 + 2*log(pow2(Msq)/pow2(Mst1)))*pow2(Msq))/
        3. + (2*(-98250*Dmglst1 - 22642*Mst1 + 1500*(21*Dmglst1 + 5*Mst1)*log(
        pow2(Mst2)/pow2(Mst1)) - 15*log(pow2(Msq)/pow2(Mst1))*(6000*Dmglst1 +
        1282*Mst1 - 75*Mst1*log(pow2(Mst2)/pow2(Mst1))) + 900*Mst1*pow2(log(
        pow2(Msq)/pow2(Mst1))))*pow3(Mst1))/(2025.*pow2(Msq)) - (pow2(Dmglst1)*
        (-486235 + 864*B4 + 432*DN + 110202*log(pow2(Mst2)/pow2(Mst1)) + 4860*
        log(pow2(Msq)/pow2(Mst1))*(1 + 2*log(pow2(Mst2)/pow2(Mst1))) - 9720*
        pow2(log(pow2(Msq)/pow2(Mst1))) - 29124*pow2(log(pow2(Mst2)/pow2(Mst1))
        ) - 1440*pow3(log(pow2(Mst2)/pow2(Mst1)))))/486. - (pow2(Mst1)*(
        32565277500*Dmglst1*Mst1*pow2(Msq) + 2235759750*pow2(Dmglst1)*pow2(Msq)
        + 29787188414*pow2(Msq)*pow2(Mst1) + 14586075000*pow2(Dmglst1)*pow2(
        Mst2) - 33562974291*pow2(Msq)*pow2(Mst2) + 111132000*DN*pow2(Msq)*pow2(
        Mst2) + 1260*log(pow2(Mst2)/pow2(Mst1))*(90799450*Dmglst1*Mst1*pow2(
        Msq) + 3675*pow2(Dmglst1)*(863*pow2(Msq) - 1620*pow2(Mst2)) + 9*pow2(
        Msq)*(4731149*pow2(Mst1) + 3302894*pow2(Mst2))) + 138915000*log(pow2(
        Msq)/pow2(Mst1))*(2*log(pow2(Mst2)/pow2(Mst1))*pow2(Msq)*(42*Dmglst1*
        Mst1 + 81*pow2(Dmglst1) + 6*pow2(Mst1) + 16*pow2(Mst2)) - 3*(18*
        Dmglst1*Mst1*pow2(Msq) + 5*pow2(Dmglst1)*(17*pow2(Msq) - 8*pow2(Mst2))
        + pow2(Msq)*(-4*pow2(Msq) - 5*pow2(Mst1) + 10*pow2(Mst2)))) +
        416745000*pow2(Msq)*(3*(2*Dmglst1*Mst1 + pow2(Dmglst1) + pow2(Mst1) -
        2*pow2(Mst2)) + log(pow2(Mst2)/pow2(Mst1))*pow2(Mst2))*pow2(log(pow2(
        Msq)/pow2(Mst1))) + 132300*pow2(Msq)*(379120*Dmglst1*Mst1 + 872760*
        pow2(Dmglst1) + 46219*pow2(Mst1) - 96621*pow2(Mst2))*pow2(log(pow2(
        Mst2)/pow2(Mst1))) - 18522000*pow2(Msq)*(908*Dmglst1*Mst1 + 1560*pow2(
        Dmglst1) + 223*pow2(Mst1) + 60*pow2(Mst2))*pow3(log(pow2(Mst2)/pow2(
        Mst1))) + 833490000*pow4(Msq)))/(1.250235e8*pow2(Msq)*pow2(Mst2)) +
        Mst1*((Dmglst1*(156781 - 864*B4 - 432*DN - 131946*log(pow2(Mst2)/pow2(
        Mst1)) - 540*log(pow2(Msq)/pow2(Mst1))*(-19 + 18*log(pow2(Mst2)/pow2(
        Mst1))) + 9720*pow2(log(pow2(Msq)/pow2(Mst1))) + 36900*pow2(log(pow2(
        Mst2)/pow2(Mst1))) + 1440*pow3(log(pow2(Mst2)/pow2(Mst1)))))/243. - (
        8575*Mst1*pow2(Dmglst1)*pow2(MuSUSY)*(2117 - 1074*log(pow2(Mst2)/pow2(
        Mst1)) + 180*(-(log(pow2(Msq)/pow2(Mst1))*(3 + 2*log(pow2(Mst2)/pow2(
        Mst1)))) + pow2(log(pow2(Mst2)/pow2(Mst1))))) + 2*(334425*Dmglst1 +
        219004*Mst1 - 84*(3675*Dmglst1 + 934*Mst1)*log(pow2(Mst2)/pow2(Mst1)) +
        210*log(pow2(Msq)/pow2(Mst1))*(1470*Dmglst1 + 991*Mst1 + 231*Mst1*log(
        pow2(Mst2)/pow2(Mst1))) + 22050*Mst1*pow2(log(pow2(Msq)/pow2(Mst1))) -
        70560*Mst1*pow2(log(pow2(Mst2)/pow2(Mst1))))*pow4(Mst2))/(185220.*pow4(
        Msq))) - (pow2(Mst2)*(15435000*pow2(Dmglst1)*pow2(Msq)*pow2(Mst1) +
        5016375*pow2(Dmglst1)*pow2(Mst1)*pow2(Mst2) - 14016352*pow2(Msq)*pow2(
        Mst1)*pow2(Mst2) + 308700*pow2(Msq)*(18*Dmglst1*Mst1*pow2(Msq) + 9*
        pow2(Dmglst1)*pow2(Msq) - 2*pow2(Mst1)*(-191*pow2(Msq) + 6*pow2(Mst1) +
        9*pow2(Mst2)))*pow2(log(pow2(Mst2)/pow2(Mst1))) - 43218000*Dmglst1*
        pow2(Msq)*pow3(Mst1) - 32516400*Dmglst1*Mst1*pow4(Msq) + 485379300*
        pow2(Dmglst1)*pow4(Msq) + 38596075*pow2(Mst1)*pow4(Msq) - 2469600*DN*
        pow2(Mst1)*pow4(Msq) + 1852200*pow2(Mst2)*pow4(Msq) - 40542600*pow2(
        Mst1)*pow3(log(pow2(Mst2)/pow2(Mst1)))*pow4(Msq) + 68299875*pow2(
        Dmglst1)*pow4(Mst1) + 8481704*pow2(Msq)*pow4(Mst1) + 44100*pow2(log(
        pow2(Msq)/pow2(Mst1)))*(1260*Dmglst1*Mst1*pow4(Msq) + 630*pow2(Dmglst1)
        *pow4(Msq) - 210*log(pow2(Mst2)/pow2(Mst1))*pow2(Mst1)*pow4(Msq) +
        pow2(Mst1)*(56*pow2(Msq)*(2*pow2(Mst1) - pow2(Mst2)) + 630*pow4(Msq) +
        15*pow4(Mst1))) - 4373250*Dmglst1*pow5(Mst1) + 18522000*pow6(Msq) +
        20580*log(pow2(Mst2)/pow2(Mst1))*((-7330*pow2(Mst1) + 180*pow2(Mst2))*
        pow4(Msq) + 45*pow2(Dmglst1)*(-20*pow2(Msq)*pow2(Mst1) - 5*pow2(Mst1)*(
        9*pow2(Mst1) + pow2(Mst2)) + 28*pow4(Msq)) + pow2(Msq)*(836*pow2(Mst1)*
        pow2(Mst2) - 1336*pow4(Mst1)) + 30*Dmglst1*(-60*pow2(Msq)*pow3(Mst1) +
        32*Mst1*pow4(Msq) - 35*pow5(Mst1)) - 339*pow6(Mst1)) - 7835880*pow6(
        Mst1) - 420*log(pow2(Msq)/pow2(Mst1))*(21364*pow2(Msq)*pow2(Mst1)*pow2(
        Mst2) - 22050*pow2(Mst1)*pow4(Msq) - 44100*pow2(Mst1)*pow2(log(pow2(
        Mst2)/pow2(Mst1)))*pow4(Msq) + 30772*pow2(Msq)*pow4(Mst1) + 735*log(
        pow2(Mst2)/pow2(Mst1))*pow2(Mst1)*(2*pow2(Msq)*(8*pow2(Mst1) - 13*pow2(
        Mst2)) + 140*pow4(Msq) + 9*pow4(Mst1)) + 11025*pow2(Dmglst1)*(40*pow2(
        Msq)*pow2(Mst1) - pow2(Mst1)*pow2(Mst2) + 24*pow4(Msq) + 37*pow4(Mst1))
        + 7350*Dmglst1*(40*pow2(Msq)*pow3(Mst1) + 16*Mst1*pow4(Msq) + 27*pow5(
        Mst1)) - 88200*pow6(Msq) + 43935*pow6(Mst1))))/(2.7783e6*pow2(Mst1)*
        pow4(Msq)) + (-3*pow2(Mst1)*pow2(Mst2)*((40*OepS2 - 4941*S2)*pow2(Mst2)
        + 8*(136*OepS2 - 3969*S2)*pow2(MuSUSY)) + 2*((8*OepS2 - 999*S2)*pow2(
        Mst2) + 2*(-1480*OepS2 + 55863*S2)*pow2(MuSUSY))*pow4(Mst1) + 27*(-8*
        OepS2 + 81*S2)*pow6(Mst2) + 162*S2*log(pow2(Mst2)/pow2(Mst1))*(-2*(
        pow2(Mst2) - 370*pow2(MuSUSY))*pow4(Mst1) + 3*pow2(Mst1)*(136*pow2(
        Mst2)*pow2(MuSUSY) + 5*pow4(Mst2)) + 27*pow6(Mst2)))/(729.*pow4(Mst2))
        + (pow2(MuSUSY)*(166698000000*pow2(Dmglst1)*pow2(Msq)*pow2(Mst1)*pow2(
        Mst2) + 81496800000*Dmglst1*pow2(Msq)*pow2(Mst2)*pow3(Mst1) -
        258631175250*pow2(Dmglst1)*pow2(Mst1)*pow4(Msq) + 11557728000*B4*pow2(
        Dmglst1)*pow2(Mst1)*pow4(Msq) + 4889808000*DN*pow2(Dmglst1)*pow2(Mst1)*
        pow4(Msq) + 53161877328*pow2(Mst1)*pow2(Mst2)*pow4(Msq) - 16892064000*
        B4*pow2(Mst1)*pow2(Mst2)*pow4(Msq) + 444528000*DN*pow2(Mst1)*pow2(Mst2)
        *pow4(Msq) + 225024292500*Dmglst1*pow3(Mst1)*pow4(Msq) - 5334336000*B4*
        Dmglst1*pow3(Mst1)*pow4(Msq) - 4445280000*Dmglst1*DN*pow3(Mst1)*pow4(
        Msq) + 74088000*pow2(Mst1)*(2418*Dmglst1*Mst1 + 1873*pow2(Dmglst1) +
        56*(70*pow2(Mst1) + 31*pow2(Mst2)))*pow3(log(pow2(Mst2)/pow2(Mst1)))*
        pow4(Msq) + 22596840000*pow2(Msq)*pow2(Mst2)*pow4(Mst1) + 531403689547*
        pow4(Msq)*pow4(Mst1) - 16892064000*B4*pow4(Msq)*pow4(Mst1) + 444528000*
        DN*pow4(Msq)*pow4(Mst1) - 12387513600*pow2(Msq)*pow2(Mst1)*pow4(Mst2) -
        44915850000*Dmglst1*pow3(Mst1)*pow4(Mst2) - 7419770460*pow4(Mst1)*pow4(
        Mst2) + 31752000*pow2(log(pow2(Msq)/pow2(Mst1)))*(-315*pow2(Mst1)*pow2(
        Mst2)*pow4(Msq) + 210*log(pow2(Mst2)/pow2(Mst1))*pow2(Mst1)*(6*Dmglst1*
        Mst1 + 3*pow2(Dmglst1) + 4*(pow2(Mst1) + pow2(Mst2)))*pow4(Msq) + pow2(
        pow2(Mst1) - pow2(Mst2))*pow4(Mst2)) + 3500599320*pow2(Mst1)*pow6(Mst2)
        - 75600*log(pow2(Msq)/pow2(Mst1))*(176400*pow2(Dmglst1)*pow2(Msq)*pow2(
        Mst1)*(19*pow2(Msq) - 5*pow2(Mst2)) - 14700*Dmglst1*pow2(Mst2)*(40*
        pow2(Msq) + pow2(Mst2))*pow3(Mst1) + 308700*pow2(Mst1)*pow2(Mst2)*pow4(
        Msq) + 88200*pow2(Mst1)*(6*Dmglst1*Mst1 + 3*pow2(Dmglst1) + 12*pow2(
        Mst1) + 8*pow2(Mst2))*pow2(log(pow2(Mst2)/pow2(Mst1)))*pow4(Msq) -
        147000*pow2(Msq)*pow2(Mst2)*pow4(Mst1) - 330750*pow4(Msq)*pow4(Mst1) -
        41160*pow2(Msq)*pow2(Mst1)*pow4(Mst2) + 12479*pow4(Mst1)*pow4(Mst2) -
        210*log(pow2(Mst2)/pow2(Mst1))*(2520*pow2(Dmglst1)*pow2(Msq)*pow2(Mst1)
        *(8*pow2(Msq) - 3*pow2(Mst2)) - 112*pow2(Msq)*pow2(Mst1)*pow2(Mst2)*(
        11*pow2(Mst1) + pow2(Mst2)) + 70*pow4(Msq)*(94*pow2(Mst1)*pow2(Mst2) +
        121*pow4(Mst1)) + 140*Dmglst1*pow3(Mst1)*(-28*pow2(Msq)*pow2(Mst2) +
        108*pow4(Msq) + 5*pow4(Mst2)) - pow4(Mst2)*(104*pow2(Mst1)*pow2(Mst2) +
        27*pow4(Mst1) + 18*pow4(Mst2)) + 1680*pow2(Mst1)*pow6(Msq)) - 26218*
        pow2(Mst1)*pow6(Mst2) - 9571*pow8(Mst2)) + 1331874540*pow8(Mst2) -
        264600*pow2(log(pow2(Mst2)/pow2(Mst1)))*(210*pow2(Dmglst1)*pow2(Msq)*
        pow2(Mst1)*(8947*pow2(Msq) - 1080*pow2(Mst2)) + pow4(Msq)*(1487556*
        pow2(Mst1)*pow2(Mst2) + 2258059*pow4(Mst1)) + 140*Dmglst1*pow3(Mst1)*(-
        840*pow2(Msq)*pow2(Mst2) + 9491*pow4(Msq) + 150*pow4(Mst2)) - 840*pow2(
        Msq)*(104*pow2(Mst2)*pow4(Mst1) + 29*pow2(Mst1)*pow4(Mst2)) - 30*(321*
        pow4(Mst1)*pow4(Mst2) + 216*pow2(Mst1)*pow6(Mst2) + 32*pow8(Mst2))) +
        1260*log(pow2(Mst2)/pow2(Mst1))*(36750*pow2(Dmglst1)*pow2(Msq)*pow2(
        Mst1)*(14027*pow2(Msq) - 3024*pow2(Mst2)) + 4900*Dmglst1*pow3(Mst1)*(-
        6960*pow2(Msq)*pow2(Mst2) + 68551*pow4(Msq) + 5010*pow4(Mst2)) + 3*(19*
        pow4(Msq)*(7684376*pow2(Mst1)*pow2(Mst2) + 11467149*pow4(Mst1)) - 3920*
        pow2(Msq)*(1928*pow2(Mst2)*pow4(Mst1) + 353*pow2(Mst1)*pow4(Mst2)) +
        3528000*pow2(Mst1)*pow6(Msq) - 10*(5151*pow4(Mst1)*pow4(Mst2) + 126132*
        pow2(Mst1)*pow6(Mst2) + 31294*pow8(Mst2))))))/(1.000188e9*pow4(Msq)*
        pow4(Mst2)))))/3888.)/pow4(Mt)/pow2(Sbeta)*12.;

   return result;
}

/**
 * 	@return returns the susy log^1 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double himalaya::H5::calc_coef_at_as2_no_sm_logs_log1(){

   const double result =
      ((20956800*s2t*pow2(Dmglst1)*pow2(Mst2)*pow2(Sbeta)*pow3(Mst1)*pow3(Mt)*
        pow4(Msq) + 11145600*s2t*z2*pow2(Dmglst1)*pow2(Mst2)*pow2(Sbeta)*pow3(
        Mst1)*pow3(Mt)*pow4(Msq) + 302400*Cbeta*Dmglst1*MuSUSY*Sbeta*pow2(Mst2)
        *pow2(Mt)*pow2(s2t)*pow4(Msq)*pow4(Mst1) + 3196800*Cbeta*Dmglst1*
        MuSUSY*Sbeta*z2*pow2(Mst2)*pow2(Mt)*pow2(s2t)*pow4(Msq)*pow4(Mst1) +
        1540800*pow2(Dmglst1)*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*pow4(Msq)*pow4(
        Mst1) + 619200*z2*pow2(Dmglst1)*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*pow4(
        Msq)*pow4(Mst1) + 1296000*pow2(Mst2)*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*
        pow4(Msq)*pow4(Mst1) + 28800*z2*pow2(Mst2)*pow2(Mt)*pow2(MuSUSY)*pow2(
        s2t)*pow4(Msq)*pow4(Mst1) - 5479200*pow2(Dmglst1)*pow2(Mst2)*pow2(Mt)*
        pow2(s2t)*pow2(Sbeta)*pow4(Msq)*pow4(Mst1) - 6912000*z2*pow2(Dmglst1)*
        pow2(Mst2)*pow2(Mt)*pow2(s2t)*pow2(Sbeta)*pow4(Msq)*pow4(Mst1) -
        1540800*pow2(Dmglst1)*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*pow2(Sbeta)*pow4(
        Msq)*pow4(Mst1) - 619200*z2*pow2(Dmglst1)*pow2(Mt)*pow2(MuSUSY)*pow2(
        s2t)*pow2(Sbeta)*pow4(Msq)*pow4(Mst1) - 1296000*pow2(Mst2)*pow2(Mt)*
        pow2(MuSUSY)*pow2(s2t)*pow2(Sbeta)*pow4(Msq)*pow4(Mst1) - 28800*z2*
        pow2(Mst2)*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*pow2(Sbeta)*pow4(Msq)*pow4(
        Mst1) + 244800*Cbeta*MuSUSY*s2t*Sbeta*pow2(Dmglst1)*pow3(Mt)*pow4(Msq)*
        pow4(Mst1) + 12441600*Cbeta*MuSUSY*s2t*Sbeta*z2*pow2(Dmglst1)*pow3(Mt)*
        pow4(Msq)*pow4(Mst1) - 6753600*Cbeta*MuSUSY*s2t*Sbeta*pow2(Mst2)*pow3(
        Mt)*pow4(Msq)*pow4(Mst1) - 1382400*Cbeta*MuSUSY*s2t*Sbeta*z2*pow2(Mst2)
        *pow3(Mt)*pow4(Msq)*pow4(Mst1) + 11232000*Dmglst1*s2t*pow2(MuSUSY)*
        pow3(Mt)*pow4(Msq)*pow4(Mst1) - 3283200*Dmglst1*s2t*z2*pow2(MuSUSY)*
        pow3(Mt)*pow4(Msq)*pow4(Mst1) + 14745600*Dmglst1*s2t*pow2(Mst2)*pow2(
        Sbeta)*pow3(Mt)*pow4(Msq)*pow4(Mst1) + 11347200*Dmglst1*s2t*z2*pow2(
        Mst2)*pow2(Sbeta)*pow3(Mt)*pow4(Msq)*pow4(Mst1) - 11232000*Dmglst1*s2t*
        pow2(MuSUSY)*pow2(Sbeta)*pow3(Mt)*pow4(Msq)*pow4(Mst1) + 3283200*
        Dmglst1*s2t*z2*pow2(MuSUSY)*pow2(Sbeta)*pow3(Mt)*pow4(Msq)*pow4(Mst1) +
        10800*Cbeta*Mt*MuSUSY*Sbeta*pow2(Dmglst1)*pow2(Mst2)*pow3(s2t)*pow4(
        Msq)*pow4(Mst1) - 4032000*s2t*pow2(Dmglst1)*pow2(Msq)*pow2(Sbeta)*pow3(
        Mst1)*pow3(Mt)*pow4(Mst2) + 5060400*pow2(Dmglst1)*pow2(Mst1)*pow2(Mt)*
        pow2(s2t)*pow2(Sbeta)*pow4(Msq)*pow4(Mst2) + 691200*z2*pow2(Dmglst1)*
        pow2(Mst1)*pow2(Mt)*pow2(s2t)*pow2(Sbeta)*pow4(Msq)*pow4(Mst2) +
        7922400*Dmglst1*pow2(Mt)*pow2(s2t)*pow2(Sbeta)*pow3(Mst1)*pow4(Msq)*
        pow4(Mst2) + 1382400*Dmglst1*z2*pow2(Mt)*pow2(s2t)*pow2(Sbeta)*pow3(
        Mst1)*pow4(Msq)*pow4(Mst2) - 4009600*Mst1*s2t*pow2(Dmglst1)*pow2(Sbeta)
        *pow3(Mt)*pow4(Msq)*pow4(Mst2) - 2009600*Dmglst1*s2t*pow2(Mst1)*pow2(
        Sbeta)*pow3(Mt)*pow4(Msq)*pow4(Mst2) + 403200*Dmglst1*s2t*z2*pow2(Mst1)
        *pow2(Sbeta)*pow3(Mt)*pow4(Msq)*pow4(Mst2) + 662400*s2t*pow2(Sbeta)*
        pow3(Mst1)*pow3(Mt)*pow4(Msq)*pow4(Mst2) + 403200*s2t*z2*pow2(Sbeta)*
        pow3(Mst1)*pow3(Mt)*pow4(Msq)*pow4(Mst2) - 2215200*Mt*pow2(Dmglst1)*
        pow2(Sbeta)*pow3(Mst1)*pow3(s2t)*pow4(Msq)*pow4(Mst2) - 1555200*Mt*z2*
        pow2(Dmglst1)*pow2(Sbeta)*pow3(Mst1)*pow3(s2t)*pow4(Msq)*pow4(Mst2) -
        2304000*Dmglst1*s2t*pow2(Msq)*pow2(Sbeta)*pow3(Mt)*pow4(Mst1)*pow4(
        Mst2) + 3636000*pow2(Mt)*pow2(s2t)*pow2(Sbeta)*pow4(Msq)*pow4(Mst1)*
        pow4(Mst2) + 691200*z2*pow2(Mt)*pow2(s2t)*pow2(Sbeta)*pow4(Msq)*pow4(
        Mst1)*pow4(Mst2) + 561600*Cbeta*Mt*MuSUSY*Sbeta*pow3(s2t)*pow4(Msq)*
        pow4(Mst1)*pow4(Mst2) - 162000*Cbeta*Mt*MuSUSY*Sbeta*z2*pow3(s2t)*pow4(
        Msq)*pow4(Mst1)*pow4(Mst2) - 3009600*Dmglst1*Mt*pow2(Sbeta)*pow3(s2t)*
        pow4(Msq)*pow4(Mst1)*pow4(Mst2) - 1310400*Dmglst1*Mt*z2*pow2(Sbeta)*
        pow3(s2t)*pow4(Msq)*pow4(Mst1)*pow4(Mst2) + 2083360*pow2(Dmglst1)*pow2(
        Mst1)*pow2(Mst2)*pow2(Sbeta)*pow4(Msq)*pow4(Mt) + 1036800*z2*pow2(
        Dmglst1)*pow2(Mst1)*pow2(Mst2)*pow2(Sbeta)*pow4(Msq)*pow4(Mt) +
        2947200*Dmglst1*pow2(Mst2)*pow2(Sbeta)*pow3(Mst1)*pow4(Msq)*pow4(Mt) +
        2073600*Dmglst1*z2*pow2(Mst2)*pow2(Sbeta)*pow3(Mst1)*pow4(Msq)*pow4(Mt)
        + 2832000*Cbeta*Dmglst1*MuSUSY*Sbeta*pow2(Msq)*pow2(Mst2)*pow4(Mst1)*
        pow4(Mt) + 2160000*pow2(Dmglst1)*pow2(Msq)*pow2(Mst2)*pow2(Sbeta)*pow4(
        Mst1)*pow4(Mt) - 11531200*Cbeta*Dmglst1*MuSUSY*Sbeta*pow4(Msq)*pow4(
        Mst1)*pow4(Mt) - 11750400*Cbeta*Dmglst1*MuSUSY*Sbeta*z2*pow4(Msq)*pow4(
        Mst1)*pow4(Mt) + 3456000*pow2(MuSUSY)*pow4(Msq)*pow4(Mst1)*pow4(Mt) +
        921600*z2*pow2(MuSUSY)*pow4(Msq)*pow4(Mst1)*pow4(Mt) - 9313760*pow2(
        Dmglst1)*pow2(Sbeta)*pow4(Msq)*pow4(Mst1)*pow4(Mt) - 2160000*z2*pow2(
        Dmglst1)*pow2(Sbeta)*pow4(Msq)*pow4(Mst1)*pow4(Mt) + 201600*pow2(Mst2)*
        pow2(Sbeta)*pow4(Msq)*pow4(Mst1)*pow4(Mt) + 1065600*z2*pow2(Mst2)*pow2(
        Sbeta)*pow4(Msq)*pow4(Mst1)*pow4(Mt) - 3456000*pow2(MuSUSY)*pow2(Sbeta)
        *pow4(Msq)*pow4(Mst1)*pow4(Mt) - 921600*z2*pow2(MuSUSY)*pow2(Sbeta)*
        pow4(Msq)*pow4(Mst1)*pow4(Mt) + 165600*pow2(Dmglst1)*pow2(Msq)*pow2(
        Mst1)*pow2(Sbeta)*pow4(Mst2)*pow4(Mt) - 216000*Dmglst1*pow2(Msq)*pow2(
        Sbeta)*pow3(Mst1)*pow4(Mst2)*pow4(Mt) - 551200*Dmglst1*Mst1*pow2(Sbeta)
        *pow4(Msq)*pow4(Mst2)*pow4(Mt) - 615552*pow2(Dmglst1)*pow2(Sbeta)*pow4(
        Msq)*pow4(Mst2)*pow4(Mt) - 604400*pow2(Mst1)*pow2(Sbeta)*pow4(Msq)*
        pow4(Mst2)*pow4(Mt) + 597600*z2*pow2(Mst1)*pow2(Sbeta)*pow4(Msq)*pow4(
        Mst2)*pow4(Mt) + 172800*z3*pow2(Mst1)*pow2(Sbeta)*pow4(Msq)*pow4(Mst2)*
        pow4(Mt) + 546000*Cbeta*Dmglst1*MuSUSY*Sbeta*pow4(Mst1)*pow4(Mst2)*
        pow4(Mt) + 473400*pow2(Dmglst1)*pow2(Sbeta)*pow4(Mst1)*pow4(Mst2)*pow4(
        Mt) - 396000*pow2(Msq)*pow2(Sbeta)*pow4(Mst1)*pow4(Mst2)*pow4(Mt) -
        200400*pow2(Dmglst1)*pow2(Sbeta)*pow4(Msq)*pow4(Mst1)*pow4(Mst2)*pow4(
        s2t) - 77400*z2*pow2(Dmglst1)*pow2(Sbeta)*pow4(Msq)*pow4(Mst1)*pow4(
        Mst2)*pow4(s2t) - 10656000*s2t*pow2(Dmglst1)*pow2(Msq)*pow2(Mst2)*pow2(
        Sbeta)*pow3(Mt)*pow5(Mst1) + 820800*Cbeta*MuSUSY*Sbeta*pow2(Dmglst1)*
        pow2(Mt)*pow2(s2t)*pow4(Msq)*pow5(Mst1) + 4665600*Cbeta*MuSUSY*Sbeta*
        z2*pow2(Dmglst1)*pow2(Mt)*pow2(s2t)*pow4(Msq)*pow5(Mst1) + 43200*Cbeta*
        MuSUSY*Sbeta*pow2(Mst2)*pow2(Mt)*pow2(s2t)*pow4(Msq)*pow5(Mst1) +
        86400*Cbeta*MuSUSY*Sbeta*z2*pow2(Mst2)*pow2(Mt)*pow2(s2t)*pow4(Msq)*
        pow5(Mst1) - 1382400*Dmglst1*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*pow4(Msq)*
        pow5(Mst1) + 1238400*Dmglst1*z2*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*pow4(
        Msq)*pow5(Mst1) - 2980800*Dmglst1*pow2(Mst2)*pow2(Mt)*pow2(s2t)*pow2(
        Sbeta)*pow4(Msq)*pow5(Mst1) - 2764800*Dmglst1*z2*pow2(Mst2)*pow2(Mt)*
        pow2(s2t)*pow2(Sbeta)*pow4(Msq)*pow5(Mst1) + 1382400*Dmglst1*pow2(Mt)*
        pow2(MuSUSY)*pow2(s2t)*pow2(Sbeta)*pow4(Msq)*pow5(Mst1) - 1238400*
        Dmglst1*z2*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*pow2(Sbeta)*pow4(Msq)*pow5(
        Mst1) - 9705600*Cbeta*Dmglst1*MuSUSY*s2t*Sbeta*pow3(Mt)*pow4(Msq)*pow5(
        Mst1) + 2764800*Cbeta*Dmglst1*MuSUSY*s2t*Sbeta*z2*pow3(Mt)*pow4(Msq)*
        pow5(Mst1) + 9561600*s2t*pow2(MuSUSY)*pow3(Mt)*pow4(Msq)*pow5(Mst1) +
        864000*s2t*z2*pow2(MuSUSY)*pow3(Mt)*pow4(Msq)*pow5(Mst1) + 60729600*
        s2t*pow2(Dmglst1)*pow2(Sbeta)*pow3(Mt)*pow4(Msq)*pow5(Mst1) + 38563200*
        s2t*z2*pow2(Dmglst1)*pow2(Sbeta)*pow3(Mt)*pow4(Msq)*pow5(Mst1) +
        1843200*s2t*pow2(Mst2)*pow2(Sbeta)*pow3(Mt)*pow4(Msq)*pow5(Mst1) +
        3916800*s2t*z2*pow2(Mst2)*pow2(Sbeta)*pow3(Mt)*pow4(Msq)*pow5(Mst1) -
        9561600*s2t*pow2(MuSUSY)*pow2(Sbeta)*pow3(Mt)*pow4(Msq)*pow5(Mst1) -
        864000*s2t*z2*pow2(MuSUSY)*pow2(Sbeta)*pow3(Mt)*pow4(Msq)*pow5(Mst1) +
        1447200*Cbeta*Dmglst1*Mt*MuSUSY*Sbeta*pow2(Mst2)*pow3(s2t)*pow4(Msq)*
        pow5(Mst1) + 1900800*Mt*pow2(Dmglst1)*pow2(Mst2)*pow2(Sbeta)*pow3(s2t)*
        pow4(Msq)*pow5(Mst1) - 3204000*s2t*pow2(Dmglst1)*pow2(Sbeta)*pow3(Mt)*
        pow4(Mst2)*pow5(Mst1) - 192000*s2t*pow2(Msq)*pow2(Sbeta)*pow3(Mt)*pow4(
        Mst2)*pow5(Mst1) - 2419200*Mt*pow2(Sbeta)*pow3(s2t)*pow4(Msq)*pow4(
        Mst2)*pow5(Mst1) - 273600*Mt*z2*pow2(Sbeta)*pow3(s2t)*pow4(Msq)*pow4(
        Mst2)*pow5(Mst1) + 14688000*Cbeta*MuSUSY*Sbeta*pow2(Dmglst1)*pow2(Msq)*
        pow4(Mt)*pow5(Mst1) + 3672000*Cbeta*MuSUSY*Sbeta*pow2(Dmglst1)*pow2(
        Mst2)*pow4(Mt)*pow5(Mst1) + 720000*Cbeta*MuSUSY*Sbeta*pow2(Msq)*pow2(
        Mst2)*pow4(Mt)*pow5(Mst1) + 1440000*Dmglst1*pow2(Msq)*pow2(Mst2)*pow2(
        Sbeta)*pow4(Mt)*pow5(Mst1) - 2203200*Cbeta*MuSUSY*Sbeta*pow4(Msq)*pow4(
        Mt)*pow5(Mst1) - 4320000*Cbeta*MuSUSY*Sbeta*z2*pow4(Msq)*pow4(Mt)*pow5(
        Mst1) - 11539200*Dmglst1*pow2(Sbeta)*pow4(Msq)*pow4(Mt)*pow5(Mst1) -
        345600*Dmglst1*z2*pow2(Sbeta)*pow4(Msq)*pow4(Mt)*pow5(Mst1) + 234000*
        Cbeta*MuSUSY*Sbeta*pow4(Mst2)*pow4(Mt)*pow5(Mst1) + 114000*Dmglst1*
        pow2(Sbeta)*pow4(Mst2)*pow4(Mt)*pow5(Mst1) - 534600*Dmglst1*pow2(Sbeta)
        *pow4(Msq)*pow4(Mst2)*pow4(s2t)*pow5(Mst1) - 154800*Dmglst1*z2*pow2(
        Sbeta)*pow4(Msq)*pow4(Mst2)*pow4(s2t)*pow5(Mst1) + 3600*Sbeta*log(pow2(
        Msq)/pow2(Mst1))*pow3(Mt)*(-20*Cbeta*Mt*MuSUSY*pow4(Mst1)*(118*Dmglst1*
        pow4(Msq) + 30*Mst1*pow4(Msq) - 7*Dmglst1*pow4(Mst2) - 3*Mst1*pow4(
        Mst2)) - 60*log(pow2(Mst2)/pow2(Mst1))*pow2(Mst1)*pow4(Msq)*(-4*Cbeta*(
        7*Dmglst1 + 3*Mst1)*Mt*MuSUSY*pow2(Mst1) + Sbeta*(6*Mst1*pow2(Dmglst1)*
        (-5*Mst1*Mt + 20*s2t*pow2(Mst1) + 4*s2t*pow2(Mst2)) + 8*s2t*pow2(Mst2)*
        pow3(Mst1) - 5*Mt*pow4(Mst1) + Mt*pow4(Mst2) + 4*Mst1*s2t*pow4(Mst2) +
        4*Dmglst1*(6*s2t*pow2(Mst1)*pow2(Mst2) - 5*Mt*pow3(Mst1) + 15*s2t*pow4(
        Mst1) + s2t*pow4(Mst2)) + 12*s2t*pow5(Mst1))) + Sbeta*(5*pow2(Mst1)*(-
        8*pow2(Msq)*(4*Mt*pow2(Mst1) + Mt*pow2(Mst2) + 2*Mst1*s2t*pow2(Mst2) -
        2*s2t*pow3(Mst1))*pow4(Mst2) - pow4(Mst2)*(6*Mt*pow2(Mst1)*pow2(Mst2) +
        8*s2t*pow2(Mst2)*pow3(Mst1) + 19*Mt*pow4(Mst1) + 3*Mt*pow4(Mst2) + 4*
        Mst1*s2t*pow4(Mst2) - 12*s2t*pow5(Mst1)) + 6*pow4(Msq)*(-6*Mt*pow2(
        Mst1)*pow2(Mst2) + 8*s2t*pow2(Mst2)*pow3(Mst1) - 14*Mt*pow4(Mst1) - 7*
        Mt*pow4(Mst2) + 12*Mst1*s2t*pow4(Mst2) + 8*s2t*pow5(Mst1))) + 20*
        Dmglst1*Mst1*(4*Mst1*pow2(Msq)*(-3*Mst1*Mt + s2t*pow2(Mst1) - s2t*pow2(
        Mst2))*pow4(Mst2) + Mst1*pow4(Mst2)*(-3*Mst1*Mt*pow2(Mst2) - 6*s2t*
        pow2(Mst1)*pow2(Mst2) - 13*Mt*pow3(Mst1) + 7*s2t*pow4(Mst1) - s2t*pow4(
        Mst2)) + 2*pow4(Msq)*(-9*Mt*pow2(Mst1)*pow2(Mst2) + 42*s2t*pow2(Mst2)*
        pow3(Mst1) - 57*Mt*pow4(Mst1) - 2*Mt*pow4(Mst2) + 17*Mst1*s2t*pow4(
        Mst2) + 66*s2t*pow5(Mst1))) + 6*pow2(Dmglst1)*(-20*Mt*pow2(Msq)*pow2(
        Mst1)*pow4(Mst2) + 5*pow2(Mst1)*(-11*Mt*pow2(Mst1) - Mt*pow2(Mst2) - 4*
        Mst1*s2t*pow2(Mst2) + 4*s2t*pow3(Mst1))*pow4(Mst2) + 2*pow4(Msq)*(-15*
        Mt*pow2(Mst1)*pow2(Mst2) + 260*s2t*pow2(Mst2)*pow3(Mst1) - 385*Mt*pow4(
        Mst1) + Mt*pow4(Mst2) + 20*Mst1*s2t*pow4(Mst2) + 740*s2t*pow5(Mst1)))))
        + 216000*pow2(Mst2)*pow2(Mt)*pow2(s2t)*pow2(Sbeta)*pow4(Mst1)*pow6(Msq)
        - 108000*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*pow3(s2t)*pow4(Mst1)*pow6(
        Msq) - 432000*pow2(Mst1)*pow2(Mt)*pow2(s2t)*pow2(Sbeta)*pow4(Mst2)*
        pow6(Msq) - 432000*pow2(Mst1)*pow2(Mst2)*pow2(Sbeta)*pow4(Mt)*pow6(Msq)
        - 432000*pow2(Sbeta)*pow4(Mst2)*pow4(Mt)*pow6(Msq) + 27000*pow2(Sbeta)*
        pow4(Mst1)*pow4(Mst2)*pow4(s2t)*pow6(Msq) - 3168000*Dmglst1*s2t*pow2(
        Msq)*pow2(Mst2)*pow2(Sbeta)*pow3(Mt)*pow6(Mst1) - 2754000*Cbeta*
        Dmglst1*MuSUSY*Sbeta*pow2(Mt)*pow2(s2t)*pow4(Msq)*pow6(Mst1) + 3196800*
        Cbeta*Dmglst1*MuSUSY*Sbeta*z2*pow2(Mt)*pow2(s2t)*pow4(Msq)*pow6(Mst1) -
        70200*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*pow4(Msq)*pow6(Mst1) + 28800*z2*
        pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*pow4(Msq)*pow6(Mst1) - 1576800*pow2(
        Mst2)*pow2(Mt)*pow2(s2t)*pow2(Sbeta)*pow4(Msq)*pow6(Mst1) + 70200*pow2(
        Mt)*pow2(MuSUSY)*pow2(s2t)*pow2(Sbeta)*pow4(Msq)*pow6(Mst1) - 28800*z2*
        pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*pow2(Sbeta)*pow4(Msq)*pow6(Mst1) -
        3600000*Cbeta*MuSUSY*s2t*Sbeta*pow3(Mt)*pow4(Msq)*pow6(Mst1) - 1382400*
        Cbeta*MuSUSY*s2t*Sbeta*z2*pow3(Mt)*pow4(Msq)*pow6(Mst1) + 8769600*
        Dmglst1*s2t*pow2(Sbeta)*pow3(Mt)*pow4(Msq)*pow6(Mst1) + 19756800*
        Dmglst1*s2t*z2*pow2(Sbeta)*pow3(Mt)*pow4(Msq)*pow6(Mst1) + 683100*
        Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*pow3(s2t)*pow4(Msq)*pow6(Mst1) +
        1018800*Dmglst1*Mt*pow2(Mst2)*pow2(Sbeta)*pow3(s2t)*pow4(Msq)*pow6(
        Mst1) - 954000*Dmglst1*s2t*pow2(Sbeta)*pow3(Mt)*pow4(Mst2)*pow6(Mst1) +
        6000000*Cbeta*Dmglst1*MuSUSY*Sbeta*pow2(Msq)*pow4(Mt)*pow6(Mst1) +
        1500000*Cbeta*Dmglst1*MuSUSY*Sbeta*pow2(Mst2)*pow4(Mt)*pow6(Mst1) +
        360000*pow2(Msq)*pow2(Mst2)*pow2(Sbeta)*pow4(Mt)*pow6(Mst1) - 3720600*
        pow2(Sbeta)*pow4(Msq)*pow4(Mt)*pow6(Mst1) + 345600*z2*pow2(Sbeta)*pow4(
        Msq)*pow4(Mt)*pow6(Mst1) - 61500*pow2(Sbeta)*pow4(Mst2)*pow4(Mt)*pow6(
        Mst1) - 30375*pow2(Sbeta)*pow4(Msq)*pow4(Mst2)*pow4(s2t)*pow6(Mst1) -
        40500*z2*pow2(Sbeta)*pow4(Msq)*pow4(Mst2)*pow4(s2t)*pow6(Mst1) - 27000*
        pow2(Mst2)*pow2(Sbeta)*pow4(s2t)*pow6(Msq)*pow6(Mst1) - 528000*Dmglst1*
        s2t*pow2(Msq)*pow2(Mst1)*pow2(Sbeta)*pow3(Mt)*pow6(Mst2) - 468000*s2t*
        pow2(Dmglst1)*pow2(Sbeta)*pow3(Mst1)*pow3(Mt)*pow6(Mst2) - 528000*s2t*
        pow2(Msq)*pow2(Sbeta)*pow3(Mst1)*pow3(Mt)*pow6(Mst2) - 88800*Dmglst1*
        Mst1*pow2(Mt)*pow2(s2t)*pow2(Sbeta)*pow4(Msq)*pow6(Mst2) + 296400*pow2(
        Dmglst1)*pow2(Mt)*pow2(s2t)*pow2(Sbeta)*pow4(Msq)*pow6(Mst2) - 280800*
        pow2(Mst1)*pow2(Mt)*pow2(s2t)*pow2(Sbeta)*pow4(Msq)*pow6(Mst2) + 40800*
        Mst1*Mt*pow2(Dmglst1)*pow2(Sbeta)*pow3(s2t)*pow4(Msq)*pow6(Mst2) +
        2908800*Dmglst1*Mt*pow2(Mst1)*pow2(Sbeta)*pow3(s2t)*pow4(Msq)*pow6(
        Mst2) + 244800*Dmglst1*Mt*z2*pow2(Mst1)*pow2(Sbeta)*pow3(s2t)*pow4(Msq)
        *pow6(Mst2) + 2404800*Mt*pow2(Sbeta)*pow3(Mst1)*pow3(s2t)*pow4(Msq)*
        pow6(Mst2) + 244800*Mt*z2*pow2(Sbeta)*pow3(Mst1)*pow3(s2t)*pow4(Msq)*
        pow6(Mst2) - 468000*Dmglst1*s2t*pow2(Sbeta)*pow3(Mt)*pow4(Mst1)*pow6(
        Mst2) - 63000*pow2(Dmglst1)*pow2(Mst1)*pow2(Sbeta)*pow4(Mt)*pow6(Mst2)
        - 192000*pow2(Msq)*pow2(Mst1)*pow2(Sbeta)*pow4(Mt)*pow6(Mst2) - 126000*
        Dmglst1*pow2(Sbeta)*pow3(Mst1)*pow4(Mt)*pow6(Mst2) - 43200*pow2(Sbeta)*
        pow4(Msq)*pow4(Mt)*pow6(Mst2) - 63000*pow2(Sbeta)*pow4(Mst1)*pow4(Mt)*
        pow6(Mst2) + 200100*pow2(Dmglst1)*pow2(Mst1)*pow2(Sbeta)*pow4(Msq)*
        pow4(s2t)*pow6(Mst2) + 77400*z2*pow2(Dmglst1)*pow2(Mst1)*pow2(Sbeta)*
        pow4(Msq)*pow4(s2t)*pow6(Mst2) + 156600*Dmglst1*pow2(Sbeta)*pow3(Mst1)*
        pow4(Msq)*pow4(s2t)*pow6(Mst2) + 154800*Dmglst1*z2*pow2(Sbeta)*pow3(
        Mst1)*pow4(Msq)*pow4(s2t)*pow6(Mst2) - 445500*pow2(Sbeta)*pow4(Msq)*
        pow4(Mst1)*pow4(s2t)*pow6(Mst2) + 77400*z2*pow2(Sbeta)*pow4(Msq)*pow4(
        Mst1)*pow4(s2t)*pow6(Mst2) - 156000*s2t*pow2(Sbeta)*pow3(Mt)*pow5(Mst1)
        *pow6(Mst2) + 216000*pow2(Mt)*pow2(s2t)*pow2(Sbeta)*pow6(Msq)*pow6(
        Mst2) + 27000*pow2(Mst1)*pow2(Sbeta)*pow4(s2t)*pow6(Msq)*pow6(Mst2) +
        900*pow2(Mst1)*pow2(log(pow2(Mst2)/pow2(Mst1)))*pow4(Msq)*(pow4(Mst1)*(
        -24*pow2(Mt)*pow2(s2t)*(102*Cbeta*Dmglst1*MuSUSY*Sbeta + 231*pow2(
        MuSUSY)*(-1 + pow2(Sbeta)) + 16*pow2(Mst2)*pow2(Sbeta)) + 16*s2t*Sbeta*
        (81*Cbeta*MuSUSY + 968*Dmglst1*Sbeta)*pow3(Mt) + 16*Mt*Sbeta*(-73*
        Cbeta*MuSUSY + 32*Dmglst1*Sbeta)*pow2(Mst2)*pow3(s2t) + 48*pow2(Sbeta)*
        pow4(Mt) - 27*pow2(Sbeta)*pow4(Mst2)*pow4(s2t)) - 2*pow3(Mst1)*(-8*
        pow2(Mt)*pow2(s2t)*(162*Cbeta*MuSUSY*Sbeta*pow2(Dmglst1) - 165*Cbeta*
        MuSUSY*Sbeta*pow2(Mst2) + Dmglst1*(-283*pow2(MuSUSY)*(-1 + pow2(Sbeta))
        - 192*pow2(Mst2)*pow2(Sbeta))) + 8*s2t*(-368*Cbeta*Dmglst1*MuSUSY*Sbeta
        + 334*pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 37*(51*pow2(Dmglst1) + 4*pow2(
        Mst2))*pow2(Sbeta))*pow3(Mt) + 4*Mt*Sbeta*pow2(Mst2)*(128*Cbeta*
        Dmglst1*MuSUSY - 53*Sbeta*pow2(Mst2))*pow3(s2t) + 64*Sbeta*(17*Cbeta*
        MuSUSY + 6*Dmglst1*Sbeta)*pow4(Mt) + 27*Dmglst1*pow2(Sbeta)*pow4(Mst2)*
        pow4(s2t)) - pow2(Mst1)*(4*pow2(Mt)*pow2(s2t)*(228*Cbeta*Dmglst1*
        MuSUSY*Sbeta*pow2(Mst2) + pow2(Mst2)*(802*pow2(MuSUSY)*(-1 + pow2(
        Sbeta)) - 57*pow2(Mst2)*pow2(Sbeta)) + 2*pow2(Dmglst1)*(283*pow2(
        MuSUSY)*(-1 + pow2(Sbeta)) + 768*pow2(Mst2)*pow2(Sbeta))) - 16*s2t*(
        760*Cbeta*MuSUSY*Sbeta*pow2(Dmglst1) + 33*Cbeta*MuSUSY*Sbeta*pow2(Mst2)
        + 2*Dmglst1*(-95*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 231*pow2(Mst2)*pow2(
        Sbeta)))*pow3(Mt) + 4*Mt*Sbeta*pow2(Mst2)*(38*Dmglst1*Sbeta*pow2(Mst2)
        + Cbeta*MuSUSY*(128*pow2(Dmglst1) + 319*pow2(Mst2)))*pow3(s2t) + 16*(
        450*Cbeta*Dmglst1*MuSUSY*Sbeta + 16*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + (
        129*pow2(Dmglst1) - 34*pow2(Mst2))*pow2(Sbeta))*pow4(Mt) + 3*(9*pow2(
        Dmglst1) - 79*pow2(Mst2))*pow2(Sbeta)*pow4(Mst2)*pow4(s2t)) + 2*Mst1*
        pow2(Mst2)*pow2(Sbeta)*(-96*s2t*pow2(Mst2)*pow3(Mt) + 24*pow2(Dmglst1)*
        (157*s2t*pow3(Mt) - 18*Mt*pow2(Mst2)*pow3(s2t)) + 228*Mt*pow3(s2t)*
        pow4(Mst2) + Dmglst1*(64*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 528*pow4(Mt) +
        155*pow4(Mst2)*pow4(s2t))) + 16*Mt*s2t*Sbeta*(-261*Cbeta*Mt*MuSUSY*s2t
        + 190*Sbeta*pow2(Mt) + 32*Sbeta*pow2(Mst2)*pow2(s2t))*pow5(Mst1) +
        pow2(Mst2)*pow2(Sbeta)*(-492*pow2(Mt)*pow2(s2t)*pow4(Mst2) + 24*
        Dmglst1*(-8*s2t*pow2(Mst2)*pow3(Mt) + 19*Mt*pow3(s2t)*pow4(Mst2)) -
        212*pow2(Mst2)*pow4(Mt) + pow2(Dmglst1)*(64*pow2(Mst2)*pow2(Mt)*pow2(
        s2t) + 528*pow4(Mt) + 155*pow4(Mst2)*pow4(s2t)) + 82*pow4(s2t)*pow6(
        Mst2))) - 288000*s2t*pow2(Msq)*pow2(Mst2)*pow2(Sbeta)*pow3(Mt)*pow7(
        Mst1) + 97200*Cbeta*MuSUSY*Sbeta*pow2(Mt)*pow2(s2t)*pow4(Msq)*pow7(
        Mst1) + 86400*Cbeta*MuSUSY*Sbeta*z2*pow2(Mt)*pow2(s2t)*pow4(Msq)*pow7(
        Mst1) - 2692800*s2t*pow2(Sbeta)*pow3(Mt)*pow4(Msq)*pow7(Mst1) +
        4032000*s2t*z2*pow2(Sbeta)*pow3(Mt)*pow4(Msq)*pow7(Mst1) - 18000*Mt*
        pow2(Mst2)*pow2(Sbeta)*pow3(s2t)*pow4(Msq)*pow7(Mst1) - 18000*s2t*pow2(
        Sbeta)*pow3(Mt)*pow4(Mst2)*pow7(Mst1) + 1008000*Cbeta*MuSUSY*Sbeta*
        pow2(Msq)*pow4(Mt)*pow7(Mst1) + 252000*Cbeta*MuSUSY*Sbeta*pow2(Mst2)*
        pow4(Mt)*pow7(Mst1) - 30*log(pow2(Mst2)/pow2(Mst1))*(-600*Sbeta*pow2(
        Mst1)*pow2(Mst2)*pow3(Mt)*(-4*Cbeta*Mt*MuSUSY*pow2(Mst1)*(36*Mst1*pow2(
        Dmglst1) + 22*Dmglst1*pow2(Mst1) + 7*Dmglst1*pow2(Mst2) + 3*Mst1*pow2(
        Mst2) + 6*pow3(Mst1)) + Sbeta*pow2(Mst2)*(6*Mt*pow2(Mst1)*pow2(Mst2) +
        6*(Mt + 4*Mst1*s2t)*pow2(Dmglst1)*(5*pow2(Mst1) + pow2(Mst2)) + 8*s2t*
        pow2(Mst2)*pow3(Mst1) + 5*Mt*pow4(Mst1) + 3*Mt*pow4(Mst2) + 4*Mst1*s2t*
        pow4(Mst2) + 4*Dmglst1*(3*Mst1*Mt*pow2(Mst2) + 6*s2t*pow2(Mst1)*pow2(
        Mst2) + 5*Mt*pow3(Mst1) + 15*s2t*pow4(Mst1) + s2t*pow4(Mst2)) + 12*s2t*
        pow5(Mst1))) + 2400*Sbeta*pow2(Msq)*pow2(Mst1)*pow3(Mt)*(4*Cbeta*Mt*
        MuSUSY*pow2(Mst1)*(36*Mst1*pow2(Dmglst1) + 22*Dmglst1*pow2(Mst1) + 7*
        Dmglst1*pow2(Mst2) + 3*Mst1*pow2(Mst2) + 6*pow3(Mst1)) - Sbeta*pow2(
        Mst2)*(3*Mt*pow2(Mst1)*pow2(Mst2) + 8*s2t*pow2(Mst2)*pow3(Mst1) + 3*
        pow2(Dmglst1)*((Mt + 8*Mst1*s2t)*pow2(Mst2) + 40*s2t*pow3(Mst1)) + 2*
        Mt*pow4(Mst2) + 4*Mst1*s2t*pow4(Mst2) + Dmglst1*(6*Mst1*Mt*pow2(Mst2) +
        24*s2t*pow2(Mst1)*pow2(Mst2) + 60*s2t*pow4(Mst1) + 4*s2t*pow4(Mst2)) +
        12*s2t*pow5(Mst1))) + 1800*pow2(Mst1)*pow2(s2t)*pow6(Msq)*(pow2(Mst1)*(
        8*pow2(Mt)*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + pow2(s2t)*pow2(Sbeta)*
        pow4(Mst2)) - pow2(s2t)*pow2(Sbeta)*pow6(Mst2)) + pow4(Msq)*(160*
        Dmglst1*Mst1*pow2(Sbeta)*pow4(Mst2)*(-48*pow2(Mst2)*pow2(Mt)*pow2(s2t)
        + 9*Dmglst1*s2t*pow3(Mt) + 22*pow4(Mt) + 12*pow4(Mst2)*pow4(s2t)) + 48*
        pow2(Dmglst1)*pow2(Sbeta)*pow4(Mst2)*(-80*pow2(Mst2)*pow2(Mt)*pow2(s2t)
        - 3*pow4(Mt) + 20*pow4(Mst2)*pow4(s2t)) - 10*pow4(Mst1)*(24*pow2(Mt)*
        pow2(s2t)*(1140*Cbeta*Dmglst1*MuSUSY*Sbeta*pow2(Mst2) + pow2(Mst2)*(
        727*pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 415*pow2(Mst2)*pow2(Sbeta)) +
        pow2(Dmglst1)*(671*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 901*pow2(Mst2)*
        pow2(Sbeta))) - 16*s2t*(2022*Cbeta*MuSUSY*Sbeta*pow2(Dmglst1) - 945*
        Cbeta*MuSUSY*Sbeta*pow2(Mst2) + Dmglst1*(-3864*pow2(MuSUSY)*(-1 + pow2(
        Sbeta)) + 4418*pow2(Mst2)*pow2(Sbeta)))*pow3(Mt) + 24*Mt*Sbeta*pow2(
        Mst2)*(-116*Dmglst1*Sbeta*pow2(Mst2) + Cbeta*MuSUSY*(173*pow2(Dmglst1)
        + 206*pow2(Mst2)))*pow3(s2t) + 8*(8876*Cbeta*Dmglst1*MuSUSY*Sbeta +
        960*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + (2153*pow2(Dmglst1) - 813*pow2(
        Mst2))*pow2(Sbeta))*pow4(Mt) - 3*(53*pow2(Dmglst1) + 97*pow2(Mst2))*
        pow2(Sbeta)*pow4(Mst2)*pow4(s2t)) + 20*pow2(Mst2)*pow2(Sbeta)*pow3(
        Mst1)*(392*pow2(Dmglst1)*(91*s2t*pow3(Mt) - 12*Mt*pow2(Mst2)*pow3(s2t))
        + 96*(17*s2t*pow2(Mst2)*pow3(Mt) + 33*Mt*pow3(s2t)*pow4(Mst2)) +
        Dmglst1*(7368*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 6696*pow4(Mt) + 783*pow4(
        Mst2)*pow4(s2t))) - 20*(-24*pow2(Mt)*pow2(s2t)*(1281*Cbeta*MuSUSY*
        Sbeta*pow2(Dmglst1) - 690*Cbeta*MuSUSY*Sbeta*pow2(Mst2) + Dmglst1*(-
        401*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 263*pow2(Mst2)*pow2(Sbeta))) + 8*
        s2t*(3324*Cbeta*Dmglst1*MuSUSY*Sbeta + 4344*pow2(MuSUSY)*(-1 + pow2(
        Sbeta)) - (24875*pow2(Dmglst1) + 1464*pow2(Mst2))*pow2(Sbeta))*pow3(Mt)
        + 24*Mt*Sbeta*pow2(Mst2)*(38*Cbeta*Dmglst1*MuSUSY + 231*Sbeta*pow2(
        Dmglst1) - 98*Sbeta*pow2(Mst2))*pow3(s2t) + 16*Sbeta*(801*Cbeta*MuSUSY
        + 200*Dmglst1*Sbeta)*pow4(Mt) + 651*Dmglst1*pow2(Sbeta)*pow4(Mst2)*
        pow4(s2t))*pow5(Mst1) + 5*(-24*pow2(Mt)*pow2(s2t)*(1218*Cbeta*Dmglst1*
        MuSUSY*Sbeta + 1471*pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 836*pow2(Mst2)*
        pow2(Sbeta)) + 32*s2t*Sbeta*(-2199*Cbeta*MuSUSY + 12656*Dmglst1*Sbeta)*
        pow3(Mt) - 12*Mt*Sbeta*(17*Cbeta*MuSUSY + 708*Dmglst1*Sbeta)*pow2(Mst2)
        *pow3(s2t) + 888*pow2(Sbeta)*pow4(Mt) - 2421*pow2(Sbeta)*pow4(Mst2)*
        pow4(s2t))*pow6(Mst1) + 10*pow2(Mst1)*pow2(Mst2)*pow2(Sbeta)*(32*
        Dmglst1*(103*s2t*pow2(Mst2)*pow3(Mt) + 198*Mt*pow3(s2t)*pow4(Mst2)) +
        pow2(Dmglst1)*(5832*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 6824*pow4(Mt) +
        783*pow4(Mst2)*pow4(s2t)) + 15*(-160*pow2(Mt)*pow2(s2t)*pow4(Mst2) +
        72*pow2(Mst2)*pow4(Mt) + 63*pow4(s2t)*pow6(Mst2))) + 240*Mt*s2t*Sbeta*(
        -1317*Cbeta*Mt*MuSUSY*s2t + 1724*Sbeta*pow2(Mt) - 21*Sbeta*pow2(Mst2)*
        pow2(s2t))*pow7(Mst1))) - 78000*Dmglst1*s2t*pow2(Mst1)*pow2(Sbeta)*
        pow3(Mt)*pow8(Mst2) - 78000*s2t*pow2(Sbeta)*pow3(Mst1)*pow3(Mt)*pow8(
        Mst2) + 21600*pow2(Mt)*pow2(s2t)*pow2(Sbeta)*pow4(Msq)*pow8(Mst2) +
        10800*Cbeta*Mt*MuSUSY*Sbeta*pow3(s2t)*pow4(Msq)*pow8(Mst2) - 40500*
        pow2(Mst1)*pow2(Sbeta)*pow4(Mt)*pow8(Mst2) + 16200*Dmglst1*Mst1*pow2(
        Sbeta)*pow4(Msq)*pow4(s2t)*pow8(Mst2) - 2400*pow2(Dmglst1)*pow2(Sbeta)*
        pow4(Msq)*pow4(s2t)*pow8(Mst2) + 307800*pow2(Mst1)*pow2(Sbeta)*pow4(
        Msq)*pow4(s2t)*pow8(Mst2) - 36900*z2*pow2(Mst1)*pow2(Sbeta)*pow4(Msq)*
        pow4(s2t)*pow8(Mst2) - 27000*pow2(Sbeta)*pow4(s2t)*pow6(Msq)*pow8(Mst2)
        )/(16200.*pow2(Mst1)*pow4(Msq)*pow4(Mst2)))/pow4(Mt)/pow2(Sbeta)*12.; 

   return result;
}

/**
 * 	@return returns the susy log^2 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double himalaya::H5::calc_coef_at_as2_no_sm_logs_log2(){

   const double result =
      ((263040*s2t*pow2(Dmglst1)*pow2(Mst2)*pow2(Sbeta)*pow3(Mst1)*pow3(Mt) -
        34560*Cbeta*Dmglst1*MuSUSY*Sbeta*pow2(Mst2)*pow2(Mt)*pow2(s2t)*pow4(
        Mst1) - 11520*pow2(Dmglst1)*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*pow4(Mst1)
        + 14760*pow2(Mst2)*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*pow4(Mst1) - 22440*
        pow2(Dmglst1)*pow2(Mst2)*pow2(Mt)*pow2(s2t)*pow2(Sbeta)*pow4(Mst1) +
        11520*pow2(Dmglst1)*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*pow2(Sbeta)*pow4(
        Mst1) - 14760*pow2(Mst2)*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*pow2(Sbeta)*
        pow4(Mst1) - 23040*Cbeta*MuSUSY*s2t*Sbeta*pow2(Dmglst1)*pow3(Mt)*pow4(
        Mst1) - 67920*Cbeta*MuSUSY*s2t*Sbeta*pow2(Mst2)*pow3(Mt)*pow4(Mst1) +
        170880*Dmglst1*s2t*pow2(MuSUSY)*pow3(Mt)*pow4(Mst1) + 136320*Dmglst1*
        s2t*pow2(Mst2)*pow2(Sbeta)*pow3(Mt)*pow4(Mst1) - 170880*Dmglst1*s2t*
        pow2(MuSUSY)*pow2(Sbeta)*pow3(Mt)*pow4(Mst1) + 22740*Cbeta*Mt*MuSUSY*
        Sbeta*pow2(Dmglst1)*pow2(Mst2)*pow3(s2t)*pow4(Mst1) + 33360*pow2(
        Dmglst1)*pow2(Mst1)*pow2(Mt)*pow2(s2t)*pow2(Sbeta)*pow4(Mst2) + 82080*
        Dmglst1*pow2(Mt)*pow2(s2t)*pow2(Sbeta)*pow3(Mst1)*pow4(Mst2) - 12160*
        Mst1*s2t*pow2(Dmglst1)*pow2(Sbeta)*pow3(Mt)*pow4(Mst2) + 36160*Dmglst1*
        s2t*pow2(Mst1)*pow2(Sbeta)*pow3(Mt)*pow4(Mst2) + 25920*s2t*pow2(Sbeta)*
        pow3(Mst1)*pow3(Mt)*pow4(Mst2) + 42960*pow2(Mt)*pow2(s2t)*pow2(Sbeta)*
        pow4(Mst1)*pow4(Mst2) + 6960*Cbeta*Mt*MuSUSY*Sbeta*pow3(s2t)*pow4(Mst1)
        *pow4(Mst2) - 19680*Dmglst1*Mt*pow2(Sbeta)*pow3(s2t)*pow4(Mst1)*pow4(
        Mst2) + 84400*pow2(Dmglst1)*pow2(Mst1)*pow2(Mst2)*pow2(Sbeta)*pow4(Mt)
        + 86880*Dmglst1*pow2(Mst2)*pow2(Sbeta)*pow3(Mst1)*pow4(Mt) - 172480*
        Cbeta*Dmglst1*MuSUSY*Sbeta*pow4(Mst1)*pow4(Mt) + 30720*pow2(MuSUSY)*
        pow4(Mst1)*pow4(Mt) - 57440*pow2(Dmglst1)*pow2(Sbeta)*pow4(Mst1)*pow4(
        Mt) + 28080*pow2(Mst2)*pow2(Sbeta)*pow4(Mst1)*pow4(Mt) - 30720*pow2(
        MuSUSY)*pow2(Sbeta)*pow4(Mst1)*pow4(Mt) + 29600*Dmglst1*Mst1*pow2(
        Sbeta)*pow4(Mst2)*pow4(Mt) - 9744*pow2(Dmglst1)*pow2(Sbeta)*pow4(Mst2)*
        pow4(Mt) + 29440*pow2(Mst1)*pow2(Sbeta)*pow4(Mst2)*pow4(Mt) + 7200*log(
        pow2(Msq)/pow2(Mst1))*pow2(Mst1)*pow2(Sbeta)*pow4(Mst2)*pow4(Mt) -
        10005*pow2(Dmglst1)*pow2(Sbeta)*pow4(Mst1)*pow4(Mst2)*pow4(s2t) -
        11520*Cbeta*MuSUSY*Sbeta*pow2(Mst2)*pow2(Mt)*pow2(s2t)*pow5(Mst1) -
        7680*Dmglst1*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*pow5(Mst1) - 29520*
        Dmglst1*pow2(Mst2)*pow2(Mt)*pow2(s2t)*pow2(Sbeta)*pow5(Mst1) + 7680*
        Dmglst1*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*pow2(Sbeta)*pow5(Mst1) - 76800*
        Cbeta*Dmglst1*MuSUSY*s2t*Sbeta*pow3(Mt)*pow5(Mst1) + 109440*s2t*pow2(
        MuSUSY)*pow3(Mt)*pow5(Mst1) + 416640*s2t*pow2(Dmglst1)*pow2(Sbeta)*
        pow3(Mt)*pow5(Mst1) + 17280*s2t*pow2(Mst2)*pow2(Sbeta)*pow3(Mt)*pow5(
        Mst1) - 109440*s2t*pow2(MuSUSY)*pow2(Sbeta)*pow3(Mt)*pow5(Mst1) +
        22440*Cbeta*Dmglst1*Mt*MuSUSY*Sbeta*pow2(Mst2)*pow3(s2t)*pow5(Mst1) -
        11520*Mt*pow2(Dmglst1)*pow2(Mst2)*pow2(Sbeta)*pow3(s2t)*pow5(Mst1) -
        19680*Mt*pow2(Sbeta)*pow3(s2t)*pow4(Mst2)*pow5(Mst1) - 43200*Cbeta*
        MuSUSY*Sbeta*pow4(Mt)*pow5(Mst1) - 18240*Dmglst1*pow2(Sbeta)*pow4(Mt)*
        pow5(Mst1) - 8490*Dmglst1*pow2(Sbeta)*pow4(Mst2)*pow4(s2t)*pow5(Mst1) -
        1920*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*pow6(Mst1) - 12840*pow2(Mst2)*
        pow2(Mt)*pow2(s2t)*pow2(Sbeta)*pow6(Mst1) + 1920*pow2(Mt)*pow2(MuSUSY)*
        pow2(s2t)*pow2(Sbeta)*pow6(Mst1) - 42240*Cbeta*MuSUSY*s2t*Sbeta*pow3(
        Mt)*pow6(Mst1) + 97920*Dmglst1*s2t*pow2(Sbeta)*pow3(Mt)*pow6(Mst1) +
        8340*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*pow3(s2t)*pow6(Mst1) - 11520*
        Dmglst1*Mt*pow2(Mst2)*pow2(Sbeta)*pow3(s2t)*pow6(Mst1) + 480*pow2(
        Sbeta)*pow4(Mt)*pow6(Mst1) - 345*pow2(Sbeta)*pow4(Mst2)*pow4(s2t)*pow6(
        Mst1) - 14160*Dmglst1*Mst1*pow2(Mt)*pow2(s2t)*pow2(Sbeta)*pow6(Mst2) +
        600*pow2(Dmglst1)*pow2(Mt)*pow2(s2t)*pow2(Sbeta)*pow6(Mst2) - 9000*
        pow2(Mst1)*pow2(Mt)*pow2(s2t)*pow2(Sbeta)*pow6(Mst2) + 11520*Mst1*Mt*
        pow2(Dmglst1)*pow2(Sbeta)*pow3(s2t)*pow6(Mst2) + 31200*Dmglst1*Mt*pow2(
        Mst1)*pow2(Sbeta)*pow3(s2t)*pow6(Mst2) + 23520*Mt*pow2(Sbeta)*pow3(
        Mst1)*pow3(s2t)*pow6(Mst2) + 4395*pow2(Dmglst1)*pow2(Mst1)*pow2(Sbeta)*
        pow4(s2t)*pow6(Mst2) + 1110*Dmglst1*pow2(Sbeta)*pow3(Mst1)*pow4(s2t)*
        pow6(Mst2) - 5325*pow2(Sbeta)*pow4(Mst1)*pow4(s2t)*pow6(Mst2) - 30*log(
        pow2(Mst2)/pow2(Mst1))*pow2(Mst1)*(4*pow4(Mst1)*(-2*pow2(Mt)*pow2(s2t)*
        (534*Cbeta*Dmglst1*MuSUSY*Sbeta + 173*pow2(MuSUSY)*(-1 + pow2(Sbeta)) -
        96*pow2(Mst2)*pow2(Sbeta)) + 16*s2t*Sbeta*(-36*Cbeta*MuSUSY + 115*
        Dmglst1*Sbeta)*pow3(Mt) + 172*pow2(Sbeta)*pow4(Mt) - 33*pow2(Sbeta)*
        pow4(Mst2)*pow4(s2t)) - 2*pow3(Mst1)*(8*pow2(Mt)*pow2(s2t)*(144*Cbeta*
        MuSUSY*Sbeta*pow2(Dmglst1) + 171*Cbeta*MuSUSY*Sbeta*pow2(Mst2) +
        Dmglst1*(91*pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 96*pow2(Mst2)*pow2(Sbeta)
        )) + 16*s2t*(144*Cbeta*Dmglst1*MuSUSY*Sbeta + 155*pow2(MuSUSY)*(-1 +
        pow2(Sbeta)) - 4*(115*pow2(Dmglst1) + 11*pow2(Mst2))*pow2(Sbeta))*pow3(
        Mt) - 292*Mt*pow2(Sbeta)*pow3(s2t)*pow4(Mst2) + 16*Sbeta*(29*Cbeta*
        MuSUSY - 54*Dmglst1*Sbeta)*pow4(Mt) + 91*Dmglst1*pow2(Sbeta)*pow4(Mst2)
        *pow4(s2t)) + pow2(Mst1)*(-8*pow2(Mt)*pow2(s2t)*(534*Cbeta*Dmglst1*
        MuSUSY*Sbeta*pow2(Mst2) + pow2(Dmglst1)*(91*pow2(MuSUSY)*(-1 + pow2(
        Sbeta)) - 96*pow2(Mst2)*pow2(Sbeta)) + pow2(Mst2)*(173*pow2(MuSUSY)*(-1
        + pow2(Sbeta)) - 89*pow2(Mst2)*pow2(Sbeta))) - 96*s2t*(24*Cbeta*MuSUSY*
        Sbeta*pow2(Dmglst1) + 8*Cbeta*MuSUSY*Sbeta*pow2(Mst2) + Dmglst1*(73*
        pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 44*pow2(Mst2)*pow2(Sbeta)))*pow3(Mt)
        + 8*Mt*Sbeta*(-66*Cbeta*MuSUSY + 137*Dmglst1*Sbeta)*pow3(s2t)*pow4(
        Mst2) - 32*(117*Cbeta*Dmglst1*MuSUSY*Sbeta + 16*pow2(MuSUSY)*(-1 +
        pow2(Sbeta)) - (49*pow2(Dmglst1) + 16*pow2(Mst2))*pow2(Sbeta))*pow4(Mt)
        + 91*(-pow2(Dmglst1) + pow2(Mst2))*pow2(Sbeta)*pow4(Mst2)*pow4(s2t)) +
        2*Mst1*pow2(Mst2)*pow2(Sbeta)*(4*Mt*s2t*pow2(Mst2)*(-60*pow2(Mt) + 41*
        pow2(Mst2)*pow2(s2t)) + 192*pow2(Dmglst1)*(11*s2t*pow3(Mt) + 2*Mt*pow2(
        Mst2)*pow3(s2t)) + Dmglst1*(384*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 512*
        pow4(Mt) + 91*pow4(Mst2)*pow4(s2t))) + 16*s2t*Sbeta*(-171*Cbeta*MuSUSY*
        s2t + 92*Mt*Sbeta)*pow2(Mt)*pow5(Mst1) + pow2(Mst2)*pow2(Sbeta)*(-328*
        pow2(Mt)*pow2(s2t)*pow4(Mst2) + 8*Dmglst1*(-60*s2t*pow2(Mst2)*pow3(Mt)
        + 41*Mt*pow3(s2t)*pow4(Mst2)) - 408*pow2(Mst2)*pow4(Mt) + pow2(Dmglst1)
        *(384*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 512*pow4(Mt) + 91*pow4(Mst2)*
        pow4(s2t)) + 41*pow4(s2t)*pow6(Mst2))) + 1920*s2t*pow2(Sbeta)*pow3(Mt)*
        pow7(Mst1) - 3840*Mt*pow2(Mst2)*pow2(Sbeta)*pow3(s2t)*pow7(Mst1) +
        1770*Dmglst1*Mst1*pow2(Sbeta)*pow4(s2t)*pow8(Mst2) - 75*pow2(Dmglst1)*
        pow2(Sbeta)*pow4(s2t)*pow8(Mst2) + 3585*pow2(Mst1)*pow2(Sbeta)*pow4(
        s2t)*pow8(Mst2))/(540.*pow2(Mst1)*pow4(Mst2)))/pow4(Mt)/pow2(Sbeta)*12.; 

   return result;
}

/**
 * 	@return returns the susy log^3 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double himalaya::H5::calc_coef_at_as2_no_sm_logs_log3(){

   const double result =
      ((-224*pow2(Sbeta)*pow4(Mt))/9.)/pow4(Mt)/pow2(Sbeta)*12.; 

   return result;
}

