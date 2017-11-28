#define Pi M_PI

#include "H5g1.hpp"
#include "HierarchyCalculator.hpp"
#include "Utils.hpp"
#include <type_traits>
#include <cmath>

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
himalaya::H5g1::H5g1(std::map<unsigned int, unsigned int> flagMap, double Al4p, double beta, double Dmglst1,
		 double lmMt, double lmMst1, double lmMst2, double lmMsq,
		 double Mgl, double Mt, double Mst1, double Mst2, double Msq, double MuSUSY,
		 double s2t,
		 int mdrFlag, int oneLoopFlag, int twoLoopFlag, int threeLoopFlag){
   // abbrev for tan(beta) and sin(beta)
   Tbeta = tan(beta);
   Sbeta = sin(beta);
   Cbeta = cos(beta);
   this -> Dmglst1 = Dmglst1;
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
   // zeta functions
   z2 = pow2(Pi) / 6.;
   z3 = 1.202056903159594;
   z4 = pow4(Pi) / 90.;
   // poly logs
   double pl412 = 0.51747906167389934317668576113647; // PolyLog[4,1/2]
   std::complex<double> pl2expPi3(0.27415567780803773941206919444, 1.014941606409653625021202554275); // PolyLog[2, Exp[I Pi / 3]]
   std::complex<double> pl3expPi6sqrt3(0.51928806536375962552715984277228, - 0.33358157526196370641686908633664); // PolyLog[3, Exp[- I Pi / 6] / Sqrt[3]]

   // polylog functions, checked
   B4 = (-4 * z2 * pow2(log(2)) + 2 / 3.* pow4(log(2)) - 13 / 2. * z4 + 16. * pl412);
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
   // expansion flags
   xDmglst1 = flagMap.at(HierarchyCalculator::xxDmglst1);
   xMsq = flagMap.at(HierarchyCalculator::xxMsq);
   
   s1 = 
   #include "../hierarchies/h5g1/sigS1Full.inc"
   ;
   s2 = 
   #include "../hierarchies/h5g1/sigS2Full.inc"
   ;
   s12 = 
   #include "../hierarchies/h5g1/sigS12Full.inc"
   ;
}

/**
 * 	@return The diagonal (1, 1) matrix element of the Higgs mass matrix as a double for the hierarchy 'H5g1'
 */
double himalaya::H5g1::getS1(){
   return s1;
}

/**
 * 	@return The diagonal (2, 2) matrix element of the Higgs mass matrix as a double for the hierarchy 'H5g1'
 */
double himalaya::H5g1::getS2(){
   return s2;
}

/**
 * 	@return The off-diagonal (1, 2) = (2, 1) matrix element of the Higgs mass matrix as a double for the hierarchy 'H5g1'
 */
double himalaya::H5g1::getS12(){
   return s12;
}

/**
 * 	@return returns the constant term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double himalaya::H5g1::calc_at_as2_no_logs(){

   const double result =
      (pow2(Mt)*pow2(MuSUSY)*((Mt*s2t*pow3(Mst1)*((3*(535 + (1255*Dmglst1*(
        Dmglst1 + Mgl))/pow2(Mgl)))/pow4(Msq) + (214598 - 38448*B4 - 216*DN + (
        3*Dmglst1*(27938 - 23760*B4 + 216*DN)*(Dmglst1 + Mgl))/pow2(Mgl))/pow4(
        Mst2) + (90*(-2*(409*Dmglst1*Mgl + 1392*pow2(Dmglst1) + 67*pow2(Mgl))*
        pow2(Mst1)*pow2(Mst2) - 8*pow2(Msq)*(4*pow2(Mgl)*(14*pow2(Mst1) + 17*
        pow2(Mst2)) + 6*Dmglst1*Mgl*(29*pow2(Mst1) + 41*pow2(Mst2)) + pow2(
        Dmglst1)*(203*pow2(Mst1) + 246*pow2(Mst2)))))/(pow2(Mgl)*pow4(Msq)*
        pow4(Mst2))))/243. + pow2(s2t)*((1672*pow2(Mst1))/(135.*pow2(Msq)) + (-
        (pow2(Mst1)*(225000*Dmglst1*(67*Dmglst1 + 22*Mgl)*pow2(Mst1) + pow2(
        Mgl)*((3228977 - 1026000*B4 + 27000*DN)*pow2(Msq) + 1372500*pow2(Mst1))
        ))/(60750.*pow2(Msq)*pow2(Mst2)) + (257250*Dmglst1*(3087*Dmglst1 + 970*
        Mgl)*pow4(Mst1) + pow2(Mgl)*(-19447774*pow2(Mst1)*pow2(Mst2) +
        41220947*pow4(Mst1) - 7399303*pow4(Mst2)))/(5.5566e6*pow4(Msq)))/pow2(
        Mgl) + ((-4*S2*pow2(Mst1)*(2069*pow2(Mst1) + 882*pow2(Mst2)))/27. - (
        531.3038044317668 - (152*B4)/9. + (4*DN)/9. + (Dmglst1*(
        224.98199588477365 - (16*B4)/3. - (40*DN)/9.))/Mgl + ((-130639 + 24192*
        B4 + 1728*DN)*pow2(Dmglst1))/(3888.*pow2(Mgl)))*pow4(Mst1) + (32*OepS2*
        (102*pow2(Mst1)*pow2(Mst2) + 185*pow4(Mst1)))/729.)/pow4(Mst2)) - (16*(
        -139 + 2*B4 + DN)*pow2(Mst1)*pow2(Mt))/(9.*pow4(Mst2))) + z4*((4*Mt*
        MuSUSY*pow2(Mst1)*(-54*Dmglst1*Mgl*Mst1*Mt*(Mt*(-507*MuSUSY*s2t + 60*
        Cbeta*Mst1*s2t*Sbeta) - 4*Cbeta*Sbeta*pow2(Mt) + (265*Mst1*MuSUSY +
        561*Cbeta*Sbeta*pow2(Mst1) + 561*Cbeta*Sbeta*pow2(Mst2))*pow2(s2t)) -
        27*Mst1*pow2(Dmglst1)*(-6*s2t*(169*MuSUSY - 30*Cbeta*Mst1*Sbeta)*pow2(
        Mt) + 3*Mt*(185*Mst1*MuSUSY + 989*Cbeta*Sbeta*pow2(Mst1) + 374*Cbeta*
        Sbeta*pow2(Mst2))*pow2(s2t) - 8*Cbeta*Sbeta*pow3(Mt) + 120*Cbeta*Mst1*
        Sbeta*pow2(Mst2)*pow3(s2t)) + pow2(Mgl)*(-2*s2t*(8451*Mst1*MuSUSY + 20*
        Cbeta*Sbeta*pow2(Mst1) + 21*Cbeta*Sbeta*pow2(Mst2))*pow2(Mt) + Mt*pow2(
        s2t)*(-127*MuSUSY*pow2(Mst1) + 39*MuSUSY*pow2(Mst2) + 2916*Cbeta*Mst1*
        Sbeta*pow2(Mst2) + 2916*Cbeta*Sbeta*pow3(Mst1)) + 108*(15*MuSUSY + 2*
        Cbeta*Mst1*Sbeta)*pow3(Mt) + Cbeta*Sbeta*pow2(Mst2)*(83*pow2(Mst1) +
        1614*pow2(Mst2))*pow3(s2t))))/(243.*pow2(Mgl)*pow4(Mst2)) - pow2(Sbeta)
        *(pow2(Mt)*pow2(s2t)*(2*pow2(Mst2) - pow2(Mst1)*((80*Dmglst1)/(3.*Mgl)
        + (40*pow2(Dmglst1))/pow2(Mgl) + (2*(95 - (26*pow2(MuSUSY))/pow2(Mst2))
        )/81.) - (4*(-(pow2(Mgl)*(pow2(Mst2) - 127*pow2(MuSUSY))) + 135*
        Dmglst1*(111*Dmglst1 + 106*Mgl)*pow2(MuSUSY))*pow4(Mst1))/(243.*pow2(
        Mgl)*pow4(Mst2))) + (80*pow2(Mst1)*pow2(MuSUSY)*pow4(Mt))/(3.*pow4(
        Mst2)) + (Mst1*((-2*Mt*(pow2(Mgl)*(169*pow2(Mst1) - 241*pow2(Mst2)) +
        Dmglst1*Mgl*(989*pow2(Mst1) - 241*pow2(Mst2)) + pow2(Dmglst1)*(2219*
        pow2(Mst1) - 241*pow2(Mst2)))*pow3(s2t))/9. + (8*s2t*pow3(Mt)*(4*(
        Dmglst1*Mgl + pow2(Dmglst1) + pow2(Mgl)) + ((507*Dmglst1*(Dmglst1 +
        Mgl) - 313*pow2(Mgl))*pow2(Mst1)*pow2(MuSUSY))/pow4(Mst2)))/9.) - ((
        14310*Dmglst1*Mgl*pow2(Mst1)*(pow2(Mst1) - pow2(Mst2)) + 405*pow2(
        Dmglst1)*(-53*pow2(Mst1)*pow2(Mst2) + 69*pow4(Mst1)) + pow2(Mgl)*(-
        6495*pow2(Mst1)*pow2(Mst2) + 3062*pow4(Mst1) + 3267*pow4(Mst2)))*pow4(
        s2t))/486.)/pow2(Mgl))) - (pow2(Mgl)*(-6912*Mst1*Mt*s2t*(pow2(Mt) +
        pow2(Mst2)*pow2(s2t))*pow2(Sbeta)*pow2(z2)*pow4(Mst2) + 3456*Mt*pow2(
        z2)*pow3(Mst1)*(27*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*pow2(s2t) + 44*s2t*
        pow2(Mt)*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 2*Cbeta*MuSUSY*Sbeta*pow3(
        Mt) - 7*pow2(Sbeta)*pow3(s2t)*pow4(Mst2)) + 6*s2t*pow2(Mst1)*pow2(Mst2)
        *(-2*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*pow2(s2t)*(236*T1ep + 537*pow2(
        z2)) - 2*s2t*pow2(Mt)*(pow2(Mst2)*pow2(Sbeta)*(20*T1ep + 87*pow2(z2)) +
        8*pow2(MuSUSY)*(-1 + pow2(Sbeta))*(68*T1ep + 213*pow2(z2))) + 56*Cbeta*
        MuSUSY*Sbeta*(4*T1ep + 3*pow2(z2))*pow3(Mt) + pow2(Sbeta)*(100*T1ep +
        111*pow2(z2))*pow3(s2t)*pow4(Mst2)) + s2t*pow4(Mst1)*(-664*Cbeta*Mt*
        MuSUSY*Sbeta*pow2(Mst2)*pow2(s2t)*(4*T1ep + 3*pow2(z2)) - 8*s2t*pow2(
        Mt)*(-(pow2(Mst2)*pow2(Sbeta)*(4*T1ep + 3*pow2(z2))) + 2*pow2(MuSUSY)*(
        -1 + pow2(Sbeta))*(740*T1ep + 1527*pow2(z2))) + 320*Cbeta*MuSUSY*Sbeta*
        (4*T1ep + 3*pow2(z2))*pow3(Mt) - pow2(Sbeta)*(44*T1ep + 1113*pow2(z2))*
        pow3(s2t)*pow4(Mst2)) + 93312*Cbeta*MuSUSY*Sbeta*pow2(Mt)*pow2(s2t)*
        pow2(z2)*pow5(Mst1) + 27*pow2(s2t)*pow2(Sbeta)*(-4*pow2(Mt)*(4*T1ep -
        5*pow2(z2)) + pow2(Mst2)*pow2(s2t)*(4*T1ep + 35*pow2(z2)))*pow6(Mst2))
        + 432*Dmglst1*Mgl*Mst1*pow2(z2)*(456*Cbeta*MuSUSY*Sbeta*pow2(Mt)*pow2(
        s2t)*pow4(Mst1) - 16*Mt*s2t*(pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow2(
        Sbeta)*pow4(Mst2) - pow2(s2t)*pow3(Mst1)*(8*pow2(Mt)*pow2(MuSUSY)*(-1 +
        pow2(Sbeta)) + pow2(s2t)*pow2(Sbeta)*pow4(Mst2)) + 8*Mt*pow2(Mst1)*(57*
        Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*pow2(s2t) + 84*s2t*pow2(Mt)*pow2(
        MuSUSY)*(-1 + pow2(Sbeta)) + 2*Cbeta*MuSUSY*Sbeta*pow3(Mt) - 17*pow2(
        Sbeta)*pow3(s2t)*pow4(Mst2)) + Mst1*pow2(Sbeta)*pow4(s2t)*pow6(Mst2)) +
        216*Mst1*pow2(Dmglst1)*pow2(z2)*(1632*Cbeta*MuSUSY*Sbeta*pow2(Mt)*pow2(
        s2t)*pow4(Mst1) - 32*Mt*s2t*(pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow2(
        Sbeta)*pow4(Mst2) - 3*pow2(s2t)*pow3(Mst1)*(8*pow2(Mt)*pow2(MuSUSY)*(-1
        + pow2(Sbeta)) + pow2(s2t)*pow2(Sbeta)*pow4(Mst2)) + 16*Mt*pow2(Mst1)*(
        57*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*pow2(s2t) + 84*s2t*pow2(Mt)*pow2(
        MuSUSY)*(-1 + pow2(Sbeta)) + 2*Cbeta*MuSUSY*Sbeta*pow3(Mt) - 32*pow2(
        Sbeta)*pow3(s2t)*pow4(Mst2)) + 3*Mst1*pow2(Sbeta)*pow4(s2t)*pow6(Mst2))
        )/(972.*pow2(Mgl)*pow4(Mst2)) + (z3*((32*pow2(Mst1)*pow2(Mt)*pow2(
        MuSUSY)*(216*Dmglst1*Mgl*Mst1*s2t*(-231*Mt + 382*Mst1*s2t) + 27*Mst1*
        s2t*(-1848*Mt + 6491*Mst1*s2t)*pow2(Dmglst1) + pow2(Mgl)*(-62100*Mst1*
        Mt*s2t + 4968*pow2(Mt) + (67606*pow2(Mst1) + 26349*pow2(Mst2))*pow2(
        s2t))))/pow4(Mst2) - (16*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst1)*(27*Dmglst1*
        Mgl*Mst1*(320*Mst1*s2t*pow2(Mt) + 3*Mt*(2183*pow2(Mst1) - 808*pow2(
        Mst2))*pow2(s2t) - 3760*pow3(Mt) + 2846*Mst1*pow2(Mst2)*pow3(s2t)) +
        27*Mst1*pow2(Dmglst1)*(2464*Mst1*s2t*pow2(Mt) + 3*Mt*(9453*pow2(Mst1) -
        808*pow2(Mst2))*pow2(s2t) - 3760*pow3(Mt) + 6197*Mst1*pow2(Mst2)*pow3(
        s2t)) - pow2(Mgl)*(4260*s2t*pow2(Mst2)*pow2(Mt) + 29565*Mt*pow2(s2t)*
        pow3(Mst1) + 324*Mst1*(221*Mt*pow2(Mst2)*pow2(s2t) + 86*pow3(Mt)) +
        pow2(Mst1)*(2936*s2t*pow2(Mt) - 41257*pow2(Mst2)*pow3(s2t)) - 21921*
        pow3(s2t)*pow4(Mst2))))/pow4(Mst2) + pow2(Sbeta)*((27*pow2(Dmglst1)*(
        4242*pow2(Mst1)*pow2(Mst2) + 22079*pow4(Mst1) - 1533*pow4(Mst2)) + 54*
        Dmglst1*Mgl*(810*pow2(Mst1)*pow2(Mst2) + 5077*pow4(Mst1) - 195*pow4(
        Mst2)) + 4*pow2(Mgl)*(17493*pow2(Mst1)*pow2(Mst2) + 19336*pow4(Mst1) +
        4428*pow4(Mst2)))*pow4(s2t) + (4*Mt*(108*Mst1*pow2(Mst2)*pow3(s2t)*(
        pow2(Mgl)*(-618*pow2(Mst1)*pow2(Mst2) + 519*pow4(Mst1) - 266*pow4(Mst2)
        ) + Dmglst1*Mgl*(-692*pow2(Mst1)*pow2(Mst2) + 2991*pow4(Mst1) - 116*
        pow4(Mst2)) + pow2(Dmglst1)*(-665*pow2(Mst1)*pow2(Mst2) + 9559*pow4(
        Mst1) + 559*pow4(Mst2))) + Mt*(54*pow2(Mt)*(pow2(Dmglst1)*(18608*pow2(
        Mst1)*pow2(Mst2) + 26736*pow4(Mst1) - 2019*pow4(Mst2)) + 8*pow2(Mgl)*(
        4*pow2(Mst1)*(36*pow2(Mst2) - 23*pow2(MuSUSY)) + 317*pow4(Mst1) - 163*
        pow4(Mst2)) + 2*Dmglst1*Mgl*(3088*pow2(Mst1)*pow2(Mst2) + 5264*pow4(
        Mst1) - 107*pow4(Mst2))) - 216*Mst1*Mt*s2t*(4*Mgl*(Mgl*(pow2(Mst1)*(
        183*pow2(Mst2) - 575*pow2(MuSUSY)) + 448*pow4(Mst1) - 54*pow4(Mst2)) +
        Dmglst1*(pow2(Mst1)*(507*pow2(Mst2) - 462*pow2(MuSUSY)) + 2196*pow4(
        Mst1) - 37*pow4(Mst2))) + pow2(Dmglst1)*(84*pow2(Mst1)*(45*pow2(Mst2) -
        22*pow2(MuSUSY)) + 25672*pow4(Mst1) + 1695*pow4(Mst2))) - pow2(s2t)*(
        27*pow2(Dmglst1)*(-8*(2189*pow2(Mst2) - 6491*pow2(MuSUSY))*pow4(Mst1) +
        16063*pow2(Mst1)*pow4(Mst2) - 3479*pow6(Mst2)) + 216*Dmglst1*Mgl*((-
        398*pow2(Mst2) + 3056*pow2(MuSUSY))*pow4(Mst1) + 373*pow2(Mst1)*pow4(
        Mst2) - 55*pow6(Mst2)) + 4*pow2(Mgl)*((-662*pow2(Mst2) + 135212*pow2(
        MuSUSY))*pow4(Mst1) + pow2(Mst1)*(52698*pow2(Mst2)*pow2(MuSUSY) + 969*
        pow4(Mst2)) + 1161*pow6(Mst2))))))/pow4(Mst2))))/(1944.*pow2(Mgl)) + (
        pow2(Sbeta)*(296352000*Mt*z2*pow2(Mst2)*pow3(Mst1)*pow3(s2t)*(3*pow2(
        Dmglst1)*(-450*pow4(Mst1)*pow4(Mst2) + pow4(Msq)*(-8630*pow2(Mst1)*
        pow2(Mst2) + 445*pow4(Mst1) + 1254*pow4(Mst2)) + 240*pow2(Msq)*(13*
        pow2(Mst2)*pow4(Mst1) - 2*pow2(Mst1)*pow4(Mst2))) + pow2(Mgl)*(-90*
        pow4(Mst1)*pow4(Mst2) + pow4(Msq)*(-2710*pow2(Mst1)*pow2(Mst2) + 737*
        pow4(Mst1) + 3762*pow4(Mst2)) + 240*pow2(Msq)*(2*pow2(Mst2)*pow4(Mst1)
        - pow2(Mst1)*pow4(Mst2))) + Dmglst1*Mgl*(-450*pow4(Mst1)*pow4(Mst2) +
        pow4(Msq)*(-11646*pow2(Mst1)*pow2(Mst2) + 1397*pow4(Mst1) + 3762*pow4(
        Mst2)) + 720*pow2(Msq)*(4*pow2(Mst2)*pow4(Mst1) - pow2(Mst1)*pow4(Mst2)
        ))) + pow2(Mst2)*pow4(s2t)*(257250*pow2(Dmglst1)*pow2(Mst1)*pow2(Mst2)*
        (-2880*pow2(Msq)*pow2(Mst1)*pow2(Mst2)*((659 + 300*z2)*pow2(Mst1) - 27*
        (6 + z2)*pow2(Mst2)) + 2700*(213 + 40*z2)*pow4(Mst1)*pow4(Mst2) + pow4(
        Msq)*(144*(-3583 + 216*B4 + 180*DN + 10476*z2)*pow2(Mst1)*pow2(Mst2) -
        (232169 + 86400*B4 + 53568*DN + 4204944*z2)*pow4(Mst1) + 792*(587 +
        324*z2)*pow4(Mst2))) + 514500*Dmglst1*Mgl*pow2(Mst1)*pow2(Mst2)*(-2880*
        pow2(Msq)*pow2(Mst1)*pow2(Mst2)*((137 + 60*z2)*pow2(Mst1) - (41 + 9*z2)
        *pow2(Mst2)) + 180*(521 + 120*z2)*pow4(Mst1)*pow4(Mst2) + pow4(Msq)*(
        144*(1151 + 72*B4 + 60*DN + 3492*z2)*pow2(Mst1)*pow2(Mst2) - (773389 +
        10368*B4 + 8640*DN + 1030032*z2)*pow4(Mst1) + 216*(7 + 396*z2)*pow4(
        Mst2))) + pow2(Mgl)*(-29635200*pow2(Msq)*((2068 + 750*z2)*pow2(Mst1) -
        (1357 + 225*z2)*pow2(Mst2))*pow4(Mst1)*pow4(Mst2) + pow2(Mst1)*pow2(
        Mst2)*pow4(Msq)*(65856*(-5876588 + 13500*B4 + 20250*DN + 50000*OepS2 -
        3989250*S2 + 2151375*z2)*pow2(Mst1)*pow2(Mst2) - (257823187891 +
        8890560000*B4 + 444528000*DN + 241472000*OepS2 - 20818728000*S2 +
        202706826000*z2)*pow4(Mst1) + 3087000*(53749 + 2592*B4 - 288*DN + 192*
        OepS2 + 21384*S2 - 23544*z2)*pow4(Mst2)) - 6667920000*(-1 + 2*z2)*(
        pow2(Mst1) + pow2(Mst2))*pow2(pow2(Mst1) - pow2(Mst2))*pow6(Msq) +
        12960*(1009961 + 214375*z2)*pow6(Mst1)*pow6(Mst2))) - 192*pow4(Mt)*(
        171500*Dmglst1*Mgl*pow2(Mst1)*(360*pow2(Msq)*pow2(Mst1)*pow2(Mst2)*(
        110*pow2(Mst1) + 9*(9 - 8*z2)*pow2(Mst2)) - 15*pow2(Mst1)*((-2095 +
        1872*z2)*pow2(Mst1) + 6*(-115 + 72*z2)*pow2(Mst2))*pow4(Mst2) + 2*pow4(
        Msq)*(2*(64093 - 10872*z2)*pow2(Mst1)*pow2(Mst2) + 3*(99985 + 67176*z2)
        *pow4(Mst1) + 3*(2215 - 3564*z2)*pow4(Mst2))) + 3430*pow2(Dmglst1)*
        pow2(Mst1)*(720*pow2(Msq)*pow2(Mst1)*pow2(Mst2)*(6875*pow2(Mst1) + (
        3349 - 2700*z2)*pow2(Mst2)) + 15*pow2(Mst1)*((215363 - 212400*z2)*pow2(
        Mst1) + 450*(115 - 72*z2)*pow2(Mst2))*pow4(Mst2) + 2*pow4(Msq)*(12*(
        1477147 - 135900*z2)*pow2(Mst1)*pow2(Mst2) + (37614091 + 27183600*z2)*
        pow4(Mst1) - 12*(60647 + 66825*z2)*pow4(Mst2))) + pow2(Mgl)*(82320*
        pow2(Msq)*pow2(Mst1)*pow2(Mst2)*(-9*(-4113 + 4000*z2)*pow2(Mst1)*pow2(
        Mst2) + 20625*pow4(Mst1) + (13517 - 9000*z2)*pow4(Mst2)) - 15*pow2(
        Mst1)*pow4(Mst2)*(514500*(-115 + 72*z2)*pow2(Mst1)*pow2(Mst2) + (-
        160842737 + 117306000*z2)*pow4(Mst1) + 3*(-13597579 + 6174000*z2)*pow4(
        Mst2)) + 555660000*(-1 + 2*z2)*pow2(Mst2)*(pow2(Mst1) + pow2(Mst2))*
        pow6(Msq) + 4*pow4(Msq)*(-2744*((-1229794 + 325125*z2)*pow2(Mst2) +
        6750*(-139 + 2*B4 + DN - 110*z2)*pow2(MuSUSY))*pow4(Mst1) - 128625*(
        10501 + 4704*z2)*pow2(Mst1)*pow4(Mst2) + (6336529126 + 3678160500*z2)*
        pow6(Mst1) + 13891500*(-1 + 2*z2)*pow6(Mst2)))) + 32928000*s2t*pow3(
        Mst1)*pow3(Mt)*(pow2(Mgl)*(2*pow4(Msq)*(pow2(Mst1)*(6*(83 + 16668*z2)*
        pow2(Mst2) + (-107299 + 19224*B4 + 108*DN - 86652*z2)*pow2(MuSUSY)) +
        26*(-1561 + 6354*z2)*pow4(Mst1) + 18*(371 + 52*z2)*pow4(Mst2)) - 120*
        pow2(Msq)*(4*pow2(Mst1)*pow2(Mst2)*((-55 + 18*z2)*pow2(Mst2) - 6*(17 +
        3*z2)*pow2(MuSUSY)) + 6*((3 + 36*z2)*pow2(Mst2) - 56*pow2(MuSUSY))*
        pow4(Mst1) + (67 - 72*z2)*pow6(Mst2)) + 15*pow2(Mst2)*(12*(-4*(-25 + 9*
        z2)*pow2(Mst2) + (67 + 18*z2)*pow2(MuSUSY))*pow4(Mst1) + pow2(Mst1)*(-
        107*pow2(Mst2)*pow2(MuSUSY) + 32*(-14 + 9*z2)*pow4(Mst2)) + 2*(-169 +
        72*z2)*pow6(Mst2))) + Dmglst1*Mgl*(2*pow4(Msq)*(pow2(Mst1)*((49442 +
        298800*z2)*pow2(Mst2) + 3*(-13969 + 11880*B4 - 108*DN + 24732*z2)*pow2(
        MuSUSY)) + 2*(54037 + 398142*z2)*pow4(Mst1) + 78*(539 + 12*z2)*pow4(
        Mst2)) - 120*pow2(Msq)*(12*pow2(Mst1)*pow2(Mst2)*((-64 + 6*z2)*pow2(
        Mst2) - 3*(41 + 6*z2)*pow2(MuSUSY)) + 18*((11 + 60*z2)*pow2(Mst2) + 2*(
        -29 + 12*z2)*pow2(MuSUSY))*pow4(Mst1) + (67 - 72*z2)*pow6(Mst2)) + 15*
        pow2(Mst2)*(-12*(4*(-92 + 21*z2)*pow2(Mst2) - (409 + 90*z2)*pow2(
        MuSUSY))*pow4(Mst1) + pow2(Mst1)*(-251*pow2(Mst2)*pow2(MuSUSY) + 96*(-
        14 + 9*z2)*pow4(Mst2)) + 2*(-169 + 72*z2)*pow6(Mst2))) + pow2(Dmglst1)*
        (2*pow4(Msq)*(pow2(Mst1)*(2*(57151 + 298494*z2)*pow2(Mst2) + 3*(-13969
        + 11880*B4 - 108*DN + 24732*z2)*pow2(MuSUSY)) + 72*(11194 + 32437*z2)*
        pow4(Mst1) + 24*(9208 + 39*z2)*pow4(Mst2)) - 120*pow2(Msq)*(6*pow2(
        Mst1)*pow2(Mst2)*((-331 + 12*z2)*pow2(Mst2) - 6*(41 + 6*z2)*pow2(
        MuSUSY)) + 6*(36*(4 + 15*z2)*pow2(Mst2) + (-203 + 324*z2)*pow2(MuSUSY))
        *pow4(Mst1) + (67 - 72*z2)*pow6(Mst2)) + 15*pow2(Mst2)*(-72*((-185 +
        26*z2)*pow2(Mst2) - (232 + 45*z2)*pow2(MuSUSY))*pow4(Mst1) + pow2(Mst1)
        *(-251*pow2(Mst2)*pow2(MuSUSY) + 192*(-14 + 9*z2)*pow4(Mst2)) + 2*(-169
        + 72*z2)*pow6(Mst2)))) + 8*pow2(Mt)*pow2(s2t)*(-257250*pow2(Dmglst1)*
        pow2(Mst1)*(1440*pow2(Msq)*pow2(Mst1)*pow2(Mst2)*(pow2(Mst1)*(577*pow2(
        Mst2) - 2*(335 + 246*z2)*pow2(MuSUSY)) - 27*pow4(Mst2)) + 180*pow2(
        Mst1)*pow4(Mst2)*(pow2(Mst1)*(497*pow2(Mst2) + 3*(1029 + 200*z2)*pow2(
        MuSUSY)) + 117*pow4(Mst2)) + pow4(Msq)*((8*(135281 + 534888*z2)*pow2(
        Mst2) - (-130639 + 24192*B4 + 1728*DN + 418032*z2)*pow2(MuSUSY))*pow4(
        Mst1) + 24*(-266599 + 864*B4 + 432*DN - 110664*z2)*pow2(Mst1)*pow4(
        Mst2) + 3888*(163 + 66*z2)*pow6(Mst2))) - 514500*Dmglst1*Mgl*pow2(Mst1)
        *(1440*pow2(Msq)*pow2(Mst1)*pow2(Mst2)*(pow2(Mst1)*(131*pow2(Mst2) - 2*
        (55 + 42*z2)*pow2(MuSUSY)) - 21*pow4(Mst2)) - 180*pow2(Mst1)*(pow2(
        Mst1)*(17*pow2(Mst2) - 5*(97 + 24*z2)*pow2(MuSUSY)) - 39*pow4(Mst2))*
        pow4(Mst2) + pow4(Msq)*((8*(63295 + 66456*z2)*pow2(Mst2) + (-437365 +
        10368*B4 + 8640*DN + 232272*z2)*pow2(MuSUSY))*pow4(Mst1) + 8*(-156781 +
        864*B4 + 432*DN - 110664*z2)*pow2(Mst1)*pow4(Mst2) + 288*(-79 + 297*z2)
        *pow6(Mst2))) + pow2(Mgl)*(180*pow2(Mst1)*pow4(Mst2)*((15671760*pow2(
        Mst2) - (41220947 + 15435000*z2)*pow2(MuSUSY))*pow4(Mst1) + 7399303*
        pow2(MuSUSY)*pow4(Mst2) - 2*pow2(Mst1)*(-9723887*pow2(Mst2)*pow2(
        MuSUSY) + 6570120*pow4(Mst2))) + 6667920000*(-1 + 2*z2)*pow2(-(Mst2*
        pow2(Mst1)) + pow3(Mst2))*pow6(Msq) - 987840*pow2(Msq)*pow2(Mst1)*pow2(
        Mst2)*(pow2(Mst1)*pow2(Mst2)*(3091*pow2(Mst2) + 30*(418 + 225*z2)*pow2(
        MuSUSY)) + (22642*pow2(Mst2) - 375*(61 + 24*z2)*pow2(MuSUSY))*pow4(
        Mst1) - 5108*pow6(Mst2)) + pow4(Msq)*(-8232*pow2(Mst2)*((-32617079 +
        108000*DN + 20000*OepS2 - 2470500*S2 - 30213750*z2)*pow2(Mst2) + 2*(-
        3228977 + 1026000*B4 - 27000*DN + 272000*OepS2 - 7938000*S2 - 466500*
        z2)*pow2(MuSUSY))*pow4(Mst1) + (16*(-14893594207 + 1372000*OepS2 -
        171328500*S2 + 75888750*z2)*pow2(Mst2) + (531403689547 - 16892064000*B4
        + 444528000*DN - 8122240000*OepS2 + 306576144000*S2 + 142720242000*z2)*
        pow2(MuSUSY))*pow6(Mst1) + 3087000*(-4501 + 288*DN - 96*OepS2 + 972*S2
        - 10998*z2)*pow2(Mst1)*pow6(Mst2) + 666792000*(-1 + 2*z2)*pow8(Mst2))))
        ))/(8.001504e9*pow2(Mgl)*pow2(Mst1)*pow4(Msq)*pow4(Mst2)) + Mt*(-(Mst1*
        pow2(Sbeta)*pow3(s2t)*(8*(132407 + 11880*B4 - 108*DN + (3*Dmglst1*(
        Dmglst1*(57101 + 17640*B4 - 396*DN) + (54053 + 9432*B4 - 180*DN)*Mgl))/
        pow2(Mgl))*pow2(Mst1) + ((pow2(Msq)*(pow2(Dmglst1)*(9306949*pow2(Msq) -
        2632320*pow2(Mst2)) + Dmglst1*Mgl*(3058187*pow2(Msq) - 924480*pow2(
        Mst2)) + pow2(Mgl)*(622151*pow2(Msq) - 239040*pow2(Mst2)))*pow4(Mst1))/
        pow2(Mst2) - 36*pow2(Mst2)*(-(pow2(Mgl)*(5920*pow2(Msq)*pow2(Mst1) + 6*
        (-4439 + 136*B4 + 4*DN)*pow4(Msq) + 1715*pow4(Mst1))) - Dmglst1*Mgl*(
        20160*pow2(Msq)*pow2(Mst1) + 6*(-3779 + 136*B4 + 4*DN)*pow4(Msq) +
        9035*pow4(Mst1)) - pow2(Dmglst1)*(45040*pow2(Msq)*pow2(Mst1) + (826 +
        816*B4 + 24*DN)*pow4(Msq) + 29415*pow4(Mst1))) - 180*(pow2(Mgl)*(48*
        pow2(Msq) + 43*pow2(Mst1)) + Dmglst1*Mgl*(48*pow2(Msq) + 91*pow2(Mst1))
        + pow2(Dmglst1)*(48*pow2(Msq) + 163*pow2(Mst1)))*pow4(Mst2))/(pow2(Mgl)
        *pow4(Msq))))/3888. + MuSUSY*(-(z2*(36*Mt*MuSUSY*(-27*pow2(s2t)*(pow2(
        Mst1)*(20/(3.*pow2(Msq)) - (2*(311 + (180*(14*Dmglst1*Mgl + 41*pow2(
        Dmglst1) + 2*pow2(Mgl))*pow2(Mst1))/(pow2(Mgl)*pow2(Msq))))/(81.*pow2(
        Mst2))) + (25*(4*Dmglst1*Mgl + 10*pow2(Dmglst1) + pow2(Mgl))*pow4(Mst1)
        )/(9.*pow2(Mgl)*pow4(Msq)) - ((-58068*Dmglst1*Mgl + 52254*pow2(Dmglst1)
        + 69349*pow2(Mgl))*pow4(Mst1))/(486.*pow2(Mgl)*pow4(Mst2))) + (8*Mt*
        pow2(Mst1)*(-660*Mt - Mst1*s2t*(2407 - (2061*Dmglst1*(Dmglst1 + Mgl))/
        pow2(Mgl)) + (15*Mst1*s2t*(pow2(Mgl)*(8*pow2(Msq) + 3*pow2(Mst1))*pow2(
        Mst2) + 3*pow2(Dmglst1)*(15*pow2(Mst1)*pow2(Mst2) + 8*pow2(Msq)*(-9*
        pow2(Mst1) + pow2(Mst2))) + 3*Dmglst1*Mgl*(5*pow2(Mst1)*pow2(Mst2) + 8*
        pow2(Msq)*(-2*pow2(Mst1) + pow2(Mst2)))))/(pow2(Mgl)*pow4(Msq))))/pow4(
        Mst2)) + (Cbeta*Sbeta*((8*pow2(Mt)*(s2t*(33516*Dmglst1*Mgl*pow2(Mst1) -
        117486*pow2(Dmglst1)*pow2(Mst1) + pow2(Mgl)*(52798*pow2(Mst1) + 52503*
        pow2(Mst2))) - (12960*Mt*(5*Dmglst1*Mgl + 15*pow2(Dmglst1) + pow2(Mgl))
        *pow3(Mst1))/pow2(Msq))*pow4(Mst1))/(pow2(Mgl)*pow4(Mst2)) + (-108*Mt*
        pow2(s2t)*((-4*(263 + (-1971*Dmglst1*(Dmglst1 + Mgl) + (60*(9*Dmglst1*
        Mgl + 33*pow2(Dmglst1) + pow2(Mgl))*pow2(Mst1))/pow2(Msq))/pow2(Mgl)))/
        pow2(Mst2) + ((240*(3*Dmglst1*Mgl + 3*pow2(Dmglst1) + pow2(Mgl)))/pow2(
        Msq) + pow2(Mst1)*((90*(5*Dmglst1*Mgl + 15*pow2(Dmglst1) + pow2(Mgl)))/
        pow4(Msq) + (6487*Dmglst1*Mgl + 20793*pow2(Dmglst1) - 1789*pow2(Mgl))/
        pow4(Mst2)))/pow2(Mgl)) + 576*pow3(Mt)*((15*(7*Dmglst1*Mgl + 7*pow2(
        Dmglst1) + 3*pow2(Mgl)))/(pow2(Mgl)*pow4(Msq)) + (1381 + (4163*Dmglst1*
        (Dmglst1 + Mgl))/pow2(Mgl))/pow4(Mst2)))*pow5(Mst1) - pow3(s2t)*(120*(-
        274 - (54*pow2(Msq))/pow2(Mst2) - (27*pow2(Mst2))/pow2(Msq))*pow4(Mst1)
        + 648*pow4(Mst2) + ((pow2(Mgl)*(7560*pow2(Msq)*pow2(Mst2) + 65617*pow4(
        Msq) - 1350*pow4(Mst2)) + 18*pow2(Dmglst1)*(5460*pow2(Msq)*pow2(Mst2) +
        16943*pow4(Msq) - 750*pow4(Mst2)) + 36*Dmglst1*Mgl*(1020*pow2(Msq)*
        pow2(Mst2) + 3067*pow4(Msq) - 150*pow4(Mst2)))*pow6(Mst1))/(pow2(Mgl)*
        pow2(Mst2)*pow4(Msq)))))/pow2(Mst1)))/972. + Cbeta*Sbeta*((pow3(Mst1)*
        pow3(Mt)*(81*((5*(669 + (1409*Dmglst1*(Dmglst1 + Mgl))/pow2(Mgl)))/(81.
        *pow4(Msq)) + (4*(5241 - (16303*Dmglst1*(Dmglst1 + Mgl))/pow2(Mgl)))/(
        243.*pow4(Mst2))) + (20*(-((871*Dmglst1*Mgl + 3220*pow2(Dmglst1) + 135*
        pow2(Mgl))*pow2(Mst1)*pow2(Mst2)) - 4*pow2(Msq)*(439*Dmglst1*Mgl*(pow2(
        Mst1) + pow2(Mst2)) + 99*pow2(Mgl)*(pow2(Mst1) + pow2(Mst2)) + pow2(
        Dmglst1)*(1174*pow2(Mst1) + 439*pow2(Mst2)))))/(pow2(Mgl)*pow4(Msq)*
        pow4(Mst2))))/81. + pow3(s2t)*(pow2(Mst1)*(110.18859670781893 - (40*B4)
        /9. - (2*DN)/9. - ((85750*Dmglst1*(3141*Dmglst1 + 1006*Mgl)*pow2(Mst1)
        + pow2(Mgl)*(48706000*pow2(Msq) + 20222907*pow2(Mst1)))*pow2(Mst2))/(
        3.7044e6*pow2(Mgl)*pow4(Msq))) + ((9600*Dmglst1*Mgl + 24850*pow2(
        Dmglst1) + 2361*pow2(Mgl))*pow4(Mst1))/(135.*pow2(Mgl)*pow2(Msq)) + ((
        pow2(Mst1)*(90*pow2(Msq) + 2374*S2*pow2(Mst1) + 2655*S2*pow2(Mst2)))/
        27. + (239.07595982905215 + (606133*Dmglst1)/(3888.*Mgl) + ((
        36.42193930041152 + (64*B4)/9. + (32*DN)/9.)*pow2(Dmglst1))/pow2(Mgl))*
        pow4(Mst1) - (8*OepS2*(177*pow2(Mst1)*pow2(Mst2) + 166*pow4(Mst1)))/
        729.)/pow2(Mst2) + (313/(45.*pow2(Msq)) - 1/(3.*pow2(Mst1)) + ((
        1.0841585681891803 + (5*Dmglst1)/(6.*Mgl) + (5*pow2(Dmglst1))/(4.*pow2(
        Mgl)))*pow2(Mst1))/pow4(Msq))*pow4(Mst2)) + Mt*pow2(s2t)*(((20*(249*
        Dmglst1*(Dmglst1 + Mgl) + 71*pow2(Mgl)))/(9.*pow2(Mgl)*pow2(Msq)) + (
        77.49382716049382 + 96*B4 + ((Dmglst1*(10021 + 5328*B4 - 72*DN)*(
        Dmglst1 + Mgl))/27. - (40*(12*Dmglst1*Mgl + 59*pow2(Dmglst1) + 2*pow2(
        Mgl))*pow2(Mst1))/(3.*pow2(Msq)))/pow2(Mgl))/pow2(Mst2))*pow3(Mst1) - (
        5*Mst1*pow2(Mst2)*(2*(59 + (131*Dmglst1*(Dmglst1 + Mgl))/pow2(Mgl))*
        pow2(Mst1) - (11*(Dmglst1*Mgl + pow2(Dmglst1) + pow2(Mgl))*pow2(Mst2))/
        pow2(Mgl)))/(108.*pow4(Msq)) + (((60*(5159*Dmglst1*Mgl + 17171*pow2(
        Dmglst1) + 911*pow2(Mgl)))/pow4(Msq) + (Dmglst1*(3539195 + 255744*B4 -
        3456*DN)*Mgl + (10707109 + 452736*B4 - 8640*DN)*pow2(Dmglst1) + 9*(
        80287 + 13824*B4)*pow2(Mgl))/pow4(Mst2))*pow5(Mst1))/(1296.*pow2(Mgl)))
        + (s2t*pow2(Mt)*(1543500*Dmglst1*Mgl*pow4(Mst1)*(6600*pow2(Msq)*pow2(
        Mst2) + 2*(-17783 + 144*B4 + 72*DN)*pow4(Msq) + 165*pow4(Mst2)) +
        257250*pow2(Dmglst1)*pow4(Mst1)*(99000*pow2(Msq)*pow2(Mst2) + 2*(-
        333689 + 1296*B4 + 648*DN)*pow4(Msq) + 13815*pow4(Mst2)) + pow2(Mgl)*(
        2*pow2(Mst1)*((-1311202751 + 27440000*OepS2 - 1287279000*S2)*pow2(Mst1)
        + 2058*(-7874051 + 14000*OepS2 - 708750*S2)*pow2(Mst2))*pow4(Msq) +
        123480*pow2(Msq)*(20625*pow2(Mst2)*pow4(Mst1) - 2017*pow2(Mst1)*pow4(
        Mst2)) - 675*(187967*pow4(Mst1)*pow4(Mst2) - 334425*pow2(Mst1)*pow6(
        Mst2) + 103583*pow8(Mst2)))))/(6.251175e7*pow2(Mgl)*pow4(Msq)*pow4(
        Mst2))))))/pow4(Mt)/pow2(Sbeta)*12.;
 
   return result;
}

/**
 * 	@return returns the susy log^0 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double himalaya::H5g1::calc_coef_at_as2_no_sm_logs_log0(){

   const double result =
      (pow2(Mt)*pow2(MuSUSY)*((Mt*s2t*pow3(Mst1)*((3*(535 + 60*log(pow2(Msq)/
        pow2(Mst1)) - 60*log(pow2(Mst2)/pow2(Mst1)) + (5*Dmglst1*(Dmglst1 +
        Mgl)*(251 + 12*log(pow2(Msq)/pow2(Mst1)) - 12*log(pow2(Mst2)/pow2(Mst1)
        )))/pow2(Mgl)))/pow4(Msq) + (214598 - 38448*B4 - 216*DN + 64800*log(
        pow2(Msq)/pow2(Mst1)) + 12960*pow2(log(pow2(Msq)/pow2(Mst1))) - 48*log(
        pow2(Mst2)/pow2(Mst1))*(10441 + 2295*log(pow2(Msq)/pow2(Mst1)) + 405*
        pow2(log(pow2(Msq)/pow2(Mst1)))) + 36*(5641 - 180*log(pow2(Msq)/pow2(
        Mst1)))*pow2(log(pow2(Mst2)/pow2(Mst1))) - (6*Dmglst1*(Dmglst1 + Mgl)*(
        -13969 + 11880*B4 - 108*DN + 49128*log(pow2(Mst2)/pow2(Mst1)) + 1080*(-
        2 + 3*log(pow2(Mst2)/pow2(Mst1)))*pow2(log(pow2(Msq)/pow2(Mst1))) -
        10494*pow2(log(pow2(Mst2)/pow2(Mst1))) + 1080*log(pow2(Msq)/pow2(Mst1))
        *(-4 + log(pow2(Mst2)/pow2(Mst1)) + 9*pow2(log(pow2(Mst2)/pow2(Mst1))))
        - 5040*pow3(log(pow2(Mst2)/pow2(Mst1)))))/pow2(Mgl) - 15408*pow3(log(
        pow2(Mst2)/pow2(Mst1))))/pow4(Mst2) + (90*(log(pow2(Mst2)/pow2(Mst1))*(
        pow2(Mgl)*(81*pow2(Mst1)*pow2(Mst2) + 32*pow2(Msq)*(27*pow2(Mst1) + 10*
        pow2(Mst2))) + 3*Dmglst1*Mgl*(159*pow2(Mst1)*pow2(Mst2) + 32*pow2(Msq)*
        (24*pow2(Mst1) + 11*pow2(Mst2))) + 3*pow2(Dmglst1)*(505*pow2(Mst1)*
        pow2(Mst2) + 16*pow2(Msq)*(61*pow2(Mst1) + 22*pow2(Mst2)))) - 2*(pow2(
        Mgl)*(67*pow2(Mst1)*pow2(Mst2) + 16*pow2(Msq)*(14*pow2(Mst1) + 17*pow2(
        Mst2))) + Dmglst1*Mgl*(409*pow2(Mst1)*pow2(Mst2) + 24*pow2(Msq)*(29*
        pow2(Mst1) + 41*pow2(Mst2))) + 4*pow2(Dmglst1)*(348*pow2(Mst1)*pow2(
        Mst2) + pow2(Msq)*(203*pow2(Mst1) + 246*pow2(Mst2)))) + 12*log(pow2(
        Msq)/pow2(Mst1))*(-2*(pow2(Mgl)*(3*pow2(Mst1)*pow2(Mst2) + 8*pow2(Msq)*
        (pow2(Mst1) + pow2(Mst2))) + 3*Dmglst1*Mgl*(5*pow2(Mst1)*pow2(Mst2) +
        8*pow2(Msq)*(pow2(Mst1) + pow2(Mst2))) + 3*pow2(Dmglst1)*(15*pow2(Mst1)
        *pow2(Mst2) + 8*pow2(Msq)*(2*pow2(Mst1) + pow2(Mst2)))) + log(pow2(
        Mst2)/pow2(Mst1))*(pow2(Mgl)*(3*pow2(Mst1)*pow2(Mst2) + 8*pow2(Msq)*(3*
        pow2(Mst1) + pow2(Mst2))) + 3*Dmglst1*Mgl*(5*pow2(Mst1)*pow2(Mst2) + 8*
        pow2(Msq)*(3*pow2(Mst1) + pow2(Mst2))) + 3*pow2(Dmglst1)*(15*pow2(Mst1)
        *pow2(Mst2) + 8*pow2(Msq)*(6*pow2(Mst1) + pow2(Mst2))))) + 6*(-(pow2(
        Mgl)*(8*pow2(Msq) + 3*pow2(Mst1))*pow2(Mst2)) + 3*pow2(Dmglst1)*(8*
        pow2(Msq)*(9*pow2(Mst1) - pow2(Mst2)) - 15*pow2(Mst1)*pow2(Mst2)) + 3*
        Dmglst1*Mgl*(8*pow2(Msq)*(2*pow2(Mst1) - pow2(Mst2)) - 5*pow2(Mst1)*
        pow2(Mst2)))*pow2(log(pow2(Mst2)/pow2(Mst1)))))/(pow2(Mgl)*pow4(Msq)*
        pow4(Mst2))))/243. - (16*pow2(Mst1)*pow2(Mt)*(-417 + 6*B4 + 3*DN + 522*
        log(pow2(Mst2)/pow2(Mst1)) - 144*pow2(log(pow2(Mst2)/pow2(Mst1))) + 2*
        pow3(log(pow2(Mst2)/pow2(Mst1)))))/(27.*pow4(Mst2)) + pow2(s2t)*((2*
        pow2(Mst1)*(836 + 353*log(pow2(Mst2)/pow2(Mst1)) + 30*log(pow2(Msq)/
        pow2(Mst1))*(-7 + 4*log(pow2(Mst2)/pow2(Mst1))) - 435*pow2(log(pow2(
        Mst2)/pow2(Mst1)))))/(135.*pow2(Msq)) + ((-4*pow2(Mst1)*(S2*(2069*pow2(
        Mst1) + 882*pow2(Mst2)) + 6*log(pow2(Mst2)/pow2(Mst1))*(15*pow2(Msq) +
        30*log(pow2(Msq)/pow2(Mst1))*pow2(Msq) + 185*S2*pow2(Mst1) + 102*S2*
        pow2(Mst2))))/27. - ((225024292500*Dmglst1*Mgl - 5334336000*B4*Dmglst1*
        Mgl - 4445280000*Dmglst1*DN*Mgl - 33606882750*pow2(Dmglst1) +
        6223392000*B4*pow2(Dmglst1) + 444528000*DN*pow2(Dmglst1) +
        531403689547*pow2(Mgl) - 16892064000*B4*pow2(Mgl) + 444528000*DN*pow2(
        Mgl) + 1260*log(pow2(Mst2)/pow2(Mst1))*(335899900*Dmglst1*Mgl +
        851392150*pow2(Dmglst1) + 653627493*pow2(Mgl)) + 6667920000*log(pow2(
        Mst2)/pow2(Mst1))*(6*Dmglst1*Mgl + 9*pow2(Dmglst1) + 4*pow2(Mgl))*pow2(
        log(pow2(Msq)/pow2(Mst1))) - 264600*(1328740*Dmglst1*Mgl + 3207610*
        pow2(Dmglst1) + 2258059*pow2(Mgl))*pow2(log(pow2(Mst2)/pow2(Mst1))) -
        555660000*log(pow2(Msq)/pow2(Mst1))*(456*pow2(Dmglst1) - 45*pow2(Mgl) -
        2*log(pow2(Mst2)/pow2(Mst1))*(216*Dmglst1*Mgl + 504*pow2(Dmglst1) +
        121*pow2(Mgl)) + 36*(2*Dmglst1*Mgl + 3*pow2(Dmglst1) + 4*pow2(Mgl))*
        pow2(log(pow2(Mst2)/pow2(Mst1)))) + 74088000*(2418*Dmglst1*Mgl + 4291*
        pow2(Dmglst1) + 3920*pow2(Mgl))*pow3(log(pow2(Mst2)/pow2(Mst1))))*pow4(
        Mst1))/(1.000188e9*pow2(Mgl)) + (32*OepS2*(102*pow2(Mst1)*pow2(Mst2) +
        185*pow4(Mst1)))/729.)/pow4(Mst2) + (-(pow2(Mst1)*(45000*Dmglst1*pow2(
        Mst1)*(5*(67*Dmglst1 + 22*Mgl) - (247*Dmglst1 + 58*Mgl)*log(pow2(Mst2)/
        pow2(Mst1)) - 6*log(pow2(Msq)/pow2(Mst1))*(-5*(5*Dmglst1 + 2*Mgl) + (
        41*Dmglst1 + 14*Mgl)*log(pow2(Mst2)/pow2(Mst1))) + 3*(41*Dmglst1 + 14*
        Mgl)*pow2(log(pow2(Mst2)/pow2(Mst1)))) + pow2(Mgl)*(900*pow2(Mst1)*(
        1525 - 1928*log(pow2(Mst2)/pow2(Mst1)) - 30*log(pow2(Msq)/pow2(Mst1))*(
        -25 + 44*log(pow2(Mst2)/pow2(Mst1))) + 1560*pow2(log(pow2(Mst2)/pow2(
        Mst1)))) + pow2(Msq)*(3228977 - 1026000*B4 + 27000*DN + 33521130*log(
        pow2(Mst2)/pow2(Mst1)) + 202500*(-3 + 8*log(pow2(Mst2)/pow2(Mst1)))*
        pow2(log(pow2(Msq)/pow2(Mst1))) - 23907150*pow2(log(pow2(Mst2)/pow2(
        Mst1))) - 67500*log(pow2(Msq)/pow2(Mst1))*(21 - 94*log(pow2(Mst2)/pow2(
        Mst1)) + 48*pow2(log(pow2(Mst2)/pow2(Mst1)))) + 7812000*pow3(log(pow2(
        Mst2)/pow2(Mst1)))))))/(60750.*pow2(Msq)*pow2(Mst2)) + (-19447774*pow2(
        Mgl)*pow2(Mst1)*pow2(Mst2) - 176400*pow2(Mgl)*pow2(log(pow2(Msq)/pow2(
        Mst1)))*pow2(pow2(Mst1) - pow2(Mst2)) + 249532500*Dmglst1*Mgl*pow4(
        Mst1) + 794130750*pow2(Dmglst1)*pow4(Mst1) + 41220947*pow2(Mgl)*pow4(
        Mst1) - 7399303*pow2(Mgl)*pow4(Mst2) + 44100*pow2(log(pow2(Mst2)/pow2(
        Mst1)))*(700*Dmglst1*Mgl*pow4(Mst1) + 1750*pow2(Dmglst1)*pow4(Mst1) -
        pow2(Mgl)*(216*pow2(Mst1)*pow2(Mst2) + 321*pow4(Mst1) + 32*pow4(Mst2)))
        - 210*log(pow2(Mst2)/pow2(Mst1))*(818300*Dmglst1*Mgl*pow4(Mst1) +
        2133950*pow2(Dmglst1)*pow4(Mst1) - pow2(Mgl)*(126132*pow2(Mst1)*pow2(
        Mst2) + 5151*pow4(Mst1) + 31294*pow4(Mst2))) - 420*log(pow2(Msq)/pow2(
        Mst1))*(14700*Dmglst1*Mgl*pow4(Mst1) + 345450*pow2(Dmglst1)*pow4(Mst1)
        + pow2(Mgl)*(26218*pow2(Mst1)*pow2(Mst2) - 12479*pow4(Mst1) + 9571*
        pow4(Mst2)) + 210*log(pow2(Mst2)/pow2(Mst1))*(700*Dmglst1*Mgl*pow4(
        Mst1) + 1750*pow2(Dmglst1)*pow4(Mst1) - pow2(Mgl)*(104*pow2(Mst1)*pow2(
        Mst2) + 27*pow4(Mst1) + 18*pow4(Mst2)))))/(5.5566e6*pow4(Msq)))/pow2(
        Mgl))) + z4*((4*Mt*MuSUSY*pow2(Mst1)*(-54*Dmglst1*Mgl*Mst1*Mt*(Mt*(-
        507*MuSUSY*s2t + 60*Cbeta*Mst1*s2t*Sbeta) - 4*Cbeta*Sbeta*pow2(Mt) + (
        265*Mst1*MuSUSY + 561*Cbeta*Sbeta*pow2(Mst1) + 561*Cbeta*Sbeta*pow2(
        Mst2))*pow2(s2t)) - 27*Mst1*pow2(Dmglst1)*(-6*s2t*(169*MuSUSY - 30*
        Cbeta*Mst1*Sbeta)*pow2(Mt) + 3*Mt*(185*Mst1*MuSUSY + 989*Cbeta*Sbeta*
        pow2(Mst1) + 374*Cbeta*Sbeta*pow2(Mst2))*pow2(s2t) - 8*Cbeta*Sbeta*
        pow3(Mt) + 120*Cbeta*Mst1*Sbeta*pow2(Mst2)*pow3(s2t)) + pow2(Mgl)*(-2*
        s2t*(8451*Mst1*MuSUSY + 20*Cbeta*Sbeta*pow2(Mst1) + 21*Cbeta*Sbeta*
        pow2(Mst2))*pow2(Mt) + Mt*pow2(s2t)*(-127*MuSUSY*pow2(Mst1) + 39*
        MuSUSY*pow2(Mst2) + 2916*Cbeta*Mst1*Sbeta*pow2(Mst2) + 2916*Cbeta*
        Sbeta*pow3(Mst1)) + 108*(15*MuSUSY + 2*Cbeta*Mst1*Sbeta)*pow3(Mt) +
        Cbeta*Sbeta*pow2(Mst2)*(83*pow2(Mst1) + 1614*pow2(Mst2))*pow3(s2t))))/(
        243.*pow2(Mgl)*pow4(Mst2)) - pow2(Sbeta)*(pow2(Mt)*pow2(s2t)*(2*pow2(
        Mst2) - pow2(Mst1)*((80*Dmglst1)/(3.*Mgl) + (40*pow2(Dmglst1))/pow2(
        Mgl) + (2*(95 - (26*pow2(MuSUSY))/pow2(Mst2)))/81.) - (4*(-(pow2(Mgl)*(
        pow2(Mst2) - 127*pow2(MuSUSY))) + 135*Dmglst1*(111*Dmglst1 + 106*Mgl)*
        pow2(MuSUSY))*pow4(Mst1))/(243.*pow2(Mgl)*pow4(Mst2))) + (80*pow2(Mst1)
        *pow2(MuSUSY)*pow4(Mt))/(3.*pow4(Mst2)) + (Mst1*((-2*Mt*(pow2(Mgl)*(
        169*pow2(Mst1) - 241*pow2(Mst2)) + Dmglst1*Mgl*(989*pow2(Mst1) - 241*
        pow2(Mst2)) + pow2(Dmglst1)*(2219*pow2(Mst1) - 241*pow2(Mst2)))*pow3(
        s2t))/9. + (8*s2t*pow3(Mt)*(4*(Dmglst1*Mgl + pow2(Dmglst1) + pow2(Mgl))
        + ((507*Dmglst1*(Dmglst1 + Mgl) - 313*pow2(Mgl))*pow2(Mst1)*pow2(
        MuSUSY))/pow4(Mst2)))/9.) - ((14310*Dmglst1*Mgl*pow2(Mst1)*(pow2(Mst1)
        - pow2(Mst2)) + 405*pow2(Dmglst1)*(-53*pow2(Mst1)*pow2(Mst2) + 69*pow4(
        Mst1)) + pow2(Mgl)*(-6495*pow2(Mst1)*pow2(Mst2) + 3062*pow4(Mst1) +
        3267*pow4(Mst2)))*pow4(s2t))/486.)/pow2(Mgl))) - (pow2(Mgl)*(-6912*
        Mst1*Mt*s2t*(pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow2(Sbeta)*pow2(z2)*
        pow4(Mst2) + 3456*Mt*pow2(z2)*pow3(Mst1)*(27*Cbeta*Mt*MuSUSY*Sbeta*
        pow2(Mst2)*pow2(s2t) + 44*s2t*pow2(Mt)*pow2(MuSUSY)*(-1 + pow2(Sbeta))
        + 2*Cbeta*MuSUSY*Sbeta*pow3(Mt) - 7*pow2(Sbeta)*pow3(s2t)*pow4(Mst2)) +
        6*s2t*pow2(Mst1)*pow2(Mst2)*(-2*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*pow2(
        s2t)*(236*T1ep + 537*pow2(z2)) - 2*s2t*pow2(Mt)*(pow2(Mst2)*pow2(Sbeta)
        *(20*T1ep + 87*pow2(z2)) + 8*pow2(MuSUSY)*(-1 + pow2(Sbeta))*(68*T1ep +
        213*pow2(z2))) + 56*Cbeta*MuSUSY*Sbeta*(4*T1ep + 3*pow2(z2))*pow3(Mt) +
        pow2(Sbeta)*(100*T1ep + 111*pow2(z2))*pow3(s2t)*pow4(Mst2)) + s2t*pow4(
        Mst1)*(-664*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*pow2(s2t)*(4*T1ep + 3*
        pow2(z2)) - 8*s2t*pow2(Mt)*(-(pow2(Mst2)*pow2(Sbeta)*(4*T1ep + 3*pow2(
        z2))) + 2*pow2(MuSUSY)*(-1 + pow2(Sbeta))*(740*T1ep + 1527*pow2(z2))) +
        320*Cbeta*MuSUSY*Sbeta*(4*T1ep + 3*pow2(z2))*pow3(Mt) - pow2(Sbeta)*(
        44*T1ep + 1113*pow2(z2))*pow3(s2t)*pow4(Mst2)) + 93312*Cbeta*MuSUSY*
        Sbeta*pow2(Mt)*pow2(s2t)*pow2(z2)*pow5(Mst1) + 27*pow2(s2t)*pow2(Sbeta)
        *(-4*pow2(Mt)*(4*T1ep - 5*pow2(z2)) + pow2(Mst2)*pow2(s2t)*(4*T1ep +
        35*pow2(z2)))*pow6(Mst2)) + 432*Dmglst1*Mgl*Mst1*pow2(z2)*(456*Cbeta*
        MuSUSY*Sbeta*pow2(Mt)*pow2(s2t)*pow4(Mst1) - 16*Mt*s2t*(pow2(Mt) +
        pow2(Mst2)*pow2(s2t))*pow2(Sbeta)*pow4(Mst2) - pow2(s2t)*pow3(Mst1)*(8*
        pow2(Mt)*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + pow2(s2t)*pow2(Sbeta)*pow4(
        Mst2)) + 8*Mt*pow2(Mst1)*(57*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*pow2(s2t)
        + 84*s2t*pow2(Mt)*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 2*Cbeta*MuSUSY*
        Sbeta*pow3(Mt) - 17*pow2(Sbeta)*pow3(s2t)*pow4(Mst2)) + Mst1*pow2(
        Sbeta)*pow4(s2t)*pow6(Mst2)) + 216*Mst1*pow2(Dmglst1)*pow2(z2)*(1632*
        Cbeta*MuSUSY*Sbeta*pow2(Mt)*pow2(s2t)*pow4(Mst1) - 32*Mt*s2t*(pow2(Mt)
        + pow2(Mst2)*pow2(s2t))*pow2(Sbeta)*pow4(Mst2) - 3*pow2(s2t)*pow3(Mst1)
        *(8*pow2(Mt)*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + pow2(s2t)*pow2(Sbeta)*
        pow4(Mst2)) + 16*Mt*pow2(Mst1)*(57*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*
        pow2(s2t) + 84*s2t*pow2(Mt)*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 2*Cbeta*
        MuSUSY*Sbeta*pow3(Mt) - 32*pow2(Sbeta)*pow3(s2t)*pow4(Mst2)) + 3*Mst1*
        pow2(Sbeta)*pow4(s2t)*pow6(Mst2)))/(972.*pow2(Mgl)*pow4(Mst2)) + Mt*(
        MuSUSY*(Cbeta*Sbeta*(pow3(s2t)*(pow2(Mst1)*(110.18859670781893 - (40*
        B4)/9. - (2*DN)/9. + (80*log(pow2(Msq)/pow2(Mst1)))/9. + (10*pow2(log(
        pow2(Msq)/pow2(Mst1))))/3. + (log(pow2(Mst2)/pow2(Mst1))*(393289 +
        88500*log(pow2(Msq)/pow2(Mst1)) + 31500*pow2(log(pow2(Msq)/pow2(Mst1)))
        ))/2700. - (13*(627 + 100*log(pow2(Msq)/pow2(Mst1)))*pow2(log(pow2(
        Mst2)/pow2(Mst1))))/60. + (2833*pow3(log(pow2(Mst2)/pow2(Mst1))))/54. +
        (pow2(Mst2)*(-48706000*pow2(Mgl)*pow2(Msq) - 86264500*Dmglst1*Mgl*pow2(
        Mst1) - 269340750*pow2(Dmglst1)*pow2(Mst1) - 20222907*pow2(Mgl)*pow2(
        Mst1) + 3430*log(pow2(Mst2)/pow2(Mst1))*(16700*Dmglst1*Mgl*pow2(Mst1) +
        43550*pow2(Dmglst1)*pow2(Mst1) + pow2(Mgl)*(3256*pow2(Msq) + 2469*pow2(
        Mst1))) + 420*log(pow2(Msq)/pow2(Mst1))*(4900*Dmglst1*Mgl*pow2(Mst1) +
        115150*pow2(Dmglst1)*pow2(Mst1) - pow2(Mgl)*(9800*pow2(Msq) + 12899*
        pow2(Mst1)) + 490*log(pow2(Mst2)/pow2(Mst1))*(100*Dmglst1*Mgl*pow2(
        Mst1) + 250*pow2(Dmglst1)*pow2(Mst1) + pow2(Mgl)*(4*pow2(Msq) + 11*
        pow2(Mst1)))) + 176400*pow2(Mgl)*pow2(Mst1)*pow2(log(pow2(Msq)/pow2(
        Mst1))) + 102900*(-100*Dmglst1*Mgl*pow2(Mst1) - 250*pow2(Dmglst1)*pow2(
        Mst1) + pow2(Mgl)*(76*pow2(Msq) + 15*pow2(Mst1)))*pow2(log(pow2(Mst2)/
        pow2(Mst1)))))/(3.7044e6*pow2(Mgl)*pow4(Msq))) + ((9600*Dmglst1*Mgl +
        24850*pow2(Dmglst1) + 2361*pow2(Mgl) - 25*log(pow2(Mst2)/pow2(Mst1))*(
        254*Dmglst1*Mgl + 701*pow2(Dmglst1) + 63*pow2(Mgl)) - 30*log(pow2(Msq)/
        pow2(Mst1))*(-6*(10*Dmglst1*Mgl + 15*pow2(Dmglst1) + 3*pow2(Mgl)) + 5*
        log(pow2(Mst2)/pow2(Mst1))*(34*Dmglst1*Mgl + 91*pow2(Dmglst1) + 8*pow2(
        Mgl))) + 75*(34*Dmglst1*Mgl + 91*pow2(Dmglst1) + 15*pow2(Mgl))*pow2(
        log(pow2(Mst2)/pow2(Mst1))))*pow4(Mst1))/(135.*pow2(Mgl)*pow2(Msq)) + (
        (pow2(Mst1)*(90*(1 + 2*log(pow2(Msq)/pow2(Mst1)))*pow2(Msq) + S2*(2374*
        pow2(Mst1) + 2655*pow2(Mst2) + 6*log(pow2(Mst2)/pow2(Mst1))*(166*pow2(
        Mst1) + 177*pow2(Mst2)))))/27. + ((311855428500*Dmglst1*Mgl +
        72857573250*pow2(Dmglst1) + 14224896000*B4*pow2(Dmglst1) + 7112448000*
        DN*pow2(Dmglst1) + 478241812219*pow2(Mgl) + 1260*log(pow2(Mst2)/pow2(
        Mst1))*(255108700*Dmglst1*Mgl + 708331750*pow2(Dmglst1) + 215618061*
        pow2(Mgl)) + 10001880000*(2*Dmglst1*Mgl + 3*pow2(Dmglst1) + pow2(Mgl))*
        pow2(log(pow2(Msq)/pow2(Mst1))) - 264600*(703780*Dmglst1*Mgl + 2270170*
        pow2(Dmglst1) + 770503*pow2(Mgl))*pow2(log(pow2(Mst2)/pow2(Mst1))) +
        1666980000*log(pow2(Msq)/pow2(Mst1))*(44*Dmglst1*Mgl - 22*pow2(Dmglst1)
        + 29*pow2(Mgl) + 6*log(pow2(Mst2)/pow2(Mst1))*(16*Dmglst1*Mgl + 44*
        pow2(Dmglst1) + 3*pow2(Mgl)) - 16*pow2(Mgl)*pow2(log(pow2(Mst2)/pow2(
        Mst1)))) + 592704000*(190*Dmglst1*Mgl + 368*pow2(Dmglst1) + 273*pow2(
        Mgl))*pow3(log(pow2(Mst2)/pow2(Mst1))))*pow4(Mst1))/(2.000376e9*pow2(
        Mgl)) - (8*OepS2*(177*pow2(Mst1)*pow2(Mst2) + 166*pow4(Mst1)))/729.)/
        pow2(Mst2) + (-(1 + 2*log(pow2(Mst2)/pow2(Mst1)))/(3.*pow2(Mst1)) + (
        939 - 760*log(pow2(Mst2)/pow2(Mst1)) - 30*log(pow2(Msq)/pow2(Mst1))*(-
        12 + 5*log(pow2(Mst2)/pow2(Mst1))) + 150*pow2(log(pow2(Mst2)/pow2(Mst1)
        )))/(135.*pow2(Msq)) + (pow2(Mst1)*(1.0841585681891803 + (5*Dmglst1)/(
        6.*Mgl) + log(pow2(Msq)/pow2(Mst1))*(0.6291383219954648 - (43*log(pow2(
        Mst2)/pow2(Mst1)))/63.) - (47419*log(pow2(Mst2)/pow2(Mst1)))/26460. + (
        5*pow2(Dmglst1))/(4.*pow2(Mgl)) - pow2(log(pow2(Msq)/pow2(Mst1)))/21. +
        (46*pow2(log(pow2(Mst2)/pow2(Mst1))))/63.))/pow4(Msq))*pow4(Mst2)) + (
        pow3(Mst1)*pow3(Mt)*(-105360*Dmglst1*Mgl*pow2(Msq)*pow2(Mst1) - 281760*
        pow2(Dmglst1)*pow2(Msq)*pow2(Mst1) - 23760*pow2(Mgl)*pow2(Msq)*pow2(
        Mst1) - 105360*Dmglst1*Mgl*pow2(Msq)*pow2(Mst2) - 105360*pow2(Dmglst1)*
        pow2(Msq)*pow2(Mst2) - 23760*pow2(Mgl)*pow2(Msq)*pow2(Mst2) - 52260*
        Dmglst1*Mgl*pow2(Mst1)*pow2(Mst2) - 193200*pow2(Dmglst1)*pow2(Mst1)*
        pow2(Mst2) - 8100*pow2(Mgl)*pow2(Mst1)*pow2(Mst2) - 72*pow2(Msq)*(pow2(
        Dmglst1)*(3037*pow2(Msq) - 2700*pow2(Mst1)) + Dmglst1*Mgl*(3037*pow2(
        Msq) - 900*pow2(Mst1)) + 3*pow2(Mgl)*(293*pow2(Msq) - 60*pow2(Mst1)))*
        pow2(log(pow2(Mst2)/pow2(Mst1))) - 65212*Dmglst1*Mgl*pow4(Msq) - 65212*
        pow2(Dmglst1)*pow4(Msq) + 20964*pow2(Mgl)*pow4(Msq) + 72*(831*Dmglst1*
        Mgl + 831*pow2(Dmglst1) + 265*pow2(Mgl))*pow3(log(pow2(Mst2)/pow2(Mst1)
        ))*pow4(Msq) + 6*log(pow2(Mst2)/pow2(Mst1))*(3*pow2(Mgl)*(120*pow2(Msq)
        *(5*pow2(Mst1) + 7*pow2(Mst2)) - 15*pow2(Mst2)*(2*pow2(Mst1) + 15*pow2(
        Mst2)) + 1286*pow4(Msq)) + Dmglst1*Mgl*(-120*pow2(Msq)*(pow2(Mst1) -
        65*pow2(Mst2)) + 15*(106*pow2(Mst1) - 89*pow2(Mst2))*pow2(Mst2) +
        31778*pow4(Msq)) + pow2(Dmglst1)*(-120*pow2(Msq)*(157*pow2(Mst1) - 65*
        pow2(Mst2)) + 15*(430*pow2(Mst1) - 89*pow2(Mst2))*pow2(Mst2) + 31778*
        pow4(Msq))) + 21135*Dmglst1*Mgl*pow4(Mst2) + 21135*pow2(Dmglst1)*pow4(
        Mst2) + 10035*pow2(Mgl)*pow4(Mst2) - 1080*pow2(log(pow2(Msq)/pow2(Mst1)
        ))*(6*log(pow2(Mst2)/pow2(Mst1))*(Dmglst1*Mgl + pow2(Dmglst1) + pow2(
        Mgl))*pow4(Msq) + (7*Dmglst1*Mgl + 7*pow2(Dmglst1) + 3*pow2(Mgl))*pow4(
        Mst2)) - 180*log(pow2(Msq)/pow2(Mst1))*(-604*Dmglst1*Mgl*pow4(Msq) -
        604*pow2(Dmglst1)*pow4(Msq) - 252*pow2(Mgl)*pow4(Msq) + 144*(3*Dmglst1*
        Mgl + 3*pow2(Dmglst1) + pow2(Mgl))*pow2(log(pow2(Mst2)/pow2(Mst1)))*
        pow4(Msq) + Dmglst1*Mgl*pow4(Mst2) + pow2(Dmglst1)*pow4(Mst2) - 3*pow2(
        Mgl)*pow4(Mst2) - 6*log(pow2(Mst2)/pow2(Mst1))*(3*pow2(Mgl)*(2*pow2(
        Mst1)*pow2(Mst2) + 4*pow2(Msq)*(pow2(Mst1) + pow2(Mst2)) - 2*pow4(Msq)
        + pow4(Mst2)) + Dmglst1*Mgl*(22*pow2(Mst1)*pow2(Mst2) + 28*pow2(Msq)*(
        pow2(Mst1) + pow2(Mst2)) + 66*pow4(Msq) + 7*pow4(Mst2)) + pow2(Dmglst1)
        *(58*pow2(Mst1)*pow2(Mst2) + 4*pow2(Msq)*(13*pow2(Mst1) + 7*pow2(Mst2))
        + 66*pow4(Msq) + 7*pow4(Mst2))))))/(243.*pow2(Mgl)*pow4(Msq)*pow4(Mst2)
        ) - (s2t*pow2(Mt)*((4116*pow2(Mst1)*(30*(2017 + 3750*log(pow2(Mst2)/
        pow2(Mst1)) - 15*log(pow2(Msq)/pow2(Mst1))*(-532 + 75*log(pow2(Mst2)/
        pow2(Mst1))) - 900*pow2(log(pow2(Msq)/pow2(Mst1))) + 3375*pow2(log(
        pow2(Mst2)/pow2(Mst1)))) - (112500*Dmglst1*(5*Dmglst1 + 2*Mgl)*(11 + 6*
        log(pow2(Msq)/pow2(Mst1)) - 6*log(pow2(Mst2)/pow2(Mst1)))*pow2(Mst1) -
        pow2(Mgl)*(-11250*pow2(Mst1)*(55 + 30*log(pow2(Msq)/pow2(Mst1)) - 30*
        log(pow2(Mst2)/pow2(Mst1)) - 9*pow2(log(pow2(Mst2)/pow2(Mst1)))) +
        pow2(Msq)*(7874051 - 7815060*log(pow2(Mst2)/pow2(Mst1)) + 303750*pow2(
        log(pow2(Msq)/pow2(Mst1))) + 1978425*pow2(log(pow2(Mst2)/pow2(Mst1))) -
        101250*log(pow2(Msq)/pow2(Mst1))*(-9 + 6*log(pow2(Mst2)/pow2(Mst1)) +
        2*pow2(log(pow2(Mst2)/pow2(Mst1)))) + 713250*pow3(log(pow2(Mst2)/pow2(
        Mst1))))))/(pow2(Mgl)*pow2(Mst2))))/pow2(Msq) + (343000*pow2(Mst1)*((-
        160*OepS2 + 54*S2*(139 + 60*log(pow2(Mst2)/pow2(Mst1))))*pow2(Mst1) +
        21*(-8*OepS2 + 81*S2*(5 + 2*log(pow2(Mst2)/pow2(Mst1))))*pow2(Mst2)))/
        pow4(Mst2) + (-675*pow2(Mst2)*(334425*pow2(Mst1) - 28*log(pow2(Mst2)/
        pow2(Mst1))*(245*pow2(Mst1) - 5359*pow2(Mst2)) - 103583*pow2(Mst2) +
        420*log(pow2(Msq)/pow2(Mst1))*(735*pow2(Mst1) - 256*pow2(Mst2) + 7*log(
        pow2(Mst2)/pow2(Mst1))*(70*pow2(Mst1) + 37*pow2(Mst2))) - 44100*pow2(
        Mst2)*pow2(log(pow2(Msq)/pow2(Mst1))) - 5880*(35*pow2(Mst1) + 11*pow2(
        Mst2))*pow2(log(pow2(Mst2)/pow2(Mst1)))) + (pow4(Mst1)*(1543500*
        Dmglst1*Mgl*(2*(17783 - 144*B4 - 72*DN + 4860*log(pow2(Msq)/pow2(Mst1))
        - 60*(1061 + 90*log(pow2(Msq)/pow2(Mst1)))*log(pow2(Mst2)/pow2(Mst1)) -
        8019*pow2(log(pow2(Mst2)/pow2(Mst1))) + 5688*pow3(log(pow2(Mst2)/pow2(
        Mst1))))*pow4(Msq) + 15*(-11 + 144*log(pow2(Msq)/pow2(Mst1)) + 60*log(
        pow2(Mst2)/pow2(Mst1)))*pow4(Mst2)) + 257250*pow2(Dmglst1)*(2*(333689 -
        1296*B4 - 648*DN - 498072*log(pow2(Mst2)/pow2(Mst1)) - 8100*log(pow2(
        Msq)/pow2(Mst1))*(-13 + 10*log(pow2(Mst2)/pow2(Mst1))) - 241371*pow2(
        log(pow2(Mst2)/pow2(Mst1))) + 91008*pow3(log(pow2(Mst2)/pow2(Mst1))))*
        pow4(Msq) + 45*(-307 + 936*log(pow2(Msq)/pow2(Mst1)) + 300*log(pow2(
        Mst2)/pow2(Mst1)))*pow4(Mst2)) + pow2(Mgl)*(2*(1311202751 -
        43575800310*log(pow2(Mst2)/pow2(Mst1)) + 1681003800*pow2(log(pow2(Mst2)
        /pow2(Mst1))) - 416745000*log(pow2(Msq)/pow2(Mst1))*(-2 + 5*log(pow2(
        Mst2)/pow2(Mst1)) + pow2(log(pow2(Mst2)/pow2(Mst1)))) + 3533071500*
        pow3(log(pow2(Mst2)/pow2(Mst1))))*pow4(Msq) + 675*(187967 + 471968*log(
        pow2(Mst2)/pow2(Mst1)) - 420*log(pow2(Msq)/pow2(Mst1))*(-2194 + 49*log(
        pow2(Mst2)/pow2(Mst1))) - 44100*pow2(log(pow2(Msq)/pow2(Mst1))) +
        205800*pow2(log(pow2(Mst2)/pow2(Mst1))))*pow4(Mst2))))/(pow2(Mgl)*pow4(
        Mst2)))/pow4(Msq)))/6.251175e7 + Mt*pow2(s2t)*(pow3(Mst1)*((-20*(-249*
        Dmglst1*Mgl - 249*pow2(Dmglst1) - 71*pow2(Mgl) + 12*log(pow2(Msq)/pow2(
        Mst1))*(-2 + log(pow2(Mst2)/pow2(Mst1)))*(3*Dmglst1*Mgl + 3*pow2(
        Dmglst1) + pow2(Mgl)) + 4*log(pow2(Mst2)/pow2(Mst1))*(33*Dmglst1*Mgl +
        33*pow2(Dmglst1) + 10*pow2(Mgl)) - 6*(3*Dmglst1*Mgl + 3*pow2(Dmglst1) +
        pow2(Mgl))*pow2(log(pow2(Mst2)/pow2(Mst1)))))/(9.*pow2(Mgl)*pow2(Msq))
        + (77.49382716049382 + 96*B4 + 40*log(pow2(Msq)/pow2(Mst1)) + log(pow2(
        Mst2)/pow2(Mst1))*(835.1481481481482 + 200*log(pow2(Msq)/pow2(Mst1)) +
        40*pow2(log(pow2(Msq)/pow2(Mst1)))) + (10*(-325 + 36*log(pow2(Msq)/
        pow2(Mst1)))*pow2(log(pow2(Mst2)/pow2(Mst1))))/9. + ((-40*pow2(Mst1)*(
        36*Dmglst1*Mgl + 177*pow2(Dmglst1) + 6*pow2(Mgl) + 2*log(pow2(Mst2)/
        pow2(Mst1))*(39*Dmglst1*Mgl + 24*pow2(Dmglst1) + 17*pow2(Mgl) + 6*log(
        pow2(Msq)/pow2(Mst1))*(3*Dmglst1*Mgl + 6*pow2(Dmglst1) + pow2(Mgl))) +
        3*(9*Dmglst1*Mgl + 33*pow2(Dmglst1) + pow2(Mgl))*pow2(log(pow2(Mst2)/
        pow2(Mst1)))))/(9.*pow2(Msq)) + (Dmglst1*(Dmglst1 + Mgl)*(10021 + 5328*
        B4 - 72*DN + 7560*log(pow2(Msq)/pow2(Mst1)) + 15*log(pow2(Mst2)/pow2(
        Mst1))*(239 - 216*log(pow2(Msq)/pow2(Mst1)) + 72*pow2(log(pow2(Msq)/
        pow2(Mst1)))) + 18*(107 + 300*log(pow2(Msq)/pow2(Mst1)))*pow2(log(pow2(
        Mst2)/pow2(Mst1))) - 3396*pow3(log(pow2(Mst2)/pow2(Mst1)))))/27.)/pow2(
        Mgl) + (136*pow3(log(pow2(Mst2)/pow2(Mst1))))/9.)/pow2(Mst2)) - (5*
        Mst1*pow2(Mst2)*(2*(59 + 12*log(pow2(Msq)/pow2(Mst1)) - 12*log(pow2(
        Mst2)/pow2(Mst1)) + (Dmglst1*(Dmglst1 + Mgl)*(131 + 12*log(pow2(Msq)/
        pow2(Mst1)) - 12*log(pow2(Mst2)/pow2(Mst1))))/pow2(Mgl))*pow2(Mst1) + (
        (-11 - 12*log(pow2(Msq)/pow2(Mst1)) + 12*log(pow2(Mst2)/pow2(Mst1)))*(
        Dmglst1*Mgl + pow2(Dmglst1) + pow2(Mgl))*pow2(Mst2))/pow2(Mgl)))/(108.*
        pow4(Msq)) + ((3539195*Dmglst1*Mgl*pow4(Msq) + 255744*B4*Dmglst1*Mgl*
        pow4(Msq) - 3456*Dmglst1*DN*Mgl*pow4(Msq) + 10707109*pow2(Dmglst1)*
        pow4(Msq) + 452736*B4*pow2(Dmglst1)*pow4(Msq) - 8640*DN*pow2(Dmglst1)*
        pow4(Msq) + 722583*pow2(Mgl)*pow4(Msq) + 124416*B4*pow2(Mgl)*pow4(Msq)
        + 51840*log(pow2(Mst2)/pow2(Mst1))*(Dmglst1*Mgl + pow2(Dmglst1) + pow2(
        Mgl))*pow2(log(pow2(Msq)/pow2(Mst1)))*pow4(Msq) - 288*(889*Dmglst1*Mgl
        + 2206*pow2(Dmglst1) - 69*pow2(Mgl))*pow3(log(pow2(Mst2)/pow2(Mst1)))*
        pow4(Msq) + 12*log(pow2(Mst2)/pow2(Mst1))*(pow2(Dmglst1)*(95783*pow4(
        Msq) - 45510*pow4(Mst2)) + 3*Dmglst1*Mgl*(31467*pow4(Msq) - 4790*pow4(
        Mst2)) + 83*pow2(Mgl)*(1223*pow4(Msq) - 30*pow4(Mst2))) + 309540*
        Dmglst1*Mgl*pow4(Mst2) + 1030260*pow2(Dmglst1)*pow4(Mst2) + 54660*pow2(
        Mgl)*pow4(Mst2) + 72*pow2(log(pow2(Mst2)/pow2(Mst1)))*(pow2(Mgl)*(-
        9079*pow4(Msq) + 90*pow4(Mst2)) + Dmglst1*Mgl*(6853*pow4(Msq) + 450*
        pow4(Mst2)) + 3*pow2(Dmglst1)*(14057*pow4(Msq) + 450*pow4(Mst2))) +
        720*log(pow2(Msq)/pow2(Mst1))*(72*(5*Dmglst1*Mgl + 11*pow2(Dmglst1) +
        pow2(Mgl))*pow2(log(pow2(Mst2)/pow2(Mst1)))*pow4(Msq) + pow2(Mgl)*(-
        126*pow4(Msq) + 37*pow4(Mst2)) + Dmglst1*Mgl*(-270*pow4(Msq) + 181*
        pow4(Mst2)) + pow2(Dmglst1)*(1098*pow4(Msq) + 541*pow4(Mst2)) - 18*log(
        pow2(Mst2)/pow2(Mst1))*(pow2(Mgl)*(-10*pow4(Msq) + pow4(Mst2)) +
        Dmglst1*Mgl*(78*pow4(Msq) + 5*pow4(Mst2)) + pow2(Dmglst1)*(314*pow4(
        Msq) + 15*pow4(Mst2)))))*pow5(Mst1))/(1296.*pow2(Mgl)*pow4(Msq)*pow4(
        Mst2)))) - (z2*(36*Mt*MuSUSY*(-27*pow2(s2t)*(pow2(Mst1)*(20/(3.*pow2(
        Msq)) - (2*(311 - 1080*log(pow2(Msq)/pow2(Mst1)) + 3348*log(pow2(Mst2)/
        pow2(Mst1)) + (180*(14*Dmglst1*Mgl + 41*pow2(Dmglst1) + 2*pow2(Mgl))*
        pow2(Mst1))/(pow2(Mgl)*pow2(Msq))))/(81.*pow2(Mst2))) + (25*(4*Dmglst1*
        Mgl + 10*pow2(Dmglst1) + pow2(Mgl))*pow4(Mst1))/(9.*pow2(Mgl)*pow4(Msq)
        ) + (-(pow2(Mst1)*((-58068*Dmglst1*Mgl + 52254*pow2(Dmglst1) + 69349*
        pow2(Mgl))*pow2(Mst1) + 216*log(pow2(Mst2)/pow2(Mst1))*(966*Dmglst1*
        Mgl*pow2(Mst1) + 1441*pow2(Dmglst1)*pow2(Mst1) + pow2(Mgl)*(-60*pow2(
        Msq) + 498*pow2(Mst1))))) + 6480*log(pow2(Msq)/pow2(Mst1))*(6*Dmglst1*
        Mgl + 9*pow2(Dmglst1) + 2*pow2(Mgl))*pow4(Mst1))/(486.*pow2(Mgl)*pow4(
        Mst2))) + (8*Mt*pow2(Mst1)*(12*Mt*(-55 + 16*log(pow2(Mst2)/pow2(Mst1)))
        - Mst1*s2t*(2407 - 180*log(pow2(Msq)/pow2(Mst1)) - 408*log(pow2(Mst2)/
        pow2(Mst1)) - (3*Dmglst1*(Dmglst1 + Mgl)*(687 + 540*log(pow2(Msq)/pow2(
        Mst1)) - 860*log(pow2(Mst2)/pow2(Mst1))))/pow2(Mgl)) + (15*Mst1*s2t*(
        pow2(Mgl)*(8*pow2(Msq) + 3*pow2(Mst1))*pow2(Mst2) + 3*pow2(Dmglst1)*(
        15*pow2(Mst1)*pow2(Mst2) + 8*pow2(Msq)*(-9*pow2(Mst1) + pow2(Mst2))) +
        3*Dmglst1*Mgl*(5*pow2(Mst1)*pow2(Mst2) + 8*pow2(Msq)*(-2*pow2(Mst1) +
        pow2(Mst2)))))/(pow2(Mgl)*pow4(Msq))))/pow4(Mst2)) + (Cbeta*Sbeta*((8*
        pow2(Mt)*(s2t*(33516*Dmglst1*Mgl*pow2(Mst1) - 117486*pow2(Dmglst1)*
        pow2(Mst1) + 108*log(pow2(Mst2)/pow2(Mst1))*(828*Dmglst1*Mgl*pow2(Mst1)
        + 2618*pow2(Dmglst1)*pow2(Mst1) + pow2(Mgl)*(128*pow2(Mst1) - 37*pow2(
        Mst2))) + pow2(Mgl)*(52798*pow2(Mst1) + 52503*pow2(Mst2))) - (12960*Mt*
        (5*Dmglst1*Mgl + 15*pow2(Dmglst1) + pow2(Mgl))*pow3(Mst1))/pow2(Msq))*
        pow4(Mst1))/(pow2(Mgl)*pow4(Mst2)) + (576*pow3(Mt)*((15*(7*Dmglst1*Mgl
        + 7*pow2(Dmglst1) + 3*pow2(Mgl)))/(pow2(Mgl)*pow4(Msq)) + (1381 + 360*
        log(pow2(Msq)/pow2(Mst1)) - 795*log(pow2(Mst2)/pow2(Mst1)) + (Dmglst1*(
        Dmglst1 + Mgl)*(4163 + 1080*log(pow2(Msq)/pow2(Mst1)) - 2421*log(pow2(
        Mst2)/pow2(Mst1))))/pow2(Mgl))/pow4(Mst2)) - 108*Mt*pow2(s2t)*((-4*(263
        - 180*log(pow2(Msq)/pow2(Mst1)) + 84*log(pow2(Mst2)/pow2(Mst1)) + (-3*
        Dmglst1*(Dmglst1 + Mgl)*(657 + 300*log(pow2(Msq)/pow2(Mst1)) - 526*log(
        pow2(Mst2)/pow2(Mst1))) + (60*(9*Dmglst1*Mgl + 33*pow2(Dmglst1) + pow2(
        Mgl))*pow2(Mst1))/pow2(Msq))/pow2(Mgl)))/pow2(Mst2) + ((240*(3*Dmglst1*
        Mgl + 3*pow2(Dmglst1) + pow2(Mgl)))/pow2(Msq) + (pow2(Mst1)*(6487*
        Dmglst1*Mgl*pow4(Msq) + 20793*pow2(Dmglst1)*pow4(Msq) - 1789*pow2(Mgl)*
        pow4(Msq) + 720*log(pow2(Msq)/pow2(Mst1))*(5*Dmglst1*Mgl + 11*pow2(
        Dmglst1) + pow2(Mgl))*pow4(Msq) - 12*log(pow2(Mst2)/pow2(Mst1))*(953*
        Dmglst1*Mgl + 2402*pow2(Dmglst1) + 99*pow2(Mgl))*pow4(Msq) + 450*
        Dmglst1*Mgl*pow4(Mst2) + 1350*pow2(Dmglst1)*pow4(Mst2) + 90*pow2(Mgl)*
        pow4(Mst2)))/(pow4(Msq)*pow4(Mst2)))/pow2(Mgl)))*pow5(Mst1) - pow3(s2t)
        *(-12*(2740 + 1350*log(pow2(Msq)/pow2(Mst1)) - 4689*log(pow2(Mst2)/
        pow2(Mst1)) + (540*pow2(Msq))/pow2(Mst2) + (270*pow2(Mst2))/pow2(Msq))*
        pow4(Mst1) + 648*pow4(Mst2) + ((1728*log(pow2(Mst2)/pow2(Mst1))*(42*
        Dmglst1*Mgl + 62*pow2(Dmglst1) + 39*pow2(Mgl))*pow4(Msq) + pow2(Mgl)*(
        7560*pow2(Msq)*pow2(Mst2) + 65617*pow4(Msq) - 1350*pow4(Mst2)) + 18*
        pow2(Dmglst1)*(5460*pow2(Msq)*pow2(Mst2) + 16943*pow4(Msq) - 750*pow4(
        Mst2)) + 36*Dmglst1*Mgl*(1020*pow2(Msq)*pow2(Mst2) + 3067*pow4(Msq) -
        150*pow4(Mst2)))*pow6(Mst1))/(pow2(Mgl)*pow2(Mst2)*pow4(Msq)))))/pow2(
        Mst1)))/972.) + (Mst1*pow2(Sbeta)*pow3(s2t)*(-1297272*Dmglst1*Mgl*pow2(
        Mst1)*pow2(Mst2)*pow4(Msq) - 226368*B4*Dmglst1*Mgl*pow2(Mst1)*pow2(
        Mst2)*pow4(Msq) + 4320*Dmglst1*DN*Mgl*pow2(Mst1)*pow2(Mst2)*pow4(Msq) -
        1370424*pow2(Dmglst1)*pow2(Mst1)*pow2(Mst2)*pow4(Msq) - 423360*B4*pow2(
        Dmglst1)*pow2(Mst1)*pow2(Mst2)*pow4(Msq) + 9504*DN*pow2(Dmglst1)*pow2(
        Mst1)*pow2(Mst2)*pow4(Msq) - 1059256*pow2(Mgl)*pow2(Mst1)*pow2(Mst2)*
        pow4(Msq) - 95040*B4*pow2(Mgl)*pow2(Mst1)*pow2(Mst2)*pow4(Msq) + 864*
        DN*pow2(Mgl)*pow2(Mst1)*pow2(Mst2)*pow4(Msq) - 25920*(Dmglst1*Mgl +
        pow2(Dmglst1) + pow2(Mgl))*pow2(Mst2)*(2*(pow2(Mst1) - pow2(Mst2)) +
        log(pow2(Mst2)/pow2(Mst1))*(pow2(Mst1) + pow2(Mst2)))*pow2(log(pow2(
        Msq)/pow2(Mst1)))*pow4(Msq) + 924480*Dmglst1*Mgl*pow2(Msq)*pow2(Mst2)*
        pow4(Mst1) + 2632320*pow2(Dmglst1)*pow2(Msq)*pow2(Mst2)*pow4(Mst1) +
        239040*pow2(Mgl)*pow2(Msq)*pow2(Mst2)*pow4(Mst1) - 3058187*Dmglst1*Mgl*
        pow4(Msq)*pow4(Mst1) - 9306949*pow2(Dmglst1)*pow4(Msq)*pow4(Mst1) -
        622151*pow2(Mgl)*pow4(Msq)*pow4(Mst1) - 725760*Dmglst1*Mgl*pow2(Msq)*
        pow2(Mst1)*pow4(Mst2) - 1621440*pow2(Dmglst1)*pow2(Msq)*pow2(Mst1)*
        pow4(Mst2) - 213120*pow2(Mgl)*pow2(Msq)*pow2(Mst1)*pow4(Mst2) + 816264*
        Dmglst1*Mgl*pow4(Msq)*pow4(Mst2) - 29376*B4*Dmglst1*Mgl*pow4(Msq)*pow4(
        Mst2) - 864*Dmglst1*DN*Mgl*pow4(Msq)*pow4(Mst2) - 29736*pow2(Dmglst1)*
        pow4(Msq)*pow4(Mst2) - 29376*B4*pow2(Dmglst1)*pow4(Msq)*pow4(Mst2) -
        864*DN*pow2(Dmglst1)*pow4(Msq)*pow4(Mst2) + 958824*pow2(Mgl)*pow4(Msq)*
        pow4(Mst2) - 29376*B4*pow2(Mgl)*pow4(Msq)*pow4(Mst2) - 864*DN*pow2(Mgl)
        *pow4(Msq)*pow4(Mst2) - 325260*Dmglst1*Mgl*pow4(Mst1)*pow4(Mst2) -
        1058940*pow2(Dmglst1)*pow4(Mst1)*pow4(Mst2) - 61740*pow2(Mgl)*pow4(
        Mst1)*pow4(Mst2) + 288*pow3(log(pow2(Mst2)/pow2(Mst1)))*pow4(Msq)*(
        Dmglst1*Mgl*(712*pow2(Mst1)*pow2(Mst2) + 323*pow4(Mst1) - 146*pow4(
        Mst2)) + pow2(Dmglst1)*(1663*pow2(Mst1)*pow2(Mst2) + 689*pow4(Mst1) -
        146*pow4(Mst2)) - pow2(Mgl)*(-78*pow2(Mst1)*pow2(Mst2) + pow4(Mst1) +
        146*pow4(Mst2))) - 72*pow2(log(pow2(Mst2)/pow2(Mst1)))*(3*pow2(Dmglst1)
        *(pow4(Msq)*(6738*pow2(Mst1)*pow2(Mst2) + 8913*pow4(Mst1) - 1594*pow4(
        Mst2)) + 450*pow4(Mst1)*pow4(Mst2) - 240*pow2(Msq)*(13*pow2(Mst2)*pow4(
        Mst1) - 2*pow2(Mst1)*pow4(Mst2))) - pow2(Mgl)*(-90*pow4(Mst1)*pow4(
        Mst2) + pow4(Msq)*(1718*pow2(Mst1)*pow2(Mst2) + 2579*pow4(Mst1) + 4782*
        pow4(Mst2)) + 240*pow2(Msq)*(2*pow2(Mst2)*pow4(Mst1) - pow2(Mst1)*pow4(
        Mst2))) + Dmglst1*Mgl*(pow4(Msq)*(6066*pow2(Mst1)*pow2(Mst2) + 5569*
        pow4(Mst1) - 4782*pow4(Mst2)) + 450*pow4(Mst1)*pow4(Mst2) - 720*pow2(
        Msq)*(4*pow2(Mst2)*pow4(Mst1) - pow2(Mst1)*pow4(Mst2)))) + 8640*
        Dmglst1*Mgl*pow2(Msq)*pow6(Mst2) + 8640*pow2(Dmglst1)*pow2(Msq)*pow6(
        Mst2) + 8640*pow2(Mgl)*pow2(Msq)*pow6(Mst2) + 16380*Dmglst1*Mgl*pow2(
        Mst1)*pow6(Mst2) + 29340*pow2(Dmglst1)*pow2(Mst1)*pow6(Mst2) + 7740*
        pow2(Mgl)*pow2(Mst1)*pow6(Mst2) + 2160*log(pow2(Msq)/pow2(Mst1))*(12*
        pow2(Mst2)*(pow2(Dmglst1)*(-23*pow2(Mst1) + pow2(Mst2)) + Dmglst1*Mgl*(
        -11*pow2(Mst1) + pow2(Mst2)) + pow2(Mgl)*(-3*pow2(Mst1) + pow2(Mst2)))*
        pow2(log(pow2(Mst2)/pow2(Mst1)))*pow4(Msq) + 2*log(pow2(Mst2)/pow2(
        Mst1))*(pow2(Mgl)*(8*pow2(Msq)*pow2(Mst1)*pow2(Mst2)*(pow2(Mst1) +
        pow2(Mst2)) + 6*pow4(Msq)*(-3*pow2(Mst1)*pow2(Mst2) + 5*pow4(Mst1) - 7*
        pow4(Mst2)) + 3*pow4(Mst1)*pow4(Mst2)) + 3*Dmglst1*Mgl*(8*pow2(Msq)*
        pow2(Mst1)*pow2(Mst2)*(pow2(Mst1) + pow2(Mst2)) + 2*pow4(Msq)*(13*pow2(
        Mst1)*pow2(Mst2) + 33*pow4(Mst1) - 7*pow4(Mst2)) + 5*pow4(Mst1)*pow4(
        Mst2)) + pow2(Dmglst1)*(48*pow2(Msq)*pow2(Mst1)*pow2(Mst2)*(pow2(Mst1)
        + pow2(Mst2)) + 6*pow4(Msq)*(53*pow2(Mst1)*pow2(Mst2) + 111*pow4(Mst1)
        - 7*pow4(Mst2)) + 45*pow4(Mst1)*pow4(Mst2))) + pow2(Mgl)*(32*pow2(Msq)*
        pow2(Mst1)*(pow2(Mst1) - pow2(Mst2))*pow2(Mst2) - 13*pow4(Mst1)*pow4(
        Mst2) + 6*pow4(Msq)*(-28*pow2(Mst1)*pow2(Mst2) + 11*pow4(Mst1) + 24*
        pow4(Mst2)) + pow2(Mst1)*pow6(Mst2)) + Dmglst1*Mgl*(96*pow2(Msq)*pow2(
        Mst1)*(pow2(Mst1) - pow2(Mst2))*pow2(Mst2) - 61*pow4(Mst1)*pow4(Mst2) +
        6*pow4(Msq)*(-64*pow2(Mst1)*pow2(Mst2) + 43*pow4(Mst1) + 36*pow4(Mst2))
        + pow2(Mst1)*pow6(Mst2)) + pow2(Dmglst1)*(192*pow2(Msq)*pow2(Mst1)*(
        pow2(Mst1) - pow2(Mst2))*pow2(Mst2) - 181*pow4(Mst1)*pow4(Mst2) + 2*
        pow4(Msq)*(-502*pow2(Mst1)*pow2(Mst2) + 177*pow4(Mst1) + 142*pow4(Mst2)
        ) + pow2(Mst1)*pow6(Mst2))) - 12*log(pow2(Mst2)/pow2(Mst1))*(pow2(
        Dmglst1)*(90*pow2(Mst1)*(-507*pow2(Mst1) + 2*pow2(Mst2))*pow4(Mst2) +
        pow4(Msq)*(-250152*pow2(Mst1)*pow2(Mst2) + 254699*pow4(Mst1) + 91236*
        pow4(Mst2)) + 1440*pow2(Msq)*(29*pow2(Mst2)*pow4(Mst1) - 45*pow2(Mst1)*
        pow4(Mst2))) + pow2(Mgl)*(90*pow2(Mst1)*(-29*pow2(Mst1) + 2*pow2(Mst2))
        *pow4(Mst2) + pow4(Msq)*(13336*pow2(Mst1)*pow2(Mst2) + 11313*pow4(Mst1)
        + 76860*pow4(Mst2)) - 960*pow2(Msq)*(7*pow2(Mst2)*pow4(Mst1) + 10*pow2(
        Mst1)*pow4(Mst2))) + 3*Dmglst1*Mgl*(-4830*pow4(Mst1)*pow4(Mst2) + pow4(
        Msq)*(-23192*pow2(Mst1)*pow2(Mst2) + 26687*pow4(Mst1) + 27972*pow4(
        Mst2)) - 960*pow2(Msq)*(2*pow2(Mst2)*pow4(Mst1) + 11*pow2(Mst1)*pow4(
        Mst2)) + 60*pow2(Mst1)*pow6(Mst2)))))/(3888.*pow2(Mgl)*pow2(Mst2)*pow4(
        Msq))) + (z3*((32*pow2(Mst1)*pow2(Mt)*pow2(MuSUSY)*(216*Dmglst1*Mgl*
        Mst1*s2t*(-231*Mt + 382*Mst1*s2t) + 27*Mst1*s2t*(-1848*Mt + 6491*Mst1*
        s2t)*pow2(Dmglst1) + pow2(Mgl)*(-62100*Mst1*Mt*s2t + 4968*pow2(Mt) + (
        67606*pow2(Mst1) + 26349*pow2(Mst2))*pow2(s2t)) - 162*log(pow2(Mst2)/
        pow2(Mst1))*(2*Dmglst1*Mgl*Mst1*s2t*(225*Mt - 23*Mst1*s2t) + Mst1*s2t*(
        450*Mt - 53*Mst1*s2t)*pow2(Dmglst1) + pow2(Mgl)*(138*Mst1*Mt*s2t + 4*
        pow2(Mt) - 21*(pow2(Mst1) + pow2(Mst2))*pow2(s2t)))))/pow4(Mst2) + (16*
        Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst1)*(-27*Dmglst1*Mgl*Mst1*(320*Mst1*s2t*
        pow2(Mt) + 3*Mt*(2183*pow2(Mst1) - 808*pow2(Mst2))*pow2(s2t) - 3760*
        pow3(Mt) + 2846*Mst1*pow2(Mst2)*pow3(s2t)) - 27*Mst1*pow2(Dmglst1)*(
        2464*Mst1*s2t*pow2(Mt) + 3*Mt*(9453*pow2(Mst1) - 808*pow2(Mst2))*pow2(
        s2t) - 3760*pow3(Mt) + 6197*Mst1*pow2(Mst2)*pow3(s2t)) + pow2(Mgl)*(
        4260*s2t*pow2(Mst2)*pow2(Mt) + 29565*Mt*pow2(s2t)*pow3(Mst1) + 324*
        Mst1*(221*Mt*pow2(Mst2)*pow2(s2t) + 86*pow3(Mt)) + pow2(Mst1)*(2936*
        s2t*pow2(Mt) - 41257*pow2(Mst2)*pow3(s2t)) - 21921*pow3(s2t)*pow4(Mst2)
        ) + 162*log(pow2(Mst2)/pow2(Mst1))*(4*Dmglst1*Mgl*Mst1*Mt*(4*Mst1*Mt*
        s2t + 4*pow2(Mt) + 171*(pow2(Mst1) + pow2(Mst2))*pow2(s2t)) + 2*Mst1*
        pow2(Dmglst1)*(12*Mst1*s2t*pow2(Mt) + 9*Mt*(77*pow2(Mst1) + 38*pow2(
        Mst2))*pow2(s2t) + 8*pow3(Mt) + 8*Mst1*pow2(Mst2)*pow3(s2t)) + pow2(
        Mgl)*(216*Mt*pow2(s2t)*pow3(Mst1) + 8*Mst1*(27*Mt*pow2(Mst2)*pow2(s2t)
        + 2*pow3(Mt)) - 21*pow3(s2t)*pow4(Mst2)))))/pow4(Mst2) + pow2(Sbeta)*((
        -648*log(pow2(Mst2)/pow2(Mst1))*pow2(Mst1)*(pow2(Dmglst1)*(85*pow2(
        Mst1) - 69*pow2(Mst2)) + 46*Dmglst1*Mgl*(pow2(Mst1) - pow2(Mst2)) + 21*
        pow2(Mgl)*(pow2(Mst1) - pow2(Mst2))) + 27*pow2(Dmglst1)*(4242*pow2(
        Mst1)*pow2(Mst2) + 22079*pow4(Mst1) - 1533*pow4(Mst2)) + 54*Dmglst1*
        Mgl*(810*pow2(Mst1)*pow2(Mst2) + 5077*pow4(Mst1) - 195*pow4(Mst2)) + 4*
        pow2(Mgl)*(17493*pow2(Mst1)*pow2(Mst2) + 19336*pow4(Mst1) + 4428*pow4(
        Mst2)))*pow4(s2t) + (4*Mt*(108*Mst1*pow2(Mst2)*pow3(s2t)*(18*log(pow2(
        Mst2)/pow2(Mst1))*pow2(Mst2)*(pow2(Dmglst1)*(-155*pow2(Mst1) + pow2(
        Mst2)) + Dmglst1*Mgl*(-77*pow2(Mst1) + pow2(Mst2)) + pow2(Mgl)*(-25*
        pow2(Mst1) + pow2(Mst2))) + pow2(Mgl)*(-618*pow2(Mst1)*pow2(Mst2) +
        519*pow4(Mst1) - 266*pow4(Mst2)) + Dmglst1*Mgl*(-692*pow2(Mst1)*pow2(
        Mst2) + 2991*pow4(Mst1) - 116*pow4(Mst2)) + pow2(Dmglst1)*(-665*pow2(
        Mst1)*pow2(Mst2) + 9559*pow4(Mst1) + 559*pow4(Mst2))) + Mt*(54*pow2(Mt)
        *(pow2(Dmglst1)*(18608*pow2(Mst1)*pow2(Mst2) + 26736*pow4(Mst1) - 2019*
        pow4(Mst2)) + 2*Dmglst1*Mgl*(3088*pow2(Mst1)*pow2(Mst2) + 5264*pow4(
        Mst1) - 107*pow4(Mst2)) + 8*pow2(Mgl)*(4*pow2(Mst1)*(36*pow2(Mst2) -
        23*pow2(MuSUSY)) + 317*pow4(Mst1) - 163*pow4(Mst2) - 6*log(pow2(Mst2)/
        pow2(Mst1))*(-2*pow2(Mst1)*pow2(MuSUSY) + pow4(Mst2)))) - 216*Mst1*Mt*
        s2t*(4*pow2(Mgl)*(pow2(Mst1)*(183*pow2(Mst2) - 575*pow2(MuSUSY)) + 448*
        pow4(Mst1) - 54*pow4(Mst2)) + 4*Dmglst1*Mgl*(pow2(Mst1)*(507*pow2(Mst2)
        - 462*pow2(MuSUSY)) + 2196*pow4(Mst1) - 37*pow4(Mst2)) + pow2(Dmglst1)*
        (84*pow2(Mst1)*(45*pow2(Mst2) - 22*pow2(MuSUSY)) + 25672*pow4(Mst1) +
        1695*pow4(Mst2)) + 12*log(pow2(Mst2)/pow2(Mst1))*(Dmglst1*Mgl*(-225*
        pow2(Mst1)*pow2(MuSUSY) + 4*pow4(Mst2)) + pow2(Dmglst1)*(-225*pow2(
        Mst1)*pow2(MuSUSY) + 4*pow4(Mst2)) + pow2(Mgl)*(-69*pow2(Mst1)*pow2(
        MuSUSY) + 4*pow4(Mst2)))) - pow2(s2t)*(1296*log(pow2(Mst2)/pow2(Mst1))*
        pow2(Mst1)*(21*pow2(Mgl)*(pow2(Mst1) + pow2(Mst2))*pow2(MuSUSY) +
        Dmglst1*Mgl*(46*pow2(Mst1)*pow2(MuSUSY) + 4*pow4(Mst2)) + pow2(Dmglst1)
        *(53*pow2(Mst1)*pow2(MuSUSY) + 6*pow4(Mst2))) - 216*Dmglst1*Mgl*((398*
        pow2(Mst2) - 3056*pow2(MuSUSY))*pow4(Mst1) - 373*pow2(Mst1)*pow4(Mst2)
        + 55*pow6(Mst2)) + 4*pow2(Mgl)*((-662*pow2(Mst2) + 135212*pow2(MuSUSY))
        *pow4(Mst1) + pow2(Mst1)*(52698*pow2(Mst2)*pow2(MuSUSY) + 969*pow4(
        Mst2)) + 1161*pow6(Mst2)) - 27*pow2(Dmglst1)*(8*(2189*pow2(Mst2) -
        6491*pow2(MuSUSY))*pow4(Mst1) - 16063*pow2(Mst1)*pow4(Mst2) + 3479*
        pow6(Mst2))))))/pow4(Mst2))))/(1944.*pow2(Mgl)) + pow2(Sbeta)*(-(((5*(1
        + 2*log(pow2(Msq)/pow2(Mst1)))*(1 - 2*log(pow2(Mst2)/pow2(Mst1)))*pow2(
        Msq)*pow2(Mst1))/6. + ((-6667920000*pow2(Mgl)*pow2(Msq) + 397908640500*
        Dmglst1*Mgl*pow2(Mst2) + 5334336000*B4*Dmglst1*Mgl*pow2(Mst2) +
        4445280000*Dmglst1*DN*Mgl*pow2(Mst2) + 59725475250*pow2(Dmglst1)*pow2(
        Mst2) + 22226400000*B4*pow2(Dmglst1)*pow2(Mst2) + 13780368000*DN*pow2(
        Dmglst1)*pow2(Mst2) + 1260*log(pow2(Mst2)/pow2(Mst1))*(154384300*
        Dmglst1*Mgl + 526022350*pow2(Dmglst1) - 15635871*pow2(Mgl))*pow2(Mst2)
        + 257823187891*pow2(Mgl)*pow2(Mst2) + 8890560000*B4*pow2(Mgl)*pow2(
        Mst2) + 444528000*DN*pow2(Mgl)*pow2(Mst2) - 3333960000*(log(pow2(Mst2)/
        pow2(Mst1))*(12*Dmglst1*Mgl + 18*pow2(Dmglst1) + 7*pow2(Mgl)) - pow2(3*
        Dmglst1 + Mgl))*pow2(Mst2)*pow2(log(pow2(Msq)/pow2(Mst1))) - 264600*(
        86380*Dmglst1*Mgl + 1344070*pow2(Dmglst1) - 256523*pow2(Mgl))*pow2(
        Mst2)*pow2(log(pow2(Mst2)/pow2(Mst1))) - 555660000*log(pow2(Msq)/pow2(
        Mst1))*(pow2(Mgl)*(24*pow2(Msq) - 55*pow2(Mst2)) - 324*Dmglst1*Mgl*
        pow2(Mst2) - 458*pow2(Dmglst1)*pow2(Mst2) - 16*log(pow2(Mst2)/pow2(
        Mst1))*(9*Dmglst1*Mgl + 36*pow2(Dmglst1) - 4*pow2(Mgl))*pow2(Mst2) - 6*
        (12*Dmglst1*Mgl + 18*pow2(Dmglst1) + 5*pow2(Mgl))*pow2(Mst2)*pow2(log(
        pow2(Mst2)/pow2(Mst1)))) + 37044000*(1244*Dmglst1*Mgl + 3194*pow2(
        Dmglst1) + 1535*pow2(Mgl))*pow2(Mst2)*pow3(log(pow2(Mst2)/pow2(Mst1))))
        *pow4(Mst1))/(8.001504e9*pow2(Mgl)*pow2(Mst2)) + pow2(Mst2)*((5*(1 + 2*
        log(pow2(Msq)/pow2(Mst1)))*(1 + 2*log(pow2(Mst2)/pow2(Mst1)))*pow2(Msq)
        )/6. - (pow2(Mst1)*(2589750*Dmglst1*Mgl + 162000*B4*Dmglst1*Mgl +
        135000*Dmglst1*DN*Mgl - 4030875*pow2(Dmglst1) + 243000*B4*pow2(Dmglst1)
        + 202500*DN*pow2(Dmglst1) - 11753176*pow2(Mgl) + 27000*B4*pow2(Mgl) +
        40500*DN*pow2(Mgl) - 540*log(pow2(Mst2)/pow2(Mst1))*(8550*Dmglst1*Mgl +
        15700*pow2(Dmglst1) + 1661*pow2(Mgl)) - 101250*(6*Dmglst1*Mgl + 9*pow2(
        Dmglst1) + 7*pow2(Mgl) + 6*log(pow2(Mst2)/pow2(Mst1))*(2*Dmglst1*Mgl +
        3*pow2(Dmglst1) + pow2(Mgl)))*pow2(log(pow2(Msq)/pow2(Mst1))) + 4050*(
        1210*Dmglst1*Mgl + 1815*pow2(Dmglst1) + 1124*pow2(Mgl))*pow2(log(pow2(
        Mst2)/pow2(Mst1))) + 33750*log(pow2(Msq)/pow2(Mst1))*(126*Dmglst1*Mgl +
        329*pow2(Dmglst1) - 53*pow2(Mgl) - 12*log(pow2(Mst2)/pow2(Mst1))*(6*
        Dmglst1*Mgl + 9*pow2(Dmglst1) + 2*pow2(Mgl)) + 6*(6*Dmglst1*Mgl + 9*
        pow2(Dmglst1) + 5*pow2(Mgl))*pow2(log(pow2(Mst2)/pow2(Mst1)))) - 2250*(
        898*Dmglst1*Mgl + 1347*pow2(Dmglst1) + 1097*pow2(Mgl))*pow3(log(pow2(
        Mst2)/pow2(Mst1)))))/(243000.*pow2(Mgl)) + ((6850*Dmglst1*Mgl + 16475*
        pow2(Dmglst1) + 2068*pow2(Mgl) - log(pow2(Mst2)/pow2(Mst1))*(4900*
        Dmglst1*Mgl + 11350*pow2(Dmglst1) + 991*pow2(Mgl)) - 15*log(pow2(Msq)/
        pow2(Mst1))*(-20*Dmglst1*Mgl + 70*pow2(Dmglst1) - 23*pow2(Mgl) + log(
        pow2(Mst2)/pow2(Mst1))*(200*Dmglst1*Mgl + 500*pow2(Dmglst1) + 41*pow2(
        Mgl))) + 30*(50*Dmglst1*Mgl + 125*pow2(Dmglst1) + 14*pow2(Mgl))*pow2(
        log(pow2(Mst2)/pow2(Mst1))))*pow4(Mst1))/(270.*pow2(Mgl)*pow2(Msq))) +
        (2*OepS2*(-150*pow2(Mst1)*pow2(Mst2) + 11*pow4(Mst1) - 27*pow4(Mst2)))/
        729. - ((53749 + 2592*B4 - 288*DN + 13320*log(pow2(Msq)/pow2(Mst1)) + (
        2160*(1 + 2*log(pow2(Msq)/pow2(Mst1)))*pow2(Msq))/pow2(Mst1) + 5400*
        pow2(log(pow2(Msq)/pow2(Mst1))) + 54*(731 + 60*log(pow2(Msq)/pow2(Mst1)
        ))*pow2(log(pow2(Mst2)/pow2(Mst1))) + (36*Dmglst1*(7 - 300*log(pow2(
        Msq)/pow2(Mst1)) + 226*log(pow2(Mst2)/pow2(Mst1)) + 180*pow2(log(pow2(
        Msq)/pow2(Mst1))) + 18*pow2(log(pow2(Mst2)/pow2(Mst1)))))/Mgl + 6*(-(
        log(pow2(Mst2)/pow2(Mst1))*(14209 + 2100*log(pow2(Msq)/pow2(Mst1)) +
        180*pow2(log(pow2(Msq)/pow2(Mst1))))) + (pow2(Dmglst1)*(6457 - 4020*
        log(pow2(Msq)/pow2(Mst1)) + 2670*log(pow2(Mst2)/pow2(Mst1)) + 1620*
        pow2(log(pow2(Msq)/pow2(Mst1))) + 162*pow2(log(pow2(Mst2)/pow2(Mst1))))
        )/pow2(Mgl)) - 7668*pow3(log(pow2(Mst2)/pow2(Mst1))) - (3*pow2(Mst1)*(-
        56252000*Dmglst1*Mgl*pow2(Msq) - 111132000*pow2(Dmglst1)*pow2(Msq) -
        37236080*pow2(Mgl)*pow2(Msq) - 44675750*Dmglst1*Mgl*pow2(Mst1) -
        136985625*pow2(Dmglst1)*pow2(Mst1) - 12119532*pow2(Mgl)*pow2(Mst1) +
        35*log(pow2(Mst2)/pow2(Mst1))*(4900*Dmglst1*Mgl*(276*pow2(Msq) + 167*
        pow2(Mst1)) + 2450*pow2(Dmglst1)*(828*pow2(Msq) + 871*pow2(Mst1)) +
        pow2(Mgl)*(457464*pow2(Msq) + 215819*pow2(Mst1))) + 420*log(pow2(Msq)/
        pow2(Mst1))*(2450*Dmglst1*Mgl*(16*pow2(Msq) + pow2(Mst1)) + 1225*pow2(
        Dmglst1)*(128*pow2(Msq) + 47*pow2(Mst1)) - 4*pow2(Mgl)*(4165*pow2(Msq)
        + 2306*pow2(Mst1)) + 35*log(pow2(Mst2)/pow2(Mst1))*(140*Dmglst1*Mgl*(6*
        pow2(Msq) + 5*pow2(Mst1)) + 70*pow2(Dmglst1)*(18*pow2(Msq) + 25*pow2(
        Mst1)) + pow2(Mgl)*(168*pow2(Msq) + 163*pow2(Mst1)))) + 176400*pow2(
        Mgl)*pow2(Mst1)*pow2(log(pow2(Msq)/pow2(Mst1))) - 7350*(140*Dmglst1*
        Mgl*(6*pow2(Msq) + 5*pow2(Mst1)) + 70*pow2(Dmglst1)*(18*pow2(Msq) + 25*
        pow2(Mst1)) + pow2(Mgl)*(-252*pow2(Msq) + 79*pow2(Mst1)))*pow2(log(
        pow2(Mst2)/pow2(Mst1)))))/(8575.*pow2(Mgl)*pow4(Msq)))*pow4(Mst2))/
        2592. + (S2*(3546*pow2(Mst1)*pow2(Mst2) - 281*pow4(Mst1) - 891*pow4(
        Mst2) + log(pow2(Mst2)/pow2(Mst1))*(900*pow2(Mst1)*pow2(Mst2) - 66*
        pow4(Mst1) + 162*pow4(Mst2))))/108.)*pow4(s2t)) + (Mst1*s2t*pow3(Mt)*(
        177120*Dmglst1*Mgl*pow2(Msq)*pow2(Mst1)*pow2(Mst2)*pow2(MuSUSY) +
        25920*Dmglst1*Mgl*z2*pow2(Msq)*pow2(Mst1)*pow2(Mst2)*pow2(MuSUSY) +
        177120*pow2(Dmglst1)*pow2(Msq)*pow2(Mst1)*pow2(Mst2)*pow2(MuSUSY) +
        25920*z2*pow2(Dmglst1)*pow2(Msq)*pow2(Mst1)*pow2(Mst2)*pow2(MuSUSY) +
        48960*pow2(Mgl)*pow2(Msq)*pow2(Mst1)*pow2(Mst2)*pow2(MuSUSY) + 8640*z2*
        pow2(Mgl)*pow2(Msq)*pow2(Mst1)*pow2(Mst2)*pow2(MuSUSY) + 98884*Dmglst1*
        Mgl*pow2(Mst1)*pow2(Mst2)*pow4(Msq) + 597600*Dmglst1*Mgl*z2*pow2(Mst1)*
        pow2(Mst2)*pow4(Msq) + 228604*pow2(Dmglst1)*pow2(Mst1)*pow2(Mst2)*pow4(
        Msq) + 1193976*z2*pow2(Dmglst1)*pow2(Mst1)*pow2(Mst2)*pow4(Msq) + 996*
        pow2(Mgl)*pow2(Mst1)*pow2(Mst2)*pow4(Msq) + 200016*z2*pow2(Mgl)*pow2(
        Mst1)*pow2(Mst2)*pow4(Msq) - 83814*Dmglst1*Mgl*pow2(Mst1)*pow2(MuSUSY)*
        pow4(Msq) + 71280*B4*Dmglst1*Mgl*pow2(Mst1)*pow2(MuSUSY)*pow4(Msq) -
        648*Dmglst1*DN*Mgl*pow2(Mst1)*pow2(MuSUSY)*pow4(Msq) + 148392*Dmglst1*
        Mgl*z2*pow2(Mst1)*pow2(MuSUSY)*pow4(Msq) - 83814*pow2(Dmglst1)*pow2(
        Mst1)*pow2(MuSUSY)*pow4(Msq) + 71280*B4*pow2(Dmglst1)*pow2(Mst1)*pow2(
        MuSUSY)*pow4(Msq) - 648*DN*pow2(Dmglst1)*pow2(Mst1)*pow2(MuSUSY)*pow4(
        Msq) + 148392*z2*pow2(Dmglst1)*pow2(Mst1)*pow2(MuSUSY)*pow4(Msq) -
        214598*pow2(Mgl)*pow2(Mst1)*pow2(MuSUSY)*pow4(Msq) + 38448*B4*pow2(Mgl)
        *pow2(Mst1)*pow2(MuSUSY)*pow4(Msq) + 216*DN*pow2(Mgl)*pow2(Mst1)*pow2(
        MuSUSY)*pow4(Msq) - 173304*z2*pow2(Mgl)*pow2(Mst1)*pow2(MuSUSY)*pow4(
        Msq) - 23760*Dmglst1*Mgl*pow2(Msq)*pow2(Mst2)*pow4(Mst1) - 129600*
        Dmglst1*Mgl*z2*pow2(Msq)*pow2(Mst2)*pow4(Mst1) - 103680*pow2(Dmglst1)*
        pow2(Msq)*pow2(Mst2)*pow4(Mst1) - 388800*z2*pow2(Dmglst1)*pow2(Msq)*
        pow2(Mst2)*pow4(Mst1) - 2160*pow2(Mgl)*pow2(Msq)*pow2(Mst2)*pow4(Mst1)
        - 25920*z2*pow2(Mgl)*pow2(Msq)*pow2(Mst2)*pow4(Mst1) + 125280*Dmglst1*
        Mgl*pow2(Msq)*pow2(MuSUSY)*pow4(Mst1) - 51840*Dmglst1*Mgl*z2*pow2(Msq)*
        pow2(MuSUSY)*pow4(Mst1) + 146160*pow2(Dmglst1)*pow2(Msq)*pow2(MuSUSY)*
        pow4(Mst1) - 233280*z2*pow2(Dmglst1)*pow2(Msq)*pow2(MuSUSY)*pow4(Mst1)
        + 40320*pow2(Mgl)*pow2(Msq)*pow2(MuSUSY)*pow4(Mst1) + 73620*Dmglst1*
        Mgl*pow2(Mst2)*pow2(MuSUSY)*pow4(Mst1) + 16200*Dmglst1*Mgl*z2*pow2(
        Mst2)*pow2(MuSUSY)*pow4(Mst1) + 250560*pow2(Dmglst1)*pow2(Mst2)*pow2(
        MuSUSY)*pow4(Mst1) + 48600*z2*pow2(Dmglst1)*pow2(Mst2)*pow2(MuSUSY)*
        pow4(Mst1) + 12060*pow2(Mgl)*pow2(Mst2)*pow2(MuSUSY)*pow4(Mst1) + 3240*
        z2*pow2(Mgl)*pow2(Mst2)*pow2(MuSUSY)*pow4(Mst1) + 216148*Dmglst1*Mgl*
        pow4(Msq)*pow4(Mst1) + 1592568*Dmglst1*Mgl*z2*pow4(Msq)*pow4(Mst1) +
        1611936*pow2(Dmglst1)*pow4(Msq)*pow4(Mst1) + 4670928*z2*pow2(Dmglst1)*
        pow4(Msq)*pow4(Mst1) - 81172*pow2(Mgl)*pow4(Msq)*pow4(Mst1) + 330408*
        z2*pow2(Mgl)*pow4(Msq)*pow4(Mst1) + 36*pow2(log(pow2(Mst2)/pow2(Mst1)))
        *(pow2(Mgl)*(-120*pow2(Msq)*pow2(Mst1)*pow2(Mst2)*(3*pow2(Mst1) - pow2(
        MuSUSY)) + 45*pow2(Mst2)*pow2(MuSUSY)*pow4(Mst1) + pow4(Msq)*(pow2(
        Mst1)*(2094*pow2(Mst2) - 5641*pow2(MuSUSY)) + 4649*pow4(Mst1) - 96*
        pow4(Mst2))) + Dmglst1*Mgl*(225*pow2(Mst2)*pow2(MuSUSY)*pow4(Mst1) -
        360*pow2(Msq)*(-(pow2(Mst1)*pow2(Mst2)*pow2(MuSUSY)) + (5*pow2(Mst2) +
        2*pow2(MuSUSY))*pow4(Mst1)) + pow4(Msq)*(pow2(Mst1)*(6536*pow2(Mst2) -
        1749*pow2(MuSUSY)) + 23091*pow4(Mst1) - 96*pow4(Mst2))) + pow2(Dmglst1)
        *(675*pow2(Mst2)*pow2(MuSUSY)*pow4(Mst1) - 360*pow2(Msq)*(-(pow2(Mst1)*
        pow2(Mst2)*pow2(MuSUSY)) + 3*(5*pow2(Mst2) + 3*pow2(MuSUSY))*pow4(Mst1)
        ) + pow4(Msq)*(pow2(Mst1)*(13007*pow2(Mst2) - 1749*pow2(MuSUSY)) +
        66686*pow4(Mst1) - 96*pow4(Mst2)))) - 72*pow3(log(pow2(Mst2)/pow2(Mst1)
        ))*pow4(Msq)*(pow2(Mgl)*(pow2(Mst1)*(279*pow2(Mst2) - 214*pow2(MuSUSY))
        + 384*pow4(Mst1) - 14*pow4(Mst2)) + Dmglst1*Mgl*(5*pow2(Mst1)*(169*
        pow2(Mst2) + 84*pow2(MuSUSY)) + 1830*pow4(Mst1) - 14*pow4(Mst2)) +
        pow2(Dmglst1)*(14*pow2(Mst1)*(121*pow2(Mst2) + 30*pow2(MuSUSY)) + 5223*
        pow4(Mst1) - 14*pow4(Mst2))) + 92160*Dmglst1*Mgl*pow2(Msq)*pow2(Mst1)*
        pow4(Mst2) - 8640*Dmglst1*Mgl*z2*pow2(Msq)*pow2(Mst1)*pow4(Mst2) +
        238320*pow2(Dmglst1)*pow2(Msq)*pow2(Mst1)*pow4(Mst2) - 8640*z2*pow2(
        Dmglst1)*pow2(Msq)*pow2(Mst1)*pow4(Mst2) + 26400*pow2(Mgl)*pow2(Msq)*
        pow2(Mst1)*pow4(Mst2) - 8640*z2*pow2(Mgl)*pow2(Msq)*pow2(Mst1)*pow4(
        Mst2) - 3765*Dmglst1*Mgl*pow2(Mst1)*pow2(MuSUSY)*pow4(Mst2) - 3765*
        pow2(Dmglst1)*pow2(Mst1)*pow2(MuSUSY)*pow4(Mst2) - 1605*pow2(Mgl)*pow2(
        Mst1)*pow2(MuSUSY)*pow4(Mst2) + 84084*Dmglst1*Mgl*pow4(Msq)*pow4(Mst2)
        + 1872*Dmglst1*Mgl*z2*pow4(Msq)*pow4(Mst2) + 441984*pow2(Dmglst1)*pow4(
        Msq)*pow4(Mst2) + 1872*z2*pow2(Dmglst1)*pow4(Msq)*pow4(Mst2) + 13356*
        pow2(Mgl)*pow4(Msq)*pow4(Mst2) + 1872*z2*pow2(Mgl)*pow4(Msq)*pow4(Mst2)
        + 66240*Dmglst1*Mgl*pow4(Mst1)*pow4(Mst2) - 15120*Dmglst1*Mgl*z2*pow4(
        Mst1)*pow4(Mst2) + 199800*pow2(Dmglst1)*pow4(Mst1)*pow4(Mst2) - 28080*
        z2*pow2(Dmglst1)*pow4(Mst1)*pow4(Mst2) + 18000*pow2(Mgl)*pow4(Mst1)*
        pow4(Mst2) - 6480*z2*pow2(Mgl)*pow4(Mst1)*pow4(Mst2) - 8040*Dmglst1*
        Mgl*pow2(Msq)*pow6(Mst2) + 8640*Dmglst1*Mgl*z2*pow2(Msq)*pow6(Mst2) -
        8040*pow2(Dmglst1)*pow2(Msq)*pow6(Mst2) + 8640*z2*pow2(Dmglst1)*pow2(
        Msq)*pow6(Mst2) - 8040*pow2(Mgl)*pow2(Msq)*pow6(Mst2) + 8640*z2*pow2(
        Mgl)*pow2(Msq)*pow6(Mst2) - 20160*Dmglst1*Mgl*pow2(Mst1)*pow6(Mst2) +
        12960*Dmglst1*Mgl*z2*pow2(Mst1)*pow6(Mst2) - 40320*pow2(Dmglst1)*pow2(
        Mst1)*pow6(Mst2) + 25920*z2*pow2(Dmglst1)*pow2(Mst1)*pow6(Mst2) - 6720*
        pow2(Mgl)*pow2(Mst1)*pow6(Mst2) + 4320*z2*pow2(Mgl)*pow2(Mst1)*pow6(
        Mst2) - 6*log(pow2(Mst2)/pow2(Mst1))*(pow2(Mgl)*(2*pow4(Msq)*(pow2(
        Mst1)*(3*(2197 + 3396*z2)*pow2(Mst2) - 4*(10441 + 612*z2)*pow2(MuSUSY))
        + 6*(2903 + 2196*z2)*pow4(Mst1) + 36*(13 - 18*z2)*pow4(Mst2)) - 120*
        pow2(Msq)*(3*(5*pow2(Mst2) - 36*pow2(MuSUSY))*pow4(Mst1) - 8*pow2(Mst1)
        *(5*pow2(Mst2)*pow2(MuSUSY) + 2*pow4(Mst2)) + 4*pow6(Mst2)) + 15*(-2*
        pow2(Mst1)*(17*pow2(Mst2) + pow2(MuSUSY))*pow4(Mst2) + 3*pow4(Mst1)*(
        27*pow2(Mst2)*pow2(MuSUSY) + 7*pow4(Mst2)) - 29*pow8(Mst2))) + pow2(
        Dmglst1)*(2*pow4(Msq)*(4*pow2(Mst1)*((13499 + 14742*z2)*pow2(Mst2) + 3*
        (-2047 + 1290*z2)*pow2(MuSUSY)) + 13*(20965 + 14292*z2)*pow4(Mst1) - (
        1067 + 648*z2)*pow4(Mst2)) - 120*pow2(Msq)*((423*pow2(Mst2) - 366*pow2(
        MuSUSY))*pow4(Mst1) - 12*pow2(Mst1)*(11*pow2(Mst2)*pow2(MuSUSY) + 8*
        pow4(Mst2)) + 4*pow6(Mst2)) + 15*(-2*pow2(Mst1)*(102*pow2(Mst2) + pow2(
        MuSUSY))*pow4(Mst2) + 15*pow4(Mst1)*(101*pow2(Mst2)*pow2(MuSUSY) + 21*
        pow4(Mst2)) - 29*pow8(Mst2))) + Dmglst1*Mgl*(2*pow4(Msq)*(3*pow2(Mst1)*
        (9*(947 + 1100*z2)*pow2(Mst2) + 4*(-2047 + 1290*z2)*pow2(MuSUSY)) + 16*
        (5696 + 3987*z2)*pow4(Mst1) - 2*(145 + 324*z2)*pow4(Mst2)) - 120*pow2(
        Msq)*(3*(37*pow2(Mst2) - 96*pow2(MuSUSY))*pow4(Mst1) - 12*pow2(Mst1)*(
        11*pow2(Mst2)*pow2(MuSUSY) + 4*pow4(Mst2)) + 4*pow6(Mst2)) + 15*(-2*
        pow2(Mst1)*(51*pow2(Mst2) + pow2(MuSUSY))*pow4(Mst2) + 3*pow4(Mst1)*(
        159*pow2(Mst2)*pow2(MuSUSY) + 35*pow4(Mst2)) - 29*pow8(Mst2)))) - 5070*
        Dmglst1*Mgl*pow8(Mst2) + 2160*Dmglst1*Mgl*z2*pow8(Mst2) - 5070*pow2(
        Dmglst1)*pow8(Mst2) + 2160*z2*pow2(Dmglst1)*pow8(Mst2) - 5070*pow2(Mgl)
        *pow8(Mst2) + 2160*z2*pow2(Mgl)*pow8(Mst2) + 1080*pow2(log(pow2(Msq)/
        pow2(Mst1)))*(6*log(pow2(Mst2)/pow2(Mst1))*(Dmglst1*Mgl + pow2(Dmglst1)
        + pow2(Mgl))*pow4(Msq)*(3*pow2(Mst1)*pow2(MuSUSY) + pow4(Mst2)) + pow2(
        Mgl)*(-12*pow2(Mst1)*pow2(MuSUSY)*pow4(Msq) - 3*pow4(Mst1)*pow4(Mst2) +
        2*pow2(Mst1)*pow6(Mst2) + pow2(Msq)*(-4*pow2(Mst1)*pow4(Mst2) + 4*pow6(
        Mst2)) + pow8(Mst2)) + Dmglst1*Mgl*(-12*pow2(Mst1)*pow2(MuSUSY)*pow4(
        Msq) - 7*pow4(Mst1)*pow4(Mst2) + 6*pow2(Mst1)*pow6(Mst2) + pow2(Msq)*(-
        4*pow2(Mst1)*pow4(Mst2) + 4*pow6(Mst2)) + pow8(Mst2)) + pow2(Dmglst1)*(
        -12*pow2(Mst1)*pow2(MuSUSY)*pow4(Msq) - 13*pow4(Mst1)*pow4(Mst2) + 12*
        pow2(Mst1)*pow6(Mst2) + pow2(Msq)*(-4*pow2(Mst1)*pow4(Mst2) + 4*pow6(
        Mst2)) + pow8(Mst2))) + 180*log(pow2(Msq)/pow2(Mst1))*(36*pow2(Mst1)*(
        pow2(Mgl)*(6*pow2(Mst1) + 4*pow2(Mst2) + pow2(MuSUSY)) + 3*Dmglst1*Mgl*
        (10*pow2(Mst1) + 4*pow2(Mst2) + 3*pow2(MuSUSY)) + 3*pow2(Dmglst1)*(30*
        pow2(Mst1) + 8*pow2(Mst2) + 3*pow2(MuSUSY)))*pow2(log(pow2(Mst2)/pow2(
        Mst1)))*pow4(Msq) + pow2(Mgl)*(18*pow4(Msq)*(2*pow2(Mst1)*((1 + 8*z2)*
        pow2(Mst2) + 2*(-5 + z2)*pow2(MuSUSY)) + 2*(-7 + 12*z2)*pow4(Mst1) -
        11*pow4(Mst2)) - pow2(Mst1)*(4*pow2(Mst2) + pow2(MuSUSY))*pow4(Mst2) +
        12*pow4(Mst1)*(3*pow2(Mst2)*pow2(MuSUSY) + pow4(Mst2)) + 4*pow2(Msq)*(
        24*pow2(MuSUSY)*pow4(Mst1) + pow2(Mst1)*(24*pow2(Mst2)*pow2(MuSUSY) -
        7*pow4(Mst2)) + 7*pow6(Mst2)) - 8*pow8(Mst2)) + Dmglst1*Mgl*(2*pow4(
        Msq)*(18*pow2(Mst1)*((7 + 24*z2)*pow2(Mst2) + 2*(-2 + 9*z2)*pow2(
        MuSUSY)) + 90*(-1 + 12*z2)*pow4(Mst1) - 251*pow4(Mst2)) - pow2(Mst1)*(
        12*pow2(Mst2) + pow2(MuSUSY))*pow4(Mst2) + 20*pow4(Mst1)*(9*pow2(Mst2)*
        pow2(MuSUSY) + pow4(Mst2)) + 4*pow2(Msq)*(72*pow2(MuSUSY)*pow4(Mst1) +
        pow2(Mst1)*(72*pow2(Mst2)*pow2(MuSUSY) - 7*pow4(Mst2)) + 7*pow6(Mst2))
        - 8*pow8(Mst2)) + pow2(Dmglst1)*(2*pow4(Msq)*(36*pow2(Mst1)*(2*(5 + 12*
        z2)*pow2(Mst2) + (-2 + 9*z2)*pow2(MuSUSY)) + 1080*(1 + 3*z2)*pow4(Mst1)
        - 419*pow4(Mst2)) - pow2(Mst1)*(24*pow2(Mst2) + pow2(MuSUSY))*pow4(
        Mst2) + 4*pow4(Mst1)*(135*pow2(Mst2)*pow2(MuSUSY) + 8*pow4(Mst2)) + 4*
        pow2(Msq)*(144*pow2(MuSUSY)*pow4(Mst1) + pow2(Mst1)*(72*pow2(Mst2)*
        pow2(MuSUSY) - 7*pow4(Mst2)) + 7*pow6(Mst2)) - 8*pow8(Mst2)) - 6*log(
        pow2(Mst2)/pow2(Mst1))*(pow2(Mgl)*(3*pow2(Mst2)*(pow2(Mst2) + pow2(
        MuSUSY))*pow4(Mst1) + 6*pow4(Msq)*(pow2(Mst1)*(4*pow2(Mst2) - 17*pow2(
        MuSUSY)) + 17*pow4(Mst1) - 2*pow4(Mst2)) + 2*pow2(Mst1)*pow6(Mst2) + 4*
        pow2(Msq)*(2*pow2(Mst1)*pow2(Mst2)*(pow2(Mst2) + pow2(MuSUSY)) + 6*
        pow2(MuSUSY)*pow4(Mst1) + pow6(Mst2)) + pow8(Mst2)) + Dmglst1*Mgl*(15*
        pow2(Mst2)*(pow2(Mst2) + pow2(MuSUSY))*pow4(Mst1) + 6*pow4(Msq)*(pow2(
        Mst1)*(20*pow2(Mst2) - pow2(MuSUSY)) + 97*pow4(Mst1) - 2*pow4(Mst2)) +
        6*pow2(Mst1)*pow6(Mst2) + 4*pow2(Msq)*(6*pow2(Mst1)*pow2(Mst2)*(pow2(
        Mst2) + pow2(MuSUSY)) + 18*pow2(MuSUSY)*pow4(Mst1) + pow6(Mst2)) +
        pow8(Mst2)) + pow2(Dmglst1)*(45*pow2(Mst2)*(pow2(Mst2) + pow2(MuSUSY))*
        pow4(Mst1) + 6*pow4(Msq)*(pow2(Mst1)*(52*pow2(Mst2) - pow2(MuSUSY)) +
        321*pow4(Mst1) - 2*pow4(Mst2)) + 12*pow2(Mst1)*pow6(Mst2) + 4*pow2(Msq)
        *(6*pow2(Mst1)*pow2(Mst2)*(2*pow2(Mst2) + pow2(MuSUSY)) + 36*pow2(
        MuSUSY)*pow4(Mst1) + pow6(Mst2)) + pow8(Mst2))))))/(243.*pow2(Mgl)*
        pow4(Msq)*pow4(Mst2)) + (pow4(Mt)*(10501 + 4920*log(pow2(Msq)/pow2(
        Mst1)) + (27*(40*(1 + 2*log(pow2(Msq)/pow2(Mst1)))*pow2(Msq) + 4*(1 +
        2*log(pow2(Mst2)/pow2(Mst1)))*pow2(Mst2)))/pow2(Mst1) + 1980*pow2(log(
        pow2(Msq)/pow2(Mst1))) + 18*log(pow2(Mst2)/pow2(Mst1))*(-149 + 50*log(
        pow2(Msq)/pow2(Mst1)) + log(pow2(Mst2)/pow2(Mst1)) + 30*pow2(log(pow2(
        Msq)/pow2(Mst1)))) + (2*Dmglst1*(-2215 + 300*log(pow2(Msq)/pow2(Mst1))
        + 254*log(pow2(Mst2)/pow2(Mst1)) + 1620*pow2(log(pow2(Msq)/pow2(Mst1)))
        + 162*pow2(log(pow2(Mst2)/pow2(Mst1)))))/Mgl - 360*pow3(log(pow2(Msq)/
        pow2(Mst1))) + 204*pow3(log(pow2(Mst2)/pow2(Mst1))) + (12*(3*z2*pow2(
        Dmglst1)*pow2(Mst1)*(360*pow2(Msq)*pow2(Mst1)*pow4(Mst2) + 590*pow4(
        Mst1)*pow4(Mst2) + pow4(Msq)*(4*(151 - 108*log(pow2(Mst2)/pow2(Mst1)))*
        pow2(Mst1)*pow2(Mst2) - 4*(2517 + 750*log(pow2(Msq)/pow2(Mst1)) - 1198*
        log(pow2(Mst2)/pow2(Mst1)))*pow4(Mst1) + 297*pow4(Mst2)) + 90*pow2(
        Mst1)*pow6(Mst2)) + 2*Dmglst1*Mgl*z2*pow2(Mst1)*(360*pow2(Msq)*pow2(
        Mst1)*pow4(Mst2) + 390*pow4(Mst1)*pow4(Mst2) + pow4(Msq)*(4*(151 - 108*
        log(pow2(Mst2)/pow2(Mst1)))*pow2(Mst1)*pow2(Mst2) - 18*(311 + 100*log(
        pow2(Msq)/pow2(Mst1)) - 154*log(pow2(Mst2)/pow2(Mst1)))*pow4(Mst1) +
        297*pow4(Mst2)) + 90*pow2(Mst1)*pow6(Mst2)) + pow2(Mgl)*(120*z2*pow2(
        Msq)*pow2(Mst1)*(4*pow2(Mst1) + pow2(Mst2))*pow4(Mst2) + 15*z2*pow2(
        Mst1)*pow4(Mst2)*(6*pow2(Mst1)*pow2(Mst2) + 19*pow4(Mst1) + 3*pow4(
        Mst2)) - 180*z2*pow2(Mst2)*(pow2(Mst1) + pow2(Mst2))*pow6(Msq) + pow4(
        Msq)*(2*(z2*(289 - 228*log(pow2(Mst2)/pow2(Mst1)))*pow2(Mst2) + 2*pow2(
        MuSUSY)*(3*(-139 + 2*B4 + DN - 110*z2) + 6*(87 + 16*z2)*log(pow2(Mst2)/
        pow2(Mst1)) - 144*pow2(log(pow2(Mst2)/pow2(Mst1))) + 2*pow3(log(pow2(
        Mst2)/pow2(Mst1)))))*pow4(Mst1) - 8*z2*(-49 + 45*log(pow2(Msq)/pow2(
        Mst1)))*pow2(Mst1)*pow4(Mst2) + z2*(-2383 - 900*log(pow2(Msq)/pow2(
        Mst1)) + 1230*log(pow2(Mst2)/pow2(Mst1)))*pow6(Mst1) - 18*z2*pow6(Mst2)
        ))))/(pow2(Mgl)*pow2(Mst1)*pow4(Msq)*pow4(Mst2)) - (43967798000*
        Dmglst1*Mgl*pow2(Mst1)*pow2(Mst2)*pow4(Msq) + 121598741040*pow2(
        Dmglst1)*pow2(Mst1)*pow2(Mst2)*pow4(Msq) + 13498218944*pow2(Mgl)*pow2(
        Mst1)*pow2(Mst2)*pow4(Msq) - 6174000*pow2(Mst1)*(2*pow2(Dmglst1)*(898*
        pow2(Mst1) - 99*pow2(Mst2)) + 5*pow2(Mgl)*(29*pow2(Mst1) - 14*pow2(
        Mst2)) + 12*Dmglst1*Mgl*(57*pow2(Mst1) - 11*pow2(Mst2)))*pow3(log(pow2(
        Mst2)/pow2(Mst1)))*pow4(Msq) + 6791400000*Dmglst1*Mgl*pow2(Msq)*pow2(
        Mst2)*pow4(Mst1) + 16978500000*pow2(Dmglst1)*pow2(Msq)*pow2(Mst2)*pow4(
        Mst1) + 1697850000*pow2(Mgl)*pow2(Msq)*pow2(Mst2)*pow4(Mst1) +
        102884565000*Dmglst1*Mgl*pow4(Msq)*pow4(Mst1) + 258032664260*pow2(
        Dmglst1)*pow4(Msq)*pow4(Mst1) + 25346116504*pow2(Mgl)*pow4(Msq)*pow4(
        Mst1) + 5000940000*Dmglst1*Mgl*pow2(Msq)*pow2(Mst1)*pow4(Mst2) +
        8270690400*pow2(Dmglst1)*pow2(Msq)*pow2(Mst1)*pow4(Mst2) + 3047239440*
        pow2(Mgl)*pow2(Msq)*pow2(Mst1)*pow4(Mst2) - 4992461040*pow2(Dmglst1)*
        pow4(Msq)*pow4(Mst2) + 5389387500*Dmglst1*Mgl*pow4(Mst1)*pow4(Mst2) +
        11080426350*pow2(Dmglst1)*pow4(Mst1)*pow4(Mst2) + 2412641055*pow2(Mgl)*
        pow4(Mst1)*pow4(Mst2) - 555660000*pow2(Mgl)*pow2(Mst2)*pow6(Msq) +
        1112719440*pow2(Mgl)*pow2(Msq)*pow6(Mst2) + 1775025000*Dmglst1*Mgl*
        pow2(Mst1)*pow6(Mst2) + 2662537500*pow2(Dmglst1)*pow2(Mst1)*pow6(Mst2)
        + 887512500*pow2(Mgl)*pow2(Mst1)*pow6(Mst2) - 1323000*pow2(Mst2)*pow2(
        log(pow2(Msq)/pow2(Mst1)))*(140*Dmglst1*Mgl*pow2(Mst1)*(12*pow2(Msq)*
        pow2(Mst2) + 13*pow2(Mst1)*pow2(Mst2) + 9*pow4(Msq) + 3*pow4(Mst2)) +
        70*pow2(Dmglst1)*(36*pow2(Msq)*pow2(Mst1)*pow2(Mst2) + 27*(pow2(Mst1) +
        pow2(Mst2))*pow4(Msq) + 59*pow2(Mst2)*pow4(Mst1) + 9*pow2(Mst1)*pow4(
        Mst2)) + pow2(Mgl)*(630*pow2(Mst1)*pow4(Msq) + 617*pow2(Mst2)*pow4(
        Mst1) + 210*pow2(Mst1)*pow4(Mst2) + 28*pow2(Msq)*(37*pow2(Mst1)*pow2(
        Mst2) + 7*pow4(Mst2)) + 57*pow6(Mst2))) + 611891055*pow2(Mgl)*pow8(
        Mst2) + 88200*pow2(log(pow2(Mst2)/pow2(Mst1)))*(140*Dmglst1*Mgl*pow2(
        Mst1)*(1449*pow2(Mst1) - 353*pow2(Mst2))*pow4(Msq) + 105*pow2(Dmglst1)*
        pow4(Msq)*(-706*pow2(Mst1)*pow2(Mst2) + 4944*pow4(Mst1) - 27*pow4(Mst2)
        ) + pow2(Mgl)*(pow4(Msq)*(-23044*pow2(Mst1)*pow2(Mst2) + 44291*pow4(
        Mst1)) + 1260*pow2(Msq)*pow6(Mst2) + 720*pow8(Mst2))) + 6300*log(pow2(
        Msq)/pow2(Mst1))*(441000*(4*Dmglst1*Mgl + 10*pow2(Dmglst1) + pow2(Mgl))
        *pow2(log(pow2(Mst2)/pow2(Mst1)))*pow4(Msq)*pow4(Mst1) + 4900*Dmglst1*
        Mgl*pow2(Mst1)*(162*(2*pow2(Mst1) - pow2(Mst2))*pow4(Msq) + 137*pow2(
        Mst1)*pow4(Mst2) + 24*pow2(Msq)*(5*pow2(Mst1)*pow2(Mst2) + 3*pow4(Mst2)
        ) + 24*pow6(Mst2)) + 98*pow2(Dmglst1)*(2*pow4(Msq)*(-6075*pow2(Mst1)*
        pow2(Mst2) + 29250*pow4(Mst1) + 2006*pow4(Mst2)) + 600*pow2(Msq)*(25*
        pow2(Mst2)*pow4(Mst1) + 19*pow2(Mst1)*pow4(Mst2)) + 25*(763*pow4(Mst1)*
        pow4(Mst2) + 72*pow2(Mst1)*pow6(Mst2))) - 210*log(pow2(Mst2)/pow2(Mst1)
        )*(70*pow2(Dmglst1)*pow2(Mst1)*(990*pow2(Mst1)*pow4(Msq) - 18*pow2(Msq)
        *pow4(Mst2) - 25*pow2(Mst1)*pow4(Mst2) - 9*pow6(Mst2)) + 140*Dmglst1*
        Mgl*pow2(Mst1)*(174*pow2(Mst1)*pow4(Msq) - 6*pow2(Msq)*pow4(Mst2) - 5*
        pow2(Mst1)*pow4(Mst2) - 3*pow6(Mst2)) + pow2(Mgl)*(5040*pow4(Msq)*pow4(
        Mst1) - pow4(Mst2)*(210*pow2(Mst1)*pow2(Mst2) + 175*pow4(Mst1) + 9*
        pow4(Mst2)) - 28*pow2(Msq)*(15*pow2(Mst1)*pow4(Mst2) + 4*pow6(Mst2))))
        + pow2(Mgl)*(44100*pow4(Msq)*(-9*pow2(Mst1)*pow2(Mst2) + 4*pow4(Mst1))
        + 220709*pow4(Mst1)*pow4(Mst2) - 176400*pow2(Mst2)*pow6(Msq) + 58800*
        pow2(Mst1)*pow6(Mst2) + 392*pow2(Msq)*(375*pow2(Mst2)*pow4(Mst1) + 116*
        pow2(Mst1)*pow4(Mst2) + 41*pow6(Mst2)) + 49209*pow8(Mst2))) - 210*log(
        pow2(Mst2)/pow2(Mst1))*(4900*Dmglst1*Mgl*pow2(Mst1)*(4*(8608*pow2(Mst1)
        - 2069*pow2(Mst2))*pow4(Msq) + 180*pow2(Msq)*pow4(Mst2) + 15*(55*pow2(
        Mst1) + 69*pow2(Mst2))*pow4(Mst2)) + 98*pow2(Dmglst1)*(13500*pow2(Msq)*
        pow2(Mst1)*pow4(Mst2) + pow4(Msq)*(-634300*pow2(Mst1)*pow2(Mst2) +
        4044750*pow4(Mst1) + 6364*pow4(Mst2)) + 375*(275*pow4(Mst1)*pow4(Mst2)
        + 207*pow2(Mst1)*pow6(Mst2))) + pow2(Mgl)*(4*pow4(Msq)*(-5313854*pow2(
        Mst1)*pow2(Mst2) + 9874311*pow4(Mst1)) + 5880*pow2(Msq)*(75*pow2(Mst1)*
        pow4(Mst2) + 482*pow6(Mst2)) + 15*(67375*pow4(Mst1)*pow4(Mst2) +
        169050*pow2(Mst1)*pow6(Mst2) + 131493*pow8(Mst2)))))/(514500.*pow2(Mgl)
        *pow4(Msq)*pow4(Mst2))))/81. + pow2(s2t)*(pow2(Mt)*((40*(1 + 2*log(
        pow2(Msq)/pow2(Mst1)))*pow2(Msq))/3. - (pow2(Mst1)*(16282638750*
        Dmglst1*Mgl*pow2(Msq)*pow2(Mst1) + 17400518625*pow2(Dmglst1)*pow2(Msq)*
        pow2(Mst1) + 14893594207*pow2(Mgl)*pow2(Msq)*pow2(Mst1) + 630*log(pow2(
        Mst2)/pow2(Mst1))*pow2(Mst1)*(3*pow2(Mgl)*(14193447*pow2(Msq) - 245000*
        pow2(Mst2)) + 2450*Dmglst1*Mgl*(37061*pow2(Msq) - 1260*pow2(Mst2)) +
        50225*pow2(Dmglst1)*(1871*pow2(Msq) - 180*pow2(Mst2))) + 6065955000*
        Dmglst1*Mgl*pow2(Mst1)*pow2(Mst2) + 13358992500*pow2(Dmglst1)*pow2(
        Mst1)*pow2(Mst2) + 1397917080*pow2(Mgl)*pow2(Mst1)*pow2(Mst2) +
        13891500*pow2(Mst1)*(90*Dmglst1*Mgl*pow2(Msq) + 135*pow2(Dmglst1)*pow2(
        Msq) + pow2(Mgl)*(45*pow2(Msq) - 4*pow2(Mst2)))*pow2(log(pow2(Msq)/
        pow2(Mst1))) + 66150*(379120*Dmglst1*Mgl + 1251880*pow2(Dmglst1) +
        46219*pow2(Mgl))*pow2(Msq)*pow2(Mst1)*pow2(log(pow2(Mst2)/pow2(Mst1)))
        - 9261000*(908*Dmglst1*Mgl + 2468*pow2(Dmglst1) + 223*pow2(Mgl))*pow2(
        Msq)*pow2(Mst1)*pow3(log(pow2(Mst2)/pow2(Mst1))) + 416745000*pow2(Mgl)*
        pow4(Msq) + 926100*log(pow2(Msq)/pow2(Mst1))*(75*log(pow2(Mst2)/pow2(
        Mst1))*pow2(Mst1)*(84*Dmglst1*Mgl*pow2(Msq) + 246*pow2(Dmglst1)*pow2(
        Msq) + pow2(Mgl)*(12*pow2(Msq) - pow2(Mst2))) - 150*Dmglst1*Mgl*pow2(
        Mst1)*(27*pow2(Msq) - 40*pow2(Mst2)) + 75*pow2(Dmglst1)*pow2(Mst1)*(-
        309*pow2(Msq) + 200*pow2(Mst2)) + pow2(Mgl)*(1125*pow2(Msq)*pow2(Mst1)
        + 1282*pow2(Mst1)*pow2(Mst2) + 900*pow4(Msq)))))/(6.251175e7*pow2(Mgl)*
        pow2(Msq)*pow2(Mst2)) - (pow2(Mst2)*(4501 - 288*DN + 1080*log(pow2(Msq)
        /pow2(Mst1)) + (2160*(1 + 2*log(pow2(Msq)/pow2(Mst1)))*pow2(Msq))/pow2(
        Mst1) + 3240*pow2(log(pow2(Msq)/pow2(Mst1))) - 24*log(pow2(Mst2)/pow2(
        Mst1))*(733 + 210*log(pow2(Msq)/pow2(Mst1)) + 45*pow2(log(pow2(Msq)/
        pow2(Mst1)))) + 72*(191 + 30*log(pow2(Msq)/pow2(Mst1)))*pow2(log(pow2(
        Mst2)/pow2(Mst1))) + (24*Dmglst1*(-158 - 240*log(pow2(Msq)/pow2(Mst1))
        + 96*log(pow2(Mst2)/pow2(Mst1)) + 270*pow2(log(pow2(Msq)/pow2(Mst1))) +
        27*pow2(log(pow2(Mst2)/pow2(Mst1)))))/Mgl - 4728*pow3(log(pow2(Mst2)/
        pow2(Mst1))) + (324*((pow2(Dmglst1)*(1467 - 520*log(pow2(Msq)/pow2(
        Mst1)) + 148*log(pow2(Mst2)/pow2(Mst1)) + 270*pow2(log(pow2(Msq)/pow2(
        Mst1))) + 27*pow2(log(pow2(Mst2)/pow2(Mst1)))))/9. - (2*pow2(Mst1)*(
        15750*Dmglst1*Mgl + 10125*pow2(Dmglst1) - 3091*pow2(Mgl) + 30*log(pow2(
        Mst2)/pow2(Mst1))*(450*Dmglst1*Mgl + 675*pow2(Dmglst1) + 334*pow2(Mgl))
        + 30*log(pow2(Msq)/pow2(Mst1))*(1500*Dmglst1*Mgl + 3750*pow2(Dmglst1) +
        157*pow2(Mgl) + 60*log(pow2(Mst2)/pow2(Mst1))*pow2(Mgl)) - 1800*pow2(
        Mgl)*pow2(log(pow2(Msq)/pow2(Mst1))) + 1350*pow2(Mgl)*pow2(log(pow2(
        Mst2)/pow2(Mst1)))))/(2025.*pow2(Msq)) - ((291550*Dmglst1*Mgl -
        4261775*pow2(Dmglst1) + 522392*pow2(Mgl) + 4116*log(pow2(Mst2)/pow2(
        Mst1))*(350*Dmglst1*Mgl + 1025*pow2(Dmglst1) + 113*pow2(Mgl)) + 420*
        log(pow2(Msq)/pow2(Mst1))*(13230*Dmglst1*Mgl + 40425*pow2(Dmglst1) +
        2929*pow2(Mgl) + 441*log(pow2(Mst2)/pow2(Mst1))*pow2(Mgl)) - 44100*
        pow2(Mgl)*pow2(log(pow2(Msq)/pow2(Mst1))))*pow4(Mst1))/(185220.*pow4(
        Msq))))/pow2(Mgl)))/324. + (2*(-675*(1 + 2*log(pow2(Mst2)/pow2(Mst1)))*
        pow2(Msq) + pow2(Mst1)*(5108 + log(pow2(Msq)/pow2(Mst1))*(3270 - 2925*
        log(pow2(Mst2)/pow2(Mst1))) - 6270*log(pow2(Mst2)/pow2(Mst1)) + 900*
        pow2(log(pow2(Msq)/pow2(Mst1))) + 2025*pow2(log(pow2(Mst2)/pow2(Mst1)))
        ))*pow4(Mst2))/(2025.*pow2(Msq)*pow2(Mst1)) + pow2(Mst1)*(
        268.4533251028807 - (8*DN)/9. + (100*log(pow2(Msq)/pow2(Mst1)))/3. +
        20*pow2(log(pow2(Msq)/pow2(Mst1))) - (2*log(pow2(Mst2)/pow2(Mst1))*(
        33703 + 4000*log(pow2(Msq)/pow2(Mst1)) + 375*pow2(log(pow2(Msq)/pow2(
        Mst1)))))/225. - (2*pow2(MuSUSY)*(836 + 353*log(pow2(Mst2)/pow2(Mst1))
        + 30*log(pow2(Msq)/pow2(Mst1))*(-7 + 4*log(pow2(Mst2)/pow2(Mst1))) -
        435*pow2(log(pow2(Mst2)/pow2(Mst1)))))/(135.*pow2(Msq)) + (4601*pow2(
        log(pow2(Mst2)/pow2(Mst1))))/45. + (80*pow3(log(pow2(Mst2)/pow2(Mst1)))
        )/9. + (Dmglst1*(156781 - 864*B4 - 432*DN - 131946*log(pow2(Mst2)/pow2(
        Mst1)) - 540*log(pow2(Msq)/pow2(Mst1))*(-19 + 18*log(pow2(Mst2)/pow2(
        Mst1))) + 9720*pow2(log(pow2(Msq)/pow2(Mst1))) + 36900*pow2(log(pow2(
        Mst2)/pow2(Mst1))) + 1440*pow3(log(pow2(Mst2)/pow2(Mst1)))))/(243.*Mgl)
        + (pow2(Dmglst1)*(log(pow2(Msq)/pow2(Mst1))*(32.22222222222222 - 60*
        log(pow2(Mst2)/pow2(Mst1))) + 60*pow2(log(pow2(Msq)/pow2(Mst1))) + (
        266599 - 864*B4 - 432*DN - 124698*log(pow2(Mst2)/pow2(Mst1)) + 34308*
        pow2(log(pow2(Mst2)/pow2(Mst1))) + 1440*pow3(log(pow2(Mst2)/pow2(Mst1))
        ))/162.))/pow2(Mgl) - ((668850*Dmglst1*Mgl + 1003275*pow2(Dmglst1) +
        438008*pow2(Mgl) - 84*log(pow2(Mst2)/pow2(Mst1))*(7350*Dmglst1*Mgl +
        11025*pow2(Dmglst1) + 1868*pow2(Mgl)) + 420*log(pow2(Msq)/pow2(Mst1))*(
        1470*Dmglst1*Mgl + 2205*pow2(Dmglst1) + 991*pow2(Mgl) + 231*log(pow2(
        Mst2)/pow2(Mst1))*pow2(Mgl)) + 44100*pow2(Mgl)*pow2(log(pow2(Msq)/pow2(
        Mst1))) - 141120*pow2(Mgl)*pow2(log(pow2(Mst2)/pow2(Mst1))))*pow4(Mst2)
        )/(185220.*pow2(Mgl)*pow4(Msq))) + (-3*pow2(Mst1)*pow2(Mst2)*((40*OepS2
        - 4941*S2)*pow2(Mst2) + 8*(136*OepS2 - 3969*S2)*pow2(MuSUSY)) + 2*((8*
        OepS2 - 999*S2)*pow2(Mst2) + 2*(-1480*OepS2 + 55863*S2)*pow2(MuSUSY))*
        pow4(Mst1) + 27*(-8*OepS2 + 81*S2)*pow6(Mst2) + 162*S2*log(pow2(Mst2)/
        pow2(Mst1))*(-2*(pow2(Mst2) - 370*pow2(MuSUSY))*pow4(Mst1) + 3*pow2(
        Mst1)*(136*pow2(Mst2)*pow2(MuSUSY) + 5*pow4(Mst2)) + 27*pow6(Mst2)))/(
        729.*pow4(Mst2)) + (pow2(MuSUSY)*(pow2(Mst1)*((-514500*Dmglst1*(-437365
        + 10368*B4 + 8640*DN)*Mgl + 257250*(-130639 + 24192*B4 + 1728*DN)*pow2(
        Dmglst1) + (531403689547 - 16892064000*B4 + 444528000*DN)*pow2(Mgl))*
        pow2(Mst1) + 1260*log(pow2(Mst2)/pow2(Mst1))*(335899900*Dmglst1*Mgl*
        pow2(Mst1) + 851392150*pow2(Dmglst1)*pow2(Mst1) + 9*pow2(Mgl)*(1176000*
        pow2(Msq) + 72625277*pow2(Mst1))) + 6667920000*log(pow2(Mst2)/pow2(
        Mst1))*(6*Dmglst1*Mgl + 9*pow2(Dmglst1) + 4*pow2(Mgl))*pow2(Mst1)*pow2(
        log(pow2(Msq)/pow2(Mst1))) - 264600*(1328740*Dmglst1*Mgl + 3207610*
        pow2(Dmglst1) + 2258059*pow2(Mgl))*pow2(Mst1)*pow2(log(pow2(Mst2)/pow2(
        Mst1))) - 555660000*log(pow2(Msq)/pow2(Mst1))*(3*(152*pow2(Dmglst1) -
        15*pow2(Mgl))*pow2(Mst1) - 2*log(pow2(Mst2)/pow2(Mst1))*(216*Dmglst1*
        Mgl*pow2(Mst1) + 504*pow2(Dmglst1)*pow2(Mst1) + pow2(Mgl)*(24*pow2(Msq)
        + 121*pow2(Mst1))) + 36*(2*Dmglst1*Mgl + 3*pow2(Dmglst1) + 4*pow2(Mgl))
        *pow2(Mst1)*pow2(log(pow2(Mst2)/pow2(Mst1)))) + 74088000*(2418*Dmglst1*
        Mgl + 4291*pow2(Dmglst1) + 3920*pow2(Mgl))*pow2(Mst1)*pow3(log(pow2(
        Mst2)/pow2(Mst1)))) + (12*pow2(Mst2)*(4430156444*pow2(Mgl)*pow2(Mst1)*
        pow4(Msq) - 1407672000*B4*pow2(Mgl)*pow2(Mst1)*pow4(Msq) + 37044000*DN*
        pow2(Mgl)*pow2(Mst1)*pow4(Msq) + 10718064000*pow2(Mgl)*pow2(Mst1)*pow3(
        log(pow2(Mst2)/pow2(Mst1)))*pow4(Msq) + 2646000*pow2(Mgl)*pow2(log(
        pow2(Msq)/pow2(Mst1)))*(pow2(-(Mst2*pow2(Mst1)) + pow3(Mst2)) - 315*
        pow2(Mst1)*pow4(Msq) + 840*log(pow2(Mst2)/pow2(Mst1))*pow2(Mst1)*pow4(
        Msq)) + 6791400000*Dmglst1*Mgl*pow2(Msq)*pow4(Mst1) + 20682900000*pow2(
        Dmglst1)*pow2(Msq)*pow4(Mst1) + 1883070000*pow2(Mgl)*pow2(Msq)*pow4(
        Mst1) - 3742987500*Dmglst1*Mgl*pow2(Mst2)*pow4(Mst1) - 11911961250*
        pow2(Dmglst1)*pow2(Mst2)*pow4(Mst1) - 618314205*pow2(Mgl)*pow2(Mst2)*
        pow4(Mst1) + 291716610*pow2(Mgl)*pow2(Mst1)*pow4(Mst2) - 6300*log(pow2(
        Msq)/pow2(Mst1))*(705600*pow2(Mgl)*pow2(Mst1)*pow2(log(pow2(Mst2)/pow2(
        Mst1)))*pow4(Msq) - 14700*Dmglst1*Mgl*(40*pow2(Msq) + pow2(Mst2))*pow4(
        Mst1) - 7350*pow2(Dmglst1)*(200*pow2(Msq) + 47*pow2(Mst2))*pow4(Mst1) -
        210*log(pow2(Mst2)/pow2(Mst1))*(-140*Dmglst1*Mgl*(28*pow2(Msq) - 5*
        pow2(Mst2))*pow4(Mst1) + 70*pow2(Dmglst1)*(-164*pow2(Msq) + 25*pow2(
        Mst2))*pow4(Mst1) + pow2(Mgl)*(6580*pow2(Mst1)*pow4(Msq) - 1232*pow2(
        Msq)*pow4(Mst1) - 27*pow2(Mst2)*pow4(Mst1) - 104*pow2(Mst1)*pow4(Mst2)
        - 18*pow6(Mst2))) + pow2(Mgl)*(308700*pow2(Mst1)*pow4(Msq) - 147000*
        pow2(Msq)*pow4(Mst1) + 12479*pow2(Mst2)*pow4(Mst1) - 26218*pow2(Mst1)*
        pow4(Mst2) - 9571*pow6(Mst2))) + 110989545*pow2(Mgl)*pow6(Mst2) -
        132300*pow2(log(pow2(Mst2)/pow2(Mst1)))*(-700*Dmglst1*Mgl*(28*pow2(Msq)
        - 5*pow2(Mst2))*pow4(Mst1) + 350*pow2(Dmglst1)*(-164*pow2(Msq) + 25*
        pow2(Mst2))*pow4(Mst1) + pow2(Mgl)*(247926*pow2(Mst1)*pow4(Msq) -
        14560*pow2(Msq)*pow4(Mst1) - 5*(321*pow2(Mst2)*pow4(Mst1) + 216*pow2(
        Mst1)*pow4(Mst2) + 32*pow6(Mst2)))) + 630*log(pow2(Mst2)/pow2(Mst1))*(-
        24500*Dmglst1*Mgl*(232*pow2(Msq) - 167*pow2(Mst2))*pow4(Mst1) - 159250*
        pow2(Dmglst1)*(152*pow2(Msq) - 67*pow2(Mst2))*pow4(Mst1) + pow2(Mgl)*(
        73001572*pow2(Mst1)*pow4(Msq) - 3778880*pow2(Msq)*pow4(Mst1) - 5*(5151*
        pow2(Mst2)*pow4(Mst1) + 126132*pow2(Mst1)*pow4(Mst2) + 31294*pow6(Mst2)
        )))))/pow4(Msq)))/(1.000188e9*pow2(Mgl)*pow4(Mst2))) + z2*(pow2(s2t)*((
        -5*pow2(Msq)*(pow2(Mst1) - pow2(Mst2))*(2*log(pow2(Mst2)/pow2(Mst1))*
        pow2(Mst1)*pow2(Mst2) + pow4(Mst1) - pow4(Mst2)))/(3.*pow2(Mst1)*pow2(
        Mst2)) - ((1080*pow2(Mst1)*(pow2(Dmglst1)*(100*pow2(Mst1) - 9*pow2(
        Mst2)) + Dmglst1*Mgl*(40*pow2(Mst1) - 6*pow2(Mst2)) + pow2(Mgl)*(10*
        pow2(Mst1) - 3*pow2(Mst2)))*pow2(Mst2))/pow2(Msq) + 18*pow2(Dmglst1)*(-
        108*(97 + 30*log(pow2(Msq)/pow2(Mst1)) - 105*log(pow2(Mst2)/pow2(Mst1))
        )*pow2(Mst1)*pow2(Mst2) + (29201 + 3240*log(pow2(Msq)/pow2(Mst1)) -
        5388*log(pow2(Mst2)/pow2(Mst1)))*pow4(Mst1) - 1782*pow4(Mst2)) + 36*
        Dmglst1*Mgl*(-3492*pow2(Mst1)*pow2(Mst2) + 7153*pow4(Mst1) + 1080*log(
        pow2(Msq)/pow2(Mst1))*(-(pow2(Mst1)*pow2(Mst2)) + pow4(Mst1)) - 252*
        log(pow2(Mst2)/pow2(Mst1))*(-15*pow2(Mst1)*pow2(Mst2) + 7*pow4(Mst1)) -
        594*pow4(Mst2)) - (1350*(4*Dmglst1*Mgl + 10*pow2(Dmglst1) + pow2(Mgl))*
        pow4(Mst1)*pow4(Mst2))/pow4(Msq) + pow2(Mgl)*(-68844*pow2(Mst1)*pow2(
        Mst2) + 98497*pow4(Mst1) + 108*log(pow2(Mst2)/pow2(Mst1))*(670*pow2(
        Mst1)*pow2(Mst2) + 103*pow4(Mst1) - 149*pow4(Mst2)) + 35316*pow4(Mst2)
        + 3240*log(pow2(Msq)/pow2(Mst1))*(-6*pow2(Mst1)*pow2(Mst2) + 5*pow4(
        Mst1) + pow4(Mst2))))/(3888.*pow2(Mgl))) + (Mt*(-18*s2t*pow2(Mst2)*
        pow3(Mst1)*(-3*pow2(Dmglst1)*(-450*pow4(Mst1)*pow4(Mst2) + pow4(Msq)*(-
        8630*pow2(Mst1)*pow2(Mst2) + 445*pow4(Mst1) + 4*log(pow2(Mst2)/pow2(
        Mst1))*(1369*pow2(Mst1)*pow2(Mst2) + 1129*pow4(Mst1) - 96*pow4(Mst2)) -
        120*log(pow2(Msq)/pow2(Mst1))*(23*pow2(Mst1)*pow2(Mst2) - pow4(Mst2)) +
        1254*pow4(Mst2)) + 240*pow2(Msq)*(13*pow2(Mst2)*pow4(Mst1) - 2*pow2(
        Mst1)*pow4(Mst2))) + pow2(Mgl)*(pow4(Msq)*(2710*pow2(Mst1)*pow2(Mst2) -
        737*pow4(Mst1) - 12*log(pow2(Mst2)/pow2(Mst1))*(124*pow2(Mst1)*pow2(
        Mst2) + 71*pow4(Mst1) - 96*pow4(Mst2)) + 360*log(pow2(Msq)/pow2(Mst1))*
        (3*pow2(Mst1)*pow2(Mst2) - pow4(Mst2)) - 3762*pow4(Mst2)) + 90*pow4(
        Mst1)*pow4(Mst2) - 240*pow2(Msq)*(2*pow2(Mst2)*pow4(Mst1) - pow2(Mst1)*
        pow4(Mst2))) + Dmglst1*Mgl*(pow4(Msq)*(11646*pow2(Mst1)*pow2(Mst2) -
        1397*pow4(Mst1) - 12*log(pow2(Mst2)/pow2(Mst1))*(622*pow2(Mst1)*pow2(
        Mst2) + 427*pow4(Mst1) - 96*pow4(Mst2)) + 360*log(pow2(Msq)/pow2(Mst1))
        *(11*pow2(Mst1)*pow2(Mst2) - pow4(Mst2)) - 3762*pow4(Mst2)) + 450*pow4(
        Mst1)*pow4(Mst2) - 720*pow2(Msq)*(4*pow2(Mst2)*pow4(Mst1) - pow2(Mst1)*
        pow4(Mst2)))) + Mt*(18*Dmglst1*pow2(Mst1)*(Dmglst1*(30*pow2(MuSUSY)*
        pow4(Mst1)*(164*pow2(Msq)*pow2(Mst2) - 25*pow4(Mst2)) + pow4(Msq)*((4*(
        -7429 + 8538*log(pow2(Mst2)/pow2(Mst1)))*pow2(Mst2) + (2903 - 3240*log(
        pow2(Msq)/pow2(Mst1)) + 17292*log(pow2(Mst2)/pow2(Mst1)))*pow2(MuSUSY))
        *pow4(Mst1) + 12*(1537 - 228*log(pow2(Mst2)/pow2(Mst1)))*pow2(Mst1)*
        pow4(Mst2) - 1782*pow6(Mst2))) + 2*Mgl*(30*pow2(MuSUSY)*pow4(Mst1)*(28*
        pow2(Msq)*pow2(Mst2) - 5*pow4(Mst2)) + pow4(Msq)*((-3692*pow2(Mst2) -
        1613*pow2(MuSUSY) - 1080*log(pow2(Msq)/pow2(Mst1))*pow2(MuSUSY) + 84*
        log(pow2(Mst2)/pow2(Mst1))*(70*pow2(Mst2) + 69*pow2(MuSUSY)))*pow4(
        Mst1) + 4*(1537 - 228*log(pow2(Mst2)/pow2(Mst1)))*pow2(Mst1)*pow4(Mst2)
        - 594*pow6(Mst2)))) - pow2(Mgl)*(-1080*pow2(Msq)*(4*pow2(Mst1) - 3*
        pow2(Mst2))*pow2(Mst2)*pow2(MuSUSY)*pow4(Mst1) - 6480*(pow2(-(Mst2*
        pow2(Mst1)) + pow3(Mst2)) - 2*log(pow2(Mst2)/pow2(Mst1))*pow2(MuSUSY)*
        pow4(Mst1))*pow6(Msq) + 1350*pow2(MuSUSY)*pow4(Mst2)*pow6(Mst1) - pow4(
        Msq)*(3732*pow2(Mst2)*pow2(MuSUSY)*pow4(Mst1) - 12960*log(pow2(Msq)/
        pow2(Mst1))*(pow2(Mst1) + pow2(Mst2))*pow2(MuSUSY)*pow4(Mst1) + 120855*
        pow4(Mst1)*pow4(Mst2) + 590*pow2(Mst2)*pow6(Mst1) + 69349*pow2(MuSUSY)*
        pow6(Mst1) - 16497*pow2(Mst1)*pow6(Mst2) + 216*log(pow2(Mst2)/pow2(
        Mst1))*(pow4(Mst1)*(186*pow2(Mst2)*pow2(MuSUSY) - 40*pow4(Mst2)) + 3*(
        55*pow2(Mst2) + 166*pow2(MuSUSY))*pow6(Mst1) + 3*pow2(Mst1)*pow6(Mst2))
        + 648*pow8(Mst2))))))/(486.*pow2(Mgl)*pow2(Mst1)*pow4(Msq)*pow4(Mst2)))
        )))/pow4(Mt)/pow2(Sbeta)*12.;

   return result;
}

/**
 * 	@return returns the susy log^1 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double himalaya::H5g1::calc_coef_at_as2_no_sm_logs_log1(){

   const double result =
      ((1296000*pow2(Mgl)*pow2(Mst2)*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*pow4(Msq)*
        pow4(Mst1) + 28800*z2*pow2(Mgl)*pow2(Mst2)*pow2(Mt)*pow2(MuSUSY)*pow2(
        s2t)*pow4(Msq)*pow4(Mst1) - 1296000*pow2(Mgl)*pow2(Mst2)*pow2(Mt)*pow2(
        MuSUSY)*pow2(s2t)*pow2(Sbeta)*pow4(Msq)*pow4(Mst1) - 28800*z2*pow2(Mgl)
        *pow2(Mst2)*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*pow2(Sbeta)*pow4(Msq)*pow4(
        Mst1) - 6753600*Cbeta*MuSUSY*s2t*Sbeta*pow2(Mgl)*pow2(Mst2)*pow3(Mt)*
        pow4(Msq)*pow4(Mst1) - 1382400*Cbeta*MuSUSY*s2t*Sbeta*z2*pow2(Mgl)*
        pow2(Mst2)*pow3(Mt)*pow4(Msq)*pow4(Mst1) - 2009600*Dmglst1*Mgl*s2t*
        pow2(Sbeta)*pow3(Mst1)*pow3(Mt)*pow4(Msq)*pow4(Mst2) + 403200*Dmglst1*
        Mgl*s2t*z2*pow2(Sbeta)*pow3(Mst1)*pow3(Mt)*pow4(Msq)*pow4(Mst2) -
        6019200*s2t*pow2(Dmglst1)*pow2(Sbeta)*pow3(Mst1)*pow3(Mt)*pow4(Msq)*
        pow4(Mst2) + 403200*s2t*z2*pow2(Dmglst1)*pow2(Sbeta)*pow3(Mst1)*pow3(
        Mt)*pow4(Msq)*pow4(Mst2) + 662400*s2t*pow2(Mgl)*pow2(Sbeta)*pow3(Mst1)*
        pow3(Mt)*pow4(Msq)*pow4(Mst2) + 403200*s2t*z2*pow2(Mgl)*pow2(Sbeta)*
        pow3(Mst1)*pow3(Mt)*pow4(Msq)*pow4(Mst2) + 7922400*Dmglst1*Mgl*pow2(Mt)
        *pow2(s2t)*pow2(Sbeta)*pow4(Msq)*pow4(Mst1)*pow4(Mst2) + 1382400*
        Dmglst1*Mgl*z2*pow2(Mt)*pow2(s2t)*pow2(Sbeta)*pow4(Msq)*pow4(Mst1)*
        pow4(Mst2) + 12982800*pow2(Dmglst1)*pow2(Mt)*pow2(s2t)*pow2(Sbeta)*
        pow4(Msq)*pow4(Mst1)*pow4(Mst2) + 2073600*z2*pow2(Dmglst1)*pow2(Mt)*
        pow2(s2t)*pow2(Sbeta)*pow4(Msq)*pow4(Mst1)*pow4(Mst2) + 3636000*pow2(
        Mgl)*pow2(Mt)*pow2(s2t)*pow2(Sbeta)*pow4(Msq)*pow4(Mst1)*pow4(Mst2) +
        691200*z2*pow2(Mgl)*pow2(Mt)*pow2(s2t)*pow2(Sbeta)*pow4(Msq)*pow4(Mst1)
        *pow4(Mst2) + 561600*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mgl)*pow3(s2t)*pow4(
        Msq)*pow4(Mst1)*pow4(Mst2) - 162000*Cbeta*Mt*MuSUSY*Sbeta*z2*pow2(Mgl)*
        pow3(s2t)*pow4(Msq)*pow4(Mst1)*pow4(Mst2) + 3456000*pow2(Mgl)*pow2(
        MuSUSY)*pow4(Msq)*pow4(Mst1)*pow4(Mt) + 921600*z2*pow2(Mgl)*pow2(
        MuSUSY)*pow4(Msq)*pow4(Mst1)*pow4(Mt) + 2947200*Dmglst1*Mgl*pow2(Mst2)*
        pow2(Sbeta)*pow4(Msq)*pow4(Mst1)*pow4(Mt) + 2073600*Dmglst1*Mgl*z2*
        pow2(Mst2)*pow2(Sbeta)*pow4(Msq)*pow4(Mst1)*pow4(Mt) + 5030560*pow2(
        Dmglst1)*pow2(Mst2)*pow2(Sbeta)*pow4(Msq)*pow4(Mst1)*pow4(Mt) +
        3110400*z2*pow2(Dmglst1)*pow2(Mst2)*pow2(Sbeta)*pow4(Msq)*pow4(Mst1)*
        pow4(Mt) + 201600*pow2(Mgl)*pow2(Mst2)*pow2(Sbeta)*pow4(Msq)*pow4(Mst1)
        *pow4(Mt) + 1065600*z2*pow2(Mgl)*pow2(Mst2)*pow2(Sbeta)*pow4(Msq)*pow4(
        Mst1)*pow4(Mt) - 3456000*pow2(Mgl)*pow2(MuSUSY)*pow2(Sbeta)*pow4(Msq)*
        pow4(Mst1)*pow4(Mt) - 921600*z2*pow2(Mgl)*pow2(MuSUSY)*pow2(Sbeta)*
        pow4(Msq)*pow4(Mst1)*pow4(Mt) - 551200*Dmglst1*Mgl*pow2(Mst1)*pow2(
        Sbeta)*pow4(Msq)*pow4(Mst2)*pow4(Mt) - 1166752*pow2(Dmglst1)*pow2(Mst1)
        *pow2(Sbeta)*pow4(Msq)*pow4(Mst2)*pow4(Mt) - 604400*pow2(Mgl)*pow2(
        Mst1)*pow2(Sbeta)*pow4(Msq)*pow4(Mst2)*pow4(Mt) + 597600*z2*pow2(Mgl)*
        pow2(Mst1)*pow2(Sbeta)*pow4(Msq)*pow4(Mst2)*pow4(Mt) + 172800*z3*pow2(
        Mgl)*pow2(Mst1)*pow2(Sbeta)*pow4(Msq)*pow4(Mst2)*pow4(Mt) - 216000*
        Dmglst1*Mgl*pow2(Msq)*pow2(Sbeta)*pow4(Mst1)*pow4(Mst2)*pow4(Mt) -
        50400*pow2(Dmglst1)*pow2(Msq)*pow2(Sbeta)*pow4(Mst1)*pow4(Mst2)*pow4(
        Mt) - 396000*pow2(Mgl)*pow2(Msq)*pow2(Sbeta)*pow4(Mst1)*pow4(Mst2)*
        pow4(Mt) + 302400*Cbeta*Dmglst1*Mgl*MuSUSY*Sbeta*pow2(Mst2)*pow2(Mt)*
        pow2(s2t)*pow4(Msq)*pow5(Mst1) + 3196800*Cbeta*Dmglst1*Mgl*MuSUSY*
        Sbeta*z2*pow2(Mst2)*pow2(Mt)*pow2(s2t)*pow4(Msq)*pow5(Mst1) + 302400*
        Cbeta*MuSUSY*Sbeta*pow2(Dmglst1)*pow2(Mst2)*pow2(Mt)*pow2(s2t)*pow4(
        Msq)*pow5(Mst1) + 3196800*Cbeta*MuSUSY*Sbeta*z2*pow2(Dmglst1)*pow2(
        Mst2)*pow2(Mt)*pow2(s2t)*pow4(Msq)*pow5(Mst1) + 43200*Cbeta*MuSUSY*
        Sbeta*pow2(Mgl)*pow2(Mst2)*pow2(Mt)*pow2(s2t)*pow4(Msq)*pow5(Mst1) +
        86400*Cbeta*MuSUSY*Sbeta*z2*pow2(Mgl)*pow2(Mst2)*pow2(Mt)*pow2(s2t)*
        pow4(Msq)*pow5(Mst1) + 11232000*Dmglst1*Mgl*s2t*pow2(MuSUSY)*pow3(Mt)*
        pow4(Msq)*pow5(Mst1) - 3283200*Dmglst1*Mgl*s2t*z2*pow2(MuSUSY)*pow3(Mt)
        *pow4(Msq)*pow5(Mst1) + 11232000*s2t*pow2(Dmglst1)*pow2(MuSUSY)*pow3(
        Mt)*pow4(Msq)*pow5(Mst1) - 3283200*s2t*z2*pow2(Dmglst1)*pow2(MuSUSY)*
        pow3(Mt)*pow4(Msq)*pow5(Mst1) + 9561600*s2t*pow2(Mgl)*pow2(MuSUSY)*
        pow3(Mt)*pow4(Msq)*pow5(Mst1) + 864000*s2t*z2*pow2(Mgl)*pow2(MuSUSY)*
        pow3(Mt)*pow4(Msq)*pow5(Mst1) + 14745600*Dmglst1*Mgl*s2t*pow2(Mst2)*
        pow2(Sbeta)*pow3(Mt)*pow4(Msq)*pow5(Mst1) + 11347200*Dmglst1*Mgl*s2t*
        z2*pow2(Mst2)*pow2(Sbeta)*pow3(Mt)*pow4(Msq)*pow5(Mst1) + 35702400*s2t*
        pow2(Dmglst1)*pow2(Mst2)*pow2(Sbeta)*pow3(Mt)*pow4(Msq)*pow5(Mst1) +
        22492800*s2t*z2*pow2(Dmglst1)*pow2(Mst2)*pow2(Sbeta)*pow3(Mt)*pow4(Msq)
        *pow5(Mst1) + 1843200*s2t*pow2(Mgl)*pow2(Mst2)*pow2(Sbeta)*pow3(Mt)*
        pow4(Msq)*pow5(Mst1) + 3916800*s2t*z2*pow2(Mgl)*pow2(Mst2)*pow2(Sbeta)*
        pow3(Mt)*pow4(Msq)*pow5(Mst1) - 11232000*Dmglst1*Mgl*s2t*pow2(MuSUSY)*
        pow2(Sbeta)*pow3(Mt)*pow4(Msq)*pow5(Mst1) + 3283200*Dmglst1*Mgl*s2t*z2*
        pow2(MuSUSY)*pow2(Sbeta)*pow3(Mt)*pow4(Msq)*pow5(Mst1) - 11232000*s2t*
        pow2(Dmglst1)*pow2(MuSUSY)*pow2(Sbeta)*pow3(Mt)*pow4(Msq)*pow5(Mst1) +
        3283200*s2t*z2*pow2(Dmglst1)*pow2(MuSUSY)*pow2(Sbeta)*pow3(Mt)*pow4(
        Msq)*pow5(Mst1) - 9561600*s2t*pow2(Mgl)*pow2(MuSUSY)*pow2(Sbeta)*pow3(
        Mt)*pow4(Msq)*pow5(Mst1) - 864000*s2t*z2*pow2(Mgl)*pow2(MuSUSY)*pow2(
        Sbeta)*pow3(Mt)*pow4(Msq)*pow5(Mst1) - 2304000*Dmglst1*Mgl*s2t*pow2(
        Msq)*pow2(Sbeta)*pow3(Mt)*pow4(Mst2)*pow5(Mst1) - 6336000*s2t*pow2(
        Dmglst1)*pow2(Msq)*pow2(Sbeta)*pow3(Mt)*pow4(Mst2)*pow5(Mst1) - 192000*
        s2t*pow2(Mgl)*pow2(Msq)*pow2(Sbeta)*pow3(Mt)*pow4(Mst2)*pow5(Mst1) -
        3009600*Dmglst1*Mgl*Mt*pow2(Sbeta)*pow3(s2t)*pow4(Msq)*pow4(Mst2)*pow5(
        Mst1) - 1310400*Dmglst1*Mgl*Mt*z2*pow2(Sbeta)*pow3(s2t)*pow4(Msq)*pow4(
        Mst2)*pow5(Mst1) - 5224800*Mt*pow2(Dmglst1)*pow2(Sbeta)*pow3(s2t)*pow4(
        Msq)*pow4(Mst2)*pow5(Mst1) - 2865600*Mt*z2*pow2(Dmglst1)*pow2(Sbeta)*
        pow3(s2t)*pow4(Msq)*pow4(Mst2)*pow5(Mst1) - 2419200*Mt*pow2(Mgl)*pow2(
        Sbeta)*pow3(s2t)*pow4(Msq)*pow4(Mst2)*pow5(Mst1) - 273600*Mt*z2*pow2(
        Mgl)*pow2(Sbeta)*pow3(s2t)*pow4(Msq)*pow4(Mst2)*pow5(Mst1) + 2832000*
        Cbeta*Dmglst1*Mgl*MuSUSY*Sbeta*pow2(Msq)*pow2(Mst2)*pow4(Mt)*pow5(Mst1)
        + 2832000*Cbeta*MuSUSY*Sbeta*pow2(Dmglst1)*pow2(Msq)*pow2(Mst2)*pow4(
        Mt)*pow5(Mst1) + 720000*Cbeta*MuSUSY*Sbeta*pow2(Mgl)*pow2(Msq)*pow2(
        Mst2)*pow4(Mt)*pow5(Mst1) - 11531200*Cbeta*Dmglst1*Mgl*MuSUSY*Sbeta*
        pow4(Msq)*pow4(Mt)*pow5(Mst1) - 11750400*Cbeta*Dmglst1*Mgl*MuSUSY*
        Sbeta*z2*pow4(Msq)*pow4(Mt)*pow5(Mst1) - 11531200*Cbeta*MuSUSY*Sbeta*
        pow2(Dmglst1)*pow4(Msq)*pow4(Mt)*pow5(Mst1) - 11750400*Cbeta*MuSUSY*
        Sbeta*z2*pow2(Dmglst1)*pow4(Msq)*pow4(Mt)*pow5(Mst1) - 2203200*Cbeta*
        MuSUSY*Sbeta*pow2(Mgl)*pow4(Msq)*pow4(Mt)*pow5(Mst1) - 4320000*Cbeta*
        MuSUSY*Sbeta*z2*pow2(Mgl)*pow4(Msq)*pow4(Mt)*pow5(Mst1) + 546000*Cbeta*
        Dmglst1*Mgl*MuSUSY*Sbeta*pow4(Mst2)*pow4(Mt)*pow5(Mst1) + 546000*Cbeta*
        MuSUSY*Sbeta*pow2(Dmglst1)*pow4(Mst2)*pow4(Mt)*pow5(Mst1) + 234000*
        Cbeta*MuSUSY*Sbeta*pow2(Mgl)*pow4(Mst2)*pow4(Mt)*pow5(Mst1) - 3600*
        Sbeta*log(pow2(Msq)/pow2(Mst1))*pow2(Mst1)*pow3(Mt)*(20*Cbeta*Mt*
        MuSUSY*pow3(Mst1)*(Dmglst1*Mgl*(118*pow4(Msq) - 7*pow4(Mst2)) + pow2(
        Dmglst1)*(118*pow4(Msq) - 7*pow4(Mst2)) + 3*pow2(Mgl)*(10*pow4(Msq) -
        pow4(Mst2))) + Sbeta*(pow2(Dmglst1)*(-40*Mst1*pow2(Msq)*(-9*Mst1*Mt +
        2*s2t*pow2(Mst1) - 2*s2t*pow2(Mst2))*pow4(Mst2) + 10*Mst1*pow4(Mst2)*(
        9*Mst1*Mt*pow2(Mst2) + 24*s2t*pow2(Mst1)*pow2(Mst2) + 59*Mt*pow3(Mst1)
        - 26*s2t*pow4(Mst1) + 2*s2t*pow4(Mst2)) + pow4(Msq)*(540*Mt*pow2(Mst1)*
        pow2(Mst2) - 4800*s2t*pow2(Mst2)*pow3(Mst1) + 6900*Mt*pow4(Mst1) + 68*
        Mt*pow4(Mst2) - 920*Mst1*s2t*pow4(Mst2) - 11520*s2t*pow5(Mst1))) - 5*
        pow2(Mgl)*(-8*pow2(Msq)*(4*Mt*pow2(Mst1) + Mt*pow2(Mst2) + 2*Mst1*s2t*
        pow2(Mst2) - 2*s2t*pow3(Mst1))*pow4(Mst2) - pow4(Mst2)*(6*Mt*pow2(Mst1)
        *pow2(Mst2) + 8*s2t*pow2(Mst2)*pow3(Mst1) + 19*Mt*pow4(Mst1) + 3*Mt*
        pow4(Mst2) + 4*Mst1*s2t*pow4(Mst2) - 12*s2t*pow5(Mst1)) + 6*pow4(Msq)*(
        -6*Mt*pow2(Mst1)*pow2(Mst2) + 8*s2t*pow2(Mst2)*pow3(Mst1) - 14*Mt*pow4(
        Mst1) - 7*Mt*pow4(Mst2) + 12*Mst1*s2t*pow4(Mst2) + 8*s2t*pow5(Mst1))) -
        20*Dmglst1*Mgl*(4*Mst1*pow2(Msq)*(-3*Mst1*Mt + s2t*pow2(Mst1) - s2t*
        pow2(Mst2))*pow4(Mst2) + Mst1*pow4(Mst2)*(-3*Mst1*Mt*pow2(Mst2) - 6*
        s2t*pow2(Mst1)*pow2(Mst2) - 13*Mt*pow3(Mst1) + 7*s2t*pow4(Mst1) - s2t*
        pow4(Mst2)) + 2*pow4(Msq)*(-9*Mt*pow2(Mst1)*pow2(Mst2) + 42*s2t*pow2(
        Mst2)*pow3(Mst1) - 57*Mt*pow4(Mst1) - 2*Mt*pow4(Mst2) + 17*Mst1*s2t*
        pow4(Mst2) + 66*s2t*pow5(Mst1)))) - 60*log(pow2(Mst2)/pow2(Mst1))*pow4(
        Msq)*(4*Cbeta*Mt*MuSUSY*(7*Dmglst1*Mgl + 7*pow2(Dmglst1) + 3*pow2(Mgl))
        *pow3(Mst1) - Sbeta*(4*Dmglst1*Mgl*Mst1*(6*s2t*pow2(Mst1)*pow2(Mst2) -
        5*Mt*pow3(Mst1) + 15*s2t*pow4(Mst1) + s2t*pow4(Mst2)) + pow2(Mgl)*(8*
        s2t*pow2(Mst2)*pow3(Mst1) - 5*Mt*pow4(Mst1) + Mt*pow4(Mst2) + 4*Mst1*
        s2t*pow4(Mst2) + 12*s2t*pow5(Mst1)) + 2*pow2(Dmglst1)*(24*s2t*pow2(
        Mst2)*pow3(Mst1) - 25*Mt*pow4(Mst1) + 2*Mst1*s2t*pow4(Mst2) + 90*s2t*
        pow5(Mst1))))) + 216000*pow2(Mgl)*pow2(Mst2)*pow2(Mt)*pow2(s2t)*pow2(
        Sbeta)*pow4(Mst1)*pow6(Msq) - 108000*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mgl)*
        pow2(Mst2)*pow3(s2t)*pow4(Mst1)*pow6(Msq) - 432000*pow2(Mgl)*pow2(Mst1)
        *pow2(Mt)*pow2(s2t)*pow2(Sbeta)*pow4(Mst2)*pow6(Msq) - 432000*pow2(Mgl)
        *pow2(Mst1)*pow2(Mst2)*pow2(Sbeta)*pow4(Mt)*pow6(Msq) - 432000*pow2(
        Mgl)*pow2(Sbeta)*pow4(Mst2)*pow4(Mt)*pow6(Msq) + 27000*pow2(Mgl)*pow2(
        Sbeta)*pow4(Mst1)*pow4(Mst2)*pow4(s2t)*pow6(Msq) - 1382400*Dmglst1*Mgl*
        pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*pow4(Msq)*pow6(Mst1) + 1238400*Dmglst1*
        Mgl*z2*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*pow4(Msq)*pow6(Mst1) + 158400*
        pow2(Dmglst1)*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*pow4(Msq)*pow6(Mst1) +
        1857600*z2*pow2(Dmglst1)*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*pow4(Msq)*
        pow6(Mst1) - 70200*pow2(Mgl)*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*pow4(Msq)*
        pow6(Mst1) + 28800*z2*pow2(Mgl)*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*pow4(
        Msq)*pow6(Mst1) - 2980800*Dmglst1*Mgl*pow2(Mst2)*pow2(Mt)*pow2(s2t)*
        pow2(Sbeta)*pow4(Msq)*pow6(Mst1) - 2764800*Dmglst1*Mgl*z2*pow2(Mst2)*
        pow2(Mt)*pow2(s2t)*pow2(Sbeta)*pow4(Msq)*pow6(Mst1) - 8460000*pow2(
        Dmglst1)*pow2(Mst2)*pow2(Mt)*pow2(s2t)*pow2(Sbeta)*pow4(Msq)*pow6(Mst1)
        - 9676800*z2*pow2(Dmglst1)*pow2(Mst2)*pow2(Mt)*pow2(s2t)*pow2(Sbeta)*
        pow4(Msq)*pow6(Mst1) - 1576800*pow2(Mgl)*pow2(Mst2)*pow2(Mt)*pow2(s2t)*
        pow2(Sbeta)*pow4(Msq)*pow6(Mst1) + 1382400*Dmglst1*Mgl*pow2(Mt)*pow2(
        MuSUSY)*pow2(s2t)*pow2(Sbeta)*pow4(Msq)*pow6(Mst1) - 1238400*Dmglst1*
        Mgl*z2*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*pow2(Sbeta)*pow4(Msq)*pow6(Mst1)
        - 158400*pow2(Dmglst1)*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*pow2(Sbeta)*
        pow4(Msq)*pow6(Mst1) - 1857600*z2*pow2(Dmglst1)*pow2(Mt)*pow2(MuSUSY)*
        pow2(s2t)*pow2(Sbeta)*pow4(Msq)*pow6(Mst1) + 70200*pow2(Mgl)*pow2(Mt)*
        pow2(MuSUSY)*pow2(s2t)*pow2(Sbeta)*pow4(Msq)*pow6(Mst1) - 28800*z2*
        pow2(Mgl)*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*pow2(Sbeta)*pow4(Msq)*pow6(
        Mst1) - 9705600*Cbeta*Dmglst1*Mgl*MuSUSY*s2t*Sbeta*pow3(Mt)*pow4(Msq)*
        pow6(Mst1) + 2764800*Cbeta*Dmglst1*Mgl*MuSUSY*s2t*Sbeta*z2*pow3(Mt)*
        pow4(Msq)*pow6(Mst1) - 9460800*Cbeta*MuSUSY*s2t*Sbeta*pow2(Dmglst1)*
        pow3(Mt)*pow4(Msq)*pow6(Mst1) + 15206400*Cbeta*MuSUSY*s2t*Sbeta*z2*
        pow2(Dmglst1)*pow3(Mt)*pow4(Msq)*pow6(Mst1) - 3600000*Cbeta*MuSUSY*s2t*
        Sbeta*pow2(Mgl)*pow3(Mt)*pow4(Msq)*pow6(Mst1) - 1382400*Cbeta*MuSUSY*
        s2t*Sbeta*z2*pow2(Mgl)*pow3(Mt)*pow4(Msq)*pow6(Mst1) + 1447200*Cbeta*
        Dmglst1*Mgl*Mt*MuSUSY*Sbeta*pow2(Mst2)*pow3(s2t)*pow4(Msq)*pow6(Mst1) +
        1458000*Cbeta*Mt*MuSUSY*Sbeta*pow2(Dmglst1)*pow2(Mst2)*pow3(s2t)*pow4(
        Msq)*pow6(Mst1) + 683100*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mgl)*pow2(Mst2)*
        pow3(s2t)*pow4(Msq)*pow6(Mst1) + 1440000*Dmglst1*Mgl*pow2(Msq)*pow2(
        Mst2)*pow2(Sbeta)*pow4(Mt)*pow6(Mst1) + 3600000*pow2(Dmglst1)*pow2(Msq)
        *pow2(Mst2)*pow2(Sbeta)*pow4(Mt)*pow6(Mst1) + 360000*pow2(Mgl)*pow2(
        Msq)*pow2(Mst2)*pow2(Sbeta)*pow4(Mt)*pow6(Mst1) - 11539200*Dmglst1*Mgl*
        pow2(Sbeta)*pow4(Msq)*pow4(Mt)*pow6(Mst1) - 345600*Dmglst1*Mgl*z2*pow2(
        Sbeta)*pow4(Msq)*pow4(Mt)*pow6(Mst1) - 20852960*pow2(Dmglst1)*pow2(
        Sbeta)*pow4(Msq)*pow4(Mt)*pow6(Mst1) - 2505600*z2*pow2(Dmglst1)*pow2(
        Sbeta)*pow4(Msq)*pow4(Mt)*pow6(Mst1) - 3720600*pow2(Mgl)*pow2(Sbeta)*
        pow4(Msq)*pow4(Mt)*pow6(Mst1) + 345600*z2*pow2(Mgl)*pow2(Sbeta)*pow4(
        Msq)*pow4(Mt)*pow6(Mst1) + 114000*Dmglst1*Mgl*pow2(Sbeta)*pow4(Mst2)*
        pow4(Mt)*pow6(Mst1) + 587400*pow2(Dmglst1)*pow2(Sbeta)*pow4(Mst2)*pow4(
        Mt)*pow6(Mst1) - 61500*pow2(Mgl)*pow2(Sbeta)*pow4(Mst2)*pow4(Mt)*pow6(
        Mst1) - 534600*Dmglst1*Mgl*pow2(Sbeta)*pow4(Msq)*pow4(Mst2)*pow4(s2t)*
        pow6(Mst1) - 154800*Dmglst1*Mgl*z2*pow2(Sbeta)*pow4(Msq)*pow4(Mst2)*
        pow4(s2t)*pow6(Mst1) - 735000*pow2(Dmglst1)*pow2(Sbeta)*pow4(Msq)*pow4(
        Mst2)*pow4(s2t)*pow6(Mst1) - 232200*z2*pow2(Dmglst1)*pow2(Sbeta)*pow4(
        Msq)*pow4(Mst2)*pow4(s2t)*pow6(Mst1) - 30375*pow2(Mgl)*pow2(Sbeta)*
        pow4(Msq)*pow4(Mst2)*pow4(s2t)*pow6(Mst1) - 40500*z2*pow2(Mgl)*pow2(
        Sbeta)*pow4(Msq)*pow4(Mst2)*pow4(s2t)*pow6(Mst1) - 27000*pow2(Mgl)*
        pow2(Mst2)*pow2(Sbeta)*pow4(s2t)*pow6(Msq)*pow6(Mst1) - 528000*Dmglst1*
        Mgl*s2t*pow2(Msq)*pow2(Sbeta)*pow3(Mst1)*pow3(Mt)*pow6(Mst2) - 528000*
        s2t*pow2(Dmglst1)*pow2(Msq)*pow2(Sbeta)*pow3(Mst1)*pow3(Mt)*pow6(Mst2)
        - 528000*s2t*pow2(Mgl)*pow2(Msq)*pow2(Sbeta)*pow3(Mst1)*pow3(Mt)*pow6(
        Mst2) - 88800*Dmglst1*Mgl*pow2(Mst1)*pow2(Mt)*pow2(s2t)*pow2(Sbeta)*
        pow4(Msq)*pow6(Mst2) + 207600*pow2(Dmglst1)*pow2(Mst1)*pow2(Mt)*pow2(
        s2t)*pow2(Sbeta)*pow4(Msq)*pow6(Mst2) - 280800*pow2(Mgl)*pow2(Mst1)*
        pow2(Mt)*pow2(s2t)*pow2(Sbeta)*pow4(Msq)*pow6(Mst2) + 2908800*Dmglst1*
        Mgl*Mt*pow2(Sbeta)*pow3(Mst1)*pow3(s2t)*pow4(Msq)*pow6(Mst2) + 244800*
        Dmglst1*Mgl*Mt*z2*pow2(Sbeta)*pow3(Mst1)*pow3(s2t)*pow4(Msq)*pow6(Mst2)
        + 2949600*Mt*pow2(Dmglst1)*pow2(Sbeta)*pow3(Mst1)*pow3(s2t)*pow4(Msq)*
        pow6(Mst2) + 244800*Mt*z2*pow2(Dmglst1)*pow2(Sbeta)*pow3(Mst1)*pow3(
        s2t)*pow4(Msq)*pow6(Mst2) + 2404800*Mt*pow2(Mgl)*pow2(Sbeta)*pow3(Mst1)
        *pow3(s2t)*pow4(Msq)*pow6(Mst2) + 244800*Mt*z2*pow2(Mgl)*pow2(Sbeta)*
        pow3(Mst1)*pow3(s2t)*pow4(Msq)*pow6(Mst2) - 192000*pow2(Mgl)*pow2(Msq)*
        pow2(Mst1)*pow2(Sbeta)*pow4(Mt)*pow6(Mst2) - 43200*pow2(Mgl)*pow2(
        Sbeta)*pow4(Msq)*pow4(Mt)*pow6(Mst2) - 126000*Dmglst1*Mgl*pow2(Sbeta)*
        pow4(Mst1)*pow4(Mt)*pow6(Mst2) - 189000*pow2(Dmglst1)*pow2(Sbeta)*pow4(
        Mst1)*pow4(Mt)*pow6(Mst2) - 63000*pow2(Mgl)*pow2(Sbeta)*pow4(Mst1)*
        pow4(Mt)*pow6(Mst2) + 156600*Dmglst1*Mgl*pow2(Sbeta)*pow4(Msq)*pow4(
        Mst1)*pow4(s2t)*pow6(Mst2) + 154800*Dmglst1*Mgl*z2*pow2(Sbeta)*pow4(
        Msq)*pow4(Mst1)*pow4(s2t)*pow6(Mst2) + 356700*pow2(Dmglst1)*pow2(Sbeta)
        *pow4(Msq)*pow4(Mst1)*pow4(s2t)*pow6(Mst2) + 232200*z2*pow2(Dmglst1)*
        pow2(Sbeta)*pow4(Msq)*pow4(Mst1)*pow4(s2t)*pow6(Mst2) - 445500*pow2(
        Mgl)*pow2(Sbeta)*pow4(Msq)*pow4(Mst1)*pow4(s2t)*pow6(Mst2) + 77400*z2*
        pow2(Mgl)*pow2(Sbeta)*pow4(Msq)*pow4(Mst1)*pow4(s2t)*pow6(Mst2) -
        468000*Dmglst1*Mgl*s2t*pow2(Sbeta)*pow3(Mt)*pow5(Mst1)*pow6(Mst2) -
        936000*s2t*pow2(Dmglst1)*pow2(Sbeta)*pow3(Mt)*pow5(Mst1)*pow6(Mst2) -
        156000*s2t*pow2(Mgl)*pow2(Sbeta)*pow3(Mt)*pow5(Mst1)*pow6(Mst2) +
        216000*pow2(Mgl)*pow2(Mt)*pow2(s2t)*pow2(Sbeta)*pow6(Msq)*pow6(Mst2) +
        27000*pow2(Mgl)*pow2(Mst1)*pow2(Sbeta)*pow4(s2t)*pow6(Msq)*pow6(Mst2) +
        900*pow2(Mst1)*pow2(log(pow2(Mst2)/pow2(Mst1)))*pow4(Msq)*(2*Dmglst1*
        Mgl*Mst1*(8*Mt*s2t*Sbeta*(-153*Cbeta*Mt*MuSUSY*s2t + 968*Sbeta*pow2(Mt)
        + 32*Sbeta*pow2(Mst2)*pow2(s2t))*pow4(Mst1) + 12*Mt*s2t*(-8*pow2(Mt) +
        19*pow2(Mst2)*pow2(s2t))*pow2(Sbeta)*pow4(Mst2) - 4*Mt*pow2(Mst1)*(114*
        Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*pow2(s2t) + 4*s2t*pow2(Mt)*(95*pow2(
        MuSUSY)*(-1 + pow2(Sbeta)) - 231*pow2(Mst2)*pow2(Sbeta)) + 900*Cbeta*
        MuSUSY*Sbeta*pow3(Mt) + 19*pow2(Sbeta)*pow3(s2t)*pow4(Mst2)) + Mst1*
        pow2(Mst2)*pow2(Sbeta)*(64*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 528*pow4(Mt)
        + 155*pow4(Mst2)*pow4(s2t)) - pow3(Mst1)*(8*pow2(Mt)*pow2(s2t)*(283*
        pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 192*pow2(Mst2)*pow2(Sbeta)) - 2944*
        Cbeta*MuSUSY*s2t*Sbeta*pow3(Mt) + 512*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*
        pow3(s2t) + 384*pow2(Sbeta)*pow4(Mt) + 27*pow2(Sbeta)*pow4(Mst2)*pow4(
        s2t))) + Mst1*pow2(Dmglst1)*(16*Mt*s2t*Sbeta*(9*Cbeta*Mt*MuSUSY*s2t +
        2855*Sbeta*pow2(Mt) + 32*Sbeta*pow2(Mst2)*pow2(s2t))*pow4(Mst1) + 24*
        Mt*s2t*(-8*pow2(Mt) + 19*pow2(Mst2)*pow2(s2t))*pow2(Sbeta)*pow4(Mst2) -
        8*Mt*pow2(Mst1)*(114*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*pow2(s2t) + 2*
        s2t*pow2(Mt)*(190*pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 933*pow2(Mst2)*
        pow2(Sbeta)) + 900*Cbeta*MuSUSY*Sbeta*pow3(Mt) + 127*pow2(Sbeta)*pow3(
        s2t)*pow4(Mst2)) + 3*Mst1*pow2(Mst2)*pow2(Sbeta)*(64*pow2(Mst2)*pow2(
        Mt)*pow2(s2t) + 528*pow4(Mt) + 155*pow4(Mst2)*pow4(s2t)) - 3*pow3(Mst1)
        *(8*pow2(Mt)*pow2(s2t)*(283*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 384*pow2(
        Mst2)*pow2(Sbeta)) - 6016*Cbeta*MuSUSY*s2t*Sbeta*pow3(Mt) + 512*Cbeta*
        Mt*MuSUSY*Sbeta*pow2(Mst2)*pow3(s2t) + 944*pow2(Sbeta)*pow4(Mt) + 27*
        pow2(Sbeta)*pow4(Mst2)*pow4(s2t))) + pow2(Mgl)*(24*Mst1*Mt*s2t*(-8*
        pow2(Mt) + 19*pow2(Mst2)*pow2(s2t))*pow2(Sbeta)*pow4(Mst2) - 8*Mt*pow3(
        Mst1)*(330*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*pow2(s2t) + 4*s2t*pow2(Mt)*
        (167*pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 74*pow2(Mst2)*pow2(Sbeta)) +
        272*Cbeta*MuSUSY*Sbeta*pow3(Mt) - 53*pow2(Sbeta)*pow3(s2t)*pow4(Mst2))
        + 2*pow2(Sbeta)*pow4(Mst2)*(-246*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 106*
        pow4(Mt) + 41*pow4(Mst2)*pow4(s2t)) + pow4(Mst1)*(-24*pow2(Mt)*pow2(
        s2t)*(231*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 16*pow2(Mst2)*pow2(Sbeta))
        + 1296*Cbeta*MuSUSY*s2t*Sbeta*pow3(Mt) - 1168*Cbeta*Mt*MuSUSY*Sbeta*
        pow2(Mst2)*pow3(s2t) + 48*pow2(Sbeta)*pow4(Mt) - 27*pow2(Sbeta)*pow4(
        Mst2)*pow4(s2t)) + 16*Mt*s2t*Sbeta*(-261*Cbeta*Mt*MuSUSY*s2t + 190*
        Sbeta*pow2(Mt) + 32*Sbeta*pow2(Mst2)*pow2(s2t))*pow5(Mst1) + pow2(Mst1)
        *(4*pow2(Mst2)*pow2(Mt)*pow2(s2t)*(-802*pow2(MuSUSY)*(-1 + pow2(Sbeta))
        + 57*pow2(Mst2)*pow2(Sbeta)) + 528*Cbeta*MuSUSY*s2t*Sbeta*pow2(Mst2)*
        pow3(Mt) - 1276*Cbeta*Mt*MuSUSY*Sbeta*pow3(s2t)*pow4(Mst2) - 32*(8*
        pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 17*pow2(Mst2)*pow2(Sbeta))*pow4(Mt) +
        237*pow2(Sbeta)*pow4(s2t)*pow6(Mst2)))) - 30*log(pow2(Mst2)/pow2(Mst1))
        *pow2(Mst1)*(2*pow2(Dmglst1)*(-600*Mst1*Sbeta*pow2(Mst2)*pow3(Mt)*(-2*
        Cbeta*Mt*MuSUSY*pow2(Mst1)*(58*pow2(Mst1) + 7*pow2(Mst2)) + Sbeta*pow2(
        Mst2)*(9*Mst1*Mt*pow2(Mst2) + 24*s2t*pow2(Mst1)*pow2(Mst2) + 25*Mt*
        pow3(Mst1) + 90*s2t*pow4(Mst1) + 2*s2t*pow4(Mst2))) + 1200*Mst1*Sbeta*
        pow2(Msq)*pow3(Mt)*(4*Cbeta*Mt*MuSUSY*pow2(Mst1)*(58*pow2(Mst1) + 7*
        pow2(Mst2)) - Sbeta*pow2(Mst2)*(9*Mst1*Mt*pow2(Mst2) + 48*s2t*pow2(
        Mst1)*pow2(Mst2) + 180*s2t*pow4(Mst1) + 4*s2t*pow4(Mst2))) + pow4(Msq)*
        (80*Mst1*Mt*s2t*(215*pow2(Mt) + 396*pow2(Mst2)*pow2(s2t))*pow2(Sbeta)*
        pow4(Mst2) - 80*Mt*pow3(Mst1)*(1710*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*
        pow2(s2t) + 3*s2t*pow2(Mt)*(1288*pow2(MuSUSY)*(-1 + pow2(Sbeta)) -
        2959*pow2(Mst2)*pow2(Sbeta)) + 4438*Cbeta*MuSUSY*Sbeta*pow3(Mt) + 414*
        pow2(Sbeta)*pow3(s2t)*pow4(Mst2)) + 8*pow2(Sbeta)*pow4(Mst2)*(-720*
        pow2(Mst2)*pow2(Mt)*pow2(s2t) + 211*pow4(Mt) + 180*pow4(Mst2)*pow4(s2t)
        ) + 5*pow2(Mst1)*pow2(Mst2)*pow2(Sbeta)*(20568*pow2(Mst2)*pow2(Mt)*
        pow2(s2t) + 20216*pow4(Mt) + 2349*pow4(Mst2)*pow4(s2t)) - 5*pow4(Mst1)*
        (72*pow2(Mt)*pow2(s2t)*(491*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 125*pow2(
        Mst2)*pow2(Sbeta)) + 20832*Cbeta*MuSUSY*s2t*Sbeta*pow3(Mt) + 5976*
        Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*pow3(s2t) + 23624*pow2(Sbeta)*pow4(Mt)
        + 1143*pow2(Sbeta)*pow4(Mst2)*pow4(s2t)) + 40*Mt*s2t*Sbeta*(5859*Cbeta*
        Mt*MuSUSY*s2t + 75062*Sbeta*pow2(Mt) - 1917*Sbeta*pow2(Mst2)*pow2(s2t))
        *pow5(Mst1))) + 20*Dmglst1*Mgl*(-120*Mst1*Sbeta*pow2(Mst2)*pow3(Mt)*(-(
        Cbeta*Mt*MuSUSY*pow2(Mst1)*(22*pow2(Mst1) + 7*pow2(Mst2))) + Sbeta*
        pow2(Mst2)*(3*Mst1*Mt*pow2(Mst2) + 6*s2t*pow2(Mst1)*pow2(Mst2) + 5*Mt*
        pow3(Mst1) + 15*s2t*pow4(Mst1) + s2t*pow4(Mst2))) + 240*Mst1*Sbeta*
        pow2(Msq)*pow3(Mt)*(2*Cbeta*Mt*MuSUSY*pow2(Mst1)*(22*pow2(Mst1) + 7*
        pow2(Mst2)) - Sbeta*pow2(Mst2)*(3*Mst1*Mt*pow2(Mst2) + 12*s2t*pow2(
        Mst1)*pow2(Mst2) + 30*s2t*pow4(Mst1) + 2*s2t*pow4(Mst2))) + pow4(Msq)*(
        16*Mst1*Mt*s2t*(103*pow2(Mt) + 198*pow2(Mst2)*pow2(s2t))*pow2(Sbeta)*
        pow4(Mst2) - 16*Mt*pow3(Mst1)*(855*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*
        pow2(s2t) + s2t*pow2(Mt)*(1932*pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 2209*
        pow2(Mst2)*pow2(Sbeta)) + 2219*Cbeta*MuSUSY*Sbeta*pow3(Mt) - 87*pow2(
        Sbeta)*pow3(s2t)*pow4(Mst2)) + 16*pow2(Sbeta)*pow4(Mst2)*(-24*pow2(
        Mst2)*pow2(Mt)*pow2(s2t) + 11*pow4(Mt) + 6*pow4(Mst2)*pow4(s2t)) + 3*
        pow2(Mst1)*pow2(Mst2)*pow2(Sbeta)*(2456*pow2(Mst2)*pow2(Mt)*pow2(s2t) +
        2232*pow4(Mt) + 261*pow4(Mst2)*pow4(s2t)) - pow4(Mst1)*(24*pow2(Mt)*
        pow2(s2t)*(401*pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 263*pow2(Mst2)*pow2(
        Sbeta)) + 26592*Cbeta*MuSUSY*s2t*Sbeta*pow3(Mt) + 912*Cbeta*Mt*MuSUSY*
        Sbeta*pow2(Mst2)*pow3(s2t) + 3200*pow2(Sbeta)*pow4(Mt) + 651*pow2(
        Sbeta)*pow4(Mst2)*pow4(s2t)) + 4*Mt*s2t*Sbeta*(-1827*Cbeta*Mt*MuSUSY*
        s2t + 25312*Sbeta*pow2(Mt) - 531*Sbeta*pow2(Mst2)*pow2(s2t))*pow5(Mst1)
        )) + 15*pow2(Mgl)*(-160*Sbeta*pow2(Msq)*pow3(Mt)*(-12*Cbeta*Mt*MuSUSY*(
        2*pow2(Mst1) + pow2(Mst2))*pow3(Mst1) + Sbeta*pow2(Mst2)*(3*Mt*pow2(
        Mst1)*pow2(Mst2) + 8*s2t*pow2(Mst2)*pow3(Mst1) + 2*Mt*pow4(Mst2) + 4*
        Mst1*s2t*pow4(Mst2) + 12*s2t*pow5(Mst1))) - 40*Sbeta*pow2(Mst2)*pow3(
        Mt)*(-12*Cbeta*Mt*MuSUSY*(2*pow2(Mst1) + pow2(Mst2))*pow3(Mst1) +
        Sbeta*pow2(Mst2)*(6*Mt*pow2(Mst1)*pow2(Mst2) + 8*s2t*pow2(Mst2)*pow3(
        Mst1) + 5*Mt*pow4(Mst1) + 3*Mt*pow4(Mst2) + 4*Mst1*s2t*pow4(Mst2) + 12*
        s2t*pow5(Mst1))) + 120*pow2(s2t)*pow6(Msq)*(pow2(Mst1)*(8*pow2(Mt)*
        pow2(MuSUSY)*(-1 + pow2(Sbeta)) + pow2(s2t)*pow2(Sbeta)*pow4(Mst2)) -
        pow2(s2t)*pow2(Sbeta)*pow6(Mst2)) + pow4(Msq)*(128*Mst1*Mt*s2t*(17*
        pow2(Mt) + 33*pow2(Mst2)*pow2(s2t))*pow2(Sbeta)*pow4(Mst2) - 64*Mt*
        pow3(Mst1)*(345*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mst2)*pow2(s2t) + 4*s2t*
        pow2(Mt)*(181*pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 61*pow2(Mst2)*pow2(
        Sbeta)) + 267*Cbeta*MuSUSY*Sbeta*pow3(Mt) - 49*pow2(Sbeta)*pow3(s2t)*
        pow4(Mst2)) + 10*pow2(Sbeta)*pow4(Mst2)*(-160*pow2(Mst2)*pow2(Mt)*pow2(
        s2t) + 72*pow4(Mt) + 63*pow4(Mst2)*pow4(s2t)) + pow4(Mst1)*(-8*pow2(Mt)
        *pow2(s2t)*(1471*pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 836*pow2(Mst2)*pow2(
        Sbeta)) - 23456*Cbeta*MuSUSY*s2t*Sbeta*pow3(Mt) - 68*Cbeta*Mt*MuSUSY*
        Sbeta*pow2(Mst2)*pow3(s2t) + 296*pow2(Sbeta)*pow4(Mt) - 807*pow2(Sbeta)
        *pow4(Mst2)*pow4(s2t)) + 16*Mt*s2t*Sbeta*(-1317*Cbeta*Mt*MuSUSY*s2t +
        1724*Sbeta*pow2(Mt) - 21*Sbeta*pow2(Mst2)*pow2(s2t))*pow5(Mst1) - 2*
        pow2(Mst1)*(-8*pow2(Mst2)*pow2(Mt)*pow2(s2t)*(-727*pow2(MuSUSY)*(-1 +
        pow2(Sbeta)) + 415*pow2(Mst2)*pow2(Sbeta)) + 5040*Cbeta*MuSUSY*s2t*
        Sbeta*pow2(Mst2)*pow3(Mt) + 1648*Cbeta*Mt*MuSUSY*Sbeta*pow3(s2t)*pow4(
        Mst2) + 8*(320*pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 271*pow2(Mst2)*pow2(
        Sbeta))*pow4(Mt) - 97*pow2(Sbeta)*pow4(s2t)*pow6(Mst2))))) - 3168000*
        Dmglst1*Mgl*s2t*pow2(Msq)*pow2(Mst2)*pow2(Sbeta)*pow3(Mt)*pow7(Mst1) -
        13824000*s2t*pow2(Dmglst1)*pow2(Msq)*pow2(Mst2)*pow2(Sbeta)*pow3(Mt)*
        pow7(Mst1) - 288000*s2t*pow2(Mgl)*pow2(Msq)*pow2(Mst2)*pow2(Sbeta)*
        pow3(Mt)*pow7(Mst1) - 2754000*Cbeta*Dmglst1*Mgl*MuSUSY*Sbeta*pow2(Mt)*
        pow2(s2t)*pow4(Msq)*pow7(Mst1) + 3196800*Cbeta*Dmglst1*Mgl*MuSUSY*
        Sbeta*z2*pow2(Mt)*pow2(s2t)*pow4(Msq)*pow7(Mst1) - 1933200*Cbeta*
        MuSUSY*Sbeta*pow2(Dmglst1)*pow2(Mt)*pow2(s2t)*pow4(Msq)*pow7(Mst1) +
        7862400*Cbeta*MuSUSY*Sbeta*z2*pow2(Dmglst1)*pow2(Mt)*pow2(s2t)*pow4(
        Msq)*pow7(Mst1) + 97200*Cbeta*MuSUSY*Sbeta*pow2(Mgl)*pow2(Mt)*pow2(s2t)
        *pow4(Msq)*pow7(Mst1) + 86400*Cbeta*MuSUSY*Sbeta*z2*pow2(Mgl)*pow2(Mt)*
        pow2(s2t)*pow4(Msq)*pow7(Mst1) + 8769600*Dmglst1*Mgl*s2t*pow2(Sbeta)*
        pow3(Mt)*pow4(Msq)*pow7(Mst1) + 19756800*Dmglst1*Mgl*s2t*z2*pow2(Sbeta)
        *pow3(Mt)*pow4(Msq)*pow7(Mst1) + 69499200*s2t*pow2(Dmglst1)*pow2(Sbeta)
        *pow3(Mt)*pow4(Msq)*pow7(Mst1) + 58320000*s2t*z2*pow2(Dmglst1)*pow2(
        Sbeta)*pow3(Mt)*pow4(Msq)*pow7(Mst1) - 2692800*s2t*pow2(Mgl)*pow2(
        Sbeta)*pow3(Mt)*pow4(Msq)*pow7(Mst1) + 4032000*s2t*z2*pow2(Mgl)*pow2(
        Sbeta)*pow3(Mt)*pow4(Msq)*pow7(Mst1) + 1018800*Dmglst1*Mgl*Mt*pow2(
        Mst2)*pow2(Sbeta)*pow3(s2t)*pow4(Msq)*pow7(Mst1) + 2919600*Mt*pow2(
        Dmglst1)*pow2(Mst2)*pow2(Sbeta)*pow3(s2t)*pow4(Msq)*pow7(Mst1) - 18000*
        Mt*pow2(Mgl)*pow2(Mst2)*pow2(Sbeta)*pow3(s2t)*pow4(Msq)*pow7(Mst1) -
        954000*Dmglst1*Mgl*s2t*pow2(Sbeta)*pow3(Mt)*pow4(Mst2)*pow7(Mst1) -
        4158000*s2t*pow2(Dmglst1)*pow2(Sbeta)*pow3(Mt)*pow4(Mst2)*pow7(Mst1) -
        18000*s2t*pow2(Mgl)*pow2(Sbeta)*pow3(Mt)*pow4(Mst2)*pow7(Mst1) +
        6000000*Cbeta*Dmglst1*Mgl*MuSUSY*Sbeta*pow2(Msq)*pow4(Mt)*pow7(Mst1) +
        20688000*Cbeta*MuSUSY*Sbeta*pow2(Dmglst1)*pow2(Msq)*pow4(Mt)*pow7(Mst1)
        + 1008000*Cbeta*MuSUSY*Sbeta*pow2(Mgl)*pow2(Msq)*pow4(Mt)*pow7(Mst1) +
        1500000*Cbeta*Dmglst1*Mgl*MuSUSY*Sbeta*pow2(Mst2)*pow4(Mt)*pow7(Mst1) +
        5172000*Cbeta*MuSUSY*Sbeta*pow2(Dmglst1)*pow2(Mst2)*pow4(Mt)*pow7(Mst1)
        + 252000*Cbeta*MuSUSY*Sbeta*pow2(Mgl)*pow2(Mst2)*pow4(Mt)*pow7(Mst1) -
        78000*Dmglst1*Mgl*s2t*pow2(Sbeta)*pow3(Mst1)*pow3(Mt)*pow8(Mst2) -
        78000*s2t*pow2(Dmglst1)*pow2(Sbeta)*pow3(Mst1)*pow3(Mt)*pow8(Mst2) -
        78000*s2t*pow2(Mgl)*pow2(Sbeta)*pow3(Mst1)*pow3(Mt)*pow8(Mst2) + 21600*
        pow2(Mgl)*pow2(Mt)*pow2(s2t)*pow2(Sbeta)*pow4(Msq)*pow8(Mst2) + 10800*
        Cbeta*Mt*MuSUSY*Sbeta*pow2(Mgl)*pow3(s2t)*pow4(Msq)*pow8(Mst2) - 40500*
        pow2(Mgl)*pow2(Mst1)*pow2(Sbeta)*pow4(Mt)*pow8(Mst2) + 16200*Dmglst1*
        Mgl*pow2(Mst1)*pow2(Sbeta)*pow4(Msq)*pow4(s2t)*pow8(Mst2) + 13800*pow2(
        Dmglst1)*pow2(Mst1)*pow2(Sbeta)*pow4(Msq)*pow4(s2t)*pow8(Mst2) +
        307800*pow2(Mgl)*pow2(Mst1)*pow2(Sbeta)*pow4(Msq)*pow4(s2t)*pow8(Mst2)
        - 36900*z2*pow2(Mgl)*pow2(Mst1)*pow2(Sbeta)*pow4(Msq)*pow4(s2t)*pow8(
        Mst2) - 27000*pow2(Mgl)*pow2(Sbeta)*pow4(s2t)*pow6(Msq)*pow8(Mst2))/(
        16200.*pow2(Mgl)*pow2(Mst1)*pow4(Msq)*pow4(Mst2)))/pow4(Mt)/pow2(Sbeta)*
        12.; 

   return result;
}

/**
 * 	@return returns the susy log^2 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double himalaya::H5g1::calc_coef_at_as2_no_sm_logs_log2(){

   const double result =
      ((14760*pow2(Mgl)*pow2(Mst1)*pow2(Mst2)*pow2(Mt)*pow2(MuSUSY)*pow2(s2t) -
        14760*pow2(Mgl)*pow2(Mst1)*pow2(Mst2)*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*
        pow2(Sbeta) - 34560*Cbeta*Dmglst1*Mgl*MuSUSY*Sbeta*pow2(Mst2)*pow2(Mt)*
        pow2(s2t)*pow3(Mst1) - 34560*Cbeta*MuSUSY*Sbeta*pow2(Dmglst1)*pow2(
        Mst2)*pow2(Mt)*pow2(s2t)*pow3(Mst1) - 11520*Cbeta*MuSUSY*Sbeta*pow2(
        Mgl)*pow2(Mst2)*pow2(Mt)*pow2(s2t)*pow3(Mst1) - 67920*Cbeta*MuSUSY*s2t*
        Sbeta*pow2(Mgl)*pow2(Mst1)*pow2(Mst2)*pow3(Mt) + 170880*Dmglst1*Mgl*
        s2t*pow2(MuSUSY)*pow3(Mst1)*pow3(Mt) + 170880*s2t*pow2(Dmglst1)*pow2(
        MuSUSY)*pow3(Mst1)*pow3(Mt) + 109440*s2t*pow2(Mgl)*pow2(MuSUSY)*pow3(
        Mst1)*pow3(Mt) + 136320*Dmglst1*Mgl*s2t*pow2(Mst2)*pow2(Sbeta)*pow3(
        Mst1)*pow3(Mt) + 399360*s2t*pow2(Dmglst1)*pow2(Mst2)*pow2(Sbeta)*pow3(
        Mst1)*pow3(Mt) + 17280*s2t*pow2(Mgl)*pow2(Mst2)*pow2(Sbeta)*pow3(Mst1)*
        pow3(Mt) - 170880*Dmglst1*Mgl*s2t*pow2(MuSUSY)*pow2(Sbeta)*pow3(Mst1)*
        pow3(Mt) - 170880*s2t*pow2(Dmglst1)*pow2(MuSUSY)*pow2(Sbeta)*pow3(Mst1)
        *pow3(Mt) - 109440*s2t*pow2(Mgl)*pow2(MuSUSY)*pow2(Sbeta)*pow3(Mst1)*
        pow3(Mt) - 7680*Dmglst1*Mgl*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*pow4(Mst1)
        - 19200*pow2(Dmglst1)*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*pow4(Mst1) -
        1920*pow2(Mgl)*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*pow4(Mst1) - 29520*
        Dmglst1*Mgl*pow2(Mst2)*pow2(Mt)*pow2(s2t)*pow2(Sbeta)*pow4(Mst1) -
        51960*pow2(Dmglst1)*pow2(Mst2)*pow2(Mt)*pow2(s2t)*pow2(Sbeta)*pow4(
        Mst1) - 12840*pow2(Mgl)*pow2(Mst2)*pow2(Mt)*pow2(s2t)*pow2(Sbeta)*pow4(
        Mst1) + 7680*Dmglst1*Mgl*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*pow2(Sbeta)*
        pow4(Mst1) + 19200*pow2(Dmglst1)*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*pow2(
        Sbeta)*pow4(Mst1) + 1920*pow2(Mgl)*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*
        pow2(Sbeta)*pow4(Mst1) - 76800*Cbeta*Dmglst1*Mgl*MuSUSY*s2t*Sbeta*pow3(
        Mt)*pow4(Mst1) - 99840*Cbeta*MuSUSY*s2t*Sbeta*pow2(Dmglst1)*pow3(Mt)*
        pow4(Mst1) - 42240*Cbeta*MuSUSY*s2t*Sbeta*pow2(Mgl)*pow3(Mt)*pow4(Mst1)
        + 22440*Cbeta*Dmglst1*Mgl*Mt*MuSUSY*Sbeta*pow2(Mst2)*pow3(s2t)*pow4(
        Mst1) + 45180*Cbeta*Mt*MuSUSY*Sbeta*pow2(Dmglst1)*pow2(Mst2)*pow3(s2t)*
        pow4(Mst1) + 8340*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mgl)*pow2(Mst2)*pow3(s2t)*
        pow4(Mst1) + 82080*Dmglst1*Mgl*pow2(Mst1)*pow2(Mt)*pow2(s2t)*pow2(
        Sbeta)*pow4(Mst2) + 115440*pow2(Dmglst1)*pow2(Mst1)*pow2(Mt)*pow2(s2t)*
        pow2(Sbeta)*pow4(Mst2) + 42960*pow2(Mgl)*pow2(Mst1)*pow2(Mt)*pow2(s2t)*
        pow2(Sbeta)*pow4(Mst2) + 36160*Dmglst1*Mgl*Mst1*s2t*pow2(Sbeta)*pow3(
        Mt)*pow4(Mst2) + 24000*Mst1*s2t*pow2(Dmglst1)*pow2(Sbeta)*pow3(Mt)*
        pow4(Mst2) + 25920*Mst1*s2t*pow2(Mgl)*pow2(Sbeta)*pow3(Mt)*pow4(Mst2) +
        6960*Cbeta*Mt*MuSUSY*Sbeta*pow2(Mgl)*pow2(Mst1)*pow3(s2t)*pow4(Mst2) -
        19680*Dmglst1*Mgl*Mt*pow2(Sbeta)*pow3(Mst1)*pow3(s2t)*pow4(Mst2) -
        19680*Mt*pow2(Dmglst1)*pow2(Sbeta)*pow3(Mst1)*pow3(s2t)*pow4(Mst2) -
        19680*Mt*pow2(Mgl)*pow2(Sbeta)*pow3(Mst1)*pow3(s2t)*pow4(Mst2) + 30720*
        pow2(Mgl)*pow2(Mst1)*pow2(MuSUSY)*pow4(Mt) + 86880*Dmglst1*Mgl*pow2(
        Mst1)*pow2(Mst2)*pow2(Sbeta)*pow4(Mt) + 171280*pow2(Dmglst1)*pow2(Mst1)
        *pow2(Mst2)*pow2(Sbeta)*pow4(Mt) + 28080*pow2(Mgl)*pow2(Mst1)*pow2(
        Mst2)*pow2(Sbeta)*pow4(Mt) - 30720*pow2(Mgl)*pow2(Mst1)*pow2(MuSUSY)*
        pow2(Sbeta)*pow4(Mt) - 172480*Cbeta*Dmglst1*Mgl*MuSUSY*Sbeta*pow3(Mst1)
        *pow4(Mt) - 172480*Cbeta*MuSUSY*Sbeta*pow2(Dmglst1)*pow3(Mst1)*pow4(Mt)
        - 43200*Cbeta*MuSUSY*Sbeta*pow2(Mgl)*pow3(Mst1)*pow4(Mt) - 18240*
        Dmglst1*Mgl*pow2(Sbeta)*pow4(Mst1)*pow4(Mt) - 75680*pow2(Dmglst1)*pow2(
        Sbeta)*pow4(Mst1)*pow4(Mt) + 480*pow2(Mgl)*pow2(Sbeta)*pow4(Mst1)*pow4(
        Mt) + 29600*Dmglst1*Mgl*pow2(Sbeta)*pow4(Mst2)*pow4(Mt) + 19856*pow2(
        Dmglst1)*pow2(Sbeta)*pow4(Mst2)*pow4(Mt) + 29440*pow2(Mgl)*pow2(Sbeta)*
        pow4(Mst2)*pow4(Mt) + 7200*log(pow2(Msq)/pow2(Mst1))*pow2(Mgl)*pow2(
        Sbeta)*pow4(Mst2)*pow4(Mt) - 8490*Dmglst1*Mgl*pow2(Sbeta)*pow4(Mst1)*
        pow4(Mst2)*pow4(s2t) - 18495*pow2(Dmglst1)*pow2(Sbeta)*pow4(Mst1)*pow4(
        Mst2)*pow4(s2t) - 345*pow2(Mgl)*pow2(Sbeta)*pow4(Mst1)*pow4(Mst2)*pow4(
        s2t) + 97920*Dmglst1*Mgl*s2t*pow2(Sbeta)*pow3(Mt)*pow5(Mst1) + 514560*
        s2t*pow2(Dmglst1)*pow2(Sbeta)*pow3(Mt)*pow5(Mst1) + 1920*s2t*pow2(Mgl)*
        pow2(Sbeta)*pow3(Mt)*pow5(Mst1) - 11520*Dmglst1*Mgl*Mt*pow2(Mst2)*pow2(
        Sbeta)*pow3(s2t)*pow5(Mst1) - 23040*Mt*pow2(Dmglst1)*pow2(Mst2)*pow2(
        Sbeta)*pow3(s2t)*pow5(Mst1) - 3840*Mt*pow2(Mgl)*pow2(Mst2)*pow2(Sbeta)*
        pow3(s2t)*pow5(Mst1) - 14160*Dmglst1*Mgl*pow2(Mt)*pow2(s2t)*pow2(Sbeta)
        *pow6(Mst2) - 13560*pow2(Dmglst1)*pow2(Mt)*pow2(s2t)*pow2(Sbeta)*pow6(
        Mst2) - 9000*pow2(Mgl)*pow2(Mt)*pow2(s2t)*pow2(Sbeta)*pow6(Mst2) +
        31200*Dmglst1*Mgl*Mst1*Mt*pow2(Sbeta)*pow3(s2t)*pow6(Mst2) + 42720*
        Mst1*Mt*pow2(Dmglst1)*pow2(Sbeta)*pow3(s2t)*pow6(Mst2) + 23520*Mst1*Mt*
        pow2(Mgl)*pow2(Sbeta)*pow3(s2t)*pow6(Mst2) + 1110*Dmglst1*Mgl*pow2(
        Mst1)*pow2(Sbeta)*pow4(s2t)*pow6(Mst2) + 5505*pow2(Dmglst1)*pow2(Mst1)*
        pow2(Sbeta)*pow4(s2t)*pow6(Mst2) - 5325*pow2(Mgl)*pow2(Mst1)*pow2(
        Sbeta)*pow4(s2t)*pow6(Mst2) - 30*log(pow2(Mst2)/pow2(Mst1))*(Mst1*pow2(
        Dmglst1)*(48*s2t*Sbeta*(-137*Cbeta*MuSUSY*s2t + 460*Mt*Sbeta)*pow2(Mt)*
        pow4(Mst1) + 8*Mt*s2t*(-60*pow2(Mt) + 41*pow2(Mst2)*pow2(s2t))*pow2(
        Sbeta)*pow4(Mst2) - 8*Mt*pow2(Mst1)*(534*Cbeta*Mt*MuSUSY*Sbeta*pow2(
        Mst2)*pow2(s2t) + 12*s2t*pow2(Mt)*(73*pow2(MuSUSY)*(-1 + pow2(Sbeta)) -
        88*pow2(Mst2)*pow2(Sbeta)) + 468*Cbeta*MuSUSY*Sbeta*pow3(Mt) - 233*
        pow2(Sbeta)*pow3(s2t)*pow4(Mst2)) + 3*Mst1*pow2(Mst2)*pow2(Sbeta)*(384*
        pow2(Mst2)*pow2(Mt)*pow2(s2t) + 512*pow4(Mt) + 91*pow4(Mst2)*pow4(s2t))
        + pow3(Mst1)*(-24*pow2(Mt)*pow2(s2t)*(91*pow2(MuSUSY)*(-1 + pow2(Sbeta)
        ) - 96*pow2(Mst2)*pow2(Sbeta)) - 6912*Cbeta*MuSUSY*s2t*Sbeta*pow3(Mt) +
        3296*pow2(Sbeta)*pow4(Mt) - 273*pow2(Sbeta)*pow4(Mst2)*pow4(s2t))) + 2*
        Dmglst1*Mgl*Mst1*(8*s2t*Sbeta*(-267*Cbeta*MuSUSY*s2t + 460*Mt*Sbeta)*
        pow2(Mt)*pow4(Mst1) + 4*Mt*s2t*(-60*pow2(Mt) + 41*pow2(Mst2)*pow2(s2t))
        *pow2(Sbeta)*pow4(Mst2) - 4*Mt*pow2(Mst1)*(534*Cbeta*Mt*MuSUSY*Sbeta*
        pow2(Mst2)*pow2(s2t) + 12*s2t*pow2(Mt)*(73*pow2(MuSUSY)*(-1 + pow2(
        Sbeta)) - 44*pow2(Mst2)*pow2(Sbeta)) + 468*Cbeta*MuSUSY*Sbeta*pow3(Mt)
        - 137*pow2(Sbeta)*pow3(s2t)*pow4(Mst2)) + Mst1*pow2(Mst2)*pow2(Sbeta)*(
        384*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 512*pow4(Mt) + 91*pow4(Mst2)*pow4(
        s2t)) + pow3(Mst1)*(-8*pow2(Mt)*pow2(s2t)*(91*pow2(MuSUSY)*(-1 + pow2(
        Sbeta)) - 96*pow2(Mst2)*pow2(Sbeta)) - 2304*Cbeta*MuSUSY*s2t*Sbeta*
        pow3(Mt) + 864*pow2(Sbeta)*pow4(Mt) - 91*pow2(Sbeta)*pow4(Mst2)*pow4(
        s2t))) + pow2(Mgl)*(8*Mst1*Mt*s2t*(-60*pow2(Mt) + 41*pow2(Mst2)*pow2(
        s2t))*pow2(Sbeta)*pow4(Mst2) - 8*Mt*pow3(Mst1)*(342*Cbeta*Mt*MuSUSY*
        Sbeta*pow2(Mst2)*pow2(s2t) + 4*s2t*pow2(Mt)*(155*pow2(MuSUSY)*(-1 +
        pow2(Sbeta)) - 44*pow2(Mst2)*pow2(Sbeta)) + 116*Cbeta*MuSUSY*Sbeta*
        pow3(Mt) - 73*pow2(Sbeta)*pow3(s2t)*pow4(Mst2)) + pow2(Sbeta)*pow4(
        Mst2)*(-328*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 408*pow4(Mt) + 41*pow4(
        Mst2)*pow4(s2t)) + 4*pow4(Mst1)*(2*pow2(Mt)*pow2(s2t)*(-173*pow2(
        MuSUSY)*(-1 + pow2(Sbeta)) + 96*pow2(Mst2)*pow2(Sbeta)) - 576*Cbeta*
        MuSUSY*s2t*Sbeta*pow3(Mt) + 172*pow2(Sbeta)*pow4(Mt) - 33*pow2(Sbeta)*
        pow4(Mst2)*pow4(s2t)) + 16*s2t*Sbeta*(-171*Cbeta*MuSUSY*s2t + 92*Mt*
        Sbeta)*pow2(Mt)*pow5(Mst1) + pow2(Mst1)*(8*pow2(Mst2)*pow2(Mt)*pow2(
        s2t)*(-173*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 89*pow2(Mst2)*pow2(Sbeta))
        - 768*Cbeta*MuSUSY*s2t*Sbeta*pow2(Mst2)*pow3(Mt) - 528*Cbeta*Mt*MuSUSY*
        Sbeta*pow3(s2t)*pow4(Mst2) + 512*(-(pow2(MuSUSY)*(-1 + pow2(Sbeta))) +
        pow2(Mst2)*pow2(Sbeta))*pow4(Mt) + 91*pow2(Sbeta)*pow4(s2t)*pow6(Mst2))
        )) + 1770*Dmglst1*Mgl*pow2(Sbeta)*pow4(s2t)*pow8(Mst2) + 1695*pow2(
        Dmglst1)*pow2(Sbeta)*pow4(s2t)*pow8(Mst2) + 3585*pow2(Mgl)*pow2(Sbeta)*
        pow4(s2t)*pow8(Mst2))/(540.*pow2(Mgl)*pow4(Mst2)))/pow4(Mt)/pow2(Sbeta)*
        12.; 

   return result;
}

/**
 * 	@return returns the susy log^3 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double himalaya::H5g1::calc_coef_at_as2_no_sm_logs_log3(){

   const double result =
      ((-224*pow2(Sbeta)*pow4(Mt))/9.)/pow4(Mt)/pow2(Sbeta)*12.; 

   return result;
}



