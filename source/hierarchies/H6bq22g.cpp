// ====================================================================
// This file is part of Himalaya.
//
// Himalaya is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================

#include "H6bq22g.hpp"
#include "Hierarchies.hpp"
#include "Constants.hpp"
#include "Utils.hpp"
#include <cmath>
#include <type_traits>

namespace himalaya{
namespace hierarchies{

/**
 * Constructor
 * @param flagMap the flagMap for the truncation of expansion variables
 * @param Al4p a double alpha_s/4/Pi
 * @param beta a double which is the mixing angle beta
 * @param Dmglst2 a double Mgl - Mst2
 * @param Dmsqst2 a double Msq - Mst2
 * @param lmMt a double log((<renormalization scale> / Mt)^2)
 * @param lmMst1 a double log((<renormalization scale> / Mst1)^2)
 * @param lmMst2 a double log((<renormalization scale> / Mst2)^2)
 * @param Mt a double top/bottom quark mass
 * @param Mst1 a double stop 1 mass
 * @param Mst2 a double stop 2 mass
 * @param Msq a double average squark mass w/o the stop quark
 * @param MuSUSY a double mu parameter
 * @param s2t a double 2 times the sine of the stop/sbottom quark mixing angle
 * @param mdrFlag an int 0 for DR and 1 for MDR scheme
 * @param oneLoopFlag an int flag to consider the one-loop expansion terms
 * @param twoLoopFlag an int flag to consider the two-loop expansion terms
 * @param threeLoopFlag an int flag to consider the three-loop expansion terms
 */
H6bq22g::H6bq22g(const std::map<unsigned int, unsigned int>& flagMap, double Al4p, double beta, double Dmglst2,
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
   this -> oneLoopFlag = oneLoopFlag;
   this -> twoLoopFlag = twoLoopFlag;
   this -> threeLoopFlag = threeLoopFlag;
   this -> Al4p = Al4p;
   // mdr flags, indicates if one wants to shift the dr stop mass to the mdr stop mass
   shiftst1 = mdrFlag;
   shiftst2 = mdrFlag;
   shiftst3 = mdrFlag;
   xDR2DRMOD = mdrFlag;
   // expansion flags
   xDmglst2 = flagMap.at(ExpansionDepth::xxDmglst2);
   xDmsqst2 = flagMap.at(ExpansionDepth::xxDmsqst2);
   xMst = flagMap.at(ExpansionDepth::xxMst);
}

/**
 * @return The diagonal (1, 1) matrix element of the Higgs mass matrix as a double for the hierarchy 'H6bq22g'
 */
double H6bq22g::getS1() const {
   return -(threeLoopFlag*pow2(Al4p)*(pow2(Mt)*pow2(s2t)*((5*Dmglst2*Dmsqst2*(-1081
        + 165*lmMst1 - 165*lmMst2)*(-1 + pow2(Sbeta))*pow2(Sbeta)*pow4(Mst1))/(
        9.*pow2(Msq)*pow3(Mst2)) + pow2(MuSUSY)*(53.385802469135804 + (40*B4)/
        9. - (4*D3)/9. + (2*DN)/9. + (1672*lmMst1)/27. + (53*pow2(lmMst1))/9. -
        lmMst2*(129.92592592592592 - 72*lmMst1 + (13*pow2(lmMst1))/3.) + (2*(-
        470 + 147*lmMst1)*pow2(lmMst2))/9. + (Dmglst2*Mst2*(4.444444444444445 -
        (4*(99 + 16*lmMst1)*lmMst2)/9. - (82*pow2(lmMst2))/3. - (40*Dmsqst2)/
        pow2(Msq)))/pow2(Mst1) + ((5*Dmsqst2*(1141 + 48*lmMst1*(7 - 3*lmMst2) -
        264*lmMst2 + 144*pow2(lmMst2)))/54. - ((Dmsqst2*(30 - 60*lmMst2) + (103
        + 186*lmMst2 + 32*lmMst1*(1 + lmMst2) + 91*pow2(lmMst2))*pow2(Msq))*
        pow2(Mst2))/(9.*pow2(Mst1)))/pow2(Msq) - (271*pow3(lmMst2))/9. - (
        Dmglst2*(4320*Dmsqst2*(4 + lmMst1 - lmMst2) + pow2(Msq)*(14267 - 432*B4
        + 576*D3 - 360*DN + 4752*lmMst1 - 1404*pow2(lmMst1) + 72*lmMst2*(-63 -
        232*lmMst1 + 16*pow2(lmMst1)) - 72*(-281 + 123*lmMst1)*pow2(lmMst2) +
        7704*pow3(lmMst2))))/(162.*Mst2*pow2(Msq)) - (pow2(Mst1)*(-(Mst2*(
        434.270658436214 + (76*B4)/9. - (2*DN)/9. + (69088*lmMst1)/405. - (
        1313*pow2(lmMst1))/27. - (4*lmMst2*(16192 - 26430*lmMst1 + 3465*pow2(
        lmMst1)))/405. + ((-5735 + 3072*lmMst1)*pow2(lmMst2))/27. + (Dmsqst2*(
        201.74098765432097 - (622*lmMst2)/15. + (2*lmMst1*(311 + 10*lmMst2))/
        15. - (22*pow2(lmMst1))/3. + 6*pow2(lmMst2)))/pow2(Msq) - (62*pow3(
        lmMst1))/27. - (2086*pow3(lmMst2))/27.)) + (2*Dmglst2*(2695042 - 40500*
        B4 + 54000*D3 - 33750*DN - 326895*lmMst1 - 324900*pow2(lmMst1) + 15*
        lmMst2*(-19607 - 129030*lmMst1 + 62550*pow2(lmMst1)) + 450*(5023 -
        5610*lmMst1)*pow2(lmMst2) + (50625*Dmsqst2*(-23 + 42*lmMst1 - 42*
        lmMst2))/pow2(Msq) + 11250*pow3(lmMst1) + 1575000*pow3(lmMst2)))/30375.
        ))/pow3(Mst2) + (16*(pow3(lmMst1) + ((1 + lmMst2)*(4*Dmglst2*lmMst2 +
        Mst2 + lmMst2*Mst2)*pow3(Mst2))/pow4(Mst1)))/9. - ((-(Mst2*(
        628.1736268201578 + (76*B4)/9. - (2*DN)/9. + (6317839*lmMst1)/396900. -
        (66307*pow2(lmMst1))/315. - lmMst2*(12.52907281431091 - (182909*lmMst1)
        /315. + (274*pow2(lmMst1))/3.) + (2*(-58301 + 37135*lmMst1)*pow2(
        lmMst2))/315. + (Dmsqst2*(237.28785508324435 + (16526*lmMst2)/3969. + (
        2*lmMst1*(-8263 + 71820*lmMst2))/3969. - (520*pow2(lmMst1))/21. - (80*
        pow2(lmMst2))/7.))/pow2(Msq) - (44*pow3(lmMst1))/9. - (1256*pow3(
        lmMst2))/9.)) + Dmglst2*(585.1892843532082 - (8*B4)/3. + (32*D3)/9. - (
        20*DN)/9. - (20109937*lmMst1)/297675. - (15886*pow2(lmMst1))/945. +
        lmMst2*(17.112243218274966 - (144628*lmMst1)/945. + (1180*pow2(lmMst1))
        /9.) + (169.85608465608465 - (2632*lmMst1)/9.)*pow2(lmMst2) - (92*pow3(
        lmMst1))/27. + (4448*pow3(lmMst2))/27.))*pow4(Mst1))/pow5(Mst2) + (2*
        Dmglst2*pow2(Msq)*(54*(344*OepS2 + 9*(15643 - 774*lmMst1 + 774*lmMst2)*
        S2)*pow2(Mst1)*pow2(Mst2) + 4*(17308*OepS2 + 27*(93919 - 12981*lmMst1 +
        12981*lmMst2)*S2)*pow4(Mst1) + 27*(56*OepS2 - 81*(-1677 + 14*lmMst1 -
        14*lmMst2)*S2)*pow4(Mst2)) + (7290*Mst2*shiftst1*(Dmsqst2*(3 - 2*
        lmMst2) + (1 - 2*lmMst2)*pow2(Msq))*((pow2(Mst1) + pow2(Mst2))*pow4(
        Mst2) - 2*(lmMst1 - lmMst2)*pow2(Mst1)*(pow2(Mst1)*pow2(Mst2) + pow4(
        Mst1) + pow4(Mst2))))/pow2(Mst1) - 3*Mst2*(pow2(Msq)*(3*(8456*OepS2 -
        81*(11243 + 2114*lmMst1 - 2114*lmMst2)*S2)*pow2(Mst1)*pow2(Mst2) + (
        52948*OepS2 - 27*(194357 + 39711*lmMst1 - 39711*lmMst2)*S2)*pow4(Mst1)
        + 27*(184*OepS2 - 81*(307 + 46*lmMst1 - 46*lmMst2)*S2)*pow4(Mst2)) +
        10*Dmsqst2*(9*(104*OepS2 + 27*(17 - 78*lmMst1 + 78*lmMst2)*S2)*pow2(
        Mst1)*pow2(Mst2) + 13*(136*OepS2 - 27*(95 + 102*lmMst1 - 102*lmMst2)*
        S2)*pow4(Mst1) + 27*(8*OepS2 - 81*(-15 + 2*lmMst1 - 2*lmMst2)*S2)*pow4(
        Mst2))))/(2187.*pow2(Msq)*pow5(Mst2)) + ((1 - 2*lmMst2)*shiftst3*(2*(1
        - 2*lmMst1 + 2*lmMst2)*pow2(Mst2)*pow4(Mst1) + 2*(1 - lmMst1 + lmMst2)*
        pow2(Mst1)*pow4(Mst2) + (2 - 6*lmMst1 + 6*lmMst2)*pow6(Mst1) + pow6(
        Mst2)))/(3.*pow2(Mst1)*pow4(Mst2)))) + (Mt*(-5145*z3*pow4(Mst1)*(3*s2t*
        pow2(Mst2)*(5*xDmsqst2*pow2(Dmsqst2)*(Mt*pow2(MuSUSY)*(48*Mst2*(3024*Mt
        - 1783*Mst2*s2t)*pow2(Mst1) + 9*(8064*Mt - 5507*Mst2*s2t)*pow3(Mst2)) +
        8*s2t*(3542*Mt*pow2(MuSUSY) + 243*Mt*pow2(Mst2)*(-1 + pow2(Sbeta))*
        pow2(Sbeta) - 810*s2t*(-1 + pow2(Sbeta))*pow2(Sbeta)*pow3(Mst2))*pow4(
        Mst1)) + Mt*s2t*pow2(MuSUSY)*(80*Dmsqst2*pow2(Msq)*(315*pow2(Mst1)*
        pow2(Mst2) + 1771*pow4(Mst1) - 621*pow4(Mst2)) + 16*pow4(Msq)*(3*(35719
        + 108*lmMst1 - 108*lmMst2)*pow2(Mst1)*pow2(Mst2) + 2*(91963 + 162*
        lmMst1 - 162*lmMst2)*pow4(Mst1) + 27*(1319 + 6*lmMst1 - 6*lmMst2)*pow4(
        Mst2)))) + 4*Dmglst2*Mt*(-3*Dmglst2*xDmglst2*pow2(MuSUSY)*(8*pow2(Mst1)
        *(-161918*Mst2*Mt*s2t + 77334*pow2(Mt) - 114571*pow2(Mst2)*pow2(s2t)) -
        3*pow2(Mst2)*(262184*Mst2*Mt*s2t - 87552*pow2(Mt) + 122917*pow2(Mst2)*
        pow2(s2t)))*pow4(Msq) + 2*Mst2*s2t*(-4860*Mst2*xDmsqst2*pow2(Dmsqst2)*(
        pow2(MuSUSY)*(8*(5*Mt - 2*Mst2*s2t)*pow2(Mst1) + 8*Mt*pow2(Mst2) - 6*
        s2t*pow3(Mst2)) + 5*Mst2*s2t*(-1 + pow2(Sbeta))*pow2(Sbeta)*pow4(Mst1))
        + s2t*(4860*Dmsqst2*pow2(Msq)*pow2(Mst2)*(2*(8*pow2(Mst1) + 3*pow2(
        Mst2))*pow2(MuSUSY) - 5*(-1 + pow2(Sbeta))*pow2(Sbeta)*pow4(Mst1)) +
        pow2(MuSUSY)*pow4(Msq)*(1052676*pow2(Mst1)*pow2(Mst2) + 1412464*pow4(
        Mst1) + 475605*pow4(Mst2)))))) - 6*s2t*xDmsqst2*pow2(Dmsqst2)*pow2(
        Mst1)*pow2(Mst2)*(s2t*(694575*Mst2*pow2(Sbeta)*(-6*(11 + 39*lmMst1 -
        39*lmMst2)*s2t*pow2(Mst2)*(-1 + pow2(Sbeta)) + (Dmglst2*(4324 - 660*
        lmMst1 + 660*lmMst2) + 3*(-167 + 5*lmMst1 - 5*lmMst2)*Mst2)*Mt*pow2(
        Sbeta))*pow6(Mst1) + Mt*(1029*Mst2*pow2(MuSUSY)*(75*(1728*Dmglst2*(4 +
        lmMst1 - lmMst2) + Mst2*(-8771 + 128*OepS2 + 72*lmMst1*(-29 + 12*lmMst2
        - 36*S2) + 57024*S2 + 72*lmMst2*(23 + 36*S2) - 864*pow2(lmMst2)))*pow2(
        Mst1)*pow2(Mst2) - 4*(4050*Dmglst2*(23 - 42*lmMst1 + 42*lmMst2) + Mst2*
        (212566 - 69615*lmMst2 - 10000*OepS2 - 1350*(947 + 150*lmMst2)*S2 + 45*
        lmMst1*(1547 - 180*lmMst2 + 4500*S2) - 4050*pow2(lmMst1) + 12150*pow2(
        lmMst2)))*pow4(Mst1) - 16200*(-12*Dmglst2 + Mst2 + 2*lmMst2*Mst2)*pow4(
        Mst2)) - (2*(593331163 - 60642400*OepS2 + 1260*lmMst2*(8263 - 974610*
        S2) + 1143733500*S2 + 1260*lmMst1*(-8263 + 71820*lmMst2 + 974610*S2) -
        61916400*pow2(lmMst1) - 28576800*pow2(lmMst2))*pow2(MuSUSY) + 694575*
        Mst2*(Dmglst2*(4324 - 660*lmMst1 + 660*lmMst2) + 3*(-167 + 5*lmMst1 -
        5*lmMst2)*Mst2)*pow2(Sbeta))*pow6(Mst1))) + 411600*Mt*pow2(MuSUSY)*(
        pow2(Mst1)*(648*Dmglst2*Mt*((8 - 15*lmMst1 + 15*lmMst2)*pow2(Mst1) + (-
        2 - 5*lmMst1 + 5*lmMst2)*pow2(Mst2)) + 9*(36*(5 + 4*lmMst1 - 4*lmMst2)*
        Mt - Mst2*s2t*(9*(-1 + 2*lmMst1 - 2*lmMst2)*(-2 + lmMst2)*shiftst1 + 4*
        T1ep))*pow3(Mst2)) + 6*Mst2*(108*(1 + 3*lmMst1 - 3*lmMst2)*Mt + Mst2*
        s2t*(-54*lmMst2*shiftst1 - 25*T1ep + 27*shiftst1*(-(lmMst1*(-2 +
        lmMst2)) + pow2(lmMst2))))*pow4(Mst1) + s2t*(-((81*(lmMst1 - lmMst2)*(-
        3 + 2*lmMst2)*shiftst1 + 442*T1ep)*pow6(Mst1)) + 81*(-2 + lmMst2)*
        shiftst1*pow6(Mst2)))) - 686*Mt*pow2(MuSUSY)*(30*pow4(Mst1)*(4*Mt*z3*
        pow2(Msq)*(19440*Dmsqst2*s2t*pow2(Mst2)*(7*Mst2*pow2(Mst1) - Dmglst2*(
        5*pow2(Mst1) + pow2(Mst2)) + 3*pow3(Mst2)) + pow2(Msq)*(36*Mst2*pow2(
        Mst1)*(-750*Dmglst2*Mt + 981*Mst2*Mt + 23728*Dmglst2*Mst2*s2t + 7318*
        s2t*pow2(Mst2)) + 27*(-664*Dmglst2*Mt + 920*Mst2*Mt + 15137*Dmglst2*
        Mst2*s2t + 5965*s2t*pow2(Mst2))*pow3(Mst2) + 2*(6759*Mt + 791162*
        Dmglst2*s2t + 149448*Mst2*s2t)*pow4(Mst1))) + 6*pow2(Mt)*pow4(Msq)*(3*
        Mst2*pow2(Mst1)*(3*Mst2*(10667 - 96*B4 + 96*D3 - 48*DN - 3072*lmMst1 -
        384*(-25 + 9*lmMst1)*lmMst2 - 384*lmMst1*lmMt + 384*(1 + lmMst2)*lmMt -
        384*(-13 + 4*lmMst1)*pow2(lmMst2) + 1536*pow3(lmMst2)) + 2*Dmglst2*(
        28405 - 288*B4 + 288*D3 - 144*DN + 5472*lmMst1 - 288*(-89 + 2*lmMst1)*
        lmMst2 + 576*(-2 + lmMst1 - lmMst2)*lmMt - 576*(-9 + 8*lmMst1)*pow2(
        lmMst2) + 4608*pow3(lmMst2))) + 432*Dmglst2*(180 - 2*B4 + 2*D3 - DN +
        144*lmMst2 - 16*(-2 + lmMst1)*pow2(lmMst2) + 16*(lmMst1 + pow3(lmMst2))
        )*pow3(Mst2) + OepS2*(672*(2*Dmglst2 - 3*Mst2)*Mst2*pow2(Mst1) - 8720*
        pow4(Mst1)) + ((383185 - 2592*B4 + 2592*D3 - 1296*DN - 187704*lmMst1 +
        216*lmMst2*(1733 + 859*lmMst2 - 6*lmMst1*(105 + 41*lmMst2) - 26*pow2(
        lmMst1)) - 7992*pow2(lmMst1) + 3456*lmMt*(3 + 6*lmMst2 - 2*lmMst1*(3 +
        lmMst2) + pow2(lmMst1) + pow2(lmMst2)) + 720*pow3(lmMst1) + 58032*pow3(
        lmMst2))*pow4(Mst1))/2. - 243*S2*(4*Mst2*(2*Dmglst2*(65 + 14*lmMst1 -
        14*lmMst2) + 3*(43 - 14*lmMst1 + 14*lmMst2)*Mst2)*pow2(Mst1) + 96*(4*
        Dmglst2 + 3*Mst2)*pow3(Mst2) + (4*(57 - 545*lmMst1 + 545*lmMst2)*pow4(
        Mst1))/3.) + 72*(436 - 6*B4 + 6*D3 - 3*DN - 48*lmMst1 + (408 - 96*
        lmMst1)*lmMst2 + 24*lmMt - 48*(-4 + lmMst1)*pow2(lmMst2) + 48*pow3(
        lmMst2))*pow4(Mst2)) - 3*pow2(z2)*(-(pow4(Msq)*(-6*pow2(Mst1)*(12*
        Dmglst2*Mst2*(2845*Mst2*Mt*s2t + 42*pow2(Mt) - 237*pow2(Mst2)*pow2(s2t)
        ) + 3*pow2(Mst2)*(3184*Mst2*Mt*s2t - 252*pow2(Mt) + 841*pow2(Mst2)*
        pow2(s2t)) + pow2(Dmglst2)*(24958*Mst2*Mt*s2t*xDmglst2 - 84*xDmglst2*
        pow2(Mt) + 530*xDmglst2*pow2(Mst2)*pow2(s2t))) + 18*s2t*(3*Dmglst2*
        Mst2*(-1763*Mt + 158*Mst2*s2t) + (-5072*Mt*xDmglst2 + 209*Mst2*s2t*
        xDmglst2)*pow2(Dmglst2) + 9*(-185*Mt + Mst2*s2t)*pow2(Mst2))*pow3(Mst2)
        + (-20*(17281*Dmglst2 + 4101*Mst2)*Mt*s2t + 19620*pow2(Mt) + (42392*
        Dmglst2 - 35823*Mst2)*Mst2*pow2(s2t))*pow4(Mst1))) + 60*Dmsqst2*pow2(
        Mst2)*pow2(s2t)*(Dmsqst2*xDmsqst2*(75*pow2(Mst1)*pow2(Mst2) + 221*pow4(
        Mst1) + 18*pow4(Mst2)) + pow2(Msq)*(117*pow2(Mst1)*pow2(Mst2) + 221*
        pow4(Mst1) + 27*pow4(Mst2)))) + 2*z4*(Mst2*pow2(s2t)*(4*Dmglst2*pow4(
        Msq)*(4995*pow2(Mst1)*pow2(Mst2) + 11327*pow4(Mst1) + 2862*pow4(Mst2))
        - 3*Mst2*(pow4(Msq)*(6828*pow2(Mst1)*pow2(Mst2) + 13723*pow4(Mst1) +
        1080*pow4(Mst2)) + 20*Dmsqst2*(Dmsqst2*xDmsqst2*(75*pow2(Mst1)*pow2(
        Mst2) + 221*pow4(Mst1) + 18*pow4(Mst2)) + pow2(Msq)*(117*pow2(Mst1)*
        pow2(Mst2) + 221*pow4(Mst1) + 27*pow4(Mst2))))) + pow4(Msq)*(6*
        xDmglst2*pow2(Dmglst2)*(6*pow2(Mst2)*(218*Mst2*Mt*s2t + 162*pow2(Mt) +
        793*pow2(Mst2)*pow2(s2t)) + pow2(Mst1)*(8090*Mst2*Mt*s2t + 2028*pow2(
        Mt) + 7489*pow2(Mst2)*pow2(s2t))) + 36*pow2(Mt)*(6*Mst2*(94*Dmglst2 +
        75*Mst2)*pow2(Mst1) + 162*(2*Dmglst2 + Mst2)*pow3(Mst2) + 1031*pow4(
        Mst1)) + 2*Mt*s2t*(18*(709*Dmglst2 + 1135*Mst2)*pow2(Mst1)*pow2(Mst2) +
        4*(659*Dmglst2 + 8823*Mst2)*pow4(Mst1) + 7317*Dmglst2*pow4(Mst2) +
        6885*pow5(Mst2))))) + pow2(Msq)*(-120*T1ep*pow4(Mst1)*(60*Dmsqst2*pow2(
        Mst2)*pow2(s2t)*(117*pow2(Mst1)*pow2(Mst2) + 221*pow4(Mst1) + 27*pow4(
        Mst2)) + pow2(Msq)*(18*Mst2*pow2(Mst1)*(4*Dmglst2*(253*Mst2*Mt*s2t +
        42*pow2(Mt) - 129*pow2(Mst2)*pow2(s2t)) + Mst2*(-272*Mst2*Mt*s2t - 252*
        pow2(Mt) + 1057*pow2(Mst2)*pow2(s2t))) + (4*(16421*Dmglst2 - 2823*Mst2)
        *Mt*s2t - 19620*pow2(Mt) + Mst2*(-34616*Dmglst2 + 39711*Mst2)*pow2(s2t)
        )*pow4(Mst1) + 54*s2t*(7*Dmglst2*(5*Mt - 2*Mst2*s2t) + 3*Mst2*(-7*Mt +
        23*Mst2*s2t))*pow4(Mst2))) + s2t*(Mt*pow2(Mst1)*(36*pow2(Mst2)*(64800*
        Dmglst2*Dmsqst2*(8 - 15*lmMst1 + 15*lmMst2) + 64800*Dmsqst2*(1 + 3*
        lmMst1 - 3*lmMst2)*Mst2 + 10*Mst2*pow2(Msq)*(75569 + 13716*B4 - 54*DN -
        33426*lmMst1 - 2376*pow2(lmMst1) + 54*lmMst2*(1427 - 1012*lmMst1 + 16*
        pow2(lmMst1)) - 108*(-642 + 203*lmMst1)*pow2(lmMst2) + 21060*pow3(
        lmMst2)) + Dmglst2*pow2(Msq)*(66761 + 301320*B4 - 4860*DN - 205380*
        lmMst1 + 23760*pow2(lmMst1) + 180*lmMst2*(4993 - 1956*lmMst1 + 48*pow2(
        lmMst1)) - 1080*(-482 + 331*lmMst1)*pow2(lmMst2) + 348840*pow3(lmMst2))
        )*pow4(Mst1) + pow2(Mst1)*(1215*(1920*Dmsqst2*(1 + lmMst1 - lmMst2) +
        pow2(Msq)*(9561 + 1760*B4 - 16*DN - 1984*lmMst1 - 256*pow2(lmMst1) -
        64*lmMst2*(-214 + 73*lmMst1 + 4*pow2(lmMst1)) - 32*(-268 + 57*lmMst1)*
        pow2(lmMst2) + 2080*pow3(lmMst2)))*pow5(Mst2) + pow2(Msq)*(160*OepS2*(
        36*(253*Dmglst2 - 68*Mst2)*pow2(Mst1)*pow2(Mst2) + 2*(16421*Dmglst2 -
        2823*Mst2)*pow4(Mst1) + 945*Dmglst2*pow4(Mst2) - 567*pow5(Mst2)) +
        43740*S2*((4*(Dmglst2*(-4978 - 7590*lmMst1 + 7590*lmMst2) + 15*(169 +
        136*lmMst1 - 136*lmMst2)*Mst2)*pow2(Mst1)*pow2(Mst2))/45. - (2*(
        Dmglst2*(123113 + 98526*lmMst1 - 98526*lmMst2) - 3*(9185 + 5646*lmMst1
        - 5646*lmMst2)*Mst2)*pow4(Mst1))/81. + Dmglst2*(90.6 - 70*lmMst1 + 70*
        lmMst2)*pow4(Mst2) + 3*(-1 + 14*lmMst1 - 14*lmMst2)*pow5(Mst2)))) - 27*
        Dmglst2*(pow2(Mst1)*(86400*Dmsqst2*(2 + 5*lmMst1 - 5*lmMst2) - pow2(
        Msq)*(23917 + 188640*B4 - 3600*DN - 37440*lmMst1 + 11520*pow2(lmMst1) -
        2880*lmMst2*(-237 + 55*lmMst1 + 4*pow2(lmMst1)) - 1440*(-280 + 121*
        lmMst1)*pow2(lmMst2) + 185760*pow3(lmMst2)))*pow4(Mst2) + 180*pow2(Msq)
        *((5707.039094650206 - 3416*B4 + 52*DN + (111356*lmMst1)/27. - (294416*
        lmMst2)/27. + 3384*lmMst1*lmMst2 + 8*pow2(lmMst1) - 512*lmMst2*pow2(
        lmMst1) - 4816*pow2(lmMst2) + 5160*lmMst1*pow2(lmMst2) - 64*pow3(
        lmMst1) - 4584*pow3(lmMst2))*pow6(Mst1) - 128*(-1 + 2*lmMst2 + 3*pow2(
        lmMst2))*pow6(Mst2))) + pow2(Msq)*(30*Mst2*(1702429 + 257904*B4 - 648*
        DN - 748656*lmMst1 + 41904*pow2(lmMst1) + 216*lmMst2*(5971 - 6106*
        lmMst1 + 576*pow2(lmMst1)) - 41904*(-34 + 15*lmMst1)*pow2(lmMst2) -
        3456*pow3(lmMst1) + 507600*pow3(lmMst2))*pow6(Mst1) + 622080*pow2(1 +
        lmMst2)*pow7(Mst2))) + 4860*xDR2DRMOD*(64*Mt*pow2(Msq)*pow2(Mst1)*(2*
        pow2(Mst2)*((1 + lmMst2)*Mst2*(1 - 10*lmMst2 + 4*lmMst1*(2 + lmMst2) -
        4*pow2(lmMst2)) + 2*Dmglst2*(8 + 13*lmMst2 - 8*pow2(lmMst2) + lmMst1*(-
        5 + 5*lmMst2 + 6*pow2(lmMst2)) - 6*pow3(lmMst2)))*pow4(Mst1) + (1 +
        lmMst2)*(-2*(1 + 5*lmMst2 - lmMst1*(3 + 2*lmMst2) + 2*pow2(lmMst2))*
        pow2(Mst1)*pow5(Mst2) + Mst2*(7 - 32*lmMst2 + 4*lmMst1*(7 + 3*lmMst2) -
        12*pow2(lmMst2))*pow6(Mst1)) + Dmglst2*((49 + 103*lmMst2 - 36*(1 +
        lmMst2)*pow2(lmMst2) + 4*lmMst1*(-11 + 6*lmMst2 + 9*pow2(lmMst2)))*
        pow6(Mst1) + 2*(pow2(Mst1)*(8 + 7*lmMst2 - 11*pow2(lmMst2) + lmMst1*(-3
        + 5*lmMst2 + 6*pow2(lmMst2)) - 6*pow3(lmMst2))*pow4(Mst2) + (1 - 2*
        lmMst2 - 3*pow2(lmMst2))*pow6(Mst2))) - 2*pow2(1 + lmMst2)*pow7(Mst2))
        - Mst2*s2t*((75*Dmsqst2*(1 - 2*lmMst1 + 2*lmMst2) + pow2(Msq)*(189 +
        726*lmMst2 + 32*(1 + lmMst2)*pow2(lmMst1) + 707*pow2(lmMst2) - 2*
        lmMst1*(253 + 332*lmMst2 + 123*pow2(lmMst2)) + 214*pow3(lmMst2)))*pow4(
        Mst1)*pow5(Mst2) + 2*pow2(Mst2)*(30*Dmglst2*Dmsqst2*(lmMst1 - lmMst2) -
        75*Dmsqst2*(lmMst1 - lmMst2)*Mst2 + 2*Dmglst2*pow2(Msq)*(48 + 4*lmMst2*
        (31 + 36*pow2(lmMst1)) + 278*pow2(lmMst2) - lmMst1*(44 + 278*lmMst2 +
        379*pow2(lmMst2)) + 235*pow3(lmMst2)) + Mst2*pow2(Msq)*(32 + 285*lmMst2
        + 144*(1 + lmMst2)*pow2(lmMst1) + 444*pow2(lmMst2) - lmMst1*(253 + 588*
        lmMst2 + 379*pow2(lmMst2)) + 235*pow3(lmMst2)))*pow6(Mst1) + pow2(Mst1)
        *(-(Dmglst2*(30*Dmsqst2*((1 - 2*lmMst1 + 2*lmMst2)*pow2(Mst1) + pow2(
        Mst2)) - 2*pow2(Msq)*((-20 + 2*(99 + 16*lmMst1)*lmMst2 + 123*pow2(
        lmMst2))*pow2(Mst2) + pow2(Mst1)*(60 + 206*lmMst2 + 32*lmMst2*pow2(
        lmMst1) + lmMst1*(8 - 460*lmMst2 - 246*pow2(lmMst2)) + 519*pow2(lmMst2)
        + 214*pow3(lmMst2))))*pow4(Mst2)) - 9*Mst2*shiftst1*(Dmsqst2 + pow2(
        Msq))*(15*(pow2(Mst1) + pow2(Mst2))*pow4(Mst2) - 30*(lmMst1 - lmMst2)*
        pow2(Mst1)*(pow2(Mst1)*pow2(Mst2) + pow4(Mst1) + pow4(Mst2))) + (75*
        Dmsqst2 + (205 + 252*lmMst2 + 32*lmMst1*(1 + lmMst2) + 91*pow2(lmMst2))
        *pow2(Msq))*pow7(Mst2)) + (4*Dmglst2*pow2(Msq)*(48 + 4*lmMst2*(45 + 68*
        pow2(lmMst1)) + 310*pow2(lmMst2) - lmMst1*(92 + 310*lmMst2 + 635*pow2(
        lmMst2)) + 363*pow3(lmMst2)) + 2*Mst2*(-75*Dmsqst2*(lmMst1 - lmMst2) +
        pow2(Msq)*(40 + 277*lmMst2 + 272*(1 + lmMst2)*pow2(lmMst1) + 556*pow2(
        lmMst2) - lmMst1*(237 + 828*lmMst2 + 635*pow2(lmMst2)) + 363*pow3(
        lmMst2))))*pow8(Mst1) - 16*(1 + lmMst2)*(4*Dmglst2*lmMst2 + Mst2 +
        lmMst2*Mst2)*pow2(Msq)*pow8(Mst2))))))))/(3.000564e7*pow4(Msq)*pow4(
        Mst1)*pow6(Mst2)))) + pow2(Mt)*((pow2(MuSUSY)*(91125*oneLoopFlag*pow2(
        s2t)*(2 - lmMst1 + lmMst2 - (2*(lmMst1 - lmMst2)*pow2(Mst1)*(pow2(Mst1)
        + pow2(Mst2)))/pow4(Mst2)) + (Al4p*(13500*Mst2*s2t*(-(xDR2DRMOD*((9*
        Mst2*s2t*(4*(2*Dmglst2*lmMst2 + Mst2 + lmMst2*Mst2)*twoLoopFlag*((pow2(
        Mst1) + pow2(Mst2))*pow4(Mst2) - 2*(lmMst1 - lmMst2)*pow2(Mst1)*(pow2(
        Mst1)*pow2(Mst2) + pow4(Mst1) + pow4(Mst2))) - (5*Al4p*Mst2*
        threeLoopFlag*xDmsqst2*pow2(Dmsqst2)*(pow2(Mst2)*(pow2(Mst1) + pow2(
        Mst2))*(-10*Dmglst2*Mst2 + 7*xDmglst2*pow2(Dmglst2) + (-11 + 18*
        shiftst1)*pow2(Mst2)) + 2*lmMst2*pow2(Mst1)*(pow2(Mst2)*(-10*Dmglst2*
        Mst2 + 7*xDmglst2*pow2(Dmglst2) + (-11 + 18*shiftst1)*pow2(Mst2)) +
        pow2(Mst1)*(4*Dmglst2*Mst2 + 26*xDmglst2*pow2(Dmglst2) + (-11 + 18*
        shiftst1)*pow2(Mst2)) + 2*(-5 + 9*shiftst1)*pow4(Mst1)) - 2*lmMst1*(
        pow2(Mst1)*pow2(Mst2)*(-10*Dmglst2*Mst2 + 7*xDmglst2*pow2(Dmglst2) + (-
        11 + 18*shiftst1)*pow2(Mst2)) + (4*Dmglst2*Mst2 + 26*xDmglst2*pow2(
        Dmglst2) + (-11 + 18*shiftst1)*pow2(Mst2))*pow4(Mst1) + 2*(-5 + 9*
        shiftst1)*pow6(Mst1))))/pow4(Msq)))/pow2(Mst1) + 2*xDmglst2*pow2(
        Dmglst2)*((-18*(2 - lmMst2)*s2t*twoLoopFlag*((pow2(Mst1) + pow2(Mst2))*
        pow4(Mst2) - 2*(lmMst1 - lmMst2)*pow2(Mst1)*(pow2(Mst1)*pow2(Mst2) +
        pow4(Mst1) + pow4(Mst2))))/pow2(Mst1) + Al4p*Mst2*threeLoopFlag*(Mt*(-
        64*pow2(Mst2)*(5 + 149*lmMst2 + 12*pow2(lmMst2) + 6*lmMst1*(-7 - 8*
        lmMst2 + 6*pow2(lmMst2)) - 36*pow3(lmMst2)) + 64*pow2(Mst1)*(85 - 215*
        lmMst2 + lmMst1*(66 + 90*lmMst2 - 72*pow2(lmMst2)) - 54*pow2(lmMst2) +
        72*pow3(lmMst2))) + (3*(128*Mt*(-3 - 4*lmMst2 + 3*pow2(lmMst2)) + Mst2*
        s2t*(-380 + 32*lmMst1*(-2 + lmMst2) - 38*lmMst2 + 251*pow2(lmMst2)))*
        pow4(Mst2))/pow2(Mst1) + (Mst2*s2t*((585*Dmsqst2*pow2(Mst1)*(2*lmMst1*
        pow2(Mst1) - 2*lmMst2*pow2(Mst1) - pow2(Mst2))*(pow2(Mst1) + pow2(Mst2)
        ))/pow2(Msq) + pow2(Mst2)*(-1540 - 2506*lmMst2 + 96*(-2 + lmMst2)*pow2(
        lmMst1) + lmMst1*(2760 + 36*lmMst2 - 738*pow2(lmMst2)) + 141*pow2(
        lmMst2) + 642*pow3(lmMst2))*pow4(Mst1) + 6*(-336 - 716*lmMst2 + 144*(-2
        + lmMst2)*pow2(lmMst1) + lmMst1*(668 + 534*lmMst2 - 379*pow2(lmMst2)) -
        246*pow2(lmMst2) + 235*pow3(lmMst2))*pow6(Mst1) + 96*(2 + lmMst2 - 3*
        pow2(lmMst2))*pow6(Mst2)))/pow4(Mst1))))) + (9*twoLoopFlag*(4*(4*Mt*(5
        + 6*lmMst2 - lmMst1*(4 + 3*lmMst2) + 3*pow2(lmMst2)) + Mst2*s2t*(-1 +
        13*lmMst2 - lmMst1*(13 + 8*lmMst2) + pow2(lmMst1) + 7*pow2(lmMst2)))*
        pow3(Mst2)*pow4(Mst1) + 2*(-(Mst2*s2t*(-14 - 20*lmMst2 + 2*lmMst1*(5 +
        3*lmMst2) + pow2(lmMst1) - 7*pow2(lmMst2))) + 8*Mt*(4 + 3*lmMst2 -
        lmMst1*(1 + lmMst2) + pow2(lmMst2)))*pow2(Mst1)*pow5(Mst2) + 4*Dmglst2*
        (Mst2*s2t*(-5 + 8*lmMst2 - 4*lmMst1*(2 + lmMst2) + 4*pow2(lmMst2)) +
        Mt*(65 + lmMst1*(34 - 20*lmMst2) - 26*lmMst2 + 20*pow2(lmMst2)))*pow6(
        Mst1) + Mst2*(Mst2*s2t*(-1 + 50*lmMst2 - 2*lmMst1*(25 + 32*lmMst2) +
        20*pow2(lmMst1) + 44*pow2(lmMst2)) + Mt*(84 + 152*lmMst2 - 40*lmMst1*(3
        + 2*lmMst2) + 80*pow2(lmMst2)))*pow6(Mst1) + 8*Dmglst2*((Mst2*s2t*(-2 +
        3*lmMst2 - lmMst1*(3 + 2*lmMst2) + 2*pow2(lmMst2)) + Mt*(22 + lmMst1*(8
        - 6*lmMst2) - 4*lmMst2 + 6*pow2(lmMst2)))*pow2(Mst2)*pow4(Mst1) + (-(
        Mst2*s2t*(1 + lmMst1 - 2*lmMst2 + 2*lmMst1*lmMst2 - 2*pow2(lmMst2))) +
        2*Mt*(6 + lmMst1 + lmMst2 - lmMst1*lmMst2 + pow2(lmMst2)))*pow2(Mst1)*
        pow4(Mst2) + lmMst2*s2t*pow7(Mst2)) + 4*(1 + lmMst2)*s2t*pow8(Mst2)))/
        pow2(Mst1)) + (2*xDmglst2*pow2(Dmglst2)*(-121500*s2t*twoLoopFlag*pow2(
        Mst1)*(2*(-8*(10 + lmMst1 - lmMst2)*Mt + Mst2*s2t*(8 - 9*lmMst1 + 9*
        lmMst2 + 2*lmMst1*lmMst2 - 2*pow2(lmMst2)))*pow2(Mst2)*pow4(Mst1) + 2*(
        -4*(8 + lmMst1 - lmMst2)*Mt + Mst2*s2t*(8 - 5*lmMst1 + 4*lmMst2 + 2*
        lmMst1*lmMst2 - 2*pow2(lmMst2)))*pow2(Mst1)*pow4(Mst2) + (-268*Mt +
        Mst2*s2t*(17 + 4*lmMst1*(-7 + lmMst2) + 28*lmMst2 - 4*pow2(lmMst2)))*
        pow6(Mst1) - 2*(-2 + lmMst2)*s2t*pow7(Mst2)) + Al4p*Mst2*threeLoopFlag*
        (75*pow2(Mst2)*(-8*Mst2*Mt*s2t*(334138 - 61560*B4 + 1620*DN - 236520*
        lmMst1 - 180*(-823 + 438*lmMst1)*lmMst2 + 2240*OepS2 - 81*(-1373 + 560*
        lmMst1 - 560*lmMst2)*S2 - 3360*T1ep - 4320*pow2(lmMst1) + 1080*(29 +
        48*lmMst1)*pow2(lmMst2) - 51840*pow3(lmMst2)) + pow2(Mst2)*pow2(s2t)*(
        13169 - 41040*B4 + 43200*D3 - 22680*DN + 282960*lmMst1 + 1120*OepS2 -
        324*(65819 + 70*lmMst1 - 70*lmMst2)*S2 - 1680*T1ep - 13500*pow2(lmMst1)
        + 720*lmMst2*(-775 + 12*lmMst1 + 24*pow2(lmMst1)) - 1080*(-2 + 123*
        lmMst1)*pow2(lmMst2) + 115560*pow3(lmMst2)) + 96*pow2(Mt)*(33934 - 90*
        B4 + 90*D3 - 45*DN + 120*(163 + 24*lmMst1)*lmMst2 - 120*lmMt - 10206*S2
        - 720*(2 + lmMst1)*pow2(lmMst2) + 720*(lmMst1 + pow3(lmMst2))))*pow4(
        Mst1) - 40500*s2t*(Mst2*s2t*(812 - 32*lmMst1*(-2 + lmMst2) + 38*lmMst2
        - 251*pow2(lmMst2)) + 128*Mt*(3 + 4*lmMst2 - 3*pow2(lmMst2)))*pow2(
        Mst1)*pow5(Mst2) + 2*(-2*pow2(Mst2)*pow2(s2t)*(2011073 + 1417500*B4 -
        1458000*D3 + 749250*DN + 934245*lmMst1 - 1178000*OepS2 + 1350*(620417 +
        17670*lmMst1 - 17670*lmMst2)*S2 + 1767000*T1ep + 3150900*pow2(lmMst1) -
        45*lmMst2*(-124139 + 189090*lmMst1 + 71550*pow2(lmMst1)) + 4050*(1323 +
        1970*lmMst1)*pow2(lmMst2) + 101250*pow3(lmMst1) - 4860000*pow3(lmMst2))
        - 125*Mst2*Mt*s2t*(996211 - 295488*B4 + 7776*DN - 1030176*lmMst1 - 144*
        (-1883 + 3618*lmMst1)*lmMst2 + 98336*OepS2 - 756*(259 + 2634*lmMst1 -
        2634*lmMst2)*S2 - 147504*T1ep + 49248*pow2(lmMst1) + 5184*(67 + 48*
        lmMst1)*pow2(lmMst2) - 248832*pow3(lmMst2)) + 150*pow2(Mt)*(2199511 -
        4320*B4 + 4320*D3 - 2160*DN + 140160*lmMst1 + 960*(1426 + 303*lmMst1)*
        lmMst2 - 2880*(-16 + 5*lmMst1 - 5*lmMst2)*lmMt - 1120*OepS2 + 324*(-
        3411 + 70*lmMst1 - 70*lmMst2)*S2 + 1680*T1ep - 2880*(77 + 24*lmMst1)*
        pow2(lmMst2) + 69120*pow3(lmMst2)))*pow6(Mst1) + 1296000*(2 + lmMst2 -
        3*pow2(lmMst2))*pow2(s2t)*pow8(Mst2))))/pow4(Mst1)))/pow7(Mst2)))/
        729000. - (s2t*xMst*(27*(lmMst1 - lmMst2)*oneLoopFlag*s2t*pow2(MuSUSY)*
        pow3(Mst2) + Al4p*twoLoopFlag*(-8*Mst2*(-18*(-2 + lmMst2)*(-lmMst1 +
        lmMst2)*s2t*xDmglst2*xDR2DRMOD*pow2(Dmglst2) + Dmglst2*Mt*(785 + 6*
        lmMst1*(85 - 42*lmMst2) - 438*lmMst2 + 252*pow2(lmMst2)) + Mst2*Mt*(193
        + 474*lmMst2 - 6*lmMst1*(67 + 42*lmMst2) + 252*pow2(lmMst2)) - Dmglst2*
        Mst2*s2t*(49 - 84*lmMst2 + lmMst1*(84 - 36*lmMst2*(-1 + xDR2DRMOD)) +
        36*(-1 + xDR2DRMOD)*pow2(lmMst2)) - s2t*(1 + 3*lmMst2*(-37 + 6*
        xDR2DRMOD) - 3*lmMst1*(-37 + 6*lmMst2*(-12 + xDR2DRMOD) + 6*xDR2DRMOD)
        - 81*pow2(lmMst1) + 9*(-15 + 2*xDR2DRMOD)*pow2(lmMst2))*pow2(Mst2))*
        pow2(MuSUSY) + (2*xDmglst2*pow2(Dmglst2)*(48*(-143 + 18*lmMst1 - 18*
        lmMst2)*Mt*Tbeta*pow2(MuSUSY) + 2*Mst2*MuSUSY*(MuSUSY*s2t*Tbeta*(157 -
        348*lmMst1 + 348*lmMst2 + 36*lmMst1*lmMst2 - 36*pow2(lmMst2)) + 30*(-43
        + 60*lmMst1 - 60*lmMst2)*Mt*pow2(Sbeta)) + 15*(-43 + 60*lmMst1 - 60*
        lmMst2)*s2t*Tbeta*(-1 + pow2(Sbeta))*pow2(Sbeta)*pow3(Mst2)))/Tbeta))*
        pow6(Mst1))/(108.*pow9(Mst2))) + (Al4p*Mt*z2*(Al4p*threeLoopFlag*((-75*
        s2t*xDmsqst2*pow2(Dmsqst2)*(s2t*(81*Mst2*pow2(Sbeta)*(-2*(11 + 10*
        lmMst1 - 10*lmMst2)*s2t*pow2(Mst2)*(-1 + pow2(Sbeta)) + (Dmglst2*(372 -
        40*lmMst1 + 40*lmMst2) + 3*(-17 + 2*lmMst1 - 2*lmMst2)*Mst2)*Mt*pow2(
        Sbeta))*pow6(Mst1) + Mt*(6*Mst2*pow2(MuSUSY)*(-6*(25 + 9*lmMst1 - 9*
        lmMst2)*pow2(Mst1)*pow3(Mst2) + (1607 - 432*lmMst1 + 432*lmMst2)*Mst2*
        pow4(Mst1) + 216*Dmglst2*(3*(2 + lmMst1 - lmMst2)*pow2(Mst1)*pow2(Mst2)
        + (3 + 8*lmMst1 - 8*lmMst2)*pow4(Mst1) + pow4(Mst2)) - 162*pow5(Mst2))
        + ((41402 - 1296*lmMst1 + 1296*lmMst2)*pow2(MuSUSY) + 81*Mst2*(4*
        Dmglst2*(-93 + 10*lmMst1 - 10*lmMst2) + 3*(17 - 2*lmMst1 + 2*lmMst2)*
        Mst2)*pow2(Sbeta))*pow6(Mst1))) + 648*Mt*pow2(MuSUSY)*(pow2(Mst1)*(8*
        Dmglst2*Mt*((1 - 9*lmMst1 + 9*lmMst2)*pow2(Mst1) + (-2 - 3*lmMst1 + 3*
        lmMst2)*pow2(Mst2)) + (4*(7 + 3*lmMst1 - 3*lmMst2)*Mt + (-1 + 2*lmMst1
        - 2*lmMst2)*Mst2*s2t*shiftst1)*pow3(Mst2)) + 2*Mst2*(4*(5 + 3*lmMst1 -
        3*lmMst2)*Mt + (lmMst1 - lmMst2)*Mst2*s2t*shiftst1)*pow4(Mst1) + s2t*
        shiftst1*(2*(lmMst1 - lmMst2)*pow6(Mst1) - pow6(Mst2)))))/(pow2(Mst1)*
        pow4(Msq)*pow4(Mst2)) + Mt*(7290*pow2(s2t)*((10*Dmglst2*Dmsqst2*(-93 +
        10*lmMst1 - 10*lmMst2)*(-1 + pow2(Sbeta))*pow2(Sbeta)*pow4(Mst1))/(3.*
        pow2(Msq)*pow3(Mst2)) - pow2(MuSUSY)*(103.64814814814815 + (94*lmMst1)/
        9. - (184*lmMst2)/9. - (Dmglst2*(60.407407407407405 - 24*lmMst1 + (560*
        lmMst2)/9. - (40*Dmsqst2*(2 + lmMst1 - lmMst2))/pow2(Msq)))/Mst2 + (8*
        Dmglst2*Mst2*(7 + (15*Dmsqst2)/pow2(Msq)))/(9.*pow2(Mst1)) - ((
        4.222222222222222 + (20*Dmsqst2)/(3.*pow2(Msq)))*pow2(Mst2))/pow2(Mst1)
        + ((265*Dmsqst2)/9. + ((10800*Dmglst2*Dmsqst2*(3 + 8*lmMst1 - 8*lmMst2)
        + 150*Dmsqst2*(1097 - 72*lmMst1 + 72*lmMst2)*Mst2 + 36*Dmglst2*(-6803 +
        1810*lmMst1 - 2670*lmMst2)*pow2(Msq) + (533629 - 900*lmMst1 + 180*
        lmMst2)*Mst2*pow2(Msq))*pow2(Mst1))/(810.*pow3(Mst2)))/pow2(Msq) - (20*
        shiftst1*(1 + Dmsqst2/pow2(Msq))*(1 + pow2(Mst2)/pow2(Mst1) - (2*(
        lmMst1 - lmMst2)*(pow2(Mst1)*pow2(Mst2) + pow4(Mst1) + pow4(Mst2)))/
        pow4(Mst2)))/3. - ((Dmglst2*(1167.8951989026064 - (1520*lmMst1)/9. + (
        1864*lmMst2)/9.) - (Mst2*(6649153 - 81540*lmMst1 + 77220*lmMst2 + (100*
        Dmsqst2*(20701 - 648*lmMst1 + 648*lmMst2))/pow2(Msq)))/4860.)*pow4(
        Mst1))/pow5(Mst2) - (2*shiftst3*(2*(1 - 2*lmMst1 + 2*lmMst2)*pow2(Mst2)
        *pow4(Mst1) + 2*(1 - lmMst1 + lmMst2)*pow2(Mst1)*pow4(Mst2) + (2 - 6*
        lmMst1 + 6*lmMst2)*pow6(Mst1) + pow6(Mst2)))/(3.*pow2(Mst1)*pow4(Mst2))
        )) + (pow2(MuSUSY)*(90*pow2(Mt)*(-6*Mst2*(Dmglst2*(7286 - 240*lmMst1 +
        6384*lmMst2) + (-353 + 72*lmMst1 + 696*lmMst2)*Mst2)*pow2(Mst1) - 96*(
        24*Dmglst2*(5 + 6*lmMst2) + (53 + 24*lmMst2)*Mst2)*pow3(Mst2) + (38401
        + 1080*lmMst1 - 7992*lmMst2)*pow4(Mst1)) + s2t*(Mt*((-388800*Dmsqst2*
        pow2(Mst2)*(Dmglst2*(1 - 9*lmMst1 + 9*lmMst2)*pow2(Mst1) + (5 + 3*
        lmMst1 - 3*lmMst2)*Mst2*pow2(Mst1) + Dmglst2*(-2 - 3*lmMst1 + 3*lmMst2)
        *pow2(Mst2) + (2 + lmMst1 - lmMst2)*pow3(Mst2)))/pow2(Msq) - 3*Mst2*(
        120*(3937 + 2988*lmMst1 - 2232*lmMst2)*pow2(Mst1)*pow2(Mst2) + (533726
        + 792720*lmMst1 - 693360*lmMst2)*pow4(Mst1) + 135*(2269 + 664*lmMst1 -
        56*lmMst2)*pow4(Mst2)) - Dmglst2*(36*(364291 - 88560*lmMst1 + 147960*
        lmMst2)*pow2(Mst1)*pow2(Mst2) + 2*(15057833 - 4014360*lmMst1 + 5563080*
        lmMst2)*pow4(Mst1) + 27*(173947 - 25080*lmMst1 + 68760*lmMst2)*pow4(
        Mst2))) - (6480*xDR2DRMOD*(2*(32*(1 + lmMst2)*Mt + (-4 - 13*lmMst1 + 9*
        lmMst2)*Mst2*s2t)*pow3(Mst2)*pow4(Mst1) + (32*(1 + lmMst2)*Mt + (5 -
        26*lmMst1 + 18*lmMst2)*Mst2*s2t)*pow2(Mst1)*pow5(Mst2) - 4*Mst2*
        xDmglst2*pow2(Dmglst2)*((-8*(-2 + 3*lmMst2)*Mt + (1 - 10*lmMst1 + 12*
        lmMst2)*Mst2*s2t)*pow2(Mst1)*pow2(Mst2) - 2*(8*(-2 + 3*lmMst2)*Mt + (2
        + 5*lmMst1 - 6*lmMst2)*Mst2*s2t)*pow4(Mst1) + 5*s2t*pow5(Mst2)) + 2*
        Mst2*(48*(1 + lmMst2)*Mt + (-4 - 13*lmMst1 + 9*lmMst2)*Mst2*s2t)*pow6(
        Mst1) + Dmglst2*(2*(32*(1 + 3*lmMst2)*Mt + (7*lmMst1 - 15*lmMst2)*Mst2*
        s2t)*pow2(Mst2)*pow4(Mst1) + (32*(1 + 3*lmMst2)*Mt + (-7 + 14*lmMst1 -
        30*lmMst2)*Mst2*s2t)*pow2(Mst1)*pow4(Mst2) + 2*(48*(1 + 3*lmMst2)*Mt +
        (7*lmMst1 - 15*lmMst2)*Mst2*s2t)*pow6(Mst1) - 7*s2t*pow7(Mst2)) + 13*
        s2t*pow8(Mst2)))/pow2(Mst1))))/pow6(Mst2))) - (3*Mt*pow2(MuSUSY)*(1620*
        s2t*twoLoopFlag*pow2(Mst1)*pow2(Mst2)*(4*Mt*pow2(Mst2)*pow4(Mst1) + 4*
        Mt*pow2(Mst1)*pow4(Mst2) + 4*Mt*xMst*pow6(Mst1) + (4*Mt - Mst2*s2t)*
        pow6(Mst2)) - 6480*Dmglst2*Mst2*s2t*twoLoopFlag*pow2(Mst1)*((-13*Mt +
        Mst2*s2t)*pow2(Mst2)*pow4(Mst1) + (-9*Mt + Mst2*s2t)*pow2(Mst1)*pow4(
        Mst2) + (-17*Mt*xMst + Mst2*s2t*xMst)*pow6(Mst1) + (-5*Mt + Mst2*s2t)*
        pow6(Mst2)) + xDmglst2*pow2(Dmglst2)*(-3240*s2t*twoLoopFlag*pow2(Mst1)*
        ((-36*Mt + Mst2*s2t)*pow2(Mst2)*pow4(Mst1) + (-24*Mt + Mst2*s2t)*pow2(
        Mst1)*pow4(Mst2) + (-48*Mt*xMst + Mst2*s2t*xMst)*pow6(Mst1) + (-12*Mt +
        Mst2*s2t)*pow6(Mst2)) + Al4p*threeLoopFlag*pow3(Mst2)*(3*pow2(Mst1)*
        pow2(Mst2)*(8*(46987 + 15390*lmMst1 + 4050*lmMst2)*Mst2*Mt*s2t + 192*(
        17 + 1320*lmMst2)*pow2(Mt) + (-40001 - 52560*lmMst1 + 37080*lmMst2)*
        pow2(Mst2)*pow2(s2t)) + 2*(5*(-184517 + 74952*lmMst1 + 18360*lmMst2)*
        Mst2*Mt*s2t + 6*(12761 + 5400*lmMst1 + 178920*lmMst2)*pow2(Mt) + (
        261247 - 163890*lmMst1 + 140670*lmMst2)*pow2(Mst2)*pow2(s2t))*pow4(
        Mst1) - 28080*pow2(s2t)*pow6(Mst2)))))/(pow2(Mst1)*pow9(Mst2))))/7290.;
}

/**
 * @return The diagonal (2, 2) matrix element of the Higgs mass matrix as a double for the hierarchy 'H6bq22g'
 */
double H6bq22g::getS2() const {
   return -(oneLoopFlag*((4*Mt*MuSUSY*s2t*(4*(lmMst1 - lmMst2)*pow2(Mt) + ((-2 -
        lmMst1 + lmMst2)*pow2(Mst1) + (2 - lmMst1 + lmMst2)*pow2(Mst2))*pow2(
        s2t)))/Tbeta + 4*pow2(Mt)*pow2(s2t)*(2*(lmMst1 - lmMst2)*(pow2(Mst1) -
        pow2(Mst2)) + pow2(MuSUSY)*(2 - lmMst1 + lmMst2 - (2*(lmMst1 - lmMst2)*
        pow2(Mst1)*(pow2(Mst1) + pow2(Mst2)))/pow4(Mst2))) - (4*pow2(Mt)*pow2(
        MuSUSY)*pow2(s2t)*(2 - lmMst1 + lmMst2 - (2*(lmMst1 - lmMst2)*pow2(
        Mst1)*(pow2(Mst1) + pow2(Mst2)))/pow4(Mst2)))/pow2(Sbeta) + 16*(lmMst1
        + lmMst2 - 2*lmMt)*pow4(Mt) - (-4*pow2(Mst1)*pow2(Mst2) + (2 + lmMst1 -
        lmMst2)*pow4(Mst1) + (2 - lmMst1 + lmMst2)*pow4(Mst2))*pow4(s2t)))/32.
         - Al4p*((xDR2DRMOD*(2700*Mst2*(2*Dmglst2*lmMst2 + Mst2 + lmMst2*Mst2)*
        twoLoopFlag*pow2(Mst1)*pow4(Msq)*(-(pow2(Mst2)*pow2(s2t)*pow4(Mst1)*(-
        8*(lmMst1 - lmMst2)*Tbeta*pow2(Mt)*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 4*
        Mt*(-(MuSUSY*s2t) + 2*Mt*Tbeta)*pow2(Mst2)*pow2(Sbeta) + (1 - 2*lmMst1
        + 2*lmMst2)*Tbeta*pow2(s2t)*pow2(Sbeta)*pow4(Mst2))) - pow2(Mst1)*pow4(
        Mst2)*(-4*Tbeta*pow2(Mt)*pow2(s2t)*((-1 + 2*lmMst1 - 2*lmMst2)*pow2(
        MuSUSY)*(-1 + pow2(Sbeta)) + 4*pow2(Mst2)*pow2(Sbeta)) + pow2(Sbeta)*(
        16*MuSUSY*s2t*pow3(Mt) + 8*(-lmMst1 + lmMst2)*Mt*MuSUSY*pow2(Mst2)*
        pow3(s2t) - 16*Tbeta*pow4(Mt) + (1 + 2*lmMst1 - 2*lmMst2)*Tbeta*pow4(
        Mst2)*pow4(s2t))) + Tbeta*pow2(s2t)*(8*(lmMst1 - lmMst2)*pow2(Mt)*pow2(
        MuSUSY)*(-1 + pow2(Sbeta)) + pow2(s2t)*pow2(Sbeta)*pow4(Mst2))*pow6(
        Mst1) + (-4*Tbeta*pow2(Mt)*pow2(s2t)*(pow2(MuSUSY)*(-1 + pow2(Sbeta)) +
        2*pow2(Mst2)*pow2(Sbeta)) + pow2(Sbeta)*(16*MuSUSY*s2t*pow3(Mt) - 4*Mt*
        MuSUSY*pow2(Mst2)*pow3(s2t) + 16*Tbeta*pow4(Mt) + Tbeta*pow4(Mst2)*
        pow4(s2t)))*pow6(Mst2)) - 3375*Al4p*threeLoopFlag*xDmsqst2*pow2(
        Dmsqst2)*pow2(Mst1)*pow2(Mst2)*(xDmglst2*pow2(Dmglst2)*(-(pow2(s2t)*
        pow4(Mst1)*(-208*(lmMst1 - lmMst2)*Tbeta*pow2(Mt)*pow2(MuSUSY)*(-1 +
        pow2(Sbeta)) + 28*Mt*(-(MuSUSY*s2t) + 2*Mt*Tbeta)*pow2(Mst2)*pow2(
        Sbeta) + 7*(1 - 2*lmMst1 + 2*lmMst2)*Tbeta*pow2(s2t)*pow2(Sbeta)*pow4(
        Mst2))) + 7*pow4(Mst2)*(-4*Tbeta*pow2(Mt)*pow2(s2t)*(pow2(MuSUSY)*(-1 +
        pow2(Sbeta)) + 2*pow2(Mst2)*pow2(Sbeta)) + pow2(Sbeta)*(16*MuSUSY*s2t*
        pow3(Mt) - 4*Mt*MuSUSY*pow2(Mst2)*pow3(s2t) + 16*Tbeta*pow4(Mt) +
        Tbeta*pow4(Mst2)*pow4(s2t))) - 7*pow2(Mst1)*pow2(Mst2)*(-4*Tbeta*pow2(
        Mt)*pow2(s2t)*((-1 + 2*lmMst1 - 2*lmMst2)*pow2(MuSUSY)*(-1 + pow2(
        Sbeta)) + 4*pow2(Mst2)*pow2(Sbeta)) + pow2(Sbeta)*(16*MuSUSY*s2t*pow3(
        Mt) + 8*(-lmMst1 + lmMst2)*Mt*MuSUSY*pow2(Mst2)*pow3(s2t) - 16*Tbeta*
        pow4(Mt) + (1 + 2*lmMst1 - 2*lmMst2)*Tbeta*pow4(Mst2)*pow4(s2t)))) +
        Tbeta*pow2(s2t)*(16*(lmMst1 - lmMst2)*(-5 + 9*shiftst1)*pow2(Mt)*pow2(
        MuSUSY)*(-1 + pow2(Sbeta)) + (-11 + 18*shiftst1)*pow2(s2t)*pow2(Sbeta)*
        pow4(Mst2))*pow6(Mst1) - 2*Dmglst2*(-(Mst2*pow2(s2t)*pow4(Mst1)*(16*(
        lmMst1 - lmMst2)*Tbeta*pow2(Mt)*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 20*
        Mt*(-(MuSUSY*s2t) + 2*Mt*Tbeta)*pow2(Mst2)*pow2(Sbeta) + 5*(1 - 2*
        lmMst1 + 2*lmMst2)*Tbeta*pow2(s2t)*pow2(Sbeta)*pow4(Mst2))) - 5*pow2(
        Mst1)*pow3(Mst2)*(-4*Tbeta*pow2(Mt)*pow2(s2t)*((-1 + 2*lmMst1 - 2*
        lmMst2)*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 4*pow2(Mst2)*pow2(Sbeta)) +
        pow2(Sbeta)*(16*MuSUSY*s2t*pow3(Mt) + 8*(-lmMst1 + lmMst2)*Mt*MuSUSY*
        pow2(Mst2)*pow3(s2t) - 16*Tbeta*pow4(Mt) + (1 + 2*lmMst1 - 2*lmMst2)*
        Tbeta*pow4(Mst2)*pow4(s2t))) + 5*((-4*Tbeta*pow2(Mt)*pow2(s2t)*(pow2(
        MuSUSY)*(-1 + pow2(Sbeta)) + 2*pow2(Mst2)*pow2(Sbeta)) + pow2(Sbeta)*(
        16*MuSUSY*s2t*pow3(Mt) - 4*Mt*MuSUSY*pow2(Mst2)*pow3(s2t) + 16*Tbeta*
        pow4(Mt) + Tbeta*pow4(Mst2)*pow4(s2t)))*pow5(Mst2) + Tbeta*pow2(Sbeta)*
        pow3(Mst2)*pow4(s2t)*pow6(Mst1))) + (11 - 18*shiftst1)*(pow2(Mst2)*
        pow2(s2t)*pow4(Mst1)*(-8*(lmMst1 - lmMst2)*Tbeta*pow2(Mt)*pow2(MuSUSY)*
        (-1 + pow2(Sbeta)) + 4*Mt*(-(MuSUSY*s2t) + 2*Mt*Tbeta)*pow2(Mst2)*pow2(
        Sbeta) + (1 - 2*lmMst1 + 2*lmMst2)*Tbeta*pow2(s2t)*pow2(Sbeta)*pow4(
        Mst2)) + pow2(Mst1)*pow4(Mst2)*(-4*Tbeta*pow2(Mt)*pow2(s2t)*((-1 + 2*
        lmMst1 - 2*lmMst2)*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 4*pow2(Mst2)*pow2(
        Sbeta)) + pow2(Sbeta)*(16*MuSUSY*s2t*pow3(Mt) + 8*(-lmMst1 + lmMst2)*
        Mt*MuSUSY*pow2(Mst2)*pow3(s2t) - 16*Tbeta*pow4(Mt) + (1 + 2*lmMst1 - 2*
        lmMst2)*Tbeta*pow4(Mst2)*pow4(s2t))) - (-4*Tbeta*pow2(Mt)*pow2(s2t)*(
        pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 2*pow2(Mst2)*pow2(Sbeta)) + pow2(
        Sbeta)*(16*MuSUSY*s2t*pow3(Mt) - 4*Mt*MuSUSY*pow2(Mst2)*pow3(s2t) + 16*
        Tbeta*pow4(Mt) + Tbeta*pow4(Mst2)*pow4(s2t)))*pow6(Mst2))) + 2*
        xDmglst2*pow2(Dmglst2)*pow2(Msq)*(-1350*(2 - lmMst2)*twoLoopFlag*pow2(
        Msq)*pow2(Mst1)*(-(pow2(Mst2)*pow2(s2t)*pow4(Mst1)*(-8*(lmMst1 -
        lmMst2)*Tbeta*pow2(Mt)*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 4*Mt*(-(
        MuSUSY*s2t) + 2*Mt*Tbeta)*pow2(Mst2)*pow2(Sbeta) + (1 - 2*lmMst1 + 2*
        lmMst2)*Tbeta*pow2(s2t)*pow2(Sbeta)*pow4(Mst2))) - pow2(Mst1)*pow4(
        Mst2)*(-4*Tbeta*pow2(Mt)*pow2(s2t)*((-1 + 2*lmMst1 - 2*lmMst2)*pow2(
        MuSUSY)*(-1 + pow2(Sbeta)) + 4*pow2(Mst2)*pow2(Sbeta)) + pow2(Sbeta)*(
        16*MuSUSY*s2t*pow3(Mt) + 8*(-lmMst1 + lmMst2)*Mt*MuSUSY*pow2(Mst2)*
        pow3(s2t) - 16*Tbeta*pow4(Mt) + (1 + 2*lmMst1 - 2*lmMst2)*Tbeta*pow4(
        Mst2)*pow4(s2t))) + Tbeta*pow2(s2t)*(8*(lmMst1 - lmMst2)*pow2(Mt)*pow2(
        MuSUSY)*(-1 + pow2(Sbeta)) + pow2(s2t)*pow2(Sbeta)*pow4(Mst2))*pow6(
        Mst1) + (-4*Tbeta*pow2(Mt)*pow2(s2t)*(pow2(MuSUSY)*(-1 + pow2(Sbeta)) +
        2*pow2(Mst2)*pow2(Sbeta)) + pow2(Sbeta)*(16*MuSUSY*s2t*pow3(Mt) - 4*Mt*
        MuSUSY*pow2(Mst2)*pow3(s2t) + 16*Tbeta*pow4(Mt) + Tbeta*pow4(Mst2)*
        pow4(s2t)))*pow6(Mst2)) + Al4p*Mst2*threeLoopFlag*(100*Mt*MuSUSY*pow2(
        Sbeta)*(128*pow2(Msq)*((215 + 37*lmMst2 + 3*(-5 + 8*lmMst2)*lmMt - 18*
        lmMst1*(1 - 2*lmMst2 + lmMt) - 42*pow2(lmMst2))*pow2(Mst2) + pow2(Mst1)
        *(1097 + lmMst2*(1531 - 246*lmMt) - 339*lmMt + (444 - 36*lmMt)*pow2(
        lmMst2) - 18*lmMst1*(39 - 2*lmMst2*(-9 + lmMt) - 7*lmMt + 2*pow2(
        lmMst2)) + 36*pow3(lmMst2)))*pow3(Mt)*pow4(Mst1) - 36*(195*Dmsqst2*
        Mst2*s2t + (Mst2*s2t*(220 - 32*lmMst1*(-2 + lmMst2) + 70*lmMst2 - 59*
        pow2(lmMst2)) + 64*Mt*(3 + 4*lmMst2 - 3*pow2(lmMst2)))*pow2(Msq))*pow2(
        Mst1)*pow2(Mt)*pow4(Mst2) - 288*Mt*pow2(Msq)*pow2(Mst1)*pow2(Mst2)*
        pow2(s2t)*(pow2(Mst1)*pow2(Mst2)*(13 - 125*lmMst2 + 6*lmMst1*(7 + 8*
        lmMst2 - 6*pow2(lmMst2)) - 30*pow2(lmMst2) + 36*pow3(lmMst2)) + 6*(15 -
        11*lmMst2 + lmMst1*(4 + 7*lmMst2 - 6*pow2(lmMst2)) - 7*pow2(lmMst2) +
        6*pow3(lmMst2))*pow4(Mst1) + 6*(-3 - 4*lmMst2 + 3*pow2(lmMst2))*pow4(
        Mst2)) - 3*pow3(Mst2)*pow3(s2t)*(585*Dmsqst2*(2*(lmMst1 - lmMst2)*pow2(
        Mst2)*pow4(Mst1) - pow2(Mst1)*pow4(Mst2) + pow6(Mst1)) - pow2(Msq)*(2*
        pow2(Mst2)*(200 + 1196*lmMst2 - 48*(-2 + lmMst2)*pow2(lmMst1) + 306*
        pow2(lmMst2) + 3*lmMst1*(-492 + 10*lmMst2 + 123*pow2(lmMst2)) - 321*
        pow3(lmMst2))*pow4(Mst1) + 3*(444 - 32*lmMst1*(-2 + lmMst2) + 70*lmMst2
        - 347*pow2(lmMst2))*pow2(Mst1)*pow4(Mst2) + (476 + 1790*lmMst2 - 768*(-
        2 + lmMst2)*pow2(lmMst1) + 1617*pow2(lmMst2) + 96*lmMst1*(-13 - 33*
        lmMst2 + 16*pow2(lmMst2)) - 768*pow3(lmMst2))*pow6(Mst1) - 96*(2 +
        lmMst2 - 3*pow2(lmMst2))*pow6(Mst2))) + 12*Mst2*s2t*pow2(Mt)*(585*
        Dmsqst2*pow2(Mst2)*pow4(Mst1) + pow2(Msq)*(pow2(Mst2)*(964 + 306*lmMst2
        + 96*(-2 + lmMst2)*pow2(lmMst1) - 81*pow2(lmMst2) + 96*(lmMst1 + 3*
        lmMst1*lmMst2 - 2*lmMst1*pow2(lmMst2) + pow3(lmMst2)))*pow4(Mst1) + 96*
        (26 + 29*lmMst2 + (-2 + lmMst2)*pow2(lmMst1) + 4*pow2(lmMst2) - 2*
        lmMst1*(8 + lmMst2 + pow2(lmMst2)) + pow3(lmMst2))*pow6(Mst1) + 96*(2 +
        lmMst2 - 3*pow2(lmMst2))*pow6(Mst2)))) + Tbeta*(300*s2t*pow2(Mt)*pow2(
        MuSUSY)*(-3*(Mst2*s2t*(380 - 32*lmMst1*(-2 + lmMst2) + 38*lmMst2 - 251*
        pow2(lmMst2)) + 128*Mt*(3 + 4*lmMst2 - 3*pow2(lmMst2)))*pow2(Msq)*pow2(
        Mst1)*pow4(Mst2) + Mt*pow2(Msq)*(-64*pow2(Mst2)*(5 + 149*lmMst2 + 12*
        pow2(lmMst2) + 6*lmMst1*(-7 - 8*lmMst2 + 6*pow2(lmMst2)) - 36*pow3(
        lmMst2))*pow4(Mst1) + 64*(85 - 215*lmMst2 + lmMst1*(66 + 90*lmMst2 -
        72*pow2(lmMst2)) - 54*pow2(lmMst2) + 72*pow3(lmMst2))*pow6(Mst1)) +
        Mst2*s2t*(585*Dmsqst2*pow2(Mst1)*(2*lmMst1*pow2(Mst1) - 2*lmMst2*pow2(
        Mst1) - pow2(Mst2))*(pow2(Mst1) + pow2(Mst2)) + pow2(Msq)*(pow2(Mst2)*(
        -1540 - 2506*lmMst2 + 96*(-2 + lmMst2)*pow2(lmMst1) + lmMst1*(2760 +
        36*lmMst2 - 738*pow2(lmMst2)) + 141*pow2(lmMst2) + 642*pow3(lmMst2))*
        pow4(Mst1) + 6*(-336 - 716*lmMst2 + 144*(-2 + lmMst2)*pow2(lmMst1) +
        lmMst1*(668 + 534*lmMst2 - 379*pow2(lmMst2)) - 246*pow2(lmMst2) + 235*
        pow3(lmMst2))*pow6(Mst1) + 96*(2 + lmMst2 - 3*pow2(lmMst2))*pow6(Mst2))
        )) + pow2(Sbeta)*(9600*Mt*pow2(Msq)*pow2(Mst1)*pow3(s2t)*pow4(Mst2)*(
        pow2(Mst1)*pow2(Mst2)*(31 - 101*lmMst2 + 6*lmMst1*(7 + 8*lmMst2 - 6*
        pow2(lmMst2)) - 48*pow2(lmMst2) + 36*pow3(lmMst2)) + (77 + 59*lmMst2 -
        6*lmMst1*(3 + lmMst2) - 12*pow2(lmMst2))*pow4(Mst1) + 6*(-3 - 4*lmMst2
        + 3*pow2(lmMst2))*pow4(Mst2)) + 6400*s2t*pow2(Msq)*pow2(Mst1)*pow3(Mt)*
        (-3*(6*pow2(Mst2)*(116 + lmMst2*(181 - 30*lmMt) - 36*lmMt - 4*(-14 +
        lmMt)*pow2(lmMst2) - 2*lmMst1*(42 + lmMst2*(21 - 2*lmMt) - 8*lmMt + 2*
        pow2(lmMst2)) + 4*pow3(lmMst2)) + pow2(MuSUSY)*(85 - 215*lmMst2 +
        lmMst1*(66 + 90*lmMst2 - 72*pow2(lmMst2)) - 54*pow2(lmMst2) + 72*pow3(
        lmMst2)))*pow4(Mst1) + 18*(3 + 4*lmMst2 - 3*pow2(lmMst2))*(2*pow2(Mst2)
        + pow2(MuSUSY))*pow4(Mst2) - pow2(Mst1)*(-3*pow2(Mst2)*pow2(MuSUSY)*(5
        + 149*lmMst2 + 12*pow2(lmMst2) + 6*lmMst1*(-7 - 8*lmMst2 + 6*pow2(
        lmMst2)) - 36*pow3(lmMst2)) + (553 + 194*lmMst2 + 18*lmMst1*(-1 + 4*
        lmMst2 - 2*lmMt) - 30*lmMt + 48*lmMst2*lmMt - 192*pow2(lmMst2))*pow4(
        Mst2))) + 75*pow4(s2t)*pow5(Mst2)*(-585*Dmsqst2*pow2(Mst1)*(-((1 + 2*
        lmMst1 - 2*lmMst2)*pow2(Mst1)*pow2(Mst2)) + (-1 + 2*lmMst1 - 2*lmMst2)*
        pow4(Mst1) + pow4(Mst2)) + pow2(Msq)*(pow2(Mst2)*(932 - 2182*lmMst2 +
        96*(-2 + lmMst2)*pow2(lmMst1) - 1653*pow2(lmMst2) - 6*lmMst1*(-524 +
        26*lmMst2 + 123*pow2(lmMst2)) + 642*pow3(lmMst2))*pow4(Mst1) + 3*(-508
        + 32*lmMst1*(-2 + lmMst2) - 102*lmMst2 + 443*pow2(lmMst2))*pow2(Mst1)*
        pow4(Mst2) + (-76 + 602*lmMst2 + 672*(-2 + lmMst2)*pow2(lmMst1) - 1005*
        pow2(lmMst2) - 6*lmMst1*(284 - 538*lmMst2 + 133*pow2(lmMst2)) + 126*
        pow3(lmMst2))*pow6(Mst1) + 96*(2 + lmMst2 - 3*pow2(lmMst2))*pow6(Mst2))
        ) - 16*Mst2*pow4(Mt)*(43875*Dmsqst2*pow2(Mst1)*pow2(Mst2)*(pow2(Mst1) +
        pow2(Mst2)) - pow2(Msq)*((394748 - 14400*lmMst1*(8 + 2*lmMst2 - lmMt) -
        102480*lmMt - 6*lmMst2*(-50363 + 4680*lmMt) + 41355*pow2(lmMst2))*pow2(
        Mst2)*pow4(Mst1) + 225*(-188 + 32*lmMst1*(-2 + lmMst2) - 166*lmMst2 +
        59*pow2(lmMst2))*pow2(Mst1)*pow4(Mst2) + 7200*(283 + lmMst2*(489 - 93*
        lmMt) - 94*lmMt + 9*(21 - 2*lmMt)*pow2(lmMst2) - 2*lmMst1*(115 - 9*
        lmMst2*(-8 + lmMt) - 24*lmMt + 9*pow2(lmMst2)) + 18*pow3(lmMst2))*pow6(
        Mst1) + 7200*(2 + lmMst2 - 3*pow2(lmMst2))*pow6(Mst2))) + 300*Mst2*
        pow2(Mt)*pow2(s2t)*(585*Dmsqst2*(pow4(Mst1)*((1 - 2*lmMst1 + 2*lmMst2)*
        pow2(Mst2)*pow2(MuSUSY) - 4*pow4(Mst2)) + pow2(Mst1)*(2*pow2(Mst2) +
        pow2(MuSUSY))*pow4(Mst2) + 2*(pow2(Mst2) + (-lmMst1 + lmMst2)*pow2(
        MuSUSY))*pow6(Mst1)) - pow2(Msq)*(pow4(Mst1)*(pow2(Mst2)*pow2(MuSUSY)*(
        -1540 - 2506*lmMst2 + 96*(-2 + lmMst2)*pow2(lmMst1) + lmMst1*(2760 +
        36*lmMst2 - 738*pow2(lmMst2)) + 141*pow2(lmMst2) + 642*pow3(lmMst2)) +
        4*(812 + 258*lmMst2 + 48*(-2 + lmMst2)*pow2(lmMst1) - 129*pow2(lmMst2)
        + 48*(lmMst1*(3 + 2*lmMst2 - 2*pow2(lmMst2)) + pow3(lmMst2)))*pow4(
        Mst2)) + 2*((1532 + 2478*lmMst2 - 96*lmMst1*(17 + 5*lmMst2) + 465*pow2(
        lmMst2))*pow2(Mst2) + 3*pow2(MuSUSY)*(-336 - 716*lmMst2 + 144*(-2 +
        lmMst2)*pow2(lmMst1) + lmMst1*(668 + 534*lmMst2 - 379*pow2(lmMst2)) -
        246*pow2(lmMst2) + 235*pow3(lmMst2)))*pow6(Mst1) + 96*(2 + lmMst2 - 3*
        pow2(lmMst2))*(2*pow2(Mst2) + pow2(MuSUSY))*pow6(Mst2) - 3*pow2(Mst1)*(
        (380 - 32*lmMst1*(-2 + lmMst2) + 38*lmMst2 - 251*pow2(lmMst2))*pow2(
        MuSUSY)*pow4(Mst2) + (568 - 64*lmMst1*(-2 + lmMst2) + 204*lmMst2 - 310*
        pow2(lmMst2))*pow6(Mst2))))))))))/(16200.*Tbeta*pow2(Sbeta)*pow4(Msq)*
        pow4(Mst1)*pow6(Mst2)) + (twoLoopFlag*((-72*Mt*pow3(s2t)*(-(Mst2*(4*(3
        + lmMst1*(2 + lmMst2) - pow2(lmMst2))*pow2(Mst1)*pow2(Mst2) + (3 + 2*
        lmMst1 - 2*lmMst2)*pow4(Mst1) + 4*(-4 + lmMst1 - 3*lmMst2 + lmMst1*
        lmMst2 - pow2(lmMst2))*pow4(Mst2))) + Dmglst2*(-4*(1 + lmMst1*(-2 +
        lmMst2) + 4*lmMst2 - pow2(lmMst2))*pow2(Mst1)*pow2(Mst2) + (1 + 6*
        lmMst1 - 6*lmMst2)*pow4(Mst1) + 4*(6 + lmMst1 + lmMst2 - lmMst1*lmMst2
        + pow2(lmMst2))*pow4(Mst2))))/pow2(Mst2) - (144*s2t*pow3(Mt)*(4*pow2(
        Mst1)*((4 - 5*lmMst2 + lmMst1*(5 + 2*lmMst2 - 2*lmMt) + 2*lmMst2*lmMt -
        2*pow2(lmMst2))*pow2(Mst2) + (-5 - 6*lmMst2 + lmMst1*(4 + 3*lmMst2) -
        3*pow2(lmMst2))*pow2(MuSUSY))*pow3(Mst2) + Mst2*(4*(5 - 7*lmMst2 +
        lmMst1*(7 + 2*lmMst2 - 2*lmMt) + 2*lmMst2*lmMt - 2*pow2(lmMst2))*pow2(
        Mst2) + (-21 - 38*lmMst2 + 10*lmMst1*(3 + 2*lmMst2) - 20*pow2(lmMst2))*
        pow2(MuSUSY))*pow4(Mst1) + 4*(-4 + lmMst1 - 3*lmMst2 + lmMst1*lmMst2 -
        pow2(lmMst2))*pow2(MuSUSY)*pow5(Mst2) + Dmglst2*(4*pow2(Mst1)*pow2(
        Mst2)*((6 + 13*lmMst2 - 4*lmMt - 2*lmMst2*lmMt + lmMst1*(-9 - 2*lmMst2
        + 2*lmMt) + 2*pow2(lmMst2))*pow2(Mst2) + (-11 + 2*lmMst2 + lmMst1*(-4 +
        3*lmMst2) - 3*pow2(lmMst2))*pow2(MuSUSY)) - (4*(1 - 29*lmMst2 + lmMst1*
        (25 + 6*lmMst2 - 6*lmMt) + 4*lmMt + 6*lmMst2*lmMt - 6*pow2(lmMst2))*
        pow2(Mst2) + (65 + lmMst1*(34 - 20*lmMst2) - 26*lmMst2 + 20*pow2(
        lmMst2))*pow2(MuSUSY))*pow4(Mst1) - 4*(6 + lmMst1 + lmMst2 - lmMst1*
        lmMst2 + pow2(lmMst2))*pow2(MuSUSY)*pow4(Mst2) + 8*(5 + lmMst1*(-1 +
        lmMst2) + lmMst2 - pow2(lmMst2))*pow6(Mst2)) + 4*(3 - 4*lmMst2 + 2*(
        lmMst1 + lmMst1*lmMst2 + lmMt) - 2*pow2(lmMst2))*pow7(Mst2)))/pow6(
        Mst2) + (32*pow4(Mt)*(2*Dmglst2*((53 - 18*lmMst1 + 24*lmMst2 - 24*lmMt)
        *pow2(Mst1)*pow4(Mst2) + 18*((5 + lmMst2*(11 - 2*lmMt) - 2*lmMst1*(4 +
        lmMst2 - lmMt) - 3*lmMt + 2*pow2(lmMst2))*pow2(Mst2)*pow4(Mst1) + (3 -
        6*lmMst2*(-5 + lmMt) - 5*lmMt + lmMst1*(-25 - 6*lmMst2 + 6*lmMt) + 6*
        pow2(lmMst2))*pow6(Mst1)) - 18*lmMst2*pow6(Mst2)) + 9*(2*(2 + 2*lmMst1*
        (3 + lmMst2 - lmMt) + lmMt + lmMst2*(-7 + 2*lmMt) - 2*pow2(lmMst2))*
        pow3(Mst2)*pow4(Mst1) + (3 + 2*lmMst1*(3 + lmMst2) + lmMt + lmMst2*(-5
        + 4*lmMt) + pow2(lmMst1) - pow2(lmMst2) - 6*pow2(lmMt))*pow2(Mst1)*
        pow5(Mst2) + Mst2*(9 + lmMst1*(22 + 6*lmMst2 - 6*lmMt) + 6*lmMst2*(-4 +
        lmMt) + 2*lmMt - 6*pow2(lmMst2))*pow6(Mst1) - 2*(1 + lmMst2)*pow7(Mst2)
        )))/(pow2(Mst1)*pow5(Mst2)) + (36*Mt*MuSUSY*(16*pow3(Mt)*(2*(4 +
        lmMst1*(3 + 2*lmMst2 - lmMt) + lmMst2*(-4 + lmMt) + lmMt - 2*pow2(
        lmMst2))*pow2(Mst1)*pow3(Mst2) + Mst2*(13 + 2*lmMst1*(6 + 3*lmMst2 - 2*
        lmMt) + 2*lmMt + 2*lmMst2*(-7 + 2*lmMt) - 6*pow2(lmMst2))*pow4(Mst1) +
        Dmglst2*(11 - 6*lmMst1*lmMst2 - 8*lmMst2*(-5 + lmMt) + 8*lmMst1*(-4 +
        lmMt) - 8*lmMt + 6*pow2(lmMst2))*pow4(Mst1) + 2*Dmglst2*((7 + 7*lmMst2
        + lmMst1*(-5 + lmMt) - 2*lmMt - lmMst2*lmMt)*pow2(Mst1)*pow2(Mst2) + (5
        + lmMst1*(-1 + lmMst2) + lmMst2 - pow2(lmMst2))*pow4(Mst2)) + 2*(2 +
        lmMst1 + (-2 + lmMst1)*lmMst2 + lmMt - pow2(lmMst2))*pow5(Mst2)) + 6*
        Mt*pow2(Mst2)*pow2(s2t)*(4*(1 + 3*lmMst2 - lmMst1*(3 + 2*lmMst2) + 2*
        pow2(lmMst2))*pow2(Mst1)*pow3(Mst2) + Mst2*(1 + 14*lmMst2 - 2*lmMst1*(7
        + 4*lmMst2) + 8*pow2(lmMst2))*pow4(Mst1) + Dmglst2*(4*(5 + lmMst1*(3 -
        2*lmMst2) - 3*lmMst2 + 2*pow2(lmMst2))*pow2(Mst1)*pow2(Mst2) + (21 +
        lmMst1*(18 - 8*lmMst2) - 18*lmMst2 + 8*pow2(lmMst2))*pow4(Mst1) + 4*(6
        + lmMst1 + lmMst2 - lmMst1*lmMst2 + pow2(lmMst2))*pow4(Mst2)) + 4*(4 +
        3*lmMst2 - lmMst1*(1 + lmMst2) + pow2(lmMst2))*pow5(Mst2)) + (8*Mst2*
        s2t*pow2(Mt)*(Dmglst2*(4*(1 - lmMst1 + lmMst2)*pow2(Mst2)*pow4(Mst1) +
        2*(1 + 2*lmMst2)*pow2(Mst1)*pow4(Mst2) + (4 - 8*lmMst1 + 8*lmMst2)*
        pow6(Mst1) - 4*lmMst2*pow6(Mst2)) - Mst2*(-2*lmMst1*((1 + 2*lmMst2)*(
        pow2(Mst1) + pow2(Mst2))*pow4(Mst1) + (3 + lmMst2)*pow2(Mst1)*pow4(
        Mst2)) + pow2(lmMst1)*(2*pow2(Mst2)*pow4(Mst1) - pow2(Mst1)*pow4(Mst2)
        + 2*pow6(Mst1)) + pow2(lmMst2)*(2*pow2(Mst2)*pow4(Mst1) + 3*pow2(Mst1)*
        pow4(Mst2) + 2*pow6(Mst1)) + 2*pow6(Mst2) + 2*lmMst2*(pow2(Mst2)*pow4(
        Mst1) + 2*pow2(Mst1)*pow4(Mst2) + pow6(Mst1) + pow6(Mst2)))))/pow2(
        Mst1) + (pow3(Mst2)*pow3(s2t)*(2*(-16 + 6*lmMst2 - 2*lmMst1*(8 + 5*
        lmMst2) + 3*pow2(lmMst1) + 7*pow2(lmMst2))*pow3(Mst2)*pow4(Mst1) - 2*(-
        12 - 18*lmMst2 + 2*lmMst1*(5 + 3*lmMst2) + pow2(lmMst1) - 7*pow2(
        lmMst2))*pow2(Mst1)*pow5(Mst2) + Mst2*(3 + lmMst1*(2 - 32*lmMst2) - 2*
        lmMst2 + 16*pow2(lmMst1) + 16*pow2(lmMst2))*pow6(Mst1) - 4*Dmglst2*(2*(
        1 + 2*lmMst1 - lmMst2)*pow2(Mst2)*pow4(Mst1) + 2*(1 + lmMst1 - lmMst2 +
        2*lmMst1*lmMst2 - 2*pow2(lmMst2))*pow2(Mst1)*pow4(Mst2) + (1 + 2*lmMst1
        - 2*lmMst2)*pow6(Mst1) - 2*lmMst2*pow6(Mst2)) + 4*(1 + lmMst2)*pow7(
        Mst2)))/pow2(Mst1)))/(Tbeta*pow6(Mst2)) + 36*pow2(Mt)*pow2(s2t)*(4*(-2*
        lmMst1*(-2 + lmMst2) - lmMst2*(2 + lmMst2) + 3*pow2(lmMst1))*pow2(Mst1)
        - 4*(2 - 2*lmMst2 + 2*lmMst1*(3 + lmMst2) + pow2(lmMst1) - 3*pow2(
        lmMst2))*pow2(Mst2) + (8*(2*Dmglst2*lmMst2 + Mst2 + lmMst2*Mst2)*pow3(
        Mst2))/pow2(Mst1) - (8*Dmglst2*(pow2(Mst1)*pow2(Mst2) - 2*lmMst1*pow2(
        Mst1)*(pow2(Mst1) + pow2(Mst2)) + 2*lmMst2*pow4(Mst1) + pow4(Mst2) + 4*
        lmMst2*pow4(Mst2)))/pow3(Mst2) + (pow2(MuSUSY)*(4*(-1 - 13*lmMst1 + 13*
        lmMst2 - 8*lmMst1*lmMst2 + pow2(lmMst1) + 7*pow2(lmMst2))*pow3(Mst2)*
        pow4(Mst1) - 2*(-14 + 10*lmMst1 - 20*lmMst2 + 6*lmMst1*lmMst2 + pow2(
        lmMst1) - 7*pow2(lmMst2))*pow2(Mst1)*pow5(Mst2) + Mst2*(-1 - 50*lmMst1
        + 50*lmMst2 - 64*lmMst1*lmMst2 + 20*pow2(lmMst1) + 44*pow2(lmMst2))*
        pow6(Mst1) - 4*Dmglst2*(2*(2 - 3*lmMst2 + lmMst1*(3 + 2*lmMst2) - 2*
        pow2(lmMst2))*pow2(Mst2)*pow4(Mst1) + 2*(1 + lmMst1 - 2*lmMst2 + 2*
        lmMst1*lmMst2 - 2*pow2(lmMst2))*pow2(Mst1)*pow4(Mst2) + (5 - 8*lmMst2 +
        4*lmMst1*(2 + lmMst2) - 4*pow2(lmMst2))*pow6(Mst1) - 2*lmMst2*pow6(
        Mst2)) + 4*(1 + lmMst2)*pow7(Mst2)))/(pow2(Mst1)*pow5(Mst2))) + ((-9*
        pow4(s2t)*(4*(-14 - 3*lmMst1 - 6*lmMst2 - 2*lmMst1*lmMst2 + 2*pow2(
        lmMst1))*pow3(Mst2)*pow4(Mst1) - 2*(-10 + 10*lmMst1 - 16*lmMst2 + 6*
        lmMst1*lmMst2 + pow2(lmMst1) - 7*pow2(lmMst2))*pow2(Mst1)*pow5(Mst2) +
        Mst2*(35 + 34*lmMst1 - 14*lmMst2 - 12*lmMst1*lmMst2 + 10*pow2(lmMst1) +
        2*pow2(lmMst2))*pow6(Mst1) + 4*Dmglst2*(-2*(lmMst1 - 2*lmMst1*lmMst2 +
        2*pow2(lmMst2))*pow2(Mst2)*pow4(Mst1) - 2*(1 + lmMst1 + 2*lmMst1*lmMst2
        - 2*pow2(lmMst2))*pow2(Mst1)*pow4(Mst2) + (1 + 2*lmMst1)*pow6(Mst1) +
        2*lmMst2*pow6(Mst2)) + 4*(1 + lmMst2)*pow7(Mst2)))/Mst2 - (36*s2t*pow2(
        Mt)*pow2(MuSUSY)*(4*(4*Mt*(5 + 6*lmMst2 - lmMst1*(4 + 3*lmMst2) + 3*
        pow2(lmMst2)) + Mst2*s2t*(-1 + 13*lmMst2 - lmMst1*(13 + 8*lmMst2) +
        pow2(lmMst1) + 7*pow2(lmMst2)))*pow3(Mst2)*pow4(Mst1) + 2*(-(Mst2*s2t*(
        -14 - 20*lmMst2 + 2*lmMst1*(5 + 3*lmMst2) + pow2(lmMst1) - 7*pow2(
        lmMst2))) + 8*Mt*(4 + 3*lmMst2 - lmMst1*(1 + lmMst2) + pow2(lmMst2)))*
        pow2(Mst1)*pow5(Mst2) + 4*Dmglst2*(Mst2*s2t*(-5 + 8*lmMst2 - 4*lmMst1*(
        2 + lmMst2) + 4*pow2(lmMst2)) + Mt*(65 + lmMst1*(34 - 20*lmMst2) - 26*
        lmMst2 + 20*pow2(lmMst2)))*pow6(Mst1) + Mst2*(Mst2*s2t*(-1 + 50*lmMst2
        - 2*lmMst1*(25 + 32*lmMst2) + 20*pow2(lmMst1) + 44*pow2(lmMst2)) + Mt*(
        84 + 152*lmMst2 - 40*lmMst1*(3 + 2*lmMst2) + 80*pow2(lmMst2)))*pow6(
        Mst1) + 8*Dmglst2*((Mst2*s2t*(-2 + 3*lmMst2 - lmMst1*(3 + 2*lmMst2) +
        2*pow2(lmMst2)) + Mt*(22 + lmMst1*(8 - 6*lmMst2) - 4*lmMst2 + 6*pow2(
        lmMst2)))*pow2(Mst2)*pow4(Mst1) + (-(Mst2*s2t*(1 + lmMst1 - 2*lmMst2 +
        2*lmMst1*lmMst2 - 2*pow2(lmMst2))) + 2*Mt*(6 + lmMst1 + lmMst2 -
        lmMst1*lmMst2 + pow2(lmMst2)))*pow2(Mst1)*pow4(Mst2) + lmMst2*s2t*pow7(
        Mst2)) + 4*(1 + lmMst2)*s2t*pow8(Mst2)))/(pow2(Sbeta)*pow6(Mst2)))/
        pow2(Mst1)))/216. + xDmglst2*pow2(Dmglst2)*(-(twoLoopFlag*(72*Mt*pow3(
        s2t)*pow4(Mst2)*(8*pow2(Mst1)*pow2(Mst2) + (3 - 6*lmMst1 + 6*lmMst2)*
        pow4(Mst1) + 2*(8 + lmMst1 - lmMst2)*pow4(Mst2)) + (96*Mst2*pow4(Mt)*(
        3*(37 + 41*lmMst2 - 13*lmMt - 6*lmMst2*lmMt + lmMst1*(-28 - 6*lmMst2 +
        6*lmMt) + 6*pow2(lmMst2))*pow2(Mst2)*pow4(Mst1) + (17 - 3*lmMst1 + 9*
        lmMst2 - 3*lmMt)*pow2(Mst1)*pow4(Mst2) + 3*(75 + 174*lmMst2 - 37*lmMt -
        30*lmMst2*lmMt + lmMst1*(-137 - 30*lmMst2 + 30*lmMt) + 30*pow2(lmMst2))
        *pow6(Mst1) + 3*(-2 + lmMst2)*pow6(Mst2)))/pow2(Mst1) - 16*s2t*pow3(Mt)
        *(9*(2*(27 + 70*lmMst2 - 4*lmMst1*(14 + 3*lmMst2 - 3*lmMt) - 2*(7 + 6*
        lmMst2)*lmMt + 12*pow2(lmMst2))*pow2(Mst2) + 67*pow2(MuSUSY))*pow4(
        Mst1) + 18*(8 + lmMst1 - lmMst2)*pow2(MuSUSY)*pow4(Mst2) + 18*pow2(
        Mst1)*(2*(10 + lmMst1 - lmMst2)*pow2(Mst2)*pow2(MuSUSY) + (19 + 17*
        lmMst2 - 6*lmMt - 2*lmMst2*lmMt + lmMst1*(-11 - 2*lmMst2 + 2*lmMt) + 2*
        pow2(lmMst2))*pow4(Mst2)) + (19 + 36*lmMst1 - 42*lmMst2 + 6*lmMt)*pow6(
        Mst2)) - (12*Mst2*pow2(Mt)*pow2(s2t)*(-2*pow4(Mst1)*(3*(8 - 9*lmMst1 +
        9*lmMst2 + 2*lmMst1*lmMst2 - 2*pow2(lmMst2))*pow2(Mst2)*pow2(MuSUSY) +
        (-13 + 18*lmMst1 - 24*lmMst2)*pow4(Mst2)) - 3*(4*(-4 + 7*lmMst1 - 7*
        lmMst2)*pow2(Mst2) + (17 + 4*lmMst1*(-7 + lmMst2) + 28*lmMst2 - 4*pow2(
        lmMst2))*pow2(MuSUSY))*pow6(Mst1) + 6*(-2 + lmMst2)*(2*pow2(Mst2) +
        pow2(MuSUSY))*pow6(Mst2) - 2*pow2(Mst1)*(3*(8 - 5*lmMst1 + 4*lmMst2 +
        2*lmMst1*lmMst2 - 2*pow2(lmMst2))*pow2(MuSUSY)*pow4(Mst2) + (-29 + 12*
        lmMst2)*pow6(Mst2))))/pow2(Mst1) + (-9*pow4(s2t)*pow5(Mst2)*(2*(-6 +
        lmMst1 - 2*lmMst1*lmMst2 + 2*pow2(lmMst2))*pow2(Mst2)*pow4(Mst1) + 2*(4
        + 6*lmMst2 + lmMst1*(-5 + 2*lmMst2) - 2*pow2(lmMst2))*pow2(Mst1)*pow4(
        Mst2) + (1 - 2*lmMst1)*pow6(Mst1) - 2*(-2 + lmMst2)*pow6(Mst2)) + (36*
        s2t*pow2(Mt)*pow2(MuSUSY)*(-2*(-8*(10 + lmMst1 - lmMst2)*Mt + Mst2*s2t*
        (8 - 9*lmMst1 + 9*lmMst2 + 2*lmMst1*lmMst2 - 2*pow2(lmMst2)))*pow2(
        Mst2)*pow4(Mst1) - 2*(-4*(8 + lmMst1 - lmMst2)*Mt + Mst2*s2t*(8 - 5*
        lmMst1 + 4*lmMst2 + 2*lmMst1*lmMst2 - 2*pow2(lmMst2)))*pow2(Mst1)*pow4(
        Mst2) + (268*Mt + Mst2*s2t*(-17 - 4*lmMst1*(-7 + lmMst2) - 28*lmMst2 +
        4*pow2(lmMst2)))*pow6(Mst1) + 2*(-2 + lmMst2)*s2t*pow7(Mst2)))/pow2(
        Sbeta) - (4*Mt*MuSUSY*(-2*(pow2(Mst2)*(36*(5 - 3*lmMst1 + 3*lmMst2)*
        Mst2*s2t*pow2(Mt) - 54*(12 + lmMst1 - lmMst2)*Mt*pow2(Mst2)*pow2(s2t) +
        4*(155 + 3*lmMst2*(41 - 6*lmMt) - 18*lmMst1*(4 + lmMst2 - lmMt) - 51*
        lmMt + 18*pow2(lmMst2))*pow3(Mt) - 9*(4*lmMst1 - 5*lmMst2)*pow3(Mst2)*
        pow3(s2t))*pow4(Mst1) + pow2(Mst1)*(6*(17 - 6*lmMst2)*Mst2*s2t*pow2(Mt)
        - 54*(8 + lmMst1 - lmMst2)*Mt*pow2(Mst2)*pow2(s2t) + 4*(11 + 18*lmMst1
        - 21*lmMst2 + 3*lmMt)*pow3(Mt) + 9*(6 + 5*lmMst2 + lmMst1*(-5 + 2*
        lmMst2) - 2*pow2(lmMst2))*pow3(Mst2)*pow3(s2t))*pow4(Mst2)) - (72*(9 -
        10*lmMst1 + 10*lmMst2)*Mst2*s2t*pow2(Mt) + 54*(-27 + 4*lmMst1 - 4*
        lmMst2)*Mt*pow2(Mst2)*pow2(s2t) + 8*(335 + 699*lmMst2 - 18*lmMst1*(29 +
        7*lmMst2 - 7*lmMt) - 3*(59 + 42*lmMst2)*lmMt + 126*pow2(lmMst2))*pow3(
        Mt) + 9*(1 - 10*lmMst1 + 10*lmMst2)*pow3(Mst2)*pow3(s2t))*pow6(Mst1) +
        18*(-2 + lmMst2)*s2t*(-4*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow7(Mst2)))/
        Tbeta)/pow2(Mst1)))/(108.*pow7(Mst2)) + (Al4p*threeLoopFlag*((175*Mt*
        pow3(s2t)*(12*pow2(Mst1)*pow2(Mst2)*(282298 - 61560*B4 + 1620*DN -
        236520*lmMst1 - 180*(-439 + 438*lmMst1)*lmMst2 + 2240*OepS2 - 81*(-1373
        + 560*lmMst1 - 560*lmMst2)*S2 - 4320*pow2(lmMst1) + 1080*(77 + 48*
        lmMst1)*pow2(lmMst2) - 51840*pow3(lmMst2)) - (2727217 - 437920*OepS2 +
        360*lmMst2*(4958 - 24633*S2) + 3648132*S2 + 360*lmMst1*(-1460 + 1980*
        lmMst2 + 24633*S2) - 349920*pow2(lmMst1) - 673920*pow2(lmMst2))*pow4(
        Mst1) + 103680*(3 + 4*lmMst2 - 3*pow2(lmMst2))*pow4(Mst2)))/(Mst2*pow2(
        Mst1)) - (50*s2t*pow3(Mt)*((8*pow2(Mst2)*(7384247 - 12322800*lmMt +
        1183840*OepS2 + 5821092*S2 + 22680*(41 + 4*lmMst2 - 4*lmMt)*pow2(
        lmMst1) + 1890*(-2489 + 48*lmMt)*pow2(lmMst2) + 315*lmMst1*(1073 +
        7734*lmMst2 + 690*lmMt - 76104*S2 + 864*pow2(lmMst2) - 864*pow2(lmMt))
        + 816480*pow2(lmMt) + 315*lmMst2*(11587 + 966*lmMt + 76104*S2 + 864*
        pow2(lmMt)) - 362880*pow3(lmMst2)) + 35*pow2(MuSUSY)*(996211 - 295488*
        B4 + 7776*DN - 1030176*lmMst1 - 144*(-1883 + 3618*lmMst1)*lmMst2 +
        98336*OepS2 - 756*(259 + 2634*lmMst1 - 2634*lmMst2)*S2 + 49248*pow2(
        lmMst1) + 5184*(67 + 48*lmMst1)*pow2(lmMst2) - 248832*pow3(lmMst2)))*
        pow4(Mst1) + 725760*(3 + 4*lmMst2 - 3*pow2(lmMst2))*(2*pow2(Mst2) +
        pow2(MuSUSY))*pow4(Mst2) + 12*pow2(Mst1)*(7*pow2(Mst2)*pow2(MuSUSY)*(
        334138 - 61560*B4 + 1620*DN - 236520*lmMst1 - 180*(-823 + 438*lmMst1)*
        lmMst2 + 2240*OepS2 - 81*(-1373 + 560*lmMst1 - 560*lmMst2)*S2 - 4320*
        pow2(lmMst1) + 1080*(29 + 48*lmMst1)*pow2(lmMst2) - 51840*pow3(lmMst2))
        - 2*(1243547 + 1157940*lmMt - 15680*OepS2 + 1680*lmMst2*(-1091 + 42*
        lmMt - 189*S2) - 953208*S2 + 7560*lmMst1*(135 + 33*lmMst2 - 8*lmMt +
        42*S2) + 30240*pow2(lmMst1) - 486360*pow2(lmMst2) + 15120*pow2(lmMt))*
        pow4(Mst2))))/(pow2(Mst1)*pow5(Mst2)) + 2551500*pow4(s2t)*((-2*OepS2*(
        1136*pow2(Mst1) + 21*pow2(Mst2)))/729. + (S2*(4*(14023 + 8520*lmMst1 -
        8520*lmMst2)*pow2(Mst1) + 9*(65819 + 70*lmMst1 - 70*lmMst2)*pow2(Mst2))
        )/540. + pow2(Mst1)*(29.427737997256514 - B4/3. + (4*D3)/9. - (5*DN)/
        18. + (270961*lmMst1)/8100. + (653*pow2(lmMst1))/90. - (lmMst2*(332311
        + 189090*lmMst1 + 57150*pow2(lmMst1)))/8100. + (7.95 + (74*lmMst1)/9.)*
        pow2(lmMst2) + (5*pow3(lmMst1))/18. - (13*pow3(lmMst2))/9.) - pow2(
        Mst2)*(47.56630658436214 - (19*B4)/9. + (20*D3)/9. - (7*DN)/6. + (163*
        lmMst1)/9. - (25*pow2(lmMst1))/36. + (2*lmMst2*(-347 - 18*lmMst1 + 12*
        pow2(lmMst1)))/27. - ((99 + 41*lmMst1)*pow2(lmMst2))/6. + (107*pow3(
        lmMst2))/18.) + ((940 - 32*lmMst1*(-2 + lmMst2) + 102*lmMst2 - 443*
        pow2(lmMst2))*pow4(Mst2))/(36.*pow2(Mst1)) - (8*(2 + lmMst2 - 3*pow2(
        lmMst2))*pow6(Mst2))/(9.*pow4(Mst1))) + (pow2(Mt)*pow2(s2t)*(-75*pow4(
        Mst1)*(-7*pow2(Mst2)*pow2(MuSUSY)*(13169 - 41040*B4 + 43200*D3 - 22680*
        DN + 282960*lmMst1 + 1120*OepS2 - 324*(65819 + 70*lmMst1 - 70*lmMst2)*
        S2 - 13500*pow2(lmMst1) + 720*lmMst2*(-775 + 12*lmMst1 + 24*pow2(
        lmMst1)) - 1080*(-2 + 123*lmMst1)*pow2(lmMst2) + 115560*pow3(lmMst2)) +
        18*(2164753 - 3360*B4 + 3360*D3 - 1680*DN + 89040*lmMst1 - 6720*lmMt +
        7840*OepS2 - 324*(523 + 490*lmMst1 - 490*lmMst2)*S2 + 17220*pow2(
        lmMst1) - 140*lmMst2*(-6485 - 1098*lmMst1 + 96*pow2(lmMst1)) - 420*(129
        + 32*lmMst1)*pow2(lmMst2) + 26880*pow3(lmMst2))*pow4(Mst2)) - 2*(14*
        pow2(MuSUSY)*(2011073 + 1417500*B4 - 1458000*D3 + 749250*DN + 934245*
        lmMst1 - 1178000*OepS2 + 1350*(620417 + 17670*lmMst1 - 17670*lmMst2)*S2
        + 3150900*pow2(lmMst1) - 45*lmMst2*(-124139 + 189090*lmMst1 + 71550*
        pow2(lmMst1)) + 4050*(1323 + 1970*lmMst1)*pow2(lmMst2) + 101250*pow3(
        lmMst1) - 4860000*pow3(lmMst2)) + 9*pow2(Mst2)*(68526791 + 9072000*lmMt
        + 2772000*OepS2 - 107403300*S2 + 420*lmMst2*(102713 + 6000*lmMt +
        133650*S2) + 6300*(131 - 30*lmMst2)*pow2(lmMst1) - 35116200*pow2(
        lmMst2) - 840*lmMst1*(-47431 - 41010*lmMst2 + 3000*lmMt + 66825*S2 +
        6975*pow2(lmMst2)) + 63000*pow3(lmMst1) + 5985000*pow3(lmMst2)))*pow6(
        Mst1) + 9072000*(2 + lmMst2 - 3*pow2(lmMst2))*(2*pow2(Mst2) + pow2(
        MuSUSY))*pow6(Mst2) - 283500*pow2(Mst1)*((812 - 32*lmMst1*(-2 + lmMst2)
        + 38*lmMst2 - 251*pow2(lmMst2))*pow2(MuSUSY)*pow4(Mst2) + (1432 - 64*
        lmMst1*(-2 + lmMst2) + 204*lmMst2 - 310*pow2(lmMst2))*pow6(Mst2))))/(
        pow4(Mst1)*pow4(Mst2)) + (2551500*MuSUSY*(-(pow2(Mt)*pow2(s2t)*((64*
        Mst2*(3 + 4*lmMst2 - 3*pow2(lmMst2)))/(3.*pow2(Mst1)) + (
        761.0320987654321 - 152*B4 + 4*DN - 584*lmMst1 + (4*(631 - 438*lmMst1)*
        lmMst2)/9. - (32*pow2(lmMst1))/3. + (8*(53 + 48*lmMst1)*pow2(lmMst2))/
        3. - 128*pow3(lmMst2))/Mst2 + ((56*OepS2*(415*pow2(Mst1) + 24*pow2(
        Mst2)))/243. - (S2*((21422 + 87150*lmMst1 - 87150*lmMst2)*pow2(Mst1) +
        9*(-1373 + 560*lmMst1 - 560*lmMst2)*pow2(Mst2)))/45.)/pow3(Mst2) + (
        pow2(Mst1)*(199.87633744855967 - 152*B4 + 4*DN - (12848*lmMst1)/27. - (
        8*(293 + 1152*lmMst1)*lmMst2)/27. + (184*pow2(lmMst1))/3. + 8*(35 + 16*
        lmMst1)*pow2(lmMst2) - 128*pow3(lmMst2)))/pow3(Mst2))) + Mt*(pow3(s2t)*
        (92.93189300411522 - (76*B4)/9. + (80*D3)/9. - (14*DN)/3. + (196*
        lmMst1)/3. - (25*pow2(lmMst1))/9. + (2*lmMst2*(-1493 - 24*lmMst1 + 48*
        pow2(lmMst1)))/27. - ((247 + 246*lmMst1)*pow2(lmMst2))/9. + (8*OepS2*(
        21 + (1157*pow2(Mst1))/pow2(Mst2)))/729. + ((-876 + 32*lmMst1*(-2 +
        lmMst2) - 70*lmMst2 + 347*pow2(lmMst2))*pow2(Mst2))/(9.*pow2(Mst1)) - (
        S2*((648463 + 34710*lmMst1 - 34710*lmMst2)*pow2(Mst1) + 9*(65819 + 70*
        lmMst1 - 70*lmMst2)*pow2(Mst2)))/(135.*pow2(Mst2)) - (pow2(Mst1)*(
        24.779058984910836 + (64*B4)/9. - (64*D3)/9. + (32*DN)/9. + (138661*
        lmMst1)/2025. + (159*pow2(lmMst1))/5. - lmMst2*(53.5116049382716 + (
        458*lmMst1)/5. + (286*pow2(lmMst1))/9.) + (2*(1333 + 1355*lmMst1)*pow2(
        lmMst2))/45. + (10*pow3(lmMst1))/9. - (266*pow3(lmMst2))/9.))/pow2(
        Mst2) + (214*pow3(lmMst2))/9. + (32*(2 + lmMst2 - 3*pow2(lmMst2))*pow4(
        Mst2))/(9.*pow4(Mst1))) - (4*T1ep*(21*pow2(Mst2)*(36*Mst2*s2t*pow2(Mt)
        - 24*Mt*pow2(Mst2)*pow2(s2t) + 32*pow3(Mt) + pow3(Mst2)*pow3(s2t)) +
        pow2(Mst1)*(4320*Mst2*s2t*pow2(Mt) - 8715*Mt*pow2(Mst2)*pow2(s2t) +
        17584*pow3(Mt) + 1157*pow3(Mst2)*pow3(s2t))))/(243.*pow5(Mst2))) + (4*(
        -3*(1648637 + 1173060*lmMt - 15680*OepS2 + 1680*lmMst2*(-1214 + 42*lmMt
        - 189*S2) - 953208*S2 + 7560*lmMst1*(131 + 33*lmMst2 - 8*lmMt + 42*S2)
        + 30240*pow2(lmMst1) - 304920*pow2(lmMst2) + 15120*pow2(lmMt))*pow2(
        Mst1)*pow2(Mst2) + 2*(-1657412 - 7512750*lmMt + 615440*OepS2 + 4340358*
        S2 + 11340*(35 + 4*lmMst2 - 4*lmMt)*pow2(lmMst1) + 945*(-1981 + 48*
        lmMt)*pow2(lmMst2) + 315*lmMst1*(-4223 + 2679*lmMst2 + 201*lmMt -
        39564*S2 + 432*pow2(lmMst2) - 432*pow2(lmMt)) + 385560*pow2(lmMt) +
        945*lmMst2*(5321 + 193*lmMt + 13188*S2 + 144*pow2(lmMt)) - 181440*pow3(
        lmMst2))*pow4(Mst1) - 181440*(-3 - 4*lmMst2 + 3*pow2(lmMst2))*pow4(
        Mst2))*pow4(Mt))/(25515.*pow2(Mst1)*pow5(Mst2)) + (s2t*pow3(Mt)*(75*
        pow2(Mst2)*(2435233 - 3360*B4 + 3360*D3 - 1680*DN + 115920*lmMst1 -
        6720*lmMt + 7840*OepS2 - 324*(523 + 490*lmMst1 - 490*lmMst2)*S2 +
        17220*pow2(lmMst1) - 140*lmMst2*(-6695 - 1002*lmMst1 + 96*pow2(lmMst1))
        - 1680*(47 + 8*lmMst1)*pow2(lmMst2) + 26880*pow3(lmMst2))*pow4(Mst1) +
        31500*(652 - 32*lmMst1*(-2 + lmMst2) + 70*lmMst2 - 59*pow2(lmMst2))*
        pow2(Mst1)*pow4(Mst2) + 2*(127852633 - 126000*B4 + 126000*D3 - 63000*DN
        + 23638020*lmMst1 - 252000*(-17 + 5*lmMst1 - 5*lmMst2)*lmMt + 1680000*
        OepS2 - 2700*(22243 + 12600*lmMst1 - 12600*lmMst2)*S2 + 1058400*pow2(
        lmMst1) - 420*lmMst2*(-136544 - 53535*lmMst1 + 1425*pow2(lmMst1)) -
        6300*(3257 + 545*lmMst1)*pow2(lmMst2) + 31500*pow3(lmMst1) + 4000500*
        pow3(lmMst2))*pow6(Mst1) - 1008000*(2 + lmMst2 - 3*pow2(lmMst2))*pow6(
        Mst2)))/(70875.*pow4(Mst1)*pow4(Mst2))))/Tbeta + (4*pow4(Mt)*(15*pow4(
        Mst1)*(840*pow2(Mst2)*pow2(MuSUSY)*(33934 - 90*B4 + 90*D3 - 45*DN +
        720*lmMst1 + 120*(163 + 24*lmMst1)*lmMst2 - 120*lmMt - 10206*S2 - 720*(
        2 + lmMst1)*pow2(lmMst2) + 720*pow3(lmMst2)) - (46632377 + 17278380*
        lmMt - 744800*OepS2 - 29679210*S2 - 168*lmMst2*(-56837 + 7170*lmMt +
        89775*S2) + 132300*pow2(lmMst1) + 10427760*pow2(lmMst2) + 37800*lmMst1*
        (66 - 240*lmMst2 + 399*S2 + 32*(lmMt + pow2(lmMst2))) - 1587600*pow2(
        lmMt) - 1209600*pow3(lmMst2))*pow4(Mst2)) + (pow2(Mst2)*(648916519 -
        1795878000*lmMt + 158732000*OepS2 + 66753450*S2 + 113400*(1403 + 285*
        lmMst2 - 225*lmMt)*pow2(lmMst1) + 94500*(-15839 + 558*lmMt)*pow2(
        lmMst2) + 3780*lmMst1*(132381 + lmMst2*(272035 - 7200*lmMt) - 21415*
        lmMt - 850350*S2 + 22200*pow2(lmMst2) - 21600*pow2(lmMt)) + 176904000*
        pow2(lmMt) + 1260*lmMst2*(-245893 + 169395*lmMt + 2551050*S2 + 64800*
        pow2(lmMt)) - 2268000*pow3(lmMst1) - 113967000*pow3(lmMst2)) + 525*
        pow2(MuSUSY)*(2199511 - 4320*B4 + 4320*D3 - 2160*DN + 140160*lmMst1 +
        960*(1426 + 303*lmMst1)*lmMst2 - 2880*(-16 + 5*lmMst1 - 5*lmMst2)*lmMt
        - 1120*OepS2 + 324*(-3411 + 70*lmMst1 - 70*lmMst2)*S2 - 2880*(77 + 24*
        lmMst1)*pow2(lmMst2) + 69120*pow3(lmMst2)))*pow6(Mst1) + 283500*(620 -
        32*lmMst1*(-2 + lmMst2) + 166*lmMst2 - 59*pow2(lmMst2))*pow2(Mst1)*
        pow6(Mst2) - 9072000*(2 + lmMst2 - 3*pow2(lmMst2))*pow8(Mst2)))/(pow4(
        Mst1)*pow6(Mst2)) + (-10500*T1ep*(4*pow2(Mst1)*(-2*pow2(Mst2)*pow2(Mt)*
        pow2(s2t)*(-589*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 891*pow2(Mst2)*pow2(
        Sbeta)) - 14*Mst2*s2t*(439*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 1208*pow2(
        Mst2)*pow2(Sbeta))*pow3(Mt) + 4*(-21*pow2(MuSUSY)*(-1 + pow2(Sbeta)) +
        5669*pow2(Mst2)*pow2(Sbeta))*pow4(Mt) + (2737*Mt - 284*Mst2*s2t)*pow2(
        Sbeta)*pow3(s2t)*pow5(Mst2)) + 21*pow3(Mst2)*(-4*Mst2*pow2(Mt)*pow2(
        s2t)*(-(pow2(MuSUSY)*(-1 + pow2(Sbeta))) + 18*pow2(Mst2)*pow2(Sbeta)) -
        64*s2t*(pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 2*pow2(Mst2)*pow2(Sbeta))*
        pow3(Mt) + pow2(Sbeta)*(32*Mt*pow3(s2t)*pow4(Mst2) + 304*Mst2*pow4(Mt)
        - pow4(s2t)*pow5(Mst2)))) - (7*pow2(Mt)*pow2(MuSUSY)*(75*pow2(Mst2)*(-
        8*Mst2*Mt*s2t*(334138 - 61560*B4 + 1620*DN - 236520*lmMst1 - 180*(-823
        + 438*lmMst1)*lmMst2 + 2240*OepS2 - 81*(-1373 + 560*lmMst1 - 560*
        lmMst2)*S2 - 4320*pow2(lmMst1) + 1080*(29 + 48*lmMst1)*pow2(lmMst2) -
        51840*pow3(lmMst2)) + 96*pow2(Mt)*(33934 - 90*B4 + 90*D3 - 45*DN + 720*
        lmMst1 + 120*(163 + 24*lmMst1)*lmMst2 - 120*lmMt - 10206*S2 - 720*(2 +
        lmMst1)*pow2(lmMst2) + 720*pow3(lmMst2)) + pow2(Mst2)*pow2(s2t)*(13169
        - 41040*B4 + 43200*D3 - 22680*DN + 282960*lmMst1 + 1120*OepS2 - 324*(
        65819 + 70*lmMst1 - 70*lmMst2)*S2 - 13500*pow2(lmMst1) + 720*lmMst2*(-
        775 + 12*lmMst1 + 24*pow2(lmMst1)) - 1080*(-2 + 123*lmMst1)*pow2(
        lmMst2) + 115560*pow3(lmMst2)))*pow4(Mst1) - 40500*s2t*(Mst2*s2t*(812 -
        32*lmMst1*(-2 + lmMst2) + 38*lmMst2 - 251*pow2(lmMst2)) + 128*Mt*(3 +
        4*lmMst2 - 3*pow2(lmMst2)))*pow2(Mst1)*pow5(Mst2) + 2*(-2*pow2(Mst2)*
        pow2(s2t)*(2011073 + 1417500*B4 - 1458000*D3 + 749250*DN + 934245*
        lmMst1 - 1178000*OepS2 + 1350*(620417 + 17670*lmMst1 - 17670*lmMst2)*S2
        + 3150900*pow2(lmMst1) - 45*lmMst2*(-124139 + 189090*lmMst1 + 71550*
        pow2(lmMst1)) + 4050*(1323 + 1970*lmMst1)*pow2(lmMst2) + 101250*pow3(
        lmMst1) - 4860000*pow3(lmMst2)) - 125*Mst2*Mt*s2t*(996211 - 295488*B4 +
        7776*DN - 1030176*lmMst1 - 144*(-1883 + 3618*lmMst1)*lmMst2 + 98336*
        OepS2 - 756*(259 + 2634*lmMst1 - 2634*lmMst2)*S2 + 49248*pow2(lmMst1) +
        5184*(67 + 48*lmMst1)*pow2(lmMst2) - 248832*pow3(lmMst2)) + 150*pow2(
        Mt)*(2199511 - 4320*B4 + 4320*D3 - 2160*DN + 140160*lmMst1 + 960*(1426
        + 303*lmMst1)*lmMst2 - 2880*(-16 + 5*lmMst1 - 5*lmMst2)*lmMt - 1120*
        OepS2 + 324*(-3411 + 70*lmMst1 - 70*lmMst2)*S2 - 2880*(77 + 24*lmMst1)*
        pow2(lmMst2) + 69120*pow3(lmMst2)))*pow6(Mst1) + 1296000*(2 + lmMst2 -
        3*pow2(lmMst2))*pow2(s2t)*pow8(Mst2)))/pow4(Mst1))/(pow2(Sbeta)*pow6(
        Mst2))))/2.5515e6)) - (Al4p*((35*Al4p*threeLoopFlag*(12*xDmglst2*pow2(
        Dmglst2)*pow4(Msq)*(-2*pow2(Mst1)*(-8*pow2(Mst2)*pow2(Mt)*(-(Tbeta*
        pow2(MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta))*(-229142*z3 - 7489*z4 + 795*
        pow2(z2))) + pow2(Sbeta)*(-72*Mt*MuSUSY*s2t*(2534*z3 + 33*z4 + 90*pow2(
        z2)) - 2*Tbeta*pow2(Mt)*(-619510*z3 + 11338*z4 + 17007*pow2(z2)))) - 4*
        s2t*pow2(Mt)*pow2(Sbeta)*(3*MuSUSY*s2t*(127198*z3 + 6782*z4 - 14613*
        pow2(z2)) + 16*Mt*Tbeta*(-61328*z3 + 2114*z4 + 3171*pow2(z2)))*pow3(
        Mst2) - 8*Mst2*MuSUSY*(-2*Mt*pow2(Sbeta)*(-342281*z3 + 8792*z4 + 13188*
        pow2(z2)) - MuSUSY*s2t*Tbeta*(-1 + pow2(Sbeta))*(-323836*z3 - 8090*z4 +
        37437*pow2(z2)))*pow3(Mt) - 2*Mt*pow2(s2t)*pow2(Sbeta)*(MuSUSY*s2t*(
        547817*z3 + 10924*z4 - 6942*pow2(z2)) + 9*Mt*Tbeta*(-15053*z3 + 792*z4
        + 1188*pow2(z2)))*pow4(Mst2) + Tbeta*(-48*pow2(MuSUSY)*(-1 + pow2(
        Sbeta))*(-25778*z3 + 338*z4 + 21*pow2(z2))*pow4(Mt) + pow2(Sbeta)*(
        Mst2*s2t*(89533*z3 - 4054*z4 - 5352*pow2(z2)) + 28*Mt*(-9920*z3 + 782*
        z4 + 1173*pow2(z2)))*pow3(s2t)*pow5(Mst2))) + 3*pow2(Mst2)*(-4*pow2(
        Mst2)*pow2(Mt)*(-(Tbeta*pow2(MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta))*(
        122917*z3 + 6344*z4 + 1254*pow2(z2))) + 2*Mt*pow2(Sbeta)*(7*Mt*Tbeta*(-
        29483*z3 + 152*z4 + 228*pow2(z2)) + 3*MuSUSY*s2t*(55597*z3 - 264*z4 +
        252*pow2(z2)))) + 4*Mt*s2t*pow2(Sbeta)*(4*Mt*(3*MuSUSY*s2t*(32773*z3 +
        218*z4 - 3804*pow2(z2)) + 2*Mt*Tbeta*(-32323*z3 + 112*z4 + 168*pow2(z2)
        )) + Mst2*s2t*(3*Mt*Tbeta*(55597*z3 - 264*z4 + 252*pow2(z2)) + MuSUSY*
        s2t*(122917*z3 + 6344*z4 + 1254*pow2(z2))))*pow3(Mst2) + 32*Mst2*
        MuSUSY*(Mt*pow2(Sbeta)*(32323*z3 - 112*z4 - 168*pow2(z2)) - MuSUSY*s2t*
        Tbeta*(-1 + pow2(Sbeta))*(-32773*z3 - 218*z4 + 3804*pow2(z2)))*pow3(Mt)
        + Tbeta*(-576*(608*z3 - 9*z4)*pow2(MuSUSY)*(-1 + pow2(Sbeta))*pow4(Mt)
        + pow2(Sbeta)*(-(Mst2*s2t*(122917*z3 + 6344*z4 + 1254*pow2(z2))) + 16*
        Mt*(-32773*z3 - 218*z4 + 3804*pow2(z2)))*pow3(s2t)*pow5(Mst2)))) + 8*
        Dmglst2*(9720*Dmsqst2*z3*(Dmsqst2*xDmsqst2 + pow2(Msq))*pow2(Mst2)*(
        s2t*pow2(Sbeta)*(-2*s2t*pow2(Mt)*(48*MuSUSY + 5*Mst2*Tbeta*(-2 + pow2(
        Sbeta))*pow2(Sbeta)) + 704*Tbeta*pow3(Mt) + 5*Tbeta*pow3(Mst2)*pow3(
        s2t))*pow4(Mst1) - 2*pow2(Mst1)*(16*Mst2*MuSUSY*pow2(Mt)*pow2(s2t)*(-(
        MuSUSY*Tbeta*(-1 + pow2(Sbeta))) + 3*Mst2*pow2(Sbeta)) + 8*s2t*Tbeta*(
        5*pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 28*pow2(Mst2)*pow2(Sbeta))*pow3(Mt)
        - 2*Mt*(5*MuSUSY + 6*Mst2*Tbeta)*pow2(Sbeta)*pow3(Mst2)*pow3(s2t) + 32*
        (10*MuSUSY + 9*Mst2*Tbeta)*pow2(Sbeta)*pow4(Mt) + Tbeta*pow2(Sbeta)*
        pow4(s2t)*pow5(Mst2)) - pow2(Mst2)*(12*Mst2*MuSUSY*pow2(Mt)*pow2(s2t)*(
        -(MuSUSY*Tbeta*(-1 + pow2(Sbeta))) + 2*Mst2*pow2(Sbeta)) + 16*s2t*
        Tbeta*(pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 12*pow2(Mst2)*pow2(Sbeta))*
        pow3(Mt) - 4*Mt*(3*MuSUSY + 2*Mst2*Tbeta)*pow2(Sbeta)*pow3(Mst2)*pow3(
        s2t) + 32*(6*MuSUSY + 5*Mst2*Tbeta)*pow2(Sbeta)*pow4(Mt) + 3*Tbeta*
        pow2(Sbeta)*pow4(s2t)*pow5(Mst2))) + pow4(Msq)*(-(pow4(Mst1)*(-4*Mst2*
        pow2(Mt)*pow2(s2t)*(4*Tbeta*pow2(MuSUSY)*(-1 + pow2(Sbeta))*(353116*z3
        + 11327*z4 + 15897*pow2(z2)) + pow2(Sbeta)*(36*Tbeta*pow2(Mst2)*(-734*
        z3 + 1538*z4 + 2307*pow2(z2)) - 3*Mst2*MuSUSY*(-728116*z3 + 10126*z4 +
        105585*pow2(z2)))) - 8*s2t*(-(Tbeta*pow2(MuSUSY)*(-1 + pow2(Sbeta))*(-
        4*(395581*z3 + 659*z4) + 259215*pow2(z2))) + 4*Mst2*pow2(Sbeta)*(Mst2*
        Tbeta*(653582*z3 - 25730*z4 - 38595*pow2(z2)) - 9*MuSUSY*(2806*z3 +
        1514*z4 + 2757*pow2(z2))))*pow3(Mt) + pow2(Sbeta)*(-2*Mt*(8*MuSUSY*(
        89947*z3 + 6332*z4 + 9498*pow2(z2)) + Mst2*Tbeta*(-565214*z3 + 31142*z4
        + 46713*pow2(z2)))*pow3(Mst2)*pow3(s2t) - 32*(2*Mst2*Tbeta*(-591386*z3
        + 5582*z4 + 8373*pow2(z2)) + MuSUSY*(-834482*z3 + 28538*z4 + 48639*
        pow2(z2)))*pow4(Mt)) + Tbeta*pow2(Sbeta)*(-217283*z3 + 16796*z4 +
        25194*pow2(z2))*pow4(s2t)*pow5(Mst2))) + 18*Mst2*pow2(Mst1)*(pow2(Mst2)
        *pow2(Mt)*(16*Tbeta*pow2(Mt)*pow2(Sbeta)*(-28846*z3 + 322*z4 + 483*
        pow2(z2)) + 3*MuSUSY*pow2(s2t)*(Mst2*pow2(Sbeta)*(99002*z3 + 1210*z4 -
        18273*pow2(z2)) + 8*MuSUSY*Tbeta*(-1 + pow2(Sbeta))*(9747*z3 + 185*z4 +
        237*pow2(z2))) - 12*Mt*s2t*pow2(Sbeta)*(8*MuSUSY*(590*z3 - 4*z4 + 75*
        pow2(z2)) + Mst2*Tbeta*(-26782*z3 + 922*z4 + 1383*pow2(z2)))) + 8*Mst2*
        MuSUSY*(24*Mt*pow2(Sbeta)*(-1675*z3 + 26*z4 + 93*pow2(z2)) - MuSUSY*
        s2t*Tbeta*(-1 + pow2(Sbeta))*(-47456*z3 - 709*z4 + 8535*pow2(z2)))*
        pow3(Mt) + 6*Mt*pow2(s2t)*pow2(Sbeta)*(4*Mt*Tbeta*(741*z3 + 100*z4 +
        150*pow2(z2)) + MuSUSY*s2t*(21373*z3 + 316*z4 + 474*pow2(z2)))*pow4(
        Mst2) + Tbeta*(-48*pow2(MuSUSY)*(-1 + pow2(Sbeta))*(250*z3 - 94*z4 +
        21*pow2(z2))*pow4(Mt) + pow2(Sbeta)*(3*Mst2*s2t*(-1879*z3 + 54*z4) +
        Mt*(-8180*z3 + 416*z4 + 2406*pow2(z2)))*pow3(s2t)*pow5(Mst2))) + 27*
        pow3(Mst2)*(4*Mt*pow2(Mst2)*(pow2(Mst2)*pow2(s2t)*pow2(Sbeta)*(4*Mt*
        Tbeta*(439*z3 - 108*z4) + MuSUSY*s2t*(17615*z3 + 424*z4 + 474*pow2(z2))
        ) + Mt*(Tbeta*pow2(MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta))*(17615*z3 +
        424*z4 + 474*pow2(z2)) + 2*Mt*pow2(Sbeta)*(4*MuSUSY*s2t*(-439*z3 + 108*
        z4) + Mt*Tbeta*(9017*z3 + 56*z4 + 84*pow2(z2))))) - 2*s2t*pow2(Mt)*
        pow2(Sbeta)*(4*Mt*Tbeta*(-18*z3 - 506*z4 + 105*pow2(z2)) + 3*MuSUSY*
        s2t*(-30274*z3 - 542*z4 + 5289*pow2(z2)))*pow3(Mst2) + 4*Mst2*MuSUSY*(
        2*Mt*pow2(Sbeta)*(-18*z3 - 506*z4 + 105*pow2(z2)) - MuSUSY*s2t*Tbeta*(-
        1 + pow2(Sbeta))*(-30274*z3 - 542*z4 + 5289*pow2(z2)))*pow3(Mt) +
        Tbeta*(-64*(83*z3 - 27*z4)*pow2(MuSUSY)*(-1 + pow2(Sbeta))*pow4(Mt) +
        pow2(Sbeta)*(-(Mst2*s2t*(17615*z3 + 424*z4 + 474*pow2(z2))) + 2*Mt*(-
        30274*z3 - 542*z4 + 5289*pow2(z2)))*pow3(s2t)*pow5(Mst2))))) - 3*(5*
        Dmsqst2*pow2(Mst2)*(16*pow2(Msq)*(-2*pow4(Mst1)*(2*pow2(Mt)*pow2(s2t)*(
        8*Mst2*pow2(Sbeta)*(972*MuSUSY*z3 + Mst2*Tbeta*(14*z3 - 2*z4 - 3*pow2(
        z2))) - Tbeta*pow2(MuSUSY)*(-1 + pow2(Sbeta))*(-1771*z3 + 442*z4 + 663*
        pow2(z2))) + pow2(Sbeta)*(4*s2t*(-2592*Mst2*Tbeta*z3 + MuSUSY*(-1930*z3
        + 106*z4 + 159*pow2(z2)))*pow3(Mt) - 208*Mt*MuSUSY*pow2(Mst2)*(2*(-7*z3
        + z4) + 3*pow2(z2))*pow3(s2t) + 15552*Tbeta*z3*pow4(Mt) + Tbeta*(-260*
        z3 + 14*z4 + 21*pow2(z2))*pow4(Mst2)*pow4(s2t))) + 9*Mst2*pow2(Mst1)*(
        288*Mt*z3*(s2t*pow2(Mst2)*(-12*Mt*MuSUSY*s2t + 8*Tbeta*pow2(Mt) +
        Tbeta*pow2(Mst2)*pow2(s2t))*pow2(Sbeta) - 2*MuSUSY*pow2(Mt)*(7*MuSUSY*
        s2t*Tbeta*(-1 + pow2(Sbeta)) + 8*Mt*pow2(Sbeta))) - 4*Mst2*pow2(Mt)*(-(
        Tbeta*pow2(MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta))*(-35*z3 + 26*z4 + 39*
        pow2(z2))) + 2*Mt*pow2(Sbeta)*(288*Mt*Tbeta*z3 + MuSUSY*s2t*(-202*z3 +
        10*z4 + 15*pow2(z2)))) + pow2(Sbeta)*(8*Mt*pow2(s2t)*(Mt*Tbeta*(-14*z3
        + 2*z4 + 3*pow2(z2)) + MuSUSY*s2t*(-52*z3 + 10*z4 + 15*pow2(z2)))*pow3(
        Mst2) + Tbeta*(173*z3 - 14*z4 - 21*pow2(z2))*pow4(s2t)*pow5(Mst2))) +
        27*pow3(Mst2)*(96*Mt*z3*(s2t*pow2(Mst2)*(-9*Mt*MuSUSY*s2t + 8*Tbeta*
        pow2(Mt) + 3*Tbeta*pow2(Mst2)*pow2(s2t))*pow2(Sbeta) - 2*MuSUSY*pow2(
        Mt)*(3*MuSUSY*s2t*Tbeta*(-1 + pow2(Sbeta)) + 4*Mt*pow2(Sbeta))) - 4*
        Mst2*pow2(Mt)*(-(Tbeta*pow2(MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta))*(23*z3
        + 2*z4 + 3*pow2(z2))) + 2*Mt*pow2(Sbeta)*(48*Mt*Tbeta*z3 + MuSUSY*s2t*(
        -58*z3 + 2*z4 + 3*pow2(z2)))) + pow2(Sbeta)*(4*Mt*pow2(s2t)*(Mt*Tbeta*(
        -58*z3 + 2*z4 + 3*pow2(z2)) + MuSUSY*s2t*(23*z3 + 2*z4 + 3*pow2(z2)))*
        pow3(Mst2) - Tbeta*(23*z3 + 2*z4 + 3*pow2(z2))*pow4(s2t)*pow5(Mst2))))
        + Dmsqst2*xDmsqst2*(-(pow4(Mst1)*(pow2(Sbeta)*(128*s2t*(-2592*Mst2*
        Tbeta*z3 + MuSUSY*(-1930*z3 + 106*z4 + 159*pow2(z2)))*pow3(Mt) +
        497664*Tbeta*z3*pow4(Mt) - 10773*Tbeta*z3*pow4(Mst2)*pow4(s2t) - 192*
        Mt*pow2(Mst2)*pow3(s2t)*(MuSUSY*(526*z3 + 38*z4 + 57*pow2(z2)) + 135*
        Mst2*Tbeta*z3*(-2 - 2*pow2(Sbeta) + pow4(Sbeta)))) + 32*pow2(Mt)*pow2(
        s2t)*(-2*Tbeta*pow2(MuSUSY)*(-1 + pow2(Sbeta))*(-1771*z3 + 442*z4 +
        663*pow2(z2)) + Mst2*pow2(Sbeta)*(15552*MuSUSY*z3 + Mst2*Tbeta*(28*z4 +
        42*pow2(z2) + z3*(-196 - 486*pow2(Sbeta) + 243*pow4(Sbeta))))))) + 6*
        Mst2*pow2(Mst1)*(3456*Mt*z3*(s2t*pow2(Mst2)*(-30*Mt*MuSUSY*s2t + 32*
        Tbeta*pow2(Mt) + 3*Tbeta*pow2(Mst2)*pow2(s2t))*pow2(Sbeta) - 4*MuSUSY*
        pow2(Mt)*(7*MuSUSY*s2t*Tbeta*(-1 + pow2(Sbeta)) + 8*Mt*pow2(Sbeta))) -
        32*Mst2*pow2(Mt)*(-(Tbeta*pow2(MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta))*(
        1783*z3 + 50*z4 + 75*pow2(z2))) + 2*Mt*pow2(Sbeta)*(1944*Mt*Tbeta*z3 +
        MuSUSY*s2t*(-494*z3 + 14*z4 + 21*pow2(z2)))) + 2*Mt*pow2(s2t)*pow2(
        Sbeta)*(2*Mt*Tbeta*(-679*z3 + 16*z4 + 24*pow2(z2)) + MuSUSY*s2t*(12007*
        z3 + 608*z4 + 912*pow2(z2)))*pow3(Mst2) - Tbeta*pow2(Sbeta)*(-2257*z3 +
        208*z4 + 312*pow2(z2))*pow4(s2t)*pow5(Mst2)) + 9*pow3(Mst2)*(2304*Mt*
        z3*(s2t*pow2(Mst2)*(-21*Mt*MuSUSY*s2t + 24*Tbeta*pow2(Mt) + 7*Tbeta*
        pow2(Mst2)*pow2(s2t))*pow2(Sbeta) - 2*MuSUSY*pow2(Mt)*(7*MuSUSY*s2t*
        Tbeta*(-1 + pow2(Sbeta)) + 12*Mt*pow2(Sbeta))) - 4*Mst2*pow2(Mt)*(-(
        Tbeta*pow2(MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta))*(5507*z3 + 64*z4 + 96*
        pow2(z2))) + 4*Mt*pow2(Sbeta)*(1827*Mt*Tbeta*z3 + MuSUSY*s2t*(-1091*z3
        + 32*z4 + 48*pow2(z2)))) + pow2(Sbeta)*(4*Mt*pow2(s2t)*(2*Mt*Tbeta*(-
        1091*z3 + 32*z4 + 48*pow2(z2)) + MuSUSY*s2t*(5507*z3 + 64*z4 + 96*pow2(
        z2)))*pow3(Mst2) - Tbeta*(5507*z3 + 64*z4 + 96*pow2(z2))*pow4(s2t)*
        pow5(Mst2))))) + 4*pow4(Msq)*(6*pow2(Mst1)*pow2(Mst2)*(-4*pow2(Mst2)*
        pow2(Mt)*(2*Tbeta*pow2(Mt)*pow2(Sbeta)*(69050*z3 + 418*z4 + 627*pow2(
        z2)) - Tbeta*pow2(MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta))*((-71438 - 216*
        lmMst1 + 216*lmMst2)*z3 + 2276*z4 + 2523*pow2(z2)) + Mt*MuSUSY*s2t*
        pow2(Sbeta)*(3718*z3 + 3470*z4 + 5205*pow2(z2))) + 2*Mt*s2t*pow2(Sbeta)
        *(-4*Tbeta*pow2(Mt)*(-48454*z3 + 754*z4 + 1131*pow2(z2)) + 4*Mst2*
        MuSUSY*pow2(s2t)*((-23848 - 54*lmMst1 + 54*lmMst2)*z3 + 958*z4 + 1275*
        pow2(z2)) + 2*Mst2*Mt*s2t*Tbeta*(4*(-55 + 54*lmMst1 - 54*lmMst2)*z3 +
        826*z4 + 1887*pow2(z2)) + 3*Mt*MuSUSY*s2t*(-22754*z3 - 3010*z4 + 4557*
        pow2(z2)))*pow3(Mst2) - 16*Mst2*MuSUSY*(-4*Mt*pow2(Sbeta)*(-8219*z3 +
        326*z4 + 165*pow2(z2)) - MuSUSY*s2t*Tbeta*(-1 + pow2(Sbeta))*(-14636*z3
        - 1135*z4 + 2388*pow2(z2)))*pow3(Mt) + Tbeta*(-144*pow2(MuSUSY)*(-1 +
        pow2(Sbeta))*(218*z3 + 50*z4 + 21*pow2(z2))*pow4(Mt) + pow2(Sbeta)*(
        Mst2*s2t*(23954*z3 - 1556*z4 - 2577*pow2(z2)) + 4*Mt*(-6518*z3 + 740*z4
        + 219*pow2(z2)))*pow3(s2t)*pow5(Mst2))) + pow4(Mst1)*(-4*pow2(Mst2)*
        pow2(Mt)*(Tbeta*pow2(MuSUSY)*pow2(s2t)*(1 - pow2(Sbeta))*(-8*(91963 +
        162*lmMst1 - 162*lmMst2)*z3 + 27446*z4 + 35823*pow2(z2)) + 4*Mt*pow2(
        Sbeta)*(16*Mt*Tbeta*(24241*z3 + 146*z4 + 219*pow2(z2)) + MuSUSY*s2t*(-
        27086*z3 + 12062*z4 + 18093*pow2(z2)))) + 4*Mt*s2t*pow2(Sbeta)*(-2*Mt*(
        3*MuSUSY*s2t*(11816*z3 + 4954*z4 - 6177*pow2(z2)) + 8*Mt*Tbeta*(-57758*
        z3 + 1370*z4 + 2055*pow2(z2))) + Mst2*s2t*(7*MuSUSY*s2t*(-43868*z3 +
        1970*z4 + 2955*pow2(z2)) + Mt*Tbeta*(-65326*z3 + 13714*z4 + 20571*pow2(
        z2))))*pow3(Mst2) - 16*Mst2*MuSUSY*(-4*Mt*pow2(Sbeta)*(-107072*z3 +
        3326*z4 + 3045*pow2(z2)) - MuSUSY*s2t*Tbeta*(-1 + pow2(Sbeta))*(-99632*
        z3 - 11764*z4 + 20505*pow2(z2)))*pow3(Mt) + Tbeta*(-48*pow2(MuSUSY)*(-1
        + pow2(Sbeta))*(1502*z3 + 2062*z4 + 1635*pow2(z2))*pow4(Mt) + 4*Mt*
        pow2(Sbeta)*(-44630*z3 + 878*z4 + 1317*pow2(z2))*pow3(s2t)*pow5(Mst2) +
        pow2(Sbeta)*((20900 - 648*lmMst1 + 648*lmMst2)*z3 - 2294*z4 - 5385*
        pow2(z2))*pow4(s2t)*pow6(Mst2))) + 18*pow4(Mst2)*(12*pow2(Mst2)*pow2(
        Mt)*(-(Mt*MuSUSY*s2t*pow2(Sbeta)*(6*(77 - 8*lmMst1 + 8*lmMst2)*z3 +
        202*z4 + 159*pow2(z2))) + Tbeta*(pow2(MuSUSY)*pow2(s2t)*(1 - pow2(
        Sbeta))*(2*(1319 + 6*lmMst1 - 6*lmMst2)*z3 - 40*z4 + 3*pow2(z2)) - 2*
        pow2(Mt)*pow2(Sbeta)*(8*(514 - 3*(lmMst1 + lmMst2) + 6*lmMt)*z3 + 14*z4
        + 21*pow2(z2)))) + 6*s2t*pow2(Mt)*pow2(Sbeta)*(-4*Mt*Tbeta*(-1922*z3 +
        206*z4 + 21*pow2(z2)) + 5*MuSUSY*s2t*(-2386*z3 - 102*z4 + 333*pow2(z2))
        )*pow3(Mst2) - 4*Mst2*MuSUSY*(-6*Mt*pow2(Sbeta)*(-1922*z3 + 206*z4 +
        21*pow2(z2)) - 5*MuSUSY*s2t*Tbeta*(-1 + pow2(Sbeta))*(-2386*z3 - 102*z4
        + 333*pow2(z2)))*pow3(Mt) - 6*Mt*pow2(s2t)*pow2(Sbeta)*(2*MuSUSY*s2t*(
        2*(1319 + 6*lmMst1 - 6*lmMst2)*z3 - 40*z4 + 3*pow2(z2)) - Mt*Tbeta*(6*(
        77 - 8*lmMst1 + 8*lmMst2)*z3 + 202*z4 + 159*pow2(z2)))*pow4(Mst2) +
        Tbeta*(-32*(230*z3 + 27*z4)*pow2(MuSUSY)*(-1 + pow2(Sbeta))*pow4(Mt) +
        10*Mt*pow2(Sbeta)*(2386*z3 + 102*z4 - 333*pow2(z2))*pow3(s2t)*pow5(
        Mst2) + 3*pow2(Sbeta)*(2*(1319 + 6*lmMst1 - 6*lmMst2)*z3 - 40*z4 + 3*
        pow2(z2))*pow4(s2t)*pow6(Mst2)))))))/(pow4(Msq)*pow6(Mst2)) - 2*z2*(6*
        xDmglst2*pow2(Dmglst2)*(Al4p*threeLoopFlag*((196560*(-4*Tbeta*pow2(Mt)*
        pow2(s2t)*(pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 2*pow2(Mst2)*pow2(Sbeta))
        + pow2(Sbeta)*(16*MuSUSY*s2t*pow3(Mt) - 4*Mt*MuSUSY*pow2(Mst2)*pow3(
        s2t) + 16*Tbeta*pow4(Mt) + Tbeta*pow4(Mst2)*pow4(s2t))))/pow2(Mst1) - (
        3*(-4*pow2(Mst2)*pow2(Mt)*(-7*(40001 + 52560*lmMst1 - 37080*lmMst2)*
        Tbeta*pow2(MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta)) + 4*Mt*(3*(255749 +
        2100*lmMst1 + 441420*lmMst2)*MuSUSY*s2t + (3572731 - 486360*lmMst1 +
        912240*lmMst2 - 587160*lmMt)*Mt*Tbeta)*pow2(Sbeta)) + 4*Mt*s2t*(Mst2*
        s2t*(7*(30641 + 52560*lmMst1 - 37080*lmMst2)*MuSUSY*s2t + 6*(233909 +
        2100*lmMst1 + 441420*lmMst2)*Mt*Tbeta) - 4*Mt*(21*(46987 + 15390*lmMst1
        + 4050*lmMst2)*MuSUSY*s2t + 4*(-497573 + 107730*lmMst1 - 233100*lmMst2
        + 125370*lmMt)*Mt*Tbeta))*pow2(Sbeta)*pow3(Mst2) - 32*Mst2*MuSUSY*(7*(
        46987 + 15390*lmMst1 + 4050*lmMst2)*MuSUSY*s2t*Tbeta*(-1 + pow2(Sbeta))
        + 4*(217444 - 53865*lmMst1 + 116550*lmMst2 - 62685*lmMt)*Mt*pow2(Sbeta)
        )*pow3(Mt) + Tbeta*(-5376*(17 + 1320*lmMst2)*pow2(MuSUSY)*(-1 + pow2(
        Sbeta))*pow4(Mt) + 7*(16*(46987 + 15390*lmMst1 + 4050*lmMst2)*Mt + (-
        21281 - 52560*lmMst1 + 37080*lmMst2)*Mst2*s2t)*pow2(Sbeta)*pow3(s2t)*
        pow5(Mst2))))/pow4(Mst2) + (4*pow2(Mst1)*(2*pow2(Mst2)*pow2(Mt)*(-7*(-
        261247 + 163890*lmMst1 - 140670*lmMst2)*Tbeta*pow2(MuSUSY)*pow2(s2t)*(-
        1 + pow2(Sbeta)) + 2*Mt*(36*(297316 - 20265*lmMst1 + 231945*lmMst2)*
        MuSUSY*s2t + (71735177 - 8595720*lmMst1 + 18272520*lmMst2 - 6940080*
        lmMt)*Mt*Tbeta)*pow2(Sbeta)) - s2t*(-21*(-1486429 + 190080*lmMst1 +
        43200*lmMst2)*MuSUSY*s2t + 32*(4999046 - 615195*lmMst1 + 1302210*lmMst2
        - 566055*lmMt)*Mt*Tbeta)*pow2(Mt)*pow2(Sbeta)*pow3(Mst2) + 2*Mst2*
        MuSUSY*(35*(-184517 + 74952*lmMst1 + 18360*lmMst2)*MuSUSY*s2t*Tbeta*(-1
        + pow2(Sbeta)) + 8*(10736701 - 1553580*lmMst1 + 3303720*lmMst2 -
        1508220*lmMt)*Mt*pow2(Sbeta))*pow3(Mt) + Mt*(7*(642497 - 170100*lmMst1
        + 170100*lmMst2)*MuSUSY*s2t + 90*(-186703 + 16632*lmMst1 - 97272*
        lmMst2)*Mt*Tbeta)*pow2(s2t)*pow2(Sbeta)*pow4(Mst2) + 84*(12761 + 5400*
        lmMst1 + 178920*lmMst2)*Tbeta*pow2(MuSUSY)*(-1 + pow2(Sbeta))*pow4(Mt)
        + 7*((2050273 - 5400*lmMst1 + 5400*lmMst2)*Mt + 5*(-36721 + 621*lmMst1
        - 2943*lmMst2)*Mst2*s2t)*Tbeta*pow2(Sbeta)*pow3(s2t)*pow5(Mst2)))/pow6(
        Mst2)) + (22680*twoLoopFlag*(-4*pow2(Mt)*(-8*(7*MuSUSY + 15*Mst2*Tbeta)
        *pow2(Mt)*pow2(Sbeta) + Mst2*MuSUSY*pow2(s2t)*(MuSUSY*Tbeta*(-1 + pow2(
        Sbeta)) - 18*Mst2*pow2(Sbeta)) + 12*Mt*s2t*Tbeta*(-3*pow2(MuSUSY)*(-1 +
        pow2(Sbeta)) + 4*pow2(Mst2)*pow2(Sbeta)))*pow4(Mst1) + s2t*(4*(12*Mt -
        Mst2*s2t)*Tbeta*pow2(Mt)*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + s2t*pow2(
        Mst2)*pow2(Sbeta)*(-4*Mst2*Mt*s2t*(MuSUSY + 6*Mst2*Tbeta) + 72*MuSUSY*
        pow2(Mt) + Tbeta*pow2(s2t)*pow3(Mst2)))*pow4(Mst2) - pow2(Mst1)*pow2(
        Mst2)*(-4*Mst2*MuSUSY*pow2(Mt)*pow2(s2t)*(-(MuSUSY*Tbeta*(-1 + pow2(
        Sbeta))) + 18*Mst2*pow2(Sbeta)) + 32*s2t*Tbeta*(-3*pow2(MuSUSY)*(-1 +
        pow2(Sbeta)) + pow2(Mst2)*pow2(Sbeta))*pow3(Mt) - 32*(MuSUSY + 3*Mst2*
        Tbeta)*pow2(Sbeta)*pow4(Mt) + Tbeta*pow2(Sbeta)*pow4(s2t)*pow5(Mst2))))
        /pow7(Mst2)) - (420*Al4p*threeLoopFlag*(108*xDR2DRMOD*(8*Mt*MuSUSY*
        pow2(Sbeta)*(pow2(Mst1)*pow2(Mst2)*pow2(s2t)*(-(s2t*pow3(Mst2)*((7*
        Dmglst2 - 13*Mst2)*pow2(Mst1) + 2*Dmglst2*(7*lmMst1 - 15*lmMst2)*pow2(
        Mst2) - 2*(4 + 13*lmMst1 - 9*lmMst2)*pow3(Mst2))) - 48*(Dmglst2 + 3*
        Dmglst2*lmMst2 + Mst2 + lmMst2*Mst2)*Mt*(pow2(Mst1)*pow2(Mst2) + pow4(
        Mst1) + pow4(Mst2))) + 64*pow2(Mst1)*pow3(Mt)*(-(Dmglst2*(7 + lmMst2)*
        pow2(Mst1)*pow2(Mst2)) - 2*Dmglst2*(11 + 5*lmMst2)*pow4(Mst1) +
        Dmglst2*(-1 + lmMst2)*pow4(Mst2) + (1 + lmMst2)*Mst2*(3*pow2(Mst1)*
        pow2(Mst2) + 6*pow4(Mst1) + pow4(Mst2))) + 4*(7*Dmglst2 - 13*Mst2)*s2t*
        (pow2(Mst1) - pow2(Mst2))*pow2(Mt)*pow5(Mst2)) + 8*Mt*pow2(Sbeta)*(-80*
        Tbeta*xDmglst2*pow2(Dmglst2)*pow3(Mt) + (7*Dmglst2 - 13*Mst2)*MuSUSY*
        pow3(Mst2)*pow3(s2t))*pow6(Mst2) + Mst2*xDmglst2*pow2(Dmglst2)*(32*Mt*
        MuSUSY*pow2(Sbeta)*(12*(2 - 3*lmMst2)*Mt*pow2(Mst1)*pow2(Mst2)*(pow2(
        Mst1) + pow2(Mst2))*pow2(s2t) + 20*s2t*(pow2(Mst1) - pow2(Mst2))*pow2(
        Mt)*pow3(Mst2) + 16*pow2(Mst1)*((7 + 2*lmMst2)*pow2(Mst1) - pow2(Mst2))
        *pow3(Mt) - pow3(Mst2)*pow3(s2t)*(2*(2 + 5*lmMst1 - 6*lmMst2)*pow2(
        Mst1)*pow2(Mst2) + 5*pow4(Mst1) - 5*pow4(Mst2))) + Tbeta*(-32*s2t*pow2(
        Mt)*pow2(MuSUSY)*((-8*(-2 + 3*lmMst2)*Mt + (1 - 10*lmMst1 + 12*lmMst2)*
        Mst2*s2t)*pow2(Mst1)*pow2(Mst2) - 2*(8*(-2 + 3*lmMst2)*Mt + (2 + 5*
        lmMst1 - 6*lmMst2)*Mst2*s2t)*pow4(Mst1) + 5*s2t*pow5(Mst2)) + pow2(
        Sbeta)*(128*pow2(Mst1)*(pow2(Mst2)*(-(Mst2*Mt) + 4*s2t*pow2(Mst2) + 2*(
        2 - 3*lmMst2)*s2t*pow2(MuSUSY)) - 4*pow2(Mst1)*(-3*(8 + 3*lmMst2)*Mst2*
        Mt + 2*(4 + lmMst2)*s2t*pow2(Mst2) + (-2 + 3*lmMst2)*s2t*pow2(MuSUSY)))
        *pow3(Mt) - 32*Mst2*pow2(Mt)*pow2(s2t)*((2*(2 + 5*lmMst1 - 6*lmMst2)*
        pow2(Mst1) - 5*pow2(Mst2))*(pow2(Mst1) + pow2(Mst2))*pow2(MuSUSY) + 20*
        pow2(Mst1)*pow4(Mst2) - 10*pow2(Mst2)*(pow4(Mst1) + pow4(Mst2))) + (8*(
        9 + 10*lmMst1 - 12*lmMst2)*pow2(Mst1)*pow2(Mst2) + (8 - 80*lmMst1 + 96*
        lmMst2)*pow4(Mst1) - 40*pow4(Mst2))*pow4(s2t)*pow5(Mst2) + 128*(-2 + 3*
        lmMst2)*Mt*pow2(Mst1)*pow3(s2t)*pow6(Mst2)))) + Tbeta*(s2t*pow2(Mt)*
        pow2(MuSUSY)*(256*(Dmglst2 + 3*Dmglst2*lmMst2 + Mst2 + lmMst2*Mst2)*Mt*
        pow2(Mst1)*(2*pow2(Mst1)*pow2(Mst2) + 3*pow4(Mst1) + pow4(Mst2)) + 8*
        Mst2*s2t*(2*(Dmglst2*(7*lmMst1 - 15*lmMst2) + (-4 - 13*lmMst1 + 9*
        lmMst2)*Mst2)*(pow2(Mst1) + pow2(Mst2))*pow4(Mst1) - Dmglst2*((7 - 14*
        lmMst1 + 30*lmMst2)*pow2(Mst1) + 7*pow2(Mst2))*pow4(Mst2) + (5 - 26*
        lmMst1 + 18*lmMst2)*pow2(Mst1)*pow5(Mst2) + 13*pow7(Mst2))) + pow2(
        Sbeta)*(2*(pow2(Mst1) - pow2(Mst2))*pow4(s2t)*(2*(4 + 13*lmMst1 - 9*
        lmMst2)*pow2(Mst1)*pow3(Mst2) + 13*Mst2*pow4(Mst1) + Dmglst2*(2*(-7*
        lmMst1 + 15*lmMst2)*pow2(Mst1)*pow2(Mst2) - 7*pow4(Mst1) + 7*pow4(Mst2)
        ) - 13*pow5(Mst2))*pow5(Mst2) + 32*Mst2*pow4(Mt)*(48*(1 + lmMst2)*Mst2*
        (2*pow2(Mst1) + pow2(Mst2))*pow4(Mst1) + (29 + 16*lmMst2)*pow2(Mst1)*
        pow5(Mst2) - Dmglst2*(96*(2 + lmMst2)*pow2(Mst2)*pow4(Mst1) + 39*pow2(
        Mst1)*pow4(Mst2) + 192*(3 + 2*lmMst2)*pow6(Mst1) + 7*pow6(Mst2)) + 13*
        pow7(Mst2)) - 256*s2t*pow2(Mst1)*pow3(Mt)*(-6*Dmglst2*(5 + 3*lmMst2)*
        pow2(Mst2)*pow4(Mst1) + 4*(-(Dmglst2*(3 + lmMst2)) + (1 + lmMst2)*Mst2)
        *pow2(Mst1)*pow4(Mst2) + (Dmglst2 + 3*Dmglst2*lmMst2 + Mst2 + lmMst2*
        Mst2)*pow2(MuSUSY)*(2*pow2(Mst1)*pow2(Mst2) + 3*pow4(Mst1) + pow4(Mst2)
        ) + 2*Dmglst2*(-1 + lmMst2)*pow6(Mst2) + 2*(1 + lmMst2)*(3*pow3(Mst2)*
        pow4(Mst1) + pow7(Mst2))) + 128*(Dmglst2 + 3*Dmglst2*lmMst2 + Mst2 +
        lmMst2*Mst2)*Mt*pow2(Mst1)*pow3(s2t)*pow8(Mst2) + 8*Mst2*pow2(Mt)*pow2(
        s2t)*(2*(7*Dmglst2 - 13*Mst2)*pow4(Mst1)*pow4(Mst2) + 4*(-7*Dmglst2 +
        13*Mst2)*pow2(Mst1)*pow6(Mst2) + pow2(MuSUSY)*(2*(4 + 13*lmMst1 - 9*
        lmMst2)*Mst2*(pow2(Mst1) + pow2(Mst2))*pow4(Mst1) + 7*Dmglst2*(pow2(
        Mst1) + pow2(Mst2))*pow4(Mst2) + Dmglst2*(-14*lmMst1 + 30*lmMst2)*pow2(
        Mst1)*(pow2(Mst1)*pow2(Mst2) + pow4(Mst1) + pow4(Mst2)) + (-5 + 26*
        lmMst1 - 18*lmMst2)*pow2(Mst1)*pow5(Mst2) - 13*pow7(Mst2)) + 2*(7*
        Dmglst2 - 13*Mst2)*pow8(Mst2))))) + (5*xDmsqst2*pow2(Dmsqst2)*pow2(
        Mst2)*(4*Mt*MuSUSY*pow2(Sbeta)*(-18*pow2(Mst1)*pow3(Mst2)*(2*Mst2*s2t*(
        25 - 36*shiftst1)*pow2(Mt) + 108*(7 + 3*lmMst1 - 3*lmMst2)*Mt*pow2(
        Mst2)*pow2(s2t) + 144*(2 - 3*lmMst1 + 6*lmMst2 - 3*lmMt)*pow3(Mt) + (2
        + lmMst2*(9 - 36*shiftst1) + 9*lmMst1*(-1 + 4*shiftst1))*pow3(Mst2)*
        pow3(s2t)) - 3*Mst2*(530*Mst2*s2t*pow2(Mt) + 1296*(4 + 3*lmMst1 - 3*
        lmMst2)*Mt*pow2(Mst2)*pow2(s2t) - 3456*(-1 + lmMst1 - 2*lmMst2 + lmMt)*
        pow3(Mt) + (1757 - 378*lmMst1 + 378*lmMst2 + 108*shiftst1)*pow3(Mst2)*
        pow3(s2t))*pow4(Mst1) + 2*s2t*(-972*(5 + 4*lmMst1 - 4*lmMst2)*Mst2*Mt*
        s2t + 4057*pow2(Mt) - 3324*pow2(Mst2)*pow2(s2t))*pow6(Mst1) + 162*s2t*(
        3 + 2*shiftst1)*(-4*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow6(Mst2) + 324*
        Dmglst2*(2*pow2(Mst1)*pow2(Mst2)*(20*Mst2*s2t*pow2(Mt) + 6*(2 + 3*
        lmMst1 - 3*lmMst2)*Mt*pow2(Mst2)*pow2(s2t) - 8*(-8 + 3*lmMst1 - 6*
        lmMst2 + 3*lmMt)*pow3(Mt) + (-5 - 3*lmMst1 + 3*lmMst2)*pow3(Mst2)*pow3(
        s2t)) - 2*(-36*Mst2*s2t*pow2(Mt) + 18*(1 - 2*lmMst1 + 2*lmMst2)*Mt*
        pow2(Mst2)*pow2(s2t) + 80*(-3 + lmMst1 - 2*lmMst2 + lmMt)*pow3(Mt) + (-
        3 + 5*lmMst1 - 5*lmMst2)*pow3(Mst2)*pow3(s2t))*pow4(Mst1) + 8*s2t*pow2(
        Mt)*pow5(Mst2) + (6*(-17 + 12*lmMst1 - 12*lmMst2)*Mt + 13*Mst2*s2t)*
        pow2(s2t)*pow6(Mst1) - 2*pow3(s2t)*pow7(Mst2))) + Tbeta*(5184*s2t*pow2(
        Mst1)*pow2(MuSUSY)*(2*Dmglst2*(1 - 9*lmMst1 + 9*lmMst2)*pow2(Mst1) + 2*
        (5 + 3*lmMst1 - 3*lmMst2)*Mst2*pow2(Mst1) - 2*Dmglst2*(2 + 3*lmMst1 -
        3*lmMst2)*pow2(Mst2) + (7 + 3*lmMst1 - 3*lmMst2)*pow3(Mst2))*pow3(Mt) +
        2*Mt*pow2(s2t)*(-81*Mst2*(-2*(11 + 10*lmMst1 - 10*lmMst2)*s2t*pow2(
        Mst2)*(-2 + pow2(Sbeta)) + (Dmglst2*(372 - 40*lmMst1 + 40*lmMst2) + 3*(
        -17 + 2*lmMst1 - 2*lmMst2)*Mst2)*Mt*pow2(Sbeta))*pow4(Sbeta)*pow6(Mst1)
        + 2*Mt*(3*Mst2*pow2(MuSUSY)*(-6*(25 + 9*lmMst1 - 9*lmMst2)*pow2(Mst1)*
        pow3(Mst2) + (1607 - 432*lmMst1 + 432*lmMst2)*Mst2*pow4(Mst1) + 216*
        Dmglst2*(3*(2 + lmMst1 - lmMst2)*pow2(Mst1)*pow2(Mst2) + (3 + 8*lmMst1
        - 8*lmMst2)*pow4(Mst1) + pow4(Mst2)) - 162*pow5(Mst2)) + ((20701 - 648*
        lmMst1 + 648*lmMst2)*pow2(MuSUSY) + 81*Mst2*(Dmglst2*(372 - 40*lmMst1 +
        40*lmMst2) + 3*(-17 + 2*lmMst1 - 2*lmMst2)*Mst2)*pow4(Sbeta))*pow6(
        Mst1))) + pow2(Sbeta)*(-5184*s2t*pow2(Mst1)*pow3(Mt)*(-2*Dmglst2*((2 +
        3*lmMst1 - 3*lmMst2)*pow2(Mst2)*pow2(MuSUSY) + pow2(Mst1)*((-51 + 14*
        lmMst1 - 28*lmMst2 + 14*lmMt)*pow2(Mst2) + (-1 + 9*lmMst1 - 9*lmMst2)*
        pow2(MuSUSY)) + (-87 + 22*lmMst1 - 44*lmMst2 + 22*lmMt)*pow4(Mst1) + (-
        19 + 6*lmMst1 - 12*lmMst2 + 6*lmMt)*pow4(Mst2)) + Mst2*((7 + 3*lmMst1 -
        3*lmMst2)*pow2(Mst2)*pow2(MuSUSY) + 2*pow2(Mst1)*(4*(-3 + lmMst1 - 2*
        lmMst2 + lmMt)*pow2(Mst2) + (5 + 3*lmMst1 - 3*lmMst2)*pow2(MuSUSY)) +
        2*(-5 + 2*lmMst1 - 4*lmMst2 + 2*lmMt)*pow4(Mst1) + (-7 + 6*lmMst1 - 12*
        lmMst2 + 6*lmMt)*pow4(Mst2))) + 648*Mt*pow2(Mst1)*pow2(Mst2)*pow3(s2t)*
        (4*(Dmglst2*(10 - 6*lmMst1 + 6*lmMst2) + (1 + 3*lmMst1 - 3*lmMst2)*
        Mst2)*pow2(Mst1)*pow2(Mst2) + (44*Dmglst2 + (-21 - 10*lmMst1 + 10*
        lmMst2)*Mst2)*pow4(Mst1) - 8*Dmglst2*(2 + 3*lmMst1 - 3*lmMst2)*pow4(
        Mst2) + 4*(7 + 3*lmMst1 - 3*lmMst2)*pow5(Mst2)) - 2592*pow4(Mt)*((4*
        Dmglst2*(-23 + 5*lmMst1 - 10*lmMst2 + 5*lmMt) + (25 - 7*lmMst1 + 14*
        lmMst2 - 7*lmMt)*Mst2)*pow2(Mst1)*pow3(Mst2) + Mst2*(4*Dmglst2*(-77 +
        18*lmMst1 - 36*lmMst2 + 18*lmMt) + 3*(23 - 6*lmMst1 + 12*lmMst2 - 6*
        lmMt)*Mst2)*pow4(Mst1) + (-4*Dmglst2 + 3*Mst2)*pow5(Mst2) - 4*(-7 + 3*
        lmMst1 - 6*lmMst2 + 3*lmMt)*pow6(Mst1)) + 4*pow2(Mt)*pow2(s2t)*(3*pow2(
        Mst2)*(115*pow2(Mst2) + (-1607 + 432*lmMst1 - 432*lmMst2)*pow2(MuSUSY))
        *pow4(Mst1) - 18*pow2(Mst1)*(29*pow2(Mst2) + (-25 - 9*lmMst1 + 9*
        lmMst2)*pow2(MuSUSY))*pow4(Mst2) + (1355*pow2(Mst2) + (-20701 + 648*
        lmMst1 - 648*lmMst2)*pow2(MuSUSY))*pow6(Mst1) - 648*Dmglst2*Mst2*(pow2(
        Mst1)*pow2(Mst2)*(8*pow2(Mst2) + 3*(2 + lmMst1 - lmMst2)*pow2(MuSUSY))
        + (8*pow2(Mst2) + (3 + 8*lmMst1 - 8*lmMst2)*pow2(MuSUSY))*pow4(Mst1) +
        (2*pow2(Mst2) + pow2(MuSUSY))*pow4(Mst2) + 8*pow6(Mst1)) + 486*(2*pow2(
        Mst2) + pow2(MuSUSY))*pow6(Mst2)) + 3*pow3(Mst2)*pow4(s2t)*((1745 -
        324*lmMst1 + 324*lmMst2)*pow3(Mst2)*pow4(Mst1) + 6*(29 - 9*lmMst1 + 9*
        lmMst2)*pow2(Mst1)*pow5(Mst2) + 27*(17 + 14*lmMst1 - 14*lmMst2)*Mst2*
        pow6(Mst1) - 108*Dmglst2*(4*(4 - lmMst1 + lmMst2)*pow2(Mst2)*pow4(Mst1)
        + (7 + 10*lmMst1 - 10*lmMst2)*pow6(Mst1) - 2*((4 + 3*lmMst1 - 3*lmMst2)
        *pow2(Mst1)*pow4(Mst2) + pow6(Mst2))) - 162*pow7(Mst2))) - 324*
        shiftst1*(4*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*(pow2(Mst1)*pow4(Mst2) - 2*
        (lmMst1 - lmMst2)*pow2(Mst1)*(pow2(Mst1)*pow2(Mst2) + pow4(Mst1) +
        pow4(Mst2)) + pow6(Mst2)) + pow2(Sbeta)*(16*(pow2(Mst1) + pow2(Mst2))*
        pow4(Mst2)*pow4(Mt) + (pow2(Mst1) - pow2(Mst2))*(2*(lmMst1 - lmMst2)*
        pow2(Mst1)*pow2(Mst2) + pow4(Mst1) - pow4(Mst2))*pow4(Mst2)*pow4(s2t) -
        4*pow2(Mt)*pow2(s2t)*(2*pow4(Mst1)*pow4(Mst2) - 4*pow2(Mst1)*pow6(Mst2)
        + pow2(MuSUSY)*(pow2(Mst1)*pow4(Mst2) - 2*(lmMst1 - lmMst2)*pow2(Mst1)*
        (pow2(Mst1)*pow2(Mst2) + pow4(Mst1) + pow4(Mst2)) + pow6(Mst2)) + 2*
        pow8(Mst2)))))))/pow4(Msq)))/(pow2(Mst1)*pow6(Mst2)) + (68040*
        twoLoopFlag*(pow5(Mst2)*(4*Mst2*MuSUSY*pow2(Mt)*pow2(s2t)*(-(MuSUSY*
        Tbeta*(-1 + pow2(Sbeta))) + 6*Mst2*pow2(Sbeta)) + 16*s2t*Tbeta*(pow2(
        MuSUSY)*(-1 + pow2(Sbeta)) - 4*pow2(Mst2)*pow2(Sbeta))*pow3(Mt) - 4*Mt*
        (MuSUSY + 2*Mst2*Tbeta)*pow2(Sbeta)*pow3(Mst2)*pow3(s2t) + 32*(2*MuSUSY
        + Mst2*Tbeta)*pow2(Sbeta)*pow4(Mt) + Tbeta*pow2(Sbeta)*pow4(s2t)*pow5(
        Mst2)) + Mst2*pow4(Mst1)*(16*s2t*Tbeta*(pow2(MuSUSY)*(-1 + pow2(Sbeta))
        - 4*pow2(Mst2)*pow2(Sbeta))*pow3(Mt) + 96*(2*MuSUSY + Mst2*Tbeta)*pow2(
        Sbeta)*pow4(Mt) + Tbeta*pow2(Sbeta)*pow4(s2t)*pow5(Mst2)) - 2*pow2(
        Mst1)*pow3(Mst2)*(pow2(Sbeta)*(-2*Mt*(MuSUSY + 2*Mst2*Tbeta)*pow3(Mst2)
        *pow3(s2t) - 32*(2*MuSUSY + Mst2*Tbeta)*pow4(Mt)) + Tbeta*(8*s2t*(pow2(
        MuSUSY) + (4*pow2(Mst2) - pow2(MuSUSY))*pow2(Sbeta))*pow3(Mt) + pow2(
        Sbeta)*pow4(s2t)*pow5(Mst2))) + 4*Dmglst2*(4*pow2(Mt)*(-12*(MuSUSY + 2*
        Mst2*Tbeta)*pow2(Mt)*pow2(Sbeta) + Mst2*MuSUSY*pow2(s2t)*(-(MuSUSY*
        Tbeta*(-1 + pow2(Sbeta))) + 6*Mst2*pow2(Sbeta)) + Mt*s2t*Tbeta*(13*
        pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 12*pow2(Mst2)*pow2(Sbeta)))*pow4(
        Mst1) + pow2(Mst1)*pow2(Mst2)*(8*s2t*(3*MuSUSY*s2t + 2*Mt*Tbeta)*pow2(
        Mst2)*pow2(Mt)*pow2(Sbeta) + Tbeta*(-4*Mst2*pow2(Mt)*pow2(MuSUSY)*pow2(
        s2t)*(-1 + pow2(Sbeta)) + 36*s2t*pow2(MuSUSY)*(-1 + pow2(Sbeta))*pow3(
        Mt) + (2*Mt - Mst2*s2t)*pow2(Sbeta)*pow3(s2t)*pow4(Mst2) - 32*Mst2*
        pow2(Sbeta)*pow4(Mt))) + pow4(Mst2)*(2*Mst2*MuSUSY*pow2(Mt)*pow2(s2t)*(
        -2*MuSUSY*Tbeta*(-1 + pow2(Sbeta)) + 15*Mst2*pow2(Sbeta)) + 4*s2t*
        Tbeta*(5*pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 4*pow2(Mst2)*pow2(Sbeta))*
        pow3(Mt) - 2*Mt*(2*MuSUSY + 5*Mst2*Tbeta)*pow2(Sbeta)*pow3(Mst2)*pow3(
        s2t) + 16*MuSUSY*pow2(Sbeta)*pow4(Mt) + Tbeta*pow2(Sbeta)*pow4(s2t)*
        pow5(Mst2)))) + (Al4p*threeLoopFlag*(28*Mt*MuSUSY*pow2(Sbeta)*(3*Mt*
        pow2(Mst1)*pow2(Mst2)*pow2(s2t)*(3*Mst2*pow2(Msq)*(15*(11075 + 17928*
        lmMst1 - 17352*lmMst2)*pow2(Mst1)*pow2(Mst2) + (61286 + 434160*lmMst1 -
        425520*lmMst2)*pow4(Mst1) + 135*(2269 + 664*lmMst1 - 56*lmMst2)*pow4(
        Mst2)) + 194400*Dmglst2*Dmsqst2*(6*(1 - 2*lmMst1 + 2*lmMst2)*pow2(Mst1)
        *pow2(Mst2) + (17 - 12*lmMst1 + 12*lmMst2)*pow4(Mst1) - 2*(2 + 3*lmMst1
        - 3*lmMst2)*pow4(Mst2)) + 194400*Dmsqst2*Mst2*(2*(3 + 2*lmMst1 - 2*
        lmMst2)*pow2(Mst1)*pow2(Mst2) + (5 + 4*lmMst1 - 4*lmMst2)*pow4(Mst1) +
        2*(2 + lmMst1 - lmMst2)*pow4(Mst2)) + Dmglst2*pow2(Msq)*(9*(935323 -
        279000*lmMst1 + 385560*lmMst2)*pow2(Mst1)*pow2(Mst2) + 10*(1700119 -
        484056*lmMst1 + 579960*lmMst2)*pow4(Mst1) + 27*(173947 - 25080*lmMst1 +
        68760*lmMst2)*pow4(Mst2))) + 4*pow2(Mst1)*pow3(Mt)*(388800*Dmsqst2*
        pow2(Mst2)*(10*Dmglst2*(-3 + lmMst1 - 2*lmMst2 + lmMt)*pow2(Mst1) - 2*(
        -1 + lmMst1 - 2*lmMst2 + lmMt)*Mst2*pow2(Mst1) + Dmglst2*(-8 + 3*lmMst1
        - 6*lmMst2 + 3*lmMt)*pow2(Mst2) - (lmMst1 - 2*lmMst2 + lmMt)*pow3(Mst2)
        ) + pow2(Msq)*(Dmglst2*(216*(-67843 + 10050*lmMst1 - 17490*lmMst2 +
        7560*lmMt)*pow2(Mst1)*pow2(Mst2) + 4*(-13463149 + 1492020*lmMst1 -
        2713500*lmMst2 + 625320*lmMt)*pow4(Mst1) + 405*(-3245 + 1672*lmMst1 -
        1448*lmMst2 + 1888*lmMt)*pow4(Mst2)) - 3*(48*(-41936 + 7245*lmMst1 -
        19665*lmMst2 + 720*lmMt)*pow2(Mst1)*pow3(Mst2) + 4*(-1077991 + 184140*
        lmMst1 - 400140*lmMst2 + 8640*lmMt)*Mst2*pow4(Mst1) + 9*(-66773 + 9960*
        lmMst1 - 45480*lmMst2 + 3840*lmMt)*pow5(Mst2)))) + 9720*s2t*(shiftst3*
        pow2(Msq)*pow2(Mst2)*(4*pow2(Mt) + ((1 + lmMst1 - lmMst2)*pow2(Mst1) -
        pow2(Mst2))*pow2(s2t)) + 10*shiftst1*(Dmsqst2 + pow2(Msq))*(-4*pow2(
        Mst1)*pow2(Mt) + 4*pow2(Mst2)*pow2(Mt) + 2*(lmMst1 - lmMst2)*pow2(Mst1)
        *pow2(Mst2)*pow2(s2t) + pow2(s2t)*(pow4(Mst1) - pow4(Mst2))))*pow6(
        Mst2) - pow3(Mst2)*pow3(s2t)*(4*Dmglst2*pow2(Msq)*(27*(32663 - 7620*
        lmMst1 + 7620*lmMst2)*pow2(Mst2)*pow4(Mst1) - 135*(648*lmMst1 - 7*(257
        + 240*lmMst2))*pow2(Mst1)*pow4(Mst2) + 4*(788723 - 80595*lmMst1 +
        80595*lmMst2)*pow6(Mst1) - 22680*pow6(Mst2)) - 3*Mst2*pow2(Msq)*(12*(
        224837 - 4680*lmMst1 + 8370*lmMst2)*pow2(Mst2)*pow4(Mst1) + 90*(5825 +
        564*lmMst1 - 1104*lmMst2)*pow2(Mst1)*pow4(Mst2) + (3447379 - 76140*
        lmMst1 + 76140*lmMst2)*pow6(Mst1) - 20520*pow6(Mst2)) + Dmsqst2*(-300*
        Mst2*(18*(469 - 36*lmMst1 + 36*lmMst2)*pow2(Mst2)*pow4(Mst1) + 1755*
        pow2(Mst1)*pow4(Mst2) + 10828*pow6(Mst1) - 324*pow6(Mst2)) + 97200*
        Dmglst2*(2*(3 - 5*lmMst1 + 5*lmMst2)*pow2(Mst2)*pow4(Mst1) + 13*pow6(
        Mst1) - 2*((5 + 3*lmMst1 - 3*lmMst2)*pow2(Mst1)*pow4(Mst2) + pow6(Mst2)
        )))) - 6*Mst2*s2t*pow2(Mt)*(129600*Dmglst2*Dmsqst2*pow2(Mst2)*(5*pow2(
        Mst1)*pow2(Mst2) + 9*pow4(Mst1) + pow4(Mst2)) + 100*Dmsqst2*Mst2*(3249*
        pow2(Mst2)*pow4(Mst1) + 1431*pow2(Mst1)*pow4(Mst2) + 4057*pow6(Mst1) -
        648*pow6(Mst2)) + pow2(Msq)*(3*(623599 + 65160*lmMst1 - 134280*lmMst2)*
        pow3(Mst2)*pow4(Mst1) + 9*(1367 + 13560*lmMst1 - 36600*lmMst2)*pow2(
        Mst1)*pow5(Mst2) + (5161826 + 406080*lmMst1 - 613440*lmMst2)*Mst2*pow6(
        Mst1) - 24*Dmglst2*(6*(19579 + 1770*lmMst1 + 12630*lmMst2)*pow2(Mst2)*
        pow4(Mst1) + 6*(4429 - 270*lmMst1 + 8910*lmMst2)*pow2(Mst1)*pow4(Mst2)
        + (510139 + 44280*lmMst1 + 76680*lmMst2)*pow6(Mst1) - 2520*pow6(Mst2))
        - 41040*pow7(Mst2)))) + Tbeta*(-56*pow2(Mst1)*pow2(MuSUSY)*pow3(Mt)*(
        388800*Dmsqst2*s2t*pow2(Mst2)*((Dmglst2*(1 - 9*lmMst1 + 9*lmMst2) + (5
        + 3*lmMst1 - 3*lmMst2)*Mst2)*pow2(Mst1) - Dmglst2*(2 + 3*lmMst1 - 3*
        lmMst2)*pow2(Mst2) + (2 + lmMst1 - lmMst2)*pow3(Mst2)) - pow2(Msq)*(36*
        Mst2*(Dmglst2*(30*(-3643 + 120*lmMst1 - 3192*lmMst2)*Mt + (-364291 +
        88560*lmMst1 - 147960*lmMst2)*Mst2*s2t) - 5*Mst2*(3*(-353 + 72*lmMst1 +
        696*lmMst2)*Mt + 2*(3937 + 2988*lmMst1 - 2232*lmMst2)*Mst2*s2t))*pow2(
        Mst1) - 27*(5*Mst2*(64*(53 + 24*lmMst2)*Mt + 3*(2269 + 664*lmMst1 - 56*
        lmMst2)*Mst2*s2t) + Dmglst2*(7680*(5 + 6*lmMst2)*Mt + (173947 - 25080*
        lmMst1 + 68760*lmMst2)*Mst2*s2t))*pow3(Mst2) + (90*(38401 + 1080*lmMst1
        - 7992*lmMst2)*Mt + 2*Dmglst2*(-15057833 + 4014360*lmMst1 - 5563080*
        lmMst2)*s2t - 6*(266863 + 396360*lmMst1 - 346680*lmMst2)*Mst2*s2t)*
        pow4(Mst1))) - 28*Mst2*pow2(Mt)*pow2(s2t)*(48600*Dmglst2*Dmsqst2*pow2(
        Mst2)*(4*pow2(MuSUSY)*(3*(2 + lmMst1 - lmMst2)*pow2(Mst1)*pow2(Mst2) +
        (3 + 8*lmMst1 - 8*lmMst2)*pow4(Mst1) + pow4(Mst2)) + (-93 + 10*lmMst1 -
        10*lmMst2)*(-2 + pow2(Sbeta))*pow4(Sbeta)*pow6(Mst1)) + 3*Mst2*pow2(
        MuSUSY)*(pow2(Msq)*(6*(533629 - 900*lmMst1 + 180*lmMst2)*pow2(Mst2)*
        pow4(Mst1) + 90*(5597 + 564*lmMst1 - 1104*lmMst2)*pow2(Mst1)*pow4(Mst2)
        + (6649153 - 81540*lmMst1 + 77220*lmMst2)*pow6(Mst1) - 20520*pow6(Mst2)
        ) + 100*Dmsqst2*(9*(1097 - 72*lmMst1 + 72*lmMst2)*pow2(Mst2)*pow4(Mst1)
        + 1431*pow2(Mst1)*pow4(Mst2) + (20701 - 648*lmMst1 + 648*lmMst2)*pow6(
        Mst1) - 324*pow6(Mst2))) - 4*Dmglst2*pow2(Msq)*pow2(MuSUSY)*(162*(6803
        - 1810*lmMst1 + 2670*lmMst2)*pow2(Mst2)*pow4(Mst1) - 135*(648*lmMst1 -
        7*(233 + 240*lmMst2))*pow2(Mst1)*pow4(Mst2) + (4256978 - 615600*lmMst1
        + 754920*lmMst2)*pow6(Mst1) - 22680*pow6(Mst2))) + pow2(Sbeta)*(28*Mt*
        pow2(Mst1)*pow3(s2t)*pow4(Mst2)*(3*Mst2*pow2(Msq)*(30*(4673 - 5976*
        lmMst1 + 8424*lmMst2)*pow2(Mst1)*pow2(Mst2) + 17*(6167 - 9720*lmMst1 +
        9720*lmMst2)*pow4(Mst1) - 135*(2269 + 664*lmMst1 - 56*lmMst2)*pow4(
        Mst2)) - Dmglst2*pow2(Msq)*(18*(206741 - 101880*lmMst1 + 89640*lmMst2)*
        pow2(Mst1)*pow2(Mst2) + (8583283 - 2329560*lmMst1 + 2329560*lmMst2)*
        pow4(Mst1) + 27*(173947 - 25080*lmMst1 + 68760*lmMst2)*pow4(Mst2)) +
        194400*Dmsqst2*(2*Dmglst2*(-5 + 3*lmMst1 - 3*lmMst2)*pow2(Mst1)*pow2(
        Mst2) - 2*(1 + lmMst1 - lmMst2)*pow2(Mst1)*pow3(Mst2) - 11*Dmglst2*
        pow4(Mst1) + Mst2*pow4(Mst1) + 2*Dmglst2*(2 + 3*lmMst1 - 3*lmMst2)*
        pow4(Mst2) - 2*(2 + lmMst1 - lmMst2)*pow5(Mst2))) - 7*pow4(s2t)*pow5(
        Mst2)*(-18*(24*Dmglst2*(900*Dmsqst2*(4 - lmMst1 + lmMst2) + (5917 -
        1095*lmMst1 - 195*lmMst2)*pow2(Msq)) - Mst2*(150*Dmsqst2*(743 - 72*
        lmMst1 + 72*lmMst2) + (362299 - 17820*lmMst1 + 33300*lmMst2)*pow2(Msq))
        )*pow2(Mst2)*pow4(Mst1) + 270*(720*Dmglst2*Dmsqst2*(4 + 3*lmMst1 - 3*
        lmMst2) + 2310*Dmsqst2*Mst2 + 2*Dmglst2*(648*lmMst1 - 7*(281 + 240*
        lmMst2))*pow2(Msq) + (6053 + 564*lmMst1 - 1104*lmMst2)*Mst2*pow2(Msq))*
        pow2(Mst1)*pow4(Mst2) - (-15*Mst2*(40*Dmsqst2*(1193 + 324*lmMst1 - 324*
        lmMst2) + (149867 - 3996*lmMst1 - 4860*lmMst2)*pow2(Msq)) + Dmglst2*(
        97200*Dmsqst2*(7 + 10*lmMst1 - 10*lmMst2) + 4*(2272991 - 116640*lmMst1
        + 116640*lmMst2)*pow2(Msq)))*pow6(Mst1) - 3240*(-60*Dmglst2*Dmsqst2 +
        30*Dmsqst2*Mst2 - 28*Dmglst2*pow2(Msq) + 19*Mst2*pow2(Msq))*pow6(Mst2))
        - 28*Mst2*pow2(Mt)*pow2(s2t)*(-3*Mst2*pow2(Msq)*(6*pow2(Mst2)*((309749
        + 12240*lmMst1 - 12240*lmMst2)*pow2(Mst2) + (533629 - 900*lmMst1 + 180*
        lmMst2)*pow2(MuSUSY))*pow4(Mst1) + 9*pow2(Mst1)*((5927 + 13560*lmMst1 -
        36600*lmMst2)*pow2(Mst2) + 10*(5597 + 564*lmMst1 - 1104*lmMst2)*pow2(
        MuSUSY))*pow4(Mst2) + ((3291029 + 210600*lmMst1 - 210600*lmMst2)*pow2(
        Mst2) + (6649153 - 81540*lmMst1 + 77220*lmMst2)*pow2(MuSUSY))*pow6(
        Mst1) - 20520*(2*pow2(Mst2) + pow2(MuSUSY))*pow6(Mst2)) - 300*Dmsqst2*
        Mst2*(9*pow2(Mst2)*(202*pow2(Mst2) + (1097 - 72*lmMst1 + 72*lmMst2)*
        pow2(MuSUSY))*pow4(Mst1) + (808*pow2(Mst2) + (20701 - 648*lmMst1 + 648*
        lmMst2)*pow2(MuSUSY))*pow6(Mst1) - 324*(2*pow2(Mst2) + pow2(MuSUSY))*
        pow6(Mst2) + 27*pow2(Mst1)*(53*pow2(MuSUSY)*pow4(Mst2) + 77*pow6(Mst2))
        ) + 4*Dmglst2*(-48600*Dmsqst2*pow2(Mst2)*(pow2(Mst1)*pow2(Mst2)*(8*
        pow2(Mst2) + 3*(2 + lmMst1 - lmMst2)*pow2(MuSUSY)) + (8*pow2(Mst2) + (3
        + 8*lmMst1 - 8*lmMst2)*pow2(MuSUSY))*pow4(Mst1) + (2*pow2(Mst2) + pow2(
        MuSUSY))*pow4(Mst2) + 8*pow6(Mst1)) + pow2(Msq)*(162*pow2(Mst2)*(20*(
        505 + 68*lmMst1 + 124*lmMst2)*pow2(Mst2) + (6803 - 1810*lmMst1 + 2670*
        lmMst2)*pow2(MuSUSY))*pow4(Mst1) - 27*pow2(Mst1)*(4*(-4849 + 270*lmMst1
        - 8910*lmMst2)*pow2(Mst2) + 5*(-1631 + 648*lmMst1 - 1680*lmMst2)*pow2(
        MuSUSY))*pow4(Mst2) + 2*(45*(78533 + 6732*lmMst1 + 180*lmMst2)*pow2(
        Mst2) + (2128489 - 307800*lmMst1 + 377460*lmMst2)*pow2(MuSUSY))*pow6(
        Mst1) - 22680*(2*pow2(Mst2) + pow2(MuSUSY))*pow6(Mst2)))) + 56*s2t*
        pow2(Mst1)*pow3(Mt)*(388800*Dmsqst2*pow2(Mst2)*(Dmglst2*((-2 - 3*lmMst1
        + 3*lmMst2)*pow2(Mst2)*pow2(MuSUSY) + pow2(Mst1)*((51 - 14*lmMst1 + 28*
        lmMst2 - 14*lmMt)*pow2(Mst2) + (1 - 9*lmMst1 + 9*lmMst2)*pow2(MuSUSY))
        + (87 - 22*lmMst1 + 44*lmMst2 - 22*lmMt)*pow4(Mst1) + (19 - 6*lmMst1 +
        12*lmMst2 - 6*lmMt)*pow4(Mst2)) + Mst2*((2 + lmMst1 - lmMst2)*pow2(
        Mst2)*pow2(MuSUSY) + pow2(Mst1)*((-5 + 2*lmMst1 - 4*lmMst2 + 2*lmMt)*
        pow2(Mst2) + (5 + 3*lmMst1 - 3*lmMst2)*pow2(MuSUSY)) + (-5 + 2*lmMst1 -
        4*lmMst2 + 2*lmMt)*pow4(Mst1) + (-1 + 2*lmMst1 - 4*lmMst2 + 2*lmMt)*
        pow4(Mst2))) + pow2(Msq)*(Dmglst2*(2*(4*(9908167 - 949320*lmMst1 +
        1769040*lmMst2 - 217080*lmMt)*pow2(Mst2) + (15057833 - 4014360*lmMst1 +
        5563080*lmMst2)*pow2(MuSUSY))*pow4(Mst1) + 18*pow2(Mst1)*(2*(364291 -
        88560*lmMst1 + 147960*lmMst2)*pow2(Mst2)*pow2(MuSUSY) + 3*(510149 -
        55320*lmMst1 + 118200*lmMst2 - 32160*lmMt)*pow4(Mst2)) + 27*((173947 -
        25080*lmMst1 + 68760*lmMst2)*pow2(MuSUSY)*pow4(Mst2) + 30*(-1672*lmMst1
        + 1448*lmMst2 + 59*(71 - 32*lmMt))*pow6(Mst2))) - 3*Mst2*(2*((2299036 -
        388800*lmMst1 + 656640*lmMst2)*pow2(Mst2) + (-266863 - 396360*lmMst1 +
        346680*lmMst2)*pow2(MuSUSY))*pow4(Mst1) + 6*pow2(Mst1)*(-20*(3937 +
        2988*lmMst1 - 2232*lmMst2)*pow2(Mst2)*pow2(MuSUSY) + (470657 - 86040*
        lmMst1 + 178200*lmMst2)*pow4(Mst2)) + 9*(-15*(2269 + 664*lmMst1 - 56*
        lmMst2)*pow2(MuSUSY)*pow4(Mst2) + 2*(68693 - 9960*lmMst1 + 45480*lmMst2
        - 3840*lmMt)*pow6(Mst2))))) + 16*pow4(Mt)*(680400*Dmsqst2*pow2(Mst2)*(-
        2*(-3 + lmMst1 - 2*lmMst2 + lmMt)*pow2(Mst1)*pow2(Mst2)*(2*pow2(Mst1) +
        pow2(Mst2)) + 2*(7 - 3*lmMst1 + 6*lmMst2 - 3*lmMt)*pow6(Mst1) + pow6(
        Mst2)) - 4*Dmglst2*(340200*Dmsqst2*pow3(Mst2)*((23 - 5*lmMst1 + 10*
        lmMst2 - 5*lmMt)*pow2(Mst1)*pow2(Mst2) + (77 - 18*lmMst1 + 36*lmMst2 -
        18*lmMt)*pow4(Mst1) + pow4(Mst2)) + Mst2*pow2(Msq)*(27*pow2(Mst1)*pow2(
        Mst2)*((209341 - 26880*lmMst1 + 82320*lmMst2 - 1680*lmMt)*pow2(Mst2) -
        6720*(5 + 6*lmMst2)*pow2(MuSUSY)) + 9*((3184397 - 476280*lmMst1 +
        1156680*lmMst2 - 161280*lmMt)*pow2(Mst2) + 105*(-3643 + 120*lmMst1 -
        3192*lmMst2)*pow2(MuSUSY))*pow4(Mst1) + 4*(19152737 - 3674160*lmMst1 +
        6872040*lmMst2 - 850500*lmMt)*pow6(Mst1) + 158760*pow6(Mst2))) + 21*
        pow2(Msq)*(3*pow4(Mst1)*(30*(-353 + 72*lmMst1 + 696*lmMst2)*pow2(Mst2)*
        pow2(MuSUSY) + (310039 - 91080*lmMst1 + 145080*lmMst2 - 43920*lmMt)*
        pow4(Mst2)) + (2*(626869 - 300240*lmMst1 + 395280*lmMst2 - 78840*lmMt)*
        pow2(Mst2) - 15*(38401 + 1080*lmMst1 - 7992*lmMst2)*pow2(MuSUSY))*pow6(
        Mst1) + 9*pow2(Mst1)*(160*(53 + 24*lmMst2)*pow2(MuSUSY)*pow4(Mst2) + (
        48983 - 12480*lmMst1 + 26820*lmMst2 - 12180*lmMt)*pow6(Mst2)) + 20520*
        pow8(Mst2)))) + 68040*pow2(Mst2)*(shiftst3*pow2(Msq)*(-4*pow2(Mst1)*
        pow2(s2t)*pow4(Mst2)*(2*(-1 + lmMst1 - lmMst2)*pow2(Mt)*(-(pow2(MuSUSY)
        *(-1 + pow2(Sbeta))) + pow2(Mst2)*pow2(Sbeta)) + pow2(s2t)*pow2(Sbeta)*
        pow4(Mst2)) + 8*(-1 + 3*lmMst1 - 3*lmMst2)*pow2(Mt)*pow2(MuSUSY)*pow2(
        s2t)*(-1 + pow2(Sbeta))*pow6(Mst1) + (-4*pow2(Mt)*pow2(s2t)*(pow2(
        MuSUSY)*(-1 + pow2(Sbeta)) + 2*pow2(Mst2)*pow2(Sbeta)) + pow2(Sbeta)*(
        16*pow4(Mt) + pow4(Mst2)*pow4(s2t)))*pow6(Mst2) + pow4(Mst1)*(8*(-1 +
        2*lmMst1 - 2*lmMst2)*pow2(Mst2)*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*(-1 +
        pow2(Sbeta)) + (3 + 2*lmMst1 - 2*lmMst2)*pow2(Sbeta)*pow4(s2t)*pow6(
        Mst2))) + 10*shiftst1*(Dmsqst2 + pow2(Msq))*(4*pow2(Mt)*pow2(MuSUSY)*
        pow2(s2t)*(pow2(Mst1)*pow4(Mst2) - 2*(lmMst1 - lmMst2)*pow2(Mst1)*(
        pow2(Mst1)*pow2(Mst2) + pow4(Mst1) + pow4(Mst2)) + pow6(Mst2)) + pow2(
        Sbeta)*(16*(pow2(Mst1) + pow2(Mst2))*pow4(Mst2)*pow4(Mt) + (pow2(Mst1)
        - pow2(Mst2))*(2*(lmMst1 - lmMst2)*pow2(Mst1)*pow2(Mst2) + pow4(Mst1) -
        pow4(Mst2))*pow4(Mst2)*pow4(s2t) - 4*pow2(Mt)*pow2(s2t)*(2*pow4(Mst1)*
        pow4(Mst2) - 4*pow2(Mst1)*pow6(Mst2) + pow2(MuSUSY)*(pow2(Mst1)*pow4(
        Mst2) - 2*(lmMst1 - lmMst2)*pow2(Mst1)*(pow2(Mst1)*pow2(Mst2) + pow4(
        Mst1) + pow4(Mst2)) + pow6(Mst2)) + 2*pow8(Mst2))))))))/(pow2(Msq)*
        pow2(Mst1)))/pow6(Mst2) + (544320*twoLoopFlag*xMst*pow2(Mt)*(xDmglst2*
        pow2(Dmglst2)*(16*(11*MuSUSY + 21*Mst2*Tbeta)*pow2(Mt)*pow2(Sbeta) +
        Mst2*MuSUSY*pow2(s2t)*(-(MuSUSY*Tbeta*(-1 + pow2(Sbeta))) + 18*Mst2*
        pow2(Sbeta)) + 24*Mt*s2t*Tbeta*(2*pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 5*
        pow2(Mst2)*pow2(Sbeta))) + 2*Mt*pow2(Mst2)*(8*Mt*(2*MuSUSY + Mst2*
        Tbeta)*pow2(Sbeta) + s2t*Tbeta*(pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 4*
        pow2(Mst2)*pow2(Sbeta))) + 2*Dmglst2*Mst2*(-16*(2*MuSUSY + 3*Mst2*
        Tbeta)*pow2(Mt)*pow2(Sbeta) + Mst2*MuSUSY*pow2(s2t)*(-(MuSUSY*Tbeta*(-1
        + pow2(Sbeta))) + 6*Mst2*pow2(Sbeta)) + Mt*s2t*Tbeta*(17*pow2(MuSUSY)*(
        -1 + pow2(Sbeta)) + 20*pow2(Mst2)*pow2(Sbeta))))*pow6(Mst1))/pow9(Mst2)
        )))/(816480.*Tbeta*pow2(Sbeta)) - (xMst*pow6(Mst1)*((-27*(lmMst1 -
        lmMst2)*oneLoopFlag*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta))*
        pow3(Mst2))/pow2(Sbeta) - Al4p*twoLoopFlag*(8*Mst2*s2t*pow2(Mt)*(Mst2*(
        Mt*(-193 - 474*lmMst2 + 6*lmMst1*(67 + 42*lmMst2) - 252*pow2(lmMst2)) +
        Dmglst2*s2t*(49 - 84*lmMst2 + 12*lmMst1*(7 + 3*lmMst2) - 36*pow2(
        lmMst2)))*pow2(MuSUSY) + Dmglst2*Mt*(-785 + 438*lmMst2 + 6*lmMst1*(-85
        + 42*lmMst2) - 252*pow2(lmMst2))*pow2(MuSUSY) + s2t*(1 - 111*lmMst2 +
        3*lmMst1*(37 + 72*lmMst2) - 81*pow2(lmMst1) - 135*pow2(lmMst2))*pow2(
        Mst2)*pow2(MuSUSY) + 4*(9*Dmglst2*(-lmMst1 + lmMst2)*s2t + Mt*(49 - 75*
        lmMst2 + 3*lmMst1*(25 + 6*lmMst2 - 6*lmMt) + 18*lmMst2*lmMt - 18*pow2(
        lmMst2)))*pow3(Mst2)) - 32*(Dmglst2*s2t*(95 - 447*lmMst2 + lmMst1*(411
        + 90*lmMst2 - 90*lmMt) + 36*lmMt + 90*lmMst2*lmMt - 90*pow2(lmMst2)) +
        Mt*(65 - 159*lmMst2 + 6*lmMst1*(25 + 6*lmMst2 - 6*lmMt) + 9*lmMt + 36*
        lmMst2*lmMt - 36*pow2(lmMst2)))*pow3(Mst2)*pow3(Mt) + 2*Mst2*Mt*MuSUSY*
        ((4*Mt*MuSUSY*s2t*(Dmglst2*Mst2*s2t*(-49 + 84*lmMst2 - 12*lmMst1*(7 +
        3*lmMst2) + 36*pow2(lmMst2)) + Dmglst2*Mt*(785 - 438*lmMst2 - 6*lmMst1*
        (-85 + 42*lmMst2) + 252*pow2(lmMst2)) + Mst2*Mt*(193 + 474*lmMst2 - 6*
        lmMst1*(67 + 42*lmMst2) + 252*pow2(lmMst2)) + s2t*(-1 + 111*lmMst2 - 3*
        lmMst1*(37 + 72*lmMst2) + 81*pow2(lmMst1) + 135*pow2(lmMst2))*pow2(
        Mst2)))/pow2(Sbeta) + (8*Dmglst2*(36*(-1 + 3*lmMst1 - 3*lmMst2)*Mst2*
        s2t*pow2(Mt) + 3*Mt*(-50 + 51*lmMst2 + 3*lmMst1*(-17 + 6*lmMst2) - 18*
        pow2(lmMst2))*pow2(Mst2)*pow2(s2t) + 4*(7 - 381*lmMst2 + lmMst1*(327 +
        72*lmMst2 - 81*lmMt) + 54*lmMt + 81*lmMst2*lmMt - 72*pow2(lmMst2))*
        pow3(Mt) + 2*(1 + 3*lmMst1 - 3*lmMst2)*pow3(Mst2)*pow3(s2t)) - Mst2*(-
        144*Mst2*s2t*(-lmMst1 + lmMst2 - 2*lmMst1*lmMst2 + pow2(lmMst1) + pow2(
        lmMst2))*pow2(Mt) + 24*Mt*(1 + 33*lmMst2 - 3*lmMst1*(11 + 6*lmMst2) +
        18*pow2(lmMst2))*pow2(Mst2)*pow2(s2t) + 32*(83 - 96*lmMst2 + 3*lmMst1*(
        29 + 12*lmMst2 - 9*lmMt) + 9*lmMt + 27*lmMst2*lmMt - 36*pow2(lmMst2))*
        pow3(Mt) + (5 + lmMst1*(6 - 288*lmMst2) - 6*lmMst2 + 144*pow2(lmMst1) +
        144*pow2(lmMst2))*pow3(Mst2)*pow3(s2t)))/Tbeta) + 576*Dmglst2*(4 + 6*
        lmMst1*(9 + 2*lmMst2 - 2*lmMt) + 7*lmMt + lmMst2*(-61 + 12*lmMt) - 12*
        pow2(lmMst2))*pow2(Mst2)*pow4(Mt) + 4*Dmglst2*(11 + 42*lmMst1 - 42*
        lmMst2)*Mt*pow3(s2t)*pow5(Mst2) + (-144*(lmMst1 - lmMst2)*Mst2*
        xDR2DRMOD*(2*Dmglst2*lmMst2*Mst2 + (-2 + lmMst2)*xDmglst2*pow2(Dmglst2)
        + (1 + lmMst2)*pow2(Mst2))*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*(-1 + pow2(
        Sbeta)) + (xDmglst2*pow2(Dmglst2)*(24*s2t*(4*(-143 + 18*lmMst1 - 18*
        lmMst2)*Tbeta*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + Mst2*pow2(Sbeta)*(-12*
        Mst2*Tbeta*(18 + lmMst2*(169 - 30*lmMt) - 3*lmMst1*(49 + 10*lmMst2 -
        10*lmMt) - 22*lmMt + 30*pow2(lmMst2)) + MuSUSY*(156 + 252*lmMst2 - 5*(
        43 + 60*lmMst2)*pow2(Sbeta) + 12*lmMst1*(-21 + 25*pow2(Sbeta)))))*pow3(
        Mt) - 2*Mst2*pow2(Mt)*pow2(s2t)*(-2*Tbeta*(157 + 348*lmMst2 + 12*
        lmMst1*(-29 + 3*lmMst2) - 36*pow2(lmMst2))*pow2(MuSUSY)*(-1 + pow2(
        Sbeta)) + 36*(85 - 36*lmMst1 + 36*lmMst2)*Mst2*MuSUSY*pow2(Sbeta) + 3*
        Tbeta*pow2(Mst2)*pow2(Sbeta)*(96 - 430*pow2(Sbeta) + 215*pow4(Sbeta) -
        12*(lmMst1 - lmMst2)*(22 - 50*pow2(Sbeta) + 25*pow4(Sbeta)))) + pow2(
        Sbeta)*(16*Mt*((1 - 24*lmMst1 + 24*lmMst2)*MuSUSY + 6*(1 - 6*lmMst1 +
        6*lmMst2)*Mst2*Tbeta)*pow3(Mst2)*pow3(s2t) + 32*(9*Mst2*Tbeta*(94 +
        475*lmMst2 - 6*lmMst1*(67 + 14*lmMst2 - 14*lmMt) - 73*lmMt - 84*lmMst2*
        lmMt + 84*pow2(lmMst2)) + MuSUSY*(398 + 2085*lmMst2 - 18*lmMst1*(95 +
        22*lmMst2 - 22*lmMt) - 375*lmMt - 396*lmMst2*lmMt + 396*pow2(lmMst2)))*
        pow4(Mt) + (5 + 6*lmMst1 - 6*lmMst2)*Tbeta*pow4(s2t)*pow5(Mst2))))/
        Tbeta)/pow2(Sbeta) + 2*(5 + 6*lmMst1 - 6*lmMst2)*(-2*Mt + Dmglst2*s2t)*
        pow3(s2t)*pow6(Mst2) - (11 + 6*lmMst1 - 6*lmMst2)*pow4(s2t)*pow7(Mst2))
        ))/(108.*pow9(Mst2)) + pow2(Al4p)*(threeLoopFlag*(-(pow2(Mt)*pow2(s2t)*
        ((-5*Dmglst2*Dmsqst2*(-1081 + 165*lmMst1 - 165*lmMst2)*(-2 + pow2(
        Sbeta))*pow2(Sbeta)*pow4(Mst1))/(9.*pow2(Msq)*pow3(Mst2)) + (pow2(
        MuSUSY)*(53.385802469135804 + (40*B4)/9. - (4*D3)/9. + (2*DN)/9. + (
        1672*lmMst1)/27. + (53*pow2(lmMst1))/9. - lmMst2*(129.92592592592592 -
        72*lmMst1 + (13*pow2(lmMst1))/3.) + (2*(-470 + 147*lmMst1)*pow2(lmMst2)
        )/9. + (Dmglst2*Mst2*(4.444444444444445 - (4*(99 + 16*lmMst1)*lmMst2)/
        9. - (82*pow2(lmMst2))/3. - (40*Dmsqst2)/pow2(Msq)))/pow2(Mst1) + ((5*
        Dmsqst2*(1141 + 48*lmMst1*(7 - 3*lmMst2) - 264*lmMst2 + 144*pow2(
        lmMst2)))/54. - ((Dmsqst2*(30 - 60*lmMst2) + (103 + 186*lmMst2 + 32*
        lmMst1*(1 + lmMst2) + 91*pow2(lmMst2))*pow2(Msq))*pow2(Mst2))/(9.*pow2(
        Mst1)))/pow2(Msq) + (16*pow3(lmMst1))/9. - (271*pow3(lmMst2))/9. - (
        Dmglst2*(4320*Dmsqst2*(4 + lmMst1 - lmMst2) + pow2(Msq)*(14267 - 432*B4
        + 576*D3 - 360*DN + 4752*lmMst1 - 1404*pow2(lmMst1) + 72*lmMst2*(-63 -
        232*lmMst1 + 16*pow2(lmMst1)) - 72*(-281 + 123*lmMst1)*pow2(lmMst2) +
        7704*pow3(lmMst2))))/(162.*Mst2*pow2(Msq)) + (pow2(Mst1)*(Mst2*(
        434.270658436214 + (76*B4)/9. - (2*DN)/9. + (69088*lmMst1)/405. - (
        1313*pow2(lmMst1))/27. - (4*lmMst2*(16192 - 26430*lmMst1 + 3465*pow2(
        lmMst1)))/405. + ((-5735 + 3072*lmMst1)*pow2(lmMst2))/27. + (Dmsqst2*(
        201.74098765432097 - (622*lmMst2)/15. + (2*lmMst1*(311 + 10*lmMst2))/
        15. - (22*pow2(lmMst1))/3. + 6*pow2(lmMst2)))/pow2(Msq) - (62*pow3(
        lmMst1))/27. - (2086*pow3(lmMst2))/27.) + Dmglst2*((8*B4)/3. + (10*
        Dmsqst2*(23 - 42*lmMst1 + 42*lmMst2))/(3.*pow2(Msq)) - (2*(2695042 +
        54000*D3 - 33750*DN - 326895*lmMst1 - 324900*pow2(lmMst1) + 15*lmMst2*(
        -19607 - 129030*lmMst1 + 62550*pow2(lmMst1)) - 450*(-5023 + 5610*
        lmMst1)*pow2(lmMst2) + 11250*pow3(lmMst1) + 1575000*pow3(lmMst2)))/
        30375.)))/pow3(Mst2) + (16*(1 + lmMst2)*(4*Dmglst2*lmMst2 + Mst2 +
        lmMst2*Mst2)*pow3(Mst2))/(9.*pow4(Mst1)) - ((-(Mst2*(628.1736268201578
         + (76*B4)/
        9. - (2*DN)/9. + (6317839*lmMst1)/396900. - (66307*pow2(lmMst1))/315. -
        lmMst2*(12.52907281431091 - (182909*lmMst1)/315. + (274*pow2(lmMst1))/
        3.) + (2*(-58301 + 37135*lmMst1)*pow2(lmMst2))/315. + (Dmsqst2*(
        237.28785508324435 + (16526*lmMst2)/3969. + (2*lmMst1*(-8263 + 71820*
        lmMst2))/3969. - (520*pow2(lmMst1))/21. - (80*pow2(lmMst2))/7.))/pow2(
        Msq) - (44*pow3(lmMst1))/9. - (1256*pow3(lmMst2))/9.)) + Dmglst2*(
        585.1892843532082 - (8*B4)/3. + (32*D3)/9. - (20*DN)/9. - (20109937*
        lmMst1)/297675. - (15886*pow2(lmMst1))/945. + lmMst2*(
        17.112243218274966 - (144628*lmMst1)/945. + (1180*pow2(lmMst1))/9.) + (
        2*(80257 - 138180*lmMst1)*pow2(lmMst2))/945. - (92*pow3(lmMst1))/27. +
        (4448*pow3(lmMst2))/27.))*pow4(Mst1))/pow5(Mst2) + (2*Dmglst2*pow2(Msq)
        *(54*(344*OepS2 + 9*(15643 - 774*lmMst1 + 774*lmMst2)*S2)*pow2(Mst1)*
        pow2(Mst2) + 4*(17308*OepS2 + 27*(93919 - 12981*lmMst1 + 12981*lmMst2)*
        S2)*pow4(Mst1) + 27*(56*OepS2 - 81*(-1677 + 14*lmMst1 - 14*lmMst2)*S2)*
        pow4(Mst2)) - 3*Mst2*(pow2(Msq)*(3*(8456*OepS2 - 81*(11243 + 2114*
        lmMst1 - 2114*lmMst2)*S2)*pow2(Mst1)*pow2(Mst2) + (52948*OepS2 - 27*(
        194357 + 39711*lmMst1 - 39711*lmMst2)*S2)*pow4(Mst1) + 27*(184*OepS2 -
        81*(307 + 46*lmMst1 - 46*lmMst2)*S2)*pow4(Mst2)) + 10*Dmsqst2*(9*(104*
        OepS2 + 27*(17 - 78*lmMst1 + 78*lmMst2)*S2)*pow2(Mst1)*pow2(Mst2) + 13*
        (136*OepS2 - 27*(95 + 102*lmMst1 - 102*lmMst2)*S2)*pow4(Mst1) + 27*(8*
        OepS2 - 81*(-15 + 2*lmMst1 - 2*lmMst2)*S2)*pow4(Mst2))))/(2187.*pow2(
        Msq)*pow5(Mst2))))/pow2(Sbeta))) + pow4(s2t)*(-(Mst2*pow2(Mst1)*(Mst2*(
        79.01365226337448 - B4/9. + (2*D3)/9. - DN/6. + (4372*lmMst1)/405. + (
        lmMst2*(16051 + 22980*lmMst1 - 5175*pow2(lmMst1)))/810. - (1631*pow2(
        lmMst1))/108. + ((-92 + 327*lmMst1)*pow2(lmMst2))/27. - (Dmsqst2*(
        3.2221604938271606 + lmMst1*(5.188888888888889 - 7*lmMst2) - (317*
        lmMst2)/90. + (11*pow2(lmMst1))/6. + (31*pow2(lmMst2))/6.))/pow2(Msq) -
        (79*pow3(lmMst1))/54. - (115*pow3(lmMst2))/27.) + Dmglst2*(
        0.7822304526748971 - (2*B4)/3. + (8*D3)/9. - (5*DN)/9. + (81193*lmMst1)
        /4050. + (137*pow2(lmMst1))/135. - (lmMst2*(81643 + 86970*lmMst1 +
        48150*pow2(lmMst1)))/4050. + ((4969 + 3840*lmMst1)*pow2(lmMst2))/270. +
        (5*Dmsqst2*(75 - 26*lmMst1 + 26*lmMst2))/(6.*pow2(Msq)) - (5*pow3(
        lmMst1))/27. - (58*pow3(lmMst2))/27.))) + (Dmglst2*(15707 - 432*B4 +
        576*D3 - 360*DN + 4752*lmMst1 - 1404*pow2(lmMst1) + 72*lmMst2*(-277 -
        264*lmMst1 + 16*pow2(lmMst1)) + 72*(142 - 123*lmMst1)*pow2(lmMst2) + (
        4320*Dmsqst2*(1 + lmMst1 - lmMst2))/pow2(Msq) + 7704*pow3(lmMst2))*
        pow3(Mst2))/648. + (46.745471895783595 + B4 + D3/9. - DN/9. + (
        34838747*lmMst1)/529200. + (50723*pow2(lmMst1))/1890. + (lmMst2*(-
        23468297 - 17276980*lmMst1 + 3601500*pow2(lmMst1)))/529200. + (4*(2941
        - 2415*lmMst1)*pow2(lmMst2))/945. + (Dmsqst2*(15.13649301931237 +
        lmMst1*(13.99649785840262 - (106*lmMst2)/21.) - (621671*lmMst2)/39690.
         + (53*pow2(lmMst1))/
        21. + (53*pow2(lmMst2))/21.))/pow2(Msq) - (10*pow3(lmMst1))/27. + (409*
        pow3(lmMst2))/108. + (Dmglst2*(79.58863384550371 + (1436147*lmMst1)/
        1.1907e6 + (1363*pow2(lmMst1))/315. + lmMst2*(6.96052994037121 - (101*
        lmMst1)/315. + (11*pow2(lmMst1))/3.) - (0.7285714285714285 + (11*
        lmMst1)/3.)*pow2(lmMst2) + (5*Dmsqst2*(76 - 249*lmMst1 + 249*lmMst2))/(
        54.*pow2(Msq)) - (11*pow3(lmMst1))/9. + (11*pow3(lmMst2))/9.))/Mst2)*
        pow4(Mst1) - ((30*Dmsqst2*(1213 + 48*lmMst1*(7 - 3*lmMst2) - 408*lmMst2
        + 144*pow2(lmMst2)) + pow2(Msq)*(25289 + 1440*B4 - 144*D3 + 72*DN +
        22368*lmMst1 + 1908*pow2(lmMst1) - 12*lmMst2*(2296 - 2136*lmMst1 + 117*
        pow2(lmMst1)) + 504*(-53 + 21*lmMst1)*pow2(lmMst2) + 576*pow3(lmMst1) -
        9756*pow3(lmMst2)))*pow4(Mst2))/(1296.*pow2(Msq)) - S2*(Dmglst2*(
        61.111111111111114 - 72*lmMst1 + 72*lmMst2)*Mst2*pow2(Mst1) + ((5717 +
        1286*lmMst1 - 1286*lmMst2 + (10*Dmsqst2*(253 + 42*lmMst1 - 42*lmMst2))/
        pow2(Msq))*pow2(Mst1)*pow2(Mst2))/12. + Dmglst2*(838.5 - 7*lmMst1 + 7*
        lmMst2)*pow3(Mst2) - (Dmglst2*(51635 + 25194*lmMst1 - 25194*lmMst2)*
        pow4(Mst1))/(162.*Mst2) + (5*(4*Dmsqst2*(163 + 42*lmMst1 - 42*lmMst2) +
        (3370 + 1077*lmMst1 - 1077*lmMst2)*pow2(Msq))*pow4(Mst1))/(108.*pow2(
        Msq)) + (3*(307 + 46*lmMst1 - 46*lmMst2)*pow4(Mst2))/4. + (15*Dmsqst2*(
        -15 + 2*lmMst1 - 2*lmMst2)*pow4(Mst2))/(2.*pow2(Msq))) + (((360*
        Dmglst2*Dmsqst2 + Dmsqst2*(30 - 60*lmMst2)*Mst2 + Mst2*(135 + 250*
        lmMst2 + 32*lmMst1*(1 + lmMst2) + 123*pow2(lmMst2))*pow2(Msq) + 2*
        Dmglst2*(-20 + (262 + 32*lmMst1)*lmMst2 + 187*pow2(lmMst2))*pow2(Msq))*
        pow5(Mst2))/(36.*pow2(Mst1)) + (OepS2*(60*Dmsqst2*Mst2*(63*pow2(Mst1)*
        pow2(Mst2) + 14*pow4(Mst1) + 27*pow4(Mst2)) + pow2(Msq)*(11574*pow2(
        Mst1)*pow3(Mst2) + 5385*Mst2*pow4(Mst1) - 4*Dmglst2*(1944*pow2(Mst1)*
        pow2(Mst2) + 4199*pow4(Mst1) + 189*pow4(Mst2)) + 3726*pow5(Mst2))))/(
        2187.*Mst2))/pow2(Msq)) - (pow2(MuSUSY)*pow3(Mt)*((-36*Mst2*pow2(Mst1)*
        (64800*Dmsqst2*Mst2*(Dmglst2*(8 - 15*lmMst1 + 15*lmMst2) + (1 + 3*
        lmMst1 - 3*lmMst2)*Mst2)*s2t + pow2(Msq)*(5*Mst2*(9*Mt*(10667 - 96*B4 +
        96*D3 - 48*DN - 3072*lmMst1 - 384*(-25 + 9*lmMst1)*lmMst2 - 384*(-1 +
        lmMst1 - lmMst2)*lmMt - 224*OepS2 + 324*(-43 + 14*lmMst1 - 14*lmMst2)*
        S2 - 384*(-13 + 4*lmMst1)*pow2(lmMst2) + 1536*pow3(lmMst2)) + 2*Mst2*
        s2t*(75569 + 13716*B4 - 54*DN - 33426*lmMst1 - 1088*OepS2 + 162*(169 +
        136*lmMst1 - 136*lmMst2)*S2 - 2376*pow2(lmMst1) + 54*lmMst2*(1427 -
        1012*lmMst1 + 16*pow2(lmMst1)) - 108*(-642 + 203*lmMst1)*pow2(lmMst2) +
        21060*pow3(lmMst2))) + Dmglst2*(30*Mt*(28405 - 288*B4 + 288*D3 - 144*DN
        + 5472*lmMst1 - 288*(-89 + 2*lmMst1)*lmMst2 + 576*(-2 + lmMst1 -
        lmMst2)*lmMt + 224*OepS2 - 324*(65 + 14*lmMst1 - 14*lmMst2)*S2 - 576*(-
        9 + 8*lmMst1)*pow2(lmMst2) + 4608*pow3(lmMst2)) + Mst2*s2t*(66761 +
        301320*B4 - 4860*DN - 205380*lmMst1 + 40480*OepS2 - 216*(2489 + 3795*
        lmMst1 - 3795*lmMst2)*S2 + 23760*pow2(lmMst1) + 180*lmMst2*(4993 -
        1956*lmMst1 + 48*pow2(lmMst1)) - 1080*(-482 + 331*lmMst1)*pow2(lmMst2)
        + 348840*pow3(lmMst2))))) + 27*(86400*Dmsqst2*Mst2*(Dmglst2*(2 + 5*
        lmMst1 - 5*lmMst2) + (-1 - lmMst1 + lmMst2)*Mst2)*s2t + pow2(Msq)*(-15*
        Mst2*(32*Mt*(436 - 6*B4 + 6*D3 - 3*DN - 48*lmMst1 + (408 - 96*lmMst1)*
        lmMst2 + 24*lmMt - 972*S2 - 48*(-4 + lmMst1)*pow2(lmMst2) + 48*pow3(
        lmMst2)) + Mst2*s2t*(28683 + 5280*B4 - 48*DN - 5952*lmMst1 - 224*OepS2
        + 324*(-1 + 14*lmMst1 - 14*lmMst2)*S2 - 768*pow2(lmMst1) - 192*lmMst2*(
        -214 + 73*lmMst1 + 4*pow2(lmMst1)) - 96*(-268 + 57*lmMst1)*pow2(lmMst2)
        + 6240*pow3(lmMst2))) - Dmglst2*(2880*Mt*(180 - 2*B4 + 2*D3 - DN + 16*
        lmMst1 + 144*lmMst2 - 216*S2 - 16*(-2 + lmMst1)*pow2(lmMst2) + 16*pow3(
        lmMst2)) + Mst2*s2t*(23917 + 188640*B4 - 3600*DN - 37440*lmMst1 + 5600*
        OepS2 - 324*(-453 + 350*lmMst1 - 350*lmMst2)*S2 + 11520*pow2(lmMst1) -
        2880*lmMst2*(-237 + 55*lmMst1 + 4*pow2(lmMst1)) - 1440*(-280 + 121*
        lmMst1)*pow2(lmMst2) + 185760*pow3(lmMst2)))))*pow3(Mst2))/pow2(Msq) -
        10*(9*Mt*(383185 - 2592*B4 + 2592*D3 - 1296*DN - 187704*lmMst1 - 17440*
        OepS2 + 648*(-57 + 545*lmMst1 - 545*lmMst2)*S2 - 7992*pow2(lmMst1) -
        216*lmMst2*(-1733 + 630*lmMst1 + 26*pow2(lmMst1)) - 216*(-859 + 246*
        lmMst1)*pow2(lmMst2) + 3456*lmMt*(3 + 6*lmMst2 - 2*lmMst1*(3 + lmMst2)
        + pow2(lmMst1) + pow2(lmMst2)) + 720*pow3(lmMst1) + 58032*pow3(lmMst2))
        + s2t*(-(Dmglst2*(2773621 - 1660176*B4 + 25272*DN + 2004408*lmMst1 -
        525472*OepS2 + 108*(123113 + 98526*lmMst1 - 98526*lmMst2)*S2 + 3888*
        pow2(lmMst1) - 144*lmMst2*(36802 - 11421*lmMst1 + 1728*pow2(lmMst1)) +
        167184*(-14 + 15*lmMst1)*pow2(lmMst2) - 31104*pow3(lmMst1) - 2227824*
        pow3(lmMst2))) + 3*Mst2*(1702429 + 257904*B4 - 648*DN - 748656*lmMst1 -
        30112*OepS2 + 108*(9185 + 5646*lmMst1 - 5646*lmMst2)*S2 + 41904*pow2(
        lmMst1) + 216*lmMst2*(5971 - 6106*lmMst1 + 576*pow2(lmMst1)) - 41904*(-
        34 + 15*lmMst1)*pow2(lmMst2) - 3456*pow3(lmMst1) + 507600*pow3(lmMst2))
        ))*pow4(Mst1) - (622080*(1 + lmMst2)*(-Dmglst2 + 3*Dmglst2*lmMst2 +
        Mst2 + lmMst2*Mst2)*s2t*pow6(Mst2))/pow2(Mst1)))/(43740.*pow2(Sbeta)*
        pow6(Mst2)) + pow2(Mt)*pow2(s2t)*(pow2(Mst2)*(312.7354938271605 + (8*
        D3)/9. - (8*DN)/9. + (806*lmMst1)/27. + (32*lmMt)/3. + (2*lmMst2*(2267
        + 90*lmMst1 - 117*pow2(lmMst1)))/27. + (124*pow2(lmMst1))/9. + (4*(133
        + 8*lmMst1)*pow2(lmMst2))/9. + (5*Dmsqst2*(391 + 32*lmMst1 - 128*
        lmMst2))/(18.*pow2(Msq)) + (32*pow3(lmMst1))/9. + (14*pow3(lmMst2))/9.)
        + pow2(Mst1)*(97.4135390946502 - (8*D3)/9. + (8*DN)/9. - (88984*lmMst1)
        /405. - (64*(lmMst1 - lmMst2)*lmMt)/3. - (1738*pow2(lmMst1))/27. + (2*
        lmMst2*(52322 - 4710*lmMst1 + 315*pow2(lmMst1)))/405. + (4*(479 - 237*
        lmMst1)*pow2(lmMst2))/27. - (Dmsqst2*(17.999506172839506 + lmMst1*(
        32.62222222222222 - 16*lmMst2) - (2068*lmMst2)/45. + 8*pow2(lmMst1) +
        8*pow2(lmMst2)))/pow2(Msq) - (260*pow3(lmMst1))/27. + (1166*pow3(
        lmMst2))/27.) + ((-4*Dmglst2*(16875*Dmsqst2*((113 - 6*lmMst1 + 6*
        lmMst2)*pow2(Mst1) + 46*pow2(Mst2)) - pow2(Msq)*(225*pow2(Mst2)*(11674
        - 120*B4 + 120*D3 - 60*DN + 690*lmMst1 + 345*pow2(lmMst1) - 5*lmMst2*(-
        3427 - 54*lmMst1 + 96*pow2(lmMst1)) + (4515 - 480*lmMst1)*pow2(lmMst2)
        + 960*pow3(lmMst2)) + pow2(Mst1)*(1216808 + 1164855*lmMst2 - 162000*(2
        + lmMst2)*lmMt + 225*(-341 + 30*lmMst2)*pow2(lmMst1) + 30*lmMst1*(34484
        - 16260*lmMst2 + 5400*lmMt - 21825*pow2(lmMst2)) + 365400*pow2(lmMst2)
        - 2250*pow3(lmMst1) + 650250*pow3(lmMst2)))))/(30375.*Mst2) + (2*(-360*
        Dmglst2*Dmsqst2 + 30*Dmsqst2*(-1 + 2*lmMst2)*Mst2 - Mst2*(87 + 154*
        lmMst2 + 32*lmMst1*(1 + lmMst2) + 75*pow2(lmMst2))*pow2(Msq) - 2*
        Dmglst2*(-52 + 2*(67 + 16*lmMst1)*lmMst2 + 91*pow2(lmMst2))*pow2(Msq))*
        pow3(Mst2))/(9.*pow2(Mst1)))/pow2(Msq) - ((Dmglst2*(175.16754355781114
         - (17578814*lmMst1)/
        33075. + lmMst2*(257.7056386999244 + (105592*lmMst1)/315. - (208*pow2(
        lmMst1))/9.) - (35576*pow2(lmMst1))/315. + (16*(-4376 + 2555*lmMst1)*
        pow2(lmMst2))/315. + (64*lmMt*(2 + 5*lmMst2 - lmMst1*(5 + 2*lmMst2) +
        pow2(lmMst1) + pow2(lmMst2)))/3. + (160*Dmsqst2*(41 - 12*lmMst1 + 12*
        lmMst2))/(27.*pow2(Msq)) + (16*pow3(lmMst1))/27. - (2896*pow3(lmMst2))/
        27.) - Mst2*(199.98139767323744 - (18614063*lmMst1)/66150. + lmMst2*(
        293.39173091458804 - (25514*lmMst1)/945. - (286*pow2(lmMst1))/9.) - (
        48143*pow2(lmMst1))/945. + (32*lmMt*(-2*lmMst1*(1 + lmMst2) + lmMst2*(2
        + lmMst2) + pow2(lmMst1)))/3. + (77.94391534391535 - (2*lmMst1)/9.)*
        pow2(lmMst2) + (Dmsqst2*(3.5881655848700444 + (470824*lmMst2)/19845. +
        (4*lmMst1*(-117706 + 34965*lmMst2))/19845. - (74*pow2(lmMst1))/21. - (
        74*pow2(lmMst2))/21.))/pow2(Msq) + (190*pow3(lmMst1))/27. + (674*pow3(
        lmMst2))/27.))*pow4(Mst1))/pow3(Mst2) + pow2(MuSUSY)*(
        53.385802469135804 + (40*B4)/9. - (4*D3)/9. + (2*DN)/9. + (1672*lmMst1)
        /27. + (53*pow2(lmMst1))/9. - lmMst2*(129.92592592592592 - 72*lmMst1 +
        (13*pow2(lmMst1))/3.) + (2*(-470 + 147*lmMst1)*pow2(lmMst2))/9. + (
        Dmglst2*Mst2*(4.444444444444445 - (4*(99 + 16*lmMst1)*lmMst2)/9. - (82*
        pow2(lmMst2))/3. - (40*Dmsqst2)/pow2(Msq)))/pow2(Mst1) + ((5*Dmsqst2*(
        1141 + 48*lmMst1*(7 - 3*lmMst2) - 264*lmMst2 + 144*pow2(lmMst2)))/54. -
        ((Dmsqst2*(30 - 60*lmMst2) + (103 + 186*lmMst2 + 32*lmMst1*(1 + lmMst2)
        + 91*pow2(lmMst2))*pow2(Msq))*pow2(Mst2))/(9.*pow2(Mst1)))/pow2(Msq) +
        (16*pow3(lmMst1))/9. - (271*pow3(lmMst2))/9. - (Dmglst2*(4320*Dmsqst2*(
        4 + lmMst1 - lmMst2) + pow2(Msq)*(14267 - 432*B4 + 576*D3 - 360*DN +
        4752*lmMst1 - 1404*pow2(lmMst1) + 72*lmMst2*(-63 - 232*lmMst1 + 16*
        pow2(lmMst1)) - 72*(-281 + 123*lmMst1)*pow2(lmMst2) + 7704*pow3(lmMst2)
        )))/(162.*Mst2*pow2(Msq)) - (pow2(Mst1)*(-(Mst2*(434.270658436214 + (
        76*B4)/9. - (2*DN)/9. + (69088*lmMst1)/405. - (1313*pow2(lmMst1))/27. -
        (4*lmMst2*(16192 - 26430*lmMst1 + 3465*pow2(lmMst1)))/405. + ((-5735 +
        3072*lmMst1)*pow2(lmMst2))/27. + (Dmsqst2*(201.74098765432097 - (622*
        lmMst2)/15. + (2*lmMst1*(311 + 10*lmMst2))/15. - (22*pow2(lmMst1))/3. +
        6*pow2(lmMst2)))/pow2(Msq) - (62*pow3(lmMst1))/27. - (2086*pow3(lmMst2)
        )/27.)) + (2*Dmglst2*(2695042 - 40500*B4 + 54000*D3 - 33750*DN -
        326895*lmMst1 - 324900*pow2(lmMst1) + 15*lmMst2*(-19607 - 129030*lmMst1
        + 62550*pow2(lmMst1)) + 450*(5023 - 5610*lmMst1)*pow2(lmMst2) + (50625*
        Dmsqst2*(-23 + 42*lmMst1 - 42*lmMst2))/pow2(Msq) + 11250*pow3(lmMst1) +
        1575000*pow3(lmMst2)))/30375.))/pow3(Mst2) + (16*(1 + lmMst2)*(4*
        Dmglst2*lmMst2 + Mst2 + lmMst2*Mst2)*pow3(Mst2))/(9.*pow4(Mst1)) - ((-(
        Mst2*(628.1736268201578 + (76*B4)/9. - (2*DN)/9. + (6317839*lmMst1)/
        396900. - (66307*pow2(lmMst1))/315. - lmMst2*(12.52907281431091 - (
        182909*lmMst1)/315. + (274*pow2(lmMst1))/3.) + (2*(-58301 + 37135*
        lmMst1)*pow2(lmMst2))/315. + (Dmsqst2*(237.28785508324435 + (16526*
        lmMst2)/3969. + (2*lmMst1*(-8263 + 71820*lmMst2))/3969. - (520*pow2(
        lmMst1))/21. - (80*pow2(lmMst2))/7.))/pow2(Msq) - (44*pow3(lmMst1))/9.
         - (1256*pow3(lmMst2))/
        9.)) + Dmglst2*(585.1892843532082 - (8*B4)/3. + (32*D3)/9. - (20*DN)/9.
         - (20109937*lmMst1)/297675. - (15886*pow2(lmMst1))/945. + lmMst2*(
        17.112243218274966 - (144628*lmMst1)/945. + (1180*pow2(lmMst1))/9.) + (
        169.85608465608465 - (2632*lmMst1)/9.)*pow2(lmMst2) - (92*pow3(lmMst1))
        /27. + (4448*pow3(lmMst2))/27.))*pow4(Mst1))/pow5(Mst2)) + (32*(1 +
        lmMst2)*(4*Dmglst2*lmMst2 + Mst2 + lmMst2*Mst2)*pow5(Mst2))/(9.*pow4(
        Mst1)) + (-3*Mst2*(pow2(Msq)*(6*pow2(Mst1)*pow2(Mst2)*((25160*OepS2 -
        81*(9191 + 6290*lmMst1 - 6290*lmMst2)*S2)*pow2(Mst2) + 5*(8456*OepS2 -
        81*(11243 + 2114*lmMst1 - 2114*lmMst2)*S2)*pow2(MuSUSY)) + ((274280*
        OepS2 - 27*(399127 + 205710*lmMst1 - 205710*lmMst2)*S2)*pow2(Mst2) +
        10*(52948*OepS2 - 27*(194357 + 39711*lmMst1 - 39711*lmMst2)*S2)*pow2(
        MuSUSY))*pow4(Mst1) + 27*((2120*OepS2 - 81*(-141 + 530*lmMst1 - 530*
        lmMst2)*S2)*pow2(Mst2) + 10*(184*OepS2 - 81*(307 + 46*lmMst1 - 46*
        lmMst2)*S2)*pow2(MuSUSY))*pow4(Mst2)) + 100*Dmsqst2*(9*pow2(Mst1)*pow2(
        Mst2)*(2*(8*OepS2 - 27*(23 + 6*lmMst1 - 6*lmMst2)*S2)*pow2(Mst2) + (
        104*OepS2 + 27*(17 - 78*lmMst1 + 78*lmMst2)*S2)*pow2(MuSUSY)) + (4*(16*
        OepS2 - 27*(35 + 12*lmMst1 - 12*lmMst2)*S2)*pow2(Mst2) + 13*(136*OepS2
        - 27*(95 + 102*lmMst1 - 102*lmMst2)*S2)*pow2(MuSUSY))*pow4(Mst1) + 27*(
        (8*OepS2 - 81*(9 + 2*lmMst1 - 2*lmMst2)*S2)*pow2(Mst2) + (8*OepS2 - 81*
        (-15 + 2*lmMst1 - 2*lmMst2)*S2)*pow2(MuSUSY))*pow4(Mst2))) + 4*Dmglst2*
        pow2(Msq)*(2*(9*(30760*OepS2 - 27*(28283 + 23070*lmMst1 - 23070*lmMst2)
        *S2)*pow2(Mst2) + 10*(17308*OepS2 + 27*(93919 - 12981*lmMst1 + 12981*
        lmMst2)*S2)*pow2(MuSUSY))*pow4(Mst1) + 135*(56*OepS2 - 81*(-1677 + 14*
        lmMst1 - 14*lmMst2)*S2)*pow2(MuSUSY)*pow4(Mst2) + 54*pow2(Mst1)*(5*(
        344*OepS2 + 9*(15643 - 774*lmMst1 + 774*lmMst2)*S2)*pow2(Mst2)*pow2(
        MuSUSY) + (2000*OepS2 - 162*(31 + 250*lmMst1 - 250*lmMst2)*S2)*pow4(
        Mst2)) - 2768742*S2*pow6(Mst2)))/(21870.*pow2(Msq)*pow5(Mst2))) + (Mt*
        pow3(s2t)*((583200*Dmsqst2*(4*(-1 + lmMst1 - lmMst2)*pow2(Mst1)*pow3(
        Mst2) + (-4 + lmMst1 - lmMst2)*Mst2*pow4(Mst1) + Dmglst2*(4*(12 - 5*
        lmMst1 + 5*lmMst2)*pow2(Mst1)*pow2(Mst2) + (40 - 3*lmMst1 + 3*lmMst2)*
        pow4(Mst1) - 4*(2 + 5*lmMst1 - 5*lmMst2)*pow4(Mst2)) + 4*(1 + lmMst1 -
        lmMst2)*pow5(Mst2)))/pow2(Msq) + (Dmglst2*(18*pow2(Mst2)*(27211 +
        36720*B4 + 1080*DN - 298440*lmMst1 + 64160*OepS2 - 108*(14033 + 12030*
        lmMst1 - 12030*lmMst2)*S2 + 12960*pow2(lmMst1) + 360*lmMst2*(-503 -
        636*lmMst1 + 144*pow2(lmMst1)) - 2160*(30 + 89*lmMst1)*pow2(lmMst2) +
        140400*pow3(lmMst2))*pow4(Mst1) + 27*pow2(Mst1)*(69997 + 188640*B4 -
        3600*DN - 37440*lmMst1 + 5600*OepS2 - 324*(-453 + 350*lmMst1 - 350*
        lmMst2)*S2 + 11520*pow2(lmMst1) - 2880*lmMst2*(-205 + 55*lmMst1 + 4*
        pow2(lmMst1)) - 1440*(-184 + 121*lmMst1)*pow2(lmMst2) + 185760*pow3(
        lmMst2))*pow4(Mst2) - (31897243 - 2491360*OepS2 + 90290268*S2 - 360*
        lmMst2*(18652 + 140139*S2) + 38880*(37 - 40*lmMst2)*pow2(lmMst1) +
        3188160*pow2(lmMst2) + 360*lmMst1*(17410 - 12852*lmMst2 + 140139*S2 +
        11232*pow2(lmMst2)) - 311040*pow3(lmMst1) - 2177280*pow3(lmMst2))*pow6(
        Mst1) + 622080*(-1 + 2*lmMst2 + 3*pow2(lmMst2))*pow6(Mst2)) + 15*Mst2*(
        6*pow2(Mst2)*(51041 + 7344*B4 + 216*DN - 80136*lmMst1 - 2336*OepS2 +
        324*(347 + 146*lmMst1 - 146*lmMst2)*S2 - 2592*pow2(lmMst1) + 216*
        lmMst2*(-221 - 428*lmMst1 + 48*pow2(lmMst1)) - 432*(-122 + 89*lmMst1)*
        pow2(lmMst2) + 28080*pow3(lmMst2))*pow4(Mst1) + 27*pow2(Mst1)*(25611 +
        5280*B4 - 48*DN - 5952*lmMst1 - 224*OepS2 + 324*(-1 + 14*lmMst1 - 14*
        lmMst2)*S2 - 768*pow2(lmMst1) - 192*lmMst2*(-182 + 73*lmMst1 + 4*pow2(
        lmMst1)) - 96*(-236 + 57*lmMst1)*pow2(lmMst2) + 6240*pow3(lmMst2))*
        pow4(Mst2) + (551987 - 14048*OepS2 + 661068*S2 - 216*lmMst2*(46 + 1317*
        S2) + 864*(205 + 216*lmMst2)*pow2(lmMst1) + 216000*pow2(lmMst2) - 216*
        lmMst1*(248 + 1820*lmMst2 - 1317*S2 + 1632*pow2(lmMst2)) - 6912*pow3(
        lmMst1) + 172800*pow3(lmMst2))*pow6(Mst1) + 41472*pow2(1 + lmMst2)*
        pow6(Mst2)))/pow2(Mst1)))/(87480.*pow2(Mst2)) + pow4(Mt)*(
        480.98395061728394 - (640*lmMst1)/9. - (392*pow2(lmMst1))/9. + (4*
        lmMst2*(-553 - 1224*lmMst1 + 63*pow2(lmMst1)))/81. + (32*(121 - 18*
        lmMst1)*pow2(lmMst2))/27. + (4*lmMt*(926 - 749*lmMst2 + 6*lmMst1*(27 +
        32*lmMst2) + 9*pow2(lmMst1) + 57*pow2(lmMst2)))/27. - (4*(-12 + lmMst1
        + 55*lmMst2)*pow2(lmMt))/3. + ((40*Dmsqst2*(-103 + 42*lmMst1 + 78*lmMt
        - 6*lmMst2*(23 + 3*lmMt) + 9*pow2(lmMst2) + 9*pow2(lmMt)))/27. + (4*(
        Dmsqst2*(30 - 60*lmMst2) + (71 + 122*lmMst2 + 32*lmMst1*(1 + lmMst2) +
        59*pow2(lmMst2))*pow2(Msq))*pow2(Mst2))/(9.*pow2(Mst1)))/pow2(Msq) - (
        68*pow3(lmMst1))/9. + (8*pow3(lmMst2))/9. + (184*pow3(lmMt))/3. - (64*(
        1 + lmMst2)*(4*Dmglst2*lmMst2 + Mst2 + lmMst2*Mst2)*pow3(Mst2))/(9.*
        pow4(Mst1)) + (600*Dmglst2*pow2(Mst2)*(350*Dmsqst2*((2954 - 864*lmMst1
        + 1713*lmMst2 - 849*lmMt)*pow2(Mst1) + 324*pow2(Mst2)) + pow2(Msq)*(
        630*(-84 + (70 + 32*lmMst1)*lmMst2 + 59*pow2(lmMst2))*pow2(Mst2) +
        pow2(Mst1)*(773533 - 229915*lmMst2 + 525*(251 + 148*lmMst2)*lmMt -
        4410*pow2(lmMst1) - 199920*pow2(lmMst2) + 5040*lmMst1*(13 + 28*lmMst2 -
        8*lmMt + 8*pow2(lmMst2)) + 63000*pow2(lmMt) - 40320*pow3(lmMst2)))) + (
        2520*Dmsqst2*Mst2*(-91321 + 30*lmMst1*(1717 + 120*lmMst2 - 150*lmMt) +
        40500*lmMt + 30*lmMst2*(-3067 + 150*lmMt) + 450*pow2(lmMst1) - 4050*
        pow2(lmMst2)) + 5670000*Dmglst2*Dmsqst2*(401 + lmMst2*(274 - 4*lmMt) -
        2*lmMst1*(71 + 2*lmMst2 - 2*lmMt) - 132*lmMt + 4*pow2(lmMst2)) + 4*
        Dmglst2*pow2(Msq)*(98884301 - 57865500*lmMt + 12600*(1213 + 285*lmMst2
        - 225*lmMt)*pow2(lmMst1) + 18900*(-5953 + 310*lmMt)*pow2(lmMst2) +
        13608000*pow2(lmMt) + 1260*(lmMst1*(28194 + lmMst2*(61415 - 2400*lmMt)
        - 1075*lmMt + 13800*pow2(lmMst2) - 7200*pow2(lmMt)) + lmMst2*(8581 +
        6025*lmMt + 7200*pow2(lmMt))) - 252000*pow3(lmMst1) - 20727000*pow3(
        lmMst2)) + 35*Mst2*pow2(Msq)*(7319011 + 3223800*lmMt - 1800*(391 + 114*
        lmMst2 - 90*lmMt)*pow2(lmMst1) + 3600*(568 + 99*lmMt)*pow2(lmMst2) -
        259200*pow2(lmMt) - 120*lmMst2*(22352 - 2835*lmMt + 4320*pow2(lmMt)) +
        120*lmMst1*(4217 + 1215*lmMt - 15*lmMst2*(871 + 288*lmMt) + 3240*pow2(
        lmMst2) + 4320*pow2(lmMt)) + 14400*pow3(lmMst1) - 198000*pow3(lmMst2)))
        *pow4(Mst1))/(425250.*pow2(Msq)*pow2(Mst1)*pow3(Mst2)) - ((-(Mst2*(
        1500.0856244066115 + (43574647*lmMst1)/99225. - (139432*pow2(lmMst1))/
        945. - (2*lmMst2*(35585111 - 2612085*lmMst1 + 2072700*pow2(lmMst1)))/
        99225. + (2*(34339 + 46620*lmMst1)*pow2(lmMst2))/945. + (lmMt*(3155 +
        1418*lmMst2 - 2*lmMst1*(513 + 572*lmMst2) + 312*pow2(lmMst1) + 832*
        pow2(lmMst2)))/9. + (64*(-1 + 3*lmMst1 - 3*lmMst2)*pow2(lmMt))/3. + (8*
        Dmsqst2*(-26331136 + 105*lmMst1*(236317 + 35070*lmMst2 - 36750*lmMt) +
        11962125*lmMt + 210*lmMst2*(-175121 + 18375*lmMt) + 88200*pow2(lmMst1)
        - 3770550*pow2(lmMst2)))/(231525.*pow2(Msq)) + (64*pow3(lmMst1))/27. -
        (1600*pow3(lmMst2))/27.)) + Dmglst2*(2929.938520304849 + (55510684*
        lmMst1)/59535. - (126272*pow2(lmMst1))/189. - (4*lmMst2*(42300121 +
        12578580*lmMst1 + 2487240*pow2(lmMst1)))/59535. + (32*(10166 - 693*
        lmMst1)*pow2(lmMst2))/189. + (8*lmMt*(5695 + 1974*lmMst2 - 12*lmMst1*(
        163 + 47*lmMst2) + 468*pow2(lmMst1) + 96*pow2(lmMst2)))/27. + (128*(-5
        + 6*lmMst1 - 6*lmMst2)*pow2(lmMt))/3. + (256*pow3(lmMst1))/27. + (7424*
        pow3(lmMst2))/27.))*pow4(Mst1))/pow5(Mst2) + (-(pow2(MuSUSY)*(6*Mst2*
        pow2(Mst1)*(3*Mst2*(10667 - 96*B4 + 96*D3 - 48*DN - 3072*lmMst1 - 384*(
        -25 + 9*lmMst1)*lmMst2 - 384*(-1 + lmMst1 - lmMst2)*lmMt - 384*(-13 +
        4*lmMst1)*pow2(lmMst2) + 1536*pow3(lmMst2)) + 2*Dmglst2*(28405 - 288*B4
        + 288*D3 - 144*DN + 5472*lmMst1 - 288*(-89 + 2*lmMst1)*lmMst2 + 576*(-2
        + lmMst1 - lmMst2)*lmMt - 576*(-9 + 8*lmMst1)*pow2(lmMst2) + 4608*pow3(
        lmMst2))) + 144*(6*Dmglst2*(180 - 2*B4 + 2*D3 - DN + 16*lmMst1 + 144*
        lmMst2 - 16*(-2 + lmMst1)*pow2(lmMst2) + 16*pow3(lmMst2)) + Mst2*(436 -
        6*B4 + 6*D3 - 3*DN - 48*lmMst1 + (408 - 96*lmMst1)*lmMst2 + 24*lmMt -
        48*(-4 + lmMst1)*pow2(lmMst2) + 48*pow3(lmMst2)))*pow3(Mst2) + (383185
        - 2592*B4 + 2592*D3 - 1296*DN - 187704*lmMst1 - 7992*pow2(lmMst1) -
        216*lmMst2*(-1733 + 630*lmMst1 + 26*pow2(lmMst1)) - 216*(-859 + 246*
        lmMst1)*pow2(lmMst2) + 3456*lmMt*(3 + 6*lmMst2 - 2*lmMst1*(3 + lmMst2)
        + pow2(lmMst1) + pow2(lmMst2)) + 720*pow3(lmMst1) + 58032*pow3(lmMst2))
        *pow4(Mst1)))/486. + (16*OepS2*(4*Dmglst2*Mst2*(63*pow2(Mst1)*(23*pow2(
        Mst2) - 3*pow2(MuSUSY)) + 5582*pow4(Mst1) + 189*pow4(Mst2)) + 3*((1168*
        pow2(Mst2) + 1635*pow2(MuSUSY))*pow4(Mst1) + pow2(Mst1)*(378*pow2(Mst2)
        *pow2(MuSUSY) + 627*pow4(Mst2)) + 189*pow6(Mst2))))/2187. - (S2*(4*
        Dmglst2*Mst2*(45*pow2(Mst1)*(23*(115 + 588*lmMst1 - 588*lmMst2)*pow2(
        Mst2) - 126*(65 + 14*lmMst1 - 14*lmMst2)*pow2(MuSUSY)) + (2001242 +
        2344440*lmMst1 - 2344440*lmMst2)*pow4(Mst1) - 81*(3360*pow2(Mst2)*pow2(
        MuSUSY) + (1593 - 980*lmMst1 + 980*lmMst2)*pow4(Mst2))) + 42*((8*(8653
        + 4380*lmMst1 - 4380*lmMst2)*pow2(Mst2) + 90*(-57 + 545*lmMst1 - 545*
        lmMst2)*pow2(MuSUSY))*pow4(Mst1) + 9*pow2(Mst1)*(90*(-43 + 14*lmMst1 -
        14*lmMst2)*pow2(Mst2)*pow2(MuSUSY) + (4163 + 2090*lmMst1 - 2090*lmMst2)
        *pow4(Mst2)) + 81*(-240*pow2(MuSUSY)*pow4(Mst2) + (109 + 70*lmMst1 -
        70*lmMst2)*pow6(Mst2)))))/2835.)/pow6(Mst2)) - (s2t*pow3(Mt)*(18*(
        1036800*Dmglst2*Dmsqst2*(7 - 3*lmMst1 + 6*lmMst2 - 3*lmMt) + 259200*
        Dmsqst2*(-3 + 2*lmMst1 - 4*lmMst2 + 2*lmMt)*Mst2 - 3*Dmglst2*pow2(Msq)*
        (322421 + 414720*lmMt - 2880*(21 + 12*lmMst2 - 4*lmMt)*pow2(lmMst1) +
        480*(1097 - 24*lmMt)*pow2(lmMst2) - 69120*pow2(lmMt) - 160*lmMst2*(-
        2852 + 231*lmMt + 216*pow2(lmMt)) + 160*lmMst1*(-2969 - 2211*lmMst2 -
        39*lmMt + 72*pow2(lmMst2) + 216*pow2(lmMt)) + 23040*pow3(lmMst2)) +
        Mst2*pow2(Msq)*(1156193 + 198720*lmMt - 6480*lmMst2*(21 + 16*lmMt*(1 +
        lmMt)) + 8640*(-13 + 4*lmMst2 + 4*lmMt)*pow2(lmMst1) + 17280*(51 - 2*
        lmMt)*pow2(lmMst2) - 2160*lmMst1*(295 + 360*lmMst2 - 52*lmMt + 112*
        pow2(lmMst2) - 48*pow2(lmMt)) + 207360*pow3(lmMst2)))*pow4(Mst1)*pow4(
        Mst2) + 2*Dmglst2*pow2(Mst1)*pow2(Mst2)*(32400*Dmsqst2*(9*(465 + 374*
        lmMst2 - 4*lmMst1*(52 + 3*lmMst2 - 3*lmMt) - 2*(83 + 6*lmMst2)*lmMt +
        12*pow2(lmMst2))*pow4(Mst1) + 2*(289 - 180*lmMst1 + 354*lmMst2 - 174*
        lmMt)*pow4(Mst2)) - pow2(Msq)*(4*(20964397 + 4563540*lmMt - 9720*(266 +
        74*lmMst2 - 55*lmMt)*pow2(lmMst1) + 3240*(1493 + 39*lmMt)*pow2(lmMst2)
        - 180*lmMst1*(32857 + 11844*lmMt + 18*lmMst2*(485 + 204*lmMt) - 1674*
        pow2(lmMst2) - 3888*pow2(lmMt)) + 1440*lmMst2*(1181 + 1332*lmMt - 486*
        pow2(lmMt)) - 466560*pow2(lmMt) + 9720*pow3(lmMst1) + 408240*pow3(
        lmMst2))*pow4(Mst1) + 135*(30183 + lmMst2*(67184 - 4704*lmMt) + 45408*
        lmMt + 2304*(-1 + lmMst2)*pow2(lmMst1) + 96*(199 + 48*lmMt)*pow2(
        lmMst2) + 288*lmMst1*(-30 + 16*lmMt - 2*lmMst2*(5 + 8*lmMt) + 41*pow2(
        lmMst2)) - 14112*pow3(lmMst2))*pow4(Mst2))) + pow2(Msq)*pow2(Mst1)*(
        108*S2*(54*pow2(Mst1)*((4367 + 3770*lmMst1 - 3770*lmMst2)*pow2(Mst2) +
        10*(169 + 136*lmMst1 - 136*lmMst2)*pow2(MuSUSY))*pow3(Mst2) + 6*Mst2*(
        8*(16297 + 10275*lmMst1 - 10275*lmMst2)*pow2(Mst2) + 5*(9185 + 5646*
        lmMst1 - 5646*lmMst2)*pow2(MuSUSY))*pow4(Mst1) - Dmglst2*(18*pow2(Mst1)
        *pow2(Mst2)*((15691 + 41490*lmMst1 - 41490*lmMst2)*pow2(Mst2) + 4*(2489
        + 3795*lmMst1 - 3795*lmMst2)*pow2(MuSUSY)) + 2*(8*(220117 + 192975*
        lmMst1 - 192975*lmMst2)*pow2(Mst2) + 5*(123113 + 98526*lmMst1 - 98526*
        lmMst2)*pow2(MuSUSY))*pow4(Mst1) - 81*(50*(27 - 14*lmMst1 + 14*lmMst2)*
        pow2(Mst2) + (453 - 350*lmMst1 + 350*lmMst2)*pow2(MuSUSY))*pow4(Mst2))
        - 243*(2*(79 - 70*lmMst1 + 70*lmMst2)*pow2(Mst2) + 5*(1 - 14*lmMst1 +
        14*lmMst2)*pow2(MuSUSY))*pow5(Mst2)) + 160*OepS2*(Dmglst2*(2*(51460*
        pow2(Mst2) + 16421*pow2(MuSUSY))*pow4(Mst1) + 945*(2*pow2(Mst2) + pow2(
        MuSUSY))*pow4(Mst2) + 18*pow2(Mst1)*(506*pow2(Mst2)*pow2(MuSUSY) +
        1383*pow4(Mst2))) - 3*(2*(941*Mst2*pow2(MuSUSY) + 2740*pow3(Mst2))*
        pow4(Mst1) + 189*(2*pow2(Mst2) + pow2(MuSUSY))*pow5(Mst2) + 6*pow2(
        Mst1)*(136*pow2(MuSUSY)*pow3(Mst2) + 377*pow5(Mst2))))) - 24*(24300*
        Dmsqst2*(35 + lmMst2*(50 - 4*lmMt) - 4*lmMst1*(8 + lmMst2 - lmMt) - 18*
        lmMt + 4*pow2(lmMst2)) - pow2(Msq)*(2016907 + 110160*lmMt + 1080*(565 -
        3*lmMt)*pow2(lmMst2) - 1080*((158 - 6*lmMst2 - 39*lmMt)*pow2(lmMst1) +
        lmMst1*(421 + 60*lmMt + lmMst2*(413 + 36*lmMt) + 129*pow2(lmMst2) - 72*
        pow2(lmMt)) + lmMst2*(167 - 66*lmMt + 72*pow2(lmMt))) + 1080*pow3(
        lmMst1) + 131760*pow3(lmMst2)))*pow3(Mst2)*pow6(Mst1) + 622080*(1 +
        lmMst2)*pow2(Msq)*(Dmglst2*(-1 + 3*lmMst2)*(2*pow2(Mst2) + pow2(MuSUSY)
        ) + 2*(1 + lmMst2)*pow3(Mst2))*pow6(Mst2) + pow2(MuSUSY)*(Dmglst2*pow2(
        Mst1)*(2332800*Dmsqst2*pow2(Mst2)*((8 - 15*lmMst1 + 15*lmMst2)*pow2(
        Mst1) + (-2 - 5*lmMst1 + 5*lmMst2)*pow2(Mst2)) - pow2(Msq)*(-36*pow2(
        Mst1)*pow2(Mst2)*(66761 + 301320*B4 - 4860*DN - 205380*lmMst1 + 23760*
        pow2(lmMst1) + 180*lmMst2*(4993 - 1956*lmMst1 + 48*pow2(lmMst1)) -
        1080*(-482 + 331*lmMst1)*pow2(lmMst2) + 348840*pow3(lmMst2)) + 10*(
        2773621 - 1660176*B4 + 25272*DN + 2004408*lmMst1 + 3888*pow2(lmMst1) -
        144*lmMst2*(36802 - 11421*lmMst1 + 1728*pow2(lmMst1)) + 167184*(-14 +
        15*lmMst1)*pow2(lmMst2) - 31104*pow3(lmMst1) - 2227824*pow3(lmMst2))*
        pow4(Mst1) - 27*(23917 + 188640*B4 - 3600*DN - 37440*lmMst1 + 11520*
        pow2(lmMst1) - 2880*lmMst2*(-237 + 55*lmMst1 + 4*pow2(lmMst1)) - 1440*(
        -280 + 121*lmMst1)*pow2(lmMst2) + 185760*pow3(lmMst2))*pow4(Mst2))) +
        15*(155520*Dmsqst2*pow2(Mst1)*((1 + 3*lmMst1 - 3*lmMst2)*pow2(Mst1) + (
        1 + lmMst1 - lmMst2)*pow2(Mst2))*pow3(Mst2) + Mst2*pow2(Msq)*(24*pow2(
        Mst2)*(75569 + 13716*B4 - 54*DN - 33426*lmMst1 - 2376*pow2(lmMst1) +
        54*lmMst2*(1427 - 1012*lmMst1 + 16*pow2(lmMst1)) - 108*(-642 + 203*
        lmMst1)*pow2(lmMst2) + 21060*pow3(lmMst2))*pow4(Mst1) + 81*pow2(Mst1)*(
        9561 + 1760*B4 - 16*DN - 1984*lmMst1 - 256*pow2(lmMst1) - 64*lmMst2*(-
        214 + 73*lmMst1 + 4*pow2(lmMst1)) - 32*(-268 + 57*lmMst1)*pow2(lmMst2)
        + 2080*pow3(lmMst2))*pow4(Mst2) + 2*(1702429 + 257904*B4 - 648*DN -
        748656*lmMst1 + 41904*pow2(lmMst1) + 216*lmMst2*(5971 - 6106*lmMst1 +
        576*pow2(lmMst1)) - 41904*(-34 + 15*lmMst1)*pow2(lmMst2) - 3456*pow3(
        lmMst1) + 507600*pow3(lmMst2))*pow6(Mst1) + 41472*pow2(1 + lmMst2)*
        pow6(Mst2)))) + 486*pow2(Mst1)*(2400*Dmsqst2*(7 + 4*lmMst1 - 10*lmMst2
        + 6*lmMt) + pow2(Msq)*(20253 - 4160*lmMt - 80*lmMst2*(13 + 66*lmMt) -
        1280*(1 + lmMst2)*pow2(lmMst1) + 160*(155 - 16*lmMt)*pow2(lmMst2) -
        160*lmMst1*(54 + 2*lmMst2*(61 - 8*lmMt) - 16*lmMt + 41*pow2(lmMst2)) -
        3840*pow2(lmMt) + 7840*pow3(lmMst2)))*pow7(Mst2)))/(43740.*pow2(Msq)*
        pow2(Mst1)*pow6(Mst2)) + (MuSUSY*(-(s2t*pow3(Mt)*(604.5820987654321 + (
        16*D3)/9. - (16*DN)/9. + (1228*lmMst1)/27. + (64*lmMt)/3. + (248*pow2(
        lmMst1))/9. - (4*lmMst2*(-1901 + 6*lmMst1 + 117*pow2(lmMst1)))/27. + (
        92 + (64*lmMst1)/9.)*pow2(lmMst2) + ((5*Dmsqst2*(367 + 32*lmMst1 - 80*
        lmMst2))/9. - (4*(Dmsqst2*(30 - 60*lmMst2) + (71 + 122*lmMst2 + 32*
        lmMst1*(1 + lmMst2) + 59*pow2(lmMst2))*pow2(Msq))*pow2(Mst2))/(9.*pow2(
        Mst1)))/pow2(Msq) + (64*pow3(lmMst1))/9. + (28*pow3(lmMst2))/9. + (64*(
        1 + lmMst2)*(4*Dmglst2*lmMst2 + Mst2 + lmMst2*Mst2)*pow3(Mst2))/(9.*
        pow4(Mst1)) + ((Mst2*(1199.3719723012073 - (99165049*lmMst1)/99225. +
        lmMst2*(1427.8402519526328 - (115988*lmMst1)/945. - (700*pow2(lmMst1))/
        9.) - (181826*pow2(lmMst1))/945. + (400.4804232804233 - (572*lmMst1)/9.
        )*pow2(lmMst2) + (64*lmMt*(1 + 4*lmMst2 - 2*lmMst1*(2 + lmMst2) + pow2(
        lmMst1) + pow2(lmMst2)))/3. + (Dmsqst2*(175.06620771294996 + (1883624*
        lmMst2)/19845. + (8*lmMst1*(-235453 + 114345*lmMst2))/19845. - (484*
        pow2(lmMst1))/21. - (484*pow2(lmMst2))/21.))/pow2(Msq) + (52*pow3(
        lmMst1))/27. + (3764*pow3(lmMst2))/27.) + Dmglst2*(622.8225754358181 -
        (64*B4)/9. + (64*D3)/9. - (32*DN)/9. + (29853268*lmMst1)/19845. + (
        38704*pow2(lmMst1))/189. + (4*lmMst2*(2917823 - 3813600*lmMst1 + 97020*
        pow2(lmMst1)))/19845. + (16*(8677 - 5439*lmMst1)*pow2(lmMst2))/189. - (
        128*lmMt*(4 + 6*lmMst2 - 2*lmMst1*(3 + lmMst2) + pow2(lmMst1) + pow2(
        lmMst2)))/3. - (16*pow3(lmMst1))/9. + (1328*pow3(lmMst2))/3.))*pow4(
        Mst1))/pow5(Mst2) + (-3*pow2(Mst1)*pow2(Mst2)*(32*Dmglst2*(50625*
        Dmsqst2*(65 - 2*lmMst1 + 2*lmMst2) - 2*pow2(Msq)*(1928479 - 13500*B4 +
        13500*D3 - 6750*DN + 635385*lmMst1 + 81000*(-2 + lmMst1 - lmMst2)*lmMt
        + 450*pow2(lmMst1) - 15*lmMst2*(-153166 + 17835*lmMst1 + 3375*pow2(
        lmMst1)) - 225*(-2627 + 1695*lmMst1)*pow2(lmMst2) - 1125*pow3(lmMst1) +
        433125*pow3(lmMst2))) - 5*Mst2*(12*Dmsqst2*(339977 + 96120*lmMst2 +
        1080*lmMst1*(-89 + 60*lmMst2) - 32400*pow2(lmMst1) - 32400*pow2(lmMst2)
        ) + pow2(Msq)*(19425643 + 518400*lmMt + 240*lmMst2*(82997 + 4320*lmMt)
        - 3600*(683 + 96*lmMst2)*pow2(lmMst1) + 5684400*pow2(lmMst2) - 240*
        lmMst1*(42047 + 4800*lmMst2 + 4320*lmMt + 6390*pow2(lmMst2)) - 295200*
        pow3(lmMst1) + 2174400*pow3(lmMst2)))) - (21600*Dmglst2*(150*Dmsqst2*(
        41*pow2(Mst1) + 18*pow2(Mst2)) - pow2(Msq)*(-15*(-52 + 2*(51 + 16*
        lmMst1)*lmMst2 + 59*pow2(lmMst2))*pow2(Mst2) + pow2(Mst1)*(12454 - 120*
        B4 + 120*D3 - 60*DN + 690*lmMst1 + 345*pow2(lmMst1) - 5*lmMst2*(-3121 +
        42*lmMst1 + 96*pow2(lmMst1)) + (3630 - 480*lmMst1)*pow2(lmMst2) + 960*
        pow3(lmMst2))))*pow4(Mst2))/pow2(Mst1) + 100*(24*Dmglst2*pow2(Msq)*(15*
        (400*OepS2 - 81*(139 + 100*lmMst1 - 100*lmMst2)*S2)*pow2(Mst1)*pow2(
        Mst2) + (36760*OepS2 - 54*(17269 + 13785*lmMst1 - 13785*lmMst2)*S2)*
        pow4(Mst1) - 153819*S2*pow4(Mst2)) - Mst2*(pow2(Msq)*(3*(69400*OepS2 -
        81*(17113 + 17350*lmMst1 - 17350*lmMst2)*S2)*pow2(Mst1)*pow2(Mst2) + 4*
        (120620*OepS2 - 27*(138286 + 90465*lmMst1 - 90465*lmMst2)*S2)*pow4(
        Mst1) + 27*(2120*OepS2 - 81*(-141 + 530*lmMst1 - 530*lmMst2)*S2)*pow4(
        Mst2)) + 100*Dmsqst2*(9*(40*OepS2 - 27*(127 + 30*lmMst1 - 30*lmMst2)*
        S2)*pow2(Mst1)*pow2(Mst2) + (424*OepS2 - 27*(1283 + 318*lmMst1 - 318*
        lmMst2)*S2)*pow4(Mst1) + 27*(8*OepS2 - 81*(9 + 2*lmMst1 - 2*lmMst2)*S2)
        *pow4(Mst2)))))/(364500.*pow2(Msq)*pow5(Mst2)))) - pow2(Mt)*pow2(s2t)*(
        (64*(1 + lmMst2)*(-Dmglst2 + 3*Dmglst2*lmMst2 + Mst2 + lmMst2*Mst2)*
        pow2(Mst2))/(3.*pow2(Mst1)) + Mst2*(377.0416666666667 + (220*B4)/3. - (
        2*DN)/3. - (248*lmMst1)/3. - (32*pow2(lmMst1))/3. - (8*lmMst2*(-198 +
        73*lmMst1 + 4*pow2(lmMst1)))/3. + (336 - 76*lmMst1)*pow2(lmMst2) + (80*
        Dmsqst2*(1 + lmMst1 - lmMst2))/pow2(Msq) + (260*pow3(lmMst2))/3.) + (
        pow2(Mst1)*(Mst2*(534.5756172839506 + 96*B4 + (1142*lmMst2)/3. + (8*(-7
        + 8*lmMst2)*pow2(lmMst1))/3. + (1496*pow2(lmMst2))/3. - (2*lmMst1*(495
        + 720*lmMst2 + 292*pow2(lmMst2)))/3. + (160*Dmsqst2*(lmMst1 - lmMst2))/
        pow2(Msq) + (520*pow3(lmMst2))/3.) + Dmglst2*(60.275617283950616 + (
        592*B4)/3. - (8*DN)/3. - (1970*lmMst1)/9. + (56*pow2(lmMst1))/3. + (2*
        lmMst2*(2149 - 1296*lmMst1 + 96*pow2(lmMst1)))/9. + (269.3333333333333
         - 280*lmMst1)*pow2(lmMst2) + (800*Dmsqst2*(1 - lmMst1 + lmMst2))/pow2(
        Msq) + (776*pow3(lmMst2))/3.)))/pow2(Mst2) + ((818.5195473251028 + 96*
        B4 - (3218*lmMst1)/9. + (652*pow2(lmMst1))/9. + (4*lmMst2*(845 - 1535*
        lmMst1 + 264*pow2(lmMst1)))/9. + (609.7777777777778 - 376*lmMst1)*pow2(
        lmMst2) + (20*Dmsqst2*(-4 + 9*lmMst1 - 9*lmMst2))/pow2(Msq) - (32*pow3(
        lmMst1))/9. + (2360*pow3(lmMst2))/9.)*pow4(Mst1))/pow3(Mst2) + S2*((9*(
        -1 + 14*lmMst1 - 14*lmMst2)*Mst2)/2. - ((Dmglst2*(799.6333333333333 +
        907*lmMst1 - 907*lmMst2) - ((685 + 418*lmMst1 - 418*lmMst2)*Mst2)/2.)*
        pow2(Mst1))/pow2(Mst2) + ((6143 + 3198*lmMst1 - 3198*lmMst2)*pow4(Mst1)
        )/(9.*pow3(Mst2)) + Dmglst2*(135.9 - 105*lmMst1 + 105*lmMst2 - ((525961
        + 356010*lmMst1 - 356010*lmMst2)*pow4(Mst1))/(135.*pow4(Mst2)))) +
        Dmglst2*(43.4787037037037 + (524*B4)/3. - (10*DN)/3. - (104*lmMst1)/3.
         + (32*pow2(lmMst1))/
        3. - (8*lmMst2*(-221 + 55*lmMst1 + 4*pow2(lmMst1)))/3. + (4*(232 - 121*
        lmMst1)*pow2(lmMst2))/3. - (80*Dmsqst2*(2 + 5*lmMst1 - 5*lmMst2))/pow2(
        Msq) + 172*pow3(lmMst2) - ((1033.594170096022 - (592*B4)/3. + (8*DN)/3.
         + (35140*lmMst1)/81. + (92*pow2(lmMst1))/3. - (2*lmMst2*(28667 - 5238*
        lmMst1 + 3024*pow2(lmMst1)))/81. + (8*(-60 + 157*lmMst1)*pow2(lmMst2))/
        3. + (20*Dmsqst2*(-80 + 43*lmMst1 - 43*lmMst2))/pow2(Msq) - (32*pow3(
        lmMst1))/3. - (1000*pow3(lmMst2))/3.)*pow4(Mst1))/pow4(Mst2)) + (4*
        OepS2*(Dmglst2*(8163*pow2(Mst1)*pow2(Mst2) + 23734*pow4(Mst1) + 945*
        pow4(Mst2)) - 3*(627*pow2(Mst1)*pow3(Mst2) + 1066*Mst2*pow4(Mst1) +
        189*pow5(Mst2))))/(729.*pow4(Mst2))) + (pow4(Mt)*((129600*Dmsqst2*pow2(
        Mst2)*(Dmglst2*(497 - 306*lmMst1 + 609*lmMst2 - 303*lmMt)*pow2(Mst1) +
        9*(1 + 6*lmMst1 - 13*lmMst2 + 7*lmMt)*Mst2*pow2(Mst1) + Dmglst2*(101 -
        90*lmMst1 + 177*lmMst2 - 87*lmMt)*pow2(Mst2) + 9*(5 + 2*lmMst1 - 5*
        lmMst2 + 3*lmMt)*pow3(Mst2)))/pow2(Msq) + (144*(110219 - 1080*lmMt -
        4400*OepS2 + 74034*S2 + 540*(-15 + 4*lmMt)*pow2(lmMst1) + 810*(121 - 8*
        lmMt)*pow2(lmMst2) - 135*lmMst1*(337 + lmMst2*(588 - 32*lmMt) - 132*
        lmMt - 660*S2 + 194*pow2(lmMst2) - 48*pow2(lmMt)) - 6480*pow2(lmMt) -
        405*lmMst2*(31 + 54*lmMt + 220*S2 + 16*pow2(lmMt)) + 26190*pow3(lmMst2)
        )*pow3(Mst2)*pow4(Mst1) + 81*pow2(Mst1)*(56439 - 24000*lmMt - 1120*
        OepS2 - 25596*S2 - 360*lmMst2*(-12 + 44*lmMt + 63*S2) - 3840*(1 +
        lmMst2)*pow2(lmMst1) + 480*(163 - 16*lmMt)*pow2(lmMst2) - 120*lmMst1*(
        184 + 8*lmMst2*(57 - 8*lmMt) - 64*lmMt - 189*S2 + 164*pow2(lmMst2)) -
        11520*pow2(lmMt) + 23520*pow3(lmMst2))*pow5(Mst2) + 60*Mst2*(678923 +
        19440*lmMt - 32480*OepS2 + 881712*S2 + 108*(-457 + 12*lmMst2 + 126*
        lmMt)*pow2(lmMst1) + 540*(661 - 30*lmMt)*pow2(lmMst2) - 108*lmMst1*(
        1841 + lmMst2*(2626 - 24*lmMt) - 420*lmMt - 6090*S2 + 840*pow2(lmMst2)
        - 288*pow2(lmMt)) - 15552*pow2(lmMt) - 108*lmMst2*(619 + 498*lmMt +
        6090*S2 + 288*pow2(lmMt)) + 216*pow3(lmMst1) + 89208*pow3(lmMst2))*
        pow6(Mst1) - Dmglst2*(216*pow2(Mst2)*(97837 + 71580*lmMt - 9920*OepS2 +
        43272*S2 - 360*(23 + 8*lmMst2 - 4*lmMt)*pow2(lmMst1) + 720*(97 + 2*
        lmMt)*pow2(lmMst2) - 8640*pow2(lmMt) - 30*lmMst2*(-2911 + 396*lmMt +
        6696*S2 + 144*pow2(lmMt)) + 10*lmMst1*(-6157 + 642*lmMt - 6*lmMst2*(791
        + 48*lmMt) + 20088*S2 + 882*pow2(lmMst2) + 432*pow2(lmMt)) - 5940*pow3(
        lmMst2))*pow4(Mst1) + 135*pow2(Mst1)*(57495 + 45408*lmMt - 1120*OepS2 -
        43740*S2 - 8*lmMst2*(-6952 + 588*lmMt + 2835*S2) + 2304*(-1 + lmMst2)*
        pow2(lmMst1) + 96*(79 + 48*lmMt)*pow2(lmMst2) + 72*lmMst1*(-88 + 64*
        lmMt - 8*lmMst2*(9 + 8*lmMt) + 315*S2 + 164*pow2(lmMst2)) - 14112*pow3(
        lmMst2))*pow4(Mst2) + 20*(5659217 + 1592460*lmMt - 518816*OepS2 +
        9976392*S2 - 972*(569 + 180*lmMst2 - 126*lmMt)*pow2(lmMst1) + 324*(5353
        + 126*lmMt)*pow2(lmMst2) - 186624*pow2(lmMt) + 72*lmMst1*(-27653 -
        3015*lmMt - 18*lmMst2*(689 + 126*lmMt) + 145917*S2 + 2160*pow2(lmMst2)
        + 2592*pow2(lmMt)) - 36*lmMst2*(-39031 - 3204*lmMt + 291834*S2 + 5184*
        pow2(lmMt)) + 1944*pow3(lmMst1) + 17496*pow3(lmMst2))*pow6(Mst1) -
        622080*(-1 + 2*lmMst2 + 3*pow2(lmMst2))*pow6(Mst2)) + 622080*pow2(1 +
        lmMst2)*pow7(Mst2))/pow2(Mst1)))/(21870.*pow6(Mst2)) + Mt*(-(pow2(Mst2)
        *(-2*s2t*(5*Dmsqst2*(3 - 2*lmMst2)*shiftst1*pow2(Mst2)*pow2(s2t) + (1 -
        2*lmMst2)*pow2(Msq)*(-2*shiftst3*pow2(Mt) + 5*shiftst1*pow2(Mst2)*pow2(
        s2t))) + (1 - 2*lmMst2)*shiftst3*pow2(Msq)*((1 + lmMst1 - lmMst2)*pow2(
        Mst1) - pow2(Mst2))*pow3(s2t)) + 10*s2t*shiftst1*(Dmsqst2*(3 - 2*
        lmMst2) + (1 - 2*lmMst2)*pow2(Msq))*(4*pow2(Mst2)*pow2(Mt) - 2*pow2(
        Mst1)*(2*pow2(Mt) + (-lmMst1 + lmMst2)*pow2(Mst2)*pow2(s2t)) + pow2(
        s2t)*pow4(Mst1)))/(3.*pow2(Msq)*pow2(Mst1)) + pow3(s2t)*((pow2(Mst2)*(
        30*Dmsqst2*(1177 + 48*lmMst1*(7 - 3*lmMst2) - 336*lmMst2 + 144*pow2(
        lmMst2)) + pow2(Msq)*(21005 + 1440*B4 - 144*D3 + 72*DN + 21216*lmMst1 +
        1908*pow2(lmMst1) - 12*lmMst2*(2950 - 2040*lmMst1 + 117*pow2(lmMst1)) +
        108*(-283 + 98*lmMst1)*pow2(lmMst2) + 576*pow3(lmMst1) - 9756*pow3(
        lmMst2))))/(324.*pow2(Msq)) + (pow2(Mst1)*(4627751 + 48600*B4 + 5400*D3
        - 5400*DN + 1320240*lmMst1 - 662400*pow2(lmMst1) - 30*lmMst2*(12148 -
        76560*lmMst1 + 12105*pow2(lmMst1)) + 2250*(-583 + 438*lmMst1)*pow2(
        lmMst2) + (12*Dmsqst2*(97294 - 17235*lmMst2 + 45*lmMst1*(233 + 330*
        lmMst2) - 7425*pow2(lmMst1) - 7425*pow2(lmMst2)))/pow2(Msq) - 49500*
        pow3(lmMst1) - 572850*pow3(lmMst2)))/12150. - (Dmglst2*((202500*
        Dmsqst2*((-55 + 34*lmMst1 - 34*lmMst2)*pow2(Mst1) + 4*(5 + 2*lmMst1 -
        2*lmMst2)*pow2(Mst2)))/pow2(Msq) + 375*pow2(Mst2)*(14987 - 432*B4 +
        576*D3 - 360*DN + 4752*lmMst1 - 1404*pow2(lmMst1) + 144*lmMst2*(-81 -
        124*lmMst1 + 8*pow2(lmMst1)) - 36*(-439 + 246*lmMst1)*pow2(lmMst2) +
        7704*pow3(lmMst2)) + pow2(Mst1)*(5430043 + 524580*lmMst2 + 900*(-859 +
        3690*lmMst2)*pow2(lmMst1) + 1454400*pow2(lmMst2) - 60*lmMst1*(51493 +
        24630*lmMst2 + 112950*pow2(lmMst2)) + 45000*pow3(lmMst1) + 3411000*
        pow3(lmMst2))))/(60750.*Mst2) - ((360*Dmglst2*Dmsqst2 + Dmsqst2*(30 -
        60*lmMst2)*Mst2 + Mst2*(119 + 218*lmMst2 + 32*lmMst1*(1 + lmMst2) +
        107*pow2(lmMst2))*pow2(Msq) + 2*Dmglst2*(-20 + (230 + 32*lmMst1)*lmMst2
        + 155*pow2(lmMst2))*pow2(Msq))*pow3(Mst2))/(9.*pow2(Msq)*pow2(Mst1)) -
        ((-(Mst2*(193.90296838394383 - (61388401*lmMst1)/396900. + lmMst2*(
        147.3919148400101 + (302047*lmMst1)/945. - (514*pow2(lmMst1))/9.) - (
        152966*pow2(lmMst1))/945. - (157.75767195767196 - 122*lmMst1)*pow2(
        lmMst2) + (Dmsqst2*(35.54686742892336 + (905536*lmMst2)/19845. +
        lmMst1*(-45.630435878054925 + (244*lmMst2)/7.) - (122*pow2(lmMst1))/7.
         - (122*pow2(lmMst2))/
        7.))/pow2(Msq) - (70*pow3(lmMst1))/27. - (1682*pow3(lmMst2))/27.)) +
        Dmglst2*(407.7379592503276 - (2740559*lmMst1)/59535. + (866*pow2(
        lmMst1))/189. + lmMst2*(36.477181489879904 - (4840*lmMst1)/189. + (208*
        pow2(lmMst1))/3.) + (21.026455026455025 - (1136*lmMst1)/9.)*pow2(
        lmMst2) + (10*Dmsqst2*(-419 + 57*lmMst1 - 57*lmMst2))/(27.*pow2(Msq)) -
        (112*pow3(lmMst1))/27. + (1648*pow3(lmMst2))/27.))*pow4(Mst1))/pow3(
        Mst2) + (2*Dmglst2*pow2(Msq)*(27*(632*OepS2 + 9*(16193 - 1422*lmMst1 +
        1422*lmMst2)*S2)*pow2(Mst1)*pow2(Mst2) + 2*(25328*OepS2 + 27*(47051 -
        18996*lmMst1 + 18996*lmMst2)*S2)*pow4(Mst1) + 27*(56*OepS2 - 81*(-1677
        + 14*lmMst1 - 14*lmMst2)*S2)*pow4(Mst2)) - 3*Mst2*(pow2(Msq)*(60*(340*
        OepS2 - 81*(424 + 85*lmMst1 - 85*lmMst2)*S2)*pow2(Mst1)*pow2(Mst2) +
        35*(788*OepS2 - 27*(2662 + 591*lmMst1 - 591*lmMst2)*S2)*pow4(Mst1) +
        27*(184*OepS2 - 81*(307 + 46*lmMst1 - 46*lmMst2)*S2)*pow4(Mst2)) + 10*
        Dmsqst2*(18*(40*OepS2 - 27*(59 + 30*lmMst1 - 30*lmMst2)*S2)*pow2(Mst1)*
        pow2(Mst2) + 4*(208*OepS2 - 27*(347 + 156*lmMst1 - 156*lmMst2)*S2)*
        pow4(Mst1) + 27*(8*OepS2 - 81*(-15 + 2*lmMst1 - 2*lmMst2)*S2)*pow4(
        Mst2))))/(2187.*pow2(Msq)*pow3(Mst2)) + (16*(1 + lmMst2)*(4*Dmglst2*
        lmMst2 + Mst2 + lmMst2*Mst2)*pow5(Mst2))/(9.*pow4(Mst1))) - (2*T1ep*(3*
        Mst2*(-(pow2(Msq)*(3*pow2(Mst1)*pow2(Mst2)*(-3470*Mst2*s2t*pow2(Mt) -
        627*Mt*pow2(Mst2)*pow2(s2t) + 1760*pow3(Mt) + 1700*pow3(Mst2)*pow3(s2t)
        ) + (-24124*Mst2*s2t*pow2(Mt) - 3198*Mt*pow2(Mst2)*pow2(s2t) + 16240*
        pow3(Mt) + 6895*pow3(Mst2)*pow3(s2t))*pow4(Mst1) + 27*(-106*Mst2*s2t*
        pow2(Mt) - 21*Mt*pow2(Mst2)*pow2(s2t) + 28*pow3(Mt) + 46*pow3(Mst2)*
        pow3(s2t))*pow4(Mst2))) + 20*Dmsqst2*Mst2*s2t*(2*(53*pow2(Mt) - 52*
        pow2(Mst2)*pow2(s2t))*pow4(Mst1) + 54*pow2(Mt)*pow4(Mst2) + 90*pow2(
        Mst1)*(pow2(Mst2)*pow2(Mt) - pow2(s2t)*pow4(Mst2)) - 27*pow2(s2t)*pow6(
        Mst2))) + Dmglst2*pow2(Msq)*(27*pow2(Mst1)*pow2(Mst2)*(-800*Mst2*s2t*
        pow2(Mt) - 907*Mt*pow2(Mst2)*pow2(s2t) + 1984*pow3(Mt) + 316*pow3(Mst2)
        *pow3(s2t)) + 2*(-66168*Mst2*s2t*pow2(Mt) - 35601*Mt*pow2(Mst2)*pow2(
        s2t) + 129704*pow3(Mt) + 12664*pow3(Mst2)*pow3(s2t))*pow4(Mst1) + 189*(
        20*pow3(Mt)*pow4(Mst2) - 15*Mt*pow2(s2t)*pow6(Mst2) + 4*pow3(s2t)*pow7(
        Mst2)))))/(729.*pow2(Msq)*pow6(Mst2)))))/Tbeta - (T1ep*(pow2(Msq)*(
        pow4(Mst1)*(4*Mst2*pow2(Mt)*pow2(s2t)*(34616*Dmglst2*pow2(MuSUSY)*(-1 +
        pow2(Sbeta)) - 39711*Mst2*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 55368*
        Dmglst2*pow2(Mst2)*pow2(Sbeta) - 20571*pow2(Sbeta)*pow3(Mst2)) + 16*
        s2t*(-16421*Dmglst2*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 2823*Mst2*pow2(
        MuSUSY)*(-1 + pow2(Sbeta)) - 51460*Dmglst2*pow2(Mst2)*pow2(Sbeta) +
        8220*pow2(Sbeta)*pow3(Mst2))*pow3(Mt) + 4*(15571*Dmglst2 - 1317*Mst2)*
        Mt*pow2(Sbeta)*pow3(s2t)*pow4(Mst2) + 16*(4905*pow2(MuSUSY)*(-1 + pow2(
        Sbeta)) + 8*Mst2*(2791*Dmglst2 + 438*Mst2)*pow2(Sbeta))*pow4(Mt) + (-
        16796*Dmglst2 + 5385*Mst2)*pow2(Sbeta)*pow4(s2t)*pow5(Mst2)) + 54*pow4(
        Mst2)*(14*Dmglst2*(4*Mst2*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*(-1 + pow2(
        Sbeta)) - 10*s2t*(pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 2*pow2(Mst2)*pow2(
        Sbeta))*pow3(Mt) + 5*Mt*pow2(Sbeta)*pow3(s2t)*pow4(Mst2) + 16*Mst2*
        pow2(Sbeta)*pow4(Mt) - pow2(Sbeta)*pow4(s2t)*pow5(Mst2)) + 3*Mst2*(-2*
        Mst2*pow2(Mt)*pow2(s2t)*(46*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 53*pow2(
        Mst2)*pow2(Sbeta)) + 28*s2t*(pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 2*pow2(
        Mst2)*pow2(Sbeta))*pow3(Mt) - 14*Mt*pow2(Sbeta)*pow3(s2t)*pow4(Mst2) +
        56*Mst2*pow2(Sbeta)*pow4(Mt) + 23*pow2(Sbeta)*pow4(s2t)*pow5(Mst2))) +
        18*Mst2*pow2(Mst1)*(-4*Dmglst2*(-12*pow2(Mst2)*pow2(Mt)*pow2(s2t)*(43*
        pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 50*pow2(Mst2)*pow2(Sbeta)) + 2*Mst2*
        s2t*(506*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 1383*pow2(Mst2)*pow2(Sbeta))
        *pow3(Mt) + 56*(3*pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 23*pow2(Mst2)*pow2(
        Sbeta))*pow4(Mt) - 401*Mt*pow2(Sbeta)*pow3(s2t)*pow5(Mst2) + 108*pow2(
        Sbeta)*pow4(s2t)*pow6(Mst2)) + Mst2*(-4*pow2(Mst2)*pow2(Mt)*pow2(s2t)*(
        1057*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 629*pow2(Mst2)*pow2(Sbeta)) + 8*
        Mst2*s2t*(136*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 377*pow2(Mst2)*pow2(
        Sbeta))*pow3(Mt) + 8*(126*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 209*pow2(
        Mst2)*pow2(Sbeta))*pow4(Mt) - 292*Mt*pow2(Sbeta)*pow3(s2t)*pow5(Mst2) +
        643*pow2(Sbeta)*pow4(s2t)*pow6(Mst2)))) + 60*Dmsqst2*pow2(Mst2)*pow2(
        s2t)*(pow4(Mst1)*(-4*pow2(Mt)*(221*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 8*
        pow2(Mst2)*pow2(Sbeta)) + 14*pow2(s2t)*pow2(Sbeta)*pow4(Mst2)) + 9*
        pow2(Mst1)*(pow2(Mst2)*pow2(Mt)*(-52*pow2(MuSUSY)*(-1 + pow2(Sbeta)) -
        8*pow2(Mst2)*pow2(Sbeta)) + 7*pow2(s2t)*pow2(Sbeta)*pow6(Mst2)) + 27*(-
        4*pow2(Mt)*(pow2(MuSUSY)*(-1 + pow2(Sbeta)) + pow2(Mst2)*pow2(Sbeta))*
        pow4(Mst2) + pow2(s2t)*pow2(Sbeta)*pow8(Mst2)))))/(1458.*pow2(Msq)*
        pow2(Sbeta)*pow6(Mst2)) + ((5*shiftst1*(Dmsqst2*(3 - 2*lmMst2) + (1 -
        2*lmMst2)*pow2(Msq))*(-(pow4(Mst2)*(pow2(s2t)*(-8*pow2(Mt) + (pow2(
        Mst1) - pow2(Mst2))*pow2(s2t))*pow2(Sbeta)*pow4(Mst1) + pow2(Mst1)*(4*
        pow2(Mt)*pow2(s2t)*(-(pow2(MuSUSY)*(-1 + pow2(Sbeta))) + 4*pow2(Mst2)*
        pow2(Sbeta)) + pow2(Sbeta)*(16*pow4(Mt) - pow4(Mst2)*pow4(s2t))) +
        pow2(Mst2)*(-4*pow2(Mt)*pow2(s2t)*(pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 2*
        pow2(Mst2)*pow2(Sbeta)) + pow2(Sbeta)*(16*pow4(Mt) + pow4(Mst2)*pow4(
        s2t))))) + 2*(lmMst1 - lmMst2)*pow2(Mst1)*pow2(s2t)*(-4*pow2(Mt)*pow2(
        MuSUSY)*(-1 + pow2(Sbeta))*(pow4(Mst1) + pow4(Mst2)) - pow2(Mst1)*pow2(
        Mst2)*(4*pow2(Mt)*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + pow2(s2t)*pow2(
        Sbeta)*pow4(Mst2)) + pow2(s2t)*pow2(Sbeta)*pow8(Mst2))))/(6.*pow2(Msq)*
        pow2(Sbeta)) + ((1 - 2*lmMst2)*shiftst3*((8*(-1 + lmMst1 - lmMst2)*
        pow2(Mst1)*pow2(Mt)*pow2(s2t) - 16*pow4(Mt) - (-4*pow2(Mst1)*pow2(Mst2)
        + (3 + 2*lmMst1 - 2*lmMst2)*pow4(Mst1) + pow4(Mst2))*pow4(s2t))*pow6(
        Mst2) - (4*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*(2*(1 - 2*lmMst1 + 2*lmMst2)
        *pow2(Mst2)*pow4(Mst1) + 2*(1 - lmMst1 + lmMst2)*pow2(Mst1)*pow4(Mst2)
        + (2 - 6*lmMst1 + 6*lmMst2)*pow6(Mst1) + pow6(Mst2)))/pow2(Sbeta) + 4*
        pow2(Mt)*pow2(s2t)*(pow2(MuSUSY)*pow6(Mst2) + 2*(pow2(MuSUSY)*((1 - 2*
        lmMst1 + 2*lmMst2)*pow2(Mst2)*pow4(Mst1) + (1 - lmMst1 + lmMst2)*pow2(
        Mst1)*pow4(Mst2) + (1 - 3*lmMst1 + 3*lmMst2)*pow6(Mst1)) + pow8(Mst2)))
        ))/12.)/(pow2(Mst1)*pow4(Mst2))) + (threeLoopFlag*(xDmsqst2*pow2(
        Dmsqst2)*pow2(Mst1)*pow2(Mst2)*(4*Mt*MuSUSY*pow2(Sbeta)*(37044000*pow2(
        Mst1)*(4*Dmglst2*(497 - 306*lmMst1 + 609*lmMst2 - 303*lmMt)*pow2(Mst1)
        + 36*(1 + 6*lmMst1 - 13*lmMst2 + 7*lmMt)*Mst2*pow2(Mst1) + 4*Dmglst2*(
        101 - 90*lmMst1 + 177*lmMst2 - 87*lmMt)*pow2(Mst2) + 9*(16 + 16*lmMst1
        - 35*lmMst2 + 19*lmMt)*pow3(Mst2))*pow3(Mt) + 500094000*Mt*pow2(Mst1)*
        pow2(s2t)*(-2*(1 + 8*lmMst1 - 8*lmMst2)*pow2(Mst1)*pow3(Mst2) + (4 - 9*
        lmMst1 + 9*lmMst2)*Mst2*pow4(Mst1) + Dmglst2*(40*(-1 + lmMst1 - lmMst2)
        *pow2(Mst1)*pow2(Mst2) + (-80 + 43*lmMst1 - 43*lmMst2)*pow4(Mst1) + 4*(
        2 + 5*lmMst1 - 5*lmMst2)*pow4(Mst2)) - 2*(5 + 4*lmMst1 - 4*lmMst2)*
        pow5(Mst2)) + 4*s2t*pow2(Mt)*(385875*(5904*Dmglst2 - Mst2*(3655 - 64*
        OepS2 + 8424*S2 + 144*lmMst1*(1 + 9*S2) - 144*lmMst2*(4 + 9*S2)))*pow2(
        Mst1)*pow3(Mst2) + 5145*Mst2*(16200*Dmglst2*(65 - 2*lmMst1 + 2*lmMst2)
        - Mst2*(279089 - 5600*OepS2 + 1260*lmMst2*(29 - 90*S2) + 812700*S2 +
        180*lmMst1*(-203 + 90*lmMst2 + 630*S2) - 8100*(pow2(lmMst1) + pow2(
        lmMst2))))*pow4(Mst1) - 83349000*(-12*Dmglst2 + Mst2 + 2*lmMst2*Mst2)*
        pow5(Mst2) - (1094369501 - 72716000*OepS2 + 2520*lmMst2*(235453 -
        584325*S2) + 5940931500*S2 + 2520*lmMst1*(-235453 + 114345*lmMst2 +
        584325*S2) - 144074700*(pow2(lmMst1) + pow2(lmMst2)))*pow6(Mst1)) +
        2058000*s2t*(-81*(2 - lmMst2)*shiftst1*(-4*pow2(Mst1)*pow2(Mt) + 4*
        pow2(Mst2)*pow2(Mt) + 2*(lmMst1 - lmMst2)*pow2(Mst1)*pow2(Mst2)*pow2(
        s2t) + pow2(s2t)*(pow4(Mst1) - pow4(Mst2)))*pow4(Mst2) - 2*T1ep*pow2(
        Mst1)*((106*pow2(Mt) - 57*pow2(Mst2)*pow2(s2t))*pow4(Mst1) + 36*pow2(
        Mt)*pow4(Mst2) + pow2(Mst1)*(42*pow2(Mst2)*pow2(Mt) - 57*pow2(s2t)*
        pow4(Mst2)) - 18*pow2(s2t)*pow6(Mst2))) + 3*Mst2*pow3(s2t)*(1715*(
        192439 - 30400*OepS2 - 837000*S2 - 180*lmMst2*(857 + 3420*S2) + 180*
        lmMst1*(677 + 180*lmMst2 + 3420*S2) - 16200*(pow2(lmMst1) + pow2(
        lmMst2)))*pow3(Mst2)*pow4(Mst1) + 128625*(8555 - 128*OepS2 - 57024*S2 -
        72*lmMst2*(29 + 36*S2) + 72*lmMst1*(29 - 12*lmMst2 + 36*S2) + 864*pow2(
        lmMst2))*pow2(Mst1)*pow5(Mst2) - 2*Mst2*(55313478 + 26068000*OepS2 +
        563377500*S2 + 315*lmMst2*(-423553 + 1675800*S2) - 315*lmMst1*(-423553
        + 166320*lmMst2 + 1675800*S2) + 26195400*(pow2(lmMst1) + pow2(lmMst2)))
        *pow6(Mst1) + 3087000*Dmglst2*(9*(55 - 34*lmMst1 + 34*lmMst2)*pow2(
        Mst2)*pow4(Mst1) - 36*(5 + 2*lmMst1 - 2*lmMst2)*pow2(Mst1)*pow4(Mst2) +
        (419 - 57*lmMst1 + 57*lmMst2)*pow6(Mst1) - 108*pow6(Mst2)) + 27783000*(
        1 + 2*lmMst2)*pow7(Mst2))) + Tbeta*(2667168000*s2t*pow2(Mst1)*pow2(
        MuSUSY)*(2*Dmglst2*(8 - 15*lmMst1 + 15*lmMst2)*pow2(Mst1) + (2 + 6*
        lmMst1 - 6*lmMst2)*Mst2*pow2(Mst1) - 2*Dmglst2*(2 + 5*lmMst1 - 5*
        lmMst2)*pow2(Mst2) + (5 + 4*lmMst1 - 4*lmMst2)*pow3(Mst2))*pow3(Mt) +
        20*Mt*pow2(s2t)*(-694575*Mst2*(-6*(11 + 39*lmMst1 - 39*lmMst2)*s2t*
        pow2(Mst2)*(-2 + pow2(Sbeta)) + (Dmglst2*(4324 - 660*lmMst1 + 660*
        lmMst2) + 3*(-167 + 5*lmMst1 - 5*lmMst2)*Mst2)*Mt*pow2(Sbeta))*pow4(
        Sbeta)*pow6(Mst1) + Mt*(1029*Mst2*pow2(MuSUSY)*(75*(1728*Dmglst2*(4 +
        lmMst1 - lmMst2) + Mst2*(-8771 + 128*OepS2 + 72*lmMst1*(-29 + 12*lmMst2
        - 36*S2) + 57024*S2 + 72*lmMst2*(23 + 36*S2) - 864*pow2(lmMst2)))*pow2(
        Mst1)*pow2(Mst2) - 4*(4050*Dmglst2*(23 - 42*lmMst1 + 42*lmMst2) + Mst2*
        (212566 - 69615*lmMst2 - 10000*OepS2 - 1350*(947 + 150*lmMst2)*S2 + 45*
        lmMst1*(1547 - 180*lmMst2 + 4500*S2) - 4050*pow2(lmMst1) + 12150*pow2(
        lmMst2)))*pow4(Mst1) - 16200*(-12*Dmglst2 + Mst2 + 2*lmMst2*Mst2)*pow4(
        Mst2)) - 2*((593331163 - 60642400*OepS2 + 1260*lmMst2*(8263 - 974610*
        S2) + 1143733500*S2 + 1260*lmMst1*(-8263 + 71820*lmMst2 + 974610*S2) -
        61916400*pow2(lmMst1) - 28576800*pow2(lmMst2))*pow2(MuSUSY) - 694575*
        Mst2*(Dmglst2*(4324 - 660*lmMst1 + 660*lmMst2) + 3*(-167 + 5*lmMst1 -
        5*lmMst2)*Mst2)*pow4(Sbeta))*pow6(Mst1))) + 4116000*pow2(s2t)*(-(T1ep*
        pow2(Mst1)*(3*pow2(s2t)*pow2(Sbeta)*pow4(Mst2)*(13*pow2(Mst1)*pow2(
        Mst2) + 6*pow4(Mst2)) + 4*pow2(Mt)*pow2(MuSUSY)*(75*pow2(Mst1)*pow2(
        Mst2) + 221*pow4(Mst1) + 18*pow4(Mst2)))) - 2*pow2(Mt)*(-2*T1ep*pow2(
        Mst1)*pow2(Sbeta)*((-7*pow2(Mst2) + 221*pow2(MuSUSY))*pow4(Mst1) + 18*(
        pow2(Mst2) + pow2(MuSUSY))*pow4(Mst2) + 3*pow2(Mst1)*(25*pow2(Mst2)*
        pow2(MuSUSY) + pow4(Mst2))) + 81*shiftst1*pow2(MuSUSY)*(2*lmMst1*(-2 +
        lmMst2)*pow2(Mst1)*pow2(Mst2)*(pow2(Mst1) + pow2(Mst2)) + 2*(pow2(Mst1)
        + pow2(Mst2))*pow4(Mst2) - 2*pow2(lmMst2)*pow2(Mst1)*(pow2(Mst1)*pow2(
        Mst2) + pow4(Mst1) + pow4(Mst2)) + lmMst1*(-3 + 2*lmMst2)*pow6(Mst1) +
        lmMst2*(4*pow2(Mst2)*pow4(Mst1) + 3*pow2(Mst1)*pow4(Mst2) + 3*pow6(
        Mst1) - pow6(Mst2))))) + 166698000*shiftst1*pow2(Sbeta)*(pow2(s2t)*((2
        - lmMst2)*pow2(Mst2)*pow4(Mst1)*(8*pow2(Mt)*(pow2(Mst2) + (-lmMst1 +
        lmMst2)*pow2(MuSUSY)) + (1 - 2*lmMst1 + 2*lmMst2)*pow2(s2t)*pow4(Mst2))
        + (4*(lmMst1 - lmMst2)*(-3 + 2*lmMst2)*pow2(Mt)*pow2(MuSUSY) + (-2 +
        lmMst2)*pow2(s2t)*pow4(Mst2))*pow6(Mst1)) - (2 - lmMst2)*(pow2(Mst1)*
        pow4(Mst2)*(4*pow2(Mt)*(4*pow2(Mst2) + (-1 + 2*lmMst1 - 2*lmMst2)*pow2(
        MuSUSY))*pow2(s2t) + 16*pow4(Mt) + (-1 - 2*lmMst1 + 2*lmMst2)*pow4(
        Mst2)*pow4(s2t)) + (-4*pow2(Mt)*(2*pow2(Mst2) + pow2(MuSUSY))*pow2(s2t)
        + 16*pow4(Mt) + pow4(Mst2)*pow4(s2t))*pow6(Mst2))) + pow2(Sbeta)*(-
        74088000*s2t*pow2(Mst1)*pow3(Mt)*(9*Mst2*(4*(5 + 4*lmMst1 - 4*lmMst2)*
        pow2(Mst2)*pow2(MuSUSY) + 4*pow2(Mst1)*((-23 + 14*lmMst1 - 28*lmMst2 +
        14*lmMt)*pow2(Mst2) + 2*(1 + 3*lmMst1 - 3*lmMst2)*pow2(MuSUSY)) + 2*(-
        35 - 50*lmMst2 + 4*lmMst1*(8 + lmMst2 - lmMt) + 18*lmMt + 4*lmMst2*lmMt
        - 4*pow2(lmMst2))*pow4(Mst1) + (13 + 32*lmMst1 - 70*lmMst2 + 38*lmMt)*
        pow4(Mst2)) - 2*Dmglst2*(36*(2 + 5*lmMst1 - 5*lmMst2)*pow2(Mst2)*pow2(
        MuSUSY) + 36*pow2(Mst1)*(8*(-7 + 3*lmMst1 - 6*lmMst2 + 3*lmMt)*pow2(
        Mst2) + (-8 + 15*lmMst1 - 15*lmMst2)*pow2(MuSUSY)) + 9*(-465 - 374*
        lmMst2 + 4*lmMst1*(52 + 3*lmMst2 - 3*lmMt) + 166*lmMt + 12*lmMst2*lmMt
        - 12*pow2(lmMst2))*pow4(Mst1) + (-578 + 360*lmMst1 - 708*lmMst2 + 348*
        lmMt)*pow4(Mst2))) - 166698000*Mt*pow2(Mst1)*pow2(Mst2)*pow3(s2t)*(32*(
        1 - lmMst1 + lmMst2)*pow2(Mst1)*pow3(Mst2) + (45 + 37*lmMst1 - 37*
        lmMst2)*Mst2*pow4(Mst1) - 4*Dmglst2*(4*(12 - 5*lmMst1 + 5*lmMst2)*pow2(
        Mst1)*pow2(Mst2) + (40 - 3*lmMst1 + 3*lmMst2)*pow4(Mst1) - 4*(2 + 5*
        lmMst1 - 5*lmMst2)*pow4(Mst2)) - 8*(5 + 4*lmMst1 - 4*lmMst2)*pow5(Mst2)
        ) - 144*pow4(Mt)*(-42875*(8*Dmglst2*(2954 - 864*lmMst1 + 1713*lmMst2 -
        849*lmMt) + Mst2*(-5567 + 2232*lmMst1 + 2838*lmMt - 6*lmMst2*(917 + 36*
        lmMt) + 108*pow2(lmMst2) + 108*pow2(lmMt)))*pow2(Mst1)*pow3(Mst2) -
        30870*Mst2*(300*Dmglst2*(401 + 274*lmMst2 - 2*lmMst1*(71 + 2*lmMst2 -
        2*lmMt) - 4*(33 + lmMst2)*lmMt + 4*pow2(lmMst2)) - Mst2*(23941 + 20010*
        lmMst2 - 30*lmMst1*(347 + 10*lmMst2 - 10*lmMt) - 300*(32 + lmMst2)*lmMt
        + 300*pow2(lmMst2)))*pow4(Mst1) + 9261000*(-12*Dmglst2 + Mst2 + 2*
        lmMst2*Mst2)*pow5(Mst2) + 24*(26331136 - 105*lmMst1*(236317 + 35070*
        lmMst2 - 36750*lmMt) + 210*lmMst2*(175121 - 18375*lmMt) - 11962125*lmMt
        - 88200*pow2(lmMst1) + 3770550*pow2(lmMst2))*pow6(Mst1)) + 4*pow2(Mt)*
        pow2(s2t)*(20580*pow2(Mst2)*((2482 - 400*OepS2 + 90*lmMst2*(443 - 90*
        S2) + 90450*S2 - 90*lmMst1*(263 - 90*(lmMst2 + S2)) - 4050*(pow2(
        lmMst1) + pow2(lmMst2)))*pow2(Mst2) - (45*lmMst1*(-1547 + 180*lmMst2 -
        4500*S2) + 45*lmMst2*(1547 + 4500*S2) + 2*(-106283 + 5000*OepS2 +
        639225*S2) + 4050*pow2(lmMst1) - 12150*pow2(lmMst2))*pow2(MuSUSY))*
        pow4(Mst1) + 385875*pow2(Mst1)*(2*(3439 - 64*OepS2 + 8424*S2 + 144*
        lmMst1*(1 + 9*S2) - 144*lmMst2*(7 + 9*S2))*pow2(Mst2) + (8771 - 128*
        OepS2 - 57024*S2 - 72*lmMst2*(23 + 36*S2) + 72*lmMst1*(29 - 12*lmMst2 +
        36*S2) + 864*pow2(lmMst2))*pow2(MuSUSY))*pow4(Mst2) + 2*((194011877 +
        9604000*OepS2 - 319504500*S2 + 1260*lmMst2*(60437 + 154350*S2) - 1260*
        lmMst1*(60437 + 15120*lmMst2 + 154350*S2) + 9525600*(pow2(lmMst1) +
        pow2(lmMst2)))*pow2(Mst2) + 5*(593331163 - 60642400*OepS2 + 1260*
        lmMst2*(8263 - 974610*S2) + 1143733500*S2 + 1260*lmMst1*(-8263 + 71820*
        lmMst2 + 974610*S2) - 61916400*pow2(lmMst1) - 28576800*pow2(lmMst2))*
        pow2(MuSUSY))*pow6(Mst1) - 9261000*Dmglst2*Mst2*(((678 - 36*lmMst1 +
        36*lmMst2)*pow2(Mst2) + 9*(-23 + 42*lmMst1 - 42*lmMst2)*pow2(MuSUSY))*
        pow4(Mst1) + 108*(2*pow2(Mst2) + pow2(MuSUSY))*pow4(Mst2) + 12*pow2(
        Mst1)*(6*(4 + lmMst1 - lmMst2)*pow2(Mst2)*pow2(MuSUSY) + 23*pow4(Mst2))
        + (656 - 192*lmMst1 + 192*lmMst2)*pow6(Mst1)) + 83349000*(1 + 2*lmMst2)
        *(2*pow2(Mst2) + pow2(MuSUSY))*pow6(Mst2)) - 3*pow3(Mst2)*pow4(s2t)*(-
        3430*(224593 + 10400*OepS2 - 1719900*S2 + 1170*lmMst2*(-1 + 180*S2) -
        90*lmMst1*(-193 + 540*lmMst2 + 2340*S2) + 8100*pow2(lmMst1) + 40500*
        pow2(lmMst2))*pow3(Mst2)*pow4(Mst1) + 128625*(8339 - 128*OepS2 - 57024*
        S2 - 72*lmMst2*(35 + 36*S2) + 72*lmMst1*(29 - 12*lmMst2 + 36*S2) + 864*
        pow2(lmMst2))*pow2(Mst1)*pow5(Mst2) - Mst2*(440659841 + 1890*lmMst1*(
        251761 - 26040*lmMst2) - 531394290*lmMst2 - 308700000*S2 + 24607800*
        pow2(lmMst1) + 24607800*pow2(lmMst2))*pow6(Mst1) - 3087000*Dmglst2*(9*(
        -75 + 26*lmMst1 - 26*lmMst2)*pow2(Mst2)*pow4(Mst1) + 72*(1 + lmMst1 -
        lmMst2)*pow2(Mst1)*pow4(Mst2) + (76 - 249*lmMst1 + 249*lmMst2)*pow6(
        Mst1) + 108*pow6(Mst2)) + 27783000*(1 + 2*lmMst2)*pow7(Mst2))))) +
        308700*xDR2DRMOD*pow2(Msq)*(4*MuSUSY*pow2(Sbeta)*(-1215*shiftst1*(
        Dmsqst2 + pow2(Msq))*pow2(Mst1)*(4*s2t*(pow2(Mst1) - pow2(Mst2))*pow3(
        Mt) - Mt*pow3(s2t)*(2*(lmMst1 - lmMst2)*pow2(Mst1)*pow2(Mst2) + pow4(
        Mst1) - pow4(Mst2)))*pow6(Mst2) + 128*pow2(Msq)*pow2(Mst1)*pow4(Mt)*(9*
        (1 + lmMst2)*Mst2*((3 - 21*lmMst2 + 2*lmMst1*(8 + 3*lmMst2 - 3*lmMt) +
        5*lmMt + 6*lmMst2*lmMt - 6*pow2(lmMst2))*pow2(Mst2)*pow4(Mst1) + (-1 -
        7*lmMst2 + 2*lmMst1*(2 + lmMst2 - lmMt) + 3*lmMt + 2*lmMst2*lmMt - 2*
        pow2(lmMst2))*pow2(Mst1)*pow4(Mst2) + (12 - 45*lmMst2 + 2*lmMst1*(19 +
        6*lmMst2 - 6*lmMt) + 7*lmMt + 12*lmMst2*lmMt - 12*pow2(lmMst2))*pow6(
        Mst1) - 2*(1 + lmMst2)*pow6(Mst2)) + Dmglst2*(pow2(Mst2)*(263 + lmMst2*
        (962 - 213*lmMt) - 177*lmMt + (393 - 18*lmMt)*pow2(lmMst2) - 18*lmMst1*
        (26 - lmMst2*(-17 + lmMt) - 7*lmMt + pow2(lmMst2)) + 18*pow3(lmMst2))*
        pow4(Mst1) - pow2(Mst1)*(17*(-7 + 3*lmMt) + lmMst2*(-224 + 15*lmMt) -
        3*(5 + 6*lmMt)*pow2(lmMst2) - 18*lmMst1*(-4 + lmMt - lmMst2*(1 + lmMt)
        + pow2(lmMst2)) + 18*pow3(lmMst2))*pow4(Mst2) + (290 + lmMst2*(2447 -
        645*lmMt) - 375*lmMt - 3*(-509 + 60*lmMt)*pow2(lmMst2) - 18*lmMst1*(87
        + lmMst2*(71 - 10*lmMt) - 22*lmMt + 10*pow2(lmMst2)) + 180*pow3(lmMst2)
        )*pow6(Mst1) - 18*(-1 + 2*lmMst2 + 3*pow2(lmMst2))*pow6(Mst2))) - 864*
        pow2(Msq)*pow2(Mst1)*pow2(Mst2)*pow2(Mt)*pow2(s2t)*((1 + lmMst2)*Mst2*(
        2*(2 - 5*lmMst2 + lmMst1*(5 + 2*lmMst2) - 2*pow2(lmMst2))*pow2(Mst2)*
        pow4(Mst1) + 2*(-2*lmMst2*(2 + lmMst2) + lmMst1*(3 + 2*lmMst2))*pow2(
        Mst1)*pow4(Mst2) + (5 - 12*lmMst2 + 4*lmMst1*(3 + lmMst2) - 4*pow2(
        lmMst2))*pow6(Mst1) - 2*(1 + lmMst2)*pow6(Mst2)) + Dmglst2*(2*pow2(
        Mst2)*(8 + 19*lmMst2 - 5*pow2(lmMst2) + lmMst1*(-7 + 5*lmMst2 + 6*pow2(
        lmMst2)) - 6*pow3(lmMst2))*pow4(Mst1) + 2*pow2(Mst1)*(7 + 9*lmMst2 - 8*
        pow2(lmMst2) + lmMst1*(-3 + 5*lmMst2 + 6*pow2(lmMst2)) - 6*pow3(lmMst2)
        )*pow4(Mst2) + (17 + 51*lmMst2 - 4*pow2(lmMst2) + 4*lmMst1*(-6 + lmMst2
        + 3*pow2(lmMst2)) - 12*pow3(lmMst2))*pow6(Mst1) - 2*(-1 + 2*lmMst2 + 3*
        pow2(lmMst2))*pow6(Mst2))) - 9*Mt*pow3(Mst2)*pow3(s2t)*(75*Dmsqst2*
        pow2(Mst1)*pow3(Mst2)*(2*(lmMst1 - lmMst2)*pow2(Mst1)*pow2(Mst2) +
        pow4(Mst1) - pow4(Mst2)) - 2*Dmglst2*(pow2(Mst1)*pow2(Mst2)*(4*pow2(
        Msq)*(20*pow2(Mst1)*pow2(Mst2) + 9*pow4(Mst1) - 5*pow4(Mst2)) + 15*
        Dmsqst2*(pow4(Mst1) - pow4(Mst2))) - 30*Dmsqst2*lmMst2*pow4(Mst1)*pow4(
        Mst2) + pow4(Mst1)*(30*Dmsqst2*lmMst1*pow4(Mst2) + pow2(Msq)*(32*
        lmMst2*pow2(lmMst1)*(8*pow2(Mst1)*pow2(Mst2) + 8*pow4(Mst1) + pow4(
        Mst2)) + 2*pow3(lmMst2)*(128*pow2(Mst1)*(pow2(Mst1) + pow2(Mst2)) +
        107*pow4(Mst2)))) + pow2(Msq)*(-2*lmMst1*((-4 + 246*lmMst2 + 123*pow2(
        lmMst2))*pow4(Mst1)*pow4(Mst2) + 16*(3 + 3*lmMst2 + 16*pow2(lmMst2))*
        pow2(Mst2)*pow6(Mst1) - 16*lmMst2*pow2(Mst1)*pow6(Mst2) + 16*(3 + 2*
        lmMst2 + 16*pow2(lmMst2))*pow8(Mst1)) + pow2(lmMst2)*(396*pow4(Mst1)*
        pow4(Mst2) + 37*pow2(Mst2)*pow6(Mst1) + 155*pow2(Mst1)*pow6(Mst2) + 64*
        pow8(Mst1) - 32*pow8(Mst2))) + 2*lmMst2*pow2(Msq)*(4*pow4(Mst1)*pow4(
        Mst2) + 21*pow2(Mst2)*pow6(Mst1) + 115*pow2(Mst1)*pow6(Mst2) + 56*pow8(
        Mst1) - 16*pow8(Mst2))) - Mst2*pow2(Msq)*(2*(-8 + 237*lmMst2 + 16*(1 +
        lmMst2)*pow2(lmMst1) + 308*pow2(lmMst2) - lmMst1*(269 + 348*lmMst2 +
        123*pow2(lmMst2)) + 107*pow3(lmMst2))*pow4(Mst1)*pow4(Mst2) + pow2(
        Mst2)*(-125 - 156*lmMst2 - 512*lmMst1*lmMst2*(1 + lmMst2) + 256*(1 +
        lmMst2)*pow2(lmMst1) + 181*pow2(lmMst2) + 256*pow3(lmMst2))*pow6(Mst1)
        + (221 + 284*lmMst2 + 32*lmMst1*(1 + lmMst2) + 107*pow2(lmMst2))*pow2(
        Mst1)*pow6(Mst2) + 16*(1 + lmMst2)*(1 + lmMst1*(2 - 32*lmMst2) - 2*
        lmMst2 + 16*pow2(lmMst1) + 16*pow2(lmMst2))*pow8(Mst1) - 16*pow2(1 +
        lmMst2)*pow8(Mst2))) + 12*Mst2*s2t*pow3(Mt)*(-90*Dmglst2*Dmsqst2*pow2(
        Mst1)*(pow2(Mst1) - pow2(Mst2))*pow4(Mst2) + 225*Dmsqst2*pow2(Mst1)*(
        pow2(Mst1) - pow2(Mst2))*pow5(Mst2) + 2*Dmglst2*pow2(Msq)*((52 + 370*
        lmMst2 - 96*lmMst2*pow2(lmMst1) + 81*pow2(lmMst2) + 96*lmMst1*(-1 +
        lmMst2 + 2*pow2(lmMst2)) - 96*pow3(lmMst2))*pow4(Mst1)*pow4(Mst2) - 96*
        (lmMst2*pow2(lmMst1) + lmMst1*(4 + 2*lmMst2 - 2*pow2(lmMst2)) + (-4 +
        lmMst2)*pow2(1 + lmMst2))*pow2(Mst2)*pow6(Mst1) - 3*(-52 + 2*(51 + 16*
        lmMst1)*lmMst2 + 59*pow2(lmMst2))*pow2(Mst1)*pow6(Mst2) - 96*(-6 - 14*
        lmMst2 + lmMst2*pow2(lmMst1) + lmMst1*(9 + 6*lmMst2 - 2*pow2(lmMst2)) -
        6*pow2(lmMst2) + pow3(lmMst2))*pow8(Mst1) + 96*lmMst2*(1 + lmMst2)*
        pow8(Mst2)) - 3*Mst2*pow2(Msq)*(-((109 + 76*lmMst2 - 21*pow2(lmMst2) +
        64*lmMst1*pow2(1 + lmMst2) - 32*((1 + lmMst2)*pow2(lmMst1) + pow3(
        lmMst2)))*pow4(Mst1)*pow4(Mst2)) + (173 + 188*lmMst2 + 32*lmMst1*(1 +
        lmMst2) + 59*pow2(lmMst2))*pow2(Mst1)*pow6(Mst2) + 32*(1 + lmMst2)*(
        pow2(1 - lmMst1 + lmMst2)*pow2(Mst2)*pow6(Mst1) + (1 + 3*lmMst2 -
        lmMst1*(3 + 2*lmMst2) + pow2(lmMst1) + pow2(lmMst2))*pow8(Mst1)) - 16*
        pow2(1 + lmMst2)*pow8(Mst2)))) + Tbeta*(pow2(Mst1)*(-9*pow6(Mst2)*(-30*
        Dmglst2*Dmsqst2*pow2(Sbeta)*pow4(s2t)*pow5(Mst2) + 2*Dmglst2*s2t*pow2(
        Msq)*(256*(-1 + 2*lmMst2 + 3*pow2(lmMst2))*pow2(MuSUSY)*pow3(Mt) + (-20
        + (262 + 32*lmMst1)*lmMst2 + 187*pow2(lmMst2))*pow2(Sbeta)*pow3(s2t)*
        pow5(Mst2)) + (75*Dmsqst2 + (237 + 316*lmMst2 + 32*lmMst1*(1 + lmMst2)
        + 123*pow2(lmMst2))*pow2(Msq))*pow2(Sbeta)*pow4(s2t)*pow6(Mst2)) +
        1215*shiftst1*(Dmsqst2 + pow2(Msq))*pow2(Mst2)*(4*pow2(Mt)*pow2(MuSUSY)
        *pow2(s2t)*(pow2(Mst1)*pow4(Mst2) - 2*(lmMst1 - lmMst2)*pow2(Mst1)*(
        pow2(Mst1)*pow2(Mst2) + pow4(Mst1) + pow4(Mst2)) + pow6(Mst2)) + pow2(
        Sbeta)*(16*(pow2(Mst1) + pow2(Mst2))*pow4(Mst2)*pow4(Mt) + (pow2(Mst1)
        - pow2(Mst2))*(2*(lmMst1 - lmMst2)*pow2(Mst1)*pow2(Mst2) + pow4(Mst1) -
        pow4(Mst2))*pow4(Mst2)*pow4(s2t) - 4*pow2(Mt)*pow2(s2t)*(2*pow4(Mst1)*
        pow4(Mst2) - 4*pow2(Mst1)*pow6(Mst2) + pow2(MuSUSY)*(pow2(Mst1)*pow4(
        Mst2) - 2*(lmMst1 - lmMst2)*pow2(Mst1)*(pow2(Mst1)*pow2(Mst2) + pow4(
        Mst1) + pow4(Mst2)) + pow6(Mst2)) + 2*pow8(Mst2))))) + pow2(Sbeta)*(9*
        pow4(Mst1)*(Dmsqst2*(30*Dmglst2 - 75*Mst2)*((-1 + 2*lmMst1 - 2*lmMst2)*
        pow2(Mst1)*pow2(Mst2) + pow4(Mst1) + (-1 - 2*lmMst1 + 2*lmMst2)*pow4(
        Mst2)) + pow2(Msq)*(-2*Dmglst2*(pow2(Mst1)*pow2(Mst2)*(-44 + 34*lmMst2
        + 224*lmMst2*pow2(lmMst1) - 359*pow2(lmMst2) - 2*lmMst1*(52 - 198*
        lmMst2 + 133*pow2(lmMst2)) + 42*pow3(lmMst2)) + (-36 + (70 + 32*lmMst1)
        *lmMst2 + 27*pow2(lmMst2))*pow4(Mst1) + (100 - 222*lmMst2 + 32*lmMst2*
        pow2(lmMst1) + lmMst1*(8 - 524*lmMst2 - 246*pow2(lmMst2)) + 241*pow2(
        lmMst2) + 214*pow3(lmMst2))*pow4(Mst2)) - Mst2*(pow2(Mst1)*pow2(Mst2)*(
        -109 - 630*lmMst2 + 224*(1 + lmMst2)*pow2(lmMst1) + lmMst1*(538 + 184*
        lmMst2 - 266*pow2(lmMst2)) - 435*pow2(lmMst2) + 42*pow3(lmMst2)) + (141
        + 140*lmMst2 + 32*lmMst1*(1 + lmMst2) + 43*pow2(lmMst2))*pow4(Mst1) + (
        -237 + 190*lmMst2 + 32*(1 + lmMst2)*pow2(lmMst1) + 509*pow2(lmMst2) -
        2*lmMst1*(285 + 364*lmMst2 + 123*pow2(lmMst2)) + 214*pow3(lmMst2))*
        pow4(Mst2))))*pow4(s2t)*pow5(Mst2) + 1152*Mt*pow2(Msq)*pow2(Mst1)*pow3(
        s2t)*pow4(Mst2)*((1 + lmMst2)*Mst2*(2*(2 + 2*lmMst1 - lmMst2)*pow2(
        Mst2)*pow4(Mst1) + 2*(1 - 3*lmMst2 + lmMst1*(3 + 2*lmMst2) - 2*pow2(
        lmMst2))*pow2(Mst1)*pow4(Mst2) + (1 + 2*lmMst1 - 2*lmMst2)*pow6(Mst1) -
        2*(1 + lmMst2)*pow6(Mst2)) + Dmglst2*(2*(1 - 4*lmMst1 + 10*lmMst2 + 3*
        pow2(lmMst2))*pow2(Mst2)*pow4(Mst1) + 2*pow2(Mst1)*(6 + 11*lmMst2 - 5*
        pow2(lmMst2) + lmMst1*(-3 + 5*lmMst2 + 6*pow2(lmMst2)) - 6*pow3(lmMst2)
        )*pow4(Mst2) + (1 + 13*lmMst2 - 2*lmMst1*(5 + 3*lmMst2) + 6*pow2(
        lmMst2))*pow6(Mst1) - 2*(-1 + 2*lmMst2 + 3*pow2(lmMst2))*pow6(Mst2))) +
        256*s2t*pow2(Msq)*pow2(Mst1)*pow3(Mt)*(9*(1 + lmMst2)*Mst2*(-2*pow2(
        Mst2)*((3 + 2*lmMst1*(7 + 2*lmMst2 - 2*lmMt) + 4*lmMst2*(-4 + lmMt) +
        2*lmMt - 4*pow2(lmMst2))*pow2(Mst2) + (1 - 10*lmMst2 + 4*lmMst1*(2 +
        lmMst2) - 4*pow2(lmMst2))*pow2(MuSUSY))*pow4(Mst1) + pow2(Mst1)*((1 -
        2*lmMst1*(5 + 2*lmMst2 - 2*lmMt) - 4*lmMst2*(-3 + lmMt) - 6*lmMt + 4*
        pow2(lmMst2))*pow2(Mst2) + 2*(1 + 5*lmMst2 - lmMst1*(3 + 2*lmMst2) + 2*
        pow2(lmMst2))*pow2(MuSUSY))*pow4(Mst2) - (2*(8 - 27*lmMst2 + lmMst1*(25
        + 6*lmMst2 - 6*lmMt) + 2*lmMt + 6*lmMst2*lmMt - 6*pow2(lmMst2))*pow2(
        Mst2) + (7 - 32*lmMst2 + 4*lmMst1*(7 + 3*lmMst2) - 12*pow2(lmMst2))*
        pow2(MuSUSY))*pow6(Mst1) + 2*(1 + lmMst2)*(2*pow2(Mst2) + pow2(MuSUSY))
        *pow6(Mst2)) - Dmglst2*(pow2(Mst1)*(pow2(Mst2)*(253 + 5*lmMst2*(107 -
        6*lmMt) - 102*lmMt - 18*lmMst1*(9 + lmMst2 - 2*lmMt + 2*lmMst2*lmMt -
        2*pow2(lmMst2)) + 12*(10 + 3*lmMt)*pow2(lmMst2) - 36*pow3(lmMst2)) +
        18*pow2(MuSUSY)*(8 + 7*lmMst2 - 11*pow2(lmMst2) + lmMst1*(-3 + 5*lmMst2
        + 6*pow2(lmMst2)) - 6*pow3(lmMst2)))*pow4(Mst2) + 9*((49 + 103*lmMst2 -
        36*(1 + lmMst2)*pow2(lmMst2) + 4*lmMst1*(-11 + 6*lmMst2 + 9*pow2(
        lmMst2)))*pow2(MuSUSY) + pow2(Mst2)*(28 + lmMst2*(378 - 96*lmMt) - 44*
        lmMt + lmMst1*(-274 + 60*lmMt + 18*lmMst2*(-13 + 2*lmMt) - 36*pow2(
        lmMst2)) + (270 - 36*lmMt)*pow2(lmMst2) + 36*pow3(lmMst2)))*pow6(Mst1)
        + 18*(pow2(Mst2)*(2*pow2(MuSUSY)*(8 + 13*lmMst2 - 8*pow2(lmMst2) +
        lmMst1*(-5 + 5*lmMst2 + 6*pow2(lmMst2)) - 6*pow3(lmMst2)) + pow2(Mst2)*
        (23 + lmMst2*(93 - 22*lmMt) - 14*lmMt - 4*(-11 + lmMt)*pow2(lmMst2) -
        2*lmMst1*(25 + lmMst2*(17 - 2*lmMt) - 6*lmMt + 2*pow2(lmMst2)) + 4*
        pow3(lmMst2)))*pow4(Mst1) + (1 - 2*lmMst2 - 3*pow2(lmMst2))*(2*pow2(
        Mst2) + pow2(MuSUSY))*pow6(Mst2)))) - 16*Mst2*pow4(Mt)*(675*Dmsqst2*
        pow2(Mst1)*(pow2(Mst1) + pow2(Mst2))*pow5(Mst2) - Mst2*pow2(Msq)*(-((
        2117 + lmMst2*(4796 - 1248*lmMt) - 576*lmMst1*(1 + lmMst2)*(3 + lmMst2
        - lmMt) - 672*lmMt + (3651 - 576*lmMt)*pow2(lmMst2) + 576*pow3(lmMst2))
        *pow4(Mst1)*pow4(Mst2)) + 288*(1 + lmMst2)*((5 - 57*lmMst2 + 2*lmMst1*(
        25 + 6*lmMst2 - 6*lmMt) + 7*lmMt + 12*lmMst2*lmMt - 12*pow2(lmMst2))*
        pow2(Mst1) + (-2 - 27*lmMst2 + lmMst1*(22 + 6*lmMst2 - 6*lmMt) + 5*lmMt
        + 6*lmMst2*lmMt - 6*pow2(lmMst2))*pow2(Mst2))*pow6(Mst1) - 9*(173 +
        188*lmMst2 + 32*lmMst1*(1 + lmMst2) + 59*pow2(lmMst2))*pow2(Mst1)*pow6(
        Mst2) + 144*pow2(1 + lmMst2)*pow8(Mst2)) - 2*Dmglst2*(135*Dmsqst2*pow2(
        Mst1)*(pow2(Mst1) + pow2(Mst2))*pow4(Mst2) + pow2(Msq)*((4220 + lmMst2*
        (6114 - 1296*lmMt) - 576*lmMst1*(4 + 2*lmMst2 - lmMt) - 1392*lmMt +
        1341*pow2(lmMst2))*pow4(Mst1)*pow4(Mst2) - 9*(-84 + (70 + 32*lmMst1)*
        lmMst2 + 59*pow2(lmMst2))*pow2(Mst1)*pow6(Mst2) + 288*(pow2(Mst2)*(31 +
        lmMst2*(95 - 23*lmMt) - 16*lmMt + (51 - 6*lmMt)*pow2(lmMst2) - 2*
        lmMst1*(25 + lmMst2*(20 - 3*lmMt) - 6*lmMt + 3*pow2(lmMst2)) + 6*pow3(
        lmMst2))*pow6(Mst1) + (42 + lmMst2*(242 - 62*lmMt) - 33*lmMt + 6*(29 -
        4*lmMt)*pow2(lmMst2) - 2*lmMst1*(81 + lmMst2*(74 - 12*lmMt) - 18*lmMt +
        12*pow2(lmMst2)) + 24*pow3(lmMst2))*pow8(Mst1)) + 288*lmMst2*(1 +
        lmMst2)*pow8(Mst2)))) - 12*Mst2*pow2(Mt)*pow2(s2t)*(2*Dmglst2*(45*
        Dmsqst2*pow2(Mst1)*pow2(Mst2)*(2*(pow2(Mst2) + (-lmMst1 + lmMst2)*pow2(
        MuSUSY))*pow4(Mst1) + pow2(Mst1)*((1 - 2*lmMst1 + 2*lmMst2)*pow2(Mst2)*
        pow2(MuSUSY) - 4*pow4(Mst2)) + (2*pow2(Mst2) + pow2(MuSUSY))*pow4(Mst2)
        ) + pow2(Msq)*(2*(-3*pow2(Mst2)*pow2(MuSUSY)*(48 + 4*lmMst2*(31 + 36*
        pow2(lmMst1)) + 278*pow2(lmMst2) - lmMst1*(44 + 278*lmMst2 + 379*pow2(
        lmMst2)) + 235*pow3(lmMst2)) + (332 + 302*lmMst2 - 288*lmMst1*(1 +
        lmMst2) + 111*pow2(lmMst2))*pow4(Mst2))*pow6(Mst1) - pow4(Mst1)*(3*
        pow2(MuSUSY)*(60 + 206*lmMst2 + 32*lmMst2*pow2(lmMst1) + lmMst1*(8 -
        460*lmMst2 - 246*pow2(lmMst2)) + 519*pow2(lmMst2) + 214*pow3(lmMst2))*
        pow4(Mst2) + 4*(52 + lmMst2*(-338 + 48*pow2(lmMst1)) - 129*pow2(lmMst2)
        + 48*(lmMst1 - 2*lmMst1*lmMst2*(1 + lmMst2) + pow3(lmMst2)))*pow6(Mst2)
        ) + 6*(32*(2 + 7*lmMst2 - lmMst1*(5 + 4*lmMst2) + 4*pow2(lmMst2))*pow2(
        Mst2) - pow2(MuSUSY)*(48 + 4*lmMst2*(45 + 68*pow2(lmMst1)) + 310*pow2(
        lmMst2) - lmMst1*(92 + 310*lmMst2 + 635*pow2(lmMst2)) + 363*pow3(
        lmMst2)))*pow8(Mst1) + 96*lmMst2*(1 + lmMst2)*(2*pow2(Mst2) + pow2(
        MuSUSY))*pow8(Mst2) + 3*pow2(Mst1)*(-((-20 + 2*(99 + 16*lmMst1)*lmMst2
        + 123*pow2(lmMst2))*pow2(MuSUSY)*pow6(Mst2)) - 2*(-52 + 2*(67 + 16*
        lmMst1)*lmMst2 + 91*pow2(lmMst2))*pow8(Mst2)))) - 3*Mst2*(6*(25*Dmsqst2
        + 47*pow2(Msq))*pow4(Mst2)*pow6(Mst1) - 12*(25*Dmsqst2 + 47*pow2(Msq))*
        pow4(Mst1)*pow6(Mst2) + pow2(MuSUSY)*(3*(25*Dmsqst2 + 63*pow2(Msq))*
        pow4(Mst1)*pow4(Mst2) + 5*(15*Dmsqst2 + 41*pow2(Msq))*pow2(Mst1)*pow6(
        Mst2)) - 2*lmMst1*(75*Dmsqst2*pow2(MuSUSY)*pow4(Mst1)*(pow2(Mst1)*pow2(
        Mst2) + pow4(Mst1) + pow4(Mst2)) + pow2(Msq)*((32*(3 + 5*lmMst2 + 2*
        pow2(lmMst2))*pow2(Mst2) + (253 + 332*lmMst2 + 123*pow2(lmMst2))*pow2(
        MuSUSY))*pow4(Mst1)*pow4(Mst2) + (253 + 588*lmMst2 + 379*pow2(lmMst2))*
        pow2(Mst2)*pow2(MuSUSY)*pow6(Mst1) - 16*(1 + lmMst2)*pow2(Mst1)*(2*
        pow2(Mst2) + pow2(MuSUSY))*pow6(Mst2) + (32*(1 + lmMst2)*pow2(Mst2) + (
        237 + 828*lmMst2 + 635*pow2(lmMst2))*pow2(MuSUSY))*pow8(Mst1))) + 6*(
        25*Dmsqst2 + 63*pow2(Msq))*pow2(Mst1)*pow8(Mst2) + 2*lmMst2*(75*
        Dmsqst2*pow2(MuSUSY)*pow4(Mst1)*(pow2(Mst1)*pow2(Mst2) + pow4(Mst1) +
        pow4(Mst2)) + pow2(Msq)*((285*pow2(Mst2)*pow2(MuSUSY) + 172*pow4(Mst2))
        *pow6(Mst1) - 33*pow4(Mst1)*(-11*pow2(MuSUSY)*pow4(Mst2) + 8*pow6(Mst2)
        ) + (32*pow2(Mst2) + 277*pow2(MuSUSY))*pow8(Mst1) - 16*(2*pow2(Mst2) +
        pow2(MuSUSY))*pow8(Mst2) + 2*pow2(Mst1)*(63*pow2(MuSUSY)*pow6(Mst2) +
        110*pow8(Mst2)))) + pow2(Msq)*(pow4(Mst1)*(32*(1 + lmMst2)*pow2(lmMst1)
        *(pow2(MuSUSY)*(9*pow2(Mst1)*pow2(Mst2) + 17*pow4(Mst1) + pow4(Mst2)) +
        2*pow6(Mst2)) + 2*pow3(lmMst2)*(pow2(MuSUSY)*(235*pow2(Mst1)*pow2(Mst2)
        + 363*pow4(Mst1) + 107*pow4(Mst2)) + 32*pow6(Mst2))) - 16*pow2(MuSUSY)*
        pow8(Mst2) - pow2(lmMst2)*(-6*(148*pow2(Mst2)*pow2(MuSUSY) + 25*pow4(
        Mst2))*pow6(Mst1) + pow4(Mst1)*(-707*pow2(MuSUSY)*pow4(Mst2) + 76*pow6(
        Mst2)) - 8*(8*pow2(Mst2) + 139*pow2(MuSUSY))*pow8(Mst1) + 16*(2*pow2(
        Mst2) + pow2(MuSUSY))*pow8(Mst2) - pow2(Mst1)*(91*pow2(MuSUSY)*pow6(
        Mst2) + 150*pow8(Mst2)))) + pow2(Msq)*(16*pow2(MuSUSY)*(4*pow2(Mst2)*
        pow6(Mst1) + 5*pow8(Mst1)) - 32*power10(Mst2))))) - 36*s2t*pow2(Mt)*
        pow2(MuSUSY)*(15*Dmsqst2*s2t*pow2(Mst1)*pow2(Mst2)*(-2*Dmglst2*Mst2*(
        pow2(Mst1) + pow2(Mst2))*(-2*lmMst1*pow2(Mst1) + 2*lmMst2*pow2(Mst1) +
        pow2(Mst2)) + 5*(pow2(Mst1) + pow2(Mst2))*pow4(Mst2) - 10*(lmMst1 -
        lmMst2)*pow2(Mst1)*(pow2(Mst1)*pow2(Mst2) + pow4(Mst1) + pow4(Mst2))) +
        pow2(Msq)*(-2*Dmglst2*(64*Mt*(8 + 7*lmMst2 - 11*pow2(lmMst2) + lmMst1*(
        -3 + 5*lmMst2 + 6*pow2(lmMst2)) - 6*pow3(lmMst2)) - Mst2*s2t*(60 + 206*
        lmMst2 + 32*lmMst2*pow2(lmMst1) + lmMst1*(8 - 460*lmMst2 - 246*pow2(
        lmMst2)) + 519*pow2(lmMst2) + 214*pow3(lmMst2)))*pow4(Mst1)*pow4(Mst2)
        + (128*(1 + lmMst2)*Mt*(1 + 5*lmMst2 - lmMst1*(3 + 2*lmMst2) + 2*pow2(
        lmMst2)) + Mst2*s2t*(189 + 726*lmMst2 + 32*(1 + lmMst2)*pow2(lmMst1) +
        707*pow2(lmMst2) - 2*lmMst1*(253 + 332*lmMst2 + 123*pow2(lmMst2)) +
        214*pow3(lmMst2)))*pow4(Mst1)*pow5(Mst2) + (Mst2*s2t*(205 + 252*lmMst2
        + 32*lmMst1*(1 + lmMst2) + 91*pow2(lmMst2)) + 128*Mt*pow2(1 + lmMst2))*
        pow2(Mst1)*pow7(Mst2) - 4*Dmglst2*(pow2(Mst2)*(64*Mt*(8 + 13*lmMst2 -
        8*pow2(lmMst2) + lmMst1*(-5 + 5*lmMst2 + 6*pow2(lmMst2)) - 6*pow3(
        lmMst2)) - Mst2*s2t*(48 + 4*lmMst2*(31 + 36*pow2(lmMst1)) + 278*pow2(
        lmMst2) - lmMst1*(44 + 278*lmMst2 + 379*pow2(lmMst2)) + 235*pow3(
        lmMst2)))*pow6(Mst1) + (16*Mt*(49 + 103*lmMst2 - 36*(1 + lmMst2)*pow2(
        lmMst2) + 4*lmMst1*(-11 + 6*lmMst2 + 9*pow2(lmMst2))) - Mst2*s2t*(48 +
        4*lmMst2*(45 + 68*pow2(lmMst1)) + 310*pow2(lmMst2) - lmMst1*(92 + 310*
        lmMst2 + 635*pow2(lmMst2)) + 363*pow3(lmMst2)))*pow8(Mst1)) - 2*Mst2*(
        pow2(Mst2)*(64*(1 + lmMst2)*Mt*(1 - 10*lmMst2 + 4*lmMst1*(2 + lmMst2) -
        4*pow2(lmMst2)) - Mst2*s2t*(32 + 285*lmMst2 + 144*(1 + lmMst2)*pow2(
        lmMst1) + 444*pow2(lmMst2) - lmMst1*(253 + 588*lmMst2 + 379*pow2(
        lmMst2)) + 235*pow3(lmMst2)))*pow6(Mst1) + (32*(1 + lmMst2)*Mt*(7 - 32*
        lmMst2 + 4*lmMst1*(7 + 3*lmMst2) - 12*pow2(lmMst2)) - Mst2*s2t*(40 +
        277*lmMst2 + 272*(1 + lmMst2)*pow2(lmMst1) + 556*pow2(lmMst2) - lmMst1*
        (237 + 828*lmMst2 + 635*pow2(lmMst2)) + 363*pow3(lmMst2)))*pow8(Mst1))
        - 2*Dmglst2*s2t*((20 - 2*(99 + 16*lmMst1)*lmMst2 - 123*pow2(lmMst2))*
        pow2(Mst1)*pow7(Mst2) + 32*lmMst2*(1 + lmMst2)*pow9(Mst2)) - 16*s2t*
        pow2(1 + lmMst2)*power10(Mst2)))))))/(1.000188e8*Tbeta*pow2(Sbeta)*
        pow4(Msq)*pow4(Mst1)*pow6(Mst2)));
}

/**
 * @return The off-diagonal (1, 2) = (2, 1) matrix element of the Higgs mass matrix as a double for the hierarchy 'H6bq22g'
 */
double H6bq22g::getS12() const {
   return (Al4p*Mt*MuSUSY*xDmglst2*z2*pow2(Dmglst2)*(-22680*twoLoopFlag*pow2(Mst1)*
        (2*Mt*pow2(Mst1)*pow2(Mst2)*(-24*Mt*MuSUSY*s2t + 4*Tbeta*pow2(Mt) +
        Mst2*(MuSUSY + 9*Mst2*Tbeta)*pow2(s2t)) + 2*Mt*(-36*Mt*MuSUSY*s2t + 28*
        Tbeta*pow2(Mt) + Mst2*(MuSUSY + 9*Mst2*Tbeta)*pow2(s2t))*pow4(Mst1) +
        s2t*(2*Mst2*Mt*s2t*(MuSUSY + 9*Mst2*Tbeta) - 24*MuSUSY*pow2(Mt) -
        Tbeta*pow2(s2t)*pow3(Mst2))*pow4(Mst2)) + Al4p*Mst2*threeLoopFlag*(
        pow4(Mst1)*(4*Mst2*s2t*(35*(-184517 + 74952*lmMst1 + 18360*lmMst2)*
        MuSUSY + 36*(-297316 + 20265*lmMst1 - 231945*lmMst2)*Mst2*Tbeta)*pow2(
        Mt) - 7*Mt*(4*(-261247 + 163890*lmMst1 - 140670*lmMst2)*MuSUSY + 3*(-
        1486429 + 190080*lmMst1 + 43200*lmMst2)*Mst2*Tbeta)*pow2(Mst2)*pow2(
        s2t) + 8*(21*(12761 + 5400*lmMst1 + 178920*lmMst2)*MuSUSY + 2*(-
        10736701 + 1553580*lmMst1 - 3303720*lmMst2 + 1508220*lmMt)*Mst2*Tbeta)*
        pow3(Mt) + 7*(-642497 + 170100*lmMst1 - 170100*lmMst2)*Tbeta*pow3(s2t)*
        pow4(Mst2)) - 3*pow2(Mst1)*pow2(Mst2)*(4*Mst2*s2t*(-28*(46987 + 15390*
        lmMst1 + 4050*lmMst2)*MuSUSY + 3*(255749 + 2100*lmMst1 + 441420*lmMst2)
        *Mst2*Tbeta)*pow2(Mt) + 14*Mt*((40001 + 52560*lmMst1 - 37080*lmMst2)*
        MuSUSY + 6*(46987 + 15390*lmMst1 + 4050*lmMst2)*Mst2*Tbeta)*pow2(Mst2)*
        pow2(s2t) - 32*(84*(17 + 1320*lmMst2)*MuSUSY + (-217444 + 53865*lmMst1
        - 116550*lmMst2 + 62685*lmMt)*Mst2*Tbeta)*pow3(Mt) - 7*(30641 + 52560*
        lmMst1 - 37080*lmMst2)*Tbeta*pow3(s2t)*pow4(Mst2)) + 196560*s2t*(-2*Mt*
        MuSUSY*s2t - 4*Tbeta*pow2(Mt) + Tbeta*pow2(Mst2)*pow2(s2t))*pow6(Mst2))
        ))/(34020.*Tbeta*pow2(Mst1)*pow7(Mst2)) + threeLoopFlag*z2*pow2(Al4p)*(
        MuSUSY*(-(pow2(Mt)*pow2(s2t)*((pow2(Mst1)*(Dmglst2*(866.0398148148148 -
        (775*lmMst1)/3. + 357*lmMst2 + (120*Dmsqst2*(1 - 2*lmMst1 + 2*lmMst2)*(
        pow3(Dmsqst2) + pow6(Msq)))/pow8(Msq)) + Mst2*(51.273148148148145 + 83*
        lmMst1 - (241*lmMst2)/3. + (40*((6 + 5*lmMst1 - 5*lmMst2)*pow4(Dmsqst2)
        + Dmsqst2*(3 + 2*lmMst1 - 2*lmMst2)*pow6(Msq)))/pow8(Msq))))/pow2(Mst2)
        + (Mst2*(2269 + 664*lmMst1 - 56*lmMst2 + (960*Dmsqst2*(2 + lmMst1 -
        lmMst2))/pow2(Msq) + (480*(13 + 5*lmMst1 - 5*lmMst2)*pow4(Dmsqst2))/
        pow8(Msq)))/24. + (pow4(Mst1)*(18.915432098765432 + 134*lmMst1 - (394*
        lmMst2)/3. + (20*Dmsqst2*(5 + 4*lmMst1 - 4*lmMst2)*(pow3(Dmsqst2) +
        pow6(Msq)))/pow8(Msq)))/pow3(Mst2) + Dmglst2*(483.18611111111113 - (
        209*lmMst1)/3. + 191*lmMst2 + (pow4(Mst1)*(1749.093621399177 - 498*
        lmMst1 + (1790*lmMst2)/3. + (20*Dmsqst2*(17 - 12*lmMst1 + 12*lmMst2)*(
        pow3(Dmsqst2) + pow6(Msq)))/pow8(Msq)))/pow4(Mst2) - (40*Dmsqst2*(2 +
        3*lmMst1 - 3*lmMst2)*(pow3(Dmsqst2) + pow6(Msq)))/pow8(Msq)))) + (s2t*
        pow3(Mt)*(129600*Dmglst2*Dmsqst2*pow2(Mst2)*(5*pow2(Mst1)*pow2(Mst2) +
        9*pow4(Mst1) + pow4(Mst2))*(pow3(Dmsqst2) + pow6(Msq)) + 100*Mst2*(
        pow4(Dmsqst2)*(-8883*pow2(Mst2)*pow4(Mst1) - 4212*pow2(Mst1)*pow4(Mst2)
        + 4057*pow6(Mst1) - 1620*pow6(Mst2)) + Dmsqst2*pow6(Msq)*(3249*pow2(
        Mst2)*pow4(Mst1) + 1431*pow2(Mst1)*pow4(Mst2) + 4057*pow6(Mst1) - 648*
        pow6(Mst2))) + Mst2*(3*(623599 + 65160*lmMst1 - 134280*lmMst2)*pow2(
        Mst2)*pow4(Mst1) + 9*(1367 + 13560*lmMst1 - 36600*lmMst2)*pow2(Mst1)*
        pow4(Mst2) + (5161826 + 406080*lmMst1 - 613440*lmMst2)*pow6(Mst1) -
        41040*pow6(Mst2))*pow8(Msq) - 24*Dmglst2*(6*(19579 + 1770*lmMst1 +
        12630*lmMst2)*pow2(Mst2)*pow4(Mst1) + 6*(4429 - 270*lmMst1 + 8910*
        lmMst2)*pow2(Mst1)*pow4(Mst2) + (510139 + 44280*lmMst1 + 76680*lmMst2)*
        pow6(Mst1) - 2520*pow6(Mst2))*pow8(Msq)))/(4860.*pow2(Mst1)*pow5(Mst2)*
        pow8(Msq)) + (pow4(Mt)*(4*(Dmglst2*(13463149 - 1492020*lmMst1 +
        2713500*lmMst2 - 625320*lmMt) + 3*(-1077991 + 184140*lmMst1 - 400140*
        lmMst2 + 8640*lmMt)*Mst2)*pow4(Mst1) + 405*Dmglst2*pow4(Mst2)*(3245 -
        1672*lmMst1 + 1448*lmMst2 - 1888*lmMt + (960*Dmsqst2*(8 - 3*lmMst1 + 6*
        lmMst2 - 3*lmMt)*(pow3(Dmsqst2) + pow6(Msq)))/pow8(Msq)) + (27*pow5(
        Mst2)*(7200*(-6 + 5*lmMst1 - 10*lmMst2 + 5*lmMt)*pow4(Dmsqst2) + 14400*
        Dmsqst2*(lmMst1 - 2*lmMst2 + lmMt)*pow6(Msq) + (-66773 + 9960*lmMst1 -
        45480*lmMst2 + 3840*lmMt)*pow8(Msq)))/pow8(Msq) + (72*pow2(Mst1)*pow2(
        Mst2)*(-54000*Dmglst2*Dmsqst2*(-3 + lmMst1 - 2*lmMst2 + lmMt)*(pow3(
        Dmsqst2) + pow6(Msq)) + 10800*Dmsqst2*(-1 + lmMst1 - 2*lmMst2 + lmMt)*
        Mst2*(pow3(Dmsqst2) + pow6(Msq)) + 3*Dmglst2*(67843 - 10050*lmMst1 +
        17490*lmMst2 - 7560*lmMt)*pow8(Msq) + 2*(-41936 + 7245*lmMst1 - 19665*
        lmMst2 + 720*lmMt)*Mst2*pow8(Msq)))/pow8(Msq)))/(7290.*pow6(Mst2)) +
        Mt*((s2t*(shiftst3*pow2(Mst2)*(-4*pow2(Mt) + ((-1 - lmMst1 + lmMst2)*
        pow2(Mst1) + pow2(Mst2))*pow2(s2t)) + (10*shiftst1*(-4*pow2(Mst2)*pow2(
        Mt) + 2*pow2(Mst1)*(2*pow2(Mt) + (-lmMst1 + lmMst2)*pow2(Mst2)*pow2(
        s2t)) - pow2(s2t)*pow4(Mst1) + pow2(s2t)*pow4(Mst2))*(pow4(Dmsqst2) +
        Dmsqst2*pow6(Msq) + pow8(Msq)))/pow8(Msq)))/(3.*pow2(Mst1)) + pow3(s2t)
        *((pow4(Mst1)*(Mst2*((47*lmMst1)/6. + (-3447379 - 76140*lmMst2 - (
        1082800*Dmsqst2)/pow2(Msq) + (171200*pow4(Dmsqst2))/pow8(Msq))/9720.) +
        Dmglst2*(432.76982167352537 - (398*lmMst1)/9. + (398*lmMst2)/9. + (130*
        Dmsqst2*(pow3(Dmsqst2) + pow6(Msq)))/(3.*pow8(Msq)))))/pow3(Mst2) +
        Dmglst2*(Mst2*(33.31481481481482 - 12*lmMst1 + (280*lmMst2)/9. - (20*
        Dmsqst2*(5 + 3*lmMst1 - 3*lmMst2)*(pow3(Dmsqst2) + pow6(Msq)))/(3.*
        pow8(Msq))) + (pow2(Mst1)*(120.97407407407407 - (254*lmMst1)/9. + (254*
        lmMst2)/9. + (20*Dmsqst2*(3 - 5*lmMst1 + 5*lmMst2)*(pow3(Dmsqst2) +
        pow6(Msq)))/(3.*pow8(Msq))))/Mst2) - pow2(Mst2)*(53.93518518518518 + (
        47*lmMst1)/9. - (92*lmMst2)/9. + (325*Dmsqst2)/(18.*pow2(Msq)) - (5*(7
        + lmMst1 - lmMst2)*pow4(Dmsqst2))/pow8(Msq)) + (pow2(Mst1)*(75*(119 +
        234*lmMst1 - 234*lmMst2)*pow4(Dmsqst2) + 150*Dmsqst2*(-469 + 36*lmMst1
        - 36*lmMst2)*pow6(Msq) + (-224837 + 4680*lmMst1 - 8370*lmMst2)*pow8(
        Msq)))/(810.*pow8(Msq)) + (pow3(Mst2)*(-4*Dmglst2*(15*pow4(Dmsqst2) +
        15*Dmsqst2*pow6(Msq) + 7*pow8(Msq)) + Mst2*(75*pow4(Dmsqst2) + 30*
        Dmsqst2*pow6(Msq) + 19*pow8(Msq))))/(9.*pow2(Mst1)*pow8(Msq))))) + (-(
        pow2(Mt)*pow2(MuSUSY)*(-(pow2(Mt)*(6*Mst2*(Dmglst2*(7286 - 240*lmMst1 +
        6384*lmMst2) + (-353 + 72*lmMst1 + 696*lmMst2)*Mst2)*pow2(Mst1) + 96*(
        24*Dmglst2*(5 + 6*lmMst2) + (53 + 24*lmMst2)*Mst2)*pow3(Mst2) + (-38401
        - 1080*lmMst1 + 7992*lmMst2)*pow4(Mst1)))/(81.*pow6(Mst2)) + (2*pow2(
        s2t)*(shiftst3*(2*(1 - lmMst1 + lmMst2) + (2*(1 - 2*lmMst1 + 2*lmMst2)*
        pow2(Mst1))/pow2(Mst2) + pow2(Mst2)/pow2(Mst1) + ((2 - 6*lmMst1 + 6*
        lmMst2)*pow4(Mst1))/pow4(Mst2)) + 10*shiftst1*(1 + pow2(Mst2)/pow2(
        Mst1) - (2*(lmMst1 - lmMst2)*(pow2(Mst1)*pow2(Mst2) + pow4(Mst1) +
        pow4(Mst2)))/pow4(Mst2))*(1 + Dmsqst2/pow2(Msq) + pow4(Dmsqst2)/pow8(
        Msq))))/3. - Mt*s2t*(((Dmglst2*(15057833 - 4014360*lmMst1 + 5563080*
        lmMst2) + 3*(266863 + 396360*lmMst1 - 346680*lmMst2)*Mst2)*pow4(Mst1))/
        (3645.*pow6(Mst2)) + (pow2(Mst1)*((4*Mst2*(3937 + 2988*lmMst1 - 2232*
        lmMst2 + (1080*Dmsqst2*(5 + 3*lmMst1 - 3*lmMst2)*(pow3(Dmsqst2) + pow6(
        Msq)))/pow8(Msq)))/81. + Dmglst2*(1798.967901234568 - (1312*lmMst1)/3.
         + (2192*lmMst2)/
        3. + (160*Dmsqst2*(1 - 9*lmMst1 + 9*lmMst2)*(pow3(Dmsqst2) + pow6(Msq))
        )/(3.*pow8(Msq)))))/pow4(Mst2) + (Dmglst2*(644.2481481481482 - (836*
        lmMst1)/9. + (764*lmMst2)/3. - (160*Dmsqst2*(2 + 3*lmMst1 - 3*lmMst2)*(
        pow3(Dmsqst2) + pow6(Msq)))/(3.*pow8(Msq))))/pow2(Mst2) + (480*(13 + 5*
        lmMst1 - 5*lmMst2)*pow4(Dmsqst2) + 960*Dmsqst2*(2 + lmMst1 - lmMst2)*
        pow6(Msq) + (2269 + 664*lmMst1 - 56*lmMst2)*pow8(Msq))/(18.*Mst2*pow8(
        Msq))))) + pow2(s2t)*((5*Dmsqst2*Mt*pow2(Sbeta)*pow4(Mst1)*((-6*(11 +
        10*lmMst1 - 10*lmMst2)*s2t*pow2(Mst2)*(-1 + pow2(Sbeta)) + (Dmglst2*(
        372 - 40*lmMst1 + 40*lmMst2) + 9*(-17 + 2*lmMst1 - 2*lmMst2)*Mst2)*Mt*
        pow2(Sbeta))*pow3(Dmsqst2) + 4*Dmglst2*(93 - 10*lmMst1 + 10*lmMst2)*Mt*
        pow2(Sbeta)*pow6(Msq)))/(6.*pow3(Mst2)*pow8(Msq)) + pow2(Mt)*((5*pow2(
        Sbeta)*pow4(Mst1)*(9*(17 - 2*lmMst1 + 2*lmMst2)*Mst2*pow4(Dmsqst2) + 4*
        Dmglst2*Dmsqst2*(-93 + 10*lmMst1 - 10*lmMst2)*(pow3(Dmsqst2) + pow6(
        Msq))))/(6.*pow3(Mst2)*pow8(Msq)) + pow2(MuSUSY)*(103.64814814814815 +
        (94*lmMst1)/9. - (184*lmMst2)/9. + (265*Dmsqst2)/(9.*pow2(Msq)) - (
        pow4(Mst1)*(Dmglst2*(1167.8951989026064 - (1520*lmMst1)/9. + (1864*
        lmMst2)/9.) - Mst2*(1368.138477366255 - (151*lmMst1)/9. + (143*lmMst2)/
        9. + (5*Dmsqst2*(20701 - 648*lmMst1 + 648*lmMst2)*(pow3(Dmsqst2) +
        pow6(Msq)))/(243.*pow8(Msq)))))/pow5(Mst2) + (8*Dmglst2*Mst2*(7 + (15*
        Dmsqst2*(pow3(Dmsqst2) + pow6(Msq)))/pow8(Msq)))/(9.*pow2(Mst1)) - (
        Dmglst2*(60.407407407407405 - 24*lmMst1 + (560*lmMst2)/9. - (40*
        Dmsqst2*(2 + lmMst1 - lmMst2)*(pow3(Dmsqst2) + pow6(Msq)))/pow8(Msq)))/
        Mst2 - (10*(26 + 3*lmMst1 - 3*lmMst2)*pow4(Dmsqst2))/(3.*pow8(Msq)) + (
        (-2*pow2(Mst2)*(75*pow4(Dmsqst2) + 30*Dmsqst2*pow6(Msq) + 19*pow8(Msq))
        )/(9.*pow2(Mst1)) + (pow2(Mst1)*(-150*(587 + 288*lmMst1 - 288*lmMst2)*
        Mst2*pow4(Dmsqst2) + 150*Dmsqst2*(1097 - 72*lmMst1 + 72*lmMst2)*Mst2*
        pow6(Msq) + 10800*Dmglst2*Dmsqst2*(3 + 8*lmMst1 - 8*lmMst2)*(pow3(
        Dmsqst2) + pow6(Msq)) + 36*Dmglst2*(-6803 + 1810*lmMst1 - 2670*lmMst2)*
        pow8(Msq) + (533629 - 900*lmMst1 + 180*lmMst2)*Mst2*pow8(Msq)))/(810.*
        pow3(Mst2)))/pow8(Msq)))))/Tbeta + (Mt*((5*xDmsqst2*pow2(Dmsqst2)*((
        s2t*(s2t*(81*Mst2*pow2(Sbeta)*(-2*(11 + 10*lmMst1 - 10*lmMst2)*s2t*
        pow2(Mst2)*(-1 + pow2(Sbeta)) + (Dmglst2*(372 - 40*lmMst1 + 40*lmMst2)
        + 3*(-17 + 2*lmMst1 - 2*lmMst2)*Mst2)*Mt*pow2(Sbeta))*pow6(Mst1) + Mt*(
        6*Mst2*pow2(MuSUSY)*(-6*(25 + 9*lmMst1 - 9*lmMst2)*pow2(Mst1)*pow3(
        Mst2) + (1607 - 432*lmMst1 + 432*lmMst2)*Mst2*pow4(Mst1) + 216*Dmglst2*
        (3*(2 + lmMst1 - lmMst2)*pow2(Mst1)*pow2(Mst2) + (3 + 8*lmMst1 - 8*
        lmMst2)*pow4(Mst1) + pow4(Mst2)) - 162*pow5(Mst2)) + ((41402 - 1296*
        lmMst1 + 1296*lmMst2)*pow2(MuSUSY) + 81*Mst2*(4*Dmglst2*(-93 + 10*
        lmMst1 - 10*lmMst2) + 3*(17 - 2*lmMst1 + 2*lmMst2)*Mst2)*pow2(Sbeta))*
        pow6(Mst1))) + 648*Mt*pow2(MuSUSY)*(8*Dmglst2*Mt*pow2(Mst1)*((1 - 9*
        lmMst1 + 9*lmMst2)*pow2(Mst1) + (-2 - 3*lmMst1 + 3*lmMst2)*pow2(Mst2))
        + (4*(7 + 3*lmMst1 - 3*lmMst2)*Mt + (-1 + 2*lmMst1 - 2*lmMst2)*Mst2*
        s2t*shiftst1)*pow2(Mst1)*pow3(Mst2) + 2*Mst2*(4*(5 + 3*lmMst1 - 3*
        lmMst2)*Mt + (lmMst1 - lmMst2)*Mst2*s2t*shiftst1)*pow4(Mst1) + s2t*
        shiftst1*(2*(lmMst1 - lmMst2)*pow6(Mst1) - pow6(Mst2)))))/Tbeta +
        MuSUSY*(-18*pow2(Mst1)*pow3(Mst2)*(2*Mst2*s2t*(25 - 36*shiftst1)*pow2(
        Mt) + 108*(7 + 3*lmMst1 - 3*lmMst2)*Mt*pow2(Mst2)*pow2(s2t) + 144*(2 -
        3*lmMst1 + 6*lmMst2 - 3*lmMt)*pow3(Mt) + (2 + lmMst2*(9 - 36*shiftst1)
        + 9*lmMst1*(-1 + 4*shiftst1))*pow3(Mst2)*pow3(s2t)) - 3*Mst2*(530*Mst2*
        s2t*pow2(Mt) + 1296*(4 + 3*lmMst1 - 3*lmMst2)*Mt*pow2(Mst2)*pow2(s2t) -
        3456*(-1 + lmMst1 - 2*lmMst2 + lmMt)*pow3(Mt) + (1757 - 378*lmMst1 +
        378*lmMst2 + 108*shiftst1)*pow3(Mst2)*pow3(s2t))*pow4(Mst1) + 2*s2t*(-
        972*(5 + 4*lmMst1 - 4*lmMst2)*Mst2*Mt*s2t + 4057*pow2(Mt) - 3324*pow2(
        Mst2)*pow2(s2t))*pow6(Mst1) + 162*s2t*(3 + 2*shiftst1)*(-4*pow2(Mt) +
        pow2(Mst2)*pow2(s2t))*pow6(Mst2) + 324*Dmglst2*(2*pow2(Mst1)*pow2(Mst2)
        *(20*Mst2*s2t*pow2(Mt) + 6*(2 + 3*lmMst1 - 3*lmMst2)*Mt*pow2(Mst2)*
        pow2(s2t) - 8*(-8 + 3*lmMst1 - 6*lmMst2 + 3*lmMt)*pow3(Mt) + (-5 - 3*
        lmMst1 + 3*lmMst2)*pow3(Mst2)*pow3(s2t)) - 2*(-36*Mst2*s2t*pow2(Mt) +
        18*(1 - 2*lmMst1 + 2*lmMst2)*Mt*pow2(Mst2)*pow2(s2t) + 80*(-3 + lmMst1
        - 2*lmMst2 + lmMt)*pow3(Mt) + (-3 + 5*lmMst1 - 5*lmMst2)*pow3(Mst2)*
        pow3(s2t))*pow4(Mst1) + 8*s2t*pow2(Mt)*pow5(Mst2) + (6*(-17 + 12*lmMst1
        - 12*lmMst2)*Mt + 13*Mst2*s2t)*pow2(s2t)*pow6(Mst1) - 2*pow3(s2t)*pow7(
        Mst2)))))/(486.*pow4(Msq)*pow4(Mst2)) + (4*MuSUSY*xDR2DRMOD*(pow3(Mst2)
        *(64*(1 + lmMst2)*(2*MuSUSY*s2t + 3*Mt*Tbeta)*pow2(Mt) - 4*Mst2*Mt*((4
        + 13*lmMst1 - 9*lmMst2)*MuSUSY + 12*(1 + lmMst2)*Mst2*Tbeta)*pow2(s2t)
        + 13*Tbeta*pow3(Mst2)*pow3(s2t))*pow4(Mst1) + 2*pow2(Mst1)*((32*(1 +
        lmMst2)*MuSUSY*s2t - 26*Mst2*s2t*Tbeta)*pow2(Mt) + Mst2*Mt*((5 - 26*
        lmMst1 + 18*lmMst2)*MuSUSY - 24*(1 + lmMst2)*Mst2*Tbeta)*pow2(s2t) +
        32*(1 + lmMst2)*Tbeta*pow3(Mt) + (4 + 13*lmMst1 - 9*lmMst2)*Tbeta*pow3(
        Mst2)*pow3(s2t))*pow5(Mst2) - 4*Mst2*xDmglst2*pow2(Dmglst2)*(2*pow2(
        Mst1)*pow2(Mst2)*(-2*s2t*(4*(-2 + 3*lmMst2)*MuSUSY + 5*Mst2*Tbeta)*
        pow2(Mt) + Mst2*Mt*((1 - 10*lmMst1 + 12*lmMst2)*MuSUSY + 6*(-2 + 3*
        lmMst2)*Mst2*Tbeta)*pow2(s2t) + 8*Tbeta*pow3(Mt) + (2 + 5*lmMst1 - 6*
        lmMst2)*Tbeta*pow3(Mst2)*pow3(s2t)) + (32*(2 - 3*lmMst2)*MuSUSY*s2t*
        pow2(Mt) + 4*Mst2*Mt*((-2 - 5*lmMst1 + 6*lmMst2)*MuSUSY + 3*(-2 + 3*
        lmMst2)*Mst2*Tbeta)*pow2(s2t) - 16*(7 + 2*lmMst2)*Tbeta*pow3(Mt) + 5*
        Tbeta*pow3(Mst2)*pow3(s2t))*pow4(Mst1) - 5*s2t*(-2*Mt*MuSUSY*s2t - 4*
        Tbeta*pow2(Mt) + Tbeta*pow2(Mst2)*pow2(s2t))*pow5(Mst2)) + 4*Mst2*Mt*(
        48*(1 + lmMst2)*Mt*(MuSUSY*s2t + 2*Mt*Tbeta) - Mst2*((4 + 13*lmMst1 -
        9*lmMst2)*MuSUSY + 12*(1 + lmMst2)*Mst2*Tbeta)*pow2(s2t))*pow6(Mst1) +
        Dmglst2*(-(pow2(Mst2)*(-128*(1 + 3*lmMst2)*MuSUSY*s2t*pow2(Mt) + 4*
        Mst2*Mt*(-7*lmMst1*MuSUSY + 15*lmMst2*MuSUSY + 12*Mst2*Tbeta + 36*
        lmMst2*Mst2*Tbeta)*pow2(s2t) + 64*(7 + lmMst2)*Tbeta*pow3(Mt) + 7*
        Tbeta*pow3(Mst2)*pow3(s2t))*pow4(Mst1)) + 2*pow2(Mst1)*(2*s2t*(16*(1 +
        3*lmMst2)*MuSUSY + 7*Mst2*Tbeta)*pow2(Mt) - Mst2*Mt*((7 - 14*lmMst1 +
        30*lmMst2)*MuSUSY + 24*(1 + 3*lmMst2)*Mst2*Tbeta)*pow2(s2t) + 32*(-1 +
        lmMst2)*Tbeta*pow3(Mt) + (-7*lmMst1 + 15*lmMst2)*Tbeta*pow3(Mst2)*pow3(
        s2t))*pow4(Mst2) - 4*Mt*(-48*(1 + 3*lmMst2)*Mt*MuSUSY*s2t + 32*(11 + 5*
        lmMst2)*Tbeta*pow2(Mt) + Mst2*(-7*lmMst1*MuSUSY + 15*lmMst2*MuSUSY +
        12*Mst2*Tbeta + 36*lmMst2*Mst2*Tbeta)*pow2(s2t))*pow6(Mst1) + 7*s2t*(-
        2*Mt*MuSUSY*s2t - 4*Tbeta*pow2(Mt) + Tbeta*pow2(Mst2)*pow2(s2t))*pow7(
        Mst2)) - 13*s2t*(-2*Mt*MuSUSY*s2t - 4*Tbeta*pow2(Mt) + Tbeta*pow2(Mst2)
        *pow2(s2t))*pow8(Mst2)))/(9.*Tbeta*pow6(Mst2))))/pow2(Mst1)) + MuSUSY*(
        (Mt*oneLoopFlag*s2t*(4*(lmMst1 - lmMst2)*pow2(Mt) - ((2 + lmMst1 -
        lmMst2)*pow2(Mst1) + (-2 + lmMst1 - lmMst2)*pow2(Mst2))*pow2(s2t) - (2*
        Mt*MuSUSY*s2t*(2 - lmMst1 + lmMst2 - (2*(lmMst1 - lmMst2)*pow2(Mst1)*(
        pow2(Mst1) + pow2(Mst2)))/pow4(Mst2)))/Tbeta))/16. - (Al4p*Mt*
        xDR2DRMOD*(-270*Al4p*s2t*threeLoopFlag*xDmsqst2*pow2(Dmsqst2)*pow2(
        Mst1)*pow2(Mst2)*pow4(Msq)*((-2*lmMst1 + 2*lmMst2)*s2t*pow2(Mst1)*(-2*
        Mt*MuSUSY*pow2(Mst1)*(4*Dmglst2*Mst2 + 26*xDmglst2*pow2(Dmglst2) + (-11
        + 18*shiftst1)*pow2(Mst2)) + pow2(Mst2)*(-10*Dmglst2*Mst2 + 7*xDmglst2*
        pow2(Dmglst2) + (-11 + 18*shiftst1)*pow2(Mst2))*(-2*Mt*MuSUSY + s2t*
        Tbeta*pow2(Mst2)) - 4*Mt*MuSUSY*(-5 + 9*shiftst1)*pow4(Mst1)) + pow2(
        Mst2)*(10*Dmglst2*Mst2 - 7*xDmglst2*pow2(Dmglst2) + (11 - 18*shiftst1)*
        pow2(Mst2))*(2*Mt*(MuSUSY*s2t - 2*Mt*Tbeta)*pow2(Mst1) + pow2(Mst2)*(2*
        Mt*MuSUSY*s2t + 4*Tbeta*pow2(Mt) - Tbeta*pow2(Mst2)*pow2(s2t)) + Tbeta*
        pow2(s2t)*pow4(Mst1))) + 216*Mst2*(2*Dmglst2*lmMst2 + Mst2 + lmMst2*
        Mst2)*s2t*twoLoopFlag*pow2(Mst1)*(-((2*Mt*(MuSUSY*s2t - 2*Mt*Tbeta)*
        pow2(Mst1) + pow2(Mst2)*(2*Mt*MuSUSY*s2t + 4*Tbeta*pow2(Mt) - Tbeta*
        pow2(Mst2)*pow2(s2t)) + Tbeta*pow2(s2t)*pow4(Mst1))*pow4(Mst2)) + 2*(
        lmMst1 - lmMst2)*s2t*pow2(Mst1)*(2*Mt*MuSUSY*(pow2(Mst1)*pow2(Mst2) +
        pow4(Mst1) + pow4(Mst2)) - s2t*Tbeta*pow6(Mst2)))*pow8(Msq) - xDmglst2*
        pow2(Dmglst2)*(-216*(2 - lmMst2)*s2t*twoLoopFlag*pow2(Mst1)*(Tbeta*(-4*
        pow2(Mst1)*pow2(Mt) + 4*pow2(Mst2)*pow2(Mt) + 2*(lmMst1 - lmMst2)*pow2(
        Mst1)*pow2(Mst2)*pow2(s2t) + pow2(s2t)*(pow4(Mst1) - pow4(Mst2)))*pow4(
        Mst2) + 2*Mt*MuSUSY*s2t*((pow2(Mst1) + pow2(Mst2))*pow4(Mst2) - 2*(
        lmMst1 - lmMst2)*pow2(Mst1)*(pow2(Mst1)*pow2(Mst2) + pow4(Mst1) + pow4(
        Mst2))))*pow8(Msq) + Al4p*Mst2*threeLoopFlag*(12*Mt*MuSUSY*s2t*(128*Mt*
        (-(pow2(Mst2)*(5 + 149*lmMst2 + 12*pow2(lmMst2) + 6*lmMst1*(-7 - 8*
        lmMst2 + 6*pow2(lmMst2)) - 36*pow3(lmMst2))) + pow2(Mst1)*(85 - 215*
        lmMst2 + lmMst1*(66 + 90*lmMst2 - 72*pow2(lmMst2)) - 54*pow2(lmMst2) +
        72*pow3(lmMst2)))*pow4(Mst1)*pow8(Msq) + 3*pow2(Mst1)*pow4(Mst2)*(195*
        Mst2*s2t*(pow4(Dmsqst2) - 2*Dmsqst2*pow6(Msq)) - 2*(Mst2*s2t*(380 - 32*
        lmMst1*(-2 + lmMst2) + 38*lmMst2 - 251*pow2(lmMst2)) + 128*Mt*(3 + 4*
        lmMst2 - 3*pow2(lmMst2)))*pow8(Msq)) + Mst2*s2t*(45*(52*(lmMst1 -
        lmMst2)*pow2(Mst1) + (31 - 62*lmMst1 + 62*lmMst2)*pow2(Mst2))*pow4(
        Dmsqst2)*pow4(Mst1) - 1170*Dmsqst2*(pow2(Mst2) - 2*(lmMst1 - lmMst2)*(
        pow2(Mst1) + pow2(Mst2)))*pow4(Mst1)*pow6(Msq) - 2*(pow2(Mst2)*(1540 +
        2506*lmMst2 - 96*(-2 + lmMst2)*pow2(lmMst1) - 141*pow2(lmMst2) + 6*
        lmMst1*(-460 - 6*lmMst2 + 123*pow2(lmMst2)) - 642*pow3(lmMst2))*pow4(
        Mst1) + 6*(336 + 716*lmMst2 - 144*(-2 + lmMst2)*pow2(lmMst1) + 246*
        pow2(lmMst2) + lmMst1*(-668 - 534*lmMst2 + 379*pow2(lmMst2)) - 235*
        pow3(lmMst2))*pow6(Mst1) - 96*(2 + lmMst2 - 3*pow2(lmMst2))*pow6(Mst2))
        *pow8(Msq))) + Tbeta*(512*((215 + 37*lmMst2 + 3*(-5 + 8*lmMst2)*lmMt -
        18*lmMst1*(1 - 2*lmMst2 + lmMt) - 42*pow2(lmMst2))*pow2(Mst2) + pow2(
        Mst1)*(1097 + lmMst2*(1531 - 246*lmMt) - 339*lmMt + (444 - 36*lmMt)*
        pow2(lmMst2) - 18*lmMst1*(39 - 2*lmMst2*(-9 + lmMt) - 7*lmMt + 2*pow2(
        lmMst2)) + 36*pow3(lmMst2)))*pow3(Mt)*pow4(Mst1)*pow8(Msq) - 1152*Mt*
        pow2(Mst1)*pow2(Mst2)*pow2(s2t)*(pow2(Mst1)*pow2(Mst2)*(13 - 125*lmMst2
        + 6*lmMst1*(7 + 8*lmMst2 - 6*pow2(lmMst2)) - 30*pow2(lmMst2) + 36*pow3(
        lmMst2)) + 6*(15 - 11*lmMst2 + lmMst1*(4 + 7*lmMst2 - 6*pow2(lmMst2)) -
        7*pow2(lmMst2) + 6*pow3(lmMst2))*pow4(Mst1) + 6*(-3 - 4*lmMst2 + 3*
        pow2(lmMst2))*pow4(Mst2))*pow8(Msq) + 72*pow2(Mst1)*pow2(Mt)*pow4(Mst2)
        *(195*Mst2*s2t*(pow4(Dmsqst2) - 2*Dmsqst2*pow6(Msq)) - 2*(Mst2*s2t*(220
        - 32*lmMst1*(-2 + lmMst2) + 70*lmMst2 - 59*pow2(lmMst2)) + 64*Mt*(3 +
        4*lmMst2 - 3*pow2(lmMst2)))*pow8(Msq)) - 3*pow3(Mst2)*pow3(s2t)*(2340*
        Dmsqst2*pow2(Mst1)*(2*(lmMst1 - lmMst2)*pow2(Mst1)*pow2(Mst2) + pow4(
        Mst1) - pow4(Mst2))*pow6(Msq) - 45*pow4(Dmsqst2)*(52*(lmMst1 - lmMst2)*
        pow2(Mst2)*pow4(Mst1) - 29*pow2(Mst1)*pow4(Mst2) + 62*pow6(Mst1)) - 4*(
        2*pow2(Mst2)*(200 + 1196*lmMst2 - 48*(-2 + lmMst2)*pow2(lmMst1) + 306*
        pow2(lmMst2) + 3*lmMst1*(-492 + 10*lmMst2 + 123*pow2(lmMst2)) - 321*
        pow3(lmMst2))*pow4(Mst1) + 3*(444 - 32*lmMst1*(-2 + lmMst2) + 70*lmMst2
        - 347*pow2(lmMst2))*pow2(Mst1)*pow4(Mst2) + (476 + 1790*lmMst2 - 768*(-
        2 + lmMst2)*pow2(lmMst1) + 1617*pow2(lmMst2) + 96*lmMst1*(-13 - 33*
        lmMst2 + 16*pow2(lmMst2)) - 768*pow3(lmMst2))*pow6(Mst1) - 96*(2 +
        lmMst2 - 3*pow2(lmMst2))*pow6(Mst2))*pow8(Msq)) - 24*Mst2*s2t*pow2(Mt)*
        (45*pow2(Mst2)*pow4(Mst1)*(31*pow4(Dmsqst2) - 26*Dmsqst2*pow6(Msq)) -
        2*(pow2(Mst2)*(964 + 306*lmMst2 + 96*(-2 + lmMst2)*pow2(lmMst1) - 81*
        pow2(lmMst2) + 96*(lmMst1 + 3*lmMst1*lmMst2 - 2*lmMst1*pow2(lmMst2) +
        pow3(lmMst2)))*pow4(Mst1) + 96*(26 + 29*lmMst2 + (-2 + lmMst2)*pow2(
        lmMst1) + 4*pow2(lmMst2) - 2*lmMst1*(8 + lmMst2 + pow2(lmMst2)) + pow3(
        lmMst2))*pow6(Mst1) + 96*(2 + lmMst2 - 3*pow2(lmMst2))*pow6(Mst2))*
        pow8(Msq)))))))/(648.*Tbeta*pow4(Mst1)*pow6(Mst2)*pow8(Msq)) + (Al4p*
        Mt*twoLoopFlag*(12*Mt*pow2(s2t)*(2*(Mst2*(4 + 3*lmMst2 - lmMst1*(1 +
        lmMst2) + pow2(lmMst2)) + ((Dmglst2*(5 + lmMst1*(3 - 2*lmMst2) - 3*
        lmMst2 + 2*pow2(lmMst2)) + Mst2*(1 + 3*lmMst2 - lmMst1*(3 + 2*lmMst2) +
        2*pow2(lmMst2)))*pow2(Mst1))/pow2(Mst2)) + ((0.5 + 7*lmMst2 - lmMst1*(7
        + 4*lmMst2) + 4*pow2(lmMst2))*pow4(Mst1))/pow3(Mst2) + Dmglst2*(12 - 2*
        lmMst1*lmMst2 + 2*(lmMst1 + lmMst2 + pow2(lmMst2)) + ((10.5 + lmMst1*(9
        - 4*lmMst2) - 9*lmMst2 + 4*pow2(lmMst2))*pow4(Mst1))/pow4(Mst2))) + 4*
        s2t*pow2(Mt)*(12*lmMst1 + 4*(-2 + lmMst1)*lmMst2 + 2*pow2(lmMst1) - 6*
        pow2(lmMst2) - (4*(1 + lmMst2)*pow2(Mst2))/pow2(Mst1) + (Dmglst2*(4 +
        lmMst2*(8 - (8*pow2(Mst2))/pow2(Mst1))))/Mst2 + (4*(1 - lmMst1 +
        lmMst2)*(2*Dmglst2 + (lmMst1 - lmMst2)*Mst2)*pow2(Mst1))/pow3(Mst2) - (
        4*(Dmglst2*(-2 + 4*lmMst1 - 4*lmMst2) + Mst2*(-lmMst1 + lmMst2 - 2*
        lmMst1*lmMst2 + pow2(lmMst1) + pow2(lmMst2)))*pow4(Mst1))/pow5(Mst2)) -
        (16*pow3(Mt)*(Dmglst2*(2*(-7 + lmMst2*(-7 + lmMt) - lmMst1*(-5 + lmMt)
        + 2*lmMt)*pow2(Mst1)*pow2(Mst2) + (-11 + lmMst1*(32 + 6*lmMst2 - 8*
        lmMt) + 8*lmMst2*(-5 + lmMt) + 8*lmMt - 6*pow2(lmMst2))*pow4(Mst1) - 2*
        (5 + lmMst1*(-1 + lmMst2) + lmMst2 - pow2(lmMst2))*pow4(Mst2)) - Mst2*(
        2*(4 + lmMst1*(3 + 2*lmMst2 - lmMt) + lmMst2*(-4 + lmMt) + lmMt - 2*
        pow2(lmMst2))*pow2(Mst1)*pow2(Mst2) + (13 + 2*lmMst1*(6 + 3*lmMst2 - 2*
        lmMt) + 2*lmMt + 2*lmMst2*(-7 + 2*lmMt) - 6*pow2(lmMst2))*pow4(Mst1) +
        2*(2 + lmMst1 - 2*lmMst2 + lmMst1*lmMst2 + lmMt - pow2(lmMst2))*pow4(
        Mst2))))/pow6(Mst2) + ((pow3(s2t)*(2*(-16 + 6*lmMst2 - 2*lmMst1*(8 + 5*
        lmMst2) + 3*pow2(lmMst1) + 7*pow2(lmMst2))*pow3(Mst2)*pow4(Mst1) - 2*(-
        12 - 18*lmMst2 + 2*lmMst1*(5 + 3*lmMst2) + pow2(lmMst1) - 7*pow2(
        lmMst2))*pow2(Mst1)*pow5(Mst2) + Mst2*(3 + lmMst1*(2 - 32*lmMst2) - 2*
        lmMst2 + 16*pow2(lmMst1) + 16*pow2(lmMst2))*pow6(Mst1) - 4*Dmglst2*(2*(
        1 + 2*lmMst1 - lmMst2)*pow2(Mst2)*pow4(Mst1) + 2*(1 + lmMst1 - lmMst2 +
        2*lmMst1*lmMst2 - 2*pow2(lmMst2))*pow2(Mst1)*pow4(Mst2) + (1 + 2*lmMst1
        - 2*lmMst2)*pow6(Mst1) - 2*lmMst2*pow6(Mst2)) + 4*(1 + lmMst2)*pow7(
        Mst2)))/pow3(Mst2) - (2*Mt*MuSUSY*s2t*(4*(4*Mt*(5 + 6*lmMst2 - lmMst1*(
        4 + 3*lmMst2) + 3*pow2(lmMst2)) + Mst2*s2t*(-1 + 13*lmMst2 - lmMst1*(13
        + 8*lmMst2) + pow2(lmMst1) + 7*pow2(lmMst2)))*pow3(Mst2)*pow4(Mst1) +
        2*(-(Mst2*s2t*(-14 - 20*lmMst2 + 2*lmMst1*(5 + 3*lmMst2) + pow2(lmMst1)
        - 7*pow2(lmMst2))) + 8*Mt*(4 + 3*lmMst2 - lmMst1*(1 + lmMst2) + pow2(
        lmMst2)))*pow2(Mst1)*pow5(Mst2) + 4*Dmglst2*(Mst2*s2t*(-5 + 8*lmMst2 -
        4*lmMst1*(2 + lmMst2) + 4*pow2(lmMst2)) + Mt*(65 + lmMst1*(34 - 20*
        lmMst2) - 26*lmMst2 + 20*pow2(lmMst2)))*pow6(Mst1) + Mst2*(Mst2*s2t*(-1
        + 50*lmMst2 - 2*lmMst1*(25 + 32*lmMst2) + 20*pow2(lmMst1) + 44*pow2(
        lmMst2)) + Mt*(84 + 152*lmMst2 - 40*lmMst1*(3 + 2*lmMst2) + 80*pow2(
        lmMst2)))*pow6(Mst1) + 8*Dmglst2*((Mst2*s2t*(-2 + 3*lmMst2 - lmMst1*(3
        + 2*lmMst2) + 2*pow2(lmMst2)) + Mt*(22 + lmMst1*(8 - 6*lmMst2) - 4*
        lmMst2 + 6*pow2(lmMst2)))*pow2(Mst2)*pow4(Mst1) + (-(Mst2*s2t*(1 +
        lmMst1 - 2*lmMst2 + 2*lmMst1*lmMst2 - 2*pow2(lmMst2))) + 2*Mt*(6 +
        lmMst1 + lmMst2 - lmMst1*lmMst2 + pow2(lmMst2)))*pow2(Mst1)*pow4(Mst2)
        + lmMst2*s2t*pow7(Mst2)) + 4*(1 + lmMst2)*s2t*pow8(Mst2)))/(Tbeta*pow6(
        Mst2)))/pow2(Mst1)))/12. + Al4p*xDmglst2*pow2(Dmglst2)*((Mt*
        twoLoopFlag*(-2*(pow2(Mst2)*(36*(5 - 3*lmMst1 + 3*lmMst2)*Mst2*s2t*
        pow2(Mt) - 54*(12 + lmMst1 - lmMst2)*Mt*pow2(Mst2)*pow2(s2t) + 4*(155 +
        3*lmMst2*(41 - 6*lmMt) - 18*lmMst1*(4 + lmMst2 - lmMt) - 51*lmMt + 18*
        pow2(lmMst2))*pow3(Mt) - 9*(4*lmMst1 - 5*lmMst2)*pow3(Mst2)*pow3(s2t))*
        pow4(Mst1) + pow2(Mst1)*(6*(17 - 6*lmMst2)*Mst2*s2t*pow2(Mt) - 54*(8 +
        lmMst1 - lmMst2)*Mt*pow2(Mst2)*pow2(s2t) + 4*(11 + 18*lmMst1 - 21*
        lmMst2 + 3*lmMt)*pow3(Mt) + 9*(6 + 5*lmMst2 + lmMst1*(-5 + 2*lmMst2) -
        2*pow2(lmMst2))*pow3(Mst2)*pow3(s2t))*pow4(Mst2)) - (72*(9 - 10*lmMst1
        + 10*lmMst2)*Mst2*s2t*pow2(Mt) + 54*(-27 + 4*lmMst1 - 4*lmMst2)*Mt*
        pow2(Mst2)*pow2(s2t) + 8*(335 + 699*lmMst2 - 18*lmMst1*(29 + 7*lmMst2 -
        7*lmMt) - 3*(59 + 42*lmMst2)*lmMt + 126*pow2(lmMst2))*pow3(Mt) + 9*(1 -
        10*lmMst1 + 10*lmMst2)*pow3(Mst2)*pow3(s2t))*pow6(Mst1) + 18*(-2 +
        lmMst2)*s2t*(-4*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow7(Mst2) + (18*Mt*
        MuSUSY*s2t*(2*(-8*(10 + lmMst1 - lmMst2)*Mt + Mst2*s2t*(8 - 9*lmMst1 +
        9*lmMst2 + 2*lmMst1*lmMst2 - 2*pow2(lmMst2)))*pow2(Mst2)*pow4(Mst1) +
        2*(-4*(8 + lmMst1 - lmMst2)*Mt + Mst2*s2t*(8 - 5*lmMst1 + 4*lmMst2 + 2*
        lmMst1*lmMst2 - 2*pow2(lmMst2)))*pow2(Mst1)*pow4(Mst2) + (-268*Mt +
        Mst2*s2t*(17 + 4*lmMst1*(-7 + lmMst2) + 28*lmMst2 - 4*pow2(lmMst2)))*
        pow6(Mst1) - 2*(-2 + lmMst2)*s2t*pow7(Mst2)))/Tbeta))/(54.*pow2(Mst1)*
        pow7(Mst2)) + Al4p*threeLoopFlag*(-(pow2(Mt)*pow2(s2t)*((32*Mst2*(3 +
        4*lmMst2 - 3*pow2(lmMst2)))/(3.*pow2(Mst1)) + (380.51604938271606 - 76*
        B4 + 2*DN - 292*lmMst1 + (2*(631 - 438*lmMst1)*lmMst2)/9. - (16*pow2(
        lmMst1))/3. + (4*(53 + 48*lmMst1)*pow2(lmMst2))/3. - 64*pow3(lmMst2))/
        Mst2 + ((28*OepS2*(415*pow2(Mst1) + 24*pow2(Mst2)))/243. - (S2*((21422
        + 87150*lmMst1 - 87150*lmMst2)*pow2(Mst1) + 9*(-1373 + 560*lmMst1 -
        560*lmMst2)*pow2(Mst2)))/90.)/pow3(Mst2) + (pow2(Mst1)*(
        99.93816872427983 - 76*B4 + 2*DN - (6424*lmMst1)/27. - (4*(293 + 1152*
        lmMst1)*lmMst2)/27. + (92*pow2(lmMst1))/3. + 4*(35 + 16*lmMst1)*pow2(
        lmMst2) - 64*pow3(lmMst2)))/pow3(Mst2))) + s2t*pow3(Mt)*((2*(652 - 32*
        lmMst1*(-2 + lmMst2) + 70*lmMst2 - 59*pow2(lmMst2)))/(9.*pow2(Mst1)) +
        (1288.4830687830688 - (16*B4)/9. + (16*D3)/9. - (8*DN)/9. + (184*
        lmMst1)/3. - (32*lmMt)/9. + (2*lmMst2*(6695 + 1002*lmMst1 - 96*pow2(
        lmMst1)))/27. + (82*pow2(lmMst1))/9. - (8*(47 + 8*lmMst1)*pow2(lmMst2))
        /9. + (128*pow3(lmMst2))/9.)/pow2(Mst2) - (64*(2 + lmMst2 - 3*pow2(
        lmMst2))*pow2(Mst2))/(9.*pow4(Mst1)) - S2*((89.65714285714286 + 84*
        lmMst1 - 84*lmMst2)/pow2(Mst2) + ((847.3523809523809 + 480*lmMst1 -
        480*lmMst2)*pow2(Mst1))/pow4(Mst2)) + (16*OepS2*(40*pow2(Mst1) + 7*
        pow2(Mst2)))/(27.*pow4(Mst2)) + (pow2(Mst1)*(1803.9172204585539 - (16*
        B4)/9. + (16*D3)/9. - (8*DN)/9. + (225124*lmMst1)/675. + (32*(17 - 5*
        lmMst1 + 5*lmMst2)*lmMt)/9. + (4*lmMst2*(136544 + 53535*lmMst1 - 1425*
        pow2(lmMst1)))/675. + (224*pow2(lmMst1))/15. - (4*(3257 + 545*lmMst1)*
        pow2(lmMst2))/45. + (4*pow3(lmMst1))/9. + (508*pow3(lmMst2))/9.))/pow4(
        Mst2)) + Mt*(pow3(s2t)*(46.46594650205761 - (38*B4)/9. + (40*D3)/9. - (
        7*DN)/3. + (98*lmMst1)/3. - (25*pow2(lmMst1))/18. + (lmMst2*(-1493 -
        24*lmMst1 + 48*pow2(lmMst1)))/27. - ((247 + 246*lmMst1)*pow2(lmMst2))/
        18. + (4*OepS2*(21 + (1157*pow2(Mst1))/pow2(Mst2)))/729. + ((-876 + 32*
        lmMst1*(-2 + lmMst2) - 70*lmMst2 + 347*pow2(lmMst2))*pow2(Mst2))/(18.*
        pow2(Mst1)) - (S2*((648463 + 34710*lmMst1 - 34710*lmMst2)*pow2(Mst1) +
        9*(65819 + 70*lmMst1 - 70*lmMst2)*pow2(Mst2)))/(270.*pow2(Mst2)) - (
        pow2(Mst1)*(12.389529492455418 + (32*B4)/9. - (32*D3)/9. + (16*DN)/9. +
        (138661*lmMst1)/4050. + (159*pow2(lmMst1))/10. - lmMst2*(
        26.7558024691358 + (229*lmMst1)/5. + (143*pow2(lmMst1))/9.) + ((1333 +
        1355*lmMst1)*pow2(lmMst2))/45. + (5*pow3(lmMst1))/9. - (133*pow3(
        lmMst2))/9.))/pow2(Mst2) + (107*pow3(lmMst2))/9. + (16*(2 + lmMst2 - 3*
        pow2(lmMst2))*pow4(Mst2))/(9.*pow4(Mst1))) - (2*T1ep*(21*pow2(Mst2)*(
        36*Mst2*s2t*pow2(Mt) - 24*Mt*pow2(Mst2)*pow2(s2t) + 32*pow3(Mt) + pow3(
        Mst2)*pow3(s2t)) + pow2(Mst1)*(4320*Mst2*s2t*pow2(Mt) - 8715*Mt*pow2(
        Mst2)*pow2(s2t) + 17584*pow3(Mt) + 1157*pow3(Mst2)*pow3(s2t))))/(243.*
        pow5(Mst2))) + (2*(-3*(1648637 + 1173060*lmMt - 15680*OepS2 + 1680*
        lmMst2*(-1214 + 42*lmMt - 189*S2) - 953208*S2 + 7560*lmMst1*(131 + 33*
        lmMst2 - 8*lmMt + 42*S2) + 30240*pow2(lmMst1) - 304920*pow2(lmMst2) +
        15120*pow2(lmMt))*pow2(Mst1)*pow2(Mst2) + 2*(-1657412 - 7512750*lmMt +
        615440*OepS2 + 4340358*S2 + 11340*(35 + 4*lmMst2 - 4*lmMt)*pow2(lmMst1)
        + 945*(-1981 + 48*lmMt)*pow2(lmMst2) + 315*lmMst1*(-4223 + 2679*lmMst2
        + 201*lmMt - 39564*S2 + 432*pow2(lmMst2) - 432*pow2(lmMt)) + 385560*
        pow2(lmMt) + 945*lmMst2*(5321 + 193*lmMt + 13188*S2 + 144*pow2(lmMt)) -
        181440*pow3(lmMst2))*pow4(Mst1) - 181440*(-3 - 4*lmMst2 + 3*pow2(
        lmMst2))*pow4(Mst2))*pow4(Mt))/(25515.*pow2(Mst1)*pow5(Mst2)) - (
        MuSUSY*pow2(Mt)*(75*pow2(Mst2)*(-8*Mst2*Mt*s2t*(334138 - 61560*B4 +
        1620*DN - 236520*lmMst1 - 180*(-823 + 438*lmMst1)*lmMst2 + 2240*OepS2 -
        81*(-1373 + 560*lmMst1 - 560*lmMst2)*S2 - 3360*T1ep - 4320*pow2(lmMst1)
        + 1080*(29 + 48*lmMst1)*pow2(lmMst2) - 51840*pow3(lmMst2)) + 96*pow2(
        Mt)*(33934 - 90*B4 + 90*D3 - 45*DN + 720*lmMst1 + 120*(163 + 24*lmMst1)
        *lmMst2 - 120*lmMt - 10206*S2 - 720*(2 + lmMst1)*pow2(lmMst2) + 720*
        pow3(lmMst2)) + pow2(Mst2)*pow2(s2t)*(13169 - 41040*B4 + 43200*D3 -
        22680*DN + 282960*lmMst1 + 1120*OepS2 - 324*(65819 + 70*lmMst1 - 70*
        lmMst2)*S2 - 1680*T1ep - 13500*pow2(lmMst1) + 720*lmMst2*(-775 + 12*
        lmMst1 + 24*pow2(lmMst1)) - 1080*(-2 + 123*lmMst1)*pow2(lmMst2) +
        115560*pow3(lmMst2)))*pow4(Mst1) - 40500*s2t*(Mst2*s2t*(812 - 32*
        lmMst1*(-2 + lmMst2) + 38*lmMst2 - 251*pow2(lmMst2)) + 128*Mt*(3 + 4*
        lmMst2 - 3*pow2(lmMst2)))*pow2(Mst1)*pow5(Mst2) + 2*(-2*pow2(Mst2)*
        pow2(s2t)*(2011073 + 1417500*B4 - 1458000*D3 + 749250*DN + 934245*
        lmMst1 - 1178000*OepS2 + 1350*(620417 + 17670*lmMst1 - 17670*lmMst2)*S2
        + 1767000*T1ep + 3150900*pow2(lmMst1) - 45*lmMst2*(-124139 + 189090*
        lmMst1 + 71550*pow2(lmMst1)) + 4050*(1323 + 1970*lmMst1)*pow2(lmMst2) +
        101250*pow3(lmMst1) - 4860000*pow3(lmMst2)) - 125*Mst2*Mt*s2t*(996211 -
        295488*B4 + 7776*DN - 1030176*lmMst1 - 144*(-1883 + 3618*lmMst1)*lmMst2
        + 98336*OepS2 - 756*(259 + 2634*lmMst1 - 2634*lmMst2)*S2 - 147504*T1ep
        + 49248*pow2(lmMst1) + 5184*(67 + 48*lmMst1)*pow2(lmMst2) - 248832*
        pow3(lmMst2)) + 150*pow2(Mt)*(2199511 - 4320*B4 + 4320*D3 - 2160*DN +
        140160*lmMst1 + 960*(1426 + 303*lmMst1)*lmMst2 - 2880*(-16 + 5*lmMst1 -
        5*lmMst2)*lmMt - 1120*OepS2 + 324*(-3411 + 70*lmMst1 - 70*lmMst2)*S2 +
        1680*T1ep - 2880*(77 + 24*lmMst1)*pow2(lmMst2) + 69120*pow3(lmMst2)))*
        pow6(Mst1) + 1296000*(2 + lmMst2 - 3*pow2(lmMst2))*pow2(s2t)*pow8(Mst2)
        ))/(364500.*Tbeta*pow4(Mst1)*pow6(Mst2))))) - threeLoopFlag*pow2(Al4p)*
        ((Mt*(-8*MuSUSY*z4*(60*Dmsqst2*s2t*pow2(Mst2)*(pow3(Dmsqst2)*(9*pow2(
        Mst1)*pow2(Mst2)*(-2*Mt*MuSUSY*s2t - 6*Tbeta*pow2(Mt) + Tbeta*pow2(
        Mst2)*pow2(s2t)) + (442*Mt*MuSUSY*s2t + 106*Tbeta*pow2(Mt) + 37*Tbeta*
        pow2(Mst2)*pow2(s2t))*pow4(Mst1)) + (18*pow2(Mst1)*pow2(Mst2)*(13*Mt*
        MuSUSY*s2t + 5*Tbeta*pow2(Mt) - 5*Tbeta*pow2(Mst2)*pow2(s2t)) + 2*(221*
        Mt*MuSUSY*s2t + 53*Tbeta*pow2(Mt) - 52*Tbeta*pow2(Mst2)*pow2(s2t))*
        pow4(Mst1) + 27*(2*Mt*MuSUSY*s2t + 2*Tbeta*pow2(Mt) - Tbeta*pow2(Mst2)*
        pow2(s2t))*pow4(Mst2))*pow6(Msq)) + 60*s2t*xDmsqst2*pow2(Dmsqst2)*pow2(
        Mst2)*pow4(Msq)*(3*pow2(Mst1)*pow2(Mst2)*(50*Mt*MuSUSY*s2t + 14*Tbeta*
        pow2(Mt) - 19*Tbeta*pow2(Mst2)*pow2(s2t)) + (442*Mt*MuSUSY*s2t + 106*
        Tbeta*pow2(Mt) - 57*Tbeta*pow2(Mst2)*pow2(s2t))*pow4(Mst1) + 36*Mt*(
        MuSUSY*s2t + Mt*Tbeta)*pow4(Mst2) - 18*Tbeta*pow2(s2t)*pow6(Mst2)) + (-
        ((4*s2t*(2636*Dmglst2*MuSUSY + 35292*Mst2*MuSUSY + 27252*Dmglst2*Mst2*
        Tbeta - 18093*Tbeta*pow2(Mst2))*pow2(Mt) + 2*Mst2*Mt*(45308*Dmglst2*
        MuSUSY - 41169*Mst2*MuSUSY + 15189*Dmglst2*Mst2*Tbeta - 22293*Tbeta*
        pow2(Mst2))*pow2(s2t) + 8*(9279*MuSUSY - 28538*Dmglst2*Tbeta + 9978*
        Mst2*Tbeta)*pow3(Mt) + (-25328*Dmglst2 + 20685*Mst2)*Tbeta*pow3(Mst2)*
        pow3(s2t))*pow4(Mst1)) - 9*pow2(Mst2)*(4*xDmglst2*pow2(Dmglst2)*(2*
        Mst2*s2t*(218*MuSUSY - 99*Mst2*Tbeta)*pow2(Mt) + Mt*(1586*MuSUSY - 327*
        Mst2*Tbeta)*pow2(Mst2)*pow2(s2t) + 4*(81*MuSUSY + 28*Mst2*Tbeta)*pow3(
        Mt) - 793*Tbeta*pow3(s2t)*pow4(Mst2)) + 9*pow2(Mst2)*(2*Mst2*s2t*(170*
        MuSUSY - 101*Mst2*Tbeta)*pow2(Mt) - 5*Mt*(16*MuSUSY + 51*Mst2*Tbeta)*
        pow2(Mst2)*pow2(s2t) + 4*(36*MuSUSY + 103*Mst2*Tbeta)*pow3(Mt) + 40*
        Tbeta*pow3(s2t)*pow4(Mst2)) - 3*Dmglst2*Mst2*(4*Mst2*s2t*(-271*MuSUSY +
        216*Mst2*Tbeta)*pow2(Mt) + Mt*(-848*MuSUSY + 813*Mst2*Tbeta)*pow2(Mst2)
        *pow2(s2t) - 4*(216*MuSUSY + 253*Mst2*Tbeta)*pow3(Mt) + 424*Tbeta*pow3(
        s2t)*pow4(Mst2))) + 3*pow2(Mst1)*(-2*xDmglst2*pow2(Dmglst2)*(4*Mst2*
        s2t*(4045*MuSUSY + 594*Mst2*Tbeta)*pow2(Mt) + Mt*(14978*MuSUSY - 10173*
        Mst2*Tbeta)*pow2(Mst2)*pow2(s2t) + 8*(507*MuSUSY + 2198*Mst2*Tbeta)*
        pow3(Mt) - 2731*Tbeta*pow3(s2t)*pow4(Mst2)) + 3*Dmglst2*Mst2*(8*Mst2*
        s2t*(-709*MuSUSY + 24*Mst2*Tbeta)*pow2(Mt) + 15*Mt*(-296*MuSUSY + 121*
        Mst2*Tbeta)*pow2(Mst2)*pow2(s2t) - 96*(47*MuSUSY - 26*Mst2*Tbeta)*pow3(
        Mt) + 948*Tbeta*pow3(s2t)*pow4(Mst2)) - 3*pow2(Mst2)*(10*Mst2*s2t*(908*
        MuSUSY - 347*Mst2*Tbeta)*pow2(Mt) - Mt*(4552*MuSUSY + 4515*Mst2*Tbeta)*
        pow2(Mst2)*pow2(s2t) + 16*(225*MuSUSY + 326*Mst2*Tbeta)*pow3(Mt) +
        1916*Tbeta*pow3(s2t)*pow4(Mst2))))*pow8(Msq)) - z3*(15*xDmsqst2*pow2(
        Dmsqst2)*pow2(Mst2)*pow4(Msq)*(MuSUSY*(3*Mst2*pow2(Mst1)*(-64*s2t*(
        1512*MuSUSY + 247*Mst2*Tbeta)*pow2(Mt) + 32*Mst2*Mt*(1783*MuSUSY +
        1620*Mst2*Tbeta)*pow2(s2t) + 55296*Tbeta*pow3(Mt) - 12007*Tbeta*pow3(
        Mst2)*pow3(s2t)) + 9*pow3(Mst2)*(-4*s2t*(4032*MuSUSY + 1091*Mst2*Tbeta)
        *pow2(Mt) + 2*Mst2*Mt*(5507*MuSUSY + 6048*Mst2*Tbeta)*pow2(s2t) +
        13824*Tbeta*pow3(Mt) - 5507*Tbeta*pow3(Mst2)*pow3(s2t))) - 16*s2t*(
        3860*MuSUSY*Tbeta*pow2(Mt) + 6*pow2(Mst2)*pow2(s2t)*(263*MuSUSY*Tbeta -
        135*Mst2*(-1 + pow2(Sbeta))*pow2(Sbeta)) + Mt*s2t*(-7776*Mst2*MuSUSY*
        Tbeta + 3542*pow2(MuSUSY) + 243*pow2(Mst2)*(-1 + pow2(Sbeta))*pow2(
        Sbeta)))*pow4(Mst1) - 5184*Dmglst2*(MuSUSY*pow2(Mst1)*(-40*MuSUSY*s2t*
        pow2(Mt) + 8*Mst2*Mt*(2*MuSUSY + 3*Mst2*Tbeta)*pow2(s2t) + 160*Tbeta*
        pow3(Mt) - 5*Tbeta*pow3(Mst2)*pow3(s2t)) + MuSUSY*pow2(Mst2)*(-8*
        MuSUSY*s2t*pow2(Mt) + 6*Mst2*Mt*(MuSUSY + Mst2*Tbeta)*pow2(s2t) + 48*
        Tbeta*pow3(Mt) - 3*Tbeta*pow3(Mst2)*pow3(s2t)) + Mt*pow2(s2t)*(24*
        MuSUSY*Tbeta - 5*Mst2*(-1 + pow2(Sbeta))*pow2(Sbeta))*pow4(Mst1))) +
        pow2(Mst2)*(15*pow4(Dmsqst2)*(9*Mst2*MuSUSY*pow2(Mst1)*(-576*s2t*(56*
        MuSUSY + 5*Mst2*Tbeta)*pow2(Mt) + 32*Mst2*Mt*(1853*MuSUSY + 756*Mst2*
        Tbeta)*pow2(s2t) + 18432*Tbeta*pow3(Mt) - 15335*Tbeta*pow3(Mst2)*pow3(
        s2t)) + 27*MuSUSY*pow3(Mst2)*(-4*s2t*(1728*MuSUSY + 163*Mst2*Tbeta)*
        pow2(Mt) + 2*Mst2*Mt*(4771*MuSUSY + 2592*Mst2*Tbeta)*pow2(s2t) + 7680*
        Tbeta*pow3(Mt) - 4771*Tbeta*pow3(Mst2)*pow3(s2t)) - 16*s2t*(3860*
        MuSUSY*Tbeta*pow2(Mt) + 2*pow2(Mst2)*pow2(s2t)*(3823*MuSUSY*Tbeta -
        1215*Mst2*(-1 + pow2(Sbeta))*pow2(Sbeta)) + Mt*s2t*(-7776*Mst2*MuSUSY*
        Tbeta + 3542*pow2(MuSUSY) + 729*pow2(Mst2)*(-1 + pow2(Sbeta))*pow2(
        Sbeta)))*pow4(Mst1) - 5184*Dmglst2*(MuSUSY*pow2(Mst1)*(-40*MuSUSY*s2t*
        pow2(Mt) + 8*Mst2*Mt*(2*MuSUSY + 3*Mst2*Tbeta)*pow2(s2t) + 160*Tbeta*
        pow3(Mt) - 5*Tbeta*pow3(Mst2)*pow3(s2t)) + MuSUSY*pow2(Mst2)*(-8*
        MuSUSY*s2t*pow2(Mt) + 6*Mst2*Mt*(MuSUSY + Mst2*Tbeta)*pow2(s2t) + 48*
        Tbeta*pow3(Mt) - 3*Tbeta*pow3(Mst2)*pow3(s2t)) + Mt*pow2(s2t)*(24*
        MuSUSY*Tbeta - 5*Mst2*(-1 + pow2(Sbeta))*pow2(Sbeta))*pow4(Mst1))) -
        240*Dmsqst2*(MuSUSY*(27*pow3(Mst2)*(4*s2t*(72*MuSUSY + 29*Mst2*Tbeta)*
        pow2(Mt) - 2*Mst2*Mt*(23*MuSUSY + 108*Mst2*Tbeta)*pow2(s2t) - 192*
        Tbeta*pow3(Mt) + 23*Tbeta*pow3(Mst2)*pow3(s2t)) - 18*Mst2*pow2(Mst1)*(-
        2*s2t*(504*MuSUSY + 101*Mst2*Tbeta)*pow2(Mt) + Mst2*Mt*(-35*MuSUSY +
        432*Mst2*Tbeta)*pow2(s2t) + 576*Tbeta*pow3(Mt) + 52*Tbeta*pow3(Mst2)*
        pow3(s2t)) + 2*s2t*(Mt*s2t*(1771*MuSUSY - 3888*Mst2*Tbeta) + 1930*
        Tbeta*pow2(Mt) - 728*Tbeta*pow2(Mst2)*pow2(s2t))*pow4(Mst1)) + 324*
        Dmglst2*(MuSUSY*pow2(Mst1)*(-40*MuSUSY*s2t*pow2(Mt) + 8*Mst2*Mt*(2*
        MuSUSY + 3*Mst2*Tbeta)*pow2(s2t) + 160*Tbeta*pow3(Mt) - 5*Tbeta*pow3(
        Mst2)*pow3(s2t)) + MuSUSY*pow2(Mst2)*(-8*MuSUSY*s2t*pow2(Mt) + 6*Mst2*
        Mt*(MuSUSY + Mst2*Tbeta)*pow2(s2t) + 48*Tbeta*pow3(Mt) - 3*Tbeta*pow3(
        Mst2)*pow3(s2t)) + Mt*pow2(s2t)*(24*MuSUSY*Tbeta - 5*Mst2*(-1 + pow2(
        Sbeta))*pow2(Sbeta))*pow4(Mst1)))*pow6(Msq)) + 4*MuSUSY*(4*(-2*Mt*s2t*(
        Mst2*s2t*(Dmglst2*(706232*MuSUSY - 546087*Mst2*Tbeta) + 6*Mst2*((91963
        + 162*lmMst1 - 162*lmMst2)*MuSUSY - 4431*Mst2*Tbeta)) + Mt*(1582324*
        Dmglst2*MuSUSY + 298896*Mst2*MuSUSY + 50508*Dmglst2*Mst2*Tbeta + 40629*
        Tbeta*pow2(Mst2))) - 4*(6759*MuSUSY + 834482*Dmglst2*Tbeta - 321216*
        Mst2*Tbeta)*pow3(Mt) + 11*(16354*Dmglst2 + 20937*Mst2)*Tbeta*pow3(Mst2)
        *pow3(s2t))*pow4(Mst1) + 3*pow2(Mst1)*(xDmglst2*pow2(Dmglst2)*(-32*
        Mst2*s2t*(80959*MuSUSY + 22806*Mst2*Tbeta)*pow2(Mt) + 4*Mt*(-458284*
        MuSUSY + 190797*Mst2*Tbeta)*pow2(Mst2)*pow2(s2t) + 8*(154668*MuSUSY +
        342281*Mst2*Tbeta)*pow3(Mt) + 547817*Tbeta*pow3(s2t)*pow4(Mst2)) + 6*
        Mst2*(Dmglst2*(-32*Mst2*s2t*(11864*MuSUSY + 885*Mst2*Tbeta)*pow2(Mt) +
        3*Mt*(-77976*MuSUSY + 49501*Mst2*Tbeta)*pow2(Mst2)*pow2(s2t) + 2400*(5*
        MuSUSY - 67*Mst2*Tbeta)*pow3(Mt) + 64119*Tbeta*pow3(s2t)*pow4(Mst2)) +
        Mst2*(2*Mst2*s2t*(-58544*MuSUSY + 1859*Mst2*Tbeta)*pow2(Mt) + Mt*(-4*(
        35719 + 108*lmMst1 - 108*lmMst2)*MuSUSY + 34131*Mst2*Tbeta)*pow2(Mst2)*
        pow2(s2t) - 16*(981*MuSUSY - 8219*Mst2*Tbeta)*pow3(Mt) + 4*(11924 + 27*
        lmMst1 - 27*lmMst2)*Tbeta*pow3(s2t)*pow4(Mst2)))) + 9*pow2(Mst2)*(
        xDmglst2*pow2(Dmglst2)*(-2*Mst2*Mt*s2t*(Mst2*s2t*(122917*MuSUSY -
        196638*Mst2*Tbeta) + Mt*(262184*MuSUSY + 166791*Mst2*Tbeta)) + 8*(
        21888*MuSUSY + 32323*Mst2*Tbeta)*pow3(Mt) + 122917*Tbeta*pow3(s2t)*
        pow4(Mst2)) + 6*Mst2*(Dmglst2*(-4*Mst2*s2t*(15137*MuSUSY + 878*Mst2*
        Tbeta)*pow2(Mt) + Mt*(-35230*MuSUSY + 45411*Mst2*Tbeta)*pow2(Mst2)*
        pow2(s2t) + 4*(664*MuSUSY - 9*Mst2*Tbeta)*pow3(Mt) + 17615*Tbeta*pow3(
        s2t)*pow4(Mst2)) + Mst2*(-2*Mst2*s2t*(11930*MuSUSY + 9*(-77 + 8*lmMst1
        - 8*lmMst2)*Mst2*Tbeta)*pow2(Mt) + 3*Mt*(-4*(1319 + 6*lmMst1 - 6*
        lmMst2)*MuSUSY + 5965*Mst2*Tbeta)*pow2(Mst2)*pow2(s2t) - 4*(920*MuSUSY
        - 2883*Mst2*Tbeta)*pow3(Mt) + 6*(1319 + 6*lmMst1 - 6*lmMst2)*Tbeta*
        pow3(s2t)*pow4(Mst2)))))*pow8(Msq))))/(11664.*Tbeta*pow6(Mst2)*pow8(
        Msq)) + MuSUSY*(-(pow2(Mt)*pow2(s2t)*((32*(1 + lmMst2)*(-Dmglst2 + 3*
        Dmglst2*lmMst2 + Mst2 + lmMst2*Mst2)*pow2(Mst2))/(3.*pow2(Mst1)) + ((2*
        OepS2*(Dmglst2*(8163*pow2(Mst1)*pow2(Mst2) + 23734*pow4(Mst1) + 945*
        pow4(Mst2)) - 3*(627*pow2(Mst1)*pow3(Mst2) + 1066*Mst2*pow4(Mst1) +
        189*pow5(Mst2))))/729. - (S2*(Dmglst2*(9*(23989 + 27210*lmMst1 - 27210*
        lmMst2)*pow2(Mst1)*pow2(Mst2) + 2*(525961 + 356010*lmMst1 - 356010*
        lmMst2)*pow4(Mst1) + 81*(-453 + 350*lmMst1 - 350*lmMst2)*pow4(Mst2)) -
        15*(9*(685 + 418*lmMst1 - 418*lmMst2)*pow2(Mst1)*pow3(Mst2) + 2*(6143 +
        3198*lmMst1 - 3198*lmMst2)*Mst2*pow4(Mst1) + 81*(-1 + 14*lmMst1 - 14*
        lmMst2)*pow5(Mst2))))/540.)/pow4(Mst2) + (pow2(Mst1)*(Mst2*(
        267.2878086419753 + 48*B4 + (571*lmMst2)/3. + (4*(-7 + 8*lmMst2)*pow2(
        lmMst1))/3. + (748*pow2(lmMst2))/3. - (lmMst1*(495 + 720*lmMst2 + 292*
        pow2(lmMst2)))/3. + (80*Dmsqst2*(lmMst1 - lmMst2))/pow2(Msq) + (260*
        pow3(lmMst2))/3. + (20*(3 + 16*lmMst1 - 16*lmMst2)*pow4(Dmsqst2))/pow8(
        Msq)) + Dmglst2*(30.137808641975308 + (296*B4)/3. - (4*DN)/3. - (985*
        lmMst1)/9. + (28*pow2(lmMst1))/3. + (lmMst2*(2149 - 1296*lmMst1 + 96*
        pow2(lmMst1)))/9. + (134.66666666666666 - 140*lmMst1)*pow2(lmMst2) + (
        388*pow3(lmMst2))/3. + (400*Dmsqst2*(1 - lmMst1 + lmMst2)*(pow3(
        Dmsqst2) + pow6(Msq)))/pow8(Msq))))/pow2(Mst2) + Mst2*(
        188.52083333333334 + (110*B4)/3. - DN/3. - (124*lmMst1)/3. - (16*pow2(
        lmMst1))/3. - (4*lmMst2*(-198 + 73*lmMst1 + 4*pow2(lmMst1)))/3. + (168
        - 38*lmMst1)*pow2(lmMst2) + (40*Dmsqst2*(1 + lmMst1 - lmMst2))/pow2(
        Msq) + (130*pow3(lmMst2))/3. + (20*(11 + 8*lmMst1 - 8*lmMst2)*pow4(
        Dmsqst2))/pow8(Msq)) + (pow4(Mst1)*(409.2597736625514 + 48*B4 - (1609*
        lmMst1)/9. + (326*pow2(lmMst1))/9. + (2*lmMst2*(845 - 1535*lmMst1 +
        264*pow2(lmMst1)))/9. + (304.8888888888889 - 188*lmMst1)*pow2(lmMst2) -
        (16*pow3(lmMst1))/9. + (1180*pow3(lmMst2))/9. + (10*Dmsqst2*(-4 + 9*
        lmMst1 - 9*lmMst2)*(pow3(Dmsqst2) + pow6(Msq)))/pow8(Msq)))/pow3(Mst2)
        + Dmglst2*(21.73935185185185 + (262*B4)/3. - (5*DN)/3. - (52*lmMst1)/3.
         + (16*pow2(lmMst1))/3. - (4*lmMst2*(-221 + 55*lmMst1 + 4*pow2(lmMst1))
        )/3. + (2*(232 - 121*lmMst1)*pow2(lmMst2))/3. + 86*pow3(lmMst2) - (
        pow4(Mst1)*(516.797085048011 - (296*B4)/3. + (4*DN)/3. + (17570*lmMst1)
        /81. + (46*pow2(lmMst1))/3. - (lmMst2*(28667 - 5238*lmMst1 + 3024*pow2(
        lmMst1)))/81. + (4*(-60 + 157*lmMst1)*pow2(lmMst2))/3. - (16*pow3(
        lmMst1))/3. - (500*pow3(lmMst2))/3. + (10*Dmsqst2*(-80 + 43*lmMst1 -
        43*lmMst2)*(pow3(Dmsqst2) + pow6(Msq)))/pow8(Msq)))/pow4(Mst2) - (40*
        Dmsqst2*(2 + 5*lmMst1 - 5*lmMst2)*(pow3(Dmsqst2) + pow6(Msq)))/pow8(
        Msq)))) + ((2*(-1 + 2*lmMst2)*s2t*shiftst3*pow2(Mst2)*pow3(Mt))/3. + (
        5*Mt*shiftst1*pow3(s2t)*pow4(Mst2)*(1 - 2*lmMst2 + (Dmsqst2*(3 - 2*
        lmMst2)*(pow3(Dmsqst2) + pow6(Msq)))/pow8(Msq)))/3.)/pow2(Mst1) - (Mt*
        T1ep*(Dmglst2*(27*pow2(Mst1)*pow2(Mst2)*(-800*Mst2*s2t*pow2(Mt) - 907*
        Mt*pow2(Mst2)*pow2(s2t) + 1984*pow3(Mt) + 316*pow3(Mst2)*pow3(s2t)) +
        2*(-66168*Mst2*s2t*pow2(Mt) - 35601*Mt*pow2(Mst2)*pow2(s2t) + 129704*
        pow3(Mt) + 12664*pow3(Mst2)*pow3(s2t))*pow4(Mst1) + 189*(20*pow3(Mt)*
        pow4(Mst2) - 15*Mt*pow2(s2t)*pow6(Mst2) + 4*pow3(s2t)*pow7(Mst2)))*
        pow8(Msq) + 3*Mst2*(20*Mst2*s2t*pow2(Mst1)*(9*pow2(Mst2)*(-6*pow2(Mt) +
        pow2(Mst2)*pow2(s2t)) + pow2(Mst1)*(106*pow2(Mt) + 37*pow2(Mst2)*pow2(
        s2t)))*pow4(Dmsqst2) - 20*Dmsqst2*Mst2*s2t*(-2*(53*pow2(Mt) - 52*pow2(
        Mst2)*pow2(s2t))*pow4(Mst1) + 27*(-2*pow2(Mt) + pow2(Mst2)*pow2(s2t))*
        pow4(Mst2) + 90*pow2(Mst1)*(-(pow2(Mst2)*pow2(Mt)) + pow2(s2t)*pow4(
        Mst2)))*pow6(Msq) + (-3*pow2(Mst1)*pow2(Mst2)*(-3470*Mst2*s2t*pow2(Mt)
        - 627*Mt*pow2(Mst2)*pow2(s2t) + 1760*pow3(Mt) + 1700*pow3(Mst2)*pow3(
        s2t)) + (24124*Mst2*s2t*pow2(Mt) + 3198*Mt*pow2(Mst2)*pow2(s2t) -
        16240*pow3(Mt) - 6895*pow3(Mst2)*pow3(s2t))*pow4(Mst1) - 27*(-106*Mst2*
        s2t*pow2(Mt) - 21*Mt*pow2(Mst2)*pow2(s2t) + 28*pow3(Mt) + 46*pow3(Mst2)
        *pow3(s2t))*pow4(Mst2))*pow8(Msq))))/(729.*pow6(Mst2)*pow8(Msq)) + (
        pow4(Mt)*(pow2(Mst1)*(-12*Mst2*(24*(2200*OepS2 - 81*(457 + 550*lmMst1 -
        550*lmMst2)*S2)*pow2(Mst1)*pow2(Mst2) + 40*(4060*OepS2 - 27*(4082 +
        3045*lmMst1 - 3045*lmMst2)*S2)*pow4(Mst1) + 27*(280*OepS2 - 81*(-79 +
        70*lmMst1 - 70*lmMst2)*S2)*pow4(Mst2)) + 4*Dmglst2*(432*(1240*OepS2 -
        9*(601 + 2790*lmMst1 - 2790*lmMst2)*S2)*pow2(Mst1)*pow2(Mst2) + 40*(
        64852*OepS2 - 27*(46187 + 48639*lmMst1 - 48639*lmMst2)*S2)*pow4(Mst1) +
        675*(56*OepS2 - 81*(-27 + 14*lmMst1 - 14*lmMst2)*S2)*pow4(Mst2)))*pow8(
        Msq) + 60*Mst2*(678923 + 19440*lmMt + 540*(661 - 30*lmMt)*pow2(lmMst2)
        - 15552*pow2(lmMt) - 108*((457 - 12*lmMst2 - 126*lmMt)*pow2(lmMst1) +
        lmMst1*(1841 + lmMst2*(2626 - 24*lmMt) - 420*lmMt + 840*pow2(lmMst2) -
        288*pow2(lmMt)) + lmMst2*(619 + 498*lmMt + 288*pow2(lmMt))) + 216*pow3(
        lmMst1) + 89208*pow3(lmMst2))*pow6(Mst1)*pow8(Msq) + 622080*pow2(1 +
        lmMst2)*pow7(Mst2)*pow8(Msq) + 243*pow2(Mst1)*pow5(Mst2)*(1200*(8 + 32*
        lmMst1 - 65*lmMst2 + 33*lmMt)*pow4(Dmsqst2) + 4800*Dmsqst2*(5 + 2*
        lmMst1 - 5*lmMst2 + 3*lmMt)*pow6(Msq) + (18813 + 480*lmMst2*(3 - 11*
        lmMt) - 8000*lmMt - 1280*(1 + lmMst2)*pow2(lmMst1) + 160*(163 - 16*
        lmMt)*pow2(lmMst2) - 160*lmMst1*(46 + 2*lmMst2*(57 - 8*lmMt) - 16*lmMt
        + 41*pow2(lmMst2)) - 3840*pow2(lmMt) + 7840*pow3(lmMst2))*pow8(Msq)) +
        72*pow2(Mst2)*pow4(Mst1)*(1800*Dmglst2*Dmsqst2*(497 - 306*lmMst1 + 609*
        lmMst2 - 303*lmMt)*(pow3(Dmsqst2) + pow6(Msq)) + 16200*Dmsqst2*(1 + 6*
        lmMst1 - 13*lmMst2 + 7*lmMt)*Mst2*(pow3(Dmsqst2) + pow6(Msq)) - 3*
        Dmglst2*(97837 + 71580*lmMt - 360*(23 + 8*lmMst2 - 4*lmMt)*pow2(lmMst1)
        + 720*(97 + 2*lmMt)*pow2(lmMst2) - 8640*pow2(lmMt) - 30*lmMst2*(-2911 +
        396*lmMt + 144*pow2(lmMt)) + 10*lmMst1*(-6157 + 642*lmMt - 6*lmMst2*(
        791 + 48*lmMt) + 882*pow2(lmMst2) + 432*pow2(lmMt)) - 5940*pow3(lmMst2)
        )*pow8(Msq) + 2*Mst2*(110219 - 1080*lmMt + 540*(-15 + 4*lmMt)*pow2(
        lmMst1) + 810*(121 - 8*lmMt)*pow2(lmMst2) - 135*lmMst1*(337 + lmMst2*(
        588 - 32*lmMt) - 132*lmMt + 194*pow2(lmMst2) - 48*pow2(lmMt)) - 6480*
        pow2(lmMt) - 405*lmMst2*(31 + 54*lmMt + 16*pow2(lmMt)) + 26190*pow3(
        lmMst2))*pow8(Msq)) - 5*Dmglst2*(4*(5659217 + 1592460*lmMt - 972*(569 +
        180*lmMst2 - 126*lmMt)*pow2(lmMst1) + 324*(5353 + 126*lmMt)*pow2(
        lmMst2) + 36*lmMst2*(39031 + 3204*lmMt - 5184*pow2(lmMt)) - 72*lmMst1*(
        27653 + 3015*lmMt + 18*lmMst2*(689 + 126*lmMt) - 2160*pow2(lmMst2) -
        2592*pow2(lmMt)) - 186624*pow2(lmMt) + 1944*pow3(lmMst1) + 17496*pow3(
        lmMst2))*pow6(Mst1)*pow8(Msq) - 124416*(-1 + 2*lmMst2 + 3*pow2(lmMst2))
        *pow6(Mst2)*pow8(Msq) - 27*pow2(Mst1)*pow4(Mst2)*(960*Dmsqst2*(101 -
        90*lmMst1 + 177*lmMst2 - 87*lmMt)*(pow3(Dmsqst2) + pow6(Msq)) - (57495
        + lmMst2*(55616 - 4704*lmMt) + 45408*lmMt + 2304*(-1 + lmMst2)*pow2(
        lmMst1) + 96*(79 + 48*lmMt)*pow2(lmMst2) + 288*lmMst1*(-22 + 16*lmMt -
        2*lmMst2*(9 + 8*lmMt) + 41*pow2(lmMst2)) - 14112*pow3(lmMst2))*pow8(
        Msq)))))/(43740.*pow2(Mst1)*pow6(Mst2)*pow8(Msq)) - s2t*pow3(Mt)*(
        302.29104938271604 + (8*D3)/9. - (8*DN)/9. + (614*lmMst1)/27. + (32*
        lmMt)/3. + (124*pow2(lmMst1))/9. - (2*lmMst2*(-1901 + 6*lmMst1 + 117*
        pow2(lmMst1)))/27. + (46 + (32*lmMst1)/9.)*pow2(lmMst2) + (5*Dmsqst2*(
        367 + 32*lmMst1 - 80*lmMst2))/(18.*pow2(Msq)) + (32*pow3(lmMst1))/9. +
        (14*pow3(lmMst2))/9. + (32*(1 + lmMst2)*(4*Dmglst2*lmMst2 + Mst2 +
        lmMst2*Mst2)*pow3(Mst2))/(9.*pow4(Mst1)) + (pow2(Mst1)*(Mst2*(
        399.70458847736626 - (84094*lmMst1)/405. + (32*(1 - 2*lmMst1 + 2*
        lmMst2)*lmMt)/3. - (1366*pow2(lmMst1))/27. - (2*lmMst2*(-82997 + 4800*
        lmMst1 + 1440*pow2(lmMst1)))/405. + (2*(1579 - 426*lmMst1)*pow2(lmMst2)
        )/27. + (Dmsqst2*(83.94493827160494 - (356*lmMst1)/15. + (4*(89 + 60*
        lmMst1)*lmMst2)/15. - 8*(pow2(lmMst1) + pow2(lmMst2))))/pow2(Msq) - (
        164*pow3(lmMst1))/27. + (1208*pow3(lmMst2))/27. + ((176.66444444444446
         + (106*lmMst1)/
        45. - (2*(53 + 270*lmMst1)*lmMst2)/45. + 6*(pow2(lmMst1) + pow2(lmMst2)
        ))*pow4(Dmsqst2))/pow8(Msq)) + (4*Dmglst2*(3856958 - 27000*B4 + 27000*
        D3 - 13500*DN + 1270770*lmMst1 + 162000*(-2 + lmMst1 - lmMst2)*lmMt +
        900*pow2(lmMst1) - 30*lmMst2*(-153166 + 17835*lmMst1 + 3375*pow2(
        lmMst1)) + 450*(2627 - 1695*lmMst1)*pow2(lmMst2) - 2250*pow3(lmMst1) +
        866250*pow3(lmMst2) + (50625*Dmsqst2*(-65 + 2*lmMst1 - 2*lmMst2)*(pow3(
        Dmsqst2) + pow6(Msq)))/pow8(Msq)))/30375.))/pow3(Mst2) + (pow4(Mst1)*(
        Dmglst2*(311.41128771790903 - (32*B4)/9. + (32*D3)/9. - (16*DN)/9. + (
        14926634*lmMst1)/19845. + (19352*pow2(lmMst1))/189. + (2*lmMst2*(
        2917823 - 3813600*lmMst1 + 97020*pow2(lmMst1)))/19845. + (8*(8677 -
        5439*lmMst1)*pow2(lmMst2))/189. - (64*lmMt*(4 + 6*lmMst2 - 2*lmMst1*(3
        + lmMst2) + pow2(lmMst1) + pow2(lmMst2)))/3. - (8*pow3(lmMst1))/9. + (
        664*pow3(lmMst2))/3.) + Mst2*(599.6859861506036 - (99165049*lmMst1)/
        198450. + lmMst2*(713.9201259763164 - (57994*lmMst1)/945. - (350*pow2(
        lmMst1))/9.) - (90913*pow2(lmMst1))/945. + (200.24021164021164 - (286*
        lmMst1)/9.)*pow2(lmMst2) + (32*lmMt*(1 + 4*lmMst2 - 2*lmMst1*(2 +
        lmMst2) + pow2(lmMst1) + pow2(lmMst2)))/3. + (26*pow3(lmMst1))/27. + (
        1882*pow3(lmMst2))/27. + (Dmsqst2*(87.53310385647498 - (941812*lmMst1)/
        19845. + (47.45840262030738 + (484*lmMst1)/21.)*lmMst2 - (242*(pow2(
        lmMst1) + pow2(lmMst2)))/21.)*(pow3(Dmsqst2) + pow6(Msq)))/pow8(Msq))))
        /pow5(Mst2) - S2*(42.3 - 159*lmMst1 + 159*lmMst2 + (2532*Dmglst2)/(5.*
        Mst2) - (30*Dmsqst2*(9 + 2*lmMst1 - 2*lmMst2))/pow2(Msq) + (pow2(Mst1)*
        (Dmglst2*(556 + 400*lmMst1 - 400*lmMst2) - Mst2*(570.4333333333333 + (
        1735*lmMst1)/3. - (1735*lmMst2)/3. + (10*Dmsqst2*(127 + 30*lmMst1 - 30*
        lmMst2))/(3.*pow2(Msq)) + ((156.66666666666666 - 60*lmMst1 + 60*lmMst2)
        *pow4(Dmsqst2))/pow8(Msq))))/pow3(Mst2) + (2*pow4(Mst1)*(12*Dmglst2*(
        17269 + 13785*lmMst1 - 13785*lmMst2) - Mst2*(138286 + 90465*lmMst1 -
        90465*lmMst2 + (25*Dmsqst2*(1283 + 318*lmMst1 - 318*lmMst2)*(pow3(
        Dmsqst2) + pow6(Msq)))/pow8(Msq))))/(135.*pow5(Mst2)) - (240*pow4(
        Dmsqst2))/pow8(Msq)) - (5*(-1453 + 48*lmMst1 + 96*lmMst2)*pow4(Dmsqst2)
        )/(54.*pow8(Msq)) - (2*pow2(Mst2)*(-30*(5 + 2*lmMst2)*pow4(Dmsqst2) +
        30*Dmsqst2*(1 - 2*lmMst2)*pow6(Msq) + (71 + 122*lmMst2 + 32*lmMst1*(1 +
        lmMst2) + 59*pow2(lmMst2))*pow8(Msq)))/(9.*pow2(Mst1)*pow8(Msq)) + ((-
        4*Dmglst2*(150*Dmsqst2*(41*pow2(Mst1) + 18*pow2(Mst2))*(pow3(Dmsqst2) +
        pow6(Msq)) - (-15*(-52 + 2*(51 + 16*lmMst1)*lmMst2 + 59*pow2(lmMst2))*
        pow2(Mst2) + pow2(Mst1)*(12454 - 120*B4 + 120*D3 - 60*DN + 690*lmMst1 +
        345*pow2(lmMst1) - 5*lmMst2*(-3121 + 42*lmMst1 + 96*pow2(lmMst1)) + (
        3630 - 480*lmMst1)*pow2(lmMst2) + 960*pow3(lmMst2)))*pow8(Msq)))/(135.*
        Mst2*pow2(Mst1)) + (OepS2*(96*Dmglst2*pow2(Mst1)*(919*pow2(Mst1) + 150*
        pow2(Mst2))*pow8(Msq) - 4*Mst2*(20*pow4(Dmsqst2)*(-27*pow2(Mst1)*pow2(
        Mst2) + 53*pow4(Mst1)) + 20*Dmsqst2*(45*pow2(Mst1)*pow2(Mst2) + 53*
        pow4(Mst1) + 27*pow4(Mst2))*pow6(Msq) + (5205*pow2(Mst1)*pow2(Mst2) +
        12062*pow4(Mst1) + 1431*pow4(Mst2))*pow8(Msq))))/(729.*pow5(Mst2)))/
        pow8(Msq)) + Mt*((-5*s2t*shiftst1*(4*pow2(Mst2)*pow2(Mt) - 2*pow2(Mst1)
        *(2*pow2(Mt) + (-lmMst1 + lmMst2)*pow2(Mst2)*pow2(s2t)) + pow2(s2t)*
        pow4(Mst1))*(Dmsqst2*(3 - 2*lmMst2)*(pow3(Dmsqst2) + pow6(Msq)) + (1 -
        2*lmMst2)*pow8(Msq)))/(3.*pow2(Mst1)*pow8(Msq)) + pow3(s2t)*(-((-1 + 2*
        lmMst2)*shiftst3*pow2(Mst2)*((-1 - lmMst1 + lmMst2)*pow2(Mst1) + pow2(
        Mst2)))/(6.*pow2(Mst1)) + (8*(1 + lmMst2)*(4*Dmglst2*lmMst2 + Mst2 +
        lmMst2*Mst2)*pow5(Mst2))/(9.*pow4(Mst1)) - (pow4(Mst1)*(-(Mst2*(
        96.95148419197191 - (61388401*lmMst1)/793800. + lmMst2*(
        73.69595742000504 + (302047*lmMst1)/1890. - (257*pow2(lmMst1))/9.) - (
        76483*pow2(lmMst1))/945. - (78.87883597883598 - 61*lmMst1)*pow2(lmMst2)
        + (Dmsqst2*(17.77343371446168 - (452768*lmMst1)/19845. + (
        22.815217939027463 + (122*lmMst1)/7.)*lmMst2 - (61*(pow2(lmMst1) +
        pow2(lmMst2)))/7.))/pow2(Msq) - (35*pow3(lmMst1))/27. - (841*pow3(
        lmMst2))/27. - ((55.45597659639988 + (27119*lmMst1)/11340. + (-
        2.3914462081128747 + 16*lmMst1)*lmMst2 - 8*(pow2(lmMst1) + pow2(lmMst2)
        ))*pow4(Dmsqst2))/pow8(Msq))) + Dmglst2*(203.8689796251638 - (2740559*
        lmMst1)/119070. + (433*pow2(lmMst1))/189. + lmMst2*(18.238590744939952
         - (2420*lmMst1)/
        189. + (104*pow2(lmMst1))/3.) + (10.513227513227513 - (568*lmMst1)/9.)*
        pow2(lmMst2) - (56*pow3(lmMst1))/27. + (824*pow3(lmMst2))/27. + (5*
        Dmsqst2*(-419 + 57*lmMst1 - 57*lmMst2)*(pow3(Dmsqst2) + pow6(Msq)))/(
        27.*pow8(Msq)))))/pow3(Mst2) - Dmglst2*((pow2(Mst1)*(5430043 - 3089580*
        lmMst1 - 773100*pow2(lmMst1) + 60*lmMst2*(8743 - 24630*lmMst1 + 55350*
        pow2(lmMst1)) + 1800*(808 - 3765*lmMst1)*pow2(lmMst2) + 45000*pow3(
        lmMst1) + 3411000*pow3(lmMst2) + (202500*Dmsqst2*(-55 + 34*lmMst1 - 34*
        lmMst2)*(pow3(Dmsqst2) + pow6(Msq)))/pow8(Msq)))/(121500.*Mst2) + Mst2*
        (46.25617283950617 - (4*B4)/3. + (16*D3)/9. - (10*DN)/9. + (44*lmMst1)/
        3. - (13*pow2(lmMst1))/3. + (4*lmMst2*(-81 - 124*lmMst1 + 8*pow2(
        lmMst1)))/9. + (48.77777777777778 - (82*lmMst1)/3.)*pow2(lmMst2) + (
        214*pow3(lmMst2))/9. + (20*Dmsqst2*(5 + 2*lmMst1 - 2*lmMst2)*(pow3(
        Dmsqst2) + pow6(Msq)))/(3.*pow8(Msq)))) + pow2(Mst1)*(190.4424279835391
         + 2*B4 + (2*D3)/9. - (2*DN)/9. + (22004*lmMst1)/405. - (736*pow2(
        lmMst1))/27. - (lmMst2*(12148 - 76560*lmMst1 + 12105*pow2(lmMst1)))/
        810. + (5*(-583 + 438*lmMst1)*pow2(lmMst2))/54. + (Dmsqst2*(97294 +
        10485*lmMst1 + 45*(-383 + 330*lmMst1)*lmMst2 - 7425*(pow2(lmMst1) +
        pow2(lmMst2))))/(2025.*pow2(Msq)) - (55*pow3(lmMst1))/27. - (1273*pow3(
        lmMst2))/54. - ((36.69808641975309 - (2453*lmMst1)/90. + ((2753 + 420*
        lmMst1)*lmMst2)/90. - (7*(pow2(lmMst1) + pow2(lmMst2)))/3.)*pow4(
        Dmsqst2))/pow8(Msq)) + (pow2(Mst2)*(15*(3847 + lmMst1*(744 - 288*
        lmMst2) - 744*lmMst2 + 288*pow2(lmMst2))*pow4(Dmsqst2) + 30*Dmsqst2*(
        1177 + 48*lmMst1*(7 - 3*lmMst2) - 336*lmMst2 + 144*pow2(lmMst2))*pow6(
        Msq) + (21005 + 1440*B4 - 144*D3 + 72*DN + 21216*lmMst1 + 1908*pow2(
        lmMst1) - 12*lmMst2*(2950 - 2040*lmMst1 + 117*pow2(lmMst1)) + 108*(-283
        + 98*lmMst1)*pow2(lmMst2) + 576*pow3(lmMst1) - 9756*pow3(lmMst2))*pow8(
        Msq)))/(648.*pow8(Msq)) - (pow3(Mst2)*(Mst2*(-30*(5 + 2*lmMst2)*pow4(
        Dmsqst2) + 30*Dmsqst2*(1 - 2*lmMst2)*pow6(Msq) + (119 + 218*lmMst2 +
        32*lmMst1*(1 + lmMst2) + 107*pow2(lmMst2))*pow8(Msq)) + Dmglst2*(360*
        pow4(Dmsqst2) + 360*Dmsqst2*pow6(Msq) + 2*(-20 + 230*lmMst2 + 32*
        lmMst1*lmMst2 + 155*pow2(lmMst2))*pow8(Msq))))/(18.*pow2(Mst1)*pow8(
        Msq)) + (2*Dmglst2*(27*(632*OepS2 + 9*(16193 - 1422*lmMst1 + 1422*
        lmMst2)*S2)*pow2(Mst1)*pow2(Mst2) + 2*(25328*OepS2 + 27*(47051 - 18996*
        lmMst1 + 18996*lmMst2)*S2)*pow4(Mst1) + 27*(56*OepS2 - 81*(-1677 + 14*
        lmMst1 - 14*lmMst2)*S2)*pow4(Mst2))*pow8(Msq) - 3*Mst2*(-10*pow4(
        Dmsqst2)*(9*(8*OepS2 - 27*(391 + 6*lmMst1 - 6*lmMst2)*S2)*pow2(Mst1)*
        pow2(Mst2) + (296*OepS2 - 27*(3871 + 222*lmMst1 - 222*lmMst2)*S2)*pow4(
        Mst1) - 126846*S2*pow4(Mst2)) + 10*Dmsqst2*(18*(40*OepS2 - 27*(59 + 30*
        lmMst1 - 30*lmMst2)*S2)*pow2(Mst1)*pow2(Mst2) + 4*(208*OepS2 - 27*(347
        + 156*lmMst1 - 156*lmMst2)*S2)*pow4(Mst1) + 27*(8*OepS2 - 81*(-15 + 2*
        lmMst1 - 2*lmMst2)*S2)*pow4(Mst2))*pow6(Msq) + (60*(340*OepS2 - 81*(424
        + 85*lmMst1 - 85*lmMst2)*S2)*pow2(Mst1)*pow2(Mst2) + 35*(788*OepS2 -
        27*(2662 + 591*lmMst1 - 591*lmMst2)*S2)*pow4(Mst1) + 27*(184*OepS2 -
        81*(307 + 46*lmMst1 - 46*lmMst2)*S2)*pow4(Mst2))*pow8(Msq)))/(4374.*
        pow3(Mst2)*pow8(Msq))))) + ((xDmsqst2*pow2(Dmsqst2)*(10*Mt*s2t*(s2t*(
        694575*Mst2*pow2(Sbeta)*(-6*(11 + 39*lmMst1 - 39*lmMst2)*s2t*pow2(Mst2)
        *(-1 + pow2(Sbeta)) + (Dmglst2*(4324 - 660*lmMst1 + 660*lmMst2) + 3*(-
        167 + 5*lmMst1 - 5*lmMst2)*Mst2)*Mt*pow2(Sbeta))*pow6(Mst1) + Mt*(1029*
        Mst2*pow2(MuSUSY)*(75*(1728*Dmglst2*(4 + lmMst1 - lmMst2) + Mst2*(-8771
        + 128*OepS2 + 72*lmMst1*(-29 + 12*lmMst2 - 36*S2) + 57024*S2 + 72*
        lmMst2*(23 + 36*S2) - 864*pow2(lmMst2)))*pow2(Mst1)*pow2(Mst2) - 4*(
        4050*Dmglst2*(23 - 42*lmMst1 + 42*lmMst2) + Mst2*(212566 - 69615*lmMst2
        - 10000*OepS2 - 1350*(947 + 150*lmMst2)*S2 + 45*lmMst1*(1547 - 180*
        lmMst2 + 4500*S2) - 4050*pow2(lmMst1) + 12150*pow2(lmMst2)))*pow4(Mst1)
        - 16200*(-12*Dmglst2 + Mst2 + 2*lmMst2*Mst2)*pow4(Mst2)) - (2*(
        593331163 - 60642400*OepS2 + 1260*lmMst2*(8263 - 974610*S2) +
        1143733500*S2 + 1260*lmMst1*(-8263 + 71820*lmMst2 + 974610*S2) -
        61916400*pow2(lmMst1) - 28576800*pow2(lmMst2))*pow2(MuSUSY) + 694575*
        Mst2*(Dmglst2*(4324 - 660*lmMst1 + 660*lmMst2) + 3*(-167 + 5*lmMst1 -
        5*lmMst2)*Mst2)*pow2(Sbeta))*pow6(Mst1))) + 411600*Mt*pow2(MuSUSY)*(
        648*Dmglst2*Mt*pow2(Mst1)*((8 - 15*lmMst1 + 15*lmMst2)*pow2(Mst1) + (-2
        - 5*lmMst1 + 5*lmMst2)*pow2(Mst2)) + 9*(36*(5 + 4*lmMst1 - 4*lmMst2)*Mt
        - Mst2*s2t*(9*(-1 + 2*lmMst1 - 2*lmMst2)*(-2 + lmMst2)*shiftst1 + 4*
        T1ep))*pow2(Mst1)*pow3(Mst2) + 6*Mst2*(108*(1 + 3*lmMst1 - 3*lmMst2)*Mt
        + Mst2*s2t*(-27*(lmMst1 - lmMst2)*(-2 + lmMst2)*shiftst1 - 25*T1ep))*
        pow4(Mst1) + s2t*(-((81*(lmMst1 - lmMst2)*(-3 + 2*lmMst2)*shiftst1 +
        442*T1ep)*pow6(Mst1)) + 81*(-2 + lmMst2)*shiftst1*pow6(Mst2)))) +
        MuSUSY*Tbeta*(-4116000*T1ep*pow2(Mst1)*(-3*Mt*pow2(Mst2)*pow3(s2t)*(19*
        pow2(Mst1)*(pow2(Mst1) + pow2(Mst2)) + 6*pow4(Mst2)) + 2*s2t*pow3(Mt)*(
        21*pow2(Mst1)*pow2(Mst2) + 53*pow4(Mst1) + 18*pow4(Mst2))) + 37044000*
        pow2(Mst1)*(4*Dmglst2*(497 - 306*lmMst1 + 609*lmMst2 - 303*lmMt)*pow2(
        Mst1) + 36*(1 + 6*lmMst1 - 13*lmMst2 + 7*lmMt)*Mst2*pow2(Mst1) + 4*
        Dmglst2*(101 - 90*lmMst1 + 177*lmMst2 - 87*lmMt)*pow2(Mst2) + 9*(16 +
        16*lmMst1 - 35*lmMst2 + 19*lmMt)*pow3(Mst2))*pow4(Mt) + 500094000*pow2(
        Mst1)*pow2(Mt)*pow2(s2t)*(-2*(1 + 8*lmMst1 - 8*lmMst2)*pow2(Mst1)*pow3(
        Mst2) + (4 - 9*lmMst1 + 9*lmMst2)*Mst2*pow4(Mst1) + Dmglst2*(40*(-1 +
        lmMst1 - lmMst2)*pow2(Mst1)*pow2(Mst2) + (-80 + 43*lmMst1 - 43*lmMst2)*
        pow4(Mst1) + 4*(2 + 5*lmMst1 - 5*lmMst2)*pow4(Mst2)) - 2*(5 + 4*lmMst1
        - 4*lmMst2)*pow5(Mst2)) + 4*s2t*pow3(Mt)*(385875*(5904*Dmglst2 - Mst2*(
        3655 - 64*OepS2 + 8424*S2 + 144*lmMst1*(1 + 9*S2) - 144*lmMst2*(4 + 9*
        S2)))*pow2(Mst1)*pow3(Mst2) + 5145*Mst2*(16200*Dmglst2*(65 - 2*lmMst1 +
        2*lmMst2) - Mst2*(279089 - 5600*OepS2 + 1260*lmMst2*(29 - 90*S2) +
        812700*S2 + 180*lmMst1*(-203 + 90*lmMst2 + 630*S2) - 8100*(pow2(lmMst1)
        + pow2(lmMst2))))*pow4(Mst1) - 83349000*(-12*Dmglst2 + Mst2 + 2*lmMst2*
        Mst2)*pow5(Mst2) - (1094369501 - 72716000*OepS2 + 2520*lmMst2*(235453 -
        584325*S2) + 5940931500*S2 + 2520*lmMst1*(-235453 + 114345*lmMst2 +
        584325*S2) - 144074700*(pow2(lmMst1) + pow2(lmMst2)))*pow6(Mst1)) + Mt*
        (-166698000*(2 - lmMst2)*s2t*shiftst1*(-4*pow2(Mst1)*pow2(Mt) + 4*pow2(
        Mst2)*pow2(Mt) + 2*(lmMst1 - lmMst2)*pow2(Mst1)*pow2(Mst2)*pow2(s2t) +
        pow2(s2t)*(pow4(Mst1) - pow4(Mst2)))*pow4(Mst2) + 3*Mst2*pow3(s2t)*(
        1715*(192439 - 30400*OepS2 - 837000*S2 - 180*lmMst2*(857 + 3420*S2) +
        180*lmMst1*(677 + 180*lmMst2 + 3420*S2) - 16200*(pow2(lmMst1) + pow2(
        lmMst2)))*pow3(Mst2)*pow4(Mst1) + 128625*(8555 - 128*OepS2 - 57024*S2 -
        72*lmMst2*(29 + 36*S2) + 72*lmMst1*(29 - 12*lmMst2 + 36*S2) + 864*pow2(
        lmMst2))*pow2(Mst1)*pow5(Mst2) - 2*Mst2*(55313478 + 26068000*OepS2 +
        563377500*S2 + 315*lmMst2*(-423553 + 1675800*S2) - 315*lmMst1*(-423553
        + 166320*lmMst2 + 1675800*S2) + 26195400*(pow2(lmMst1) + pow2(lmMst2)))
        *pow6(Mst1) + 3087000*Dmglst2*(9*(55 - 34*lmMst1 + 34*lmMst2)*pow2(
        Mst2)*pow4(Mst1) - 36*(5 + 2*lmMst1 - 2*lmMst2)*pow2(Mst1)*pow4(Mst2) +
        (419 - 57*lmMst1 + 57*lmMst2)*pow6(Mst1) - 108*pow6(Mst2)) + 27783000*(
        1 + 2*lmMst2)*pow7(Mst2))))))/(5.00094e7*pow2(Mst1)*pow4(Msq)*pow4(
        Mst2)) + (pow2(Mt)*pow2(MuSUSY)*((3*pow2(Mt)*(6*Mst2*pow2(Mst1)*(3*
        Mst2*(10667 - 96*B4 + 96*D3 - 48*DN - 3072*lmMst1 - 384*(-25 + 9*
        lmMst1)*lmMst2 - 384*(-1 + lmMst1 - lmMst2)*lmMt - 224*OepS2 + 324*(-43
        + 14*lmMst1 - 14*lmMst2)*S2 - 384*(-13 + 4*lmMst1)*pow2(lmMst2) + 1536*
        pow3(lmMst2)) + 2*Dmglst2*(28405 - 288*B4 + 288*D3 - 144*DN + 5472*
        lmMst1 - 288*(-89 + 2*lmMst1)*lmMst2 + 576*(-2 + lmMst1 - lmMst2)*lmMt
        + 224*OepS2 - 324*(65 + 14*lmMst1 - 14*lmMst2)*S2 - 576*(-9 + 8*lmMst1)
        *pow2(lmMst2) + 4608*pow3(lmMst2))) + 144*(6*Dmglst2*(180 - 2*B4 + 2*D3
        - DN + 16*lmMst1 + 144*lmMst2 - 216*S2 - 16*(-2 + lmMst1)*pow2(lmMst2)
        + 16*pow3(lmMst2)) + Mst2*(436 - 6*B4 + 6*D3 - 3*DN - 48*lmMst1 + (408
        - 96*lmMst1)*lmMst2 + 24*lmMt - 972*S2 - 48*(-4 + lmMst1)*pow2(lmMst2)
        + 48*pow3(lmMst2)))*pow3(Mst2) + (383185 - 2592*B4 + 2592*D3 - 1296*DN
        - 187704*lmMst1 - 17440*OepS2 + 648*(-57 + 545*lmMst1 - 545*lmMst2)*S2
        - 7992*pow2(lmMst1) - 216*lmMst2*(-1733 + 630*lmMst1 + 26*pow2(lmMst1))
        - 216*(-859 + 246*lmMst1)*pow2(lmMst2) + 3456*lmMt*(3 + 6*lmMst2 - 2*
        lmMst1*(3 + lmMst2) + pow2(lmMst1) + pow2(lmMst2)) + 720*pow3(lmMst1) +
        58032*pow3(lmMst2))*pow4(Mst1)))/pow6(Mst2) + (4*T1ep*(-18*Mst2*pow2(
        Mst1)*(4*Dmglst2*(253*Mst2*Mt*s2t + 42*pow2(Mt) - 129*pow2(Mst2)*pow2(
        s2t)) + Mst2*(-272*Mst2*Mt*s2t - 252*pow2(Mt) + 1057*pow2(Mst2)*pow2(
        s2t))) + (Mt*(-65684*Dmglst2*s2t + 11292*Mst2*s2t) + 19620*pow2(Mt) + (
        34616*Dmglst2 - 39711*Mst2)*Mst2*pow2(s2t))*pow4(Mst1) - 54*s2t*(7*
        Dmglst2*(5*Mt - 2*Mst2*s2t) + 3*Mst2*(-7*Mt + 23*Mst2*s2t))*pow4(Mst2)
        - (60*pow2(Mst2)*pow2(s2t)*(pow4(Dmsqst2)*(-9*pow2(Mst1)*pow2(Mst2) +
        221*pow4(Mst1)) + Dmsqst2*(117*pow2(Mst1)*pow2(Mst2) + 221*pow4(Mst1) +
        27*pow4(Mst2))*pow6(Msq)))/pow8(Msq)))/pow6(Mst2) - 1458*Mt*s2t*((-128*
        Mst2*pow2(1 + lmMst2))/(9.*pow2(Mst1)) - ((1702429 + 257904*B4 - 648*DN
        - 748656*lmMst1 + 41904*pow2(lmMst1) + 216*lmMst2*(5971 - 6106*lmMst1 +
        576*pow2(lmMst1)) - 41904*(-34 + 15*lmMst1)*pow2(lmMst2) - 3456*pow3(
        lmMst1) + 507600*pow3(lmMst2))*pow4(Mst1))/(1458.*pow5(Mst2)) + ((8*
        OepS2*(2448*pow2(Mst1)*pow3(Mst2) + 5646*Mst2*pow4(Mst1) - Dmglst2*(
        9108*pow2(Mst1)*pow2(Mst2) + 32842*pow4(Mst1) + 945*pow4(Mst2)) + 567*
        pow5(Mst2)))/2187. + (S2*(Dmglst2*(72*(2489 + 3795*lmMst1 - 3795*
        lmMst2)*pow2(Mst1)*pow2(Mst2) + 10*(123113 + 98526*lmMst1 - 98526*
        lmMst2)*pow4(Mst1) + 81*(-453 + 350*lmMst1 - 350*lmMst2)*pow4(Mst2)) -
        15*(36*(169 + 136*lmMst1 - 136*lmMst2)*pow2(Mst1)*pow3(Mst2) + 2*(9185
        + 5646*lmMst1 - 5646*lmMst2)*Mst2*pow4(Mst1) + 81*(-1 + 14*lmMst1 - 14*
        lmMst2)*pow5(Mst2))))/405.)/pow6(Mst2) + Dmglst2*((-128*(-1 + 2*lmMst2
        + 3*pow2(lmMst2)))/(9.*pow2(Mst1)) + ((634.115454961134 - (3416*B4)/9.
         + (52*DN)/
        9. + (111356*lmMst1)/243. + (8*pow2(lmMst1))/9. - (8*lmMst2*(36802 -
        11421*lmMst1 + 1728*pow2(lmMst1)))/243. + (344*(-14 + 15*lmMst1)*pow2(
        lmMst2))/9. - (64*pow3(lmMst1))/9. - (1528*pow3(lmMst2))/3.)*pow4(Mst1)
        )/pow6(Mst2) - (14.76358024691358 + (1048*B4)/9. - (20*DN)/9. - (208*
        lmMst1)/9. + (64*pow2(lmMst1))/9. - (16*lmMst2*(-237 + 55*lmMst1 + 4*
        pow2(lmMst1)))/9. + (8*(280 - 121*lmMst1)*pow2(lmMst2))/9. + (344*pow3(
        lmMst2))/3. - (160*Dmsqst2*(2 + 5*lmMst1 - 5*lmMst2)*(pow3(Dmsqst2) +
        pow6(Msq)))/(3.*pow8(Msq)))/pow2(Mst2)) - (pow2(Mst1)*((2*Mst2*(75569 +
        13716*B4 - 54*DN - 33426*lmMst1 - 2376*pow2(lmMst1) + 54*lmMst2*(1427 -
        1012*lmMst1 + 16*pow2(lmMst1)) + 108*(642 - 203*lmMst1)*pow2(lmMst2) +
        21060*pow3(lmMst2) + (6480*Dmsqst2*(1 + 3*lmMst1 - 3*lmMst2)*(pow3(
        Dmsqst2) + pow6(Msq)))/pow8(Msq)))/243. + Dmglst2*(54.94732510288066 +
        248*B4 - 4*DN - (4564*lmMst1)/27. + (176*pow2(lmMst1))/9. + (4*lmMst2*(
        4993 - 1956*lmMst1 + 48*pow2(lmMst1)))/27. + (8*(482 - 331*lmMst1)*
        pow2(lmMst2))/9. + (2584*pow3(lmMst2))/9. + (160*Dmsqst2*(8 - 15*lmMst1
        + 15*lmMst2)*(pow3(Dmsqst2) + pow6(Msq)))/(3.*pow8(Msq)))))/pow4(Mst2)
        - (960*(11 + 8*lmMst1 - 8*lmMst2)*pow4(Dmsqst2) + 1920*Dmsqst2*(1 +
        lmMst1 - lmMst2)*pow6(Msq) + (9561 + 1760*B4 - 16*DN - 1984*lmMst1 -
        256*pow2(lmMst1) - 64*lmMst2*(-214 + 73*lmMst1 + 4*pow2(lmMst1)) - 32*(
        -268 + 57*lmMst1)*pow2(lmMst2) + 2080*pow3(lmMst2))*pow8(Msq))/(36.*
        Mst2*pow8(Msq))) - (486*pow2(s2t)*((1 - 2*lmMst2)*shiftst3*(2*(1 - 2*
        lmMst1 + 2*lmMst2)*pow2(Mst2)*pow4(Mst1) + 2*(1 - lmMst1 + lmMst2)*
        pow2(Mst1)*pow4(Mst2) + (2 - 6*lmMst1 + 6*lmMst2)*pow6(Mst1) + pow6(
        Mst2)) + (10*shiftst1*(pow4(Dmsqst2)*(3*(pow2(Mst1) + pow2(Mst2))*pow4(
        Mst2) - 4*pow2(lmMst2)*pow2(Mst1)*(pow2(Mst1)*pow2(Mst2) + pow4(Mst1) +
        pow4(Mst2)) + 2*lmMst1*(2*(-3 + lmMst2)*pow2(Mst2)*pow4(Mst1) + (-3 +
        2*lmMst2)*pow2(Mst1)*pow4(Mst2) + (-3 + 2*lmMst2)*pow6(Mst1)) + 2*
        lmMst2*(6*pow2(Mst2)*pow4(Mst1) + 2*pow2(Mst1)*pow4(Mst2) + 3*pow6(
        Mst1) - pow6(Mst2))) + ((pow2(Mst1) + pow2(Mst2))*pow4(Mst2) - 2*(
        lmMst1 - lmMst2)*pow2(Mst1)*(pow2(Mst1)*pow2(Mst2) + pow4(Mst1) + pow4(
        Mst2)))*(Dmsqst2*(3 - 2*lmMst2)*pow6(Msq) + (1 - 2*lmMst2)*pow8(Msq))))
        /pow8(Msq)))/(pow2(Mst1)*pow4(Mst2))))/1458. - pow2(s2t)*((-5*Dmsqst2*
        Mt*pow2(Sbeta)*pow4(Mst1)*((-18*(11 + 39*lmMst1 - 39*lmMst2)*s2t*pow2(
        Mst2)*(-1 + pow2(Sbeta)) + (Dmglst2*(4324 - 660*lmMst1 + 660*lmMst2) +
        9*(-167 + 5*lmMst1 - 5*lmMst2)*Mst2)*Mt*pow2(Sbeta))*pow3(Dmsqst2) + 4*
        Dmglst2*(1081 - 165*lmMst1 + 165*lmMst2)*Mt*pow2(Sbeta)*pow6(Msq)))/(
        36.*pow3(Mst2)*pow8(Msq)) + pow2(Mt)*((5*pow2(Sbeta)*pow4(Mst1)*(9*(-
        167 + 5*lmMst1 - 5*lmMst2)*Mst2*pow4(Dmsqst2) + 4*Dmglst2*Dmsqst2*(1081
        - 165*lmMst1 + 165*lmMst2)*(pow3(Dmsqst2) + pow6(Msq))))/(36.*pow3(
        Mst2)*pow8(Msq)) + pow2(MuSUSY)*(53.385802469135804 + (40*B4)/9. - (4*
        D3)/9. + (2*DN)/9. + (1672*lmMst1)/27. + (53*pow2(lmMst1))/9. - lmMst2*
        (129.92592592592592 - 72*lmMst1 + (13*pow2(lmMst1))/3.) + (2*(-470 +
        147*lmMst1)*pow2(lmMst2))/9. + (5*Dmsqst2*(1141 + 48*lmMst1*(7 - 3*
        lmMst2) - 264*lmMst2 + 144*pow2(lmMst2)))/(54.*pow2(Msq)) + (16*pow3(
        lmMst1))/9. - (271*pow3(lmMst2))/9. + (16*(1 + lmMst2)*(4*Dmglst2*
        lmMst2 + Mst2 + lmMst2*Mst2)*pow3(Mst2))/(9.*pow4(Mst1)) - (pow2(Mst1)*
        (-(Mst2*(434.270658436214 + (76*B4)/9. - (2*DN)/9. + (69088*lmMst1)/
        405. - (1313*pow2(lmMst1))/27. - (4*lmMst2*(16192 - 26430*lmMst1 +
        3465*pow2(lmMst1)))/405. + ((-5735 + 3072*lmMst1)*pow2(lmMst2))/27. + (
        Dmsqst2*(201.74098765432097 - (622*lmMst2)/15. + (2*lmMst1*(311 + 10*
        lmMst2))/15. - (22*pow2(lmMst1))/3. + 6*pow2(lmMst2)))/pow2(Msq) - (62*
        pow3(lmMst1))/27. - (2086*pow3(lmMst2))/27. + ((121.37234567901234 +
        lmMst1*(88.95555555555555 - (68*lmMst2)/3.) - (4003*lmMst2)/45. + (14*
        pow2(lmMst1))/3. + 18*pow2(lmMst2))*pow4(Dmsqst2))/pow8(Msq))) + (2*
        Dmglst2*(2695042 - 40500*B4 + 54000*D3 - 33750*DN - 326895*lmMst1 -
        324900*pow2(lmMst1) + 15*lmMst2*(-19607 - 129030*lmMst1 + 62550*pow2(
        lmMst1)) + 450*(5023 - 5610*lmMst1)*pow2(lmMst2) + 11250*pow3(lmMst1) +
        1575000*pow3(lmMst2) + (50625*Dmsqst2*(-23 + 42*lmMst1 - 42*lmMst2)*(
        pow3(Dmsqst2) + pow6(Msq)))/pow8(Msq)))/30375.))/pow3(Mst2) - (pow4(
        Mst1)*(Dmglst2*(585.1892843532082 - (8*B4)/3. + (32*D3)/9. - (20*DN)/9.
         - (20109937*lmMst1)/297675. - (15886*pow2(lmMst1))/945. + lmMst2*(
        17.112243218274966 - (144628*lmMst1)/945. + (1180*pow2(lmMst1))/9.) + (
        2*(80257 - 138180*lmMst1)*pow2(lmMst2))/945. - (92*pow3(lmMst1))/27. +
        (4448*pow3(lmMst2))/27.) - Mst2*(628.1736268201578 + (76*B4)/9. - (2*
        DN)/9. + (6317839*lmMst1)/396900. - (66307*pow2(lmMst1))/315. - lmMst2*
        (12.52907281431091 - (182909*lmMst1)/315. + (274*pow2(lmMst1))/3.) + (
        2*(-58301 + 37135*lmMst1)*pow2(lmMst2))/315. - (44*pow3(lmMst1))/9. - (
        1256*pow3(lmMst2))/9. + (Dmsqst2*(237.28785508324435 + (16526*lmMst2)/
        3969. + (2*lmMst1*(-8263 + 71820*lmMst2))/3969. - (520*pow2(lmMst1))/
        21. - (80*pow2(lmMst2))/7.)*(pow3(Dmsqst2) + pow6(Msq)))/pow8(Msq))))/
        pow5(Mst2) + (Dmglst2*(-88.0679012345679 + (8*B4)/3. - (32*D3)/9. + (
        20*DN)/9. - (88*lmMst1)/3. + lmMst2*(28 + (928*lmMst1)/9. - (64*pow2(
        lmMst1))/9.) + (26*pow2(lmMst1))/3. + (4*(-281 + 123*lmMst1)*pow2(
        lmMst2))/9. - (428*pow3(lmMst2))/9. + (pow2(Mst2)*(4.444444444444445 -
        (4*(99 + 16*lmMst1)*lmMst2)/9. - (82*pow2(lmMst2))/3. - (40*Dmsqst2*(
        pow3(Dmsqst2) + pow6(Msq)))/pow8(Msq)))/pow2(Mst1) - (80*Dmsqst2*(4 +
        lmMst1 - lmMst2)*(pow3(Dmsqst2) + pow6(Msq)))/(3.*pow8(Msq))))/Mst2 + (
        5*(4207 + lmMst1*(744 - 288*lmMst2) - 600*lmMst2 + 288*pow2(lmMst2))*
        pow4(Dmsqst2))/(108.*pow8(Msq)) - (pow2(Mst2)*(-30*(5 + 2*lmMst2)*pow4(
        Dmsqst2) + 30*Dmsqst2*(1 - 2*lmMst2)*pow6(Msq) + (103 + 186*lmMst2 +
        32*lmMst1*(1 + lmMst2) + 91*pow2(lmMst2))*pow8(Msq)))/(9.*pow2(Mst1)*
        pow8(Msq)) + (4*OepS2*(4*Dmglst2*(2322*pow2(Mst1)*pow2(Mst2) + 8654*
        pow4(Mst1) + 189*pow4(Mst2))*pow8(Msq) - 3*Mst2*(20*pow4(Dmsqst2)*(-9*
        pow2(Mst1)*pow2(Mst2) + 221*pow4(Mst1)) + 20*Dmsqst2*(117*pow2(Mst1)*
        pow2(Mst2) + 221*pow4(Mst1) + 27*pow4(Mst2))*pow6(Msq) + (6342*pow2(
        Mst1)*pow2(Mst2) + 13237*pow4(Mst1) + 1242*pow4(Mst2))*pow8(Msq))))/(
        2187.*pow5(Mst2)*pow8(Msq)) + S2*(921 + 138*lmMst1 - 138*lmMst2 + (
        Dmglst2*(3354 - 28*lmMst1 + 28*lmMst2))/Mst2 + (30*Dmsqst2*(-15 + 2*
        lmMst1 - 2*lmMst2))/pow2(Msq) + (pow4(Mst1)*(8*Dmglst2*(93919 - 12981*
        lmMst1 + 12981*lmMst2) + 3*Mst2*(194357 + 39711*lmMst1 - 39711*lmMst2 +
        (130*Dmsqst2*(95 + 102*lmMst1 - 102*lmMst2)*(pow3(Dmsqst2) + pow6(Msq))
        )/pow8(Msq))))/(81.*pow5(Mst2)) - (1740*pow4(Dmsqst2))/pow8(Msq) + (
        pow2(Mst1)*(4*Dmglst2*(15643 - 774*lmMst1 + 774*lmMst2)*pow8(Msq) - 3*
        Mst2*(10*(913 + 6*lmMst1 - 6*lmMst2)*pow4(Dmsqst2) + 10*Dmsqst2*(17 -
        78*lmMst1 + 78*lmMst2)*pow6(Msq) + (-11243 - 2114*lmMst1 + 2114*lmMst2)
        *pow8(Msq))))/(9.*pow3(Mst2)*pow8(Msq)))))))/Tbeta + MuSUSY*(-(Mt*pow2(
        z2)*(60*s2t*pow2(Mst1)*pow2(Mst2)*(9*pow2(Mst2)*(-2*Mt*MuSUSY*s2t - 6*
        Tbeta*pow2(Mt) + Tbeta*pow2(Mst2)*pow2(s2t)) + pow2(Mst1)*(442*Mt*
        MuSUSY*s2t + 106*Tbeta*pow2(Mt) + 37*Tbeta*pow2(Mst2)*pow2(s2t)))*pow4(
        Dmsqst2) - 60*s2t*xDmsqst2*pow2(Dmsqst2)*pow2(Mst2)*pow4(Msq)*(3*pow2(
        Mst1)*pow2(Mst2)*(-50*Mt*MuSUSY*s2t - 14*Tbeta*pow2(Mt) + 19*Tbeta*
        pow2(Mst2)*pow2(s2t)) + (-442*Mt*MuSUSY*s2t - 106*Tbeta*pow2(Mt) + 57*
        Tbeta*pow2(Mst2)*pow2(s2t))*pow4(Mst1) + 18*(-2*Mt*MuSUSY*s2t - 2*
        Tbeta*pow2(Mt) + Tbeta*pow2(Mst2)*pow2(s2t))*pow4(Mst2)) - 60*Dmsqst2*
        s2t*pow2(Mst2)*(18*pow2(Mst1)*pow2(Mst2)*(-13*Mt*MuSUSY*s2t - 5*Tbeta*
        pow2(Mt) + 5*Tbeta*pow2(Mst2)*pow2(s2t)) - 2*(221*Mt*MuSUSY*s2t + 53*
        Tbeta*pow2(Mt) - 52*Tbeta*pow2(Mst2)*pow2(s2t))*pow4(Mst1) + 27*(-2*Mt*
        MuSUSY*s2t - 2*Tbeta*pow2(Mt) + Tbeta*pow2(Mst2)*pow2(s2t))*pow4(Mst2))
        *pow6(Msq) - (-9*pow3(Mst2)*(-2*xDmglst2*pow2(Dmglst2)*(-4*s2t*(2536*
        MuSUSY - 63*Mst2*Tbeta)*pow2(Mt) + 2*Mst2*Mt*(209*MuSUSY + 3804*Mst2*
        Tbeta)*pow2(s2t) + 224*Tbeta*pow3(Mt) - 209*Tbeta*pow3(Mst2)*pow3(s2t))
        + 9*pow2(Mst2)*(2*s2t*(370*MuSUSY + 53*Mst2*Tbeta)*pow2(Mt) - Mst2*Mt*(
        4*MuSUSY + 555*Mst2*Tbeta)*pow2(s2t) - 28*Tbeta*pow3(Mt) + 2*Tbeta*
        pow3(Mst2)*pow3(s2t)) + 3*Dmglst2*Mst2*(7052*MuSUSY*s2t*pow2(Mt) -
        Mst2*Mt*(632*MuSUSY + 5289*Mst2*Tbeta)*pow2(s2t) + 140*Tbeta*pow3(Mt) +
        316*Tbeta*pow3(Mst2)*pow3(s2t))) + (-4*s2t*(172810*Dmglst2*MuSUSY +
        41010*Mst2*MuSUSY - 33084*Dmglst2*Mst2*Tbeta + 18093*Tbeta*pow2(Mst2))*
        pow2(Mt) + 2*Mst2*Mt*(42392*Dmglst2*MuSUSY - 35823*Mst2*MuSUSY +
        105585*Dmglst2*Mst2*Tbeta + 18531*Tbeta*pow2(Mst2))*pow2(s2t) + 8*(
        4905*MuSUSY - 32426*Dmglst2*Tbeta + 6090*Mst2*Tbeta)*pow3(Mt) + (-
        25328*Dmglst2 + 20685*Mst2)*Tbeta*pow3(Mst2)*pow3(s2t))*pow4(Mst1) + 3*
        pow2(Mst1)*(-3*Dmglst2*Mst2*(-80*Mst2*s2t*(-569*MuSUSY + 30*Mst2*Tbeta)
        *pow2(Mt) - 3*Mt*(1264*MuSUSY + 6091*Mst2*Tbeta)*pow2(Mst2)*pow2(s2t) +
        96*(7*MuSUSY + 62*Mst2*Tbeta)*pow3(Mt) + 948*Tbeta*pow3(s2t)*pow4(Mst2)
        ) + 2*xDmglst2*pow2(Dmglst2)*(4*Mst2*s2t*(-12479*MuSUSY + 1080*Mst2*
        Tbeta)*pow2(Mt) + Mt*(-1060*MuSUSY + 14613*Mst2*Tbeta)*pow2(Mst2)*pow2(
        s2t) + 56*(3*MuSUSY + 314*Mst2*Tbeta)*pow3(Mt) + 1157*Tbeta*pow3(s2t)*
        pow4(Mst2)) + 3*pow2(Mst2)*(-2*Mst2*s2t*(6368*MuSUSY + 1735*Mst2*Tbeta)
        *pow2(Mt) + Mt*(-3364*MuSUSY + 4557*Mst2*Tbeta)*pow2(Mst2)*pow2(s2t) +
        16*(63*MuSUSY + 110*Mst2*Tbeta)*pow3(Mt) + 1700*Tbeta*pow3(s2t)*pow4(
        Mst2))))*pow8(Msq)))/(972.*Tbeta*pow6(Mst2)*pow8(Msq)) + xDR2DRMOD*((
        64*pow4(Mt)*(9*(1 + lmMst2)*Mst2*((3 - 21*lmMst2 + 2*lmMst1*(8 + 3*
        lmMst2 - 3*lmMt) + 5*lmMt + 6*lmMst2*lmMt - 6*pow2(lmMst2))*pow2(Mst2)*
        pow4(Mst1) + (-1 - 7*lmMst2 + 2*lmMst1*(2 + lmMst2 - lmMt) + 3*lmMt +
        2*lmMst2*lmMt - 2*pow2(lmMst2))*pow2(Mst1)*pow4(Mst2) + (12 - 45*lmMst2
        + 2*lmMst1*(19 + 6*lmMst2 - 6*lmMt) + 7*lmMt + 12*lmMst2*lmMt - 12*
        pow2(lmMst2))*pow6(Mst1) - 2*(1 + lmMst2)*pow6(Mst2)) + Dmglst2*(pow2(
        Mst2)*(263 + lmMst2*(962 - 213*lmMt) - 177*lmMt + (393 - 18*lmMt)*pow2(
        lmMst2) - 18*lmMst1*(26 - lmMst2*(-17 + lmMt) - 7*lmMt + pow2(lmMst2))
        + 18*pow3(lmMst2))*pow4(Mst1) - pow2(Mst1)*(17*(-7 + 3*lmMt) + lmMst2*(
        -224 + 15*lmMt) - 3*(5 + 6*lmMt)*pow2(lmMst2) - 18*lmMst1*(-4 + lmMt -
        lmMst2*(1 + lmMt) + pow2(lmMst2)) + 18*pow3(lmMst2))*pow4(Mst2) + (290
        + lmMst2*(2447 - 645*lmMt) - 375*lmMt - 3*(-509 + 60*lmMt)*pow2(lmMst2)
        - 18*lmMst1*(87 + lmMst2*(71 - 10*lmMt) - 22*lmMt + 10*pow2(lmMst2)) +
        180*pow3(lmMst2))*pow6(Mst1) - 18*(-1 + 2*lmMst2 + 3*pow2(lmMst2))*
        pow6(Mst2))))/(81.*pow2(Mst1)*pow6(Mst2)) - (16*pow2(Mt)*pow2(s2t)*((1
        + lmMst2)*Mst2*(2*(2 - 5*lmMst2 + lmMst1*(5 + 2*lmMst2) - 2*pow2(
        lmMst2))*pow2(Mst2)*pow4(Mst1) + 2*(-2*lmMst2*(2 + lmMst2) + lmMst1*(3
        + 2*lmMst2))*pow2(Mst1)*pow4(Mst2) + (5 - 12*lmMst2 + 4*lmMst1*(3 +
        lmMst2) - 4*pow2(lmMst2))*pow6(Mst1) - 2*(1 + lmMst2)*pow6(Mst2)) +
        Dmglst2*(2*pow2(Mst2)*(8 + 19*lmMst2 - 5*pow2(lmMst2) + lmMst1*(-7 + 5*
        lmMst2 + 6*pow2(lmMst2)) - 6*pow3(lmMst2))*pow4(Mst1) + 2*pow2(Mst1)*(7
        + 9*lmMst2 - 8*pow2(lmMst2) + lmMst1*(-3 + 5*lmMst2 + 6*pow2(lmMst2)) -
        6*pow3(lmMst2))*pow4(Mst2) + (17 + 51*lmMst2 - 4*pow2(lmMst2) + 4*
        lmMst1*(-6 + lmMst2 + 3*pow2(lmMst2)) - 12*pow3(lmMst2))*pow6(Mst1) -
        2*(-1 + 2*lmMst2 + 3*pow2(lmMst2))*pow6(Mst2))))/(3.*pow2(Mst1)*pow4(
        Mst2)) - (15*shiftst1*(4*s2t*(pow2(Mst1) - pow2(Mst2))*pow3(Mt) - Mt*
        pow3(s2t)*(2*(lmMst1 - lmMst2)*pow2(Mst1)*pow2(Mst2) + pow4(Mst1) -
        pow4(Mst2)))*(1 + Dmsqst2/pow2(Msq) + pow4(Dmsqst2)/pow8(Msq)))/(2.*
        pow2(Mst1)) + (MuSUSY*s2t*pow2(Mt)*((64*Mt*((1 + lmMst2)*Mst2*(2*(1 -
        10*lmMst2 + 4*lmMst1*(2 + lmMst2) - 4*pow2(lmMst2))*pow2(Mst2)*pow4(
        Mst1) - 2*(1 + 5*lmMst2 - lmMst1*(3 + 2*lmMst2) + 2*pow2(lmMst2))*pow2(
        Mst1)*pow4(Mst2) + (7 - 32*lmMst2 + 4*lmMst1*(7 + 3*lmMst2) - 12*pow2(
        lmMst2))*pow6(Mst1) - 2*(1 + lmMst2)*pow6(Mst2)) + Dmglst2*(4*pow2(
        Mst2)*(8 + 13*lmMst2 - 8*pow2(lmMst2) + lmMst1*(-5 + 5*lmMst2 + 6*pow2(
        lmMst2)) - 6*pow3(lmMst2))*pow4(Mst1) + (49 + 103*lmMst2 - 36*(1 +
        lmMst2)*pow2(lmMst2) + 4*lmMst1*(-11 + 6*lmMst2 + 9*pow2(lmMst2)))*
        pow6(Mst1) + 2*(pow2(Mst1)*(8 + 7*lmMst2 - 11*pow2(lmMst2) + lmMst1*(-3
        + 5*lmMst2 + 6*pow2(lmMst2)) - 6*pow3(lmMst2))*pow4(Mst2) + (1 - 2*
        lmMst2 - 3*pow2(lmMst2))*pow6(Mst2)))))/(9.*pow2(Mst1)*pow6(Mst2)) - (
        s2t*(756 - 2024*lmMst1 + 128*pow2(lmMst1) + 8*lmMst2*(363 - 332*lmMst1
        + 16*pow2(lmMst1)) + 4*(707 - 246*lmMst1)*pow2(lmMst2) + 856*pow3(
        lmMst2) - (64*(1 + lmMst2)*(4*Dmglst2*lmMst2 + Mst2 + lmMst2*Mst2)*
        pow3(Mst2))/pow4(Mst1) - 540*shiftst1*(1 + pow2(Mst2)/pow2(Mst1) - (2*(
        lmMst1 - lmMst2)*(pow2(Mst1)*pow2(Mst2) + pow4(Mst1) + pow4(Mst2)))/
        pow4(Mst2))*(1 + Dmsqst2/pow2(Msq) + pow4(Dmsqst2)/pow8(Msq)) + (150*(-
        1 + 2*lmMst1 - 2*lmMst2)*(pow4(Dmsqst2) - 2*Dmsqst2*pow6(Msq)))/pow8(
        Msq) + (pow2(Mst2)*(-105*pow4(Dmsqst2) + 300*Dmsqst2*pow6(Msq) + 4*(205
        + 252*lmMst2 + 32*lmMst1*(1 + lmMst2) + 91*pow2(lmMst2))*pow8(Msq)))/(
        pow2(Mst1)*pow8(Msq)) + (4*(pow2(Mst2)*pow4(Mst1)*(60*Dmglst2*Dmsqst2*(
        lmMst1 - lmMst2)*(pow3(Dmsqst2) + pow6(Msq)) - 15*(lmMst1 - lmMst2)*
        Mst2*(13*pow4(Dmsqst2) + 10*Dmsqst2*pow6(Msq)) + 4*Dmglst2*(48 + 4*
        lmMst2*(31 + 36*pow2(lmMst1)) + 278*pow2(lmMst2) - lmMst1*(44 + 278*
        lmMst2 + 379*pow2(lmMst2)) + 235*pow3(lmMst2))*pow8(Msq) + 2*Mst2*(32 +
        285*lmMst2 + 144*(1 + lmMst2)*pow2(lmMst1) + 444*pow2(lmMst2) - lmMst1*
        (253 + 588*lmMst2 + 379*pow2(lmMst2)) + 235*pow3(lmMst2))*pow8(Msq)) +
        2*pow6(Mst1)*(-75*Dmsqst2*(lmMst1 - lmMst2)*Mst2*(pow3(Dmsqst2) + pow6(
        Msq)) + 2*Dmglst2*(48 + 4*lmMst2*(45 + 68*pow2(lmMst1)) + 310*pow2(
        lmMst2) - lmMst1*(92 + 310*lmMst2 + 635*pow2(lmMst2)) + 363*pow3(
        lmMst2))*pow8(Msq) + Mst2*(40 + 277*lmMst2 + 272*(1 + lmMst2)*pow2(
        lmMst1) + 556*pow2(lmMst2) - lmMst1*(237 + 828*lmMst2 + 635*pow2(
        lmMst2)) + 363*pow3(lmMst2))*pow8(Msq)) - Dmglst2*pow4(Mst2)*(-15*((19
        - 38*lmMst1 + 38*lmMst2)*pow2(Mst1) + pow2(Mst2))*pow4(Dmsqst2) + 30*
        Dmsqst2*((1 - 2*lmMst1 + 2*lmMst2)*pow2(Mst1) + pow2(Mst2))*pow6(Msq) -
        2*((-20 + 2*(99 + 16*lmMst1)*lmMst2 + 123*pow2(lmMst2))*pow2(Mst2) +
        pow2(Mst1)*(60 + 206*lmMst2 + 32*lmMst2*pow2(lmMst1) + lmMst1*(8 - 460*
        lmMst2 - 246*pow2(lmMst2)) + 519*pow2(lmMst2) + 214*pow3(lmMst2)))*
        pow8(Msq))))/(pow2(Mst1)*pow5(Mst2)*pow8(Msq))))/36.))/Tbeta + (Mt*
        pow3(s2t)*(pow2(Mst1)*pow3(Mst2)*(15*pow4(Dmsqst2)*(14*(lmMst1 -
        lmMst2)*pow2(Mst1)*pow2(Mst2) + 10*pow4(Mst1) - 7*pow4(Mst2)) - 300*
        Dmsqst2*(2*(lmMst1 - lmMst2)*pow2(Mst1)*pow2(Mst2) + pow4(Mst1) - pow4(
        Mst2))*pow6(Msq)) + 2*Dmglst2*(60*lmMst2*pow4(Mst1)*pow4(Mst2)*(pow4(
        Dmsqst2) - 2*Dmsqst2*pow6(Msq)) + pow2(Mst1)*pow2(Mst2)*(pow4(Dmsqst2)*
        (-570*pow4(Mst1) + 75*pow4(Mst2)) + 60*Dmsqst2*(pow4(Mst1) - pow4(Mst2)
        )*pow6(Msq) + 16*(20*pow2(Mst1)*pow2(Mst2) + 9*pow4(Mst1) - 5*pow4(
        Mst2))*pow8(Msq)) - 4*lmMst1*(15*pow4(Mst1)*pow4(Mst2)*(pow4(Dmsqst2) -
        2*Dmsqst2*pow6(Msq)) + pow8(Msq)*(2*(-4 + 246*lmMst2 + 123*pow2(lmMst2)
        )*pow4(Mst1)*pow4(Mst2) + 32*(3 + 3*lmMst2 + 16*pow2(lmMst2))*pow2(
        Mst2)*pow6(Mst1) - 32*lmMst2*pow2(Mst1)*pow6(Mst2) + 32*(3 + 2*lmMst2 +
        16*pow2(lmMst2))*pow8(Mst1))) + 4*lmMst2*pow8(Msq)*(32*pow2(lmMst1)*
        pow4(Mst1)*(8*pow2(Mst1)*pow2(Mst2) + 8*pow4(Mst1) + pow4(Mst2)) +
        lmMst2*(2*(198 + 107*lmMst2)*pow4(Mst1)*pow4(Mst2) + (37 + 256*lmMst2)*
        pow2(Mst2)*pow6(Mst1) + 155*pow2(Mst1)*pow6(Mst2) + 64*(1 + 4*lmMst2)*
        pow8(Mst1) - 32*pow8(Mst2))) + 8*lmMst2*pow8(Msq)*(4*pow4(Mst1)*pow4(
        Mst2) + 21*pow2(Mst2)*pow6(Mst1) + 115*pow2(Mst1)*pow6(Mst2) + 56*pow8(
        Mst1) - 16*pow8(Mst2))) + 4*Mst2*pow8(Msq)*(2*(-8 + 237*lmMst2 + 16*(1
        + lmMst2)*pow2(lmMst1) + 308*pow2(lmMst2) - lmMst1*(269 + 348*lmMst2 +
        123*pow2(lmMst2)) + 107*pow3(lmMst2))*pow4(Mst1)*pow4(Mst2) + pow2(
        Mst2)*(-125 - 156*lmMst2 - 512*lmMst1*lmMst2*(1 + lmMst2) + 256*(1 +
        lmMst2)*pow2(lmMst1) + 181*pow2(lmMst2) + 256*pow3(lmMst2))*pow6(Mst1)
        + (221 + 284*lmMst2 + 32*lmMst1*(1 + lmMst2) + 107*pow2(lmMst2))*pow2(
        Mst1)*pow6(Mst2) + 16*(1 + lmMst2)*(1 + lmMst1*(2 - 32*lmMst2) - 2*
        lmMst2 + 16*pow2(lmMst1) + 16*pow2(lmMst2))*pow8(Mst1) - 16*pow2(1 +
        lmMst2)*pow8(Mst2))))/(72.*pow3(Mst2)*pow4(Mst1)*pow8(Msq)) - (s2t*
        pow3(Mt)*(pow2(Mst1)*pow5(Mst2)*(45*(10*pow2(Mst1) - 7*pow2(Mst2))*
        pow4(Dmsqst2) - 900*Dmsqst2*(pow2(Mst1) - pow2(Mst2))*pow6(Msq)) + 12*
        Mst2*pow8(Msq)*(-((109 + 76*lmMst2 - 21*pow2(lmMst2) + 64*lmMst1*pow2(1
        + lmMst2) - 32*((1 + lmMst2)*pow2(lmMst1) + pow3(lmMst2)))*pow4(Mst1)*
        pow4(Mst2)) + (173 + 188*lmMst2 + 32*lmMst1*(1 + lmMst2) + 59*pow2(
        lmMst2))*pow2(Mst1)*pow6(Mst2) + 32*(1 + lmMst2)*(pow2(1 - lmMst1 +
        lmMst2)*pow2(Mst2)*pow6(Mst1) + (1 + 3*lmMst2 - lmMst1*(3 + 2*lmMst2) +
        pow2(lmMst1) + pow2(lmMst2))*pow8(Mst1)) - 16*pow2(1 + lmMst2)*pow8(
        Mst2)) + 4*Dmglst2*(pow2(Mst1)*pow4(Mst2)*(45*(-19*pow2(Mst1) + pow2(
        Mst2))*pow4(Dmsqst2) + 90*Dmsqst2*(pow2(Mst1) - pow2(Mst2))*pow6(Msq))
        - 2*pow8(Msq)*((52 + 370*lmMst2 - 96*lmMst2*pow2(lmMst1) + 81*pow2(
        lmMst2) + 96*lmMst1*(-1 + lmMst2 + 2*pow2(lmMst2)) - 96*pow3(lmMst2))*
        pow4(Mst1)*pow4(Mst2) - 96*(lmMst2*pow2(lmMst1) + lmMst1*(4 + 2*lmMst2
        - 2*pow2(lmMst2)) + (-4 + lmMst2)*pow2(1 + lmMst2))*pow2(Mst2)*pow6(
        Mst1) - 3*(-52 + 2*(51 + 16*lmMst1)*lmMst2 + 59*pow2(lmMst2))*pow2(
        Mst1)*pow6(Mst2) - 96*(-6 - 14*lmMst2 + lmMst2*pow2(lmMst1) + lmMst1*(9
        + 6*lmMst2 - 2*pow2(lmMst2)) - 6*pow2(lmMst2) + pow3(lmMst2))*pow8(
        Mst1) + 96*lmMst2*(1 + lmMst2)*pow8(Mst2)))))/(54.*pow4(Mst1)*pow5(
        Mst2)*pow8(Msq))))) + (Mt*xMst*(27*(lmMst1 - lmMst2)*Mt*oneLoopFlag*
        pow2(MuSUSY)*pow2(s2t)*pow3(Mst2) - Al4p*twoLoopFlag*(8*Dmglst2*Mst2*
        MuSUSY*(s2t*(36*(-1 + 3*lmMst1 - 3*lmMst2)*Mst2*Tbeta + MuSUSY*(785 +
        6*lmMst1*(85 - 42*lmMst2) - 438*lmMst2 + 252*pow2(lmMst2)))*pow2(Mt) -
        Mst2*Mt*(3*Mst2*Tbeta*(50 + lmMst1*(51 - 18*lmMst2) - 51*lmMst2 + 18*
        pow2(lmMst2)) + MuSUSY*(49 - 84*lmMst2 + lmMst1*(84 - 36*lmMst2*(-1 +
        xDR2DRMOD)) + 36*(-1 + xDR2DRMOD)*pow2(lmMst2)))*pow2(s2t) + Tbeta*(4*(
        7 - 381*lmMst2 + lmMst1*(327 + 72*lmMst2 - 81*lmMt) + 54*lmMt + 81*
        lmMst2*lmMt - 72*pow2(lmMst2))*pow3(Mt) + 2*(1 + 3*lmMst1 - 3*lmMst2)*
        pow3(Mst2)*pow3(s2t))) - MuSUSY*pow2(Mst2)*(-8*s2t*(18*Mst2*Tbeta*(-
        lmMst1 + lmMst2 - 2*lmMst1*lmMst2 + pow2(lmMst1) + pow2(lmMst2)) +
        MuSUSY*(193 + 474*lmMst2 - 6*lmMst1*(67 + 42*lmMst2) + 252*pow2(lmMst2)
        ))*pow2(Mt) + 8*Mst2*Mt*(3*Mst2*Tbeta*(1 + 33*lmMst2 - 3*lmMst1*(11 +
        6*lmMst2) + 18*pow2(lmMst2)) + MuSUSY*(1 + 3*lmMst2*(-37 + 6*xDR2DRMOD)
        - 3*lmMst1*(-37 + 6*lmMst2*(-12 + xDR2DRMOD) + 6*xDR2DRMOD) - 81*pow2(
        lmMst1) + 9*(-15 + 2*xDR2DRMOD)*pow2(lmMst2)))*pow2(s2t) + Tbeta*(32*(
        83 - 96*lmMst2 + 3*lmMst1*(29 + 12*lmMst2 - 9*lmMt) + 9*lmMt + 27*
        lmMst2*lmMt - 36*pow2(lmMst2))*pow3(Mt) + (5 + lmMst1*(6 - 288*lmMst2)
        - 6*lmMst2 + 144*pow2(lmMst1) + 144*pow2(lmMst2))*pow3(Mst2)*pow3(s2t))
        ) + 2*xDmglst2*pow2(Dmglst2)*(-(Mst2*Mt*pow2(s2t)*(18*(85 - 36*lmMst1 +
        36*lmMst2)*Mst2*MuSUSY*Tbeta + (314 + 24*lmMst2*(29 - 6*xDR2DRMOD) -
        24*lmMst1*(29 + 3*lmMst2*(-1 + xDR2DRMOD) - 6*xDR2DRMOD) + 72*(-1 +
        xDR2DRMOD)*pow2(lmMst2))*pow2(MuSUSY) + 15*(-43 + 60*lmMst1 - 60*
        lmMst2)*pow2(Mst2)*(-1 + pow2(Sbeta))*pow2(Sbeta))) + MuSUSY*(6*s2t*
        pow2(Mt)*(8*(143 - 18*lmMst1 + 18*lmMst2)*MuSUSY + Mst2*Tbeta*(371 +
        lmMst2*(552 - 600*pow2(Sbeta)) - 430*pow2(Sbeta) + 24*lmMst1*(-23 + 25*
        pow2(Sbeta)))) + Tbeta*(8*(398 + lmMst2*(2085 - 396*lmMt) - 18*lmMst1*(
        95 + 22*lmMst2 - 22*lmMt) - 375*lmMt + 396*pow2(lmMst2))*pow3(Mt) + 4*(
        1 - 24*lmMst1 + 24*lmMst2)*pow3(Mst2)*pow3(s2t))))))*pow6(Mst1))/(108.*
        Tbeta*pow9(Mst2)) - (Al4p*Mt*MuSUSY*twoLoopFlag*z2*(8*(-(MuSUSY*s2t) +
        6*Mt*Tbeta)*pow2(Mt)*pow4(Mst1)*pow4(Mst2) + 8*(-(MuSUSY*s2t) + 8*Mt*
        Tbeta)*xMst*pow2(Mst2)*pow2(Mt)*pow6(Mst1) + 4*Mt*xDmglst2*xMst*pow2(
        Dmglst2)*(-48*Mt*MuSUSY*s2t + 88*Tbeta*pow2(Mt) + Mst2*(MuSUSY + 9*
        Mst2*Tbeta)*pow2(s2t))*pow6(Mst1) + pow2(Mst1)*(-8*MuSUSY*s2t*pow2(Mt)
        + 32*Tbeta*pow3(Mt) + Tbeta*pow3(Mst2)*pow3(s2t))*pow6(Mst2) + 2*
        Dmglst2*(4*Mt*(-13*Mt*MuSUSY*s2t - 6*Tbeta*pow2(Mt) + Mst2*(MuSUSY + 3*
        Mst2*Tbeta)*pow2(s2t))*pow3(Mst2)*pow4(Mst1) + 4*Mt*s2t*(-9*Mt*MuSUSY +
        Mst2*s2t*(MuSUSY + 3*Mst2*Tbeta))*pow2(Mst1)*pow5(Mst2) + 4*Mst2*Mt*
        xMst*(-17*Mt*MuSUSY*s2t - 16*Tbeta*pow2(Mt) + Mst2*(MuSUSY + 3*Mst2*
        Tbeta)*pow2(s2t))*pow6(Mst1) + (-20*MuSUSY*s2t*pow2(Mt) + Mst2*Mt*(4*
        MuSUSY + 15*Mst2*Tbeta)*pow2(s2t) + 8*Tbeta*pow3(Mt) - 2*Tbeta*pow3(
        Mst2)*pow3(s2t))*pow7(Mst2)) + (-8*MuSUSY*s2t*pow2(Mt) + 2*Mst2*Mt*(
        MuSUSY + 3*Mst2*Tbeta)*pow2(s2t) + 16*Tbeta*pow3(Mt) - Tbeta*pow3(Mst2)
        *pow3(s2t))*pow8(Mst2)))/(3.*Tbeta*pow9(Mst2));
}

/**
 * @return returns the susy log^0 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double H6bq22g::calc_coef_at_as2_no_sm_logs_log0() const {

   const double result =
      (588*Mst2*pow2(Dmglst2)*pow4(Msq)*(1134000*pow2(Mst1)*pow4(Mst2)*(-8*(-
        179 + 26*z2)*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 768*Mst2*s2t*pow3(Mt) -
        192*Mt*pow3(Mst2)*pow3(s2t) + 16*(-155 + 26*z2)*pow4(Mt) + (-235 + 26*
        z2)*pow4(Mst2)*pow4(s2t)) + 15*pow2(Mst2)*pow4(Mst1)*(-360*pow2(Mst2)*
        pow2(Mt)*pow2(s2t)*(-2164753 + 3360*B4 - 3360*D3 + 1680*DN - 7840*OepS2
        + 169452*S2 + 11760*T1ep + 467818*z2 + 1945895*z3 - 9240*z4 + 8820*
        pow2(z2)) + 320*Mst2*s2t*(-1243547 + 15680*OepS2 + 953208*S2 - 23520*
        T1ep - 2985438*z2 + 3393915*z3 - 11760*z4 - 17640*pow2(z2))*pow3(Mt) +
        560*Mt*(-282298 + 61560*B4 - 1620*DN - 2240*OepS2 - 111213*S2 + 3360*
        T1ep - 281922*z2 + 983190*z3 + 6540*z4 - 114120*pow2(z2))*pow3(Mst2)*
        pow3(s2t) - 16*(-46632377 + 744800*OepS2 + 29679210*S2 - 1117200*T1ep -
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
        s2t))*pow6(Mst2)) + 4*Dmglst2*pow2(Mst1)*(-46305000*Mst2*pow2(Dmsqst2)*
        (4*pow2(Mst1)*pow3(Mst2)*(36*(23 - 24*z2)*pow2(Mst2)*pow2(Mt)*pow2(s2t)
        + 24*Mst2*s2t*(289 - 342*z2 + 216*z3)*pow3(Mt) - 216*Mt*(-2 + 2*z2 -
        z3)*pow3(Mst2)*pow3(s2t) + 8*(-1477 + 1242*z2 - 540*z3)*pow4(Mt) + 27*(
        -2 + 4*z2 - 3*z3)*pow4(Mst2)*pow4(s2t)) - 9*Mst2*pow4(Mst1)*(8*(-113 +
        48*z2)*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 192*Mst2*s2t*(51*z2 - 28*(2 +
        z3))*pow3(Mt) - 96*Mt*(5*z2 + 3*(-4 + z3))*pow3(Mst2)*pow3(s2t) - 48*(-
        401 + 308*z2 - 144*z3)*pow4(Mt) + 3*(-75 + 32*z2 + 8*z3)*pow4(Mst2)*
        pow4(s2t)) + 108*(-3 + z2)*pow2(-4*pow2(Mt) + pow2(Mst2)*pow2(s2t))*
        pow5(Mst2) - 6*s2t*(Mst2*s2t*(1931 - 1098*z2 - 180*z3)*pow2(Mt) - 72*
        Mt*(-20 + 11*z2)*pow2(Mst2)*pow2(s2t) + 72*(-465 + 348*z2 - 176*z3)*
        pow3(Mt) + (38 + 63*z2 - 90*z3)*pow3(Mst2)*pow3(s2t))*pow6(Mst1)) -
        46305000*Dmsqst2*Mst2*pow2(Msq)*(4*pow2(Mst1)*pow3(Mst2)*(36*(23 - 24*
        z2)*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 24*Mst2*s2t*(289 - 342*z2 + 216*z3)
        *pow3(Mt) - 216*Mt*(-2 + 2*z2 - z3)*pow3(Mst2)*pow3(s2t) + 8*(-1477 +
        1242*z2 - 540*z3)*pow4(Mt) + 27*(-2 + 4*z2 - 3*z3)*pow4(Mst2)*pow4(s2t)
        ) - 9*Mst2*pow4(Mst1)*(8*(-113 + 48*z2)*pow2(Mst2)*pow2(Mt)*pow2(s2t) +
        192*Mst2*s2t*(51*z2 - 28*(2 + z3))*pow3(Mt) - 96*Mt*(5*z2 + 3*(-4 + z3)
        )*pow3(Mst2)*pow3(s2t) - 48*(-401 + 308*z2 - 144*z3)*pow4(Mt) + 3*(-75
        + 32*z2 + 8*z3)*pow4(Mst2)*pow4(s2t)) + 108*(-3 + z2)*pow2(-4*pow2(Mt)
        + pow2(Mst2)*pow2(s2t))*pow5(Mst2) - 6*s2t*(Mst2*s2t*(1931 - 1098*z2 -
        180*z3)*pow2(Mt) - 72*Mt*(-20 + 11*z2)*pow2(Mst2)*pow2(s2t) + 72*(-465
        + 348*z2 - 176*z3)*pow3(Mt) + (38 + 63*z2 - 90*z3)*pow3(Mst2)*pow3(s2t)
        )*pow6(Mst1)) + pow4(Msq)*(-1764*pow2(Mst2)*pow4(Mst1)*(-112*pow2(Mst2)
        *pow2(Mt)*pow2(s2t)*(1216808 + 150000*OepS2 - 376650*S2 - 225000*T1ep -
        3408750*z2 - 833625*z3 - 112500*z4 - 168750*pow2(z2)) + 1050*Mst2*s2t*(
        -322421 + 73760*OepS2 - 564876*S2 - 110640*T1ep - 3060894*z2 + 1606920*
        z3 - 55320*z4 - 82980*pow2(z2))*pow3(Mt) - 175*Mt*(27211 + 36720*B4 +
        1080*DN + 64160*OepS2 - 1515564*S2 - 96240*T1ep - 1240446*z2 + 245400*
        z3 - 12480*z4 - 72180*pow2(z2))*pow3(Mst2)*pow3(s2t) - 8*(98884301 +
        4508000*OepS2 - 17853750*S2 - 6762000*T1ep - 477659550*z2 + 302883000*
        z3 - 3381000*z4 - 5071500*pow2(z2))*pow4(Mt) - 7*(-95041 + 81000*B4 -
        108000*D3 + 67500*DN - 432000*OepS2 - 7425000*S2 + 648000*T1ep +
        5325300*z2 + 4227750*z3 - 121500*z4)*pow4(Mst2)*pow4(s2t)) - 66150*
        pow2(Mst1)*pow4(Mst2)*(672*(-11674 + 120*B4 - 120*D3 + 60*DN + 17091*S2
        + 4849*z2 + 2195*z3 - 540*z4)*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 140*Mst2*
        s2t*(1120*OepS2 + 3*(-10061 + 14580*S2 - 560*T1ep - 25134*z2 + 72*z3 +
        2024*z4 - 420*pow2(z2)))*pow3(Mt) - 7*Mt*(69997 + 188640*B4 - 3600*DN +
        5600*OepS2 + 146772*S2 - 8400*T1ep - 1043682*z2 + 1816440*z3 + 32520*z4
        - 317340*pow2(z2))*pow3(Mst2)*pow3(s2t) - 16*(1547066 + 7840*OepS2 +
        258066*S2 - 11760*T1ep - 1256046*z2 - 946785*z3 - 5880*z4 - 8820*pow2(
        z2))*pow4(Mt) + 35*(-15707 + 432*B4 - 576*D3 + 360*DN + 224*OepS2 +
        543348*S2 - 336*T1ep - 11802*z2 - 105690*z3 - 2544*z4 - 2844*pow2(z2))*
        pow4(Mst2)*pow4(s2t)) + (72*pow2(Mst2)*pow2(Mt)*pow2(s2t)*(-3650009897
        + 2110136000*OepS2 - 52385772600*S2 - 3165204000*T1ep - 40405228500*z2
        + 755286000*z3 - 1582602000*z4 - 2373903000*pow2(z2)) - 274400*Mst2*
        s2t*(-20964397 + 2058400*OepS2 - 47545272*S2 - 3087600*T1ep - 59449002*
        z2 + 39214920*z3 - 1543800*z4 - 2315700*pow2(z2))*pow3(Mt) + 17150*Mt*(
        -31897243 + 2491360*OepS2 - 90290268*S2 - 3737040*T1ep - 51499698*z2 +
        33912840*z3 - 1868520*z4 - 2802780*pow2(z2))*pow3(Mst2)*pow3(s2t) + 80*
        (-54946675289 + 3063401600*OepS2 - 52952863320*S2 - 4595102400*T1ep -
        225236187120*z2 + 243414477600*z3 - 2297551200*z4 - 3446326800*pow2(z2)
        )*pow4(Mt) + (119405394763 - 11522056000*OepS2 + 478191735000*S2 +
        17283084000*T1ep + 233890773900*z2 - 111792103500*z3 + 8641542000*z4 +
        12962313000*pow2(z2))*pow4(Mst2)*pow4(s2t))*pow6(Mst1) - 333396000*(2*
        Mt + Mst2*s2t)*(14*Mt*(3 + z2) + Mst2*s2t*(5 + 7*z2))*pow2(-2*Mt +
        Mst2*s2t)*pow6(Mst2))) + 55566000*pow3(log(pow2(Mst1)/pow2(Mst2)))*
        pow4(Msq)*pow4(Mst1)*(4*pow2(Mst1)*pow3(Mst2)*(1166*pow2(Mst2)*pow2(Mt)
        *pow2(s2t) - 2304*Mst2*s2t*pow3(Mt) + 780*Mt*pow3(Mst2)*pow3(s2t) -
        440*pow4(Mt) + 115*pow4(Mst2)*pow4(s2t)) + Mst2*pow4(Mst1)*(2696*pow2(
        Mst2)*pow2(Mt)*pow2(s2t) - 7808*Mst2*s2t*pow3(Mt) + 3200*Mt*pow3(Mst2)*
        pow3(s2t) - 6400*pow4(Mt) + 409*pow4(Mst2)*pow4(s2t)) - 4*Dmglst2*(
        pow4(Mst1)*(-2896*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 2016*Mst2*s2t*pow3(
        Mt) - 672*Mt*pow3(Mst2)*pow3(s2t) + 7424*pow4(Mt) - 33*pow4(Mst2)*pow4(
        s2t)) - 2*pow2(Mst1)*pow2(Mst2)*(1156*pow2(Mst2)*pow2(Mt)*pow2(s2t) +
        384*Mst2*s2t*pow3(Mt) + 390*Mt*pow3(Mst2)*pow3(s2t) - 2632*pow4(Mt) +
        29*pow4(Mst2)*pow4(s2t)) - 3*pow4(Mst2)*(256*pow2(Mst2)*pow2(Mt)*pow2(
        s2t) - 784*Mst2*s2t*pow3(Mt) + 516*Mt*pow3(Mst2)*pow3(s2t) - 512*pow4(
        Mt) + 107*pow4(Mst2)*pow4(s2t))) + 3*(56*pow2(Mst2)*pow2(Mt)*pow2(s2t)
        - 3136*Mst2*s2t*pow3(Mt) + 1040*Mt*pow3(Mst2)*pow3(s2t) + 32*pow4(Mt) +
        271*pow4(Mst2)*pow4(s2t))*pow5(Mst2) + 6*pow2(Dmglst2)*(-512*pow3(Mst2)
        *pow4(Mt) + 2*Mst2*pow2(Mst1)*(380*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 512*
        Mst2*s2t*pow3(Mt) + 1608*pow4(Mt) + 13*pow4(Mst2)*pow4(s2t)) + 256*
        pow2(Mt)*pow2(s2t)*pow5(Mst2) + 768*Mt*pow3(s2t)*pow6(Mst2) + 107*pow4(
        s2t)*pow7(Mst2))) + 3*Mst2*(40*Dmsqst2*pow2(Msq)*pow2(Mst1)*(3087*pow2(
        Mst2)*pow4(Mst1)*(-8*pow2(Mst2)*pow2(Mt)*pow2(s2t)*(36449 + 4000*OepS2
        - 310500*S2 - 6000*T1ep - 75750*z2 + 21000*z3 - 3000*z4 - 4500*pow2(z2)
        ) - 864000*Mst2*s2t*(-6 + 5*z2 - 4*z3)*pow3(Mt) - 432000*Mt*(1 + z2 -
        z3)*pow3(Mst2)*pow3(s2t) + 96*(-91321 + 54000*z2 - 36000*z3)*pow4(Mt) +
        (52199 + 28000*OepS2 - 3415500*S2 - 42000*T1ep - 557250*z2 + 259500*z3
        - 21000*z4 - 31500*pow2(z2))*pow4(Mst2)*pow4(s2t)) + 1157625*pow2(Mst1)
        *pow4(Mst2)*(-4*pow2(Mst2)*pow2(Mt)*pow2(s2t)*(32*OepS2 - 3*(391 + 972*
        S2 + 16*T1ep + 154*z2 - 232*z3 + 8*z4 + 12*pow2(z2))) - 1152*Mst2*s2t*(
        7 + 2*z2 - 8*z3)*pow3(Mt) - 1152*Mt*(-1 + 2*z2 - 3*z3)*pow3(Mst2)*pow3(
        s2t) + 64*(-103 + 108*z2 - 72*z3)*pow4(Mt) + (-1213 + 32*OepS2 + 4860*
        S2 - 48*T1ep - 462*z2 - 276*z3 - 24*z4 - 36*pow2(z2))*pow4(Mst2)*pow4(
        s2t)) + 2*(-4*pow2(Mst2)*pow2(Mt)*pow2(s2t)*(-22430251 + 5488000*OepS2
        - 324135000*S2 - 8232000*T1ep - 103929000*z2 + 28812000*z3 - 4116000*z4
        - 6174000*pow2(z2)) - 333396000*Mst2*s2t*(-35 + 20*z2 - 16*z3)*pow3(Mt)
        + 333396000*Mt*(-2 + z2)*pow3(Mst2)*pow3(s2t) + 3456*(-6582784 +
        2701125*z2 - 2315250*z3)*pow4(Mt) + (378483467 + 9604000*OepS2 -
        754771500*S2 - 14406000*T1ep - 306899250*z2 + 133770000*z3 - 7203000*z4
        - 10804500*pow2(z2))*pow4(Mst2)*pow4(s2t))*pow6(Mst1) + 41674500*(1 +
        2*z2)*pow2(-4*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow6(Mst2)) + 10*pow2(
        Dmsqst2)*pow2(Mst1)*(10290*pow2(Mst2)*pow4(Mst1)*(-4*pow2(Mst2)*pow2(
        Mt)*pow2(s2t)*(-9928 + 1600*OepS2 - 361800*S2 - 2400*T1ep + 34500*z2 +
        50925*z3 - 1200*z4 - 1800*pow2(z2)) - 518400*Mst2*s2t*(-23 + 24*z2 -
        16*z3)*pow3(Mt) - 259200*Mt*(4 + z2 - 3*z3)*pow3(Mst2)*pow3(s2t) + 864*
        (-23941 + 20700*z2 - 10800*z3)*pow4(Mt) + (449186 + 20800*OepS2 -
        3439800*S2 - 31200*T1ep - 523500*z2 + 169275*z3 - 15600*z4 - 23400*
        pow2(z2))*pow4(Mst2)*pow4(s2t)) + 385875*pow2(Mst1)*pow4(Mst2)*(-8*
        pow2(Mst2)*pow2(Mt)*pow2(s2t)*(-6878 + 128*OepS2 - 16848*S2 - 192*T1ep
        - 696*z2 + 3273*z3 - 96*z4 - 144*pow2(z2)) - 3456*Mst2*s2t*(13 + 28*z2
        - 48*z3)*pow3(Mt) - 6912*Mt*(-5 + 7*z2 - 7*z3)*pow3(Mst2)*pow3(s2t) +
        16*(-11134 + 10800*z2 - 5481*z3)*pow4(Mt) + (-16678 + 256*OepS2 +
        114048*S2 - 384*T1ep - 1392*z2 - 16521*z3 - 192*z4 - 288*pow2(z2))*
        pow4(Mst2)*pow4(s2t)) + (8*pow2(Mst2)*pow2(Mt)*pow2(s2t)*(-1351886621 +
        19208000*OepS2 - 639009000*S2 - 28812000*T1ep + 365552250*z2 +
        225865500*z3 - 14406000*z4 - 21609000*pow2(z2)) - 2667168000*Mst2*s2t*(
        -35 + 20*z2 - 16*z3)*pow3(Mt) + 166698000*Mt*(-101 + 106*z2 - 60*z3)*
        pow3(Mst2)*pow3(s2t) + 27648*(-6582784 + 2701125*z2 - 2315250*z3)*pow4(
        Mt) - 3*(-881319682 + 617400000*S2 + 472311000*z2 - 461892375*z3)*pow4(
        Mst2)*pow4(s2t))*pow6(Mst1) + 166698000*(-1 + 3*z2)*pow2(-4*pow2(Mt) +
        pow2(Mst2)*pow2(s2t))*pow6(Mst2)) + pow4(Msq)*(308700*pow4(Mst1)*pow4(
        Mst2)*(2*pow2(Mst2)*pow2(Mt)*pow2(s2t)*(1013263 + 2880*D3 - 2880*DN -
        25440*OepS2 - 137052*S2 + 38160*T1ep + 35562*z2 + 83160*z3 + 36360*z4 +
        28620*pow2(z2)) + 24*Mst2*s2t*(-60759 + 1120*OepS2 + 25596*S2 - 1680*
        T1ep - 137386*z2 + 115320*z3 - 12360*z4 - 1260*pow2(z2))*pow3(Mt) + 30*
        Mt*(25611 + 5280*B4 - 48*DN - 224*OepS2 - 324*S2 + 336*T1ep - 13614*z2
        + 47720*z3 + 2040*z4 - 6660*pow2(z2))*pow3(Mst2)*pow3(s2t) + 8*(389597
        + 3360*OepS2 - 105948*S2 - 5040*T1ep + 293898*z2 - 740160*z3 - 2520*z4
        - 3780*pow2(z2))*pow4(Mt) - 5*(25289 + 1440*B4 - 144*D3 + 72*DN - 2208*
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
        pow2(s2t))*pow8(Mst2))) - 1260*log(pow2(Mst1)/pow2(Mst2))*(441000*
        Dmglst2*Mst2*pow2(Dmsqst2)*pow4(Mst1)*(-3*s2t*(2*Mst2*s2t*(37 - 60*z2)*
        pow2(Mt) + 72*Mt*pow2(Mst2)*pow2(s2t) + 1056*(-17 + 8*z2)*pow3(Mt) + (
        83 - 60*z2)*pow3(Mst2)*pow3(s2t))*pow4(Mst1) + 4*pow3(Mst2)*(-48*Mst2*
        s2t*(-59 + 36*z2)*pow3(Mt) + 72*Mt*(-5 + 3*z2)*pow3(Mst2)*pow3(s2t) +
        4*(-571 + 360*z2)*pow4(Mt) + 9*(2 - 3*z2)*pow4(Mst2)*pow4(s2t)) - 18*
        Mst2*pow2(Mst1)*(-8*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 128*Mst2*s2t*(-12 +
        7*z2)*pow3(Mt) + 16*Mt*(5 - 3*z2)*pow3(Mst2)*pow3(s2t) - 16*(-137 + 72*
        z2)*pow4(Mt) + (-13 + 4*z2)*pow4(Mst2)*pow4(s2t))) + 441000*Dmglst2*
        Dmsqst2*Mst2*pow2(Msq)*pow4(Mst1)*(-3*s2t*(2*Mst2*s2t*(37 - 60*z2)*
        pow2(Mt) + 72*Mt*pow2(Mst2)*pow2(s2t) + 1056*(-17 + 8*z2)*pow3(Mt) + (
        83 - 60*z2)*pow3(Mst2)*pow3(s2t))*pow4(Mst1) + 4*pow3(Mst2)*(-48*Mst2*
        s2t*(-59 + 36*z2)*pow3(Mt) + 72*Mt*(-5 + 3*z2)*pow3(Mst2)*pow3(s2t) +
        4*(-571 + 360*z2)*pow4(Mt) + 9*(2 - 3*z2)*pow4(Mst2)*pow4(s2t)) - 18*
        Mst2*pow2(Mst1)*(-8*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 128*Mst2*s2t*(-12 +
        7*z2)*pow3(Mt) + 16*Mt*(5 - 3*z2)*pow3(Mst2)*pow3(s2t) - 16*(-137 + 72*
        z2)*pow4(Mt) + (-13 + 4*z2)*pow4(Mst2)*pow4(s2t))) - 588*Mst2*pow2(
        Dmglst2)*pow4(Msq)*(-450*pow2(Mst1)*pow4(Mst2)*(-408*pow2(Mst2)*pow2(
        Mt)*pow2(s2t) - 2048*Mst2*s2t*pow3(Mt) + 512*Mt*pow3(Mst2)*pow3(s2t) +
        1328*pow4(Mt) + 51*pow4(Mst2)*pow4(s2t)) + 2*pow2(Mst2)*pow4(Mst1)*(
        300*(6485 + 1134*S2 - 6306*z2)*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 1600*
        Mst2*s2t*(2182 + 378*S2 - 1665*z2)*pow3(Mt) - 600*Mt*(439 + 252*S2 +
        135*z2)*pow3(Mst2)*pow3(s2t) - 16*(-56837 + 89775*S2 - 162900*z2)*pow4(
        Mt) + 75*(-1388 + 63*S2 - 618*z2)*pow4(Mst2)*pow4(s2t)) + (24*(102713 +
        133650*S2 - 173700*z2)*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 400*Mst2*s2t*(
        11587 + 76104*S2 - 49608*z2)*pow3(Mt) - 200*Mt*(-4958 + 24633*S2 - 90*
        z2)*pow3(Mst2)*pow3(s2t) + (3934288 - 40816800*S2 + 34804800*z2)*pow4(
        Mt) + (332311 + 511200*S2 - 49050*z2)*pow4(Mst2)*pow4(s2t))*pow6(Mst1)
        + 7200*pow2(-4*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow6(Mst2)) - 4*
        Dmglst2*pow4(Msq)*(-14700*pow4(Mst1)*pow4(Mst2)*(12*(-3427 + 1782*z2)*
        pow2(Mst2)*pow2(Mt)*pow2(s2t) + 4*Mst2*s2t*(-8398 + 2835*S2 - 3258*z2)*
        pow3(Mt) - 9*Mt*(1640 + 315*S2 - 1146*z2)*pow3(Mst2)*pow3(s2t) + (26276
        - 9072*S2 + 28224*z2)*pow4(Mt) + 9*(277 + 63*S2 - 140*z2)*pow4(Mst2)*
        pow4(s2t)) - 294*pow2(Mst2)*(-8*(77657 + 202500*S2 - 55800*z2)*pow2(
        Mst2)*pow2(Mt)*pow2(s2t) + 200*Mst2*s2t*(-11408 + 37341*S2 - 17730*z2)*
        pow3(Mt) - 300*Mt*(-503 + 3609*S2 - 1494*z2)*pow3(Mst2)*pow3(s2t) - 48*
        (8581 + 72450*S2 - 137700*z2)*pow4(Mt) + (-81643 + 291600*S2 + 5850*z2)
        *pow4(Mst2)*pow4(s2t))*pow6(Mst1) + 132300*pow2(Mst1)*(-536*pow2(Mst2)*
        pow2(Mt)*pow2(s2t) - 512*Mst2*s2t*pow3(Mt) + 128*Mt*pow3(Mst2)*pow3(
        s2t) + 560*pow4(Mt) + 131*pow4(Mst2)*pow4(s2t))*pow6(Mst2) + (72*(-
        4261807 + 33912900*S2 - 73500*z2)*pow2(Mst2)*pow2(Mt)*pow2(s2t) -
        78400*Mst2*s2t*(-4724 + 115785*S2 - 29484*z2)*pow3(Mt) + 4900*Mt*(18652
        + 140139*S2 - 38826*z2)*pow3(Mst2)*pow3(s2t) + 80*(42300121 + 49233240*
        S2 - 64139040*z2)*pow4(Mt) + (8287903 - 185175900*S2 + 9525600*z2)*
        pow4(Mst2)*pow4(s2t))*pow8(Mst1) + 16934400*pow2(Mt)*(-2*pow2(Mt) +
        pow2(Mst2)*pow2(s2t))*pow8(Mst2)) - 3*Mst2*(10*pow2(Dmsqst2)*pow2(Mst1)
        *(1470*pow2(Mst2)*pow4(Mst1)*(-8*(-443 + 90*S2)*pow2(Mst2)*pow2(Mt)*
        pow2(s2t) - 11520*Mst2*s2t*(-7 + 4*z2)*pow3(Mt) + 1440*Mt*(-4 + 3*z2)*
        pow3(Mst2)*pow3(s2t) + 144*(-667 + 360*z2)*pow4(Mt) + (-13 + 2340*S2 -
        540*z2)*pow4(Mst2)*pow4(s2t)) + 14700*pow2(Mst1)*pow4(Mst2)*(-48*(7 +
        9*S2)*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 144*Mst2*s2t*(-35 + 24*z2)*pow3(
        Mt) + 144*Mt*(-4 + 3*z2)*pow3(Mst2)*pow3(s2t) + 28*(-131 + 72*z2)*pow4(
        Mt) + 3*(35 + 36*S2 - 3*z2)*pow4(Mst2)*pow4(s2t)) + (2*(318121 +
        1234800*S2 + 396900*z2)*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 4233600*Mst2*
        s2t*(-25 + 8*z2)*pow3(Mt) - 132300*Mt*(-113 + 60*z2)*pow3(Mst2)*pow3(
        s2t) + 1152*(-175121 + 44100*z2)*pow4(Mt) + 9*(-281161 + 102900*z2)*
        pow4(Mst2)*pow4(s2t))*pow6(Mst1) - 264600*pow2(-4*pow2(Mt) + pow2(Mst2)
        *pow2(s2t))*pow6(Mst2)) + 40*Dmsqst2*pow2(Msq)*pow2(Mst1)*(22050*pow2(
        Mst1)*pow4(Mst2)*(-4*(16 + 27*S2)*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 96*
        Mst2*s2t*(-5 + 4*z2)*pow3(Mt) + 48*Mt*(-1 + z2)*pow3(Mst2)*pow3(s2t) +
        16*(-23 + 12*z2)*pow4(Mt) + (17 + 27*S2)*pow4(Mst2)*pow4(s2t)) + 441*
        pow2(Mst2)*pow4(Mst1)*(-8*(-517 + 450*S2)*pow2(Mst2)*pow2(Mt)*pow2(s2t)
        - 19200*Mst2*s2t*(-2 + z2)*pow3(Mt) + 2400*Mt*(-1 + z2)*pow3(Mst2)*
        pow3(s2t) + 16*(-3067 + 1200*z2)*pow4(Mt) + (-317 + 3150*S2 - 300*z2)*
        pow4(Mst2)*pow4(s2t)) + (16*(58853 - 44100*S2)*pow2(Mst2)*pow2(Mt)*
        pow2(s2t) - 1058400*Mst2*s2t*(-25 + 8*z2)*pow3(Mt) - 264600*Mt*pow3(
        Mst2)*pow3(s2t) + 288*(-175121 + 44100*z2)*pow4(Mt) + (-621671 +
        308700*S2 + 132300*z2)*pow4(Mst2)*pow4(s2t))*pow6(Mst1) - 66150*pow2(-
        4*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow6(Mst2)) + pow4(Msq)*(9800*pow4(
        Mst1)*pow4(Mst2)*(-6*(-4534 + 4293*S2 + 1830*z2 - 72*z3)*pow2(Mst2)*
        pow2(Mt)*pow2(s2t) + 72*Mst2*s2t*(26 + 189*S2 - 758*z2)*pow3(Mt) - 126*
        Mt*(-208 + 27*S2 - 2*z2)*pow3(Mst2)*pow3(s2t) + 8*(-553 + 1701*S2 +
        4023*z2 + 108*z3)*pow4(Mt) + 3*(1148 + 1863*S2 + 276*z2 - 18*z3)*pow4(
        Mst2)*pow4(s2t)) + 1960*pow2(Mst2)*(-4*(-52322 + 84915*S2 + 3060*z2 +
        540*z3)*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 1080*Mst2*s2t*(42 + 377*S2 -
        330*z2)*pow3(Mt) - 180*Mt*(221 + 219*S2 - 234*z2)*pow3(Mst2)*pow3(s2t)
        + 8*(-22352 + 28215*S2 + 36270*z2)*pow4(Mt) + (-16051 + 86805*S2 -
        8325*z2)*pow4(Mst2)*pow4(s2t))*pow6(Mst1) + 88200*pow2(Mst1)*(-616*
        pow2(Mst2)*pow2(Mt)*pow2(s2t) - 1024*Mst2*s2t*pow3(Mt) + 256*Mt*pow3(
        Mst2)*pow3(s2t) + 976*pow4(Mt) + 125*pow4(Mst2)*pow4(s2t))*pow6(Mst2) +
        (-24*(-19407863 + 50398950*S2 + 2866500*z2)*pow2(Mst2)*pow2(Mt)*pow2(
        s2t) + 940800*Mst2*s2t*(167 + 2055*S2 - 912*z2)*pow3(Mt) - 58800*Mt*(46
        + 1317*S2 - 918*z2)*pow3(Mst2)*pow3(s2t) + 32*(-35585111 + 25754400*S2
        + 32281200*z2)*pow4(Mt) + 3*(-23468297 + 26386500*S2 + 661500*z2 +
        176400*z3)*pow4(Mst2)*pow4(s2t))*pow8(Mst1) + 11289600*pow2(Mt)*(-2*
        pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow8(Mst2)))) + 1587600*pow2(log(pow2(
        Mst1)/pow2(Mst2)))*(-7*Mst2*pow2(Dmglst2)*pow4(Msq)*(-15*pow2(Mst1)*
        pow4(Mst2)*(-1240*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 3072*Mst2*s2t*pow3(
        Mt) + 768*Mt*pow3(Mst2)*pow3(s2t) + 944*pow4(Mt) + 443*pow4(Mst2)*pow4(
        s2t)) + 2*pow2(Mst2)*pow4(Mst1)*(7740*pow2(Mst2)*pow2(Mt)*pow2(s2t) -
        61760*Mst2*s2t*pow3(Mt) + 18480*Mt*pow3(Mst2)*pow3(s2t) - 66208*pow4(
        Mt) + 4455*pow4(Mst2)*pow4(s2t)) + (133776*pow2(Mst2)*pow2(Mt)*pow2(
        s2t) + 398240*Mst2*s2t*pow3(Mt) + 24960*Mt*pow3(Mst2)*pow3(s2t) -
        1267120*pow4(Mt) + 4293*pow4(Mst2)*pow4(s2t))*pow6(Mst1) + 1440*pow2(-
        4*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow6(Mst2)) - 2*Dmglst2*(100800*
        Mst2*pow2(Dmsqst2)*(-(Mst2*Mt) + 3*s2t*pow2(Mst1))*pow3(Mt)*pow6(Mst1)
        + 100800*Dmsqst2*Mst2*pow2(Msq)*(-(Mst2*Mt) + 3*s2t*pow2(Mst1))*pow3(
        Mt)*pow6(Mst1) + pow4(Msq)*(-140*pow4(Mst1)*pow4(Mst2)*(1806*pow2(Mst2)
        *pow2(Mt)*pow2(s2t) + 1592*Mst2*s2t*pow3(Mt) + 1104*Mt*pow3(Mst2)*pow3(
        s2t) - 3808*pow4(Mt) + 213*pow4(Mst2)*pow4(s2t)) + 7*pow2(Mst2)*(-
        12992*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 175520*Mst2*s2t*pow3(Mt) + 3600*
        Mt*pow3(Mst2)*pow3(s2t) + 285744*pow4(Mt) + 4969*pow4(Mst2)*pow4(s2t))*
        pow6(Mst1) - 105*pow2(Mst1)*(-728*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 1536*
        Mst2*s2t*pow3(Mt) + 384*Mt*pow3(Mst2)*pow3(s2t) + 944*pow4(Mt) + 187*
        pow4(Mst2)*pow4(s2t))*pow6(Mst2) + (-420096*pow2(Mst2)*pow2(Mt)*pow2(
        s2t) - 1672160*Mst2*s2t*pow3(Mt) + 68880*Mt*pow3(Mst2)*pow3(s2t) +
        3253120*pow4(Mt) + 1377*pow4(Mst2)*pow4(s2t))*pow8(Mst1) - 26880*pow2(
        Mt)*(-2*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow8(Mst2))) + Mst2*(pow4(Msq)
        *(70*pow4(Mst1)*pow4(Mst2)*(3192*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 14880*
        Mst2*s2t*pow3(Mt) + 5664*Mt*pow3(Mst2)*pow3(s2t) + 7744*pow4(Mt) +
        1113*pow4(Mst2)*pow4(s2t)) + 560*pow2(Mst2)*(479*pow2(Mst2)*pow2(Mt)*
        pow2(s2t) - 2448*Mst2*s2t*pow3(Mt) + 366*Mt*pow3(Mst2)*pow3(s2t) +
        1136*pow4(Mt) + 23*pow4(Mst2)*pow4(s2t))*pow6(Mst1) + 105*pow2(Mst1)*(-
        600*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 1024*Mst2*s2t*pow3(Mt) + 256*Mt*
        pow3(Mst2)*pow3(s2t) + 944*pow4(Mt) + 123*pow4(Mst2)*pow4(s2t))*pow6(
        Mst2) + 4*(73657*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 316400*Mst2*s2t*pow3(
        Mt) + 35000*Mt*pow3(Mst2)*pow3(s2t) + 68678*pow4(Mt) + 11764*pow4(Mst2)
        *pow4(s2t))*pow8(Mst1) + 13440*pow2(Mt)*(-2*pow2(Mt) + pow2(Mst2)*pow2(
        s2t))*pow8(Mst2)) - 90*Dmsqst2*pow2(Msq)*pow4(Mst1)*(2*pow4(Mst1)*(74*
        pow2(Mst2)*pow2(Mt)*pow2(s2t) - 1120*Mst2*s2t*pow3(Mt) + 2736*pow4(Mt)
        - 53*pow4(Mst2)*pow4(s2t)) - 7*pow2(Mst1)*(-48*pow2(Mt)*pow2(s2t)*pow4(
        Mst2) - 144*pow2(Mst2)*pow4(Mt) + 31*pow4(s2t)*pow6(Mst2)) + 140*(-4*
        pow4(Mst2)*pow4(Mt) + pow4(s2t)*pow8(Mst2))) - 90*pow2(Dmsqst2)*pow4(
        Mst1)*(pow4(Mst1)*(-32*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 2240*Mst2*s2t*
        pow3(Mt) + 5472*pow4(Mt) - 31*pow4(Mst2)*pow4(s2t)) - 35*pow2(Mst1)*(-
        4*pow2(Mt)*pow2(s2t)*pow4(Mst2) - 16*pow2(Mst2)*pow4(Mt) + 5*pow4(s2t)*
        pow6(Mst2)) + 140*(-2*pow4(Mst2)*pow4(Mt) + pow4(s2t)*pow8(Mst2))))))/(
        5.00094e8*pow4(Msq)*pow4(Mst1)*pow4(Mt)*pow5(Mst2));

   return result;
}

/**
 * @return returns the susy log^1 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double H6bq22g::calc_coef_at_as2_no_sm_logs_log1() const {

   const double result =
      (4*Mst2*pow2(Dmglst2)*pow4(Msq)*(-225*pow2(Mst1)*pow4(Mst2)*(-664*pow2(
        Mst2)*pow2(Mt)*pow2(s2t) - 2048*Mst2*s2t*pow3(Mt) + 512*Mt*pow3(Mst2)*
        pow3(s2t) + 1840*pow4(Mt) + 83*pow4(Mst2)*pow4(s2t)) + 2*pow2(Mst2)*
        pow4(Mst1)*(-150*(-7121 + 6336*z2)*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 400*
        Mst2*s2t*(-1934 + 1791*z2)*pow3(Mt) - 300*Mt*(-875 + 648*z2)*pow3(Mst2)
        *pow3(s2t) + 8*(71687 + 76050*z2)*pow4(Mt) + 75*(-205 + 129*z2)*pow4(
        Mst2)*pow4(s2t)) + 25*(12*(7903 - 5760*z2)*pow2(Mst2)*pow2(Mt)*pow2(
        s2t) - 96*Mst2*s2t*(-1055 + 2181*z2)*pow3(Mt) + 13992*Mt*pow3(Mst2)*
        pow3(s2t) + 80*(-605 + 4608*z2)*pow4(Mt) + 3*(409 - 258*z2)*pow4(Mst2)*
        pow4(s2t))*pow6(Mst1) + 3600*pow2(-4*pow2(Mt) + pow2(Mst2)*pow2(s2t))*
        pow6(Mst2)) + 450*pow2(log(pow2(Mst1)/pow2(Mst2)))*pow4(Msq)*pow4(Mst1)
        *(8*pow2(Mst1)*pow3(Mst2)*(425*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 704*
        Mst2*s2t*pow3(Mt) + 212*Mt*pow3(Mst2)*pow3(s2t) - 76*pow4(Mt) + 3*pow4(
        Mst2)*pow4(s2t)) + Mst2*pow4(Mst1)*(2688*pow2(Mst2)*pow2(Mt)*pow2(s2t)
        - 5056*Mst2*s2t*pow3(Mt) + 1024*Mt*pow3(Mst2)*pow3(s2t) - 2848*pow4(Mt)
        + 41*pow4(Mst2)*pow4(s2t)) - 8*Dmglst2*(8*Mt*(-157*Mst2*s2t*pow2(Mt) -
        108*Mt*pow2(Mst2)*pow2(s2t) + 398*pow3(Mt) - 16*pow3(Mst2)*pow3(s2t))*
        pow4(Mst1) + pow4(Mst2)*(-320*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 848*Mst2*
        s2t*pow3(Mt) - 532*Mt*pow3(Mst2)*pow3(s2t) + 512*pow4(Mt) - 99*pow4(
        Mst2)*pow4(s2t)) + pow2(Mst1)*pow2(Mst2)*(-768*pow2(Mst2)*pow2(Mt)*
        pow2(s2t) - 448*Mst2*s2t*pow3(Mt) - 212*Mt*pow3(Mst2)*pow3(s2t) + 1896*
        pow4(Mt) + 35*pow4(Mst2)*pow4(s2t))) + (296*pow2(Mst2)*pow2(Mt)*pow2(
        s2t) - 6784*Mst2*s2t*pow3(Mt) + 2208*Mt*pow3(Mst2)*pow3(s2t) - 672*
        pow4(Mt) + 519*pow4(Mst2)*pow4(s2t))*pow5(Mst2) + 4*Mst2*pow2(Dmglst2)*
        (320*pow2(Mt)*pow2(s2t)*pow4(Mst2) - 512*pow2(Mst2)*pow4(Mt) + pow2(
        Mst1)*(768*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 1152*Mst2*s2t*pow3(Mt) +
        3640*pow4(Mt) - 35*pow4(Mst2)*pow4(s2t)) + 768*Mt*pow3(s2t)*pow5(Mst2)
        + 99*pow4(s2t)*pow6(Mst2))) - 100*Dmglst2*(-240*Mst2*pow2(Dmsqst2)*
        pow3(Mt)*pow4(Mst1)*(36*Mst2*(Mt*(33 - 18*z2) + 2*Mst2*s2t*(-12 + 7*z2)
        )*pow2(Mst1) + (Mt*(283 - 180*z2) + 12*Mst2*s2t*(-29 + 18*z2))*pow3(
        Mst2) + 18*s2t*(-83 + 44*z2)*pow4(Mst1)) - 240*Dmsqst2*Mst2*pow2(Msq)*
        pow3(Mt)*pow4(Mst1)*(36*Mst2*(Mt*(33 - 18*z2) + 2*Mst2*s2t*(-12 + 7*z2)
        )*pow2(Mst1) + (Mt*(283 - 180*z2) + 12*Mst2*s2t*(-29 + 18*z2))*pow3(
        Mst2) + 18*s2t*(-83 + 44*z2)*pow4(Mst1)) + pow4(Msq)*(2*pow4(Mst1)*
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
        pow2(Mst2)*pow2(s2t))*pow8(Mst2))) + 25*Mst2*(360*Dmsqst2*pow2(Msq)*
        pow2(Mst1)*(3*pow2(Mst2)*pow4(Mst1)*(8*pow2(Mst2)*pow2(Mt)*pow2(s2t) -
        64*Mst2*s2t*(-2 + z2)*pow3(Mt) + 16*(-9 + 4*z2)*pow4(Mt) + pow4(Mst2)*
        pow4(s2t)) + pow2(Mst1)*pow4(Mst2)*(-48*pow2(Mst2)*pow2(Mt)*pow2(s2t) -
        96*Mst2*s2t*(-3 + 2*z2)*pow3(Mt) + 32*(-8 + 3*z2)*pow4(Mt) + 3*pow4(
        Mst2)*pow4(s2t)) + 3*(-16*Mst2*s2t*(-9 + 4*z2)*pow3(Mt) + 8*(-31 + 12*
        z2)*pow4(Mt) - pow4(Mst2)*pow4(s2t))*pow6(Mst1) - 3*pow2(-4*pow2(Mt) +
        pow2(Mst2)*pow2(s2t))*pow6(Mst2)) + 120*pow2(Dmsqst2)*(pow4(Mst1)*pow4(
        Mst2)*(-144*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 72*Mst2*s2t*(-19 + 12*z2)*
        pow3(Mt) + 2*(-545 + 252*z2)*pow4(Mt) + 9*pow4(Mst2)*pow4(s2t)) + 9*
        pow2(Mst2)*(8*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 32*Mst2*s2t*(-7 + 4*z2)*
        pow3(Mt) + 16*(-16 + 9*z2)*pow4(Mt) + pow4(Mst2)*pow4(s2t))*pow6(Mst1)
        - 9*pow2(Mst1)*pow2(-4*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow6(Mst2) + 9*
        (-16*Mst2*s2t*(-9 + 4*z2)*pow3(Mt) + 8*(-31 + 12*z2)*pow4(Mt) - pow4(
        Mst2)*pow4(s2t))*pow8(Mst1)) + pow4(Msq)*(4*pow4(Mst1)*pow4(Mst2)*(-72*
        (-445 + 96*z2)*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 144*Mst2*s2t*(-121 +
        296*z2)*pow3(Mt) - 144*Mt*(-151 + 19*z2)*pow3(Mst2)*pow3(s2t) + 8*(-
        1993 + 2151*z2 + 216*z3)*pow4(Mt) + 81*(8 + 5*z2)*pow4(Mst2)*pow4(s2t))
        - 36*pow2(Mst2)*(-696*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 32*Mst2*s2t*(-179
        + 128*z2)*pow3(Mt) + 16*Mt*(148 - 17*z2)*pow3(Mst2)*pow3(s2t) - 8*(-403
        + 300*z2)*pow4(Mt) + (551 + 86*z2)*pow4(Mst2)*pow4(s2t))*pow6(Mst1) +
        36*pow2(Mst1)*(-744*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 1024*Mst2*s2t*pow3(
        Mt) + 256*Mt*pow3(Mst2)*pow3(s2t) + 1232*pow4(Mt) + 141*pow4(Mst2)*
        pow4(s2t))*pow6(Mst2) + 9*(864*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 512*
        Mst2*s2t*(-49 + 31*z2)*pow3(Mt) - 784*Mt*pow3(Mst2)*pow3(s2t) + 8*(-
        2503 + 1408*z2)*pow4(Mt) + (1547 + 164*z2)*pow4(Mst2)*pow4(s2t))*pow8(
        Mst1) + 4608*pow2(Mt)*(-2*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow8(Mst2)))
        + 30*log(pow2(Mst1)/pow2(Mst2))*(2*Mst2*pow2(Dmglst2)*pow4(Msq)*(2*
        pow2(Mst2)*pow4(Mst1)*(3240*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 45920*Mst2*
        s2t*pow3(Mt) - 9720*Mt*pow3(Mst2)*pow3(s2t) + 37408*pow4(Mt) - 4635*
        pow4(Mst2)*pow4(s2t)) + 45*pow2(Mst1)*pow4(Mst2)*(-456*pow2(Mst2)*pow2(
        Mt)*pow2(s2t) - 1024*Mst2*s2t*pow3(Mt) + 256*Mt*pow3(Mst2)*pow3(s2t) +
        400*pow4(Mt) + 153*pow4(Mst2)*pow4(s2t)) + 2*(-34080*pow2(Mst2)*pow2(
        Mt)*pow2(s2t) - 147560*Mst2*s2t*pow3(Mt) - 5880*Mt*pow3(Mst2)*pow3(s2t)
        + 415932*pow4(Mt) + 1005*pow4(Mst2)*pow4(s2t))*pow6(Mst1) - 1440*pow2(-
        4*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow6(Mst2)) - 20*Dmglst2*(1440*Mst2*
        pow2(Dmsqst2)*(-(Mst2*Mt) + 3*s2t*pow2(Mst1))*pow3(Mt)*pow6(Mst1) +
        1440*Dmsqst2*Mst2*pow2(Msq)*(-(Mst2*Mt) + 3*s2t*pow2(Mst1))*pow3(Mt)*
        pow6(Mst1) + pow4(Msq)*(-4*pow4(Mst1)*pow4(Mst2)*(1860*pow2(Mst2)*pow2(
        Mt)*pow2(s2t) + 1472*Mst2*s2t*pow3(Mt) + 774*Mt*pow3(Mst2)*pow3(s2t) -
        2464*pow4(Mt) + 15*pow4(Mst2)*pow4(s2t)) + 2*pow2(Mst2)*(-432*pow2(
        Mst2)*pow2(Mt)*pow2(s2t) - 11656*Mst2*s2t*pow3(Mt) + 996*Mt*pow3(Mst2)*
        pow3(s2t) + 18748*pow4(Mt) + 207*pow4(Mst2)*pow4(s2t))*pow6(Mst1) - 3*
        pow2(Mst1)*(-856*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 1536*Mst2*s2t*pow3(Mt)
        + 384*Mt*pow3(Mst2)*pow3(s2t) + 1200*pow4(Mt) + 203*pow4(Mst2)*pow4(
        s2t))*pow6(Mst2) + 4*(-738*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 10004*Mst2*
        s2t*pow3(Mt) + 135*Mt*pow3(Mst2)*pow3(s2t) + 17532*pow4(Mt) + 12*pow4(
        Mst2)*pow4(s2t))*pow8(Mst1) - 768*pow2(Mt)*(-2*pow2(Mt) + pow2(Mst2)*
        pow2(s2t))*pow8(Mst2))) + 5*Mst2*(pow4(Msq)*(4*pow4(Mst1)*pow4(Mst2)*(
        3372*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 9024*Mst2*s2t*pow3(Mt) + 3912*Mt*
        pow3(Mst2)*pow3(s2t) + 6112*pow4(Mt) + 579*pow4(Mst2)*pow4(s2t)) - 24*
        pow2(Mst2)*(-534*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 1824*Mst2*s2t*pow3(Mt)
        - 60*Mt*pow3(Mst2)*pow3(s2t) - 934*pow4(Mt) + 97*pow4(Mst2)*pow4(s2t))*
        pow6(Mst1) + 6*pow2(Mst1)*(-728*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 1024*
        Mst2*s2t*pow3(Mt) + 256*Mt*pow3(Mst2)*pow3(s2t) + 1200*pow4(Mt) + 139*
        pow4(Mst2)*pow4(s2t))*pow6(Mst2) + 3*(4640*pow2(Mst2)*pow2(Mt)*pow2(
        s2t) - 15296*Mst2*s2t*pow3(Mt) + 240*Mt*pow3(Mst2)*pow3(s2t) + 7128*
        pow4(Mt) - 279*pow4(Mst2)*pow4(s2t))*pow8(Mst1) + 768*pow2(Mt)*(-2*
        pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow8(Mst2)) - 360*pow2(Dmsqst2)*pow4(
        Mst1)*(8*(5*Mt - 2*Mst2*s2t)*pow3(Mt)*pow4(Mst1) - 4*pow4(Mst2)*pow4(
        Mt) + pow2(Mst1)*(4*pow2(Mst2)*pow4(Mt) - pow4(s2t)*pow6(Mst2)) + pow4(
        s2t)*pow8(Mst2)) - 360*Dmsqst2*pow2(Msq)*pow4(Mst1)*(8*(5*Mt - 2*Mst2*
        s2t)*pow3(Mt)*pow4(Mst1) - 8*pow4(Mst2)*pow4(Mt) + pow2(Mst1)*(8*pow2(
        Mst2)*pow4(Mt) - pow4(s2t)*pow6(Mst2)) + pow4(s2t)*pow8(Mst2)))))/(
        1350.*pow4(Msq)*pow4(Mst1)*pow4(Mt)*pow5(Mst2));

   return result;
}

/**
 * @return returns the susy log^2 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double H6bq22g::calc_coef_at_as2_no_sm_logs_log2() const {

   const double result =
      (Mst2*pow2(Dmglst2)*pow4(Msq)*(pow2(Mst2)*pow4(Mst1)*(33360*pow2(Mst2)*
        pow2(Mt)*pow2(s2t) + 52480*Mst2*s2t*pow3(Mt) + 18896*pow4(Mt) - 10005*
        pow4(Mst2)*pow4(s2t)) + 15*pow2(Mst1)*pow4(Mst2)*(-1496*pow2(Mst2)*
        pow2(Mt)*pow2(s2t) - 3072*Mst2*s2t*pow3(Mt) + 768*Mt*pow3(Mst2)*pow3(
        s2t) + 1456*pow4(Mt) + 475*pow4(Mst2)*pow4(s2t)) + 5*(120*pow2(Mst2)*
        pow2(Mt)*pow2(s2t) - 22656*Mst2*s2t*pow3(Mt) - 2304*Mt*pow3(Mst2)*pow3(
        s2t) + 52384*pow4(Mt) + 879*pow4(Mst2)*pow4(s2t))*pow6(Mst1) - 1440*
        pow2(-4*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow6(Mst2)) + 30*log(pow2(
        Mst1)/pow2(Mst2))*pow4(Msq)*pow4(Mst1)*(pow2(Mst1)*pow3(Mst2)*(1096*
        pow2(Mst2)*pow2(Mt)*pow2(s2t) - 1280*Mst2*s2t*pow3(Mt) + 328*Mt*pow3(
        Mst2)*pow3(s2t) - 32*pow4(Mt) - 91*pow4(Mst2)*pow4(s2t)) - Mst2*pow4(
        Mst1)*(-768*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 1216*Mst2*s2t*pow3(Mt) +
        400*pow4(Mt) + 41*pow4(Mst2)*pow4(s2t)) - 2*Dmglst2*(32*pow2(Mt)*(-57*
        Mst2*Mt*s2t + 119*pow2(Mt) - 24*pow2(Mst2)*pow2(s2t))*pow4(Mst1) +
        pow4(Mst2)*(-384*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 912*Mst2*s2t*pow3(Mt)
        - 548*Mt*pow3(Mst2)*pow3(s2t) + 512*pow4(Mt) - 91*pow4(Mst2)*pow4(s2t))
        + pow2(Mst1)*pow2(Mst2)*(-768*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 640*Mst2*
        s2t*pow3(Mt) - 164*Mt*pow3(Mst2)*pow3(s2t) + 2016*pow4(Mt) + 91*pow4(
        Mst2)*pow4(s2t))) + 4*(14*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 456*Mst2*s2t*
        pow3(Mt) + 146*Mt*pow3(Mst2)*pow3(s2t) - 166*pow4(Mt) + 33*pow4(Mst2)*
        pow4(s2t))*pow5(Mst2) + Mst2*pow2(Dmglst2)*(384*pow2(Mt)*pow2(s2t)*
        pow4(Mst2) - 512*pow2(Mst2)*pow4(Mt) + pow2(Mst1)*(768*pow2(Mst2)*pow2(
        Mt)*pow2(s2t) - 1280*Mst2*s2t*pow3(Mt) + 4000*pow4(Mt) - 91*pow4(Mst2)*
        pow4(s2t)) + 768*Mt*pow3(s2t)*pow5(Mst2) + 91*pow4(s2t)*pow6(Mst2))) -
        10*Dmglst2*pow4(Msq)*(pow4(Mst1)*pow4(Mst2)*(-8208*pow2(Mst2)*pow2(Mt)*
        pow2(s2t) - 4640*Mst2*s2t*pow3(Mt) - 1968*Mt*pow3(Mst2)*pow3(s2t) +
        4816*pow4(Mt) + 849*pow4(Mst2)*pow4(s2t)) - 3*pow2(Mst2)*(-472*pow2(
        Mst2)*pow2(Mt)*pow2(s2t) + 2496*Mst2*s2t*pow3(Mt) - 1040*Mt*pow3(Mst2)*
        pow3(s2t) - 3360*pow4(Mt) + 37*pow4(Mst2)*pow4(s2t))*pow6(Mst1) - 3*
        pow2(Mst1)*(-984*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 1536*Mst2*s2t*pow3(Mt)
        + 384*Mt*pow3(Mst2)*pow3(s2t) + 1456*pow4(Mt) + 219*pow4(Mst2)*pow4(
        s2t))*pow6(Mst2) + 3*(-2240*Mst2*s2t*pow3(Mt) + 3744*pow4(Mt) - 59*
        pow4(Mst2)*pow4(s2t))*pow8(Mst1) - 768*pow2(Mt)*(-2*pow2(Mt) + pow2(
        Mst2)*pow2(s2t))*pow8(Mst2)) + 5*Mst2*(720*pow2(Dmsqst2)*pow4(Mst1)*
        pow4(Mst2)*pow4(Mt) + 1440*Dmsqst2*pow2(Msq)*pow4(Mst1)*pow4(Mst2)*
        pow4(Mt) + pow4(Msq)*(pow4(Mst1)*pow4(Mst2)*(8592*pow2(Mst2)*pow2(Mt)*
        pow2(s2t) - 4800*Mst2*s2t*pow3(Mt) + 3936*Mt*pow3(Mst2)*pow3(s2t) +
        4256*pow4(Mt) - 69*pow4(Mst2)*pow4(s2t)) - 3*pow2(Mst2)*(600*pow2(Mst2)
        *pow2(Mt)*pow2(s2t) - 128*Mst2*s2t*pow3(Mt) + 1568*Mt*pow3(Mst2)*pow3(
        s2t) + 672*pow4(Mt) + 355*pow4(Mst2)*pow4(s2t))*pow6(Mst1) + 3*pow2(
        Mst1)*(-856*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 1024*Mst2*s2t*pow3(Mt) +
        256*Mt*pow3(Mst2)*pow3(s2t) + 1456*pow4(Mt) + 155*pow4(Mst2)*pow4(s2t))
        *pow6(Mst2) + (384*Mst2*s2t*pow3(Mt) - 2400*pow4(Mt) + 717*pow4(Mst2)*
        pow4(s2t))*pow8(Mst1) + 384*pow2(Mt)*(-2*pow2(Mt) + pow2(Mst2)*pow2(
        s2t))*pow8(Mst2))))/(45.*pow4(Msq)*pow4(Mst1)*pow4(Mt)*pow5(Mst2));

   return result;
}

/**
 * @return returns the susy log^3 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double H6bq22g::calc_coef_at_as2_no_sm_logs_log3() const {

   const double result =
      -298.6666666666667;

   return result;
}

}        // hierarchies
}        // himalaya
