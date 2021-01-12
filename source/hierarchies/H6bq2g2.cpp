// ====================================================================
// This file is part of Himalaya.
//
// Himalaya is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================

#include "H6bq2g2.hpp"
#include "Enums.hpp"
#include "Constants.hpp"
#include "Powers.hpp"
#include <cmath>

namespace himalaya{
namespace hierarchies{

/**
 * Constructor
 * @param expansionDepth the flagMap for the truncation of expansion variables
 * @param Al4p a double alpha_s/4/Pi
 * @param beta a double which is the mixing angle beta
 * @param Dmglst2 a double Mgl - Mst2
 * @param Dmsqst2 a double Msq - Mst2
 * @param lmMt a double log((<renormalization scale> / Mt)^2)
 * @param lmMst1 a double log((<renormalization scale> / Mst1)^2)
 * @param lmMst2 a double log((<renormalization scale> / Mst2)^2)
 * @param Mgl a double gluino mass
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
H6bq2g2::H6bq2g2(const ExpansionFlags_t& expansionDepth, double Al4p, double beta, double Dmglst2,
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
   xDmglst2 = expansionDepth.at(ExpansionDepth::Dmglst2);
   xDmsqst2 = expansionDepth.at(ExpansionDepth::Dmsqst2);
   xMst = expansionDepth.at(ExpansionDepth::Mst);
}

/**
 * @return The diagonal (1, 1) matrix element of the Higgs mass matrix as a double for the hierarchy 'H6bq2g2'
 */
double H6bq2g2::getS1() const {
   return -(Al4p*Mt*z2*((Mt*xDmglst2*pow2(Dmglst2)*(9720*s2t*twoLoopFlag*pow2(
        MuSUSY)*(-3*(-14*Mt + Mst2*s2t)*pow2(Mst1)*pow2(Mst2) + (62*Mt - 3*
        Mst2*s2t)*pow4(Mst1) + (22*Mt - 3*Mst2*s2t)*pow4(Mst2)) + (Al4p*
        threeLoopFlag*(Mt*pow2(MuSUSY)*(388800*Dmsqst2*s2t*pow2(Mst2)*((1 - 9*
        lmMst1 + 9*lmMst2)*pow2(Mst1) + (-2 - 3*lmMst1 + 3*lmMst2)*pow2(Mst2))
        + pow2(Msq)*(6*Mst2*(6*(122051 + 1800*lmMst1 + 274680*lmMst2)*Mt + (
        1263161 - 156600*lmMst1 + 979560*lmMst2)*Mst2*s2t)*pow2(Mst1) + 9*(192*
        (617 + 2040*lmMst2)*Mt + (897737 + 47880*lmMst1 + 238680*lmMst2)*Mst2*
        s2t)*pow3(Mst2) + 2*(15057833 - 4014360*lmMst1 + 5563080*lmMst2)*s2t*
        pow4(Mst1))) - (Mst2*pow2(s2t)*(-24300*Dmsqst2*pow2(Mst2)*(4*pow2(
        MuSUSY)*(3*(2 + lmMst1 - lmMst2)*pow2(Mst1)*pow2(Mst2) + (3 + 8*lmMst1
        - 8*lmMst2)*pow4(Mst1) + pow4(Mst2)) - (-93 + 10*lmMst1 - 10*lmMst2)*(-
        1 + pow2(Sbeta))*pow2(Sbeta)*pow6(Mst1)) + pow2(Msq)*pow2(MuSUSY)*(30*(
        21223 + 13230*lmMst1 + 702*lmMst2)*pow2(Mst2)*pow4(Mst1) + 9*(88931 +
        33120*lmMst1 + 13320*lmMst2)*pow2(Mst1)*pow4(Mst2) + (8513956 -
        1231200*lmMst1 + 1509840*lmMst2)*pow6(Mst1) + 38880*pow6(Mst2))))/pow2(
        Mst1)))/pow2(Msq)))/(pow2(Mgl)*pow5(Mst2)) + (4860*Mt*s2t*twoLoopFlag*
        pow2(MuSUSY)*(2*(82*Mt - 3*Mst2*s2t)*xDmglst2*xMst*pow2(Dmglst2)*pow6(
        Mst1) + pow2(Mgl)*(4*Mt*pow2(Mst2)*pow4(Mst1) + 4*Mt*pow2(Mst1)*pow4(
        Mst2) + 4*Mt*xMst*pow6(Mst1) + (4*Mt - Mst2*s2t)*pow6(Mst2)) - 4*
        Dmglst2*Mgl*((-13*Mt + Mst2*s2t)*pow2(Mst2)*pow4(Mst1) + (-9*Mt + Mst2*
        s2t)*pow2(Mst1)*pow4(Mst2) + (-17*Mt*xMst + Mst2*s2t*xMst)*pow6(Mst1) +
        (-5*Mt + Mst2*s2t)*pow6(Mst2))))/(pow2(Mgl)*pow7(Mst2)) - Al4p*
        threeLoopFlag*(7290*Mt*pow2(s2t)*((10*Dmglst2*Dmsqst2*(-93 + 10*lmMst1
        - 10*lmMst2)*(-1 + pow2(Sbeta))*pow2(Sbeta)*pow4(Mst1))/(3.*Mgl*pow2(
        Msq)*pow2(Mst2)) - pow2(MuSUSY)*(103.64814814814815 + (94*lmMst1)/9. -
        (184*lmMst2)/9. + (Dmglst2*(24*lmMst1 - (7*(233 + 240*lmMst2))/27.) + (
        5*Dmsqst2*(72*Dmglst2*(2 + lmMst1 - lmMst2) + 53*Mgl))/(9.*pow2(Msq)))/
        Mgl + ((10800*Dmglst2*Dmsqst2*(3 + 8*lmMst1 - 8*lmMst2) + 150*Dmsqst2*(
        1097 - 72*lmMst1 + 72*lmMst2)*Mgl + 36*Dmglst2*(-6803 + 1810*lmMst1 -
        2670*lmMst2)*pow2(Msq) + (533629 - 900*lmMst1 + 180*lmMst2)*Mgl*pow2(
        Msq))*pow2(Mst1))/(810.*Mgl*pow2(Msq)*pow2(Mst2)) - (2*(19 + (-28*
        Dmglst2 + (30*Dmsqst2*(-2*Dmglst2 + Mgl))/pow2(Msq))/Mgl)*pow2(Mst2))/(
        9.*pow2(Mst1)) + ((19947459 - 244620*lmMst1 + 231660*lmMst2 + (8*
        Dmglst2*(-2128489 + 307800*lmMst1 - 377460*lmMst2))/Mgl + (300*Dmsqst2*
        (20701 - 648*lmMst1 + 648*lmMst2))/pow2(Msq))*pow4(Mst1))/(14580.*pow4(
        Mst2)) - (20*shiftst1*(1 + Dmsqst2/pow2(Msq))*(1 + pow2(Mst2)/pow2(
        Mst1) - (2*(lmMst1 - lmMst2)*(pow2(Mst1)*pow2(Mst2) + pow4(Mst1) +
        pow4(Mst2)))/pow4(Mst2)))/3. - (2*shiftst3*(2*(1 - 2*lmMst1 + 2*lmMst2)
        *pow2(Mst2)*pow4(Mst1) + 2*(1 - lmMst1 + lmMst2)*pow2(Mst1)*pow4(Mst2)
        + (2 - 6*lmMst1 + 6*lmMst2)*pow6(Mst1) + pow6(Mst2)))/(3.*pow2(Mst1)*
        pow4(Mst2)))) + ((-75*s2t*xDmsqst2*pow2(Dmsqst2)*(324*Mst2*Mt*xDmglst2*
        pow2(Dmglst2)*(pow2(MuSUSY)*(-4*(4*(2 + 3*lmMst1 - 3*lmMst2)*Mt - 3*(2
        + lmMst1 - lmMst2)*Mst2*s2t)*pow2(Mst1)*pow2(Mst2) + 4*(4*(1 - 9*lmMst1
        + 9*lmMst2)*Mt + (3 + 8*lmMst1 - 8*lmMst2)*Mst2*s2t)*pow4(Mst1) + 4*
        s2t*pow5(Mst2)) - (-93 + 10*lmMst1 - 10*lmMst2)*Mst2*s2t*(-1 + pow2(
        Sbeta))*pow2(Sbeta)*pow6(Mst1)) + Mgl*(s2t*(81*pow2(Mst2)*pow2(Sbeta)*(
        -2*(11 + 10*lmMst1 - 10*lmMst2)*Mgl*Mst2*s2t*(-1 + pow2(Sbeta)) + (
        Dmglst2*(372 - 40*lmMst1 + 40*lmMst2) + 3*(-17 + 2*lmMst1 - 2*lmMst2)*
        Mgl)*Mt*pow2(Sbeta))*pow6(Mst1) + Mt*(81*(4*Dmglst2*(-93 + 10*lmMst1 -
        10*lmMst2) + 3*(17 - 2*lmMst1 + 2*lmMst2)*Mgl)*pow2(Mst2)*pow2(Sbeta)*
        pow6(Mst1) + 2*pow2(MuSUSY)*(648*Dmglst2*pow2(Mst2)*(3*(2 + lmMst1 -
        lmMst2)*pow2(Mst1)*pow2(Mst2) + (3 + 8*lmMst1 - 8*lmMst2)*pow4(Mst1) +
        pow4(Mst2)) + Mgl*(3*(1607 - 432*lmMst1 + 432*lmMst2)*pow2(Mst2)*pow4(
        Mst1) - 18*(25 + 9*lmMst1 - 9*lmMst2)*pow2(Mst1)*pow4(Mst2) + (20701 -
        648*lmMst1 + 648*lmMst2)*pow6(Mst1) - 486*pow6(Mst2))))) + 648*Mt*pow2(
        MuSUSY)*(8*Dmglst2*Mst2*Mt*pow2(Mst1)*((1 - 9*lmMst1 + 9*lmMst2)*pow2(
        Mst1) + (-2 - 3*lmMst1 + 3*lmMst2)*pow2(Mst2)) + Mgl*((4*(7 + 3*lmMst1
        - 3*lmMst2)*Mt + (-1 + 2*lmMst1 - 2*lmMst2)*Mst2*s2t*shiftst1)*pow2(
        Mst1)*pow3(Mst2) + 2*Mst2*(4*(5 + 3*lmMst1 - 3*lmMst2)*Mt + (lmMst1 -
        lmMst2)*Mst2*s2t*shiftst1)*pow4(Mst1) + 2*(lmMst1 - lmMst2)*s2t*
        shiftst1*pow6(Mst1) - s2t*shiftst1*pow6(Mst2))))))/(pow2(Mst1)*pow4(
        Msq)*pow4(Mst2)) + (Mt*pow2(MuSUSY)*(Mgl*Mt*(540*(Dmglst2*(-7286 + 240*
        lmMst1 - 6384*lmMst2) + (353 - 72*lmMst1 - 696*lmMst2)*Mgl)*Mt*pow2(
        Mst1)*pow2(Mst2) - (388800*Dmsqst2*s2t*(Dmglst2*(1 - 9*lmMst1 + 9*
        lmMst2)*pow2(Mst1) + (5 + 3*lmMst1 - 3*lmMst2)*Mgl*pow2(Mst1) -
        Dmglst2*(2 + 3*lmMst1 - 3*lmMst2)*pow2(Mst2) + (2 + lmMst1 - lmMst2)*
        Mgl*pow2(Mst2))*pow3(Mst2))/pow2(Msq) + 90*(38401 + 1080*lmMst1 - 7992*
        lmMst2)*Mgl*Mt*pow4(Mst1) - 8640*(24*Dmglst2*(5 + 6*lmMst2) + (53 + 24*
        lmMst2)*Mgl)*Mt*pow4(Mst2) - Mst2*s2t*(3*Mgl*(120*(3937 + 2988*lmMst1 -
        2232*lmMst2)*pow2(Mst1)*pow2(Mst2) + (533726 + 792720*lmMst1 - 693360*
        lmMst2)*pow4(Mst1) + 135*(2269 + 664*lmMst1 - 56*lmMst2)*pow4(Mst2)) +
        Dmglst2*(36*(364291 - 88560*lmMst1 + 147960*lmMst2)*pow2(Mst1)*pow2(
        Mst2) + (30115666 - 8028720*lmMst1 + 11126160*lmMst2)*pow4(Mst1) + 27*(
        173947 - 25080*lmMst1 + 68760*lmMst2)*pow4(Mst2)))) - (6480*Mst2*s2t*
        xDR2DRMOD*(Mgl*(32*(Dmglst2 + 3*Dmglst2*lmMst2 + Mgl + lmMst2*Mgl)*Mt*
        pow2(Mst1)*(2*pow2(Mst1)*pow2(Mst2) + 3*pow4(Mst1) + pow4(Mst2)) -
        Mst2*s2t*(-2*(Dmglst2*(7*lmMst1 - 15*lmMst2) + (-4 - 13*lmMst1 + 9*
        lmMst2)*Mgl)*(pow2(Mst1) + pow2(Mst2))*pow4(Mst1) + (Dmglst2*(7 - 14*
        lmMst1 + 30*lmMst2) + (-5 + 26*lmMst1 - 18*lmMst2)*Mgl)*pow2(Mst1)*
        pow4(Mst2) + (7*Dmglst2 - 13*Mgl)*pow6(Mst2))) - xDmglst2*pow2(Dmglst2)
        *(2*(-32*(-1 + 6*lmMst2)*Mt + (-8 - 27*lmMst1 + 39*lmMst2)*Mst2*s2t)*
        pow2(Mst2)*pow4(Mst1) + (32*(1 - 6*lmMst2)*Mt + (11 - 54*lmMst1 + 78*
        lmMst2)*Mst2*s2t)*pow2(Mst1)*pow4(Mst2) - 2*(48*(1 + 3*lmMst2)*Mt + (7*
        lmMst1 - 15*lmMst2)*Mst2*s2t)*pow6(Mst1) + 27*s2t*pow7(Mst2))))/pow2(
        Mst1)))/pow6(Mst2))/pow2(Mgl))))/7290. - threeLoopFlag*pow2(Al4p)*(
        pow2(Mt)*pow2(s2t)*((5*Dmglst2*Dmsqst2*(-1081 + 165*lmMst1 - 165*
        lmMst2)*(-1 + pow2(Sbeta))*pow2(Sbeta)*pow4(Mst1))/(9.*Mgl*pow2(Msq)*
        pow2(Mst2)) + pow2(MuSUSY)*(53.385802469135804 + (40*B4)/9. - (4*D3)/9.
         + (2*DN)/9. + (1672*lmMst1)/27. + (53*pow2(lmMst1))/9. - lmMst2*(
        129.92592592592592 - 72*lmMst1 + (13*pow2(lmMst1))/3.) + (2*(-470 +
        147*lmMst1)*pow2(lmMst2))/9. + (16*pow3(lmMst1))/9. - (271*pow3(lmMst2)
        )/9. + ((-5*Dmsqst2*(288*Dmglst2*(4 + lmMst1 - lmMst2) + Mgl*(-1141 +
        264*lmMst2 + 48*lmMst1*(-7 + 3*lmMst2) - 144*pow2(lmMst2))))/(54.*pow2(
        Msq)) - (Dmglst2*(14267 - 432*B4 + 576*D3 - 360*DN + 4752*lmMst1 -
        1404*pow2(lmMst1) + 72*lmMst2*(-63 - 232*lmMst1 + 16*pow2(lmMst1)) +
        72*(281 - 123*lmMst1)*pow2(lmMst2) + 7704*pow3(lmMst2)))/162.)/Mgl + (
        pow2(Mst1)*(434.270658436214 + (76*B4)/9. - (2*DN)/9. + (69088*lmMst1)/
        405. - (1313*pow2(lmMst1))/27. - (4*lmMst2*(16192 - 26430*lmMst1 +
        3465*pow2(lmMst1)))/405. + ((-5735 + 3072*lmMst1)*pow2(lmMst2))/27. + (
        Dmsqst2*(201.74098765432097 + (622*lmMst1)/15. + (2*(-311 + 10*lmMst1)*
        lmMst2)/15. + (10*Dmglst2*(23 - 42*lmMst1 + 42*lmMst2))/(3.*Mgl) - (22*
        pow2(lmMst1))/3. + 6*pow2(lmMst2)))/pow2(Msq) - (62*pow3(lmMst1))/27. -
        (2086*pow3(lmMst2))/27. - (2*Dmglst2*(2695042 - 40500*B4 + 54000*D3 -
        33750*DN - 326895*lmMst1 - 324900*pow2(lmMst1) + 15*lmMst2*(-19607 -
        129030*lmMst1 + 62550*pow2(lmMst1)) - 450*(-5023 + 5610*lmMst1)*pow2(
        lmMst2) + 11250*pow3(lmMst1) + 1575000*pow3(lmMst2)))/(30375.*Mgl)))/
        pow2(Mst2) + ((628.1736268201578 + (76*B4)/9. - (2*DN)/9. + (6317839*
        lmMst1)/396900. - (66307*pow2(lmMst1))/315. - lmMst2*(12.52907281431091
         - (182909*lmMst1)/315. + (274*pow2(lmMst1))/3.) + (2*(-58301 + 37135*
        lmMst1)*pow2(lmMst2))/315. + (Dmsqst2*(237.28785508324435 + (16526*
        lmMst2)/3969. + (2*lmMst1*(-8263 + 71820*lmMst2))/3969. - (520*pow2(
        lmMst1))/21. - (80*pow2(lmMst2))/7.))/pow2(Msq) - (44*pow3(lmMst1))/9.
         - (1256*pow3(lmMst2))/
        9. - (Dmglst2*(585.1892843532082 - (8*B4)/3. + (32*D3)/9. - (20*DN)/9.
         - (20109937*lmMst1)/
        297675. - (15886*pow2(lmMst1))/945. + lmMst2*(17.112243218274966 - (
        144628*lmMst1)/945. + (1180*pow2(lmMst1))/9.) + (169.85608465608465 - (
        2632*lmMst1)/9.)*pow2(lmMst2) - (92*pow3(lmMst1))/27. + (4448*pow3(
        lmMst2))/27.))/Mgl)*pow4(Mst1))/pow4(Mst2) + (-((360*Dmglst2*Dmsqst2 +
        Dmsqst2*(30 - 60*lmMst2)*Mgl + Mgl*(103 + 186*lmMst2 + 32*lmMst1*(1 +
        lmMst2) + 91*pow2(lmMst2))*pow2(Msq) + 2*Dmglst2*(-20 + 2*(99 + 16*
        lmMst1)*lmMst2 + 123*pow2(lmMst2))*pow2(Msq))*pow2(Mst2))/(9.*pow2(Msq)
        *pow2(Mst1)) + (16*(1 + lmMst2)*(4*Dmglst2*lmMst2 + Mgl + lmMst2*Mgl)*
        pow4(Mst2))/(9.*pow4(Mst1)) + (-30*Dmsqst2*Mgl*(9*(104*OepS2 + 27*(17 -
        78*lmMst1 + 78*lmMst2)*S2)*pow2(Mst1)*pow2(Mst2) + 13*(136*OepS2 - 27*(
        95 + 102*lmMst1 - 102*lmMst2)*S2)*pow4(Mst1) + 27*(8*OepS2 - 81*(-15 +
        2*lmMst1 - 2*lmMst2)*S2)*pow4(Mst2)) + pow2(Msq)*(3*Mgl*(-3*(8456*OepS2
        - 81*(11243 + 2114*lmMst1 - 2114*lmMst2)*S2)*pow2(Mst1)*pow2(Mst2) + (-
        52948*OepS2 + 27*(194357 + 39711*lmMst1 - 39711*lmMst2)*S2)*pow4(Mst1)
        + 27*(-184*OepS2 + 81*(307 + 46*lmMst1 - 46*lmMst2)*S2)*pow4(Mst2)) +
        2*Dmglst2*(54*(344*OepS2 + 9*(15643 - 774*lmMst1 + 774*lmMst2)*S2)*
        pow2(Mst1)*pow2(Mst2) + 4*(17308*OepS2 + 27*(93919 - 12981*lmMst1 +
        12981*lmMst2)*S2)*pow4(Mst1) + 27*(56*OepS2 - 81*(-1677 + 14*lmMst1 -
        14*lmMst2)*S2)*pow4(Mst2))))/(2187.*pow2(Msq)*pow4(Mst2)))/Mgl + ((10*
        shiftst1*(Dmsqst2*(3 - 2*lmMst2) + (1 - 2*lmMst2)*pow2(Msq))*((pow2(
        Mst1) + pow2(Mst2))*pow4(Mst2) - 2*(lmMst1 - lmMst2)*pow2(Mst1)*(pow2(
        Mst1)*pow2(Mst2) + pow4(Mst1) + pow4(Mst2))))/pow2(Msq) + (1 - 2*
        lmMst2)*shiftst3*(2*(1 - 2*lmMst1 + 2*lmMst2)*pow2(Mst2)*pow4(Mst1) +
        2*(1 - lmMst1 + lmMst2)*pow2(Mst1)*pow4(Mst2) + (2 - 6*lmMst1 + 6*
        lmMst2)*pow6(Mst1) + pow6(Mst2)))/(3.*pow2(Mst1)*pow4(Mst2)))) + (Mt*(
        3*Mst2*pow2(Mst1)*(-1715*z3*pow2(Mst1)*(-15*Mst2*s2t*xDmsqst2*pow2(
        Dmsqst2)*(2592*Dmglst2*Mgl*Mst2*Mt*(pow2(MuSUSY)*(8*(5*Mt - 2*Mst2*s2t)
        *pow2(Mst1) + 8*Mt*pow2(Mst2) - 6*s2t*pow3(Mst2)) + 5*Mst2*s2t*(-1 +
        pow2(Sbeta))*pow2(Sbeta)*pow4(Mst1)) + 2592*Mst2*Mt*xDmglst2*pow2(
        Dmglst2)*(pow2(MuSUSY)*(8*(5*Mt - 2*Mst2*s2t)*pow2(Mst1) + 8*Mt*pow2(
        Mst2) - 6*s2t*pow3(Mst2)) + 5*Mst2*s2t*(-1 + pow2(Sbeta))*pow2(Sbeta)*
        pow4(Mst1)) - pow2(Mgl)*(Mt*pow2(MuSUSY)*(48*Mst2*(3024*Mt - 1783*Mst2*
        s2t)*pow2(Mst1) + 9*(8064*Mt - 5507*Mst2*s2t)*pow3(Mst2)) + 8*s2t*(
        3542*Mt*pow2(MuSUSY) + 243*Mt*pow2(Mst2)*(-1 + pow2(Sbeta))*pow2(Sbeta)
        - 810*s2t*(-1 + pow2(Sbeta))*pow2(Sbeta)*pow3(Mst2))*pow4(Mst1))) + Mt*
        (240*Dmsqst2*Mst2*s2t*pow2(Msq)*(-162*Mst2*xDmglst2*pow2(Dmglst2)*(
        pow2(MuSUSY)*(8*(5*Mt - 2*Mst2*s2t)*pow2(Mst1) + 8*Mt*pow2(Mst2) - 6*
        s2t*pow3(Mst2)) + 5*Mst2*s2t*(-1 + pow2(Sbeta))*pow2(Sbeta)*pow4(Mst1))
        + Mgl*s2t*(162*Dmglst2*pow2(Mst2)*(2*(8*pow2(Mst1) + 3*pow2(Mst2))*
        pow2(MuSUSY) - 5*(-1 + pow2(Sbeta))*pow2(Sbeta)*pow4(Mst1)) + Mgl*pow2(
        MuSUSY)*(315*pow2(Mst1)*pow2(Mst2) + 1771*pow4(Mst1) - 621*pow4(Mst2)))
        ) + 4*pow2(MuSUSY)*pow4(Msq)*(xDmglst2*pow2(Dmglst2)*(48*Mst2*pow2(
        Mst1)*(152143*Mst2*Mt*s2t - 40917*pow2(Mt) + 101147*pow2(Mst2)*pow2(
        s2t)) + 9*(443828*Mst2*Mt*s2t - 95520*pow2(Mt) + 228607*pow2(Mst2)*
        pow2(s2t))*pow3(Mst2) + 16*s2t*(395581*Mt + 176558*Mst2*s2t)*pow4(Mst1)
        ) + 2*Mgl*Mst2*pow2(s2t)*(Dmglst2*(1052676*pow2(Mst1)*pow2(Mst2) +
        1412464*pow4(Mst1) + 475605*pow4(Mst2)) + 6*Mgl*(3*(35719 + 108*lmMst1
        - 108*lmMst2)*pow2(Mst1)*pow2(Mst2) + 2*(91963 + 162*lmMst1 - 162*
        lmMst2)*pow4(Mst1) + 27*(1319 + 6*lmMst1 - 6*lmMst2)*pow4(Mst2)))))) +
        2*Mst2*s2t*xDmsqst2*pow2(Dmsqst2)*(-2778300*Mst2*Mt*xDmglst2*pow2(
        Dmglst2)*(pow2(MuSUSY)*(48*(-2*(2 + 5*lmMst1 - 5*lmMst2)*Mt + (4 +
        lmMst1 - lmMst2)*Mst2*s2t)*pow2(Mst1)*pow2(Mst2) - 6*(16*(-8 + 15*
        lmMst1 - 15*lmMst2)*Mt + (23 - 42*lmMst1 + 42*lmMst2)*Mst2*s2t)*pow4(
        Mst1) + 72*s2t*pow5(Mst2)) - (-1081 + 165*lmMst1 - 165*lmMst2)*Mst2*
        s2t*(-1 + pow2(Sbeta))*pow2(Sbeta)*pow6(Mst1)) + Mgl*(-411600*Mt*pow2(
        MuSUSY)*(648*Dmglst2*Mst2*Mt*pow2(Mst1)*((8 - 15*lmMst1 + 15*lmMst2)*
        pow2(Mst1) + (-2 - 5*lmMst1 + 5*lmMst2)*pow2(Mst2)) + Mgl*(9*(36*(5 +
        4*lmMst1 - 4*lmMst2)*Mt - Mst2*s2t*(9*(-1 + 2*lmMst1 - 2*lmMst2)*(-2 +
        lmMst2)*shiftst1 + 4*T1ep))*pow2(Mst1)*pow3(Mst2) + 6*Mst2*(108*(1 + 3*
        lmMst1 - 3*lmMst2)*Mt + Mst2*s2t*(-27*(lmMst1 - lmMst2)*(-2 + lmMst2)*
        shiftst1 - 25*T1ep))*pow4(Mst1) - s2t*(81*(lmMst1 - lmMst2)*(-3 + 2*
        lmMst2)*shiftst1 + 442*T1ep)*pow6(Mst1) + 81*(-2 + lmMst2)*s2t*
        shiftst1*pow6(Mst2))) + s2t*(-694575*pow2(Mst2)*pow2(Sbeta)*(-6*(11 +
        39*lmMst1 - 39*lmMst2)*Mgl*Mst2*s2t*(-1 + pow2(Sbeta)) + (Dmglst2*(4324
        - 660*lmMst1 + 660*lmMst2) + 3*(-167 + 5*lmMst1 - 5*lmMst2)*Mgl)*Mt*
        pow2(Sbeta))*pow6(Mst1) + Mt*(694575*(Dmglst2*(4324 - 660*lmMst1 + 660*
        lmMst2) + 3*(-167 + 5*lmMst1 - 5*lmMst2)*Mgl)*pow2(Mst2)*pow2(Sbeta)*
        pow6(Mst1) + pow2(MuSUSY)*(16669800*Dmglst2*pow2(Mst2)*(-8*(4 + lmMst1
        - lmMst2)*pow2(Mst1)*pow2(Mst2) + (23 - 42*lmMst1 + 42*lmMst2)*pow4(
        Mst1) - 12*pow4(Mst2)) + Mgl*(-4116*(45*lmMst1*(-1547 + 180*lmMst2 -
        4500*S2) + 45*lmMst2*(1547 + 4500*S2) + 2*(-106283 + 5000*OepS2 +
        639225*S2) + 4050*pow2(lmMst1) - 12150*pow2(lmMst2))*pow2(Mst2)*pow4(
        Mst1) + 77175*(8771 - 128*OepS2 - 57024*S2 - 72*lmMst2*(23 + 36*S2) +
        72*lmMst1*(29 - 12*lmMst2 + 36*S2) + 864*pow2(lmMst2))*pow2(Mst1)*pow4(
        Mst2) + 2*(593331163 - 60642400*OepS2 + 1260*lmMst2*(8263 - 974610*S2)
        + 1143733500*S2 + 1260*lmMst1*(-8263 + 71820*lmMst2 + 974610*S2) -
        61916400*pow2(lmMst1) - 28576800*pow2(lmMst2))*pow6(Mst1) + 16669800*(1
        + 2*lmMst2)*pow6(Mst2)))))))) - 686*Mt*pow2(MuSUSY)*(30*pow4(Mst1)*(2*
        Mst2*xDmglst2*pow2(Dmglst2)*pow4(Msq)*(6*Mst2*pow2(Mst1)*(Mst2*Mt*s2t*(
        12344*z4 - 88647*pow2(z2)) + pow2(Mt)*(5412*z4 - 630*pow2(z2)) + pow2(
        Mst2)*pow2(s2t)*(10819*z4 + 3471*pow2(z2))) + 9*(1944*z4*pow2(Mt) +
        Mst2*Mt*s2t*(2498*z4 - 31083*pow2(z2)) + pow2(Mst2)*pow2(s2t)*(4444*z4
        + 2049*pow2(z2)))*pow3(Mst2) + 2*s2t*(2636*Mt*z4 + 22654*Mst2*s2t*z4 -
        259215*Mt*pow2(z2) + 31794*Mst2*s2t*pow2(z2))*pow4(Mst1)) + 2*Dmglst2*
        Mgl*Mst2*pow2(Msq)*(-38880*Dmsqst2*Mt*s2t*z3*pow2(Mst2)*(5*pow2(Mst1) +
        pow2(Mst2)) + pow2(Msq)*(36*Mst2*pow2(Mst1)*(Mst2*Mt*s2t*(47456*z3 +
        709*z4 - 8535*pow2(z2)) - 6*pow2(Mt)*(250*z3 - 94*z4 + 21*pow2(z2)) +
        3*pow2(Mst2)*pow2(s2t)*(185*z4 + 237*pow2(z2))) + 27*(16*(-83*z3 + 27*
        z4)*pow2(Mt) + Mst2*Mt*s2t*(30274*z3 + 542*z4 - 5289*pow2(z2)) + 2*
        pow2(Mst2)*pow2(s2t)*(212*z4 + 237*pow2(z2)))*pow3(Mst2) + 2*s2t*(
        1582324*Mt*z3 + 2636*Mt*z4 + 22654*Mst2*s2t*z4 - 259215*Mt*pow2(z2) +
        31794*Mst2*s2t*pow2(z2))*pow4(Mst1))) + 3*pow2(Mgl)*(-20*Dmsqst2*s2t*
        pow2(Msq)*pow2(Mst2)*(9*Mst2*pow2(Mst1)*(-1008*Mt*z3 + 26*Mst2*s2t*z4 +
        39*Mst2*s2t*pow2(z2)) + 27*(-144*Mt*z3 + 2*Mst2*s2t*z4 + 3*Mst2*s2t*
        pow2(z2))*pow3(Mst2) + 221*s2t*(2*z4 + 3*pow2(z2))*pow4(Mst1)) - 20*
        xDmsqst2*pow2(Dmsqst2)*pow2(Mst2)*pow2(s2t)*(2*z4 + 3*pow2(z2))*(75*
        pow2(Mst1)*pow2(Mst2) + 221*pow4(Mst1) + 18*pow4(Mst2)) + pow4(Msq)*(-
        6*pow2(Mst1)*pow2(Mst2)*(-36*pow2(Mt)*(218*z3 + 50*z4 + 21*pow2(z2)) +
        4*Mst2*Mt*s2t*(-14636*z3 - 1135*z4 + 2388*pow2(z2)) + pow2(Mst2)*pow2(
        s2t)*(2276*z4 + 2523*pow2(z2))) + (4*Mst2*Mt*s2t*(99632*z3 + 11764*z4 -
        20505*pow2(z2)) + 12*pow2(Mt)*(1502*z3 + 2062*z4 + 1635*pow2(z2)) -
        pow2(Mst2)*pow2(s2t)*(27446*z4 + 35823*pow2(z2)))*pow4(Mst1) + 18*(8*(
        230*z3 + 27*z4)*pow2(Mt) + 5*Mst2*Mt*s2t*(2386*z3 + 102*z4 - 333*pow2(
        z2)) + 3*pow2(Mst2)*pow2(s2t)*(-40*z4 + 3*pow2(z2)))*pow4(Mst2)))) +
        Mgl*pow2(Msq)*(180*pow2(Msq)*pow2(Mt)*pow4(Mst1)*(9*pow2(Mst1)*pow2(
        Mst2)*(Mgl*(10667 - 96*B4 + 96*D3 - 48*DN - 3072*lmMst1 - 384*(-25 + 9*
        lmMst1)*lmMst2 - 384*(-1 + lmMst1 - lmMst2)*lmMt - 384*(-13 + 4*lmMst1)
        *pow2(lmMst2) + 1536*pow3(lmMst2)) + Dmglst2*(18936.666666666668 - 192*
        B4 + 192*D3 - 96*DN + 3648*lmMst1 - 192*(-89 + 2*lmMst1)*lmMst2 + 384*(
        -2 + lmMst1 - lmMst2)*lmMt - 384*(-9 + 8*lmMst1)*pow2(lmMst2) + 3072*
        pow3(lmMst2))) + (Mgl*(383185 - 2592*B4 + 2592*D3 - 1296*DN - 187704*
        lmMst1 + 216*lmMst2*(1733 + 859*lmMst2 - 6*lmMst1*(105 + 41*lmMst2) -
        26*pow2(lmMst1)) - 7992*pow2(lmMst1) + 3456*lmMt*(3 + 6*lmMst2 - 2*
        lmMst1*(3 + lmMst2) + pow2(lmMst1) + pow2(lmMst2)) + 720*pow3(lmMst1) +
        58032*pow3(lmMst2))*pow4(Mst1))/2. + OepS2*(672*(2*Dmglst2 - 3*Mgl)*
        pow2(Mst1)*pow2(Mst2) - 8720*Mgl*pow4(Mst1)) + 72*(6*Dmglst2*(180 - 2*
        B4 + 2*D3 - DN + 16*lmMst1 + 144*lmMst2 - 16*(-2 + lmMst1)*pow2(lmMst2)
        + 16*pow3(lmMst2)) + Mgl*(436 - 6*B4 + 6*D3 - 3*DN - 48*lmMst1 + (408 -
        96*lmMst1)*lmMst2 + 24*lmMt - 48*(-4 + lmMst1)*pow2(lmMst2) + 48*pow3(
        lmMst2)))*pow4(Mst2) - 243*S2*(4*(2*Dmglst2*(65 + 14*lmMst1 - 14*
        lmMst2) + 3*(43 - 14*lmMst1 + 14*lmMst2)*Mgl)*pow2(Mst1)*pow2(Mst2) + (
        4*(57 - 545*lmMst1 + 545*lmMst2)*Mgl*pow4(Mst1))/3. + 96*(4*Dmglst2 +
        3*Mgl)*pow4(Mst2))) + 120*T1ep*pow4(Mst1)*(-(pow2(Mst2)*pow2(s2t)*(60*
        Dmsqst2*Mgl*(117*pow2(Mst1)*pow2(Mst2) + 221*pow4(Mst1) + 27*pow4(Mst2)
        ) + pow2(Msq)*(-4*Dmglst2*(2322*pow2(Mst1)*pow2(Mst2) + 8654*pow4(Mst1)
        + 189*pow4(Mst2)) + 3*Mgl*(6342*pow2(Mst1)*pow2(Mst2) + 13237*pow4(
        Mst1) + 1242*pow4(Mst2))))) + pow2(Msq)*(-2*Dmglst2*Mst2*Mt*(36*Mst2*(
        42*Mt + 253*Mst2*s2t)*pow2(Mst1) + 32842*s2t*pow4(Mst1) + 945*s2t*pow4(
        Mst2)) + 6*Mgl*Mt*(12*(63*Mt + 68*Mst2*s2t)*pow2(Mst1)*pow2(Mst2) + 2*(
        1635*Mt + 941*Mst2*s2t)*pow4(Mst1) + 189*s2t*pow5(Mst2)))) + Mst2*Mt*
        s2t*pow2(Mst1)*(2332800*Dmsqst2*pow2(Mst1)*pow2(Mst2)*(Dmglst2*(8 - 15*
        lmMst1 + 15*lmMst2)*pow2(Mst1) + (1 + 3*lmMst1 - 3*lmMst2)*Mgl*pow2(
        Mst1) - Dmglst2*(2 + 5*lmMst1 - 5*lmMst2)*pow2(Mst2) + (1 + lmMst1 -
        lmMst2)*Mgl*pow2(Mst2)) + pow2(Msq)*(Dmglst2*(36*pow2(Mst2)*(66761 +
        301320*B4 - 4860*DN - 205380*lmMst1 + 40480*OepS2 - 216*(2489 + 3795*
        lmMst1 - 3795*lmMst2)*S2 + 23760*pow2(lmMst1) + 180*lmMst2*(4993 -
        1956*lmMst1 + 48*pow2(lmMst1)) - 1080*(-482 + 331*lmMst1)*pow2(lmMst2)
        + 348840*pow3(lmMst2))*pow4(Mst1) + 27*pow2(Mst1)*(23917 + 188640*B4 -
        3600*DN - 37440*lmMst1 + 5600*OepS2 - 324*(-453 + 350*lmMst1 - 350*
        lmMst2)*S2 + 11520*pow2(lmMst1) - 2880*lmMst2*(-237 + 55*lmMst1 + 4*
        pow2(lmMst1)) - 1440*(-280 + 121*lmMst1)*pow2(lmMst2) + 185760*pow3(
        lmMst2))*pow4(Mst2) - 10*(2773621 - 1660176*B4 + 25272*DN + 2004408*
        lmMst1 - 525472*OepS2 + 108*(123113 + 98526*lmMst1 - 98526*lmMst2)*S2 +
        3888*pow2(lmMst1) - 144*lmMst2*(36802 - 11421*lmMst1 + 1728*pow2(
        lmMst1)) + 167184*(-14 + 15*lmMst1)*pow2(lmMst2) - 31104*pow3(lmMst1) -
        2227824*pow3(lmMst2))*pow6(Mst1) + 622080*(-1 + 2*lmMst2 + 3*pow2(
        lmMst2))*pow6(Mst2)) + 15*Mgl*(24*pow2(Mst2)*(75569 + 13716*B4 - 54*DN
        - 33426*lmMst1 - 1088*OepS2 + 162*(169 + 136*lmMst1 - 136*lmMst2)*S2 -
        2376*pow2(lmMst1) + 54*lmMst2*(1427 - 1012*lmMst1 + 16*pow2(lmMst1)) -
        108*(-642 + 203*lmMst1)*pow2(lmMst2) + 21060*pow3(lmMst2))*pow4(Mst1) +
        27*pow2(Mst1)*(28683 + 5280*B4 - 48*DN - 5952*lmMst1 - 224*OepS2 + 324*
        (-1 + 14*lmMst1 - 14*lmMst2)*S2 - 768*pow2(lmMst1) - 192*lmMst2*(-214 +
        73*lmMst1 + 4*pow2(lmMst1)) - 96*(-268 + 57*lmMst1)*pow2(lmMst2) +
        6240*pow3(lmMst2))*pow4(Mst2) + 2*(1702429 + 257904*B4 - 648*DN -
        748656*lmMst1 - 30112*OepS2 + 108*(9185 + 5646*lmMst1 - 5646*lmMst2)*S2
        + 41904*pow2(lmMst1) + 216*lmMst2*(5971 - 6106*lmMst1 + 576*pow2(
        lmMst1)) - 41904*(-34 + 15*lmMst1)*pow2(lmMst2) - 3456*pow3(lmMst1) +
        507600*pow3(lmMst2))*pow6(Mst1) + 41472*pow2(1 + lmMst2)*pow6(Mst2))))
        + 4860*Mst2*s2t*xDR2DRMOD*(64*Mt*pow2(Msq)*pow2(Mst1)*(((1 + lmMst2)*
        Mgl*(7 - 32*lmMst2 + 4*lmMst1*(7 + 3*lmMst2) - 12*pow2(lmMst2)) +
        Dmglst2*(49 + 103*lmMst2 - 36*(1 + lmMst2)*pow2(lmMst2) + 4*lmMst1*(-11
        + 6*lmMst2 + 9*pow2(lmMst2))))*pow6(Mst1) + 2*(pow2(Mst2)*((1 + lmMst2)
        *Mgl*(1 - 10*lmMst2 + 4*lmMst1*(2 + lmMst2) - 4*pow2(lmMst2)) + 2*
        Dmglst2*(8 + 13*lmMst2 - 8*pow2(lmMst2) + lmMst1*(-5 + 5*lmMst2 + 6*
        pow2(lmMst2)) - 6*pow3(lmMst2)))*pow4(Mst1) + pow2(Mst1)*((1 + lmMst2)*
        Mgl*(-1 - 5*lmMst2 + lmMst1*(3 + 2*lmMst2) - 2*pow2(lmMst2)) + Dmglst2*
        (8 + 7*lmMst2 - 11*pow2(lmMst2) + lmMst1*(-3 + 5*lmMst2 + 6*pow2(
        lmMst2)) - 6*pow3(lmMst2)))*pow4(Mst2) - (1 + lmMst2)*(-Dmglst2 + 3*
        Dmglst2*lmMst2 + Mgl + lmMst2*Mgl)*pow6(Mst2))) + Mst2*s2t*((15*
        Dmsqst2*(1 - 2*lmMst1 + 2*lmMst2)*(2*Dmglst2 - 5*Mgl) + pow2(Msq)*(-2*
        Dmglst2*(60 + 206*lmMst2 + 32*lmMst2*pow2(lmMst1) + lmMst1*(8 - 460*
        lmMst2 - 246*pow2(lmMst2)) + 519*pow2(lmMst2) + 214*pow3(lmMst2)) -
        Mgl*(189 + 726*lmMst2 + 32*(1 + lmMst2)*pow2(lmMst1) + 707*pow2(lmMst2)
        - 2*lmMst1*(253 + 332*lmMst2 + 123*pow2(lmMst2)) + 214*pow3(lmMst2))))*
        pow4(Mst1)*pow4(Mst2) + 9*Mgl*shiftst1*(Dmsqst2 + pow2(Msq))*pow2(Mst1)
        *(15*(pow2(Mst1) + pow2(Mst2))*pow4(Mst2) - 30*(lmMst1 - lmMst2)*pow2(
        Mst1)*(pow2(Mst1)*pow2(Mst2) + pow4(Mst1) + pow4(Mst2))) - pow2(Mst2)*(
        30*Dmsqst2*(lmMst1 - lmMst2)*(2*Dmglst2 - 5*Mgl) + pow2(Msq)*(4*
        Dmglst2*(48 + 4*lmMst2*(31 + 36*pow2(lmMst1)) + 278*pow2(lmMst2) -
        lmMst1*(44 + 278*lmMst2 + 379*pow2(lmMst2)) + 235*pow3(lmMst2)) + 2*
        Mgl*(32 + 285*lmMst2 + 144*(1 + lmMst2)*pow2(lmMst1) + 444*pow2(lmMst2)
        - lmMst1*(253 + 588*lmMst2 + 379*pow2(lmMst2)) + 235*pow3(lmMst2))))*
        pow6(Mst1) + (15*Dmsqst2*(2*Dmglst2 - 5*Mgl) + (-(Mgl*(205 + 252*lmMst2
        + 32*lmMst1*(1 + lmMst2) + 91*pow2(lmMst2))) - 2*Dmglst2*(-20 + 2*(99 +
        16*lmMst1)*lmMst2 + 123*pow2(lmMst2)))*pow2(Msq))*pow2(Mst1)*pow6(Mst2)
        + (150*Dmsqst2*(lmMst1 - lmMst2)*Mgl - 2*pow2(Msq)*(2*Dmglst2*(48 + 4*
        lmMst2*(45 + 68*pow2(lmMst1)) + 310*pow2(lmMst2) - lmMst1*(92 + 310*
        lmMst2 + 635*pow2(lmMst2)) + 363*pow3(lmMst2)) + Mgl*(40 + 277*lmMst2 +
        272*(1 + lmMst2)*pow2(lmMst1) + 556*pow2(lmMst2) - lmMst1*(237 + 828*
        lmMst2 + 635*pow2(lmMst2)) + 363*pow3(lmMst2))))*pow8(Mst1) + 16*(1 +
        lmMst2)*(4*Dmglst2*lmMst2 + Mgl + lmMst2*Mgl)*pow2(Msq)*pow8(Mst2))))))
        )/(3.000564e7*pow2(Mgl)*pow4(Msq)*pow4(Mst1)*pow6(Mst2))) - pow2(Mt)*((
        Al4p*xDmglst2*pow2(Dmglst2)*(Al4p*threeLoopFlag*((5*Dmsqst2*(-1081 +
        165*lmMst1 - 165*lmMst2)*pow2(s2t)*(-1 + pow2(Sbeta))*pow2(Sbeta)*pow4(
        Mst1))/(9.*pow2(Msq)*pow2(Mst2)) - (pow2(MuSUSY)*((36*pow2(Mt)*(24*
        pow2(Mst2)*(50134 - 270*B4 + 270*D3 - 135*DN + 2160*lmMst1 + 120*(271 +
        24*lmMst1)*lmMst2 - 120*lmMt - 29646*S2 - 720*(-2 + 3*lmMst1)*pow2(
        lmMst2) + 2160*pow3(lmMst2)) + pow2(Mst1)*(3051661 - 12960*B4 + 12960*
        D3 - 6480*DN + 304320*lmMst1 + 960*(2227 + 285*lmMst1)*lmMst2 + 2880*(4
        + lmMst1 - lmMst2)*lmMt + 5600*OepS2 - 324*(5361 + 350*lmMst1 - 350*
        lmMst2)*S2 - 2880*(23 + 72*lmMst1)*pow2(lmMst2) + 207360*pow3(lmMst2)))
        )/pow4(Mst2) + 43740*pow2(s2t)*(90.77757201646091 - (100*B4)/9. + (112*
        D3)/9. - (62*DN)/9. + (788*lmMst1)/9. - (103*pow2(lmMst1))/9. + (16*
        lmMst2*(-241 - 171*lmMst1 + 18*pow2(lmMst1)))/27. + (125.33333333333333
         - 82*lmMst1)*pow2(lmMst2) + (80*Dmsqst2*(4 + lmMst1 - lmMst2))/(3.*
        pow2(Msq)) + ((360*Dmsqst2 + (-852 + 358*lmMst2 + 32*lmMst1*(-2 + 3*
        lmMst2) + 497*pow2(lmMst2))*pow2(Msq))*pow2(Mst2))/(9.*pow2(Msq)*pow2(
        Mst1)) + (214*pow3(lmMst2))/3. + (pow2(Mst1)*(155.38193689986284 - (
        164*B4)/9. + (176*D3)/9. - (94*DN)/9. - (21449*lmMst1)/675. - (7556*
        pow2(lmMst1))/135. + (lmMst2*(-54451 - 22990*lmMst1 + 65550*pow2(
        lmMst1)))/675. + (2*(6077 - 17130*lmMst1)*pow2(lmMst2))/135. + (10*
        Dmsqst2*(-23 + 42*lmMst1 - 42*lmMst2))/(3.*pow2(Msq)) - (10*pow3(
        lmMst1))/27. + (4240*pow3(lmMst2))/27.))/pow2(Mst2) + ((
        585.1892843532082 - (8*B4)/3. + (32*D3)/9. - (20*DN)/9. - (20109937*
        lmMst1)/297675. - (15886*pow2(lmMst1))/945. + lmMst2*(
        17.112243218274966 - (144628*lmMst1)/945. + (1180*pow2(lmMst1))/9.) + (
        169.85608465608465 - (2632*lmMst1)/9.)*pow2(lmMst2) - (92*pow3(lmMst1))
        /27. + (4448*pow3(lmMst2))/27.)*pow4(Mst1))/pow4(Mst2) - (32*(-2 +
        lmMst2 + 5*pow2(lmMst2))*pow4(Mst2))/(9.*pow4(Mst1)) - (6*(7400*OepS2 +
        27*(1089707 - 5550*lmMst1 + 5550*lmMst2)*S2)*pow2(Mst1)*pow2(Mst2) +
        40*(17308*OepS2 + 27*(93919 - 12981*lmMst1 + 12981*lmMst2)*S2)*pow4(
        Mst1) + 9*(1400*OepS2 + 81*(116129 - 350*lmMst1 + 350*lmMst2)*S2)*pow4(
        Mst2))/(10935.*pow4(Mst2))) + (240*T1ep*(30*Mst2*pow2(Mst1)*(311*Mst2*
        Mt*s2t - 42*pow2(Mt) + 37*pow2(Mst2)*pow2(s2t)) + 2*s2t*(-16421*Mt +
        8654*Mst2*s2t)*pow4(Mst1) + 63*s2t*(Mt + 5*Mst2*s2t)*pow4(Mst2)))/pow5(
        Mst2) - (Mt*s2t*(pow2(Mst2)*((2332800*Dmsqst2*((-8 + 15*lmMst1 - 15*
        lmMst2)*pow2(Mst1) + (2 + 5*lmMst1 - 5*lmMst2)*pow2(Mst2)))/pow2(Msq) +
        6*pow2(Mst1)*(4580489 - 3285360*B4 + 68040*DN - 3918600*lmMst1 +
        248800*OepS2 - 108*(-20803 + 46650*lmMst1 - 46650*lmMst2)*S2 + 103680*
        pow2(lmMst1) - 360*lmMst2*(11213 + 1368*lmMst1 + 144*pow2(lmMst1)) +
        6480*(-214 + 523*lmMst1)*pow2(lmMst2) - 3337200*pow3(lmMst2))) + 10*(
        2773621 - 1660176*B4 + 25272*DN + 2004408*lmMst1 - 525472*OepS2 + 108*(
        123113 + 98526*lmMst1 - 98526*lmMst2)*S2 + 3888*pow2(lmMst1) - 144*
        lmMst2*(36802 - 11421*lmMst1 + 1728*pow2(lmMst1)) + 167184*(-14 + 15*
        lmMst1)*pow2(lmMst2) - 31104*pow3(lmMst1) - 2227824*pow3(lmMst2))*pow4(
        Mst1) + 9*(2601353 - 1058400*B4 + 23760*DN - 1779840*lmMst1 + 1120*
        OepS2 - 324*(-1387 + 70*lmMst1 - 70*lmMst2)*S2 - 69120*pow2(lmMst1) +
        1440*lmMst2*(-599 - 108*lmMst1 + 24*pow2(lmMst1)) + 4320*(-222 + 217*
        lmMst1)*pow2(lmMst2) - 972000*pow3(lmMst2))*pow4(Mst2) + (1244160*(2 +
        lmMst2 - 3*pow2(lmMst2))*pow6(Mst2))/pow2(Mst1)))/pow5(Mst2)))/43740.)
        + (s2t*twoLoopFlag*pow2(MuSUSY)*(2*(4*Mt*(-31 + 3*lmMst1*(-2 + lmMst2)
        + 4*lmMst2 - 3*pow2(lmMst2)) + 3*Mst2*s2t*(4 + lmMst2 + lmMst1*(-1 + 2*
        lmMst2) - 2*pow2(lmMst2)))*pow2(Mst2)*pow4(Mst1) + 2*(Mst2*s2t*(10 +
        lmMst1*(-3 + 6*lmMst2) - 6*pow2(lmMst2)) + 4*Mt*(-14 + lmMst1*(-2 +
        lmMst2) - pow2(lmMst2)))*pow2(Mst1)*pow4(Mst2) + (Mt*(-398 + 52*lmMst2
        + lmMst1*(-68 + 40*lmMst2) - 40*pow2(lmMst2)) + 3*Mst2*s2t*(9 + 4*
        lmMst1*(-1 + lmMst2) + 4*lmMst2 - 4*pow2(lmMst2)))*pow6(Mst1) + 2*(2 -
        3*lmMst2)*s2t*pow7(Mst2)))/(3.*pow2(Mst1)*pow5(Mst2))))/pow2(Mgl) + (
        s2t*((2*xMst*(27*(lmMst1 - lmMst2)*Mst2*oneLoopFlag*s2t*pow2(MuSUSY) -
        (2*Al4p*twoLoopFlag*(4*(18*(lmMst1 - lmMst2)*(-2 + 3*lmMst2)*Mst2*s2t*
        xDmglst2*xDR2DRMOD*pow2(Dmglst2) + Dmglst2*Mgl*(Mt*(785 + 6*lmMst1*(85
        - 42*lmMst2) - 438*lmMst2 + 252*pow2(lmMst2)) - Mst2*s2t*(49 - 84*
        lmMst2 + lmMst1*(84 - 36*lmMst2*(-1 + xDR2DRMOD)) + 36*(-1 + xDR2DRMOD)
        *pow2(lmMst2))) + (Mt*(193 + 474*lmMst2 - 6*lmMst1*(67 + 42*lmMst2) +
        252*pow2(lmMst2)) - Mst2*s2t*(1 + 3*lmMst2*(-37 + 6*xDR2DRMOD) - 3*
        lmMst1*(-37 + 6*lmMst2*(-12 + xDR2DRMOD) + 6*xDR2DRMOD) - 81*pow2(
        lmMst1) + 9*(-15 + 2*xDR2DRMOD)*pow2(lmMst2)))*pow2(Mgl))*pow2(MuSUSY)
        - (xDmglst2*pow2(Dmglst2)*(4*Mt*Tbeta*(-2501 + 222*lmMst2 + 42*lmMst1*(
        -7 + 6*lmMst2) - 252*pow2(lmMst2))*pow2(MuSUSY) + 6*Mst2*MuSUSY*(
        MuSUSY*s2t*Tbeta*(85 - 60*lmMst1 + 60*lmMst2 + 36*lmMst1*lmMst2 - 36*
        pow2(lmMst2)) + 10*(-43 + 60*lmMst1 - 60*lmMst2)*Mt*pow2(Sbeta)) + 15*(
        -43 + 60*lmMst1 - 60*lmMst2)*s2t*Tbeta*(-1 + pow2(Sbeta))*pow2(Sbeta)*
        pow3(Mst2)))/Tbeta))/pow2(Mgl))*pow6(Mst1))/pow7(Mst2) - pow2(MuSUSY)*(
        27*oneLoopFlag*s2t*(2 - lmMst1 + lmMst2 - (2*(lmMst1 - lmMst2)*pow2(
        Mst1)*(pow2(Mst1) + pow2(Mst2)))/pow4(Mst2)) + Al4p*((36*twoLoopFlag*(
        Mgl*(4*(4*Mt*(5 + 6*lmMst2 - lmMst1*(4 + 3*lmMst2) + 3*pow2(lmMst2)) +
        Mst2*s2t*(-1 + 13*lmMst2 - lmMst1*(13 + 8*lmMst2) + pow2(lmMst1) + 7*
        pow2(lmMst2)))*pow2(Mst2)*pow4(Mst1) + 2*(-(Mst2*s2t*(-14 - 20*lmMst2 +
        2*lmMst1*(5 + 3*lmMst2) + pow2(lmMst1) - 7*pow2(lmMst2))) + 8*Mt*(4 +
        3*lmMst2 - lmMst1*(1 + lmMst2) + pow2(lmMst2)))*pow2(Mst1)*pow4(Mst2) +
        (Mst2*s2t*(-1 + 50*lmMst2 - 2*lmMst1*(25 + 32*lmMst2) + 20*pow2(lmMst1)
        + 44*pow2(lmMst2)) + Mt*(84 + 152*lmMst2 - 40*lmMst1*(3 + 2*lmMst2) +
        80*pow2(lmMst2)))*pow6(Mst1) + 4*(1 + lmMst2)*s2t*pow7(Mst2)) + 4*
        Dmglst2*((Mst2*s2t*(-5 + 8*lmMst2 - 4*lmMst1*(2 + lmMst2) + 4*pow2(
        lmMst2)) + Mt*(65 + lmMst1*(34 - 20*lmMst2) - 26*lmMst2 + 20*pow2(
        lmMst2)))*pow6(Mst1) + 2*((Mst2*s2t*(-2 + 3*lmMst2 - lmMst1*(3 + 2*
        lmMst2) + 2*pow2(lmMst2)) + Mt*(22 + lmMst1*(8 - 6*lmMst2) - 4*lmMst2 +
        6*pow2(lmMst2)))*pow2(Mst2)*pow4(Mst1) + (-(Mst2*s2t*(1 + lmMst1 - 2*
        lmMst2 + 2*lmMst1*lmMst2 - 2*pow2(lmMst2))) + 2*Mt*(6 + lmMst1 + lmMst2
        - lmMst1*lmMst2 + pow2(lmMst2)))*pow2(Mst1)*pow4(Mst2) + lmMst2*s2t*
        pow7(Mst2)))))/(Mgl*pow2(Mst1)*pow5(Mst2)) - 4*xDR2DRMOD*(9*s2t*((4*(2*
        Dmglst2*lmMst2 + Mgl + lmMst2*Mgl)*twoLoopFlag*((pow2(Mst1) + pow2(
        Mst2))*pow4(Mst2) - 2*(lmMst1 - lmMst2)*pow2(Mst1)*(pow2(Mst1)*pow2(
        Mst2) + pow4(Mst1) + pow4(Mst2))))/(Mgl*pow2(Mst1)*pow4(Mst2)) + (Al4p*
        threeLoopFlag*xDmsqst2*pow2(Dmsqst2)*(55 - 110*lmMst1 + 110*lmMst2 + (
        50*Dmglst2*(1 - 2*lmMst1 + 2*lmMst2))/Mgl + (5*(10*Dmglst2 + 11*Mgl)*
        pow2(Mst2))/(Mgl*pow2(Mst1)) - (6*xDmglst2*pow2(Dmglst2)*(5*lmMst1 - 5*
        lmMst2 - (50*(lmMst1 - lmMst2)*pow2(Mst1))/pow2(Mst2) - (5*(1 + pow2(
        Mst2)/pow2(Mst1)))/2.))/pow2(Mgl) + (10*(lmMst1 - lmMst2)*pow2(Mst1)*(-
        10*pow2(Mst1) + ((4*Dmglst2 - 11*Mgl)*pow2(Mst2))/Mgl))/pow4(Mst2) -
        90*shiftst1*(1 + pow2(Mst2)/pow2(Mst1) - (2*(lmMst1 - lmMst2)*(pow2(
        Mst1)*pow2(Mst2) + pow4(Mst1) + pow4(Mst2)))/pow4(Mst2))))/pow4(Msq)) +
        (2*xDmglst2*pow2(Dmglst2)*(-18*(2 - 3*lmMst2)*Mst2*s2t*twoLoopFlag*
        pow2(Mst1)*((pow2(Mst1) + pow2(Mst2))*pow4(Mst2) - 2*(lmMst1 - lmMst2)*
        pow2(Mst1)*(pow2(Mst1)*pow2(Mst2) + pow4(Mst1) + pow4(Mst2))) + Al4p*
        threeLoopFlag*(-((Mst2*s2t*(1180 + 1270*lmMst2 + 96*(2 - 3*lmMst2)*
        pow2(lmMst1) - 3255*pow2(lmMst2) + 6*lmMst1*(-468 + 454*lmMst2 + 369*
        pow2(lmMst2)) - 1926*pow3(lmMst2)) + 64*Mt*(53 + 191*lmMst2 - 54*pow2(
        lmMst2) + 6*lmMst1*(-10 - 3*lmMst2 + 12*pow2(lmMst2)) - 72*pow3(lmMst2)
        ))*pow4(Mst1)*pow4(Mst2)) - 2*pow2(Mst2)*(3*Mst2*s2t*(240 + 468*lmMst2
        + 144*(2 - 3*lmMst2)*pow2(lmMst1) - 310*pow2(lmMst2) + lmMst1*(-580 +
        22*lmMst2 + 1137*pow2(lmMst2)) - 705*pow3(lmMst2)) + 32*Mt*(11 + 371*
        lmMst2 - 42*pow2(lmMst2) + 6*lmMst1*(-21 - 5*lmMst2 + 24*pow2(lmMst2))
        - 144*pow3(lmMst2)))*pow6(Mst1) - 3*(Mst2*s2t*(420 + lmMst1*(64 - 96*
        lmMst2) - 358*lmMst2 - 497*pow2(lmMst2)) + 256*Mt*(2 + lmMst2 - 3*pow2(
        lmMst2)))*pow2(Mst1)*pow6(Mst2) - 12*(16*Mt*(49 + 103*lmMst2 - 36*(1 +
        lmMst2)*pow2(lmMst2) + 4*lmMst1*(-11 + 6*lmMst2 + 9*pow2(lmMst2))) -
        Mst2*s2t*(48 + 4*lmMst2*(45 + 68*pow2(lmMst1)) + 310*pow2(lmMst2) -
        lmMst1*(92 + 310*lmMst2 + 635*pow2(lmMst2)) + 363*pow3(lmMst2)))*pow8(
        Mst1) + s2t*((675*Dmsqst2*pow2(Mst1)*(2*(lmMst1 - lmMst2)*pow2(Mst1) -
        pow2(Mst2))*(pow2(Mst1) + pow2(Mst2))*pow3(Mst2))/pow2(Msq) - 96*(-2 +
        lmMst2 + 5*pow2(lmMst2))*pow9(Mst2)))))/(pow2(Mgl)*pow4(Mst1)*pow5(
        Mst2)))))))/216.);
}

/**
 * @return The diagonal (2, 2) matrix element of the Higgs mass matrix as a double for the hierarchy 'H6bq2g2'
 */
double H6bq2g2::getS2() const {
   return -(oneLoopFlag*((4*Mt*MuSUSY*s2t*(4*(lmMst1 - lmMst2)*pow2(Mt) + ((-2 -
        lmMst1 + lmMst2)*pow2(Mst1) + (2 - lmMst1 + lmMst2)*pow2(Mst2))*pow2(
        s2t)))/Tbeta + 4*pow2(Mt)*pow2(s2t)*(2*(lmMst1 - lmMst2)*(pow2(Mst1) -
        pow2(Mst2)) + pow2(MuSUSY)*(2 - lmMst1 + lmMst2 - (2*(lmMst1 - lmMst2)*
        pow2(Mst1)*(pow2(Mst1) + pow2(Mst2)))/pow4(Mst2))) - (4*pow2(Mt)*pow2(
        MuSUSY)*pow2(s2t)*(2 - lmMst1 + lmMst2 - (2*(lmMst1 - lmMst2)*pow2(
        Mst1)*(pow2(Mst1) + pow2(Mst2)))/pow4(Mst2)))/pow2(Sbeta) + 16*(lmMst1
        + lmMst2 - 2*lmMt)*pow4(Mt) - (-4*pow2(Mst1)*pow2(Mst2) + (2 + lmMst1 -
        lmMst2)*pow4(Mst1) + (2 - lmMst1 + lmMst2)*pow4(Mst2))*pow4(s2t)))/32.
         - (Al4p*twoLoopFlag*((-72*Mt*pow3(s2t)*(-(Mgl*(4*(3 + lmMst1*(2 +
        lmMst2) - pow2(lmMst2))*pow2(Mst1)*pow2(Mst2) + (3 + 2*lmMst1 - 2*
        lmMst2)*pow4(Mst1) + 4*(-4 + lmMst1 - 3*lmMst2 + lmMst1*lmMst2 - pow2(
        lmMst2))*pow4(Mst2))) + Dmglst2*(-4*(1 + lmMst1*(-2 + lmMst2) + 4*
        lmMst2 - pow2(lmMst2))*pow2(Mst1)*pow2(Mst2) + (1 + 6*lmMst1 - 6*
        lmMst2)*pow4(Mst1) + 4*(6 + lmMst1 + lmMst2 - lmMst1*lmMst2 + pow2(
        lmMst2))*pow4(Mst2))))/Mst2 + (32*pow4(Mt)*(2*Dmglst2*((53 - 18*lmMst1
        + 24*lmMst2 - 24*lmMt)*pow2(Mst1)*pow4(Mst2) + 18*((5 + lmMst2*(11 - 2*
        lmMt) - 2*lmMst1*(4 + lmMst2 - lmMt) - 3*lmMt + 2*pow2(lmMst2))*pow2(
        Mst2)*pow4(Mst1) + (3 - 6*lmMst2*(-5 + lmMt) - 5*lmMt + lmMst1*(-25 -
        6*lmMst2 + 6*lmMt) + 6*pow2(lmMst2))*pow6(Mst1)) - 18*lmMst2*pow6(Mst2)
        ) + 9*Mgl*(2*(2 + 2*lmMst1*(3 + lmMst2 - lmMt) + lmMt + lmMst2*(-7 + 2*
        lmMt) - 2*pow2(lmMst2))*pow2(Mst2)*pow4(Mst1) + (3 + 2*lmMst1*(3 +
        lmMst2) + lmMt + lmMst2*(-5 + 4*lmMt) + pow2(lmMst1) - pow2(lmMst2) -
        6*pow2(lmMt))*pow2(Mst1)*pow4(Mst2) + (9 + lmMst1*(22 + 6*lmMst2 - 6*
        lmMt) + 6*lmMst2*(-4 + lmMt) + 2*lmMt - 6*pow2(lmMst2))*pow6(Mst1) - 2*
        (1 + lmMst2)*pow6(Mst2))))/(pow2(Mst1)*pow4(Mst2)) - (144*s2t*pow3(Mt)*
        (Dmglst2*(4*pow2(Mst1)*pow2(Mst2)*((6 + 13*lmMst2 - 4*lmMt - 2*lmMst2*
        lmMt + lmMst1*(-9 - 2*lmMst2 + 2*lmMt) + 2*pow2(lmMst2))*pow2(Mst2) + (
        -11 + 2*lmMst2 + lmMst1*(-4 + 3*lmMst2) - 3*pow2(lmMst2))*pow2(MuSUSY))
        - (4*(1 - 29*lmMst2 + lmMst1*(25 + 6*lmMst2 - 6*lmMt) + 4*lmMt + 6*
        lmMst2*lmMt - 6*pow2(lmMst2))*pow2(Mst2) + (65 + lmMst1*(34 - 20*
        lmMst2) - 26*lmMst2 + 20*pow2(lmMst2))*pow2(MuSUSY))*pow4(Mst1) - 4*(6
        + lmMst1 + lmMst2 - lmMst1*lmMst2 + pow2(lmMst2))*pow2(MuSUSY)*pow4(
        Mst2) + 8*(5 + lmMst1*(-1 + lmMst2) + lmMst2 - pow2(lmMst2))*pow6(Mst2)
        ) + Mgl*((4*(5 - 7*lmMst2 + lmMst1*(7 + 2*lmMst2 - 2*lmMt) + 2*lmMst2*
        lmMt - 2*pow2(lmMst2))*pow2(Mst2) + (-21 - 38*lmMst2 + 10*lmMst1*(3 +
        2*lmMst2) - 20*pow2(lmMst2))*pow2(MuSUSY))*pow4(Mst1) + 4*((-4 + lmMst1
        - 3*lmMst2 + lmMst1*lmMst2 - pow2(lmMst2))*pow2(MuSUSY)*pow4(Mst2) +
        pow2(Mst1)*((-5 - 6*lmMst2 + lmMst1*(4 + 3*lmMst2) - 3*pow2(lmMst2))*
        pow2(Mst2)*pow2(MuSUSY) + (4 - 5*lmMst2 + lmMst1*(5 + 2*lmMst2 - 2*
        lmMt) + 2*lmMst2*lmMt - 2*pow2(lmMst2))*pow4(Mst2)) + (3 - 4*lmMst2 +
        2*(lmMst1 + lmMst1*lmMst2 + lmMt) - 2*pow2(lmMst2))*pow6(Mst2)))))/
        pow5(Mst2) + (36*Mt*MuSUSY*(16*pow3(Mt)*(Mgl*(13 + 2*lmMst1*(6 + 3*
        lmMst2 - 2*lmMt) + 2*lmMt + 2*lmMst2*(-7 + 2*lmMt) - 6*pow2(lmMst2))*
        pow4(Mst1) + Dmglst2*(11 - 6*lmMst1*lmMst2 - 8*lmMst2*(-5 + lmMt) + 8*
        lmMst1*(-4 + lmMt) - 8*lmMt + 6*pow2(lmMst2))*pow4(Mst1) + 2*Dmglst2*((
        7 + 7*lmMst2 + lmMst1*(-5 + lmMt) - 2*lmMt - lmMst2*lmMt)*pow2(Mst1)*
        pow2(Mst2) + (5 + lmMst1*(-1 + lmMst2) + lmMst2 - pow2(lmMst2))*pow4(
        Mst2)) + 2*Mgl*((4 + lmMst1*(3 + 2*lmMst2 - lmMt) + lmMst2*(-4 + lmMt)
        + lmMt - 2*pow2(lmMst2))*pow2(Mst1)*pow2(Mst2) + (2 + lmMst1 + (-2 +
        lmMst1)*lmMst2 + lmMt - pow2(lmMst2))*pow4(Mst2))) + 6*Mt*pow2(Mst2)*
        pow2(s2t)*(Dmglst2*(4*(5 + lmMst1*(3 - 2*lmMst2) - 3*lmMst2 + 2*pow2(
        lmMst2))*pow2(Mst1)*pow2(Mst2) + (21 + lmMst1*(18 - 8*lmMst2) - 18*
        lmMst2 + 8*pow2(lmMst2))*pow4(Mst1) + 4*(6 + lmMst1 + lmMst2 - lmMst1*
        lmMst2 + pow2(lmMst2))*pow4(Mst2)) + Mgl*(4*(1 + 3*lmMst2 - lmMst1*(3 +
        2*lmMst2) + 2*pow2(lmMst2))*pow2(Mst1)*pow2(Mst2) + (1 + 14*lmMst2 - 2*
        lmMst1*(7 + 4*lmMst2) + 8*pow2(lmMst2))*pow4(Mst1) + 4*(4 + 3*lmMst2 -
        lmMst1*(1 + lmMst2) + pow2(lmMst2))*pow4(Mst2))) + (pow3(Mst2)*pow3(
        s2t)*(-4*Dmglst2*(2*(1 + 2*lmMst1 - lmMst2)*pow2(Mst2)*pow4(Mst1) + 2*(
        1 + lmMst1 - lmMst2 + 2*lmMst1*lmMst2 - 2*pow2(lmMst2))*pow2(Mst1)*
        pow4(Mst2) + (1 + 2*lmMst1 - 2*lmMst2)*pow6(Mst1) - 2*lmMst2*pow6(Mst2)
        ) + Mgl*(2*(-16 + 6*lmMst2 - 2*lmMst1*(8 + 5*lmMst2) + 3*pow2(lmMst1) +
        7*pow2(lmMst2))*pow2(Mst2)*pow4(Mst1) - 2*(-12 - 18*lmMst2 + 2*lmMst1*(
        5 + 3*lmMst2) + pow2(lmMst1) - 7*pow2(lmMst2))*pow2(Mst1)*pow4(Mst2) +
        (3 + lmMst1*(2 - 32*lmMst2) - 2*lmMst2 + 16*pow2(lmMst1) + 16*pow2(
        lmMst2))*pow6(Mst1) + 4*(1 + lmMst2)*pow6(Mst2))))/pow2(Mst1) + (8*
        Mst2*s2t*pow2(Mt)*(Dmglst2*(4*(1 - lmMst1 + lmMst2)*pow2(Mst2)*pow4(
        Mst1) + 2*(1 + 2*lmMst2)*pow2(Mst1)*pow4(Mst2) + (4 - 8*lmMst1 + 8*
        lmMst2)*pow6(Mst1) - 4*lmMst2*pow6(Mst2)) - Mgl*(-2*lmMst1*((1 + 2*
        lmMst2)*(pow2(Mst1) + pow2(Mst2))*pow4(Mst1) + (3 + lmMst2)*pow2(Mst1)*
        pow4(Mst2)) + pow2(lmMst1)*(2*pow2(Mst2)*pow4(Mst1) - pow2(Mst1)*pow4(
        Mst2) + 2*pow6(Mst1)) + pow2(lmMst2)*(2*pow2(Mst2)*pow4(Mst1) + 3*pow2(
        Mst1)*pow4(Mst2) + 2*pow6(Mst1)) + 2*pow6(Mst2) + 2*lmMst2*(pow2(Mst2)*
        pow4(Mst1) + 2*pow2(Mst1)*pow4(Mst2) + pow6(Mst1) + pow6(Mst2)))))/
        pow2(Mst1)))/(Tbeta*pow5(Mst2)) - 36*pow2(Mt)*pow2(s2t)*(4*(Dmglst2*(2
        - 4*lmMst1) + Mgl*(2*lmMst1*(-2 + lmMst2) + lmMst2*(2 + lmMst2) - 3*
        pow2(lmMst1)))*pow2(Mst1) + 4*(Dmglst2*(2 + 8*lmMst2) + Mgl*(2 - 2*
        lmMst2 + 2*lmMst1*(3 + lmMst2) + pow2(lmMst1) - 3*pow2(lmMst2)))*pow2(
        Mst2) - (16*Dmglst2*(lmMst1 - lmMst2)*pow4(Mst1))/pow2(Mst2) + (-8*(2*
        Dmglst2*lmMst2 + Mgl + lmMst2*Mgl)*pow4(Mst2) + (pow2(MuSUSY)*(4*
        Dmglst2*(2*(2 - 3*lmMst2 + lmMst1*(3 + 2*lmMst2) - 2*pow2(lmMst2))*
        pow2(Mst2)*pow4(Mst1) + 2*(1 + lmMst1 - 2*lmMst2 + 2*lmMst1*lmMst2 - 2*
        pow2(lmMst2))*pow2(Mst1)*pow4(Mst2) + (5 - 8*lmMst2 + 4*lmMst1*(2 +
        lmMst2) - 4*pow2(lmMst2))*pow6(Mst1) - 2*lmMst2*pow6(Mst2)) + Mgl*(-4*(
        -1 - 13*lmMst1 + 13*lmMst2 - 8*lmMst1*lmMst2 + pow2(lmMst1) + 7*pow2(
        lmMst2))*pow2(Mst2)*pow4(Mst1) + 2*(-14 + 10*lmMst1 - 20*lmMst2 + 6*
        lmMst1*lmMst2 + pow2(lmMst1) - 7*pow2(lmMst2))*pow2(Mst1)*pow4(Mst2) +
        (1 + 50*lmMst1 - 50*lmMst2 + 64*lmMst1*lmMst2 - 20*pow2(lmMst1) - 44*
        pow2(lmMst2))*pow6(Mst1) - 4*(1 + lmMst2)*pow6(Mst2))))/pow4(Mst2))/
        pow2(Mst1)) + (-9*pow4(s2t)*(4*Dmglst2*(-2*(lmMst1 - 2*lmMst1*lmMst2 +
        2*pow2(lmMst2))*pow2(Mst2)*pow4(Mst1) - 2*(1 + lmMst1 + 2*lmMst1*lmMst2
        - 2*pow2(lmMst2))*pow2(Mst1)*pow4(Mst2) + (1 + 2*lmMst1)*pow6(Mst1) +
        2*lmMst2*pow6(Mst2)) + Mgl*(-4*(14 + 6*lmMst2 + lmMst1*(3 + 2*lmMst2) -
        2*pow2(lmMst1))*pow2(Mst2)*pow4(Mst1) - 2*(-10 - 16*lmMst2 + 2*lmMst1*(
        5 + 3*lmMst2) + pow2(lmMst1) - 7*pow2(lmMst2))*pow2(Mst1)*pow4(Mst2) +
        (35 + lmMst1*(34 - 12*lmMst2) - 14*lmMst2 + 10*pow2(lmMst1) + 2*pow2(
        lmMst2))*pow6(Mst1) + 4*(1 + lmMst2)*pow6(Mst2))) - (36*s2t*pow2(Mt)*
        pow2(MuSUSY)*(Mgl*(4*(4*Mt*(5 + 6*lmMst2 - lmMst1*(4 + 3*lmMst2) + 3*
        pow2(lmMst2)) + Mst2*s2t*(-1 + 13*lmMst2 - lmMst1*(13 + 8*lmMst2) +
        pow2(lmMst1) + 7*pow2(lmMst2)))*pow2(Mst2)*pow4(Mst1) + 2*(-(Mst2*s2t*(
        -14 - 20*lmMst2 + 2*lmMst1*(5 + 3*lmMst2) + pow2(lmMst1) - 7*pow2(
        lmMst2))) + 8*Mt*(4 + 3*lmMst2 - lmMst1*(1 + lmMst2) + pow2(lmMst2)))*
        pow2(Mst1)*pow4(Mst2) + (Mst2*s2t*(-1 + 50*lmMst2 - 2*lmMst1*(25 + 32*
        lmMst2) + 20*pow2(lmMst1) + 44*pow2(lmMst2)) + Mt*(84 + 152*lmMst2 -
        40*lmMst1*(3 + 2*lmMst2) + 80*pow2(lmMst2)))*pow6(Mst1) + 4*(1 +
        lmMst2)*s2t*pow7(Mst2)) + 4*Dmglst2*((Mst2*s2t*(-5 + 8*lmMst2 - 4*
        lmMst1*(2 + lmMst2) + 4*pow2(lmMst2)) + Mt*(65 + lmMst1*(34 - 20*
        lmMst2) - 26*lmMst2 + 20*pow2(lmMst2)))*pow6(Mst1) + 2*((Mst2*s2t*(-2 +
        3*lmMst2 - lmMst1*(3 + 2*lmMst2) + 2*pow2(lmMst2)) + Mt*(22 + lmMst1*(8
        - 6*lmMst2) - 4*lmMst2 + 6*pow2(lmMst2)))*pow2(Mst2)*pow4(Mst1) + (-(
        Mst2*s2t*(1 + lmMst1 - 2*lmMst2 + 2*lmMst1*lmMst2 - 2*pow2(lmMst2))) +
        2*Mt*(6 + lmMst1 + lmMst2 - lmMst1*lmMst2 + pow2(lmMst2)))*pow2(Mst1)*
        pow4(Mst2) + lmMst2*s2t*pow7(Mst2)))))/(pow2(Sbeta)*pow5(Mst2)))/pow2(
        Mst1)))/(216.*Mgl) + (Al4p*xDmglst2*pow2(Dmglst2)*(Al4p*threeLoopFlag*(
        (2187.9919341563786 + (1360*lmMst1)/9. - (4*(50551 + 7200*lmMst1)*
        lmMst2)/2025. + (16*(4997 - 240*lmMst1 + 686*lmMst2)*lmMt)/135. - (28*
        pow2(lmMst1))/9. + (16*(-311 + 720*lmMst1)*pow2(lmMst2))/135. + (464*
        pow2(lmMt))/9. + (40*Dmsqst2*(2954 - 864*lmMst1 + 1713*lmMst2 - 849*
        lmMt))/(81.*pow2(Msq)) + ((4*(360*Dmsqst2 + (-788 - 26*lmMst2 + 32*
        lmMst1*(-2 + 3*lmMst2) + 177*pow2(lmMst2))*pow2(Msq))*pow2(Mst2))/(9.*
        pow2(Mst1)) + (pow2(Mst1)*(8505000*Dmsqst2*(401 + lmMst2*(274 - 4*lmMt)
        - 2*lmMst1*(71 + 2*lmMst2 - 2*lmMt) - 132*lmMt + 4*pow2(lmMst2)) +
        pow2(Msq)*(-55610713 + 1448685000*lmMt - 37800*(1783 + 285*lmMst2 -
        225*lmMt)*pow2(lmMst1) - 18900*(-43477 + 930*lmMt)*pow2(lmMst2) -
        95256000*pow2(lmMt) + 3780*lmMst1*(-75993 + 19265*lmMt + 15*lmMst2*(-
        9947 + 160*lmMt) + 5400*pow2(lmMst2) + 7200*pow2(lmMt)) - 1260*lmMst2*(
        -297379 + 133245*lmMt + 21600*pow2(lmMt)) + 756000*pow3(lmMst1) -
        10395000*pow3(lmMst2))))/(637875.*pow2(Mst2)))/pow2(Msq) - (256*pow3(
        lmMst2))/3. - ((2929.938520304849 + (55510684*lmMst1)/59535. - (126272*
        pow2(lmMst1))/189. - (4*lmMst2*(42300121 + 12578580*lmMst1 + 2487240*
        pow2(lmMst1)))/59535. + (32*(10166 - 693*lmMst1)*pow2(lmMst2))/189. + (
        8*lmMt*(5695 + 1974*lmMst2 - 12*lmMst1*(163 + 47*lmMst2) + 468*pow2(
        lmMst1) + 96*pow2(lmMst2)))/27. + (128*(-5 + 6*lmMst1 - 6*lmMst2)*pow2(
        lmMt))/3. + (256*pow3(lmMst1))/27. + (7424*pow3(lmMst2))/27.)*pow4(
        Mst1))/pow4(Mst2) - (128*(-2 + lmMst2 + 5*pow2(lmMst2))*pow4(Mst2))/(9.
        *pow4(Mst1)) + (-63*pow2(MuSUSY)*(24*pow2(Mst2)*(50134 - 270*B4 + 270*
        D3 - 135*DN + 2160*lmMst1 + 120*(271 + 24*lmMst1)*lmMst2 - 120*lmMt -
        720*(-2 + 3*lmMst1)*pow2(lmMst2) + 2160*pow3(lmMst2)) + pow2(Mst1)*(
        3051661 - 12960*B4 + 12960*D3 - 6480*DN + 304320*lmMst1 + 960*(2227 +
        285*lmMst1)*lmMst2 + 2880*(4 + lmMst1 - lmMst2)*lmMt - 2880*(23 + 72*
        lmMst1)*pow2(lmMst2) + 207360*pow3(lmMst2))) + 1120*OepS2*(-3*pow2(
        Mst1)*(4703*pow2(Mst2) + 105*pow2(MuSUSY)) + 11164*pow4(Mst1) - 819*
        pow4(Mst2)) - 54*S2*(3*pow2(Mst1)*((128797 - 1975260*lmMst1 + 1975260*
        lmMst2)*pow2(Mst2) - 126*(5361 + 350*lmMst1 - 350*lmMst2)*pow2(MuSUSY))
        + (4002484 + 4688880*lmMst1 - 4688880*lmMst2)*pow4(Mst1) + 189*(-4392*
        pow2(Mst2)*pow2(MuSUSY) + (3869 - 1820*lmMst1 + 1820*lmMst2)*pow4(Mst2)
        )))/(76545.*pow4(Mst2)))*pow4(Mt) + pow4(s2t)*(-(pow2(Mst1)*pow2(Mst2)*
        (30.209968449931413 - B4 + (4*D3)/3. - (5*DN)/6. + (144449*lmMst1)/
        2700. + (2233*pow2(lmMst1))/270. - (lmMst2*(165199 + 121010*lmMst1 +
        51150*pow2(lmMst1)))/2700. + (26.353703703703705 + (202*lmMst1)/9.)*
        pow2(lmMst2) + (5*Dmsqst2*(75 - 26*lmMst1 + 26*lmMst2))/(6.*pow2(Msq))
        + (5*pow3(lmMst1))/54. - (97*pow3(lmMst2))/27.)) + (79.58863384550371 +
        (1436147*lmMst1)/1.1907e6 + (1363*pow2(lmMst1))/315. + lmMst2*(
        6.96052994037121 - (101*lmMst1)/315. + (11*pow2(lmMst1))/3.) - (
        0.7285714285714285 + (11*lmMst1)/3.)*pow2(lmMst2) + (5*Dmsqst2*(76 -
        249*lmMst1 + 249*lmMst2))/(54.*pow2(Msq)) - (11*pow3(lmMst1))/9. + (11*
        pow3(lmMst2))/9.)*pow4(Mst1) + (71.80550411522634 - (25*B4)/9. + (28*
        D3)/9. - (31*DN)/18. + (229*lmMst1)/9. - (103*pow2(lmMst1))/36. + (
        lmMst2*(-1525 - 828*lmMst1 + 72*pow2(lmMst1)))/27. - ((13 + 369*lmMst1)
        *pow2(lmMst2))/18. + (20*Dmsqst2*(1 + lmMst1 - lmMst2))/(3.*pow2(Msq))
        + (107*pow3(lmMst2))/6.)*pow4(Mst2) - (2*OepS2*(480*pow2(Mst1)*pow2(
        Mst2) + 8398*pow4(Mst1) + 315*pow4(Mst2)))/2187. + (S2*(12*(-22273 +
        1200*lmMst1 - 1200*lmMst2)*pow2(Mst1)*pow2(Mst2) + 10*(51635 + 25194*
        lmMst1 - 25194*lmMst2)*pow4(Mst1) + 27*(-116129 + 350*lmMst1 - 350*
        lmMst2)*pow4(Mst2)))/1620. + ((360*Dmsqst2 + (-980 + 422*lmMst2 + 32*
        lmMst1*(-2 + 3*lmMst2) + 817*pow2(lmMst2))*pow2(Msq))*pow6(Mst2))/(36.*
        pow2(Msq)*pow2(Mst1))) + (Mt*pow3(s2t)*(583200*Dmsqst2*(4*(12 - 5*
        lmMst1 + 5*lmMst2)*pow2(Mst2)*pow4(Mst1) - 4*(2 + 5*lmMst1 - 5*lmMst2)*
        pow2(Mst1)*pow4(Mst2) + (40 - 3*lmMst1 + 3*lmMst2)*pow6(Mst1)) - pow2(
        Msq)*(-60*pow2(Mst2)*(280885 + 11016*B4 + 324*DN - 142092*lmMst1 -
        24544*OepS2 + 2808*(-32 + 177*lmMst1 - 177*lmMst2)*S2 - 31104*pow2(
        lmMst1) + 36*lmMst2*(3449 + 72*lmMst1 + 432*pow2(lmMst1)) - 648*(134 +
        89*lmMst1)*pow2(lmMst2) + 42120*pow3(lmMst2))*pow4(Mst1) + 9*pow2(Mst1)
        *(2048393 - 1058400*B4 + 23760*DN - 1779840*lmMst1 + 1120*OepS2 - 324*(
        -1387 + 70*lmMst1 - 70*lmMst2)*S2 - 69120*pow2(lmMst1) + 1440*lmMst2*(-
        791 - 108*lmMst1 + 24*pow2(lmMst1)) + 4320*(-30 + 217*lmMst1)*pow2(
        lmMst2) - 972000*pow3(lmMst2))*pow4(Mst2) + (31897243 - 2491360*OepS2 +
        90290268*S2 - 360*lmMst2*(18652 + 140139*S2) + 38880*(37 - 40*lmMst2)*
        pow2(lmMst1) + 3188160*pow2(lmMst2) + 360*lmMst1*(17410 - 12852*lmMst2
        + 140139*S2 + 11232*pow2(lmMst2)) - 311040*pow3(lmMst1) - 2177280*pow3(
        lmMst2))*pow6(Mst1) + 1244160*(2 + lmMst2 - 3*pow2(lmMst2))*pow6(Mst2))
        ))/(87480.*Mst2*pow2(Msq)*pow2(Mst1)) + (T1ep*(-2*pow4(Mst1)*(8*Mst2*
        pow2(Mt)*pow2(s2t)*(4327*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 6921*pow2(
        Mst2)*pow2(Sbeta)) - 4*s2t*(16421*pow2(MuSUSY)*(-1 + pow2(Sbeta)) +
        51460*pow2(Mst2)*pow2(Sbeta))*pow3(Mt) + 15571*Mt*pow2(Sbeta)*pow3(s2t)
        *pow4(Mst2) + 89312*Mst2*pow2(Sbeta)*pow4(Mt) - 4199*pow2(Sbeta)*pow4(
        s2t)*pow5(Mst2)) + 63*pow4(Mst2)*(-4*Mst2*pow2(Mt)*pow2(s2t)*(5*pow2(
        MuSUSY)*(-1 + pow2(Sbeta)) + 18*pow2(Mst2)*pow2(Sbeta)) - 4*s2t*(pow2(
        MuSUSY)*(-1 + pow2(Sbeta)) + 2*pow2(Mst2)*pow2(Sbeta))*pow3(Mt) + 2*Mt*
        pow2(Sbeta)*pow3(s2t)*pow4(Mst2) + 208*Mst2*pow2(Sbeta)*pow4(Mt) + 5*
        pow2(Sbeta)*pow4(s2t)*pow5(Mst2)) + 24*Mst2*pow2(Mst1)*(pow2(Mst2)*
        pow2(Mt)*pow2(s2t)*(-185*pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 1791*pow2(
        Mst2)*pow2(Sbeta)) + Mst2*s2t*(-1555*pow2(MuSUSY)*(-1 + pow2(Sbeta)) -
        4307*pow2(Mst2)*pow2(Sbeta))*pow3(Mt) + (210*pow2(MuSUSY)*(-1 + pow2(
        Sbeta)) + 9406*pow2(Mst2)*pow2(Sbeta))*pow4(Mt) + 767*Mt*pow2(Sbeta)*
        pow3(s2t)*pow5(Mst2) + 20*pow2(Sbeta)*pow4(s2t)*pow6(Mst2))))/(729.*
        pow2(Sbeta)*pow5(Mst2)) + pow2(Mt)*pow2(s2t)*(pow2(Mst2)*(
        1491.268253968254 - (16*B4)/3. + (16*D3)/3. - (8*DN)/3. + (608*lmMst1)/
        9. - (32*lmMt)/9. + (2*lmMst2*(13339 + 1206*lmMst1 - 288*pow2(lmMst1)))
        /27. + (58*pow2(lmMst1))/3. + (2*(473 - 96*lmMst1)*pow2(lmMst2))/9. - (
        920*Dmsqst2)/(9.*pow2(Msq)) + (128*pow3(lmMst2))/3.) + (pow2(Mst1)*(
        273721621 + 420*lmMst2*(463453 - 3600*lmMt) + 9072000*lmMt - 6300*(289
        + 30*lmMst2)*pow2(lmMst1) + 840*lmMst1*(211261 + 90510*lmMst2 + 1800*
        lmMt - 64575*pow2(lmMst2)) - 84886200*pow2(lmMst2) + (945000*Dmsqst2*(-
        113 + 6*lmMst1 - 6*lmMst2))/pow2(Msq) + 63000*pow3(lmMst1) + 54369000*
        pow3(lmMst2)))/425250. - ((175.16754355781114 - (17578814*lmMst1)/
        33075. + lmMst2*(257.7056386999244 + (105592*lmMst1)/315. - (208*pow2(
        lmMst1))/9.) - (35576*pow2(lmMst1))/315. + (16*(-4376 + 2555*lmMst1)*
        pow2(lmMst2))/315. + (64*lmMt*(2 + 5*lmMst2 - lmMst1*(5 + 2*lmMst2) +
        pow2(lmMst1) + pow2(lmMst2)))/3. + (160*Dmsqst2*(41 - 12*lmMst1 + 12*
        lmMst2))/(27.*pow2(Msq)) + (16*pow3(lmMst1))/27. - (2896*pow3(lmMst2))/
        27.)*pow4(Mst1))/pow2(Mst2) - (2*(360*Dmsqst2 + (-820 + 166*lmMst2 +
        32*lmMst1*(-2 + 3*lmMst2) + 337*pow2(lmMst2))*pow2(Msq))*pow4(Mst2))/(
        9.*pow2(Msq)*pow2(Mst1)) - pow2(MuSUSY)*(90.77757201646091 - (100*B4)/
        9. + (112*D3)/9. - (62*DN)/9. + (788*lmMst1)/9. - (103*pow2(lmMst1))/9.
         + (16*lmMst2*(-241 - 171*lmMst1 + 18*pow2(lmMst1)))/27. + (
        125.33333333333333 - 82*lmMst1)*pow2(lmMst2) + (80*Dmsqst2*(4 + lmMst1
        - lmMst2))/(3.*pow2(Msq)) + ((360*Dmsqst2 + (-852 + 358*lmMst2 + 32*
        lmMst1*(-2 + 3*lmMst2) + 497*pow2(lmMst2))*pow2(Msq))*pow2(Mst2))/(9.*
        pow2(Msq)*pow2(Mst1)) + (214*pow3(lmMst2))/3. + (pow2(Mst1)*(
        155.38193689986284 - (164*B4)/9. + (176*D3)/9. - (94*DN)/9. - (21449*
        lmMst1)/675. - (7556*pow2(lmMst1))/135. + (lmMst2*(-54451 - 22990*
        lmMst1 + 65550*pow2(lmMst1)))/675. + (2*(6077 - 17130*lmMst1)*pow2(
        lmMst2))/135. + (10*Dmsqst2*(-23 + 42*lmMst1 - 42*lmMst2))/(3.*pow2(
        Msq)) - (10*pow3(lmMst1))/27. + (4240*pow3(lmMst2))/27.))/pow2(Mst2) +
        ((585.1892843532082 - (8*B4)/3. + (32*D3)/9. - (20*DN)/9. - (20109937*
        lmMst1)/297675. - (15886*pow2(lmMst1))/945. + lmMst2*(
        17.112243218274966 - (144628*lmMst1)/945. + (1180*pow2(lmMst1))/9.) + (
        169.85608465608465 - (2632*lmMst1)/9.)*pow2(lmMst2) - (92*pow3(lmMst1))
        /27. + (4448*pow3(lmMst2))/27.)*pow4(Mst1))/pow4(Mst2) - (32*(-2 +
        lmMst2 + 5*pow2(lmMst2))*pow4(Mst2))/(9.*pow4(Mst1))) - S2*((
        807.2952380952381 + 796*lmMst1 - 796*lmMst2)*pow2(Mst1) + (
        596.0571428571428 + 84*lmMst1 - 84*lmMst2)*pow2(Mst2) + (4*(28283 +
        23070*lmMst1 - 23070*lmMst2)*pow4(Mst1))/(45.*pow2(Mst2)) - (pow2(
        MuSUSY)*(6*(1089707 - 5550*lmMst1 + 5550*lmMst2)*pow2(Mst1)*pow2(Mst2)
        + 40*(93919 - 12981*lmMst1 + 12981*lmMst2)*pow4(Mst1) + 27*(116129 -
        350*lmMst1 + 350*lmMst2)*pow4(Mst2)))/(405.*pow4(Mst2))) + (64*(-2 +
        lmMst2 + 5*pow2(lmMst2))*pow6(Mst2))/(9.*pow4(Mst1)) + (8*OepS2*(4*(
        6921*pow2(Mst2) + 4327*pow2(MuSUSY))*pow4(Mst1) + 6*pow2(Mst1)*(185*
        pow2(Mst2)*pow2(MuSUSY) + 1791*pow4(Mst2)) + 63*(5*pow2(MuSUSY)*pow4(
        Mst2) + 18*pow6(Mst2))))/(2187.*pow4(Mst2))) - (MuSUSY*(s2t*pow3(Mt)*(
        3314.980952380952 - (32*B4)/3. + (32*D3)/3. - (16*DN)/3. + (1472*
        lmMst1)/9. - (64*lmMt)/9. + lmMst2*(1916.5925925925926 + 136*lmMst1 - (
        128*pow2(lmMst1))/3.) + (116*pow2(lmMst1))/3. + (32*(37 - 12*lmMst1)*
        pow2(lmMst2))/9. - (3280*Dmsqst2)/(9.*pow2(Msq)) - (4*(360*Dmsqst2 + (-
        756 + 134*lmMst2 + 32*lmMst1*(-2 + 3*lmMst2) + 177*pow2(lmMst2))*pow2(
        Msq))*pow2(Mst2))/(9.*pow2(Msq)*pow2(Mst1)) + (256*pow3(lmMst2))/3. + (
        pow2(Mst1)*(4623.658770135215 - (32*B4)/3. + (32*D3)/3. - (16*DN)/3. +
        (2028488*lmMst1)/2025. + (64*(5 + lmMst1 - lmMst2)*lmMt)/9. + (8*
        lmMst2*(715964 + 124935*lmMst1 - 11025*pow2(lmMst1)))/2025. + (4064*
        pow2(lmMst1))/135. - (8*(4517 + 5025*lmMst1)*pow2(lmMst2))/135. + (40*
        Dmsqst2*(-65 + 2*lmMst1 - 2*lmMst2))/(3.*pow2(Msq)) + (8*pow3(lmMst1))/
        27. + (9208*pow3(lmMst2))/27.))/pow2(Mst2) + ((622.8225754358181 - (64*
        B4)/9. + (64*D3)/9. - (32*DN)/9. + (29853268*lmMst1)/19845. + (38704*
        pow2(lmMst1))/189. + (4*lmMst2*(2917823 - 3813600*lmMst1 + 97020*pow2(
        lmMst1)))/19845. + (16*(8677 - 5439*lmMst1)*pow2(lmMst2))/189. - (128*
        lmMt*(4 + 6*lmMst2 - 2*lmMst1*(3 + lmMst2) + pow2(lmMst1) + pow2(
        lmMst2)))/3. - (16*pow3(lmMst1))/9. + (1328*pow3(lmMst2))/3.)*pow4(
        Mst1))/pow4(Mst2) + (128*(-2 + lmMst2 + 5*pow2(lmMst2))*pow4(Mst2))/(9.
        *pow4(Mst1)) + ((32*OepS2*(660*pow2(Mst1)*pow2(Mst2) + 1838*pow4(Mst1)
        + 63*pow4(Mst2)))/243. - (4*S2*(12*(18419 + 11550*lmMst1 - 11550*
        lmMst2)*pow2(Mst1)*pow2(Mst2) + 28*(17269 + 13785*lmMst1 - 13785*
        lmMst2)*pow4(Mst1) + 27*(3477 + 490*lmMst1 - 490*lmMst2)*pow4(Mst2)))/
        315.)/pow4(Mst2)) - pow2(Mt)*pow2(s2t)*((pow2(Mst1)*(139.60072016460904
         - (1048*B4)/3. + (20*DN)/3. - (6938*lmMst1)/27. + (128*pow2(lmMst1))/
        3. - (2*lmMst2*(7619 + 720*lmMst1 + 288*pow2(lmMst1)))/27. + (8*(4 +
        153*lmMst1)*pow2(lmMst2))/3. + (800*Dmsqst2*(-1 + lmMst1 - lmMst2))/
        pow2(Msq) - (1160*pow3(lmMst2))/3.))/Mst2 + Mst2*(717.5533950617285 - (
        980*B4)/3. + (22*DN)/3. - (1648*lmMst1)/3. - (64*pow2(lmMst1))/3. + (4*
        lmMst2*(-695 - 108*lmMst1 + 24*pow2(lmMst1)))/9. + (28*(-18 + 31*
        lmMst1)*pow2(lmMst2))/3. + (80*Dmsqst2*(2 + 5*lmMst1 - 5*lmMst2))/pow2(
        Msq) - 300*pow3(lmMst2)) + (128*(2 + lmMst2 - 3*pow2(lmMst2))*pow3(
        Mst2))/(3.*pow2(Mst1)) + ((1033.594170096022 - (592*B4)/3. + (8*DN)/3.
         + (35140*lmMst1)/
        81. + (92*pow2(lmMst1))/3. - (2*lmMst2*(28667 - 5238*lmMst1 + 3024*
        pow2(lmMst1)))/81. + (8*(-60 + 157*lmMst1)*pow2(lmMst2))/3. + (20*
        Dmsqst2*(-80 + 43*lmMst1 - 43*lmMst2))/pow2(Msq) - (32*pow3(lmMst1))/3.
         - (1000*pow3(lmMst2))/3.)*pow4(Mst1))/pow3(Mst2) + ((-4*OepS2*(-9267*
        pow2(Mst1)*pow2(Mst2) + 23734*pow4(Mst1) - 63*pow4(Mst2)))/729. + (S2*(
        3*(29123 - 92670*lmMst1 + 92670*lmMst2)*pow2(Mst1)*pow2(Mst2) + 2*(
        525961 + 356010*lmMst1 - 356010*lmMst2)*pow4(Mst1) + 27*(1387 - 70*
        lmMst1 + 70*lmMst2)*pow4(Mst2)))/270.)/pow3(Mst2)) + (pow4(Mt)*((
        907200*Dmsqst2*pow2(Mst2)*((-497 + 306*lmMst1 - 609*lmMst2 + 303*lmMt)*
        pow2(Mst1) + (-101 + 90*lmMst1 - 177*lmMst2 + 87*lmMt)*pow2(Mst2)))/
        pow2(Msq) + 24*pow2(Mst1)*pow2(Mst2)*(2848907 - 10515960*lmMt + 605920*
        OepS2 + 11406852*S2 + 1890*lmMst2*(8232 - 203*lmMt + 6492*S2) - 90720*(
        -3 + lmMst2)*pow2(lmMst1) - 1890*lmMst1*(3460 - 281*lmMt + lmMst2*(689
        + 96*lmMt) + 6492*S2 - 438*pow2(lmMst2)) + 1890*(347 + 96*lmMt)*pow2(
        lmMst2) + 226800*pow2(lmMt) - 737100*pow3(lmMst2)) + 140*(5659217 +
        1592460*lmMt - 518816*OepS2 + 9976392*S2 - 972*(569 + 180*lmMst2 - 126*
        lmMt)*pow2(lmMst1) + 324*(5353 + 126*lmMt)*pow2(lmMst2) - 186624*pow2(
        lmMt) + 72*lmMst1*(-27653 - 3015*lmMt - 18*lmMst2*(689 + 126*lmMt) +
        145917*S2 + 2160*pow2(lmMst2) + 2592*pow2(lmMt)) - 36*lmMst2*(-39031 -
        3204*lmMt + 291834*S2 + 5184*pow2(lmMt)) + 1944*pow3(lmMst1) + 17496*
        pow3(lmMst2))*pow4(Mst1) - 9*(7152121 + 4616640*lmMt - 7840*OepS2 +
        17640*lmMst2*(-1256 + 60*lmMt - 9*S2) - 3032964*S2 - 241920*(-2 +
        lmMst2)*pow2(lmMst1) + 7560*lmMst1*(1136 - 128*lmMt + 16*lmMst2*(21 +
        4*lmMt) + 21*S2 - 164*pow2(lmMst2)) - 30240*(107 + 16*lmMt)*pow2(
        lmMst2) + 120960*pow2(lmMt) + 1481760*pow3(lmMst2))*pow4(Mst2) + (
        8709120*(2 + lmMst2 - 3*pow2(lmMst2))*pow6(Mst2))/pow2(Mst1)))/(153090.
        *pow5(Mst2)) + Mt*((2*T1ep*(3*pow2(Mst1)*pow2(Mst2)*(-15840*Mst2*s2t*
        pow2(Mt) + 9267*Mt*pow2(Mst2)*pow2(s2t) - 17312*pow3(Mt) + 530*pow3(
        Mst2)*pow3(s2t)) + 2*(-66168*Mst2*s2t*pow2(Mt) - 35601*Mt*pow2(Mst2)*
        pow2(s2t) + 129704*pow3(Mt) + 12664*pow3(Mst2)*pow3(s2t))*pow4(Mst1) +
        63*(-72*Mst2*s2t*pow2(Mt) + 3*Mt*pow2(Mst2)*pow2(s2t) - 4*pow3(Mt) +
        10*pow3(Mst2)*pow3(s2t))*pow4(Mst2)))/(729.*pow5(Mst2)) + pow3(s2t)*(
        pow2(Mst2)*(185.44423868312757 - (100*B4)/9. + (112*D3)/9. - (62*DN)/9.
         + (284*lmMst1)/3. - (103*pow2(lmMst1))/9. + (2*lmMst2*(-2465 - 1512*
        lmMst1 + 144*pow2(lmMst1)))/27. + (70.11111111111111 - 82*lmMst1)*pow2(
        lmMst2) + (40*Dmsqst2*(5 + 2*lmMst1 - 2*lmMst2))/(3.*pow2(Msq)) + (214*
        pow3(lmMst2))/3.) + pow2(Mst1)*(64.60436488340191 - (64*B4)/9. + (64*
        D3)/9. - (32*DN)/9. - (80549*lmMst1)/675. - (6011*pow2(lmMst1))/135. +
        (lmMst2*(41949 + 45410*lmMst1 + 58350*pow2(lmMst1)))/675. - (2*(2383 +
        11595*lmMst1)*pow2(lmMst2))/135. + (10*Dmsqst2*(-55 + 34*lmMst1 - 34*
        lmMst2))/(3.*pow2(Msq)) - (10*pow3(lmMst1))/27. + (2314*pow3(lmMst2))/
        27.) + ((407.7379592503276 - (2740559*lmMst1)/59535. + (866*pow2(
        lmMst1))/189. + lmMst2*(36.477181489879904 - (4840*lmMst1)/189. + (208*
        pow2(lmMst1))/3.) + (21.026455026455025 - (1136*lmMst1)/9.)*pow2(
        lmMst2) + (10*Dmsqst2*(-419 + 57*lmMst1 - 57*lmMst2))/(27.*pow2(Msq)) -
        (112*pow3(lmMst1))/27. + (1648*pow3(lmMst2))/27.)*pow4(Mst1))/pow2(
        Mst2) - ((101.77777777777777 + (32*lmMst1*(2 - 3*lmMst2))/9. - (130*
        lmMst2)/3. - 73*pow2(lmMst2) - (40*Dmsqst2)/pow2(Msq))*pow4(Mst2))/
        pow2(Mst1) - (159*(200*OepS2 + 27*(21401 - 150*lmMst1 + 150*lmMst2)*S2)
        *pow2(Mst1)*pow2(Mst2) + 20*(25328*OepS2 + 27*(47051 - 18996*lmMst1 +
        18996*lmMst2)*S2)*pow4(Mst1) + 9*(1400*OepS2 + 81*(116129 - 350*lmMst1
        + 350*lmMst2)*S2)*pow4(Mst2))/(10935.*pow2(Mst2)) - (32*(-2 + lmMst2 +
        5*pow2(lmMst2))*pow6(Mst2))/(9.*pow4(Mst1))))))/Tbeta - (s2t*pow3(Mt)*(
        9*Mst2*(6779161 + 4495680*lmMt + 5040*lmMst2*(-4309 + 210*lmMt) -
        241920*(-2 + lmMst2)*pow2(lmMst1) - 30240*(195 + 16*lmMt)*pow2(lmMst2)
        - 30240*lmMst1*(-4*(75 - 8*lmMt + lmMst2*(19 + 4*lmMt)) + 41*pow2(
        lmMst2)) + 120960*pow2(lmMt) + (50400*Dmsqst2*(289 - 180*lmMst1 + 354*
        lmMst2 - 174*lmMt))/pow2(Msq) + 1481760*pow3(lmMst2)) - (8709120*(2 +
        lmMst2 - 3*pow2(lmMst2))*pow3(Mst2))/pow2(Mst1) + ((-3*pow2(Mst1)*(
        21772800*Dmsqst2*(-7 + 3*lmMst1 - 6*lmMst2 + 3*lmMt) + pow2(Msq)*(
        79386499 - 72455040*lmMt + 52920*lmMst2*(1095 + 2*lmMt) - 725760*(-5 +
        2*lmMst2)*pow2(lmMst1) - 4460400*pow2(lmMst2) + 7560*lmMst1*(-3601 -
        370*lmMst2 + 178*lmMt + 384*pow2(lmMst2)) + 2177280*pow2(lmMt) -
        1451520*pow3(lmMst2))))/Mst2 + (28*(72900*Dmsqst2*(465 + lmMst2*(374 -
        12*lmMt) - 4*lmMst1*(52 + 3*lmMst2 - 3*lmMt) - 166*lmMt + 12*pow2(
        lmMst2)) - pow2(Msq)*(20964397 + 4563540*lmMt - 9720*(266 + 74*lmMst2 -
        55*lmMt)*pow2(lmMst1) + 3240*(1493 + 39*lmMt)*pow2(lmMst2) - 180*
        lmMst1*(32857 + 11844*lmMt + 18*lmMst2*(485 + 204*lmMt) - 1674*pow2(
        lmMst2) - 3888*pow2(lmMt)) + 1440*lmMst2*(1181 + 1332*lmMt - 486*pow2(
        lmMt)) - 466560*pow2(lmMt) + 9720*pow3(lmMst1) + 408240*pow3(lmMst2)))*
        pow4(Mst1))/pow3(Mst2))/pow2(Msq) - 153090*pow2(MuSUSY)*((256*Mst2*(2 +
        lmMst2 - 3*pow2(lmMst2)))/(9.*pow2(Mst1)) + (535.2578189300411 - (1960*
        B4)/9. + (44*DN)/9. - (3296*lmMst1)/9. - (128*pow2(lmMst1))/9. + (8*
        lmMst2*(-599 - 108*lmMst1 + 24*pow2(lmMst1)))/27. + (8*(-222 + 217*
        lmMst1)*pow2(lmMst2))/9. + (160*Dmsqst2*(2 + 5*lmMst1 - 5*lmMst2))/(3.*
        pow2(Msq)) - 200*pow3(lmMst2))/Mst2 + (pow2(Mst1)*(628.3249657064472 -
        (1352*B4)/3. + (28*DN)/3. - (43540*lmMst1)/81. + (128*pow2(lmMst1))/9.
         - (4*lmMst2*(11213 + 1368*lmMst1 + 144*pow2(lmMst1)))/
        81. + (8*(-214 + 523*lmMst1)*pow2(lmMst2))/9. + (160*Dmsqst2*(-8 + 15*
        lmMst1 - 15*lmMst2))/(3.*pow2(Msq)) - (4120*pow3(lmMst2))/9.))/pow3(
        Mst2) + ((634.115454961134 - (3416*B4)/9. + (52*DN)/9. + (111356*
        lmMst1)/243. + (8*pow2(lmMst1))/9. - (8*lmMst2*(36802 - 11421*lmMst1 +
        1728*pow2(lmMst1)))/243. + (344*(-14 + 15*lmMst1)*pow2(lmMst2))/9. - (
        64*pow3(lmMst1))/9. - (1528*pow3(lmMst2))/3.)*pow4(Mst1))/pow5(Mst2)) +
        (-560*OepS2*(-2*(51460*pow2(Mst2) + 16421*pow2(MuSUSY))*pow4(Mst1) +
        63*(2*pow2(Mst2) + pow2(MuSUSY))*pow4(Mst2) + 6*pow2(Mst1)*(1555*pow2(
        Mst2)*pow2(MuSUSY) + 4307*pow4(Mst2))) - 54*S2*(14*(8*(220117 + 192975*
        lmMst1 - 192975*lmMst2)*pow2(Mst2) + 5*(123113 + 98526*lmMst1 - 98526*
        lmMst2)*pow2(MuSUSY))*pow4(Mst1) + 6*pow2(Mst1)*(7*(20803 - 46650*
        lmMst1 + 46650*lmMst2)*pow2(Mst2)*pow2(MuSUSY) + (760703 - 904470*
        lmMst1 + 904470*lmMst2)*pow4(Mst2)) + 27*(7*(1387 - 70*lmMst1 + 70*
        lmMst2)*pow2(MuSUSY)*pow4(Mst2) + (18722 - 980*lmMst1 + 980*lmMst2)*
        pow6(Mst2))))/pow5(Mst2)))/153090. + pow2(Mt)*((5*Dmsqst2*(-1081 + 165*
        lmMst1 - 165*lmMst2)*pow2(s2t)*(-2 + pow2(Sbeta))*pow2(Sbeta)*pow4(
        Mst1))/(9.*pow2(Msq)*pow2(Mst2)) + (pow2(MuSUSY)*((pow2(Mt)*(24*pow2(
        Mst2)*(50134 - 270*B4 + 270*D3 - 135*DN + 2160*lmMst1 + 120*(271 + 24*
        lmMst1)*lmMst2 - 120*lmMt - 29646*S2 - 720*(-2 + 3*lmMst1)*pow2(lmMst2)
        + 2160*pow3(lmMst2)) + pow2(Mst1)*(3051661 - 12960*B4 + 12960*D3 -
        6480*DN + 304320*lmMst1 + 960*(2227 + 285*lmMst1)*lmMst2 + 2880*(4 +
        lmMst1 - lmMst2)*lmMt + 5600*OepS2 - 324*(5361 + 350*lmMst1 - 350*
        lmMst2)*S2 - 2880*(23 + 72*lmMst1)*pow2(lmMst2) + 207360*pow3(lmMst2)))
        )/(1215.*pow4(Mst2)) + pow2(s2t)*(90.77757201646091 - (100*B4)/9. + (
        112*D3)/9. - (62*DN)/9. + (788*lmMst1)/9. - (103*pow2(lmMst1))/9. + (
        16*lmMst2*(-241 - 171*lmMst1 + 18*pow2(lmMst1)))/27. + (
        125.33333333333333 - 82*lmMst1)*pow2(lmMst2) + (80*Dmsqst2*(4 + lmMst1
        - lmMst2))/(3.*pow2(Msq)) + ((360*Dmsqst2 + (-852 + 358*lmMst2 + 32*
        lmMst1*(-2 + 3*lmMst2) + 497*pow2(lmMst2))*pow2(Msq))*pow2(Mst2))/(9.*
        pow2(Msq)*pow2(Mst1)) + (214*pow3(lmMst2))/3. + (pow2(Mst1)*(
        155.38193689986284 - (164*B4)/9. + (176*D3)/9. - (94*DN)/9. - (21449*
        lmMst1)/675. - (7556*pow2(lmMst1))/135. + (lmMst2*(-54451 - 22990*
        lmMst1 + 65550*pow2(lmMst1)))/675. + (2*(6077 - 17130*lmMst1)*pow2(
        lmMst2))/135. + (10*Dmsqst2*(-23 + 42*lmMst1 - 42*lmMst2))/(3.*pow2(
        Msq)) - (10*pow3(lmMst1))/27. + (4240*pow3(lmMst2))/27.))/pow2(Mst2) +
        ((585.1892843532082 - (8*B4)/3. + (32*D3)/9. - (20*DN)/9. - (20109937*
        lmMst1)/297675. - (15886*pow2(lmMst1))/945. + lmMst2*(
        17.112243218274966 - (144628*lmMst1)/945. + (1180*pow2(lmMst1))/9.) + (
        2*(80257 - 138180*lmMst1)*pow2(lmMst2))/945. - (92*pow3(lmMst1))/27. +
        (4448*pow3(lmMst2))/27.)*pow4(Mst1))/pow4(Mst2) - (32*(-2 + lmMst2 + 5*
        pow2(lmMst2))*pow4(Mst2))/(9.*pow4(Mst1)) - (6*(7400*OepS2 + 27*(
        1089707 - 5550*lmMst1 + 5550*lmMst2)*S2)*pow2(Mst1)*pow2(Mst2) + 40*(
        17308*OepS2 + 27*(93919 - 12981*lmMst1 + 12981*lmMst2)*S2)*pow4(Mst1) +
        9*(1400*OepS2 + 81*(116129 - 350*lmMst1 + 350*lmMst2)*S2)*pow4(Mst2))/(
        10935.*pow4(Mst2))) + Mt*s2t*((160*Dmsqst2*((8 - 15*lmMst1 + 15*lmMst2)
        *pow2(Mst1) + (-2 - 5*lmMst1 + 5*lmMst2)*pow2(Mst2)))/(3.*pow2(Msq)*
        pow3(Mst2)) - (6*pow2(Mst2)*(4580489 - 3285360*B4 + 68040*DN - 3918600*
        lmMst1 + 248800*OepS2 - 108*(-20803 + 46650*lmMst1 - 46650*lmMst2)*S2 +
        103680*pow2(lmMst1) - 360*lmMst2*(11213 + 1368*lmMst1 + 144*pow2(
        lmMst1)) + 6480*(-214 + 523*lmMst1)*pow2(lmMst2) - 3337200*pow3(lmMst2)
        )*pow4(Mst1) + 9*pow2(Mst1)*(2601353 - 1058400*B4 + 23760*DN - 1779840*
        lmMst1 + 1120*OepS2 - 324*(-1387 + 70*lmMst1 - 70*lmMst2)*S2 - 69120*
        pow2(lmMst1) + 1440*lmMst2*(-599 - 108*lmMst1 + 24*pow2(lmMst1)) +
        4320*(-222 + 217*lmMst1)*pow2(lmMst2) - 972000*pow3(lmMst2))*pow4(Mst2)
        + 10*(2773621 - 1660176*B4 + 25272*DN + 2004408*lmMst1 - 525472*OepS2 +
        108*(123113 + 98526*lmMst1 - 98526*lmMst2)*S2 + 3888*pow2(lmMst1) -
        144*lmMst2*(36802 - 11421*lmMst1 + 1728*pow2(lmMst1)) + 167184*(-14 +
        15*lmMst1)*pow2(lmMst2) - 31104*pow3(lmMst1) - 2227824*pow3(lmMst2))*
        pow6(Mst1) + 1244160*(2 + lmMst2 - 3*pow2(lmMst2))*pow6(Mst2))/(43740.*
        pow2(Mst1)*pow5(Mst2)))))/pow2(Sbeta))) + (twoLoopFlag*((36*Mt*pow3(
        s2t)*(4*(3 - lmMst1*(-2 + lmMst2) - 4*lmMst2 + pow2(lmMst2))*pow2(Mst1)
        *pow2(Mst2) + (7 - 6*lmMst1 + 6*lmMst2)*pow4(Mst1) + 4*(14 - lmMst1*(-2
        + lmMst2) + pow2(lmMst2))*pow4(Mst2)))/Mst2 + (32*pow4(Mt)*(-9*(-27 -
        19*lmMst2 + 2*lmMst1*(6 + lmMst2 - lmMt) + 7*lmMt + 2*lmMst2*lmMt - 2*
        pow2(lmMst2))*pow2(Mst2)*pow4(Mst1) + (-2 + 9*lmMst1 + 3*lmMst2 + 15*
        lmMt)*pow2(Mst1)*pow4(Mst2) + 27*(23 + 38*lmMst2 - 9*lmMt - 6*lmMst2*
        lmMt + lmMst1*(-29 - 6*lmMst2 + 6*lmMt) + 6*pow2(lmMst2))*pow6(Mst1) +
        9*(-2 + 3*lmMst2)*pow6(Mst2)))/(pow2(Mst1)*pow4(Mst2)) + (8*s2t*pow3(
        Mt)*(9*(4*(-28 - 41*lmMst2 + lmMst1*(31 + 6*lmMst2 - 6*lmMt) + 10*lmMt
        + 6*lmMst2*lmMt - 6*pow2(lmMst2))*pow2(Mst2) + (-199 + 26*lmMst2 +
        lmMst1*(-34 + 20*lmMst2) - 20*pow2(lmMst2))*pow2(MuSUSY))*pow4(Mst1) +
        36*(-14 + lmMst1*(-2 + lmMst2) - pow2(lmMst2))*pow2(MuSUSY)*pow4(Mst2)
        + 36*pow2(Mst1)*((-31 + 3*lmMst1*(-2 + lmMst2) + 4*lmMst2 - 3*pow2(
        lmMst2))*pow2(Mst2)*pow2(MuSUSY) + (-13 + 2*lmMst1 - 4*lmMst2 + 2*lmMt)
        *pow4(Mst2)) + 2*(161 + 36*lmMst1*(-2 + lmMst2) + 78*lmMst2 - 6*lmMt -
        36*pow2(lmMst2))*pow6(Mst2)))/pow5(Mst2) + (12*pow2(Mt)*pow2(s2t)*(2*
        pow4(Mst1)*(9*(4 - lmMst1 + lmMst2 + 2*lmMst1*lmMst2 - 2*pow2(lmMst2))*
        pow2(Mst2)*pow2(MuSUSY) + (-7 + 6*lmMst1 - 24*lmMst2)*pow4(Mst2)) + (
        12*(-4 + 5*lmMst1 - 5*lmMst2)*pow2(Mst2) + 9*(9 + 4*lmMst1*(-1 +
        lmMst2) + 4*lmMst2 - 4*pow2(lmMst2))*pow2(MuSUSY))*pow6(Mst1) - 6*(-2 +
        3*lmMst2)*(2*pow2(Mst2) + pow2(MuSUSY))*pow6(Mst2) + 2*pow2(Mst1)*(3*(
        10 - 3*lmMst1 + 6*lmMst1*lmMst2 - 6*pow2(lmMst2))*pow2(MuSUSY)*pow4(
        Mst2) + (-23 + 36*lmMst2)*pow6(Mst2))))/(pow2(Mst1)*pow4(Mst2)) + (9*
        pow4(s2t)*(-6*(-2 + lmMst1 - 2*lmMst1*lmMst2 + 2*pow2(lmMst2))*pow2(
        Mst2)*pow4(Mst1) + 6*(-2 + lmMst1 - 2*lmMst2 - 2*lmMst1*lmMst2 + 2*
        pow2(lmMst2))*pow2(Mst1)*pow4(Mst2) + (1 + 6*lmMst1)*pow6(Mst1) + 2*(-2
        + 3*lmMst2)*pow6(Mst2)) + (4*Mt*MuSUSY*((9*Mt*MuSUSY*s2t*(-2*(4*Mt*(-31
        + 3*lmMst1*(-2 + lmMst2) + 4*lmMst2 - 3*pow2(lmMst2)) + 3*Mst2*s2t*(4 +
        lmMst2 + lmMst1*(-1 + 2*lmMst2) - 2*pow2(lmMst2)))*pow2(Mst2)*pow4(
        Mst1) - 2*(Mst2*s2t*(10 + lmMst1*(-3 + 6*lmMst2) - 6*pow2(lmMst2)) + 4*
        Mt*(-14 + lmMst1*(-2 + lmMst2) - pow2(lmMst2)))*pow2(Mst1)*pow4(Mst2) +
        (-3*Mst2*s2t*(9 + 4*lmMst1*(-1 + lmMst2) + 4*lmMst2 - 4*pow2(lmMst2)) +
        Mt*(398 + lmMst1*(68 - 40*lmMst2) - 52*lmMst2 + 40*pow2(lmMst2)))*pow6(
        Mst1) + 2*(-2 + 3*lmMst2)*s2t*pow7(Mst2)))/pow2(Sbeta) - (-2*pow2(Mst2)
        *(36*(3 - lmMst1 + lmMst2)*Mst2*s2t*pow2(Mt) + 54*Mt*(-17 + 2*lmMst1*(-
        2 + lmMst2) + 4*lmMst2 - 2*pow2(lmMst2))*pow2(Mst2)*pow2(s2t) + 4*(29 -
        18*lmMst1*(-1 + lmMst2) - 3*lmMst2 - 15*lmMt + 18*pow2(lmMst2))*pow3(
        Mt) + 9*(2 + 3*lmMst2)*pow3(Mst2)*pow3(s2t))*pow4(Mst1) + 2*pow2(Mst1)*
        (6*(-11 + 18*lmMst2)*Mst2*s2t*pow2(Mt) + 54*Mt*(14 - lmMst1*(-2 +
        lmMst2) + pow2(lmMst2))*pow2(Mst2)*pow2(s2t) + 4*(79 + 18*lmMst1*(-2 +
        lmMst2) + 39*lmMst2 - 3*lmMt - 18*pow2(lmMst2))*pow3(Mt) - 9*(8 + 3*
        lmMst2 + lmMst1*(-3 + 6*lmMst2) - 6*pow2(lmMst2))*pow3(Mst2)*pow3(s2t))
        *pow4(Mst2) - (72*(7 - 6*lmMst1 + 6*lmMst2)*Mst2*s2t*pow2(Mt) - 27*Mt*(
        75 + lmMst1*(10 - 8*lmMst2) - 10*lmMst2 + 8*pow2(lmMst2))*pow2(Mst2)*
        pow2(s2t) + 8*(236 + 339*lmMst2 - 18*lmMst1*(13 + 4*lmMst2 - 3*lmMt) -
        3*(35 + 18*lmMst2)*lmMt + 72*pow2(lmMst2))*pow3(Mt) + 27*(1 - 2*lmMst1
        + 2*lmMst2)*pow3(Mst2)*pow3(s2t))*pow6(Mst1) + 18*(-2 + 3*lmMst2)*s2t*(
        -4*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow7(Mst2))/Tbeta))/pow5(Mst2))/
        pow2(Mst1)))/108.))/pow2(Mgl) + threeLoopFlag*pow2(Al4p)*(-(pow2(Mt)*
        pow2(s2t)*((-5*Dmglst2*Dmsqst2*(-1081 + 165*lmMst1 - 165*lmMst2)*(-2 +
        pow2(Sbeta))*pow2(Sbeta)*pow4(Mst1))/(9.*Mgl*pow2(Msq)*pow2(Mst2)) + (
        pow2(MuSUSY)*(53.385802469135804 + (40*B4)/9. - (4*D3)/9. + (2*DN)/9. +
        (1672*lmMst1)/27. + (53*pow2(lmMst1))/9. - lmMst2*(129.92592592592592 -
        72*lmMst1 + (13*pow2(lmMst1))/3.) + (2*(-470 + 147*lmMst1)*pow2(lmMst2)
        )/9. + (16*pow3(lmMst1))/9. + ((-15*Dmsqst2*(288*Dmglst2*(4 + lmMst1 -
        lmMst2) + Mgl*(-1141 + 264*lmMst2 + 48*lmMst1*(-7 + 3*lmMst2) - 144*
        pow2(lmMst2))))/pow2(Msq) + Dmglst2*(-14267 + 432*B4 - 576*D3 + 360*DN
        - 4752*lmMst1 + 4536*lmMst2 + 16704*lmMst1*lmMst2 + 1404*pow2(lmMst1) -
        1152*lmMst2*pow2(lmMst1) - 20232*pow2(lmMst2) + 8856*lmMst1*pow2(
        lmMst2) - 7704*pow3(lmMst2)))/(162.*Mgl) - (271*pow3(lmMst2))/9. + (
        pow2(Mst1)*(434.270658436214 + (76*B4)/9. - (2*DN)/9. + (69088*lmMst1)/
        405. - (1313*pow2(lmMst1))/27. - (4*lmMst2*(16192 - 26430*lmMst1 +
        3465*pow2(lmMst1)))/405. + ((-5735 + 3072*lmMst1)*pow2(lmMst2))/27. + (
        Dmsqst2*(201.74098765432097 + (622*lmMst1)/15. + (2*(-311 + 10*lmMst1)*
        lmMst2)/15. + (Dmglst2*(76.66666666666667 - 140*lmMst1 + 140*lmMst2))/
        Mgl - (22*pow2(lmMst1))/3. + 6*pow2(lmMst2)))/pow2(Msq) - (62*pow3(
        lmMst1))/27. - (2086*pow3(lmMst2))/27. - (2*Dmglst2*(2695042 - 40500*B4
        + 54000*D3 - 33750*DN - 326895*lmMst1 - 324900*pow2(lmMst1) + 15*
        lmMst2*(-19607 - 129030*lmMst1 + 62550*pow2(lmMst1)) - 450*(-5023 +
        5610*lmMst1)*pow2(lmMst2) + 11250*pow3(lmMst1) + 1575000*pow3(lmMst2)))
        /(30375.*Mgl)))/pow2(Mst2) + ((628.1736268201578 + (76*B4)/9. - (2*DN)/
        9. + (6317839*lmMst1)/396900. - (66307*pow2(lmMst1))/315. - lmMst2*(
        12.52907281431091 - (182909*lmMst1)/315. + (274*pow2(lmMst1))/3.) + (2*
        (-58301 + 37135*lmMst1)*pow2(lmMst2))/315. + (Dmsqst2*(
        237.28785508324435 + (16526*lmMst2)/3969. + (2*lmMst1*(-8263 + 71820*
        lmMst2))/3969. - (520*pow2(lmMst1))/21. - (80*pow2(lmMst2))/7.))/pow2(
        Msq) - (44*pow3(lmMst1))/9. - (1256*pow3(lmMst2))/9. - (Dmglst2*(
        585.1892843532082 - (8*B4)/3. + (32*D3)/9. - (20*DN)/9. - (20109937*
        lmMst1)/297675. - (15886*pow2(lmMst1))/945. + lmMst2*(
        17.112243218274966 - (144628*lmMst1)/945. + (1180*pow2(lmMst1))/9.) + (
        2*(80257 - 138180*lmMst1)*pow2(lmMst2))/945. - (92*pow3(lmMst1))/27. +
        (4448*pow3(lmMst2))/27.))/Mgl)*pow4(Mst1))/pow4(Mst2) + (-((360*
        Dmglst2*Dmsqst2 + Dmsqst2*(30 - 60*lmMst2)*Mgl + Mgl*(103 + 186*lmMst2
        + 32*lmMst1*(1 + lmMst2) + 91*pow2(lmMst2))*pow2(Msq) + 2*Dmglst2*(-20
        + 2*(99 + 16*lmMst1)*lmMst2 + 123*pow2(lmMst2))*pow2(Msq))*pow2(Mst2))/
        (9.*pow2(Msq)*pow2(Mst1)) + (16*(1 + lmMst2)*(4*Dmglst2*lmMst2 + Mgl +
        lmMst2*Mgl)*pow4(Mst2))/(9.*pow4(Mst1)) + (-30*Dmsqst2*Mgl*(9*(104*
        OepS2 + 27*(17 - 78*lmMst1 + 78*lmMst2)*S2)*pow2(Mst1)*pow2(Mst2) + 13*
        (136*OepS2 - 27*(95 + 102*lmMst1 - 102*lmMst2)*S2)*pow4(Mst1) + 27*(8*
        OepS2 - 81*(-15 + 2*lmMst1 - 2*lmMst2)*S2)*pow4(Mst2)) + pow2(Msq)*(3*
        Mgl*(-3*(8456*OepS2 - 81*(11243 + 2114*lmMst1 - 2114*lmMst2)*S2)*pow2(
        Mst1)*pow2(Mst2) + (-52948*OepS2 + 27*(194357 + 39711*lmMst1 - 39711*
        lmMst2)*S2)*pow4(Mst1) + 27*(-184*OepS2 + 81*(307 + 46*lmMst1 - 46*
        lmMst2)*S2)*pow4(Mst2)) + 2*Dmglst2*(54*(344*OepS2 + 9*(15643 - 774*
        lmMst1 + 774*lmMst2)*S2)*pow2(Mst1)*pow2(Mst2) + 4*(17308*OepS2 + 27*(
        93919 - 12981*lmMst1 + 12981*lmMst2)*S2)*pow4(Mst1) + 27*(56*OepS2 -
        81*(-1677 + 14*lmMst1 - 14*lmMst2)*S2)*pow4(Mst2))))/(2187.*pow2(Msq)*
        pow4(Mst2)))/Mgl))/pow2(Sbeta))) + pow4(s2t)*(-(pow2(Mst1)*pow2(Mst2)*(
        79.01365226337448 - B4/9. + (2*D3)/9. - DN/6. + (4372*lmMst1)/405. + (
        lmMst2*(16051 + 22980*lmMst1 - 5175*pow2(lmMst1)))/810. - (1631*pow2(
        lmMst1))/108. + ((-92 + 327*lmMst1)*pow2(lmMst2))/27. - (Dmsqst2*(
        3.2221604938271606 + (467*lmMst1)/90. - (3.522222222222222 + 7*lmMst1)*
        lmMst2 + (5*Dmglst2*(-75 + 26*lmMst1 - 26*lmMst2))/(6.*Mgl) + (11*pow2(
        lmMst1))/6. + (31*pow2(lmMst2))/6.))/pow2(Msq) - (79*pow3(lmMst1))/54.
         + (Dmglst2*(
        0.7822304526748971 - (2*B4)/3. + (8*D3)/9. - (5*DN)/9. + (81193*lmMst1)
        /4050. + (137*pow2(lmMst1))/135. - (lmMst2*(81643 + 86970*lmMst1 +
        48150*pow2(lmMst1)))/4050. + ((4969 + 3840*lmMst1)*pow2(lmMst2))/270. -
        (5*pow3(lmMst1))/27. - (58*pow3(lmMst2))/27.))/Mgl - (115*pow3(lmMst2))
        /27.)) + (46.745471895783595 + B4 + D3/9. - DN/9. + (34838747*lmMst1)/
        529200. + (50723*pow2(lmMst1))/1890. + (lmMst2*(-23468297 - 17276980*
        lmMst1 + 3601500*pow2(lmMst1)))/529200. + (4*(2941 - 2415*lmMst1)*pow2(
        lmMst2))/945. + (Dmsqst2*(15.13649301931237 + (555521*lmMst1)/39690. -
        ((621671 + 200340*lmMst1)*lmMst2)/39690. + (5*Dmglst2*(76 - 249*lmMst1
        + 249*lmMst2))/(54.*Mgl) + (53*pow2(lmMst1))/21. + (53*pow2(lmMst2))/
        21.))/pow2(Msq) - (10*pow3(lmMst1))/27. + (409*pow3(lmMst2))/108. + (
        Dmglst2*(79.58863384550371 + (8287903*lmMst2)/1.1907e6 + (
        4.326984126984127 + (11*lmMst2)/3.)*pow2(lmMst1) + lmMst1*(
        1.2061367262954565 - (101*lmMst2)/315. - (11*pow2(lmMst2))/3.) - (51*
        pow2(lmMst2))/70. - (11*pow3(lmMst1))/9. + (11*pow3(lmMst2))/9.))/Mgl)*
        pow4(Mst1) - (S2*(3*(9*(5717 + 1286*lmMst1 - 1286*lmMst2)*pow2(Mst1)*
        pow2(Mst2) + 5*(3370 + 1077*lmMst1 - 1077*lmMst2)*pow4(Mst1) + 81*(307
        + 46*lmMst1 - 46*lmMst2)*pow4(Mst2)) - (2*Dmglst2*(36*(-275 + 324*
        lmMst1 - 324*lmMst2)*pow2(Mst1)*pow2(Mst2) + (51635 + 25194*lmMst1 -
        25194*lmMst2)*pow4(Mst1) + 81*(-1677 + 14*lmMst1 - 14*lmMst2)*pow4(
        Mst2)))/Mgl + (30*Dmsqst2*(9*(253 + 42*lmMst1 - 42*lmMst2)*pow2(Mst1)*
        pow2(Mst2) + (326 + 84*lmMst1 - 84*lmMst2)*pow4(Mst1) + 81*(-15 + 2*
        lmMst1 - 2*lmMst2)*pow4(Mst2)))/pow2(Msq)))/324. + (27*(Dmsqst2*(8640*
        Dmglst2*(1 + lmMst1 - lmMst2) + 30*Mgl*(-1213 + 408*lmMst2 + 48*lmMst1*
        (-7 + 3*lmMst2) - 144*pow2(lmMst2))) + pow2(Msq)*(-(Mgl*(25289 + 1440*
        B4 - 144*D3 + 72*DN + 22368*lmMst1 + 1908*pow2(lmMst1) - 12*lmMst2*(
        2296 - 2136*lmMst1 + 117*pow2(lmMst1)) + 504*(-53 + 21*lmMst1)*pow2(
        lmMst2) + 576*pow3(lmMst1) - 9756*pow3(lmMst2))) + 2*Dmglst2*(15707 -
        432*B4 + 576*D3 - 360*DN + 4752*lmMst1 - 1404*pow2(lmMst1) + 72*lmMst2*
        (-277 - 264*lmMst1 + 16*pow2(lmMst1)) - 72*(-142 + 123*lmMst1)*pow2(
        lmMst2) + 7704*pow3(lmMst2))))*pow4(Mst2) + 16*OepS2*(60*Dmsqst2*Mgl*(
        63*pow2(Mst1)*pow2(Mst2) + 14*pow4(Mst1) + 27*pow4(Mst2)) + pow2(Msq)*(
        -4*Dmglst2*(1944*pow2(Mst1)*pow2(Mst2) + 4199*pow4(Mst1) + 189*pow4(
        Mst2)) + 3*Mgl*(3858*pow2(Mst1)*pow2(Mst2) + 1795*pow4(Mst1) + 1242*
        pow4(Mst2)))) + (972*(360*Dmglst2*Dmsqst2 + Dmsqst2*(30 - 60*lmMst2)*
        Mgl + Mgl*(135 + 250*lmMst2 + 32*lmMst1*(1 + lmMst2) + 123*pow2(lmMst2)
        )*pow2(Msq) + 2*Dmglst2*(-20 + (262 + 32*lmMst1)*lmMst2 + 187*pow2(
        lmMst2))*pow2(Msq))*pow6(Mst2))/pow2(Mst1))/(34992.*Mgl*pow2(Msq))) - (
        s2t*pow3(Mt)*(((622080*(1 + lmMst2)*(-Dmglst2 + 3*Dmglst2*lmMst2 + Mgl
        + lmMst2*Mgl)*pow3(Mst2))/pow2(Mst1) + (27*Mst2*(2400*Dmglst2*Dmsqst2*(
        289 - 180*lmMst1 + 354*lmMst2 - 174*lmMt) + 21600*Dmsqst2*(7 + 4*lmMst1
        - 10*lmMst2 + 6*lmMt)*Mgl - 5*Dmglst2*pow2(Msq)*(30183 + lmMst2*(67184
        - 4704*lmMt) + 45408*lmMt + 2304*(-1 + lmMst2)*pow2(lmMst1) + 96*(199 +
        48*lmMt)*pow2(lmMst2) + 288*lmMst1*(-30 + 16*lmMt - 2*lmMst2*(5 + 8*
        lmMt) + 41*pow2(lmMst2)) - 14112*pow3(lmMst2)) + 9*Mgl*pow2(Msq)*(20253
        - 4160*lmMt - 80*lmMst2*(13 + 66*lmMt) - 1280*(1 + lmMst2)*pow2(lmMst1)
        + 160*(155 - 16*lmMt)*pow2(lmMst2) - 160*lmMst1*(54 + 2*lmMst2*(61 - 8*
        lmMt) - 16*lmMt + 41*pow2(lmMst2)) - 3840*pow2(lmMt) + 7840*pow3(
        lmMst2))) - (9*pow2(Mst1)*(1036800*Dmglst2*Dmsqst2*(-7 + 3*lmMst1 - 6*
        lmMst2 + 3*lmMt) + 259200*Dmsqst2*(3 - 2*lmMst1 + 4*lmMst2 - 2*lmMt)*
        Mgl + 3*Dmglst2*pow2(Msq)*(322421 + 414720*lmMt - 2880*(21 + 12*lmMst2
        - 4*lmMt)*pow2(lmMst1) + 480*(1097 - 24*lmMt)*pow2(lmMst2) - 69120*
        pow2(lmMt) - 160*lmMst2*(-2852 + 231*lmMt + 216*pow2(lmMt)) + 160*
        lmMst1*(-2969 - 2211*lmMst2 - 39*lmMt + 72*pow2(lmMst2) + 216*pow2(
        lmMt)) + 23040*pow3(lmMst2)) - Mgl*pow2(Msq)*(1156193 + 198720*lmMt -
        6480*lmMst2*(21 + 16*lmMt*(1 + lmMt)) + 8640*(-13 + 4*lmMst2 + 4*lmMt)*
        pow2(lmMst1) + 17280*(51 - 2*lmMt)*pow2(lmMst2) - 2160*lmMst1*(295 +
        360*lmMst2 - 52*lmMt + 112*pow2(lmMst2) - 48*pow2(lmMt)) + 207360*pow3(
        lmMst2))))/Mst2 + (4*(72900*Dmsqst2*(Mgl*(-35 - 50*lmMst2 + 4*lmMst1*(8
        + lmMst2 - lmMt) + 18*lmMt + 4*lmMst2*lmMt - 4*pow2(lmMst2)) + Dmglst2*
        (465 + 374*lmMst2 - 4*lmMst1*(52 + 3*lmMst2 - 3*lmMt) - 166*lmMt - 12*
        lmMst2*lmMt + 12*pow2(lmMst2))) + pow2(Msq)*(3*Mgl*(2016907 + 110160*
        lmMt + 1080*(565 - 3*lmMt)*pow2(lmMst2) - 1080*((158 - 6*lmMst2 - 39*
        lmMt)*pow2(lmMst1) + lmMst1*(421 + 60*lmMt + lmMst2*(413 + 36*lmMt) +
        129*pow2(lmMst2) - 72*pow2(lmMt)) + lmMst2*(167 - 66*lmMt + 72*pow2(
        lmMt))) + 1080*pow3(lmMst1) + 131760*pow3(lmMst2)) - Dmglst2*(20964397
        + 4563540*lmMt - 9720*(266 + 74*lmMst2 - 55*lmMt)*pow2(lmMst1) + 3240*(
        1493 + 39*lmMt)*pow2(lmMst2) - 180*lmMst1*(32857 + 11844*lmMt + 18*
        lmMst2*(485 + 204*lmMt) - 1674*pow2(lmMst2) - 3888*pow2(lmMt)) + 1440*
        lmMst2*(1181 + 1332*lmMt - 486*pow2(lmMt)) - 466560*pow2(lmMt) + 9720*
        pow3(lmMst1) + 408240*pow3(lmMst2))))*pow4(Mst1))/pow3(Mst2))/pow2(Msq)
        )/Mgl + 21870*pow2(MuSUSY)*(((128*(1 + lmMst2)*(-Dmglst2 + 3*Dmglst2*
        lmMst2 + Mgl + lmMst2*Mgl)*Mst2)/(9.*pow2(Mst1)) - (86400*Dmsqst2*(
        Dmglst2*(2 + 5*lmMst1 - 5*lmMst2) + (-1 - lmMst1 + lmMst2)*Mgl) + pow2(
        Msq)*(-45*Mgl*(9561 + 1760*B4 - 16*DN - 1984*lmMst1 - 256*pow2(lmMst1)
        - 64*lmMst2*(-214 + 73*lmMst1 + 4*pow2(lmMst1)) - 32*(-268 + 57*lmMst1)
        *pow2(lmMst2) + 2080*pow3(lmMst2)) - Dmglst2*(23917 + 188640*B4 - 3600*
        DN - 37440*lmMst1 + 11520*pow2(lmMst1) - 2880*lmMst2*(-237 + 55*lmMst1
        + 4*pow2(lmMst1)) - 1440*(-280 + 121*lmMst1)*pow2(lmMst2) + 185760*
        pow3(lmMst2))))/(1620.*Mst2*pow2(Msq)))/Mgl + (pow2(Mst1)*(
        621.9670781893004 + (1016*B4)/9. - (4*DN)/9. - (2476*lmMst1)/9. - (176*
        pow2(lmMst1))/9. + (4*lmMst2*(1427 - 1012*lmMst1 + 16*pow2(lmMst1)))/9.
         + (8*(642 - 203*lmMst1)*pow2(lmMst2))/9. + (520*pow3(lmMst2))/3. + ((
        160*Dmsqst2*(Dmglst2*(8 - 15*lmMst1 + 15*lmMst2) + (1 + 3*lmMst1 - 3*
        lmMst2)*Mgl))/(3.*pow2(Msq)) + Dmglst2*(248*B4 - 4*DN + (66761 +
        898740*lmMst2 + 2160*(11 + 4*lmMst2)*pow2(lmMst1) + 520560*pow2(lmMst2)
        - 180*lmMst1*(1141 + 1956*lmMst2 + 1986*pow2(lmMst2)) + 348840*pow3(
        lmMst2))/1215.))/Mgl))/pow3(Mst2) + ((1702429 + 257904*B4 - 648*DN -
        748656*lmMst1 + 41904*pow2(lmMst1) + 216*lmMst2*(5971 - 6106*lmMst1 +
        576*pow2(lmMst1)) + 41904*(34 - 15*lmMst1)*pow2(lmMst2) - 3456*pow3(
        lmMst1) - (1458*Dmglst2*(634.115454961134 - (3416*B4)/9. + (52*DN)/9. +
        (111356*lmMst1)/243. + (8*pow2(lmMst1))/9. - (8*lmMst2*(36802 - 11421*
        lmMst1 + 1728*pow2(lmMst1)))/243. + (344*(-14 + 15*lmMst1)*pow2(lmMst2)
        )/9. - (64*pow3(lmMst1))/9. - (1528*pow3(lmMst2))/3.))/Mgl + 507600*
        pow3(lmMst2))*pow4(Mst1))/(1458.*pow5(Mst2))) + (2*Dmglst2*(2*(8*(
        257300*OepS2 - 27*(220117 + 192975*lmMst1 - 192975*lmMst2)*S2)*pow2(
        Mst2) + 5*(131368*OepS2 - 27*(123113 + 98526*lmMst1 - 98526*lmMst2)*S2)
        *pow2(MuSUSY))*pow4(Mst1) + 27*(50*(56*OepS2 - 81*(-27 + 14*lmMst1 -
        14*lmMst2)*S2)*pow2(Mst2) + (1400*OepS2 - 81*(-453 + 350*lmMst1 - 350*
        lmMst2)*S2)*pow2(MuSUSY))*pow4(Mst2) + 18*pow2(Mst1)*(4*(5060*OepS2 -
        27*(2489 + 3795*lmMst1 - 3795*lmMst2)*S2)*pow2(Mst2)*pow2(MuSUSY) + 3*(
        18440*OepS2 - 9*(15691 + 41490*lmMst1 - 41490*lmMst2)*S2)*pow4(Mst2)))
        - 6*Mgl*(2*(8*(13700*OepS2 - 27*(16297 + 10275*lmMst1 - 10275*lmMst2)*
        S2)*pow2(Mst2) + 5*(7528*OepS2 - 27*(9185 + 5646*lmMst1 - 5646*lmMst2)*
        S2)*pow2(MuSUSY))*pow4(Mst1) + 135*(56*OepS2 + 81*(1 - 14*lmMst1 + 14*
        lmMst2)*S2)*pow2(MuSUSY)*pow4(Mst2) + 6*pow2(Mst1)*(10*(544*OepS2 - 81*
        (169 + 136*lmMst1 - 136*lmMst2)*S2)*pow2(Mst2)*pow2(MuSUSY) + (15080*
        OepS2 - 81*(4367 + 3770*lmMst1 - 3770*lmMst2)*S2)*pow4(Mst2)) + 54*(
        280*OepS2 - 81*(-79 + 70*lmMst1 - 70*lmMst2)*S2)*pow6(Mst2)))/(Mgl*
        pow5(Mst2))))/21870. + (Mt*pow3(s2t)*(583200*Dmsqst2*pow2(Mst1)*(
        Dmglst2*(4*(12 - 5*lmMst1 + 5*lmMst2)*pow2(Mst1)*pow2(Mst2) + (40 - 3*
        lmMst1 + 3*lmMst2)*pow4(Mst1) - 4*(2 + 5*lmMst1 - 5*lmMst2)*pow4(Mst2))
        + Mgl*(4*(-1 + lmMst1 - lmMst2)*pow2(Mst1)*pow2(Mst2) + (-4 + lmMst1 -
        lmMst2)*pow4(Mst1) + 4*(1 + lmMst1 - lmMst2)*pow4(Mst2))) + pow2(Msq)*(
        Dmglst2*(18*pow2(Mst2)*(27211 + 36720*B4 + 1080*DN - 298440*lmMst1 +
        64160*OepS2 - 108*(14033 + 12030*lmMst1 - 12030*lmMst2)*S2 + 12960*
        pow2(lmMst1) + 360*lmMst2*(-503 - 636*lmMst1 + 144*pow2(lmMst1)) -
        2160*(30 + 89*lmMst1)*pow2(lmMst2) + 140400*pow3(lmMst2))*pow4(Mst1) +
        27*pow2(Mst1)*(69997 + 188640*B4 - 3600*DN - 37440*lmMst1 + 5600*OepS2
        - 324*(-453 + 350*lmMst1 - 350*lmMst2)*S2 + 11520*pow2(lmMst1) - 2880*
        lmMst2*(-205 + 55*lmMst1 + 4*pow2(lmMst1)) - 1440*(-184 + 121*lmMst1)*
        pow2(lmMst2) + 185760*pow3(lmMst2))*pow4(Mst2) - (31897243 - 2491360*
        OepS2 + 90290268*S2 - 360*lmMst2*(18652 + 140139*S2) + 38880*(37 - 40*
        lmMst2)*pow2(lmMst1) + 3188160*pow2(lmMst2) + 360*lmMst1*(17410 -
        12852*lmMst2 + 140139*S2 + 11232*pow2(lmMst2)) - 311040*pow3(lmMst1) -
        2177280*pow3(lmMst2))*pow6(Mst1) + 622080*(-1 + 2*lmMst2 + 3*pow2(
        lmMst2))*pow6(Mst2)) + 15*Mgl*(6*pow2(Mst2)*(51041 + 7344*B4 + 216*DN -
        80136*lmMst1 - 2336*OepS2 + 324*(347 + 146*lmMst1 - 146*lmMst2)*S2 -
        2592*pow2(lmMst1) + 216*lmMst2*(-221 - 428*lmMst1 + 48*pow2(lmMst1)) -
        432*(-122 + 89*lmMst1)*pow2(lmMst2) + 28080*pow3(lmMst2))*pow4(Mst1) +
        27*pow2(Mst1)*(25611 + 5280*B4 - 48*DN - 5952*lmMst1 - 224*OepS2 + 324*
        (-1 + 14*lmMst1 - 14*lmMst2)*S2 - 768*pow2(lmMst1) - 192*lmMst2*(-182 +
        73*lmMst1 + 4*pow2(lmMst1)) - 96*(-236 + 57*lmMst1)*pow2(lmMst2) +
        6240*pow3(lmMst2))*pow4(Mst2) + (551987 - 14048*OepS2 + 661068*S2 -
        216*lmMst2*(46 + 1317*S2) + 864*(205 + 216*lmMst2)*pow2(lmMst1) +
        216000*pow2(lmMst2) - 216*lmMst1*(248 + 1820*lmMst2 - 1317*S2 + 1632*
        pow2(lmMst2)) - 6912*pow3(lmMst1) + 172800*pow3(lmMst2))*pow6(Mst1) +
        41472*pow2(1 + lmMst2)*pow6(Mst2)))))/(87480.*Mgl*Mst2*pow2(Msq)*pow2(
        Mst1)) + (pow2(MuSUSY)*pow3(Mt)*(1080*Dmglst2*Mt*pow2(Mst2)*(72*pow2(
        Mst2)*(180 - 2*B4 + 2*D3 - DN + 16*lmMst1 + 144*lmMst2 - 216*S2 - 16*(-
        2 + lmMst1)*pow2(lmMst2) + 16*pow3(lmMst2)) + pow2(Mst1)*(28405 - 288*
        B4 + 288*D3 - 144*DN + 5472*lmMst1 - 288*(-89 + 2*lmMst1)*lmMst2 + 576*
        (-2 + lmMst1 - lmMst2)*lmMt + 224*OepS2 - 324*(65 + 14*lmMst1 - 14*
        lmMst2)*S2 - 576*(-9 + 8*lmMst1)*pow2(lmMst2) + 4608*pow3(lmMst2))) + (
        2332800*Dmsqst2*s2t*(Dmglst2*(8 - 15*lmMst1 + 15*lmMst2)*pow2(Mst1) + (
        1 + 3*lmMst1 - 3*lmMst2)*Mgl*pow2(Mst1) - Dmglst2*(2 + 5*lmMst1 - 5*
        lmMst2)*pow2(Mst2) + (1 + lmMst1 - lmMst2)*Mgl*pow2(Mst2))*pow3(Mst2))/
        pow2(Msq) + 90*Mgl*Mt*(18*pow2(Mst1)*pow2(Mst2)*(10667 - 96*B4 + 96*D3
        - 48*DN - 3072*lmMst1 - 384*(-25 + 9*lmMst1)*lmMst2 - 384*(-1 + lmMst1
        - lmMst2)*lmMt - 224*OepS2 + 324*(-43 + 14*lmMst1 - 14*lmMst2)*S2 -
        384*(-13 + 4*lmMst1)*pow2(lmMst2) + 1536*pow3(lmMst2)) + (383185 -
        2592*B4 + 2592*D3 - 1296*DN - 187704*lmMst1 - 17440*OepS2 + 648*(-57 +
        545*lmMst1 - 545*lmMst2)*S2 - 7992*pow2(lmMst1) - 216*lmMst2*(-1733 +
        630*lmMst1 + 26*pow2(lmMst1)) - 216*(-859 + 246*lmMst1)*pow2(lmMst2) +
        3456*lmMt*(3 + 6*lmMst2 - 2*lmMst1*(3 + lmMst2) + pow2(lmMst1) + pow2(
        lmMst2)) + 720*pow3(lmMst1) + 58032*pow3(lmMst2))*pow4(Mst1) + 144*(436
        - 6*B4 + 6*D3 - 3*DN - 48*lmMst1 + (408 - 96*lmMst1)*lmMst2 + 24*lmMt -
        972*S2 - 48*(-4 + lmMst1)*pow2(lmMst2) + 48*pow3(lmMst2))*pow4(Mst2)) +
        (Mst2*s2t*(Dmglst2*(36*pow2(Mst2)*(66761 + 301320*B4 - 4860*DN -
        205380*lmMst1 + 40480*OepS2 - 216*(2489 + 3795*lmMst1 - 3795*lmMst2)*S2
        + 23760*pow2(lmMst1) + 180*lmMst2*(4993 - 1956*lmMst1 + 48*pow2(lmMst1)
        ) - 1080*(-482 + 331*lmMst1)*pow2(lmMst2) + 348840*pow3(lmMst2))*pow4(
        Mst1) + 27*pow2(Mst1)*(23917 + 188640*B4 - 3600*DN - 37440*lmMst1 +
        5600*OepS2 - 324*(-453 + 350*lmMst1 - 350*lmMst2)*S2 + 11520*pow2(
        lmMst1) - 2880*lmMst2*(-237 + 55*lmMst1 + 4*pow2(lmMst1)) - 1440*(-280
        + 121*lmMst1)*pow2(lmMst2) + 185760*pow3(lmMst2))*pow4(Mst2) - 10*(
        2773621 - 1660176*B4 + 25272*DN + 2004408*lmMst1 - 525472*OepS2 + 108*(
        123113 + 98526*lmMst1 - 98526*lmMst2)*S2 + 3888*pow2(lmMst1) - 144*
        lmMst2*(36802 - 11421*lmMst1 + 1728*pow2(lmMst1)) + 167184*(-14 + 15*
        lmMst1)*pow2(lmMst2) - 31104*pow3(lmMst1) - 2227824*pow3(lmMst2))*pow6(
        Mst1) + 622080*(-1 + 2*lmMst2 + 3*pow2(lmMst2))*pow6(Mst2)) + 15*Mgl*(
        24*pow2(Mst2)*(75569 + 13716*B4 - 54*DN - 33426*lmMst1 - 1088*OepS2 +
        162*(169 + 136*lmMst1 - 136*lmMst2)*S2 - 2376*pow2(lmMst1) + 54*lmMst2*
        (1427 - 1012*lmMst1 + 16*pow2(lmMst1)) - 108*(-642 + 203*lmMst1)*pow2(
        lmMst2) + 21060*pow3(lmMst2))*pow4(Mst1) + 27*pow2(Mst1)*(28683 + 5280*
        B4 - 48*DN - 5952*lmMst1 - 224*OepS2 + 324*(-1 + 14*lmMst1 - 14*lmMst2)
        *S2 - 768*pow2(lmMst1) - 192*lmMst2*(-214 + 73*lmMst1 + 4*pow2(lmMst1))
        - 96*(-268 + 57*lmMst1)*pow2(lmMst2) + 6240*pow3(lmMst2))*pow4(Mst2) +
        2*(1702429 + 257904*B4 - 648*DN - 748656*lmMst1 - 30112*OepS2 + 108*(
        9185 + 5646*lmMst1 - 5646*lmMst2)*S2 + 41904*pow2(lmMst1) + 216*lmMst2*
        (5971 - 6106*lmMst1 + 576*pow2(lmMst1)) - 41904*(-34 + 15*lmMst1)*pow2(
        lmMst2) - 3456*pow3(lmMst1) + 507600*pow3(lmMst2))*pow6(Mst1) + 41472*
        pow2(1 + lmMst2)*pow6(Mst2))))/pow2(Mst1)))/(43740.*Mgl*pow2(Sbeta)*
        pow6(Mst2)) + pow2(Mt)*pow2(s2t)*(pow2(Mst1)*(97.4135390946502 - (8*D3)
        /9. + (8*DN)/9. - (88984*lmMst1)/405. - (64*(lmMst1 - lmMst2)*lmMt)/3.
         - (1738*pow2(lmMst1))/
        27. + (2*lmMst2*(52322 - 4710*lmMst1 + 315*pow2(lmMst1)))/405. + (4*(
        479 - 237*lmMst1)*pow2(lmMst2))/27. - (Dmsqst2*(17.999506172839506 + (
        1468*lmMst1)/45. - (4*(517 + 180*lmMst1)*lmMst2)/45. + (20*Dmglst2*(113
        - 6*lmMst1 + 6*lmMst2))/(9.*Mgl) + 8*pow2(lmMst1) + 8*pow2(lmMst2)))/
        pow2(Msq) - (260*pow3(lmMst1))/27. - (4*Dmglst2*(15*lmMst2*(-77657 +
        10800*lmMt) + 8*(-152101 + 40500*lmMt) - 225*(-341 + 30*lmMst2)*pow2(
        lmMst1) - 365400*pow2(lmMst2) + 30*lmMst1*(-34484 + 16260*lmMst2 -
        5400*lmMt + 21825*pow2(lmMst2)) + 2250*pow3(lmMst1) - 650250*pow3(
        lmMst2)))/(30375.*Mgl) + (1166*pow3(lmMst2))/27.) + ((
        199.98139767323744 - (18614063*lmMst1)/66150. + lmMst2*(
        293.39173091458804 - (25514*lmMst1)/945. - (286*pow2(lmMst1))/9.) - (
        48143*pow2(lmMst1))/945. + (32*lmMt*(-2*lmMst1*(1 + lmMst2) + lmMst2*(2
        + lmMst2) + pow2(lmMst1)))/3. + (77.94391534391535 - (2*lmMst1)/9.)*
        pow2(lmMst2) + (Dmsqst2*(3.5881655848700444 - (470824*lmMst1)/19845. +
        (4*(117706 + 34965*lmMst1)*lmMst2)/19845. + (160*Dmglst2*(-41 + 12*
        lmMst1 - 12*lmMst2))/(27.*Mgl) - (74*pow2(lmMst1))/21. - (74*pow2(
        lmMst2))/21.))/pow2(Msq) + (190*pow3(lmMst1))/27. - (Dmglst2*(
        175.16754355781114 - (17578814*lmMst1)/33075. + lmMst2*(
        257.7056386999244 + (105592*lmMst1)/315. - (208*pow2(lmMst1))/9.) - (
        35576*pow2(lmMst1))/315. + (16*(-4376 + 2555*lmMst1)*pow2(lmMst2))/315.
         + (64*lmMt*(2 + 5*lmMst2 - lmMst1*(5 + 2*lmMst2) + pow2(lmMst1) +
        pow2(lmMst2)))/3. + (16*pow3(lmMst1))/27. - (2896*pow3(lmMst2))/27.))/
        Mgl + (674*pow3(lmMst2))/27.)*pow4(Mst1))/pow2(Mst2) + pow2(MuSUSY)*(
        53.385802469135804 + (40*B4)/9. - (4*D3)/9. + (2*DN)/9. + (1672*lmMst1)
        /27. + (53*pow2(lmMst1))/9. - lmMst2*(129.92592592592592 - 72*lmMst1 +
        (13*pow2(lmMst1))/3.) + (2*(-470 + 147*lmMst1)*pow2(lmMst2))/9. + (16*
        pow3(lmMst1))/9. - (271*pow3(lmMst2))/9. + ((-5*Dmsqst2*(288*Dmglst2*(4
        + lmMst1 - lmMst2) + Mgl*(-1141 + 264*lmMst2 + 48*lmMst1*(-7 + 3*
        lmMst2) - 144*pow2(lmMst2))))/(54.*pow2(Msq)) - (Dmglst2*(14267 - 432*
        B4 + 576*D3 - 360*DN + 4752*lmMst1 - 1404*pow2(lmMst1) + 72*lmMst2*(-63
        - 232*lmMst1 + 16*pow2(lmMst1)) + 72*(281 - 123*lmMst1)*pow2(lmMst2) +
        7704*pow3(lmMst2)))/162.)/Mgl + (pow2(Mst1)*(434.270658436214 + (76*B4)
        /9. - (2*DN)/9. + (69088*lmMst1)/405. - (1313*pow2(lmMst1))/27. - (4*
        lmMst2*(16192 - 26430*lmMst1 + 3465*pow2(lmMst1)))/405. + ((-5735 +
        3072*lmMst1)*pow2(lmMst2))/27. + (Dmsqst2*(201.74098765432097 + (622*
        lmMst1)/15. + (2*(-311 + 10*lmMst1)*lmMst2)/15. + (10*Dmglst2*(23 - 42*
        lmMst1 + 42*lmMst2))/(3.*Mgl) - (22*pow2(lmMst1))/3. + 6*pow2(lmMst2)))
        /pow2(Msq) - (62*pow3(lmMst1))/27. - (2086*pow3(lmMst2))/27. - (2*
        Dmglst2*(2695042 - 40500*B4 + 54000*D3 - 33750*DN - 326895*lmMst1 -
        324900*pow2(lmMst1) + 15*lmMst2*(-19607 - 129030*lmMst1 + 62550*pow2(
        lmMst1)) - 450*(-5023 + 5610*lmMst1)*pow2(lmMst2) + 11250*pow3(lmMst1)
        + 1575000*pow3(lmMst2)))/(30375.*Mgl)))/pow2(Mst2) + (pow2(Mst2)*(-(((
        360*Dmglst2*Dmsqst2 + Dmsqst2*(30 - 60*lmMst2)*Mgl + Mgl*(103 + 186*
        lmMst2 + 32*lmMst1*(1 + lmMst2) + 91*pow2(lmMst2))*pow2(Msq) + 2*
        Dmglst2*(-20 + 2*(99 + 16*lmMst1)*lmMst2 + 123*pow2(lmMst2))*pow2(Msq))
        *pow2(Mst1))/pow2(Msq)) + 16*(1 + lmMst2)*(4*Dmglst2*lmMst2 + Mgl +
        lmMst2*Mgl)*pow2(Mst2)))/(9.*Mgl*pow4(Mst1)) + ((628.1736268201578 + (
        76*B4)/9. - (2*DN)/9. + (6317839*lmMst1)/396900. - (66307*pow2(lmMst1))
        /315. - lmMst2*(12.52907281431091 - (182909*lmMst1)/315. + (274*pow2(
        lmMst1))/3.) + (2*(-58301 + 37135*lmMst1)*pow2(lmMst2))/315. + (
        Dmsqst2*(237.28785508324435 + (16526*lmMst2)/3969. + (2*lmMst1*(-8263 +
        71820*lmMst2))/3969. - (520*pow2(lmMst1))/21. - (80*pow2(lmMst2))/7.))/
        pow2(Msq) - (44*pow3(lmMst1))/9. - (1256*pow3(lmMst2))/9. - (Dmglst2*(
        585.1892843532082 - (8*B4)/3. + (32*D3)/9. - (20*DN)/9. - (20109937*
        lmMst1)/297675. - (15886*pow2(lmMst1))/945. + lmMst2*(
        17.112243218274966 - (144628*lmMst1)/945. + (1180*pow2(lmMst1))/9.) + (
        169.85608465608465 - (2632*lmMst1)/9.)*pow2(lmMst2) - (92*pow3(lmMst1))
        /27. + (4448*pow3(lmMst2))/27.))/Mgl)*pow4(Mst1))/pow4(Mst2)) + ((pow2(
        Mst2)*((-720*(360*Dmglst2*Dmsqst2 + Dmsqst2*(30 - 60*lmMst2)*Mgl + Mgl*
        (87 + 154*lmMst2 + 32*lmMst1*(1 + lmMst2) + 75*pow2(lmMst2))*pow2(Msq)
        + 2*Dmglst2*(-52 + 2*(67 + 16*lmMst1)*lmMst2 + 91*pow2(lmMst2))*pow2(
        Msq))*pow2(Mst2))/pow2(Mst1) - 96*Dmglst2*(3450*Dmsqst2 - pow2(Msq)*(
        11674 - 120*B4 + 120*D3 - 60*DN + 690*lmMst1 + 345*pow2(lmMst1) - 5*
        lmMst2*(-3427 - 54*lmMst1 + 96*pow2(lmMst1)) + (4515 - 480*lmMst1)*
        pow2(lmMst2) + 960*pow3(lmMst2))) + Mgl*(900*Dmsqst2*(391 + 32*lmMst1 -
        128*lmMst2) + pow2(Msq)*(1013263 + 2880*D3 - 2880*DN + 96720*lmMst1 +
        34560*lmMt + 44640*pow2(lmMst1) - 240*lmMst2*(-2267 - 90*lmMst1 + 117*
        pow2(lmMst1)) + 1440*(133 + 8*lmMst1)*pow2(lmMst2) + 11520*pow3(lmMst1)
        + 5040*pow3(lmMst2)))))/(3240.*pow2(Msq)) + (32*(1 + lmMst2)*(4*
        Dmglst2*lmMst2 + Mgl + lmMst2*Mgl)*pow6(Mst2))/(9.*pow4(Mst1)))/Mgl + (
        -300*Dmsqst2*Mgl*((4*(16*OepS2 - 27*(35 + 12*lmMst1 - 12*lmMst2)*S2)*
        pow2(Mst2) + 13*(136*OepS2 - 27*(95 + 102*lmMst1 - 102*lmMst2)*S2)*
        pow2(MuSUSY))*pow4(Mst1) + 27*((8*OepS2 - 81*(9 + 2*lmMst1 - 2*lmMst2)*
        S2)*pow2(Mst2) + (8*OepS2 - 81*(-15 + 2*lmMst1 - 2*lmMst2)*S2)*pow2(
        MuSUSY))*pow4(Mst2) + 9*pow2(Mst1)*((104*OepS2 + 27*(17 - 78*lmMst1 +
        78*lmMst2)*S2)*pow2(Mst2)*pow2(MuSUSY) + 2*(8*OepS2 - 27*(23 + 6*lmMst1
        - 6*lmMst2)*S2)*pow4(Mst2))) + pow2(Msq)*(-3*Mgl*(6*pow2(Mst1)*pow2(
        Mst2)*((25160*OepS2 - 81*(9191 + 6290*lmMst1 - 6290*lmMst2)*S2)*pow2(
        Mst2) + 5*(8456*OepS2 - 81*(11243 + 2114*lmMst1 - 2114*lmMst2)*S2)*
        pow2(MuSUSY)) + ((274280*OepS2 - 27*(399127 + 205710*lmMst1 - 205710*
        lmMst2)*S2)*pow2(Mst2) + 10*(52948*OepS2 - 27*(194357 + 39711*lmMst1 -
        39711*lmMst2)*S2)*pow2(MuSUSY))*pow4(Mst1) + 27*((2120*OepS2 - 81*(-141
        + 530*lmMst1 - 530*lmMst2)*S2)*pow2(Mst2) + 10*(184*OepS2 - 81*(307 +
        46*lmMst1 - 46*lmMst2)*S2)*pow2(MuSUSY))*pow4(Mst2)) + 4*Dmglst2*(54*
        pow2(Mst1)*pow2(Mst2)*((2000*OepS2 - 162*(31 + 250*lmMst1 - 250*lmMst2)
        *S2)*pow2(Mst2) + 5*(344*OepS2 + 9*(15643 - 774*lmMst1 + 774*lmMst2)*
        S2)*pow2(MuSUSY)) + 2*(9*(30760*OepS2 - 27*(28283 + 23070*lmMst1 -
        23070*lmMst2)*S2)*pow2(Mst2) + 10*(17308*OepS2 + 27*(93919 - 12981*
        lmMst1 + 12981*lmMst2)*S2)*pow2(MuSUSY))*pow4(Mst1) + 135*(56*OepS2 -
        81*(-1677 + 14*lmMst1 - 14*lmMst2)*S2)*pow2(MuSUSY)*pow4(Mst2) -
        2768742*S2*pow6(Mst2))))/(21870.*Mgl*pow2(Msq)*pow4(Mst2))) + (4*
        Dmglst2*Mst2*T1ep*pow2(Msq)*(189*pow4(Mst2)*(-4*Mst2*pow2(Mt)*pow2(
        MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta)) + 10*s2t*(pow2(MuSUSY)*(-1 + pow2(
        Sbeta)) + 2*pow2(Mst2)*pow2(Sbeta))*pow3(Mt) - 5*Mt*pow2(Sbeta)*pow3(
        s2t)*pow4(Mst2) - 16*Mst2*pow2(Sbeta)*pow4(Mt) + pow2(Sbeta)*pow4(s2t)*
        pow5(Mst2)) + pow4(Mst1)*(-8*Mst2*pow2(Mt)*pow2(s2t)*(4327*pow2(MuSUSY)
        *(-1 + pow2(Sbeta)) + 6921*pow2(Mst2)*pow2(Sbeta)) + 4*s2t*(16421*pow2(
        MuSUSY)*(-1 + pow2(Sbeta)) + 51460*pow2(Mst2)*pow2(Sbeta))*pow3(Mt) -
        15571*Mt*pow2(Sbeta)*pow3(s2t)*pow4(Mst2) - 89312*Mst2*pow2(Sbeta)*
        pow4(Mt) + 4199*pow2(Sbeta)*pow4(s2t)*pow5(Mst2)) + 18*Mst2*pow2(Mst1)*
        (-12*pow2(Mst2)*pow2(Mt)*pow2(s2t)*(43*pow2(MuSUSY)*(-1 + pow2(Sbeta))
        + 50*pow2(Mst2)*pow2(Sbeta)) + 2*Mst2*s2t*(506*pow2(MuSUSY)*(-1 + pow2(
        Sbeta)) + 1383*pow2(Mst2)*pow2(Sbeta))*pow3(Mt) + 56*(3*pow2(MuSUSY)*(-
        1 + pow2(Sbeta)) - 23*pow2(Mst2)*pow2(Sbeta))*pow4(Mt) - 401*Mt*pow2(
        Sbeta)*pow3(s2t)*pow5(Mst2) + 108*pow2(Sbeta)*pow4(s2t)*pow6(Mst2))) -
        3*Mgl*T1ep*(20*Dmsqst2*pow2(Mst2)*pow2(s2t)*(27*pow4(Mst2)*(-4*pow2(Mt)
        *(pow2(MuSUSY)*(-1 + pow2(Sbeta)) + pow2(Mst2)*pow2(Sbeta)) + pow2(s2t)
        *pow2(Sbeta)*pow4(Mst2)) + 9*pow2(Mst1)*pow2(Mst2)*(pow2(Mt)*(-52*pow2(
        MuSUSY)*(-1 + pow2(Sbeta)) - 8*pow2(Mst2)*pow2(Sbeta)) + 7*pow2(s2t)*
        pow2(Sbeta)*pow4(Mst2)) + pow4(Mst1)*(-4*pow2(Mt)*(221*pow2(MuSUSY)*(-1
        + pow2(Sbeta)) + 8*pow2(Mst2)*pow2(Sbeta)) + 14*pow2(s2t)*pow2(Sbeta)*
        pow4(Mst2))) + pow2(Msq)*(54*pow5(Mst2)*(-2*Mst2*pow2(Mt)*pow2(s2t)*(
        46*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 53*pow2(Mst2)*pow2(Sbeta)) + 28*
        s2t*(pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 2*pow2(Mst2)*pow2(Sbeta))*pow3(
        Mt) - 14*Mt*pow2(Sbeta)*pow3(s2t)*pow4(Mst2) + 56*Mst2*pow2(Sbeta)*
        pow4(Mt) + 23*pow2(Sbeta)*pow4(s2t)*pow5(Mst2)) + 6*pow2(Mst1)*pow2(
        Mst2)*(-4*pow2(Mst2)*pow2(Mt)*pow2(s2t)*(1057*pow2(MuSUSY)*(-1 + pow2(
        Sbeta)) + 629*pow2(Mst2)*pow2(Sbeta)) + 8*Mst2*s2t*(136*pow2(MuSUSY)*(-
        1 + pow2(Sbeta)) + 377*pow2(Mst2)*pow2(Sbeta))*pow3(Mt) + 8*(126*pow2(
        MuSUSY)*(-1 + pow2(Sbeta)) + 209*pow2(Mst2)*pow2(Sbeta))*pow4(Mt) -
        292*Mt*pow2(Sbeta)*pow3(s2t)*pow5(Mst2) + 643*pow2(Sbeta)*pow4(s2t)*
        pow6(Mst2)) + pow4(Mst1)*(-4*pow2(Mst2)*pow2(Mt)*pow2(s2t)*(13237*pow2(
        MuSUSY)*(-1 + pow2(Sbeta)) + 6857*pow2(Mst2)*pow2(Sbeta)) + 16*Mst2*
        s2t*(941*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 2740*pow2(Mst2)*pow2(Sbeta))
        *pow3(Mt) + 16*(1635*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 1168*pow2(Mst2)*
        pow2(Sbeta))*pow4(Mt) - 1756*Mt*pow2(Sbeta)*pow3(s2t)*pow5(Mst2) +
        1795*pow2(Sbeta)*pow4(s2t)*pow6(Mst2)))))/(1458.*Mgl*pow2(Msq)*pow2(
        Sbeta)*pow6(Mst2)) + pow4(Mt)*(480.98395061728394 - (640*lmMst1)/9. - (
        392*pow2(lmMst1))/9. + (4*lmMst2*(-553 - 1224*lmMst1 + 63*pow2(lmMst1))
        )/81. + (32*(121 - 18*lmMst1)*pow2(lmMst2))/27. + (4*lmMt*(926 - 749*
        lmMst2 + 6*lmMst1*(27 + 32*lmMst2) + 9*pow2(lmMst1) + 57*pow2(lmMst2)))
        /27. - (4*(-12 + lmMst1 + 55*lmMst2)*pow2(lmMt))/3. - (68*pow3(lmMst1))
        /9. + (4*((-350*Dmsqst2*(Dmglst2*(-2954 + 864*lmMst1 - 1713*lmMst2 +
        849*lmMt) - 3*Mgl*(-103 + 42*lmMst1 + 78*lmMt - 6*lmMst2*(23 + 3*lmMt)
        + 9*pow2(lmMst2) + 9*pow2(lmMt))))/pow2(Msq) + Dmglst2*(773533 +
        131775*lmMt + 35*lmMst2*(-6569 + 2220*lmMt) - 4410*pow2(lmMst1) -
        199920*pow2(lmMst2) + 5040*lmMst1*(13 + 28*lmMst2 - 8*lmMt + 8*pow2(
        lmMst2)) + 63000*pow2(lmMt) - 40320*pow3(lmMst2))))/(2835.*Mgl) + (8*
        pow3(lmMst2))/9. + (184*pow3(lmMt))/3. + ((1500.0856244066115 + (
        43574647*lmMst1)/99225. - (139432*pow2(lmMst1))/945. - (2*lmMst2*(
        35585111 - 2612085*lmMst1 + 2072700*pow2(lmMst1)))/99225. + (2*(34339 +
        46620*lmMst1)*pow2(lmMst2))/945. + (lmMt*(3155 + 1418*lmMst2 - 2*
        lmMst1*(513 + 572*lmMst2) + 312*pow2(lmMst1) + 832*pow2(lmMst2)))/9. +
        (64*(-1 + 3*lmMst1 - 3*lmMst2)*pow2(lmMt))/3. + (8*Dmsqst2*(-26331136 +
        105*lmMst1*(236317 + 35070*lmMst2 - 36750*lmMt) + 11962125*lmMt + 210*
        lmMst2*(-175121 + 18375*lmMt) + 88200*pow2(lmMst1) - 3770550*pow2(
        lmMst2)))/(231525.*pow2(Msq)) + (64*pow3(lmMst1))/27. - (1600*pow3(
        lmMst2))/27. - (Dmglst2*(2929.938520304849 + (55510684*lmMst1)/59535. -
        (126272*pow2(lmMst1))/189. - (4*lmMst2*(42300121 + 12578580*lmMst1 +
        2487240*pow2(lmMst1)))/59535. + (32*(10166 - 693*lmMst1)*pow2(lmMst2))/
        189. + (8*lmMt*(5695 + 1974*lmMst2 - 12*lmMst1*(163 + 47*lmMst2) + 468*
        pow2(lmMst1) + 96*pow2(lmMst2)))/27. + (128*(-5 + 6*lmMst1 - 6*lmMst2)*
        pow2(lmMt))/3. + (256*pow3(lmMst1))/27. + (7424*pow3(lmMst2))/27.))/
        Mgl)*pow4(Mst1))/pow4(Mst2) + (((4*(360*Dmglst2*Dmsqst2 + Dmsqst2*(30 -
        60*lmMst2)*Mgl + 2*Dmglst2*(-84 + (70 + 32*lmMst1)*lmMst2 + 59*pow2(
        lmMst2))*pow2(Msq) + Mgl*(71 + 122*lmMst2 + 32*lmMst1*(1 + lmMst2) +
        59*pow2(lmMst2))*pow2(Msq))*pow2(Mst2))/(9.*pow2(Mst1)) + (pow2(Mst1)*(
        2520*Dmsqst2*Mgl*(-91321 + 30*lmMst1*(1717 + 120*lmMst2 - 150*lmMt) +
        40500*lmMt + 30*lmMst2*(-3067 + 150*lmMt) + 450*pow2(lmMst1) - 4050*
        pow2(lmMst2)) + 5670000*Dmglst2*Dmsqst2*(401 + lmMst2*(274 - 4*lmMt) -
        2*lmMst1*(71 + 2*lmMst2 - 2*lmMt) - 132*lmMt + 4*pow2(lmMst2)) + 4*
        Dmglst2*pow2(Msq)*(98884301 - 57865500*lmMt + 12600*(1213 + 285*lmMst2
        - 225*lmMt)*pow2(lmMst1) + 18900*(-5953 + 310*lmMt)*pow2(lmMst2) +
        13608000*pow2(lmMt) + 1260*(lmMst1*(28194 + lmMst2*(61415 - 2400*lmMt)
        - 1075*lmMt + 13800*pow2(lmMst2) - 7200*pow2(lmMt)) + lmMst2*(8581 +
        6025*lmMt + 7200*pow2(lmMt))) - 252000*pow3(lmMst1) - 20727000*pow3(
        lmMst2)) + 35*Mgl*pow2(Msq)*(7319011 + 3223800*lmMt - 1800*(391 + 114*
        lmMst2 - 90*lmMt)*pow2(lmMst1) + 3600*(568 + 99*lmMt)*pow2(lmMst2) -
        259200*pow2(lmMt) - 120*lmMst2*(22352 - 2835*lmMt + 4320*pow2(lmMt)) +
        120*lmMst1*(4217 + 1215*lmMt - 15*lmMst2*(871 + 288*lmMt) + 3240*pow2(
        lmMst2) + 4320*pow2(lmMt)) + 14400*pow3(lmMst1) - 198000*pow3(lmMst2)))
        )/(425250.*pow2(Mst2)))/pow2(Msq) - (64*(1 + lmMst2)*(4*Dmglst2*lmMst2
        + Mgl + lmMst2*Mgl)*pow4(Mst2))/(9.*pow4(Mst1)) + (-(pow2(MuSUSY)*(12*
        Dmglst2*pow2(Mst2)*(72*pow2(Mst2)*(180 - 2*B4 + 2*D3 - DN + 16*lmMst1 +
        144*lmMst2 - 16*(-2 + lmMst1)*pow2(lmMst2) + 16*pow3(lmMst2)) + pow2(
        Mst1)*(28405 - 288*B4 + 288*D3 - 144*DN + 5472*lmMst1 - 288*(-89 + 2*
        lmMst1)*lmMst2 + 576*(-2 + lmMst1 - lmMst2)*lmMt - 576*(-9 + 8*lmMst1)*
        pow2(lmMst2) + 4608*pow3(lmMst2))) + Mgl*(18*pow2(Mst1)*pow2(Mst2)*(
        10667 - 96*B4 + 96*D3 - 48*DN - 3072*lmMst1 - 384*(-25 + 9*lmMst1)*
        lmMst2 - 384*(-1 + lmMst1 - lmMst2)*lmMt - 384*(-13 + 4*lmMst1)*pow2(
        lmMst2) + 1536*pow3(lmMst2)) + (383185 - 2592*B4 + 2592*D3 - 1296*DN -
        187704*lmMst1 - 7992*pow2(lmMst1) - 216*lmMst2*(-1733 + 630*lmMst1 +
        26*pow2(lmMst1)) - 216*(-859 + 246*lmMst1)*pow2(lmMst2) + 3456*lmMt*(3
        + 6*lmMst2 - 2*lmMst1*(3 + lmMst2) + pow2(lmMst1) + pow2(lmMst2)) +
        720*pow3(lmMst1) + 58032*pow3(lmMst2))*pow4(Mst1) + 144*(436 - 6*B4 +
        6*D3 - 3*DN - 48*lmMst1 + (408 - 96*lmMst1)*lmMst2 + 24*lmMt - 48*(-4 +
        lmMst1)*pow2(lmMst2) + 48*pow3(lmMst2))*pow4(Mst2))))/486. + (16*OepS2*
        (4*Dmglst2*pow2(Mst2)*(63*pow2(Mst1)*(23*pow2(Mst2) - 3*pow2(MuSUSY)) +
        5582*pow4(Mst1) + 189*pow4(Mst2)) + 3*Mgl*((1168*pow2(Mst2) + 1635*
        pow2(MuSUSY))*pow4(Mst1) + pow2(Mst1)*(378*pow2(Mst2)*pow2(MuSUSY) +
        627*pow4(Mst2)) + 189*pow6(Mst2))))/2187. - (S2*(4*Dmglst2*pow2(Mst2)*(
        45*pow2(Mst1)*(23*(115 + 588*lmMst1 - 588*lmMst2)*pow2(Mst2) - 126*(65
        + 14*lmMst1 - 14*lmMst2)*pow2(MuSUSY)) + (2001242 + 2344440*lmMst1 -
        2344440*lmMst2)*pow4(Mst1) - 81*(3360*pow2(Mst2)*pow2(MuSUSY) + (1593 -
        980*lmMst1 + 980*lmMst2)*pow4(Mst2))) + 42*Mgl*((8*(8653 + 4380*lmMst1
        - 4380*lmMst2)*pow2(Mst2) + 90*(-57 + 545*lmMst1 - 545*lmMst2)*pow2(
        MuSUSY))*pow4(Mst1) + 9*pow2(Mst1)*(90*(-43 + 14*lmMst1 - 14*lmMst2)*
        pow2(Mst2)*pow2(MuSUSY) + (4163 + 2090*lmMst1 - 2090*lmMst2)*pow4(Mst2)
        ) + 81*(-240*pow2(MuSUSY)*pow4(Mst2) + (109 + 70*lmMst1 - 70*lmMst2)*
        pow6(Mst2)))))/2835.)/pow6(Mst2))/Mgl) + (MuSUSY*(-(pow2(Mt)*pow2(s2t)*
        (Mst2*(377.0416666666667 + (220*B4)/3. - (2*DN)/3. - (248*lmMst1)/3. -
        (32*pow2(lmMst1))/3. - (8*lmMst2*(-198 + 73*lmMst1 + 4*pow2(lmMst1)))/
        3. + (336 - 76*lmMst1)*pow2(lmMst2) + (260*pow3(lmMst2))/3. + ((-80*
        Dmsqst2*(Dmglst2*(2 + 5*lmMst1 - 5*lmMst2) + (-1 - lmMst1 + lmMst2)*
        Mgl))/pow2(Msq) + Dmglst2*(43.4787037037037 + (524*B4)/3. - (10*DN)/3.
         - (104*lmMst1)/
        3. + (32*pow2(lmMst1))/3. - (8*lmMst2*(-221 + 55*lmMst1 + 4*pow2(
        lmMst1)))/3. + (4*(232 - 121*lmMst1)*pow2(lmMst2))/3. + 172*pow3(
        lmMst2)))/Mgl) + (pow2(Mst1)*(534.5756172839506 + 96*B4 - 330*lmMst1 -
        (56*pow2(lmMst1))/3. + (2*lmMst2*(571 - 720*lmMst1 + 32*pow2(lmMst1)))/
        3. + (8*(187 - 73*lmMst1)*pow2(lmMst2))/3. + (520*pow3(lmMst2))/3. + ((
        160*Dmsqst2*(Dmglst2*(5 - 5*lmMst1 + 5*lmMst2) + (lmMst1 - lmMst2)*Mgl)
        )/pow2(Msq) + Dmglst2*(60.275617283950616 + (592*B4)/3. - (8*DN)/3. - (
        1970*lmMst1)/9. + (56*pow2(lmMst1))/3. + (2*lmMst2*(2149 - 1296*lmMst1
        + 96*pow2(lmMst1)))/9. + (269.3333333333333 - 280*lmMst1)*pow2(lmMst2)
        + (776*pow3(lmMst2))/3.))/Mgl))/Mst2 + ((818.5195473251028 + 96*B4 - (
        3218*lmMst1)/9. + (652*pow2(lmMst1))/9. + (4*lmMst2*(845 - 1535*lmMst1
        + 264*pow2(lmMst1)))/9. + (609.7777777777778 - 376*lmMst1)*pow2(lmMst2)
        - (32*pow3(lmMst1))/9. + ((20*Dmsqst2*(Dmglst2*(80 - 43*lmMst1 + 43*
        lmMst2) + (-4 + 9*lmMst1 - 9*lmMst2)*Mgl))/pow2(Msq) - Dmglst2*(
        1033.594170096022 - (592*B4)/3. + (8*DN)/3. + (35140*lmMst1)/81. + (92*
        pow2(lmMst1))/3. - (2*lmMst2*(28667 - 5238*lmMst1 + 3024*pow2(lmMst1)))
        /81. + (8*(-60 + 157*lmMst1)*pow2(lmMst2))/3. - (32*pow3(lmMst1))/3. -
        (1000*pow3(lmMst2))/3.))/Mgl + (2360*pow3(lmMst2))/9.)*pow4(Mst1))/
        pow3(Mst2) - S2*((4.5 - 63*lmMst1 + 63*lmMst2 + (3*Dmglst2*(-453 + 350*
        lmMst1 - 350*lmMst2))/(10.*Mgl))*Mst2 - ((342.5 + 209*lmMst1 - 209*
        lmMst2 - (Dmglst2*(799.6333333333333 + 907*lmMst1 - 907*lmMst2))/Mgl)*
        pow2(Mst1))/Mst2 + ((Dmglst2*(525961 + 356010*lmMst1 - 356010*lmMst2) -
        15*(6143 + 3198*lmMst1 - 3198*lmMst2)*Mgl)*pow4(Mst1))/(135.*Mgl*pow3(
        Mst2))) + ((64*(1 + lmMst2)*(-Dmglst2 + 3*Dmglst2*lmMst2 + Mgl +
        lmMst2*Mgl)*pow3(Mst2))/(3.*pow2(Mst1)) + (4*OepS2*(-3*Mgl*(627*pow2(
        Mst1)*pow2(Mst2) + 1066*pow4(Mst1) + 189*pow4(Mst2)) + Dmglst2*(8163*
        pow2(Mst1)*pow2(Mst2) + 23734*pow4(Mst1) + 945*pow4(Mst2))))/(729.*
        pow3(Mst2)))/Mgl)) - s2t*pow3(Mt)*(604.5820987654321 + (16*D3)/9. - (
        16*DN)/9. + (1228*lmMst1)/27. + (64*lmMt)/3. + (248*pow2(lmMst1))/9. -
        (4*lmMst2*(-1901 + 6*lmMst1 + 117*pow2(lmMst1)))/27. + (92 + (64*
        lmMst1)/9.)*pow2(lmMst2) + (64*pow3(lmMst1))/9. + (28*pow3(lmMst2))/9.
         + ((-5*Dmsqst2*(656*Dmglst2 + (-367 - 32*lmMst1 + 80*lmMst2)*Mgl))/(
        9.*pow2(Msq)) + (8*Dmglst2*(12454 - 120*B4 + 120*D3 - 60*DN + 690*
        lmMst1 + 345*pow2(lmMst1) - 5*lmMst2*(-3121 + 42*lmMst1 + 96*pow2(
        lmMst1)) + (3630 - 480*lmMst1)*pow2(lmMst2) + 960*pow3(lmMst2)))/135.)/
        Mgl + ((1199.3719723012073 - (99165049*lmMst1)/99225. + lmMst2*(
        1427.8402519526328 - (115988*lmMst1)/945. - (700*pow2(lmMst1))/9.) - (
        181826*pow2(lmMst1))/945. + (400.4804232804233 - (572*lmMst1)/9.)*pow2(
        lmMst2) + (64*lmMt*(1 + 4*lmMst2 - 2*lmMst1*(2 + lmMst2) + pow2(lmMst1)
        + pow2(lmMst2)))/3. + (Dmsqst2*(175.06620771294996 + (1883624*lmMst2)/
        19845. + (8*lmMst1*(-235453 + 114345*lmMst2))/19845. - (484*pow2(
        lmMst1))/21. - (484*pow2(lmMst2))/21.))/pow2(Msq) + (52*pow3(lmMst1))/
        27. + (3764*pow3(lmMst2))/27. + (Dmglst2*(622.8225754358181 - (64*B4)/
        9. + (64*D3)/9. - (32*DN)/9. + (29853268*lmMst1)/19845. + (38704*pow2(
        lmMst1))/189. + (4*lmMst2*(2917823 - 3813600*lmMst1 + 97020*pow2(
        lmMst1)))/19845. + (16*(8677 - 5439*lmMst1)*pow2(lmMst2))/189. - (128*
        lmMt*(4 + 6*lmMst2 - 2*lmMst1*(3 + lmMst2) + pow2(lmMst1) + pow2(
        lmMst2)))/3. - (16*pow3(lmMst1))/9. + (1328*pow3(lmMst2))/3.))/Mgl)*
        pow4(Mst1))/pow4(Mst2) + (((-4*(360*Dmglst2*Dmsqst2 + Dmsqst2*(30 - 60*
        lmMst2)*Mgl + 2*Dmglst2*(-52 + 2*(51 + 16*lmMst1)*lmMst2 + 59*pow2(
        lmMst2))*pow2(Msq) + Mgl*(71 + 122*lmMst2 + 32*lmMst1*(1 + lmMst2) +
        59*pow2(lmMst2))*pow2(Msq))*pow2(Mst2))/(9.*pow2(Mst1)) - (pow2(Mst1)*(
        32*Dmglst2*(50625*Dmsqst2*(65 - 2*lmMst1 + 2*lmMst2) - 2*pow2(Msq)*(
        1928479 - 13500*B4 + 13500*D3 - 6750*DN + 635385*lmMst1 + 81000*(-2 +
        lmMst1 - lmMst2)*lmMt + 450*pow2(lmMst1) - 15*lmMst2*(-153166 + 17835*
        lmMst1 + 3375*pow2(lmMst1)) - 225*(-2627 + 1695*lmMst1)*pow2(lmMst2) -
        1125*pow3(lmMst1) + 433125*pow3(lmMst2))) - 5*Mgl*(12*Dmsqst2*(339977 +
        96120*lmMst2 + 1080*lmMst1*(-89 + 60*lmMst2) - 32400*pow2(lmMst1) -
        32400*pow2(lmMst2)) + pow2(Msq)*(19425643 + 518400*lmMt + 240*lmMst2*(
        82997 + 4320*lmMt) - 3600*(683 + 96*lmMst2)*pow2(lmMst1) + 5684400*
        pow2(lmMst2) - 240*lmMst1*(42047 + 4800*lmMst2 + 4320*lmMt + 6390*pow2(
        lmMst2)) - 295200*pow3(lmMst1) + 2174400*pow3(lmMst2)))))/(121500.*
        pow2(Mst2)))/pow2(Msq) + (64*(1 + lmMst2)*(4*Dmglst2*lmMst2 + Mgl +
        lmMst2*Mgl)*pow4(Mst2))/(9.*pow4(Mst1)) + (-100*Dmsqst2*Mgl*(9*(40*
        OepS2 - 27*(127 + 30*lmMst1 - 30*lmMst2)*S2)*pow2(Mst1)*pow2(Mst2) + (
        424*OepS2 - 27*(1283 + 318*lmMst1 - 318*lmMst2)*S2)*pow4(Mst1) + 27*(8*
        OepS2 - 81*(9 + 2*lmMst1 - 2*lmMst2)*S2)*pow4(Mst2)) + pow2(Msq)*(24*
        Dmglst2*(15*(400*OepS2 - 81*(139 + 100*lmMst1 - 100*lmMst2)*S2)*pow2(
        Mst1)*pow2(Mst2) + (36760*OepS2 - 54*(17269 + 13785*lmMst1 - 13785*
        lmMst2)*S2)*pow4(Mst1) - 153819*S2*pow4(Mst2)) - Mgl*(3*(69400*OepS2 -
        81*(17113 + 17350*lmMst1 - 17350*lmMst2)*S2)*pow2(Mst1)*pow2(Mst2) + 4*
        (120620*OepS2 - 27*(138286 + 90465*lmMst1 - 90465*lmMst2)*S2)*pow4(
        Mst1) + 27*(2120*OepS2 - 81*(-141 + 530*lmMst1 - 530*lmMst2)*S2)*pow4(
        Mst2))))/(3645.*pow2(Msq)*pow4(Mst2)))/Mgl) + (pow4(Mt)*(129600*
        Dmglst2*Dmsqst2*pow2(Mst1)*pow2(Mst2)*((497 - 306*lmMst1 + 609*lmMst2 -
        303*lmMt)*pow2(Mst1) + (101 - 90*lmMst1 + 177*lmMst2 - 87*lmMt)*pow2(
        Mst2)) + 1166400*Dmsqst2*Mgl*pow2(Mst1)*pow2(Mst2)*((1 + 6*lmMst1 - 13*
        lmMst2 + 7*lmMt)*pow2(Mst1) + (5 + 2*lmMst1 - 5*lmMst2 + 3*lmMt)*pow2(
        Mst2)) - Dmglst2*pow2(Msq)*(216*pow2(Mst2)*(97837 + 71580*lmMt - 9920*
        OepS2 + 43272*S2 - 360*(23 + 8*lmMst2 - 4*lmMt)*pow2(lmMst1) + 720*(97
        + 2*lmMt)*pow2(lmMst2) - 8640*pow2(lmMt) - 30*lmMst2*(-2911 + 396*lmMt
        + 6696*S2 + 144*pow2(lmMt)) + 10*lmMst1*(-6157 + 642*lmMt - 6*lmMst2*(
        791 + 48*lmMt) + 20088*S2 + 882*pow2(lmMst2) + 432*pow2(lmMt)) - 5940*
        pow3(lmMst2))*pow4(Mst1) + 135*pow2(Mst1)*(57495 + 45408*lmMt - 1120*
        OepS2 - 43740*S2 - 8*lmMst2*(-6952 + 588*lmMt + 2835*S2) + 2304*(-1 +
        lmMst2)*pow2(lmMst1) + 96*(79 + 48*lmMt)*pow2(lmMst2) + 72*lmMst1*(-88
        + 64*lmMt - 8*lmMst2*(9 + 8*lmMt) + 315*S2 + 164*pow2(lmMst2)) - 14112*
        pow3(lmMst2))*pow4(Mst2) + 20*(5659217 + 1592460*lmMt - 518816*OepS2 +
        9976392*S2 - 972*(569 + 180*lmMst2 - 126*lmMt)*pow2(lmMst1) + 324*(5353
        + 126*lmMt)*pow2(lmMst2) - 186624*pow2(lmMt) + 72*lmMst1*(-27653 -
        3015*lmMt - 18*lmMst2*(689 + 126*lmMt) + 145917*S2 + 2160*pow2(lmMst2)
        + 2592*pow2(lmMt)) - 36*lmMst2*(-39031 - 3204*lmMt + 291834*S2 + 5184*
        pow2(lmMt)) + 1944*pow3(lmMst1) + 17496*pow3(lmMst2))*pow6(Mst1) -
        622080*(-1 + 2*lmMst2 + 3*pow2(lmMst2))*pow6(Mst2)) + 3*Mgl*pow2(Msq)*(
        48*pow2(Mst2)*(110219 - 1080*lmMt - 4400*OepS2 + 74034*S2 + 540*(-15 +
        4*lmMt)*pow2(lmMst1) + 810*(121 - 8*lmMt)*pow2(lmMst2) - 135*lmMst1*(
        337 + lmMst2*(588 - 32*lmMt) - 132*lmMt - 660*S2 + 194*pow2(lmMst2) -
        48*pow2(lmMt)) - 6480*pow2(lmMt) - 405*lmMst2*(31 + 54*lmMt + 220*S2 +
        16*pow2(lmMt)) + 26190*pow3(lmMst2))*pow4(Mst1) + 27*pow2(Mst1)*(56439
        - 24000*lmMt - 1120*OepS2 - 25596*S2 - 360*lmMst2*(-12 + 44*lmMt + 63*
        S2) - 3840*(1 + lmMst2)*pow2(lmMst1) + 480*(163 - 16*lmMt)*pow2(lmMst2)
        - 120*lmMst1*(184 + 8*lmMst2*(57 - 8*lmMt) - 64*lmMt - 189*S2 + 164*
        pow2(lmMst2)) - 11520*pow2(lmMt) + 23520*pow3(lmMst2))*pow4(Mst2) + 20*
        (678923 + 19440*lmMt - 32480*OepS2 + 881712*S2 + 108*(-457 + 12*lmMst2
        + 126*lmMt)*pow2(lmMst1) + 540*(661 - 30*lmMt)*pow2(lmMst2) - 108*
        lmMst1*(1841 + lmMst2*(2626 - 24*lmMt) - 420*lmMt - 6090*S2 + 840*pow2(
        lmMst2) - 288*pow2(lmMt)) - 15552*pow2(lmMt) - 108*lmMst2*(619 + 498*
        lmMt + 6090*S2 + 288*pow2(lmMt)) + 216*pow3(lmMst1) + 89208*pow3(
        lmMst2))*pow6(Mst1) + 207360*pow2(1 + lmMst2)*pow6(Mst2))))/(21870.*
        Mgl*pow2(Msq)*pow2(Mst1)*pow5(Mst2)) + Mt*(-(pow2(Mst2)*(-2*s2t*(5*
        Dmsqst2*(3 - 2*lmMst2)*shiftst1*pow2(Mst2)*pow2(s2t) + (1 - 2*lmMst2)*
        pow2(Msq)*(-2*shiftst3*pow2(Mt) + 5*shiftst1*pow2(Mst2)*pow2(s2t))) + (
        1 - 2*lmMst2)*shiftst3*pow2(Msq)*((1 + lmMst1 - lmMst2)*pow2(Mst1) -
        pow2(Mst2))*pow3(s2t)) + 10*s2t*shiftst1*(Dmsqst2*(3 - 2*lmMst2) + (1 -
        2*lmMst2)*pow2(Msq))*(4*pow2(Mst2)*pow2(Mt) - 2*pow2(Mst1)*(2*pow2(Mt)
        + (-lmMst1 + lmMst2)*pow2(Mst2)*pow2(s2t)) + pow2(s2t)*pow4(Mst1)))/(3.
        *pow2(Msq)*pow2(Mst1)) + pow3(s2t)*((pow2(Mst1)*(202500*Dmglst2*
        Dmsqst2*(55 - 34*lmMst1 + 34*lmMst2) + 60*Dmsqst2*Mgl*(97294 - 17235*
        lmMst2 + 45*lmMst1*(233 + 330*lmMst2) - 7425*pow2(lmMst1) - 7425*pow2(
        lmMst2)) + 5*Mgl*pow2(Msq)*(4627751 + 48600*B4 + 5400*D3 - 5400*DN +
        1320240*lmMst1 - 662400*pow2(lmMst1) - 30*lmMst2*(12148 - 76560*lmMst1
        + 12105*pow2(lmMst1)) + 2250*(-583 + 438*lmMst1)*pow2(lmMst2) - 49500*
        pow3(lmMst1) - 572850*pow3(lmMst2)) - Dmglst2*pow2(Msq)*(5430043 +
        524580*lmMst2 + 900*(-859 + 3690*lmMst2)*pow2(lmMst1) + 1454400*pow2(
        lmMst2) - 60*lmMst1*(51493 + 24630*lmMst2 + 112950*pow2(lmMst2)) +
        45000*pow3(lmMst1) + 3411000*pow3(lmMst2))))/(60750.*Mgl*pow2(Msq)) + (
        (193.90296838394383 - (61388401*lmMst1)/396900. + lmMst2*(
        147.3919148400101 + (302047*lmMst1)/945. - (514*pow2(lmMst1))/9.) - (
        152966*pow2(lmMst1))/945. - (157.75767195767196 - 122*lmMst1)*pow2(
        lmMst2) + (Dmsqst2*(35.54686742892336 - (905536*lmMst1)/19845. + (
        45.630435878054925 + (244*lmMst1)/7.)*lmMst2 + (10*Dmglst2*(419 - 57*
        lmMst1 + 57*lmMst2))/(27.*Mgl) - (122*pow2(lmMst1))/7. - (122*pow2(
        lmMst2))/7.))/pow2(Msq) - (70*pow3(lmMst1))/27. - (1682*pow3(lmMst2))/
        27. - (Dmglst2*(407.7379592503276 + (2171669*lmMst2)/59535. + (2*(433 +
        6552*lmMst2)*pow2(lmMst1))/189. + (3974*pow2(lmMst2))/189. - (lmMst1*(
        2740559 + 1524600*lmMst2 + 7514640*pow2(lmMst2)))/59535. - (112*pow3(
        lmMst1))/27. + (1648*pow3(lmMst2))/27.))/Mgl)*pow4(Mst1))/pow2(Mst2) +
        ((-(pow2(Mst2)*(Dmsqst2*(4320*Dmglst2*(5 + 2*lmMst1 - 2*lmMst2) + 30*
        Mgl*(-1177 + 336*lmMst2 + 48*lmMst1*(-7 + 3*lmMst2) - 144*pow2(lmMst2))
        ) + pow2(Msq)*(-(Mgl*(21005 + 1440*B4 - 144*D3 + 72*DN + 21216*lmMst1 +
        1908*pow2(lmMst1) - 12*lmMst2*(2950 - 2040*lmMst1 + 117*pow2(lmMst1)) +
        108*(-283 + 98*lmMst1)*pow2(lmMst2) + 576*pow3(lmMst1) - 9756*pow3(
        lmMst2))) + 2*Dmglst2*(14987 - 432*B4 + 576*D3 - 360*DN + 4752*lmMst1 -
        1404*pow2(lmMst1) + 144*lmMst2*(-81 - 124*lmMst1 + 8*pow2(lmMst1)) -
        36*(-439 + 246*lmMst1)*pow2(lmMst2) + 7704*pow3(lmMst2)))))/324. - ((
        360*Dmglst2*Dmsqst2 + Dmsqst2*(30 - 60*lmMst2)*Mgl + Mgl*(119 + 218*
        lmMst2 + 32*lmMst1*(1 + lmMst2) + 107*pow2(lmMst2))*pow2(Msq) + 2*
        Dmglst2*(-20 + (230 + 32*lmMst1)*lmMst2 + 155*pow2(lmMst2))*pow2(Msq))*
        pow4(Mst2))/(9.*pow2(Mst1)))/pow2(Msq) + (-30*Dmsqst2*Mgl*(18*(40*OepS2
        - 27*(59 + 30*lmMst1 - 30*lmMst2)*S2)*pow2(Mst1)*pow2(Mst2) + 4*(208*
        OepS2 - 27*(347 + 156*lmMst1 - 156*lmMst2)*S2)*pow4(Mst1) + 27*(8*OepS2
        - 81*(-15 + 2*lmMst1 - 2*lmMst2)*S2)*pow4(Mst2)) + pow2(Msq)*(3*Mgl*(-
        60*(340*OepS2 - 81*(424 + 85*lmMst1 - 85*lmMst2)*S2)*pow2(Mst1)*pow2(
        Mst2) - 35*(788*OepS2 - 27*(2662 + 591*lmMst1 - 591*lmMst2)*S2)*pow4(
        Mst1) + 27*(-184*OepS2 + 81*(307 + 46*lmMst1 - 46*lmMst2)*S2)*pow4(
        Mst2)) + 2*Dmglst2*(27*(632*OepS2 + 9*(16193 - 1422*lmMst1 + 1422*
        lmMst2)*S2)*pow2(Mst1)*pow2(Mst2) + 2*(25328*OepS2 + 27*(47051 - 18996*
        lmMst1 + 18996*lmMst2)*S2)*pow4(Mst1) + 27*(56*OepS2 - 81*(-1677 + 14*
        lmMst1 - 14*lmMst2)*S2)*pow4(Mst2))))/(2187.*pow2(Msq)*pow2(Mst2)) + (
        16*(1 + lmMst2)*(4*Dmglst2*lmMst2 + Mgl + lmMst2*Mgl)*pow6(Mst2))/(9.*
        pow4(Mst1)))/Mgl) - (2*T1ep*(-3*Mgl*pow2(Msq)*(3*pow2(Mst1)*pow2(Mst2)*
        (-3470*Mst2*s2t*pow2(Mt) - 627*Mt*pow2(Mst2)*pow2(s2t) + 1760*pow3(Mt)
        + 1700*pow3(Mst2)*pow3(s2t)) + (-24124*Mst2*s2t*pow2(Mt) - 3198*Mt*
        pow2(Mst2)*pow2(s2t) + 16240*pow3(Mt) + 6895*pow3(Mst2)*pow3(s2t))*
        pow4(Mst1) + 27*(-106*Mst2*s2t*pow2(Mt) - 21*Mt*pow2(Mst2)*pow2(s2t) +
        28*pow3(Mt) + 46*pow3(Mst2)*pow3(s2t))*pow4(Mst2)) + 60*Dmsqst2*Mgl*
        Mst2*s2t*(2*(53*pow2(Mt) - 52*pow2(Mst2)*pow2(s2t))*pow4(Mst1) + 54*
        pow2(Mt)*pow4(Mst2) + 90*pow2(Mst1)*(pow2(Mst2)*pow2(Mt) - pow2(s2t)*
        pow4(Mst2)) - 27*pow2(s2t)*pow6(Mst2)) + Dmglst2*pow2(Msq)*(27*pow2(
        Mst1)*pow2(Mst2)*(-800*Mst2*s2t*pow2(Mt) - 907*Mt*pow2(Mst2)*pow2(s2t)
        + 1984*pow3(Mt) + 316*pow3(Mst2)*pow3(s2t)) + 2*(-66168*Mst2*s2t*pow2(
        Mt) - 35601*Mt*pow2(Mst2)*pow2(s2t) + 129704*pow3(Mt) + 12664*pow3(
        Mst2)*pow3(s2t))*pow4(Mst1) + 189*(20*pow3(Mt)*pow4(Mst2) - 15*Mt*pow2(
        s2t)*pow6(Mst2) + 4*pow3(s2t)*pow7(Mst2)))))/(729.*Mgl*pow2(Msq)*pow5(
        Mst2)))))/Tbeta + ((5*shiftst1*(Dmsqst2*(3 - 2*lmMst2) + (1 - 2*lmMst2)
        *pow2(Msq))*(-(pow4(Mst2)*(pow2(s2t)*(-8*pow2(Mt) + (pow2(Mst1) - pow2(
        Mst2))*pow2(s2t))*pow2(Sbeta)*pow4(Mst1) + pow2(Mst1)*(4*pow2(Mt)*pow2(
        s2t)*(-(pow2(MuSUSY)*(-1 + pow2(Sbeta))) + 4*pow2(Mst2)*pow2(Sbeta)) +
        pow2(Sbeta)*(16*pow4(Mt) - pow4(Mst2)*pow4(s2t))) + pow2(Mst2)*(-4*
        pow2(Mt)*pow2(s2t)*(pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 2*pow2(Mst2)*
        pow2(Sbeta)) + pow2(Sbeta)*(16*pow4(Mt) + pow4(Mst2)*pow4(s2t))))) + 2*
        (lmMst1 - lmMst2)*pow2(Mst1)*pow2(s2t)*(-4*pow2(Mt)*pow2(MuSUSY)*(-1 +
        pow2(Sbeta))*(pow4(Mst1) + pow4(Mst2)) - pow2(Mst1)*pow2(Mst2)*(4*pow2(
        Mt)*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + pow2(s2t)*pow2(Sbeta)*pow4(Mst2))
        + pow2(s2t)*pow2(Sbeta)*pow8(Mst2))))/(6.*pow2(Msq)*pow2(Sbeta)) + ((1
        - 2*lmMst2)*shiftst3*((8*(-1 + lmMst1 - lmMst2)*pow2(Mst1)*pow2(Mt)*
        pow2(s2t) - 16*pow4(Mt) - (-4*pow2(Mst1)*pow2(Mst2) + (3 + 2*lmMst1 -
        2*lmMst2)*pow4(Mst1) + pow4(Mst2))*pow4(s2t))*pow6(Mst2) - (4*pow2(Mt)*
        pow2(MuSUSY)*pow2(s2t)*(2*(1 - 2*lmMst1 + 2*lmMst2)*pow2(Mst2)*pow4(
        Mst1) + 2*(1 - lmMst1 + lmMst2)*pow2(Mst1)*pow4(Mst2) + (2 - 6*lmMst1 +
        6*lmMst2)*pow6(Mst1) + pow6(Mst2)))/pow2(Sbeta) + 4*pow2(Mt)*pow2(s2t)*
        (pow2(MuSUSY)*pow6(Mst2) + 2*(pow2(MuSUSY)*((1 - 2*lmMst1 + 2*lmMst2)*
        pow2(Mst2)*pow4(Mst1) + (1 - lmMst1 + lmMst2)*pow2(Mst1)*pow4(Mst2) + (
        1 - 3*lmMst1 + 3*lmMst2)*pow6(Mst1)) + pow8(Mst2)))))/12.)/(pow2(Mst1)*
        pow4(Mst2))) - Al4p*xDR2DRMOD*(((2*Dmglst2*lmMst2 + Mgl + lmMst2*Mgl)*
        twoLoopFlag*(-(pow2(Mst2)*pow2(s2t)*pow4(Mst1)*(-8*(lmMst1 - lmMst2)*
        Tbeta*pow2(Mt)*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 4*Mt*(-(MuSUSY*s2t) +
        2*Mt*Tbeta)*pow2(Mst2)*pow2(Sbeta) + (1 - 2*lmMst1 + 2*lmMst2)*Tbeta*
        pow2(s2t)*pow2(Sbeta)*pow4(Mst2))) - pow2(Mst1)*pow4(Mst2)*(-4*Tbeta*
        pow2(Mt)*pow2(s2t)*((-1 + 2*lmMst1 - 2*lmMst2)*pow2(MuSUSY)*(-1 + pow2(
        Sbeta)) + 4*pow2(Mst2)*pow2(Sbeta)) + pow2(Sbeta)*(16*MuSUSY*s2t*pow3(
        Mt) + 8*(-lmMst1 + lmMst2)*Mt*MuSUSY*pow2(Mst2)*pow3(s2t) - 16*Tbeta*
        pow4(Mt) + (1 + 2*lmMst1 - 2*lmMst2)*Tbeta*pow4(Mst2)*pow4(s2t))) +
        Tbeta*pow2(s2t)*(8*(lmMst1 - lmMst2)*pow2(Mt)*pow2(MuSUSY)*(-1 + pow2(
        Sbeta)) + pow2(s2t)*pow2(Sbeta)*pow4(Mst2))*pow6(Mst1) + (-4*Tbeta*
        pow2(Mt)*pow2(s2t)*(pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 2*pow2(Mst2)*
        pow2(Sbeta)) + pow2(Sbeta)*(16*MuSUSY*s2t*pow3(Mt) - 4*Mt*MuSUSY*pow2(
        Mst2)*pow3(s2t) + 16*Tbeta*pow4(Mt) + Tbeta*pow4(Mst2)*pow4(s2t)))*
        pow6(Mst2)))/(6.*Mgl*Tbeta*pow2(Mst1)*pow2(Sbeta)*pow4(Mst2)) + (Al4p*
        threeLoopFlag*xDmsqst2*pow2(Dmsqst2)*((20*pow2(Mt)*pow2(MuSUSY)*pow2(
        s2t)*(11 - 22*lmMst1 + 22*lmMst2 + (10*Dmglst2*(1 - 2*lmMst1 + 2*
        lmMst2))/Mgl + ((10*Dmglst2 + 11*Mgl)*pow2(Mst2))/(Mgl*pow2(Mst1)) + (
        lmMst1 - lmMst2)*((2*(4*Dmglst2 - 11*Mgl)*pow2(Mst1))/(Mgl*pow2(Mst2))
        - (20*pow4(Mst1))/pow4(Mst2))))/pow2(Sbeta) + (20*Mt*MuSUSY*s2t*(10*
        Dmglst2 + Mgl*(11 - 18*shiftst1))*(-4*pow2(Mst1)*pow2(Mt) + 4*pow2(
        Mst2)*pow2(Mt) + 2*(lmMst1 - lmMst2)*pow2(Mst1)*pow2(Mst2)*pow2(s2t) +
        pow2(s2t)*(pow4(Mst1) - pow4(Mst2))))/(Mgl*Tbeta*pow2(Mst1)) + 80*(11 +
        (10*Dmglst2)/Mgl)*(1 + pow2(Mst2)/pow2(Mst1))*pow4(Mt) - 90*shiftst1*((
        4*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*(1 + pow2(Mst2)/pow2(Mst1) - (2*(
        lmMst1 - lmMst2)*(pow2(Mst1)*pow2(Mst2) + pow4(Mst1) + pow4(Mst2)))/
        pow4(Mst2)))/pow2(Sbeta) - 4*pow2(Mt)*pow2(s2t)*(2*pow2(Mst1) - 4*pow2(
        Mst2) + (2*pow4(Mst2))/pow2(Mst1) + pow2(MuSUSY)*(1 + pow2(Mst2)/pow2(
        Mst1) - (2*(lmMst1 - lmMst2)*(pow2(Mst1)*pow2(Mst2) + pow4(Mst1) +
        pow4(Mst2)))/pow4(Mst2))) + 16*(1 + pow2(Mst2)/pow2(Mst1))*pow4(Mt) + (
        (pow2(Mst1) - pow2(Mst2))*(2*(lmMst1 - lmMst2)*pow2(Mst1)*pow2(Mst2) +
        pow4(Mst1) - pow4(Mst2))*pow4(s2t))/pow2(Mst1)) + (24*xDmglst2*pow2(
        Dmglst2)*((5*(1 - 2*lmMst1 + 2*lmMst2 + (20*(lmMst1 - lmMst2)*pow2(
        Mst1))/pow2(Mst2) + pow2(Mst2)/pow2(Mst1))*pow2(Mt)*pow2(MuSUSY)*pow2(
        s2t))/(2.*pow2(Sbeta)) - (5*MuSUSY*(4*s2t*(pow2(Mst1) - pow2(Mst2))*
        pow3(Mt) - Mt*pow3(s2t)*(2*(lmMst1 - lmMst2)*pow2(Mst1)*pow2(Mst2) +
        pow4(Mst1) - pow4(Mst2))))/(2.*Tbeta*pow2(Mst1)) - pow2(Mt)*pow2(s2t)*(
        5*pow2(Mst1) - 10*pow2(Mst2) + (5*(1 - 2*lmMst1 + 2*lmMst2 + (20*(
        lmMst1 - lmMst2)*pow2(Mst1))/pow2(Mst2) + pow2(Mst2)/pow2(Mst1))*pow2(
        MuSUSY))/2. + (5*pow4(Mst2))/pow2(Mst1)) + 10*(1 + pow2(Mst2)/pow2(
        Mst1))*pow4(Mt) + (5*pow4(s2t)*(3*(-1 + 2*lmMst1 - 2*lmMst2)*pow2(Mst2)
        *pow4(Mst1) - 3*(1 + 2*lmMst1 - 2*lmMst2)*pow2(Mst1)*pow4(Mst2) + 10*
        pow6(Mst1) + 3*pow6(Mst2)))/(24.*pow2(Mst1))))/pow2(Mgl) + (5*pow2(s2t)
        *((10*Dmglst2 + 11*Mgl)*(pow2(Mst1) - pow2(Mst2))*pow2(s2t)*(2*(lmMst1
        - lmMst2)*pow2(Mst1)*pow2(Mst2) + pow4(Mst1) - pow4(Mst2)) - (4*pow2(
        Mt)*(2*Dmglst2*(2*pow2(Mst2)*(5*pow2(Mst2) + 2*(lmMst1 - lmMst2)*pow2(
        MuSUSY))*pow4(Mst1) - 5*pow2(Mst1)*(4*pow2(Mst2) + (-1 + 2*lmMst1 - 2*
        lmMst2)*pow2(MuSUSY))*pow4(Mst2) + 5*(2*pow2(Mst2) + pow2(MuSUSY))*
        pow6(Mst2)) + Mgl*(22*pow2(Mst2)*(pow2(Mst2) + (-lmMst1 + lmMst2)*pow2(
        MuSUSY))*pow4(Mst1) - 11*pow2(Mst1)*(4*pow2(Mst2) + (-1 + 2*lmMst1 - 2*
        lmMst2)*pow2(MuSUSY))*pow4(Mst2) - 20*(lmMst1 - lmMst2)*pow2(MuSUSY)*
        pow6(Mst1) + 11*(2*pow2(Mst2) + pow2(MuSUSY))*pow6(Mst2))))/pow4(Mst2))
        )/(Mgl*pow2(Mst1))))/(24.*pow4(Msq)) + (xDmglst2*pow2(Dmglst2)*(-((2 -
        3*lmMst2)*twoLoopFlag*(4*Mt*MuSUSY*s2t*((-4*pow2(Mst1)*pow2(Mt) + 4*
        pow2(Mst2)*pow2(Mt) + 2*(lmMst1 - lmMst2)*pow2(Mst1)*pow2(Mst2)*pow2(
        s2t) + pow2(s2t)*(pow4(Mst1) - pow4(Mst2)))/Tbeta + (Mt*MuSUSY*s2t*((
        pow2(Mst1) + pow2(Mst2))*pow4(Mst2) - 2*(lmMst1 - lmMst2)*pow2(Mst1)*(
        pow2(Mst1)*pow2(Mst2) + pow4(Mst1) + pow4(Mst2))))/(pow2(Sbeta)*pow4(
        Mst2))) + 16*(pow2(Mst1) + pow2(Mst2))*pow4(Mt) + (pow2(Mst1) - pow2(
        Mst2))*(2*(lmMst1 - lmMst2)*pow2(Mst1)*pow2(Mst2) + pow4(Mst1) - pow4(
        Mst2))*pow4(s2t) - (4*pow2(Mt)*pow2(s2t)*(2*pow4(Mst1)*((-lmMst1 +
        lmMst2)*pow2(Mst2)*pow2(MuSUSY) + pow4(Mst2)) + 2*(-lmMst1 + lmMst2)*
        pow2(MuSUSY)*pow6(Mst1) + pow2(Mst1)*((1 - 2*lmMst1 + 2*lmMst2)*pow2(
        MuSUSY)*pow4(Mst2) - 4*pow6(Mst2)) + (2*pow2(Mst2) + pow2(MuSUSY))*
        pow6(Mst2)))/pow4(Mst2)))/(6.*pow2(Mst1)) + (Al4p*threeLoopFlag*((6400*
        s2t*pow3(Mt)*(-3*pow4(Mst1)*(pow2(Mst2)*pow2(MuSUSY)*(-11 - 371*lmMst2
        + 42*pow2(lmMst2) - 6*lmMst1*(-21 - 5*lmMst2 + 24*pow2(lmMst2)) + 144*
        pow3(lmMst2)) + 6*(93 + 88*lmMst2 - 22*lmMt - 8*lmMst2*lmMt + lmMst1*(-
        34 - 8*lmMst2 + 4*lmMt) + 12*pow2(lmMst2))*pow4(Mst2)) + 9*((49 + 103*
        lmMst2 - 36*(1 + lmMst2)*pow2(lmMst2) + 4*lmMst1*(-11 + 6*lmMst2 + 9*
        pow2(lmMst2)))*pow2(MuSUSY) + pow2(Mst2)*(28 + lmMst2*(378 - 96*lmMt) -
        44*lmMt + lmMst1*(-274 + 60*lmMt + 18*lmMst2*(-13 + 2*lmMt) - 36*pow2(
        lmMst2)) + (270 - 36*lmMt)*pow2(lmMst2) + 36*pow3(lmMst2)))*pow6(Mst1)
        + 36*(2 + lmMst2 - 3*pow2(lmMst2))*(2*pow2(Mst2) + pow2(MuSUSY))*pow6(
        Mst2) - pow2(Mst1)*(3*pow2(MuSUSY)*(-53 - 191*lmMst2 + lmMst1*(60 + 18*
        lmMst2 - 72*pow2(lmMst2)) + 54*pow2(lmMst2) + 72*pow3(lmMst2))*pow4(
        Mst2) + (300 + 72*lmMt + lmMst2*(-341 + 78*lmMt) + 18*lmMst1*(8 - 4*
        lmMt + lmMst2*(5 + 2*lmMt) - 2*pow2(lmMst2)) - 12*(26 + 3*lmMt)*pow2(
        lmMst2) + 36*pow3(lmMst2))*pow6(Mst2))))/(pow2(Mst1)*pow5(Mst2)) + (75*
        pow3(s2t)*((128*Mt*(-((-71 + 6*lmMst1*(-1 + lmMst2) + lmMst2 + 30*pow2(
        lmMst2))*pow2(Mst2)*pow4(Mst1)) + pow2(Mst1)*(-5 - 167*lmMst2 + lmMst1*
        (60 + 18*lmMst2 - 72*pow2(lmMst2)) - 18*pow2(lmMst2) + 72*pow3(lmMst2))
        *pow4(Mst2) - 3*(1 + 13*lmMst2 - 2*lmMst1*(5 + 3*lmMst2) + 6*pow2(
        lmMst2))*pow6(Mst1) + 12*(-2 - lmMst2 + 3*pow2(lmMst2))*pow6(Mst2)))/
        Mst2 + s2t*(pow2(Mst2)*(-340 + 806*lmMst2 + 672*(-2 + 3*lmMst2)*pow2(
        lmMst1) - 3159*pow2(lmMst2) - 6*lmMst1*(388 - 934*lmMst2 + 399*pow2(
        lmMst2)) + 378*pow3(lmMst2))*pow4(Mst1) + pow2(Mst1)*(1532 - 3514*
        lmMst2 + 96*(-2 + 3*lmMst2)*pow2(lmMst1) - 207*pow2(lmMst2) - 6*lmMst1*
        (-532 + 550*lmMst2 + 369*pow2(lmMst2)) + 1926*pow3(lmMst2))*pow4(Mst2)
        + 6*(-36 + (70 + 32*lmMst1)*lmMst2 + 27*pow2(lmMst2))*pow6(Mst1) + 3*(-
        548 + 422*lmMst2 + 32*lmMst1*(-2 + 3*lmMst2) + 817*pow2(lmMst2))*pow6(
        Mst2) - (45*Dmsqst2*(15*(-1 + 2*lmMst1 - 2*lmMst2)*pow2(Mst2)*pow4(
        Mst1) - 15*(1 + 2*lmMst1 - 2*lmMst2)*pow2(Mst1)*pow4(Mst2) + 2*pow6(
        Mst1) + 15*pow6(Mst2)))/pow2(Msq))))/pow2(Mst1) + (16*pow4(Mt)*(225*(-
        356 - 26*lmMst2 + 32*lmMst1*(-2 + 3*lmMst2) + 177*pow2(lmMst2))*pow2(
        Mst1)*pow2(Mst2) - (50625*Dmsqst2*pow2(Mst1)*(pow2(Mst1) + pow2(Mst2)))
        /pow2(Msq) + (183748 - 240*(137 + 60*lmMst1)*lmMt + 6*lmMst2*(-587 +
        4800*lmMst1 + 6120*lmMt) - 25695*pow2(lmMst2))*pow4(Mst1) - 7200*(-2 +
        lmMst2 + 5*pow2(lmMst2))*pow4(Mst2) + (7200*(221 + lmMst2*(299 - 47*
        lmMt) - 62*lmMt + (87 - 6*lmMt)*pow2(lmMst2) - 2*lmMst1*(65 + lmMst2*(
        32 - 3*lmMt) - 12*lmMt + 3*pow2(lmMst2)) + 6*pow3(lmMst2))*pow6(Mst1))/
        pow2(Mst2) + (14400*(-42 + 33*lmMt + lmMst2*(-242 + 62*lmMt) + 6*(-29 +
        4*lmMt)*pow2(lmMst2) + 2*lmMst1*(81 + lmMst2*(74 - 12*lmMt) - 18*lmMt +
        12*pow2(lmMst2)) - 24*pow3(lmMst2))*pow8(Mst1))/pow4(Mst2)))/pow4(Mst1)
        - (100*Mt*MuSUSY*((128*pow3(Mt)*(pow2(Mst2)*(-834 + 162*lmMt + lmMst2*(
        -569 + 33*lmMt) + 3*(-17 + 6*lmMt)*pow2(lmMst2) + 18*lmMst1*(13 +
        lmMst2 - lmMst2*lmMt + pow2(lmMst2)) - 18*pow3(lmMst2))*pow4(Mst1) +
        pow2(Mst1)*(lmMst2*(187 - 39*lmMt) - 12*(8 + 3*lmMt) + 3*(19 + 6*lmMt)*
        pow2(lmMst2) + 18*lmMst1*(-3 + 2*lmMt - lmMst2*(3 + lmMt) + pow2(
        lmMst2)) - 18*pow3(lmMst2))*pow4(Mst2) + (290 + lmMst2*(2447 - 645*
        lmMt) - 375*lmMt - 3*(-509 + 60*lmMt)*pow2(lmMst2) - 18*lmMst1*(87 +
        lmMst2*(71 - 10*lmMt) - 22*lmMt + 10*pow2(lmMst2)) + 180*pow3(lmMst2))*
        pow6(Mst1) + 36*(2 + lmMst2 - 3*pow2(lmMst2))*pow6(Mst2)))/(pow2(Mst1)*
        pow5(Mst2)) + (288*Mt*pow2(s2t)*(6*pow2(Mst2)*(7 - 30*lmMst2 + lmMst1*(
        11 + 2*lmMst2 - 12*pow2(lmMst2)) - 2*pow2(lmMst2) + 12*pow3(lmMst2))*
        pow4(Mst1) + pow2(Mst1)*(-29 - 179*lmMst2 + lmMst1*(60 + 18*lmMst2 -
        72*pow2(lmMst2)) + 18*pow2(lmMst2) + 72*pow3(lmMst2))*pow4(Mst2) - 3*(
        17 + 51*lmMst2 - 4*pow2(lmMst2) + 4*lmMst1*(-6 + lmMst2 + 3*pow2(
        lmMst2)) - 12*pow3(lmMst2))*pow6(Mst1) + 12*(-2 - lmMst2 + 3*pow2(
        lmMst2))*pow6(Mst2)))/(pow2(Mst1)*pow3(Mst2)) - 3*pow3(s2t)*(pow2(Mst1)
        *(260 + 1538*lmMst2 + 768*(2 - 3*lmMst2)*pow2(lmMst1) + 1395*pow2(
        lmMst2) + 96*lmMst1*(-7 - 27*lmMst2 + 48*pow2(lmMst2)) - (675*Dmsqst2)/
        pow2(Msq) - 2304*pow3(lmMst2)) - 2*pow2(Mst2)*(40 - 1172*lmMst2 + 48*(-
        2 + 3*lmMst2)*pow2(lmMst1) + 882*pow2(lmMst2) - 3*lmMst1*(-500 + 502*
        lmMst2 + 369*pow2(lmMst2)) + (675*Dmsqst2*(lmMst1 - lmMst2))/pow2(Msq)
        + 963*pow3(lmMst2)) - (96*(16*lmMst2*pow2(lmMst1) - 2*lmMst1*(3 + 2*
        lmMst2 + 16*pow2(lmMst2)) + lmMst2*(7 + 4*lmMst2 + 16*pow2(lmMst2)))*
        pow4(Mst1))/pow2(Mst2) + (27*(53.77777777777778 + (64*lmMst1)/9. - (2*(
        65 + 16*lmMst1)*lmMst2)/3. - 73*pow2(lmMst2) + (25*Dmsqst2)/pow2(Msq))*
        pow4(Mst2))/pow2(Mst1) + (96*(-2 + lmMst2 + 5*pow2(lmMst2))*pow6(Mst2))
        /pow4(Mst1)) - (12*s2t*pow2(Mt)*((675*Dmsqst2*pow2(Mst1)*(pow2(Mst1) -
        pow2(Mst2)))/pow2(Msq) + 3*(-324 + 134*lmMst2 + 32*lmMst1*(-2 + 3*
        lmMst2) + 177*pow2(lmMst2))*pow2(Mst1)*pow2(Mst2) + (860 - 434*lmMst2 +
        96*(-2 + 3*lmMst2)*pow2(lmMst1) + 96*lmMst1*(3 + lmMst2 - 6*pow2(
        lmMst2)) - 243*pow2(lmMst2) + 288*pow3(lmMst2))*pow4(Mst1) - 96*(-2 +
        lmMst2 + 5*pow2(lmMst2))*pow4(Mst2) + (96*((-2 + 3*lmMst2)*pow2(lmMst1)
        + lmMst1*(-8 + 2*lmMst2 - 6*pow2(lmMst2)) + 3*(6 + 5*lmMst2 + pow3(
        lmMst2)))*pow6(Mst1))/pow2(Mst2) + (192*(-6 - 14*lmMst2 + lmMst2*pow2(
        lmMst1) + lmMst1*(9 + 6*lmMst2 - 2*pow2(lmMst2)) - 6*pow2(lmMst2) +
        pow3(lmMst2))*pow8(Mst1))/pow4(Mst2)))/pow4(Mst1)))/Tbeta - (300*pow2(
        Mt)*pow2(s2t)*(-675*Dmsqst2*pow2(Mst1)*pow2(Mst2)*(2*(pow2(Mst2) + (-
        lmMst1 + lmMst2)*pow2(MuSUSY))*pow4(Mst1) + pow2(Mst1)*((1 - 2*lmMst1 +
        2*lmMst2)*pow2(Mst2)*pow2(MuSUSY) - 4*pow4(Mst2)) + (2*pow2(Mst2) +
        pow2(MuSUSY))*pow4(Mst2)) + pow2(Msq)*((4*pow2(Mst2)*(916 - 418*lmMst2
        + 48*(-2 + 3*lmMst2)*pow2(lmMst1) - 387*pow2(lmMst2) - 48*lmMst1*(-5 +
        2*lmMst2 + 6*pow2(lmMst2)) + 144*pow3(lmMst2)) + pow2(MuSUSY)*(-1180 -
        1270*lmMst2 + 96*(-2 + 3*lmMst2)*pow2(lmMst1) + 3255*pow2(lmMst2) - 6*
        lmMst1*(-468 + 454*lmMst2 + 369*pow2(lmMst2)) + 1926*pow3(lmMst2)))*
        pow4(Mst1)*pow4(Mst2) + 2*pow2(Mst2)*((868 + 96*lmMst1*(-11 + lmMst2) +
        1874*lmMst2 + 243*pow2(lmMst2))*pow2(Mst2) + 3*pow2(MuSUSY)*(-240 -
        468*lmMst2 + 144*(-2 + 3*lmMst2)*pow2(lmMst1) + lmMst1*(580 - 22*lmMst2
        - 1137*pow2(lmMst2)) + 310*pow2(lmMst2) + 705*pow3(lmMst2)))*pow6(Mst1)
        - 3*pow2(Mst1)*(-2*(-388 + 166*lmMst2 + 32*lmMst1*(-2 + 3*lmMst2) +
        337*pow2(lmMst2))*pow2(Mst2) + (420 + lmMst1*(64 - 96*lmMst2) - 358*
        lmMst2 - 497*pow2(lmMst2))*pow2(MuSUSY))*pow6(Mst2) - 12*(32*(2 + 7*
        lmMst2 - lmMst1*(5 + 4*lmMst2) + 4*pow2(lmMst2))*pow2(Mst2) - pow2(
        MuSUSY)*(48 + 4*lmMst2*(45 + 68*pow2(lmMst1)) + 310*pow2(lmMst2) -
        lmMst1*(92 + 310*lmMst2 + 635*pow2(lmMst2)) + 363*pow3(lmMst2)))*pow8(
        Mst1) - 96*(-2 + lmMst2 + 5*pow2(lmMst2))*(2*pow2(Mst2) + pow2(MuSUSY))
        *pow8(Mst2))))/(pow2(Msq)*pow4(Mst1)*pow4(Mst2)) + (300*s2t*pow2(Mt)*
        pow2(MuSUSY)*(675*Dmsqst2*s2t*pow2(Mst1)*(2*(lmMst1 - lmMst2)*pow2(
        Mst1) - pow2(Mst2))*(pow2(Mst1) + pow2(Mst2))*pow3(Mst2) + pow2(Msq)*(-
        ((Mst2*s2t*(1180 + 1270*lmMst2 + 96*(2 - 3*lmMst2)*pow2(lmMst1) - 3255*
        pow2(lmMst2) + 6*lmMst1*(-468 + 454*lmMst2 + 369*pow2(lmMst2)) - 1926*
        pow3(lmMst2)) + 64*Mt*(53 + 191*lmMst2 - 54*pow2(lmMst2) + 6*lmMst1*(-
        10 - 3*lmMst2 + 12*pow2(lmMst2)) - 72*pow3(lmMst2)))*pow4(Mst1)*pow4(
        Mst2)) - 2*pow2(Mst2)*(3*Mst2*s2t*(240 + 468*lmMst2 + 144*(2 - 3*
        lmMst2)*pow2(lmMst1) - 310*pow2(lmMst2) + lmMst1*(-580 + 22*lmMst2 +
        1137*pow2(lmMst2)) - 705*pow3(lmMst2)) + 32*Mt*(11 + 371*lmMst2 - 42*
        pow2(lmMst2) + 6*lmMst1*(-21 - 5*lmMst2 + 24*pow2(lmMst2)) - 144*pow3(
        lmMst2)))*pow6(Mst1) - 3*(Mst2*s2t*(420 + lmMst1*(64 - 96*lmMst2) -
        358*lmMst2 - 497*pow2(lmMst2)) + 256*Mt*(2 + lmMst2 - 3*pow2(lmMst2)))*
        pow2(Mst1)*pow6(Mst2) - 12*(16*Mt*(49 + 103*lmMst2 - 36*(1 + lmMst2)*
        pow2(lmMst2) + 4*lmMst1*(-11 + 6*lmMst2 + 9*pow2(lmMst2))) - Mst2*s2t*(
        48 + 4*lmMst2*(45 + 68*pow2(lmMst1)) + 310*pow2(lmMst2) - lmMst1*(92 +
        310*lmMst2 + 635*pow2(lmMst2)) + 363*pow3(lmMst2)))*pow8(Mst1) - 96*
        s2t*(-2 + lmMst2 + 5*pow2(lmMst2))*pow9(Mst2))))/(pow2(Msq)*pow2(Sbeta)
        *pow4(Mst1)*pow5(Mst2))))/8100.))/pow2(Mgl)) - ((threeLoopFlag*pow2(
        Al4p)*(12*pow2(z2)*(-4*Mst2*xDmglst2*pow2(Dmglst2)*pow4(Msq)*(3*Mst2*
        pow2(Mst1)*(4*pow2(Mst2)*pow2(Mt)*(-1157*Tbeta*pow2(MuSUSY)*pow2(s2t)*(
        -1 + pow2(Sbeta)) + 2*Mt*(1980*MuSUSY*s2t + 4703*Mt*Tbeta)*pow2(Sbeta))
        + s2t*(84045*MuSUSY*s2t - 17228*Mt*Tbeta)*pow2(Mt)*pow2(Sbeta)*pow3(
        Mst2) + 4*Mst2*MuSUSY*(29549*MuSUSY*s2t*Tbeta*(-1 + pow2(Sbeta)) +
        4328*Mt*pow2(Sbeta))*pow3(Mt) - 2*Mt*(265*MuSUSY*s2t + 3582*Mt*Tbeta)*
        pow2(s2t)*pow2(Sbeta)*pow4(Mst2) + 840*Tbeta*pow2(MuSUSY)*(-1 + pow2(
        Sbeta))*pow4(Mt) + 4*(767*Mt - 223*Mst2*s2t)*Tbeta*pow2(Sbeta)*pow3(
        s2t)*pow5(Mst2)) + (9*pow4(Mst2)*(2*Mst2*pow2(Mt)*pow2(s2t)*(-1366*
        Tbeta*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 3*Mst2*(10361*MuSUSY - 84*Mst2*
        Tbeta)*pow2(Sbeta)) + 4*s2t*(10361*Tbeta*pow2(MuSUSY)*(-1 + pow2(Sbeta)
        ) - 14*Mst2*(-18*MuSUSY + Mst2*Tbeta)*pow2(Sbeta))*pow3(Mt) + pow2(
        Sbeta)*(-2*Mt*(1366*MuSUSY + 10361*Mst2*Tbeta)*pow3(Mst2)*pow3(s2t) +
        56*(MuSUSY + 26*Mst2*Tbeta)*pow4(Mt) + 683*Tbeta*pow4(s2t)*pow5(Mst2)))
        )/2. + pow4(Mst1)*(2*Mst2*pow2(Mt)*pow2(s2t)*(-21196*Tbeta*pow2(MuSUSY)
        *(-1 + pow2(Sbeta)) - 3*Mst2*(-35195*MuSUSY + 9228*Mst2*Tbeta)*pow2(
        Sbeta)) + 4*s2t*(86405*Tbeta*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 4*Mst2*(
        8271*MuSUSY + 12865*Mst2*Tbeta)*pow2(Sbeta))*pow3(Mt) + pow2(Sbeta)*(-(
        Mt*(25328*MuSUSY + 15571*Mst2*Tbeta)*pow3(Mst2)*pow3(s2t)) - 16*(16213*
        MuSUSY + 5582*Mst2*Tbeta)*pow4(Mt) + 4199*Tbeta*pow4(s2t)*pow5(Mst2))))
        + 60*xDmsqst2*pow2(Dmsqst2)*pow2(Mgl)*pow2(Mst2)*(4*MuSUSY*pow2(Sbeta)*
        (-3*Mt*pow2(Mst2)*pow3(s2t)*(19*pow2(Mst1)*(pow2(Mst1) + pow2(Mst2)) +
        6*pow4(Mst2)) + 2*s2t*pow3(Mt)*(21*pow2(Mst1)*pow2(Mst2) + 53*pow4(
        Mst1) + 18*pow4(Mst2))) + 3*Tbeta*pow2(Sbeta)*pow4(Mst2)*(13*pow2(Mst1)
        *pow2(Mst2) + 6*pow4(Mst2))*pow4(s2t) + 4*Tbeta*pow2(Mt)*pow2(s2t)*((-
        221*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 7*pow2(Mst2)*pow2(Sbeta))*pow4(
        Mst1) - 3*pow2(Mst1)*(25*pow2(Mst2)*pow2(MuSUSY)*(-1 + pow2(Sbeta)) +
        pow2(Sbeta)*pow4(Mst2)) - 18*(pow2(MuSUSY)*(-1 + pow2(Sbeta))*pow4(
        Mst2) + pow2(Sbeta)*pow6(Mst2)))) + Mgl*pow2(Msq)*(4*Mt*MuSUSY*(Mst2*
        pow2(Sbeta)*(4*pow2(Msq)*pow3(Mt)*(72*(186*Dmglst2 - 55*Mgl)*pow2(Mst1)
        *pow2(Mst2) + 4*(16213*Dmglst2 - 3045*Mgl)*pow4(Mst1) + 189*(5*Dmglst2
        - 3*Mgl)*pow4(Mst2)) - 3*Mt*pow2(Msq)*pow2(Mst2)*pow2(s2t)*(9*(6091*
        Dmglst2 + 1519*Mgl)*pow2(Mst1)*pow2(Mst2) + 2*(35195*Dmglst2 + 6177*
        Mgl)*pow4(Mst1) + 27*(1763*Dmglst2 + 555*Mgl)*pow4(Mst2)) + 6*Mst2*s2t*
        pow2(Mt)*(15*(60*Dmsqst2*Mgl + (-240*Dmglst2 + 347*Mgl)*pow2(Msq))*
        pow2(Mst1)*pow2(Mst2) + 2*(530*Dmsqst2*Mgl + (-11028*Dmglst2 + 6031*
        Mgl)*pow2(Msq))*pow4(Mst1) + 27*Mgl*(20*Dmsqst2 + 53*pow2(Msq))*pow4(
        Mst2)) - pow3(Mst2)*pow3(s2t)*(36*(150*Dmsqst2*Mgl + (-237*Dmglst2 +
        425*Mgl)*pow2(Msq))*pow2(Mst1)*pow2(Mst2) + (6240*Dmsqst2*Mgl + (-
        25328*Dmglst2 + 20685*Mgl)*pow2(Msq))*pow4(Mst1) + 54*(30*Dmsqst2*Mgl -
        (158*Dmglst2 + 3*Mgl)*pow2(Msq))*pow4(Mst2))) - Mt*MuSUSY*Tbeta*(-(
        pow2(Mst2)*pow2(s2t)*(3*Mgl*pow2(Msq)*(5046*pow2(Mst1)*pow2(Mst2) +
        11941*pow4(Mst1) - 54*pow4(Mst2)) + 60*Dmsqst2*Mgl*(117*pow2(Mst1)*
        pow2(Mst2) + 221*pow4(Mst1) + 27*pow4(Mst2)) - 4*Dmglst2*pow2(Msq)*(
        4266*pow2(Mst1)*pow2(Mst2) + 10598*pow4(Mst1) + 2133*pow4(Mst2)))) +
        pow2(Msq)*(36*pow2(Mt)*(-42*(2*Dmglst2 - 3*Mgl)*pow2(Mst1)*pow2(Mst2) +
        545*Mgl*pow4(Mst1)) - 2*Mst2*Mt*s2t*(36*(2845*Dmglst2 + 796*Mgl)*pow2(
        Mst1)*pow2(Mst2) + 10*(17281*Dmglst2 + 4101*Mgl)*pow4(Mst1) + 27*(1763*
        Dmglst2 + 555*Mgl)*pow4(Mst2))))) + Tbeta*pow2(Sbeta)*(4*Mt*pow2(Msq)*
        pow3(s2t)*(18*(401*Dmglst2 - 73*Mgl)*pow2(Mst1)*pow2(Mst2) + (15571*
        Dmglst2 - 1317*Mgl)*pow4(Mst1) + 27*(1763*Dmglst2 + 555*Mgl)*pow4(Mst2)
        )*pow5(Mst2) + (18*Mgl*(210*Dmsqst2 + 859*pow2(Msq))*pow2(Mst1)*pow2(
        Mst2) + (840*Dmsqst2*Mgl + (-16796*Dmglst2 + 5385*Mgl)*pow2(Msq))*pow4(
        Mst1) + 54*(30*Dmsqst2*Mgl - (158*Dmglst2 + 3*Mgl)*pow2(Msq))*pow4(
        Mst2))*pow4(s2t)*pow6(Mst2) - 4*pow2(Mst2)*pow2(Mt)*pow2(s2t)*(3*(160*
        Dmsqst2*Mgl + (-18456*Dmglst2 + 6857*Mgl)*pow2(Msq))*pow2(Mst2)*pow4(
        Mst1) + 18*(60*Dmsqst2*Mgl + (-600*Dmglst2 + 629*Mgl)*pow2(Msq))*pow2(
        Mst1)*pow4(Mst2) + pow2(MuSUSY)*(3*Mgl*pow2(Msq)*(5046*pow2(Mst1)*pow2(
        Mst2) + 11941*pow4(Mst1) - 54*pow4(Mst2)) + 60*Dmsqst2*Mgl*(117*pow2(
        Mst1)*pow2(Mst2) + 221*pow4(Mst1) + 27*pow4(Mst2)) - 4*Dmglst2*pow2(
        Msq)*(4266*pow2(Mst1)*pow2(Mst2) + 10598*pow4(Mst1) + 2133*pow4(Mst2)))
        + 81*Mgl*(20*Dmsqst2 + 53*pow2(Msq))*pow6(Mst2)) + 8*pow2(Msq)*pow3(Mt)
        *(-(Mst2*s2t*(40*(2573*Dmglst2 - 411*Mgl)*pow2(Mst2)*pow4(Mst1) + 18*(
        1383*Dmglst2 - 377*Mgl)*pow2(Mst1)*pow4(Mst2) + pow2(MuSUSY)*(36*(2845*
        Dmglst2 + 796*Mgl)*pow2(Mst1)*pow2(Mst2) + 10*(17281*Dmglst2 + 4101*
        Mgl)*pow4(Mst1) + 27*(1763*Dmglst2 + 555*Mgl)*pow4(Mst2)) + 378*(5*
        Dmglst2 - 3*Mgl)*pow6(Mst2))) + 2*Mt*(8*(2791*Dmglst2 + 438*Mgl)*pow2(
        Mst2)*pow4(Mst1) + 9*pow2(MuSUSY)*(-42*(2*Dmglst2 - 3*Mgl)*pow2(Mst1)*
        pow2(Mst2) + 545*Mgl*pow4(Mst1)) + 9*(644*Dmglst2 + 209*Mgl)*pow2(Mst1)
        *pow4(Mst2) + 189*(4*Dmglst2 + 3*Mgl)*pow6(Mst2)))))) + 8*z4*(-4*Mst2*
        xDmglst2*pow2(Dmglst2)*pow4(Msq)*((3*Mst2*pow2(Mst1)*(4*pow2(Mst2)*
        pow2(Mt)*(-10819*Tbeta*pow2(MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta)) + 4*
        Mt*(522*MuSUSY*s2t + 4703*Mt*Tbeta)*pow2(Sbeta)) - 2*Mt*s2t*(3*Mt*s2t*(
        8597*MuSUSY + 2388*Mst2*Tbeta) + 17228*Tbeta*pow2(Mt) + 8306*Mst2*
        MuSUSY*pow2(s2t))*pow2(Sbeta)*pow3(Mst2) + 32*Mst2*MuSUSY*(-1543*
        MuSUSY*s2t*Tbeta*(-1 + pow2(Sbeta)) + 1730*Mt*pow2(Sbeta))*pow3(Mt) +
        Tbeta*(-21648*pow2(MuSUSY)*(-1 + pow2(Sbeta))*pow4(Mt) + (9700*Mt -
        2513*Mst2*s2t)*pow2(Sbeta)*pow3(s2t)*pow5(Mst2))))/2. + 9*pow3(Mst2)*(-
        4*pow2(Mst2)*pow2(Mt)*(1111*Tbeta*pow2(MuSUSY)*pow2(s2t)*(-1 + pow2(
        Sbeta)) + 2*Mt*(423*MuSUSY*s2t - 91*Mt*Tbeta)*pow2(Sbeta)) + Mt*s2t*(-
        3747*Mt*MuSUSY*s2t + 1692*Mst2*Mt*s2t*Tbeta - 3484*Tbeta*pow2(Mt) -
        4444*Mst2*MuSUSY*pow2(s2t))*pow2(Sbeta)*pow3(Mst2) + 2*Mst2*MuSUSY*(-
        1249*MuSUSY*s2t*Tbeta*(-1 + pow2(Sbeta)) + 1742*Mt*pow2(Sbeta))*pow3(
        Mt) + Tbeta*(-1944*pow2(MuSUSY)*(-1 + pow2(Sbeta))*pow4(Mt) + (1249*Mt
        + 1111*Mst2*s2t)*pow2(Sbeta)*pow3(s2t)*pow5(Mst2))) + pow4(Mst1)*(2*
        Mst2*pow2(Mt)*pow2(s2t)*(-22654*Tbeta*pow2(MuSUSY)*(-1 + pow2(Sbeta)) -
        3*Mst2*(-5063*MuSUSY + 9228*Mst2*Tbeta)*pow2(Sbeta)) + 8*s2t*(-659*
        Tbeta*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 2*Mst2*(6813*MuSUSY + 12865*
        Mst2*Tbeta)*pow2(Sbeta))*pow3(Mt) + pow2(Sbeta)*(-(Mt*(25328*MuSUSY +
        15571*Mst2*Tbeta)*pow3(Mst2)*pow3(s2t)) - 16*(14269*MuSUSY + 5582*Mst2*
        Tbeta)*pow4(Mt) + 4199*Tbeta*pow4(s2t)*pow5(Mst2)))) + 60*xDmsqst2*
        pow2(Dmsqst2)*pow2(Mgl)*pow2(Mst2)*(4*MuSUSY*pow2(Sbeta)*(-3*Mt*pow2(
        Mst2)*pow3(s2t)*(19*pow2(Mst1)*(pow2(Mst1) + pow2(Mst2)) + 6*pow4(Mst2)
        ) + 2*s2t*pow3(Mt)*(21*pow2(Mst1)*pow2(Mst2) + 53*pow4(Mst1) + 18*pow4(
        Mst2))) + 3*Tbeta*pow2(Sbeta)*pow4(Mst2)*(13*pow2(Mst1)*pow2(Mst2) + 6*
        pow4(Mst2))*pow4(s2t) + 4*Tbeta*pow2(Mt)*pow2(s2t)*((-221*pow2(MuSUSY)*
        (-1 + pow2(Sbeta)) + 7*pow2(Mst2)*pow2(Sbeta))*pow4(Mst1) - 3*pow2(
        Mst1)*(25*pow2(Mst2)*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + pow2(Sbeta)*
        pow4(Mst2)) - 18*(pow2(MuSUSY)*(-1 + pow2(Sbeta))*pow4(Mst2) + pow2(
        Sbeta)*pow6(Mst2)))) + Mgl*pow2(Msq)*(4*Mt*MuSUSY*(Mst2*pow2(Sbeta)*(-
        3*Mt*pow2(Msq)*pow2(Mst2)*pow2(s2t)*(-45*(121*Dmglst2 + 301*Mgl)*pow2(
        Mst1)*pow2(Mst2) + 2*(5063*Dmglst2 - 7431*Mgl)*pow4(Mst1) - 27*(271*
        Dmglst2 + 255*Mgl)*pow4(Mst2)) + 4*pow2(Msq)*pow3(Mt)*(72*(78*Dmglst2 -
        163*Mgl)*pow2(Mst1)*pow2(Mst2) + 4*(14269*Dmglst2 - 4989*Mgl)*pow4(
        Mst1) - 27*(253*Dmglst2 + 309*Mgl)*pow4(Mst2)) - pow3(Mst2)*pow3(s2t)*(
        36*(150*Dmsqst2*Mgl + (-237*Dmglst2 + 479*Mgl)*pow2(Msq))*pow2(Mst1)*
        pow2(Mst2) + (6240*Dmsqst2*Mgl + (-25328*Dmglst2 + 20685*Mgl)*pow2(Msq)
        )*pow4(Mst1) + 108*(15*Dmsqst2*Mgl + 2*(-53*Dmglst2 + 15*Mgl)*pow2(Msq)
        )*pow4(Mst2)) + 6*Mst2*s2t*pow2(Mt)*(3*(300*Dmsqst2*Mgl + (96*Dmglst2 +
        1735*Mgl)*pow2(Msq))*pow2(Mst1)*pow2(Mst2) + 2*(530*Dmsqst2*Mgl + (-
        9084*Dmglst2 + 6031*Mgl)*pow2(Msq))*pow4(Mst1) + 27*(20*Dmsqst2*Mgl + (
        144*Dmglst2 + 101*Mgl)*pow2(Msq))*pow4(Mst2))) - Mt*MuSUSY*Tbeta*(2*Mt*
        pow2(Msq)*(18*Mt*(6*(94*Dmglst2 + 75*Mgl)*pow2(Mst1)*pow2(Mst2) + 1031*
        Mgl*pow4(Mst1) + 162*(2*Dmglst2 + Mgl)*pow4(Mst2)) + Mst2*s2t*(18*(709*
        Dmglst2 + 1135*Mgl)*pow2(Mst1)*pow2(Mst2) + 4*(659*Dmglst2 + 8823*Mgl)*
        pow4(Mst1) + 27*(271*Dmglst2 + 255*Mgl)*pow4(Mst2))) - pow2(Mst2)*pow2(
        s2t)*(60*Dmsqst2*Mgl*(117*pow2(Mst1)*pow2(Mst2) + 221*pow4(Mst1) + 27*
        pow4(Mst2)) + pow2(Msq)*(3*Mgl*(6828*pow2(Mst1)*pow2(Mst2) + 13723*
        pow4(Mst1) + 1080*pow4(Mst2)) - 4*Dmglst2*(4995*pow2(Mst1)*pow2(Mst2) +
        11327*pow4(Mst1) + 2862*pow4(Mst2)))))) + Tbeta*pow2(Sbeta)*(4*Mt*pow2(
        Msq)*pow3(s2t)*(36*(52*Dmglst2 - 185*Mgl)*pow2(Mst1)*pow2(Mst2) + (
        15571*Dmglst2 - 1317*Mgl)*pow4(Mst1) - 27*(271*Dmglst2 + 255*Mgl)*pow4(
        Mst2))*pow5(Mst2) + (36*(105*Dmsqst2*Mgl + (81*Dmglst2 + 389*Mgl)*pow2(
        Msq))*pow2(Mst1)*pow2(Mst2) + (840*Dmsqst2*Mgl + (-16796*Dmglst2 +
        3441*Mgl)*pow2(Msq))*pow4(Mst1) + 108*(15*Dmsqst2*Mgl + 2*(-53*Dmglst2
        + 15*Mgl)*pow2(Msq))*pow4(Mst2))*pow4(s2t)*pow6(Mst2) - 4*pow2(Mst2)*
        pow2(Mt)*pow2(s2t)*(3*(160*Dmsqst2*Mgl + (-18456*Dmglst2 + 6857*Mgl)*
        pow2(Msq))*pow2(Mst2)*pow4(Mst1) + 18*(60*Dmsqst2*Mgl + (-600*Dmglst2 +
        413*Mgl)*pow2(Msq))*pow2(Mst1)*pow4(Mst2) + pow2(MuSUSY)*(60*Dmsqst2*
        Mgl*(117*pow2(Mst1)*pow2(Mst2) + 221*pow4(Mst1) + 27*pow4(Mst2)) +
        pow2(Msq)*(3*Mgl*(6828*pow2(Mst1)*pow2(Mst2) + 13723*pow4(Mst1) + 1080*
        pow4(Mst2)) - 4*Dmglst2*(4995*pow2(Mst1)*pow2(Mst2) + 11327*pow4(Mst1)
        + 2862*pow4(Mst2)))) + 81*(20*Dmsqst2*Mgl + (144*Dmglst2 + 101*Mgl)*
        pow2(Msq))*pow6(Mst2)) + 8*pow2(Msq)*pow3(Mt)*(2*Mt*(8*(2791*Dmglst2 +
        438*Mgl)*pow2(Mst2)*pow4(Mst1) + 9*(644*Dmglst2 + 209*Mgl)*pow2(Mst1)*
        pow4(Mst2) + 9*pow2(MuSUSY)*(6*(94*Dmglst2 + 75*Mgl)*pow2(Mst1)*pow2(
        Mst2) + 1031*Mgl*pow4(Mst1) + 162*(2*Dmglst2 + Mgl)*pow4(Mst2)) + 189*(
        4*Dmglst2 + 3*Mgl)*pow6(Mst2)) + Mst2*s2t*(-40*(2573*Dmglst2 - 411*Mgl)
        *pow2(Mst2)*pow4(Mst1) - 18*(1383*Dmglst2 - 377*Mgl)*pow2(Mst1)*pow4(
        Mst2) + pow2(MuSUSY)*(18*(709*Dmglst2 + 1135*Mgl)*pow2(Mst1)*pow2(Mst2)
        + 4*(659*Dmglst2 + 8823*Mgl)*pow4(Mst1) + 27*(271*Dmglst2 + 255*Mgl)*
        pow4(Mst2)) + 54*(253*Dmglst2 + 309*Mgl)*pow6(Mst2)))))) - z3*(15*
        xDmsqst2*pow2(Dmsqst2)*pow2(Mst2)*(-5184*Mst2*xDmglst2*pow2(Dmglst2)*(
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
        pow2(Sbeta)*pow4(s2t)*pow5(Mst2))) + Mgl*(4*Mt*MuSUSY*pow2(Sbeta)*(
        41472*Mst2*(4*(5*Dmglst2 - Mgl)*pow2(Mst1) + 3*(2*Dmglst2 - Mgl)*pow2(
        Mst2))*pow3(Mt) + 4*Mgl*s2t*pow2(Mt)*(11856*pow2(Mst1)*pow2(Mst2) +
        15440*pow4(Mst1) + 9819*pow4(Mst2)) + 15552*Mst2*Mt*pow2(s2t)*(2*
        Dmglst2*pow2(2*pow2(Mst1) + pow2(Mst2)) - Mgl*(10*pow2(Mst1)*pow2(Mst2)
        + 8*pow4(Mst1) + 7*pow4(Mst2))) - 3*pow2(Mst2)*pow3(s2t)*(1728*Dmglst2*
        pow2(Mst2)*(5*pow2(Mst1) + 3*pow2(Mst2)) - Mgl*(12007*pow2(Mst1)*pow2(
        Mst2) + 8416*pow4(Mst1) + 16521*pow4(Mst2)))) + Tbeta*(41472*Mst2*s2t*(
        7*Mgl*(2*pow2(Mst1) + pow2(Mst2)) - 2*Dmglst2*(5*pow2(Mst1) + pow2(
        Mst2)))*pow2(MuSUSY)*pow3(Mt) + 4*Mt*pow2(s2t)*(648*pow2(Mst2)*(10*Mgl*
        Mst2*s2t*(-2 + pow2(Sbeta)) + (20*Dmglst2 - 3*Mgl)*Mt*pow2(Sbeta))*
        pow4(Mst1)*pow4(Sbeta) + Mt*(pow2(MuSUSY)*(48*(864*Dmglst2 - 1783*Mgl)*
        pow2(Mst1)*pow2(Mst2) + 28336*Mgl*pow4(Mst1) + 9*(1728*Dmglst2 - 5507*
        Mgl)*pow4(Mst2)) - 1296*(20*Dmglst2 - 3*Mgl)*pow2(Mst2)*pow4(Mst1)*
        pow4(Sbeta))) + pow2(Sbeta)*(-10368*Mt*pow3(Mst2)*pow3(s2t)*(6*(2*
        Dmglst2 - Mgl)*pow2(Mst1)*pow2(Mst2) + 5*Mgl*pow4(Mst1) + 2*(2*Dmglst2
        - 7*Mgl)*pow4(Mst2)) - 41472*Mst2*s2t*pow3(Mt)*((7*Mgl*(2*pow2(Mst1) +
        pow2(Mst2)) - 2*Dmglst2*(5*pow2(Mst1) + pow2(Mst2)))*pow2(MuSUSY) + 8*(
        (7*Dmglst2 - 2*Mgl)*pow2(Mst1)*pow2(Mst2) + (11*Dmglst2 - Mgl)*pow4(
        Mst1)) + 12*(2*Dmglst2 - Mgl)*pow4(Mst2)) - 1296*(-576*(4*Dmglst2 -
        Mgl)*pow2(Mst1)*pow2(Mst2) + 384*Mgl*pow4(Mst1) + (-640*Dmglst2 + 203*
        Mgl)*pow4(Mst2))*pow4(Mt) - 3*pow4(Mst2)*(1728*Dmglst2*(-2*pow2(Mst1)*
        pow2(Mst2) + 5*pow4(Mst1) - 3*pow4(Mst2)) + Mgl*(-4514*pow2(Mst1)*pow2(
        Mst2) - 3591*pow4(Mst1) + 16521*pow4(Mst2)))*pow4(s2t) + 4*pow2(Mt)*
        pow2(s2t)*(-(pow2(MuSUSY)*(48*(864*Dmglst2 - 1783*Mgl)*pow2(Mst1)*pow2(
        Mst2) + 28336*Mgl*pow4(Mst1) + 9*(1728*Dmglst2 - 5507*Mgl)*pow4(Mst2)))
        + 2*Mgl*(784*pow2(Mst2)*pow4(Mst1) - 2037*pow2(Mst1)*pow4(Mst2) - 9819*
        pow6(Mst2))))))) + pow2(Msq)*(-4*Mst2*xDmglst2*pow2(Dmglst2)*(19440*
        Dmsqst2*pow2(Mst2)*(s2t*pow2(Sbeta)*(-2*s2t*pow2(Mt)*(48*MuSUSY + 5*
        Mst2*Tbeta*(-2 + pow2(Sbeta))*pow2(Sbeta)) + 704*Tbeta*pow3(Mt) + 5*
        Tbeta*pow3(Mst2)*pow3(s2t))*pow4(Mst1) - 2*pow2(Mst1)*(16*Mst2*MuSUSY*
        pow2(Mt)*pow2(s2t)*(-(MuSUSY*Tbeta*(-1 + pow2(Sbeta))) + 3*Mst2*pow2(
        Sbeta)) + 8*s2t*Tbeta*(5*pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 28*pow2(
        Mst2)*pow2(Sbeta))*pow3(Mt) - 2*Mt*(5*MuSUSY + 6*Mst2*Tbeta)*pow2(
        Sbeta)*pow3(Mst2)*pow3(s2t) + 32*(10*MuSUSY + 9*Mst2*Tbeta)*pow2(Sbeta)
        *pow4(Mt) + Tbeta*pow2(Sbeta)*pow4(s2t)*pow5(Mst2)) - pow2(Mst2)*(12*
        Mst2*MuSUSY*pow2(Mt)*pow2(s2t)*(-(MuSUSY*Tbeta*(-1 + pow2(Sbeta))) + 2*
        Mst2*pow2(Sbeta)) + 16*s2t*Tbeta*(pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 12*
        pow2(Mst2)*pow2(Sbeta))*pow3(Mt) - 4*Mt*(3*MuSUSY + 2*Mst2*Tbeta)*pow2(
        Sbeta)*pow3(Mst2)*pow3(s2t) + 32*(6*MuSUSY + 5*Mst2*Tbeta)*pow2(Sbeta)*
        pow4(Mt) + 3*Tbeta*pow2(Sbeta)*pow4(s2t)*pow5(Mst2))) + pow2(Msq)*(6*
        Mst2*pow2(Mst1)*(-32*pow2(Mst2)*pow2(Mt)*(-101147*Tbeta*pow2(MuSUSY)*
        pow2(s2t)*(-1 + pow2(Sbeta)) + Mt*(56232*MuSUSY*s2t - 223217*Mt*Tbeta)*
        pow2(Sbeta)) + 2*Mt*s2t*(2*Mt*(827103*MuSUSY*s2t - 499172*Mt*Tbeta) +
        Mst2*s2t*(932531*MuSUSY*s2t - 82125*Mt*Tbeta))*pow2(Sbeta)*pow3(Mst2) +
        16*Mst2*MuSUSY*(304286*MuSUSY*s2t*Tbeta*(-1 + pow2(Sbeta)) + 221681*Mt*
        pow2(Sbeta))*pow3(Mt) + Tbeta*(-1309344*pow2(MuSUSY)*(-1 + pow2(Sbeta))
        *pow4(Mt) - 5*(-45736*Mt + 24671*Mst2*s2t)*pow2(Sbeta)*pow3(s2t)*pow5(
        Mst2))) - 9*pow3(Mst2)*(4*pow2(Mst2)*pow2(Mt)*(-228607*Tbeta*pow2(
        MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta)) + (354654*Mt*MuSUSY*s2t - 520966*
        Tbeta*pow2(Mt))*pow2(Sbeta)) - 4*Mt*s2t*(3*Mt*s2t*(221914*MuSUSY +
        59109*Mst2*Tbeta) - 258368*Tbeta*pow2(Mt) + 228607*Mst2*MuSUSY*pow2(
        s2t))*pow2(Sbeta)*pow3(Mst2) - 176*Mst2*MuSUSY*(10087*MuSUSY*s2t*Tbeta*
        (-1 + pow2(Sbeta)) + 5872*Mt*pow2(Sbeta))*pow3(Mt) + Tbeta*(382080*
        pow2(MuSUSY)*(-1 + pow2(Sbeta))*pow4(Mt) + (887656*Mt + 228607*Mst2*
        s2t)*pow2(Sbeta)*pow3(s2t)*pow5(Mst2))) + 2*pow4(Mst1)*(16*Mst2*pow2(
        Mt)*pow2(s2t)*(353116*Tbeta*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + (546087*
        Mst2*MuSUSY - 6606*Tbeta*pow2(Mst2))*pow2(Sbeta)) - 32*s2t*(-395581*
        Tbeta*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + (25254*Mst2*MuSUSY - 653582*
        Tbeta*pow2(Mst2))*pow2(Sbeta))*pow3(Mt) + pow2(Sbeta)*(52*Mt*(27676*
        MuSUSY - 21739*Mst2*Tbeta)*pow3(Mst2)*pow3(s2t) - 64*(417241*MuSUSY +
        591386*Mst2*Tbeta)*pow4(Mt) + 217283*Tbeta*pow4(s2t)*pow5(Mst2))))) +
        Mgl*(32*Mst2*Mt*MuSUSY*pow2(Sbeta)*(4*pow3(Mt)*(36*(2160*Dmsqst2*(5*
        Dmglst2 - Mgl) + (10050*Dmglst2 - 8219*Mgl)*pow2(Msq))*pow2(Mst1)*pow2(
        Mst2) + 4*(417241*Dmglst2 - 160608*Mgl)*pow2(Msq)*pow4(Mst1) + 81*(480*
        Dmsqst2*(3*Dmglst2 - Mgl) + (3*Dmglst2 - 961*Mgl)*pow2(Msq))*pow4(Mst2)
        ) + 3*Mt*pow2(Mst2)*pow2(s2t)*(9*(8640*Dmsqst2*(Dmglst2 - Mgl) - (
        49501*Dmglst2 + 11377*Mgl)*pow2(Msq))*pow2(Mst1)*pow2(Mst2) + 4*(19440*
        Dmsqst2*(Dmglst2 - Mgl) - (182029*Dmglst2 + 8862*Mgl)*pow2(Msq))*pow4(
        Mst1) + 27*(720*Dmsqst2*(Dmglst2 - 3*Mgl) - (15137*Dmglst2 + 5965*Mgl)*
        pow2(Msq))*pow4(Mst2)) - pow3(Mst2)*pow3(s2t)*(9*(Dmglst2*(5400*Dmsqst2
        + 64119*pow2(Msq)) + 4*Mgl*(780*Dmsqst2 + (11924 + 27*lmMst1 - 27*
        lmMst2)*pow2(Msq)))*pow2(Mst1)*pow2(Mst2) + 2*(21840*Dmsqst2*Mgl + 11*(
        16354*Dmglst2 + 20937*Mgl)*pow2(Msq))*pow4(Mst1) - 27*(690*Dmsqst2*Mgl
        - 6*(1319 + 6*lmMst1 - 6*lmMst2)*Mgl*pow2(Msq) - 5*Dmglst2*(216*Dmsqst2
        + 3523*pow2(Msq)))*pow4(Mst2)) + 6*Mst2*s2t*pow2(Mt)*(20*Dmsqst2*Mgl*(
        909*pow2(Mst1)*pow2(Mst2) + 965*pow4(Mst1) + 783*pow4(Mst2)) + pow2(
        Msq)*(12*Dmglst2*(3540*pow2(Mst1)*pow2(Mst2) + 2806*pow4(Mst1) + 1317*
        pow4(Mst2)) + Mgl*(-5577*pow2(Mst1)*pow2(Mst2) + 27086*pow4(Mst1) + 81*
        (-77 + 8*lmMst1 - 8*lmMst2)*pow4(Mst2))))) + 8*Tbeta*(4*pow2(Mst2)*
        pow2(Mt)*pow2(s2t)*(pow2(Msq)*pow2(MuSUSY)*(Dmglst2*(1052676*pow2(Mst1)
        *pow2(Mst2) + 1412464*pow4(Mst1) + 475605*pow4(Mst2)) + 6*Mgl*(3*(35719
        + 108*lmMst1 - 108*lmMst2)*pow2(Mst1)*pow2(Mst2) + 2*(91963 + 162*
        lmMst1 - 162*lmMst2)*pow4(Mst1) + 27*(1319 + 6*lmMst1 - 6*lmMst2)*pow4(
        Mst2))) + 30*Dmsqst2*(Mgl*pow2(MuSUSY)*(315*pow2(Mst1)*pow2(Mst2) +
        1771*pow4(Mst1) - 621*pow4(Mst2)) + 162*Dmglst2*pow2(Mst2)*(2*(8*pow2(
        Mst1) + 3*pow2(Mst2))*pow2(MuSUSY) + 5*(-2 + pow2(Sbeta))*pow4(Mst1)*
        pow4(Sbeta)))) + Mst2*s2t*pow2(Sbeta)*(4*Mt*pow2(s2t)*pow4(Mst2)*(-18*(
        1080*Dmsqst2*(3*Dmglst2 - Mgl) + (-2045*Dmglst2 + 3259*Mgl)*pow2(Msq))*
        pow2(Mst1)*pow2(Mst2) + (282607*Dmglst2 - 66945*Mgl)*pow2(Msq)*pow4(
        Mst1) - 27*(720*Dmsqst2*(Dmglst2 - 3*Mgl) - (15137*Dmglst2 + 5965*Mgl)*
        pow2(Msq))*pow4(Mst2)) - pow3(s2t)*(-18*(1080*Dmglst2*Dmsqst2 + 2595*
        Dmsqst2*Mgl + 5637*Dmglst2*pow2(Msq) + 11977*Mgl*pow2(Msq))*pow2(Mst1)*
        pow2(Mst2) + (Dmglst2*(48600*Dmsqst2 + 217283*pow2(Msq)) - 6*Mgl*(2600*
        Dmsqst2 + (5225 - 162*lmMst1 + 162*lmMst2)*pow2(Msq)))*pow4(Mst1) + 27*
        (690*Dmsqst2*Mgl - 6*(1319 + 6*lmMst1 - 6*lmMst2)*Mgl*pow2(Msq) - 5*
        Dmglst2*(216*Dmsqst2 + 3523*pow2(Msq)))*pow4(Mst2))*pow5(Mst2) - 8*
        pow3(Mt)*(8*(9720*Dmsqst2*(11*Dmglst2 - Mgl) + (326791*Dmglst2 - 86637*
        Mgl)*pow2(Msq))*pow2(Mst2)*pow4(Mst1) + 18*(4320*Dmsqst2*(7*Dmglst2 -
        Mgl) + 7*(5739*Dmglst2 - 3461*Mgl)*pow2(Msq))*pow2(Mst1)*pow4(Mst2) +
        pow2(MuSUSY)*(-72*(270*Dmsqst2*(5*Dmglst2 - 7*Mgl) - (11864*Dmglst2 +
        3659*Mgl)*pow2(Msq))*pow2(Mst1)*pow2(Mst2) + 4*(395581*Dmglst2 + 74724*
        Mgl)*pow2(Msq)*pow4(Mst1) - 27*(720*Dmsqst2*(Dmglst2 - 3*Mgl) - (15137*
        Dmglst2 + 5965*Mgl)*pow2(Msq))*pow4(Mst2)) + 162*(480*Dmsqst2*(3*
        Dmglst2 - Mgl) + (3*Dmglst2 - 961*Mgl)*pow2(Msq))*pow6(Mst2)) - 4*Mst2*
        s2t*pow2(Mt)*(3*(1120*Dmsqst2*Mgl + 367*(-24*Dmglst2 + 89*Mgl)*pow2(
        Msq))*pow2(Mst2)*pow4(Mst1) + 36*(210*Dmsqst2*Mgl + (2223*Dmglst2 + (55
        - 54*lmMst1 + 54*lmMst2)*Mgl)*pow2(Msq))*pow2(Mst1)*pow4(Mst2) + pow2(
        MuSUSY)*(30*Dmsqst2*(324*Dmglst2*pow2(Mst2)*(8*pow2(Mst1) + 3*pow2(
        Mst2)) + Mgl*(315*pow2(Mst1)*pow2(Mst2) + 1771*pow4(Mst1) - 621*pow4(
        Mst2))) + pow2(Msq)*(Dmglst2*(1052676*pow2(Mst1)*pow2(Mst2) + 1412464*
        pow4(Mst1) + 475605*pow4(Mst2)) + 6*Mgl*(3*(35719 + 108*lmMst1 - 108*
        lmMst2)*pow2(Mst1)*pow2(Mst2) + 2*(91963 + 162*lmMst1 - 162*lmMst2)*
        pow4(Mst1) + 27*(1319 + 6*lmMst1 - 6*lmMst2)*pow4(Mst2)))) + 27*(1740*
        Dmsqst2*Mgl + (1756*Dmglst2 + 9*(-77 + 8*lmMst1 - 8*lmMst2)*Mgl)*pow2(
        Msq))*pow6(Mst2))) + 8*pow3(Mt)*(pow2(MuSUSY)*(18*Mt*pow2(Msq)*(-6*(
        250*Dmglst2 - 327*Mgl)*pow2(Mst1)*pow2(Mst2) + 751*Mgl*pow4(Mst1) - 12*
        (83*Dmglst2 - 115*Mgl)*pow4(Mst2)) + Mst2*s2t*(-72*(270*Dmsqst2*(5*
        Dmglst2 - 7*Mgl) - (11864*Dmglst2 + 3659*Mgl)*pow2(Msq))*pow2(Mst1)*
        pow2(Mst2) + 4*(395581*Dmglst2 + 74724*Mgl)*pow2(Msq)*pow4(Mst1) - 27*(
        720*Dmsqst2*(Dmglst2 - 3*Mgl) - (15137*Dmglst2 + 5965*Mgl)*pow2(Msq))*
        pow4(Mst2))) + Mt*pow2(Sbeta)*(38880*Dmsqst2*(Dmglst2*(18*pow2(Mst1) +
        5*pow2(Mst2))*pow4(Mst2) - Mgl*(3*pow2(Mst2)*pow4(Mst1) + 2*pow2(Mst1)*
        pow4(Mst2) + pow6(Mst2))) + pow2(Msq)*(Dmglst2*pow2(Mst2)*(72*pow2(
        Mst1)*(14423*pow2(Mst2) + 375*pow2(MuSUSY)) + 27*pow2(Mst2)*(-9017*
        pow2(Mst2) + 664*pow2(MuSUSY)) + 4731088*pow4(Mst1)) - 6*Mgl*((193928*
        pow2(Mst2) + 2253*pow2(MuSUSY))*pow4(Mst1) + 4140*pow2(MuSUSY)*pow4(
        Mst2) + 3*pow2(Mst1)*(1962*pow2(Mst2)*pow2(MuSUSY) + 34525*pow4(Mst2))
        + 108*(514 - 3*(lmMst1 + lmMst2) + 6*lmMt)*pow6(Mst2)))))))))))/(23328.
        *Tbeta*pow2(Sbeta)*pow4(Msq)*pow6(Mst2)) - (Al4p*z2*(-8*xDmglst2*pow2(
        Dmglst2)*((-17010*twoLoopFlag*(4*pow2(Mt)*(8*(4*MuSUSY + 9*Mst2*Tbeta)*
        pow2(Mt)*pow2(Sbeta) + 3*Mst2*MuSUSY*pow2(s2t)*(-(MuSUSY*Tbeta*(-1 +
        pow2(Sbeta))) + 10*Mst2*pow2(Sbeta)) + 2*Mt*s2t*Tbeta*(31*pow2(MuSUSY)*
        (-1 + pow2(Sbeta)) - 12*pow2(Mst2)*pow2(Sbeta)))*pow4(Mst1) + pow2(
        Mst1)*pow2(Mst2)*(12*Mst2*MuSUSY*pow2(Mt)*pow2(s2t)*(MuSUSY*Tbeta + 10*
        Mst2*pow2(Sbeta) - MuSUSY*Tbeta*pow2(Sbeta)) + 168*s2t*Tbeta*pow2(
        MuSUSY)*(-1 + pow2(Sbeta))*pow3(Mt) + 4*Mt*Tbeta*pow2(Sbeta)*pow3(s2t)*
        pow4(Mst2) + 32*(MuSUSY + Mst2*Tbeta)*pow2(Sbeta)*pow4(Mt) - 3*Tbeta*
        pow2(Sbeta)*pow4(s2t)*pow5(Mst2)) + pow4(Mst2)*(12*Mst2*MuSUSY*pow2(Mt)
        *pow2(s2t)*(-(MuSUSY*Tbeta*(-1 + pow2(Sbeta))) + 11*Mst2*pow2(Sbeta)) +
        8*s2t*Tbeta*(11*pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 4*pow2(Mst2)*pow2(
        Sbeta))*pow3(Mt) - 4*Mt*(3*MuSUSY + 11*Mst2*Tbeta)*pow2(Sbeta)*pow3(
        Mst2)*pow3(s2t) + 32*MuSUSY*pow2(Sbeta)*pow4(Mt) + 3*Tbeta*pow2(Sbeta)*
        pow4(s2t)*pow5(Mst2))))/(Tbeta*pow2(Sbeta)*pow5(Mst2)) + Al4p*
        threeLoopFlag*((7*Mt*pow3(s2t)*(194400*Dmsqst2*(2*(5 - 3*lmMst1 + 3*
        lmMst2)*pow2(Mst1)*pow2(Mst2) + 11*pow4(Mst1) - 2*(2 + 3*lmMst1 - 3*
        lmMst2)*pow4(Mst2)) + pow2(Msq)*(-60*(143005 + 30024*lmMst1 - 26352*
        lmMst2)*pow2(Mst1)*pow2(Mst2) + (8583283 - 2329560*lmMst1 + 2329560*
        lmMst2)*pow4(Mst1) + 9*(897737 + 47880*lmMst1 + 238680*lmMst2)*pow4(
        Mst2))))/(2.*Mst2*pow2(Msq)) + (4*pow4(Mt)*(680400*Dmsqst2*pow2(Mst2)*(
        (23 - 5*lmMst1 + 10*lmMst2 - 5*lmMt)*pow2(Mst1)*pow2(Mst2) + (77 - 18*
        lmMst1 + 36*lmMst2 - 18*lmMt)*pow4(Mst1) + pow4(Mst2)) + pow2(Msq)*(-
        63*pow2(Mst1)*pow2(Mst2)*(5*(66191 - 9288*lmMst1 + 11952*lmMst2 -
        16488*lmMt)*pow2(Mst2) + 48*(617 + 2040*lmMst2)*pow2(MuSUSY)) - 3*((
        52628795 - 5738040*lmMst1 + 11332440*lmMst2 - 5972400*lmMt)*pow2(Mst2)
        + 21*(122051 + 1800*lmMst1 + 274680*lmMst2)*pow2(MuSUSY))*pow4(Mst1) +
        8*(19152737 - 3674160*lmMst1 + 6872040*lmMst2 - 850500*lmMt)*pow6(Mst1)
        - 272160*pow6(Mst2))))/(pow2(Msq)*pow2(Mst1)*pow4(Mst2)) - (s2t*pow3(
        Mt)*(14*(4*(9908167 - 949320*lmMst1 + 1769040*lmMst2 - 217080*lmMt)*
        pow2(Mst2) + (15057833 - 4014360*lmMst1 + 5563080*lmMst2)*pow2(MuSUSY))
        *pow4(Mst1) + 63*(897737 + 47880*lmMst1 + 238680*lmMst2)*pow2(MuSUSY)*
        pow4(Mst2) - 6*pow2(Mst1)*(7*(-1263161 + 156600*lmMst1 - 979560*lmMst2)
        *pow2(Mst2)*pow2(MuSUSY) + (47845349 - 6357960*lmMst1 + 13388760*lmMst2
        - 7030800*lmMt)*pow4(Mst2)) + (2721600*Dmsqst2*pow2(Mst2)*((-2 - 3*
        lmMst1 + 3*lmMst2)*pow2(Mst2)*pow2(MuSUSY) + pow2(Mst1)*((51 - 14*
        lmMst1 + 28*lmMst2 - 14*lmMt)*pow2(Mst2) + (1 - 9*lmMst1 + 9*lmMst2)*
        pow2(MuSUSY)) + (87 - 22*lmMst1 + 44*lmMst2 - 22*lmMt)*pow4(Mst1) + (19
        - 6*lmMst1 + 12*lmMst2 - 6*lmMt)*pow4(Mst2)))/pow2(Msq) + 18*(-2661049
        + 335160*lmMst1 - 1408680*lmMst2 + 408240*lmMt)*pow6(Mst2)))/pow5(Mst2)
        - 51030*pow4(s2t)*(-((31.72798353909465 + (41*lmMst1)/6. + (15*lmMst2)/
        2. + (20*Dmsqst2*(-4 + lmMst1 - lmMst2))/(3.*pow2(Msq)))*pow2(Mst1)*
        pow2(Mst2)) + (155.89787379972566 - 8*lmMst1 + 8*lmMst2 + (5*Dmsqst2*(7
        + 10*lmMst1 - 10*lmMst2))/(3.*pow2(Msq)))*pow4(Mst1) + (
        24.781172839506173 + (92*lmMst1)/9. + (37*lmMst2)/9. - (10*Dmsqst2*(4 +
        3*lmMst1 - 3*lmMst2))/(3.*pow2(Msq)))*pow4(Mst2) + (2*(2 - (5*Dmsqst2)/
        pow2(Msq))*pow6(Mst2))/(3.*pow2(Mst1))) - (pow2(Mt)*pow2(s2t)*(680400*
        Dmsqst2*pow2(Mst2)*((8*pow2(Mst2) + (3 + 8*lmMst1 - 8*lmMst2)*pow2(
        MuSUSY))*pow4(Mst1) + (2*pow2(Mst2) + pow2(MuSUSY))*pow4(Mst2) + pow2(
        Mst1)*(3*(2 + lmMst1 - lmMst2)*pow2(Mst2)*pow2(MuSUSY) + 8*pow4(Mst2))
        + 8*pow6(Mst1)) - pow2(Msq)*(-30*pow2(Mst2)*(9*(5208*lmMst1 - 37*(7339
        + 3192*lmMst2))*pow2(Mst2) - 7*(21223 + 13230*lmMst1 + 702*lmMst2)*
        pow2(MuSUSY))*pow4(Mst1) + 9*pow2(Mst1)*(18*(123227 - 1820*lmMst1 +
        230300*lmMst2)*pow2(Mst2) + 7*(88931 + 33120*lmMst1 + 13320*lmMst2)*
        pow2(MuSUSY))*pow4(Mst2) + 28*(45*(78533 + 6732*lmMst1 + 180*lmMst2)*
        pow2(Mst2) + (2128489 - 307800*lmMst1 + 377460*lmMst2)*pow2(MuSUSY))*
        pow6(Mst1) + 272160*(2*pow2(Mst2) + pow2(MuSUSY))*pow6(Mst2))))/(pow2(
        Msq)*pow2(Mst1)*pow4(Mst2)) + (7*pow2(Mt)*(Mt*pow2(MuSUSY)*(388800*
        Dmsqst2*s2t*pow2(Mst2)*((1 - 9*lmMst1 + 9*lmMst2)*pow2(Mst1) + (-2 - 3*
        lmMst1 + 3*lmMst2)*pow2(Mst2)) + pow2(Msq)*(6*Mst2*(6*(122051 + 1800*
        lmMst1 + 274680*lmMst2)*Mt + (1263161 - 156600*lmMst1 + 979560*lmMst2)*
        Mst2*s2t)*pow2(Mst1) + 9*(192*(617 + 2040*lmMst2)*Mt + (897737 + 47880*
        lmMst1 + 238680*lmMst2)*Mst2*s2t)*pow3(Mst2) + 2*(15057833 - 4014360*
        lmMst1 + 5563080*lmMst2)*s2t*pow4(Mst1))) - (Mst2*pow2(s2t)*(-24300*
        Dmsqst2*pow2(Mst2)*(4*pow2(MuSUSY)*(3*(2 + lmMst1 - lmMst2)*pow2(Mst1)*
        pow2(Mst2) + (3 + 8*lmMst1 - 8*lmMst2)*pow4(Mst1) + pow4(Mst2)) + (-93
        + 10*lmMst1 - 10*lmMst2)*(-2 + pow2(Sbeta))*pow4(Sbeta)*pow6(Mst1)) +
        pow2(Msq)*pow2(MuSUSY)*(30*(21223 + 13230*lmMst1 + 702*lmMst2)*pow2(
        Mst2)*pow4(Mst1) + 9*(88931 + 33120*lmMst1 + 13320*lmMst2)*pow2(Mst1)*
        pow4(Mst2) + 4*(2128489 - 307800*lmMst1 + 377460*lmMst2)*pow6(Mst1) +
        38880*pow6(Mst2))))/pow2(Mst1)))/(pow2(Msq)*pow2(Sbeta)*pow5(Mst2)) + (
        51030*Mt*MuSUSY*(-(Mt*pow2(s2t)*(Mst2*(1662.4759259259258 + (266*
        lmMst1)/3. + 442*lmMst2 - (80*Dmsqst2*(2 + 3*lmMst1 - 3*lmMst2))/pow2(
        Msq)) - ((103.0179012345679 + 282*lmMst1 - (2302*lmMst2)/3. + (240*
        Dmsqst2*(-1 + 2*lmMst1 - 2*lmMst2))/pow2(Msq))*pow2(Mst1))/Mst2 + ((
        3498.187242798354 - 996*lmMst1 + (3580*lmMst2)/3. + (40*Dmsqst2*(17 -
        12*lmMst1 + 12*lmMst2))/pow2(Msq))*pow4(Mst1))/pow3(Mst2))) - s2t*pow2(
        Mt)*(803.7269841269841 - (104*lmMst1)/9. + (13160*lmMst2)/9. - (800*
        Dmsqst2)/(3.*pow2(Msq)) + ((3677.1978835978834 - (200*lmMst1)/3. +
        2712*lmMst2 - (480*Dmsqst2)/pow2(Msq))*pow2(Mst1))/pow2(Mst2) + (32*(2
        - (5*Dmsqst2)/pow2(Msq))*pow2(Mst2))/(3.*pow2(Mst1)) + (4*(510139 +
        44280*lmMst1 + 76680*lmMst2)*pow4(Mst1))/(405.*pow4(Mst2))) + (pow3(Mt)
        *(48*(-3231296 + 460215*lmMst1 - 1100925*lmMst2 + 515970*lmMt)*pow2(
        Mst1)*pow2(Mst2) - (2721600*Dmsqst2*pow2(Mst2)*(10*(-3 + lmMst1 - 2*
        lmMst2 + lmMt)*pow2(Mst1) + (-8 + 3*lmMst1 - 6*lmMst2 + 3*lmMt)*pow2(
        Mst2)))/pow2(Msq) + 28*(13463149 - 1492020*lmMst1 + 2713500*lmMst2 -
        625320*lmMt)*pow4(Mst1) + 9*(-2456929 + 335160*lmMst1 - 1408680*lmMst2
        + 408240*lmMt)*pow4(Mst2)))/(25515.*pow5(Mst2)) + (pow3(s2t)*(pow2(Msq)
        *(3*(-54563 + 32940*lmMst1 - 32940*lmMst2)*pow2(Mst2)*pow4(Mst1) + 9*(
        84611 + 33120*lmMst1 + 13320*lmMst2)*pow2(Mst1)*pow4(Mst2) + 8*(788723
        - 80595*lmMst1 + 80595*lmMst2)*pow6(Mst1) + 38880*pow6(Mst2)) + 48600*
        Dmsqst2*(2*(3 - 5*lmMst1 + 5*lmMst2)*pow2(Mst2)*pow4(Mst1) + 13*pow6(
        Mst1) - 2*((5 + 3*lmMst1 - 3*lmMst2)*pow2(Mst1)*pow4(Mst2) + pow6(Mst2)
        ))))/(7290.*pow2(Msq)*pow2(Mst1)*pow2(Mst2))))/Tbeta)) + (420*(-1296*
        twoLoopFlag*xMst*pow2(Mt)*(-(xDmglst2*pow2(Dmglst2)*(16*(7*MuSUSY + 15*
        Mst2*Tbeta)*pow2(Mt)*pow2(Sbeta) + 3*Mst2*MuSUSY*pow2(s2t)*(-(MuSUSY*
        Tbeta*(-1 + pow2(Sbeta))) + 10*Mst2*pow2(Sbeta)) + 2*Mt*s2t*Tbeta*(41*
        pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 40*pow2(Mst2)*pow2(Sbeta)))) + 2*
        Dmglst2*Mgl*(16*(2*MuSUSY + 3*Mst2*Tbeta)*pow2(Mt)*pow2(Sbeta) + Mst2*
        MuSUSY*pow2(s2t)*(MuSUSY*Tbeta*(-1 + pow2(Sbeta)) - 6*Mst2*pow2(Sbeta))
        + Mt*s2t*Tbeta*(-17*pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 20*pow2(Mst2)*
        pow2(Sbeta))) - 2*Mt*pow2(Mgl)*(8*Mt*(2*MuSUSY + Mst2*Tbeta)*pow2(
        Sbeta) + s2t*Tbeta*(pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 4*pow2(Mst2)*
        pow2(Sbeta))))*pow8(Mst1) - (Al4p*threeLoopFlag*pow2(Mst2)*(108*
        xDR2DRMOD*pow4(Msq)*(-(xDmglst2*pow2(Dmglst2)*(Mt*MuSUSY*pow2(Sbeta)*(
        384*Mt*pow2(Mst1)*pow2(Mst2)*pow2(s2t)*((-1 + 6*lmMst2)*pow2(Mst2)*(
        pow2(Mst1) + pow2(Mst2)) + (1 + 3*lmMst2)*pow4(Mst1)) + 512*pow2(Mst1)*
        pow3(Mt)*(-(lmMst2*pow2(Mst1)*pow2(Mst2)) + 2*(11 + 5*lmMst2)*pow4(
        Mst1) - (-2 + lmMst2)*pow4(Mst2)) + 864*s2t*(-pow2(Mst1) + pow2(Mst2))*
        pow2(Mt)*pow5(Mst2) + 9*pow3(s2t)*((16*(8 + 27*lmMst1 - 39*lmMst2)*
        pow2(Mst1)*pow2(Mst2))/9. + 24*pow4(Mst1) - 24*pow4(Mst2))*pow5(Mst2))
        + Tbeta*(8*s2t*pow2(Mt)*pow2(MuSUSY)*(2*(-32*(-1 + 6*lmMst2)*Mt + (-8 -
        27*lmMst1 + 39*lmMst2)*Mst2*s2t)*pow2(Mst2)*pow4(Mst1) + (32*(1 - 6*
        lmMst2)*Mt + (11 - 54*lmMst1 + 78*lmMst2)*Mst2*s2t)*pow2(Mst1)*pow4(
        Mst2) - 2*(48*(1 + 3*lmMst2)*Mt + (7*lmMst1 - 15*lmMst2)*Mst2*s2t)*
        pow6(Mst1) + 27*s2t*pow7(Mst2)) + pow2(Sbeta)*(32*Mst2*pow4(Mt)*(-48*(4
        + lmMst2)*pow2(Mst2)*pow4(Mst1) + 43*pow2(Mst1)*pow4(Mst2) + 192*(3 +
        2*lmMst2)*pow6(Mst1) + 27*pow6(Mst2)) + pow4(s2t)*pow5(Mst2)*(2*(-11 +
        54*lmMst1 - 78*lmMst2)*pow2(Mst2)*pow4(Mst1) - 2*(43 + 54*lmMst1 - 78*
        lmMst2)*pow2(Mst1)*pow4(Mst2) + 14*pow6(Mst1) + 54*pow6(Mst2)) - 256*
        s2t*pow2(Mst1)*pow3(Mt)*(3*(2*(5 + 3*lmMst2)*pow2(Mst2) - (1 + 3*
        lmMst2)*pow2(MuSUSY))*pow4(Mst1) + (1 - 6*lmMst2)*pow2(MuSUSY)*pow4(
        Mst2) - 2*pow2(Mst1)*((-1 + 6*lmMst2)*pow2(Mst2)*pow2(MuSUSY) + 2*pow4(
        Mst2)) - 2*(-2 + lmMst2)*pow6(Mst2)) + 8*Mst2*pow2(Mt)*pow2(s2t)*(pow4(
        Mst1)*(2*(8 + 27*lmMst1 - 39*lmMst2)*pow2(Mst2)*pow2(MuSUSY) - 54*pow4(
        Mst2)) + 2*(7*lmMst1 - 15*lmMst2)*pow2(MuSUSY)*pow6(Mst1) - 27*(2*pow2(
        Mst2) + pow2(MuSUSY))*pow6(Mst2) + pow2(Mst1)*((-11 + 54*lmMst1 - 78*
        lmMst2)*pow2(MuSUSY)*pow4(Mst2) + 108*pow6(Mst2))) + 128*(1 - 6*lmMst2)
        *Mt*pow2(Mst1)*pow3(s2t)*pow8(Mst2))))) + Mgl*(Mt*MuSUSY*pow2(Sbeta)*(-
        384*(Dmglst2 + 3*Dmglst2*lmMst2 + Mgl + lmMst2*Mgl)*Mt*pow2(Mst1)*pow2(
        Mst2)*pow2(s2t)*(pow2(Mst1)*pow2(Mst2) + pow4(Mst1) + pow4(Mst2)) +
        512*pow2(Mst1)*pow3(Mt)*(-(Dmglst2*(7 + lmMst2)*pow2(Mst1)*pow2(Mst2))
        - 2*Dmglst2*(11 + 5*lmMst2)*pow4(Mst1) + Dmglst2*(-1 + lmMst2)*pow4(
        Mst2) + (1 + lmMst2)*Mgl*(3*pow2(Mst1)*pow2(Mst2) + 6*pow4(Mst1) +
        pow4(Mst2))) + 32*(7*Dmglst2 - 13*Mgl)*s2t*(pow2(Mst1) - pow2(Mst2))*
        pow2(Mt)*pow5(Mst2) - pow3(s2t)*(16*(Dmglst2*(7*lmMst1 - 15*lmMst2) + (
        -4 - 13*lmMst1 + 9*lmMst2)*Mgl)*pow2(Mst1)*pow2(Mst2) + 8*(7*Dmglst2 -
        13*Mgl)*pow4(Mst1) - 8*(7*Dmglst2 - 13*Mgl)*pow4(Mst2))*pow5(Mst2)) +
        Tbeta*(s2t*pow2(Mt)*pow2(MuSUSY)*(256*(Dmglst2 + 3*Dmglst2*lmMst2 + Mgl
        + lmMst2*Mgl)*Mt*pow2(Mst1)*(2*pow2(Mst1)*pow2(Mst2) + 3*pow4(Mst1) +
        pow4(Mst2)) - 8*Mst2*s2t*(-2*(Dmglst2*(7*lmMst1 - 15*lmMst2) + (-4 -
        13*lmMst1 + 9*lmMst2)*Mgl)*(pow2(Mst1) + pow2(Mst2))*pow4(Mst1) + (
        Dmglst2*(7 - 14*lmMst1 + 30*lmMst2) + (-5 + 26*lmMst1 - 18*lmMst2)*Mgl)
        *pow2(Mst1)*pow4(Mst2) + (7*Dmglst2 - 13*Mgl)*pow6(Mst2))) + pow2(
        Sbeta)*(-2*(pow2(Mst1) - pow2(Mst2))*(Dmglst2*(2*(7*lmMst1 - 15*lmMst2)
        *pow2(Mst1)*pow2(Mst2) + 7*pow4(Mst1) - 7*pow4(Mst2)) + Mgl*(-2*(4 +
        13*lmMst1 - 9*lmMst2)*pow2(Mst1)*pow2(Mst2) - 13*pow4(Mst1) + 13*pow4(
        Mst2)))*pow4(s2t)*pow5(Mst2) + 256*s2t*pow2(Mst1)*pow3(Mt)*(6*(Dmglst2*
        (5 + 3*lmMst2) - (1 + lmMst2)*Mgl)*pow2(Mst2)*pow4(Mst1) + 4*(Dmglst2*(
        3 + lmMst2) - (1 + lmMst2)*Mgl)*pow2(Mst1)*pow4(Mst2) - (Dmglst2 + 3*
        Dmglst2*lmMst2 + Mgl + lmMst2*Mgl)*pow2(MuSUSY)*(2*pow2(Mst1)*pow2(
        Mst2) + 3*pow4(Mst1) + pow4(Mst2)) - 2*(Dmglst2*(-1 + lmMst2) + (1 +
        lmMst2)*Mgl)*pow6(Mst2)) - 32*Mst2*pow4(Mt)*(Dmglst2*(96*(2 + lmMst2)*
        pow2(Mst2)*pow4(Mst1) + 39*pow2(Mst1)*pow4(Mst2) + 192*(3 + 2*lmMst2)*
        pow6(Mst1) + 7*pow6(Mst2)) - Mgl*(48*(1 + lmMst2)*(2*pow2(Mst1) + pow2(
        Mst2))*pow4(Mst1) + (29 + 16*lmMst2)*pow2(Mst1)*pow4(Mst2) + 13*pow6(
        Mst2))) + 8*Mst2*pow2(Mt)*pow2(s2t)*(2*(7*Dmglst2 - 13*Mgl)*pow2(pow2(
        Mst1) - pow2(Mst2))*pow4(Mst2) + pow2(MuSUSY)*(-2*(Dmglst2*(7*lmMst1 -
        15*lmMst2) + (-4 - 13*lmMst1 + 9*lmMst2)*Mgl)*(pow2(Mst1) + pow2(Mst2))
        *pow4(Mst1) + (Dmglst2*(7 - 14*lmMst1 + 30*lmMst2) + (-5 + 26*lmMst1 -
        18*lmMst2)*Mgl)*pow2(Mst1)*pow4(Mst2) + (7*Dmglst2 - 13*Mgl)*pow6(Mst2)
        )) + 128*(Dmglst2 + 3*Dmglst2*lmMst2 + Mgl + lmMst2*Mgl)*Mt*pow2(Mst1)*
        pow3(s2t)*pow8(Mst2))))) + Mst2*xDmsqst2*pow2(Dmsqst2)*(324*Mst2*
        xDmglst2*pow2(Dmglst2)*(Mt*MuSUSY*pow2(Sbeta)*(320*pow2(Mst1)*(-10*(-3
        + lmMst1 - 2*lmMst2 + lmMt)*pow2(Mst1) + (8 - 3*lmMst1 + 6*lmMst2 - 3*
        lmMt)*pow2(Mst2))*pow3(Mt) + 160*Mst2*s2t*pow2(Mt)*(5*pow2(Mst1)*pow2(
        Mst2) + 9*pow4(Mst1) + pow4(Mst2)) - 3*Mt*pow2(Mst1)*pow2(s2t)*(240*(1
        - 2*lmMst1 + 2*lmMst2)*pow2(Mst1)*pow2(Mst2) + (680 - 480*lmMst1 + 480*
        lmMst2)*pow4(Mst1) - 80*(2 + 3*lmMst1 - 3*lmMst2)*pow4(Mst2)) + 20*
        Mst2*pow3(s2t)*(2*(3 - 5*lmMst1 + 5*lmMst2)*pow2(Mst2)*pow4(Mst1) + 13*
        pow6(Mst1) - 2*((5 + 3*lmMst1 - 3*lmMst2)*pow2(Mst1)*pow4(Mst2) + pow6(
        Mst2)))) + Tbeta*(160*pow3(Mt)*(s2t*pow2(Mst1)*((1 - 9*lmMst1 + 9*
        lmMst2)*pow2(Mst1) + (-2 - 3*lmMst1 + 3*lmMst2)*pow2(Mst2))*pow2(
        MuSUSY) + Mst2*Mt*pow2(Sbeta)*((23 - 5*lmMst1 + 10*lmMst2 - 5*lmMt)*
        pow2(Mst1)*pow2(Mst2) + (77 - 18*lmMst1 + 36*lmMst2 - 18*lmMt)*pow4(
        Mst1) + pow4(Mst2))) + Mst2*pow2(Mt)*pow2(s2t)*(40*pow2(MuSUSY)*(3*(2 +
        lmMst1 - lmMst2)*pow2(Mst1)*pow2(Mst2) + (3 + 8*lmMst1 - 8*lmMst2)*
        pow4(Mst1) + pow4(Mst2)) + 10*(-93 + 10*lmMst1 - 10*lmMst2)*(-2 + pow2(
        Sbeta))*pow4(Sbeta)*pow6(Mst1)) + s2t*pow2(Sbeta)*(Mt*pow2(Mst1)*pow2(
        Mst2)*pow2(s2t)*(80*(5 - 3*lmMst1 + 3*lmMst2)*pow2(Mst1)*pow2(Mst2) +
        440*pow4(Mst1) - 80*(2 + 3*lmMst1 - 3*lmMst2)*pow4(Mst2)) - 160*pow2(
        Mst1)*pow3(Mt)*((-2 - 3*lmMst1 + 3*lmMst2)*pow2(Mst2)*pow2(MuSUSY) +
        pow2(Mst1)*((51 - 14*lmMst1 + 28*lmMst2 - 14*lmMt)*pow2(Mst2) + (1 - 9*
        lmMst1 + 9*lmMst2)*pow2(MuSUSY)) + (87 - 22*lmMst1 + 44*lmMst2 - 22*
        lmMt)*pow4(Mst1) + (19 - 6*lmMst1 + 12*lmMst2 - 6*lmMt)*pow4(Mst2)) -
        40*Mst2*s2t*pow2(Mt)*(8*pow2(Mst1)*(pow2(Mst1)*pow2(Mst2) + pow4(Mst1)
        + pow4(Mst2)) + pow2(MuSUSY)*(3*(2 + lmMst1 - lmMst2)*pow2(Mst1)*pow2(
        Mst2) + (3 + 8*lmMst1 - 8*lmMst2)*pow4(Mst1) + pow4(Mst2)) + 2*pow6(
        Mst2)) - 5*pow3(Mst2)*pow3(s2t)*(4*(4 - lmMst1 + lmMst2)*pow2(Mst2)*
        pow4(Mst1) + (7 + 10*lmMst1 - 10*lmMst2)*pow6(Mst1) - 2*((4 + 3*lmMst1
        - 3*lmMst2)*pow2(Mst1)*pow4(Mst2) + pow6(Mst2)))))) - 1620*shiftst1*
        Tbeta*pow2(Mgl)*(4*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*(pow2(Mst1)*pow4(
        Mst2) - 2*(lmMst1 - lmMst2)*pow2(Mst1)*(pow2(Mst1)*pow2(Mst2) + pow4(
        Mst1) + pow4(Mst2)) + pow6(Mst2)) + pow2(Sbeta)*(16*(pow2(Mst1) + pow2(
        Mst2))*pow4(Mst2)*pow4(Mt) + (pow2(Mst1) - pow2(Mst2))*(2*(lmMst1 -
        lmMst2)*pow2(Mst1)*pow2(Mst2) + pow4(Mst1) - pow4(Mst2))*pow4(Mst2)*
        pow4(s2t) - 4*pow2(Mt)*pow2(s2t)*(2*pow4(Mst1)*pow4(Mst2) - 4*pow2(
        Mst1)*pow6(Mst2) + pow2(MuSUSY)*(pow2(Mst1)*pow4(Mst2) - 2*(lmMst1 -
        lmMst2)*pow2(Mst1)*(pow2(Mst1)*pow2(Mst2) + pow4(Mst1) + pow4(Mst2)) +
        pow6(Mst2)) + 2*pow8(Mst2)))) + Mgl*(4*MuSUSY*pow2(Sbeta)*(81*Mgl*
        shiftst1*(80*s2t*(pow2(Mst1) - pow2(Mst2))*pow3(Mt) - 20*Mt*pow3(s2t)*(
        2*(lmMst1 - lmMst2)*pow2(Mst1)*pow2(Mst2) + pow4(Mst1) - pow4(Mst2)))*
        pow4(Mst2) - 9720*Mst2*pow2(Mst1)*pow2(Mt)*pow2(s2t)*(2*(Dmglst2*(3 -
        6*lmMst1 + 6*lmMst2) + (4 + 3*lmMst1 - 3*lmMst2)*Mgl)*pow2(Mst1)*pow2(
        Mst2) + (Dmglst2*(17 - 12*lmMst1 + 12*lmMst2) + (5 + 4*lmMst1 - 4*
        lmMst2)*Mgl)*pow4(Mst1) - (Dmglst2*(4 + 6*lmMst1 - 6*lmMst2) + (-7 - 3*
        lmMst1 + 3*lmMst2)*Mgl)*pow4(Mst2)) + 12960*Mst2*pow2(Mst1)*(4*(-5*
        Dmglst2*(-3 + lmMst1 - 2*lmMst2 + lmMt) + (-1 + lmMst1 - 2*lmMst2 +
        lmMt)*Mgl)*pow2(Mst1) + (2*Dmglst2*(8 - 3*lmMst1 + 6*lmMst2 - 3*lmMt) +
        (-2 + 3*lmMst1 - 6*lmMst2 + 3*lmMt)*Mgl)*pow2(Mst2))*pow4(Mt) + 10*s2t*
        pow3(Mt)*(1296*Dmglst2*pow2(Mst2)*(5*pow2(Mst1)*pow2(Mst2) + 9*pow4(
        Mst1) + pow4(Mst2)) + Mgl*(-795*pow2(Mst2)*pow4(Mst1) - 450*pow2(Mst1)*
        pow4(Mst2) + 4057*pow6(Mst1) - 972*pow6(Mst2))) + 15*Mt*pow2(Mst2)*
        pow3(s2t)*((216*Dmglst2*(3 - 5*lmMst1 + 5*lmMst2) + 7*(-251 + 54*lmMst1
        - 54*lmMst2)*Mgl)*pow2(Mst2)*pow4(Mst1) - 6*(36*Dmglst2*(5 + 3*lmMst1 -
        3*lmMst2) + (2 - 9*lmMst1 + 9*lmMst2)*Mgl)*pow2(Mst1)*pow4(Mst2) + 4*(
        351*Dmglst2 - 554*Mgl)*pow6(Mst1) - 54*(4*Dmglst2 - 3*Mgl)*pow6(Mst2)))
        + Tbeta*(324*Mst2*s2t*pow2(Mst1)*(160*(Dmglst2*(1 - 9*lmMst1 + 9*
        lmMst2) + (5 + 3*lmMst1 - 3*lmMst2)*Mgl)*pow2(Mst1) - 80*(Dmglst2*(4 +
        6*lmMst1 - 6*lmMst2) + (-7 - 3*lmMst1 + 3*lmMst2)*Mgl)*pow2(Mst2))*
        pow2(MuSUSY)*pow3(Mt) + 6*Mt*pow2(s2t)*(-135*pow2(Mst2)*(-2*(11 + 10*
        lmMst1 - 10*lmMst2)*Mgl*Mst2*s2t*(-2 + pow2(Sbeta)) + (Dmglst2*(372 -
        40*lmMst1 + 40*lmMst2) + 3*(-17 + 2*lmMst1 - 2*lmMst2)*Mgl)*Mt*pow2(
        Sbeta))*pow4(Sbeta)*pow6(Mst1) + 2*Mt*(135*(Dmglst2*(372 - 40*lmMst1 +
        40*lmMst2) + 3*(-17 + 2*lmMst1 - 2*lmMst2)*Mgl)*pow2(Mst2)*pow4(Sbeta)*
        pow6(Mst1) + (5*pow2(MuSUSY)*(648*Dmglst2*pow2(Mst2)*(3*(2 + lmMst1 -
        lmMst2)*pow2(Mst1)*pow2(Mst2) + (3 + 8*lmMst1 - 8*lmMst2)*pow4(Mst1) +
        pow4(Mst2)) + Mgl*(3*(1607 - 432*lmMst1 + 432*lmMst2)*pow2(Mst2)*pow4(
        Mst1) - 18*(25 + 9*lmMst1 - 9*lmMst2)*pow2(Mst1)*pow4(Mst2) + (20701 -
        648*lmMst1 + 648*lmMst2)*pow6(Mst1) - 486*pow6(Mst2))))/3.)) + pow2(
        Sbeta)*(3240*Mt*pow2(Mst1)*pow3(Mst2)*pow3(s2t)*(4*(Dmglst2*(10 - 6*
        lmMst1 + 6*lmMst2) + (1 + 3*lmMst1 - 3*lmMst2)*Mgl)*pow2(Mst1)*pow2(
        Mst2) + (44*Dmglst2 + (-21 - 10*lmMst1 + 10*lmMst2)*Mgl)*pow4(Mst1) -
        4*(Dmglst2*(4 + 6*lmMst1 - 6*lmMst2) + (-7 - 3*lmMst1 + 3*lmMst2)*Mgl)*
        pow4(Mst2)) - 324*Mst2*s2t*pow2(Mst1)*pow3(Mt)*((160*(Dmglst2*(1 - 9*
        lmMst1 + 9*lmMst2) + (5 + 3*lmMst1 - 3*lmMst2)*Mgl)*pow2(Mst1) - 80*(
        Dmglst2*(4 + 6*lmMst1 - 6*lmMst2) + (-7 - 3*lmMst1 + 3*lmMst2)*Mgl)*
        pow2(Mst2))*pow2(MuSUSY) + 160*((Dmglst2*(51 - 14*lmMst1 + 28*lmMst2 -
        14*lmMt) + 4*(-3 + lmMst1 - 2*lmMst2 + lmMt)*Mgl)*pow2(Mst1)*pow2(Mst2)
        + (Dmglst2*(87 - 22*lmMst1 + 44*lmMst2 - 22*lmMt) + (-5 + 2*lmMst1 - 4*
        lmMst2 + 2*lmMt)*Mgl)*pow4(Mst1)) + 80*(2*Dmglst2*(19 - 6*lmMst1 + 12*
        lmMst2 - 6*lmMt) + (-7 + 6*lmMst1 - 12*lmMst2 + 6*lmMt)*Mgl)*pow4(Mst2)
        ) + 12960*pow4(Mt)*(4*Dmglst2*pow2(Mst2)*((23 - 5*lmMst1 + 10*lmMst2 -
        5*lmMt)*pow2(Mst1)*pow2(Mst2) + (77 - 18*lmMst1 + 36*lmMst2 - 18*lmMt)*
        pow4(Mst1) + pow4(Mst2)) + Mgl*(3*(-23 + 6*lmMst1 - 12*lmMst2 + 6*lmMt)
        *pow2(Mst2)*pow4(Mst1) + (-25 + 7*lmMst1 - 14*lmMst2 + 7*lmMt)*pow2(
        Mst1)*pow4(Mst2) + 4*(-7 + 3*lmMst1 - 6*lmMst2 + 3*lmMt)*pow6(Mst1) -
        3*pow6(Mst2))) - 15*pow4(Mst2)*pow4(s2t)*((432*Dmglst2*(4 - lmMst1 +
        lmMst2) + (-1745 + 324*lmMst1 - 324*lmMst2)*Mgl)*pow2(Mst2)*pow4(Mst1)
        - 6*(36*Dmglst2*(4 + 3*lmMst1 - 3*lmMst2) + (29 - 9*lmMst1 + 9*lmMst2)*
        Mgl)*pow2(Mst1)*pow4(Mst2) + 27*(4*Dmglst2*(7 + 10*lmMst1 - 10*lmMst2)
        + (-17 - 14*lmMst1 + 14*lmMst2)*Mgl)*pow6(Mst1) - 54*(4*Dmglst2 - 3*
        Mgl)*pow6(Mst2)) - 20*pow2(Mt)*pow2(s2t)*(3*(1728*Dmglst2 - 115*Mgl)*
        pow4(Mst1)*pow4(Mst2) + (5184*Dmglst2 - 1355*Mgl)*pow2(Mst2)*pow6(Mst1)
        + pow2(MuSUSY)*(648*Dmglst2*pow2(Mst2)*(3*(2 + lmMst1 - lmMst2)*pow2(
        Mst1)*pow2(Mst2) + (3 + 8*lmMst1 - 8*lmMst2)*pow4(Mst1) + pow4(Mst2)) +
        Mgl*(3*(1607 - 432*lmMst1 + 432*lmMst2)*pow2(Mst2)*pow4(Mst1) - 18*(25
        + 9*lmMst1 - 9*lmMst2)*pow2(Mst1)*pow4(Mst2) + (20701 - 648*lmMst1 +
        648*lmMst2)*pow6(Mst1) - 486*pow6(Mst2))) + 18*(288*Dmglst2 + 29*Mgl)*
        pow2(Mst1)*pow6(Mst2) + 324*(4*Dmglst2 - 3*Mgl)*pow8(Mst2))))))))/pow4(
        Msq)))/(Tbeta*pow2(Mst1)*pow2(Sbeta)*pow7(Mst2)) + (Mgl*(68040*Mst2*
        twoLoopFlag*(4*Dmglst2*(4*pow2(Mt)*(-12*(MuSUSY + 2*Mst2*Tbeta)*pow2(
        Mt)*pow2(Sbeta) + Mst2*MuSUSY*pow2(s2t)*(-(MuSUSY*Tbeta*(-1 + pow2(
        Sbeta))) + 6*Mst2*pow2(Sbeta)) + Mt*s2t*Tbeta*(13*pow2(MuSUSY)*(-1 +
        pow2(Sbeta)) + 12*pow2(Mst2)*pow2(Sbeta)))*pow4(Mst1) + pow2(Mst1)*
        pow2(Mst2)*(8*s2t*(3*MuSUSY*s2t + 2*Mt*Tbeta)*pow2(Mst2)*pow2(Mt)*pow2(
        Sbeta) + Tbeta*(-4*Mst2*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*(-1 + pow2(
        Sbeta)) + 36*s2t*pow2(MuSUSY)*(-1 + pow2(Sbeta))*pow3(Mt) + (2*Mt -
        Mst2*s2t)*pow2(Sbeta)*pow3(s2t)*pow4(Mst2) - 32*Mst2*pow2(Sbeta)*pow4(
        Mt))) + pow4(Mst2)*(2*Mst2*MuSUSY*pow2(Mt)*pow2(s2t)*(-2*MuSUSY*Tbeta*(
        -1 + pow2(Sbeta)) + 15*Mst2*pow2(Sbeta)) + 4*s2t*Tbeta*(5*pow2(MuSUSY)*
        (-1 + pow2(Sbeta)) - 4*pow2(Mst2)*pow2(Sbeta))*pow3(Mt) - 2*Mt*(2*
        MuSUSY + 5*Mst2*Tbeta)*pow2(Sbeta)*pow3(Mst2)*pow3(s2t) + 16*MuSUSY*
        pow2(Sbeta)*pow4(Mt) + Tbeta*pow2(Sbeta)*pow4(s2t)*pow5(Mst2))) + Mgl*(
        pow4(Mst2)*(4*Mst2*MuSUSY*pow2(Mt)*pow2(s2t)*(-(MuSUSY*Tbeta*(-1 +
        pow2(Sbeta))) + 6*Mst2*pow2(Sbeta)) + 16*s2t*Tbeta*(pow2(MuSUSY)*(-1 +
        pow2(Sbeta)) - 4*pow2(Mst2)*pow2(Sbeta))*pow3(Mt) - 4*Mt*(MuSUSY + 2*
        Mst2*Tbeta)*pow2(Sbeta)*pow3(Mst2)*pow3(s2t) + 32*(2*MuSUSY + Mst2*
        Tbeta)*pow2(Sbeta)*pow4(Mt) + Tbeta*pow2(Sbeta)*pow4(s2t)*pow5(Mst2)) +
        pow4(Mst1)*(16*s2t*Tbeta*(pow2(MuSUSY)*(-1 + pow2(Sbeta)) - 4*pow2(
        Mst2)*pow2(Sbeta))*pow3(Mt) + 96*(2*MuSUSY + Mst2*Tbeta)*pow2(Sbeta)*
        pow4(Mt) + Tbeta*pow2(Sbeta)*pow4(s2t)*pow5(Mst2)) - 2*pow2(Mst1)*pow2(
        Mst2)*(pow2(Sbeta)*(-2*Mt*(MuSUSY + 2*Mst2*Tbeta)*pow3(Mst2)*pow3(s2t)
        - 32*(2*MuSUSY + Mst2*Tbeta)*pow4(Mt)) + Tbeta*(8*s2t*(pow2(MuSUSY) + (
        4*pow2(Mst2) - pow2(MuSUSY))*pow2(Sbeta))*pow3(Mt) + pow2(Sbeta)*pow4(
        s2t)*pow5(Mst2))))) + (Al4p*threeLoopFlag*(28*Mst2*Mt*MuSUSY*pow2(
        Sbeta)*(9*Mt*pow2(Mst1)*pow2(Mst2)*pow2(s2t)*(pow2(Msq)*(3*Dmglst2*(
        935323 - 279000*lmMst1 + 385560*lmMst2)*pow2(Mst1)*pow2(Mst2) + 15*(
        11075 + 17928*lmMst1 - 17352*lmMst2)*Mgl*pow2(Mst1)*pow2(Mst2) + (10*
        Dmglst2*(1700119 - 484056*lmMst1 + 579960*lmMst2)*pow4(Mst1))/3. + (
        61286 + 434160*lmMst1 - 425520*lmMst2)*Mgl*pow4(Mst1) + 9*Dmglst2*(
        173947 - 25080*lmMst1 + 68760*lmMst2)*pow4(Mst2) + 135*(2269 + 664*
        lmMst1 - 56*lmMst2)*Mgl*pow4(Mst2)) + 64800*Dmsqst2*(Dmglst2*(6*(1 - 2*
        lmMst1 + 2*lmMst2)*pow2(Mst1)*pow2(Mst2) + (17 - 12*lmMst1 + 12*lmMst2)
        *pow4(Mst1) - 2*(2 + 3*lmMst1 - 3*lmMst2)*pow4(Mst2)) + Mgl*(2*(3 + 2*
        lmMst1 - 2*lmMst2)*pow2(Mst1)*pow2(Mst2) + (5 + 4*lmMst1 - 4*lmMst2)*
        pow4(Mst1) + 2*(2 + lmMst1 - lmMst2)*pow4(Mst2)))) + 4*pow2(Mst1)*pow3(
        Mt)*(388800*Dmsqst2*pow2(Mst2)*(10*Dmglst2*(-3 + lmMst1 - 2*lmMst2 +
        lmMt)*pow2(Mst1) - 2*(-1 + lmMst1 - 2*lmMst2 + lmMt)*Mgl*pow2(Mst1) +
        Dmglst2*(-8 + 3*lmMst1 - 6*lmMst2 + 3*lmMt)*pow2(Mst2) - (lmMst1 - 2*
        lmMst2 + lmMt)*Mgl*pow2(Mst2)) + pow2(Msq)*(Dmglst2*(216*(-67843 +
        10050*lmMst1 - 17490*lmMst2 + 7560*lmMt)*pow2(Mst1)*pow2(Mst2) + 4*(-
        13463149 + 1492020*lmMst1 - 2713500*lmMst2 + 625320*lmMt)*pow4(Mst1) +
        405*(-3245 + 1672*lmMst1 - 1448*lmMst2 + 1888*lmMt)*pow4(Mst2)) - 3*
        Mgl*(48*(-41936 + 7245*lmMst1 - 19665*lmMst2 + 720*lmMt)*pow2(Mst1)*
        pow2(Mst2) + 4*(-1077991 + 184140*lmMst1 - 400140*lmMst2 + 8640*lmMt)*
        pow4(Mst1) + 9*(-66773 + 9960*lmMst1 - 45480*lmMst2 + 3840*lmMt)*pow4(
        Mst2)))) + 9720*Mgl*s2t*(shiftst3*pow2(Msq)*pow2(Mst2)*(4*pow2(Mt) + ((
        1 + lmMst1 - lmMst2)*pow2(Mst1) - pow2(Mst2))*pow2(s2t)) + 10*shiftst1*
        (Dmsqst2 + pow2(Msq))*(-4*pow2(Mst1)*pow2(Mt) + 4*pow2(Mst2)*pow2(Mt) +
        2*(lmMst1 - lmMst2)*pow2(Mst1)*pow2(Mst2)*pow2(s2t) + pow2(s2t)*(pow4(
        Mst1) - pow4(Mst2))))*pow5(Mst2) - 6*Mst2*s2t*pow2(Mt)*(pow2(Msq)*(Mgl*
        (3*(623599 + 65160*lmMst1 - 134280*lmMst2)*pow2(Mst2)*pow4(Mst1) + 9*(
        1367 + 13560*lmMst1 - 36600*lmMst2)*pow2(Mst1)*pow4(Mst2) + (5161826 +
        406080*lmMst1 - 613440*lmMst2)*pow6(Mst1) - 41040*pow6(Mst2)) - 24*
        Dmglst2*(6*(19579 + 1770*lmMst1 + 12630*lmMst2)*pow2(Mst2)*pow4(Mst1) +
        6*(4429 - 270*lmMst1 + 8910*lmMst2)*pow2(Mst1)*pow4(Mst2) + (510139 +
        44280*lmMst1 + 76680*lmMst2)*pow6(Mst1) - 2520*pow6(Mst2))) + 100*
        Dmsqst2*(1296*Dmglst2*pow2(Mst2)*(5*pow2(Mst1)*pow2(Mst2) + 9*pow4(
        Mst1) + pow4(Mst2)) + Mgl*(3249*pow2(Mst2)*pow4(Mst1) + 1431*pow2(Mst1)
        *pow4(Mst2) + 4057*pow6(Mst1) - 648*pow6(Mst2)))) + 3*pow3(Mst2)*pow3(
        s2t)*(-12*(150*Dmsqst2*(36*Dmglst2*(3 - 5*lmMst1 + 5*lmMst2) + (-469 +
        36*lmMst1 - 36*lmMst2)*Mgl) + (3*Dmglst2*(32663 - 7620*lmMst1 + 7620*
        lmMst2) + (-224837 + 4680*lmMst1 - 8370*lmMst2)*Mgl)*pow2(Msq))*pow2(
        Mst2)*pow4(Mst1) + 90*(Dmsqst2*(720*Dmglst2*(5 + 3*lmMst1 - 3*lmMst2) +
        1950*Mgl) + (2*Dmglst2*(648*lmMst1 - 7*(257 + 240*lmMst2)) + (5825 +
        564*lmMst1 - 1104*lmMst2)*Mgl)*pow2(Msq))*pow2(Mst1)*pow4(Mst2) + ((3*
        Mgl*(1082800*Dmsqst2 + (3447379 - 76140*lmMst1 + 76140*lmMst2)*pow2(
        Msq)) - 16*Dmglst2*(78975*Dmsqst2 + (788723 - 80595*lmMst1 + 80595*
        lmMst2)*pow2(Msq)))*pow6(Mst1))/3. + 1080*(60*Dmglst2*Dmsqst2 - 30*
        Dmsqst2*Mgl + 28*Dmglst2*pow2(Msq) - 19*Mgl*pow2(Msq))*pow6(Mst2))) +
        Tbeta*(-56*pow2(Mst1)*pow2(MuSUSY)*pow3(Mt)*(388800*Dmsqst2*s2t*(
        Dmglst2*(1 - 9*lmMst1 + 9*lmMst2)*pow2(Mst1) + (5 + 3*lmMst1 - 3*
        lmMst2)*Mgl*pow2(Mst1) - Dmglst2*(2 + 3*lmMst1 - 3*lmMst2)*pow2(Mst2) +
        (2 + lmMst1 - lmMst2)*Mgl*pow2(Mst2))*pow3(Mst2) + pow2(Msq)*(Dmglst2*
        Mst2*(36*Mst2*(30*(3643 - 120*lmMst1 + 3192*lmMst2)*Mt + (364291 -
        88560*lmMst1 + 147960*lmMst2)*Mst2*s2t)*pow2(Mst1) + 27*(7680*(5 + 6*
        lmMst2)*Mt + (173947 - 25080*lmMst1 + 68760*lmMst2)*Mst2*s2t)*pow3(
        Mst2) + 2*(15057833 - 4014360*lmMst1 + 5563080*lmMst2)*s2t*pow4(Mst1))
        - 3*Mgl*(-60*(3*(-353 + 72*lmMst1 + 696*lmMst2)*Mt + 2*(3937 + 2988*
        lmMst1 - 2232*lmMst2)*Mst2*s2t)*pow2(Mst1)*pow2(Mst2) + 2*(15*(38401 +
        1080*lmMst1 - 7992*lmMst2)*Mt + (-266863 - 396360*lmMst1 + 346680*
        lmMst2)*Mst2*s2t)*pow4(Mst1) - 45*(64*(53 + 24*lmMst2)*Mt + 3*(2269 +
        664*lmMst1 - 56*lmMst2)*Mst2*s2t)*pow4(Mst2)))) + pow2(Sbeta)*(7*pow4(
        s2t)*pow6(Mst2)*(pow2(Msq)*(3*Mgl*(6*(-362299 + 17820*lmMst1 - 33300*
        lmMst2)*pow2(Mst2)*pow4(Mst1) - 90*(6053 + 564*lmMst1 - 1104*lmMst2)*
        pow2(Mst1)*pow4(Mst2) + 5*(-149867 + 3996*lmMst1 + 4860*lmMst2)*pow6(
        Mst1) + 20520*pow6(Mst2)) - 4*Dmglst2*(108*(-5917 + 1095*lmMst1 + 195*
        lmMst2)*pow2(Mst2)*pow4(Mst1) + 135*(-1967 + 648*lmMst1 - 1680*lmMst2)*
        pow2(Mst1)*pow4(Mst2) + (-2272991 + 116640*lmMst1 - 116640*lmMst2)*
        pow6(Mst1) + 22680*pow6(Mst2))) + Dmsqst2*(-300*Mgl*(9*(743 - 72*lmMst1
        + 72*lmMst2)*pow2(Mst2)*pow4(Mst1) + 2079*pow2(Mst1)*pow4(Mst2) + (2386
        + 648*lmMst1 - 648*lmMst2)*pow6(Mst1) - 324*pow6(Mst2)) + 97200*
        Dmglst2*(4*(4 - lmMst1 + lmMst2)*pow2(Mst2)*pow4(Mst1) + (7 + 10*lmMst1
        - 10*lmMst2)*pow6(Mst1) - 2*((4 + 3*lmMst1 - 3*lmMst2)*pow2(Mst1)*pow4(
        Mst2) + pow6(Mst2))))) + 56*Mst2*s2t*pow2(Mst1)*pow3(Mt)*(388800*
        Dmsqst2*pow2(Mst2)*(Dmglst2*((-2 - 3*lmMst1 + 3*lmMst2)*pow2(Mst2)*
        pow2(MuSUSY) + pow2(Mst1)*((51 - 14*lmMst1 + 28*lmMst2 - 14*lmMt)*pow2(
        Mst2) + (1 - 9*lmMst1 + 9*lmMst2)*pow2(MuSUSY)) + (87 - 22*lmMst1 + 44*
        lmMst2 - 22*lmMt)*pow4(Mst1) + (19 - 6*lmMst1 + 12*lmMst2 - 6*lmMt)*
        pow4(Mst2)) + Mgl*((2 + lmMst1 - lmMst2)*pow2(Mst2)*pow2(MuSUSY) +
        pow2(Mst1)*((-5 + 2*lmMst1 - 4*lmMst2 + 2*lmMt)*pow2(Mst2) + (5 + 3*
        lmMst1 - 3*lmMst2)*pow2(MuSUSY)) + (-5 + 2*lmMst1 - 4*lmMst2 + 2*lmMt)*
        pow4(Mst1) + (-1 + 2*lmMst1 - 4*lmMst2 + 2*lmMt)*pow4(Mst2))) + pow2(
        Msq)*(Dmglst2*(2*(4*(9908167 - 949320*lmMst1 + 1769040*lmMst2 - 217080*
        lmMt)*pow2(Mst2) + (15057833 - 4014360*lmMst1 + 5563080*lmMst2)*pow2(
        MuSUSY))*pow4(Mst1) + 18*pow2(Mst1)*(2*(364291 - 88560*lmMst1 + 147960*
        lmMst2)*pow2(Mst2)*pow2(MuSUSY) + 3*(510149 - 55320*lmMst1 + 118200*
        lmMst2 - 32160*lmMt)*pow4(Mst2)) + 27*((173947 - 25080*lmMst1 + 68760*
        lmMst2)*pow2(MuSUSY)*pow4(Mst2) + 30*(-1672*lmMst1 + 1448*lmMst2 + 59*(
        71 - 32*lmMt))*pow6(Mst2))) - 3*Mgl*(2*((2299036 - 388800*lmMst1 +
        656640*lmMst2)*pow2(Mst2) + (-266863 - 396360*lmMst1 + 346680*lmMst2)*
        pow2(MuSUSY))*pow4(Mst1) + 6*pow2(Mst1)*(-20*(3937 + 2988*lmMst1 -
        2232*lmMst2)*pow2(Mst2)*pow2(MuSUSY) + (470657 - 86040*lmMst1 + 178200*
        lmMst2)*pow4(Mst2)) + 9*(-15*(2269 + 664*lmMst1 - 56*lmMst2)*pow2(
        MuSUSY)*pow4(Mst2) + 2*(68693 - 9960*lmMst1 + 45480*lmMst2 - 3840*lmMt)
        *pow6(Mst2))))) + 28*Mt*pow2(Mst2)*pow2(s2t)*(s2t*pow2(Mst1)*pow3(Mst2)
        *(3*Mgl*pow2(Msq)*(30*(4673 - 5976*lmMst1 + 8424*lmMst2)*pow2(Mst1)*
        pow2(Mst2) + 17*(6167 - 9720*lmMst1 + 9720*lmMst2)*pow4(Mst1) - 135*(
        2269 + 664*lmMst1 - 56*lmMst2)*pow4(Mst2)) - Dmglst2*pow2(Msq)*(18*(
        206741 - 101880*lmMst1 + 89640*lmMst2)*pow2(Mst1)*pow2(Mst2) + (8583283
        - 2329560*lmMst1 + 2329560*lmMst2)*pow4(Mst1) + 27*(173947 - 25080*
        lmMst1 + 68760*lmMst2)*pow4(Mst2)) + 194400*Dmsqst2*(Dmglst2*(2*(-5 +
        3*lmMst1 - 3*lmMst2)*pow2(Mst1)*pow2(Mst2) - 11*pow4(Mst1) + 2*(2 + 3*
        lmMst1 - 3*lmMst2)*pow4(Mst2)) + Mgl*(-2*(1 + lmMst1 - lmMst2)*pow2(
        Mst1)*pow2(Mst2) + pow4(Mst1) - 2*(2 + lmMst1 - lmMst2)*pow4(Mst2)))) +
        Mt*(pow2(Msq)*(-4*Dmglst2*(162*pow2(Mst2)*(20*(505 + 68*lmMst1 + 124*
        lmMst2)*pow2(Mst2) + (6803 - 1810*lmMst1 + 2670*lmMst2)*pow2(MuSUSY))*
        pow4(Mst1) - 27*pow2(Mst1)*(4*(-4849 + 270*lmMst1 - 8910*lmMst2)*pow2(
        Mst2) + 5*(-1631 + 648*lmMst1 - 1680*lmMst2)*pow2(MuSUSY))*pow4(Mst2) +
        2*(45*(78533 + 6732*lmMst1 + 180*lmMst2)*pow2(Mst2) + (2128489 -
        307800*lmMst1 + 377460*lmMst2)*pow2(MuSUSY))*pow6(Mst1) - 22680*(2*
        pow2(Mst2) + pow2(MuSUSY))*pow6(Mst2)) + 3*Mgl*(6*pow2(Mst2)*((309749 +
        12240*lmMst1 - 12240*lmMst2)*pow2(Mst2) + (533629 - 900*lmMst1 + 180*
        lmMst2)*pow2(MuSUSY))*pow4(Mst1) + 9*pow2(Mst1)*((5927 + 13560*lmMst1 -
        36600*lmMst2)*pow2(Mst2) + 10*(5597 + 564*lmMst1 - 1104*lmMst2)*pow2(
        MuSUSY))*pow4(Mst2) + ((3291029 + 210600*lmMst1 - 210600*lmMst2)*pow2(
        Mst2) + (6649153 - 81540*lmMst1 + 77220*lmMst2)*pow2(MuSUSY))*pow6(
        Mst1) - 20520*(2*pow2(Mst2) + pow2(MuSUSY))*pow6(Mst2))) + 300*Dmsqst2*
        (648*Dmglst2*pow2(Mst2)*(pow2(Mst1)*pow2(Mst2)*(8*pow2(Mst2) + 3*(2 +
        lmMst1 - lmMst2)*pow2(MuSUSY)) + (8*pow2(Mst2) + (3 + 8*lmMst1 - 8*
        lmMst2)*pow2(MuSUSY))*pow4(Mst1) + (2*pow2(Mst2) + pow2(MuSUSY))*pow4(
        Mst2) + 8*pow6(Mst1)) + Mgl*(9*pow2(Mst2)*(202*pow2(Mst2) + (1097 - 72*
        lmMst1 + 72*lmMst2)*pow2(MuSUSY))*pow4(Mst1) + (808*pow2(Mst2) + (20701
        - 648*lmMst1 + 648*lmMst2)*pow2(MuSUSY))*pow6(Mst1) - 324*(2*pow2(Mst2)
        + pow2(MuSUSY))*pow6(Mst2) + 27*pow2(Mst1)*(53*pow2(MuSUSY)*pow4(Mst2)
        + 77*pow6(Mst2)))))) - 16*pow4(Mt)*(680400*Dmsqst2*(2*Dmglst2*pow4(
        Mst2)*((23 - 5*lmMst1 + 10*lmMst2 - 5*lmMt)*pow2(Mst1)*pow2(Mst2) + (77
        - 18*lmMst1 + 36*lmMst2 - 18*lmMt)*pow4(Mst1) + pow4(Mst2)) - Mgl*(-2*(
        -3 + lmMst1 - 2*lmMst2 + lmMt)*pow2(Mst1)*(2*pow2(Mst1) + pow2(Mst2))*
        pow4(Mst2) + 2*(7 - 3*lmMst1 + 6*lmMst2 - 3*lmMt)*pow2(Mst2)*pow6(Mst1)
        + pow8(Mst2))) - pow2(Msq)*(-4*Dmglst2*pow2(Mst2)*(27*pow2(Mst1)*pow2(
        Mst2)*((209341 - 26880*lmMst1 + 82320*lmMst2 - 1680*lmMt)*pow2(Mst2) -
        6720*(5 + 6*lmMst2)*pow2(MuSUSY)) + 9*((3184397 - 476280*lmMst1 +
        1156680*lmMst2 - 161280*lmMt)*pow2(Mst2) + 105*(-3643 + 120*lmMst1 -
        3192*lmMst2)*pow2(MuSUSY))*pow4(Mst1) + 4*(19152737 - 3674160*lmMst1 +
        6872040*lmMst2 - 850500*lmMt)*pow6(Mst1) + 158760*pow6(Mst2)) + 21*Mgl*
        (3*pow2(Mst2)*((310039 - 91080*lmMst1 + 145080*lmMst2 - 43920*lmMt)*
        pow2(Mst2) + 30*(-353 + 72*lmMst1 + 696*lmMst2)*pow2(MuSUSY))*pow4(
        Mst1) + 9*pow2(Mst1)*((48983 - 12480*lmMst1 + 26820*lmMst2 - 12180*
        lmMt)*pow2(Mst2) + 160*(53 + 24*lmMst2)*pow2(MuSUSY))*pow4(Mst2) + (2*(
        626869 - 300240*lmMst1 + 395280*lmMst2 - 78840*lmMt)*pow2(Mst2) - 15*(
        38401 + 1080*lmMst1 - 7992*lmMst2)*pow2(MuSUSY))*pow6(Mst1) + 20520*
        pow8(Mst2))))) + pow2(Mst2)*(-28*pow2(Mt)*pow2(s2t)*(pow2(Msq)*pow2(
        MuSUSY)*(-4*Dmglst2*(162*(6803 - 1810*lmMst1 + 2670*lmMst2)*pow2(Mst2)*
        pow4(Mst1) + 135*(1631 - 648*lmMst1 + 1680*lmMst2)*pow2(Mst1)*pow4(
        Mst2) + (4256978 - 615600*lmMst1 + 754920*lmMst2)*pow6(Mst1) - 22680*
        pow6(Mst2)) + 3*Mgl*(6*(533629 - 900*lmMst1 + 180*lmMst2)*pow2(Mst2)*
        pow4(Mst1) + 90*(5597 + 564*lmMst1 - 1104*lmMst2)*pow2(Mst1)*pow4(Mst2)
        + (6649153 - 81540*lmMst1 + 77220*lmMst2)*pow6(Mst1) - 20520*pow6(Mst2)
        )) + 300*Dmsqst2*(162*Dmglst2*pow2(Mst2)*(4*pow2(MuSUSY)*(3*(2 + lmMst1
        - lmMst2)*pow2(Mst1)*pow2(Mst2) + (3 + 8*lmMst1 - 8*lmMst2)*pow4(Mst1)
        + pow4(Mst2)) + (-93 + 10*lmMst1 - 10*lmMst2)*(-2 + pow2(Sbeta))*pow4(
        Sbeta)*pow6(Mst1)) + Mgl*pow2(MuSUSY)*(9*(1097 - 72*lmMst1 + 72*lmMst2)
        *pow2(Mst2)*pow4(Mst1) + 1431*pow2(Mst1)*pow4(Mst2) + (20701 - 648*
        lmMst1 + 648*lmMst2)*pow6(Mst1) - 324*pow6(Mst2)))) + 68040*Mgl*(
        shiftst3*pow2(Msq)*(-4*pow2(Mst1)*pow2(s2t)*pow4(Mst2)*(2*(-1 + lmMst1
        - lmMst2)*pow2(Mt)*(-(pow2(MuSUSY)*(-1 + pow2(Sbeta))) + pow2(Mst2)*
        pow2(Sbeta)) + pow2(s2t)*pow2(Sbeta)*pow4(Mst2)) + 8*(-1 + 3*lmMst1 -
        3*lmMst2)*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta))*pow6(Mst1)
        + (-4*pow2(Mt)*pow2(s2t)*(pow2(MuSUSY)*(-1 + pow2(Sbeta)) + 2*pow2(
        Mst2)*pow2(Sbeta)) + pow2(Sbeta)*(16*pow4(Mt) + pow4(Mst2)*pow4(s2t)))*
        pow6(Mst2) + pow4(Mst1)*(8*(-1 + 2*lmMst1 - 2*lmMst2)*pow2(Mst2)*pow2(
        Mt)*pow2(MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta)) + (3 + 2*lmMst1 - 2*
        lmMst2)*pow2(Sbeta)*pow4(s2t)*pow6(Mst2))) + 10*shiftst1*(Dmsqst2 +
        pow2(Msq))*(4*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*(pow2(Mst1)*pow4(Mst2) -
        2*(lmMst1 - lmMst2)*pow2(Mst1)*(pow2(Mst1)*pow2(Mst2) + pow4(Mst1) +
        pow4(Mst2)) + pow6(Mst2)) + pow2(Sbeta)*(16*(pow2(Mst1) + pow2(Mst2))*
        pow4(Mst2)*pow4(Mt) + (pow2(Mst1) - pow2(Mst2))*(2*(lmMst1 - lmMst2)*
        pow2(Mst1)*pow2(Mst2) + pow4(Mst1) - pow4(Mst2))*pow4(Mst2)*pow4(s2t) -
        4*pow2(Mt)*pow2(s2t)*(2*pow4(Mst1)*pow4(Mst2) - 4*pow2(Mst1)*pow6(Mst2)
        + pow2(MuSUSY)*(pow2(Mst1)*pow4(Mst2) - 2*(lmMst1 - lmMst2)*pow2(Mst1)*
        (pow2(Mst1)*pow2(Mst2) + pow4(Mst1) + pow4(Mst2)) + pow6(Mst2)) + 2*
        pow8(Mst2)))))))))/(pow2(Msq)*pow2(Mst1))))/(Tbeta*pow2(Sbeta)*pow6(
        Mst2))))/408240. + (-(xMst*(Mgl*(27*(lmMst1 - lmMst2)*Mgl*Mst2*
        oneLoopFlag*Tbeta*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*(-1 + pow2(Sbeta)) +
        2*Al4p*Mt*MuSUSY*twoLoopFlag*pow2(Sbeta)*(8*Dmglst2*(36*(-1 + 3*lmMst1
        - 3*lmMst2)*Mst2*s2t*pow2(Mt) - 3*Mt*(50 + lmMst1*(51 - 18*lmMst2) -
        51*lmMst2 + 18*pow2(lmMst2))*pow2(Mst2)*pow2(s2t) + 4*(7 - 381*lmMst2 +
        lmMst1*(327 + 72*lmMst2 - 81*lmMt) + 54*lmMt + 81*lmMst2*lmMt - 72*
        pow2(lmMst2))*pow3(Mt) + 2*(1 + 3*lmMst1 - 3*lmMst2)*pow3(Mst2)*pow3(
        s2t)) - Mgl*(-144*Mst2*s2t*(-lmMst1 + lmMst2 - 2*lmMst1*lmMst2 + pow2(
        lmMst1) + pow2(lmMst2))*pow2(Mt) + 24*Mt*(1 + 33*lmMst2 - 3*lmMst1*(11
        + 6*lmMst2) + 18*pow2(lmMst2))*pow2(Mst2)*pow2(s2t) + 32*(83 - 96*
        lmMst2 + 3*lmMst1*(29 + 12*lmMst2 - 9*lmMt) + 9*lmMt + 27*lmMst2*lmMt -
        36*pow2(lmMst2))*pow3(Mt) + (5 + lmMst1*(6 - 288*lmMst2) - 6*lmMst2 +
        144*pow2(lmMst1) + 144*pow2(lmMst2))*pow3(Mst2)*pow3(s2t))) + Al4p*s2t*
        Tbeta*twoLoopFlag*(8*pow2(Mt)*((Dmglst2*Mst2*s2t*(-49 + 84*lmMst2 - 12*
        lmMst1*(7 + 3*lmMst2) + 36*pow2(lmMst2)) + Mgl*Mst2*s2t*(-1 + 111*
        lmMst2 - 3*lmMst1*(37 + 72*lmMst2) + 81*pow2(lmMst1) + 135*pow2(lmMst2)
        ) + Dmglst2*Mt*(785 - 438*lmMst2 - 6*lmMst1*(-85 + 42*lmMst2) + 252*
        pow2(lmMst2)) + Mgl*Mt*(193 + 474*lmMst2 - 6*lmMst1*(67 + 42*lmMst2) +
        252*pow2(lmMst2)))*pow2(MuSUSY) - Mt*((4*Dmglst2*(95 - 447*lmMst2 +
        lmMst1*(411 + 90*lmMst2 - 90*lmMt) + 36*lmMt + 90*lmMst2*lmMt - 90*
        pow2(lmMst2)) - 4*Mgl*(49 - 75*lmMst2 + 3*lmMst1*(25 + 6*lmMst2 - 6*
        lmMt) + 18*lmMst2*lmMt - 18*pow2(lmMst2)))*pow2(Mst2) + (Dmglst2*(785 +
        6*lmMst1*(85 - 42*lmMst2) - 438*lmMst2 + 252*pow2(lmMst2)) + Mgl*(193 +
        474*lmMst2 - 6*lmMst1*(67 + 42*lmMst2) + 252*pow2(lmMst2)))*pow2(
        MuSUSY))*pow2(Sbeta)) + Mst2*s2t*pow2(Sbeta)*(8*pow2(Mt)*(Mgl*(1 - 111*
        lmMst2 + 3*lmMst1*(37 + 72*lmMst2) - 81*pow2(lmMst1) - 135*pow2(lmMst2)
        )*pow2(MuSUSY) - Dmglst2*(36*(lmMst1 - lmMst2)*pow2(Mst2) + (-49 + 84*
        lmMst2 - 12*lmMst1*(7 + 3*lmMst2) + 36*pow2(lmMst2))*pow2(MuSUSY))) +
        4*(Dmglst2*(11 + 42*lmMst1 - 42*lmMst2) + (-5 - 6*lmMst1 + 6*lmMst2)*
        Mgl)*Mt*s2t*pow3(Mst2) + (2*Dmglst2*(5 + 6*lmMst1 - 6*lmMst2) + (-11 -
        6*lmMst1 + 6*lmMst2)*Mgl)*pow2(s2t)*pow4(Mst2))) + 32*Al4p*Mst2*Tbeta*
        twoLoopFlag*(-(Mgl*(65 - 159*lmMst2 + 6*lmMst1*(25 + 6*lmMst2 - 6*lmMt)
        + 9*lmMt + 36*lmMst2*lmMt - 36*pow2(lmMst2))) + 18*Dmglst2*(4 - 61*
        lmMst2 + 6*lmMst1*(9 + 2*lmMst2 - 2*lmMt) + 7*lmMt + 12*lmMst2*lmMt -
        12*pow2(lmMst2)))*pow2(Sbeta)*pow4(Mt))*power10(Mst1) + Al4p*
        twoLoopFlag*(-144*(lmMst1 - lmMst2)*Mgl*(2*Dmglst2*lmMst2 + Mgl +
        lmMst2*Mgl)*Mst2*Tbeta*xDR2DRMOD*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*(-1 +
        pow2(Sbeta)) + xDmglst2*pow2(Dmglst2)*(8*s2t*(Tbeta*(2501 + 42*lmMst1*(
        7 - 6*lmMst2) - 222*lmMst2 + 252*pow2(lmMst2))*pow2(MuSUSY)*(1 - pow2(
        Sbeta)) + Mst2*pow2(Sbeta)*(-4*Mst2*Tbeta*(257 + 6*lmMst2*(179 - 30*
        lmMt) - 12*lmMst1*(76 + 15*lmMst2 - 15*lmMt) - 162*lmMt + 180*pow2(
        lmMst2)) + 3*MuSUSY*(132 - 215*pow2(Sbeta) + 60*(lmMst1 - lmMst2)*(-3 +
        5*pow2(Sbeta)))))*pow3(Mt) - 6*Mst2*pow2(Mt)*pow2(s2t)*(4*Mst2*MuSUSY*(
        355 + 6*lmMst2 - 6*lmMst1*(1 + 6*lmMst2) + 36*pow2(lmMst2))*pow2(Sbeta)
        + Tbeta*(-2*(85 + 60*lmMst2 + 12*lmMst1*(-5 + 3*lmMst2) - 36*pow2(
        lmMst2))*pow2(MuSUSY)*(-1 + pow2(Sbeta)) + pow2(Mst2)*pow2(Sbeta)*(96 -
        430*pow2(Sbeta) + 215*pow4(Sbeta) - 12*(lmMst1 - lmMst2)*(18 - 50*pow2(
        Sbeta) + 25*pow4(Sbeta))))) + pow2(Sbeta)*(4*Mt*(12*(1 - 6*lmMst1 + 6*
        lmMst2)*MuSUSY + (35 - 102*lmMst1 + 102*lmMst2)*Mst2*Tbeta)*pow3(Mst2)*
        pow3(s2t) + 32*(9*Mst2*Tbeta*(102 + lmMst2*(353 - 60*lmMt) - 6*lmMst1*(
        49 + 10*lmMst2 - 10*lmMt) - 59*lmMt + 60*pow2(lmMst2)) + MuSUSY*(412 -
        6*lmMst1*(176 + 42*lmMst2 - 39*lmMt) + 9*lmMst2*(147 - 26*lmMt) - 267*
        lmMt + 252*pow2(lmMst2)))*pow4(Mt) + 3*(5 + 6*lmMst1 - 6*lmMst2)*Tbeta*
        pow4(s2t)*pow5(Mst2))))*power10(Mst1) + 48*Al4p*Mst2*Tbeta*xDmglst2*
        pow2(Dmglst2)*pow2(s2t)*(2*Al4p*threeLoopFlag*pow14(Mst2)*(2 + lmMst2 -
        3*pow2(lmMst2))*pow2(s2t)*pow2(Sbeta) - 2*Al4p*threeLoopFlag*xDR2DRMOD*
        pow14(Mst2)*(2 + lmMst2 - 3*pow2(lmMst2))*pow2(s2t)*pow2(Sbeta) + 3*(-
        lmMst1 + lmMst2)*(-2 + 3*lmMst2)*twoLoopFlag*xDR2DRMOD*pow2(Mt)*pow2(
        MuSUSY)*(-1 + pow2(Sbeta))*power10(Mst1))))/(108.*pow7(Mst2)) - (
        threeLoopFlag*pow2(Al4p)*(Mst2*xDmsqst2*pow2(Dmsqst2)*pow2(Mst1)*(
        2058000*Tbeta*pow2(Mgl)*(2*pow2(s2t)*(-(T1ep*pow2(Mst1)*(3*pow2(s2t)*
        pow2(Sbeta)*pow4(Mst2)*(13*pow2(Mst1)*pow2(Mst2) + 6*pow4(Mst2)) + 4*
        pow2(Mt)*pow2(MuSUSY)*(75*pow2(Mst1)*pow2(Mst2) + 221*pow4(Mst1) + 18*
        pow4(Mst2)))) - 2*pow2(Mt)*(-2*T1ep*pow2(Mst1)*pow2(Sbeta)*((-7*pow2(
        Mst2) + 221*pow2(MuSUSY))*pow4(Mst1) + 18*(pow2(Mst2) + pow2(MuSUSY))*
        pow4(Mst2) + 3*pow2(Mst1)*(25*pow2(Mst2)*pow2(MuSUSY) + pow4(Mst2))) +
        81*shiftst1*pow2(MuSUSY)*(2*lmMst1*(-2 + lmMst2)*pow2(Mst1)*pow2(Mst2)*
        (pow2(Mst1) + pow2(Mst2)) + 2*(pow2(Mst1) + pow2(Mst2))*pow4(Mst2) - 2*
        pow2(lmMst2)*pow2(Mst1)*(pow2(Mst1)*pow2(Mst2) + pow4(Mst1) + pow4(
        Mst2)) + lmMst1*(-3 + 2*lmMst2)*pow6(Mst1) + lmMst2*(4*pow2(Mst2)*pow4(
        Mst1) + 3*pow2(Mst1)*pow4(Mst2) + 3*pow6(Mst1) - pow6(Mst2))))) + 81*
        shiftst1*pow2(Sbeta)*(pow2(s2t)*((2 - lmMst2)*pow2(Mst2)*pow4(Mst1)*(8*
        pow2(Mt)*(pow2(Mst2) + (-lmMst1 + lmMst2)*pow2(MuSUSY)) + (1 - 2*lmMst1
        + 2*lmMst2)*pow2(s2t)*pow4(Mst2)) + (4*(lmMst1 - lmMst2)*(-3 + 2*
        lmMst2)*pow2(Mt)*pow2(MuSUSY) + (-2 + lmMst2)*pow2(s2t)*pow4(Mst2))*
        pow6(Mst1)) - (2 - lmMst2)*(pow2(Mst1)*pow4(Mst2)*(4*pow2(Mt)*(4*pow2(
        Mst2) + (-1 + 2*lmMst1 - 2*lmMst2)*pow2(MuSUSY))*pow2(s2t) + 16*pow4(
        Mt) + (-1 - 2*lmMst1 + 2*lmMst2)*pow4(Mst2)*pow4(s2t)) + (-4*pow2(Mt)*(
        2*pow2(Mst2) + pow2(MuSUSY))*pow2(s2t) + 16*pow4(Mt) + pow4(Mst2)*pow4(
        s2t))*pow6(Mst2)))) + Mgl*(4*Mt*MuSUSY*pow2(Sbeta)*(37044000*Mst2*pow2(
        Mst1)*(4*Dmglst2*((497 - 306*lmMst1 + 609*lmMst2 - 303*lmMt)*pow2(Mst1)
        + (101 - 90*lmMst1 + 177*lmMst2 - 87*lmMt)*pow2(Mst2)) + 9*Mgl*(4*(1 +
        6*lmMst1 - 13*lmMst2 + 7*lmMt)*pow2(Mst1) + (16 + 16*lmMst1 - 35*lmMst2
        + 19*lmMt)*pow2(Mst2)))*pow3(Mt) - 500094000*Mst2*Mt*pow2(Mst1)*pow2(
        s2t)*(40*Dmglst2*(1 - lmMst1 + lmMst2)*pow2(Mst1)*pow2(Mst2) + Dmglst2*
        (80 - 43*lmMst1 + 43*lmMst2)*pow4(Mst1) + (-4 + 9*lmMst1 - 9*lmMst2)*
        Mgl*pow4(Mst1) - 4*Dmglst2*(2 + 5*lmMst1 - 5*lmMst2)*pow4(Mst2) + 2*
        Mgl*((1 + 8*lmMst1 - 8*lmMst2)*pow2(Mst1)*pow2(Mst2) + (5 + 4*lmMst1 -
        4*lmMst2)*pow4(Mst2))) + 3*pow2(Mst2)*pow3(s2t)*(3087000*Dmglst2*(9*(55
        - 34*lmMst1 + 34*lmMst2)*pow2(Mst2)*pow4(Mst1) - 36*(5 + 2*lmMst1 - 2*
        lmMst2)*pow2(Mst1)*pow4(Mst2) + (419 - 57*lmMst1 + 57*lmMst2)*pow6(
        Mst1) - 108*pow6(Mst2)) + Mgl*(1715*(192439 - 30400*OepS2 - 837000*S2 -
        180*lmMst2*(857 + 3420*S2) + 180*lmMst1*(677 + 180*lmMst2 + 3420*S2) -
        16200*(pow2(lmMst1) + pow2(lmMst2)))*pow2(Mst2)*pow4(Mst1) + 128625*(
        8555 - 128*OepS2 - 57024*S2 - 72*lmMst2*(29 + 36*S2) + 72*lmMst1*(29 -
        12*lmMst2 + 36*S2) + 864*pow2(lmMst2))*pow2(Mst1)*pow4(Mst2) - 2*(
        55313478 + 26068000*OepS2 + 563377500*S2 + 315*lmMst2*(-423553 +
        1675800*S2) - 315*lmMst1*(-423553 + 166320*lmMst2 + 1675800*S2) +
        26195400*(pow2(lmMst1) + pow2(lmMst2)))*pow6(Mst1) + 27783000*(1 + 2*
        lmMst2)*pow6(Mst2))) + 4*s2t*pow2(Mt)*(27783000*Dmglst2*pow2(Mst2)*(82*
        pow2(Mst1)*pow2(Mst2) + (195 - 6*lmMst1 + 6*lmMst2)*pow4(Mst1) + 36*
        pow4(Mst2)) - Mgl*(5145*(279089 - 5600*OepS2 + 1260*lmMst2*(29 - 90*S2)
        + 812700*S2 + 180*lmMst1*(-203 + 90*lmMst2 + 630*S2) - 8100*(pow2(
        lmMst1) + pow2(lmMst2)))*pow2(Mst2)*pow4(Mst1) + 385875*(3655 - 64*
        OepS2 + 8424*S2 + 144*lmMst1*(1 + 9*S2) - 144*lmMst2*(4 + 9*S2))*pow2(
        Mst1)*pow4(Mst2) + (1094369501 - 72716000*OepS2 + 2520*lmMst2*(235453 -
        584325*S2) + 5940931500*S2 + 2520*lmMst1*(-235453 + 114345*lmMst2 +
        584325*S2) - 144074700*(pow2(lmMst1) + pow2(lmMst2)))*pow6(Mst1) +
        83349000*(1 + 2*lmMst2)*pow6(Mst2))) + 2058000*Mgl*s2t*(-81*(2 -
        lmMst2)*shiftst1*(-4*pow2(Mst1)*pow2(Mt) + 4*pow2(Mst2)*pow2(Mt) + 2*(
        lmMst1 - lmMst2)*pow2(Mst1)*pow2(Mst2)*pow2(s2t) + pow2(s2t)*(pow4(
        Mst1) - pow4(Mst2)))*pow4(Mst2) - 2*T1ep*pow2(Mst1)*((106*pow2(Mt) -
        57*pow2(Mst2)*pow2(s2t))*pow4(Mst1) + 36*pow2(Mt)*pow4(Mst2) + pow2(
        Mst1)*(42*pow2(Mst2)*pow2(Mt) - 57*pow2(s2t)*pow4(Mst2)) - 18*pow2(s2t)
        *pow6(Mst2)))) + Tbeta*(2667168000*Mst2*s2t*pow2(Mst1)*(2*Dmglst2*(8 -
        15*lmMst1 + 15*lmMst2)*pow2(Mst1) + (2 + 6*lmMst1 - 6*lmMst2)*Mgl*pow2(
        Mst1) - 2*Dmglst2*(2 + 5*lmMst1 - 5*lmMst2)*pow2(Mst2) + (5 + 4*lmMst1
        - 4*lmMst2)*Mgl*pow2(Mst2))*pow2(MuSUSY)*pow3(Mt) + 20*Mt*pow2(s2t)*(
        1389150*(Dmglst2*(4324 - 660*lmMst1 + 660*lmMst2) + 3*(-167 + 5*lmMst1
        - 5*lmMst2)*Mgl)*Mt*pow2(Mst2)*pow4(Sbeta)*pow6(Mst1) - 694575*pow2(
        Mst2)*(-6*(11 + 39*lmMst1 - 39*lmMst2)*Mgl*Mst2*s2t*(-2 + pow2(Sbeta))
        + (Dmglst2*(4324 - 660*lmMst1 + 660*lmMst2) + 3*(-167 + 5*lmMst1 - 5*
        lmMst2)*Mgl)*Mt*pow2(Sbeta))*pow4(Sbeta)*pow6(Mst1) - Mt*pow2(MuSUSY)*(
        16669800*Dmglst2*pow2(Mst2)*(-8*(4 + lmMst1 - lmMst2)*pow2(Mst1)*pow2(
        Mst2) + (23 - 42*lmMst1 + 42*lmMst2)*pow4(Mst1) - 12*pow4(Mst2)) + Mgl*
        (-4116*(45*lmMst1*(-1547 + 180*lmMst2 - 4500*S2) + 45*lmMst2*(1547 +
        4500*S2) + 2*(-106283 + 5000*OepS2 + 639225*S2) + 4050*pow2(lmMst1) -
        12150*pow2(lmMst2))*pow2(Mst2)*pow4(Mst1) + 77175*(8771 - 128*OepS2 -
        57024*S2 - 72*lmMst2*(23 + 36*S2) + 72*lmMst1*(29 - 12*lmMst2 + 36*S2)
        + 864*pow2(lmMst2))*pow2(Mst1)*pow4(Mst2) + 2*(593331163 - 60642400*
        OepS2 + 1260*lmMst2*(8263 - 974610*S2) + 1143733500*S2 + 1260*lmMst1*(-
        8263 + 71820*lmMst2 + 974610*S2) - 61916400*pow2(lmMst1) - 28576800*
        pow2(lmMst2))*pow6(Mst1) + 16669800*(1 + 2*lmMst2)*pow6(Mst2)))) +
        pow2(Sbeta)*(166698000*Mt*pow2(Mst1)*pow3(Mst2)*pow3(s2t)*(4*Dmglst2*(
        4*(12 - 5*lmMst1 + 5*lmMst2)*pow2(Mst1)*pow2(Mst2) + (40 - 3*lmMst1 +
        3*lmMst2)*pow4(Mst1) - 4*(2 + 5*lmMst1 - 5*lmMst2)*pow4(Mst2)) - Mgl*(
        32*(1 - lmMst1 + lmMst2)*pow2(Mst1)*pow2(Mst2) + (45 + 37*lmMst1 - 37*
        lmMst2)*pow4(Mst1) - 8*(5 + 4*lmMst1 - 4*lmMst2)*pow4(Mst2))) -
        74088000*Mst2*s2t*pow2(Mst1)*pow3(Mt)*(9*Mgl*(4*(5 + 4*lmMst1 - 4*
        lmMst2)*pow2(Mst2)*pow2(MuSUSY) + 4*pow2(Mst1)*((-23 + 14*lmMst1 - 28*
        lmMst2 + 14*lmMt)*pow2(Mst2) + 2*(1 + 3*lmMst1 - 3*lmMst2)*pow2(MuSUSY)
        ) + 2*(-35 - 50*lmMst2 + 4*lmMst1*(8 + lmMst2 - lmMt) + 18*lmMt + 4*
        lmMst2*lmMt - 4*pow2(lmMst2))*pow4(Mst1) + (13 + 32*lmMst1 - 70*lmMst2
        + 38*lmMt)*pow4(Mst2)) - 2*Dmglst2*(36*(2 + 5*lmMst1 - 5*lmMst2)*pow2(
        Mst2)*pow2(MuSUSY) + 36*pow2(Mst1)*(8*(-7 + 3*lmMst1 - 6*lmMst2 + 3*
        lmMt)*pow2(Mst2) + (-8 + 15*lmMst1 - 15*lmMst2)*pow2(MuSUSY)) + 9*(-465
        - 374*lmMst2 + 4*lmMst1*(52 + 3*lmMst2 - 3*lmMt) + 166*lmMt + 12*
        lmMst2*lmMt - 12*pow2(lmMst2))*pow4(Mst1) + (-578 + 360*lmMst1 - 708*
        lmMst2 + 348*lmMt)*pow4(Mst2))) + 3*pow4(Mst2)*pow4(s2t)*(3087000*
        Dmglst2*(9*(-75 + 26*lmMst1 - 26*lmMst2)*pow2(Mst2)*pow4(Mst1) + 72*(1
        + lmMst1 - lmMst2)*pow2(Mst1)*pow4(Mst2) + (76 - 249*lmMst1 + 249*
        lmMst2)*pow6(Mst1) + 108*pow6(Mst2)) + Mgl*(3430*(224593 + 10400*OepS2
        - 1719900*S2 + 1170*lmMst2*(-1 + 180*S2) - 90*lmMst1*(-193 + 540*lmMst2
        + 2340*S2) + 8100*pow2(lmMst1) + 40500*pow2(lmMst2))*pow2(Mst2)*pow4(
        Mst1) - 128625*(8339 - 128*OepS2 - 57024*S2 - 72*lmMst2*(35 + 36*S2) +
        72*lmMst1*(29 - 12*lmMst2 + 36*S2) + 864*pow2(lmMst2))*pow2(Mst1)*pow4(
        Mst2) + (440659841 + 1890*lmMst1*(251761 - 26040*lmMst2) - 531394290*
        lmMst2 - 308700000*S2 + 24607800*pow2(lmMst1) + 24607800*pow2(lmMst2))*
        pow6(Mst1) - 27783000*(1 + 2*lmMst2)*pow6(Mst2))) + 144*pow4(Mt)*(
        343000*Dmglst2*pow2(Mst2)*((2954 - 864*lmMst1 + 1713*lmMst2 - 849*lmMt)
        *pow2(Mst1)*pow2(Mst2) + 27*(401 + 274*lmMst2 - 2*lmMst1*(71 + 2*lmMst2
        - 2*lmMt) - 4*(33 + lmMst2)*lmMt + 4*pow2(lmMst2))*pow4(Mst1) + 324*
        pow4(Mst2)) - Mgl*(30870*(23941 + 30*lmMst2*(667 - 10*lmMt) - 30*
        lmMst1*(347 + 10*lmMst2 - 10*lmMt) - 9600*lmMt + 300*pow2(lmMst2))*
        pow2(Mst2)*pow4(Mst1) + 42875*(5567 - 2232*lmMst1 - 2838*lmMt + 6*
        lmMst2*(917 + 36*lmMt) - 108*pow2(lmMst2) - 108*pow2(lmMt))*pow2(Mst1)*
        pow4(Mst2) + 24*(26331136 - 105*lmMst1*(236317 + 35070*lmMst2 - 36750*
        lmMt) + 210*lmMst2*(175121 - 18375*lmMt) - 11962125*lmMt - 88200*pow2(
        lmMst1) + 3770550*pow2(lmMst2))*pow6(Mst1) + 9261000*(1 + 2*lmMst2)*
        pow6(Mst2))) - 4*pow2(Mt)*pow2(s2t)*(9261000*Dmglst2*pow2(Mst2)*(((678
        - 36*lmMst1 + 36*lmMst2)*pow2(Mst2) + 9*(-23 + 42*lmMst1 - 42*lmMst2)*
        pow2(MuSUSY))*pow4(Mst1) + 108*(2*pow2(Mst2) + pow2(MuSUSY))*pow4(Mst2)
        + 12*pow2(Mst1)*(6*(4 + lmMst1 - lmMst2)*pow2(Mst2)*pow2(MuSUSY) + 23*
        pow4(Mst2)) + (656 - 192*lmMst1 + 192*lmMst2)*pow6(Mst1)) - Mgl*(20580*
        pow2(Mst2)*((2482 - 400*OepS2 + 90*lmMst2*(443 - 90*S2) + 90450*S2 -
        90*lmMst1*(263 - 90*(lmMst2 + S2)) - 4050*(pow2(lmMst1) + pow2(lmMst2))
        )*pow2(Mst2) - (45*lmMst1*(-1547 + 180*lmMst2 - 4500*S2) + 45*lmMst2*(
        1547 + 4500*S2) + 2*(-106283 + 5000*OepS2 + 639225*S2) + 4050*pow2(
        lmMst1) - 12150*pow2(lmMst2))*pow2(MuSUSY))*pow4(Mst1) + 385875*pow2(
        Mst1)*(2*(3439 - 64*OepS2 + 8424*S2 + 144*lmMst1*(1 + 9*S2) - 144*
        lmMst2*(7 + 9*S2))*pow2(Mst2) + (8771 - 128*OepS2 - 57024*S2 - 72*
        lmMst2*(23 + 36*S2) + 72*lmMst1*(29 - 12*lmMst2 + 36*S2) + 864*pow2(
        lmMst2))*pow2(MuSUSY))*pow4(Mst2) + 2*((194011877 + 9604000*OepS2 -
        319504500*S2 + 1260*lmMst2*(60437 + 154350*S2) - 1260*lmMst1*(60437 +
        15120*lmMst2 + 154350*S2) + 9525600*(pow2(lmMst1) + pow2(lmMst2)))*
        pow2(Mst2) + 5*(593331163 - 60642400*OepS2 + 1260*lmMst2*(8263 -
        974610*S2) + 1143733500*S2 + 1260*lmMst1*(-8263 + 71820*lmMst2 +
        974610*S2) - 61916400*pow2(lmMst1) - 28576800*pow2(lmMst2))*pow2(
        MuSUSY))*pow6(Mst1) + 83349000*(1 + 2*lmMst2)*(2*pow2(Mst2) + pow2(
        MuSUSY))*pow6(Mst2)))))) - 3087000*Mst2*xDmglst2*pow2(Dmglst2)*(Tbeta*(
        -1728*s2t*pow2(Mst1)*((8 - 15*lmMst1 + 15*lmMst2)*pow2(Mst1) + (-2 - 5*
        lmMst1 + 5*lmMst2)*pow2(Mst2))*pow2(MuSUSY)*pow3(Mt) + 18*Mst2*pow2(Mt)
        *pow2(s2t)*(pow2(MuSUSY)*(-48*(4 + lmMst1 - lmMst2)*pow2(Mst1)*pow2(
        Mst2) + 6*(23 - 42*lmMst1 + 42*lmMst2)*pow4(Mst1) - 72*pow4(Mst2)) - (-
        1081 + 165*lmMst1 - 165*lmMst2)*(-2 + pow2(Sbeta))*pow4(Sbeta)*pow6(
        Mst1)) + pow2(Sbeta)*(-216*Mt*pow2(Mst1)*pow2(Mst2)*pow3(s2t)*(4*(12 -
        5*lmMst1 + 5*lmMst2)*pow2(Mst1)*pow2(Mst2) + (40 - 3*lmMst1 + 3*lmMst2)
        *pow4(Mst1) - 4*(2 + 5*lmMst1 - 5*lmMst2)*pow4(Mst2)) + 48*s2t*pow2(
        Mst1)*pow3(Mt)*(-36*(2 + 5*lmMst1 - 5*lmMst2)*pow2(Mst2)*pow2(MuSUSY) +
        36*pow2(Mst1)*(8*(7 - 3*lmMst1 + 6*lmMst2 - 3*lmMt)*pow2(Mst2) + (8 -
        15*lmMst1 + 15*lmMst2)*pow2(MuSUSY)) + 9*(465 + 374*lmMst2 - 4*lmMst1*(
        52 + 3*lmMst2 - 3*lmMt) - 2*(83 + 6*lmMst2)*lmMt + 12*pow2(lmMst2))*
        pow4(Mst1) + (578 - 360*lmMst1 + 708*lmMst2 - 348*lmMt)*pow4(Mst2)) -
        16*Mst2*((2954 - 864*lmMst1 + 1713*lmMst2 - 849*lmMt)*pow2(Mst1)*pow2(
        Mst2) + 27*(401 + 274*lmMst2 - 2*lmMst1*(71 + 2*lmMst2 - 2*lmMt) - 4*(
        33 + lmMst2)*lmMt + 4*pow2(lmMst2))*pow4(Mst1) + 324*pow4(Mst2))*pow4(
        Mt) + 12*Mst2*pow2(Mt)*pow2(s2t)*(3*(2*(113 - 6*lmMst1 + 6*lmMst2)*
        pow2(Mst2) + 3*(-23 + 42*lmMst1 - 42*lmMst2)*pow2(MuSUSY))*pow4(Mst1) +
        108*(2*pow2(Mst2) + pow2(MuSUSY))*pow4(Mst2) + 12*pow2(Mst1)*(6*(4 +
        lmMst1 - lmMst2)*pow2(Mst2)*pow2(MuSUSY) + 23*pow4(Mst2)) + 16*(41 -
        12*lmMst1 + 12*lmMst2)*pow6(Mst1)) - 3*pow3(Mst2)*pow4(s2t)*(9*(-75 +
        26*lmMst1 - 26*lmMst2)*pow2(Mst2)*pow4(Mst1) + 72*(1 + lmMst1 - lmMst2)
        *pow2(Mst1)*pow4(Mst2) + (76 - 249*lmMst1 + 249*lmMst2)*pow6(Mst1) +
        108*pow6(Mst2)))) + 12*Mt*MuSUSY*pow2(Sbeta)*(-4*pow2(Mst1)*pow2(Mst2)*
        (246*Mst2*s2t*pow2(Mt) + 54*(2 + 5*lmMst1 - 5*lmMst2)*Mt*pow2(Mst2)*
        pow2(s2t) + 4*(101 - 90*lmMst1 + 177*lmMst2 - 87*lmMt)*pow3(Mt) - 9*(5
        + 2*lmMst1 - 2*lmMst2)*pow3(Mst2)*pow3(s2t)) - (36*(65 - 2*lmMst1 + 2*
        lmMst2)*Mst2*s2t*pow2(Mt) + 2160*(-1 + lmMst1 - lmMst2)*Mt*pow2(Mst2)*
        pow2(s2t) + 16*(497 - 306*lmMst1 + 609*lmMst2 - 303*lmMt)*pow3(Mt) + 9*
        (55 - 34*lmMst1 + 34*lmMst2)*pow3(Mst2)*pow3(s2t))*pow4(Mst1) + (54*(80
        - 43*lmMst1 + 43*lmMst2)*Mt + (-419 + 57*lmMst1 - 57*lmMst2)*Mst2*s2t)*
        pow2(s2t)*pow6(Mst1) + 108*(-4*s2t*pow2(Mt)*pow5(Mst2) + pow3(s2t)*
        pow7(Mst2))))) + 308700*Mgl*xDR2DRMOD*pow2(Msq)*(4*MuSUSY*pow2(Sbeta)*(
        -1215*Mgl*shiftst1*(Dmsqst2 + pow2(Msq))*pow2(Mst1)*(4*s2t*(pow2(Mst1)
        - pow2(Mst2))*pow3(Mt) - Mt*pow3(s2t)*(2*(lmMst1 - lmMst2)*pow2(Mst1)*
        pow2(Mst2) + pow4(Mst1) - pow4(Mst2)))*pow5(Mst2) + 128*pow2(Msq)*pow2(
        Mst1)*pow4(Mt)*(9*(1 + lmMst2)*Mgl*((3 - 21*lmMst2 + 2*lmMst1*(8 + 3*
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
        pow6(Mst2))) - 864*pow2(Msq)*pow2(Mst1)*pow2(Mst2)*pow2(Mt)*pow2(s2t)*(
        (1 + lmMst2)*Mgl*(2*(2 - 5*lmMst2 + lmMst1*(5 + 2*lmMst2) - 2*pow2(
        lmMst2))*pow2(Mst2)*pow4(Mst1) + 2*(-2*lmMst2*(2 + lmMst2) + lmMst1*(3
        + 2*lmMst2))*pow2(Mst1)*pow4(Mst2) + (5 - 12*lmMst2 + 4*lmMst1*(3 +
        lmMst2) - 4*pow2(lmMst2))*pow6(Mst1) - 2*(1 + lmMst2)*pow6(Mst2)) +
        Dmglst2*(2*pow2(Mst2)*(8 + 19*lmMst2 - 5*pow2(lmMst2) + lmMst1*(-7 + 5*
        lmMst2 + 6*pow2(lmMst2)) - 6*pow3(lmMst2))*pow4(Mst1) + 2*pow2(Mst1)*(7
        + 9*lmMst2 - 8*pow2(lmMst2) + lmMst1*(-3 + 5*lmMst2 + 6*pow2(lmMst2)) -
        6*pow3(lmMst2))*pow4(Mst2) + (17 + 51*lmMst2 - 4*pow2(lmMst2) + 4*
        lmMst1*(-6 + lmMst2 + 3*pow2(lmMst2)) - 12*pow3(lmMst2))*pow6(Mst1) -
        2*(-1 + 2*lmMst2 + 3*pow2(lmMst2))*pow6(Mst2))) - 9*Mt*pow3(Mst2)*pow3(
        s2t)*(-2*Dmglst2*(pow2(Mst1)*pow2(Mst2)*(4*pow2(Msq)*(20*pow2(Mst1)*
        pow2(Mst2) + 9*pow4(Mst1) - 5*pow4(Mst2)) + 15*Dmsqst2*(pow4(Mst1) -
        pow4(Mst2))) - 30*Dmsqst2*lmMst2*pow4(Mst1)*pow4(Mst2) + pow4(Mst1)*(
        30*Dmsqst2*lmMst1*pow4(Mst2) + pow2(Msq)*(32*lmMst2*pow2(lmMst1)*(8*
        pow2(Mst1)*pow2(Mst2) + 8*pow4(Mst1) + pow4(Mst2)) + 2*pow3(lmMst2)*(
        128*pow2(Mst1)*(pow2(Mst1) + pow2(Mst2)) + 107*pow4(Mst2)))) + pow2(
        Msq)*(-2*lmMst1*((-4 + 246*lmMst2 + 123*pow2(lmMst2))*pow4(Mst1)*pow4(
        Mst2) + 16*(3 + 3*lmMst2 + 16*pow2(lmMst2))*pow2(Mst2)*pow6(Mst1) - 16*
        lmMst2*pow2(Mst1)*pow6(Mst2) + 16*(3 + 2*lmMst2 + 16*pow2(lmMst2))*
        pow8(Mst1)) + pow2(lmMst2)*(396*pow4(Mst1)*pow4(Mst2) + 37*pow2(Mst2)*
        pow6(Mst1) + 155*pow2(Mst1)*pow6(Mst2) + 64*pow8(Mst1) - 32*pow8(Mst2))
        ) + 2*lmMst2*pow2(Msq)*(4*pow4(Mst1)*pow4(Mst2) + 21*pow2(Mst2)*pow6(
        Mst1) + 115*pow2(Mst1)*pow6(Mst2) + 56*pow8(Mst1) - 16*pow8(Mst2))) +
        Mgl*(75*Dmsqst2*pow2(Mst1)*pow2(Mst2)*(2*(lmMst1 - lmMst2)*pow2(Mst1)*
        pow2(Mst2) + pow4(Mst1) - pow4(Mst2)) - pow2(Msq)*(2*(-8 + 237*lmMst2 +
        16*(1 + lmMst2)*pow2(lmMst1) + 308*pow2(lmMst2) - lmMst1*(269 + 348*
        lmMst2 + 123*pow2(lmMst2)) + 107*pow3(lmMst2))*pow4(Mst1)*pow4(Mst2) +
        pow2(Mst2)*(-125 - 156*lmMst2 - 512*lmMst1*lmMst2*(1 + lmMst2) + 256*(1
        + lmMst2)*pow2(lmMst1) + 181*pow2(lmMst2) + 256*pow3(lmMst2))*pow6(
        Mst1) + (221 + 284*lmMst2 + 32*lmMst1*(1 + lmMst2) + 107*pow2(lmMst2))*
        pow2(Mst1)*pow6(Mst2) + 16*(1 + lmMst2)*(1 + lmMst1*(2 - 32*lmMst2) -
        2*lmMst2 + 16*pow2(lmMst1) + 16*pow2(lmMst2))*pow8(Mst1) - 16*pow2(1 +
        lmMst2)*pow8(Mst2)))) - 12*Mst2*s2t*pow3(Mt)*(45*Dmsqst2*(2*Dmglst2 -
        5*Mgl)*pow2(Mst1)*(pow2(Mst1) - pow2(Mst2))*pow4(Mst2) + pow2(Msq)*(-2*
        Dmglst2*((52 + 370*lmMst2 - 96*lmMst2*pow2(lmMst1) + 81*pow2(lmMst2) +
        96*lmMst1*(-1 + lmMst2 + 2*pow2(lmMst2)) - 96*pow3(lmMst2))*pow4(Mst1)*
        pow4(Mst2) - 96*(lmMst2*pow2(lmMst1) + lmMst1*(4 + 2*lmMst2 - 2*pow2(
        lmMst2)) + (-4 + lmMst2)*pow2(1 + lmMst2))*pow2(Mst2)*pow6(Mst1) - 3*(-
        52 + 2*(51 + 16*lmMst1)*lmMst2 + 59*pow2(lmMst2))*pow2(Mst1)*pow6(Mst2)
        - 96*(-6 - 14*lmMst2 + lmMst2*pow2(lmMst1) + lmMst1*(9 + 6*lmMst2 - 2*
        pow2(lmMst2)) - 6*pow2(lmMst2) + pow3(lmMst2))*pow8(Mst1) + 96*lmMst2*(
        1 + lmMst2)*pow8(Mst2)) + 3*Mgl*(-((109 + 76*lmMst2 - 21*pow2(lmMst2) +
        64*lmMst1*pow2(1 + lmMst2) - 32*((1 + lmMst2)*pow2(lmMst1) + pow3(
        lmMst2)))*pow4(Mst1)*pow4(Mst2)) + (173 + 188*lmMst2 + 32*lmMst1*(1 +
        lmMst2) + 59*pow2(lmMst2))*pow2(Mst1)*pow6(Mst2) + 32*(1 + lmMst2)*(
        pow2(1 - lmMst1 + lmMst2)*pow2(Mst2)*pow6(Mst1) + (1 + 3*lmMst2 -
        lmMst1*(3 + 2*lmMst2) + pow2(lmMst1) + pow2(lmMst2))*pow8(Mst1)) - 16*
        pow2(1 + lmMst2)*pow8(Mst2))))) + Tbeta*(1215*Mgl*Mst2*shiftst1*(
        Dmsqst2 + pow2(Msq))*pow2(Mst1)*(4*pow2(Mt)*pow2(MuSUSY)*pow2(s2t)*(
        pow2(Mst1)*pow4(Mst2) - 2*(lmMst1 - lmMst2)*pow2(Mst1)*(pow2(Mst1)*
        pow2(Mst2) + pow4(Mst1) + pow4(Mst2)) + pow6(Mst2)) + pow2(Sbeta)*(16*(
        pow2(Mst1) + pow2(Mst2))*pow4(Mst2)*pow4(Mt) + (pow2(Mst1) - pow2(Mst2)
        )*(2*(lmMst1 - lmMst2)*pow2(Mst1)*pow2(Mst2) + pow4(Mst1) - pow4(Mst2))
        *pow4(Mst2)*pow4(s2t) - 4*pow2(Mt)*pow2(s2t)*(2*pow4(Mst1)*pow4(Mst2) -
        4*pow2(Mst1)*pow6(Mst2) + pow2(MuSUSY)*(pow2(Mst1)*pow4(Mst2) - 2*(
        lmMst1 - lmMst2)*pow2(Mst1)*(pow2(Mst1)*pow2(Mst2) + pow4(Mst1) + pow4(
        Mst2)) + pow6(Mst2)) + 2*pow8(Mst2)))) - 36*s2t*pow2(Mt)*pow2(MuSUSY)*(
        Dmsqst2*s2t*pow2(Mst1)*(-30*Dmglst2*(pow2(Mst1) + pow2(Mst2))*(-2*
        lmMst1*pow2(Mst1) + 2*lmMst2*pow2(Mst1) + pow2(Mst2))*pow3(Mst2) + 75*
        Mgl*Mst2*((pow2(Mst1) + pow2(Mst2))*pow4(Mst2) - 2*(lmMst1 - lmMst2)*
        pow2(Mst1)*(pow2(Mst1)*pow2(Mst2) + pow4(Mst1) + pow4(Mst2)))) + pow2(
        Msq)*(-2*Dmglst2*((64*Mt*(8 + 7*lmMst2 - 11*pow2(lmMst2) + lmMst1*(-3 +
        5*lmMst2 + 6*pow2(lmMst2)) - 6*pow3(lmMst2)) - Mst2*s2t*(60 + 206*
        lmMst2 + 32*lmMst2*pow2(lmMst1) + lmMst1*(8 - 460*lmMst2 - 246*pow2(
        lmMst2)) + 519*pow2(lmMst2) + 214*pow3(lmMst2)))*pow4(Mst1)*pow4(Mst2)
        + (-64*Mt*(-1 + 2*lmMst2 + 3*pow2(lmMst2)) - Mst2*s2t*(-20 + 2*(99 +
        16*lmMst1)*lmMst2 + 123*pow2(lmMst2)))*pow2(Mst1)*pow6(Mst2) + 2*(pow2(
        Mst2)*(64*Mt*(8 + 13*lmMst2 - 8*pow2(lmMst2) + lmMst1*(-5 + 5*lmMst2 +
        6*pow2(lmMst2)) - 6*pow3(lmMst2)) - Mst2*s2t*(48 + 4*lmMst2*(31 + 36*
        pow2(lmMst1)) + 278*pow2(lmMst2) - lmMst1*(44 + 278*lmMst2 + 379*pow2(
        lmMst2)) + 235*pow3(lmMst2)))*pow6(Mst1) + (16*Mt*(49 + 103*lmMst2 -
        36*(1 + lmMst2)*pow2(lmMst2) + 4*lmMst1*(-11 + 6*lmMst2 + 9*pow2(
        lmMst2))) - Mst2*s2t*(48 + 4*lmMst2*(45 + 68*pow2(lmMst1)) + 310*pow2(
        lmMst2) - lmMst1*(92 + 310*lmMst2 + 635*pow2(lmMst2)) + 363*pow3(
        lmMst2)))*pow8(Mst1)) + 32*lmMst2*(1 + lmMst2)*s2t*pow9(Mst2)) + Mgl*((
        128*(1 + lmMst2)*Mt*(1 + 5*lmMst2 - lmMst1*(3 + 2*lmMst2) + 2*pow2(
        lmMst2)) + Mst2*s2t*(189 + 726*lmMst2 + 32*(1 + lmMst2)*pow2(lmMst1) +
        707*pow2(lmMst2) - 2*lmMst1*(253 + 332*lmMst2 + 123*pow2(lmMst2)) +
        214*pow3(lmMst2)))*pow4(Mst1)*pow4(Mst2) + (Mst2*s2t*(205 + 252*lmMst2
        + 32*lmMst1*(1 + lmMst2) + 91*pow2(lmMst2)) + 128*Mt*pow2(1 + lmMst2))*
        pow2(Mst1)*pow6(Mst2) - 2*(pow2(Mst2)*(64*(1 + lmMst2)*Mt*(1 - 10*
        lmMst2 + 4*lmMst1*(2 + lmMst2) - 4*pow2(lmMst2)) - Mst2*s2t*(32 + 285*
        lmMst2 + 144*(1 + lmMst2)*pow2(lmMst1) + 444*pow2(lmMst2) - lmMst1*(253
        + 588*lmMst2 + 379*pow2(lmMst2)) + 235*pow3(lmMst2)))*pow6(Mst1) + (32*
        (1 + lmMst2)*Mt*(7 - 32*lmMst2 + 4*lmMst1*(7 + 3*lmMst2) - 12*pow2(
        lmMst2)) - Mst2*s2t*(40 + 277*lmMst2 + 272*(1 + lmMst2)*pow2(lmMst1) +
        556*pow2(lmMst2) - lmMst1*(237 + 828*lmMst2 + 635*pow2(lmMst2)) + 363*
        pow3(lmMst2)))*pow8(Mst1)) - 16*s2t*pow2(1 + lmMst2)*pow9(Mst2)))) +
        pow2(Sbeta)*(256*s2t*pow2(Msq)*pow2(Mst1)*pow3(Mt)*(9*(1 + lmMst2)*Mgl*
        (-2*pow2(Mst2)*((3 + 2*lmMst1*(7 + 2*lmMst2 - 2*lmMt) + 4*lmMst2*(-4 +
        lmMt) + 2*lmMt - 4*pow2(lmMst2))*pow2(Mst2) + (1 - 10*lmMst2 + 4*
        lmMst1*(2 + lmMst2) - 4*pow2(lmMst2))*pow2(MuSUSY))*pow4(Mst1) + pow2(
        Mst1)*((1 - 2*lmMst1*(5 + 2*lmMst2 - 2*lmMt) - 4*lmMst2*(-3 + lmMt) -
        6*lmMt + 4*pow2(lmMst2))*pow2(Mst2) + 2*(1 + 5*lmMst2 - lmMst1*(3 + 2*
        lmMst2) + 2*pow2(lmMst2))*pow2(MuSUSY))*pow4(Mst2) - (2*(8 - 27*lmMst2
        + lmMst1*(25 + 6*lmMst2 - 6*lmMt) + 2*lmMt + 6*lmMst2*lmMt - 6*pow2(
        lmMst2))*pow2(Mst2) + (7 - 32*lmMst2 + 4*lmMst1*(7 + 3*lmMst2) - 12*
        pow2(lmMst2))*pow2(MuSUSY))*pow6(Mst1) + 2*(1 + lmMst2)*(2*pow2(Mst2) +
        pow2(MuSUSY))*pow6(Mst2)) - Dmglst2*(pow2(Mst1)*(pow2(Mst2)*(253 + 5*
        lmMst2*(107 - 6*lmMt) - 102*lmMt - 18*lmMst1*(9 + lmMst2 - 2*lmMt + 2*
        lmMst2*lmMt - 2*pow2(lmMst2)) + 12*(10 + 3*lmMt)*pow2(lmMst2) - 36*
        pow3(lmMst2)) + 18*pow2(MuSUSY)*(8 + 7*lmMst2 - 11*pow2(lmMst2) +
        lmMst1*(-3 + 5*lmMst2 + 6*pow2(lmMst2)) - 6*pow3(lmMst2)))*pow4(Mst2) +
        9*((49 + 103*lmMst2 - 36*(1 + lmMst2)*pow2(lmMst2) + 4*lmMst1*(-11 + 6*
        lmMst2 + 9*pow2(lmMst2)))*pow2(MuSUSY) + pow2(Mst2)*(28 + lmMst2*(378 -
        96*lmMt) - 44*lmMt + lmMst1*(-274 + 60*lmMt + 18*lmMst2*(-13 + 2*lmMt)
        - 36*pow2(lmMst2)) + (270 - 36*lmMt)*pow2(lmMst2) + 36*pow3(lmMst2)))*
        pow6(Mst1) + 18*(pow2(Mst2)*(2*pow2(MuSUSY)*(8 + 13*lmMst2 - 8*pow2(
        lmMst2) + lmMst1*(-5 + 5*lmMst2 + 6*pow2(lmMst2)) - 6*pow3(lmMst2)) +
        pow2(Mst2)*(23 + lmMst2*(93 - 22*lmMt) - 14*lmMt - 4*(-11 + lmMt)*pow2(
        lmMst2) - 2*lmMst1*(25 + lmMst2*(17 - 2*lmMt) - 6*lmMt + 2*pow2(lmMst2)
        ) + 4*pow3(lmMst2)))*pow4(Mst1) + (1 - 2*lmMst2 - 3*pow2(lmMst2))*(2*
        pow2(Mst2) + pow2(MuSUSY))*pow6(Mst2)))) + 9*pow2(Mst1)*pow3(s2t)*pow4(
        Mst2)*(128*Mt*pow2(Msq)*((1 + lmMst2)*Mgl*(2*(2 + 2*lmMst1 - lmMst2)*
        pow2(Mst2)*pow4(Mst1) + 2*(1 - 3*lmMst2 + lmMst1*(3 + 2*lmMst2) - 2*
        pow2(lmMst2))*pow2(Mst1)*pow4(Mst2) + (1 + 2*lmMst1 - 2*lmMst2)*pow6(
        Mst1) - 2*(1 + lmMst2)*pow6(Mst2)) + Dmglst2*(2*(1 - 4*lmMst1 + 10*
        lmMst2 + 3*pow2(lmMst2))*pow2(Mst2)*pow4(Mst1) + 2*pow2(Mst1)*(6 + 11*
        lmMst2 - 5*pow2(lmMst2) + lmMst1*(-3 + 5*lmMst2 + 6*pow2(lmMst2)) - 6*
        pow3(lmMst2))*pow4(Mst2) + (1 + 13*lmMst2 - 2*lmMst1*(5 + 3*lmMst2) +
        6*pow2(lmMst2))*pow6(Mst1) - 2*(-1 + 2*lmMst2 + 3*pow2(lmMst2))*pow6(
        Mst2))) + Mst2*s2t*(Dmsqst2*(30*Dmglst2 - 75*Mgl)*(pow2(Mst1) - pow2(
        Mst2))*(2*(lmMst1 - lmMst2)*pow2(Mst1)*pow2(Mst2) + pow4(Mst1) - pow4(
        Mst2)) + pow2(Msq)*(-(Mgl*(pow2(Mst2)*(-109 - 630*lmMst2 + 224*(1 +
        lmMst2)*pow2(lmMst1) + lmMst1*(538 + 184*lmMst2 - 266*pow2(lmMst2)) -
        435*pow2(lmMst2) + 42*pow3(lmMst2))*pow4(Mst1) + pow2(Mst1)*(-237 +
        190*lmMst2 + 32*(1 + lmMst2)*pow2(lmMst1) + 509*pow2(lmMst2) - 2*
        lmMst1*(285 + 364*lmMst2 + 123*pow2(lmMst2)) + 214*pow3(lmMst2))*pow4(
        Mst2) + (141 + 140*lmMst2 + 32*lmMst1*(1 + lmMst2) + 43*pow2(lmMst2))*
        pow6(Mst1) + (237 + 316*lmMst2 + 32*lmMst1*(1 + lmMst2) + 123*pow2(
        lmMst2))*pow6(Mst2))) - 2*Dmglst2*(pow2(Mst2)*(-44 + 34*lmMst2 + 224*
        lmMst2*pow2(lmMst1) - 359*pow2(lmMst2) - 2*lmMst1*(52 - 198*lmMst2 +
        133*pow2(lmMst2)) + 42*pow3(lmMst2))*pow4(Mst1) + pow2(Mst1)*(100 -
        222*lmMst2 + 32*lmMst2*pow2(lmMst1) + lmMst1*(8 - 524*lmMst2 - 246*
        pow2(lmMst2)) + 241*pow2(lmMst2) + 214*pow3(lmMst2))*pow4(Mst2) + (-36
        + (70 + 32*lmMst1)*lmMst2 + 27*pow2(lmMst2))*pow6(Mst1) + (-20 + (262 +
        32*lmMst1)*lmMst2 + 187*pow2(lmMst2))*pow6(Mst2))))) + 16*Mst2*pow4(Mt)
        *(135*Dmsqst2*(2*Dmglst2 - 5*Mgl)*pow2(Mst1)*(pow2(Mst1) + pow2(Mst2))*
        pow4(Mst2) + pow2(Msq)*(2*Dmglst2*((4220 + lmMst2*(6114 - 1296*lmMt) -
        576*lmMst1*(4 + 2*lmMst2 - lmMt) - 1392*lmMt + 1341*pow2(lmMst2))*pow4(
        Mst1)*pow4(Mst2) - 9*(-84 + (70 + 32*lmMst1)*lmMst2 + 59*pow2(lmMst2))*
        pow2(Mst1)*pow6(Mst2) + 288*(pow2(Mst2)*(31 + lmMst2*(95 - 23*lmMt) -
        16*lmMt + (51 - 6*lmMt)*pow2(lmMst2) - 2*lmMst1*(25 + lmMst2*(20 - 3*
        lmMt) - 6*lmMt + 3*pow2(lmMst2)) + 6*pow3(lmMst2))*pow6(Mst1) + (42 +
        lmMst2*(242 - 62*lmMt) - 33*lmMt + 6*(29 - 4*lmMt)*pow2(lmMst2) - 2*
        lmMst1*(81 + lmMst2*(74 - 12*lmMt) - 18*lmMt + 12*pow2(lmMst2)) + 24*
        pow3(lmMst2))*pow8(Mst1)) + 288*lmMst2*(1 + lmMst2)*pow8(Mst2)) + Mgl*(
        -((2117 + lmMst2*(4796 - 1248*lmMt) - 576*lmMst1*(1 + lmMst2)*(3 +
        lmMst2 - lmMt) - 672*lmMt + (3651 - 576*lmMt)*pow2(lmMst2) + 576*pow3(
        lmMst2))*pow4(Mst1)*pow4(Mst2)) + 288*(1 + lmMst2)*((5 - 57*lmMst2 + 2*
        lmMst1*(25 + 6*lmMst2 - 6*lmMt) + 7*lmMt + 12*lmMst2*lmMt - 12*pow2(
        lmMst2))*pow2(Mst1) + (-2 - 27*lmMst2 + lmMst1*(22 + 6*lmMst2 - 6*lmMt)
        + 5*lmMt + 6*lmMst2*lmMt - 6*pow2(lmMst2))*pow2(Mst2))*pow6(Mst1) - 9*(
        173 + 188*lmMst2 + 32*lmMst1*(1 + lmMst2) + 59*pow2(lmMst2))*pow2(Mst1)
        *pow6(Mst2) + 144*pow2(1 + lmMst2)*pow8(Mst2)))) - 12*Mst2*pow2(Mt)*
        pow2(s2t)*(2*Dmglst2*(45*Dmsqst2*pow2(Mst1)*pow2(Mst2)*(2*(pow2(Mst2) +
        (-lmMst1 + lmMst2)*pow2(MuSUSY))*pow4(Mst1) + pow2(Mst1)*((1 - 2*lmMst1
        + 2*lmMst2)*pow2(Mst2)*pow2(MuSUSY) - 4*pow4(Mst2)) + (2*pow2(Mst2) +
        pow2(MuSUSY))*pow4(Mst2)) + pow2(Msq)*(2*(-3*pow2(Mst2)*pow2(MuSUSY)*(
        48 + 4*lmMst2*(31 + 36*pow2(lmMst1)) + 278*pow2(lmMst2) - lmMst1*(44 +
        278*lmMst2 + 379*pow2(lmMst2)) + 235*pow3(lmMst2)) + (332 + 302*lmMst2
        - 288*lmMst1*(1 + lmMst2) + 111*pow2(lmMst2))*pow4(Mst2))*pow6(Mst1) -
        pow4(Mst1)*(3*pow2(MuSUSY)*(60 + 206*lmMst2 + 32*lmMst2*pow2(lmMst1) +
        lmMst1*(8 - 460*lmMst2 - 246*pow2(lmMst2)) + 519*pow2(lmMst2) + 214*
        pow3(lmMst2))*pow4(Mst2) + 4*(52 + lmMst2*(-338 + 48*pow2(lmMst1)) -
        129*pow2(lmMst2) + 48*(lmMst1 - 2*lmMst1*lmMst2*(1 + lmMst2) + pow3(
        lmMst2)))*pow6(Mst2)) + 6*(32*(2 + 7*lmMst2 - lmMst1*(5 + 4*lmMst2) +
        4*pow2(lmMst2))*pow2(Mst2) - pow2(MuSUSY)*(48 + 4*lmMst2*(45 + 68*pow2(
        lmMst1)) + 310*pow2(lmMst2) - lmMst1*(92 + 310*lmMst2 + 635*pow2(
        lmMst2)) + 363*pow3(lmMst2)))*pow8(Mst1) + 96*lmMst2*(1 + lmMst2)*(2*
        pow2(Mst2) + pow2(MuSUSY))*pow8(Mst2) + 3*pow2(Mst1)*(-((-20 + 2*(99 +
        16*lmMst1)*lmMst2 + 123*pow2(lmMst2))*pow2(MuSUSY)*pow6(Mst2)) - 2*(-52
        + 2*(67 + 16*lmMst1)*lmMst2 + 91*pow2(lmMst2))*pow8(Mst2)))) - 3*Mgl*(
        6*(25*Dmsqst2 + 47*pow2(Msq))*pow4(Mst2)*pow6(Mst1) - 12*(25*Dmsqst2 +
        47*pow2(Msq))*pow4(Mst1)*pow6(Mst2) + pow2(MuSUSY)*(3*(25*Dmsqst2 + 63*
        pow2(Msq))*pow4(Mst1)*pow4(Mst2) + 5*(15*Dmsqst2 + 41*pow2(Msq))*pow2(
        Mst1)*pow6(Mst2)) - 2*lmMst1*(75*Dmsqst2*pow2(MuSUSY)*pow4(Mst1)*(pow2(
        Mst1)*pow2(Mst2) + pow4(Mst1) + pow4(Mst2)) + pow2(Msq)*((32*(3 + 5*
        lmMst2 + 2*pow2(lmMst2))*pow2(Mst2) + (253 + 332*lmMst2 + 123*pow2(
        lmMst2))*pow2(MuSUSY))*pow4(Mst1)*pow4(Mst2) + (253 + 588*lmMst2 + 379*
        pow2(lmMst2))*pow2(Mst2)*pow2(MuSUSY)*pow6(Mst1) - 16*(1 + lmMst2)*
        pow2(Mst1)*(2*pow2(Mst2) + pow2(MuSUSY))*pow6(Mst2) + (32*(1 + lmMst2)*
        pow2(Mst2) + (237 + 828*lmMst2 + 635*pow2(lmMst2))*pow2(MuSUSY))*pow8(
        Mst1))) + 6*(25*Dmsqst2 + 63*pow2(Msq))*pow2(Mst1)*pow8(Mst2) + 2*
        lmMst2*(75*Dmsqst2*pow2(MuSUSY)*pow4(Mst1)*(pow2(Mst1)*pow2(Mst2) +
        pow4(Mst1) + pow4(Mst2)) + pow2(Msq)*((285*pow2(Mst2)*pow2(MuSUSY) +
        172*pow4(Mst2))*pow6(Mst1) - 33*pow4(Mst1)*(-11*pow2(MuSUSY)*pow4(Mst2)
        + 8*pow6(Mst2)) + (32*pow2(Mst2) + 277*pow2(MuSUSY))*pow8(Mst1) - 16*(
        2*pow2(Mst2) + pow2(MuSUSY))*pow8(Mst2) + 2*pow2(Mst1)*(63*pow2(MuSUSY)
        *pow6(Mst2) + 110*pow8(Mst2)))) + pow2(Msq)*(pow4(Mst1)*(32*(1 +
        lmMst2)*pow2(lmMst1)*(pow2(MuSUSY)*(9*pow2(Mst1)*pow2(Mst2) + 17*pow4(
        Mst1) + pow4(Mst2)) + 2*pow6(Mst2)) + 2*pow3(lmMst2)*(pow2(MuSUSY)*(
        235*pow2(Mst1)*pow2(Mst2) + 363*pow4(Mst1) + 107*pow4(Mst2)) + 32*pow6(
        Mst2))) - 16*pow2(MuSUSY)*pow8(Mst2) - pow2(lmMst2)*(-6*(148*pow2(Mst2)
        *pow2(MuSUSY) + 25*pow4(Mst2))*pow6(Mst1) + pow4(Mst1)*(-707*pow2(
        MuSUSY)*pow4(Mst2) + 76*pow6(Mst2)) - 8*(8*pow2(Mst2) + 139*pow2(
        MuSUSY))*pow8(Mst1) + 16*(2*pow2(Mst2) + pow2(MuSUSY))*pow8(Mst2) -
        pow2(Mst1)*(91*pow2(MuSUSY)*pow6(Mst2) + 150*pow8(Mst2)))) + pow2(Msq)*
        (16*pow2(MuSUSY)*(4*pow2(Mst2)*pow6(Mst1) + 5*pow8(Mst1)) - 32*power10(
        Mst2)))))))))/(1.000188e8*pow4(Msq)*pow5(Mst2)))/(Tbeta*pow2(Sbeta)*
        pow4(Mst1)))/pow2(Mgl);
}

/**
 * @return The off-diagonal (1, 2) = (2, 1) matrix element of the Higgs mass matrix as a double for the hierarchy 'H6bq2g2'
 */
double H6bq2g2::getS12() const {
   return Mt*((xMst*(27*(lmMst1 - lmMst2)*Mst2*Mt*oneLoopFlag*pow2(Mgl)*pow2(
        MuSUSY)*pow2(s2t) - Al4p*twoLoopFlag*(8*Mt*s2t*(18*(lmMst1 - lmMst2)*(-
        2 + 3*lmMst2)*Mst2*s2t*xDmglst2*xDR2DRMOD*pow2(Dmglst2) + Dmglst2*Mgl*(
        Mt*(785 + 6*lmMst1*(85 - 42*lmMst2) - 438*lmMst2 + 252*pow2(lmMst2)) -
        Mst2*s2t*(49 - 84*lmMst2 + lmMst1*(84 - 36*lmMst2*(-1 + xDR2DRMOD)) +
        36*(-1 + xDR2DRMOD)*pow2(lmMst2))) + (Mt*(193 + 474*lmMst2 - 6*lmMst1*(
        67 + 42*lmMst2) + 252*pow2(lmMst2)) - Mst2*s2t*(1 + 3*lmMst2*(-37 + 6*
        xDR2DRMOD) - 3*lmMst1*(-37 + 6*lmMst2*(-12 + xDR2DRMOD) + 6*xDR2DRMOD)
        - 81*pow2(lmMst1) + 9*(-15 + 2*xDR2DRMOD)*pow2(lmMst2)))*pow2(Mgl))*
        pow2(MuSUSY) + Mgl*MuSUSY*Tbeta*(8*Dmglst2*(36*(-1 + 3*lmMst1 - 3*
        lmMst2)*Mst2*s2t*pow2(Mt) - 3*Mt*(50 + lmMst1*(51 - 18*lmMst2) - 51*
        lmMst2 + 18*pow2(lmMst2))*pow2(Mst2)*pow2(s2t) + 4*(7 - 381*lmMst2 +
        lmMst1*(327 + 72*lmMst2 - 81*lmMt) + 54*lmMt + 81*lmMst2*lmMt - 72*
        pow2(lmMst2))*pow3(Mt) + 2*(1 + 3*lmMst1 - 3*lmMst2)*pow3(Mst2)*pow3(
        s2t)) - Mgl*(-144*Mst2*s2t*(-lmMst1 + lmMst2 - 2*lmMst1*lmMst2 + pow2(
        lmMst1) + pow2(lmMst2))*pow2(Mt) + 24*Mt*(1 + 33*lmMst2 - 3*lmMst1*(11
        + 6*lmMst2) + 18*pow2(lmMst2))*pow2(Mst2)*pow2(s2t) + 32*(83 + 3*
        lmMst1*(29 + 12*lmMst2 - 9*lmMt) + 9*lmMt + 3*lmMst2*(-32 + 9*lmMt) -
        36*pow2(lmMst2))*pow3(Mt) + (5 + lmMst1*(6 - 288*lmMst2) - 6*lmMst2 +
        144*(pow2(lmMst1) + pow2(lmMst2)))*pow3(Mst2)*pow3(s2t))) + 2*xDmglst2*
        pow2(Dmglst2)*(-3*Mst2*Mt*pow2(s2t)*(2*Mst2*MuSUSY*Tbeta*(355 + 6*
        lmMst2 - 6*lmMst1*(1 + 6*lmMst2) + 36*pow2(lmMst2)) + 2*(85 + 60*lmMst2
        + 12*lmMst1*(-5 + 3*lmMst2) - 36*pow2(lmMst2))*pow2(MuSUSY) + 5*(-43 +
        60*lmMst1 - 60*lmMst2)*pow2(Mst2)*(-1 + pow2(Sbeta))*pow2(Sbeta)) +
        MuSUSY*(2*s2t*pow2(Mt)*(MuSUSY*(5002 + 84*lmMst1*(7 - 6*lmMst2) - 444*
        lmMst2 + 504*pow2(lmMst2)) + 3*Mst2*Tbeta*(347 - 430*pow2(Sbeta) + 120*
        (lmMst1 - lmMst2)*(-4 + 5*pow2(Sbeta)))) + Tbeta*(8*(412 - 6*lmMst1*(
        176 + 42*lmMst2 - 39*lmMt) + 9*lmMst2*(147 - 26*lmMt) - 267*lmMt + 252*
        pow2(lmMst2))*pow3(Mt) + 12*(1 - 6*lmMst1 + 6*lmMst2)*pow3(Mst2)*pow3(
        s2t))))))*pow6(Mst1))/(108.*Tbeta*pow2(Mgl)*pow7(Mst2)) + (MuSUSY*(3*
        oneLoopFlag*s2t*(4*(lmMst1 - lmMst2)*pow2(Mt) - ((2 + lmMst1 - lmMst2)*
        pow2(Mst1) + (-2 + lmMst1 - lmMst2)*pow2(Mst2))*pow2(s2t) - (2*Mt*
        MuSUSY*s2t*(2 - lmMst1 + lmMst2 - (2*(lmMst1 - lmMst2)*pow2(Mst1)*(
        pow2(Mst1) + pow2(Mst2)))/pow4(Mst2)))/Tbeta) + 4*Al4p*twoLoopFlag*(12*
        Mt*pow2(s2t)*(2*Mst2*(4 + 3*lmMst2 - lmMst1*(1 + lmMst2) + pow2(lmMst2)
        + (Dmglst2*(6 + lmMst1 + lmMst2 - lmMst1*lmMst2 + pow2(lmMst2)))/Mgl) +
        ((2 - 6*lmMst1 + (6 - 4*lmMst1)*lmMst2 + 4*pow2(lmMst2) + (Dmglst2*(10
        + lmMst1*(6 - 4*lmMst2) - 6*lmMst2 + 4*pow2(lmMst2)))/Mgl)*pow2(Mst1))/
        Mst2 + ((0.5 - 7*lmMst1 + (7 - 4*lmMst1)*lmMst2 + 4*pow2(lmMst2) + (
        Dmglst2*(10.5 + lmMst1*(9 - 4*lmMst2) - 9*lmMst2 + 4*pow2(lmMst2)))/
        Mgl)*pow4(Mst1))/pow3(Mst2)) + 8*s2t*pow2(Mt)*(6*lmMst1 + 2*(-2 +
        lmMst1)*lmMst2 + pow2(lmMst1) - 3*pow2(lmMst2) + (2*(Dmglst2 + 2*
        Dmglst2*lmMst2 + ((1 - lmMst1 + lmMst2)*(2*Dmglst2 + (lmMst1 - lmMst2)*
        Mgl)*pow2(Mst1))/pow2(Mst2) - ((2*Dmglst2*lmMst2 + Mgl + lmMst2*Mgl)*
        pow2(Mst2))/pow2(Mst1) - ((Dmglst2*(-2 + 4*lmMst1 - 4*lmMst2) + Mgl*(-
        lmMst1 + lmMst2 - 2*lmMst1*lmMst2 + pow2(lmMst1) + pow2(lmMst2)))*pow4(
        Mst1))/pow4(Mst2)))/Mgl) + (16*Tbeta*pow2(Mst1)*pow3(Mt)*(Mgl*(13 + 2*
        lmMst1*(6 + 3*lmMst2 - 2*lmMt) + 2*lmMt + 2*lmMst2*(-7 + 2*lmMt) - 6*
        pow2(lmMst2))*pow4(Mst1) + Dmglst2*(11 - 6*lmMst1*lmMst2 - 8*lmMst2*(-5
        + lmMt) + 8*lmMst1*(-4 + lmMt) - 8*lmMt + 6*pow2(lmMst2))*pow4(Mst1) +
        2*Dmglst2*((7 + 7*lmMst2 + lmMst1*(-5 + lmMt) - 2*lmMt - lmMst2*lmMt)*
        pow2(Mst1)*pow2(Mst2) + (5 + lmMst1*(-1 + lmMst2) + lmMst2 - pow2(
        lmMst2))*pow4(Mst2)) + 2*Mgl*((4 + lmMst1*(3 + 2*lmMst2 - lmMt) +
        lmMst2*(-4 + lmMt) + lmMt - 2*pow2(lmMst2))*pow2(Mst1)*pow2(Mst2) + (2
        + lmMst1 + (-2 + lmMst1)*lmMst2 + lmMt - pow2(lmMst2))*pow4(Mst2))) +
        Tbeta*pow3(Mst2)*pow3(s2t)*(-4*Dmglst2*(2*(1 + 2*lmMst1 - lmMst2)*pow2(
        Mst2)*pow4(Mst1) + 2*(1 + lmMst1 - lmMst2 + 2*lmMst1*lmMst2 - 2*pow2(
        lmMst2))*pow2(Mst1)*pow4(Mst2) + (1 + 2*lmMst1 - 2*lmMst2)*pow6(Mst1) -
        2*lmMst2*pow6(Mst2)) + Mgl*(2*(-16 + 6*lmMst2 - 2*lmMst1*(8 + 5*lmMst2)
        + 3*pow2(lmMst1) + 7*pow2(lmMst2))*pow2(Mst2)*pow4(Mst1) - 2*(-12 - 18*
        lmMst2 + 2*lmMst1*(5 + 3*lmMst2) + pow2(lmMst1) - 7*pow2(lmMst2))*pow2(
        Mst1)*pow4(Mst2) + (3 + lmMst1*(2 - 32*lmMst2) - 2*lmMst2 + 16*(pow2(
        lmMst1) + pow2(lmMst2)))*pow6(Mst1) + 4*(1 + lmMst2)*pow6(Mst2))) - 2*
        Mt*MuSUSY*s2t*(2*Mgl*(-(Mst2*s2t*(-14 - 20*lmMst2 + 2*lmMst1*(5 + 3*
        lmMst2) + pow2(lmMst1) - 7*pow2(lmMst2))) + 8*Mt*(4 + 3*lmMst2 -
        lmMst1*(1 + lmMst2) + pow2(lmMst2)))*pow2(Mst1)*pow4(Mst2) + 4*Dmglst2*
        (Mst2*s2t*(-5 + 8*lmMst2 - 4*lmMst1*(2 + lmMst2) + 4*pow2(lmMst2)) +
        Mt*(65 + lmMst1*(34 - 20*lmMst2) - 26*lmMst2 + 20*pow2(lmMst2)))*pow6(
        Mst1) + Mgl*(Mst2*s2t*(-1 + 50*lmMst2 - 2*lmMst1*(25 + 32*lmMst2) + 20*
        pow2(lmMst1) + 44*pow2(lmMst2)) + Mt*(84 + 152*lmMst2 - 40*lmMst1*(3 +
        2*lmMst2) + 80*pow2(lmMst2)))*pow6(Mst1) + 8*Dmglst2*((Mst2*s2t*(-2 +
        3*lmMst2 - lmMst1*(3 + 2*lmMst2) + 2*pow2(lmMst2)) + Mt*(22 + lmMst1*(8
        - 6*lmMst2) - 4*lmMst2 + 6*pow2(lmMst2)))*pow2(Mst2)*pow4(Mst1) + (-(
        Mst2*s2t*(1 + lmMst1 - 2*lmMst2 + 2*lmMst1*lmMst2 - 2*pow2(lmMst2))) +
        2*Mt*(6 + lmMst1 + lmMst2 - lmMst1*lmMst2 + pow2(lmMst2)))*pow2(Mst1)*
        pow4(Mst2) + lmMst2*s2t*pow7(Mst2)) + 4*Mgl*((4*Mt*(5 + 6*lmMst2 -
        lmMst1*(4 + 3*lmMst2) + 3*pow2(lmMst2)) + Mst2*s2t*(-1 + 13*lmMst2 -
        lmMst1*(13 + 8*lmMst2) + pow2(lmMst1) + 7*pow2(lmMst2)))*pow2(Mst2)*
        pow4(Mst1) + (1 + lmMst2)*s2t*pow7(Mst2))))/(Mgl*Tbeta*pow2(Mst1)*pow5(
        Mst2)))))/48. + (Al4p*z2*(9720*xDmglst2*pow2(Dmglst2)*((-2*MuSUSY*
        twoLoopFlag*(2*Mt*pow2(Mst1)*pow2(Mst2)*(-42*Mt*MuSUSY*s2t + 4*Tbeta*
        pow2(Mt) + 3*Mst2*(MuSUSY + 5*Mst2*Tbeta)*pow2(s2t)) + 2*Mt*(-62*Mt*
        MuSUSY*s2t + 16*Tbeta*pow2(Mt) + 3*Mst2*(MuSUSY + 5*Mst2*Tbeta)*pow2(
        s2t))*pow4(Mst1) + (-44*MuSUSY*s2t*pow2(Mt) + 3*Mst2*Mt*(2*MuSUSY + 11*
        Mst2*Tbeta)*pow2(s2t) + 8*Tbeta*pow3(Mt) - 3*Tbeta*pow3(Mst2)*pow3(s2t)
        )*pow4(Mst2)))/(Tbeta*pow5(Mst2)) + 3*Al4p*threeLoopFlag*((Mt*((Mt*
        pow2(MuSUSY)*((2*(122051 + 1800*lmMst1 + 274680*lmMst2)*Mst2*Mt*pow2(
        Mst1))/405. + (32*(617 + 2040*lmMst2)*Mt*pow3(Mst2))/135. + (
        4131.092729766804 - (3304*lmMst1)/3. + (13736*lmMst2)/9.)*s2t*pow4(
        Mst1) + s2t*pow4(Mst2)*(1108.3172839506174 + (532*lmMst1)/9. + (884*
        lmMst2)/3. - (160*Dmsqst2*(2 + 3*lmMst1 - 3*lmMst2)*(pow3(Dmsqst2) +
        pow6(Msq)))/(3.*pow8(Msq))) + s2t*pow2(Mst1)*pow2(Mst2)*(
        1039.638683127572 - (1160*lmMst1)/9. + (7256*lmMst2)/9. + (160*Dmsqst2*
        (1 - 9*lmMst1 + 9*lmMst2)*(pow3(Dmsqst2) + pow6(Msq)))/(3.*pow8(Msq))))
        )/pow5(Mst2) - pow2(s2t)*(pow2(MuSUSY)*(109.79135802469136 + (368*
        lmMst1)/9. + (148*lmMst2)/9. + ((1167.8951989026064 - (1520*lmMst1)/9.
         + (1864*lmMst2)/
        9.)*pow4(Mst1))/pow4(Mst2) + (pow2(Mst1)*(87.33744855967078 + (490*
        lmMst1)/9. + (26*lmMst2)/9. - (40*Dmsqst2*(3 + 8*lmMst1 - 8*lmMst2)*(
        pow3(Dmsqst2) + pow6(Msq)))/(3.*pow8(Msq))))/pow2(Mst2) + (-40*Dmsqst2*
        (2 + lmMst1 - lmMst2)*(pow3(Dmsqst2) + pow6(Msq)) - (8*pow2(Mst2)*(5*
        pow4(Dmsqst2) + 5*Dmsqst2*pow6(Msq) - 2*pow8(Msq)))/(3.*pow2(Mst1)))/
        pow8(Msq)) + (10*Dmsqst2*(-93 + 10*lmMst1 - 10*lmMst2)*(-1 + pow2(
        Sbeta))*pow2(Sbeta)*pow4(Mst1)*(pow3(Dmsqst2) + pow6(Msq)))/(3.*pow2(
        Mst2)*pow8(Msq)))))/Tbeta + MuSUSY*(-(pow3(s2t)*((-4*pow4(Mst2)*(2 - (
        5*Dmsqst2*(pow3(Dmsqst2) + pow6(Msq)))/pow8(Msq)))/(3.*pow2(Mst1)) - (
        pow4(Mst1)*(432.76982167352537 - (398*lmMst1)/9. + (398*lmMst2)/9. + (
        130*Dmsqst2*(pow3(Dmsqst2) + pow6(Msq)))/(3.*pow8(Msq))))/pow2(Mst2) +
        pow2(Mst1)*(11.226954732510288 - (61*lmMst1)/9. + (61*lmMst2)/9. + (20*
        Dmsqst2*(-3 + 5*lmMst1 - 5*lmMst2)*(pow3(Dmsqst2) + pow6(Msq)))/(3.*
        pow8(Msq))) - pow2(Mst2)*(52.22901234567901 + (184*lmMst1)/9. + (74*
        lmMst2)/9. - (20*Dmsqst2*(5 + 3*lmMst1 - 3*lmMst2)*(pow3(Dmsqst2) +
        pow6(Msq)))/(3.*pow8(Msq))))) - Mt*pow2(s2t)*(Mst2*(831.2379629629629 +
        (133*lmMst1)/3. + 221*lmMst2 - (40*Dmsqst2*(2 + 3*lmMst1 - 3*lmMst2)*(
        pow3(Dmsqst2) + pow6(Msq)))/pow8(Msq)) - (pow2(Mst1)*(51.50895061728395
         + 141*lmMst1 - (1151*lmMst2)/3. + (120*Dmsqst2*(-1 + 2*lmMst1 - 2*
        lmMst2)*(pow3(Dmsqst2) + pow6(Msq)))/pow8(Msq)))/Mst2 + (pow4(Mst1)*(
        1749.093621399177 - 498*lmMst1 + (1790*lmMst2)/3. + (20*Dmsqst2*(17 -
        12*lmMst1 + 12*lmMst2)*(pow3(Dmsqst2) + pow6(Msq)))/pow8(Msq)))/pow3(
        Mst2)) + (pow3(Mt)*((7387.187379972565 - (2456*lmMst1)/3. + (13400*
        lmMst2)/9. - (3088*lmMt)/9.)*pow4(Mst1) - pow2(Mst1)*pow2(Mst2)*(
        3039.4318636096414 - (3896*lmMst1)/9. + (9320*lmMst2)/9. - (1456*lmMt)/
        3. + (1600*Dmsqst2*(-3 + lmMst1 - 2*lmMst2 + lmMt)*(pow3(Dmsqst2) +
        pow6(Msq)))/(3.*pow8(Msq))) - pow4(Mst2)*(433.32081128747797 - (532*
        lmMst1)/9. + (2236*lmMst2)/9. - 72*lmMt + (160*Dmsqst2*(-8 + 3*lmMst1 -
        6*lmMst2 + 3*lmMt)*(pow3(Dmsqst2) + pow6(Msq)))/(3.*pow8(Msq)))))/pow5(
        Mst2) - s2t*pow2(Mt)*(401.86349206349206 - (52*lmMst1)/9. + (6580*
        lmMst2)/9. + (2*(510139 + 44280*lmMst1 + 76680*lmMst2)*pow4(Mst1))/(
        405.*pow4(Mst2)) + (pow2(Mst1)*(1838.5989417989417 - (100*lmMst1)/3. +
        1356*lmMst2 - (240*Dmsqst2*(pow3(Dmsqst2) + pow6(Msq)))/pow8(Msq)))/
        pow2(Mst2) - (400*Dmsqst2*(pow3(Dmsqst2) + pow6(Msq)))/(3.*pow8(Msq)) +
        (16*pow2(Mst2)*(-5*pow4(Dmsqst2) - 5*Dmsqst2*pow6(Msq) + 2*pow8(Msq)))/
        (3.*pow2(Mst1)*pow8(Msq)))))) + (9720*MuSUSY*twoLoopFlag*(4*Mt*
        xDmglst2*xMst*pow2(Dmglst2)*(82*Mt*MuSUSY*s2t - 56*Tbeta*pow2(Mt) - 3*
        Mst2*(MuSUSY + 5*Mst2*Tbeta)*pow2(s2t))*pow6(Mst1) + pow2(Mgl)*(-(pow2(
        Mst1)*(-8*MuSUSY*s2t*pow2(Mt) + 32*Tbeta*pow3(Mt) + Tbeta*pow3(Mst2)*
        pow3(s2t))*pow4(Mst2)) + 8*pow2(Mt)*((MuSUSY*s2t - 6*Mt*Tbeta)*pow2(
        Mst2)*pow4(Mst1) + (MuSUSY*s2t - 8*Mt*Tbeta)*xMst*pow6(Mst1)) + (8*
        MuSUSY*s2t*pow2(Mt) - 2*Mst2*Mt*(MuSUSY + 3*Mst2*Tbeta)*pow2(s2t) - 16*
        Tbeta*pow3(Mt) + Tbeta*pow3(Mst2)*pow3(s2t))*pow6(Mst2)) + 2*Dmglst2*
        Mgl*(4*Mt*(pow2(Mst2)*(13*Mt*MuSUSY*s2t + 6*Tbeta*pow2(Mt) - Mst2*(
        MuSUSY + 3*Mst2*Tbeta)*pow2(s2t))*pow4(Mst1) + s2t*(9*Mt*MuSUSY - Mst2*
        s2t*(MuSUSY + 3*Mst2*Tbeta))*pow2(Mst1)*pow4(Mst2) + xMst*(17*Mt*
        MuSUSY*s2t + 16*Tbeta*pow2(Mt) - Mst2*(MuSUSY + 3*Mst2*Tbeta)*pow2(s2t)
        )*pow6(Mst1)) + (20*MuSUSY*s2t*pow2(Mt) - Mst2*Mt*(4*MuSUSY + 15*Mst2*
        Tbeta)*pow2(s2t) - 8*Tbeta*pow3(Mt) + 2*Tbeta*pow3(Mst2)*pow3(s2t))*
        pow6(Mst2))) + (Al4p*Mst2*threeLoopFlag*(300*xDmsqst2*pow2(Dmsqst2)*
        pow2(Mst2)*pow4(Msq)*(-324*Mst2*xDmglst2*pow2(Dmglst2)*(-2*MuSUSY*(
        pow2(Mst1)*pow2(Mst2)*(-4*s2t*((4 + 6*lmMst1 - 6*lmMst2)*MuSUSY - 5*
        Mst2*Tbeta)*pow2(Mt) + 6*Mst2*Mt*((2 + lmMst1 - lmMst2)*MuSUSY + (2 +
        3*lmMst1 - 3*lmMst2)*Mst2*Tbeta)*pow2(s2t) + Tbeta*(8*(8 + 6*lmMst2 -
        3*(lmMst1 + lmMt))*pow3(Mt) - (5 + 3*lmMst1 - 3*lmMst2)*pow3(Mst2)*
        pow3(s2t))) + (4*s2t*((2 - 18*lmMst1 + 18*lmMst2)*MuSUSY + 9*Mst2*
        Tbeta)*pow2(Mt) + 2*Mst2*Mt*((3 + 8*lmMst1 - 8*lmMst2)*MuSUSY + 9*(-1 +
        2*lmMst1 - 2*lmMst2)*Mst2*Tbeta)*pow2(s2t) - 80*(-3 + lmMst1 - 2*lmMst2
        + lmMt)*Tbeta*pow3(Mt) + (3 - 5*lmMst1 + 5*lmMst2)*Tbeta*pow3(Mst2)*
        pow3(s2t))*pow4(Mst1) + s2t*(2*Mt*MuSUSY*s2t + 4*Tbeta*pow2(Mt) -
        Tbeta*pow2(Mst2)*pow2(s2t))*pow5(Mst2)) + pow2(s2t)*(6*(17 - 12*lmMst1
        + 12*lmMst2)*Mt*MuSUSY*Tbeta - 13*Mst2*MuSUSY*s2t*Tbeta + (-93 + 10*
        lmMst1 - 10*lmMst2)*Mst2*Mt*(-1 + pow2(Sbeta))*pow2(Sbeta))*pow6(Mst1))
        + Mgl*(s2t*(s2t*(81*pow2(Mst2)*pow2(Sbeta)*(-2*(11 + 10*lmMst1 - 10*
        lmMst2)*Mgl*Mst2*s2t*(-1 + pow2(Sbeta)) + (Dmglst2*(372 - 40*lmMst1 +
        40*lmMst2) + 3*(-17 + 2*lmMst1 - 2*lmMst2)*Mgl)*Mt*pow2(Sbeta))*pow6(
        Mst1) + Mt*(81*(4*Dmglst2*(-93 + 10*lmMst1 - 10*lmMst2) + 3*(17 - 2*
        lmMst1 + 2*lmMst2)*Mgl)*pow2(Mst2)*pow2(Sbeta)*pow6(Mst1) + 2*pow2(
        MuSUSY)*(648*Dmglst2*pow2(Mst2)*(3*(2 + lmMst1 - lmMst2)*pow2(Mst1)*
        pow2(Mst2) + (3 + 8*lmMst1 - 8*lmMst2)*pow4(Mst1) + pow4(Mst2)) + Mgl*(
        3*(1607 - 432*lmMst1 + 432*lmMst2)*pow2(Mst2)*pow4(Mst1) - 18*(25 + 9*
        lmMst1 - 9*lmMst2)*pow2(Mst1)*pow4(Mst2) + (20701 - 648*lmMst1 + 648*
        lmMst2)*pow6(Mst1) - 486*pow6(Mst2))))) + 648*Mt*pow2(MuSUSY)*(8*
        Dmglst2*Mst2*Mt*pow2(Mst1)*((1 - 9*lmMst1 + 9*lmMst2)*pow2(Mst1) + (-2
        - 3*lmMst1 + 3*lmMst2)*pow2(Mst2)) + Mgl*((4*(7 + 3*lmMst1 - 3*lmMst2)*
        Mt + (-1 + 2*lmMst1 - 2*lmMst2)*Mst2*s2t*shiftst1)*pow2(Mst1)*pow3(
        Mst2) + 2*Mst2*(4*(5 + 3*lmMst1 - 3*lmMst2)*Mt + (lmMst1 - lmMst2)*
        Mst2*s2t*shiftst1)*pow4(Mst1) + 2*(lmMst1 - lmMst2)*s2t*shiftst1*pow6(
        Mst1) - s2t*shiftst1*pow6(Mst2)))) + MuSUSY*Tbeta*(Mgl*(-18*pow2(Mst1)*
        pow3(Mst2)*(2*Mst2*s2t*(25 - 36*shiftst1)*pow2(Mt) + 108*(7 + 3*lmMst1
        - 3*lmMst2)*Mt*pow2(Mst2)*pow2(s2t) + 144*(2 + 6*lmMst2 - 3*(lmMst1 +
        lmMt))*pow3(Mt) + (2 + lmMst2*(9 - 36*shiftst1) + 9*lmMst1*(-1 + 4*
        shiftst1))*pow3(Mst2)*pow3(s2t)) - 3*Mst2*(530*Mst2*s2t*pow2(Mt) +
        1296*(4 + 3*lmMst1 - 3*lmMst2)*Mt*pow2(Mst2)*pow2(s2t) - 3456*(-1 +
        lmMst1 - 2*lmMst2 + lmMt)*pow3(Mt) + (1757 - 378*lmMst1 + 378*lmMst2 +
        108*shiftst1)*pow3(Mst2)*pow3(s2t))*pow4(Mst1) + s2t*((-1944*(5 + 4*
        lmMst1 - 4*lmMst2)*Mst2*Mt*s2t + 8114*pow2(Mt) - 6648*pow2(Mst2)*pow2(
        s2t))*pow6(Mst1) + 162*(3 + 2*shiftst1)*(-4*pow2(Mt) + pow2(Mst2)*pow2(
        s2t))*pow6(Mst2))) + 324*Dmglst2*Mst2*(2*(pow2(Mst1)*pow2(Mst2)*(20*
        Mst2*s2t*pow2(Mt) + 6*(2 + 3*lmMst1 - 3*lmMst2)*Mt*pow2(Mst2)*pow2(s2t)
        + 8*(8 + 6*lmMst2 - 3*(lmMst1 + lmMt))*pow3(Mt) - (5 + 3*lmMst1 - 3*
        lmMst2)*pow3(Mst2)*pow3(s2t)) + (36*Mst2*s2t*pow2(Mt) + 18*(-1 + 2*
        lmMst1 - 2*lmMst2)*Mt*pow2(Mst2)*pow2(s2t) - 80*(-3 + lmMst1 - 2*lmMst2
        + lmMt)*pow3(Mt) + (3 - 5*lmMst1 + 5*lmMst2)*pow3(Mst2)*pow3(s2t))*
        pow4(Mst1)) + 8*s2t*pow2(Mt)*pow5(Mst2) + (6*(-17 + 12*lmMst1 - 12*
        lmMst2)*Mt + 13*Mst2*s2t)*pow2(s2t)*pow6(Mst1) - 2*pow3(s2t)*pow7(Mst2)
        )))) + 3240*Mst2*MuSUSY*xDR2DRMOD*(Mgl*(Tbeta*(64*Mt*pow2(Mst1)*(-3*(
        Dmglst2 + 3*Dmglst2*lmMst2 + Mgl + lmMst2*Mgl)*pow2(Mst2)*pow2(s2t)*(
        pow2(Mst1)*pow2(Mst2) + pow4(Mst1) + pow4(Mst2)) + 4*pow2(Mt)*(-(
        Dmglst2*(7 + lmMst2)*pow2(Mst1)*pow2(Mst2)) - 2*Dmglst2*(11 + 5*lmMst2)
        *pow4(Mst1) + Dmglst2*(-1 + lmMst2)*pow4(Mst2) + (1 + lmMst2)*Mgl*(3*
        pow2(Mst1)*pow2(Mst2) + 6*pow4(Mst1) + pow4(Mst2)))) + (16*(7*Dmglst2 -
        13*Mgl)*s2t*(pow2(Mst1) - pow2(Mst2))*pow2(Mt) - pow3(s2t)*(8*(Dmglst2*
        (7*lmMst1 - 15*lmMst2) + (-4 - 13*lmMst1 + 9*lmMst2)*Mgl)*pow2(Mst1)*
        pow2(Mst2) + 4*(7*Dmglst2 - 13*Mgl)*(pow4(Mst1) - pow4(Mst2))))*pow5(
        Mst2)) + Mt*MuSUSY*s2t*(256*(Dmglst2 + 3*Dmglst2*lmMst2 + Mgl + lmMst2*
        Mgl)*Mt*pow2(Mst1)*(2*pow2(Mst1)*pow2(Mst2) + 3*pow4(Mst1) + pow4(Mst2)
        ) - 8*Mst2*s2t*(-2*(Dmglst2*(7*lmMst1 - 15*lmMst2) + (-4 - 13*lmMst1 +
        9*lmMst2)*Mgl)*(pow2(Mst1) + pow2(Mst2))*pow4(Mst1) + (Dmglst2*(7 - 14*
        lmMst1 + 30*lmMst2) + (-5 + 26*lmMst1 - 18*lmMst2)*Mgl)*pow2(Mst1)*
        pow4(Mst2) + (7*Dmglst2 - 13*Mgl)*pow6(Mst2)))) - xDmglst2*pow2(
        Dmglst2)*(Tbeta*(64*Mt*pow2(Mst1)*(3*pow2(Mst2)*pow2(s2t)*((-1 + 6*
        lmMst2)*pow2(Mst2)*(pow2(Mst1) + pow2(Mst2)) + (1 + 3*lmMst2)*pow4(
        Mst1)) + 4*pow2(Mt)*(-(lmMst2*pow2(Mst1)*pow2(Mst2)) + 2*(11 + 5*
        lmMst2)*pow4(Mst1) - (-2 + lmMst2)*pow4(Mst2))) + (-432*s2t*(pow2(Mst1)
        - pow2(Mst2))*pow2(Mt) + 9*pow3(s2t)*((8*(8 + 27*lmMst1 - 39*lmMst2)*
        pow2(Mst1)*pow2(Mst2))/9. + 12*pow4(Mst1) - 12*pow4(Mst2)))*pow5(Mst2))
        + 8*Mt*MuSUSY*s2t*(2*(-32*(-1 + 6*lmMst2)*Mt + (-8 - 27*lmMst1 + 39*
        lmMst2)*Mst2*s2t)*pow2(Mst2)*pow4(Mst1) + (32*(1 - 6*lmMst2)*Mt + (11 -
        54*lmMst1 + 78*lmMst2)*Mst2*s2t)*pow2(Mst1)*pow4(Mst2) - 2*(48*(1 + 3*
        lmMst2)*Mt + (7*lmMst1 - 15*lmMst2)*Mst2*s2t)*pow6(Mst1) + 27*s2t*pow7(
        Mst2))))*pow8(Msq) + Mgl*(Mst2*MuSUSY*Tbeta*(pow3(Mst2)*pow3(s2t)*(-
        300*Dmsqst2*Mgl*pow6(Msq)*(18*(469 - 36*lmMst1 + 36*lmMst2)*pow2(Mst2)*
        pow4(Mst1) + 1755*pow2(Mst1)*pow4(Mst2) + 10828*pow6(Mst1) - 324*pow6(
        Mst2)) + 300*Mgl*pow4(Dmsqst2)*(9*(119 + 234*lmMst1 - 234*lmMst2)*pow2(
        Mst2)*pow4(Mst1) + 486*(7 + lmMst1 - lmMst2)*pow2(Mst1)*pow4(Mst2) +
        1712*pow6(Mst1) + 810*pow6(Mst2)) + 97200*Dmglst2*Dmsqst2*(pow3(
        Dmsqst2) + pow6(Msq))*(2*(3 - 5*lmMst1 + 5*lmMst2)*pow2(Mst2)*pow4(
        Mst1) + 13*pow6(Mst1) - 2*((5 + 3*lmMst1 - 3*lmMst2)*pow2(Mst1)*pow4(
        Mst2) + pow6(Mst2))) + 4*Dmglst2*(27*(32663 - 7620*lmMst1 + 7620*
        lmMst2)*pow2(Mst2)*pow4(Mst1) + 135*(1799 - 648*lmMst1 + 1680*lmMst2)*
        pow2(Mst1)*pow4(Mst2) + 4*(788723 - 80595*lmMst1 + 80595*lmMst2)*pow6(
        Mst1) - 22680*pow6(Mst2))*pow8(Msq) - 3*Mgl*(12*(224837 - 4680*lmMst1 +
        8370*lmMst2)*pow2(Mst2)*pow4(Mst1) + 90*(5825 + 564*lmMst1 - 1104*
        lmMst2)*pow2(Mst1)*pow4(Mst2) + (3447379 - 76140*lmMst1 + 76140*lmMst2)
        *pow6(Mst1) - 20520*pow6(Mst2))*pow8(Msq)) + s2t*(6*Mst2*pow2(Mt)*(100*
        (pow4(Dmsqst2)*(1296*Dmglst2*pow2(Mst2)*(5*pow2(Mst1)*pow2(Mst2) + 9*
        pow4(Mst1) + pow4(Mst2)) + Mgl*(-8883*pow2(Mst2)*pow4(Mst1) - 4212*
        pow2(Mst1)*pow4(Mst2) + 4057*pow6(Mst1) - 1620*pow6(Mst2))) + Dmsqst2*
        pow6(Msq)*(1296*Dmglst2*pow2(Mst2)*(5*pow2(Mst1)*pow2(Mst2) + 9*pow4(
        Mst1) + pow4(Mst2)) + Mgl*(3249*pow2(Mst2)*pow4(Mst1) + 1431*pow2(Mst1)
        *pow4(Mst2) + 4057*pow6(Mst1) - 648*pow6(Mst2)))) + (Mgl*(3*(623599 +
        65160*lmMst1 - 134280*lmMst2)*pow2(Mst2)*pow4(Mst1) + 9*(1367 + 13560*
        lmMst1 - 36600*lmMst2)*pow2(Mst1)*pow4(Mst2) + (5161826 + 406080*lmMst1
        - 613440*lmMst2)*pow6(Mst1) - 41040*pow6(Mst2)) - 24*Dmglst2*(6*(19579
        + 1770*lmMst1 + 12630*lmMst2)*pow2(Mst2)*pow4(Mst1) + 6*(4429 - 270*
        lmMst1 + 8910*lmMst2)*pow2(Mst1)*pow4(Mst2) + (510139 + 44280*lmMst1 +
        76680*lmMst2)*pow6(Mst1) - 2520*pow6(Mst2)))*pow8(Msq)) + 9720*Mgl*
        pow5(Mst2)*(-(shiftst3*pow2(Mst2)*(4*pow2(Mt) + ((1 + lmMst1 - lmMst2)*
        pow2(Mst1) - pow2(Mst2))*pow2(s2t))*pow8(Msq)) - 10*shiftst1*(-4*pow2(
        Mst1)*pow2(Mt) + 4*pow2(Mst2)*pow2(Mt) + 2*(lmMst1 - lmMst2)*pow2(Mst1)
        *pow2(Mst2)*pow2(s2t) + pow2(s2t)*(pow4(Mst1) - pow4(Mst2)))*(pow4(
        Dmsqst2) + Dmsqst2*pow6(Msq) + pow8(Msq)))) + pow2(Mst1)*(-9*Mt*pow2(
        Mst2)*pow2(s2t)*(2*pow4(Mst1)*(32400*Dmsqst2*(Dmglst2*(17 - 12*lmMst1 +
        12*lmMst2) + (5 + 4*lmMst1 - 4*lmMst2)*Mgl)*(pow3(Dmsqst2) + pow6(Msq))
        + (Dmglst2*(2.8335316666666665e6 - 806760*lmMst1 + 966600*lmMst2) + (
        30643 + 217080*lmMst1 - 212760*lmMst2)*Mgl)*pow8(Msq)) + 15*pow2(Mst1)*
        pow2(Mst2)*(8640*((6 + 5*lmMst1 - 5*lmMst2)*Mgl*pow4(Dmsqst2) +
        Dmsqst2*(3 + 2*lmMst1 - 2*lmMst2)*Mgl*pow6(Msq) - 3*Dmglst2*Dmsqst2*(-1
        + 2*lmMst1 - 2*lmMst2)*(pow3(Dmsqst2) + pow6(Msq))) + (Dmglst2*(
        187064.6 - 55800*lmMst1 + 77112*lmMst2) + (11075 + 17928*lmMst1 -
        17352*lmMst2)*Mgl)*pow8(Msq)) - 135*pow4(Mst2)*(480*(Dmglst2*(4 + 6*
        lmMst1 - 6*lmMst2) + (-13 - 5*lmMst1 + 5*lmMst2)*Mgl)*pow4(Dmsqst2) +
        960*Dmsqst2*(Dmglst2*(2 + 3*lmMst1 - 3*lmMst2) + (-2 - lmMst1 + lmMst2)
        *Mgl)*pow6(Msq) + (Dmglst2*(-11596.466666666667 + 1672*lmMst1 - 4584*
        lmMst2) + (-2269 - 664*lmMst1 + 56*lmMst2)*Mgl)*pow8(Msq))) + 12*pow3(
        Mt)*((4*(Dmglst2*(13463149 - 1492020*lmMst1 + 2713500*lmMst2 - 625320*
        lmMt) + 3*(-1077991 + 184140*lmMst1 - 400140*lmMst2 + 8640*lmMt)*Mgl)*
        pow4(Mst1)*pow8(Msq))/3. + 24*pow2(Mst1)*pow2(Mst2)*(10800*Dmsqst2*(-5*
        Dmglst2*(-3 + lmMst1 - 2*lmMst2 + lmMt) + (-1 + lmMst1 - 2*lmMst2 +
        lmMt)*Mgl)*(pow3(Dmsqst2) + pow6(Msq)) + (3*Dmglst2*(67843 - 10050*
        lmMst1 + 17490*lmMst2 - 7560*lmMt) + 2*(-41936 + 7245*lmMst1 - 19665*
        lmMst2 + 720*lmMt)*Mgl)*pow8(Msq)) - 9*pow4(Mst2)*(-15*Dmglst2*(960*
        Dmsqst2*(8 + 6*lmMst2 - 3*(lmMst1 + lmMt))*(pow3(Dmsqst2) + pow6(Msq))
        + (-1672*lmMst1 + 1448*lmMst2 + 59*(55 - 32*lmMt))*pow8(Msq)) + Mgl*(
        7200*(6 + 10*lmMst2 - 5*(lmMst1 + lmMt))*pow4(Dmsqst2) - 14400*Dmsqst2*
        (lmMst1 - 2*lmMst2 + lmMt)*pow6(Msq) + (66773 - 9960*lmMst1 + 45480*
        lmMst2 - 3840*lmMt)*pow8(Msq)))))) + 2*(pow2(s2t)*(12150*Dmsqst2*pow2(
        Sbeta)*pow4(Mst2)*((-6*(11 + 10*lmMst1 - 10*lmMst2)*Mgl*Mst2*s2t*(-1 +
        pow2(Sbeta)) + (Dmglst2*(372 - 40*lmMst1 + 40*lmMst2) + 9*(-17 + 2*
        lmMst1 - 2*lmMst2)*Mgl)*Mt*pow2(Sbeta))*pow3(Dmsqst2) + 4*Dmglst2*(93 -
        10*lmMst1 + 10*lmMst2)*Mt*pow2(Sbeta)*pow6(Msq))*pow6(Mst1) + Mt*pow2(
        Mst2)*(12150*pow2(Mst2)*pow2(Sbeta)*(9*(17 - 2*lmMst1 + 2*lmMst2)*Mgl*
        pow4(Dmsqst2) + 4*Dmglst2*Dmsqst2*(-93 + 10*lmMst1 - 10*lmMst2)*(pow3(
        Dmsqst2) + pow6(Msq)))*pow6(Mst1) + pow2(MuSUSY)*(300*(Dmsqst2*pow6(
        Msq)*(648*Dmglst2*pow2(Mst2)*(3*(2 + lmMst1 - lmMst2)*pow2(Mst1)*pow2(
        Mst2) + (3 + 8*lmMst1 - 8*lmMst2)*pow4(Mst1) + pow4(Mst2)) + Mgl*(9*(
        1097 - 72*lmMst1 + 72*lmMst2)*pow2(Mst2)*pow4(Mst1) + 1431*pow2(Mst1)*
        pow4(Mst2) + (20701 - 648*lmMst1 + 648*lmMst2)*pow6(Mst1) - 324*pow6(
        Mst2))) + pow4(Dmsqst2)*(648*Dmglst2*pow2(Mst2)*(3*(2 + lmMst1 -
        lmMst2)*pow2(Mst1)*pow2(Mst2) + (3 + 8*lmMst1 - 8*lmMst2)*pow4(Mst1) +
        pow4(Mst2)) - Mgl*(9*(587 + 288*lmMst1 - 288*lmMst2)*pow2(Mst2)*pow4(
        Mst1) + 162*(26 + 3*lmMst1 - 3*lmMst2)*pow2(Mst1)*pow4(Mst2) + (-20701
        + 648*lmMst1 - 648*lmMst2)*pow6(Mst1) + 810*pow6(Mst2)))) + (-4*
        Dmglst2*(162*(6803 - 1810*lmMst1 + 2670*lmMst2)*pow2(Mst2)*pow4(Mst1) +
        135*(1631 - 648*lmMst1 + 1680*lmMst2)*pow2(Mst1)*pow4(Mst2) + (4256978
        - 615600*lmMst1 + 754920*lmMst2)*pow6(Mst1) - 22680*pow6(Mst2)) + 3*
        Mgl*(6*(533629 - 900*lmMst1 + 180*lmMst2)*pow2(Mst2)*pow4(Mst1) + 90*(
        5597 + 564*lmMst1 - 1104*lmMst2)*pow2(Mst1)*pow4(Mst2) + (6649153 -
        81540*lmMst1 + 77220*lmMst2)*pow6(Mst1) - 20520*pow6(Mst2)))*pow8(Msq))
        )) + 6*Mt*pow2(MuSUSY)*(-1620*Mgl*pow2(Mst2)*pow2(s2t)*(shiftst3*(2*(1
        - 2*lmMst1 + 2*lmMst2)*pow2(Mst2)*pow4(Mst1) + 2*(1 - lmMst1 + lmMst2)*
        pow2(Mst1)*pow4(Mst2) + (2 - 6*lmMst1 + 6*lmMst2)*pow6(Mst1) + pow6(
        Mst2))*pow8(Msq) + 10*shiftst1*(pow2(Mst1)*(pow4(Mst2) - 2*(lmMst1 -
        lmMst2)*(pow2(Mst1)*pow2(Mst2) + pow4(Mst1) + pow4(Mst2))) + pow6(Mst2)
        )*(pow4(Dmsqst2) + Dmsqst2*pow6(Msq) + pow8(Msq))) + Mt*pow2(Mst1)*(30*
        Mt*(12*Dmglst2*pow2(Mst2)*((3643 - 120*lmMst1 + 3192*lmMst2)*pow2(Mst1)
        + 192*(5 + 6*lmMst2)*pow2(Mst2)) + Mgl*(6*(-353 + 72*lmMst1 + 696*
        lmMst2)*pow2(Mst1)*pow2(Mst2) + (-38401 - 1080*lmMst1 + 7992*lmMst2)*
        pow4(Mst1) + 96*(53 + 24*lmMst2)*pow4(Mst2)))*pow8(Msq) + Mst2*s2t*((2*
        (Dmglst2*(15057833 - 4014360*lmMst1 + 5563080*lmMst2) + 3*(266863 +
        396360*lmMst1 - 346680*lmMst2)*Mgl)*pow4(Mst1)*pow8(Msq))/3. + 9*pow4(
        Mst2)*(7200*(13 + 5*lmMst1 - 5*lmMst2)*Mgl*pow4(Dmsqst2) + 14400*
        Dmsqst2*(2 + lmMst1 - lmMst2)*Mgl*pow6(Msq) - 14400*Dmglst2*Dmsqst2*(2
        + 3*lmMst1 - 3*lmMst2)*(pow3(Dmsqst2) + pow6(Msq)) + Dmglst2*(173947 -
        25080*lmMst1 + 68760*lmMst2)*pow8(Msq) + 15*(2269 + 664*lmMst1 - 56*
        lmMst2)*Mgl*pow8(Msq)) + 30*pow2(Mst1)*pow2(Mst2)*(4320*Dmsqst2*(
        Dmglst2*(1 - 9*lmMst1 + 9*lmMst2) + (5 + 3*lmMst1 - 3*lmMst2)*Mgl)*(
        pow3(Dmsqst2) + pow6(Msq)) + (Dmglst2*(145716.4 - 35424*lmMst1 + 59184*
        lmMst2) + 4*(3937 + 2988*lmMst1 - 2232*lmMst2)*Mgl)*pow8(Msq)))))))))/(
        pow2(Mst1)*pow8(Msq)))/(Tbeta*pow7(Mst2))))/(29160.*pow2(Mgl))) +
        threeLoopFlag*pow2(Al4p)*(-(xDmsqst2*pow2(Dmsqst2)*((-9261000*Mst2*Mt*
        xDmglst2*pow2(Dmglst2)*(MuSUSY*(-4*pow2(Mst1)*pow2(Mst2)*(-6*s2t*(12*(2
        + 5*lmMst1 - 5*lmMst2)*MuSUSY - 41*Mst2*Tbeta)*pow2(Mt) + 18*Mst2*Mt*(
        2*(4 + lmMst1 - lmMst2)*MuSUSY + 3*(2 + 5*lmMst1 - 5*lmMst2)*Mst2*
        Tbeta)*pow2(s2t) - 4*(-101 + 90*lmMst1 - 177*lmMst2 + 87*lmMt)*Tbeta*
        pow3(Mt) - 9*(5 + 2*lmMst1 - 2*lmMst2)*Tbeta*pow3(Mst2)*pow3(s2t)) - (
        36*s2t*(8*(8 - 15*lmMst1 + 15*lmMst2)*MuSUSY + (65 - 2*lmMst1 + 2*
        lmMst2)*Mst2*Tbeta)*pow2(Mt) - 18*Mst2*Mt*((23 - 42*lmMst1 + 42*lmMst2)
        *MuSUSY + 120*(1 - lmMst1 + lmMst2)*Mst2*Tbeta)*pow2(s2t) + Tbeta*(16*(
        497 - 306*lmMst1 + 609*lmMst2 - 303*lmMt)*pow3(Mt) + 9*(55 - 34*lmMst1
        + 34*lmMst2)*pow3(Mst2)*pow3(s2t)))*pow4(Mst1) + 108*s2t*(-2*Mt*MuSUSY*
        s2t - 4*Tbeta*pow2(Mt) + Tbeta*pow2(Mst2)*pow2(s2t))*pow5(Mst2)) +
        pow2(s2t)*(54*(80 - 43*lmMst1 + 43*lmMst2)*Mt*MuSUSY*Tbeta + (-419 +
        57*lmMst1 - 57*lmMst2)*Mst2*MuSUSY*s2t*Tbeta + 3*(-1081 + 165*lmMst1 -
        165*lmMst2)*Mst2*Mt*(-1 + pow2(Sbeta))*pow2(Sbeta))*pow6(Mst1)))/Tbeta
        + Mgl*(Mt*MuSUSY*(-166698000*(2 - lmMst2)*Mgl*s2t*shiftst1*(-4*pow2(
        Mst1)*pow2(Mt) + 4*pow2(Mst2)*pow2(Mt) + 2*(lmMst1 - lmMst2)*pow2(Mst1)
        *pow2(Mst2)*pow2(s2t) + pow2(s2t)*(pow4(Mst1) - pow4(Mst2)))*pow4(Mst2)
        + 2058000*pow2(Mst1)*(18*Mst2*(4*Dmglst2*((497 - 306*lmMst1 + 609*
        lmMst2 - 303*lmMt)*pow2(Mst1) + (101 - 90*lmMst1 + 177*lmMst2 - 87*
        lmMt)*pow2(Mst2)) + 9*Mgl*(4*(1 + 6*lmMst1 - 13*lmMst2 + 7*lmMt)*pow2(
        Mst1) + (16 + 16*lmMst1 - 35*lmMst2 + 19*lmMt)*pow2(Mst2)))*pow3(Mt) +
        6*Mgl*T1ep*pow2(Mst2)*pow3(s2t)*(19*pow2(Mst1)*(pow2(Mst1) + pow2(Mst2)
        ) + 6*pow4(Mst2)) - 4*Mgl*s2t*T1ep*pow2(Mt)*(21*pow2(Mst1)*pow2(Mst2) +
        53*pow4(Mst1) + 18*pow4(Mst2)) - 243*Mst2*Mt*pow2(s2t)*(40*Dmglst2*(1 -
        lmMst1 + lmMst2)*pow2(Mst1)*pow2(Mst2) + Dmglst2*(80 - 43*lmMst1 + 43*
        lmMst2)*pow4(Mst1) + (-4 + 9*lmMst1 - 9*lmMst2)*Mgl*pow4(Mst1) - 4*
        Dmglst2*(2 + 5*lmMst1 - 5*lmMst2)*pow4(Mst2) + 2*Mgl*((1 + 8*lmMst1 -
        8*lmMst2)*pow2(Mst1)*pow2(Mst2) + (5 + 4*lmMst1 - 4*lmMst2)*pow4(Mst2))
        )) + 3*pow2(Mst2)*pow3(s2t)*(3087000*Dmglst2*(9*(55 - 34*lmMst1 + 34*
        lmMst2)*pow2(Mst2)*pow4(Mst1) - 36*(5 + 2*lmMst1 - 2*lmMst2)*pow2(Mst1)
        *pow4(Mst2) + (419 - 57*lmMst1 + 57*lmMst2)*pow6(Mst1) - 108*pow6(Mst2)
        ) + Mgl*(1715*(192439 - 30400*OepS2 - 837000*S2 - 180*lmMst2*(857 +
        3420*S2) + 180*lmMst1*(677 + 180*lmMst2 + 3420*S2) - 16200*(pow2(
        lmMst1) + pow2(lmMst2)))*pow2(Mst2)*pow4(Mst1) + 128625*(8555 - 128*
        OepS2 - 57024*S2 - 72*lmMst2*(29 + 36*S2) + 72*lmMst1*(29 - 12*lmMst2 +
        36*S2) + 864*pow2(lmMst2))*pow2(Mst1)*pow4(Mst2) - 2*(55313478 +
        26068000*OepS2 + 563377500*S2 + 315*lmMst2*(-423553 + 1675800*S2) -
        315*lmMst1*(-423553 + 166320*lmMst2 + 1675800*S2) + 26195400*(pow2(
        lmMst1) + pow2(lmMst2)))*pow6(Mst1) + 27783000*(1 + 2*lmMst2)*pow6(
        Mst2))) + 4*s2t*pow2(Mt)*(27783000*Dmglst2*pow2(Mst2)*(82*pow2(Mst1)*
        pow2(Mst2) + (195 - 6*lmMst1 + 6*lmMst2)*pow4(Mst1) + 36*pow4(Mst2)) -
        Mgl*(5145*(279089 - 5600*OepS2 + 1260*lmMst2*(29 - 90*S2) + 812700*S2 +
        180*lmMst1*(-203 + 90*lmMst2 + 630*S2) - 8100*(pow2(lmMst1) + pow2(
        lmMst2)))*pow2(Mst2)*pow4(Mst1) + 385875*(3655 - 64*OepS2 + 8424*S2 +
        144*lmMst1*(1 + 9*S2) - 144*lmMst2*(4 + 9*S2))*pow2(Mst1)*pow4(Mst2) +
        (1094369501 - 72716000*OepS2 + 2520*lmMst2*(235453 - 584325*S2) +
        5940931500*S2 + 2520*lmMst1*(-235453 + 114345*lmMst2 + 584325*S2) -
        144074700*(pow2(lmMst1) + pow2(lmMst2)))*pow6(Mst1) + 83349000*(1 + 2*
        lmMst2)*pow6(Mst2)))) - (10*Mt*s2t*(-411600*Mt*pow2(MuSUSY)*(648*
        Dmglst2*Mst2*Mt*pow2(Mst1)*((8 - 15*lmMst1 + 15*lmMst2)*pow2(Mst1) + (-
        2 - 5*lmMst1 + 5*lmMst2)*pow2(Mst2)) + Mgl*(9*(36*(5 + 4*lmMst1 - 4*
        lmMst2)*Mt - Mst2*s2t*(9*(-1 + 2*lmMst1 - 2*lmMst2)*(-2 + lmMst2)*
        shiftst1 + 4*T1ep))*pow2(Mst1)*pow3(Mst2) + 6*Mst2*(108*(1 + 3*lmMst1 -
        3*lmMst2)*Mt + Mst2*s2t*(-54*lmMst2*shiftst1 - 25*T1ep + 27*shiftst1*(-
        (lmMst1*(-2 + lmMst2)) + pow2(lmMst2))))*pow4(Mst1) - s2t*(81*(lmMst1 -
        lmMst2)*(-3 + 2*lmMst2)*shiftst1 + 442*T1ep)*pow6(Mst1) + 81*(-2 +
        lmMst2)*s2t*shiftst1*pow6(Mst2))) + s2t*(-694575*pow2(Mst2)*pow2(Sbeta)
        *(-6*(11 + 39*lmMst1 - 39*lmMst2)*Mgl*Mst2*s2t*(-1 + pow2(Sbeta)) + (
        Dmglst2*(4324 - 660*lmMst1 + 660*lmMst2) + 3*(-167 + 5*lmMst1 - 5*
        lmMst2)*Mgl)*Mt*pow2(Sbeta))*pow6(Mst1) + Mt*(694575*(Dmglst2*(4324 -
        660*lmMst1 + 660*lmMst2) + 3*(-167 + 5*lmMst1 - 5*lmMst2)*Mgl)*pow2(
        Mst2)*pow2(Sbeta)*pow6(Mst1) + pow2(MuSUSY)*(16669800*Dmglst2*pow2(
        Mst2)*(-8*(4 + lmMst1 - lmMst2)*pow2(Mst1)*pow2(Mst2) + (23 - 42*lmMst1
        + 42*lmMst2)*pow4(Mst1) - 12*pow4(Mst2)) + Mgl*(-4116*(45*lmMst1*(-1547
        + 180*lmMst2 - 4500*S2) + 45*lmMst2*(1547 + 4500*S2) + 2*(-106283 +
        5000*OepS2 + 639225*S2) + 4050*pow2(lmMst1) - 12150*pow2(lmMst2))*pow2(
        Mst2)*pow4(Mst1) + 77175*(8771 - 128*OepS2 - 57024*S2 - 72*lmMst2*(23 +
        36*S2) + 72*lmMst1*(29 - 12*lmMst2 + 36*S2) + 864*pow2(lmMst2))*pow2(
        Mst1)*pow4(Mst2) + 2*(593331163 - 60642400*OepS2 + 1260*lmMst2*(8263 -
        974610*S2) + 1143733500*S2 + 1260*lmMst1*(-8263 + 71820*lmMst2 +
        974610*S2) - 61916400*pow2(lmMst1) - 28576800*pow2(lmMst2))*pow6(Mst1)
        + 16669800*(1 + 2*lmMst2)*pow6(Mst2)))))))/Tbeta)))/(5.00094e7*pow2(
        Mgl)*pow2(Mst1)*pow4(Msq)*pow4(Mst2)) + ((Mt*((15*(8*MuSUSY*z4*(
        Dmglst2*Mst2*(Mgl*((-16*s2t*(659*MuSUSY + 6813*Mst2*Tbeta)*pow2(Mt) -
        2*Mst2*Mt*(45308*MuSUSY + 15189*Mst2*Tbeta)*pow2(s2t) + 228304*Tbeta*
        pow3(Mt) + 25328*Tbeta*pow3(Mst2)*pow3(s2t))*pow4(Mst1) + 27*pow3(Mst2)
        *(-4*(Mst2*s2t*(271*MuSUSY - 216*Mst2*Tbeta) + Mt*(216*MuSUSY + 253*
        Mst2*Tbeta))*pow2(Mt) + Mt*(-848*MuSUSY + 813*Mst2*Tbeta)*pow2(Mst2)*
        pow2(s2t) + 424*Tbeta*pow3(s2t)*pow4(Mst2)) + 9*Mst2*pow2(Mst1)*(8*
        Mst2*s2t*(-709*MuSUSY + 24*Mst2*Tbeta)*pow2(Mt) + 15*Mt*(-296*MuSUSY +
        121*Mst2*Tbeta)*pow2(Mst2)*pow2(s2t) - 96*(47*MuSUSY - 26*Mst2*Tbeta)*
        pow3(Mt) + 948*Tbeta*pow3(s2t)*pow4(Mst2))) + Dmglst2*xDmglst2*((-16*
        s2t*(659*MuSUSY + 6813*Mst2*Tbeta)*pow2(Mt) - 2*Mst2*Mt*(45308*MuSUSY +
        15189*Mst2*Tbeta)*pow2(s2t) + 228304*Tbeta*pow3(Mt) + 25328*Tbeta*pow3(
        Mst2)*pow3(s2t))*pow4(Mst1) + 9*pow3(Mst2)*(-4*(Mst2*s2t*(1249*MuSUSY -
        846*Mst2*Tbeta) + Mt*(972*MuSUSY + 871*Mst2*Tbeta))*pow2(Mt) + Mt*(-
        8888*MuSUSY + 3747*Mst2*Tbeta)*pow2(Mst2)*pow2(s2t) + 4444*Tbeta*pow3(
        s2t)*pow4(Mst2)) + 3*Mst2*pow2(Mst1)*(-16*(Mst2*s2t*(3086*MuSUSY + 261*
        Mst2*Tbeta) + Mt*(1353*MuSUSY + 1730*Mst2*Tbeta))*pow2(Mt) + Mt*(-
        43276*MuSUSY + 25791*Mst2*Tbeta)*pow2(Mst2)*pow2(s2t) + 8306*Tbeta*
        pow3(s2t)*pow4(Mst2))))*pow8(Msq) + 3*pow2(Mgl)*(20*Dmsqst2*s2t*pow2(
        Mst2)*(pow2(Mst1)*(9*pow2(Mst2)*(-2*Mt*MuSUSY*s2t - 6*Tbeta*pow2(Mt) +
        Tbeta*pow2(Mst2)*pow2(s2t)) + pow2(Mst1)*(442*Mt*MuSUSY*s2t + 106*
        Tbeta*pow2(Mt) + 37*Tbeta*pow2(Mst2)*pow2(s2t)))*pow3(Dmsqst2) + pow6(
        Msq)*(18*pow2(Mst1)*pow2(Mst2)*(13*Mt*MuSUSY*s2t + 5*Tbeta*pow2(Mt) -
        5*Tbeta*pow2(Mst2)*pow2(s2t)) + 2*(221*Mt*MuSUSY*s2t + 53*Tbeta*pow2(
        Mt) - 52*Tbeta*pow2(Mst2)*pow2(s2t))*pow4(Mst1) + 54*Mt*(MuSUSY*s2t +
        Mt*Tbeta)*pow4(Mst2) - 27*Tbeta*pow2(s2t)*pow6(Mst2)) + Dmsqst2*
        xDmsqst2*pow4(Msq)*(3*pow2(Mst1)*pow2(Mst2)*(50*Mt*MuSUSY*s2t + 14*
        Tbeta*pow2(Mt) - 19*Tbeta*pow2(Mst2)*pow2(s2t)) + (442*Mt*MuSUSY*s2t +
        106*Tbeta*pow2(Mt) - 57*Tbeta*pow2(Mst2)*pow2(s2t))*pow4(Mst1) + 36*Mt*
        (MuSUSY*s2t + Mt*Tbeta)*pow4(Mst2) - 18*Tbeta*pow2(s2t)*pow6(Mst2))) -
        (27*pow4(Mst2)*(2*Mst2*s2t*(170*MuSUSY - 101*Mst2*Tbeta)*pow2(Mt) - 5*
        Mt*(16*MuSUSY + 51*Mst2*Tbeta)*pow2(Mst2)*pow2(s2t) + 4*(36*MuSUSY +
        103*Mst2*Tbeta)*pow3(Mt) + 40*Tbeta*pow3(s2t)*pow4(Mst2)) + 3*pow2(
        Mst1)*pow2(Mst2)*(10*Mst2*s2t*(908*MuSUSY - 347*Mst2*Tbeta)*pow2(Mt) -
        Mt*(4552*MuSUSY + 4515*Mst2*Tbeta)*pow2(Mst2)*pow2(s2t) + 16*(225*
        MuSUSY + 326*Mst2*Tbeta)*pow3(Mt) + 1916*Tbeta*pow3(s2t)*pow4(Mst2)) +
        pow4(Mst1)*(4*Mst2*s2t*(11764*MuSUSY - 6031*Mst2*Tbeta)*pow2(Mt) - 2*
        Mt*(13723*MuSUSY + 7431*Mst2*Tbeta)*pow2(Mst2)*pow2(s2t) + 8*(3093*
        MuSUSY + 3326*Mst2*Tbeta)*pow3(Mt) + 6895*Tbeta*pow3(s2t)*pow4(Mst2)))*
        pow8(Msq))) - z3*(3*pow2(Mgl)*(5*xDmsqst2*pow2(Dmsqst2)*pow2(Mst2)*
        pow4(Msq)*(9*MuSUSY*pow3(Mst2)*(4*s2t*(4032*MuSUSY + 1091*Mst2*Tbeta)*
        pow2(Mt) - 2*Mst2*Mt*(5507*MuSUSY + 6048*Mst2*Tbeta)*pow2(s2t) - 13824*
        Tbeta*pow3(Mt) + 5507*Tbeta*pow3(Mst2)*pow3(s2t)) + 3*Mst2*MuSUSY*pow2(
        Mst1)*(64*s2t*(1512*MuSUSY + 247*Mst2*Tbeta)*pow2(Mt) - 32*Mst2*Mt*(
        1783*MuSUSY + 1620*Mst2*Tbeta)*pow2(s2t) - 55296*Tbeta*pow3(Mt) +
        12007*Tbeta*pow3(Mst2)*pow3(s2t)) + 16*s2t*(3860*MuSUSY*Tbeta*pow2(Mt)
        + 6*pow2(Mst2)*pow2(s2t)*(263*MuSUSY*Tbeta - 135*Mst2*(-1 + pow2(Sbeta)
        )*pow2(Sbeta)) + Mt*s2t*(-7776*Mst2*MuSUSY*Tbeta + 3542*pow2(MuSUSY) +
        243*pow2(Mst2)*(-1 + pow2(Sbeta))*pow2(Sbeta)))*pow4(Mst1)) + 5*pow2(
        Mst2)*pow4(Dmsqst2)*(27*MuSUSY*pow3(Mst2)*((6912*MuSUSY*s2t + 652*Mst2*
        s2t*Tbeta)*pow2(Mt) - 2*Mst2*Mt*(4771*MuSUSY + 2592*Mst2*Tbeta)*pow2(
        s2t) - 7680*Tbeta*pow3(Mt) + 4771*Tbeta*pow3(Mst2)*pow3(s2t)) + 9*Mst2*
        MuSUSY*pow2(Mst1)*(576*s2t*(56*MuSUSY + 5*Mst2*Tbeta)*pow2(Mt) - 32*
        Mst2*Mt*(1853*MuSUSY + 756*Mst2*Tbeta)*pow2(s2t) - 18432*Tbeta*pow3(Mt)
        + 15335*Tbeta*pow3(Mst2)*pow3(s2t)) + 16*s2t*(3860*MuSUSY*Tbeta*pow2(
        Mt) + 2*pow2(Mst2)*pow2(s2t)*(3823*MuSUSY*Tbeta - 1215*Mst2*(-1 + pow2(
        Sbeta))*pow2(Sbeta)) + Mt*s2t*(-7776*Mst2*MuSUSY*Tbeta + 3542*pow2(
        MuSUSY) + 729*pow2(Mst2)*(-1 + pow2(Sbeta))*pow2(Sbeta)))*pow4(Mst1)) -
        80*Dmsqst2*MuSUSY*pow2(Mst2)*(-27*pow3(Mst2)*(4*s2t*(72*MuSUSY + 29*
        Mst2*Tbeta)*pow2(Mt) - 2*Mst2*Mt*(23*MuSUSY + 108*Mst2*Tbeta)*pow2(s2t)
        - 192*Tbeta*pow3(Mt) + 23*Tbeta*pow3(Mst2)*pow3(s2t)) + 18*Mst2*pow2(
        Mst1)*(-2*s2t*(504*MuSUSY + 101*Mst2*Tbeta)*pow2(Mt) + Mst2*Mt*(-35*
        MuSUSY + 432*Mst2*Tbeta)*pow2(s2t) + 576*Tbeta*pow3(Mt) + 52*Tbeta*
        pow3(Mst2)*pow3(s2t)) + 2*s2t*(-1771*Mt*MuSUSY*s2t + 3888*Mst2*Mt*s2t*
        Tbeta - 1930*Tbeta*pow2(Mt) + 728*Tbeta*pow2(Mst2)*pow2(s2t))*pow4(
        Mst1))*pow6(Msq) + 8*MuSUSY*(2*pow4(Mst1)*(2*Mst2*s2t*(99632*MuSUSY +
        13543*Mst2*Tbeta)*pow2(Mt) - 4*Mt*((-91963 - 162*lmMst1 + 162*lmMst2)*
        MuSUSY + 4431*Mst2*Tbeta)*pow2(Mst2)*pow2(s2t) + (9012*MuSUSY - 428288*
        Mst2*Tbeta)*pow3(Mt) - 76769*Tbeta*pow3(s2t)*pow4(Mst2)) - 3*pow2(Mst1)
        *pow2(Mst2)*(2*Mst2*s2t*(-58544*MuSUSY + 1859*Mst2*Tbeta)*pow2(Mt) +
        Mt*(-4*(35719 + 108*lmMst1 - 108*lmMst2)*MuSUSY + 34131*Mst2*Tbeta)*
        pow2(Mst2)*pow2(s2t) - 16*(981*MuSUSY - 8219*Mst2*Tbeta)*pow3(Mt) + 4*(
        11924 + 27*lmMst1 - 27*lmMst2)*Tbeta*pow3(s2t)*pow4(Mst2)) - 9*pow4(
        Mst2)*(-2*Mst2*s2t*(11930*MuSUSY + 9*(-77 + 8*lmMst1 - 8*lmMst2)*Mst2*
        Tbeta)*pow2(Mt) + 3*Mt*(-4*(1319 + 6*lmMst1 - 6*lmMst2)*MuSUSY + 5965*
        Mst2*Tbeta)*pow2(Mst2)*pow2(s2t) - 4*(920*MuSUSY - 2883*Mst2*Tbeta)*
        pow3(Mt) + 6*(1319 + 6*lmMst1 - 6*lmMst2)*Tbeta*pow3(s2t)*pow4(Mst2)))*
        pow8(Msq)) + 4*Dmglst2*Mst2*(2*Mgl*(9720*Dmsqst2*pow2(Mst2)*(MuSUSY*
        pow2(Mst1)*(-40*MuSUSY*s2t*pow2(Mt) + 8*Mst2*Mt*(2*MuSUSY + 3*Mst2*
        Tbeta)*pow2(s2t) + 160*Tbeta*pow3(Mt) - 5*Tbeta*pow3(Mst2)*pow3(s2t)) +
        MuSUSY*pow2(Mst2)*(-8*MuSUSY*s2t*pow2(Mt) + 6*Mst2*Mt*(MuSUSY + Mst2*
        Tbeta)*pow2(s2t) + 48*Tbeta*pow3(Mt) - 3*Tbeta*pow3(Mst2)*pow3(s2t)) +
        Mt*pow2(s2t)*(24*MuSUSY*Tbeta - 5*Mst2*(-1 + pow2(Sbeta))*pow2(Sbeta))*
        pow4(Mst1))*(pow3(Dmsqst2) + Dmsqst2*xDmsqst2*pow4(Msq) + pow6(Msq)) +
        MuSUSY*(4*(4*s2t*(395581*MuSUSY + 12627*Mst2*Tbeta)*pow2(Mt) + Mst2*Mt*
        (706232*MuSUSY - 546087*Mst2*Tbeta)*pow2(s2t) + 1668964*Tbeta*pow3(Mt)
        - 89947*Tbeta*pow3(Mst2)*pow3(s2t))*pow4(Mst1) - 27*pow3(Mst2)*(-4*
        Mst2*s2t*(15137*MuSUSY + 878*Mst2*Tbeta)*pow2(Mt) + Mt*(-35230*MuSUSY +
        45411*Mst2*Tbeta)*pow2(Mst2)*pow2(s2t) + 4*(664*MuSUSY - 9*Mst2*Tbeta)*
        pow3(Mt) + 17615*Tbeta*pow3(s2t)*pow4(Mst2)) - 9*Mst2*pow2(Mst1)*(-32*
        Mst2*s2t*(11864*MuSUSY + 885*Mst2*Tbeta)*pow2(Mt) + 3*Mt*(-77976*MuSUSY
        + 49501*Mst2*Tbeta)*pow2(Mst2)*pow2(s2t) + 2400*(5*MuSUSY - 67*Mst2*
        Tbeta)*pow3(Mt) + 64119*Tbeta*pow3(s2t)*pow4(Mst2)))*pow8(Msq)) +
        Dmglst2*xDmglst2*(19440*Dmsqst2*pow2(Mst2)*(MuSUSY*pow2(Mst1)*(-40*
        MuSUSY*s2t*pow2(Mt) + 8*Mst2*Mt*(2*MuSUSY + 3*Mst2*Tbeta)*pow2(s2t) +
        160*Tbeta*pow3(Mt) - 5*Tbeta*pow3(Mst2)*pow3(s2t)) + MuSUSY*pow2(Mst2)*
        (-8*MuSUSY*s2t*pow2(Mt) + 6*Mst2*Mt*(MuSUSY + Mst2*Tbeta)*pow2(s2t) +
        48*Tbeta*pow3(Mt) - 3*Tbeta*pow3(Mst2)*pow3(s2t)) + Mt*pow2(s2t)*(24*
        MuSUSY*Tbeta - 5*Mst2*(-1 + pow2(Sbeta))*pow2(Sbeta))*pow4(Mst1))*(
        pow3(Dmsqst2) + Dmsqst2*xDmsqst2*pow4(Msq) + pow6(Msq)) + MuSUSY*(8*(4*
        s2t*(395581*MuSUSY + 12627*Mst2*Tbeta)*pow2(Mt) + Mst2*Mt*(706232*
        MuSUSY - 546087*Mst2*Tbeta)*pow2(s2t) + 1668964*Tbeta*pow3(Mt) - 89947*
        Tbeta*pow3(Mst2)*pow3(s2t))*pow4(Mst1) - 9*pow3(Mst2)*(-2*Mst2*Mt*s2t*(
        Mst2*s2t*(228607*MuSUSY - 332871*Mst2*Tbeta) + Mt*(443828*MuSUSY +
        177327*Mst2*Tbeta)) + 64*(2985*MuSUSY + 4037*Mst2*Tbeta)*pow3(Mt) +
        228607*Tbeta*pow3(s2t)*pow4(Mst2)) - 3*Mst2*pow2(Mst1)*(-32*Mst2*s2t*(
        152143*MuSUSY + 28116*Mst2*Tbeta)*pow2(Mt) + 2*Mt*(-1618352*MuSUSY +
        827103*Mst2*Tbeta)*pow2(Mst2)*pow2(s2t) + 8*(163668*MuSUSY + 221681*
        Mst2*Tbeta)*pow3(Mt) + 932531*Tbeta*pow3(s2t)*pow4(Mst2)))*pow8(Msq))))
        ))/(pow2(Mgl)*pow8(Msq)) - 4*Mt*pow2(MuSUSY)*((90*pow2(Mt)*(12*Dmglst2*
        pow2(Mst2)*(pow2(Mst1)*(28405 - 288*B4 + 288*D3 - 144*DN + 5472*lmMst1
        - 288*(-89 + 2*lmMst1)*lmMst2 + 576*(-2 + lmMst1 - lmMst2)*lmMt + 224*
        OepS2 - 324*(65 + 14*lmMst1 - 14*lmMst2)*S2 - 576*(-9 + 8*lmMst1)*pow2(
        lmMst2) + 4608*pow3(lmMst2)) + 72*pow2(Mst2)*(180 - 2*B4 + 2*D3 - DN +
        144*lmMst2 - 216*S2 - 16*(-2 + lmMst1)*pow2(lmMst2) + 16*(lmMst1 +
        pow3(lmMst2)))) + Mgl*(18*pow2(Mst1)*pow2(Mst2)*(10667 - 96*B4 + 96*D3
        - 48*DN - 3072*lmMst1 - 384*(-25 + 9*lmMst1)*lmMst2 - 384*lmMst1*lmMt +
        384*(1 + lmMst2)*lmMt - 224*OepS2 + 324*(-43 + 14*lmMst1 - 14*lmMst2)*
        S2 - 384*(-13 + 4*lmMst1)*pow2(lmMst2) + 1536*pow3(lmMst2)) + (383185 -
        2592*B4 + 2592*D3 - 1296*DN - 187704*lmMst1 - 17440*OepS2 + 648*(-57 +
        545*lmMst1 - 545*lmMst2)*S2 - 7992*pow2(lmMst1) - 216*lmMst2*(-1733 +
        630*lmMst1 + 26*pow2(lmMst1)) - 216*(-859 + 246*lmMst1)*pow2(lmMst2) +
        3456*lmMt*(3 + 6*lmMst2 - 2*lmMst1*(3 + lmMst2) + pow2(lmMst1) + pow2(
        lmMst2)) + 720*pow3(lmMst1) + 58032*pow3(lmMst2))*pow4(Mst1) + 144*(436
        - 6*B4 + 6*D3 - 3*DN - 48*lmMst1 + (408 - 96*lmMst1)*lmMst2 + 24*lmMt -
        972*S2 - 48*(-4 + lmMst1)*pow2(lmMst2) + 48*pow3(lmMst2))*pow4(Mst2))))
        /Mgl - 120*T1ep*((2*Dmglst2*Mst2*(36*Mst2*pow2(Mst1)*(253*Mst2*Mt*s2t +
        42*pow2(Mt) - 129*pow2(Mst2)*pow2(s2t)) + 2*s2t*(16421*Mt - 8654*Mst2*
        s2t)*pow4(Mst1) + 189*s2t*(5*Mt - 2*Mst2*s2t)*pow4(Mst2)))/Mgl - 3*(6*
        pow2(Mst1)*pow2(Mst2)*(272*Mst2*Mt*s2t + 252*pow2(Mt) - 1057*pow2(Mst2)
        *pow2(s2t)) + (3764*Mst2*Mt*s2t + 6540*pow2(Mt) - 13237*pow2(Mst2)*
        pow2(s2t))*pow4(Mst1) + 54*s2t*(7*Mt - 23*Mst2*s2t)*pow5(Mst2)) + (60*
        pow2(Mst2)*pow2(s2t)*(pow4(Dmsqst2)*(-9*pow2(Mst1)*pow2(Mst2) + 221*
        pow4(Mst1)) + Dmsqst2*(117*pow2(Mst1)*pow2(Mst2) + 221*pow4(Mst1) + 27*
        pow4(Mst2))*pow6(Msq)))/pow8(Msq)) - (14580*pow2(Mst2)*pow2(s2t)*((1 -
        2*lmMst2)*shiftst3*(2*(1 - 2*lmMst1 + 2*lmMst2)*pow2(Mst2)*pow4(Mst1) +
        2*(1 - lmMst1 + lmMst2)*pow2(Mst1)*pow4(Mst2) + (2 - 6*lmMst1 + 6*
        lmMst2)*pow6(Mst1) + pow6(Mst2)) + (10*shiftst1*(pow4(Dmsqst2)*(4*
        lmMst1*(-3 + lmMst2)*pow2(Mst2)*pow4(Mst1) + 3*(pow2(Mst1) + pow2(Mst2)
        )*pow4(Mst2) + 2*lmMst1*(-3 + 2*lmMst2)*pow2(Mst1)*(pow4(Mst1) + pow4(
        Mst2)) - 4*pow2(lmMst2)*pow2(Mst1)*(pow2(Mst1)*pow2(Mst2) + pow4(Mst1)
        + pow4(Mst2)) + 2*lmMst2*(6*pow2(Mst2)*pow4(Mst1) + 2*pow2(Mst1)*pow4(
        Mst2) + 3*pow6(Mst1) - pow6(Mst2))) + ((pow2(Mst1) + pow2(Mst2))*pow4(
        Mst2) - 2*(lmMst1 - lmMst2)*pow2(Mst1)*(pow2(Mst1)*pow2(Mst2) + pow4(
        Mst1) + pow4(Mst2)))*(Dmsqst2*(3 - 2*lmMst2)*pow6(Msq) + (1 - 2*lmMst2)
        *pow8(Msq))))/pow8(Msq)))/pow2(Mst1) + Mst2*Mt*s2t*(30*(1702429 +
        257904*B4 - 648*DN - 748656*lmMst1 + 216*lmMst2*(5971 - 6106*lmMst1 +
        576*pow2(lmMst1)) + 41904*(pow2(lmMst1) + (34 - 15*lmMst1)*pow2(lmMst2)
        ) - 3456*pow3(lmMst1) - (1458*Dmglst2*(634.115454961134 - (3416*B4)/9.
         + (52*DN)/
        9. + (111356*lmMst1)/243. + (8*pow2(lmMst1))/9. - (8*lmMst2*(36802 -
        11421*lmMst1 + 1728*pow2(lmMst1)))/243. + (344*(-14 + 15*lmMst1)*pow2(
        lmMst2))/9. - (64*pow3(lmMst1))/9. - (1528*pow3(lmMst2))/3.))/Mgl +
        507600*pow3(lmMst2))*pow4(Mst1) + 43740*pow2(Mst1)*pow2(Mst2)*(
        621.9670781893004 + (1016*B4)/9. - (4*DN)/9. - (2476*lmMst1)/9. - (176*
        pow2(lmMst1))/9. + (4*lmMst2*(1427 - 1012*lmMst1 + 16*pow2(lmMst1)))/9.
         + (8*(642 - 203*lmMst1)*pow2(lmMst2))/9. + (520*pow3(lmMst2))/3. + (
        Dmglst2*(248*B4 - 4*DN + (66761 + 898740*lmMst2 + 2160*(11 + 4*lmMst2)*
        pow2(lmMst1) + 520560*pow2(lmMst2) - 180*lmMst1*(1141 + 1956*lmMst2 +
        1986*pow2(lmMst2)) + 348840*pow3(lmMst2))/1215.) + (160*Dmsqst2*(
        Dmglst2*(8 - 15*lmMst1 + 15*lmMst2) + (1 + 3*lmMst1 - 3*lmMst2)*Mgl)*(
        pow3(Dmsqst2) + pow6(Msq)))/(3.*pow8(Msq)))/Mgl) + (-160*OepS2*(3*Mgl*(
        816*pow2(Mst1)*pow2(Mst2) + 1882*pow4(Mst1) + 189*pow4(Mst2)) -
        Dmglst2*(9108*pow2(Mst1)*pow2(Mst2) + 32842*pow4(Mst1) + 945*pow4(Mst2)
        )) - 108*S2*(Dmglst2*(72*(2489 + 3795*lmMst1 - 3795*lmMst2)*pow2(Mst1)*
        pow2(Mst2) + 10*(123113 + 98526*lmMst1 - 98526*lmMst2)*pow4(Mst1) + 81*
        (-453 + 350*lmMst1 - 350*lmMst2)*pow4(Mst2)) - 15*Mgl*(36*(169 + 136*
        lmMst1 - 136*lmMst2)*pow2(Mst1)*pow2(Mst2) + 2*(9185 + 5646*lmMst1 -
        5646*lmMst2)*pow4(Mst1) + 81*(-1 + 14*lmMst1 - 14*lmMst2)*pow4(Mst2)))
        + (622080*(1 + lmMst2)*(-Dmglst2 + 3*Dmglst2*lmMst2 + Mgl + lmMst2*Mgl)
        *pow6(Mst2))/pow2(Mst1) - (27*pow4(Mst2)*(-45*Mgl*(960*(11 + 8*lmMst1 -
        8*lmMst2)*pow4(Dmsqst2) + 1920*Dmsqst2*(1 + lmMst1 - lmMst2)*pow6(Msq)
        + (9561 + 1760*B4 - 16*DN - 1984*lmMst1 - 256*pow2(lmMst1) - 64*lmMst2*
        (-214 + 73*lmMst1 + 4*pow2(lmMst1)) - 32*(-268 + 57*lmMst1)*pow2(
        lmMst2) + 2080*pow3(lmMst2))*pow8(Msq)) + Dmglst2*(86400*Dmsqst2*(2 +
        5*lmMst1 - 5*lmMst2)*(pow3(Dmsqst2) + pow6(Msq)) - (23917 + 188640*B4 -
        3600*DN - 37440*lmMst1 + 11520*pow2(lmMst1) - 2880*lmMst2*(-237 + 55*
        lmMst1 + 4*pow2(lmMst1)) - 1440*(-280 + 121*lmMst1)*pow2(lmMst2) +
        185760*pow3(lmMst2))*pow8(Msq))))/pow8(Msq))/Mgl))))/(174960.*pow6(
        Mst2)) + pow2(s2t)*((-5*Dmsqst2*Mt*pow2(Sbeta)*pow4(Mst1)*((4*Dmglst2*(
        1081 - 165*lmMst1 + 165*lmMst2)*Mt*pow2(Sbeta) - 9*Mgl*(2*(11 + 39*
        lmMst1 - 39*lmMst2)*Mst2*s2t*(-1 + pow2(Sbeta)) + (167 - 5*lmMst1 + 5*
        lmMst2)*Mt*pow2(Sbeta)))*pow3(Dmsqst2) + 4*Dmglst2*(1081 - 165*lmMst1 +
        165*lmMst2)*Mt*pow2(Sbeta)*pow6(Msq)))/(36.*Mgl*pow2(Mst2)*pow8(Msq)) +
        pow2(Mt)*((-5*pow2(Sbeta)*pow4(Mst1)*(9*(167 - 5*lmMst1 + 5*lmMst2)*
        Mgl*pow4(Dmsqst2) + 4*Dmglst2*Dmsqst2*(-1081 + 165*lmMst1 - 165*lmMst2)
        *(pow3(Dmsqst2) + pow6(Msq))))/(36.*Mgl*pow2(Mst2)*pow8(Msq)) + pow2(
        MuSUSY)*(53.385802469135804 + (40*B4)/9. - (4*D3)/9. + (2*DN)/9. + (
        1672*lmMst1)/27. + (53*pow2(lmMst1))/9. - lmMst2*(129.92592592592592 -
        72*lmMst1 + (13*pow2(lmMst1))/3.) + (2*(-470 + 147*lmMst1)*pow2(lmMst2)
        )/9. + (16*pow3(lmMst1))/9. - (271*pow3(lmMst2))/9. + S2*(921 + 138*
        lmMst1 - 138*lmMst2 + (Dmglst2*(3354 - 28*lmMst1 + 28*lmMst2))/Mgl + (
        30*Dmsqst2*(-15 + 2*lmMst1 - 2*lmMst2))/pow2(Msq) + (pow4(Mst1)*(583071
        + 119133*lmMst1 - 119133*lmMst2 + (8*Dmglst2*(93919 - 12981*lmMst1 +
        12981*lmMst2))/Mgl + (390*Dmsqst2*(95 + 102*lmMst1 - 102*lmMst2)*(pow3(
        Dmsqst2) + pow6(Msq)))/pow8(Msq)))/(81.*pow4(Mst2)) + (pow2(Mst1)*(
        11243 + 2114*lmMst1 - 2114*lmMst2 + (3*Dmglst2*(6952.444444444444 -
        344*lmMst1 + 344*lmMst2))/Mgl - (10*((913 + 6*lmMst1 - 6*lmMst2)*pow4(
        Dmsqst2) + Dmsqst2*(17 - 78*lmMst1 + 78*lmMst2)*pow6(Msq)))/pow8(Msq)))
        /(3.*pow2(Mst2)) - (1740*pow4(Dmsqst2))/pow8(Msq)) + (pow2(Mst1)*(
        434.270658436214 + (76*B4)/9. - (2*DN)/9. + (69088*lmMst1)/405. - (
        1313*pow2(lmMst1))/27. - (4*lmMst2*(16192 - 26430*lmMst1 + 3465*pow2(
        lmMst1)))/405. + ((-5735 + 3072*lmMst1)*pow2(lmMst2))/27. + (Dmsqst2*(
        201.74098765432097 + (622*lmMst1)/15. + (2*(-311 + 10*lmMst1)*lmMst2)/
        15. + (Dmglst2*(76.66666666666667 - 140*lmMst1 + 140*lmMst2))/Mgl - (
        22*pow2(lmMst1))/3. + 6*pow2(lmMst2)))/pow2(Msq) - (62*pow3(lmMst1))/
        27. - (2086*pow3(lmMst2))/27. - (2*Dmglst2*(2695042 - 40500*B4 + 54000*
        D3 - 33750*DN - 326895*lmMst1 - 324900*pow2(lmMst1) + 15*lmMst2*(-19607
        - 129030*lmMst1 + 62550*pow2(lmMst1)) - 450*(-5023 + 5610*lmMst1)*pow2(
        lmMst2) + 11250*pow3(lmMst1) + 1575000*pow3(lmMst2)))/(30375.*Mgl) + ((
        121.37234567901234 + (4003*lmMst1)/45. - ((4003 + 1020*lmMst1)*lmMst2)/
        45. + (Dmglst2*(76.66666666666667 - 140*lmMst1 + 140*lmMst2))/Mgl + (
        14*pow2(lmMst1))/3. + 18*pow2(lmMst2))*pow4(Dmsqst2))/pow8(Msq)))/pow2(
        Mst2) + (pow4(Mst1)*(628.1736268201578 + (76*B4)/9. - (2*DN)/9. + (
        6317839*lmMst1)/396900. - (66307*pow2(lmMst1))/315. - lmMst2*(
        12.52907281431091 - (182909*lmMst1)/315. + (274*pow2(lmMst1))/3.) + (2*
        (-58301 + 37135*lmMst1)*pow2(lmMst2))/315. - (44*pow3(lmMst1))/9. - (
        1256*pow3(lmMst2))/9. - (Dmglst2*(585.1892843532082 - (8*B4)/3. + (32*
        D3)/9. - (20*DN)/9. - (20109937*lmMst1)/297675. - (15886*pow2(lmMst1))/
        945. + lmMst2*(17.112243218274966 - (144628*lmMst1)/945. + (1180*pow2(
        lmMst1))/9.) + (2*(80257 - 138180*lmMst1)*pow2(lmMst2))/945. - (92*
        pow3(lmMst1))/27. + (4448*pow3(lmMst2))/27.))/Mgl + (Dmsqst2*(
        237.28785508324435 + (16526*lmMst2)/3969. + (2*lmMst1*(-8263 + 71820*
        lmMst2))/3969. - (520*pow2(lmMst1))/21. - (80*pow2(lmMst2))/7.)*(pow3(
        Dmsqst2) + pow6(Msq)))/pow8(Msq)))/pow4(Mst2) + ((-5*Dmsqst2*(288*
        Dmglst2*(4 + lmMst1 - lmMst2) + Mgl*(-1141 + 264*lmMst2 + 48*lmMst1*(-7
        + 3*lmMst2) - 144*pow2(lmMst2))))/(54.*pow2(Msq)) + (Dmglst2*(-14267 +
        432*B4 - 576*D3 + 360*DN - 4752*lmMst1 + 4536*lmMst2 + 16704*lmMst1*
        lmMst2 + 1404*pow2(lmMst1) - 1152*lmMst2*pow2(lmMst1) - 20232*pow2(
        lmMst2) + 8856*lmMst1*pow2(lmMst2) - 7704*pow3(lmMst2)))/162. + (16*(1
        + lmMst2)*(4*Dmglst2*lmMst2 + Mgl + lmMst2*Mgl)*pow4(Mst2))/(9.*pow4(
        Mst1)) - ((5*(576*Dmglst2*(4 + lmMst1 - lmMst2) + Mgl*(-4207 + 600*
        lmMst2 + 24*lmMst1*(-31 + 12*lmMst2) - 288*pow2(lmMst2)))*pow4(Dmsqst2)
        )/108. + (4*OepS2*(60*Mgl*pow2(Mst1)*(221*pow2(Mst1) - 9*pow2(Mst2))*
        pow4(Dmsqst2) + 60*Dmsqst2*Mgl*(117*pow2(Mst1)*pow2(Mst2) + 221*pow4(
        Mst1) + 27*pow4(Mst2))*pow6(Msq) + (-4*Dmglst2*(2322*pow2(Mst1)*pow2(
        Mst2) + 8654*pow4(Mst1) + 189*pow4(Mst2)) + 3*Mgl*(6342*pow2(Mst1)*
        pow2(Mst2) + 13237*pow4(Mst1) + 1242*pow4(Mst2)))*pow8(Msq)))/(2187.*
        pow4(Mst2)) + (pow2(Mst2)*(Mgl*(-30*(5 + 2*lmMst2)*pow4(Dmsqst2) + 30*
        Dmsqst2*(1 - 2*lmMst2)*pow6(Msq) + (103 + 186*lmMst2 + 32*lmMst1*(1 +
        lmMst2) + 91*pow2(lmMst2))*pow8(Msq)) + Dmglst2*(360*pow4(Dmsqst2) +
        360*Dmsqst2*pow6(Msq) + 2*(-20 + 198*lmMst2 + 32*lmMst1*lmMst2 + 123*
        pow2(lmMst2))*pow8(Msq))))/(9.*pow2(Mst1)))/pow8(Msq))/Mgl))))/Tbeta +
        MuSUSY*(pow2(Mt)*pow2(s2t)*(((32*(1 + lmMst2)*(-Dmglst2 + 3*Dmglst2*
        lmMst2 + Mgl + lmMst2*Mgl)*pow3(Mst2))/(3.*pow2(Mst1)) + (Dmglst2*(9*(
        36280*OepS2 - 27*(23989 + 27210*lmMst1 - 27210*lmMst2)*S2)*pow2(Mst1)*
        pow2(Mst2) + 2*(474680*OepS2 - 27*(525961 + 356010*lmMst1 - 356010*
        lmMst2)*S2)*pow4(Mst1) + 27*(1400*OepS2 - 81*(-453 + 350*lmMst1 - 350*
        lmMst2)*S2)*pow4(Mst2)) - 15*Mgl*(3*(1672*OepS2 - 81*(685 + 418*lmMst1
        - 418*lmMst2)*S2)*pow2(Mst1)*pow2(Mst2) + 2*(4264*OepS2 - 27*(6143 +
        3198*lmMst1 - 3198*lmMst2)*S2)*pow4(Mst1) + 27*(56*OepS2 + 81*(1 - 14*
        lmMst1 + 14*lmMst2)*S2)*pow4(Mst2)))/(14580.*pow3(Mst2)))/Mgl + Mst2*(
        188.52083333333334 + (110*B4)/3. - DN/3. - (124*lmMst1)/3. - (16*pow2(
        lmMst1))/3. - (4*lmMst2*(-198 + 73*lmMst1 + 4*pow2(lmMst1)))/3. + (168
        - 38*lmMst1)*pow2(lmMst2) + (130*pow3(lmMst2))/3. + ((-40*Dmsqst2*(
        Dmglst2*(2 + 5*lmMst1 - 5*lmMst2) + (-1 - lmMst1 + lmMst2)*Mgl))/pow2(
        Msq) + Dmglst2*(21.73935185185185 + (262*B4)/3. - (5*DN)/3. - (52*
        lmMst1)/3. + (16*pow2(lmMst1))/3. - (4*lmMst2*(-221 + 55*lmMst1 + 4*
        pow2(lmMst1)))/3. + (2*(232 - 121*lmMst1)*pow2(lmMst2))/3. + 86*pow3(
        lmMst2)) - (20*(2*Dmglst2*(2 + 5*lmMst1 - 5*lmMst2) + (-11 - 8*lmMst1 +
        8*lmMst2)*Mgl)*pow4(Dmsqst2))/pow8(Msq))/Mgl) + (pow4(Mst1)*(
        409.2597736625514 + 48*B4 - (1609*lmMst1)/9. + (326*pow2(lmMst1))/9. +
        (2*lmMst2*(845 - 1535*lmMst1 + 264*pow2(lmMst1)))/9. + (
        304.8888888888889 - 188*lmMst1)*pow2(lmMst2) - (16*pow3(lmMst1))/9. + (
        1180*pow3(lmMst2))/9. + (-(Dmglst2*(516.797085048011 - (296*B4)/3. + (
        4*DN)/3. + (17570*lmMst1)/81. + (46*pow2(lmMst1))/3. - (lmMst2*(28667 -
        5238*lmMst1 + 3024*pow2(lmMst1)))/81. + (4*(-60 + 157*lmMst1)*pow2(
        lmMst2))/3. - (16*pow3(lmMst1))/3. - (500*pow3(lmMst2))/3.)) + (10*
        Dmsqst2*(Dmglst2*(80 - 43*lmMst1 + 43*lmMst2) + (-4 + 9*lmMst1 - 9*
        lmMst2)*Mgl)*(pow3(Dmsqst2) + pow6(Msq)))/pow8(Msq))/Mgl))/pow3(Mst2) +
        (pow2(Mst1)*(267.2878086419753 + 48*B4 - 165*lmMst1 - (28*pow2(lmMst1))
        /3. + (lmMst2*(571 - 720*lmMst1 + 32*pow2(lmMst1)))/3. + (4*(187 - 73*
        lmMst1)*pow2(lmMst2))/3. + (80*Dmsqst2*(Dmglst2*(5 - 5*lmMst1 + 5*
        lmMst2) + (lmMst1 - lmMst2)*Mgl))/(Mgl*pow2(Msq)) + (260*pow3(lmMst2))/
        3. + (Dmglst2*(30.137808641975308 + (296*B4)/3. - (4*DN)/3. - (985*
        lmMst1)/9. + (28*pow2(lmMst1))/3. + (lmMst2*(2149 - 1296*lmMst1 + 96*
        pow2(lmMst1)))/9. + (134.66666666666666 - 140*lmMst1)*pow2(lmMst2) + (
        388*pow3(lmMst2))/3.))/Mgl + (20*(3 + 16*lmMst1 - 16*lmMst2 + (20*
        Dmglst2*(1 - lmMst1 + lmMst2))/Mgl)*pow4(Dmsqst2))/pow8(Msq)))/Mst2) -
        ((2*(-1 + 2*lmMst2)*s2t*shiftst3*pow2(Mst2)*pow3(Mt))/3. + (5*Mt*
        shiftst1*pow3(s2t)*pow4(Mst2)*(1 - 2*lmMst2 + (Dmsqst2*(3 - 2*lmMst2)*(
        pow3(Dmsqst2) + pow6(Msq)))/pow8(Msq)))/3.)/pow2(Mst1) + (Mt*pow2(z2)*(
        Mst2*xDmglst2*pow2(Dmglst2)*(2*(4*s2t*(86405*MuSUSY - 16542*Mst2*Tbeta)
        *pow2(Mt) - Mst2*Mt*(42392*MuSUSY + 105585*Mst2*Tbeta)*pow2(s2t) +
        129704*Tbeta*pow3(Mt) + 12664*Tbeta*pow3(Mst2)*pow3(s2t))*pow4(Mst1) +
        9*((41444*MuSUSY*s2t - 504*Mst2*s2t*Tbeta)*pow2(Mt) - Mst2*Mt*(2732*
        MuSUSY + 31083*Mst2*Tbeta)*pow2(s2t) - 28*Tbeta*pow3(Mt) + 1366*Tbeta*
        pow3(Mst2)*pow3(s2t))*pow4(Mst2) + 3*Mst2*pow2(Mst1)*(-8*Mst2*s2t*(-
        29549*MuSUSY + 1980*Mst2*Tbeta)*pow2(Mt) - 13*Mt*(712*MuSUSY + 6465*
        Mst2*Tbeta)*pow2(Mst2)*pow2(s2t) + 16*(105*MuSUSY - 1082*Mst2*Tbeta)*
        pow3(Mt) + 530*Tbeta*pow3(s2t)*pow4(Mst2)))*pow8(Msq) + Dmglst2*Mgl*
        Mst2*(2*(4*s2t*(86405*MuSUSY - 16542*Mst2*Tbeta)*pow2(Mt) - Mst2*Mt*(
        42392*MuSUSY + 105585*Mst2*Tbeta)*pow2(s2t) + 129704*Tbeta*pow3(Mt) +
        12664*Tbeta*pow3(Mst2)*pow3(s2t))*pow4(Mst1) + 27*(7052*MuSUSY*s2t*
        pow2(Mt) - Mst2*Mt*(632*MuSUSY + 5289*Mst2*Tbeta)*pow2(s2t) + 140*
        Tbeta*pow3(Mt) + 316*Tbeta*pow3(Mst2)*pow3(s2t))*pow4(Mst2) + 9*Mst2*
        pow2(Mst1)*(-80*Mst2*s2t*(-569*MuSUSY + 30*Mst2*Tbeta)*pow2(Mt) - 3*Mt*
        (1264*MuSUSY + 6091*Mst2*Tbeta)*pow2(Mst2)*pow2(s2t) + 96*(7*MuSUSY +
        62*Mst2*Tbeta)*pow3(Mt) + 948*Tbeta*pow3(s2t)*pow4(Mst2)))*pow8(Msq) -
        3*pow2(Mgl)*(-20*s2t*pow2(Mst1)*pow2(Mst2)*(9*pow2(Mst2)*(-2*Mt*MuSUSY*
        s2t - 6*Tbeta*pow2(Mt) + Tbeta*pow2(Mst2)*pow2(s2t)) + pow2(Mst1)*(442*
        Mt*MuSUSY*s2t + 106*Tbeta*pow2(Mt) + 37*Tbeta*pow2(Mst2)*pow2(s2t)))*
        pow4(Dmsqst2) + 20*s2t*xDmsqst2*pow2(Dmsqst2)*pow2(Mst2)*pow4(Msq)*(3*
        pow2(Mst1)*pow2(Mst2)*(-50*Mt*MuSUSY*s2t - 14*Tbeta*pow2(Mt) + 19*
        Tbeta*pow2(Mst2)*pow2(s2t)) + (-442*Mt*MuSUSY*s2t - 106*Tbeta*pow2(Mt)
        + 57*Tbeta*pow2(Mst2)*pow2(s2t))*pow4(Mst1) + 18*(-2*Mt*MuSUSY*s2t - 2*
        Tbeta*pow2(Mt) + Tbeta*pow2(Mst2)*pow2(s2t))*pow4(Mst2)) + 20*Dmsqst2*
        s2t*pow2(Mst2)*(18*pow2(Mst1)*pow2(Mst2)*(-13*Mt*MuSUSY*s2t - 5*Tbeta*
        pow2(Mt) + 5*Tbeta*pow2(Mst2)*pow2(s2t)) - 2*(221*Mt*MuSUSY*s2t + 53*
        Tbeta*pow2(Mt) - 52*Tbeta*pow2(Mst2)*pow2(s2t))*pow4(Mst1) + 27*(-2*Mt*
        MuSUSY*s2t - 2*Tbeta*pow2(Mt) + Tbeta*pow2(Mst2)*pow2(s2t))*pow4(Mst2))
        *pow6(Msq) + (3*pow2(Mst1)*pow2(Mst2)*(-2*Mst2*s2t*(6368*MuSUSY + 1735*
        Mst2*Tbeta)*pow2(Mt) + Mt*(-3364*MuSUSY + 4557*Mst2*Tbeta)*pow2(Mst2)*
        pow2(s2t) + 16*(63*MuSUSY + 110*Mst2*Tbeta)*pow3(Mt) + 1700*Tbeta*pow3(
        s2t)*pow4(Mst2)) + pow4(Mst1)*(-4*Mst2*s2t*(13670*MuSUSY + 6031*Mst2*
        Tbeta)*pow2(Mt) + 2*Mt*(-11941*MuSUSY + 6177*Mst2*Tbeta)*pow2(Mst2)*
        pow2(s2t) + 40*(327*MuSUSY + 406*Mst2*Tbeta)*pow3(Mt) + 6895*Tbeta*
        pow3(s2t)*pow4(Mst2)) - 27*(2*s2t*(370*MuSUSY + 53*Mst2*Tbeta)*pow2(Mt)
        - Mst2*Mt*(4*MuSUSY + 555*Mst2*Tbeta)*pow2(s2t) - 28*Tbeta*pow3(Mt) +
        2*Tbeta*pow3(Mst2)*pow3(s2t))*pow5(Mst2))*pow8(Msq))))/(972.*Tbeta*
        pow2(Mgl)*pow6(Mst2)*pow8(Msq)) + (Mt*(60*T1ep*(-3*Mgl*(3*pow2(Mst1)*
        pow2(Mst2)*(-3470*Mst2*s2t*pow2(Mt) - 627*Mt*pow2(Mst2)*pow2(s2t) +
        1760*pow3(Mt) + 1700*pow3(Mst2)*pow3(s2t)) + (-24124*Mst2*s2t*pow2(Mt)
        - 3198*Mt*pow2(Mst2)*pow2(s2t) + 16240*pow3(Mt) + 6895*pow3(Mst2)*pow3(
        s2t))*pow4(Mst1) + 27*(-106*Mst2*s2t*pow2(Mt) - 21*Mt*pow2(Mst2)*pow2(
        s2t) + 28*pow3(Mt) + 46*pow3(Mst2)*pow3(s2t))*pow4(Mst2)) - (60*
        Dmsqst2*Mgl*Mst2*s2t*(-2*(53*pow2(Mt) - 52*pow2(Mst2)*pow2(s2t))*pow4(
        Mst1) + 27*(-2*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow4(Mst2) + 90*pow2(
        Mst1)*(-(pow2(Mst2)*pow2(Mt)) + pow2(s2t)*pow4(Mst2))))/pow2(Msq) +
        Dmglst2*(27*pow2(Mst1)*pow2(Mst2)*(-800*Mst2*s2t*pow2(Mt) - 907*Mt*
        pow2(Mst2)*pow2(s2t) + 1984*pow3(Mt) + 316*pow3(Mst2)*pow3(s2t)) + 2*(-
        66168*Mst2*s2t*pow2(Mt) - 35601*Mt*pow2(Mst2)*pow2(s2t) + 129704*pow3(
        Mt) + 12664*pow3(Mst2)*pow3(s2t))*pow4(Mst1) + 189*(20*pow3(Mt)*pow4(
        Mst2) - 15*Mt*pow2(s2t)*pow6(Mst2) + 4*pow3(s2t)*pow7(Mst2))) + (60*
        Mgl*Mst2*s2t*pow2(Mst1)*(9*pow2(Mst2)*(-6*pow2(Mt) + pow2(Mst2)*pow2(
        s2t)) + pow2(Mst1)*(106*pow2(Mt) + 37*pow2(Mst2)*pow2(s2t)))*pow4(
        Dmsqst2))/pow8(Msq)) + pow3(Mt)*(20*(Dmglst2*(5659217 + 1592460*lmMt -
        972*(569 + 180*lmMst2 - 126*lmMt)*pow2(lmMst1) + 324*(5353 + 126*lmMt)*
        pow2(lmMst2) + 36*lmMst2*(39031 + 3204*lmMt - 5184*pow2(lmMt)) - 72*
        lmMst1*(27653 + 3015*lmMt + 18*lmMst2*(689 + 126*lmMt) - 2160*pow2(
        lmMst2) - 2592*pow2(lmMt)) - 186624*pow2(lmMt) + 1944*pow3(lmMst1) +
        17496*pow3(lmMst2)) - 3*Mgl*(678923 + 19440*lmMt + 540*(661 - 30*lmMt)*
        pow2(lmMst2) - 15552*pow2(lmMt) - 108*((457 - 12*lmMst2 - 126*lmMt)*
        pow2(lmMst1) + lmMst1*(1841 + lmMst2*(2626 - 24*lmMt) - 420*lmMt + 840*
        pow2(lmMst2) - 288*pow2(lmMt)) + lmMst2*(619 + 498*lmMt + 288*pow2(
        lmMt))) + 216*pow3(lmMst1) + 89208*pow3(lmMst2)))*pow4(Mst1) + 160*
        OepS2*(3*Mgl*(1320*pow2(Mst1)*pow2(Mst2) + 4060*pow4(Mst1) + 189*pow4(
        Mst2)) - Dmglst2*(13392*pow2(Mst1)*pow2(Mst2) + 64852*pow4(Mst1) + 945*
        pow4(Mst2))) - 108*S2*(3*Mgl*(72*(457 + 550*lmMst1 - 550*lmMst2)*pow2(
        Mst1)*pow2(Mst2) + 40*(4082 + 3045*lmMst1 - 3045*lmMst2)*pow4(Mst1) +
        81*(-79 + 70*lmMst1 - 70*lmMst2)*pow4(Mst2)) - Dmglst2*(144*(601 +
        2790*lmMst1 - 2790*lmMst2)*pow2(Mst1)*pow2(Mst2) + 40*(46187 + 48639*
        lmMst1 - 48639*lmMst2)*pow4(Mst1) + 2025*(-27 + 14*lmMst1 - 14*lmMst2)*
        pow4(Mst2))) - (622080*(1 + lmMst2)*(-Dmglst2 + 3*Dmglst2*lmMst2 + Mgl
        + lmMst2*Mgl)*pow6(Mst2))/pow2(Mst1) + (-72*pow2(Mst1)*pow2(Mst2)*(
        1800*Dmglst2*Dmsqst2*(497 - 306*lmMst1 + 609*lmMst2 - 303*lmMt)*(pow3(
        Dmsqst2) + pow6(Msq)) + 16200*Dmsqst2*(1 + 6*lmMst1 - 13*lmMst2 + 7*
        lmMt)*Mgl*(pow3(Dmsqst2) + pow6(Msq)) - 3*Dmglst2*(97837 + 71580*lmMt -
        360*(23 + 8*lmMst2 - 4*lmMt)*pow2(lmMst1) + 720*(97 + 2*lmMt)*pow2(
        lmMst2) - 8640*pow2(lmMt) - 30*lmMst2*(-2911 + 396*lmMt + 144*pow2(
        lmMt)) + 10*lmMst1*(-6157 + 642*lmMt - 6*lmMst2*(791 + 48*lmMt) + 882*
        pow2(lmMst2) + 432*pow2(lmMt)) - 5940*pow3(lmMst2))*pow8(Msq) + 2*Mgl*(
        110219 - 1080*lmMt + 540*(-15 + 4*lmMt)*pow2(lmMst1) + 810*(121 - 8*
        lmMt)*pow2(lmMst2) - 135*lmMst1*(337 + lmMst2*(588 - 32*lmMt) - 132*
        lmMt + 194*pow2(lmMst2) - 48*pow2(lmMt)) - 6480*pow2(lmMt) - 405*
        lmMst2*(31 + 54*lmMt + 16*pow2(lmMt)) + 26190*pow3(lmMst2))*pow8(Msq))
        - 27*pow4(Mst2)*(5*Dmglst2*(960*Dmsqst2*(101 - 90*lmMst1 + 177*lmMst2 -
        87*lmMt)*(pow3(Dmsqst2) + pow6(Msq)) - (57495 + lmMst2*(55616 - 4704*
        lmMt) + 45408*lmMt + 2304*(-1 + lmMst2)*pow2(lmMst1) + 96*(79 + 48*
        lmMt)*pow2(lmMst2) + 288*lmMst1*(-22 + 16*lmMt - 2*lmMst2*(9 + 8*lmMt)
        + 41*pow2(lmMst2)) - 14112*pow3(lmMst2))*pow8(Msq)) + 9*Mgl*(1200*(8 +
        32*lmMst1 - 65*lmMst2 + 33*lmMt)*pow4(Dmsqst2) + 4800*Dmsqst2*(5 + 2*
        lmMst1 - 5*lmMst2 + 3*lmMt)*pow6(Msq) + (18813 + 480*lmMst2*(3 - 11*
        lmMt) - 8000*lmMt - 1280*(1 + lmMst2)*pow2(lmMst1) + 160*(163 - 16*
        lmMt)*pow2(lmMst2) - 160*lmMst1*(46 + 2*lmMst2*(57 - 8*lmMt) - 16*lmMt
        + 41*pow2(lmMst2)) - 3840*pow2(lmMt) + 7840*pow3(lmMst2))*pow8(Msq))))/
        pow8(Msq))))/(43740.*Mgl*pow5(Mst2)) + s2t*pow3(Mt)*(302.29104938271604
         + (8*D3)/9. - (8*DN)/9. + (614*lmMst1)/27. + (32*lmMt)/3. + (124*pow2(
        lmMst1))/9. - (2*lmMst2*(-1901 + 6*lmMst1 + 117*pow2(lmMst1)))/27. + (
        46 + (32*lmMst1)/9.)*pow2(lmMst2) + (32*pow3(lmMst1))/9. + (14*pow3(
        lmMst2))/9. - S2*(42.3 - 159*lmMst1 + 159*lmMst2 + (2532*Dmglst2)/(5.*
        Mgl) - (30*Dmsqst2*(9 + 2*lmMst1 - 2*lmMst2))/pow2(Msq) - (pow2(Mst1)*(
        570.4333333333333 + (1735*lmMst1)/3. - (1735*lmMst2)/3. - (4*Dmglst2*(
        139 + 100*lmMst1 - 100*lmMst2))/Mgl + (10*Dmsqst2*(127 + 30*lmMst1 -
        30*lmMst2))/(3.*pow2(Msq)) + ((156.66666666666666 - 60*lmMst1 + 60*
        lmMst2)*pow4(Dmsqst2))/pow8(Msq)))/pow2(Mst2) - (2*pow4(Mst1)*(138286 +
        90465*lmMst1 - 90465*lmMst2 - (12*Dmglst2*(17269 + 13785*lmMst1 -
        13785*lmMst2))/Mgl + (25*Dmsqst2*(1283 + 318*lmMst1 - 318*lmMst2)*(
        pow3(Dmsqst2) + pow6(Msq)))/pow8(Msq)))/(135.*pow4(Mst2)) - (240*pow4(
        Dmsqst2))/pow8(Msq)) + ((-75*Dmsqst2*(656*Dmglst2 + (-367 - 32*lmMst1 +
        80*lmMst2)*Mgl))/pow2(Msq) + 8*Dmglst2*(12454 - 120*B4 + 120*D3 - 60*DN
        + 690*lmMst1 + 345*pow2(lmMst1) - 5*lmMst2*(-3121 + 42*lmMst1 + 96*
        pow2(lmMst1)) + (3630 - 480*lmMst1)*pow2(lmMst2) + 960*pow3(lmMst2)) -
        (25*(1968*Dmglst2 + (-1453 + 48*lmMst1 + 96*lmMst2)*Mgl)*pow4(Dmsqst2))
        /pow8(Msq))/(270.*Mgl) + (pow4(Mst1)*(599.6859861506036 - (99165049*
        lmMst1)/198450. + lmMst2*(713.9201259763164 - (57994*lmMst1)/945. - (
        350*pow2(lmMst1))/9.) - (90913*pow2(lmMst1))/945. + (200.24021164021164
         - (286*lmMst1)/9.)*pow2(lmMst2) + (32*lmMt*(1 + 4*lmMst2 - 2*lmMst1*(2
        + lmMst2) + pow2(lmMst1) + pow2(lmMst2)))/3. + (26*pow3(lmMst1))/27. +
        (1882*pow3(lmMst2))/27. + (Dmglst2*(311.41128771790903 - (32*B4)/9. + (
        32*D3)/9. - (16*DN)/9. + (14926634*lmMst1)/19845. + (19352*pow2(lmMst1)
        )/189. + (2*lmMst2*(2917823 - 3813600*lmMst1 + 97020*pow2(lmMst1)))/
        19845. + (8*(8677 - 5439*lmMst1)*pow2(lmMst2))/189. - (64*lmMt*(4 + 6*
        lmMst2 - 2*lmMst1*(3 + lmMst2) + pow2(lmMst1) + pow2(lmMst2)))/3. - (8*
        pow3(lmMst1))/9. + (664*pow3(lmMst2))/3.))/Mgl + (Dmsqst2*(
        87.53310385647498 - (941812*lmMst1)/19845. + (47.45840262030738 + (484*
        lmMst1)/21.)*lmMst2 - (242*(pow2(lmMst1) + pow2(lmMst2)))/21.)*(pow3(
        Dmsqst2) + pow6(Msq)))/pow8(Msq)))/pow4(Mst2) + ((32*(1 + lmMst2)*(4*
        Dmglst2*lmMst2 + Mgl + lmMst2*Mgl)*pow4(Mst2))/(9.*pow4(Mst1)) - (4*
        OepS2*(20*Mgl*pow2(Mst1)*(53*pow2(Mst1) - 27*pow2(Mst2))*pow4(Dmsqst2)
        + 20*Dmsqst2*Mgl*(45*pow2(Mst1)*pow2(Mst2) + 53*pow4(Mst1) + 27*pow4(
        Mst2))*pow6(Msq) + (-24*Dmglst2*(150*pow2(Mst1)*pow2(Mst2) + 919*pow4(
        Mst1)) + Mgl*(5205*pow2(Mst1)*pow2(Mst2) + 12062*pow4(Mst1) + 1431*
        pow4(Mst2)))*pow8(Msq)))/(729.*pow4(Mst2)*pow8(Msq)) + ((-2*pow2(Mst2)*
        (2*Dmglst2*(180*pow4(Dmsqst2) + 180*Dmsqst2*pow6(Msq) + (-52 + 102*
        lmMst2 + 32*lmMst1*lmMst2 + 59*pow2(lmMst2))*pow8(Msq)) + Mgl*(-30*(5 +
        2*lmMst2)*pow4(Dmsqst2) + 30*Dmsqst2*(1 - 2*lmMst2)*pow6(Msq) + (71 +
        122*lmMst2 + 32*lmMst1*(1 + lmMst2) + 59*pow2(lmMst2))*pow8(Msq))))/(9.
        *pow2(Mst1)) - (pow2(Mst1)*(32*Dmglst2*(50625*Dmsqst2*(65 - 2*lmMst1 +
        2*lmMst2)*(pow3(Dmsqst2) + pow6(Msq)) - 2*(1928479 - 13500*B4 + 13500*
        D3 - 6750*DN + 635385*lmMst1 + 81000*(-2 + lmMst1 - lmMst2)*lmMt + 450*
        pow2(lmMst1) - 15*lmMst2*(-153166 + 17835*lmMst1 + 3375*pow2(lmMst1)) -
        225*(-2627 + 1695*lmMst1)*pow2(lmMst2) - 1125*pow3(lmMst1) + 433125*
        pow3(lmMst2))*pow8(Msq)) - 5*Mgl*(108*(79499 + 20*lmMst1*(53 - 270*
        lmMst2) - 1060*lmMst2 + 2700*pow2(lmMst1) + 2700*pow2(lmMst2))*pow4(
        Dmsqst2) + 12*Dmsqst2*(339977 + 96120*lmMst2 + 1080*lmMst1*(-89 + 60*
        lmMst2) - 32400*pow2(lmMst1) - 32400*pow2(lmMst2))*pow6(Msq) + (
        19425643 + 518400*lmMt + 240*lmMst2*(82997 + 4320*lmMt) - 3600*(683 +
        96*lmMst2)*pow2(lmMst1) + 5684400*pow2(lmMst2) - 240*lmMst1*(42047 +
        4800*lmMst2 + 4320*lmMt + 6390*pow2(lmMst2)) - 295200*pow3(lmMst1) +
        2174400*pow3(lmMst2))*pow8(Msq))))/(243000.*pow2(Mst2)))/pow8(Msq))/
        Mgl) + Mt*((5*s2t*shiftst1*(4*pow2(Mst2)*pow2(Mt) - 2*pow2(Mst1)*(2*
        pow2(Mt) + (-lmMst1 + lmMst2)*pow2(Mst2)*pow2(s2t)) + pow2(s2t)*pow4(
        Mst1))*(Dmsqst2*(3 - 2*lmMst2)*(pow3(Dmsqst2) + pow6(Msq)) + (1 - 2*
        lmMst2)*pow8(Msq)))/(3.*pow2(Mst1)*pow8(Msq)) + pow3(s2t)*(((-1 + 2*
        lmMst2)*shiftst3*pow2(Mst2)*((-1 - lmMst1 + lmMst2)*pow2(Mst1) + pow2(
        Mst2)))/(6.*pow2(Mst1)) - (pow4(Mst1)*(96.95148419197191 - (61388401*
        lmMst1)/793800. + lmMst2*(73.69595742000504 + (302047*lmMst1)/1890. - (
        257*pow2(lmMst1))/9.) - (76483*pow2(lmMst1))/945. - (78.87883597883598
         - 61*lmMst1)*pow2(lmMst2) + (Dmsqst2*(
        17.77343371446168 - (452768*lmMst1)/19845. + (22.815217939027463 + (
        122*lmMst1)/7.)*lmMst2 + (5*Dmglst2*(419 - 57*lmMst1 + 57*lmMst2))/(27.
        *Mgl) - (61*(pow2(lmMst1) + pow2(lmMst2)))/7.))/pow2(Msq) - (35*pow3(
        lmMst1))/27. - (841*pow3(lmMst2))/27. - (Dmglst2*(203.8689796251638 + (
        2171669*lmMst2)/119070. + ((433 + 6552*lmMst2)*pow2(lmMst1))/189. + (
        1987*pow2(lmMst2))/189. - (lmMst1*(2740559 + 1524600*lmMst2 + 7514640*
        pow2(lmMst2)))/119070. - (56*pow3(lmMst1))/27. + (824*pow3(lmMst2))/27.
        ))/Mgl - ((55.45597659639988 + (27119*lmMst1)/11340. + (-
        2.3914462081128747 + 16*lmMst1)*lmMst2 + (5*Dmglst2*(-419 + 57*lmMst1 -
        57*lmMst2))/(27.*Mgl) - 8*(pow2(lmMst1) + pow2(lmMst2)))*pow4(Dmsqst2))
        /pow8(Msq)))/pow2(Mst2) + (pow2(Mst1)*(15*Mgl*(594509 + 495540*lmMst2 +
        180*lmMst1*(-2453 + 420*lmMst2) - 37800*pow2(lmMst1) - 37800*pow2(
        lmMst2))*pow4(Dmsqst2) - 120*Dmsqst2*Mgl*(97294 - 17235*lmMst2 + 45*
        lmMst1*(233 + 330*lmMst2) - 7425*pow2(lmMst1) - 7425*pow2(lmMst2))*
        pow6(Msq) + 405000*Dmglst2*Dmsqst2*(-55 + 34*lmMst1 - 34*lmMst2)*(pow3(
        Dmsqst2) + pow6(Msq)) - 10*Mgl*(4627751 + 48600*B4 + 5400*D3 - 5400*DN
        + 1320240*lmMst1 - 662400*pow2(lmMst1) - 30*lmMst2*(12148 - 76560*
        lmMst1 + 12105*pow2(lmMst1)) + 2250*(-583 + 438*lmMst1)*pow2(lmMst2) -
        49500*pow3(lmMst1) - 572850*pow3(lmMst2))*pow8(Msq) + 2*Dmglst2*(
        5430043 + 524580*lmMst2 + 900*(-859 + 3690*lmMst2)*pow2(lmMst1) +
        1454400*pow2(lmMst2) - 60*lmMst1*(51493 + 24630*lmMst2 + 112950*pow2(
        lmMst2)) + 45000*pow3(lmMst1) + 3411000*pow3(lmMst2))*pow8(Msq)))/(
        243000.*Mgl*pow8(Msq)) + ((-8*(1 + lmMst2)*(4*Dmglst2*lmMst2 + Mgl +
        lmMst2*Mgl)*pow6(Mst2))/(9.*pow4(Mst1)) - (30*Mgl*pow4(Dmsqst2)*(9*(8*
        OepS2 - 27*(391 + 6*lmMst1 - 6*lmMst2)*S2)*pow2(Mst1)*pow2(Mst2) + (
        296*OepS2 - 27*(3871 + 222*lmMst1 - 222*lmMst2)*S2)*pow4(Mst1) -
        126846*S2*pow4(Mst2)) - 30*Dmsqst2*Mgl*(18*(40*OepS2 - 27*(59 + 30*
        lmMst1 - 30*lmMst2)*S2)*pow2(Mst1)*pow2(Mst2) + 4*(208*OepS2 - 27*(347
        + 156*lmMst1 - 156*lmMst2)*S2)*pow4(Mst1) + 27*(8*OepS2 - 81*(-15 + 2*
        lmMst1 - 2*lmMst2)*S2)*pow4(Mst2))*pow6(Msq) + (3*Mgl*(-60*(340*OepS2 -
        81*(424 + 85*lmMst1 - 85*lmMst2)*S2)*pow2(Mst1)*pow2(Mst2) - 35*(788*
        OepS2 - 27*(2662 + 591*lmMst1 - 591*lmMst2)*S2)*pow4(Mst1) + 27*(-184*
        OepS2 + 81*(307 + 46*lmMst1 - 46*lmMst2)*S2)*pow4(Mst2)) + 2*Dmglst2*(
        27*(632*OepS2 + 9*(16193 - 1422*lmMst1 + 1422*lmMst2)*S2)*pow2(Mst1)*
        pow2(Mst2) + 2*(25328*OepS2 + 27*(47051 - 18996*lmMst1 + 18996*lmMst2)*
        S2)*pow4(Mst1) + 27*(56*OepS2 - 81*(-1677 + 14*lmMst1 - 14*lmMst2)*S2)*
        pow4(Mst2)))*pow8(Msq))/(4374.*pow2(Mst2)*pow8(Msq)) + (pow2(Mst2)*(-(
        Mgl*(15*(3847 + lmMst1*(744 - 288*lmMst2) - 744*lmMst2 + 288*pow2(
        lmMst2))*pow4(Dmsqst2) + 30*Dmsqst2*(1177 + 48*lmMst1*(7 - 3*lmMst2) -
        336*lmMst2 + 144*pow2(lmMst2))*pow6(Msq) + (21005 + 1440*B4 - 144*D3 +
        72*DN + 21216*lmMst1 + 1908*pow2(lmMst1) - 12*lmMst2*(2950 - 2040*
        lmMst1 + 117*pow2(lmMst1)) + 108*(-283 + 98*lmMst1)*pow2(lmMst2) + 576*
        pow3(lmMst1) - 9756*pow3(lmMst2))*pow8(Msq))) + 2*Dmglst2*(2160*
        Dmsqst2*(5 + 2*lmMst1 - 2*lmMst2)*(pow3(Dmsqst2) + pow6(Msq)) + (14987
        - 432*B4 + 576*D3 - 360*DN + 4752*lmMst1 - 1404*pow2(lmMst1) + 144*
        lmMst2*(-81 - 124*lmMst1 + 8*pow2(lmMst1)) - 36*(-439 + 246*lmMst1)*
        pow2(lmMst2) + 7704*pow3(lmMst2))*pow8(Msq)) + (36*pow2(Mst2)*(Mgl*(-
        30*(5 + 2*lmMst2)*pow4(Dmsqst2) + 30*Dmsqst2*(1 - 2*lmMst2)*pow6(Msq) +
        (119 + 218*lmMst2 + 32*lmMst1*(1 + lmMst2) + 107*pow2(lmMst2))*pow8(
        Msq)) + Dmglst2*(360*pow4(Dmsqst2) + 360*Dmsqst2*pow6(Msq) + 2*(-20 +
        230*lmMst2 + 32*lmMst1*lmMst2 + 155*pow2(lmMst2))*pow8(Msq))))/pow2(
        Mst1)))/(648.*pow8(Msq)))/Mgl)) - xDR2DRMOD*((-15*shiftst1*(4*s2t*(
        pow2(Mst1) - pow2(Mst2))*pow3(Mt) - Mt*pow3(s2t)*(2*(lmMst1 - lmMst2)*
        pow2(Mst1)*pow2(Mst2) + pow4(Mst1) - pow4(Mst2)))*(1 + Dmsqst2/pow2(
        Msq) + pow4(Dmsqst2)/pow8(Msq)))/(2.*pow2(Mst1)) + (MuSUSY*s2t*pow2(Mt)
        *((64*Mt*((1 + lmMst2)*Mgl*(2*(1 - 10*lmMst2 + 4*lmMst1*(2 + lmMst2) -
        4*pow2(lmMst2))*pow2(Mst2)*pow4(Mst1) - 2*(1 + 5*lmMst2 - lmMst1*(3 +
        2*lmMst2) + 2*pow2(lmMst2))*pow2(Mst1)*pow4(Mst2) + (7 - 32*lmMst2 + 4*
        lmMst1*(7 + 3*lmMst2) - 12*pow2(lmMst2))*pow6(Mst1) - 2*(1 + lmMst2)*
        pow6(Mst2)) + Dmglst2*(4*pow2(Mst2)*(8 + 13*lmMst2 - 8*pow2(lmMst2) +
        lmMst1*(-5 + 5*lmMst2 + 6*pow2(lmMst2)) - 6*pow3(lmMst2))*pow4(Mst1) +
        (49 + 103*lmMst2 - 36*(1 + lmMst2)*pow2(lmMst2) + 4*lmMst1*(-11 + 6*
        lmMst2 + 9*pow2(lmMst2)))*pow6(Mst1) + 2*(pow2(Mst1)*(8 + 7*lmMst2 -
        11*pow2(lmMst2) + lmMst1*(-3 + 5*lmMst2 + 6*pow2(lmMst2)) - 6*pow3(
        lmMst2))*pow4(Mst2) + (1 - 2*lmMst2 - 3*pow2(lmMst2))*pow6(Mst2)))))/(
        9.*Mgl*pow2(Mst1)*pow5(Mst2)) - (s2t*(756 - 2024*lmMst1 + 128*pow2(
        lmMst1) + 8*lmMst2*(363 - 332*lmMst1 + 16*pow2(lmMst1)) + 4*(707 - 246*
        lmMst1)*pow2(lmMst2) + 856*pow3(lmMst2) - 540*shiftst1*(1 + pow2(Mst2)/
        pow2(Mst1) - (2*(lmMst1 - lmMst2)*(pow2(Mst1)*pow2(Mst2) + pow4(Mst1) +
        pow4(Mst2)))/pow4(Mst2))*(1 + Dmsqst2/pow2(Msq) + pow4(Dmsqst2)/pow8(
        Msq)) + (8*Dmglst2*(60 + 206*lmMst2 + 32*lmMst2*pow2(lmMst1) + lmMst1*(
        8 - 460*lmMst2 - 246*pow2(lmMst2)) + 519*pow2(lmMst2) + 214*pow3(
        lmMst2)) - (64*(1 + lmMst2)*(4*Dmglst2*lmMst2 + Mgl + lmMst2*Mgl)*pow4(
        Mst2))/pow4(Mst1) + (30*Dmsqst2*(1 - 2*lmMst1 + 2*lmMst2)*((38*Dmglst2
        - 5*Mgl)*pow3(Dmsqst2) - 2*(2*Dmglst2 - 5*Mgl)*pow6(Msq)))/pow8(Msq) +
        (4*pow2(Mst2)*pow4(Mst1)*(60*Dmglst2*Dmsqst2*(lmMst1 - lmMst2)*(pow3(
        Dmsqst2) + pow6(Msq)) - 15*(lmMst1 - lmMst2)*Mgl*(13*pow4(Dmsqst2) +
        10*Dmsqst2*pow6(Msq)) + 4*Dmglst2*(48 + 4*lmMst2*(31 + 36*pow2(lmMst1))
        + 278*pow2(lmMst2) - lmMst1*(44 + 278*lmMst2 + 379*pow2(lmMst2)) + 235*
        pow3(lmMst2))*pow8(Msq) + 2*Mgl*(32 + 285*lmMst2 + 144*(1 + lmMst2)*
        pow2(lmMst1) + 444*pow2(lmMst2) - lmMst1*(253 + 588*lmMst2 + 379*pow2(
        lmMst2)) + 235*pow3(lmMst2))*pow8(Msq)) + 8*pow6(Mst1)*(-75*Dmsqst2*(
        lmMst1 - lmMst2)*Mgl*(pow3(Dmsqst2) + pow6(Msq)) + (2*Dmglst2*(48 + 4*
        lmMst2*(45 + 68*pow2(lmMst1)) + 310*pow2(lmMst2) - lmMst1*(92 + 310*
        lmMst2 + 635*pow2(lmMst2)) + 363*pow3(lmMst2)) + Mgl*(40 + 277*lmMst2 +
        272*(1 + lmMst2)*pow2(lmMst1) + 556*pow2(lmMst2) - lmMst1*(237 + 828*
        lmMst2 + 635*pow2(lmMst2)) + 363*pow3(lmMst2)))*pow8(Msq)) + pow6(Mst2)
        *(Mgl*(-105*pow4(Dmsqst2) + 300*Dmsqst2*pow6(Msq) + 4*(205 + 252*lmMst2
        + 32*lmMst1*(1 + lmMst2) + 91*pow2(lmMst2))*pow8(Msq)) + 4*Dmglst2*(15*
        pow4(Dmsqst2) - 30*Dmsqst2*pow6(Msq) + 2*(-20 + 198*lmMst2 + 32*lmMst1*
        lmMst2 + 123*pow2(lmMst2))*pow8(Msq))))/(pow2(Mst1)*pow4(Mst2)*pow8(
        Msq)))/Mgl))/36.))/Tbeta + (Mt*(512*pow2(Mst1)*pow3(Mt)*(9*(1 + lmMst2)
        *Mgl*((3 - 21*lmMst2 + 2*lmMst1*(8 + 3*lmMst2 - 3*lmMt) + 5*lmMt + 6*
        lmMst2*lmMt - 6*pow2(lmMst2))*pow2(Mst2)*pow4(Mst1) + (-1 - 7*lmMst2 +
        2*lmMst1*(2 + lmMst2 - lmMt) + 3*lmMt + 2*lmMst2*lmMt - 2*pow2(lmMst2))
        *pow2(Mst1)*pow4(Mst2) + (12 - 45*lmMst2 + 2*lmMst1*(19 + 6*lmMst2 - 6*
        lmMt) + 7*lmMt + 12*lmMst2*lmMt - 12*pow2(lmMst2))*pow6(Mst1) - 2*(1 +
        lmMst2)*pow6(Mst2)) + Dmglst2*(pow2(Mst2)*(263 + lmMst2*(962 - 213*
        lmMt) - 177*lmMt + (393 - 18*lmMt)*pow2(lmMst2) - 18*lmMst1*(26 -
        lmMst2*(-17 + lmMt) - 7*lmMt + pow2(lmMst2)) + 18*pow3(lmMst2))*pow4(
        Mst1) - pow2(Mst1)*(17*(-7 + 3*lmMt) + lmMst2*(-224 + 15*lmMt) - 3*(5 +
        6*lmMt)*pow2(lmMst2) - 18*lmMst1*(-4 + lmMt - lmMst2*(1 + lmMt) + pow2(
        lmMst2)) + 18*pow3(lmMst2))*pow4(Mst2) + (290 + lmMst2*(2447 - 645*
        lmMt) - 375*lmMt - 3*(-509 + 60*lmMt)*pow2(lmMst2) - 18*lmMst1*(87 +
        lmMst2*(71 - 10*lmMt) - 22*lmMt + 10*pow2(lmMst2)) + 180*pow3(lmMst2))*
        pow6(Mst1) - 18*(-1 + 2*lmMst2 + 3*pow2(lmMst2))*pow6(Mst2))) - 3456*
        Mt*pow2(Mst1)*pow2(Mst2)*pow2(s2t)*((1 + lmMst2)*Mgl*(2*(2 - 5*lmMst2 +
        lmMst1*(5 + 2*lmMst2) - 2*pow2(lmMst2))*pow2(Mst2)*pow4(Mst1) + 2*(-2*
        lmMst2*(2 + lmMst2) + lmMst1*(3 + 2*lmMst2))*pow2(Mst1)*pow4(Mst2) + (5
        - 12*lmMst2 + 4*lmMst1*(3 + lmMst2) - 4*pow2(lmMst2))*pow6(Mst1) - 2*(1
        + lmMst2)*pow6(Mst2)) + Dmglst2*(2*pow2(Mst2)*(8 + 19*lmMst2 - 5*pow2(
        lmMst2) + lmMst1*(-7 + 5*lmMst2 + 6*pow2(lmMst2)) - 6*pow3(lmMst2))*
        pow4(Mst1) + 2*pow2(Mst1)*(7 + 9*lmMst2 - 8*pow2(lmMst2) + lmMst1*(-3 +
        5*lmMst2 + 6*pow2(lmMst2)) - 6*pow3(lmMst2))*pow4(Mst2) + (17 + 51*
        lmMst2 - 4*pow2(lmMst2) + 4*lmMst1*(-6 + lmMst2 + 3*pow2(lmMst2)) - 12*
        pow3(lmMst2))*pow6(Mst1) - 2*(-1 + 2*lmMst2 + 3*pow2(lmMst2))*pow6(
        Mst2))) - (12*Mst2*s2t*pow2(Mt)*(3*Mgl*pow2(Mst1)*pow4(Mst2)*(15*(10*
        pow2(Mst1) - 7*pow2(Mst2))*pow4(Dmsqst2) - 300*Dmsqst2*(pow2(Mst1) -
        pow2(Mst2))*pow6(Msq)) + 4*Dmglst2*pow2(Mst1)*pow4(Mst2)*(45*(-19*pow2(
        Mst1) + pow2(Mst2))*pow4(Dmsqst2) + 90*Dmsqst2*(pow2(Mst1) - pow2(Mst2)
        )*pow6(Msq)) - 8*Dmglst2*pow8(Msq)*((52 + 370*lmMst2 - 96*lmMst2*pow2(
        lmMst1) + 81*pow2(lmMst2) + 96*lmMst1*(-1 + lmMst2 + 2*pow2(lmMst2)) -
        96*pow3(lmMst2))*pow4(Mst1)*pow4(Mst2) - 96*(lmMst2*pow2(lmMst1) +
        lmMst1*(4 + 2*lmMst2 - 2*pow2(lmMst2)) + (-4 + lmMst2)*pow2(1 + lmMst2)
        )*pow2(Mst2)*pow6(Mst1) - 3*(-52 + 2*(51 + 16*lmMst1)*lmMst2 + 59*pow2(
        lmMst2))*pow2(Mst1)*pow6(Mst2) - 96*(-6 - 14*lmMst2 + lmMst2*pow2(
        lmMst1) + lmMst1*(9 + 6*lmMst2 - 2*pow2(lmMst2)) - 6*pow2(lmMst2) +
        pow3(lmMst2))*pow8(Mst1) + 96*lmMst2*(1 + lmMst2)*pow8(Mst2)) + 12*Mgl*
        pow8(Msq)*(-((109 + 76*lmMst2 - 21*pow2(lmMst2) + 64*lmMst1*pow2(1 +
        lmMst2) - 32*((1 + lmMst2)*pow2(lmMst1) + pow3(lmMst2)))*pow4(Mst1)*
        pow4(Mst2)) + (173 + 188*lmMst2 + 32*lmMst1*(1 + lmMst2) + 59*pow2(
        lmMst2))*pow2(Mst1)*pow6(Mst2) + 32*(1 + lmMst2)*(pow2(1 - lmMst1 +
        lmMst2)*pow2(Mst2)*pow6(Mst1) + (1 + 3*lmMst2 - lmMst1*(3 + 2*lmMst2) +
        pow2(lmMst1) + pow2(lmMst2))*pow8(Mst1)) - 16*pow2(1 + lmMst2)*pow8(
        Mst2))))/pow8(Msq) + (9*pow3(Mst2)*pow3(s2t)*(2*Dmglst2*(60*lmMst2*
        pow4(Mst1)*pow4(Mst2)*(pow4(Dmsqst2) - 2*Dmsqst2*pow6(Msq)) + pow2(
        Mst1)*pow2(Mst2)*(pow4(Dmsqst2)*(-570*pow4(Mst1) + 75*pow4(Mst2)) + 60*
        Dmsqst2*(pow4(Mst1) - pow4(Mst2))*pow6(Msq) + 16*(20*pow2(Mst1)*pow2(
        Mst2) + 9*pow4(Mst1) - 5*pow4(Mst2))*pow8(Msq)) - 4*lmMst1*(15*pow4(
        Mst1)*pow4(Mst2)*(pow4(Dmsqst2) - 2*Dmsqst2*pow6(Msq)) + pow8(Msq)*(2*(
        -4 + 246*lmMst2 + 123*pow2(lmMst2))*pow4(Mst1)*pow4(Mst2) + 32*(3 + 3*
        lmMst2 + 16*pow2(lmMst2))*pow2(Mst2)*pow6(Mst1) - 32*lmMst2*pow2(Mst1)*
        pow6(Mst2) + 32*(3 + 2*lmMst2 + 16*pow2(lmMst2))*pow8(Mst1))) + 4*
        lmMst2*pow8(Msq)*(32*pow2(lmMst1)*pow4(Mst1)*(8*pow2(Mst1)*pow2(Mst2) +
        8*pow4(Mst1) + pow4(Mst2)) + lmMst2*(2*(198 + 107*lmMst2)*pow4(Mst1)*
        pow4(Mst2) + (37 + 256*lmMst2)*pow2(Mst2)*pow6(Mst1) + 155*pow2(Mst1)*
        pow6(Mst2) + 64*(1 + 4*lmMst2)*pow8(Mst1) - 32*pow8(Mst2))) + 8*lmMst2*
        pow8(Msq)*(4*pow4(Mst1)*pow4(Mst2) + 21*pow2(Mst2)*pow6(Mst1) + 115*
        pow2(Mst1)*pow6(Mst2) + 56*pow8(Mst1) - 16*pow8(Mst2))) + Mgl*(pow2(
        Mst1)*pow2(Mst2)*(15*pow4(Dmsqst2)*(14*(lmMst1 - lmMst2)*pow2(Mst1)*
        pow2(Mst2) + 10*pow4(Mst1) - 7*pow4(Mst2)) - 300*Dmsqst2*(2*(lmMst1 -
        lmMst2)*pow2(Mst1)*pow2(Mst2) + pow4(Mst1) - pow4(Mst2))*pow6(Msq)) +
        4*pow8(Msq)*(2*(-8 + 237*lmMst2 + 16*(1 + lmMst2)*pow2(lmMst1) + 308*
        pow2(lmMst2) - lmMst1*(269 + 348*lmMst2 + 123*pow2(lmMst2)) + 107*pow3(
        lmMst2))*pow4(Mst1)*pow4(Mst2) + pow2(Mst2)*(-125 - 156*lmMst2 - 512*
        lmMst1*lmMst2*(1 + lmMst2) + 256*(1 + lmMst2)*pow2(lmMst1) + 181*pow2(
        lmMst2) + 256*pow3(lmMst2))*pow6(Mst1) + (221 + 284*lmMst2 + 32*lmMst1*
        (1 + lmMst2) + 107*pow2(lmMst2))*pow2(Mst1)*pow6(Mst2) + 16*(1 +
        lmMst2)*(1 + lmMst1*(2 - 32*lmMst2) - 2*lmMst2 + 16*pow2(lmMst1) + 16*
        pow2(lmMst2))*pow8(Mst1) - 16*pow2(1 + lmMst2)*pow8(Mst2)))))/pow8(Msq)
        ))/(648.*Mgl*pow4(Mst1)*pow5(Mst2))))) + (Al4p*(xDmglst2*pow2(Dmglst2)*
        (-(Mt*MuSUSY*twoLoopFlag*(-(Mt*pow2(s2t)*(Mst2*(28 + 4*lmMst1 - 2*
        lmMst1*lmMst2 + 2*pow2(lmMst2)) + ((34 + 8*lmMst1 - 4*(2 + lmMst1)*
        lmMst2 + 4*pow2(lmMst2))*pow2(Mst1))/Mst2 + ((37.5 + lmMst1*(5 - 4*
        lmMst2) - 5*lmMst2 + 4*pow2(lmMst2))*pow4(Mst1))/pow3(Mst2))) + (2*s2t*
        pow2(Mt)*(11 - 18*lmMst2 + (6*(3 - lmMst1 + lmMst2)*pow2(Mst1))/pow2(
        Mst2) + (6*(-2 + 3*lmMst2)*pow2(Mst2))/pow2(Mst1) + (6*(7 - 6*lmMst1 +
        6*lmMst2)*pow4(Mst1))/pow4(Mst2)))/9. + pow3(s2t)*(((2 + 3*lmMst2)*
        pow2(Mst1))/3. + (2.6666666666666665 + lmMst2 + lmMst1*(-1 + 2*lmMst2)
        - 2*pow2(lmMst2))*pow2(Mst2) + ((0.5 - lmMst1 + lmMst2)*pow4(Mst1))/
        pow2(Mst2) + ((0.6666666666666666 - lmMst2)*pow4(Mst2))/pow2(Mst1)) + (
        4*pow3(Mt)*((29 - 18*lmMst1*(-1 + lmMst2) - 3*lmMst2 - 15*lmMt + 18*
        pow2(lmMst2))*pow2(Mst1)*pow2(Mst2) + (236 + 339*lmMst2 - 18*lmMst1*(13
        + 4*lmMst2 - 3*lmMt) - 105*lmMt - 54*lmMst2*lmMt + 72*pow2(lmMst2))*
        pow4(Mst1) + (-79 - 18*lmMst1*(-2 + lmMst2) - 39*lmMst2 + 3*lmMt + 18*
        pow2(lmMst2))*pow4(Mst2)))/(27.*pow5(Mst2)))) + (s2t*twoLoopFlag*pow2(
        Mt)*pow2(MuSUSY)*(2*(4*Mt*(-31 + 3*lmMst1*(-2 + lmMst2) + 4*lmMst2 - 3*
        pow2(lmMst2)) + 3*Mst2*s2t*(4 + lmMst2 + lmMst1*(-1 + 2*lmMst2) - 2*
        pow2(lmMst2)))*pow2(Mst2)*pow4(Mst1) + 2*(Mst2*s2t*(10 + lmMst1*(-3 +
        6*lmMst2) - 6*pow2(lmMst2)) + 4*Mt*(-14 + lmMst1*(-2 + lmMst2) - pow2(
        lmMst2)))*pow2(Mst1)*pow4(Mst2) + (Mt*(-398 + 52*lmMst2 + lmMst1*(-68 +
        40*lmMst2) - 40*pow2(lmMst2)) + 3*Mst2*s2t*(9 + 4*lmMst1*(-1 + lmMst2)
        + 4*lmMst2 - 4*pow2(lmMst2)))*pow6(Mst1) + 2*(2 - 3*lmMst2)*s2t*pow7(
        Mst2)))/(3.*Tbeta*pow2(Mst1)*pow5(Mst2)) + Al4p*Mt*MuSUSY*
        threeLoopFlag*((T1ep*(3*pow2(Mst1)*pow2(Mst2)*(-15840*Mst2*s2t*pow2(Mt)
        + 9267*Mt*pow2(Mst2)*pow2(s2t) - 17312*pow3(Mt) + 530*pow3(Mst2)*pow3(
        s2t)) + 2*(-66168*Mst2*s2t*pow2(Mt) - 35601*Mt*pow2(Mst2)*pow2(s2t) +
        129704*pow3(Mt) + 12664*pow3(Mst2)*pow3(s2t))*pow4(Mst1) + 63*(-72*
        Mst2*s2t*pow2(Mt) + 3*Mt*pow2(Mst2)*pow2(s2t) - 4*pow3(Mt) + 10*pow3(
        Mst2)*pow3(s2t))*pow4(Mst2)))/(729.*pow5(Mst2)) + pow3(s2t)*((-16*(-2 +
        lmMst2 + 5*pow2(lmMst2))*pow6(Mst2))/(9.*pow4(Mst1)) + ((-4*OepS2*(795*
        pow2(Mst1)*pow2(Mst2) + 12664*pow4(Mst1) + 315*pow4(Mst2)))/2187. - (
        S2*(159*(21401 - 150*lmMst1 + 150*lmMst2)*pow2(Mst1)*pow2(Mst2) + 20*(
        47051 - 18996*lmMst1 + 18996*lmMst2)*pow4(Mst1) + 27*(116129 - 350*
        lmMst1 + 350*lmMst2)*pow4(Mst2)))/810. + pow4(Mst1)*(203.8689796251638
         - (2740559*lmMst1)/
        119070. + (433*pow2(lmMst1))/189. + lmMst2*(18.238590744939952 - (2420*
        lmMst1)/189. + (104*pow2(lmMst1))/3.) + (10.513227513227513 - (568*
        lmMst1)/9.)*pow2(lmMst2) - (56*pow3(lmMst1))/27. + (824*pow3(lmMst2))/
        27. + (5*Dmsqst2*(-419 + 57*lmMst1 - 57*lmMst2)*(pow3(Dmsqst2) + pow6(
        Msq)))/(27.*pow8(Msq))))/pow2(Mst2) - (pow4(Mst2)*(916 + 64*lmMst1 - 6*
        (65 + 16*lmMst1)*lmMst2 - 657*pow2(lmMst2) - (360*Dmsqst2*(pow3(
        Dmsqst2) + pow6(Msq)))/pow8(Msq)))/(18.*pow2(Mst1)) + pow2(Mst1)*(
        32.30218244170096 - (32*B4)/9. + (32*D3)/9. - (16*DN)/9. - (80549*
        lmMst1)/1350. - (6011*pow2(lmMst1))/270. + (lmMst2*(41949 + 45410*
        lmMst1 + 58350*pow2(lmMst1)))/1350. - ((2383 + 11595*lmMst1)*pow2(
        lmMst2))/135. - (5*pow3(lmMst1))/27. + (1157*pow3(lmMst2))/27. + (5*
        Dmsqst2*(-55 + 34*lmMst1 - 34*lmMst2)*(pow3(Dmsqst2) + pow6(Msq)))/(3.*
        pow8(Msq))) + pow2(Mst2)*(92.72211934156378 - (50*B4)/9. + (56*D3)/9. -
        (31*DN)/9. + (142*lmMst1)/3. - lmMst2*(91.29629629629629 + 56*lmMst1 -
        (16*pow2(lmMst1))/3.) - (103*pow2(lmMst1))/18. + (35.05555555555556 -
        41*lmMst1)*pow2(lmMst2) + (107*pow3(lmMst2))/3. + (20*Dmsqst2*(5 + 2*
        lmMst1 - 2*lmMst2)*(pow3(Dmsqst2) + pow6(Msq)))/(3.*pow8(Msq)))) - Mt*
        pow2(s2t)*((64*(2 + lmMst2 - 3*pow2(lmMst2))*pow3(Mst2))/(3.*pow2(Mst1)
        ) + ((-2*OepS2*(-9267*pow2(Mst1)*pow2(Mst2) + 23734*pow4(Mst1) - 63*
        pow4(Mst2)))/729. + (S2*(3*(29123 - 92670*lmMst1 + 92670*lmMst2)*pow2(
        Mst1)*pow2(Mst2) + 2*(525961 + 356010*lmMst1 - 356010*lmMst2)*pow4(
        Mst1) + 27*(1387 - 70*lmMst1 + 70*lmMst2)*pow4(Mst2)))/540. + pow4(
        Mst1)*(516.797085048011 - (296*B4)/3. + (4*DN)/3. + (17570*lmMst1)/81.
         + (46*pow2(lmMst1))/
        3. - (lmMst2*(28667 - 5238*lmMst1 + 3024*pow2(lmMst1)))/81. + (4*(-60 +
        157*lmMst1)*pow2(lmMst2))/3. - (16*pow3(lmMst1))/3. - (500*pow3(lmMst2)
        )/3. + (10*Dmsqst2*(-80 + 43*lmMst1 - 43*lmMst2)*(pow3(Dmsqst2) + pow6(
        Msq)))/pow8(Msq)))/pow3(Mst2) + Mst2*(358.7766975308642 - (490*B4)/3. +
        (11*DN)/3. - (824*lmMst1)/3. - (32*pow2(lmMst1))/3. + (2*lmMst2*(-695 -
        108*lmMst1 + 24*pow2(lmMst1)))/9. + (14*(-18 + 31*lmMst1)*pow2(lmMst2))
        /3. - 150*pow3(lmMst2) + (40*Dmsqst2*(2 + 5*lmMst1 - 5*lmMst2)*(pow3(
        Dmsqst2) + pow6(Msq)))/pow8(Msq)) + (pow2(Mst1)*(69.80036008230452 - (
        524*B4)/3. + (10*DN)/3. - (3469*lmMst1)/27. + (64*pow2(lmMst1))/3. - (
        lmMst2*(7619 + 720*lmMst1 + 288*pow2(lmMst1)))/27. + (4*(4 + 153*
        lmMst1)*pow2(lmMst2))/3. - (580*pow3(lmMst2))/3. + (400*Dmsqst2*(-1 +
        lmMst1 - lmMst2)*(pow3(Dmsqst2) + pow6(Msq)))/pow8(Msq)))/Mst2) + pow3(
        Mt)*((256*Mst2*(2 + lmMst2 - 3*pow2(lmMst2)))/(9.*pow2(Mst1)) + (24*(
        151480*OepS2 + 27*(105619 - 113610*lmMst1 + 113610*lmMst2)*S2)*pow2(
        Mst1)*pow2(Mst2) + 35*(5659217 + 1592460*lmMt - 518816*OepS2 + 9976392*
        S2 - 972*(569 + 180*lmMst2 - 126*lmMt)*pow2(lmMst1) + 324*(5353 + 126*
        lmMt)*pow2(lmMst2) - 186624*pow2(lmMt) + 72*lmMst1*(-27653 - 3015*lmMt
        - 18*lmMst2*(689 + 126*lmMt) + 145917*S2 + 2160*pow2(lmMst2) + 2592*
        pow2(lmMt)) - 36*lmMst2*(-39031 - 3204*lmMt + 291834*S2 + 5184*pow2(
        lmMt)) + 1944*pow3(lmMst1) + 17496*pow3(lmMst2))*pow4(Mst1) + 9*(1960*
        OepS2 - 81*(-9361 + 490*lmMst1 - 490*lmMst2)*S2)*pow4(Mst2))/(76545.*
        pow5(Mst2)) - (210.23283362727807 + (2272*lmMst1)/9. + (128*pow2(
        lmMst1))/9. - (16*lmMst2*(1099 - 126*lmMst1 + 12*pow2(lmMst1)))/27. + (
        8*lmMt*(458 + 48*lmMst1*(-2 + lmMst2) + 105*lmMst2 - 48*pow2(lmMst2)))/
        27. - (8*(107 + 41*lmMst1)*pow2(lmMst2))/9. + (32*pow2(lmMt))/9. + (
        392*pow3(lmMst2))/9. + (80*Dmsqst2*(101 - 90*lmMst1 + 177*lmMst2 - 87*
        lmMt)*(pow3(Dmsqst2) + pow6(Msq)))/(27.*pow8(Msq)))/Mst2 + (pow2(Mst1)*
        (223.31232608269644 - (13840*lmMst1)/27. + (4*lmMst2*(8232 + 347*lmMst2
        + lmMst1*(-689 + 438*lmMst2) - 48*pow2(lmMst1)))/27. + (64*pow2(lmMst1)
        )/3. - (4*lmMt*(5564 + 203*lmMst2 + lmMst1*(-281 + 96*lmMst2) - 96*
        pow2(lmMst2)))/27. + (160*pow2(lmMt))/9. - (520*pow3(lmMst2))/9. + (80*
        Dmsqst2*(-497 + 306*lmMst1 - 609*lmMst2 + 303*lmMt)*(pow3(Dmsqst2) +
        pow6(Msq)))/(27.*pow8(Msq))))/pow3(Mst2)) + s2t*pow2(Mt)*(
        1657.490476190476 - (16*B4)/3. + (16*D3)/3. - (8*DN)/3. + (736*lmMst1)/
        9. - (32*lmMt)/9. + lmMst2*(958.2962962962963 + 68*lmMst1 - (64*pow2(
        lmMst1))/3.) + (58*pow2(lmMst1))/3. + (16*(37 - 12*lmMst1)*pow2(lmMst2)
        )/9. + (128*pow3(lmMst2))/3. - S2*(596.0571428571428 + 84*lmMst1 - 84*
        lmMst2 + ((1403.352380952381 + 880*lmMst1 - 880*lmMst2)*pow2(Mst1))/
        pow2(Mst2) + (8*(17269 + 13785*lmMst1 - 13785*lmMst2)*pow4(Mst1))/(45.*
        pow4(Mst2))) + (64*(-2 + lmMst2 + 5*pow2(lmMst2))*pow4(Mst2))/(9.*pow4(
        Mst1)) + ((311.41128771790903 - (32*B4)/9. + (32*D3)/9. - (16*DN)/9. +
        (14926634*lmMst1)/19845. + (19352*pow2(lmMst1))/189. + (2*lmMst2*(
        2917823 - 3813600*lmMst1 + 97020*pow2(lmMst1)))/19845. + (8*(8677 -
        5439*lmMst1)*pow2(lmMst2))/189. - (64*lmMt*(4 + 6*lmMst2 - 2*lmMst1*(3
        + lmMst2) + pow2(lmMst1) + pow2(lmMst2)))/3. - (8*pow3(lmMst1))/9. + (
        664*pow3(lmMst2))/3.)*pow4(Mst1) + (16*OepS2*(660*pow2(Mst1)*pow2(Mst2)
        + 1838*pow4(Mst1) + 63*pow4(Mst2)))/243.)/pow4(Mst2) + (pow2(Mst1)*(
        2311.8293850676073 - (16*B4)/3. + (16*D3)/3. - (8*DN)/3. + (1014244*
        lmMst1)/2025. + (32*(5 + lmMst1 - lmMst2)*lmMt)/9. + (4*lmMst2*(715964
        + 124935*lmMst1 - 11025*pow2(lmMst1)))/2025. + (2032*pow2(lmMst1))/135.
         - (4*(4517 + 5025*lmMst1)*pow2(lmMst2))/135. + (4*pow3(lmMst1))/27. +
        (4604*pow3(lmMst2))/27. + (20*Dmsqst2*(-65 + 2*lmMst1 - 2*lmMst2)*(
        pow3(Dmsqst2) + pow6(Msq)))/(3.*pow8(Msq))))/pow2(Mst2) + (2*(-820*
        Dmsqst2*(pow3(Dmsqst2) + pow6(Msq)) - (pow2(Mst2)*(360*pow4(Dmsqst2) +
        360*Dmsqst2*pow6(Msq) + (-756 + 134*lmMst2 + 32*lmMst1*(-2 + 3*lmMst2)
        + 177*pow2(lmMst2))*pow8(Msq)))/pow2(Mst1)))/(9.*pow8(Msq)))) + (Al4p*
        threeLoopFlag*((5*Dmsqst2*(-1081 + 165*lmMst1 - 165*lmMst2)*pow2(Mt)*
        pow2(s2t)*(-1 + pow2(Sbeta))*pow2(Sbeta)*pow4(Mst1)*(pow3(Dmsqst2) +
        pow6(Msq)))/(9.*pow2(Mst2)*pow8(Msq)) - pow2(MuSUSY)*(((pow2(Mst1)*(
        3051661 - 12960*B4 + 12960*D3 - 6480*DN + 304320*lmMst1 + 960*(2227 +
        285*lmMst1)*lmMst2 + 2880*(4 + lmMst1 - lmMst2)*lmMt + 5600*OepS2 -
        324*(5361 + 350*lmMst1 - 350*lmMst2)*S2 - 2880*(23 + 72*lmMst1)*pow2(
        lmMst2) + 207360*pow3(lmMst2)) + 24*pow2(Mst2)*(50134 - 270*B4 + 270*D3
        - 135*DN + 120*(271 + 24*lmMst1)*lmMst2 - 120*lmMt - 29646*S2 - 720*(-2
        + 3*lmMst1)*pow2(lmMst2) + 2160*(lmMst1 + pow3(lmMst2))))*pow4(Mt))/(
        1215.*pow4(Mst2)) - s2t*pow3(Mt)*((256*Mst2*(2 + lmMst2 - 3*pow2(
        lmMst2)))/(9.*pow2(Mst1)) + ((634.115454961134 - (3416*B4)/9. + (52*DN)
        /9. + (111356*lmMst1)/243. + (8*pow2(lmMst1))/9. - (8*lmMst2*(36802 -
        11421*lmMst1 + 1728*pow2(lmMst1)))/243. + (344*(-14 + 15*lmMst1)*pow2(
        lmMst2))/9. - (64*pow3(lmMst1))/9. - (1528*pow3(lmMst2))/3.)*pow4(Mst1)
        - (8*OepS2*(-9330*pow2(Mst1)*pow2(Mst2) + 32842*pow4(Mst1) - 63*pow4(
        Mst2)))/2187. + (S2*(6*(20803 - 46650*lmMst1 + 46650*lmMst2)*pow2(Mst1)
        *pow2(Mst2) + 10*(123113 + 98526*lmMst1 - 98526*lmMst2)*pow4(Mst1) +
        27*(1387 - 70*lmMst1 + 70*lmMst2)*pow4(Mst2)))/405.)/pow5(Mst2) + (
        pow2(Mst1)*(628.3249657064472 - (1352*B4)/3. + (28*DN)/3. - (43540*
        lmMst1)/81. + (128*pow2(lmMst1))/9. - (4*lmMst2*(11213 + 1368*lmMst1 +
        144*pow2(lmMst1)))/81. + (8*(-214 + 523*lmMst1)*pow2(lmMst2))/9. - (
        4120*pow3(lmMst2))/9. + (160*Dmsqst2*(-8 + 15*lmMst1 - 15*lmMst2)*(
        pow3(Dmsqst2) + pow6(Msq)))/(3.*pow8(Msq))))/pow3(Mst2) + (
        535.2578189300411 - (1960*B4)/9. + (44*DN)/9. - (3296*lmMst1)/9. - (
        128*pow2(lmMst1))/9. + (8*lmMst2*(-599 - 108*lmMst1 + 24*pow2(lmMst1)))
        /27. + (8*(-222 + 217*lmMst1)*pow2(lmMst2))/9. - 200*pow3(lmMst2) + (
        160*Dmsqst2*(2 + 5*lmMst1 - 5*lmMst2)*(pow3(Dmsqst2) + pow6(Msq)))/(3.*
        pow8(Msq)))/Mst2) + pow2(Mt)*((4*T1ep*(30*Mst2*pow2(Mst1)*(311*Mst2*Mt*
        s2t - 42*pow2(Mt) + 37*pow2(Mst2)*pow2(s2t)) + 2*s2t*(-16421*Mt + 8654*
        Mst2*s2t)*pow4(Mst1) + 63*s2t*(Mt + 5*Mst2*s2t)*pow4(Mst2)))/(729.*
        pow5(Mst2)) + pow2(s2t)*(90.77757201646091 - (100*B4)/9. + (112*D3)/9.
         - (62*DN)/
        9. + (788*lmMst1)/9. - (103*pow2(lmMst1))/9. + (16*lmMst2*(-241 - 171*
        lmMst1 + 18*pow2(lmMst1)))/27. + (125.33333333333333 - 82*lmMst1)*pow2(
        lmMst2) + (214*pow3(lmMst2))/3. - (32*(-2 + lmMst2 + 5*pow2(lmMst2))*
        pow4(Mst2))/(9.*pow4(Mst1)) + ((585.1892843532082 - (8*B4)/3. + (32*D3)
        /9. - (20*DN)/9. - (20109937*lmMst1)/297675. - (15886*pow2(lmMst1))/
        945. + lmMst2*(17.112243218274966 - (144628*lmMst1)/945. + (1180*pow2(
        lmMst1))/9.) + (2*(80257 - 138180*lmMst1)*pow2(lmMst2))/945. - (92*
        pow3(lmMst1))/27. + (4448*pow3(lmMst2))/27.)*pow4(Mst1) - (8*OepS2*(
        1110*pow2(Mst1)*pow2(Mst2) + 17308*pow4(Mst1) + 315*pow4(Mst2)))/2187.
         - (S2*(6*(1089707 - 5550*lmMst1 + 5550*lmMst2)*pow2(Mst1)*pow2(Mst2) +
        40*(93919 - 12981*lmMst1 + 12981*lmMst2)*pow4(Mst1) + 27*(116129 - 350*
        lmMst1 + 350*lmMst2)*pow4(Mst2)))/405.)/pow4(Mst2) + (pow2(Mst1)*(
        155.38193689986284 - (164*B4)/9. + (176*D3)/9. - (94*DN)/9. - (21449*
        lmMst1)/675. - (7556*pow2(lmMst1))/135. + (lmMst2*(-54451 - 22990*
        lmMst1 + 65550*pow2(lmMst1)))/675. + (2*(6077 - 17130*lmMst1)*pow2(
        lmMst2))/135. - (10*pow3(lmMst1))/27. + (4240*pow3(lmMst2))/27. + (10*
        Dmsqst2*(-23 + 42*lmMst1 - 42*lmMst2)*(pow3(Dmsqst2) + pow6(Msq)))/(3.*
        pow8(Msq))))/pow2(Mst2) + ((80*Dmsqst2*(4 + lmMst1 - lmMst2)*(pow3(
        Dmsqst2) + pow6(Msq)))/3. + (pow2(Mst2)*(360*pow4(Dmsqst2) + 360*
        Dmsqst2*pow6(Msq) + (-852 + 358*lmMst2 + 32*lmMst1*(-2 + 3*lmMst2) +
        497*pow2(lmMst2))*pow8(Msq)))/(9.*pow2(Mst1)))/pow8(Msq))))))/Tbeta) -
        (Mt*MuSUSY*xDR2DRMOD*(54*Mst2*s2t*pow2(Mst1)*pow4(Msq)*(-5*Al4p*
        threeLoopFlag*xDmsqst2*pow2(Dmsqst2)*(3*xDmglst2*pow2(Dmglst2)*pow2(
        Mst2)*(4*Tbeta*pow2(Mst2)*(-pow2(Mst1) + pow2(Mst2))*pow2(Mt) + Tbeta*
        pow2(Mst2)*pow2(s2t)*(2*(lmMst1 - lmMst2)*pow2(Mst1)*pow2(Mst2) + pow4(
        Mst1) - pow4(Mst2)) + 2*Mt*MuSUSY*s2t*(pow2(Mst1)*pow2(Mst2) - 2*
        lmMst1*pow2(Mst1)*pow2(Mst2) + 2*lmMst2*pow2(Mst1)*pow2(Mst2) + 20*(
        lmMst1 - lmMst2)*pow4(Mst1) + pow4(Mst2))) + Mgl*((10*Dmglst2 + Mgl*(11
        - 18*shiftst1))*Tbeta*(-4*pow2(Mst1)*pow2(Mt) + 4*pow2(Mst2)*pow2(Mt) +
        2*(lmMst1 - lmMst2)*pow2(Mst1)*pow2(Mst2)*pow2(s2t) + pow2(s2t)*(pow4(
        Mst1) - pow4(Mst2)))*pow4(Mst2) + 2*Mt*MuSUSY*s2t*((10*Dmglst2 + Mgl*(
        11 - 18*shiftst1))*(pow2(Mst1) + pow2(Mst2))*pow4(Mst2) + 2*(lmMst1 -
        lmMst2)*pow2(Mst1)*(2*Dmglst2*(2*pow2(Mst1) - 5*pow2(Mst2))*pow2(Mst2)
        + Mgl*((-11 + 18*shiftst1)*pow2(Mst1)*pow2(Mst2) + 2*(-5 + 9*shiftst1)*
        pow4(Mst1) + (-11 + 18*shiftst1)*pow4(Mst2)))))) + 4*Mgl*(2*Dmglst2*
        lmMst2 + Mgl + lmMst2*Mgl)*twoLoopFlag*pow4(Msq)*(-((2*Mt*(MuSUSY*s2t -
        2*Mt*Tbeta)*pow2(Mst1) + pow2(Mst2)*(2*Mt*MuSUSY*s2t + 4*Tbeta*pow2(Mt)
        - Tbeta*pow2(Mst2)*pow2(s2t)) + Tbeta*pow2(s2t)*pow4(Mst1))*pow4(Mst2))
        + 2*(lmMst1 - lmMst2)*s2t*pow2(Mst1)*(2*Mt*MuSUSY*(pow2(Mst1)*pow2(
        Mst2) + pow4(Mst1) + pow4(Mst2)) - s2t*Tbeta*pow6(Mst2)))) - xDmglst2*
        pow2(Dmglst2)*(-216*(2 - 3*lmMst2)*Mst2*s2t*twoLoopFlag*pow2(Mst1)*(
        Tbeta*(-4*pow2(Mst1)*pow2(Mt) + 4*pow2(Mst2)*pow2(Mt) + 2*(lmMst1 -
        lmMst2)*pow2(Mst1)*pow2(Mst2)*pow2(s2t) + pow2(s2t)*(pow4(Mst1) - pow4(
        Mst2)))*pow4(Mst2) + 2*Mt*MuSUSY*s2t*((pow2(Mst1) + pow2(Mst2))*pow4(
        Mst2) - 2*(lmMst1 - lmMst2)*pow2(Mst1)*(pow2(Mst1)*pow2(Mst2) + pow4(
        Mst1) + pow4(Mst2))))*pow8(Msq) + Al4p*threeLoopFlag*(-128*Tbeta*pow2(
        Mst1)*(-9*Mt*pow2(Mst2)*pow2(s2t)*(6*pow2(Mst2)*(-7 + 30*lmMst2 + 2*
        pow2(lmMst2) + lmMst1*(-11 - 2*lmMst2 + 12*pow2(lmMst2)) - 12*pow3(
        lmMst2))*pow4(Mst1) + pow2(Mst1)*(29 + 179*lmMst2 - 18*pow2(lmMst2) +
        6*lmMst1*(-10 - 3*lmMst2 + 12*pow2(lmMst2)) - 72*pow3(lmMst2))*pow4(
        Mst2) + 3*(17 + 51*lmMst2 - 4*pow2(lmMst2) + 4*lmMst1*(-6 + lmMst2 + 3*
        pow2(lmMst2)) - 12*pow3(lmMst2))*pow6(Mst1) + 12*(2 + lmMst2 - 3*pow2(
        lmMst2))*pow6(Mst2)) + 4*pow3(Mt)*(-(pow2(Mst2)*(834 + lmMst2*(569 -
        33*lmMt) - 162*lmMt + (51 - 18*lmMt)*pow2(lmMst2) - 18*lmMst1*(13 +
        lmMst2 - lmMst2*lmMt + pow2(lmMst2)) + 18*pow3(lmMst2))*pow4(Mst1)) -
        pow2(Mst1)*(96 + 36*lmMt + lmMst2*(-187 + 39*lmMt) - 3*(19 + 6*lmMt)*
        pow2(lmMst2) - 18*lmMst1*(-3 + 2*lmMt - lmMst2*(3 + lmMt) + pow2(
        lmMst2)) + 18*pow3(lmMst2))*pow4(Mst2) + (290 + lmMst2*(2447 - 645*
        lmMt) - 375*lmMt - 3*(-509 + 60*lmMt)*pow2(lmMst2) - 18*lmMst1*(87 +
        lmMst2*(71 - 10*lmMt) - 22*lmMt + 10*pow2(lmMst2)) + 180*pow3(lmMst2))*
        pow6(Mst1) + 36*(2 + lmMst2 - 3*pow2(lmMst2))*pow6(Mst2)))*pow8(Msq) +
        3*Tbeta*pow3(Mst2)*pow3(s2t)*(-3*pow2(Mst1)*pow6(Mst2)*(585*pow4(
        Dmsqst2) - 900*Dmsqst2*pow6(Msq) + 4*(-484 + 390*lmMst2 + 32*lmMst1*(-2
        + 3*lmMst2) + 657*pow2(lmMst2))*pow8(Msq)) + 2*pow2(Mst2)*pow6(Mst1)*(
        3105*pow4(Dmsqst2) - 1350*Dmsqst2*pow6(Msq) + 2*(260 + 1538*lmMst2 +
        768*(2 - 3*lmMst2)*pow2(lmMst1) + 1395*pow2(lmMst2) + 96*lmMst1*(-7 -
        27*lmMst2 + 48*pow2(lmMst2)) - 2304*pow3(lmMst2))*pow8(Msq)) - 4*pow4(
        Mst1)*pow4(Mst2)*(-675*(lmMst1 - lmMst2)*(pow4(Dmsqst2) - 2*Dmsqst2*
        pow6(Msq)) + 2*(40 - 1172*lmMst2 + 48*(-2 + 3*lmMst2)*pow2(lmMst1) +
        882*pow2(lmMst2) - 3*lmMst1*(-500 + 502*lmMst2 + 369*pow2(lmMst2)) +
        963*pow3(lmMst2))*pow8(Msq)) + 384*pow8(Msq)*((2*lmMst1*(3 + 2*lmMst2 +
        16*pow2(lmMst2)) - lmMst2*(7 + 4*lmMst2 + 16*(pow2(lmMst1) + pow2(
        lmMst2))))*pow8(Mst1) + (-2 + lmMst2 + 5*pow2(lmMst2))*pow8(Mst2))) +
        12*Mt*s2t*(-2*Mst2*Mt*Tbeta*(-2*(860 - 434*lmMst2 + 96*(-2 + 3*lmMst2)*
        pow2(lmMst1) + 96*lmMst1*(3 + lmMst2 - 6*pow2(lmMst2)) - 243*pow2(
        lmMst2) + 288*pow3(lmMst2))*pow4(Mst1)*pow4(Mst2)*pow8(Msq) + 192*pow2(
        Mst2)*((2 - 3*lmMst2)*pow2(lmMst1) + lmMst1*(8 - 2*lmMst2 + 6*pow2(
        lmMst2)) - 3*(6 + 5*lmMst2 + pow3(lmMst2)))*pow6(Mst1)*pow8(Msq) + 3*
        pow2(Mst1)*pow4(Mst2)*(45*(23*pow2(Mst1) - 5*pow2(Mst2))*pow4(Dmsqst2)
        - 450*Dmsqst2*(pow2(Mst1) - pow2(Mst2))*pow6(Msq) - 2*(-324 + 134*
        lmMst2 + 32*lmMst1*(-2 + 3*lmMst2) + 177*pow2(lmMst2))*pow2(Mst2)*pow8(
        Msq)) - 384*(-6 - 14*lmMst2 + lmMst2*pow2(lmMst1) + lmMst1*(9 + 6*
        lmMst2 - 2*pow2(lmMst2)) - 6*pow2(lmMst2) + pow3(lmMst2))*pow8(Msq)*
        pow8(Mst1) + 192*(-2 + lmMst2 + 5*pow2(lmMst2))*pow8(Msq)*pow8(Mst2)) +
        MuSUSY*(135*Dmsqst2*s2t*pow2(Mst1)*pow3(Mst2)*(pow3(Dmsqst2)*((23 - 46*
        lmMst1 + 46*lmMst2)*pow2(Mst1)*pow2(Mst2) + 20*(lmMst1 - lmMst2)*pow4(
        Mst1) + 5*pow4(Mst2)) + 10*(2*(lmMst1 - lmMst2)*pow2(Mst1) - pow2(Mst2)
        )*(pow2(Mst1) + pow2(Mst2))*pow6(Msq)) + 2*pow8(Msq)*(-((Mst2*s2t*(1180
        + 1270*lmMst2 + 96*(2 - 3*lmMst2)*pow2(lmMst1) - 3255*pow2(lmMst2) + 6*
        lmMst1*(-468 + 454*lmMst2 + 369*pow2(lmMst2)) - 1926*pow3(lmMst2)) +
        64*Mt*(53 + 191*lmMst2 - 54*pow2(lmMst2) + 6*lmMst1*(-10 - 3*lmMst2 +
        12*pow2(lmMst2)) - 72*pow3(lmMst2)))*pow4(Mst1)*pow4(Mst2)) - 2*pow2(
        Mst2)*(3*Mst2*s2t*(240 + 468*lmMst2 + 144*(2 - 3*lmMst2)*pow2(lmMst1) -
        310*pow2(lmMst2) + lmMst1*(-580 + 22*lmMst2 + 1137*pow2(lmMst2)) - 705*
        pow3(lmMst2)) + 32*Mt*(11 + 371*lmMst2 - 42*pow2(lmMst2) + 6*lmMst1*(-
        21 - 5*lmMst2 + 24*pow2(lmMst2)) - 144*pow3(lmMst2)))*pow6(Mst1) - 3*(
        Mst2*s2t*(420 + lmMst1*(64 - 96*lmMst2) - 358*lmMst2 - 497*pow2(lmMst2)
        ) + 256*Mt*(2 + lmMst2 - 3*pow2(lmMst2)))*pow2(Mst1)*pow6(Mst2) - 12*(
        16*Mt*(49 + 103*lmMst2 - 36*(1 + lmMst2)*pow2(lmMst2) + 4*lmMst1*(-11 +
        6*lmMst2 + 9*pow2(lmMst2))) - Mst2*s2t*(48 + 4*lmMst2*(45 + 68*pow2(
        lmMst1)) + 310*pow2(lmMst2) - lmMst1*(92 + 310*lmMst2 + 635*pow2(
        lmMst2)) + 363*pow3(lmMst2)))*pow8(Mst1) - 96*s2t*(-2 + lmMst2 + 5*
        pow2(lmMst2))*pow9(Mst2))))))))/(648.*Tbeta*pow4(Mst1)*pow5(Mst2)*pow8(
        Msq))))/pow2(Mgl);
}

/**
 * @return returns the susy log^0 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double H6bq2g2::calc_coef_at_as2_no_sm_logs_log0() const {

   const double result =
      (-55566000*pow3(log(pow2(Mst1)/pow2(Mst2)))*pow4(Msq)*pow4(Mst1)*(4*
        Dmglst2*Mgl*(pow4(Mst1)*(-2896*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 2016*
        Mst2*s2t*pow3(Mt) - 672*Mt*pow3(Mst2)*pow3(s2t) + 7424*pow4(Mt) - 33*
        pow4(Mst2)*pow4(s2t)) - 2*pow2(Mst1)*pow2(Mst2)*(1156*pow2(Mst2)*pow2(
        Mt)*pow2(s2t) + 384*Mst2*s2t*pow3(Mt) + 390*Mt*pow3(Mst2)*pow3(s2t) -
        2632*pow4(Mt) + 29*pow4(Mst2)*pow4(s2t)) - 3*pow4(Mst2)*(256*pow2(Mst2)
        *pow2(Mt)*pow2(s2t) - 784*Mst2*s2t*pow3(Mt) + 516*Mt*pow3(Mst2)*pow3(
        s2t) - 512*pow4(Mt) + 107*pow4(Mst2)*pow4(s2t))) + pow2(Mgl)*(pow4(
        Mst1)*(-2696*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 7808*Mst2*s2t*pow3(Mt) -
        3200*Mt*pow3(Mst2)*pow3(s2t) + 6400*pow4(Mt) - 409*pow4(Mst2)*pow4(s2t)
        ) - 4*pow2(Mst1)*pow2(Mst2)*(1166*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 2304*
        Mst2*s2t*pow3(Mt) + 780*Mt*pow3(Mst2)*pow3(s2t) - 440*pow4(Mt) + 115*
        pow4(Mst2)*pow4(s2t)) - 3*pow4(Mst2)*(56*pow2(Mst2)*pow2(Mt)*pow2(s2t)
        - 3136*Mst2*s2t*pow3(Mt) + 1040*Mt*pow3(Mst2)*pow3(s2t) + 32*pow4(Mt) +
        271*pow4(Mst2)*pow4(s2t))) + 2*pow2(Dmglst2)*(2*pow4(Mst1)*(-2896*pow2(
        Mst2)*pow2(Mt)*pow2(s2t) - 2016*Mst2*s2t*pow3(Mt) - 672*Mt*pow3(Mst2)*
        pow3(s2t) + 7424*pow4(Mt) - 33*pow4(Mst2)*pow4(s2t)) - 2*pow2(Mst1)*
        pow2(Mst2)*(3452*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 768*Mst2*s2t*pow3(Mt)
        + 780*Mt*pow3(Mst2)*pow3(s2t) - 440*pow4(Mt) + 97*pow4(Mst2)*pow4(s2t))
        - 3*pow4(Mst2)*(768*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 1568*Mst2*s2t*pow3(
        Mt) + 1800*Mt*pow3(Mst2)*pow3(s2t) - 1536*pow4(Mt) + 321*pow4(Mst2)*
        pow4(s2t)))) + 120*Dmsqst2*pow2(Msq)*pow2(Mst1)*(-1543500*Dmglst2*Mgl*
        Mst2*(4*pow2(Mst1)*pow3(Mst2)*(36*(23 - 24*z2)*pow2(Mst2)*pow2(Mt)*
        pow2(s2t) + 24*Mst2*s2t*(289 - 342*z2 + 216*z3)*pow3(Mt) - 216*Mt*(-2 +
        2*z2 - z3)*pow3(Mst2)*pow3(s2t) + 8*(-1477 + 1242*z2 - 540*z3)*pow4(Mt)
        + 27*(-2 + 4*z2 - 3*z3)*pow4(Mst2)*pow4(s2t)) - 9*Mst2*pow4(Mst1)*(8*(-
        113 + 48*z2)*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 192*Mst2*s2t*(51*z2 - 28*(
        2 + z3))*pow3(Mt) - 96*Mt*(5*z2 + 3*(-4 + z3))*pow3(Mst2)*pow3(s2t) -
        48*(-401 + 308*z2 - 144*z3)*pow4(Mt) + 3*(-75 + 32*z2 + 8*z3)*pow4(
        Mst2)*pow4(s2t)) + 108*(-3 + z2)*pow2(-4*pow2(Mt) + pow2(Mst2)*pow2(
        s2t))*pow5(Mst2) - 6*s2t*(Mst2*s2t*(1931 - 1098*z2 - 180*z3)*pow2(Mt) -
        72*Mt*(-20 + 11*z2)*pow2(Mst2)*pow2(s2t) + 72*(-465 + 348*z2 - 176*z3)*
        pow3(Mt) + (38 + 63*z2 - 90*z3)*pow3(Mst2)*pow3(s2t))*pow6(Mst1)) -
        1543500*Mst2*pow2(Dmglst2)*(4*pow2(Mst1)*pow3(Mst2)*(36*(23 - 24*z2)*
        pow2(Mst2)*pow2(Mt)*pow2(s2t) + 24*Mst2*s2t*(289 - 342*z2 + 216*z3)*
        pow3(Mt) - 216*Mt*(-2 + 2*z2 - z3)*pow3(Mst2)*pow3(s2t) + 8*(-1477 +
        1242*z2 - 540*z3)*pow4(Mt) + 27*(-2 + 4*z2 - 3*z3)*pow4(Mst2)*pow4(s2t)
        ) - 9*Mst2*pow4(Mst1)*(8*(-113 + 48*z2)*pow2(Mst2)*pow2(Mt)*pow2(s2t) +
        192*Mst2*s2t*(51*z2 - 28*(2 + z3))*pow3(Mt) - 96*Mt*(5*z2 + 3*(-4 + z3)
        )*pow3(Mst2)*pow3(s2t) - 48*(-401 + 308*z2 - 144*z3)*pow4(Mt) + 3*(-75
        + 32*z2 + 8*z3)*pow4(Mst2)*pow4(s2t)) + 108*(-3 + z2)*pow2(-4*pow2(Mt)
        + pow2(Mst2)*pow2(s2t))*pow5(Mst2) - 6*s2t*(Mst2*s2t*(1931 - 1098*z2 -
        180*z3)*pow2(Mt) - 72*Mt*(-20 + 11*z2)*pow2(Mst2)*pow2(s2t) + 72*(-465
        + 348*z2 - 176*z3)*pow3(Mt) + (38 + 63*z2 - 90*z3)*pow3(Mst2)*pow3(s2t)
        )*pow6(Mst1)) + pow2(Mgl)*(3087*pow2(Mst2)*pow4(Mst1)*(-8*pow2(Mst2)*
        pow2(Mt)*pow2(s2t)*(36449 + 4000*OepS2 - 310500*S2 - 6000*T1ep - 75750*
        z2 + 21000*z3 - 3000*z4 - 4500*pow2(z2)) - 864000*Mst2*s2t*(-6 + 5*z2 -
        4*z3)*pow3(Mt) - 432000*Mt*(1 + z2 - z3)*pow3(Mst2)*pow3(s2t) + 96*(-
        91321 + 54000*z2 - 36000*z3)*pow4(Mt) + (52199 + 28000*OepS2 - 3415500*
        S2 - 42000*T1ep - 557250*z2 + 259500*z3 - 21000*z4 - 31500*pow2(z2))*
        pow4(Mst2)*pow4(s2t)) + 1157625*pow2(Mst1)*pow4(Mst2)*(-4*pow2(Mst2)*
        pow2(Mt)*pow2(s2t)*(32*OepS2 - 3*(391 + 972*S2 + 16*T1ep + 154*z2 -
        232*z3 + 8*z4 + 12*pow2(z2))) - 1152*Mst2*s2t*(7 + 2*z2 - 8*z3)*pow3(
        Mt) - 1152*Mt*(-1 + 2*z2 - 3*z3)*pow3(Mst2)*pow3(s2t) + 64*(-103 + 108*
        z2 - 72*z3)*pow4(Mt) + (-1213 + 32*OepS2 + 4860*S2 - 48*T1ep - 462*z2 -
        276*z3 - 24*z4 - 36*pow2(z2))*pow4(Mst2)*pow4(s2t)) + 2*(-4*pow2(Mst2)*
        pow2(Mt)*pow2(s2t)*(-22430251 + 5488000*OepS2 - 324135000*S2 - 8232000*
        T1ep - 103929000*z2 + 28812000*z3 - 4116000*z4 - 6174000*pow2(z2)) -
        333396000*Mst2*s2t*(-35 + 20*z2 - 16*z3)*pow3(Mt) + 333396000*Mt*(-2 +
        z2)*pow3(Mst2)*pow3(s2t) + 3456*(-6582784 + 2701125*z2 - 2315250*z3)*
        pow4(Mt) + (378483467 + 9604000*OepS2 - 754771500*S2 - 14406000*T1ep -
        306899250*z2 + 133770000*z3 - 7203000*z4 - 10804500*pow2(z2))*pow4(
        Mst2)*pow4(s2t))*pow6(Mst1) + 41674500*(1 + 2*z2)*pow2(-4*pow2(Mt) +
        pow2(Mst2)*pow2(s2t))*pow6(Mst2))) + 30*pow2(Dmsqst2)*pow2(Mst1)*(-
        6174000*Dmglst2*Mgl*Mst2*(4*pow2(Mst1)*pow3(Mst2)*(36*(23 - 24*z2)*
        pow2(Mst2)*pow2(Mt)*pow2(s2t) + 24*Mst2*s2t*(289 - 342*z2 + 216*z3)*
        pow3(Mt) - 216*Mt*(-2 + 2*z2 - z3)*pow3(Mst2)*pow3(s2t) + 8*(-1477 +
        1242*z2 - 540*z3)*pow4(Mt) + 27*(-2 + 4*z2 - 3*z3)*pow4(Mst2)*pow4(s2t)
        ) - 9*Mst2*pow4(Mst1)*(8*(-113 + 48*z2)*pow2(Mst2)*pow2(Mt)*pow2(s2t) +
        192*Mst2*s2t*(51*z2 - 28*(2 + z3))*pow3(Mt) - 96*Mt*(5*z2 + 3*(-4 + z3)
        )*pow3(Mst2)*pow3(s2t) - 48*(-401 + 308*z2 - 144*z3)*pow4(Mt) + 3*(-75
        + 32*z2 + 8*z3)*pow4(Mst2)*pow4(s2t)) + 108*(-3 + z2)*pow2(-4*pow2(Mt)
        + pow2(Mst2)*pow2(s2t))*pow5(Mst2) - 6*s2t*(Mst2*s2t*(1931 - 1098*z2 -
        180*z3)*pow2(Mt) - 72*Mt*(-20 + 11*z2)*pow2(Mst2)*pow2(s2t) + 72*(-465
        + 348*z2 - 176*z3)*pow3(Mt) + (38 + 63*z2 - 90*z3)*pow3(Mst2)*pow3(s2t)
        )*pow6(Mst1)) - 6174000*Mst2*pow2(Dmglst2)*(4*pow2(Mst1)*pow3(Mst2)*(
        36*(23 - 24*z2)*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 24*Mst2*s2t*(289 - 342*
        z2 + 216*z3)*pow3(Mt) - 216*Mt*(-2 + 2*z2 - z3)*pow3(Mst2)*pow3(s2t) +
        8*(-1477 + 1242*z2 - 540*z3)*pow4(Mt) + 27*(-2 + 4*z2 - 3*z3)*pow4(
        Mst2)*pow4(s2t)) - 9*Mst2*pow4(Mst1)*(8*(-113 + 48*z2)*pow2(Mst2)*pow2(
        Mt)*pow2(s2t) + 192*Mst2*s2t*(51*z2 - 28*(2 + z3))*pow3(Mt) - 96*Mt*(5*
        z2 + 3*(-4 + z3))*pow3(Mst2)*pow3(s2t) - 48*(-401 + 308*z2 - 144*z3)*
        pow4(Mt) + 3*(-75 + 32*z2 + 8*z3)*pow4(Mst2)*pow4(s2t)) + 108*(-3 + z2)
        *pow2(-4*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow5(Mst2) - 6*s2t*(Mst2*s2t*
        (1931 - 1098*z2 - 180*z3)*pow2(Mt) - 72*Mt*(-20 + 11*z2)*pow2(Mst2)*
        pow2(s2t) + 72*(-465 + 348*z2 - 176*z3)*pow3(Mt) + (38 + 63*z2 - 90*z3)
        *pow3(Mst2)*pow3(s2t))*pow6(Mst1)) + pow2(Mgl)*(10290*pow2(Mst2)*pow4(
        Mst1)*(-4*pow2(Mst2)*pow2(Mt)*pow2(s2t)*(-9928 + 1600*OepS2 - 361800*S2
        - 2400*T1ep + 34500*z2 + 50925*z3 - 1200*z4 - 1800*pow2(z2)) - 518400*
        Mst2*s2t*(-23 + 24*z2 - 16*z3)*pow3(Mt) - 259200*Mt*(4 + z2 - 3*z3)*
        pow3(Mst2)*pow3(s2t) + 864*(-23941 + 20700*z2 - 10800*z3)*pow4(Mt) + (
        449186 + 20800*OepS2 - 3439800*S2 - 31200*T1ep - 523500*z2 + 169275*z3
        - 15600*z4 - 23400*pow2(z2))*pow4(Mst2)*pow4(s2t)) + 385875*pow2(Mst1)*
        pow4(Mst2)*(-8*pow2(Mst2)*pow2(Mt)*pow2(s2t)*(-6878 + 128*OepS2 -
        16848*S2 - 192*T1ep - 696*z2 + 3273*z3 - 96*z4 - 144*pow2(z2)) - 3456*
        Mst2*s2t*(13 + 28*z2 - 48*z3)*pow3(Mt) - 6912*Mt*(-5 + 7*z2 - 7*z3)*
        pow3(Mst2)*pow3(s2t) + 16*(-11134 + 10800*z2 - 5481*z3)*pow4(Mt) + (-
        16678 + 256*OepS2 + 114048*S2 - 384*T1ep - 1392*z2 - 16521*z3 - 192*z4
        - 288*pow2(z2))*pow4(Mst2)*pow4(s2t)) + (8*pow2(Mst2)*pow2(Mt)*pow2(
        s2t)*(-1351886621 + 19208000*OepS2 - 639009000*S2 - 28812000*T1ep +
        365552250*z2 + 225865500*z3 - 14406000*z4 - 21609000*pow2(z2)) -
        2667168000*Mst2*s2t*(-35 + 20*z2 - 16*z3)*pow3(Mt) + 166698000*Mt*(-101
        + 106*z2 - 60*z3)*pow3(Mst2)*pow3(s2t) + 27648*(-6582784 + 2701125*z2 -
        2315250*z3)*pow4(Mt) - 3*(-881319682 + 617400000*S2 + 472311000*z2 -
        461892375*z3)*pow4(Mst2)*pow4(s2t))*pow6(Mst1) + 166698000*(-1 + 3*z2)*
        pow2(-4*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow6(Mst2))) + pow4(Msq)*(4*
        Dmglst2*Mgl*pow2(Mst1)*(-1764*pow2(Mst2)*pow4(Mst1)*(-112*pow2(Mst2)*
        pow2(Mt)*pow2(s2t)*(1216808 + 150000*OepS2 - 376650*S2 - 225000*T1ep -
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
        Mst2*s2t)*pow6(Mst2)) + 3*pow2(Mgl)*(308700*pow4(Mst1)*pow4(Mst2)*(2*
        pow2(Mst2)*pow2(Mt)*pow2(s2t)*(1013263 + 2880*D3 - 2880*DN - 25440*
        OepS2 - 137052*S2 + 38160*T1ep + 35562*z2 + 83160*z3 + 36360*z4 +
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
        pow2(s2t))*pow8(Mst2)) + 4*pow2(Dmglst2)*(-2205*pow4(Mst1)*pow4(Mst2)*(
        360*pow2(Mst2)*pow2(Mt)*pow2(s2t)*(-2818497 + 10080*B4 - 10080*D3 +
        5040*DN - 7840*OepS2 + 1126548*S2 + 11760*T1ep + 739362*z2 + 2068815*z3
        - 39480*z4 + 8820*pow2(z2)) - 40*Mst2*s2t*(-6779161 + 7840*OepS2 +
        3032964*S2 - 11760*T1ep - 15966294*z2 + 27128640*z3 - 731640*z4 - 8820*
        pow2(z2))*pow3(Mt) - 70*Mt*(-2048393 + 1058400*B4 - 23760*DN - 1120*
        OepS2 - 449388*S2 + 1680*T1ep - 5386422*z2 + 13314840*z3 + 149880*z4 -
        1864980*pow2(z2))*pow3(Mst2)*pow3(s2t) + 112*(-13292051 + 72800*OepS2 +
        3133890*S2 - 109200*T1ep - 9928650*z2 + 19536225*z3 - 54600*z4 - 81900*
        pow2(z2))*pow4(Mt) + 35*(-1395899 + 54000*B4 - 60480*D3 + 33480*DN +
        5600*OepS2 + 37625796*S2 - 8400*T1ep - 481746*z2 - 6858210*z3 - 266640*
        z4 - 122940*pow2(z2))*pow4(Mst2)*pow4(s2t)) - 294*pow2(Mst2)*(-12*pow2(
        Mst2)*pow2(Mt)*pow2(s2t)*(273721621 + 16716000*OepS2 - 343302300*S2 -
        25074000*T1ep - 610971750*z2 + 71859375*z3 - 12537000*z4 - 18805500*
        pow2(z2)) - 100*Mst2*s2t*(79386499 + 4823840*OepS2 + 82155924*S2 -
        7235760*T1ep - 287072094*z2 + 104826120*z3 - 3617880*z4 - 5426820*pow2(
        z2))*pow3(Mt) - 3500*Mt*(280885 + 11016*B4 + 324*DN - 24544*OepS2 -
        89856*S2 + 36816*T1ep + 858030*z2 - 343020*z3 + 29100*z4 + 27612*pow2(
        z2))*pow3(Mst2)*pow3(s2t) + 8*(55610713 + 131684000*OepS2 + 173875950*
        S2 - 197526000*T1ep - 7894319250*z2 + 4687557000*z3 - 98763000*z4 -
        148144500*pow2(z2))*pow4(Mt) - 7*(-22023067 + 729000*B4 - 972000*D3 +
        607500*DN - 320000*OepS2 - 120274200*S2 + 480000*T1ep - 23129700*z2 +
        92516250*z3 - 3769500*z4 - 4014000*pow2(z2))*pow4(Mst2)*pow4(s2t))*
        pow6(Mst1) + 166698000*pow2(Mst1)*(8*(205 - 12*z2)*pow2(Mst2)*pow2(Mt)*
        pow2(s2t) + 1024*Mst2*s2t*pow3(Mt) - 256*Mt*pow3(Mst2)*pow3(s2t) + 16*(
        -197 + 12*z2)*pow4(Mt) + (-245 + 12*z2)*pow4(Mst2)*pow4(s2t))*pow6(
        Mst2) + (72*pow2(Mst2)*pow2(Mt)*pow2(s2t)*(-3650009897 + 2110136000*
        OepS2 - 52385772600*S2 - 3165204000*T1ep - 40405228500*z2 + 755286000*
        z3 - 1582602000*z4 - 2373903000*pow2(z2)) - 274400*Mst2*s2t*(-20964397
        + 2058400*OepS2 - 47545272*S2 - 3087600*T1ep - 59449002*z2 + 39214920*
        z3 - 1543800*z4 - 2315700*pow2(z2))*pow3(Mt) + 17150*Mt*(-31897243 +
        2491360*OepS2 - 90290268*S2 - 3737040*T1ep - 51499698*z2 + 33912840*z3
        - 1868520*z4 - 2802780*pow2(z2))*pow3(Mst2)*pow3(s2t) + 80*(-
        54946675289 + 3063401600*OepS2 - 52952863320*S2 - 4595102400*T1ep -
        225236187120*z2 + 243414477600*z3 - 2297551200*z4 - 3446326800*pow2(z2)
        )*pow4(Mt) + (119405394763 - 11522056000*OepS2 + 478191735000*S2 +
        17283084000*T1ep + 233890773900*z2 - 111792103500*z3 + 8641542000*z4 +
        12962313000*pow2(z2))*pow4(Mst2)*pow4(s2t))*pow8(Mst1) + 2667168000*
        pow2(-4*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow8(Mst2))) + 1260*log(pow2(
        Mst1)/pow2(Mst2))*(30*pow2(Dmsqst2)*pow2(Mst1)*(-14700*Dmglst2*Mgl*
        Mst2*pow2(Mst1)*(-3*s2t*(2*Mst2*s2t*(37 - 60*z2)*pow2(Mt) + 72*Mt*pow2(
        Mst2)*pow2(s2t) + 1056*(-17 + 8*z2)*pow3(Mt) + (83 - 60*z2)*pow3(Mst2)*
        pow3(s2t))*pow4(Mst1) + 4*pow3(Mst2)*(-48*Mst2*s2t*(-59 + 36*z2)*pow3(
        Mt) + 72*Mt*(-5 + 3*z2)*pow3(Mst2)*pow3(s2t) + 4*(-571 + 360*z2)*pow4(
        Mt) + 9*(2 - 3*z2)*pow4(Mst2)*pow4(s2t)) - 18*Mst2*pow2(Mst1)*(-8*pow2(
        Mst2)*pow2(Mt)*pow2(s2t) + 128*Mst2*s2t*(-12 + 7*z2)*pow3(Mt) + 16*Mt*(
        5 - 3*z2)*pow3(Mst2)*pow3(s2t) - 16*(-137 + 72*z2)*pow4(Mt) + (-13 + 4*
        z2)*pow4(Mst2)*pow4(s2t))) - 14700*Mst2*pow2(Dmglst2)*pow2(Mst1)*(-3*
        s2t*(2*Mst2*s2t*(37 - 60*z2)*pow2(Mt) + 72*Mt*pow2(Mst2)*pow2(s2t) +
        1056*(-17 + 8*z2)*pow3(Mt) + (83 - 60*z2)*pow3(Mst2)*pow3(s2t))*pow4(
        Mst1) + 4*pow3(Mst2)*(-48*Mst2*s2t*(-59 + 36*z2)*pow3(Mt) + 72*Mt*(-5 +
        3*z2)*pow3(Mst2)*pow3(s2t) + 4*(-571 + 360*z2)*pow4(Mt) + 9*(2 - 3*z2)*
        pow4(Mst2)*pow4(s2t)) - 18*Mst2*pow2(Mst1)*(-8*pow2(Mst2)*pow2(Mt)*
        pow2(s2t) + 128*Mst2*s2t*(-12 + 7*z2)*pow3(Mt) + 16*Mt*(5 - 3*z2)*pow3(
        Mst2)*pow3(s2t) - 16*(-137 + 72*z2)*pow4(Mt) + (-13 + 4*z2)*pow4(Mst2)*
        pow4(s2t))) + pow2(Mgl)*(1470*pow2(Mst2)*pow4(Mst1)*(-8*(-443 + 90*S2)*
        pow2(Mst2)*pow2(Mt)*pow2(s2t) - 11520*Mst2*s2t*(-7 + 4*z2)*pow3(Mt) +
        1440*Mt*(-4 + 3*z2)*pow3(Mst2)*pow3(s2t) + 144*(-667 + 360*z2)*pow4(Mt)
        + (-13 + 2340*S2 - 540*z2)*pow4(Mst2)*pow4(s2t)) + 14700*pow2(Mst1)*
        pow4(Mst2)*(-48*(7 + 9*S2)*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 144*Mst2*
        s2t*(-35 + 24*z2)*pow3(Mt) + 144*Mt*(-4 + 3*z2)*pow3(Mst2)*pow3(s2t) +
        28*(-131 + 72*z2)*pow4(Mt) + 3*(35 + 36*S2 - 3*z2)*pow4(Mst2)*pow4(s2t)
        ) + (2*(318121 + 1234800*S2 + 396900*z2)*pow2(Mst2)*pow2(Mt)*pow2(s2t)
        - 4233600*Mst2*s2t*(-25 + 8*z2)*pow3(Mt) - 132300*Mt*(-113 + 60*z2)*
        pow3(Mst2)*pow3(s2t) + 1152*(-175121 + 44100*z2)*pow4(Mt) + 9*(-281161
        + 102900*z2)*pow4(Mst2)*pow4(s2t))*pow6(Mst1) - 264600*pow2(-4*pow2(Mt)
        + pow2(Mst2)*pow2(s2t))*pow6(Mst2))) + 120*Dmsqst2*pow2(Msq)*pow2(Mst1)
        *(-3675*Dmglst2*Mgl*Mst2*pow2(Mst1)*(-3*s2t*(2*Mst2*s2t*(37 - 60*z2)*
        pow2(Mt) + 72*Mt*pow2(Mst2)*pow2(s2t) + 1056*(-17 + 8*z2)*pow3(Mt) + (
        83 - 60*z2)*pow3(Mst2)*pow3(s2t))*pow4(Mst1) + 4*pow3(Mst2)*(-48*Mst2*
        s2t*(-59 + 36*z2)*pow3(Mt) + 72*Mt*(-5 + 3*z2)*pow3(Mst2)*pow3(s2t) +
        4*(-571 + 360*z2)*pow4(Mt) + 9*(2 - 3*z2)*pow4(Mst2)*pow4(s2t)) - 18*
        Mst2*pow2(Mst1)*(-8*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 128*Mst2*s2t*(-12 +
        7*z2)*pow3(Mt) + 16*Mt*(5 - 3*z2)*pow3(Mst2)*pow3(s2t) - 16*(-137 + 72*
        z2)*pow4(Mt) + (-13 + 4*z2)*pow4(Mst2)*pow4(s2t))) - 3675*Mst2*pow2(
        Dmglst2)*pow2(Mst1)*(-3*s2t*(2*Mst2*s2t*(37 - 60*z2)*pow2(Mt) + 72*Mt*
        pow2(Mst2)*pow2(s2t) + 1056*(-17 + 8*z2)*pow3(Mt) + (83 - 60*z2)*pow3(
        Mst2)*pow3(s2t))*pow4(Mst1) + 4*pow3(Mst2)*(-48*Mst2*s2t*(-59 + 36*z2)*
        pow3(Mt) + 72*Mt*(-5 + 3*z2)*pow3(Mst2)*pow3(s2t) + 4*(-571 + 360*z2)*
        pow4(Mt) + 9*(2 - 3*z2)*pow4(Mst2)*pow4(s2t)) - 18*Mst2*pow2(Mst1)*(-8*
        pow2(Mst2)*pow2(Mt)*pow2(s2t) + 128*Mst2*s2t*(-12 + 7*z2)*pow3(Mt) +
        16*Mt*(5 - 3*z2)*pow3(Mst2)*pow3(s2t) - 16*(-137 + 72*z2)*pow4(Mt) + (-
        13 + 4*z2)*pow4(Mst2)*pow4(s2t))) + pow2(Mgl)*(22050*pow2(Mst1)*pow4(
        Mst2)*(-4*(16 + 27*S2)*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 96*Mst2*s2t*(-5
        + 4*z2)*pow3(Mt) + 48*Mt*(-1 + z2)*pow3(Mst2)*pow3(s2t) + 16*(-23 + 12*
        z2)*pow4(Mt) + (17 + 27*S2)*pow4(Mst2)*pow4(s2t)) + 441*pow2(Mst2)*
        pow4(Mst1)*(-8*(-517 + 450*S2)*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 19200*
        Mst2*s2t*(-2 + z2)*pow3(Mt) + 2400*Mt*(-1 + z2)*pow3(Mst2)*pow3(s2t) +
        16*(-3067 + 1200*z2)*pow4(Mt) + (-317 + 3150*S2 - 300*z2)*pow4(Mst2)*
        pow4(s2t)) + (16*(58853 - 44100*S2)*pow2(Mst2)*pow2(Mt)*pow2(s2t) -
        1058400*Mst2*s2t*(-25 + 8*z2)*pow3(Mt) - 264600*Mt*pow3(Mst2)*pow3(s2t)
        + 288*(-175121 + 44100*z2)*pow4(Mt) + (-621671 + 308700*S2 + 132300*z2)
        *pow4(Mst2)*pow4(s2t))*pow6(Mst1) - 66150*pow2(-4*pow2(Mt) + pow2(Mst2)
        *pow2(s2t))*pow6(Mst2))) + pow4(Msq)*(3*pow2(Mgl)*(9800*pow4(Mst1)*
        pow4(Mst2)*(-6*(-4534 + 4293*S2 + 1830*z2 - 72*z3)*pow2(Mst2)*pow2(Mt)*
        pow2(s2t) + 72*Mst2*s2t*(26 + 189*S2 - 758*z2)*pow3(Mt) - 126*Mt*(-208
        + 27*S2 - 2*z2)*pow3(Mst2)*pow3(s2t) + 8*(-553 + 1701*S2 + 4023*z2 +
        108*z3)*pow4(Mt) + 3*(1148 + 1863*S2 + 276*z2 - 18*z3)*pow4(Mst2)*pow4(
        s2t)) + 1960*pow2(Mst2)*(-4*(-52322 + 84915*S2 + 3060*z2 + 540*z3)*
        pow2(Mst2)*pow2(Mt)*pow2(s2t) + 1080*Mst2*s2t*(42 + 377*S2 - 330*z2)*
        pow3(Mt) - 180*Mt*(221 + 219*S2 - 234*z2)*pow3(Mst2)*pow3(s2t) + 8*(-
        22352 + 28215*S2 + 36270*z2)*pow4(Mt) + (-16051 + 86805*S2 - 8325*z2)*
        pow4(Mst2)*pow4(s2t))*pow6(Mst1) + 88200*pow2(Mst1)*(-616*pow2(Mst2)*
        pow2(Mt)*pow2(s2t) - 1024*Mst2*s2t*pow3(Mt) + 256*Mt*pow3(Mst2)*pow3(
        s2t) + 976*pow4(Mt) + 125*pow4(Mst2)*pow4(s2t))*pow6(Mst2) + (-24*(-
        19407863 + 50398950*S2 + 2866500*z2)*pow2(Mst2)*pow2(Mt)*pow2(s2t) +
        940800*Mst2*s2t*(167 + 2055*S2 - 912*z2)*pow3(Mt) - 58800*Mt*(46 +
        1317*S2 - 918*z2)*pow3(Mst2)*pow3(s2t) + 32*(-35585111 + 25754400*S2 +
        32281200*z2)*pow4(Mt) + 3*(-23468297 + 26386500*S2 + 661500*z2 +
        176400*z3)*pow4(Mst2)*pow4(s2t))*pow8(Mst1) + 11289600*pow2(Mt)*(-2*
        pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow8(Mst2)) + 4*Dmglst2*Mgl*(-14700*
        pow4(Mst1)*pow4(Mst2)*(12*(-3427 + 1782*z2)*pow2(Mst2)*pow2(Mt)*pow2(
        s2t) + 4*Mst2*s2t*(-8398 + 2835*S2 - 3258*z2)*pow3(Mt) - 9*Mt*(1640 +
        315*S2 - 1146*z2)*pow3(Mst2)*pow3(s2t) + (26276 - 9072*S2 + 28224*z2)*
        pow4(Mt) + 9*(277 + 63*S2 - 140*z2)*pow4(Mst2)*pow4(s2t)) - 294*pow2(
        Mst2)*(-8*(77657 + 202500*S2 - 55800*z2)*pow2(Mst2)*pow2(Mt)*pow2(s2t)
        + 200*Mst2*s2t*(-11408 + 37341*S2 - 17730*z2)*pow3(Mt) - 300*Mt*(-503 +
        3609*S2 - 1494*z2)*pow3(Mst2)*pow3(s2t) - 48*(8581 + 72450*S2 - 137700*
        z2)*pow4(Mt) + (-81643 + 291600*S2 + 5850*z2)*pow4(Mst2)*pow4(s2t))*
        pow6(Mst1) + 132300*pow2(Mst1)*(-536*pow2(Mst2)*pow2(Mt)*pow2(s2t) -
        512*Mst2*s2t*pow3(Mt) + 128*Mt*pow3(Mst2)*pow3(s2t) + 560*pow4(Mt) +
        131*pow4(Mst2)*pow4(s2t))*pow6(Mst2) + (72*(-4261807 + 33912900*S2 -
        73500*z2)*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 78400*Mst2*s2t*(-4724 +
        115785*S2 - 29484*z2)*pow3(Mt) + 4900*Mt*(18652 + 140139*S2 - 38826*z2)
        *pow3(Mst2)*pow3(s2t) + 80*(42300121 + 49233240*S2 - 64139040*z2)*pow4(
        Mt) + (8287903 - 185175900*S2 + 9525600*z2)*pow4(Mst2)*pow4(s2t))*pow8(
        Mst1) + 16934400*pow2(Mt)*(-2*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow8(
        Mst2)) + 4*pow2(Dmglst2)*(-294*pow4(Mst1)*pow4(Mst2)*(-300*(13339 +
        1134*S2 - 9870*z2)*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 600*Mst2*s2t*(8618 +
        63*S2 - 3354*z2)*pow3(Mt) + 150*Mt*(-3164 + 63*S2 + 3978*z2)*pow3(Mst2)
        *pow3(s2t) + 8*(50551 + 122850*S2 - 149400*z2)*pow4(Mt) + 75*(3050 +
        315*S2 - 222*z2)*pow4(Mst2)*pow4(s2t)) - 147*pow2(Mst2)*(-8*(463453 +
        805950*S2 - 632700*z2)*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 3600*Mst2*s2t*(
        2555 + 4307*S2 - 3542*z2)*pow3(Mt) + 200*Mt*(-3449 + 13806*S2 + 4392*
        z2)*pow3(Mst2)*pow3(s2t) + 16*(-297379 + 2116350*S2 - 1349100*z2)*pow4(
        Mt) + 3*(-165199 + 24000*S2 + 20250*z2)*pow4(Mst2)*pow4(s2t))*pow6(
        Mst1) + 66150*pow2(Mst1)*(-664*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 1024*
        Mst2*s2t*pow3(Mt) - 256*Mt*pow3(Mst2)*pow3(s2t) - 208*pow4(Mt) + 211*
        pow4(Mst2)*pow4(s2t))*pow6(Mst2) + (72*(-4261807 + 33912900*S2 - 73500*
        z2)*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 78400*Mst2*s2t*(-4724 + 115785*S2 -
        29484*z2)*pow3(Mt) + 4900*Mt*(18652 + 140139*S2 - 38826*z2)*pow3(Mst2)*
        pow3(s2t) + 80*(42300121 + 49233240*S2 - 64139040*z2)*pow4(Mt) + (
        8287903 - 185175900*S2 + 9525600*z2)*pow4(Mst2)*pow4(s2t))*pow8(Mst1) +
        1058400*(8*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 16*pow4(Mt) + pow4(Mst2)*
        pow4(s2t))*pow8(Mst2)))) - 1587600*pow2(log(pow2(Mst1)/pow2(Mst2)))*(
        pow4(Msq)*(2*Dmglst2*Mgl*(-140*pow4(Mst1)*pow4(Mst2)*(1806*pow2(Mst2)*
        pow2(Mt)*pow2(s2t) + 1592*Mst2*s2t*pow3(Mt) + 1104*Mt*pow3(Mst2)*pow3(
        s2t) - 3808*pow4(Mt) + 213*pow4(Mst2)*pow4(s2t)) + 7*pow2(Mst2)*(-
        12992*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 175520*Mst2*s2t*pow3(Mt) + 3600*
        Mt*pow3(Mst2)*pow3(s2t) + 285744*pow4(Mt) + 4969*pow4(Mst2)*pow4(s2t))*
        pow6(Mst1) - 105*pow2(Mst1)*(-728*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 1536*
        Mst2*s2t*pow3(Mt) + 384*Mt*pow3(Mst2)*pow3(s2t) + 944*pow4(Mt) + 187*
        pow4(Mst2)*pow4(s2t))*pow6(Mst2) + (-420096*pow2(Mst2)*pow2(Mt)*pow2(
        s2t) - 1672160*Mst2*s2t*pow3(Mt) + 68880*Mt*pow3(Mst2)*pow3(s2t) +
        3253120*pow4(Mt) + 1377*pow4(Mst2)*pow4(s2t))*pow8(Mst1) - 26880*pow2(
        Mt)*(-2*pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow8(Mst2)) - pow2(Mgl)*(70*
        pow4(Mst1)*pow4(Mst2)*(3192*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 14880*Mst2*
        s2t*pow3(Mt) + 5664*Mt*pow3(Mst2)*pow3(s2t) + 7744*pow4(Mt) + 1113*
        pow4(Mst2)*pow4(s2t)) + 560*pow2(Mst2)*(479*pow2(Mst2)*pow2(Mt)*pow2(
        s2t) - 2448*Mst2*s2t*pow3(Mt) + 366*Mt*pow3(Mst2)*pow3(s2t) + 1136*
        pow4(Mt) + 23*pow4(Mst2)*pow4(s2t))*pow6(Mst1) + 105*pow2(Mst1)*(-600*
        pow2(Mst2)*pow2(Mt)*pow2(s2t) - 1024*Mst2*s2t*pow3(Mt) + 256*Mt*pow3(
        Mst2)*pow3(s2t) + 944*pow4(Mt) + 123*pow4(Mst2)*pow4(s2t))*pow6(Mst2) +
        4*(73657*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 316400*Mst2*s2t*pow3(Mt) +
        35000*Mt*pow3(Mst2)*pow3(s2t) + 68678*pow4(Mt) + 11764*pow4(Mst2)*pow4(
        s2t))*pow8(Mst1) + 13440*pow2(Mt)*(-2*pow2(Mt) + pow2(Mst2)*pow2(s2t))*
        pow8(Mst2)) + pow2(Dmglst2)*(14*pow4(Mst1)*pow4(Mst2)*(-28380*pow2(
        Mst2)*pow2(Mt)*pow2(s2t) - 93600*Mst2*s2t*pow3(Mt) - 3600*Mt*pow3(Mst2)
        *pow3(s2t) + 9952*pow4(Mt) + 195*pow4(Mst2)*pow4(s2t)) + 7*pow2(Mst2)*(
        107792*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 47200*Mst2*s2t*pow3(Mt) + 32160*
        Mt*pow3(Mst2)*pow3(s2t) - 695632*pow4(Mt) + 14231*pow4(Mst2)*pow4(s2t))
        *pow6(Mst1) - 105*pow2(Mst1)*(-2696*pow2(Mst2)*pow2(Mt)*pow2(s2t) -
        6144*Mst2*s2t*pow3(Mt) + 1536*Mt*pow3(Mst2)*pow3(s2t) + 2832*pow4(Mt) +
        817*pow4(Mst2)*pow4(s2t))*pow6(Mst2) + 2*(-420096*pow2(Mst2)*pow2(Mt)*
        pow2(s2t) - 1672160*Mst2*s2t*pow3(Mt) + 68880*Mt*pow3(Mst2)*pow3(s2t) +
        3253120*pow4(Mt) + 1377*pow4(Mst2)*pow4(s2t))*pow8(Mst1) + 3360*(-40*
        pow2(Mst2)*pow2(Mt)*pow2(s2t) + 80*pow4(Mt) + 3*pow4(Mst2)*pow4(s2t))*
        pow8(Mst2))) + 90*Dmsqst2*pow2(Msq)*pow4(Mst1)*(2240*Dmglst2*Mgl*Mst2*
        pow2(Mst1)*(-(Mst2*Mt) + 3*s2t*pow2(Mst1))*pow3(Mt) + 2240*Mst2*pow2(
        Dmglst2)*pow2(Mst1)*(-(Mst2*Mt) + 3*s2t*pow2(Mst1))*pow3(Mt) + pow2(
        Mgl)*(2*pow4(Mst1)*(74*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 1120*Mst2*s2t*
        pow3(Mt) + 2736*pow4(Mt) - 53*pow4(Mst2)*pow4(s2t)) - 7*pow2(Mst1)*(-
        48*pow2(Mt)*pow2(s2t)*pow4(Mst2) - 144*pow2(Mst2)*pow4(Mt) + 31*pow4(
        s2t)*pow6(Mst2)) + 140*(-4*pow4(Mst2)*pow4(Mt) + pow4(s2t)*pow8(Mst2)))
        ) + 90*pow2(Dmsqst2)*(2240*Dmglst2*Mgl*Mst2*(-(Mst2*Mt) + 3*s2t*pow2(
        Mst1))*pow3(Mt)*pow6(Mst1) + 2240*Mst2*pow2(Dmglst2)*(-(Mst2*Mt) + 3*
        s2t*pow2(Mst1))*pow3(Mt)*pow6(Mst1) + pow2(Mgl)*pow4(Mst1)*(pow4(Mst1)*
        (-32*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 2240*Mst2*s2t*pow3(Mt) + 5472*
        pow4(Mt) - 31*pow4(Mst2)*pow4(s2t)) - 35*pow2(Mst1)*(-4*pow2(Mt)*pow2(
        s2t)*pow4(Mst2) - 16*pow2(Mst2)*pow4(Mt) + 5*pow4(s2t)*pow6(Mst2)) +
        140*(-2*pow4(Mst2)*pow4(Mt) + pow4(s2t)*pow8(Mst2))))))/(5.00094e8*
        pow2(Mgl)*pow4(Msq)*pow4(Mst1)*pow4(Mst2)*pow4(Mt));

   return result;
}

/**
 * @return returns the susy log^1 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double H6bq2g2::calc_coef_at_as2_no_sm_logs_log1() const {

   const double result =
      (-450*pow2(log(pow2(Mst1)/pow2(Mst2)))*pow4(Msq)*pow4(Mst1)*(pow2(Mgl)*(
        pow4(Mst2)*(-296*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 6784*Mst2*s2t*pow3(Mt)
        - 2208*Mt*pow3(Mst2)*pow3(s2t) + 672*pow4(Mt) - 519*pow4(Mst2)*pow4(
        s2t)) + pow4(Mst1)*(-2688*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 5056*Mst2*
        s2t*pow3(Mt) - 1024*Mt*pow3(Mst2)*pow3(s2t) + 2848*pow4(Mt) - 41*pow4(
        Mst2)*pow4(s2t)) - 8*pow2(Mst1)*pow2(Mst2)*(425*pow2(Mst2)*pow2(Mt)*
        pow2(s2t) - 704*Mst2*s2t*pow3(Mt) + 212*Mt*pow3(Mst2)*pow3(s2t) - 76*
        pow4(Mt) + 3*pow4(Mst2)*pow4(s2t))) + 8*Dmglst2*Mgl*(8*Mt*(-157*Mst2*
        s2t*pow2(Mt) - 108*Mt*pow2(Mst2)*pow2(s2t) + 398*pow3(Mt) - 16*pow3(
        Mst2)*pow3(s2t))*pow4(Mst1) + pow4(Mst2)*(-320*pow2(Mst2)*pow2(Mt)*
        pow2(s2t) + 848*Mst2*s2t*pow3(Mt) - 532*Mt*pow3(Mst2)*pow3(s2t) + 512*
        pow4(Mt) - 99*pow4(Mst2)*pow4(s2t)) + pow2(Mst1)*pow2(Mst2)*(-768*pow2(
        Mst2)*pow2(Mt)*pow2(s2t) - 448*Mst2*s2t*pow3(Mt) - 212*Mt*pow3(Mst2)*
        pow3(s2t) + 1896*pow4(Mt) + 35*pow4(Mst2)*pow4(s2t))) + 4*pow2(Dmglst2)
        *(16*Mt*(-157*Mst2*s2t*pow2(Mt) - 108*Mt*pow2(Mst2)*pow2(s2t) + 398*
        pow3(Mt) - 16*pow3(Mst2)*pow3(s2t))*pow4(Mst1) + pow4(Mst2)*(-960*pow2(
        Mst2)*pow2(Mt)*pow2(s2t) + 1696*Mst2*s2t*pow3(Mt) - 1832*Mt*pow3(Mst2)*
        pow3(s2t) + 1536*pow4(Mt) - 297*pow4(Mst2)*pow4(s2t)) + pow2(Mst1)*
        pow2(Mst2)*(-2304*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 256*Mst2*s2t*pow3(Mt)
        - 424*Mt*pow3(Mst2)*pow3(s2t) + 152*pow4(Mt) + 105*pow4(Mst2)*pow4(s2t)
        ))) + 3000*pow2(Dmsqst2)*pow2(Mst1)*(8*Dmglst2*Mgl*Mst2*pow2(Mst1)*
        pow3(Mt)*(36*Mst2*(Mt*(33 - 18*z2) + 2*Mst2*s2t*(-12 + 7*z2))*pow2(
        Mst1) + (Mt*(283 - 180*z2) + 12*Mst2*s2t*(-29 + 18*z2))*pow3(Mst2) +
        18*s2t*(-83 + 44*z2)*pow4(Mst1)) + 8*Mst2*pow2(Dmglst2)*pow2(Mst1)*
        pow3(Mt)*(36*Mst2*(Mt*(33 - 18*z2) + 2*Mst2*s2t*(-12 + 7*z2))*pow2(
        Mst1) + (Mt*(283 - 180*z2) + 12*Mst2*s2t*(-29 + 18*z2))*pow3(Mst2) +
        18*s2t*(-83 + 44*z2)*pow4(Mst1)) + pow2(Mgl)*(9*pow2(Mst2)*pow4(Mst1)*(
        8*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 32*Mst2*s2t*(-7 + 4*z2)*pow3(Mt) +
        16*(-16 + 9*z2)*pow4(Mt) + pow4(Mst2)*pow4(s2t)) + pow2(Mst1)*pow4(
        Mst2)*(-144*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 72*Mst2*s2t*(-19 + 12*z2)*
        pow3(Mt) + 2*(-545 + 252*z2)*pow4(Mt) + 9*pow4(Mst2)*pow4(s2t)) + 9*(-
        16*Mst2*s2t*(-9 + 4*z2)*pow3(Mt) + 8*(-31 + 12*z2)*pow4(Mt) - pow4(
        Mst2)*pow4(s2t))*pow6(Mst1) - 9*pow2(-4*pow2(Mt) + pow2(Mst2)*pow2(s2t)
        )*pow6(Mst2))) + 3000*Dmsqst2*pow2(Msq)*pow2(Mst1)*(8*Dmglst2*Mgl*Mst2*
        pow2(Mst1)*pow3(Mt)*(36*Mst2*(Mt*(33 - 18*z2) + 2*Mst2*s2t*(-12 + 7*z2)
        )*pow2(Mst1) + (Mt*(283 - 180*z2) + 12*Mst2*s2t*(-29 + 18*z2))*pow3(
        Mst2) + 18*s2t*(-83 + 44*z2)*pow4(Mst1)) + 8*Mst2*pow2(Dmglst2)*pow2(
        Mst1)*pow3(Mt)*(36*Mst2*(Mt*(33 - 18*z2) + 2*Mst2*s2t*(-12 + 7*z2))*
        pow2(Mst1) + (Mt*(283 - 180*z2) + 12*Mst2*s2t*(-29 + 18*z2))*pow3(Mst2)
        + 18*s2t*(-83 + 44*z2)*pow4(Mst1)) + 3*pow2(Mgl)*(3*pow2(Mst2)*pow4(
        Mst1)*(8*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 64*Mst2*s2t*(-2 + z2)*pow3(Mt)
        + 16*(-9 + 4*z2)*pow4(Mt) + pow4(Mst2)*pow4(s2t)) + pow2(Mst1)*pow4(
        Mst2)*(-48*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 96*Mst2*s2t*(-3 + 2*z2)*
        pow3(Mt) + 32*(-8 + 3*z2)*pow4(Mt) + 3*pow4(Mst2)*pow4(s2t)) + 3*(-16*
        Mst2*s2t*(-9 + 4*z2)*pow3(Mt) + 8*(-31 + 12*z2)*pow4(Mt) - pow4(Mst2)*
        pow4(s2t))*pow6(Mst1) - 3*pow2(-4*pow2(Mt) + pow2(Mst2)*pow2(s2t))*
        pow6(Mst2))) + pow4(Msq)*(-100*Dmglst2*Mgl*(2*pow4(Mst1)*pow4(Mst2)*(
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
        Mst2)*pow4(s2t))*pow8(Mst2))) - 30*log(pow2(Mst1)/pow2(Mst2))*(1800*
        pow2(Dmsqst2)*(16*Dmglst2*Mgl*Mst2*(-(Mst2*Mt) + 3*s2t*pow2(Mst1))*
        pow3(Mt)*pow6(Mst1) + 16*Mst2*pow2(Dmglst2)*(-(Mst2*Mt) + 3*s2t*pow2(
        Mst1))*pow3(Mt)*pow6(Mst1) + pow2(Mgl)*pow4(Mst1)*(8*(5*Mt - 2*Mst2*
        s2t)*pow3(Mt)*pow4(Mst1) - 4*pow4(Mst2)*pow4(Mt) + pow2(Mst1)*(4*pow2(
        Mst2)*pow4(Mt) - pow4(s2t)*pow6(Mst2)) + pow4(s2t)*pow8(Mst2))) + 1800*
        Dmsqst2*pow2(Msq)*pow4(Mst1)*(16*Dmglst2*Mgl*Mst2*pow2(Mst1)*(-(Mst2*
        Mt) + 3*s2t*pow2(Mst1))*pow3(Mt) + 16*Mst2*pow2(Dmglst2)*pow2(Mst1)*(-(
        Mst2*Mt) + 3*s2t*pow2(Mst1))*pow3(Mt) + pow2(Mgl)*(8*(5*Mt - 2*Mst2*
        s2t)*pow3(Mt)*pow4(Mst1) - 8*pow4(Mst2)*pow4(Mt) + pow2(Mst1)*(8*pow2(
        Mst2)*pow4(Mt) - pow4(s2t)*pow6(Mst2)) + pow4(s2t)*pow8(Mst2))) + pow4(
        Msq)*(20*Dmglst2*Mgl*(-4*pow4(Mst1)*pow4(Mst2)*(1860*pow2(Mst2)*pow2(
        Mt)*pow2(s2t) + 1472*Mst2*s2t*pow3(Mt) + 774*Mt*pow3(Mst2)*pow3(s2t) -
        2464*pow4(Mt) + 15*pow4(Mst2)*pow4(s2t)) + 2*pow2(Mst2)*(-432*pow2(
        Mst2)*pow2(Mt)*pow2(s2t) - 11656*Mst2*s2t*pow3(Mt) + 996*Mt*pow3(Mst2)*
        pow3(s2t) + 18748*pow4(Mt) + 207*pow4(Mst2)*pow4(s2t))*pow6(Mst1) - 3*
        pow2(Mst1)*(-856*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 1536*Mst2*s2t*pow3(Mt)
        + 384*Mt*pow3(Mst2)*pow3(s2t) + 1200*pow4(Mt) + 203*pow4(Mst2)*pow4(
        s2t))*pow6(Mst2) + 4*(-738*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 10004*Mst2*
        s2t*pow3(Mt) + 135*Mt*pow3(Mst2)*pow3(s2t) + 17532*pow4(Mt) + 12*pow4(
        Mst2)*pow4(s2t))*pow8(Mst1) - 768*pow2(Mt)*(-2*pow2(Mt) + pow2(Mst2)*
        pow2(s2t))*pow8(Mst2)) - 5*pow2(Mgl)*(4*pow4(Mst1)*pow4(Mst2)*(3372*
        pow2(Mst2)*pow2(Mt)*pow2(s2t) - 9024*Mst2*s2t*pow3(Mt) + 3912*Mt*pow3(
        Mst2)*pow3(s2t) + 6112*pow4(Mt) + 579*pow4(Mst2)*pow4(s2t)) - 24*pow2(
        Mst2)*(-534*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 1824*Mst2*s2t*pow3(Mt) -
        60*Mt*pow3(Mst2)*pow3(s2t) - 934*pow4(Mt) + 97*pow4(Mst2)*pow4(s2t))*
        pow6(Mst1) + 6*pow2(Mst1)*(-728*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 1024*
        Mst2*s2t*pow3(Mt) + 256*Mt*pow3(Mst2)*pow3(s2t) + 1200*pow4(Mt) + 139*
        pow4(Mst2)*pow4(s2t))*pow6(Mst2) + 3*(4640*pow2(Mst2)*pow2(Mt)*pow2(
        s2t) - 15296*Mst2*s2t*pow3(Mt) + 240*Mt*pow3(Mst2)*pow3(s2t) + 7128*
        pow4(Mt) - 279*pow4(Mst2)*pow4(s2t))*pow8(Mst1) + 768*pow2(Mt)*(-2*
        pow2(Mt) + pow2(Mst2)*pow2(s2t))*pow8(Mst2)) + 2*pow2(Dmglst2)*(2*pow4(
        Mst1)*pow4(Mst2)*(-40440*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 75360*Mst2*
        s2t*pow3(Mt) - 5760*Mt*pow3(Mst2)*pow3(s2t) + 11872*pow4(Mt) + 4335*
        pow4(Mst2)*pow4(s2t)) + 2*pow2(Mst2)*(29760*pow2(Mst2)*pow2(Mt)*pow2(
        s2t) + 31000*Mst2*s2t*pow3(Mt) + 15840*Mt*pow3(Mst2)*pow3(s2t) -
        228452*pow4(Mt) + 1065*pow4(Mst2)*pow4(s2t))*pow6(Mst1) - 15*pow2(Mst1)
        *(-3080*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 6144*Mst2*s2t*pow3(Mt) + 1536*
        Mt*pow3(Mst2)*pow3(s2t) + 3600*pow4(Mt) + 865*pow4(Mst2)*pow4(s2t))*
        pow6(Mst2) + 40*(-738*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 10004*Mst2*s2t*
        pow3(Mt) + 135*Mt*pow3(Mst2)*pow3(s2t) + 17532*pow4(Mt) + 12*pow4(Mst2)
        *pow4(s2t))*pow8(Mst1) + 480*(-40*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 80*
        pow4(Mt) + 3*pow4(Mst2)*pow4(s2t))*pow8(Mst2)))))/(1350.*pow2(Mgl)*
        pow4(Msq)*pow4(Mst1)*pow4(Mst2)*pow4(Mt));

   return result;
}

/**
 * @return returns the susy log^2 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double H6bq2g2::calc_coef_at_as2_no_sm_logs_log2() const {

   const double result =
      (-30*log(pow2(Mst1)/pow2(Mst2))*pow4(Msq)*pow4(Mst1)*(pow2(Mgl)*(-4*pow4(
        Mst2)*(14*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 456*Mst2*s2t*pow3(Mt) + 146*
        Mt*pow3(Mst2)*pow3(s2t) - 166*pow4(Mt) + 33*pow4(Mst2)*pow4(s2t)) +
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
        pow6(Mst2)))) - 10*Dmglst2*Mgl*pow4(Msq)*(pow4(Mst1)*pow4(Mst2)*(-8208*
        pow2(Mst2)*pow2(Mt)*pow2(s2t) - 4640*Mst2*s2t*pow3(Mt) - 1968*Mt*pow3(
        Mst2)*pow3(s2t) + 4816*pow4(Mt) + 849*pow4(Mst2)*pow4(s2t)) - 3*pow2(
        Mst2)*(-472*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 2496*Mst2*s2t*pow3(Mt) -
        1040*Mt*pow3(Mst2)*pow3(s2t) - 3360*pow4(Mt) + 37*pow4(Mst2)*pow4(s2t))
        *pow6(Mst1) - 3*pow2(Mst1)*(-984*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 1536*
        Mst2*s2t*pow3(Mt) + 384*Mt*pow3(Mst2)*pow3(s2t) + 1456*pow4(Mt) + 219*
        pow4(Mst2)*pow4(s2t))*pow6(Mst2) + 3*(-2240*Mst2*s2t*pow3(Mt) + 3744*
        pow4(Mt) - 59*pow4(Mst2)*pow4(s2t))*pow8(Mst1) - 768*pow2(Mt)*(-2*pow2(
        Mt) + pow2(Mst2)*pow2(s2t))*pow8(Mst2)) + pow2(Dmglst2)*pow4(Msq)*(
        pow4(Mst1)*pow4(Mst2)*(115440*pow2(Mst2)*pow2(Mt)*pow2(s2t) + 98880*
        Mst2*s2t*pow3(Mt) + 19680*Mt*pow3(Mst2)*pow3(s2t) - 29264*pow4(Mt) -
        18495*pow4(Mst2)*pow4(s2t)) + 5*pow2(Mst2)*(-2712*pow2(Mst2)*pow2(Mt)*
        pow2(s2t) - 7680*Mst2*s2t*pow3(Mt) - 8544*Mt*pow3(Mst2)*pow3(s2t) +
        32224*pow4(Mt) + 1101*pow4(Mst2)*pow4(s2t))*pow6(Mst1) + 15*pow2(Mst1)*
        (-3464*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 6144*Mst2*s2t*pow3(Mt) + 1536*
        Mt*pow3(Mst2)*pow3(s2t) + 4368*pow4(Mt) + 913*pow4(Mst2)*pow4(s2t))*
        pow6(Mst2) - 30*(-2240*Mst2*s2t*pow3(Mt) + 3744*pow4(Mt) - 59*pow4(
        Mst2)*pow4(s2t))*pow8(Mst1) - 480*(-40*pow2(Mst2)*pow2(Mt)*pow2(s2t) +
        80*pow4(Mt) + 3*pow4(Mst2)*pow4(s2t))*pow8(Mst2)) + 5*pow2(Mgl)*(720*
        pow2(Dmsqst2)*pow4(Mst1)*pow4(Mst2)*pow4(Mt) + 1440*Dmsqst2*pow2(Msq)*
        pow4(Mst1)*pow4(Mst2)*pow4(Mt) + pow4(Msq)*(pow4(Mst1)*pow4(Mst2)*(
        8592*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 4800*Mst2*s2t*pow3(Mt) + 3936*Mt*
        pow3(Mst2)*pow3(s2t) + 4256*pow4(Mt) - 69*pow4(Mst2)*pow4(s2t)) - 3*
        pow2(Mst2)*(600*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 128*Mst2*s2t*pow3(Mt) +
        1568*Mt*pow3(Mst2)*pow3(s2t) + 672*pow4(Mt) + 355*pow4(Mst2)*pow4(s2t))
        *pow6(Mst1) + 3*pow2(Mst1)*(-856*pow2(Mst2)*pow2(Mt)*pow2(s2t) - 1024*
        Mst2*s2t*pow3(Mt) + 256*Mt*pow3(Mst2)*pow3(s2t) + 1456*pow4(Mt) + 155*
        pow4(Mst2)*pow4(s2t))*pow6(Mst2) + (384*Mst2*s2t*pow3(Mt) - 2400*pow4(
        Mt) + 717*pow4(Mst2)*pow4(s2t))*pow8(Mst1) + 384*pow2(Mt)*(-2*pow2(Mt)
        + pow2(Mst2)*pow2(s2t))*pow8(Mst2))))/(45.*pow2(Mgl)*pow4(Msq)*pow4(
        Mst1)*pow4(Mst2)*pow4(Mt));

   return result;
}

/**
 * @return returns the susy log^3 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
 */
double H6bq2g2::calc_coef_at_as2_no_sm_logs_log3() const {

   const double result =
      -298.6666666666667;

   return result;
}

}        // hierarchies
}        // himalaya
