// ====================================================================
// This file is part of Himalaya.
//
// Himalaya is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================

#include "himalaya/HierarchyCalculator.hpp"
#include "himalaya/HierarchyObject.hpp"
#include "himalaya/version.hpp"

#include "Enums.hpp"
#include "DSZHiggs.hpp"
#include "Flags.hpp"
#include "Hierarchies.hpp"
#include "Mh2EFTCalculator.hpp"
#include "Constants.hpp"
#include "linalg2.hpp"
#include "Logger.hpp"
#include "ThresholdCalculator.hpp"
#include "MSSM_mass_eigenstates.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>

#include <Eigen/Eigenvalues>

/**
 * @file HierarchyCalculator.cpp
 * @brief Implementation of the HierarchyCalculator.
 */

namespace himalaya {

static bool isInfoPrinted; ///< If this bool is true, than no info will be printed in further runs

/**
 * Define static variables
 */
namespace {

/**
 * Returns sqrt of smallest eigenvalue of given 2x2 matrix.
 * @param matrix the matrix whose eigenvalues should be computed
 * @return sqrt of smallest eigenvalue
 */
double calcSmallestEigenvalue(const Eigen::Matrix2d& matrix)
{
   Eigen::EigenSolver<Eigen::Matrix2d> solver(matrix, false);
   const auto eigenvalues = solver.eigenvalues().real().cwiseSqrt();

   return eigenvalues(0) < eigenvalues(1) ? eigenvalues(0) : eigenvalues(1);
}

/// set flags to omit all corrections, except O(at*as^n)
void disable_non_as_terms(himalaya::mh2_eft::Mh2EFTCalculator& mhc)
{
   using namespace himalaya::mh2_eft;

   mhc.setCorrectionFlag(CouplingOrders::G12G22, 0);
   mhc.setCorrectionFlag(CouplingOrders::G12YB2, 0);
   mhc.setCorrectionFlag(CouplingOrders::G14, 0);
   mhc.setCorrectionFlag(CouplingOrders::G24, 0);
   mhc.setCorrectionFlag(CouplingOrders::G12YB2, 0);
   mhc.setCorrectionFlag(CouplingOrders::G22YB2, 0);
   mhc.setCorrectionFlag(CouplingOrders::YB4, 0);
   mhc.setCorrectionFlag(CouplingOrders::G12YTAU2, 0);
   mhc.setCorrectionFlag(CouplingOrders::G22YTAU2, 0);
   mhc.setCorrectionFlag(CouplingOrders::YTAU4, 0);
   mhc.setCorrectionFlag(CouplingOrders::G12YT2, 0);
   mhc.setCorrectionFlag(CouplingOrders::G22YT2, 0);
   mhc.setCorrectionFlag(CouplingOrders::G32YB4, 0);
   mhc.setCorrectionFlag(CouplingOrders::YB6, 0);
   mhc.setCorrectionFlag(CouplingOrders::YT6, 0);
   mhc.setCorrectionFlag(CouplingOrders::YTAU2YB4, 0);
   mhc.setCorrectionFlag(CouplingOrders::YTAU6, 0);
   mhc.setCorrectionFlag(CouplingOrders::YT2YB4, 0);
   mhc.setCorrectionFlag(CouplingOrders::YB2YT4, 0);
   mhc.setCorrectionFlag(CouplingOrders::YTAU4YB2, 0);
}

} // anonymous namespace

/**
 * Constructor
 * @param p_ a HimalayaInterface struct
 * @param verbose_ suppress informative output during the calculation, if set to false
 */
HierarchyCalculator::HierarchyCalculator(const Parameters& p_, bool verbose_)
   : p(p_), verbose(verbose_)
{
   if (!isInfoPrinted && verbose) {
      printInfo();
      isInfoPrinted = true;
   }

   p.validate(verbose);

   expansionDepth.fill(1);
}

/**
 * Calculates the 3-loop mass matrix and other information of the hierarchy selection process.
 * @param isAlphab a bool which determines if the returned object is proportinal to alpha_b.
 * @return A HierarchyObject which holds all information of the calculation.
 */
himalaya::HierarchyObject HierarchyCalculator::calculateDMh3L(bool isAlphab)
{
   HierarchyObject ho (isAlphab);

   if (isAlphab)
      INFO_MSG("3-loop threshold correction Δλ not available for O(αb*αs^2)!");

   const int mdrFlag = 0;

   // set mdrFlag
   ho.setMDRFlag(mdrFlag);

   // compare hierarchies and get the best fitting hierarchy
   compareHierarchies(ho);

   // calculate the 3-loop Higgs mass matrix for the obtained hierachy in the (M)DRbar' scheme
   ho.setDMh(3, calculateHierarchy(ho, 0, 0, 1) + shiftH3mToDRbarPrime(ho));

   // set the alpha_x contributions
   ho.setDMh(1, getMt41L(ho, mdrFlag, mdrFlag));

   // set the alpha_x*alpha_s contributions
   ho.setDMh(2, getMt42L(ho, mdrFlag, mdrFlag));

   // estimate the uncertainty of the expansion at 3-loop level
   ho.setDMhExpUncertainty(3, getExpansionUncertainty(ho,
                                                   ho.getDMh(0) + ho.getDMh(1)
                                                   + ho.getDMh(2), 0, 0, 1));

   // set the uncertainty of the expansion at 1-loop level to 0 by default,
   // if the user needs this value getExpansionUncertainty should be called
   ho.setDMhExpUncertainty(1, 0.);

   // calculate shifts needed to convert DR' to other renormalization schemes,
   // here one needs it with the minus sign to convert DR -> H3m
   ho.setDMhDRbarPrimeToH3mShift(-shiftH3mToDRbarPrime(ho));

   // to obtain delta_lambda one has to divide the difference of the two calculations by v^2
   const double v2 = calcV2();
   const double gt = sqrt2*p.Mt/std::sqrt(v2);

   // calculate delta_lambda @ 3-loop level
   calcDeltaLambda3L(ho, false);

   // create a modified parameters struct and construct
   // Mh2EFTCalculator and ThresholdCalculator
   auto p_mass_ES = p;
   p_mass_ES.mu2(2,2) = pow2(p.MSt(0));
   p_mass_ES.mq2(2,2) = pow2(p.MSt(1));
   himalaya::mh2_eft::Mh2EFTCalculator mh2EFTCalculator(p_mass_ES);
   himalaya::mh2_eft::ThresholdCalculator tc (p_mass_ES);
   disable_non_as_terms(mh2EFTCalculator); // omit all but O(at*as^n)

   // mh_eft^2
   const double mh2_eft = mh2EFTCalculator.getDeltaMh2EFT0Loop();
   // 1-Loop prefactor at
   const double pref_1L = oneLoop * pow2(p.Mt * gt);
   // 2-Loop prefactor at*as
   const double pref_2L = twoLoop * pow2(p.Mt * gt * p.g3);

   ho.setDLambda(0, mh2_eft/v2);
   ho.setDLambda(1, pref_1L*(tc.getThresholdCorrection(
      mh2_eft::ThresholdCouplingOrders::LAMBDA_AT, mh2_eft::RenSchemes::DRBARPRIME, 1))/v2);
   ho.setDLambda(2, pref_2L*(tc.getThresholdCorrection(
      mh2_eft::ThresholdCouplingOrders::LAMBDA_AT_AS, mh2_eft::RenSchemes::DRBARPRIME, 1))/v2);
   ho.setDLambda(3, ho.getDLambdaEFT());
   ho.setDLambdaDRbarPrimeToMSbarShift(0, 0.);
   ho.setDLambdaDRbarPrimeToMSbarShift(1, 0.);
   ho.setDLambdaDRbarPrimeToMSbarShift(2, pref_2L*(-4*ho.getDLambda(1)
      *tc.getThresholdCorrection(mh2_eft::ThresholdCouplingOrders::YT_AS,
                                 mh2_eft::RenSchemes::DRBARPRIME, 1))/v2);

   himalaya::mh2_fo::MSSM_mass_eigenstates mfo(p);
   const auto dmh_fo = mfo.calculate_Mh2();
   ho.setDMh2FO(0, std::get<0>(dmh_fo));
   ho.setDMh2FO(1, std::get<1>(dmh_fo));
   ho.setDMh2FO(2, std::get<2>(dmh_fo));
   ho.setDMh2FO(3, 0); // filled later, see below

   himalaya::mh2_fo::MSSM_mass_eigenstates mfo_atas(p, true);
   const auto dmh_fo_atas = mfo_atas.calculate_Mh2();
   ho.setDMh2FOAt(0, 0.);
   ho.setDMh2FOAt(1, std::get<1>(dmh_fo_atas));
   ho.setDMh2FOAt(2, std::get<2>(dmh_fo_atas));
   ho.setDMh2FOAt(3, 0.); // filled later, see below

   // fill in results of EFT calculation
   himalaya::mh2_eft::Mh2EFTCalculator mh2EFTCalculator_full(p);
   ho.setDMh2EFT(0, mh2EFTCalculator_full.getDeltaMh2EFT0Loop());
   ho.setDMh2EFT(1, mh2EFTCalculator_full.getDeltaMh2EFT1Loop(1, 1));
   ho.setDMh2EFT(2, mh2EFTCalculator_full.getDeltaMh2EFT2Loop(1, 1));
   ho.setDMh2EFT(3, mh2EFTCalculator_full.getDeltaMh2EFT3Loop(1, 1, 0)
      + ho.getDLambdaEFT()*v2);

   {
      auto ho_mdr = ho;
      ho_mdr.setMDRFlag(1);
      // calculate the DR to MDR shift with the obtained hierarchy
      ho_mdr.setDMhDRbarPrimeToMDRbarPrimeShift(calcDRbarToMDRbarShift(ho_mdr, true, true));
      ho_mdr.setDMh(3, calculateHierarchy(ho_mdr, 0, 0, 1) + shiftH3mToDRbarPrime(ho_mdr));
      Eigen::Vector2d mdrMasses;
      mdrMasses(0) = ho_mdr.getMDRMasses()(0);
      mdrMasses(1) = ho_mdr.getMDRMasses()(1);
      ho.setMDRMasses(mdrMasses);
      ho.setDMhDRbarPrimeToMDRbarPrimeShift(ho_mdr.getDMhDRbarPrimeToMDRbarPrimeShift()
                                            + ho_mdr.getDMh(3) - ho.getDMh(3));
   }

   // perturbatively diagonalize
   {
      const auto DMh2 = flexiblesusy::fs_diagonalize_hermitian_perturbatively(
         ho.getDMh(0), ho.getDMh(1), ho.getDMh(2), ho.getDMh(3));
      ho.setDMh2(0, std::get<0>(DMh2)(0));
      ho.setDMh2(1, std::get<1>(DMh2)(0));
      ho.setDMh2(2, std::get<2>(DMh2)(0));
      ho.setDMh2(3, std::get<3>(DMh2)(0));
      // set 3L FO corrections to perturbatively diagonalized O(at*as^n) corrections
      ho.setDMh2FO(3,std::get<3>(DMh2)(0));
      ho.setDMh2FOAt(3,std::get<3>(DMh2)(0));
   }

   return ho;
}

/**
 * Compares deviation of all hierarchies with the exact two-loop result and returns the hierarchy which minimizes the error.
 * @param ho a HierarchyObject with constant isAlphab.
 * @return An integer which is identified with the suitable hierarchy.
 */
int HierarchyCalculator::compareHierarchies(himalaya::HierarchyObject& ho)
{
   using namespace himalaya::hierarchies;

   // set flags to truncate the expansion
   expansionDepth.at(ExpansionDepth::threeLoop) = 0;
   expansionDepth.at(ExpansionDepth::Mst) = 0;
   double error = -1.;
   int suitableHierarchy = -1;

   // sine of 2 times beta
   const double s2b = std::sin(2*calcBeta());
   const double tbeta = calcTanBeta();

   // tree level Higgs mass matrix
   Eigen::Matrix2d treelvl;
   treelvl (0,0) = s2b/2.*(pow2(p.MZ) / tbeta + pow2(p.MA) * tbeta);
   treelvl (1,0) = s2b/2.*(-pow2(p.MZ) - pow2(p.MA));
   treelvl (0,1) = treelvl (1,0);
   treelvl (1,1) = s2b/2.*(pow2(p.MZ) * tbeta + pow2(p.MA) / tbeta);

   ho.setDMh(0, treelvl);

   // compare the exact higgs mass at 2-loop level with the expanded expressions to find a suitable hierarchy
   for (int hierarchy = Hierarchies::FIRST; hierarchy < Hierarchies::NUMBER_OF_HIERARCHIES; hierarchy++) {
      // first, check if the hierarchy is suitable to the mass spectrum
      ho.setSuitableHierarchy(hierarchy);

      if(isHierarchySuitable(ho)){
         // calculate the exact 1-loop result (only alpha_t/b)
         const Eigen::Matrix2d Mt41L = getMt41L(ho, ho.getMDRFlag(), 0);

         // call the routine of Pietro Slavich to get the alpha_s alpha_t/b corrections with the MDRbar masses
         const Eigen::Matrix2d Mt42L = getMt42L(ho, ho.getMDRFlag(), 0);

         // calculate the exact Higgs mass at 2-loop (only up to alpha_s alpha_t/b)
         const double Mh2l = calcSmallestEigenvalue(treelvl + Mt41L + Mt42L);

         // calculate the expanded 2-loop expression with the specific hierarchy
         const double Mh2LExpanded = calcSmallestEigenvalue(treelvl + Mt41L + calculateHierarchy(ho, 0, 1, 0));

         // estimate the error
         const double twoLoopError = std::abs((Mh2l - Mh2LExpanded));

         // estimate the uncertainty of the expansion at 2L
         const double expUncertainty2L = getExpansionUncertainty(ho, treelvl
            + Mt41L, 0, 1, 0);

         // estimate the uncertainty of the expansion at 3L
         const double expUncertainty3L = getExpansionUncertainty(ho, treelvl
            + Mt41L + Mt42L, 0, 0, 1);

         // estimate the uncertainty of the difference of delta_lambda @ 3-loop
         // normalized to the exact logarithms
         calcDeltaLambda3L(ho, true);
         /*const double deltaLambdaUncertainty = std::abs((ho.getDLambdaEFT()
          - ho.getDLambdaH3m())/(ho.getDLambdaEFT() - ho.getDLambdaNonLog()));*/

         // add these errors to include the error of the expansion in the comparison
         const double currError = std::sqrt(pow2(twoLoopError/Mh2l)
            + pow2(expUncertainty2L/Mh2LExpanded));

         // if the error is negative, it is the first iteration and there is no hierarchy which fits better
         if(error < 0){
            error = currError;
            suitableHierarchy = hierarchy;
            ho.setAbsDiff2L(twoLoopError);
            ho.setRelDiff2L(twoLoopError/Mh2l);
            ho.setDMhExpUncertainty(2, expUncertainty2L);
            ho.setDMhExpUncertainty(3, expUncertainty3L);
         }
         // compare the current error with the last error and choose the hierarchy which fits best (lowest error)
         else if(currError < error){
            error = currError;
            suitableHierarchy = hierarchy;
            ho.setAbsDiff2L(twoLoopError);
            ho.setRelDiff2L(twoLoopError/Mh2l);
            ho.setDMhExpUncertainty(2, expUncertainty2L);
            ho.setDMhExpUncertainty(3, expUncertainty3L);
         }
      }
   }
   ho.setSuitableHierarchy(suitableHierarchy);

   // reset the flags
   expansionDepth.at(ExpansionDepth::threeLoop) = 1;
   expansionDepth.at(ExpansionDepth::Mst) = 1;
   return suitableHierarchy;
}

/**
 * Checks if a hierarchy is suitable to the given mass spectrum.
 * @param ho a HierarchyObject with constant isAlphab and a hierarchy candidate.
 * @returns A bool if the hierarchy candidate is suitable to the given mass spectrum.
 */
bool HierarchyCalculator::isHierarchySuitable(const himalaya::HierarchyObject& ho) const
{
   using namespace himalaya::hierarchies;

   double Mst1, Mst2;
   if (!ho.getIsAlphab()) {
      Mst1 = p.MSt(0);
      Mst2 = p.MSt(1);
   } else {
      Mst1 = p.MSb(0);
      Mst2 = p.MSb(1);
   }

   // check if the squark mass is the heaviest mass
   const double delta = 1.3;        // allow for an offset of 30%
   const double Mgl = p.MG;
   const double Msq = calcMeanMsq();

   if(Mst2 > delta*Msq) return false;
   if(Mgl > delta*Msq) return false;

   switch (ho.getSuitableHierarchy()){
      case Hierarchies::h3:
         return Mgl > Mst2;
      case Hierarchies::h32q2g:
         return (Mst2 >= Msq) && (Mst2 > Mgl);
      case Hierarchies::h3q22g:
         return (Msq > Mst2) && (Mst2 > Mgl);
      case Hierarchies::h4:
         return (Mst1 < Msq) && (Mst1 >= Mgl);
      case Hierarchies::h5:
         return (Mst2 - Mst1 > 0.1*Mst1) && ((Mgl - Mst1) < std::abs(Mgl - Mst2)) && (Mst2 < Msq) && (Mst1 >= Mgl);
      case Hierarchies::h5g1:
         return (Mst2 - Mst1 > 0.1*Mst1) && ((Mgl - Mst1) < std::abs(Mgl - Mst2)) && (Mst2 < Msq) && (Mgl > Mst1);
      case Hierarchies::h6:
         return (Mst2 - Mst1 > 0.1*Mst1) && ((Mst2 - Mgl) < std::abs(Mgl - Mst1)) && (Mst2*1.3 < Msq) && (Mst2 >= Mgl);
      case Hierarchies::h6g2:
         return (Mst2 - Mst1 > 0.1*Mst1) && ((Mst2 - Mgl) < std::abs(Mgl - Mst1)) && (Mst2*1.3 < Msq) && (Mgl > Mst2);
      case Hierarchies::h6b:
         return (Mst2 - Mst1 > 0.1*Mst1) && ((Mst2 - Mgl) < std::abs(Mgl - Mst1)) && (Mst2 >= Msq) && (Mst2 >= Mgl);
      case Hierarchies::h6b2qg2:
         return (Mst2 - Mst1 > 0.1*Mst1) && ((Mst2 - Mgl) < std::abs(Mgl - Mst1)) && (Mst2 >= Msq) && (Mgl > Mst2);
      case Hierarchies::h6bq22g:
         return (Mst2 - Mst1 > 0.1*Mst1) && ((Mst2 - Mgl) < std::abs(Mgl - Mst1)) && (Msq > Mst2) && (Mst2 >= Mgl);
      case Hierarchies::h6bq2g2:
         return (Mst2 - Mst1 > 0.1*Mst1) && ((Mst2 - Mgl) < std::abs(Mgl - Mst1)) && (Msq > Mst2) && (Mgl > Mst2);
      case Hierarchies::h9:
         return (Mst2 >= Msq) && ((Mst2 - Mst1) < (Mst1 - Mgl));
      case Hierarchies::h9q2:
         return (Msq > Mst2) && ((Mst1 - Mst1) < (Mst1 - Mgl));
   }
   return false;
}

// TODO(avoigt): if one is interested in the expansion at one- and two-loop choose a unified choice for the MDR scheme
/**
 * Calculates the hierarchy contributions for a specific hierarchy at a specific loop order.
 * @param ho a HierarchyObject with constant isAlphab.
 * @param oneLoopFlagIn an integer flag which is 0 or 1 in order to add or omit the expanded one-loop results to the returned value, respectivley.
 * @param twoLoopFlagIn an integer flag which is 0 or 1 in order to add or omit the expanded two-loop results to the returned value, respectivley.
 * @param threeLoopFlagIn an integer flag which is 0 or 1 in order to add or omit the expanded three-loop results to the returned value, respectivley.
 * @throws runtime_error Throws a runtime_error if the tree-level is requested in terms of hierarchies.
 * @return The loop corrected Higgs mass matrix which contains the expanded corrections at the given order.
 */
Eigen::Matrix2d HierarchyCalculator::calculateHierarchy(
   himalaya::HierarchyObject& ho, int oneLoopFlagIn,
   int twoLoopFlagIn, int threeLoopFlagIn) const
{
   using namespace himalaya::hierarchies;

   // get the hierarchy
   const int hierarchy = ho.getSuitableHierarchy();

   // the hierarchy files containing 1-, 2- and 3-loop terms (alpha_s^0 alpha_t/b, alpha_s alpha_t/b, alpha_s^2 alpha_t/b)
   double sigS1Full = 0., sigS2Full = 0., sigS12Full = 0.;

   // common variables
   double At, Mt, s2t, Mst1 = 0., Mst2 = 0.;
   if (!ho.getIsAlphab()) {
      At = p.Au(2,2);
      Mt = p.Mt;
      s2t = p.s2t;
   } else {
      At = p.Ad(2,2);
      Mt = p.Mb;
      s2t = p.s2b;
   }

   const double beta = calcBeta();
   const double lmMt = std::log(pow2(p.scale / Mt));
   const double Mgl = p.MG;
   const double lmMgl = std::log(pow2(p.scale / p.MG));

   // this loop is needed to calculate the suitable mass shift order by order
   for(int currentLoopOrder = 1; currentLoopOrder <= 3; currentLoopOrder ++){
      bool runThisOrder;
      double curSig1 = 0., curSig2 = 0., curSig12 = 0.;
      int oneLoopFlag = 0, twoLoopFlag = 0, threeLoopFlag = 0;
      switch (currentLoopOrder){
         case 1:
            oneLoopFlag = 1;
            runThisOrder = oneLoopFlag == oneLoopFlagIn;
         break;
         case 2:
            twoLoopFlag = 1;
            runThisOrder = twoLoopFlag == twoLoopFlagIn;
         break;
         case 3:
            threeLoopFlag = 1;
            runThisOrder = threeLoopFlag == threeLoopFlagIn;
         break;
      }
      if(runThisOrder){
         const double Al4p = calcAsOver4Pi();
         const double Msq = calcMeanMsq();
         // set the Msx masses according to MDRFlag
         if(oneLoopFlag == 1){
            Mst1 = shiftMst1ToMDR(ho, 0, 0);
            Mst2 = shiftMst2ToMDR(ho, 0, 0);
         }
         else if(twoLoopFlag == 1){
            Mst1 = shiftMst1ToMDR(ho, ho.getMDRFlag(), 0);
            Mst2 = shiftMst2ToMDR(ho, ho.getMDRFlag(), 0);
         }
         else if(threeLoopFlag == 1){
            Mst1 = shiftMst1ToMDR(ho, ho.getMDRFlag(), ho.getMDRFlag());
            Mst2 = shiftMst2ToMDR(ho, ho.getMDRFlag(), ho.getMDRFlag());
         }
         else{
            throw std::runtime_error("There are no tree-level hierarchies included!");
         }
         // select the suitable hierarchy for the specific hierarchy and set variables
         switch(getMotherHierarchy(hierarchy)){
            case Hierarchies::h3:{
               const double Dmglst1 = Mgl - Mst1;
               const double Dmsqst1 = pow2(Msq) - pow2(Mst1);
               const double Dmst12 = pow2(Mst1) - pow2(Mst2);
               const double lmMst1 = std::log(pow2(p.scale / Mst1));
               switch(hierarchy){
                  case Hierarchies::h3:{
                     const H3 hierarchy3(expansionDepth, Al4p, beta,
                        Dmglst1, Dmst12, Dmsqst1, lmMt, lmMst1,
                        Mgl, Mt, Mst1, Mst2, Msq, p.mu,
                        s2t,
                        ho.getMDRFlag(), oneLoopFlag, twoLoopFlag, threeLoopFlag);
                     curSig1 = hierarchy3.getS1();
                     curSig2 = hierarchy3.getS2();
                     curSig12 = hierarchy3.getS12();
                     if(oneLoopFlagIn == 0 && twoLoopFlagIn == 0 && threeLoopFlagIn == 1){
                        const double c = hierarchy3.calc_coef_at_as2_no_sm_logs_log0();
                        ho.setDLambdaH3m(c
                           + lmMst1 * hierarchy3.calc_coef_at_as2_no_sm_logs_log1()
                           + pow2(lmMst1) * hierarchy3.calc_coef_at_as2_no_sm_logs_log2()
                           + pow3(lmMst1) * hierarchy3.calc_coef_at_as2_no_sm_logs_log3());
                        ho.setDLambdaNonLog(c);
                     }
                  }
                  break;
                  case Hierarchies::h32q2g:{
                     const H32q2g hierarchy32q2g(expansionDepth, Al4p, beta,
                        Dmglst1, Dmst12, Dmsqst1, lmMt, lmMst1,
                        Mt, Mst1, Mst2, p.mu,
                        s2t,
                        ho.getMDRFlag(), oneLoopFlag, twoLoopFlag, threeLoopFlag);
                     curSig1 = hierarchy32q2g.getS1();
                     curSig2 = hierarchy32q2g.getS2();
                     curSig12 = hierarchy32q2g.getS12();
                     if(oneLoopFlagIn == 0 && twoLoopFlagIn == 0 && threeLoopFlagIn == 1){
                        const double c = hierarchy32q2g.calc_coef_at_as2_no_sm_logs_log0();
                        ho.setDLambdaH3m(c
                           + lmMst1 * hierarchy32q2g.calc_coef_at_as2_no_sm_logs_log1()
                           + pow2(lmMst1) * hierarchy32q2g.calc_coef_at_as2_no_sm_logs_log2()
                           + pow3(lmMst1) * hierarchy32q2g.calc_coef_at_as2_no_sm_logs_log3());
                        ho.setDLambdaNonLog(c);
                     }
                  }
                  break;
                  case Hierarchies::h3q22g:{
                     const H3q22g hierarchy3q22g(expansionDepth, Al4p, beta,
                        Dmglst1, Dmst12, Dmsqst1, lmMt, lmMst1,
                        Mt, Mst1, Mst2, Msq, p.mu,
                        s2t,
                        ho.getMDRFlag(), oneLoopFlag, twoLoopFlag, threeLoopFlag);
                     curSig1 = hierarchy3q22g.getS1();
                     curSig2 = hierarchy3q22g.getS2();
                     curSig12 = hierarchy3q22g.getS12();
                     if(oneLoopFlagIn == 0 && twoLoopFlagIn == 0 && threeLoopFlagIn == 1){
                        const double c = hierarchy3q22g.calc_coef_at_as2_no_sm_logs_log0();
                        ho.setDLambdaH3m(c
                           + lmMst1 * hierarchy3q22g.calc_coef_at_as2_no_sm_logs_log1()
                           + pow2(lmMst1) * hierarchy3q22g.calc_coef_at_as2_no_sm_logs_log2()
                           + pow3(lmMst1) * hierarchy3q22g.calc_coef_at_as2_no_sm_logs_log3());
                        ho.setDLambdaNonLog(c);
                     }
                  }
                  break;
               }
            }
            break;
            case Hierarchies::h4:{
               const double Msusy = (Mst1 + Mst2 + Mgl) / 3.;
               const double lmMsusy = std::log(pow2(p.scale / Msusy));
               const double lmMst1 = std::log(pow2(p.scale / Mst1));
               const double lmMsq = std::log(pow2(p.scale / calcMeanMsq()));
               const H4 hierarchy4(expansionDepth, Al4p, At, beta,
                  lmMt, lmMsq, lmMsusy, Mt, Msusy, Msq,
                  ho.getMDRFlag(), oneLoopFlag, twoLoopFlag, threeLoopFlag);
               curSig1 = hierarchy4.getS1();
               curSig2 = hierarchy4.getS2();
               curSig12 = hierarchy4.getS12();
               if(oneLoopFlagIn == 0 && twoLoopFlagIn == 0 && threeLoopFlagIn == 1){
                  const double c = hierarchy4.calc_coef_at_as2_no_sm_logs_log0();
                  ho.setDLambdaH3m(c
                     + lmMst1 * hierarchy4.calc_coef_at_as2_no_sm_logs_log1()
                     + pow2(lmMst1) * hierarchy4.calc_coef_at_as2_no_sm_logs_log2()
                     + pow3(lmMst1) * hierarchy4.calc_coef_at_as2_no_sm_logs_log3());
                  ho.setDLambdaNonLog(c);
               }
            }
            break;
            case Hierarchies::h5:{
               const double Dmglst1 = Mgl - Mst1;
               const double lmMst1 = std::log(pow2(p.scale / Mst1));
               const double lmMst2 = std::log(pow2(p.scale / Mst2));
               const double lmMsq = std::log(pow2(p.scale / calcMeanMsq()));
               switch(hierarchy){
                  case Hierarchies::h5:{
                     const H5 hierarchy5(expansionDepth, Al4p, beta, Dmglst1,
                        lmMt, lmMst1, lmMst2, lmMsq, Mt, Mst1,
                        Mst2, Msq, p.mu,
                        s2t,
                        ho.getMDRFlag(), oneLoopFlag, twoLoopFlag, threeLoopFlag);
                     curSig1 = hierarchy5.getS1();
                     curSig2 = hierarchy5.getS2();
                     curSig12 = hierarchy5.getS12();
                     if(oneLoopFlagIn == 0 && twoLoopFlagIn == 0 && threeLoopFlagIn == 1){
                        const double c = hierarchy5.calc_coef_at_as2_no_sm_logs_log0();
                        ho.setDLambdaH3m(c
                           + lmMst1 * hierarchy5.calc_coef_at_as2_no_sm_logs_log1()
                           + pow2(lmMst1) * hierarchy5.calc_coef_at_as2_no_sm_logs_log2()
                           + pow3(lmMst1) * hierarchy5.calc_coef_at_as2_no_sm_logs_log3());
                        ho.setDLambdaNonLog(c);
                     }
                  }
                  break;
                  case Hierarchies::h5g1:{
                     const H5g1 hierarchy5g1(expansionDepth, Al4p, beta, Dmglst1,
                        lmMt, lmMst1, lmMst2, lmMsq, Mgl, Mt, Mst1,
                        Mst2, Msq, p.mu,
                        s2t,
                        ho.getMDRFlag(), oneLoopFlag, twoLoopFlag, threeLoopFlag);
                     curSig1 = hierarchy5g1.getS1();
                     curSig2 = hierarchy5g1.getS2();
                     curSig12 = hierarchy5g1.getS12();
                     if(oneLoopFlagIn == 0 && twoLoopFlagIn == 0 && threeLoopFlagIn == 1){
                        const double c = hierarchy5g1.calc_coef_at_as2_no_sm_logs_log0();
                        ho.setDLambdaH3m(c
                           + lmMst1 * hierarchy5g1.calc_coef_at_as2_no_sm_logs_log1()
                           + pow2(lmMst1) * hierarchy5g1.calc_coef_at_as2_no_sm_logs_log2()
                           + pow3(lmMst1) * hierarchy5g1.calc_coef_at_as2_no_sm_logs_log3());
                        ho.setDLambdaNonLog(c);
                     }
                  }
                  break;
               }
            }
            break;
            case Hierarchies::h6:{
               const double Dmglst2 = Mgl - Mst2;
               const double lmMst1 = std::log(pow2(p.scale / Mst1));
               const double lmMst2 = std::log(pow2(p.scale / Mst2));
               const double lmMsq = std::log(pow2(p.scale / calcMeanMsq()));
               switch(hierarchy){
                  case Hierarchies::h6:{
                     const H6 hierarchy6(expansionDepth, Al4p, beta, Dmglst2,
                        lmMt, lmMst1, lmMst2, lmMsq,
                        Mt, Mst1, Mst2, Msq, p.mu,
                        s2t,
                        ho.getMDRFlag(), oneLoopFlag, twoLoopFlag, threeLoopFlag);
                     curSig1 = hierarchy6.getS1();
                     curSig2 = hierarchy6.getS2();
                     curSig12 = hierarchy6.getS12();
                     if(oneLoopFlagIn == 0 && twoLoopFlagIn == 0 && threeLoopFlagIn == 1){
                        const double c = hierarchy6.calc_coef_at_as2_no_sm_logs_log0();
                        ho.setDLambdaH3m(c
                           + lmMst1 * hierarchy6.calc_coef_at_as2_no_sm_logs_log1()
                           + pow2(lmMst1) * hierarchy6.calc_coef_at_as2_no_sm_logs_log2()
                           + pow3(lmMst1) * hierarchy6.calc_coef_at_as2_no_sm_logs_log3());
                        ho.setDLambdaNonLog(c);
                     };
                  }
                  break;
                  case Hierarchies::h6g2:{
                     const H6g2 hierarchy6g2(expansionDepth, Al4p, beta, Dmglst2,
                        lmMt, lmMst1, lmMst2, lmMsq,
                        Mgl, Mt, Mst1, Mst2, Msq, p.mu,
                        s2t,
                        ho.getMDRFlag(), oneLoopFlag, twoLoopFlag, threeLoopFlag);
                     curSig1 = hierarchy6g2.getS1();
                     curSig2 = hierarchy6g2.getS2();
                     curSig12 = hierarchy6g2.getS12();
                     if(oneLoopFlagIn == 0 && twoLoopFlagIn == 0 && threeLoopFlagIn == 1){
                        const double c = hierarchy6g2.calc_coef_at_as2_no_sm_logs_log0();
                        ho.setDLambdaH3m(c
                           + lmMst1 * hierarchy6g2.calc_coef_at_as2_no_sm_logs_log1()
                           + pow2(lmMst1) * hierarchy6g2.calc_coef_at_as2_no_sm_logs_log2()
                           + pow3(lmMst1) * hierarchy6g2.calc_coef_at_as2_no_sm_logs_log3());
                        ho.setDLambdaNonLog(c);
                     }
                  }
                  break;
               }
            }
            break;
            case Hierarchies::h6b:{
               const double Dmglst2 = Mgl - Mst2;
               const double Dmsqst2 = Msq - Mst2;
               const double lmMst1 = std::log(pow2(p.scale / Mst1));
               const double lmMst2 = std::log(pow2(p.scale / Mst2));
               switch(hierarchy){
                  case Hierarchies::h6b:{
                     const H6b hierarchy6b(expansionDepth, Al4p, beta, Dmglst2,
                        Dmsqst2, lmMt, lmMst1, lmMst2,
                        Mt, Mst1, Mst2, p.mu,
                        s2t,
                        ho.getMDRFlag(), oneLoopFlag, twoLoopFlag, threeLoopFlag);
                     curSig1 = hierarchy6b.getS1();
                     curSig2 = hierarchy6b.getS2();
                     curSig12 = hierarchy6b.getS12();
                     if(oneLoopFlagIn == 0 && twoLoopFlagIn == 0 && threeLoopFlagIn == 1){
                        const double c = hierarchy6b.calc_coef_at_as2_no_sm_logs_log0();
                        ho.setDLambdaH3m(c
                           + lmMst1 * hierarchy6b.calc_coef_at_as2_no_sm_logs_log1()
                           + pow2(lmMst1) * hierarchy6b.calc_coef_at_as2_no_sm_logs_log2()
                           + pow3(lmMst1) * hierarchy6b.calc_coef_at_as2_no_sm_logs_log3());
                        ho.setDLambdaNonLog(c);
                     }
                  }
                  break;
                  case Hierarchies::h6b2qg2:{
                     const H6b2qg2 hierarchy6b2qg2(expansionDepth, Al4p, beta, Dmglst2,
                        Dmsqst2, lmMt, lmMst1, lmMst2,
                        Mgl, Mt, Mst1, Mst2, p.mu,
                        s2t,
                        ho.getMDRFlag(), oneLoopFlag, twoLoopFlag, threeLoopFlag);
                     curSig1 = hierarchy6b2qg2.getS1();
                     curSig2 = hierarchy6b2qg2.getS2();
                     curSig12 = hierarchy6b2qg2.getS12();
                     if(oneLoopFlagIn == 0 && twoLoopFlagIn == 0 && threeLoopFlagIn == 1){
                        const double c = hierarchy6b2qg2.calc_coef_at_as2_no_sm_logs_log0();
                        ho.setDLambdaH3m(c
                           + lmMst1 * hierarchy6b2qg2.calc_coef_at_as2_no_sm_logs_log1()
                           + pow2(lmMst1) * hierarchy6b2qg2.calc_coef_at_as2_no_sm_logs_log2()
                           + pow3(lmMst1) * hierarchy6b2qg2.calc_coef_at_as2_no_sm_logs_log3());
                        ho.setDLambdaNonLog(c);
                     }
                  }
                  break;
                  case Hierarchies::h6bq22g:{
                     const H6bq22g hierarchy6bq22g(expansionDepth, Al4p, beta, Dmglst2,
                        Dmsqst2, lmMt, lmMst1, lmMst2,
                        Mt, Mst1, Mst2, Msq, p.mu,
                        s2t,
                        ho.getMDRFlag(), oneLoopFlag, twoLoopFlag, threeLoopFlag);
                     curSig1 = hierarchy6bq22g.getS1();
                     curSig2 = hierarchy6bq22g.getS2();
                     curSig12 = hierarchy6bq22g.getS12();
                     if(oneLoopFlagIn == 0 && twoLoopFlagIn == 0 && threeLoopFlagIn == 1){
                        const double c = hierarchy6bq22g.calc_coef_at_as2_no_sm_logs_log0();
                        ho.setDLambdaH3m(c
                           + lmMst1 * hierarchy6bq22g.calc_coef_at_as2_no_sm_logs_log1()
                           + pow2(lmMst1) * hierarchy6bq22g.calc_coef_at_as2_no_sm_logs_log2()
                           + pow3(lmMst1) * hierarchy6bq22g.calc_coef_at_as2_no_sm_logs_log3());
                        ho.setDLambdaNonLog(c);
                     }
                  }
                  break;
                  case Hierarchies::h6bq2g2:{
                     const H6bq2g2 hierarchy6bq2g2(expansionDepth, Al4p, beta, Dmglst2,
                        Dmsqst2, lmMt, lmMst1, lmMst2,
                        Mgl, Mt, Mst1,Mst2, Msq, p.mu,
                        s2t,
                        ho.getMDRFlag(), oneLoopFlag, twoLoopFlag, threeLoopFlag);
                     curSig1 = hierarchy6bq2g2.getS1();
                     curSig2 = hierarchy6bq2g2.getS2();
                     curSig12 = hierarchy6bq2g2.getS12();
                     if(oneLoopFlagIn == 0 && twoLoopFlagIn == 0 && threeLoopFlagIn == 1){
                        const double c = hierarchy6bq2g2.calc_coef_at_as2_no_sm_logs_log0();
                        ho.setDLambdaH3m(c
                           + lmMst1 * hierarchy6bq2g2.calc_coef_at_as2_no_sm_logs_log1()
                           + pow2(lmMst1) * hierarchy6bq2g2.calc_coef_at_as2_no_sm_logs_log2()
                           + pow3(lmMst1) * hierarchy6bq2g2.calc_coef_at_as2_no_sm_logs_log3());
                        ho.setDLambdaNonLog(c);
                     }
                  }
                  break;
               }
            }
            break;
            case Hierarchies::h9:{
               const double lmMst1 = std::log(pow2(p.scale / Mst1));
               const double Dmst12 = pow2(Mst1) - pow2(Mst2);
               const double Dmsqst1 = pow2(Msq) - pow2(Mst1);
               switch(hierarchy){
                  case Hierarchies::h9:{
                     const H9 hierarchy9(expansionDepth, Al4p, beta, Dmst12, Dmsqst1,
                        lmMt, lmMgl, lmMst1,
                        Mgl, Mt, Mst1, Mst2, p.mu,
                        s2t,
                        ho.getMDRFlag(), oneLoopFlag, twoLoopFlag, threeLoopFlag);
                     curSig1 = hierarchy9.getS1();
                     curSig2 = hierarchy9.getS2();
                     curSig12 = hierarchy9.getS12();
                     if(oneLoopFlagIn == 0 && twoLoopFlagIn == 0 && threeLoopFlagIn == 1){
                        const double c = hierarchy9.calc_coef_at_as2_no_sm_logs_log0();
                        ho.setDLambdaH3m(c
                           + lmMst1 * hierarchy9.calc_coef_at_as2_no_sm_logs_log1()
                           + pow2(lmMst1) * hierarchy9.calc_coef_at_as2_no_sm_logs_log2()
                           + pow3(lmMst1) * hierarchy9.calc_coef_at_as2_no_sm_logs_log3());
                        ho.setDLambdaNonLog(c);
                     }
                  }
                  break;
                  case Hierarchies::h9q2:{
                     const H9q2 hierarchy9q2(expansionDepth, Al4p, beta, Dmst12, Dmsqst1,
                        lmMt, lmMgl, lmMst1,
                        Mgl, Mt, Mst1, Mst2, Msq, p.mu,
                        s2t,
                        ho.getMDRFlag(), oneLoopFlag, twoLoopFlag, threeLoopFlag);
                     curSig1 = hierarchy9q2.getS1();
                     curSig2 = hierarchy9q2.getS2();
                     curSig12 = hierarchy9q2.getS12();
                     if(oneLoopFlagIn == 0 && twoLoopFlagIn == 0 && threeLoopFlagIn == 1){
                        const double c = hierarchy9q2.calc_coef_at_as2_no_sm_logs_log0();
                        ho.setDLambdaH3m(c
                           + lmMst1 * hierarchy9q2.calc_coef_at_as2_no_sm_logs_log1()
                           + pow2(lmMst1) * hierarchy9q2.calc_coef_at_as2_no_sm_logs_log2()
                           + pow3(lmMst1) * hierarchy9q2.calc_coef_at_as2_no_sm_logs_log3());
                        ho.setDLambdaNonLog(c);
                     }
                  }
                  break;
               }
            }
            break;
         }
      }
      sigS1Full += curSig1;
      sigS2Full += curSig2;
      sigS12Full += curSig12;
   }

   // add the MDR masses to the hierarchy object only if a 3-loop calculation has to be done, otherwise let the user decide
   if (oneLoopFlagIn == 0 && twoLoopFlagIn == 0 && threeLoopFlagIn == 1) {
      Eigen::Vector2d mdrMasses;
      mdrMasses(0) = Mst1;
      mdrMasses(1) = Mst2;
      ho.setMDRMasses(mdrMasses);
   }

   Eigen::Matrix2d higgsMassMatrix;
   higgsMassMatrix(0, 0) = sigS1Full;
   higgsMassMatrix(0, 1) = sigS12Full;
   higgsMassMatrix(1, 0) = higgsMassMatrix(0, 1);
   higgsMassMatrix(1, 1) = sigS2Full;

   return calcHiggsMassMatrixPrefactor() * higgsMassMatrix;
}

/**
 * Shifts Msx1 according to the hierarchy to the MDR scheme.
 * @param ho a HierarchyObject with constant isAlphab.
 * @param oneLoopFlag an integer flag which is 0 or 1 in order to shift the order O(alpha_s).
 * @param twoLoopFlag an integer flag which is 0 or 1 in order to shift the order O(alpha_s^2).
 * @return A double which is the MDR sx_1 mass.
 */
double HierarchyCalculator::shiftMst1ToMDR(const himalaya::HierarchyObject& ho,
                                           unsigned oneLoopFlag,
                                           unsigned twoLoopFlag) const
{
   using namespace himalaya::hierarchies;

   double Mst1mod = 0., Mst1, Mst2;

   if (!ho.getIsAlphab()) {
      Mst1 = p.MSt(0);
      Mst2 = p.MSt(1);
   } else {
      Mst1 = p.MSb(0);
      Mst2 = p.MSb(1);
   }

   const double Al4p = calcAsOver4Pi();
   const double Msq = calcMeanMsq();
   const double Mgl = p.MG;
   const double lmMst2 = std::log(pow2(p.scale) / pow2(Mst2));
   const double lmMgl = std::log(pow2(p.scale / p.MG));
   const double lmMsq = std::log(pow2(p.scale / calcMeanMsq()));
   const double Dmglst2 = Mgl - Mst2;
   const double mdr2mst1ka = (-8. * twoLoopFlag * pow2(Al4p)
      * (10 * pow2(Msq) * (-1 + 2 * lmMsq + 2 * z2) + pow2(Mst2)
      * (-1 + 2 * lmMst2 + 2 * z2))) / (3. * pow2(Mst1));
   switch (getMotherHierarchy(ho.getSuitableHierarchy())) {
   case Hierarchies::h3:
      Mst1mod = (1 + mdr2mst1ka);
      break;
   case Hierarchies::h4:
      Mst1mod = (1 + mdr2mst1ka);
      break;
   case Hierarchies::h5:
      Mst1mod = (1 + mdr2mst1ka);
      break;
   case Hierarchies::h6:
      Mst1mod = (144 * oneLoopFlag * Al4p * (1 + lmMgl) * pow2(Mgl) * pow4(Msq) + 27 * (1 + mdr2mst1ka) * pow4(Msq) * pow2(Mst1) +
         twoLoopFlag * pow2(Al4p) * Mgl * (-5 * (67 + 84 * lmMgl - 84 * lmMsq) * pow5(Mgl) - 40 * (43 + 30 * lmMgl - 30 * lmMsq) * pow3(Mgl) * pow2(Msq) +
            288 * Dmglst2 * pow4(Msq) * (1 - 2 * z2) + 12 * Mgl * pow4(Msq) * (79 + 144 * pow2(lmMgl) - 150 * lmMsq +
               90 * pow2(lmMsq) - 90 * lmMgl * (-3 + 2 * lmMsq) + 208 * z2))) / (27. * pow4(Msq) * pow2(Mst1));
      break;
   case Hierarchies::h6b:
      Mst1mod = (48 * oneLoopFlag * Al4p * (1 + lmMgl) * pow2(Mgl) + 9 * (1 + mdr2mst1ka) * pow2(Mst1) +
         8 * twoLoopFlag * pow2(Al4p) * (-135 * pow2(Msq) + 12 * Dmglst2 * Mgl * (1 - 22 * z2) +
            pow2(Mgl) * (77 + 135 * lmMgl + 72 * pow2(lmMgl) - 75 * lmMsq -
               90 * lmMgl * lmMsq + 45 * pow2(lmMsq) + 104 * z2))) / (9. * pow2(Mst1));
      break;
   case Hierarchies::h9:
      Mst1mod = (1 + mdr2mst1ka);
      break;
   }
   return Mst1 * std::sqrt(Mst1mod);
}

/**
 * Shifts Msx2 according to the hierarchy to the MDR scheme.
 * @param ho a HierarchyObject with constant isAlphab.
 * @param oneLoopFlag an integer flag which is 0 or 1 in order to shift the order O(alpha_s).
 * @param twoLoopFlag an integer flag which is 0 or 1 in order to shift the order O(alpha_s^2).
 * @return A double which is the MDR sx_2 mass.
 */
double HierarchyCalculator::shiftMst2ToMDR(const himalaya::HierarchyObject& ho,
                                           unsigned oneLoopFlag,
                                           unsigned twoLoopFlag) const
{
   using namespace himalaya::hierarchies;

   double Mst2mod = 0., Mst2;
   if (!ho.getIsAlphab()) {
      Mst2 = p.MSt(1);
   } else {
      Mst2 = p.MSb(1);
   }
   const double Al4p = calcAsOver4Pi();
   const double Mgl = p.MG;
   const double Msq = calcMeanMsq();
   const double lmMgl = std::log(pow2(p.scale / p.MG));
   const double lmMsq = std::log(pow2(p.scale / calcMeanMsq()));
   const double Dmglst2 = Mgl - Mst2;
   const double mdr2mst2ka = (-80. * twoLoopFlag * pow2(Al4p)
      * pow2(Msq) * (-1 + 2 * lmMsq + 2 * z2)) / (3. * pow2(Mst2));
   switch (getMotherHierarchy(ho.getSuitableHierarchy())) {
   case Hierarchies::h3:
      Mst2mod = (1 + mdr2mst2ka);
      break;
   case Hierarchies::h4:
      Mst2mod = (1 + mdr2mst2ka);
      break;
   case Hierarchies::h5:
      Mst2mod = (1 + mdr2mst2ka);
      break;
   case Hierarchies::h6:
      Mst2mod = (144 * oneLoopFlag * Al4p * (1 + lmMgl) * pow2(Mgl) * pow4(Msq) + 27 * (1 + mdr2mst2ka) * pow4(Msq) * pow2(Mst2) +
         twoLoopFlag * pow2(Al4p) * Mgl * (-5 * (67 + 84 * lmMgl - 84 * lmMsq) * pow5(Mgl) - 40 * (43 + 30 * lmMgl - 30 * lmMsq) * pow3(Mgl) * pow2(Msq) +
            288 * Dmglst2 * pow4(Msq) * (1 - 2 * z2) + 12 * Mgl * pow4(Msq) * (79 + 144 * pow2(lmMgl) - 150 * lmMsq +
               90 * pow2(lmMsq) - 90 * lmMgl * (-3 + 2 * lmMsq) + 208 * z2))) / (27. * pow4(Msq) * pow2(Mst2));
      break;
   case Hierarchies::h6b:
      Mst2mod = (48 * oneLoopFlag * Al4p * (1 + lmMgl) * pow2(Mgl) + 9 * (1 + mdr2mst2ka) * pow2(Mst2) +
         8 * twoLoopFlag * pow2(Al4p) * (-135 * pow2(Msq) + 12 * Dmglst2 * Mgl * (1 - 22 * z2) +
            pow2(Mgl) * (77 + 135 * lmMgl + 72 * pow2(lmMgl) - 75 * lmMsq -
               90 * lmMgl * lmMsq + 45 * pow2(lmMsq) + 104 * z2))) / (9. * pow2(Mst2));
      break;
   case Hierarchies::h9:
      Mst2mod = (1 + mdr2mst2ka);
      break;
   }
   return Mst2 * std::sqrt(Mst2mod);
}

/**
 * Shifts the H3m renormalization scheme to DR' scheme. This shift has to be added to the H3m result!
 * @param ho a HierarchyObject with constant isAlphab.
 * @return A matrix which shifts the H3m scheme to the DR' scheme at three-loop level
 *
 */
Eigen::Matrix2d HierarchyCalculator::shiftH3mToDRbarPrime(
   const himalaya::HierarchyObject& ho) const
{
   Eigen::Matrix2d shift;

   // truncate shift at O(Xt^2) to be consistent with H3m result
   int truncateXt = 1;

   const int suitableHierarchy = ho.getSuitableHierarchy();
   if(suitableHierarchy == himalaya::hierarchies::Hierarchies::h3
      || suitableHierarchy == himalaya::hierarchies::Hierarchies::h32q2g
      || suitableHierarchy == himalaya::hierarchies::Hierarchies::h3q22g
      || suitableHierarchy == himalaya::hierarchies::Hierarchies::h9
      || suitableHierarchy == himalaya::hierarchies::Hierarchies::h9q2) truncateXt = 0;

   // pre-factor of shift -> checked normalization against H3m normalization and they coincide
   const double yt = sqrt2 * p.Mt / p.vu;
   const double prefac = threeLoop * pow4(p.g3) * pow2(p.Mt * yt);

   // tanbeta
   const double Tbeta = calcTanBeta();

   // stop masses
   const double Mst1 = shiftMst1ToMDR(ho, ho.getMDRFlag(), ho.getMDRFlag());
   const double Mst2 = shiftMst2ToMDR(ho, ho.getMDRFlag(), ho.getMDRFlag());
   double Xt = p.Au(2,2) - p.mu / Tbeta;
   // Hierarchy h4 only covers O(Xt^0)
   if(suitableHierarchy == himalaya::hierarchies::Hierarchies::h4) Xt = 0;

   // threshold for degenerate squark mass case is 1% of the stop mass
   const double eps = Mst1 * 0.01;

   // squared masses
   const double Mst12 = pow2(Mst1);
   const double Mgl = p.MG;
   const double Mgl2 = pow2(p.MG);
   const double Msq2 = pow2(calcMeanMsq());
   const double scale2 = pow2(p.scale);
   const double Xt2 = pow2(Xt);
   const double Mst22 = pow2(Mst2);
   const double Dmst12 = Mst12 - Mst22;

   // logarithms
   const double lmMst1 = std::log(scale2 / Mst12);
   const double lmMst2 = std::log(scale2 / Mst22);
   const double lmMsq = std::log(scale2 / Msq2);
   const double lmMgl = std::log(scale2 / Mgl2);

   // degenerate mass case flag
   bool isDegen = false;

   // check for degenerate squark masses
   if(std::abs(Mst1 - Mst2) < eps){
      const double Mst2shift = Mst1 + std::sqrt(std::abs(Dmst12))/2.;
      const double lmMst2shift = std::log(scale2 / pow2(Mst2shift));
      // limit
      const double lim = (32*Xt2*(-3*(1 + lmMgl)*pow2(Mgl) + 5*(1 + lmMsq)*Msq2 + (1 +
        lmMst1)*pow2(Mst1))*pow2(p.mu))/(3.*pow6(Mst1));
      // exact result
      const double exact = (16*Xt2*(-6*(1 + lmMgl)*pow2(Mgl) + 10*(1 + lmMsq)*Msq2 + (1 +
        lmMst1)*pow2(Mst1) + (1 + lmMst2)*pow2(Mst2))*pow2(p.mu)*(-4*std::log(Mst1/
        Mst2)*pow2(Mst1)*pow2(Mst2) + pow4(Mst1) - pow4(Mst2)))/(pow2(Mst1)*
        pow2(Mst2)*pow3(pow2(Mst1) - pow2(Mst2)));
      // exact result with shifted stop_2 mass
      const double exactShifted = (16*Xt2*(-6*(1 + lmMgl)*pow2(Mgl) + 10*(1 + lmMsq)*Msq2 + (1 +
        lmMst1)*pow2(Mst1) + (1 + lmMst2shift)*pow2(Mst2shift))*pow2(p.mu)*(-4*std::log(Mst1/
        Mst2shift)*pow2(Mst1)*pow2(Mst2shift) + pow4(Mst1) - pow4(Mst2shift)))/(pow2(Mst1)*
        pow2(Mst2shift)*pow3(pow2(Mst1) - pow2(Mst2shift)));

      isDegen = std::abs(exactShifted - lim) >= std::abs(exact - lim)
         || !std::isfinite(exact) || !std::isfinite(exactShifted);
   }

   if(isDegen){
      // matrix elements
      shift(0, 0) = (32*Xt2*(-3*(1 + lmMgl)*pow2(Mgl) + 5*(1 + lmMsq)*Msq2 + (1 +
        lmMst1)*pow2(Mst1))*pow2(p.mu))/(3.*pow6(Mst1));
      shift(1, 0) = (-32*p.mu*(-3*(1 + lmMgl)*pow2(Mgl) + 5*(1 + lmMsq)*Msq2 + (1 +
        lmMst1)*pow2(Mst1))*(p.mu*Xt2 - 3*Tbeta*Xt*pow2(Mst1) + Tbeta*pow3(Xt)))/
        (3.*Tbeta*pow6(Mst1));
      shift(0, 1) = shift(1,0);
      shift(1, 1) = (32*(-3*(1 + lmMgl)*pow2(Mgl) + 5*(1 + lmMsq)*Msq2 + (1 + lmMst1)*
        pow2(Mst1))*(-6*Tbeta*(p.mu*Xt + Tbeta*Xt2)*pow2(Mst1) + Xt2*pow2(p.mu) +
        truncateXt*pow2(Tbeta)*pow2(Xt2) + 2*p.mu*Tbeta*pow3(Xt) + 6*pow2(Tbeta)*
        pow4(Mst1)))/(3.*pow2(Tbeta)*pow6(Mst1));
   }
   else{
      // matrix elements
      shift(0, 0) = (16*Xt2*(-6*(1 + lmMgl)*pow2(Mgl) + 10*(1 + lmMsq)*Msq2 + (1 +
        lmMst1)*pow2(Mst1) + (1 + lmMst2)*pow2(Mst2))*pow2(p.mu)*(-4*std::log(Mst1/
        Mst2)*pow2(Mst1)*pow2(Mst2) + pow4(Mst1) - pow4(Mst2)))/(pow2(Mst1)*
        pow2(Mst2)*pow3(pow2(Mst1) - pow2(Mst2)));
      shift(1, 0) = (16*p.mu*(-6*(1 + lmMgl)*pow2(Mgl) + 10*(1 + lmMsq)*Msq2 + (1 +
        lmMst1)*pow2(Mst1) + (1 + lmMst2)*pow2(Mst2))*(4*std::log(Mst1/Mst2)*pow2(
        Mst1)*pow2(Mst2)*(p.mu*Xt2 + Tbeta*pow3(Xt)) + p.mu*Xt2*(-pow4(Mst1) +
        pow4(Mst2)) + Tbeta*Xt*(pow3(pow2(Mst1) - pow2(Mst2)) + pow2(Xt)*(-
        pow4(Mst1) + pow4(Mst2)))))/(Tbeta*pow2(Mst1)*pow2(Mst2)*pow3(pow2(
        Mst1) - pow2(Mst2)));
      shift(0, 1) = shift(1,0);
      shift(1, 1) = (16*(-6*(1 + lmMgl)*pow2(Mgl) + 10*(1 + lmMsq)*Msq2 + (1 + lmMst1)*
        pow2(Mst1) + (1 + lmMst2)*pow2(Mst2))*(-4*std::log(Mst1/Mst2)*pow2(Mst1)*
        pow2(Mst2)*(Xt2*pow2(p.mu) + truncateXt*pow2(Tbeta)*pow2(Xt2) + 2*p.mu*
        Tbeta*pow3(Xt)) + Xt2*(-2*pow2(Tbeta)*pow3(pow2(Mst1) - pow2(Mst2)) +
        pow2(p.mu)*(pow4(Mst1) - pow4(Mst2))) + Tbeta*(-2*p.mu*Xt*pow3(pow2(Mst1) -
        pow2(Mst2)) + Tbeta*(pow2(Mst1) + pow2(Mst2))*pow3(pow2(Mst1) - pow2(
        Mst2)) + 2*p.mu*pow3(Xt)*(pow4(Mst1) - pow4(Mst2))) + truncateXt*pow2(
        Tbeta)*pow2(Xt2)*(pow4(Mst1) - pow4(Mst2))))/(pow2(Mst1)*pow2(Mst2)*
        pow2(Tbeta)*pow3(pow2(Mst1) - pow2(Mst2)));
   }

   return prefac * shift;
}

/**
 * Shifts the H3m renormalization scheme to DR' scheme. This shift has to be added to the H3m result!
 * Note: This shift is WITHOUT the three-loop pre-factor g3^4*k^3*Mt^2*yt^2*Sin[beta]^2 with k = 1 / (16 Pi^2)
 * @param ho a HierarchyObject with constant isAlphab.
 * @return A double which shifts the H3m scheme to the DR' scheme at three-loop level
 *
 */
double HierarchyCalculator::shiftH3mToDRbarPrimeMh2(
   const himalaya::HierarchyObject& ho, int omitLogs) const
{
   double shift;

   // truncate shift at O(Xt^2) to be consistent with H3m result
   int truncateXt = 1;
   const int suitableHierarchy = ho.getSuitableHierarchy();
   if(suitableHierarchy == himalaya::hierarchies::Hierarchies::h3
      || suitableHierarchy == himalaya::hierarchies::Hierarchies::h32q2g
      || suitableHierarchy == himalaya::hierarchies::Hierarchies::h3q22g
      || suitableHierarchy == himalaya::hierarchies::Hierarchies::h9
      || suitableHierarchy == himalaya::hierarchies::Hierarchies::h9q2) truncateXt = 0;

   // tanbeta
   const double Tbeta = calcTanBeta();

   // stop masses
   const double Mst1 = p.MSt(0);
   const double Mst2 = p.MSt(1);

   double Xt = p.Au(2,2) - p.mu / Tbeta;
   // Hierarchy h4 only covers O(Xt^0)
   if(suitableHierarchy == himalaya::hierarchies::Hierarchies::h4) Xt = 0;

   // threshold for degenerate squark mass case is 1% of the stop mass
   const double eps = Mst1 * 0.01;

   // squared masses
   const double Mst12 = pow2(Mst1);
   const double Mgl = p.MG;
   const double Mgl2 = pow2(p.MG);
   const double Msq2 = pow2(calcMeanMsq());
   const double scale2 = pow2(p.scale);
   const double Xt2 = pow2(Xt);
   const double Mst22 = pow2(Mst2);
   const double Dmst12 = Mst12 - Mst22;

   // logarithms
   const double lmMst1 = omitLogs * std::log(scale2 / Mst12);
   const double lmMst2 = omitLogs * std::log(scale2 / Mst12) - std::log(Mst22 / Mst12);
   const double lmMsq = omitLogs * std::log(scale2 / Mst12) - std::log(Msq2 / Mst12);
   const double lmMgl = omitLogs * std::log(scale2 / Mst12) - std::log(Mgl2 / Mst12);

   // degenerate mass case flag
   bool isDegen = false;

   // check for degenerate squark masses
   if(std::abs(Mst1 - Mst2) < eps){
      const double Mst2shift = Mst1 + std::sqrt(std::abs(Dmst12))/2.;
      const double lmMst2shift = std::log(scale2 / pow2(Mst2shift));
      // limit
      const double lim = (32*(-3*(1 + lmMgl)*pow2(Mgl) + 5*(1 + lmMsq)*Msq2 + (1 + lmMst1)*
        pow2(Mst1))*(-6*Xt2*pow2(Mst1) + truncateXt*pow2(Xt2) + 6*pow4(Mst1)))/
        (3.*pow6(Mst1));
      // exact result
      const double exact =  (16*(-6*(1 + lmMgl)*pow2(Mgl) + 10*(1 + lmMsq)*Msq2 + (1 + lmMst1)*
        pow2(Mst1) + (1 + lmMst2)*pow2(Mst2))*(4*truncateXt*std::log(Mst2/Mst1)*
        pow2(Mst1)*pow2(Mst2)*pow2(Xt2) - 2*Xt2*pow3(pow2(Mst1) - pow2(Mst2)) +
        (pow2(Mst1) + pow2(Mst2))*pow3(pow2(Mst1) - pow2(Mst2)) + truncateXt*
        pow2(Xt2)*(pow4(Mst1) - pow4(Mst2))))/(pow2(Mst1)*pow2(Mst2)*pow3(pow2(
        Mst1) - pow2(Mst2)));
      // exact result with shifted stop_2 mass
      const double exactShifted =  (16*(-6*(1 + lmMgl)*pow2(Mgl) + 10*(1 + lmMsq)*Msq2 + (1 + lmMst1)*
        pow2(Mst1) + (1 + lmMst2shift)*pow2(Mst2shift))*(4*truncateXt*std::log(Mst2shift/Mst1)*
        pow2(Mst1)*pow2(Mst2shift)*pow2(Xt2) - 2*Xt2*pow3(pow2(Mst1) - pow2(Mst2shift)) +
        (pow2(Mst1) + pow2(Mst2shift))*pow3(pow2(Mst1) - pow2(Mst2shift)) + truncateXt*
        pow2(Xt2)*(pow4(Mst1) - pow4(Mst2shift))))/(pow2(Mst1)*pow2(Mst2shift)*pow3(pow2(
        Mst1) - pow2(Mst2shift)));

      isDegen = std::abs(exactShifted - lim) >= std::abs(exact - lim)
         || !std::isfinite(exact) || !std::isfinite(exactShifted);
   }

   if(isDegen){
      shift = (32*(-3*(1 + lmMgl)*pow2(Mgl) + 5*(1 + lmMsq)*Msq2 + (1 + lmMst1)*
        pow2(Mst1))*(-6*Xt2*pow2(Mst1) + truncateXt*pow2(Xt2) + 6*pow4(Mst1)))/
        (3.*pow6(Mst1));
   }
   else{
      shift = (16*(-6*(1 + lmMgl)*pow2(Mgl) + 10*(1 + lmMsq)*Msq2 + (1 + lmMst1)*
        pow2(Mst1) + (1 + lmMst2)*pow2(Mst2))*(4*truncateXt*std::log(Mst2/Mst1)*
        pow2(Mst1)*pow2(Mst2)*pow2(Xt2) - 2*Xt2*pow3(pow2(Mst1) - pow2(Mst2)) +
        (pow2(Mst1) + pow2(Mst2))*pow3(pow2(Mst1) - pow2(Mst2)) + truncateXt*
        pow2(Xt2)*(pow4(Mst1) - pow4(Mst2))))/(pow2(Mst1)*pow2(Mst2)*pow3(pow2(
        Mst1) - pow2(Mst2)));
   }

   return shift;
}


/**
 * Calculates the loop corrected Higgs mass matrix at the order O(alpha_x). Here, x can be t or b.
 * @param ho a HierarchyObject with constant isAlphab.
 * @param shiftOneLoop An integer flag which is 0 or 1 in order to shift the one-loop terms to the MDR scheme.
 * @param shiftTwoLoop An integer flag which is 0 or 1 in order to shift the two-loop terms to the MDR scheme.
 * @return The loop corrected Higgs mass matrix at the order O(alpha_x).
 */
Eigen::Matrix2d HierarchyCalculator::getMt41L(
   const himalaya::HierarchyObject& ho,
   unsigned shiftOneLoop,
   unsigned shiftTwoLoop) const
{
   using std::log;

   Eigen::Matrix2d Mt41L;
   const double pi2 = Pi*Pi;
   const double GF = 1/(sqrt2 * calcV2());
   const double beta = calcBeta();
   const double Mst1 = shiftMst1ToMDR(ho, shiftOneLoop, shiftTwoLoop);
   const double Mst2 = shiftMst2ToMDR(ho, shiftOneLoop, shiftTwoLoop);
   const double sbeta = std::sin(beta);
   const double cbeta = std::cos(beta);
   double Mt, s2t;
   if (!ho.getIsAlphab()) {
      s2t = p.s2t;
      Mt = p.Mt;
   } else {
      s2t = p.s2b;
      Mt = p.Mb;
   }

   Mt41L(0, 0) = (-3 * GF * pow2(Mt) * pow2(p.mu) * pow2(1 / sbeta) *
      (-pow2(Mst1) + pow2(Mst2) + pow2(Mst1) * log(Mst1) +
      pow2(Mst2) * log(Mst1) - pow2(Mst1) * log(Mst2) -
      pow2(Mst2) * log(Mst2)) * pow2(s2t)) /
      (4. * sqrt2 * (pow2(Mst1) - pow2(Mst2)) * pi2);

   Mt41L(0, 1) = (3 * GF * pow2(1 / sbeta) *
      (-(pow3(Mt) * p.mu * (log(Mst1) - log(Mst2)) * s2t) / 2. +
      (pow2(Mt) * pow2(p.mu) * 1 / tan(beta) *
      (-pow2(Mst1) + pow2(Mst2) + pow2(Mst1) * log(Mst1) +
      pow2(Mst2) * log(Mst1) - pow2(Mst1) * log(Mst2) -
      pow2(Mst2) * log(Mst2)) * pow2(s2t)) /
      (4. * (pow2(Mst1) - pow2(Mst2))) +
      (Mt * p.mu * (-pow2(Mst1) + pow2(Mst2) + pow2(Mst1) * log(Mst1) +
      pow2(Mst2) * log(Mst1) - pow2(Mst1) * log(Mst2) -
      pow2(Mst2) * log(Mst2)) * pow3(s2t)) / 8.)) /
      (sqrt2 * pi2);

   Mt41L (1,0) = Mt41L(0,1);

   Mt41L(1, 1) = (3 * GF * pow2(1 / sbeta) *
      (pow4(Mt) * (log(Mst1) + log(Mst2) - 2 * log(Mt)) +
      pow3(Mt) * p.mu * 1 / tan(beta) * (log(Mst1) - log(Mst2)) * s2t +
      (pow2(Mt) * pow2(1 / sbeta) *
      (pow2(Mst1) * pow2(p.mu) * pow2(cbeta) -
      pow2(Mst2) * pow2(p.mu) * pow2(cbeta) -
      pow2(Mst1) * pow2(p.mu) * pow2(cbeta) * log(Mst1) -
      pow2(Mst2) * pow2(p.mu) * pow2(cbeta) * log(Mst1) +
      pow2(Mst1) * pow2(p.mu) * pow2(cbeta) * log(Mst2) +
      pow2(Mst2) * pow2(p.mu) * pow2(cbeta) * log(Mst2) +
      2 * pow4(Mst1) * log(Mst1) * pow2(sbeta) -
      4 * pow2(Mst1) * pow2(Mst2) * log(Mst1) * pow2(sbeta) +
      2 * pow4(Mst2) * log(Mst1) * pow2(sbeta) -
      2 * pow4(Mst1) * log(Mst2) * pow2(sbeta) +
      4 * pow2(Mst1) * pow2(Mst2) * log(Mst2) * pow2(sbeta) -
      2 * pow4(Mst2) * log(Mst2) * pow2(sbeta)) * pow2(s2t)) /
      (4. * (pow2(Mst1) - pow2(Mst2))) -
      (Mt * p.mu * 1 / tan(beta) * (-pow2(Mst1) + pow2(Mst2) +
      pow2(Mst1) * log(Mst1) + pow2(Mst2) * log(Mst1) -
      pow2(Mst1) * log(Mst2) - pow2(Mst2) * log(Mst2)) *
      pow3(s2t)) / 4. -
      ((pow2(Mst1) - pow2(Mst2)) *
      (-pow2(Mst1) + pow2(Mst2) + pow2(Mst1) * log(Mst1) +
      pow2(Mst2) * log(Mst1) - pow2(Mst1) * log(Mst2) -
      pow2(Mst2) * log(Mst2)) * pow4(s2t)) / 16.)) /
      (sqrt2 * pi2);

    return Mt41L;
}

/**
 * Calculates the loop corrected Higgs mass matrix at the order O(alpha_x*alpha_s). Here, x can be t or b.
 * @param ho a HierarchyObject with constant isAlphab.
 * @param shiftOneLoop An integer flag which is 0 or 1 in order to shift the one-loop terms to the MDR scheme.
 * @param shiftTwoLoop An integer flag which is 0 or 1 in order to shift the two-loop terms to the MDR scheme.
 * @return The loop corrected Higgs mass matrix at the order O(alpha_x*alpha_s).
 */
Eigen::Matrix2d
HierarchyCalculator::getMt42L(const himalaya::HierarchyObject& ho,
                              unsigned shiftOneLoop,
                              unsigned shiftTwoLoop) const
{
   using namespace himalaya::mssm_twoloophiggs;

   double Mt2, st, ct;

   if (!ho.getIsAlphab()) {
      const double theta = p.theta_t;
      Mt2 = pow2(p.Mt);
      st = std::sin(theta);
      ct = std::cos(theta);
   } else {
      const double theta = p.theta_b;
      Mt2 = pow2(p.Mb);
      st = std::sin(theta);
      ct = std::cos(theta);
   }

   const double Mst12 = pow2(shiftMst1ToMDR(ho, shiftOneLoop, shiftTwoLoop));
   const double Mst22 = pow2(shiftMst2ToMDR(ho, shiftOneLoop, shiftTwoLoop));
   const double MG = p.MG;
   const double scale2 = pow2(p.scale);
   const double mu = - p.mu; // note the sign difference in mu
   const double tanb = calcTanBeta();
   const double v2 = calcV2();
   const double gs = p.g3;
   const int include_heavy_higgs = 0;

   const Eigen::Matrix2d Mt42L = delta_mh2_2loop_at_as(
      Mt2, MG, Mst12, Mst22, st, ct, scale2, mu, tanb, v2, gs,
      include_heavy_higgs);

   return Mt42L;
}

/**
 * Calculates the contribution to the order (alpha_x) and (alpha_s alpha_x) as the difference
 * of the Higgs mass matrices of the MDR and DR scheme. Here, x can be t or b.
 * @param ho a HierarchyObject with constant isAlphab.
 * @param shiftOneLoop a bool to shift the terms at one-loop level.
 * @param shiftTwoLoop a bool to shift the terms at two-loop level.
 * @return The loop corrected Higgs mass matrix difference of the MDR and DR scheme at the given order.
 */
Eigen::Matrix2d
HierarchyCalculator::calcDRbarToMDRbarShift(const himalaya::HierarchyObject& ho,
                                            bool shiftOneLoop,
                                            bool shiftTwoLoop) const
{
   if(shiftOneLoop && shiftTwoLoop){
      return getMt41L(ho, 1, 1) + getMt42L(ho, 1, 1) - getMt41L(ho, 0, 0) - getMt42L(ho, 0, 0);
   }
   if(shiftOneLoop){
      return getMt41L(ho, 1, 1) - getMt41L(ho, 0, 0);
   }
   if(shiftTwoLoop){
      return getMt42L(ho, 1, 1) - getMt42L(ho, 0, 0);
   }

   return Eigen::Matrix2d::Zero();
}

double HierarchyCalculator::calcTanBeta() const
{
   return p.vu / p.vd;
}


double HierarchyCalculator::calcBeta() const
{
   return std::atan(calcTanBeta());
}


double HierarchyCalculator::calcV2() const
{
   return pow2(p.vu) + pow2(p.vd);
}


double HierarchyCalculator::calcHiggsMassMatrixPrefactor() const
{
   // GF = 1/(sqrt(2) * (vu^2 + vd^2)) is calculated in the DR'-bar scheme
   return 3. / (2. * calcV2() * Pi * Pi * pow2(std::sin(calcBeta())));
}


double HierarchyCalculator::calcAsOver4Pi() const
{
   return oneLoop * pow2(p.g3);
}


double HierarchyCalculator::calcMeanMsq() const
{
   return std::sqrt(std::abs(p.calculateMsq2()));
}


/**
 * Fills in delta_lambda @ 3L to the given HierarchyObject
 * @param ho a HierrachyObject
 */
void HierarchyCalculator::calcDeltaLambda3L(himalaya::HierarchyObject& ho, bool omitXtOrders) const
{
   // set Xt order truncation for EFT contribution to be consistent with H3m
   const int suitableHierarchy = ho.getSuitableHierarchy();
   const int xtOrder =
      (suitableHierarchy == himalaya::hierarchies::Hierarchies::h3
       || suitableHierarchy == himalaya::hierarchies::Hierarchies::h32q2g
       || suitableHierarchy == himalaya::hierarchies::Hierarchies::h3q22g
       || suitableHierarchy == himalaya::hierarchies::Hierarchies::h9
       || suitableHierarchy == himalaya::hierarchies::Hierarchies::h9q2) ? 3 : 4;

   // to obtain delta_lambda one has to divide the difference of the two calculations by v^2
   const double v2 = calcV2();
   const double gt = sqrt2*p.Mt/std::sqrt(v2);
   const double pref = threeLoop * pow2(p.Mt * gt * pow2(p.g3));

   // create a modified parameters struct and construct
   // Mh2EFTCalculator and ThresholdCalculator
   auto p_mass_ES = p;
   p_mass_ES.mu2(2,2) = pow2(p.MSt(0));
   p_mass_ES.mq2(2,2) = pow2(p.MSt(1));
   himalaya::mh2_eft::Mh2EFTCalculator mh2EFTCalculator(p_mass_ES);
   himalaya::mh2_eft::ThresholdCalculator tc (p_mass_ES);
   disable_non_as_terms(mh2EFTCalculator); // omit all but O(at*as^n)

   // calculate the (non-)logarithmic part of Mh2 without delta_lambda_3L
   //
   // The first line is equivalent to 64 * dytas - 84 * pow2(dytas) -
   // 24 * dytas2 + catas2 including all log(mu^2/Mst1^2) which are
   // also included in the H3m result, and the second line omits these
   // log(mu^2/Mst1^2) terms since they originate from SM
   // contributions. The non-logarithmic part contains only Xt orders
   // up to O(Xt^2) which are also included in H3m so no subtraction
   // is needed. Checked.
   const double subtractionTermH3m = mh2EFTCalculator.getDeltaMh2EFT3Loop(0,1,0);
   const double subtractionTermEFT = mh2EFTCalculator.getDeltaMh2EFT3Loop(0,0,0);

   if(omitXtOrders) tc.setXtOrderOfDeltaLambdaAtAs2(xtOrder);

   // calculate the EFT logs. In the first call we calculate the full reconstructed contribution to delta_lambda_3L
   // including all logarithmic contributions. In the second line we subtract all non-logarithmic contributions
   // to isolate the logarithmic ones. Checked.
   const double eftLogs = pref*(
      tc.getThresholdCorrection(
            mh2_eft::ThresholdCouplingOrders::LAMBDA_AT_AS2, mh2_eft::RenSchemes::DRBARPRIME, 1)
      - tc.getThresholdCorrection(
            mh2_eft::ThresholdCouplingOrders::LAMBDA_AT_AS2, mh2_eft::RenSchemes::DRBARPRIME, 0));

   // calculate the non-logarithmic part of delta_lambda @ 3L
   const double deltaLambda3LNonLog = pref*(ho.getDLambdaNonLog()
      + shiftH3mToDRbarPrimeMh2(ho,0)) - subtractionTermEFT;

   // calculate delta_lambda_H3m
   ho.setDLambdaH3m((pref*(ho.getDLambdaH3m()
      + shiftH3mToDRbarPrimeMh2(ho,1)) - subtractionTermH3m)/v2);

   // caluclate delta_lambda_EFT
   ho.setDLambdaEFT((deltaLambda3LNonLog + eftLogs)/v2);

   // save the non-logarithmic part of delta_lambda @ 3L
   ho.setDLambdaNonLog(deltaLambda3LNonLog/v2);

   // calculate DR' -> MS shift for delta_lambda 3L
   // this shift generates Xt^5*Log(mu) terms for the EFT expression
   ho.setDLambdaDRbarPrimeToMSbarShift(3, pref*tc.getDRbarPrimeToMSbarShift(xtOrder,1,0)/v2);
   // If orders of Xt are omitted, we subtract them from delta_lambda_EFT to be at the same order
   // as delta_lambda_H3m. This ensures that in the hierarchy selection process we don't compare
   // wrong orders of Xt.
   if(omitXtOrders){
      ho.setDLambdaEFT(ho.getDLambdaEFT() + pref*(tc.getDRbarPrimeToMSbarShift(xtOrder,1,0)
         - tc.getDRbarPrimeToMSbarShift(xtOrder,1,1))/v2);
   }

   // set the uncertainty of delta_lambda due to missing Xt terms
   // to summarize: delta_lambda_H3m is always of the x_t order of the suitable
   // hierachy (h3, h9 ~ xt^3, h5, h6, h6b ~ xt^4), whereas delta_lambda_EFT
   // the non-logarithmic term is of order of the hierarchy but the logarithmic
   // part is of order O(xt^4) and using the shift O(xt^5). Note that the shift
   // for delta_lambda_H3m is always of order of the hierarchy as well.
   // the logarithms are omitted since one would double count the error when
   // combining it with delta_exp
   const int xt4Flag = xtOrder == 3 ? 1 : 0;
   ho.setDLambdaXtUncertainty(
      std::abs(xt4Flag*pref/v2*tc.getDRbarPrimeToMSbarXtTerms(tc.getLimit(), 4, 1)));
}


/**
 * Estimates the uncertainty of the expansion at a given order.
 * @param ho a HierarchyObject with constant isAlphab.
 * @param massMatrix the CP-even Higgs mass matrix without the corrections whose uncertainty should be estimated.
 * @param oneLoopFlag an integer flag which is 0 or 1 in order to estimate the uncertainty of the one-loop expansion terms.
 * @param twoLoopFlag an integer flag which is 0 or 1 in order to estimate the uncertainty of the two-loop expansion terms.
 * @param threeLoopFlag an integer flag which is 0 or 1 in order to estimte the uncertainty of the three-loop expansion terms.
 * @return A double which is the estimated uncertainty.
 */
double HierarchyCalculator::getExpansionUncertainty(
   himalaya::HierarchyObject& ho, const Eigen::Matrix2d& massMatrix,
   unsigned oneLoopFlag, unsigned twoLoopFlag, unsigned threeLoopFlag)
{
   using namespace himalaya::hierarchies;

   // re-computes the Higgs mass eigenvalues (modifies its arguments)
   const auto recomputeMh =
      [this, &massMatrix, oneLoopFlag, twoLoopFlag, threeLoopFlag]
      (HierarchyObject& ho) {
         return calcSmallestEigenvalue(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag));
      };

   // reset flags
   expansionDepth.at(ExpansionDepth::Mst) = 1;

   const double Mh = recomputeMh(ho);
   double uncertainty{};

   switch (getMotherHierarchy(ho.getSuitableHierarchy())) {
   case Hierarchies::h3:
      // truncate the expansion at all variables with one order lower than the expansion depth and evaluate the expansion uncertainty
      expansionDepth.at(ExpansionDepth::Dmglst1) = 0;
      uncertainty += pow2(Mh - recomputeMh(ho));
      expansionDepth.at(ExpansionDepth::Dmglst1) = 1;
      expansionDepth.at(ExpansionDepth::Dmsqst1) = 0;
      uncertainty += pow2(Mh - recomputeMh(ho));
      expansionDepth.at(ExpansionDepth::Dmsqst1) = 1;
      expansionDepth.at(ExpansionDepth::Dmst12) = 0;
      uncertainty += pow2(Mh - recomputeMh(ho));
      expansionDepth.at(ExpansionDepth::Dmst12) = 1;
      break;
   case Hierarchies::h4:
      expansionDepth.at(ExpansionDepth::At) = 0;
      uncertainty += pow2(Mh - recomputeMh(ho));
      expansionDepth.at(ExpansionDepth::At) = 1;
      expansionDepth.at(ExpansionDepth::lmMsusy) = 0;
      uncertainty += pow2(Mh - recomputeMh(ho));
      expansionDepth.at(ExpansionDepth::lmMsusy) = 1;
      expansionDepth.at(ExpansionDepth::Msq) = 0;
      uncertainty += pow2(Mh - recomputeMh(ho));
      expansionDepth.at(ExpansionDepth::Msq) = 1;
      expansionDepth.at(ExpansionDepth::Msusy) = 0;
      uncertainty += pow2(Mh - recomputeMh(ho));
      expansionDepth.at(ExpansionDepth::Msusy) = 1;
      break;
   case Hierarchies::h5:
      expansionDepth.at(ExpansionDepth::Dmglst1) = 0;
      uncertainty += pow2(Mh - recomputeMh(ho));
      expansionDepth.at(ExpansionDepth::Dmglst1) = 1;
      expansionDepth.at(ExpansionDepth::Msq) = 0;
      uncertainty += pow2(Mh - recomputeMh(ho));
      expansionDepth.at(ExpansionDepth::Msq) = 1;
      break;
   case Hierarchies::h6:
      expansionDepth.at(ExpansionDepth::Dmglst2) = 0;
      uncertainty += pow2(Mh - recomputeMh(ho));
      expansionDepth.at(ExpansionDepth::Dmglst2) = 1;
      expansionDepth.at(ExpansionDepth::Msq) = 0;
      uncertainty += pow2(Mh - recomputeMh(ho));
      expansionDepth.at(ExpansionDepth::Msq) = 1;
      break;
   case Hierarchies::h6b:
      expansionDepth.at(ExpansionDepth::Dmglst2) = 0;
      uncertainty += pow2(Mh - recomputeMh(ho));
      expansionDepth.at(ExpansionDepth::Dmglst2) = 1;
      expansionDepth.at(ExpansionDepth::Dmsqst2) = 0;
      uncertainty += pow2(Mh - recomputeMh(ho));
      expansionDepth.at(ExpansionDepth::Dmsqst2) = 1;
      break;
   case Hierarchies::h9:
      expansionDepth.at(ExpansionDepth::Dmsqst1) = 0;
      uncertainty += pow2(Mh - recomputeMh(ho));
      expansionDepth.at(ExpansionDepth::Dmsqst1) = 1;
      expansionDepth.at(ExpansionDepth::Dmst12) = 0;
      uncertainty += pow2(Mh - recomputeMh(ho));
      expansionDepth.at(ExpansionDepth::Dmst12) = 1;
      expansionDepth.at(ExpansionDepth::Mgl) = 0;
      uncertainty += pow2(Mh - recomputeMh(ho));
      expansionDepth.at(ExpansionDepth::Mgl) = 1;
      break;
   }

   // set the expansion depth for the next comparison
   expansionDepth.at(ExpansionDepth::Mst) = 0;
   expansionDepth.at(ExpansionDepth::threeLoop) = 0;

   return std::sqrt(uncertainty);
}

/**
 * Maps a hierarchy to it's mother hierarchy.
 * @param hierarchy the key to a hierarchy.
 * @throws runtime_error Throws a runtime_error if the given hierarchy is not included.
 * @returns The key of the mother hierarchy.
 */
int HierarchyCalculator::getMotherHierarchy(int hierarchy) const
{
   using namespace himalaya::hierarchies;

   if (hierarchy < Hierarchies::Hierarchies::FIRST ||
       hierarchy >= Hierarchies::Hierarchies::NUMBER_OF_HIERARCHIES) {
      if (hierarchy == -1) {
         throw std::runtime_error("No suitable hierarchy found!");
      } else {
         throw std::runtime_error("Hierarchy " + std::to_string(hierarchy) +
                                  " not included!");
      }
   }

   // maps all hierarchies to their mother hierarchy
   static const auto hierarchyMap = [] {
      std::array<int, Hierarchies::NUMBER_OF_HIERARCHIES> a{};
      a[hierarchies::Hierarchies::h3]      = hierarchies::Hierarchies::h3;
      a[hierarchies::Hierarchies::h32q2g]  = hierarchies::Hierarchies::h3;
      a[hierarchies::Hierarchies::h3q22g]  = hierarchies::Hierarchies::h3;
      a[hierarchies::Hierarchies::h4]      = hierarchies::Hierarchies::h4;
      a[hierarchies::Hierarchies::h5]      = hierarchies::Hierarchies::h5;
      a[hierarchies::Hierarchies::h5g1]    = hierarchies::Hierarchies::h5;
      a[hierarchies::Hierarchies::h6]      = hierarchies::Hierarchies::h6;
      a[hierarchies::Hierarchies::h6g2]    = hierarchies::Hierarchies::h6;
      a[hierarchies::Hierarchies::h6b]     = hierarchies::Hierarchies::h6b;
      a[hierarchies::Hierarchies::h6b2qg2] = hierarchies::Hierarchies::h6b;
      a[hierarchies::Hierarchies::h6bq22g] = hierarchies::Hierarchies::h6b;
      a[hierarchies::Hierarchies::h6bq2g2] = hierarchies::Hierarchies::h6b;
      a[hierarchies::Hierarchies::h9]      = hierarchies::Hierarchies::h9;
      a[hierarchies::Hierarchies::h9q2]    = hierarchies::Hierarchies::h9;
      return a;
   }();

   return hierarchyMap.at(hierarchy);
}

/**
 * Prints out some information about Himalaya.
 */
void HierarchyCalculator::printInfo() const
{
   std::cerr << "........................................................................\n";
   std::cerr << "Himalaya " << Himalaya_VERSION_MAJOR << "." << Himalaya_VERSION_MINOR << "." << Himalaya_VERSION_RELEASE << "\tѧѦѧ \n";
   std::cerr << "Uses code by\n";
   std::cerr << "  P. Slavich et al. (2-loop αt*αs) [hep-ph/0105096]\n";
   std::cerr << "  P. Slavich et al. (2-loop αt^2) [hep-ph/0305127]\n";
   std::cerr << "  P. Slavich et al. (2-loop αb*ατ) [hep-ph/0406166]\n";
   std::cerr << "  P. Slavich et al. (2-loop ατ^2) [hep-ph/0112177]\n";
   std::cerr << "Uses contributions\n";
   std::cerr << "  3-loop αt*αs^2 to Mh^2 of Kant et al. [arXiv:1005.5709]\n";
   std::cerr << "  2-loop αt*αs to λ of Bagnaschi et.al. [arXiv:1407.4081]\n";
   std::cerr << "  2-loop αt^2 to λ of Bagnaschi et.al. [arXiv:1703.08166]\n";
   std::cerr << "  2-loop αs^2 to yt of Bednyakov et.al. [hep-ph/0210258, hep-ph/0507139]\n";
   std::cerr << "........................................................................\n";
}

} // namespace himalaya
