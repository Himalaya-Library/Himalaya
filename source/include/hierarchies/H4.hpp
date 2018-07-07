// ====================================================================
// This file is part of Himalaya.
//
// Himalaya is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================

#pragma once

#include <map>

namespace himalaya{
namespace hierarchies{

   class H4{
   public:
      /**
       *         Constructor
       *         @param flagMap the flagMap for the truncation of expansion variables
       *         @param Al4p a double alpha_s/4/Pi
       *         @param At a double tri-linear breaking term
       *         @param beta a double which is the mixing angle beta
       *         @param lmMt a double log((<renormalization scale> / Mt)^2)
       *         @param lmMsq a double log((<renormalization scale> / Msq)^2)
       *         @param Mt a double top/bottom quark mass
       *         @param Msusy a double (Mst1 + Mst2 + Mgl) / 3.
       *         @param Msq a double the average squark mass w/o the top squark
       *         @param mdrFlag an int 0 for DR and 1 for MDR scheme
       *         @param oneLoopFlag an int flag to consider the one-loop expansion terms
       *         @param twoLoopFlag an int flag to consider the two-loop expansion terms
       *         @param threeLoopFlag an int flag to consider the three-loop expansion terms
       */
      H4(const std::map<unsigned int, unsigned int>& flagMap, double Al4p, double At, double beta,
                 double lmMt, double lmMsq, double lmMsusy, double Mt, double Msusy, double Msq,
                 int mdrFlag, int oneLoopFlag, int twoLoopFlag, int threeLoopFlag);
      /**
       *         @return The diagonal (1, 1) matrix element of the Higgs mass matrix as a double for the hierarchy 'H4'
       */
      double getS1() const;
      /**
       *         @return The diagonal (2, 2) matrix element of the Higgs mass matrix as a double for the hierarchy 'H4'
       */
      double getS2() const;
      /**
       *         @return The off-diagonal (1, 2) = (2, 1) matrix element of the Higgs mass matrix as a double for the hierarchy 'H4'
       */
      double getS12() const;
      /**
       *         @return returns the susy log^0 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
       */
      double calc_coef_at_as2_no_sm_logs_log0() const;
      /**
       *         @return returns the susy log^1 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
       */
      double calc_coef_at_as2_no_sm_logs_log1() const;
      /**
       *         @return returns the susy log^2 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
       */
      double calc_coef_at_as2_no_sm_logs_log2() const;
      /**
       *         @return returns the susy log^3 term of Mh^2 @ O(at*as^2) without any log(mu^2) terms normalized to DO (H3m*12/Mt^4/Sbeta^2)
       */
      double calc_coef_at_as2_no_sm_logs_log3() const;
   private:
      double At{}, lmMt{}, lmMsq{}, lmMsusy{}, Msusy{}, Mt{}, Msq{}, Cbeta{}, Sbeta{}, Al4p{}; /**< common variables */
      double shiftst1{}, shiftst2{}, shiftst3{}, xAt{}, xMsq{}, xlmMsusy{}, xMsusy{}; /**< MDR and truncation flags */
      int oneLoopFlag{}, twoLoopFlag{}, threeLoopFlag{}; /**< loop flags */
   };

}        // hierarchies
}        // himalaya
