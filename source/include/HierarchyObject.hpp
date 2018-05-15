#pragma once

#include "version.hpp"
#include "Hierarchies.hpp"
#include <iosfwd>
#include <Eigen/Eigenvalues>
#include <vector>
#include <map>

namespace himalaya{
   /**
    * 	The HierarchyObject class.
    */
   class HierarchyObject{
   public:
      /**
       * 	A constructor.
       * 	@param isAlphab the boolean which determines wether the members are proportinal to alpha_b or alpha_t.
       */
      HierarchyObject(bool isAlphab);
      /**
       * 	@return the value of isAlphab.
       */
      bool getIsAlphab() const;
      /**
       * 	@return The key to the suitable hierarchy
       */
      int getSuitableHierarchy() const;
      /**
       * 	@return The absolute difference of the exact and expanded Higgs masses at two-loop level at the order O(alpha_x + alpha_x*alpha_s).
       */
      double getAbsDiff2L() const;
      /**
       * 	@return The relative difference of the exact and expanded Higgs masses at two-loop level at the order O(alpha_x + alpha_x*alpha_s).
       */
      double getRelDiff2L() const;
      /**
       * 	@param loops an integer which can be 1, 2 or 3.
       * 	@return A double which is the expansion uncertainty for the given loop order.
       */
      double getExpUncertainty(int loops) const;
      /**
       * 	@param loops an integer which can be 0, 1, 2, 3. Here, 0 corresponds to the tree-level matrix.
       * 	@return The CP-even Higgs mass matrix at the given loop order.
       */
      Eigen::Matrix2d getDMh(int loops) const;
      /**
       * 	@return The matrix M(MDR) - M(DR) at the order O(alpha_x + alpha_x*alpha_s)
       */
      Eigen::Matrix2d getDRToMDRShift() const;
      /**
       * 	@return A vector of the MDR stop/sbottom masses. The 0th entry corresponds to the lighter particle.
       */
      Eigen::Matrix<double, 2, 1> getMDRMasses() const;
      /**
       * 	@return the MDRFlag integer.
       */
      int getMDRFlag() const;
      /**
       * 	@return Renormalization scheme key
       */
      int getRenormalizationScheme() const;
      /**
       * 	@return 3-loop delta_lambda with Himalaya logs
       */
      double getDeltaLambdaHimalaya() const;
      /**
       * 	@return 3-loop delta_lambda with EFT logs
       */
      double getDeltaLambdaEFT() const;
      /**
       *        @return 3-loop non-logarithmic part of delta_lambda
       */
      double getDeltaLambdaNonLog() const;
      /**
       *        @return uncertainty of 3-loop delta_lambda_Himalaya
       */
      double getDeltaLambdaUncertaintyHimalaya() const;
      /**
       *        @return uncertainty of 3-loop delta_lambda_EFT
       */
      double getDeltaLambdaUncertaintyEFT() const;
      /**
       * 	@return the DR' -> MS shift for delta_lambda_Himalaya which should be added to the DR' result
       */
      double getDRbarPrimeToMSbarShiftHimalaya() const;
      /**
       * 	@return the DR' -> MS shift for delta_lambda_EFT which should be added to the DR' result
       */
      double getDRbarPrimeToMSbarShiftEFT() const;
      /**
       * 	Sets the suitable hierarchy
       * 	@param hierarchy the integer key of the hierarchy.
       */
      void setSuitableHierarchy(int hierarchy);
      /**
       * 	Sets the absolute difference of the Higgs masses at two-loop level
       * 	@param absDiff2L the absolute difference of the Higgs masses as a double.
       */
      void setAbsDiff2L(double absDiff2L);
      /**
       * 	Sets the relative difference ot the Higgs masses at two-loop level
       * 	@param relDiff2L the relative difference of the Higgs masses as a double.
       */
      void setRelDiff2L(double relDiff2L);
      /**
       * 	Sets the uncertainty of the expansion at a given loop level.
       * 	@param loops the integer value of the corresponding loops. Can be 1, 2 or 3.
       * 	@param uncertainty the expansion untertainty at the given loop order as a double.
       */
      void setExpUncertainty(int loops, double uncertainty);
      /**
       * 	Sets the DR -> MDR shift
       * 	@param mdrShift the DR -> MDR shiftet matrix of the form M(MDR) - M(DR).
       */
      void setDRToMDRShift(const Eigen::Matrix2d& mdrShift);
      /**
       * 	Sets the MDR masses
       * 	@param mdrMasses a vector containting the MDR masses with the lightest particle at position 0.
       */
      void setMDRMasses(Eigen::Matrix<double, 2, 1>& mdrMasses);
      /**
       * 	Sets the delta of the CP-even Higgs mass matrix
       * 	@param loops the integer value of the corresponding loops. Can be 0, 1, 2 or 3. 0 corresponds to the tree-level.
       * 	@param dMh the delta of the mass matrix.
       */
      void setDMh(int loops, const Eigen::Matrix2d& dMh);
      /**
       * 	Sets the mdrFlag to calculate the corretions in the without MDR (0) or with MDR (1) shifts
       * 	@param mdrFlag an int. (0) for H3m (DR')- and (1) for MDR-scheme.
       * 	@throws runtime_exception if the flag is neither 0 or 1 an exception is thrown.
       */
      void setMDRFlag(int mdrFlag);
      /**
       * 	Sets the renormalization scheme accodring to the RenScheme enum
       * 	@param renScheme an int according to the RenScheme enum.
       * 	@throws runtime_exception if the flag is not in {0,1,2,3} an exception is thrown
       */
      void setRenormalizationScheme(int renScheme);
      /**
       * 	Sets the delta_lambda at 3-loop order with Himalaya logs.
       * 	@param deltaLambda delta_lambda at 3-loop order.
       */
      void setDeltaLambdaHimalaya(double deltaLambda);
      /**
       * 	Sets the delta_lambda at 3-loop order with EFT logs
       * 	@param deltaLambda delta_lambda at 3-loop order.
       */
      void setDeltaLambdaEFT(double deltaLambda);
      /**
       * 	Sets the non-logarithmic part of delta_lambda at 3-loop order.
       * 	@param deltaLambda constant part of delta_lambda at 3-loop order.
       */
      void setDeltaLambdaNonLog(double deltaLambda);
      /**
       * 	Sets the Xt parts of the uncertainty of delta_lambda_Himalaya
       *        @param uncertainty of 3-loop delta_lambda
       */
      void setDeltaLambdaXtUncertaintyHimalaya(double uncertainty);
      /**
       * 	Sets the Xt parts of the uncertainty of delta_lambda_EFT
       *        @param uncertainty of 3-loop delta_lambda
       */
      void setDeltaLambdaXtUncertaintyEFT(double uncertainty);
      /**
       * 	Sets the DR' -> MS shift for delta_lambda_Himalaya which should be added to the DR' result
       * 	@param shift the DR' -> MS shift which should be added to the 3-loop threshold correction
       */
      void setDRbarPrimeToMSbarShiftHimalaya(double shift);
      /**
       * 	Sets the DR' -> MS shift for delta_lambda_EFT which should be added to the DR' result
       * 	@param shift the DR' -> MS shift which should be added to the 3-loop threshold correction
       */
      void setDRbarPrimeToMSbarShiftEFT(double shift);
      /**
       *	Returns the H3m notation of a given hierarchy.
       *	@param hierarchy An integer of a Himalaya hierarchy.
       *	@return Returns the corresponding H3m notation of the given hierarchy as a string.
       */
      std::string getH3mHierarchyNotation(int hierarchy) const;
   private:
      bool isAlphab{false};								/**< the bool isAlphab */
      int hierarchy{};									/**< the suitable hierarchy */
      int mdrFlag{0};									/**< the MDR-scheme flag */
      int renormalizationScheme{RenSchemes::DRBARPRIME};				/**< the renormalization scheme flag */
      double absDiff2L{};								/**< the absolute difference of the two loop Higgs masses */
      double relDiff2L{};								/**< the relative difference of the two loop Higgs masses */
      std::map<int, double> expUncertainties{};						/**< the map which holds the expansion uncertainties, the keys are the loop order: 1, 2, 3 */
      std::map<int, Eigen::Matrix2d> dMhMap{};						/**< the map which holds all mass matrices at the given loop order */
      Eigen::Matrix2d mdrShift{};							/**< the mass matrix of the difference of the MDR - DR contributions of the order alpha_x + alpha_x*alpha_s */
      Eigen::Matrix<double, 2, 1> mdrMasses{};						/**< the 'vector' which holds the MDR masses */
      double deltaLambdaHimalaya{};							/**< delta_lambda 3-loop with Himalaya logs */
      double deltaLambdaEFT{};								/**< delta_lambda 3-loop with EFT logs */
      double deltaLambdaNonLog{};                                                       /**< delta_lambda 3-loop, non-logarithmic part */
      double drBarPrimeToMSbarShiftHimalaya{};						/**< The shift to convert the DR' delta_lambda_Himalaya to the MS scheme */
      double drBarPrimeToMSbarShiftEFT{};						/**< The shift to convert the DR' delta_lambda_EFT to the MS scheme */
      double deltaLambdaXtUncertaintyHimalaya{};					/**< The uncertainty of delta_lambda_Himalaya due to mising Xt terms */
      double deltaLambdaXtUncertaintyEFT{};						/**< The uncertainty of delta_lambda_EFT due to mising Xt terms */
      /**
       * 	Sorts a vector.
       * 	@param vector The vector which should be sorted.
       *	@return Returns a vector the lightest entry at position 0.
       */
      Eigen::Matrix<double, 2, 1> sortVector(Eigen::Matrix<double, 2, 1>& vector);
   };
   /**
    * 	Prints out all information of the HierarchyObject
    */
   std::ostream& operator<<(std::ostream&, const HierarchyObject&);
}	// himalaya
