#ifndef HierarchyObject_HPP
#define HierarchyObject_HPP

#include <Eigen>
#include <vector>
#include <map>

namespace himalaya{
   
   class HierarchyObject{
   public:
      HierarchyObject(const bool& isAlphab);						// constructor
      bool getIsAlphab() const;								// returns if the corresponding object contains contributions proportional to alpha_b or alpha_t
      int getSuitableHierarchy() const;							// returns the suitable hierarchy for the parameter point
      double getAbsDiff2L() const;							// returns the absolute difference of the Higgs masses at 2-loop level of the exact and the expanded terms
      double getRelDiff2L() const;							// returns the relative difference of the Higgs masses at 2-loop level of the exact and the expanded terms |Mh(2l,exact) - Mh(2l, expanded)|/Mh(2l, exact)
      double getExpUncertainty(const int& loops) const;					// returns the uncertainty of the expansion at the given loop order
      Eigen::Matrix2d getDMh(const int& loops) const;					// returns the @loops mass matrix for the given hierarchy and parameter point
      Eigen::Matrix2d getDRToMDRShift() const;						// returns the difference of the MDR - DR contributions of the order alpha_x + alpha_x*alpha_s, with x is b,t with respect to isAlphab
      Eigen::Matrix<double, 2, 1> getMDRMasses() const;					// returns a vector of the MDR masses. The first entry is Msx1 and the second Msx2, with x is b or t
      void setSuitableHierarchy(const int& hierarchy);					// sets the value of the suitable hierarchy
      void setAbsDiff2L(const double& absDiff2L);					// sets the absolute difference of the Higgs masses at 2-loop level
      void setRelDiff2L(const double& relDiff2L);					// sets the relative difference of the Higgs masses at 2-loop level
      void setExpUncertainty(const int& loops, const double& uncertainty);		// sets the expansion uncertainties
      void setDRToMDRShift(const Eigen::Matrix2d& dMh2L);				// sets the difference of the MDR - DR contributions of the order alpha_x + alpha_x*alpha_s, with x is b or t 
      void setMDRMasses(Eigen::Matrix<double, 2, 1>& mdrMasses);			// sets the MDR masses
      void setDMh(const int& loops, const Eigen::Matrix2d& dMh);				// sets the @loops mass matrix for the given hierarchy and parameter point
   private:
      bool isAlphab;									// the bool isAlphab
      int hierarchy;									// the suitable hierarchy
      double absDiff2L;									// the absolute difference of the two loop Higgs masses
      double relDiff2L;									// the relative difference of the two loop Higgs masses
      std::map<int, double> expUncertainties;						// the map which holds the expansion uncertainties, the keys are the loop order: 1, 2, 3
      std::map<int, Eigen::Matrix2d> dMhMap;						// the map which holds all mass matrices at the given loop order
      Eigen::Matrix2d mdrShift;								// the mass matrix of the difference of the MDR - DR contributions of the order alpha_x + alpha_x*alpha_s
      Eigen::Matrix<double, 2, 1> mdrMasses;						// the 'vector' which holds the MDR masses
      Eigen::Matrix<double, 2, 1> sortVector(Eigen::Matrix<double, 2, 1>& vector);	// sorts the vector
   };
}	// himalaya
#endif	// HierarchyObject_HPP