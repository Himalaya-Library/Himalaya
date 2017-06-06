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
      double getAbsDiff2l() const;							// returns the absolute difference of the Higgs masses at 2-loop level of the exact and the expanded terms
      double getRelDiff2l() const;							// returns the relative difference of the Higgs masses at 2-loop level of the exact and the expanded terms |Mh(2l,exact) - Mh(2l, expanded)|/Mh(2l, exact)
      double getExpUncertainty(const int& loops) const;					// returns the uncertainty of the expansion at the given loop order
      Eigen::Matrix2d getDMh3L() const;							// returns the 3-loop mass matrix for the given hierarchy and parameter point
      Eigen::Matrix2d getDRToMDRShift() const;						// returns the difference of the MDR - DR contributions of the order alpha_x + alpha_x*alpha_s, with x is b,t with respect to isAlphab
      Eigen::Matrix2d getDMh2L() const;							// returns the 2-loop Higgs mass matrix with only alpha_x*alpha_s contributions
      Eigen::Matrix2d getDMh1l() const;							// returns the 1-loop Higgs mass matrix with only alpha_x contributions
      Eigen::Matrix2d getDMh0l() const;							// returns the 0-loop (tree level) Higgs mass matrix
      Eigen::Matrix<double,2,1> getMDRMasses() const;					// returns a vector of the MDR masses. The first entry is Msx1 and the second Msx2, with x is b or t
      void setSuitableHierarchy(const int& hierarchy);					// sets the value of the suitable hierarchy
      void setAbsDiff2l(const double& absDiff2l);					// sets the absolute difference of the Higgs masses at 2-loop level
      void setRelDiff2l(const double& relDiff2l);					// sets the relative difference of the Higgs masses at 2-loop level
      void setExpUncertainty(const int& loops, const double& uncertainty);		// sets the expansion uncertainties
      void setDMh3l(const Eigen::Matrix2d& dMh3l);					// sets the 3-loop Higgs mass matrix
      void setDRToMDRShift(const Eigen::Matrix2d& dMh2l);				// sets the difference of the MDR - DR contributions of the order alpha_x + alpha_x*alpha_s, with x is b or t 
      void setDMh2l(const Eigen::Matrix2d& dMh2l);					// sets the 2-loop Higgs mass matrix with only alpha_x*alpha_s contributions
      void setDMh1l(const Eigen::Matrix2d& dMh1l);					// sets the 1-loop Higgs mass matrix with only alpha_x contributions
      void setDMh0l(const Eigen::Matrix2d& dMh0l);					// sets the 0-loop (tree level) Higgs mass matrix
      void setMDRMasses(Eigen::Matrix< double, int(2), int(1) >& mdrMasses);			// sets the MDR masses
   private:
      bool isAlphab;									// the bool isAlphab
      int hierarchy;									// the suitable hierarchy
      double absDiff2l;									// the absolute difference of the two loop Higgs masses
      double relDiff2l;									// the relative difference of the two loop Higgs masses
      std::map<int, double> expUncertainties;						// the map which holds the expansion uncertainties, the keys are the loop order: 1, 2, 3
      Eigen::Matrix2d dMh3l;								// the 3-loop mass matrix corresponding to the chosen hierarchy
      Eigen::Matrix2d mdrShift;								// the mass matrix of the difference of the MDR - DR contributions of the order alpha_x + alpha_x*alpha_s
      Eigen::Matrix2d dMh2l;								// the two loop mass matrix with only alpha_x*alpha_s contributions
      Eigen::Matrix2d dMh1l;								// the one loop mass matrix with only alpha_x contributions
      Eigen::Matrix2d dMh0l;								// the tree level mass matrix
      Eigen::Matrix<double,2,1> mdrMasses;						// the 'vector' which holds the MDR masses
      Eigen::Matrix<double,2,1> sortVector(Eigen::Matrix< double, int(2), int(1) >& vector);	// sorts the vector
   };
}	// himalaya
#endif	// HierarchyObject_HPP