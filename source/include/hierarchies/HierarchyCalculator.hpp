#ifndef HierarchyCalculator_HPP
#define HierarchyCalculator_HPP

#include <Eigen>
#include <Himalaya_interface.hpp>
#include <HierarchyObject.hpp>
#include <map>
#include <vector>

namespace himalaya{
   class HierarchyCalculator{
   public:
      HierarchyCalculator(const Parameters& p);												// constructor
      HierarchyObject calculateDMh3L(const bool& isAlphab);										// calculates the 3-loop mass matrix and other information of the hierarchy selection process
      int compareHierarchies(HierarchyObject& ho);											// compares deviation of all hierarchies with the exact two-loop result and returns the hierarchy which minimizes the error
      Eigen::Matrix2d calculateHierarchy(HierarchyObject& ho, const unsigned int oneLoopFlagIn, const unsigned int twoLoopFlagIn,
			      const unsigned int threeLoopFlagIn);									// calculates the hierarchy contributions for a specific hierarchy (hierarchy) and a specific loop order
      Eigen::Matrix2d calcDRbarToMDRbarShift(const HierarchyObject& ho, const bool shiftOneLoop, const bool shiftTwoLoop);		// calculates the contribution to the order (alpha_t) and (alpha_s alpha_t) in the MDRbar scheme
      Eigen::Matrix2d getMt41L(const HierarchyObject& ho, const unsigned int shiftOneLoop, const unsigned int shiftTwoLoop);		// calculates the one loop higgs mass matrix of the order alpha_t/b
      Eigen::Matrix2d getMt42L(const HierarchyObject& ho, const unsigned int shiftOneLoop, const unsigned int shiftTwoLoop);		// calculates the two loop higgs mass matrix of the order alpha_s * alpha_t
      double shiftMst1ToMDR(const HierarchyObject& ho, const unsigned int oneLoopFlag, const unsigned int twoLoopFlag);			// shifts Mst1 according to the hierarchy to the MDRbar scheme
      double shiftMst2ToMDR(const HierarchyObject& ho, const unsigned int oneLoopFlag, const unsigned int twoLoopFlag);			// shifts Mst2 according to the hierarchy to the MDRbar scheme
      double getExpansionUncertainty(himalaya::HierarchyObject& ho, const Eigen::Matrix2d& massMatrix, const unsigned int oneLoopFlag, const unsigned int twoLoopFlag, const unsigned int threeLoopFlag);
      void checkTerms();														// checks the expansion terms
   private:
      Parameters p;
      double Al4p,lmMgl, lmMsq, Mgl, Msq, prefac, z2, z3, z4, deltaDSZ = 0.;
      double B4, D3, DN, OepS2, S2, T1ep;
      void init();
      bool isHierarchySuitable(const HierarchyObject& ho);										// checks if the hierarchy is suitable to the spectrum
      std::vector<double> sortEigenvalues(const Eigen::EigenSolver<Eigen::Matrix2d> es);						// sorts the eigenvalus of a 2x2 matrix. The lower index is the lower eigenvalue
      Eigen::Matrix2d getShift(const HierarchyObject& ho);										// shifts the 1-loop terms to the MDRbar scheme to compare top level with different hierarchies
      std::map<unsigned int, unsigned int> flagMap;
      std::string tf(const bool tf);													// returns "true" or "false" with respect to the bool tf
      int getCorrectHierarchy(const int hierarchy);											// gets the correct hierarchy according to hierarchy map
      void printInfo();															// prints the information of Himalaya
      //hierarchy keys TODO: use an enum instead?
      static const int h3 		= 0;
      static const int h32q2g 		= 1;
      static const int h3q22g 		= 2;
      static const int h4 		= 3;
      static const int h5 		= 4;
      static const int h5g1 		= 5;
      static const int h6 		= 6;
      static const int h6b 		= 7;
      static const int h6b2qg2 		= 8;
      static const int h6bq22g		= 9;
      static const int h6bq2g2		= 10;
      static const int h6g2		= 11;
      static const int h9 		= 12;
      static const int h9q2		= 13;
      const std::map<int, int> hierarchyMap = {{h3, h3}, {h32q2g, h3}, {h3q22g, h3}, {h4, h4}, {h5, h5}, {h5g1, h5},
	 {h6, h6}, {h6g2, h6}, {h6b, h6b}, {h6b2qg2, h6b}, {h6bq22g, h6b}, {h6bq2g2, h6b}, {h9, h9}, {h9q2, h9}};
      // expansion depth flags
      const unsigned int xx			= 14;
      const unsigned int xxMst			= 15;
      const unsigned int xxDmglst1		= 16;
      const unsigned int xxDmsqst1		= 17;
      const unsigned int xxDmst12		= 18;
      const unsigned int xxAt			= 19;
      const unsigned int xxlmMsusy		= 20;
      const unsigned int xxMsq			= 21;
      const unsigned int xxMsusy		= 22;
      const unsigned int xxDmglst2		= 23;
      const unsigned int xxDmsqst2		= 24;
      const unsigned int xxMgl			= 25;
  };
}	// himalaya
#endif	// HierarchyCalculator_HPP