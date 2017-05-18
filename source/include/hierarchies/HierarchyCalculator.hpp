#ifndef HierarchyCalculator_HPP
#define HierarchyCalculator_HPP

#include <Eigen>
#include <complex>
#include <H3m_interface.hpp>
#include <map>
#include <vector>

namespace h3m{
   class HierarchyCalculator{
   public:
      HierarchyCalculator(const Parameters& p);														// constructor
      int compareHierarchies(const bool isBottom, const double maxError);													// compares deviation of all hierarchies with the exact two-loop result and returns the hierarchy which minimizes the error
      Eigen::Matrix2d calculateHierarchy(const unsigned int tag, const bool isbottom, const unsigned int oneLoopFlagIn, const unsigned int twoLoopFlagIn,
			      const unsigned int threeLoopFlagIn);											// calculates the hierarchy contributions for a specific hierarchy (tag) and a specific loop order
      Eigen::Matrix2d calcDRbarToMDRbarShift(const unsigned int tag, const bool isBottom, const bool shiftOneLoop, const bool shiftTwoLoop);	// calculates the contribution to the order (alpha_t) and (alpha_s alpha_t) in the MDRbar scheme
      double shiftMst1ToMDR(const unsigned int tag, const bool isBottom, const unsigned int twoLoopFlag, const unsigned int threeLoopFlag);		// shifts Mst1 according to the hierarchy to the MDRbar scheme
      double shiftMst2ToMDR(const unsigned int tag, const bool isBottom, const unsigned int twoLoopFlag, const unsigned int threeLoopFlag);		// shifts Mst2 according to the hierarchy to the MDRbar scheme
      double getExpansionError(const unsigned int tag, const bool isBottom, const Eigen::Matrix2d& twoLoopMassMatrix);
      void checkTerms();																// checks the expansion terms
   private:
      Parameters p;
      double Al4p,lmMgl, lmMsq, Mgl, Msq, prefac, z2, z3, z4, deltaDSZ = 0.;
      double B4, D3, DN, OepS2, S2, T1ep;
      void init();
      bool isHierarchySuitable(const unsigned int tag, const bool isBottom);										// checks if the hierarchy is suitable to the spectrum
      std::vector<double> sortEigenvalues(const Eigen::EigenSolver<Eigen::Matrix2d> es);								// sorts the eigenvalus of a 2x2 matrix. The lower index is the lower eigenvalue
      Eigen::Matrix2d getMt41L(const unsigned int tag, const bool isBottom, const unsigned int shiftOneLoop, const unsigned int shiftTwoLoop);		// calculates the one loop higgs mass matrix of the order alpha_t/b
      Eigen::Matrix2d getMt42L(const unsigned int tag, const bool isBottom, const unsigned int shiftOneLoop, const unsigned int shiftTwoLoop);		// calculates the two loop higgs mass matrix of the order alpha_s * alpha_t
      Eigen::Matrix2d getShift(const unsigned int tag, const bool isBottom);										// shifts the 1-loop terms to the MDRbar scheme to compare top level with different hierarchies
      std::map<unsigned int, unsigned int> flagMap;
      std::string tf(const bool tf);															// returns "true" or "false" with respect to the bool tf
      //hierarchy keys
      static const unsigned int h3 		= 0;
      static const unsigned int h32q2g 		= 1;
      static const unsigned int h3q22g 		= 2;
      static const unsigned int h4 		= 3;
      static const unsigned int h5 		= 4;
      static const unsigned int h5g1 		= 5;
      static const unsigned int h6 		= 6;
      static const unsigned int h6b 		= 7;
      static const unsigned int h6b2qg2 	= 8;
      static const unsigned int h6bq22g		= 9;
      static const unsigned int h6bq2g2		= 10;
      static const unsigned int h6g2		= 11;
      static const unsigned int h9 		= 12;
      static const unsigned int h9q2		= 13;
      const std::map<unsigned int, unsigned int> hierarchyMap = {{h3, h3}, {h32q2g, h3}, {h3q22g, h3}, {h4, h4}, {h5, h5}, {h5g1, h5},
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
}	// h3m
#endif	// HierarchyCalculator_HPP