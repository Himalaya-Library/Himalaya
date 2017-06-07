#ifndef Himalaya_interface_HPP
#define Himalaya_interface_HPP

#include <complex>
#include <Eigen>

namespace himalaya {

typedef Eigen::Matrix<double,2,1> V2;
typedef Eigen::Matrix<double,2,2> RM22;
typedef Eigen::Matrix<double,3,3> RM33;

/**
 * 	The Himalaya interface struct
 */
struct Parameters {
   // DR-bar parameters
   double scale{};         /** renormalization scale */
   double mu{};            /** mu parameter */
   double g3{};            /** gauge coupling g3 SU(3) */
   double vd{};            /** VEV of down Higgs */
   double vu{};            /** VEV of up Higgs */
   RM33 mq2{RM33::Zero()}; /** soft-breaking squared left-handed squark mass parameters */
   RM33 md2{RM33::Zero()}; /** soft-breaking squared right-handed down-squark mass parameters */
   RM33 mu2{RM33::Zero()}; /** soft-breaking squared right-handed up-squark mass parameters */
   double At{};            /** soft-breaking trilinear coupling */
   double Ab{};            /** soft-breaking trilinear coupling */

   // DR-bar masses
   double MG{};            /** gluino */
   double MW{};            /** W */
   double MZ{};            /** Z */
   double Mt{};            /** top-quark */
   double Mb{};            /** down-quark */
   double MA{};
   V2 MSt{V2::Zero()};     /** stops */
   V2 MSb{V2::Zero()};     /** sbottoms */

   // DR-bar mixing angles
   double s2t{};	   /** sine of 2 times the stop mixing angle */
   double s2b{};	   /** sine of 2 times the sbot mixing angle */
   
   /**
    * 	Checks if the stob/sbottom masses are ordered in the right way. If these masses are wrongly ordered
    * 	the right ordering will be introduced.
    */
   void validate(){
      if (MSt(0) > MSt(1)) {
	 std::swap(MSt(0), MSt(1));
	 s2t *= -1;
      }

      if (MSb(0) > MSb(1)) {
	 std::swap(MSb(0), MSb(1));
	 s2b *= -1;
      }
   };
};

}	//	himalaya

#endif	//	Himalaya_interface_HPP
