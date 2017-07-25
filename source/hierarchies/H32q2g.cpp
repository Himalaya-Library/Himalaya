#define Pi M_PI

#include "H32q2g.hpp"
#include "HierarchyCalculator.hpp"
#include <type_traits>
#include <math.h>

// some templates to perform operations between int's and complex<double>
template< typename T, typename SCALAR > inline
typename std::enable_if< !std::is_same<T,SCALAR>::value, std::complex<T> >::type
operator* ( const std::complex<T>& c, SCALAR n ) { return c * T(n) ; }

template< typename T, typename SCALAR > inline
typename std::enable_if< !std::is_same<T,SCALAR>::value, std::complex<T> >::type
operator* ( SCALAR n, const std::complex<T>& c ) { return T(n) * c ; }

template< typename T, typename SCALAR > inline
typename std::enable_if< !std::is_same<T,SCALAR>::value, std::complex<T> >::type
operator/ ( const std::complex<T>& c, SCALAR n ) { return c / T(n) ; }

template< typename T, typename SCALAR > inline
typename std::enable_if< !std::is_same<T,SCALAR>::value, std::complex<T> >::type
operator/ ( SCALAR n, const std::complex<T>& c ) { return T(n) / c ; }

template< typename T, typename SCALAR > inline
typename std::enable_if< !std::is_same<T,SCALAR>::value, std::complex<T> >::type
operator+ ( const std::complex<T>& c, SCALAR n ) { return c + T(n) ; }

template< typename T, typename SCALAR > inline
typename std::enable_if< !std::is_same<T,SCALAR>::value, std::complex<T> >::type
operator+ ( SCALAR n, const std::complex<T>& c ) { return T(n) + c ; }

template< typename T, typename SCALAR > inline
typename std::enable_if< !std::is_same<T,SCALAR>::value, std::complex<T> >::type
operator- ( const std::complex<T>& c, SCALAR n ) { return c - T(n) ; }

template< typename T, typename SCALAR > inline
typename std::enable_if< !std::is_same<T,SCALAR>::value, std::complex<T> >::type
operator- ( SCALAR n, const std::complex<T>& c ) { return T(n) - c ; }

template <typename T> T pow2(T x)  { return x*x; }
template <typename T> T pow3(T x)  { return x*x*x; }
template <typename T> T pow4(T x)  { return x*x*x*x; }
template <typename T> T pow5(T x)  { return x*x*x*x*x; }
template <typename T> T pow6(T x)  { return x*x*x*x*x*x; }
template <typename T> T pow7(T x)  { return x*x*x*x*x*x*x; }
template <typename T> T pow8(T x)  { return x*x*x*x*x*x*x*x; }
template <typename T> T pow9(T x)  { return x*x*x*x*x*x*x*x*x; }
template <typename T> T power10(T x) { return x*x*x*x*x*x*x*x*x*x; }
template <typename T> T pow11(T x) { return x*x*x*x*x*x*x*x*x*x*x; }
template <typename T> T pow12(T x) { return x*x*x*x*x*x*x*x*x*x*x*x; }

/**
 * 	Constuctor
 * 	@param flagMap the flagMap for the truncation of expansion variables
 * 	@param Al4p a double alpha_s/4/Pi
 * 	@param beta a double which is the mixing angle beta
 * 	@param Dmglst1 a double Mgl - Mst1
 * 	@param Dmst12 a double Mst1^2 - Mst2^2
 * 	@param Dmsqst1 a double Msq^2 - Mst1^2
 * 	@param lmMt a double log((<renormalization scale> / Mt)^2)
 * 	@param lmMst1 a double log((<renormalization scale> / Mst1)^2)
 * 	@param Mt a double top/bottom quark mass
 * 	@param Mst1 a double stop 1 mass
 * 	@param Mst2 a double stop 2 mass
 * 	@param MuSUSY a double mu parameter
 * 	@param s2t a double 2 times the sine of the stop/sbottom quark mixing angle
 * 	@param mdrFlag an int 0 for DR and 1 for MDR scheme
 * 	@param oneLoopFlag an int flag to consider the one-loop expansion terms
 * 	@param twoLoopFlag an int flag to consider the two-loop expansion terms
 * 	@param threeLoopFlag an int flag to consider the three-loop expansion terms
 */
himalaya::H32q2g::H32q2g(std::map<unsigned int, unsigned int> flagMap, double Al4p, double beta,
		 double Dmglst1, double Dmst12, double Dmsqst1, double lmMt, double lmMst1,
		 double Mt, double Mst1, double Mst2, double MuSUSY,
		 double s2t,
		 int mdrFlag, int oneLoopFlag, int twoLoopFlag, int threeLoopFlag){
   // abbrev for tan(beta) and sin(beta0
   double Tbeta = tan(beta);
   double Sbeta = sin(beta);
   // zeta functions
   double z2 = pow2(Pi)/6.;
   double z3 = 1.202056903159594;
   // mdr flags, indicates if one wants to shift the dr stop mass to the mdr stop mass
   int shiftst1 = mdrFlag;
   int shiftst2 = mdrFlag;
   int shiftst3 = mdrFlag;
   // expansion flags
   int xDmglst1 = flagMap.at(HierarchyCalculator::xxDmglst1);
   int xDmst12 = flagMap.at(HierarchyCalculator::xxDmglst1);
   int xDmsqst1 = flagMap.at(HierarchyCalculator::xxDmsqst1);
   s1 = 
   #include "../hierarchies/h32q2g//sigS1Full.inc"
   ;
   s2 = 
   #include "../hierarchies/h32q2g/sigS2Full.inc"
   ;
   s12 = 
   #include "../hierarchies/h32q2g/sigS12Full.inc"
   ;
}

/**
 * 	@return The diagonal (1, 1) matrix element of the Higgs mass matrix as a double for the hierarchy 'H32q2g'
 */
double himalaya::H32q2g::getS1(){
   return s1;
}

/**
 * 	@return The diagonal (2, 2) matrix element of the Higgs mass matrix as a double for the hierarchy 'H32q2g'
 */
double himalaya::H32q2g::getS2(){
   return s2;
}

/**
 * 	@return The off-diagonal (1, 2) = (2, 1) matrix element of the Higgs mass matrix as a double for the hierarchy 'H32q2g'
 */
double himalaya::H32q2g::getS12(){
   return s12;
}


