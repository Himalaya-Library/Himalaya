#define Pi M_PI

#include <HierarchyCalculator.hpp>
#include <iostream>
#include <stdexcept>
#include <iomanip>
#include <type_traits>

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

extern "C" void DSZHiggs_(double *t, double *mg, double *T1, double *T2, double *st, double *ct, double *q, double *mu, double *tanb,
      double *v2, double *gs, int *OS, double *S11, double *S22, double *S12);

/*
 * 	constructor
 */
himalaya::HierarchyCalculator::HierarchyCalculator(const Parameters& p){
   printInfo();
   this -> p = p;
   this -> p.validate();
   // init constants
   // imaginary unit
   const std::complex<double> I(0., 1.);

   // Riemann-Zeta
   z2 = pow2(Pi)/6.;
   z3 = 1.202056903159594;
   z4 = pow4(Pi)/90.;

   // polylogs
   double pl412 = 0.51747906167389934317668576113647; // PolyLog[4,1/2]
   std::complex<double> pl2expPi3 (0.27415567780803773941206919444, 1.014941606409653625021202554275); // PolyLog[2, Exp[I Pi / 3]]
   std::complex<double> pl3expPi6sqrt3 (0.51928806536375962552715984277228, - 0.33358157526196370641686908633664); // PolyLog[3, Exp[- I Pi / 6] / Sqrt[3]]

   // polylog functions, checked
   B4 = (-4 * z2 * pow2(log(2)) + 2 / 3.* pow4(log(2)) - 13 / 2. * z4 + 16. * pl412);
   D3 = 6 * z3 - 15 / 4. * z4 - 6. * pow2(std::imag(pl2expPi3));
   DN = 6 * z3 - 4 * z2 * pow2(log(2)) + 2 / 3. * pow4(log(2)) - 21 / 2. * z4 + 16. * pl412;
   OepS2 = - 763 / 32. - (9 * Pi * sqrt(3) * pow2(log(3))) / 16. - (35 * pow3(Pi) * sqrt(3)) / 48.
      + 195 / 16. * z2 - 15 / 4. * z3 + 57 / 16. * z4 + 45 * sqrt(3) / 2. * std::imag(pl2expPi3)
      - 27 * sqrt(3) * std::imag(pl3expPi6sqrt3);
   S2 = 4 * std::imag(pl2expPi3) / (9. * sqrt(3));
   T1ep = - 45 / 2. - (Pi * sqrt(3) * pow2(log(3))) / 8. - (35 * pow3(Pi) * sqrt(3)) / 216. - 9 / 2. * z2 + z3 
      + 6. * sqrt(3) * std::imag(pl2expPi3) - 6. * sqrt(3) * std::imag(pl3expPi6sqrt3);

   // init common variables
   init();
}

/*
 *  	init common variables
 */
void himalaya::HierarchyCalculator::init(){
   // fill flag list
   flagMap.clear();
   for(int i = xx; i <= xxMgl; i++){
      flagMap.insert(std::pair<unsigned int, unsigned int> (i, 1));
   }
   // beta
   const double beta = atan(p.vu / p.vd);

   //sw2
   const double sw2 = 1 - pow2(p.MW / p.MZ);

   // Al4p
   Al4p = pow2(p.g3 / (4 * Pi));

   // MGl
   Mgl = p.MG;

   // Msq, checked
   Msq = (2 * sqrt(p.mq2(0, 0)) + sqrt(p.mu2(0, 0)) + sqrt(p.md2(0, 0))	// sup and sdown
      + 2 * sqrt(p.mq2(1, 1)) + sqrt(p.mu2(1, 1)) + sqrt(p.md2(1, 1))	// scharm and sstrange
      // sbottom
      + sqrt(p.mq2(2, 2) + pow2(p.Mb) - (1 / 2. - 1 / 3. * sw2) * pow2(p.MZ) * cos(2 * beta))
      + sqrt(p.md2(2, 2) + pow2(p.Mb) - 1 / 3. * sw2 * pow2(p.MZ) * cos(2 * beta))) / 10.;

   // lmMsq, checked
   lmMsq = log(pow2(p.scale / Msq));

   // lmMgl, checked
   lmMgl = log(pow2(p.scale / Mgl));

   // prefactor, GF = 1/(sqrt(2) * (vu^2 + vd^2)) (here, GF is calculated in the DRbar scheme, checked)
   prefac = (3. / (sqrt(2) * (pow2(p.vu) + pow2(p.vd)) * sqrt(2) * pow2(Pi) * pow2(sin(beta))));
}

/*
 * 	calculates everything that is needed for the 3-loop Higgs mass matrix and returns these results in the hierarchy object ho
 */
himalaya::HierarchyObject himalaya::HierarchyCalculator::calculateDMh3L(bool isAlphab){
   HierarchyObject ho (isAlphab);
   
   // compare hierarchies and get the best fitting hierarchy
   compareHierarchies(ho);
   
   // calculate the DR to MDR shift with the obtained hierarchy
   ho.setDRToMDRShift(calcDRbarToMDRbarShift(ho, true, true));
   
   // calculate the 3-loop Higgs mass matrix for the obtained hierarhy
   ho.setDMh(3, calculateHierarchy(ho, 0, 0, 1));
   
   // set the alpha_x contributions
   ho.setDMh(1, getMt41L(ho, 1, 1));
   
   // set the alpha_x*alpha_s contributions
   ho.setDMh(2, getMt42L(ho, 1, 1));
   
   // estimate the uncertainty of the expansion at 3-loop level
   ho.setExpUncertainty(3, getExpansionUncertainty(ho,
						   ho.getDMh(0) + ho.getDMh(1) + ho.getDMh(2), 0, 0, 1));
   
   // set the uncertainty of the expansion at 1-loop level to 0 by convention, if the user needs this value getExpansionUncertainty should be called
   ho.setExpUncertainty(1, 0.);
   
   return ho;
}



/*
 * 	compares deviation of all hierarchies with the exact 2-loop result and returns the hierarchy which minimizes the error
 */
int himalaya::HierarchyCalculator::compareHierarchies(himalaya::HierarchyObject& ho){
   // set flags to truncate the expansion
   flagMap.at(xx) = 0;
   flagMap.at(xxMst) = 0;
   double error = -1.;
   int suitableHierarchy = -1;
   // sine of 2 times beta
   const double s2b = sin(2*atan(p.vu/p.vd));
   const double tbeta = p.vu/p.vd;
   // tree level Higgs mass matrix
   Eigen::Matrix2d treelvl;
   treelvl (0,0) = s2b/2.*(pow2(p.MZ) / tbeta + pow2(p.MA) * tbeta);
   treelvl (1,0) = s2b/2.*(-pow2(p.MZ) - pow2(p.MA));
   treelvl (0,1) = treelvl (1,0);
   treelvl (1,1) = s2b/2.*(pow2(p.MZ) * tbeta + pow2(p.MA) / tbeta);
   
   ho.setDMh(0, treelvl);

   // compare the exact higgs mass at 2-loop level with the expanded expressions to find a suitable hierarchy
   for(int hierarchy = h3; hierarchy <= h9q2; hierarchy ++){
      // first, check if the hierarchy is suitable to the mass spectrum
      ho.setSuitableHierarchy(hierarchy);
      if(isHierarchySuitable(ho)){
	 // calculate the exact 1-loop result (only alpha_t/b)
	 Eigen::Matrix2d Mt41L = getMt41L(ho, 1, 0);
	 
	 // call the routine of Pietro Slavich to get the alpha_s alpha_t/b corrections with the MDRbar masses
	 Eigen::Matrix2d Mt42L = getMt42L(ho, 1, 0);
	 
	 // check for spurious poles. If this is the case slightly change Mst2
	 if(std::isnan(Mt42L(0,0)) || std::isnan(Mt42L(1,0)) || std::isnan(Mt42L(1,1))){
	    deltaDSZ = 1.0E-6;
	    Mt42L = getMt42L(ho, 1, 0);
	 }
	 
	 //DEPRECATED calc 1-loop shift for DRbar -> MDRbar
	 //calc difference of Mt41L or Mt41L in the MDRbar scheme directly
	 //it seems that in H3m the sign of the function getShift is wrong as well(Mt4LDRbar - Mt4LMDRbar)????
	 //to be consistent everything should be calculated in the MDRbar-scheme so we should subtract Mt4LDRbar, shouldn't we?
	 //Eigen::Matrix2d shift = getShift(hierarchyMap.at(hierarchy), isAlphab);
	 
	 //calculate the exact higgs mass at 2-loop (only up to alpha_s alpha_t/b)
	 Eigen::EigenSolver<Eigen::Matrix2d> es2L (treelvl + Mt41L + Mt42L);
	 double Mh2l = sortEigenvalues(es2L).at(0);

	 // calculate the expanded 2-loop expression with the specific hierarchy
	 Eigen::EigenSolver<Eigen::Matrix2d> esExpanded (treelvl + Mt41L + calculateHierarchy(ho, 0, 1, 0));
	 
	 // calculate the higgs mass in the given mass hierarchy and compare the result to estimate the error
	 double Mh2LExpanded = sortEigenvalues(esExpanded).at(0);

	 // estimate the error
	 double twoLoopError = fabs((Mh2l - Mh2LExpanded));

	 // estimate the uncertainty of the expansion
	 double expUncertainty = getExpansionUncertainty(ho, treelvl + Mt41L, 0, 1, 0);
	 
	 // add these errors to include the error of the expansion in the comparison
	 double currError = sqrt(pow2(twoLoopError) + pow2(expUncertainty));
	 
	 // if the error is negative, it is the first iteration and there is no hierarchy which fits better
	 if(error < 0){
	    error = currError;
	    suitableHierarchy = hierarchy;
	    ho.setAbsDiff2L(twoLoopError);
	    ho.setRelDiff2L(twoLoopError/Mh2l);
	    ho.setExpUncertainty(2, expUncertainty);
	 }
	 // compare the current error with the last error and choose the hierarchy which fits best (lowest error)
	 else if(currError < error){
	    error = currError;
	    suitableHierarchy = hierarchy;
	    ho.setAbsDiff2L(twoLoopError);
	    ho.setRelDiff2L(twoLoopError/Mh2l);
	    ho.setExpUncertainty(2, expUncertainty);
	 }
      }
      else{
	 ho.setSuitableHierarchy(suitableHierarchy);	// captures the case if no hierarchy is fitting
      }
   }
   ho.setSuitableHierarchy(suitableHierarchy);
   // reset the flags
   flagMap.at(xx) = 1;
   flagMap.at(xxMst) = 1;
   return suitableHierarchy;
}

/*
 * 	calculates the expanded self-energy to a given order and for a specific hierarchy
 */
Eigen::Matrix2d himalaya::HierarchyCalculator::calculateHierarchy(himalaya::HierarchyObject& ho, const unsigned int oneLoopFlagIn, const unsigned int twoLoopFlagIn, const unsigned int threeLoopFlagIn) {
   // get the hierarchy
   const int hierarchy = ho.getSuitableHierarchy();

   // the hierarchy files containing 1-, 2- and 3-loop terms (alpha_s^0 alpha_t/b, alpha_s alpha_t/b, alpha_s^2 alpha_t/b)
   double sigS1Full = 0., sigS2Full = 0., sigS12Full = 0.;

   // common variables
   double At, Mt, s2t, Mst1, Mst2;
   if (!ho.getIsAlphab()) {
      At = p.At;
      Mt = p.Mt;
      s2t = p.s2t;
   }
   else {
      At = p.Ab;
      Mt = p.Mb;
      s2t = p.s2b;
   }

   const double Tbeta = p.vu / p.vd;
   const double Cbeta = cos(atan(Tbeta));
   const double Sbeta = sin(atan(Tbeta));
   const double scale = p.scale;
   const double lmMt = log(pow2(scale / Mt));
   const double MuSUSY = p.mu;

   // these shifts are needed to eliminate huge corrections (cf. arXiv:1005.5709 [hep-ph] eq. (46) - (49))
   int shiftst1 = 1, shiftst2 = 1, shiftst3 = 1;
   // specific variables for hierarchies
   double Dmglst1, Dmglst2, Dmsqst1, Dmsqst2, Dmst12, lmMst1, lmMst2, Msusy, lmMsusy;
   int xDR2DRMOD;
   
   // flags to truncate the expansion while comparing at 2-loop level or to estimate the error
   int x, xMst, xDmglst1, xDmsqst1, xDmst12, xAt, xlmMsusy, xMsq, xMsusy, xDmglst2, xDmsqst2, xMgl;
   x = flagMap.at(xx);
   xMst = flagMap.at(xxMst);
   xDmglst1 = flagMap.at(xxDmglst1);
   xDmsqst1 = flagMap.at(xxDmsqst1);
   xDmst12 = flagMap.at(xxDmst12);
   xAt = flagMap.at(xxAt);
   xlmMsusy = flagMap.at(xxlmMsusy);
   xMsq = flagMap.at(xxMsq);
   xMsusy = flagMap.at(xxMsusy);
   xDmglst2 = flagMap.at(xxDmglst2);
   xDmsqst2 = flagMap.at(xxDmsqst2);
   xMgl = flagMap.at(xxMgl);
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
	 // select the suitable hierarchy for the specific hierarchy and set variables
	 switch(getCorrectHierarchy(hierarchy)){
	    case h3:
	       if(threeLoopFlag == 0){
		  Mst1 = shiftMst1ToMDR(ho, twoLoopFlag, 0);
		  Mst2 = shiftMst2ToMDR(ho, twoLoopFlag, 0);
	       }
	       else{
		  Mst1 = shiftMst1ToMDR(ho, 1, 1);
		  Mst2 = shiftMst2ToMDR(ho, 1, 1);
	       }
	       Dmglst1 = Mgl - Mst1;
	       Dmsqst1 = pow2(Msq) - pow2(Mst1);
	       Dmst12 = pow2(Mst1) - pow2(Mst2);
	       lmMst1 = log(pow2(scale / Mst1));
	       lmMsusy = log(pow2(scale / ((Mst1 + Mst2 + Mgl + 10*Msq) / 13.)));
	       switch(hierarchy){
		  case h3:
		     curSig1 = 
		     #include "../hierarchies/h3/sigS1Full.inc"
		     ;
		     curSig2 = 
		     #include "../hierarchies/h3/sigS2Full.inc"
		     ;
		     curSig12 = 
		     #include "../hierarchies/h3/sigS12Full.inc"
		     ;
		     break;
		  case h32q2g:
		     curSig1 = 
		     #include "../hierarchies/h32q2g/sigS1Full.inc"
		     ;
		     curSig2 =
		     #include "../hierarchies/h32q2g/sigS2Full.inc"
		     ;
		     curSig12 = 
		     #include "../hierarchies/h32q2g/sigS12Full.inc"
		     ;
		     break;
		  case h3q22g:
		     curSig1 = 
		     #include "../hierarchies/h3q22g/sigS1Full.inc"
		     ;
		     curSig2 = 
		     #include "../hierarchies/h3q22g/sigS2Full.inc"
		     ;
		     curSig12 = 
		     #include "../hierarchies/h3q22g/sigS12Full.inc"
		     ;
		     break;
	       }
	       break;
	    case h4:
	       if(threeLoopFlag == 0){
		  Mst1 = shiftMst1ToMDR(ho, twoLoopFlag, 0);
		  Mst2 = shiftMst2ToMDR(ho, twoLoopFlag, 0);
	       }
	       else{
		  Mst1 = shiftMst1ToMDR(ho, 1, 1);
		  Mst2 = shiftMst2ToMDR(ho, 1, 1);
	       }
	       Msusy = (Mst1 + Mst2 + Mgl) / 3.;
	       lmMsusy = log(pow2(scale / Msusy));
	       curSig1 = 
	       #include "../hierarchies/h4/sigS1Full.inc"
	       ;
	       curSig2 = 
	       #include "../hierarchies/h4/sigS2Full.inc"
	       ;
	       curSig12 = 
	       #include "../hierarchies/h4/sigS12Full.inc"
	       ;
	       break;
	    case h5:
	       if(threeLoopFlag == 0){
		  Mst1 = shiftMst1ToMDR(ho, twoLoopFlag, 0);
		  Mst2 = shiftMst2ToMDR(ho, twoLoopFlag, 0);
	       }
	       else{
		  Mst1 = shiftMst1ToMDR(ho, 1, 1);
		  Mst2 = shiftMst2ToMDR(ho, 1, 1);
	       }
	       Dmglst1 = Mgl - Mst1;
	       lmMst1 = log(pow2(scale / Mst1));
	       lmMst2 = log(pow2(scale / Mst2));
	       switch(hierarchy){
		  case h5:
		     curSig1 = 
		     #include "../hierarchies/h5/sigS1Full.inc"
		     ;
		     curSig2 = 
		     #include "../hierarchies/h5/sigS2Full.inc"
		     ;
		     curSig12 = 
		     #include "../hierarchies/h5/sigS12Full.inc"
		     ;
		     break;
		  case h5g1:
		     curSig1 = 
		     #include "../hierarchies/h5g1/sigS1Full.inc"
		     ;
		     curSig2 = 
		     #include "../hierarchies/h5g1/sigS2Full.inc"
		     ;
		     curSig12 = 
		     #include "../hierarchies/h5g1/sigS12Full.inc"
		     ;
		     break;
	       }
	       break;
	    case h6:
	       if(threeLoopFlag == 0){
		  Mst1 = shiftMst1ToMDR(ho, twoLoopFlag, 0);
		  Mst2 = shiftMst2ToMDR(ho, twoLoopFlag, 0);
	       }
	       else{
		  Mst1 = shiftMst1ToMDR(ho, 1, 1);
		  Mst2 = shiftMst2ToMDR(ho, 1, 1);
	       }
	       Dmglst2 = Mgl - Mst2;
	       lmMst1 = log(pow2(scale / Mst1));
	       lmMst2 = log(pow2(scale / Mst2));
	       xDR2DRMOD = 1;
	       switch(hierarchy){
		  case h6:
		     curSig1 = 
		     #include "../hierarchies/h6/sigS1Full.inc"
		     ;
		     curSig2 = 
		     #include "../hierarchies/h6/sigS2Full.inc"
		     ;
		     curSig12 = 
		     #include "../hierarchies/h6/sigS12Full.inc"
		     ;
		     break;
		  case h6g2:
		     curSig1 = 
		     #include "../hierarchies/h6g2/sigS1Full.inc"
		     ;
		     curSig2 = 
		     #include "../hierarchies/h6g2/sigS2Full.inc"
		     ;
		     curSig12 = 
		     #include "../hierarchies/h6g2/sigS12Full.inc"
		     ;
		     break;
	       }
	       break;
	    case h6b:
	       if(threeLoopFlag == 0){
		  Mst1 = shiftMst1ToMDR(ho, twoLoopFlag, 0);
		  Mst2 = shiftMst2ToMDR(ho, twoLoopFlag, 0);
	       }
	       else{
		  Mst1 = shiftMst1ToMDR(ho, 1, 1);
		  Mst2 = shiftMst2ToMDR(ho, 1, 1);
	       }
	       Dmglst2 = Mgl - Mst2;
	       Dmsqst2 = Msq - Mst2;
	       lmMst1 = log(pow2(scale / Mst1));
	       lmMst2 = log(pow2(scale / Mst2));
	       xDR2DRMOD = 1;
	       switch(hierarchy){
		  case h6b:
		     curSig1 = 
		     #include "../hierarchies/h6b/sigS1Full.inc"
		     ;
		     curSig2 = 
		     #include "../hierarchies/h6b/sigS2Full.inc"
		     ;
		     curSig12 = 
		     #include "../hierarchies/h6b/sigS12Full.inc"
		     ;
		     break;
		  case h6b2qg2:
		     curSig1 = 
		     #include "../hierarchies/h6b2qg2/sigS1Full.inc"
		     ;
		     curSig2 = 
		     #include "../hierarchies/h6b2qg2/sigS2Full.inc"
		     ;
		     curSig12 = 
		     #include "../hierarchies/h6b2qg2/sigS12Full.inc"
		     ;
		     break;
		  case h6bq22g:
		     curSig1 = 
		     #include "../hierarchies/h6bq22g/sigS1Full.inc"
		     ;
		     curSig2 = 
		     #include "../hierarchies/h6bq22g/sigS2Full.inc"
		     ;
		     curSig12 = 
		     #include "../hierarchies/h6bq22g/sigS12Full.inc"
		     ;
		     break;
		  case h6bq2g2:
		     curSig1 = 
		     #include "../hierarchies/h6bq2g2/sigS1Full.inc"
		     ;
		     curSig2 = 
		     #include "../hierarchies/h6bq2g2/sigS2Full.inc"
		     ;
		     curSig12 = 
		     #include "../hierarchies/h6bq2g2/sigS12Full.inc"
		     ;
		     break;
	       }
	       break;
	    case h9:
	       if(threeLoopFlag == 0){
		  Mst1 = shiftMst1ToMDR(ho, twoLoopFlag, 0);
		  Mst2 = shiftMst2ToMDR(ho, twoLoopFlag, 0);
	       }
	       else{
		  Mst1 = shiftMst1ToMDR(ho, 1, 1);
		  Mst2 = shiftMst2ToMDR(ho, 1, 1);
	       }
	       lmMst1 = log(pow2(scale / Mst1));
	       Dmst12 = pow2(Mst1) - pow2(Mst2);
	       Dmsqst1 = pow2(Msq) - pow2(Mst1);
	       switch(hierarchy){
		  case h9:
		     curSig1 = 
		     #include "../hierarchies/h9/sigS1Full.inc"
		     ;
		     curSig2 = 
		     #include "../hierarchies/h9/sigS2Full.inc"
		     ;
		     curSig12 = 
		     #include "../hierarchies/h9/sigS12Full.inc"
		     ;
		     break;
		  case h9q2:
		     curSig1 = 
		     #include "../hierarchies/h9q2/sigS1Full.inc"
		     ;
		     curSig2 = 
		     #include "../hierarchies/h9q2/sigS2Full.inc"
		     ;
		     curSig12 = 
		     #include "../hierarchies/h9q2/sigS12Full.inc"
		     ;
		     break;
	       }
	       break;
	 }
      }
      sigS1Full += curSig1;
      sigS2Full += curSig2;
      sigS12Full += curSig12;
   }
   // add the MDR masses to the hierarchy object only if a 3-loop calculation has to be done, otherwise let the user decide
   if(oneLoopFlagIn == 0 && twoLoopFlagIn == 0 && threeLoopFlagIn == 1){
      Eigen::Matrix<double,2,1> mdrMasses;
      mdrMasses(0) = Mst1;
      mdrMasses(1) = Mst2;
      ho.setMDRMasses(mdrMasses);
   }
   Eigen::Matrix2d higgsMassMatrix;
   higgsMassMatrix(0, 0) = prefac * sigS1Full;
   higgsMassMatrix(0, 1) = prefac * sigS12Full;
   higgsMassMatrix(1, 0) = higgsMassMatrix(0, 1);
   higgsMassMatrix(1, 1) = prefac * sigS2Full;
   return higgsMassMatrix;
}

/*
 * 	checks if hierarchy is suitable to the given mass spectrum
 */
bool himalaya::HierarchyCalculator::isHierarchySuitable(const himalaya::HierarchyObject& ho){
   double Mst1, Mst2;
   if(!ho.getIsAlphab()){
      Mst1 = p.MSt(0);
      Mst2 = p.MSt(1);
   }
   else{
      Mst1 = p.MSb(0);
      Mst2 = p.MSb(1);
   }
   switch (ho.getSuitableHierarchy()){
      case h3:
	 return Mgl > Mst2;
      case h32q2g:
	 return (Mst2 >= Msq) && (Mst2 > Mgl);
      case h3q22g:
	 return (Msq > Mst2) && (Mst2 > Mgl);
      case h4:
	 return (Mst1 < Msq) && (Mst1 >= Mgl);
      case h5:
	 return (Mst2 - Mst1 > 0.1*Mst1) && ((Mgl - Mst1) < std::abs(Mgl - Mst2)) && (Mst2 < Msq) && (Mst1 >= Mgl);
      case h5g1:
	 return (Mst2 - Mst1 > 0.1*Mst1) && ((Mgl - Mst1) < std::abs(Mgl - Mst2)) && (Mst2 < Msq) && (Mgl > Mst1);
      case h6:
	 return (Mst2 - Mst1 > 0.1*Mst1) && ((Mst2 - Mgl) < std::abs(Mgl - Mst1)) && (Mst2 < Msq) && (Mst2 >= Mgl);
      case h6g2:
	 return (Mst2 - Mst1 > 0.1*Mst1) && ((Mst2 - Mgl) < std::abs(Mgl - Mst1)) && (Mst2 < Msq) && (Mgl > Mst2);
      case h6b:
	 return (Mst2 - Mst1 > 0.1*Mst1) && ((Mst2 - Mgl) < std::abs(Mgl - Mst1)) && (Mst2 >= Msq) && (Mst2 >= Mgl);
      case h6b2qg2:
	 return (Mst2 - Mst1 > 0.1*Mst1) && ((Mst2 - Mgl) < std::abs(Mgl - Mst1)) && (Mst2 >= Msq) && (Mgl > Mst2);
      case h6bq22g:
	 return (Mst2 - Mst1 > 0.1*Mst1) && ((Mst2 - Mgl) < std::abs(Mgl - Mst1)) && (Msq > Mst2) && (Mst2 >= Mgl);
      case h6bq2g2:
	 return (Mst2 - Mst1 > 0.1*Mst1) && ((Mst2 - Mgl) < std::abs(Mgl - Mst1)) && (Msq > Mst2) && (Mgl > Mst2);
      case h9:
	 return (Mst2 >= Msq) && ((Mst2 - Mst1) < (Mst1 - Mgl));
      case h9q2:
	 return (Msq > Mst2) && ((Mst1 - Mst1) < (Mst1 - Mgl));
   }
}

/*
 * 	shifts Mst1/Msb1 according to the hierarchy to the MDRbar scheme, checked
 */
double himalaya::HierarchyCalculator::shiftMst1ToMDR(const himalaya::HierarchyObject& ho, const unsigned int oneLoopFlag, const unsigned int twoLoopFlag) {
   double Mst1mod = 0., Mst1, Mst2;
   if(!ho.getIsAlphab()){
      Mst1 = p.MSt(0);
      Mst2 = p.MSt(1);
   }
   else{
      Mst1 = p.MSb(0);
      Mst2 = p.MSb(1);
   }
   double lmMst2 = log(pow2(p.scale) / pow2(Mst2));
   double Dmglst2 = Mgl - Mst2;
   double mdr2mst1ka = (-8. * twoLoopFlag * pow2(Al4p) * (10 * pow2(Msq) * (-1 + 2 * lmMsq + 2 * z2) + pow2(Mst2) * (-1 + 2 * lmMst2 + 2 * z2))) / (3. * pow2(Mst1));
   switch (getCorrectHierarchy(ho.getSuitableHierarchy())) {
   case h3:
      Mst1mod = (1 + mdr2mst1ka);
      break;
   case h4:
      Mst1mod = (1 + mdr2mst1ka);
      break;
   case h5:
      Mst1mod = (1 + mdr2mst1ka);
      break;
   case h6:
      Mst1mod = (144 * oneLoopFlag * Al4p * (1 + lmMgl) * pow2(Mgl) * pow4(Msq) + 27 * (1 + mdr2mst1ka) * pow4(Msq) * pow2(Mst1) +
         twoLoopFlag * pow2(Al4p) * Mgl * (-5 * (67 + 84 * lmMgl - 84 * lmMsq) * pow5(Mgl) - 40 * (43 + 30 * lmMgl - 30 * lmMsq) * pow3(Mgl) * pow2(Msq) +
            288 * Dmglst2 * pow4(Msq) * (1 - 2 * z2) + 12 * Mgl * pow4(Msq) * (79 + 144 * pow2(lmMgl) - 150 * lmMsq +
               90 * pow2(lmMsq) - 90 * lmMgl * (-3 + 2 * lmMsq) + 208 * z2))) / (27. * pow4(Msq) * pow2(Mst1));
      break;
   case h6b:
      Mst1mod = (48 * oneLoopFlag * Al4p * (1 + lmMgl) * pow2(Mgl) + 9 * (1 + mdr2mst1ka) * pow2(Mst1) +
         8 * twoLoopFlag * pow2(Al4p) * (-135 * pow2(Msq) + 12 * Dmglst2 * Mgl * (1 - 22 * z2) +
            pow2(Mgl) * (77 + 135 * lmMgl + 72 * pow2(lmMgl) - 75 * lmMsq -
               90 * lmMgl * lmMsq + 45 * pow2(lmMsq) + 104 * z2))) / (9. * pow2(Mst1));
      break;
   case h9:
      Mst1mod = (1 + mdr2mst1ka);
      break;
   }
   return Mst1 * sqrt(Mst1mod);
}

/*
 * 	shifts Mst2/Msb2 according to the hierarchy to the MDRbar scheme, checked
 */
double himalaya::HierarchyCalculator::shiftMst2ToMDR(const himalaya::HierarchyObject& ho, const unsigned int oneLoopFlag, const unsigned int twoLoopFlag) {
   double Mst2mod;
   double Mst1;
   double Mst2;
   if(!ho.getIsAlphab()){
      Mst1 = p.MSt(0);
      Mst2 = p.MSt(1);
   }
   else{
      Mst1 = p.MSb(0);
      Mst2 = p.MSb(1);
   }
   double Dmglst2 = Mgl - Mst2;
   double mdr2mst2ka = (-80. * twoLoopFlag * pow2(Al4p) * pow2(Msq) * (-1 + 2 * lmMsq + 2 * z2)) / (3. * pow2(Mst2));
   switch (getCorrectHierarchy(ho.getSuitableHierarchy())) {
   case h3:
      Mst2mod = (1 + mdr2mst2ka);
      break;
   case h4:
      Mst2mod = (1 + mdr2mst2ka);
      break;
   case h5:
      Mst2mod = (1 + mdr2mst2ka);
      break;
   case h6:
      Mst2mod = (144 * oneLoopFlag * Al4p * (1 + lmMgl) * pow2(Mgl) * pow4(Msq) + 27 * (1 + mdr2mst2ka) * pow4(Msq) * pow2(Mst2) +
         twoLoopFlag * pow2(Al4p) * Mgl * (-5 * (67 + 84 * lmMgl - 84 * lmMsq) * pow5(Mgl) - 40 * (43 + 30 * lmMgl - 30 * lmMsq) * pow3(Mgl) * pow2(Msq) +
            288 * Dmglst2 * pow4(Msq) * (1 - 2 * z2) + 12 * Mgl * pow4(Msq) * (79 + 144 * pow2(lmMgl) - 150 * lmMsq +
               90 * pow2(lmMsq) - 90 * lmMgl * (-3 + 2 * lmMsq) + 208 * z2))) / (27. * pow4(Msq) * pow2(Mst2));
      break;
   case h6b:
      Mst2mod = (48 * oneLoopFlag * Al4p * (1 + lmMgl) * pow2(Mgl) + 9 * (1 + mdr2mst2ka) * pow2(Mst2) +
         8 * twoLoopFlag * pow2(Al4p) * (-135 * pow2(Msq) + 12 * Dmglst2 * Mgl * (1 - 22 * z2) +
            pow2(Mgl) * (77 + 135 * lmMgl + 72 * pow2(lmMgl) - 75 * lmMsq -
               90 * lmMgl * lmMsq + 45 * pow2(lmMsq) + 104 * z2))) / (9. * pow2(Mst2));
      break;
   case h9:
      Mst2mod = (1 + mdr2mst2ka);
      break;
   }
   return Mst2 * sqrt(Mst2mod);
}


/*
 * 	sorts the eigenvalues of a 2x2 matrix and returns them in a vector
 */
std::vector<double> himalaya::HierarchyCalculator::sortEigenvalues(const Eigen::EigenSolver<Eigen::Matrix2d> es){
  std::vector<double> sortedEigenvalues = {sqrt(std::real(es.eigenvalues()(0))), sqrt(std::real(es.eigenvalues()(1)))};
  std::sort(sortedEigenvalues.begin(), sortedEigenvalues.end());
  return sortedEigenvalues;
}

/*
 * 	calculates the 1-loop alpha_t/b Higgs mass matrix
 */
Eigen::Matrix2d himalaya::HierarchyCalculator::getMt41L(const himalaya::HierarchyObject& ho, const unsigned int shiftOneLoop, const unsigned int shiftTwoLoop){
   Eigen::Matrix2d Mt41L;
   double GF = 1/(sqrt(2) * (pow2(p.vu) + pow2(p.vd)));
   double Mst1;
   double Mst2;
   double Mt;
   double s2t;
   const double beta = atan(p.vu/p.vd);
   Mst1 = shiftMst1ToMDR(ho, shiftOneLoop, shiftTwoLoop);
   Mst2 = shiftMst2ToMDR(ho, shiftOneLoop, shiftTwoLoop);
   if(!ho.getIsAlphab()){
      s2t = p.s2t;
      Mt = p.Mt;
   }
   else{
      s2t = p.s2b;
      Mt = p.Mb;
   }
   Mt41L (0,0) = (-3*GF*pow2(Mt)*pow2(p.mu)*pow2(1/sin(beta))*
      (-pow2(Mst1) + pow2(Mst2) + pow2(Mst1)*log(Mst1) + 
        pow2(Mst2)*log(Mst1) - pow2(Mst1)*log(Mst2) - 
        pow2(Mst2)*log(Mst2))*pow2(s2t))/
    (4.*sqrt(2)*(pow2(Mst1) - pow2(Mst2))*pow2(Pi));
   Mt41L (0,1) = (3*GF*pow2(1/sin(beta))*
      (-(pow3(Mt)*p.mu*(log(Mst1) - log(Mst2))*s2t)/2. + 
        (pow2(Mt)*pow2(p.mu)*1/tan(beta)*
           (-pow2(Mst1) + pow2(Mst2) + pow2(Mst1)*log(Mst1) + 
             pow2(Mst2)*log(Mst1) - pow2(Mst1)*log(Mst2) - 
             pow2(Mst2)*log(Mst2))*pow2(s2t))/
         (4.*(pow2(Mst1) - pow2(Mst2))) + 
        (Mt*p.mu*(-pow2(Mst1) + pow2(Mst2) + pow2(Mst1)*log(Mst1) + 
             pow2(Mst2)*log(Mst1) - pow2(Mst1)*log(Mst2) - 
             pow2(Mst2)*log(Mst2))*pow3(s2t))/8.))/
    (sqrt(2)*pow2(Pi));
   Mt41L (1,0) = Mt41L(0,1);
   Mt41L (1,1) =  (3*GF*pow2(1/sin(beta))*
      (pow4(Mt)*(log(Mst1) + log(Mst2) - 2*log(Mt)) + 
        pow3(Mt)*p.mu*1/tan(beta)*(log(Mst1) - log(Mst2))*s2t + 
        (pow2(Mt)*pow2(1/sin(beta))*
           (pow2(Mst1)*pow2(p.mu)*pow2(cos(beta)) - 
             pow2(Mst2)*pow2(p.mu)*pow2(cos(beta)) - 
             pow2(Mst1)*pow2(p.mu)*pow2(cos(beta))*log(Mst1) - 
             pow2(Mst2)*pow2(p.mu)*pow2(cos(beta))*log(Mst1) + 
             pow2(Mst1)*pow2(p.mu)*pow2(cos(beta))*log(Mst2) + 
             pow2(Mst2)*pow2(p.mu)*pow2(cos(beta))*log(Mst2) + 
             2*pow4(Mst1)*log(Mst1)*pow2(sin(beta)) - 
             4*pow2(Mst1)*pow2(Mst2)*log(Mst1)*pow2(sin(beta)) + 
             2*pow4(Mst2)*log(Mst1)*pow2(sin(beta)) - 
             2*pow4(Mst1)*log(Mst2)*pow2(sin(beta)) + 
             4*pow2(Mst1)*pow2(Mst2)*log(Mst2)*pow2(sin(beta)) - 
             2*pow4(Mst2)*log(Mst2)*pow2(sin(beta)))*pow2(s2t))/
         (4.*(pow2(Mst1) - pow2(Mst2))) - 
        (Mt*p.mu*1/tan(beta)*(-pow2(Mst1) + pow2(Mst2) + 
             pow2(Mst1)*log(Mst1) + pow2(Mst2)*log(Mst1) - 
             pow2(Mst1)*log(Mst2) - pow2(Mst2)*log(Mst2))*
           pow3(s2t))/4. - 
        ((pow2(Mst1) - pow2(Mst2))*
           (-pow2(Mst1) + pow2(Mst2) + pow2(Mst1)*log(Mst1) + 
             pow2(Mst2)*log(Mst1) - pow2(Mst1)*log(Mst2) - 
             pow2(Mst2)*log(Mst2))*pow4(s2t))/16.))/
    (sqrt(2)*pow2(Pi));
    return Mt41L;
}

/*
 * 	DEPRECATED! Just use differences of Mt41L
 * 	calculates the DRbar -> MDRbar shift for the one loop contribution
 */
Eigen::Matrix2d himalaya::HierarchyCalculator::getShift(const himalaya::HierarchyObject& ho){
   Eigen::Matrix2d shift;
   double GF = 1/(sqrt(2) * (pow2(p.vu) + pow2(p.vd)));
   double Mst1;
   double Mst2;
   double Mt;
   double s2t;
   const double beta = atan(p.vu/p.vd);
   double deltamst1, deltamst2;
   if(!ho.getIsAlphab()){
      Mst1 = p.MSt(0);
      Mst2 = p.MSt(1);
      s2t = p.s2t;
      Mt = p.Mt;
      deltamst1 = -(Mst1 - shiftMst1ToMDR(ho, 1, 0));
      deltamst2 = -(Mst2 - shiftMst2ToMDR(ho, 1, 0));
   }
   else{
      Mst1 = p.MSb(0);
      Mst2 = p.MSb(1);
      s2t = p.s2b;
      Mt = p.Mb;
      deltamst1 = -(Mst1 - shiftMst1ToMDR(ho, 1, 0));
      deltamst2 = -(Mst2 - shiftMst2ToMDR(ho, 1, 0));
   }
   shift(0, 0) = (3 * GF * (deltamst2 * Mst1 - deltamst1 * Mst2) * pow2(Mt) * pow2(p.mu) * pow2(1 / Pi) *
      pow2(1/sin(beta)) * pow2(1 / (pow2(Mst1) - pow2(Mst2))) * pow2(s2t) *
      (4 * (log(Mst1) - log(Mst2)) * pow2(Mst1) * pow2(Mst2) - pow4(Mst1) +
	 pow4(Mst2))) / (4. * sqrt(2) * Mst1 * Mst2);
   shift(0, 1) = (-3 * GF * Mt * p.mu * pow2(1 / Pi) * pow2(1/sin(beta)) *
      pow2(1 / (pow2(Mst1) - pow2(Mst2))) * s2t *
      (-(pow2(pow2(Mst1) - pow2(Mst2)) *
	    (4 * (-(deltamst2 * Mst1) + deltamst1 * Mst2) * pow2(Mt) +
	       (-2 * Mst1 * Mst2 * (deltamst1 * Mst1 + deltamst2 * Mst2) *
		  (log(Mst1) - log(Mst2)) +
		  (deltamst2 * Mst1 + deltamst1 * Mst2) * pow2(Mst1) -
		  (deltamst2 * Mst1 + deltamst1 * Mst2) * pow2(Mst2)) * pow2(s2t))) + 2 * (deltamst2 * Mst1 - deltamst1 * Mst2) * Mt * p.mu * 1/tan(beta) *
	 (4 * (log(Mst1) - log(Mst2)) * pow2(Mst1) * pow2(Mst2) - pow4(Mst1) +
	    pow4(Mst2)) * s2t)) / (8. * sqrt(2) * Mst1 * Mst2);
   shift(1, 0) = shift(0, 1);
   shift(1, 1) = (3 * GF * pow2(1 / Pi) * pow2(1/sin(beta)) *
      ((Mt * p.mu * 1/tan(beta) * (-(deltamst1 * Mst1) + deltamst2 * Mst2 +
	       2 * deltamst1 * Mst1 * log(Mst1) + 2 * deltamst2 * Mst2 * log(Mst1) -
	       2 * deltamst1 * Mst1 * log(Mst2) - 2 * deltamst2 * Mst2 * log(Mst2) -
	       (deltamst2 * pow2(Mst1)) / Mst2 + (deltamst1 * pow2(Mst2)) / Mst1) *
	    pow3(s2t)) / 4. + (pow2(Mt) * pow2(1 / (pow2(Mst1) - pow2(Mst2))) *
	    pow2(s2t) * (2 * deltamst2 * pow7(Mst1) -
	       2 * deltamst1 * pow6(Mst1) * Mst2 -
	       2 * deltamst2 * Mst1 * pow6(Mst2) + 2 * deltamst1 * pow7(Mst2) -
	       4 * deltamst1 * pow6(Mst1) * Mst2 * log(Mst1) +
	       4 * deltamst2 * Mst1 * pow6(Mst2) * log(Mst1) +
	       4 * deltamst1 * pow6(Mst1) * Mst2 * log(Mst2) -
	       4 * deltamst2 * Mst1 * pow6(Mst2) * log(Mst2) -
	       deltamst2 * pow5(Mst1) * pow2(p.mu) * pow2(1/tan(beta)) -
	       deltamst1 * pow5(Mst2) * pow2(p.mu) * pow2(1/tan(beta)) +
	       2 * deltamst2 * pow2(Mst2) *
	       (pow5(Mst1) * (-3 + 2 * log(Mst1) - 2 * log(Mst2)) +
		  2 * (log(Mst1) - log(Mst2)) * pow2(p.mu) * pow2(1/tan(beta)) *
		  pow3(Mst1)) - 2 * deltamst1 * pow2(Mst1) *
	       (pow5(Mst2) * (3 + 2 * log(Mst1) - 2 * log(Mst2)) +
		  2 * (log(Mst1) - log(Mst2)) * pow2(p.mu) * pow2(1/tan(beta)) *
		  pow3(Mst2)) + deltamst1 * Mst2 * pow2(p.mu) * pow2(1/tan(beta)) *
	       pow4(Mst1) + 6 * deltamst1 * pow3(Mst2) * pow4(Mst1) +
	       8 * deltamst1 * log(Mst1) * pow3(Mst2) * pow4(Mst1) -
	       8 * deltamst1 * log(Mst2) * pow3(Mst2) * pow4(Mst1) +
	       deltamst2 * Mst1 * pow2(p.mu) * pow2(1/tan(beta)) * pow4(Mst2) +
	       6 * deltamst2 * pow3(Mst1) * pow4(Mst2) -
	       8 * deltamst2 * log(Mst1) * pow3(Mst1) * pow4(Mst2) +
	       8 * deltamst2 * log(Mst2) * pow3(Mst1) * pow4(Mst2))) / (4. * Mst1 * Mst2) -
	 ((deltamst2 * Mst1 + deltamst1 * Mst2) * pow4(Mt)) / (Mst1 * Mst2) -
	 ((deltamst2 * pow5(Mst1) + deltamst1 * pow5(Mst2) -
	    4 * deltamst2 * pow2(Mst2) * pow3(Mst1) -
	    4 * deltamst1 * pow2(Mst1) * pow3(Mst2) +
	    3 * deltamst1 * Mst2 * pow4(Mst1) -
	    4 * deltamst1 * Mst2 * log(Mst1) * pow4(Mst1) +
	    4 * deltamst1 * Mst2 * log(Mst2) * pow4(Mst1) +
	    3 * deltamst2 * Mst1 * pow4(Mst2) +
	    4 * deltamst2 * Mst1 * log(Mst1) * pow4(Mst2) -
	    4 * deltamst2 * Mst1 * log(Mst2) * pow4(Mst2)) * pow4(s2t)) /
	 (16. * Mst1 * Mst2) + (-(deltamst1 / Mst1) + deltamst2 / Mst2) * p.mu *
	 1/tan(beta) * pow3(Mt) * s2t)) / sqrt(2);
   return shift;
}


/*
 * 	calculates the 2-loop higgs mass matrix according to Pietro Slavich
 */
Eigen::Matrix2d himalaya::HierarchyCalculator::getMt42L(const himalaya::HierarchyObject& ho, const unsigned int shiftOneLoop, const unsigned int shiftTwoLoop){
   const int hierarchy = getCorrectHierarchy(ho.getSuitableHierarchy());
   Eigen::Matrix2d Mt42L;
   double S11, S12, S22;
   double Mt2;
   double MG = p.MG;
   double Mst12;
   double Mst22;
   double st;
   double ct;
   Mst12 = pow2(shiftMst1ToMDR(ho, shiftOneLoop, shiftTwoLoop));
   Mst22 = pow2(shiftMst2ToMDR(ho, shiftOneLoop, shiftTwoLoop) + deltaDSZ);
   if(!ho.getIsAlphab()){
      const double theta = asin(p.s2t)/2.;
      Mt2 = pow2(p.Mt);
      st = sin(theta);
      ct = cos(theta);
   }
   else{
      const double theta = asin(p.s2b)/2.;
      Mt2 = pow2(p.Mb);
      st = sin(theta);
      ct = cos(theta);
   }
   double scale2 = pow2(p.scale);
   // note the sign difference in mu
   double mu = - p.mu;
   double tanb = p.vu/p.vd;
   double v2 = pow2(p.vu) + pow2(p.vd);
   double gs = p.g3;
   int os = 0;
   DSZHiggs_(&Mt2, &MG, &Mst12, &Mst22, &st, &ct, &scale2, &mu, &tanb, &v2, &gs, &os, &S11, &S22, &S12);
   Mt42L(0, 0) = S11;
   Mt42L(1, 0) = S12;
   Mt42L(0, 1) = S12;
   Mt42L(1, 1) = S22;
   return Mt42L;
}


/*
 * 	calculates the contribution to the order (alpha_t) and (alpha_s alpha_t) in the MDRbar scheme
 * 	the result is the difference in the Higgs mass matrix of the MDRbar and DRbar scheme
 */
Eigen::Matrix2d himalaya::HierarchyCalculator::calcDRbarToMDRbarShift(const himalaya::HierarchyObject& ho, const bool shiftOneLoop, const bool shiftTwoLoop){
   if(shiftOneLoop && shiftTwoLoop){
      return getMt41L(ho, 1, 1) + getMt42L(ho, 1, 1) - getMt41L(ho, 0, 0) - getMt42L(ho, 0, 0);
   }
   else if(shiftOneLoop){
      return getMt41L(ho, 1, 1) - getMt41L(ho, 0, 0);
   }
   else if(shiftTwoLoop){
      return getMt42L(ho, 1, 1) - getMt42L(ho, 0, 0);
   }
}


/*
 * 	evaluates the error due to the expansion in mass ratios and differences
 * 	first calc the full 3-loop contribution and than truncate the expansion at different variable orders to estimate the expansion error
 */
double himalaya::HierarchyCalculator::getExpansionUncertainty(himalaya::HierarchyObject& ho, const Eigen::Matrix2d& massMatrix, const unsigned int oneLoopFlag, const unsigned int twoLoopFlag, const unsigned int threeLoopFlag){
   double Mh;
   double Mhcut;
   std::vector<double> errors;
   // reset flags
   flagMap.at(xxMst) = 1;
   Eigen::EigenSolver<Eigen::Matrix2d> es;
   switch (getCorrectHierarchy(ho.getSuitableHierarchy())) {
   case h3:
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mh = sortEigenvalues(es).at(0);
      // truncate the expansion at all variables with one order lower than the expansion depth and evaluate the expansion error 
      flagMap.at(xxDmglst1) = 0;
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mhcut = sortEigenvalues(es).at(0);
      errors.push_back(fabs(Mh - Mhcut));
      flagMap.at(xxDmglst1) = 1;
      flagMap.at(xxDmsqst1) = 0;
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mhcut = sortEigenvalues(es).at(0);
      errors.push_back(fabs(Mh - Mhcut));
      flagMap.at(xxDmsqst1) = 1;
      flagMap.at(xxDmst12) = 0;
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mhcut = sortEigenvalues(es).at(0);
      errors.push_back(fabs(Mh - Mhcut));
      flagMap.at(xxDmst12) = 1;
      break;
   case h4:
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mh = sortEigenvalues(es).at(0);
      flagMap.at(xxAt) = 0;
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mhcut = sortEigenvalues(es).at(0);
      errors.push_back(fabs(Mh - Mhcut));
      flagMap.at(xxAt) = 1;
      flagMap.at(xxlmMsusy) = 0;
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mhcut = sortEigenvalues(es).at(0);
      errors.push_back(fabs(Mh - Mhcut));
      flagMap.at(xxlmMsusy) = 1;
      flagMap.at(xxMsq) = 0;
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mhcut = sortEigenvalues(es).at(0);
      errors.push_back(fabs(Mh - Mhcut));
      flagMap.at(xxMsq) = 1;
      flagMap.at(xxMsusy) = 0;
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mhcut = sortEigenvalues(es).at(0);
      errors.push_back(fabs(Mh - Mhcut));
      flagMap.at(xxMsusy) = 1;
      break;
   case h5:
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mh = sortEigenvalues(es).at(0);
      flagMap.at(xxDmglst1) = 0;
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mhcut = sortEigenvalues(es).at(0);
      errors.push_back(fabs(Mh - Mhcut));
      flagMap.at(xxDmglst1) = 1;
      flagMap.at(xxMsq) = 0;
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mhcut = sortEigenvalues(es).at(0);
      errors.push_back(fabs(Mh - Mhcut));
      flagMap.at(xxMsq) = 1;
      break;
   case h6:
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mh = sortEigenvalues(es).at(0);
      flagMap.at(xxDmglst2) = 0;
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mhcut = sortEigenvalues(es).at(0);
      errors.push_back(fabs(Mh - Mhcut));
      flagMap.at(xxDmglst2) = 1;
      flagMap.at(xxMsq) = 0;
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mhcut = sortEigenvalues(es).at(0);
      errors.push_back(fabs(Mh - Mhcut));
      flagMap.at(xxMsq) = 1;
      break;
   case h6b:
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mh = sortEigenvalues(es).at(0);
      flagMap.at(xxDmglst2) = 0;
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mhcut = sortEigenvalues(es).at(0);
      errors.push_back(fabs(Mh - Mhcut));
      flagMap.at(xxDmglst2) = 1;
      flagMap.at(xxDmsqst2) = 0;
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mhcut = sortEigenvalues(es).at(0);
      errors.push_back(fabs(Mh - Mhcut));
      flagMap.at(xxDmsqst2) = 1;
      break;
   case h9:
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mh = sortEigenvalues(es).at(0);
      flagMap.at(xxDmsqst1) = 0;
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mhcut = sortEigenvalues(es).at(0);
      errors.push_back(fabs(Mh - Mhcut));
      flagMap.at(xxDmsqst1) = 1;
      flagMap.at(xxDmst12) = 0;
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mhcut = sortEigenvalues(es).at(0);
      errors.push_back(fabs(Mh - Mhcut));
      flagMap.at(xxDmst12) = 1;
      flagMap.at(xxMgl) = 0;
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mhcut = sortEigenvalues(es).at(0);
      errors.push_back(fabs(Mh - Mhcut));
      flagMap.at(xxMgl) = 1;
      break;
   }
   // evalue the sqrt of the squared errors
   double squaredErrorSum = 0.;
   for(auto const& error: errors){
      squaredErrorSum = squaredErrorSum + pow2(error);
   }
   // set the expansion depth for the next comparison
   flagMap.at(xxMst) = 0;
   flagMap.at(xx) = 0;
   return sqrt(squaredErrorSum);
}

/*
 * 	input parameters for the check
 */
himalaya::Parameters checkTermsXt33(){

   himalaya::Parameters pars;

   pars.scale = 1973.75;
   pars.mu = 1999.82;
   pars.g3 =  1.02907;
   pars.vd = 49.5751;
   pars.vu = 236.115;
   pars.mq2 <<  4.00428e+06 , 0, 0,
               0, 4.00428e+06, 0,
               0, 0, 3.99786e+06;
   pars.md2 << 4.00361e+06, 0, 0,
               0, 4.00361e+06, 0,
               0, 0, 4.00346e+06;
   pars.mu2 << 4.00363e+06 , 0, 0,
               0, 4.00363e+06, 0,
               0, 0, 3.99067e+06;
   pars.Ab = 9996.81;
   pars.At = 6992.34;

   pars.MA = 1992.14;
   pars.MG = 2000.96;
   pars.MW = 76.7777;
   pars.MZ = 88.4219;
   pars.Mt = 147.295;
   pars.Mb = 2.23149;
   pars.MSt << 1745.3 , 2232.1;
   pars.MSb <<  2000.14, 2001.09;
   pars.s2t = -0.999995;
   pars.s2b =-0.550527;

   return pars;
}

/*
 * 	check the expansion terms
 */
void himalaya::HierarchyCalculator::checkTerms(){
   p = checkTermsXt33();
   init();
   himalaya::HierarchyObject ho (false);
   for(int i = h3; i <= h9q2; i++){
      ho.setSuitableHierarchy(i);
      Eigen::Matrix2d oloMat = calculateHierarchy(ho, 1,0,0);
      Eigen::Matrix2d twloMat = calculateHierarchy(ho, 0,1,0);
      Eigen::Matrix2d thloMat = calculateHierarchy(ho, 0,0,1);
      bool ck1LPassed, ck2LPassed, ck3LPassed;
      switch(i){
	 case h3:
	    ck1LPassed = false;
	    ck2LPassed = false;
	    ck3LPassed = false;
	    ck1LPassed = oloMat(0,0) - (-1033.437882123761) < 1e-06 && oloMat(1,0) - (-394.3521101999062) < 1e-06 && oloMat(1,1) - 17633.47392819223 < 1e-06;
	    ck2LPassed = twloMat(0,0) - (-13.4248798129061) < 1e-06 && twloMat(1,0) - 10.91323060388626 < 1e-06 && twloMat(1,1) - 1477.154147584478 < 1e-06;
	    ck3LPassed = thloMat(0,0) - 1.163582002655875 < 1e-06 && thloMat(1,0) - 9.897241079348351 < 1e-06 && thloMat(1,1) - 369.9741236956309 < 1e-06;
	    std::cout << "Hierarchy " << i << " passed checks 1L: " << tf(ck1LPassed) << " 2L: " << tf(ck2LPassed) << " 3L: "<< tf(ck3LPassed) << "." << std::endl;
	 break;
	 case h32q2g:
	    ck1LPassed = false;
	    ck2LPassed = false;
	    ck3LPassed = false;
	    ck1LPassed = oloMat(0,0) - (-1033.437882123761) < 1e-06 && oloMat(1,0) - (-394.3521101999062) < 1e-06 && oloMat(1,1) - 17633.47392819223 < 1e-06;
	    ck2LPassed = twloMat(0,0) - (-13.66052379180129) < 1e-06 && twloMat(1,0) - 11.26755617866339 < 1e-06 && twloMat(1,1) - 1477.465656153518 < 1e-06;
	    ck3LPassed = thloMat(0,0) - 1.113051431370291 < 1e-06 && thloMat(1,0) - 9.903809573970422 < 1e-06 && thloMat(1,1) - 369.7408109643386 < 1e-06;
	    std::cout << "Hierarchy " << i << " passed checks 1L: " << tf(ck1LPassed) << " 2L: " << tf(ck2LPassed) << " 3L: "<< tf(ck3LPassed) << "." << std::endl;
	 break;
	 case h3q22g:
	    ck1LPassed = false;
	    ck2LPassed = false;
	    ck3LPassed = false;
	    ck1LPassed = oloMat(0,0) - (-1033.437882123761) < 1e-06 && oloMat(1,0) - (-394.3521101999062) < 1e-06 && oloMat(1,1) - 17633.47392819223 < 1e-06;
	    ck2LPassed = twloMat(0,0) - (-13.66052379180129) < 1e-06 && twloMat(1,0) - 11.26755617866339 < 1e-06 && twloMat(1,1) - 1477.465656153518 < 1e-06;
	    ck3LPassed = thloMat(0,0) - 1.058450932536496 < 1e-06 && thloMat(1,0) - 10.0141272838662 < 1e-06 && thloMat(1,1) - 370.3301180635573 < 1e-06;
	    std::cout << "Hierarchy " << i << " passed checks 1L: " << tf(ck1LPassed) << " 2L: " << tf(ck2LPassed) << " 3L: "<< tf(ck3LPassed) << "." << std::endl;
	 break;
	 case h4:
	    ck1LPassed = false;
	    ck2LPassed = false;
	    ck3LPassed = false;
	    ck1LPassed = oloMat(0,0) - 0 < 1e-06 && oloMat(1,0) - 0 < 1e-06 && oloMat(1,1) - 6685.123085628641 < 1e-06;
	    ck2LPassed = twloMat(0,0) - 0 < 1e-06 && twloMat(1,0) - 1183.325484493686 < 1e-06 && twloMat(1,1) - 1458.970501474495 < 1e-06;
	    ck3LPassed = thloMat(0,0) - 162.1379208650191 < 1e-06 && thloMat(1,0) - 326.0219627343553 < 1e-06 && thloMat(1,1) - 431.6926278454841 < 1e-06;
	    std::cout << "Hierarchy " << i << " passed checks 1L: " << tf(ck1LPassed) << " 2L: " << tf(ck2LPassed) << " 3L: "<< tf(ck3LPassed) << "." << std::endl;
	 break;
	 case h5:
	    ck1LPassed = false;
	    ck2LPassed = false;
	    ck3LPassed = false;
	    ck1LPassed = oloMat(0,0) - 15921.69462848581 < 1e-06 && oloMat(1,0) - (-388569.2043081555) < 1e-06 && oloMat(1,1) - 7874.401574063407 < 1e-06;
	    ck2LPassed = twloMat(0,0) - (-86.77887344841422) < 1e-06 && twloMat(1,0) - (-20625.63783863484) < 1e-06 && twloMat(1,1) - (-42446.62009872038) < 1e-06;
	    ck3LPassed = thloMat(0,0) - 2442.115080578889 < 1e-06 && thloMat(1,0) - (-3859.942907446577) < 1e-06 && thloMat(1,1) - 60593.055768119 < 1e-06;
	    std::cout << "Hierarchy " << i << " passed checks 1L: " << tf(ck1LPassed) << " 2L: " << tf(ck2LPassed) << " 3L: "<< tf(ck3LPassed) << "." << std::endl;
	 break;
	 case h5g1:
	    ck1LPassed = false;
	    ck2LPassed = false;
	    ck3LPassed = false;
	    ck1LPassed = oloMat(0,0) - 15921.69462848581 < 1e-06 && oloMat(1,0) - (-388569.2043081556) < 1e-06 && oloMat(1,1) - 7874.401574063407 < 1e-06;
	    ck2LPassed = twloMat(0,0) - (-114.6037388932203) < 1e-06 && twloMat(1,0) - (-20341.84471909946) < 1e-06 && twloMat(1,1) - (-42843.48046642416) < 1e-06;
	    ck3LPassed = thloMat(0,0) - 2415.507513838155 < 1e-06 && thloMat(1,0) - (-3766.750163753644) < 1e-06 && thloMat(1,1) - 59380.34497121828 < 1e-06;
	    std::cout << "Hierarchy " << i << " passed checks 1L: " << tf(ck1LPassed) << " 2L: " << tf(ck2LPassed) << " 3L: "<< tf(ck3LPassed) << "." << std::endl;
	 break;
	 case h6:
	    ck1LPassed = false;
	    ck2LPassed = false;
	    ck3LPassed = false;
	    ck1LPassed = oloMat(0,0) - 9272.477351702315 < 1e-06 && oloMat(1,0) - (-184.7601614832763) < 1e-06 && oloMat(1,1) - 7581.278122072418 < 1e-06;
	    ck2LPassed = twloMat(0,0) - (-1078.578574572312) < 1e-06 && twloMat(1,0) - 7096.529601647042 < 1e-06 && twloMat(1,1) - (-1927.791631086123) < 1e-06;
	    ck3LPassed = thloMat(0,0) - 245.4412216221288 < 1e-06 && thloMat(1,0) - 573.1296253278389 < 1e-06 && thloMat(1,1) - 8448.4582538127 < 1e-06;
	    std::cout << "Hierarchy " << i << " passed checks 1L: " << tf(ck1LPassed) << " 2L: " << tf(ck2LPassed) << " 3L: "<< tf(ck3LPassed) << "." << std::endl;
	 break;
	 case h6b:
	    ck1LPassed = false;
	    ck2LPassed = false;
	    ck3LPassed = false;
	    ck1LPassed = oloMat(0,0) - 9272.477351702311 < 1e-06 && oloMat(1,0) - (-184.7601614832763) < 1e-06 && oloMat(1,1) - 7581.278122072418 < 1e-06;
	    ck2LPassed = twloMat(0,0) - (-1078.578574572312) < 1e-06 && twloMat(1,0) - 7096.52960164704 < 1e-06 && twloMat(1,1) - (-1900.197036824461) < 1e-06;
	    ck3LPassed = thloMat(0,0) - 283.0253770519464 < 1e-06 && thloMat(1,0) - 566.2182257407396 < 1e-06 && thloMat(1,1) - 10093.33785879814 < 1e-06;
	    std::cout << "Hierarchy " << i << " passed checks 1L: " << tf(ck1LPassed) << " 2L: " << tf(ck2LPassed) << " 3L: "<< tf(ck3LPassed) << "." << std::endl;
	 break;
	 case h6b2qg2:
	    ck1LPassed = false;
	    ck2LPassed = false;
	    ck3LPassed = false;
	    ck1LPassed = oloMat(0,0) - 9272.477351702311 < 1e-06 && oloMat(1,0) - (-184.7601614832759) < 1e-06 && oloMat(1,1) - 7581.278122072418 < 1e-06;
	    ck2LPassed = twloMat(0,0) - (-1089.201418061661) < 1e-06 && twloMat(1,0) - 7145.267026465748 < 1e-06 && twloMat(1,1) - (-2077.345120153528) < 1e-06;
	    ck3LPassed = thloMat(0,0) - 285.3154791763894 < 1e-06 && thloMat(1,0) - 544.3654284413091 < 1e-06 && thloMat(1,1) - 10336.22756889787 < 1e-06;
	    std::cout << "Hierarchy " << i << " passed checks 1L: " << tf(ck1LPassed) << " 2L: " << tf(ck2LPassed) << " 3L: "<< tf(ck3LPassed) << "." << std::endl;
	 break;
	 case h6bq22g:
	    ck1LPassed = false;
	    ck2LPassed = false;
	    ck3LPassed = false;
	    ck1LPassed = oloMat(0,0) - 9272.477351702315 < 1e-06 && oloMat(1,0) - (-184.7601614832763) < 1e-06 && oloMat(1,1) - 7581.278122072418 < 1e-06;
	    ck2LPassed = twloMat(0,0) - (-1078.578574572311) < 1e-06 && twloMat(1,0) - 7096.529601647042 < 1e-06 && twloMat(1,1) - (-1900.197036824461) < 1e-06;
	    ck3LPassed = thloMat(0,0) - 283.0220052455883 < 1e-06 && thloMat(1,0) - 566.2190953470737 < 1e-06 && thloMat(1,1) - 10093.33986048966 < 1e-06;
	    std::cout << "Hierarchy " << i << " passed checks 1L: " << tf(ck1LPassed) << " 2L: " << tf(ck2LPassed) << " 3L: "<< tf(ck3LPassed) << "." << std::endl;;
	 break;
	 case h6bq2g2:
	    ck1LPassed = false;
	    ck2LPassed = false;
	    ck3LPassed = false;
	    ck1LPassed = oloMat(0,0) - 9272.477351702315 < 1e-06 && oloMat(1,0) - (-184.7601614832759) < 1e-06 && oloMat(1,1) - 7581.278122072418 < 1e-06;
	    ck2LPassed = twloMat(0,0) - (-1089.201418061661) < 1e-06 && twloMat(1,0) - 7145.267026465748 < 1e-06 && twloMat(1,1) - (-2077.345120153528) < 1e-06;
	    ck3LPassed = thloMat(0,0) - 285.3120881213721 < 1e-06 && thloMat(1,0) - 544.3662758149513 < 1e-06 && thloMat(1,1) - 10336.23012077387 < 1e-06;
	    std::cout << "Hierarchy " << i << " passed checks 1L: " << tf(ck1LPassed) << " 2L: " << tf(ck2LPassed) << " 3L: "<< tf(ck3LPassed) << "." << std::endl;
	 break;
	 case h6g2:
	    ck1LPassed = false;
	    ck2LPassed = false;
	    ck3LPassed = false;
	    ck1LPassed = oloMat(0,0) - 9272.477351702315 < 1e-06 && oloMat(1,0) - (-184.7601614832761) < 1e-06 && oloMat(1,1) - 7581.278122072418 < 1e-06;
	    ck2LPassed = twloMat(0,0) - (-1089.201418061661) < 1e-06 && twloMat(1,0) - 7145.267026465748 < 1e-06 && twloMat(1,1) - (-2112.642999123034) < 1e-06;
	    ck3LPassed = thloMat(0,0) - 246.0217489966267 < 1e-06 && thloMat(1,0) - 557.451210096066 < 1e-06 && thloMat(1,1) - 8628.076480526881 < 1e-06;
	    std::cout << "Hierarchy " << i << " passed checks 1L: " << tf(ck1LPassed) << " 2L: " << tf(ck2LPassed) << " 3L: "<< tf(ck3LPassed) << "." << std::endl;
	 break;
	 case h9:
	    ck1LPassed = false;
	    ck2LPassed = false;
	    ck3LPassed = false;
	    ck1LPassed = oloMat(0,0) - (-1033.437882123761) < 1e-06 && oloMat(1,0) - (-394.352110199906) < 1e-06 && oloMat(1,1) - 17633.47392819223 < 1e-06;
	    ck2LPassed = twloMat(0,0) - 420.2050380976995 < 1e-06 && twloMat(1,0) - (-554.6021924866435) < 1e-06 && twloMat(1,1) - (-797.8089039452509) < 1e-06;
	    ck3LPassed = thloMat(0,0) - 132.8584579769461 < 1e-06 && thloMat(1,0) - (-171.9326869339159) < 1e-06 && thloMat(1,1) - (-800.8408283898472) < 1e-06;
	    std::cout << "Hierarchy " << i << " passed checks 1L: " << tf(ck1LPassed) << " 2L: " << tf(ck2LPassed) << " 3L: "<< tf(ck3LPassed) << "." << std::endl;
	 break;
	 case h9q2:
	    ck1LPassed = false;
	    ck2LPassed = false;
	    ck3LPassed = false;
	    ck1LPassed = oloMat(0,0) - (-1033.437882123761) < 1e-06 && oloMat(1,0) - (-394.3521101999065) < 1e-06 && oloMat(1,1) - 17633.47392819223 < 1e-06;
	    ck2LPassed = twloMat(0,0) - 420.2050380976993 < 1e-06 && twloMat(1,0) - (-554.6021924866436) < 1e-06 && twloMat(1,1) - (-797.8089039452487) < 1e-06;
	    ck3LPassed = thloMat(0,0) - 132.6358855624267 < 1e-06 && thloMat(1,0) - (-171.4711818838455) < 1e-06 && thloMat(1,1) - (-800.9569014303727) < 1e-06;
	    std::cout << "Hierarchy " << i << " passed checks 1L: " << tf(ck1LPassed) << " 2L: " << tf(ck2LPassed) << " 3L: "<< tf(ck3LPassed) << "." << std::endl;
	 break;
      }
   }
}

/*
 * 	returns "true" or "false" with respect to the bool tf
 */
std::string himalaya::HierarchyCalculator::tf(const bool tf){
   return tf ? "true" : "false";
}

int himalaya::HierarchyCalculator::getCorrectHierarchy(const int hierarchy){
   if(hierarchy < 0 || hierarchy > 13){
      throw std::runtime_error("Hierarchy " + std::to_string(hierarchy) + " not included!");
   }
   return hierarchyMap.at(hierarchy);
}

void himalaya::HierarchyCalculator::printInfo(){
   //std::cout << ".  .         .          \n|__|*._ _  _.| _.  . _.\n|  ||[ | )(_]|(_]\\_|(_]\n                 ._|   " << std::endl;
   std::cout << "Himalaya contains code by: P. Slavich et al. (2-loop rMSSM Higgs self-energies)" << std::endl;
   std::cout << "Himalaya uses the 3-loop contributions of Kant et al." << std::endl;
}