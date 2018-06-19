// ====================================================================
// This file is part of Himalaya.
//
// Himalaya is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================

#include "HierarchyObject.hpp"
#include "Hierarchies.hpp"
#include <cmath>
#include <iostream>

/**
 * 	A constructor.
 * 	@param isAlphab the boolean which determines wether the members are proportinal to alpha_b or alpha_t.
 */
himalaya::HierarchyObject::HierarchyObject(bool isAlphab)
   : isAlphab(isAlphab)
{
}

/**
 * 	@return The value of isAlphab.
 */
bool himalaya::HierarchyObject::getIsAlphab() const{
   return isAlphab;
}

/**
 * 	Sets the suitable hierarchy
 * 	@param hierarchy the integer key of the hierarchy.
 */
void himalaya::HierarchyObject::setSuitableHierarchy(int hierarchy){
   this -> hierarchy = hierarchy;
}

/**
 * 	@return The key to the suitable hierarchy
 */
int himalaya::HierarchyObject::getSuitableHierarchy() const{
   return hierarchy;
}

/**
 * 	Sets the absolute difference of the Higgs masses at two-loop level
 * 	@param absDiff2L the absolute difference of the Higgs masses as a double.
 */
void himalaya::HierarchyObject::setAbsDiff2L(double absDiff2L){
   this -> absDiff2L = absDiff2L;
}

/**
 * 	@return The absolute difference of the exact and expanded Higgs masses at two-loop level at the order O(alpha_x + alpha_x*alpha_s).
 */
double himalaya::HierarchyObject::getAbsDiff2L() const{
   return absDiff2L;
}

/**
 * 	Sets the relative difference ot the Higgs masses at two-loop level
 * 	@param relDiff2L the relative difference of the Higgs masses as a double.
 */
void himalaya::HierarchyObject::setRelDiff2L(double relDiff2L){
   this -> relDiff2L = relDiff2L;
}

/**
 * 	@return The relative difference of the exact and expanded Higgs masses at two-loop level at the order O(alpha_x + alpha_x*alpha_s).
 */
double himalaya::HierarchyObject::getRelDiff2L() const{
   return relDiff2L;
}


/**
 * 	Sets the uncertainty of the expansion at a given loop level.
 * 	@param loops the integer value of the corresponding loops. Can be 1, 2 or 3.
 * 	@param uncertainty the expansion untertainty at the given loop order as a double.
 */
void himalaya::HierarchyObject::setExpUncertainty(int loops, double uncertainty){
   if(loops > 0 && loops <= 3){
      expUncertainties[loops] = uncertainty;
   }
   else {
      throw std::runtime_error("Expansion uncertainty for " + std::to_string(loops) + " loop(s) is not available.");
   }
}

/**
 * 	@param loops an integer which can be 1, 2 or 3.
 * 	@return A double which is the expansion uncertainty for the given loop order.
 */
double himalaya::HierarchyObject::getExpUncertainty(int loops) const{
   if(loops > 0 && loops <= 3){
      return expUncertainties.at(loops);
   }
   
   throw std::runtime_error("Expansion uncertainty for " + std::to_string(loops) + " loop(s) is not available.");
   
}

/**
 * 	Sets the delta of the CP-even Higgs mass matrix
 * 	@param loops the integer value of the corresponding loops. Can be 0, 1, 2 or 3. 0 corresponds to the tree-level.
 * 	@param dMh the delta of the mass matrix.
 */
void himalaya::HierarchyObject::setDMh(int loops, const Eigen::Matrix2d& dMh){
   if(loops >= 0 && loops <= 3){
      dMhMap[loops] = dMh;
   }
   else {
      throw std::runtime_error("Higgs mass matrix for " + std::to_string(loops) + " loop(s) is not available.");
   }
}

/**
 * 	@param loops an integer which can be 0, 1, 2, 3. Here, 0 corresponds to the tree-level matrix.
 * 	@return The CP-even Higgs mass matrix at the given loop order.
 */
Eigen::Matrix2d himalaya::HierarchyObject::getDMh(int loops) const{
   if(loops >= 0 && loops <= 3){
      return dMhMap.at(loops);
   }
   
   throw std::runtime_error("Higgs mass matrix for " + std::to_string(loops) + " loop(s) is not available.");
   
}

/**
 * 	Sets the DR' -> MDR' shift
 * 	@param mdrShift the DR' -> MDR' shiftet matrix of the form M(MDR') - M(DR').
 */
void himalaya::HierarchyObject::setDRbarPrimeToMDRbarPrimeShift(const Eigen::Matrix2d& mdrShift){
   this -> mdrShift = mdrShift;
}

/**
 * 	@return The matrix M(MDR') - M(DR') at the order O(alpha_x + alpha_x*alpha_s)
 */
Eigen::Matrix2d himalaya::HierarchyObject::getDRbarPrimeToMDRbarPrimeShift() const{
   return mdrShift;
}

void himalaya::HierarchyObject::setDRbarPrimeToH3mShift(const Eigen::Matrix2d shift){
   h3mShift = shift;
}

Eigen::Matrix2d himalaya::HierarchyObject::getDRbarPrimeToH3mShift() const{
   return h3mShift;
}


/**
 * 	Sets the MDR masses
 * 	@param mdrMasses a vector containting the MDR masses with the lightest particle at position 0.
 */
void himalaya::HierarchyObject::setMDRMasses(Eigen::Vector2d& mdrMasses){
   this -> mdrMasses = sortVector(mdrMasses);
}

/**
 * 	@return A vector of the MDR stop/sbottom masses. The 0th entry corresponds to the lighter particle.
 */
Eigen::Vector2d himalaya::HierarchyObject::getMDRMasses() const{
   return mdrMasses;
}

/**
 * 	Sets the mdrFlag to calculate the corretions in without MDR (0) or with MDR (1) shifts.
 * 	@param mdrFlag an int. (0) for H3m (DR')- and (1) for MDR-scheme.
 * 	@throws runtime_exception if the flag is neither 0 or 1 an exception is thrown.
 */
void himalaya::HierarchyObject::setMDRFlag(int mdrFlag){
   if(mdrFlag != 0 && mdrFlag != 1) {
      throw std::runtime_error("The MDR-flag has to be 0 (DR-scheme) or 1 (MDR-scheme). Input: " + std::to_string(mdrFlag) + ".");
}
   this -> mdrFlag = mdrFlag;
}


/**
 * 	@return the MDRFlag integer.
 */
int himalaya::HierarchyObject::getMDRFlag() const{
   return mdrFlag;
}

/**
 * 	Sets the renormalization scheme accodring to the RenScheme enum
 * 	@param renScheme an int according to the RenScheme enum.
 * 	@throws runtime_exception if the flag is not in {0,1,2,3} an exception is thrown
 */
void himalaya::HierarchyObject::setRenormalizationScheme(int renScheme){
   if(renScheme < 0 || renScheme > 3) {
      throw std::runtime_error("The renormalization scheme has to be 0 (H3m), 1 (DR'), 2 (H3m with MDR), 3 (MDR'). Input: " + std::to_string(renScheme) + ".");
   }
   renormalizationScheme = renScheme;
}

/**
 * 	@return Renormalization scheme key
 */
int himalaya::HierarchyObject::getRenormalizationScheme() const{
   return renormalizationScheme;
}


/**
 * 	Sets the delta_lambda at 3-loop order with H3m logs.
 * 	@param deltaLambda delta_lambda at 3-loop order.
 */
void himalaya::HierarchyObject::setDeltaLambdaH3m(double deltaLambda){
   this -> deltaLambdaH3m = deltaLambda;
}

/**
 * 	@return 3-loop delta_lambda with H3m logs
 */
double himalaya::HierarchyObject::getDeltaLambdaH3m() const{
   return deltaLambdaH3m;
}

/**
 * 	Sets the delta_lambda at 3-loop order with EFT logs.
 * 	@param deltaLambda delta_lambda at 3-loop order.
 */
void himalaya::HierarchyObject::setDeltaLambdaEFT(double deltaLambda){
   this -> deltaLambdaEFT = deltaLambda;
}

/**
 * 	@return 3-loop delta_lambda with EFT logs
 */
double himalaya::HierarchyObject::getDeltaLambdaEFT() const{
   return deltaLambdaEFT;
}

/**
 * 	Sets the constant part of delta_lambda at 3-loop order.
 * 	@param deltaLambda constant part of delta_lambda at 3-loop order.
 */
void himalaya::HierarchyObject::setDeltaLambdaNonLog(double deltaLambda){
   deltaLambdaNonLog = deltaLambda;
}

/**
 *        @return 3-loop delta_lambda for the degenerated mass case with EFT logs
 */
double himalaya::HierarchyObject::getDeltaLambdaNonLog() const{
   return deltaLambdaNonLog;
}

/**
 * 	Sets the Xt parts of the uncertainty of delta_lambda_H3m
 *        @param uncertainty of 3-loop delta_lambda
 */
void himalaya::HierarchyObject::setDeltaLambdaXtUncertaintyH3m(double uncertainty){
   deltaLambdaXtUncertaintyH3m = uncertainty;
}

/**
 *        @return uncertainty of 3-loop delta_lambda_H3m
 */
double himalaya::HierarchyObject::getDeltaLambdaUncertaintyH3m() const {
   return std::abs(getDeltaLambdaEFT() - getDeltaLambdaH3m()) 
      + std::abs(deltaLambdaXtUncertaintyH3m)
      + std::abs(expansionUncertaintyDeltaLambda);
}

/**
 * 	Sets the Xt parts of the uncertainty of delta_lambda_EFT
 *      @param uncertainty of 3-loop delta_lambda
 */
void himalaya::HierarchyObject::setDeltaLambdaXtUncertaintyEFT(double uncertainty){
   deltaLambdaXtUncertaintyEFT = uncertainty;
}

/**
 *	@return uncertainty of 3-loop delta_lambda_EFT
 *	TODO: the difference should only be taken up to O(Xt^3/Xt^4) depending on the chosen hierarchy
 */
double himalaya::HierarchyObject::getDeltaLambdaUncertaintyEFT() const {
   return std::abs(deltaLambdaXtUncertaintyEFT)
      + std::abs(expansionUncertaintyDeltaLambda);
}

/**
 * 	Sets the DR' -> MS shift for delta_lambda_himalaya which should be added to the DR' result
 * 	@param shift the DR' -> MS shift which should be added to the 3-loop threshold correction
 */
void himalaya::HierarchyObject::setDRbarPrimeToMSbarShiftH3m(double shift){
   drBarPrimeToMSbarShiftH3m = shift;
}

/**
 * 	@return the DR' -> MS shift for delta_lambda_himalaya which should be added to the DR' result
 */
double himalaya::HierarchyObject::getDRbarPrimeToMSbarShiftH3m() const{
   return drBarPrimeToMSbarShiftH3m;
}

/**
 * 	Sets the DR' -> MS shift for delta_lambda_EFT which should be added to the DR' result
 * 	@param shift the DR' -> MS shift which should be added to the 3-loop threshold correction
 */
void himalaya::HierarchyObject::setDRbarPrimeToMSbarShiftEFT(double shift){
   drBarPrimeToMSbarShiftEFT = shift;
}

/**
 * 	@return the DR' -> MS shift for delta_lambda_EFT which should be added to the DR' result
 */
double himalaya::HierarchyObject::getDRbarPrimeToMSbarShiftEFT() const{
   return drBarPrimeToMSbarShiftEFT;
}

/**
 * 	Set the expasion uncertainty for delta_lambda
 * 	@param expUncertLambda the expansion uncertainty for delta_lambda
 */
void himalaya::HierarchyObject::setExpUncertaintyDeltaLambda(double expUncertLambda){
   expansionUncertaintyDeltaLambda = expUncertLambda;
}


/**
 * 	Sorts a vector.
 * 	@param vector The vector which should be sorted.
 * 	@return Returns a vector the lightest entry at position 0.
 */
Eigen::Vector2d himalaya::HierarchyObject::sortVector(Eigen::Vector2d& vector){
   // checks if all variables are ordered in the right way
   if (vector(0) > vector(1)) {
      std::swap(vector(0), vector(1));
   }
   return vector;
}

double himalaya::HierarchyObject::getDeltaLambda0L() const{
   return deltaLambda0L;
}

double himalaya::HierarchyObject::getDeltaLambda1L() const{
   return deltaLambda1L;
}

void himalaya::HierarchyObject::setDeltaLambda0L(double deltaLambda){
   deltaLambda0L = deltaLambda;
}

void himalaya::HierarchyObject::setDeltaLambda1L(double deltaLambda){
   deltaLambda1L = deltaLambda;
}

double himalaya::HierarchyObject::getDeltaLambda2L() const{
   return deltaLambda2L;
}

void himalaya::HierarchyObject::setDeltaLambda2L(double deltaLambda){
   deltaLambda2L = deltaLambda;
}

double himalaya::HierarchyObject::getDRbarPrimeToMSbarShiftDeltaLambda1L() const{
   return drToMSDL1L;
}

void himalaya::HierarchyObject::setDRbarPrimeToMSbarShiftDeltaLambda1L(double shift){
   drToMSDL1L = shift;
}

double himalaya::HierarchyObject::getDRbarPrimeToMSbarShiftDeltaLambda2L() const{
   return drToMSDL2L;
}

void himalaya::HierarchyObject::setDRbarPrimeToMSbarShiftDeltaLambda2L(double shift){
   drToMSDL2L = shift;
}

double himalaya::HierarchyObject::getDeltaMh2EFT0L() const{
   return mh2EFT0L;
}

void himalaya::HierarchyObject::setDeltaMh2EFT0L(double deltaMh2){
   mh2EFT0L = deltaMh2;
}

double himalaya::HierarchyObject::getDeltaMh2EFT1L() const{
   return mh2EFT1L;
}

void himalaya::HierarchyObject::setDeltaMh2EFT1L(double deltaMh2){
   mh2EFT1L = deltaMh2;
}

double himalaya::HierarchyObject::getDeltaMh2EFT2L() const{
   return mh2EFT2L;
}

void himalaya::HierarchyObject::setDeltaMh2EFT2L(double deltaMh2){
   mh2EFT2L = deltaMh2;
}

double himalaya::HierarchyObject::getDeltaMh2EFT3L() const{
   return mh2EFT3L;
}

void himalaya::HierarchyObject::setDeltaMh2EFT3L(double deltaMh2){
   mh2EFT3L = deltaMh2;
}

/**
 *      Returns the H3m notation of a given hierarchy.
 *      @param hierarchy An integer of a Himalaya hierarchy.
 *      @return Returns the corresponding H3m notation of the given hierarchy as a string.
 */
std::string himalaya::HierarchyObject::getH3mHierarchyNotation(int hierarchy) const{
   switch (hierarchy){
      case Hierarchies::h3:
	 return "h3";
      case Hierarchies::h32q2g:
	 return "h32q2g";
      case Hierarchies::h3q22g:
	 return "h3q22g";
      case Hierarchies::h4:
	 return "h4";
      case Hierarchies::h5:
	 return "h5";
      case Hierarchies::h5g1:
	 return "h5g1";
      case Hierarchies::h6:
	 return "h6";
      case Hierarchies::h6b:
	 return "h6b";
      case Hierarchies::h6b2qg2:
	 return "h6b2qg2";
      case Hierarchies::h6bq22g:
	 return "h6bq22g";
      case Hierarchies::h6bq2g2:
	 return "h6bq2g2";
      case Hierarchies::h6g2:
	 return "h6g2";
      case Hierarchies::h9:
	 return "h9";
      case Hierarchies::h9q2:
	 return "h9q2";
      default:
	 return "Hierarchy " + std::to_string(hierarchy) + " not included";
   }
}

/**
 * 	Prints out all information of the HierarchyObject
 */
std::ostream& himalaya::operator<<(std::ostream& ostr, himalaya::HierarchyObject const &ho){
   const int suitableHierarchy = ho.getSuitableHierarchy();
   const std::string renSchemeString = (ho.getRenormalizationScheme() == RenSchemes::H3m 
      || ho.getRenormalizationScheme() == RenSchemes::H3mMDRBAR) ? "H3m scheme" : "DR'";
   const std::string massString = ho.getIsAlphab() ? "Msbottom" : "Mstop";
   const std::string spaces = ho.getIsAlphab() == true ? "             " : "                ";
   ostr << "===================================\n"
	<< "Himalaya HierarchyObject parameters\n"
        << "===================================\n"
	<< "Ren. scheme            =  " << renSchemeString << "\n"
        << "Hierarchy              =  " << suitableHierarchy << " (" << ho.getH3mHierarchyNotation(suitableHierarchy) << ")\n"
	<< massString << "_1" << spaces << "=  " << ho.getMDRMasses()(0) << " GeV (MDR')\n"
	<< massString << "_2" << spaces << "=  " << ho.getMDRMasses()(1) << " GeV (MDR')\n"
        << "Abs. diff 2L           =  " << ho.getAbsDiff2L() << " GeV\n"
        << "Rel. diff 2L           =  " << ho.getRelDiff2L()*100 << " %\n"
        << "Mh^2_0L                =  {{" << ho.getDMh(0).row(0)(0) << ", " << ho.getDMh(0).row(0)(1)
		   << "}, {" << ho.getDMh(0).row(1)(0) << ", " << ho.getDMh(0).row(1)(1) << "}} GeV^2\n"
        << "ΔMh^2_1L               =  {{" << ho.getDMh(1).row(0)(0) << ", " << ho.getDMh(1).row(0)(1)
		   << "}, {" << ho.getDMh(1).row(1)(0) << ", " << ho.getDMh(1).row(1)(1) << "}} GeV^2\n"
        << "ΔMh^2_2L               =  {{" << ho.getDMh(2).row(0)(0) << ", " << ho.getDMh(2).row(0)(1)
		   << "}, {" << ho.getDMh(2).row(1)(0) << ", " << ho.getDMh(2).row(1)(1) << "}} GeV^2\n"
        << "ΔMh^2_3L               =  {{" << ho.getDMh(3).row(0)(0) << ", " << ho.getDMh(3).row(0)(1)
		   << "}, {" << ho.getDMh(3).row(1)(0) << ", " << ho.getDMh(3).row(1)(1) << "}} GeV^2\n"
        << "Exp. uncert. 1L        =  " << ho.getExpUncertainty(1) << " GeV\n"
        << "Exp. uncert. 2L        =  " << ho.getExpUncertainty(2) << " GeV\n"
        << "Exp. uncert. 3L        =  " << ho.getExpUncertainty(3) << " GeV\n"
	<< "DR' -> MDR' shift      =  {{" << ho.getDRbarPrimeToMDRbarPrimeShift().row(0)(0) << ", " << ho.getDRbarPrimeToMDRbarPrimeShift().row(0)(1)
		   << "}, {" << ho.getDRbarPrimeToMDRbarPrimeShift().row(1)(0) << ", " << ho.getDRbarPrimeToMDRbarPrimeShift().row(1)(1)  << "}} GeV^2\n"
	<< "DR' -> H3m shift       =  {{" << ho.getDRbarPrimeToH3mShift().row(0)(0) << ", " << ho.getDRbarPrimeToH3mShift().row(0)(1)
		   << "}, {" << ho.getDRbarPrimeToH3mShift().row(1)(0) << ", " << ho.getDRbarPrimeToH3mShift().row(1)(1) << "}} GeV^2\n"
	<< "Δλ_0L                  =  " << ho.getDeltaLambda0L() << " O(g_1^2, g_2^2)\n"
	<< "Δλ_0L DR' -> MS shift  =  " << 0. << "\n"
	<< "Δλ_1L                  =  " << ho.getDeltaLambda1L() << " O(α_t)\n"
	<< "Δλ_1L DR' -> MS shift  =  " << ho.getDRbarPrimeToMSbarShiftDeltaLambda1L() << "\n"
	<< "Δλ_2L                  =  " << ho.getDeltaLambda2L() << " O(α_t*α_s)\n"
	<< "Δλ_2L DR' -> MS shift  =  " << ho.getDRbarPrimeToMSbarShiftDeltaLambda2L() << "\n"
	<< "Δλ_H3m                 =  " << ho.getDeltaLambdaH3m() << " +/- " << ho.getDeltaLambdaUncertaintyH3m() << " (expanded coefficients of logarithms)\n"
	<< "Δλ_H3m DR' -> MS shift =  " << ho.getDRbarPrimeToMSbarShiftH3m() << " (should be added to Δλ to convert it to MS)\n"
        << "Δλ_EFT                 =  " << ho.getDeltaLambdaEFT() << " +/- " << ho.getDeltaLambdaUncertaintyEFT() << " (exact mass dependence of coefficients of logarithms)\n"
	<< "Δλ_EFT DR' -> MS shift =  " << ho.getDRbarPrimeToMSbarShiftEFT() << " (should be added to Δλ to convert it to MS)\n"
	<< "Mh^2_EFT_0L            =  " << ho.getDeltaMh2EFT0L() << " GeV^2\n"
	<< "ΔMh^2_EFT_1L           =  " << ho.getDeltaMh2EFT1L() << " GeV^2\n"
	<< "ΔMh^2_EFT_2L           =  " << ho.getDeltaMh2EFT2L() << " GeV^2\n"
	<< "ΔMh^2_EFT_3L           =  " << ho.getDeltaMh2EFT3L() << " GeV^2";

   return ostr;
}

