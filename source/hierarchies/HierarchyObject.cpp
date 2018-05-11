#include "HierarchyObject.hpp"
#include "Hierarchies.hpp"
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
   if(loops > 0 && loops <=3){
      expUncertainties.emplace(loops, uncertainty);
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
   if(loops > 0 && loops <=3){
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
   if(loops >= 0 && loops <=3){
      dMhMap.emplace(loops, dMh);
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
   if(loops >= 0 && loops <=3){
      return dMhMap.at(loops);
   }
   
   throw std::runtime_error("Higgs mass matrix for " + std::to_string(loops) + " loop(s) is not available.");
   
}

/**
 * 	Sets the DR -> MDR shift
 * 	@param mdrShift the DR -> MDR shiftet matrix of the form M(MDR) - M(DR).
 */
void himalaya::HierarchyObject::setDRToMDRShift(const Eigen::Matrix2d& mdrShift){
   this -> mdrShift = mdrShift;
}

/**
 * 	@return The matrix M(MDR) - M(DR) at the order O(alpha_x + alpha_x*alpha_s)
 */
Eigen::Matrix2d himalaya::HierarchyObject::getDRToMDRShift() const{
   return mdrShift;
}

/**
 * 	Sets the MDR masses
 * 	@param mdrMasses a vector containting the MDR masses with the lightest particle at position 0.
 */
void himalaya::HierarchyObject::setMDRMasses(Eigen::Matrix<double, 2, 1>& mdrMasses){
   this -> mdrMasses = sortVector(mdrMasses);
}

/**
 * 	@return A vector of the MDR stop/sbottom masses. The 0th entry corresponds to the lighter particle.
 */
Eigen::Matrix<double, 2, 1> himalaya::HierarchyObject::getMDRMasses() const{
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
 * 	Sets the delta_lambda at 3-loop order with Himalaya logs.
 * 	@param deltaLambda delta_lambda at 3-loop order.
 */
void himalaya::HierarchyObject::setDeltaLambdaHimalaya(double deltaLambda){
   this -> deltaLambdaHimalaya = deltaLambda;
}

/**
 * 	@return 3-loop delta_lambda with Himalaya logs
 */
double himalaya::HierarchyObject::getDeltaLambdaHimalaya() const{
   return deltaLambdaHimalaya;
}

/**
 * 	Sets the delta_lambda at 3-loop order with EFT logs.
 * 	@param zeta delta_lambda at 3-loop order.
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
 * 	@param zeta constant part of delta_lambda at 3-loop order.
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
 * 	Sets the DR' -> MS shift for the 3-loop threshold correction which should be added to the DR' result
 * 	@param shift the DR' -> MS shift which should be added to the 3-loop threshold correction
 */
void himalaya::HierarchyObject::setDRbarPrimeToMSbarShift(double shift){
   drBarPrimeToMSbarShift = shift;
}

/**
 * 	@return the DR' -> MS shift for the 3-loop threshold correction which should be added to the DR' result
 */
double himalaya::HierarchyObject::getDRbarPrimeToMSbarShift() const{
   return drBarPrimeToMSbarShift;
}



/**
 * 	Sorts a vector.
 * 	@param vector The vector which should be sorted.
 * 	@return Returns a vector the lightest entry at position 0.
 */
Eigen::Matrix<double, 2, 1> himalaya::HierarchyObject::sortVector(Eigen::Matrix<double, 2, 1>& vector){
   // checks if all variables are ordered in the right way
   if (vector(0) > vector(1)) {
      std::swap(vector(0), vector(1));
   }
   return vector;
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
   const std::string mdrString = ho.getMDRFlag() == 0 ? "No MDR shifts applied" : "MDR shifts applied";
   const int suitableHierarchy = ho.getSuitableHierarchy();
   std::string renSchemeString = (ho.getRenormalizationScheme() == RenSchemes::H3m 
      || ho.getRenormalizationScheme() == RenSchemes::H3mMDRBAR) ? "H3m scheme" : "DR'";
   ostr << "===================================\n"
	<< "Himalaya HierarchyObject parameters\n"
        << "===================================\n"
	<< "Ren. scheme       =  " << renSchemeString << "\n"
	<< "MDR shifts?       =  " << mdrString << "\n"
        << "Hierarchy         =  " << suitableHierarchy << " (" << ho.getH3mHierarchyNotation(suitableHierarchy) << ")\n"
	<< "Mstop_1           =  " << ho.getMDRMasses()(0) << " GeV (" << renSchemeString << ")\n"
	<< "Mstop_2           =  " << ho.getMDRMasses()(1) << " GeV (" << renSchemeString << ")\n"
        << "Abs. diff 2L      =  " << ho.getAbsDiff2L() << " GeV\n"
        << "Rel. diff 2L      =  " << ho.getRelDiff2L()*100 << " %\n"
        << "Mh^2_tree         =  {{" << ho.getDMh(0).row(0)(0) << ", " << ho.getDMh(0).row(0)(1)
		   << "}, {" << ho.getDMh(0).row(1)(0) << ", " << ho.getDMh(0).row(1)(1) << "}} GeV^2\n"
        << "Mh^2_1L           =  {{" << ho.getDMh(1).row(0)(0) << ", " << ho.getDMh(1).row(0)(1)
		   << "}, {" << ho.getDMh(1).row(1)(0) << ", " << ho.getDMh(1).row(1)(1) << "}} GeV^2\n"
        << "Mh^2_2L           =  {{" << ho.getDMh(2).row(0)(0) << ", " << ho.getDMh(2).row(0)(1)
		   << "}, {" << ho.getDMh(2).row(1)(0) << ", " << ho.getDMh(2).row(1)(1) << "}} GeV^2\n"
        << "Mh^2_3L           =  {{" << ho.getDMh(3).row(0)(0) << ", " << ho.getDMh(3).row(0)(1)
		   << "}, {" << ho.getDMh(3).row(1)(0) << ", " << ho.getDMh(3).row(1)(1) << "}} GeV^2\n"
        << "Exp. uncert. 1L   =  " << ho.getExpUncertainty(1) << " GeV\n"
        << "Exp. uncert. 2L   =  " << ho.getExpUncertainty(2) << " GeV\n"
        << "Exp. uncert. 3L   =  " << ho.getExpUncertainty(3) << " GeV\n"
	<< "DR -> MDR shift   =  {{" << ho.getDRToMDRShift().row(0)(0) << ", " << ho.getDRToMDRShift().row(0)(1)
		   << "}, {" << ho.getDRToMDRShift().row(1)(0) << ", " << ho.getDRToMDRShift().row(1)(1)  << "}} GeV^2\n"
	<< "Δλ 3L Himalaya    =  " << ho.getDeltaLambdaHimalaya() << " (expanded coefficients of logarithms)\n"
        << "Δλ 3L EFT         =  " << ho.getDeltaLambdaEFT() << " (exact mass dependence of coefficients of logarithms)\n"
	<< "DR' -> MS shift   =  " << ho.getDRbarPrimeToMSbarShift();

   return ostr;
}

