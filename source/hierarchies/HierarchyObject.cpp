#include <HierarchyObject.hpp>

himalaya::HierarchyObject::HierarchyObject(const bool& isAlphab){
   this -> isAlphab = isAlphab;
}

bool himalaya::HierarchyObject::getIsAlphab() const{
   return isAlphab;
}

void himalaya::HierarchyObject::setSuitableHierarchy(const int& hierarchy){
   this -> hierarchy = hierarchy;
}

int himalaya::HierarchyObject::getSuitableHierarchy() const{
   return hierarchy;
}

void himalaya::HierarchyObject::setAbsDiff2L(const double& absDiff2L){
   this -> absDiff2L = absDiff2L;
}

double himalaya::HierarchyObject::getAbsDiff2L() const{
   return absDiff2L;
}

void himalaya::HierarchyObject::setRelDiff2L(const double& relDiff2L){
   this -> relDiff2L = relDiff2L;
}

double himalaya::HierarchyObject::getRelDiff2L() const{
   return relDiff2L;
}

void himalaya::HierarchyObject::setExpUncertainty(const int& loops, const double& uncertainty){
   if(loops > 0 && loops <=3){
      expUncertainties.insert(std::pair<int, double> (loops, uncertainty));
   }
   else {
      throw std::runtime_error("Expansion uncertainty for " + std::to_string(loops) + " loop(s) is not available.");
   }
}

double himalaya::HierarchyObject::getExpUncertainty(const int& loops) const{
   if(loops > 0 && loops <=3){
      return expUncertainties.at(loops);
   }
   else {
      throw std::runtime_error("Expansion uncertainty for " + std::to_string(loops) + " loop(s) is not available.");
   }
}

void himalaya::HierarchyObject::setDMh(const int& loops, const Eigen::Matrix2d& dMh){
   if(loops >= 0 && loops <=3){
      dMhMap.insert(std::pair<int, Eigen::Matrix2d> (loops, dMh));
   }
   else {
      throw std::runtime_error("Higgs mass matrix for " + std::to_string(loops) + " loop(s) is not available.");
   }
}

Eigen::Matrix2d himalaya::HierarchyObject::getDMh(const int& loops) const{
   if(loops >= 0 && loops <=3){
      return dMhMap.at(loops);
   }
   else {
      throw std::runtime_error("Higgs mass matrix for " + std::to_string(loops) + " loop(s) is not available.");
   }
}

void himalaya::HierarchyObject::setDRToMDRShift(const Eigen::Matrix2d& mdrShift){
   this -> mdrShift = mdrShift;
}

Eigen::Matrix2d himalaya::HierarchyObject::getDRToMDRShift() const{
   return this -> mdrShift;
}


void himalaya::HierarchyObject::setMDRMasses(Eigen::Matrix<double,2,1>& mdrMasses){
   this -> mdrMasses = sortVector(mdrMasses);
}

Eigen::Matrix<double,2,1> himalaya::HierarchyObject::getMDRMasses() const{
   return this -> mdrMasses;
}

Eigen::Matrix<double,2,1> himalaya::HierarchyObject::sortVector(Eigen::Matrix<double,2,1>& vector){
   // checks if all variables are ordered in the right way
   if (vector(0) > vector(1)) {
      std::swap(vector(0), vector(1));
   }
   return vector;
}
