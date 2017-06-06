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

void himalaya::HierarchyObject::setAbsDiff2l(const double& absDiff2l){
   this -> absDiff2l = absDiff2l;
}

double himalaya::HierarchyObject::getAbsDiff2l() const{
   return absDiff2l;
}

void himalaya::HierarchyObject::setRelDiff2l(const double& relDiff2l){
   this -> relDiff2l = relDiff2l;
}

double himalaya::HierarchyObject::getRelDiff2l() const{
   return relDiff2l;
}

void himalaya::HierarchyObject::setExpUncertainty(const int& loops, const double& uncertainty){
   expUncertainties.insert(std::pair<int, double> (loops, uncertainty));
}

double himalaya::HierarchyObject::getExpUncertainty(const int& loops) const{
   if(loops > 0 && loops <=3){
      return expUncertainties.at(loops);
   }
   else {
      throw std::runtime_error("Expansion uncertainty for " + std::to_string(loops) + " loop(s) is not available.");
   }
}

void himalaya::HierarchyObject::setDMh3l(const Eigen::Matrix2d& dMh3l){
   this -> dMh3l = dMh3l;
}

Eigen::Matrix2d himalaya::HierarchyObject::getDMh3L() const{
   return this -> dMh3l;
}

void himalaya::HierarchyObject::setDRToMDRShift(const Eigen::Matrix2d& mdrShift){
   this -> mdrShift = mdrShift;
}

Eigen::Matrix2d himalaya::HierarchyObject::getDRToMDRShift() const{
   return this -> mdrShift;
}

void himalaya::HierarchyObject::setDMh2l(const Eigen::Matrix2d& dMh2l){
   this -> dMh2l = dMh2l;
}

Eigen::Matrix2d himalaya::HierarchyObject::getDMh2L() const{
   return dMh2l;
}

void himalaya::HierarchyObject::setDMh1l(const Eigen::Matrix2d& dMh1l){
   this -> dMh1l = dMh1l;
}

Eigen::Matrix2d himalaya::HierarchyObject::getDMh1l() const{
   return dMh1l;
}

void himalaya::HierarchyObject::setDMh0l(const Eigen::Matrix2d& dMh0l){
   this -> dMh0l = dMh0l;
}

Eigen::Matrix2d himalaya::HierarchyObject::getDMh0l() const{
   return dMh0l;
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
