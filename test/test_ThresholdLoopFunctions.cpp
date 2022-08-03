#include "doctest.h"
#include "himalaya/mh2_eft/ThresholdLoopFunctions.hpp"
#include <algorithm>
#include <fstream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>


#define CHECK_CLOSE(a,b,eps) CHECK((a) == doctest::Approx(b).epsilon(eps))


const char PATH_SEPARATOR =
#ifdef _WIN32
   '\\';
#else
   '/';
#endif


struct Phi {
   double x{}, y{}, z{}, phi{};
};


std::ostream& operator<<(std::ostream& ostr, const Phi& phi)
{
   ostr << std::setprecision(std::numeric_limits<double>::digits10)
        << "Phi(x=" << phi.x << ", y=" << phi.y << ", z=" << phi.z << ") = " << phi.phi;
   return ostr;
}


/// tranform vector of elements of type A -> B
template <class B, class A, class F>
std::vector<B> map(F f, const std::vector<A>& in)
{
   std::vector<B> out;
   std::transform(in.cbegin(), in.cend(), std::back_inserter(out), f);
   return out;
}


/// read data line-wise from file
template <class T>
std::vector<std::vector<T>> read_data(const std::string& filename)
{
   std::ifstream fstr(filename);
   std::string line;
   std::vector<std::vector<T>> data;

   while (std::getline(fstr, line)) {
      std::istringstream isstr(line);
      data.emplace_back((std::istream_iterator<T>(isstr)),
                        std::istream_iterator<T>());
   }

   return data;
}


/// read Phi(x,y,z) function data
std::vector<Phi> read_phi(const std::string& filename)
{
   return map<Phi>(
      [](const auto& d) {
         return Phi{d.at(0), d.at(1), d.at(2), d.at(3)};
      },
      read_data<double>(filename));
}


/// test Phi(x,y,z)
TEST_CASE("test_Phi")
{
   const auto filename = std::string(TEST_DATA_DIR) + PATH_SEPARATOR + "data" +
                         PATH_SEPARATOR + "Phi.dat";
   const auto data = read_phi(filename);
   const double eps = 1e-13;

   for (auto d: data) {
      const auto phi = himalaya::mh2_eft::threshold_loop_functions::phi_xyz(d.x, d.y, d.z);
      INFO("expected: " << d);
      INFO("observed: " << phi);
      CHECK_CLOSE(d.phi, phi, eps);
   }
}
