#include "doctest.h"
#include "himalaya/mh2_fo/PV.hpp"
#include <fstream>
#include <limits>
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


struct A0 {
   double m2{}, q2{}, a0{};
};


struct B0 {
   double p2{}, m12{}, m22{}, q2{}, b0{};
};


std::ostream& operator<<(std::ostream& ostr, const A0& a0)
{
   ostr << std::setprecision(std::numeric_limits<double>::digits10)
        << "A0(m2=" << a0.m2 << ", q2=" << a0.q2 << ") = " << a0.a0;
   return ostr;
}


std::ostream& operator<<(std::ostream& ostr, const B0& b0)
{
   ostr << std::setprecision(std::numeric_limits<double>::digits10)
        << "B0(p2=" << b0.p2 << ", m12=" << b0.m12
        << ", m22=" << b0.m22 << ", q2=" << b0.q2 << ") = " << b0.b0;
   return ostr;
}


/// read A0 function data
std::vector<A0> read_a0(const std::string& filename)
{
   std::ifstream fstr(filename);
   std::string line;
   std::vector<A0> data;

   while (std::getline(fstr, line)) {
      double m2{}, q2{}, a0{};

      std::istringstream isstr(line);
      isstr >> m2 >> q2 >> a0;

      data.emplace_back(A0{m2, q2, a0});
   }

   return data;
}


/// read B0 function data
std::vector<B0> read_b0(const std::string& filename)
{
   std::ifstream fstr(filename);
   std::string line;
   std::vector<B0> data;

   while (std::getline(fstr, line)) {
      double p2{}, m12{}, m22{}, q2{}, b0{};

      std::istringstream isstr(line);
      isstr >> p2 >> m12 >> m22 >> q2 >> b0;

      data.emplace_back(B0{p2, m12, m22, q2, b0});
   }

   return data;
}


/// test A0
TEST_CASE("test_A0")
{
   const auto filename = std::string(TEST_DATA_DIR) + PATH_SEPARATOR + "data" +
                         PATH_SEPARATOR + "A0.dat";
   const auto data = read_a0(filename);
   const double eps = 1e-12;

   for (auto d: data) {
      const auto a0 = himalaya::mh2_fo::a0(d.m2, d.q2);
      INFO("expected: " << d);
      INFO("observed: " << a0);
      CHECK_CLOSE(d.a0, a0, eps);
   }
}


/// test B0 for m1^2 == m2^2
TEST_CASE("test_B0xx")
{
   const auto filename = std::string(TEST_DATA_DIR) + PATH_SEPARATOR + "data" +
                         PATH_SEPARATOR + "B0xx.dat";
   const auto data = read_b0(filename);
   const double eps = 1e-12;

   for (auto d: data) {
      const auto b0 = himalaya::mh2_fo::b0xx(d.p2, d.m12, d.q2);
      INFO("expected: " << d);
      INFO("observed: " << b0);
      CHECK_CLOSE(d.b0, b0, eps);
   }
}


/// test B0
TEST_CASE("test_B0")
{
   const auto filename = std::string(TEST_DATA_DIR) + PATH_SEPARATOR + "data" +
                         PATH_SEPARATOR + "B0.dat";
   const auto data = read_b0(filename);
   const double eps = 1e-12;

   for (auto d: data) {
      const auto b0 = himalaya::mh2_fo::b0(d.p2, d.m12, d.m22, d.q2);
      INFO("expected: " << d);
      INFO("observed: " << b0);
      CHECK_CLOSE(d.b0, b0, eps);
   }
}
