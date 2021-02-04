#include "doctest.h"
#include "himalaya/mh2_fo/PV.hpp"
#include <algorithm>
#include <fstream>
#include <iterator>
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


struct DB0 {
   double m12{}, m22{}, db0{};
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


std::ostream& operator<<(std::ostream& ostr, const DB0& db0)
{
   ostr << std::setprecision(std::numeric_limits<double>::digits10)
        << "DB0(m12=" << db0.m12 << ", m22=" << db0.m22 << ") = " << db0.db0;
   return ostr;
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


/// tranform vector of elements of type A -> B
template <class B, class A, class F>
std::vector<B> transform_to(const std::vector<A>& in, F f)
{
   std::vector<B> out;
   std::transform(in.cbegin(), in.cend(), std::back_inserter(out), f);
   return out;
}


/// read A0 function data
std::vector<A0> read_a0(const std::string& filename)
{
   return transform_to<A0>(read_data<double>(filename), [](const auto& d) {
      return A0{d.at(0), d.at(1), d.at(2)};
   });
}


/// read B0 function data
std::vector<B0> read_b0(const std::string& filename)
{
   return transform_to<B0>(read_data<double>(filename), [](const auto& d) {
      return B0{d.at(0), d.at(1), d.at(2), d.at(3), d.at(4)};
   });
}


/// read DB0 function data
std::vector<DB0> read_db0(const std::string& filename)
{
   return transform_to<DB0>(read_data<double>(filename), [](const auto& d) {
      return DB0{d.at(0), d.at(1), d.at(2)};
   });
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


/// test DB0
TEST_CASE("test_DB0")
{
   const auto filename = std::string(TEST_DATA_DIR) + PATH_SEPARATOR + "data" +
                         PATH_SEPARATOR + "DB0.dat";
   const auto data = read_db0(filename);
   const double eps = 1e-12;

   for (auto d: data) {
      const auto db0 = himalaya::mh2_fo::d1_b0(d.m12, d.m22);
      INFO("expected: " << d);
      INFO("observed: " << db0);
      CHECK_CLOSE(d.db0, db0, eps);
   }
}
