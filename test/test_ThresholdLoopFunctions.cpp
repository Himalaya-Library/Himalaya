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


template <class Fn>
void test_3(const std::string& function_name, Fn fn, double eps)
{
   const auto filename = std::string(TEST_DATA_DIR) + PATH_SEPARATOR + "data" +
                         PATH_SEPARATOR + function_name + ".dat";

   struct Data { double x{}, y{}, z{}, f{}; };

   const auto data = read_data<double>(filename);

   for (auto d: data) {
      const double x = d.at(0), y = d.at(1), z = d.at(2), ex = d.at(3);
      const auto f = fn(x, y, z);
      INFO("expected: " << ex << ", observed: " << f);
      CHECK_CLOSE(ex, f, eps);
   }
}

TEST_CASE("test_3")
{
   test_3("Phi", himalaya::mh2_eft::threshold_loop_functions::phi_xyz, 1e-13);
}
