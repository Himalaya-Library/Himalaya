#include "doctest.h"
#include "himalaya/mh2_eft/ThresholdLoopFunctions.hpp"
#include <algorithm>
#include <cmath>
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


double rel_diff(double x, double y)
{
   return (x - y)/std::max(std::fabs(x), std::fabs(y));
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


std::string create_filename(const std::string& function_name)
{
   return std::string(TEST_DATA_DIR) + PATH_SEPARATOR + "data" +
                      PATH_SEPARATOR + function_name + ".txt";
}


template <class Fn>
void test_1(const std::string& function_name, Fn fn, double eps)
{
   const auto data = read_data<double>(create_filename(function_name));

   for (auto d: data) {
      const auto ex = d.at(1);
      const auto f = fn(d.at(0));
      INFO(function_name << "( " << d.at(0) << " ): expected: " << ex << ", observed: " << f << ", rel. dev.: " << rel_diff(ex, f));
      CHECK_CLOSE(ex, f, eps);
   }
}


template <class Fn>
void test_2(const std::string& function_name, Fn fn, double eps)
{
   const auto data = read_data<double>(create_filename(function_name));

   for (auto d: data) {
      const auto ex = d.at(2);
      const auto f = fn(d.at(0), d.at(1));
      INFO(function_name << "( " << d.at(0) << ", " << d.at(1) << "): expected: " << ex << ", observed: " << f << ", rel. dev.: " << rel_diff(ex, f));
      CHECK_CLOSE(ex, f, eps);
   }
}


template <class Fn>
void test_3(const std::string& function_name, Fn fn, double eps)
{
   const auto data = read_data<double>(create_filename(function_name));

   for (auto d: data) {
      const auto ex = d.at(3);
      const auto f = fn(d.at(0), d.at(1), d.at(2));
      INFO(function_name << "( " << d.at(0) << ", " << d.at(1) << ", " << d.at(2) << "): expected: " << ex << ", observed: " << f << ", rel. dev.: " << rel_diff(ex, f));
      CHECK_CLOSE(ex, f, eps);
   }
}


TEST_CASE("test_1")
{
   test_1("f1_", himalaya::mh2_eft::threshold_loop_functions::f1, 1e-11);
   test_1("f2_", himalaya::mh2_eft::threshold_loop_functions::f2, 1e-11);
   test_1("f3_", himalaya::mh2_eft::threshold_loop_functions::f3, 1e-11);
   test_1("f4_", himalaya::mh2_eft::threshold_loop_functions::f4, 1e-11);

   test_1("F1", himalaya::mh2_eft::threshold_loop_functions::F1, 1e-11);
   test_1("F2", himalaya::mh2_eft::threshold_loop_functions::F2, 1e-11);
   test_1("F3", himalaya::mh2_eft::threshold_loop_functions::F3, 1e-11);
   test_1("F4", himalaya::mh2_eft::threshold_loop_functions::F4, 1e-11);
   test_1("F5", himalaya::mh2_eft::threshold_loop_functions::F5, 1e-11);
   test_1("F6", himalaya::mh2_eft::threshold_loop_functions::F6, 1e-11);
   test_1("F7", himalaya::mh2_eft::threshold_loop_functions::F7, 1e-11);
}


TEST_CASE("test_2")
{
   test_2("f5_", himalaya::mh2_eft::threshold_loop_functions::f5, 1e-8);
   test_2("f6_", himalaya::mh2_eft::threshold_loop_functions::f6, 1e-10);
   test_2("f7_", himalaya::mh2_eft::threshold_loop_functions::f7, 1e-7);
   test_2("f8_", himalaya::mh2_eft::threshold_loop_functions::f8, 1e-9);

   test_2("F8", himalaya::mh2_eft::threshold_loop_functions::F8, 1e-8);
   test_2("F9", himalaya::mh2_eft::threshold_loop_functions::F9, 1e-11);
}


TEST_CASE("test_3")
{
   test_3("Phi", himalaya::mh2_eft::threshold_loop_functions::phi_xyz, 1e-13);
}
