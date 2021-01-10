#include "doctest.h"
#include "himalaya/HierarchyCalculator.hpp"
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <utility>
#include <vector>
#include <Eigen/Core>


#define CHECK_CLOSE(a,b,eps) CHECK((a) == doctest::Approx(b).epsilon(eps))


namespace {

const char PATH_SEPARATOR =
#ifdef _WIN32
   '\\';
#else
   '/';
#endif

const std::string DATA_FILE = std::string(TEST_DATA_DIR) + PATH_SEPARATOR + "test_points.txt";
const bool VERBOSE = false;
const bool IS_ALPHA_B = false; // calculate O(as^2*ab^2) corrections
const int N_DIGITS = 6; // number of digits to test


struct Point {
   double MS;
   double xt;
   double tb;
};


struct Data {
   double MhFO{};
   double MhEFT{};
};


std::ostream& operator<<(std::ostream& ostr, const Point& point)
{
   ostr << point.MS << '\t' << point.xt << '\t' << point.tb;
   return ostr;
}


std::istream& operator>>(std::istream& istr, Point& point)
{
   istr >> point.MS >> point.xt >> point.tb;
   return istr;
}


std::ostream& operator<<(std::ostream& ostr, const Data& data)
{
   ostr << data.MhFO << '\t' << data.MhEFT;
   return ostr;
}


std::istream& operator>>(std::istream& istr, Data& data)
{
   istr >> data.MhFO >> data.MhEFT;
   return istr;
}


himalaya::Parameters make_point(const Point& point)
{
   himalaya::Parameters pars;

   const double MS = point.MS;
   const double xt = point.xt;
   const double tb = point.tb;
   const double MS2 = MS*MS;
   const double Xt = xt*MS;
   const double beta = std::atan(tb);

   pars.scale = MS;
   pars.mu = MS;
   pars.g1 = 0.46;
   pars.g2 = 0.65;
   pars.g3 = 1.166;
   pars.vd = 246*std::cos(beta);
   pars.vu = 246*std::sin(beta);
   pars.mq2 << MS2, 0, 0,
               0, MS2, 0,
               0, 0, MS2;
   pars.md2 << MS2, 0, 0,
               0, MS2, 0,
               0, 0, MS2;
   pars.mu2 << MS2, 0, 0,
               0, MS2, 0,
               0, 0, MS2;
   pars.ml2 << MS2, 0, 0,
               0, MS2, 0,
               0, 0, MS2;
   pars.me2 << MS2, 0, 0,
               0, MS2, 0,
               0, 0, MS2;
   pars.Au(2,2) = Xt + pars.mu/tb;
   pars.Yu(2,2) = 0.862;
   pars.Yd(2,2) = 0.133;
   pars.Ye(2,2) = 0.101;
   pars.MA = MS;
   pars.M1 = MS;
   pars.M2 = MS;
   pars.MG = MS;

   return pars;
}


Data calculate_all(const himalaya::Parameters& point)
{
   Data data;

   try {
      himalaya::HierarchyCalculator hc(point, VERBOSE);
      const auto ho = hc.calculateDMh3L(IS_ALPHA_B);

      data.MhFO  = std::sqrt(ho.getDMh2FO(0)
                             + ho.getDMh2FOAt(1)
                             + ho.getDMh2FOAt(2)
                             + ho.getDMh2FOAt(3));
      data.MhEFT = std::sqrt(ho.getDMh2EFTAt(0)
                             + ho.getDMh2EFTAt(1)
                             + ho.getDMh2EFTAt(2)
                             + ho.getDMh2EFTAt(3));
   } catch (const std::exception& e) {
      std::cerr << e.what() << '\n';
   }

   return data;
}


void write_points_with_xt_tb(std::ostream& ostr, double xt, double tb)
{
   const unsigned N = 100;
   const double MS_start = 200;
   const double MS_stop  = 10000 + MS_start;

   for (unsigned i = 0; i < N; i++) {
      const double MS = MS_start + i*(MS_stop - MS_start)/N;
      Point point{MS, xt, tb};

      const auto data = calculate_all(make_point(point));
      ostr << std::setprecision(N_DIGITS+1) << point << '\t' << data << '\n';
   }
}


void write_points(const std::string& filename)
{
   std::ofstream ostr(filename);
   write_points_with_xt_tb(ostr,             0.0, 20.0);
   write_points_with_xt_tb(ostr, -std::sqrt(6.0), 20.0);
   write_points_with_xt_tb(ostr,  std::sqrt(6.0), 20.0);
}


std::vector<std::pair<Point, Data>> read_points(const std::string& filename)
{
   std::ifstream istr(filename);
   std::vector<std::pair<Point, Data>> vec;
   std::string line;

   while (std::getline(istr, line)) {
      Point point;
      Data data;

      std::istringstream isstr(line);
      isstr >> point >> data;

      vec.push_back({point, data});
   }

   return vec;
}


} // anonymous namespace


// TEST_CASE("write_points")
// {
//    write_points(DATA_FILE);
// }


TEST_CASE("test_points")
{
   const double eps  = std::pow(10.0, -N_DIGITS);
   const auto points = read_points(DATA_FILE);

   for (const auto& p: points) {
      const auto point = make_point(p.first);
      const auto data = calculate_all(point);
      INFO("point =\n" << point);
      CHECK_CLOSE(p.second.MhFO , data.MhFO , eps);
      CHECK_CLOSE(p.second.MhEFT, data.MhEFT, eps);
   }
}
