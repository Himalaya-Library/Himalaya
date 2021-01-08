#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "himalaya/HierarchyCalculator.hpp"
#include <cmath>
#include <limits>
#include <utility>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>

#define CHECK_CLOSE(a,b,eps) CHECK((a) == doctest::Approx(b).epsilon(eps))


namespace {

const bool verbose = false;


struct Point {
   double MS;
   double xt;
   double tb;
};


struct Data {
   double Mh;
};


std::ostream& operator<<(std::ostream& ostr, const Data& data)
{
   ostr << "Data["
        << data.Mh << "]";
   return ostr;
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

   pars.validate(verbose);

   return pars;
}


const std::pair<Point, Data> points[] = {
   { {1000.0, 0.0, 20.0}, {96.8024604087} }
};


Data calculate_all(const himalaya::Parameters& point)
{
   himalaya::HierarchyCalculator hc(point);
   const auto ho = hc.calculateDMh3L(false);
   const auto dMh_0L = ho.getDMh(0);
   const auto dMh_1L = ho.getDMh(3);
   const auto dMh_2L = ho.getDMh(3);
   const auto dMh_3L = ho.getDMh(3);

   Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> solver(dMh_0L + dMh_1L + dMh_2L + dMh_3L);

   Data data{ std::sqrt(solver.eigenvalues()(0)) };

   return data;
}


} // anonymous namespace


TEST_CASE("test_points")
{
   const double eps  = 1e-12;

   for (const auto& p: points) {
      const auto point = make_point(p.first);
      const auto data = calculate_all(point);
      INFO("point = " << point);
      CHECK_CLOSE(p.second.Mh, data.Mh, eps);
   }
}
