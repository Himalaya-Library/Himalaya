#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"

#include "himalaya/HierarchyCalculator.hpp"
#include "stopwatch.hpp"

#include <iostream>
#include <random>
#include <utility>

namespace {

const bool verbose = false;


double sqr(double x) noexcept { return x*x; }


/// random number generator
double random(double start, double stop)
{
   static std::minstd_rand gen;
   std::uniform_real_distribution<double> dist(start, stop);
   return dist(gen);
}


/// random value for SUSY scale MS
auto rMS = [] {
   const double MS_min = 400;
   const double MS_max = 4000;
   return random(MS_min, MS_max);
};


/// random value for tan(beta)
auto rTB = [] { return random(2, 100); };


/// random value for x_t
auto rXT = [] { return random(-3.0, 3.0); };


/// generate random parameter point of FO calculation
himalaya::Parameters setup_point(double MS, double tb, double xt)
{
   himalaya::Parameters pars;

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
   pars.Au << 0, 0, 0,
              0, 0, 0,
              0, 0, Xt + pars.mu/tb;
   pars.Ad << 0, 0, 0, 0, 0, 0, 0, 0, 0;
   pars.Ae << 0, 0, 0, 0, 0, 0, 0, 0, 0;
   pars.Yu << 0, 0, 0, 0, 0, 0, 0, 0, 0.862;
   pars.Yd << 0, 0, 0, 0 ,0 ,0 ,0 ,0, 0.133;
   pars.Ye << 0, 0, 0, 0, 0, 0, 0, 0, 0.101;
   pars.MA = MS;
   pars.M1 = MS;
   pars.M2 = MS;
   pars.MG = MS;

   pars.validate(verbose);

   return pars;
}


himalaya::Parameters make_point()
{
   return setup_point(rMS(), rTB(), rXT());
}


/// calculate amu and uncertainty
void calculate(const himalaya::Parameters& point)
{
   try {
      himalaya::HierarchyCalculator hc(point, verbose);
      const auto ho = hc.calculateDMh3L(false);
   } catch (const std::exception& e) {
      std::cerr << e.what() << '\n';
   }
}


/// measure time for calling the function f()
template <class F>
double time_in_milliseconds(F&& f)
{
   himalaya::Stopwatch sw;
   sw.start();
   f();
   sw.stop();
   return sw.get_time_in_milliseconds();
}


/// measure time for calling the function f() N times
template <class F>
double time_in_milliseconds(unsigned N, F&& f)
{
   auto loop = [N, f]() {
      for (unsigned i = 0; i < N; i++) {
         try {
            (void) f();
         } catch (...) {
         }
      }
   };

   return time_in_milliseconds(loop);
}

} // anonymous namespace


TEST_CASE("benchmark fixed-order")
{
   const unsigned N = 1000;
   const auto time_in_ms = time_in_milliseconds(
      N, [] { return calculate(make_point()); });

   std::cout << "Average time per point: " << time_in_ms/N << " ms\n";
}
