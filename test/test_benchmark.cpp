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


/// random number for variation of MS
auto rVar = [] { return random(0.9, 1.1); };


/// generate random parameter point of FO calculation
himalaya::Parameters setup_point(double MS, double tb, double xt)
{
   himalaya::Parameters pars;

   const double Xt = xt*MS;
   const double beta = std::atan(tb);
   pars.scale = MS*rVar();
   pars.mu = MS*rVar();
   pars.g1 = 0.46;
   pars.g2 = 0.65;
   pars.g3 = 1.166;
   pars.vd = 246*std::cos(beta);
   pars.vu = 246*std::sin(beta);
   pars.mq2 << sqr(MS*rVar()), 0, 0,
               0, sqr(MS*rVar()), 0,
               0, 0, sqr(MS*rVar());
   pars.md2 << sqr(MS*rVar()), 0, 0,
               0, sqr(MS*rVar()), 0,
               0, 0, sqr(MS*rVar());
   pars.mu2 << sqr(MS*rVar()), 0, 0,
               0, sqr(MS*rVar()), 0,
               0, 0, sqr(MS*rVar());
   pars.ml2 << sqr(MS*rVar()), 0, 0,
               0, sqr(MS*rVar()), 0,
               0, 0, sqr(MS*rVar());
   pars.me2 << sqr(MS*rVar()), 0, 0,
               0, sqr(MS*rVar()), 0,
               0, 0, sqr(MS*rVar());
   pars.Au << 0, 0, 0,
              0, 0, 0,
              0, 0, Xt + pars.mu/tb;
   pars.Ad << 0, 0, 0, 0, 0, 0, 0, 0, 0;
   pars.Ae << 0, 0, 0, 0, 0, 0, 0, 0, 0;
   pars.Yu << 0, 0, 0, 0, 0, 0, 0, 0, 0.862;
   pars.Yd << 0, 0, 0, 0 ,0 ,0 ,0 ,0, 0.133;
   pars.Ye << 0, 0, 0, 0, 0, 0, 0, 0, 0.101;
   pars.MA = MS*rVar();
   pars.M1 = MS*rVar();
   pars.M2 = MS*rVar();
   pars.MG = MS*rVar();

   pars.validate(verbose);

   return pars;
}


himalaya::Parameters make_point()
{
   return setup_point(rMS(), rTB(), rXT());
}


/// calculate all Higgs mass corrections
himalaya::HierarchyObject calculate_all(const himalaya::Parameters& point, bool isAlphab)
{
   himalaya::HierarchyObject ho(isAlphab);

   try {
      himalaya::HierarchyCalculator hc(point, verbose);
      ho = hc.calculateDMh3L(isAlphab);
   } catch (const std::exception& e) {
   }

   return ho;
}


/// avoid removal of value
template <class T>
inline void do_not_optimize(const T& value) {
   asm volatile("" : : "r,m"(value) : "memory");
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
   return time_in_milliseconds([N, f] {
      for (unsigned i = 0; i < N; i++) {
         try {
            do_not_optimize(f());
         } catch (...) {
         }
      }
   });
}

} // anonymous namespace


TEST_CASE("benchmark all loop corrections O(as^2*at^2)")
{
   const unsigned N = 1000;
   const auto time_in_ms = time_in_milliseconds(
      N, [] { return calculate_all(make_point(), false); });

   std::cout << "Average time per point: " << time_in_ms/N << " ms O(as^2*at^2)\n";
}


TEST_CASE("benchmark all loop corrections O(as^2*ab^2)")
{
   const unsigned N = 1000;
   const auto time_in_ms = time_in_milliseconds(
      N, [] { return calculate_all(make_point(), false); });

   std::cout << "Average time per point: " << time_in_ms/N << " ms O(as^2*ab^2)\n";
}
