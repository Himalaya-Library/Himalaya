// ====================================================================
// This file is part of FlexibleSUSY.
//
// FlexibleSUSY is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// FlexibleSUSY is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with FlexibleSUSY.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

#include "threshold_loop_functions.hpp"
#include "dilog.hpp"

#include <cmath>
#include <limits>

/**
 * @file threshold_loop_functions.cpp
 * @brief Implementation of the loop functions from [arXiv:1407.4081]
 * @note The file has been taken from FlexibleSUSY.
 */

namespace himalaya {
namespace threshold_loop_functions {

namespace {
   const double Pi = 3.1415926535897932384626433832795;

   template <typename T> T sqr(T x) noexcept { return x*x; }
   template <typename T> T cube(T x) noexcept { return x*x*x; }
   template <typename T> T quad(T x) noexcept { return x*x*x*x; }
   template <typename T> T pow5(T x) noexcept { return x*x*x*x*x; }
   template <typename T> T pow6(T x) noexcept { return x*x*x*x*x*x; }
   template <typename T> T pow7(T x) noexcept { return x*x*x*x*x*x*x; }
   template <typename T> T pow8(T x) noexcept { return x*x*x*x*x*x*x*x; }
   template <typename T> T pow9(T x) noexcept { return x*x*x*x*x*x*x*x*x; }
   template <typename T> T power10(T x) noexcept { return x*x*x*x*x*x*x*x*x*x; }

   template <typename T>
   bool is_zero(T a, T prec = std::numeric_limits<T>::epsilon()) noexcept
   {
      return std::fabs(a) < prec;
   }

   template <typename T>
   bool is_equal(T a, T b, T prec = std::numeric_limits<T>::epsilon()) noexcept
   {
      return is_zero(a - b, prec);
   }

   template <typename T>
   bool is_equal_rel(T a, T b, T prec = std::numeric_limits<T>::epsilon()) noexcept
   {
      if (is_equal(a, b, std::numeric_limits<T>::epsilon()))
         return true;

      if (std::fabs(a) < std::numeric_limits<T>::epsilon())
         return is_equal(a, b, prec);

      return std::fabs((a - b)/a) < prec;
   }

} // anonymous namespace

double F1(double x) noexcept
{
   const double x2 = sqr(x);

   if (is_equal(x, 0.))
      return 0.;

   if (is_equal(x, 1., 0.01))
      return (16 + 41*x - 44*x2 + 21*cube(x) - 4*quad(x))/30.;

   return x*std::log(x2)/(x2-1);
}

double F2(double x) noexcept
{
   const double x2 = sqr(x);

   if (is_equal(x, 0.))
      return 0.;

   if (is_equal(x, 1., 0.01))
      return (-5 + 216*x - 226*x2 + 104*cube(x) - 19*quad(x))/70.;

   return 6*x2*(2-2*x2+(1+x2)*std::log(x2))/cube(x2-1);
}

double F3(double x) noexcept
{
   const double x2 = sqr(x);

   if (is_equal(x, 0.))
      return 0.;

   if (is_equal(x, 1., 0.01))
      return (-27 + 218*x - 142*x2 + 48*cube(x) - 7*quad(x))/90.;

   return 2*x*(5*(1-x2)+(1+4*x2)*std::log(x2))/(3*sqr(x2-1));
}

double F4(double x) noexcept
{
   const double x2 = sqr(x);

   if (is_equal(x, 0.))
      return 0.;

   if (is_equal(x, 1., 0.01))
      return (31 + 22*x - 42*x2 + 24*cube(x) - 5*quad(x))/30.;

   return 2*x*(x2-1-std::log(x2))/sqr(x2-1);
}

double F5(double x) noexcept
{
   const double x2 = sqr(x);
   const double x4 = quad(x);

   if (is_equal(x, 0.))
      return 0.;

   if (is_equal(x, 1., 0.01))
      return (13 + 165*x - 174*x2 + 81*cube(x) - 15*x4)/70.;

   if (is_equal(x, -1., 0.01))
     return (-13 + 165*x + 174*x2 + 81*cube(x) + 15*x4)/70.;

   return 3*x*(1-x4+2*x2*std::log(x2))/cube(1-x2);
}

double F6(double x) noexcept
{
   const double x2 = sqr(x);

   if (is_equal(x, 0.))
      return -0.75;

   if (is_equal(x, 1., 0.01))
      return (-103 + 128*x - 26*x2 + quad(x))/120.;

   if (is_equal(x, -1., 0.01))
      return (-103 - 128*x - 26*x2 + quad(x))/120.;

   return (x2-3)/(4*(1-x2)) + x2*(x2-2)/(2*sqr(1.-x2))*std::log(x2);
}

double F7(double x) noexcept
{
   const double x2 = sqr(x);
   const double x4 = quad(x);

   if (is_equal(x, 0.))
      return -1.5;

   if (is_equal(x, 1., 0.01))
      return -1.8642857142857139 + 2.057142857142856*x
         + 1.7142857142857153*x2 - 1.1428571428571432*cube(x)
         + 0.2357142857142858*x4;

   if (is_equal(x, -1., 0.01))
      return -1.8642857142857139 - 2.057142857142856*x
         + 1.7142857142857153*x2 + 1.1428571428571432*cube(x)
         + 0.2357142857142858*x4;

   return (-3*(x4-6*x2+1.))/(2*sqr(x2-1))
      + (3*x4*(x2-3.))/(cube(x2-1.))*std::log(x2);
}

/// F8(x1,x2) in the limit x1 -> 1 and x2 -> 1
static double F8_1_1(double x1, double x2) noexcept
{
   return 1. + 1.333333333333333*(-1 + x1) +
      (1.333333333333333 - 0.6666666666666661*(-1 + x1))*(-1 + x2);
}

/// F8(x1,x2) in the limit x1 -> 1
static double F8_1_x2(double x1, double x2) noexcept
{
   const double lx22 = std::log(sqr(x2));

   return -2. + (4.*(-1. + x1)*(-0.5 + 2.*sqr(x2) - 1.5*quad(x2)
                                + quad(x2)*lx22))
      /cube(-1. + sqr(x2))
      + (2. - 2.*sqr(x2) + 2.*quad(x2)*lx22)
      /sqr(-1. + sqr(x2))
      - (1.3333333333333333*cube(-1. + x1)*(
            -sqr(x2) - 9.*quad(x2) + 9.*pow6(x2)
            + pow8(x2)
            + (-6.*quad(x2) - 6.*pow6(x2))*lx22))
      /pow5(-1. + sqr(x2))
      + (2.*sqr(-1. + x1)*(
            -0.16666666666666666 + 1.5*sqr(x2) + 1.5*quad(x2)
            - 2.8333333333333335*pow6(x2)
            + (3.*quad(x2) + pow6(x2))*lx22))
      /quad(-1. + sqr(x2))
      + (0.26666666666666666*quad(-1. + x1)*(
            0.25 + 2.5*sqr(x2) + 80.*quad(x2) - 47.5*pow6(x2)
            - 36.25*pow8(x2) + power10(x2)
            + (37.5*quad(x2) + 75.*pow6(x2)
               + 7.5*pow8(x2))*lx22))
      /pow6(-1. + sqr(x2));
}

/// F8(x1,x2) in the limit x1 -> 0
static double F8_0_x2(double x1, double x2) noexcept
{
   const double lx22 = std::log(sqr(x2));

   return -2. + (2.*sqr(x1)*lx22)/(-1. + sqr(x2)) +
      (2.*sqr(x2)*lx22)/(-1. + sqr(x2));
}

// F8(x1,x2) in the limit x1 -> x2
static double F8_x1_x2(double x1, double x2) noexcept
{
   const double lx22 = std::log(sqr(x2));

   return -2. + (2.*(x1 - x2)*(3.*x2 - 4.*cube(x2) + pow5(x2) +
        2.*x2*lx22))/cube(-1. + sqr(x2)) +
   (2.*sqr(x2)*(-1. + sqr(x2) - 2.*lx22 +
        sqr(x2)*lx22))/sqr(-1. + sqr(x2)) -
   (0.33333333333333326*sqr(x1 - x2)*
      (17.000000000000007 - 9.000000000000007*sqr(x2) -
        9.*quad(x2) + pow6(x2) +
        6.000000000000002*lx22 +
        18.000000000000004*sqr(x2)*lx22))/
    quad(-1. + sqr(x2)) -
   (1.333333333333333*cube(x1 - x2)*
      (-1. - 9.000000000000002*sqr(x2) +
        9.000000000000002*quad(x2) + pow6(x2) -
        6.*sqr(x2)*lx22 -
        6.000000000000003*quad(x2)*lx22))/
    (x2*pow5(-1. + sqr(x2)));
}

double F8(double x1, double x2) noexcept
{
   if (is_equal(x1, 0.) && is_equal(x2, 0.))
      return -2.;

   if (is_equal(std::fabs(x1), 1., 0.01) && is_equal(std::fabs(x2), 1., 0.01))
      return F8_1_1(std::fabs(x1), std::fabs(x2));

   if (is_equal(std::fabs(x1), 1., 0.01)) {
      if (is_equal(x2, 0., 0.01)) {
         return -2.333333333333332 + 5.66666666666667*sqr(x2) +
            std::fabs(x1)*(2.6666666666666643 - 5.333333333333339*sqr(x2)) +
            sqr(x1)*(-0.33333333333333215 + 1.6666666666666696*sqr(x2));
      }

      return F8_1_x2(std::fabs(x1), x2);
   }

   if (is_equal(std::fabs(x2), 1., 0.01)) {
      if (is_equal(x1, 0., 0.01)) {
         return -2.3333333333333335 + 2.6666666666666665*std::fabs(x2) -
            0.3333333333333333*sqr(x2) +
            sqr(x1)*(5.666666666666667 - 5.333333333333334*std::fabs(x2) +
                     1.6666666666666667*sqr(x2));
      }

      return F8_1_x2(std::fabs(x2), x1);
   }

   if (is_equal(x1, 0., 0.0001))
      return F8_0_x2(x1, x2);

   if (is_equal(x2, 0., 0.0001))
      return F8_0_x2(x2, x1);

   if (is_equal(std::fabs(x1), std::fabs(x2), 0.00001))
      return F8_x1_x2(std::fabs(x1), std::fabs(x2));

   const double x12 = sqr(x1);
   const double x22 = sqr(x2);

   return -2. + 2./(x12-x22)
      *(quad(x1)/(x12-1.)*std::log(x12)
        -quad(x2)/(x22-1.)*std::log(x22));
}

/// F9(x1,x2) in the limit x1 -> 1 and x2 -> 1
static double F9_1_1(double x1, double x2) noexcept
{
   return 8.223809523809523 - 12.863492063492064*x2
      + 10.580952380952382*sqr(x2) - 4.609523809523809*cube(x2)
      + 0.8349206349206348*quad(x2)
      + x1*(-12.863492063492064 + 26.260317460317456*x2
            - 24.609523809523807*sqr(x2) + 11.530158730158728*cube(x2)
            - 2.184126984126984*quad(x2))
      + cube(x1)*(-4.60952380952381 + 11.53015873015873*x2
                        - 11.961904761904762*sqr(x2)
                        + 5.942857142857143*cube(x2)
                        - 1.1682539682539683*quad(x2))
      + quad(x1)*(0.8349206349206351 - 2.1841269841269844*x2
                        + 2.3190476190476197*sqr(x2)
                        - 1.1682539682539685*cube(x2)
                        + 0.23174603174603178*quad(x2))
      + sqr(x1)*(10.580952380952379 - 24.609523809523804*x2
                 + 24.6047619047619*sqr(x2)
                 - 11.96190476190476*cube(x2)
                 + 2.319047619047619*quad(x2));
}

/// F9(x1,x2) in the limit x1 -> 1
static double F9_1_x2(double x1, double x2) noexcept
{
   const double lx22 = std::log(sqr(x2));

   return (-2.*(-1. + x1)*cube(-1. + sqr(x2))*(
              -1. + quad(x2) - 2.*sqr(x2)*lx22)
           + 2.*quad(-1. + sqr(x2))*(
              1. - sqr(x2) + sqr(x2)*lx22)
           - 1.3333333333333333*cube(-1. + x1)*(
              -1. + sqr(x2))*(
                 -1. - 9.*sqr(x2) + 9.*quad(x2) + pow6(x2)
                 + (-6.*sqr(x2) - 6.*quad(x2))*lx22)
           + 0.3333333333333333*sqr(-1. + x1)*sqr(-1. + sqr(x2))
           *(5. + 15.*sqr(x2) - 21.*quad(x2) + pow6(x2)
             + (18.*sqr(x2) + 6.*quad(x2))*lx22)
           - 0.06666666666666667*quad(-1. + x1)*(
              -16. - 305.*sqr(x2) + 170.*quad(x2) + 160.*pow6(x2)
              - 10.*pow8(x2) + power10(x2)
              + (-150.*sqr(x2) - 300.*quad(x2)
                 - 30.*pow6(x2))*lx22))
      /pow6(-1. + sqr(x2));
}

/// F9(x1,x2) in the limit x1 -> 0
static double F9_0_x2(double, double x2) noexcept
{
   return (2.*std::log(sqr(x2)))/(-1. + sqr(x2));
}

/// F9(x1,x2) in the limit x1 -> x2
static double F9_x1_x2(double x1, double x2) noexcept
{
   const double lx22 = std::log(sqr(x2));

   return (2.*(-1. + sqr(x2) - lx22))/
      sqr(-1. + sqr(x2)) -
      (2.*(x1 - x2)*(-1. + quad(x2) -
                     2.*sqr(x2)*lx22))/
      (x2*cube(-1. + sqr(x2))) -
      (1.3333333333333333*cube(x1 - x2)*
       (-1.0000000000000002 - 9.*sqr(x2) + 9.*quad(x2) +
        pow6(x2) - 6.*sqr(x2)*lx22 -
        6.*quad(x2)*lx22))/
      (x2*pow5(-1. + sqr(x2))) +
      (1.6666666666666665*sqr(x1 - x2)*
       (0.2000000000000001 - 4.200000000000001*sqr(x2) +
        3.000000000000001*quad(x2) + pow6(x2) -
        1.2000000000000002*sqr(x2)*lx22 -
        3.6000000000000005*quad(x2)*lx22))/
      (sqr(x2)*quad(-1. + sqr(x2)));
}

double F9(double x1, double x2) noexcept
{
   if (is_equal(std::fabs(x1), 1., 0.01) && is_equal(std::fabs(x2), 1., 0.01))
      return F9_1_1(std::fabs(x1), std::fabs(x2));

   if (is_equal(std::fabs(x1), 1., 0.01)) {
      if (is_equal(x2, 0., 0.01))
         return 2. - 2.*(-1 + std::fabs(x1)) + 5./3.*sqr(-1 + std::fabs(x1));

      return F9_1_x2(std::fabs(x1), x2);
   }

   if (is_equal(std::fabs(x2), 1., 0.01)) {
      if (is_equal(x1, 0., 0.01))
         return 2. - 2.*(-1 + std::fabs(x2)) + 5./3.*sqr(-1 + std::fabs(x2));

      return F9_1_x2(std::fabs(x2), x1);
   }

   if (is_equal(x1, 0., 0.0001))
      return F9_0_x2(x1, x2);

   if (is_equal(x2, 0., 0.0001))
      return F9_0_x2(x2, x1);

   if (is_equal(std::fabs(x1), std::fabs(x2), 0.00001))
      return F9_x1_x2(std::fabs(x1), std::fabs(x2));

   const double x12 = sqr(x1);
   const double x22 = sqr(x2);

   return 2./(x12-x22)*(x12/(x12-1.)*std::log(x12)-x22/(x22-1.)*std::log(x22));
}

double f(double r) noexcept
{
   return F5(r);
}

double g(double r) noexcept
{
   return F7(r);
}

double f1(double r) noexcept
{
   if (is_equal(r, 0., 0.01))
      return 18./7.*sqr(r);

   // The function is even under x-> -x
   if (is_equal(r, 1., 0.01))
      return (-81 + 464*r + 270*sqr(r)
              - 208*cube(r) + 45*quad(r))/490.;

   // Notice the flipped sign for the odd terms
   if (is_equal(r, -1., 0.01))
      return (-81 - 464*r + 270*sqr(r)
              + 208*cube(r) + 45*quad(r))/490.;

   const double r2 = sqr(r);

   return (6*(r2+3)*r2)/(7*sqr(r2-1))
      + (6*(r2-5)*quad(r)*std::log(r2))/(7*cube(r2-1));
}

double f2(double r) noexcept
{
   if (is_equal(r, 0., 0.01))
      return 22./9.*sqr(r);

   // The function is even under x -> -x
   if (is_equal(r, 1., 0.01))
      return (-285 + 1616*r + 1230*sqr(r)
              - 848*cube(r) + 177*quad(r))/1890.;

   // Notice the flipped sign for the odd terms
   if (is_equal(r, -1., 0.01))
      return (-285 - 1616*r + 1230*sqr(r)
              + 848*cube(r) + 177*quad(r))/1890.;

   const double r2 = sqr(r);

   return (2*(r2+11)*r2)/(9*sqr(r2-1))
      + (2*(5*r2-17)*quad(r)*std::log(r2))/(9*cube(r2-1));
}

double f3(double r) noexcept
{
   if (is_equal(r, 0., 0.001))
      return 4./3.;

   if (is_equal(r, 1., 0.01))
      return (849 - 1184*r + 1566*sqr(r)
              - 736*cube(r) + 135*quad(r))/630.;

   if (is_equal(r, -1., 0.01))
      return (849 + 1184*r + 1566*sqr(r)
              + 736*cube(r) + 135*quad(r))/630.;

   const double r2 = sqr(r);
   const double r4 = quad(r);

   return (2*(r4+9*r2+2))/(3*sqr(r2-1))
      + (2*(r4-7*r2-6)*r2*std::log(r2))/(3*cube(r2-1));
}

double f4(double r) noexcept
{
   if (is_equal(r, 0., 0.001))
      return 12./7.;

   const double r2 = sqr(r);
   const double r4 = quad(r);

   if (is_equal(r, 1., 0.01))
      return (2589 - 3776*r + 4278*r2 - 1984*cube(r) + 363*r4)/1470.;

   if (is_equal(r, -1., 0.01))
      return (2589 + 3776*r + 4278*r2 + 1984*cube(r) + 363*r4)/1470.;

   return (2*(5*r4+25*r2+6))/(7*sqr(r2-1))
      + (2*(r4-19*r2-18)*r2*std::log(r2))/(7*cube(r2-1));
}

/// f5(r1,r2) in the limit r1 -> 1 and r2 -> 1
static double f5_1_1(double r1, double r2) noexcept
{
   return 0.772943722943723
      - 0.5524891774891774*r2
      + 0.7870670995670994*sqr(r2)
      - 0.3316558441558441*cube(r2)
      + 0.056277056277056266*quad(r2)
      + r1*(-0.5524891774891774
            + 1.0700757575757573*r2
            - 0.6625541125541123*sqr(r2)
            + 0.22483766233766228*cube(r2)
            - 0.03344155844155843*quad(r2))
      + cube(r1)*(-0.33165584415584404
                        + 0.22483766233766223*r2
                        - 0.08755411255411245*sqr(r2)
                        + 0.01650432900432896*cube(r2)
                        - 0.0007034632034631958*quad(r2))
      + quad(r1)*(0.05627705627705626
                        - 0.03344155844155841*r2
                        + 0.010281385281385256*sqr(r2)
                        - 0.0007034632034631921*cube(r2)
                        - 0.0002705627705627725*quad(r2))
      + sqr(r1)*(0.7870670995670994 - 0.6625541125541123*r2
                 + 0.32061688311688297*sqr(r2)
                 - 0.08755411255411248*cube(r2)
                 + 0.01028138528138527*quad(r2));
}

/// f5(r1,r2) in the limit r1 -> 1
static double f5_1_r2(double r1, double r2) noexcept
{
   const double lr22 = std::log(sqr(r2));

   return (-0.025*cube(-1. + r1)*(
              4. - 17.*r2 + 4.*sqr(r2)
              - 25.*cube(r2)
              - 20.*quad(r2)
              + 41.*pow5(r2)
              + 12.*pow6(r2)
              + pow7(r2)
              + (-30.*cube(r2) - 30.*pow5(r2))
              *lr22))/(pow6(-1. + r2)*sqr(1. + r2))
      - (0.125*sqr(-1. + r1)*(
            1. - 4.*r2 + sqr(r2)
            - 4.*cube(r2)
            - 5.*quad(r2)
            + 8.*pow5(r2)
            + 3.*pow6(r2)
            + (-6.*cube(r2) - 6.*pow5(r2))
            *lr22))/(pow5(-1. + r2)*sqr(1. + r2))
      + (0.75*(-1 + r2 + 2*sqr(r2)
               - quad(r2)
               - pow5(r2)
               + (cube(r2) + pow5(r2))
               *lr22))/(cube(-1 + r2)*sqr(1 + r2))
      + (0.25*(-1. + r1)*(
            1. - r2 - 2.*sqr(r2) + 8.*cube(r2)
            + quad(r2) - 7.*pow5(r2)
            + (3.*cube(r2) + 3.*pow5(r2))
            *lr22))/(quad(-1. + r2)*sqr(1. + r2))
      + (0.05*quad(-1. + r1)*(
            -1. + 4.5*r2 + 2.*sqr(r2)
            + 16.5*cube(r2) - 16.5*pow5(r2) - 2.*pow6(r2)
            - 4.5*pow7(r2) + pow8(r2)
            + (15.*cube(r2) + 15.*pow5(r2))
            *lr22))/(pow7(-1. + r2)*sqr(1. + r2));
}

/// f5(r1,r2) in the limit r1 -> 0
static double f5_0_r2(double r1, double r2) noexcept
{
   const double r22 = sqr(r2);
   const double lr22 = std::log(r22);

   return ((1 + r22)*(1 - r22 + r22*lr22))/sqr(-1 + r22) +
      (r1*r2*(2 - 2*r22 + lr22 +
              r22*lr22))/sqr(-1 + r22) +
      sqr(r1)*(-2/(-1 + r22) + ((1 + r22)*lr22)/sqr(-1 + r22));
}

/// f5(r1,r2) in the limit r1 -> 0 and r2 -> 1
static double f5_0_1(double, double r2) noexcept
{
   return 0.75*(1 + (-1 + r2)/3. + sqr(-1 + r2)/6.);
}

/// f5(r1,r2) in the limit r1 -> r2
static double f5_r1_r2(double r1, double r2) noexcept
{
   const double r22 = sqr(r2);
   const double lr22 = std::log(r22);

   return ((r1 - r2)*(11*r2 + 3*cube(r2) - 15*pow5(r2) + pow7(r2) +
        3*r2*lr22 + 18*cube(r2)*lr22 +
        3*pow5(r2)*lr22))/quad(-1 + r22)\
    + (sqr(r1 - r2)*(-17 - 116*r22 + 90*quad(r2) +
        44*pow6(r2) - pow8(r2) - 3*lr22 -
        75*r22*lr22 -
        105*quad(r2)*lr22 -
        9*pow6(r2)*lr22))/
    (3.*pow5(-1 + r22)) +
   (-1 - 5*r22 + 5*quad(r2) + pow6(r2) -
      3*r22*lr22 -
      6*quad(r2)*lr22 + pow6(r2)*lr22)
     /cube(-1 + r22) +
   (cube(r1 - r2)*(3 + 273*r22 + 314*quad(r2) -
        498*pow6(r2) - 93*pow8(r2) + power10(r2) +
        90*r22*lr22 +
        510*quad(r2)*lr22 +
        342*pow6(r2)*lr22 +
        18*pow8(r2)*lr22))/
    (6.*r2*pow6(-1 + r22));
}

double f5(double r1, double r2) noexcept
{
   if (is_equal(r1, 0., 0.0001) && is_equal(r2, 0., 0.0001))
      return 0.75;

   if (is_equal(r1, 1., 0.01) && is_equal(r2, 1., 0.01))
      return f5_1_1(r1, r2);

   if (is_equal(r1, -1., 0.01) && is_equal(r2, -1., 0.01))
      return f5_1_1(-r1, -r2);

   if (is_equal(r1, 1., 0.01)) {
      if (is_equal(r2, 0., 0.01))
         return f5_0_1(r2, r1);

      return f5_1_r2(r1, r2);
   }

   if (is_equal(r2, 1., 0.01)) {
      if (is_equal(r1, 0., 0.01))
         return f5_0_1(r1, r2);

      return f5_1_r2(r2, r1);
   }

   if (is_equal(r1, 0., 0.0001))
      return 0.75 * f5_0_r2(r1, r2);

   if (is_equal(r2, 0., 0.0001))
      return 0.75 * f5_0_r2(r2, r1);

   if (is_equal(r1, r2, 0.0001))
      return 0.75 * f5_r1_r2(r2, r1);

   const double r12 = sqr(r1);
   const double r22 = sqr(r2);

   const double result
      = (1+sqr(r1+r2)-r12*r22)/((r12-1)*(r22-1))
      + (cube(r1)*(r12+1)*std::log(r12))/(sqr(r12-1)*(r1-r2))
      - (cube(r2)*(r22+1)*std::log(r22))/((r1-r2)*sqr(r22-1));

   return 0.75 * result;
}

/// f6(r1,r2) in the limit r1 -> 1 and r2 -> 1
static double f6_1_1(double r1, double r2) noexcept
{
   return 1 + (4*(-1 + r1))/7. - (2*sqr(-1 + r1))/35.
      - cube(-1 + r1)/70. + (9*quad(-1 + r1))/490.
      + (0.5714285714285714 - (2*(-1 + r1))/35. - sqr(-1 + r1)/70.
         + (9*cube(-1 + r1))/490.
         - (3*quad(-1 + r1))/245.)*(-1 + r2)
      + (-0.05714285714285714 + (1 - r1)/70. + (9*sqr(-1 + r1))/490.
         - (3*cube(-1 + r1))/245.
         + quad(-1 + r1)/147.)*sqr(-1 + r2)
      + (-0.014285714285714285 + (9*(-1 + r1))/490.
         - (3*sqr(-1 + r1))/245. + cube(-1 + r1)/147.
         - quad(-1 + r1)/294.)*cube(-1 + r2)
      + (0.018367346938775512 - (3*(-1 + r1))/245. + sqr(-1 + r1)/147.
         - cube(-1 + r1)/294.
         + (5*quad(-1 + r1))/3234.)*quad(-1 + r2);
}

/// f6(r1,r2) in the limit r1 -> 1
static double f6_1_r2(double r1, double r2) noexcept
{
   const double r22 = sqr(r2);

   return ((-1 + r22)*(
              -3 + 16*r2 + 33*r22 - 332*cube(r2)
              + 573*quad(r2) - 584*pow5(r2) + 297*pow6(r2)
              - 60*pow7(r2) - 2*cube(r1)*(
                 8 - 49*r2 + 121*r22 - 129*cube(r2)
                 - 99*quad(r2) + 28*pow5(r2))
              + quad(r1)*(
                 3 - 18*r2 + 43*r22 - 42*cube(r2)
                 - 47*quad(r2) + pow6(r2))
              - 2*r1*(-8 + 10*r2 + 106*r22 - 359*cube(r2)
                      + 211*quad(r2) - 101*pow5(r2)
                      + 21*pow6(r2))
              - 2*sqr(r1)*(-15 + 98*r2 - 264*r22
                                  + 331*cube(r2) + 106*quad(r2)
                                  - 99*pow5(r2) + 23*pow6(r2)))
           + 60*pow5(r2)*(5 + quad(r1) + cube(r1)*(-5 + r2)
                                - 10*r2 + 10*r22 - 5*cube(r2)
                                + quad(r2) + sqr(r1)*(
                                   10 - 5*r2 + r22)
                                + r1*(-10 + 10*r2 - 5*r22
                                      + cube(r2)))*std::log(r22))
      /(70.*pow7(-1 + r2)*sqr(1 + r2));
}

/// f6(r1,r2) in the limit r1 -> 0
static double f6_0_r2(double r1, double r2) noexcept
{
   const double r22 = sqr(r2);
   const double lr22 = std::log(r22);

   return (sqr(r1)*(1 - r22 + r22*lr22))/sqr(-1 + r22) +
      (r1*(r2 - cube(r2) + cube(r2)*lr22))/sqr(-1 + r22) +
      (r22 - quad(r2) + quad(r2)*lr22)/sqr(-1 + r22);
}

// f6(r1,r2) in the limit r1 -> r2
static double f6_r1_r2(double r1, double r2) noexcept
{
   const double r22 = sqr(r2);
   const double lr22 = std::log(r22);

   return ((r1 - r2)*(3*r2 + 7*cube(r2) - 11*pow5(r2) + pow7(r2) +
        10*cube(r2)*lr22 +
        2*pow5(r2)*lr22))/quad(-1 + r22)\
    + (sqr(r1 - r2)*(-3 - 62*r22 + 36*quad(r2) +
        30*pow6(r2) - pow8(r2) -
        30*r22*lr22 -
        60*quad(r2)*lr22 -
        6*pow6(r2)*lr22))/
    (3.*pow5(-1 + r22)) +
   (-3*r22 + 2*quad(r2) + pow6(r2) -
      5*quad(r2)*lr22 + pow6(r2)*lr22)
     /cube(-1 + r22) +
   (cube(r1 - r2)*(107*r2 + 206*cube(r2) - 252*pow5(r2) -
        62*pow7(r2) + pow9(r2) + 30*r2*lr22 +
        240*cube(r2)*lr22 +
        198*pow5(r2)*lr22 +
        12*pow7(r2)*lr22))/
    (6.*pow6(-1 + r22));
}

/// f6(r1,r2) in the limit r1 -> 0 and r2 -> 1
static double f6_0_1(double, double r2) noexcept
{
   return 6./7.*(0.5 + (2*(-1 + r2))/3. - cube(-1 + r2)/15.
                 + quad(-1 + r2)/20.);
}

double f6(double r1, double r2) noexcept
{
   if (is_equal(r1, 0., 0.0001) && is_equal(r2, 0., 0.0001))
      return 0.;

   if (is_equal(r1, 1., 0.01) && is_equal(r2, 1., 0.01))
      return f6_1_1(r1, r2);

   if (is_equal(r1, -1., 0.01) && is_equal(r2, -1., 0.01))
      return f6_1_1(-r1, -r2);

   if (is_equal(r1, 1., 0.01)) {
      if (is_equal(r2, 0., 0.0001))
         return f6_0_1(r2, r1);

      return f6_1_r2(r1, r2);
   }

   if (is_equal(r2, 1., 0.01)) {
      if (is_equal(r1, 0., 0.0001))
         return f6_0_1(r1, r2);

      return f6_1_r2(r2, r1);
   }

   if (is_equal(r1, 0., 0.0001))
      return 6./7. * f6_0_r2(r1, r2);

   if (is_equal(r2, 0., 0.0001))
      return 6./7. * f6_0_r2(r2, r1);

   if (is_equal(r1, r2, 0.0001))
      return 6./7. * f6_r1_r2(r2, r1);

   const double r12 = sqr(r1);
   const double r22 = sqr(r2);

   const double result
      = (r12+r22+r1*r2-r12*r22)/((r12-1)*(r22-1))
      + (pow5(r1)*std::log(r12))/(sqr(r12-1)*(r1-r2))
      - (pow5(r2)*std::log(r22))/((r1-r2)*sqr(r22-1));

   return 6./7. * result;
}

/// f7(r1,r2) in the limit r1 -> 1 and r2 -> 1
static double f7_1_1(double r1, double r2) noexcept
{
   const double r22 = sqr(r2);

   return (15700 - 14411*r2 + 7850*r22 - 2498*cube(r2)
           + 355*quad(r2)
           + sqr(r1)*(7850 - 2558*r2 - 750*r22 + 940*cube(r2)
                      - 235*quad(r2))
           + quad(r1)*(355 + 65*r2 - 235*r22 + 142*cube(r2)
                             - 30*quad(r2))
           + r1*(-14411 + 8375*r2 - 2558*r22 + 180*cube(r2)
                 + 65*quad(r2))
           + cube(r1)*(-2498 + 180*r2 + 940*r22
                             - 645*cube(r2) + 142*quad(r2)))
      /2310.;
}

/// f7(r1,r2) in the limit r1 -> 1
static double f7_1_r2(double r1, double r2) noexcept
{
   const double r22 = sqr(r2);
   const double lr22 = std::log(r22);

   return (-10*(-1 + r1)*cube(-1 + r2)*(
              2 - 5*r2 - 4*r22 + 4*cube(r2) + 2*quad(r2)
              + pow5(r2) - 6*cube(r2)*lr22)
           - 30*quad(-1 + r2)*(
              1 - 2*r2 - 2*r22 + 2*cube(r2) + quad(r2)
              - 2*cube(r2)*lr22)
           + 10*sqr(-1 + r1)*sqr(-1 + r2)*(
              -1 + 3*r2 + 3*r22 - 3*quad(r2) - 3*pow5(r2)
              + pow6(r2) + 6*cube(r2)*lr22)
           + 2*cube(-1 + r1)*(-1 + r2)*(
              -2 + 6*r2 + 18*r22 + 15*cube(r2) - 30*quad(r2)
              - 18*pow5(r2) + 14*pow6(r2) - 3*pow7(r2)
              + 30*cube(r2)*lr22)
           + quad(-1 + r1)*(
              -1 + 48*r22 + 42*cube(r2) - 90*quad(r2)
              - 24*pow5(r2) + 40*pow6(r2) - 18*pow7(r2)
              + 3*pow8(r2) + 60*cube(r2)*lr22))
      /(10.*pow7(-1 + r2)*sqr(1 + r2));
}

/// f7(r1,r2) in the limit r1 -> 0
static double f7_0_r2(double r1, double r2) noexcept
{
   const double r22 = sqr(r2);
   const double lr22 = std::log(r22);

   return -((r1*r2*(-1 + r22 - lr22))/sqr(-1 + r22)) +
      (sqr(r1)*(1 - r22 + lr22))/sqr(-1 + r22) +
      (1 - r22 + r22*lr22)/sqr(-1 + r22);
}

/// f7(r1,r2) in the limit r1 -> 0 and r2 -> 1
static double f7_0_1(double, double r2) noexcept
{
   return 6.*(0.5 + (1 - r2)/3. + sqr(-1 + r2)/6. - cube(-1 + r2)/15.
              + quad(-1 + r2)/60.);
}

/// f7(r1,r2) in the limit r1 -> r2
static double f7_r1_r2(double r1, double r2) noexcept
{
   const double r22 = sqr(r2);
   const double lr22 = std::log(r22);

   return (-1 - 2*r22 + 3*quad(r2) -
      3*r22*lr22 - quad(r2)*lr22)
     /cube(-1 + r22) +
   ((r1 - r2)*(8*r2 - 4*cube(r2) - 4*pow5(r2) +
        3*r2*lr22 + 8*cube(r2)*lr22 +
        pow5(r2)*lr22))/quad(-1 + r22) +
   (sqr(r1 - r2)*(-14 - 54*r22 + 54*quad(r2) +
        14*pow6(r2) - 3*lr22 -
        45*r22*lr22 -
        45*quad(r2)*lr22 -
        3*pow6(r2)*lr22))/
    (3.*pow5(-1 + r22)) +
   (cube(r1 - r2)*(3 + 166*r22 + 108*quad(r2) -
        246*pow6(r2) - 31*pow8(r2) +
        60*r22*lr22 +
        270*quad(r2)*lr22 +
        144*pow6(r2)*lr22 +
        6*pow8(r2)*lr22))/
    (6.*r2*pow6(-1 + r22));
}

double f7(double r1, double r2) noexcept
{
   if (is_equal(r1, 0., 0.0001) && is_equal(r2, 0., 0.0001))
      return 6.;

   if (is_equal(r1, 1., 0.01) && is_equal(r2, 1., 0.01))
      return f7_1_1(r1, r2);

   if (is_equal(r1, -1., 0.01) && is_equal(r2, -1., 0.01))
      return f7_1_1(-r1, -r2);

   if (is_equal(r1, 1., 0.01)) {
      if (is_equal(r2, 0., 0.0001))
         return f7_0_1(r2, r1);

      return f7_1_r2(r1, r2);
   }

   if (is_equal(r2, 1., 0.01)) {
      if (is_equal(r1, 0., 0.0001))
         return f7_0_1(r1, r2);

      return f7_1_r2(r2, r1);
   }

   if (is_equal(r1, 0., 0.0001))
      return 6. * f7_0_r2(r1, r2);

   if (is_equal(r2, 0., 0.0001))
      return 6. * f7_0_r2(r2, r1);

   if (is_equal(r1, r2, 0.0001))
      return 6. * f7_r1_r2(r2, r1);

   const double r12 = sqr(r1);
   const double r22 = sqr(r2);

   const double result
      = (1+r1*r2)/((r12-1)*(r22-1))
      + (cube(r1)*std::log(r12))/(sqr(r12-1)*(r1-r2))
      - (cube(r2)*std::log(r22))/((r1-r2)*sqr(r22-1));

   return 6. * result;
}

/// f8(r1,r2) in the limit r1 -> 1 and r2 -> 1
static double f8_1_1(double r1, double r2) noexcept
{
   return 1 - sqr(-1 + r1)/10. + (3*cube(-1 + r1))/40.
      - (3*quad(-1 + r1))/70.
      + ((1 - r1)/10. + (3*sqr(-1 + r1))/40.
         - (3*cube(-1 + r1))/70.
         + (3*quad(-1 + r1))/140.)*(-1 + r2)
      + (-0.1 + (3*(-1 + r1))/40. - (3*sqr(-1 + r1))/70.
         + (3*cube(-1 + r1))/140.
         - quad(-1 + r1)/105.)*sqr(-1 + r2)
      + (0.075 - (3*(-1 + r1))/70. + (3*sqr(-1 + r1))/140.
         - cube(-1 + r1)/105.
         + quad(-1 + r1)/280.)*cube(-1 + r2)
      + (-0.04285714285714286 + (3*(-1 + r1))/140.
         - sqr(-1 + r1)/105. + cube(-1 + r1)/280.
         - quad(-1 + r1)/1155.)*quad(-1 + r2);
}

/// f8(r1,r2) in the limit r1 -> 1
static double f8_1_r2(double r1, double r2) noexcept
{
   const double r22 = sqr(r2);
   const double lr22 = std::log(r22);

   return (30*quad(-1 + r2)*(
              -1 + 4*r22 - 3*quad(r2)
              + 2*quad(r2)*lr22)
           + 10*sqr(-1 + r1)*sqr(-1 + r2)*(
              1 - 4*r2 + 4*r22 + 8*cube(r2)
              - 5*quad(r2) - 4*pow5(r2)
              + 6*quad(r2)*lr22)
           + 10*(-1 + r1)*cube(-1 + r2)*(
              1 - 4*r2 + 4*r22 + 8*cube(r2)
              - 5*quad(r2) - 4*pow5(r2)
              + 6*quad(r2)*lr22)
           + 2*cube(-1 + r1)*(-1 + r2)*(
              3 - 14*r2 + 18*r22 + 30*cube(r2)
              - 15*quad(r2) - 18*pow5(r2) - 6*pow6(r2)
              + 2*pow7(r2) + 30*quad(r2)*lr22)
           + quad(-1 + r1)*(
              3 - 16*r2 + 24*r22 + 48*cube(r2)
              - 48*pow5(r2) - 24*pow6(r2) + 16*pow7(r2)
              - 3*pow8(r2) + 60*quad(r2)*lr22))
      /(40.*pow7(-1 + r2)*sqr(1 + r2));
}

/// f8(r1,r2) in the limit r1 -> 0
static double f8_0_r2(double r1, double r2) noexcept
{
   const double r22 = sqr(r2);
   const double lr22 = std::log(r22);

   return -((sqr(r1)*r2*(-1 + r22 - lr22))/sqr(-1 + r22)) +
      (r1*(1 - r22 + r22*lr22))/sqr(-1 + r22) +
      (r2 - cube(r2) + cube(r2)*lr22)/sqr(-1 + r22);
}

/// f8(r1,r2) in the limit r1 -> 0 and r2 -> 1
static double f8_0_1(double, double r2) noexcept
{
   return 1.5*(0.5 + (-1 + r2)/6. - sqr(-1 + r2)/6. + cube(-1 + r2)/10.
               - quad(-1 + r2)/20.);
}

/// f8(r1,r2) in the limit r1 -> r2
static double f8_r1_r2(double r1, double r2) noexcept
{
   const double r22 = sqr(r2);
   const double lr22 = std::log(r22);

   return (2*(-r2 + pow5(r2) - 2*cube(r2)*lr22))/
    cube(-1 + r22) +
   ((r1 - r2)*(1 + 9*r22 - 9*quad(r2) - pow6(r2) +
        6*r22*lr22 +
        6*quad(r2)*lr22))/quad(-1 + r22)\
    + (2*sqr(r1 - r2)*(-19*r2 - 9*cube(r2) +
        27*pow5(r2) + pow7(r2) - 6*r2*lr22 -
        30*cube(r2)*lr22 -
        12*pow5(r2)*lr22))/
    (3.*pow5(-1 + r22)) +
   (cube(r1 - r2)*(31 + 246*r22 - 108*quad(r2) -
        166*pow6(r2) - 3*pow8(r2) + 6*lr22 +
        144*r22*lr22 +
        270*quad(r2)*lr22 +
        60*pow6(r2)*lr22))/
    (6.*pow6(-1 + r22));
}

double f8(double r1, double r2) noexcept
{
   if (is_equal(r1, 0., 0.0001) && is_equal(r2, 0., 0.0001))
      return 0.;

   if (is_equal(r1, 1., 0.01) && is_equal(r2, 1., 0.01))
      return f8_1_1(r1, r2);

   if (is_equal(r1, -1., 0.01) && is_equal(r2, -1., 0.01))
      return -1.;

   if (is_equal(r1, 1., 0.01)) {
      if (is_equal(r2, 0., 0.0001))
         return f8_0_1(r2, r1);

      return f8_1_r2(r1, r2);
   }

   if (is_equal(r2, 1., 0.01)) {
      if (is_equal(r1, 0., 0.0001))
         return f8_0_1(r1, r2);

      return f8_1_r2(r2, r1);
   }

   if (is_equal(r1, 0., 0.0001))
      return 1.5 * f8_0_r2(r1, r2);

   if (is_equal(r2, 0., 0.0001))
      return 1.5 * f8_0_r2(r2, r1);

   if (is_equal(r1, r2, 0.0001))
      return 1.5 * f8_r1_r2(r2, r1);

   const double r12 = sqr(r1);
   const double r22 = sqr(r2);

   const double result
      = (r1+r2)/((r12-1)*(r22-1))
      + (quad(r1)*std::log(r12))/(sqr(r12-1)*(r1-r2))
      - (quad(r2)*std::log(r22))/((r1-r2)*sqr(r22-1));

   return 1.5 * result;
}

double fth1(double y) noexcept
{
   using std::log;

   if (is_zero(y))
      return 0.;

   if (is_equal(std::abs(y), 1.))
      return -1.;

   if (!is_zero(y) && !is_equal(std::abs(y), 1.)) {
      const double y2 = sqr(y);

      return y2*log(y2) / (1. - y2);
   }

   return 0.;
}

double fth2(double y) noexcept
{
   using std::log;

   if (is_zero(y))
      return 0.;

   if (is_equal(std::abs(y), 1.))
      return 0.5;

   if (!is_zero(y) && !is_equal(std::abs(y), 1.)) {
      const double y2 = sqr(y);

      return (1. + y2*log(y2) / (1. - y2)) / (1 - y2);
   }

   return 0.;
}

double fth3(double y) noexcept
{
   using std::log;
   const double z2 = sqr(Pi) / 6.;

   if (is_zero(y))
      return z2;

   if (is_equal(std::abs(y), 1.))
      return -9./4.;

   if (!is_zero(y) && !is_equal(std::abs(y), 1.)) {
      const double y2 = sqr(y);
      const double y4 = sqr(y2);

      return (-1. + 2*y2 + 2*y4)
         *(-z2 - y2*log(y2) + log(y2)*log(1. - y2) + dilog(y2))
         / sqr(1 - y2);
   }

   return 0.;
}

/// Iabc(a,a,a)
static double Iaaa(double a, double b, double c) noexcept
{
   return (151.*quad(a) + 13.*sqr(b)*sqr(c) - 128.*cube(a)*(b + c) - 40.*a*b*c*(b + c)
           + sqr(a)*(37.*sqr(b) + 128.*b*c + 37.*sqr(c))) / (60.*pow6(a));
}

/// Iabc(a,a,c)
static double Iaac(double a, double b, double c) noexcept
{
   return ((sqr(a) - sqr(c))
           * (17.*pow6(a) - 16.*pow5(a)*b - 40.*cube(a)*b*sqr(c)
              + 8.*a*b*quad(c) - sqr(b)*quad(c) + quad(a)*(5.*sqr(b) + 8.*sqr(c))
              + sqr(a)*(20.*sqr(b)*sqr(c) - quad(c)))
           - 6.*sqr(a)*sqr(c) * log(sqr(a)/sqr(c))
           * (6.*quad(a) - 8.*cube(a)*b + 3.*sqr(a)*(sqr(b) - sqr(c)) + sqr(c)*(sqr(b) + sqr(c))))
      / (6.*sqr(a)*quad(sqr(a) - sqr(c)));
}

/// Iabc(a,a,0)
double Iaa0(double a, double b) noexcept
{
   return (17.*sqr(a) - 16.*a*b + 5.*sqr(b)) / (6.*quad(a));
}

/// Iabc(0,b,c)
double I0bc(double b, double c)
{
   return log(sqr(b/c))/(sqr(b) - sqr(c));
}

double Iabc(double a, double b, double c) noexcept {
   if ((is_zero(a) && is_zero(b) && is_zero(c)) ||
       (is_zero(a) && is_zero(b)) ||
       (is_zero(a) && is_zero(c)) ||
       (is_zero(b) && is_zero(c)))
      return 0.;

   if (is_equal_rel(std::abs(a), std::abs(b), 0.01) && is_equal_rel(std::abs(a), std::abs(c), 0.01))
      return Iaaa(std::abs(a),std::abs(b),std::abs(c));

   if (is_equal_rel(std::abs(a), std::abs(b), 0.01)) {
      if (is_zero(c))
         return Iaa0(std::abs(a),std::abs(b));
      return Iaac(std::abs(a),std::abs(b),c);
   }

   if (is_equal_rel(std::abs(b), std::abs(c), 0.01)) {
      if (is_zero(a))
         return Iaa0(std::abs(b),std::abs(c));
      return Iaac(std::abs(b),std::abs(c),a);
   }

   if (is_equal_rel(std::abs(a), std::abs(c), 0.01)) {
      if (is_zero(b))
         return Iaa0(std::abs(a),std::abs(c));
      return Iaac(std::abs(a),std::abs(c),b);
   }

   if (is_zero(a))
      return I0bc(b,c);

   if (is_zero(b))
      return I0bc(c,a);

   if (is_zero(c))
      return I0bc(a,b);

   return ( (sqr(a * b) * log(sqr(a / b))
           + sqr(b * c) * log(sqr(b / c))
           + sqr(c * a) * log(sqr(c / a)))
           / ((sqr(a) - sqr(b)) * (sqr(b) - sqr(c)) * (sqr(a) - sqr(c))) );
}

/// Delta function from hep-ph/0907.47682v1
double delta_xyz(double x, double y, double z) noexcept
{
   return sqr(x)+sqr(y)+sqr(z)-2*(x*y+x*z+y*z);
}

namespace {
   /// lambda^2(u,v)
   double lambda_2(double u, double v) noexcept
   {
      return sqr(1 - u - v) - 4*u*v;
   }

   /// u < 1 && v < 1, lambda^2(u,v) > 0
   double phi_pos(double u, double v) noexcept
   {
      using std::log;
      const auto lambda = std::sqrt(lambda_2(u,v));

      return (-(log(u)*log(v))
              + 2*log((1 - lambda + u - v)/2.)*log((1 - lambda - u + v)/2.)
              - 2*dilog((1 - lambda + u - v)/2.)
              - 2*dilog((1 - lambda - u + v)/2.)
              + sqr(Pi)/3.)/lambda;
   }

   /// lambda^2(u,v) < 0
   double phi_neg(double u, double v) noexcept
   {
      using std::acos;
      using std::sqrt;
      const auto lambda = std::sqrt(-lambda_2(u,v));

      return 2*(+ clausen_2(2*acos((1 + u - v)/(2.*sqrt(u))))
                + clausen_2(2*acos((1 - u + v)/(2.*sqrt(v))))
                + clausen_2(2*acos((-1 + u + v)/(2.*sqrt(u*v)))))/lambda;
   }

   /**
    * Phi(u,v) with u = x/z, v = y/z.
    *
    * The following identities hold:
    * Phi(u,v) = Phi(v,u) = Phi(1/u,v/u)/u = Phi(1/v,u/v)/v
    */
   double phi_uv(double u, double v) noexcept
   {
      const auto lambda = lambda_2(u,v);

      if (lambda > 0.) {
         if (u <= 1 && v <= 1)
            return phi_pos(u,v);
         if (u >= 1 && v/u <= 1)
            return phi_pos(1./u,v/u)/u;
         // v >= 1 && u/v <= 1
         return phi_pos(1./v,u/v)/v;
      }

      return phi_neg(u,v);
   }
} // anonymous namespace

/**
 * \f$\Phi(x,y,z)\f$ function.  The arguments x, y and z are
 * interpreted as squared masses.
 *
 * Davydychev and Tausk, Nucl. Phys. B397 (1993) 23
 *
 * @param x squared mass
 * @param y squared mass
 * @param z squared mass
 *
 * @return \f$\Phi(x,y,z)\f$
 */
double phi_xyz(double x, double y, double z) noexcept
{
   const auto u = x/z, v = y/z;
   return phi_uv(u,v);
}

} // namespace threshold_loop_functions
} // namespace himalaya
