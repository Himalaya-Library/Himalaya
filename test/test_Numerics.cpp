#include "doctest.h"
#include "himalaya/misc/Numerics.hpp"
#include <limits>


TEST_CASE("test_dabs")
{
   using namespace himalaya;

   CHECK(dabs(   0) ==   0);
   CHECK(dabs( 0.0) == 0.0);
   CHECK(dabs(-0.0) == 0.0);

   CHECK(dabs(   1) ==   1);
   CHECK(dabs(  -1) ==   1);
   CHECK(dabs( 1.0) == 1.0);
   CHECK(dabs(-1.0) == 1.0);
}


TEST_CASE("test_is_zero")
{
   using namespace himalaya;

   const double deps = std::numeric_limits<double>::epsilon();

   CHECK( is_zero(0));
   CHECK( is_zero(0.0));
   CHECK( is_zero(deps, deps));
   CHECK(!is_zero(2*deps, deps));
   CHECK(!is_zero(1.0, deps));
}


TEST_CASE("test_is_equal")
{
   using namespace himalaya;

   const double deps = std::numeric_limits<double>::epsilon();

   CHECK( is_equal(0, 0, 0));
   CHECK( is_equal(0.0, 0.0, deps));
   CHECK( is_equal(deps, deps));
   CHECK( is_equal(2*deps, deps, deps));
   CHECK(!is_equal(3*deps, deps, deps));
   CHECK( is_equal(1.0, 2.0, 1.0));
   CHECK(!is_equal(1.0, 3.0, 1.0));
}


TEST_CASE("test_is_equal_rel")
{
   using namespace himalaya;

   const double deps = std::numeric_limits<double>::epsilon();
   const double small = deps/10;

   CHECK( is_equal_rel(0, 0, 0));
   CHECK( is_equal_rel(0.0, 0.0, deps));
   CHECK( is_equal_rel(0.0, deps, deps));
   CHECK( is_equal_rel(deps, deps));
   CHECK( is_equal_rel(2*deps, deps, deps));
   CHECK(!is_equal_rel(3*deps, deps, deps));
   CHECK( is_equal_rel(small, deps, deps));
   CHECK( is_equal_rel(small, 2*small, deps));
   CHECK( is_equal_rel(1.0, 1.0 + deps, deps));
   CHECK(!is_equal_rel(1.0, 1.0 + 2*deps, deps));
   CHECK( is_equal_rel(1.0, 2.0, 0.5));
   CHECK( is_equal_rel(1.0, 2.0 + deps, 0.5));
   CHECK(!is_equal_rel(1.0, 2.0 + 10*deps, 0.5));
   CHECK( is_equal_rel(1.0, 3.0, 2./3));
   CHECK(!is_equal_rel(1.0, 3.0 + 10*deps, 2./3));
}


TEST_CASE("test_is_close")
{
   using namespace himalaya;

   const double deps = std::numeric_limits<double>::epsilon();
   const double small = deps/10;

   CHECK( is_close(0, 0, 0));
   CHECK( is_close(0.0, 0.0, deps));
   CHECK( is_close(0.0, deps, deps));
   CHECK( is_close(deps, deps));
   CHECK( is_close(2*deps, deps, deps));
   CHECK(!is_close(3*deps, deps, deps));
   CHECK( is_close(small, deps, deps));
   CHECK( is_close(small, 2*small, deps));
   CHECK( is_close(1.0, 1.0 + deps, deps));
   CHECK( is_close(1.0, 1.0 + 2*deps, deps));
   CHECK( is_close(1.0, 2.0, 0.5));
   CHECK( is_close(1.0, 2.0 + deps, 0.5));
   CHECK( is_close(1.0, 2.0 + 10*deps, 0.5));
   CHECK( is_close(1.0, 3.0, 2./3));
   CHECK( is_close(1.0, 3.0 + 10*deps, 2./3));
}
