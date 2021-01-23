// ====================================================================
// This file is part of Himalaya.
//
// Himalaya is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================

#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <type_traits>

namespace himalaya {
namespace {

/// absolute
template <typename T>
constexpr T dabs(T a) noexcept
{
   return a >= T{0} ? a : -a;
}

/// compares a number for being close to zero
template <typename T>
constexpr typename std::enable_if<!std::is_unsigned<T>::value, bool>::type
is_zero(T a, T prec = std::numeric_limits<T>::epsilon()) noexcept
{
   return dabs(a) <= prec;
}

/// compares two numbers for absolute equality
template <typename T>
constexpr bool is_equal(T a, T b, T prec = std::numeric_limits<T>::epsilon()) noexcept
{
   return is_zero(a - b, prec);
}

/// compares two numbers for relative equality, treating numbers with
/// small differences as equal
template <typename T>
constexpr bool is_equal_rel(T a, T b, T prec = std::numeric_limits<T>::epsilon()) noexcept
{
   if (is_equal(a, b, std::numeric_limits<T>::epsilon()))
      return true;

   const T min = std::min(dabs(a), dabs(b));

   if (min < std::numeric_limits<T>::epsilon()) {
      return is_equal(a, b, prec);
   }

   const T max = std::max(dabs(a), dabs(b));

   return is_equal(a, b, prec*max);
}

/// compares two numbers for relative equality
template <typename T>
constexpr bool is_close(T a, T b, T prec = std::numeric_limits<T>::epsilon()) noexcept
{
   const T max = std::max(dabs(a), dabs(b));
   return is_zero(a - b, prec*(1 + max));
}

} // anonymous namespace
} // namespace himalaya
