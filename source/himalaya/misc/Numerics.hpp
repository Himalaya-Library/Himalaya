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

/// compares a number for being close to zero
template <typename T>
typename std::enable_if<!std::is_unsigned<T>::value, bool>::type
is_zero(T a, T prec = std::numeric_limits<T>::epsilon()) noexcept
{
   return std::abs(a) <= prec;
}

/// compares two numbers for absolute equality
template <typename T>
bool is_equal(T a, T b, T prec = std::numeric_limits<T>::epsilon()) noexcept
{
   return is_zero(a - b, prec);
}

/// compares two numbers for relative equality, treating numbers with
/// small differences as equal
template <typename T>
bool is_equal_rel(T a, T b, T prec = std::numeric_limits<T>::epsilon()) noexcept
{
   if (is_equal(a, b, std::numeric_limits<T>::epsilon()))
      return true;

   const T min = std::min(std::abs(a), std::abs(b));

   if (min < std::numeric_limits<T>::epsilon()) {
      return is_equal(a, b, prec);
   }

   const T max = std::max(std::abs(a), std::abs(b));

   return is_equal(a, b, prec*max);
}

/// compares two numbers for relative equality
template <typename T>
bool is_close(T a, T b, T prec = std::numeric_limits<T>::epsilon()) noexcept
{
   const T max = std::max(std::abs(a), std::abs(b));
   return is_zero(a - b, prec*(1.0 + max));
}

} // anonymous namespace
} // namespace himalaya
