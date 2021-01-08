// ====================================================================
// This file is part of Himalaya.
//
// Himalaya is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================

#pragma once

namespace himalaya {
namespace {

template <typename T> T pow2(T x)    noexcept { return x * x; }
template <typename T> T pow3(T x)    noexcept { return x * x * x; }
template <typename T> T pow4(T x)    noexcept { return pow2(pow2(x)); }
template <typename T> T pow5(T x)    noexcept { return pow4(x) * x; }
template <typename T> T pow6(T x)    noexcept { return pow2(pow2(x) * x); }
template <typename T> T pow7(T x)    noexcept { return pow6(x) * x; }
template <typename T> T pow8(T x)    noexcept { return pow2(pow4(x)); }
template <typename T> T pow9(T x)    noexcept { return pow8(x) * x; }
template <typename T> T power10(T x) noexcept { return pow2(pow5(x)); }
template <typename T> T pow11(T x)   noexcept { return power10(x) * x; }
template <typename T> T pow12(T x)   noexcept { return pow2(pow6(x)); }
template <typename T> T pow13(T x)   noexcept { return pow4(x) * pow9(x); }
template <typename T> T pow14(T x)   noexcept { return pow2(pow7(x)); }
template <typename T> T pow15(T x)   noexcept { return pow6(x) * pow9(x); }
template <typename T> T pow16(T x)   noexcept { return pow2(pow8(x)); }
template <typename T> T pow17(T x)   noexcept { return pow16(x) * x; }
template <typename T> T pow18(T x)   noexcept { return pow2(pow9(x)); }
template <typename T> T pow19(T x)   noexcept { return pow18(x) * x; }
template <typename T> T pow20(T x)   noexcept { return pow2(power10(x)); }

} // anonymous namespace
} // himalaya namespace
