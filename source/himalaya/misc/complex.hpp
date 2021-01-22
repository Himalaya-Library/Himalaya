// ====================================================================
// This file is part of Himalaya.
//
// Himalaya is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================

#pragma once

#include <cmath>
#include <complex>

namespace himalaya {

template <typename T>
struct Complex {
   constexpr Complex(T re_ = T{}, T im_ = T{}) : re(re_), im(im_) {}
   operator std::complex<T>() const noexcept { return std::complex<T>(re, im); }
   T re{};
   T im{};
};

/// converts -0.0 to 0.0
template <typename T>
Complex<T> pos(const Complex<T>& z) noexcept
{
   const T rz = z.re == T(0) ? std::abs(z.re) : z.re;
   const T iz = z.im == T(0) ? std::abs(z.im) : z.im;
   return { rz, iz };
}

template <typename T>
constexpr T arg(const Complex<T>& z) noexcept
{
   return std::atan2(z.im, z.re);
}

template <typename T>
constexpr Complex<T> conj(const Complex<T>& z) noexcept
{
   return { z.re, -z.im };
}

template <typename T>
Complex<T> log(const Complex<T>& z_) noexcept
{
   const Complex<T> z = pos(z_);
   return { 0.5*std::log(norm_sqr(z)), arg(z) };
}

template <typename T>
constexpr T norm_sqr(const Complex<T>& z) noexcept
{
   return z.re*z.re + z.im*z.im;
}

template <typename T>
constexpr Complex<T> operator+(const Complex<T>& a, const Complex<T>& b) noexcept
{
   return { a.re + b.re, a.im + b.im };
}

template <typename T>
constexpr Complex<T> operator+(const Complex<T>& z, T x) noexcept
{
   return { z.re + x, z.im };
}

template <typename T>
constexpr Complex<T> operator+(T x, const Complex<T>& z) noexcept
{
   return { x + z.re, z.im };
}

template <typename T>
constexpr Complex<T> operator-(const Complex<T>& a, const Complex<T>& b) noexcept
{
   return { a.re - b.re, a.im - b.im };
}

template <typename T>
constexpr Complex<T> operator-(T x, const Complex<T>& z) noexcept
{
   return { x - z.re, -z.im };
}

template <typename T>
constexpr Complex<T> operator-(const Complex<T>& z, T x) noexcept
{
   return { z.re - x, z.im };
}

template <typename T>
constexpr Complex<T> operator-(const Complex<T>& z) noexcept
{
   return { -z.re, -z.im };
}

template <typename T>
constexpr Complex<T> operator*(const Complex<T>& a, const Complex<T>& b) noexcept
{
   return { a.re*b.re - a.im*b.im, a.re*b.im + a.im*b.re };
}

template <typename T>
constexpr Complex<T> operator*(T x, const Complex<T>& z) noexcept
{
   return { x*z.re, x*z.im };
}

template <typename T>
constexpr Complex<T> operator*(const Complex<T>& z, T x) noexcept
{
   return x*z;
}

template <typename T>
constexpr Complex<T> operator/(T x, const Complex<T>& z) noexcept
{
   return x*conj(z)/norm_sqr(z);
}

template <typename T>
constexpr Complex<T> operator/(const Complex<T>& z, T x) noexcept
{
   return { z.re/x, z.im/x };
}

} // namespace himalaya
