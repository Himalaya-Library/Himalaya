// ====================================================================
// This file is part of Himalaya.
//
// Himalaya is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================

#pragma once

/**
 * @file sum.hpp
 *
 * @brief Contains the definition of the \a SUM macro.
 *
 * Usage example of the \a SUM macro:
 *
 * Calculate \f$\sum_{n=1}^{10} n^2\f$
 *
 * @code{.cpp}
const double s = SUM(n,1,10,n*n);
   @endcode
 */

#include <type_traits>
#include <cstddef>
#include <Eigen/Core>

namespace himalaya {
namespace mh2_fo {

#define SUM(...) (get_sum(__VA_ARGS__)(__VA_ARGS__))

#define get_sum(...) get_sum_macro(__VA_ARGS__, sum_user_t, sum_ptrdiff_t,)

#define get_sum_macro(_1, _2, _3, _4, _5, name, ...) name

#define sum_ptrdiff_t(idx, ini, fin, expr)      \
    sum<std::ptrdiff_t>((ini), (fin), [&](std::ptrdiff_t (idx)) { return (expr); })

#define sum_user_t(type, idx, ini, fin, expr)	\
    sum<type>((ini), (fin), [&](type (idx)) { return (expr); })

template<typename T>
struct is_eigen_type
{
    static constexpr auto value =
	std::is_base_of<Eigen::EigenBase<T>, T>::value;
};

template<typename Idx, typename Function, bool isEigenType>
struct EvalEigenXprImpl {
    static auto eval(Idx i, Function f) -> decltype(f(i)) {
	return f(i);
    }
};

template<typename Idx, typename Function>
struct EvalEigenXprImpl<Idx, Function, true> {
    static auto eval(Idx i, Function f) ->
	typename std::remove_reference<decltype(f(i).eval())>::type
    {
	return f(i).eval();
    }
};

template<typename Idx, typename Function>
auto EvalEigenXpr(Idx i, Function f) ->
    decltype(
	EvalEigenXprImpl<Idx, Function, is_eigen_type<decltype(f(i))>::value>::
	eval(i, f))
{
    return
	EvalEigenXprImpl<Idx, Function, is_eigen_type<decltype(f(i))>::value>::
	eval(i, f);
}

template<typename T, bool isEigenType>
struct create_zero {
    static const T zero() {
	return T();
    }
};

template<typename T>
struct create_zero<T, true> {
    static const T zero() {
	T z;
	z.setZero();
	return z;
    }
};

template<class Idx, class Function>
auto sum(Idx ini, Idx fin, Function f) -> decltype(EvalEigenXpr<Idx>(ini, f))
{
    using Evaled = decltype(EvalEigenXpr<Idx>(ini, f));
    using Acc = typename std::remove_cv<Evaled>::type;
    Acc s = create_zero<Acc, is_eigen_type<Evaled>::value>::zero();
    for (Idx i = ini; i <= fin; i++) s += f(i);
    return s;
}

} // namespace mh2_fo
} // namespace himalaya
