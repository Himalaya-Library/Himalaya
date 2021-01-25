// ====================================================================
// This file is part of Himalaya.
//
// Himalaya is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================

#pragma once

namespace himalaya {

/// temporarily sets a variable to a new value, and resets the value
/// to the old one when destoyed
template <class T>
class Temporarily_set {
public:
   Temporarily_set(T& variable_, T new_value_) noexcept
      : variable(variable_)
      , old_value(variable_)
      {
         variable_ = new_value_;
      }
   Temporarily_set(const Temporarily_set&) = delete;
   Temporarily_set(Temporarily_set&&) noexcept = delete;
   ~Temporarily_set() { variable = old_value; }
   Temporarily_set& operator=(const Temporarily_set&) = delete;
   Temporarily_set& operator=(Temporarily_set&& other) noexcept = delete;
private:
   T& variable;
   T old_value;
};

} // namespace himalaya
