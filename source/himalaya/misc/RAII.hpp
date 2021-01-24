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
class RAII_tmp_set {
public:
   RAII_tmp_set(T& variable_, T new_value_) noexcept
      : variable(variable_)
      , old_value(variable_)
      {
         variable_ = new_value_;
      }
   RAII_tmp_set(const RAII_tmp_set&) = delete;
   RAII_tmp_set(RAII_tmp_set&&) noexcept = delete;
   ~RAII_tmp_set() { variable = old_value; }
   RAII_tmp_set& operator=(const RAII_tmp_set&) = delete;
   RAII_tmp_set& operator=(RAII_tmp_set&& other) noexcept = delete;
private:
   T& variable;
   T old_value;
};

} // namespace himalaya
