//---------------------------------Spheral++----------------------------------//
// Bisection search of an array to find the indicies bracketing the given value.
// The return value is the lower bracketing index; -1 is returned if the value
// is off the lower end of the scale, and j = table.size() - 1 if the value is
// off the upper end.
// We assume here that the input table is sorted.
//----------------------------------------------------------------------------//
#ifndef __Spheral_bisectSearch__
#define __Spheral_bisectSearch__

#include <vector>

namespace Spheral {

template<typename DataType>
inline
int
bisectSearch(const std::vector<DataType>& table, 
             const DataType& val) {
  const int n = table.size();
  const bool ascnd = (table[n - 1] >= table[0]);
  int jl = -1;
  int ju = n;
  while (ju - jl > 1) {
    const int jm = (ju + jl)/2;
    if ((val >= table[jm]) == ascnd) {
      jl = jm;
    } else {
      ju = jm;
    }
  }

  // Post conditions.
  ENSURE(ju - jl == 1);
  ENSURE((jl == -1 && (val <= table[0]) == ascnd) ||
         (jl == n - 1 && (val >= table[n - 1]) == ascnd) ||
         (((val >= table[jl]) == ascnd) &&
          ((val <= table[ju]) == ascnd)));
  return jl;
}

}

#endif
