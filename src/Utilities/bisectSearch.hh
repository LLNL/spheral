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
#include <algorithm>

namespace Spheral {

// General iterator based method.
template<typename DataType, typename IteratorType>
inline
int
bisectSearch(const IteratorType& begin,
             const IteratorType& end,
             const DataType& val) {
  const int n = std::distance(begin, end);
  if (n > 1) {
    const bool ascnd = (*(end - 1) >= *begin);
    int jl = -1;
    int ju = n;
    while (ju - jl > 1) {
      const int jm = (ju + jl)/2;
      if ((val >= *(begin + jm)) == ascnd) {
        jl = jm;
      } else {
        ju = jm;
      }
    }

    // Post conditions.
    ENSURE(ju - jl == 1);
    ENSURE((jl == -1 and (val <= *begin) == ascnd) or
           (jl == n - 1 and (val >= *(end - 1)) == ascnd) or
           (((val >= *(begin + jl)) == ascnd) and
            ((val <= *(begin + ju)) == ascnd)));
    return jl;

  } else if (n == 1) {
    if (val <= *begin) {
      return -1;
    } else {
      return 0;
    }
  } else {
    return -1;
  }
}

// Specialized for a std::vector for backwards compatibility.
template<typename DataType>
inline
int
bisectSearch(const std::vector<DataType>& table, 
             const DataType& val) {
  return bisectSearch(table.begin(), table.end(), val);
}

}

#endif
