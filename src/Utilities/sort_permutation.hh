//------------------------------------------------------------------------------
// Sort a vector by the ordering of another.
// 
// sort_permutation: computes the permutation that would sort a given vector
// apply_permutation: returns the result of applying a permutation 
// apply_permutation_inplace: same as above in place
//
// Based on the examples from
// http://stackoverflow.com/questions/17074324/how-can-i-sort-two-vectors-in-the-same-way-with-criteria-that-uses-only-one-of
//------------------------------------------------------------------------------
#ifndef __Spheral_sort_permutation__
#define __Spheral_sort_permutation__

#include <vector>
#include <algorithm>

namespace Spheral {

//------------------------------------------------------------------------------
// Find a permutation.
//------------------------------------------------------------------------------
template <typename T, typename Compare>
std::vector<std::size_t> sort_permutation(const std::vector<T>& vec,
                                          Compare& compare) {
  std::vector<std::size_t> p(vec.size());
  std::iota(p.begin(), p.end(), 0);
  std::sort(p.begin(), p.end(),
            [&](std::size_t i, std::size_t j){ return compare(vec[i], vec[j]); });
  return p;
}

//------------------------------------------------------------------------------
// Apply a permutation (returns a new copy).
//------------------------------------------------------------------------------
template <typename T>
std::vector<T> apply_permutation(const std::vector<T>& vec,
                                 const std::vector<std::size_t>& p) {
  std::vector<T> sorted_vec(vec.size());
  std::transform(p.begin(), p.end(), sorted_vec.begin(),
                 [&](std::size_t i){ return vec[i]; });
  return sorted_vec;
}

//------------------------------------------------------------------------------
// Apply a permutation (in-place version).
//------------------------------------------------------------------------------
template <typename T>
void apply_permutation_in_place(std::vector<T>& vec,
                                const std::vector<std::size_t>& p) {
  std::vector<bool> done(vec.size());
  for (std::size_t i = 0; i < vec.size(); ++i) {
    if (done[i]) continue;
    done[i] = true;
    std::size_t prev_j = i;
    std::size_t j = p[i];
    while (i != j) {
      std::swap(vec[prev_j], vec[j]);
      done[j] = true;
      prev_j = j;
      j = p[j];
    }
  }
}

}

#endif
