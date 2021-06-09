//---------------------------------Spheral++----------------------------------//
// uniform_random
//
// Encapsulate a uniform random number generator for an arbitrary range of
// values: defaults to range [0, 1).
//
// We also require this generator be able to serialize/deserialize, both for
// restart and communication.
//
// Created by JMO, Mon May 10 16:02:11 PDT 2021
//----------------------------------------------------------------------------//
#ifndef __Spheral_uniform_random__
#define __Spheral_uniform_random__

#include <vector>
#include <random>

namespace Spheral {

class uniform_random {
public:
  //--------------------------- Public Interface ---------------------------//
  // Constructor, destructor.
  uniform_random(const size_t seed = std::random_device()(),
                 const double minVal = 0.0,
                 const double maxVal = 1.0);
  uniform_random(const uniform_random& rhs);
  uniform_random& operator=(const uniform_random& rhs);
  ~uniform_random();

  // Generate a random number.
  double operator()();

  // Poke the internal state.
  size_t seed() const;
  size_t numCalls() const;
  double min() const;
  double max() const;
  void seed(const size_t val);                // Set the seed value
  void advance(const size_t n);               // Advance n times in random sequence
  void range(const double a, const double b); // Set the possible range of values [a,b)

  // Comparison
  bool operator==(const uniform_random& rhs) const;
  bool operator!=(const uniform_random& rhs) const;

  // Methods for serializing our state to/from vector<char> buffers.
  void serialize(std::vector<char>& buffer) const;
  void deserialize(std::vector<char>::const_iterator& itr,
                   const std::vector<char>::const_iterator endItr);

private:
  //--------------------------- Private Interface ---------------------------//
  std::mt19937 mGen;
  std::uniform_real_distribution<double> mRan;
  size_t mSeed, mNumCalls;
  double mMin, mMax;
};

}

#include "uniform_random_Inline.hh"

#else

namespace Spheral {
  class uniform_random;
}

#endif
