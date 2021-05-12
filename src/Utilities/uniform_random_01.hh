//---------------------------------Spheral++----------------------------------//
// uniform_random_01
//
// Encapsulate a random number generator to generate numbers in [0,1).
//
// We also require this generator be able to serialize/deserialize, both for
// restart and communication.
//
// Created by JMO, Mon May 10 16:02:11 PDT 2021
//----------------------------------------------------------------------------//
#ifndef __Spheral_uniform_random_01__
#define __Spheral_uniform_random_01__

#include <vector>
#include <random>

namespace Spheral {

class uniform_random_01 {
public:
  //--------------------------- Public Interface ---------------------------//
  // Constructor, destructor.
  uniform_random_01();
  uniform_random_01(const size_t seed);
  uniform_random_01(const uniform_random_01& rhs);
  uniform_random_01& operator=(const uniform_random_01& rhs);
  ~uniform_random_01();

  // Generate a random number.
  double operator()();

  // Poke the internal state.
  size_t seed() const;
  size_t numCalls() const;
  void seed(const size_t val);   // Set the seed value
  void advance(const size_t n);  // Advance n times in random sequence

  // Comparison
  bool operator==(const uniform_random_01& rhs) const;
  bool operator!=(const uniform_random_01& rhs) const;

  // Methods for serializing our state to/from vector<char> buffers.
  void serialize(std::vector<char>& buffer) const;
  void deserialize(std::vector<char>::const_iterator& itr,
                   const std::vector<char>::const_iterator endItr);

private:
  //--------------------------- Private Interface ---------------------------//
  std::mt19937 mGen;
  std::uniform_real_distribution<double> mRan;
  size_t mSeed, mNumCalls;
};

}

#include "uniform_random_01_Inline.hh"

#else

namespace Spheral {
  class uniform_random_01;
}

#endif
