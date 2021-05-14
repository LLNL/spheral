//---------------------------------Spheral++----------------------------------//
// uniform_random
//
// Encapsulate a random number generator to generate numbers in [0,1).
//
// We also require this generator be able to serialize/deserialize, both for
// restart and communication.
//
// Created by JMO, Mon May 10 16:02:11 PDT 2021
//----------------------------------------------------------------------------//

#include "uniform_random.hh"
#include "Utilities/packElement.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructors
//------------------------------------------------------------------------------
uniform_random::uniform_random(const size_t seed,
                               const double minValue,
                               const double maxValue):
  mGen(seed),
  mRan(minValue, maxValue),
  mSeed(seed),
  mNumCalls(0u),
  mMin(minValue),
  mMax(maxValue) {
}

uniform_random::uniform_random(const uniform_random& rhs):
  mGen(rhs.mGen),
  mRan(rhs.mRan),
  mSeed(rhs.mSeed),
  mNumCalls(rhs.mNumCalls),
  mMin(rhs.mMin),
  mMax(rhs.mMax) {
}

uniform_random&
uniform_random::operator=(const uniform_random& rhs) {
  if (this != &rhs) {
    mGen = rhs.mGen;
    mRan = rhs.mRan;
    mSeed = rhs.mSeed;
    mNumCalls = rhs.mNumCalls;
    mMin = rhs.mMin;
    mMax = rhs.mMax;
  }
  return *this;
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
uniform_random::~uniform_random() {
}


//------------------------------------------------------------------------------
// seed
//------------------------------------------------------------------------------
size_t
uniform_random::seed() const {
  return mSeed;
}

//------------------------------------------------------------------------------
// numCalls
//------------------------------------------------------------------------------
size_t
uniform_random::numCalls() const {
  return mNumCalls;
}

//------------------------------------------------------------------------------
// Set the seed
//------------------------------------------------------------------------------
void
uniform_random::seed(const size_t val) {
  mSeed = val;
  mGen.seed(val);
}


//------------------------------------------------------------------------------
// Set the range
//------------------------------------------------------------------------------
void
uniform_random::range(const double a, const double b) {
  mRan = std::uniform_real_distribution<double>(a, b);
}

//------------------------------------------------------------------------------
// Advance the state number of steps, as though we generated that many random
// values.
//------------------------------------------------------------------------------
void
uniform_random::advance(const size_t n) {
  // mGen.discard(n * mGen.state_size);
  // mNumCalls += n;
  for (auto i = 0u; i < n; ++i) (*this)();
}

//------------------------------------------------------------------------------
// Serialize to a buffer
//------------------------------------------------------------------------------
void
uniform_random::serialize(std::vector<char>& buffer) const {
  packElement(mSeed, buffer);
  packElement(mNumCalls, buffer);
  packElement(mMin, buffer);
  packElement(mMax, buffer);
}

//------------------------------------------------------------------------------
// Deserialize from a buffer
//------------------------------------------------------------------------------
void
uniform_random::deserialize(std::vector<char>::const_iterator& itr,
                            const std::vector<char>::const_iterator endItr) {
  unpackElement(mSeed, itr, endItr);
  unpackElement(mNumCalls, itr, endItr);
  unpackElement(mMin, itr, endItr);
  unpackElement(mMax, itr, endItr);
  this->seed(mSeed);
  this->range(mMin, mMax);
  this->advance(mNumCalls);
}

}
