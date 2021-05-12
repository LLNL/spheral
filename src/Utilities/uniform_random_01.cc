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

#include "uniform_random_01.hh"
#include "Utilities/packElement.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructors
//------------------------------------------------------------------------------
uniform_random_01::uniform_random_01():
  mGen(0u),
  mRan(0.0, 1.0),
  mSeed(0u),
  mNumCalls(0u) {
}

uniform_random_01::uniform_random_01(const size_t seed):
  mGen(seed),
  mRan(0.0, 1.0),
  mSeed(seed),
  mNumCalls(0u) {
}

uniform_random_01::uniform_random_01(const uniform_random_01& rhs):
  mGen(rhs.mGen),
  mRan(rhs.mRan),
  mSeed(rhs.mSeed),
  mNumCalls(rhs.mNumCalls) {
}

uniform_random_01&
uniform_random_01::operator=(const uniform_random_01& rhs) {
  if (this != &rhs) {
    mGen = rhs.mGen;
    mRan = rhs.mRan;
    mSeed = rhs.mSeed;
    mNumCalls = rhs.mNumCalls;
  }
  return *this;
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
uniform_random_01::~uniform_random_01() {
}


//------------------------------------------------------------------------------
// seed
//------------------------------------------------------------------------------
size_t
uniform_random_01::seed() const {
  return mSeed;
}

//------------------------------------------------------------------------------
// numCalls
//------------------------------------------------------------------------------
size_t
uniform_random_01::numCalls() const {
  return mNumCalls;
}

//------------------------------------------------------------------------------
// Set the seed
//------------------------------------------------------------------------------
void
uniform_random_01::seed(const size_t val) {
  mSeed = val;
  mGen.seed(val);
}

//------------------------------------------------------------------------------
// Advance the state number of steps, as though we generated that many random
// values.
//------------------------------------------------------------------------------
void
uniform_random_01::advance(const size_t n) {
  // mGen.discard(n * mGen.state_size);
  // mNumCalls += n;
  for (auto i = 0u; i < n; ++i) (*this)();
}

//------------------------------------------------------------------------------
// Serialize to a buffer
//------------------------------------------------------------------------------
void
uniform_random_01::serialize(std::vector<char>& buffer) const {
  packElement(mSeed, buffer);
  packElement(mNumCalls, buffer);
}

//------------------------------------------------------------------------------
// Deserialize from a buffer
//------------------------------------------------------------------------------
void
uniform_random_01::deserialize(std::vector<char>::const_iterator& itr,
                               const std::vector<char>::const_iterator endItr) {
  unpackElement(mSeed, itr, endItr);
  unpackElement(mNumCalls, itr, endItr);
  this->seed(mSeed);
  this->advance(mNumCalls);
}

}
