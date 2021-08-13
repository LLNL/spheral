namespace Spheral {

//------------------------------------------------------------------------------
// Return a random value
//------------------------------------------------------------------------------
inline
double
uniform_random::operator()() {
  ++mNumCalls;
  return mRan(mGen);
}

//------------------------------------------------------------------------------
// ==
//------------------------------------------------------------------------------
inline
bool
uniform_random::operator==(const uniform_random& rhs) const {
  return (mGen == rhs.mGen) && (mRan == rhs.mRan);
}

//------------------------------------------------------------------------------
// !=
//------------------------------------------------------------------------------
inline
bool
uniform_random::operator!=(const uniform_random& rhs) const {
  return (mGen != rhs.mGen) || (mRan != rhs.mRan);
}

}
