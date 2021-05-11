namespace Spheral {

//------------------------------------------------------------------------------
// Return a random value
//------------------------------------------------------------------------------
inline
double
uniform_random_01::operator()() {
  return mRan(mGen);
}

//------------------------------------------------------------------------------
// ==
//------------------------------------------------------------------------------
inline
bool
uniform_random_01::operator==(const uniform_random_01& rhs) const {
  return (mGen == rhs.mGen) && (mRan == rhs.mRan);
}

//------------------------------------------------------------------------------
// !=
//------------------------------------------------------------------------------
inline
bool
uniform_random_01::operator!=(const uniform_random_01& rhs) const {
  return (mGen != rhs.mGen) || (mRan == rhs.mRan);
}

}
