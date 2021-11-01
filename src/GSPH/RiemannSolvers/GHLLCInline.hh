namespace Spheral {


//--------------------------------------------------------
// setter/getter for the gravitational acceleration
//--------------------------------------------------------
template<typename Dimension>
typename Dimension::Vector
GHLLC<Dimension>::
gravitationalAcceleration() const{
  return mGravitationalAcceleration;
}

template<typename Dimension>
void
GHLLC<Dimension>::
gravitationalAcceleration(const typename Dimension::Vector g) {
  mGravitationalAcceleration = g;
}

}