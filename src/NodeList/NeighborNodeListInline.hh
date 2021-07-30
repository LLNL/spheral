namespace Spheral {

//------------------------------------------------------------------------------
// Access the maximum number of neighbors to allow.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
unsigned 
NeighborNodeList<Dimension>::maxNumNeighbors() const {
  return mMaxNumNeighbors;
}

template<typename Dimension>
void
NeighborNodeList<Dimension>::maxNumNeighbors(unsigned val) {
  mMaxNumNeighbors = val;
}


}
