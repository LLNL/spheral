namespace Spheral {

//------------------------------------------------------------------------------
// Singleton instance method.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
NestedGridDistributedBoundary<Dimension>&
NestedGridDistributedBoundary<Dimension>::
instance() {
  static NestedGridDistributedBoundary<Dimension> theInstance;
  return theInstance;
}

}
