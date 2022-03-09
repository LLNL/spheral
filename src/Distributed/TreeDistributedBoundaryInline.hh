namespace Spheral {

//------------------------------------------------------------------------------
// Singleton instance method.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
TreeDistributedBoundary<Dimension>&
TreeDistributedBoundary<Dimension>::
instance() {
  static TreeDistributedBoundary<Dimension> theInstance;
  return theInstance;
}

}
