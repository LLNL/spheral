namespace Spheral {

//------------------------------------------------------------------------------
// Singleton instance method.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
BoundingVolumeDistributedBoundary<Dimension>&
BoundingVolumeDistributedBoundary<Dimension>::
instance() {
  static BoundingVolumeDistributedBoundary<Dimension> theInstance;
  return theInstance;
}

}
