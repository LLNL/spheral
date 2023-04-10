namespace Spheral {

template<typename Dimension>
inline
const std::vector<SolidBoundary<Dimension>*>&
DEMBoundaryPolicy<Dimension>::solidBoundaryConditions() const {
  return mSolidBoundaries;
}

}
