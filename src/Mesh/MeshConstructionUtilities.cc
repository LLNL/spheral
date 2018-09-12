#include <limits>

#include "MeshConstructionUtilities.hh"
#include "Geometry/Dimension.hh"

namespace Spheral {

#ifdef __IBMCPP__
  template<> const uint64_t MeshTraits<Dim<1>::Vector>::ncells = 1UL << 32;
  template<> const uint64_t MeshTraits<Dim<1>::Vector>::ncoarse = 1UL << 14;
  template<> const uint64_t MeshTraits<Dim<1>::Vector>::ncellsperbin = 1UL << 18;
  template<> const uint64_t MeshTraits<Dim<2>::Vector>::ncells = 1UL << 32;
  template<> const uint64_t MeshTraits<Dim<2>::Vector>::ncoarse = 1UL << 14;
  template<> const uint64_t MeshTraits<Dim<2>::Vector>::ncellsperbin = 1UL << 18;
  template<> const uint64_t MeshTraits<Dim<3>::Vector>::ncells = 1UL << 32;
  template<> const uint64_t MeshTraits<Dim<3>::Vector>::ncoarse = 1UL << 14;
  template<> const uint64_t MeshTraits<Dim<3>::Vector>::ncellsperbin = 1UL << 18;
#else
  template<typename Vector> const uint64_t MeshTraits<Vector>::ncells = 1UL << 32;
  template<typename Vector> const uint64_t MeshTraits<Vector>::ncoarse = 1UL << 14;
  template<typename Vector> const uint64_t MeshTraits<Vector>::ncellsperbin = 1UL << 18;
#endif
  template struct MeshTraits<Dim<1>::Vector>;
  template struct MeshTraits<Dim<2>::Vector>;
  template struct MeshTraits<Dim<3>::Vector>;

}
