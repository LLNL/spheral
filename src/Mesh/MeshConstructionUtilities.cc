#include <limits>

#include "MeshConstructionUtilities.hh"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace MeshSpace {
    template<typename Vector> const uint64_t MeshTraits<Vector>::ncells = 1UL << 32;
    template<typename Vector> const uint64_t MeshTraits<Vector>::ncoarse = 1UL << 14;
    template<typename Vector> const uint64_t MeshTraits<Vector>::ncellsperbin = 1UL << 18;
  }
}
