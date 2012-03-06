#include "ExtendJacobi.hh"
#include "Geometry/Jacobi2.hh"
#include "Geometry/Dimension.hh"

namespace Spheral {

int 
jacobiDiagonalize3d(const GeomSymmetricTensor<3>& A,
                    GeomTensor<3>& eigenVectors,
                    GeomVector<3>& eigenValues,
                    const double convergenceThreshold,
                    const int maxSweeps) {
  return Spheral::Geometry::jacobiDiagonalize< Spheral::Dim<3> >(A,
                                                        eigenVectors,
                                                        eigenValues,
                                                        convergenceThreshold,
                                                        maxSweeps);
}

}
