#include "Geometry/GeomVector.hh"
#include "Geometry/GeomTensor.hh"
#include "Geometry/GeomSymmetricTensor.hh"

namespace Spheral {

int 
jacobiDiagonalize3d(const GeomSymmetricTensor<3>& A,
                    GeomTensor<3>& eigenVectors,
                    GeomVector<3>& eigenValues,
                    const double convergenceThreshold = 1.0e-15,
                    const int maxSweeps = 50);

}
