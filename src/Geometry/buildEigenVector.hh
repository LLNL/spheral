// Helper method to compute the eigen vector corresponding to a given 
// eigen value for a 3x3 symmetric matrix.
#include "Dimension.hh"

namespace Spheral {

Dim<3>::Vector
buildEigenVector(const Dim<3>::SymTensor& A,
                 const double lambda);

Dim<3>::Vector
buildUniqueEigenVector(const Dim<3>::SymTensor& A,
                       const double lambda,
                       const Dim<3>::Vector& U0,
                       const Dim<3>::Vector& U1);

}
