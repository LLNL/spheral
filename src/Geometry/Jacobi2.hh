// A simple home brewed implementation of the Jacobi transformation method
// for diagonalizing a real symmetric matrix.

namespace Spheral {

template<typename Dimension>
int jacobiDiagonalize(const typename Dimension::SymTensor& A,
                      typename Dimension::Tensor& eigenVectors,
                      typename Dimension::Vector& eigenValues,
                      const double convergenceThreshold = 1.0e-15,
                      const int maxSweeps = 50);

}
