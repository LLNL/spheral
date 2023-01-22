//---------------------------------Spheral++----------------------------------//
// SPHIntegrationKernel
//----------------------------------------------------------------------------//
#include "SPHIntegrationKernel.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template<typename Dimension>
SPHIntegrationKernel<Dimension>::
SPHIntegrationKernel(const TableKernel<Dimension>& kernel):
  mKernel(kernel) {
}

//------------------------------------------------------------------------------
// Evaluate the functions
//------------------------------------------------------------------------------
template<typename Dimension>
void
SPHIntegrationKernel<Dimension>::
evaluate(const Vector& xp,
         const std::vector<std::pair<int, int>>& indices,
         const FieldList<Dimension, Vector>& position,
         const FieldList<Dimension, SymTensor>& H,
         const FieldList<Dimension, Scalar>& volume,
         const Scalar Hmult,
         std::vector<Scalar>& values,
         std::vector<Vector>& dvalues) const {
  const auto numIndices = indices.size();
  CHECK(values.size() == numIndices
        && dvalues.size() == numIndices);
  for (auto j = 0u; j < numIndices; ++j) {
    auto nodeListj = indices[j].first;
    auto nodej = indices[j].second;
    const auto xj = position(nodeListj, nodej);
    const auto Hj = H(nodeListj, nodej) * Hmult;
    const auto xpj = xp - xj;
    const auto eta = Hj * xpj;
    const auto etaMag = eta.magnitude();
    const auto etaUnit = eta.unitVector();
    const auto Hdet = Hj.Determinant();
    const auto HetaUnit = Hj * etaUnit;
    
    values[j] = mKernel.kernelValue(etaMag, Hdet);
    dvalues[j] = HetaUnit * mKernel.gradValue(etaMag, Hdet);
  }
}

} // end namespace Spheral
