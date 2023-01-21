//---------------------------------Spheral++----------------------------------//
// SPHIntegrationKernel
//----------------------------------------------------------------------------//
#ifndef __Spheral_SPHIntegrationKernel__
#define __Spheral_SPHIntegrationKernel__

#include "Field/FieldList.hh"
#include "Geometry/Dimension.hh"
#include "IntegrationKernel.hh"

#include <utility>
#include <vector>

namespace Spheral {

template<typename Dimension>
class SPHIntegrationKernel : public IntegrationKernel<Dimension> {
public:
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::SymTensor SymTensor;
  
  SPHIntegrationKernel(const TableKernel<Dimension>& kernel);

  virtual double extent(const Scalar Hmult) const { return mKernel.kernelExtent() / Hmult; }
  
  virtual void evaluate(const Vector& xp,
                        const std::vector<std::pair<int, int>>& indices,
                        const FieldList<Dimension, Vector>& position,
                        const FieldList<Dimension, SymTensor>& H,
                        const FieldList<Dimension, Scalar>& volume,
                        const Scalar Hmult,
                        std::vector<Scalar>& values,
                        std::vector<Vector>& dvalues) const override;
  
protected:
  const TableKernel<Dimension>& mKernel;
};

} // end namespace Spheral

#endif
