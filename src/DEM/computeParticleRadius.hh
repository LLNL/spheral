//---------------------------------Spheral++------------------------------------
// Compute the SPH mass density summation.
//------------------------------------------------------------------------------
#ifndef __Spheral__computeParticleRadius__
#define __Spheral__computeParticleRadius__

namespace Spheral {

  // Forward declarations.
  template<typename Dimension, typename DataType> class FieldList;

  template<typename Dimension>
  void
  computeParticleRadius(const FieldList<Dimension, typename Dimension::SymTensor>& H,
                              FieldList<Dimension, typename Dimension::Scalar>& particleRadius);

}

#endif
