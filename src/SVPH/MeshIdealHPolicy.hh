//---------------------------------Spheral++----------------------------------//
// MeshIdealHPolicy
//
// Created by JMO, Fri Aug 16 14:47:57 PDT 2013
//----------------------------------------------------------------------------//
#ifndef __Spheral_MeshIdealHPolicy_hh__
#define __Spheral_MeshIdealHPolicy_hh__

#include <string>

#include "DataBase/UpdatePolicyBase.hh"
#include "NodeList/SmoothingScaleBase.hh"

namespace Spheral {

template<typename Dimension>
class MeshIdealHPolicy: 
    public UpdatePolicyBase<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  using KeyType = typename UpdatePolicyBase<Dimension>::KeyType;

  // Constructors, destructor.
  MeshIdealHPolicy(const SmoothingScaleBase<Dimension>& smoothingScaleBase,
                   const Scalar hmin,
                   const Scalar hmax,
                   const Scalar hminratio,
                   const Scalar nPerh);
  virtual ~MeshIdealHPolicy();
  
  // Overload the methods describing how to update Fields.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt) override;

  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const override;

private:
  //--------------------------- Private Interface ---------------------------//
  static bool mFired;
  const SmoothingScaleBase<Dimension>& mSmoothingScaleBase;
  Scalar mhmin, mhmax, mhminratio, mnPerh;

  MeshIdealHPolicy(const MeshIdealHPolicy& rhs);
  MeshIdealHPolicy& operator=(const MeshIdealHPolicy& rhs);
};

}

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class MeshIdealHPolicy;
}

#endif
