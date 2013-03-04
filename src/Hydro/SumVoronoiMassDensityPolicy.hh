//---------------------------------Spheral++----------------------------------//
// SumVoronoiMassDensityPolicy -- An implementation of UpdatePolicyBase 
// specialized for the updating the mass density according to the specific 
// volume from the Voronoi tesselation.  This version uses a weighted sum of the
// local mass and volume to do the deed.
//
// Created by JMO, Mon Aug  1 10:48:03 PDT 2011
//----------------------------------------------------------------------------//
#ifndef __Spheral_SumVoronoiMassDensityPolicy_hh__
#define __Spheral_SumVoronoiMassDensityPolicy_hh__

#include <string>

#include "DataBase/ReplaceState.hh"
#include "Kernel/TableKernel.hh"
#include "Physics/Physics.hh"

namespace Spheral {

template<typename Dimension>
class SumVoronoiMassDensityPolicy: public ReplaceState<Dimension, typename Dimension::Scalar> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename ReplaceState<Dimension, Scalar>::KeyType KeyType;

  // Constructors, destructor.
  SumVoronoiMassDensityPolicy(const KernelSpace::TableKernel<Dimension>& W,
                              const PhysicsSpace::Physics<Dimension>& package,
                              const double rhoMin,
                              const double rhoMax);
  virtual ~SumVoronoiMassDensityPolicy();
  
  // Overload the methods describing how to update Fields.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt);

  virtual void updateAsIncrement(const KeyType& key,
                                 State<Dimension>& state,
                                 StateDerivatives<Dimension>& derivs,
                                 const double multiplier,
                                 const double t,
                                 const double dt);

  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const;

private:
  //--------------------------- Private Interface ---------------------------//
  bool mEnforceBoundaries;
  const KernelSpace::TableKernel<Dimension>& mW;
  const PhysicsSpace::Physics<Dimension>& mPackage;
  double mRhoMin, mRhoMax;
  SumVoronoiMassDensityPolicy(const SumVoronoiMassDensityPolicy& rhs);
  SumVoronoiMassDensityPolicy& operator=(const SumVoronoiMassDensityPolicy& rhs);
};

}

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class SumVoronoiMassDensityPolicy;
}

#endif
