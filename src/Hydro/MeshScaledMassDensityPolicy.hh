//---------------------------------Spheral++----------------------------------//
// MeshScaledMassDensityPolicy --
// Implements the algorithm for updating the mass density by computing the
// Voronoi tesselation of the point positions and using the cell volumes.
//
// Created by JMO, Mon Jan 24 22:47:32 PST 2011
//----------------------------------------------------------------------------//
#ifndef __Spheral_PressurePolicy_hh__
#define __Spheral_PressurePolicy_hh__

#include "DataBase/UpdatePolicyBase.hh"

#include <string>

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class FluidNodeList;
template<typename Dimension, typename DataType> class Field;

template<typename Dimension>
class MeshScaledMassDensityPolicy: public UpdatePolicyBase<Dimension, typename Dimension::Scalar> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  typedef typename Dimension::Scalar Scalar;
  typedef typename UpdatePolicyBase<Dimension, Scalar>::KeyType KeyType;

  // Constructors, destructor.
  MeshScaledMassDensityPolicy();
  virtual ~MeshScaledMassDensityPolicy();
  
  // Overload the methods describing how to update Fields.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt);

  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension, Scalar>& rhs) const;

private:
  //--------------------------- Private Interface ---------------------------//
  MeshScaledMassDensityPolicy(const MeshScaledMassDensityPolicy& rhs);
  MeshScaledMassDensityPolicy& operator=(const MeshScaledMassDensityPolicy& rhs);
};

}

#endif
