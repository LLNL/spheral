//---------------------------------Spheral++----------------------------------//
// MeshPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the Mesh in the state.
//
// Created by JMO, Sat Feb 12 14:37:57 PST 2011
//----------------------------------------------------------------------------//
#ifndef __Spheral_MeshPolicy_hh__
#define __Spheral_MeshPolicy_hh__

#include <string>

#include "DataBase/UpdatePolicyBase.hh"
#include "Physics/Physics.hh"

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;

template<typename Dimension>
class MeshPolicy: public UpdatePolicyBase<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  typedef typename Dimension::Vector Vector;
  typedef typename UpdatePolicyBase<Dimension>::KeyType KeyType;

  // Constructors, destructor.
  MeshPolicy(const PhysicsSpace::Physics<Dimension>& package,
             const double voidThreshold = 2.0,
             const bool meshGhostNodes = false);
  MeshPolicy(const PhysicsSpace::Physics<Dimension>& package,
             const Vector& xmin,
             const Vector& xmax,
             const double voidThreshold = 2.0,
             const bool meshGhostNodes = false);
  virtual ~MeshPolicy();
  
  // Overload the methods describing how to update Fields.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt);

  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const;

private:
  //--------------------------- Private Interface ---------------------------//
  const PhysicsSpace::Physics<Dimension>& mPackage;
  double mVoidThreshold;
  bool mComputeBounds, mMeshGhostNodes;
  Vector mXmin, mXmax;

  MeshPolicy();
  MeshPolicy(const MeshPolicy& rhs);
  MeshPolicy& operator=(const MeshPolicy& rhs);
};

}

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class MeshPolicy;
}

#endif
