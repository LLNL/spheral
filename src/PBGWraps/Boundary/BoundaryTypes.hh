#ifndef __PBGWRAPS_BOUNDARYTYPES__
#define __PBGWRAPS_BOUNDARYTYPES__

#include <vector>
#include "Boundary/Boundary.hh"
#include "Boundary/PlanarBoundary.hh"
#include "Boundary/ReflectingBoundary.hh"
#include "Boundary/RigidBoundary.hh"
#include "Boundary/PeriodicBoundary.hh"
#include "Boundary/ConstantVelocityBoundary.hh"
#include "Boundary/ConstantXVelocityBoundary.hh"
#include "Boundary/ConstantYVelocityBoundary.hh"
#include "Boundary/ConstantZVelocityBoundary.hh"
#include "Boundary/ConstantRVelocityBoundary.hh"
#include "Boundary/ConstantBoundary.hh"
#include "Boundary/CRKSPHVoidBoundary.hh"
#include "Boundary/SphericalBoundary.hh"
#include "Boundary/CylindricalBoundary.hh"
#include "Boundary/AxialSymmetryBoundary.hh"
#include "Boundary/AxisBoundaryRZ.hh"

#include "PBGWraps/referenceAsPointer.hh"

//------------------------------------------------------------------------------
// Names!
//------------------------------------------------------------------------------
namespace Spheral {

typedef Boundary<Dim<1> > Boundary1d;
typedef Boundary<Dim<2> > Boundary2d;
typedef Boundary<Dim<3> > Boundary3d;

typedef PlanarBoundary<Dim<1> > PlanarBoundary1d;
typedef PlanarBoundary<Dim<2> > PlanarBoundary2d;
typedef PlanarBoundary<Dim<3> > PlanarBoundary3d;

typedef ReflectingBoundary<Dim<1> > ReflectingBoundary1d;
typedef ReflectingBoundary<Dim<2> > ReflectingBoundary2d;
typedef ReflectingBoundary<Dim<3> > ReflectingBoundary3d;

typedef RigidBoundary<Dim<1> > RigidBoundary1d;
typedef RigidBoundary<Dim<2> > RigidBoundary2d;
typedef RigidBoundary<Dim<3> > RigidBoundary3d;

typedef PeriodicBoundary<Dim<1> > PeriodicBoundary1d;
typedef PeriodicBoundary<Dim<2> > PeriodicBoundary2d;
typedef PeriodicBoundary<Dim<3> > PeriodicBoundary3d;

typedef ConstantVelocityBoundary<Dim<1> > ConstantVelocityBoundary1d;
typedef ConstantVelocityBoundary<Dim<2> > ConstantVelocityBoundary2d;
typedef ConstantVelocityBoundary<Dim<3> > ConstantVelocityBoundary3d;

typedef ConstantXVelocityBoundary<Dim<1> > ConstantXVelocityBoundary1d;
typedef ConstantXVelocityBoundary<Dim<2> > ConstantXVelocityBoundary2d;
typedef ConstantXVelocityBoundary<Dim<3> > ConstantXVelocityBoundary3d;

typedef ConstantYVelocityBoundary<Dim<2> > ConstantYVelocityBoundary2d;
typedef ConstantYVelocityBoundary<Dim<3> > ConstantYVelocityBoundary3d;

typedef ConstantZVelocityBoundary<Dim<3> > ConstantZVelocityBoundary3d;

typedef ConstantRVelocityBoundary<Dim<1> > ConstantRVelocityBoundary1d;
typedef ConstantRVelocityBoundary<Dim<2> > ConstantRVelocityBoundary2d;
typedef ConstantRVelocityBoundary<Dim<3> > ConstantRVelocityBoundary3d;

typedef ConstantBoundary<Dim<1> > ConstantBoundary1d;
typedef ConstantBoundary<Dim<2> > ConstantBoundary2d;
typedef ConstantBoundary<Dim<3> > ConstantBoundary3d;

typedef CRKSPHVoidBoundary<Dim<1> > CRKSPHVoidBoundary1d;
typedef CRKSPHVoidBoundary<Dim<2> > CRKSPHVoidBoundary2d;
typedef CRKSPHVoidBoundary<Dim<3> > CRKSPHVoidBoundary3d;

//------------------------------------------------------------------------------
// Extract the BoundaryNodes for the given NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Boundary<Dimension>::BoundaryNodes*
accessBoundaryNodesFromBoundary(Boundary<Dimension>& self,
                                NodeList<Dimension>& nodes) {
  return &(self.accessBoundaryNodes(nodes));
}

//------------------------------------------------------------------------------
// Extract the various sets from a BoundaryNode instance.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
std::vector<int>*
controlNodesFromBoundaryNodes(typename Boundary<Dimension>::BoundaryNodes& bn) {
  return &(bn.controlNodes);
}

template<typename Dimension>
inline
std::vector<int>*
ghostNodesFromBoundaryNodes(typename Boundary<Dimension>::BoundaryNodes& bn) {
  return &(bn.ghostNodes);
}

template<typename Dimension>
inline
std::vector<int>*
violationNodesFromBoundaryNodes(typename Boundary<Dimension>::BoundaryNodes& bn) {
  return &(bn.violationNodes);
}

//------------------------------------------------------------------------------
// Support dynamic_casting from a Boundary to one of descendent types.  Kind of
// an ugly collision between having C++ vector's of boundaries and needing
// dynamic_cast at the python level.
//------------------------------------------------------------------------------
template<typename Boundary1, typename Boundary2>
Boundary2*
dynamicCastBoundary(Boundary1* boundPtr) {
  if (Boundary2* result = dynamic_cast<Boundary2*>(boundPtr)) {
    return result;
  } else {
    PyErr_SetString(PyExc_ValueError, "Failed attempt to dynamic_cast a Boundary class");
    return NULL;
  }
}

}

typedef std::vector<Spheral::Boundary1d*> vector_of_Boundary1d;
typedef std::vector<Spheral::Boundary2d*> vector_of_Boundary2d;
typedef std::vector<Spheral::Boundary3d*> vector_of_Boundary3d;

#endif
