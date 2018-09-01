//---------------------------------Spheral++----------------------------------//
// A simple form for the artificial viscosity due to Monaghan & Gingold.
//----------------------------------------------------------------------------//
#include "NohViscosity.hh"
#include "Field/FieldList.hh"
#include "DataBase/DataBase.hh"
#include "Boundary/Boundary.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Default constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
NohViscosity<Dimension>::
NohViscosity():
  MonaghanGingoldViscosity<Dimension>(),
  mNodePositions(),
  mHfield(),
  mCurrentTime(0.0) {
}

//------------------------------------------------------------------------------
// Construct with the given value for the linear and quadratic coefficients.
//------------------------------------------------------------------------------
template<typename Dimension>
NohViscosity<Dimension>::
NohViscosity(Scalar Clinear, Scalar Cquadratic):
  MonaghanGingoldViscosity<Dimension>(Clinear, Cquadratic),
  mNodePositions(),
  mHfield(),
  mCurrentTime(0.0) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
NohViscosity<Dimension>::
~NohViscosity() {
}

//------------------------------------------------------------------------------
// Initialization method for the beginning of a cycle.
//------------------------------------------------------------------------------
template<typename Dimension>
void
NohViscosity<Dimension>::
initialize(const DataBase<Dimension>& dataBase,
           typename ArtificialViscosity<Dimension>::ConstBoundaryIterator boundaryBegin,
           typename ArtificialViscosity<Dimension>::ConstBoundaryIterator boundaryEnd,
           const typename Dimension::Scalar time,
           const typename Dimension::Scalar dt,
           const TableKernel<Dimension>& W) {

  // Invoke the base class initialization.
  MonaghanGingoldViscosity<Dimension>::initialize(dataBase, 
                                                  boundaryBegin,
                                                  boundaryEnd,
                                                  time,
                                                  dt,
                                                  W);

  // Store the time.
  mCurrentTime = time;

  // Copy the node positions and H tensors.
  mNodePositions = dataBase.fluidPosition();
  mHfield = dataBase.fluidHfield();
}

//------------------------------------------------------------------------------
// Determine the contribution that mimics the artificial internal energy.
// Just call the MonaghanGingold base method if Q should be active.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
NohViscosity<Dimension>::
viscousInternalEnergy(const NodeIDIterator<Dimension>& nodeI,
                      const NodeIDIterator<Dimension>& nodeJ,
                      const Vector& rij, const Vector& vij,
                      const Vector& etai, const Vector& etaj,
                      const Scalar ci, const Scalar cj) const {

  const Vector ri = mNodePositions(nodeI);
  const Vector rj = mNodePositions(nodeJ);
  const Vector rjiNorm = (rj - ri).unitVector();
  const Scalar hi = 1.0/((mHfield(nodeI)*rjiNorm).magnitude());
  const Scalar rshock = mCurrentTime/3.0;
  if (ri.magnitude() <= rshock + 2.0*hi) {
    return MonaghanGingoldViscosity<Dimension>::viscousInternalEnergy(nodeI, nodeJ,
                                                                      rij, vij,
                                                                      etai, etaj,
                                                                      ci, cj);
  } else {
    return 0.0;
  }
}

}

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
namespace Spheral {
  template class NohViscosity< Dim<1> >;
  template class NohViscosity< Dim<2> >;
  template class NohViscosity< Dim<3> >;
}
