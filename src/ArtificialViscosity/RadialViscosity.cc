//---------------------------------Spheral++----------------------------------//
// A specialized form of the Artificial viscosity to restrict the viscous 
// interactions to the radial direction.
//----------------------------------------------------------------------------//
#include "RadialViscosity.hh"
#include "Field/FieldList.hh"
#include "DataBase/DataBase.hh"
#include "Boundary/Boundary.hh"
#include "Geometry/Dimension.hh"

#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Default constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
RadialViscosity<Dimension>::
RadialViscosity():
  MonaghanGingoldViscosity<Dimension>() {
}

//------------------------------------------------------------------------------
// Construct with the given value for the linear and quadratic coefficients.
//------------------------------------------------------------------------------
template<typename Dimension>
RadialViscosity<Dimension>::
RadialViscosity(Scalar Clinear, Scalar Cquadratic):
  MonaghanGingoldViscosity<Dimension>(Clinear, Cquadratic) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
RadialViscosity<Dimension>::
~RadialViscosity() {
}

//------------------------------------------------------------------------------
// Determine the contribution that mimics the artificial internal energy.
// Just call the MonaghanGingold base method if Q should be active.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
RadialViscosity<Dimension>::
viscousInternalEnergy(const NodeIDIterator<Dimension>& nodeI,
                      const NodeIDIterator<Dimension>& nodeJ,
		      const Vector& rij, const Vector& vij,
                      const Vector& etai, const Vector& etaj,
                      const Scalar ci, const Scalar cj) const {

  const Vector ri = nodeI.nodeListPtr()->positions()(nodeI);
  const Vector rj = nodeJ.nodeListPtr()->positions()(nodeJ);
  const Vector riNorm = ri.unitVector();
  const Vector rjiNorm = (rj - ri).unitVector();
  const double scale = pow(riNorm.dot(rjiNorm), 10);
  return scale*MonaghanGingoldViscosity<Dimension>::viscousInternalEnergy(nodeI,
                                                                          nodeJ,
                                                                          rij,
                                                                          vij,
                                                                          etai,
                                                                          etaj,
                                                                          ci,
                                                                          cj);
}

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
template class RadialViscosity< Dim<1> >;
template class RadialViscosity< Dim<2> >;
template class RadialViscosity< Dim<3> >;
}
