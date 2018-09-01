//---------------------------------Spheral++----------------------------------//
// A specialized form of the Artificial viscosity to restrict the viscous 
// interactions to the radial direction.
//----------------------------------------------------------------------------//
#include "MonaghanGingoldKurapatenkoViscosity.hh"
#include "Field/FieldList.hh"
#include "Material/EquationOfState.hh"
#include "DataBase/DataBase.hh"
#include "Boundary/Boundary.hh"
#include "Geometry/Dimension.hh"

#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Default constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
MonaghanGingoldKurapatenkoViscosity<Dimension>::
MonaghanGingoldKurapatenkoViscosity():
  MonaghanGingoldViscosity<Dimension>() {
}

//------------------------------------------------------------------------------
// Construct with the given value for the linear and quadratic coefficients.
//------------------------------------------------------------------------------
template<typename Dimension>
MonaghanGingoldKurapatenkoViscosity<Dimension>::
MonaghanGingoldKurapatenkoViscosity(Scalar Clinear, Scalar Cquadratic):
  MonaghanGingoldViscosity<Dimension>(Clinear, Cquadratic) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
MonaghanGingoldKurapatenkoViscosity<Dimension>::
~MonaghanGingoldKurapatenkoViscosity() {
}

//------------------------------------------------------------------------------
// Determine the contribution that mimics the artificial internal energy.
// Just call the MonaghanGingold base method if Q should be active.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
MonaghanGingoldKurapatenkoViscosity<Dimension>::
viscousInternalEnergy(const NodeIDIterator<Dimension>& nodeI,
                      const NodeIDIterator<Dimension>& nodeJ,
		      const Vector& rij, const Vector& vij,
                      const Vector& etai, const Vector& etaj,
                      const Scalar ci, const Scalar cj) const {

  const FluidNodeList<Dimension>* fluidPtr = nodeI.fluidNodeListPtr();
  CHECK(fluidPtr != 0);
  const EquationOfState<Dimension>& eosi = fluidPtr->equationOfState();
  const Scalar rhoi = fluidPtr->massDensity()(nodeI);
  const Scalar epsi = fluidPtr->specificThermalEnergy()(nodeI);
  const Scalar gamma = eosi.gamma(rhoi, epsi);
  const Scalar vijrel = max(0.0, -vij.dot(rij.unitVector()));
  const Scalar thpt = 0.25*(gamma + 1);
  return (Cq()*thpt*vijrel + 
	  sqrt(FastMath::square(Cq()*thpt*vijrel) + FastMath::square(Cl()*ci)))*vijrel;
}

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
template class MonaghanGingoldKurapatenkoViscosity< Dim<1> >;
template class MonaghanGingoldKurapatenkoViscosity< Dim<2> >;
template class MonaghanGingoldKurapatenkoViscosity< Dim<3> >;
}
