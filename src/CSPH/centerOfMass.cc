//---------------------------------Spheral++------------------------------------
// Compute the center of mass of a FacetedVolume assuming a linear mass density
// field.
//------------------------------------------------------------------------------
#include "centerOfMass.hh"
#include "Geometry/Dimension.hh"
#include "Utilities/DBC.hh"
#include "Utilities/FastMath.hh"

namespace Spheral {
namespace CSPHSpace {

using FastMath::pow3;

//------------------------------------------------------------------------------
// 1D
//------------------------------------------------------------------------------
GeomVector<1>
centerOfMass(const Dim<1>::FacetedVolume& polyvol,
             const Dim<1>::Vector& gradRhoi) {
  typedef Dim<1>::Vector Vector;
  const double grhoi = gradRhoi.x(),
               x0 = polyvol.xmin().x(),
               x1 = polyvol.xmax().x(),
               x01 = x1 - x0;
  const double rho0 = (grhoi > 0.0 ? 
                       1.0 :
                       1.0 - grhoi*x01);
  const double num = 0.5*(rho0 - grhoi*x0)*(x1*x1 - x0*x0) + grhoi*(pow3(x1) - pow3(x0))/3.0;
  const double den = (rho0 - x0*grhoi)*x01 + 0.5*grhoi*(x1*x1 - x0*x0);
  CHECK2(den > 0.0, den << " " << x0 << " " << x1 << " " << grhoi);
  const double result = num/den;
  ENSURE(result >= x0 and result <= x1);
  return Vector(result);
}

//------------------------------------------------------------------------------
// 2D
//------------------------------------------------------------------------------
Dim<2>::Vector
centerOfMass(const Dim<2>::FacetedVolume& polyvol,
             const Dim<2>::Vector& gradRhoi) {
  VERIFY2("Implement me!", false);
}

//------------------------------------------------------------------------------
// 3D
//------------------------------------------------------------------------------
Dim<3>::Vector
centerOfMass(const Dim<3>::FacetedVolume& polyvol,
             const Dim<3>::Vector& gradRhoi) {
  VERIFY2("Implement me!", false);
}

}
}

