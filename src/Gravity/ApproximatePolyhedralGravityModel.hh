//---------------------------------Spheral++----------------------------------//
// Analytic Polhedral Model
///----------------------------------------------------------------------------//
//   Werner, R.A., "The Gravitational Potential of a Homogeneous Polyhedron
//   or Don't Cut Corners," Celestial Mechanics and Dynamical Astronomy, Vol.
//   59, pp. 253-278, 1994.
//
//   Werner, R.A., Scheeres, D.J., Exterior "Gravitation of a Polyhedron
//   Dervied and Compared with Harmonic and Mascon Gravitation
//   Representations of Aeroid 4769 Castalia," Celestial Mechanics and
//   Dynamical Astronomy, Vol. 65, pp. 313–344, 1996.
//
//   Werner, R.A., "The Solid Angle Hidden in Polyhedron Gravitation
//   Formulations," Journal of Geodesy, Vol. 97, pp. 307-328, 2017.
//----------------------------------------------------------------------------//
// J. M. Pearl 2022
//----------------------------------------------------------------------------//
#ifndef __Spheral_approximatePolyhedralGravityModel__
#define __Spheral_approximatePolyhedralGravityModel__

#include "Geometry/Dimension.hh"
#include "Geometry/GeomPolyhedron.hh"

#include <vector>
#include <string>

namespace Spheral {

typedef typename Dim<3>::Scalar Scalar;
typedef typename Dim<3>::Vector Vector;

class ApproximatePolyhedralGravityModel{

  public:

    ApproximatePolyhedralGravityModel(const GeomPolyhedron & poly, const double Mass, const double G);
    ~ApproximatePolyhedralGravityModel();


    Scalar potential(const Vector& position) const;
    Vector acceleration(const Vector& position) const;

    unsigned int numQuadraturePoints() const;

    const std::vector<Vector>& quadraturePoints() const;
    const std::vector<Vector>& values() const;
    const std::vector<Scalar>& resolutions() const;

  private:
    unsigned int mNumQuadraturePoints;      
    std::vector<Vector> mQuadraturePoints;  // quadrature point positions
    std::vector<Vector> mValues;            // quadrature point area vectors
    std::vector<Scalar> mResolutions;       // quadrature point length scale
};

}

#include "ApproximatePolyhedralGravityModelInline.hh"

#endif
