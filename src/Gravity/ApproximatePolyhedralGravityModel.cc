//---------------------------------Spheral++----------------------------------//
// Approximate Polhedral Model
///----------------------------------------------------------------------------//
//   Model used to calculate the gravitational fields of an irregular shaped 
//   body by apply numerical quadrature to a discetization of its surface.
//
//   Pearl, Hitt, A Fast Quadrature-Based Gravity Model for the Homogeneous 
//   Polyhedron, MNRAS, 492(1), 2020.
//----------------------------------------------------------------------------//
// J. M. Pearl 2022
//----------------------------------------------------------------------------//

#include "ApproximatePolyhedralGravityModel.hh"

namespace Spheral {

//--------------------------------------------------------------------------------------------
// Constructor 
//--------------------------------------------------------------------------------------------
ApproximatePolyhedralGravityModel::
ApproximatePolyhedralGravityModel(const GeomPolyhedron & poly, const double Mass, const double G):
  mNumQuadraturePoints(poly.facets().size()),
  mQuadraturePoints(poly.facetCentroids()),
  mValues(poly.facetAreaVectors()){
    
	  auto rho = Mass/poly.volume();

    mResolutions.reserve(this->numQuadraturePoints());

	  for(unsigned int i=0; i < this->numQuadraturePoints(); i++){
      Scalar A = mValues[i].magnitude(); 
      mResolutions[i] = std::sqrt(A);
      mValues[i] *= G*rho;
    }

}

//--------------------------------------------------------------------------------------------
// Destructor 
//--------------------------------------------------------------------------------------------
ApproximatePolyhedralGravityModel::
~ApproximatePolyhedralGravityModel() {
}

//--------------------------------------------------------------------------------------------
// acceleration
//--------------------------------------------------------------------------------------------
Dim<3>::Vector
ApproximatePolyhedralGravityModel::
acceleration(const Dim<3>::Vector& position) const {
  
  Dim<3>::Vector acceleration = Vector::zero;

  const std::vector<Vector>& GrhoAn = this-> values();
  const std::vector<Vector>& quadPoints = this->quadraturePoints();
  const std::vector<Scalar>& res = this->resolutions();

  for(unsigned int i=0; i < this->numQuadraturePoints(); i++){
    Scalar r = (position - quadPoints[i]).magnitude();
    acceleration -= GrhoAn[i]/std::max(r,res[i]);
  }

  return acceleration;
}

//--------------------------------------------------------------------------------------------
// potential
//--------------------------------------------------------------------------------------------
Dim<3>::Scalar
ApproximatePolyhedralGravityModel::
potential(const Dim<3>::Vector& position) const {
  
  Scalar potential = 0.0;

  const std::vector<Vector>& GrhoAn = this-> values();
  const std::vector<Vector>& quadPoints = this->quadraturePoints();
  const std::vector<Scalar>& res = this->resolutions();

  for(unsigned int i=0; i < this->numQuadraturePoints(); i++){
    Vector r = (position - quadPoints[i]);
    potential -= GrhoAn[i].dot(r)/std::max(r.magnitude(),0.001*res[i]);
  }

  return 0.5*potential;
}

}


