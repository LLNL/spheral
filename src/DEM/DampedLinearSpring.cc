//---------------------------------Spheral++----------------------------------//
// Physics -- root abstract class for all DEM contact models in Spheral++
//----------------------------------------------------------------------------//
#include "DampedLinearSpring.hh"


namespace Spheral {

//------------------------------------------------------------------------------
// Default constructor
//------------------------------------------------------------------------------
template<typename Dimension>
DampedLinearSpring<Dimension>::
DampedLinearSpring(Scalar YoungsModulus,
                   Scalar restitutionCoefficient):
                   mRestitutionCoefficient(restitutionCoefficient),
                   mYoungsModulus(YoungsModulus) {
      mBeta = 3.14159/std::log(restitutionCoefficient);
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
DampedLinearSpring<Dimension>::
~DampedLinearSpring() {}


//------------------------------------------------------------------------------
// get our acceleration and other things
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Vector
DampedLinearSpring<Dimension>::
force(const typename Dimension::Scalar mi, 
      const typename Dimension::Scalar mj,
      const typename Dimension::Vector ri, 
      const typename Dimension::Vector rj,
      const typename Dimension::Vector vi, 
      const typename Dimension::Vector vj,
      const typename Dimension::Scalar hi, 
      const typename Dimension::Scalar hj) const{

};

template<typename Dimension>
typename Dimension::Vector
DampedLinearSpring<Dimension>:: 
torque(const typename Dimension::Scalar mi, 
       const typename Dimension::Scalar mj,
       const typename Dimension::Vector ri, 
       const typename Dimension::Vector rj,
       const typename Dimension::Vector vi, 
       const typename Dimension::Vector vj,
       const typename Dimension::Scalar hi, 
       const typename Dimension::Scalar hj) const{

};


}
