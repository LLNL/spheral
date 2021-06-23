//---------------------------------Spheral++----------------------------------//
// Physics -- dampled linear spring contact model Spheral++
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



template<typename Dimension>
typename Dimension::Scalar
DampedLinearSpring<Dimension>::
timeStep(const DataBase<Dimension>& dataBase,
         const State<Dimension>& state,
         const StateDerivatives<Dimension>& derivs,
               typename Dimension::Scalar /*currentTime*/) const{
  const auto& mask = state.fields(HydroFieldNames::timeStepMask, 1);
  const auto& mass = state.fields(HydroFieldNames::mass, 0.0); 
  const auto& position = state.fields(HydroFieldNames::position, Vector::zero);
  const auto& velocity = state.fields(HydroFieldNames::velocity, Vector::zero);  
  const auto& angularVelocity = state.fields("angularVelocity", Vector::zero);  

  const auto& connectivityMap = dataBase.connectivityMap(this->requireGhostConnectivity(),
                                                         this->requireOverlapConnectivity());

  const float pi = 3.14159;
  const float c1 = 1.0;
  const float c2 = 1.0;

  auto minContactTime = 1e30;

  #pragma omp parallel
  {

    int i, j, nodeListi, nodeListj;

  #pragma omp for
    for (auto kk = 0u; kk < npairs; ++kk) {

      const auto mi = mass(nodeListi,i);
      const auto mj = mass(nodeListj,j);
      const auto mij = (mi*mj)/(mi+mj);

      const auto stiffness = c1/mij;
      const auto dissipation = c2/(2.0*mij);
      const auto contactFrequency = std::sqrt(stiffness - dissipation*dissipation);
      const auto contactTime = pi/contactFrequency;
      minContactTime = min(contactTime,minContactTime);
    }
  }
  return minContactTime;
};

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
