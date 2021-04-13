//---------------------------------Spheral++----------------------------------//
// Physics -- root abstract class for all DEM contact models in Spheral++
//----------------------------------------------------------------------------//
#include "ContactModelBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "DataBase/DataBase.hh"
#include "Field/FieldList.hh"
#include "Field/NodeIterators.hh"
#include "Boundary/Boundary.hh"
#include "Neighbor/ConnectivityMap.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Default constructor
//------------------------------------------------------------------------------
template<typename Dimension>
ContactModelBase<Dimension>::
ContactModelBase(){}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
ContactModelBase<Dimension>::
~ContactModelBase() {}

template<typename Dimension>
typename Dimension::Vector
ContactModelBase<Dimension>::
force(const typename Dimension::Scalar mi, 
      const typename Dimension::Scalar mj,
      const typename Dimension::Vector ri, 
      const typename Dimension::Vector rj,
      const typename Dimension::Vector vi, 
      const typename Dimension::Vector vj,
      const typename Dimension::Scalar hi, 
      const typename Dimension::Scalar hj) const {

}

template<typename Dimension>
typename Dimension::Vector
ContactModelBase<Dimension>::
torque(const typename Dimension::Scalar mi, 
       const typename Dimension::Scalar mj,
       const typename Dimension::Vector ri, 
       const typename Dimension::Vector rj,
       const typename Dimension::Vector vi, 
       const typename Dimension::Vector vj,
       const typename Dimension::Scalar hi, 
       const typename Dimension::Scalar hj) const{

}



}
