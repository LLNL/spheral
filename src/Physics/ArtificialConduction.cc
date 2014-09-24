//---------------------------------Spheral++----------------------------------//
// ArtificialConduction -- Artificial smoothing of energy discontinuities
//
//
// Created by CDR, 9/24/2014
//----------------------------------------------------------------------------//

#include "ArtificialConduction.hh"
#include "Field/Field.hh"
#include "Hydro/HydroFieldNames.hh"

namespace Spheral {
namespace PhysicsSpace {

using DataBaseSpace::DataBase;
using FieldSpace::Field;

//------------------------------------------------------------------------------
// Default constructor
//------------------------------------------------------------------------------
template<typename Dimension>
ArtificialConduction<Dimension>::
ArtificialConduction(const bool gradPMode):
    Physics<Dimension>(),
    mDepsDt(FieldSpace::Copy),
    mGradPMode(gradPMode){
    
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
ArtificialConduction<Dimension>::~ArtificialConduction() {
}
    
//------------------------------------------------------------------------------
// Accessor Fns
//------------------------------------------------------------------------------
template<typename Dimension>
void
ArtificialConduction<Dimension>::
gradPMode(bool val)
{
    mGradPMode = val;
}
template<typename Dimension>
bool
ArtificialConduction<Dimension>::
gradPMode() const{ return mGradPMode;}

template<typename Dimension>
typename Dimension::Scalar
ArtificialConduction<Dimension>::
vsig() const{ return mVsig;}

template<typename Dimension>
const FieldList<Dimension, typename Dimension::Vector>&
ArtificialConduction<Dimension>::
DepsDt() const {
    return mDepsDt;
}

//------------------------------------------------------------------------------
// Create and register the conduction
//------------------------------------------------------------------------------
template<typename Dimension>
void
ArtificialConduction<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {
    typedef typename State<Dimension>::PolicyPointer PolicyPointer;


}





}
}