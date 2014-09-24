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
ArtificialConduction(const bool gradPMode,
                     const Scalar alphaArCond):
    Physics<Dimension>(),
    mDepsDt(FieldSpace::Copy),
    mGradPMode(gradPMode)
    mAlphaArCond(alphaArCond){
    
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
gradPMode(bool val) { mGradPMode = val;}

template<typename Dimension>
bool
ArtificialConduction<Dimension>::
gradPMode() const { return mGradPMode;}

template<typename Dimension>
const FieldList<Dimension, typename Dimension::Scalar>&
ArtificialConduction<Dimension>::
vsig() const { return mVsig;}

template<typename Dimension>
const FieldList<Dimension, typename Dimension::Scalar>&
ArtificialConduction<Dimension>::
DepsDt() const { return mDepsDt;}
    
//------------------------------------------------------------------------------
// On problem start up, we need to initialize our internal data.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ArtificialConduction<Dimension>::
initializeProblemStartup(DataBase<Dimension>& dataBase) {
    mVsig = dataBase.newFluidFieldList(0.0, FieldList<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::artiCondSignalVelocity);
    
    // I don't want an incrementBoundedFieldList here. just a fieldlist
}
    
//------------------------------------------------------------------------------
// Meat and potatoes
//------------------------------------------------------------------------------
template<typename Dimension>
void
ArtificialConduction<Dimension>::
computeConduction(const DataBase<Dimension>& dataBase,
                  const State<Dimension>& state) {
    /* PSEUDO CODE ---------
     
    grab eps
    grab P
    grab rho
    grab m
    
    if !gradPMode
        vsig = sqrt(abs(Pi-Pj)/(0.5*(rhoi+rhoj)))
    else
        vsig = sqrt(abs(gradP*dr)/(0.5*(rhoi+rhoj)))
     
     DepsDti += mj/rhoj * (0.5 * alphaArCond * vsig * (uj-ui) * dotp(rij,gradWij)
    
    */
    
    const FieldList<Dimension, Scalar> pressure = state.fields(HydroFieldNames::pressure, 0.0);
    const FieldList<Dimension, SymTensor> Hsmooth = state.fields(HydroFieldNames::H, SymTensor::zero);
    const FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
    const FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
    
    if (!mGradPMode)
    {
        
    }
    else
    {
        
    }
}


}
}