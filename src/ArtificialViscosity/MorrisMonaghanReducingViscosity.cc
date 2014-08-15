//---------------------------------Spheral++----------------------------------//
// A simple form of the reducing artificial viscosity from Morris & Monaghan.
//----------------------------------------------------------------------------//
#include "MorrisMonaghanReducingViscosity.hh"
#include "DataOutput/Restart.hh"
#include "Field/FieldList.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "NodeList/FluidNodeList.hh"
#include "Neighbor/Neighbor.hh"
#include "Material/EquationOfState.hh"
#include "Boundary/Boundary.hh"
#include "Hydro/HydroFieldNames.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/IncrementBoundedFieldList.hh"

namespace Spheral {
namespace ArtificialViscositySpace {
    
using namespace std;
using std::min;
using std::max;
using std::abs;

using PhysicsSpace::Physics;
using DataOutput::Restart;
using FieldSpace::Field;
using FieldSpace::FieldList;
using DataBaseSpace::DataBase;
using NodeSpace::NodeList;
using NodeSpace::FluidNodeList;

//------------------------------------------------------------------------------
// Construct with the given value for the linear and quadratic coefficients.
//------------------------------------------------------------------------------
template<typename Dimension>
MorrisMonaghanReducingViscosity<Dimension>::
    MorrisMonaghanReducingViscosity(ArtificialViscosity<Dimension>& q,
                                    const Scalar nh,
                                    const Scalar aMin,
                                    const Scalar aMax):
    Physics<Dimension>(),
    mDrvAlphaDt(FieldSpace::Copy),
    mnh(nh),
    maMin(aMin),
    maMax(aMax),
    myq(q){
        myq.reducingViscosityCorrection(true);
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
MorrisMonaghanReducingViscosity<Dimension>::
~MorrisMonaghanReducingViscosity() {
}

    
// Accessor Fns
    template<typename Dimension>
    void
    MorrisMonaghanReducingViscosity<Dimension>::
    aMin(Scalar val)
    {
        maMin = val;
    }
    template<typename Dimension>
    void
    MorrisMonaghanReducingViscosity<Dimension>::
    aMax(Scalar val)
    {
        maMax = val;
    }
    template<typename Dimension>
    void
    MorrisMonaghanReducingViscosity<Dimension>::
    nh(Scalar val)
    {
        mnh = val;
    }
    
    template<typename Dimension>
    typename Dimension::Scalar
    MorrisMonaghanReducingViscosity<Dimension>::
    aMin() const{ return maMin;}
    
    template<typename Dimension>
    typename Dimension::Scalar
    MorrisMonaghanReducingViscosity<Dimension>::
    aMax() const{ return maMax;}
    
    template<typename Dimension>
    typename Dimension::Scalar
    MorrisMonaghanReducingViscosity<Dimension>::
    nh() const{ return mnh;}
    
    template<typename Dimension>
    const FieldList<Dimension, typename Dimension::Scalar>&
    MorrisMonaghanReducingViscosity<Dimension>::
    DrvAlphaDt() const{ return mDrvAlphaDt;}
    
//------------------------------------------------------------------------------
// On problem start up, we need to initialize our internal data.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MorrisMonaghanReducingViscosity<Dimension>::
initializeProblemStartup(DataBase<Dimension>& dataBase) {
    mDrvAlphaDt = dataBase.newFluidFieldList(0.0, IncrementBoundedFieldList<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::reducingViscosityMultiplier);
    dataBase.resizeFluidFieldList(myq.reducingViscosityMultiplier(), 1.0, HydroFieldNames::reducingViscosityMultiplier);
}
    

//------------------------------------------------------------------------------
// Register the state we need/are going to evolve.
//------------------------------------------------------------------------------
template<typename Dimension>
void
    MorrisMonaghanReducingViscosity<Dimension>::
    registerState(DataBase<Dimension>& dataBase,
                  State<Dimension>& state) {
        typedef typename State<Dimension>::PolicyPointer PolicyPointer;
        
        // Create local storage for rvAlpha
        FieldList<Dimension, Scalar>& rvAlpha = myq.reducingViscosityMultiplier();
        PolicyPointer reducingViscosityMultiplierPolicy(new IncrementBoundedFieldList<Dimension, Scalar>(maMin,maMax));
        state.enroll(rvAlpha, reducingViscosityMultiplierPolicy);
    }
    
//------------------------------------------------------------------------------
// Register the state derivative fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MorrisMonaghanReducingViscosity<Dimension>::
registerDerivatives(DataBase<Dimension>& dataBase,
                    StateDerivatives<Dimension>& derivs) {
    derivs.enroll(mDrvAlphaDt);
}

//------------------------------------------------------------------------------
// Determine the principle derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MorrisMonaghanReducingViscosity<Dimension>::
    evaluateDerivatives(const typename Dimension::Scalar time,
                        const typename Dimension::Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivs) const {
    // Get Qtys for Derivs
    const FieldList<Dimension, Scalar> soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
    const FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
    const FieldList<Dimension, Scalar> pressure = state.fields(HydroFieldNames::pressure, 0.0);
    const FieldList<Dimension, SymTensor> Hsmooth = state.fields(HydroFieldNames::H, SymTensor::zero);
    FieldList<Dimension, Scalar> reducingViscosityMultiplier = state.fields(HydroFieldNames::reducingViscosityMultiplier, 0.0);
    
    // Derivative FieldLists
    const FieldList<Dimension, Tensor> DvDx = derivs.fields(HydroFieldNames::velocityGradient, Tensor::zero);
    FieldList<Dimension, Scalar> DrvAlphaDt = derivs.fields(IncrementBoundedFieldList<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::reducingViscosityMultiplier, 0.0);
    
    unsigned int numNodeLists = pressure.numFields();
    CHECK(soundSpeed.size() == numNodeLists);
    CHECK(reducingViscosityMultiplier.size() == numNodeLists);
    CHECK(DrvAlphaDt.size() == numNodeLists);
    CHECK(DvDx.size() == numNodeLists);
    CHECK(massDensity.size() == numNodeLists);
    CHECK(Hsmooth.size() == numNodeLists);
    
    //Walk the nodes
    for (unsigned i = 0; i < numNodeLists; ++i){
        unsigned int numNodes = pressure[i]->numInternalElements();
        for (unsigned int j = 0; j < numNodes; j++){
            Scalar rv = reducingViscosityMultiplier(i,j);
            const Scalar source = max(-DvDx(i,j).Trace(),0.0);
            const Scalar adiabatIndex = soundSpeed(i,j)*soundSpeed(i,j)*massDensity(i,j)/pressure(i,j);
            const Scalar decayConst = (1.0/mnh)*sqrt((adiabatIndex-1.0)/(2.0*adiabatIndex));
            const Scalar h = 1.0/(Dimension::rootnu(Hsmooth(i,j).Determinant()));
            const Scalar decayTime = h/(decayConst*soundSpeed(i,j));
            

            
            DrvAlphaDt(i,j) = (maMax-rv)*source - (rv - maMin)/decayTime;
            
            // shock detection switch
            if (DrvAlphaDt(i,j) > 0){
                //myq.reducingViscosityMultiplier(i,j) = maMax;
            }
        }
    }
}
    
//------------------------------------------------------------------------------
// Vote on a time step.
//------------------------------------------------------------------------------
template<typename Dimension>
typename MorrisMonaghanReducingViscosity<Dimension>::TimeStepType
MorrisMonaghanReducingViscosity<Dimension>::
dt(const DataBase<Dimension>& dataBase,
   const State<Dimension>& state,
   const StateDerivatives<Dimension>& derivs,
   const Scalar currentTime) const {
    return TimeStepType(1.0e100, "Rate of viscosity change -- NO VOTE.");
}

    
}
}
