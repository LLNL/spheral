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
#include "FileIO/FileIO.hh"

namespace Spheral {
namespace ArtificialViscositySpace {
    
using namespace std;
using std::min;
using std::max;
using std::abs;
using FileIOSpace::FileIO;

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
                                    const Scalar nhQ,
                                    const Scalar nhL,
                                    const Scalar aMin,
                                    const Scalar aMax):
    Physics<Dimension>(),
    mDrvAlphaDtQ(FieldSpace::Copy),
    mDrvAlphaDtL(FieldSpace::Copy),
    mnhQ(nhQ),
    mnhL(nhL),
    maMin(aMin),
    maMax(aMax),
    myq(q),
    mRestart(DataOutput::registerWithRestart(*this)){
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
    nhQ(Scalar val)
    {
        mnhQ = val;
    }
    template<typename Dimension>
    void
    MorrisMonaghanReducingViscosity<Dimension>::
    nhL(Scalar val)
    {
        mnhL = val;
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
    nhQ() const{ return mnhQ;}
    
    template<typename Dimension>
    typename Dimension::Scalar
    MorrisMonaghanReducingViscosity<Dimension>::
    nhL() const{ return mnhL;}
    
    template<typename Dimension>
    const FieldList<Dimension, typename Dimension::Scalar>&
    MorrisMonaghanReducingViscosity<Dimension>::
    DrvAlphaDtQ() const{ return mDrvAlphaDtQ;}
    
    template<typename Dimension>
    const FieldList<Dimension, typename Dimension::Scalar>&
    MorrisMonaghanReducingViscosity<Dimension>::
    DrvAlphaDtL() const{ return mDrvAlphaDtL;}
    
//------------------------------------------------------------------------------
// On problem start up, we need to initialize our internal data.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MorrisMonaghanReducingViscosity<Dimension>::
initializeProblemStartup(DataBase<Dimension>& dataBase) {
    mDrvAlphaDtQ = dataBase.newFluidFieldList(0.0, IncrementBoundedFieldList<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::reducingViscosityMultiplierQ);
    mDrvAlphaDtL = dataBase.newFluidFieldList(0.0, IncrementBoundedFieldList<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::reducingViscosityMultiplierL);
    dataBase.resizeFluidFieldList(myq.reducingViscosityMultiplierQ(), 1.0, HydroFieldNames::reducingViscosityMultiplierQ);
    dataBase.resizeFluidFieldList(myq.reducingViscosityMultiplierL(), 1.0, HydroFieldNames::reducingViscosityMultiplierL);
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
        FieldList<Dimension, Scalar>& rvAlphaQ = myq.reducingViscosityMultiplierQ();
        FieldList<Dimension, Scalar>& rvAlphaL = myq.reducingViscosityMultiplierL();
        PolicyPointer reducingViscosityMultiplierPolicy(new IncrementBoundedFieldList<Dimension, Scalar>(maMin,maMax));
        state.enroll(rvAlphaQ, reducingViscosityMultiplierPolicy);
        state.enroll(rvAlphaL, reducingViscosityMultiplierPolicy);
    }
    
//------------------------------------------------------------------------------
// Register the state derivative fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MorrisMonaghanReducingViscosity<Dimension>::
registerDerivatives(DataBase<Dimension>& dataBase,
                    StateDerivatives<Dimension>& derivs) {
    derivs.enroll(mDrvAlphaDtQ);
    derivs.enroll(mDrvAlphaDtL);
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
    FieldList<Dimension, Scalar> reducingViscosityMultiplierQ = state.fields(HydroFieldNames::reducingViscosityMultiplierQ, 0.0);
    FieldList<Dimension, Scalar> reducingViscosityMultiplierL = state.fields(HydroFieldNames::reducingViscosityMultiplierL, 0.0);
    
    // Derivative FieldLists
    const FieldList<Dimension, Tensor> DvDx = derivs.fields(HydroFieldNames::velocityGradient, Tensor::zero);
    FieldList<Dimension, Scalar> DrvAlphaDtQ = derivs.fields(IncrementBoundedFieldList<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::reducingViscosityMultiplierQ, 0.0);
    FieldList<Dimension, Scalar> DrvAlphaDtL = derivs.fields(IncrementBoundedFieldList<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::reducingViscosityMultiplierL, 0.0);
    
    unsigned int numNodeLists = pressure.numFields();
    CHECK(soundSpeed.size() == numNodeLists);
    CHECK(reducingViscosityMultiplierQ.size() == numNodeLists);
    CHECK(reducingViscosityMultiplierL.size() == numNodeLists);
    CHECK(DrvAlphaDtQ.size() == numNodeLists);
    CHECK(DrvAlphaDtL.size() == numNodeLists);
    CHECK(DvDx.size() == numNodeLists);
    CHECK(massDensity.size() == numNodeLists);
    CHECK(Hsmooth.size() == numNodeLists);
        
    const Scalar ngs = myq.negligibleSoundSpeed();
    
    //Walk the nodes
    for (unsigned i = 0; i < numNodeLists; ++i){
        unsigned int numNodes = pressure[i]->numInternalElements();
        for (unsigned int j = 0; j < numNodes; j++){
            Scalar rvQ = reducingViscosityMultiplierQ(i,j);
            Scalar rvL = reducingViscosityMultiplierL(i,j);
            const Scalar cs = soundSpeed(i,j);
            const Scalar csSafe = (cs*cs + ngs*ngs)/cs;
            const Scalar pmin = max(ngs*ngs*massDensity(i,j),abs(pressure(i,j)));
            const Scalar source = max(-DvDx(i,j).Trace(),0.0);
            const Scalar adiabatIndex = max(cs*cs*massDensity(i,j)/pmin,1.0+ngs);
            const Scalar decayConstQ = (1.0/mnhQ)*sqrt((adiabatIndex-1.0)/(2.0*adiabatIndex));
            const Scalar decayConstL = (1.0/mnhL)*sqrt((adiabatIndex-1.0)/(2.0*adiabatIndex));
            const Scalar h = 1.0/(Dimension::rootnu(Hsmooth(i,j).Determinant()));
            const Scalar decayTimeQ = h/(decayConstQ*csSafe);
            const Scalar decayTimeL = h/(decayConstL*csSafe);
            
            // safeInverse instead of forcing cs != 0

            
            DrvAlphaDtQ(i,j) = (maMax-rvQ)*source - (rvQ - maMin)/decayTimeQ;
            DrvAlphaDtL(i,j) = (maMax-rvL)*source - (rvL - maMin)/decayTimeL;
            
            // shock detection switch

        }
    }
}

//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MorrisMonaghanReducingViscosity<Dimension>::
dumpState(FileIO& file, string pathName) const {
    file.write(mDrvAlphaDtQ, pathName + "/DrvAlphaDtQ");
    file.write(mDrvAlphaDtL, pathName + "/DrvAlphaDtL");
}
    
//------------------------------------------------------------------------------
// Restore the state from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MorrisMonaghanReducingViscosity<Dimension>::
restoreState(const FileIO& file, string pathName) {
    file.read(mDrvAlphaDtQ, pathName + "/DrvAlphaDtQ");
    file.read(mDrvAlphaDtL, pathName + "/DrvAlphaDtL");
}
    
//------------------------------------------------------------------------------
// Boundaries
//------------------------------------------------------------------------------
template<typename Dimension>
void
MorrisMonaghanReducingViscosity<Dimension>::
    applyGhostBoundaries(State<Dimension>& state,
                         StateDerivatives<Dimension>& derivs)
{
    FieldList<Dimension, Scalar> reducingViscosityMultiplierQ = state.fields(HydroFieldNames::reducingViscosityMultiplierQ, 0.0);
    FieldList<Dimension, Scalar> reducingViscosityMultiplierL = state.fields(HydroFieldNames::reducingViscosityMultiplierL, 0.0);
    
    for (ConstBoundaryIterator boundaryItr = this->boundaryBegin();
         boundaryItr != this->boundaryEnd();
         ++boundaryItr) {
        (*boundaryItr)->applyFieldListGhostBoundary(reducingViscosityMultiplierQ);
        (*boundaryItr)->applyFieldListGhostBoundary(reducingViscosityMultiplierL);
    }
}
    
template<typename Dimension>
void
MorrisMonaghanReducingViscosity<Dimension>::
enforceBoundaries(State<Dimension>& state,
                     StateDerivatives<Dimension>& derivs)
{
    FieldList<Dimension, Scalar> reducingViscosityMultiplierQ = state.fields(HydroFieldNames::reducingViscosityMultiplierQ, 0.0);
    FieldList<Dimension, Scalar> reducingViscosityMultiplierL = state.fields(HydroFieldNames::reducingViscosityMultiplierL, 0.0);
    
    for (ConstBoundaryIterator boundaryItr = this->boundaryBegin();
         boundaryItr != this->boundaryEnd();
         ++boundaryItr) {
        (*boundaryItr)->enforceFieldListBoundary(reducingViscosityMultiplierQ);
        (*boundaryItr)->enforceFieldListBoundary(reducingViscosityMultiplierL);
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
