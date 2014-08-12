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
using NeighborSpace::Neighbor;
using Material::EquationOfState;
using BoundarySpace::Boundary;
using NeighborSpace::ConnectivityMap;
using KernelSpace::TableKernel;

//------------------------------------------------------------------------------
// Construct with the given value for the linear and quadratic coefficients.
//------------------------------------------------------------------------------
template<typename Dimension>
MorrisMonaghanReducingViscosity<Dimension>::
    MorrisMonaghanReducingViscosity(const Scalar nh,
                            const Scalar aStar,
                            const Scalar Clinear,
                            const Scalar Cquadratic,
                            const bool linearInExpansion,
                            const bool quadraticInExpansion):
    Physics<Dimension>(),
    ArtificialViscosity<Dimension>(Clinear, Cquadratic),
    mLinearInExpansion(linearInExpansion),
    mQuadraticInExpansion(quadraticInExpansion),
    mrvAlpha(FieldSpace::Copy),
    mDrvAlphaDt(FieldSpace::Copy),
    mnh(nh),
    maStar(aStar){
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
MorrisMonaghanReducingViscosity<Dimension>::
~MorrisMonaghanReducingViscosity() {
}

//------------------------------------------------------------------------------
// On problem start up, we need to initialize our internal data.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MorrisMonaghanReducingViscosity<Dimension>::
initializeProblemStartup(DataBase<Dimension>& dataBase) {
    mrvAlpha = dataBase.newFluidFieldList(1.0, HydroFieldNames::rvAlpha);
    mDrvAlphaDt = dataBase.newFluidFieldList(0.0, IncrementBoundedFieldList<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::rvAlpha);
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
        dataBase.resizeFluidFieldList(mrvAlpha, 1.0, HydroFieldNames::rvAlpha, false);
        PolicyPointer rvAlphaPolicy(new IncrementBoundedFieldList<Dimension, Scalar>(maStar,1.0));
        state.enroll(mrvAlpha, rvAlphaPolicy);
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
    const FieldList<Dimension, Scalar> Hsmooth = state.fields(HydroFieldNames::Hsmooth, 0.0);
    const FieldList<Dimension, Scalar> rvAlpha = state.fields(HydroFieldNames::rvAlpha, 0.0);
    
    // Derivative FieldLists
    const FieldList<Dimension, Tensor> DvDx = derivs.fields(HydroFieldNames::velocityGradient, Tensor::zero);
    FieldList<Dimension, Scalar> DrvAlphaDt = derivs.fields(HydroFieldNames::rvAlpha, 0.0);
    
    unsigned int numNodeLists = pressure.numFields();
    CHECK(soundSpeed.size() == numNodeLists);
    CHECK(rvAlpha.size() == numNodeLists);
    CHECK(DrvAlphaDt.size() == numNodeLists);
    CHECK(DvDx.size() == numNodeLists);
    CHECK(massDensity.size() == numNodeLists);
    CHECK(Hsmooth.size() == numNodeLists);
    
    //Walk the nodes
    for (unsigned i = 0; i < numNodeLists; i++){
        unsigned int numNodes = pressure[i]->numInternalElements();
        for (unsigned int j = 0; j < numNodes; j++){
            const Scalar source = max(-DvDx(i,j).Trace(),0.0);
            const Scalar adiabatIndex = soundSpeed(i,j)*soundSpeed(i,j)*massDensity(i,j)/pressure(i,j);
            const Scalar decayConst = (1.0/mnh)*sqrt((adiabatIndex-1.0)/(2.0*adiabatIndex));
            const Scalar decayTime = Hsmooth(i,j)/(decayConst*soundSpeed(i,j));
            
            DrvAlphaDt(i,j) = source - (rvAlpha(i,j) - maStar)/decayTime;
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
    
//------------------------------------------------------------------------------
// The required method to compute the artificial viscous P/rho^2.
//------------------------------------------------------------------------------
template<typename Dimension>
pair<typename Dimension::Tensor,
typename Dimension::Tensor>
MorrisMonaghanReducingViscosity<Dimension>::
Piij(const unsigned nodeListi, const unsigned i,
     const unsigned nodeListj, const unsigned j,
     const Vector& xi,
     const Vector& etai,
     const Vector& vi,
     const Scalar rhoi,
     const Scalar csi,
     const SymTensor& Hi,
     const Vector& xj,
     const Vector& etaj,
     const Vector& vj,
     const Scalar rhoj,
     const Scalar csj,
     const SymTensor& Hj) const {
    
    const double Cl = this->mClinear;
    const double Cq = this->mCquadratic;
    const double eps2 = this->mEpsilon2;
    const bool balsaraShearCorrection = this->mBalsaraShearCorrection;
    
    // Are we applying the shear corrections?
    const Vector vij = vi - vj;
    Scalar fshear = 1.0;
    if (balsaraShearCorrection) {
        fshear = abs(vij.unitVector().dot(vij.unitVector()));
    }
    
    // Compute mu.
    const Scalar mui = vij.dot(etai)/(etai.magnitude2() + eps2);
    const Scalar muj = vij.dot(etaj)/(etaj.magnitude2() + eps2);
    
    // The artificial internal energy.
    const Scalar ei = fshear*(-Cl*csi*(mLinearInExpansion    ? mui                : min(0.0, mui)) +
                              Cq     *(mQuadraticInExpansion ? -sgn(mui)*mui*mui  : FastMath::square(min(0.0, mui))));
    const Scalar ej = fshear*(-Cl*csj*(mLinearInExpansion    ? muj                : min(0.0, muj)) +
                              Cq     *(mQuadraticInExpansion ? -sgn(muj)*muj*muj  : FastMath::square(min(0.0, muj))));
    CHECK(ei >= 0.0 or (mLinearInExpansion or mQuadraticInExpansion));
    CHECK(ej >= 0.0 or (mLinearInExpansion or mQuadraticInExpansion));
    
    // Now compute the symmetrized artificial viscous pressure.
    return make_pair(ei/rhoi*Tensor::one * mrvAlpha(nodeListi,i),
                     ej/rhoj*Tensor::one * mrvAlpha(nodeListj,j));
}

//------------------------------------------------------------------------------
// linearInExpansion
//------------------------------------------------------------------------------
template<typename Dimension>
bool
MorrisMonaghanReducingViscosity<Dimension>::
linearInExpansion() const {
    return mLinearInExpansion;
}

template<typename Dimension>
void
MorrisMonaghanReducingViscosity<Dimension>::
linearInExpansion(const bool x) {
    mLinearInExpansion = x;
}

//------------------------------------------------------------------------------
// quadraticInExpansion
//------------------------------------------------------------------------------
template<typename Dimension>
bool
MorrisMonaghanReducingViscosity<Dimension>::
quadraticInExpansion() const {
    return mQuadraticInExpansion;
}

template<typename Dimension>
void
MorrisMonaghanReducingViscosity<Dimension>::
quadraticInExpansion(const bool x) {
    mQuadraticInExpansion = x;
}
    
}
}
