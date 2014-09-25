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
ArtificialConduction(const Scalar alphaArCond):
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
    mVsig = dataBase.newFluidFieldList(0.0, "Artificial Conduction vsig");
    mGradP = dataBase.newFluidFieldList(0.0, "Pressure Gradient");
}
    //register derivs for gradP
//------------------------------------------------------------------------------
// Register gradP
//------------------------------------------------------------------------------
template<typename Dimension>
void
SPHHydroBase<Dimension>::
registerDerivatives(DataBase<Dimension>& dataBase,
                    StateDerivatives<Dimension>& derivs) {
    derivs.enroll(mGradP);
    // do i want to do this?? is it sufficient merely to calculate it below in eval?
}

    
//------------------------------------------------------------------------------
// Meat and potatoes
//------------------------------------------------------------------------------
template<typename Dimension>
void
SPHHydroBase<Dimension>::
evaluateDerivatives(const typename Dimension::Scalar time,
                    const typename Dimension::Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivatives) const {
    
    const TableKernel<Dimension>& W = this->kernel();
    
    bool CSPHisOn = 0;
    
    // The connectivity map
    const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
    const vector<const NodeList<Dimension>*>& nodeLists = connectivityMap.nodeLists();
    const size_t numNodeLists = nodeLists.size();
    
    // The relevant fields
    const FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
    const FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero);
    const FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
    const FieldList<Dimension, Scalar> specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
    const FieldList<Dimension, SymTensor> H = state.fields(HydroFieldNames::H, SymTensor::zero);
    const FieldList<Dimension, Scalar> pressure = state.fields(HydroFieldNames::pressure, 0.0);
    CHECK(mass.size() == numNodeLists);
    CHECK(position.size() == numNodeLists);
    CHECK(massDensity.size() == numNodeLists);
    CHECK(specificThermalEnergy.size() == numNodeLists);
    CHECK(H.size() == numNodeLists);
    CHECK(pressure.size() == numNodeLists);
    
    // The relevant derivatives
    FieldList<Dimension, Scalar> DepsDt = derivatives.fields(IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::specificThermalEnergy, 0.0);
    FieldList<Dimension, Scalar> gradP = derivatives.fields("Pressure Gradient", 0.0);
    CHECK(DepsDt.size() == numNodeLists);
    CHECK(gradP.size() == numNodeLists);
    
    // Now check if CSPH is active
    if (state.registered(A_CSPH))
    {
        CSPHisOn = 1;
        const FieldList<Dimension, Scalar> A = state.fields(HydroFieldNames::A_CSPH, 0.0);
        const FieldList<Dimension, Vector> B = state.fields(HydroFieldNames::B_CSPH, Vector::zero);
        const FieldList<Dimension, Vector> C = state.fields(HydroFieldNames::C_CSPH, Vector::zero);
        const FieldList<Dimension, Tensor> D = state.fields(HydroFieldNames::D_CSPH, Tensor::zero);
        const FieldList<Dimension, Vector> gradA0 = state.fields(HydroFieldNames::gradA0_CSPH, Vector::zero);
        const FieldList<Dimension, Vector> gradA = state.fields(HydroFieldNames::gradA_CSPH, Vector::zero);
        const FieldList<Dimension, Tensor> gradB = state.fields(HydroFieldNames::gradB_CSPH, Tensor::zero);
        CHECK(A.size() == numNodeLists);
        CHECK(B.size() == numNodeLists);
        CHECK(C.size() == numNodeLists);
        CHECK(D.size() == numNodeLists);
        CHECK(gradA0.size() == numNodeLists);
        CHECK(gradA.size() == numNodeLists);
        CHECK(gradB.size() == numNodeLists);
        //gradP = gradientCSPH(pressure, position, mass, H, A, B, C, D, gradA, gradB, connectivityMap, W);
    }

}


}
}