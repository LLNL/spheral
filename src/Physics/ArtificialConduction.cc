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

//------------------------------------------------------------------------------
// Register gradP
//------------------------------------------------------------------------------
template<typename Dimension>
void
ArtificialConduction<Dimension>::
registerDerivatives(DataBase<Dimension>& dataBase,
                    StateDerivatives<Dimension>& derivs) {
    derivs.enroll(mGradP);
    // do i want to do this?? is it sufficient merely to calculate it below in eval?
}

//------------------------------------------------------------------------------
// Register vsig
//------------------------------------------------------------------------------
template<typename Dimension>
void
ArtificialConduction<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {
    typedef typename State<Dimension>::PolicyPointer PolicyPointer;
    
    dataBase.resizeFluidFieldList(mVsig, 0.0,"Artificial Conduction vsig", false);
    state.enroll(mVsig);
}
    
//------------------------------------------------------------------------------
// Meat and potatoes
//------------------------------------------------------------------------------
template<typename Dimension>
void
ArtificialConduction<Dimension>::
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
    const FieldList<Dimension, Scalar> vsig = state.fields("Artificial Conduction vsig", 0.0);
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
    
    // Start our big loop over all FluidNodeLists.
    size_t nodeListi = 0;
    for (typename DataBase<Dimension>::ConstFluidNodeListIterator itr = dataBase.fluidNodeListBegin();
         itr != dataBase.fluidNodeListEnd();
         ++itr, ++nodeListi) {
        const NodeList<Dimension>& nodeList = **itr;
        
        // Iterate over the internal nodes in this NodeList.
        for (typename ConnectivityMap<Dimension>::const_iterator iItr = connectivityMap.begin(nodeListi);
             iItr != connectivityMap.end(nodeListi);
             ++iItr) {
            const int i = *iItr;
            
            // Get the state for node i.
            const Vector& ri = position(nodeListi, i);
            const Scalar& mi = mass(nodeListi, i);
            const Scalar& rhoi = massDensity(nodeListi, i);
            const Scalar& epsi = specificThermalEnergy(nodeListi, i);
            const Scalar& Pi = pressure(nodeListi, i);
            const Scalar& vsigi = pressure(nodeListi, i);
            const SymTensor& Hi = H(nodeListi, i);
            if (CSPHisOn)
            {
            const Scalar& Ai = A(nodeListi, i);
            const Vector& Bi = B(nodeListi, i);
            const Vector& gradA0i = gradA0(nodeListi, i);
            const Vector& gradAi = gradA(nodeListi, i);
            const Tensor& gradBi = gradB(nodeListi, i);
            CHECK(Ai > 0.0);
            }
            CHECK(mi > 0.0);
            CHECK(rhoi > 0.0);
            
            Scalar& DepsDti = DepsDt(nodeListi, i);
            Scalar& gradPi = gradP(nodeListi, i);
            
            // Get the connectivity info for this node.
            const vector< vector<int> >& fullConnectivity = connectivityMap.connectivityForNode(&nodeList, i);
            
            // Iterate over the NodeLists.
            for (size_t nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {
                
                // Connectivity of this node with this NodeList.  We only need to proceed if
                // there are some nodes in this list.
                const vector<int>& connectivity = fullConnectivity[nodeListj];
                if (connectivity.size() > 0) {
                    
                    // Loop over the neighbors.
#pragma vector always
                    for (vector<int>::const_iterator jItr = connectivity.begin();
                         jItr != connectivity.end();
                         ++jItr) {
                        const int j = *jItr;
                        
                        // Only proceed if this node pair has not been calculated yet.
                        // not sure about this..... <--------------------------------
                        if (connectivityMap.calculatePairInteraction(nodeListi, i,
                                                                     nodeListj, j,
                                                                     firstGhostNodej)) {
                            // Get the state for node j
                            const Vector& rj = position(nodeListj, j);
                            const Scalar& mj = mass(nodeListj, j);
                            const Scalar& rhoj = massDensity(nodeListj, j);
                            const Scalar& epsj = specificThermalEnergy(nodeListj, j);
                            const Scalar& Pj = pressure(nodeListj, j);
                            const Scalar& vsigj = pressure(nodeListj, j);
                            const SymTensor& Hj = H(nodeListj, j);
                            if (CSPHisOn)
                            {
                            const Scalar& Aj = A(nodeListj, j);
                            const Vector& Bj = B(nodeListj, j);
                            const Vector& gradA0j = gradA0(nodeListj, j);
                            const Vector& gradAj = gradA(nodeListj, j);
                            const Tensor& gradBj = gradB(nodeListj, j);
                            CHECK(Aj > 0.0 or j >= firstGhostNodej);
                            }
                            CHECK(mj > 0.0);
                            CHECK(rhoj > 0.0);
                            
                            Scalar& DepsDtj = DepsDt(nodeListj, j);
                            Scalar& gradPj = gradP(nodeListj, j);
                            
                            // Node displacement.
                            const Vector rij = ri - rj;
                            
                        }
                    }
                }
            }
        }
    }
}


}
}