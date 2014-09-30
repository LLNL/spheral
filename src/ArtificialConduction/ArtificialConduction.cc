//---------------------------------Spheral++----------------------------------//
// ArtificialConduction -- Artificial smoothing of energy discontinuities
//
//
// Created by CDR, 9/24/2014
//----------------------------------------------------------------------------//
#include <stdio.h>
#include "ArtificialConduction.hh"
#include "Field/Field.hh"
#include "Hydro/HydroFieldNames.hh"
#include "FieldOperations/FieldListFunctions.hh"
#include "CSPH/gradientCSPH.hh"
#include "DataBase/IncrementFieldList.hh"
#include "Neighbor/ConnectivityMap.hh"

namespace Spheral {
namespace PhysicsSpace {

using std::vector;
using DataBaseSpace::DataBase;
using FieldSpace::Field;
using KernelSpace::TableKernel;
using NodeSpace::NodeList;
using NodeSpace::FluidNodeList;
using FieldSpace::FieldList;
using NeighborSpace::ConnectivityMap;
    
using CSPHSpace::gradientCSPH;
using FieldSpace::gradient;

//------------------------------------------------------------------------------
// Default constructor
//------------------------------------------------------------------------------
template<typename Dimension>
ArtificialConduction<Dimension>::
ArtificialConduction(const TableKernel<Dimension>& W,
                     const Scalar alphaArCond):
    Physics<Dimension>(),
    mKernel(W),
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
    
//------------------------------------------------------------------------------
// On problem start up, we need to initialize our internal data.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ArtificialConduction<Dimension>::
initializeProblemStartup(DataBase<Dimension>& dataBase) {
    mGradP = dataBase.newFluidFieldList(Vector::zero, "Pressure Gradient");
    mDepsDtArty = dataBase.newFluidFieldList(Scalar, "Artificial Cond. DepsDt");
}

//------------------------------------------------------------------------------
// Register the state (no op here)
//------------------------------------------------------------------------------
template<typename Dimension>
void
ArtificialConduction<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {
    typedef typename State<Dimension>::PolicyPointer PolicyPointer;
    
    // get the eps policy
    FieldList<Dimension, Scalar> specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
    PolicyPointer energyPolicy = state.policy(specificThermalEnergy.name()); /* this needs to be the key */
    PolicyPointer artificialConductionPolicy(new ArtificialConductionPolicy<Dimension, Scalar>(energyPolicy));
    state.enroll(specificThermalEnergy, artificialConductionPolicy);
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
    derivs.enroll(mDepsDtArty);
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
    
    const TableKernel<Dimension>& W = mKernel;
    
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
    FieldList<Dimension, Scalar> DepsDt = derivatives.fields("Artificial Cond. DepsDt", 0.0);
    FieldList<Dimension, Vector> gradP = derivatives.fields("Pressure Gradient", Vector::zero);
    CHECK(DepsDt.size() == numNodeLists);
    CHECK(gradP.size() == numNodeLists);
    
    // Now check if CSPH is active
    if (state.registered(HydroFieldNames::A_CSPH))
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
        
        gradP = gradientCSPH(pressure, position, mass, H, A, B, C, D, gradA, gradB, connectivityMap, W);
    }
    else { gradP = gradient(pressure,position,mass,mass,massDensity,H,W); }
    
    
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
            const SymTensor& Hi = H(nodeListi, i);
            CHECK(mi > 0.0);
            CHECK(rhoi > 0.0);
            
            Scalar& DepsDti = DepsDt(nodeListi, i);
            Vector& gradPi = gradP(nodeListi, i);
                      
            
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
                        const int firstGhostNodej = nodeLists[nodeListj]->firstGhostNode();
                        
                        // Only proceed if this node pair has not been calculated yet.
                        if (connectivityMap.calculatePairInteraction(nodeListi, i,
                                                                     nodeListj, j,
                                                                     firstGhostNodej)) {
                            
                            
                            // get the state for node j
                            const Vector& rj        = position(nodeListj, j);
                            const Scalar& mj        = mass(nodeListj, j);
                            const Scalar& rhoj      = massDensity(nodeListj, j);
                            const Scalar& epsj      = specificThermalEnergy(nodeListj, j);
                            const Scalar& Pj        = pressure(nodeListj, j);
                            const SymTensor& Hj     = H(nodeListj, j);
                            CHECK(mj > 0.0);
                            CHECK(rhoj > 0.0);
                            
                            Scalar& DepsDtj         = DepsDt(nodeListj, j);
                            Vector& gradPj          = gradP(nodeListj, j);
                            
                            // get some differentials
                            const Vector rij        = ri - rj;
                            const Vector rji        = rj - ri;
                            const Vector etai       = Hi*rij;
                            const Vector etaj       = Hj*rij;
                            const Vector etaiNorm   = etai.unitVector();
                            const Vector etajNorm   = etaj.unitVector();
                            const Scalar rhoij      = 0.5 * (rhoi + rhoj);
                            const Scalar uij        = epsi - epsj;
                            const Scalar Pij        = Pi - Pj;
                            const Scalar DPij       = 0.5 * (gradPi.dot(rji) -
                                                             gradPj.dot(rij));
                            
                            // start a-calculatin' all the things
                            const Scalar deltaPij   = min(fabs(Pij),fabs(Pij+DPij));
                            
                            const Scalar vsigij     = sqrt(deltaPij/rhoij);
                            //const Scalar vsigij     = sqrt(fabs(Pij)/rhoij); /* this is purely for testing against Price */
                            const Vector gradWij    = 0.5*(Hi*etaiNorm*W.grad(etai, Hi) +
                                                           Hj*etajNorm*W.grad(etaj, Hj));
                            
                            const Scalar deltaU     = (mj/rhoij) * (mAlphaArCond) * vsigij * uij * rij.dot(gradWij)/rij.magnitude();
                            
                            DepsDti += deltaU;
                            DepsDtj += -deltaU;
                            //if(i==50 || j==50) printf("%02d->%02d %3.2e: rhoij=%3.2e mj=%3.2e vsigij=%3.2e uij=%3.2e gradWij=%3.2e dU/U=%3.2e DuDt=%3.2e\n",j,i,(i==50? deltaU : -deltaU),rhoij,mj,vsigij,uij,gradWij.magnitude(),deltaU/epsi,DepsDti);
                        }
                    }
                }
            }
        }
    }
}
    
//------------------------------------------------------------------------------
// Vote on a time step.
//------------------------------------------------------------------------------
template<typename Dimension>
typename ArtificialConduction<Dimension>::TimeStepType
ArtificialConduction<Dimension>::
dt(const DataBase<Dimension>& dataBase,
   const State<Dimension>& state,
   const StateDerivatives<Dimension>& derivs,
   const Scalar currentTime) const {
    return TimeStepType(1.0e100, "Rate of artificial conduction change -- NO VOTE.");
}


}
}
