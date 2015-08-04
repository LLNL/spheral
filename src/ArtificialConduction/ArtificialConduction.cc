//---------------------------------Spheral++----------------------------------//
// ArtificialConduction -- Artificial smoothing of energy discontinuities
//
//
// Created by CDR, 9/24/2014
//----------------------------------------------------------------------------//
#include <stdio.h>
#include "ArtificialConduction.hh"
#include "ArtificialConductionPolicy.hh"
#include "Field/Field.hh"
#include "Hydro/HydroFieldNames.hh"
#include "FieldOperations/FieldListFunctions.hh"
#include "CRKSPH/gradientCRKSPH.hh"
#include "DataBase/IncrementFieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Utilities/safeInv.hh"

namespace Spheral {
namespace PhysicsSpace {

using std::vector;
using std::min;
using std::max;
using std::abs;
using std::string;
using DataBaseSpace::DataBase;
using FieldSpace::Field;
using KernelSpace::TableKernel;
using NodeSpace::NodeList;
using NodeSpace::FluidNodeList;
using FieldSpace::FieldList;
using NeighborSpace::ConnectivityMap;
    
using CRKSPHSpace::gradientCRKSPH;
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
    mDepsDtArty = dataBase.newFluidFieldList(0.0, "Artificial Cond DepsDt");
    mVsigMax = dataBase.newFluidFieldList(0.0, "Maximum Artificial Cond Signal Speed");
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
    
    // get the eps policy and enroll the new one passing the old one as arg
    FieldList<Dimension, Scalar> specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
    PolicyPointer energyPolicy = state.policy(state.key(specificThermalEnergy)); /* this needs to be the key */
    PolicyPointer artificialConductionPolicy(new ArtificialConductionPolicy<Dimension>(energyPolicy));
    state.enroll(specificThermalEnergy, artificialConductionPolicy);
    state.enroll(mVsigMax);
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
    
    bool CRKSPHisOn = 0;
    
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
    FieldList<Dimension, Scalar> vsigMax = state.fields("Maximum Artificial Cond Signal Speed", 0.0);
    CHECK(mass.size() == numNodeLists);
    CHECK(position.size() == numNodeLists);
    CHECK(massDensity.size() == numNodeLists);
    CHECK(specificThermalEnergy.size() == numNodeLists);
    CHECK(H.size() == numNodeLists);
    CHECK(pressure.size() == numNodeLists);
    
    // The relevant derivatives
    FieldList<Dimension, Scalar> DepsDt = derivatives.fields("Artificial Cond DepsDt", 0.0);
    FieldList<Dimension, Vector> gradP = derivatives.fields("Pressure Gradient", Vector::zero);
    CHECK(DepsDt.size() == numNodeLists);
    CHECK(gradP.size() == numNodeLists);
    
    // Now check if CRKSPH is active
    if (state.registered(HydroFieldNames::A_CRKSPH))
    {
        CRKSPHisOn = 1;
        const FieldList<Dimension, Scalar> A = state.fields(HydroFieldNames::A_CRKSPH, 0.0);
        const FieldList<Dimension, Vector> B = state.fields(HydroFieldNames::B_CRKSPH, Vector::zero);
        const FieldList<Dimension, Vector> gradA0 = state.fields(HydroFieldNames::gradA0_CRKSPH, Vector::zero);
        const FieldList<Dimension, Vector> gradA = state.fields(HydroFieldNames::gradA_CRKSPH, Vector::zero);
        const FieldList<Dimension, Tensor> gradB = state.fields(HydroFieldNames::gradB_CRKSPH, Tensor::zero);
        CHECK(A.size() == numNodeLists);
        CHECK(B.size() == numNodeLists);
        CHECK(gradA0.size() == numNodeLists);
        CHECK(gradA.size() == numNodeLists);
        CHECK(gradB.size() == numNodeLists);
        
        gradP = gradientCRKSPH(pressure, position, mass, H, A, B, gradA, gradB, connectivityMap, W, NodeCoupling());
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
            const Vector& ri    = position(nodeListi, i);
            const Scalar& mi    = mass(nodeListi, i);
            const Scalar& rhoi  = massDensity(nodeListi, i);
            const Scalar& epsi  = specificThermalEnergy(nodeListi, i);
            const Scalar& Pi    = pressure(nodeListi, i);
            const SymTensor& Hi = H(nodeListi, i);
            const Scalar hhi    = Dimension::nDim/(Hi.Trace());
            
            CHECK(mi > 0.0);
            CHECK(rhoi > 0.0);
            
            Scalar& DepsDti = DepsDt(nodeListi, i);
            Vector& gradPi = gradP(nodeListi, i);
            Scalar& vsigi = vsigMax(nodeListi, i);
            
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
                            const Scalar hhj        = Dimension::nDim/(Hj.Trace());
                            const Scalar hAve       = 0.5*(hhi+hhj);
                            
                            CHECK(mj > 0.0);
                            CHECK(rhoj > 0.0);
                            
                            Scalar& DepsDtj         = DepsDt(nodeListj, j);
                            Vector& gradPj          = gradP(nodeListj, j);
                            Scalar& vsigj           = vsigMax(nodeListj, j);
                            
                            // get some differentials
                            const Vector rij        = ri - rj; /* this is sign flipped but it's ok! */
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
                            //const Scalar deltaPij   = min(fabs(Pij),fabs(Pij+DPij)); 
                            const Scalar deltaPij   = abs(Pij);
                            const Scalar vsigij     = sqrt(deltaPij/rhoij);
                            const Vector gradWij    = 0.5*(Hi*etaiNorm*W.grad(etai, Hi) +
                                                           Hj*etajNorm*W.grad(etaj, Hj));

                            // store max vsig back to i and j
                            vsigi = max(vsigi,vsigij);
                            vsigj = max(vsigj,vsigij);
                                        
                            // calc and add change in energy
                            const Scalar deltaU     = (mj/rhoij) * (mAlphaArCond) * vsigij * uij * rij.dot(gradWij) * safeInv(rij.magnitude(),0.01*hAve);
                        
                            //if (((abs(ri.magnitude()-0.5)<=0.006 || abs(rj.magnitude()-0.5)<=0.006)) && abs(uij)>0)
			    //printf("%02d->%02d %0.2d %3.2e: vsigij=%3.2e ui,j=(%3.2e,%3.2e,%3.2e) gradWij=%3.2e ri,j=(%3.2e,%3.2e,%3.2e) deltaPij=%3.2e\n",
			    //         j,i,firstGhostNodej,deltaU,vsigij,epsi,epsj,uij,gradWij.magnitude(),ri.magnitude(),rj.magnitude(),rij.magnitude(),deltaPij);

                                    
                            DepsDti += deltaU;
                            DepsDtj += -deltaU;


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
    
    
    // Get some useful fluid variables from the DataBase.
    const FieldList<Dimension, int> mask = state.fields(HydroFieldNames::timeStepMask, 1);
    const FieldList<Dimension, Scalar> rho = state.fields(HydroFieldNames::massDensity, 0.0);
    const FieldList<Dimension, SymTensor> H = state.fields(HydroFieldNames::H, SymTensor::zero);
    const FieldList<Dimension, Scalar> eps = state.fields(HydroFieldNames::specificThermalEnergy, Scalar());
    const FieldList<Dimension, Scalar> soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
    const FieldList<Dimension, Scalar> vsigMax = state.fields("Maximum Artificial Cond Signal Speed", 0.0);
    const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap(this->requireGhostConnectivity());
    const int numNodeLists = connectivityMap.nodeLists().size();

    // Initialize the return value to some impossibly high value.
    Scalar minDt = FLT_MAX;
    // Set up some history variables to track what set the minimum Dt.
    Scalar lastMinDt = minDt;
    int lastNodeID;
    string lastNodeListName, reason;
    Scalar lastNodeScale, lastEps, lastVsig, lastRho;
    
    // Loop over every fluid node.
    size_t nodeListi = 0;
    for (typename DataBase<Dimension>::ConstFluidNodeListIterator nodeListItr = dataBase.fluidNodeListBegin();
         nodeListItr != dataBase.fluidNodeListEnd();
         ++nodeListItr, ++nodeListi) {
        const FluidNodeList<Dimension>& fluidNodeList = **nodeListItr;
        const Scalar kernelExtent = fluidNodeList.neighbor().kernelExtent();
        CHECK(kernelExtent > 0.0);
        
        for (typename ConnectivityMap<Dimension>::const_iterator iItr = connectivityMap.begin(nodeListi);
             iItr != connectivityMap.end(nodeListi);
             ++iItr) {
            const int i = *iItr;
            
            // If this node is masked, don't worry about it.
            if (mask(nodeListi, i) == 1) {
                
                // Get this nodes minimum characteristic smoothing scale.
                CHECK2(H(nodeListi, i).Determinant() >  0.0,
                       "Bad H tensor : " << H(nodeListi, i) << " : " << fluidNodeList.name() << " " << i << " " << fluidNodeList.firstGhostNode());
                const Scalar nodeScale = 1.0/Dimension::rootnu(H(nodeListi, i).Determinant());

                // Artificial conduction effective signal speed.
                CHECK(rho(nodeListi, i) > 0.0);
                const Scalar csq = vsigMax(nodeListi, i);
                const double csqDt = 0.3 * nodeScale/(csq + FLT_MIN); //ideally this should be Cfl() instead of 0.3
                if (csqDt < minDt) {
                    minDt = csqDt;
                    reason = "artificial conduction signal velocity limit";
                }
                
                if (minDt < lastMinDt) {
                    lastMinDt = minDt;
                    lastNodeID = i;
                    lastNodeListName = fluidNodeList.name();
                    lastRho = rho(nodeListi, i);
                    lastEps = eps(nodeListi, i);
                    lastVsig = vsigMax(nodeListi, i);
                }
            }
        }
    }
    
    std::stringstream reasonStream;
    reasonStream << "mindt = " << minDt << std::endl
    << reason << std::endl
    << "  (nodeList, node) = (" << lastNodeListName << ", " << lastNodeID << ") | "
    << "  h = " << lastNodeScale << std::endl
    << "  vsig = " << lastVsig << std::endl
    << "  rho = " << lastRho << std::endl
    << "  eps = " << lastEps << std::endl
    << std::ends;
    
    // Now build the result.
    TimeStepType result(minDt, reasonStream.str());
    
    return result;
}


}
}
