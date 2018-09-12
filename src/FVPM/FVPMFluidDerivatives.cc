//---------------------------------Spheral++----------------------------------//
// FVPMFluidDerivatives -- The Finite Volume Point Method form of the the 
// fluid derivativies.
//----------------------------------------------------------------------------//
#include "FVPMFluidDerivatives.hh"
#include "NodeList/FluidNodeList.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "Material/EquationOfState.hh"
#include "Kernel/TableKernel.hh"
#include "Neighbor/Neighbor.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Geometry/EigenStruct.hh"
#include "Utilities/timingUtilities.hh"
#include "Utilities/safeInv.hh"

#include "Hydro/HydroFieldNames.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/ReplaceState.hh"

#include "NodeList/secondMomentUtilities.hh"

#include "Utilities/DBC.hh"

namespace Spheral {


//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
FVPMFluidDerivatives<Dimension>::
FVPMFluidDerivatives():
  FluidDerivativeProducer<Dimension>() {
}

//------------------------------------------------------------------------------
// Copy constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
FVPMFluidDerivatives<Dimension>::
FVPMFluidDerivatives(const FVPMFluidDerivatives& rhs):
  FluidDerivativeProducer<Dimension>(rhs) {
}

//------------------------------------------------------------------------------
// Assignment.
//------------------------------------------------------------------------------
template<typename Dimension>
FVPMFluidDerivatives<Dimension>&
FVPMFluidDerivatives<Dimension>::
operator=(const FVPMFluidDerivatives& rhs) {
  FluidDerivativeProducer<Dimension>::operator=(rhs);
  return *this;
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
FVPMFluidDerivatives<Dimension>::
~FVPMFluidDerivatives<Dimension>() {
}

//------------------------------------------------------------------------------
template<typename Dimension>
void
FVPMFluidDerivatives<Dimension>::
updateMassDensity(const State<Dimension>& state,
                  const TableKernel<Dimension>& W,
                  const ConnectivityMap<Dimension>& connectivityMap,
                  Field<Dimension, Scalar>& massDensity) const {
  // Nothing to do here for now.
  VERIFY2(false, "FVPMFluidDerivatives::updateMassDensity undefined!");
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template<typename Dimension>
void
FVPMFluidDerivatives<Dimension>::
updateWeight(const State<Dimension>& state,
             const TableKernel<Dimension>& W,
             const ConnectivityMap<Dimension>& connectivityMap,
             Field<Dimension, typename Dimension::Scalar>& weight) const {
  VERIFY2(false, "FVPMFluidDerivatives::updateWeight undefined!");
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template<typename Dimension>
void
FVPMFluidDerivatives<Dimension>::
updateCorrections(const State<Dimension>& state,
                  const TableKernel<Dimension>& W,
                  const ConnectivityMap<Dimension>& connectivityMap,
                  Field<Dimension, Scalar>& omegaGradh,
                  Field<Dimension, Scalar>& A,
                  Field<Dimension, Vector>& B,
                  Field<Dimension, Vector>& C,
                  Field<Dimension, Tensor>& D,
                  Field<Dimension, Vector>& gradA,
                  Field<Dimension, Tensor>& gradB,
                  const bool useGradhCorrections) const {
  VERIFY2(false, "FVPMFluidDerivatives::updateCorrections undefined!");
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template<typename Dimension>
void
FVPMFluidDerivatives<Dimension>::
calculateDerivatives(const typename Dimension::Scalar time,
                     const typename Dimension::Scalar dt,
                     const FluidNodeList<Dimension>& nodeList,
                     const typename Dimension::Scalar nPerh,
                     const bool XSPH,
                     const bool compatibleEnergyEvolution,
                     const typename Dimension::Scalar epsTensile,
                     const TableKernel<Dimension>& W,
                     const ConnectivityMap<Dimension>& connectivityMap,
                     const State<Dimension>& state,
                     StateDerivatives<Dimension>& derivatives) const 
{
  REQUIRE(dt >= 0.0);
  REQUIRE(nPerh > 0.0);
  REQUIRE(epsTensile >= 0.0);

  typedef typename Timing::Time Time;
  const vector<const NodeList<Dimension>*>& nodeLists = connectivityMap.nodeLists();
  const int numNodeLists = nodeLists.size();
  double kernelExtent = W.kernelExtent();

  // Get the state and derivative FieldLists.
  // State FieldLists.
  const FieldList<Dimension, Scalar> mass = state.scalarFields(HydroFieldNames::mass);
  const FieldList<Dimension, Vector> position = state.vectorFields(HydroFieldNames::position);
  const FieldList<Dimension, Vector> velocity = state.vectorFields(HydroFieldNames::velocity);
  const FieldList<Dimension, Scalar> massDensity = state.scalarFields(HydroFieldNames::massDensity);
  const FieldList<Dimension, Scalar> specificThermalEnergy = state.scalarFields(HydroFieldNames::specificThermalEnergy);
  const FieldList<Dimension, SymTensor> H = state.symTensorFields(HydroFieldNames::H);
  const FieldList<Dimension, Scalar> pressure = state.scalarFields(HydroFieldNames::pressure);
  const FieldList<Dimension, Scalar> soundSpeed = state.scalarFields(HydroFieldNames::soundSpeed);
  const FieldList<Dimension, Scalar> positionWeight = state.scalarFields(HydroFieldNames::positionWeight);
  CHECK(mass.size() == numNodeLists);
  CHECK(position.size() == numNodeLists);
  CHECK(velocity.size() == numNodeLists);
  CHECK(massDensity.size() == numNodeLists);
  CHECK(specificThermalEnergy.size() == numNodeLists);
  CHECK(H.size() == numNodeLists);
  CHECK(pressure.size() == numNodeLists);
  CHECK(soundSpeed.size() == numNodeLists);
  CHECK(positionWeight.size() == numNodeLists);

  // Derivative FieldLists.
  FieldList<Dimension, Vector> DxDt = derivatives.vectorFields(IncrementState<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::position);
  FieldList<Dimension, Scalar> DrhoDt = derivatives.scalarFields(IncrementState<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::massDensity);
  FieldList<Dimension, Vector> DvDt = derivatives.vectorFields(IncrementState<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::velocity);
  FieldList<Dimension, Scalar> DepsDt = derivatives.scalarFields(IncrementState<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::specificThermalEnergy);
  FieldList<Dimension, Tensor> DvDx = derivatives.tensorFields(HydroFieldNames::velocityGradient);
  FieldList<Dimension, Tensor> localDvDx = derivatives.tensorFields(HydroFieldNames::internalVelocityGradient);
  FieldList<Dimension, SymTensor> DHDt = derivatives.symTensorFields(IncrementState<Dimension, Field<Dimension, SymTensor> >::prefix() + HydroFieldNames::H);
  FieldList<Dimension, SymTensor> Hideal = derivatives.symTensorFields(ReplaceState<Dimension, Field<Dimension, SymTensor> >::prefix() + HydroFieldNames::H);
  FieldList<Dimension, vector<Vector> > pairAccelerations = derivatives.vectorVectorFields(HydroFieldNames::pairAccelerations);
  FieldList<Dimension, Scalar> weightedNeighborSum = derivatives.scalarFields(HydroFieldNames::weightedNeighborSum);
  FieldList<Dimension, SymTensor> massSecondMoment = derivatives.symTensorFields(HydroFieldNames::massSecondMoment);
  CHECK(rhoSum.size() == numNodeLists);
  CHECK(DxDt.size() == numNodeLists);
  CHECK(DrhoDt.size() == numNodeLists);
  CHECK(DvDt.size() == numNodeLists);
  CHECK(DepsDt.size() == numNodeLists);
  CHECK(DvDx.size() == numNodeLists);
  CHECK(localDvDx.size() == numNodeLists);
  CHECK(DHDt.size() == numNodeLists);
  CHECK(Hideal.size() == numNodeLists);
  CHECK(pairAccelerations.size() == numNodeLists);
  CHECK(weightedNeighborSum.size() == numNodeLists);
  CHECK(massSecondMoment.size() == numNodeLists);

  // Figure out which NodeList we are.
  const int nodeListi = distance(nodeLists.begin(), find(nodeLists.begin(), nodeLists.end(), &nodeList));
  CHECK(nodeListi < numNodeLists);
  const int firstGhostNodei = nodeLists[nodeListi]->firstGhostNode();

  // Get the work field for this NodeList.
  Field<Dimension, Scalar>& workFieldi = nodeLists[nodeListi]->work();

  // Iterate over the internal nodes in this NodeList.
  for (typename ConnectivityMap<Dimension>::const_iterator iItr = connectivityMap.begin(nodeListi);
       iItr != connectivityMap.end(nodeListi);
       ++iItr) 
  {
    const int i = *iItr;

    // Prepare to accumulate the time.
    const Time start = Timing::currentTime();
    size_t ncalc = 0;

    // Get the state for node i.
    const Vector& xi = position(nodeListi, i);
    const Vector& vi = velocity(nodeListi, i);
    const Scalar& rhoi = massDensity(nodeListi, i);
    const Scalar& epsi = specificThermalEnergy(nodeListi, i);
    const Scalar& Pi = pressure(nodeListi, i);
    const SymTensor& Hi = H(nodeListi, i);
    const Scalar& ci = soundSpeed(nodeListi, i);
    const Scalar Hdeti = Hi.Determinant();
    CHECK(rhoi > 0.0);
    CHECK(Hdeti > 0.0);

    Vector& DxDti = DxDt(nodeListi, i);
    Scalar& DrhoDti = DrhoDt(nodeListi, i);
    Vector& DvDti = DvDt(nodeListi, i);
    Scalar& DepsDti = DepsDt(nodeListi, i);
    Tensor& DvDxi = DvDx(nodeListi, i);
    Tensor& localDvDxi = localDvDx(nodeListi, i);
    SymTensor& DHDti = DHDt(nodeListi, i);
    SymTensor& Hideali = Hideal(nodeListi, i);
    vector<Vector>& pairAccelerationsi = pairAccelerations(nodeListi, i);
    Scalar& weightedNeighborSumi = weightedNeighborSum(nodeListi, i);
    SymTensor& massSecondMomenti = massSecondMoment(nodeListi, i);
    Scalar& worki = workFieldi(i);

    // Get the connectivity info for this node.
    const vector< vector<int> >& fullConnectivity = connectivityMap.connectivityForNode(&nodeList, i);

    // Assemble the positions {xk} and effective radii {rk} of all nodes 
    // within the support domain of i (including i itself).
    vector<Vector> xks;
    vector<double> rks;
    xks.push_back(xi);
    double detHi = Hi.det();
    rks.push_back(kernelExtent/Dimension::rootnu(detHi));
    for (int nodeListk = 0; nodeListk != numNodeLists; ++nodeListk)
    {
      for (int nodek = 0; nodek < fullConnectivity[nodeListk].size(); ++nodek)
      {
        int k = fullConnectivity[nodeListk][nodek];
        xks.push_back(position(nodeListk, k));
        double detHk = H(nodeListk, k).det();
        rks.push_back(kernelExtent/Dimension::rootnu(detHk));
      }
    }

    // Iterate over the NodeLists.
    for (int nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) 
    {
      // Connectivity of this node with this NodeList.  We only need to proceed if
      // there are some nodes in this list.
      const vector<int>& connectivity = fullConnectivity[nodeListj];
      if (connectivity.size() > 0) 
      {
        const int firstGhostNodej = nodeLists[nodeListj]->firstGhostNode();

        // Loop over the neighbors.
#pragma vector always
        for (vector<int>::const_iterator jItr = connectivity.begin();
             jItr != connectivity.end();
             ++jItr) 
        {
          const int j = *jItr;

          // Only proceed if this node pair has not been calculated yet.
          if (connectivityMap.calculatePairInteraction(nodeListi, i, 
                                                       nodeListj, j,
                                                       firstGhostNodej)) 
          {
            ++ncalc;

            // Get the state for node j
            const Vector& xj = position(nodeListj, j);
            const Scalar& mj = mass(nodeListj, j);
            const Vector& vj = velocity(nodeListj, j);
            const Scalar& rhoj = massDensity(nodeListj, j);
            const Scalar& epsj = specificThermalEnergy(nodeListj, j);
            const Scalar& Pj = pressure(nodeListj, j);
            const SymTensor& Hj = H(nodeListj, j);
            const Scalar& cj = soundSpeed(nodeListj, j);
            const Scalar& pwj = positionWeight(nodeListj, j);
            const Scalar Hdetj = Hj.Determinant();
            CHECK(mj > 0.0);
            CHECK(rhoj > 0.0);
            CHECK(Hdetj > 0.0);

            Vector& DxDtj = DxDt(nodeListj, j);
            Scalar& DrhoDtj = DrhoDt(nodeListj, j);
            Vector& DvDtj = DvDt(nodeListj, j);
            Scalar& DepsDtj = DepsDt(nodeListj, j);
            Tensor& DvDxj = DvDx(nodeListj, j);
            Tensor& localDvDxj = localDvDx(nodeListj, j);
            SymTensor& DHDtj = DHDt(nodeListj, j);
            SymTensor& Hidealj = Hideal(nodeListj, j);
            vector<Vector>& pairAccelerationsj = pairAccelerations(nodeListj, j);
            Scalar& weightedNeighborSumj = weightedNeighborSum(nodeListj, j);
            SymTensor& massSecondMomentj = massSecondMoment(nodeListj, j);

            // Compute the interaction vector beta_ij on the intersection 
            // of support domains i and j.
            Vector betaij = testFunction.interactionVector(xi, Hi, xj, Hj, xks, Hks);

            // Compute the approximate flux between the nodes i and j.
            // FIXME
          }
        }
      }
    }

    // Now treat the boundary nodes.
    // FIXME
  }
}
//------------------------------------------------------------------------------

}
}

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
namespace Spheral {
  template class FVPMFluidDerivatives< Dim<1> >;
  template class FVPMFluidDerivatives< Dim<2> >;
  template class FVPMFluidDerivatives< Dim<3> >;
}
