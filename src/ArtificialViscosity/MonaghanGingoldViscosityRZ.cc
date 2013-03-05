//---------------------------------Spheral++----------------------------------//
// A simple form for the artificial viscosity due to Monaghan & Gingold.
//----------------------------------------------------------------------------//
#include "MonaghanGingoldViscosityRZ.hh"
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

namespace Spheral {
namespace ArtificialViscositySpace {

using namespace std;

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
MonaghanGingoldViscosityRZ::
MonaghanGingoldViscosityRZ(Scalar Clinear, Scalar Cquadratic):
  MonaghanGingoldViscosity<Dim<2> >(Clinear, Cquadratic) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
MonaghanGingoldViscosityRZ::
~MonaghanGingoldViscosityRZ() {
}

//------------------------------------------------------------------------------
// Method to apply the viscous acceleration, work, and pressure, to the derivatives
// all in one step (efficiency and all).
//------------------------------------------------------------------------------
void
MonaghanGingoldViscosityRZ::
viscousEffects(const DataBase<Dim<2> >& dataBase,
               const ConnectivityMap<Dim<2> >& connectivityMap,
               const State<Dim<2> >& state,
               StateDerivatives<Dim<2> >& derivatives) const {

  // Get the set of NodeLists.
  const vector<const NodeList<Dimension >*>& nodeLists = connectivityMap.nodeLists();
  const int numNodeLists = nodeLists.size();

  const double Cl = this->Cl();
  const double Cq = this->Cq();
  const double eps2 = this->epsilon2();
  const bool balsaraShearCorrection = this->balsaraShearCorrection();

  // Get the full FieldLists for the state.
  const FieldList<Dimension, Scalar> mass = state.scalarFields(HydroFieldNames::mass);
  const FieldList<Dimension, Vector> position = state.vectorFields(HydroFieldNames::position);
  const FieldList<Dimension, Vector> velocity = state.vectorFields(HydroFieldNames::velocity);
  const FieldList<Dimension, SymTensor> H = state.symTensorFields(HydroFieldNames::H);
  const FieldList<Dimension, Scalar> rho = state.scalarFields(HydroFieldNames::massDensity);
  const FieldList<Dimension, Scalar> soundSpeed = state.scalarFields(HydroFieldNames::soundSpeed);

  // Derivative FieldLists.
  FieldList<Dimension, Scalar> maxQPressure = derivatives.scalarFields(HydroFieldNames::maxViscousPressure);
  FieldList<Dimension, Vector> DvDt = derivatives.vectorFields(IncrementState<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::velocity);
  FieldList<Dimension, Scalar> DepsDt = derivatives.scalarFields(IncrementState<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::specificThermalEnergy);

  // Prepare an empty FieldList to hold the pair-wise accelerations if we're using
  // the compatible energy discretization.
  FieldList<Dimension, vector<Vector> > QpairAccelerations = dataBase.newFluidFieldList(vector<Vector>(), "Q accelerations");

  // Iterate over the FluidNodeLists.
  int nodeListi = 0;
  for (DataBase<Dimension>::ConstFluidNodeListIterator itr = dataBase.fluidNodeListBegin();
       itr != dataBase.fluidNodeListEnd();
       ++itr, ++nodeListi) {
    const bool compatibleEnergyEvolution = (*itr)->compatibleEnergyEvolution();
    const TableKernel<Dimension>& W = (*itr)->Qkernel();

    // Grab the state for the NodeList we're working on.
    const Field<Dimension, Scalar>& massThis = *(mass[nodeListi]);
    const Field<Dimension, Vector>& positionThis = *(position[nodeListi]);
    const Field<Dimension, Vector>& velocityThis = *(velocity[nodeListi]);
    const Field<Dimension, SymTensor>& HThis = *(H[nodeListi]);
    const Field<Dimension, Scalar>& rhoThis = *(rho[nodeListi]);
    const Field<Dimension, Scalar>& soundSpeedThis = *(soundSpeed[nodeListi]);

    // Now grab the derivative fields for this NodeList.
    Field<Dimension, Scalar>& maxQPressureThis = *maxQPressure[nodeListi];
    Field<Dimension, Vector>& DvDtThis = *DvDt[nodeListi];
    Field<Dimension, Scalar>& DepsDtThis = *DepsDt[nodeListi];
    Field<Dimension, vector<Vector> >& QpairAccelerationsThis = *QpairAccelerations[nodeListi];

    // Iterate over the internal nodes.
    const int numNodes = NodeListRegistrar<Dimension>::instance().domainDecompositionIndependent() ? (*itr)->numNodes() : (*itr)->numInternalNodes();
    for (int i = 0; i != numNodes; ++i) {

      // State of node I.
      const Scalar& mi = massThis(i);
      const Vector& xi = positionThis(i);
      const Vector& vi = velocityThis(i);
      const SymTensor& Hi = HThis(i);
      const Scalar& rhoi = rhoThis(i);
      const Scalar& ci = soundSpeedThis(i);
      const Scalar Hdeti = Hi.Determinant();
      const Scalar ri = xi.x();
      const Scalar zi = xi.y();
      CHECK(ri >= 0.0);

      Scalar& maxQPressurei = maxQPressureThis(i);
      Vector& DvDti = DvDtThis(i);
      Scalar& DepsDti = DepsDtThis(i);
      vector<Vector>& QpairAccelerationsi = QpairAccelerationsThis(i);

      // Connectivity for this node.
      const vector< vector<int> >& fullConnectivity = connectivityMap.connectivityForNode(*itr, i);
      CHECK(fullConnectivity.size() == numNodeLists);

      // Are we using the shear correction?
      Scalar fsheari = 1.0;
      if (balsaraShearCorrection) fsheari = (*(this->mShearMultiplier[nodeListi]))(i);
      CHECK(fsheari >= 0.0 && fsheari <= 1.0);

      // Iterate over the neighboring NodeLists.
      for (int nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {

        // Connectivity of this node with this NodeList.  We only need to proceed if
        // there are some nodes in this list.
        const vector<int>& connectivity = fullConnectivity[nodeListj];
        if (connectivity.size() > 0) {

          // Grab the state for the neighbor NodeList.
          const Field<Dimension, Scalar>& massThem = *(mass[nodeListj]);
          const Field<Dimension, Vector>& positionThem = *(position[nodeListj]);
          const Field<Dimension, Vector>& velocityThem = *(velocity[nodeListj]);
          const Field<Dimension, SymTensor>& HThem = *(H[nodeListj]);
          const Field<Dimension, Scalar>& rhoThem = *(rho[nodeListj]);
          const Field<Dimension, Scalar>& soundSpeedThem = *(soundSpeed[nodeListj]);

          Field<Dimension, Scalar>& maxQPressureThem = *maxQPressure[nodeListj];
          Field<Dimension, Vector>& DvDtThem = *DvDt[nodeListj];
          Field<Dimension, Scalar>& DepsDtThem = *DepsDt[nodeListj];
          Field<Dimension, vector<Vector> >& QpairAccelerationsThem = *QpairAccelerations[nodeListj];

          const int firstGhostNodej = nodeLists[nodeListj]->firstGhostNode();

          // Iterate over the neighbors in this NodeList.
          for (vector<int>::const_iterator jItr = connectivity.begin();
               jItr != connectivity.end();
               ++jItr) {
            const int j = *jItr;
            CHECK(j < nodeLists[nodeListj]->numNodes());

            // Only proceed if this node pair has not been calculated yet.
            if (connectivityMap.calculatePairInteraction(nodeListi, i, 
                                                         nodeListj, j,
                                                         firstGhostNodej)) {

              // State for node J.
              const Scalar& mj = massThem(j);
              const Vector& xj = positionThem(j);
              const Vector& vj = velocityThem(j);
              const SymTensor& Hj = HThem(j);
              const Scalar& rhoj = rhoThem(j);
              const Scalar& cj = soundSpeedThem(j);
              const Scalar Hdetj = Hj.Determinant();
              const Scalar rj = abs(xj.x());
              const Scalar zj = xj.y();
              CHECK(rj >= 0.0);

              Scalar& maxQPressurej = maxQPressureThem(j);
              Vector& DvDtj = DvDtThem(j);
              Scalar& DepsDtj = DepsDtThem(j);
              vector<Vector>& QpairAccelerationsj = QpairAccelerationsThem(j);

              // Are we using the shear correction?
              Scalar fshearj = 1.0;
              if (balsaraShearCorrection) fshearj = (*(this->mShearMultiplier[nodeListj]))(j);
              CHECK(fshearj >= 0.0 && fshearj <= 1.0);

              // Compute the displacement.
              const Vector xij = xi - xj;
              const Vector vij = vi - vj;

              const Vector etai = Hi*xij;
              const Vector etaj = Hj*xij;
              const Scalar etaMagi = etai.magnitude();
              const Scalar etaMagj = etaj.magnitude();

              // Symmetrized gradient.
              const Vector Hetai = Hi*etai.unitVector();
              const Vector gradWi = Hetai*W.grad(etaMagi, Hdeti);
        
              const Vector Hetaj = Hj*etaj.unitVector();
              const Vector gradWj = Hetaj*W.grad(etaMagj, Hdetj);

              const Vector gradWij = 0.5*(gradWi + gradWj);

              // The effective smoothing scale along this node pairs direction.
              const Scalar hi = sqrt(xij.magnitude2()/(etaMagi*etaMagi + 1.0e-30));
              const Scalar hj = sqrt(xij.magnitude2()/(etaMagj*etaMagj + 1.0e-30));
              const Scalar hij = 0.5*(hi + hj);

              const Scalar riInv = ri/(ri*ri + eps2*hi*hi);
              const Scalar rjInv = rj/(rj*rj + eps2*hj*hj);

        const Scalar etair = etai.x();
        const Scalar etajr = etaj.x();
        const Scalar etar12 = 0.5*(etair + etajr);
        const Scalar etairInv = etair/(etair*etair + 1.0e-30);
        const Scalar etajrInv = etajr/(etajr*etajr + 1.0e-30);
        const Scalar etar12Inv = etar12/(etar12*etar12 + 1.0e-30);

//           // Compute mu.
//           const Scalar rc = 0.5*(ri + rj);
//           const Scalar rij = ri - rj;
//           const Scalar zij = zi - zj;
//           const Scalar rcInv = rc/(rc*rc + eps2*hij*hij);
//           const Scalar rijInv = rij/(rij*rij + eps2*hij*hij);
//           const Scalar zijInv = zij/(zij*zij + eps2*hij*hij);
//           const Scalar thpt = rcInv*(ri*vi.x() - rj*vj.x())*rijInv + (vi.y() - vj.y())*zijInv;
//           const Scalar mui = min(0.0, hi*thpt);
//           const Scalar muj = min(0.0, hj*thpt);

//              // Compute mu.
//              const Scalar mui = min(0.0, (vij.dot(etai))/(etai.magnitude2() + eps2));
//              const Scalar muj = min(0.0, (vij.dot(etaj))/(etaj.magnitude2() + eps2));
//         const Scalar mui = min(0.0, (vij.x()*etai.x()*riInv + vij.y()*etai.y()*rjInv)/(etai.magnitude2() + this->epsilon2()));
//         const Scalar muj = min(0.0, (vij.x()*etaj.x()*riInv + vij.y()*etaj.y()*rjInv)/(etaj.magnitude2() + this->epsilon2()));

        const Scalar rc = 0.5*(ri + rj);
        const Scalar rij = ri - rj;
        const Scalar zij = zi - zj;
        const Scalar ack = this->epsilon2() * hij*hij;
        const Scalar thpt = (ri*vi.x() - rj*vj.x())*rij/(rij*rij + ack)*rc/(rc*rc + ack) + (vi.y() - vj.y())*zij/(zij*zij + ack);
        const Scalar mui = min(0.0, hi*thpt);
        const Scalar muj = min(0.0, hj*thpt);
        // The artificial internal energy.
        const Scalar ei = fsheari*(-Cl*ci + Cq*mui)*mui;
        const Scalar ej = fshearj*(-Cl*cj + Cq*muj)*muj;
        CHECK(ei >= 0.0);
        CHECK(ej >= 0.0);

              // Now we can compute the acceleration, work, and artificial pressure.
              const Scalar QPi = 0.5*(ei/rhoi + ej/rhoj);
              CHECK(QPi >= 0.0);
              const Scalar work = QPi*(vij.x()*riInv*gradWij.x() + vij.y()*rjInv*gradWij.y());
              const Vector acceleration = QPi*gradWij;
              const Scalar pressure = 0.5*rhoi*(ei + ej);
              if (!(work >= 0.0)) {
                cerr << work << " "
                     << i << " "
                     << j << " "
                     << xi << " "
                     << xj << " "
                     << vij << " "
                     << acceleration << " "
                     << pressure << " "
                     << ei << " "
                     << ej << " "
                     << mui << " "
                     << muj << " "
                     << vij.x()*riInv*gradWij.x() << " "
                     << vij.y()*rjInv*gradWij.y() << " "
                     << endl;
              }
              CHECK(work >= 0.0);

              // Apply them!
              DvDti -= mj*acceleration;
              DvDtj += mi*acceleration;
              if (compatibleEnergyEvolution) {
                QpairAccelerationsi.push_back(-mj*acceleration);
                QpairAccelerationsj.push_back(mi*acceleration);
              } else {
                DepsDti += 0.5*mj*work;
                DepsDtj += 0.5*mi*work;
              }
              maxQPressurei = max(maxQPressurei, pressure);
              maxQPressurej = max(maxQPressurej, pressure);

            }
          }
        }
      }
    }
  }

  // If we're doing compatible energy evolution, we have to add the viscous
  // pair-wise accelerations to the full pair-wise values.
  FieldList<Dimension, vector<Vector> > pairAccelerations = derivatives.vectorVectorFields(HydroFieldNames::pairAccelerations);
  nodeListi = 0;
  for (DataBase<Dimension>::ConstFluidNodeListIterator itr = dataBase.fluidNodeListBegin();
       itr != dataBase.fluidNodeListEnd();
       ++itr, ++nodeListi) {
    if ((*itr)->compatibleEnergyEvolution()) {
      for (int i = 0; i != (*itr)->numInternalNodes(); ++i) {
        vector<Vector>& pa = (*pairAccelerations[nodeListi])(i);
        const vector<Vector>& qpa = (*QpairAccelerations[nodeListi])(i);
        CHECK(pa.size() == qpa.size());
        for (int j = 0; j != pa.size(); ++j) pa[j] += qpa[j];
      }
    }
  }


}

}
}
