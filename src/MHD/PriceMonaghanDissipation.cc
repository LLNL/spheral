//---------------------------------Spheral++----------------------------------//
// A simple form for the artificial viscosity due to Monaghan & Gingold.
//----------------------------------------------------------------------------//
#include "MHD/PriceMonaghanDissipation.hh"
#include "MHD/MHDFieldNames.hh"
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

using std::abs;
using std::min;
using std::max;

//------------------------------------------------------------------------------
// Construct with the given value for the linear and quadratic coefficients.
//------------------------------------------------------------------------------
PriceMonaghanDissipation::
PriceMonaghanDissipation(Scalar alpha, 
                         Scalar alphaU,
                         Scalar alphaB,
                         Scalar beta, 
                         Scalar mu0):
  ArtificialViscosity<Dim<3> >(alpha, beta),
  mAlpha(alpha),
  mAlphaU(alphaU),
  mAlphaB(alphaB),
  mBeta(beta),
  mMu0(mu0),
  mMinDt(FLT_MAX) {
  setLimiter(false);
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
PriceMonaghanDissipation::
~PriceMonaghanDissipation() {
}

//------------------------------------------------------------------------------
// Apply the "viscous" effects to the derivatives.
//------------------------------------------------------------------------------
void
PriceMonaghanDissipation::
viscousEffects(const DataBase<Dim<3> >& dataBase,
               const ConnectivityMap<Dim<3> >& connectivityMap,
               const State<Dim<3> >& state,
               StateDerivatives<Dim<3> >& derivatives) const {

  // Get the set of NodeLists.
  const vector<const NodeList<Dim<3> >*>& nodeLists = connectivityMap.nodeLists();
  const int numNodeLists = nodeLists.size();

  const double Cl = this->Cl();
  const double Cq = this->Cq();
  const double eps2 = this->epsilon2();
  const bool balsaraShearCorrection = this->balsaraShearCorrection();

  // Get the full FieldLists for the state.
  const FieldList<Dim<3>, Scalar> mass = state.scalarFields(HydroFieldNames::mass);
  const FieldList<Dim<3>, Vector> position = state.vectorFields(HydroFieldNames::position);
  const FieldList<Dim<3>, Vector> velocity = state.vectorFields(HydroFieldNames::velocity);
  const FieldList<Dim<3>, SymTensor> H = state.symTensorFields(HydroFieldNames::H);
  const FieldList<Dim<3>, Scalar> rho = state.scalarFields(HydroFieldNames::massDensity);
  const FieldList<Dim<3>, Scalar> soundSpeed = state.scalarFields(HydroFieldNames::soundSpeed);
  const FieldList<Dim<3>, Scalar> u = state.scalarFields(HydroFieldNames::specificThermalEnergy);
  const FieldList<Dim<3>, Scalar> e = state.scalarFields(MHDFieldNames::totalSpecificEnergy);
  const FieldList<Dim<3>, Vector> B = state.vectorFields(MHDFieldNames::magneticInduction);

  // Derivative FieldLists.
  FieldList<Dim<3>, Scalar> maxQPressure = derivatives.scalarFields(HydroFieldNames::maxViscousPressure);
  FieldList<Dim<3>, Vector> DvDt = derivatives.vectorFields(IncrementState<Dim<3>, Field<Dim<3>, Vector> >::prefix() + HydroFieldNames::velocity);
  FieldList<Dim<3>, Scalar> DuDt, DeDt;
  if (derivatives.scalarFieldNameRegistered(IncrementState<Dim<3>, Field<Dim<3>, Scalar> >::prefix() + HydroFieldNames::specificThermalEnergy))
     DuDt = derivatives.scalarFields(IncrementState<Dim<3>, Field<Dim<3>, Scalar> >::prefix() + HydroFieldNames::specificThermalEnergy);
  if (derivatives.scalarFieldNameRegistered(IncrementState<Dim<3>, Field<Dim<3>, Scalar> >::prefix() + MHDFieldNames::totalSpecificEnergy))
     DeDt = derivatives.scalarFields(IncrementState<Dim<3>, Field<Dim<3>, Scalar> >::prefix() + MHDFieldNames::totalSpecificEnergy);
  FieldList<Dim<3>, Vector> DBDt = derivatives.vectorFields(IncrementState<Dim<3>, Field<Dim<3>, Vector> >::prefix() + MHDFieldNames::magneticInduction);

  // Prepare empty FieldLists to hold the pair-wise accelerations and 
  // heat flows if we're using the compatible energy discretization.
  FieldList<Dim<3>, vector<Vector> > QpairAccelerations = dataBase.newFluidFieldList(vector<Vector>(), "Q accelerations");
  FieldList<Dim<3>, vector<Scalar> > QpairHeatFlows = dataBase.newFluidFieldList(vector<Scalar>(), "Q heat flows");

  // Some useful pure numbers.
  const Scalar invSqrt2 = 1.0/sqrt(2.0);
  const Scalar oneThird = 1.0/3.0;

  // Reset the minimum timestep.
  mMinDt = FLT_MAX;

  // Iterate over the FluidNodeLists.
  int nodeListi = 0;
  for (DataBase<Dim<3> >::ConstFluidNodeListIterator itr = dataBase.fluidNodeListBegin();
       itr != dataBase.fluidNodeListEnd();
       ++itr, ++nodeListi) {
    const bool compatibleEnergyEvolution = (*itr)->compatibleEnergyEvolution();
    const TableKernel<Dim<3> >& W = (*itr)->Qkernel();

    // Grab the state for the NodeList we're working on.
    const Field<Dim<3>, Scalar>& massThis = *(mass[nodeListi]);
    const Field<Dim<3>, Vector>& positionThis = *(position[nodeListi]);
    const Field<Dim<3>, Vector>& velocityThis = *(velocity[nodeListi]);
    const Field<Dim<3>, SymTensor>& HThis = *(H[nodeListi]);
    const Field<Dim<3>, Scalar>& rhoThis = *(rho[nodeListi]);
    const Field<Dim<3>, Scalar>& soundSpeedThis = *(soundSpeed[nodeListi]);
    const Field<Dim<3>, Scalar>& uThis = *(u[nodeListi]);
    const Field<Dim<3>, Vector>& BThis = *(B[nodeListi]);

    // Now grab the derivative fields for this NodeList.
    Field<Dim<3>, Scalar>& maxQPressureThis = *maxQPressure[nodeListi];
    Field<Dim<3>, Vector>& DvDtThis = *DvDt[nodeListi];
    Field<Dim<3>, Scalar>* DuDtThis = (DuDt.numFields() > 0) ? DuDt[nodeListi] : 0;
    Field<Dim<3>, Scalar>* DeDtThis = (DeDt.numFields() > 0) ? DeDt[nodeListi] : 0;
    Field<Dim<3>, Vector>& DBDtThis = *DBDt[nodeListi];
    Field<Dim<3>, vector<Vector> >& QpairAccelerationsThis = *QpairAccelerations[nodeListi];
    Field<Dim<3>, vector<Scalar> >& QpairHeatFlowsThis = *QpairHeatFlows[nodeListi];

    // Iterate over the internal nodes.
    for (int i = 0; i != (*itr)->numInternalNodes(); ++i) {

      // State of node I.
      const Scalar& mi = massThis(i);
      const Vector& ri = positionThis(i);
      const Vector& vi = velocityThis(i);
      const Scalar vi2 = vi.magnitude2();
      const SymTensor& Hi = HThis(i);
      const Scalar& rhoi = rhoThis(i);
      const Scalar& ci = soundSpeedThis(i);
      const Scalar ci2 = ci*ci;
      const Scalar Hdeti = Hi.Determinant();
      const Scalar hi = 1.0/pow(Hdeti, oneThird);
      const Scalar& ui = uThis(i);
      const Vector& Bi = BThis(i);
      const Scalar Bi2 = Bi.magnitude2();
      const Scalar invRhoiMu0 = 1.0 / (rhoi*mMu0);
      const Scalar Vfi2 = ci2 + Bi2*invRhoiMu0;

      Scalar& maxQPressurei = maxQPressureThis(i);
      Vector& DvDti = DvDtThis(i);
      Vector& DBDti = DBDtThis(i);
      vector<Vector>& QpairAccelerationsi = QpairAccelerationsThis(i);
      //vector<Scalar>& QpairHeatFlowsi = QpairHeatFlowsThis(i);

      // Maximum signal velocity involving this node.
      Scalar maxVsig = -FLT_MAX;

      // Connectivity for this node.
      const vector< vector<int> >& fullConnectivity = connectivityMap.connectivityForNode(*itr, i);
      CHECK(fullConnectivity.size() == numNodeLists);

      // Iterate over the neighboring NodeLists.
      for (int nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {

        // Connectivity of this node with this NodeList.  We only need to proceed if
        // there are some nodes in this list.
        const vector<int>& connectivity = fullConnectivity[nodeListj];
        if (connectivity.size() > 0) {

          // Grab the state for the neighbor NodeList.
          const Field<Dim<3>, Scalar>& massThem = *(mass[nodeListj]);
          const Field<Dim<3>, Vector>& positionThem = *(position[nodeListj]);
          const Field<Dim<3>, Vector>& velocityThem = *(velocity[nodeListj]);
          const Field<Dim<3>, SymTensor>& HThem = *(H[nodeListj]);
          const Field<Dim<3>, Scalar>& rhoThem = *(rho[nodeListj]);
          const Field<Dim<3>, Scalar>& soundSpeedThem = *(soundSpeed[nodeListj]);
          const Field<Dim<3>, Scalar>& uThem = *(u[nodeListj]);
          const Field<Dim<3>, Vector>& BThem = *(B[nodeListj]);

          Field<Dim<3>, Scalar>& maxQPressureThem = *maxQPressure[nodeListj];
          Field<Dim<3>, Vector>& DvDtThem = *DvDt[nodeListj];
          Field<Dim<3>, Scalar>* DuDtThem = (DuDt.numFields() > 0) ? DuDt[nodeListj] : 0;
          Field<Dim<3>, Scalar>* DeDtThem = (DeDt.numFields() > 0) ? DeDt[nodeListj] : 0;
          Field<Dim<3>, Vector>& DBDtThem = *DBDt[nodeListj];
          Field<Dim<3>, vector<Vector> >& QpairAccelerationsThem = *QpairAccelerations[nodeListj];
          //Field<Dim<3>, vector<Scalar> >& QpairHeatFlowsThem = *QpairHeatFlows[nodeListj];

          const int firstGhostNodej = nodeLists[nodeListj]->firstGhostNode();

          // Iterate over the neighbors in this NodeList.
          for (vector<int>::const_iterator jItr = connectivity.begin();
               jItr != connectivity.end();
               ++jItr) {
            const int j = *jItr;
            CHECK(j < nodeLists[nodeListj]->numNodes());

            // Only proceed if this node pair has not been calculated yet.
            if ((nodeListj > nodeListi) ||
                ((nodeListj == nodeListi) && (j > i)) ||
                ((nodeListj < nodeListi) && (j >= firstGhostNodej))) {

              // State for node J.
              const Scalar& mj = massThem(j);
              const Vector& rj = positionThem(j);
              const Vector& vj = velocityThem(j);
              const Scalar vj2 = vj.magnitude2();
              const SymTensor& Hj = HThem(j);
              const Scalar& rhoj = rhoThem(j);
              const Scalar& cj = soundSpeedThem(j);
              const Scalar& cj2 = cj*cj;
              const Scalar Hdetj = Hj.Determinant();
              const Scalar hj = 1.0/pow(Hdetj, oneThird);
              const Scalar& uj = uThem(j);
              const Vector& Bj = BThem(j);
              const Scalar Bj2 = Bj.magnitude2();
              const Scalar invRhojMu0 = 1.0 / (rhoj*mMu0);
              const Scalar Vfj2 = cj2 + Bj2*invRhojMu0;

              Scalar& maxQPressurej = maxQPressureThem(j);
              Vector& DvDtj = DvDtThem(j);
              Vector& DBDtj = DBDtThem(j);
              vector<Vector>& QpairAccelerationsj = QpairAccelerationsThem(j);
              //vector<Scalar>& QpairHeatFlowsj = QpairHeatFlowsThem(j);

              // Compute the displacement.
              const Vector rij = ri - rj;
              const Vector vij = vi - vj;
              const Vector etai = Hi*rij;
              const Vector etaj = Hj*rij;
              const Scalar etaMagi = etai.magnitude();
              const Scalar etaMagj = etaj.magnitude();

              // Symmetrized gradient.
              const Vector Hetai = Hi*etai.unitVector();
              const Vector gradWi = Hetai*W.grad(etaMagi, Hdeti);

              const Vector Hetaj = Hj*etaj.unitVector();
              const Vector gradWj = Hetaj*W.grad(etaMagj, Hdetj);
              const Vector gradWij = 0.5 * (gradWi + gradWj);

              // Compute the relative signal velocity for nodes i and j 
              // approaching one another.
              const Vector rhat = rij.unitVector();
              const Scalar rhatogradWij = rhat.dot(gradWij);
              const Scalar vijorhat = vij.dot(rhat);
              const Scalar Biorhat = Bi.dot(rhat);
              const Scalar Bjorhat = Bj.dot(rhat);
              const Scalar vsi = invSqrt2 * sqrt(Vfi2 + sqrt(Vfi2*Vfi2 - 
                                    4.0*ci2*Biorhat*Biorhat*invRhoiMu0));
              const Scalar vsj = invSqrt2 * sqrt(Vfj2 + sqrt(Vfj2*Vfj2 - 
                                    4.0*cj2*Bjorhat*Bjorhat*invRhojMu0));
              const Scalar vsig = max(0.0, vsi + vsj - mBeta*vijorhat);
#if 0
if ((i == 1753) && (abs(gradWij.x()) == 0.) && (abs(gradWij.y()) > 1e-8) && (gradWij.z() > 0.))
{
   gradWijy += gradWij.y();
   cout << "rij = " << rij << ", gradWij = " << gradWij << "(" << gradWijy << ")\n";
   cout << "vsi = " << vsi << ", vsj = " << vsj << ", vsig = " << vsig << endl;
}
else if ((j == 1753) && (abs(gradWij.x()) == 0.) && (abs(gradWij.y()) > 1e-8) && (gradWij.z() > 0.))
{
   gradWijy -= gradWij.y();
   cout << "rij = " << -rij << ", gradWij = " << -gradWij << "(" << gradWijy << ")\n";
   cout << "vsi = " << vsj << ", vsj = " << vsi << ", vsig = " << vsig << endl;
}
#endif

              // Revise the minimum timestep if necessary.
              const Scalar cflvsig = max(0.0, vsi + vsj + mBeta*abs(vijorhat));
              const Scalar dtQ = 0.25 * (hi + hj) / (cflvsig + FLT_MIN);
              if (dtQ < mMinDt)
                 mMinDt = dtQ;

              // Compute the artificial dissipation terms.
              const Scalar rhoij = 0.5 * (rhoi + rhoj);
              const Scalar rhoij2 = rhoij*rhoij;
              const Vector Bij = Bi - Bj;
              Vector Fij;
              Scalar eiStar, ejStar, pressure = 0.0;
              if (vijorhat < 0.0)
              {
                // Viscous accelerations.
                Fij = 0.5 * mi * mj * mAlpha * vsig * vijorhat * gradWij / rhoij; // (Eq. 4.80, Price thesis)
//                Fij = 0.5 * mi * mj * mAlpha * vsig * vij * rij.unitVector().dot(gradWij) / rhoij; // (Eq. 4.92, Price thesis)
                DvDti += Fij / mi;
                CHECK(DvDti.magnitude() == DvDti.magnitude()); // NaN?
                DvDtj -= Fij / mj;
                CHECK(DvDtj.magnitude() == DvDtj.magnitude()); // NaN?

                // Viscous heat transfers.
                if (DuDtThis != 0)
                {
                   CHECK(DuDtThem != 0);
                   (*DuDtThis)(i) -= mj * (vsig/rhoij) * rhatogradWij * 
                                     (0.5 * mAlpha * vij.magnitude2());
                   (*DuDtThem)(j) -= mi * (vsig/rhoij) * rhatogradWij * 
                                     (0.5 * mAlpha * vij.magnitude2());
                }
                else if (DeDtThis != 0)
                {
                   CHECK(DeDtThem != 0);
                   (*DeDtThis)(i) -= mj * (vsig/rhoij) * rhatogradWij * 
                                     (0.5 * mAlpha * vij.magnitude2());
                   (*DeDtThem)(j) -= mi * (vsig/rhoij) * rhatogradWij * 
                                     (0.5 * mAlpha * vij.magnitude2());
                }

                // Dissipative terms for the total specific energy.
                eiStar = 0.5 * mAlpha * vi2 + mAlphaU * ui;
                ejStar = 0.5 * mAlpha * vj2 + mAlphaU * uj;

                // Viscous pressure.
                pressure = -rhoi * mAlpha * vsig * vijorhat / rhoij;
              }
              else
              {
                // The nodes aren't approaching each other, so there 
                // are no viscous forces.  There is, however, conductive 
                // heat transfer.
                eiStar = mAlphaU * ui;
                ejStar = mAlphaU * uj;
              }

              // Resistive inductions.
              DBDti += 0.5 * rhoi * mj * mAlphaB * (vsig/rhoij2) * 
                       Bij * rhatogradWij;
              DBDtj -= 0.5 * rhoj * mi * mAlphaB * (vsig/rhoij2) * 
                       Bij * rhatogradWij;

              // Resistive heat transfers.
              if (DuDtThis != 0)
              {
                 CHECK(DuDtThem != 0);
                 (*DuDtThis)(i) -= 0.25 * mj * (vsig/rhoij2) * rhatogradWij * 
                       (mAlphaB/mMu0) * Bij.magnitude2();
                 (*DuDtThem)(j) -= 0.25 * mi * (vsig/rhoij2) * rhatogradWij * 
                       (mAlphaB/mMu0) * Bij.magnitude2();
              }
              else if (DeDtThis != 0)
              {
                 CHECK(DeDtThem != 0);
                 (*DeDtThis)(i) -= 0.25 * mj * (vsig/rhoij2) * rhatogradWij * 
                       (mAlphaB/mMu0) * Bij.magnitude2();
                 (*DeDtThem)(j) -= 0.25 * mi * (vsig/rhoij2) * rhatogradWij * 
                       (mAlphaB/mMu0) * Bij.magnitude2();
              }

              // Magnetic terms in total specific energy dissipation.
              eiStar += 0.5*mAlphaB*Bi2/(rhoij*mMu0);
              ejStar += 0.5*mAlphaB*Bj2/(rhoij*mMu0);

              //DeDti -= 0.5 * mj * (vsig*(eiStar - ejStar)/rhoij) * 
              //         rhatogradWij;
              //DeDtj += 0.5 * mi * (vsig*(eiStar - ejStar)/rhoij) * 
              //         rhatogradWij;

              // Finally, apply artificial conductive heating.
              if (DuDtThis != 0)
              {
                 (*DuDtThis)(i) += 0.5 * mj * (vsig/rhoij) * rhatogradWij * 
                                   mAlphaU * (ui - uj);
                 (*DuDtThem)(j) += 0.5 * mi * (vsig/rhoij) * rhatogradWij * 
                                   mAlphaU * (uj - ui);
              }
              else if (DeDtThis != 0)
              {
                 (*DeDtThis)(i) += 0.5 * mj * (vsig/rhoij) * rhatogradWij * 
                                   mAlphaU * (ui - uj);
                 (*DeDtThem)(j) += 0.5 * mi * (vsig/rhoij) * rhatogradWij * 
                                   mAlphaU * (uj - ui);
              }

              // Now place the pairwise accelerations into the array for 
              // the compatible energy update (if we're using it).
              if (compatibleEnergyEvolution) {
                QpairAccelerationsi.push_back(Fij/mi);
                QpairAccelerationsj.push_back(-Fij/mj);

                // Clear out the thermal energy derivatives.
                if (DuDtThis != 0)
                {
                   (*DuDtThis)(i) = 0.0;
                   (*DuDtThem)(j) = 0.0;
                }
              }

              // Max viscous pressure.
              maxQPressurei = max(maxQPressurei, pressure);
              maxQPressurej = max(maxQPressurej, pressure);
            }
          }
        }
      }
    }
  }

  // If we're doing compatible energy evolution, we have to add the 
  // pair-wise accelerations and heat flows to the full pair-wise values.
  FieldList<Dim<3>, vector<Vector> > pairAccelerations = derivatives.vectorVectorFields(HydroFieldNames::pairAccelerations);
  nodeListi = 0;
  for (DataBase<Dim<3> >::ConstFluidNodeListIterator itr = dataBase.fluidNodeListBegin();
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

//------------------------------------------------------------------------------
// Overridden timestep method.
//------------------------------------------------------------------------------
ArtificialViscosity<Dim<3> >::TimeStepType
PriceMonaghanDissipation::
dt(const DataBase<Dim<3> >& dataBase, 
   const State<Dim<3> >& state,
   const StateDerivatives<Dim<3> >& derivs,
   const Scalar currentTime) const
{
   return TimeStepType(mMinDt, "Artificial dissipation signal speed CFL");
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Test if the ArtificialViscosity is valid, i.e., ready to use.
//------------------------------------------------------------------------------
bool
PriceMonaghanDissipation::valid() const {
  return (ArtificialViscosity<Dim<3> >::valid() &&
          this->epsilon2() > 0.0);
}

}
