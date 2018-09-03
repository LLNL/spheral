#include "MHD/CurrentDensityUpdatePolicy.hh"
#include "Hydro/HydroFieldNames.hh"
#include "MHD/ConductingFluidNodeList.hh"
#include "MHD/MHDFieldNames.hh"
#include "Geometry/Dimension.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/UpdatePolicyBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "Utilities/DBC.hh"

namespace Spheral {


//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
CurrentDensityUpdatePolicy::
CurrentDensityUpdatePolicy(const TableKernel<Dim<3> >& kernel,
                           const DataBase<Dim<3> >& dataBase,
                           double mu0):
   UpdatePolicyBase<Dim<3>, Field<Dim<3>, Dim<3>::Vector> >(HydroFieldNames::position,
                                                            MHDFieldNames::magneticInduction),
   mKernel(kernel),
   mDataBase(dataBase),
   mMu0(mu0)
{
} // end constructor
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
CurrentDensityUpdatePolicy::
~CurrentDensityUpdatePolicy()
{
} // end destructor
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
void
CurrentDensityUpdatePolicy::
update(const KeyType& key,
       State<Dim<3> >& state,
       StateDerivatives<Dim<3> >& derivs,
       const double multiplier,
       const double t,
       const double dt)
{
  typedef Dim<3>::Scalar Scalar;
  typedef Dim<3>::Vector Vector;
  typedef Dim<3>::SymTensor SymTensor;

  // Get access to data.
  const FieldList<Dim<3>, Dim<3>::Scalar>& m = state.scalarFields(HydroFieldNames::mass);
  const FieldList<Dim<3>, Dim<3>::Scalar>& rho = state.scalarFields(HydroFieldNames::massDensity);
  const FieldList<Dim<3>, Dim<3>::Vector>& x = state.vectorFields(HydroFieldNames::position);
  const FieldList<Dim<3>, Dim<3>::SymTensor>& H = state.symTensorFields(HydroFieldNames::H);
  const FieldList<Dim<3>, Dim<3>::Vector>& B = state.vectorFields(MHDFieldNames::magneticInduction);
  const FieldList<Dim<3>, Dim<3>::Scalar>& Omega = state.scalarFields(HydroFieldNames::omegaGradh);
  FieldList<Dim<3>, Dim<3>::Vector> J = state.vectorFields(MHDFieldNames::currentDensity);

  // Compute the curl of the magnetic induction to obtain the current density.
  // Now compute the gravitational acceleration, which is the negative gradient of the 
  // potential, and add it to the nodal accelerations.
  const ConnectivityMap<Dim<3> >& connectivityMap = mDataBase.connectivityMap();
  int numNodeLists = J.numFields();
  for (int iNodeList = 0; iNodeList < numNodeLists; ++iNodeList) {
    const NodeList<Dim<3> >& nodeList = J[iNodeList]->nodeList(); 
    for (int iNode = 0; iNode < nodeList.numInternalNodes(); ++iNode)
    {
      // Get data for the ith node.
      Scalar mi = (*m[iNodeList])[iNode];
      Scalar rhoi = (*m[iNodeList])[iNode];
      const Vector& Bi = (*B[iNodeList])[iNode];
      const SymTensor& Hi = (*H[iNodeList])[iNode];
      Scalar detHi = Hi.Determinant();
      Scalar Omegai = (*Omega[iNodeList])[iNode];
      Vector Ji;
    
      // Iterate over the neighboring NodeLists.
      const vector< vector<int> >& fullConnectivity = 
         connectivityMap.connectivityForNode(&nodeList, iNode);
      for (int jNodeList = 0; jNodeList != numNodeLists; ++jNodeList) {
        const vector<int>& connectivity = fullConnectivity[jNodeList];
        if (connectivity.size() > 0) {
          // Loop over the neighbors.
          for (vector<int>::const_iterator jItr = connectivity.begin();
               jItr != connectivity.end();
               ++jItr) {
            int jNode = *jItr;
            // Get data for the jth node.
            Vector xij = (*x[iNodeList])[iNode] - (*x[jNodeList])[jNode];
            CHECK(xij.magnitude2() != 0.0);
            Scalar mj = (*m[jNodeList])[jNode];
            Scalar rhoj = (*rho[jNodeList])[jNode];
            const Vector& Bj = (*B[jNodeList])[jNode];
            Vector Bij = Bi - Bj;

            const SymTensor& Hj = (*H[jNodeList])[jNode];
            Scalar detHj = Hj.Determinant();

            Vector etai = Hi.dot(xij);
            Scalar dWi = mKernel.gradValue(etai.magnitude(), detHi);
            Vector etaiHat = etai.unitVector();
            Vector gradWi = dWi * Hi.dot(etaiHat);

            Vector etaj = Hj.dot(xij);
            Scalar dWj = mKernel.gradValue(etaj.magnitude(), detHj);
            Vector etajHat = etaj.unitVector();
            Vector gradWj = dWj * Hj.dot(etajHat);

            // Compute the curl of B.
            Vector gradWij = 0.5 * (gradWi + gradWj);
            Ji += (mj/rhoj) * Bij.cross(gradWij) / mMu0;
          }
        }
      }

      // Etch the computed current density into its field.
      (*J[iNodeList])[iNode] = Ji;
    } 
  } // end for
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
bool
CurrentDensityUpdatePolicy::
operator==(const Spheral::UpdatePolicyBase<Dim<3>, Field<Dim<3>, Dim<3>::Vector> >& rhs) const
{
   // We're only equal if the other guy is the same type.
   const CurrentDensityUpdatePolicy* rhsPtr = 
      dynamic_cast<const CurrentDensityUpdatePolicy*>(&rhs);
   if (rhsPtr == 0) {
      return false;
   } else {
      return true;
   }
}
//----------------------------------------------------------------------------

}
