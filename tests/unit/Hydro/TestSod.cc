#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "Geometry/Dimension.hh"
#include "Neighbor/Neighbor.hh"
#include "Kernel/GaussianKernel.hh"
#include "Kernel/TableKernel.hh"
#include "Material/GammaLawGas.hh"
#include "Material/MKSUnits.hh"
#include "Material/PhysicalConstants.hh"
#include "NodeList/SphNodeList.hh"
#include "Neighbor/NestedGridNeighbor.hh"
#include "Boundary/ReflectingBoundary.hh"
#include "DataBase/DataBase.hh"
#include "ArtificialViscosity/MonaghanGingoldViscosity.hh"
#include "Hydro/Hydro.hh"
#include "Integrator/SynchronousRK2.hh"


// Define a C++ stripped down version of the controller.
template<typename Dimension, typename FluidNodeListType, typename KernelType>
class SpheralController {
public:
  SpheralController(Integrator<Dimension, FluidNodeListType>& integrator, 
		    KernelType& kernel,
		    const vector<Boundary<Dimension, FluidNodeListType>*>& bcs):
    mIntegratorPtr(&integrator),
    mKernelPtr(&kernel),
    mTotalSteps(0),
    mBoundaryConditions(bcs) {
    applyBoundaryConditions();
    mIntegratorPtr->initialize();
  }

  double time() const {
    return mIntegratorPtr->currentTime();
  }

  double lastDt() const {
    return mIntegratorPtr->lastDt();
  }

  void applyBoundaryConditions() {
    DataBase<Dimension, FluidNodeListType>& db = 
      const_cast<DataBase<Dimension, FluidNodeListType>&>(mIntegratorPtr->dataBase());
    for (DataBase<Dimension, FluidNodeListType>::FluidNodeListIterator 
	   nodeListItr = db.fluidNodeListBegin();
	 nodeListItr < db.fluidNodeListEnd();
	 ++nodeListItr) {
      (*nodeListItr)->numGhostNodes(0);
      (*nodeListItr)->neighborPtr()->updateNodes();
      for (vector<Boundary<Dimension, FluidNodeListType>*>::iterator
	     bcItr = mBoundaryConditions.begin();
	   bcItr < mBoundaryConditions.end(); ++bcItr) {
	(*bcItr)->setGhostNodes(**nodeListItr);
	(*bcItr)->applyBoundary((*nodeListItr)->mass());
	(*bcItr)->applyBoundary((*nodeListItr)->Hfield());
	(*bcItr)->applyBoundary((*nodeListItr)->velocity());
	(*bcItr)->applyBoundary((*nodeListItr)->specificThermalEnergy());
      }
      (*nodeListItr)->neighborPtr()->updateNodes();
      if (mIntegratorPtr->sumForMassDensity()) 
	(*nodeListItr)->updateMassDensity(*mKernelPtr);
      for (vector<Boundary<Dimension, FluidNodeListType>*>::iterator
	     bcItr = mBoundaryConditions.begin();
	   bcItr < mBoundaryConditions.end(); ++bcItr) {
	(*bcItr)->applyBoundary((*nodeListItr)->massDensity());
      }
      (*nodeListItr)->updateWeight();
      for (vector<Boundary<Dimension, FluidNodeListType>*>::iterator
	     bcItr = mBoundaryConditions.begin();
	   bcItr < mBoundaryConditions.end(); ++bcItr) {
	(*bcItr)->applyBoundary((*nodeListItr)->weight());
      }
    }
  }

  void appendBoundary(Boundary<Dimension, FluidNodeListType>& boundary) {
    mBoundaryConditions.push_back(&boundary);
  }

//     # Smooth the physical variables.
//     def smoothState(self, smoothIters=1):
//         for iter in xrange(smoothIters):
//             db = self.integrator.dataBase
//             smoothedVelocity = db.fluidVelocity.smoothFields(db.fluidHfield, self.kernel)
//             smoothedSpecificThermalEnergy = db.fluidSpecificThermalEnergy.smoothFields(db.fluidHfield, self.kernel)
//             smoothedHfield = db.fluidHfield.smoothFields(db.fluidHfield, self.kernel)
//             fieldID = 0
//             for nodeList in db.fluidNodeLists:
//                 nodeList.velocity[:] = smoothedVelocity[fieldID][:]
//                 nodeList.specificThermalEnergy[:] = smoothedSpecificThermalEnergy[fieldID][:]
//                 nodeList.Hfield[:] = smoothedHfield[fieldID][:]
            
//             # After everything is done, we have to redo the boundary conditions.
//             self.applyBoundaryConditions(self.kernel)

  // Advance the system to the given simulation time.  The user can also
  // specify a max number of steps to take.
  void advance(double goalTime, int maxSteps=0) {
    int currentSteps = 0;
    while (this->time() < goalTime &&
	   (maxSteps == 0 || currentSteps < maxSteps)) {
      mIntegratorPtr->step(goalTime);
      currentSteps += 1;
      mTotalSteps += 1;
      cout << "Cycle: " << mTotalSteps
	   << ", Time: " << this->time()
	   << ", TimeStep: " << this->lastDt()
	   << endl;
    }
  }

private:
  Integrator<Dimension, FluidNodeListType>* mIntegratorPtr;
  KernelType* mKernelPtr;
  int mTotalSteps;
  vector<Boundary<Dimension, FluidNodeListType>*> mBoundaryConditions;
};

int main() {

  // Run parameters.
  int nx1 = 100;
  int nx2 = 25;
  int nx = nx1 + nx2;

  double rho1 = 1.0;
  double rho2 = 0.25;

  double P1 = 1.0;
  double P2 = 0.1795;

  double x0 = -0.5;
  double x1 = 0.0;
  double x2 = 0.5;

  double gamma = 1.4;
  double mu = 1.0;

  double Cl = 0.75;
  double Cq = 1.5;
  double epsilon2 = 1e-2;

  double HsmoothMin = 0.0001;
  double HsmoothMax = 0.1;
  double cfl = 0.1;

  int numGridLevels = 10;
  double topGridCellSize = 0.25;
  Dim<1>::Vector origin(0.0);

  Neighbor<Dim<1>, Dim<1>::Scalar>::NeighborSearchType neighborSearchType = Neighbor<Dim<1>, Dim<1>::Scalar>::GatherScatter;

  double goalTime = 0.15;
  int maxSteps = 500;
  double dt = 0.0001;
  int sumForMassDensity = 1;

  double m1 = (x1 - x0)*rho1/nx1;
  double m2 = (x2 - x1)*rho2/nx2;

  double eps1 = P1/((gamma - 1.0)*rho1);
  double eps2 = P2/((gamma - 1.0)*rho2);

  double dx1 = (x1 - x0)/nx1;
  double dx2 = (x2 - x1)/nx2;

  double h1 = 1.0/(2.0*dx1);
  double h2 = 1.0/(2.0*dx2);

  // Set up the Kernel.
  GaussianKernel< Dim<1> > W;
  TableKernel< Dim<1> > WT;
  WT.setTableData(W, 100);

  // Set up the material.
  GammaLawGas<Dim<1>, PhysicalConstants<MKSUnits> > eos(gamma, mu);

  // Set up the node list.
  SphNodeList< Dim<1> > nodes1(nx, eos);

  // Set node properties.
  Field<Dim<1>, Dim<1>::Scalar>& mass = nodes1.mass();
  Field<Dim<1>, Dim<1>::Vector>& positions = nodes1.positions();
  Field<Dim<1>, Dim<1>::Scalar>& eps = nodes1.specificThermalEnergy();
  Field<Dim<1>, Dim<1>::Vector>& vel = nodes1.velocity();
  Field<Dim<1>, Dim<1>::Scalar>& H = nodes1.Hfield();
  Field<Dim<1>, Dim<1>::Scalar>& rho = nodes1.massDensity();
  for (int ix = 0; ix < nx1; ++ix) {
    int nodeID = ix;
    mass(nodeID) = m1;
    positions(nodeID) = x0 + (ix + 0.5)*dx1;
    eps(nodeID) = eps1;
    vel(nodeID) = Dim<1>::Vector(0.0);
    H(nodeID) = h1;
    if (!sumForMassDensity) rho(nodeID) = rho1;
  }
  for (int ix = 0; ix < nx2; ++ix) {
    int nodeID = ix + nx1;
    mass(nodeID) = m2;
    positions(nodeID) = x1 + (ix + 0.5)*dx2;
    eps(nodeID) = eps2;
    vel(nodeID) = Dim<1>::Vector(0.0);
    H(nodeID) = h2;
    if (!sumForMassDensity) rho(nodeID) = rho2;
  }
  cout << "Created Sph node list at " << &nodes1 << endl
       << "\tnumNodes = " << nodes1.numNodes() << endl;
  
  // Construct a neighbor object.
  NestedGridNeighbor<Dim<1>, Dim<1>::Scalar> neighbor(nodes1,
						      neighborSearchType,
						      numGridLevels,
						      topGridCellSize,
						      origin,
						      WT.kernelExtent());
  nodes1.registerNeighbor(neighbor);

  // Create boundary conditions.
  GeomPlane< Dim<1> > xPlane0(Dim<1>::Vector(x0), Dim<1>::Vector(1.0));
  GeomPlane< Dim<1> > xPlane1(Dim<1>::Vector(x2), Dim<1>::Vector(-1.0));
  ReflectingBoundary<Dim<1>, SphNodeList< Dim<1> > > xbc0(xPlane0);
  ReflectingBoundary<Dim<1>, SphNodeList< Dim<1> > > xbc1(xPlane1);

  // Construct a data base.
  DataBase<Dim<1>, SphNodeList< Dim<1> > > db;
  db.appendNodeList(nodes1);
  cout << "DataBase at " << &db << endl
       << "\tnumNodeLists: " << db.numNodeLists() << endl;

  // Construct a standard Monaghan-Gingold viscosity.
  MonaghanGingoldViscosity<Dim<1>, SphNodeList< Dim<1> > > q(Cl, Cq);

  // Construct the hydro object.
  Hydro<Dim<1>, SphNodeList< Dim<1> >, TableKernel< Dim<1> >, TableKernel< Dim<1> > > hydro(q);
  hydro.setKernel(WT);
  hydro.setPiKernel(WT);
  hydro.setCfl(cfl);
  cout << "Created Sph hydro object at " << &hydro << endl
       << "\thydro.cfl: " << hydro.cfl() << endl
       << "\thydro.valid: " << hydro.valid() << endl;

  // Construct a synchronous RK2 integrator.
  SynchronousRK2<Dim<1>, SphNodeList< Dim<1> > > integrator(db);
  integrator.appendPhysicsPackage(hydro);
  integrator.appendBoundary(xbc0);
  integrator.appendBoundary(xbc1);
  integrator.setHsmoothMin(HsmoothMin);
  integrator.setHsmoothMax(HsmoothMax);
  integrator.setLastDt(dt);
  integrator.setSumForMassDensity(sumForMassDensity);

  // Build our little controller.
  vector<Boundary<Dim<1>, SphNodeList< Dim<1> > >*> bcVector;
  bcVector.push_back(&xbc0);
  bcVector.push_back(&xbc1);
  SpheralController<Dim<1>, SphNodeList< Dim<1> >, TableKernel< Dim<1> > >
    control(integrator, WT, bcVector);

  // Advance to the end time.
  control.advance(goalTime, maxSteps);

  // At the end, write out the internal field values.
  ofstream os("TestSod.out");
  Field<Dim<1>, Dim<1>::Scalar> P = nodes1.pressure();
  for (int nodeID = 0; nodeID < nodes1.numInternalNodes(); ++nodeID) {
    os << positions(nodeID).x() << "  "
       << mass(nodeID) << "  "
       << vel(nodeID).x() << "  "
       << rho(nodeID) << "  "
       << eps(nodeID) << "  "
       << P(nodeID) << "  "
       << 1.0/H(nodeID) << endl;
  }

  return 0;
}
