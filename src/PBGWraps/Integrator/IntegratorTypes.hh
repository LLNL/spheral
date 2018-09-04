#ifndef __PBGWRAPS_INTEGRATORTYPES__
#define __PBGWRAPS_INTEGRATORTYPES__

#include "Geometry/Dimension.hh"
#include "Integrator/Integrator.hh"
#include "Integrator/PredictorCorrector.hh"
#include "Integrator/SynchronousRK1.hh"
#include "Integrator/SynchronousRK2.hh"
#include "Integrator/SynchronousRK4.hh"
#include "Integrator/CheapSynchronousRK2.hh"
#include "Integrator/Verlet.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Names!
//------------------------------------------------------------------------------
typedef Integrator<Dim<1> > Integrator1d;
typedef Integrator<Dim<2> > Integrator2d;
typedef Integrator<Dim<3> > Integrator3d;

typedef PredictorCorrector<Dim<1> > PredictorCorrectorIntegrator1d;
typedef PredictorCorrector<Dim<2> > PredictorCorrectorIntegrator2d;
typedef PredictorCorrector<Dim<3> > PredictorCorrectorIntegrator3d;

typedef SynchronousRK1<Dim<1> > SynchronousRK1Integrator1d;
typedef SynchronousRK1<Dim<2> > SynchronousRK1Integrator2d;
typedef SynchronousRK1<Dim<3> > SynchronousRK1Integrator3d;

typedef SynchronousRK2<Dim<1> > SynchronousRK2Integrator1d;
typedef SynchronousRK2<Dim<2> > SynchronousRK2Integrator2d;
typedef SynchronousRK2<Dim<3> > SynchronousRK2Integrator3d;

typedef SynchronousRK4<Dim<1> > SynchronousRK4Integrator1d;
typedef SynchronousRK4<Dim<2> > SynchronousRK4Integrator2d;
typedef SynchronousRK4<Dim<3> > SynchronousRK4Integrator3d;

typedef CheapSynchronousRK2<Dim<1> > CheapSynchronousRK2Integrator1d;
typedef CheapSynchronousRK2<Dim<2> > CheapSynchronousRK2Integrator2d;
typedef CheapSynchronousRK2<Dim<3> > CheapSynchronousRK2Integrator3d;

typedef Verlet<Dim<1> > VerletIntegrator1d;
typedef Verlet<Dim<2> > VerletIntegrator2d;
typedef Verlet<Dim<3> > VerletIntegrator3d;

}

#endif
