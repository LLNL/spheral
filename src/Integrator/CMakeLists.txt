include_directories(.)
set(Integrator_inst
    Integrator
    ImplicitIntegrator
    PredictorCorrector 
    ForwardEuler
    SynchronousRK2
    SynchronousRK4
    CheapSynchronousRK2
    Verlet
    )

set(Integrator_sources )

instantiate(Integrator_inst Integrator_sources)

set(Integrator_headers
    Integrator.hh
    ImplicitIntegrator.hh
    CheapSynchronousRK2.hh
    PredictorCorrector.hh
    ForwardEuler.hh
    SynchronousRK2.hh
    SynchronousRK4.hh
    Verlet.hh
    )

spheral_add_obj_library(Integrator SPHERAL_OBJ_LIBS)

