include_directories(.)

set(Integrator_sources )

# Add explicit instantiations
list(APPEND Integrator_sources
     IntegratorInst.cc
     PredictorCorrectorInst.cc
     SynchronousRK1Inst.cc
     SynchronousRK2Inst.cc
     SynchronousRK4Inst.cc
     CheapSynchronousRK2Inst.cc
     VerletInst.cc
     )

set(Integrator_headers
    CheapSynchronousRK2.hh
    Integrator.hh
    IntegratorInline.hh
    PredictorCorrector.hh
    SynchronousRK1.hh
    SynchronousRK2.hh
    SynchronousRK4.hh
    Verlet.hh
    )

spheral_add_obj_library(Integrator SPHERAL_OBJ_LIBS)

