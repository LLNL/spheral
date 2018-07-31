set(Integrator_sources2d
    CheapSynchronousRK2Inst.cc.2.cc
    IntegratorInst.cc.2.cc
    PredictorCorrectorInst.cc.2.cc
    SynchronousRK1Inst.cc.2.cc
    SynchronousRK2Inst.cc.2.cc
    SynchronousRK4Inst.cc.2.cc
    VerletInst.cc.2.cc
    )

set(Integrator_sources3d
    CheapSynchronousRK2Inst.cc.3.cc
    IntegratorInst.cc.3.cc
    PredictorCorrectorInst.cc.3.cc
    SynchronousRK1Inst.cc.3.cc
    SynchronousRK2Inst.cc.3.cc
    SynchronousRK4Inst.cc.3.cc
    VerletInst.cc.3.cc
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

INCLUDE_DIRECTORIES(../)

set(spheral_depends BOOST Eigen qhull silo polytope mpi)

blt_add_library(NAME        Spheral_Integrator
                SOURCES     ${Integrator_sources2d} ${Integrator_sources3d}
                HEADERS     ${Integrator_headers}
                DEPENDS_ON  ${spheral_depends}
                )

install(TARGETS              Spheral_Integrator
        EXPORT               ${PROJECT_NAME}-targets
        DESTINATION          lib
        INCLUDES DESTINATION include
        )

install(FILES       ${Integrator_headers}
        DESTINATION include/Integrator
        )
