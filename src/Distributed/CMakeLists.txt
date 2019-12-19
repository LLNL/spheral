set(Distributed_sources
    Communicator.cc
   )

if (ENABLE_MPI)
   list(APPEND Distributed_sources
        waitAllWithDeadlockDetection.cc
       )

   set(Distributed_inst
       DistributedBoundary
       NestedGridDistributedBoundary
       RedistributeNodes
       TreeDistributedBoundary
    )
endif()

instantiate(Distributed_inst Distributed_sources)

set(Distributed_headers
    Communicator.hh
    CompareDomainNodesByPosition.hh
    DistributedBoundary.hh
    DistributedBoundaryInline.hh
    NestedGridDistributedBoundary.hh
    NestedGridUtilities.hh
    RedistributeNodes.hh
    RedistributeNodesInline.hh
    TreeDistributedBoundary.hh
    waitAllWithDeadlockDetection.hh
    )

blt_add_library(NAME        Spheral_Distributed
                SOURCES     ${Distributed_sources}
                HEADERS     ${Distributed_headers}
                DEPENDS_ON  ${spheral_depends}
                )

install(TARGETS              Spheral_Distributed
        EXPORT               ${PROJECT_NAME}-targets
        DESTINATION          lib
        INCLUDES DESTINATION include
        )

install(FILES       ${Distributed_headers}
        DESTINATION include/Distributed
        )
