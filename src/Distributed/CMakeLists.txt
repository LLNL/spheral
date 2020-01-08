include_directories(.)
set(Distributed_sources
    Communicator.cc
    SortAndDivideRedistributeNodes1d.cc
    SortAndDivideRedistributeNodes2d.cc
    SortAndDivideRedistributeNodes3d.cc
   )

if (ENABLE_MPI)
   list(APPEND Distributed_sources
        waitAllWithDeadlockDetection.cc
       )

endif()
   set(Distributed_inst
     DistributedBoundary
     NestedGridDistributedBoundary
     TreeDistributedBoundary
     BoundingVolumeDistributedBoundary
     RedistributeNodes
     DistributeByXPosition
     SortAndDivideRedistributeNodes
     SpaceFillingCurveRedistributeNodes
     MortonOrderRedistributeNodes
     PeanoHilbertOrderRedistributeNodes
     VoronoiRedistributeNodes
    )

instantiate(Distributed_inst Distributed_sources)

set(Distributed_headers
    Communicator.hh
    CompareDomainNodesByPosition.hh
    DistributedBoundary.hh
    DistributedBoundaryInline.hh
    #DomainNode.hh
    #DomainNodeInline.hh
    NestedGridDistributedBoundary.hh
    NestedGridUtilities.hh
    RedistributeNodes.hh
    RedistributeNodesInline.hh
    TreeDistributedBoundary.hh
    waitAllWithDeadlockDetection.hh
    SortAndDivideRedistributeNodes.hh
    SortAndDivideRedistributeNodes1d.hh
    )

spheral_add_cxx_library(Distributed)

install(TARGETS              Spheral_Distributed
        EXPORT               ${PROJECT_NAME}-targets
        DESTINATION          lib
        INCLUDES DESTINATION include
        )

install(FILES       ${Distributed_headers}
        DESTINATION include/Distributed
        )

list(APPEND SPHERAL_PYTHON_DEPENDS Spheral_Distributed)
set(SPHERAL_PYTHON_DEPENDS "${SPHERAL_PYTHON_DEPENDS}" PARENT_SCOPE)
