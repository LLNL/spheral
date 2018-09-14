set(Distributed_sources
    Communicator.cc
   )

if (ENABLE_MPI)
   list(APPEND Distributed_sources
        waitAllWithDeadlockDetection.cc
       )

   if(ENABLE_2D)
     list(APPEND Distributed_sources
          SortAndDivideRedistributeNodes2d.cc
         )
   endif()
   if(ENABLE_3D)
     list(APPEND Distributed_sources
          SortAndDivideRedistributeNodes3d.cc
         )
   endif()

   set(Distributed_inst
       BoundingVolumeDistributedBoundary
       DistributeByXPosition
       DistributedBoundary
       MortonOrderRedistributeNodes
       NestedGridDistributedBoundary
       NestedGridRedistributeNodes
       PeanoHilbertOrderRedistributeNodes
       RedistributeNodes
       SortAndDivideRedistributeNodes
       SpaceFillingCurveRedistributeNodes
       TreeDistributedBoundary
       VoronoiRedistributeNodes
    )
endif()

instantiate(Distributed_inst Distributed_sources)

set(Distributed_headers
    BoundingVolumeDistributedBoundary.hh
    Communicator.hh
    CompareDomainNodesByPosition.hh
    DistributeByXPosition.hh
    DistributedBoundary.hh
    DistributedBoundaryInline.hh
    DomainNode.hh
    DomainNodeInline.hh
    MortonOrderRedistributeNodes.hh
    NestedGridDistributedBoundary.hh
    NestedGridRedistributeNodes.hh
    NestedGridRedistributeNodesInline.hh
    NestedGridUtilities.hh
    ParmetisRedistributeNodes.hh
    ParmetisRedistributeNodesInline.hh
    PeanoHilbertOrderRedistributeNodes.hh
    RedistributeNodes.hh
    RedistributeNodesInline.hh
    SortAndDivideRedistributeNodes.hh
    SortAndDivideRedistributeNodes1d.hh
    SortAndDivideRedistributeNodes2d.hh
    SortAndDivideRedistributeNodes3d.hh
    SortAndDivideRedistributeNodesInline.hh
    SpaceFillingCurveRedistributeNodes.hh
    TreeDistributedBoundary.hh
    VoronoiRedistributeNodes.hh
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
