if (ENABLE_MPI)
   set(Distributed_sources2d
       BoundingVolumeDistributedBoundaryInst.cc.2.cc
       DistributeByXPositionInst.cc.2.cc
       DistributedBoundaryInst.cc.2.cc
       MortonOrderRedistributeNodesInst.cc.2.cc
       NestedGridDistributedBoundaryInst.cc.2.cc
       NestedGridRedistributeNodesInst.cc.2.cc
       PeanoHilbertOrderRedistributeNodesInst.cc.2.cc
       RedistributeNodesInst.cc.2.cc
       SortAndDivideRedistributeNodesInst.cc.2.cc
       SortAndDivideRedistributeNodes2d.cc
       SpaceFillingCurveRedistributeNodesInst.cc.2.cc
       TreeDistributedBoundaryInst.cc.2.cc
       VoronoiRedistributeNodesInst.cc.2.cc
       waitAllWithDeadlockDetection.cc
       )

   set(Distributed_sources3d
       BoundingVolumeDistributedBoundaryInst.cc.3.cc
       DistributeByXPositionInst.cc.3.cc
       DistributedBoundaryInst.cc.3.cc
       MortonOrderRedistributeNodesInst.cc.3.cc
       NestedGridDistributedBoundaryInst.cc.3.cc
       NestedGridRedistributeNodesInst.cc.3.cc
       PeanoHilbertOrderRedistributeNodesInst.cc.3.cc
       RedistributeNodesInst.cc.3.cc
       SortAndDivideRedistributeNodesInst.cc.3.cc
       SortAndDivideRedistributeNodes3d.cc
       SpaceFillingCurveRedistributeNodesInst.cc.3.cc
       TreeDistributedBoundaryInst.cc.3.cc
       VoronoiRedistributeNodesInst.cc.3.cc
       waitAllWithDeadlockDetection.cc
       )
endif()

set(Distributed_sources
    Communicator.cc
    )

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

INCLUDE_DIRECTORIES(../)

set(spheral_depends BOOST Eigen qhull silo polytope mpi)

blt_add_library(NAME        Spheral_Distributed
                SOURCES     ${Distributed_sources2d} ${Distributed_sources3d} ${Distributed_sources}
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
