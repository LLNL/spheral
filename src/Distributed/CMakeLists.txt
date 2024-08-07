include_directories(.)

set(Distributed_sources
    Communicator.cc
    Process.cc
    RegisterMPIDataTypes.cc
    )

set(Distributed_inst )

set(Distributed_headers
    allReduce.hh
    Communicator.hh
    Process.hh
    RegisterMPIDataTypes.hh
    RegisterMPIDataTypesInline.hh
    )

if (ENABLE_MPI)
  #----------------------------------------------------------------------------
  # MPI parallel on
  #----------------------------------------------------------------------------
  list(APPEND Distributed_sources waitAllWithDeadlockDetection.cc)

  if (ENABLE_1D)
    list(APPEND Distributed_sources SortAndDivideRedistributeNodes1d.cc)
  endif()

  if (ENABLE_2D)
    list(APPEND Distributed_sources SortAndDivideRedistributeNodes2d.cc)
  endif()

  if (ENABLE_3D)
    list(APPEND Distributed_sources SortAndDivideRedistributeNodes3d.cc)
  endif()

  list(APPEND Distributed_inst
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

 list(APPEND Distributed_headers
      BoundingVolumeDistributedBoundary.hh
      BoundingVolumeDistributedBoundaryInline.hh
      CompareDomainNodesByPosition.hh
      DistributeByXPosition.hh
      DistributedBoundary.hh
      DistributedBoundaryInline.hh
      MortonOrderRedistributeNodes.hh
      NestedGridDistributedBoundary.hh
      NestedGridDistributedBoundaryInline.hh
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
      TreeDistributedBoundaryInline.hh
      VoronoiRedistributeNodes.hh
      waitAllWithDeadlockDetection.hh
      )

  set(MPIPY_FILE_NAME "mpi_mpi4py.py")

else()
  #----------------------------------------------------------------------------
  # MPI parallel off
  #----------------------------------------------------------------------------
  set(MPIPY_FILE_NAME "fakempi.py")

endif()

if (NOT ENABLE_CXXONLY)
  install(FILES       ${MPIPY_FILE_NAME}
          DESTINATION ${SPHERAL_SITE_PACKAGES_PATH}/Spheral
          RENAME      mpi.py
          )
endif()
instantiate(Distributed_inst Distributed_sources)

spheral_add_obj_library(Distributed SPHERAL_OBJ_LIBS)
