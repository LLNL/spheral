include_directories(.)
set(NodeList_inst
    ASPHSmoothingScale
    FluidNodeList
    NodeList
    NodeListRegistrar
    SPHSmoothingScale
    SmoothingScaleBase
    SolidNodeList
    generateVoidNodes
    FixedSmoothingScale
    nthNodalMoment
   )

if(SPHERAL_ENABLE_DEM)
  list(APPEND NodeList_inst DEMNodeList)
endif()

set(NodeList_sources )

instantiate(NodeList_inst NodeList_sources)

set(NodeList_headers
    ASPHSmoothingScale.hh
    FixedSmoothingScale.hh
    FluidNodeList.hh
    FluidNodeListInline.hh
    FluidNodeTraits.hh
    NodeList.hh
    NodeListInline.hh
    NodeListRegistrar.hh
    NodeListRegistrarInline.hh
    SPHSmoothingScale.hh
    SmoothingScaleBase.hh
    SmoothingScaleBaseInline.hh
    SolidNodeList.hh
    SolidNodeListInline.hh
    generateVoidNodes.hh
    nthNodalMoment.hh
    secondMomentUtilities.hh
    )

if(SPHERAL_ENABLE_DEM)
  list(APPEND NodeList_headers DEMNodeList.hh DEMNodeListInline.hh)
endif()

spheral_install_python_files(
  FluidNodeLists.py
  SolidNodeLists.py
  VoidNodeLists.py
  )

if(SPHERAL_ENABLE_DEM)
  spheral_install_python_files(DEMNodeLists.py)
endif()

spheral_add_obj_library(NodeList SPHERAL_OBJ_LIBS)

