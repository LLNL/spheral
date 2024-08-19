include_directories(.)

set(NodeList_sources )

# Add explicit instantiations
list(APPEND NodeList_sources
    ASPHSmoothingScaleInst.cc
    FluidNodeListInst.cc
    NodeListInst.cc
    NodeListRegistrarInst.cc
    SPHSmoothingScaleInst.cc
    SmoothingScaleBaseInst.cc
    SolidNodeListInst.cc
    generateVoidNodesInst.cc
    FixedSmoothingScaleInst.cc
    nthNodalMomentInst.cc
    DEMNodeListInst.cc
   )

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
    DEMNodeList.hh
    DEMNodeListInline.hh
    )

spheral_install_python_files(
  FluidNodeLists.py
  SolidNodeLists.py
  VoidNodeLists.py
  DEMNodeLists.py
  )

spheral_add_obj_library(NodeList SPHERAL_OBJ_LIBS)

