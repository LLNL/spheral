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

set(NodeList_sources )

instantiate(NodeList_inst NodeList_sources)

set(NodeList_headers
    ASPHSmoothingScale.hh
    FluidNodeList.hh
    FluidNodeListInline.hh
    FluidNodeTraits.hh
    NodeList.hh
    NodeListInline.hh
    NodeListRegistrar.hh
    NodeListRegistrarInline.hh
    SPHSmoothingScale.hh
    SmoothingScaleBase.hh
    SolidNodeList.hh
    SolidNodeListInline.hh
    generateVoidNodes.hh
    secondMomentUtilities.hh
    )

set(NodeList_python_targets
    FluidNodeLists.py
    SolidNodeLists.py
    VoidNodeLists.py
    )

spheral_add_cxx_library(NodeList)

if(NOT ENABLE_CXXONLY)
  install(FILES       ${NodeList_python_targets}
          DESTINATION .
          )
endif()

list(APPEND SPHERAL_PYTHON_DEPENDS Spheral_NodeList)
set(SPHERAL_PYTHON_DEPENDS "${SPHERAL_PYTHON_DEPENDS}" PARENT_SCOPE)
