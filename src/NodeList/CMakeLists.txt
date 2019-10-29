set(NodeList_inst
  ASPHSmoothingScale
  FixedSmoothingScale
  #FluidNodeList
  generateVoidNodes
  NodeList
  NodeListRegistrar
  nthNodalMoment
  SmoothingScaleBase
  #SolidNodeList
  SPHSmoothingScale
  )

set(NodeList_sources )
instantiate(NodeList_inst NodeList_sources)

set(NodeList_headers
  ASPHSmoothingScale.hh
  FixedSmoothingScale.hh
  FluidNodeList.hh
  FluidNodeListInline.hh
  FluidNodeTraits.hh
  generateVoidNodes.hh
  NodeList.hh
  NodeListInline.hh
  NodeListRegistrar.hh
  NodeListRegistrarInline.hh
  nthNodalMoment.hh
  secondMomentUtilities.hh
  SmoothingScaleBase.hh
  SolidNodeList.hh
  SolidNodeListInline.hh
  SPHSmoothingScale.hh
  )

set(NodeList_includes
  ${SPHERAL_SOURCE_DIR}
  ${SPHERAL_SOURCE_DIR}/NodeList
  ${EIGEN_LIB_LOCATION}
  )


blt_add_library(NAME        Spheral_NodeList
                HEADERS     ${NodeList_headers}
                SOURCES     ${NodeList_sources}
                DEPENDS_ON  ${spheral_depends}
                INCLUDES    ${NodeList_includes}
                )

install(TARGETS             Spheral_NodeList
        EXPORT              ${PROJECT_NAME}-targets
        DESTINATION         lib
        INCLUDES DESTINATION include
        )

install(FILES       ${NodeList_headers}
        DESTINATION include/NodeList
        )
