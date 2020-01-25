set(NodeList_inst
    ASPHSmoothingScale
    FluidNodeList
    NodeList
    NodeListRegistrar
    SPHSmoothingScale
    SmoothingScaleBase
    SolidNodeList
    generateVoidNodes
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

blt_add_library(NAME        Spheral_NodeList
                SOURCES     ${NodeList_sources}
                HEADERS     ${NodeList_headers}
                DEPENDS_ON  ${spheral_depends}
                )

install(TARGETS              Spheral_NodeList
        EXPORT               ${PROJECT_NAME}-targets
        DESTINATION          lib
        INCLUDES DESTINATION include
        )

install(FILES       ${NodeList_headers}
        DESTINATION include/NodeList
        )
