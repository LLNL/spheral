set(NodeList_sources2d
    ASPHSmoothingScaleInst.cc.2.cc
    FixedSmoothingScaleInst.cc.2.cc
    FluidNodeListInst.cc.2.cc
    NodeListInst.cc.2.cc
    NodeListRegistrarInst.cc.2.cc
    SPHSmoothingScaleInst.cc.2.cc
    SmoothingScaleBaseInst.cc.2.cc
    SolidNodeListInst.cc.2.cc
    generateVoidNodesInst.cc.2.cc
    nthNodalMomentInst.cc.2.cc
    )

set(NodeList_sources3d
    ASPHSmoothingScaleInst.cc.3.cc
    FixedSmoothingScaleInst.cc.3.cc
    FluidNodeListInst.cc.3.cc
    NodeListInst.cc.3.cc
    NodeListRegistrarInst.cc.3.cc
    SPHSmoothingScaleInst.cc.3.cc
    SmoothingScaleBaseInst.cc.3.cc
    SolidNodeListInst.cc.3.cc
    generateVoidNodesInst.cc.3.cc
    nthNodalMomentInst.cc.3.cc
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
    SolidNodeList.hh
    SolidNodeListInline.hh
    generateVoidNodes.hh
    nthNodalMoment.hh
    secondMomentUtilities.hh
    )

INCLUDE_DIRECTORIES(../)

set(spheral_depends BOOST Eigen qhull silo polytope mpi)

blt_add_library(NAME        Spheral_NodeList
                SOURCES     ${NodeList_sources2d} ${NodeList_sources3d}
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
