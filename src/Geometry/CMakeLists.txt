set(Geometry_inst
  invertRankNTensor
  )

set(Geometry_sources
    buildEigenVector.cc
    #clipConvexFacetedVolumeByPlanes2d.cc
    #clipConvexFacetedVolumeByPlanes3d.cc
    #clipFacetedVolumeByPlanes1d.cc
    #clipFacetedVolumeByPlanes3d.cc
    computeEigenValues.cc
    Dimension.cc
    EigenStruct.cc
    fakemain.cc
    Geom3Vector.cc
    GeomFacet2d.cc
    GeomFacet3d.cc
    GeomFifthRankTensor.cc
    GeomFourthRankTensor.cc
    GeomPlane.cc
    GeomPolygon.cc
    GeomPolyhedron.cc
    GeomSymmetricTensor.cc
    GeomSymmetricTensor_default.cc
    GeomSymmetricTensor_eigen.cc
    GeomTensor.cc
    GeomTensor_default.cc
    GeomTensor_eigen.cc
    GeomThirdRankTensor.cc
    GeomVector.cc
    invertRankNTensor.cc
    Jacobi.cc
    Jacobi2.cc
    polyclipper2d.cc
    polyclipper3d.cc
    #TestGeometry.cc
    #TestMTL.cc
    ../Utilities/SpheralTimers.cc
    )

instantiate(Geometry_inst Geometry_sources)

set(Geometry_headers
    aggregateFacetedVolumes.hh
    Box1d.hh
    Box1dInline.hh
    buildEigenVector.hh
    CellFaceFlag.hh
    clipConvexFacetedVolumeByPlanes.hh
    clipFacetedVolumeByPlanes.hh
    computeEigenValues.hh
    Dimension.hh
    EigenStruct.hh
    EigenStruct_fwd.hh
    EigenStructInline.hh
    FacetedVolumeUtilities.hh
    findEigenValues3.hh
    Geom3Vector.hh
    Geom3VectorInline.hh
    GeomFacet2d.hh
    GeomFacet2dInline.hh
    GeomFacet3d.hh
    GeomFacet3dInline.hh
    GeomFifthRankTensor.hh
    GeomFifthRankTensor_fwd.hh
    GeomFifthRankTensorInline.hh
    GeomFourthRankTensor.hh
    GeomFourthRankTensor_fwd.hh
    GeomFourthRankTensorInline.hh
    GeomPlane.hh
    GeomPlaneInline.hh
    GeomPolygon.hh
    GeomPolygonInline.hh
    GeomPolyhedron.hh
    GeomPolyhedronInline.hh
    GeomSymmetricTensor.hh
    GeomSymmetricTensor_default.hh
    GeomSymmetricTensor_eigen.hh
    GeomSymmetricTensor_fwd.hh
    GeomSymmetricTensorBase.hh
    GeomSymmetricTensorInline_default.hh
    GeomSymmetricTensorInline_eigen.hh
    GeomTensor.hh
    GeomTensor_default.hh
    GeomTensor_eigen.hh
    GeomTensor_fwd.hh
    GeomTensorBase.hh
    GeomTensorInline_default.hh
    GeomTensorInline_eigen.hh
    GeomThirdRankTensor.hh
    GeomThirdRankTensor_fwd.hh
    GeomThirdRankTensorInline.hh
    GeomVector.hh
    GeomVector_array.hh
    GeomVector_default.hh
    GeomVector_eigen.hh
    GeomVector_fwd.hh
    GeomVectorBase_default.hh
    GeomVectorInline_array.hh
    GeomVectorInline_default.hh
    GeomVectorInline_eigen.hh
    innerDoubleProduct.hh
    innerProduct.hh
    invertRankNTensor.hh
    Jacobi.hh
    Jacobi2.hh
    makefile.in
    MathTraits.hh
    outerProduct.hh
    polyclipper.hh
    polyclipper_utilities.hh
    RankNTensor.hh
    RankNTensorInline.hh
    TensorType.hh
    )

set(Geometry_includes
  ${SPHERAL_SOURCE_DIR}
  ${SPHERAL_SOURCE_DIR}/Geometry
  ${SPHERAL_SOURCE_DIR}/Utilities
  ${EIGEN_LIB_LOCATION}
  ${QHULL_LIB_LOCATION}
  )

blt_add_library(NAME        Spheral_Geometry
                HEADERS     ${Geometry_headers}
                SOURCES     ${Geometry_sources}
                DEPENDS_ON  ${spheral_depends}
                INCLUDES    ${Geometry_includes}
                )

install(TARGETS             Spheral_Geometry
        EXPORT              ${PROJECT_NAME}-targets
        DESTINATION         lib
        INCLUDES DESTINATION include
        )

install(FILES       ${Geometry_headers}
        DESTINATION include/Geometry
        )


list(APPEND SPHERAL_PYTHON_DEPENDS Spheral_Geometry)
set(SPHERAL_PYTHON_DEPENDS "${SPHERAL_PYTHON_DEPENDS}" PARENT_SCOPE)
