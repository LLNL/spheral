set(Geometry_inst
    buildEigenVector
    clipConvexFacetedVolumeByPlanes2d
    clipConvexFacetedVolumeByPlanes3d
    clipFacetedVolumeByPlanes1d
    clipFacetedVolumeByPlanes3d
    computeEigenValues
    Dimension
    EigenStruct
    fakemain
    Geom3Vector
    GeomFacet2d
    GeomFacet3d
    GeomFifthRankTensor
    GeomFourthRankTensor
    GeomPlane
    GeomPolygon
    GeomPolyhedron
    GeomSymmetricTensor
    GeomSymmetricTensor_default
    GeomSymmetricTensor_eigen
    GeomTensor
    GeomTensor_default
    GeomTensor_eigen
    GeomThirdRankTensor
    GeomVector
    invertRankNTensor
    Jacobi
    Jacobi2
    polyclipper2d
    polyclipper3d
    TestGeometry
    TestMTL
    )

set(Geometry_sources )
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

blt_add_library(NAME	    Spheral_Geometry
		HEADERS     ${Geometry_headers}
		SOURCES     ${Geometry_sources}
		DEPENDS_ON  ${spheral_depends}
		)

install(TARGETS		    Spheral_Geometry
	EXPORT		    ${PROJECT_NAME}-targets
	DESTINATION	    lib
	INCLUDES DESTINATION include
	)

install(FILES	    ${Geometry_headers}
	DESTINATION include/Geometry
	)


