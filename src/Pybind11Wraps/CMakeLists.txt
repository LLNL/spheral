set(PYTHONPKGS
  "polytope Distributed OpenMP CXXTypes Geometry PolyClipper Silo DataOutput NodeList Field Kernel Neighbor Material FileIO DataBase Boundary Physics Hydro ExternalForce Gravity Integrator Utilities NodeGenerators FieldOperations SPH CRKSPH ArtificialViscosity SVPH Mesh Damage SolidMaterial Strength ArtificialConduction"
  )

configure_file(
  SpheralCompiledPackages.py.in
  ${PROJECT_BINARY_DIR}/Spheral/SpheralCompiledPackages.py)

set(Pybind11Wraps_python_targets
  SpheralPickle.py
  PolyhedronFileUtilities.py
  )

install(FILES       ${Pybind11Wraps_python_targets}
        DESTINATION .
        )
