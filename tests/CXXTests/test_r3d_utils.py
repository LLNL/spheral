#-------------------------------------------------------------------------------
# Exercise the C++ unit tests of the r3d utilities.
#-------------------------------------------------------------------------------
#ATS:test(SELF, "", label="r2d/r3d unit tests.")
from SpheralModules import Spheral as sph
if "Testing" in dir(sph):
    for method in (sph.Testing.test_polygon_to_r2d_poly,
                   sph.Testing.test_r2d_poly_to_polygon,
                   sph.Testing.test_polyhedron_to_r3d_poly,
                   sph.Testing.test_r3d_poly_to_polyhedron,
                   sph.Testing.test_clip_polygon,
                   sph.Testing.test_orphan_polygon,
                   sph.Testing.test_clip_polyhedron):
        print "Testing ", str(method), " : ", method()
        assert method() == "OK"
    print "All r2d/r3d unit tests passed."
else:
    print "C++ testing module not built, skipping."
