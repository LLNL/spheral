#-------------------------------------------------------------------------------
# Exercise the C++ unit tests of the r3d utilities.
#-------------------------------------------------------------------------------
#ATS:test(SELF, "", label="r2d/r3d unit tests.")
import Spheral
for method in ("test_polygon_to_r2d_poly",
               "test_r2d_poly_to_polygon",
               "test_polyhedron_to_r3d_poly",
               "test_r3d_poly_to_polyhedron",
               "test_clip_polygon",
               "test_orphan_polygon",
               "test_clip_polyhedron"):
    if method in dir(Spheral):
        print("Testing ", method, " : ", eval("Spheral.%s()" % method))
        assert eval("Spheral.%s()" % method) == "OK"
    else:
        print("Skipping ", method)
