from PYB11Generator import *

@PYB11template("Dimension")
class BilinearIndex:
    PYB11typedefs = """
    typedef typename %(Dimension)s::Vector Vector;
    typedef typename std::array<int, %(Dimension)s::nDim> ArrayDim;
"""
    def pyinit(self):
        "Indexing for bilinear form"
    
    @PYB11const
    def numElementsForRow(self,
                          ni = "const std::pair<int, int>&"):
        "Number of overlap elements for this i"
        return "int"
        
    @PYB11const
    def flatIndex(self,
                  ni = "const std::pair<int, int>&",
                  nj = "const std::pair<int, int>&"):
        "Return the flat index"
        return "int"

    @PYB11const
    def nodeIndex(self,
                  ni = "const std::pair<int, int>&",
                  j = "const int&"):
        "Return the node index"
        return "std::pair<int, int>"

    @PYB11const
    def numSurfacesForRow(self,
                          ni = "const std::pair<int, int>&"):
        "Number of unique surfaces for this row"
        return "int"
    
    @PYB11const
    @PYB11pycppname("surfaceIndex")
    def surfaceIndex1(self,
                      ni = "const std::pair<int, int>&",
                      normal = "const Vector&"):
        "For the row ni, get the surface index s with the given normal"
        return "int"

    @PYB11const
    @PYB11pycppname("surfaceIndex")
    def surfaceIndex2(self,
                      ni = "const std::pair<int, int>&",
                      normal = "const ArrayDim&"):
        "For the row ni, get the surface index s with the given normal"
        return "int"
    
    @PYB11const
    def normal(self,
               ni = "const std::pair<int, int>&"):
        "For the row ni, for the surface s, get the normal"
        return "const std::vector<Vector>&"
    
    @PYB11const
    def surfaceFlags(self,
                     ni = "const std::pair<int, int>&"):
        "Get the facet indices for this cell that we have flagged for integration"
        return "const std::vector<int>&"
