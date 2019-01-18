from PYB11Generator import *

#-------------------------------------------------------------------------------
# Tessellation template
#-------------------------------------------------------------------------------
@PYB11template("Dimension", "RealType")
class Tessellation:
    "Tessellation - A basic descriptor class for a topologically-consistent arbitrary poly(gonal/hedral) mesh."

    def pyinit(self):
        "Default constructor."
        return

    #...........................................................................
    # Methods
    def clear(self):
        "Clears the tessellation, emptying it of all data."
        return "void"

    @PYB11const
    def empty(self):
        "Returns true if the tessellation is empty (not defined), false otherwise."
        return "bool"

    def computeNodeCells(self):
        "Find the set of cells that touch each mesh node."
        return "std::vector<std::set<unsigned> >"

    def computeCellToNodes(self):
        "Collect the nodes around each cell"
        return "std::vector<std::set<unsigned> >"

    @PYB11implementation("[](Tessellation<%(Dimension)s, %(RealType)s>& self) -> std::string { std::stringstream os; os << self; return os.str(); }")
    def __str__(self):
        "output operator."
        return

    #...........................................................................
    # Attributes
    nodes = PYB11readwrite(returnpolicy = "reference_internal",
                           doc = """An array of (%(Dimension)s*numNodes) values containing components of 
node positions. The components are stored in node-major order and 
the 0th component of the ith node appears in nodes[%(Dimension)s*i].""")

    cells = PYB11readwrite(returnpolicy = "reference_internal",
                           doc = """This two-dimensional array defines the cell-face topology of the 
mesh. A cell has an arbitrary number of faces in 2D and 3D.
cells[i][j] gives the index of the jth face of the ith cell.
A negative face index indicates the actual face index is the 1's 
complement of the value (~cells[i][j]) and the face is oriented
witih an inward pointing normal for cells[i].""")

    faces = PYB11readwrite(returnpolicy = "reference",
                           doc = """This two-dimensional array defines the topology of the faces of the 
mesh. A face has an arbitrary number of nodes in 3D and 2 nodes in 2D. 
faces[i][j] gives the index of the jth node of the ith face.
Nodes for a given face are arranged counterclockwise around the face
viewed from the "positive" (outside) direction. """)

    #boundaryNodes = PYB11readwrite(doc = "Indices of all nodes that are on the boundary of the tessellation.")

    #boundaryFaces = PYB11readwrite(doc = "Indices of all faces on the boundary of the tessellation.")
        
    faceCells = PYB11readwrite(returnpolicy = "reference",
                               doc = """An array of cell indices for each face, i.e., the cells that share
the face.
For a given cell there will be either 1 or 2 cells -- the cases with 1
cell indicate a face on a boundary of the tessellation.""")

    neighborDomains = PYB11readwrite(returnpolicy = "reference",
                                     doc = """Parallel data structure: the set of neighbor domains this portion of
the tessellation is in contact with.""")

    sharedNodes = PYB11readwrite(returnpolicy = "reference",
                                 doc = """Parallel data structure: the nodes and faces this domain shares with
each neighbor domain.
NOTE: we implicitly assume that any domains of rank less than ours we
      are receiving from, while any domains of greater rank we send
      to.""")

    sharedFaces = PYB11readwrite(returnpolicy = "reference",
                                 doc = """Parallel data structure: the nodes and faces this domain shares with
each neighbor domain.
NOTE: we implicitly assume that any domains of rank less than ours we
      are receiving from, while any domains of greater rank we send
      to.""")

    #...........................................................................
    # A few handy properties that implement transformations on the Tessellation data.
    xnodes = PYB11property(getterraw="""[](const Tessellation<%(Dimension)s, %(RealType)s>& self) -> std::vector<%(RealType)s> { 
                                          const auto n = self.nodes.size()/%(Dimension)s;
                                          std::vector<double> result(n);
                                          for (auto i = 0; i < n; ++i) result[i] = self.nodes[%(Dimension)s*i];
                                          return result;
                                        }""",
                           doc = "Extract the X coordinates of the nodes")

    ynodes = PYB11property(getterraw="""[](const Tessellation<%(Dimension)s, %(RealType)s>& self) -> std::vector<%(RealType)s> { 
                                          const auto n = self.nodes.size()/%(Dimension)s;
                                          std::vector<double> result(n);
                                          for (auto i = 0; i < n; ++i) result[i] = self.nodes[%(Dimension)s*i + 1];
                                          return result;
                                        }""",
                           doc = "Extract the Y coordinates of the nodes")

    znodes = PYB11property(getterraw="""[](const Tessellation<%(Dimension)s, %(RealType)s>& self) -> std::vector<%(RealType)s> { 
                                          if (%(Dimension)s != 3) throw py::type_error("Cannot extract z component from 2D Tessellation");
                                          const auto n = self.nodes.size()/%(Dimension)s;
                                          std::vector<double> result(n);
                                          for (auto i = 0; i < n; ++i) result[i] = self.nodes[%(Dimension)s*i + 2];
                                          return result;
                                        }""",
                           doc = "Extract the Z coordinates of the nodes")

    facesAsInts = PYB11property(getterraw="""[](const Tessellation<%(Dimension)s, %(RealType)s>& self) -> std::vector<std::vector<int>> {
                                               std::vector<std::vector<int>> result;
                                               for (const auto& inds: self.faces) { result.push_back(std::vector<int>(inds.begin(), inds.end())); }
                                               return result;
                                             }""",
                                doc = "Same as 'faces' attribute, but returns vector<vector<int>> rather than vector<vector<unsigned>>")

    zoneNodes = PYB11property(getterraw="""[](const Tessellation<%(Dimension)s, %(RealType)s>& self) -> std::vector<std::vector<int>> {
                                             const auto nzones = self.cells.size();
                                             std::vector<std::vector<int>> result(nzones);
                                             if (%(Dimension)s == 2) {
                                               // In 2D we read the points out ordered counterclockwise.
                                               for (auto izone = 0; izone < nzones; ++izone) {
                                                 std::transform(self.cells[izone].begin(), self.cells[izone].end(), std::back_inserter(result[izone]),
                                                                [&](const int iface) { return iface < 0 ? self.faces[~iface][1] : self.faces[iface][0]; });
                                               }
                                             } else {
                                               // In 3D we just return the unique set of nodes for each zone.
                                               for (auto izone = 0; izone < nzones; ++izone) {
                                                 for (auto iface: self.cells[izone]) {
                                                   iface = iface < 0 ? ~iface : iface;
                                                   std::copy(self.faces[iface].begin(), self.faces[iface].end(), std::back_inserter(result[izone]));
                                                 }
                                                 std::sort(result[izone].begin(), result[izone].end());
                                                 result[izone].erase(std::unique(result[izone].begin(), result[izone].end()), result[izone].end());
                                               }
                                             }
                                             return result;
                                           }""",
                              doc = "Return the unique node IDs for each zone.  In 2D these are arranged counterclockwise around the zone.")
                                                   
