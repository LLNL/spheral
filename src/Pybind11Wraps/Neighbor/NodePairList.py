from PYB11Generator import *

class NodePairIdxType:
  def pyinit(self,
             i_n = "int",
             i_l = "int",
             j_n = "int",
             j_l = "int"):
    "Constructor"

    i_node = PYB11readwrite()
    j_node = PYB11readwrite()
    i_list = PYB11readwrite()
    j_list = PYB11readwrite()

class NodePairList:
  def pyinit(self):
    "Default Constructor"

  def push_back(self,
                nodePair = "NodePairIdxType"):
    "Push new Node Idx Data onto vector"
    return "void"

  def clear(self):
    "Clears all data from Node pair List."
    return "void"

  def size(self):
    "Returns the number of Node Pairs in the lsit"
    return "unsigned int"

  @PYB11implementation("[](const NodePairList& self) { return self.size(); }")
  def __len__(self):
    return

  @PYB11cppname("operator[]")
  @PYB11returnpolicy("reference_internal")
  @PYB11implementation('[](NodePairList& self, int i) { const int n = self.size(); return &self[(i %% n + n) %% n]; }')
  def __getitem__(self):
      return
