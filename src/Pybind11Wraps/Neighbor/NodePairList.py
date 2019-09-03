from PYB11Generator import *

class NodePairIdxType:
  def pyinit(self,
	     i_n = "int",
	     i_l = "int",
	     j_n = "int",
	     j_l = "int"):
    "Constructor"

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

  @PyB11const
  def size(self):
    "Returns the number of Node Pairs in the lsit"
    return "unsigned int"
