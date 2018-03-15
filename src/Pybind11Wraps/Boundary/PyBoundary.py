import sys
sys.path.append("..")
from TrampolineGenerator import TrampolineGenerator

class PyBoundary(TrampolineGenerator):

    def __init__(self):
        TrampolineGenerator.__init__(self)
        self.namespaces = ["Spheral"]
        self.templates = ["Dimension"]
        return

    def setGhostNodes(self,
                      returnType = "void",
                      args = (("NodeList<Dimension>&", "nodeList"))):
        return

if __name__ == "__main__":
    PyBoundary()()
