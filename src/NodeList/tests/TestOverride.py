from SpheralTestUtilities import *
import NodeList

class DummySphNodeList1d(SphNodeList1d):
    def __init__(self,
                 numInternal = 100,
                 numGhost = 0):
        SphNodeList1d(numInternal, numGhost)
        print "Instantiating dummy sph node list."
        return

nodes = SphNodeList1d(100)
output("nodes")

