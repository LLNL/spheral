from SpheralCompiledPackages import *

from spheralDimensions import spheralDimensions
dims = spheralDimensions()
#-------------------------------------------------------------------------------
# Convienent constructor for slide surfaces.
#-------------------------------------------------------------------------------
SlideSurfaceFactoryString = """
def makeSlideSurfaces%(dim)s(W,
                             dataBase,
                             slideSurfaces=None):

        contactTypes = [0]*(dataBase.numNodeLists**2)
        
        if slideSurfaces:
            
            # create the map nodelist --> index
            nodeLists = dataBase.nodeLists()
            nodeListMap = {}
            for i in range(dataBase.numNodeLists):        
                nodeListMap[nodeLists[i]]=i

            # dumb table this should be fixed laters
            for slide in slideSurfaces:
                nodeListi = nodeListMap[slide[0]]
                nodeListj = nodeListMap[slide[1]]
                contactTypes[dataBase.numNodeLists*nodeListi+nodeListj]=1
                contactTypes[nodeListi+dataBase.numNodeLists*nodeListj]=1

        result = SlideSurface%(dim)s(W, vector_of_int(contactTypes))

        result.numNodeLists=dataBase.numNodeLists

        return result
"""

for dim in dims:
    exec(SlideSurfaceFactoryString % {"dim": "%id" % dim})