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
        for slide in slideSurfaces:
            contactTypes[dataBase.numNodeLists*slide[0]+slide[1]]=1
            contactTypes[slide[0]+dataBase.numNodeLists*slide[1]]=1

        result = SlideSurface%(dim)s(W, vector_of_int(contactTypes))

        result.numNodeLists=dataBase.numNodeLists

        return result
"""

for dim in dims:
    exec(SlideSurfaceFactoryString % {"dim": "%id" % dim})