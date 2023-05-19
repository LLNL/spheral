from PYB11Generator import *

@PYB11template("Dimension")
class KernelIntegrationData:
    def pyinit(self):
        "Stores the data for integration"

    weight = PYB11readonly()
    ordinate = PYB11readonly()
    
    values = PYB11readonly()
    dvalues = PYB11readonly()

    index0 = PYB11readonly()
    nodeIndex0 = PYB11readonly()
    
    indices = PYB11readonly()
    nodeIndices = PYB11readonly()
    
    localIndex = PYB11readonly()
    
    normal = PYB11readonly()
    
    surfaceIndex0 = PYB11readonly()
    surfaceIndex = PYB11readonly()

    numSurfaces = PYB11readonly()

    time = PYB11readonly()
