from GenerateNodeDistribution1d import GenerateNodeDistribution1d

#-------------------------------------------------------------------------------
# Specialize the 1d Node generator for spherical coordinates
#-------------------------------------------------------------------------------
class GenerateSphericalNodeDistribution1d(GenerateNodeDistribution1d):

    #---------------------------------------------------------------------------
    # Constructor
    #---------------------------------------------------------------------------
    def __init__(self,
                 nr,
                 rho,
                 rmin,
                 rmax,
                 nNodePerh = 2.01,
                 offset = 0.0,
                 rejecter = None):
        GenerateNodeDistribution1d.__init__(self,
                                            n = nr,
                                            rho = rho,
                                            xmin = rmin,
                                            xmax = rmax,
                                            nNodePerh = nNodePerh,
                                            offset = offset,
                                            rejecter = rejecter)

        # Convert the mass for spherical shells
        self.m = [mi * ri * ri for mi, ri in zip(self.m, self.x)]

        return
