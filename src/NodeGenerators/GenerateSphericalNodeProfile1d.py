from GenerateNodeProfile import GenerateNodeProfile1d

#-------------------------------------------------------------------------------
# Specialize the 1d profile Node generator for spherical coordinates
#-------------------------------------------------------------------------------
class GenerateSphericalNodeProfile1d(GenerateNodeProfile1d):

    #---------------------------------------------------------------------------
    # Constructor
    #---------------------------------------------------------------------------
    def __init__(self,
                 nr,
                 rho,
                 rmin,
                 rmax,
                 nNodePerh = 2.01,
                 numbins = 10000):
        GenerateNodeProfile1d.__init__(self,
                                       nx = nr,
                                       rho = rho,
                                       xmin = rmin,
                                       xmax = rmax,
                                       nNodePerh = nNodePerh,
                                       numbins = numbins)

        # Convert the mass for spherical shells
        self.m = [mi * ri * ri for mi, ri in zip(self.m, self.x)]

        return
