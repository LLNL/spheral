#-------------------------------------------------------------------------------
# Measure the cloud mass fraction for the blob test, as in Figure 6 from
# 
# 1.	Agertz O, Moore B, Stadel J, Potter D, Miniati F, Read J, et al. 
#       Fundamental differences between SPH and grid methods. Monthly Notices of
#       the Royal Astronomical Society. 2007; 380(3):963-978.
#       doi:10.1111/j.1365-2966.2007.12183.x.
#
# This method is generalized for both the 2D and 3D tests.
#-------------------------------------------------------------------------------
from NodeHistory import NodeHistory
import Spheral
import mpi

class CloudMassFraction(NodeHistory):

    #---------------------------------------------------------------------------
    # Constructor
    #---------------------------------------------------------------------------
    def __init__(self,
                 r0,               # initial blob radius
                 rhoThreshold,     # density cutoff
                 epsThreshold,     # specific thermal energy cutoff
                 nodes,            # blob NodeList
                 filename):        # file to write the results to
        self.r0 = r0
        self.rho0 = rhoThreshold
        self.eps0 = epsThreshold
        NodeHistory.__init__(self,
                             nodeList = nodes,
                             nodeIndicies = [],
                             sampleMethod = self.measureCloudFraction,
                             filename = filename,
                             labels = ("mfrac",
                                       "mass",
                                       "volume"))

        # Check our dimensionality
        if isinstance(nodes, Spheral.NodeList2d):
            self.ndim = 2
        elif isinstance(nodes, Spheral.NodeList3d):
            self.ndim = 3
        else:
            raise RuntimeError, "What the heck is %s?" % nodes

        # Find the starting mass of the cloud.
        self.M0 = nodes.mass().sumElements()
        return

    #---------------------------------------------------------------------------
    # Do our measurements.
    #---------------------------------------------------------------------------
    def measureCloudFraction(self, nodes, indices):
        mass = nodes.mass()
        rho = nodes.massDensity()
        eps = nodes.specificThermalEnergy()
        msum, volsum = 0.0, 0.0
        for i in xrange(nodes.numInternalNodes):
            if rho[i] > self.rho0 and eps[i] < self.eps0:
                msum += mass[i]
                volsum += mass[i]/rho[i]
        msum = mpi.allreduce(msum, mpi.SUM)
        volsum = mpi.allreduce(volsum, mpi.SUM)
        return msum/self.M0, msum, volsum
