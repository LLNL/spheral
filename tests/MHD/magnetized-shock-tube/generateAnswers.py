#-------------------------------------------------------------------------------
# This script generates answers for the magnetized shock tube configurations 
# using Jim Stone's Athena code.
#-------------------------------------------------------------------------------
from math import *
from Spheral import commandLine
#from SpheralTestUtilities import *
#from SpheralVisitDump import dumpPhysicsState

# Load the mpi module if we"re parallel.
import loadmpi
mpi, procID, numProcs = loadmpi.loadmpi()

from Spasmos import X, XYZ, Mesh, Materials, Athena

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(n = 64,    # Resolution of answer.

            # Shock tube initial configuration.
            # The tests are those suite mentioned by 
            # Dai and Woodward (DW, JCP 111, 354, 1994), 
            # Ryu and Jones (RJ, ApJ 442, 228, 1995), and 
            # Brio and Wu (BW, JCP 75, 500, 1988) and are as follows: 
            # 1. Magnetized Noh problem (DW), 
            # 2. 2D discontinuities (DW), 
            # 3. Weak 3D discontinuities (DW), 
            # 4. 3D discontinuities (RJ), 
            # 5. Vt = Bn = 0 discontinuities (DW),
            # 6. Magnetosonic rarefactions (RJ), 
            # 7. Switch-on fast shock (RJ), 
            # 8. Switch-off fast rarefaction (RJ), 
            # 9. Switch-off slow shock (RJ), 
            # 10. Switch-on slow rarefaction (RJ), 
            # 11. Slow compound structures (BW, gamma = 2), 
            # 12. Slow compound structures (RJ, gamma = 5/3), 
            # 13. Fast compound structures (RJ)
            config = None,
            goalTime = None)

# Shock tube initial conditions. 
class ShockTubeIC(object):
   """ShockTubeIC(name, eos, leftState, rightState,
                  Bx, endTime, solution = None)
   Creates a set of initial conditions for a magnetized shock tube
   problem.  leftState and rightState are tuples containing 
   (rho, P, Vx, Vy, Vz, By, Bz) for their respective sides."""

   def __init__(self, name, eos, leftState, rightState, Bx, 
                endTime, solution = None):
      self.name = name
      self.eos = eos
      (self.rhoL, self.pL, VxL, VyL, VzL, ByL, BzL) = leftState
      (self.rhoR, self.pR, VxR, VyR, VzR, ByR, BzR) = rightState
      self.vL = XYZ.Vector(VxL, VyL, VzL)
      self.vR = XYZ.Vector(VxR, VyR, VzR)
      self.BL = XYZ.Vector(Bx, ByL, BzL)
      self.BR = XYZ.Vector(Bx, ByR, BzR)
      self.endTime = endTime
      self.BCs = (0, 2, 2)
      self.solution = solution

      # Specific thermal energies.
      if type(eos) is Materials.GammaLawGasEOS:
         self.uL = self.pL / ((eos.gamma - 1)*self.rhoL)
         self.uR = self.pR / ((eos.gamma - 1)*self.rhoR)
      elif type(eos) is Materials.IsothermalEOS:
         self.uL = self.uR = 0.0
      else:
         assert(False)

   def AthenaFluid(self, N):
      """Initializes an Athena fluid according to this set of initial conditions."""

      xL = X.Vector(-0.5)
      xR = X.Vector( 0.5)
      mesh = Mesh.SMesh.Line(N, xL, xR, Nghost = Athena.Nghost)
      Nc = mesh.numCells
      fluid = Athena.Fluid('Athena gas', mesh, self.eos.gamma)

      # Fluid properties.
      xC = (xL + xR) / 2.0
      NL = Nc/2
      NR = Nc - NL
      rho = fluid.massDensity
      v = fluid.velocity
      B = fluid.magneticInduction
      p = fluid.pressure

      # Set up the left half of the shock tube.
      dx = (xC - xL) / NL
      for i in xrange(NL):
         v[i] = self.vL
         rho[i] = self.rhoL
         B[i] = self.BL
         p[i] = self.pL

      # Set up the right half.
      dx = (xR - xC) / NR
      for i in xrange(NR):
         j = i+NL
         v[j] = self.vR
         rho[j] = self.rhoR
         B[j] = self.BR
         p[j] = self.pR

      return fluid

# Configurations for this problem.
sqrt4pi = sqrt(4*pi)
gamma53 = Materials.GammaLawGasEOS('numeric', 5.0/3.0)
gamma2 = Materials.GammaLawGasEOS('numeric', 2.0)
shockTubeICs = \
   [ ShockTubeIC('Magnetized Inflow problem (DW)', gamma53,
                 (1., 20., 10., 0., 0., 1.4105, 0.), 
                 (1., 1., -10., 0., 0., 1.4105, 0.),
                 1.4105, 0.08),
     ShockTubeIC('2D discontinuities (DW)', gamma53,
                 (1., 1., 0., 0., 0., 1.4105, 0.), 
                 (0.1, 10., 0., 0., 0., 0.56419, 0.),
                 0.84628, 0.03),
     ShockTubeIC('3D weak discontinuities (DW)', gamma53,
                 (1.08, 0.95, 1.2, 0.01, 0.5, 1.0155, 0.56419), 
                 (1., 1., 0., 0., 0., 1.1284, 0.56419),
                 0.56419, 0.2),
     ShockTubeIC('3D discontinuities (RJ)', gamma53,
                 (1., 1., 0., 0., 0., 1.6926, 0.), 
                 (0.1, 10., 0., 2., 1., 0.28209, 0.),
                 0.84628, 0.035),
     ShockTubeIC('vt = Bn = 0 discontinuities (DW)', gamma53,
                 (0.1, 0.4, 50., 0., 0., -0.28209, -0.56419), 
                 (0.1, 0.2, 0., 0., 0., 0.28209, 0.56419),
                 0., 0.01),
     ShockTubeIC('Magnetosonic rarefactions (RJ)', gamma53,
                 (1, 1, -1, 0, 0, 1, 0), 
                 (1, 1, 1, 0, 0, 1, 0), 
                 0, 0.1),
     ShockTubeIC('Switch-on fast shock (RJ)', gamma53,
                 (1, 1, 0, 0, 0, 1, 0), 
                 (0.2, 0.1, 0, 0, 0, 0, 0), 
                 1, 0.15),
     ShockTubeIC('Switch-off fast rarefaction (RJ)', gamma53,
                 (0.4, 0.52467, -0.66991, 0.98263, 0.0, 0.0025293, 0.0), 
                 (1, 1, 0, 0, 0, 1, 0), 
                 1.3, 0.15),
     ShockTubeIC('Switch-off slow shock (RJ)', gamma53,
                 (0.65, 0.5, 0.667, -0.257, 0, 0.55, 0), 
                 (1, 0.75, 0.4, -0.94, 0, 0, 0), 
                 0.75, 0.15),
     ShockTubeIC('Switch-on slow rarefaction (RJ)', gamma53,
                 (1, 1, 0, 0, 0, 0, 0), 
                 (0.3, 0.2, 0, 0, 1, 1, 0), 
                 0.7, 0.16),
     ShockTubeIC('Slow compound structures (BW, gamma = 2)', gamma2,
                 (1., 1., 0., 0., 0., 1., 0.), 
                 (0.125, 0.1, 0., 0., 0., -1., 0.),
                 0.75, 0.1),
     ShockTubeIC('Slow compound structures (RJ, gamma = 5/3)', gamma53,
                 (1., 1., 0., 0., 0., 1., 0.), 
                 (0.125, 0.1, 0., 0., 0., -1., 0.),
                 0.75, 0.1),
     ShockTubeIC('Fast compound structures (RJ, gamma = 5/3)', gamma53,
                 (1., 1., 0., 0., 0., 1., 0.), 
                 (0.4, 0.4, 0., 0., 0., -1., 0.),
                 1.3, 0.16),
   ]


def runrunrun(config, t): 
   # Set things up. 
   ICs = shockTubeICs[config-1]
   fluid = ICs.AthenaFluid(n)
    
   print t
   if t is None:
      t = ICs.endTime

   # Run the problem.
   print 'Running %s to t = %g...'%(ICs.name, t)
   MHD = Athena.MHD(X.Vector)
   MHD.setBoundaries(*ICs.BCs)
   MHD.run([fluid], 0.0, endTime = t)

   # Write the data.
   xns = [u.x for u in fluid.nodePositions]
   xs = [0.5 * (xns[j] + xns[j+1]) for j in fluid.cells]
   rhos = fluid.massDensity[:]
   vxs = [u.x for u in fluid.velocity]
   vys = [u.y for u in fluid.velocity]
   vzs = [u.z for u in fluid.velocity]
   us = fluid.specificThermalEnergy[:]
   Bxs = [u.x for u in fluid.magneticInduction]
   Bys = [u.y for u in fluid.magneticInduction]
   Bzs = [u.z for u in fluid.magneticInduction]

   import pickle
   file = open('shock-tube-answers-%d-%g.dat'%(config, t), 'w')
   pickle.dump([xs, rhos, vxs, vys, vzs, us, Bxs, Bys, Bzs], file)
   file.close()

   # Only one Athena fluid may exist at a time because of silly globals.
   del fluid

if config is None:
   # Step through the myriad configurations and run Athena to generate answers 
   # that we can check against.
   for i in xrange(len(shockTubeICs)):
      runrunrun(i+1, goalTime)
else:
   runrunrun(config, goalTime)
    
