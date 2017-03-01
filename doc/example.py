import os, sys, shutil
from math import *
from Spheral2d import *
from GenerateNodeDistribution2d import *

title("Sedov 2d Test")

commandLine(nRadial = 50,	# number of radial bins/particles
            nTheta = 50,	# number of theta bins/particles
            rmin = 0.0,		# minimum radius
            rmax = 1.0,		# maximum radius
            nPerh = 1.51,	# particle neighbor count
            
            KernelConstructor = NBSplineKernel,
            order = 5,
            
            rho0 = 1.0,
            eps0 = 0.0,
            gamma = 5.0/3.0
            mu = 1.0,
            
            Cl = 1.0,
            Cq = 0.75,
            hmin = 1e-15,
            hmax = 1.0,
            cfl = 0.5,
            rhomin = 1e-10,
            
            HydroConstructor = SPHHydro,
            Qconstructor = MonaghanGingoldViscosity,
            balsaraCorrection = False,
            linearInExpansion = False,
            densityUpdate = RigorousSumDensity,
            HEvolution = IdealH,
            compatibleEnergy = True,
            
            IntegratorConstructor = CheapSynchronousRK2Integrator,
            steps = None,
            goalTime = 1.0,
            dt = 1e-8,
            dtMin = 1e-8,
            dtMax = None,
            dtGrowth = 2.0,
            restoreCycle = None,
            restartStep = 1000,
            
            dataDir = "dumps-sedov2d")

dataDir = os.path.join(dataDir,
                       HydroConstructor.__name__,
                       Qconstructor.__name__,
                       "nperh=%4.2f" % nPerh,
                       "compatibleEnergy=\%s" \% compatibleEnergy,
                       "nr=%i_nt=%i" % (nRadial, nTheta))

restartDir = os.path.join(dataDir, "restarts")
vizDir = os.path.join(dataDir, "visit")
restartBaseName = os.path.join(restartDir, "sedov-2d-%i" % nRadial)

if not os.path.exists(restartDir):
    os.makedirs(restartDir)
if not os.path.exists(vizDir):
    os.makedirs(vizDir)
