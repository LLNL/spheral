
# Solid FSISPH
#
#ATS:fsisph1 = test(           SELF, "--fsisph True --solid True --nx1 500 --cfl 0.45 --graphics None --clearDirectories True  --restartStep 20 --steps 40", label="Copper plastic-wave problem with FSISPH -- 1-D (serial)")
#ATS:fsisph2 = testif(fsisph1, SELF, "--fsisph True --solid True --nx1 500 --cfl 0.45 --graphics None --clearDirectories False --restartStep 20 --steps 20 --restoreCycle 20 --checkRestart True", label="Copper plastic-wave problem with FSISPH -- 1-D (serial) RESTART CHECK")


import os, sys, shutil
from SolidSpheral1d import *
from SpheralTestUtilities import *

from GenerateNodeProfile import GenerateNodeProfile1d
from VoronoiDistributeNodes import distributeNodes1d

title("1-D Elastic-Plastic copy piston problem")

#-------------------------------------------------------------------------------
# Maire, et. al., A nominally second-order cell-centered Lagrangian scheme for 
# simulating elastic–plastic flows on two-dimensional unstructured grids,
# Journal of Computational Physics, 235, 626-665, 2013, 10.1016/j.jcp.2012.10.017
#-------------------------------------------------------------------------------
            
commandLine(# materials properties
            rho0 = 8930.0,  # refence density
            etamin = 0.0,   # min density ratio
            etamax = 100.0, # max density ratio
            c0 = 3940.0,    # reference sound speed
            gamma0 = 2.0,   # gruniesen coefficient
            S1 = 1.49,      # some other parameter
            S2 = 0.00,      # some other parameter
            S3 = 0.00,      # some other parameter
            b = 0.0,        # some other parameter
            MW = 64.0,      # molecular weight
            P0 = 10e5,

            mu = 4.5e10,    # shear modulus
            Y0 = 9e7,       # yield strength

            v0 = 20.0,      # velocity

            # domain
            x0 = 0.0,  
            x1 = 1.0,      
            nx1 = 100,    

            # kernel things
            KernelConstructor = WendlandC2Kernel,  # (NBSplineKernel, WendlandC2Kernel, WendlandC4Kernel)
            order = 5,                             # spline order for NBSplineKernel
            nPerh = 3.01,                          # neighbors per smoothing scale
            HUpdate = IdealH,                      # (IdealH, IntegrateH) how do we update smoothing scale?
            iterateInitialH = True,                # do we want to find an ideal H to start 
            gradhCorrection = True,                # correct for adaptive-h (not implemented in FSISPH)
            hmin = 1e-10,                    
            hmax = 1.0,

            # hydro type (if all false default to SPH)
            crksph = False,   # based on conservative formulation w/ repoducing kernels
            psph = False,     # pressure based sph
            fsisph = False,   # multimaterial patching method
            gsph = False,     # Convolution-free godunov-SPH
            mfm = False,      # moving finite mass -- hopkins 2015

            # FSI parameters
            fsiRhoStabilizeCoeff = 0.0,         # diffusion operating through the vel-gradient
            fsiEpsDiffuseCoeff = 0.0,           # diffusion coeff for specific thermal energy
            fsiXSPHCoeff = 0.0,                 # ramps xsph up
            fsiInterfaceMethod = HLLCInterface, # (HLLCInterface, ModulusInterface, NoInterface)

            # CRK parameters
            correctionOrder = LinearOrder,
            volumeType = RKSumVolume,

            # GSPH parameters
            gsphReconstructionGradient = RiemannGradient,
            linearReconstruction = True,

            # general hydro parameters
            densityUpdate = IntegrateDensity,     # (RigorousSumDensity, IntegrateDensity, CorrectedSumDensity)
            correctVelocityGradient = True,       # linear correction to velocity gradient
            compatibleEnergy = True,              # conservative specific energy updat method (reccomended)
            evolveTotalEnergy = False,            # evolve total rather than specific energy
            solid = False,                        # turns on solid nodelist & hydro objects
            XSPH = False,                         # position updated w/ smoothed velocity
            epsilonTensile = 0.0,
            nTensile = 8,
            rhoMin = 0.01,
            filter = 0.00,

            # artificial viscosity
            Cl = None,
            Cq = None,
            boolReduceViscosity = False,
            nh = 5.0,
            aMin = 0.1,
            aMax = 2.0,
            Qlimiter = None,
            epsilon2 = None,
            etaCritFrac = None,
            etaFoldFrac = None,
            linearInExpansion = None,
            quadraticInExpansion = None,
            boolCullenViscosity = False,
            alphMax = 2.0,
            alphMin = 0.02,
            betaC = 0.7,
            betaD = 0.05,
            betaE = 1.0,
            fKern = 1.0/3.0,
            boolHopkinsCorrection = True,
            
            # integrator settings
            IntegratorConstructor = CheapSynchronousRK2Integrator,            
            cfl = 0.25,
            goalTime = 150e-6,
            dt = 1.0e-13,
            dtMin = 1.0e-13,
            dtMax = 0.1,
            dtGrowth = 2.0,
            steps = None,
            maxSteps = None,
            statsStep = 10,
            rigorousBoundaries = False,
            dtverbose = False,
            
            # IO settings
            clearDirectories = False,
            restoreCycle = -1,
            restartStep = 10000,
            dataDirBase = "dumps-Piston-1d-Cu",
            restartBaseName = "Piston-1d-Cu-restart",
            outputFile = "None",
            checkRestart = False,
            graphics = True,
            )

assert not(boolReduceViscosity and boolCullenViscosity)
assert not (fsisph and not solid)              # only implemented for solid
assert (mfm + fsisph + crksph + psph + gsph <= 1)    # only one hydro selection

if crksph:
    hydroname = os.path.join("CRKSPH",
                             str(volumeType),
                             str(correctionOrder))
elif psph:
    hydroname = "PSPH"
elif fsisph:
    hydroname = "FSISPH"
elif gsph:
    hydroname = os.path.join("GSPH",str(gsphReconstructionGradient))
elif mfm:
    hydroname = os.path.join("MFM",str(gsphReconstructionGradient))
else:
    hydroname = "SPH"
if solid:
    hydroname = "Solid" + hydroname

dataDir = os.path.join(dataDirBase, 
                       hydroname,
                       "nPerh=%f" % nPerh,
                       "compatibleEnergy=%s" % compatibleEnergy,
                       "Cullen=%s" % boolCullenViscosity,
                       "%i" % (nx1))
restartDir = os.path.join(dataDir, "restarts")
restartBaseName = os.path.join(restartDir, "Piston-1d-Cu-%i" % (nx1))

#-------------------------------------------------------------------------------
# Check if the necessary output directories exist.  If not, create them.
#-------------------------------------------------------------------------------
if mpi.rank == 0:
    if clearDirectories and os.path.exists(dataDir):
        shutil.rmtree(dataDir)
    if not os.path.exists(restartDir):
        os.makedirs(restartDir)
mpi.barrier()

#-------------------------------------------------------------------------------
# Material properties.
#-------------------------------------------------------------------------------
units = MKS()

eos = GruneisenEquationOfState(referenceDensity = rho0,
                               etamin = etamin,
                               etamax = etamax,
                               C0 = c0,
                               S1 = S1,
                               S2 = S2,
                               S3 = S3,
                               gamma0 = gamma0,
                               b = b,
                               atomicWeight = MW,
                               externalPressure = 0.0,
                               constants=units)

# gammaStiff = 6.0
# Pstiff = 1.6e10
# eos = StiffenedGas(gamma=gammaStiff, 
#                     P0=Pstiff, 
#                     Cv=3000,
#                     constants=units)

#eps0 = (P0+gammaStiff*Pstiff)/((gammaStiff - 1.0)*rho0)
# eos = TillotsonEquationOfState("copper",
#                                                etamin,
#                                                etamax,
#                                                units)
eps0 = eos.specificThermalEnergyForPressure(Ptarget = P0,
                                         rho =rho0,
                                         epsMin = -1e10,
                                         epsMax = 1e10,
                                         epsTol = 1e-5,
                                         Ptol = 1e-5,
                                         maxIterations = 100)

print eos.soundSpeed(rho0,eps0)
print sqrt(eos.computeDPDrho(rho0,eps0))
strengthModel = ConstantStrength(mu0 = mu,
                                 Y0 = Y0,
                                 muD = mu,
                                 YD = Y0)

#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
if KernelConstructor==NBSplineKernel:
    assert order in (3,5,7)
    WT = TableKernel(KernelConstructor(order), 1000)
else:
    WT = TableKernel(KernelConstructor(), 1000) 
kernelExtent = WT.kernelExtent
output("WT")

#-------------------------------------------------------------------------------
# Make the NodeLists.
#-------------------------------------------------------------------------------

nodes = makeSolidNodeList("nodes", eos, 
                           strength=strengthModel,
                           hmin = hmin,
                           hmax = hmax,
                           nPerh = nPerh,
                           kernelExtent = kernelExtent,
                           rhoMin = rhoMin)
nodeSet = [nodes]

#-------------------------------------------------------------------------------
# Set the node properties.
#-------------------------------------------------------------------------------
gen = GenerateNodeProfile1d(nx = nx1,
                             rho = rho0,
                             xmin = x0,
                             xmax = x1,
                             nNodePerh = nPerh)
distributeNodes1d((nodes, gen))
output("nodes.numNodes")

# Set node specific thermal energies
pos = nodes.positions()
eps = nodes.specificThermalEnergy()
rho = nodes.massDensity()
vel = nodes.velocity()
for i in xrange(nodes.numInternalNodes):
    vel[i].x = v0
    eps[i] = eps0
#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase()
output("db")
for nodes in nodeSet:
    output("db.appendNodeList(nodes)")
output("db.numNodeLists")
output("db.numFluidNodeLists")

#-------------------------------------------------------------------------------
# Construct the hydro physics object.
#-------------------------------------------------------------------------------

if crksph:
    hydro = CRKSPH(dataBase = db,
                   order = correctionOrder,
                   filter = filter,
                   cfl = cfl,
                   compatibleEnergyEvolution = compatibleEnergy,
                   evolveTotalEnergy = evolveTotalEnergy,
                   XSPH = XSPH,
                   densityUpdate = RigorousSumDensity,#densityUpdate,
                   HUpdate = HUpdate)
elif psph:
    hydro = PSPH(dataBase = db,
                 W = WT,
                 cfl = cfl,
                 compatibleEnergyEvolution = compatibleEnergy,
                 evolveTotalEnergy = evolveTotalEnergy,
                 densityUpdate = densityUpdate,
                 HUpdate = HUpdate,
                 XSPH = XSPH,
                 correctVelocityGradient = correctVelocityGradient)
elif fsisph:
    hydro = FSISPH(dataBase = db,
                   W = WT,
                   filter = filter,
                   cfl = cfl,
                   sumDensityNodeLists=[],                       
                   densityStabilizationCoefficient = fsiRhoStabilizeCoeff,
                   specificThermalEnergyDiffusionCoefficient = fsiEpsDiffuseCoeff,
                   xsphCoefficient = fsiXSPHCoeff,
                   interfaceMethod = fsiInterfaceMethod,
                   compatibleEnergyEvolution = compatibleEnergy,
                   evolveTotalEnergy = evolveTotalEnergy,
                   correctVelocityGradient = correctVelocityGradient,
                   HUpdate = HUpdate)
else:
    hydro = SPH(dataBase = db,
                W = WT,
                cfl = cfl,
                compatibleEnergyEvolution = compatibleEnergy,
                evolveTotalEnergy = evolveTotalEnergy,
                gradhCorrection = gradhCorrection,
                densityUpdate = densityUpdate,
                HUpdate = HUpdate,
                XSPH = XSPH,
                correctVelocityGradient = correctVelocityGradient,
                epsTensile = epsilonTensile,
                nTensile = nTensile)
output("hydro")
output("hydro.compatibleEnergyEvolution")
packages = [hydro]

#-------------------------------------------------------------------------------
# Tweak the artificial viscosity.
#-------------------------------------------------------------------------------
q = hydro.Q
if not Cl is None:
    q.Cl = Cl
if not Cq is None:
    q.Cq = Cq
if not linearInExpansion is None:
    q.linearInExpansion = linearInExpansion
if not quadraticInExpansion is None:
    q.quadraticInExpansion = quadraticInExpansion
if not Qlimiter is None:
    q.limiter = Qlimiter
if not epsilon2 is None:
    q.epsilon2 = epsilon2
if not etaCritFrac is None:
    q.etaCritFrac = etaCritFrac
if not etaFoldFrac is None:
    q.etaFoldFrac = etaFoldFrac
output("q")
output("q.Cl")
output("q.Cq")
output("q.limiter")
output("q.epsilon2")
output("q.linearInExpansion")
output("q.quadraticInExpansion")

#-------------------------------------------------------------------------------
# Construct the MMRV physics object.
#-------------------------------------------------------------------------------
if boolReduceViscosity:
    evolveReducingViscosityMultiplier = MorrisMonaghanReducingViscosity(q,nh,nh,aMin,aMax)
    packages.append(evolveReducingViscosityMultiplier)
elif boolCullenViscosity:
    evolveCullenViscosityMultiplier = CullenDehnenViscosity(q,WT,alphMax,alphMin,betaC,betaD,betaE,fKern,boolHopkinsCorrection)
    packages.append(evolveCullenViscosityMultiplier)



#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
xPlane0 = Plane(Vector(x0), Vector( 1.0))
xPlane1 = Plane(Vector(x1), Vector(-1.0))
xbc0 = InflowOutflowBoundary(db,xPlane0)
xbc1 = ReflectingBoundary(xPlane1)

for p in packages:
    p.appendBoundary(xbc0)
    p.appendBoundary(xbc1)

#-------------------------------------------------------------------------------
# Construct an integrator, and add the one physics package.
#-------------------------------------------------------------------------------
integrator = IntegratorConstructor(db)
for p in packages:
    integrator.appendPhysicsPackage(p)
integrator.lastDt = dt
integrator.dtGrowth = dtGrowth
if dtMin:
    integrator.dtMin = dtMin
if dtMax:
    integrator.dtMax = dtMax
integrator.rigorousBoundaries = rigorousBoundaries
integrator.verbose = dtverbose
output("integrator")
output("integrator.havePhysicsPackage(hydro)")
output("integrator.lastDt")
output("integrator.dtMin")
output("integrator.dtMax")
output("integrator.rigorousBoundaries")

#-------------------------------------------------------------------------------
# Make the problem controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator,
                            iterateInitialH=iterateInitialH,
                            kernel = WT,
                            volumeType = volumeType,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            restoreCycle = restoreCycle)
output("control")


#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
if not steps is None:
    if checkRestart:
        control.setRestartBaseName(restartBaseName + "_CHECK")
    control.step(steps)
    if checkRestart:
        control.setRestartBaseName(restartBaseName)

    # Are we doing the restart test?
    if checkRestart:
        state0 = State(db, integrator.physicsPackages())
        state0.copyState()
        control.loadRestartFile(control.totalSteps)
        state1 = State(db, integrator.physicsPackages())
        if not state1 == state0:
            raise ValueError, "The restarted state does not match!"
        else:
            print "Restart check PASSED."

else:
    control.advance(goalTime, maxSteps)
    control.dropRestartFile()

#-------------------------------------------------------------------------------
# Plot the final state.
#-------------------------------------------------------------------------------
dx1 = (x1 - x0)/nx1
h1 = 1.0/(nPerh*dx1)

# answer = SodSolution(nPoints=nx1,
#                      gamma1 = gamma1,
#                      gamma2 = gamma2,
#                      rho1 = rho1,
#                      P1 = P1,
#                      Po1 = Po1,
#                      rho2 = rho2,
#                      P2 = P2,
#                      Po2 =Po2,
#                      x0 = x0,
#                      x1 = x1,
#                      x2 = x2,
#                      h1 = 1.0/h1,
#                      h2 = 1.0/h2)



# Make a flat list from a FieldList
def createList(x):
    result = []
    for i in xrange(len(x)):
        for j in xrange(x[i].numInternalElements):
            result.append(x(i,j))
    return mpi.allreduce(result, mpi.SUM)

# Compute the simulated specific entropy.
Ss = db.solidDeviatoricStress
Sxx = db.newSolidScalarFieldList(0.0, "deviatoricStress")
Y = hydro.yieldStrength
mu = hydro.shearModulus
nodeLists = db.nodeLists()
numNodeLists = db.numNodeLists
for nodeListi in range(numNodeLists):
    
    for i in range(nodeLists[nodeListi].numInternalNodes):
        print Ss(nodeListi,i)
        print Y(nodeListi,i)
        print mu(nodeListi,i)
        Sxx[nodeListi][i]=Ss(nodeListi,i).xx

print Sxx
SxxList = createList(Sxx)
#gammaList = createList(gamma)

cs = hydro.soundSpeed
csList = createList(cs)

#rho = createList(db.fluidMassDensity)
#P = createList(hydro.pressure)
#A = [Pi/rhoi**gammai for (Pi, rhoi,gammai) in zip(P, rho,gammaList)]

# The analytic solution for the simulated entropy.
xprof = [x.x for x in createList(db.fluidPosition)]
#xans, vans, uans, rhoans, Pans, hans, gammaans = answer.solution(control.time(), xprof)
#Aans = [Pi/rhoi**gammai for (Pi, rhoi, gammai) in zip(Pans,  rhoans, gammaans)]
#csAns = [sqrt(gammai*Pi/rhoi) for (Pi, rhoi, gammai) in zip(Pans,  rhoans, gammaans)]

if graphics:
    from SpheralMatplotlib import *

    rhoPlot, velPlot, epsPlot, PPlot, HPlot = plotState(db)
    # plotAnswer(answer, control.time(),
    #            rhoPlot = rhoPlot,
    #            velPlot = velPlot,
    #            epsPlot = epsPlot,
    #            PPlot = PPlot,
    #            HPlot = HPlot)
    pE = plotEHistory(control.conserve)

    
    SPlot = plotFieldList(Sxx, winTitle="Sigmaxx")
    csPlot = plotFieldList(cs, winTitle="Sigmaxx")
    #            label = "Analytic")

    #APlot = newFigure()
    #APlot.plot(xprof, A, "ro", label="Simulation")
    #APlot.plot(xans, Aans, "k-", label="Analytic")
    #plt.title("A entropy")

    plots = [(rhoPlot, "Sod-planar-rho.png"),
             (velPlot, "Sod-planar-vel.png"),
             (epsPlot, "Sod-planar-eps.png"),
             (PPlot, "Sod-planar-P.png"),
             (HPlot, "Sod-planar-h.png"),
             (SPlot, "Sod-planar-S.png"),
             (csPlot, "Sod-planar-cs.png")]
             #(APlot, "Sod-planar-entropy.png")]
    
    # if crksph:
    #     volPlot = plotFieldList(control.RKCorrections.volume, 
    #                             winTitle = "volume",
    #                             colorNodeLists = False, plotGhosts = False)
    #     splot = plotFieldList(control.RKCorrections.surfacePoint,
    #                           winTitle = "surface point",
    #                           colorNodeLists = False)
    #     plots += [(volPlot, "Sod-planar-vol.png"),
    #                (splot, "Sod-planar-surfacePoint.png")]
    # if not (gsph or mfm):
    #     viscPlot = plotFieldList(hydro.maxViscousPressure,
    #                          winTitle = "max($\\rho^2 \pi_{ij}$)",
    #                          colorNodeLists = False)
    #     plots.append((viscPlot, "Sod-planar-viscosity.png"))
    
    #     if boolCullenViscosity:
    #         cullAlphaPlot = plotFieldList(q.ClMultiplier,
    #                                   winTitle = "Cullen alpha")
    #         cullDalphaPlot = plotFieldList(evolveCullenViscosityMultiplier.DalphaDt,
    #                                    winTitle = "Cullen DalphaDt")
    #         plots += [(cullAlphaPlot, "Sod-planar-Cullen-alpha.png"),
    #               (cullDalphaPlot, "Sod-planar-Cullen-DalphaDt.png")]

    #     if boolReduceViscosity:
    #         alphaPlot = plotFieldList(q.ClMultiplier,
    #         alphaPlot = plotFieldList(q.ClMultiplier,
    #                               winTitle = "rvAlpha",
    #                               colorNodeLists = False)

    # Make hardcopies of the plots.
    for p, filename in plots:
        savefig(p, os.path.join(dataDir, filename))

print "Energy conservation: original=%g, final=%g, error=%g" % (control.conserve.EHistory[0],
                                                                control.conserve.EHistory[-1],
                                                                (control.conserve.EHistory[-1] - control.conserve.EHistory[0])/control.conserve.EHistory[0])

'''
#-------------------------------------------------------------------------------
# If requested, write out the state in a global ordering to a file.
#-------------------------------------------------------------------------------
from SpheralGnuPlotUtilities import multiSort
mof = mortonOrderIndices(db)
mo = createList(mof)
rhoprof = createList(db.fluidMassDensity)
Pprof = createList(hydro.pressure)
vprof = [v.x for v in createList(db.fluidVelocity)]
epsprof = createList(db.fluidSpecificThermalEnergy)
hprof = [1.0/Hi.xx for Hi in createList(db.fluidHfield)]

rmin = x0
rmax = x1
if mpi.rank == 0:
    multiSort(mo, xprof, rhoprof, Pprof, vprof, epsprof, hprof)
    if outputFile != "None":
        outputFile = os.path.join(dataDir, outputFile)
        f = open(outputFile, "w")
        f.write(("#  " + 17*"'%s' " + "\n") % ("x", "rho", "P", "v", "eps", "h", "mo",
                                               "rhoans", "Pans", "vans", "hans",
                                               "x_UU", "rho_UU", "P_UU", "v_UU", "eps_UU", "h_UU"))
        for (xi, rhoi, Pi, vi, epsi, hi, mi,
             rhoansi, Pansi, vansi, hansi) in zip(xprof, rhoprof, Pprof, vprof, epsprof, hprof, mo,
                                                         rhoans, Pans, vans, hans):
            f.write((6*"%16.12e " + "%i " + 4*"%16.12e " + 6*"%i " + '\n') % 
                    (xi, rhoi, Pi, vi, epsi, hi, mi,
                     rhoansi, Pansi, vansi,  hansi,
                     unpackElementUL(packElementDouble(xi)),
                     unpackElementUL(packElementDouble(rhoi)),
                     unpackElementUL(packElementDouble(Pi)),
                     unpackElementUL(packElementDouble(vi)),
                     unpackElementUL(packElementDouble(epsi)),
                     unpackElementUL(packElementDouble(hi))))
        f.close()

    import Pnorm
    print "\tQuantity \t\tL1 \t\t\tL2 \t\t\tLinf"
    failure = False
    hD = []
    for (name, data, ans) in [("Mass Density", rhoprof, rhoans),
                              ("Pressure", Pprof, Pans),
                              ("Velocity", vprof, vans),
                              ("Thermal E", epsprof, uans),
                              ("h       ", hprof, hans)]:
        assert len(data) == len(ans)
        error = [data[i] - ans[i] for i in xrange(len(data))]
        Pn = Pnorm.Pnorm(error, xprof)
        L1 = Pn.gridpnorm(1, rmin, rmax)
        L2 = Pn.gridpnorm(2, rmin, rmax)
        Linf = Pn.gridpnorm("inf", rmin, rmax)
        print "\t%s \t\t%g \t\t%g \t\t%g" % (name, L1, L2, Linf)
        #f.write(("\t\t%g") % (L1))
        hD.append([L1,L2,Linf])
    #f.write("\n")

    print "%d\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t" % (nx1,hD[0][0],hD[1][0],hD[2][0],hD[3][0],
                                                                                hD[0][1],hD[1][1],hD[2][1],hD[3][1],
                                                                                hD[0][2],hD[1][2],hD[2][2],hD[3][2])
'''