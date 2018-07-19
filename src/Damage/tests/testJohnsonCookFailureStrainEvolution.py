#-------------------------------------------------------------------------------
# Unit test of the Johnson-Cook failure strain evolution.  Makes sure we follow
# the exected transition through linear interpolation and the minimum value
# as a function of P/sigma.
#-------------------------------------------------------------------------------
from SolidSpheral1d import *
from SpheralTestUtilities import *
from math import *

#-------------------------------------------------------------------------------
# Identify ourselves!
#-------------------------------------------------------------------------------
title("Johnson-Cook failure strain evolution test")

#-------------------------------------------------------------------------------
# Material properties.
#-------------------------------------------------------------------------------
units = CGuS()
eos = TillotsonEquationOfState(materialName = "aluminum",
                               etamin = 0.1,
                               etamax = 10.0,
                               units = units)
strength = ConstantStrength(materialName = "aluminum",
                            units = units)
rho0 = eos.referenceDensity

#-------------------------------------------------------------------------------
# Create the NodeLists.
#-------------------------------------------------------------------------------
nodes = makeSolidNodeList("nodes", eos, strength)
from DistributeNodes import distributeNodesInRange1d
distributeNodesInRange1d([(nodes, 1000, rho0, (0.0, 1.0))],
                         nPerh = nPerh)

#-------------------------------------------------------------------------------
# The JC damage model.
#-------------------------------------------------------------------------------
D1, D2, D3, D4, D5 = 1.0, 2.0, 0.5, 0.25, 0.0
epsilondot0 = 0.1
Tcrit = -2.0
sigmamax = -4.0
efailmin = 0.01
damage = JohnsonCookDamage(nodeList = nodes,
                           D1 = D1,
                           D2 = D2,
                           D3 = D3,
                           D4 = D4,
                           D5 = D5,
                           aD1 = 0.0,       # D1 constant
                           bD1 = 0.0,
                           eps0D1 = 0.0,
                           aD2 = 0.0,       # D2 constant
                           bD2 = 0.0,
                           eps0D2 = 0.0,
                           epsilondot0 = epsilondot0,
                           Tcrit = Tcrit,
                           sigmamax = sigmamax,
                           efailmin = efailmin,
                           seed = 104899110)

# The analytic expection for the JW failure strain.
# Note we are neglecting the melt term here, which is why we set D5=0 above.
def JWfailStrain(P, S, psr):
    sigVM = sqrt(1.5*S.doubledot(S))
    sigVMinv = sigVM/(1.0e-20 + sigVM*sigVM)
    if -P*sigVMinv < Tcrit:
        return (D1 + D2*exp(D3*P*sigVMinv))*(1.0 + D4*log(psr/epsilondot0))
    elif -P < -sigmamax:
        e0 = (D1 + D2*exp(D3*Tcrit))*(1.0 + D4*log(psr/epsilondot0))
        chi = (-P*sigmaVMinv + Tcrit)/(-sigmamax*sigmaVMinv + Tcrit)
        return (1.0 - chi)*e0 + chi*efailmin
    else:
        return efailmin

#-------------------------------------------------------------------------------
# Generate the state and derivatives objects.
#-------------------------------------------------------------------------------
db = DataBase()
db.appendNodes(nodes)
packages = vector_of_Physics()
packages.append(damage)
state = State(db, packages)
derivs = StateDerivatives(db, packages)
