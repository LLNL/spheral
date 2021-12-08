from Spheral import *
from SpheralTestUtilities import *

################################################################################
n1 = 10
n2 = 20
n3 = 30

m1 = 1.0
m2 = 2.0
m3 = 3.0

rho1 = 1.0
rho2 = 2.0
rho3 = 3.0

eps1 = 10.0
eps2 = 20.0
eps3 = 30.0

v1 = Vector1d(1.0)
v2 = Vector1d(-1.0)
v3 = Vector1d(2.0)

gamma = 5.0/3.0
mu = 1.0

################################################################################
title('1-D DataBase Test')

eos = GammaLawGasMKS1d(gamma, mu)

nodes1 = SphNodeList1d(eos, n1)
nodes2 = SphNodeList1d(eos, n2)
nodes3 = SphNodeList1d(eos, n3)

output('nodes1.numNodes')
output('nodes2.numNodes')
output('nodes3.numNodes')

nodes1.setMass(ScalarField1d(nodes1, m1))
nodes2.setMass(ScalarField1d(nodes2, m2))
nodes3.setMass(ScalarField1d(nodes3, m3))

nodes1.setMassDensity(ScalarField1d(nodes1, rho1))
nodes2.setMassDensity(ScalarField1d(nodes2, rho2))
nodes3.setMassDensity(ScalarField1d(nodes3, rho3))

nodes1.setSpecificThermalEnergy(ScalarField1d(nodes1, eps1))
nodes2.setSpecificThermalEnergy(ScalarField1d(nodes2, eps2))
nodes3.setSpecificThermalEnergy(ScalarField1d(nodes3, eps3))

nodes1.setVelocity(VectorField1d(nodes1, v1))
nodes2.setVelocity(VectorField1d(nodes2, v2))
nodes3.setVelocity(VectorField1d(nodes3, v3))

db = DataBase1d()
output('db')
output('db.appendNodeList(nodes1)')
output('db.appendNodeList(nodes2)')
output('db.appendNodeList(nodes3)')
output('db.numNodeLists')
output('db.numFluidNodeLists')

output('db.globalMass')
output('db.fluidMass')
output('db.globalPosition')
output('db.fluidMassDensity')
output('db.fluidSpecificThermalEnergy')
output('db.fluidVelocity')
output('db.fluidWeight')
#output('db.fluidHfield')

pressure = db.fluidPressure
output('pressure')
temperature = db.fluidTemperature
output('temperature')
cs = db.fluidSoundSpeed
output('cs')
pmom = db.fluidLinearMomentum
output('pmom')
