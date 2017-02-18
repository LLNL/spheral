from Spheral import *
from SpheralTestUtilities import *

################################################################################
title("Testing MKS units")
mks = MKSUnits()

output("mks.unitL")
output("mks.unitM")
output("mks.unitT")
output("mks.ProtonMass")
output("mks.ElectronMass")
output("mks.G")
output("mks.c")

################################################################################
title("Testing CGS units")
cgs = CGSUnits()

output("cgs.unitL")
output("cgs.unitM")
output("cgs.unitT")
output("cgs.ProtonMass")
output("cgs.ElectronMass")
output("cgs.G")
output("cgs.c")

################################################################################
title("Testing Cosmological units")
cosmo = CosmologicalUnits()

output("cosmo.unitL")
output("cosmo.unitM")
output("cosmo.unitT")
output("cosmo.ProtonMass")
output("cosmo.ElectronMass")
output("cosmo.G")
output("cosmo.c")
