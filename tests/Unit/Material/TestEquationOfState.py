from Spheral import *
from SpheralTestUtilities import *

################################################################################
title("Testing GammaLawGas EOS")
eos = GammaLawGasMKS3d(5.0/3.0, 1.0)

output("eos")
output("eos.gamma")
output("eos.molecularWeight")
output("eos.valid")

nodes = SphNodeList3d(10)
output("nodes")

rho = ScalarField3d(nodes, 1.0)
u = ScalarField3d(nodes, range(10))
output("rho")
output("rho[:]")
output("u")
output("u[:]")

P = ScalarField3d(nodes)
eos.setPressure(P, rho, u)
output("P")
output("P[:]")
