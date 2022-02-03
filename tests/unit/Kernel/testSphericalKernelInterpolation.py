#-------------------------------------------------------------------------------
# Test spherical interpolation of a constant field using the spherical kernel
#-------------------------------------------------------------------------------
from Spheral import *
from SpheralTestUtilities import fuzzyEqual
from GenerateSphericalNodeProfile1d import *
from SortAndDivideRedistributeNodes import *
import matplotlib.pyplot as plt

# Command line arguments
commandLine(n = 100,
            nPerh = 4.0,
            F0 = 2.0,
            r0 = 0.0,
            r1 = 1.0)

# Build the SphericalKernel
WT1 = TableKernel3d(BSplineKernel3d(), 500)
WT2 = TableKernel3d(WendlandC4Kernel3d(), 500)
t0 = time.time()
W1 = SphericalTableKernel(WT1)
print("Required %0.4f sec to construct SphericalTableKernel(Cubic B spline)"% (time.time() - t0))
t0 = time.time()
W2 = SphericalTableKernel(WT2)
print("Required %0.4f sec to construct SphericalTableKernel(Wendland C4)"% (time.time() - t0))

for W in (W1,):

    # Generate some points
    eos = GammaLawGasMKS1d(2.0, 1.0)
    nodes = makeFluidNodeList1d("nodes", eos, 
                                hmin = 1e-10,
                                hmax = 1e10,
                                nPerh = nPerh,
                                kernelExtent = W.etamax)
    gen = GenerateSphericalNodeProfile1d(nr = n,
                                         rho = 1.0,
                                         xmin = r0,
                                         xmax = r1,
                                         nNodePerh = nPerh)
    distributeNodes1d((nodes, gen))
    pos = nodes.positions()
    mass = nodes.mass()
    rho = nodes.massDensity()
    H = nodes.Hfield()

    # Build the DataBase and connectivity
    db = DataBase1d()
    db.appendNodeList(nodes)
    db.updateConnectivityMap(False, False, False)
    cm = db.connectivityMap()
    pairs = cm.nodePairList

    # Create our constant Field
    field0 = ScalarField1d("Initial field", nodes, F0)

    # Interpolate to the new field
    field1 = ScalarField1d("Interpolated field (gather)", nodes, 0.0)
    field2 = ScalarField1d("Interpolated field (scatter)", nodes, 0.0)
    rho1 = ScalarField1d("Interpolated rho (gather)", nodes, 0.0)
    rho2 = ScalarField1d("Interpolated rho (scatter)", nodes, 0.0)

    def pair_sum(i, j):
        ri, mi, rhoi, Hi = pos[i], mass[i], rho[i], H[i]
        rj, mj, rhoj, Hj = pos[j], mass[j], rho[j], H[j]
        Wiji = W(Hi*rj, Hi*ri, Hi.Determinant())
        Wijj = W(Hj*rj, Hj*ri, Hj.Determinant())
        field1[i] += mj/rhoj * field0[j] * Wiji
        field2[i] += mj/rhoj * field0[j] * Wijj
        rho1[i] += mj * Wiji
        rho2[i] += mj * Wijj

    for pair in pairs:
        pair_sum(pair.i_node, pair.j_node)
    for i in xrange(nodes.numInternalNodes):
        pair_sum(i, i)

    # Plot the results
    fig1 = plt.figure()
    rvals = [posi.x for posi in pos.internalValues()]
    plt.plot(rvals, field0.internalValues(), label="Actual")
    plt.plot(rvals, field1.internalValues(), label="Gather interpolated")
    plt.plot(rvals, field2.internalValues(), label="Scatter interpolated")
    plt.legend()

plt.show()

    
