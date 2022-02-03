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
    field1 = ScalarField1d("Interpolated field", nodes)
    grad_field1 = VectorField1d("Interpolated gradient", nodes)
    rho1 = ScalarField1d("Interpolated rho (gather)", nodes)

    def pair_sum(i, j):
        ri, mi, rhoi, Hi = pos[i], mass[i], rho[i], H[i]
        rj, mj, rhoj, Hj = pos[j], mass[j], rho[j], H[j]
        Wijj = W(Hj*rj, Hj*ri, Hj.Determinant())
        gradWijj = W.grad(Hj*rj, Hj*ri, Hj.Determinant())
        field1[i] += mj/rhoj * field0[j] * Wijj
        grad_field1[i] += mj/rhoj * field0[j] * gradWijj
        rho1[i] += mj * Wijj

    for pair in pairs:
        pair_sum(pair.i_node, pair.j_node)
        pair_sum(pair.j_node, pair.i_node)
    for i in xrange(nodes.numInternalNodes):
        pair_sum(i, i)

    # Plot the results
    fig1 = plt.figure()
    rvals = [posi.x for posi in pos.internalValues()]
    plt.plot(rvals, field0.internalValues(), label="Actual")
    plt.plot(rvals, field1.internalValues(), label="Scatter interpolated")
    plt.legend()

    fig2 = plt.figure()
    plt.plot(rvals, [grad_field1[i].x for i in xrange(nodes.numInternalNodes)], label="Numerical gradient")
    plt.legend()

plt.show()

    
