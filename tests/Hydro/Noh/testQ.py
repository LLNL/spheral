import Gnuplot
from Spheral import *

################################################################################
def testQlimit(Q, nodes, nodeID):
    neighbor = nodes.neighbor
    neighbor.setMasterList(nodeID)
    neighbor.setRefineNeighborList(nodeID)

    pos = nodes.positions
    vel = nodes.velocity
    DvDx = nodes.DvelocityDx
    H = nodes.Hfield

    x = Numeric.array([0.0]*neighbor.numRefine)
    y = Numeric.array([0.0]*neighbor.numRefine)
    vx = Numeric.array([0.0]*neighbor.numRefine)
    vy = Numeric.array([0.0]*neighbor.numRefine)
    Qvx = Numeric.array([0.0]*neighbor.numRefine)
    Qvy = Numeric.array([0.0]*neighbor.numRefine)
    velFraction = Numeric.array([0.0]*neighbor.numRefine)

    i = 0
    for neighborID in neighbor.refineNeighborList:
        ri = pos[nodeID]
        rj = pos[neighborID]

        vi = vel[nodeID]
        vj = vel[neighborID]

        rij = Vector2d(ri.x - rj.x, ri.y - rj.y)
        rijNorm = rij.unitVector()

        rji = Vector2d(rj.x - ri.x, rj.y - ri.y)
        rjiNorm = rji.unitVector()

        vij = Vector2d(vi.x - vj.x, vi.y - vj.y)

        sigi = DvDx[nodeID]
        sigj = DvDx[neighborID]

        Hi = H[nodeID]
        Hj = H[neighborID]

        DelVeli = (sigi.dotvec(rijNorm)).magnitude() / (Hi.dotvec(rijNorm)).magnitude();
        DelVelj = (sigj.dotvec(rijNorm)).magnitude() / (Hj.dotvec(rijNorm)).magnitude();

        fi = Q.limitSigma(vi, vj, rij, rijNorm, sigi, sigj);
        fj = Q.limitSigma(vj, vi, rji, rjiNorm, sigj, sigi);

        dvi = sigi.dotvec(rij)
        dvj = sigj.dotvec(rij)

        correctedVij = Vector2d(vij.x - 0.5*(fi*dvi.x + fj*dvj.x),
                                vij.y - 0.5*(fi*dvi.y + fj*dvj.y))

        x[i] = rj.x
        y[i] = rj.y
        vx[i] = -(vij.x)
        vy[i] = -(vij.y)
        Qvx[i] = -(correctedVij.x)
        Qvy[i] = -(correctedVij.y)
        velFraction[i] = correctedVij.magnitude() / (vij.magnitude() + 1e-30)

        i = i + 1

    return x, y, vx, vy, Qvx, Qvy, velFraction

################################################################################
def plotQvij(x, y, vx, vy,
             plot = Gnuplot.Gnuplot()):
    data = Gnuplot.Data(x, y, vx, vy,
                        with = 'vector')
    plot.replot(data)
    plot('set size square')
    return plot
