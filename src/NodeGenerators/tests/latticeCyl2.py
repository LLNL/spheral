
from math import *
import matplotlib.pyplot as plt

class Vector2d:
    def __init__(self,x,y):
        self.x = float(x)
        self.y = float(y)
        self.items = [x,y]
    def __getitem__(self,i):
        return items[i]
    def magnitude(self):
        return sqrt(self.x**2+self.y**2)\

def vectorfromtuple(tup):
    return Vector2d(tup[0],tup[1])

def density(pos):
    x = pos.x
    y = pos.y
    return pos.magnitude()**(-0.5)

class constantDensity:
    def __init__(self,d):
        self.d = d
    def __call__(self,pos):
        return self.d



def create(nx, ny, rho,
            xmin = (0.0, 0.0),
            xmax = (1.0, 1.0),
            centroid = (0.5,0.5),
            weightFunction = None,
            nNodePerh = 2.01):
    x = []
    y = []
    m = []

    rmax = vectorfromtuple(xmax).magnitude()

    # create temporary lists for the cylindrical and lattice positions
    xc = []
    yc = []
    mc = []
    xl = []
    yl = []
    ml = []

    if type(rho) == float:
        rho = constantDensity(rho)
    
    centroid = vectorfromtuple(centroid)
    
    # first fill the lattice positions
    dx = (xmax[0] - xmin[0])/nx
    dy = (xmax[1] - xmin[1])/ny
    hx = 1.0/(nNodePerh*dx)

    for j in range(ny):
        for i in range(nx):
            xx = xmin[0] + (i + 0.5)*dx
            yy = xmin[1] + (j + 0.5)*dy
            r = sqrt(xx*xx + yy*yy)
            m0 = dx*dy*rho(Vector2d(xx, yy))
            xl.append(xx)
            yl.append(yy)
            ml.append(m0)

    # now fill the cylindrical positions
    theta = 2*3.14159
    ri = rmax+2.0*nNodePerh/nx

    #import random
    #random.seed()

    while ri > 0:   
        rhoi = rho(Vector2d(ri, 0.0))
        dr = sqrt(m0/rhoi)
        arclength = theta*ri
        arcmass = arclength*dr*rhoi
        nTheta = max(1, int(arcmass/m0))
        print(m0,nTheta,dr,ri)
        dTheta = theta/nTheta
        mi = arcmass/nTheta
        #hi = nNodePerh*0.5*(dr + ri*dTheta)
        #Hi = SymTensor2d(1.0/hi, 0.0,
        #                    0.0, 1.0/hi)

        # Now assign the nodes for this radius.
        for i in range(nTheta):
            thetai = (i + 0.5)*dTheta
            xi = ri*cos(thetai)
            yi = ri*sin(thetai)
            if xi >= xmin[0] and xi <= xmax[0] and yi >= xmin[1] and yi <= xmax[1]:
                xc.append(ri*cos(thetai))
                yc.append(ri*sin(thetai))
                mc.append(mi)
            #Hc.append(Hi)

        ri = max(0.0, ri - dr)

    plt.plot(xc,yc,".",color="red")
    #plt.plot(xl,yl,".",color="k")
    plt.show()

    return x, y, m

x,y,m = create(20,20,rho=1.0,xmin=(-1.0,-1.0),xmax=(1.0,1.0))

