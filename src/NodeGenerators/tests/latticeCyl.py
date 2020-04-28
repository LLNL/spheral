
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
            rmin = None,
            rmax = None,
            centroid = (0.0,0.0),
            nNodePerh = 2.01):
    if type(rho) == float:
        rho = constantDensity(rho)
    
    centroid = Vector2d(centroid[0],centroid[1])

    k = 0
    np = 0
    if rmin is None:
        rmin = 0
    if rmax is None:
        rmax = Vector2d(xmax[0],xmax[1]).magnitude()
    deltar = rmax - rmin
    dx = (xmax[0] - xmin[0])/nx
    dy = (xmax[1] - xmin[1])/ny

    hx = 1.0/(nNodePerh*dx)
    #hy = 1.0/(nNodePerh*dy)
    #H0 = SymTensor2d(hx, 0.0, 0.0, hy)

    x = []
    y = []
    m = []
    #H = []


    ml = []
    #Hl = []
    xl = []
    yl = []

    xc = []
    yc = []
    mc = []
    #Hc = []

    for j in range(ny):
        for i in range(nx):
            xx = xmin[0] + (i + 0.5)*dx
            yy = xmin[1] + (j + 0.5)*dy
            r = sqrt(xx*xx + yy*yy)
            m0 = dx*dy*rho(Vector2d(xx, yy))
            if (r>=rmin*0.8):
                xl.append(xx)
                yl.append(yy)
                ml.append(m0)
                #Hl.append(H0)
                k = k + 1
            if (r>=rmax):
                x.append(xx)
                y.append(yy)
                m.append(m0)
                #H.append(H0)
                np = np + 1

    # Start at the outermost radius, and work our way inward.
    theta = 2*3.14159
    ri = rmax+2.0*nNodePerh/nx

    #import random
    #random.seed()

    while ri > 0:
        
        # Get the nominal delta r, delta theta, number of nodes, and mass per
        # node at this radius.
        rhoi = rho(Vector2d(ri, 0.0))
        dr = sqrt(m0/rhoi)
        arclength = theta*ri
        arcmass = arclength*dr*rhoi
        nTheta = max(1, int(arcmass/m0))
        dTheta = theta/nTheta
        mi = arcmass/nTheta
        #hi = nNodePerh*0.5*(dr + ri*dTheta)
        #Hi = SymTensor2d(1.0/hi, 0.0,
        #                    0.0, 1.0/hi)

        # Now assign the nodes for this radius.
        for i in range(nTheta):
            thetai = (i + 0.5)*dTheta
            xc.append(ri*cos(thetai))
            yc.append(ri*sin(thetai))
            mc.append(mi)
            #Hc.append(Hi)
            xi = ri*cos(thetai)
            yi = ri*sin(thetai)
            
            if(ri < rmin):
                x.append(xi)
                y.append(yi)
                m.append(mi)
                #H.append(Hi)
                np = np + 1
            elif(ri>=rmin):
                #eps = random.random()
                #func = ((ri-rmin)/deltar)**2
                func = 1-ri/(rmin-rmax) - rmax/(rmax-rmin)
                if(func>1.0):
                    func = 1.0
                if(func<0.0):
                    func = 0.0
                #if (eps <= func):
                #x.append(ri*cos(thetai))
                #y.append(ri*sin(thetai))
                #m.append(mi)
                #H.append(Hi)
                #else:
                minddr = nx
                mink = 2*k
                for j in range(k):
                    ddr = sqrt((xl[j]-xi)**2+(yl[j]-yi)**2)
                    if (ddr < minddr):
                        minddr = ddr
                        mink = j
                xi = xi+(xl[mink]-xi)*func
                yi = yi+(yl[mink]-yi)*func
                
                minddr = nx
                for j in range(np):
                    ddr = sqrt((x[j]-xi)**2 + (y[j]-yi)**2)
                    if (ddr < minddr):
                        minddr = ddr
                if(minddr > (1.0/hx)*0.5):
                    x.append(xi+(xl[mink]-xi)*func)
                    y.append(yi+(yl[mink]-yi)*func)
                    m.append(ml[mink])
                    #H.append(Hl[mink])


        # Decrement to the next radial bin inward.
        ri = max(0.0, ri - dr)
    return x, y, m

x,y,m = create(20,20,rho=1.0,xmin=(-1.0,-1.0),xmax=(1.0,1.0),rmax=0.8,rmin=0.5)

plt.plot(x,y,".")
plt.show()