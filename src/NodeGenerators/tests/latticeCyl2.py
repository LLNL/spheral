
from math import *
import matplotlib.pyplot as plt
from vectormath import *

import numpy as np
import scipy.spatial as spatial

def density(pos):
    x = pos.x
    y = pos.y
    return pos.magnitude()**(-0.5)

class constantDensity:
    def __init__(self,d):
        self.d = d
    def __call__(self,pos):
        return self.d

class weightFunction:
    def __init__(self,minr,maxr):
        self.minr = minr
        self.maxr = maxr
    def __call__(self,d):
        if d > self.maxr:
            return 1.0
        elif d > self.minr:
            return (d-self.minr)/(self.maxr-self.minr) # just linear for now
        else: 
            return 0.0

nx = ny = 50
centroid = (0.5,0.5)
xmin = (-1,-1)
xmax = (1,1)
rho = 1.0
nNodePerh = 2.01

x = []
y = []
m = []

centroid = vectorfromtuple(centroid)
maxd = max((vectorfromtuple(xmax)-centroid).magnitude(),(centroid-vectorfromtuple(xmin)).magnitude())
rmax = maxd

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

# first fill the lattice positions  ------------------------
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

lattice = np.column_stack((xl,yl))
lattree = spatial.KDTree(lattice)


# now fill the cylindrical positions -----------------------
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
    dTheta = theta/nTheta
    mi = arcmass/nTheta
    #hi = nNodePerh*0.5*(dr + ri*dTheta)
    #Hi = SymTensor2d(1.0/hi, 0.0,
    #                    0.0, 1.0/hi)

    # Now assign the nodes for this radius.
    for i in range(nTheta):
        thetai = (i + 0.5)*dTheta
        xi = ri*cos(thetai) + centroid.x
        yi = ri*sin(thetai) + centroid.y
        if xi >= xmin[0] and xi <= xmax[0] and yi >= xmin[1] and yi <= xmax[1]:
            xc.append(xi)
            yc.append(yi)
            mc.append(mi)
        #Hc.append(Hi)

    ri = max(0.0, ri - dr)

cylinder = np.column_stack((xc,yc))
cyltree = spatial.KDTree(cylinder)

minr = 0.7
maxr = 0.9

wt = weightFunction(minr,maxr)

# first grab all lattice points that are outside maxr
for i in range(len(lattice)):
    pos = Vector2d(xl[i],yl[i])
    if (pos-centroid).magnitude() > maxr:
        x.append(xl[i])
        y.append(yl[i])


n = len(cylinder)
for i in range(n):
    pos = cylinder[i]
    idx = lattree.query(pos)[1]
    lpos = lattice[idx]
    vpos = Vector2d(pos[0],pos[1])
    d = (vpos-centroid).magnitude()
    w = wt(d)

    vlpos = Vector2d(lpos[0],lpos[1])
    d = (vlpos - vpos).magnitude()
    v = (vlpos - vpos).norm()*w*d
    x.append(pos[0]+v[0])
    y.append(pos[1]+v[1])

#plt.plot(xl,yl,".",color="k")
#plt.plot(xc,yc,".",color="red")
plt.plot(x,y,"+",color="blue")   
plt.show()






#x,y,m = create(50,50,rho=1.0,xmin=(-1.0,-1.0),xmax=(1.0,1.0),centroid=(0.5,0.5))

