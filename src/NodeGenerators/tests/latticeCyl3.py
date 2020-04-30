
from math import *
import matplotlib.pyplot as plt
from vectormath import *

def density(pos):
    x = pos.x
    y = pos.y
    return pos.magnitude()**(-0.5)

class constantDensity:
    def __init__(self,d):
        self.d = d
    def __call__(self,pos):
        return self.d



def latticeDistribution(nr, rho0,
                            xmin,
                            xmax,
                            nNodePerh = 2.01):
        
    assert rho0 > 0
    nx = 2*nr+1
    ny = 2*nr+1

    dx = (xmax[0] - xmin[0])/nx
    dy = (xmax[1] - xmin[1])/ny

    m0 = dx*dy*rho0

    x = []
    y = []
    m = []

    for j in range(ny):
        for i in range(nx):
            xx = xmin[0] + (i + 0.5)*dx
            yy = xmin[1] + (j + 0.5)*dy
            x.append(xx)
            y.append(yy)
            m.append(m0)
    
    return x, y, m


def profileMethod(r):
    if r == 0.0:
        r = 1e-5
    return r**(-1.5)

nx = ny = 50
centroid = (0.25,0.5)
xmin = (-1,-1)
xmax = (1,1)
rho0 = 1.0
nNodePerh = 2.01

x = []
y = []
m = []

centroid = vectorfromtuple(centroid)
maxd = max((vectorfromtuple(xmax)-centroid).magnitude(),(centroid-vectorfromtuple(xmin)).magnitude())
#rmax = maxd
rmax = 1.0

#if type(rho) == float:
#    rho = constantDensity(rho)

centroid = vectorfromtuple(centroid)

# first fill the lattice positions  ------------------------
nr = int(nx*0.5)
xl,yl,ml = latticeDistribution(nr,rho0,xmin,xmax,nNodePerh)

for i in range(len(xl)):
    xx = xl[i]
    yy = yl[i]
    if abs(xx-centroid[0])>rmax or abs(yy-centroid[1])>rmax:
        x.append(xx)
        y.append(yy)

# now fill the stretched positions -----------------------
xxmin0 = max(centroid[0]-rmax,xmin[0])
xxmin1 = max(centroid[1]-rmax,xmin[1])
xxmin = (xxmin0,xxmin1)
xxmax0 = min(centroid[0]+rmax,xmax[0])
xxmax1 = min(centroid[1]+rmax,xmax[1])
xxmax = (xxmax0,xxmax1)

nr = int(nr*((xxmax[0]-xxmin[0])/(xmax[0]-xmin[0])))

xc,yc,mc = latticeDistribution(nr,rho0,xxmin,xxmax,nNodePerh)
rc = []
for i in range(len(xc)):
    rc.append(sqrt((xc[i]-centroid[0])**2+(yc[i]-centroid[0])**2))

zipped = zip(rc,xc,yc,mc)
zipped = sorted(zipped)
rc,xc,yc,mc = zip(*zipped)

nxx  = 2*nr+1
eta = (xxmax[0] - xxmin[0])/nxx
        
print("Stretching lattice...")

dr  = eta * 0.01    # this will essentially be the error in the new dumb way
r0p = 0
rp  = 0
rn  = 0
for i in range(1,len(rc)):
    #print "%d / %d" % (i,len(self.rl))
    r0 = rc[i]
    if (abs(r0-r0p)/r0>1e-10):
        sol     = r0**2*rho0/2.0
        iter    = int(10*rmax // dr)
        fn      = 0
        for j in range(iter+1):
            rj  = dr*j
            rjj = dr*(j+1)
            fj  = rj * profileMethod(rj)
            fjj = rjj * profileMethod(rjj)
            fn  = fn + 0.5*(fj+fjj)*dr
            if (fn>=sol):
                rn = rj
                break
    r0p = r0
    if (rn <= rmax):
        x.append(xc[i] * rn/r0)
        y.append(yc[i] * rn/r0)
        m.append(mc[i])

#seededMass = sum(m)

#mAdj = targetMass / seededMass
#for i in xrange(len(m)):
#    m[i] = m[i] * mAdj

#plt.plot(xl,yl,".",color="k")
#plt.plot(xc,yc,".",color="red")
plt.plot(x,y,"+",color="blue")   
plt.show()






#x,y,m = create(50,50,rho=1.0,xmin=(-1.0,-1.0),xmax=(1.0,1.0),centroid=(0.5,0.5))

