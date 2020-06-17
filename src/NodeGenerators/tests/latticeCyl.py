
from math import *
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

class SquareCurve:
    def __init__(self,sideLength):
        self.sideLength = sideLength
    def __call__(self,t):
        a = t*4.0
        b = a
        if a < 1.0:
            Q = Vector2d(-1.0,1.0)
            P = Vector2d(1.0,1.0)
        elif a < 2.0 and a >= 1.0:
            b = a-1.0
            Q = Vector2d(1.0,1.0)
            P = Vector2d(1.0,-1.0)
        elif a < 3.0 and a >= 2.0:
            b = a-2.0
            Q = Vector2d(1.0,-1.0)
            P = Vector2d(-1.0,-1.0)
        else:
            b = a-3.0
            Q = Vector2d(-1.0,-1.0)
            P = Vector2d(-1.0,1.0)
        P = P*self.sideLength
        Q = Q*self.sideLength
        return P*b+Q*(1.0-b)

class CircleCurve:
    def __init__(self,radius):
        self.radius = radius
    def __call__(self,t):
        t = 1.0-t
        t = t*2*3.14159
        t = t - 3.14159 - 3.14159/4.0
        return Vector2d(cos(t),sin(t))*self.radius

def latticeDistribution3d(nx, ny, nz, rhofunc,
                        xmin,
                        xmax,
                        nNodePerh = 2.01):

    assert nx > 0
    assert ny > 0
    assert nz > 0

    dx = (xmax[0] - xmin[0])/nx
    dy = (xmax[1] - xmin[1])/ny
    dz = (xmax[2] - xmin[2])/nz

    n = nx*ny*nz
    volume = ((xmax[0] - xmin[0])*
                (xmax[1] - xmin[1])*
                (xmax[2] - xmin[2]))

    hx = 1.0/(nNodePerh*dx)
    hy = 1.0/(nNodePerh*dy)
    hz = 1.0/(nNodePerh*dz)
    #H0 = SymTensor3d(hx, 0.0, 0.0,
                       # 0.0, hy, 0.0,
                       # 0.0, 0.0, hz)

    x = []
    y = []
    z = []
    m = []
    H = []

    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                xx = xmin[0] + (i + 0.5)*dx
                yy = xmin[1] + (j + 0.5)*dy
                zz = xmin[2] + (k + 0.5)*dz
                x.append(xx)
                y.append(yy)
                z.append(zz)
                m.append(dx*dy*dz * rhofunc(Vector3d(xx, yy, zz)))
                H.append(0)

    return x, y, z, m, H

def profileMethod(rad,maxrad):
    ri = rad/maxrad
    ri = 1.0+(1.0-ri)
    return ri

def latticeCylindrcalGenerator3d(nx,ny,nz,
                                xmin = (-1.0, -1.0, -1.0),
                                xmax = (1.0, 1.0, 1.0),
                                centroid = (0.0, 0.0), # centroid remains a 2d vector in x-y plane
                                rho0 = 1.0,
                                nNodePerh = 2.01,
                                refineMethod = None,
                                rmax = 1.0,):

    x = []
    y = []
    z = []
    m = []
    H = []

    centroid = vectorfromtuple(centroid)
    maxd = max((vectorfromtuple((xmax[0],xmax[1]))-centroid).magnitude(),(centroid-vectorfromtuple((xmin[0],xmin[1]))).magnitude())

    # first fill the regular lattice positions  ------------------------
    nr = int(nx*0.5)
    xl,yl,zl,ml,Hl, = latticeDistribution3d(nx=nx,ny=ny,nz=nz,rhofunc=constantDensity(rho0),xmin=xmin,xmax=xmax,nNodePerh=nNodePerh)

    for i in range(len(xl)):
        xx = xl[i]
        yy = yl[i]
        zz = zl[i]
        hh = Hl[i]
        if abs(xx-centroid[0])>rmax or abs(yy-centroid[1])>rmax:
            x.append(xx)
            y.append(yy)
            z.append(zz)
            H.append(hh)


    # now fill the transfinite positions  -------------------------------
    xxmin0 = max(centroid[0]-rmax,xmin[0])
    xxmin1 = max(centroid[1]-rmax,xmin[1])
    xxmin = (xxmin0,xxmin1,xmin[2])
    xxmax0 = min(centroid[0]+rmax,xmax[0])
    xxmax1 = min(centroid[1]+rmax,xmax[1])
    xxmax = (xxmax0,xxmax1,xmax[2])

    def interpT(rad,minrad,maxrad):
        if rad > minrad and rad < maxrad:
            return 1.0-(rad-minrad)/(maxrad-minrad)
        elif rad >= maxrad:
            return 0.0
        else:
            return 1.0

    dx = (xmax[0] - xmin[0])/nx
    xxmin = (xxmin[0]-dx*0.5,xxmin[1],xmin[2])

    nx0 = int(nx*(2.0*rmax)/(xmax[0]-xmin[0]))
    dx = (xmax[0] - xmin[0])/nx
    dr = dx
    dz = dr
    r = rmax-dr*0.5
    thcount = 1e5
    while thcount > 1:
        prevv = Vector2d(0,0)
        dist = 0.0
        sq = SquareCurve(r)
        cir = CircleCurve(r)
        ntot = nx0*4
        thcount = max(int(ntot*(r/rmax)*refineMethod(r,rmax)),1)
        ti = interpT(r,rmax*0.4,rmax*0.9) # linear interp from square to circle
        for i in range(thcount):
            th = 1.0/thcount*i
            rs = sq(th)
            rc = cir(th)
            v = rs.lerp(rc,ti)+centroid
            if prevv.magnitude() > 0.0:
                dist += (v-prevv).magnitude()
            xx = v[0]
            yy = v[1]
            if xx < xxmax[0] and yy < xxmax[1] and xx >= xxmin[0] and yy >= xxmin[1]:
                nzi = int((xmax[2]-xmin[2])/dz)
                zi = xmin[2]+dz/2.0
                while zi < xmax[2]:
                    x.append(xx)
                    y.append(yy)
                    z.append(zi)
                    zi += dz
            prevv = v
        dist = dist / float(thcount)
        dr = dist
        dz = dr
        r = r-dr
    
    return x,y,z,m

x,y,z,m = latticeCylindrcalGenerator3d(nx=12,ny=12,nz=12,centroid=(0.5,0.5),refineMethod=profileMethod,rmax=0.8)
print("Created %d particles"%len(x))

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x, y,z, zdir='z', s=20, c=None, depthshade=True)
plt.show()