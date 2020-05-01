from math import *
import matplotlib.pyplot as plt
from vectormath import *

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

x=[]
y=[]
xc = []
yc = []
xs = []
ys = []
thcount = 100
radcount = 20
sq = SquareCurve(8.0)
cir = CircleCurve(8.0)

nx0 = 20

def profileMethod(rad,maxrad):
    ri = rad/maxrad
    ri = 1.0+(1.0-ri)
    return ri

def interpT(rad,minrad,maxrad):
    if rad > minrad and rad < maxrad:
        return 1.0-(rad-minrad)/(maxrad-minrad)
    elif rad >= maxrad:
        return 0.0
    else:
        return 1.0

for i in range(thcount):
    th = float(1.0/thcount*i)
    rs = sq(th)
    rc = cir(th)
    #v = interpVector(rs,rc,ti)
    xc.append(rc[0])
    yc.append(rc[1])
    xs.append(rs[0])
    ys.append(rs[1])

for j in range(radcount):
    r = radcount - j
    sq = SquareCurve(r)
    cir = CircleCurve(r)
    ntot = nx0*4
    thcount = int(ntot*(r/radcount)*profileMethod(r,radcount))
    ti = interpT(r,8.0,radcount-3) # linear interp from square to circle
    for i in range(thcount):
        th = 1.0/thcount*i
        rs = sq(th)
        rc = cir(th)
        v = rs.lerp(rc,ti)
        x.append(v[0])
        y.append(v[1])



plt.plot(x,y,".","k")
#plt.plot(xs,ys,".",color="red") 
#plt.plot(xc,yc,"+",color="blue")   
plt.show()

