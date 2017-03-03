import numpy as np
from math import *

class quaqua():

    def __init__(self):
        self.maxLevel = 6

        self.sph = []
        self.corners = []
        self.qua = []

        x = []
        y = []
        z = []
        m = []
        H = []

        self.mm = np.array([
                            [[1.0,0.0,0.0,0.0],
                             [0.5,0.5,0.0,0.0],
                             [0.5,0.0,0.5,0.0],
                             [0.5,0.0,0.0,0.5]],
                            [[0.0,1.0,0.0,0.0],
                             [0.5,0.5,0.0,0.0],
                             [0.5,0.5,-0.5,0.5],
                             [0.5,0.0,0.0,0.5]],
                            [[0.5, 0.5, -0.5, 0.5],
                             [0., 1., -0.5, 0.5],
                             [0., 1., 0., 0.],
                             [0., 0.5, 0.5, 0.]],
                            [[0.5, 0., 0.5, 0.],
                             [0., 0.5, 0.5, 0.],
                             [0., 0., 1., 0.],
                             [0., 0., 0.5, 0.5]],
                            [[1., 0., -0.5, 0.5],
                             [0.5, 0.5, -0.5, 0.5],
                             [0.5, 0., 0., 0.5],
                             [0.5, 0., -0.5, 1.]],
                            [[0., 1., -1., 1.],
                             [0.5, 0.5, -1., 1.],
                             [0.5, 0., -0.5, 1.],
                             [0.5, 0., 0., 0.5]],
                            [[0., 0., 0., 1.],
                             [0.25, 0.5, -0.75, 1.],
                             [0.5, 0., -0.5, 1.],
                             [0.5, 0., 0., 0.5]],
                            [[0., 0., 0.5, 0.5],
                             [0.25, 0.5, -0.25, 0.5],
                             [0., 1., -0.5, 0.5],
                             [0., 1., -1., 1.]]])
        self.fixpt = np.array([0.4285714, 0.2857143, 0.1428571])

        vec = np.array([0., - sqrt(3.), 0.,
                        0., 0., 0.,
                        0., 0., 1.,
                        -1., 0., 1.])
            
        level = 0
        ii = 0

        self.checkNorm(level,ii)
        vec = self.moveCenter(vec)

        vec = self.recurse(level,ii,vec)
        '''
        print self.sph
        print self.corners
        print self.qua
        '''
        
        #self.qua = self.sph
        
        #self.qua = self.qua + self.sph
        
        print "X Y Z value"
        for i in xrange(len(self.qua)):
            print self.qua[i][0], self.qua[i][1], self.qua[i][2], 1
        
        return

    def checkNorm(self,level,ii):
        b = pow(8,self.maxLevel)
        #print "This will produce %e points" % b
        
        for i in xrange(8):
            for j in xrange(4):
                b = 0
                for k in xrange(4):
                    b = b + self.mm[i][j][k]
        #if (b!=1):
        #print "b = %f,%i,%i" %(b,i,j)
        return
            
    def moveCenter(self,vec):
        for i in xrange(4):
            j = i*3
        #print "%f %f %f" % (vec[j],vec[j+1],vec[j+2])
        
        for j in xrange(4):
            i = j*3
            vec[i] += (1.0-self.fixpt[0])
            vec[i+1] += sqrt(3.0)*self.fixpt[1]
            vec[i+2] += -self.fixpt[2]
        
        return vec

    def savePositions(self,vec):
        for i in xrange(4):
            self.corners.append([vec[i*3],vec[i*3+1],vec[i*3+2]])
        self.corners.append([vec[6],vec[7],vec[8]])
        self.corners.append([vec[0],vec[1],vec[2]])
        return

    def writeCenter(self,vec):
        cntr = np.zeros(3)
        for i in xrange(3):
            cntr[i] = vec[3+i] + (1.0-self.fixpt[0])*(vec[9+i]-vec[6+i])
            cntr[i] += self.fixpt[1]*(vec[i]-vec[3+i])
            cntr[i] += self.fixpt[2]*(vec[6+i]-vec[3+i])
        
        range = 0.14
        
        if ((abs(cntr[0]) < range) and (abs(cntr[1]) < range) and (abs(cntr[2]) < range)):
            self.qua.append([cntr[0],cntr[1],cntr[2]])
        
        radius = 0.0
        
        for i in xrange(3):
            radius += cntr[i]*cntr[i]
        radius = sqrt(radius)
        
        if (radius < range):
            self.sph.append([cntr[0],cntr[1],cntr[2]])
        return

    def recurse(self,level,ii,vec):
        vec2 = np.zeros(12)
        
        if (level == self.maxLevel):
            self.writeCenter(vec)
            return vec
        else:
            level += 1
            for ii in xrange(8):
                for i in xrange(12):
                    vec2[i] = 0
                for i in xrange(4):
                    for k in xrange(3):
                        for j in xrange(4):
                            vec2[k+3*i] += self.mm[ii][i][j] * vec[k+3*j]
                vec2 = self.recurse(level,ii,vec2)
            
            return vec2

qvt = quaqua()


