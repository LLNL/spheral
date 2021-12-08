#-------------------------------------------------------------------------
# herzian solution for the stress inside a disk subjected to a point load
# y is oriented along the load-line, x perpendicular.
#-------------------------------------------------------------------------
class HerzianSolution:
    def __init__(self,
                 force,  # load
                 R):     # radius of disk

        self.sigma0 = 2.0*force/3.1415
        self.R = R

    def sigmaxx(self,X,Y):
        sigmax = [0]*len(X)
        for i in range(len(X)):
            sigmax[i] = self.sigma0*(0.5/self.R - ( (X[i]*X[i]*(self.R-Y[i]))/(X[i]*X[i]+(self.R-Y[i])**2)**2.0 + 
                                                    (X[i]*X[i]*(self.R+Y[i]))/(X[i]*X[i]+(self.R+Y[i])**2)**2.0 ))
        return sigmax

    def sigmayy(self,X,Y):
        sigmay = [0]*len(X)
        for i in range(len(X)):
            sigmay[i] = self.sigma0*(0.5/self.R - ( ((self.R-Y[i]))**3.0/(X[i]*X[i]+(self.R-Y[i])**2)**2.0 + 
                                                    ((self.R+Y[i]))**3.0/(X[i]*X[i]+(self.R+Y[i])**2)**2.0 ))
        return sigmay

    def sigmaxy(self,X,Y):
        sigmaxy = [0]*len(X)
        for i in range(len(X)):
            sigmaxy[i] = self.sigma0*( ( (X[i]*(self.R-Y[i]))**2.0/(X[i]*X[i]+(self.R-Y[i])**2) + 
                                         (X[i]*(self.R+Y[i]))**2.0/(X[i]*X[i]+(self.R+Y[i])**2) ))
        return sigmaxy