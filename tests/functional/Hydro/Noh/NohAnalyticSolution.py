################################################################################
# A class to provide the analytic solution for the Noh problem
################################################################################
class NohSolution:

    # Constructor
    def __init__(self, nDim,
                 r = None,
                 nPoints = 101,
                 gamma = 5.0/3.0,
                 R0 = 1.0,
                 rho0 = 1.0,
                 v0 = -1.0,
                 h0 = 2.01*0.02):
        self.nDim = nDim
        self.r = r
        self.nPoints = nPoints
        self.gamma = gamma
        self.R0 = R0
        self.rho0 = rho0
        self.v0 = v0
        self.h0 = h0
        return

    # Method to provide the analytic solution of the Noh problem at the
    # requested time.
    def solution(self, time,
                 r = None):

        # If the user has not specified the desired r coordinates, compute
        # them evenly between [0, R0].
        if r is None:
            if self.r is None:
                assert self.nPoints > 1
                R0 = self.R0 + self.v0*time
                dr = R0/(self.nPoints - 1)
                r = [i*dr for i in xrange(self.nPoints)]
            else:
                r = self.r
        assert not r is None

        # Prepare arrays for the values we're going to compute.
        v = []
        u = []
        rho = []
        P = []
        h = []

        # The current postion of the shock.
        rshock = -time/3.0*self.v0

        # Fill in the state values for each value of r.
        for ri in r:
            if abs(ri) <= rshock:
                v.append(0.0)
                u.append(0.5*self.v0**2)
                rho.append(self.rho0*2.0**(2*self.nDim))
                h.append(self.h0*(self.rho0/rho[-1])**(1.0/self.nDim))
            else:
                v.append(self.v0)
                u.append(0.0)
                rho.append(self.rho0*(1.0 - self.v0 * time/(abs(ri) + 1.0e-10))**(self.nDim - 1))
                h.append(self.h0)
            P.append((self.gamma - 1.0)*u[-1]*rho[-1])

        return r, v, u, rho, P, h

    # For higher dimensional cases, compute the ideal radial and
    # tangential smoothing scales.
    def hrtsolution(self, time):

        # First get the standard solution.
        r, v, u, rho, P, h = self.solution(time)

        # The current postion of the shock.
        rshock = -time/3.0*self.v0

        # For each position, compute them h's.
        hr = []
        ht = []
        for ri, rhoi in zip(r, rho):
            if ri <= rshock:
                hs = self.h0*(self.rho0/rhoi)**(1.0/self.nDim)
                hr.append(hs)
                ht.append(hs)
            else:
                hr.append(self.h0)
                ht.append(self.h0*(self.rho0/rhoi)**(1.0/(self.nDim - 1)))

        # Post-conditions.
        assert len(hr) == len(r)
        assert len(ht) == len(r)
        return r, hr, ht
