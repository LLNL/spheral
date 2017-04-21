# Test out integrating kernels that are half-filled.
from math import *
from Spheral import *

class Wintegral(ScalarFunctor):
    def __init__(self, W, ndim, useGradientAsKernel):
        assert ndim in (1, 2, 3)
        self.W = W
        self.ndim = ndim
        self.useGradientAsKernel = useGradientAsKernel
        ScalarFunctor.__init__(self)
        return
    def __call__(self, x):
        if self.useGradientAsKernel:
            result = abs(W.gradValue(x, 1.0))
        else:
            result = W.kernelValue(x, 1.0)
        if self.ndim == 1:
            return result
        elif self.ndim == 2:
            return pi*x*result
        else:
            return 2.0*pi*x*x*result

nperh = 2.0
deta = 1.0/nperh
neta = 5
etas1d, etas2d, etas3d = [], [], []
for ix in xrange(neta):
    etas1d.append(Vector1d((ix + 0.5)*deta))
    for iy in xrange(-neta + 1, neta):
        etas2d.append(Vector2d((ix + 0.5)*deta, (iy + 0.5)*deta))
        for iz in xrange(-neta + 1, neta):
            etas3d.append(Vector3d((ix + 0.5)*deta, (iy + 0.5)*deta, (iz + 0.5)*deta))

for (W, ndim, etas, zero) in ((TableKernel1d(BSplineKernel1d(), 1000), 1, etas1d, Vector1d.zero),
                              (TableKernel2d(BSplineKernel2d(), 1000), 2, etas2d, Vector2d.zero),
                              (TableKernel3d(BSplineKernel3d(), 1000), 3, etas3d, Vector3d.zero)):
    result = simpsonsIntegrationDouble(Wintegral(W, ndim, True), 0.0, W.kernelExtent, 1000)
    print "Expected half zeroth moment in %i dimensions:  %g" % (ndim, result)

    Wsum = 0.0
    W1sum = zero
    for eta in etas:
        Wi = abs(W.gradValue(eta.magnitude(), 1.0))
        Wsum += Wi
        W1sum += Wi*eta
    W1sum /= Wsum
    print "Result of summing W: ", Wsum, Wsum**(1.0/ndim), W1sum.magnitude() # , (Wsum/W.volumeNormalization)**(1.0/ndim), Wsum**(1.0/ndim)/W.volumeNormalization


