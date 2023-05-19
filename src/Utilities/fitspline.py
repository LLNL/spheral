#-------------------------------------------------------------------------------
# Pair of methods to fit spline curves to points.  These methods are direct
# adaptations of the Numerical Recipes algorithms, translated here to python.
# I chose the interface to match a previous spline fitting package I was trying
# out, but that just didn't work.
#
# Usage:
#
# Get the coefficients for the spline fit, then you can use those to evaluate the
# spline at a given point.
# >>> coefs = fitspline(yknots, xknots)
# >>> yfit = evalspline(coefs, x)
#-------------------------------------------------------------------------------
def fitspline(y, x,
              yp0 = 1.0e30,
              yp1 = 1.0e30):
    n = len(x)
    assert len(y) == n
    y2 = [0.0]*n
    u = [0.0]*(n - 1)
    if yp0 > 0.99e30:
        y2[0] = 0.0
        u[0] = 0.0
    else:
        y2[0] = -0.5
        u[0] = 3.0/(x[1] - x[0]) * (y[1] - y[0])/(x[1] - x[0])
    for i in range(1, n-1):
        i0 = i - 1
        i1 = i + 1
        sig = (x[i] - x[i0])/(x[i1] - x[i0])
        p = sig*y2[i0] + 2.0
        assert abs(p) > 0.0
        y2[i] = (sig - 1.0)/p
        u[i] = (y[i1] - y[i])/(x[i1] - x[i]) - (y[i] - y[i0])/(x[i] - x[i0])
        u[i] = (6.0*u[i]/(x[i1] - x[i0]) - sig*u[i0])/p
    if yp1 > 0.99e30:
        qn = 0.0
        un = 0.0
    else:
        qn = 0.5
        un = 3.0/(x[n - 1] - x[n - 2]) * (yp1 - (y[n - 1] - y[n - 2])/(x[n - 1] - x[n - 2]))
    y2[n - 1] = (un - qn*u[n - 2])/(qn*y2[n - 2] + 1.0)
    for k in range(n - 2, 0, -1):
        y2[k] = y2[k]*y[k + 1] + u[k]
    return (x, y, y2)

def evalspline(coefs, x,
               dummyinterval = None):
    assert len(coefs) == 3
    xa = coefs[0]
    ya = coefs[1]
    y2a = coefs[2]
    n = len(xa)
    assert len(ya) == n
    assert len(y2a) == n
    klo = 0
    khi = n - 1
    while khi - klo > 1:
        k = (khi + klo) >> 1
        if xa[k] > x:
            khi = k
        else:
            klo = k
    h = xa[khi] - xa[klo]
    assert abs(h) > 0.0
    a = (xa[khi] - x)/h
    b = (x - xa[klo])/h
    return a*ya[klo] + b*ya[khi] + ((a**3 - a)*y2a[klo] + (b**3 - b)*y2a[khi])*h*h/6.0
