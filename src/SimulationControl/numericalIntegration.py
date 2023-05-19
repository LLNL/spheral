#-------------------------------------------------------------------------------
# trapezoidalIntegration
# 
# Implements the trapezoid rule for evaluating definite integrals.
#-------------------------------------------------------------------------------
def trapezoidalIntegration(func,       # Function to be evaluated.
                           a,          # Beginning point.
                           b,          # Ending point.
                           numPoints,  # Number of points to use.
                           ):
    assert numPoints > 2

    dx = (b - a)/(numPoints - 1)

    result = 0.0
    for i in range(1, numPoints - 1):
        result += func(a + i*dx)

    return dx*(result + 0.5*(func(a) + func(b)))

#-------------------------------------------------------------------------------
# simpsonsRuleIntegration
# 
# Implements the Simpson's rule for evaluating definite integrals.
#-------------------------------------------------------------------------------
def simpsonsRuleIntegration(func,       # Function to be evaluated.
                            a,          # Beginning point.
                            b,          # Ending point.
                            numPoints,  # Number of points to use.
                            ):
    assert numPoints > 4 and numPoints % 2 == 1

    dx = (b - a)/(numPoints - 1)

    result1 = 0.0
    for i in range(1, numPoints, 2):
        x = a + i*dx
        result1 += func(x)

    result2 = 0.0
    for i in range(2, numPoints, 2):
        x = a + i*dx
        result2 += func(x)
        
    return dx*(func(a) + 4.0*result1 + 2.0*result2 + func(b))/3.0
