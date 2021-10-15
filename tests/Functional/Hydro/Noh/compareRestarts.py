def readRestartData(file):
    f = open(file, 'r')
    data = f.readlines()
    f.close

    m = data[1].split()[2:]
    r = data[2].replace('(', '').replace(')', '').split()[2:]
    v = data[3].replace('(', '').replace(')', '').split()[2:]
    rho = data[4].split()[2:]
    eps = data[5].split()[2:]
    H = data[6].replace('(', '').replace(')', '').split()[2:]

    return m, r, v, rho, eps, H

def testLists(list1, list2, message):
    if list1 != list2:
        print 'Test ', message, ' FAILED'
        print list1
        print list2
        return 0
    else:
        return 1

def compareRestarts(file1, file2):
    m1, r1, v1, rho1, eps1, H1 = readRestartData(file1)
    m2, r2, v2, rho2, eps2, H2 = readRestartData(file2)

    for list in [m2, r2, v2, rho2, eps2, H2]:
        list.reverse()

    testMass = testLists(m1, m2, 'mass')
    testPosition = testLists(r1, r2, 'positions')
    testVelocity = testLists(v1, v2, 'velocity')
    testDensity = testLists(rho1, rho2, 'density')
    testEnergy = testLists(eps1, eps2, 'specific thermal energy')
    testH = testLists(H1, H2, 'H')

    if (testMass and testPosition and testVelocity and testDensity and
        testEnergy and testH):
        return 1
    else:
        print file1, ' and ', file2, ' DO NOT MATCH'
        return 0
